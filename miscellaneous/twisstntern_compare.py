#!/usr/bin/env python3
"""Generate the composite comparison figure used in the Dasha notebook.

This script compares two topology-weight datasets (CSV or tree sequence) and
produces the same 3×2 panel figure rendered in cell 37 of
`Examples/Figures_Notebook_4_DASHA.ipynb`.

Example
-------
python twisstntern_compare.py \
    --reference Examples/data_files/ne0.25.csv \
    --comparison Examples/data_files/ne0.15.csv \
    --alpha 0.05 \
    --reference-label "TRUTH" \
    --comparison-label "MODEL" \
    --output comparison_figure.pdf
"""

from __future__ import annotations

import argparse
import contextlib
import io
import gzip
import shutil
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde, chi2

from twisstntern.pipeline import detect_file_type
from twisstntern.tree_processing import trees_to_twisst_weights_unified
from twisstntern.utils import (
    dump_data,
    T1_lim, T2_lim, T3_lim,
    T2 as T2_line,
    T3 as T3_line,
    h, return_triangle_coord,
)

# ---------------------------------------------------------------------------
# Configuration defaults (mirroring notebook globals)
# ---------------------------------------------------------------------------
DATA_MODEL_COLORMAP = "viridis"
RESIDUALS_COLORMAP = "RdBu_r"
L2_COLORMAP = "magma"
HISTOGRAM_COLOR = "#A2C5F2"
KDE_COLOR = "#22223b"
HATCH_PATTERN = "+++"
HATCH_DATA = True
HATCH_RESIDUALS = False
HATCH_L2 = False


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def load_dataset(path: Path, *, axis_order: Optional[list[str]] = None,
                 taxon_names: Optional[list[str]] = None,
                 outgroup: Optional[str] = None) -> pd.DataFrame:
    """Load CSV or tree sequence into a (T1, T2, T3) DataFrame."""
    file_type = detect_file_type(str(path))

    if file_type == "csv":
        data = dump_data(str(path), axis_order=axis_order)
    elif file_type == "tree":
        # Use the built-in tree processor to convert any supported tree format
        # (ts/trees/newick/nexus) into a topology-weight dataframe.  We reuse
        # twisstntern's pipeline so the CLI stays aligned with the library.
        df = trees_to_twisst_weights_unified(
            str(path),
            taxon_names=taxon_names,
            outgroup=outgroup,
            output_file=None,
            verbose=False,
        )
        # numeric conversion / cleanup to mimic dump_data()
        df = df.apply(pd.to_numeric, errors="coerce").dropna().reset_index(drop=True)
        expected_cols = [col for col in df.columns if col.upper().startswith("T")]
        if len(expected_cols) < 3:
            raise ValueError(f"Tree sequence conversion did not yield topology columns: {df.columns}")
        data = df[[expected_cols[0], expected_cols[1], expected_cols[2]]].copy()
        data.columns = ["T1", "T2", "T3"]
        row_sums = data.sum(axis=1).replace(0, np.nan)
        data = data.div(row_sums, axis=0).dropna()
        data = data.loc[data["T2"] != data["T3"]]
    else:
        raise ValueError(f"Unsupported file type for {path}")

    return data[["T1", "T2", "T3"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Grid analysis utilities (from notebook)
# ---------------------------------------------------------------------------

def create_triangular_grid(alpha: float) -> list[dict[str, tuple[float, float]]]:
    triangles = []
    steps = int(1 / alpha)
    for k in range(steps):
        a1 = round(k * alpha, 10)
        b1 = round((k + 1) * alpha, 10)
        t2_steps = round((1 - k * alpha) / alpha)
        a3 = round(1 - (k + 1) * alpha, 10)
        b3 = round(1 - k * alpha, 10)
        for t2_step in range(t2_steps):
            a2 = round(t2_step * alpha, 10)
            b2 = round((t2_step + 1) * alpha, 10)
            if a3 >= 0:
                triangles.append({"T1": (a1, b1), "T2": (a2, b2), "T3": (a3, b3)})
            a3_second = round(a3 - alpha, 10)
            b3_second = round(b3 - alpha, 10)
            if a3_second >= 0:
                triangles.append({"T1": (a1, b1), "T2": (a2, b2), "T3": (a3_second, b3_second)})
            a3, b3 = a3_second, b3_second
    return triangles


def triangle_count(a1, b1, a2, b2, a3, b3, data: pd.DataFrame) -> int:
    cond_a1 = (a1 <= data.T1) if a1 == 0 else (a1 < data.T1)
    cond_a2 = (a2 <= data.T2) if a2 == 0 else (a2 < data.T2)
    cond_a3 = (a3 <= data.T3) if a3 == 0 else (a3 < data.T3)
    mask = (
        cond_a1 & (data.T1 <= b1) &
        cond_a2 & (data.T2 <= b2) &
        cond_a3 & (data.T3 <= b3)
    )
    return int(mask.sum())


def perform_enhanced_grid_analysis(data1: pd.DataFrame, data2: pd.DataFrame, alpha: float):
    triangles = create_triangular_grid(alpha)
    results = []
    for tri in triangles:
        count_data = triangle_count(*tri["T1"], *tri["T2"], *tri["T3"], data1)
        count_model = triangle_count(*tri["T1"], *tri["T2"], *tri["T3"], data2)
        prop_data = count_data / len(data1)
        prop_model = count_model / len(data2)
        prop_residual = prop_data - prop_model
        results.append({
            # store bin bounds plus all derived quantities for plotting
            "T1_bounds": tri["T1"],
            "T2_bounds": tri["T2"],
            "T3_bounds": tri["T3"],
            "count_data": count_data,
            "count_model": count_model,
            "prop_data": prop_data,
            "prop_model": prop_model,
            "prop_residual": prop_residual,
            "residual_squared": prop_residual ** 2,
            "count_residual": count_data - count_model,
        })
    results_df = pd.DataFrame(results)
    meaningful = (results_df["count_data"] > 0) | (results_df["count_model"] > 0)
    filtered = results_df[meaningful]

    if "l2" not in results_df.columns:
        results_df["l2"] = np.sqrt(results_df["residual_squared"])

    # Reproduce notebook metrics (L², chi², etc.)
    prop1 = results_df["prop_data"]
    prop2 = results_df["prop_model"]
    l2_distance = np.linalg.norm(prop1 - prop2) * (1 / len(triangles))

    counts1 = results_df["count_data"].astype(float)
    counts2 = results_df["count_model"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_stat = np.nansum((counts1 - counts2) ** 2 / (counts1 + 1e-8))
    dof = max(int(meaningful.sum()) - 1, 1)
    p_value = 1 - chi2.cdf(chi2_stat, dof)

    statistics = {
        "L2_distance": float(l2_distance),
        "chi2_statistic": float(chi2_stat),
        "p_value": float(p_value),
        "total_data_points": int(len(data1)),
        "total_model_points": int(len(data2)),
        "num_triangles": int(len(results_df)),
        "num_meaningful_triangles": int(meaningful.sum()),
        "num_empty_triangles": int((~meaningful).sum()),
        "degrees_freedom": int(dof),
        "mean_residual": float(filtered["prop_residual"].mean() if not filtered.empty else 0.0),
        "std_residual": float(filtered["prop_residual"].std() if not filtered.empty else 0.0),
        "max_residual": float(filtered["prop_residual"].abs().max() if not filtered.empty else 0.0),
        "mean_squared_residual": float(filtered["residual_squared"].mean() if not filtered.empty else 0.0),
    }

    return results_df, statistics


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def draw_empty_triangle(ax, xcoords: Iterable[float], ycoords: Iterable[float], *, hatch: bool) -> None:
    tri = Polygon(list(zip(xcoords, ycoords)), closed=True,
                  facecolor="white",
                  edgecolor="grey" if hatch else "none",
                  hatch=HATCH_PATTERN if hatch else None,
                  linewidth=0.5 if hatch else 0)
    ax.add_patch(tri)


def plot_ternary_base(ax, alpha: float) -> None:
    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)
    ax.plot(x_side_T2, T2_line(0, x_side_T2), color="grey", lw=1)
    ax.plot(x_side_T3, T3_line(0, x_side_T3), color="grey", lw=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="grey", lw=1)
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        ax.hlines(y=y * h, xmin=T1_lim(y)[0], xmax=T1_lim(y)[1], color="grey", lw=1)
        x2 = np.linspace(T2_lim(y)[0], T2_lim(y)[1], 100)
        ax.plot(x2, T2_line(y, x2), color="grey", lw=1)
        x3 = np.linspace(T3_lim(y)[0], T3_lim(y)[1], 100)
        ax.plot(x3, T3_line(y, x3), color="grey", lw=1)
    ax.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_heatmap_data(ax, df: pd.DataFrame, alpha: float, column: str,
                      title: str, cmap, norm: Optional[Normalize] = None,
                      value_transform=None, hatch_flag: bool = True) -> Normalize:
    plot_ternary_base(ax, alpha)
    data_values = df[column].to_numpy()
    if value_transform is not None:
        data_values = value_transform(data_values)
    if norm is None:
        vmax = data_values.max() if len(data_values) else 1
        vmin = 1 if vmax > 1 else 0
        norm = Normalize(vmin=vmin, vmax=vmax)
    df_reset = df.reset_index(drop=True)
    for idx, row in df_reset.iterrows():
        xtri, ytri, _ = return_triangle_coord(
            row['T1_bounds'][0], row['T1_bounds'][1],
            row['T2_bounds'][0], row['T2_bounds'][1],
            row['T3_bounds'][0], row['T3_bounds'][1],
        )
        has_data = (row['count_data'] > 0) or (row['count_model'] > 0)
        if not has_data:
            draw_empty_triangle(ax, xtri, ytri, hatch=hatch_flag)
        else:
            ax.fill(xtri, ytri, color=cmap(norm(data_values[idx])),
                    edgecolor='none', alpha=0.85)
    ax.set_title(title, fontsize=12, fontweight='bold')
    sns.despine(ax=ax)
    return norm


def plot_composite_comparison(results_df: pd.DataFrame,
                              statistics: dict[str, float],
                              alpha: float,
                              *,
                              reference_label: str,
                              comparison_label: str,
                              output_path: Optional[Path] = None) -> Path:
    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.3], hspace=0.45, wspace=0.55)
    ax_counts_ref = fig.add_subplot(gs[0, 0])
    ax_counts_cmp = fig.add_subplot(gs[0, 1])
    ax_counts_resid = fig.add_subplot(gs[0, 2])
    ax_l2 = fig.add_subplot(gs[1, 0])
    ax_hist = fig.add_subplot(gs[1, 1:3])
    ax_hist.yaxis.set_ticks_position('right')
    ax_hist.yaxis.set_label_position('right')

    def add_horizontal_colorbar(ax, norm, cmap, label, offset=0.05):
        fig.canvas.draw_idle()
        bbox = ax.get_position()
        cax = fig.add_axes([bbox.x0, bbox.y0 - offset, bbox.width, 0.02])
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label(label, fontsize=10)
        return cbar

    def add_vertical_colorbar(ax, norm, cmap, label, width=0.015):
        fig.canvas.draw_idle()
        bbox = ax.get_position()
        cax = fig.add_axes([bbox.x1 + 0.03, bbox.y0, width, bbox.height])
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=9)
        cbar.set_label(label, fontsize=10)
        return cbar

    # Row 1: counts with shared normalization
    cmap_counts = plt.get_cmap(DATA_MODEL_COLORMAP)
    vmax_counts = max(results_df['count_data'].max(), results_df['count_model'].max())
    norm_counts = Normalize(vmin=1 if vmax_counts > 1 else 0, vmax=vmax_counts)
    plot_heatmap_data(ax_counts_ref, results_df, alpha, 'count_data', reference_label,
                      cmap=cmap_counts, norm=norm_counts, hatch_flag=HATCH_DATA)
    plot_heatmap_data(ax_counts_cmp, results_df, alpha, 'count_model', comparison_label,
                      cmap=cmap_counts, norm=norm_counts, hatch_flag=HATCH_DATA)

    # Residuals heatmap (top right)
    residual_norm = Normalize(vmin=-results_df['count_residual'].abs().max(),
                              vmax=results_df['count_residual'].abs().max())
    cmap_resid = sns.color_palette(RESIDUALS_COLORMAP, as_cmap=True)
    plot_heatmap_data(ax_counts_resid, results_df, alpha, 'count_residual', 'Residuals',
                      cmap=cmap_resid, norm=residual_norm, hatch_flag=HATCH_RESIDUALS)

    if 'l2' not in results_df.columns:
        results_df['l2'] = np.sqrt(results_df['residual_squared'])
    vmax_l2 = results_df['l2'].max()
    norm_l2 = Normalize(vmin=1 if vmax_l2 > 1 else 0, vmax=vmax_l2)
    cmap_l2 = plt.get_cmap(L2_COLORMAP)
    plot_heatmap_data(ax_l2, results_df, alpha, 'l2', 'L² Distance',
                      cmap=cmap_l2, norm=norm_l2, hatch_flag=HATCH_L2)
    add_horizontal_colorbar(ax_counts_ref, norm_counts, cmap_counts, "Count")
    add_horizontal_colorbar(ax_counts_cmp, norm_counts, cmap_counts, "Count")
    add_horizontal_colorbar(ax_counts_resid, residual_norm, cmap_resid, "Count")
    add_vertical_colorbar(ax_l2, norm_l2, cmap_l2, r"$L^2$")

    # Residual histogram (bottom middle)
    meaningful = (results_df['count_data'] > 0) | (results_df['count_model'] > 0)
    residuals = results_df.loc[meaningful, 'count_residual'].to_numpy()
    hist_vals, bin_edges, _ = ax_hist.hist(
        residuals,
        bins=30,
        color=HISTOGRAM_COLOR,
        alpha=0.7,
        edgecolor='white',
        density=True,
    )
    if residuals.size > 1 and residuals.var() > 0:
        xg = np.linspace(residuals.min(), residuals.max(), 200)
        try:
            kde = gaussian_kde(residuals)
            ax_hist.plot(xg, kde(xg), color=KDE_COLOR, linewidth=2, label="KDE fit")
            ax_hist.legend(loc='upper right', fontsize=11)
        except np.linalg.LinAlgError:
            pass
    ax_hist.set_title('Residuals Distribution (Count)', fontsize=13, fontweight='bold')
    ax_hist.set_xlabel('Count Residual (Reference - Comparison)')
    ax_hist.set_ylabel('Density')
    sns.despine(ax=ax_hist)
    ax_hist.grid(False)
    ax_hist.yaxis.set_ticks_position('right')
    ax_hist.yaxis.set_label_position('right')

    # Summary stats annotation (same text block as the notebook figure)
    summary_lines = [
        f"L² distance: {statistics['L2_distance']:.6f}",
        f"χ² statistic: {statistics['chi2_statistic']:.6f}",
        f"χ² p-value: {statistics['p_value']:.2e}",
        f"Mean residual: {statistics['mean_residual']:.6f}",
        f"Std residual: {statistics['std_residual']:.6f}",
        f"Max |residual|: {statistics['max_residual']:.6f}",
        f"Mean squared residual: {statistics['mean_squared_residual']:.6f}",
        f"Data points: {statistics['total_data_points']}",
        f"Model points: {statistics['total_model_points']}",
        f"Grid triangles: {statistics['num_triangles']}",
        f"Meaningful triangles: {statistics['num_meaningful_triangles']}",
        f"Empty triangles: {statistics['num_empty_triangles']}",
    ]
    ax_hist.text(
        0.7,
        0.95,
        "\n".join(summary_lines),
        transform=ax_hist.transAxes,
        fontsize=9,
        va='top',
        ha='left',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, edgecolor='grey'),
    )

    fig.suptitle(f'Residual Analysis: {comparison_label} vs {reference_label}',
                 fontsize=20, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if output_path:
        output_path = output_path.resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight')
        print(f"Figure saved to {output_path}")
    return output_path or Path("twisstntern_compare.pdf")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True, type=Path,
                        help='Reference dataset (CSV/trees)')
    parser.add_argument('--comparison', required=True, type=Path,
                        help='Comparison dataset (CSV/trees)')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Grid granularity (default: 0.05)')
    parser.add_argument('--reference-label', default='Reference',
                        help='Label for reference dataset')
    parser.add_argument('--comparison-label', default='Comparison',
                        help='Label for comparison dataset')
    parser.add_argument('--output', type=Path, default=Path('twisstntern_compare.pdf'),
                        help='Output figure path (default: twisstntern_compare.pdf)')
    parser.add_argument('--taxon-names', nargs='*', help='Optional taxon names (tree inputs)')
    parser.add_argument('--outgroup', help='Optional outgroup (tree inputs)')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sns.set_context('notebook')
    sns.set_style('white')

    data_ref = load_dataset(args.reference, taxon_names=args.taxon_names, outgroup=args.outgroup)
    data_cmp = load_dataset(args.comparison, taxon_names=args.taxon_names, outgroup=args.outgroup)

    results_df, statistics = perform_enhanced_grid_analysis(data_ref, data_cmp, args.alpha)

    plot_composite_comparison(
        results_df,
        statistics,
        args.alpha,
        reference_label=args.reference_label,
        comparison_label=args.comparison_label,
        output_path=args.output,
    )


if __name__ == '__main__':
    # Example usage:
    # python twisstntern_compare.py \
    #   --reference Examples/data_files/ne0.25.csv \
    #   --comparison Examples/data_files/ne0.15.csv \
    #   --alpha 0.05 \
    #   --reference-label "TRUTH" \
    #   --comparison-label "MODEL" \
    #   --output comparison_figure.pdf
    main()
