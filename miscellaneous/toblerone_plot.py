#!/usr/bin/env python3
"""Plot ternary weights through time in a Toblerone (triangle × line) volume."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from twisstntern.utils import cartizian
TRI_HEIGHT = np.sqrt(3) / 2.0
BASE_VERTICES = np.array(
    [
        [0.0, TRI_HEIGHT],  # T1 apex (matches twisstntern)
        [-0.5, 0.0],  # T2 left base
        [0.5, 0.0],  # T3 right base
    ]
)
ALIAS_MAP = {
    "T1": {"t1", "topo1", "topology1", "topo_1"},
    "T2": {"t2", "topo2", "topology2", "topo_2"},
    "T3": {"t3", "topo3", "topology3", "topo_3"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a Toblerone-style ternary plot extruded through time. "
            "Defaults to Examples/data_files/Littorina_data.csv"
        )
    )
    parser.add_argument(
        "--file",
        default="Examples/data_files/Littorina_data.csv",
        help="CSV with three topology columns (relative path from repo root).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Figure DPI when saving/showing (default: 150).",
    )
    parser.add_argument(
        "--save",
        metavar="PATH",
        help="Optional path to save the figure instead of displaying.",
    )
    parser.add_argument(
        "--start-row",
        type=int,
        default=None,
        help="Optional zero-based row index (inclusive) to start reading from.",
    )
    parser.add_argument(
        "--end-row",
        type=int,
        default=None,
        help="Optional zero-based row index (exclusive) to stop reading.",
    )
    return parser.parse_args()


def find_columns(df: pd.DataFrame) -> list[str]:
    lower_map = {c.lower(): c for c in df.columns}
    ordering: list[str] = []
    used: set[str] = set()
    for target in ("T1", "T2", "T3"):
        found = next(
            (
                lower_map[name]
                for name in lower_map
                if name in ALIAS_MAP[target] and lower_map[name] not in used
            ),
            None,
        )
        if found is None:
            fallback = next((col for col in df.columns if col not in used), None)
            found = fallback
        if found is None:
            raise ValueError("Need at least three numeric columns for weights")
        ordering.append(found)
        used.add(found)
    return ordering


def load_weights(
    csv_path: Path, start_row: int | None, end_row: int | None
) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path, comment="#")
    cols = find_columns(df)
    numeric = df[cols].apply(pd.to_numeric, errors="coerce")
    valid_numeric = numeric.dropna()
    if start_row is not None or end_row is not None:
        valid_numeric = valid_numeric.iloc[start_row:end_row]
    if valid_numeric.empty:
        raise ValueError("No numeric data available in the requested row range")
    arr = valid_numeric.to_numpy(dtype=float)
    sums = arr.sum(axis=1, keepdims=True)
    mask = sums[:, 0] > 0
    if not np.any(mask):
        raise ValueError("All rows summed to zero; cannot normalize")
    normalized = arr[mask] / sums[mask]
    valid_rows = valid_numeric.index.to_numpy(dtype=float)[mask]
    return normalized, valid_rows


def barycentric_to_xy(weights: np.ndarray) -> np.ndarray:
    t1, t2, t3 = weights.T
    x, y = cartizian(t1, t2, t3)
    return np.column_stack((x, y))


def draw_prism(ax, height: float) -> None:
    tri = BASE_VERTICES
    for i in range(3):
        j = (i + 1) % 3
        ax.plot([tri[i, 0], tri[j, 0]], [tri[i, 1], tri[j, 1]], [0, 0], color="0.35", lw=1.2)
        ax.plot(
            [tri[i, 0], tri[j, 0]],
            [tri[i, 1], tri[j, 1]],
            [height, height],
            color="0.55",
            lw=1.0,
        )
        ax.plot(
            [tri[i, 0], tri[i, 0]],
            [tri[i, 1], tri[i, 1]],
            [0, height],
            color="0.85",
            lw=1.0,
            ls="--",
        )


def plot_toblerone(
    xy: np.ndarray,
    time_axis: np.ndarray,
    title: str,
    dpi: int,
    save_path: str | None,
) -> None:
    cmap = cm.plasma
    t_min, t_max = time_axis.min(), time_axis.max()
    if np.isclose(t_min, t_max):
        t_max = t_min + 1.0
    norm = Normalize(t_min, t_max)
    colors = cmap(norm(time_axis))

    fig = plt.figure(figsize=(9, 8), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")

    draw_prism(ax, max(float(time_axis.max()), 0.5))
    ax.scatter(xy[:, 0], xy[:, 1], time_axis, c=colors, s=14, depthshade=False)

    ax.text(0.0, TRI_HEIGHT + 0.05, 0, r"$\mathbf{T}_1$", ha="center", va="bottom", fontsize=12)
    ax.text(-0.55, -0.02, 0, r"$\mathbf{T}_2$", ha="center", va="top", fontsize=12)
    ax.text(0.55, -0.02, 0, r"$\mathbf{T}_3$", ha="center", va="top", fontsize=12)
    ax.text(0.65, -0.08, time_axis.max(), "Genome coordinates", ha="left", va="center", fontsize=11, rotation=90)

    x_bounds = (-0.55, 0.55)
    y_bounds = (-0.05, TRI_HEIGHT + 0.05)
    z_bounds = (time_axis.min(), max(time_axis.max(), time_axis.min() + 1e-6))
    ax.set(xlim=x_bounds, ylim=y_bounds, zlim=z_bounds, title=title)
    z_span = max(np.ptp(time_axis), 1.0)
    x_span = x_bounds[1] - x_bounds[0]
    y_span = y_bounds[1] - y_bounds[0]
    try:
        ax.set_box_aspect((x_span, y_span, z_span))
    except AttributeError:
        pass  # Matplotlib <3.4
    ax.view_init(elev=20, azim=220)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.grid(False)
    ax.set_axis_off()

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, pad=0.05, label="Genome coordinates (compressed)")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=dpi)
    else:
        plt.show()


def main() -> None:
    args = parse_args()
    csv_arg = Path(args.file)
    csv_path = csv_arg if csv_arg.is_absolute() else REPO_ROOT / csv_arg
    weights, row_ids = load_weights(csv_path, args.start_row, args.end_row)
    xy = barycentric_to_xy(weights)
    n = len(weights)
    if n == 0:
        raise ValueError("Selected row range produced no data points")
    if n == 1:
        time_axis = np.zeros(1, dtype=float)
    else:
        raw = row_ids - row_ids.min()
        span = raw.max() - raw.min()
        if span == 0:
            time_axis = np.zeros_like(raw)
        else:
            target_height = TRI_HEIGHT * 2.0  # elongate prism (~2× triangle height)
            time_axis = (raw / span) * target_height
    title = f"Toblerone ternary evolution — {csv_path.name}"
    plot_toblerone(xy, time_axis, title, args.dpi, args.save)


if __name__ == "__main__":
    main()
