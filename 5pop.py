import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objs as go
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pandas as pd
from twisstntern.visualization import plot
from twisstntern.visualization import plot_density_colored_radcount
from sklearn.neighbors import NearestNeighbors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from twisstntern.visualization import cartizian, h
from twisstntern.utils import (
    cartizian,
    return_triangle_coord,
    T2,
    T3,
    T1_lim,
    T2_lim,
    T3_lim,
    T3_lim_symm,
    h,
    mid_point_triangle,
    right_triangle_coordinates_list,
)

import msprime
from twisstntern_simulate.ts_processing import ts_to_twisst_weights
import random



# Hellp dear user (Dasha:), the simulation stuff are at the end of the file ~ line 650, 
# rest are helpers, and plotting functions
# have fun! :D




########################################################
# nice to have
def count_binary_tree_topologies(n, unrooted=True):
    """
    Count the number of distinct labeled binary tree topologies
    for n leaves. Set unrooted=False for rooted trees.
    """
    if n < 2:
        return 0  # No binary tree with < 2 leaves
    if n == 2:
        return 1 if not unrooted else 0  # Only rooted case possible

    k = 2 * n - 3 if not unrooted else 2 * n - 5

    if k <= 1:
        return 1  # For n = 3, unrooted => (2*3 - 5)!! = 1

    result = 1
    while k > 1:
        result *= k
        k -= 2
    return result


########################################################
########################################################
# For myself: good to know, behind the machinery used to produce the fake data- dirichlet distrbutions in the simplex

# explaining the alpha_value
# The Dirichlet distribution PDF for K-dimensional vector x = (x₁, ..., x_K), where each xᵢ > 0 and ∑xᵢ = 1, is:
#
#     f(x₁, ..., x_K; α₁, ..., α_K) = (1 / B(α)) * ∏_{i=1}^K xᵢ^{αᵢ - 1}
#
# where:
#     α = (α₁, ..., α_K) are the concentration parameters (αᵢ > 0),
#     B(α) is the multivariate Beta function:
#
#     B(α) = ∏_{i=1}^K Γ(αᵢ) / Γ(∑_{i=1}^K αᵢ)
#
# and Γ is the Gamma function (generalization of factorial).
# Dirichlet alpha value explanation:
# 
# When using alpha = np.ones(K) * alpha_value, the distribution is symmetric,
# but the *shape* of the distribution depends on the value of alpha_value:
#
#   - alpha_value = 1:   Uniform distribution over the simplex (all points equally likely)
#   - alpha_value < 1:   Points cluster near the corners (sparse vectors, one value ≈ 1)
#   - alpha_value > 1:   Points cluster near the center (balanced vectors, all values ≈ 1/K)
#
# This allows you to control how "peaked" or "flat" the distribution is across the simplex.
#
# Note: Regardless of alpha, the sampled vectors always have positive entries that sum to 1.





## Here is the function that produces the fake data
def generate_simplex_data(n_samples, n_dimensions, alpha_value=1.0):
    """
    Generate a dataset of Dirichlet-distributed points on a simplex.
    Args:
        n_samples: number of samples (rows)
        n_dimensions: number of dimensions (columns)
        alpha_value: Dirichlet concentration parameter (default 1.0 for uniform)
    Returns:
        np.ndarray of shape (n_samples, n_dimensions)
    """
    alpha = np.ones(n_dimensions) * alpha_value
    return np.random.dirichlet(alpha, size=n_samples)

# Number of samples and dimensions
n_samples = 5000
n_dimensions = 15

# Generate fake data
fake_data = generate_simplex_data(n_samples, n_dimensions)

# Check shape and first few rows
#print(fake_data.shape)  # should be (5000, 15)
#print(fake_data[:5])


##################################################################################################
##################################################################################################
# PLOTTING!!!!
##################################################################################################
##################################################################################################
  



## A helper function for the 3D plot -- not much use
def group_and_stack(data, groups):
    """
    Groups columns of data according to the provided groups.
    Args:
        data: np.ndarray of shape (n_samples, n_dimensions)
        groups: list of lists, each sublist contains column indices for a group
    Returns:
        np.ndarray of shape (n_samples, len(groups)), where each column is the sum of the columns in that group
    """
    return np.stack([data[:, group].sum(axis=1) for group in groups], axis=1)

# Example: group into 4 groups (A, B, C, D)
group_A = [0, 1, 2, 3]
group_B = [4, 5, 6, 7]
group_C = [8, 9, 10]
group_D = [11, 12, 13, 14]
groups = [group_A, group_B, group_C, group_D]

simplex_4d = group_and_stack(fake_data, groups)
#print(simplex_4d[:5])  # shape should be (5000, 4)


########################################################
## Plotting in 3D -- just for fun basically -- not much use

def plot_topology_weights_tetrahedron(topology_weights, groups, interactive=False):
    """
    Plots topology weights (points on a simplex) in a 3D tetrahedron.
    - topology_weights: np.ndarray of shape (n_samples, n_topologies)
    - groups: list of 4 lists, each containing column indices to sum for each vertex
      (e.g., [[0,2,5,6], [1,3], [4,7,8], [9,10,11,12,13,14]])
    - interactive: if True, uses Plotly; else uses matplotlib
    Each vertex of the tetrahedron is labeled with the string of the group indices.
    A point [1,0,0,0] will be plotted at the first vertex (group 0), etc.
    """
    assert len(groups) == 4, "You must provide exactly 4 groups."
    # Sum columns in each group
    grouped = np.stack([topology_weights[:, group].sum(axis=1) for group in groups], axis=1)
    # Create labels for each vertex as the string of the group indices
    vertex_labels = [str(g) for g in groups]
    # Use regular tetrahedron vertices
    V = np.array([
        [1, 1, 1],
        [-1, -1, 1],
        [-1, 1, -1],
        [1, -1, -1],
    ])
    def project(x):
        # Projects a 4D simplex point to 3D using barycentric coordinates
        # [1,0,0,0] -> V[0], [0,1,0,0] -> V[1], etc.
        return x[0]*V[0] + x[1]*V[1] + x[2]*V[2] + x[3]*V[3]
    pts3d = np.array([project(x) for x in grouped])
    if interactive:
        # --- Interactive Plotly Plot ---
        import plotly.graph_objs as go
        # Tetrahedron edges
        edge_indices = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
        edge_traces = []
        for i,j in edge_indices:
            edge_traces.append(go.Scatter3d(
                x=[V[i,0], V[j,0]], y=[V[i,1], V[j,1]], z=[V[i,2], V[j,2]],
                mode='lines', line=dict(color='black', width=3), showlegend=False
            ))
        # Tetrahedron faces (semi-transparent)
        face_indices = [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]
        face_traces = []
        for i, j, k in face_indices:
            face_traces.append(go.Mesh3d(
                x=[V[i,0], V[j,0], V[k,0]],
                y=[V[i,1], V[j,1], V[k,1]],
                z=[V[i,2], V[j,2], V[k,2]],
                color='gray', opacity=0.08, showscale=False, hoverinfo='skip', name='Tetrahedron'
            ))
        # Vertex labels
        label_traces = []
        for i, label in enumerate(vertex_labels):
            label_traces.append(go.Scatter3d(
                x=[V[i,0]], y=[V[i,1]], z=[V[i,2]],
                mode='text',
                text=[label],
                textposition='top center',
                showlegend=False,
                textfont=dict(color='black', size=18, family='Arial')
            ))
        scatter = go.Scatter3d(
            x=pts3d[:,0], y=pts3d[:,1], z=pts3d[:,2],
            mode='markers', marker=dict(size=3, opacity=0.5), name='Samples'
        )
        fig = go.Figure(data=[scatter]+edge_traces+face_traces+label_traces)
        fig.update_layout(
            scene=dict(
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showline=False, backgroundcolor='rgba(0,0,0,0)'),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showline=False, backgroundcolor='rgba(0,0,0,0)'),
                zaxis=dict(showgrid=False, zeroline=False, showticklabels=False, showline=False, backgroundcolor='rgba(0,0,0,0)'),
                aspectmode='cube',
            ),
            title='Interactive 3D Simplex (Tetrahedron) Projection',
            margin=dict(l=0, r=0, b=0, t=30)
        )
        fig.show()
    else:
        # --- Static 3D Tetrahedron Plot (matplotlib) ---
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(pts3d[:,0], pts3d[:,1], pts3d[:,2], alpha=0.5, s=10)
        # Draw tetrahedron faces
        faces = [[V[0], V[1], V[2]], [V[0], V[1], V[3]], [V[0], V[2], V[3]], [V[1], V[2], V[3]]]
        poly3d = Poly3DCollection(faces, alpha=0.02, facecolor='gray', edgecolor='k')
        ax.add_collection3d(poly3d)
        # Draw tetrahedron edges
        for i in range(4):
            for j in range(i+1, 4):
                ax.plot([V[i,0], V[j,0]], [V[i,1], V[j,1]], [V[i,2], V[j,2]], 'k-', lw=1)
        # Vertex labels
        for i, label in enumerate(vertex_labels):
            ax.text(V[i,0], V[i,1], V[i,2], label, color='black', fontsize=14, fontweight='bold', ha='center', va='center')
        # Hide axes, grid, and ticks for a clean look
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig("simplex_3d.png", dpi=150)
        # plt.show()  # Uncomment to display interactively


########################################################
########################################################
# Now the real deal: Projecting to ternary coordinates
# Thats the relevant plot!
def project_to_ternary(data, t1_idx, t_focal_idx):
    """
    Project n-dimensional simplex data to ternary coordinates for plotting.
    - data: np.ndarray or pd.DataFrame, shape (n_samples, n_topologies)
    - t1_idx: int, index of T1 (top vertex)
    - t_focal_idx: int, index of T2 (left vertex)
    Returns: pd.DataFrame with columns ['T1', 'T2', 'T3']
    """
    T1 = data[:, t1_idx]
    T2 = data[:, t_focal_idx]
    mask = [i for i in range(data.shape[1]) if i not in (t1_idx, t_focal_idx)]
    T3 = data[:, mask].sum(axis=1)
    df = pd.DataFrame({'T1': T1, 'T2': T2, 'T3': T3})
    # Normalize (should already sum to 1, but just in case)
    df = df.div(df.sum(axis=1), axis=0)
    return df

# Example usage:
# ternary_data = project_to_ternary(your_data, t1_idx=0, t_focal_idx=1)
# plot(ternary_data, granularity='fine', file_name='my_ternary_plot')


### you can ignore these two functions
def plot_ternary_projection_old(data, t1_idx, t_focal_idx, file_name="ternary_projection"):
    """
    Project n-dimensional simplex data to ternary coordinates and plot with density coloring.
    - data: np.ndarray or pd.DataFrame, shape (n_samples, n_topologies)
    - t1_idx: int, index of T1 (top vertex)
    - t_focal_idx: int, index of T2 (left vertex)
    - file_name: output file name prefix
    """
    T1 = data[:, t1_idx]
    T2 = data[:, t_focal_idx]
    mask = [i for i in range(data.shape[1]) if i not in (t1_idx, t_focal_idx)]
    T3 = data[:, mask].sum(axis=1)
    df = pd.DataFrame({'T1': T1, 'T2': T2, 'T3': T3})
    df = df.div(df.sum(axis=1), axis=0)
    plot_density_colored_radcount(df, file_name)



def plot_ternary_projection_with_expectation_old(data, t1_idx, t_focal_idx, expectation, file_name="manifold_projection"):
    """
    Project n-dimensional simplex data to ternary coordinates with a nonlinear T_focal axis based on expectation.
    - data: np.ndarray or pd.DataFrame, shape (n_samples, n_topologies)
    - t1_idx: int, index of T1 (top vertex)
    - t_focal_idx: int, index of T2 (left vertex)
    - expectation: float, expected value for T_focal (0 < expectation < 1)
    - file_name: output file name prefix
    The median vertical line in the triangle will correspond to the expectation value.
    Points with T_focal = expectation will be on this line; others are mapped piecewise linearly.
    """
    import matplotlib.pyplot as plt
    from twisstntern.visualization import cartizian, h
    T1 = data[:, t1_idx]
    T_focal = data[:, t_focal_idx]
    mask = [i for i in range(data.shape[1]) if i not in (t1_idx, t_focal_idx)]
    T_other = data[:, mask].sum(axis=1)
    # Piecewise linear transform for T_focal axis
    T2 = np.zeros_like(T_focal)
    for i, val in enumerate(T_focal):
        if val <= expectation:
            T2[i] = 0.5 * val / expectation if expectation > 0 else 0.0
        else:
            T2[i] = 0.5 + 0.5 * (val - expectation) / (1 - expectation) if expectation < 1 else 1.0
    # T1 remains, T2 is transformed, T3 = 1 - T1 - T2
    T3 = 1 - T1 - T2
    df = pd.DataFrame({'T1': T1, 'T2': T2, 'T3': T3})
    df = df.div(df.sum(axis=1), axis=0)
    # Plot with density coloring
    fig = plot_density_colored_radcount(df, file_name)
    # Get the correct axis from the figure
    ax = fig.gca()
    # The median vertical line in the triangle is at x=0
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, zorder=5)
    # Place the annotation just below the triangle, using figure coordinates
    fig.text(0.5, 0.04, f"Expectation: {expectation:.2f}", color='red', fontsize=13, ha='center', va='top', fontweight='bold')
    # Save without tight_layout to avoid cutting off annotation
    fig.savefig(f"{file_name}_expectation_{expectation:.2f}.png", dpi=150, bbox_inches='tight')
    # plt.show()  # Uncomment to display interactively
    return fig



# helper function for plotting
def draw_grey_grid_lines(ax, alpha=0.1):
    """
    Draw grey grid lines on ternary plot, copied from simple_density_plot.py
    Uses fixed granularity of 0.1 as requested.
    Grid lines are drawn with zorder=1 to ensure they appear under data points.
    """
    # Draw grid lines using twisstntern functions (all grey, under data points)
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        # T1 lines (horizontal)
        ax.hlines(
            y=y * h,
            xmin=T1_lim(y)[0],
            xmax=T1_lim(y)[1],
            color="grey",
            linewidth=1,
            zorder=1,
        )
        # T2 lines
        x2 = np.linspace(T2_lim(y)[0], T2_lim(y)[1], 100)
        ax.plot(x2, T2(y, x2), color="grey", linewidth=1, zorder=1)
        # T3 lines
        x3 = np.linspace(T3_lim(y)[0], T3_lim(y)[1], 100)
        ax.plot(x3, T3(y, x3), color="grey", linewidth=1, zorder=1)

    # Note: intentionally no central vertical median line

def plot_ternary_projection(
    data, t1_idx, t_focal_idx, expectation=None, colormap=None, file_name=None, ax=None
):
    """
    Plot ternary coordinate grid with data points colored by local density,
    and add a median line for the expectation value of T_focal.
    Also adds a legend for T1, T2, T3 indices.

    Parameters:
    - data: np.ndarray (n_samples, n_topologies) or DataFrame (n_samples, n_topologies) or DataFrame with columns T1, T2, T3
    - file_name: output file name prefix
    - t1_idx: index of T1 in the original data
    - t_focal_idx: index of T2 (focal) in the original data
    - expectation: float or None, expected value for T_focal (default None, uses 0.5 and hides annotation)
    - colormap: string or matplotlib colormap, e.g. 'viridis_r', 'plasma', etc.
    """
    import pandas as pd
    import numpy as np
    
    # Handle default expectation
    show_expectation_legend = True
    if expectation is None:
        expectation = 0.5
        show_expectation_legend = False
    # Always extract T1, T2, T3 from the input, regardless of format
    if isinstance(data, pd.DataFrame) and set(['T1', 'T2', 'T3']).issubset(data.columns):
        arr = data[['T1', 'T2', 'T3']].values
        T1 = arr[:, 0]
        T2 = arr[:, 1]
    else:
        if isinstance(data, pd.DataFrame):
            arr = data.values
        else:
            arr = data
        T1 = arr[:, t1_idx]
        T2 = arr[:, t_focal_idx]
    # Piecewise transform for T2 axis
    T2_new = np.zeros_like(T2)
    for i, val in enumerate(T2):
        if val <= expectation:
            T2_new[i] = 0.5 * val / expectation if expectation > 0 else 0.0
        else:
            T2_new[i] = 0.5 + 0.5 * (val - expectation) / (1 - expectation) if expectation < 1 else 1.0
    T3 = 1 - T1 - T2_new
    ternary_df = pd.DataFrame({'T1': T1, 'T2': T2_new, 'T3': T3})
    ternary_df = ternary_df.div(ternary_df.sum(axis=1), axis=0)

    if file_name is None:
        import inspect
        callers_locals = inspect.currentframe().f_back.f_locals
        for var_name, val in callers_locals.items():
            if val is data:
                file_name = var_name
                break
        else:
            file_name = "unnamed"

    if colormap is None:
        colormap = "viridis_r"

    alpha = 0.1
    grid = True
    point_alpha = 0.8
    density_method = "neighbors"
    bandwidth = 0.02

    import matplotlib.pyplot as plt
    
    # If no axis provided, create a new figure
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes()
        create_new_figure = True
    else:
        fig = ax.figure
        create_new_figure = False

    if grid:
        draw_grey_grid_lines(ax, alpha=0.1)

    x_data, y_data = cartizian(ternary_df["T1"], ternary_df["T2"], ternary_df["T3"])

    if density_method == "neighbors":
        points = np.column_stack([x_data, y_data])
        nn = NearestNeighbors(radius=bandwidth)
        nn.fit(points)
        density = nn.radius_neighbors(points, return_distance=False)
        density = np.array([len(neighbors) for neighbors in density])

    scatter = plt.scatter(
        x_data,
        y_data,
        c=density,
        cmap=colormap,
        alpha=point_alpha,
        s=15,
        edgecolors="none",
        zorder=2,
    )

    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)

    # No median line

    label_color = "black"
    label_size = 12
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)
    ax.text(0.51, -0.005, r"$\mathbf{T}_{sum}$", size=label_size, color=label_color)
    ax.text(-0.535, -0.005, r"$\mathbf{T}_{focal}$", size=label_size, color=label_color)

    for spine in ax.spines.values():
        spine.set_color("none")
    ax.set_yticks([])
    ax.set_xticks([])

    # Add colorbar and legend only if creating a new figure
    if create_new_figure:
        sm = plt.cm.ScalarMappable(
            cmap=colormap,
            norm=plt.cm.colors.Normalize(vmin=density.min(), vmax=density.max()),
        )
        sm.set_array([])
        cax = inset_axes(
            ax,
            width="3%",
            height="25%",
            loc="upper right",
            bbox_to_anchor=(0.05, 0.05, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=1,
        )
        cbar = plt.colorbar(sm, cax=cax)
        cbar.ax.set_title("Count", fontsize=10, pad=6)

        legend_text = r"$T_1$ = [{}], $T_2$ = [{}], $T_3$ = all others".format(t1_idx, t_focal_idx)
        fig.text(0.02, 0.98, legend_text, ha='left', va='top', fontsize=12, color='black')

    if show_expectation_legend and create_new_figure:
        fig.text(0.5, 0.12, rf"Expectation: $T_2$ = {expectation:.4f}", color='darkblue', fontsize=13, ha='center', va='top')

    # Only save and return figure if we created a new one
    if create_new_figure:
        title = f"{file_name}_expectation_{expectation:.4f}.png"
        fig.savefig(title, dpi=150, bbox_inches='tight')
        return fig
    else:
        return ax

########################################################
########################################################
# PLOTTING THE PANEL OF PLOTS!!! -- FINAL PLOT
def all_plots(data, t1_idx, colormap="plasma", expectation=None):
    """
    Create a single panel with all ternary projection plots for all possible T2 (focal) topologies
    given a fixed T1 index. Uses the existing plot_ternary_projection function.
    
    Parameters:
    - data: np.ndarray, data to plot
    - t1_idx: int, index of T1 (top vertex)
    - colormap: str, colormap to use
    - expectation: float or None, expectation value for T2 (default None)
    
    Returns:
    - matplotlib figure with subplots arranged in a panel
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    # Get all remaining indices (excluding t1_idx)
    remaining_indices = [i for i in range(data.shape[1]) if i != t1_idx]
    
    # Create subplot layout: 3 rows of 5 plots each (15 total, but we'll use 14)
    fig, axes = plt.subplots(3, 5, figsize=(20, 12))
    axes = axes.flatten()  # Flatten for easier indexing
    
    for i, t2_idx in enumerate(remaining_indices):
        if i < 14:  # Only plot first 14
            ax = axes[i]
            
            # Debug: print what we're working with
            print(f"Plotting T1={t1_idx}, T2={t2_idx}")
            
            # Process data manually (same as in plot_ternary_projection)
            T1 = data[:, t1_idx]
            T2 = data[:, t2_idx]
            mask = [j for j in range(data.shape[1]) if j not in (t1_idx, t2_idx)]
            T3 = data[:, mask].sum(axis=1)
            
            # Apply expectation transformation if provided
            if expectation is not None:
                T2_new = np.zeros_like(T2)
                for j, val in enumerate(T2):
                    if val <= expectation:
                        T2_new[j] = 0.5 * val / expectation if expectation > 0 else 0.0
                    else:
                        T2_new[j] = 0.5 + 0.5 * (val - expectation) / (1 - expectation) if expectation < 1 else 1.0
                T2 = T2_new
            
            T3 = 1 - T1 - T2
            ternary_df = pd.DataFrame({'T1': T1, 'T2': T2, 'T3': T3})
            ternary_df = ternary_df.div(ternary_df.sum(axis=1), axis=0)
            
            # Convert to cartesian coordinates
            x_data, y_data = cartizian(ternary_df["T1"], ternary_df["T2"], ternary_df["T3"])
            
            # Calculate density
            points = np.column_stack([x_data, y_data])
            nn = NearestNeighbors(radius=0.02)
            nn.fit(points)
            density = nn.radius_neighbors(points, return_distance=False)
            density = np.array([len(neighbors) for neighbors in density])
            
            # Add grid lines first
            draw_grey_grid_lines(ax, alpha=0.1)
            
            # Plot on the subplot
            scatter = ax.scatter(
                x_data, y_data,
                c=density,
                cmap=colormap,
                alpha=0.8,
                s=8,
                edgecolors="none"
            )
            
            # Draw triangle
            triangle_x = [0, -0.5, 0.5, 0]
            triangle_y = [h, 0, 0, h]
            ax.plot(triangle_x, triangle_y, color="k", linewidth=0.8, zorder=3)
            
            # No median line
            
            # Labels
            ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=8, color="black")
            ax.text(0.51, -0.005, r"$\mathbf{T}_{sum}$", size=8, color="black")
            ax.text(-0.535, -0.005, r"$\mathbf{T}_{focal}$", size=8, color="black")
            
            # Clean up axes
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            
            # Add title to subplot
            ax.set_title(f"T1={t1_idx}, T2={t2_idx}", fontsize=10)
    
    # Add colorbar to the empty subplot (15th subplot)
    if len(remaining_indices) >= 14:
        cbar_ax = axes[14]  # The 15th subplot (index 14)
        cbar_ax.set_visible(True)
        
        # Create a colorbar using the density values from the last plot
        sm = plt.cm.ScalarMappable(
            cmap=colormap,
            norm=plt.cm.colors.Normalize(vmin=density.min(), vmax=density.max()),
        )
        sm.set_array([])
        
        # Add colorbar to the subplot
        cbar = plt.colorbar(sm, ax=cbar_ax, orientation='vertical')
        cbar.set_label('Density', fontsize=12)
        cbar_ax.set_title('Colorbar', fontsize=12)
        
        # Remove axes from colorbar subplot
        cbar_ax.set_xticks([])
        cbar_ax.set_yticks([])
        for spine in cbar_ax.spines.values():
            spine.set_visible(False)
    
    plt.tight_layout()
    filename = f"all_plots_t1_{t1_idx}"
    if expectation is not None:
        filename += f"_exp_{expectation:.2f}"
    plt.savefig(f"{filename}.png", dpi=150, bbox_inches='tight')
    return fig



print("########################################################")
print("PLOTTING FAKE DATA-- as an example")
print("########################################################")
fig = all_plots(fake_data, t1_idx=0, colormap="plasma")
print("########################################################")
print("please proceed to simulating syntatic data")
print("########################################################")
##########################################################################################
##########################################################################################
# SIMULATING DATA
##########################################################################################
##########################################################################################
""" Now we move on to supply simulated data to our specification instead of a random fake_data
you might need a touch of debegging but this should mostly work.
here are some istructions to what we wish to simulatefrom our meeting 4.9.:

one needs to pick a nice set of parameters in which one can see a difference. Suggestions: Supplements: page 66. Section h, line  I5
 (ne = 500 for each, 100 gen between each split, migration– where to put it..? Maybe P4 → P2 for 10% of the genome and m=0.1. 90 % of loci are simulated without geneflow. And 10# are simulated with the m=0.1. For 9000 loci – m=0. And then for 1000 loci with the same model but m=0.1 P4 --> P2. Then we talk about all those topology weights and stick them together.
 [reminder, we want complex patterns], population model (the “real”history = (1,2),3),4),0) 
In the plots on pages ~50, the loci that have no gene flow are colored grey, and the m != 0 are colored red.
And then, after simulating this, we run the 5 pop visualisation, do a figure where we see all T3=[x], x= 1,2,...15 = focal topology, put it in a panel and see if we can see sth interesting– a panel with these 14.

"""



def simulate_locus_simple(
    populations,
    splits=None,
    migration=None,
    n_loci=10,
    locus_length=1,
    ploidy=1,
    seed=None,
):
    """
    Simulate independent non-recombining loci with msprime

    Args:
        populations: list of dicts with keys:
            - name (str): population name (e.g., 'p0')
            - Ne (float): effective population size
            - sample_size (int): number of haploid samples from this population
            - growth_rate (float, optional): per-generation growth rate (default 0)
        splits: list of dicts, each with either:
            - {time, derived_pop1, derived_pop2, ancestral}
            or
            - {time, derived: [name1, name2], ancestral}
        migration: dict mapping 'source>dest' -> rate (float). Only non-zero rates are applied.
        n_loci: number of independent loci (num_replicates)
        locus_length: length per locus (bp). For pure topology, 1 is fine
        ploidy: 1 for haploid (default), 2 for diploid
        seed: random seed (int). If None, a random seed is drawn

    Returns:
        generator of tskit.TreeSequence objects (one per locus)
    """


    # Build demography
    demography = msprime.Demography()

    # Add populations
    for pop in populations:
        name = pop["name"]
        Ne = pop["Ne"]
        growth_rate = pop.get("growth_rate", 0.0)
        demography.add_population(name=name, initial_size=Ne, growth_rate=growth_rate)

    # Add splits
    if splits:
        for split in splits:
            time = split["time"]
            if "derived" in split:
                derived = split["derived"]
                if len(derived) != 2:
                    raise ValueError("Each split must have exactly two derived populations")
                demography.add_population_split(
                    time=time, derived=derived, ancestral=split["ancestral"]
                )
            else:
                demography.add_population_split(
                    time=time,
                    derived=[split["derived_pop1"], split["derived_pop2"]],
                    ancestral=split["ancestral"],
                )

    # Add migration routes
    if migration:
        for route, rate in migration.items():
            if rate and rate > 0:
                source, dest = route.split(">")
                demography.set_migration_rate(source=source, dest=dest, rate=rate)

    # Seed
    if seed is None:
        seed = random.randint(0, 2**32 - 1)
        print(f"Using random seed: {seed}")

    # Compose samples dict
    samples = {
        pop["name"]: int(pop.get("sample_size", 0)) for pop in populations if int(pop.get("sample_size", 0)) > 0
    }

    ts_iter = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        num_replicates=n_loci,
        sequence_length=locus_length,
        ploidy=ploidy,
        random_seed=seed,
        recombination_rate=0.0,
    )
    return ts_iter


def simulate_chromosome_simple(
    populations,
    splits=None,
    migration=None,
    chromosome_length=1e6,
    rec_rate=1e-8,
    ploidy=1,
    seed=None,
    return_weights=False,
    outgroup=None,
    twisst_verbose=False,
    topology_mapping=None,
    population_labels=None,
):
    """
    Simulate a recombining chromosome with msprime (no external config).

    Args:
        populations: list of dicts with keys:
            - name (str): population name (e.g., 'p0')
            - Ne (float): effective population size
            - sample_size (int): number of haploid samples from this population
            - growth_rate (float, optional): per-generation growth rate (default 0)
        splits: list of dicts, each with either:
            - {time, derived_pop1, derived_pop2, ancestral}
            or
            - {time, derived: [name1, name2], ancestral}
        migration: dict mapping 'source>dest' -> rate (float). Only non-zero rates are applied.
        chromosome_length: total length in bp
        rec_rate: recombination rate per bp per generation
        ploidy: 1 for haploid (default), 2 for diploid
        seed: random seed (int). If None, a random seed is drawn
        return_weights: if True, also compute twisst weights and return (ts, df)
        outgroup: None, int id (e.g., 0), or string ('0' or 'p0'). If None, auto-picked
        twisst_verbose: pass-through to twisst
        topology_mapping: optional mapping for topology reordering
        population_labels: optional labels for logging

    Returns:
        tskit.TreeSequence if return_weights=False
        (tskit.TreeSequence, pandas.DataFrame) if return_weights=True
    """


    # Build demography
    demography = msprime.Demography()

    # Add populations
    for pop in populations:
        name = pop["name"]
        Ne = pop["Ne"]
        growth_rate = pop.get("growth_rate", 0.0)
        demography.add_population(name=name, initial_size=Ne, growth_rate=growth_rate)

    # Add splits
    if splits:
        for split in splits:
            time = split["time"]
            if "derived" in split:
                derived = split["derived"]
                if len(derived) != 2:
                    raise ValueError("Each split must have exactly two derived populations")
                demography.add_population_split(
                    time=time, derived=derived, ancestral=split["ancestral"]
                )
            else:
                demography.add_population_split(
                    time=time,
                    derived=[split["derived_pop1"], split["derived_pop2"]],
                    ancestral=split["ancestral"],
                )

    # Add migration routes
    if migration:
        for route, rate in migration.items():
            if rate and rate > 0:
                source, dest = route.split(">")
                demography.set_migration_rate(source=source, dest=dest, rate=rate)

    # Seed
    if seed is None:
        seed = random.randint(0, 2**32 - 1)
        print(f"Using random seed: {seed}")

    # Compose samples dict
    samples = {
        pop["name"]: int(pop.get("sample_size", 0)) for pop in populations if int(pop.get("sample_size", 0)) > 0
    }

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=chromosome_length,
        recombination_rate=rec_rate,
        ploidy=ploidy,
        random_seed=seed,
    )

    if not return_weights:
        return ts

    # Optional: compute twisst weights
    from twisstntern_simulate.ts_processing import ts_to_twisst_weights

    # Normalize outgroup input to the format expected by ts_to_twisst_weights (string id like '0')
    og = None
    if outgroup is None:
        og = None
    elif isinstance(outgroup, int):
        og = str(outgroup)
    elif isinstance(outgroup, str):
        if outgroup.isdigit():
            og = outgroup
        elif outgroup.startswith("p") and outgroup[1:].isdigit():
            og = outgroup[1:]
        else:
            # Unknown naming; defer to auto-pick
            og = None

    df = ts_to_twisst_weights(
        ts,
        outgroup=og,
        verbose=True,
        twisst_verbose=twisst_verbose,
        topology_mapping=topology_mapping,
        population_labels=population_labels,
    )
    return ts, df

### SIMULATE DATA
pops = [{'name':'p0','Ne':1e4,'sample_size':10}, {'name':'p1','Ne':1e4,'sample_size':10},
        {'name':'p2','Ne':1e4,'sample_size':10}, {'name':'p3','Ne':1e4,'sample_size':10},
        {'name':'p4','Ne':1e4,'sample_size':10}]

# The splits are the population splits that we want to simulate, 
# e.g.
splits = [{'time':500,'derived':['p3','p4'],'ancestral':'p2'},
          {'time':1000,'derived':['p1','p2'],'ancestral':'p0'}]

migration = {
    'p1>p2': 1e-4,  # one-way p1 -> p2
    'p2>p1': 1e-4,  # add this to make it symmetric
    'p3>p4': 5e-5,  # one-way p3 -> p4
    # no other routes will be used
}

          
### UNCOMMENT THIS TO SIMULATE DATA
# i guess we only wanna have loci mode, not chromosome.
# notice that the df for chromome as an extra "Position" column that needs to e erased beforre plotting
ts, df = simulate_chromosome_simple(pops, splits, return_weights=True, outgroup='p0') 

df = df.drop(columns=['Position'])

# now we can plot the data
all_plots(df, 0, colormap="plasma")

# we can also plot the data in a ternary plot
plot_ternary_projection(df, 0, 1, colormap="plasma")
####################
# LOCI MODE -- FINALLY!
routes = [('p1','p2',1e-4), ('p3','p4',5e-5)]  # (source, dest, rate)
migration = {f'{src}>{dst}': rate for src, dst, rate in routes}


ts_iter = simulate_locus_simple(
    populations=pops,
    splits=splits,
    migration=migration,
    n_loci=50,
    locus_length=1,
    ploidy=1,
    seed=42,
)

ts_list = list(ts_iter)

# IMPORTANT: outgroup must match population IDs inside the ts, which are strings '0','1',...
# If unsure, pass outgroup=None (auto-picks the first sampled population).
df = ts_to_twisst_weights(
    ts_iter,
    outgroup=None,         # or '0' for the first pop
    verbose=True,
    twisst_verbose=False,
)
print(df.head())

all_plots(df, 0, colormap="plasma")

plot_ternary_projection(df, 0, 1, colormap="plasma")