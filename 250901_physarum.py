import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# ------------------------------
# Parameters
# ------------------------------
mu = 0.1             # decay coefficient
delta_t = 0.1        # time step
iterations = 100     # number of iterations
I0 = 1.0             # unit flow

# ------------------------------
# Read edges from CSV
# ------------------------------
# CSV format: source,target,length,D_max
edge_df = pd.read_csv("edges.csv")

# Min-Max normalization of edge lengths
L_min = edge_df["length"].min()
L_max = edge_df["length"].max()
if L_max > L_min:
    edge_df["length_norm"] = (edge_df["length"] - L_min) / (L_max - L_min)
else:
    edge_df["length_norm"] = 1.0  # all lengths are the same

edges = list(edge_df.itertuples(index=False, name=None))

# Define start and end nodes manually
start = 0
end = 63

# ------------------------------
# Create the graph
# ------------------------------
G = nx.Graph()
for u, v, l, dmax, l_norm in edges:
    G.add_edge(u, v,
               length=l,           # original length
               length_norm=l_norm, # normalized length
               D=1.0,
               D_max=dmax,
               Q=0.0)

nodes = list(G.nodes)
node_index = {node: i for i, node in enumerate(nodes)}
n = len(nodes)

# ------------------------------
# Main loop
# ------------------------------
for step in range(iterations):
    # --- Build matrix A and vector b for linear system Ap = b ---
    A = lil_matrix((n, n))
    b = np.zeros(n)

    for u in G.nodes:
        i = node_index[u]
        for v in G.neighbors(u):
            j = node_index[v]
            edge = G[u][v]
            D = edge['D']
            L = edge['length']
            c = D / L
            A[i, i] += c
            A[i, j] -= c

    # Boundary conditions (source and sink)
    b[node_index[start]] = I0
    b[node_index[end]] = -I0

    # --- Fix reference node potential to 0 to avoid singular matrix ---
    ref_node = n - 1
    A[ref_node, :] = 0
    A[:, ref_node] = 0
    A[ref_node, ref_node] = 1
    b[ref_node] = 0

    # Solve for potentials
    p = spsolve(A.tocsr(), b)

    # --- Compute flow Q and update conductivity D ---
    for u, v in G.edges:
        i = node_index[u]
        j = node_index[v]
        edge = G[u][v]
        L = edge['length']
        D = edge['D']
        Q = (D / L) * (p[i] - p[j])
        edge['Q'] = Q
        new_D = D + delta_t * (abs(Q) - mu * D)
        edge['D'] = max(min(new_D, edge['D_max']), 1e-6)

# ------------------------------
# Save edge results to CSV
# ------------------------------
output_data = []
for u, v in G.edges:
    edge = G[u][v]
    output_data.append({
        "source": u,
        "target": v,
        "length": edge["length"],
        "length_norm": edge["length_norm"],
        "D_final": edge["D"],
        "D_max": edge["D_max"],
        "Q_final": edge["Q"]
    })

output_df = pd.DataFrame(output_data)
output_df.to_csv("physarum_results.csv", index=False)
print("Saved results to 'physarum_results.csv'")

# ------------------------------
# Extract path by following max-Q edges
# ------------------------------
path_edges = []
visited = set()
current = start

while current != end and current not in visited:
    visited.add(current)
    neighbors = [(nbr, G[current][nbr]['Q']) for nbr in G.neighbors(current) if nbr not in visited]
    if not neighbors:
        print("Could not reach the goal node.")
        break
    # pick neighbor with maximum |Q|
    next_node, _ = max(neighbors, key=lambda x: abs(x[1]))
    path_edges.append((current, next_node, G[current][next_node]['Q']))
    current = next_node

# ------------------------------
# Save path to CSV
# ------------------------------
if current == end:
    path_data = []
    for u, v, Q in path_edges:
        edge = G[u][v]
        path_data.append({
            "source": u,
            "target": v,
            "length": edge["length"],
            "length_norm": edge["length_norm"],
            "D_final": edge["D"],
            "D_max": edge["D_max"],
            "Q_final": Q
        })
    path_df = pd.DataFrame(path_data)
    path_df.to_csv("physarum_path.csv", index=False)
    print("Saved max-Q path to 'physarum_path.csv'")

# 4x4x4 cube positions
pos_3d = {i*16 + j*4 + k: (i, j, k) for i in range(4) for j in range(4) for k in range(4)}

def draw_3d_network(G, pos_3d, edge_attr, filename, cmap=plt.cm.plasma, highlight_edges=None, node_color="skyblue", node_size=50):
    """
    Draw a 3D network and save as PNG.
    
    Parameters:
    -----------
    G : networkx.Graph
        The graph containing nodes and edges.
    pos_3d : dict
        Mapping {node: (x,y,z)} for 3D positions.
    edge_attr : str
        Edge attribute to color edges by (e.g., 'D' or 'D_max').
    filename : str
        Output PNG file path.
    cmap : matplotlib colormap
        Colormap to use for edges.
    highlight_edges : list of (u,v) tuples
        Edges to highlight in red (optional).
    node_color : str
        Node color.
    node_size : int
        Node marker size.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection="3d")
    
    # Prepare edge colors for colormap if needed
    if highlight_edges is None:
        edge_values = [G[u][v][edge_attr] for u, v in G.edges]
        norm = plt.Normalize(vmin=min(edge_values), vmax=max(edge_values))
    else:
        norm = None  # not used for highlighted edges

    # Draw edges
    for u, v in G.edges:
        if highlight_edges is not None and ((u, v) in highlight_edges or (v, u) in highlight_edges):
            color = "red"
            linewidth = 3
        else:
            if highlight_edges is None:
                color = cmap(norm(G[u][v][edge_attr]))  # <-- float -> RGBA
            else:
                color = "lightgray"
            linewidth = 2

        x = [pos_3d[u][0], pos_3d[v][0]]
        y = [pos_3d[u][1], pos_3d[v][1]]
        z = [pos_3d[u][2], pos_3d[v][2]]
        ax.plot(x, y, z, color=color, linewidth=linewidth)
    
    # Draw nodes
    xs, ys, zs = zip(*[pos_3d[n] for n in G.nodes])
    ax.scatter(xs, ys, zs, c=node_color, s=node_size, edgecolors="k")
    
    # Draw labels
    for n, (x, y, z) in pos_3d.items():
        ax.text(x, y, z, str(n), fontsize=8)
    
    # Colorbar for non-highlighted figures
    if highlight_edges is None:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.1)
        cbar.set_label(edge_attr)
    
    # Remove axes
    ax.set_axis_off()
    
    # Save figure
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved 3D visualization to '{filename}'")

# ------------------------------
# 1) Full network colored by D_final
# ------------------------------
draw_3d_network(G, pos_3d, edge_attr="D", filename="physarum_network.png")

# ------------------------------
# 2) Max-Q path only (highlight in red)
# ------------------------------
# Prepare highlight_edges as list of (u,v)
highlight_edges = [(int(row["source"]), int(row["target"])) for idx, row in path_df.iterrows()]
draw_3d_network(G, pos_3d, edge_attr="D", filename="physarum_path.png", highlight_edges=highlight_edges, node_color="lightgray")

# ------------------------------
# 3) Full network colored by D_max
# ------------------------------
draw_3d_network(G, pos_3d, edge_attr="D_max", filename="physarum_dmax.png")