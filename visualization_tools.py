import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import imageio


# ---------------------------
# hyperplane basis
# ---------------------------
def plane_basis(n):
    n = n / np.linalg.norm(n)
    U, _, _ = np.linalg.svd(np.eye(4) - np.outer(n, n))
    return U[:, :3]


# ---------------------------
# edges
# ---------------------------
def extract_edges(hull):
    edges = set()
    for simplex in hull.simplices:
        for i in range(len(simplex)):
            for j in range(i + 1, len(simplex)):
                edges.add(tuple(sorted((simplex[i], simplex[j]))))
    return list(edges)


# ---------------------------
# intersections
# ---------------------------
def intersect_polyhedron(vertices, edges, n, d):
    pts = []

    for i, j in edges:
        v1, v2 = vertices[i], vertices[j]

        denom = np.dot(n, v2 - v1)
        if abs(denom) < 1e-10:
            continue

        t = -(np.dot(n, v1) + d) / denom

        if 0 <= t <= 1:
            pts.append(v1 + t * (v2 - v1))

    if len(pts) < 4:
        return None

    return np.array(pts)


# ---------------------------
# plotting pictures for GIF images
# ---------------------------
def plot_section(ax, pts3d, color):
    if pts3d is None or len(pts3d) < 4:
        return

    hull = ConvexHull(pts3d)

    faces = [pts3d[s] for s in hull.simplices]

    poly = Poly3DCollection(
        faces,
        alpha=0.2,                
        facecolor=color,
        edgecolor='black',        
        linewidth=0.1
    )

    poly.set_zsort('average')  

    ax.add_collection3d(poly)


# ---------------------------
def make_4d_sections_gif(
    hulls,
    normal_vector,
    center0,
    gif_path="sections.gif",
    gif_duration=6,
    num_frames=120,
    show=False
):
    assert len(hulls) == 8

    d_values = np.linspace(-1, 1, num_frames)

    d_values -= (normal_vector * center0).sum()

    colors = [
        "#d62728",
        "#1f77b4",
        "#2ca02c",
        "#9467bd",
        "#ff7f0e",
        "#bcbd22",
        "#8c564b",
        "#e377c2",
    ]
    basis = plane_basis(normal_vector)

    frames = []

    for frame_idx, d in enumerate(d_values):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        empty = True

        for hull, color in zip(hulls, colors):
            vertices = hull.points
            edges = extract_edges(hull)

            pts4d = intersect_polyhedron(vertices, edges, normal_vector, d)

            if pts4d is not None:
                empty = False
                pts3d = pts4d @ basis
                plot_section(ax, pts3d, color)

        if not empty:
            # --- ВАЖНО: геометрия вида ---
            ax.set_box_aspect([1, 1, 1])

            ax.set_axis_off()

            # --- ВРАЩЕНИЕ (даёт ощущение 3D) ---
            ax.view_init(elev=20, azim=frame_idx * 3)

            if empty:
                ax.set_facecolor('white')


            # --- buffer_rgba ---
            fig.canvas.draw()
            buf = np.asarray(fig.canvas.buffer_rgba())
            image = buf[:, :, :3]

            frames.append(image)
            plt.close(fig)

  
    # save gif
    imageio.mimsave(gif_path, frames, fps=len(frames)/gif_duration, loop=0)

    if show:
        from IPython.display import Image, display
        display(Image(filename=gif_path))

    return gif_path
