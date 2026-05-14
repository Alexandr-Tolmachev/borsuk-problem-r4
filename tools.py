import numpy as np
import itertools
import json
from scipy.spatial import ConvexHull, HalfspaceIntersection
from sklearn.metrics import pairwise_distances
from scipy.optimize import minimize

import torch
import torch.nn as nn
import torch.optim as optim

import math
import argparse

np.math = math



def spheres_params(n):
    r0 = np.sqrt(n / (2 * (n + 1)))
    r1 = 1.0
    a = 0.5 * np.sqrt(2 / (n * (n + 1)))

    c0 = np.zeros(n)
    c1 = np.zeros(n)
    c0[-1] = a
    c1[-1] = a + r0
    return c0, c1, r0, r1


def generate_projections_on_U(points, n, c0, c1, r0, r1):  # U = Lassak cover
    polyhedron_vertices = []
    for point in points:
        d = (point - c0) / np.linalg.norm(point - c0)

        if point[-1] >= 0:
            polyhedron_vertices.append(c0 + d * r0)
        else:
            A = np.sum(d ** 2)  # = 1
            B = 2 * np.dot(c0 - c1, d)
            C = np.sum((c0 - c1) ** 2) - r1 ** 2
            t = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)

            if t < r0:
                polyhedron_vertices.append(c0 + d * t)
            else:
                polyhedron_vertices.append(c0 + d * r0)

    return np.array(polyhedron_vertices)

def generate_points_on_hypercube_surface(n, c0, r0, grid_size=4): # an uniform grid on the surface of the n-dimensional hypercube
    hypercube_vertices = [                                        # centered at the origin, edge length = 2r0   
        np.linspace(c0[i] - r0, c0[i] + r0, num=grid_size)
        for i in range(n)
    ]
    hypercube_vertices_copy = hypercube_vertices.copy()

    points = []
    for i in range(n):
        hypercube_vertices_copy[i] = [c0[i] - r0, c0[i] + r0]
        points.extend(list(itertools.product(*hypercube_vertices_copy)))
        hypercube_vertices_copy[i] = hypercube_vertices[i][1:-1]

    points = np.array(points)
    assert len(points) == grid_size ** n - (grid_size - 2) ** n
    return points


def construct_polyhedron(polyhedron_vertices, n, c0, c1, r0, r1, additional_planes=None):
    all_planes = []

    for v in polyhedron_vertices:
        if v[-1] >= 0:
            normal = (c0 - v) / np.linalg.norm(c0 - v)
        else:
            normal = (c1 - v) / np.linalg.norm(c1 - v)

        offset = -np.dot(normal, v)
        if np.dot(normal, c0) + offset >= 0:
            normal, offset = -normal, -offset

        all_planes.append(np.hstack((normal, offset)))

    all_planes = np.array(all_planes)

    if additional_planes is not None:
        all_planes = np.vstack((all_planes, additional_planes))

    hs = HalfspaceIntersection(all_planes, c0)
    return ConvexHull(hs.intersections)


def get_projection_on_the_polyhedron_all(directions, convex_hull_polyhedron, c0):
    directions = directions / np.linalg.norm(directions, axis=1)[:, None]

    A = convex_hull_polyhedron.equations[:, :-1]
    b = convex_hull_polyhedron.equations[:, -1]

    denom = A @ directions.T
    denom[denom < 1e-12] = 1e-12

    t = -(A @ c0 + b)[:, None] / denom
    t_min = np.min(t, axis=0)

    curr_points = c0 + t_min[:, None] * directions
    assert np.all(abs((curr_points @ A.T + b).max(axis=1)) < 1e-9), 'wrong projection'

    return curr_points
    

def select_polyhedron_vertices(curr_convex_hull, convex_hull_polyhedron, c0):
    c0_index = np.argmin(((curr_convex_hull.points - c0) ** 2).sum(axis=1))
    assert c0_index == 0

    facet_mask = [c0_index in f for f in curr_convex_hull.simplices]
    curr_equations = curr_convex_hull.equations[facet_mask]

    pts = convex_hull_polyhedron.points
    pts1 = np.hstack([pts, np.ones((len(pts), 1))])

    mask = np.all(pts1 @ curr_equations.T <= 1e-9, axis=1)
    return pts[mask], mask


def find_polyhedron_diameter(vertices):
    d = pairwise_distances(vertices)
    return np.max(d), np.argmax(d)



# def compute_diameters_from_facets(vertices, facets, convex_hull_polyhedron, c0):
#     polyhedron_vertex_masks = np.zeros((len(facets), len(convex_hull_polyhedron.points)))
#     diameters = np.zeros(len(facets))

#     for i, facet in enumerate(facets):
#             curr_cone_vertices = vertices[facet]
#         #try:
#             curr_hull = ConvexHull(np.vstack((c0, curr_cone_vertices)))
#             poly_v, mask = select_polyhedron_vertices(
#                 curr_hull, convex_hull_polyhedron, c0
#             )
#             polyhedron_vertex_masks[i] = mask


#             center_index = np.argmin(((curr_hull.points - c0) ** 2).sum())
#             assert center_index == 0, 'ConvexHull permutes points!'

#             # выделим плоскости, которые проходят через с0 (границы конуса)
#             facet_mask = [center_index in facet for facet in curr_hull.simplices]
#             curr_equations = curr_hull.equations[facet_mask]

#             # ищем пересечение границ конуса с плоскостями, задающими описанный вокруг покрышки многогранник
#             interior_point = np.mean(curr_hull.points, axis=0)  # точка, строго внутри этого многранника

#             # cчитаем диаметр части разбиения как диаметр описанного вокруг нее многогранника
#             hs = HalfspaceIntersection(np.vstack((convex_hull_polyhedron.equations, curr_equations)), interior_point)
#             part_vertices = ConvexHull(hs.intersections).points

#             diameters[i], ij = find_polyhedron_diameter(part_vertices)

#         #except Exception:
#         #    diameters[i] = 2.0
        

#     #assert np.all(polyhedron_vertex_masks.sum(axis=0) >= 1), 'Not all polyhedron vertices have been covered!'
#     return diameters


def find_polyhedron_hulls_from_facets(vertices, facets, convex_hull_polyhedron, c0):
    polyhedron_vertex_masks = np.zeros((len(facets), len(convex_hull_polyhedron.points)))
    diameters = np.zeros(len(facets))

    polyhedron_hulls = []

    for i, facet in enumerate(facets):
            curr_cone_vertices = vertices[facet]
        #try:
            curr_hull = ConvexHull(np.vstack((c0, curr_cone_vertices)))
            poly_v, mask = select_polyhedron_vertices(
                curr_hull, convex_hull_polyhedron, c0
            )
            polyhedron_vertex_masks[i] = mask


            center_index = np.argmin(((curr_hull.points - c0) ** 2).sum())
            assert center_index == 0, 'ConvexHull permutes points!'

            # выделим плоскости, которые проходят через с0 (границы конуса)
            facet_mask = [center_index in facet for facet in curr_hull.simplices]
            curr_equations = curr_hull.equations[facet_mask]

            # ищем пересечение границ конуса с плоскостями, задающими описанный вокруг покрышки многогранник
            interior_point = np.mean(curr_hull.points, axis=0)  # точка, строго внутри этого многранника

            # cчитаем диаметр части разбиения как диаметр описанного вокруг нее многогранника
            hs = HalfspaceIntersection(np.vstack((convex_hull_polyhedron.equations, curr_equations)), interior_point)

            polyhedron_hulls.append(ConvexHull(hs.intersections))

        #except Exception:
        #    diameters[i] = 2.0


    #assert np.all(polyhedron_vertex_masks.sum(axis=0) >= 1), 'Not all polyhedron vertices have been covered!'
    return polyhedron_hulls


def compute_diameters_from_facets(vertices, facets, convex_hull_polyhedron, center0):
    polyhedron_hulls = find_polyhedron_hulls_from_facets(vertices, facets, convex_hull_polyhedron, center0)
    diameters = np.zeros(len(facets))

    for i, curr_part in enumerate(polyhedron_hulls):
        diameters[i], ij = find_polyhedron_diameter(curr_part.points)
        
    return diameters

    
# ============================================================
#  TORCH approximate OBJECTIVE
# ============================================================


def np_project_on_polyhedron(directions, center, A, b):
    denom = A @ directions.T    
    numer = (A @ center + b)[:, None]             # (n_planes, 1)
    denom = np.clip(denom, a_min=1e-15, a_max = 1e9)
    t = -numer / denom      
    t_min = np.min(t, axis = 0)
    curr_points = center + t_min[:, None] * directions
    return curr_points

def torch_project_on_polyhedron(directions, center, A, b):
    directions = directions / torch.norm(directions, dim=1, keepdim=True)
    denom = A @ directions.T
    denom = torch.clamp(denom, min=1e-15)

    numer = (A @ center + b)[:, None]             # (n_planes, 1)
    t = -numer / denom                            # (n_planes, n_dirs)
    t_min, _ = t.min(dim=0)

    curr_points = center + t_min[:, None] * directions

    return curr_points


def approximate_diameters(vertices, center, facets): # silmple lower estimate, only base vertices are used
    d_all = []
    for f in facets:
        pts = torch.cat([vertices[f], center[None, :]], dim=0)
        d = torch.cdist(pts, pts).max()
        d_all.append(d)
    return torch.stack(d_all)


def ray_trace_two_spheres_torch(points, center, n, c0, c1, r0, r1, EQ, EQ_b):

    directions = points.clone()
    
    norms = torch.norm(directions, dim=1, keepdim=True) # normalization
    epsilon = 1e-10
    d = directions / (norms + epsilon)
    
    # quadratic equation, sphere (c0, r0)
    A = torch.ones(points.shape[0], device=points.device, dtype=torch.double)
    B = 2 * torch.sum((center - c0) * d, dim=1)
    C = torch.sum((center - c0) ** 2) - r0 ** 2
    
    discriminant = B * B - 4 * A * C
    sqrt_disc = torch.sqrt(discriminant)
    t_sph0 = (-B + sqrt_disc) / 2
    
    # quadratic equation, sphere (c1, r1)
    A = torch.ones(points.shape[0], device=points.device)  
    B = 2 * torch.sum((center - c1) * d, dim=1)
    C = torch.sum((center - c1) ** 2) - r1 ** 2
    
    discriminant = B * B - 4 * A * C
    sqrt_disc = torch.sqrt(discriminant)
    t_sph1 = (-B + sqrt_disc) / 2

    # hyperplanes
    denom = EQ @ d.T
    denom = torch.clamp(denom, min=epsilon)

    numer = (EQ @ center + EQ_b)[:, None]             # (n_planes, 1)
    t = -numer / denom                            # (n_planes, n_dirs)
    t_planes, _ = t.min(dim=0)
    
    
    distances, _ = torch.min(torch.stack((t_sph0, t_sph1, t_planes)), 0)
    
    return center + d * distances.unsqueeze(-1)

def approximate_diameters_intermediate(vertices, center, facets, c0, c1, A, b, additional_points2d = 1, additional_grid3d = 5):
    d_all = []
    diam_vert = []
    denom = additional_points2d + 1
    gg = additional_grid3d
    for f in facets:  # convex linear  	
        addpts = []
        for i in range(len(f)-1):
            for j in range(i+1,len(f)):
                for k in range(additional_points2d):
                    p = (k+1)/denom * vertices[f[i]] + (denom-k-1)/denom * vertices[f[j]] #- center
                    addpts.append(p)
        if additional_grid3d:                    
            for i in range(len(f)-2):
                for j in range(i+1,len(f)-1):            
                    for k in range(j+1,len(f)):
                        for aa in range(1, gg-1):
                            for bb in range(1, gg-aa):
                                cc = gg-aa-bb
                                p =  (aa*vertices[f[i]] + bb*vertices[f[j]]/3 + cc*vertices[f[k]])/gg
                                addpts.append(p)

        pts1 = torch.stack(addpts, dim = 0)
        pts12 = torch.cat([vertices[f], pts1], dim=0)
        pts3 = ray_trace_two_spheres_torch(pts12, center, n, c0, c1, r0, r1, A, b)
        pts = torch.cat([pts3,  center[None, :]], dim=0)                    	                    	                    	                    	
        
        x2 = torch.square(pts)
        x2s = torch.sum(x2,1)
        distm =  - 2* pts.mm(pts.t()) + x2s + x2s.reshape((-1,1))
        d = torch.max(distm)
        d_all.append(d)
    return torch.stack(d_all)


def validity_approximate_loss(directions, eps=1e-3):
    # log-barrier for origin inside convex hull of directions
    norms = torch.norm(directions, dim=1)
    loss = torch.relu(eps - norms).sum()
    return loss
    
    
def fun(x): # wrapper for a black-box optimization    
    global convex_hull_pylyhedron, facets
    x = np.array(x)
    mn = x.shape[0]
    m = mn // 4
    xx = x.reshape((m,4))
    directions = xx[:-1,:]
    center = xx[-1,:]
    A = convex_hull_polyhedron.equations[:, :-1]
    b = convex_hull_polyhedron.equations[:, -1]        
    verts = np_project_on_polyhedron(directions, center, A, b)                    
    diam_exact, _ = compute_diameters_from_facets(
        verts, facets, convex_hull_polyhedron, center
    )    
    return max(diam_exact)

