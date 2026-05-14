import numpy as np
import math
import argparse

np.math = math

from tools import *
from argparse import ArgumentParser   
import sys
    
if __name__ == "__main__":
    
    parser = ArgumentParser(
        prog = 'valdation.py',
        description = 'Check that a partition of a truncated Lassak cover in R4 into 8 parts of a smaller diameter is correct.'
    )    
    
    n = 4
    c0, c1, r0, r1 = spheres_params(n)

    parser.add_argument('filepath')        
    parser.add_argument('-g', '--gridsize')   
    
    args = parser.parse_args()

    try:
        with open(args.filepath, 'r') as f:
            ff = f.read()    
        print(f'CHECKING THE PARTITION <{args.filepath}> OF THE TRUNCATED LASSAK COVER...')     
        exec(ff)                    
        print(f'hyperplanes = {additional_planes}')            
        print(f'directions = {init_dirs}')        
        print(f'center = {center0}')
        print(f'grid size = {int(args.gridsize)}')
    except:
        print('Requires a description of the polyhedral cones and a grid size. Usage: validation.py <filepath> -g <gridsize>.')
        sys.exit()

    facets = np.array([
        [0,1,2,3,4,5,6,7],
        [8,9,10,11,12,13,14,15],
        [0,1,2,3,8,9,10,11],
        [4,5,6,7,12,13,14,15],
        [0,1,4,5,8,9,12,13],
        [2,3,6,7,10,11,14,15],
        [0,2,4,6,8,10,12,14],
        [1,3,5,7,9,11,13,15],
    ])  	    
            
    hypercube_points = generate_points_on_hypercube_surface(n, c0, r0, grid_size=int(args.gridsize))

    polyhedron_vertices = generate_projections_on_U(
        hypercube_points, n, c0, c1, r0, r1
    )

    
		
    convex_hull_polyhedron = construct_polyhedron(
        polyhedron_vertices, n, c0, c1, r0, r1, additional_planes
    )
    
    verts_np = get_projection_on_the_polyhedron_all(init_dirs, convex_hull_polyhedron, center0)
          
    diam = compute_diameters_from_facets(verts_np, facets, convex_hull_polyhedron, center0)
    
    print('Upper estimates of diameters:')
    for i,d in enumerate(diam):
        print(f'part {i}: {d}')
    print(f'Maximal diameter: {diam.max()}')
    if diam.max() < 1.0:
        print('---CORRECT---')
    else:    
        print('---PROBABLY ERRONEOUS---')

