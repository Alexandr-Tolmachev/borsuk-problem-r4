import numpy as np
import math
import argparse

np.math = math

from tools import *
from visualization_tools import *

from argparse import ArgumentParser   
import sys
    
if __name__ == "__main__":
    
    parser = ArgumentParser(
        prog = 'visualization.py',
        description = 'Visualize the projections of 8 parts into moving plane parameterized by its normal vector'
    )    
    
    n = 4
    c0, c1, r0, r1 = spheres_params(n)

    parser.add_argument('filename')        
    parser.add_argument('-g', '--gridsize', default='4')  # default grid parameter
    parser.add_argument('-normal', default='[1,1,1,1]')    # default normal vector
    parser.add_argument('-gif_duration', default='6')      # default value (seconds)
    parser.add_argument('-gif_path', default='animation.gif')  # default path to save the animation 
    parser.add_argument('-num_frames', default='120')      # default num_frames value 
    
    args = parser.parse_args()

    try:
        with open(args.filename, 'r') as f:
            ff = f.read()    
        print(f'CHECKING THE PARTITION <{args.filename}> OF THE TRUNCATED LASSAK COVER...')     
        exec(ff)                    
        print(f'hyperplanes = {additional_planes}')            
        print(f'directions = {init_dirs}')        
        print(f'center = {center0}')
        print(f'grid size = {int(args.gridsize)}')
        print(f'gif duration = {int(args.gif_duration)} sec')
        print(f'num frames = {int(args.num_frames)}')
        print(f'gif path = {args.gif_path}')
    except Exception as e:
        print(f'Requires a description of GIF. Usage: visualization.py <filename> -g <gridsize> -normal [1, 1, 1, 1]. Err or: {e}')
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
    
    # find parts
    hulls = find_polyhedron_hulls_from_facets(verts_np, facets, convex_hull_polyhedron, center0)

    # transform normal_vector from str to list
    normal_vector_str = args.normal
    if normal_vector_str.startswith('[') and normal_vector_str.endswith(']'):
        normal_vector = eval(normal_vector_str)
    else:
        normal_vector = [float(x) for x in normal_vector_str.split(',')]
    
    # plot GIF
    make_4d_sections_gif(
        hulls,
        normal_vector,
        center0=center0,
        gif_path=str(args.gif_path),
        gif_duration=int(args.gif_duration),
        num_frames=int(args.num_frames),
		show=False
    )
