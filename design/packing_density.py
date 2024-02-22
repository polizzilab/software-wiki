from scipy.spatial import Delaunay
import prody as pr
import numpy as np
from numba import jit, prange
import copy
from numpy import dot
from math import sqrt
from numba import jit
from concurrent.futures import ThreadPoolExecutor
import argparse
import os

"""
Example command line usage:

> pdb=/Users/nickpolizzi/projects/Combs2/designs/apx/paper_designs/LABLE_design.pdb
> outdir=/Users/nickpolizzi/projects/Combs2/designs/scripts
> filename='packdensity.txt'
> python packing_density.py --p $pdb -o $outdir -f $filename

Example usage in a script:
from packing_density import calc_packing_density

pack_density = calc_packing_density(pdb, n_workers=1)
"""

def is_inside_spheres_no_thread(points, centers, radii, chunk_size=100000):
    n_points = len(points)
    inside_flags = np.zeros(n_points, dtype=bool)
    
    for i in range(0, n_points, chunk_size):
        end = min(i + chunk_size, n_points)
        chunk = points[i:end]
        
        # Expand dimensions for broadcasting, but only for the chunk
        chunk_exp = np.expand_dims(chunk, axis=1)
        centers_exp = np.expand_dims(centers, axis=0)
        
        # Calculate distances from centers to points in the chunk
        distances = np.linalg.norm(chunk_exp - centers_exp, axis=2)
        
        # Identify points that are inside any of the spheres
        inside_chunk = np.any(distances < radii, axis=1)
        
        # Store the results
        inside_flags[i:end] = inside_chunk
    
    return inside_flags

def chunkwise_processing(chunk, centers, radii):
    chunk_exp = np.expand_dims(chunk, axis=1)
    centers_exp = np.expand_dims(centers, axis=0)
    distances = np.linalg.norm(chunk_exp - centers_exp, axis=2)
    inside_chunk = np.any(distances < radii, axis=1)
    return inside_chunk

def is_inside_spheres(points, centers, radii, chunk_size=100000, n_workers=8):
    n_points = len(points)
    inside_flags = np.zeros(n_points, dtype=bool)
    
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for i in range(0, n_points, chunk_size):
            end = min(i + chunk_size, n_points)
            chunk = points[i:end]
            futures.append(executor.submit(chunkwise_processing, chunk, centers, radii))

        for i, future in enumerate(futures):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, n_points)
            inside_flags[start:end] = future.result()

    return inside_flags


def set_elements_sel(pdb):
    sel = pdb.select('protein')
    names = sel.getNames()
    elems = []
    for name in names:
        if name[0].isalpha():
            elems.append(name[0].upper())
        elif name[0].isdigit() and name[1].isalpha():
            elems.append(name[1].upper())
        else:
            elems.append('X')
    sel.setElements(elems)


def calc_packing_density(path_to_pdb: str, 
                        #  outdir='/.', 
                        #  filename='packingdensity.txt',
                        #  overwrite=False,
                         selection: str = None, 
                         n_points: int = 1000000, 
                         chunk_size: int = 100000, 
                         n_workers: int = 4):
    '''
    Calculate the packing density of a protein structure.
    Inputs:
        path_to_pdb: str, path to the PDB file
        selection: str, prody selection string for the atoms to include in the calculation
        n_points: int, number of random points to use in the Monte Carlo simulation
        chunk_size: int, number of points to process at a time
        n_workers: int, number of threads to use for parallel processing
    '''
    
    pdb = pr.parsePDB(path_to_pdb)
    if selection is not None:
        pdb = pdb.select(selection)
    ahull = AlphaHull()
    ahull.set_coords(pdb)
    ahull.calc_hull()
    vol = ahull.get_volume()
    
    sel = pdb.select('protein')
    if (sel.getElements() == '').all():
        set_elements_sel(pdb)

    ex, inter, bur = partition_res_by_burial(sel)
    resinds = " ".join(str(s) for s in bur | inter)
    sel_C = pdb.select(f'sidechain and protein and element C and resindex {resinds}')
    sel_S = pdb.select(f'sidechain and protein and element S and resindex {resinds}')
    sel_O = pdb.select(f'sidechain and protein and element O and resindex {resinds}')
    sel_N = pdb.select(f'sidechain and protein and element N and resindex {resinds}')

    protein_coords = sel.select('backbone or name CB').getCoords()
    bbox_min = np.min(protein_coords, axis=0) 
    bbox_max = np.max(protein_coords, axis=0)

    # print("Bounding box min coordinates:", bbox_min)
    # print("Bounding box max coordinates:", bbox_max)
    centers_C, centers_S, centers_O, centers_N = [np.array([]).reshape(0,3)] * 4
    radii_C, radii_S, radii_O, radii_N = [np.array([])] * 4
    if sel_C is not None:
        centers_C = sel_C.getCoords()
        radii_C = 1.7*np.ones(len(centers_C))
    if sel_S is not None:
        centers_S = sel_S.getCoords()
        radii_S = 1.8*np.ones(len(centers_S))
    if sel_N is not None:
        centers_N = sel_N.getCoords()
        radii_N = 1.5*np.ones(len(centers_N))
    if sel_O is not None:
        centers_O = sel_O.getCoords()
        radii_O = 1.4*np.ones(len(centers_O))
    if sel_C is None and sel_S is None and sel_O is None and sel_N is None:
        raise ValueError("No C/S or O/N sidechain atoms found.")

    centers = np.vstack((centers_C, centers_S, centers_N, centers_O))
    radii = np.hstack((radii_C, radii_S, radii_N, radii_O))

    # # Monte Carlo simulation
    # points = np.random.rand(n_points, 3) * (bbox_max - bbox_min) + bbox_min
    # # Check which points are inside any of the spheres
    # inside_flags = is_inside_spheres(points, centers, radii, chunk_size=chunk_size, n_workers=n_workers)
    # # Count how many points are inside
    # inside_count = np.sum(inside_flags)
    # # Volume of the bounding box
    # bbox_volume = np.prod(bbox_max - bbox_min)
    # # Estimate the collective volume of the spheres
    # estimated_volume = (inside_count / n_points) * bbox_volume
    # # print("Estimated collective volume:", estimated_volume)

    # Monte Carlo simulation
    n_points_orig = n_points
    points = np.random.rand(n_points, 3) * (bbox_max - bbox_min) + bbox_min
    # print("Number of points:", n_points_orig)
    points = points[ahull.pnts_in_hull(points)]
    if points.shape[0] < n_points_orig / 2:
        new_points = np.random.rand(n_points_orig, 3) * (bbox_max - bbox_min) + bbox_min
        new_points = points[ahull.pnts_in_hull(points)]
        points = np.vstack((points, new_points))
    n_points = points.shape[0]
    # print("Number of points inside hull:", n_points)

    # Check which points are inside any of the spheres
    inside_flags = is_inside_spheres(points, centers, radii, chunk_size=chunk_size, n_workers=n_workers)
    # Count how many points are inside
    inside_count = np.sum(inside_flags)
    # Estimate the collective volume of the spheres
    estimated_volume = (inside_count / n_points) * vol

    # flag = 'a'
    # if overwrite:
    #     print(f'Overwriting {filename}')
    #     flag = 'w'

    # with open(os.path.join(outdir, filename), flag) as f:
    #     f.write(f'{path_to_pdb} {estimated_volume / vol} \n')

    return estimated_volume / vol


@jit("f8(f8[:], f8[:,:])", nopython=True, cache=True)
def distance_(P, TRI):
    # function [dist,PP0] = pointTriangleDistance(TRI,P)
    # calculate distance between a point and a triangle in 3D
    # SYNTAX
    #   dist = pointTriangleDistance(TRI,P)
    #   [dist,PP0] = pointTriangleDistance(TRI,P)
    #
    # DESCRIPTION
    #   Calculate the distance of a given point P from a triangle TRI.
    #   Point P is a row vector of the form 1x3. The triangle is a matrix
    #   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
    #   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
    #   to the triangle TRI.
    #   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
    #   closest point PP0 to P on the triangle TRI.
    #
    # Author: Gwolyn Fischer
    # Release: 1.0
    # Release date: 09/02/02
    # Release: 1.1 Fixed Bug because of normalization
    # Release: 1.2 Fixed Bug because of typo in region 5 20101013
    # Release: 1.3 Fixed Bug because of typo in region 2 20101014

    # Possible extention could be a version tailored not to return the distance
    # and additionally the closest point, but instead return only the closest
    # point. Could lead to a small speed gain.

    # Example:
    # %% The Problem
    # P0 = [0.5 -0.3 0.5]
    #
    # P1 = [0 -1 0]
    # P2 = [1  0 0]
    # P3 = [0  0 0]
    #
    # vertices = [P1; P2; P3]
    # faces = [1 2 3]
    #
    # %% The Engine
    # [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0)
    #
    # %% Visualization
    # [x,y,z] = sphere(20)
    # x = dist*x+P0(1)
    # y = dist*y+P0(2)
    # z = dist*z+P0(3)
    #
    # figure
    # hold all
    # patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8)
    # plot3(P0(1),P0(2),P0(3),'b*')
    # plot3(PP0(1),PP0(2),PP0(3),'*g')
    # surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
    # view(3)

    # The algorithm is based on
    # "David Eberly, 'Distance Between Point and Triangle in 3D',
    # Geometric Tools, LLC, (1999)"
    # http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    #
    #        ^t
    #  \     |
    #   \reg2|
    #    \   |
    #     \  |
    #      \ |
    #       \|
    #        *P2
    #        |\
    #        | \
    #  reg3  |  \ reg1
    #        |   \
    #        |reg0\
    #        |     \
    #        |      \ P1
    # -------*-------*------->s
    #        |P0      \
    #  reg4  | reg5    \ reg6
    # rewrite triangle in normal form
    B = TRI[0, :]
    E0 = TRI[1, :] - B
    # E0 = E0/sqrt(sum(E0.^2)); %normalize vector
    E1 = TRI[2, :] - B
    # E1 = E1/sqrt(sum(E1.^2)); %normalize vector
    D = B - P
    a = dot(E0, E0)
    b = dot(E0, E1)
    c = dot(E1, E1)
    d = dot(E0, D)
    e = dot(E1, D)
    f = dot(D, D)

    #print "{0} {1} {2} ".format(B,E1,E0)
    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e

    # Terible tree of conditionals to determine in which region of the diagram
    # shown above the projection of the point into the triangle-plane lies.
    if (s + t) <= det:
        if s < 0.0:
            if t < 0.0:
                # region4
                if d < 0:
                    t = 0.0
                    if -d >= a:
                        s = 1.0
                        sqrdistance = a + 2.0 * d + f
                    else:
                        s = -d / a
                        sqrdistance = d * s + f
                else:
                    s = 0.0
                    if e >= 0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        if -e >= c:
                            t = 1.0
                            sqrdistance = c + 2.0 * e + f
                        else:
                            t = -e / c
                            sqrdistance = e * t + f

                            # of region 4
            else:
                ## region 3
                s = 0
                if e >= 0:
                    t = 0
                    sqrdistance = f
                else:
                    if -e >= c:
                        t = 1
                        sqrdistance = c + 2.0 * e + f
                    else:
                        t = -e / c
                        sqrdistance = e * t + f
                        # of region 3
        else:
            if t < 0:
                ## region 5
                t = 0
                if d >= 0:
                    s = 0
                    sqrdistance = f
                else:
                    if -d >= a:
                        s = 1
                        sqrdistance = a + 2.0 * d + f;  # GF 20101013 fixed typo d*s ->2*d
                    else:
                        s = -d / a
                        sqrdistance = d * s + f
            else:
                ## region 0
                invDet = 1.0 / det
                s = s * invDet
                t = t * invDet
                sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f
    else:
        if s < 0.0:
            ## region 2
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0:  # minimum on edge s+t=1
                numer = tmp1 - tmp0
                denom = a - 2.0 * b + c
                if numer >= denom:
                    s = 1.0
                    t = 0.0
                    sqrdistance = a + 2.0 * d + f;  # GF 20101014 fixed typo 2*b -> 2*d
                else:
                    s = numer / denom
                    t = 1 - s
                    sqrdistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f

            else:  # minimum on edge s=0
                s = 0.0
                if tmp1 <= 0.0:
                    t = 1
                    sqrdistance = c + 2.0 * e + f
                else:
                    if e >= 0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        t = -e / c
                        sqrdistance = e * t + f
                        # of region 2
        else:
            if t < 0.0:
                # region6
                tmp0 = b + e
                tmp1 = a + d
                if tmp1 > tmp0:
                    numer = tmp1 - tmp0
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        t = 1.0
                        s = 0
                        sqrdistance = c + 2.0 * e + f
                    else:
                        t = numer / denom
                        s = 1 - t
                        sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f

                else:
                    t = 0.0
                    if tmp1 <= 0.0:
                        s = 1
                        sqrdistance = a + 2.0 * d + f
                    else:
                        if d >= 0.0:
                            s = 0.0
                            sqrdistance = f
                        else:
                            s = -d / a
                            sqrdistance = d * s + f
            else:
                ## region 1
                numer = c + e - b - d
                if numer <= 0:
                    s = 0.0
                    t = 1.0
                    sqrdistance = c + 2.0 * e + f
                else:
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        s = 1.0
                        t = 0.0
                        sqrdistance = a + 2.0 * d + f
                    else:
                        s = numer / denom
                        t = 1 - s
                        sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f

    # account for numerical round-off error
    if sqrdistance < 0:
        sqrdistance = 0

    dist = sqrt(sqrdistance)

    #PP0 = B + s * E0 + t * E1
    return dist  #, PP0


@jit("f4(f4[:], f8[:,:])", nopython=True, cache=True)
def distance(P, TRI):
    # function [dist,PP0] = pointTriangleDistance(TRI,P)
    # calculate distance between a point and a triangle in 3D
    # SYNTAX
    #   dist = pointTriangleDistance(TRI,P)
    #   [dist,PP0] = pointTriangleDistance(TRI,P)
    #
    # DESCRIPTION
    #   Calculate the distance of a given point P from a triangle TRI.
    #   Point P is a row vector of the form 1x3. The triangle is a matrix
    #   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
    #   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
    #   to the triangle TRI.
    #   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
    #   closest point PP0 to P on the triangle TRI.
    #
    # Author: Gwolyn Fischer
    # Release: 1.0
    # Release date: 09/02/02
    # Release: 1.1 Fixed Bug because of normalization
    # Release: 1.2 Fixed Bug because of typo in region 5 20101013
    # Release: 1.3 Fixed Bug because of typo in region 2 20101014

    # Possible extention could be a version tailored not to return the distance
    # and additionally the closest point, but instead return only the closest
    # point. Could lead to a small speed gain.

    # Example:
    # %% The Problem
    # P0 = [0.5 -0.3 0.5]
    #
    # P1 = [0 -1 0]
    # P2 = [1  0 0]
    # P3 = [0  0 0]
    #
    # vertices = [P1; P2; P3]
    # faces = [1 2 3]
    #
    # %% The Engine
    # [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0)
    #
    # %% Visualization
    # [x,y,z] = sphere(20)
    # x = dist*x+P0(1)
    # y = dist*y+P0(2)
    # z = dist*z+P0(3)
    #
    # figure
    # hold all
    # patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8)
    # plot3(P0(1),P0(2),P0(3),'b*')
    # plot3(PP0(1),PP0(2),PP0(3),'*g')
    # surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
    # view(3)

    # The algorithm is based on
    # "David Eberly, 'Distance Between Point and Triangle in 3D',
    # Geometric Tools, LLC, (1999)"
    # http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    #
    #        ^t
    #  \     |
    #   \reg2|
    #    \   |
    #     \  |
    #      \ |
    #       \|
    #        *P2
    #        |\
    #        | \
    #  reg3  |  \ reg1
    #        |   \
    #        |reg0\
    #        |     \
    #        |      \ P1
    # -------*-------*------->s
    #        |P0      \
    #  reg4  | reg5    \ reg6
    # rewrite triangle in normal form
    B = TRI[0, :]
    E0 = TRI[1, :] - B
    # E0 = E0/sqrt(sum(E0.^2)); %normalize vector
    E1 = TRI[2, :] - B
    # E1 = E1/sqrt(sum(E1.^2)); %normalize vector
    D = B - P
    a = dot(E0, E0)
    b = dot(E0, E1)
    c = dot(E1, E1)
    d = dot(E0, D)
    e = dot(E1, D)
    f = dot(D, D)

    #print "{0} {1} {2} ".format(B,E1,E0)
    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e

    # Terible tree of conditionals to determine in which region of the diagram
    # shown above the projection of the point into the triangle-plane lies.
    if (s + t) <= det:
        if s < 0.0:
            if t < 0.0:
                # region4
                if d < 0:
                    t = 0.0
                    if -d >= a:
                        s = 1.0
                        sqrdistance = a + 2.0 * d + f
                    else:
                        s = -d / a
                        sqrdistance = d * s + f
                else:
                    s = 0.0
                    if e >= 0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        if -e >= c:
                            t = 1.0
                            sqrdistance = c + 2.0 * e + f
                        else:
                            t = -e / c
                            sqrdistance = e * t + f

                            # of region 4
            else:
                ## region 3
                s = 0
                if e >= 0:
                    t = 0
                    sqrdistance = f
                else:
                    if -e >= c:
                        t = 1
                        sqrdistance = c + 2.0 * e + f
                    else:
                        t = -e / c
                        sqrdistance = e * t + f
                        # of region 3
        else:
            if t < 0:
                ## region 5
                t = 0
                if d >= 0:
                    s = 0
                    sqrdistance = f
                else:
                    if -d >= a:
                        s = 1
                        sqrdistance = a + 2.0 * d + f;  # GF 20101013 fixed typo d*s ->2*d
                    else:
                        s = -d / a
                        sqrdistance = d * s + f
            else:
                ## region 0
                invDet = 1.0 / det
                s = s * invDet
                t = t * invDet
                sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f
    else:
        if s < 0.0:
            ## region 2
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0:  # minimum on edge s+t=1
                numer = tmp1 - tmp0
                denom = a - 2.0 * b + c
                if numer >= denom:
                    s = 1.0
                    t = 0.0
                    sqrdistance = a + 2.0 * d + f;  # GF 20101014 fixed typo 2*b -> 2*d
                else:
                    s = numer / denom
                    t = 1 - s
                    sqrdistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f

            else:  # minimum on edge s=0
                s = 0.0
                if tmp1 <= 0.0:
                    t = 1
                    sqrdistance = c + 2.0 * e + f
                else:
                    if e >= 0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        t = -e / c
                        sqrdistance = e * t + f
                        # of region 2
        else:
            if t < 0.0:
                # region6
                tmp0 = b + e
                tmp1 = a + d
                if tmp1 > tmp0:
                    numer = tmp1 - tmp0
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        t = 1.0
                        s = 0
                        sqrdistance = c + 2.0 * e + f
                    else:
                        t = numer / denom
                        s = 1 - t
                        sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f

                else:
                    t = 0.0
                    if tmp1 <= 0.0:
                        s = 1
                        sqrdistance = a + 2.0 * d + f
                    else:
                        if d >= 0.0:
                            s = 0.0
                            sqrdistance = f
                        else:
                            s = -d / a
                            sqrdistance = d * s + f
            else:
                ## region 1
                numer = c + e - b - d
                if numer <= 0:
                    s = 0.0
                    t = 1.0
                    sqrdistance = c + 2.0 * e + f
                else:
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        s = 1.0
                        t = 0.0
                        sqrdistance = a + 2.0 * d + f
                    else:
                        s = numer / denom
                        t = 1 - s
                        sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f

    # account for numerical round-off error
    if sqrdistance < 0:
        sqrdistance = 0

    dist = sqrt(sqrdistance)

    #PP0 = B + s * E0 + t * E1
    return dist  #, PP0


@jit("f8(f8[:],f8[:],f8[:],f8[:])", nopython=True, cache=True)
def vol(a, b, c, d):
    M = np.zeros((3, 3))
    M[0, :] = np.subtract(a, d)
    M[1, :] = np.subtract(b, d)
    M[2, :] = np.subtract(c, d)
    return np.abs(np.linalg.det(M)) / 6


@jit("f8(f8[:,:])", nopython=True, cache=True)
def get_radius(points):
    a = np.linalg.norm(points[0] - points[1])
    a1 = np.linalg.norm(points[2] - points[3])
    b = np.linalg.norm(points[0] - points[2])
    b1 = np.linalg.norm(points[1] - points[3])
    c = np.linalg.norm(points[0] - points[3])
    c1 = np.linalg.norm(points[1] - points[2])
    p = (a * a1 + b * b1 + c * c1) / 2
    V = vol(points[0], points[1], points[2], points[3])
    if V > 0:
        return 1 / (6 * V) * np.sqrt(p * (p - a * a1) * (p - b * b1) * (p - c * c1))
    else:
        return np.inf


@jit("i4[:,:](f8[:,:], i4[:,:], f8)", nopython=True, cache=True)
def _calc_alpha_simplex(C, S, a):
    M = S.shape[0]
    N = S.shape[1]
    Result = np.zeros((M, N))
    j = 0
    for i in range(M):
        s = S[i, :]
        ps = C[s]
        r = get_radius(ps)
        if r < a:
            Result[j, :] = s
            j += 1
    return Result[:j, :].astype(np.int32)


combos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])


@jit("i4[:,:](i4[:,:], i8[:,:])", nopython=True, cache=True)
def make_simplex_set(S, combos):
    M = S.shape[0] * 4
    N = S.shape[1] - 1
    R = np.zeros((M, N), dtype=np.int32)
    for k, s in enumerate(range(0, M, 4)):
        for i in range(4):
            for j in range(3):
                R[s + i, j] = S[k, combos[i, j]]
    return R


@jit("f8[:](f4[:], i4[:,:], f8[:,:])", nopython=True, cache=True, parallel=True)
def get_distances(pnt, hull, coords):
    S = hull.shape[0]
    distances = np.zeros(S)
    for i in prange(S):
        distances[i] = distance(pnt, coords[hull[i]])
    return distances


@jit("f8[:](f8[:], i4[:,:], f8[:,:])", nopython=True, cache=True, parallel=True)
def get_distances_(pnt, hull, coords):
    S = hull.shape[0]
    distances = np.zeros(S)
    for i in prange(S):
        distances[i] = distance_(pnt, coords[hull[i]])
    return distances


class AlphaHull:

    def __init__(self, alpha=9):
        self.alpha = alpha
        self.hull = None
        self._hull = None
        self.tri = None
        self._tri = None
        self.coords = None
        self.simplices = None
        self.resindices = None

    def set_coords(self, pdb):
        type1 = isinstance(pdb, pr.atomic.selection.Selection)
        type2 = isinstance(pdb, pr.atomic.atomgroup.AtomGroup)
        if type1 or type2:
            self._set_coords(pdb)
        elif isinstance(pdb, np.ndarray):
            self.coords = pdb
        else:
            raise ValueError('*pdb* must be prody instance or numpy array')

    def _set_coords(self, pdb):
        """pdb is a prody object. pdb should have CB atoms where appropriate."""
        self.coords = pdb.select('name CB or (resname GLY and name CA)').getCoords()

    def set_tri(self):
        self.tri = Delaunay(self.coords)
        self._tri = copy.deepcopy(self.tri)

    def set_resindices(self, pdb):
        """pdb is a prody object. pdb should have CB atoms where appropriate."""
        self.resindices = pdb.select('name CB or (resname GLY and name CA)').getResindices()

    def calc_alpha_simplices(self):
        if self.tri is None:
            self.set_tri()
        self.tri.simplices.sort()
        self.tri.simplices = self.tri.simplices[self.tri.simplices[:, 0].argsort()]
        self.simplices = _calc_alpha_simplex(self.coords, self.tri.simplices, self.alpha)
        self._tri.simplices = self.simplices
        self._tri.neighbors = self.simplices

    def calc_hull(self):
        if self.simplices is None:
            self.calc_alpha_simplices()
        simpl_set = make_simplex_set(self.simplices, combos)
        un, ind, co = np.unique(simpl_set, axis=0,
                                return_counts=True, return_index=True)
        self.hull = np.array([simpl_set[i] for i in ind[co == 1]], dtype=np.int32)

    def pnts_in_hull(self, pnts):
        return self._tri.find_simplex(pnts) >= 0

    def get_pnt_distance(self, pnt):
        distances = get_distances(pnt, self.hull, self.coords)
        inout = self.pnts_in_hull(pnt)
        if inout:
            return distances.min()
        else:
            return (-1 * distances).max()

    def get_pnt_distance_(self, pnt):
        distances = get_distances_(pnt, self.hull, self.coords)
        inout = self.pnts_in_hull(pnt)
        if inout:
            return distances.min()
        else:
            return (-1 * distances).max()

    def get_pnts_distance(self, pnts):
        if pnts.dtype == 'float32':
            return [self.get_pnt_distance(pnt) for pnt in pnts]
        elif pnts.dtype == 'float64':
            return [self.get_pnt_distance_(pnt) for pnt in pnts]
        else:
            raise TypeError('Coordinates must have dtype float32 or float64.')

    def get_volume(self):
        if self.hull is None:
            self.calc_hull()
        vertices = np.unique(self.hull.flatten())
        # calculate hull centroid
        d = np.mean(self.coords[vertices], axis=0)
        volume = 0.
        for face in self.hull:
            a, b, c = (self.coords[face[0]], 
                       self.coords[face[1]], 
                       self.coords[face[2]])
            mat = np.zeros((3, 3))
            mat[0], mat[1], mat[2] = b - a, c - a, d - a
            volume += np.abs(np.linalg.det(mat)) / 6.
        return volume

    def get_surface_area(self):
        if self.hull is None:
            self.calc_hull()
        area = 0.
        for face in self.hull:
            a, b, c = (self.coords[face[0]], 
                       self.coords[face[1]], 
                       self.coords[face[2]])
            area += np.linalg.norm(np.cross(b - a, c - a)) / 2.
        return area


def partition_res_by_burial(pdb_ala, alpha=9, ahull_ca=None, ahull_cb=None, 
                            assign_intermediate_by_distance=False,
                            distance_threshold=-1.0,):
    """Returns residue indices of exposed, intermediate, and buried residues
    based on CA hull and CB hull."""
    if ahull_ca is None:
        ahull_ca = AlphaHull(alpha=alpha)
        ahull_ca.coords = pdb_ala.select('name CA').getCoords()
        ahull_ca.calc_hull()
    if ahull_cb is None:
        ahull_cb = AlphaHull(alpha=alpha)
        ahull_cb.set_coords(pdb_ala)
        ahull_cb.calc_hull()
    ahull_cb.set_resindices(pdb_ala)
    cb_in_ca_hull = ahull_ca.pnts_in_hull(ahull_cb.coords)
    resindices_cb_in_ca_hull = set(ahull_cb.resindices[cb_in_ca_hull])
    resindices_cb_hull = set(ahull_cb.resindices[np.unique(ahull_cb.hull)])
    resindices_not_cb_hull = set(ahull_cb.resindices) - resindices_cb_hull
    resindices_exposed = resindices_cb_hull - resindices_cb_in_ca_hull
    resindices_intermediate = resindices_cb_in_ca_hull - resindices_not_cb_hull # CB is part of CB hull, but CB in CA hull.
    resindices_buried = resindices_cb_in_ca_hull - resindices_intermediate
    res_ = resindices_not_cb_hull - resindices_buried
    resindices_intermediate |= res_

    if assign_intermediate_by_distance:
        for resindex in resindices_exposed.copy():
            if ahull_ca.get_pnt_distance_(ahull_cb.coords[resindex]) > distance_threshold:
                resindices_intermediate.add(resindex)
                resindices_exposed.remove(resindex)

    return resindices_exposed, resindices_intermediate, resindices_buried



def main():
    par = argparse.ArgumentParser()
    par.add_argument('-p', '--path_to_pdb', required=True, help='path to pdb file or pdb accession code (for fetching)')
    par.add_argument('-s', '--sel', default=None, help='prody selection string within pdb file')
    par.add_argument('-o', '--outdir', default='./', help='output directory')
    par.add_argument('-f', '--filename', default='packdensity.txt', help='filename for output')
    par.add_argument('-n', '--numthreads', default=4, help='number of threads for calculation')
    par.add_argument('-q', '--numpoints', default=1000000, help='number of sample points for MC integration')
    par.add_argument('-c', '--chunksize', default=100000, help='chunk size for sample points distance calc')
    par.add_argument('-w', '--overwrite', action='store_true', help='overwrite previous output file. Default is to append.')
    args = par.parse_args()

    path_to_pdb = args.path_to_pdb
    sel = args.sel
    outdir = args.outdir
    filename = args.filename
    numthreads = int(args.numthreads)
    numpoints = int(args.numpoints)
    overwrite = args.overwrite
    chunk_size = args.chunksize

    calc_packing_density(path_to_pdb=path_to_pdb, 
                         outdir=outdir, 
                         filename=filename,
                         overwrite=overwrite,
                         selection=sel, 
                         n_points=numpoints, 
                         chunk_size=chunk_size, 
                         n_workers=numthreads)


if __name__ == '__main__':
    main()