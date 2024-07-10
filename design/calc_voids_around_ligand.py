import prody as pr
from sklearn.neighbors import NearestNeighbors
import numpy as np
from scipy.spatial import Delaunay
from numba import jit
import copy
import argparse
import os
from sklearn.cluster import DBSCAN

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

def sample_points_on_sphere(center, radius, num_points):
    points = []
    offset = 2 / num_points
    increment = np.pi * (3 - np.sqrt(5))

    for i in range(num_points):
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - y*y)

        phi = ((i + 1) % num_points) * increment

        x = np.cos(phi) * r
        z = np.sin(phi) * r

        # Convert to Cartesian coordinates and adjust for center
        point = (x * radius + center[0], y * radius + center[1], z * radius + center[2])
        points.append(point)

    # # add gaussian noise to points
    # points = np.array(points)
    # noise = np.random.normal(0, 0.1, points.shape)
    # points += noise

    return points 


def _print_grid_pts(grid_pts, filename):
    ag = pr.AtomGroup('ligand')
    ag.setCoords(grid_pts)
    ag.setNames(['C'] * grid_pts.shape[0])
    ag.setResnames(['ALA'] * grid_pts.shape[0])
    ag.setResnums(np.arange(grid_pts.shape[0]))
    ag.setChids(['A'] * grid_pts.shape[0])
    pr.writePDB(filename, ag)


def get_grid_points(coords, radius, num_points=300, print_grid_pts=False, grid_pts_filename='grid_pts.pdb'):
    grid_pts = []
    for coord in coords:
        r = radius
        grid_pts.append(coord)
        while r > 0:
            grid_pts.extend(sample_points_on_sphere(coord, radius=r, num_points=num_points))
            r -= 1
    grid_pts = np.array(grid_pts)
    if print_grid_pts:
        _print_grid_pts(grid_pts, grid_pts_filename)
    return grid_pts


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
    

def calc_void_vol_around_lig(path_to_pdb,
                            lig_resname,
                            lig_chain,
                            outdir='./', 
                            selection=None,
                            print_cluster_volumes=False,
                            print_cluster_pts=False):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    if outdir[-1] != '/':
        outdir += '/'
    
    pdb = pr.parsePDB(path_to_pdb)
    if selection is not None:
        pdb = pdb.select(selection)

    if (pdb.getElements() == '').all():
        set_elements_sel(pdb)

    lig = pdb.select(f'resname {lig_resname} and chain {lig_chain}')
    lig_noH = lig.select('not element H')
    lig_noH_coords = lig_noH.getCoords()

    ahull_prot = AlphaHull()
    ahull_prot.set_coords(pdb.select('protein'))
    ahull_prot.calc_hull()

    grid_points = get_grid_points(lig_noH_coords, radius=8.0)
    grid_points = grid_points[ahull_prot.pnts_in_hull(grid_points)]

    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(pdb.select('not element H').getCoords())
    inds = nbrs.radius_neighbors(grid_points, radius=2.8, return_distance=False)
    grid_points = grid_points[[len(i) == 0 for i in inds]]

    if print_cluster_pts:
        _print_grid_pts(grid_points, f'{outdir}gridpts_void_unclustered.pdb')

    clustering = DBSCAN(eps=1.0, min_samples=6).fit(grid_points)
    cluster_labels = set(clustering.labels_) - {-1}
    vol = 0
    idx = 0
    for label in cluster_labels:
        cluster = grid_points[clustering.labels_ == label]
        if cluster.shape[0] < 12: # twice min_samples
            continue
        ahull = AlphaHull(alpha=3.0)
        ahull.set_coords(cluster)
        ahull.calc_hull()
        vol_ = ahull.get_volume()
        vol += vol_
        if print_cluster_volumes:
            print(idx, vol_)
        if print_cluster_pts:
            _print_grid_pts(cluster, f'{outdir}gridpts_void_{idx}.pdb')
        idx += 1

    #if overwrite:
        #with open(f'{outdir}{filename}', 'w') as f:
            #f.write(f'{path_to_pdb} {vol}\n')
    #else:
        #with open(f'{outdir}{filename}', 'a') as f:
            #f.write(f'{path_to_pdb} {vol}\n')

    return vol
