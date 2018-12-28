#!/usr/bin/python3
import os
import sys
import numpy as np
import spatial
import neighbour_list

np.random.seed(100)

r1 = np.zeros(3, dtype=np.float64)
r2 = np.ones(3, dtype=np.float64)
box = np.identity(3, dtype=np.float64)*30
box_inv = np.linalg.inv(box).astype(np.float64)

print(spatial.distance_pbc(r1, r2, box, box_inv, 0))

n_atoms = 4000
r_cut = 4.0
spositions = np.random.random((n_atoms, 3))
positions = spositions.dot(box)
atom_labels = np.arange(n_atoms)
# build cell list
celllist_1d = neighbour_list.construct_cell_list(positions.T, atom_labels, n_atoms, box.T, box_inv.T, r_cut)
cell_index_map = neighbour_list.build_cell_index_map(celllist_1d, n_atoms)
cell_index_map = cell_index_map.T
nbs = neighbour_list.find_neighbours_for_all(positions.T, celllist_1d, n_atoms, r_cut, box.T, box_inv.T)
nbs = nbs.T

dimensions = np.flip(celllist_1d[0:4])
celllist = celllist_1d[4:].reshape(dimensions)
print('---cell_list components---')
print(celllist.shape)
i = 6
j = 6
k = 6
na, nb, nc, nmax = dimensions
print(celllist[i, j, k, :])
ib = nmax * (k + j*nc + i*nc*nb)
ie = nmax * (k + j*nc + i*nc*nb + 1)
print(celllist_1d[4+ib:4+ie])
print('---cell_index_map_test---')
print(positions[392])
print(box[0]/dimensions[0], box[1]/dimensions[1], box[2]/dimensions[2])
print(cell_index_map[392])
print('---neighbours---')
i0 = 392
i1 = 618
print('indices')
print(i0, i1)
print(nbs[i0])
print(nbs[i1])
print('positions')
r0 = positions[i0]
r1 = positions[i1]
print('distances')
print(spatial.distance_pbc(r0, r1, box, box_inv, 0))
print('cell_index_map')
print(cell_index_map[i0])
print(cell_index_map[i1])

################################
# systematic test neighbour list
################################
# for i_atom0, neighbours in enumerate(nbs):
#      r0 = positions[i_atom0]
#      print('Testing particle', i_atom0)
#      for i_atom1 in range(n_atoms):
#          if i_atom1 == i_atom0:
#              continue 
#          r1 = positions[i_atom1]
#          dr = spatial.distance_pbc(r0, r1, box, box_inv, 0)
#          if dr >= r_cut and i_atom1 in neighbours:
#              print('Wrong interacting pair:', i_atom0, i_atom1, dr)
#          elif dr < r_cut and i_atom1 not in neighbours:
#              print('Missing interacting pair:', i_atom0, i_atom1, dr)

print('---number of neighbours---')
print('max', 'min', 'aver')
n_nbs = (((nbs>0)*nbs)>0).sum(axis=1)
print(n_nbs.max(), n_nbs.min(), np.average(n_nbs))
