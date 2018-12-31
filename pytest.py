#!/usr/bin/python3
import os
import sys
import numpy as np
import spatial
import neighbour_list

print('OMP_NUM_THREADS:', os.environ["OMP_NUM_THREADS"])

# ALERT: must keep it consistent with include/constants.h
MAX_N_NEIGHBOURS = 80

np.random.seed(100)

r1 = np.zeros(3)
r2 = np.ones(3)
# box = np.array([
#      [19.0, 0.00, 0.00],
#      [0.0, 19.0, 0.0],
#      [0.0, 0.0, 19.0]
#      ])
box = np.array([
     [30.0, 0.00, 0.00],
     [0.0, 40.0, 0.0],
     [0.0, 0.0, 35.0]
     ])
# box = np.array([
#      [30.0, 0.00, 0.00],
#      [-17.320508075688775, 30.0, 0.0],
#      [0.0, 0.0, 35.0]
#      ])
box_inv = np.linalg.inv(box)
#box_inv = np.identity(3)/30.0

print(spatial.distance_pbc(r1, r2, box, box_inv, 0))

# construct inputs
n_atoms = 4000
r_cut = 4.0
spositions = np.random.random((n_atoms, 3))
positions = spositions.dot(box)
atom_labels = np.arange(n_atoms)

# build cell list
celllist_1d = neighbour_list.construct_cell_list(positions.T, atom_labels, n_atoms, box.T, box_inv.T, r_cut)
cell_index_map = neighbour_list.build_cell_index_map(celllist_1d, n_atoms)
cell_index_map = cell_index_map.T
# find all neighbours
distances2 = np.zeros((MAX_N_NEIGHBOURS, n_atoms), order='F')
for i in range(1000):
    nbs = neighbour_list.find_neighbours_for_all(positions.T, celllist_1d, n_atoms, r_cut, box.T, box_inv.T, distances2)
nbs = nbs.T
distances2 = distances2.T
dimensions = np.flip(celllist_1d[0:4])
celllist = celllist_1d[4:].reshape(dimensions)


################################
# check data in the cell list
################################
# print('*****, printing some randomly selected elements in neighbour list')
# print('---cell_list components---')
# print(celllist.shape)
# i = 5
# j = 6
# k = 6
# na, nb, nc, nmax = dimensions
# print(celllist[i, j, k, :])
# ib = nmax * (k + j*nc + i*nc*nb)
# ie = nmax * (k + j*nc + i*nc*nb + 1)
# print(celllist_1d[4+ib:4+ie])
# # print('---cell_index_map_test---')
# # print(positions[392])
# # print(box[0]/dimensions[0], box[1]/dimensions[1], box[2]/dimensions[2])
# # print(cell_index_map[392])
# print('---neighbours---')
# i0 = 1
# i1 = 1809
# print('indices')
# print(i0, i1)
# print(nbs[i0])
# print(np.sqrt(distances2)[i0,:])
# print(nbs[i1])
# # print(neighbour_list.find_neighbours_for_atom(positions.T, celllist_1d, i0, cell_index_map.T, r_cut, box.T, box_inv.T))
# print('positions')
# r0 = positions[i0]
# r1 = positions[i1]
# print(r0);
# print(r1);
# print('distances')
# print(spatial.distance_pbc(r0, r1, box.T, box_inv.T, 0))
# print('cell_index_map')
# print(cell_index_map[i0])
# print(cell_index_map[i1])

################################
# systematic test neighbour list
################################
flag = 0
print('*****, Start testing nbs and distances2...')
for i_atom0, neighbours in enumerate(nbs):
    r0 = positions[i_atom0]
    print('Testing particle', i_atom0)
    for i_atom1 in range(n_atoms):
        if i_atom1 == i_atom0:
            continue 
        r1 = positions[i_atom1]
        dr = spatial.distance_pbc(r0, r1, box.T, box_inv.T, 0)
        if dr >= r_cut and i_atom1 in neighbours:
            print('ERROR: Wrong interacting pair:', i_atom0, i_atom1, dr)
            flag = 1
        elif dr < r_cut:
            index = np.where(neighbours==i_atom1)[0]
            if len(index) == 0:
                print('ERROR: Missing interacting pair:', i_atom0, i_atom1, dr)
                flag = 1
            elif len(index) > 1:
                print('ERROR: Duplication in neighbour list:', i_atom0, i_atom1, dr)
                flag = 1
            else:
                if abs(np.sqrt(distances2[i_atom0, index[0]]) - dr) > 1e-6:
                    print('ERROR: Distances mismatch:', i_atom0, i_atom1, dr, np.sqrt(distances2[i_atom0, index[0]]))
                    flag = 1
if flag == 0:
   print('Neighbour list checks out.')
        

print('---number of neighbours---')
print('max', 'min', 'aver')
n_nbs = ((nbs>0)>0).sum(axis=1)
print(n_nbs.max(), n_nbs.min(), np.average(n_nbs))
