#!/usr/bin/python3
import os
import sys
import numpy as np
import spatial
import neighbour_list

r1 = np.zeros(3, dtype=np.float64)
r2 = np.ones(3, dtype=np.float64)
box = np.identity(3, dtype=np.float64)*30
box_inv = np.linalg.inv(box).astype(np.float64)

print(spatial.distance_pbc(r1, r2, box, box_inv, 0))

n_atoms = 1000
r_cut = 4.0
spositions = np.random.random((n_atoms, 3))
positions = spositions.dot(box)
atom_labels = np.arange(n_atoms)
celllist = neighbour_list.construct_cell_list_pywrap(positions, atom_labels, n_atoms, box, box_inv, r_cut)
dimensions = celllist[0:4]
celllist = celllist[4:].reshape(dimensions)
print(box)
print(celllist.shape)
print(celllist[6, 6, 6, :])
print(celllist[5, 6, 3, :])
print(celllist[2, 2, 1, :])
print(celllist[1, 1, 4, :])
print(positions[celllist[6, 6, 6, 0], :])
