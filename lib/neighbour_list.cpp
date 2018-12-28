/*
 * =====================================================================================
 *
 *       Filename:  neighbour_list.cpp
 *
 *    Description:  library that takes care of neighbour list search
 *
 *        Version:  1.0
 *        Created:  12/26/2018 01:59:46 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (Kuang Yu), 
 *   Organization:  TBSI
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <mylibs.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

using namespace std;
using namespace Eigen;

// this function returns a cell neighbour list
// positions: n_atoms*3 array, storing the positions of the atoms you want to neighbour list
//              Here, positions MUST BE dense.
// atom_labels: n_atoms*1 array, just integers to label the atoms, usually just 0, 1, ... n_atoms-1
//              but can be different if we are constructing a neighbour list for only part of the system
//              for example, if we only want to neighbour list the real atoms, but not the v-sites, then
//              the v-site indices should be missing from atom_labels
//              that is, atom_labels must NOT be dense.
// n_atoms: integer, must be dense
// box, box_inv: box dimensions and its inversion
// r_cut: realspace cutoff distance, in Angstrom
Tensor<int, 4> construct_cell_list(const MatrixXd & positions, const ArrayXi & atom_labels,
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut)
{
  Tensor<int, 4> celllist;
  double vol = volume(box);
  double ab, bc, ca;
  int na, nb, nc, n_max;
  Vector3d a, b, c;
  a = box.col(0);
  b = box.col(1);
  c = box.col(2);
  ab = (a.cross(b)).norm();
  bc = (b.cross(c)).norm();
  ca = (c.cross(a)).norm();
  // perpendicular heights along each dimension
  double ha, hb, hc;
  ha = vol / bc;
  hb = vol / ca;
  hc = vol / ab;
  // number of cells in each dimension
  na = floor(ha/r_cut);
  nb = floor(hb/r_cut);
  nc = floor(hc/r_cut);
  double density = (double)n_atoms/na/nb/nc;
  // Note here +N is the magic number, we assume the density fluctuation is small enough
  // so +N is large enough buffer.
  n_max = floor(density) + CELL_LIST_BUFFER;
  // create celllist, initialize with -1
  // -1 signals the end of the list
  celllist.resize(n_max, nc, nb, na);
  celllist.setConstant(-1);

  MatrixXd cell=box, cell_inv;
  cell.col(0) /= na;
  cell.col(1) /= nb;
  cell.col(2) /= nc;
  cell_inv = cell.inverse();

  // cell indices of each particle
  MatrixXd spositions = cell_inv * positions;
  MatrixXi indices = spositions.array().floor().cast<int>();
  for (int i_atom=0; i_atom<n_atoms; i_atom++) {
    // pbc shifts
    int label = atom_labels(i_atom);
    int ia = indices(0, i_atom) % na;
    int ib = indices(1, i_atom) % nb;
    int ic = indices(2, i_atom) % nc;
    if (indices(0, i_atom) < 0) ia += na;
    if (indices(0, i_atom) < 0) ib += nb;
    if (indices(0, i_atom) < 0) ic += nc;
    // assign atom to a particular cell
    int i;
    for (i=0; i<n_max; i++) {
      if (celllist(i, ic, ib, ia) < 0) {
        celllist(i, ic, ib, ia) = label;
        break;
      }
    }
    // overflow
    if (i == n_max) {
      cout << "ERROR: cell list overflow, increase CELL_LIST_BUFFER" << endl;
      throw exception();
    }
  }

  return celllist;
}


// construct a n_atoms*3 array, stores the cell indices in which each particle resides
// Watchout: for a system with virtual sites, n_atoms should be the number of all sites
// including both atoms (which you want to neighbour list), and the vsites (which you do
// not want to neighbour list.
// Hence, the returned cell_index_map will contain -1 for virtual sites, indicating that 
// they are not assigned to a particular cell.
// n_atoms: integer, must be not dense
// return values not dense
MatrixXi build_cell_index_map(const Tensor<int, 4> & celllist, int n_atoms) 
{
  MatrixXi cell_index_map;
  cell_index_map.resize(3, n_atoms);
  cell_index_map.setConstant(-1);
  int na = celllist.dimension(3);
  int nb = celllist.dimension(2);
  int nc = celllist.dimension(1);
  int nmax = celllist.dimension(0);
  for (int ia=0; ia<na; ia++){
    for (int ib=0; ib<nb; ib++) {
      for (int ic=0; ic<nc; ic++) {
        for (int n=0; n<nmax; n++) {
          if (celllist(n, ic, ib, ia) < 0) break;
          int iatom = celllist(n, ic, ib, ia);
          cell_index_map(0, iatom) = ia;
          cell_index_map(1, iatom) = ib;
          cell_index_map(2, iatom) = ic;
        }
      }
    }
  }
  return cell_index_map;
}

// Here, positions must NOT be dense
// cell_index_map must NOT be dense
// i_atom must NOT be dense
ArrayXi find_neighbours(const MatrixXd & positions, const Tensor<int, 4> & celllist,
    int i_atom, int i0, int j0, int k0, double r_cut, const MatrixXd & box, const MatrixXd & box_inv)
{
  // original cell indices
  ArrayXi neighbours;
  neighbours.resize(MAX_N_NEIGHBOURS);
  neighbours.setConstant(-1);
  int n_neighbour = 0;
  int n1 = celllist.dimension(3);
  int n2 = celllist.dimension(2);
  int n3 = celllist.dimension(1);
  int nmax = celllist.dimension(0);
  double r_cut2 = r_cut*r_cut;
  VectorXd r0 = positions.col(i_atom);
  // loop over the 3*3*3=27 neighbouring cells
  for (int di=-1; di<=1; di++) {
    int i = i0 + di;
    // pbc shifts
    i = i % n1;
    if (i < 0) i += n1;
    for (int dj=-1; dj<=1; dj++) {
      int j = j0 + dj;
      j = j % n2;
      if (j < 0) j += n2;
      for (int dk=-1; dk<=1; dk++) {
        int k = k0 + dk;
        k = k % n3;
        if (k < 0) k += n3;
        for (int n=0; n<nmax; n++) {
          // end of cell list
          if (celllist(n, k, j, i) < 0) break;
          int i_atom1 = celllist(n, k, j, i);
          if (i_atom1 == i_atom) continue;
          // if r_cut < 0, meaning we do not do distance judge
          if (r_cut < 0) {
            if (n_neighbour == MAX_N_NEIGHBOURS) {
              cout << "ERROR: overflow in neighbour search, increase MAX_N_NEIGHBOURS" << endl;
              throw exception();
            }
            neighbours(n_neighbour) = i_atom1;
            n_neighbour += 1;
          }
          else {
            VectorXd r1 = positions.col(i_atom1);
            VectorXd dr = dr_vec_pbc(r0, r1, box, box_inv, 0);
            if (abs(dr(0))>r_cut or abs(dr(1))>r_cut or abs(dr(2))>r_cut) continue;
            if (dr.squaredNorm() < r_cut2) {
              if (n_neighbour == MAX_N_NEIGHBOURS) {
                cout << "ERROR: overflow in neighbour search, increase MAX_N_NEIGHBOURS" << endl;
                throw exception();
              }
              neighbours(n_neighbour) = i_atom1;
              n_neighbour += 1;
            }
          }
        }
      }
    }
  }
  return neighbours;
}


ArrayXi find_neighbours_for_atom(const MatrixXd & positions, const Tensor<int, 4> & celllist,
    int i_atom0, const MatrixXi & cell_index_map, double r_cut, const MatrixXd & box, const MatrixXd & box_inv)
{
  int i0 = cell_index_map(0, i_atom0);
  int j0 = cell_index_map(1, i_atom0);
  int k0 = cell_index_map(2, i_atom0);
  return find_neighbours(positions, celllist, i_atom0, i0, j0, k0, r_cut, box, box_inv);
}


// Positions here are not dense
// n_atoms here is not dense
MatrixXi find_neighbours_for_all(const MatrixXd & positions, const Tensor<int, 4> & celllist,
    int n_atoms, double r_cut, const MatrixXd & box, const MatrixXd & box_inv)
{
  MatrixXi nb_list;
  nb_list.resize(MAX_N_NEIGHBOURS, n_atoms);
  nb_list.setConstant(-1);
  ArrayXi n_nbs;
  n_nbs.resize(n_atoms);
  n_nbs.setConstant(0);
  int n1 = celllist.dimension(3);
  int n2 = celllist.dimension(2);
  int n3 = celllist.dimension(1);
  int nmax = celllist.dimension(0);
  double r_cut2 = r_cut * r_cut;
  for (int i0=0; i0<n1; i0++) {
    for (int j0=0; j0<n2; j0++) {
      for (int k0=0; k0<n3; k0++) {
        for (int n0=0; n0<nmax; n0++) {
          // finish looping over this cell
          if (celllist(n0, k0, j0, i0) < 0) break;
          int i_atom0 = celllist(n0, k0, j0, i0);
          VectorXd r0 = positions.col(i_atom0);
          // loop over the 3*3*3=27 neighbouring cells
          for (int di=-1; di<=1; di++) {
            int i = i0 + di;
            i = i % n1;
            if (i < 0) i += n1;
            for (int dj=-1; dj<=1; dj++) {
              int j = j0 + dj;
              j = j % n2;
              if (j < 0) j += n2;
              for (int dk=-1; dk<=1; dk++) {
                int k = k0 + dk;
                k = k % n3;
                if (k < 0) k += n3;
                for (int n=0; n<nmax; n++) {
                  if (celllist(n, k, j, i) < 0) break;
                  int i_atom1 = celllist(n, k, j, i);
                  if (i_atom1 <= i_atom0) continue;
                  // do not judge distance, just put atom0, atom1 be neighbours
                  if (r_cut < 0.0) { 
                    if (n_nbs(i_atom0) == MAX_N_NEIGHBOURS or n_nbs(i_atom1) == MAX_N_NEIGHBOURS) {
                      cout << "ERROR: overflow in neighbour search, increase MAX_N_NEIGHBOURS" << endl;
                      throw exception();
                    }
                    // add neighbour to atom0
                    nb_list(n_nbs(i_atom0), i_atom0) = i_atom1;
                    n_nbs(i_atom0) += 1;
                    // add neighbour to atom1
                    nb_list(n_nbs(i_atom1), i_atom1) = i_atom0;
                    n_nbs(i_atom1) += 1;
                  }
                  else {
                    VectorXd r1 = positions.col(i_atom1);
                    VectorXd dr = dr_vec_pbc(r0, r1, box, box_inv, 0);
                    // pre-screening
                    if (abs(dr(0))>r_cut or abs(dr(1))>r_cut or abs(dr(2))>r_cut) continue;
                    if (dr.squaredNorm() < r_cut2) {
                      if (n_nbs(i_atom0) == MAX_N_NEIGHBOURS or n_nbs(i_atom1) == MAX_N_NEIGHBOURS) {
                        cout << "ERROR: overflow in neighbour search, increase MAX_N_NEIGHBOURS" << endl;
                        throw exception();
                      }
                      // add neighbour to atom0
                      nb_list(n_nbs(i_atom0), i_atom0) = i_atom1;
                      n_nbs(i_atom0) += 1;
                      // add neighbour to atom1
                      nb_list(n_nbs(i_atom1), i_atom1) = i_atom0;
                      n_nbs(i_atom1) += 1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return nb_list;
}


/* ************************************************ */
/* Python wrappers for the neighbour list functions */
/* ************************************************ */

// pybind11 does not recognize the multidimensional Tensor datatype, hence the beaucracy

// convert a 4d tensor to 1d array (assume column based layout)
// T(n1, n2, n3, n4) -> A = {n1, n2, n3, n4, T{1111}, T{2111}, ..., T{n1,n2,n3,n4}}
ArrayXi Tensor_4d_to_1d(const Tensor<int, 4> & t_4d)
{
  ArrayXi results;
  Map<const ArrayXi> data(t_4d.data(), t_4d.size());
  results.resize(t_4d.size()+4);
  // copy the dimension information
  results(0) = t_4d.dimension(0);
  results(1) = t_4d.dimension(1);
  results(2) = t_4d.dimension(2);
  results(3) = t_4d.dimension(3);
  // copy the real data
  results.segment(4, t_4d.size()) = data;
  return results;
}

// convert a 1d array to 4d tensor
Tensor<int, 4> Tensor_1d_to_4d(const ArrayXi & A_1d)
{
  int n1 = A_1d(0);
  int n2 = A_1d(1);
  int n3 = A_1d(2);
  int n4 = A_1d(3);
  TensorMap<Tensor<const int, 4>> T_4d(A_1d.segment(4, A_1d.size()-4).data(), n1, n2, n3, n4);
  return T_4d;
}

ArrayXi construct_cell_list_py(const MatrixXd & positions, const ArrayXi & atom_labels,
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut)
{
  Tensor<int, 4> celllist;
  celllist = construct_cell_list(positions, atom_labels, n_atoms, box, box_inv, r_cut);
  return Tensor_4d_to_1d(celllist);
}


MatrixXi build_cell_index_map_py(const ArrayXi & celllist_py, int n_atoms) 
{
  Tensor<int, 4> celllist;
  celllist = Tensor_1d_to_4d(celllist_py);
  return build_cell_index_map(celllist, n_atoms);
}


ArrayXi find_neighbours_for_atom_py(const MatrixXd & positions, const ArrayXi & celllist_py,
    int i_atom0, const MatrixXi & cell_index_map, double r_cut, const MatrixXd & box, const MatrixXd & box_inv)
{
  Tensor<int, 4> celllist;
  celllist = Tensor_1d_to_4d(celllist_py);
  return find_neighbours_for_atom(positions, celllist, i_atom0, cell_index_map, r_cut, box, box_inv);
}


MatrixXi find_neighbours_for_all_py(const MatrixXd & positions, const ArrayXi & celllist_py,
    int n_atoms, double r_cut, const MatrixXd & box, const MatrixXd & box_inv)
{
  Tensor<int, 4> celllist;
  celllist = Tensor_1d_to_4d(celllist_py);
  return find_neighbours_for_all(positions, celllist, n_atoms, r_cut, box, box_inv);
}  


PYBIND11_MODULE(neighbour_list, m)
{
  m.doc() = "Neighbour cell list module";
  m.def("construct_cell_list", &construct_cell_list_py);
  m.def("build_cell_index_map", &build_cell_index_map_py);
  m.def("find_neighbours_for_atom", &find_neighbours_for_atom_py);
  m.def("find_neighbours_for_all", &find_neighbours_for_all_py);
}
