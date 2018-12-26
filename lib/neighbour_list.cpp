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
// positions: n_atoms*3 array, storing the position information of the atoms
// atom_labels: n_atoms*1 array, just integers to label the atoms, usually just 0, 1, ... n_atoms-1
//              but can be different if we are constructing a neighbour list for only part of atoms
// n_atoms: integer
// box, box_inv: box dimensions and its inversion
// r_cut: realspace cutoff distance, in Angstrom
Tensor<int, 4, RowMajor> construct_cell_list(const MatrixXd & positions, const ArrayXi & atom_labels,
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut)
{
  Tensor<int, 4, RowMajor> celllist;
  double vol = volume(box);
  double ab, bc, ca;
  int na, nb, nc, n_max;
  Vector3d a, b, c;
  a = box.row(0);
  b = box.row(1);
  c = box.row(2);
  ab = (a.cross(b)).norm();
  bc = (b.cross(c)).norm();
  ca = (c.cross(a)).norm();
  // perpendicular heights along each dimension
  double ha, hb, hc;
  ha = vol / bc;
  hb = vol / ca;
  hc = vol / ab;
  // number of cells in each dimension
  na = floor(ha/r_cut) + 1;
  nb = floor(hb/r_cut) + 1;
  nc = floor(hc/r_cut) + 1;
  double density = (double)n_atoms/na/nb/nc;
  // Note here +N is the magic number, we assume the density fluctuation is small enough
  // so +N is large enough buffer.
  n_max = floor(density) + 10;
  // create celllist, initialize with -1
  // -1 signals the end of the list
  celllist.resize(na, nb, nc, n_max);
  celllist.setConstant(-1);

  MatrixXd cell=box, cell_inv;
  cell.row(0) /= na;
  cell.row(1) /= nb;
  cell.row(2) /= nc;
  cell_inv = cell.inverse();

  // cell indices of each particle
  MatrixXd spositions = positions * cell_inv;
  MatrixXi indices = spositions.array().floor().cast<int>();
  for (int i_atom=0; i_atom<n_atoms; i_atom++) {
    // pbc shifts
    int label = atom_labels(i_atom);
    int ia = indices(i_atom, 0) % na;
    int ib = indices(i_atom, 1) % nb;
    int ic = indices(i_atom, 2) % nc;
    if (indices(i_atom, 0) < 0) ia += na;
    if (indices(i_atom, 1) < 0) ib += nb;
    if (indices(i_atom, 2) < 0) ic += nc;
    // assign atom to a particular cell
    for (int i=0; i<n_max; i++) {
      if (celllist(ia, ib, ic, i) < 0) {
        celllist(ia, ib, ic, i) = label;
        break;
      }
    }
  }

  return celllist;
}

// pybind11 does not recognize the multidimensional Tensor datatype, hence the beaucracy
ArrayXi construct_cell_list_pywrap(const MatrixXd & positions, const ArrayXi & atom_labels,
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut)
{
  Tensor<int, 4, RowMajor> celllist;
  celllist = construct_cell_list(positions, atom_labels, n_atoms, box, box_inv, r_cut);
  ArrayXi results;
  Map<ArrayXi> data(celllist.data(), celllist.size());
  // the first four elements are dimensions
  results.resize(celllist.size() + 4);
  results(0) = celllist.dimension(0);
  results(1) = celllist.dimension(1);
  results(2) = celllist.dimension(2);
  results(3) = celllist.dimension(3);
  results.segment(4, celllist.size()) = data;
  return results;
}

PYBIND11_MODULE(neighbour_list, m)
{
  m.doc() = "Neighbour cell list module";
  m.def("construct_cell_list_pywrap", &construct_cell_list_pywrap);
}
