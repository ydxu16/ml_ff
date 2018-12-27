/*
 * =====================================================================================
 *
 *       Filename:  neighbour_list.h
 *
 *    Description:  header for neighbour_list.cpp
 *
 *        Version:  1.0
 *        Created:  12/26/2018 02:51:04 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

using namespace std;
using namespace Eigen;


Tensor<int, 4> construct_cell_list(const MatrixXd & positions, const ArrayXi & atom_labels,
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut);


MatrixXi build_cell_index_map(const Tensor<int, 4, RowMajor> & celllist, int n_atoms);


ArrayXi find_neighbours(const MatrixXd & positions, const Tensor<int, 4, RowMajor> & celllist,
    int i_atom, int i0, int j0, int k0, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);


ArrayXi find_neighbours_for_atom(const MatrixXd & positions, const Tensor<int, 4, RowMajor> & celllist,
    int i_atom0, const MatrixXi & cell_index_map, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);


MatrixXi find_neighbours_for_all(const MatrixXd & positions, const Tensor<int, 4, RowMajor> & celllist,
    int n_atoms, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);


/* python wrapper beaucracy  */
ArrayXi Tensor_4d_to_1d(const Tensor<int, 4, RowMajor> & t_4d);

Tensor<int, 4, RowMajor> Tensor_1d_to_4d(const ArrayXi & A_1d);

ArrayXi construct_cell_list_py(const MatrixXd & positions, const ArrayXi & atom_labels,                                                                                  
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut);

MatrixXi build_cell_index_map_py(const ArrayXi * celllist_py, int n_atoms);

ArrayXi find_neighbours_for_atom_py(const MatrixXd & positions, const ArrayXi & celllist_py,
    int i_atom0, const MatrixXi & cell_index_map, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);

MatrixXi find_neighbours_for_all_py(const MatrixXd & positions, const ArrayXi & celllist_py,
    int n_atoms, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);
