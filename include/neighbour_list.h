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


ArrayXi ind1d_to_ind3d(int n, int na, int nb, int nc);


int ind3d_to_ind1d(int i, int j, int k, int na, int nb, int nc);


bool is_neighbour_cell(ArrayXi ind0, ArrayXi ind1, int n1, int n2, int n3);


Tensor<int, 4> construct_cell_list(const MatrixXd & positions, const ArrayXi & atom_labels,
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut);


MatrixXi build_cell_index_map(const Tensor<int, 4> & celllist, int n_atoms);


ArrayXi find_neighbours(const MatrixXd & positions, const Tensor<int, 4> & celllist,
    int i_atom, const VectorXi & ind0, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);


ArrayXi find_neighbours_for_atom(const MatrixXd & positions, const Tensor<int, 4> & celllist,
    int i_atom0, const MatrixXi & cell_index_map, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);


MatrixXi find_neighbours_for_all(const MatrixXd & positions, const Tensor<int, 4> & celllist,
    int n_atoms, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);


/* python wrapper beaucracy  */
ArrayXi Tensor_4d_to_1d(const Tensor<int, 4> & t_4d);

Tensor<int, 4> Tensor_1d_to_4d(const ArrayXi & A_1d);

ArrayXi construct_cell_list_py(const MatrixXd & positions, const ArrayXi & atom_labels,                                                                                  
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut);

MatrixXi build_cell_index_map_py(const ArrayXi & celllist_py, int n_atoms);

ArrayXi find_neighbours_for_atom_py(const MatrixXd & positions, const ArrayXi & celllist_py,
    int i_atom0, const MatrixXi & cell_index_map, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);

MatrixXi find_neighbours_for_all_py(const MatrixXd & positions, const ArrayXi & celllist_py,
    int n_atoms, double r_cut, const MatrixXd & box, const MatrixXd & box_inv);
//
