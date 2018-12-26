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

ArrayXi construct_cell_list_pywrap(const MatrixXd & positions, const ArrayXi & atom_labels,                                                                                  
    int n_atoms, const MatrixXd & box, const MatrixXd & box_inv, double r_cut);
