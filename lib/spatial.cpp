/*
 * =====================================================================================
 *
 *       Filename:  spatial.cpp
 *
 *    Description:  library taking care of maths related to spatial geometries
 *
 *        Version:  1.0
 *        Created:  12/26/2018 02:17:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (Kuang Yu), 
 *   Organization:  TBSI
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <mylibs.h>
#include <stdlib.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

using namespace std;
using namespace Eigen;

// find out the vector point from r1->r2, under pbc
// flag == 0: feed in real coordinates
// flag == 1: feed in direct coordinates
// if flag == 1: then the box_inv value does not matter
VectorXd dr_vec_pbc(const VectorXd & r1, const VectorXd & r2, const MatrixXd & box, const MatrixXd & box_inv, int flag)
{
  VectorXd dr, ds;
  if (flag==0){
    dr = r2 - r1;
    ds = box_inv * dr;
  }
  else {
    ds = r2 - r1;
  }
  // pbc shift
  ds -= (ds.array()+0.5).floor().matrix();
  //for (int dd=0; dd<3; dd++) ds(dd) -= floor(ds(dd));
  return (box * ds);
}

// return distance between r1 & r2, under pbc
double distance_pbc(const VectorXd & r1, const VectorXd & r2, const MatrixXd & box, const MatrixXd & box_inv, int flag)
{
  VectorXd dr;
  dr = dr_vec_pbc(r1, r2, box, box_inv, flag);
  return dr.norm();
}

// compute the interatomic distances in batch
// pos1 is a 3*na1 position matirx
// pos2 is a 3*na2 position matrix
// calculate the interatomic distances between these two
// results is a na2*na1 distance^2 matrix
// padding structure:
// pos1_pad:
// [r1_1, r1_1, ... r1_1  | r1_2, r1_2, ... r1_2  | ... | r1_n1, r1_n1, ... r1_n1]
// pos2_pad:
// [r2_1, r2_2, ... r2_n2 | r2_1, r2_2, ... r2_n2 | ... | r2_1,  r2_2, ...  r2_n2]
MatrixXd distance2_pbc_batch(const MatrixXd & pos1, const MatrixXd & pos2, 
    const MatrixXd & box, const MatrixXd & box_inv, int flag)
{
  int na1 = pos1.cols();
  int na2 = pos2.cols();
  MatrixXd pos1_pad, pos2_pad;
  pos1_pad.resize(3, na1*na2);
  pos2_pad.resize(3, na1*na2);
  // pad pos1
  for (int i=0; i<na1; i++) {
    pos1_pad.block(0, i*na2, 3, na2) = pos1.col(i).replicate(1, na2);
  }
  // pad pos2
  for (int i=0; i<na1; i++) {
    pos2_pad.block(0, i*na2, 3, na2) = pos2;
  }
  // the eigen native replicate function: really slow, why???
  //pos2_pad = pos2.replicate(1, na1);
  MatrixXd dr_vecs, ds_vecs;
  if (flag==0) {
    dr_vecs = pos2_pad - pos1_pad;
    ds_vecs = box_inv * dr_vecs;
  }
  else {
    ds_vecs = pos2_pad - pos1_pad;
  }
  // pbc shifts
  ds_vecs -= (ds_vecs.array()+0.5).floor().matrix();
  dr_vecs = box * ds_vecs;
  // find out distances
  VectorXd distances2;
  distances2 = dr_vecs.colwise().squaredNorm();
  // cast results into matrix
  Map<MatrixXd> results(distances2.data(), na2, na1);

  // bypass the main body of the function
//  MatrixXd results2(na2, na1);
//  results2.setConstant(100.0);
  
  return results;
}

// return angle between r1 & r2, in radian
double angle_vec(const VectorXd & r1, const VectorXd & r2)
{
  float cos_A = r1.dot(r2) / r1.norm() / r2.norm();
  return acos(cos_A);
}

// return 1-2-3 angle in radian values
double angle_pbc(const VectorXd & r1, const VectorXd & r2, const VectorXd & r3, const MatrixXd & box, const MatrixXd & box_inv, int flag)
{
  Vector3d dr21, dr23;
  dr21 = dr_vec_pbc(r2, r1, box, box_inv, flag);
  dr23 = dr_vec_pbc(r2, r3, box, box_inv, flag);
  return angle_vec(dr21, dr23);
}

// return 1-2-3-4 dihedral in radian values
double dihedral_pbc(const VectorXd & r1, const VectorXd & r2, const VectorXd & r3,
    const VectorXd & r4, const MatrixXd & box, const MatrixXd & box_inv, int flag)
{
  Vector3d dr21, dr23, dr32, dr34;
  Vector3d n1, n2;
  dr21 = dr_vec_pbc(r2, r1, box, box_inv, flag);
  dr23 = dr_vec_pbc(r2, r3, box, box_inv, flag);
  dr32 = -dr23;
  dr34 = dr_vec_pbc(r3, r4, box, box_inv, flag);
  n1 = dr23.cross(dr21);
  n2 = dr34.cross(dr32);
  return angle_vec(n1, n2);
}

double volume(const MatrixXd & box) 
{
  return box.determinant();
}

PYBIND11_MODULE(spatial, m)
{
  m.doc() = "Spatial operations regarding 3d vectors";
  m.def("distance_pbc", &distance_pbc, "distance between two points under pbc.");
  m.def("angle_pbc", &angle_pbc, "angle between r1, r2, and r3 under pbc.");
}
