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
  return (ds.transpose() * box);
}

// return distance between r1 & r2, under pbc
double distance_pbc(const VectorXd & r1, const VectorXd & r2, const MatrixXd & box, const MatrixXd & box_inv, int flag)
{
  VectorXd dr;
  dr = dr_vec_pbc(r1, r2, box, box_inv, flag);
  return dr.norm();
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
