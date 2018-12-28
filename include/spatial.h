#include <string>
using namespace std;
using namespace Eigen;

VectorXd dr_vec_pbc(const VectorXd & r1, const VectorXd & r2, const MatrixXd & box, const MatrixXd & box_inv, int flag);

double distance_pbc(const VectorXd & r1, const VectorXd & r2, const MatrixXd & box, const MatrixXd & box_inv, int flag);

double angle_vec(const VectorXd & r1, const VectorXd & r2);

double angle_pbc(const VectorXd & r1, const VectorXd & r2, const VectorXd & r3, const MatrixXd & box, const MatrixXd & box_inv, int flag);

double dihedral_pbc(const VectorXd & r1, const VectorXd & r2, const VectorXd & r3,
    const VectorXd & r4, const MatrixXd & box, const MatrixXd & box_inv, int flag);

double volume(const MatrixXd & box); 
