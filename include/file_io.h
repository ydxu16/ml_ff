#include <string>

using namespace std;

void read_dmatrix2d_from_file(string ifn, int * size, Eigen::MatrixXd & matrix);

template <typename Vectype>
void read_vector_from_file(string ifn, int * size, Vectype & vector)
{
  string line;
  ifstream ifile;
  ifile.open(ifn, ios::in);
  if (ifile.is_open())
  {
    getline(ifile, line);
    stringstream stream(line);
    stream >> *size;
    vector.resize(*size);
    for (int i=0; i<(*size); i++)
    {
      getline(ifile, line);
      stringstream stream(line);
      stream >> vector(i);
    }
  }
  else
  {
    cout << "ERROR: can't open file: " << ifn << endl;
  }
  ifile.close();
}
