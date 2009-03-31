#include <iostream>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace std;

USING_PART_OF_NAMESPACE_EIGEN

typedef double num;
typedef complex<num> cnum;
typedef Eigen::SparseMatrix< cnum , Eigen::RowMajor> esm;

int main(int argc, char** argv) {
	MatrixXcd foo(3, 3);
	foo(0, 2) = cnum(2, 3);
	cout << foo << endl;

	foo = foo.adjoint().eval();
	cout << "adjoint: \n";

	cout << foo << endl;
}
