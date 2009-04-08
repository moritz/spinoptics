#include <iostream>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace std;

using namespace Eigen;

typedef double num;
typedef complex<num> cnum;
typedef Eigen::SparseMatrix< cnum , Eigen::RowMajor> esm;
typedef Eigen::RandomSetter< esm > ers;
typedef Eigen::SparseLU<esm, Eigen::SuperLU> eslu;

int main(int argc, char** argv) {
    esm a(2, 2);
    {
        ers s(a);
        s(0, 0) = cnum(2, 0);
        s(0, 1) = cnum(1, 1);
        s(1, 0) = cnum(2, 1);
    }
    {
        eslu s(a);
        VectorXcd rhs(2);
        VectorXcd solution(2);
        rhs(0) = cnum(1, 0);
        rhs(1) = cnum(2, 0);

        s.solve(rhs, &solution);
        cout << solution << "\n\n";
        cout << a * solution << "\n\n";

        s.solve(rhs, &solution, SvTranspose);
        cout << solution << "\n\n";
        cout << a.transpose() * solution << "\n\n";

        s.solve(rhs, &solution, SvAdjoint);
        cout << solution << "\n\n";
        cout << a.adjoint() * solution << "\n\n";
    }
}
