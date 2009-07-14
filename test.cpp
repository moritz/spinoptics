#include <iostream>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stdlib.h>
using namespace std;

using namespace Eigen;

typedef std::complex<double> cnum;
typedef Eigen::SparseMatrix< cnum , Eigen::RowMajor> esm;
typedef Eigen::RandomSetter< esm > ers;
typedef Eigen::SparseLU<esm, Eigen::SuperLU> eslu;


void pseudo_sparse_solve(const eslu * const slu,
                         const esm &rhs,
                         esm &result,
                         const bool adjoint = false) {
    assert( rhs.cols() == rhs.rows() );
    int n = rhs.cols();

    ers setter(result);
    Eigen::VectorXcd *invCol = new Eigen::VectorXcd(n);
    int transpose_flag;
    if (adjoint) {
        transpose_flag = Eigen::SvAdjoint;
    } else {
        transpose_flag = Eigen::SvNoTrans;
    }
    for (int k=0; k<rhs.outerSize(); ++k) {
        Eigen::VectorXcd base(n);
        int i = 0;
        for (esm::InnerIterator it(rhs,k); it; ++it) {
            base(it.col()) = it.value();
            i++;
        }
        if (i != 0) {
            slu->solve(base, invCol, transpose_flag);
            for (int j = 0; j < n; j++) {
                if (abs((*invCol)[j]) > 1e-18) {
                    setter(j, k) = (*invCol)[j];
                }
            }
        }
    }
    delete invCol;
}


int main(int argc, char** argv) {
    srand(1);
    esm a(100, 100);
    {
        ers s(a);
        for (int x = 0; x < 100; x++) {
            for (int y = 0; y < 100; y++) {
                if (rand() < RAND_MAX / 2)
                    s(x, y) = 3;
            }
        }
    }
    esm b(100, 100);
    eslu slu(a);
    pseudo_sparse_solve(&slu, a, b, false);

}
 // vim: ts=4 sw=4 expandtab
