#ifndef __MATH_UTILS_H
#define __MATH_UTILS_H

#include <complex>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>

namespace ub = boost::numeric::ublas;

typedef double num;
typedef std::complex<num> cnum;
typedef ub::matrix<cnum> cmatrix;
typedef ub::compressed_matrix<cnum, ub::row_major> sparse_cm;
typedef Eigen::SparseMatrix< cnum , Eigen::RowMajor> esm;
typedef Eigen::RandomSetter< esm > ers;
typedef Eigen::SparseLU<esm, Eigen::SuperLU> eslu;

void sparse_inverse(const esm &m, cmatrix &inv) {
    assert(m.rows() == m.cols());
    int size = m.rows();
    eslu slu(m);
    Eigen::VectorXcd base(size), invCol(size);
    for (int i=0; i<size; ++i) {
        base = Eigen::VectorXcd::Unit(size, i);
        slu.solve(base, &invCol);
        for (int j=0; j<size; ++j) {
            inv(j, i) = invCol[j];
        }
    }
}

void ublas_to_eigen(const sparse_cm &m, esm &result) {
    result.setZero();
    ers setter(result);
    sparse_cm::const_iterator1 x = m.begin1();
    sparse_cm::const_iterator1 x_end = m.end1();
    for (; x != x_end; ++x) {
        sparse_cm::const_iterator2 y = x.begin();
        sparse_cm::const_iterator2 y_end = x.end();
        for (; y != y_end; y++) {
            setter(y.index1(), y.index2()) = *y;
        }
    }
    return;
}

void eigen_to_ublas(const esm &m, sparse_cm &result) {
    result.clear();
    for (int k=0; k<m.outerSize(); ++k) {
        for (esm::InnerIterator it(m,k); it; ++it) {
            result(it.row(), it.col()) = it.value();
        }
    }
}

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


#endif /* __MATH_UTILS_H */
