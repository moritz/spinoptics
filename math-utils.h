#ifndef __MATH_UTILS_H
#define __MATH_UTILS_H

#include <complex>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include "zmumps_c.h"

namespace ub = boost::numeric::ublas;

typedef double num;
typedef std::complex<num> cnum;
typedef ub::matrix<cnum> cmatrix;
typedef ub::compressed_matrix<cnum, ub::row_major> sparse_cm;
typedef Eigen::SparseMatrix< cnum , Eigen::RowMajor> esm;
typedef Eigen::RandomSetter< esm > ers;
typedef Eigen::SparseLU<esm, Eigen::SuperLU> eslu;

// calculate the (dense) inverse of a sparse, complex matrix
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

// convert a sparse boost::numeric::ublas to an Eigen::SparseMatrix
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

// convert an Eigen::SparseMatrix to a boost::numeric::ublas matrix
void eigen_to_ublas(const esm &m, sparse_cm &result) {
    result.clear();
    for (int k=0; k<m.outerSize(); ++k) {
        for (esm::InnerIterator it(m,k); it; ++it) {
            result(it.row(), it.col()) = it.value();
        }
    }
}

// given a sparse LU decompostion slu of a matrix A and a right-hand side
// rhs, solve the equation A * x = rhs or A^\dagger * x = rhs (if the 
// `adjoint' flag is true
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

// real(trace(a * b))
num r_prod_trace(const esm &a, const esm &b) {
    num sum = 0.0;
    cnum null = cnum(0, 0);
    for (int k = 0; k < a.outerSize(); ++k) {
        for (esm::InnerIterator it(a,k); it; ++it) {
            cnum y = b.coeff(it.col(), it.row());
            sum += real(it.value() * y);
        }
    }

    return sum;
}

typedef struct {
    int             size;
    int             nz; // non-zeros
    int             *row_ptr;
    int             *col_ptr;
    ZMUMPS_COMPLEX  *values;
} MUMPS_complex_sparse;

MUMPS_complex_sparse* eigen_to_mumps(const esm &e) {
    assert(e.rows() == e.cols());
    int nz = e.nonZeros();
    MUMPS_complex_sparse *m = new MUMPS_complex_sparse;
    m->size     = e.rows();
    m->nz       = nz;
    m->row_ptr  = new int[nz];
    m->col_ptr  = new int[nz];
    m->values   = new ZMUMPS_COMPLEX[nz];
    int i = 0 ;
    for (int k=0; k < e.outerSize(); ++k) {
        for (esm::InnerIterator it(e,k); it; ++it) {
            m->row_ptr[i]   = 1 + it.row();
            m->col_ptr[i]   = 1 + it.col();
            m->values[i].r  = real(it.value());
            m->values[i].i  = imag(it.value());
        }
    }
    return m;
}

ZMUMPS_STRUC_C* MUMPS_lr(const esm &e) {
    assert(e.rows() == e.cols());
    ZMUMPS_STRUC_C *id = new ZMUMPS_STRUC_C;
    id->job = -1;           // initialize
    id->par = 1;            // just one job in parallel
    id->sym = 0;            // no symmetries
    id->comm_fortran = -987654; // black magic
//    zmumps_c(id);           // init
    MUMPS_complex_sparse *m = eigen_to_mumps(e);
    id->n   = m->size;
    id->nz  = m->nz;
    id->irn = m->row_ptr;
    id->jcn = m->col_ptr;

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    id->ICNTL(1) = -1;      // be silent
    id->ICNTL(2) = -1;      // be silent
    id->ICNTL(3) = -1;      // be silent
    id->ICNTL(4) =  0;      // be silent

    id->job      = 4;       // analyze + decompose
//    zmumps_c(id);
    return id;
}

void MUMPS_solve(ZMUMPS_STRUC_C *id, const esm &e, esm &result) {
    assert(id->n    == e.rows());
    assert(e.rows() == e.cols());
    assert(e.rows() == result.rows());
    assert(e.rows() == result.cols());

}

#endif /* __MATH_UTILS_H */
