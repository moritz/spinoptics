#include <iostream>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <stdlib.h>
#include "invert-matrix.hpp"
#include "sparse_io.hpp"
#include <getopt.h>
#include <assert.h>
#include <time.h>
// for getpid();
#include <sys/types.h>
#include <unistd.h>
// Eigen libs

#include <Eigen/Core>
#include <Eigen/Sparse>


#define IDX(x, y) ((x) + Nx * (y))

using namespace boost::numeric::ublas;
using namespace std;
USING_PART_OF_NAMESPACE_EIGEN

typedef double num;
typedef complex<num> cnum;
typedef matrix<cnum> cmatrix;
typedef compressed_matrix<cnum, row_major> sparse_cm;
typedef unsigned int idx_t;

int Nx               = 10;
int Ny               = Nx;
const int N_leads    = 8;

const num epsilon    = 1e-5;

const num pi         = 3.14159265358979323846264;
const num h_bar      = 6.582122E-16;     // [eV*s]
const num h_planck   = 4.135669E-15;     // [eV*s]
const num electron_mass
                     = 9.109389E-31;     // [kg]
const num mass       = 0.381 * electron_mass;
const num e_charge   = 1.60217653E-19;   // [C = A*s]
const num bohr_magneton
                     = 9.27400915E-24;   // [A * m^2]

const num width_lead = 30.0;             // [nm]
const num g_factor   = 20.0;

num global_gauge     = 1.0;

// XXX is the +1 correct?
const num a_lead     = width_lead / (num) (Nx + 1);
const int size       = Nx * Ny * 2;      // `Nfin'
const num V          = 1.0;              // hopping term

// e_tot is our choice of energy zero-level.
const num e_tot      = -2.0 * V;
const num width_disorder  = 0.0;

num alpha = -0.02 / a_lead / 2.0;

void log_tick(const char* desc, bool end = false) {
    static time_t prev = time(NULL);
    time_t t = time(NULL);
    printf("[Tick] %06ld %s\n", t-prev, desc);
    prev = t;
}

num flux_from_field(const num B) {
    return 2.0 * pi * B 
        * (a_lead * a_lead) * 1e-18 // a_lead is in nm
        / h_planck; 
}

cnum b_factor(const num flux, const int n) {
    return exp(cnum(0.0, -1) * flux * (num) n);
}

void correct_phase(sparse_cm &m, num flux) {
    sparse_cm::iterator1 x = m.begin1();
    sparse_cm::iterator1 x_end = m.end1();
    for (; x != x_end; ++x) {
        sparse_cm::iterator2 y = x.begin();
        sparse_cm::iterator2 y_end = x.end();
        for (; y != y_end; y++) {
            int x1 = y.index1() % Nx;
            int y1 = (y.index1()/Nx) % Ny;

            int x2 = y.index2() % Nx;
            int y2 = (y.index2()/Nx) % Ny;

            cnum phi = b_factor(x2 * y2 - x1 * y1, flux);
            (*y) *= phi;
        }
    }

}

template <class T>
void set_zero(matrix<T> &m) {
    for (unsigned int x = 0; x < m.size1(); x++){
        for (unsigned int y = 0; y < m.size2(); y++){
            m(x, y) = 0.0;
        }
    }
}

template <class T>
idx_t count_nonzero(matrix<T> &m) {
    idx_t i = 0;
    for (idx_t x = 0; x < m.size1(); x++){
        for (idx_t y = 0; y < m.size2(); y++){
            if (m(x, y) != 0.0) 
                i++;
        }
    }
    return i;
}

idx_t count_nonzero(sparse_cm &m) {
    idx_t i = 0;
    sparse_cm::const_iterator1 x = m.begin1();
    sparse_cm::const_iterator1 x_end = m.end1();
    for (; x != x_end; ++x) {
        sparse_cm::const_iterator2 y = x.begin();
        sparse_cm::const_iterator2 y_end = x.end();
        for (; y != y_end; y++) {
            if ( (*y) != (num) 0.0) 
                i++;
        }
    }

    return i;
}

void sparse_inverse(Eigen::SparseMatrix< cnum > &m, cmatrix &inv) {
    Eigen::SparseLU<Eigen::SparseMatrix< cnum >,Eigen::SuperLU> slu(m);
    Eigen::VectorXcd base(size), invCol(size);
    for (int i=0; i<size; ++i) {
        base = Eigen::VectorXcd::Unit(size, i);
        slu.solve(base, &invCol);
        for (int j=0; j<size; ++j) {
            inv(j, i) = invCol[j];
//            (*inv)(j,i) = invCol[j];
        }
    }
//    return inv;
}

// sparse_product(a, b, c) computes the matrix product
// a * b where a is a sparse matrix and b is a full one,
// and assumes that c has a different storage location than
// a and b
void sparse_product(const sparse_cm &m1, const cmatrix &m2, sparse_cm &r) {
    sparse_cm::const_iterator1 x = m1.begin1();
    sparse_cm::const_iterator1 x_end = m1.end1();
    r.clear();
    idx_t s = m2.size1();
    for (; x != x_end; ++x) {
        sparse_cm::const_iterator2 y = x.begin();
        sparse_cm::const_iterator2 y_end = x.end();
        for (; y != y_end; y++) {
            for (idx_t j = 0; j < s; j++){
                idx_t xi = y.index1();
                idx_t yi = y.index2();
                r(xi, j) += (*y) * m2(yi, j);
            }
        }
    }
}

// sparse_herm_product(a, b, c) computes the matrix product
// a * herm(b) where a is a sparse matrix and b is a full one,
// (herm(b) is the hermitian conjugate, ie complex conjugation
// and transposition).
// It assumes that c has a different storage location than
// a and b
void sparse_herm_product(const sparse_cm &m1, const cmatrix &m2, sparse_cm &r) {
    sparse_cm::const_iterator1 x = m1.begin1();
    sparse_cm::const_iterator1 x_end = m1.end1();
    r.clear();
    idx_t s = m2.size1();
    for (; x != x_end; ++x) {
        sparse_cm::const_iterator2 y = x.begin();
        sparse_cm::const_iterator2 y_end = x.end();
        for (; y != y_end; y++) {
            for (idx_t j = 0; j < s; j++){
                idx_t xi = y.index1();
                idx_t yi = y.index2();
                r(xi, j) += (*y) * conj(m2(j, yi));
            }
        }
    }
}

inline num rashba(const num alpha) {
    return 2.0 * alpha * a_lead;
}

inline num mods(const int n, const int nle) {
    return 2.0 * V * (cos(pi * (num) (n+1) / ((num) nle + 1.0)) 
            - 1.0);
}

cnum findk(const num Emod) {
    return sqrt(
            cnum(2 * mass / (h_bar * h_bar)
            * (e_tot - Emod) * 10.0 / 1.60219, 0)
        );
}

sparse_cm* hamiltonian(const num rashb, const num B) {
    sparse_cm* Hnn = new sparse_cm(size, size, 4 * size);

    // division by e_charge to convert from electron volt to Joule
    num zeeman = 0.5 * g_factor * bohr_magneton * B / e_charge;
//    num zeeman = 0.0;
    num flux = flux_from_field(B);
    num gauge = 1.0;
    num xflux = gauge * flux;
    num yflux = (1.0 - gauge) * flux;
    cout << "Zeeman term: " << zeeman << endl;

    // diagonal elements
    for (int i = 0; i < size / 2; i++) {
        // later we might want to add random disorder,
        // in which case these items might be different per
        // iteration, but every two diagonal items with distance
        // (size/2) must still have the same value
        cnum energy = 4.0 * V + e_tot;
        (*Hnn)(i, i)                     = energy - zeeman;
        (*Hnn)(i + size/2, i + size/2)   = energy + zeeman;
    }

    /*                Nx
     *       +----------------+
     *      /                  \
     *      X   X   X ...   X   X  -
     *      X   X   X ...   X   X    \
     *      .   .   .       .   .    |
     *      .   .   .       .   .    |  Ny
     *      .   .   .       .   .    |
     *      X   X   X ...   X   X    /
     *      X   X   X ...   X   X  -
     *
     *      Numbering scheme:
     *      0   1     ...       (Nx-1)
     *      Nx  Nx+1  ...       (2*Nx -1)
     *
     *      if the index i is given, 
     *      y = i / Nx
     *      x = i % Nx
     */

    // interaction in x direction
    int s = size / 2;

    cnum r = rashba(rashb);
    for (int x = 0; x < Nx - 1; x++) {
        for (int y = 0; y < Ny; y++) {
            cnum h = -V * conj(b_factor(xflux, y));
            // kinetic energy
            (*Hnn)(    IDX(x,   y),     IDX(x+1, y)) = h;
            (*Hnn)(    IDX(x+1, y),     IDX(x,   y)) = conj(h);
            (*Hnn)(s + IDX(x,   y), s + IDX(x+1, y)) = h;
            (*Hnn)(s + IDX(x+1, y), s + IDX(x,   y)) = conj(h);
            // Rashba terms
            // "1 and 102"
            // with spin flip
            cnum b = b_factor(xflux, y);
            (*Hnn)(IDX(x, y)      , IDX(x+1, y) + s) = -r * conj(b);
            (*Hnn)(IDX(x+1, y) + s, IDX(x,  y)     ) = -r * b;
            // "101 and 2"
            (*Hnn)(IDX(x+1, y)    , IDX(x, y) + s  ) = r * b;
            (*Hnn)(IDX(x, y) + s  , IDX(x+1, y)    ) = r * conj(b);
        }
    }

    // kinetic energy in y direction
    // spin up

    for (int x = 0; x < Nx; x++){
        for (int y = 0; y < Ny - 1; y++) {
            cnum b = b_factor(yflux, x);
            cnum h = -V * b;
            (*Hnn)(IDX(x, y)      , IDX(x, y+1)    ) = h;
            (*Hnn)(IDX(x, y+1)    , IDX(x, y  )    ) = conj(h);
            (*Hnn)(IDX(x, y)   + s, IDX(x, y+1) + s) = h;
            (*Hnn)(IDX(x, y+1) + s, IDX(x, y  ) + s) = conj(h);
            // Rashba terms
            // "11 and 101"
            h = cnum(0, 1) * b * r;
            (*Hnn)(IDX(x, y+1)    , IDX(x, y  ) + s) = conj(h);
            (*Hnn)(IDX(x, y  ) + s, IDX(x, y+1)    ) = h;
            // "1 and 111"
            (*Hnn)(IDX(x, y  )    , IDX(x, y+1) + s) = h;
            (*Hnn)(IDX(x, y+1) + s, IDX(x, y  )    ) = conj(h);
        }
    }

    log_tick("hamiltonian");
    return Hnn;
};

sparse_cm** self_energy(num flux, num gauge) {
    // Glp1lp1n = G_{l+1, l+1}n
    cmatrix Glp1lp1n = cmatrix(Nx, Nx);
    Glp1lp1n.clear();
    for (int r = 0; r < Nx; r++) {
        num x = (e_tot - mods(r, Nx))
            / (2.0 * V) + 1.0;
        cnum theta;
        if (x > 1.0) {
            // evanescent mode, calculate cosh^-1
            theta = cnum(0, 1)
                * log((cnum) (x + sqrt(x*x - 1.0)));
        } else if (x < -1.0) {
            theta = cnum(0, 1)
                * log((cnum) (x - sqrt(x*x - 1.0)));
        } else {
            theta = acos(x);
        }

        cnum unit = cnum(1.0, 0.0);
        num  f    = (num) (Nx + 1);
        cnum tmpp = exp(cnum(0, 1) * (num) (2.0 * pi 
                    * (((num) ((r+1) * Nx )) / f)));
        cnum tmpm = exp(cnum(0, -1) *(num)  (2.0 * pi 
                    * ((num) (r+1) / f)));

        // "AnorN1(mm)" in nano0903c.f
        num y = 1.0 / sqrt(0.5 * Nx +
            real((unit - tmpp) / (unit-tmpp) * cnum(0.5, 0.0)));

        for (int p = 0; p < Nx; p++) {
            // "psiN1(ii)" in nano0903c.f
            cnum y1 = y * sin(pi * (num) ((p+1) * (r+1))/(1.0 + Nx));
            for (int q = 0; q < Nx; q++) {
                // "psiN1(jj)" in nano0903c.f
                cnum y2 = y * sin(pi * (num) ((q+1) * (r+1))/(1.0 + Nx));
                Glp1lp1n(p, q) += exp(cnum(0.0,1.0) * theta)/V * y1 * y2;
            }
        }
    }

    int s = size / 2;

    sparse_cm *G_left_up     = new sparse_cm(size, size, s);
    sparse_cm *G_left_down   = new sparse_cm(size, size, s);
    sparse_cm *G_right_up    = new sparse_cm(size, size, s);
    sparse_cm *G_right_down  = new sparse_cm(size, size, s);

    sparse_cm *G_top_up      = new sparse_cm(size, size, s);
    sparse_cm *G_top_down    = new sparse_cm(size, size, s);
    sparse_cm *G_bottom_up   = new sparse_cm(size, size, s);
    sparse_cm *G_bottom_down = new sparse_cm(size, size, s);


    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            cnum g = Glp1lp1n(i, j);
            (*G_left_up)    (IDX(0, i)       , IDX(0, j)       ) = g;
            (*G_left_down)  (IDX(0, i) + s   , IDX(0, j) + s   ) = g;

            (*G_right_up)   (IDX(Nx-1, i)    , IDX(Nx-1, j)    ) = g;
            (*G_right_down) (IDX(Nx-1, i) + s, IDX(Nx-1, j) + s) = g;

            (*G_top_up)     (IDX(i, 0)       , IDX(j, 0)       ) = g;
            (*G_top_down)   (IDX(i, 0) + s   , IDX(j, 0) + s   ) = g;

            (*G_bottom_up)  (IDX(i, Ny-1)    , IDX(j, Ny-1)    ) = g;
            (*G_bottom_down)(IDX(i, Ny-1) + s, IDX(j, Ny-1) + s) = g;
        }
    }

    sparse_cm** sr = new sparse_cm*[N_leads];
    sr[0] = G_left_up;
    sr[2] = G_left_down;

    sr[4] = G_top_up;
    // 5, not 6
    sr[5] = G_top_down;

    sr[1] = G_right_up;
    sr[3] = G_right_down;
    // 6, not 5
    sr[6] = G_bottom_up;
    sr[7] = G_bottom_down;

    for (int i = 0; i < N_leads; i++) {
        switch(i) {
            case 0:
            case 1:
            case 2:
            case 3:
                correct_phase(*sr[i], -(1.0 - gauge) * flux);
                break;
            case 4:
            case 5:
            case 6:
            case 7:
                cout << "gauge scaling factor: " << flux << endl;
                correct_phase(*sr[i], gauge * flux);
                break;
        }
    }

    log_tick("self-energy");
    return sr;
}

matrix<num>* greenji(sparse_cm* Hnn, num flux, num gauge) {
    sparse_cm **sigma_r   = self_energy(flux, gauge);

    for (int k = 0; k < N_leads; k++){
        noalias(*Hnn) -= *(sigma_r[k]);
    }
    MatrixXcd I = MatrixXcd::Identity(size, size), e_green(size, size);
    Eigen::SparseMatrix< cnum > e_green_inv(size, size);

    {
        Eigen::RandomSetter<Eigen::SparseMatrix< cnum > > setter(e_green_inv);
        sparse_cm::const_iterator1 x = Hnn->begin1();
        sparse_cm::const_iterator1 x_end = Hnn->end1();
        for (; x != x_end; ++x) {
            sparse_cm::const_iterator2 y = x.begin();
            sparse_cm::const_iterator2 y_end = x.end();
            for (; y != y_end; y++) {
                setter(y.index1(), y.index2()) = *y;
            }
        }
    }
    log_tick("green_inv");
    cmatrix *green = new cmatrix(size, size);
    sparse_inverse(e_green_inv, *green);
//    cout << "Green: " << *green << endl;
//    cout << green->size1() * green->size2() << "\t";
//    cout << count_nonzero(*green) << endl;
    log_tick("matrix inversion");


    matrix<num> *tpq = new matrix<num>(N_leads, N_leads);
    tpq->clear();
    sparse_cm **gamma_g_adv = new sparse_cm*[N_leads];
    sparse_cm **gamma_g_ret = new sparse_cm*[N_leads];

    // T_{p, q} = Trace( \Gamma_p * G^R * \Gamma_q * G^A )
    // where G^A = (G^R)^\dagger
    //
    //
    // \Gamma = -2 * Im(\Sigma_r)   
    // (calculate slices of \Gamma (aka gamm_i) on the fly
    // to save memory
    // first carry out the two products 
    // \Gamma_p * G^R and \Gamma_q * G^A
    sparse_cm * gamm_i = new sparse_cm(size, size);
    for (int i = 0; i < N_leads; i++) {
        sparse_cm *g_adv = new sparse_cm(size, size);
        sparse_cm *g_ret = new sparse_cm(size, size);
        assert(V * V == 1);
        noalias(*gamm_i) = -2 * imag(*sigma_r[i]);
        delete sigma_r[i];
        sigma_r[i] = NULL;
        sparse_product     (*gamm_i, *green, *g_ret);
        sparse_herm_product(*gamm_i, *green, *g_adv);
        gamma_g_adv[i] = g_adv;
        gamma_g_ret[i] = g_ret;
    }
    log_tick("products");
    delete gamm_i;      gamm_i      = NULL;
    delete green;       green       = NULL;

    delete[] sigma_r;
    sigma_r = NULL;

    // Now calculate the trace
    for (int i = 0; i < N_leads; i++){
        for (int j = 0; j < N_leads; j++){
            for (int n = 0; n < size; n++){
                for (int m = 0; m < size; m++){
                    cnum x = (*(gamma_g_ret[i]))(n, m);
                    cnum y = (*gamma_g_adv[j])(m, n);
                    (*tpq)(i, j) += real(x * y);
                }
            }
        }
    }
    log_tick("trace");

    for (int i = 0; i < N_leads; i++){
        for (int n = 0; n < size; n++){
            cnum x = (*gamma_g_ret[i])(n, n);
            cnum y = (*gamma_g_adv[i])(n, n);
            (*tpq)(i, i) += real(cnum(0, 1) * y - cnum(0, 1) * x);

        }
        delete gamma_g_adv[i];
        delete gamma_g_ret[i];
    }
    for (int i = 0; i < Nx; i++){
        cnum k = findk(mods(i, Nx));
        if (imag(k) == 0.0) {
            for (int j = 0; j < N_leads; j++){
                (*tpq)(j, j) += 1.0;
            }
        }

    }
    log_tick("tpq corrections");

    delete[] gamma_g_adv;
    delete[] gamma_g_ret;

    return tpq;
}


int main (int argc, char** argv) {
    log_tick("start");
    num Bz = 6;

    int opt;
    while ((opt = getopt(argc, argv, "r:b:s:")) != -1) {
       switch (opt) {
            case 'r':
                alpha = atof(optarg);
                break;
            case 's':
                Nx = atoi(optarg);
                Ny = Nx;
                break;
            case 'b':
                Bz = atof(optarg);
                break;
            default:
                cerr << "Error while processing command line args\n";
                exit(1);

        }
    }

    cout << "PID:     " << getpid() << endl;
    cout << "Size:    " << Nx << "x" << Ny << endl;
    cout << "Bz:      " << Bz << endl;
    sparse_cm *Hnn = hamiltonian(alpha, Bz);
#ifndef NDEBUG
    {
        sparse_cm HH = herm(*Hnn);
        int errors = 0;
        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++){
                cnum c1 = (*Hnn)(x, y);
                cnum c2 = HH(x, y);
                if (abs(c1 - c2 ) > 1e-5) {
                    errors++;
                }
            }
        }
        if (errors > 0) {
            cerr << "ERROR: Hamiltonian is not hermitian ("
                << errors << " differences)\n";
            exit(1);
        }
    }
#endif

//    cout << "Hamiltonian: " << *Hnn << "\n";
//    cout << io::sparse(*Hnn) << endl;
    num flux = flux_from_field(Bz);
    matrix<num> *tpq = greenji(Hnn, flux, global_gauge);
    delete Hnn;
    cout << "final tpq" << *tpq << endl;
    boost::numeric::ublas::vector<num> r;
    boost::numeric::ublas::vector<num> c;
    bool is_first = true;
    num ref = 0.0;
    for (int i = 0; i < N_leads; i++) {
        num r_sum = 0.0;
        num c_sum = 0.0;
        r = row(*tpq, i);
        c = column(*tpq, i);
        for (int j = 0; j < N_leads; j++) {
            r_sum += r(j);
            c_sum += c(j);
        }
        if (is_first) {
            ref = r_sum;
            is_first = false;
        }
        if (abs(r_sum - ref) > epsilon) {
            cout << "ERROR: sum rule violated for row " << i << endl;
        }
        if (abs(c_sum - ref) > epsilon) {
            cout << "ERROR: sum rule violated for column " << i << endl;
        }
    }
    delete tpq;
}

// vim: ft=cpp sw=4 ts=4 expandtab
