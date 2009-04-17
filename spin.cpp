#include <iostream>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <stdlib.h>
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

namespace ub = boost::numeric::ublas;
using namespace std;

#include "math-utils.h"
typedef unsigned int idx_t;

const int Nx               = 10;
const int Ny               = Nx;
const int Spin_idx         = Nx * Ny;

const int N_leads    = 8;

// width of leads in units of lattice sites
const int lead_sites       = Nx;

/*   Numbering  scheme for the sites
 *
 *                Nx
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
 *      The same is repeated once more for spin down electrons
 *
 *      IDX(x, y, spin) returns an index into the site matrices
 *      (where `x' and `y' are zero based, and `spin' is either 0 or 1
 *
 *      X_IDX(i) returns the x index for site index i
 *      Y_IDX(i) returns the y index for site index i
 *      S_IDX(i) returns the spin index for site index i
 */

#define IDX(x, y, s) ((x) + Nx * (y) + (s) * Spin_idx)
#define X_IDX(i) (((i) % Spin_idx) % Nx)
#define Y_IDX(i) (((i) % Spin_idx) / Nx)
#define S_IDX(i) ((int) (i) / Spin_idx)

int lead_offset[N_leads];

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

const num width_sample
                     = 30.0;             // [nm]
const num g_factor   = 20.0;

num global_gauge     = 1.0;

// XXX is the +1 correct?
const num a_sample   = width_sample / (num) (Nx + 1);
const int size       = Nx * Ny * 2;      // `Nfin'
const num V          = 1.0;              // hopping term

// e_tot is our choice of energy zero-level.
const num e_tot      = -2.0 * V;
const num width_disorder  = 0.0;

num alpha = -0.02 / a_sample / 2.0;

void log_tick(const char* desc) {
    static time_t prev = time(NULL);
    time_t t = time(NULL);
    printf("[Tick] %06ld %s\n", t-prev, desc);
    prev = t;
}

num flux_from_field(const num B) {
    return 2.0 * pi * B
        * (a_sample * a_sample) * 1e-18 // a_sample is in nm
        / h_planck;
}

cnum b_factor(const num flux, const int n) {
    return exp(cnum(0.0, -1) * flux * (num) n);
}

void correct_phase(esm &m, const num flux) {
    if (flux == num(0))
        return;
//    cout << "correcting a phase...\n";
    for (int k = 0; k < m.outerSize(); ++k) {
        for (esm::InnerIterator it(m,k); it; ++it) {
            int x1 = X_IDX(it.row());
            int y1 = Y_IDX(it.row());
            assert(it.row() == IDX(x1, y1, S_IDX(it.row())));

            int x2 = X_IDX(it.col());
            int y2 = Y_IDX(it.col());
            assert(it.col() == IDX(x2, y2, S_IDX(it.col())));

//            cout << "product: " << x2 * y2 - x1 * y1 << endl;
            cnum phi = b_factor(flux, x2 * y2 - x1 * y1);
//            cout << "phi: " << phi << endl;
//            cout << "before: " << *y;
//            (*y) *= phi;
            it.value() *= phi;
//            cout << " after: " << *y << endl;
        }
    }
}

inline num rashba(const num alpha) {
    return 2.0 * alpha * a_sample;
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

esm* hamiltonian(const num rashb, const num B) {
    esm *H = new esm(size, size);
    ers Hnn( *H );

    // division by e_charge to convert from electron volt to Joule
    num zeeman = 0.5 * g_factor * bohr_magneton * B / e_charge;
//    num zeeman = 0.0;
    num flux = flux_from_field(B);
    num gauge = global_gauge;
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
        Hnn(i, i)                     = energy - zeeman;
        Hnn(i + size/2, i + size/2)   = energy + zeeman;
    }


    // interaction in x direction
    cnum r = rashba(rashb);
    for (int x = 0; x < Nx - 1; x++) {
        for (int y = 0; y < Ny; y++) {
            cnum h = -V * conj(b_factor(xflux, y));
            // kinetic energy
            (Hnn)(IDX(x,   y, 0), IDX(x+1, y, 0)) = h;
            (Hnn)(IDX(x+1, y, 0), IDX(x,   y, 0)) = conj(h);
            (Hnn)(IDX(x,   y, 1), IDX(x+1, y, 1)) = h;
            (Hnn)(IDX(x+1, y, 1), IDX(x,   y, 1)) = conj(h);
            // Rashba terms
            // "1 and 102"
            // with spin flip
            cnum b = b_factor(xflux, y);
            (Hnn)(IDX(x,   y, 0), IDX(x+1, y, 1)) = -r * conj(b);
            (Hnn)(IDX(x+1, y, 1), IDX(x,   y, 0)) = -r * b;
            // "101 and 2"
            (Hnn)(IDX(x+1, y, 0), IDX(x,   y, 1)) = r * b;
            (Hnn)(IDX(x,   y, 1), IDX(x+1, y, 0)) = r * conj(b);
        }
    }

    // kinetic energy in y direction
    // spin up

    for (int x = 0; x < Nx; x++){
        for (int y = 0; y < Ny - 1; y++) {
            cnum b = b_factor(yflux, x);
            cnum h = -V * b;
            Hnn(IDX(x, y, 0)  , IDX(x, y+1, 0)) = h;
            Hnn(IDX(x, y+1, 0), IDX(x, y  , 0)) = conj(h);
            Hnn(IDX(x, y, 1)  , IDX(x, y+1, 1)) = h;
            Hnn(IDX(x, y+1, 1), IDX(x, y  , 1)) = conj(h);
            // Rashba terms
            // "11 and 101"
            h = cnum(0, 1) * b * r;
            Hnn(IDX(x, y+1, 0), IDX(x, y  , 1)) = conj(h);
            Hnn(IDX(x, y  , 1), IDX(x, y+1, 0)) = h;
            // "1 and 111"
            Hnn(IDX(x, y  , 0), IDX(x, y+1, 1)) = h;
            Hnn(IDX(x, y+1, 1), IDX(x, y  , 0)) = conj(h);
        }
    }

    log_tick("hamiltonian");

    return H;

};

esm** self_energy(const num flux, const num gauge) {
    // analytical green's function in the leads
    // gl = G_{l+1, l+1}n
    Eigen::MatrixXcd gl(Nx, Nx);
    gl.setZero();


    for (int r = 0; r < Nx; r++) {
        num x = (e_tot - mods(r, lead_sites)) / (2.0 * V) + 1.0;
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
            real((unit - tmpp) / (unit-tmpm) * cnum(0.5, 0.0)));

        for (int p = 0; p < Nx; p++) {
            // "psiN1(ii)" in nano0903c.f
            cnum y1 = y * sin(pi * (num) ((p+1) * (r+1))/(1.0 + Nx));
            for (int q = 0; q < Nx; q++) {
                // "psiN1(jj)" in nano0903c.f
                cnum y2 = y * sin(pi * (num) ((q+1) * (r+1))/(1.0 + Nx));
                gl(p, q) += exp(cnum(0.0,1.0) * theta)/V * y1 * y2;
            }
        }
    }

    esm** e = new esm*[N_leads];
    ers** s = new ers*[N_leads];

    for (int i = 0; i < N_leads; i++) {
        assert(lead_sites + lead_offset[i] <= Nx);
        e[i]  = new esm(size, size);
        s[i]  = new ers( *e[i] );
    }

    for (int i = 0; i < lead_sites; i++){
        for (int j = 0; j < lead_sites; j++){
            cnum g = gl(i, j);

            /* left */
            (*s[0])(IDX(0, i+lead_offset[0], 0),
                    IDX(0, j+lead_offset[0], 0))      = g;
            (*s[2])(IDX(0, i+lead_offset[2], 1),
                    IDX(0, j+lead_offset[2], 1))      = g;

            /* right */
            (*s[1])(IDX(Nx-1, i+lead_offset[1], 0),
                    IDX(Nx-1, j+lead_offset[1], 0))   = g;
            (*s[3])(IDX(Nx-1, i+lead_offset[3], 1),
                    IDX(Nx-1, j+lead_offset[3], 1))   = g;

            /* top */
            (*s[4])(IDX(i+lead_offset[4], 0, 0),
                    IDX(j+lead_offset[4], 0, 0))      = g;
            /*  5 (sic) */
            (*s[5])(IDX(i+lead_offset[5], 0, 1),
                    IDX(j+lead_offset[5], 0, 1))      = g;

            /* bottom */
            /*   6 (sic) */
            (*s[6])(IDX(i+lead_offset[6], Ny-1, 0),
                    IDX(j+lead_offset[6], Ny-1, 0))   = g;
            (*s[7])(IDX(i+lead_offset[7], Ny-1, 1),
                    IDX(j+lead_offset[7], Ny-1, 1))   = g;
        }
    }

    for (int i = 0; i < N_leads; i++) {
        delete s[i];

        switch(i) {
            case 0:
            case 1:
            case 2:
            case 3:
                correct_phase(*e[i], -(1.0 - gauge) * flux);
                break;
            case 4:
            case 5:
            case 6:
            case 7:
//                cout << "gauge scaling factor: " << flux << endl;
                correct_phase(*e[i], gauge * flux);
                break;
        }
    }
    delete[] s;

    log_tick("self-energy");
    return e;
}

ub::matrix<num>* transmission(esm *H, const num flux, const num gauge) {
    esm **sigma_r   = self_energy(flux, gauge);

    for (int k = 0; k < N_leads; k++){
        *H -= *sigma_r[k];
    }
    esm e_green_inv(size, size);

    log_tick("hamiltonian + self-energy");
    // the second parameter is the ordering method that the solver uses
    // internally. Doesn't change results, only execution time
    eslu *slu = new eslu(*H, Eigen::MinimumDegree_ATA);
    delete H;
    H = NULL;
    log_tick("LU decomposition");


    ub::matrix<num> *tpq = new ub::matrix<num>(N_leads, N_leads);
    tpq->clear();
    esm **gamma_g_adv = new esm*[N_leads];
    esm **gamma_g_ret = new esm*[N_leads];

    // T_{p, q} = Trace( \Gamma_p * G^R * \Gamma_q * G^A )
    // where G^A = (G^R)^\dagger
    //
    //
    // \Gamma = -2 * Im(\Sigma_r)
    // (calculate slices of \Gamma (aka gamm_i) on the fly
    // to save memory
    // first carry out the two products
    // \Gamma_p * G^R and \Gamma_q * G^A
    for (int i = 0; i < N_leads; i++) {
        cout << "working on lead " << i << endl;
        esm *g_adv = new esm(size, size);
        esm *g_ret = new esm(size, size);


        assert(V * V == 1);
        *sigma_r[i] = (-2) * sigma_r[i]->imag();

        esm m1(size, size);
        esm result(size, size);

        // since *sigma_r[i] is a real matrix by now (although declared
        // complex) we can use transpose() instead of adjoint();
        pseudo_sparse_solve(slu, sigma_r[i]->transpose(), result, true);
        *g_ret = result.adjoint();

        pseudo_sparse_solve(slu, *sigma_r[i], result, false);
        delete sigma_r[i];
        sigma_r[i] = NULL;


        *g_adv = result.adjoint();

        gamma_g_adv[i] = g_adv;
        gamma_g_ret[i] = g_ret;
    }

    log_tick("solving");

    delete slu;
    slu = NULL;
    delete[] sigma_r;
    sigma_r = NULL;

    cout << "trace...\n";
    // Now calculate the trace
    cnum null = cnum(0, 0);
    for (int i = 0; i < N_leads; i++){
        for (int j = 0; j < N_leads; j++){
            (*tpq)(i, j) = r_prod_trace(*gamma_g_ret[i], *gamma_g_adv[j]);
        }
    }
    log_tick("trace");


    for (int i = 0; i < N_leads; i++){
        for (int n = 0; n < size; n++){
            cnum x = gamma_g_ret[i]->coeff(n, n);
            cnum y = gamma_g_adv[i]->coeff(n, n);
            (*tpq)(i, i) += real(cnum(0, 1) * y - cnum(0, 1) * x);

        }
        delete gamma_g_adv[i];
        delete gamma_g_ret[i];
    }
    for (int i = 0; i < lead_sites; i++){
        cnum k = findk(mods(i, lead_sites));
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
    num Bz = +6;

    int opt;
    while ((opt = getopt(argc, argv, "r:b:s:")) != -1) {
       switch (opt) {
            case 'r':
                alpha = atof(optarg);
                break;
            case 'b':
                Bz = atof(optarg);
                break;
            default:
                cerr << "Error while processing command line args\n";
                exit(1);

        }
    }

    for (int i = 0; i < N_leads; i++) {
        lead_offset[i] = 0;
    }

    cout << "PID:        " << getpid() << endl;
    cout << "Size:       " << Nx << "x" << Ny << endl;
    cout << "lead width: " << lead_sites << endl;
    cout << "Bz:         " << Bz << endl;
    esm *H = hamiltonian(alpha, Bz);
#ifndef NDEBUG
    {
        esm Hcheck(size, size);
        Hcheck = *H - esm(H->adjoint()).eval();
        num x = Hcheck.cwise().abs().sum();
        if (x > 0.01) {
            cerr << "ERROR: Hamiltonian is not hermitian ("
                 << x << ")\n";
            exit(1);
        }
    }
#endif

    num flux = flux_from_field(Bz);
    ub::matrix<num> *tpq = transmission(H, flux, global_gauge);
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
            cout << "Number of modes: " << ref << endl;
            is_first = false;
        }
        if (abs(r_sum - ref) > epsilon) {
            cout << "ERROR: sum rule violated for row " << i
                 << "  (" << r_sum << ")" << endl;
        }
        if (abs(c_sum - ref) > epsilon) {
            cout << "ERROR: sum rule violated for column " << i
                 << "  (" << c_sum << ")" << endl;
        }
    }
    log_tick("Done");
    delete tpq;
    return 0;
}

// vim: ft=cpp sw=4 ts=4 expandtab
