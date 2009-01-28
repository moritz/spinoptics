#include <iostream>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <stdlib.h>
#include <time.h>
#include "invert-matrix.hpp"

#define IDX(x, y) ((x) + Nx * (y))

using namespace boost::numeric::ublas;
using namespace std;

typedef double num;
typedef complex<num> cnum;
typedef matrix<cnum> cmatrix;
typedef compressed_matrix<cnum> sparse_cmatrix;

const int Nx         = 4;
const int Ny         = Nx;
const int N_leads    = 8;

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
const num g_factor   = 2.0;

// XXX is the +1 correct?
const num a_lead     = width_lead / (double) (Nx + 1);
const int size       = Nx * Ny * 2;      // `Nfin'
const num V          = -1.0;             // hopping term

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

cnum b_field(const num B, const int y) {
    return exp(cnum(0.0, e_charge / h_bar * B * a_lead * ((num) y)) );
}

template <class T>
void set_zero(matrix<T>* m) {
    for (unsigned int x = 0; x < m->size1(); x++){
        for (unsigned int y = 0; y < m->size2(); y++){
            (*m)(x, y) = 0.0;
        }
    }
}

inline num rashba(const num alpha) {
    return 2.0 * alpha * a_lead;
}

inline num mods(const int n, const int nle) {
    return 2.0 * V * (cos(pi * (double) (n+1) / ((double) nle + 1.0)) 
            - 1.0);
}

cnum findk(const num Emod) {
    return sqrt(
            cnum(2 * mass / (h_bar * h_bar)
            * (e_tot - Emod) * 10.0 / 1.60219, 0)
        );
}

sparse_cmatrix* hamiltonian(const num rashb, const num B) {
    sparse_cmatrix* Hnn = new sparse_cmatrix(size, size, 4 * size);
//    set_zero(Hnn);

    // division by e_charge to convert from electron volt to Joule
    num wb = 0.5 * B * g_factor * bohr_magneton / e_charge;

    // diagonal elements
    for (int i = 0; i < size / 2; i++) {
        // later we might want to add random disorder,
        // in which case these items might be different per
        // iteration, but every two diagonal items with distance
        // (size/2) must still have the same value
        cnum energy = 4.0 * V + e_tot;
        (*Hnn)(i, i)                     = energy + wb;
        (*Hnn)(i + size/2, i + size/2)   = energy - wb;
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

    // kinetic energy in x direction
    // No distinction between spin up and spin down needed
    int s = size / 2;
    for (int i = 0; i < size; i++) {
        if ((i+1) % Nx != 0) {
            (*Hnn)(i, i+1) = -V * b_field(B, i / Nx);
            (*Hnn)(i+1, i) = -V * conj(b_field(B, i / Nx));
        }
    }

    // kinetic energy in y direction
    // spin up
    for (int i = 0; i < size / 2 - Nx; i++) {
        (*Hnn)(i, i + Nx) = -V;
        (*Hnn)(i + Nx, i) = -V;
    }
    // spin down
    for (int i = size / 2; i < size - Nx; i++) {
        (*Hnn)(i, i + Nx) = -V;
        (*Hnn)(i + Nx, i) = -V;
    }

    // Rashba terms
    // in x-direction
    for (int i = 0; i < size/2; i++) {
        if ((i+1) % Nx != 0) {
            // "1 and 102"
            // with spin flip
            (*Hnn)(i, i + s + 1) = rashba(alpha) * b_field(B, i / Nx);
            (*Hnn)(i + s + 1, i) = rashba(alpha) * conj( b_field(B, i / Nx));
            // "101 and 2"
            (*Hnn)(i + 1, i + s) = - rashba(alpha);
            (*Hnn)(i + s, i + 1) = - rashba(alpha);
        }
    }

    // in y-direction
    for (int i = 0; i < Nx * (Ny -1); i++) {
        // "11 and 101"
        // with spin flip
        (*Hnn)(i + Nx, i + s) = cnum(0,  1) * rashba(alpha);
        (*Hnn)(i + s, i + Nx) = cnum(0, -1) * rashba(alpha);
        // "1 and 111"
        (*Hnn)(i, i + s + Nx) = cnum(0, -1) * rashba(alpha);
        (*Hnn)(i + s + Nx, i) = cnum(0,  1) * rashba(alpha);
    }
    log_tick("hamiltonian");
    return Hnn;
};

sparse_cmatrix** self_energy(void) {
    // Glp1lp1n = G_{l+1, l+1}n
    cmatrix Glp1lp1n = cmatrix(Nx, Nx);
    set_zero(&Glp1lp1n);
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
        cnum tmpp = exp(cnum(0, 1) * (2.0 * pi * ((num) ((r+1) * Nx )
                        / (num) (Nx + 1))));
        cnum tmpm = exp(cnum(0, -1) * (2.0 * pi * ((num) (r+1)
                        / (num) (Nx + 1))));

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
//    cout << "Glp1lp1n: " << Glp1lp1n << endl;
    // Glp1lp1n identical with that of nano0903c.f

    sparse_cmatrix *G_lp1_lp1_up   = new sparse_cmatrix(size, size, size / 2);
    sparse_cmatrix *G_lp1_lp1_down = new sparse_cmatrix(size, size, size / 2);
    sparse_cmatrix *G_lm1_lm1_up   = new sparse_cmatrix(size, size, size / 2);
    sparse_cmatrix *G_lm1_lm1_down = new sparse_cmatrix(size, size, size / 2);

    sparse_cmatrix *G_xp1_xp1_up   = new sparse_cmatrix(size, size, size / 2);
    sparse_cmatrix *G_xp1_xp1_down = new sparse_cmatrix(size, size, size / 2);
    sparse_cmatrix *G_xm1_xm1_up   = new sparse_cmatrix(size, size, size / 2);
    sparse_cmatrix *G_xm1_xm1_down = new sparse_cmatrix(size, size, size / 2);

    int s = size / 2;
    for (int i = 0; i < Nx; i++){
        int n = Nx * i;
        for (int j = 0; j < Ny; j++){
            int m = Ny * j;
            cnum g = Glp1lp1n(i, j);
            (*G_lp1_lp1_up)(m, n)                   = g;
            (*G_lm1_lm1_up)(m + Nx - 1, n + Nx - 1) = g;
            (*G_lp1_lp1_down)(m + s, n + s)         = g;
            (*G_lm1_lm1_down)(m+s+Nx-1, n+s+Nx- 1)  = g;

            (*G_xp1_xp1_up)(i, j)                   = g;
            (*G_xp1_xp1_down)(i + s, j + s)         = g;
            (*G_xm1_xm1_up)(i + s - Nx, j + s - Nx) = g;
            (*G_xm1_xm1_down)(i + size - Nx, j+size-Nx) = g;
        }
    }
    sparse_cmatrix** sr = new sparse_cmatrix*[N_leads];
    sr[0] = G_lp1_lp1_up;
    sr[2] = G_lp1_lp1_down;
    sr[4] = G_xp1_xp1_up;
    // 5, not 6
    sr[5] = G_xp1_xp1_down;

    sr[1] = G_lm1_lm1_up;
    sr[3] = G_lm1_lm1_down;
    // 6, not 5
    sr[6] = G_xm1_xm1_up;
    sr[7] = G_xm1_xm1_down;
//    cout << "G_xp1_xp1_up: " << *G_xp1_xp1_up << endl;
    // re-checked: G_xp1_xp1_*
    log_tick("self-energy");
    return sr;
}

matrix<num>* greenji(sparse_cmatrix* Hnn) {
    sparse_cmatrix **sigma_r   = self_energy();

    cmatrix *green_inv  = new cmatrix(size, size);
    (*green_inv) = *Hnn;
    for (int k = 0; k < N_leads; k++){
        noalias(*green_inv) -= *(sigma_r[k]);
    }
    cmatrix *green = new cmatrix(size, size);
    log_tick("green_inv");
    InvertMatrix(*green_inv, *green);
//    cout << "Green: " << *green << endl;
    log_tick("matrix inversion");
    delete green_inv;
    green_inv = NULL;

    cmatrix *green_herm = new cmatrix(size, size);
    *green_herm = herm(*green);

    matrix<num> *tpq = new matrix<num>(N_leads, N_leads);
    set_zero(tpq);
    cmatrix **gamma_g_adv = new cmatrix*[N_leads];
    cmatrix **gamma_g_ret = new cmatrix*[N_leads];

    // T_{p, q} = Trace( \Gamma_p G^R \Gamma_q G^A )
    // where G^A = (G^R)^*
    //
    //
    // \Gamma = -2 * Im(\Sigma_r)   
    // (calculate slices of \Gamma (aka gamm_i) on the fly
    // to save memory
    // first carry out the two products \Gamma_p * G^R and \Gamma_q * G^A
    cmatrix * gamm_i = new cmatrix(size, size);
    for (int i = 0; i < N_leads; i++) {
        cmatrix *g_adv = new cmatrix(size, size);
        cmatrix *g_ret = new cmatrix(size, size);
        noalias(*gamm_i) = -2 * imag(*sigma_r[i]);
        // axpy_prod writes its result into the third argument
        axpy_prod(*gamm_i, *green, *g_ret, true);
        axpy_prod(*gamm_i, *green_herm, *g_adv, true);
        gamma_g_adv[i] = g_adv;
        gamma_g_ret[i] = g_ret;
    }
    log_tick("products");
    delete gamm_i;      gamm_i      = NULL;
    delete green;       green       = NULL;
    delete green_herm;  green_herm  = NULL;

    // we don't need sigma_r any more
    for (int i = 0; i < N_leads; i++){
        delete sigma_r[i];
    }
    delete[] sigma_r;
    sigma_r = NULL;

    // Now calculate the trace
    for (int i = 0; i < N_leads; i++){
        for (int j = 0; j < N_leads; j++){
            for (int n = 0; n < size; n++){
                for (int m = 0; m < size; m++){
                    (*tpq)(i, j) += real((*gamma_g_ret[i])(n, m)
                                    * (*gamma_g_adv[j])(m, n));
                }
            }
        }
    }
    log_tick("trace");

    for (int i = 0; i < N_leads; i++){
        for (int n = 0; n < size; n++){
            (*tpq)(i, i) += real(cnum(0, 1) * (*(gamma_g_adv[i]))(n, n)
                          - cnum(0, 1) * (*(gamma_g_ret[i]))(n, n));

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
    sparse_cmatrix *Hnn = hamiltonian(0.3, 0.0);
    cout << "Hamiltonian: " << *Hnn << "\n";
    matrix<num> *tpq = greenji(Hnn);
    delete Hnn;
    cout << "final tpq" << *tpq << endl;
    boost::numeric::ublas::vector<num> r;
    boost::numeric::ublas::vector<num> c;
    for (int i = 0; i < N_leads; i++) {
        num r_sum = 0.0;
        num c_sum = 0.0;
        r = row(*tpq, i);
        c = column(*tpq, i);
        for (int j = 0; j < N_leads; j++) {
            r_sum += r(j);
            c_sum += c(j);
        }
        cout << i <<  "\t" <<  r_sum << "\t" << c_sum << "\n";
    }
    delete tpq;
}

// vim: ft=cpp sw=4 ts=4 expandtab
