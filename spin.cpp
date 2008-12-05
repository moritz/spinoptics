#include <iostream>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <stdlib.h>
#include "invert-matrix.hpp"

using namespace boost::numeric::ublas;
using namespace std;

typedef long double num;
typedef complex<num> cnum;
typedef matrix<cnum> cmatrix;

const int Nx         = 10;
const int Ny         = Nx;
const int N_leads    = 8;

const num pi         = 3.14159265358979323846264;
const num h_bar      = 6.582122E-16;     // [eV*s]
const num h_planck   = 4.135669E-15;     // [eV*s]
const num electron_mass
                     = 9.109389E-31;     // [kg]
const num mass       = 0.381 * electron_mass; 
const num e_charge   = 1.60217653E-19;   // [C = A*s]

const num width_lead = 30.0;             // [nm]

// XXX is the +1 correct?
const num a_lead     = width_lead / (double) (Nx + 1); 
const int size       = Nx * Ny * 2;      // `Nfin'
const num V          = 1.0;              // hopping term

const num e_tot      = -2.0 * V;
const num width_disorder  = 0.0;

num alpha = -0.02;

template <class T>
void set_zero(matrix<T>* m) {
    for (unsigned int x = 0; x < m->size1(); x++){
        for (unsigned int y = 0; y < m->size2(); y++){
            (*m)(x, y) = 0.0;
        }
    }
}

inline num rashba(num alpha) {
    return 2.0 * alpha * a_lead;
}

inline num mods(int n, int nle) {
    return 2.0 * V * (cos(pi * (double) n / ((double) nle + 1)) -1);
}

num findk(num Emod) {
    num Etot = -2.0 * V; // ???
    return sqrt(
            2 * mass / (h_bar * h_bar) 
            * (Etot - Emod) * 10.0 / e_charge
        );
}

cmatrix* hamiltonian(num rashb) {
    cmatrix* Hnn = new cmatrix(size, size);
    cmatrix H = (*Hnn);
    set_zero(Hnn);

    // diagonal elements
    for (int i = 0; i < size / 2; i++) {
        // later we might want to add random disorder,
        // in which case these items might be different per
        // iteration, but every two diagonal items with distance 
        // (size/2) must still have the same value
        cnum energy = -4.0 * V;
        H(i, i)                     = energy;
        H(i + size/2, i + size/2)   = energy;
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
     *      0   1               (Nx-1)
     *      Nx  Nx+1            (2*Nx -1)    
     */

    // kinetic energy in x direction
    // No distinction between spin up and spin down needed
    for (int i = 0; i < size/2; i++) {
        if ((i+1) % Nx != 0) {
            H(i, i+1) = V;
            H(i+1, i) = V;
        }
    }

    // kinetic energy in y direction
    // spin up
    for (int i = 0; i < size / 2 - Nx; i++) {
        H(i, i + Nx) = V;
        H(i + Nx, i) = V;
    }
    // spin down
    for (int i = size / 2; i < size - Nx; i++) {
        H(i, i + Nx) = V;
        H(i + Nx, i) = V;
    }

    // Rashba terms
    // in x-direction
    for (int i = 0; i < size/2; i++) {
        if ((i+1) % Nx != 0) {
            // "1 and 102"
            // with spin flip
            H(i, i + Nx * Ny + 1) = rashba(alpha);
            H(i + Nx * Ny + 1, i) = rashba(alpha);
            // "101 and 2"
            H(i + 1, i + Nx * Ny) = - rashba(alpha);
            H(i + Nx * Ny, i + 1) = - rashba(alpha);
        }
    }

    // in y-direction
    for (int i = 0; i < Nx * (Ny -1); i++) {
        // "11 and 101"
        // with spin flip
        H(i + Nx, i) = cnum(0, 1) * rashba(alpha);
        // "1 and 111"
        H(i, i + Nx) = cnum(0, -1) * rashba(alpha);
    }

    std::cout << H;
    return Hnn;
};

cmatrix** self_energy(void) {
    // Glp1lp1n = G_{l+1, l+1}n
    cmatrix Glp1lp1n = cmatrix(N_leads, N_leads);
    set_zero(&Glp1lp1n);
    for (int p = 0; p < N_leads; p++) {
        for (int q = 0; q < N_leads; q++) {
            for (int r = 0; r < N_leads; r++) {
                num x = (e_tot - mods(r+1, N_leads)) / (2.0 * V) + 1.0;
                cnum theta;
                if (x > 1.0) {
                    theta = cnum(0, 1)
                        * log((cnum) (x + sqrt(x*x - 1.0)));
                } else if (x < -1.0) {
                    theta = cnum(0, 1)
                        * log((cnum) (x - sqrt(x*x - 1.0)));
                } else {
                    theta = acos(x);
                }
                cnum unit = cnum(1.0, 0.0);
                cnum tmpp = exp(cnum(0, 2.0) * pi * 
                        ((cnum) N_leads) / ((cnum) (N_leads + 1)));
                cnum tmpm = exp(cnum(0, -2.0) * pi * 
                        ((cnum) N_leads) / ((cnum) (N_leads + 1)));
                // "AnorN1(mm)" in nano0903c.f
                cnum y = sqrt(cnum(0.5 * N_leads, 0.0) + 
                    (unit - tmpp) / (unit-tmpp) * cnum(0.5, 0.0));
                // "psiN1(ii)" in nano0903c.f
                cnum y1 = y * sin(pi * (num) (p * r)/(1.0 + N_leads));
                // "psiN1(jj)" in nano0903c.f
                cnum y2 = y * sin(pi * (num) (q * r)/(1.0 + N_leads));
                Glp1lp1n(p, q) += exp(cnum(0.0,1.0) * theta)/V * y1 * y2;
            }
        }
    }

    cmatrix *G_lp1_lp1_up   = new cmatrix(size, size); set_zero(G_lp1_lp1_up);
    cmatrix *G_lp1_lp1_down = new cmatrix(size, size); set_zero(G_lp1_lp1_down);
    cmatrix *G_lm1_lm1_up   = new cmatrix(size, size); set_zero(G_lm1_lm1_up);
    cmatrix *G_lm1_lm1_down = new cmatrix(size, size); set_zero(G_lm1_lm1_down);

    cmatrix *G_xp1_xp1_up   = new cmatrix(size, size); set_zero(G_xp1_xp1_up);
    cmatrix *G_xp1_xp1_down = new cmatrix(size, size); set_zero(G_xp1_xp1_down);
    cmatrix *G_xm1_xm1_up   = new cmatrix(size, size); set_zero(G_xm1_xm1_up);
    cmatrix *G_xm1_xm1_down = new cmatrix(size, size); set_zero(G_xm1_xm1_down);

    int s = size / 2;
    for (int i = 1; i <= Nx; i++){
        int n = Nx * i;
        for (int j = 1; j <= Ny; j++){
            int m = Ny * j;
            cnum g = V * V * Glp1lp1n(i, j);
            (*G_lp1_lp1_up)(m, n)                   = g;
            (*G_lm1_lm1_up)(m + Nx - 1, n + Nx - 1) = g;
            (*G_lp1_lp1_down)(m + s, n + s)         = g;
            (*G_lm1_lm1_down)(m+s+Nx-1, n+s+Nx- 1)  = g;

            (*G_xp1_xp1_up)(i, j)                   = g;
            (*G_xp1_xp1_down)(i + s, j + s)         = g;
            (*G_xm1_xm1_up)(i + s - Nx, i + s - Nx) = g;
            (*G_xm1_xm1_down)(i+size-Nx, i+size-Nx) = g;
        }
    }
    cmatrix** sr = new cmatrix*[8];
//    std::vector<cmatrix*> *sr = new std::vector<matrix*>();
    sr[0] = G_lp1_lp1_up;
    sr[2] = G_lp1_lp1_down;
    sr[4] = G_xp1_xp1_up;
    sr[6] = G_xp1_xp1_down;

    sr[1] = G_lm1_lm1_up;
    sr[3] = G_lm1_lm1_down;
    sr[5] = G_xm1_xm1_up;
    sr[7] = G_xm1_xm1_down;
    return sr;
}

cmatrix* greenji(num rashba) {
    cmatrix *Hnn        = hamiltonian(rashba);
    cmatrix **sigma_r   = self_energy();
    cmatrix *green_inv  = new cmatrix(size, size);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            (*green_inv)(i, j) = (*Hnn)(i, j);
            for (int k = 0; k < N_leads; k++){
                (*green_inv)(i, j) -= (*(sigma_r[k]))(i, j);
            }
        }
    }
    cmatrix *green = new cmatrix(size, size);
    InvertMatrix(*green_inv, *green);
    delete green_inv;

    cmatrix *tpq = new cmatrix(N_leads, N_leads);
    set_zero(tpq);
    cmatrix **gamma_g_adv = new cmatrix*[N_leads];
    cmatrix **gamma_g_ret = new cmatrix*[N_leads];
    cnum two = cnum(2, 0);

    // T_{p, q} = Trace( \Gamma_p G^R \Gamma_q G^A )
    // where G^A = (G^R)^* 
    //
    // first carry out the first two products
    for (int i = 0; i < N_leads; i++) {
        cmatrix *g_adv = new cmatrix(size, size);
        cmatrix *g_ret = new cmatrix(size, size);
        set_zero(g_adv);
        set_zero(g_ret);
        for (int n = 0; n < N_leads; i++){
            for (int m = 0; n < N_leads; i++){
                for (int nn = 0; nn < N_leads; nn++){
                    (*g_adv)(n, m) -= two * (*green)(nn, m) 
                                      * imag((*(sigma_r[i]))(m, n));
                    (*g_ret)(n, m) -= two * conj((*green)(nn, m)) 
                                      * imag((*(sigma_r[i]))(m, n));
                }
            }
        }
        gamma_g_adv[i] = g_adv;
        gamma_g_ret[i] = g_ret;
    }
    // Now calculate the trace
    for (int i = 0; i < N_leads; i++){
        for (int j = 0; j < N_leads; j++){
            for (int n = 0; n < size; n++){
                for (int m = 0; m < size; m++){
                    (*tpq)(i, j) += (*gamma_g_ret[i])(n, m)
                                    * (*gamma_g_adv[j])(m, n);
                }
            }
        }
    }

    for (int i = 0; i < N_leads; i++){
        for (int n = 0; n < size; n++){
            (*tpq)(i, i) += cnum(0, 1) * (*(gamma_g_adv[i]))(n, n)
                          - cnum(0, 1) * (*(gamma_g_ret[i]))(n, n);

        }
    }
    // clean up temporary variables
    for (int i = 0; i < N_leads; i++){
        delete gamma_g_adv[i];
        delete gamma_g_ret[i];
    }
    delete[] gamma_g_adv;
    delete[] gamma_g_ret;

    return tpq;
}


int main (int argc, char** argv) {
    return 0;
}

// vim: ft=cpp sw=4 ts=4 expandtab
