#include <iostream>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <stdlib.h>

using namespace boost::numeric::ublas;
using namespace std;

typedef long double num;
typedef complex<num> cnum;
typedef matrix<cnum> cmatrix;

const int    Nx         = 10;
const int    Ny         = Nx;
const int    N_leads    = Nx;

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

inline num rashba(num alpha) {
    return 2.0 * alpha * a_lead;
}

inline num mods(int n, int nle) {
    return 2.0 * V * (cos(pi * (double) n / ((double) N_leads + 1)) -1);
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
    for (int x = 0; x < size; x++) {
        for (int y = 0; y < size; y++) {
            H(x, y) = 0.0;
        }
    }

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
            H(i, i + Nx * Ny + 1) = rashba_v;
            H(i + Nx * Ny + 1, i) = rashba_v;
            // "101 and 2"
            H(i + 1, i + Nx * Ny) = - rashba_v;
            H(i + Nx * Ny, i + 1) = - rashba_v;
        }
    }

    // in y-direction
    for (int i = 0; i < Nx * (Ny -1); i++) {
        // "11 and 101"
        H(i + Nx, i) = complex(0, 1) * rashba_v;
        // "1 and 111"
        H(i, i + Nx) = complex(0, -1) * rashba_v;
    }

    std::cout << H;
    return Hnn;
}

matrix<num>* greenji(num rashba) {
    cmatrix *Hnn = hamiltonian(rashba);

}


int main (int argc, char** argv) {
/*    for (num alpha = 0.02; fabs(alpha) <= fabs(V); alpha *= sqrt(10.0)) {
        num r = rashba(alpha);

    } */
    return 0;
}

// vim: ft=cpp sw=4 ts=4 expandtab
