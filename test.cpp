#include <iostream>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stdlib.h>
#include "math-utils.h"
using namespace std;


int main(int argc, char** argv) {
    esm a(3, 3);
    esm b(3, 3);
    {
        ers s(a);
        s(0, 0) = cnum(1.0, 0.0);
        s(0, 2) = cnum(2.0, 0.0);

        s(1, 1) = cnum(4.0, 0.0);
        s(1, 2) = cnum(3.0, 0.0);

        s(2, 0) = cnum(2.0, 0.0);
        s(2, 1) = cnum(2.0, 0.0);
    }

    {
        ers s(b);
        s(0, 0) = cnum(1.0, 0.0);
        s(0, 1) = cnum(1.0, 0.0);
        s(0, 2) = cnum(1.0, 0.0);

        s(1, 0) = cnum(2.0, 0.0);
        s(1, 1) = cnum(2.0, 0.0);
        s(1, 2) = cnum(2.0, 0.0);

        s(2, 0) = cnum(3.0, 0.0);
        s(2, 1) = cnum(3.0, 0.0);
        s(2, 2) = cnum(3.0, 0.0);

    }

    ZMUMPS_STRUC_C * mumps = MUMPS_lr(a);
    esm c(3, 3);
    MUMPS_solve(mumps, b, c, 5);
    cout << c;
    MUMPS_free(mumps);

}
 // vim: ts=4 sw=4 expandtab
