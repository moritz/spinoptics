#include <stdio.h>
#include <math.h>


#define     PLUS        1.0
#define     MINUS       -1.0

#define     UP          2
#define     DOWN        3

#define     pi          3.1415926535897932384626433832795

double alpha    = 0.1;

double chi_n_norm(double pm, double phi) {
    return sqrt(sin(phi) * sin(phi) + pow(pm - cos(phi), 2));
}

double p_so(double pm) {
    return sqrt(1 + alpha*alpha) - pm * alpha;
}

double p_x_so(double pm, double phi) {
    return sqrt(pow(p_so(pm),2) - pow(sin(phi), 2));
}

double chi_n(double pm, int component, double phi) {
    if (component == UP) {
        return (-cos(phi) + pm) / chi_n_norm(pm, phi);
    } else {
        return sin(phi) / chi_n_norm(pm, phi);
    }
}

double chi_so_norm(double pm, double phi) {
    return sqrt(
              pow(-p_x_so(pm, phi) + pm * p_so(pm), 2)
            + pow(sin(phi), 2)
           );
}

double chi_so(double pm, int component, double phi) {
    if (component == UP) {
        return (-p_x_so(pm, phi) + pm * p_so(pm)) / chi_so_norm(pm, phi);
    } else {
        return sin(phi) / chi_so_norm(pm, phi);
    }
}

inline double A(double phi) {
    return p_x_so(PLUS, phi) / cos(phi) + 1.0;
}

inline double C(double phi) {
    return p_x_so(MINUS, phi) / cos(phi) + 1.0;
}

inline double B(double phi) {
    return alpha / cos(phi);
}

int main(int argc, char** argv) {
    for (double phi = 0.01; phi < pi / 2.0; phi += 0.01) {
        double nom = 2.0 * ( 
                chi_n(PLUS, UP,   phi) * (C(phi) + B(phi)) * chi_so(MINUS, DOWN, phi)
              - chi_n(PLUS, DOWN, phi) * (C(phi) - B(phi)) * chi_so(MINUS, UP,   phi) 
            );
        double denom = 
            (A(phi) - B(phi)) * (C(phi) + B(phi)) * chi_so(MINUS, DOWN, phi)
                                                  * chi_so(PLUS,  UP,   phi)
           -(A(phi) + B(phi)) * (C(phi) - B(phi)) * chi_so(MINUS, UP,   phi)
                                                  * chi_so(PLUS,  DOWN, phi);
        if (isnan(nom/denom))
            break;
        printf("%.5f\t%.5f\t%.5f\t%.5f\n", phi, nom/denom, 
                chi_n(PLUS, UP,   phi) * (C(phi) + B(phi)) * chi_so(MINUS, DOWN, phi),
              - chi_n(PLUS, DOWN, phi) * (C(phi) - B(phi)) * chi_so(MINUS, UP,   phi) 
              );

    }
}

/* vim: ts=4 sw=4 expandtab
 */
