#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

typedef struct {
    double real;
    double imag;
} Root;

int solve_polynomial(const double *coefficients, int degree, double tolerance,
                     Root *roots);
void print_roots(const Root *roots, int count);

#endif
