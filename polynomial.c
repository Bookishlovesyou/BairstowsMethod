#include "polynomial.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_ITERATIONS 2000

static void remainders(const double *a, int degree, double r, double s,
                       double *b, double *rem_r, double *rem_s) {
    b[0] = a[0];
    b[1] = a[1] - r * b[0];
    for (int i = 2; i <= degree; ++i)
        b[i] = a[i] - r * b[i - 1] - s * b[i - 2];
    *rem_r = b[degree - 1];
    *rem_s = b[degree];
}

static int factor_quadratic(const double *a, int degree, double tolerance,
                            double *r, double *s, double *quotient) {
    double *b = malloc((size_t)(degree + 1) * sizeof(*b));
    double scale = 1.0;
    if (!b) return 0;
    for (int i = 0; i <= degree; ++i)
        if (fabs(a[i]) > scale) scale = fabs(a[i]);

    for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        double f, g;
        remainders(a, degree, *r, *s, b, &f, &g);
        printf("%4d  r=% .10g  s=% .10g  remainder=(% .3e, % .3e)\n",
               iteration + 1, *r, *s, f, g);
        if (fabs(f) <= tolerance * scale && fabs(g) <= tolerance * scale) {
            for (int i = 0; i <= degree - 2; ++i) quotient[i] = b[i];
            free(b);
            return 1;
        }

        
        double h = 1e-6 * (1.0 + fabs(*r) + fabs(*s));
        double fp, gp, fq, gq;
        remainders(a, degree, *r + h, *s, b, &fp, &gp);
        remainders(a, degree, *r, *s + h, b, &fq, &gq);
        double j11 = (fp - f) / h, j21 = (gp - g) / h;
        double j12 = (fq - f) / h, j22 = (gq - g) / h;
        double determinant = j11 * j22 - j12 * j21;
        if (!isfinite(determinant) || fabs(determinant) < 1e-18 * scale) {
            free(b);
            return 0;
        }
        double dr = (-f * j22 + g * j12) / determinant;
        double ds = (-g * j11 + f * j21) / determinant;
        if (!isfinite(dr) || !isfinite(ds)) { free(b); return 0; }

       
        double factor = 1.0, old_norm = hypot(f, g);
        int accepted = 0;
        for (int trial = 0; trial < 24; ++trial, factor *= 0.5) {
            double nf, ng;
            remainders(a, degree, *r + factor * dr, *s + factor * ds,
                       b, &nf, &ng);
            if (isfinite(nf) && isfinite(ng) && hypot(nf, ng) < old_norm) {
                *r += factor * dr;
                *s += factor * ds;
                accepted = 1;
                break;
            }
        }
        if (!accepted) { free(b); return 0; }
    }
    free(b);
    return 0;
}

static void quadratic_roots(double r, double s, Root *roots) {
    double discriminant = r * r - 4.0 * s;
    if (discriminant >= 0.0) {
        double d = sqrt(discriminant);
        roots[0] = (Root){(-r + d) / 2.0, 0.0};
        roots[1] = (Root){(-r - d) / 2.0, 0.0};
    } else {
        roots[0] = (Root){-r / 2.0, sqrt(-discriminant) / 2.0};
        roots[1] = (Root){-r / 2.0, -sqrt(-discriminant) / 2.0};
    }
}

int solve_polynomial(const double *coefficients, int degree, double tolerance,
                     Root *roots) {
    double *work = malloc((size_t)(degree + 1) * sizeof(*work));
    if (!work) return 0;
    for (int i = 0; i <= degree; ++i) work[i] = coefficients[i] / coefficients[0];
    int remaining = degree, root_count = 0;

    while (remaining > 2) {
        double r = 0.0, s = 0.0;
        double *quotient = malloc((size_t)(remaining - 1) * sizeof(*quotient));
        if (!quotient) { free(work); return 0; }
        int found = 0;
         
        for (int attempt = 0; attempt < 12 && !found; ++attempt) {
            if (attempt) {
                r = (attempt % 2 ? 1.0 : -1.0) * 0.5 * attempt;
                s = (attempt % 3 ? 1.0 : -1.0) * 0.25 * attempt;
            }
            printf("\nBairstow iteration table (degree %d, start %d)\n",
                   remaining, attempt + 1);
            found = factor_quadratic(work, remaining, tolerance, &r, &s,
                                     quotient);
        }
        if (!found) { free(quotient); free(work); return 0; }
        quadratic_roots(r, s, roots + root_count);
        root_count += 2;
        for (int i = 0; i <= remaining - 2; ++i) work[i] = quotient[i];
        free(quotient);
        remaining -= 2;
    }

    if (remaining == 2) {
        double disc = work[1] * work[1] - 4.0 * work[2];
        if (disc >= 0.0) {
            roots[root_count++] = (Root){(-work[1] + sqrt(disc)) / 2.0, 0.0};
            roots[root_count++] = (Root){(-work[1] - sqrt(disc)) / 2.0, 0.0};
        } else {
            roots[root_count++] = (Root){-work[1] / 2.0, sqrt(-disc) / 2.0};
            roots[root_count++] = (Root){-work[1] / 2.0, -sqrt(-disc) / 2.0};
        }
    } else if (remaining == 1) {
        roots[root_count++] = (Root){-work[1] / work[0], 0.0};
    }
    free(work);
    return root_count;
}

void print_roots(const Root *roots, int count) {
    puts("\nRoots:");
    for (int i = 0; i < count; ++i) {
        if (fabs(roots[i].imag) < 1e-12)
            printf("x%d = %.10g\n", i + 1, roots[i].real);
        else
            printf("x%d = %.10g %c %.10gi\n", i + 1, roots[i].real,
                   roots[i].imag < 0 ? '-' : '+', fabs(roots[i].imag));
    }
}
