#include "polynomial.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    int degree;
    double tolerance;
    printf("Enter degree of the polynomial (2-5): ");
    if (scanf("%d", &degree) != 1 || degree < 2 || degree > 5) {
        fprintf(stderr, "Degree must be between 2 and 5.\n");
        return EXIT_FAILURE;
    }
    double *coefficients = malloc((size_t)(degree + 1) * sizeof(*coefficients));
    Root *roots = malloc((size_t)degree * sizeof(*roots));
    if (!coefficients || !roots) {
        fprintf(stderr, "Unable to allocate memory.\n");
        free(coefficients); free(roots);
        return EXIT_FAILURE;
    }
    for (int i = 0; i <= degree; ++i) {
        printf("Enter coefficient of x^%d: ", degree - i);
        if (scanf("%lf", &coefficients[i]) != 1 || !isfinite(coefficients[i])) {
            fprintf(stderr, "Invalid coefficient.\n");
            free(coefficients); free(roots);
            return EXIT_FAILURE;
        }
    }
    if (coefficients[0] == 0.0) {
        fprintf(stderr, "The leading coefficient cannot be zero.\n");
        free(coefficients); free(roots);
        return EXIT_FAILURE;
    }
    printf("Enter positive tolerance: ");
    if (scanf("%lf", &tolerance) != 1 || !isfinite(tolerance) || tolerance <= 0.0) {
        fprintf(stderr, "Tolerance must be a positive number.\n");
        free(coefficients); free(roots);
        return EXIT_FAILURE;
    }
    int count = solve_polynomial(coefficients, degree, tolerance, roots);
    if (!count) {
        fprintf(stderr, "Bairstow iteration did not converge. Try a different tolerance.\n");
        free(coefficients); free(roots);
        return EXIT_FAILURE;
    }
    print_roots(roots, count);
    int verified =
        verify_reconstruction(coefficients, degree, roots, count, tolerance);
    free(coefficients); free(roots);
    return verified ? EXIT_SUCCESS : EXIT_FAILURE;
}
