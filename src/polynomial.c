#include "polynomial.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_ITERATIONS 2000
#define MAX_STARTS 12

typedef struct {
    int iteration;
    double r;
    double s;
    double remainder_r;
    double remainder_s;
} IterationRow;

typedef struct {
    int degree;
    int start;
    IterationRow *rows;
    size_t row_count;
    size_t row_capacity;
} IterationTable;

typedef struct {
    IterationTable *tables;
    size_t count;
    size_t capacity;
} TableCollection;

static void free_tables(TableCollection *collection) {
    for (size_t i = 0; i < collection->count; ++i)
        free(collection->tables[i].rows);
    free(collection->tables);
    collection->tables = NULL;
    collection->count = 0;
    collection->capacity = 0;
}

static IterationTable *new_table(TableCollection *collection, int degree,
                                 int start) {
    if (collection->count == collection->capacity) {
        size_t capacity = collection->capacity ? collection->capacity * 2 : 4;
        IterationTable *tables =
            realloc(collection->tables, capacity * sizeof(*tables));
        if (!tables) return NULL;
        collection->tables = tables;
        collection->capacity = capacity;
    }
    IterationTable *table = &collection->tables[collection->count++];
    *table = (IterationTable){0};
    table->degree = degree;
    table->start = start;
    return table;
}

static int append_row(IterationTable *table, IterationRow row) {
    if (table->row_count == table->row_capacity) {
        size_t capacity = table->row_capacity ? table->row_capacity * 2 : 16;
        IterationRow *rows = realloc(table->rows, capacity * sizeof(*rows));
        if (!rows) return 0;
        table->rows = rows;
        table->row_capacity = capacity;
    }
    table->rows[table->row_count++] = row;
    return 1;
}

static void print_table(const IterationTable *table, size_t number,
                        const char *position) {
    printf("\n%s Bairstow iteration table %zu (degree %d, start %d)\n",
           position, number, table->degree, table->start);
    for (size_t i = 0; i < table->row_count; ++i) {
        const IterationRow *row = &table->rows[i];
        printf("%4d  r=% .10g  s=% .10g  remainder=(% .3e, % .3e)\n",
               row->iteration, row->r, row->s, row->remainder_r,
               row->remainder_s);
    }
}

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
                            double *r, double *s, double *quotient,
                            IterationTable *table) {
    double *b = malloc((size_t)(degree + 1) * sizeof(*b));
    double scale = 1.0;
    if (!b) {
        fprintf(stderr, "Unable to allocate Bairstow workspace.\n");
        return 0;
    }
    for (int i = 0; i <= degree; ++i)
        if (fabs(a[i]) > scale) scale = fabs(a[i]);

    for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        double f, g;
        remainders(a, degree, *r, *s, b, &f, &g);
        if (!append_row(table, (IterationRow){iteration + 1, *r, *s, f, g})) {
            fprintf(stderr, "Unable to store Bairstow iteration table.\n");
            free(b);
            return -1;
        }
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
    if (!work) {
        fprintf(stderr, "Unable to allocate polynomial workspace.\n");
        return 0;
    }
    TableCollection tables = {0};
    for (int i = 0; i <= degree; ++i) work[i] = coefficients[i] / coefficients[0];
    int remaining = degree, root_count = 0;

    while (remaining > 2) {
        double r = 0.0, s = 0.0;
        double *quotient = malloc((size_t)(remaining - 1) * sizeof(*quotient));
        if (!quotient) {
            fprintf(stderr, "Unable to allocate quotient workspace.\n");
            free(work);
            free_tables(&tables);
            return 0;
        }
        int found = 0;
         
        for (int attempt = 0; attempt < MAX_STARTS && !found; ++attempt) {
            if (attempt) {
                r = (attempt % 2 ? 1.0 : -1.0) * 0.5 * attempt;
                s = (attempt % 3 ? 1.0 : -1.0) * 0.25 * attempt;
            }
            IterationTable *table = new_table(&tables, remaining, attempt + 1);
            if (!table) {
                fprintf(stderr, "Unable to allocate Bairstow table.\n");
                free(quotient);
                free(work);
                free_tables(&tables);
                return 0;
            }
            int factor_result = factor_quadratic(work, remaining, tolerance,
                                                  &r, &s, quotient, table);
            if (factor_result < 0) {
                free(quotient);
                free(work);
                free_tables(&tables);
                return 0;
            }
            found = factor_result;
        }
        if (!found) {
            free(quotient);
            free(work);
            free_tables(&tables);
            return 0;
        }
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
    printf("\nNumber of Bairstow tables: %zu\n", tables.count);
    if (tables.count) {
        size_t middle = tables.count / 2;
        print_table(&tables.tables[0], 1,
                    tables.count == 1 ? "First and only" : "First");
        if (middle != 0 && middle != tables.count - 1)
            print_table(&tables.tables[middle], middle + 1, "Middle");
        if (tables.count > 1)
            print_table(&tables.tables[tables.count - 1], tables.count,
                        "Last");
    }
    free_tables(&tables);
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

int verify_reconstruction(const double *coefficients, int degree,
                          const Root *roots, int count, double tolerance) {
    Root *current = calloc((size_t)degree + 1, sizeof(*current));
    Root *next = calloc((size_t)degree + 1, sizeof(*next));
    if (!current || !next) {
        fprintf(stderr, "Unable to allocate reconstruction workspace.\n");
        free(current);
        free(next);
        return 0;
    }

    current[0].real = coefficients[0];
    int current_degree = 0;
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j <= degree; ++j) next[j] = (Root){0.0, 0.0};
        for (int j = 0; j <= current_degree; ++j) {
            next[j].real += current[j].real;
            next[j].imag += current[j].imag;
            next[j + 1].real -= current[j].real * roots[i].real -
                                current[j].imag * roots[i].imag;
            next[j + 1].imag -= current[j].real * roots[i].imag +
                                current[j].imag * roots[i].real;
        }
        Root *swap = current;
        current = next;
        next = swap;
        ++current_degree;
    }

    double scale = 1.0;
    for (int i = 0; i <= degree; ++i)
        if (fabs(coefficients[i]) > scale) scale = fabs(coefficients[i]);
    double max_error = 0.0;
    printf("\nPolynomial reconstructed from roots:\n");
    for (int i = 0; i <= degree; ++i) {
        double error = hypot(current[i].real - coefficients[i], current[i].imag);
        if (error > max_error) max_error = error;
        printf("x^%d: reconstructed=% .10g %+.3ei, input=% .10g, error=%.3e\n",
               degree - i, current[i].real, current[i].imag,
               coefficients[i], error);
    }
    int verified = max_error <= tolerance * scale;
    printf("Maximum coefficient error: %.3e\n", max_error);
    printf("Reconstruction verification: %s\n", verified ? "PASS" : "FAIL");
    free(current);
    free(next);
    return verified;
}
