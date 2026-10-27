# Bairstow's Method

This program finds the roots of a polynomial using Bairstow's method. It accepts polynomials of degree 2 through 5, prints the roots, and checks that they reconstruct the original polynomial.

## Example

For the polynomial `x^2 - 3x + 2`, enter the degree, coefficients from highest to lowest power, and a positive tolerance:

```text
2
1
-3
2
0.00000001
```

Expected roots and verification:

```text
Roots:
+------+------------------------------+
| Root | Value                        |
+------+------------------------------+
| x1   | 2                            |
| x2   | 1                            |
+------+------------------------------+

Reconstruction verification: PASS
```
