# Bairstow's Method

A command-line C program for finding the roots of a polynomial using Bairstow's method. The program accepts polynomial degrees from 2 through 5, prints intermediate polynomial calculations, and displays the calculated roots.

## Requirements

- A C compiler such as GCC
- The math library (`libm`), linked with `-lm` when using GCC

## Build

From the project directory, compile with GCC:

```sh
gcc main.c -o bairstow -lm
```

On Windows with MinGW-w64, run the generated executable with:

```powershell
.\bairstow.exe
```

On Linux or macOS:

```sh
./bairstow
```

## Use

1. Enter the polynomial's degree (2–5).
2. Enter each coefficient when prompted, in descending power order, starting with the coefficient of the highest power and ending with the constant term.
3. Enter a positive tolerance for convergence.

For example, for the polynomial `x^2 - 3x + 2`, enter degree `2`, then coefficients `1`, `-3`, and `2`, followed by a tolerance such as `0.0001`.

## Notes

- Coefficients and tolerance are entered as decimal numbers.
- The starting values for Bairstow's iteration are `u = 1` and `v = 1`.
- The results depend on the selected tolerance and the convergence of the iteration for the entered polynomial.
