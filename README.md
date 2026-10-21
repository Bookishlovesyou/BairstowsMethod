# Bairstow's Method

A command-line C program for finding the roots of a polynomial using Bairstow's method. The program accepts polynomial degrees from 2 through 5, prints intermediate polynomial calculations, and displays the calculated roots.

## Requirements

- A C compiler such as GCC
- The math library (`libm`), linked with `-lm` when using GCC

## Build

The program has two C source files: `main.c` (the entry point) and
`polynomial.c` (the polynomial solver). Compile both from the project
directory with GCC:

```sh
gcc main.c polynomial.c -o bairstow -lm
```

You can also compile each source file separately and then link the object files:

```sh
gcc -c main.c -o main.o
gcc -c polynomial.c -o polynomial.o
gcc main.o polynomial.o -o bairstow -lm
```

`polynomial.h` is included by the source files; it is not compiled separately.

On Windows with MinGW-w64, the output executable has an `.exe` extension:

```powershell
.\bairstow.exe
```

On Linux or macOS:

```sh
./bairstow
```

## Use

1. Enter the polynomial's degree (2–5).
2. Enter each coefficient when prompted, in descending power order: from the highest power through the constant term.
3. Enter a positive tolerance for convergence.

For example, for `x^2 - 3x + 2`, enter:

```text
2
1
-3
2
0.0001
```

The program prints the Bairstow iteration tables and then lists the roots. Real roots are printed as real numbers; complex roots are printed as a real part plus or minus an imaginary part.

## Notes

- The leading coefficient must be nonzero. Coefficients and tolerance must be finite numbers, and the tolerance must be positive.
- Bairstow iteration starts with `r = 0` and `s = 0`, and tries alternate starting values if it does not converge.
- The results depend on the selected tolerance. Some polynomials may not converge with the available starting values; in that case, try a different tolerance.
