Bairstow Method - Test Cases
============================

Run the program and enter each test case's inputs one value per line. For
successful cases, the program also prints a Bairstow iteration table; the
expected output below shows the final Roots section.

Test 1: Quadratic with two real roots
-------------------------------------
Input:
2
1
-3
2
1e-8

Polynomial: x^2 - 3x + 2
Expected output:
Roots:
x1 = 2
x2 = 1

Test 2: Cubic with three real roots
-----------------------------------
Input:
3
1
-6
11
-6
1e-8

Polynomial: x^3 - 6x^2 + 11x - 6
Expected output:
Roots:
x1 = 1.999999997
x2 = 1.000000001
x3 = 3.000000002

Test 3: Quartic with real and complex roots
-------------------------------------------
Input:
4
1
0
-3
0
-4
1e-8

Polynomial: x^4 - 3x^2 - 4
Expected output:
Roots:
x1 = -0 + 1.000000001i
x2 = -0 - 1.000000001i
x3 = 2
x4 = -2

Test 4: Quintic with real and complex roots
-------------------------------------------
Input:
5
1
-1
-3
3
-4
4
1e-8

Polynomial: x^5 - x^4 - 3x^3 + 3x^2 - 4x + 4
Expected output:
Roots:
x1 = -0 + 1.000000001i
x2 = -0 - 1.000000001i
x3 = 2
x4 = -2
x5 = 1

Test 5: Degree below the supported range
----------------------------------------
Input:
1

Expected error:
Degree must be between 2 and 5.
Expected result: program exits with failure.

Test 6: Zero leading coefficient
--------------------------------
Input:
2
0
1
1
1e-8

Expected error:
The leading coefficient cannot be zero.
Expected result: program exits with failure.

Test 7: Zero tolerance
----------------------
Input:
2
1
-3
2
0

Expected error:
Tolerance must be a positive number.
Expected result: program exits with failure.
