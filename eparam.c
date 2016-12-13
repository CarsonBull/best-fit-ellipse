//
// Given coefficients B,C,D,E,F for implicit function
// for an ellipse (B^2 < 4C):
//    F(x,y; B,C,D,E,F) = x^2 + B*xy + C*y^2 + D*x + E*y + F = 0
// we find the parameters (a, b, phi, Xc, Yc) for the parameterized
// version of an ellipse:
//    X(t) = Xc + a cos(t) cos(phi) - b sin(t) sin(phi)
//    Y(t) = Yc + a cos(t) sin(phi) + b sin(t) cos(phi)
// where
//    0 <= t < 2*pi
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsm.h"

//
// Coefficients for implicit function describing ellipse;
// Implements Equation 1.
//
double B, C, D, E, F;
double G(double x, double y) {
  return x*x + B*x*y + C*y*y + D*x + E*y + F;
}

//
// Function to minimize.
// Implement Equation 6 from writeup.
// x[0],x[1] = a, b
// x[2] = phi
// x[3], x[4] = xc, yc
//
double g(const double x[]) {
  // XXX
  return 100000;
}

int main(int argc, char *argv[]) {

  //
  // (1) Get coefficients of implicite equation for ellipse.
  //
  double coeffs[5];
  for (int i = 0; i < 5; i++)
    if (scanf("%lf", &coeffs[i]) != 1) {
      fprintf(stderr, "error reading coefficient %d!\n", i);
      exit(-1);
    }

  B = coeffs[0];
  C = coeffs[1];
  D = coeffs[2];
  E = coeffs[3];
  F = coeffs[4];

  //
  // (2) Get initial guess at ellipse parameters.
  //
  // XXX

  //
  // (3) Minimize g use Downhill Simplex Method
  //
  double xmin[5] = {0,0,0,0,0};
  int iters = 0;
  // XXX

  //
  // (4) Output ellipse parameters and error info.
  //
  for (int i = 0; i < 5; i++)
    printf("%0.10f\n", xmin[i]);

  printf("iterations = %d, fmin = %0.10f\n", iters, sqrt(g(xmin)));

  return 0;
}
