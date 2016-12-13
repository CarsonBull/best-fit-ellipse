//
// Find best fitting ellipse to input sample points.
// W. Cochran   wcochran@acm.org
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsm.h"

//
// Input sample points
//
int numPoints;
int capacity;
double (*points)[2];

//
// Implicit function for general ellipse.
//    F(x,y; B,C,D,E,F) = 0 when (x,y) on ellipse
// See http://en.wikipedia.org/wiki/Ellipse
// Equation 1 in project writeup.
// We have normized the implicit function for a general
// ellipse so so that A = 1, yielding 5 coefficients.
//
double F(double x, double y, 
	 double B, double C, double D, double E, double F) {
  // XXX
  return 0;
}

//
// Functions we are trying to minimize:
//  f(B,C,D,E,F) = sum over all points (x_i,y_i) of F(x_i,y_i; A,B,C,D,E,F)^2
// Equation 2 in project writeup.
// We are trying to find the best ellipse coefficents B, C, D, E, F
//
double f(const double x[]) {  // x = 5-D vector (holds 5 coefficients)
  // XXX
  return 100000;
}

int main(int argc, char *argv[]) {
  //
  // (1) Initialize points buffer.
  //
  numPoints = 0;
  capacity = 10;
  points = (double(*)[2]) malloc(2*sizeof(double)*capacity);

  //
  // (2) Read input points into buffer
  //
  double x, y;
  while (scanf("%lf %lf", &x, &y) == 2) {
    if (numPoints >= capacity) { // out of room: double capacity
      capacity *= 2;
      points = (double(*)[2]) realloc(points, 2*sizeof(double)*capacity);
    }
    points[numPoints][0] = x;
    points[numPoints][1] = y;
    numPoints++;
  }

  //
  // (3) Find bounding rectangle of input points.
  //
  // XXX

  //
  // (4) Initial guess at optimal ellipse will be a circle
  // inscribed in the bounding rectangle of the input points:
  //    f(x,y) = (x - cx)^2 + (y - xy)^2 - R^2
  //           = x^2 - + 0*x*y - 2*cx*x + y^2 - 2*cy*y + cx^2 + cy^2 - R^2
  // B = 0, C = 1, D = -2*cx, E = -2*xy, F = cx^2 + cy^2 - R^2
  //
  // XXX

  //
  // (5) Minimize f use Downhill Simplex Method
  //
  double xmin[5] = {0,0,0,0,0};
  int iters = 0;
  // XXX

  //
  // (6) Print ellipse coefficients
  //
  for (int i = 0; i < 5; i++)
    printf("%0.10f\n", xmin[i]);

  //
  // (7) Print number of iterations used and "average error".
  //
  const double error = sqrt(f(xmin)/numPoints);
  printf("iters=%d, error=%0.12f\n", iters, error);

  return 0;
}
