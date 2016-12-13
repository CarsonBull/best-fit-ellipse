#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * Linear interpolation between 2 N-dimensional vectors:
 *    P(t) = A + (B - A)*t
 * Note:
 *    P(0) = A
 *    P(1) = B
 *    P(0.5) = (A + B)/2 = midpoint
 *    P(2) = 2B - A = B + (B - A) = point A reflected around B.
 * Safe to alter A or B "in place" (i.e., A and P as well as B and P
 * can be the same vector).
 */
void lerp(unsigned N, double t, 
          const double A[], const double B[], 
          double P[]) {
  for (unsigned k = 0; k < N; k++)
    P[k] = A[k] + (B[k] - A[k])*t;
}

/*
 * Clever, efficient hack to swap the contents of two integer variables
 * without involving any third temporary value.
 */
#define SWAP(a,b) (a ^= b, b ^= a, a ^= b)

/*
 * Implementation of the "Downhill Simplex Method" described in 
 * J. A. Nelder and R. Mead
 * A Simplex Method for Function Minimization 
 * The Computer Journal 1965 7: 308-313; doi:10.1093/comjnl/7.4.308 
 *
 * Input:
 *   N : dimension (number of arguments for f).
 *   f : function we are trying to minimize.
 *   p : initial guess = N-dimensional point used to seed algorithm.
 *   delta : "characteristic scale" = distance betweem vertices
 *     in initial simplex.
 *   ftolerance : difference in f() values used to determine convergence.
 *   maxIters: cap on number of iterations or simplexes to try.
 * Output:
 *   xmin : result = N-dimensional point that is our guess at the minimum.
 *   fmin : f(xmin)
 * Returns:
 *   number of iterations used.
 *   If iterations were exhausted it might be worth trying again
 *   using the returned xmin value as the need seed p.
 */
int downhillSimplexMethod(unsigned N, double (*f)(const double x[]),
                          const double p[], double delta, 
                          double ftolerance, unsigned maxIters,
                          double xmin[], double *fmin) {
  const double tiny = 1e-10;  /* a small positive number > machine-epsilon */
  const double Ninverse = 1.0/N;

  /*
   * x holds N+1 vertices of simplex in N dimensions.
   * y holds value of function f at simplex vertices.
   * initialize simplex as follows:
   *   x[0] = p  (initial guess)
   *   x[i] = x[0] + delta*ei  i=1..N
   * where ei is the unit vector pointing in direction of ith axis.
   */
  double x[N+1][N];
  double y[N+1];
  const unsigned bytesPerPoint = N*sizeof(double);
  memcpy(x[0],p,bytesPerPoint);
  y[0] = f(x[0]);
  for (unsigned i = 1; i <= N; i++) {
    memcpy(x[i],p,bytesPerPoint);
    x[i][i-1] += delta;
    y[i] = f(x[i]);
  }

  unsigned best, worst; // index of "best" and "worst" simplex vertex

#ifdef PLOT_SIMPLEX
  char fname[30];
  sprintf(fname, "simplex%u.plot", N);
  FILE *fp = fopen(fname, "w");
  if (fp == NULL) {perror(fname); exit(-1);}
#endif

  /*
   * Main loop.
   *   For each iteration, create a new simplex from the old one 
   *   that (hopefully) walks its way to a local minimum.
   */
  unsigned iters;
  for (iters = 0; iters < maxIters; iters++) {

#ifdef PLOT_SIMPLEX
    fprintf(fp, "# iteration %u\n", iters);
    for (unsigned i = 0; i <= N; i++) {
      for (unsigned k = 0; k < N; k++)
        fprintf(fp, "%0.10f ", x[i][k]);
      fprintf(fp, "\n");
    }
    if (N == 2) { // replicate first point for triangles
      for (unsigned k = 0; k < N; k++)  
        fprintf(fp, "%0.10f ", x[0][k]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
#endif

    /*
     * Get indices for "worst" and "best"
     * vertices in simplex:
     *   y[best] <= ... <= y[i] <= ... <= y[worst]
     */
    best = 0;
    worst = 1;
    if (y[worst] < y[best])
      SWAP(worst, best);
    for (unsigned i = 2; i <= N; i++)
      if (y[i] < y[best])
        best = i;
      else if (y[i] > y[worst])
        worst = i;

    /*
     * If the relative difference between the function values
     * a the best and worst vertices is less than the prescribed 
     * tolerance, then we're done. 'tiny' is used to make sure we
     * don't divide by zero.
     */
    const double relativeFuncDifference = 
      2*fabs(y[worst] - y[best])/(fabs(y[worst]) + fabs(y[best]) + tiny);
    if (relativeFuncDifference <= ftolerance)
      break;

    /*
     * Compute the centroid of the face of the simplex not
     * containing the worst point.
     */
    double centroid[N];
    bzero(centroid, bytesPerPoint);
    for (unsigned i = 0; i <= N; i++)
      if (i != worst)
        for (unsigned k = 0; k < N; k++)
          centroid[k] += x[i][k];
    for (unsigned k = 0; k < N; k++)
      centroid[k] *= Ninverse;

    /*
     * Reflect worst point about simplex face not
     * containing the worst point.
     */
    double reflection[N];
    lerp(N, 2, x[worst], centroid, reflection);
    double y_reflection = f(reflection);

    /*
     * Three cases:
     * (1) Major improvement : reflected point better than
     *   the best. We replace the worst vertex with the
     *   newfound best. To push our luck, we extend the
     *   reflection out farther and see if we get an
     *   even better result.
     * (2) Reflection sucks: reflection no better than the
     *   worst vertex. We look for a better point via
     *   a "one-dimensional contraction" which moves the worst
     *   point half-the-distance towards the centroid.
     *   If this still yields a poor result, we contract the
     *   entire simplex towards the best point.
     * (3) Minor improvement : reflected point not the best, but
     *   better than some of the other simplex vertices.
     *   We replace the worst vertex with refected vertex.
     */
    if (y_reflection <= y[best]) {  // major improvement!

#ifdef PLOT_SIMPLEX
      fprintf(fp, "# reflection...major improvement\n");
#endif

      /*
       * Replace our worst simplex point with reflection.
       */
      y[worst] = y_reflection;
      memcpy(x[worst],reflection,bytesPerPoint);

      /*
       * Let's extend the reflection even farther.
       * If we still get an improvement, we'll again replace
       * the worst value with the extended-reflected vertex.
       */
      lerp(N, 2, centroid, reflection, reflection);
      y_reflection = f(reflection);
      if (y_reflection <= y[best]) {
#ifdef PLOT_SIMPLEX
        fprintf(fp, "# extended reflection...\n");
#endif
        y[worst] = y_reflection;
        memcpy(x[worst],reflection,bytesPerPoint);
      }

    } else if (y_reflection >= y[worst]) {  // reflection sucks bad!

      /*
       * One dimensional contraction of worst point towards centroid,.
       */
      double contraction[N];
      lerp(N, 0.5, x[worst], centroid, contraction);
      double y_contraction = f(contraction);

      if (y_contraction >= y[worst]) {  // 1-D contraction didn't help

#ifdef PLOT_SIMPLEX
        fprintf(fp, "# contraction to best...\n");
#endif

        /*
         * Contract entire simplex towards best vertex.
         */
        for (unsigned i = 0; i <= N; i++)
          if (i != best) {
            lerp(N, 0.5, x[i], x[best], x[i]);
            y[i] = f(x[i]);
          }

      } else {  // 1-D contraction improved things

#ifdef PLOT_SIMPLEX
        fprintf(fp, "# 1-D contraction...\n");
#endif
        /*
         * Replace worst vertex with contracted vertex.
         */
        y[worst] = y_contraction;
        memcpy(x[worst],contraction,bytesPerPoint);     

      }

    } else { // minor improvement

#ifdef PLOT_SIMPLEX
        fprintf(fp, "# reflection -- minor improvement...\n");
#endif
      /*
       * Replace our worst simplex point with reflection.
       */
      y[worst] = y_reflection;
      memcpy(x[worst],reflection,bytesPerPoint);

    }
    
  }

#ifdef PLOT_SIMPLEX
      fclose(fp);
#endif

  /*
   * We give up since we have exhausted the specified
   * number of iterations.
   */
  memcpy(xmin, x[best], bytesPerPoint);
  *fmin = y[best];
  return iters;
}

#ifdef TEST_DSM

/*
 * f(x,y) = exp(|x-5| + |y-13|)
 * has a minumum at (5,13).
 */
double func(const double x[]) {
  return exp(fabs(x[0] - 5) + fabs(x[1] - 13));
}

int main(void) {
  double guess[2] = {0, 0};
  double xmin[2], fmin;
  int iters = downhillSimplexMethod(2, func, guess, 1.0, 0.0001, 100,
                                  xmin, &fmin);
  printf("iters = %d; f(xmin = (%0.5f, %0.5f)) = %0.5f\n", 
         iters, xmin[0], xmin[1], fmin);
  return 0;
}

#endif // TEST_DSM

#ifdef TEST_DSM_3D

/*
 * f(x,y,z) = exp(|x-5| + |y-13| + |z+4|)
 * has a minumum at (5,13,-4).
 */
double func(const double x[]) {
  return exp(fabs(x[0] - 5) + fabs(x[1] - 13) + fabs(x[2] + 4));
}

int main(void) {
  double guess[3] = {0, 0, 0};
  double xmin[3], fmin;
  int iters = downhillSimplexMethod(3, func, guess, 1.0, 0.0001, 1000,
                                  xmin, &fmin);
  printf("iters = %d; f(xmin = (%0.5f, %0.5f, %0.5f)) = %0.5f\n", 
         iters, xmin[0], xmin[1], xmin[2], fmin);
  return 0;
}

#endif // TEST_DSM_3D
