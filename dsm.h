#ifndef DSM_H
#define DSM_H

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
			  double xmin[], double *fmin);
#endif // DSM_H
