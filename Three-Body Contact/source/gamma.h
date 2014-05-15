#ifndef GAMMA_H
#define GAMMA_H

#include <math.h>
#define PI 3.14159265358979

// Implement log(gamma(x)) using the Lanczos algorithm, based on Numerical Recipes
double LogGamma(double x);
// Handles negative arguments, 
// TO DO: Implement efficiency gains for small gamma over exp(LogGamma(x))
double Gamma(double x);

// Iterative implementations of factorial function
double LogFactorial(int x);
double Factorial(int x);

// Binomial coefficient
double LogBinomial(int a, int b);
double Binomial(int a, int b);

#endif
