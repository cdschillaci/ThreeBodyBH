#include "gamma.h"

#include <iostream>
using namespace std;

// Uses the Lanczos approximation with g=671/128 and n=15 to calculate log(Gamma(x))
// Requires Re(x)>0
double LogGamma(double x) {

  int i;
  double y,temp,ser;
  static double cof[14]={57.1562356658629235,-59.5979603554754912,14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5}; // Coefficients in series defined by NR as c_1 through c_14

  if (x<=0) {
    cout<<"LogGamma error for x="<<x<<". Lanczos algorithm valid only for right half plane.\n\n";
    return 0;
  }

  y=x;
  temp=x;
  temp+=5.24218750000; // This is the rational number 671/128, defined in NR as \gamma+1/2
  temp=(x+0.5)*log(temp)-temp;
  
  ser=0.99999999999999709; // This is c_0 from the expansion as defined in NR
  for(i=0;i<14;i++) ser+=cof[i]/++y;
  
  return temp+log(2.5066282746310005*ser/x); // The decimal is sqrt(2*\pi)

}

// This function returns the gamma function of the argument for both negative and positive numbers.
// Gives an error for x=0 or negative integer
double Gamma(double x) {

  if (x<=0) {
    
    if( x-int(x)!=0 ) return PI / sin( PI*(1-x) ) / exp(LogGamma( (1-x) ));
    
    //Error handling for negative integers
    else {
      cout<<"Gamma("<<x<<") is infinite \n";
      return 0;
    }
  }
    
  else return exp(LogGamma(x));

}

double LogFactorial(int x) {

  double log_fact=0;

  if (x<0) {
    cout<<"("<<x<<")! is undefined.\n";
    return 0;
  }

  for (int i=1;i<=x;++i) log_fact+=log(i);

  return log_fact;

}

double Factorial(int x) {

  double fact=x;

  if (x<0) {
    cout<<"("<<x<<")! is undefined.\n";
  }

  for (int i=1;i<x;++i) fact*=i;

  return fact;

}

// Binomial coefficient
double LogBinomial(int n, int k) {
  return LogFactorial(n)-LogFactorial(k)-LogFactorial(n-k);
}
double Binomial(int n, int k) {
  return exp(LogBinomial(n,k));
}


