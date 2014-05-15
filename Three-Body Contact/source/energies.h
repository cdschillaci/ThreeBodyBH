#ifndef energies_H
#define energies_H

#define MAX_EIGENVALUE 1234

#define USE_SEED 0 // If true, will use initial energy sweep to seed future energy sweeps even if no initial seed is provided

#include <iostream>
#include <fstream>
#include <string.h>

#include "basis.h"
#include "heff.h"
#include "heff3.h"

class energies {
 public:

  void writePairs(string filename) const;

 protected:

  // Constructor 
  energies(const basis& set_basis, double set_e3Step, double set_e3Range);
  ~energies();

  const double e3Step_,e3Range_;

  const basis& basis_;

  unsigned int progressCounter_,totalStepEstimate_; // These might need to be moved into the derived classes

  heff * Heff_; // Use a pointer so that this can refer to  heff or heff3

  Matrix<double, Dynamic, Dynamic> eLambdaPairs_; // Number of columns depends on whether we sweep in a3 or just a2Inv

  void printPairs(string dataDescription) const;

  double calcPair(const MatrixXd& mm);
  void calcAllPairs(double aInvMin, double aInvMax, double aInvStep, double eStep, double eRange, const VectorXd& initSeed);
  void eSweep(const VectorXd& contactStrengths, const VectorXd& seed); 

 // Create an initial set of seeds by stepping from eMin to eMax. Note that eRange should be set to zero by constructor or this won't make sense
  VectorXd initSeed(double e3Min, double e3Max, double e3Step);

  #if USE_SEED
  VectorXd genSeed(unsigned int nEnergies, double eStep);
  #endif

};

#endif
