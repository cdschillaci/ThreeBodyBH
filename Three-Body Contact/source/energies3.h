#ifndef energies3_H
#define energies3_H

#include "energies.h"
#include "heff3.h"

class energies3 : public energies {
 public:
  // Constructor for brute force initialization
  energies3(const basis& set_basis, double set_a2InvMin, double set_a2InvMax, double set_a2InvStep, double set_a3Min, double set_a3Max, double set_a3Step, double e3Min, double e3Max, double set_e3Step );
  // Constructor with a first e Seed input
  energies3(const basis& set_basis, double set_a2InvMin, double set_a2InvMax, double set_a2InvStep, double set_a3Min, double set_a3Max, double set_a3Step, double set_e3Step, double set_e3Range, const VectorXd& seed );
  ~energies3();
  
  void printPairs() {energies::printPairs("The best eigenvalue at each inverse a2, a3, and energy E is \n( 1/a2 , a3, e , best_lambda)\n");}; // Prints the eLambdaPairs_
  
 private:
  const double a2InvMin_,a2InvMax_,a2InvStep_;
  const double a3Min_,a3Max_,a3Step_;

  void Init(); // Instantiates Heff_ and sets eLambdaPairs to the correct size

  void calcAllPairs(double a2InvMin, double a2InvMax, double a2InvStep, double a3Min, double a3Max, double a3Step, double e3Step, double e3Range, const VectorXd& initSeed);

};

#endif
