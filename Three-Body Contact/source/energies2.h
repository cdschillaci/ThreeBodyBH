#ifndef energies2_H
#define energies2_H

#include "energies.h"
#include "heff.h"

class energies2 : public energies {
 public:
  // Constructor for brute force initialization
  energies2(const basis& set_basis, double set_aInvMin, double set_aInvMax, double set_aInvStep, double eMin, double eMax, double set_eStep );
  // Constructor with a first e Seed input
  energies2(const basis& set_basis, double set_aInvMin, double set_aInvMax, double set_aInvStep, double set_eStep, double set_eRange, const VectorXd& seed );
  ~energies2();

  void printPairs() {energies::printPairs("The best eigenvalue at each inverse a2 and energy E is \n( 1/a2 , e , best_lambda)\n");}; // Prints the eLambdaPairs_

 private:
  const double a2InvMin_,a2InvMax_,a2InvStep_;

  void Init(); // Instantiates Heff_ and sets eLambdaPairs to the correct size

  void calcAllPairs(double aInvMin, double aInvMax, double aInvStep, double eStep, double eRange, const VectorXd& initSeed);

};

#endif
