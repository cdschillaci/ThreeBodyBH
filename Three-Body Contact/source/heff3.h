#ifndef heff3_H
#define heff3_H

#include "heff.h"

class heff3 : public heff {
 public:
  // Constructor with only a basis
  heff3(const basis& set_basis );
  // Constructor that initializes everything
  heff3(const basis& set_basis, double set_e3, double set_a2Inv, double set_a3Inv);
  ~heff3();
  
  virtual void set(double e3, const VectorXd& contactStrengths) {set(e3,contactStrengths(0),contactStrengths(1));}; // Updates the matrices for a new energy, a2Inv, a3
  
 protected:
 
  double a3_; // Inverse scattering length and three body energy

  MatrixXd V123_; // The three body contact interaction, independent of e3_
  MatrixXd V123eff_; // The effective three body t-matrix as in Luu eqn 34

  void setV123(double a3); // Independent of e3_
  void setV123eff(double E3, double a2Inv, double a3); // The arguments don't really accomplish anything except to remind you of the dependence, could be removed
  
  void set(double set_e3, double set_a2Inv, double set_a3); // Updates the matrices for a new energy, scattering length, and three body contact strength
  void calcHeff(double E3); // Make sure G0 and T12 have been updated before calling!

    
};

#endif
