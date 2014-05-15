#ifndef heff_H
#define heff_H

#define PI 3.14159265358979
#define ROOT2 1.41421356237310

#include <iostream>
#include <fstream>
#include <string.h>

#include <Eigen/Dense>
#include <Eigen/Core>

#include <omp.h> // Should add something so that this is not included unless threading

#include "HOstate.h"
#include "basis.h"
#include "gamma.h"

class heff{
 public:
  // Constructor with only a basis
  heff(const basis& set_basis );
  // Constructor that initializes everything
  heff(const basis& set_basis, double set_e, double set_a2Inv);
  virtual ~heff();
  
  MatrixXd getG0() const; // Returns the last calculated HO Green's function, with error handling
  MatrixXd getHeff() const; // Returns the effective Hamiltonian, with error handling
  
  virtual void set(double e3, const VectorXd& contactStrengths) {set(e3,contactStrengths(0));}; // Allows for consistent interaction with derived classes
  
 protected:
  bool initialized_; // Flag to make sure that the matrices have been initialized
  
  double a2Inv_,e3_; // Inverse scattering length and three body energy
  
  const basis& basis_; // The basis being used
  
  VectorXd G0_; // The harmonic oscillator Green's function in the P+Q space
                // This is a vector of the diagonal coefficients. It MUST be used in
                // matrix multiplication as G0_.asDiagonal() !

  MatrixXd PH0P_; // The harmonic oscillator Hamiltonian in the P space symmetric basis
  MatrixXd T12_; // The t-matrix in the P+Q space
  MatrixXd V12eff_; // The effective 2+1 interaction in the P+Q space
  MatrixXd Heff_; //  1/E*(H0+(1+\Pi)*(1-t12 G0 Q \Pi)^-1 t12) P |\psi> = Heff P |\psi>= P|\psi>
  
  virtual void set(double e3, double a2Inv); // Updates the matrices for a new energy and scattering length 
  
  void setPH0P();
  void setG0(double E3);
  void setT12(double E3, double a2Inv);
  void setV12eff(double E3, double a2Inv);
  
  double C0(double E2, double a2Inv); // E2 is the 2-body energy, E2=E-2*N-L-1.5
 
  void calcHeff(double E); // Make sure G0 and T12 have been updated before calling!

    
};

#endif
