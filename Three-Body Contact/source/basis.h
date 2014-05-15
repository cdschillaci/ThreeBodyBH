#ifndef basis_H
#define basis_H

#define CSV_ENTRY_LENGTH 25     //Space for negative sign+leading zero+period+up to four more zeros+17 digits+comma=25

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#include "HOstate.h"

using namespace std;
using namespace Eigen;

class basis {
public:

  // Constructs the physical basis using files in dir_name output from Mathematica.
  // Work in the restricted HO basis with center of mass energy E<=nmax+3
  // Also reads in a basis for the Q space (with lmax set by a Mathematica script) and the matrix Q\Pi
  // Make sure that set_dirName ends in a "/", i.e. "../states/"
  basis(int set_nmax,int set_Qnmax,int set_J,string set_dirName);

  // Destructor, no effect
  ~basis();
  
  // Changes the basis of available states
  void set(int set_nmax,int set_Qnmax,int set_J);

  // Check if the basis population worked as expected
  bool isValid() const {return isValid_;};

  // These functions return values of private class variables
  int nmax() const {return nmax_;};
  int J() const {return J_;};
  unsigned int PbasisDim() const {return basisDim_;};
  unsigned int QbasisDim() const {return QbasisDim_;};
  unsigned int statesDim() const {return statesDim_;};
  unsigned int fullBasisDim() const {return basisDim_+QbasisDim_;}; // The full basis of P + Q HO states

  HOstate getBasisState(unsigned int i) const; // Gets HO states from the P+Q basis 
  HOstate getPState(unsigned int i) const; // Gets HO states from the P basis
  HOstate getQState(unsigned int i) const; // Gets HO states from the Q basis

  const MatrixXd& getPPi() const {return PPi_;}; // Returns the matrix P\Pi
  const MatrixXd& getQPi() const {return QPi_;}; // Returns the matrix (1-P)\Pi
  const MatrixXd& fullToSym() const {return fullToSym_;}; //Projects the mixed symmetry P basis onto the symmetric P states

  VectorXd getState(unsigned int i) const {return states_.col(i);}; // Gets the representations of symmetric three body states in the P basis

  void printPBasis() const {printBasis(0);};
  void printQBasis() const {printBasis(1);};
  void printStates() const;
  void printPi() const; 
  void printPPi() const; 
  void printQPi() const; 

  MatrixXd toSymBasis(const MatrixXd& m) const; // Convert a matrix m in P basis to the symmetric basis 

protected:
  int nmax_; // POOR NOMENCLATURE. Highest energy allowed is E=nmax_+3
  int Qnmax_; // This is the maximum Q space state energy to include
  int J_; // States will all have total angular momentum J_

  string dirName_; // When reading in states, holds the directory path

  void Init(); // Calls the other functions when constructor or set are called

  bool isValid_; // Goes to zero if there is an error reading the basis

  unsigned int basisDim_; // Number of HO basis kets in the P space
  unsigned int QbasisDim_; // Number of HO basis kets we will use in the Q space
  unsigned int statesDim_; // Number of physical states in the P space with correct permutation symmetry

  Matrix<HOstate, Dynamic, 1> PBasis; // The array of basis kets is handled using the Eigen routines
  Matrix<HOstate, Dynamic, 1> QBasis; // The array of basis kets is handled using the Eigen routines
  MatrixXd PPi_; // The matrix P\Pi
  MatrixXd QPi_; // The matrix Q\Pi
  MatrixXd states_; // Contains the representation in the P basis of states with proper permutation symmetry
  MatrixXd fullToSym_;

  void readBasis(bool Q); // Reads in basis states. For Q==False, reads in P states. For Q==TRUE, reads in Q states
  void readSymBasis(); // Reads in the (anti)symmetric three particle basis states. 
  void readPi(); // Reads in the matrices Q\Pi=(1-P) and P\Pi 
  void popFullToSym(); //

  void printBasis(bool Q) const;

};

#endif
