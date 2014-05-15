/*******************************************************************************
This class creates harmonic oscillator states of the form | n l, N L, J >

-Upon construction, it verifies that the quantum numbers are physical. If they 
 are not a warning is printed.

-User visualization of state is handled via printHeader() and printKet()

-The CoM energy for a state in HO units can be obtained from the energy() 
 function

*******************************************************************************/


#ifndef HOstate_H
#define HOstate_H

#include <iostream>
#include <math.h>

using namespace std;

class HOstate {

public: 
  HOstate();
  HOstate(int set_n, int set_l,int set_N, int set_L,int set_J);
  ~HOstate();

  void Init(); // Initializer
  
  void printHeader() const; // Prints | nl, NL, J>
  void printKet() const; // Print a nicely formatted ket expression;

  // Return energy of the state
  int energy() const;

  // Print out any problems with current quantum numbers
  bool IsValid() const {return valid_;}
  void IsValidVerbose() const;

  // These functions set values of private class variables
  void setAll(int new_n, int new_l, int new_N, int new_L, int new_J);
  void set_n(int new_n);
  void set_l(int new_l);
  void set_N(int new_N);
  void set_L(int new_L);
  void set_J(int new_J);
 
  // These functions return values of private class variables
  int n() const {return n_;}
  int l() const {return l_;}
  int N() const {return N_;}
  int L() const {return L_;}
  int J() const {return J_;}

  // Comparison operators
  friend bool operator== (HOstate &state1, HOstate &state2);
  friend bool operator!= (HOstate &state1, HOstate &state2);
  
protected:
  int n_;
  int l_;
  int N_;
  int L_;
  int J_;
  
  bool valid_;

  // Utility that tests whether the triad (l1,l2,L) satisfies the triangular inequality |l1-l2|<L<l1+l2
  bool IsTriangular (int l1,int l2,int L) const;
  bool checkState();
  
};

#endif
