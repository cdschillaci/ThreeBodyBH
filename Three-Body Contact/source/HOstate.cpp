#include "HOstate.h"

// Empty constructor
HOstate::HOstate() {}

// Constructor
HOstate::HOstate(int set_n, int set_l,int set_N, int set_L,int set_J) : n_(set_n),l_(set_l),N_(set_N),L_(set_L),J_(set_J) {
  Init();
}
// Destructor
HOstate::~HOstate() {
}

// Initializer
void HOstate::Init() {
  
  valid_=checkState(); // Sanity check
  
}

// Prints | nl, NL, J>
void HOstate::printHeader() const {
  cout<<"| n l, N L, J >\n";
}

// Prints the values for | nl, NL, J>
void HOstate::printKet() const {
  cout<<"| "<<n_<<" "<<l_<<", "<<N_<<" "<<L_<<", "<<J_<<" >\n";
}

// Returns the HO energy 2n+l+2N+L
int HOstate::energy() const {
  return 2*n_+l_+2*N_+L_+3;
}

// Check that quantum numbers are positive and J is possible, displaying errors
void HOstate::IsValidVerbose() const {
  
  if (n_<0) cout<<"WARNING: n_ initialized to "<<n_<<"\n";
  if (l_<0) cout<<"WARNING: l_ initialized to "<<l_<<"\n";
  if (N_<0) cout<<"WARNING: N_ initialized to "<<N_<<"\n";
  if (L_<0) cout<<"WARNING: L_ initialized to "<<L_<<"\n";
  
  if (!IsTriangular(l_,L_,J_)) cout<<"WARNING: J_ initialized to "<<J_<<" with l_="<<l_<<" and L="<<L_<<"\n";

}

// These functions set values of private class variables and perform sanity check
void HOstate::setAll(int new_n, int new_l, int new_N, int new_L, int new_J) {
  n_=new_n;
  l_=new_l;
  N_=new_N;
  L_=new_L;
  J_=new_J;

  valid_=checkState();

};

void HOstate::set_n(int new_n) {
  n_=new_n;
  valid_=checkState();
}
void HOstate::set_l(int new_l) {
  l_=new_l;
  valid_=checkState();
}
void HOstate::set_N(int new_N) {
  N_=new_N;
  valid_=checkState();
}
void HOstate::set_L(int new_L) {
  L_=new_L;
  valid_=checkState();
}
void HOstate::set_J(int new_J) {
  J_=new_J;
  valid_=checkState();
}

// Tests whether the triad (l1,l2,L) satisfies the triangular inequality |l1-l2|<L<l1+l2
bool HOstate::IsTriangular (int l1,int l2,int L) const {
  return ( fabs(l1-l2)<=L && L<=l1+l2 ) ? 1:0;
}

// Check that quantum numbers are positive and J is posible
bool HOstate::checkState() {
  

  // Reduce error checking overhead when NDEBUG defined
  #ifndef NDEBUG

  if (!IsTriangular(l_,L_,J_)) return 0;

  if (n_<0 || l_<0 || N_<0 || L_<0 || J_<0) return 0;

  #endif 

  return 1;

}

bool operator== (HOstate &state1, HOstate &state2) {

  return (state1.n_==state2.n_ && state1.l_==state2.l_ && state1.N_==state2.N_ && state1.L_==state2.L_ && state1.J_==state2.J_ );
 
}
bool operator!= (HOstate &state1, HOstate &state2) {
  return !(state1==state2);
}
