#include "heff.h"

// // Constructor
// heff::heff() {

// }

// Constructor with only a basis
heff::heff(const basis& set_basis ) : basis_(set_basis) {
  
  // TEST
  cout<<"Using "<<Eigen::nbThreads()<<" threads out of a maximum of "<<omp_get_max_threads()<<".\n\n";

  // Set the matrix PH0P_, independendent of e_ and a2Inv_
  setPH0P();
  
  // Flag that the effective hamiltonian has not yet been calculated
  initialized_=0;

}


// Constructor with initial parameters
heff::heff(const basis& set_basis, double set_e, double set_a2Inv ) : basis_(set_basis) {

  // TEST
  //cout<<"Using "<<Eigen::nbThreads()<<" threads out of a maximum of "<<omp_get_max_threads()<<".\n\n";

  // Set the matrix PH0P_, independendent of e3_ and a2Inv_
  setPH0P();

  // Set all the matrices from the input e, a2Inv
  set(set_e,set_a2Inv);

  // Set the initialized_ flag to true
  initialized_=1;

}


// Destructor
heff::~heff() {
}

// Updates the matrices
void heff::set(double set_e, double set_a2Inv) {

  // Do nothing if no update
  if( set_e==e3_ && set_a2Inv==a2Inv_) return;

  a2Inv_=set_a2Inv;

  // Update G0_ only if E changed or if it's the first call
  if(set_e!=e3_ || !initialized_ ) {
      e3_=set_e;
      setG0(e3_); 
  }
 
  // Set T12 
  setT12(e3_,a2Inv_);

  // Set V12eff
  setV12eff(e3_, a2Inv_);

  // Calculate Heff in the symmetric basis
  calcHeff(e3_);

  // Set the initialized_ flag to true
  initialized_=1;

}


// Calculates the HO Hamiltonian matrix in the P space symmetric basis
void heff::setPH0P() {

  // Set PH0P_ to zero
  PH0P_.setZero( basis_.statesDim(),basis_.statesDim() );

  MatrixXd fullPH0P;

  // Set fullPH0P to the identity so we can just loop over the diagonal
  fullPH0P.setIdentity( basis_.PbasisDim(),basis_.PbasisDim() );

  // Loop over diagonal and set entries
  for(unsigned int i=0; i<basis_.PbasisDim(); i++) fullPH0P(i,i)= basis_.getPState(i).energy() ;

  // Calculate matrix elements in the symmetric basis,
  // the Hamiltonian should be diagonal
  for(unsigned int i=0; i<basis_.statesDim(); i++) {
      PH0P_(i,i)=basis_.getState(i).transpose()*fullPH0P*basis_.getState(i); 
  }
		 
}

// Calculates the HO Green's function matrix at energy E in the P+Q space
// and stores it as a vector of the digonal components
void heff::setG0(double E) {

  // Resize G0_ to the correct length
  G0_.resize( basis_.fullBasisDim() );

  for(unsigned int i=0; i<basis_.fullBasisDim(); i++) {
    G0_(i)= 1/(E - basis_.getBasisState(i).energy());
  }
		 
}
// Calculates the two body t-matrix for a regulated two-body contact interaction at energy E and interaction strength a
// in the P+Q space
void heff::setT12(double E3, double a2Inv) {

  HOstate state,primeState;

  double temp;

  // Ensure that the matrix is the right size, if size is correct the eigen package knows and does nothing
  T12_.resize( basis_.fullBasisDim(), basis_.fullBasisDim() );

  for(unsigned int i=0; i<basis_.fullBasisDim(); i++) {
    
    if(basis_.getBasisState(i).l()==0 ) {

      primeState=basis_.getBasisState(i);
      
      // Loop only on lower left half
      for (unsigned int j=0; j<=i; j++ ) {
	
	state=basis_.getBasisState(j);

	// If N=N',L=L',l=l'=0 use 2-body t-matrix
	if( state.N()==primeState.N() && state.L()==primeState.L() && state.l()==0 ) {
	  
	  // Matrix elements of the two-body delta function in the HO basis
	  temp=exp (( LogGamma(primeState.n()+1.5)+LogGamma(state.n()+1.5)-LogGamma(primeState.n()+1.0)-LogGamma(state.n()+1) )/2.0 );

	  // Multiply by the coefficient of the delta function
	  temp*=2*ROOT2/PI;	  

	  temp*=C0( E3-2*state.N()-state.L()-1.5, a2Inv );	  
	  
	}
	
	// Otherwise the entries are zero
	else temp=0;

	T12_(i,j)= temp;
	T12_(j,i)= temp;
	
      } // End j loop
    } // End if (l'==0)
    
    // Inserts zeros when l'!=0
    else for (unsigned int j=0; j<basis_.fullBasisDim(); j++ ) {
	T12_(i,j)=0;
	T12_(j,i)=0;
      } // End else for

  } // End i loop
  
} // End setT12


// Required by setT12, defined as in Tom's notes
// Calculate using Busch formula, not Tom's formula for t12 (by which he means v12eff_)
double heff::C0(double E2, double a2Inv) {

  double c0;

  c0= -ROOT2;

  // The gamma function arguments can be negative, provide cases 
  if(.25-E2/2 > 0) c0*=exp( LogGamma( .75 - E2/2 ) - LogGamma(.25-E2/2) );
  else if ( .75-E2/2>0 ) c0*=exp(LogGamma( .75 -E2/2 ) + LogGamma(.75+E2/2) )*sin( PI*(.25-E2/2) )/PI; 
  else c0*=exp( LogGamma(.75+E2/2) - LogGamma(.25+E2/2) )*sin( PI*(.25-E2/2) )/sin( PI*(.75-E2/2) );

  c0+=a2Inv;

  return 1.0/c0;

}// End C0

// Calculate V12eff_ according to Tom's thesis eqn. 3.40
// MAKE SURE that T12 has been calculated already!
void heff::setV12eff(double E3, double a2Inv) {

  // Ensure that the matrix is the right size, if size is correct the eigen package knows and does nothing
  V12eff_.resize( basis_.fullBasisDim(), basis_.fullBasisDim() );

  // A temp matrix to store P space stuff
  MatrixXd temp(basis_.PbasisDim(), basis_.PbasisDim() );

  //This product shows up several times, best to just calculate it once
  MatrixXd G0T12=G0_.asDiagonal() * T12_;
  
  // temp= \Gamma_\infty = G0 T12 G0
  temp.noalias() = G0T12.block(0,0,basis_.PbasisDim(), basis_.PbasisDim() )*G0_.head( basis_.PbasisDim() ).asDiagonal();
 
  // temp= \Gamma_\infty + \Gamma_0 
  for(unsigned int i=0;i<basis_.PbasisDim();++i) temp(i,i) += G0_(i);
 
  // temp= ( \Gamma_0+\Gamma_\infty )^-1 
  // Inversion is done in the P basis, then projected back into the mixed symmetry basis.
  temp=basis_.fullToSym().transpose()*basis_.toSymBasis(temp).inverse()*basis_.fullToSym();

  // V12eff_= T12_ *G0_*( \Gamma_0+\Gamma_\infty )^-1*G0_*T12_
  // Try to build the answer without doing so many multiplications by zero, at the expense of clarity

  V12eff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim()).noalias()= G0T12.transpose().topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim()) * temp* G0T12.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim());
 
  V12eff_.topRightCorner(basis_.PbasisDim(),basis_.QbasisDim()).noalias()= G0T12.transpose().topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim()) *temp*G0T12.topRightCorner(basis_.PbasisDim(),basis_.QbasisDim());
  
  V12eff_.bottomRightCorner(basis_.QbasisDim(),basis_.QbasisDim()).noalias()= G0T12.transpose().bottomLeftCorner(basis_.QbasisDim(),basis_.PbasisDim()) * temp * G0T12.topRightCorner(basis_.PbasisDim(),basis_.QbasisDim());
  
  // Use self-adjointness
  V12eff_.bottomLeftCorner(basis_.QbasisDim(),basis_.PbasisDim()).noalias()= V12eff_.topRightCorner(basis_.PbasisDim(),basis_.QbasisDim()).transpose();

  // V12eff_= T12_-T12_ *G0_*( \Gamma_0+\Gamma_\infty )^-1*G0_*T12_
  V12eff_=T12_-V12eff_;
  
} // end setV12eff
 

// Heff  =   1/E3*(H0+(1+\Pi)*(1-V12eff G0 Q \Pi)^-1 V12eff)P
void heff::calcHeff(double E3) {

  // Heff_= G0 Q \Pi =  G0 * Q \Pi , this step should automatically resize the matrix
  Heff_.noalias() = G0_.asDiagonal() * basis_.getQPi();

  // Heff_= (1-V12eff G0 Q \Pi) in the P+Q space
  // This is a slow step, could probably do explicit multiplication using the block diagonal property of Heff_
  Heff_= MatrixXd::Identity(basis_.fullBasisDim() ,basis_.fullBasisDim() )-V12eff_.selfadjointView<Upper>()*Heff_;

  // Heff_= (1-V12eff_ G0 Q \Pi)^-1 in the P+Q space
  // This is a slow step
  Heff_=Heff_.inverse();

  // Heff_= P (1-V12eff_ G0 Q \Pi)^-1 V12eff_ P, we only need P space elements
  Heff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim())= Heff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim())*V12eff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim())+Heff_.topRightCorner(basis_.PbasisDim(),basis_.QbasisDim())*V12eff_.bottomLeftCorner(basis_.QbasisDim(),basis_.PbasisDim());

  // Project Heff_=Veff into symmetric P basis
  Heff_=basis_.toSymBasis(Heff_);

  // Multiply Heff_=Veff by (1+\Pi)
  Heff_*=3;

  //TEST
  //cout<<"V_eff in the symmetric basis is \n"<<Heff_<<endl<<endl;

  // Add H0
  Heff_+=PH0P_;
  
  // Divide by E3 to get eigenvalues of 1 rather than E3
  Heff_/=E3;

}

// Return the calculated effective Hamiltonian, with error checking to make sure it has been calculated
MatrixXd heff::getHeff() const {
    if(initialized_==0) cout<<"ERROR in heff.h: The effective Hamiltonian has not been calculated\n\n"; 
    return Heff_;
}
