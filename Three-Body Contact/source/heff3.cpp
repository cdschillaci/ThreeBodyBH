#include "heff3.h"

// Constructor with only a basis
heff3::heff3(const basis& set_basis ) : heff(set_basis) {
  
  // Note that we have set the matrix PH0P_ and the basis_ using the base constructor

  // Flag that the effective hamiltonian has not yet been calculated
  initialized_=0;

}


// Constructor with initial parameters
heff3::heff3(const basis& set_basis, double set_e, double set_a2Inv, double set_a3 ) : heff(set_basis) {

  // Note that we have set the matrix PH0P_ and the basis_ using the base constructor

  // Set all the matrices from the input e, aInv
  set(set_e,set_a2Inv,set_a3);

  // Set the initialized_ flag to true
  initialized_=1;

}


// Destructor
heff3::~heff3() {
}

// Updates the matrices
void heff3::set(double set_e3, double set_a2Inv, double set_a3) {

  // Do nothing if no update
  if( set_e3==e3_ && set_a2Inv==a2Inv_ && set_a3==a3_  ) return;

  // If a3=0, don't waste time! 
  if( set_a3==0 ) {
    heff::set(set_e3,set_a2Inv); // Use inherited set function
    a3_=0.0; // Make sure that this object is keeping track of a3_
    initialized_=1; // Heff_ has been set
    return; 
  }

  // Update G0_ only if E changed or if it's the first call
  if(set_e3!=e3_ || !initialized_ ) setG0(set_e3); 
 
  // These don't update if only a3 changes
  if( set_e3!=e3_ || set_a2Inv==a2Inv_ || !initialized_ ) {
   
    // Set T12 
    setT12(set_e3,set_a2Inv);
    
    // Set V12eff
    setV12eff(set_e3, set_a2Inv);   

  } // End if
 
  // Always need to updated V123_
  setV123(set_a3);

  // Calculate V123eff in the symmetric basis
  setV123eff(set_e3, set_a2Inv, set_a3);
  
  // Calculate Heff in the symmetric basis
  calcHeff(set_e3);

  // Set the initialized_ flag to true
  initialized_=1;

  // Update energy, interaction strengths
  e3_=set_e3;
  a2Inv_=set_a2Inv;
  a3_=set_a3;
  

}

void heff3::setV123(double a3) {

  // Ensure V123_ is the correct size
  V123_.resize(basis_.fullBasisDim() ,basis_.fullBasisDim() );

  HOstate state,primeState;

  double temp;

  for(unsigned int i=0; i<basis_.fullBasisDim(); i++) {
    
    if( basis_.getBasisState(i).l()==0 && basis_.getBasisState(i).L()==0 ) {

      primeState=basis_.getBasisState(i);
      
      // Loop only on lower left half
      for (unsigned int j=0; j<=i; j++ ) {
	
	state=basis_.getBasisState(j);

	// If N=N',L=L'=l=l'=0 use 2-body t-matrix
	if( state.N()==primeState.N() && state.L()==primeState.L() && state.l()==primeState.L() ) {
	  
	  // Matrix elements of the three-body delta function in the HO basis
	  temp=exp (( LogGamma(primeState.n()+1.5)+LogGamma(state.n()+1.5)+LogGamma(primeState.N()+1.5)+LogGamma(state.N()+1.5)-LogGamma(primeState.n()+1.0)-LogGamma(state.n()+1) -LogGamma(primeState.N()+1.0)-LogGamma(state.N()+1) )/2.0 );

	  // Multiply by the coefficient of the three-body delta function
	  temp*=a3;	  	  
	  
	}
	
	// Otherwise the entries are zero
	else temp=0;

	V123_(i,j)= temp;
	V123_(j,i)= temp;
	
      } // End j loop
    } // End if (l'==0 && L'==0)
    
    // Inserts zeros when l'!=0 or L'!=0
    else for (unsigned int j=0; j<basis_.fullBasisDim(); j++ ) {
	V123_(i,j)=0;
	V123_(j,i)=0;
      } // End else for

  } // End i loop

} // End setV123
 
void heff3::setV123eff(double E3, double a2Inv, double a3) {

  MatrixXd temp;
  
  // temp = (1-V12eff_ G0_ ) V123/3
  temp.noalias()=(MatrixXd::Identity( basis_.fullBasisDim() ,basis_.fullBasisDim() )+V12eff_*G0_.asDiagonal())*V123_/3.0;
  
  // V123eff_ = 1-(1+V12eff_ G0_ ) V123/3 G0
  V123eff_.noalias()= MatrixXd::Identity( basis_.fullBasisDim() ,basis_.fullBasisDim() ) - temp*G0_.asDiagonal();

  // V123eff_ = ( 1+(1-V12eff_ G0_ ) V123/3 G0 )^-1
  V123eff_=V123eff_.inverse();

  // temp = V12eff_ + (1-V12eff_ G0_ ) V123/3
  temp+=V12eff_;

  // V123eff_ = ( 1+(1-V12eff_ G0_ ) V123/3 G0 )^-1*( V12eff_+(1+V12eff_ G0_ ) V123/3 )
  V123eff_=V123eff_*temp;  

} // End setV123eff


// Heff  =   1/E*(H0+(1+\Pi)*(1-V123eff G0 Q \Pi)^-1 V123eff)P
void heff3::calcHeff(double E) {

  // Heff_= G0 Q \Pi =  G0 * Q \Pi , this step should automatically resize the matrix
  Heff_.noalias() = G0_.asDiagonal() * basis_.getQPi();

  // Heff_= (1-V12eff G0 Q \Pi) in the P+Q space
  // This is a slow step, could probably do explicit multiplication using the block diagonal property of Heff_
  Heff_= MatrixXd::Identity(basis_.fullBasisDim() ,basis_.fullBasisDim() )-V123eff_.selfadjointView<Upper>()*Heff_;

  // Heff_= (1-V12eff_ G0 Q \Pi)^-1 in the P+Q space
  // This is a slow step
  Heff_=Heff_.inverse();

  // Heff_= P (1-V12eff_ G0 Q \Pi)^-1 V12eff_ P, we only need P space elements
  Heff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim())= Heff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim())*V123eff_.topLeftCorner(basis_.PbasisDim(),basis_.PbasisDim())+Heff_.topRightCorner(basis_.PbasisDim(),basis_.QbasisDim())*V123eff_.bottomLeftCorner(basis_.QbasisDim(),basis_.PbasisDim());

  // Project Heff_=Veff into symmetric P basis
  Heff_=basis_.toSymBasis(Heff_);

  // Multiply Heff_=Veff by (1+\Pi)
  Heff_*=3;

  //TEST
  //  cout<<"V_eff in the symmetric basis is \n"<<Heff_<<endl<<endl;

  // Add H0
  Heff_+=PH0P_;
  
  // Divide by E to get eigenvalues of 1 rather than E
  Heff_/=E;

}
