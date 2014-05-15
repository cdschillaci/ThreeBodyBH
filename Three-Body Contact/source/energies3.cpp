#include "energies3.h"

// Constructor when no seed is provided to guide the initial energy sweep
energies3::energies3(const basis& set_basis, double set_a2InvMin, double set_a2InvMax, double set_a2InvStep, double set_a3Min, double set_a3Max, double set_a3Step, double e3Min, double e3Max, double set_e3Step ) :energies( set_basis, set_e3Step, 0), a2InvMin_(set_a2InvMin),a2InvMax_(set_a2InvMax),a2InvStep_(set_a2InvStep), a3Min_(set_a3Min),a3Max_(set_a3Max),a3Step_(set_a3Step) {

  Init();

  // Create an initial set of seeds by stepping from eMin to eMax. Note that eRange is set to zero by constructor
  VectorXd startSeed=initSeed(e3Min,e3Max,set_e3Step);

  calcAllPairs(a2InvMin_, a2InvMax_, a2InvStep_, a3Min_, a3Max_, a3Step_, e3Step_, e3Range_, startSeed);
  
}

// Constructor when a seed is provided to guide the initial energy sweep
energies3::energies3(const basis& set_basis, double set_a2InvMin, double set_a2InvMax, double set_a2InvStep, double set_a3Min, double set_a3Max, double set_a3Step, double set_e3Step, double set_e3Range, const VectorXd& seed ) : energies( set_basis, set_e3Step, set_e3Range), a2InvMin_(set_a2InvMin),a2InvMax_(set_a2InvMax),a2InvStep_(set_a2InvStep), a3Min_(set_a3Min),a3Max_(set_a3Max),a3Step_(set_a3Step) {

  Init();

  calcAllPairs(a2InvMin_, a2InvMax_, a2InvStep_, a3Min_, a3Max_, a3Step_, e3Step_, e3Range_, seed);
  
}

// Destructor
energies3::~energies3() {
  
  delete Heff_; // Delete Heff_
  
}

void energies3::Init() {
  Heff_=new heff3(basis_);  // instantiate Heff_

  eLambdaPairs_.resize(0,4); // 4 columns per row store (a2Inv, a3, e, bestLambda)

}


//////////////////////////////////////////////////////////////////////////////////
// Calls calcpair for                                                           //
//      Heff  =   1/E*(H0+(1+\Pi)*(1-t12 G0 Q \Pi)^-1 t12)                      //
// at each of the points in the aInv,e grid specified                           //
// The results are saved to eLambdaPairs                                        //
//                                                                              //
// TO DO: Prevent from running with e=integer (divergent G0 matrix elts)        //
//////////////////////////////////////////////////////////////////////////////////
void energies3::calcAllPairs(double a2InvMin, double a2InvMax, double a2InvStep, double a3Min, double a3Max, double a3Step, double e3Step, double e3Range, const VectorXd& initSeed) {

  VectorXd newSeed=initSeed; // Will be used to store new seed arrays as generated. Should probably be an array
  Vector2d contactStrengths; // 

  #if USE_SEED
  unsigned int nNewEnergies=0; // Number of new energies to use in genSeed
  #endif

  totalStepEstimate_=( a2InvMax-a2InvMin+a2InvStep )/a2InvStep*( a3Max-a3Min+a3Step )/a3Step*initSeed.rows();

  // If the stepsize is zero we get infinite loops
  if(a3Step==0) cout<<"ERROR: a3Step=0.0\n";
  if(a2InvStep==0) cout<<"ERROR: a2InvStep=0.0\n";

  for(double a2Inv=a2InvMin; a2Inv<=a2InvMax; a2Inv+=a2InvStep ) {
    
    contactStrengths(0)=a2Inv;

    for(double a3=a3Min; a3<=a3Max; a3+=a3Step ) {

      contactStrengths(1)=a3;

      #if USE_SEED
      // If not using initial seed array, generate a new seed array and save the old number of good eLambdaPairs
      if(a2Inv!=a2InvMin) {
	newSeed=genSeed(nNewEnergies,e3Step);
	nNewEnergies=eLambdaPairs_.rows();
      }//End if
      #endif
      
      eSweep(contactStrengths, newSeed); 

      #if USE_SEED
      nNewEnergies=eLambdaPairs_.rows()-nNewEnergies;
      #endif

    } // End a3 loop
    
  } // end a2Inv loop

}
