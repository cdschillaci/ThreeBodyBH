#include "energies.h"

// Constructor to initialize const variables in base class
energies::energies(const basis& set_basis, double set_e3Step, double set_e3Range) :e3Step_(set_e3Step), e3Range_(set_e3Range), basis_(set_basis), progressCounter_(0) {
  
}

// Destructor
energies::~energies() {
}

// Prints the pairs using cout
void energies::printPairs(string dataDescription) const {

  cout<<dataDescription<<endl;

  for( unsigned int i=0; i<eLambdaPairs_.rows(); i++) {
    cout<<eLambdaPairs_.row(i)<<"\n";
  }

}

// Creates a csv file with data from eLambdaPairs on each row
void energies::writePairs(string filename) const {
  
  ofstream file;
  file.open(filename.c_str());

  for(unsigned int i=0; i<eLambdaPairs_.rows(); i++) {
    for(unsigned int j=0; j<eLambdaPairs_.cols(); j++) {
      file<<eLambdaPairs_(i,j);
      if(j!=eLambdaPairs_.cols()-1) file<<","; //Comma required by csv format, no comma after last data point
    } // End j loop
    file<<"\n"; // Line break after each row
  } // End i loop

  file.close();

}

//////////////////////////////////////////////////////////////////////////////////
// Calculates the eigenvalue of mm closest to 1                                 //
//  -Returns MAX_EIGENVALUE if no good eigenvalue found, best eigenvalue        //
//   if successful                                                              //
//////////////////////////////////////////////////////////////////////////////////
double energies::calcPair(const MatrixXd& mm) {

  double best=MAX_EIGENVALUE; // Start with a large number, eigenvalue must be closer to one than this
  
  // Calculate the eigenvalues
  EigenSolver<MatrixXd> eigensolver(mm);

  // Find eigenvalue closest to 1
  for( unsigned int i=0; i<basis_.statesDim(); i++ ) {
    
    if( fabs(eigensolver.eigenvalues()(i).real()-1) < fabs(best-1) ) best=eigensolver.eigenvalues()(i).real();
    
  }// end i loop
  
  // Print an error 
  if(best==MAX_EIGENVALUE) cout<<"Error, no good eigenvalues\n";
    
  // If a good eigenvalue was found, return it or return MAX_EIGENVALUE if no good eigenvalue was found
  return best;

}


#if USE_SEED
// Will someday intelligently create a vector of energies around which to sweep at the next small increment of a2 or a3
VectorXd energies::genSeed(unsigned int nEnergies, double eStep) {

  cout<<"MEGA ERROR, FUNCTIONALITY NOT IMPLEMENTED\n\n";

  return temp;

}
#endif

// Given a vector of energies, calculate the eigenvalue of Heff/e closest to 1 at each and save a triplet aInv/
void energies::eSweep(const VectorXd& contactStrengths, const VectorXd& seed) {

  double bestLambda; // Stores the best eigenvalue each time Heff is calculated

  unsigned int numPairs=eLambdaPairs_.rows(); // The total number of triplets (a2Inv,e,lambda_best) or quads (a2Inv,a3,e,lambda_best) saved. Initialize to the current number of rows in eLambdaPairs_

  unsigned int nSeeds=seed.rows(); // Number of seed energies. 

  // Loops over all seed energies
  for(unsigned int seedNumber=0; seedNumber<nSeeds; seedNumber++) {

    // Sweeps over all energies within eRange of the seed energy, with resolution eStep
    for( double e=seed(seedNumber)-e3Range_ ; e<=seed(seedNumber)+e3Range_; e+=e3Step_ ) {
      
      // Calculate Heff in the symmetric basis
      Heff_->set(e,contactStrengths);
      
      numPairs++;

      bestLambda=calcPair( Heff_->getHeff() );
      
      // Increase size of eLambdaPairs, unecessary if there was no good eigvenvalue last loop
      if( numPairs>eLambdaPairs_.rows() ) eLambdaPairs_.conservativeResize(numPairs,Eigen::NoChange);

      // If a good eigenvalue was found, add to set
      if(bestLambda!=MAX_EIGENVALUE) {
	for(unsigned int i=0; i<eLambdaPairs_.cols()-2;i++) eLambdaPairs_.row(numPairs-1)(i)=contactStrengths(i); // Loops over the contact strengths
	  
	eLambdaPairs_.row(numPairs-1)(eLambdaPairs_.cols()-2)=e;
	eLambdaPairs_.row(numPairs-1)(eLambdaPairs_.cols()-1)=bestLambda;

      }

      // Otherwise decrement numPairs. 
      else numPairs--;

      progressCounter_++;
      if( progressCounter_ % 1 == 0) cout<<"Percent complete = "<< setiosflags( ios::fixed )<<setprecision(2)<<100.0*progressCounter_/totalStepEstimate_<<"%\n";
      writePairs("../data/PairsBya2Inv.csv");

    } // end e loop

    // If the last pass did not find a good eigenvalue, remove the extra row in eLambdaPairs_
    if( numPairs<eLambdaPairs_.rows() ) eLambdaPairs_.conservativeResize(numPairs,Eigen::NoChange);

  } // end seedNumber loop 

  // If the numPairs is not the original estimate, resize eLambdaPairs
  if (numPairs!=eLambdaPairs_.rows() ) eLambdaPairs_.conservativeResize(numPairs,Eigen::NoChange);

}


 // Create an initial set of seeds by stepping from eMin to eMax. Note that eRange should be set to zero by constructor or this won't make sense
VectorXd energies::initSeed(double e3Min, double e3Max, double e3Step) { 

  VectorXd temp;
  unsigned int i=0;

  if(e3Range_!=0) cout<<"ERROR: attempting to create a full energy sweep with e3Range_ nonzero.\n";

  temp.resize( floor( (e3Max-e3Min)/e3Step_ ) +1);
  for( double e=e3Min; e<=e3Max; e+=e3Step ) {
    temp(i)=e;
    i++;
  }

  return temp;

}
