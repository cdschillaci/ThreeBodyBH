#include "basis.h"

// Constructs the physical basis using files in dir_name output from Mathematica
basis::basis(int set_nmax,int set_Qnmax,int set_J,string set_dirName) : nmax_(set_nmax),Qnmax_(set_Qnmax),J_(set_J),dirName_(set_dirName) {
  
  Init();

}

// Destructor
basis::~basis() {
}

// Changes the basis of available states
void basis::set(int set_nmax,int set_Qnmax,int set_J) {
  
  nmax_=set_nmax;
  J_=set_J;

  Init();

}

// Read in everything
void basis::Init() {

  isValid_=1; // The basis is valid by default, isValid_ will be set to 0 in case of an error

  readBasis(0); // Reads in the P basis in HO space
  readSymBasis(); // Reads in the symmetric states in the P basis
  popFullToSym(); // Creates a matrix to project from full P basis to symmetric P basis

  readBasis(1); // Reads in the Q basis in HO space, with lmax set by the mathematica script
  
  readPi(); // Reads in the matrices PPi and QPi

}

// Read in a basis prepared by Mathematica
// Q==0--> read in P basis
// Q==1--> read in Q basis
void basis::readBasis(bool Q) {

  Matrix<HOstate, Dynamic, 1> currentBasis;

  HOstate temp(0,0,0,0,J_);

  string filename; // Used to hold the full filename
  ostringstream tempstream; // Needed to insert numbers into the filename
  ifstream file; 

  stringstream ss;
  char buff[20];

  int minshell,maxshell; // Range of shells to read in
  int currentDim = 0; // Set the number of basis states to zero

  // Set range of shells to read if reading P
  if(!Q) {
    minshell=0;
    maxshell=nmax_;}

  // Set range of shells to read if reading Q
  else {
    minshell=nmax_+2;
    maxshell=Qnmax_;}

  for(int i=minshell; i<=maxshell;i+=2) {

    // Create a string with the name of the file dirName_+"Basis"+i+".csv"
    filename=dirName_;
    if(Q) filename+="Q"; //Add the Q part of the filename if necessary
    filename+="Basis";
    tempstream<<i;
    filename+=tempstream.str();
    filename+=".csv";

    tempstream.str(""); // Clear tempstream for reuse

    // Open the file 
    file.open(filename.c_str());

    // Error handling
    if( !file.is_open() ) {
      cout<<"ERROR: the file "<<filename<<" could not be opened\n\n";
      isValid_=0;
      break;
    }

    // On each loop, reads a line in and assigns to buff. Exits when reaches eof
    while( file.getline( buff, 20 ) ) {
      
      // Then assign whole line to the stream ss
      ss<<buff;
      // Set the quantum numbers of temp from the line, must be integers of digits or less
      ss.getline(buff, 6 , ','); // Read ss up until a "," or 6 digits and store in buff
      temp.set_n( atoi(buff) );  // Convert buff to integer and set n quantum number for temp
      ss.getline(buff, 6 , ','); // and so on...
      temp.set_l( atoi(buff) ); 
      ss.getline(buff, 6 , ',');
      temp.set_N( atoi(buff) ); 
      ss.getline(buff, 6 , ',');
      temp.set_L( atoi(buff) ); 
      ss.getline(buff, 6 , ',');
      temp.set_J( atoi(buff) ); 

      // Add to the HO basis
      if( temp.IsValid() ) {
	currentDim++;
	currentBasis.conservativeResize( currentDim );
	currentBasis(currentDim-1)=temp;
      }

      // If the quantum numbers don't work, warn the user that something is up
      else cout<<"Bad HO state in "<<filename<<endl;

      // This copies an empty string into ss, erasing the old contents
      ss<<" ";
      // This clears the eof flag.
      ss.clear();

    } // End while

    file.close();

  } // End loop over different files

  if(!Q) {
    basisDim_=currentDim;
    PBasis.resize( basisDim_ );
    PBasis=currentBasis;
  }
  
  else {
    QbasisDim_=currentDim;
    QBasis.resize( QbasisDim_ );
    QBasis=currentBasis;
  }

}

// This function reads in symmetric states expressed in the HO basis
void basis::readSymBasis() {
 
  string filename; // Used to hold the full filename
  ostringstream tempstream; // Needed to insert numbers into the filename
  ifstream file; 

  stringstream ss;
  char *buff; // A c-string with dynamic length for buffering data from file
  int bufflength; // Store the length here

  // Store the number of zeros needed to pad in the subspaces with E<Eshell
  int beginPad;

  ArrayXi subspaceLengths;

  statesDim_=0; // Start with zero symmetric states

  //////////////////////////////////////////////////////////////////////////////
  // In this section we populate the ArrayXi subspaceLengths with the size of //
  // each subspace in the basis                                               //
  //////////////////////////////////////////////////////////////////////////////

  // Start by sizing the array and initializing all entries to zero
  subspaceLengths.resize( nmax_/2+1 );
  for(int i=0;i<nmax_/2+1;++i) subspaceLengths(i)=0;

  // Loop over all the states, incrementing subspaceLengths at the correct shell
  // corresponding to each state
  for(unsigned int i=0;i<basisDim_;i++) subspaceLengths( (PBasis(i).energy()-3)/2 )+=1; 

  ////////////////////////////////////////////////////////////////////////////
  // In this section we process the files                                   //
  ////////////////////////////////////////////////////////////////////////////  

  // Loop over all of the files, e.g. "symBasis0.csv"
  for(int nshell=0; nshell<=nmax_;nshell+=2) {

    //  Create a string with the name of the file dirName_+"symBasis"+nshell+".csv"
    filename=dirName_;
    filename+="symBasis";
    tempstream<<nshell;
    filename+=tempstream.str();
    filename+=".csv";
    
    tempstream.str(""); // Clear tempstream for reuse
    
    // Open the file 
    file.open(filename.c_str());

    // Error handling
    if( !file.is_open() ) {
      cout<<"ERROR: the file "<<filename<<" could not be opened\n\n";
      isValid_=0;
      break;
    }

    // Store the number of states with E<Eshell
    beginPad=0;
    for(int j=0;j<nshell/2;j++) beginPad+=subspaceLengths(j);

    // Make the buffer more than long enough for negative sign+leading zero+period+up to four more zeros+17 digits+comma
    // and add a full set of extra spaces to be safe
    bufflength=(subspaceLengths(nshell/2)+1)*CSV_ENTRY_LENGTH;
    buff = (char*) malloc( bufflength*sizeof(char) );
    
    // On each loop, reads a line in and assigns to buff. Exits when reaches eof
    while( file.getline( buff, bufflength ) ) {
      
      // Then assign whole line to the stream ss
      ss<<buff;
 
      // Increment the number of states
      statesDim_++;

      // Resize the states_ matrix to hold the new state
      states_.conservativeResize(basisDim_,statesDim_);

      // Pad with initial zeros
      for(int j=0;j<beginPad;j++) states_.col(statesDim_-1)(j)=0;

      // Transfer coefficients from file to column
      for(int j=beginPad;j<beginPad+subspaceLengths(nshell/2);j++) {
	ss.getline( buff, CSV_ENTRY_LENGTH, ',');  // Read ss up until a "," or 25 digits and store in buff
	states_.col(statesDim_-1)(j)= atof( buff ); // Set column vector entries
      } //End for

      // Pad with final zeros
      for(int j=beginPad+subspaceLengths(nshell/2);j<int(basisDim_);j++) states_.col(statesDim_-1)(j)=0;
      
      ss<<" ";  // This copies an empty string into ss, erasing the old contents
      ss.clear();  // This clears the eof flag.

    } // End while loop over file
    
    free(buff);

    file.close();

  }  
} // End of readSymBasis

// This function reads in Q\Pi expressed in the QBasis
// Since QPi_ is block diagonal, this can be done one shell at a time
// for greater flexibility
// TO DO: Fix the horrible indexing for subspaceLengths!
void basis::readPi() {
 
  string filename; // Used to hold the full filename
  ostringstream tempstream; // Needed to insert numbers into the filename
  ifstream file; 

  stringstream ss;
  char *buff; // A c-string with dynamic length for buffering data from file
  int bufflength; // Store the length here
  MatrixXd temp; // Copy data here first, then to QPi_
  int counter; // Which row of the matrix is being read

  // Store the number of zeros needed to pad in the subspaces with E<Eshell
  int beginPad;

  ArrayXi subspaceLengths;
  
  QPi_.setZero( basisDim_+QbasisDim_, basisDim_+QbasisDim_ ); // Pi is a matrix of zeros of the correct size

  //////////////////////////////////////////////////////////////////////////////
  // In this section we populate the ArrayXi subspaceLengths with the size of //
  // each subspace in the Q basis                                             //
  //////////////////////////////////////////////////////////////////////////////

  // Start by sizing the array and initializing all entries to zero
  subspaceLengths.setZero( Qnmax_/2+1 );

  // Loop over all the states, incrementing subspaceLengths at the correct shell
  // corresponding to each state
  for(unsigned int i=0;i<fullBasisDim() ;i++) {
    if(i<basisDim_) subspaceLengths( (PBasis(i).energy()-3)/2 )+=1; 
    else subspaceLengths( (QBasis(i-basisDim_).energy()-3)/2 )+=1;
  }

  ////////////////////////////////////////////////////////////////////////////
  // In this section we process the files                                   //
  ////////////////////////////////////////////////////////////////////////////  

  // Loop over all of the files, e.g. "Pi52.csv"
  for(int nshell=0; nshell<=Qnmax_ ;nshell+=2) {

    //  Create a string with the name of the file dirName_+"Pi"+nshell+".csv"
    filename=dirName_;
    if(nshell<=nmax_) filename+="PPi";
    else filename+="QPi";
    tempstream<<nshell;
    filename+=tempstream.str();
    filename+=".csv";
    
    tempstream.str(""); // Clear tempstream for reuse
    
    // Open the file 
    file.open(filename.c_str());

    // Error handling
    if( !file.is_open() ) {
      cout<<"ERROR: the file "<<filename<<" could not be opened\n\n";
      isValid_=0;
      break;
    }

    // Store the number of states with E<Eshell
    beginPad=0; 
    for(int j=0; j<nshell/2; j++) beginPad+=subspaceLengths(j);

    // Make the buffer more than long enough for negative sign+leading zero+period+up to four more zeros+17 digits+comma
    // and add a full set of extra spaces to be safe
    bufflength= subspaceLengths( nshell/2 )*CSV_ENTRY_LENGTH ;
    buff = (char*) malloc( bufflength*sizeof(char) );

    // Make the matrix temp the correct size to take the input
    temp.resize(subspaceLengths( nshell/2), subspaceLengths( nshell/2 ) );
    
    counter=0; // Reset the row counter to zero

    // On each loop, reads a line in and assigns to buff. Exits when reaches eof
    while( file.getline( buff, bufflength ) ) {
      
      // Then assign whole line to the stream ss
      ss<<buff;

      // Transfer coefficients from file to temp matrix
      for(int j=0;j<subspaceLengths( nshell/2 );j++) {
	ss.getline( buff, CSV_ENTRY_LENGTH, ',');  // Read ss up until a "," or 25 digits and store in buff
	temp(counter,j)= atof( buff ); // Set temp matrix entries
      } //End for
      
      ss<<" ";  // This copies an empty string into ss, erasing the old contents
      ss.clear();  // This clears the eof flag.

      counter++; // Increment the row counter
      // Some error handling in case counter gets bigger than it should be
      if( counter>subspaceLengths(nshell/2) ) {
	cout<<"ERROR: variable counter in readPi outside of range\n";
	isValid_=0;
	break;
      }

    } // End while loop over file

    // Copy temp to a diagonal block in Pi_
    QPi_.block(beginPad,beginPad,subspaceLengths( nshell/2 ),subspaceLengths( nshell/2 ) )=temp;
    
    // Cleanup, needed either for another run or at the end of the runs
    free(buff);
    file.close();
    
  } // End loop over nshell 

  // Save P\Pi
  PPi_=QPi_.topLeftCorner(  basisDim_, basisDim_ );

  // Subtract off P\Pi
  QPi_.topLeftCorner(  basisDim_, basisDim_ ) *= MatrixXd::Identity(basisDim_, basisDim_ ) - fullToSym_.transpose()*fullToSym_;

  
} // End of readPi

// Returns an object of type HOstate corresponding to the state |i> in the P+Q basis
HOstate basis::getBasisState(unsigned int i) const {

  if(i>=basisDim_+QbasisDim_) cout<<"ERROR in basis::getBasisState! The P+Q basis is not that large\n";
  
  if(i<basisDim_) return PBasis(i);
  else return QBasis(i-basisDim_);

}

// Returns an object of type HOstate corresponding to the state |i> in the P basis
HOstate basis::getPState(unsigned int i) const {

  if(i>=basisDim_) cout<<"ERROR in basis::getPState! The P basis is not that large\n";

  return PBasis(i);

}

// Returns an object of type HOstate corresponding to the state |i> in the Q basis
HOstate basis::getQState(unsigned int i) const {

  if(i>=QbasisDim_) cout<<"ERROR in basis::getQState! The Q basis is not that large\n";

  return QBasis(i);

}


void basis::printBasis(bool Q) const {

  Matrix<HOstate, Dynamic, 1> currentBasis;

  if(Q) currentBasis=QBasis;
  else currentBasis=PBasis;

  int length=currentBasis.rows();
  
  currentBasis(0).printHeader();
  cout<<"----------------\n";

  for(int i=0;i<length;i++) {

    currentBasis(i).printKet();

  }
  cout<<"\n";
}


void basis::printPi() const { 

  MatrixXd temp=MatrixXd::Identity( fullBasisDim(), fullBasisDim()) ;

  temp.topLeftCorner( basisDim_, basisDim_ ) = PPi_;
  temp.bottomRightCorner( QbasisDim_, QbasisDim_ ) = QPi_;

  cout<<"The matrix Pi_ is \n"<<temp<<"\n";

}

void basis::printPPi() const {  cout<<"The matrix PPi_ is \n"<<PPi_<<"\n";} 

void basis::printQPi() const {  cout<<"The matrix QPi_ is \n"<<QPi_<<"\n";} 

void basis::printStates() const {

  cout<<"The basis of HO states with correct permutational symmetry contains "<<statesDim_<<" states\n";
  cout<<"----------------------------------------------------------------------------- \n\n";

  for(unsigned int i=0;i<statesDim_;i++) {
    
    for(unsigned int j=0;j<basisDim_;j++) {
      
      if(states_.col(i)(j)!=0) {
	cout<<setw(10)<<showpoint<<setprecision(6)<<states_.col(i)(j);
	PBasis(j).printKet();	
      }

    } // End j loop
    
    cout<<"\n";

  } // End i loop
	     
} // End of printStates 

void basis::popFullToSym() {
  
  fullToSym_.resize(statesDim_,basisDim_);
  for(unsigned int i=0; i<statesDim_; ++i) {
    for(unsigned int j=0; j<basisDim_; ++j) {
      fullToSym_(i,j)=states_.col(i)(j);
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////
// Converts a matrix m in the P or full basis to the symmetric basis             //
//   -Assumes the the first PbasisDim_ entries are in the P basis!               //
//   -If the matrix is of size statesDim_, assumes it's already in sym basis     //
///////////////////////////////////////////////////////////////////////////////////
MatrixXd basis::toSymBasis (const MatrixXd& m) const {

  if( m.rows()==statesDim_ ) return m;
  
  // tempM1 is the P space block of m
  MatrixXd tempM1=m.topLeftCorner( basisDim_, basisDim_ );

  // tempM2 will be m projected into the symmetric basis
  MatrixXd tempM2;
  tempM2.resize(statesDim_,statesDim_);

  VectorXd tempV; // Stores m*state for efficiency in loop below

  // Put m into states basis
  for(unsigned int j=0;j<statesDim_;j++) {  
    tempV.noalias()=tempM1*getState(j);   
    for(unsigned int i=0;i<statesDim_;i++) {
      tempM2(i,j)=getState(i).transpose()*tempV;
    } // End i loop
  } // End j loop
  
  return tempM2;

}
