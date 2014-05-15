#include <iostream>
#include <iomanip>
#include <math.h>

#include "basis.h"
#include "energies3.h"

//TEST
#include <omp.h>

using namespace std;

int main () {

  int nmax=4;
  int Qnmax=200;
  basis workingBasis(nmax,Qnmax,0,"../states/");

  //TEST
  double t1,t2;

  cout<<"nmax= "<<nmax<<" and Qnmax= "<<Qnmax<<endl;
  cout<<"The Q basis contains "<<workingBasis.QbasisDim()<<" states\n\n";

  // Abort if basis is invalid
  if(workingBasis.isValid()==0) {
    cout<<"Basis is not valid\n";
    return 0;
  }

  t1=omp_get_wtime();

  //   This object's constructor does all the work. 
  //   energies3(const basis& set_basis, double set_a2InvMin, double set_a2InvMax, double set_a2InvStep, double set_a3Min, double set_a3Max, double set_a3Step, double e3Min, double e3Max, double set_e3Step );
  energies3 newEnergies3( workingBasis, 0.0, 0.0, 1/3.0, //a2Inv
			  -0.05, 0.05, 0.005, //a3
			  7.01, 9.01, .0075*sqrt(2) ); // e3

  
  t2=omp_get_wtime();
  cout<<"energies3 ran in "<<t2-t1<<" seconds\n";

  cout<<"Writing result to file \"../data/PairsBya2Inv.csv\"\n\n";
  newEnergies3.writePairs("../data/PairsBya2Inv.csv");

  }
