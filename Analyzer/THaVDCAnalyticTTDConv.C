///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCAnalyticTTDConv                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaVDCAnalyticTTDConv.h"
#include <iostream>
#include <string>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TSpline.h"
#include <cmath>

ClassImp(THaVDCAnalyticTTDConv)


//_____________________________________________________________________________
THaVDCAnalyticTTDConv::THaVDCAnalyticTTDConv() : fIsSet(false)
{
  //Normal constructor
}

//_____________________________________________________________________________
THaVDCAnalyticTTDConv::THaVDCAnalyticTTDConv( Double_t vel ) :
  fDriftVel(vel), fIsSet(false)
{
  // Normal constructor 
}



//_____________________________________________________________________________
THaVDCAnalyticTTDConv::~THaVDCAnalyticTTDConv()
{
  // Destructor. Remove variables from global list.

}

//_____________________________________________________________________________
Double_t THaVDCAnalyticTTDConv::ConvertTimeToDist(Double_t time,
						  Double_t tanTheta,
						  Double_t *ddist)
{
  // Drift Velocity in m/s time in s
  // Return m 
  
  if( !fIsSet ) {
    Error( "THaVDCAnalyticTTDConv::ConvertTimeToDist", "Parameters not set. "
	   "Fix database." );
    return 1e38;
  }
  const Int_t max = 20000;
  const Int_t n = 200;
  
  Double_t Point_t[max];
  Double_t Point_y[max];
  Double_t Err_y[max];

  ifstream fout("$APEX/Thesis/ttd.dat");
    
  //Check the file    
  if (fout.is_open()){
        // File opened        
    cout << "It opened!" << endl;
    for(int i = 0; i < n ; i++){
      
      fout >> Point_t[i]; // read in the time component
      fout >> Point_y[i]; // read in the space component
      fout >> Err_y[i];   // read in the space error
           
      cout << Point_t[i] << " " << Point_y[i] << " "; // read in the time component
      cout << Err_y[i] << endl; // read in the space component
    }
  }

  TSpline3 *sp = new TSpline3("Cubic Spline", Point_t, Point_y, n, "b2e2", 0, 0);

  Double_t dist = sp -> Eval( time*pow(10,9) )*pow(10,-2);
  Double_t dist = fDriftVel * time
  
  Double_t unc  = sp -> Eval( fdtime*pow(10,9) )*pow(10,-2);  // watch uncertainty in the timing

  if (ddist) *ddist = unc;

  return dist;
  
}

//_____________________________________________________________________________
void THaVDCAnalyticTTDConv::SetParameters( const Double_t* A1,
					   const Double_t* A2, Double_t dtime )
{
  // Set coefficients of a1 and a2 4-th order polynomial and uncertainty
  // of drift time measurement

  for( int i=0; i<4; ++i ) {
    fA1tdcCor[i] = A1[i];
    fA2tdcCor[i] = A2[i];
  }  
  fdtime = dtime;
  fIsSet = true;
}

//_____________________________________________________________________________
void THaVDCAnalyticTTDConv::SetDefaultParam()
{
   // Set some reasonable defaults for the polynomial coefficients and 
   // drift time uncertainty. Applicable to Hall A VDCs.

   fA1tdcCor[0] = 2.12e-3;
   fA1tdcCor[1] = 0.0;
   fA1tdcCor[2] = 0.0;
   fA1tdcCor[3] = 0.0;
   fA2tdcCor[0] = -4.20e-4;
   fA2tdcCor[1] =  1.3e-3;
   fA2tdcCor[2] = 1.06e-4;
   fA2tdcCor[3] = 0.0;
  
   fdtime    = 4.e-9; // 4ns -> 200 microns
 }

///////////////////////////////////////////////////////////////////////////////
