
gROOT -> Reset();

#include <iostream>
#include <string>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"

using namespace std;

const Long64_t n = 200;

Double_t Point_t[2000];
Double_t Point_y[2000];
Double_t Err_y[2000];

void array(){

    ifstream fout("ttd.dat");
    
        // Check the file    
        if (fout.is_open())
        {
            // File opened        
            cout << "It opened!" << endl;
                
           
            for(int i = 0; i < n ; i++)
            {
        
                fout >> Point_t[i]; // read in the time component
                fout >> Point_y[i]; // read in the space component
                fout >> Err_y[i];   // read in the space error
               
                //cout << Point_t[i] << " " << Point_y[i] << " "; // read in the time component
                //cout << Err_y[i] << endl; // read in the space component
            

            
           }
           for(int j = 0; j < n ; j++)
           {
           
           cout << Point_t[j] << ", ";
           //cout << Point_y[j] << ", ";
            
           }
           
       }
}
