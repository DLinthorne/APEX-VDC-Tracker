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

//Double_t Point_t[25] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140
//                                               ,150,160,170,180,190,200,210};        
    
//Double_t Point_y[25] = {0.,0.340,0.910,1.610,2.270,2.750,3.190,3.630
//                        ,4.090,4.540,4.970,5.410,5.860,6.290,6.730,7.160,7.600,8.020,8.440
//                        ,8.870,9.290,9.730};  
                        
Double_t Point_t_old[46] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
                      105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180,
                      185, 190, 195, 200, 205, 210, 215, 220, 225};

Double_t Point_y_old[46] = {.120, .275, .550, .920, 1.305, 1.685, 2.060, 2.395, 2.665, 2.910, 3.145, 3.385,
                      3.625, 3.865, 4.110, 4.350, 4.600, 4.840, 5.065, 5.295, 5.535, 5.785, 6.030, 6.265, 6.495,
                      6.735, 6.975, 7.210, 7.430, 7.660, 7.885, 8.125, 8.355, 8.590, 8.820, 9.045, 9.275, 9.515,
                      9.740, 9.965, 10.195, 10.425, 10.650}; 
Double_t Point_yp[46];
Double_t Point_t[2000];
Double_t Point_y[2000];
Double_t Err_y[2000];

const Int_t max = 20000;         // max array size

struct data {

    Int_t n_max;            // array size 
    Int_t size = 20000;
    
    Double_t rtime[max];  // rawtime values
    Double_t time[max];   // time values 
    Double_t wire[max];   // Wire values
    Double_t n_wire[max];
    Double_t nhit;          // number of wires hit
	Double_t dis[max];           // Distance from wire 
	Double_t event;         // Events
	Double_t time_n[max];


}   treedata;

Int_t  newTTD()
{

    //Point_yp[0] = 0.;
    
    //for(Int_t i = 0; i < 44; i++){
     
      //  Point_yp[i] = Point_y[i+1] - .165;
    
        //cout << Point_yp[i];
    //}    
	
	
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
    }
	
    Double_t dis[20000];           // Distance from wire 
	Int_t size = 20000;
	
    TH1D *rawt = new TH1D("rawt","Rawtime Left Arm U1 Plane",300, 0,300);
    TH1D *Dist = new TH1D("Dist","Converted Drift Distance",100, -1, 16);
    TH2D *Clus = new TH2D("Clus","",300,0,30,20,-5,10);
    
    TSpline3 *sp = new TSpline3("Cubic Spline", Point_t, Point_y, n, "b2e2", 0, 0);
    
    //cout << sp -> Eval(1000) << endl; 
    

    
    TFile *f = new TFile("$APEX/rootfiles/apex_2060.root");
    TTree *T = (TTree*) f->Get("T");
    
    // Extract the T Tree 
    //TChain *T = new TChain("T");
    //T->Add("rootfiles/apex_2060.root");

    // Stops reading other branches
    T->SetBranchStatus("*", false);
    
    // Read array size for time array 
    T->SetBranchStatus("L.vdc.u1.rawtime", true );
    T->SetBranchAddress("L.vdc.u1.rawtime", treedata.rtime);
    
    // Read array size for time array 
    T->SetBranchStatus("Ndata.L.vdc.u1.time", true );
    T->SetBranchAddress("Ndata.L.vdc.u1.time", &treedata.n_max);
            
    // Read Times branch 
    T->SetBranchStatus("L.vdc.u1.time", true );
    T->SetBranchAddress("L.vdc.u1.time", treedata.time);
        
    // Read Number of Hits branch
    T->SetBranchStatus("L.vdc.u1.nhit", true );
    T->SetBranchAddress("L.vdc.u1.nhit", &treedata.nhit );
    
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.wire", true );
    T->SetBranchAddress("L.vdc.u1.wire", treedata.wire );
    
    
    
    TFile *ff = new TFile("D_apex_2060.root", "RECREATE");
    
    TTree *nT = new TTree("newT","");
    
    newT -> Branch("Size", &treedata.size, "size/I");
    newT -> Branch("Dist", treedata.dis,"dis[size]/D");
    newT -> Branch("Wire", treedata.wire, "wire[size]/D");
    newT -> Branch("Time", treedata.time_n, "time[size]/D");  
    newT -> Branch("Event", &treedata.event, "Event/D");  
   
    // Total ammount of Events
    Int_t nEntries = T -> GetEntries();    
    
    treedata.event = 0;         
    
    for(int i = 1; i < nEntries; i++)
    {
        if((i%10000) == 1) cout << "Event: " <<  i - 1 << endl;
    
        // Get the entry values 	                    
        T -> GetEntry(i);
           
        // Ensure that the events have hits
        if( treedata.nhit )
        {   
        
            treedata.size = treedata.nhit;
                    
            // Cycle through the hits in the event
            for(int j = 0; j < treedata.nhit ; j++)
            {
                
                if (treedata.time[j] > 0 && treedata.time[j]*pow(10,9) < 300){
                    //cout << wire[j] << endl;
                    // Fill in Histogram  [ns]
                    rawt -> Fill( treedata.time[j]*pow(10,9) );
                    Dist -> Fill(sp -> Eval( treedata.time[j]*pow(10,9)) );
                    //Clus -> Fill( treedata.wire[j], sp -> Eval( treedata.time[j]*pow(10,9)) );
                
                
                
                    treedata.time_n[j] = treedata.time[j]*pow(10,9);
        
                    //treedata.time_n[j] = treedata.rtime[j];
                
                    treedata.dis[j] = sp -> Eval( treedata.time_n[j] );
                    //cout << treedata.dis[j] << endl;
                }
                
            }
            
            
            newT -> Fill();
            
            treedata.event ++;         
            
        }	      
        else();
        
    }
    
    newT->Print();
    newT->Write();
    rawt->Write();
    Dist->Write();
    
   // delete f;
    
    TCanvas *Can = new TCanvas;
    Can -> Divide(2,2);
    Can -> cd(1);
    // Check on the Time distribution  
    //rawt -> Draw(); 
    T -> Draw("L.vdc.u1.time");
    
    Can -> cd(2);
    
    TGraph *g1 = new TGraph(n ,Point_t ,Point_y );
    
    g1->SetMarkerStyle(1);
    g1->SetMarkerSize(1.5);
   
    g1->GetXaxis()->SetLimits(0, 400);        // set real range
    g1->GetXaxis()->SetRangeUser(0, 400); // set visible range
    g1->GetXaxis()->SetTitle("Drift RawTime [ns]");
    g1->GetXaxis()->CenterTitle();

    g1->GetYaxis()->SetLimits(0, 30);
    g1->GetYaxis()->SetRangeUser(0, 30);
    g1->GetYaxis()->SetTitle("Drift Distance [mm]");
    g1->GetYaxis()->CenterTitle();
    
    g1->Draw("ap");
    
    sp -> Draw("lsame");
    
    Can -> cd(3);
    
    //Dist -> Draw();
    nT -> Draw("time");
    
    Can -> cd(4);
    
    Dist -> Draw();
    
}
    
           
