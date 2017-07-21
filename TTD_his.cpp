gROOT -> Reset();

//////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//
//                  Histrogram of Distance Offsets
//                         Dylan Linthorne
//
//  This program uses the root toolkit to extract tree data from the
//  APEX test run and create a TTD (Drift time to Drift distance conversion)
//
//
//////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#include <iostream>
#include <string>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"


using namespace std;

// prototyping 
void Print(TH1D*, TH1D*);
Double_t accedent(TH1D*);
Double_t TTD(TH1D*, TH1D*);

// Max array size 
const Int_t n_max = 20000;

struct var {

    Double_t rtime[n_max];      // Rawtime array from tree
    Double_t time[n_max];       // converted time array from tree
    Double_t wire[n_max];       // wire array from tree
    
    Double_t nhit;              // Number of hits per event
    Double_t nclust;            // Number of clusters per event 
    
    Double_t In_err_a;          // errors for left accedental in rawtime
    Double_t In_err_b;          // errors for right accedentals in rawtime

}   fvar;

ofstream fout("ttd.dat");       // output file for distances
ofstream ffout("ttd_v.dat");    // output file for drift velocities 

Int_t TTD_his(){
 
    
    TH1D *rawt = new TH1D("rawt","Rawtime Left Arm U1 Plane",3000, 0,3000);
    TH1D *clus = new TH1D("clus","Numbe of clusters",100,0,100);
    
    TFile *f = new TFile("$APEX/rootfiles/apex_2060.root");
    TTree *T = (TTree*) f->Get("T");

    // Read array size for time array 
    T->SetBranchStatus("L.vdc.u1.rawtime", true );
    T->SetBranchAddress("L.vdc.u1.rawtime", fvar.rtime);
            
    // Read Times branch 
    T->SetBranchStatus("L.vdc.u1.time", true );
    T->SetBranchAddress("L.vdc.u1.time", fvar.time);
        
    // Read Number of Hits branch
    T->SetBranchStatus("L.vdc.u1.nhit", true );
    T->SetBranchAddress("L.vdc.u1.nhit", &fvar.nhit );
    
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.wire", true );
    T->SetBranchAddress("L.vdc.u1.wire", fvar.wire );
    
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.nclust", true );
    T->SetBranchAddress("L.vdc.u1.nclust", fvar.nclust );
    
    Int_t nEntries = T -> GetEntries();             
    
    Bool_t Hit;
    
    for(int i = 1; i < nEntries; i++)
    {
        if((i%10000) == 1) cout << "Event: " <<  i - 1 << endl;
    
        // Get the entry values 	                    
        T -> GetEntry(i);
        
        Hit = kFALSE;
           
        // Ensure that the events have hits
        if( fvar.nhit )
        {  
                    
            // Cycle through the hits in the event
            for(int j = 0; j < fvar.nhit ; j++)
            {
                   // cut for the set of chosen wires
                   
                   
                   if(fvar.wire[j] < 111 && fvar.wire [j] > 96){
                        
                        
                        rawt -> Fill(fvar.rtime[j]);
                        
                        Hit = kTRUE;         
                        
                   }
            }
            
            if(Hit == kTRUE)   clus -> Fill(fvar.nclust);                  
            
            
        }	      
        else();
        
    }
    
   //Print(rawt,clus);
   cout << TTD(rawt, clus) << endl;
   
   TFile *ff = new TFile("hist.root","RECREATE");
   
   rawt -> Write();
   clus -> Write();
   
   
}

Double_t TTD ( TH1D* time, TH1D* clus ){

    
    TH1D *res= new TH1D("res","Residual",10000,-.01,.01);
    
    Double_t alpha; 
    Double_t del[20000];
    Double_t vel[20000];
    Double_t dtime[20000];    
    Double_t err_v[20000];
    Double_t err_d[20000];
    
    Double_t acc;
    Double_t er_v;
    
    alpha = 8/(clus -> Integral());     // 20 [mm]
    
    //vel[0] = (alpha/5)*( time -> Integral(600,605) - 5*accedent(time) );  
    
    Int_t bin = 1200;
    Int_t incr = 5;
    
    acc = accedent(time);
    
    vel[0] = 0; 
    dtime[0] = 1200;
    del[0] = 0;
    err_v[0];
    err_d[0];
    
    
    fout << dtime[0] << " " << del[0] << endl;
    
    for(Int_t i = 1; i < 600; i++){
    
    
        vel[i] = 2*(alpha/incr)*( time -> IntegralAndError(bin - incr, bin,er_v) - incr*acc );  
        
        del[i] = vel[i]*(incr/2) + del[i-1]; 
        
        err_v[i] = (alpha/incr)*er_v;
        err_d[i] = err_v[i]*(incr/2);
        
        res -> Fill(err_d[i]);
        
        bin = bin - incr;
        
        dtime[i] = dtime[i - 1] - incr;
        
        //dtime[i] = dtime[i -1] + incr ;
        
        fout << dtime[i] << " " << del[i] << " " << err_d[i] <<  endl;
        ffout << dtime[i] << " " << vel[i] << " " << err_v[i] <<  endl;
    
    }
    
    TGraphErrors *g = new TGraphErrors(300, dtime, vel, 0, err_v);
    TGraphErrors *h = new TGraphErrors(300, dtime, del, 0, err_d);
    
    g->SetMarkerStyle(1);
    g->SetMarkerSize(1.5);
   
    g->GetXaxis()->SetLimits(400, 1200);        // set real range
    g->GetXaxis()->SetRangeUser(400, 1200); // set visible range
    g->GetXaxis()->SetTitle("Rawtime [TDC Channel = 0.5 ns]");
    g->GetXaxis()->CenterTitle();

    g->GetYaxis()->SetLimits(0, .15);
    g->GetYaxis()->SetRangeUser(0, .15);
    g->GetYaxis()->SetTitle("Drift Velocity ");
    g->GetYaxis()->CenterTitle();
    
    
    
    TCanvas *C = new TCanvas();
    
    C -> Divide(2,2);
    C ->cd(1);
    
    g -> Draw();
    
    C -> cd(2);
    
    h -> Draw();
    
    C -> cd(3);
    
    res -> Draw();
    
    return vel[40];
    //return dis[0]*pow(10,3); 
    
    

}

Double_t accedent( TH1D* time){
    
    Double_t seg_a, seg_b, avg;

    seg_a = time -> IntegralAndError(0, 500, fvar.In_err_a);
    
    seg_b = time -> IntegralAndError(1500,3000, fvar.In_err_b);

    // Taking the average for a single bin
    avg = (seg_a + seg_b)/(2000);
    
    return avg;

}

void Print (TH1D* h1, TH1D* h2){

    TCanvas *Screen = new TCanvas();
    Screen -> Divide(2,2);
    
    Screen -> cd(1);
    h1 -> Draw();
      
    Screen -> cd(2);
    h2 -> Draw();

}
