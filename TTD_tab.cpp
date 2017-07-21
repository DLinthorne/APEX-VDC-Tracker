
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

void Print(TH1D*, TH1D*);

Double_t accedent(TH1D*);

ofstream fout("ttd.dat");
ofstream ffout("ttd_v.dat");

Int_t TTD_tab(){    


    TFile *f = new TFile("hist.root");
    TH1D *time = (TH1D*) f -> Get("rawt");
    TH1D *clus = (TH1D*) f -> Get("clus");

    TH1D *res= new TH1D("res","Residual",10000,-.01,.01);
    
    Double_t alpha; 
    Double_t del[20000];
    Double_t vel[20000];
    Double_t dtime[20000];    
    Double_t err_v[20000];
    Double_t err_d[20000];
    
    Double_t acc;
    Double_t er_v;
    
    alpha = 4.7/(clus -> Integral());     // 20 [mm]
    
    //vel[0] = (alpha/5)*( time -> Integral(600,605) - 5*accedent(time) );  
    
    Int_t bin = 1120;
    Int_t incr = 5;
    
    acc = accedent(time);
    
    vel[0] = 0; 
    dtime[0] = 0;
    del[0] = 0;
    err_v[0] = 0;
    err_d[0] = 0;
    
    
    fout << dtime[0] << " " << del[0] << " " << err_d[0] << endl;
    
    for(Int_t i = 1; i < 400; i++){
    
    
        vel[i] = 2*(alpha/incr)*( time -> IntegralAndError(bin - incr, bin ,er_v) - incr*acc );  
        
        del[i] = vel[i]*(.5*incr) + del[i-1]; 
        
        err_v[i] = (alpha/incr)*er_v;
        err_d[i] = err_v[i]*(incr/2);
        
        //res -> Fill(err_d[i]);
        
        bin = bin - incr;
        
        dtime[i] = dtime[i - 1] + .5*incr;
        
        //dtime[i] = dtime[i -1] + i ;
        
        fout << dtime[i] << " " << del[i] << " " << err_d[i] <<  endl;
        ffout << dtime[i] << " " << vel[i] << " " << err_v[i] <<  endl;
    
    }
    
    TGraphErrors *g = new TGraphErrors(400, dtime, vel, 0, err_v);
    TGraphErrors *h = new TGraphErrors(400, dtime, del, 0, err_d);
    
    g->SetMarkerStyle(1);
    g->SetMarkerSize(1.5);
   
    g->GetXaxis()->SetLimits(0, 400);        // set real range
    g->GetXaxis()->SetRangeUser(0, 400); // set visible range
    g->GetXaxis()->SetTitle("Rawtime [TDC Channel = 0.5 ns]");
    g->GetXaxis()->CenterTitle();

    g->GetYaxis()->SetLimits(0, .15);
    g->GetYaxis()->SetRangeUser(0, .15);
    g->GetYaxis()->SetTitle("Drift Velocity ");
    g->GetYaxis()->CenterTitle();
    
    h->SetMarkerStyle(1);
    h->SetMarkerSize(1.5);
   
    h->GetXaxis()->SetLimits(0, 400);        // set real range
    h->GetXaxis()->SetRangeUser(0, 400); // set visible range
    h->GetXaxis()->SetTitle("Rawtime [TDC Channel = 0.5 ns]");
    h->GetXaxis()->CenterTitle();

    h->GetYaxis()->SetLimits(0, 18);
    h->GetYaxis()->SetRangeUser(0, 18);
    h->GetYaxis()->SetTitle("Drift Distance [mm] ");
    h->GetYaxis()->CenterTitle();
    
    
    
    
    
    TCanvas *C = new TCanvas();
    
    C -> Divide(2,2);
    C ->cd(1);
    
    g -> Draw();
    
    C -> cd(2);
    
    h -> Draw();
    
    C -> cd(3);
    
    rawt -> Draw();
    
    return vel[40];
    //return dis[0]*pow(10,3); 
    
    

}

Double_t accedent( TH1D* time){
    
    Double_t seg_a, seg_b, avg;
    Double_t err_a, err_b;

    seg_a = time -> IntegralAndError(0, 500, err_a);
    
    seg_b = time -> IntegralAndError(1500,3000, err_b);

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
