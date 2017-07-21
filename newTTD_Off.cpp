gROOT -> Reset();

//////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//
//                  Histrogram of Distance Offsets
//                         Dylan Linthorne
//
//  This program uses the root toolkit to extract APEX tree data and 
//  applies cuts on the clusters for a single VDC plane. The clusters 
//  are fit with lines to find their pivots and a offset histogram 
//  is created. The line selection process is governed by weighting of 
//  the Chi2 of the tracks. 
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

// Set max array size /arbitrary 
const Int_t max = 20000;

Bool_t flag = kFALSE;
Bool_t Fitflag = kFALSE;

TH1D *h_sep = new TH1D("h_sep", "Separation", 300, -10, 10 );
TH2D *h2_sep = new TH2D("h2_sep", "2D Separation", 400, -10, 10 , 400, 0, 400);

TH1D *h_off = new TH1D("h_sep", "Offset", 300, -10, 10 );
TH2D *h2_off = new TH2D("h2_sep", "Offset", 400, -10, 10 , 400, 0, 400);

struct data{

    Double_t rtime[max];	// Rawtime (TDC) branch [channel]
    Double_t dist[max]; 	// Drfit Distance Branch [mm]
    Double_t wire[max];		// Set of wires Branch [ spacing = 4.2 mm]
    Double_t nhit;		// Hits ber Event (true array size)
    Double_t event;		// Events branch 
    Double_t nclust;		// Number of clusters per event
    Double_t pivot;		// pivot point distance [mm]
    Double_t clsiz;		// Size of each cluster 


}   dat;

Double_t LineFit ();
Double_t fint();

TF1 *f1, *f2, *fbase;   // fit functions 

Double_t Chi2 ( Double_t x_track[], Double_t y_track[], const Int_t n){

    // create graph of line to apply fit
    TGraph *line = new TGraph(n, x_track, y_track);    
    
    // Linear fit function, open parametres 
    TF1 *func = new TF1("fit", "[0] + [1]*x", x_track[0], x_track[n]);

    // fit line to function
    line -> Fit(func);
    
    // return Chi2
    return func -> GetChisquare();
}
 
Double_t finter(double *x, double*par) {
   
    // minimum function for two function intersection
    return TMath::Abs(f1->EvalPar(x,par) - f2->EvalPar(x,par));

}

Double_t xinter(double *x, double *par){
    
    // minimum function for x-intercept calculation
    return TMath::Abs(f1->EvalPar(x,par) - fbase->EvalPar(x,par)); 

}

Double_t separ(double x_0, double x_p, double y_p){

    // distance formula between two points
    Double_t sqaure =  pow(x_p - x_0, 2) +  pow(y_p, 2);
    
    return  sqrt( sqaure );

}

Int_t newTTD_Off() 
{

    Double_t Tru_dis[max];	// Distances after Cuts 
    Double_t Tru_wir[max];	// Wires after cuts
    
    // The separation histogram
    //TH1D *h_sep = new TH1D("h_sep", "Separation", 300, -5, 5 );

    // Opening file, extracting Tree leaf :FIX ME: make this 
    // More universal for my data
    TFile *f = new TFile("$APEX/rootfiles/apex_2060.root");
    TTree *T = (TTree*) f->Get("T");
    
    // Stops reading other branches
    T->SetBranchStatus("*", false);
    
    // Read array size for time array 
    T->SetBranchStatus("L.vdc.u1.rawtime", true );
    T->SetBranchAddress("L.vdc.u1.rawtime", dat.rtime);
    
    // Read Times branch 
    T->SetBranchStatus("L.vdc.u1.dist", true );
    T->SetBranchAddress("L.vdc.u1.dist", dat.dist);
        
    // Read Number of Hits branch
    T->SetBranchStatus("L.vdc.u1.nhit", true );
    T->SetBranchAddress("L.vdc.u1.nhit", &dat.nhit );
    
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.wire", true );
    T->SetBranchAddress("L.vdc.u1.wire", dat.wire );
    
    T->SetBranchStatus("g.evnum", true );
    T->SetBranchAddress("g.evnum", &dat.event );
     
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.nclust", true );
    T->SetBranchAddress("L.vdc.u1.nclust", &dat.nclust );
    
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.clpos", true );
    T->SetBranchAddress("L.vdc.u1.clpos", &dat.pivot );
    
    // Read wires branch
    T->SetBranchStatus("L.vdc.u1.clsiz", true );
    T->SetBranchAddress("L.vdc.u1.clsiz", &dat.clsiz );
    
    const Double_t min_dstr = 2.0;	// Minimum distance for start of cluster
    const Double_t max_dstr = 2000.0;	// Maximum distance for cluster :FIX ME:
    
    Bool_t badWir = kFALSE;
    Int_t accept;
    Double_t diff;
    
    // Extract ammount of events in Tree file
    Int_t nEntries = T -> GetEntries();    

    for(int i = 0; i < nEntries ; i++)
    {
        if((i%10000) == 1) cout << "Event: " <<  i - 1 << endl;
    
        // Get the entry values 	                    
        T -> GetEntry(i);
           
        // Ensure that the events have hits
        if( dat.nhit ){  
        
	    // Cut for single cluster events with atleast 4 hits 
            if( dat.nclust == 1 && dat.clsiz > 3 ){
            
            // intialize accepted cut counter
            accept = 0;   
 
                for(Int_t j = 0; j < dat.nhit; j++)
                {       
                
                        if( dat.wire[j] < 96 || dat.wire[j] > 111) continue;
                 
                        // Looking at first entry of cluster 
                        if( j == 0)
                        {
                        
		                    // Cut on the angle, ensure proper cluster hit
                            if(dat.dist[j] < min_dstr) continue;
                            
                            diff = dat.wire[j];
                            
                        
                        }
            		
			            // If any extreme hits occur, exit iteration 
                        if( dat.dist[j] > max_dstr) continue;  
                 
		                // Feed cut data to the new arrays 
                        Tru_dis[accept] = dat.dist[j];
                        Tru_wir[accept] = dat.wire[j];
                        
                        accept++;
                        
                        // if the wire spacing is more than '1' than set flag
                        if( diff - dat.wire[j] > 1) badWir = kTRUE; 
             
                        
                    
                }
                
                if( badWir == kTRUE) continue;

                
                Double_t diff = 4.26;      // [mm]
                Double_t x_dist[max];
                Double_t sum = Tru_wir[0];
                Double_t avg_w = 0;
                x_dist[0] = 0;  
   
   
                for( Int_t h = 1; h < accept; h++){
                    
                    sum = sum + Tru_wir[h];
                
                    if(Tru_wir[h] - Tru_wir[h -1] == 1 ) x_dist[h] = diff + x_dist[h - 1];
       
                    if(Tru_wir[h] - Tru_wir[h -1] == 2 ) x_dist[h] = 2.*diff + x_dist[h -1];
       
                }
                
                if( accept != 0) avg_w = sum/accept;
                
		        // Feed to linear fit algorithm 
                // h_sep -> Fill( LineFit(Tru_dis, x_dist, dat.pivot, accept ));  
                 
                LineFit(Tru_dis, x_dist, dat.pivot, accept, avg_w);
            }
	    }
    } 
    
    TCanvas *Can = new TCanvas();
    
    Can -> Divide(2,2);
    
    Can -> cd(1);
    h_sep -> Draw();
    
    Can -> cd(2);
    h2_sep -> Draw();
    
    Can -> cd(3);
    h_off -> Draw();
    
    Can -> cd(4);
    h2_off -> Draw();
    
    return 0;
    
}

Double_t LineFit ( Double_t y_dist[], Double_t x_dist[], Double_t pivot, const Int_t hits, Double_t avg){

    Double_t line_1[10];        // left track with ambigous point
    Double_t line_2[10];        // right track with ambigous point
    Double_t Lline_1[10];       // left track without ambigous point
    Double_t Lline_2[10];       // right track without ambigous point
    Double_t x_1[10];           // x/wire coord for left/amb
    Double_t x_2[10];           // x/wire coord for right/amb
    Double_t Lx_1[10];          // x/wire coord for left/unamb
    Double_t Lx_2[10];          // x/wire coord for right/unamb

    Int_t count_1 = 0;
    Int_t count_2 = 0;

    // start sorting look to isolate these tracks
    for (Int_t j = 0; j < hits; j++){
    
        if(j + 1 == hits){
    
            line_2[count_2] = y_dist[j];
            x_2[count_2] = x_dist[j];
            count_2++;
        }
        
        else if( y_dist[j] > y_dist[j+1]){
        
            line_1[count_1] = y_dist[j];
            x_1[count_1] = x_dist[j];
            
            count_1++;
            
            if(y_dist[j + 1] < y_dist[j + 2]){
            
                line_1[count_1] = y_dist[j + 1];
                x_1[count_1] = x_dist[j + 1];
                count_1++;
            }
                        
        }
        
        else if(y_dist[j] < y_dist[j +1]){
        
            line_2[count_2] = y_dist[j];
            x_2[count_2] = x_dist[j];
            
            count_2++;
        }
        
    }
    
    // create counts for Lline lines which are reduced to n-1
    Double_t  Lcount_1 = count_1 - 1;
    Double_t  Lcount_2 = count_2 - 1;
    
    // loop through the proper amount of counts
    for(Int_t i = 0; i < Lcount_1 ; i++){
    
        Lline_1[i] = line_1[i];
        Lx_1[i] = x_1[i];
    }
    
    for(Int_t k = 0; k < Lcount_2 ; k++){
    
        Lline_2[k] = line_2[k + 1];
        Lx_2[k] = x_2[k + 1];
    }
    
    // find the Chi2 for each line
    Double_t chi2_1 = Chi2(x_1, line_1, count_1);
    Double_t chi2_2 = Chi2(x_2, line_2, count_2);
    Double_t Lchi2_1 = Chi2(Lx_1, Lline_1, Lcount_1);
    Double_t Lchi2_2 = Chi2(Lx_2, Lline_2, Lcount_2);
    
    // match the chi2 to find the true set of points for each track
    if (TMath :: Abs(chi2_1 - Lchi2_1) < TMath :: Abs(chi2_2 - Lchi2_2)){
    
        // send the tracks to the fit algorithm 
        return fint(x_1, line_1, count_1, Lx_2, Lline_2, Lcount_2, avg);   
    
    }
    
    if (TMath :: Abs(chi2_1 - Lchi2_1) > TMath :: Abs(chi2_2 - Lchi2_2)){
    
        // send the tracks to the fit algorithm 
        return fint(Lx_1, Lline_1, Lcount_1, x_2, line_2, count_2, avg);  
        
    
    } 

    
}

Double_t fint(double x_1[], double y_1[], int n_1, double x_2[], double y_2[], int n_2, double avg){
   
    // Create graphs for both lower and higher tracks
    TGraph *L_line = new TGraph(n_1, x_1, y_1);    
    TGraph *H_line = new TGraph(n_2, x_2, y_2);    
    
    // create fitting polynomials for tracks and base y = 0   
    f1 = new TF1("f1","[0] + [1]*x",0 ,10);
    f2 = new TF1("f2","[0] + [1]*x",0 ,10);
    fbase = new TF1("fbase", "0", 0, 10);
   
    // fit the lines to the proper parameters 
    L_line -> Fit(f1, "QN", "", x_1[0],x_1[n_1]);
    H_line -> Fit(f2, "QN", "", x_2[0],x_2[n_2]);
   
    TFitResultPtr r = L_line -> Fit( f1 ,"S");
    Fitflag = r -> IsEmpty();
    //f1->Draw();
    //f2->Draw("same");
    //fbase -> Draw("same");
    
    // minimum functions for the intercepts
    TF1 *fint = new TF1("fint",finter,0,10,0);
    TF1 *x_base = new TF1("x_base",xinter,0,10,0);
   
    double xint = fint->GetMinimumX();
    double x_int = fint -> GetMinimumX();
    double y_int = f1 -> Eval(x_int);
    
    double base  = x_base -> GetMinimumX();
    double delta_s = separ(base, x_int, y_int);
    
    // set markers for the intercepts
    TMarker *b = new TMarker(base, 0 ,24);
    TMarker *m = new TMarker(xint,f1->Eval(xint),24);
    b->SetMarkerColor(kGreen);
    m->SetMarkerColor(kRed);
    b->SetMarkerSize(3);
    m->SetMarkerSize(3); 
    
    //b->Draw();
    //m->Draw();
    
    //printf("x-intercept = %g\n",x_int);
    //printf("y-intercept = %g\n",y_int);
    //printf("x-base = %g\n",base);
    //printf("delta-s = %g\n", delta_s);
        
    
    if( y_int < 0) delta_s = (-1)*delta_s;
    
    if(Fitflag == kTRUE) 
    {
        delta_s = 1e38;
        y_int = 1e38;
    }
    
    h_sep -> Fill(delta_s);
    h2_sep -> Fill(delta_s, avg);
    
    h_off -> Fill(y_int);
    h2_off -> Fill(y_int, avg);
    
    return delta_s;
    
}

