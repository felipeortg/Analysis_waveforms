//
//  Noise_analysis.h
//  
//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//
//

#ifndef Noise_analysis_h
#define Noise_analysis_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>


#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TVirtualPad.h>
#include <TMath.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

using namespace std;

/* global variable declaration */

const double ns = 1e-9;

/* Auxiliar functions */

struct globalArgs_t
{
    const char* data_folder;                /* -d option */
    const char* arg_pathToSetupFile;             /* -S option */
    int save_all;               /* -a  */
    //int  arg_eventplot;                        /* -E option */
    
} globalArgs;



Double_t Amplitude_calc(const char* vol_folder, Int_t data_size){
    
    Double_t pe_volt;
    
    Char_t canvas_title[200];
    sprintf(canvas_title,"Amplitude calculation %s",vol_folder);
    
    TH1D* volt_ampl= new TH1D(canvas_title, canvas_title, 150, -0.05, 0.5);
    
    cout<<"****----->Amplitude calculation of pe: "<< vol_folder << endl;
    
    //loop over every measurement on a folder
    for (int j=0; j<data_size; j++) {
        if (j%500==0) {
            cout<<"Event: "<<j<<endl;
        }
        
        Char_t datafilename[100];
        
        sprintf(datafilename,"%s/%s/C1H%05i.csv",globalArgs.data_folder,vol_folder,j);
        
        TGraph* waveform = new TGraph(datafilename,"%lg %lg","/t;,");
        
        int ROWS_DATA = waveform->GetN();
        Double_t *time = waveform->GetX();
        Double_t *volts = waveform->GetY();
        Double_t  maxV=0.0;
        
        //loop over the time window to get height of first pulse
        for (int row = 0; row<ROWS_DATA; row++) {
            if (time[row] > 0 && time[row]< 1 * ns) { //1 nanosecond window for first pulse
                if (maxV<volts[row]) {
                    maxV = volts[row];
                }
            }
        }
        
        volt_ampl->Fill(maxV);
        
        
    }
    
    
    TF1 *f1 = new TF1("f1","gaus",0,0.8); //Change range for fit of MPV
    volt_ampl->Fit("f1","R");
    pe_volt= f1->GetParameter(1);
    
    if (globalArgs.save_all==1) {
        volt_ampl->Write();
        //To save canvas instead of histogram
        /*TCanvas* c1= new TCanvas(canvas_title,canvas_title,100,100,900,700);
        volt_ampl->Draw();
        c1->SetGrid();
        c1->Write();*/
    }

    
    return pe_volt;
}



//Int_t single_plot(const char * Voltage, const Int_t event){
//
//
//    Char_t datafilename[100];
//
//    sprintf(datafilename,"%s/%s/C1H%05i.csv",data_folder.Data(),Voltage,event);
//
//    TGraph* waveform = new TGraph(datafilename,"%lg %lg","/t;,");
//
//    TCanvas* single = new TCanvas("Single signal","Single signal",100,100,900,700);
//    single->cd();
//    waveform->SetTitle(datafilename);
//    waveform->GetYaxis()->SetRangeUser(-0.1,0.5);
//    waveform->GetXaxis()->SetRangeUser(-10*ns,80*ns);
//    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
//    waveform->GetYaxis()->SetTitleOffset(1.3);
//    waveform->GetXaxis()->SetTitle("Time [s]");
//    waveform->Draw("AL");
//    single->SetGrid();
//
//    sprintf(datafilename,"C1H%05i.pdf",event);
//    single->Print(datafilename,"pdf");
//
//    return 0;
//
//
//}






#endif /* Noise_analysis_h */
