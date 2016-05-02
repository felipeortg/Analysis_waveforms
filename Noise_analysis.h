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

/* global variable declaration */

//Define name of the folder with the data measured
const TString data_folder = "ov_scan_pde_H2014";


//Define number of files on each voltage
const Int_t data_size = 2000;


//Define name of subfolders of each voltage
const Int_t vol_size = 6;
TString vol_folders[vol_size] = {"56.5V","57V","57.5V","58V","58.5V","59V"};




const double ns = 1e-9;

/* Auxiliar functions */

Double_t Amplitude_calc(int i){
    
    Double_t pe_volt;
    
    Char_t canvas_title[200];
    sprintf(canvas_title,"Amplitude calculation %s",vol_folders[i].Data());
    
    TH1D* volt_ampl= new TH1D(canvas_title, canvas_title, 150, -0.05, 0.5);
    
    cout<<"****----->Voltage amplitude analysis:"<< vol_folders[i] << endl;
    
    //loop over every measurement on a folder
    for (int j=0; j<data_size; j++) {
        if (j%500==0) {
            cout<<"Event: "<<j<<endl;
        }
        
        Char_t datafilename[100];
        
        sprintf(datafilename,"%s/%s/C1H%05i.csv",data_folder.Data(),vol_folders[i].Data(),j);
        
        TGraph* waveform = new TGraph(datafilename,"%lg %lg","/t;,");
        
        int ROWS_DATA = waveform->GetN();
        Double_t *time = waveform->GetX();
        Double_t *volts = waveform->GetY();
        Double_t  maxV=0.0;
        
        for (int row = 0; row<ROWS_DATA; row++) {
            if (time[row] > 0 && time[row]< 1 * ns) {
                if (maxV<volts[row]) {
                    maxV = volts[row];
                }
            }
        }
        
        volt_ampl->Fill(maxV);
        
        
    }
    
    
    TCanvas* c1= new TCanvas(canvas_title,canvas_title,100,100,900,700);
    
    volt_ampl->Draw();
    
    TF1 *f1 = new TF1("f1","gaus",0,0.8); //Change range for fit of MPV
    volt_ampl->Fit("f1","R");
    pe_volt= f1->GetParameter(1);
    
    c1->SetGrid();
    
    //Uncoment to save plot and check fit
    /*
    sprintf(canvas_title,"Plots/Amplitude_calculation_%s.pdf",vol_folders[i].Data());
    c1->Print(canvas_title,"pdf");
    */
     
    return pe_volt;
}



Int_t single_plot(const char * Voltage, const Int_t event){
    
    
    Char_t datafilename[100];
    
    sprintf(datafilename,"%s/%s/C1H%05i.csv",data_folder.Data(),Voltage,event);
    
    TGraph* waveform = new TGraph(datafilename,"%lg %lg","/t;,");
    
    TCanvas* single = new TCanvas("Single signal","Single signal",100,100,900,700);
    single->cd();
    waveform->SetTitle(datafilename);
    waveform->GetYaxis()->SetRangeUser(-0.1,0.5);
    waveform->GetXaxis()->SetRangeUser(-10*ns,80*ns);
    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
    waveform->GetYaxis()->SetTitleOffset(1.3);
    waveform->GetXaxis()->SetTitle("Time [s]");
    waveform->Draw("AL");
    single->SetGrid();
    
    sprintf(datafilename,"C1H%05i.pdf",event);
    single->Print(datafilename,"pdf");
    
    return 0;
    
    
}


#endif /* Noise_analysis_h */
