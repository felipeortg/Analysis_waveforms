
//
//  Noise_analysis.C
//  
//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//
//

#include "Noise_analysis.h"



TFile* Noise_analysis(int singleplot=0 ,const char * Voltage="56.5V",int event=0) {
    
    Int_t Event;
    Char_t Category[15];
    TGraph* waveform = 0;
    Double_t Amp;
    Double_t V_meas;
    
    
    const double reject_time = 4;
    
    double pe = 0.07;
    int row = 0;
    int full_n_file = 0;
    int ap_n_file = 0;
    int xtalk_n_file = 0;
    int dxtalk_n_file = 0;
    int time_dist_n_file = 0;
    
    int direct_xtalk_pulse;
    Double_t direct_xtalk_pulse_cnt=0;
    int xtalk_pulse;
    Double_t xtalk_pulse_cnt = 0;
    int after_pulse;
    Double_t after_pulse_cnt=0;
    Double_t event_cnt = 0;
    double sig_max = 0;
    double time_of_max = 0;
    double sig_max_first = 0;
    double time_of_max_first = 0;
    int max_cnt = 0;
    int max_noise_cnt = 0;
    int max_found = 0;
    
    
    //Define Voltages measured
    
    
    //Define thresholds
    
    const double after_pulse_th = 0.38;
    const double direct_xtalk_th = 1.17;
    const double xtalk_th = 0.85;
    const double time_dist_th = 0.4;
    
    
    if (singleplot) {
        single_plot(Voltage,event);
    }
    
    
    //Create a root tree with the graph of the waveform of each event and
    //classify them
    TString filename = "noiseanalysis.root";
    
    TFile *hfile = 0;
    hfile = TFile::Open(filename,"RECREATE");
    
    
    TTree *tree = new TTree("T","Noise Analysis");
    tree->Branch("Event",&Event,"Event/I");
    tree->Branch("Category",Category,"Category/C");
    tree->Branch("waveform","TGraph",&waveform);
    tree->Branch("V_meas",&V_meas,"V_meas/D");
    
    
    
    Event=0;
    
    TGraph* Correl_noise[4];
    Correl_noise[0] = new TGraph();
    Correl_noise[1] = new TGraph();
    Correl_noise[2] = new TGraph();
    Correl_noise[3] = new TGraph();
    
    TCanvas* c1[vol_size];
    TCanvas* c2[vol_size];
    TCanvas* c3[vol_size];
    TCanvas* c4[vol_size];
    
    
    cout<<"****----->Voltage Breakdown calculation"<< endl;
    
    vector <Double_t> pe_volt;
    TGraph *Vbias_ver= new TGraph();
    
    for (int i=0; i<vol_size; i++) {
        pe_volt.push_back(Amplitude_calc(i));
        V_meas = vol_folders[i].Atof();
        Vbias_ver->SetPoint(i, pe_volt.at(i), V_meas);
    }
    
    TCanvas* ca= new TCanvas("Voltage Breakdown calculation","Voltage Breakdown calculation",100,100,900,700);
    Vbias_ver->SetTitle("Voltage Breakdown calculation");
    Vbias_ver->GetYaxis()->SetTitle("Bias Volatge [V]");
    Vbias_ver->GetYaxis()->SetTitleOffset(1.2);
    Vbias_ver->GetXaxis()->SetTitle("Mean peak amplitude [V]");
    Vbias_ver->Draw("AP*");
    ca->SetGrid();
    
    TFitResultPtr  fit = Vbias_ver->Fit("pol1","S");
    Double_t VBD= fit->Value(0);
 
    ca->Print("Plots/VBD.pdf","pdf");
    
    

    /////////////////
    // Loop over all Voltages measured
    /////////////////
    for (int i=0; i<vol_size; i++) {
        
        //Important to reinitialize, use value 1 to plot axis of TGraph()
        
        int color1 = 0;
        int color2 = 0;
        int color3 = 0;
        int color4 = 0;
        
        direct_xtalk_pulse_cnt = 0;
        xtalk_pulse_cnt = 0;
        after_pulse_cnt = 0;
        event_cnt = 0;
        
        cout<<"****----->Voltage analyzed:"<< vol_folders[i] << endl;
        
        
        //Define amplitude measured
        
        Double_t pe = pe_volt.at(i);
        
        Char_t canvas_title[40];
        sprintf(canvas_title,"Direct CrossTalk %s",vol_folders[i].Data());
        c1[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        sprintf(canvas_title,"Delayed CrossTalk %s",vol_folders[i].Data());
        c2[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        sprintf(canvas_title,"After Pulse %s",vol_folders[i].Data());
        c3[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        sprintf(canvas_title,"Clean %s",vol_folders[i].Data());
        c4[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        
        //loop over every measurement on a folder
        for (int j=0; j<data_size; j++) {
            
            Char_t datafilename[100];
            
            sprintf(datafilename,"%s/%s/C1H%05i.csv",data_folder.Data(),vol_folders[i].Data(),j);
            
            waveform = new TGraph(datafilename,"%lg %lg","/t;,");
            waveform->SetName(datafilename);
            
            int ROWS_DATA = waveform->GetN();
            Double_t *time = waveform->GetX();
            Double_t *volts = waveform->GetY();
            
            Amp = waveform->GetY()[0];
            
            V_meas = vol_folders[i].Atof()-VBD;
            
            
            
            
            /////////////////////////////////////////////////////
            // Data filtering into the different type of events
            // direct x-talk  AP   delayed x-talk
            /////////////////////////////////////////////////////
            after_pulse = 0;
            xtalk_pulse = 0;
            direct_xtalk_pulse = 0;
            sig_max = 0;
            max_cnt = 0;
            max_found = 0;
            
            
            /// Print Single event
            
            bool print;
            print=false;
            
            
            // direct x-talk
            for (row = 0; row < ROWS_DATA; row++) {
                if ((time[row]>0 * ns)&(volts[row] > direct_xtalk_th * pe)) {// time larger 0ns
                    direct_xtalk_pulse++;
                }
            }
            //            if (direct_xtalk_pulse > 0) {
            //                direct_xtalk_pulse_cnt++;
            //
            //                if (print) {
            //                    cout<<"Direct x-talk cnt = %d\n"<< direct_xtalk_pulse_cnt<<endl;
            //                }
            //                //printf("Direct x-talk cnt = %d\n", direct_xtalk_pulse_cnt);
            //            }
            // after-pulse threshold
            for (row = 0; row < ROWS_DATA; row++) {
                if ((time[row]>reject_time*ns)&(volts[row] > after_pulse_th * pe)) {// time larger 4ns and ap_th
                    after_pulse++;
                }
            }
            //            if (after_pulse > 0) {
            //                after_pulse_cnt++;
            //                if (print) {
            //                    cout<<"AP and Delayed x-talk event cnt = %d\n"<<after_pulse_cnt<<endl;
            //                }
            //                //printf("AP and Delayed x-talk event cnt = %d\n", after_pulse_cnt);
            //            }
            // delayed x-talk
            for (row = 0; row < ROWS_DATA; row++) {
                if ((time[row]>reject_time*ns)&(volts[row] > xtalk_th * pe)) {// time larger 4ns and larger xtalk_th
                    xtalk_pulse++;
                }
            }
            //            if (xtalk_pulse > 0) {
            //                xtalk_pulse_cnt++;
            //                if (print) {
            //                    cout<<"X-talk event cnt = %d\n"<< xtalk_pulse_cnt<<endl;
            //                }
            //                //printf("X-talk event cnt = %d\n", xtalk_pulse_cnt);
            //            }
            //
            /////////////////////////////////////////////////////////////////////
            // Detect peaks in data after 4ns, count the number of maxima and
            // measure the time of arrival of first maxima
            /////////////////////////////////////////////////////////////////////
            max_noise_cnt = 0;
            for (row = 0; row < ROWS_DATA; row++) {
                if (time[row] > reject_time*ns) {// time larger 4ns
                    if (volts[row] > sig_max) {
                        sig_max = volts[row];        // set the max
                        time_of_max = time[row];    // time max
                        max_noise_cnt++;                   // set the histeresis cnt
                    }
                    else if (max_noise_cnt > 0) max_noise_cnt--;  //  count down if no new max is reached
                    // decide if real max or only noise, threshold has to be reached in case of a real max
                    if (max_noise_cnt>2 && sig_max > time_dist_th * pe) {
                        max_cnt++;
                        if (max_cnt == 1) {
                            sig_max_first = sig_max;  // sig max
                            time_of_max_first = time_of_max; // time max
                            max_found = 1;
                            //printf("First max found: sig=%f  time=%f ns max_noise_cnt=%d\n", sig_max, time_of_max / ns, max_noise_cnt);
                        }
                        
                        //printf("Max number is: %d   cnt=%d\n", max_cnt, max_noise_cnt);
                    }
                } // 3ns
            } //loop over data samples
            
            bool clean = true;
            char graph_title[50];
            
            if (direct_xtalk_pulse > 0){
                direct_xtalk_pulse_cnt++;
                sprintf(Category,"ImmCrosstalk");
                c1[i]->cd();
                
                color1++;
                if (color1==10 || color1==0) {
                    color1++;
                } else if (color1>49) color1=2;
                waveform->SetLineColor(color1);
                
                if (color1>1) {
                    waveform->Draw("SAME");
                }else{
                    sprintf(graph_title,"Direct CrossTalk %s",vol_folders[i].Data());
                    waveform->SetTitle(graph_title);
                    waveform->GetYaxis()->SetRangeUser(-0.1,pe*(15+3*i)/10);
                    waveform->GetXaxis()->SetRangeUser(-10*ns,80*ns);
                    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
                    waveform->GetYaxis()->SetTitleOffset(1.3);
                    waveform->GetXaxis()->SetTitle("Time [s]");
                    waveform->Draw("AL");
                    c1[i]->SetGrid();
                }
                clean = false;
            }
            
            // only delayed x-talk
            if (xtalk_pulse > 0 && direct_xtalk_pulse == 0){
                xtalk_pulse_cnt++;
                sprintf(Category,"DelCrosstalk");
                c2[i]->cd();
                
                color2++;
                if (color2==10 || color2==0) {
                    color2++;
                } else if (color2>49) color2=2;
                waveform->SetLineColor(color2);
                
                if (color2>1) {
                    waveform->Draw("SAME");
                }else{
                    sprintf(graph_title,"Delayed cross-talk %s",vol_folders[i].Data());
                    waveform->SetTitle(graph_title);
                    waveform->GetYaxis()->SetRangeUser(-0.1,pe*(15+3*i)/10);
                    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
                    waveform->GetYaxis()->SetTitleOffset(1.3);
                    waveform->GetXaxis()->SetTitle("Time [s]");
                    waveform->Draw("AL");
                    c2[i]->SetGrid();
                }
                clean = false;
            }
            
            //  after pulse
            if (after_pulse > 0 && xtalk_pulse == 0 && direct_xtalk_pulse == 0){
                after_pulse_cnt++;
                sprintf(Category,"AfterPulse");
                c3[i]->cd();
                
                color3++;
                if (color3==10 || color3==0) {
                    color3++;
                } else if (color3>49) color3=2;
                waveform->SetLineColor(color3);
                
                if (color3>1) {
                    waveform->Draw("SAME");
                }else{
                    sprintf(graph_title,"After pulse %s",vol_folders[i].Data());
                    waveform->SetTitle(graph_title);
                    waveform->GetYaxis()->SetRangeUser(-0.1,pe*(15+3*i)/10);
                    waveform->GetXaxis()->SetRangeUser(-10*ns,80*ns);
                    waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
                    waveform->GetYaxis()->SetTitleOffset(1.3);
                    waveform->GetXaxis()->SetTitle("Time [s]");
                    waveform->Draw("AL");
                    c3[i]->SetGrid();
                }
                
                clean = false;
            }
            
            if (clean){
                sprintf(Category,"Clean");
                if (j<50) {
                    c4[i]->cd();
                    
                    color4++;
                    if (color4==10 || color4==0) {
                        color4++;
                    } else if (color4>49) color4=2;
                    waveform->SetLineColor(color4);
                    
                    if (color4>1) {
                        waveform->Draw("SAME");
                    }else{
                        sprintf(graph_title,"Clean pulse %s",vol_folders[i].Data());
                        waveform->SetTitle(graph_title);
                        waveform->GetYaxis()->SetRangeUser(-0.1,pe*(15+3*i)/10);
                        waveform->GetYaxis()->SetTitle("Oscilloscope Signal [V]");
                        waveform->GetYaxis()->SetTitleOffset(1.3);
                        waveform->GetXaxis()->SetTitle("Time [s]");
                        waveform->Draw("AL");
                        c4[i]->SetGrid();
                    }
                }
                
                
            }
            
            tree->Fill();
            
            Event ++;
            if (Event%500==0) {
                cout<<"****----->Events analyzed:"<< Event << endl;
            }
            
            event_cnt++;
            
            
        }
        
        Correl_noise[0]->SetPoint(i,V_meas,direct_xtalk_pulse_cnt/event_cnt*100);
        Correl_noise[1]->SetPoint(i,V_meas,after_pulse_cnt/event_cnt*100);
        Correl_noise[2]->SetPoint(i,V_meas,xtalk_pulse_cnt/event_cnt*100);
        Correl_noise[3]->SetPoint(i,V_meas,
                                  Correl_noise[0]->GetY()[i]+Correl_noise[1]->GetY()[i]+Correl_noise[2]->GetY()[i]);
        
        
        //Uncoment to save waveforms and check them
        /*
        sprintf(canvas_title,"Plots/Immcrosstalk_%s.pdf",vol_folders[i].Data());
        c1[i]->Print(canvas_title,"pdf");
        sprintf(canvas_title,"Plots/Immcrosstalk_%s.root",vol_folders[i].Data());
        c1[i]->SaveAs(canvas_title,"root");
        sprintf(canvas_title,"Plots/Delcrosstalk_%s.pdf",vol_folders[i].Data());
        c2[i]->Print(canvas_title,"pdf");
        sprintf(canvas_title,"Plots/Afterpulse_%s.pdf",vol_folders[i].Data());
        c3[i]->Print(canvas_title,"pdf");
        sprintf(canvas_title,"Plots/Clean_%s.pdf",vol_folders[i].Data());
        c4[i]->Print(canvas_title,"pdf");
        */
        
    }
    
    tree->Write();
    
    TCanvas* c5 = new TCanvas("Correlated Noise","Correlated Noise",100,100,900,700);
    
    Correl_noise[3]->SetTitle("Total");
    Correl_noise[3]->SetMarkerColor(2);
    Correl_noise[3]->SetLineColor(2);
    Correl_noise[3]->GetYaxis()->SetRangeUser(0,30);
    Correl_noise[3]->GetYaxis()->SetTitle("Noise [%]");
    Correl_noise[3]->GetXaxis()->SetTitle("OverVoltage [V]");
    Correl_noise[3]->Draw("ALP*");
    
    Correl_noise[0]->SetTitle("Direct Cross-Talk");
    Correl_noise[1]->SetTitle("After Pulse");
    Correl_noise[2]->SetTitle("Delayed Cross-Talk");
    Correl_noise[0]->SetLineColor(3);
    Correl_noise[1]->SetLineColor(4);
    Correl_noise[2]->SetLineColor(5);
    Correl_noise[0]->SetMarkerColor(3);
    Correl_noise[1]->SetMarkerColor(4);
    Correl_noise[2]->SetMarkerColor(5);
    Correl_noise[0]->Draw("LP*");
    Correl_noise[1]->Draw("LP*");
    Correl_noise[2]->Draw("LP*");
    
    c5->BuildLegend();
    c5->SetGrid();
    c5->Print("Correlated Noise.pdf","pdf");
    c5->Write();
    
    delete hfile;
    
    return 0;

}







