
//
//  Noise_analysis.C
//  
//
//  Created by Felipe Gilberto Ortega on 25/04/16.
//
//

#include "Noise_analysis.h"



static const char *optString = "d:S:ah?";

int main(int argc, char* argv[]) {
    
    
    // initialize globalArgs
    globalArgs.data_folder = " ";
    globalArgs.arg_pathToSetupFile = " ";
    globalArgs.save_all = 0;
    
    // Get paremeter from the command
    int opt =0;
    opt = getopt(argc, argv, optString);
    if(opt == -1){
        std::cerr <<  "There is no opption in the command! Type \"output -h\" for help." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    while(opt != -1){
        switch(opt){
            case 'd':
                globalArgs.data_folder= optarg;
                //std::cout<<"-p option path= "<<globalArgs.arg_pathToData<<std::endl;
                break;
            case 'S':
                globalArgs.arg_pathToSetupFile = optarg;
                break;
            case 'a':
                globalArgs.save_all = 1;
                break;
            case 'h':
            case '?':
                std::cerr << "Usage: output -d pathToData -S pathToSetupFile [-a]" << std::endl;
                std::cerr << "----------------------------------------------------------------------------------------------------"<<std::endl;
                std::cerr << " '-d'+'-S' options are necessary!"<<std::endl;
                std::cerr << "-----------------------------------------------------------------------------------------------------"<<std::endl;
                std::cerr << " use '-a' option afterwards to save all the plots of the analysis to further check."<<std::endl;
                std::cerr << "-----------------------------------------------------------------------------------------------------"<<std::endl;
                std::cerr << "Example: ./output -d /Users/Analysis_waveforms/ov_scan_pde_H2014 -S /Users/Analysis_waveforms/config_file.txt [-a]"<<std::endl;
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
        opt = getopt(argc, argv, optString);
    }
    
    
    if((strncmp(globalArgs.data_folder," ",1) == 0|| strncmp(globalArgs.arg_pathToSetupFile," ",1) == 0)){
        std::cerr << "ERROR: -d or -S option is not set! Both of them has to be set correctly!"<<std::endl;
        exit(EXIT_FAILURE);
    }
            
    ifstream setupFile(globalArgs.arg_pathToSetupFile);
    if(!setupFile){
        std::cerr << "Failure: could not open file: \"" << globalArgs.arg_pathToSetupFile << "\"." << std::endl;
        std::cerr << "Please check if the path is correct or not!" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    
    //Read setup file:
    string s;
    vector <TString> vol_folders;
    Int_t data_size;
    while (true) {
        
        getline(setupFile, s);
        if (setupFile.eof()) break;
        
        const char* searchString = s.c_str();
        char volt [20];
        Int_t numfiles;
        
        if (s.find("#") == 0 || s=="") {
            continue; // Skip commented  or empty lines
        }
        
        //Find the voltages
        if(sscanf(searchString, "V || %s ||", volt)==1){
            vol_folders.push_back(volt);
        }
        
        //Find data size
        if(sscanf(searchString, "Files at each voltage || %d ||", &numfiles)==1){
            data_size = numfiles;
        }
    }
    
    //Initialize variables
    const Int_t vol_size = vol_folders.size();
    
    int singleplot=0;
    Int_t Event=0;
    Char_t Category[15];
    TGraph* waveform = 0;
    Double_t Amp;
    Double_t V_meas;
    
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
    
    ////////////////
    //Define thresholds
    ////////////////
    
    //in  nanoseconds
    const double reject_time = 4;
    const double time_dist_th = 0.4;
    
    //in percentage of pe
    const double after_pulse_th = 0.38;
    const double direct_xtalk_th = 1.17;
    const double xtalk_th = 0.85;
    
    
    ////////////////
    ////////////////
    
    /*const char * Voltage="56.5V";
    int event=0;
    if (singleplot) {
        single_plot(Voltage,event);
    }*/
    
    
    //Create a root tree with the graph of the waveform of each event and
    //classify them
    TString filename = "noiseanalysis.root";
    
    TFile *hfile = 0;
    hfile = TFile::Open(filename,"RECREATE");
    
    
    TTree *tree = new TTree("T","Noise Analysis");
    tree->Branch("Event",&Event,"Event/I");
    tree->Branch("Category",Category,"Category/C");
    //Uncomment if every single waveform is desired to be saved by its own on the root file
    //tree->Branch("waveform","TGraph",&waveform);
    tree->Branch("V_meas",&V_meas,"V_meas/D"); //OV of the measurement
    
   
    TGraph* Correl_noise[4];
    Correl_noise[0] = new TGraph();
    Correl_noise[1] = new TGraph();
    Correl_noise[2] = new TGraph();
    Correl_noise[3] = new TGraph();
    TGraph *Expfit_longtau[vol_size];
    TGraph *Expfit_AP[vol_size];
    
    //Fiting functions of long tau and AP recharge
    TF1 *exp_longtau= new TF1("exptau","[0]*exp(-x/[1])",0,180 * ns);
    TF1 *exp= new TF1("exp","[0]*(1-exp(-x/[1]))+[2]*exp(-x/[3])",0,180 * ns);
    
    TCanvas* c1[vol_size];
    TCanvas* c2[vol_size];
    TCanvas* c3[vol_size];
    TCanvas* c4[vol_size];
    
    TMultiGraph *Cleanwaves[vol_size];
    TCanvas* expfit_longtau_c[vol_size];
    TCanvas* expfit_AP_c[vol_size];
    
    cout<<"////////////"<< endl;
    cout<<"****----->Voltage Breakdown calculation ***"<< endl;
    
   
    vector <Double_t> pe_volt;
    TGraph *Vbias_ver= new TGraph();
    
    
        //Change to not recalculate the pe
        //pe_volt.push_back(6.87435e-02);
        /*pe_volt.push_back( 1.20426e-01);
        pe_volt.push_back(1.75262e-01);
        pe_volt.push_back(2.30936e-01);
        pe_volt.push_back(2.87958e-01);*/
        //pe_volt.push_back( 3.44156e-01);
        //Double_t VBD=55.9006;
    
    //Calculate Voltage breakdown and value of pe
    for (int i=0; i<vol_size; i++) {
        pe_volt.push_back(Amplitude_calc(vol_folders.at(i).Data(), data_size));
        V_meas = vol_folders.at(i).Atof();
        Vbias_ver->SetPoint(i, pe_volt.at(i), V_meas);
    }
    
    TCanvas* ca= new TCanvas("Voltage Breakdown calculation","Voltage Breakdown calculation",100,100,900,700);
    Vbias_ver->SetTitle("Voltage Breakdown calculation");
    Vbias_ver->GetYaxis()->SetTitle("Bias Volatge [V]");
    Vbias_ver->GetYaxis()->SetTitleOffset(1.2);
    Vbias_ver->GetXaxis()->SetTitle("Mean peak amplitude [V]");
    Vbias_ver->Draw("AP*");
    ca->SetGrid();
    
    cout<<"////////////"<< endl;
    cout<<"****----->Voltage Breakdown fit ***"<< endl;
    
    TFitResultPtr fit = Vbias_ver->Fit("pol1","S");
    Double_t VBD= fit->Value(0);
    
    if (globalArgs.save_all==1) ca->Write();
    
    
    cout<<"////////////"<< endl;
    cout<<"****----->Noise analysis ***"<< endl;
    cout<<"////////////"<< endl;
    
    /////////////////
    // Loop over all Voltages measured
    /////////////////
    for (int i=0; i<vol_size; i++) {
        
        //Important to reinitialize, the value color* = kOrange-11 is used to plot axis of TGraph()
        
        int color1 = kOrange-11;
        int color2 = kOrange-11;
        int color3 = kOrange-11;
        int color4 = kOrange-11;
        
        direct_xtalk_pulse_cnt = 0;
        xtalk_pulse_cnt = 0;
        after_pulse_cnt = 0;
        event_cnt = 0; //Events on the Voltage measured
        
        cout<<"****----->Voltage analyzed: "<< vol_folders.at(i) << endl;
        
        
        //Define amplitude measured at which OV
        
        Double_t pe = pe_volt.at(i);
        
        V_meas = vol_folders.at(i).Atof()-VBD;
        
        //Define canvases to save and check results
        
        Char_t canvas_title[40];
        sprintf(canvas_title,"Direct CrossTalk OV = %2.2f V",V_meas);
        c1[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        sprintf(canvas_title,"Delayed CrossTalk OV = %2.2f V",V_meas);
        c2[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        sprintf(canvas_title,"After Pulse OV = %2.2f V",V_meas);
        c3[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        sprintf(canvas_title,"Clean OV = %2.2f V",V_meas);
        c4[i] = new TCanvas(canvas_title,canvas_title,100,100,900,700);
        
        Cleanwaves[i]=new TMultiGraph();
        
        sprintf(canvas_title,"Exponential fit, #tau_l OV = %2.2f V",V_meas);
        expfit_longtau_c[i] = new TCanvas(canvas_title,canvas_title,300,100,900,500);
        sprintf(canvas_title,"Exponential fit OV = %2.2f V",V_meas);
        expfit_AP_c[i] = new TCanvas(canvas_title,canvas_title,300,100,900,500);
        
        Expfit_longtau[i]= new TGraph();
        Expfit_AP[i]= new TGraph();
        
        //loop over every measurement on a folder
        for (int j=0; j<data_size; j++) {
            
            Char_t datafilename[200];
            Char_t datashortfilename[100];
            
            sprintf(datafilename,"%s/%s/C1H%05i.csv",globalArgs.data_folder,vol_folders.at(i).Data(),j);
            sprintf(datashortfilename,"%s_C1H%05i",vol_folders.at(i).Data(),j);
            //Get the data of a single file:
            waveform = new TGraph(datafilename,"%lg %lg","/t;,");
            waveform->SetName(datashortfilename);
            waveform->SetTitle("");
            
            Int_t ROWS_DATA = waveform->GetN();
            Double_t *time = waveform->GetX();
            Double_t *volts = waveform->GetY();
            
            Amp = waveform->GetY()[0];
            
            
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
            
            /////////////////////////////////////////////////////
            // direct x-talk
            for (row = 0; row < ROWS_DATA; row++) {
                if ((time[row]>0 * ns)&(volts[row] > direct_xtalk_th * pe)) {// time larger 0ns
                    direct_xtalk_pulse++;
                }
            }
            
            /////////////////////////////////////////////////////
            // after-pulse threshold
            for (row = 0; row < ROWS_DATA; row++) {
                if ((time[row]>reject_time*ns)&(volts[row] > after_pulse_th * pe)) {// time larger 4ns and ap_th
                    after_pulse++;
                }
            }
            
            /////////////////////////////////////////////////////
            // delayed x-talk
            for (row = 0; row < ROWS_DATA; row++) {
                if ((time[row]>reject_time*ns)&(volts[row] > xtalk_th * pe)) {// time larger 4ns and larger xtalk_th
                    xtalk_pulse++;
                }
            }
            
            
            /////////////////////////////////////////////////////////////////////
            // Detect peaks in data after 4ns, count the number of maxima and
            // measure the time of arrival of first maxima, used later for AP exp fit
            /////////////////////////////////////////////////////////////////////
            max_noise_cnt = 0;
            for (row = 0; row < ROWS_DATA; row++) {
                if (time[row] > reject_time*ns) {// time larger 4ns
                    if (volts[row] > sig_max) {
                        sig_max = volts[row];        // set the max
                        time_of_max = time[row];    // time max
                        max_noise_cnt++;                   // set the histeresis cnt
                    }else if (max_noise_cnt > 0) max_noise_cnt--;  //  count down if no new max is reached
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
                } // 4ns
            } //loop over time
            
            bool clean = true; //The pulse is clean until the contrary can be demonstrated
            char graph_title[50];
            
            //Check for imm x-talk and plot
            if (direct_xtalk_pulse > 0){
                direct_xtalk_pulse_cnt++;
                sprintf(Category,"ImmCrosstalk");
                c1[i]->cd();
                
                //Set graph color, and counting to draw axis and title
                color1=color1+2;
                if (color1>kOrange+110) {
                    color1=kOrange-8;
                }else if (color1>kOrange+109){
                    color1=kOrange-7;
                }
                waveform->SetLineColor(color1);
                waveform->SetMarkerColor(color1);
                
                //Format the graph
                sprintf(graph_title,"Direct CrossTalk OV = %2.2f V",V_meas);
                waveform = format_graph(waveform,graph_title,pe);
                
                if (color1>kOrange-8) {
                    waveform->Draw("SAME");
                }else{
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
                
                //Set graph color, and counting to draw axis and title
                color2=color2+2;
                if (color2>kOrange+110) {
                    color2=kOrange-8;
                }else if (color2>kOrange+109){
                    color2=kOrange-7;
                }
                waveform->SetLineColor(color2);
                waveform->SetMarkerColor(color2);
                
                //Format the graph
                sprintf(graph_title,"Delayed cross-talk OV = %2.2f V",V_meas);
                waveform = format_graph(waveform,graph_title,pe);

                if (color2>kOrange-8) {
                    waveform->Draw("SAME");
                }else{
                    waveform->Draw("AL");
                    c2[i]->SetGrid();
                }
                clean = false;
            }
            
            //  Only after pulse
            if (after_pulse > 0 && xtalk_pulse == 0 && direct_xtalk_pulse == 0){
                after_pulse_cnt++;
                sprintf(Category,"AfterPulse");
                c3[i]->cd();
                
                //Set graph color, and counting to draw axis and title
                color3=color3+2;
                if (color3>kOrange+110) {
                    color3=kOrange-8;
                }else if (color3>kOrange+109){
                    color3=kOrange-7;
                }
                waveform->SetLineColor(color3);
                waveform->SetMarkerColor(color3);
                
                //Format the graph
                sprintf(graph_title,"After pulse OV = %2.2f V",V_meas);
                waveform = format_graph(waveform,graph_title,pe);
                
                if (color3>kOrange-8) {
                    waveform->Draw("SAME");
                }else{
                    waveform->Draw("AL");
                    c3[i]->SetGrid();
                }
                
                clean = false;
                
                //Fill for the exponential fit
                Expfit_AP[i]->SetPoint(after_pulse_cnt-1,time_of_max,sig_max);
            }
            
            // Only clean graphs for the sample
            if (clean){
                sprintf(Category,"Clean");
                
                if (color4 < 860 && j <100) { //Max 100 clean graphs on the plot
                    Cleanwaves[i]->Add(waveform);
                    c4[i]->cd();
                    
                    //Set graph color, and counting to draw axis and title
                    color4=color4+2;
                    if (color4>kOrange+110) {
                        color4=kOrange-8;
                    }else if (color4>kOrange+109){
                        color4=kOrange-7;
                    }
                    waveform->SetLineColor(color4);
                    waveform->SetMarkerColor(color4);
                    
                    //Format the graph
                    sprintf(graph_title,"Clean pulse OV = %2.2f V",V_meas);
                    waveform = format_graph(waveform,graph_title,pe);
                    
                    if (color4>kOrange-8) {
                        waveform->Draw("SAME");
                    }else{
                        waveform->Draw("AL");
                        c4[i]->SetGrid();
                    }
                }
                
                
            }
            
            tree->Fill();
            
            Event ++;//Total number of events analyzed on the run
            if (Event%500==0) {
                cout<<"****----->Events analyzed:"<< Event << endl;
            }
            
            event_cnt++;
            
            
        }
        
        cout<<"////////////"<< endl;
        cout<<"****----->Long tau fit ***"<< endl;
        expfit_longtau_c[i]->cd();
        Cleanwaves[i]->Draw("AP*");
        // Fit parameters and limits to calculate slow component of the pulse
        exp_longtau->SetParameter(0,pe*0.2);
        exp_longtau->SetParLimits(0,0.05*pe,0.5*pe);
        exp_longtau->SetParameter(1,80*ns);
        exp_longtau->SetParLimits(1,4*ns,200*ns);
        Cleanwaves[i]->Fit("exptau","","",4*ns,60*ns); // Fit boundaries for the slow component of the pulse
        Double_t amp0 = exp_longtau->GetParameter(0);
        Double_t tau = exp_longtau->GetParameter(1);
        if (globalArgs.save_all==1) expfit_longtau_c[i]->Write();
        
        c4[i]->cd();
        TF1* exp_tau_plot =(TF1*) exp_longtau->Clone();
        exp_tau_plot->Draw("SAME");//Draw fit-line over clean waveforms
        
        cout<<"////////////"<< endl;
        cout<<"****----->After pulse fit ***"<< endl;
        expfit_AP_c[i]->cd();
        Expfit_AP[i]->Draw("AP*");
        // Fit parameters and limits to calculate AP recharge
        exp->SetParameter(0,pe);
        exp->SetParLimits(0,0.5*pe,1.5*pe);
        exp->SetParameter(1,30*ns);
        exp->SetParLimits(1,4*ns,500*ns);
        exp->SetParameter(2,amp0);
        exp->FixParameter(2,amp0);
        exp->SetParameter(3,tau);
        exp->FixParameter(3,tau);
        Expfit_AP[i]->Fit("exp");
        if (globalArgs.save_all==1) expfit_AP_c[i]->Write();
        
        c3[i]->cd();
        TF1* exp_plot =(TF1*) exp->Clone();
        exp_plot->Draw("SAME"); //Draw fit-line over AP waveforms
        
        
        //Final result: Correlated noise
        Correl_noise[0]->SetPoint(i,V_meas,direct_xtalk_pulse_cnt/event_cnt*100);
        Correl_noise[1]->SetPoint(i,V_meas,after_pulse_cnt/event_cnt*100);
        Correl_noise[2]->SetPoint(i,V_meas,xtalk_pulse_cnt/event_cnt*100);
        Correl_noise[3]->SetPoint(i,V_meas,
                                  Correl_noise[0]->GetY()[i]+Correl_noise[1]->GetY()[i]+Correl_noise[2]->GetY()[i]);
        
        
        //Save/print reults:
        if (globalArgs.save_all==1){
            c1[i]->Write();
            c2[i]->Write();
            c3[i]->Write();
            c4[i]->Write();
        }
        
        sprintf(canvas_title,"Plots/Immcrosstalk_%s.pdf",vol_folders.at(i).Data());
        c1[i]->Print(canvas_title,"pdf");
        sprintf(canvas_title,"Plots/Delcrosstalk_%s.pdf",vol_folders.at(i).Data());
        c2[i]->Print(canvas_title,"pdf");
        sprintf(canvas_title,"Plots/Afterpulse_%s.pdf",vol_folders.at(i).Data());
        c3[i]->Print(canvas_title,"pdf");
        sprintf(canvas_title,"Plots/Clean_%s.pdf",vol_folders.at(i).Data());
        c4[i]->Print(canvas_title,"pdf");
        
        
    }
    //Save TTree with hist of noise
    //Save each event with its OV and the noise classification
    tree->Write();
    
    //Create final plot of total correlated noise
    TCanvas* c5 = new TCanvas("Correlated Noise","Correlated Noise",100,100,900,700);
    
    Correl_noise[3]->SetTitle("Correlated Noise");
    Correl_noise[3]->SetMarkerColor(kRed);
    Correl_noise[3]->SetLineColor(kRed);
    Correl_noise[3]->GetYaxis()->SetRangeUser(0,20);
    Correl_noise[3]->GetYaxis()->SetTitle("Noise [%]");
    Correl_noise[3]->GetXaxis()->SetTitle("OverVoltage [V]");
    Correl_noise[3]->Draw("ALP*");
    
    Correl_noise[0]->SetTitle("Direct Cross-Talk");
    Correl_noise[1]->SetTitle("After Pulse");
    Correl_noise[2]->SetTitle("Delayed Cross-Talk");
    Correl_noise[0]->SetLineColor(kBlue);
    Correl_noise[1]->SetLineColor(kOrange+7);
    Correl_noise[2]->SetLineColor(kGreen+2);
    Correl_noise[0]->SetMarkerColor(kBlue);
    Correl_noise[1]->SetMarkerColor(kOrange+7);
    Correl_noise[2]->SetMarkerColor(kGreen+2);
    Correl_noise[0]->Draw("LP*");
    Correl_noise[1]->Draw("LP*");
    Correl_noise[2]->Draw("LP*");
    
    TLegend* leg = new TLegend(0.15,0.65,0.47,0.87);
    leg->AddEntry(Correl_noise[3],"Total","lp");
    leg->AddEntry(Correl_noise[0],"Direct Cross-Talk","lp");
    leg->AddEntry(Correl_noise[1],"After Pulse","lp");
    leg->AddEntry(Correl_noise[2],"Delayed Cross-Talk","lp");
    leg->Draw();
    
    
    c5->SetGrid();
    c5->Print("Plots/Correlated Noise.pdf","pdf");
    c5->Write();
    
    delete hfile;
    
    return 0;

}
