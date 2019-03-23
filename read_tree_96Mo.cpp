// ************************************************************************//
// Frank's nice sorting code :)
// Modified by Cecilie, March 2019, for the 96Mo data set
// 
// Info on the experiment: 96Mo(p,p'gamma) with E_p = 16 MeV, I_p = 3.6-3.7 nA
// SiRi in backward angles, 126-140 degrees
// OSCAR in the close configuration, 16.3 cm - 17 cm from the target center position
// One OSCAR detector (no. 30) was not working properly, "jumping" in the energy
// The target was 96Mo, 1.94 mg/cm^2 and 96.7% enriched.
// The validation signal for the XiA boxes were the E detectors, 
// looking 1 us forward and 1 us backward in time, 
// so in total the time window is 2 us
//
// Data files, 96Mo: 
// sirius-20190314-084105.data and sirius-20190314-084105-big00X.data
// sirius-20190314-145712.data and sirius-20190314-145712-big-000.data
// Data files, 28Si (for calibration):
// sirius-20190313-142346.data and sirius-20190313-142346-big-00X.data  
// 
// Elog: Nd-p-2019
// ************************************************************************//


#include <iostream>
#include <stdio.h>
#include <string>
#define MAX_HITS_EDET 8    //!< Maximum number of hits in e detectors per event
#define MAX_HITS_LABR 64    //!< Maximum number of hits in labr detectors per event

void read_tree_96Mo(){
    
    Long64_t ns_to_hours = 3600000000000;
    Long64_t ns_to_minutes = 60000000000;
    Long64_t ns_to_seconds =  1000000000;
    
    TChain *t = (TChain *) new TChain("data");
    
    // Character string for the file name
    char name[1024];
    
    // Set file name, open the file and add it to the TChain t.
    sprintf(name, "sirius-20190314-084105.root");
    ifstream rootfile1_96Mo(name);
    if (rootfile1_96Mo){
        t->Add(name);
        std::cout << "Added file " << name << " to chain." << std::endl;
    } else {
        std::cout << "File not read" << std::endl;
    }
    
    // Variables for detector IDs, timestamps, energies etc
    Int_t deDet_ID=0, tel_ID=0;
    Long64_t deDet_timestamp=0;
    Double_t deDet_energy=0;
    
    Int_t deDet_mult=0, eDet_mult=0;
    Double_t eDet_time[MAX_HITS_EDET], eDet_energy[MAX_HITS_EDET];
    
    Int_t labr_mult=0, labr_ID[MAX_HITS_LABR], labr_ring[MAX_HITS_LABR];
    Double_t labr_time[MAX_HITS_LABR], labr_energy[MAX_HITS_LABR];
    
    // TO BE UPDATED FOR 96Mo
    // Conversion coefficients to go from deposited energy in SiRi
    // to excitation energy in the residual nucleus
    Double_t Ex_coeff [8][3] = {//from qkinz for 96Mo(p,p')
        {16647.790,-0.997342,-6.39e-07},
        {16650.309,-0.996624,-6.50e-07},
        {16651.622,-0.995680,-6.70e-07},
        {16653.284,-0.994807,-6.86e-07},
        {16654.474,-0.993855,-7.04e-07},
        {16655.149,-0.992823,-7.24e-07},
        {16655.258,-0.991707,-7.46e-07},
        {16654.735,-0.990502,-7.70e-07}
    };
    
    // Loop over the E_back detectors that gave a signal
    for (int i=0; i<MAX_HITS_EDET ; ++i) {
        eDet_energy[i]=0;
        eDet_time[i]=0;
    }
    
    // Loop over the LaBr3 detectors that gave a signal
    for (int i=0; i<MAX_HITS_LABR ; ++i) {
        labr_ID[i]=0;
        labr_energy[i]=0;
        labr_ring[i]=0;
        labr_time[i]=0;
    }
    
    // Set branch addresses in the TChain
    t->SetBranchAddress("deDet_ID", &deDet_ID); // Front detector ID
    t->SetBranchAddress("tel_ID", &tel_ID); // Back detector ID
    t->SetBranchAddress("deDet_energy", &deDet_energy);
    t->SetBranchAddress("deDet_timestamp", &deDet_timestamp);
    
    t->SetBranchAddress("eDet_mult", &eDet_mult);
    t->SetBranchAddress("eDet_energy", eDet_energy);
    t->SetBranchAddress("eDet_time", eDet_time);
    
    t->SetBranchAddress("labr_mult", &labr_mult);
    t->SetBranchAddress("labr_ID", labr_ID);
    t->SetBranchAddress("labr_ring", labr_ring);
    t->SetBranchAddress("labr_energy", labr_energy);
    t->SetBranchAddress("labr_time", labr_time);
    
    // Define matrices and histograms
    const int max_e = 20000, max_de = 10000; // Max energies in the E and Delta E detectors
    const int max_e_labr3 =  20000; // Max energy for the LaBr3 detectors
    
    TH1D *h_deDet_mult = new TH1D("h_eDet_mult","#Delta E_{front} multiplicity",64,0,64);
    TH1D *h_eDet_mult = new TH1D("h_eDet_mult","E_{back} multiplicity",8,0,8);
    TH1D *h_LaBr3_mult = new TH1D("h_LaBr3_mult","LaBr3 multiplicity",60,0,60); // Potentially more hits than 30, because the validation signal gives a 2 us time wndow

    TH2D *h_deDet_energy = new TH2D("h_deDet_energy","#Delta E_{front} energy",1000,0,max_de,64,0,64);
    TH2D *h_eDet_energy = new TH2D("h_eDet_energy","E_{back} energy",1000,0,max_e,8,0,8);
    TH2D *h_LaBr3_energy = new TH1D("h_LaBr3_energy","LaBr3 energy",40000,0,max_e_labr3);

    TH2D *h_deDet_time = new TH2D("h_deDet_time","#Delta E_{front} time",20000,-1000,1000,8,0,8);
    TH2D *h_eDet_time = new TH2D("h_eDet_time","E_{back} time",20000,-1000,1000,8,0,8);
    TH2D *h_LaBr3_time = new TH2D("h_LaBr3_time","LaBr3 time",20000,-1000,1000,30,0,30);

    // Cecilie, modification 22 March 2019:
    // Making Delta E - E matrices as in the old user_sort.cpp by Alexander Bürger
    // These can be used for calibration with the script peaks2D.C by Alexander Bürger
    // TO DO: 
    // 1. Make sure about the mapping of the detector IDs
    // 2. Fill the matrices and write them to the output file
    TH2D *deltaE_E_matrices[8][8] = {NULL};

    for(int b=0; b<8; ++b ) {
        for(int f=0; f<8; ++f ) {
            ostringstream histogramNameStream;
            histogramNameStream << "m_e_de_b" << b << "f" << f;
            deltaE_E_matrices[b][f] = new TH2D(histogramNameStream.str().c_str(),histogramNameStream.str().c_str(),1000,0,max_e,1000,0,max_de);
            deltaE_E_matrices[b][f]->GetXaxis()->SetTitle("E_{back} (keV)");
            deltaE_E_matrices[b][f]->GetYaxis()->SetTitle("#Delta E_{front} (keV)");
            cout << histogramNameStream.str().c_str() << endl;
        }
    }

    // ************************************************************************//
    // Histograms that are not used at present, but will probably be used later:
    //TH2D *h_eDet_energy_time = new TH2D("h_eDet_energy_time","eDet Energy-Time",8192,0,32768,10000,-1000,1000);
    
    //TH2D *h_Ex_gamma_with_bg = new TH2D("h_Ex_gamma_with_bg","Excitation Energy vs gamma Energy with Background",16384,0,16384,2048,0,16384);
    //TH2D *h_Ex_gamma_bg = new TH2D("h_Ex_gamma_bg","Excitation Energy vs gamma Energy Background",16384,0,16384,2048,0,16384);
    
    //TH2D *h_labr_energy = new TH2D("h_labr_energy","labr Energy",32768,0,32768,30,1,31);
    //TH1D *h_labr_with_bg = new TH1D("h_labr_with_bg","labr with Background",32768,0,32768);
    //TH1D *h_labr_bg = new TH1D("h_labr_bg","labr Background",32768,0,32768);
    //TH2D *h_labr_gamma_gamma = new TH2D("h_labr_gamma_gamma","labr gamma-gamma",8192,0,32762,8192,0,32762);
    //TH3D *h_labr_time = new TH3D("h_labr_time","labr Time",20000,-1000,1000,64,1,65,30,1,31);
    //TH2D *h_labr_energy_time = new TH2D("h_labr_energy_time","labr Energy-Time",1024,0,1024,40000,-1000,1000);
    //TH2D *h_labr_energy_time = new TH2D("h_labr_energy_time","labr Energy-Time",8192,0,8192,20000,-1000,1000);
    // ************************************************************************//
 
    // Possibility for reading in files with pre-defined graphical cuts - not used at present
    /*
    TFile *cut_file = new TFile("163Dy_cuts.root");
    TCutG *cut_ban = (TCutG*) cut_file->Get("p_cut");
    TCutG *cut_eDet = (TCutG*) cut_file->Get("eDet_ET_Ex_cut");//eDet_ET_p_cut_delayed //eDet_ET_p_cut
    TCutG *cut_labr = (TCutG*) cut_file->Get("labr_ET_p_cut_prompt");
    TCutG *cut_bg = (TCutG*) cut_file->Get("labr_ET_p_cut_bg");
    //TCutG *cut_Ex = (TCutG*) cut_file->Get("Ex_p_cut");
    cut_file->Close();
    */
    
    // Now preparing for reading the events in the data file
    Long64_t nrow = (Long64_t) t->GetEntries(); // Get how many events there are in the file
    Long64_t n_eDet=0, n_eDet_bg=0, n_in_cut=0;
    

    // HERE COMES THE MAIN LOOP, READING THE EVENTS
    for (Long64_t n=0; n<nrow; n++){      
        t->GetEntry(n); // Read event by event

        // Delta E multiplicity
        for (int i=0; i<deDet_mult ; ++i) {
            h_deDet_mult->Fill(eDet_mult,deDet_ID);

        
        for (int i=0; i<eDet_mult ; ++i) {
                h_de_e_energy->Fill(eDet_energy[i],deDet_energy);
                for (int j=0; j<labr_mult ; ++j) {
                    if (labr_ID[j]==1) {
                        h_labr_time->Fill(labr_time[j]);
                        if(labr_time[j]>-9 && labr_time[j]<1){
                            h_Ex_gamma_with_bg->Fill(labr_energy[j],eDet_energy[i]+deDet_energy);
                        }
                        else if(labr_time[j]>-65 && labr_time[j]<-55){
                            h_Ex_gamma_bg->Fill(labr_energy[j],eDet_energy[i]+deDet_energy);
                        }
                    }
                }
            }
        
            //h_deDet_energy->Fill(deDet_energy,deDet_ID);
            for (int i=0; i<eDet_mult ; ++i) {
                //if (cut_ban->IsInside(eDet_energy[i],deDet_energy) && cut_eDet->IsInside(eDet_energy[i],eDet_time[i]) && (Ex_coeff[(deDet_ID-tel_ID)%8][0]+(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(eDet_energy[i]+deDet_energy)*(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][2]-215.406)<6551) {
                     //h_deDet_energy->Fill(deDet_energy);
                     h_eDet_energy->Fill(eDet_energy[i],tel_ID);
                     //h_eDet_time->Fill(eDet_time[i]);
                     //h_de_e_energy->Fill(eDet_energy[i],deDet_energy);
                     //h_eDet_energy_time->Fill(eDet_energy[i],eDet_time[i]);
                     //h_eDet_mult->Fill(eDet_mult);
                     //break;
                     for (int j=0; j<labr_mult ; ++j) {
                         //h_labr_mult->Fill(labr_mult);
                         //h_labr_energy_time->Fill(labr_energy[j],labr_time[j]);
                         //h_labr_time->Fill(labr_time[j],deDet_ID,labr_ID[j]);
                         //h_labr_energy->Fill(labr_energy[j],labr_ID[j]);
                         if (labr_ID[j]==1) {
                             //h_labr_time->Fill(labr_time[j]);
                             h_labr_energy->Fill(labr_energy[j]);
                             h_labr_energy_time->Fill(labr_energy[j],labr_time[j]);
                         }
                         //if (cut_Ex->IsInside(labr_energy[j],Ex_coeff[(deDet_ID-tel_ID)%8][0]+(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(eDet_energy[i]+deDet_energy)*(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][2]-215.406)/* && cut_labr->IsInside(labr_energy[j],labr_time[j])*/) {
                             //h_labr_energy_time->Fill(labr_energy[j],labr_time[j]);
                             //h_labr_energy->Fill(labr_energy[j]);
                         //}
                         /*
                         //if (labr_time[j]>142 && labr_time[j]<172) { //for delayed_p_cut
                         if ((labr_time[j]>110 && labr_time[j]<130)||(labr_time[j]>150 && labr_time[j]<175)||(labr_time[j]>195 && labr_time[j]<220)||(labr_time[j]>240 && labr_time[j]<265)||(labr_time[j]>285 && labr_time[j]<310)||(labr_time[j]>330 && labr_time[j]<355)||(labr_time[j]>375 && labr_time[j]<400)||(labr_time[j]>420 && labr_time[j]<445)||(labr_time[j]>465 && labr_time[j]<490)) { //for delayed_d_cut
                             //h_labr_energy->Fill(labr_energy[j]);
                             h_labr_energy->Fill(labr_energy[j],labr_ID[j]);
                         }
                         */
                         
                         //if (cut_labr->IsInside(labr_energy[j],labr_time[j]) && labr_energy[j]>305 && labr_energy[j]<=(Ex_coeff[(deDet_ID-tel_ID)%8][0]+(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(eDet_energy[i]+deDet_energy)*(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][2]-215.406)) {
                             //h_Ex_gamma_with_bg->Fill(labr_energy[j],Ex_coeff[(deDet_ID-tel_ID)%8][0]+(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(eDet_energy[i]+deDet_energy)*(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][2]-215.406);
                             //h_labr_energy->Fill(labr_energy[j],labr_ID[j]);
                             //h_labr_energy->Fill(labr_energy[j]);
                             //h_labr_with_bg->Fill(labr_energy[j]);
                         //}
                         //else if (cut_bg->IsInside(labr_energy[j],labr_time[j]) && labr_energy[j]>305 && labr_energy[j]<=(Ex_coeff[(deDet_ID-tel_ID)%8][0]+(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(eDet_energy[i]+deDet_energy)*(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][2]-215.406)) {
                             //h_Ex_gamma_bg->Fill(labr_energy[j],Ex_coeff[(deDet_ID-tel_ID)%8][0]+(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(eDet_energy[i]+deDet_energy)*(eDet_energy[i]+deDet_energy)*Ex_coeff[(deDet_ID-tel_ID)%8][2]-215.406);
                         //}
                         
                     }
                    
                    n_eDet++;
                //}
                //h_eDet_energy->Fill(eDet_energy[i]);
                //h_eDet_time->Fill(eDet_time[i]);
                //h_de_e_energy->Fill(eDet_energy[i],deDet_energy);
            }
            if (n_eDet>1) {
                n_eDet_bg=n_eDet_bg+n_eDet-1;
            }
            n_in_cut=n_in_cut+n_eDet;
            n_eDet=0;
        //}
        
        /*h_deDet_ID->Fill(deDet_ID);
        h_tel_ID->Fill(tel_ID);
        h_deDet_energy->Fill(deDet_energy);
        h_deDet_timestamp->Fill(deDet_timestamp/ns_to_seconds);
 
        h_eDet_mult->Fill(eDet_mult);
        for (int i=0; i<eDet_mult ; ++i) {
            h_eDet_energy->Fill(eDet_energy[i]);
            h_eDet_time->Fill(eDet_time[i]);
        }
        
        h_labr_mult->Fill(labr_mult);
        for (int i=0; i<labr_mult ; ++i) {
            h_labr_ID->Fill(labr_ID[i]);
            h_labr_ring->Fill(labr_ring[i]);
            h_labr_energy->Fill(labr_energy[i]);
            h_labr_time->Fill(labr_time[i]);
        }*/
        
        if(n%1000000==0){
            std::cout << ".";
            std::cout.flush();
        }//*/
        
    }
    
    std::cout << endl;
    std::cout << "Total number of events:           " << nrow << endl;
    std::cout << "Number of events inside the cut:  " << n_in_cut << endl;
    std::cout << "Number of eDet background events: " << n_eDet_bg << endl;
    
    
    //TCanvas *c_de_e_energy = new TCanvas("c_de_e_energy","deDet eDet Energy",1);
    //h_de_e_energy->Draw("colz");
    
    
    TH2D *h_Ex_gamma = new TH2D(*h_Ex_gamma_with_bg);
    h_Ex_gamma->SetNameTitle("h_Ex_gamma","Excitation Energy vs gamma Energy");
    if (!(h_Ex_gamma->GetSumw2N() > 0)) h_Ex_gamma->Sumw2(kTRUE);
    h_Ex_gamma->Add(h_Ex_gamma_bg, -1.0);
    
    
    
    TFile *outputFile = new TFile("Mo96.root","update");
    //h_deDet_energy->Write("h_deDet_energy_peak4");
    h_eDet_energy->Write("h_eDet_energy",TObject::kOverwrite);
    //h_eDet_time->Write("h_eDet1_time_gs_19F_cut");
    h_labr_time->Write("h_labr_time",TObject::kOverwrite);//,TObject::kOverwrite
    h_labr_energy->Write("h_labr_energy",TObject::kOverwrite);
    //h_labr_gamma_gamma->Write("h_labr_gamma_gamma_peak2",TObject::kOverwrite);
    h_de_e_energy->Write("h_de_e_energy",TObject::kOverwrite);
    //h_eDet_energy_time->Write("h_eDet_energy_time_Ex_cut",TObject::kOverwrite);
    h_labr_energy_time->Write("h_labr_energy_time",TObject::kOverwrite);
    //h_eDet_mult->Write("h_eDet_mult_t_cut",TObject::kOverwrite);
    //h_Ex_gamma_with_bg->Write("h_Ex_gamma_with_bg",TObject::kOverwrite);
    h_Ex_gamma_bg->Write("h_Ex_gamma_bg",TObject::kOverwrite);
    h_Ex_gamma->Write("h_Ex_gamma",TObject::kOverwrite);
    //h_labr_mult->Write("h_labr_mult",TObject::kOverwrite);
    outputFile->Close();
    
}
