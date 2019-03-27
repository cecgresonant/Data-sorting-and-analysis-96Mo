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
#define MAX_HITS_LABR 64    //!< Maximum number of hits in labr detectors per event (can be more than 30 due to the long time range for the vaidation signal)

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
    Double_t deDet_energy; 
    
    Int_t deDet_mult=0, eDet_mult=0;
    Double_t eDet_time[MAX_HITS_EDET], eDet_energy[MAX_HITS_EDET];
    
    Int_t labr_mult=0, labr_ID[MAX_HITS_LABR], labr_ring[MAX_HITS_LABR];
    Double_t labr_time[MAX_HITS_LABR], labr_energy[MAX_HITS_LABR];
    
    // Conversion coefficients to go from deposited energy in SiRi
    // to excitation energy in the residual nucleus
    Double_t Ex_coeff [8][3] = {//from qkinz for 96Mo(p,p')
        {15382.792,-0.991240,-9.24e-07},
        {15383.826,-0.990678,-9.34e-07},
        {15383.509,-0.989872,-9.55e-07},
        {15383.491,-0.989155,-9.70e-07},
        {15382.895,-0.988357,-9.87e-07},
        {15381.673,-0.987474,-10.1e-07},
        {15379.764,-0.986502,-10.3e-07},
        {15377.091,-0.985435,-10.6e-07}
    };
    
    // Initializing vectors
    for (int i=0; i<MAX_HITS_EDET; ++i) {
        eDet_energy[i]=0;
        eDet_time[i]=0;
    }

    for (int i=0; i<MAX_HITS_LABR; ++i) {
        labr_ID[i]=0;
        labr_energy[i]=0;
        labr_ring[i]=0;
        labr_time[i]=0;
    }
    
    // Set branch addresses in the TChain
    t->SetBranchAddress("deDet_ID", &deDet_ID); // Front detector ID
    t->SetBranchAddress("tel_ID", &tel_ID); // Back detector ID
    t->SetBranchAddress("deDet_energy", &deDet_energy); // Delta E energy
    t->SetBranchAddress("deDet_timestamp", &deDet_timestamp); // Delta E timestamp
    
    t->SetBranchAddress("eDet_mult", &eDet_mult);   // E multiplicity
    t->SetBranchAddress("eDet_energy", eDet_energy); // E energy
    t->SetBranchAddress("eDet_time", eDet_time); // E time
    
    t->SetBranchAddress("labr_mult", &labr_mult); // LaBr3 multiplicity
    t->SetBranchAddress("labr_ID", labr_ID);    // LaBr3 ID
    t->SetBranchAddress("labr_ring", labr_ring);    // LaBr3 ring (6 different theta angles)
    t->SetBranchAddress("labr_energy", labr_energy); // LaBr3 energy
    t->SetBranchAddress("labr_time", labr_time);    // LaBr3 time
    
    // Define matrices and histograms
    const int max_e = 20000, max_de = 10000; // Max energies in the E and Delta E detectors
    const int max_e_labr3 =  20000; // Max energy for the LaBr3 detectors
    
    TH1D *h_eDet_mult = new TH1D("h_eDet_mult","E_{back} multiplicity",9,0,9);
    TH1D *h_LaBr3_mult = new TH1D("h_LaBr3_mult","LaBr3 multiplicity",60,0,60); // Potentially more hits than 30, because the validation signal gives a 2 us time wndow

    TH2D *h_deDet_energy = new TH2D("h_deDet_energy","#Delta E_{front} energy",1000,0,max_de,65,0,65);
    TH2D *h_eDet_energy = new TH2D("h_eDet_energy","E_{back} energy",1000,0,max_e,9,0,9);
    TH2D *h_LaBr3_energy = new TH2D("h_LaBr3_energy","LaBr3 energy",20000,0,max_e_labr3,32,0,32);

    TH2D *h_deDet_time = new TH2D("h_deDet_time","#Delta E_{front} time",20000,-1000,1000,65,0,65);
    TH2D *h_eDet_time = new TH2D("h_eDet_time","E_{back} time",20000,-1000,1000,9,0,9);
    TH2D *h_LaBr3_time = new TH2D("h_LaBr3_time","LaBr3 time",20000,-1000,1000,32,0,32);

    // Cecilie, modification 22 March 2019:
    // Making Delta E - E matrices as in the old user_sort.cpp by Alexander Bürger
    // These can be used for calibration with the script peaks2D.C by Alexander Bürger
    // TO DO: 
    // 1. Make sure about the mapping of the detector IDs
    // 2. Fill the matrices and write them to the output file
    TH2D *deltaE_E_matrices[64] = {NULL};

    int histo_no=0;
    for(int b=0; b<8; ++b ) {
        for(int f=0; f<8; ++f ) {
        	ostringstream histogramNameStream;
            histogramNameStream << "m_e_de_b" << b << "f" << f;
        	deltaE_E_matrices[histo_no] = new TH2D(histogramNameStream.str().c_str(),histogramNameStream.str().c_str(),1000,0,max_e,1000,0,max_de);
        	deltaE_E_matrices[histo_no]->GetXaxis()->SetTitle("E_{back} (keV)");
        	deltaE_E_matrices[histo_no]->GetYaxis()->SetTitle("#Delta E_{front} (keV)");
        	++histo_no;
        }
        //cout << histogramNameStream.str().c_str() << endl;
    }

    // ************************************************************************//
    // Histograms that are not used at present, but will probably be used later:
    //TH2D *h_eDet_energy_time = new TH2D("h_eDet_energy_time","eDet Energy-Time",8192,0,32768,10000,-1000,1000);
    
    //TH2D *h_Ex_gamma_with_bg = new TH2D("h_Ex_gamma_with_bg","Excitation Energy vs gamma Energy with Background",16384,0,16384,2048,0,16384);
    //TH2D *h_Ex_gamma_bg = new TH2D("h_Ex_gamma_bg","Excitation Energy vs gamma Energy Background",16384,0,16384,2048,0,16384);
    
    //TH1D *h_labr_with_bg = new TH1D("h_labr_with_bg","labr with Background",32768,0,32768);
    //TH1D *h_labr_bg = new TH1D("h_labr_bg","labr Background",32768,0,32768);
    //TH2D *h_labr_gamma_gamma = new TH2D("h_labr_gamma_gamma","labr gamma-gamma",8192,0,32762,8192,0,32762);
    //TH3D *h_labr_time = new TH3D("h_labr_time","labr Time",20000,-1000,1000,64,1,65,30,1,31);
    //TH2D *h_labr_energy_time = new TH2D("h_labr_energy_time","labr Energy-Time",1024,0,1024,40000,-1000,1000);
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

        h_deDet_energy->Fill(deDet_energy,deDet_ID-1);
        // E multiplicity
        h_eDet_mult->Fill(eDet_mult);
        // LaBr3 multiplicity
        h_LaBr3_mult->Fill(labr_mult);
        
        for (int i=0; i<eDet_mult; ++i) {
            n_eDet++;
            // Fill E detectors' energies and times
			h_eDet_energy->Fill(eDet_energy[i],tel_ID-1);
            h_eDet_time->Fill(eDet_time[i],tel_ID-1);

            // Fill the "banana" matrices, one for each back and front detector combination
            //std::cout << "E det ID = " << tel_ID << ", Delta E det ID = " << deDet_ID << std::endl;
            const int id_f = (deDet_ID-1) % 8; // from here, the front strip has ID 0..7
            const int id_b = (tel_ID-1);
            char *bananaName = new char[10];
            sprintf(bananaName,"m_e_de_b%df%d",id_b,id_f);

            TKey *key = gDirectory->FindKey(bananaName);
            TH2D *h =  (TH2D*)gDirectory->Get(bananaName);
            h->Fill(eDet_energy[i],deDet_energy);

            // Fill LaBr3 energies and times
            for (int j=0; j<labr_mult ; ++j){
                h_LaBr3_energy->Fill(labr_energy[j],labr_ID[j]-1);
                h_LaBr3_time->Fill(labr_time[j],labr_ID[j]-1);
            }
        
            if(n%1000000==0){
                std::cout << ".";
                std::cout.flush();
            }
        } 
    }
    
    std::cout << endl;
    std::cout << "Total number of events:           " << nrow << endl;
    //std::cout << "Number of events inside the cut:  " << n_in_cut << endl;
    //std::cout << "Number of eDet background events: " << n_eDet_bg << endl;
    
    
    //TCanvas *c_de_e_energy = new TCanvas("c_de_e_energy","deDet eDet Energy",1);
    //h_de_e_energy->Draw("colz");
    
    
    //TH2D *h_Ex_gamma = new TH2D(*h_Ex_gamma_with_bg);
    //h_Ex_gamma->SetNameTitle("h_Ex_gamma","Excitation Energy vs gamma Energy");
    //if (!(h_Ex_gamma->GetSumw2N() > 0)) h_Ex_gamma->Sumw2(kTRUE);
    //h_Ex_gamma->Add(h_Ex_gamma_bg, -1.0);
    
    
    TFile *outputFile = new TFile("Mo96_plain.root","recreate");
    //outputFile->Write(); // Write all objects to file - this takes a looong time!!

    // Also possible to only write some objects
    h_eDet_mult->Write("h_eDet_mult",TObject::kOverwrite);
    h_LaBr3_mult->Write("h_LaBr3_mult",TObject::kOverwrite);
    h_deDet_energy->Write("h_deDet_energy",TObject::kOverwrite);
    h_eDet_energy->Write("h_eDet_energy",TObject::kOverwrite);
    h_eDet_time->Write("h_eDet_time",TObject::kOverwrite);
    h_LaBr3_energy->Write("h_LaBr3_energy",TObject::kOverwrite);
    h_LaBr3_time->Write("h_LaBr3_time",TObject::kOverwrite);
    // Write all the banana plots (Delta E-E) to file
    for(int i=0;i<64;++i){
     	deltaE_E_matrices[i]->Write();
    }


    //h_eDet_time->Write("h_eDet1_time_gs_19F_cut");
    //h_labr_energy->Write("h_labr_energy",TObject::kOverwrite);
    //h_labr_gamma_gamma->Write("h_labr_gamma_gamma_peak2",TObject::kOverwrite);
    //h_de_e_energy->Write("h_de_e_energy",TObject::kOverwrite);
    //h_eDet_energy_time->Write("h_eDet_energy_time_Ex_cut",TObject::kOverwrite);
    //h_labr_energy_time->Write("h_labr_energy_time",TObject::kOverwrite);
    //h_Ex_gamma_with_bg->Write("h_Ex_gamma_with_bg",TObject::kOverwrite);
    //h_Ex_gamma_bg->Write("h_Ex_gamma_bg",TObject::kOverwrite);
    //h_Ex_gamma->Write("h_Ex_gamma",TObject::kOverwrite);
    outputFile->Close();
    std::cout << "Histograms written to file  " << nrow << endl;

    
}
