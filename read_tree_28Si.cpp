// ************************************************************************//
// Frank's nice sorting code :)
// Modified by Cecilie, March 2019, for the 96Mo data set
// This code is for the 28Si calibration run
// LAST UPDATE: 8 April 2019
// 
// Info on the experiment: 28Si(p,p'gamma) with E_p = 16 MeV, I_p = 3.6-3.7 nA
// SiRi in backward angles, 126-140 degrees
// OSCAR in the close configuration, 16.3 cm - 17 cm from the target center position
// One OSCAR detector (no. 30) was not working properly, "jumping" in the energy.
// The target was 28Si, thickness is not known but around 2-3 mg/cm^2 and ≈ 92% in 28Si (natural Si).
// The validation signal for the XiA boxes were the E detectors, the time window is 2 us
// In the sorting code, we are looking 1 us forward and 1 us backward in time, 
// so in total the time window is 2 us here also.
//
// Data files, 96Mo: 
// sirius-20190314-084105.data and sirius-20190314-084105-big00X.data
// sirius-20190314-145712.data and sirius-20190314-145712-big-000.data
// 
// Data files, 28Si (for calibration):
// sirius-20190313-142346.data and sirius-20190313-142346-big-00X.data  
// 
// Elog: Nd-p-2019
// CHECK Delta E strips 27-30, and 57, to see if there are some "shadow" 
// from neighboring strips
// ************************************************************************//


#include <iostream>
#include <stdio.h>
#include <string>
#define MAX_HITS_EDET 8    //!< Maximum number of hits in e detectors per event
#define MAX_HITS_LABR 64    //!< Maximum number of hits in labr detectors per event (can be more than 30 due to the long time range for the vaidation signal)

void read_tree_28Si(){
    
    Long64_t ns_to_hours = 3600000000000;
    Long64_t ns_to_minutes = 60000000000;
    Long64_t ns_to_seconds =  1000000000;
    
    TChain *t = (TChain *) new TChain("data");
    
    // Character string for the file name
    char name[1024];
    
    // Set file name, open the file and add it to the TChain t.
    // 96Mo
    /*sprintf(name, "sirius-20190314-084105.root");
    ifstream rootfile1_96Mo(name);
    if (rootfile1_96Mo){
        t->Add(name);
        std::cout << "Added file " << name << " to chain." << std::endl;
    } else {
        std::cout << "File not read" << std::endl;
    }*/

    // 28Si
    sprintf(name, "sirius-20190313-142346.root");
    ifstream rootfile1_28Si(name);
    if (rootfile1_28Si){
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
    Double_t Ex_coeff [8][3] = {//from qkinz for 28Si(p,p') assuming 4 mg/cm2
        {14533.325,-1.030173,-12.39e-07},
        {14532.849,-1.027640,-13.02e-07},
        {14532.560,-1.025162,-13.57e-07},
        {14529.979,-1.022221,-14.32e-07},
        {14527.384,-1.019329,-15.00e-07},
        {14523.407,-1.016213,-15.75e-07},
        {14517.837,-1.012856,-16.57e-07},
        {14510.415,-1.009234,-17.48e-07}
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
    const int max_e_labr3 =  40000; // Max energy for the LaBr3 detectors
    
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
    TH2D *deltaE_E_matrices[64] = {NULL};

    int histo_no=0;
    for(int b=0; b<8; ++b ) {
        for(int f=0; f<8; ++f ) {
        	ostringstream histogramNameStream;
            histogramNameStream << "m_e_de_b" << b << "f" << f;
        	//deltaE_E_matrices[histo_no] = new TH2D(histogramNameStream.str().c_str(),histogramNameStream.str().c_str(),500,0,max_e,100,0,max_de); // This is best for the peaks2D.C script
        	deltaE_E_matrices[histo_no] = new TH2D(histogramNameStream.str().c_str(),histogramNameStream.str().c_str(),1000,0,max_e,500,0,max_de);
            deltaE_E_matrices[histo_no]->GetXaxis()->SetTitle("E_{back} (keV)");
        	deltaE_E_matrices[histo_no]->GetYaxis()->SetTitle("#Delta E_{front} (keV)");
        	++histo_no;
        }
        //cout << histogramNameStream.str().c_str() << endl;
    }

    // Cecilie, modification 8 April 2019:
    // Making Delta E + E histograms for all 64 combinations of strips and back detectors
    // to check alignment.  
    TH1D *h_e_de_spectra[64] = {NULL};
    histo_no=0;
    for(int b=0; b<8; ++b ) {
        for(int f=0; f<8; ++f ) {
            ostringstream histogramNameStream;
            histogramNameStream << "h_e_de_b" << b << "f" << f;
            h_e_de_spectra[histo_no] = new TH1D(histogramNameStream.str().c_str(),histogramNameStream.str().c_str(),1000,0,max_e);
            h_e_de_spectra[histo_no]->GetXaxis()->SetTitle("E_{back}+#Delta E_{front} (keV)");
            h_e_de_spectra[histo_no]->GetYaxis()->SetTitle("counts");
            ++histo_no;
        }
    }

    // Cecilie, modification 8 April 2019:
    // Making h_Ex_fX histograms for all 8 rings (126-140 degrees)
    // to check alignment. Also make the sum of all, h_ex
    TH1D *h_ex_ring_spectra[8] = {NULL};
    histo_no=0;
    for(int f=0; f<8; ++f ) {
        ostringstream histogramNameStream;
        histogramNameStream << "h_ex_f" << f;
        h_ex_ring_spectra[histo_no] = new TH1D(histogramNameStream.str().c_str(),histogramNameStream.str().c_str(),1000,-1000.,max_e);
        h_ex_ring_spectra[histo_no]->GetXaxis()->SetTitle("E_{x} (keV)");
        h_ex_ring_spectra[histo_no]->GetYaxis()->SetTitle("counts");
        //std::cout << histogramNameStream.str().c_str() << endl;
        ++histo_no;
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

            // Fill the E+Delta E histograms
            char *h_e_de_Name = new char[10];
            sprintf(h_e_de_Name,"h_e_de_b%df%d",id_b,id_f);
            TKey *key2 = gDirectory->FindKey(h_e_de_Name);
            TH1D *h_2 =  (TH1D*)gDirectory->Get(h_e_de_Name);
            Double_t e_sum = eDet_energy[i]+deDet_energy;
            h_2->Fill(e_sum);

            // Fill the h_ex_f0...h_ex_f7 histograms, 
            // this is excitation energy for each angle 126-140 degrees
            char *h_ex_ring_Name = new char[10];
            sprintf(h_ex_ring_Name,"h_ex_f%d",id_f);
            TKey *key3 = gDirectory->FindKey(h_ex_ring_Name);
            TH1D *h3 =  (TH1D*)gDirectory->Get(h_ex_ring_Name);


            Double_t excitation_energy = Ex_coeff[(deDet_ID-tel_ID)%8][0]+(e_sum)*Ex_coeff[(deDet_ID-tel_ID)%8][1]+(e_sum)*(e_sum)*Ex_coeff[(deDet_ID-tel_ID)%8][2];
            //std::cout << " Ex =  " << excitation_energy << " for histogram " << h_ex_ring_Name << std::endl;
            if(excitation_energy > -1000. && excitation_energy < 20000.){
                h3->Fill(excitation_energy);
            }
            //std::cout << " DUCK DUCK " << std::endl;

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

    // When using declarations_plain.h, i.e. default calibration coefficients:
    TFile *outputFile = new TFile("Si28_plain.root","recreate");
    // When using declarations_sirical.h, i.e. SiRi calibration coefficients:
    //TFile *outputFile = new TFile("Si28_sirical.root","recreate");
    //outputFile->Write(); // Write all objects to file - this takes a looong time!!
    // Also possible to only write some objects, this goes faster
    h_eDet_mult->Write("h_eDet_mult",TObject::kOverwrite);
    h_LaBr3_mult->Write("h_LaBr3_mult",TObject::kOverwrite);
    h_deDet_energy->Write("h_deDet_energy",TObject::kOverwrite);
    h_eDet_energy->Write("h_eDet_energy",TObject::kOverwrite);
    h_eDet_time->Write("h_eDet_time",TObject::kOverwrite);
    h_LaBr3_energy->Write("h_LaBr3_energy",TObject::kOverwrite);
    h_LaBr3_time->Write("h_LaBr3_time",TObject::kOverwrite);
    // Write all the banana plots (Delta E-E) to file
    // and the Delta E + E histograms
    for(int i=0;i<64;++i){
     	deltaE_E_matrices[i]->Write();
        h_e_de_spectra[i]->Write();
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
    std::cout << "Histograms written to file. " << endl;

    
}
