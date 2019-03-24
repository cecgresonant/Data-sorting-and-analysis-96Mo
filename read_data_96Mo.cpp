// Code from Frank to convert the Sirius data files into .root files
// Modified by Cecilie, 23 March 2019
// Using declarations_plain.h with "neutral" calibration coefficients 
// to determine peak positions for calibration of SiRi and OSCAR
// To compile: 
// > c++ read_data_96Mo.cpp XIA_CFD.cpp `root-config --libs --cflags` -o read_data
// Make also Delta E multiplicity to see how many strips are firing within the event
// In declarations_plain.h, updated with MAX_HITS_DEDET so we can keep track of the
// Delta E multiplicity

//#include "declarations.h"
#include "declarations_plain.h"
#include "XIA_CFD.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <string>
#include <vector>
#include <deque>
#include <algorithm>

#include <iostream>

bool ReadHit(FILE *file, word_t &hit)
{
	// First check if EOF has been reached.
	uint32_t hitdata[4];
	if (fread(hitdata, sizeof(uint32_t), 4, file) != 4) // EOF or error...
		return false;
    
    //Check for pile-up
    hit.finishcode = ( ( hitdata[0] & 0x80000000 ) > 0 ) ? 1 : 0;
    
    //Get the address
    hit.address = ( hitdata[0] & 0x00000FFF );
    
    // Parse the timestamp
    int64_t tmp1 = hitdata[1];
    int64_t tmp2 = (hitdata[2] & 0x0000FFFF);
    hit.timestamp = tmp2 << 32;
    hit.timestamp |= tmp1;
    
    // Extract CFD value.
    hit.cfddata = (hitdata[2] & 0xFFFF0000) >> 16;
    hit.cfdfail = 1;
    
    //Get timestamp in ns
    switch (pDetector[hit.address].sfreq) {
        case f250MHz :
            hit.cfdcorr = XIA_CFD_Fraction_250MHz(hit.cfddata, &hit.cfdfail);
            hit.timestamp *= 8;
            break;
        case f500MHz :
            hit.cfdcorr = XIA_CFD_Fraction_500MHz(hit.cfddata, &hit.cfdfail);
            hit.timestamp *= 10;
            break;
        case f000MHz :
            hit.cfdcorr = 0;
            hit.timestamp *= 0;
            break;
        default :
            hit.cfdcorr = 0;
            hit.timestamp *= 0;
            break;
    }
    
    // Extract energy
    hit.adcdata = (hitdata[3] & 0xFFFF);
    
    std::fseek(file, sizeof(uint32_t)*(-4), SEEK_CUR);
    
    return true; // Done!
}

void ReadFile(std::string filename)
{
    FILE *input = fopen(filename.c_str(), "rb");

    word_t hit, hit_before, hit_after;//, hit_deDet; // Make hit_deDet into a vector like the hit_eDet
    word_t hit_deDet[MAX_HITS_DEDET], hit_eDet[MAX_HITS_EDET], hit_labr[MAX_HITS_LABR];
    
    int64_t n_tot=0, n_events=0, n_pile_up=0, n_cfd_fail=0, temp_ts=0;
    uint16_t n_before=0, n_after=0, n_while=1;
    
    int tel_ID=0, deDet_mult=0, eDet_mult=0, labr_mult=0;
    int labr_ID[MAX_HITS_LABR], labr_ring[MAX_HITS_LABR], deDet_ID[MAX_HITS_DEDET];
    
    double eDet_time[MAX_HITS_EDET], labr_time[MAX_HITS_LABR], temp_cfd=0, eDet_energy[MAX_HITS_EDET], labr_energy[MAX_HITS_LABR]; 
    double deDet_energy[MAX_HITS_DEDET], deDet_time[MAX_HITS_DEDET];
    
    Long64_t deDet_timestamp=0;
    
    std::cout << "Results of file '" << filename << "':" << std::endl;
    
    std::string base_name = filename.substr(0, filename.find(".data"));
    base_name += ".root";
    TFile *f = new TFile(base_name.c_str(),"RECREATE");
    TTree *tree = new TTree("data","data");
    
    tree->Branch("deDet_ID", &deDet_ID, "deDet_ID/I");
    tree->Branch("tel_ID", &tel_ID, "tel_ID/I");
    tree->Branch("deDet_timestamp", &deDet_timestamp, "deDet_timestamp/L");
    tree->Branch("deDet_mult", &deDet_mult, "deDet_mult/I");
    tree->Branch("deDet_energy", deDet_energy, "deDet_energy[deDet_mult]/D");
    tree->Branch("deDet_time", deDet_time, "deDet_time[deDet_mult]/D");
    
    tree->Branch("eDet_mult", &eDet_mult, "eDet_mult/I");
    tree->Branch("eDet_energy", eDet_energy, "eDet_energy[eDet_mult]/D");
    tree->Branch("eDet_time", eDet_time, "eDet_time[eDet_mult]/D");
    
    tree->Branch("labr_mult", &labr_mult, "labr_mult/I");
    tree->Branch("labr_ID", labr_ID, "labr_ID[labr_mult]/I");
    tree->Branch("labr_ring", labr_ring, "labr_ring[labr_mult]/I");
    tree->Branch("labr_energy", labr_energy, "labr_energy[labr_mult]/D");
    tree->Branch("labr_time", labr_time, "labr_time[labr_mult]/D");
    
    while (ReadHit(input, hit)){
        
        deDet_mult=0;
        eDet_mult=0;
        labr_mult=0;
        
        // Initialize E detectors' vectors
        for (int i=0; i<MAX_HITS_EDET; ++i) {
            eDet_energy[i]=0;
            eDet_time[i]=0;
        }
        // Initialize LaBr3 detectors' vectors
        for (int i=0; i<MAX_HITS_LABR; ++i) {
            labr_ID[i]=0;
            labr_ring[i]=0;
            labr_energy[i]=0;
            labr_time[i]=0;
        }
        // Initialize Delta E detectors' vectors
        for (int i=0; i<MAX_HITS_DEDET; ++i) {
            deDet_energy[i]=0;
            deDet_time[i]=0;
        }        

        
        if (!hit.finishcode) {
            
            if (pDetector[hit.address].type==2) {
                //hit_deDet = hit;
                
                while (std::ftell(input)) {
                    
                    std::fseek(input, sizeof(uint32_t)*(-4), SEEK_CUR);
                    n_before++;
                    
                    ReadHit(input, hit_before);
                    
                    if (hit.timestamp - hit_before.timestamp < 1000) {
                        if (!hit_before.finishcode) {
                            if (pDetector[hit_before.address].type==3 && pDetector[hit_before.address].telNum==pDetector[hit.address].telNum) {
                                hit_eDet[eDet_mult]=hit_before;
                                eDet_mult++;
                            } else if (pDetector[hit_before.address].type==1) {
                                hit_labr[labr_mult]=hit_before;
                                labr_mult++;
                            }
                        }
                    } else {
                        break;
                    }
                    
                }
                
                std::fseek(input, sizeof(uint32_t)*(4*n_before), SEEK_CUR);
                n_before=0;
                
                while (n_while) {
                    
                    std::fseek(input, sizeof(uint32_t)*(4), SEEK_CUR);
                    n_after++;
                    
                    if (ReadHit(input, hit_after)) {
                        if (hit_after.timestamp - hit.timestamp < 1000) {
                            if (!hit_after.finishcode) {
                                if (pDetector[hit_after.address].type==3 && pDetector[hit_after.address].telNum==pDetector[hit.address].telNum) {
                                    hit_eDet[eDet_mult]=hit_after;
                                    eDet_mult++;
                                } else if (pDetector[hit_after.address].type==1) {
                                    hit_labr[labr_mult]=hit_after;
                                    labr_mult++;
                                }
                            }
                            
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                    
                }
                
                std::fseek(input, sizeof(uint32_t)*(-4*n_after), SEEK_CUR);
                n_after=0;
                
                if (eDet_mult>0) {
                    n_events++;
                    deDet_ID=pDetector[hit_deDet.address].detectorNum;
                    tel_ID=pDetector[hit_deDet.address].telNum;
                    deDet_energy=hit_deDet.adcdata+(gRandom->Uniform())-0.5;
                    deDet_energy=deDetCal[deDet_ID-1][0]+deDetCal[deDet_ID-1][1]*deDet_energy;
                    deDet_timestamp=hit_deDet.timestamp;
                    for (int i=0; i<eDet_mult ; ++i) {
                        eDet_energy[i]=hit_eDet[i].adcdata+(gRandom->Uniform())-0.5;
                        eDet_energy[i]=eDetCal[deDet_ID-1][0]+eDetCal[deDet_ID-1][1]*eDet_energy[i];
                        temp_ts=hit_eDet[i].timestamp-hit_deDet.timestamp;
                        temp_cfd=hit_eDet[i].cfdcorr-hit_deDet.cfdcorr;
                        eDet_time[i]=temp_ts+temp_cfd+eDetTimeShift[deDet_ID-1];
                        eDet_time[i]=eDet_time[i]+eDetTimeShiftRecal[deDet_ID-1];
                    }
                    for (int i=0; i<labr_mult ; ++i) {
                        labr_ID[i]=pDetector[hit_labr[i].address].detectorNum;
                        labr_ring[i]=pDetector[hit_labr[i].address].telNum-8;
                        labr_energy[i]=hit_labr[i].adcdata+(gRandom->Uniform())-0.5;
                        labr_energy[i]=labrCal[labr_ID[i]-1][0]+labrCal[labr_ID[i]-1][1]*labr_energy[i]+labrCal[labr_ID[i]-1][2]*labr_energy[i]*labr_energy[i];
                        temp_ts=hit_labr[i].timestamp-hit_deDet.timestamp;
                        temp_cfd=hit_labr[i].cfdcorr-hit_deDet.cfdcorr;
                        labr_time[i]=temp_ts+temp_cfd+labrTimeShift[deDet_ID-1][labr_ID[i]-1];
                    }
                    tree->Fill();
                }
            }
        }
        
        n_pile_up+=hit.finishcode;
        n_cfd_fail+=hit.cfdfail;
        n_tot++;
        if(n_tot%1000000==0){
            std::cout << ".";
            std::cout.flush();
        }
        std::fseek(input, sizeof(uint32_t)*(4), SEEK_CUR);
    }

    fclose(input);

    std::cout << std::endl;
    std::cout << "Total number of hits:        " << n_tot << std::endl;
    std::cout << "Number of correlated events: " << n_events << std::endl;
    std::cout << "Total number of pile-ups:    " << n_pile_up << std::endl;
    std::cout << "Total number of cfd fails:   " << n_cfd_fail << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    
    f->Write(0,TObject::kOverwrite);
    f->Close();
}

int main(int argc, char *argv[])
{

    std::vector<std::string> input_files;

    if (argc <= 1){
        std::cout << "Usage: >" << argv[0] << "<path/to/file1> <path/to/file2> ..." << std::endl;
    } else {
        for (int i = 1 ; i < argc ; ++i){
            input_files.push_back(argv[i]);
        }
    }

	for (size_t i = 0 ; i < input_files.size() ; ++i){
            ReadFile(input_files[i]);
    }

	return 0;
}
