#ifndef DECLARATIONS_H
#define DECLARATIONS_H

// From Frank's declarations.h file, modified by Cecilie 23 March 2019
// Setting all gains and shifts to initial values
// This file is read in read_data.cpp


#define NUM_LABR_DETECTORS 30   //!< Number of LaBr detectors
#define NUM_SI_DE_DET 64        //!< Number of Si dE sections
#define NUM_SI_E_DET 8          //!< Number of E Si rings
#define NUM_SI_E_GUARD 4        //!< Number of E guard rings
#define NUM_PPAC 4              //!< Number of PPACs

#define NUM_SI_RINGS 8          //!< Number of Si rings
#define TOTAL_NUMBER_OF_MODULES 9   //!< Number of modules
#define TOTAL_NUMBER_OF_ADDRESSES 144   //! Total number of address that needs to be defined

#define MAX_HITS_EDET 8    //!< Maximum number of hits in E detectors per event
#define MAX_HITS_DEDET 128    //!< Maximum number of hits in dE detectors per event - added by Cecilie
#define MAX_HITS_LABR 64    //!< Maximum number of hits in LaBr3 detectors per event



#include <stdint.h>

enum DetectorType {
    invalid,    //!< Invalid address: type==0
    labr,       //!< Is a labr detector: type==1
    deDet,      //!< Is a Delta-E segment: type==2
    eDet,       //!< Is a E detector: type==3
    eGuard,     //!< Is a E guard ring: type==4
    ppac,       //!< Is a PPAC: type==5
    rfchan,     //!< Is a RF channel: type==6
    unused      //!< Is a unused XIA channel: type==7
};

enum ADCSamplingFreq {
    f250MHz,    //!< 250 MHz sampling frequency
    f500MHz,    //!< 500 MHz sampling frequency
    f000MHz     //!< If invalid address
};

typedef struct {
    uint16_t address;           //!< ADC address of the detector
    enum ADCSamplingFreq sfreq; //!< ADC sampling frequency
    enum DetectorType type;     //!< Type of detector
    int detectorNum;            //!< 'Linear' number of the detector
    int telNum;                 //!< Telescope number (ie. E back detector for the dE front detector).For the labr it means the ring number
} DetectorInfo_t;

typedef struct {
    uint16_t address;        //!< Holds the address of the ADC.
    uint16_t adcdata;        //!< Data read out from the ADC.
    uint16_t cfddata;       //!< Fractional difference of before/after zero-crossing.
    char finishcode;        //!< Pile-up flag.
    char cfdfail;           //!< Flag to tell if the CFD was forced or not.
    int64_t timestamp;        //!< Timestamp in [ns].
    double cfdcorr;         //!< Correction from the CFD.
} word_t;

static DetectorInfo_t pDetector[144] =
{
    {0, f000MHz, unused, 0, 0},
    {1, f000MHz, unused, 0, 0},
    {2, f000MHz, unused, 0, 0},
    {3, f000MHz, unused, 0, 0},
    {4, f000MHz, unused, 0, 0},
    {5, f000MHz, unused, 0, 0},
    {6, f000MHz, unused, 0, 0},
    {7, f000MHz, unused, 0, 0},
    {8, f000MHz, unused, 0, 0},
    {9, f000MHz, unused, 0, 0},
    {10, f000MHz, unused, 0, 0},
    {11, f000MHz, unused, 0, 0},
    {12, f000MHz, unused, 0, 0},
    {13, f000MHz, unused, 0, 0},
    {14, f000MHz, unused, 0, 0},
    {15, f000MHz, unused, 0, 0},
    {16, f000MHz, unused, 0, 0},
    {17, f000MHz, unused, 0, 0},
    {18, f000MHz, unused, 0, 0},
    {19, f000MHz, unused, 0, 0},
    {20, f000MHz, unused, 0, 0},
    {21, f000MHz, unused, 0, 0},
    {22, f000MHz, unused, 0, 0},
    {23, f000MHz, unused, 0, 0},
    {24, f000MHz, unused, 0, 0},
    {25, f000MHz, unused, 0, 0},
    {26, f000MHz, unused, 0, 0},
    {27, f000MHz, unused, 0, 0},
    {28, f000MHz, unused, 0, 0},
    {29, f000MHz, unused, 0, 0},
    {30, f000MHz, unused, 0, 0},
    {31, f000MHz, unused, 0, 0},
    {32, f500MHz, labr, 1, 9},//ring 1, theta = 37 deg
    {33, f500MHz, labr, 2, 9},
    {34, f500MHz, labr, 3, 9},
    {35, f500MHz, labr, 4, 9},
    {36, f500MHz, labr, 5, 9},
    {37, f500MHz, labr, 6, 10},//ring 2, theta = 63 deg
    {38, f500MHz, labr, 7, 11},//ring 3, theta = 79 deg
    {39, f500MHz, labr, 8, 10},
    {40, f500MHz, labr, 9, 11},
    {41, f500MHz, labr, 10, 10},
    {42, f500MHz, labr, 11, 11},
    {43, f500MHz, labr, 12, 10},
    {44, f500MHz, labr, 13, 11},
    {45, f500MHz, labr, 14, 10},
    {46, f500MHz, labr, 15, 11},
    {47, f500MHz, unused, 0, 0},
    {48, f500MHz, labr, 16, 12},//ring 4, theta = 101 deg
    {49, f500MHz, labr, 17, 13},//ring 5, theta = 117 deg
    {50, f500MHz, labr, 18, 13},
    {51, f500MHz, labr, 19, 13},
    {52, f500MHz, labr, 20, 12},
    {53, f500MHz, labr, 21, 13},
    {54, f500MHz, labr, 22, 12},
    {55, f500MHz, labr, 23, 13},
    {56, f500MHz, labr, 24, 12},
    {57, f500MHz, labr, 25, 13},
    {58, f500MHz, labr, 26, 13},
    {59, f500MHz, labr, 27, 14},//ring 5, theta = 143 deg
    {60, f500MHz, labr, 28, 14},
    {61, f500MHz, labr, 29, 14},
    {62, f500MHz, labr, 30, 12}, //This one is bad, jumping around in energy
    {63, f500MHz, unused, 0, 0},
    {64, f250MHz, deDet, 1, 1},
    {65, f250MHz, deDet, 2, 1},
    {66, f250MHz, deDet, 3, 1},
    {67, f250MHz, deDet, 4, 1},
    {68, f250MHz, deDet, 5, 1},
    {69, f250MHz, deDet, 6, 1},
    {70, f250MHz, deDet, 7, 1},
    {71, f250MHz, deDet, 8, 1},
    {72, f250MHz, deDet, 9, 2},
    {73, f250MHz, deDet, 10, 2},
    {74, f250MHz, deDet, 11, 2},
    {75, f250MHz, deDet, 12, 2},
    {76, f250MHz, deDet, 13, 2},
    {77, f250MHz, deDet, 14, 2},
    {78, f250MHz, deDet, 15, 2},
    {79, f250MHz, deDet, 16, 2},
    {80, f250MHz, deDet, 17, 3},
    {81, f250MHz, deDet, 18, 3},
    {82, f250MHz, deDet, 19, 3},
    {83, f250MHz, deDet, 20, 3},
    {84, f250MHz, deDet, 21, 3},
    {85, f250MHz, deDet, 22, 3},
    {86, f250MHz, deDet, 23, 3},
    {87, f250MHz, deDet, 24, 3},
    {88, f250MHz, deDet, 25, 4},
    {89, f250MHz, deDet, 26, 4},
    {90, f250MHz, deDet, 27, 4},
    {91, f250MHz, deDet, 28, 4},
    {92, f250MHz, deDet, 29, 4},
    {93, f250MHz, deDet, 30, 4},
    {94, f250MHz, deDet, 31, 4},
    {95, f250MHz, deDet, 32, 4},
    {96, f250MHz, deDet, 33, 5},
    {97, f250MHz, deDet, 34, 5},
    {98, f250MHz, deDet, 35, 5},
    {99, f250MHz, deDet, 36, 5},
    {100, f250MHz, deDet, 37, 5},
    {101, f250MHz, deDet, 38, 5},
    {102, f250MHz, deDet, 39, 5},
    {103, f250MHz, deDet, 40, 5},
    {104, f250MHz, deDet, 41, 6},
    {105, f250MHz, deDet, 42, 6},
    {106, f250MHz, deDet, 43, 6},
    {107, f250MHz, deDet, 44, 6},
    {108, f250MHz, deDet, 45, 6},
    {109, f250MHz, deDet, 46, 6},
    {110, f250MHz, deDet, 47, 6},
    {111, f250MHz, deDet, 48, 6},
    {112, f250MHz, deDet, 49, 7},
    {113, f250MHz, deDet, 50, 7},
    {114, f250MHz, deDet, 51, 7},
    {115, f250MHz, deDet, 52, 7},
    {116, f250MHz, deDet, 53, 7},
    {117, f250MHz, deDet, 54, 7},
    {118, f250MHz, deDet, 55, 7},
    {119, f250MHz, deDet, 56, 7},
    {120, f250MHz, deDet, 57, 8},
    {121, f250MHz, deDet, 58, 8},
    {122, f250MHz, deDet, 59, 8},
    {123, f250MHz, deDet, 60, 8},
    {124, f250MHz, deDet, 61, 8},
    {125, f250MHz, deDet, 62, 8},
    {126, f250MHz, deDet, 63, 8},
    {127, f250MHz, deDet, 64, 8},
    {128, f250MHz, eGuard, 1, 1},
    {129, f250MHz, eDet, 1, 1},
    {130, f250MHz, eGuard, 2, 2},
    {131, f250MHz, eDet, 2, 2},
    {132, f250MHz, eGuard, 3, 3},
    {133, f250MHz, eDet, 3, 3},
    {134, f250MHz, eGuard, 4, 4},
    {135, f250MHz, eDet, 4, 4},
    {136, f250MHz, eGuard, 5, 5},
    {137, f250MHz, eDet, 5, 5},
    {138, f250MHz, eGuard, 6, 6},
    {139, f250MHz, eDet, 6, 6},
    {140, f250MHz, eGuard, 7, 7},
    {141, f250MHz, eDet, 7, 7},
    {142, f250MHz, eGuard, 8, 8},
    {143, f250MHz, eDet, 8, 8},
};

double labrCal [30][3] = {
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.},
    {0.,1.000000,0.} // This one is bad in the 96Mo exp.
};


double labrTimeShift [64][30] = {
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
};


double eDetCal [64][2] = {
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//8
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//16
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//24
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//32
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//40
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//48 
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//56
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000}//64 
};

double eDetTimeShift [64] = {
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//8
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//16 
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//24
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//32 
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//40
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//48 
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//56
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.//64 
};

double eDetTimeShiftRecal [64] = {
    0.,
    0.,
    0.,//
    0.,//
    0.,//
    0.,//
    0.,
    0.,//8
    0.,//
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//16 
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//24
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//32
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//40
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//48 
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,//56
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
    0.//64 
};

double deDetCal [64][2] = {
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//8
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//16 
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//24
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//32 
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//40
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//48 
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},//56
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000},
    {0.,1.000000}//64 
};

//! Get detector method
/*! \return Detector structure containing information about the
 *  detector at address.
 */
//DetectorInfo_t GetDetector(uint16_t address   /*!< Address of the detector to get */);

//! Get sampling frequency
/*! \return The XIA module sampling frequency
 */
//enum ADCSamplingFreq GetSamplingFrequency(uint16_t address    /*!< ADC address    */);

#endif // DECLARATIONS_H
