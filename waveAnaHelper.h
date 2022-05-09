#pragma once 

#include "HistManager.h"

const Long64_t maxEvents = 50000000; 
const int useFixedGate = 1;

const int digiSamNo = 1024;
const int digiChNo = 3;
const int digiTime = 4;
const double wEx = 0.0;
const double wEy = 0.0;

const int blineSamFrom = (int)90/digiTime;
const int blineSamNo = (int)90/digiTime;
const int chargeAfterSamNo = (int)400/digiTime;
const int chargeBeforeSamNo = (int)80/digiTime;

const int peakSearchLimit = 800/digiTime;

const double timingBefore = 100;
const double timingAfter = 50;
const double timingCF = 0.10;

double chargeStart[digiChNo] = {200, 220, 220};
double chargeStop[digiChNo] = {900, 850, 850};
double baseStart = 50;
double baseStop = 150;

TString inFileName, outFileName, runName;
TFile *outFile, *inFile;
TDirectory *samples_dir;
Long64_t etp;
Long64_t evNum;
HistManager *HM;

TTree *inTree, *cryTree;
TString treeNameOut = "crytree", treeName = "Wavefull";



