#include "TChain.h"
#include "TString.h" 
#include <iostream> 
using namespace std;

#include "waveAna.h"
#include "waveAna.C"
#include "waveAnaHelper.h"
#include "waveAnaIncludes.h"

TString inFileFormat = "lysoData/data1/Lyso%d_s%d/Lyso%d__pos%dfull.root";
TString outFileFormat = "lysoData/data2/cry%d_sd%d_pos%d.root";

const int posNo = 5;
const int cryNo = 4;
const int sideNo = 1;
const int singleCrystal = 2; // -999

void goWaveAna() {

    gErrorIgnoreLevel = kFatal;
    
    if (singleCrystal < 0) {
        for (int cry = 0; cry < cryNo; cry++) {
            for (int side = 0; side < sideNo; side++) {
                for (int pos = 0; pos < posNo; pos++) {
                    inFileName = Form(inFileFormat, cry, side, cry, pos);
                    outFileName = Form(outFileFormat, cry, side, pos);
                    cout<<endl<<endl<<"       processing "<<inFileName<<" -> "<<outFileName<<endl; 
                    inFile = new TFile(inFileName);
                    inFile->GetObject(treeName, inTree);

                    waveAna *ANALYSIS = new waveAna(inTree);
                    ANALYSIS->Loop();
                }  
            }
        }
    } else {
        for (int side = 0; side < sideNo; side++) {
            for (int pos = 0; pos < posNo; pos++) {
                inFileName = Form(inFileFormat, singleCrystal, side, singleCrystal, pos);
                outFileName = Form(outFileFormat, singleCrystal, side, pos);
                cout<<endl<<endl<<"       processing "<<inFileName<<" -> "<<outFileName<<endl; 
                inFile = new TFile(inFileName);
                inFile->GetObject(treeName, inTree);

                waveAna *ANALYSIS = new waveAna(inTree);
                ANALYSIS->Loop();
            }  
        }
    }

}




