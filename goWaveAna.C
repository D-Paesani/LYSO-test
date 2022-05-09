#include "TChain.h"
#include "TString.h" 
#include <iostream> 
using namespace std;

#include "waveAna.h"
#include "waveAna.C"
#include "waveAnaHelper.h"
#include "waveAnaIncludes.h"

TString inFileFormat = "lysoData/data1/Lyso%d_s%d_May2022/Lyso%d_s%d_May2022__pos%dfull.root";
TString outFileFormat = "lysoData/data2/cry%d_sd%d_pos%d.root";

const int posNo = 6;
const int cryNo = 4;
const int sideNo = 1;

const int singleCrystal = 2; // -999
const int singleSide = 1; // -999

void goWaveAna() {

    gErrorIgnoreLevel = kFatal;

    for (int cry = 0; cry < cryNo; cry++) {
        if (!(singleCrystal < 0)) { cry = singleCrystal; }

        for (int side = 0; side < sideNo; side++) {
            if (!(singleSide < 0)) { side = singleSide; }

            for (int pos = 0; pos < posNo; pos++) {
                inFileName = Form(inFileFormat, cry, side, cry, side, pos);
                outFileName = Form(outFileFormat, cry, side, pos);
                cout<<endl<<endl<<"       processing "<<inFileName<<" -> "<<outFileName<<endl; 
                inFile = new TFile(inFileName);
                inFile->GetObject(treeName, inTree);

                waveAna *ANALYSIS = new waveAna(inTree);
                ANALYSIS->Loop();
            }  
        
            if (!(singleSide<0)) { break; }
        }

        if (!(singleCrystal<0)) { break; }
    }
}




