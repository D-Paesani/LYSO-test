#include "TChain.h"
#include "TString.h" 
#include "wfAna.h"
#include <iostream> 
using namespace std;


TString inFileFormat = "data/data1/lyso_cr%d_sd%d/lyso_cr%d_sd%d_pos%dfull.root";
TString outFileFormat = "data/data2/cry%d_sd%d_pos%d.root";
TString treeName = "Wavefull";
TString analysisName = "wfAna.C+";
TString runName = "lyso";

const int posNo = 3;
const int cryNo = 4;

void goAna() {
    for (int cry = 0; cry < cryNo; cry++) {
        for (int side = 0; side < 2; side++) {
            for (int pos = 0; pos < posNo; pos++) {
                TString inName = Form(inFileFormat, cry, side, cry, side, pos);
                TString outName = Form(outFileFormat, cry, side, pos);
                cout<<"-------------------------------> processing "<<inName<<" -> "<<outName<<endl; 
                //continue;
                TChain chain(treeName.Data());  
                chain.Add(inName); 
                chain.Process(analysisName, runName + " " + inName + " " + outName + " opt3 opt4");
            }  
        }
    }
}
