#define wfAna_cxx

#include <fstream>
#include <chrono>
#include <TLine.h>
#include <iostream>
#include <list>
#include <algorithm>
#include <unistd.h>
#include <stdio.h>

#include "TApplication.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphSmooth.h"
#include "TSpline.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TButton.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TRandom.h"
#include "TMath.h"
#include "math.h"
#include "TH2.h"
#include "TStyle.h"

#include "wfAna.h"
#include "HistManager.h"
#include "pars1.h"

using namespace std;

TFile* outFile;
TTree* cryTree;
TDirectory *samples_dir;
Long64_t etp;
HistManager HM;
Long64_t evNum;
Double_t wCharge[digiChNo], wTime[digiChNo], wTimePeak[digiChNo], wPeak[digiChNo], wBline[digiChNo], wBrms[digiChNo], waveTime[digiSamNo], wPed[digiSamNo];
Double_t wave[digiChNo][digiSamNo];

void wfAna::Begin(TTree *tree)
{
   cout<<"::::::::::::::::::::: Start analysis :::::::::::::::::::::"<<endl<<endl;

   TString opti = GetOption();
   cout<<"Options: "<<opti<<endl;
   TObjArray* opt = opti.Tokenize(" ");
   TString runName = ((TObjString*)(opt->At(0)))->String(); cout<<"RunName: "<<runName<<endl;
   TString inName = ((TObjString*)(opt->At(1)))->String(); cout<<"InFile: "<<inName<<endl;
   TString outName = ((TObjString*)(opt->At(2)))->String(); cout<<"OutFile: "<<outName<<endl<<endl;
   outFile = new TFile(outName.Data(), "RECREATE"); 
   outFile->cd();

   Init(tree);

   Long64_t en = fReader.GetEntries();
   etp = min(maxEvents, en);
   cout << "Number of events to process: " << etp << endl << endl;

   outFile->cd();
   cryTree = new TTree(treeNameOut.Data(),treeNameOut.Data());          
   cryTree->SetAutoSave(1000);
   cryTree->Branch("evNum", &evNum);
   cryTree->Branch("wave", &wave, "wave[3][1024]/D");
   cryTree->Branch("waveTime", &waveTime, "waveTime[1024]/D");
   cryTree->Branch("wCharge", &wCharge, "wCharge[3]/D");
   cryTree->Branch("wTime",   &wTime,  "wTime[3]]/D"); 
   cryTree->Branch("wBline",  &wBline,  "wBline[3]/D"); 
   cryTree->Branch("wBrms",   &wBrms,  "wBrms[3]/D"); 
   cryTree->Branch("wPed",   &wPed,  "wPed[3]/D"); 
   cryTree->Branch("wTimePeak", &wTimePeak,  "wTimePeak[3]/D"); 
   cryTree->Branch("wPeak", &wPeak,  "wPeak[3]/D"); 

   samples_dir = outFile->mkdir("randomSpecimens");
   gRandom->SetSeed();
   outFile->cd();

}

Bool_t wfAna::Process(Long64_t entry)
{
   
   double toss = gRandom->Uniform(0, etp/50);

   if (entry > etp) {return kFALSE;}
   if (!(entry%1000)) {cout << Form( "    processing evt %lld / %lld  ( %.0f%% )", entry, etp, (float)(100*entry/etp) ) << endl;} 
   fReader.SetLocalEntry(entry);

   evNum = *evnum;

   for (int iCh = 0; iCh < digiChNo; iCh++) {

      TGraphErrors waveGra = TGraphErrors(digiSamNo);
      double blineSum = 0, blineSum2 = 0, blineTmp = 0, chargeTmp = 0, peakTmp = 0, peakTimeTmp = 0, timeTmp = 0, thr = 0, pedTmp = 0;
      double valTmp = -99999, maxYTmp = -99999;
      int maxXTmp = 0;

      for (int iSam = 0; iSam < digiSamNo; iSam++) {

         if (iCh == 0) {
            valTmp = onda0[iSam];
         } else if (iCh == 1) {
            valTmp = onda1[iSam];
         } else if (iCh == 2) {
            valTmp = onda2[iSam];
         }
         valTmp =  valTmp/0.98;
         if (valTmp > maxYTmp && iSam < peakSearchLimit) { maxYTmp = valTmp; maxXTmp = iSam; }
         wave[iCh][iSam] = valTmp;
         waveTime[iSam] = time[iSam];
      }
      peakTimeTmp = (double)digiTime*maxXTmp;

      int tLow = max((int)0, maxXTmp - blineSamFrom - blineSamNo), tHi = max((int)0, maxXTmp - blineSamFrom);
      double deltat = (double)(tHi-tLow);
      for (int iSam = tLow; iSam < tHi; iSam++) {
         valTmp = wave[iCh][iSam];
         blineSum += valTmp; 
         blineSum2 += valTmp*valTmp;
      }
      pedTmp = blineTmp;
      blineTmp = blineSum/deltat;
      wBline[iCh] = blineTmp;
      wBrms[iCh] = TMath::Sqrt(TMath::Abs(blineSum2/deltat - blineTmp*blineTmp));
      peakTmp = maxYTmp - blineTmp;
      pedTmp = (double)digiTime*(pedTmp - blineTmp*deltat)/50; //non si calcola così !!!!
      wPed[iCh] = pedTmp;

      for (int iSam = 0; iSam < digiSamNo; iSam++) {
         wave[iCh][iSam] += -blineTmp;
         waveGra.SetPoint(iSam, waveTime[iSam], wave[iCh][iSam]);
         waveGra.SetPointError(iSam, wEx, wEy);
      }

      int tLow2 = max((int)0, maxXTmp - chargeBeforeSamNo), tHi2 = max((int)0, maxXTmp + chargeAfterSamNo);
      for (int iSam = tLow2; iSam < tHi2; iSam++) {
         chargeTmp += waveGra.GetPointY(iSam); // - blineTmp; //max((double)waveGra.GetPointY(iSam) - blineTmp, (double)0);
      }
      chargeTmp = (double)digiTime*chargeTmp/50;
      wCharge[iCh] = chargeTmp;
   
      double tmin = max((double)0, peakTimeTmp - timingBefore), tmax = max((double)0, peakTimeTmp + timingAfter);
      TSpline5 waveSp = TSpline5("wsp", &waveGra); 
      auto waveSpFun = [&waveSp](double *x, double *){ return waveSp.Eval(x[0]); }; 
      TF1 waveFitFun = TF1("fitf", waveSpFun, tmin, tmax, 0); 
      peakTmp = waveFitFun.GetMaximum(tmin, tmax);
      peakTimeTmp = waveFitFun.GetMaximumX(tmin, tmax);      
      thr = peakTmp*timingCF;
      timeTmp = waveFitFun.GetX(thr);   
      wTime[iCh] = timeTmp;
      wPeak[iCh] = peakTmp;
      wTimePeak[iCh] = peakTimeTmp;

      if ( toss < 1 ) {
         samples_dir->cd();
         TCanvas cc = TCanvas(Form("e%lld_ch%d", entry, iCh)); cc.cd(); 
         waveGra.SetTitle(Form("bl=%.3f br=%.3f pk=%.2f pt=%.2f rt=%.2f qq=%.2f", blineTmp, wBrms[iCh], peakTmp, peakTimeTmp, timeTmp, chargeTmp));
         waveGra.SetLineWidth(0); waveGra.SetMarkerStyle(20); waveGra.SetMarkerSize(.2); waveGra.SetMarkerColor(kBlue); waveGra.Draw(""); 
         waveFitFun.SetLineColor(kTeal); waveFitFun.Draw("same");
         waveSp.SetLineColor(kBlack); waveSp.Draw("same");
         TMarker tp = TMarker(timeTmp, waveSp.Eval(timeTmp), 2); tp.SetMarkerSize(3); tp.SetMarkerColor(kRed); tp.Draw("same"); 
         TLine t1 = TLine(tHi*digiTime, 0, tHi*digiTime, peakTmp);  t1.SetLineColor(kPink); t1.Draw("same");
         TLine t2 = TLine(tLow*digiTime, 0, tLow*digiTime, peakTmp);  t2.SetLineColor(kPink); t2.Draw("same");
         TLine t3 = TLine(tHi2*digiTime, 0, tHi2*digiTime, peakTmp);  t3.SetLineColor(kPink); t3.Draw("same");
         TLine t4 = TLine(0, peakTmp, peakTimeTmp, peakTmp);  t4.SetLineColor(kPink); t4.Draw("same");
         TLine t5 = TLine(tmin, 0, tmin, peakTmp);  t5.SetLineColor(kGreen); t5.Draw("same");
         TLine t6 = TLine(tmax, 0, tmax, peakTmp);  t6.SetLineColor(kGreen); t6.Draw("same");

         cc.Write(); 
         outFile->cd();
      }

   }
   
   cryTree->Fill();
   return kTRUE;
}

void wfAna::Terminate() {

   cout<<endl<<endl;

   HM.ProcessBoxes(); 

   cryTree->Write();
   outFile->Close();

   cout<<endl<<endl<<"::::::::::::::::::::: done :::::::::::::::::::::"<<endl<<endl;

}

void wfAna::SlaveBegin(TTree * /*tree*/) {}
void wfAna::SlaveTerminate() {}




//#undef wfAna_cxx
