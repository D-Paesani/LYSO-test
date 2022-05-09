#define waveAna_cxx
#include "waveAna.h"

#include "waveAnaHelper.h"

void waveAna::Loop() {
   if (fChain == 0) return;

   //Begin
      outFile = new TFile(outFileName.Data(), "RECREATE"); 
      outFile->cd();
      samples_dir = outFile->mkdir("randomSpecimens");

      Double_t wCharge[digiChNo], wTime[digiChNo], wTimePeak[digiChNo], wPeak[digiChNo], wBline[digiChNo], wBrms[digiChNo], wPed[digiSamNo];
      Double_t wave[digiChNo][digiSamNo], waveTime[digiSamNo];

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

      outFile->cd();

      HM = new HistManager();
      HM->SetOutFile(outFile);

      gRandom->SetSeed();
   //Begin

   //Loop
      Long64_t en = fChain->GetEntriesFast();
      etp = min(maxEvents, en);
      cout << "Number of events to process: " << etp << endl << endl;
      Long64_t nbytes{0}, nb{0};
      for (Long64_t jentry=0; jentry<etp;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry); 
         nbytes += nb;
   
         if (!(jentry%1000)) {cout << Form( "    processing evt %lld / %lld  ( %.0f%% )", jentry, etp, (float)(100*jentry/etp) ) << endl;} 
         evNum = jentry;

         for (int iCh = 0; iCh < digiChNo; iCh++) {

            TGraphErrors waveGra = TGraphErrors(digiSamNo);
            double blineSum{0}, blineSum2{0}, blineTmp{0}, pedTmp{0};
            double chargeTmp{0}, peakTmp{0}, peakTimeTmp{0}, timeTmp{0}, thr{0};
            double valTmp = -99999, tTmp;
            int samMax{0}, baseN{0}, chargeN{0};

            for (int iSam = 0; iSam < digiSamNo; iSam++) {

               if (iCh == 0) {
                  valTmp = onda0[iSam];
               } else if (iCh == 1) {
                  valTmp = onda1[iSam];
               } else if (iCh == 2) {
                  valTmp = onda2[iSam];
               }
               valTmp =  valTmp/0.98; 
               wave[iCh][iSam] = valTmp;
               waveTime[iSam] = time[iSam];
               tTmp = time[iSam];

               waveGra.SetPoint(iSam, tTmp, valTmp);
               waveGra.SetPointError(iSam, wEx, wEy);

               if (valTmp > peakTmp && iSam < peakSearchLimit) { 
                  peakTmp = valTmp; 
                  samMax = iSam;
                  peakTimeTmp = tTmp;
               }

               if (tTmp > baseStart && tTmp < baseStop) {  
                  blineSum += valTmp; 
                  blineSum2 += valTmp*valTmp; 
                  baseN++;
               }

               if (tTmp > chargeStart[iCh] && tTmp < chargeStop[iCh]) {
                  chargeTmp += valTmp; 
                  chargeN++;
               }

            }

            blineTmp = blineSum/(double)baseN;
            wBline[iCh] = blineTmp;
            waveGra.MovePoints(-blineTmp, 0);
            wBrms[iCh] = TMath::Sqrt(blineSum2/(double)baseN - blineTmp*blineTmp);
            peakTmp = peakTmp - blineTmp;
            pedTmp = -9999;
            wPed[iCh] = pedTmp;
            chargeTmp = chargeTmp - (double)chargeN*blineTmp;
            chargeTmp = chargeTmp*(double)digiTime/50;
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

            double toss = gRandom->Uniform(0, etp/50);
            if ( toss < 1 ) {
               samples_dir->cd();
               TCanvas cc = TCanvas(Form("e%lld_ch%d", jentry, iCh)); cc.cd(); 
               waveGra.SetTitle(Form("bl=%.3f br=%.3f pk=%.2f pt=%.2f rt=%.2f qq=%.2f", blineTmp, wBrms[iCh], peakTmp, peakTimeTmp, timeTmp, chargeTmp));
               waveGra.SetLineWidth(0); waveGra.SetMarkerStyle(20); waveGra.SetMarkerSize(.2); waveGra.SetMarkerColor(kBlue); waveGra.Draw(""); 
               waveFitFun.SetLineColor(kTeal); waveFitFun.Draw("same");
               waveSp.SetLineColor(kBlack); waveSp.Draw("same");
               TMarker tp = TMarker(timeTmp, waveSp.Eval(timeTmp), 2); tp.SetMarkerSize(3); tp.SetMarkerColor(kRed); tp.Draw("same"); 
               TLine t0 = TLine(baseStart*digiTime, 0, baseStart*digiTime, peakTmp);  t0.SetLineColor(kPink); t0.Draw("same");
               TLine t1 = TLine(baseStop*digiTime, 0, baseStop*digiTime, peakTmp);  t1.SetLineColor(kPink); t1.Draw("same");
               TLine t2 = TLine(chargeStart[iCh]*digiTime, 0, chargeStart[iCh]*digiTime, peakTmp);  t2.SetLineColor(kBlue); t2.Draw("same");
               TLine t3 = TLine(chargeStop[iCh]*digiTime, 0, chargeStop[iCh]*digiTime, peakTmp);  t2.SetLineColor(kBlue); t2.Draw("same");
               TLine t4 = TLine(0, peakTmp, peakTimeTmp, peakTmp);  t4.SetLineColor(kPink); t4.Draw("same");
               TLine t5 = TLine(tmin, 0, tmin, peakTmp);  t5.SetLineColor(kGreen); t5.Draw("same");
               TLine t6 = TLine(tmax, 0, tmax, peakTmp);  t6.SetLineColor(kGreen); t6.Draw("same");
               cc.Write(); 
               outFile->cd();
            }

         } //channels loop

      cryTree->Fill();
      } //entries loop
   //Loop

   //Terminate
      HM->ProcessBoxes(); 
      cryTree->Write();
      outFile->Close();
   //Terminate

}


