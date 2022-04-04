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

#include "HistManager.h"

using namespace std;

TString inFileFormat = "data/data2/cry%d_sd%d_pos%d.root";
TString outFileFormat = "data/data3/cry_%d.root";
TString outAllFileFormat = "data/data3/cryAll.root";
TString treeName = "crytree";

const TString tiDiffStr = "0.5*wTime[1] + 0.5*wTime[2] - wTime[0]";
const TString meanChaStr = "0.5*wCharge[1] + 0.5*wCharge[2]";
const TString asymStr = "(wCharge[1] - wCharge[2])/(wCharge[1] + wCharge[2])";
TString superCut = "";

const int posNo = 3;
const int cryNo = 4;
const int chNo = 3;
const int sideNo = 2;
int _cry = 0;

double fmean[cryNo][sideNo][posNo], fmeanerr[cryNo][sideNo][posNo], fsigma[cryNo][sideNo][posNo], fsigmaerr[cryNo][sideNo][posNo];

void branch2histo1d( TH1* hObj, TTree* tree, TString val, TString cut = "" ){
    TH1 *hist;
    TString Val = Form( "%s>>histo(%d, %f, %f)", val.Data(), (int)hObj->GetNbinsX(), (double)hObj->GetXaxis()->GetXmin(), (double)hObj->GetXaxis()->GetXmax());
    tree->Draw( Val, cut.Data(), "goff");
    hist = (TH1*)gDirectory->Get("histo"); 
    hObj->Add(hist);
    hist->Delete();
}
        
void namerSidePos(int hN, TString& hTag, TString& hTitleTag) {
    int side = (int)((hN+1)>posNo); int pos = (int)(hN-(side==1)*posNo); 
    hTag = Form("_sd%d_pos%d", side, pos); 
    hTitleTag = Form(" cry%d sd%d pos%d", _cry, side, pos); 
}

void namerNone(int hN, TString& hTag, TString& hTitleTag) {
    hTag = "";
    hTitleTag = "";
}

void namerCha(int hN, TString& hTag, TString& hTitleTag) {
    hTag = Form("_ch%d", hN); 
    hTitleTag = Form(" cry%d ch%d", _cry, hN); 
}

void procCharges(TH1* histObj, int histN, int& histSkipFlag) {   

    gStyle->SetOptFit(1111); 
    gStyle->SetOptStat(1111); 

    // TSpectrum tsp(1);
    // tsp.Search(histObj,  5, "noMarkov");
    // double *pks = tsp.GetPositionX();
    // double qpeak = *std::max_element(pks, pks + 2);
    double min_tmp = histObj->GetXaxis()->GetXmin();
    double max_tmp = histObj->GetXaxis()->GetXmax();
    histObj->GetXaxis()->SetRangeUser(500,1000);
    double qpeak = histObj->GetBinCenter(histObj->GetMaximumBin());
    histObj->GetXaxis()->SetRangeUser(min_tmp, max_tmp);

    double qmax = qpeak + 80, qmin = qpeak - 50; float pk,sigma;
    histObj->SetTitle(Form("%s -> [%.0f]", histObj->GetTitle(), qpeak));

    TF1 fit1 = TF1("g", "gaus", qmin, qmax); fit1.SetParameter(1, qpeak); fit1.SetParameter(2, 10); histObj->Fit(&fit1, "R"); 
    pk = fit1.GetParameter(1); sigma = fit1.GetParameter(2);
    fit1 = TF1("g", "gaus", pk - 2.5*sigma, pk + 2.5*sigma); fit1.SetParameter(1, pk); fit1.SetParameter(2, sigma); histObj->Fit(&fit1, "R"); 
    pk = fit1.GetParameter(1); sigma = fit1.GetParameter(2);

    int side = (int)((histN+1)>posNo); int pos = (int)(histN-(side==1)*posNo);
    fmean[_cry][side][pos] = fit1.GetParameter(1);
    fmeanerr[_cry][side][pos] = fit1.GetParError(1);
    fsigma[_cry][side][pos] = fit1.GetParameter(2);
    fsigmaerr[_cry][side][pos] = fit1.GetParError(2);
}

void procPeaks(TH1* histObj, int histN, int& histSkipFlag) { 
   if (!(TString(histObj->GetName()).Contains("Tag"))) {return;}  
    gStyle->SetOptFit(1); 
    TSpectrum tsp(4);
    tsp.Search(histObj,  5, "noMarkov");
    double *pks = tsp.GetPositionX();
    histObj->SetTitle(Form("%s -> pks[ %.0f, %.0f, %.0f]", histObj->GetTitle(), pks[1], pks[2], pks[3]));
}

void procAsym(TH1* histObj, int histN, int& histSkipFlag) {   
}

void Process() {
    
    for (_cry = 0; _cry < cryNo; _cry++) {

        HistManager *HM = new HistManager();

        TFile* outFile = new TFile(Form(outFileFormat, _cry), "RECREATE"); 
        HM->SetOutFile(outFile);

        const int tiBins = 2000;
        const int qBins = 1500;
        const int qTo = 6000;

        HM->AddHistBox("th1f", posNo*2, "chargeMeanRaw", "chargeMeanRaw", "Q", "pC", qBins, 0, qTo, &procPeaks, &namerSidePos); 
        HM->AddHistBox("th1f", posNo*2, "chargeCut", "chargeCut", "Q", "pC", 250, 100, 3100, &procCharges, &namerSidePos); 
        HM->AddHistBox("th1f", posNo*2, "asym", "asym", "", "", 400, -5, 5, &procAsym, &namerSidePos); 
        HM->AddHistBox("th1f", 1, "timeTagAll", "timeTag", "time", "ns", tiBins, 0, 500); 
        HM->AddHistBox("th1f", 1, "timePeakTagAll", "timePeakTag", "time", "ns", tiBins, 0, 500); 
        HM->AddHistBox("th1f", chNo, "blineAll", "blineAll", "", "mV", 500, -50, 50, HM->GetProcDef(), &namerCha); 
        HM->AddHistBox("th1f", chNo, "blineRmsAll", "blineRmsAll", "", "mV", 300, 0, 6, HM->GetProcDef(), &namerCha);
        HM->AddHistBox("th1f", chNo, "pedAll", "pedAll", "Q", "pC", 600, -100, 100, HM->GetProcDef(), &namerCha);
        HM->AddHistBox("th1f", posNo*2, "timesDiff", "timesDiff", "time", "ns", tiBins, -100, 100, HM->GetProcDef(), &namerSidePos); 
        HM->AddHistBox("th1f", posNo*2, "chargeTag", "chargeTag", "Q", "pC", qBins, 0, qTo, HM->GetProcDef(), &namerSidePos); 
        HM->AddHistBox("th1f", 1, "chargeTagAll", "chargeTagAll", "Q", "pC", qBins, 0, qTo, &procPeaks); 

        for (int _side = 0 ; _side < 2; _side++) {
            for (int _pos = 0; _pos < posNo; _pos++) {   

                int SidePos = _pos + _side*posNo;

                TString inName = Form(inFileFormat, _cry, _side, _pos); 
                TFile *inFile = new TFile(inName.Data());
                inFile->cd();
                TTree *inTree = (TTree*)inFile->Get(treeName.Data());

                for (int _cha = 0; _cha < chNo; _cha++) { 
                    int SidePosCha = 0;
                    // HM->AddHistBox("th1f", posNo, Form("times_sd%d_ch%d", _side, _cha), Form("times sd%d ch%d", _side, _cha), "time", "ns", tiBins, 0, 500, HM->GetProcDef());
                    // HM->AddHistBox("th1f", posNo, Form("timesPeak_sd%d_ch%d", _side, _cha), Form("timesPeak sd%d ch%d", _side, _cha), "time", "ns", tiBins, 0, 500, HM->GetProcDef());
                    // HM->AddHistBox("th1f", posNo, Form("chargesRaw_sd%d_ch%d", _side, _cha), Form("chargesRaw sd%d ch%d", _side, _cha), "Q", "pC", qBins, 0, qTo, HM->GetProcDef()); 
                    // branch2histo1d(HM->GetHist(Form("times_sd%d_ch%d", _side, _cha), _pos), inTree, Form("wTime[%d]",_cha));
                    // branch2histo1d(HM->GetHist(Form("timesPeak_sd%d_ch%d", _side, _cha), _pos), inTree, Form("wTimePeak[%d]", _cha));
                    // branch2histo1d(HM->GetHist(Form("chargesRaw_sd%d_ch%d", _side, _cha), _pos), inTree, Form("wCharge[%d]", _cha));
                    branch2histo1d(HM->GetHist("blineRmsAll", _cha), inTree, Form("wBrms[%d]", _cha));
                    branch2histo1d(HM->GetHist("blineAll", _cha), inTree, Form("wBline[%d]", _cha));
                    branch2histo1d(HM->GetHist("pedAll", _cha), inTree, Form("wPed[%d]", _cha));
                }

                branch2histo1d(HM->GetHist("chargeTagAll", 0), inTree, "wCharge[0]");
                branch2histo1d(HM->GetHist("chargeTag", SidePos), inTree, "wCharge[0]");
                branch2histo1d(HM->GetHist("timesDiff", SidePos), inTree, tiDiffStr.Data());

                TH1F *tagChTmp = (TH1F*)HM->GetHist("chargeTag", SidePos);
                TH1F *tiDiTmp = (TH1F*)HM->GetHist("timesDiff", SidePos);

                gStyle->SetOptFit(1); 
                double minq_ = tagChTmp->GetXaxis()->GetXmin(), maxq_ = tagChTmp->GetXaxis()->GetXmax();
                tagChTmp->GetXaxis()->SetRangeUser(800, 1400);
                double pkq = tagChTmp->GetBinCenter(tagChTmp->GetMaximumBin());
                tagChTmp->GetXaxis()->SetRangeUser(minq_, maxq_);
                double maxq = pkq + 80, minq = pkq - 70;
                TF1 fitq = TF1("g", "gaus", minq, maxq); fitq.SetParameter(1, pkq); fitq.SetParameter(2, 40); tagChTmp->Fit(&fitq, "R"); 
                fitq = TF1("g", "gaus", fitq.GetParameter(1) - 2.5*fitq.GetParameter(2), fitq.GetParameter(1) + 2.5*fitq.GetParameter(2)); fitq.SetParameter(1, fitq.GetParameter(1)); fitq.SetParameter(2, fitq.GetParameter(2)); tagChTmp->Fit(&fitq, "R"); 
                double cutqmin = fitq.GetParameter(1) - 2.5*fitq.GetParameter(2), cutqmax = fitq.GetParameter(1) + 2.5*fitq.GetParameter(2);

                double mint_ = tiDiTmp->GetXaxis()->GetXmin(), maxt_ = tiDiTmp->GetXaxis()->GetXmax();
                tiDiTmp->GetXaxis()->SetRangeUser(10, 35);
                double pkt = tiDiTmp->GetBinCenter(tiDiTmp->GetMaximumBin());
                tiDiTmp->GetXaxis()->SetRangeUser(mint_, maxt_);
                double maxt = pkt + 10, mint = pkt - 10;
                TF1 fitt = TF1("g", "gaus", mint, maxt); fitt.SetParameter(1, pkt); fitt.SetParameter(2, 2); tiDiTmp->Fit(&fitt, "R"); 
                fitt = TF1("g", "gaus", fitt.GetParameter(1) - 3*fitt.GetParameter(2), fitt.GetParameter(1) + 3*fitt.GetParameter(2)); fitt.SetParameter(1, fitt.GetParameter(1)); fitt.SetParameter(2, fitt.GetParameter(2)); tiDiTmp->Fit(&fitt, "R"); 
                double cuttmin = fitt.GetParameter(1) - 1.5*fitt.GetParameter(2), cuttmax = fitt.GetParameter(1) + 1.5*fitt.GetParameter(2);

                superCut = Form("wCharge[0] > %f && wCharge[0] < %f", cutqmin, cutqmax);
                superCut += Form( " && %s > %f && %s < %f", tiDiffStr.Data(), cuttmin, tiDiffStr.Data(), cuttmax);
                // superCut += "wBline[1] < 20 && wBline[2] < 20 && wBline[2] > - 20 && wBline[1] > -20 && ";
                // superCut += "wBrms[1] < 3 && wBrms[2] < 3";
                cout<<"--------------------------{cut}---> "<<superCut<<endl;

                branch2histo1d(HM->GetHist("asym", SidePos), inTree, asymStr.Data(), superCut.Data());
                branch2histo1d(HM->GetHist("chargeMeanRaw", SidePos), inTree, meanChaStr.Data());
                branch2histo1d(HM->GetHist("timePeakTagAll", 0), inTree, "wTimePeak[0]");
                branch2histo1d(HM->GetHist("timeTagAll", 0), inTree, "wTime[0]");
                branch2histo1d(HM->GetHist("chargeCut", SidePos), inTree, meanChaStr.Data(), superCut.Data());

            }  
        }
    
        HM->ProcessBoxes();
        HM->CloseOutFile();
        
    }

    TFile* outAllFile = new TFile(outAllFileFormat, "RECREATE"); 
    outAllFile->cd();

    TGraphErrors resoAll[2]; // = {TGraphErrors(posNo*cryNo)};
    TGraphErrors chargeAll[2]; // = {TGraphErrors(posNo*cryNo)};

    for (int iCry = 0; iCry < cryNo; iCry++) {
        for (int iPos = 0; iPos < posNo; iPos++) {    
            for (int iSd = 0; iSd < sideNo; iSd++) {               
                double sig = fsigma[iCry][iSd][iPos], esig = fsigmaerr[iCry][iSd][iPos];
                double mea = fmean[iCry][iSd][iPos], emea = fmeanerr[iCry][iSd][iPos];
                double err = TMath::Sqrt((emea/mea)*(emea/mea) + (esig/sig)*(esig/sig));
                resoAll[iSd].SetPoint(iPos + posNo*iCry, (double)iCry + 1 - 0.2 + iPos*0.4/(posNo-1), sig/mea );
                resoAll[iSd].SetPointError(iPos + posNo*iCry, 0, (sig/mea)*err);
                chargeAll[iSd].SetPoint(iPos + posNo*iCry, (double)iCry + 1 - 0.2 + iPos*0.4/(posNo-1), mea);
                chargeAll[iSd].SetPointError(iPos + posNo*iCry, 0, emea);
            }
        }
    }

    TCanvas cc; cc.cd();
    resoAll[0].SetMarkerStyle(25); resoAll[0].SetMarkerColor(kRed); resoAll[0].SetMarkerSize(1.9);
    resoAll[1].SetMarkerStyle(24); resoAll[1].SetMarkerSize(1.9); resoAll[0].SetMarkerColor(kBlue);
    resoAll[0].SetTitle("Resolution (sq->sd0, cir->sd1)"); resoAll[0].GetYaxis()->SetTitle("sigmaQ/Q"); resoAll[0].GetXaxis()->SetTitle("Crystal ID");
    resoAll[0].GetXaxis()->SetRangeUser(0, cryNo + 1); resoAll[0].GetXaxis()->SetNdivisions(cryNo+1);
    resoAll[0].GetYaxis()->SetRangeUser(0, 0.2);
    resoAll[0].Draw("P");
    resoAll[1].Draw("P same");
    cc.Write("ResolutionOverview");

    TCanvas ccc; ccc.cd();
    chargeAll[0].GetYaxis()->SetTitle("Q [pC]"); chargeAll[0].GetXaxis()->SetTitle("Crystal ID");
    chargeAll[0].SetMarkerStyle(25); 
    chargeAll[0].SetMarkerColor(kRed);
    chargeAll[1].SetMarkerStyle(24); 
    chargeAll[1].SetMarkerSize(1.9);
    chargeAll[0].SetMarkerColor(kBlue);
    chargeAll[0].SetMarkerSize(1.9);
    chargeAll[0].SetTitle("Charges (sq->sd0, cir->sd1)");
    chargeAll[0].Draw("P");
    chargeAll[1].Draw("P same");
    ccc.Write("ChargeAll");

    outAllFile->Close();


}

