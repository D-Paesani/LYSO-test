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

const int posNo = 6;
const int cryNo = 4;
const int chNo = 3;
const int sideNo = 2;
const int singleCrystal = 2; //set as negative integer to disable
const int singleSide = 1;

TString inFileFormat = "lysoData/data2/gap_cry%d_sd%d_pos%d.root";
TString outFileFormat = "lysoData/data3/gap_cry%d_sd%d.root";
// TString inFileFormat = "lysoData/data2/cry%d_sd%d_pos%d.root";
// TString outFileFormat = "lysoData/data3/cry%d_sd%d.root";
TString outAllFileFormat = "lysoData/data3/cryAll.root";

const TString treeName = "crytree";
const TString tiDiffStr = "0.5*wTime[1] + 0.5*wTime[2] - wTime[0]";
const TString tiMeanStr = "0.5*wTime[1] + 0.5*wTime[2]";
const TString meanChaStr = "0.5*wCharge[1] + 0.5*wCharge[2]";
const TString asymStr = "(wCharge[1] - wCharge[2])/(wCharge[1] + wCharge[2])";

struct _resVals { double peak, peakErr, sigma, sigmaErr, resolution, resError, npe, npeErr; };

struct _resCry {
    _resVals left[posNo], right[posNo], mean[posNo];
    TGraphErrors resoPosL, resoPosR, resoPosMean;
    TGraphErrors chargePosL, chargePosR, chargePosMean;
    TGraphErrors npePos, npeQ;
};

_resCry resCry[cryNo][sideNo];
TGraphErrors resoPosAll[sideNo], chargePosAll[sideNo], npePosAll[sideNo], npeQ[sideNo], tdiffPos[sideNo];

int _cry = 0, _side = 0; 
_resCry *thisRes;

void branch2histo1d( TH1* hObj, TTree* tree, TString val, TString cut = "" ){
    TH1 *hist;
    TString Val = Form( "%s>>histo(%d, %f, %f)", val.Data(), (int)hObj->GetNbinsX(), (double)hObj->GetXaxis()->GetXmin(), (double)hObj->GetXaxis()->GetXmax());
    tree->Draw( Val, cut.Data(), "goff");
    hist = (TH1*)gDirectory->Get("histo"); 
    hObj->Add(hist);
    hist->Delete();
}

int fillGraph(TGraphErrors *g, double x, double y, double ex, double ey) {
    g->AddPoint(x, y);
    int n = g->GetN();
    g->SetPointError(n-1, ex, ey);
    return n;    
}

int fitChargeGausStandard(TH1 *hist, double& peak, double &peakErr, double &sigma, double &sigmaErr, double sigmaToFit, TString options) {
    
    gStyle->SetOptFit(1111);
    hist->SetStats(1111); 

    peak = hist->GetBinCenter(hist->GetMaximumBin());
    sigma = hist->GetRMS(); 
    
    TF1 gaus = TF1("g", "gaus");
    gaus.SetParameter(1, peak); 
    gaus.SetParameter(2, sigma);
    if ( (hist->Fit(&gaus, "EMQ")) < 0 ) {return -1;}     
    peak = gaus.GetParameter(1); 
    sigma = gaus.GetParameter(2);
    gaus.SetParameter(1, peak); 
    gaus.SetParameter(2, sigma);
    gaus.SetRange(peak - sigmaToFit*sigma, peak + sigmaToFit*sigma); 
    if ( (hist->Fit(&gaus, "REMQ")) < 0 ) {return -1;}     
    peak = gaus.GetParameter(1); 
    sigma = gaus.GetParameter(2);
    peakErr = gaus.GetParError(1);
    sigmaErr = gaus.GetParError(2); 
    gaus.SetRange(peak - sigmaToFit*sigma, peak + sigmaToFit*sigma); 
    
    return 1;
}

void namerNone(int hN, TString& hTag, TString& hTitleTag) {
    hTag = "";
    hTitleTag = "";
}

void namerCha(int hN, TString& hTag, TString& hTitleTag) {
    hTag = Form("_ch%d", hN); 
    hTitleTag = Form(" cry%d sd%d ch%d", _cry, _side, hN); 
}

void namerPos(int hN, TString& hTag, TString& hTitleTag) {
    hTag = Form("_pos%d", hN); 
    hTitleTag = Form(" cry%d sd%d pos%d", _cry, _side, hN); 
}

void procCharges(TH1* histObj, int histN, int& histSkipFlag) {   

    TString name = histObj->GetName();
    gStyle->SetOptFit(1111); 
    gStyle->SetOptStat(1111); 

    double min_tmp = histObj->GetXaxis()->GetXmin();
    double max_tmp = histObj->GetXaxis()->GetXmax();
    histObj->GetXaxis()->SetRangeUser(500, 2000);

    double peak, epeak, sigma, esigma;
    fitChargeGausStandard(histObj, peak, epeak, sigma, esigma, 2.0, "");
    histObj->GetXaxis()->SetRangeUser(min_tmp, max_tmp);

    double resolution = sigma/peak;
    double resoError = resolution*TMath::Sqrt((epeak/peak)*(epeak/peak) + (esigma/sigma)*(esigma/sigma));
    int posGraph = histN + 1;

    if (name.Contains("chargeMeanCut")) {
        thisRes->mean[histN].peak = peak;
        thisRes->mean[histN].peakErr = epeak;
        thisRes->mean[histN].sigma = sigma;
        thisRes->mean[histN].sigmaErr = esigma;
        fillGraph( &(thisRes->resoPosMean), posGraph, resolution, 0, resoError);
        fillGraph( &(thisRes->chargePosMean), posGraph, peak, 0, epeak);
        fillGraph( &resoPosAll[_side], (double)_cry + 1 - 0.2 + ((double)posGraph-1)*0.4/((double)posNo-1), resolution, 0, resoError); 
        fillGraph( &chargePosAll[_side], (double)_cry + 1 - 0.2 + ((double)posGraph-1)*0.4/((double)posNo-1), peak, 0, epeak); 
    } else if (name.Contains("chargeCutL")) {
        thisRes->left[histN].peak = peak;
        thisRes->left[histN].peakErr = epeak;
        thisRes->left[histN].sigma = sigma;
        thisRes->left[histN].sigmaErr = esigma;
        fillGraph( &(thisRes->resoPosL), posGraph, resolution, 0, resoError);
        fillGraph( &(thisRes->chargePosL), posGraph, peak, 0, epeak);
    } else if (name.Contains("chargeCutR")){
        thisRes->right[histN].peak = peak;
        thisRes->right[histN].peakErr = epeak;
        thisRes->right[histN].sigma = sigma;
        thisRes->right[histN].sigmaErr = esigma;
        fillGraph( &(thisRes->resoPosR), posGraph, resolution, 0, resoError);
        fillGraph( &(thisRes->chargePosR), posGraph, peak, 0, epeak);
    }

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

    TString name = histObj->GetName();
    gStyle->SetOptFit(1111); 
    gStyle->SetOptStat(1111); 

    double peak, epeak, sigma, esigma;
    fitChargeGausStandard(histObj, peak, epeak, sigma, esigma, 2.2, "");
    double npe = 1/(2*sigma*sigma);
    double errNpe = npe*1.4*(esigma/sigma)*(esigma/sigma);
    int posGraph = histN + 1;

    thisRes->mean[histN].npe = npe;
    thisRes->mean[histN].npeErr = errNpe;

    fillGraph( &(thisRes->npePos), posGraph, npe, 0, errNpe);
    fillGraph( &npePosAll[_side], (double)_cry + 1 - 0.2 + ((double)posGraph-1)*0.4/((double)posNo-1), npe, 0, errNpe); 
    //fillGraph( &(thisRes->npeQ), posGraph, thisRes->mean[histN].peak, 0, errNpe);
}

void Process() {

    gErrorIgnoreLevel = kFatal;
    
    for (_cry = 0; _cry < cryNo; _cry++) {

        if (!(singleCrystal<0)) { _cry = singleCrystal; }

        for (_side = 0; _side < sideNo; _side++) { 

            if (!(singleSide<0)) { _side = singleSide; }
            
            thisRes = &(resCry[_cry][_side]);

            cout<<"-----> Creating  "<<Form(outFileFormat, _cry, _side)<<endl;
            TFile* outFile = new TFile(Form(outFileFormat, _cry, _side), "RECREATE"); 
            
            HistManager *HM = new HistManager();
            HM->SetOutFile(outFile);

            const int tiBins = 2000;
            const int qBins = 1500;
            const int qTo = 6000;

            HM->AddHistBox("th1f", chNo, "blineAll", "blineAll", "", "mV", 500, -50, 50, HM->GetProcDef(), &namerCha); //indexing on channels
            HM->AddHistBox("th1f", chNo, "blineRmsAll", "blineRmsAll", "", "mV", 300, 0, 6, HM->GetProcDef(), &namerCha);
            HM->AddHistBox("th1f", chNo, "pedAll", "pedAll", "Q", "pC", 600, -50, 50, HM->GetProcDef(), &namerCha);

            HM->AddHistBox("th1f", posNo, "chargeRawL", "chargeRawL", "Q", "pC", qBins, 0, qTo, &procPeaks, &namerPos); //indexing on position
            HM->AddHistBox("th1f", posNo, "chargeRawR", "chargeRawR", "Q", "pC", qBins, 0, qTo, &procPeaks, &namerPos); 
            HM->AddHistBox("th1f", posNo, "chargeCutL", "chargeCutL", "Q", "pC", qBins, 0, qTo, &procCharges, &namerPos); 
            HM->AddHistBox("th1f", posNo, "chargeCutR", "chargeCutR", "Q", "pC", qBins, 0, qTo, &procCharges, &namerPos); 
            HM->AddHistBox("th1f", posNo, "chargeMeanRaw", "chargeMeanRaw", "Q", "pC", qBins, 0, qTo, &procPeaks, &namerPos); 
            HM->AddHistBox("th1f", posNo, "chargeMeanCut", "chargeMeanCut", "Q", "pC", qBins, 0, qTo, &procCharges, &namerPos); 

            HM->AddHistBox("th1f", posNo, "timeRawL", "timeRawL", "Q", "pC", 2000, 0, 500, HM->GetProcDef(), &namerPos); 
            HM->AddHistBox("th1f", posNo, "timeRawR", "timeRawR", "Q", "pC", 2000, 0, 500, HM->GetProcDef(), &namerPos); 
            HM->AddHistBox("th1f", posNo, "timeMeanRaw", "timeMeanRaw", "Q", "pC", 2000, 0, 500, HM->GetProcDef(), &namerPos); 
            HM->AddHistBox("th1f", posNo, "timeMeanCut", "timeMeanCut", "Q", "pC", 2000, 0, 500, HM->GetProcDef(), &namerPos);
            HM->AddHistBox("th1f", posNo, "timeDiffRaw", "timesDiff", "time", "ns", tiBins, -100, 100, HM->GetProcDef(), &namerPos); 
            HM->AddHistBox("th1f", posNo, "timeDiffCut", "timesDiff", "time", "ns", tiBins, -100, 100, HM->GetProcDef(), &namerPos); 

            HM->AddHistBox("th1f", posNo, "asym", "asym", "", "", 600, -3, 3, &procAsym, &namerPos); 

            HM->AddHistBox("th1f", 1, "chargeTagAll", "chargeTagAll", "Q", "pC", qBins, 0, qTo, &procPeaks, &namerNone); 
            HM->AddHistBox("th1f", posNo, "chargeTag", "chargeTag", "Q", "pC", qBins, 0, qTo, HM->GetProcDef(), &namerPos); 
            HM->AddHistBox("th1f", posNo, "timeTagRaw", "timeTag", "time", "ns", tiBins, 0, 500,  HM->GetProcDef(), &namerPos);         
            HM->AddHistBox("th1f", posNo, "timeTagCut", "timeTag", "time", "ns", tiBins, 0, 500,  HM->GetProcDef(), &namerPos);      

            for (int _pos = 0; _pos < posNo; _pos++) {   

                TString inName = Form(inFileFormat, _cry, _side, _pos);
                cout<<"--------> Opening  "<<inName<<endl;
                TFile *inFile = new TFile(inName.Data());
                inFile->cd();
                TTree *inTree = (TTree*)inFile->Get(treeName.Data());

                branch2histo1d( HM->GetHist("chargeRawL", _pos), inTree, "wCharge[1]", "");
                branch2histo1d( HM->GetHist("chargeRawR", _pos), inTree, "wCharge[2]", "");
                branch2histo1d( HM->GetHist("chargeMeanRaw", _pos), inTree, meanChaStr, "");

                branch2histo1d( HM->GetHist("timeRawL", _pos), inTree, "wTime[1]", "");
                branch2histo1d( HM->GetHist("timeRawR", _pos), inTree, "wTime[2]", "");
                branch2histo1d( HM->GetHist("timeMeanRaw", _pos), inTree, tiMeanStr, "");
                branch2histo1d( HM->GetHist("timeDiffRaw", _pos), inTree, tiDiffStr, "");
                
                branch2histo1d( HM->GetHist("chargeTagAll", 0), inTree, "wCharge[0]", "");
                branch2histo1d( HM->GetHist("chargeTag", _pos), inTree, "wCharge[0]", "");
                branch2histo1d( HM->GetHist("timeTagRaw", _pos), inTree, "wTime[0]", "");

                TH1F *tagChTmp = (TH1F*)HM->GetHist("chargeTag", _pos);
                TH1F *tiDiTmp = (TH1F*)HM->GetHist("timeDiffCut", _pos);

                double sigmaCharge = 1.8;
                double sigmaTime = 2.4;

                double minq_ = tagChTmp->GetXaxis()->GetXmin(), maxq_ = tagChTmp->GetXaxis()->GetXmax();
                tagChTmp->GetXaxis()->SetRangeUser(800, 1400);
                double qPeak, qPeakErr, qSigma, qSigmaErr;
                fitChargeGausStandard(tagChTmp, qPeak, qPeakErr, qSigma, qSigmaErr, sigmaCharge, "");
                tagChTmp->GetXaxis()->SetRangeUser(minq_, maxq_);

                double cutqmin = qPeak - sigmaCharge*qSigma, cutqmax = qPeak + sigmaCharge*qSigma;
                TString selCut = Form("wCharge[0] > %f && wCharge[0] < %f", cutqmin, cutqmax);
                double qminboth = 550;
                selCut += Form(" && wCharge[1] > %f && wCharge[2] > %f", qminboth, qminboth);

                branch2histo1d( HM->GetHist("timeMeanCut", _pos), inTree, tiMeanStr, selCut);
                branch2histo1d( HM->GetHist("timeTagCut", _pos), inTree, "wTime[0]", selCut);
                branch2histo1d( tiDiTmp, inTree, tiDiffStr, selCut);
        
                double mint_ = tiDiTmp->GetXaxis()->GetXmin(), maxt_ = tiDiTmp->GetXaxis()->GetXmax();
                tiDiTmp->GetXaxis()->SetRangeUser(10, 35);
                double tpeak, tepeak, tsigma, tesigma;
                fitChargeGausStandard(tiDiTmp, tpeak, tepeak, tsigma, tesigma, sigmaTime, "");
                tiDiTmp->GetXaxis()->SetRangeUser(mint_, maxt_);
                
                double cuttmin = tpeak - sigmaTime*tsigma, cuttmax = tpeak + sigmaTime*tsigma;
                selCut = Form( "%s && %s > %f && %s < %f", selCut.Data(), tiDiffStr.Data(), cuttmin, tiDiffStr.Data(), cuttmax);
                cout<<"----> cut: "<<selCut<<endl;

                branch2histo1d( HM->GetHist("chargeCutL", _pos), inTree, "wCharge[1]", selCut);
                branch2histo1d( HM->GetHist("chargeCutR", _pos), inTree, "wCharge[2]", selCut);
                branch2histo1d( HM->GetHist("chargeMeanCut", _pos), inTree, meanChaStr, selCut);
                branch2histo1d( HM->GetHist("asym", _pos), inTree, asymStr, selCut);  

                for (int _cha = 0; _cha < chNo; _cha++) { 
                    branch2histo1d(HM->GetHist("blineRmsAll", _cha), inTree, Form("wBrms[%d]", _cha), selCut);
                    branch2histo1d(HM->GetHist("blineAll", _cha), inTree, Form("wBline[%d]", _cha), selCut);
                    branch2histo1d(HM->GetHist("pedAll", _cha), inTree, Form("wPed[%d]", _cha), selCut);
                }
            }
        
            HM->ProcessBoxes();

            TDirectory *resoDir = outFile->mkdir("results");
            resoDir->cd();

            thisRes->chargePosL.Write("chargePosL");
            thisRes->chargePosR.Write("chargePosR");
            thisRes->chargePosMean.Write("chargePosMean");
            thisRes->resoPosL.Write("resoPosL");
            thisRes->resoPosR.Write("resoPosR");
            thisRes->resoPosMean.Write("resoPosMean");
            thisRes->npePos.Write("npePos");
            thisRes->npeQ.Write("npeQ");
 
            HM->CloseOutFile();

            if (!(singleSide<0)) { break; }

        } //for side
        
        if (!(singleCrystal<0)) { break; }

    } //for cry

    if (!(singleCrystal<0)) { return; }

    TFile* outAllFile = new TFile(outAllFileFormat, "RECREATE"); 
    outAllFile->cd();

    TCanvas cc; cc.cd();
    resoPosAll[0].SetMarkerStyle(25); resoPosAll[0].SetMarkerSize(1.9); resoPosAll[0].SetMarkerColor(kRed);
    resoPosAll[1].SetMarkerStyle(24); resoPosAll[1].SetMarkerSize(1.9); resoPosAll[1].SetMarkerColor(kBlue);
    resoPosAll[0].SetTitle("Resolution (sq->sd0, cir->sd1)"); 
    resoPosAll[0].GetYaxis()->SetTitle("sigmaQ/Q"); 
    resoPosAll[0].GetXaxis()->SetTitle("Crystal ID");
    resoPosAll[0].GetXaxis()->SetRangeUser(0, cryNo + 1); 
    resoPosAll[0].GetXaxis()->SetNdivisions(cryNo+1);
    resoPosAll[0].GetYaxis()->SetRangeUser(0, 0.2);
    if (!(singleSide<0)) { 
        resoPosAll[singleSide].Draw("AP");
    } else {
        resoPosAll[0].Draw("AP");
        resoPosAll[1].Draw("P same");
    }
    cc.Write("ResolutionOverview");

    TCanvas ccc; ccc.cd();
    chargePosAll[0].SetMarkerStyle(25); chargePosAll[0].SetMarkerSize(1.9); chargePosAll[0].SetMarkerColor(kRed);
    chargePosAll[1].SetMarkerStyle(24); chargePosAll[1].SetMarkerSize(1.9); chargePosAll[1].SetMarkerColor(kBlue);
    chargePosAll[0].SetTitle("Charge (sq->sd0, cir->sd1)"); 
    chargePosAll[0].GetYaxis()->SetTitle("Q"); 
    chargePosAll[0].GetXaxis()->SetTitle("Crystal ID");
    chargePosAll[0].GetXaxis()->SetRangeUser(0, cryNo + 1); 
    chargePosAll[0].GetXaxis()->SetNdivisions(cryNo+1);
    chargePosAll[0].GetYaxis()->SetRangeUser(0, 0.2);
    if (!(singleSide<0)) { 
        chargePosAll[singleSide].Draw("AP");
    } else {
        chargePosAll[0].Draw("AP");
        chargePosAll[1].Draw("P same");
    }
    ccc.Write("ChargeOverview");

    TCanvas cccc; cccc.cd();
    npePosAll[0].SetMarkerStyle(25); npePosAll[0].SetMarkerSize(1.9); npePosAll[0].SetMarkerColor(kRed);
    npePosAll[1].SetMarkerStyle(24); npePosAll[1].SetMarkerSize(1.9); npePosAll[1].SetMarkerColor(kBlue);
    npePosAll[0].SetTitle("npe (sq->sd0, cir->sd1)"); 
    npePosAll[0].GetYaxis()->SetTitle("Q"); 
    npePosAll[0].GetXaxis()->SetTitle("Crystal ID");
    npePosAll[0].GetXaxis()->SetRangeUser(0, cryNo + 1); 
    npePosAll[0].GetXaxis()->SetNdivisions(cryNo+1);
    npePosAll[0].GetYaxis()->SetRangeUser(0, 0.2);
    if (!(singleSide<0)) { 
        npePosAll[singleSide].Draw("AP");
    } else {
        npePosAll[0].Draw("AP");
        npePosAll[1].Draw("P same");
    }
    cccc.Write("NpeOverview");

    outAllFile->Close();
}

