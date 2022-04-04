#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <vector>
#include <iostream>
#include <TApplication.h>
#include <TH1.h>
#include <stdlib.h>
#include <TFrame.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <vector>
#include <TString.h>
#include <TTimer.h>
#include <TLatex.h>
#include <TGClient.h>
#include <TRootCanvas.h>
#include <TDatime.h>
#include <time.h>

using namespace std;

//MODIFIED FOR TITLE TAG APPEND!!!!!!!!!!!!!!!!!!!!

TString _label_form = "%s [%s]";
TString _entries_form = "Entries / %.2f %s";
TString _entriesNoUnit_form = "Entries";

/* TEMPLATES
  void processHistBox(TH1* hObj, int hN, int& hSkip) {}
  void namerHistBox(int hN, TString& hTag, TString& hTitleTag) {}
*/

static void processFunc_def(TH1* hObj, int hN, int& hSkip) {};

static void namerFunc_def(int n, TString& hTag, TString& hTitleTag) { 
  hTag = Form("_%d", n); 
  hTitleTag = Form(" [%d]", n); 
};

typedef void (*ProcessFunction)(TH1*, int, int&);
typedef void (*NamerFunction)(int, TString&, TString&);


class HistBox{

  private:

    TDirectory  *_histosDir;
    TFile *_outfile;
    TH1 **_histosObj;
    TString _histType, _label1, _label2, _histosTitle, _histosName, hTag, hTitleTag;
    int _histN, hSkip;
    ProcessFunction _processFunction; 
    NamerFunction _namerFunction;

  public:
   
    HistBox(
      TFile *f,
      TString histType, int histN,
      TString histosName, TString histosTitle, TString label1, TString label2, TString label3, TString label4,
      int xBins, double xMin, double xMax, 
      int yBins, double yMin, double yMax,
      ProcessFunction processFunc,
      NamerFunction namerFunc
    ):  
      _namerFunction(namerFunc),
      _processFunction(processFunc),
      _histosTitle(histosTitle),
      _histType(histType),
      _histN(histN),
      _histosName(histosName),
      _outfile(f)
    { 
      cout<<"-----> HistManager: HistBox constructor: ["<<_histType<<"] ["<<_histN<<"] ["<<_histosName<<"] ["<<_histosTitle<<"]"<<endl;

      _histosObj = new TH1*[_histN];

      if (_histType == "th1f") { 

        if (label2 == "") {
          _label2 = _entriesNoUnit_form;
          _label1 = label1;
        } else {
          _label2 = Form(_entries_form, ((xMax-xMin)/(double)xBins), label2.Data());
          _label1 = Form(_label_form, label1.Data(), label2.Data());
        }

        for(int k=0; k<_histN; k++){
          _namerFunction(k, hTag, hTitleTag);
          _histosObj[k] = new TH1F(_histosName.Data() + hTag, _histosTitle.Data() + hTitleTag, xBins, xMin, xMax);
        }

      } else if (histType == "th2f") { 

        if (label2 == "") { _label1 = label1; } 
        else { _label1 = Form(_label_form, label1.Data(), label2.Data()); }
        
        if (label4 == "") { _label2 = label3; } 
        else { _label2 = Form(_label_form, label3.Data(), label4.Data()); }

        for(int k=0; k<_histN; k++){
          _namerFunction(k, hTag, hTitleTag);  
          _histosObj[k] = new TH2F(_histosName.Data() + hTag, _histosTitle.Data() + hTitleTag, xBins, xMin, xMax, yBins, yMin, yMax); 
        }      
      
      } else {
        cout<< "-----> HistManager: HistBox error: histType not supported"<<endl; 
        exit(EXIT_FAILURE); 
      }
    }

    ~HistBox() {}

    TH1** GetHistosObj(){ return _histosObj; }
    TH1* GetHistObj(int k){ return _histosObj[k]; }

    TDirectory* GetDirectory(){ return _histosDir; }


    void ProcessBox(){
      
      _outfile->cd(); 
      _histosDir = _outfile->mkdir(_histosName.Data(), "recreate"); 
      _histosDir->cd();

      for (int k = 0; k < _histN; k++) {

        _histosObj[k]->GetXaxis()->SetTitle(_label1);
        _histosObj[k]->GetYaxis()->SetTitle(_label2);
        
        _histosDir->cd();

        hSkip = 0; 
        _processFunction(_histosObj[k], k, hSkip);
        if (hSkip) {continue;}  

        _histosDir->cd(); 
        _histosObj[k]->Write();
      } 
    }

};


class HistManager{

  private:

    TFile *_outfile;

    ProcessFunction procUser_def = &processFunc_def;
    NamerFunction namerUser_def = &namerFunc_def;
    ProcessFunction proc_def = &processFunc_def;
    NamerFunction namer_def = &namerFunc_def;

  public:

    HistManager() {}
    ~HistManager() {}
    
    void SetNamerFun(NamerFunction namerFunc) { namerUser_def = namerFunc; }
    void SetProcessFun(ProcessFunction processFunc) { procUser_def = processFunc; }
    NamerFunction GetNamerUser() { return namerUser_def; }
    ProcessFunction GetProcUser() { return procUser_def;}
    NamerFunction GetNamerDef() { return namer_def; }
    ProcessFunction GetProcDef() { return proc_def;}

    void SetOutFile(TFile *outfile) { _outfile = outfile; }
    TFile *GetOutFIle() { return _outfile; }
    void CloseOutFile() { _outfile->Close(); }

    unordered_map<string, HistBox *> HistBoxes;

    void AddHistBox(TString, int, TString, TString, TString, TString, TString, TString, int, double, double, int, double, double, ProcessFunction, NamerFunction);
    void AddHistBox(TString, int, TString, TString, TString, TString, TString, TString, int, double, double, int, double, double, ProcessFunction);
    void AddHistBox(TString, int, TString, TString, TString, TString, TString, TString, int, double, double, int, double, double);
    void AddHistBox(TString, int, TString, TString, TString, TString, int, double, double, ProcessFunction, NamerFunction);
    void AddHistBox(TString, int, TString, TString, TString, TString, int, double, double, ProcessFunction);
    void AddHistBox(TString, int, TString, TString, TString, TString, int, double, double);

    TH1 **GetHistos(string name) { return (TH1**)HistBoxes[name]->GetHistosObj(); }
    TH1 *GetHist(string name, int n) { return (TH1*)HistBoxes[name]->GetHistObj(n); }

    void Fill2d(string name, int n, double val1, double val2) { ((TH2F*)GetHist(name, n))->Fill(val1, val2); }
    void Fill1d(string name, int n, double val1) { ((TH1F*)GetHist(name, n))->Fill(val1); }
    
    void ProcessBox(string name) { HistBoxes[name]->ProcessBox(); }

    void ProcessBoxes() {
      for (auto& HistBox: HistBoxes) {
        HistBox.second->ProcessBox(); 
      }
    };

};


void HistManager::AddHistBox( TString histType, int histN, TString histosName, TString histTitle, TString label1, TString label2, TString label3, TString label4,
  int xBins, double xMin, double xMax, int yBins, double yMin, double yMax,
  ProcessFunction processFunc, NamerFunction namerFunc
) {
  HistManager::HistBoxes.insert( {histosName.Data(), 
    new HistBox(_outfile, histType, histN, histosName, histTitle, label1, label2, label3, label4, xBins, xMin, xMax, yBins, yMin, yMax, processFunc, namerFunc)
  } );
};

void HistManager::AddHistBox( TString histType, int histN, TString histosName, TString histTitle, TString label1, TString label2, TString label3, TString label4,
  int xBins, double xMin, double xMax, int yBins, double yMin, double yMax,
  ProcessFunction processFunc
) { HistManager::AddHistBox(histType, histN, histosName, histTitle, label1, label2, label3, label4, xBins, xMin, xMax, yBins, yMin, yMax, processFunc, namerUser_def); };

void HistManager::AddHistBox( TString histType, int histN, TString histosName, TString histTitle, TString label1, TString label2, TString label3, TString label4,
  int xBins, double xMin, double xMax, int yBins, double yMin, double yMax
) { HistManager::AddHistBox(histType, histN, histosName, histTitle, label1, label2, label3, label4, xBins, xMin, xMax, yBins, yMin, yMax, procUser_def, namerUser_def); };

void HistManager::AddHistBox( TString histType, int histN, TString histosName, TString histTitle, TString label1, TString label2,
  int xBins, double xMin, double xMax,
  ProcessFunction processFunc, NamerFunction namerFunc
) { HistManager::AddHistBox(histType, histN, histosName, histTitle, label1, label2, "", "", xBins, xMin, xMax, 0, 0, 0, processFunc, namerFunc); };

void HistManager::AddHistBox( TString histType, int histN, TString histosName, TString histTitle, TString label1, TString label2,
  int xBins, double xMin, double xMax,
  ProcessFunction processFunc
) { HistManager::AddHistBox(histType, histN, histosName, histTitle, label1, label2, "", "", xBins, xMin, xMax, 0, 0, 0, processFunc, namerUser_def); };

void HistManager::AddHistBox( TString histType, int histN, TString histosName, TString histTitle, TString label1, TString label2,
  int xBins, double xMin, double xMax
) { HistManager::AddHistBox(histType, histN, histosName, histTitle, label1, label2, "", "", xBins, xMin, xMax, 0, 0, 0, procUser_def, namerUser_def); };










