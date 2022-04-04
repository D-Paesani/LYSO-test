//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 24 18:27:20 2022 by ROOT version 6.24/06
// from TTree Wavefull/Wavefull
// found on file: /Users/dp/Documents/Programmi/ROOT/LYSO/data1/lyso_cr0_sd0/lyso_cr0_sd0_pos0full.root
//////////////////////////////////////////////////////////

#ifndef wfAna_h
#define wfAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class wfAna : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> ntrig = {fReader, "ntrig"};
   TTreeReaderValue<Int_t> evnum = {fReader, "evnum"};
   TTreeReaderValue<Int_t> num = {fReader, "num"};
   TTreeReaderArray<Double_t> time = {fReader, "time"};
   TTreeReaderArray<Double_t> onda0 = {fReader, "onda0"};
   TTreeReaderValue<Double_t> bkg_charge0 = {fReader, "bkg_charge0"};
   TTreeReaderValue<Double_t> charge0 = {fReader, "charge0"};
   TTreeReaderValue<Double_t> tmean0 = {fReader, "tmean0"};
   TTreeReaderArray<Double_t> onda1 = {fReader, "onda1"};
   TTreeReaderValue<Double_t> bkg_charge1 = {fReader, "bkg_charge1"};
   TTreeReaderValue<Double_t> charge1 = {fReader, "charge1"};
   TTreeReaderValue<Double_t> tmean1 = {fReader, "tmean1"};
   TTreeReaderValue<Double_t> charge_peak = {fReader, "charge_peak"};
   TTreeReaderValue<Double_t> charge_peak2 = {fReader, "charge_peak2"};
   TTreeReaderArray<Double_t> onda2 = {fReader, "onda2"};
   TTreeReaderValue<Double_t> bkg_charge2 = {fReader, "bkg_charge2"};
   TTreeReaderValue<Double_t> charge2 = {fReader, "charge2"};
   TTreeReaderValue<Double_t> tmean2 = {fReader, "tmean2"};


   wfAna(TTree * /*tree*/ =0) { }
   virtual ~wfAna() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(wfAna,0);

};

#endif

#ifdef wfAna_cxx
void wfAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t wfAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef wfAna_cxx
