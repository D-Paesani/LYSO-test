//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May  9 11:03:59 2022 by ROOT version 6.24/06
// from TTree Wavefull/Wavefull
// found on file: lysoData/data1/Lyso2_s0/Lyso2__pos1full.root
//////////////////////////////////////////////////////////

#ifndef waveAna_h
#define waveAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class waveAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ntrig;
   Int_t           evnum;
   Int_t           num;
   Double_t        time[1024];   //[num]
   Double_t        onda0[1024];   //[num]
   Double_t        bkg_charge0;
   Double_t        charge0;
   Double_t        tmean0;
   Double_t        onda1[1024];   //[num]
   Double_t        bkg_charge1;
   Double_t        charge1;
   Double_t        tmean1;
   Double_t        charge_peak;
   Double_t        charge_peak2;
   Double_t        onda2[1024];   //[num]
   Double_t        bkg_charge2;
   Double_t        charge2;
   Double_t        tmean2;

   // List of branches
   TBranch        *b_ntrig;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_num;   //!
   TBranch        *b_time;   //!
   TBranch        *b_onda0;   //!
   TBranch        *b_bkg_charge0;   //!
   TBranch        *b_charge0;   //!
   TBranch        *b_tmean0;   //!
   TBranch        *b_onda1;   //!
   TBranch        *b_bkg_charge1;   //!
   TBranch        *b_charge1;   //!
   TBranch        *b_tmean1;   //!
   TBranch        *b_charge_peak;   //!
   TBranch        *b_charge_peak2;   //!
   TBranch        *b_onda2;   //!
   TBranch        *b_bkg_charge2;   //!
   TBranch        *b_charge2;   //!
   TBranch        *b_tmean2;   //!

   waveAna(TTree *tree=0);
   virtual ~waveAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef waveAna_cxx
waveAna::waveAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("lysoData/data1/Lyso2_s0/Lyso2__pos1full.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("lysoData/data1/Lyso2_s0/Lyso2__pos1full.root");
      }
      f->GetObject("Wavefull",tree);

   }
   Init(tree);
}

waveAna::~waveAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t waveAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t waveAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void waveAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ntrig", &ntrig, &b_ntrig);
   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("num", &num, &b_num);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("onda0", onda0, &b_onda0);
   fChain->SetBranchAddress("bkg_charge0", &bkg_charge0, &b_bkg_charge0);
   fChain->SetBranchAddress("charge0", &charge0, &b_charge0);
   fChain->SetBranchAddress("tmean0", &tmean0, &b_tmean0);
   fChain->SetBranchAddress("onda1", onda1, &b_onda1);
   fChain->SetBranchAddress("bkg_charge1", &bkg_charge1, &b_bkg_charge1);
   fChain->SetBranchAddress("charge1", &charge1, &b_charge1);
   fChain->SetBranchAddress("tmean1", &tmean1, &b_tmean1);
   fChain->SetBranchAddress("charge_peak", &charge_peak, &b_charge_peak);
   fChain->SetBranchAddress("charge_peak2", &charge_peak2, &b_charge_peak2);
   fChain->SetBranchAddress("onda2", onda2, &b_onda2);
   fChain->SetBranchAddress("bkg_charge2", &bkg_charge2, &b_bkg_charge2);
   fChain->SetBranchAddress("charge2", &charge2, &b_charge2);
   fChain->SetBranchAddress("tmean2", &tmean2, &b_tmean2);
   Notify();
}

Bool_t waveAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void waveAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t waveAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef waveAna_cxx
