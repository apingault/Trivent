//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 10 14:39:37 2017 by ROOT version 6.07/07
// from TTree EventTree/Event variables
// found on file: /eos/user/a/apingaul/CALICE/Data/SPS_12_2014/Trivent/TDHCAL_Antoine_726255POS.root
//////////////////////////////////////////////////////////

#ifndef EventDisplay_h
#define EventDisplay_h

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// Headers needed by this particular selector
#include <map>
#include <vector>

class EventDisplay : public TSelector {
public:
  TTreeReader fReader; //! the tree reader
  TTree *fChain = 0;   //! pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<UInt_t> TriggerNumber = {fReader, "TriggerNumber"};
  TTreeReaderValue<UInt_t> EventNumber = {fReader, "EventNumber"};
  TTreeReaderValue<UInt_t> NumberOfHits = {fReader, "NumberOfHits"};
  TTreeReaderArray<int> HitI = {fReader, "HitI"};
  TTreeReaderArray<int> HitJ = {fReader, "HitJ"};
  TTreeReaderArray<int> HitK = {fReader, "HitK"};
  TTreeReaderArray<int> HitThreshold = {fReader, "HitThreshold"};
  TTreeReaderValue<UInt_t> NumberOfFiredLayers = {fReader, "NumberOfFiredLayers"};
  TTreeReaderValue<UInt_t> NumberOfCerenkov1Hits = {fReader, "NumberOfCerenkov1Hits"};
  TTreeReaderValue<UInt_t> NumberOfCerenkov2Hits = {fReader, "NumberOfCerenkov2Hits"};
  TTreeReaderValue<UInt_t> NumberOfCerenkov3Hits = {fReader, "NumberOfCerenkov3Hits"};
  TTreeReaderValue<UInt_t> CerAsic = {fReader, "CerAsic"};
  TTreeReaderValue<UInt_t> CerChan = {fReader, "CerChan"};
  TTreeReaderValue<UInt_t> CerThreshold = {fReader, "CerThreshold"};
  TTreeReaderValue<Bool_t> TriggerHasTooManyCerenkov = {fReader, "TriggerHasTooManyCerenkov"};
  TTreeReaderValue<Int_t> CerenkovTime = {fReader, "CerenkovTime"};
  TTreeReaderValue<Bool_t> EventIsSelected = {fReader, "EventIsSelected"};
  TTreeReaderValue<Bool_t> EventIsNoise = {fReader, "EventIsNoise"};
  TTreeReaderValue<Bool_t> EventIsToCloseFromLast = {fReader, "EventIsToCloseFromLast"};
  TTreeReaderValue<Bool_t> EventHasNotEnoughLayers = {fReader, "EventHasNotEnoughLayers"};
  TTreeReaderValue<Bool_t> EventIsHasFullAsic = {fReader, "EventIsHasFullAsic"};

  EventDisplay(TTree * /*tree*/ = 0) : c1(nullptr), evt(-1) {}
  virtual ~EventDisplay() {}
  virtual Int_t Version() const { return 2; }
  virtual void Begin(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual void Init(TTree *tree);
  virtual Bool_t Notify();
  virtual Bool_t Process(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }
  virtual void SetOption(const char *option) { fOption = option; }
  virtual void SetObject(TObject *obj) { fObject = obj; }
  virtual void SetInputList(TList *input) { fInput = input; }
  virtual TList *GetOutputList() const { return fOutput; }
  virtual void SlaveTerminate();
  virtual void Terminate();

private:
  // std::map<threshold, std::map<histoName, histo>>
  std::map<int, std::map<std::string, TH2*>> m_mThreshTH2;

  /* TH2* ijH2; */
  /* TH2* ikH2; */
  /* TH2* jkH2; */
  /* // Threshold 2 */
  /* TH2* ijH2_1; */
  /* TH2* ikH2_1; */
  /* TH2* jkH2_1; */

  /* // Threshold 3 */
  /* TH2* ijH2_2; */
  /* TH2* ikH2_2; */
  /* TH2* jkH2_2; */

  TCanvas *c1;
  int evt;
  ClassDef(EventDisplay, 0);
};

#endif

#ifdef EventDisplay_cxx
void EventDisplay::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t EventDisplay::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef EventDisplay_cxx
