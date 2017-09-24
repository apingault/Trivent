#define EventDisplay_cxx
// The class definition in EventDisplay.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> EventTree->Process("EventDisplay.C")
// root> T->Process("EventDisplay.C","some options")
// root> T->Process("EventDisplay.C+")
//


#include "EventDisplay.h"
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TStopwatch.h"
#include <iostream>

#include <sstream>

TStopwatch w; 
void EventDisplay::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
  
   TString option = GetOption();
   
   for(int iThresh=1; iThresh<4; ++iThresh)
     {
       std::map <std::string, TH2*> mTH2;

       std::ostringstream ssName;
       ssName << "ij_"<<iThresh;
       TH2* ijTH2 = new TH2F(ssName.str().c_str(),"IvJ Display ",96,0,95,96,0,95);
       ssName.str(std::string());
       ssName << "ik_"<<iThresh;
       TH2* ikTH2 = new TH2F(ssName.str().c_str(),"IvK Display ",50,0,50,96,0,95);
       ssName.str(std::string());
       ssName << "jk_"<<iThresh;
       TH2* jkTH2 = new TH2F(ssName.str().c_str(),"JvK Display ",50,0,50,96,0,95);

       mTH2.insert({"ij", ijTH2});
       mTH2.insert({"ik", ikTH2});
       mTH2.insert({"jk", jkTH2});
       m_mThreshTH2.insert({iThresh, mTH2});
     }
     
   Double_t w = 600;
   Double_t h = 800;
   c1 = new TCanvas("c", "c", w, h);
   // c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));
   c1->Divide(1,3);
}

void EventDisplay::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
  //  ijH2 = new TH2F("ijH2","IvJ Display ",96,0,95,96,0,95);
  //  ikH2 = new TH2F("ikH2","IvK Display ",50,0,50,96,0,95);
  //  jkH2 = new TH2F("jkH2","JvK Display ",50,0,50,96,0,95);
  //  fOutput->Add(h1);
  evt=0;


}

Bool_t EventDisplay::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);
   ++evt;
  //  std::cout << "processing evt '" << evt << "'" << std::endl;
   if (!*EventIsSelected || *CerenkovTime<-10 || *CerenkovTime>-4)
      return kTRUE;
   
  //  if (*EventNumber == 100)
  //   Abort("nik");
    
  // if (evt>100)
  //   Abort("renik");


   // Threshold 1
    // ijH2 = new TH2F("ijH2","IvJ Display ",96,0,95,96,0,95);
    // ikH2 = new TH2F("ikH2","IvK Display ",50,0,50,96,0,95);
    // jkH2 = new TH2F("jkH2","JvK Display ",50,0,50,96,0,95);

    // // Threshold 2
    // ijH2_1 = new TH2F("ijH2_1","IvJ Display ",96,0,95,96,0,95);
    // ikH2_1 = new TH2F("ikH2_1","IvK Display ",50,0,50,96,0,95);
    // jkH2_1 = new TH2F("jkH2_1","JvK Display ",50,0,50,96,0,95);
    
    // // Threshold 3
    // ijH2_2 = new TH2F("ijH2_2","IvJ Display ",96,0,95,96,0,95);
    // ikH2_2 = new TH2F("ikH2_2","IvK Display ",50,0,50,96,0,95);
    // jkH2_2 = new TH2F("jkH2_2","JvK Display ",50,0,50,96,0,95);

    for ( unsigned int iHit = 0; iHit < *NumberOfHits; ++iHit)
      {
        // std::cout << "iHit " << iHit << " HitThreshold[iHit] " << HitThreshold[iHit] << std::endl;
	      //  std::map<std::string, TH2*> toto = 
          // HitThreshold goes from 
         m_mThreshTH2.at(HitThreshold[iHit]).at("ij")->Fill(HitI[iHit],HitJ[iHit]);
  
	m_mThreshTH2.at(HitThreshold[iHit]).at("ik")->Fill(HitK[iHit],HitI[iHit]);
	m_mThreshTH2.at(HitThreshold[iHit]).at("jk")->Fill(HitK[iHit],HitJ[iHit]);
	// std::cout << HitI[iHit] << " " << HitJ[iHit] << " " << HitK[iHit] << std::endl;
      }
   
    std::cout
      << " Evt: " << *EventNumber
      << "\tnHit: " << *NumberOfHits
      << "\t " << *NumberOfCerenkov1Hits
      << " " << *NumberOfCerenkov2Hits
      << " " << *NumberOfCerenkov3Hits
      << " Th : " << *CerThreshold
      << " Asic : " << *CerAsic
      << " Chan : " << *CerChan
      << " TooMany : " << *TriggerHasTooManyCerenkov
      << " CerTime : " << *CerenkovTime;
      // << std::endl;


    for(int iThresh=1; iThresh<4; ++iThresh)
      {
	Color_t color;
	switch (iThresh){
	case 1:
	  color = kGreen; break;
	case 2:
	  color = kBlue; break;
	case 3:
	  color = kRed; break;
	default:
    std::cout << "Wrong Threshold" << std::endl;
	  color = kBlack; break;
	}
	  
	c1->cd(1);
	m_mThreshTH2.at(iThresh).at("ij")->SetMarkerColor(color);
	m_mThreshTH2.at(iThresh).at("ij")->SetMarkerStyle(kFullDotLarge);
	m_mThreshTH2.at(iThresh).at("ij")->SetMarkerSize(.25);
	m_mThreshTH2.at(iThresh).at("ij")->Draw("same");

	c1->cd(2);
	m_mThreshTH2.at(iThresh).at("ik")->SetMarkerColor(color);
  m_mThreshTH2.at(iThresh).at("ik")->SetMarkerStyle(kFullDotLarge);
  m_mThreshTH2.at(iThresh).at("ik")->SetMarkerSize(.5);
	m_mThreshTH2.at(iThresh).at("ik")->Draw("same");
  
	c1->cd(3);
	m_mThreshTH2.at(iThresh).at("jk")->SetMarkerColor(color);
  m_mThreshTH2.at(iThresh).at("jk")->SetMarkerStyle(kFullDotLarge);
  m_mThreshTH2.at(iThresh).at("jk")->SetMarkerSize(.25);
	m_mThreshTH2.at(iThresh).at("jk")->Draw("same");

	//c1->WaitPrimitive();
  }
  c1->Update();
	int ch = std::cin.get();
	if (ch == 27) // ESC
	  Abort("Legraclavie");


 
    // getchar();
    //  // }

        for(int iThresh=1; iThresh<4; ++iThresh)
	  {
	    m_mThreshTH2.at(iThresh).at("ij")->Reset();
	    m_mThreshTH2.at(iThresh).at("ik")->Reset();
	    m_mThreshTH2.at(iThresh).at("jk")->Reset();
	    // m_mThreshTH2[iThresh].clear();
	  }
	// m_mThreshTH2.clear();
	return kTRUE;
}

void EventDisplay::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void EventDisplay::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
      std::cout << "Terminate.. draw histograms" << std::endl;
      gPad->WaitPrimitive();

      w.Stop();
      std::cout << "Time for query "; 
      w.Print();

     delete c1; 
      
  //  c1->Divide(1,3);
  //  c1->cd(1);
  //  TObject * oh1 = fOutput->FindObject("h1");
  //  if (oh1) oh1->Draw();
  //  c1->cd(2);
  //  TObject * oh2 = fOutput->FindObject("h2");
  //  if (oh2) oh2->Draw();
  //  c1->cd(3);
  //  TObject * oh3 = fOutput->FindObject("h3");
  //  if (oh3) oh3->Draw();
}
