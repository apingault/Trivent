#ifndef _HistogramHendler_h_
#define _HistogramHendler_h_

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TProfile.h>
#include <TClonesArray.h>
 
class HistDefTH1F{
  public:
    Text_t *name;
    Text_t *title;
    Int_t nbins;
    Float_t xlow;
    Float_t xhigh;
};
 
class HistDefTH2F{
  public: 
    Text_t *name;
    Text_t *title;
    Int_t n1bins;
    Float_t x1low;
    Float_t x1high;
    Int_t n2bins;
    Float_t x2low;
    Float_t x2high;
};
 
class HistDefTProf{
  public: 
    Text_t *name;
    Text_t *title;
    Int_t nbins;
    Float_t x1low;
    Float_t x1high;
    Float_t x2low;
    Float_t x2high;
};
//-------
HistDefTH1F histoTH1F[] = {
  { "h_noise_dist"  ,  "Noise distribution"       , 200, 0. , 200. },
  { "h_time_2_event",  "Time to previous event"   , 10001 ,0 ,10001},
  { 0,0,0,0.0,0.0}
};

HistDefTH2F histoTH2F[] = {
  {"hist2002","Number of Tracks Vs Energy",5,-0.5,4.5,100,0.0,200.0},
  { 0,0,0,0.0,0.0,0,0.0 ,0.0}
};

HistDefTProf histoTProf[] = {
  {"hist10101","   ",20, 0.0, 100.0, 0.0, 100.0},    
  {"hist10102","   ",20, -3.0, 3.0, 0.0, 100.0},
  {"hist10103","   ",20, -3.5, 3.5, 0.0, 100.0},
  {"hist10111","   ",20, 0.0, 100.0, 0.0, 100.0},
  {0,0,0,0,0,0,0}
};

class CollTree{
  void CreateHistos();
  void StoreHistos(const char *name);
  TH1F       *GetHistoTH1F(std::string name);  // returns Histogram with name "histX"
  TH2F       *GetHistoTH2F(std::string name);  // returns Histogram with name "histX"
  TProfile   *GetHistoTProf(std::string name);  // returns Histogram with name "histX"
 private:
  TClonesArray      *m_histoTH1F;  // list of histograms TH1F
  TClonesArray      *m_histoTH2F;  // list of histograms TH2F
  TClonesArray      *m_histoTProf; // list of TProfiles
  
};

#endif
