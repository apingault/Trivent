#include "HistogramHendler.h"

void CollTree::CreateHistos()
{
  //  Create histograms... the bins, names are  held in histoTH1F
 
  Int_t i=0;
  m_histoTH1F = new TClonesArray("TH1F",5);
  TClonesArray &a=*m_histoTH1F;  
  while (histoTH1F[i].nbins) {
    a[i]=new(a[i]) TH1F( histoTH1F[i].name,
        histoTH1F[i].title,
        histoTH1F[i].nbins,
        histoTH1F[i].xlow,
        histoTH1F[i].xhigh);
    i++;
  }
 
  m_histoTH2F = new TClonesArray("TH2F",8);
  TClonesArray &b=*m_histoTH2F;
  i=0;
  while (histoTH2F[i].n1bins) {
    b[i]=new(b[i]) TH2F( histoTH2F[i].name,
        histoTH2F[i].title,
        histoTH2F[i].n1bins,
        histoTH2F[i].x1low,
        histoTH2F[i].x1high,
        histoTH2F[i].n2bins,
        histoTH2F[i].x2low,
        histoTH2F[i].x2high);
    i++;
  };    
  m_histoTProf = new TClonesArray("TProfile",7);
  TClonesArray &c=*m_histoTProf;
  i=0;
  while (histoTProf[i].nbins) {
    c[i]=new(c[i]) TProfile( histoTProf[i].name,
        histoTProf[i].title,
        histoTProf[i].nbins,
        histoTProf[i].x1low,
        histoTProf[i].x1high,
        histoTProf[i].x2low,
        histoTProf[i].x2high);
    //			     histoTProf->at(i).x2high,
    //			     "s");
    i++;
  };
 
}
//_________________________________________________________________________
void CollTree::StoreHistos(const char *name)
{ // writes all the histograms individually defined in m_histoTH1F in the same root file
  TFile f(name,"RECREATE");
 
  //store histograms
 
  m_histoTH1F->Write();
  m_histoTH2F->Write();
  m_histoTProf->Write();
}
 
//_________________________________________________________________________
TH1F *CollTree::GetHistoTH1F(std::string name)
{// returns the histogram specified by the number
  return (TH1F*) m_histoTH1F->FindObject(name.c_str());
}
 
//_________________________________________________________________________
TH2F *CollTree::GetHistoTH2F(std::string name)
{
  return (TH2F*) m_histoTH2F->FindObject(name.c_str());
} 
//_________________________________________________________________________
TProfile *CollTree::GetHistoTProf(std::string name)
{
  return (TProfile*) m_histoTProf->FindObject(name.c_str());
}
