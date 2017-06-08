#ifndef _TCaloHitRoot_hh_
#define _TCaloHitRoot_hh_


#include <TObject.h>
#if (!defined(__CINT__)) 
#include <EVENT/CalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#endif

class TCaloHitRoot : public TObject
{
public:
  TCaloHitRoot() {};


#if (!defined(__CINT__)) 
  TCaloHitRoot(const EVENT::CalorimeterHit &_hit, CellIDDecoder<CalorimeterHit> *idDecoder) {
    //    std::cout<<" on debute le maping"<<std::endl;
    
    Int_t ID0,ID1;
    
    ID0=_hit.getCellID0();
    ID1=_hit.getCellID1();
    
    // ========== Decode CellID Informations =======
    I = (*idDecoder)(&_hit)["I"];
    J = (*idDecoder)(&_hit)["J"];
    K = (*idDecoder)(&_hit)["K"];
    fDif_id  = (*idDecoder)(&_hit)["Dif_id"];
    fAsic_id = (*idDecoder)(&_hit)["Asic_id"];
    fChan_id = (*idDecoder)(&_hit)["Chan_id"];
    // ========== get yhe position =======
        
    fX = _hit.getPosition()[0]/10.; // X coordinate ( cm )
    fY = _hit.getPosition()[1]/10.; // Y coordinate ( cm ) 
    fZ = _hit.getPosition()[2]/10.; // Z coordinate ( cm ) 
    
    // ======= Get Time and Thr =============
    fTime  = _hit.getTime();
    //event  = _hit.getEvent
    fThr   = int(_hit.getEnergy()); // the threshold
    
  };
#endif

  
  Double_t GetTime() const          { return fTime;};
  void     SetTime(Double_t time)   { fTime = time;};
  
  Double_t GetX() {return fX;}
  Double_t GetY() {return fY;}
  Double_t GetZ() {return fZ;}
  
protected:
 
  Double_t fX;
  Double_t fY; 
  Double_t fZ; 
  
  Int_t I;
  Int_t J;
  Int_t K;

  // for simulation
  //int event;
  unsigned int fDif_id;
  int fAsic_id;
  int fChan_id;
  Double_t fTime;
  Int_t fThr;
  
  ClassDef(TCaloHitRoot,1);
};
#endif 


