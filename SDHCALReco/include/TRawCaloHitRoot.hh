#ifndef _TRawCaloHitRoot_hh_
#define _TRawCaloHitRoot_hh_


#include <TObject.h>
#include "Mapping.h" // mapping of ASICs and DIFs
#if (!defined(__CINT__)) 
#include <EVENT/CalorimeterHit.h>
#include <EVENT/RawCalorimeterHit.h>
#include <UTIL/CellIDDecoder.h>
#endif

class TRawCaloHitRoot : public TObject
{
public:
  TRawCaloHitRoot() {};
  
#if (!defined(__CINT__)) 
  TRawCaloHitRoot(const EVENT::RawCalorimeterHit &_hit, std::map<int,LayerID> mapping) {
    //    std::cout<<" on debute le maping"<<std::endl;
    Int_t ID0,ID1;
    
    ID0=_hit.getCellID0();
    ID1=_hit.getCellID1();
    
    // ========== Decode CellID Informations =======
    fDif_id=ID0 & 0xFF; // Dif id
    fAsic_id=(ID0 & 0xFF00)>>8; // Asic id
    fChan_id=(ID0 & 0x3F0000)>>16; // Channel id
    // ========== get the position =======
    fX = (1+MapILargeHR2[fChan_id]+AsicShiftI[fAsic_id])*1.04125;
    fY = (32-(MapJLargeHR2[fChan_id]+AsicShiftJ[fAsic_id]))*1.04125;
    double DifX,DifY,DifZ;
    for(std::map<int,LayerID>::iterator itt =mapping.begin();itt!=mapping.end();itt++)
      {
	if(itt->first == fDif_id){ 
	  DifZ = itt->second.K ;
	  DifY = itt->second.DifY ;
	}
      }
     I = (1+MapILargeHR2[fChan_id]+AsicShiftI[fAsic_id]);
     J = (32-(MapJLargeHR2[fChan_id]+AsicShiftJ[fAsic_id])) + int(DifY);
     K = int(DifZ) ;
     
     fX = fX +DifY*1.0421;
     fZ = DifZ * _layer_gap;
     // ======= Get Time and Thr =============
     fThr = _hit.getAmplitude()&3;
     fTime =_hit.getTimeStamp();
    
  };
#endif

  
  Double_t GetTime() const          { return fTime;};
  void     SetTime(Double_t time)   { fTime = time;};
  
  Double_t GetX() {return fX;}
  Double_t GetY() {return fY;}
  Double_t GetZ() {return fZ;}
  void SetLayerGap(double gap=2.8) { _layer_gap = gap;}
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
  Double_t  _layer_gap; 
  ClassDef(TRawCaloHitRoot,1)
};
#endif 
