#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <EVENT/RawCalorimeterHit.h>

#include <iostream>
#include <map>
#include <cmath>
#include "Mapping.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

struct GainHist{
  TH1F *noise_on_dif;
  TH1F *gain_chan;
};

class GainCorrect{
 public:
  GainCorrect();
  ~GainCorrect();
  // methods
  uint    getCellDif_id(int cell_id);
  uint    getCellAsic_id(int cell_id);
  uint    getCellChan_id(int cell_id);
  
  void    gainCorrectionInit(std::map<int, LayerID  > _mapping);
  void    gainCorrectionParser(EVENT::LCCollection* col);
  void    gainCorrector(const string gain_file);
    
 protected:
  std::map<int, LayerID  > _mapping;
  std::map<int, GainHist > _h_gain ;
  
  TH1F *noise_dist;
  TH1F *gain_chan;
  TH1F *mean_hit_dif;
  TH1F *time_hit_dif;
};  

class MyGainCorrection{
 public:
  MyGainCorrection();
  ~MyGainCorrection();

  uint getCellDif_id(int cell_id){
    return cell_id & 0xFF;
  }
  uint getCellAsic_id(int cell_id){
    return (cell_id & 0xFF00)>>8;
  }
  uint getCellChan_id(int cell_id){
    return (cell_id & 0x3F0000)>>16;
  }

  void analyseNoise();
  void writeChannelGainCorrectionFile();
  void writeAsicGainCorrectionFile();
  void writeDifGainCorrectionFile();
  void fillHitsInMap(EVENT::RawCalorimeterHit* &raw_hit);
  void createControlHistogram();
 private : 
  float mpv;
  float meanHitNoise;
  float meanAsicNoise;
  float meanDifNoise;
  float rmsHitNoise;
  float rmsAsicNoise;
  float rmsDifNoise;
  // map : first = key (cellID, asicKey, difKey) ; second = counter for each key
  // cellID = lcio cellid0
  // asicKey = dif_id*100+asic_id
  // difKey = dif_id
  std::map<int,int> noiseHitMap;
  std::map<int,int> noiseAsicMap;
  std::map<int,int> noiseDifMap;
  TFile *ctrlFile;
  TTree *tree;
  int noiseLevel;
  int med;
  float gainFeDi;
  float gainLin;
  float gainFeDiMed;
};
