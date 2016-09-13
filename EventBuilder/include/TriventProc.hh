#ifndef _TriventProc_hh_
#define _TriventProc_hh_

// -- Asics and channels mapping for sdhcal
#include "Mapping.h"

// -- std includes
#include <string>
#include <map>

// -- marlin includes
#include <marlin/Processor.h>
#include "marlin/VerbosityLevels.h"
#include "marlin/tinyxml.h"

// -- lcio includes
#include <lcio.h>
#include "IO/LCWriter.h"
#include <EVENT/RawCalorimeterHit.h>
#include <EVENT/LCCollection.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>

// -- ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>

class TriventProc  : public marlin::Processor
{
public:

  Processor*  newProcessor() { return new TriventProc ; }

  TriventProc();

  ~TriventProc() {};

  void init();


  void    processEvent( LCEvent * evtP );
  void    processRunHeader( LCRunHeader * runH);// added by me
  void    XMLReader(std::string xmlfile);
  void    readDifGeomFile(std::string geomfile);
  void    printDifGeom();
  void    defineColors();

  uint    getCellDif_id(int cell_id);
  uint    getCellAsic_id(int cell_id);
  uint    getCellChan_id(int cell_id);

  void    getMaxTime();
  std::vector<int> getTimeSpectrum();
  std::vector<unsigned int> getPadIndex(uint dif_id, uint asic_id, uint chan_id);
  void    eventBuilder(LCCollection* col_event, int time_peak, int prev_time_peak);
  void    end();

protected:
  TH1F *noise_dist;
  TH1F *gain_chan;
  TH1F *mean_hit_dif;
  TH1F *time_hit_dif;
  // xml test
  std::map<std::string, std::string> m_parameters;

  std::vector<EVENT::RawCalorimeterHit*> _trigger_raw_hit;

  bool GAIN_CORRECTION_MODE;
  std::string _outFileName;
  std::string _noiseFileName;
  std::string _treeName;
  std::string _rootFileName;
  std::string _colName;
  std::string _fileName;
  std::string _mappingfile;
  std::string _geomXML;
  std::vector<std::string> _hcalCollections;
  int _overwrite;
  TTree *_outputTree;
  unsigned int _eventNr;
  Int_t _nHit;
  Int_t _elecNoiseCut;

  std::map<int, LayerID  > _mapping;
  std::map<int, double  > _chamberPos;//camber , pos

  // Cut parameters
  int _noiseCut;
  int _timeWin;
  int _layerCut;
  int _time2prevEventCut;
  double _layerGap;
  float _beamEnergy;
  int _trigCount;
  int _maxTime;
  int evtnum;
  int _selectedNum;
  int _rejectedNum;
  std::set<unsigned int> _firedLayersSet;
  LCWriter* _lcWriter;
  int _bcid1;
  int _bcid2;

  // Color for streamlog ouptut 
  std::string normal  ;
  std::string red     ;
  std::string green   ;
  std::string yellow  ;
  std::string blue    ;
  std::string magenta ;
  std::string white   ;


//ROOT histograms
TFile *m_rootFile; 
std::vector<TH2D*> m_vHitMapPerLayer;
Int_t m_runNumber;

};



#endif


