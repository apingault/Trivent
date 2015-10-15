#ifndef _TriventProc_hh_
#define _TriventProc_hh_

#define  HISTOGRAM_PARSER true

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/RawCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include "IO/LCWriter.h"
#include <map>
#include <algorithm>
#include "Mapping.h"
using namespace std;

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

  uint    getCellDif_id(int cell_id);
  uint    getCellAsic_id(int cell_id);
  uint    getCellChan_id(int cell_id);

  void    getMaxTime();
  std::vector<int> getTimeSpectrum();
  uint*   getPadIndex(uint dif_id, uint asic_id, uint chan_id);
  void    eventBuilder(LCCollection* col_event,int time_peak, int prev_time_peak);
  bool    peakOrNot(std::vector<int> time_spectrum, int itime ,int threshold);
  void    end();

protected:
  TH1F *noise_dist;
  TH1F *gain_chan;
  TH1F *mean_hit_dif;
  TH1F *time_hit_dif;
  // xml test
  std::map<std::string,std::string> m_parameters;

  std::vector<EVENT::RawCalorimeterHit*> _trigger_raw_hit;

  bool GAIN_CORRECTION_MODE;
  std::string _outFileName;
  std::string _noiseFileName;
  std::string _treeName;
  std::string _logrootName;
  std::string _colName;
  std::string _fileName;
  std::string _mappingfile;
  std::string _geomXML;
  std::ostream *_output;
  std::vector<std::string> _hcalCollections;
  int _overwrite;
  TTree *_outputTree;
  unsigned int _eventNr;
  Int_t _nHit;
  Int_t _elecNoiseCut;

  std::map<int, LayerID  > _mapping;
  std::map<int, double  > _chamberPos;//camber , pos

  int _noiseCut;
  int _timeWin;
  int _layerCut;
  int _time2prevEventCut;
  double _layerGap;
  float _beamEnergy;
  int _trigCount;
  int _maxTime;
  int evtnum;
  int _rejectedNum;
  uint _index[3];
  uintVec zcut;
  LCWriter* _lcWriter;
  int _bcid1;
  int _bcid2;

  std::string normal  ;
  std::string red     ;
  std::string green   ;
  std::string yellow  ;
  std::string blue    ;
  std::string magenta ;
  std::string white   ;

};



#endif


