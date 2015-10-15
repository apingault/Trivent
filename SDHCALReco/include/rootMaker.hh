#ifndef _rootMaker_hh_
#define _rootMaker_hh_


#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include "TRawCaloHitRoot.hh"
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "IO/LCWriter.h"



class rootMaker  : public marlin::Processor
{
public:

  Processor*  newProcessor() { return new rootMaker ; }
 
  rootMaker();
    
  ~rootMaker() {};

  void init();
     
  void processEvent( LCEvent * evtP );
  //void processRunHeader( LCRunHeader * runH);// added by me 

  void end();
    
protected:
  std::string _treeFileName;
  std::string _treeName;
  std::string _colName; 
  std::string _fileName;
  Int_t _elec_noise_cut;
  double _layer_gap;
  std::string _mappingfile;
  std::map<int, LayerID> _mapping;
  std::ostream *_output;
  std::vector<std::string> _hcalCollections;
  int _overwrite;
  TTree *_outputTree;  
  TClonesArray *_hitArray;
  TClonesArray *_rootArray;// added by me 
  unsigned int _eventNr;  
  unsigned int _triggerNr;  
  unsigned int _runNr;
  unsigned int _eventType;
  Int_t _Nhit; 
  Int_t _evt_id; 
  double _beamEnergy;
  
  int evtnum;

  LCWriter* _lcWriter;

};



#endif


