#include <rootMaker.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <EVENT/RawCalorimeterHit.h>
#include <limits.h>
#include <string>
#include <iostream>
#include <cmath>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>  
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TH2.h>//we can use the root histograms 
#include <TH3.h>
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include <fstream>
using namespace std;
//#include "TCalorimeterHit.hh" 
#include "TCaloHitRoot.hh" // the hit class used to build events (TObject)
//#include "RAWReader.hh"
rootMaker a_rootMaker_instance;

//=========================================================
rootMaker::rootMaker():Processor("rootMaker"),
		       _output(0),
		       _outputTree(0),
		       _rootArray(0),
		       _runNr(0),
		       _eventType(0),
		       _evt_id(0){
  
  
  _description = "DHCAL root file maker" ;
  
  _treeFileName="rootmake_test.root";
  registerProcessorParameter("TreeOutputFile" , 
			     "ROOT tree" ,
			     _treeFileName ,
			     _treeFileName);
  
  _rootArray = new TClonesArray("TCaloHitRoot");
  _colName="DHCAL hits";
  registerProcessorParameter( "MCObjects" ,
			      "The name of the PFOs" ,
			      _colName ,
			      _colName);
  
  std::vector<std::string> hcalCollections;    
  //hcalCollections.push_back(std::string(""));
  hcalCollections.push_back(std::string("SDHCAL_HIT"));
  
  //hcalCollections.push_back(std::string("HcalBarrelRegCollection"));
  //hcalCollections.push_back(std::string("HcalEndcapRingsCollection"));
  //hcalCollections.push_back(std::string("HcalEndcapsCollection"));

  registerInputCollections( LCIO::CALORIMETERHIT,
			    "HCALCollections" ,
			    "HCAL Collection Names" ,
			    _hcalCollections ,
			    hcalCollections);
  
  _treeName="DHCAL";
  registerProcessorParameter( "TreeName" , 
			      "The name of the ROOT tree" ,
			      _treeName ,
			      _treeName);
	
  _overwrite=0;
  registerProcessorParameter( "OverwriteFile" , 
			      "overwrit a file if zero." ,
			      _overwrite ,
			      _overwrite);
  
  _layer_gap = 2.8;
  registerProcessorParameter("layerGap" ,
                             "Layers Gap in (cm)",
                             _layer_gap,
                             _layer_gap);

  _elec_noise_cut = 5000;
  registerProcessorParameter("electronic_noise_cut" ,
                             "number of hit max on time stamp",
                             _elec_noise_cut,
                             _elec_noise_cut);
  
  //maping 
  _mappingfile = "mapping_ps.txt";
  registerProcessorParameter("DIFMapping" ,
                             "dif's mapping file ",
                             _mappingfile,
                             _mappingfile);
  }

//=========================================================
void rootMaker::init() {
  
  printParameters();
  // read the 
  
  std::cout << "read the mapping file .."<< std::endl;
  
  //print the maping 
  LayerID contenu;
  ifstream file(_mappingfile.c_str(), ios::in);
  if(file){ 
    while(!file.eof()){
      int Dif_id;
      char co;
      file >> Dif_id >> co 
	   >> contenu.K >> co 
	   >> contenu.DifX >> co 
	   >> contenu.DifY >> co  
	   >> contenu.IncX >> co  
	   >> contenu.IncY ;
      _mapping [Dif_id] = contenu;
    }
    file.close();
  }
  else  // sinon
    cerr << "ERROR ... maping file not correct !" << endl;
  
  
  for(std::map<int,LayerID>::iterator itt = _mapping.begin();itt!=_mapping.end();itt++)     {
    std::cout << itt->first << "\t" << itt->second.K 
	      <<"\t"<<itt->second.DifX 
	      <<"\t"<<itt->second.DifY
	      <<"\t"<<itt->second.IncX
      	      <<"\t"<<itt->second.IncY
	      << std::endl;
  }
  


  TFile *tree_file=new TFile(_treeFileName.c_str(),"RECREATE");
  
  if (!tree_file->IsOpen()) {// test the fil opening 
    delete tree_file;
    tree_file=new TFile(_treeFileName.c_str(),"NEW");
  }else{
    std::cout << "mahada !!!!" << std::endl;
  }
  
  UInt_t bufsize=64000;
  
  //create the branchies 
  _outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
  _outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB
  _outputTree->Branch("runNr",&_runNr,"runNr/I");
  _outputTree->Branch("eventNr",&_eventNr,"eventNr/I");
  _outputTree->Branch("triggerNr",&_triggerNr,"triggerNr/I");
  _outputTree->Branch("eventType",&_eventType,"eventType/I");
  _outputTree->Branch("beamEnergy",&_beamEnergy,"beamEnergy/D");
  _outputTree->Branch("Nhit",&_Nhit,"Nhit/I");
  _outputTree->Branch("DHCALEvent",&_rootArray,bufsize,99);// the same variable 
  
}
//=========================================================
void rootMaker::processEvent( LCEvent * evtP ){	
  
  _rootArray->Clear("C");
  if (evtP){	
    try{  
      _eventNr=evtP->getEventNumber();
      _runNr=evtP->getRunNumber();
      // get the LCEvent paramters ===
      std::string trigger_key = "trigger";
      const std::string energy_key = "beamEnergy";
      _triggerNr = evtP->getParameters().getIntVal(trigger_key);
      _beamEnergy = evtP->getParameters().getIntVal(energy_key);
      std::cout << "Trigger number == " << _triggerNr << std::endl;
      
      // ==============================
      
      _eventType = 0;
      _Nhit=0;
      for(unsigned int i=0; i< _hcalCollections.size(); i++){
	try{
	  LCCollection *col = evtP ->getCollection(_hcalCollections[i].c_str());
	  std::string coltype = col->getTypeName();
	  std::string str_1 ("RawCalorimeterHit");
	  int numElements = col->getNumberOfElements();
	  
	  std::cout << "The collection type is == "
		    << coltype.c_str() 
		    << "  with "
		    << numElements 
		    << " elements"
		    << std::endl;
	  
	  
	  
	  // COLLECTION OF CALORIMETER HITS
	  if( coltype.compare("CalorimeterHit") == 0){
	    CellIDDecoder<CalorimeterHit> idDecoder(col);
	    //int numElements = col->getNumberOfElements();
	    if(numElements < 100000){
	      for (int j(0); j < numElements; ++j) {
		CalorimeterHit *a_DHcalhit = 
		  dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
		//std::cout<<"Time == " <<a_DHcalhit->getTime()<<std::endl;
		
		_Nhit=numElements;
		new((*_rootArray)[j]) TCaloHitRoot(*a_DHcalhit,&idDecoder);
	      }
	    }
	  }
	  
	  // COLLECTION FO RAWCALORIMETERHIT
	  if(coltype.compare("RawCalorimeterHit") == 0){
	    std::cout<<"Trigger with == " 
		     << numElements
		     <<std::endl;// numElements << " Hits ..."<<std::endl;
	    
	    if(numElements < 100000){
	      for (int j(0); j < numElements; ++j) {
		RawCalorimeterHit *a_DHcalhit = 
		  dynamic_cast<RawCalorimeterHit*>( col->getElementAt( j ) ) ;
		
		_Nhit=numElements;
		//SetLayerGap(_layer_gap);
		new((*_rootArray)[j]) TRawCaloHitRoot(*a_DHcalhit,_mapping);
		
	      }
	    }
	  }
	  _outputTree->Fill();
      }catch (lcio::DataNotAvailableException zero) { }
    }
    }	
    catch (lcio::DataNotAvailableException err) { }
  }
  
}	
//=========================================================
void rootMaker::end(){
  
  if (_outputTree) {
    TFile *tree_file = _outputTree->GetCurrentFile(); 
    std::cout<<"test end"<<std::endl;
    tree_file->Write();
    delete tree_file;
  }
}
//=========================================================


 
