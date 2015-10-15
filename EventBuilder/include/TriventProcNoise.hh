/**
 * Yacine HADDAD
 * LLR Ecole polytechnique
 * avril 2012
 * Trivent v0.3
 */

#include <TriventProcNoise.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <limits.h>
#include <cmath>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>  
#include <Rtypes.h> 
#include <sstream>		
#include <UTIL/CellIDEncoder.h>
#include "Mapping.h"
#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include <fstream>
#include <algorithm>
#include <lcio.h>
#include "marlin/VerbosityLevels.h"
#include "marlin/tinyxml.h"

TriventProcNoise a_TriventProcNoise_instance;

//=========================================================
TriventProcNoise::TriventProcNoise()
  : Processor("TriventProcNoise"),
    _output(0),
    _outputTree(0)//,
			//_hitArray(0)
{
  
  streamlog_out( MESSAGE )<< "Trivent ... begin " << endl;
  _rejectedNum = 0;
  
  // collection 
  std::vector<std::string> hcalCollections;    
  hcalCollections.push_back(std::string("DHCALRawHits"));
  registerInputCollections( LCIO::RAWCALORIMETERHIT , 
			    "HCALCollections"       ,  
			    "HCAL Collection Names" ,  
			    _hcalCollections        , 
			    hcalCollections         ); 
  
  // Option of output file with clean events
  _outFileName="LCIO_clean_run.slcio";
  registerProcessorParameter("LCIOOutputFile" , 
			     "LCIO file" ,
			     _outFileName ,
			     _outFileName);
  // Energy
  _beamEnergy = 10;
  registerProcessorParameter("beamEnergy" ,
                             "The beam ",
                             _beamEnergy ,
                             _beamEnergy);
  // Option of output file with noise
  _noiseFileName="noise_run.slcio";
  registerProcessorParameter("NOISEutputFile" , 
			     "NOISE file" ,
			     _noiseFileName ,
			     _noiseFileName);
  
  // layer cut 
  _LayerCut = 10;
  registerProcessorParameter("LayerCut" ,
                             "cut in number of layer 10 in default",
                             _LayerCut ,
                             _LayerCut);


  // noise cut
  _noiseCut = 10;
  registerProcessorParameter("noiseCut" ,
                             "noise cut in time spectrum 10 in default",
                             _noiseCut ,
                             _noiseCut);

  // time windows
  _timeWin = 2;
  registerProcessorParameter("timeWin" ,
                             "time window = 2 in default",
                             _timeWin ,
                             _timeWin);
  //maping on XML file
  _geomXML = "setup_geometry.xml";
  registerProcessorParameter("setup_geometry" ,
                             "Dif geometry and position on the detector XML",
                             _geomXML,
                             _geomXML);
  
  //maping on txt file
  _mappingfile = "mapping_ps.txt";
  registerProcessorParameter("DIFMapping" ,
                             "dif's mapping file ",
                             _mappingfile,
                             _mappingfile);
  //inter layer
  _layer_gap = 9.;
  registerProcessorParameter("layerGap" ,
                             "Layers Gap in (cm)",
                             _layer_gap,
                             _layer_gap);
  // electronic noise cut
  _elec_noise_cut = 5000;
  registerProcessorParameter("electronic_noise_cut" ,
                             "number of hit max on time stamp",
                             _elec_noise_cut,
                             _elec_noise_cut);

  // electronic noise cut
  _time2prev_event_cut = 0;
  registerProcessorParameter("_time2prev_event_cut" ,
                             "cut on time to previous event (x 200 ns)",
                             _time2prev_event_cut,
                             _time2prev_event_cut);
  
  //log root file
  _treeName = "TEST";
  registerProcessorParameter("TreeName_logroot" ,
                             "Logroot tree name",
                             _treeName,
			     _treeName);
  // histogram control tree 
  _logrootName = "logroot.root";
  registerProcessorParameter("logroot_Name" ,
                             "Logroot name",
                             _logrootName,
			     _logrootName);
  
  GAIN_CORRECTION_MODE = false;
  registerProcessorParameter("GAIN_CORRECTION_MODE",
                             "GAIN_CORRECTION_MODE",
                             GAIN_CORRECTION_MODE,
			     GAIN_CORRECTION_MODE);
  
} 

void TriventProcNoise::XMLReader(std::string xmlfile){
  TiXmlDocument doc(xmlfile.c_str());
  bool load_key = doc.LoadFile();  
  if(load_key){
    streamlog_out( MESSAGE ) << green << "File : " << xmlfile.c_str() << normal <<std::endl;
    // tout ici 
    TiXmlHandle hDoc(&doc);
    TiXmlElement* pElem;
    TiXmlHandle hRoot(0);
    // name block
    {
      pElem=hDoc.FirstChildElement().Element();
      // should always have a valid root but handle gracefully if it does
      if (!pElem) streamlog_out( WARNING ) << red << "error elem" << normal << std::endl;
      streamlog_out( MESSAGE ) << green << pElem->Value() << normal << std::endl;
      
      // save this for later
      hRoot=TiXmlHandle(pElem);
    }
    // parameters block
    {
      m_parameters.clear();
      pElem=hRoot.FirstChild("parameter").Element();
      std::string key = pElem->Attribute("name");
      streamlog_out( MESSAGE ) << green << key.c_str() << normal << std::endl; 
      streamlog_out( DEBUG1 ) << green
			      <<"parameter : " 
			      << pElem->Attribute("name") 
			      << normal 
			      << std::endl;
    
      std::vector<std::string> lines;
      {
	std::string value = pElem->GetText() ;
	std::vector<std::string> lines;
	istringstream iss(value);
	copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter<vector<string> >(lines));
	for(unsigned int iline = 0; iline < lines.size(); iline++){
	  std::string line = lines.at(iline);
	  streamlog_out( MESSAGE ) << red << line << normal << std::endl;
	  
	  stringstream ss( line.c_str() );
	  vector<string> result;
	  
	  LayerID mapp;
	  int Dif_id;
	  while( ss.good() )
	    {
	      string substr;
	      getline( ss, substr, ',' );
	      result.push_back( substr );
	    }
	  istringstream ( result.at(0) ) >> Dif_id;
	  istringstream ( result.at(1) ) >> mapp.K;
	  istringstream ( result.at(2) ) >> mapp.DifX;
	  istringstream ( result.at(3) ) >> mapp.DifY;
	  istringstream ( result.at(4) ) >> mapp.IncX;
	  istringstream ( result.at(5) ) >> mapp.IncY;
	  _mapping[Dif_id] = mapp;
	}
      }
      pElem = pElem->NextSiblingElement();
      // ChamberGeom  Node.
      {
	streamlog_out( DEBUG1 ) << green
				<<"parameter : " 
				<< pElem->Attribute("name") 
				<< normal 
				<< std::endl;
	std::vector<std::string> lines;
	{
	  std::string value = pElem->GetText() ;
	  std::vector<std::string> lines;
	  istringstream iss(value);
	  copy(istream_iterator<string>(iss),
	       istream_iterator<string>(),
	       back_inserter<vector<string> >(lines));
	  for(unsigned int iline = 0; iline < lines.size(); iline++){
	    std::string line = lines.at(iline);
	    streamlog_out( MESSAGE ) << red << line << normal << std::endl;
	    
	    stringstream ss( line.c_str() );
	    vector<string> result;
	    
	    double position;
	    int Dif_id;
	    while( ss.good() )
	      {
		string substr;
		getline( ss, substr, ',' );
		result.push_back( substr );
	      }
	    istringstream ( result.at(0) ) >> Dif_id;
	    istringstream ( result.at(3) ) >> position;
	    
	    _chamber_pos[Dif_id] = position;
	  }
	}
      }
    }
  }else{
    streamlog_out( WARNING ) << red << "Failed to load file : " << xmlfile.c_str() << normal <<std::endl;
  }
}

void TriventProcNoise::readDifGeomFile(std::string geomfile){
  
  cout << "read the mapping file .."<< endl;
  
  LayerID contenu;
  ifstream file(geomfile.c_str(), ios::in);
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
  else
    cerr << "ERROR ... maping file not correct !" << endl;
}

void TriventProcNoise::printDifGeom(){
  
  for(std::map<int,LayerID>::iterator itt = _mapping.begin();itt!=_mapping.end();itt++)     {
    streamlog_out( MESSAGE ) << itt->first << "\t" << itt->second.K 
			     <<"\t"<<itt->second.DifX 
			     <<"\t"<<itt->second.DifY
			     <<"\t"<<itt->second.IncX
			     <<"\t"<<itt->second.IncY
			     << std::endl;
  }
}

// ============ decode the cell ids =============
uint TriventProcNoise::getCellDif_id(int cell_id){
  return cell_id & 0xFF;
}
uint TriventProcNoise::getCellAsic_id(int cell_id){
  return (cell_id & 0xFF00)>>8;
}
uint TriventProcNoise::getCellChan_id(int cell_id){
  return (cell_id & 0x3F0000)>>16;
}

uint* TriventProcNoise::getPadIndex(uint dif_id, uint asic_id, uint chan_id){
  _index[0]=_index[1]=_index[2]=0;
  double DifY = -1.,DifZ = -1.;
  DifZ = _mapping.find(dif_id)->second.K;
  DifY = _mapping.find(dif_id)->second.DifY;
  _index[0] = (1+MapILargeHR2[chan_id]+AsicShiftI[asic_id]);
  _index[1] = (32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+int(DifY);
  _index[2] = abs(int(DifZ));
  streamlog_out( DEBUG0 ) << " Dif_id == " << dif_id
			  << " Asic_id ==" << asic_id
			  << " Chan_id ==" << chan_id
			  << " I == " << _index[0]
			  << " J == " << _index[1]
			  << " K == " << _index[2]
			  << std::endl;
  return _index;
}
//===============================================
void TriventProcNoise::getMaxTime()
{
  _maxtime = 0;
  try{
    for(std::vector<EVENT::RawCalorimeterHit*>::iterator raw_hit=_trigger_raw_hit.begin();raw_hit!= _trigger_raw_hit.end();raw_hit++){
      int time =  int((*raw_hit)->getTimeStamp());
      if(time >= 0) _maxtime = max(_maxtime, time);
    }
  }catch (std::exception ec){
    streamlog_out( WARNING )<<"No hits "<<std::endl;
  }
  streamlog_out( DEBUG1 ) << " maxtime before == " << _maxtime << std::endl;
  //return maxtime;
}


std::vector<int> TriventProcNoise::getTimeSpectrum() //__attribute__((optimize(0)))
{
  std::vector<int> time_spectrum(_maxtime+1);
  try{
    for(std::vector<EVENT::RawCalorimeterHit*>::iterator raw_hit=_trigger_raw_hit.begin();raw_hit!= _trigger_raw_hit.end();raw_hit++){  
      int time =  int((*raw_hit)->getTimeStamp());
      if(time >= 0) time_spectrum[time]++;
    }
  }catch (std::exception ec){
    streamlog_out( WARNING )<<"No hits "<<std::endl;
  }
  return time_spectrum;
}

bool TriventProcNoise::peakOrNot(std::vector<int> time_spectrum ,int itime ,int threshold){
  
#if HISTOGRAM_PARSER
  noise_dist->Fill(time_spectrum[itime]);
#endif
  
  if(time_spectrum[itime] >= threshold
     && time_spectrum[itime] >  time_spectrum[itime+1]  
     && time_spectrum[itime] >= time_spectrum[itime+1]){
    return true;
  }else{
    return false;
  }
}

int IJKToKey(const int i,const int j,const int k){return 100*100*k+100*j+i;}
void TriventProcNoise::eventBuilder(LCCollection* col_event,int time_peak, int prev_time_peak){
  zcut.clear();
  col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
  col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
  CellIDEncoder<CalorimeterHitImpl> cd( "M:3,S-1:3,I:9,J:9,K-1:6" ,col_event) ;
  try{
    std::vector<int> hitKeys;
    for(std::vector<EVENT::RawCalorimeterHit*>::iterator rawhit=_trigger_raw_hit.begin();rawhit!= _trigger_raw_hit.end();rawhit++){  
      float pos[3];
      int time = (*rawhit)->getTimeStamp();
      if(fabs(time-time_peak) <=_timeWin && 
	 (time > prev_time_peak + _timeWin )){
	  
	int Dif_id  =  getCellDif_id ((*rawhit)->getCellID0());
	int Asic_id =  getCellAsic_id((*rawhit)->getCellID0());
	int Chan_id =  getCellChan_id((*rawhit)->getCellID0());
	  
	uint I = getPadIndex(Dif_id, Asic_id, Chan_id)[0];
	uint J = getPadIndex(Dif_id, Asic_id, Chan_id)[1];
	uint K = getPadIndex(Dif_id, Asic_id, Chan_id)[2];
	int aHitKey=IJKToKey(I,J,K);
	pos[0] = I*10.*1.0408;
	pos[1] = J*10.*1.0408;
	pos[2] = K*26.131;

	if(K<=0||K>64) {streamlog_out( DEBUG ) << Dif_id << std::endl; continue;}
	CalorimeterHitImpl* caloHit = new CalorimeterHitImpl();
	caloHit->setTime(float((*rawhit)->getTimeStamp())); // done !!
	
	if(float((*rawhit)->getAmplitude()&3)>2.5) caloHit->setEnergy(float((*rawhit)->getAmplitude()&3));
	else if(float((*rawhit)->getAmplitude()&3)>1.5) caloHit->setEnergy(float((*rawhit)->getAmplitude()&3)-1);
	else caloHit->setEnergy(float((*rawhit)->getAmplitude()&3)+1);
	
	//avoid two hit in the same cell
	if(std::find(hitKeys.begin(),hitKeys.end(),aHitKey)!=hitKeys.end()){
	  IMPL::CalorimeterHitImpl* hit=
	    dynamic_cast<IMPL::CalorimeterHitImpl*>(col_event->getElementAt(std::distance(hitKeys.begin(),std::find(hitKeys.begin(),hitKeys.end(),aHitKey))));
	  float hitTime=hit->getTime();
	  if( fabs(time_peak-hitTime)>fabs(time_peak-time) ){
	    hit->setEnergy(caloHit->getEnergy());
	  }
	  continue;
	}
	// set the cell id 
	cd["I"] = I ;
	cd["J"] = J ;
	cd["K-1"] = K-1 ;
	cd["M"] = 0 ;
	cd["S-1"] = 3 ;
	streamlog_out( DEBUG0 ) << " I == " << I
				<< " J == " << J
				<< " K == " << K
				<< std::endl;
	cd.setCellID( caloHit ) ;
	if(std::find(zcut.begin(), zcut.end(), K)==zcut.end())
	  zcut.push_back(K);
	caloHit->setPosition(pos);
	col_event->addElement(caloHit);
	hitKeys.push_back(aHitKey);
      }
      //}else{
      //streamlog_out( MESSAGE ) << " time peak = "<<time_peak<<" pointer --> : " << rawhit << std::endl;
      //}
    }//loop over the hit     
    hitKeys.clear();
  }catch(DataNotAvailableException &e){
    streamlog_out(WARNING) << " collection not available" << std::endl;
  }
}

//===============================================
void TriventProcNoise::init() {
  trig_count = 0;
  //========================
  //readDifGeomFile(_mappingfile.c_str());
  
  // ========================
  
  printParameters();
  // new process
  
  char cnormal[8] =  {0x1b,'[','0',';','3','9','m',0};
  char cred[8]     = {0x1b,'[','1',';','3','1','m',0};
  char cgreen[8]   = {0x1b,'[','1',';','3','2','m',0};
  char cyellow[8]  = {0x1b,'[','1',';','3','3','m',0};
  char cblue[8]    = {0x1b,'[','1',';','3','4','m',0};
  char cmagenta[8] = {0x1b,'[','1',';','3','5','m',0};
  char cwhite[8]   = {0x1b,'[','1',';','3','9','m',0};

  normal   = cnormal;
  red      = cred;
  green    = cgreen;
  yellow   = cyellow;
  blue     = cblue;
  magenta  = cmagenta;
  white    = cwhite;

  _lcWriter = LCFactory::getInstance()->createLCWriter() ;
  _lcWriter->setCompressionLevel( 0 ) ;
  _lcWriter->open(_outFileName.c_str(),LCIO::WRITE_NEW) ; 
  
 
  XMLReader(_geomXML.c_str());
  printDifGeom();
  evtnum=0;// event number

}
//==================================================================================
//void TriventProcNoise::setTriggerRawHit() {
  
//}
//==================================================================================
void TriventProcNoise::processRunHeader( LCRunHeader * runHd ) {
  
}

//==================================================================================
void TriventProcNoise::processEvent( LCEvent * evtP ) 
{	
  if (evtP != NULL){
    try{
      
      _eventNr=evtP->getEventNumber();
      for(unsigned int i=0; i< _hcalCollections.size(); i++){//!loop over collection
	try{
	  
	  LCCollection * col = NULL;
	  col = evtP ->getCollection(_hcalCollections[i].c_str());
	  int numElements = col->getNumberOfElements();// hit number 
	  
	  
	  streamlog_out( MESSAGE ) << yellow << "Trigger number == " << trig_count++ << normal << std::endl;
	  
	  if(col == NULL )  {
	    streamlog_out( WARNING )<< red << "TRIGGER SKIPED ..."<< normal <<std::endl;
	    break;
	  }
	  
	  if(numElements > _elec_noise_cut)  {
	    streamlog_out( MESSAGE ) << red << "TRIGGER SKIPED ..."<< normal <<std::endl;
	    break;
	  }
	  
	  // set raw hits 
	  _trigger_raw_hit.clear();
	  std::vector<int> vTrigger;
	  for (int ihit(0); ihit < col->getNumberOfElements(); ++ihit) {// loop over the hits
	    RawCalorimeterHit *raw_hit = 
	      dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit)) ;
	    if (NULL != raw_hit){
	      _trigger_raw_hit.push_back(raw_hit);
	      //extract abolute bcid information: 
	      if(ihit==0){
		unsigned int difid=0; 
		difid = raw_hit->getCellID0()&0xFF;
		if (difid==0) return;
		std::stringstream pname("");
		pname <<"DIF"<<difid<<"_Triggers";
		col->getParameters().getIntVals(pname.str(),vTrigger);
		if (vTrigger.size()!=0){
		  _bcid1=vTrigger[4] ;
		  _bcid2=vTrigger[3] ;
		  unsigned long long Shift=16777216ULL;//to shift the value from the 24 first bits
		  unsigned long long theBCID_=_bcid1*Shift+_bcid2; 
		  streamlog_out( DEBUG1 ) << "trigger time : " << theBCID_ << std::endl;
		}
	      }
	    }
	  }
	  getMaxTime();
	  std::vector<int> time_spectrum = getTimeSpectrum();
	  
	  //---------------------------------------------------------------
	  //! Find the condidate event
	  int ibin=0;
	  int bin_c_prev = -2 * _timeWin; //  the previous bin center
	  
	  int time_prev = 0;
	  while(ibin < (_maxtime+1)){ 
	    if(time_spectrum[ibin] < _noiseCut && 
	       time_spectrum[ibin] >  time_spectrum[ibin+1] && 
	       time_spectrum[ibin] >= time_spectrum[ibin+1] &&
	       time_spectrum[ibin] >= time_spectrum[ibin+1] && 
	       time_spectrum[ibin] >= time_spectrum[ibin+2] ){
	      LCEventImpl*  evt = new LCEventImpl() ;     // create the event
	      
	      //---------- set event paramters ------
	      const std::string parname_trigger = "trigger";
	      const std::string parname_energy  = "beamEnergy";
	      const std::string parname_bcid1 = "bcid1";
	      const std::string parname_bcid2 = "bcid2";
	      evt->parameters().setValue(parname_trigger,evtP->getEventNumber()); 
	      evt->parameters().setValue(parname_energy , _beamEnergy);
	      evt->parameters().setValue(parname_bcid1 , _bcid1);
	      evt->parameters().setValue(parname_bcid2 , _bcid2);
	      evt->setRunNumber( evtP->getRunNumber()) ;
	      //-------------------------------------
	      
	      LCCollectionVec* outcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
	      TriventProcNoise::eventBuilder(outcol,ibin,bin_c_prev);
	      
	      streamlog_out( DEBUG1 ) << "zcut.size() = " << zcut.size() << "\t _LayerCut = " << _LayerCut << std::endl;
	      if( (int)zcut.size() >_LayerCut && abs(int(ibin)-time_prev) > _time2prev_event_cut){ 
		streamlog_out( DEBUG5 ) <<green<<" Trivent find event at :==> "<< red << ibin 
					<<green<<"\t :Nhit: ==> "<< magenta
					<<outcol->getNumberOfElements() << normal <<std::endl;  		
		evt->setEventNumber( evtnum++ ) ;
		evt->addCollection(outcol, "SDHCAL_HIT");
		_lcWriter->writeEvent( evt ) ;
	      }else{
		_rejectedNum++;
		delete outcol;
	      }
	      time_prev = ibin;
	      delete evt; evt=NULL;
	      
	      bin_c_prev = ibin;
	      ibin = ibin+_timeWin;
	    }else{ibin++;}
	  }
	  
	}catch (lcio::DataNotAvailableException zero) {}
      }
    }catch (lcio::DataNotAvailableException err) {}
  }
  
	   
}	
//==============================================================
void TriventProcNoise::end()
{       
  streamlog_out( MESSAGE )<< "Trivent Select "<<_rejectedNum <<" Condidate event"<<std::endl;
  streamlog_out( MESSAGE )<< "Trivent end"<<std::endl;
  //cc.StoreHistos("test.root");
  _lcWriter->close();
  
  if (_outputTree) {
    TFile *_logroot = _outputTree->GetCurrentFile();
    _logroot->Write();
    delete _logroot;
  }
}
//==============================================================


 
