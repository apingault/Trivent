/**
 * Yacine HADDAD
 * LLR Ecole polytechnique
 * avril 2012
 * Trivent v0.3
 */

#include <TriventProc.hh>

// -- std includes
#include <fstream>    // std::stringstream
#include <iterator>   // std::next

// -- Root headers
#include <TCanvas.h>


TriventProc a_TriventProc_instance;

//=========================================================
TriventProc::TriventProc()
  : Processor("TriventProc"),
    _outputTree(0),
    m_cellSizeI(10.408),
    m_cellSizeJ(10.408),
    m_layerThickness(26.131)
{

  streamlog_out( MESSAGE ) << "Trivent ... begin " << endl;
  _rejectedNum = 0;
  _selectedNum = 0;

  // collection
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("DHCALRawHits"));
  registerInputCollections( LCIO::RAWCALORIMETERHIT ,
                            "HCALCollections"       ,
                            "HCAL Collection Names" ,
                            _hcalCollections        ,
                            hcalCollections         );

  // Option of output file with clean events
  _outFileName = "LCIO_clean_run.slcio";
  registerProcessorParameter("LCIOOutputFile" ,
                             "LCIO file" ,
                             _outFileName ,
                             _outFileName);
  // Energy
  _beamEnergy = 0;
  registerProcessorParameter("beamEnergy" ,
                             "The beam ",
                             _beamEnergy ,
                             _beamEnergy);

  // Option of output file with noise
  _noiseFileName = "noise_run.slcio";
  registerProcessorParameter("NOISEutputFile" ,
                             "NOISE file" ,
                             _noiseFileName ,
                             _noiseFileName);

  // layer cut
  _layerCut = 10;
  registerProcessorParameter("LayerCut" ,
                             "cut in number of layer 10 in default",
                             _layerCut ,
                             _layerCut);

  // noise cut
  _noiseCut = 10;
  registerProcessorParameter("NoiseCut" ,
                             "noise cut in time spectrum 10 in default",
                             _noiseCut ,
                             _noiseCut);

  // time windows
  _timeWin = 2;
  registerProcessorParameter("TimeWin" ,
                             "time window = 2 in default",
                             _timeWin ,
                             _timeWin);
  //maping on XML file
  _geomXML = "setup_geometry.xml";
  registerProcessorParameter("SetupGeometry" ,
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
  _layerGap = 9.;
  registerProcessorParameter("LayerGap" ,
                             "Layers Gap in (cm)",
                             _layerGap,
                             _layerGap);
  // electronic noise cut
  _elecNoiseCut = 5000;
  registerProcessorParameter("ElectronicNoiseCut" ,
                             "number of hit max on time stamp",
                             _elecNoiseCut,
                             _elecNoiseCut);

  // electronic noise cut
  _time2prevEventCut = 0;
  registerProcessorParameter("_time2prev_event_cut" ,
                             "cut on time to previous event (x 200 ns)",
                             _time2prevEventCut,
                             _time2prevEventCut);

  //log root file
  _treeName = "TEST";
  registerProcessorParameter("TreeName_logroot" ,
                             "Logroot tree name",
                             _treeName,
                             _treeName);
  // histogram control tree
  _rootFileName = "logroot.root";
  registerProcessorParameter("ROOTOutputFile" ,
                             "Logroot name",
                             _rootFileName,
                             _rootFileName);

  GAIN_CORRECTION_MODE = false;
  registerProcessorParameter("GAIN_CORRECTION_MODE",
                             "GAIN_CORRECTION_MODE",
                             GAIN_CORRECTION_MODE,
                             GAIN_CORRECTION_MODE);

}

//=============================================================================
void TriventProc::XMLReader(std::string xmlfile) {
  TiXmlDocument doc(xmlfile.c_str());
  bool load_key = doc.LoadFile();
  if (load_key) {
    streamlog_out( MESSAGE ) << green << "File : " << xmlfile.c_str() << normal << std::endl;
    // tout ici
    TiXmlHandle hDoc(&doc);
    TiXmlElement* pElem;
    TiXmlHandle hRoot(0);
    // name block
    {
      pElem = hDoc.FirstChildElement().Element();
      // should always have a valid root but handle gracefully if it does
      if (!pElem) streamlog_out( WARNING ) << red << "error elem" << normal << std::endl;
      streamlog_out( MESSAGE ) << green << pElem->Value() << normal << std::endl;

      // save this for later
      hRoot = TiXmlHandle(pElem);
    }
    // parameters block
    {
      m_parameters.clear();
      pElem = hRoot.FirstChild("parameter").Element();
      std::string key = pElem->Attribute("name");
      streamlog_out( MESSAGE ) << green << key.c_str() << normal << std::endl;
      streamlog_out( DEBUG1 ) << green
                              << "parameter : "
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
        for (unsigned int iline = 0; iline < lines.size(); iline++) {
          std::string line = lines.at(iline);
          streamlog_out( MESSAGE ) << red << line << normal << std::endl;

          stringstream ss( line.c_str() );
          vector<string> result;

          LayerID mapp;
          int Dif_id;
          while ( ss.good() )
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
                                << "parameter : "
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
          for (unsigned int iline = 0; iline < lines.size(); iline++) {
            std::string line = lines.at(iline);
            streamlog_out( MESSAGE ) << red << line << normal << std::endl;

            stringstream ss( line.c_str() );
            vector<string> result;

            double position;
            int Dif_id;
            while ( ss.good() )
            {
              string substr;
              getline( ss, substr, ',' );
              result.push_back( substr );
            }
            istringstream ( result.at(0) ) >> Dif_id;
            istringstream ( result.at(3) ) >> position;

            _chamberPos[Dif_id] = position;
          }
        }
      }
    }
  } else {
    streamlog_out( WARNING ) << red << "Failed to load file : " << xmlfile.c_str() << normal << std::endl;
  }
}

//=============================================================================
void TriventProc::readDifGeomFile(std::string geomfile) {

  cout << "read the mapping file .." << endl;

  LayerID contenu;
  ifstream file(geomfile.c_str(), ios::in);
  if (file) {
    while (!file.eof()) {
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

//=============================================================================
void TriventProc::printDifGeom() {

  for (std::map<int, LayerID>::iterator itt = _mapping.begin(); itt != _mapping.end(); itt++)     {
    streamlog_out( MESSAGE ) << itt->first << "\t" << itt->second.K
                             << "\t" << itt->second.DifX
                             << "\t" << itt->second.DifY
                             << "\t" << itt->second.IncX
                             << "\t" << itt->second.IncY
                             << std::endl;
  }
}

// ============ decode the cell ids =============
// bit shift & 0xFF = Apply mask 1111 1111 to binary value
// eg: Dif 1 => cellID0 = 00983297 => DifID = 1 / AsicID = 1 / ChanID = 15
uint TriventProc::getCellDif_id(int cell_id) {
  return cell_id & 0xFF;
}

//=============================================================================
//  bit shift & 0xFF00 Apply mask 1111 1111 0000 0000 then cut last 8 bits
uint TriventProc::getCellAsic_id(int cell_id) {
  return (cell_id & 0xFF00) >> 8;
}

//=============================================================================
//  bit shift & 0x3F0000 Apply mask 1111 0000 0000 0000 0000 then cut last 16 bits
uint TriventProc::getCellChan_id(int cell_id) {
  return (cell_id & 0x3F0000) >> 16;
}

// ============ ============ ============ ============ ============ ============ ============
// ============ ============ ============ ============ ============ ============ ============
// Full example with Dif 1:
// cellId0 = 00983297 -> Binaire =  1111 0000 0001 0000 0001
// binaire & 0xFF = 0000 0001 => 2^0 = 1
// binaire & 0xFF00 = 0000 0001 0000 0000 >>8 = 0000 0001 => 2^0 = 1
// binaire & 0x3F0000 =  1111 0000 0000 0000 0000 >>16 = 1111 => (2^3)+(2^2)+(2^1)+(2^0) = 15
// ============ ============ ============ ============ ============ ============ ============
// ============ ============ ============ ============ ============ ============ ============

//=============================================================================
std::vector<unsigned int> TriventProc::getPadIndex(unsigned int dif_id, unsigned int asic_id, unsigned int chan_id) {
  std::vector<unsigned int> index(3, 0);
  std::map<int, LayerID>::const_iterator findIter = _mapping.find(dif_id);

  index[0] = (1 + MapILargeHR2[chan_id] + AsicShiftI[asic_id]);
  index[1] = ( 32 - (MapJLargeHR2[chan_id] + AsicShiftJ[asic_id]) ) + findIter->second.DifY;
  index[2] = findIter->second.K;

  streamlog_out( DEBUG0 ) << " Dif_id == " << dif_id
                          << " Asic_id ==" << asic_id
                          << " Chan_id ==" << chan_id
                          << " I == " << index[0]
                          << " J == " << index[1]
                          << " K == " << index[2]
                          << std::endl;
  return index;
}

//=============================================================================
void TriventProc::getMaxTime()
{
  _maxTime = 0;
  try {
    for (std::vector<EVENT::RawCalorimeterHit*>::const_iterator raw_hit = _trigger_raw_hit.begin(); raw_hit != _trigger_raw_hit.end(); raw_hit++) {
      int time =  static_cast<int>((*raw_hit)->getTimeStamp());

      if (time >= 0)
        _maxTime = max(_maxTime, time);
    }
  } catch (std::exception ec) {
    streamlog_out( WARNING ) << "No hits " << std::endl;
  }
}

//=============================================================================
std::vector<int> TriventProc::getTimeSpectrum() //__attribute__((optimize(0)))
{
  std::vector<int> time_spectrum(_maxTime + 1);
  try {
    for (std::vector<EVENT::RawCalorimeterHit*>::const_iterator raw_hit = _trigger_raw_hit.begin(); raw_hit != _trigger_raw_hit.end(); raw_hit++) {
      int time =  static_cast<int>((*raw_hit)->getTimeStamp());
      if (time > _maxTime)
      {
        streamlog_out( WARNING ) << "\t *** WARNING *** Found Hit after _maxTime -> hitTime: " << time << " / maxTime: " << _maxTime << std::endl;
        continue;
      }
      if (time >= 0)
        ++time_spectrum.at(time);
    }
  } catch (std::exception ec) {
    streamlog_out( WARNING ) << "No hits " << std::endl;
  }
  return time_spectrum;
}

//=============================================================================
int IJKToKey(const int i, const int j, const int k) {return 100 * 100 * k + 100 * j + i;}

//=============================================================================
int findAsicKey(int i, int j, int k)
{
  if (i > 96 || i < 0 || j > 96 || j < 0) return -1;
  int jnum = (j - 1) / 8;
  int inum = (i - 1) / 8;
  int num = jnum * 12 + inum;
  return k * 1000 + num;
}

//=============================================================================
void TriventProc::eventBuilder(LCCollection* col_event, int time_peak, int prev_time_peak) {

  // reset Number of layers touched
  _firedLayersSet.clear();

  col_event->setFlag(col_event->getFlag() | ( 1 << LCIO::RCHBIT_LONG));
  col_event->setFlag(col_event->getFlag() | ( 1 << LCIO::RCHBIT_TIME));

  CellIDEncoder<CalorimeterHitImpl> cd( "M:3,S-1:3,I:9,J:9,K-1:6" , col_event);

  std::map<int, int> asicMap;
  try {
    std::vector<int> hitKeys;
    for (std::vector<EVENT::RawCalorimeterHit*>::const_iterator rawhit = _trigger_raw_hit.begin(); rawhit != _trigger_raw_hit.end(); rawhit++) {
      int time = static_cast<int>((*rawhit)->getTimeStamp());
      if (std::fabs(time - time_peak) <= _timeWin &&
          (time > prev_time_peak + _timeWin )) {

        int Dif_id  =  getCellDif_id ((*rawhit)->getCellID0());
        int Asic_id =  getCellAsic_id((*rawhit)->getCellID0());
        int Chan_id =  getCellChan_id((*rawhit)->getCellID0());

        std::vector<unsigned int> padIndex = getPadIndex(Dif_id, Asic_id, Chan_id);

        unsigned int I = padIndex[0];
        unsigned int J = padIndex[1];
        unsigned int K = padIndex[2];

        if (K <= 0 || K > 64) {
          streamlog_out( WARNING ) << " Found hit in Layer '" << K << "' DifId: " << Dif_id << std::endl;
          continue;
        }

        //find and remove square events
        int asickey = findAsicKey(I, J, K);
        if (asicMap[asickey])
          ++asicMap[asickey];
        else
          asicMap[asickey] = 1;

        if ( asicMap[asickey] == 64 ) {
          streamlog_out ( MESSAGE) << "Rejecting event with full asic. Dif " << Dif_id << " asic " << Asic_id << "' ... " << std::endl;

          _firedLayersSet.clear();
          hitKeys.clear();
          asicMap.clear();
          return;
        }


        // Creating Calorimeter Hit
        float pos[3];
        pos[0] = I * m_cellSizeI;
        pos[1] = J * m_cellSizeJ;
        pos[2] = K * m_layerThickness;

        CalorimeterHitImpl* caloHit = new CalorimeterHitImpl();
        caloHit->setTime(static_cast<float>((*rawhit)->getTimeStamp()));

        float hitShiftedAmplitude = static_cast<float>((*rawhit)->getAmplitude() & 3);
        if (hitShiftedAmplitude > 2.5)
          caloHit->setEnergy(hitShiftedAmplitude);        // 3rd treshold
        else if (hitShiftedAmplitude > 1.5)
          caloHit->setEnergy(hitShiftedAmplitude - 1);    // 2nd treshold ?
        else
          caloHit->setEnergy(hitShiftedAmplitude + 1);    // 1st treshold ?

        // Create hit Key
        int aHitKey = IJKToKey(I, J, K);

        // Avoid two hit in the same cell
        std::vector<int>::iterator findIter = std::find(hitKeys.begin(), hitKeys.end(), aHitKey);

        if (findIter != hitKeys.end())
        {
          streamlog_out( WARNING ) << " Found two hits in the same cell! " << std::endl;
          delete caloHit;
          continue;
        }


        // set the cell id
        cd["I"] = I ;
        cd["J"] = J ;
        cd["K-1"] = K - 1 ;
        cd["M"] = 0 ;
        cd["S-1"] = 3 ;


        streamlog_out( DEBUG ) << " I == " << I
                               << " J == " << J
                               << " K == " << K
                               << std::endl;

        // Fill hitsMap for each Layer, starting at K-1 = 0
        streamlog_out ( DEBUG ) << yellow << "Filling hitMap for Layer '" << K << "'..." << normal << std::endl;
        m_vHitMapPerLayer.at(K - 1)->Fill(I, J);
        streamlog_out ( DEBUG ) << blue << "Filling hitMap for Layer '" << K << "'...OK" << normal << std::endl;


        cd.setCellID( caloHit ) ;
        // add layer to list of unique touched layers
        _firedLayersSet.insert(K);
        caloHit->setPosition(pos);
        col_event->addElement(caloHit);
        hitKeys.push_back(aHitKey);
      }
      else {
        streamlog_out( MESSAGE ) << " time peak = " << time_peak << " pointer --> : " << rawhit << std::endl;
      }
    }//loop over the hit
    hitKeys.clear();
  } catch (DataNotAvailableException &e) {
    streamlog_out(WARNING) << " collection not available" << std::endl;
  }
}

//=============================================================================
void TriventProc::defineColors()
{
  char cnormal[8] =  {0x1b, '[', '0', ';', '3', '9', 'm', 0};
  char cred[8]     = {0x1b, '[', '1', ';', '3', '1', 'm', 0};
  char cgreen[8]   = {0x1b, '[', '1', ';', '3', '2', 'm', 0};
  char cyellow[8]  = {0x1b, '[', '1', ';', '3', '3', 'm', 0};
  char cblue[8]    = {0x1b, '[', '1', ';', '3', '4', 'm', 0};
  char cmagenta[8] = {0x1b, '[', '1', ';', '3', '5', 'm', 0};
  char cwhite[8]   = {0x1b, '[', '1', ';', '3', '9', 'm', 0};

  normal   = cnormal;
  red      = cred;
  green    = cgreen;
  yellow   = cyellow;
  blue     = cblue;
  magenta  = cmagenta;
  white    = cwhite;
}
//=============================================================================
void TriventProc::init() {
  _trigCount = 0;
  evtnum = 0; // event number
  // ========================
  printParameters();

  // Define colors for
  defineColors();


  // Create writer for lcio output file
  _lcWriter = LCFactory::getInstance()->createLCWriter() ;
  _lcWriter->setCompressionLevel( 0 ) ;
  _lcWriter->open(_outFileName.c_str(), LCIO::WRITE_NEW) ;


  // Read and print geometry file
  XMLReader(_geomXML.c_str());
  printDifGeom();


  /**
   * Book root histograms
   */

  m_rootFile = new TFile(_rootFileName.c_str(), "RECREATE");

  if (NULL == m_rootFile)
    return;

  // Create a list of unique Layers from geometry file
  std::set<int> layerSet;
  for (std::map<int, LayerID>::iterator itt = _mapping.begin(); itt != _mapping.end(); itt++)
    layerSet.insert(itt->second.K);

  // Find last layer and resize the vector of hitMap
  const auto maxLayer = std::max_element(layerSet.begin(), layerSet.end());
  streamlog_out( MESSAGE ) << yellow << "Max LayerId in geometryFile: '" << *maxLayer << "'" << normal << std::endl;
  m_vHitMapPerLayer.resize(*maxLayer);

  // Check if first layer is numbered 0 or 1
  // Prevent accessing non defined element in vectors...
  const auto firstLayer = std::min_element(layerSet.begin(), layerSet.end());
  bool startAt0 = false;
  if ( 0 == *firstLayer )
    startAt0 = true;

  for (std::set<int>::const_iterator layerIter = layerSet.begin(); layerIter != layerSet.end(); ++layerIter)
  {
    unsigned int iLayer = *layerIter - 1;
    if (startAt0)
      iLayer = *layerIter;

    std::stringstream oss;
    streamlog_out( MESSAGE ) << yellow << "Booking hitMap for layer '" << iLayer << "'..." << normal << std::endl;
    oss << "hitMap_Layer" << iLayer;
    m_vHitMapPerLayer.at(iLayer) = (new TH2D(oss.str().c_str(), oss.str().c_str(), 96, 0, 96, 96, 0, 96));
    m_vHitMapPerLayer.at(iLayer)->GetXaxis()->SetTitle("I");
    m_vHitMapPerLayer.at(iLayer)->GetYaxis()->SetTitle("J (DIFSide)");

    streamlog_out( MESSAGE ) << blue << "Booking hitMap for layer '" << iLayer << "'...OK" << normal << std::endl;
  }

  for (const auto histo : m_vHitMapPerLayer)
    streamlog_out( MESSAGE ) << yellow << "Booked hitMap histo for layer '" << histo << "'" << normal << std::endl;
}

//=============================================================================
void TriventProc::processRunHeader( LCRunHeader * /*runHd*/ ) {
}

//=============================================================================
void TriventProc::processEvent( LCEvent * evtP ) {
  if (evtP != NULL) {
    try {

      _eventNr = evtP->getEventNumber();
      for (unsigned int i = 0; i < _hcalCollections.size(); i++) { //!loop over collection
        try {

          LCCollection * col = NULL;
          col = evtP ->getCollection(_hcalCollections[i].c_str());
          int numElements = col->getNumberOfElements();// hit number

          _trigCount++;
          if (0 == _trigCount % 10)
            streamlog_out( MESSAGE ) << yellow << "Trigger number == " << _trigCount << normal << std::endl;

          if (col == NULL )  {
            streamlog_out( WARNING ) << red << "TRIGGER SKIPED ... col == NULL" << normal << std::endl;
            break;
          }

          if (numElements > _elecNoiseCut)  {
            streamlog_out( MESSAGE ) << red << "TRIGGER number " << _trigCount << " SKIPPED ... hitNumber > _elecNoiseCut : " << numElements << " > " << _elecNoiseCut << normal << std::endl;
            break;
          }

          // set raw hits
          _trigger_raw_hit.clear();
          std::vector<int> vTrigger;
          for (int ihit(0); ihit < col->getNumberOfElements(); ++ihit) {// loop over the hits
            RawCalorimeterHit *raw_hit =
              dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit));

            if (NULL != raw_hit) {
              _trigger_raw_hit.push_back(raw_hit);
              //extract abolute bcid information:
              if (ihit == 0) {
                unsigned int difid = 0;
                difid = raw_hit->getCellID0() & 0xFF;
                if (difid == 0) return;
                std::stringstream pname("");
                pname << "DIF" << difid << "_Triggers";
                col->getParameters().getIntVals(pname.str(), vTrigger);
                if (vTrigger.size() != 0) {
                  _bcid1 = vTrigger[4] ;
                  _bcid2 = vTrigger[3] ;
                  unsigned long long Shift = 16777216ULL; //to shift the value from the 24 first bits
                  unsigned long long theBCID_ = _bcid1 * Shift + _bcid2;
                  streamlog_out( DEBUG1 ) << "trigger time : " << theBCID_ << std::endl;
                }
              }
            }
          }
          getMaxTime();
          std::vector<int> time_spectrum = getTimeSpectrum();

          //---------------------------------------------------------------
          //! Find the condidate event
          int prevTimePeak = 0; //  the previous bin center

          // Event is built at peakTime+-TimeWindow
          // Loop on time_spectrum vector without going out of range
          // This way we might miss the first/last event but sure not to access out of range value (like time_spectrum(-2)...)
          auto beginTime = std::next(time_spectrum.begin(), _timeWin);
          const auto & endTime = std::prev(time_spectrum.end(), _timeWin);
          auto & timeIter = beginTime;

          /**
           * Method used before to find peak ( does not take into account the time window)
           */
          // while (ibin < (_maxTime + 1)) {
          //   if (time_spectrum[ibin] >= _noiseCut &&
          //       time_spectrum[ibin] >= time_spectrum[ibin + 1] &&
          //       time_spectrum[ibin] >= time_spectrum[ibin - 1] &&
          //       time_spectrum[ibin] >= time_spectrum[ibin - 2] &&
          //       time_spectrum[ibin] >= time_spectrum[ibin + 2] ) {
          while (distance(timeIter, endTime) > 0 ) // Insure that timeIter < endTime
          {
            if ( *(timeIter) >= _noiseCut ) {
              // find the bin with max hits
              const auto & maxIter = std::max_element(std::prev(timeIter, _timeWin), std::next(timeIter, _timeWin));

              //find if the current bin has the max or equal hits
              if ( maxIter == timeIter || *(maxIter) == *(timeIter)) // if bin > other bins or bin is equal to biggest bin
              {
                // Found a peak at time *(timeIter)
                // std::cout << yellow
                //           << "\tmaxIter: " << std::distance(maxIter, timeIter)
                //           << "\t_timeWin: " << _timeWin
                //           << "\t_noiseCut: " << _noiseCut
                //           << "\ttime: " << *(timeIter)
                //           << "\tnextTime: " << *(std::next(timeIter))
                //           << "\tprevTime: " << *(std::prev(timeIter))
                //           << "\t2nextTime: " << *(std::next(timeIter, 2))
                //           << "\t2prevTime: " << *(std::prev(timeIter, 2))
                //           << normal << std::endl;

                LCEventImpl*  evt = new LCEventImpl() ;     // create the event

                //---------- set event paramters ------
                const std::string parname_trigger = "trigger";
                const std::string parname_energy  = "beamEnergy";
                const std::string parname_bcid1 = "bcid1";
                const std::string parname_bcid2 = "bcid2";
                evt->parameters().setValue(parname_trigger, evtP->getEventNumber());
                evt->parameters().setValue(parname_energy , _beamEnergy);
                evt->parameters().setValue(parname_bcid1 , _bcid1);
                evt->parameters().setValue(parname_bcid2 , _bcid2);
                evt->setRunNumber( evtP->getRunNumber()) ;
                m_runNumber = evtP->getRunNumber();
                //-------------------------------------

                LCCollectionVec* outcol = new LCCollectionVec(LCIO::CALORIMETERHIT);

                int timePeak = distance (time_spectrum.begin(), timeIter);
                streamlog_out( DEBUG0 ) << yellow << " EventBuilding with timePeak '" << timePeak << "' prevTimePeak: " << prevTimePeak << normal << std::endl;
                TriventProc::eventBuilder(outcol, timePeak, prevTimePeak);
                streamlog_out( DEBUG0 ) << yellow << " EventBuilding...OK" << normal << std::endl;

                streamlog_out( DEBUG0 ) << "_firedLayersSet.size() = " << _firedLayersSet.size() << "\t _LayerCut = " << _layerCut << std::endl;

                // Apply cut on min number of firedLayer +
                if ( (int)_firedLayersSet.size() < _layerCut )
                {
                  streamlog_out( DEBUG0 ) << green << " Event rejected, too few layer hit. nLayerHit: " << _firedLayersSet.size() << " _layerCut: " << _layerCut << normal << std::endl;
                  _rejectedNum++;
                  delete outcol;
                }
                //  Apply cut on time between two events
                else if (abs(timePeak - prevTimePeak) > _time2prevEventCut)
                {
                  streamlog_out( DEBUG0 ) << green << " Trivent find event at :==> " << red << timePeak
                                          << green << "\t :Nhit: ==> " << magenta
                                          << outcol->getNumberOfElements() << normal << std::endl;
                  evt->setEventNumber( evtnum++ ) ;
                  evt->addCollection(outcol, "SDHCAL_HIT");
                  _lcWriter->writeEvent( evt ) ;
                  _selectedNum++;
                }
                else {
                  streamlog_out( DEBUG0 ) << blue << " Event rejected, Events too close. eventTime: " << timePeak << " prevEventTime: " << prevTimePeak << normal << std::endl;
                  _rejectedNum++;
                  delete outcol;
                }

                delete evt; evt = NULL;

                prevTimePeak = timePeak;
                timeIter = std::next(timeIter, _timeWin);
              } else { // is not a peak, look in next frame
                ++timeIter;
              }
            } else { // Not enough hit in frame, look in next one
              ++timeIter;
            }
          }
        } catch (lcio::DataNotAvailableException zero) {}
      }
    } catch (lcio::DataNotAvailableException err) {}
  }
}

//=============================================================================
void TriventProc::end()
{
  streamlog_out( MESSAGE ) << "Trivent rejected " << _rejectedNum << " Condidate event" << std::endl;
  streamlog_out( MESSAGE ) << "Trivent Selected " << _selectedNum << " Condidate event" << std::endl;
  streamlog_out( MESSAGE ) << "Trivent end" << std::endl;
  _lcWriter->close();

  TCanvas *c1 = new TCanvas();
  c1->SetCanvasSize(2048, 1024);
  c1->Update();
  c1->cd();
  c1->Divide(2, 1);
  c1->cd(1);
  std::cout << "Drawing for layer 48" << std::endl;
  m_vHitMapPerLayer.at(47)->Draw("colz");
  c1->cd(2);
  std::cout << "Drawing for layer 50" << std::endl;
  m_vHitMapPerLayer.at(49)->Draw("colz");
  std::stringstream ss;
  ss << "/Volumes/PagosDisk/CALICE/data/SPS_06_2016/hitMap_Layer48-50_run" << m_runNumber << ".png";
  c1->SaveAs(ss.str().c_str());
  m_rootFile->cd();
  m_rootFile->Write();
  m_rootFile->Close();
}
//==============================================================



