/**
 * Yacine HADDAD
 * LLR Ecole polytechnique
 * avril 2012
 * Trivent v0.3
 */

/*
 *  TODO : LCIO parameters
 *          - startOfSpill Time Stamp
 *          - CerTag/TimeStamp
 *
 *      Add a cut if too much hit in One layer vs before/after (mean of +-3 layers)
 *          -> Remove Grounding Event?
 *          -> Remove event if RamFull on DIF
 *      ROOT
 *        timeSpectrum + limit
 *        rejected event
 *          too few hits
 *          too few layers hit
 *          Full asic
 *          too close events
 *
 */

// -- Trivent include
#include <TriventProc.hh>

// -- std includes
#include <algorithm> // std::max_element
#include <fstream>   // std::stringstream
#include <iterator>  // std::next

// -- Root headers
#include <TCanvas.h>

TriventProc a_TriventProc_instance;

//=========================================================
TriventProc::TriventProc()
    : Processor("TriventProc"),
      m_outputCollectionName("SDHCAL_HIT"),
      m_outFileName("TDHCAL.slcio"),
      m_noiseFileName("noise_run.slcio"),
      m_rootFileName("TDHCAL.root"),
      m_treeName("EventTree"),
      m_treeDescription("Event variables"),
      m_geomXMLFile("setup_geometry"),
      m_beamEnergy(0),
      m_useGainCorrection(false),
      m_elecNoiseCut(5000),
      m_noiseCut(10),
      m_layerCut(10),
      m_timeWin(2),
      m_time2prevEventCut(0),
      m_cellSizeI(10.408),
      m_cellSizeJ(10.408),
      m_layerThickness(26.131),
      m_hasCherenkov(true),
      m_cerenkovDifId(3),
      m_cerenkovTimeWindow(25),
      m_cerAsic(0),
      m_cerChan(0),
      m_cerThreshold(0),
      m_nCerenkov1(0),
      m_nCerenkov2(0),
      m_nCerenkov3(0),
      m_nCerenkovTrigger(0),
      m_hasTooManyCerenkov(false),
      m_timeCerenkov(0),
      m_totCerenkovHits(0),
      m_cerenkovEvts(0),
      m_trigNbr(0),
      m_trigCount(0),
      m_evtNum(0),
      m_selectedNum(0),
      m_rejectedNum(0),
      m_firedLayersSet{},
      m_bcid1(0),
      m_bcid2(0),
      m_runNumber(0),
      m_plotFolder("./"),
      m_evtTrigNbr(0),
      m_evtNbr(0),
      m_nHit(0),
      m_hitI{},
      m_hitJ{},
      m_hitK{},
      m_hitThreshold(0),
      m_nFiredLayers(0),
      m_isSelected(false),
      m_isNoise(false),
      m_isTooCloseInTime(false),
      m_hasNotEnoughLayers(false),
      m_hasFullAsic(false) {

  // collection
  std::vector<std::string> hcalCollections;
  hcalCollections.push_back(std::string("DHCALRawHits"));
  registerInputCollections(LCIO::RAWCALORIMETERHIT, "InputCollectionNames", "HCAL Collection Names", m_hcalCollections,
                           hcalCollections);

  registerOutputCollection(LCIO::CALORIMETERHIT, "OutputCollectionName", "HCAL Collection Names",
                           m_outputCollectionName, m_outputCollectionName);

  // Option of output file with clean events
  registerProcessorParameter("LCIOOutputFile", "LCIO file", m_outFileName, m_outFileName);
  // Energy
  registerProcessorParameter("beamEnergy", "The beam ", m_beamEnergy, m_beamEnergy);
  // Option of output file with noise
  registerProcessorParameter("NOISEOutputFile", "NOISE file", m_noiseFileName, m_noiseFileName);
  // layer cut
  registerProcessorParameter("LayerCut", "cut in number of layer 10 in default", m_layerCut, m_layerCut);

  // noise cut
  registerProcessorParameter("NoiseCut", "noise cut in time spectrum 10 in default", m_noiseCut, m_noiseCut);

  // time windows
  registerProcessorParameter("TimeWin", "time window = 2 in default", m_timeWin, m_timeWin);
  // maping on XML file
  registerProcessorParameter("SetupGeometry", "Dif geometry and position on the detector XML", m_geomXMLFile,
                             m_geomXMLFile);

  // electronic noise cut
  registerProcessorParameter("ElectronicNoiseCut", "number of hit max on time stamp", m_elecNoiseCut, m_elecNoiseCut);

  // electronic noise cut
  registerProcessorParameter("_time2prev_event_cut", "cut on time to previous event (x 200 ns)", m_time2prevEventCut,
                             m_time2prevEventCut);

  // Root Tree
  registerProcessorParameter("TreeName", "Logroot tree name", m_treeName, m_treeName);

  // Root Tree
  registerProcessorParameter("TreeDescription", "Logroot tree name", m_treeName, m_treeName);

  // histogram control tree
  registerProcessorParameter("ROOTOutputFile", "Logroot name", m_rootFileName, m_rootFileName);

  registerProcessorParameter("GainCorrectionMode", "m_useGainCorrection", m_useGainCorrection, m_useGainCorrection);

  registerProcessorParameter("HasCerenkovDIF", "If Cerenkov dif was connected during data taking", m_hasCherenkov,
                             m_hasCherenkov);

  registerProcessorParameter("CerenkovDifId", "Dif number for cerenkov data", m_cerenkovDifId, m_cerenkovDifId);

  registerProcessorParameter("CerenkovTimeWindow", "TimeWindow around timePeak in which to look for cerenkov data",
                             m_cerenkovTimeWindow, m_cerenkovTimeWindow);

  registerProcessorParameter("PlotFolder", "Folder Path to save Plot", m_plotFolder, m_plotFolder);
}

//=============================================================================
void TriventProc::XMLReader(const std::string &xmlfile) {
  TiXmlDocument xml(xmlfile.c_str());
  bool          load_key = xml.LoadFile();

  if (load_key) {
    streamlog_out(MESSAGE) << yellow << "Found Geometry File : " << xmlfile.c_str() << normal << std::endl;

    TiXmlHandle   xmlHandle(&xml);
    TiXmlElement *pElem = xmlHandle.FirstChildElement().Element();

    // should always have a valid root but handle gracefully if it does not
    if (!pElem) {
      streamlog_out(ERROR) << red << "error No root handle Found in xml file" << normal << std::endl;
    }

    // save this for later
    TiXmlHandle rootHandle(0);
    rootHandle = TiXmlHandle(pElem);

    // parameters block
    pElem = rootHandle.FirstChild("parameter").Element();
    assert(pElem);
    streamlog_out(DEBUG0) << "Reading data for key: " << pElem->Attribute("name") << std::endl;
    streamlog_out(DEBUG0) << "parameter : " << pElem->Attribute("name") << std::endl;

    std::string value = pElem->GetText();
    streamlog_out(DEBUG0) << "value : " << value << std::endl;

    std::string              line;
    std::vector<std::string> lines;
    std::istringstream       iss(value);

    // value no longer has the formatted eol, all replaced by a space character...
    while (std::getline(iss, line, ' ')) {
      // streamlog_out( MESSAGE ) << blue << line << normal << std::endl;
      lines.push_back(line);
    }
    streamlog_out(MESSAGE) << yellow << "Found " << lines.size()
                           << " difs in geometry file corresponding to : " << (unsigned int)lines.size() / 3
                           << " layers + " << lines.size() % 3 << " difs" << normal << std::endl;
    // for (std::vector<std::string>::const_iterator lineIter = lines.begin(); lineIter != lines.end(); ++lineIter) {
    for (const auto &lineIter : lines) {
      std::stringstream        ss(lineIter.c_str());
      std::vector<std::string> result;

      while (ss.good()) {
        std::string substr;
        getline(ss, substr, ',');
        // streamlog_out( MESSAGE ) << red << substr << normal << std::endl;
        result.push_back(substr);
      }

      LayerID mapp{};
      int     difId;
      while (ss.good()) {
        std::string substr;
        getline(ss, substr, ',');
        result.push_back(substr);
      }
      std::istringstream(result.at(0)) >> difId;
      std::istringstream(result.at(1)) >> mapp.K;
      std::istringstream(result.at(2)) >> mapp.DifX;
      std::istringstream(result.at(3)) >> mapp.DifY;
      std::istringstream(result.at(4)) >> mapp.IncX;
      std::istringstream(result.at(5)) >> mapp.IncY;
      m_mDifMapping[difId] = mapp;
      m_layerSet.insert(mapp.K);
    }
    // Cerenkov layer should not be counted it the total number of layer
    if (m_hasCherenkov) {
      m_cerenkovLayerId = m_mDifMapping.find(m_cerenkovDifId)->second.K;
      m_layerSet.erase(m_cerenkovLayerId);
    }
  } else {
    std::ostringstream oss;
    oss << "Failed to load geometry file '" << xmlfile.c_str() << "'";
    streamlog_out(WARNING) << red << oss.str() << normal << std::endl;
    throw(oss.str());
  }
}

//=============================================================================
void TriventProc::printDifGeom() {
  streamlog_out(DEBUG1) << "--- Dumping geomtry File: " << std::endl;
  // for (std::map<int, LayerID>::iterator itt = m_mDifMapping.begin(); itt != m_mDifMapping.end(); ++itt) {
  for (const auto &itt : m_mDifMapping) {
    streamlog_out(DEBUG1) << itt.first << "\t" << itt.second.K << "\t" << itt.second.DifX << "\t" << itt.second.DifY
                          << "\t" << itt.second.IncX << "\t" << itt.second.IncY << std::endl;
  }
}

// ============ decode the cell ids =============
// bit shift & 0xFF = Apply mask 1111 1111 to binary value
// eg: Dif 1 => cellID0 = 00983297 => DifID = 1 / AsicID = 1 / ChanID = 15
int TriventProc::getCellDif_id(const int &cellId) { return cellId & 0xFF; }

//=============================================================================
//  bit shift & 0xFF00 Apply mask 1111 1111 0000 0000 then cut last 8 bits
int TriventProc::getCellAsic_id(const int &cellId) { return (cellId & 0xFF00) >> 8; }

//=============================================================================
//  bit shift & 0x3F0000 Apply mask 1111 0000 0000 0000 0000 then cut last 16 bits
int TriventProc::getCellChan_id(const int &cellId) { return (cellId & 0x3F0000) >> 16; }

// ============ ============ ============ ============ ============ ============ ============
// ============ ============ ============ ============ ============ ============ ============
// Full example with Dif 1:
// cellId0 = 00983297 -> Binary =  1111 0000 0001 0000 0001
// binary & 0xFF = 0000 0001 => 2^0 = 1
// binary & 0xFF00 = 0000 0001 0000 0000 >>8 = 0000 0001 => 2^0 = 1
// binary & 0x3F0000 =  1111 0000 0000 0000 0000 >>16 = 1111 => (2^3)+(2^2)+(2^1)+(2^0) = 15
// ============ ============ ============ ============ ============ ============ ============
// ============ ============ ============ ============ ============ ============ ============

//=============================================================================

bool TriventProc::checkPadLimits(const std::vector<int> &padIndex, const std::vector<int> &padLimits) {
  assert(padLimits.size() == padIndex.size() * 2);
  for (int i = 0; i < static_cast<int>(padIndex.size()); ++i) {
    if (padIndex[i] < padLimits[i * 2] || padIndex[i] > padLimits[i * 2 + 1])
      return false;
  }
  return true;
}

//=============================================================================
std::vector<int> TriventProc::getPadIndex(const int &difId, const int &asicId, const int &chanId) {
  std::map<int, LayerID>::const_iterator findIter = m_mDifMapping.find(difId);

  if (findIter == m_mDifMapping.end()) {
    streamlog_out(ERROR) << " [getPadIndex] difId '" << difId << "' not found in geometry file" << std::endl;
    return {}; // empty
  }

  std::vector<int> index{1 + MapILargeHR2[chanId] + AsicShiftI[asicId],
                         32 - (MapJLargeHR2[chanId] + AsicShiftJ[asicId] ) + findIter->second.DifY,
                         static_cast<int>(findIter->second.K)};
  std::vector<int> padLims = {1, 96, 1, 96, 0, static_cast<int>(m_layerSet.size())};
  // Cerenkov layer is not in the layerSet as it's not a physical layer, needs to account for that when checking the pad
  // limits
  //
  if (difId == m_cerenkovDifId) {
    padLims.pop_back();
    padLims.push_back(m_cerenkovLayerId);
  }
  assert(checkPadLimits(index, padLims));

  // if (difId == m_cerenkovDifId) {
  // streamlog_out(DEBUG0) << " difId== " << difId<< " asicId ==" << asicId << " chanId ==" << chanId
  // << " I == " << index[0] << " J == " << index[1] << " K == " << index[2] << std::endl;
  // }
  return index;
}

//=============================================================================
int TriventProc::getMaxTime() {
  assert(!m_triggerRawHitMap.empty());
  return m_triggerRawHitMap.rbegin()->first;
}

//=============================================================================
std::vector<int> TriventProc::getTimeSpectrum(const int &maxTime) //__attribute__((optimize(0)))
{
  std::vector<int> timeSpectrumVec(maxTime + 1, 0);
  for (const auto &mapIt : m_triggerRawHitMap) {
    int time = mapIt.first;
    for (const auto &rawIt : mapIt.second)
      assert(time == rawIt->getTimeStamp());
    assert(time <= maxTime);
    assert(time >= 0);
    if (time >= 0) {
      timeSpectrumVec.at(time) += mapIt.second.size();
    }
  }
  return timeSpectrumVec;
}

//=============================================================================
int TriventProc::IJKToKey(const std::vector<int> &padIndex) {
  return 100 * 100 * padIndex[2] + 100 * padIndex[1] + padIndex[0];
}

//=============================================================================
int TriventProc::getAsicKey(const std::vector<int> &padIndex) {
  // Not necessary to check for boundary here as already tested in getPadIndex
  std::vector<int> padLims = {1, 96, 1, 96, 0, static_cast<int>(m_layerSet.size())};
  if (padIndex[2] == m_cerenkovLayerId) {
    padLims.pop_back();
    padLims.push_back(m_cerenkovLayerId);
  }
  assert(checkPadLimits(padIndex, padLims));

  const int jnum = (padIndex[1] - 1) / 8;
  const int inum = (padIndex[0] - 1) / 8;
  const int num  = jnum * 12 + inum;
  return padIndex[2] * 1000 + num;
}

//=============================================================================
void TriventProc::resetTriggerParameters() {
  m_nCerenkovTrigger   = 0;
  m_hasTooManyCerenkov = false;
  m_triggerRawHitMap.clear();
  m_cerenkov_raw_hit.clear();
}

//=============================================================================
void TriventProc::resetEventParameters() {
  // reset Flags
  m_isSelected         = false;
  m_isNoise            = false;
  m_hasNotEnoughLayers = false;
  m_hasFullAsic        = false;
  m_isTooCloseInTime   = false;
  m_nHit               = 0;
  m_hitI.clear();
  m_hitJ.clear();
  m_hitK.clear();
  m_hitThreshold.clear();
  m_firedLayersSet.clear();

  m_nFiredLayers = 0;
  m_nCerenkov1   = 0;
  m_nCerenkov2   = 0;
  m_nCerenkov3   = 0;
  m_timeCerenkov = -2 * m_cerenkovTimeWindow;
}

//=============================================================================
void TriventProc::eventBuilder(std::unique_ptr<IMPL::LCCollectionVec> &evtCol, const int &timePeak,
                               const unsigned int &lowTimeBoundary, const unsigned int &highTimeBoundary) {

  resetEventParameters();

  evtCol->setFlag(evtCol->getFlag() | (1 << LCIO::RCHBIT_LONG));
  evtCol->setFlag(evtCol->getFlag() | (1 << LCIO::RCHBIT_TIME));

  CellIDEncoder<CalorimeterHitImpl> cellIdEncoder("M:3,S-1:3,I:9,J:9,K-1:6", evtCol.get());

  std::map<int, int> asicMap;
  std::map<int, int> hitKeys;
  // TODO: Make it nicer. Here it just ensure cerenkov hit within m_cerenkovTimeWindow are not discarded from the
  // Calorimeter hit
  // unsigned int timeWindow = m_timeWin;
  // if (getCellDif_id(rawhit->getCellID0()) == m_cerenkovDifId)
  //   timeWindow = m_cerenkovTimeWindow;

  for (unsigned int hitTime = lowTimeBoundary; hitTime <= highTimeBoundary; ++hitTime) {

    // No hit recorded at hitTime
    if (m_triggerRawHitMap.find(hitTime) == m_triggerRawHitMap.end())
      continue;

    for (const auto &rawHit : m_triggerRawHitMap.at(hitTime)) {
      assert(rawHit);

      const int difId  = getCellDif_id(rawHit->getCellID0());
      int       asicId = getCellAsic_id(
          rawHit->getCellID0()); // Can't be const to correct for cerenkovAsicId bug (in data from 2014>2016)
      const int chanId = getCellChan_id(rawHit->getCellID0());
      const int thresh = rawHit->getAmplitude();
      if ((asicId < 1) || (asicId > 48)) {
        if (m_hasCherenkov && difId == m_cerenkovDifId) {
          streamlog_out(WARNING) << yellow << "Cerenkov has a weird asicId : '" << asicId
                                 << " correcting it... to asicId = 1" << normal << std::endl;
          asicId = 1;
        } else {
          streamlog_out(ERROR) << red << "[eventBuilder] - Found a hit with weird AsicId, Dif/Asic/Chan/Thr... "
                               << difId << "/" << asicId << "/" << chanId << "/" << thresh << normal << std::endl;
          continue;
        }
      }
      if (chanId > 63) {
        streamlog_out(ERROR) << yellow << "[eventBuilder] - Found a hit with weird ChannelId, Dif/Asic/Chan/Thr... "
                             << difId << "/" << asicId << "/" << chanId << "/" << thresh << normal << std::endl;
        continue;
      }
      const std::vector<int> padIndex = getPadIndex(difId, asicId, chanId);
      if (padIndex.empty()) {
        streamlog_out(ERROR) << red << "[eventBuilder] - Dif '" << difId
                             << "' not found in geometry file...skipping hit" << normal << std::endl;
        continue;
      }

      // find and remove square events
      const int asicKey = getAsicKey(padIndex);
      assert(asicKey > 0);
      if (asicMap[asicKey]) {
        ++asicMap[asicKey];
      } else {
        asicMap[asicKey] = 1;
      }

      if (asicMap[asicKey] == 64 && difId != m_cerenkovDifId) {
        streamlog_out(WARNING) << yellow << "[eventBuilder] - Rejecting event with full asic. Dif '" << difId
                               << "' asic '" << asicId << "' at time '" << timePeak << "'" << normal << std::endl;

        m_firedLayersSet.clear();
        hitKeys.clear();
        asicMap.clear();
        m_isSelected  = false;
        m_hasFullAsic = true;
      }

      // Creating Calorimeter Hit
      float pos[3];
      pos[0] = padIndex[0] * m_cellSizeI;
      pos[1] = padIndex[1] * m_cellSizeJ;
      pos[2] = padIndex[2] * m_layerThickness;

      CalorimeterHitImpl *caloHit = new CalorimeterHitImpl();
      caloHit->setTime(static_cast<float>(rawHit->getTimeStamp()));

      const float hitShiftedAmplitude = static_cast<float>(rawHit->getAmplitude() & 3);
      if (hitShiftedAmplitude > 2.5)
        caloHit->setEnergy(hitShiftedAmplitude); // 3rd threshold
      else if (hitShiftedAmplitude > 1.5)
        caloHit->setEnergy(hitShiftedAmplitude - 1); // 2nd threshold ?
      else
        caloHit->setEnergy(hitShiftedAmplitude + 1); // 1st threshold ?

      // Create hit Key
      const int aHitKey = IJKToKey(padIndex);

      // Avoid two hit in the same cell
      std::map<int, int>::const_iterator findIter = hitKeys.find(aHitKey);

      if (findIter != hitKeys.end()) {
        if (difId != m_cerenkovDifId) {
          delete caloHit;
          continue;
        } else {
          streamlog_out(ERROR) << yellow << " So much Hit in my cherenkov " << normal << std::endl;
          abort();
        }
      }

      const int I = padIndex[0];
      const int J = padIndex[1];
      const int K = padIndex[2];

      // set the cell id
      cellIdEncoder["I"]   = I;
      cellIdEncoder["J"]   = J;
      cellIdEncoder["K-1"] = K - 1;
      cellIdEncoder["M"]   = 0;
      cellIdEncoder["S-1"] = 3;

      if (difId == m_cerenkovDifId) {
        m_totCerenkovHits++;
        streamlog_out(DEBUG0) << " m_totCerenkovHits == " << m_totCerenkovHits << " I == " << I << " J == " << J
                              << " K == " << K << " dif == " << difId << " asic == " << asicId << " chan == " << chanId
                              // << " pos[0] == " << pos[0]
                              // << " pos[1] == " << pos[1]
                              // << " pos[2] == " << pos[2]
                              << " ahitKey == " << aHitKey << " time == " << rawHit->getTimeStamp()
                              << " time - timePeak == " << rawHit->getTimeStamp() - timePeak << std::endl;
        m_vHitMapPerLayer.back()->Fill(I, J);
      } else {
        // Fill hitsMap for each Layer, starting at K-1 = 0
        // streamlog_out(DEBUG) << yellow << "Filling hitMap for Layer '" << K << "'..." << normal << std::endl;
        assert(m_vHitMapPerLayer.at(K - 1));
        m_vHitMapPerLayer.at(K - 1)->Fill(I, J);
        // streamlog_out(DEBUG) << blue << "Filling hitMap for Layer '" << K << "'...OK" << normal << std::endl;
      }

      cellIdEncoder.setCellID(caloHit);
      // add layer to list of unique touched layers
      m_firedLayersSet.insert(K);
      caloHit->setPosition(pos);
      evtCol->addElement(caloHit);
      hitKeys.insert(std::pair<int, int>(aHitKey, hitTime));
      m_hitI.push_back(I);
      m_hitJ.push_back(J);
      m_hitK.push_back(K);
      // m_hitBCID.push_back((*rawhit)->getTimeStamp());
      m_hitThreshold.push_back(thresh);
    }
  }
  m_nFiredLayers = m_firedLayersSet.size();
}

//=============================================================================
void TriventProc::defineColors() {
  char cnormal[8]  = {0x1b, '[', '0', ';', '3', '9', 'm', 0};
  char cred[8]     = {0x1b, '[', '1', ';', '3', '1', 'm', 0};
  char cgreen[8]   = {0x1b, '[', '1', ';', '3', '2', 'm', 0};
  char cyellow[8]  = {0x1b, '[', '1', ';', '3', '3', 'm', 0};
  char cblue[8]    = {0x1b, '[', '1', ';', '3', '4', 'm', 0};
  char cmagenta[8] = {0x1b, '[', '1', ';', '3', '5', 'm', 0};
  char cwhite[8]   = {0x1b, '[', '1', ';', '3', '9', 'm', 0};

  normal  = cnormal;
  red     = cred;
  green   = cgreen;
  yellow  = cyellow;
  blue    = cblue;
  magenta = cmagenta;
  white   = cwhite;
}

//=============================================================================
TTree *TriventProc::getOrCreateTree(const std::string &treeName, const std::string &treeDescription) {
  TTree *tree = static_cast<TTree *>(m_rootFile->Get(treeName.c_str()));

  if (!tree) {
    streamlog_out(DEBUG0) << "Creating tree '" << treeName << "'" << std::endl;
    tree = new TTree(treeName.c_str(), treeDescription.c_str());
  }
  assert(tree);
  return tree;
}

//=============================================================================
void TriventProc::initRootTree() {
  m_rootFile = new TFile(m_rootFileName.c_str(), "RECREATE");
  assert(m_rootFile);

  // Create Trigger tree & Branches
  // m_triggerTree = getOrCreateTree("TriggerTree", "Trigger variables");
  // m_triggerTree->Branch("TriggerNumber",                    &m_trigNbr);
  // m_triggerTree->Branch("NumberOfEvents",                   &m_nEvt);
  // m_triggerTree->Branch("TimeSpectrum", "std::vector<int>", &m_vTimeSpectrum);

  // Create Event tree & Branches
  m_eventTree = getOrCreateTree(m_treeName, m_treeDescription);
  assert(m_eventTree);
  m_eventTree->Branch("TriggerNumber", &m_evtTrigNbr);
  m_eventTree->Branch("EventNumber", &m_evtNbr);
  // m_eventTree->Branch("Hitbcid",                 &m_hitBCID);
  m_eventTree->Branch("NumberOfHits", &m_nHit);
  m_eventTree->Branch("HitI", &m_hitI);
  m_eventTree->Branch("HitJ", &m_hitJ);
  m_eventTree->Branch("HitK", &m_hitK);
  m_eventTree->Branch("HitThreshold", &m_hitThreshold);
  m_eventTree->Branch("NumberOfFiredLayers", &m_nFiredLayers);
  m_eventTree->Branch("CerAsic", &m_cerAsic);
  m_eventTree->Branch("CerChan", &m_cerChan);
  m_eventTree->Branch("CerThreshold", &m_cerThreshold);
  m_eventTree->Branch("NumberOfCerenkov1Hits", &m_nCerenkov1);
  m_eventTree->Branch("NumberOfCerenkov2Hits", &m_nCerenkov2);
  m_eventTree->Branch("NumberOfCerenkov3Hits", &m_nCerenkov3);
  m_eventTree->Branch("TriggerHasTooManyCerenkov", &m_hasTooManyCerenkov);
  m_eventTree->Branch("CerenkovTime", &m_timeCerenkov);
  m_eventTree->Branch("EventIsSelected", &m_isSelected);
  m_eventTree->Branch("EventIsNoise", &m_isNoise);
  m_eventTree->Branch("EventIsToCloseFromLast", &m_isTooCloseInTime);
  m_eventTree->Branch("EventHasNotEnoughLayers", &m_hasNotEnoughLayers);
  m_eventTree->Branch("EventIsHasFullAsic", &m_hasFullAsic);

  TDirectory *rootDir   = gDirectory;
  TDirectory *hitMapDir = rootDir->mkdir("HitMapPerLayer");
  hitMapDir->cd();

  // Check if first layer is numbered 0 or 1
  // Prevent accessing non defined element in vectors...
  const auto firstLayer = std::min_element(m_layerSet.begin(), m_layerSet.end());
  bool       startAt0   = false;
  if (0 == *firstLayer) {
    startAt0 = true;
  }

  for (const auto &layerIter : m_layerSet) {
    int iLayer = layerIter - 1;
    if (startAt0) {
      iLayer = layerIter;
    }

    std::stringstream oss;
    oss << "hitMap_Layer" << iLayer;
    streamlog_out(DEBUG0) << "Booking hitMap for layer '" << iLayer << "'..." << std::endl;
    m_vHitMapPerLayer.push_back(makeTH2(oss.str(), "I", "J (DIFSide)"));
    assert(m_vHitMapPerLayer.back());
    streamlog_out(DEBUG0) << "Booking hitMap for layer '" << iLayer << "'...OK" << std::endl;
  }

  std::stringstream oss;
  oss << "hitMap_Cerenkov";
  streamlog_out(DEBUG0) << "Booking hitMap for Cerenkov..." << std::endl;
  m_vHitMapPerLayer.push_back(makeTH2(oss.str(), "I", "J (DIFSide)"));
  assert(m_vHitMapPerLayer.back());
  streamlog_out(DEBUG0) << "Booking hitMap for Cerenkov...OK" << std::endl;

  int iLayer = 0;
  for (auto const &histo : m_vHitMapPerLayer) {
    assert(histo);
    streamlog_out(DEBUG0) << yellow << "Booked hitMap histo for layer '" << iLayer << "' at --> '" << histo << "'"
                          << normal << std::endl;
    ++iLayer;
  }
}

//=============================================================================
void TriventProc::init() {
  m_trigNbr   = 0;
  m_trigCount = 0;
  m_evtNum    = 0; // event number
  // ========================
  printParameters();
  defineColors();

  // Create writer for lcio output file
  m_lcWriter = std::unique_ptr<LCWriter>(LCFactory::getInstance()->createLCWriter());
  m_lcWriter->setCompressionLevel(0);
  m_lcWriter->open(m_outFileName.c_str(), LCIO::WRITE_NEW);

  // Read and print geometry file
  try {
    XMLReader(m_geomXMLFile);
  } catch (std::string &e) {
    std::cout << "\n------------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "\t ****** Caught Exception when parsing geometry file: " << e << std::endl;
    std::cout << "------------------------------------------------------------------------------------------\n"
              << std::endl;
    throw;
  } catch (...) {
    std::cout << "\n------------------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "\t ****** Uncaught Exception when parsing geometry file! " << std::endl;
    std::cout << "------------------------------------------------------------------------------------------\n"
              << std::endl;
    throw;
  }
  printDifGeom();

  initRootTree();
}

//=============================================================================
TH2 *TriventProc::makeTH2(const std::string &title, const std::string &xTitle, const std::string &yTitle) {
  TH2 *hMap = new TH2D(title.c_str(), title.c_str(), 96, 1, 97, 96, 1, 97);
  hMap->GetXaxis()->SetTitle(xTitle.c_str());
  hMap->GetYaxis()->SetTitle(yTitle.c_str());
  return hMap;
}

//=============================================================================
void TriventProc::findCerenkovHits(const int &timePeak) {
  for (const auto &cerHit : m_cerenkov_raw_hit) {
    assert(cerHit);
    const int bifTime = static_cast<int>(cerHit->getTimeStamp());

    if (std::fabs(bifTime - timePeak) <= m_cerenkovTimeWindow) {
      m_timeCerenkov         = bifTime - timePeak;
      const int difId        = getCellDif_id(cerHit->getCellID0());
      const int asicId       = getCellAsic_id(cerHit->getCellID0());
      const int chanId       = getCellChan_id(cerHit->getCellID0());
      const int hitThreshold = cerHit->getAmplitude();

      if (difId != m_cerenkovDifId) {
        streamlog_out(WARNING) << yellow << "[findCerenkov] - Found Cerenkov hit in wrong dif '" << difId
                               << "'... should be in dif '" << m_cerenkovDifId << "'...skipping" << normal << std::endl;
        continue;
      }

      streamlog_out(DEBUG) << "[findCerenkov] - Found Cerenkov hit at time '" << m_timeCerenkov << "'\t Asic '"
                           << asicId << "'\t Chan '" << chanId << "'\t Threshold " << hitThreshold << std::endl;

      m_cerAsic      = asicId;
      m_cerChan      = chanId;
      m_cerThreshold = hitThreshold;

      switch (hitThreshold) {
      case 1:
        m_nCerenkov1 += 1;
        break;

      case 2:
        m_nCerenkov2 += 1;
        break;

      case 3:
        m_nCerenkov3 += 1;
        break;

      default:
        streamlog_out(ERROR) << red << "[findCerenkov] - Found Cerenkov hit with weird threshold : '" << hitThreshold
                             << "'..." << normal << std::endl;
        break;
      }
    }
  }
  m_nCerenkovTrigger += (m_nCerenkov1 + m_nCerenkov2 + m_nCerenkov3);
  if (m_nCerenkovTrigger > m_cerenkov_raw_hit.size()) {
    streamlog_out(WARNING) << yellow << "[findCerenkov] - Cerenkov hit associated with multiple event in Trigger : "
                                        "Associated hit/total bif_hit in trigger : '"
                           << m_nCerenkovTrigger << "'/" << m_cerenkov_raw_hit.size() << normal << std::endl;
    m_hasTooManyCerenkov = true;
  }
}

//=============================================================================
void TriventProc::fillRawHitTrigger(const LCCollection &inputLCCol) {
  std::vector<int> vTrigger;

  for (int ihit(0); ihit < inputLCCol.getNumberOfElements(); ++ihit) // loop over the hits
  {
    RawCalorimeterHit *raw_hit = dynamic_cast<RawCalorimeterHit *>(inputLCCol.getElementAt(ihit));
    if (raw_hit) {
      // extract abolute bcid information:
      const int difId = raw_hit->getCellID0() & 0xFF;
      assert(difId > 0);
      if (ihit == 0) {
        std::stringstream pname("");
        pname << "DIF" << difId << "_Triggers";
        inputLCCol.getParameters().getIntVals(pname.str(), vTrigger);
        if (vTrigger.size() != 0) {
          m_bcid1                     = vTrigger[4];
          m_bcid2                     = vTrigger[3];
          unsigned long long Shift    = 16777216ULL; // to shift the value from the 24 first bits
          unsigned long long theBCID_ = m_bcid1 * Shift + m_bcid2;
          streamlog_out(DEBUG0) << "trigger time : " << theBCID_ << std::endl;
        }
      }

      if (difId == m_cerenkovDifId) {
        m_cerenkov_raw_hit.push_back(raw_hit);
      }
      m_triggerRawHitMap[raw_hit->getTimeStamp()].push_back(raw_hit);
    }
  }

  streamlog_out(MESSAGE) << blue << " Trigger '" << m_trigCount << "' Found " << m_cerenkov_raw_hit.size()
                         << " raw hits in BIF!" << normal << std::endl;
  streamlog_out(DEBUG1) << "at time : " << normal << std::endl;
  for (const auto &cerHit : m_cerenkov_raw_hit) {
    assert(cerHit);
    streamlog_out(DEBUG1) << blue << " \t '" << static_cast<int>(cerHit->getTimeStamp()) << " Cerenkov --> '" << cerHit
                          << normal << std::endl;
  }
}

//=============================================================================
std::vector<std::vector<int>::iterator>
TriventProc::getCandidateTimeBoundaries(std::vector<int>::iterator &beginTime, std::vector<int>::iterator &endTime,
                                        std::vector<int>::iterator &candidateTime) {
  assert(beginTime < endTime);
  assert(candidateTime >= beginTime);
  assert(candidateTime < endTime); // Can't be equal otherwise we'll access out of bound value on next ++timeIter call
  std::vector<int>::iterator lowerBound = endTime;
  std::vector<int>::iterator upperBound = beginTime;

  // Ensure there is sufficient time between two candidate + we are not looking before begining of timeSpectrumVec
  auto timeDistance = std::distance(beginTime, candidateTime);
  if (timeDistance > m_timeWin)
    lowerBound = std::prev(candidateTime, m_timeWin);
  // Don't throw a potential candidate found in the first few frames of the trigger
  else {
    streamlog_out(DEBUG0) << green << "[processEvent] - small lowerBound! m_timeWin : " << m_timeWin
                          << " distance(beginTime, candidateTime) = " << timeDistance << normal << std::endl;
    lowerBound = std::prev(candidateTime, timeDistance);
  }

  // Check we are sufficiently far from end of timeSpectrumVec
  timeDistance = std::distance(candidateTime, endTime);
  if (timeDistance > m_timeWin)
    upperBound = std::next(candidateTime, m_timeWin);
  else { // Don't throw a potential candidate found in the last few frames of the trigger
    upperBound = std::next(candidateTime, timeDistance); // distance > 0 already met in while loop
    streamlog_out(DEBUG0) << green << "[processEvent] - small upperBound! m_timeWin : " << m_timeWin
                          << " distance(candidateTime, endTime) = " << timeDistance << normal << std::endl;
  }
  streamlog_out(DEBUG0) << "low: " << distance(beginTime, lowerBound) << " up " << distance(beginTime, upperBound)
                        << std::endl;
  assert(lowerBound < upperBound);
  assert(std::distance(upperBound, lowerBound) <= 2 * m_timeWin);
  assert(std::distance(lowerBound, upperBound) > m_timeWin);
  return {lowerBound, upperBound};
}

//=============================================================================
void TriventProc::processEvent(LCEvent *evtP) {
  assert(evtP != NULL);

  m_trigNbr = evtP->getEventNumber();
  if (m_trigNbr > 1E6) {
    streamlog_out(ERROR) << red << "Too much Triggers : " << m_trigNbr << normal << std::endl;
    return;
  }

  for (unsigned int i = 0; i < m_hcalCollections.size(); i++) //! loop over collection
  {
    LCCollection *inputLCCol;
    try {
      inputLCCol = evtP->getCollection(m_hcalCollections.at(i).c_str());
    } catch (lcio::DataNotAvailableException &zero) {
      streamlog_out(ERROR) << red << "No data found in collection " << i << normal << std::endl;
    }

    if (!inputLCCol) {
      streamlog_out(WARNING) << red << "TRIGGER SKIPED ... col is nullptr" << normal << std::endl;
      continue;
    }

    ++m_trigCount;
    if (0 == m_trigCount % 100) {
      streamlog_out(MESSAGE) << yellow << "Trigger number == " << m_trigCount << normal << std::endl;
    }

    const int numElements = inputLCCol->getNumberOfElements(); // hit number in trigger
    if (numElements > m_elecNoiseCut) {
      streamlog_out(MESSAGE) << yellow << "TRIGGER number " << m_trigCount
                             << " SKIPPED ... hitNumber > m_elecNoiseCut : " << numElements << " > " << m_elecNoiseCut
                             << normal << std::endl;
      continue;
    }

    // set raw hits
    resetTriggerParameters();
    fillRawHitTrigger(*inputLCCol);
    std::vector<int> timeSpectrumVec = getTimeSpectrum(getMaxTime());

    //---------------------------------------------------------------
    //! Find the candidate event
    // Event is built at peakTime+-TimeWindow
    // Loop on timeSpectrumVec vector without going out of range
    std::vector<int>::iterator beginTimeIter = timeSpectrumVec.begin();
    std::vector<int>::iterator endTimeIter   = timeSpectrumVec.end();
    std::vector<int>::iterator timeIter      = beginTimeIter;
    std::vector<int>::iterator prevMaxIter   = beginTimeIter;

    while (std::distance(timeIter, endTimeIter) > 0) { // Ensure that timeIter < endTime
      if (*(timeIter) < m_noiseCut) {                  // Not enough hit in frame, look in next one
        ++timeIter;
        continue;
      }

      // find timeBoundaries to build the event
      const auto  boundaries = getCandidateTimeBoundaries(prevMaxIter, endTimeIter, timeIter);
      const auto &maxIter    = std::max_element(boundaries[0], boundaries[1]); // max in [lower,upper))

      if (maxIter > timeIter) {
        // timeIter is not a real peak yet, look in next frame
        ++timeIter;
        continue;
      }

      streamlog_out(DEBUG0) << yellow << "[processEvent] - upperBound '" << std::distance(beginTimeIter, boundaries[1])
                            << "' lowerBound '" << std::distance(beginTimeIter, boundaries[0]) << " max at time '"
                            << std::distance(beginTimeIter, maxIter) << "' : " << *maxIter << normal << std::endl;

      // Check we didn't already process the peak
      // if (std::distance(maxIter, prevMaxIter) < m_timeWin && std::distance(beginTimeIter, timeIter) > 0) {
      if (maxIter < prevMaxIter && std::distance(beginTimeIter, timeIter) > 0) {
        if (std::distance(maxIter, prevMaxIter) < m_timeWin) {
          streamlog_out(DEBUG0) << yellow << "[processEvent] - Found duplicate peak, at time '"
                                << std::distance(beginTimeIter, maxIter) << "' previous peak : '"
                                << std::distance(beginTimeIter, prevMaxIter) << "'..." << normal << std::endl;
          ++timeIter;
          continue;
        }
      }

      const int          timePeak     = distance(beginTimeIter, maxIter);
      const int          prevTimePeak = distance(beginTimeIter, prevMaxIter);
      const unsigned int lowBound     = distance(beginTimeIter, boundaries[0]);
      const unsigned int highBound    = distance(beginTimeIter, boundaries[1]);
      prevMaxIter                     = maxIter;

      streamlog_out(DEBUG0) << blue << "[processEvent] - Trig '" << m_trigCount << "' : Found Peak, at time '"
                            << timePeak << "' - hits : " << *maxIter << " prevTimePeak: " << prevTimePeak << normal
                            << std::endl;

      //---------- set event paramters ------
      std::unique_ptr<LCEventImpl> lcEvt           = make_unique<LCEventImpl>(); // create the event
      const std::string            parname_trigger = "trigger";
      const std::string            parname_energy  = "beamEnergy";
      const std::string            parname_bcid1   = "bcid1";
      const std::string            parname_bcid2   = "bcid2";
      lcEvt->parameters().setValue(parname_trigger, evtP->getEventNumber());
      lcEvt->parameters().setValue(parname_energy, m_beamEnergy);
      lcEvt->parameters().setValue(parname_bcid1, m_bcid1);
      lcEvt->parameters().setValue(parname_bcid2, m_bcid2);
      lcEvt->setRunNumber(evtP->getRunNumber());
      m_runNumber = evtP->getRunNumber();
      //-------------------------------------

      std::unique_ptr<LCCollectionVec> outCol = make_unique<LCCollectionVec>(LCIO::CALORIMETERHIT);

      // Event Building
      TriventProc::eventBuilder(outCol, timePeak, lowBound, highBound);
      streamlog_out(DEBUG0) << blue << " EventBuilding...OK" << normal << std::endl;

      m_evtTrigNbr = m_trigNbr;
      m_evtNbr     = m_evtNum; // dont increment here: rejected event will have same number as last accepted !
      m_nHit       = outCol->getNumberOfElements();

      if (m_hasCherenkov) {
        // streamlog_out( MESSAGE ) << green << "hit in bif: " << m_cerenkov_raw_hit.size() << normal <<
        // std::endl;

        streamlog_out(DEBUG0) << " Find Cer " << m_cerenkov_raw_hit.size() << std::endl;
        findCerenkovHits(timePeak);
        streamlog_out(DEBUG0) << " Find Cer OK" << normal << std::endl;

        if (m_timeCerenkov != -2 * m_cerenkovTimeWindow) // InitialValue
        {
          streamlog_out(DEBUG) << "[processEvent] - " << green << "Trig# " << m_evtTrigNbr << " TrigCount "
                               << m_trigCount << " Evt# " << m_evtNum << "\tFound Cerenkov!"
                               << "\t cer1 = " << m_nCerenkov1 << "\t cer2 = " << m_nCerenkov2
                               << "\t cer3 = " << m_nCerenkov3 << normal << std::endl;
        }
      }

      // Apply cut on min number of firedLayer +
      if (static_cast<int>(m_nFiredLayers) < m_layerCut) {
        streamlog_out(DEBUG0) << green << " Event rejected, too few layer hit. nLayerHit: " << m_nFiredLayers
                              << " m_layerCut: " << m_layerCut << normal << std::endl;
        ++m_rejectedNum;
        m_isSelected         = false;
        m_hasNotEnoughLayers = true;
      }

      //  Apply cut on time between two events
      assert(timePeak - prevTimePeak > m_time2prevEventCut);

      // else if (abs(timePeak - prevTimePeak) > m_time2prevEventCut) {
      streamlog_out(DEBUG0) << green << " Trivent find event at :==> " << red << timePeak << green << "\t :Nhit: ==> "
                            << magenta << outCol->getNumberOfElements() << normal << std::endl;
      lcEvt->setEventNumber(++m_evtNum);

      lcEvt->addCollection(outCol.release(), m_outputCollectionName);
      m_lcWriter->writeEvent(lcEvt.get());
      assert(lcEvt);

      ++m_selectedNum;
      m_isSelected = true;
      // } else {
      //   streamlog_out(MESSAGE) << blue << " Event rejected, Events too close. eventTime: " << timePeak
      //                          << " prevEventTime: " << prevTimePeak << normal << std::endl;
      //   ++m_rejectedNum;
      //   m_isSelected       = false;
      //   m_isTooCloseInTime = true;
      // }

      if (m_nCerenkov1 > 0 || m_nCerenkov2 > 0 || m_nCerenkov3 > 0) {
        ++m_cerenkovEvts;
      }

      m_eventTree->Fill();
      timeIter = std::next(timeIter, m_timeWin + 1);
    }
    // m_vTimeSpectrum = timeSpectrumVec;
    // m_triggerTree->Fill();
  }
}

//=============================================================================
void TriventProc::end() {

  streamlog_out(MESSAGE) << "Trivent rejected " << m_rejectedNum << " Condidate event" << std::endl;
  streamlog_out(MESSAGE) << "Trivent Selected " << m_selectedNum << " Condidate event" << std::endl;
  streamlog_out(MESSAGE) << "Cerenkov Event Selected " << m_cerenkovEvts << std::endl;

  std::string cerCut = "CerenkovTime>-" + std::to_string(m_cerenkovTimeWindow);
  m_eventTree->Draw("CerenkovTime>>hcer", cerCut.c_str());
  // std::unique_ptr<TH1> hcer(dynamic_cast<TH1 *>(gDirectory->Get("hcer")));
  TH1 *hcer = dynamic_cast<TH1 *>(gDirectory->Get("hcer"));
  streamlog_out(MESSAGE) << "Cerenkov Probable time shift at " << hcer->GetXaxis()->GetBinCenter(hcer->GetMaximumBin())
                         << " with " << hcer->GetBinContent(hcer->GetMaximumBin()) << "/" << m_cerenkovEvts
                         << " event tagged" << std::endl;
  int cerTagInTime = 0;
  for (int i = -m_timeWin; i < m_timeWin + 1; ++i) {
    streamlog_out(MESSAGE) << "Bin '" << hcer->GetXaxis()->GetBinCenter(hcer->GetMaximumBin() + i)
                           << "' :  " << hcer->GetBinContent(hcer->GetMaximumBin() + i) << std::endl;
    cerTagInTime += hcer->GetBinContent(hcer->GetMaximumBin() + i);
  }
  streamlog_out(ERROR) << "TotCerenkov in MPV+-timeWin : " << cerTagInTime << "/" << m_cerenkovEvts << " ("
                       << (float)(cerTagInTime) / (float)(m_cerenkovEvts)*100 << "%)" << std::endl;

  m_lcWriter->close();

  // TCanvas *c1 = new TCanvas();
  // std::unique_ptr<TCanvas> c1 = make_unique<TCanvas>();
  // c1->SetCanvasSize(2048, 1024);
  // c1->Update();
  // c1->cd();
  // c1->Divide(2, 1);
  // c1->cd(1);

  // std::cout << "Drawing for layer 48" << std::endl;
  // m_vHitMapPerLayer.at(47)->Draw("colz");
  // c1->cd(2);
  // std::cout << "Drawing for layer 50" << std::endl;
  // m_vHitMapPerLayer.at(49)->Draw("colz");
  // std::stringstream ss;
  // ss << m_plotFolder << "/hitMap_Layer48-50_run" << m_runNumber << ".png";
  // c1->SaveAs(ss.str().c_str());
  m_rootFile->cd();
  m_rootFile->Write();
  m_rootFile->Close();

  streamlog_out(MESSAGE) << "Trivent end" << std::endl;
}

//==============================================================
