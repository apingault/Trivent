// -- Trivent include
#include "TriventProc.hh"

// -- std includes
#include <algorithm> // std::max_element
#include <fstream>   // std::stringstream
#include <iterator>  // std::next

// -- Root headers
#include <TCanvas.h>

TriventProc triventForCommonEcal;

// -- json
#include <nlohmann/json.hpp>
using json = nlohmann::json;

//=========================================================
TriventProc::TriventProc() : Processor("TriventProc") {
  // collection
  registerInputCollections(LCIO::RAWCALORIMETERHIT, "InputCollectionNames", "HCAL Collection Names", m_inputCollections,
                           m_inputCollections);

  // layer cut
  registerProcessorParameter("LayerCut", "cut in number of layer 10 in default", m_layerCut, m_layerCut);

  // noise cut
  registerProcessorParameter("NoiseCut", "noise cut in time spectrum 10 in default", m_noiseCut, m_noiseCut);

  // time windows
  registerProcessorParameter("TimeWin", "time window = 2 in default", m_timeWin, m_timeWin);

  // dif mapping in json file
  registerProcessorParameter("SetupGeometry", "Dif geometry and position on the detector in json", m_geomFile,
                             m_geomFile);

  // electronic noise cut
  registerProcessorParameter("ElectronicNoiseCut", "number of hit max on time stamp", m_elecNoiseCut, m_elecNoiseCut);

  // electronic noise cut
  registerProcessorParameter("_time2prev_event_cut", "cut on time to previous event (x 200 ns)", m_time2prevEventCut,
                             m_time2prevEventCut);

  // Root Tree
  registerProcessorParameter("TreeName", "Logroot tree name", m_treeName, m_treeName);

  // Root Tree
  registerProcessorParameter("TreeDescription", "Logroot tree description", m_treeDescription, m_treeDescription);

  // histogram control tree
  registerProcessorParameter("ROOTOutputFile", "Logroot name", m_rootFileName, m_rootFileName);

  registerProcessorParameter("HasCerenkovDIF", "If Cerenkov dif was connected during data taking", m_hasCherenkov,
                             m_hasCherenkov);

  registerProcessorParameter("CerenkovDifId", "Dif number for cerenkov data", m_cerenkovDifId, m_cerenkovDifId);

  registerProcessorParameter("CellIdFormat",
                             "Data Format string: it could be M:3,S-1:3,I:9,J:9,K-1:6 (ILD_ENDCAP) or "
                             "I:9,J:9,K-1:6,Dif_id:8,Asic_id:6,Chan_id:7",
                             m_cellIdFormat, m_cellIdFormat);

  registerProcessorParameter("ZShift", "Shift on Z axis to compute hitPosition", m_zShift, m_zShift);
}

//=============================================================================
// insert the new dif and check if it already exists elsewhere
void TriventProc::insertDifIntoMap(int difId, difGeom &dif) {
  auto it = m_difMapping.insert({difId, dif});
  if (!it.second) {
    throw std::runtime_error("ERROR in geometry : dif " + std::to_string(difId) + " of layer " +
                             std::to_string(dif.layerId) + " already exists");
  }
}

//=============================================================================
void TriventProc::readGeometry(const std::string &geomFile) {

  std::ifstream jsonStream(geomFile);
  const auto    jsonFile  = json::parse(jsonStream);
  const auto    layerList = jsonFile.at("chambers");

  for (const auto &layer : layerList) {
    const int slotId = layer.at("slot");
    int       difId  = layer.at("left");
    difGeom   temp{slotId, 0};
    insertDifIntoMap(difId, temp);

    difId       = layer.at("center");
    temp.shiftY = 32;
    insertDifIntoMap(difId, temp);

    difId       = layer.at("right");
    temp.shiftY = 64;
    insertDifIntoMap(difId, temp);

    streamlog_out(DEBUG) << "inserting layer " << slotId << std::endl;

    // Make sure dummy layer for the bif is not in the list of layer (some geometry files implement this...)
    if (m_hasCherenkov && difId == m_cerenkovDifId) {
      throw std::runtime_error(
          "Bif should not be present in this part of the geometry, please move it to a separate 'bifId' section");
    }
    m_layerSet.insert(slotId);
  }

  const auto difsToSkipList = jsonFile.at("difsToSkip");
  for (const auto &difItem : difsToSkipList) {
    m_difsToSkip.emplace_back(difItem);
  }

  if (m_hasCherenkov) {
    // make sure we didn't already added the dif in the mapping
    const auto difIter = m_difMapping.find(m_cerenkovDifId);
    m_cerenkovLayerId  = *(m_layerSet.rbegin()) + 100; // Set dummy layerId for the bif to be 100 layers after the end
                                                       // of the prototype, so it doesn't register in the displays

    if (difIter == m_difMapping.end()) {
      try {
        const int bifId = jsonFile.at("bifId").get<int>();
        if (bifId != m_cerenkovDifId) {
          throw std::runtime_error("BifId from geometry file (" + std::to_string(bifId) +
                                   ") is different than marlin xml parameter (" + std::to_string(m_cerenkovDifId) +
                                   ")");
        }
        difGeom bifGeom{m_cerenkovLayerId, 0, 0, 1}; //
        insertDifIntoMap(bifId, bifGeom);
      } catch (std::exception &e) { // TODO: having to declare the bif id in geomFile + config file is dumb af!
        streamlog_out(WARNING) << "[" << __func__
                               << "] - No bifId (cerenkov dif) found in geometry file, using default parameter from "
                                  "marlin config file : "
                               << m_cerenkovDifId << std::endl;
        difGeom bifGeom{m_cerenkovLayerId, 0, 0, 1};
        insertDifIntoMap(m_cerenkovDifId, bifGeom);
      }
    } else { // should never happen now
      throw std::runtime_error("No reason the bif dummy layer should be present");
    }
  }

  streamlog_out(MESSAGE) << yellow << "[" << __func__ << "] - Found " << m_difMapping.size()
                         << " difs in geometry file corresponding to : "
                         << static_cast<unsigned int>(m_difMapping.size()) / 3 << " layers + "
                         << m_difMapping.size() % 3 << " difs" << normal << std::endl;
}

//=============================================================================
void TriventProc::printDifGeom() const {
  streamlog_out(DEBUG1) << "[" << __func__ << "] --- Dumping geomtry File: " << std::endl;
  for (const auto &dif : m_difMapping) {
    streamlog_out(DEBUG1) << dif.first << "\t" << dif.second.layerId << "\t" << dif.second.shiftY << "\t"
                          << dif.second.shiftX << "\t" << dif.second.nAsics << std::endl;
  }
}

// ============ decode the cell ids =============
// bit shift & 0xFF = Apply mask 1111 1111 to binary value
// eg: Dif 1 => cellID0 = 00983297 => DifID = 1 / AsicID = 1 / ChanID = 15
int TriventProc::getCellDif_id(const int cellId) const { return cellId & 0xFF; }

//=============================================================================
//  bit shift & 0xFF00 Apply mask 1111 1111 0000 0000 then cut last 8 bits
int TriventProc::getCellAsic_id(const int cellId) const { return (cellId & 0xFF00) >> 8; }

//=============================================================================
//  bit shift & 0x3F0000 Apply mask 1111 0000 0000 0000 0000 then cut last 16 bits
int TriventProc::getCellChan_id(const int cellId) const { return (cellId & 0x3F0000) >> 16; }

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

bool TriventProc::checkPadLimits(const std::vector<int> &padIndex, const std::vector<int> &padLimits) const {
  assert(padLimits.size() == padIndex.size() * 2);
  for (int i = 0; i < static_cast<int>(padIndex.size()); ++i) {
    if (padIndex[i] < padLimits[i * 2] || padIndex[i] > padLimits[i * 2 + 1]) {
      return false;
    }
  }
  return true;
}

//=============================================================================
std::vector<int> TriventProc::getPadIndex(const int difId, const int asicId, const int chanId) const {

  const auto findIter = m_difMapping.find(difId);
  if (findIter == m_difMapping.end()) {
    const auto findSkipIter = std::find(m_difsToSkip.begin(), m_difsToSkip.end(), difId);
    if (findSkipIter == m_difsToSkip.end()) {
      streamlog_out(ERROR) << "[" << __func__ << "] - difId '" << difId << "' not found in geometry file" << std::endl;
    } else
      streamlog_out(DEBUG) << "[" << __func__ << "] - difId '" << difId << "' is set to be skipped in geometry file"
                           << std::endl;
    return {}; // empty
  }

  std::vector<int> index{
      static_cast<int>(1 + MapILargeHR2.at(chanId) + AsicShiftI.at(asicId)),
      static_cast<int>(32 - (MapJLargeHR2.at(chanId) + AsicShiftJ.at(asicId)) + findIter->second.shiftY),
      static_cast<int>(findIter->second.layerId)};
  std::vector<int> padLims = {1, 96, 1, 96, 0, static_cast<int>(m_layerSet.size())};

  // Cerenkov layer is not in the layerSet as it's not a physical layer, needs to account for that when checking the pad
  // limits
  //
  if (difId == m_cerenkovDifId) {
    padLims.pop_back();
    padLims.push_back(m_cerenkovLayerId);
  }

  const bool padOk = checkPadLimits(index, padLims);
  assert(padOk);
  return index;
}

//=============================================================================
int TriventProc::getMaxTime() const {
  assert(!m_triggerRawHitMap.empty());
  return m_triggerRawHitMap.rbegin()->first;
}

//=============================================================================
std::vector<int> TriventProc::getTimeSpectrum(const int maxTime) const {
  std::vector<int> timeSpectrumVec(maxTime + 1, 0);
  for (const auto &mapIt : m_triggerRawHitMap) {
    int time = mapIt.first;
    for (const auto &rawIt : mapIt.second) {
      assert(time == rawIt->getTimeStamp());
    }
    assert(time <= maxTime);
    assert(time >= 0);
    if (time >= 0) {
      timeSpectrumVec.at(time) += mapIt.second.size();
    }
  }
  return timeSpectrumVec;
}

//=============================================================================
int TriventProc::IJKToKey(const std::vector<int> &padIndex) const {
  return 100 * 100 * padIndex[2] + 100 * padIndex[1] + padIndex[0];
}

//=============================================================================
int TriventProc::getAsicKey(const std::vector<int> &padIndex) const {
  // Not necessary to check for boundary here as already tested in getPadIndex
  const int jnum = (padIndex[1] - 1) / 8;
  const int inum = (padIndex[0] - 1) / 8;
  const int num  = jnum * 12 + inum;
  return padIndex[2] * 1000 + num;
}

//=============================================================================
void TriventProc::resetTriggerParameters() { m_triggerRawHitMap.clear(); }

//=============================================================================
void TriventProc::resetEventParameters() {
  // reset Flags
  m_nHit = 0;
  m_hitI.clear();
  m_hitJ.clear();
  m_hitK.clear();
  m_hitX.clear();
  m_hitY.clear();
  m_hitZ.clear();
  m_hitCogX = 0;
  m_hitCogY = 0;
  m_hitCogZ = 0;
  m_hitThreshold.clear();
  m_hitBcid.clear();
  m_hitRevBcid.clear();
  m_firedLayersSet.clear();
  m_nFiredLayers = 0;
}

//=============================================================================
void TriventProc::eventBuilder(const int timePeak, const int lowTimeBoundary, const int highTimeBoundary) {

  resetEventParameters();

  std::map<int, int> asicMap{};
  std::map<int, int> ramFullMap{};
  std::map<int, int> hitKeys{};

  for (int hitTime = lowTimeBoundary; hitTime <= highTimeBoundary; ++hitTime) {

    // No hit recorded at hitTime
    if (m_triggerRawHitMap.find(hitTime) == m_triggerRawHitMap.end()) {
      continue;
    }

    for (const auto &rawHit : m_triggerRawHitMap.at(hitTime)) {
      assert(rawHit);
      assert(rawHit->getTimeStamp() == hitTime);

      const int difId = getCellDif_id(rawHit->getCellID0());
      assert(difId != m_cerenkovDifId); // No hit from the bif should be present in this collection

      const int asicId = getCellAsic_id(rawHit->getCellID0());
      const int chanId = getCellChan_id(rawHit->getCellID0());
      const int thresh = rawHit->getAmplitude();
      if (asicId < 1 || asicId > 48) {
        throw std::runtime_error("Found a hit with weird AsicId, Dif/Asic/Chan/Thr... " + std::to_string(difId) + "/" +
                                 std::to_string(asicId) + "/" + std::to_string(chanId) + "/" + std::to_string(thresh));
      }
      if (chanId < 0 || chanId > 63) {
        throw std::runtime_error("Found a hit with weird ChannelId, Dif/Asic/Chan/Thr... " + std::to_string(difId) +
                                 "/" + std::to_string(asicId) + "/" + std::to_string(chanId) + "/" +
                                 std::to_string(thresh));
      }

      // find and remove ramFull events
      if (chanId == 29 || chanId == 31) {
        ramFullMap[difId]++;
      }

      if (ramFullMap[difId] > 48) {
        streamlog_out(DEBUG1) << yellow << "[" << __func__ << "] - Rejecting event with ram full. Dif : " << difId
                              << normal << std::endl;
        ++m_rejectedNum;
        hitKeys.clear();
        asicMap.clear();
        ramFullMap.clear();
        return;
      }

      // find and remove square events
      const std::vector<int> padIndex = getPadIndex(difId, asicId, chanId);
      if (padIndex.empty()) {
        continue;
      }

      const int asicKey = getAsicKey(padIndex);

      if (m_layerStartAt0)
        assert(asicKey >= 0);
      else
        assert(asicKey > 0);

      if (asicMap[asicKey] != 0) {
        ++asicMap[asicKey];
      } else {
        asicMap[asicKey] = 1;
      }

      if (asicMap[asicKey] == 64) {
        streamlog_out(DEBUG1) << yellow << "[" << __func__ << "] - Rejecting event with full asic. Dif '" << difId
                              << "' asic '" << asicId << "' at time '" << timePeak << "'" << normal << std::endl;
        ++m_rejectedNum;
        hitKeys.clear();
        asicMap.clear();
        ramFullMap.clear();
        return;
      }

      // Create hit Key
      const int aHitKey = IJKToKey(padIndex);

      // Avoid two hit in the same cell
      std::map<int, int>::const_iterator findIter = hitKeys.find(aHitKey);
      if (findIter != hitKeys.end()) {
        continue;
      }

      // Creating Calorimeter Hit
      const int            I = padIndex[0];
      const int            J = padIndex[1];
      const int            K = padIndex[2];
      std::array<float, 3> pos{static_cast<float>(J * m_cellSizeJ + m_cellSizeJ / 2),
                               static_cast<float>((96 - I) * m_cellSizeI + m_cellSizeI / 2),
                               static_cast<float>(K * m_layerThickness + m_zShift)};

      // add layer to list of unique touched layers
      m_firedLayersSet.insert(K);

      hitKeys.insert(std::pair<int, int>(aHitKey, hitTime));
      m_hitI.push_back(I);
      m_hitJ.push_back(J);
      m_hitK.push_back(K);
      m_hitX.push_back(pos.at(0));
      m_hitY.push_back(pos.at(1));
      m_hitZ.push_back(pos.at(2));
      m_hitBcid.push_back(rawHit->getTimeStamp());
      m_hitRevBcid.push_back(m_acquisitionTime - rawHit->getTimeStamp());
      m_hitThreshold.push_back(thresh);
    }
  }
  m_nFiredLayers = m_firedLayersSet.size();
}

//=============================================================================
TTree *TriventProc::getOrCreateTree(const std::string &treeName, const std::string &treeDescription) const {
  auto *tree = dynamic_cast<TTree *>(m_rootFile->Get(treeName.c_str()));

  if (tree == nullptr) {
    streamlog_out(DEBUG0) << "[" << __func__ << "] - Creating tree '" << treeName << "'" << std::endl;
    tree = new TTree(treeName.c_str(), treeDescription.c_str());
  }
  assert(tree);
  return tree;
}

//=============================================================================
void TriventProc::initRootTree() {
  m_rootFile = new TFile(m_rootFileName.c_str(), "RECREATE");
  assert(m_rootFile);

  // Create Event tree & Branches
  m_eventTree = getOrCreateTree(m_treeName, m_treeDescription);
  assert(m_eventTree);
  m_eventTree->Branch("DetId", &m_detId);
  m_eventTree->Branch("TrigNum", &m_evtTrigNbr);
  m_eventTree->Branch("TrigBcid", &m_triggerBcid);
  m_eventTree->Branch("TrigLength", &m_acquisitionTime);
  m_eventTree->Branch("EvtNum", &m_evtNbr);
  m_eventTree->Branch("EvtBcid", &m_evtBcid);
  m_eventTree->Branch("EvtRevBcid", &m_evtReversedBcid);
  m_eventTree->Branch("NHits", &m_nHit);
  m_eventTree->Branch("HitI", &m_hitI);
  m_eventTree->Branch("HitJ", &m_hitJ);
  m_eventTree->Branch("HitK", &m_hitK);
  m_eventTree->Branch("HitX", &m_hitX);
  m_eventTree->Branch("HitY", &m_hitY);
  m_eventTree->Branch("HitZ", &m_hitZ);
  m_eventTree->Branch("HitCogX", &m_hitCogX);
  m_eventTree->Branch("HitCogY", &m_hitCogY);
  m_eventTree->Branch("HitCogZ", &m_hitCogZ);
  m_eventTree->Branch("HitBcid", &m_hitBcid);
  m_eventTree->Branch("HitRevBcid", &m_hitRevBcid);
  m_eventTree->Branch("HitThresh", &m_hitThreshold);

  // Check if first layer is numbered 0 or 1
  // Prevent accessing non defined element in vectors...
  const auto firstLayer = std::min_element(m_layerSet.begin(), m_layerSet.end());
  if (0 == *firstLayer) {
    m_layerStartAt0 = true;
  }
}

//=============================================================================
void TriventProc::init() {
  printParameters();

  // Read and print geometry file
  streamlog_out(DEBUG) << "Reading geometry...\n";
  readGeometry(m_geomFile);
  streamlog_out(DEBUG) << "Reading geometry...DONE\n";
  printDifGeom();
  initRootTree();
}

//=============================================================================
void TriventProc::fillRawHitTrigger(const LCCollection &inputLCCol) {
  std::vector<int> vTrigger;

  for (int ihit(0); ihit < inputLCCol.getNumberOfElements(); ++ihit) // loop over the hits
  {
    auto *rawHit = dynamic_cast<RawCalorimeterHit *>(inputLCCol.getElementAt(ihit));
    if (!rawHit) {
      continue;
    }
    // extract abolute bcid information:
    const int difId = rawHit->getCellID0() & 0xFF;
    assert(difId > 0);
    if (ihit == 0) {
      std::stringstream pname("");
      pname << "DIF" << difId << "_Triggers";
      inputLCCol.getParameters().getIntVals(pname.str(), vTrigger);
      if (!vTrigger.empty()) {
        m_bcid1              = vTrigger[4];
        m_bcid2              = vTrigger[3];
        const uint64_t Shift = 16777216ULL; // to shift the value from the 24 first bits
        m_triggerBcid        = m_bcid1 * Shift + m_bcid2;
        m_acquisitionTime    = vTrigger[2]; // in 200ns clock

        streamlog_out(DEBUG0) << "[" << __func__ << "] - trigger time : " << m_triggerBcid << std::endl;
      }
    }

    // Weird timestamp from some dif from time to time
    if (rawHit->getTimeStamp() < 0) {
      // streamlog_out(ERROR) << red << "[" << __func__ << "] - Trig '" << m_trigNbr
      streamlog_out(DEBUG) << red << "[" << __func__ << "] - Trig '" << m_trigNbr
                           << "'Found a raw hit with negative timeStamp! : "
                           << "time: " << rawHit->getTimeStamp() << " difId: " << getCellDif_id(rawHit->getCellID0())
                           << " asicId: " << getCellAsic_id(rawHit->getCellID0())
                           << " chanId: " << getCellChan_id(rawHit->getCellID0())
                           << " thresh: " << rawHit->getAmplitude() << " removing it !" << normal << std::endl;
      continue;
    }
    if (difId != m_cerenkovDifId) {
      m_triggerRawHitMap[rawHit->getTimeStamp()].push_back(rawHit);
    }
  }
}

//=============================================================================
std::vector<std::vector<int>::iterator>
TriventProc::getCandidateTimeBoundaries(const std::vector<int>::iterator &beginTime,
                                        const std::vector<int>::iterator &endTime,
                                        const std::vector<int>::iterator &candidateTime) const {
  assert(beginTime < endTime);
  assert(candidateTime >= beginTime);
  assert(candidateTime < endTime); // Can't be equal otherwise we'll access out of bound value on next ++timeIter call
  auto lowerBound = endTime;
  auto upperBound = beginTime;

  // Ensure there is sufficient time between two candidate + we are not looking before begining of timeSpectrumVec
  auto timeDistance = std::distance(beginTime, candidateTime);
  if (timeDistance > m_timeWin) {
    lowerBound = std::prev(candidateTime, m_timeWin);
    // Don't throw a potential candidate found in the first few frames of the trigger
  } else {
    streamlog_out(DEBUG0) << green << "[" << __func__ << "] - Small lowerBound! m_timeWin : " << m_timeWin
                          << " distance(beginTime, candidateTime) = " << timeDistance << normal << std::endl;
    lowerBound = std::prev(candidateTime, timeDistance);
  }

  // Check we are sufficiently far from end of timeSpectrumVec
  timeDistance = std::distance(candidateTime, endTime);
  if (timeDistance > m_timeWin) {
    upperBound = std::next(candidateTime, m_timeWin);
  } else { // Don't throw a potential candidate found in the last few frames of the trigger
    upperBound = std::next(candidateTime, timeDistance); // distance > 0 already met in while loop
    streamlog_out(DEBUG0) << green << "[" << __func__ << "] - Small upperBound! m_timeWin : " << m_timeWin
                          << " distance(candidateTime, endTime) = " << timeDistance << normal << std::endl;
  }
  streamlog_out(DEBUG0) << "[" << __func__ << "] - low: " << distance(beginTime, lowerBound) << " up "
                        << distance(beginTime, upperBound) << std::endl;
  assert(lowerBound < upperBound);
  assert(std::distance(upperBound, lowerBound) <= 2 * m_timeWin);
  return {lowerBound, upperBound};
}

//=============================================================================
void TriventProc::processEvent(LCEvent *evtP) {
  assert(evtP != nullptr);

  m_trigNbr = evtP->getEventNumber();
  if (m_trigNbr > 1E6) {
    streamlog_out(ERROR) << red << "[" << __func__ << "] - Too much Triggers : " << m_trigNbr << normal << std::endl;
    return;
  }

  for (unsigned int i = 0; i < m_inputCollections.size(); i++) {
    LCCollection *inputLCCol{nullptr};
    try {
      inputLCCol = evtP->getCollection(m_inputCollections.at(i));
    } catch (lcio::DataNotAvailableException &zero) {
      streamlog_out(ERROR) << red << "[" << __func__ << "] - No data found in collection " << i << normal << std::endl;
    }

    if (inputLCCol == nullptr) {
      streamlog_out(WARNING) << red << "[" << __func__ << "] - TRIGGER SKIPED ... col is nullptr" << normal
                             << std::endl;
      continue;
    }

    ++m_trigCount;
    if (0 == m_trigCount % 100) {
      streamlog_out(MESSAGE) << yellow << "[" << __func__ << "] - Trigger number == " << m_trigCount << normal
                             << std::endl;
    }

    const int numElements = inputLCCol->getNumberOfElements(); // nHits in trigger
    if (numElements > m_elecNoiseCut) {
      streamlog_out(MESSAGE) << yellow << "[" << __func__ << "] - TRIGGER number " << m_trigCount
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
    auto beginTimeIter = timeSpectrumVec.begin();
    auto endTimeIter   = timeSpectrumVec.end();
    auto timeIter      = beginTimeIter;
    auto prevMaxIter   = beginTimeIter;

    while (std::distance(timeIter, endTimeIter) > 0) { // Ensure that timeIter < endTime
      if (*(timeIter) < m_noiseCut) {                  // Not enough hit in frame, look in next one
        ++timeIter;
        continue;
      }

      // find timeBoundaries to build the event
      const auto  boundaries = getCandidateTimeBoundaries(prevMaxIter, endTimeIter, timeIter);
      const auto &maxIter    = std::max_element(boundaries[0], boundaries[1]); // max in [lower,upper))

      if (maxIter > timeIter || maxIter == beginTimeIter) {
        // timeIter is not a real peak yet, look in next frame
        ++timeIter;
        continue;
      }

      streamlog_out(DEBUG0) << yellow << "[" << __func__ << "] - upperBound '"
                            << std::distance(beginTimeIter, boundaries[1]) << "' lowerBound '"
                            << std::distance(beginTimeIter, boundaries[0]) << " max at time '"
                            << std::distance(beginTimeIter, maxIter) << "' : " << *maxIter << normal << std::endl;

      // Check we didn't already process the peak
      // if (std::distance(maxIter, prevMaxIter) < m_timeWin && std::distance(beginTimeIter, timeIter) > 0) {
      if (maxIter < prevMaxIter && std::distance(beginTimeIter, timeIter) > 0) {
        if (std::distance(prevMaxIter, maxIter) < m_timeWin) {
          streamlog_out(DEBUG0) << yellow << "[" << __func__ << "] - Found duplicate peak, at time '"
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

      streamlog_out(DEBUG0) << blue << "[" << __func__ << "] - Trig '" << m_trigCount << "' : Found Peak, at time '"
                            << timePeak << "' - hits : " << *maxIter << " prevTimePeak: " << prevTimePeak << normal
                            << std::endl;

      std::unique_ptr<LCCollectionVec> outCol = std::make_unique<LCCollectionVec>(LCIO::CALORIMETERHIT);

      // Event Building
      streamlog_out(DEBUG0) << blue << "[" << __func__ << "] - EventBuilding..." << normal << std::endl;
      try {
        TriventProc::eventBuilder(timePeak, lowBound, highBound);
      } catch (const std::exception &e) {
        throw;
      }

      streamlog_out(DEBUG0) << blue << "[" << __func__ << "] - EventBuilding...OK" << normal << std::endl;
      // Apply cut on min number of firedLayer +
      if (static_cast<int>(m_nFiredLayers) < m_layerCut) {
        streamlog_out(DEBUG0) << green << "[" << __func__
                              << "] - Event rejected, too few layer hit. nLayerHit: " << m_nFiredLayers
                              << " m_layerCut: " << m_layerCut << normal << std::endl;
      }

      m_evtTrigNbr      = m_trigNbr;
      m_evtNbr          = m_evtNum; // dont increment here: rejected event will have same number as last accepted !
      m_evtBcid         = timePeak;
      m_evtReversedBcid = m_acquisitionTime - timePeak;
      m_nHit            = outCol->getNumberOfElements();
      m_hitCogX         = std::accumulate(m_hitX.begin(), m_hitX.end(), 0.0) / m_hitX.size();
      m_hitCogY         = std::accumulate(m_hitY.begin(), m_hitY.end(), 0.0) / m_hitY.size();
      m_hitCogZ         = std::accumulate(m_hitZ.begin(), m_hitZ.end(), 0.0) / m_hitZ.size();

      //  Apply cut on time between two events
      assert(timePeak - prevTimePeak > m_time2prevEventCut);

      streamlog_out(DEBUG0) << green << "[" << __func__ << "] - Trivent find event at :==> " << red << timePeak << green
                            << "\t :Nhit: ==> " << magenta << outCol->getNumberOfElements() << normal << std::endl;

      ++m_selectedNum;
      m_eventTree->Fill();
      timeIter = std::next(timeIter, m_timeWin + 1);
    }
  }
}

//=============================================================================
void TriventProc::end() {

  streamlog_out(MESSAGE) << "Trivent rejected " << m_rejectedNum << " Candidate event" << std::endl;
  streamlog_out(MESSAGE) << "Trivent Selected " << m_selectedNum << " Candidate event" << std::endl;

  m_rootFile->cd();
  m_rootFile->Write();
  m_rootFile->Close();

  streamlog_out(MESSAGE) << "Trivent end" << std::endl;
}

//==============================================================
