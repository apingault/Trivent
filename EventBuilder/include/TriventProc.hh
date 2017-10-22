#ifndef _TriventProc_hh_
#define _TriventProc_hh_

// -- Asics and channels mapping for sdhcal
#include "Mapping.h"

// -- std includes
#include <map>
#include <string>

// -- marlin includes
#include "marlin/VerbosityLevels.h"
#include "marlin/tinyxml.h"
#include <marlin/Processor.h>

// -- lcio includes
#include "IO/LCWriter.h"
#include <EVENT/LCCollection.h>
#include <EVENT/RawCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <lcio.h>

// -- ROOT includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TTree.h>

class TriventProc : public marlin::Processor {
public:
  Processor *newProcessor() { return new TriventProc; }

  TriventProc();
  ~TriventProc() { ; };

  void init();
  void initRootTree();
  void processEvent(LCEvent *evtP);
  void processRunHeader(LCRunHeader * /*runH*/){};
  void XMLReader(const std::string &xmlfile);
  void printDifGeom();
  void defineColors();

  int getCellDif_id(const int &cellId);
  int getCellAsic_id(const int &cellId);
  int getCellChan_id(const int &cellId);

  int getMaxTime();

  std::vector<int> getTimeSpectrum(const int &maxTime);

  /**
   * @brief
   *
   * @param inputLCCol
   */
  void fillRawHitTrigger(const LCCollection &inputLCCol);

  /**
   * @brief
   *
   * @param beginTime
   * @param endTime
   * @param candidateTime
   * @return std::vector<std::vector<int>::iterator>
   */
  std::vector<std::vector<int>::iterator> getCandidateTimeBoundaries(std::vector<int>::iterator &beginTime,
                                                                     std::vector<int>::iterator &endTime,
                                                                     std::vector<int>::iterator &candidateTime);

  /**
   * @brief Implementation of std::make_unique from c++14
   *
   * @tparam T
   * @tparam Args
   * @param args
   * @return std::unique_ptr<T>
   */
  template <typename T, typename... Args> std::unique_ptr<T> make_unique(Args &&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  bool checkPadLimits(const std::vector<int> &padIndex, const std::vector<int> &padLimits);

  /**
   * @brief return I, J, K for hit in given difId, asicId, chanId
   * return empty vector if dif not found in geometry
   *
   * @param difId
   * @param asicId
   * @param chanId
   * @return std::vector<int>
   */
  std::vector<int> getPadIndex(const int &difId, const int &asicId, const int &chanId);

  /**
   * @brief
   *
   */
  void resetEventParameters();

  /**
   * @brief
   *
   */
  void resetTriggerParameters();

  void eventBuilder(std::unique_ptr<IMPL::LCCollectionVec> &evtCol, const int &timePeak, const int &prevTimePeak);
  void end();

  TH2 *makeTH2(const std::string &title, const std::string &xTitle, const std::string &yTitle);

  // std::unique_ptr<TTree> getOrCreateTree(const std::string &treeName, const std::string &treeDescription);
  TTree *getOrCreateTree(const std::string &treeName, const std::string &treeDescription);
  void findCerenkovHits(const int &timePeak);
  int getAsicKey(const std::vector<int> &padIndex);
  int IJKToKey(const std::vector<int> &padIndex);

protected:
  std::unique_ptr<LCWriter>               m_lcWriter;
  std::vector<EVENT::RawCalorimeterHit *> m_trigger_raw_hit;
  std::vector<EVENT::RawCalorimeterHit *> m_cerenkov_raw_hit;
  // std::vector<std::shared_ptr<EVENT::RawCalorimeterHit>> m_trigger_raw_hit;
  // std::vector<std::shared_ptr<EVENT::RawCalorimeterHit>> m_cerenkov_raw_hit;

  std::string              m_outputCollectionName;
  std::string              m_outFileName;
  std::string              m_noiseFileName;
  std::string              m_rootFileName;
  std::string              m_treeName;
  std::string              m_treeDescription;
  std::string              m_geomXMLFile;
  std::vector<std::string> m_hcalCollections;
  float                    m_beamEnergy;
  std::map<int, LayerID> m_mDifMapping;

  // Cut parameters
  bool m_useGainCorrection;
  int  m_elecNoiseCut;
  int  m_noiseCut;
  int  m_layerCut;
  int  m_timeWin;
  int  m_time2prevEventCut;

  // Geometry paramaters
  float          m_cellSizeI;
  float          m_cellSizeJ;
  float          m_layerThickness;
  std::set<uint> m_layerSet;

  // Cerenkov
  bool         m_hasCherenkov;
  int          m_cerenkovDifId;
  int          m_cerenkovLayerId;
  int          m_cerenkovTimeWindow;
  unsigned int m_cerAsic;
  unsigned int m_cerChan;
  unsigned int m_cerThreshold;
  unsigned int m_nCerenkov1;         // Number of hit in first Cerenkov
  unsigned int m_nCerenkov2;         // Number of hit in second Cerenkov
  unsigned int m_nCerenkov3;         // Number of hit in first + second Cerenkov
  unsigned int m_nCerenkovTrigger;   // Tot number of hit in cerenkov for current trigger
  bool         m_hasTooManyCerenkov; // if m_nCerenkovTrigger > bifHit in trigger
  int          m_timeCerenkov;       // Timing between peak and Cerenkov signal
  unsigned int m_totCerenkovHits;
  unsigned int m_cerenkovEvts; // Number of events tagged with the cerenkov

  unsigned int m_trigNbr;
  unsigned int m_trigCount;
  unsigned int m_evtNum;
  unsigned int m_selectedNum;
  unsigned int m_rejectedNum;
  int          m_bcid1;
  int          m_bcid2;

  // Color for streamlog output
  std::string normal;
  std::string red;
  std::string green;
  std::string yellow;
  std::string blue;
  std::string magenta;
  std::string white;

  // ROOT histograms
  // std::unique_ptr<TFile> m_rootFile;
  TFile *m_rootFile;
  // std::vector<std::unique_ptr<TH2>> m_vHitMapPerLayer; // HitMap of selected evt for each Layer
  std::vector<TH2 *> m_vHitMapPerLayer; // HitMap of selected evt for each Layer
  unsigned int       m_runNumber;
  std::string        m_plotFolder;

  // Trees
  // TTree *m_triggerTree;
  // std::unique_ptr<TTree> m_triggerTree;
  TTree *m_eventTree;
  // std::unique_ptr<TTree> m_eventTree;

  // Trigger branches
  // unsigned int m_trigNbr;              // Current trigger number
  // int m_nEvt;                       // Number of evt in trigger
  // std::vector<unsigned int> m_vTimeSpectrum; // number of hits per time clock

  // Event branches
  unsigned int m_evtTrigNbr; // Current trigger Number
  unsigned int m_evtNbr;     // Current Evt number
  unsigned int m_nHit;       // Number of hits

  // std::vector<unsigned long int> m_hitBCID;             // Hit time
  std::vector<int> m_hitI; // Hit position
  std::vector<int> m_hitJ; // Hit position
  std::vector<int> m_hitK; // Hit position
  std::vector<int> m_hitThreshold;

  std::set<unsigned int> m_firedLayersSet;     // set of Layers touched in evt
  unsigned int           m_nFiredLayers;       // Number of Layers touched in evt = m_firedLayersSet.size()
  bool                   m_isSelected;         // Event is selected/rejected
  bool                   m_isNoise;            // If rejected, is it noise
  bool                   m_isTooCloseInTime;   // If rejected, is it too close from previous evt
  bool                   m_hasNotEnoughLayers; // If rejected, has not touched sufficient layers
  bool                   m_hasFullAsic;        // If rejected, has full asics
};

#endif
