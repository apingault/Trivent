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
  void processEvent(LCEvent *evtP);
  void processRunHeader(LCRunHeader * /*runH*/){};
  void XMLReader(std::string xmlfile);
  void readDifGeomFile(std::string geomfile);
  void printDifGeom();
  void defineColors();

  int getCellDif_id(int cell_id);
  int getCellAsic_id(int cell_id);
  int getCellChan_id(int cell_id);

  const int getMaxTime();

  const std::vector<int> getTimeSpectrum(const int &maxTime);

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

  bool checkPadLimits(std::vector<int> &padIndex, std::vector<int> &padLimits);

  std::vector<int> getPadIndex(const int &dif_id, const int &asic_id, const int &chan_id);
  void eventBuilder(std::unique_ptr<IMPL::LCCollectionVec> &col_event, int &time_peak, int &prev_time_peak);
  void end();

  TH2 *makeTH2(const std::string &title, const std::string &xTitle, const std::string &yTitle);

  // std::unique_ptr<TTree> getOrCreateTree(const std::string &treeName, const std::string &treeDescription);
  TTree *getOrCreateTree(const std::string &treeName, const std::string &treeDescription);
  void findCerenkovHits(const int timePeak);

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
  unsigned int m_cerenkovEvts;    // Number of events tagged with the cerenkov
  int          m_maxCerenkovTime; // CerenkovTime with most occurrence ->MPV for cerenkov timeShift

  unsigned int           m_trigNbr;
  unsigned int           m_trigCount;
  unsigned int           m_evtNum;
  unsigned int           m_selectedNum;
  unsigned int           m_rejectedNum;
  std::set<unsigned int> m_firedLayersSet;
  int                    m_bcid1;
  int                    m_bcid2;

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

  unsigned int m_nFiredLayers;       // Number of Layers touched in evt
  bool         m_isSelected;         // Event is selected/rejected
  bool         m_isNoise;            // If rejected, is it noise
  bool         m_isTooCloseInTime;   // If rejected, is it too close from previous evt
  bool         m_hasNotEnoughLayers; // If rejected, has not touched sufficient layers
  bool         m_hasFullAsic;        // If rejected, has full asics
};

#endif
