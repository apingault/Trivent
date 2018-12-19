#ifndef TRIVENTPROC_HH
#define TRIVENTPROC_HH

// -- Asics and channels mapping for sdhcal
#include "Mapping.hh"

// -- std includes
#include <map>
#include <string>

// -- marlin includes
#include "marlin/VerbosityLevels.h"
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
  Processor *newProcessor() override { return new TriventProc; }

  TriventProc();
  TriventProc(const TriventProc & /*proc*/)  = delete;
  TriventProc(const TriventProc && /*proc*/) = delete;
  ~TriventProc() override                    = default;

  TriventProc &operator=(const TriventProc & /*proc*/) = delete;
  TriventProc &operator=(const TriventProc && /*proc*/) = delete;

  void init() override;
  void processEvent(LCEvent *evtP) override;
  void end() override;

  void insertDifIntoMap(int difId, difGeom &dif);
  void readGeometry(const std::string &geomFile);
  void printDifGeom() const;
  void initRootTree();

  int getCellDif_id(const int cellId) const;
  int getCellAsic_id(const int cellId) const;
  int getCellChan_id(const int cellId) const;

  int getMaxTime() const;

  std::vector<int> getTimeSpectrum(const int maxTime) const;

  void fillRawHitTrigger(const LCCollection &inputLCCol);

  /**
   * @brief Return **inclusive** time boundaries for eventBuilding
   *
   * @param beginTime
   * @param endTime
   * @param candidateTime
   * @return std::vector<const std::vector<int>::iterator>
   */
  std::vector<std::vector<int>::iterator>
  getCandidateTimeBoundaries(const std::vector<int>::iterator &beginTime, const std::vector<int>::iterator &endTime,
                             const std::vector<int>::iterator &candidateTime) const;

  /**
   * @brief return I, J, K for hit in given difId, asicId, chanId
   * return empty vector if dif not found in geometry
   *
   * @param difId
   * @param asicId
   * @param chanId
   * @return std::vector<int>
   */
  std::vector<int> getPadIndex(const int difId, const int asicId, const int chanId) const;
  bool             checkPadLimits(const std::vector<int> &padIndex, const std::vector<int> &padLimits) const;

  void resetEventParameters();
  void resetTriggerParameters();

  void eventBuilder(const int timePeak, const int lowTimeBoundary, const int highTimeBoundary);

  TTree *getOrCreateTree(const std::string &treeName, const std::string &treeDescription) const;
  int    getAsicKey(const std::vector<int> &padIndex) const;
  int    IJKToKey(const std::vector<int> &padIndex) const;

private:
  // map of <hitTimeStamp, rawHit>
  std::map<int, std::vector<EVENT::RawCalorimeterHit *>> m_triggerRawHitMap{};

  std::string m_rootFileName{"TDHCAL.root"};
  std::string m_treeName{"EventTree"};
  std::string m_treeDescription{"Event variables"};
  std::string m_geomFile{"setup_geometry"};
  std::string m_cellIdFormat{"M:3,S-1:3,I:9,J:9,K-1:6"};

  std::vector<std::string> m_inputCollections{"DHCALRawHits"};
  std::map<int, difGeom>   m_difMapping{}; ///< map of <difId/difGeom>
  std::vector<int>         m_difsToSkip{};

  // Cut parameters
  int m_elecNoiseCut{5000};
  int m_noiseCut{10};
  int m_layerCut{10};
  int m_timeWin{2};
  int m_time2prevEventCut{0};

  // Geometry paramaters all in mm
  float         m_cellSizeI{10.408};
  float         m_cellSizeJ{10.408};
  float         m_layerThickness{26.131};
  float         m_zShift{225};
  std::set<int> m_layerSet{};
  bool          m_layerStartAt0{false};

  unsigned int m_trigNbr{0};
  unsigned int m_trigCount{0};
  unsigned int m_evtNum{0};
  unsigned int m_selectedNum{0};
  unsigned int m_rejectedNum{0};
  int          m_bcid1{0};
  int          m_bcid2{0};

  // Color for streamlog output
  std::string normal{0x1b, '[', '0', ';', '3', '9', 'm', 0};
  std::string red{0x1b, '[', '1', ';', '3', '1', 'm', 0};
  std::string green{0x1b, '[', '1', ';', '3', '2', 'm', 0};
  std::string yellow{0x1b, '[', '1', ';', '3', '3', 'm', 0};
  std::string blue{0x1b, '[', '1', ';', '3', '4', 'm', 0};
  std::string magenta{0x1b, '[', '1', ';', '3', '5', 'm', 0};
  std::string white{0x1b, '[', '1', ';', '3', '9', 'm', 0};

  // ROOT histograms
  TFile *m_rootFile{};
  TTree *m_eventTree{};

  // Event branches, all time length are in 5MHz clock
  unsigned int       m_detId{100};         // DetectorId (used for ecal)
  unsigned int       m_evtTrigNbr{0};      // Current trigger Number
  unsigned long long m_triggerBcid{0};     // Current trigger bcid
  unsigned long long m_acquisitionTime{0}; // Current trigger length
  unsigned int       m_evtNbr{0};          // Current Evt number
  unsigned int       m_evtBcid{0};         // Current Evt bcid
  unsigned int       m_evtReversedBcid{0}; // Current Evt reversed bcid = (acquisionTime - bcid)
  unsigned int       m_nHit{0};            // Number of hits

  std::vector<int>    m_hitBcid{};    // Hit time
  std::vector<int>    m_hitRevBcid{}; // Reversed Hit time
  std::vector<int>    m_hitI{};       // Hit position
  std::vector<int>    m_hitJ{};       // Hit position
  std::vector<int>    m_hitK{};       // Hit position
  std::vector<double> m_hitX{};       // Hit position
  std::vector<double> m_hitY{};       // Hit position
  std::vector<double> m_hitZ{};       // Hit position
  std::vector<int>    m_hitThreshold{};
  double              m_hitCogX{-1.}; // Hit Cog
  double              m_hitCogY{-1.}; // Hit Cog
  double              m_hitCogZ{-1.}; // Hit Cog

  std::set<unsigned int> m_firedLayersSet{}; // set of Layers touched in evt
  unsigned int           m_nFiredLayers{0};  // Number of Layers touched in evt = m_firedLayersSet.size()
};
#endif
