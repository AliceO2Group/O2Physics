// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
// TableMaker produces skimmed data using the DQ data model
// Events to be written are filtered using a user provided event cut and optionally the filterPP task
// Barrel and muon tracks are filtered using multiple parallel selections (currently limited to 8)
// The skimming can optionally produce just the barrel, muon, or both barrel and muon tracks
// The event filtering (filterPP), centrality, and V0Bits (from v0-selector) can be switched on/off by selecting one
//  of the process functions
// C++ includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
// other includes
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/Zorro.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FwdTrackReAlignTables.h"
#include "Common/DataModel/MftmchMatchingML.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "DetectorsVertexing/VertexTrackMatcher.h"
#include "Field/MagneticField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Primitive2D.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"

#include "TGeoGlobalMagField.h"

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

Zorro zorro;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCovOnlyStdPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullEl, aod::pidTPCFullPi,
                                                  aod::pidTPCFullKa, aod::pidTPCFullPr,
                                                  aod::pidTOFFullEl, aod::pidTOFFullPi,
                                                  aod::pidTOFFullKa, aod::pidTOFFullPr>;
using MyBarrelTracksWithV0Bits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                           aod::pidTPCFullKa, aod::pidTPCFullPr,
                                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                           aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::V0Bits>;
using MyBarrelTracksWithV0BitsForMaps = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA,
                                                  aod::pidTPCFullEl, aod::pidTPCFullPi,
                                                  aod::pidTPCFullKa, aod::pidTPCFullPr, aod::V0Bits>;
using MyBarrelTracksWithDalitzBits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                               aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                               aod::pidTPCFullKa, aod::pidTPCFullPr,
                                               aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                               aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::DalitzBits>;
using MyBarrelTracksWithV0AndDalitzBits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                                    aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                    aod::pidTPCFullKa, aod::pidTPCFullPr,
                                                    aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                                    aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::V0Bits, aod::DalitzBits>;
using MyBarrelTracksForElectronMuon = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA,
                                                aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra>;
using MyEventsWithFilter = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventFilter>;
using MyEventsWithMultsAndFilter = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::DQEventFilter>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0As, aod::CentFT0Ms>;
using MyEventsWithCentAndMults = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0As, aod::CentFT0Ms, aod::Mults, aod::MultsExtra>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA>;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsCollWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsRealignCollWithCov = soa::Join<aod::FwdTracksReAlign, aod::FwdTrksCovReAlign, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkEventFillMapWithMult = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::CollisionMultExtra;
constexpr static uint32_t gkEventFillMapWithMultsAndEventFilter = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::CollisionMultExtra | VarManager::ObjTypes::EventFilter;
constexpr static uint32_t gkEventFillMapWithCent = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkEventFillMapWithCentAndMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent | VarManager::CollisionMult | VarManager::ObjTypes::CollisionMultExtra;
// constexpr static uint32_t gkEventFillMapWithCentRun2 = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCentRun2; // Unused variable
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
constexpr static uint32_t gkTrackFillMapWithCovOnlyStdPID = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithV0Bits = gkTrackFillMap | VarManager::ObjTypes::TrackV0Bits;
constexpr static uint32_t gkTrackFillMapWithV0BitsForMaps = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackV0Bits | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackTPCPID;
constexpr static uint32_t gkTrackFillMapWithDalitzBits = gkTrackFillMap | VarManager::ObjTypes::DalitzBits;
constexpr static uint32_t gkTrackFillMapWithV0AndDalitzBits = gkTrackFillMap | VarManager::ObjTypes::TrackV0Bits | VarManager::ObjTypes::DalitzBits;
constexpr static uint32_t gkTrackFillMapForElectronMuon = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackTPCPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;
constexpr static uint32_t gkMuonFillMapWithCovAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov | VarManager::ObjTypes::AmbiMuon;
constexpr static uint32_t gkMuonRealignFillMapWithCovAmbi = VarManager::ObjTypes::MuonRealign | VarManager::ObjTypes::MuonCovRealign | VarManager::ObjTypes::AmbiMuon;
constexpr static uint32_t gkTrackFillMapWithAmbi = VarManager::ObjTypes::Track | VarManager::ObjTypes::AmbiTrack;
constexpr static uint32_t gkMFTFillMap = VarManager::ObjTypes::TrackMFT;

struct TableMaker {

  Produces<ReducedEvents> event;
  Produces<ReducedEventsExtended> eventExtended;
  Produces<ReducedEventsVtxCov> eventVtxCov;
  Produces<ReducedEventsInfo> eventInfo;
  Produces<ReducedEventsMultPV> multPV;
  Produces<ReducedEventsMultAll> multAll;
  Produces<ReducedTracksBarrelInfo> trackBarrelInfo;
  Produces<ReducedTracks> trackBasic;
  Produces<ReducedTracksBarrel> trackBarrel;
  Produces<ReducedTracksBarrelCov> trackBarrelCov;
  Produces<ReducedTracksBarrelPID> trackBarrelPID;
  Produces<ReducedMuons> muonBasic;
  Produces<ReducedMuonsExtra> muonExtra;
  Produces<ReducedMuonsCov> muonCov;
  Produces<ReducedMuonsInfo> muonInfo;
  Produces<ReducedMFTs> trackMFT;
  Produces<ReducedMFTsExtra> trackMFTExtra;

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics
  HistogramManager* fHistMan;

  struct : ConfigurableGroup {
    Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
    Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
    Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  } configCuts;
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
  } addHistoConfigurations;
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 1.0f, "Low pt cut for tracks in the barrel"};
  Configurable<float> fConfigBarrelTrackMaxAbsEta{"cfgBarrelMaxAbsEta", 0.9f, "Eta absolute value cut for tracks in the barrel"};
  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 1.0f, "Low pt cut for muons"};
  struct : ConfigurableGroup {
    Configurable<float> fConfigMinTpcSignal{"cfgMinTpcSignal", 30.0, "Minimum TPC signal"};
    Configurable<float> fConfigMaxTpcSignal{"cfgMaxTpcSignal", 300.0, "Maximum TPC signal"};
  } configTpcSignal;
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes)"};
  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};
  Configurable<bool> fIsAmbiguous{"cfgIsAmbiguous", false, "Whether we enable QA plots for ambiguous tracks"};

  struct : ConfigurableGroup {
    Configurable<bool> fConfigRunZorro{"cfgRunZorro", false, "Enable event selection with zorro"};
    Configurable<std::string> fConfigZorroTrigMask{"cfgZorroTriggerMask", "fDiMuon", "DQ Trigger masks: fSingleE,fLMeeIMR,fLMeeHMR,fDiElectron,fSingleMuLow,fSingleMuHigh,fDiMuon"};
    Configurable<bool> fConfigRunZorroSel{"cfgRunZorroSel", false, "Select events with trigger mask"};
    Configurable<uint64_t> fBcTolerance{"cfgBcTolerance", 100, "Number of BCs of margin for software triggers"};
  } useZorro;

  struct : ConfigurableGroup {
    Configurable<std::string> fConfigCcdbUrl{"useCCDBConfigurations.ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> fConfigCcdbPathTPC{"useCCDBConfigurations.ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
    Configurable<std::string> fConfigCcdbPathZorro{"useCCDBConfigurations.ccdb-path-zorro", "/Users/m/mpuccio/EventFiltering/OTS/", "base path to the ccdb object for zorro"};
  } useCCDBConfigurations;

  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas(electrons, pions, protons)"};
  Configurable<bool> fConfigComputeTPCpostCalibKaon{"cfgTPCpostCalibKaon", false, "If true, compute TPC post-calibrated n-sigmas for kaons"};
  Configurable<bool> fConfigIsOnlyforMaps{"cfgIsforMaps", false, "If true, run for postcalibration maps only"};
  Configurable<bool> fPropMuon{"cfgPropMuon", false, "Propgate muon tracks through absorber"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpmagPathRun2{"grpmagPathRun2", "GLO/GRP/GRP", "CCDB path of the GRPObject (Usage for Run 2)"};
  struct : ConfigurableGroup {
    Configurable<int> useMatCorrType{"useMatConfigurations.useMatCorrType", 1, "materialCorrType: 0: none, 1: TGeo, 2: LUT"};
    Configurable<std::string> lutPath{"useMatConfigurations.lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  } useMatConfigurations;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  o2::parameters::GRPObject* grpmagrun2 = nullptr; // for run 2, we access the GRPObject from GLO/GRP/GRP
  o2::parameters::GRPMagField* grpmag = nullptr;   // for run 3, we access GRPMagField from GLO/Config/GRPMagField
  o2::base::MatLayerCylSet* lut = nullptr;

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fMuonCuts;  //! Muon track cuts

  Preslice<MyBarrelTracks> perCollisionTracks = aod::track::collisionId;
  Preslice<MyMuons> perCollisionMuons = aod::fwdtrack::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  bool fDoDetailedQA = false; // Bool to set detailed QA true, if QA is set true
  int fCurrentRun;            // needed to detect if the run changed and trigger update of calibrations etc.

  // TODO: filter on TPC dedx used temporarily until electron PID will be improved
  Filter barrelSelectedTracks = ifnode(fIsRun2.node() == true, aod::track::trackType == uint8_t(aod::track::Run2Track), aod::track::trackType == uint8_t(aod::track::Track)) && o2::aod::track::pt >= fConfigBarrelTrackPtLow && nabs(o2::aod::track::eta) <= fConfigBarrelTrackMaxAbsEta && o2::aod::track::tpcSignal >= configTpcSignal.fConfigMinTpcSignal && o2::aod::track::tpcSignal <= configTpcSignal.fConfigMaxTpcSignal && o2::aod::track::tpcChi2NCl < 4.0f && o2::aod::track::itsChi2NCl < 36.0f;

  Filter muonFilter = o2::aod::fwdtrack::pt >= fConfigMuonPtLow;

  void init(o2::framework::InitContext& context)
  {
    DefineCuts();
    fCCDB->setURL(useCCDBConfigurations.fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    if (useMatConfigurations.useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        fCCDB->get<TGeoManager>(geoPath);
      }
    }
    if (useMatConfigurations.useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(fCCDB->get<o2::base::MatLayerCylSet>(useMatConfigurations.lutPath));
      LOGF(info, "LUT load done!");
    }
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
    // Only use detailed QA when QA is set true
    if (fConfigQA && fConfigDetailedQA) {
      fDoDetailedQA = true;
    }

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";
    if (fDoDetailedQA) {
      histClasses += "Event_BeforeCuts;";
    }
    if (fConfigQA) {
      histClasses += "Event_AfterCuts;";
    }

    bool enableBarrelHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                               context.mOptions.get<bool>("processFullWithCent") || context.mOptions.get<bool>("processFullWithCovAndEventFilter") ||
                               context.mOptions.get<bool>("processFullWithCovMultsAndEventFilter") ||
                               context.mOptions.get<bool>("processBarrelOnly") || context.mOptions.get<bool>("processBarrelOnlyWithCent") || context.mOptions.get<bool>("processBarrelOnlyWithCovWithCent") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithMults") || context.mOptions.get<bool>("processBarrelOnlyWithCentAndMults") || context.mOptions.get<bool>("processBarrelOnlyWithCovWithCentAndMults") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithCov") || context.mOptions.get<bool>("processBarrelOnlyWithCovOnlyStdPID") || context.mOptions.get<bool>("processBarrelOnlyWithEventFilter") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithMultsAndEventFilter") || context.mOptions.get<bool>("processBarrelOnlyWithCovAndEventFilter") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithDalitzBits") || context.mOptions.get<bool>("processBarrelOnlyWithV0Bits") || context.mOptions.get<bool>("processBarrelWithDalitzEvent") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithV0BitsAndMaps") || context.mOptions.get<bool>("processAmbiguousBarrelOnly")) ||
                              context.mOptions.get<bool>("processBarrelWithV0AndDalitzEvent") || context.mOptions.get<bool>("processBarrelOnlyWithV0BitsAndMults");
    bool enableMuonHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                             context.mOptions.get<bool>("processFullWithCent") || context.mOptions.get<bool>("processFullWithCovAndEventFilter") ||
                             context.mOptions.get<bool>("processFullWithCovMultsAndEventFilter") || context.mOptions.get<bool>("processMuonOnlyWithCovAndCentMults") ||
                             context.mOptions.get<bool>("processMuonOnlyWithCov") || context.mOptions.get<bool>("processMuonOnlyWithCovAndEventFilter") || context.mOptions.get<bool>("processAmbiguousMuonOnlyWithCov") ||
                             context.mOptions.get<bool>("processMuonsAndMFT") || context.mOptions.get<bool>("processMuonsAndMFTWithFilter") || context.mOptions.get<bool>("processMuonMLOnly") || context.mOptions.get<bool>("processAssociatedMuonOnlyWithCov") || context.mOptions.get<bool>("processAssociatedRealignedMuonOnlyWithCov"));

    if (enableBarrelHistos) {
      if (fDoDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
        if (fIsAmbiguous) {
          histClasses += "Ambiguous_TrackBarrel_BeforeCuts;";
          histClasses += "Orphan_TrackBarrel;";
        }
      }
      if (fConfigQA) {
        for (auto& cut : fTrackCuts) {
          histClasses += Form("TrackBarrel_%s;", cut.GetName());
          if (fIsAmbiguous) {
            histClasses += Form("Ambiguous_TrackBarrel_%s;", cut.GetName());
          }
        }
      }
      if (fConfigIsOnlyforMaps) {
        histClasses += "TrackBarrel_PostCalibElectron;";
        histClasses += "TrackBarrel_PostCalibPion;";
        histClasses += "TrackBarrel_PostCalibProton;";
      }
    }
    if (enableMuonHistos) {
      if (fDoDetailedQA) {
        histClasses += "Muons_BeforeCuts;";
        histClasses += "MftTracks;";
        if (fIsAmbiguous) {
          histClasses += "Ambiguous_Muons_BeforeCuts;";
          histClasses += "Orphan_Muons_MFTMCHMID;";
          histClasses += "Orphan_Muons_MCHMID;";
        }
      }
      if (fConfigQA) {
        for (auto& muonCut : fMuonCuts) {
          histClasses += Form("Muons_%s;", muonCut.GetName());
          if (fIsAmbiguous) {
            histClasses += Form("Ambiguous_Muons_%s;", muonCut.GetName());
          }
        }
      }
    }

    DefineHistograms(histClasses);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
    // CCDB configuration
    if (fConfigComputeTPCpostCalib) {
      fCCDB->setURL(useCCDBConfigurations.fConfigCcdbUrl.value);
      fCCDB->setCaching(true);
      fCCDB->setLocalObjectValidityChecking();
      // Not later than now objects
      fCCDB->setCreatedNotAfter(fConfigNoLaterThan.value);
    }
  }

  void DefineCuts()
  {
    // Event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = configCuts.fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    // Barrel track cuts
    TString cutNamesStr = configCuts.fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // Muon cuts
    cutNamesStr = configCuts.fConfigMuonCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, uint32_t TMFTFillMap = 0u, typename TEvent, typename TTracks, typename TMuons, typename TAmbiTracks, typename TAmbiMuons, typename TMFTTracks = std::nullptr_t>
  void fullSkimming(TEvent const& collision, aod::BCsWithTimestamps const&, TTracks const& tracksBarrel, TMuons const& tracksMuon, TAmbiTracks const& ambiTracksMid, TAmbiMuons const& ambiTracksFwd, TMFTTracks const& mftTracks = nullptr)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (fCurrentRun != bc.runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(useCCDBConfigurations.fConfigCcdbPathTPC.value, bc.timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
        if (fConfigComputeTPCpostCalibKaon) {
          VarManager::SetCalibrationObject(VarManager::kTPCKaonMean, calibList->FindObject("mean_map_kaon"));
          VarManager::SetCalibrationObject(VarManager::kTPCKaonSigma, calibList->FindObject("sigma_map_kaon"));
        }
      }
      if (fIsRun2 == true) {
        grpmagrun2 = fCCDB->getForTimeStamp<o2::parameters::GRPObject>(grpmagPathRun2, bc.timestamp());
        if (grpmagrun2 != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmagrun2);
        }
      } else {
        grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
        if (grpmag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmag);
        }
        if (fPropMuon) {
          VarManager::SetupMuonMagField();
        }
        if (useMatConfigurations.useMatCorrType == 2) {
          // setMatLUT only after magfield has been initalized
          // (setMatLUT has implicit and problematic init field call if not)
          o2::base::Propagator::Instance()->setMatLUT(lut);
        }
      }
      fCurrentRun = bc.runNumber();
    }

    // store the selection decisions
    uint64_t tag = 0;
    // store some more information in the tag
    // if the BC found by event selection does not coincide with the collision.bc()
    auto bcEvSel = collision.template foundBC_as<aod::BCsWithTimestamps>();
    if (bcEvSel.globalIndex() != bc.globalIndex()) {
      tag |= (static_cast<uint64_t>(1) << 0);
    }
    // Put the 8 first bits of the event filter in the last 8 bits of the tag
    if constexpr ((TEventFillMap & VarManager::ObjTypes::EventFilter) > 0) {
      if (!useZorro.fConfigRunZorro) {
        tag |= (collision.eventFilter() << 56);
      }
    }

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    // TODO: These variables cannot be filled in the VarManager for the moment as long as BCsWithTimestamps are used.
    //       So temporarily, we filled them here, in order to be available for eventual QA of the skimming
    VarManager::fgValues[VarManager::kRunNo] = bc.runNumber();
    VarManager::fgValues[VarManager::kBC] = bc.globalBC();
    VarManager::fgValues[VarManager::kTimestamp] = bc.timestamp();
    VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
    if (fDoDetailedQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
    }

    uint32_t triggerAliases = collision.alias_raw();
    // fill stats information, before selections
    for (int i = 0; i < kNaliases; i++) {
      if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(kNaliases));

    if (useZorro.fConfigRunZorro) {
      zorro.setBaseCCDBPath(useCCDBConfigurations.fConfigCcdbPathZorro.value);
      zorro.setBCtolerance(useZorro.fBcTolerance);
      zorro.initCCDB(fCCDB.service, fCurrentRun, bc.timestamp(), useZorro.fConfigZorroTrigMask.value);
      zorro.populateExternalHists(fCurrentRun, reinterpret_cast<TH2D*>(fStatsList->At(3)), reinterpret_cast<TH2D*>(fStatsList->At(4)));
      bool zorroSel = zorro.isSelected(bc.globalBC(), useZorro.fBcTolerance, reinterpret_cast<TH2D*>(fStatsList->At(4)));
      if (zorroSel) {
        tag |= (static_cast<uint64_t>(true) << 56); // the same bit is used for this zorro selections from ccdb
      }
      if (useZorro.fConfigRunZorroSel && (!zorroSel || !fEventCut->IsSelected(VarManager::fgValues))) {
        return;
      }
    } else {
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        return;
      }
    }

    // fill stats information, after selections
    for (int i = 0; i < kNaliases; i++) {
      if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(kNaliases));

    fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

    // create the event tables
    event(tag, bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
    if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0 && (TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                    collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                    collision.centFT0C(), collision.centFT0A(), collision.centFT0M());
    } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0) {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                    collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                    -1, -1, -1);
    } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, collision.centFT0C(), collision.centFT0A(), collision.centFT0M());
    } else {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    }
    eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
    eventInfo(collision.globalIndex());
    if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMultExtra) > 0) {
      multPV(collision.multNTracksHasITS(), collision.multNTracksHasTPC(), collision.multNTracksHasTOF(), collision.multNTracksHasTRD(),
             collision.multNTracksITSOnly(), collision.multNTracksTPCOnly(), collision.multNTracksITSTPC(),
             collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
      multAll(collision.multAllTracksTPCOnly(), collision.multAllTracksITSTPC(),
              0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0);
    }

    uint64_t trackFilteringTag = 0;
    uint8_t fwdFilteringTag = 0;
    uint8_t trackTempFilterMap = 0;
    int isAmbiguous = 0;
    if constexpr (static_cast<bool>(TTrackFillMap)) {
      trackBarrelInfo.reserve(tracksBarrel.size());
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
        trackBarrelCov.reserve(tracksBarrel.size());
      }
      trackBarrelPID.reserve(tracksBarrel.size());

      // loop over tracks
      for (auto& track : tracksBarrel) {
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::AmbiTrack) > 0) {
          if (fIsAmbiguous) {
            isAmbiguous = 0;
            for (auto& ambiTrackMid : ambiTracksMid) {
              if (ambiTrackMid.trackId() == track.globalIndex()) {
                isAmbiguous = 1;
                break;
              }
            }
          }
        }

        trackFilteringTag = static_cast<uint64_t>(0);
        trackTempFilterMap = uint8_t(0);
        VarManager::FillTrack<TTrackFillMap>(track);
        if (fDoDetailedQA) {
          fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
          if (fIsAmbiguous && isAmbiguous == 1) {
            fHistMan->FillHistClass("Ambiguous_TrackBarrel_BeforeCuts", VarManager::fgValues);
          }
        }

        // apply track cuts and fill stats histogram
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            trackTempFilterMap |= (uint8_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
              if (fIsAmbiguous && isAmbiguous == 1) {
                fHistMan->FillHistClass(Form("Ambiguous_TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
              }
            }
            (reinterpret_cast<TH1D*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
          }
        }
        if (!trackTempFilterMap) {
          continue;
        }

        // store filtering information
        /*if (track.isGlobalTrack()) {
          trackFilteringTag |= (uint64_t(1) << 0); // BIT0: global track
        }
        if (track.isGlobalTrackSDD()) {
          trackFilteringTag |= (uint64_t(1) << 1); // BIT1: global track SSD
        }*/
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackV0Bits)) { // BIT0-4: V0Bits
          trackFilteringTag = static_cast<uint64_t>(track.pidbit());
          for (int iv0 = 0; iv0 < 5; iv0++) {
            if (track.pidbit() & (uint8_t(1) << iv0)) {
              (reinterpret_cast<TH1D*>(fStatsList->At(1)))->Fill(fTrackCuts.size() + static_cast<float>(iv0));
            }
          }
          if (fConfigIsOnlyforMaps) {
            if (trackFilteringTag & (static_cast<uint64_t>(1) << VarManager::kIsConversionLeg)) { // for electron
              fHistMan->FillHistClass("TrackBarrel_PostCalibElectron", VarManager::fgValues);
            }
            if (trackFilteringTag & (static_cast<uint64_t>(1) << VarManager::kIsK0sLeg)) { // for pion
              fHistMan->FillHistClass("TrackBarrel_PostCalibPion", VarManager::fgValues);
            }
            if ((static_cast<bool>(trackFilteringTag & (static_cast<uint64_t>(1) << VarManager::kIsLambdaLeg)) * (track.sign()) > 0)) { // for proton from Lambda
              fHistMan->FillHistClass("TrackBarrel_PostCalibProton", VarManager::fgValues);
            }
            if ((static_cast<bool>(trackFilteringTag & (static_cast<uint64_t>(1) << VarManager::kIsALambdaLeg)) * (track.sign()) < 0)) { // for proton from AntiLambda
              fHistMan->FillHistClass("TrackBarrel_PostCalibProton", VarManager::fgValues);
            }
          }
        }
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
          trackFilteringTag |= (static_cast<uint64_t>(track.dalitzBits()) << VarManager::kDalitzBits); // BIT5-12: Dalitz selection bits
        }
        trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << VarManager::kBarrelUserCutsBits); // BIT13-20...:  user track filters

        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPID)) {
          if (fConfigComputeTPCpostCalib) {
            trackFilteringTag |= (static_cast<uint64_t>(1) << VarManager::kIsTPCPostcalibrated); // store the info on whether TPC pid is skimmed as postcalibrated
          }
        }

        // create the track tables
        trackBarrelInfo(track.collisionId(), collision.posX(), collision.posY(), collision.posZ(), track.globalIndex());
        trackBasic(event.lastIndex(), trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), isAmbiguous);
        trackBarrel(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                    track.tpcInnerParam(), track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                    track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.tpcChi2NCl(),
                    track.trdChi2(), track.trdPattern(), track.tofChi2(),
                    track.length(), track.dcaXY(), track.dcaZ(),
                    track.trackTime(), track.trackTimeRes(), track.tofExpMom(),
                    track.detectorMap());
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
          trackBarrelCov(track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                         track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                         track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
        }
        // NOTE: If the TPC postcalibration is switched on, then we write the postcalibrated n-sigma values directly in the skimmed data
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPID)) {
          float nSigmaEl = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaEl_Corr] : track.tpcNSigmaEl());
          float nSigmaPi = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaPi_Corr] : track.tpcNSigmaPi());
          float nSigmaKa = ((fConfigComputeTPCpostCalib && fConfigComputeTPCpostCalibKaon) ? VarManager::fgValues[VarManager::kTPCnSigmaKa_Corr] : track.tpcNSigmaKa());
          float nSigmaPr = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaPr_Corr] : track.tpcNSigmaPr());
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPIDExtra)) {
            trackBarrelPID(track.tpcSignal(),
                           nSigmaEl, track.tpcNSigmaMu(), nSigmaPi, nSigmaKa, nSigmaPr,
                           track.beta(),
                           track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                           track.trdSignal());
          } else {
            trackBarrelPID(track.tpcSignal(),
                           nSigmaEl, -1, nSigmaPi, nSigmaKa, nSigmaPr,
                           -1,
                           track.tofNSigmaEl(), -1, track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                           -1);
          }
        }
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackTPCPID)) {
          trackBarrelPID(track.tpcSignal(), track.tpcNSigmaEl(), -1, track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         -1, -1, -1, -1, -1, -1, -1);
        }
      }
    } // end if constexpr (TTrackFillMap)

    // Maps for the MFT-muon matching index
    std::map<int, int> newMFTTableSize; // key : oldMFTIndex, value: size of the table-1 at step key
    std::map<uint64_t, int> mftOffsets; // key: mftoldglobalindex, value: mft.offsets

    if constexpr (static_cast<bool>(TMFTFillMap)) {
      trackMFT.reserve(mftTracks.size());
      trackMFTExtra.reserve(mftTracks.size());
      // TODO add cuts on the MFT tracks
      // int nDel = 0;
      for (auto& mft : mftTracks) {
        if (false) // for now no cuts
        {
          // nDel++;
        } else { // it passes the cuts and will be saved in the tables
          newMFTTableSize[mft.index()] = trackMFT.lastIndex();
        }

        mftOffsets[mft.globalIndex()] = mft.offsets();

        double chi2 = mft.chi2();
        SMatrix5 tpars(mft.x(), mft.y(), mft.phi(), mft.tgl(), mft.signed1Pt());
        std::vector<double> v1;
        SMatrix55 tcovs(v1.begin(), v1.end());
        o2::track::TrackParCovFwd pars1{mft.z(), tpars, tcovs, chi2};
        pars1.propagateToZlinear(collision.posZ());

        double dcaX = (pars1.getX() - collision.posX());
        double dcaY = (pars1.getY() - collision.posY());

        VarManager::FillTrack<gkMFTFillMap>(mft);
        fHistMan->FillHistClass("MftTracks", VarManager::fgValues);

        trackMFT(event.lastIndex(), fwdFilteringTag, mft.pt(), mft.eta(), mft.phi());
        trackMFTExtra(mft.mftClusterSizesAndTrackFlags(), mft.sign(), dcaX, dcaY, mft.nClusters());
      } // end of mft : mftTracks

    } // end if constexpr (TMFTFillMap)

    if constexpr (static_cast<bool>(TMuonFillMap)) {
      // build the muon tables
      muonBasic.reserve(tracksMuon.size());
      muonExtra.reserve(tracksMuon.size());
      muonInfo.reserve(tracksMuon.size());
      if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
        muonCov.reserve(tracksMuon.size());
      }
      // loop over muons

      // first we need to get the correct indices
      int nDel = 0;
      int idxPrev = -1;
      std::map<int, int> newEntryNb;
      std::map<int, int> newMatchIndex;
      std::map<int, int> newMFTMatchIndex;

      for (auto& muon : tracksMuon) {
        fwdFilteringTag = static_cast<uint64_t>(0);
        VarManager::FillTrack<TMuonFillMap>(muon);
        if (fPropMuon) {
          VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
        }

        if (muon.index() > idxPrev + 1) { // checks if some muons are filtered even before the skimming function
          nDel += muon.index() - (idxPrev + 1);
        }
        idxPrev = muon.index();

        // check the cuts and filters
        int i = 0;
        for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues))
            trackTempFilterMap |= (uint8_t(1) << i);
        }

        if (!trackTempFilterMap) { // does not pass the cuts
          nDel++;
        } else { // it passes the cuts and will be saved in the tables
          newEntryNb[muon.index()] = muon.index() - nDel;
        }
      }

      // now let's save the muons with the correct indices and matches
      for (auto& muon : tracksMuon) {
        if constexpr ((TMuonFillMap & VarManager::ObjTypes::AmbiMuon) > 0) {
          if (fIsAmbiguous) {
            isAmbiguous = 0;
            for (auto& ambiTrackFwd : ambiTracksFwd) {
              if (ambiTrackFwd.fwdtrackId() == muon.globalIndex()) {
                isAmbiguous = 1;
                break;
              }
            }
          }
        }
        fwdFilteringTag = uint8_t(0);
        trackTempFilterMap = uint8_t(0);

        VarManager::FillTrack<TMuonFillMap>(muon);

        // recalculte pDca for global muon tracks
        VarManager::FillTrackCollision<TMuonFillMap>(muon, collision);

        if (fPropMuon) {
          VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
        }
        if (fDoDetailedQA) {
          fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
          if (fIsAmbiguous && isAmbiguous == 1) {
            fHistMan->FillHistClass("Ambiguous_Muons_BeforeCuts", VarManager::fgValues);
          }
        }
        // apply the muon selection cuts and fill the stats histogram
        int i = 0;
        for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            trackTempFilterMap |= (uint8_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(Form("Muons_%s", (*cut).GetName()), VarManager::fgValues);
              if (fIsAmbiguous && isAmbiguous == 1) {
                fHistMan->FillHistClass(Form("Ambiguous_Muons_%s", (*cut).GetName()), VarManager::fgValues);
              }
            }
            (reinterpret_cast<TH1D*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
          }
        }
        if (!trackTempFilterMap) {
          continue;
        }
        // store the cut decisions
        trackFilteringTag = trackTempFilterMap; // BIT0-7:  user selection cuts
        if (fPropMuon) {
          trackFilteringTag |= (uint8_t(1) << VarManager::kMuonIsPropagated); // store the info on whether the muon is propagated or not
        }

        // update the matching MCH/MFT index
        if (static_cast<int>(muon.trackType()) == 0 || static_cast<int>(muon.trackType()) == 2) { // MCH-MFT(2) or GLB(0) track
          int matchIdx = muon.matchMCHTrackId() - muon.offsets();                                 // simple match index, not the global index
          int matchMFTIdx = muon.matchMFTTrackId() - mftOffsets[muon.matchMFTTrackId()];

          // first for MCH matching index
          if (newEntryNb.count(matchIdx) > 0) {                                                  // if the key exists i.e the match will not get deleted
            newMatchIndex[muon.index()] = newEntryNb[matchIdx];                                  // update the match for this muon to the updated entry of the match
            newMatchIndex[muon.index()] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()]; // adding the offset of muons, muonBasic.lastIndex() start at -1

            if (static_cast<int>(muon.trackType()) == 0) {                                     // for now only do this to global tracks
              newMatchIndex[matchIdx] = newEntryNb[muon.index()];                              // add the  updated index of this muon as a match to mch track
              newMatchIndex[matchIdx] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()]; // adding the offset, muonBasic.lastIndex() start at -1
            }
          } else {
            newMatchIndex[muon.index()] = -1;
          }

          // then for MFT match index
          if (newMFTTableSize.count(matchMFTIdx) > 0) {                        // if the key exists i.e the match will not get deleted
            newMFTMatchIndex[muon.index()] = newMFTTableSize[matchMFTIdx] + 1; // adding the offset of mfts, newMFTTableSize start at -1
          } else {
            newMFTMatchIndex[muon.index()] = -1;
          }

        } else if (static_cast<int>(muon.trackType() == 4)) { // an MCH track

          newMFTMatchIndex[muon.index()] = -1;
          // in this case the matches should be filled from the other types but we need to check
          if (newMatchIndex.count(muon.index()) == 0) { // if an entry for this mch was not added it simply mean that non of the global tracks were matched to it
            newMatchIndex[muon.index()] = -1;
          }
        }

        muonBasic(event.lastIndex(), newMatchIndex.find(muon.index())->second, newMFTMatchIndex.find(muon.index())->second, trackFilteringTag, VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], muon.sign(), isAmbiguous);
        muonInfo(muon.collisionId(), collision.posX(), collision.posY(), collision.posZ());
        if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {

          if (fPropMuon) {
            muonExtra(muon.nClusters(), VarManager::fgValues[VarManager::kMuonPDca], VarManager::fgValues[VarManager::kMuonRAtAbsorberEnd],
                      muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                      muon.matchScoreMCHMFT(), muon.mchBitMap(), muon.midBitMap(),
                      muon.midBoards(), muon.trackType(), VarManager::fgValues[VarManager::kMuonDCAx], VarManager::fgValues[VarManager::kMuonDCAy],
                      muon.trackTime(), muon.trackTimeRes());
          } else {
            muonExtra(muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
                      muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                      muon.matchScoreMCHMFT(), muon.mchBitMap(), muon.midBitMap(),
                      muon.midBoards(), muon.trackType(), muon.fwdDcaX(), muon.fwdDcaY(),
                      muon.trackTime(), muon.trackTimeRes());
          }

          muonCov(VarManager::fgValues[VarManager::kX], VarManager::fgValues[VarManager::kY], VarManager::fgValues[VarManager::kZ], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kTgl], muon.sign() / VarManager::fgValues[VarManager::kPt],
                  VarManager::fgValues[VarManager::kMuonCXX], VarManager::fgValues[VarManager::kMuonCXY], VarManager::fgValues[VarManager::kMuonCYY], VarManager::fgValues[VarManager::kMuonCPhiX], VarManager::fgValues[VarManager::kMuonCPhiY], VarManager::fgValues[VarManager::kMuonCPhiPhi],
                  VarManager::fgValues[VarManager::kMuonCTglX], VarManager::fgValues[VarManager::kMuonCTglY], VarManager::fgValues[VarManager::kMuonCTglPhi], VarManager::fgValues[VarManager::kMuonCTglTgl], VarManager::fgValues[VarManager::kMuonC1Pt2X], VarManager::fgValues[VarManager::kMuonC1Pt2Y],
                  VarManager::fgValues[VarManager::kMuonC1Pt2Phi], VarManager::fgValues[VarManager::kMuonC1Pt2Tgl], VarManager::fgValues[VarManager::kMuonC1Pt21Pt2]);
        } else {
          muonExtra(muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
                    muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                    muon.matchScoreMCHMFT(), muon.mchBitMap(), muon.midBitMap(),
                    muon.midBoards(), muon.trackType(), muon.fwdDcaX(), muon.fwdDcaY(),
                    muon.trackTime(), muon.trackTimeRes());
        }
      }
    } // end if constexpr (TMuonFillMap)
  } // end fullSkimming()

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons, typename AssocTracks, typename AssocMuons>
  void fullSkimmingIndices(TEvent const& collision, aod::BCsWithTimestamps const&, TTracks const& tracksBarrel, TMuons const& tracksMuon, AssocTracks const& trackIndices, AssocMuons const& fwdtrackIndices)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (fCurrentRun != bc.runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(useCCDBConfigurations.fConfigCcdbPathTPC.value, bc.timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
        if (fConfigComputeTPCpostCalibKaon) {
          VarManager::SetCalibrationObject(VarManager::kTPCKaonMean, calibList->FindObject("mean_map_kaon"));
          VarManager::SetCalibrationObject(VarManager::kTPCKaonSigma, calibList->FindObject("sigma_map_kaon"));
        }
      }
      if (fIsRun2 == true) {
        grpmagrun2 = fCCDB->getForTimeStamp<o2::parameters::GRPObject>(grpmagPathRun2, bc.timestamp());
        if (grpmagrun2 != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmagrun2);
        }
      } else {
        grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
        if (grpmag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmag);
        }
        if constexpr (static_cast<bool>(TMuonFillMap)) {
          VarManager::SetupMuonMagField();
        }
        if (useMatConfigurations.useMatCorrType == 2) {
          // setMatLUT only after magfield has been initalized
          // (setMatLUT has implicit and problematic init field call if not)
          o2::base::Propagator::Instance()->setMatLUT(lut);
        }
      }
      fCurrentRun = bc.runNumber();
    }

    // get the trigger aliases
    uint32_t triggerAliases = collision.alias_raw();
    // store the selection decisions
    uint64_t tag = 0;
    // store some more information in the tag
    // if the BC found by event selection does not coincide with the collision.bc()
    auto bcEvSel = collision.template foundBC_as<aod::BCsWithTimestamps>();
    if (bcEvSel.globalIndex() != bc.globalIndex()) {
      tag |= (static_cast<uint64_t>(1) << 0);
    }
    // Put the 8 first bits of the event filter in the last 8 bits of the tag
    if constexpr ((TEventFillMap & VarManager::ObjTypes::EventFilter) > 0) {
      if (!useZorro.fConfigRunZorro) {
        tag |= (collision.eventFilter() << 56);
      }
    }

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    // TODO: These variables cannot be filled in the VarManager for the moment as long as BCsWithTimestamps are used.
    //       So temporarily, we filled them here, in order to be available for eventual QA of the skimming
    VarManager::fgValues[VarManager::kRunNo] = bc.runNumber();
    VarManager::fgValues[VarManager::kBC] = bc.globalBC();
    VarManager::fgValues[VarManager::kTimestamp] = bc.timestamp();
    VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
    if (fDoDetailedQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
    }

    // fill stats information, before selections
    for (int i = 0; i < kNaliases; i++) {
      if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(kNaliases));

    if (useZorro.fConfigRunZorro) {
      zorro.setBaseCCDBPath(useCCDBConfigurations.fConfigCcdbPathZorro.value);
      zorro.setBCtolerance(useZorro.fBcTolerance);
      zorro.initCCDB(fCCDB.service, fCurrentRun, bc.timestamp(), useZorro.fConfigZorroTrigMask.value);
      zorro.populateExternalHists(fCurrentRun, reinterpret_cast<TH2D*>(fStatsList->At(3)), reinterpret_cast<TH2D*>(fStatsList->At(4)));
      bool zorroSel = zorro.isSelected(bc.globalBC(), useZorro.fBcTolerance, reinterpret_cast<TH2D*>(fStatsList->At(4)));
      if (zorroSel) {
        tag |= (static_cast<uint64_t>(true) << 56); // the same bit is used for this zorro selections from ccdb
      }
      if (useZorro.fConfigRunZorroSel && (!zorroSel || !fEventCut->IsSelected(VarManager::fgValues))) {
        return;
      }
    } else {
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        return;
      }
    }

    // fill stats information, after selections
    for (int i = 0; i < kNaliases; i++) {
      if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(kNaliases));

    fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

    // create the event tables
    event(tag, bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
    if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0 && (TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                    collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                    collision.centFT0C(), collision.centFT0A(), collision.centFT0M());
    } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0) {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                    collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                    -1, -1, -1);
    } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, collision.centFT0C(), collision.centFT0A(), collision.centFT0M());
    } else {
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    }
    eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
    eventInfo(collision.globalIndex());
    if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMultExtra) > 0) {
      multPV(collision.multNTracksHasITS(), collision.multNTracksHasTPC(), collision.multNTracksHasTOF(), collision.multNTracksHasTRD(),
             collision.multNTracksITSOnly(), collision.multNTracksTPCOnly(), collision.multNTracksITSTPC(),
             collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
      multAll(collision.multAllTracksTPCOnly(), collision.multAllTracksITSTPC(),
              0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0);
    }

    uint64_t trackFilteringTag = 0;
    uint8_t trackTempFilterMap = 0;
    int isAmbiguous = 0;
    if constexpr (static_cast<bool>(TTrackFillMap)) {
      trackBarrelInfo.reserve(tracksBarrel.size());
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
        trackBarrelCov.reserve(tracksBarrel.size());
      }
      trackBarrelPID.reserve(tracksBarrel.size());

      // loop over tracks
      for (const auto& trackId : trackIndices) { // start loop over tracks
        auto track = trackId.template track_as<TTracks>();
        if constexpr ((TTrackFillMap & VarManager::ObjTypes::AmbiTrack) > 0) {
          if (fIsAmbiguous) {
            isAmbiguous = (track.compatibleCollIds().size() != 1);
          }
        }
        trackFilteringTag = static_cast<uint64_t>(0);
        trackTempFilterMap = uint8_t(0);
        VarManager::FillTrack<TTrackFillMap>(track);
        if (fDoDetailedQA) {
          fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
          if (fIsAmbiguous && isAmbiguous == 1) {
            fHistMan->FillHistClass("Ambiguous_TrackBarrel_BeforeCuts", VarManager::fgValues);
          }
        }

        // apply track cuts and fill stats histogram
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            trackTempFilterMap |= (uint8_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
              if (fIsAmbiguous && isAmbiguous == 1) {
                fHistMan->FillHistClass(Form("Ambiguous_TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
              }
            }
            (reinterpret_cast<TH1D*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
          }
        }
        if (!trackTempFilterMap) {
          continue;
        }

        // store filtering information
        if (track.isGlobalTrack()) {
          trackFilteringTag |= (static_cast<uint64_t>(1) << 0); // BIT0: global track
        }
        if (track.isGlobalTrackSDD()) {
          trackFilteringTag |= (static_cast<uint64_t>(1) << 1); // BIT1: global track SSD
        }
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackV0Bits)) { // BIT2-6: V0Bits
          trackFilteringTag |= (static_cast<uint64_t>(track.pidbit()) << 2);
          for (int iv0 = 0; iv0 < 5; iv0++) {
            if (track.pidbit() & (uint8_t(1) << iv0)) {
              (reinterpret_cast<TH1D*>(fStatsList->At(1)))->Fill(fTrackCuts.size() + static_cast<float>(iv0));
            }
          }
          if (fConfigIsOnlyforMaps) {
            if (trackFilteringTag & (static_cast<uint64_t>(1) << 2)) { // for electron
              fHistMan->FillHistClass("TrackBarrel_PostCalibElectron", VarManager::fgValues);
            }
            if (trackFilteringTag & (static_cast<uint64_t>(1) << 3)) { // for pion
              fHistMan->FillHistClass("TrackBarrel_PostCalibPion", VarManager::fgValues);
            }
            if ((static_cast<bool>(trackFilteringTag & (static_cast<uint64_t>(1) << 4)) * (track.sign()) > 0)) { // for proton from Lambda
              fHistMan->FillHistClass("TrackBarrel_PostCalibProton", VarManager::fgValues);
            }
            if ((static_cast<bool>(trackFilteringTag & (static_cast<uint64_t>(1) << 5)) * (track.sign()) < 0)) { // for proton from AntiLambda
              fHistMan->FillHistClass("TrackBarrel_PostCalibProton", VarManager::fgValues);
            }
          }
        }
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
          trackFilteringTag |= (static_cast<uint64_t>(track.dalitzBits()) << 7); // BIT7-14: Dalitz
        }
        trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << 15); // BIT15-...:  user track filters

        // create the track tables
        trackBarrelInfo(track.collisionId(), collision.posX(), collision.posY(), collision.posZ(), track.globalIndex());
        trackBasic(event.lastIndex(), trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), isAmbiguous);
        trackBarrel(track.tpcInnerParam(), track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                    track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.tpcChi2NCl(),
                    track.trdChi2(), track.trdPattern(), track.tofChi2(),
                    track.length(), track.dcaXY(), track.dcaZ(),
                    track.trackTime(), track.trackTimeRes(), track.tofExpMom(),
                    track.detectorMap());
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
          trackBarrelCov(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                         track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                         track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                         track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
        }
        // NOTE: If the TPC postcalibration is switched on, then we write the postcalibrated n-sigma values directly in the skimmed data
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPID)) {
          float nSigmaEl = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaEl_Corr] : track.tpcNSigmaEl());
          float nSigmaPi = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaPi_Corr] : track.tpcNSigmaPi());
          float nSigmaPr = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaPr_Corr] : track.tpcNSigmaPr());
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPIDExtra)) {
            trackBarrelPID(track.tpcSignal(),
                           nSigmaEl, track.tpcNSigmaMu(), nSigmaPi, track.tpcNSigmaKa(), nSigmaPr,
                           track.beta(),
                           track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                           track.trdSignal());
          } else {
            trackBarrelPID(track.tpcSignal(),
                           nSigmaEl, -1, nSigmaPi, track.tpcNSigmaKa(), nSigmaPr,
                           -1,
                           track.tofNSigmaEl(), -1, track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                           -1);
          }
        }
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackTPCPID)) {
          trackBarrelPID(track.tpcSignal(), track.tpcNSigmaEl(), -1, track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         -1, -1, -1, -1, -1, -1, -1);
        }
      }
    } // end if constexpr (TTrackFillMap)

    if constexpr (static_cast<bool>(TMuonFillMap)) {
      // build the muon tables
      muonBasic.reserve(tracksMuon.size());
      muonExtra.reserve(tracksMuon.size());
      muonInfo.reserve(tracksMuon.size());
      if constexpr (static_cast<bool>((TMuonFillMap & VarManager::ObjTypes::MuonCov) || (TMuonFillMap & VarManager::ObjTypes::MuonCovRealign))) {
        muonCov.reserve(tracksMuon.size());
      }
      // loop over muons

      // first we need to get the correct indices
      int nDel = 0;
      int idxPrev = -1;
      std::map<int, int> newEntryNb;
      std::map<int, int> newMatchIndex;

      for (const auto& muonId : fwdtrackIndices) { // start loop over tracks
        auto muon = tracksMuon.rawIteratorAt(muonId.fwdtrackId());
        trackFilteringTag = static_cast<uint64_t>(0);

        VarManager::FillTrack<TMuonFillMap>(muon);

        if (muon.index() > idxPrev + 1) { // checks if some muons are filtered even before the skimming function
          nDel += muon.index() - (idxPrev + 1);
        }
        idxPrev = muon.index();

        // check the cuts and filters
        int i = 0;
        for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues))
            trackTempFilterMap |= (uint8_t(1) << i);
        }

        if (!trackTempFilterMap) { // does not pass the cuts
          nDel++;
        } else { // it passes the cuts and will be saved in the tables
          newEntryNb[muon.index()] = muon.index() - nDel;
        }
      }

      // now let's save the muons with the correct indices and matches
      for (const auto& muonId : fwdtrackIndices) { // start loop over tracks
        auto muon = tracksMuon.rawIteratorAt(muonId.fwdtrackId());
        if constexpr ((TMuonFillMap & VarManager::ObjTypes::AmbiMuon) > 0) {
          if (fIsAmbiguous) {
            isAmbiguous = (muon.compatibleCollIds().size() != 1);
          }
        }
        trackFilteringTag = static_cast<uint64_t>(0);
        trackTempFilterMap = uint8_t(0);

        if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonRealign)) {
          if (static_cast<bool>(muon.isRemovable())) {
            continue;
          }
        }

        VarManager::FillTrack<TMuonFillMap>(muon);

        // recalculte pDca for global muon tracks
        VarManager::FillTrackCollision<TMuonFillMap>(muon, collision);

        if (fPropMuon) {
          VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
        }

        if (fDoDetailedQA) {
          fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
          if (fIsAmbiguous && isAmbiguous == 1) {
            fHistMan->FillHistClass("Ambiguous_Muons_BeforeCuts", VarManager::fgValues);
          }
        }
        // apply the muon selection cuts and fill the stats histogram
        int i = 0;
        for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            trackTempFilterMap |= (uint8_t(1) << i);
            if (fConfigQA) {
              fHistMan->FillHistClass(Form("Muons_%s", (*cut).GetName()), VarManager::fgValues);
              if (fIsAmbiguous && isAmbiguous == 1) {
                fHistMan->FillHistClass(Form("Ambiguous_Muons_%s", (*cut).GetName()), VarManager::fgValues);
              }
            }
            (reinterpret_cast<TH1D*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
          }
        }
        if (!trackTempFilterMap) {
          continue;
        }

        // store the cut decisions
        trackFilteringTag = trackTempFilterMap; // BIT0-7:  user selection cuts
        if (fPropMuon) {
          trackFilteringTag |= (uint8_t(1) << VarManager::kMuonIsPropagated); // store the info on whether the muon is propagated or not
        }

        // update the matching MCH/MFT index
        if (static_cast<int>(muon.trackType()) == 0 || static_cast<int>(muon.trackType()) == 2) { // MCH-MFT(2) or GLB(0) track
          int matchIdx = muon.matchMCHTrackId() - muon.offsets();
          if (newEntryNb.count(matchIdx) > 0) {                                                  // if the key exists i.e the match will not get deleted
            newMatchIndex[muon.index()] = newEntryNb[matchIdx];                                  // update the match for this muon to the updated entry of the match
            newMatchIndex[muon.index()] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()]; // adding the offset of muons, muonBasic.lastIndex() start at -1
            if (static_cast<int>(muon.trackType()) == 0) {                                       // for now only do this to global tracks
              newMatchIndex[matchIdx] = newEntryNb[muon.index()];                                // add the  updated index of this muon as a match to mch track
              newMatchIndex[matchIdx] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()];   // adding the offset, muonBasic.lastIndex() start at -1
            }
          } else {
            newMatchIndex[muon.index()] = -1;
          }
        } else if (static_cast<int>(muon.trackType() == 4)) { // an MCH track
          // in this case the matches should be filled from the other types but we need to check
          if (newMatchIndex.count(muon.index()) == 0) { // if an entry for this mch was not added it simply mean that non of the global tracks were matched to it
            newMatchIndex[muon.index()] = -1;
          }
        }

        muonBasic(event.lastIndex(), newMatchIndex.find(muon.index())->second, -1, trackFilteringTag, VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], muon.sign(), isAmbiguous);
        muonInfo(muon.collisionId(), collision.posX(), collision.posY(), collision.posZ());
        if constexpr (static_cast<bool>((TMuonFillMap & VarManager::ObjTypes::MuonCov) || (TMuonFillMap & VarManager::ObjTypes::MuonCovRealign))) {

          if (fPropMuon) {
            muonExtra(muon.nClusters(), VarManager::fgValues[VarManager::kMuonPDca], VarManager::fgValues[VarManager::kMuonRAtAbsorberEnd],
                      VarManager::fgValues[VarManager::kMuonChi2], muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                      muon.matchScoreMCHMFT(), muon.mchBitMap(), muon.midBitMap(),
                      muon.midBoards(), muon.trackType(), VarManager::fgValues[VarManager::kMuonDCAx], VarManager::fgValues[VarManager::kMuonDCAy],
                      muon.trackTime(), muon.trackTimeRes());
          } else {
            muonExtra(muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
                      VarManager::fgValues[VarManager::kMuonChi2], muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                      muon.matchScoreMCHMFT(), muon.mchBitMap(), muon.midBitMap(),
                      muon.midBoards(), muon.trackType(), muon.fwdDcaX(), muon.fwdDcaY(),
                      muon.trackTime(), muon.trackTimeRes());
          }

          muonCov(VarManager::fgValues[VarManager::kX], VarManager::fgValues[VarManager::kY], VarManager::fgValues[VarManager::kZ], VarManager::fgValues[VarManager::kPhi], VarManager::fgValues[VarManager::kTgl], muon.sign() / VarManager::fgValues[VarManager::kPt],
                  VarManager::fgValues[VarManager::kMuonCXX], VarManager::fgValues[VarManager::kMuonCXY], VarManager::fgValues[VarManager::kMuonCYY], VarManager::fgValues[VarManager::kMuonCPhiX], VarManager::fgValues[VarManager::kMuonCPhiY], VarManager::fgValues[VarManager::kMuonCPhiPhi],
                  VarManager::fgValues[VarManager::kMuonCTglX], VarManager::fgValues[VarManager::kMuonCTglY], VarManager::fgValues[VarManager::kMuonCTglPhi], VarManager::fgValues[VarManager::kMuonCTglTgl], VarManager::fgValues[VarManager::kMuonC1Pt2X], VarManager::fgValues[VarManager::kMuonC1Pt2Y],
                  VarManager::fgValues[VarManager::kMuonC1Pt2Phi], VarManager::fgValues[VarManager::kMuonC1Pt2Tgl], VarManager::fgValues[VarManager::kMuonC1Pt21Pt2]);
        }
      }
    } // end if constexpr (TMuonFillMap)
  } // end fullSkimming()

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      if (fConfigQA) {
        fHistMan->AddHistClass(classStr.Data());
      }

      // fill the THn histograms
      if (fConfigIsOnlyforMaps) {
        if (classStr.Contains("PostCalibElectron")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PostCalibPion")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PostCalibProton")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
      }

      TString histEventName = addHistoConfigurations.fConfigAddEventHistogram.value;
      if (classStr.Contains("Event")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", histEventName);
        }
      }

      TString histTrackName = addHistoConfigurations.fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
        }
      }

      TString histMuonName = addHistoConfigurations.fConfigAddMuonHistogram.value;
      if (classStr.Contains("Muons")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histMuonName);
        }
      }

      TString histMftName = addHistoConfigurations.fConfigAddMuonHistogram.value;
      if (classStr.Contains("Mft")) {
        if (fConfigDetailedQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histMftName);
        }
      }
    }

    // create statistics histograms
    // 0: Event statistics
    // 1: Track statistics
    // 2: Muon statistics
    // 3: Zorro information
    // 4: Zorro trigger selection
    // NOTE: Please keep the order of the histograms in the list
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH2F* histEvents = new TH2F("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, kNaliases + 1, -0.5, +kNaliases + 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    for (int ib = 1; ib <= kNaliases; ib++) {
      histEvents->GetYaxis()->SetBinLabel(ib, aliasLabels[ib - 1].data());
    }
    histEvents->GetYaxis()->SetBinLabel(kNaliases + 1, "Total");
    fStatsList->Add(histEvents); // At index 0

    // Track statistics: one bin for each track selection and 5 bins for V0 tags (gamma, K0s, Lambda, anti-Lambda, Omega)
    TH1D* histTracks = new TH1D("TrackStats", "Track statistics", fTrackCuts.size() + 5.0, -0.5, fTrackCuts.size() - 0.5 + 5.0);
    ib = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ib++) {
      histTracks->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    const char* v0TagNames[5] = {"Photon conversion", "K^{0}_{s}", "#Lambda", "#bar{#Lambda}", "#Omega"};
    for (ib = 0; ib < 5; ib++) {
      histTracks->GetXaxis()->SetBinLabel(fTrackCuts.size() + 1 + ib, v0TagNames[ib]);
    }
    fStatsList->Add(histTracks); // At index 1
    TH1D* histMuons = new TH1D("MuonStats", "Muon statistics", fMuonCuts.size(), -0.5, fMuonCuts.size() - 0.5);
    ib = 1;
    for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, ib++) {
      histMuons->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    fStatsList->Add(histMuons); // At index 2

    TH2D* histZorroInfo = new TH2D("ZorroInfo", "Zorro information", 1, -0.5, 0.5, 1, -0.5, 0.5);
    fStatsList->Add(histZorroInfo); // At index 3

    TH2D* histZorroSel = new TH2D("ZorroSel", "trigger of interested", 1, -0.5, 0.5, 1, -0.5, 0.5);
    fStatsList->Add(histZorroSel); // At index 4
  }

  // Produce barrel + muon tables -------------------------------------------------------------------------------------------------------------
  void processFull(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                   soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
  }

  void processFullWithCov(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                          soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
  }

  // Produce barrel + muon tables, with centrality --------------------------------------------------------------------------------------------
  void processFullWithCent(MyEventsWithCent::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                           soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
  }

  // Produce barrel + muon tables, with centrality and multiplicity ---------------------------------------------------------------------------
  void processFullWithCentAndMults(MyEventsWithCentAndMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                   soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
  }

  // Produce barrel + muon tables, with covariance, centrality and multiplicity ---------------------------------------------------------------------------
  void processFullWithCovCentAndMults(MyEventsWithCentAndMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                      soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
  }

  // Produce barrel + muon tables, with track covariance matrix and event filtering ----------------------------------------------------------------------------------------
  void processFullWithCovAndEventFilter(MyEventsWithFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                        soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
    }
  }

  // Produce barrel + muon tables, with track covariance matrix and event filtering, with multiplicity -------------------------------------------
  void processFullWithCovMultsAndEventFilter(MyEventsWithMultsAndFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                             soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
    }
  }

  // Produce barrel + muon tables for the eletron-muon analysis, without PIDTOF----------------------------------------------------------------------
  void processFullForElectronMuon(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                  soa::Filtered<MyBarrelTracksForElectronMuon> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapForElectronMuon, gkMuonFillMap>(collision, bcs, tracksBarrel, tracksMuon, nullptr, nullptr);
  }

  // Produce barrel only tables, with V0Bits ------------------------------------------------------------------------------------------------
  void processBarrelOnlyWithV0Bits(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                   soa::Filtered<MyBarrelTracksWithV0Bits> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithV0Bits, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with V0Bits and Mults
  void processBarrelOnlyWithV0BitsAndMults(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                           soa::Filtered<MyBarrelTracksWithV0Bits> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithV0Bits, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with V0Bits and produce maps ------------------------------------------------------------------------------
  void processBarrelOnlyWithV0BitsAndMaps(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                          soa::Filtered<MyBarrelTracksWithV0BitsForMaps> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithV0BitsForMaps, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with DalitzBits ------------------------------------------------------------------------------------------------
  void processBarrelOnlyWithDalitzBits(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                       soa::Filtered<MyBarrelTracksWithDalitzBits> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithDalitzBits, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel tables, with DalitzBits, and skip events without dalitz track  ------------------------------------------------------------------------------
  void processBarrelWithDalitzEvent(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                    soa::Filtered<MyBarrelTracksWithDalitzBits> const& tracksBarrel)
  {
    for (auto const& track : tracksBarrel) {
      if (track.dalitzBits() != uint8_t(0)) {
        fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithDalitzBits, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
        break;
      }
    }
  }

  // Produce barrel only tables, with V0Bits and DalitzBits, and skip event without dalitz/V0 electron -----------------------------------------------------------------
  void processBarrelWithV0AndDalitzEvent(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                         soa::Filtered<MyBarrelTracksWithV0AndDalitzBits> const& tracksBarrel)
  {
    for (auto const& track : tracksBarrel) {
      // check if this event has a Dalitz candidate
      if (track.dalitzBits() != uint8_t(0)) {
        fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithV0AndDalitzBits, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
        break;
      }
      // check if this event has a electron candidate from V0
      if (track.pidbit() == uint8_t(1)) {
        fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithV0AndDalitzBits, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
        break;
      }
    }
  }

  // Produce barrel only tables, with event filtering ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithEventFilter(MyEventsWithFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                        soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
    }
  }

  // Produce barrel tables only, with multiplicity ---------------------------------------------------------------------------------------------
  void processBarrelOnlyWithMults(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                  soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithMult, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with multiplicity and event filtering ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithMultsAndEventFilter(MyEventsWithMultsAndFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                                soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMapWithMultsAndEventFilter, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
    }
  }

  // Produce barrel only tables, with track covariance matrix and event filtering ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithCovAndEventFilter(MyEventsWithFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                              soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
    }
  }

  // Produce barrel only tables, with centrality -----------------------------------------------------------------------------------------------
  void processBarrelOnlyWithCent(MyEventsWithCent::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                 soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with centrality -----------------------------------------------------------------------------------------------
  void processBarrelOnlyWithCovWithCent(MyEventsWithCent::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                        soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMapWithCov, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with centrality and multiplicity -------------------------------------------------------------------
  void processBarrelOnlyWithCentAndMults(MyEventsWithCentAndMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                         soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with centrality and multiplicity -------------------------------------------------------------------
  void processBarrelOnlyWithCovWithCentAndMults(MyEventsWithCentAndMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                                soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel tables only, with track cov matrix ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithCov(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithCov, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel tables only, with track cov matrix , only std PID information used ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithCovOnlyStdPID(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                          soa::Filtered<MyBarrelTracksWithCovOnlyStdPID> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithMult, gkTrackFillMapWithCovOnlyStdPID, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce barrel tables only ----------------------------------------------------------------------------------------------------------------
  void processBarrelOnly(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                         soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr, nullptr, nullptr);
  }

  // Produce muon tables only, with centrality and muon cov matrix -------------------------------------------------------------------------------------------------
  void processMuonOnlyWithCovAndCentMults(MyEventsWithCentAndMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                          soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, 0u, gkMuonFillMapWithCov>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce muon tables only, with muon cov matrix --------------------------------------------------------------------------------------------
  void processMuonOnlyWithCov(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                              soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce muon tables only, with muon cov matrix, with event filtering --------------------------------------------------------------------------------------------
  void processMuonOnlyWithCovAndEventFilter(MyEventsWithFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                            soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
    }
  }

  // Produce muon tables only ------------------------------------------------------------------------------------------------------------------
  void processMuonMLOnly(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                         soa::Filtered<aod::FwdTracksML> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce MFT tracks tables and muons  ------------------------------------------------------------------------------------------------------------------
  void processMuonsAndMFT(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                          aod::MFTTracks const& tracksMft, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr, tracksMft);
  }

  // Produce MFT tracks tables and muons, with event filtering  ------------------------------------------------------------------------------------------------------------------
  void processMuonsAndMFTWithFilter(MyEventsWithFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                    aod::MFTTracks const& tracksMft, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    if (collision.eventFilter()) // the collision has been selected by the filterPP task for at least one of the selections
    {
      fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr, tracksMft);
    }
  }

  void processAssociatedMuonOnlyWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                        soa::Filtered<MyMuonsCollWithCov> const& tracksMuon, aod::AmbiguousFwdTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      fullSkimmingIndices<gkEventFillMap, 0u, gkMuonFillMapWithCovAmbi>(collision, bcs, nullptr, tracksMuon, nullptr, muonIdsThisCollision);
    }
  }

  void processAssociatedRealignedMuonOnlyWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                                 soa::Filtered<MyMuonsRealignCollWithCov> const& tracksMuon, aod::AmbiguousFwdTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      fullSkimmingIndices<gkEventFillMap, 0u, gkMuonRealignFillMapWithCovAmbi>(collision, bcs, nullptr, tracksMuon, nullptr, muonIdsThisCollision);
    }
  }

  void processAssociatedMuonOnlyWithCovAndCentMults(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                                    soa::Filtered<MyMuonsCollWithCov> const& tracksMuon, aod::AmbiguousFwdTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      fullSkimmingIndices<gkEventFillMapWithCentAndMults, 0u, gkMuonFillMapWithCovAmbi>(collision, bcs, nullptr, tracksMuon, nullptr, muonIdsThisCollision);
    }
  }

  void processAssociatedMuonOnlyWithCovAndMults(MyEventsWithMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                                soa::Filtered<MyMuonsCollWithCov> const& tracksMuon, aod::AmbiguousFwdTracks const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      fullSkimmingIndices<gkEventFillMapWithMult, 0u, gkMuonFillMapWithCovAmbi>(collision, bcs, nullptr, tracksMuon, nullptr, muonIdsThisCollision);
    }
  }

  void processAmbiguousMuonOnlyWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                       soa::Filtered<MyMuonsWithCov> const& tracksMuon, aod::AmbiguousFwdTracks const& ambiTracksFwd)
  {
    // Process orphan tracks
    if (fDoDetailedQA && fIsAmbiguous) {
      for (auto& ambiTrackFwd : ambiTracksFwd) {
        auto muon = ambiTrackFwd.template fwdtrack_as<MyMuonsWithCov>();
        if (muon.collisionId() < 0) {
          VarManager::FillTrack<gkMuonFillMapWithCovAmbi>(muon);
          if ((static_cast<int>(muon.trackType()) == 0)) {
            fHistMan->FillHistClass("Orphan_Muons_MFTMCHMID", VarManager::fgValues);
          } else if ((static_cast<int>(muon.trackType()) == 3)) {
            fHistMan->FillHistClass("Orphan_Muons_MCHMID", VarManager::fgValues);
          }
        }
      }
    }
    for (auto& collision : collisions) {
      auto groupedMuons = tracksMuon.sliceBy(perCollisionMuons, collision.globalIndex());
      fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCovAmbi>(collision, bcs, nullptr, groupedMuons, nullptr, ambiTracksFwd);
    }
  }

  // Produce track tables only for ambiguous tracks studies -------------------------------------------------------------------------------------
  void processAmbiguousBarrelOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                  soa::Filtered<MyBarrelTracks> const& tracksBarrel, aod::AmbiguousTracks const& ambiTracksMid)
  {
    // Process orphan tracks
    if (fDoDetailedQA && fIsAmbiguous) {
      for (auto& ambiTrack : ambiTracksMid) {
        auto trk = ambiTrack.template track_as<MyBarrelTracks>();
        if (trk.collisionId() < 0) {
          VarManager::FillTrack<gkTrackFillMapWithAmbi>(trk);
          fHistMan->FillHistClass("Orphan_TrackBarrel", VarManager::fgValues);
        }
      }
    }
    for (auto& collision : collisions) {
      auto groupedTracks = tracksBarrel.sliceBy(perCollisionTracks, collision.globalIndex());
      fullSkimming<gkEventFillMap, gkTrackFillMapWithAmbi, 0u>(collision, bcs, groupedTracks, nullptr, ambiTracksMid, nullptr);
    }
  }

  // Process the BCs and store stats for luminosity retrieval -----------------------------------------------------------------------------------
  void processOnlyBCs(soa::Join<aod::BCs, aod::BcSels>::iterator const& bc)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (bc.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(kNaliases));
  }

  // List of process functions removed because they are not using cov matrix
  /*
  // Produce muon tables only, with centrality -------------------------------------------------------------------------------------------------
  void processMuonOnlyWithCent(MyEventsWithCent::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                               soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCent, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce muon tables only, with multiplicity ---------------------------------------------------------------------------------------------
  void processMuonOnlyWithMults(MyEventsWithMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithMult, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce muon tables only, with centrality and multiplicity --------------------------------------------------------------------------------
  void processMuonOnlyWithCentAndMults(MyEventsWithCentAndMults::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                       soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce muon tables only ------------------------------------------------------------------------------------------------------------------
  void processMuonOnly(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                       soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
  }

  // Produce muon tables only, with event filtering --------------------------------------------------------------------------------------------
  void processMuonOnlyWithEventFilter(MyEventsWithFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                      soa::Filtered<MyMuons> const& tracksMuon)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
    }
  }

  // Produce muon tables only, with multiplicity and event filtering ----------------------------------------------------------------------------------------
  void processMuonOnlyWithMultsAndEventFilter(MyEventsWithMultsAndFilter::iterator const& collision, aod::BCsWithTimestamps const& bcs,
                                              soa::Filtered<MyMuons> const& tracksMuon)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias_bit(i) > 0) {
        (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2F*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMapWithMultsAndEventFilter, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon, nullptr, nullptr);
    }
  }

  // Produce muon tables only based on track-collision association tables --------------------------------------------------------------------------------------
  void processAssociatedMuonOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                 soa::Filtered<MyMuonsColl> const& tracksMuon, aod::AmbiguousTracksFwd const&, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      fullSkimmingIndices<gkEventFillMap, 0u, gkMuonFillMapWithAmbi>(collision, bcs, nullptr, tracksMuon, nullptr, muonIdsThisCollision);
    }
  }
  */

  PROCESS_SWITCH(TableMaker, processFull, "Build full DQ skimmed data model, w/o centrality", false);
  PROCESS_SWITCH(TableMaker, processFullWithCov, "Build full DQ skimmed data model, w/ track and fwdtrack covariance tables", false);
  PROCESS_SWITCH(TableMaker, processFullWithCovAndEventFilter, "Build full DQ skimmed data model, w/ track and fwdtrack covariance tables, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processFullWithCovMultsAndEventFilter, "Build full DQ skimmed data model, w/ track and fwdtrack covariance tables, w/ event filter and multiplicities", false);
  PROCESS_SWITCH(TableMaker, processFullWithCent, "Build full DQ skimmed data model, w/ centrality", false);
  PROCESS_SWITCH(TableMaker, processFullWithCentAndMults, "Build full DQ skimmed data model, w/ centrality and multiplicities", false);
  PROCESS_SWITCH(TableMaker, processFullWithCovCentAndMults, "Build full DQ skimmed data model, w/ centrality, multiplicities and track covariances", false);
  PROCESS_SWITCH(TableMaker, processFullForElectronMuon, "Build full DQ skimmed data model for electron-muon correlation analysis, w/o centrality", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithV0Bits, "Build full DQ skimmed data model, w/o centrality, w/ V0Bits", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithV0BitsAndMults, "Build full DQ skimmed data model, w/ multiplicity, w/ V0Bits", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithV0BitsAndMaps, "Build full DQ skimmed data model, w/o multiplicity, w/ V0Bits", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithDalitzBits, "Build barrel-only DQ skimmed data model, w/o centrality, w/ DalitzBits", false);
  PROCESS_SWITCH(TableMaker, processBarrelWithDalitzEvent, "Build barrel-only DQ skimmed data model, w/o centrality, w/ DalitzBits", false);
  PROCESS_SWITCH(TableMaker, processBarrelWithV0AndDalitzEvent, "Build barrel-only DQ skimmed data model, w/o centrality, w/ V0Bits and DalitzBits", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithEventFilter, "Build full DQ skimmed data model, w/o centrality, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithMults, "Build barrel-only DQ skimmed data model, w/ multiplicity", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithMultsAndEventFilter, "Build barrel-only DQ skimmed data model, w/ multiplicity, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCovAndEventFilter, "Build full DQ skimmed data model, w/ track and fwdtrack covariance tables, w/o centrality, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCent, "Build barrel-only DQ skimmed data model, w/ centrality", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCovWithCent, "Build barrel-only DQ skimmed data model, w/ centrality and w/ track covariance", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCentAndMults, "Build barrel-only DQ skimmed data model, w/ centrality and multiplicities", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCovWithCentAndMults, "Build barrel-only DQ skimmed data model, w/ centrality and multiplicities and w/ track covariance", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCov, "Build barrel-only DQ skimmed data model, w/ track cov matrix", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCovOnlyStdPID, "Build barrel-only DQ skimmed data model, w/ track cov matrix, only std PID information used", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnly, "Build barrel-only DQ skimmed data model, w/o centrality", false);
  PROCESS_SWITCH(TableMaker, processMuonOnlyWithCovAndCentMults, "Build muon-only DQ skimmed data model, w/ centrality and muon cov matrix", false);
  PROCESS_SWITCH(TableMaker, processMuonOnlyWithCov, "Build muon-only DQ skimmed data model, w/ muon cov matrix", false);
  PROCESS_SWITCH(TableMaker, processMuonOnlyWithCovAndEventFilter, "Build muon-only DQ skimmed data model, w/ muon cov matrix, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processMuonMLOnly, "Build muon-only DQ skimmed data model with global muon track by ML matching", false);
  PROCESS_SWITCH(TableMaker, processOnlyBCs, "Analyze the BCs to store sampled lumi", false);
  PROCESS_SWITCH(TableMaker, processAmbiguousMuonOnlyWithCov, "Build muon-only with cov DQ skimmed data model with QA plots for ambiguous muons", false);
  PROCESS_SWITCH(TableMaker, processAmbiguousBarrelOnly, "Build barrel-only DQ skimmed data model with QA plots for ambiguous tracks", false);
  PROCESS_SWITCH(TableMaker, processAssociatedMuonOnlyWithCov, "Build muon-only with cov DQ skimmed data model using track-collision association tables", false);
  PROCESS_SWITCH(TableMaker, processAssociatedRealignedMuonOnlyWithCov, "Build realigned muon-only with cov DQ skimmed data model using track-collision association tables", false);
  PROCESS_SWITCH(TableMaker, processAssociatedMuonOnlyWithCovAndCentMults, "Build muon-only with cov DQ skimmed data model using track-collision association tables and centrality", false);
  PROCESS_SWITCH(TableMaker, processAssociatedMuonOnlyWithCovAndMults, "Build muon-only with cov DQ skimmed data model using track-collision association tables and multiplicity", false);
  PROCESS_SWITCH(TableMaker, processMuonsAndMFT, "Build MFT and muons DQ skimmed data model", false);
  PROCESS_SWITCH(TableMaker, processMuonsAndMFTWithFilter, "Build MFT and muons DQ skimmed data model, w/ event filter", false);

  // List of process functions removed because they are not using cov matrix
  // PROCESS_SWITCH(TableMaker, processMuonOnlyWithCent, "Build muon-only DQ skimmed data model, w/ centrality", false);
  // PROCESS_SWITCH(TableMaker, processMuonOnlyWithMults, "Build muon-only DQ skimmed data model, w/ multiplicity", false);
  // PROCESS_SWITCH(TableMaker, processMuonOnlyWithMultsAndEventFilter, "Build muon-only DQ skimmed data model, w/ multiplicity, w/ event filter", false);
  // PROCESS_SWITCH(TableMaker, processMuonOnlyWithCentAndMults, "Build muon-only DQ skimmed data model, w/ centrality and multiplicities", false);
  // PROCESS_SWITCH(TableMaker, processMuonOnly, "Build muon-only DQ skimmed data model", false);
  // PROCESS_SWITCH(TableMaker, processMuonOnlyWithEventFilter, "Build muon-only DQ skimmed data model, w/ event filter", false);
  // PROCESS_SWITCH(TableMaker, processAssociatedMuonOnly, "Build muon-only DQ skimmed data model using track-collision association tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TableMaker>(cfgc)};
}
