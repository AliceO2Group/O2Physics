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
// Events to be written are filtered using a user provided event cut and optionally the filterPPwithAssociation task
// Barrel and muon (collision-track) associations are filtered using multiple parallel selections (currently limited to 8 but easily extendable up to 64)
// The skimming can optionally produce just the barrel, muon, or both barrel and muon tracks
// The event filtering, centrality, and V0Bits (from v0-selector) can be switched on/off by selecting one
//  of the process functions
#include <iostream>
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/MftmchMatchingML.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DetectorsVertexing/VertexTrackMatcher.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "MathUtils/Primitive2D.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// TODO: Since DCA depends on which collision the track is associated to, we should remove writing and subscribing to DCA tables, to optimize on CPU / memory
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithV0Bits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                           aod::pidTPCFullKa, aod::pidTPCFullPr,
                                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                           aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::V0Bits>;
using MyBarrelTracksWithV0BitsForMaps = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                                  aod::pidTPCFullEl, aod::pidTPCFullPi,
                                                  aod::pidTPCFullKa, aod::pidTPCFullPr, aod::V0Bits>;
using MyBarrelTracksWithDalitzBits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA,
                                               aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                               aod::pidTPCFullKa, aod::pidTPCFullPr,
                                               aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                               aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::DalitzBits>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyEventsWithFilter = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventFilter>;
using MyEventsWithMultsAndFilter = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::DQEventFilter>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using MyEventsWithCentAndMults = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA>;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsCollWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

namespace o2::aod
{
DECLARE_SOA_TABLE(AmbiguousTracksMid, "AOD", "AMBIGUOUSTRACK", //! Table for tracks which are not uniquely associated with a collision
                  o2::soa::Index<>, o2::aod::ambiguous::TrackId, o2::aod::ambiguous::BCIdSlice, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(AmbiguousTracksFwd, "AOD", "AMBIGUOUSFWDTR", //! Table for Fwd tracks which are not uniquely associated with a collision
                  o2::soa::Index<>, o2::aod::ambiguous::FwdTrackId, o2::aod::ambiguous::BCIdSlice, o2::soa::Marker<2>);
} // namespace o2::aod

// constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
// constexpr static uint32_t gkEventFillMapWithFilter = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::EventFilter;

// constexpr static uint32_t gkEventFillMapWithMult = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult;
constexpr static uint32_t gkEventFillMapWithMultsAndEventFilter = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult | VarManager::ObjTypes::EventFilter;
// constexpr static uint32_t gkEventFillMapWithCent = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkEventFillMapWithCentAndMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent | VarManager::CollisionMult;
//  constexpr static uint32_t gkEventFillMapWithCentRun2 = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCentRun2; // Unused variable
// constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
// constexpr static uint32_t gkTrackFillMapWithV0Bits = gkTrackFillMapWithCov | VarManager::ObjTypes::TrackV0Bits;
// constexpr static uint32_t gkTrackFillMapWithV0BitsForMaps = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackV0Bits | VarManager::ObjTypes::TrackTPCPID;
// constexpr static uint32_t gkTrackFillMapWithDalitzBits = gkTrackFillMap | VarManager::ObjTypes::DalitzBits;
// constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;
// constexpr static uint32_t gkMuonFillMapWithAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::AmbiMuon;
// constexpr static uint32_t gkMuonFillMapWithCovAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov | VarManager::ObjTypes::AmbiMuon;
// constexpr static uint32_t gkTrackFillMapWithAmbi = VarManager::ObjTypes::Track | VarManager::ObjTypes::AmbiTrack;
constexpr static uint32_t gkMFTFillMap = VarManager::ObjTypes::TrackMFT;

struct TableMaker {

  Produces<ReducedEvents> event;
  Produces<ReducedEventsExtended> eventExtended;
  Produces<ReducedEventsVtxCov> eventVtxCov;
  Produces<ReducedTracks> trackBasic;
  Produces<ReducedTracksBarrel> trackBarrel;
  Produces<ReducedTracksBarrelCov> trackBarrelCov;
  Produces<ReducedTracksBarrelPID> trackBarrelPID;
  Produces<ReducedTracksBarrelInfo> trackBarrelInfo;
  Produces<ReducedTracksAssoc> trackBarrelAssoc;
  Produces<ReducedMuons> muonBasic;
  Produces<ReducedMuonsExtra> muonExtra;
  Produces<ReducedMuonsCov> muonCov;
  Produces<ReducedMuonsInfo> muonInfo;
  Produces<ReducedMuonsAssoc> muonAssoc;
  Produces<ReducedMFTs> mftTrack;
  Produces<ReducedMFTsExtra> mftTrackExtra;
  Produces<ReducedMFTAssoc> mftAssoc;

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics
  HistogramManager* fHistMan;

  // Event and track AnalysisCut configurables
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};

  // Steer QA output
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes)"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};

  // Selections to be applied as Filter on the Track and FwdTrack
  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};
  Configurable<float> fConfigBarrelTrackMaxAbsEta{"cfgBarrelMaxAbsEta", 0.9f, "Eta absolute value cut for tracks in the barrel"};
  Configurable<float> fConfigBarrelTrackMinPt{"cfgBarrelMinPt", 0.5f, "Minimum pt for tracks in the barrel"};
  Configurable<bool> fConfigBarrelRequireTPC{"cfgBarrelRequireTPC", true, "Require TPC for tracks in the barrel"};
  Configurable<float> fConfigBarrelMinTPCncls{"cfgBarrelMinTPCncls", 50.0f, "Minimum TPC cls for tracks in the barrel"};
  Configurable<float> fConfigBarrelMaxTPCchi2{"cfgBarrelMaxTPCchi2", 10.0f, "Maximum TPC chi2/ndf for tracks in the barrel"};
  Configurable<bool> fConfigBarrelRequireITS{"cfgBarrelRequireITS", true, "Require ITS for tracks in the barrel"};
  Configurable<float> fConfigBarrelMaxITSchi2{"cfgBarrelMaxITSchi2", 36.0f, "Maximum ITS chi2/ndf for tracks in the barrel"};
  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 1.0f, "Low pt cut for muons"};

  // CCDB connection configurables
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/z/zhxiong/TPCPID/PostCalib", "base path to the ccdb object"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpmagPathRun2{"grpmagPathRun2", "GLO/GRP/GRP", "CCDB path of the GRPObject (Usage for Run 2)"};

  // TPC postcalibration related options
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas(electrons, pions, protons)"};
  Configurable<bool> fConfigComputeTPCpostCalibKaon{"cfgTPCpostCalibKaon", false, "If true, compute TPC post-calibrated n-sigmas for kaons"};
  Configurable<bool> fConfigIsOnlyforMaps{"cfgIsforMaps", false, "If true, run for postcalibration maps only"};
  Configurable<bool> fConfigSaveElectronSample{"cfgSaveElectronSample", false, "If true, only save electron sample"};
  Configurable<bool> fConfigDummyRunlist{"cfgDummyRunlist", false, "If true, use dummy runlist"};
  Configurable<int> fConfigInitRunNumber{"cfgInitRunNumber", 543215, "Initial run number used in run by run checks"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  o2::parameters::GRPObject* grpmagrun2 = nullptr; // for run 2, we access the GRPObject from GLO/GRP/GRP
  o2::parameters::GRPMagField* grpmag = nullptr;   // for run 3, we access GRPMagField from GLO/Config/GRPMagField

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fMuonCuts;  //! Muon track cuts

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  bool fDoDetailedQA = false; // Bool to set detailed QA true, if QA is set true
  int fCurrentRun;            // needed to detect if the run changed and trigger update of calibrations etc.

  // maps used to store index info; NOTE: std::map are sorted in ascending order by default (needed for track to collision indices)
  std::map<uint32_t, uint32_t> fCollIndexMap;             // key: old collision index, value: skimmed collision index
  std::map<uint32_t, uint32_t> fTrackIndexMap;            // key: old track global index, value: new track global index
  std::map<uint32_t, uint32_t> fFwdTrackIndexMap;         // key: fwd-track global index, value: new fwd-track global index
  std::map<uint32_t, uint32_t> fFwdTrackIndexMapReversed; // key: new fwd-track global index, value: fwd-track global index
  std::map<uint32_t, uint8_t> fFwdTrackFilterMap;         // key: fwd-track global index, value: fwd-track filter map
  std::map<uint32_t, uint32_t> fMftIndexMap;              // key: MFT tracklet global index, value: new MFT tracklet global index

  // FIXME: For now, the skimming is done using the Common track-collision association task, which does not allow to use
  //       our own Filtered tracks. If the filter is very selective, then it may be worth to run the association in this workflow
  //       using the Common/CollisionAssociation class
  /*Filter barrelSelectedTracks = ifnode(fIsRun2.node() == true, track::trackType == uint8_t(track::Run2Track), track::trackType == uint8_t(track::Track))
                              && track::pt > fConfigBarrelTrackMinPt
                              && nabs(track::eta) <= fConfigBarrelTrackMaxAbsEta
                              && ifnode(fConfigBarrelRequireITS.node() == true, track::itsChi2NCl < fConfigBarrelMaxITSchi2, true)
                              && ifnode(fConfigBarrelRequireTPC.node() == true, track::tpcNClsFound > fConfigBarrelMinTPCncls && track::tpcChi2NCl < fConfigBarrelMaxTPCchi2, true);

  Filter muonFilter = o2::aod::fwdtrack::pt >= fConfigMuonPtLow;*/

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::MFTTrackAssoc> mfttrackIndicesPerCollision = aod::track_association::collisionId;

  void init(o2::framework::InitContext& context)
  {
    DefineCuts();
    ccdb->setURL(fConfigCcdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
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

    bool enableBarrelHistos = (context.mOptions.get<bool>("processPPWithFilter") || context.mOptions.get<bool>("processPPWithFilterBarrelOnly") ||
                               context.mOptions.get<bool>("processPbPb") || context.mOptions.get<bool>("processPbPbBarrelOnly"));
    bool enableMuonHistos = (context.mOptions.get<bool>("processPPWithFilter") || context.mOptions.get<bool>("processPPWithFilterMuonOnly") || context.mOptions.get<bool>("processPPWithFilterMuonMFT") ||
                             context.mOptions.get<bool>("processPbPb") || context.mOptions.get<bool>("processPbPbMuonOnly") || context.mOptions.get<bool>("processPbPbMuonMFT"));

    if (enableBarrelHistos) {
      if (fDoDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
      }
      if (fConfigQA) {
        for (auto& cut : fTrackCuts) {
          histClasses += Form("TrackBarrel_%s;", cut.GetName());
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
      }
      if (fConfigQA) {
        for (auto& muonCut : fMuonCuts) {
          histClasses += Form("Muons_%s;", muonCut.GetName());
        }
      }
    }

    if (fConfigDummyRunlist) {
      VarManager::SetDummyRunlist(fConfigInitRunNumber);
    }

    DefineHistograms(histClasses);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    // CCDB configuration
    if (fConfigComputeTPCpostCalib) {
      fCCDB->setURL(fConfigCcdbUrl.value);
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
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    // Barrel track cuts
    TString cutNamesStr = fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // Muon cuts
    cutNamesStr = fConfigMuonCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

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

      TString histEventName = fConfigAddEventHistogram.value;
      if (classStr.Contains("Event")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", histEventName);
        }
      }

      TString histTrackName = fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
        }
      }

      TString histMuonName = fConfigAddMuonHistogram.value;
      if (classStr.Contains("Muons")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histMuonName);
        }
      }

      TString histMftName = fConfigAddMuonHistogram.value;
      if (classStr.Contains("Mft")) {
        if (fConfigDetailedQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histMftName);
        }
      }
    }

    // create statistics histograms (event, tracks, muons)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, o2::aod::evsel::kNsel + 1, -0.5, double(o2::aod::evsel::kNsel) + 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    for (int ib = 1; ib <= o2::aod::evsel::kNsel; ib++) {
      histEvents->GetYaxis()->SetBinLabel(ib, o2::aod::evsel::selectionLabels[ib - 1]);
    }
    histEvents->GetYaxis()->SetBinLabel(o2::aod::evsel::kNsel + 1, "Total");
    fStatsList->Add(histEvents);

    // Track statistics: one bin for each track selection and 5 bins for V0 tags (gamma, K0s, Lambda, anti-Lambda, Omega)
    TH1I* histTracks = new TH1I("TrackStats", "Track statistics", fTrackCuts.size() + 5.0, -0.5, fTrackCuts.size() - 0.5 + 5.0);
    ib = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ib++) {
      histTracks->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    const char* v0TagNames[5] = {"Photon conversion", "K^{0}_{s}", "#Lambda", "#bar{#Lambda}", "#Omega"};
    for (ib = 0; ib < 5; ib++) {
      histTracks->GetXaxis()->SetBinLabel(fTrackCuts.size() + 1 + ib, v0TagNames[ib]);
    }
    fStatsList->Add(histTracks);
    TH1I* histMuons = new TH1I("MuonStats", "Muon statistics", fMuonCuts.size(), -0.5, fMuonCuts.size() - 0.5);
    ib = 1;
    for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, ib++) {
      histMuons->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    fStatsList->Add(histMuons);
  }

  template <uint32_t TEventFillMap, typename TEvents>
  void skimCollisions(TEvents const& collisions, BCsWithTimestamps const& /*bcs*/)
  {
    // Skim collisions
    // NOTE: So far, collisions are filtered based on the user specified analysis cuts and the filterPP event filter.
    //      The collision-track associations which point to an event that is not selected for writing are discarded!

    fCollIndexMap.clear();
    int multTPC = -1.0;
    float multFV0A = -1.0;
    float multFV0C = -1.0;
    float multFT0A = -1.0;
    float multFT0C = -1.0;
    float multFDDA = -1.0;
    float multFDDC = -1.0;
    float multZNA = -1.0;
    float multZNC = -1.0;
    int multTracklets = -1.0;
    int multTracksPV = -1.0;
    float centFT0C = -1.0;

    for (const auto& collision : collisions) {

      for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
        if (collision.selection_bit(i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(o2::aod::evsel::kNsel));

      // apply the event filter computed by filter-PP
      if constexpr ((TEventFillMap & VarManager::ObjTypes::EventFilter) > 0) {
        if (!collision.eventFilter()) {
          continue;
        }
      }

      auto bc = collision.template bc_as<BCsWithTimestamps>();
      // store the selection decisions
      uint64_t tag = 0;
      // store some more information in the tag
      // if the BC found by event selection does not coincide with the collision.bc()
      auto bcEvSel = collision.template foundBC_as<BCsWithTimestamps>();
      if (bcEvSel.globalIndex() != bc.globalIndex()) {
        tag |= (uint64_t(1) << 0);
      }
      // Put the 8 first bits of the event filter in the last 8 bits of the tag
      if constexpr ((TEventFillMap & VarManager::ObjTypes::EventFilter) > 0) {
        tag |= (collision.eventFilter() << 56);
      }

      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      // TODO: These variables cannot be filled in the VarManager for the moment as long as BCsWithTimestamps are used.
      //       So temporarily, we filled them here, in order to be available for eventual QA of the skimming
      VarManager::fgValues[VarManager::kRunNo] = bc.runNumber();
      VarManager::fgValues[VarManager::kBC] = bc.globalBC();
      VarManager::fgValues[VarManager::kTimestamp] = bc.timestamp();
      VarManager::fgValues[VarManager::kRunIndex] = VarManager::GetRunIndex(bc.runNumber());
      VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      // fill stats information, before selections
      for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
        if (collision.selection_bit(i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(o2::aod::evsel::kNsel));

      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }

      // fill stats information, after selections
      for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
        if (collision.selection_bit(i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(o2::aod::evsel::kNsel));

      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

      // create the event tables
      event(tag, bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0) {
        multTPC = collision.multTPC();
        multFV0A = collision.multFV0A();
        multFV0C = collision.multFV0C();
        multFT0A = collision.multFT0A();
        multFT0C = collision.multFT0C();
        multFDDA = collision.multFDDA();
        multFDDC = collision.multFDDC();
        multZNA = collision.multZNA();
        multZNC = collision.multZNC();
        multTracklets = collision.multTracklets();
        multTracksPV = collision.multNTracksPV();
      }
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
        centFT0C = collision.centFT0C();
      }
      eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                    multTPC, multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTracksPV, centFT0C);
      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      fCollIndexMap[collision.globalIndex()] = event.lastIndex();
    }
  }

  template <uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void skimTracks(TEvent const& collision, BCsWithTimestamps const& /*bcs*/, TTracks const& /*tracks*/, TrackAssoc const& assocs)
  {
    // Skim the barrel tracks
    // Loop over the collision-track associations, retrieve the track, and apply track cuts for selection
    //     One can apply here cuts which depend on the association (e.g. DCA), which will discard (hopefully most) wrong associations.
    //     Tracks are written only once, even if they constribute to more than one association

    uint64_t trackFilteringTag = uint64_t(0);
    uint64_t trackTempFilterMap = uint8_t(0);

    for (const auto& assoc : assocs) {
      auto track = assoc.template track_as<TTracks>();

      trackFilteringTag = uint64_t(0);
      trackTempFilterMap = uint8_t(0);
      VarManager::FillTrack<TTrackFillMap>(track);
      if (fDoDetailedQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      // apply track cuts and fill stats histogram
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          trackTempFilterMap |= (uint8_t(1) << i);
          // NOTE: the QA is filled here for every (collision,track) association since the results of the track cuts can be different depending on which collision is associated (e.g. DCA cuts)
          // TODO: create a statistics histograms with unique tracks
          if (fConfigQA) {
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
          }
          (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
        }
      }
      if (!trackTempFilterMap) {
        continue;
      }

      // store selection information in the track tag
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackV0Bits)) { // BIT0-4: V0Bits
        trackFilteringTag |= uint64_t(track.pidbit());
        for (int iv0 = 0; iv0 < 5; iv0++) {
          if (track.pidbit() & (uint8_t(1) << iv0)) {
            (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(fTrackCuts.size() + static_cast<float>(iv0));
          }
        }
        if (fConfigIsOnlyforMaps) {
          if (trackFilteringTag & (uint64_t(1) << VarManager::kIsConversionLeg)) { // for electron
            fHistMan->FillHistClass("TrackBarrel_PostCalibElectron", VarManager::fgValues);
          }
          if (trackFilteringTag & (uint64_t(1) << VarManager::kIsK0sLeg)) { // for pion
            fHistMan->FillHistClass("TrackBarrel_PostCalibPion", VarManager::fgValues);
          }
          if ((static_cast<bool>(trackFilteringTag & (uint64_t(1) << VarManager::kIsLambdaLeg)) * (track.sign()) > 0)) { // for proton from Lambda
            fHistMan->FillHistClass("TrackBarrel_PostCalibProton", VarManager::fgValues);
          }
          if ((static_cast<bool>(trackFilteringTag & (uint64_t(1) << VarManager::kIsALambdaLeg)) * (track.sign()) < 0)) { // for proton from AntiLambda
            fHistMan->FillHistClass("TrackBarrel_PostCalibProton", VarManager::fgValues);
          }
        }
        if (fConfigSaveElectronSample) { // only save electron sample
          if (!(trackFilteringTag & (uint64_t(1) << VarManager::kIsConversionLeg))) {
            continue;
          }
        }
      } // end if V0Bits
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
        trackFilteringTag |= (uint64_t(track.dalitzBits()) << VarManager::kDalitzBits); // BIT5-12: Dalitz
      }
      trackFilteringTag |= (uint64_t(trackTempFilterMap) << VarManager::kBarrelUserCutsBits); // BIT13-...:  user track filters

      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPID)) {
        if (fConfigComputeTPCpostCalib) {
          trackFilteringTag |= (uint64_t(1) << VarManager::kIsTPCPostcalibrated);
        }
      }

      // write the track global index in the map for skimming (to make sure we have it just once)
      if (fTrackIndexMap.find(track.globalIndex()) == fTrackIndexMap.end()) {
        // NOTE: The collision ID that is written in the table is the one found in the first association for this track.
        //       However, in data analysis one should loop over associations, so this one should not be used.
        //      In the case of Run2-like analysis, there will be no associations, so this ID will be the one originally assigned in the AO2Ds (updated for the skims)
        uint32_t reducedEventIdx = fCollIndexMap[collision.globalIndex()];
        // NOTE: trackBarrelInfo stores the index of the collision as in AO2D (for use in some cases where the analysis on skims is done
        //   in workflows where the original AO2Ds are also present)
        trackBarrelInfo(collision.globalIndex(), collision.posX(), collision.posY(), collision.posZ());
        trackBasic(reducedEventIdx, trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), 0);
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
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackPID)) {
          float nSigmaEl = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaEl_Corr] : track.tpcNSigmaEl());
          float nSigmaPi = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaPi_Corr] : track.tpcNSigmaPi());
          float nSigmaKa = ((fConfigComputeTPCpostCalib && fConfigComputeTPCpostCalibKaon) ? VarManager::fgValues[VarManager::kTPCnSigmaKa_Corr] : track.tpcNSigmaKa());
          float nSigmaPr = (fConfigComputeTPCpostCalib ? VarManager::fgValues[VarManager::kTPCnSigmaPr_Corr] : track.tpcNSigmaPr());
          trackBarrelPID(track.tpcSignal(),
                         nSigmaEl, track.tpcNSigmaMu(), nSigmaPi, nSigmaKa, nSigmaPr,
                         track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.trdSignal());
        }
        fTrackIndexMap[track.globalIndex()] = trackBasic.lastIndex();
      }
      // write the skimmed collision - track association
      trackBarrelAssoc(fCollIndexMap[collision.globalIndex()], fTrackIndexMap[track.globalIndex()]);
    } // end loop over associations
  }   // end skimTracks

  template <uint32_t TMFTFillMap, typename TEvent>
  void skimMFT(TEvent const& collision, BCsWithTimestamps const& /*bcs*/, MFTTracks const& /*mfts*/, MFTTrackAssoc const& mftAssocs)
  {
    // Skim MFT tracks
    // So far no cuts are applied here

    for (const auto& assoc : mftAssocs) {
      auto track = assoc.template mfttrack_as<MFTTracks>();

      if (fConfigQA) {
        VarManager::FillTrack<TMFTFillMap>(track);
        fHistMan->FillHistClass("MftTracks", VarManager::fgValues);
      }

      // write the MFT track global index in the map for skimming (to make sure we have it just once)
      if (fMftIndexMap.find(track.globalIndex()) == fMftIndexMap.end()) {
        uint32_t reducedEventIdx = fCollIndexMap[collision.globalIndex()];
        mftTrack(reducedEventIdx, uint64_t(0), track.pt(), track.eta(), track.phi());
        // TODO: We are not writing the DCA at the moment, because this depend on the collision association
        mftTrackExtra(track.mftClusterSizesAndTrackFlags(), track.sign(), 0.0, 0.0, track.nClusters());

        fMftIndexMap[track.globalIndex()] = mftTrack.lastIndex();
      }
      mftAssoc(fCollIndexMap[collision.globalIndex()], fMftIndexMap[track.globalIndex()]);
    }
  }

  template <uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void skimMuons(TEvent const& collision, BCsWithTimestamps const& /*bcs*/, TMuons const& muons, FwdTrackAssoc const& muonAssocs)
  {
    // Skim the fwd-tracks (muons)
    // Loop over the collision-track associations, recompute track properties depending on the collision assigned, and apply track cuts for selection
    //     Muons are written only once, even if they constribute to more than one association,
    //         which means that in the case of multiple associations, the track parameters are wrong and should be computed again at analysis time.

    uint8_t trackFilteringTag = uint8_t(0);
    uint8_t trackTempFilterMap = uint8_t(0);
    fFwdTrackIndexMapReversed.clear();

    uint32_t offset = muonBasic.lastIndex();
    uint32_t counter = 0;
    for (const auto& assoc : muonAssocs) {
      // get the muon
      auto muon = assoc.template fwdtrack_as<TMuons>();

      trackFilteringTag = uint8_t(0);
      VarManager::FillTrack<TMuonFillMap>(muon);
      // NOTE: If a muon is associated to multiple collisions, depending on the selections,
      //       it may be accepted for some associations and rejected for other
      VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
      }
      // check the cuts and filters
      int i = 0;
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          trackTempFilterMap |= (uint8_t(1) << i);
          if (fConfigQA) {
            fHistMan->FillHistClass(Form("Muons_%s", (*cut).GetName()), VarManager::fgValues);
          }
          (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
        }
      }

      // don't skim the muon if no cut has passed
      if (!trackTempFilterMap) {
        continue;
      }
      trackFilteringTag = trackTempFilterMap; // BIT0-7:  user selection cuts

      // update the index map if this is a new muon (it can already exist in the map from a different collision association)
      if (fFwdTrackIndexMap.find(muon.globalIndex()) == fFwdTrackIndexMap.end()) {
        counter++;
        fFwdTrackIndexMap[muon.globalIndex()] = offset + counter;
        fFwdTrackIndexMapReversed[offset + counter] = muon.globalIndex();
        fFwdTrackFilterMap[muon.globalIndex()] = trackFilteringTag;                                                    // store here the filtering tag so we don't repeat the cuts in the second iteration
        if (muon.has_matchMCHTrack() && (fFwdTrackIndexMap.find(muon.matchMCHTrackId()) == fFwdTrackIndexMap.end())) { // write also the matched MCH track
          counter++;
          fFwdTrackIndexMap[muon.matchMCHTrackId()] = offset + counter;
          fFwdTrackIndexMapReversed[offset + counter] = muon.matchMCHTrackId();
          fFwdTrackFilterMap[muon.matchMCHTrackId()] = trackFilteringTag; // store here the filtering tag so we don't repeat the cuts in the second iteration
        }
      } else {
        fFwdTrackFilterMap[muon.globalIndex()] |= trackFilteringTag; // make a bitwise OR with previous existing cuts
      }
      // write the association table
      muonAssoc(fCollIndexMap[collision.globalIndex()], fFwdTrackIndexMap[muon.globalIndex()]);
    } // end loop over assocs

    // Now we have the full index map of selected muons so we can proceed with writing the muon tables
    // Special care needed for the MCH and MFT indices
    for (const auto& [skimIdx, origIdx] : fFwdTrackIndexMapReversed) {
      // get the muon
      auto muon = muons.rawIteratorAt(origIdx);
      uint32_t reducedEventIdx = fCollIndexMap[collision.globalIndex()];
      // NOTE: Currently, one writes the original AO2D momentum-vector (pt, eta and phi) in the tables because we write only one instance of the muon track,
      //       while multiple collision associations (and multiple mom vectors can exist)
      //       The momentum can be recomputed at the analysis time based on the associations written in the skims
      //   So all the information which pertains to collision association or MFT associations should not be taken from the skimmed data, but recomputed at analysis time.
      uint32_t mchIdx = -1;
      uint32_t mftIdx = -1;
      if (muon.trackType() == uint8_t(0) || muon.trackType() == uint8_t(2)) { // MCH-MID (2) or global (0)
        if (fFwdTrackIndexMap.find(muon.matchMCHTrackId()) != fFwdTrackIndexMap.end()) {
          mchIdx = fFwdTrackIndexMap[muon.matchMCHTrackId()];
        }
        if (fMftIndexMap.find(muon.matchMFTTrackId()) != fMftIndexMap.end()) {
          mftIdx = fMftIndexMap[muon.matchMFTTrackId()];
        }
      }
      muonBasic(reducedEventIdx, mchIdx, mftIdx, fFwdTrackFilterMap[muon.globalIndex()], muon.pt(), muon.eta(), muon.phi(), muon.sign(), 0);
      muonExtra(muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
                muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                muon.matchScoreMCHMFT(),
                muon.mchBitMap(), muon.midBitMap(),
                muon.midBoards(), muon.trackType(), muon.fwdDcaX(), muon.fwdDcaY(),
                muon.trackTime(), muon.trackTimeRes());
      muonInfo(muon.collisionId(), collision.posX(), collision.posY(), collision.posZ());
      if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
        muonCov(muon.x(), muon.y(), muon.z(), muon.phi(), muon.tgl(), muon.sign() / muon.pt(),
                muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(), muon.cPhiPhi(),
                muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(), muon.c1PtX(), muon.c1PtY(),
                muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2());
      }
    } // end loop over selected muons
  }   // end skimMuons

  // Produce standard barrel + muon tables with event filter (typically for pp and p-Pb) ------------------------------------------------------
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, uint32_t TMFTFillMap, typename TEvents, typename TTracks,
            typename TMuons, typename TMFTTracks, typename TTrackAssoc, typename TFwdTrackAssoc, typename TMFTTrackAssoc>
  void fullSkimming(TEvents const& collisions, BCsWithTimestamps const& bcs,
                    TTracks const& tracksBarrel, TMuons const& muons, TMFTTracks const& mftTracks,
                    TTrackAssoc const& trackAssocs, TFwdTrackAssoc const& fwdTrackAssocs, TMFTTrackAssoc const& mftAssocs)
  {

    if (fCurrentRun != bcs.begin().runNumber()) {
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, bcs.begin().timestamp());
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
        grpmagrun2 = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpmagPathRun2, bcs.begin().timestamp());
        if (grpmagrun2 != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmagrun2);
        }
      } else {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bcs.begin().timestamp());
        if (grpmag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpmag);
        }
      }
      fCurrentRun = bcs.begin().runNumber();
    }

    // skim collisions
    skimCollisions<TEventFillMap>(collisions, bcs);
    if (fCollIndexMap.size() == 0) {
      return;
    }

    if constexpr (static_cast<bool>(TTrackFillMap)) {
      fTrackIndexMap.clear();
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      trackBarrelInfo.reserve(tracksBarrel.size());
      trackBarrelCov.reserve(tracksBarrel.size());
      trackBarrelPID.reserve(tracksBarrel.size());
      trackBarrelAssoc.reserve(tracksBarrel.size());
    }

    if constexpr (static_cast<bool>(TMFTFillMap)) {
      fMftIndexMap.clear();
      mftTrack.reserve(mftTracks.size());
      mftTrackExtra.reserve(mftTracks.size());
      mftAssoc.reserve(mftTracks.size());
    }

    if constexpr (static_cast<bool>(TMuonFillMap)) {
      fFwdTrackIndexMap.clear();
      fFwdTrackFilterMap.clear();
      muonBasic.reserve(muons.size());
      muonExtra.reserve(muons.size());
      muonInfo.reserve(muons.size());
      muonCov.reserve(muons.size());
      muonAssoc.reserve(muons.size());
    }

    // loop over selected collisions and select the tracks and fwd tracks to be skimmed
    for (auto const& [origIdx, skimIdx] : fCollIndexMap) {
      auto collision = collisions.rawIteratorAt(origIdx);
      // group the tracks and muons for this collision
      if constexpr (static_cast<bool>(TTrackFillMap)) {
        auto groupedTrackIndices = trackAssocs.sliceBy(trackIndicesPerCollision, origIdx);
        skimTracks<TTrackFillMap>(collision, bcs, tracksBarrel, groupedTrackIndices);
      }
      if constexpr (static_cast<bool>(TMFTFillMap)) {
        auto groupedMFTIndices = mftAssocs.sliceBy(mfttrackIndicesPerCollision, origIdx);
        skimMFT<TMFTFillMap>(collision, bcs, mftTracks, groupedMFTIndices);
      }
      if constexpr (static_cast<bool>(TMuonFillMap)) {
        auto groupedMuonIndices = fwdTrackAssocs.sliceBy(fwdtrackIndicesPerCollision, origIdx);
        skimMuons<TMuonFillMap>(collision, bcs, muons, groupedMuonIndices);
      }
    } // end loop over skimmed collisions
  }

  // produce the full DQ skimmed data model typically for pp/p-Pb or UPC Pb-Pb (no centrality), subscribe to the DQ event filter (filter-pp or filter-PbPb)
  void processPPWithFilter(MyEventsWithMultsAndFilter const& collisions, BCsWithTimestamps const& bcs,
                           MyBarrelTracksWithCov const& tracksBarrel,
                           MyMuonsWithCov const& muons, MFTTracks const& mftTracks,
                           TrackAssoc const& trackAssocs, FwdTrackAssoc const& fwdTrackAssocs,
                           MFTTrackAssoc const& mftAssocs)
  {
    fullSkimming<gkEventFillMapWithMultsAndEventFilter, gkTrackFillMapWithCov, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, tracksBarrel, muons, mftTracks, trackAssocs, fwdTrackAssocs, mftAssocs);
  }

  // produce the barrel-only DQ skimmed data model typically for pp/p-Pb or UPC Pb-Pb (no centrality), subscribe to the DQ event filter (filter-pp or filter-PbPb)
  void processPPWithFilterBarrelOnly(MyEventsWithMultsAndFilter const& collisions, BCsWithTimestamps const& bcs,
                                     MyBarrelTracksWithCov const& tracksBarrel,
                                     TrackAssoc const& trackAssocs)
  {
    /*const int& a = 0;
    MFTTracks const& mftTracks = 0;
    FwdTrackAssoc const& fwdTrackAssocs = 0;
    MFTTrackAssoc const& mftAssocs = 0;*/
    fullSkimming<gkEventFillMapWithMultsAndEventFilter, gkTrackFillMapWithCov, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, nullptr, trackAssocs, nullptr, nullptr);
  }

  // produce the muon-only DQ skimmed data model typically for pp/p-Pb or UPC Pb-Pb (no centrality), subscribe to the DQ event filter (filter-pp or filter-PbPb)
  void processPPWithFilterMuonOnly(MyEventsWithMultsAndFilter const& collisions, BCsWithTimestamps const& bcs,
                                   MyMuonsWithCov const& muons, FwdTrackAssoc const& fwdTrackAssocs)
  {
    fullSkimming<gkEventFillMapWithMultsAndEventFilter, 0u, gkMuonFillMapWithCov, 0u>(collisions, bcs, nullptr, muons, nullptr, nullptr, fwdTrackAssocs, nullptr);
  }

  // produce the muon+mft DQ skimmed data model typically for pp/p-Pb or UPC Pb-Pb (no centrality), subscribe to the DQ event filter (filter-pp or filter-PbPb)
  void processPPWithFilterMuonMFT(MyEventsWithMultsAndFilter const& collisions, BCsWithTimestamps const& bcs,
                                  MyMuonsWithCov const& muons, MFTTracks const& mftTracks,
                                  FwdTrackAssoc const& fwdTrackAssocs, MFTTrackAssoc const& mftAssocs)
  {
    fullSkimming<gkEventFillMapWithMultsAndEventFilter, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, nullptr, muons, mftTracks, nullptr, fwdTrackAssocs, mftAssocs);
  }

  // produce the full DQ skimmed data model typically for Pb-Pb (with centrality), no subscribtion to the DQ event filter
  void processPbPb(MyEventsWithCentAndMults const& collisions, BCsWithTimestamps const& bcs,
                   MyBarrelTracksWithCov const& tracksBarrel,
                   MyMuonsWithCov const& muons, MFTTracks const& mftTracks,
                   TrackAssoc const& trackAssocs, FwdTrackAssoc const& fwdTrackAssocs,
                   MFTTrackAssoc const& mftAssocs)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, tracksBarrel, muons, mftTracks, trackAssocs, fwdTrackAssocs, mftAssocs);
  }

  // produce the barrel only DQ skimmed data model typically for Pb-Pb (with centrality), no subscribtion to the DQ event filter
  void processPbPbBarrelOnly(MyEventsWithCentAndMults const& collisions, BCsWithTimestamps const& bcs,
                             MyBarrelTracksWithCov const& tracksBarrel,
                             TrackAssoc const& trackAssocs)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, nullptr, trackAssocs, nullptr, nullptr);
  }

  // produce the muon only DQ skimmed data model typically for Pb-Pb (with centrality), no subscribtion to the DQ event filter
  void processPbPbMuonOnly(MyEventsWithCentAndMults const& collisions, BCsWithTimestamps const& bcs,
                           MyMuonsWithCov const& muons, FwdTrackAssoc const& fwdTrackAssocs)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, 0u, gkMuonFillMapWithCov, 0u>(collisions, bcs, nullptr, muons, nullptr, nullptr, fwdTrackAssocs, nullptr);
  }

  // produce the muon+mft DQ skimmed data model typically for Pb-Pb (with centrality), no subscribtion to the DQ event filter
  void processPbPbMuonMFT(MyEventsWithCentAndMults const& collisions, BCsWithTimestamps const& bcs,
                          MyMuonsWithCov const& muons, MFTTracks const& mftTracks,
                          FwdTrackAssoc const& fwdTrackAssocs, MFTTrackAssoc const& mftAssocs)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, nullptr, muons, mftTracks, nullptr, fwdTrackAssocs, mftAssocs);
  }

  // Process the BCs and store stats for luminosity retrieval -----------------------------------------------------------------------------------
  void processOnlyBCs(soa::Join<aod::BCs, aod::BcSels>::iterator const& bc)
  {
    for (int i = 0; i < o2::aod::evsel::kNsel; i++) {
      if (bc.selection_bit(i) > 0) {
        (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(o2::aod::evsel::kNsel));
  }

  PROCESS_SWITCH(TableMaker, processPPWithFilter, "Build full DQ skimmed data model typically for pp/p-Pb and UPC Pb-Pb, w/ event filtering", false);
  PROCESS_SWITCH(TableMaker, processPPWithFilterBarrelOnly, "Build barrel only DQ skimmed data model typically for pp/p-Pb and UPC Pb-Pb, w/ event filtering", false);
  PROCESS_SWITCH(TableMaker, processPPWithFilterMuonOnly, "Build muon only DQ skimmed data model typically for pp/p-Pb and UPC Pb-Pb, w/ event filtering", false);
  PROCESS_SWITCH(TableMaker, processPPWithFilterMuonMFT, "Build muon + mft DQ skimmed data model typically for pp/p-Pb and UPC Pb-Pb, w/ event filtering", false);
  PROCESS_SWITCH(TableMaker, processPbPb, "Build full DQ skimmed data model typically for Pb-Pb, w/o event filtering", false);
  PROCESS_SWITCH(TableMaker, processPbPbBarrelOnly, "Build barrel only DQ skimmed data model typically for Pb-Pb, w/o event filtering", false);
  PROCESS_SWITCH(TableMaker, processPbPbMuonOnly, "Build muon only DQ skimmed data model typically for Pb-Pb, w/o event filtering", false);
  PROCESS_SWITCH(TableMaker, processPbPbMuonMFT, "Build muon + mft DQ skimmed data model typically for Pb-Pb, w/o event filtering", false);
  PROCESS_SWITCH(TableMaker, processOnlyBCs, "Analyze the BCs to store sampled lumi", false);
};

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    TString histName = configVar.value;
    // NOTE: The level of detail for histogramming can be controlled via configurables
    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
        if (classStr.Contains("PIDCalibElectron")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PIDCalibPion")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PIDCalibProton")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
      }
      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

    if (classStr.Contains("DileptonsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "barrel");
    }

    if (classStr.Contains("HadronsSelected")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
    }

    if (classStr.Contains("DileptonHadronInvMass")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-mass");
    }

    if (classStr.Contains("DileptonHadronCorrelation")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "dilepton-hadron-correlation");
    }
  } // end loop over histogram classes
}

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TableMaker>(cfgc)};
}
