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
// Analysis task for skimming MC AODs
// Similar to tableMaker.cxx. The written skimmed data model includes, besides the reconstructed data tables, a skimmed MC stack.
//   The skimmed MC stack includes the MC truth particles corresponding to the list of user specified MC signals (see MCsignal.h)
//    and the MC truth particles corresponding to the reconstructed tracks selected by the specified track cuts on reconstructed data.

#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include "TList.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CCDB/BasicCCDBManager.h"

using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                 aod::McTrackLabels>;
using MyBarrelTracksWithDalitzBits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                               aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                               aod::pidTPCFullKa, aod::pidTPCFullPr,
                                               aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                               aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                               aod::McTrackLabels, aod::DalitzBits>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                                        aod::McTrackLabels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsCollWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTrkCompColls>;

using MyMftTracks = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsWithMults = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::McCollisionLabels>;
using MyEventsWithCentAndMults = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::MultsExtra, aod::McCollisionLabels>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkEventFillMapWithMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionMult;
constexpr static uint32_t gkEventFillMapWithCent = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkEventFillMapWithCentAndMults = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent | VarManager::CollisionMult;
constexpr static uint32_t gkEventMCFillMap = VarManager::ObjTypes::CollisionMC;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID | VarManager::ObjTypes::TrackPIDExtra;
constexpr static uint32_t gkTrackFillMapWithDalitzBits = gkTrackFillMap | VarManager::ObjTypes::DalitzBits;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;
constexpr static uint32_t gkMuonFillMapWithAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::AmbiMuon;
constexpr static uint32_t gkMuonFillMapWithCovAmbi = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov | VarManager::ObjTypes::AmbiMuon;
constexpr static uint32_t gkTrackFillMapWithAmbi = VarManager::ObjTypes::Track | VarManager::ObjTypes::AmbiTrack;
constexpr static uint32_t gkMFTFillMap = VarManager::ObjTypes::TrackMFT;

struct TableMakerMC {

  Produces<ReducedEvents> event;
  Produces<ReducedEventsExtended> eventExtended;
  Produces<ReducedEventsVtxCov> eventVtxCov;
  Produces<ReducedEventsInfo> eventInfo;
  Produces<ReducedEventsMultPV> multPV;
  Produces<ReducedEventsMultAll> multAll;
  Produces<ReducedMCEventLabels> eventMClabels;
  Produces<ReducedMCEvents> eventMC;
  Produces<ReducedTracksBarrelInfo> trackBarrelInfo;
  Produces<ReducedTracks> trackBasic;
  Produces<ReducedTracksBarrel> trackBarrel;
  Produces<ReducedTracksBarrelCov> trackBarrelCov;
  Produces<ReducedTracksBarrelPID> trackBarrelPID;
  Produces<ReducedTracksBarrelLabels> trackBarrelLabels;
  Produces<ReducedMCTracks> trackMC;
  Produces<ReducedMuons> muonBasic;
  Produces<ReducedMuonsExtra> muonExtra;
  Produces<ReducedMuonsCov> muonCov;
  Produces<ReducedMuonsLabels> muonLabels; // TODO: enable this once we have fwdtrack labels
  Produces<ReducedMFTs> trackMFT;
  Produces<ReducedMFTsExtra> trackMFTExtra;
  Produces<ReducedMFTLabels> trackMFTLabels;

  // list of MCsignal objects
  std::vector<MCSignal> fMCSignals;

  OutputObj<THashList> fOutputList{"output"};
  // TODO: add statistics histograms, similar to table-maker
  OutputObj<TList> fStatsList{"Statistics"}; //! skimming statistics
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "jpsiPID1", "barrel track cut"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
  Configurable<std::string> fConfigAddMCTruthHistogram{"cfgAddMCTruthHistogram", "", "Comma separated list of histograms"};
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 1.0f, "Low pt cut for tracks in the barrel"};
  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 1.0f, "Low pt cut for muons"};
  Configurable<float> fConfigMinTpcSignal{"cfgMinTpcSignal", 30.0, "Minimum TPC signal"};
  Configurable<float> fConfigMaxTpcSignal{"cfgMaxTpcSignal", 300.0, "Maximum TPC signal"};
  Configurable<std::string> fConfigMCSignals{"cfgMCsignals", "", "Comma separated list of MC signals"};
  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<bool> fConfigDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes)"};
  Configurable<bool> fIsAmbiguous{"cfgIsAmbiguous", false, "Whether we enable QA plots for ambiguous tracks"};
  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> fPropMuon{"cfgPropMuon", false, "Propgate muon tracks through absorber"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpmagPathRun2{"grpmagPathRun2", "GLO/GRP/GRP", "CCDB path of the GRPObject (Usage for Run 2)"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;

  o2::parameters::GRPObject* grpmagrun2 = nullptr; // for run 2, we access the GRPObject from GLO/GRP/GRP
  o2::parameters::GRPMagField* grpmag = nullptr;   // for run 3, we access GRPMagField from GLO/Config/GRPMagField

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fMuonCuts;  //! Muon track cuts

  bool fDoDetailedQA = false; // Bool to set detailed QA true, if QA is set true
  int fCurrentRun;            // needed to detect if the run changed and trigger update of calibrations etc.

  // TODO: filter on TPC dedx used temporarily until electron PID will be improved
  Filter barrelSelectedTracks = ifnode(fIsRun2.node() == true, aod::track::trackType == uint8_t(aod::track::Run2Track), aod::track::trackType == uint8_t(aod::track::Track)) && o2::aod::track::pt >= fConfigBarrelTrackPtLow && nabs(o2::aod::track::eta) <= 0.9f;

  Filter muonFilter = o2::aod::fwdtrack::pt >= fConfigMuonPtLow;

  void init(o2::framework::InitContext& context)
  {
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    if (fPropMuon) {
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        fCCDB->get<TGeoManager>(geoPath);
      }
    }

    // Define cuts --------------------------------------------------------------------------------------------
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

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Only use detailed QA when QA is set true
    if (fConfigQA && fConfigDetailedQA) {
      fDoDetailedQA = true;
    }

    TString histClasses = "";
    if (fDoDetailedQA) {
      histClasses += "Event_BeforeCuts;";
    }
    if (fConfigQA) {
      histClasses += "Event_AfterCuts;";
    }

    bool enableBarrelHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                               context.mOptions.get<bool>("processBarrelOnly") || context.mOptions.get<bool>("processBarrelOnlyWithDalitzBits") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithCent") || context.mOptions.get<bool>("processBarrelOnlyWithCov") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithCovWithCentAndMults") || context.mOptions.get<bool>("processBarrelOnlyWithCentAndMults") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithMults") || context.mOptions.get<bool>("processAmbiguousBarrelOnly"));
    bool enableMuonHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                             context.mOptions.get<bool>("processMuonOnlyWithCent") || context.mOptions.get<bool>("processMuonOnlyWithCov") ||
                             context.mOptions.get<bool>("processAmbiguousMuonOnlyWithCov") || context.mOptions.get<bool>("processAmbiguousMuonOnly") ||
                             context.mOptions.get<bool>("processAssociatedMuonOnly") || context.mOptions.get<bool>("processAssociatedMuonOnlyWithCov") ||
                             context.mOptions.get<bool>("processAssociatedMuonOnlyWithCovAndCent") || context.mOptions.get<bool>("processAssociatedMuonOnlyWithCovAndMults"));
    // TODO: switch on/off histogram classes depending on which process function we run
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
    }

    if (enableMuonHistos) {
      if (fDoDetailedQA) {
        histClasses += "Muons_BeforeCuts;";
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

    TString configNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> objArray(configNamesStr.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      for (int isig = 0; isig < objArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objArray->At(isig)->GetName());
        if (sig) {
          fMCSignals.push_back(*sig);
          if (fConfigQA) {
            histClasses += Form("MCTruth_%s;", objArray->At(isig)->GetName());
          }
        } else {
          continue;
        }
        if (fDoDetailedQA) {
          if (enableBarrelHistos) {
            for (auto& cut : fTrackCuts) {
              histClasses += Form("TrackBarrel_%s_%s;", cut.GetName(), objArray->At(isig)->GetName());
            }
          }
          if (enableMuonHistos) {
            for (auto& cut : fMuonCuts) {
              histClasses += Form("Muons_%s_%s;", cut.GetName(), objArray->At(isig)->GetName());
            }
          }
        }
      }
    }

    DefineHistograms(histClasses);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  Preslice<aod::McParticles_001> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<MyBarrelTracks> perCollisionTracks = aod::track::collisionId;
  Preslice<MyMuons> perCollisionMuons = aod::fwdtrack::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, uint32_t TMFTFillMap = 0u, typename TEvent, typename TTracks, typename TMuons, typename TAmbiTracks, typename TAmbiMuons, typename TMFTTracks = std::nullptr_t>
  void fullSkimming(TEvent const& collisions, aod::BCsWithTimestamps const& /*bcs*/, TTracks const& tracksBarrel, TMuons const& tracksMuon,
                    aod::McCollisions const& /*mcEvents*/, aod::McParticles_001 const& mcTracks, TAmbiTracks const& ambiTracksMid, TAmbiMuons const& ambiTracksFwd, TMFTTracks const& mftTracks = nullptr)
  {
    // Loop over collisions and produce skimmed data tables for:
    // 1) all selected collisions
    // 2) all MC collisions pointed to by selected collisions
    // 2.1) MC collision labels
    // 3) all selected tracks
    // 4) MC track labels
    //   Maps of indices are stored for all selected MC tracks (selected MC signals + MC truth particles for selected tracks)
    //   The MC particles table is generated in a separate loop over McParticles, after the indices for all MC particles we want
    //     to store are fed into the ordered std maps

    // temporary variables used for the indexing of the skimmed MC stack
    std::map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    std::map<uint64_t, uint16_t> fMCFlags;
    std::map<uint64_t, int> fEventIdx;
    std::map<uint64_t, int> fEventLabels;
    int fCounters[2] = {0, 0}; //! [0] - particle counter, [1] - event counter

    uint16_t mcflags = 0;
    uint64_t trackFilteringTag = 0;
    uint8_t trackTempFilterMap = 0;
    // Process orphan tracks
    if constexpr ((TTrackFillMap & VarManager::ObjTypes::AmbiTrack) > 0) {
      if (fDoDetailedQA && fIsAmbiguous) {
        for (auto& ambiTrack : ambiTracksMid) {
          auto trk = ambiTrack.template track_as<TTracks>();
          if (trk.collisionId() < 0) {
            VarManager::FillTrack<TTrackFillMap>(trk);
            fHistMan->FillHistClass("Orphan_TrackBarrel", VarManager::fgValues);
          }
        }
      }
    }

    if constexpr ((TMuonFillMap & VarManager::ObjTypes::AmbiMuon) > 0) {
      if (fDoDetailedQA && fIsAmbiguous) {
        for (auto& ambiTrackFwd : ambiTracksFwd) {
          auto muon = ambiTrackFwd.template fwdtrack_as<TMuons>();
          if (muon.collisionId() < 0) {
            VarManager::FillTrack<TMuonFillMap>(muon);
            if ((static_cast<int>(muon.trackType()) == 0)) {
              fHistMan->FillHistClass("Orphan_Muons_MFTMCHMID", VarManager::fgValues);
            } else if ((static_cast<int>(muon.trackType()) == 3)) {
              fHistMan->FillHistClass("Orphan_Muons_MCHMID", VarManager::fgValues);
            }
          }
        }
      }
    }

    for (auto& collision : collisions) {
      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (fCurrentRun != bc.runNumber()) {
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
        }
        fCurrentRun = bc.runNumber();
      }

      // get the trigger aliases
      uint32_t triggerAliases = collision.alias_raw();
      // store the selection decisions
      uint64_t tag = collision.selection_raw();
      if (collision.sel7()) {
        tag |= (static_cast<uint64_t>(1) << evsel::kNsel); //! SEL7 stored at position kNsel in the tag bit map
      }

      auto mcCollision = collision.mcCollision();
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::fgValues[VarManager::kRunNo] = bc.runNumber();
      VarManager::fgValues[VarManager::kBC] = bc.globalBC();
      VarManager::fgValues[VarManager::kTimestamp] = bc.timestamp();
      VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
      VarManager::FillEvent<gkEventMCFillMap>(mcCollision);

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }
      // fill stats information, before selections
      for (int i = 0; i < kNaliases; i++) {
        if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(kNaliases));

      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }

      // fill stats information, after selections
      for (int i = 0; i < kNaliases; i++) {
        if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(kNaliases));

      event(tag, bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0 && (TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
        eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                      collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                      collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                      collision.centFT0C());
      } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0) {
        eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                      collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                      collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                      -1);
      } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
        eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO],
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, collision.centFT0C());
      } else {
        eventExtended(bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), VarManager::fgValues[VarManager::kCentVZERO], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
      }
      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
      eventInfo(collision.globalIndex());
      // make an entry for this MC event only if it was not already added to the table
      if (!(fEventLabels.find(mcCollision.globalIndex()) != fEventLabels.end())) {
        eventMC(mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
                mcCollision.t(), mcCollision.weight(), mcCollision.impactParameter());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
      }
      eventMClabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMultExtra) > 0) {
        multPV(collision.multNTracksHasITS(), collision.multNTracksHasTPC(), collision.multNTracksHasTOF(), collision.multNTracksHasTRD(),
               collision.multNTracksITSOnly(), collision.multNTracksTPCOnly(), collision.multNTracksITSTPC(),
               collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
        multAll(collision.multAllTracksTPCOnly(), collision.multAllTracksITSTPC(),
                0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, 0.0);
      }

      // loop over the MC truth tracks and find those that need to be written
      auto groupedMcTracks = mcTracks.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (auto& mctrack : groupedMcTracks) {
        // check all the requested MC signals and fill a decision bit map
        mcflags = 0;
        int i = 0;
        for (auto& sig : fMCSignals) {
          bool checked = false;
          if constexpr (soa::is_soa_filtered_v<aod::McParticles_001>) {
            auto mctrack_raw = groupedMcTracks.rawIteratorAt(mctrack.globalIndex());
            checked = sig.CheckSignal(true, mctrack_raw);
          } else {
            checked = sig.CheckSignal(true, mctrack);
          }
          if (checked) {
            mcflags |= (static_cast<uint16_t>(1) << i);
          }
          i++;
        }
        if (mcflags == 0) {
          continue;
        }

        if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mctrack.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
          fMCFlags[mctrack.globalIndex()] = mcflags;
          fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
          fCounters[0]++;

          // if any of the MC signals was matched, then fill histograms and write that MC particle into the new stack
          // fill histograms for each of the signals, if found
          if (fConfigQA) {
            VarManager::FillTrackMC(mcTracks, mctrack);
            int j = 0;
            for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, j++) {
              if (mcflags & (static_cast<uint16_t>(1) << j)) {
                fHistMan->FillHistClass(Form("MCTruth_%s", (*signal).GetName()), VarManager::fgValues);
              }
            }
          }
        }
      } // end loop over mc stack

      int isAmbiguous = 0;
      // loop over reconstructed tracks
      if constexpr (static_cast<bool>(TTrackFillMap)) {
        trackBarrelInfo.reserve(tracksBarrel.size());
        trackBasic.reserve(tracksBarrel.size());
        trackBarrel.reserve(tracksBarrel.size());
        trackBarrelPID.reserve(tracksBarrel.size());
        trackBarrelLabels.reserve(tracksBarrel.size());
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
          trackBarrelCov.reserve(tracksBarrel.size());
        }

        auto groupedTracks = tracksBarrel.sliceBy(perCollisionTracks, collision.globalIndex());
        // loop over tracks
        for (auto& track : groupedTracks) {
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
          trackTempFilterMap = static_cast<uint8_t>(0);
          VarManager::FillTrack<TTrackFillMap>(track);
          // If no MC particle is found, skip the track
          if (!track.has_mcParticle()) {
            continue;
          }
          auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrackMC(mcTracks, mctrack);

          if (fDoDetailedQA) {
            fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
            if (fIsAmbiguous && isAmbiguous == 1) {
              fHistMan->FillHistClass("Ambiguous_TrackBarrel_BeforeCuts", VarManager::fgValues);
            }
          }
          // apply track cuts and fill stats histogram
          int i = 0;
          for (auto& cut : fTrackCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              if (fConfigQA) {
                fHistMan->FillHistClass(Form("TrackBarrel_%s", cut.GetName()), VarManager::fgValues); // fill the reconstructed truth
                if (fIsAmbiguous && isAmbiguous == 1) {
                  fHistMan->FillHistClass(Form("Ambiguous_TrackBarrel_%s", cut.GetName()), VarManager::fgValues);
                }
              }
              (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
            }
            i++;
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
              if (track.pidbit() & (static_cast<uint8_t>(1) << iv0)) {
                (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(fTrackCuts.size() + static_cast<float>(iv0));
              }
            }
          }
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
            trackFilteringTag |= (static_cast<uint64_t>(track.dalitzBits()) << 7); // BIT7-14: Dalitz
          }
          trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << 15); // BIT15-...:  user track filters

          mcflags = 0;
          i = 0;     // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mctrack)) {
              mcflags |= (static_cast<uint16_t>(1) << i);
              if (fDoDetailedQA) {
                j = 0;
                for (auto& cut : fTrackCuts) {
                  if (trackTempFilterMap & (uint8_t(1) << j)) {
                    fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", cut.GetName(), sig.GetName()), VarManager::fgValues); // fill the reconstructed truth
                  }
                  j++;
                }
              }
            }
            i++;
          }

          // if the MC truth particle corresponding to this reconstructed track is not already written,
          //   add it to the skimmed stack
          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

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
          trackBarrelPID(track.tpcSignal(),
                         track.tpcNSigmaEl(), track.tpcNSigmaMu(),
                         track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(),
                         track.tofNSigmaEl(), track.tofNSigmaMu(),
                         track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.trdSignal());
          trackBarrelLabels(fNewLabels.find(mctrack.index())->second, track.mcMask(), mcflags);
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
            trackBarrelCov(track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                           track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                           track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
          }
        } // end loop over reconstructed tracks
      }   // end if constexpr (static_cast<bool>(TTrackFillMap))

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

          // Check MC matching
          if (!mft.has_mcParticle()) {
            continue;
          }
          auto mctrack = mft.template mcParticle_as<aod::McParticles_001>();

          mcflags = 0;
          int i = 0; // runs over the MC signals
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mctrack)) {
              mcflags |= (static_cast<uint16_t>(1) << i);
            }
            i++;
          }

          // if the MC truth particle corresponding to this reconstructed track is not already written,
          //   add it to the skimmed stack
          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
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

          trackMFT(event.lastIndex(), trackFilteringTag, mft.pt(), mft.eta(), mft.phi());
          trackMFTExtra(mft.mftClusterSizesAndTrackFlags(), mft.sign(), dcaX, dcaY, mft.nClusters());
          trackMFTLabels(fNewLabels.find(mctrack.index())->second, mft.mcMask(), mcflags);
        } // end of mft : mftTracks

      } // end if constexpr (TMFTFillMap)

      if constexpr (static_cast<bool>(TMuonFillMap)) {
        // build the muon tables
        muonBasic.reserve(tracksMuon.size());
        muonExtra.reserve(tracksMuon.size());
        if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
          muonCov.reserve(tracksMuon.size());
        }
        muonLabels.reserve(tracksMuon.size()); // TODO: enable this once we have fwdtrack labels

        auto groupedMuons = tracksMuon.sliceBy(perCollisionMuons, collision.globalIndex());
        // loop over muons

        // first we need to get the correct indices
        int nDel = 0;
        int idxPrev = -1;
        std::map<int, int> newEntryNb;
        std::map<int, int> newMatchIndex;
        std::map<int, int> newMFTMatchIndex;

        for (auto& muon : groupedMuons) {
          trackFilteringTag = static_cast<uint64_t>(0);
          trackTempFilterMap = static_cast<uint8_t>(0);

          if (!muon.has_mcParticle()) {
            continue;
          }

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
          for (auto& cut : fMuonCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              if (fConfigQA) {
                fHistMan->FillHistClass(Form("Muons_%s", cut.GetName()), VarManager::fgValues);
              }
              (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
            }
            i++;
          }
          if (!trackTempFilterMap) { // does not pass the cuts
            nDel++;
          } else { // it passes the cuts and will be saved in the tables
            newEntryNb[muon.index()] = muon.index() - nDel;
          }
        }

        // now let's save the muons with the correct indices and matches
        for (auto& muon : groupedMuons) {
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

          trackFilteringTag = static_cast<uint64_t>(0);
          trackTempFilterMap = static_cast<uint8_t>(0);

          if (!muon.has_mcParticle()) {
            continue;
          }
          auto mctrack = muon.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrack<TMuonFillMap>(muon);
          if (fPropMuon) {
            VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
          }
          VarManager::FillTrackMC(mcTracks, mctrack);

          if (fDoDetailedQA) {
            fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
            if (fIsAmbiguous && isAmbiguous == 1) {
              fHistMan->FillHistClass("Ambiguous_Muons_BeforeCuts", VarManager::fgValues);
            }
          }
          // apply the muon selection cuts and fill the stats histogram
          int i = 0;
          for (auto& cut : fMuonCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              fHistMan->FillHistClass(Form("Muons_%s", cut.GetName()), VarManager::fgValues);
              if (fIsAmbiguous && isAmbiguous == 1) {
                fHistMan->FillHistClass(Form("Ambiguous_Muons_%s", cut.GetName()), VarManager::fgValues);
              }
              (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
            }
            i++;
          }
          if (!trackTempFilterMap) {
            continue;
          }
          // store the cut decisions
          trackFilteringTag |= static_cast<uint64_t>(trackTempFilterMap); // BIT0-7:  user selection cuts

          mcflags = 0;
          i = 0;     // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mctrack)) {
              mcflags |= (static_cast<uint16_t>(1) << i);
              if (fDoDetailedQA) {
                for (auto& cut : fMuonCuts) {
                  if (trackTempFilterMap & (uint8_t(1) << j)) {
                    fHistMan->FillHistClass(Form("Muons_%s_%s", cut.GetName(), sig.GetName()), VarManager::fgValues); // fill the reconstructed truth
                  }
                  j++;
                }
              }
            }
            i++;
          }

          // if the MC truth particle corresponding to this reconstructed track is not already written,
          //   add it to the skimmed stack
          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

          // update the matching MCH/MFT index

          if (static_cast<int>(muon.trackType()) == 0 || static_cast<int>(muon.trackType()) == 2) { // MCH-MFT or GLB track
            int matchIdx = muon.matchMCHTrackId() - muon.offsets();
            int matchMFTIdx = muon.matchMFTTrackId() - mftOffsets[muon.matchMFTTrackId()];

            // first for MCH matching index
            if (newEntryNb.count(matchIdx) > 0) {                                                  // if the key exists which means the match will not get deleted
              newMatchIndex[muon.index()] = newEntryNb[matchIdx];                                  // update the match for this muon to the updated entry of the match
              newMatchIndex[muon.index()] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()]; // adding the offset of muons, muonBasic.lastIndex() start at -1

              if (static_cast<int>(muon.trackType()) == 0) {                                       // for now only do this to global tracks
                newMatchIndex[matchIdx] = newEntryNb[muon.index()];                                // add the  updated index of this muon as a match to mch track
                newMatchIndex[matchIdx] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()];   // adding the offset, muonBasic.lastIndex() start at -1
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
            // in this case the matches should be filled from the other types but we need to check
            if (newMatchIndex.count(muon.index()) == 0) {
              newMatchIndex[muon.index()] = -1;
            }
          }

          muonBasic(event.lastIndex(), newMatchIndex.find(muon.index())->second, newMFTMatchIndex.find(muon.index())->second, trackFilteringTag, VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], muon.sign(), isAmbiguous);
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
          muonLabels(fNewLabels.find(mctrack.index())->second, muon.mcMask(), mcflags);
        }
      } // end if constexpr (static_cast<bool>(TMuonFillMap))
    }   // end loop over collisions

    // Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fNewLabelsReversed) {
      auto mctrack = mcTracks.iteratorAt(oldLabel);
      uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mctrack.has_mothers()) {
        for (auto& m : mctrack.mothersIds()) {
          if (m < mcTracks.size()) { // protect against bad mother indices
            if (fNewLabels.find(m) != fNewLabels.end()) {
              mothers.push_back(fNewLabels.find(m)->second);
            }
          } else {
            cout << "Mother label (" << m << ") exceeds the McParticles size (" << mcTracks.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }

      // TODO: Check that the daughter slice in the skimmed table works as expected
      //       Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcTracks.size()) { // protect against bad daughter indices
            if (fNewLabels.find(d) != fNewLabels.end()) {
              daughters.push_back(fNewLabels.find(d)->second);
            }
          } else {
            cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }
      int daughterRange[2] = {-1, -1};
      if (daughters.size() > 0) {
        daughterRange[0] = daughters[0];
        daughterRange[1] = daughters[daughters.size() - 1];
      }

      trackMC(fEventIdx.find(oldLabel)->second, mctrack.pdgCode(), mctrack.statusCode(), mctrack.flags(),
              mothers, daughterRange,
              mctrack.weight(), mctrack.pt(), mctrack.eta(), mctrack.phi(), mctrack.e(),
              mctrack.vx(), mctrack.vy(), mctrack.vz(), mctrack.vt(), mcflags);
      for (unsigned int isig = 0; isig < fMCSignals.size(); isig++) {
        if (mcflags & (static_cast<uint16_t>(1) << isig)) {
          (reinterpret_cast<TH1I*>(fStatsList->At(3)))->Fill(static_cast<float>(isig));
        }
      }
      if (mcflags == 0) {
        (reinterpret_cast<TH1I*>(fStatsList->At(3)))->Fill(static_cast<float>(fMCSignals.size()));
      }
    } // end loop over labels

    fCounters[0] = 0;
    fCounters[1] = 0;
    fNewLabels.clear();
    fNewLabelsReversed.clear();
    fMCFlags.clear();
    fEventIdx.clear();
    fEventLabels.clear();
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons, typename AssocTracks, typename AssocMuons>
  void fullSkimmingIndices(TEvent const& collisions, aod::BCsWithTimestamps const&, TTracks const& tracksBarrel, TMuons const& tracksMuon, aod::McCollisions const& /*mcEvents*/, aod::McParticles_001 const& mcTracks, AssocTracks const& trackIndices, AssocMuons const& fwdtrackIndices)
  {
    // Loop over collisions and produce skimmed data tables for:
    // 1) all selected collisions
    // 2) all MC collisions pointed to by selected collisions
    // 2.1) MC collision labels
    // 3) all selected tracks
    // 4) MC track labels
    //   Maps of indices are stored for all selected MC tracks (selected MC signals + MC truth particles for selected tracks)
    //   The MC particles table is generated in a separate loop over McParticles, after the indices for all MC particles we want
    //     to store are fed into the ordered std maps

    // temporary variables used for the indexing of the skimmed MC stack
    std::map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    std::map<uint64_t, uint16_t> fMCFlags;
    std::map<uint64_t, int> fEventIdx;
    std::map<uint64_t, int> fEventLabels;
    int fCounters[2] = {0, 0}; //! [0] - particle counter, [1] - event counter

    uint16_t mcflags = 0;
    uint64_t trackFilteringTag = 0;
    uint8_t trackTempFilterMap = 0;

    for (auto& collision : collisions) {
      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (fCurrentRun != bc.runNumber()) {
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
        }
        fCurrentRun = bc.runNumber();
      }

      // get the trigger aliases
      uint32_t triggerAliases = collision.alias_raw();
      // store the selection decisions
      uint64_t tag = collision.selection_raw();
      if (collision.sel7()) {
        tag |= (static_cast<uint64_t>(1) << evsel::kNsel); //! SEL7 stored at position kNsel in the tag bit map
      }

      auto mcCollision = collision.mcCollision();
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::fgValues[VarManager::kRunNo] = bc.runNumber();
      VarManager::fgValues[VarManager::kBC] = bc.globalBC();
      VarManager::fgValues[VarManager::kTimestamp] = bc.timestamp();
      VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
      VarManager::FillEvent<gkEventMCFillMap>(mcCollision);

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }
      // fill stats information, before selections
      for (int i = 0; i < kNaliases; i++) {
        if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(2.0, static_cast<float>(kNaliases));

      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }

      // fill stats information, after selections
      for (int i = 0; i < kNaliases; i++) {
        if (triggerAliases & (static_cast<uint32_t>(1) << i)) {
          (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(i));
        }
      }
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(kNaliases));

      event(tag, bc.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
      if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0 && (TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
        eventExtended(collision.bc().globalBC(), collision.bc().triggerMask(), 0, triggerAliases, VarManager::fgValues[VarManager::kCentVZERO],
                      collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                      collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                      collision.centFT0C());
      } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionMult) > 0) {
        eventExtended(bc.globalBC(), bc.triggerMask(), bc.timestamp(), triggerAliases, VarManager::fgValues[VarManager::kCentVZERO],
                      collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
                      collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(),
                      -1);
      } else if constexpr ((TEventFillMap & VarManager::ObjTypes::CollisionCent) > 0) {
        eventExtended(bc.globalBC(), bc.triggerMask(), bc.timestamp(), triggerAliases, VarManager::fgValues[VarManager::kCentVZERO],
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, collision.centFT0C());
      } else {
        eventExtended(bc.globalBC(), bc.triggerMask(), bc.timestamp(), triggerAliases, VarManager::fgValues[VarManager::kCentVZERO], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
      }
      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
      eventInfo(collision.globalIndex());
      // make an entry for this MC event only if it was not already added to the table
      if (!(fEventLabels.find(mcCollision.globalIndex()) != fEventLabels.end())) {
        eventMC(mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
                mcCollision.t(), mcCollision.weight(), mcCollision.impactParameter());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
      }
      eventMClabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());

      // loop over the MC truth tracks and find those that need to be written
      auto groupedMcTracks = mcTracks.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (auto& mctrack : groupedMcTracks) {
        // check all the requested MC signals and fill a decision bit map
        mcflags = 0;
        int i = 0;
        for (auto& sig : fMCSignals) {
          bool checked = false;
          if constexpr (soa::is_soa_filtered_v<aod::McParticles_001>) {
            auto mctrack_raw = groupedMcTracks.rawIteratorAt(mctrack.globalIndex());
            checked = sig.CheckSignal(true, mctrack_raw);
          } else {
            checked = sig.CheckSignal(true, mctrack);
          }
          if (checked) {
            mcflags |= (static_cast<uint16_t>(1) << i);
          }
          i++;
        }
        if (mcflags == 0) {
          continue;
        }

        if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mctrack.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
          fMCFlags[mctrack.globalIndex()] = mcflags;
          fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
          fCounters[0]++;

          // if any of the MC signals was matched, then fill histograms and write that MC particle into the new stack
          // fill histograms for each of the signals, if found
          if (fConfigQA) {
            VarManager::FillTrackMC(mcTracks, mctrack);
            int j = 0;
            for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, j++) {
              if (mcflags & (static_cast<uint16_t>(1) << j)) {
                fHistMan->FillHistClass(Form("MCTruth_%s", (*signal).GetName()), VarManager::fgValues);
              }
            }
          }
        }
      } // end loop over mc stack

      int isAmbiguous = 0;
      // loop over reconstructed tracks
      if constexpr (static_cast<bool>(TTrackFillMap)) {
        trackBarrelInfo.reserve(tracksBarrel.size());
        trackBasic.reserve(tracksBarrel.size());
        trackBarrel.reserve(tracksBarrel.size());
        trackBarrelPID.reserve(tracksBarrel.size());
        trackBarrelLabels.reserve(tracksBarrel.size());
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
          trackBarrelCov.reserve(tracksBarrel.size());
        }

        auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
        // loop over tracks
        for (auto& trackId : trackIdsThisCollision) {
          auto track = trackId.template track_as<TTracks>();
          if constexpr ((TTrackFillMap & VarManager::ObjTypes::AmbiTrack) > 0) {
            if (fIsAmbiguous) {
              isAmbiguous = (track.compatibleCollIds().size() != 1);
            }
          }
          trackFilteringTag = static_cast<uint64_t>(0);
          trackTempFilterMap = static_cast<uint8_t>(0);
          VarManager::FillTrack<TTrackFillMap>(track);
          // If no MC particle is found, skip the track
          if (!track.has_mcParticle()) {
            continue;
          }
          auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrackMC(mcTracks, mctrack);

          if (fDoDetailedQA) {
            fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
            if (fIsAmbiguous && isAmbiguous == 1) {
              fHistMan->FillHistClass("Ambiguous_TrackBarrel_BeforeCuts", VarManager::fgValues);
            }
          }
          // apply track cuts and fill stats histogram
          int i = 0;
          for (auto& cut : fTrackCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              if (fConfigQA) {
                fHistMan->FillHistClass(Form("TrackBarrel_%s", cut.GetName()), VarManager::fgValues); // fill the reconstructed truth
                if (fIsAmbiguous && isAmbiguous == 1) {
                  fHistMan->FillHistClass(Form("Ambiguous_TrackBarrel_%s", cut.GetName()), VarManager::fgValues);
                }
              }
              (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
            }
            i++;
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
              if (track.pidbit() & (static_cast<uint8_t>(1) << iv0)) {
                (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(fTrackCuts.size() + static_cast<float>(iv0));
              }
            }
          }
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
            trackFilteringTag |= (static_cast<uint64_t>(track.dalitzBits()) << 7); // BIT7-14: Dalitz
          }
          trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << 15); // BIT15-...:  user track filters

          mcflags = 0;
          i = 0;     // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mctrack)) {
              mcflags |= (static_cast<uint16_t>(1) << i);
              if (fDoDetailedQA) {
                j = 0;
                for (auto& cut : fTrackCuts) {
                  if (trackTempFilterMap & (uint8_t(1) << j)) {
                    fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", cut.GetName(), sig.GetName()), VarManager::fgValues); // fill the reconstructed truth
                  }
                  j++;
                }
              }
            }
            i++;
          }

          // if the MC truth particle corresponding to this reconstructed track is not already written,
          //   add it to the skimmed stack
          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

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
          trackBarrelPID(track.tpcSignal(),
                         track.tpcNSigmaEl(), track.tpcNSigmaMu(),
                         track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(),
                         track.tofNSigmaEl(), track.tofNSigmaMu(),
                         track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.trdSignal());
          trackBarrelLabels(fNewLabels.find(mctrack.index())->second, track.mcMask(), mcflags);
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
            trackBarrelCov(track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                           track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                           track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
          }
        } // end loop over reconstructed tracks
      }   // end if constexpr (static_cast<bool>(TTrackFillMap))

      if constexpr (static_cast<bool>(TMuonFillMap)) {
        // build the muon tables
        muonBasic.reserve(tracksMuon.size());
        muonExtra.reserve(tracksMuon.size());
        if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
          muonCov.reserve(tracksMuon.size());
        }
        muonLabels.reserve(tracksMuon.size()); // TODO: enable this once we have fwdtrack labels

        auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
        // loop over muons

        // first we need to get the correct indices
        int nDel = 0;
        int idxPrev = -1;
        std::map<int, int> newEntryNb;
        std::map<int, int> newMatchIndex;
        std::map<int, int> newMFTMatchIndex;

        for (auto& muonId : fwdtrackIdsThisCollision) {
          auto muon = muonId.template fwdtrack_as<TMuons>();
          trackFilteringTag = static_cast<uint64_t>(0);
          trackTempFilterMap = static_cast<uint8_t>(0);

          if (!muon.has_mcParticle()) {
            continue;
          }

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
          for (auto& cut : fMuonCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              if (fConfigQA) {
                fHistMan->FillHistClass(Form("Muons_%s", cut.GetName()), VarManager::fgValues);
              }
              (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
            }
            i++;
          }
          if (!trackTempFilterMap) { // does not pass the cuts
            nDel++;
          } else { // it passes the cuts and will be saved in the tables
            newEntryNb[muon.index()] = muon.index() - nDel;
          }
        }

        // now let's save the muons with the correct indices and matches
        for (auto& muonId : fwdtrackIdsThisCollision) {
          auto muon = muonId.template fwdtrack_as<TMuons>();
          if constexpr ((TMuonFillMap & VarManager::ObjTypes::AmbiMuon) > 0) {
            if (fIsAmbiguous) {
              isAmbiguous = (muon.compatibleCollIds().size() != 1);
            }
          }

          trackFilteringTag = static_cast<uint64_t>(0);
          trackTempFilterMap = static_cast<uint8_t>(0);

          if (!muon.has_mcParticle()) {
            continue;
          }
          auto mctrack = muon.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrack<TMuonFillMap>(muon);
          // recalculte DCA and pDca for global muon tracks
          VarManager::FillTrackCollision<TMuonFillMap>(muon, collision);
          if (fPropMuon) {
            VarManager::FillPropagateMuon<TMuonFillMap>(muon, collision);
          }
          VarManager::FillTrackMC(mcTracks, mctrack);

          if (fDoDetailedQA) {
            fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
            if (fIsAmbiguous && isAmbiguous == 1) {
              fHistMan->FillHistClass("Ambiguous_Muons_BeforeCuts", VarManager::fgValues);
            }
          }
          // apply the muon selection cuts and fill the stats histogram
          int i = 0;
          for (auto& cut : fMuonCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              fHistMan->FillHistClass(Form("Muons_%s", cut.GetName()), VarManager::fgValues);
              if (fIsAmbiguous && isAmbiguous == 1) {
                fHistMan->FillHistClass(Form("Ambiguous_Muons_%s", cut.GetName()), VarManager::fgValues);
              }
              (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(i));
            }
            i++;
          }
          if (!trackTempFilterMap) {
            continue;
          }
          // store the cut decisions
          trackFilteringTag |= static_cast<uint64_t>(trackTempFilterMap); // BIT0-7:  user selection cuts

          mcflags = 0;
          i = 0;     // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mctrack)) {
              mcflags |= (static_cast<uint16_t>(1) << i);
              if (fDoDetailedQA) {
                for (auto& cut : fMuonCuts) {
                  if (trackTempFilterMap & (uint8_t(1) << j)) {
                    fHistMan->FillHistClass(Form("Muons_%s_%s", cut.GetName(), sig.GetName()), VarManager::fgValues); // fill the reconstructed truth
                  }
                  j++;
                }
              }
            }
            i++;
          }

          // if the MC truth particle corresponding to this reconstructed track is not already written,
          //   add it to the skimmed stack
          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

          // update the matching MCH/MFT index

          if (static_cast<int>(muon.trackType()) == 0 || static_cast<int>(muon.trackType()) == 2) { // MCH-MFT or GLB track
            int matchIdx = muon.matchMCHTrackId() - muon.offsets();
            if (newEntryNb.count(matchIdx) > 0) {                                                  // if the key exists which means the match will not get deleted
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
            if (newMatchIndex.count(muon.index()) == 0) {
              newMatchIndex[muon.index()] = -1;
            }
          }

          muonBasic(event.lastIndex(), newMatchIndex.find(muon.index())->second, -1, trackFilteringTag, VarManager::fgValues[VarManager::kPt], VarManager::fgValues[VarManager::kEta], VarManager::fgValues[VarManager::kPhi], muon.sign(), isAmbiguous);
          if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
            if (fPropMuon) {
              muonExtra(muon.nClusters(), VarManager::fgValues[VarManager::kMuonPDca], VarManager::fgValues[VarManager::kMuonRAtAbsorberEnd],
                        muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                        muon.matchScoreMCHMFT(), muon.mchBitMap(), muon.midBitMap(),
                        muon.midBoards(), muon.trackType(), VarManager::fgValues[VarManager::kMuonDCAx], VarManager::fgValues[VarManager::kMuonDCAy],
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
          muonLabels(fNewLabels.find(mctrack.index())->second, muon.mcMask(), mcflags);
        }
      } // end if constexpr (static_cast<bool>(TMuonFillMap))
    }   // end loop over collisions

    // Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fNewLabelsReversed) {
      auto mctrack = mcTracks.iteratorAt(oldLabel);
      uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mctrack.has_mothers()) {
        for (auto& m : mctrack.mothersIds()) {
          if (m < mcTracks.size()) { // protect against bad mother indices
            if (fNewLabels.find(m) != fNewLabels.end()) {
              mothers.push_back(fNewLabels.find(m)->second);
            }
          } else {
            cout << "Mother label (" << m << ") exceeds the McParticles size (" << mcTracks.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }

      // TODO: Check that the daughter slice in the skimmed table works as expected
      //       Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcTracks.size()) { // protect against bad daughter indices
            if (fNewLabels.find(d) != fNewLabels.end()) {
              daughters.push_back(fNewLabels.find(d)->second);
            }
          } else {
            cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }
      int daughterRange[2] = {-1, -1};
      if (daughters.size() > 0) {
        daughterRange[0] = daughters[0];
        daughterRange[1] = daughters[daughters.size() - 1];
      }

      trackMC(fEventIdx.find(oldLabel)->second, mctrack.pdgCode(), mctrack.statusCode(), mctrack.flags(),
              mothers, daughterRange,
              mctrack.weight(), mctrack.pt(), mctrack.eta(), mctrack.phi(), mctrack.e(),
              mctrack.vx(), mctrack.vy(), mctrack.vz(), mctrack.vt(), mcflags);
      for (unsigned int isig = 0; isig < fMCSignals.size(); isig++) {
        if (mcflags & (static_cast<uint16_t>(1) << isig)) {
          (reinterpret_cast<TH1I*>(fStatsList->At(3)))->Fill(static_cast<float>(isig));
        }
      }
      if (mcflags == 0) {
        (reinterpret_cast<TH1I*>(fStatsList->At(3)))->Fill(static_cast<float>(fMCSignals.size()));
      }
    } // end loop over labels

    fCounters[0] = 0;
    fCounters[1] = 0;
    fNewLabels.clear();
    fNewLabelsReversed.clear();
    fMCFlags.clear();
    fEventIdx.clear();
    fEventLabels.clear();
  }

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      if (fConfigQA) {
        fHistMan->AddHistClass(classStr.Data());
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

      TString histMCTruthName = fConfigAddMCTruthHistogram.value;
      if (classStr.Contains("MCTruth")) {
        if (fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "mctruth_track", histMCTruthName);
        }
      }
    }

    // create statistics histograms (event, tracks, muons, MCsignals)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, (float)kNaliases + 1, -0.5, (float)kNaliases + 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    for (int ib = 1; ib <= kNaliases; ib++) {
      histEvents->GetYaxis()->SetBinLabel(ib, aliasLabels[ib - 1].data());
    }
    histEvents->GetYaxis()->SetBinLabel(kNaliases + 1, "Total");
    fStatsList->Add(histEvents);

    // Track statistics: one bin for each track selection and 5 bins for V0 tags (gamma, K0s, Lambda, anti-Lambda, Omega)
    TH1I* histTracks = new TH1I("TrackStats", "Track statistics", fTrackCuts.size() + 5.0, -0.5, fTrackCuts.size() - 0.5 + 5.0);
    ib = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ib++) {
      histTracks->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    const char* v0TagNames[5] = {"Photon conversion", "K^{0}_{s}", "#Lambda", "#bar{#Lambda}", "#Omega"};
    for (int ib = 0; ib < 5; ib++) {
      histTracks->GetXaxis()->SetBinLabel(fTrackCuts.size() + 1 + ib, v0TagNames[ib]);
    }
    fStatsList->Add(histTracks);
    TH1I* histMuons = new TH1I("MuonStats", "Muon statistics", fMuonCuts.size(), -0.5, fMuonCuts.size() - 0.5);
    ib = 1;
    for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, ib++) {
      histMuons->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }
    fStatsList->Add(histMuons);
    TH1I* histMCsignals = new TH1I("MCsignals", "MC signals", fMCSignals.size() + 1, -0.5, fMCSignals.size() - 0.5 + 1.0);
    ib = 1;
    for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, ib++) {
      histMCsignals->GetXaxis()->SetBinLabel(ib, (*signal).GetName());
    }
    histMCsignals->GetXaxis()->SetBinLabel(fMCSignals.size() + 1, "Others (matched to reco tracks)");
    fStatsList->Add(histMCsignals);
  }

  // Produce barrel + muon tables ------------------------------------------------------------------------------------
  void processFull(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                   soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon,
                   aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, gkMuonFillMap, 0u>(collisions, bcs, tracksBarrel, tracksMuon, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  void processFullWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                          soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon,
                          aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, gkMuonFillMapWithCov, 0u>(collisions, bcs, tracksBarrel, tracksMuon, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables ------------------------------------------------------------------------------------
  void processBarrelOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                         soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                         aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with multiplicity ------------------------------------------------------------------------------------
  void processBarrelOnlyWithMults(MyEventsWithMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                  soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                                  aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithMults, gkTrackFillMap, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with centrality ------------------------------------------------------------------------------------
  void processBarrelOnlyWithCent(MyEventsWithCent const& collisions, aod::BCsWithTimestamps const& bcs,
                                 soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                                 aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMap, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with centrality and multiplicity -------------------------------------------------------------------
  void processBarrelOnlyWithCentAndMults(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                         soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                                         aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMap, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with cov matrix-----------------------------------------------------------------------
  void processBarrelOnlyWithCov(MyEventsWithMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel,
                                aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithMults, gkTrackFillMapWithCov, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with centrality, multiplicity and cov matrix -------------------------------------------------------------------
  void processBarrelOnlyWithCovWithCentAndMults(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                                soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel,
                                                aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithCentAndMults, gkTrackFillMapWithCov, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce barrel only tables, with cov matrix and dalitz bits-----------------------------------------------------------------------
  void processBarrelOnlyWithDalitzBits(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                       soa::Filtered<MyBarrelTracksWithDalitzBits> const& tracksBarrel,
                                       aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithDalitzBits, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }

  // Produce muon only tables ------------------------------------------------------------------------------------
  void processMuonOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                       soa::Filtered<MyMuons> const& tracksMuon,
                       aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMap, 0u>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }
  // Produce muon only tables, with centrality-------------------------------------------------------------------------------
  void processMuonOnlyWithCent(MyEventsWithCent const& collisions, aod::BCsWithTimestamps const& bcs,
                               soa::Filtered<MyMuons> const& tracksMuon,
                               aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithCent, 0u, gkMuonFillMap, 0u>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }
  // Produce muon only tables, with cov matrix ------------------------------------------------------------------------------------
  void processMuonOnlyWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                              soa::Filtered<MyMuonsWithCov> const& tracksMuon,
                              aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov, 0u>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, nullptr, nullptr);
  }
  // Produce MFT tracks tables and muons  ------------------------------------------------------------------------------------------------------------------
  void processMuonsAndMFT(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                          MyMftTracks const& tracksMft, MyMuonsWithCov const& tracksMuon,
                          aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov, gkMFTFillMap>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, nullptr, tracksMft);
  }
  // Produce muon tables only based on track-collision association tables --------------------------------------------------------------------------------------
  void processAssociatedMuonOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                 soa::Filtered<MyMuonsColl> const& tracksMuon,
                                 aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    fullSkimmingIndices<gkEventFillMap, 0u, gkMuonFillMapWithAmbi>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, fwdtrackIndices);
  }

  void processAssociatedMuonOnlyWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                        soa::Filtered<MyMuonsCollWithCov> const& tracksMuon,
                                        aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    fullSkimmingIndices<gkEventFillMap, 0u, gkMuonFillMapWithCovAmbi>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, fwdtrackIndices);
  }

  void processAssociatedMuonOnlyWithCovAndCent(MyEventsWithCent const& collisions, aod::BCsWithTimestamps const& bcs,
                                               soa::Filtered<MyMuonsCollWithCov> const& tracksMuon,
                                               aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    fullSkimmingIndices<gkEventFillMap, 0u, gkMuonFillMapWithCovAmbi>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, fwdtrackIndices);
  }

  void processAssociatedMuonOnlyWithCovAndMults(MyEventsWithCentAndMults const& collisions, aod::BCsWithTimestamps const& bcs,
                                                soa::Filtered<MyMuonsCollWithCov> const& tracksMuon,
                                                aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    fullSkimmingIndices<gkEventFillMap, 0u, gkMuonFillMapWithCovAmbi>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, fwdtrackIndices);
  }
  // Produce muon tables only for ambiguous tracks studies --------------------------------------------------------------------------------------
  void processAmbiguousMuonOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                soa::Filtered<MyMuons> const& tracksMuon,
                                aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::AmbiguousFwdTracks const& ambiTracksFwd)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithAmbi, 0u>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, ambiTracksFwd, nullptr);
  }

  void processAmbiguousMuonOnlyWithCov(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                       soa::Filtered<MyMuonsWithCov> const& tracksMuon,
                                       aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::AmbiguousFwdTracks const& ambiTracksFwd)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCovAmbi, 0u>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks, nullptr, ambiTracksFwd, nullptr);
  }

  // Produce track tables only for ambiguous tracks studies -------------------------------------------------------------------------------------
  void processAmbiguousBarrelOnly(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs,
                                  soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                                  aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks, aod::AmbiguousTracks const& ambiTracksMid)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithAmbi, 0u, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks, ambiTracksMid, nullptr, nullptr);
  }
  // Process the BCs and store stats for luminosity retrieval -----------------------------------------------------------------------------------
  void processOnlyBCs(soa::Join<aod::BCs, aod::BcSels>::iterator const& bc)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (bc.alias_bit(i) > 0) {
        (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(i));
      }
    }
    (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(0.0, static_cast<float>(kNaliases));
  }

  PROCESS_SWITCH(TableMakerMC, processFull, "Produce both barrel and muon skims", false);
  PROCESS_SWITCH(TableMakerMC, processFullWithCov, "Produce both barrel and muon skims, w/ track and fwdtrack cov tables", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnly, "Produce barrel skims", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithMults, "Produce barrel skims, with mults", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithCent, "Produce barrel skims, w/ centrality", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithCentAndMults, "Produce barrel skims, w/ centrality and mults", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithCov, "Produce barrel skims, with track covariance matrix", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithCovWithCentAndMults, "Produce barrel skims, w/ centrality and mults, with track covariance matrix", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithDalitzBits, "Produce barrel skims, and dalitz bits", false);
  PROCESS_SWITCH(TableMakerMC, processMuonOnly, "Produce muon skims", false);
  PROCESS_SWITCH(TableMakerMC, processMuonOnlyWithCov, "Produce muon skims, with muon covariance matrix", false);
  PROCESS_SWITCH(TableMakerMC, processMuonsAndMFT, "Produce muon and MFT skims", false);
  PROCESS_SWITCH(TableMakerMC, processMuonOnlyWithCent, "Produce muon skims, w/ centrality", false);
  PROCESS_SWITCH(TableMakerMC, processOnlyBCs, "Analyze the BCs to store sampled lumi", false);
  PROCESS_SWITCH(TableMakerMC, processAssociatedMuonOnly, "Produce muon skims using track-collision association tables", false);
  PROCESS_SWITCH(TableMakerMC, processAssociatedMuonOnlyWithCov, "Produce muon skims using track-collision association tables, with muon covariance matrix", false);
  PROCESS_SWITCH(TableMakerMC, processAssociatedMuonOnlyWithCovAndCent, "Produce muon skims using track-collision association tables, w/ centrality", false);
  PROCESS_SWITCH(TableMakerMC, processAssociatedMuonOnlyWithCovAndMults, "Produce muon skims using track-collision association tables, w/ mults", false);
  PROCESS_SWITCH(TableMakerMC, processAmbiguousMuonOnly, "Build muon-only DQ skimmed data model with QA plots for ambiguous muons", false);
  PROCESS_SWITCH(TableMakerMC, processAmbiguousMuonOnlyWithCov, "Build muon-only with cov DQ skimmed data model with QA plots for ambiguous muons", false);
  PROCESS_SWITCH(TableMakerMC, processAmbiguousBarrelOnly, "Build barrel-only DQ skimmed data model with QA plots for ambiguous tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // TODO: For now TableMakerMC works just for PbPb (cent table is present)
  //      Implement workflow arguments for pp/PbPb and possibly merge the task with tableMaker.cxx
  return WorkflowSpec{
    adaptAnalysisTask<TableMakerMC>(cfgc)};
}
