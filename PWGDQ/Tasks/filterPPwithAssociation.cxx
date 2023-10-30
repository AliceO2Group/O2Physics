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
#include <iostream>
#include <vector>
#include <memory>
#include <cstring>
#include <TH1F.h>
#include <TH2I.h>
#include <THashList.h>
#include <TString.h>
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
#include "EventFiltering/filterTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "CommonConstants/LHCConstants.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace
{
enum DQTriggers {
  kSingleE = 0,
  kDiElectron,
  kSingleMuLow,
  kSingleMuHigh,
  kDiMuon,
  kNTriggersDQ
};
} // namespace
namespace o2::aod
{
namespace dqppfilter
{
DECLARE_SOA_COLUMN(IsDQEventSelected, isDQEventSelected, int);
DECLARE_SOA_COLUMN(IsDQBarrelSelected, isDQBarrelSelected, uint32_t);
DECLARE_SOA_COLUMN(IsDQMuonSelected, isDQMuonSelected, uint32_t);
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! Collision index
DECLARE_SOA_INDEX_COLUMN(Track, track);         //! Track index
DECLARE_SOA_INDEX_COLUMN(FwdTrack, fwdtrack);   //! FwdTrack index
} // namespace dqppfilter

DECLARE_SOA_TABLE(DQTrackAssoc, "AOD", "DQTRACKASSOC", //! Table for track-to-collision association (tracks can appear for several collisions)
                  dqppfilter::CollisionId,
                  dqppfilter::TrackId,
                  dqppfilter::IsDQBarrelSelected);
DECLARE_SOA_TABLE(DQMuonAssoc, "AOD", "DQMUONASSOC", //! Table for muon-to-collision association (tracks can appear for several collisions)
                  dqppfilter::CollisionId,
                  dqppfilter::FwdTrackId,
                  dqppfilter::IsDQMuonSelected);

DECLARE_SOA_TABLE(DQEventCuts, "AOD", "DQEVENTCUTS", dqppfilter::IsDQEventSelected);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsSelected = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventCuts>;
// TODO: subscribe to the bare needed minimum, in particular for the CEFP task
// TODO: test working with TrackIU
// TODO: remove the TOF pid if not needed (requires changes in VarManager to separate TrackPID into TPC and TOF)
using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct DQEventSelectionTask {
  Produces<aod::DQEventCuts> eventSel;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan = nullptr;
  AnalysisCompositeCut* fEventCut;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Comma separated list of event cuts; multiple cuts are applied with a logical AND"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  // TODO: configure the histogram classes to be filled by QA

  void init(o2::framework::InitContext&)
  {
    // Construct the composite cut out of the user provided event selection cut(s)
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut = new AnalysisCompositeCut(true);
    if (!eventCutStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(eventCutStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fEventCut->AddCut(dqcuts::GetAnalysisCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;"); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                 // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSelection(TEvent const& collision, aod::BCs const& bcs)
  {
    // Reset the Values array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);

    VarManager::FillEvent<gkEventFillMap>(collision);
    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
    }
    if (fEventCut->IsSelected(VarManager::fgValues)) {
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }
      eventSel(1);
    } else {
      eventSel(0);
    }
  }

  void processEventSelection(MyEvents::iterator const& collision, aod::BCs const& bcs)
  {
    runEventSelection<gkEventFillMap>(collision, bcs);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQEventSelectionTask, processEventSelection, "Run event selection", false);
  PROCESS_SWITCH(DQEventSelectionTask, processDummy, "Dummy function", false);
};

struct DQBarrelTrackSelection {
  Produces<aod::DQTrackAssoc> trackAssoc;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigCuts{"cfgBarrelTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> fConfigCcdbPathTPC{"ccdb-path-tpc", "Users/i/iarsene/Calib/TPCpostCalib", "base path to the ccdb object"};
  // Configurable<std::string> rct_path{"rct-path", "RCT/Info/RunInformation", "path to the ccdb RCT objects for the SOR timestamps"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<bool> fConfigComputeTPCpostCalib{"cfgTPCpostCalib", false, "If true, compute TPC post-calibrated n-sigmas"};
  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 0.5f, "Low pt cut for tracks in the barrel"};
  Configurable<float> fConfigMinTpcSignal{"cfgMinTpcSignal", 65.0, "Minimum TPC signal"};
  Configurable<float> fConfigMaxTpcSignal{"cfgMaxTpcSignal", 110.0, "Maximum TPC signal"};
  Configurable<int> fConfigCollisionTrackAssoc{"cfgCollisionTrackAssoc", 0, "0 - standard association, 1 - time compatibility, 2 - ambiguous"};
  Configurable<float> fConfigAssocTimeMargin{"cfgAssocTimeMargin", 0.0f, "Extra time margin to be considered when doing collision - track matching (in ns)"};
  Configurable<float> fConfigSigmaForTimeCompat{"cfgSigmaForTimeCompat", 4.0, "nSigma window when doing collision - track matching "};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  // o2::ccdb::CcdbApi fCCDB_api;                /// API to access CCDB headers

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<TString> fCutHistNames;

  int fCurrentRun; // needed to detect if the run changed and trigger update of calibrations etc.
  std::map<int64_t, uint32_t> fSelectedTracks;

  // int fTimeFrameAssociation;
  // int fTimeFrameSelection;

  Filter barrelSelectedTracks = ifnode(fIsRun2.node() == true, aod::track::trackType == uint8_t(aod::track::Run2Track), aod::track::trackType == uint8_t(aod::track::Track)) && o2::aod::track::pt >= fConfigBarrelTrackPtLow &&
                                nabs(o2::aod::track::eta) <= 0.9f &&
                                o2::aod::track::tpcSignal >= fConfigMinTpcSignal && o2::aod::track::tpcSignal <= fConfigMaxTpcSignal;

  void init(o2::framework::InitContext&)
  {
    // fTimeFrameAssociation = 0;
    // fTimeFrameSelection = 0;

    if (fConfigCollisionTrackAssoc.value < 0 || fConfigCollisionTrackAssoc.value > 2) {
      LOGF(fatal, "Invalid collision-track association option. Must be either 0, 1 or 2");
    }

    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        AnalysisCompositeCut* cut = dqcuts::GetCompositeCut(objArray->At(icut)->GetName());
        if (cut) {
          fTrackCuts.push_back(*cut);
        } else {
          LOGF(fatal, "Invalid barrel track cut provided: %s", objArray->At(icut)->GetName());
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      TString cutNames = "TrackBarrel_BeforeCuts;";
      for (auto& cut : fTrackCuts) {
        cutNames += Form("TrackBarrel_%s;", cut.GetName());
        fCutHistNames.push_back(Form("TrackBarrel_%s", cut.GetName()));
      }

      DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
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
      /*fCCDB_api.init(fConfigCcdbUrl.value);
      if (!fCCDB_api.isHostReachable()) {
        LOGF(fatal, "CCDB host %s is not reacheable, cannot go forward", fConfigCcdbUrl.value.data());
      }*/
    }
    fSelectedTracks.clear();
  }

  template <typename TTracks>
  void associateTracksToCollisionsStandard(Collisions const& collisions,
                                           TTracks const& tracksBarrel)
  {
    // This is the standard association: does nothing, just associate the already selected tracks to their own default collision
    // loop over collisions to find time-compatible tracks
    if (fSelectedTracks.size() == 0) {
      return;
    }

    for (auto const& [trackIdx, filterMap] : fSelectedTracks) {
      auto track = tracksBarrel.rawIteratorAt(trackIdx);
      trackAssoc(track.collision().globalIndex(), trackIdx, filterMap);
    }
  }

  template <typename TTracks>
  void associateTracksToCollisionsTime(Collisions const& collisions,
                                       TTracks const& tracksBarrel,
                                       BCsWithTimestamps const& bcs)
  {
    // loop over collisions to find time-compatible tracks
    if (fSelectedTracks.size() == 0) {
      return;
    }
    float timeMargin = fConfigAssocTimeMargin.value;
    float nSigmaForTimeCompat = fConfigSigmaForTimeCompat.value;

    auto trackBegin = fSelectedTracks.begin();
    const auto bOffsetMax = 241; // 6 mus (ITS)
    for (const auto& collision : collisions) {
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      uint64_t collBC = bc.globalBC();
      bool iteratorMoved = false;

      for (auto trackFiltered = trackBegin; trackFiltered != fSelectedTracks.end(); trackFiltered++) {
        auto track = tracksBarrel.rawIteratorAt(trackFiltered->first);
        auto trackBC = track.collision().template bc_as<aod::BCsWithTimestamps>();
        const int64_t bcOffset = (int64_t)trackBC.globalBC() - (int64_t)collBC;

        float trackTime{0.};
        float trackTimeRes{0.};
        if (track.isPVContributor()) {
          trackTime = track.collision().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes = 25.f;                           // 1 BC
        } else {
          trackTime = track.trackTime();
          trackTimeRes = track.trackTimeRes();
        }
        const float deltaTime = trackTime - collTime + bcOffset * 25.f;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;

        // optimization to avoid looping over all collisions in the TF (loop just for +/- bOffsetMax around the track)
        if (!iteratorMoved && bcOffset > -bOffsetMax) {
          trackBegin = trackFiltered;
          iteratorMoved = true;
        } else if (bcOffset > bOffsetMax) {
          break;
        }

        float thresholdTime = 0.;
        if (track.isPVContributor()) {
          thresholdTime = trackTimeRes;
        } else if (TESTBIT(track.flags(), o2::aod::track::TrackTimeResIsRange)) {
          thresholdTime = std::sqrt(sigmaTimeRes2) + timeMargin;
        } else {
          thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;
        }

        if (std::abs(deltaTime) < thresholdTime) {
          const auto collIdx = collision.globalIndex();
          trackAssoc(collIdx, trackFiltered->first, trackFiltered->second);
        }
      } // end loop for tracks
    }   // end for collisions
  }

  template <typename TTracks>
  void associateTracksToCollisionsAmbigous(Collisions const& collisions,
                                           TTracks const& tracksBarrel,
                                           BCsWithTimestamps const& bcs,
                                           AmbiguousTracks const& ambTracks)
  {
    if (fSelectedTracks.size() == 0) {
      return;
    }

    // map to keep all collision-track associations (ordered based on the key by construction); for each collision there is a vector of tracks
    std::map<int64_t, std::vector<int64_t>> collTrackIds;
    // map to hold collision - BC associations
    std::map<uint64_t, uint64_t> collBCmap;

    // create the standard collision - track associations (track and their default collision)
    for (auto const& [trackIdx, filterMap] : fSelectedTracks) {
      auto track = tracksBarrel.rawIteratorAt(trackIdx);
      auto collId = track.collisionId();
      if (collTrackIds.find(collId) == collTrackIds.end()) { // this colision is not in the map
        std::vector<int64_t> idxs{trackIdx};
        collTrackIds[collId] = idxs;
      } else { // collision is in the map, add track to the vector
        auto idxs = collTrackIds[collId];
        // cout << "collision in the map already with " << idxs.size() << " tracks" << endl;
        idxs.push_back(trackIdx);
        collTrackIds[collId] = idxs;
      }
      auto collision = collisions.rawIteratorAt(collId);
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      collBCmap[collId] = bc.globalBC();
    }

    // associate ambiguous tracks with compatible collisions
    for (const auto& ambTrack : ambTracks) {
      // use just the ambiguous tracks that are selected
      if (fSelectedTracks.find(ambTrack.trackId()) == fSelectedTracks.end()) {
        continue;
      }
      auto track = ambTrack.template track_as<TTracks>();

      // Loop over all collisions
      for (const auto& collision : collisions) {
        const auto collId = collision.globalIndex();
        // if this collision is the same as the default one of the track, skip since this association was added already
        if (collId == track.collisionId()) {
          continue;
        }

        // get the BC of this collision
        auto bcColl = collision.template bc_as<aod::BCsWithTimestamps>();
        const uint64_t mostProbableBc = bcColl.globalBC();
        // loop over the track BC slice and check if this collision matches
        for (const auto& bc : ambTrack.bc()) {
          if (bc.globalBC() == mostProbableBc) { // found a match, add it to the association map
            const int64_t trackId = track.globalIndex();
            if (collTrackIds.find(collId) == collTrackIds.end()) {
              std::vector<int64_t> idxs{trackId};
              collTrackIds[collId] = idxs;
            } else {
              auto idxs = collTrackIds[collId];
              idxs.push_back(trackId);
              collTrackIds[collId] = idxs;
            }
            collBCmap[collId] = bc.globalBC();
          }
        } // end loop over BCs
      }   // end loop over collisions
    }     // end loop over ambiguous tracks

    // associate tracks with in-bunch pileup collisions
    // Loop over the already existing associations
    // NOTE: to be checked if this does something, since the pileup collisions should had been already associated in the previous step
    for (auto const& [collId, assocTracks] : collTrackIds) {
      const uint64_t currentBC = collBCmap[collId];
      std::vector<uint64_t> pileupCollIds{};
      for (const auto& collision : collisions) { // NOTE: since collisions should be ordered in time, one could optimize here and not make the full loop
        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
        if (bc.globalBC() == currentBC) {
          if (collision.globalIndex() != collId) { // found a pileup collision (different collisions pointing to the same BC)
            pileupCollIds.push_back(collision.globalIndex());
          }
        }
      }

      // in-bunch pileup requires at least 2 collisions in the same bunch (collId + another one)
      if (pileupCollIds.size() < 1) {
        continue;
      }

      // here we have in bunch pileup and we will associate the tracks from collId to the other (pileup) collisions
      for (const auto& pileupCollId : pileupCollIds) {
        for (const auto& assocTrack : assocTracks) {
          if (collTrackIds.find(pileupCollId) == collTrackIds.end()) { // this colision is not in the map
            std::vector<int64_t> idxs{assocTrack};
            collTrackIds[pileupCollId] = idxs;
          } else {
            auto idxs = collTrackIds[pileupCollId];
            if (std::find(idxs.begin(), idxs.end(), assocTrack) == idxs.end()) { // check if the track is already in the list or not
              idxs.push_back(assocTrack);
            }
            collTrackIds[pileupCollId] = idxs;
          }
        } // end loop over tracks associated to collId
      }   // end loop over pileup collisions
    }

    // create the collision - track association table
    for (auto const& [collId, assocTracks] : collTrackIds) {
      for (const auto& assocTrack : assocTracks) {
        trackAssoc(collId, assocTrack, fSelectedTracks[assocTrack]);
      }
    }
  }

  // Templated function instantianed for all of the process functions
  template <uint32_t TTrackFillMap, typename TTracks>
  void runTrackSelection(aod::BCsWithTimestamps const& bcs, TTracks const& tracksBarrel)
  {
    fSelectedTracks.clear();
    auto bc = bcs.begin(); // check just the first bc to get the run number
    if (fCurrentRun != bc.runNumber()) {
      fCurrentRun = bc.runNumber();
      if (fConfigComputeTPCpostCalib) {
        auto calibList = fCCDB->getForTimeStamp<TList>(fConfigCcdbPathTPC.value, bc.timestamp());
        VarManager::SetCalibrationObject(VarManager::kTPCElectronMean, calibList->FindObject("mean_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCElectronSigma, calibList->FindObject("sigma_map_electron"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionMean, calibList->FindObject("mean_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCPionSigma, calibList->FindObject("sigma_map_pion"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonMean, calibList->FindObject("mean_map_proton"));
        VarManager::SetCalibrationObject(VarManager::kTPCProtonSigma, calibList->FindObject("sigma_map_proton"));
      }

      /*std::map<std::string, std::string> metadata, headers;
      const std::string run_path = Form("%s/%i", rct_path.value.data(), fCurrentRun);
      headers = fCCDB_api.retrieveHeaders(run_path, metadata, -1);
      if (headers.count("SOR") == 0) {
        LOGF(fatal, "Cannot find start-of-run timestamp for run number in path '%s'.", run_path.data());
      }
      if (headers.count("EOR") == 0) {
        LOGF(fatal, "Cannot find end-of-run timestamp for run number in path '%s'.", run_path.data());
      }

      int64_t sorTimestamp = atol(headers["SOR"].c_str()); // timestamp of the SOR in ms
      int64_t eorTimestamp = atol(headers["EOR"].c_str()); // timestamp of the EOR in ms
      */
    }

    uint32_t filterMap = uint32_t(0);
    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);
    for (auto& track : tracksBarrel) {
      // NOTE: Here we explicitly remove orphan tracks
      if (!track.has_collision()) {
        continue;
      }
      filterMap = uint32_t(0);

      VarManager::FillTrack<TTrackFillMap>(track);
      if (fConfigQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA) {
            fHistMan->FillHistClass(fCutHistNames[i].Data(), VarManager::fgValues);
          }
        }
      }
      if (filterMap != 0) {
        fSelectedTracks[track.globalIndex()] = filterMap;
      }
    } // end loop over tracks
  }

  void processSelection(Collisions const& collisions, aod::BCsWithTimestamps const& bcs,
                        MyBarrelTracks const& tracks, soa::Filtered<MyBarrelTracks> const& filteredTracks,
                        AmbiguousTracks const& ambTracks)
  {
    runTrackSelection<gkTrackFillMap>(bcs, filteredTracks);
    if (fConfigCollisionTrackAssoc.value == 0) {
      associateTracksToCollisionsStandard(collisions, tracks);
    } else if (fConfigCollisionTrackAssoc.value == 1) {
      associateTracksToCollisionsTime(collisions, tracks, bcs);
    } else {
      associateTracksToCollisionsAmbigous(collisions, tracks, bcs, ambTracks);
    }
  }

  void processDummy(MyBarrelTracks&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQBarrelTrackSelection, processSelection, "Run barrel track selection", false);
  PROCESS_SWITCH(DQBarrelTrackSelection, processDummy, "Dummy function", false);
};

struct DQMuonsSelection {
  Produces<aod::DQMuonAssoc> muonAssoc;
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;

  HistogramRegistry registry{
    "registry",
    {{"Association/DeltaT", "; t_vtx - t_track; counts", {HistType::kTH1F, {{2000, -50., 50.}}}},
     {"Association/AssociationTrackStatus", "; status; counts", {HistType::kTH1F, {{10, 0, 10}}}}}};

  Configurable<std::string> fConfigCuts{"cfgMuonsCuts", "muonQualityCuts", "Comma separated list of ADDITIONAL muon track cuts"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};
  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 0.5f, "Low pt cut for muons"};
  Configurable<int> fConfigCollisionMuonAssoc{"cfgCollisionMuonAssoc", 0, "0 - standard association, 1 - time compatibility, 2 - ambiguous"};
  Configurable<float> fConfigAssocTimeMargin{"cfgAssocTimeMargin", 0.0f, "Extra time margin to be considered when doing collision - muon matching (in ns)"};
  Configurable<float> fConfigSigmaForTimeCompat{"cfgSigmaForTimeCompat", 4.0, "nSigma window when doing collision - track matching "};
  Configurable<float> fSigmaTrack{"cfgSigmaTrack", 1.0, "Number of sigma for track time window"};
  Configurable<float> fSigmaVtx{"cfgSigmaVtx", 4.0, "Number of sigma for vertex time window"};
  Configurable<float> fTimeMarginTrack{"cfgTimeMarginTrack", 0.0, "Number of sigma for track time window"};
  Configurable<float> fTimeMarginVtx{"cfgTimeMarginVtx", 0.0, "Number of sigma for vertex time window"};
  Configurable<float> fTimeBias{"cfgTimeBias", 0.0, "Number of sigma for track time window"};

  // TODO: configure the histogram classes to be filled by QA

  Filter muonFilter = o2::aod::fwdtrack::pt >= fConfigMuonPtLow;

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<TString> fCutHistNames;

  std::map<int64_t, uint32_t> fSelectedMuons;
  std::map<int64_t, int> isMuonReassigned;
  std::vector<std::pair<std::pair<double, double>, int>> vtxOrdBrack;

  void init(o2::framework::InitContext&)
  {
    TString cutNamesStr = fConfigCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      TString cutNames = "Muon_BeforeCuts;";
      for (unsigned int i = 0; i < fTrackCuts.size(); i++) {
        cutNames += Form("Muon_%s;", fTrackCuts[i].GetName());
        fCutHistNames.push_back(Form("Muon_%s", fTrackCuts[i].GetName()));
      }

      DefineHistograms(fHistMan, cutNames.Data());     // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    auto hstatus2 = registry.get<TH1>(HIST("Association/AssociationTrackStatus"));
    auto* x3 = hstatus2->GetXaxis();
    x3->SetBinLabel(1, "Ambiguous tracks");
    x3->SetBinLabel(2, "Orphan tracks");
    x3->SetBinLabel(3, "Associated tracks");
    x3->SetBinLabel(4, "Ambiguous after assoc.");
    x3->SetBinLabel(5, "Orphan after assoc.");
    x3->SetBinLabel(6, "Reassociated ");
  }

  template <uint32_t TMuonFillMap, typename TMuons>
  void runMuonSelection(TMuons const& muons)
  {
    fSelectedMuons.clear();

    uint32_t filterMap = uint32_t(0);
    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
    // TODO: fill event information which might be needed in histograms or cuts that combine track and event properties

    for (auto& muon : muons) {
      filterMap = uint32_t(0);
      // NOTE: here we do not exclude orphan muon tracks
      VarManager::FillTrack<TMuonFillMap>(muon);
      if (fConfigQA) {
        fHistMan->FillHistClass("Muon_BeforeCuts", VarManager::fgValues);
      }
      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
        if ((*cut).IsSelected(VarManager::fgValues)) {
          filterMap |= (uint32_t(1) << i);
          if (fConfigQA) {
            fHistMan->FillHistClass(fCutHistNames[i].Data(), VarManager::fgValues);
          }
        }
      }
      if (filterMap != 0) {
        fSelectedMuons[muon.globalIndex()] = filterMap;
      }
    } // end loop over muons
  }

  void runCollisionMap(Collisions const& collisions, aod::BCsWithTimestamps const& bcstimestamp)
  {
    // association of time brackets to each collision
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      double t0 = bc.globalBC() * o2::constants::lhc::LHCBunchSpacingNS - collision.collisionTime();
      double err = collision.collisionTimeRes() * fSigmaVtx + fTimeMarginVtx;
      std::pair<double, double> timeBracket = {t0 - err, t0 + err};
      std::pair<std::pair<double, double>, int> timeNID = {timeBracket, collision.globalIndex()};
      vtxOrdBrack.emplace_back(timeNID);
    }
    // sorting collision according to time
    std::sort(vtxOrdBrack.begin(), vtxOrdBrack.end(), [](const std::pair<std::pair<double, double>, int>& a, const std::pair<std::pair<double, double>, int>& b) { return a.first.first < b.first.first; });
  }

  template <typename TMuons>
  void associateMuonsToCollisionsStandard(Collisions const& collisions,
                                          TMuons const& muons)
  {
    // This is the standard association: does nothing, just associate the already selected muons to their own default collision
    if (fSelectedMuons.size() == 0) {
      return;
    }

    for (auto const& [muonIdx, filterMap] : fSelectedMuons) {
      auto muon = muons.rawIteratorAt(muonIdx);
      // We do not associate orphan tracks here
      if (muon.has_collision()) {
        muonAssoc(muon.collision().globalIndex(), muonIdx, filterMap);
      }
    }
  }

  template <typename TMuons>
  void associateMuonsToCollisionsTime(Collisions const& collisions,
                                      TMuons const& muons,
                                      BCs const& bcs)
  {
    // loop over collisions to find time-compatible muons
    if (fSelectedMuons.size() == 0) {
      return;
    }
    float timeMargin = fConfigAssocTimeMargin.value;
    float nSigmaForTimeCompat = fConfigSigmaForTimeCompat.value;

    auto trackBegin = fSelectedMuons.begin();
    const auto bOffsetMax = 200; // check 200 BCs in past and future

    for (const auto& collision : collisions) {
      // get this collisions time and BC
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();

      bool iteratorMoved = false;
      // loop over filtered muons
      for (auto trackFiltered = trackBegin; trackFiltered != fSelectedMuons.end(); trackFiltered++) {
        auto muon = muons.rawIteratorAt(trackFiltered->first);
        if (!muon.has_collision()) {
          continue;
        }

        const int64_t bcOffset = (int64_t)muon.collision().bc().globalBC() - (int64_t)collBC;
        float trackTime = muon.trackTime();
        float trackTimeRes = muon.trackTimeRes();
        const float deltaTime = trackTime - collTime + bcOffset * 25.f;
        // TODO: Time uncertainty is computed in quadrature, but maybe that should just be added linearly ?
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;

        // optimization to only loop over tracks which are within a given BC window
        if (!iteratorMoved && bcOffset > -bOffsetMax) {
          trackBegin = trackFiltered;
          iteratorMoved = true;
        } else if (bcOffset > bOffsetMax) {
          break;
        }

        // If the muon and collision are compatible, create the association
        float thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;
        if (std::abs(deltaTime) < thresholdTime) {
          const auto collIdx = collision.globalIndex();
          muonAssoc(collIdx, trackFiltered->first, trackFiltered->second);
        }
      } // end for tracks
    }   // end for collisions
  }

  template <typename TMuons>
  void associateMuonsToCollisionsAmbigous(Collisions const& collisions,
                                          TMuons const& muons,
                                          BCs const& bcs,
                                          AmbiguousFwdTracks const& ambMuons)
  {
    if (fSelectedMuons.size() == 0) {
      return;
    }

    std::map<int64_t, std::vector<int64_t>> collTrackIds; // map to keep all collision-track associations (ordered based on the key by construction)
    std::map<uint64_t, uint64_t> collBCmap;               // map to hold collision - BC associations

    // first lets associate all the non-orphan muons to their primary collision Id
    for (auto const& [muonIdx, filterMap] : fSelectedMuons) {
      auto muon = muons.rawIteratorAt(muonIdx);
      if (!muon.has_collision()) {
        continue;
      }
      auto collId = muon.collisionId();
      if (collTrackIds.find(collId) == collTrackIds.end()) { // this colision is not in the map
        std::vector<int64_t> idxs{muonIdx};
        collTrackIds[collId] = idxs;
      } else { // collision is in the map, add track to the associated vector
        auto idxs = collTrackIds[collId];
        idxs.push_back(muonIdx);
        collTrackIds[collId] = idxs;
      }
      collBCmap[collId] = muon.collision().bc().globalBC();
    }

    // associate collisions with ambiguous tracks
    for (const auto& ambMuon : ambMuons) {
      // consider only the filtered muons
      if (fSelectedMuons.find(ambMuon.fwdtrackId()) == fSelectedMuons.end()) {
        continue;
      }
      auto muon = ambMuon.template fwdtrack_as<TMuons>();

      // loop over collisions in the TF
      for (const auto& collision : collisions) {
        const auto collId = collision.globalIndex();
        // If this is a non-orphan muons and is associated to this collisions, skip it (it is already taken into account)
        if (muon.has_collision()) {
          if (collId == muon.collisionId()) {
            continue;
          }
        }

        // TODO: maybe possibly dynamically expand the range of BCs on the ambMuon to take into account the time resolution of the
        //       currently checked collision ?
        const uint64_t mostProbableBc = collision.bc().globalBC(); // this collision's BC
        // loop over the BCs compatible with the muon
        for (const auto& bc : ambMuon.bc()) {
          if (bc.globalBC() == mostProbableBc) { // found a match
            const int64_t muonId = muon.globalIndex();
            if (collTrackIds.find(collId) == collTrackIds.end()) { // this colision is not in the map, adding it
              std::vector<int64_t> idxs{muonId};
              collTrackIds[collId] = idxs;
            } else {
              auto idxs = collTrackIds[collId];
              idxs.push_back(muonId);
              collTrackIds[collId] = idxs;
            }
            collBCmap[collId] = bc.globalBC();
          }
        }
      } // end loop over collisions
    }   // end loop over ambiguous tracks

    // associate tracks with in-bunch pileup collisions
    // Check only collisions that are associated with filtered muons
    for (auto const& [collId, assocTracks] : collTrackIds) {
      const uint64_t currentBC = collBCmap[collId];
      std::vector<int64_t> pileupCollIds{};
      for (const auto& collision : collisions) { // NOTE: since collisions should be ordered in time, one could optimize here and not make the full loop (similar as with the tracks in the time based association)
        if (collision.bc().globalBC() == currentBC) {
          if (collision.globalIndex() != collId) { // add the collision only if its a different one
            pileupCollIds.push_back(collision.globalIndex());
          }
        }
      }

      // in-bunch pileup requires at least 2 collisions in the same bunch (collId + another one)
      if (pileupCollIds.size() < 1) {
        continue;
      }

      // here we have in bunch pileup and we will associate the tracks from collId to the other (in-bunch pileup) collisions
      for (const auto& pileupCollId : pileupCollIds) {
        for (const auto& assocTrack : assocTracks) {
          if (collTrackIds.find(pileupCollId) == collTrackIds.end()) { // this colision is not in the map
            std::vector<int64_t> idxs{assocTrack};
            collTrackIds[pileupCollId] = idxs;
          } else {
            auto idxs = collTrackIds[pileupCollId];
            if (std::find(idxs.begin(), idxs.end(), assocTrack) == idxs.end()) { // add the muon only if it does not exist already
              idxs.push_back(assocTrack);
            }
            collTrackIds[pileupCollId] = idxs;
          }
        } // end loop over associated muons
      }   // end loop over pileup collisions
    }

    // create the collision - track association table
    for (auto const& [collId, assocTracks] : collTrackIds) {
      for (const auto& assocTrack : assocTracks) {
        muonAssoc(collId, assocTrack, fSelectedMuons[assocTrack]); // writes in the table (collId, fwdtrackId, filterMap)
      }
    }
  }

  template <typename TMuons>
  void associateMuonsToCollisionsAllTracks(Collisions const& collisions,
                                           TMuons const& muons,
                                           BCsWithTimestamps const& bcstimestamp,
                                           AmbiguousFwdTracks const& ambiTracksFwd)
  {
    // first processing tracks registered in the ambigous tracks table
    for (auto& ambiTrackFwd : ambiTracksFwd) {
      if (fSelectedMuons.find(ambiTrackFwd.fwdtrackId()) == fSelectedMuons.end()) {
        continue;
      }
      auto muon = ambiTrackFwd.template fwdtrack_as<TMuons>();
      if (muon.collisionId() < 0) {
        registry.fill(HIST("Association/AssociationTrackStatus"), 1);
      } else {
        registry.fill(HIST("Association/AssociationTrackStatus"), 0);
      }
      std::vector<int> vtxList;
      const auto& bcSlice = ambiTrackFwd.bc();
      int64_t trackBC = -1;
      if (bcSlice.size() != 0) {
        auto first = bcSlice.begin();
        trackBC = first.globalBC();
      }
      double t0 = muon.trackTime() + trackBC * o2::constants::lhc::LHCBunchSpacingNS + fTimeBias; // computing track time relative to first BC of the compatible BC slice
      double err = muon.trackTimeRes() * fSigmaTrack + fTimeMarginTrack;
      double tmin = t0 - err;
      double tmax = t0 + err;
      double vtxminOK = 0;
      double vtxmaxOK = 0;
      for (auto& vtxBracket : vtxOrdBrack) {
        double vtxmin = vtxBracket.first.first;
        double vtxmax = vtxBracket.first.second;
        if (tmax < vtxmin) {
          break; // all following collisions will be later and not compatible
        } else if (tmin > vtxmax) {
          continue; // following vertex with longer span might still match this track
        } else {
          vtxList.push_back(vtxBracket.second);
          vtxminOK = vtxmin;
          vtxmaxOK = vtxmax;
        }
      }
      isMuonReassigned[muon.globalIndex()] = -1;
      if (vtxList.size() > 1) {
        registry.fill(HIST("Association/AssociationTrackStatus"), 3); // track is still ambiguous
      } else if (vtxList.size() == 0) {
        registry.fill(HIST("Association/AssociationTrackStatus"), 4); // track is now orphan
      } else {
        isMuonReassigned[muon.globalIndex()] = vtxList.front();                             // track is non-ambiguously associated
        muonAssoc(vtxList.front(), muon.globalIndex(), fSelectedMuons[muon.globalIndex()]); // writes in the table (collId, fwdtrackId, filterMap)
        registry.fill(HIST("Association/AssociationTrackStatus"), 5);
        registry.fill(HIST("Association/DeltaT"), (vtxminOK + vtxmaxOK) / 2 - t0);
      }
    }
    auto trackBegin = fSelectedMuons.begin();
    // now processing all other tracks (which were not registered in the ambiguous table)
    for (auto trackFiltered = trackBegin; trackFiltered != fSelectedMuons.end(); trackFiltered++) {
      auto muon = muons.rawIteratorAt(trackFiltered->first);
      if (!(muon.has_collision())) {
        continue;
      }
      if (isMuonReassigned.find(muon.globalIndex()) != isMuonReassigned.end()) {
        continue;
      }
      registry.fill(HIST("Association/AssociationTrackStatus"), 2);
      std::vector<int> vtxList;
      auto collision = collisions.rawIteratorAt(muon.collisionId() - collisions.offset());
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      double t0 = muon.trackTime() + bc.globalBC() * o2::constants::lhc::LHCBunchSpacingNS - collision.collisionTime() + fTimeBias; // computing track time relative to associated collisino time
      double err = muon.trackTimeRes() * fSigmaTrack + fTimeMarginTrack;
      double tmin = t0 - err;
      double tmax = t0 + err;
      double vtxminOK = 0;
      double vtxmaxOK = 0;
      for (auto& vtxBracket : vtxOrdBrack) {
        double vtxmin = vtxBracket.first.first;
        double vtxmax = vtxBracket.first.second;
        if (tmax < vtxmin) {
          break; // all following collisions will be later and not compatible
        } else if (tmin > vtxmax) {
          continue; // following vertex with longer span might still match this track
        } else {
          vtxList.push_back(vtxBracket.second);
          vtxminOK = vtxmin;
          vtxmaxOK = vtxmax;
        }
      }
      isMuonReassigned[muon.globalIndex()] = -1;
      if (vtxList.size() > 1) {
        registry.fill(HIST("Association/AssociationTrackStatus"), 3); // track is still ambiguous
        for (auto& vtx : vtxList) {
          muonAssoc(vtx, muon.globalIndex(), fSelectedMuons[muon.globalIndex()]); // writes in the table (collId, fwdtrackId, filterMap)
        }
      } else if (vtxList.size() == 0) {
        registry.fill(HIST("Association/AssociationTrackStatus"), 4); // track is now orphan
      } else {
        isMuonReassigned[muon.globalIndex()] = vtxList.front();                             // track is non-ambiguously associated
        muonAssoc(vtxList.front(), muon.globalIndex(), fSelectedMuons[muon.globalIndex()]); // writes in the table (collId, fwdtrackId, filterMap)
        registry.fill(HIST("Association/AssociationTrackStatus"), 5);
        registry.fill(HIST("Association/DeltaT"), (vtxminOK + vtxmaxOK) / 2 - t0);
      }
    }
  }

  void processSelection(Collisions const& collisions,
                        BCsWithTimestamps const& bcstimestamp,
                        BCs const& bcs,
                        MyMuons const& muons,
                        soa::Filtered<MyMuons> const& filteredMuons,
                        AmbiguousFwdTracks const& ambFwdTracks)
  {
    runMuonSelection<gkMuonFillMap>(filteredMuons);
    if (fConfigCollisionMuonAssoc.value == 0) {
      associateMuonsToCollisionsStandard(collisions, muons);
    } else if (fConfigCollisionMuonAssoc.value == 1) {
      associateMuonsToCollisionsTime(collisions, muons, bcs);
    } else if (fConfigCollisionMuonAssoc.value == 2) {
      associateMuonsToCollisionsAmbigous(collisions, muons, bcs, ambFwdTracks);
    } else {
      runCollisionMap(collisions, bcstimestamp);
      associateMuonsToCollisionsAllTracks(collisions, muons, bcstimestamp, ambFwdTracks);
      isMuonReassigned.clear();
      vtxOrdBrack.clear();
    }
  }
  void processDummy(MyMuons&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQMuonsSelection, processSelection, "Run muon selection", false);
  PROCESS_SWITCH(DQMuonsSelection, processDummy, "Dummy function", false);
};

struct DQFilterPPTask {
  Produces<aod::DQEventFilter> eventFilter;
  Produces<aod::DqFilters> dqtable;
  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TH1I> fStats{"Statistics"};
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigBarrelSelections{"cfgBarrelSels", "jpsiPID1:pairMassLow:1", "<track-cut>:[<pair-cut>]:<n>,[<track-cut>:[<pair-cut>]:<n>],..."};
  Configurable<std::string> fConfigMuonSelections{"cfgMuonSels", "muonQualityCuts:pairNoCut:1", "<muon-cut>:[<pair-cut>]:<n>"};
  Configurable<bool> fConfigQA{"cfgWithQA", false, "If true, fill QA histograms"};

  int fNBarrelCuts;                                    // number of barrel selections
  int fNMuonCuts;                                      // number of muon selections
  std::vector<bool> fBarrelRunPairing;                 // bit map on whether the selections require pairing (barrel)
  std::vector<bool> fMuonRunPairing;                   // bit map on whether the selections require pairing (muon)
  std::vector<int> fBarrelNreqObjs;                    // minimal number of tracks/pairs required (barrel)
  std::vector<int> fMuonNreqObjs;                      // minimal number of tracks/pairs required (muon)
  std::map<int, AnalysisCompositeCut> fBarrelPairCuts; // map of barrel pair cuts
  std::map<int, AnalysisCompositeCut> fMuonPairCuts;   // map of muon pair cuts
  std::map<int, TString> fBarrelPairHistNames;         // map with names of the barrel pairing histogram directories
  std::map<int, TString> fMuonPairHistNames;           // map with names of the muon pairing histogram directories

  std::map<uint64_t, uint64_t> fFiltersMap;           // map of filters for events that passed at least one filter
  std::map<uint64_t, std::vector<bool>> fCEFPfilters; // map of CEFP filters for events that passed at least one filter

  void DefineCuts()
  {
    TString barrelSelsStr = fConfigBarrelSelections.value;
    std::unique_ptr<TObjArray> objArray(barrelSelsStr.Tokenize(","));
    fNBarrelCuts = objArray->GetEntries();
    if (fNBarrelCuts) {
      for (int icut = 0; icut < fNBarrelCuts; ++icut) {
        TString selStr = objArray->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() < 2 || sel->GetEntries() > 3) {
          continue;
        }
        // if "sel" contains 3 entries, it means the user provided both the track and pair cuts
        if (sel->GetEntries() == 3) {
          fBarrelPairCuts[icut] = (*dqcuts::GetCompositeCut(sel->At(1)->GetName()));
          fBarrelRunPairing.push_back(true);
          fBarrelNreqObjs.push_back(std::atoi(sel->At(2)->GetName()));
          fBarrelPairHistNames[icut] = Form("PairsBarrelSEPM_%s_%s", sel->At(0)->GetName(), sel->At(1)->GetName());
        } else {
          fBarrelNreqObjs.push_back(std::atoi(sel->At(1)->GetName()));
          fBarrelRunPairing.push_back(false);
        }
      }
    }
    TString muonSelsStr = fConfigMuonSelections.value;
    std::unique_ptr<TObjArray> objArray2(muonSelsStr.Tokenize(","));
    fNMuonCuts = objArray2->GetEntries();
    if (fNMuonCuts) {
      for (int icut = 0; icut < fNMuonCuts; ++icut) {
        TString selStr = objArray2->At(icut)->GetName();
        std::unique_ptr<TObjArray> sel(selStr.Tokenize(":"));
        if (sel->GetEntries() < 2 || sel->GetEntries() > 3) {
          continue;
        }
        // if "sel" contains 3 entries, it means the user provided both the track and pair cuts
        if (sel->GetEntries() == 3) {
          fMuonPairCuts[icut] = (*dqcuts::GetCompositeCut(sel->At(1)->GetName()));
          fMuonRunPairing.push_back(true);
          fMuonNreqObjs.push_back(std::atoi(sel->At(2)->GetName()));
          fMuonPairHistNames[icut] = Form("PairsForwardSEPM_%s_%s", sel->At(0)->GetName(), sel->At(1)->GetName());
        } else {
          fMuonNreqObjs.push_back(std::atoi(sel->At(1)->GetName()));
          fMuonRunPairing.push_back(false);
        }
      }
    }
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);

    // setup the Stats histogram
    fStats.setObject(new TH1I("Statistics", "Stats for DQ triggers", fNBarrelCuts + fNMuonCuts + 2, -2.5, -0.5 + fNBarrelCuts + fNMuonCuts));
    fStats->GetXaxis()->SetBinLabel(1, "Events inspected");
    fStats->GetXaxis()->SetBinLabel(2, "Events selected");
    if (fNBarrelCuts) {
      for (int ib = 3; ib < 3 + fNBarrelCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray->At(ib - 3)->GetName());
      }
    }
    if (fNMuonCuts) {
      for (int ib = 3 + fNBarrelCuts; ib < 3 + fNBarrelCuts + fNMuonCuts; ib++) {
        fStats->GetXaxis()->SetBinLabel(ib, objArray2->At(ib - 3 - fNBarrelCuts)->GetName());
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    DefineCuts();

    if (fConfigQA) {
      // initialize the variable manager
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      TString histNames = "";
      for (const auto& [key, value] : fBarrelPairHistNames) {
        histNames += value;
        histNames += ";";
      }
      for (const auto& [key, value] : fMuonPairHistNames) {
        histNames += value;
        histNames += ";";
      }
      DefineHistograms(fHistMan, histNames.Data());
      VarManager::SetUseVars(fHistMan->GetUsedVars());
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons>
  uint64_t runFilterPP(TEvent const& collision,
                       aod::BCs const& bcs,
                       TTracks const& tracksBarrel,
                       TMuons const& muons,
                       DQTrackAssoc const& barrelAssocs, DQMuonAssoc const& muonAssocs)
  {
    // Reset the values array and compute event quantities
    VarManager::ResetValues(0, VarManager::kNVars);
    VarManager::FillEvent<TEventFillMap>(collision); // event properties could be needed for cuts or histogramming

    std::vector<int> objCountersBarrel(fNBarrelCuts, 0); // init all counters to zero
    // count the number of barrel tracks fulfilling each cut
    for (auto trackAssoc : barrelAssocs) {
      for (int i = 0; i < fNBarrelCuts; ++i) {
        if (trackAssoc.isDQBarrelSelected() & (uint32_t(1) << i)) {
          objCountersBarrel[i] += 1;
        }
      }
    }

    // check which selections require pairing
    uint32_t pairingMask = 0; // in order to know which of the selections actually require pairing
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (fBarrelRunPairing[i]) {
        if (objCountersBarrel[i] > 1) { // pairing has to be enabled and at least two tracks are needed
          pairingMask |= (uint32_t(1) << i);
        }
        objCountersBarrel[i] = 0; // reset counters for selections where pairing is needed (count pairs instead)
      }
    }

    // run pairing if there is at least one selection that requires it
    uint32_t pairFilter = 0;
    if (pairingMask > 0) {
      // run pairing on the collision grouped associations
      for (auto& [a1, a2] : combinations(barrelAssocs, barrelAssocs)) {
        // check the pairing mask and that the tracks share a cut bit
        pairFilter = pairingMask & a1.isDQBarrelSelected() & a2.isDQBarrelSelected();
        if (pairFilter == 0) {
          continue;
        }

        // get the tracks from the index stored in the association
        auto t1 = a1.template track_as<TTracks>();
        auto t2 = a2.template track_as<TTracks>();
        // keep just opposite-sign pairs
        if (t1.sign() * t2.sign() > 0) {
          continue;
        }

        // construct the pair and apply pair cuts
        VarManager::FillPair<VarManager::kDecayToEE, TTrackFillMap>(t1, t2); // compute pair quantities
        for (int icut = 0; icut < fNBarrelCuts; icut++) {
          if (!(pairFilter & (uint32_t(1) << icut))) {
            continue;
          }
          if (!fBarrelPairCuts[icut].IsSelected(VarManager::fgValues)) {
            continue;
          }
          objCountersBarrel[icut] += 1; // count the pair
          if (fConfigQA) {              // fill histograms if QA is enabled
            fHistMan->FillHistClass(fBarrelPairHistNames[icut].Data(), VarManager::fgValues);
          }
        }
      }
    }

    std::vector<int> objCountersMuon(fNMuonCuts, 0); // init all counters to zero
    // count the number of muon tracks fulfilling each selection
    for (auto muon : muonAssocs) {
      for (int i = 0; i < fNMuonCuts; ++i) {
        if (muon.isDQMuonSelected() & (uint32_t(1) << i)) {
          objCountersMuon[i] += 1;
        }
      }
    }

    // check which muon selections require pairing
    pairingMask = 0; // reset the mask for the muons
    for (int i = 0; i < fNMuonCuts; i++) {
      if (fMuonRunPairing[i]) { // pairing has to be enabled and at least two tracks are needed
        if (objCountersMuon[i] > 1) {
          pairingMask |= (uint32_t(1) << i);
        }
        objCountersMuon[i] = 0; // reset counters for selections where pairing is needed (count pairs instead)
      }
    }

    // run pairing if there is at least one selection that requires it
    pairFilter = 0;
    if (pairingMask > 0) {
      // pairing is done using the collision grouped muon associations
      for (auto& [a1, a2] : combinations(muonAssocs, muonAssocs)) {

        // check the pairing mask and that the tracks share a cut bit
        pairFilter = pairingMask & a1.isDQMuonSelected() & a2.isDQMuonSelected();
        if (pairFilter == 0) {
          continue;
        }

        // get the real muon tracks
        auto t1 = a1.template fwdtrack_as<TMuons>();
        auto t2 = a2.template fwdtrack_as<TMuons>();
        // keep just opposite-sign pairs
        if (t1.sign() * t2.sign() > 0) {
          continue;
        }

        // construct the pair and apply cuts
        VarManager::FillPair<VarManager::kDecayToMuMu, TTrackFillMap>(t1, t2); // compute pair quantities
        for (int icut = 0; icut < fNMuonCuts; icut++) {
          if (!(pairFilter & (uint32_t(1) << icut))) {
            continue;
          }
          if (!fMuonPairCuts[icut].IsSelected(VarManager::fgValues)) {
            continue;
          }
          objCountersMuon[icut] += 1;
          if (fConfigQA) {
            fHistMan->FillHistClass(fMuonPairHistNames[icut].Data(), VarManager::fgValues);
          }
        }
      }
    }

    // compute the decisions and publish
    uint64_t filter = 0;
    for (int i = 0; i < fNBarrelCuts; i++) {
      if (objCountersBarrel[i] >= fBarrelNreqObjs[i]) {
        filter |= (uint64_t(1) << i);
      }
    }
    for (int i = 0; i < fNMuonCuts; i++) {
      if (objCountersMuon[i] >= fMuonNreqObjs[i]) {
        filter |= (uint64_t(1) << (i + fNBarrelCuts));
      }
    }
    return filter;
  }

  Preslice<DQTrackAssoc> trackIndicesPerCollision = aod::dqppfilter::collisionId;
  Preslice<DQMuonAssoc> muonIndicesPerCollision = aod::dqppfilter::collisionId;

  void processFilterPP(MyEventsSelected const& collisions,
                       aod::BCs const& bcs,
                       MyBarrelTracks const& tracks,
                       MyMuons const& muons,
                       DQTrackAssoc const& trackAssocs, DQMuonAssoc const& muonAssocs)
  {
    fFiltersMap.clear();
    fCEFPfilters.clear();

    uint64_t barrelMask = 0;
    for (int i = 0; i < fNBarrelCuts; i++) {
      barrelMask |= (uint64_t(1) << i);
    }
    uint64_t muonMask = 0;
    for (int i = fNBarrelCuts; i < fNBarrelCuts + fNMuonCuts; i++) {
      muonMask |= (uint64_t(1) << i);
    }

    // Loop over collisions
    int event = 0;
    int eventsFired = 0;
    for (const auto& collision : collisions) {
      // skip those that do not pass our selection
      if (!collision.isDQEventSelected()) {
        event++;
        continue;
      }
      // group the tracks and muons for this collision
      auto groupedTrackIndices = trackAssocs.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      auto groupedMuonIndices = muonAssocs.sliceBy(muonIndicesPerCollision, collision.globalIndex());

      uint64_t filter = 0;
      // if there is at least one track or muon, run the filtering function and compute triggers
      if (groupedTrackIndices.size() > 0 || groupedMuonIndices.size() > 0) {
        filter = runFilterPP<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracks, muons, groupedTrackIndices, groupedMuonIndices);
      }
      if (filter == 0) {
        event++;
        continue;
      }
      eventsFired++;
      // compute the CEPF decisions (this is done in a spacial setup with exactly kNTriggersDQ configured triggers)
      std::vector<bool> decisions(kNTriggersDQ, false); // event decisions to be transmitted to CEFP
      for (int i = 0; i < fNBarrelCuts; i++) {
        if (filter & (uint64_t(1) << i)) {
          if (i < kNTriggersDQ) {
            decisions[i] = true;
          }
        }
      }
      for (int i = fNBarrelCuts; i < fNBarrelCuts + fNMuonCuts; i++) {
        if (filter & (uint64_t(1) << i)) {
          if (i < kNTriggersDQ) {
            decisions[i] = true;
          }
        }
      }

      // if this collision fired at least one input, add it to the map, or if it is there already, update the decisions with a logical OR
      // This may happen in the case when some collisions beyond the iterator are added because they contain ambiguous tracks fired on by another collision
      if (fFiltersMap.find(collision.globalIndex()) == fFiltersMap.end()) {
        fFiltersMap[collision.globalIndex()] = filter;
        fCEFPfilters[collision.globalIndex()] = decisions;
      } else { // this collision was already fired, possible via collision - track association; add as an OR the new decisions
        fFiltersMap[collision.globalIndex()] |= filter;
        for (int i = 0; i < kNTriggersDQ; i++) {
          if (decisions[i]) {
            fCEFPfilters[collision.globalIndex()][i] = true;
          }
        }
      }

      // Now check through the associated tracks / fwdtracks and assign the same filter to their parent collisions
      // This is needed since if a collision was selected because of a track association from a neighbouring collision,
      //   then one needs to select also that collision in order to be able to redo the pairing at analysis time.
      if (filter & barrelMask) {
        for (auto& a : groupedTrackIndices) {
          auto t = a.template track_as<MyBarrelTracks>();
          auto tColl = t.collisionId();
          if (tColl == collision.globalIndex()) { // track from this collision, nothing to do
            continue;
          } else {
            if (fFiltersMap.find(tColl) == fFiltersMap.end()) {
              fFiltersMap[tColl] = filter;
              fCEFPfilters[tColl] = decisions;
            } else { // this collision was already fired, possible via collision - track association; add as an OR the new decisions
              fFiltersMap[tColl] |= filter;
              for (int i = 0; i < kNTriggersDQ; i++) {
                if (decisions[i]) {
                  fCEFPfilters[tColl][i] = true;
                }
              }
            }
          }
        }
      }
      // Do the same for muons
      if (filter & muonMask) {
        for (auto& a : groupedMuonIndices) {
          auto t = a.template fwdtrack_as<MyMuons>();
          auto tColl = t.collisionId();
          if (tColl == collision.globalIndex()) { // track from this collision, nothing to do
            continue;
          } else {
            if (fFiltersMap.find(tColl) == fFiltersMap.end()) {
              fFiltersMap[tColl] = filter;
              fCEFPfilters[tColl] = decisions;
            } else { // this collision was already fired, possible via collision - track association; add as an OR the new decisions
              fFiltersMap[tColl] |= filter;
              for (int i = 0; i < kNTriggersDQ; i++) {
                if (decisions[i]) {
                  fCEFPfilters[tColl][i] = true;
                }
              }
            }
          }
        }
      }
      event++;
    }

    // At this point, we have all the non-null decisions for all collisions.
    // we loop again over collisions and create the decision tables
    // NOTE: For the CEFP decisions, decisions are placed in a vector of bool in an ordered way:
    //       start with all configured barrel selections and then continue with those from muons
    //       The configured order has to be in sync with that implemented in the cefp task and can be done
    //       by preparing a dedicated json configuration file
    int totalEventsTriggered = 0;
    for (const auto& collision : collisions) {
      fStats->Fill(-2.0);
      if (!collision.isDQEventSelected()) {
        eventFilter(0);
        dqtable(false, false, false, false, false);
      }
      fStats->Fill(-1.0);

      if (fFiltersMap.find(collision.globalIndex()) == fFiltersMap.end()) {
        eventFilter(0);
        dqtable(false, false, false, false, false);
      } else {
        totalEventsTriggered++;
        for (int i = 0; i < fNBarrelCuts + fNMuonCuts; i++) {
          if (fFiltersMap[collision.globalIndex()] & (uint32_t(1) << i))
            fStats->Fill(static_cast<float>(i));
        }
        eventFilter(fFiltersMap[collision.globalIndex()]);
        auto dqDecisions = fCEFPfilters[collision.globalIndex()];
        dqtable(dqDecisions[0], dqDecisions[1], dqDecisions[2], dqDecisions[3], dqDecisions[4]);
      }
    }

    cout << "-------------------- In this TF, eventsFired / totalTriggered :: " << eventsFired << "/" << totalEventsTriggered << endl;
  }

  // TODO: dummy function for the case when no process function is enabled
  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPPTask, processFilterPP, "Run filter task", false);
  PROCESS_SWITCH(DQFilterPPTask, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQEventSelectionTask>(cfgc),
    adaptAnalysisTask<DQBarrelTrackSelection>(cfgc),
    adaptAnalysisTask<DQMuonsSelection>(cfgc),
    adaptAnalysisTask<DQFilterPPTask>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,vtxPbPb");
    }

    if (classStr.Contains("Track")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "its,tpcpid,dca");
    }

    if (classStr.Contains("Muon")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "muon");
    }

    if (classStr.Contains("Pairs")) {
      if (classStr.Contains("Barrel")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "vertexing-barrel");
      }
      if (classStr.Contains("Forward")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", "dimuon,vertexing-forward");
      }
    }
  }
}
