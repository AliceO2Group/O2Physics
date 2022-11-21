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
#include "TList.h"
#include <iostream>

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

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using MyEventsWithCents = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::McCollisionLabels>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkEventFillMapWithCent = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkEventMCFillMap = VarManager::ObjTypes::CollisionMC;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithDalitzBits = gkTrackFillMap | VarManager::ObjTypes::DalitzBits;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;
constexpr static uint32_t gkParticleMCFillMap = VarManager::ObjTypes::ParticleMC;

struct TableMakerMC {

  Produces<ReducedEvents> event;
  Produces<ReducedEventsExtended> eventExtended;
  Produces<ReducedEventsVtxCov> eventVtxCov;
  Produces<ReducedMCEventLabels> eventMClabels;
  Produces<ReducedMCEvents> eventMC;
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

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fMuonCuts;  //! Muon track cuts

  bool fDoDetailedQA = false; // Bool to set detailed QA true, if QA is set true

  // TODO: filter on TPC dedx used temporarily until electron PID will be improved
  Filter barrelSelectedTracks = ifnode(fIsRun2.node() == true, aod::track::trackType == uint8_t(aod::track::Run2Track), aod::track::trackType == uint8_t(aod::track::Track)) && o2::aod::track::pt >= fConfigBarrelTrackPtLow && nabs(o2::aod::track::eta) <= 0.9f;

  Filter muonFilter = o2::aod::fwdtrack::pt >= fConfigMuonPtLow;

  void init(o2::framework::InitContext& context)
  {
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
                               context.mOptions.get<bool>("processBarrelOnlyWithCent") || context.mOptions.get<bool>("processBarrelOnlyWithCov"));
    bool enableMuonHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                             context.mOptions.get<bool>("processMuonOnlyWithCent") || context.mOptions.get<bool>("processMuonOnlyWithCov"));
    // TODO: switch on/off histogram classes depending on which process function we run
    if (enableBarrelHistos) {
      if (fDoDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
      }
      if (fConfigQA) {
        for (auto& cut : fTrackCuts) {
          histClasses += Form("TrackBarrel_%s;", cut.GetName());
        }
      }
    }

    if (enableMuonHistos) {
      if (fDoDetailedQA) {
        histClasses += "Muons_BeforeCuts;";
      }
      if (fConfigQA) {
        for (auto& cut : fMuonCuts) {
          histClasses += Form("Muons_%s;", cut.GetName());
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
            histClasses += Form("TrackBarrel_BeforeCuts_%s;", objArray->At(isig)->GetName());
            for (auto& cut : fTrackCuts) {
              histClasses += Form("TrackBarrel_%s_%s;", cut.GetName(), objArray->At(isig)->GetName());
            }
          }
          if (enableMuonHistos) {
            histClasses += Form("Muons_BeforeCuts_%s;", objArray->At(isig)->GetName());
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

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons>
  void fullSkimming(TEvent const& collisions, aod::BCs const& bcs, TTracks const& tracksBarrel, TMuons const& tracksMuon,
                    aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
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
      //TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      // get the trigger aliases
      uint32_t triggerAliases = 0;
      for (int i = 0; i < kNaliases; i++) {
        if (collision.alias()[i] > 0) {
          triggerAliases |= (uint32_t(1) << i);
        }
      }
      uint64_t tag = 0;
      // store the selection decisions
      for (int i = 0; i < kNsel; i++) {
        if (collision.selection()[i] > 0) {
          tag |= (uint64_t(1) << i);
        }
      }
      if (collision.sel7()) {
        tag |= (uint64_t(1) << kNsel); //! SEL7 stored at position kNsel in the tag bit map
      }

      auto mcCollision = collision.mcCollision();
      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
      VarManager::FillEvent<gkEventMCFillMap>(mcCollision);

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }
      // fill stats information, before selections
      for (int i = 0; i < kNaliases; i++) {
        if (triggerAliases & (uint32_t(1) << i)) {
          ((TH2I*)fStatsList->At(0))->Fill(2.0, float(i));
        }
      }
      ((TH2I*)fStatsList->At(0))->Fill(2.0, float(kNaliases));

      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }

      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
      }

      // fill stats information, after selections
      for (int i = 0; i < kNaliases; i++) {
        if (triggerAliases & (uint32_t(1) << i)) {
          ((TH2I*)fStatsList->At(0))->Fill(3.0, float(i));
        }
      }
      ((TH2I*)fStatsList->At(0))->Fill(3.0, float(kNaliases));

      event(tag, collision.bc().runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
      eventExtended(collision.bc().globalBC(), collision.bc().triggerMask(), 0, triggerAliases, VarManager::fgValues[VarManager::kCentVZERO]);
      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
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
          if (sig.CheckSignal(true, mcTracks, mctrack)) {
            mcflags |= (uint16_t(1) << i);
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
            VarManager::FillTrack<gkParticleMCFillMap>(mctrack);
            int j = 0;
            for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, j++) {
              if (mcflags & (uint16_t(1) << j)) {
                fHistMan->FillHistClass(Form("MCTruth_%s", (*signal).GetName()), VarManager::fgValues);
              }
            }
          }
        }
      } // end loop over mc stack

      // loop over reconstructed tracks
      if constexpr (static_cast<bool>(TTrackFillMap)) {
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
          trackFilteringTag = uint64_t(0);
          trackTempFilterMap = uint8_t(0);
          VarManager::FillTrack<TTrackFillMap>(track);
          // If no MC particle is found, skip the track
          if (!track.has_mcParticle()) {
            continue;
          }
          auto mctrack = track.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrack<gkParticleMCFillMap>(mctrack);

          if (fDoDetailedQA) {
            fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
          }
          // apply track cuts and fill stats histogram
          int i = 0;
          for (auto& cut : fTrackCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              if (fConfigQA) {
                fHistMan->FillHistClass(Form("TrackBarrel_%s", cut.GetName()), VarManager::fgValues); // fill the reconstructed truth
              }
              ((TH1I*)fStatsList->At(1))->Fill(float(i));
            }
            i++;
          }
          if (!trackTempFilterMap) {
            continue;
          }

          // store filtering information
          if (track.isGlobalTrack()) {
            trackFilteringTag |= (uint64_t(1) << 0); // BIT0: global track
          }
          if (track.isGlobalTrackSDD()) {
            trackFilteringTag |= (uint64_t(1) << 1); // BIT1: global track SSD
          }
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackV0Bits)) { // BIT2-6: V0Bits
            trackFilteringTag |= (uint64_t(track.pidbit()) << 2);
            for (int iv0 = 0; iv0 < 5; iv0++) {
              if (track.pidbit() & (uint8_t(1) << iv0)) {
                ((TH1I*)fStatsList->At(1))->Fill(fTrackCuts.size() + float(iv0));
              }
            }
          }
          trackFilteringTag |= (uint64_t(trackTempFilterMap) << 7); // BIT7-14:  user track filters
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::DalitzBits)) {
            trackFilteringTag |= (uint64_t(track.dalitzBits()) << 15); // BIT15-...: Dalitz
          }

          mcflags = 0;
          i = 0;     // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mcTracks, mctrack)) {
              mcflags |= (uint16_t(1) << i);
              if (fDoDetailedQA) {
                j = 0;
                fHistMan->FillHistClass(Form("TrackBarrel_BeforeCuts_%s", sig.GetName()), VarManager::fgValues); // fill the reconstructed truth BeforeCuts
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

          trackBasic(event.lastIndex(), trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), 0);
          trackBarrel(track.tpcInnerParam(), track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                      track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                      track.tpcNClsShared(), track.tpcChi2NCl(),
                      track.trdChi2(), track.trdPattern(), track.tofChi2(),
                      track.length(), track.dcaXY(), track.dcaZ());
          trackBarrelPID(track.tpcSignal(),
                         track.tpcNSigmaEl(), track.tpcNSigmaMu(),
                         track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(),
                         track.tofNSigmaEl(), track.tofNSigmaMu(),
                         track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                         track.trdSignal());
          trackBarrelLabels(fNewLabels.find(mctrack.index())->second, track.mcMask(), mcflags);
          if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
            trackBarrelCov(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                           track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
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

        auto groupedMuons = tracksMuon.sliceBy(perCollisionMuons, collision.globalIndex());
        // loop over muons

        // first we need to get the correct indices
        int nDel = 0;
        int idxPrev = -1;
        std::map<int, int> newEntryNb;
        std::map<int, int> newMatchIndex;

        for (auto& muon : groupedMuons) {
          trackFilteringTag = uint64_t(0);
          trackTempFilterMap = uint8_t(0);

          if (!muon.has_mcParticle()) {
            continue;
          }

          VarManager::FillTrack<TMuonFillMap>(muon);

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
              ((TH1I*)fStatsList->At(2))->Fill(float(i));
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

          trackFilteringTag = uint64_t(0);
          trackTempFilterMap = uint8_t(0);

          if (!muon.has_mcParticle()) {
            continue;
          }
          auto mctrack = muon.template mcParticle_as<aod::McParticles_001>();
          VarManager::FillTrack<TMuonFillMap>(muon);
          VarManager::FillTrack<gkParticleMCFillMap>(mctrack);

          if (fDoDetailedQA) {
            fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
          }
          // apply the muon selection cuts and fill the stats histogram
          int i = 0;
          for (auto& cut : fMuonCuts) {
            if (cut.IsSelected(VarManager::fgValues)) {
              trackTempFilterMap |= (uint8_t(1) << i);
              fHistMan->FillHistClass(Form("Muons_%s", cut.GetName()), VarManager::fgValues);
              ((TH1I*)fStatsList->At(2))->Fill(float(i));
            }
            i++;
          }
          if (!trackTempFilterMap) {
            continue;
          }
          // store the cut decisions
          trackFilteringTag |= uint64_t(trackTempFilterMap); // BIT0-7:  user selection cuts

          mcflags = 0;
          i = 0;     // runs over the MC signals
          int j = 0; // runs over the track cuts
          // check all the specified signals and fill histograms for MC truth matched tracks
          for (auto& sig : fMCSignals) {
            if (sig.CheckSignal(true, mcTracks, mctrack)) {
              mcflags |= (uint16_t(1) << i);
              if (fDoDetailedQA) {
                fHistMan->FillHistClass(Form("Muons_BeforeCuts_%s", sig.GetName()), VarManager::fgValues); // fill the reconstructed truth BeforeCuts
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

          if (int(muon.trackType()) == 0 || int(muon.trackType()) == 2) { // MCH-MFT or GLB track
            int matchIdx = muon.matchMCHTrackId() - muon.offsets();
            if (newEntryNb.count(matchIdx) > 0) {                                                  // if the key exists which means the match will not get deleted
              newMatchIndex[muon.index()] = newEntryNb[matchIdx];                                  // update the match for this muon to the updated entry of the match
              newMatchIndex[muon.index()] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()]; // adding the offset of muons, muonBasic.lastIndex() start at -1
              if (int(muon.trackType()) == 0) {                                                    // for now only do this to global tracks
                newMatchIndex[matchIdx] = newEntryNb[muon.index()];                                // add the  updated index of this muon as a match to mch track
                newMatchIndex[matchIdx] += muonBasic.lastIndex() + 1 - newEntryNb[muon.index()];   // adding the offset, muonBasic.lastIndex() start at -1
              }
            } else {
              newMatchIndex[muon.index()] = -1;
            }
          }

          else if (int(muon.trackType() == 4)) { // an MCH track
            // in this case the matches should be filled from the other types but we need to check
            if (newMatchIndex.count(muon.index()) == 0) {
              newMatchIndex[muon.index()] = -1;
            }
          }

          muonBasic(event.lastIndex(), trackFilteringTag, muon.pt(), muon.eta(), muon.phi(), muon.sign(), 0);
          muonExtra(muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
                    muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                    muon.matchScoreMCHMFT(), newMatchIndex.find(muon.index())->second, muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), muon.trackType(), muon.fwdDcaX(), muon.fwdDcaY());
          if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
            muonCov(muon.x(), muon.y(), muon.z(), muon.phi(), muon.tgl(), muon.signed1Pt(),
                    muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(), muon.cPhiPhi(),
                    muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(), muon.c1PtX(), muon.c1PtY(),
                    muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2());
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
        if (mcflags & (uint16_t(1) << isig)) {
          ((TH1I*)fStatsList->At(3))->Fill(float(isig));
        }
      }
      if (mcflags == 0) {
        ((TH1I*)fStatsList->At(3))->Fill(float(fMCSignals.size()));
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
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "mctruth", histMCTruthName);
        }
      }
    }

    // create statistics histograms (event, tracks, muons, MCsignals)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, kNaliases + 1, -0.5, kNaliases + 0.5);
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
  void processFull(MyEvents const& collisions, aod::BCs const& bcs,
                   soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon,
                   aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collisions, bcs, tracksBarrel, tracksMuon, mcEvents, mcTracks);
  }

  void processFullWithCov(MyEvents const& collisions, aod::BCs const& bcs,
                          soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon,
                          aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(collisions, bcs, tracksBarrel, tracksMuon, mcEvents, mcTracks);
  }

  // Produce barrel only tables ------------------------------------------------------------------------------------
  void processBarrelOnly(MyEvents const& collisions, aod::BCs const& bcs,
                         soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                         aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks);
  }
  // Produce barrel only tables, with centrality ------------------------------------------------------------------------------------
  void processBarrelOnlyWithCent(MyEventsWithCents const& collisions, aod::BCs const& bcs,
                                 soa::Filtered<MyBarrelTracks> const& tracksBarrel,
                                 aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMap, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks);
  }
  // Produce barrel only tables, with cov matrix-----------------------------------------------------------------------
  void processBarrelOnlyWithCov(MyEvents const& collisions, aod::BCs const& bcs,
                                soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel,
                                aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks);
  }
  // Produce barrel only tables, with cov matrix and dalitz bits-----------------------------------------------------------------------
  void processBarrelOnlyWithDalitzBits(MyEvents const& collisions, aod::BCs const& bcs,
                                       soa::Filtered<MyBarrelTracksWithDalitzBits> const& tracksBarrel,
                                       aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithDalitzBits, 0u>(collisions, bcs, tracksBarrel, nullptr, mcEvents, mcTracks);
  }

  // Produce muon only tables ------------------------------------------------------------------------------------
  /*void processMuonOnly(MyEvents const& collisions, aod::BCs const& bcs,
                       soa::Filtered<MyMuons> const& tracksMuon,
                       aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMap>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks);
  }*/
  // Produce muon only tables, with centrality-------------------------------------------------------------------------------
  void processMuonOnlyWithCent(MyEventsWithCents const& collisions, aod::BCs const& bcs,
                               soa::Filtered<MyMuons> const& tracksMuon,
                               aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMapWithCent, 0u, gkMuonFillMap>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks);
  }
  // Produce muon only tables, with cov matrix ------------------------------------------------------------------------------------
  void processMuonOnlyWithCov(MyEvents const& collisions, aod::BCs const& bcs,
                              soa::Filtered<MyMuonsWithCov> const& tracksMuon,
                              aod::McCollisions const& mcEvents, aod::McParticles_001 const& mcTracks)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov>(collisions, bcs, nullptr, tracksMuon, mcEvents, mcTracks);
  }
  // Process the BCs and store stats for luminosity retrieval -----------------------------------------------------------------------------------
  void processOnlyBCs(soa::Join<aod::BCs, aod::BcSels>::iterator const& bc)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (bc.alias()[i] > 0) {
        ((TH2I*)fStatsList->At(0))->Fill(0.0, float(i));
      }
    }
    ((TH2I*)fStatsList->At(0))->Fill(0.0, float(kNaliases));
  }

  PROCESS_SWITCH(TableMakerMC, processFull, "Produce both barrel and muon skims", false);
  PROCESS_SWITCH(TableMakerMC, processFullWithCov, "Produce both barrel and muon skims, w/ track and fwdtrack cov tables", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnly, "Produce barrel skims", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithCent, "Produce barrel skims, w/ centrality", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithCov, "Produce barrel skims, with track covariance matrix", false);
  PROCESS_SWITCH(TableMakerMC, processBarrelOnlyWithDalitzBits, "Produce barrel skims, and dalitz bits", false);
  // PROCESS_SWITCH(TableMakerMC, processMuonOnly, "Produce muon skims", false);
  PROCESS_SWITCH(TableMakerMC, processMuonOnlyWithCov, "Produce muon skims, with muon covariance matrix", false);
  PROCESS_SWITCH(TableMakerMC, processMuonOnlyWithCent, "Produce muon skims, w/ centrality", false);
  PROCESS_SWITCH(TableMakerMC, processOnlyBCs, "Analyze the BCs to store sampled lumi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // TODO: For now TableMakerMC works just for PbPb (cent table is present)
  //      Implement workflow arguments for pp/PbPb and possibly merge the task with tableMaker.cxx
  return WorkflowSpec{
    adaptAnalysisTask<TableMakerMC>(cfgc)};
}
