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
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithCov = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended, aod::TrackSelection,
                                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                        aod::pidTPCFullKa, aod::pidTPCFullPr,
                                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                        aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyBarrelTracksWithV0Bits = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended, aod::TrackSelection,
                                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                           aod::pidTPCFullKa, aod::pidTPCFullPr,
                                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                           aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::V0Bits>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEventsWithFilter = soa::Join<aod::Collisions, aod::EvSels, aod::DQEventFilter>;
using MyEventsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
using MyMuons = aod::FwdTracks;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkEventFillMapWithCent = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision | VarManager::ObjTypes::CollisionCent;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkTrackFillMapWithV0Bits = gkTrackFillMap | VarManager::ObjTypes::TrackV0Bits;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::Muon;
constexpr static uint32_t gkMuonFillMapWithCov = VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov;

struct TableMaker {

  Produces<ReducedEvents> event;
  Produces<ReducedEventsExtended> eventExtended;
  Produces<ReducedEventsVtxCov> eventVtxCov;
  Produces<ReducedTracks> trackBasic;
  Produces<ReducedTracksBarrel> trackBarrel;
  Produces<ReducedTracksBarrelCov> trackBarrelCov;
  Produces<ReducedTracksBarrelPID> trackBarrelPID;
  Produces<ReducedMuons> muonBasic;
  Produces<ReducedMuonsExtra> muonExtra;
  Produces<ReducedMuonsCov> muonCov;

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics
  HistogramManager* fHistMan;

  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "jpsiPID1", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 1.0f, "Low pt cut for tracks in the barrel"};
  Configurable<float> fConfigMuonPtLow{"cfgMuonLowPt", 1.0f, "Low pt cut for muons"};
  Configurable<float> fConfigMinTpcSignal{"cfgMinTpcSignal", 30.0, "Minimum TPC signal"};
  Configurable<float> fConfigMaxTpcSignal{"cfgMaxTpcSignal", 300.0, "Maximum TPC signal"};
  Configurable<bool> fConfigNoQA{"cfgNoQA", false, "If true, no QA histograms"};
  Configurable<bool> fConfigDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes and more)"};
  Configurable<bool> fIsRun2{"cfgIsRun2", false, "Whether we analyze Run-2 or Run-3 data"};

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fMuonCuts;  //! Muon track cuts

  // TODO: filter on TPC dedx used temporarily until electron PID will be improved
  Filter barrelSelectedTracks = ifnode(fIsRun2.node() == true, aod::track::trackType == uint8_t(aod::track::Run2Track), aod::track::trackType == uint8_t(aod::track::Track)) && o2::aod::track::pt >= fConfigBarrelTrackPtLow && nabs(o2::aod::track::eta) <= 0.9f && o2::aod::track::tpcSignal >= fConfigMinTpcSignal && o2::aod::track::tpcSignal <= fConfigMaxTpcSignal && o2::aod::track::tpcChi2NCl < 4.0f && o2::aod::track::itsChi2NCl < 36.0f;

  Filter muonFilter = o2::aod::fwdtrack::pt >= fConfigMuonPtLow;

  void init(o2::framework::InitContext& context)
  {
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";
    if (fConfigDetailedQA) {
      histClasses += "Event_BeforeCuts;";
    }
    histClasses += "Event_AfterCuts;";

    bool enableBarrelHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                               context.mOptions.get<bool>("processFullWithCent") ||
                               context.mOptions.get<bool>("processBarrelOnly") || context.mOptions.get<bool>("processBarrelOnlyWithCent") ||
                               context.mOptions.get<bool>("processBarrelOnlyWithCov") || context.mOptions.get<bool>("processBarrelOnlyWithEventFilter"));
    bool enableMuonHistos = (context.mOptions.get<bool>("processFull") || context.mOptions.get<bool>("processFullWithCov") ||
                             context.mOptions.get<bool>("processFullWithCent") ||
                             context.mOptions.get<bool>("processMuonOnly") || context.mOptions.get<bool>("processMuonOnlyWithCent") ||
                             context.mOptions.get<bool>("processMuonOnlyWithCov") || context.mOptions.get<bool>("processMuonOnlyWithFilter"));

    if (enableBarrelHistos) {
      if (fConfigDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
      }
      for (auto& cut : fTrackCuts) {
        histClasses += Form("TrackBarrel_%s;", cut.GetName());
      }
    }
    if (enableMuonHistos) {
      if (fConfigDetailedQA) {
        histClasses += "Muons_BeforeCuts;";
      }
      for (auto& muonCut : fMuonCuts) {
        histClasses += Form("Muons_%s;", muonCut.GetName());
      }
    }

    DefineHistograms(histClasses);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
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

  // Templated function instantianed for all of the process functions
  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, uint32_t TMuonFillMap, typename TEvent, typename TTracks, typename TMuons>
  void fullSkimming(TEvent const& collision, aod::BCs const& bcs, TTracks const& tracksBarrel, TMuons const& tracksMuon)
  {
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
    // TODO: Add the event level decisions from the filtering task into the tag

    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventFillMap>(collision); // extract event information and place it in the fValues array
    if (fConfigDetailedQA) {
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
      return;
    }

    // fill stats information, after selections
    for (int i = 0; i < kNaliases; i++) {
      if (triggerAliases & (uint32_t(1) << i)) {
        ((TH2I*)fStatsList->At(0))->Fill(3.0, float(i));
      }
    }
    ((TH2I*)fStatsList->At(0))->Fill(3.0, float(kNaliases));

    if (!fConfigNoQA) {
      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);
    }

    // create the event tables
    event(tag, collision.bc().runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
    eventExtended(collision.bc().globalBC(), collision.bc().triggerMask(), 0, triggerAliases, VarManager::fgValues[VarManager::kCentVZERO]);
    eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

    uint64_t trackFilteringTag = 0;
    uint8_t trackTempFilterMap = 0;
    if constexpr (static_cast<bool>(TTrackFillMap)) {
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
        trackBarrelCov.reserve(tracksBarrel.size());
      }
      trackBarrelPID.reserve(tracksBarrel.size());

      // loop over tracks
      for (auto& track : tracksBarrel) {
        trackFilteringTag = uint64_t(0);
        trackTempFilterMap = uint8_t(0);
        VarManager::FillTrack<TTrackFillMap>(track);
        if (fConfigDetailedQA) {
          fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
        }
        // apply track cuts and fill stats histogram
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            trackTempFilterMap |= (uint8_t(1) << i);
            if (!fConfigNoQA) {
              fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut).GetName()), VarManager::fgValues);
            }
            ((TH1I*)fStatsList->At(1))->Fill(float(i));
          }
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

        // create the track tables
        trackBasic(event.lastIndex(), trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign());
        trackBarrel(track.tpcInnerParam(), track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                    track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.tpcChi2NCl(),
                    track.trdChi2(), track.trdPattern(), track.tofChi2(),
                    track.length(), track.dcaXY(), track.dcaZ());
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov)) {
          trackBarrelCov(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                         track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                         track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                         track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
        }
        trackBarrelPID(track.tpcSignal(),
                       track.tpcNSigmaEl(), track.tpcNSigmaMu(),
                       track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                       track.beta(),
                       track.tofNSigmaEl(), track.tofNSigmaMu(),
                       track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                       track.trdSignal());
      }
    } // end if constexpr (TTrackFillMap)

    if constexpr (static_cast<bool>(TMuonFillMap)) {
      // build the muon tables
      muonBasic.reserve(tracksMuon.size());
      muonExtra.reserve(tracksMuon.size());
      if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
        muonCov.reserve(tracksMuon.size());
      }
      // loop over muons
      for (auto& muon : tracksMuon) {
        trackFilteringTag = uint64_t(0);
        trackTempFilterMap = uint8_t(0);

        VarManager::FillTrack<TMuonFillMap>(muon);
        if (fConfigDetailedQA) {
          fHistMan->FillHistClass("Muons_BeforeCuts", VarManager::fgValues);
        }
        // apply the muon selection cuts and fill the stats histogram
        int i = 0;
        for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            trackTempFilterMap |= (uint8_t(1) << i);
            if (!fConfigNoQA) {
              fHistMan->FillHistClass(Form("Muons_%s", (*cut).GetName()), VarManager::fgValues);
            }
            ((TH1I*)fStatsList->At(2))->Fill(float(i));
          }
        }
        if (!trackTempFilterMap) {
          continue;
        }
        // store the cut decisions
        trackFilteringTag |= uint64_t(trackTempFilterMap); // BIT0-7:  user selection cuts

        muonBasic(event.lastIndex(), trackFilteringTag, muon.pt(), muon.eta(), muon.phi(), muon.sign());
        muonExtra(muon.nClusters(), muon.pDca(), muon.rAtAbsorberEnd(),
                  muon.chi2(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(),
                  muon.matchScoreMCHMFT(), muon.matchMFTTrackId(), muon.matchMCHTrackId(), muon.mchBitMap(), muon.midBitMap(), muon.midBoards(), muon.trackType());
        if constexpr (static_cast<bool>(TMuonFillMap & VarManager::ObjTypes::MuonCov)) {
          muonCov(muon.x(), muon.y(), muon.z(), muon.phi(), muon.tgl(), muon.signed1Pt(),
                  muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(), muon.cPhiPhi(),
                  muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(), muon.c1PtX(), muon.c1PtY(),
                  muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2());
        }
      }
    }
  } // end fullSkimming()

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      fHistMan->AddHistClass(classStr.Data());

      if (classStr.Contains("Event")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "triggerall,cent");
      }

      if (classStr.Contains("Track")) {
        if (fConfigDetailedQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "dca,its,tpcpid,tofpid");
        } else {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track");
        }
      }
      if (classStr.Contains("Muons")) {
        if (fConfigDetailedQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "muon");
        } else {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track");
        }
      }
    }

    // create statistics histograms (event, tracks, muons)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, kNaliases + 1, -0.5, kNaliases + 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    for (int ib = 1; ib <= kNaliases; ib++) {
      histEvents->GetYaxis()->SetBinLabel(ib, aliasLabels[ib - 1]);
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

  // Produce barrel + muon tables -------------------------------------------------------------------------------------------------------------
  void processFull(MyEvents::iterator const& collision, aod::BCs const& bcs,
                   soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracksBarrel, tracksMuon);
  }

  void processFullWithCov(MyEvents::iterator const& collision, aod::BCs const& bcs,
                          soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel, soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, gkMuonFillMapWithCov>(collision, bcs, tracksBarrel, tracksMuon);
  }

  // Produce barrel + muon tables, with centrality --------------------------------------------------------------------------------------------
  void processFullWithCent(MyEventsWithCent::iterator const& collision, aod::BCs const& bcs,
                           soa::Filtered<MyBarrelTracks> const& tracksBarrel, soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMap, gkMuonFillMap>(collision, bcs, tracksBarrel, tracksMuon);
  }

  // Produce barrel only tables, with V0Bits ------------------------------------------------------------------------------------------------
  void processBarrelOnlyWithV0Bits(MyEvents::iterator const& collision, aod::BCs const& bcs,
                                   soa::Filtered<MyBarrelTracksWithV0Bits> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithV0Bits, 0u>(collision, bcs, tracksBarrel, nullptr);
  }

  // Produce barrel only tables, with event filtering ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithEventFilter(MyEventsWithFilter::iterator const& collision, aod::BCs const& bcs,
                                        soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias()[i] > 0) {
        ((TH2I*)fStatsList->At(0))->Fill(1.0, float(i));
      }
    }
    ((TH2I*)fStatsList->At(0))->Fill(1.0, float(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr);
    }
  }

  // Produce barrel only tables, with centrality -----------------------------------------------------------------------------------------------
  void processBarrelOnlyWithCent(MyEventsWithCent::iterator const& collision, aod::BCs const& bcs,
                                 soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMapWithCent, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr);
  }

  // Produce barrel tables only, with track cov matrix ----------------------------------------------------------------------------------------
  void processBarrelOnlyWithCov(MyEvents::iterator const& collision, aod::BCs const& bcs,
                                soa::Filtered<MyBarrelTracksWithCov> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMapWithCov, 0u>(collision, bcs, tracksBarrel, nullptr);
  }

  // Produce barrel tables only ----------------------------------------------------------------------------------------------------------------
  void processBarrelOnly(MyEvents::iterator const& collision, aod::BCs const& bcs,
                         soa::Filtered<MyBarrelTracks> const& tracksBarrel)
  {
    fullSkimming<gkEventFillMap, gkTrackFillMap, 0u>(collision, bcs, tracksBarrel, nullptr);
  }

  // Produce muon tables only, with centrality -------------------------------------------------------------------------------------------------
  void processMuonOnlyWithCent(MyEventsWithCent::iterator const& collision, aod::BCs const& bcs,
                               soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMapWithCent, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon);
  }

  // Produce muon tables only, with muon cov matrix --------------------------------------------------------------------------------------------
  void processMuonOnlyWithCov(MyEvents::iterator const& collision, aod::BCs const& bcs,
                              soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMapWithCov>(collision, bcs, nullptr, tracksMuon);
  }

  // Produce muon tables only ------------------------------------------------------------------------------------------------------------------
  void processMuonOnly(MyEvents::iterator const& collision, aod::BCs const& bcs,
                       soa::Filtered<MyMuons> const& tracksMuon)
  {
    fullSkimming<gkEventFillMap, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon);
  }

  // Produce muon tables only, with event filtering --------------------------------------------------------------------------------------------
  void processMuonOnlyWithFilter(MyEventsWithFilter::iterator const& collision, aod::BCs const& bcs,
                                 soa::Filtered<MyMuons> const& tracksMuon)
  {
    for (int i = 0; i < kNaliases; i++) {
      if (collision.alias()[i] > 0) {
        ((TH2I*)fStatsList->At(0))->Fill(1.0, float(i));
      }
    }
    ((TH2I*)fStatsList->At(0))->Fill(1.0, float(kNaliases));
    if (collision.eventFilter()) {
      fullSkimming<gkEventFillMap, 0u, gkMuonFillMap>(collision, bcs, nullptr, tracksMuon);
    }
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

  PROCESS_SWITCH(TableMaker, processFull, "Build full DQ skimmed data model, w/o centrality", false);
  PROCESS_SWITCH(TableMaker, processFullWithCov, "Build full DQ skimmed data model, w/ track and fwdtrack covariance tables", false);
  PROCESS_SWITCH(TableMaker, processFullWithCent, "Build full DQ skimmed data model, w/ centrality", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithV0Bits, "Build full DQ skimmed data model, w/o centrality, w/ V0Bits", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithEventFilter, "Build full DQ skimmed data model, w/o centrality, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCent, "Build barrel-only DQ skimmed data model, w/ centrality", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnlyWithCov, "Build barrel-only DQ skimmed data model, w/ track cov matrix", false);
  PROCESS_SWITCH(TableMaker, processBarrelOnly, "Build barrel-only DQ skimmed data model, w/o centrality", false);
  PROCESS_SWITCH(TableMaker, processMuonOnlyWithCent, "Build muon-only DQ skimmed data model, w/ centrality", false);
  PROCESS_SWITCH(TableMaker, processMuonOnlyWithCov, "Build muon-only DQ skimmed data model, w/ muon cov matrix", false);
  PROCESS_SWITCH(TableMaker, processMuonOnly, "Build muon-only DQ skimmed data model", false);
  PROCESS_SWITCH(TableMaker, processMuonOnlyWithFilter, "Build muon-only DQ skimmed data model, w/ event filter", false);
  PROCESS_SWITCH(TableMaker, processOnlyBCs, "Analyze the BCs to store sampled lumi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TableMaker>(cfgc)};
}
