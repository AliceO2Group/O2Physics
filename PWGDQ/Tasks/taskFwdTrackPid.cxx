// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskFwdTrackPid.cxx
/// \brief Task for the analysis of forward PID with MFT
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/CCDB/EventSelectionParams.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ITSMFTBase/DPLAlpideParam.h"

#include "TGeoGlobalMagField.h"
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsMC = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedMCEventLabels>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksMC = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;
using MyMftTracks = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra>;
using MyMftTracksMC = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra, aod::ReducedMFTLabels>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMCEventFillMap = VarManager::ObjTypes::ReducedEventMC;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

void DefineHistograms(HistogramManager* histMan, TString histClasses);

struct taskFwdTrackPid {
  Produces<aod::FwdPidsAll> fwdPidAllList;

  HistogramManager* fHistMan;
  OutputObj<THashList> fOutputList{"output"};

  Configurable<float> fConfigMaxDCA{"cfgMaxDCA", 0.5f, "Manually set maximum DCA of the track"};
  Configurable<float> downSampleFactor{"downSampleFactor", 1., "Fraction of candidates to keep for ML"};
  Configurable<std::string> fConfigMCGenSignals{"cfgMCGenSignals", "", "Comma separated list of MC signals (generated)"};
  Configurable<std::string> fConfigMCRecSignals{"cfgMCRecSignals", "", "Comma separated list of MC signals (reconstructed)"};

  std::vector<TString> fGenMCSignalsNames;
  std::vector<TString> fRecMCSignalsNames;
  std::vector<MCSignal> fGenMCSignals;
  std::vector<MCSignal> fRecMCSignals;

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
    }

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    TString histNames = "";

    TString sigGenNamesStr = fConfigMCGenSignals.value;
    std::unique_ptr<TObjArray> objGenSigArray(sigGenNamesStr.Tokenize(","));
    for (int isig = 0; isig < objGenSigArray->GetEntries(); isig++) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objGenSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() == 1) { // NOTE: 1-prong signals required
          fGenMCSignals.push_back(*sig);
          histNames += Form("MCTruthGen_%s;", sig->GetName()); // TODO: Add these names to a std::vector to avoid using Form in the process function
        }
      }
    }

    DefineHistograms(fHistMan, histNames.Data());    // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    TString sigNamesStr = fConfigMCRecSignals.value;
    std::unique_ptr<TObjArray> objRecSigArray(sigNamesStr.Tokenize(","));
    if (!sigNamesStr.IsNull()) {
      for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
        if (sig) {
          if (sig->GetNProngs() == 1) {
            fRecMCSignals.push_back(*sig);
            fRecMCSignalsNames.push_back(sig->GetName());
          }
        }
      }
    }
    // Setting the MC rec signal names
    for (int isig = 0; isig < objRecSigArray->GetEntries(); ++isig) {
      MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objRecSigArray->At(isig)->GetName());
      if (sig) {
        if (sig->GetNProngs() != 1) { // NOTE: 2-prong signals required
          continue;
        }
        fRecMCSignals.push_back(*sig);
      }
    }
  }

  // Template function to pair mft tracks and muon tracks
  template <bool TMatchedOnly, uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename Muons, typename MftTracks>
  void runFwdTrackPid(TEvent const& event, Muons const& muons, MftTracks const& mftTracks)
  {
    fwdPidAllList.reserve(1);
    for (const auto& muon : muons) {
      if (muon.has_matchMFTTrack() && muon.trackType() == 0 && TMath::Abs(muon.fwdDcaX()) < fConfigMaxDCA && TMath::Abs(muon.fwdDcaY()) < fConfigMaxDCA) {
        auto mftTrack = muon.template matchMFTTrack_as<MyMftTracks>();
        fwdPidAllList(muon.trackType(), event.posX(), event.posY(), event.posZ(), event.numContrib(), muon.pt(), muon.eta(), muon.phi(), muon.sign(), mftTrack.mftClusterSizesAndTrackFlags(), muon.fwdDcaX(), muon.fwdDcaY(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(), 0);
      }
    }
    if constexpr (TMatchedOnly == false) {
      for (const auto& mftTrack : mftTracks) {
        if (TMath::Abs(mftTrack.fwdDcaX()) < fConfigMaxDCA && TMath::Abs(mftTrack.fwdDcaY()) < fConfigMaxDCA) {
          if (downSampleFactor < 1.) {
            float pseudoRndm = mftTrack.pt() * 1000. - (int64_t)(mftTrack.pt() * 1000);
            if (pseudoRndm >= downSampleFactor) {
              continue;
            }
          }
          fwdPidAllList(4, event.posX(), event.posY(), event.posZ(), event.numContrib(), mftTrack.pt(), mftTrack.eta(), mftTrack.phi(), mftTrack.sign(), mftTrack.mftClusterSizesAndTrackFlags(), mftTrack.fwdDcaX(), mftTrack.fwdDcaY(), -999, -999, 0);
        }
      }
    }
  }

  // Template function to run over reconstructed tracks
  template <bool TMatchedOnly, uint32_t TEventFillMap, uint32_t TEventMCFillMap, uint32_t TTrackFillMap, typename TEvent, typename Muons, typename MftTracks, typename TEventsMC, typename TTracksMC>
  void runFwdTrackPidMC(TEvent const& event, Muons const& muons, MftTracks const& mftTracks, TEventsMC const& /*eventsMC*/, TTracksMC const& /*tracksMC*/)
  {
    fwdPidAllList.reserve(1);
    for (const auto& muon : muons) {
      if (muon.has_matchMFTTrack() && muon.trackType() == 0 && TMath::Abs(muon.fwdDcaX()) < fConfigMaxDCA && TMath::Abs(muon.fwdDcaY()) < fConfigMaxDCA) {
        auto mftTrack = muon.template matchMFTTrack_as<MyMftTracksMC>();
        fwdPidAllList(muon.trackType(), event.posX(), event.posY(), event.posZ(), event.numContrib(), muon.pt(), muon.eta(), muon.phi(), muon.sign(), mftTrack.mftClusterSizesAndTrackFlags(), muon.fwdDcaX(), muon.fwdDcaY(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT(), muon.mcReducedFlags());
      }
    }

    if constexpr (TMatchedOnly == false) {
      for (const auto& mftTrack : mftTracks) {
        if (TMath::Abs(mftTrack.fwdDcaX()) < fConfigMaxDCA && TMath::Abs(mftTrack.fwdDcaY()) < fConfigMaxDCA) {
          fwdPidAllList(4, event.posX(), event.posY(), event.posZ(), event.numContrib(), mftTrack.pt(), mftTrack.eta(), mftTrack.phi(), mftTrack.sign(), mftTrack.mftClusterSizesAndTrackFlags(), mftTrack.fwdDcaX(), mftTrack.fwdDcaY(), -999, -999, mftTrack.mcReducedFlags());
        }
      }
    }
  }

  // Template function to run over MC tracks
  template <typename TTracksMC>
  void runMCGen(TTracksMC& groupedMCTracks)
  {
    for (auto& mctrack : groupedMCTracks) {
      VarManager::FillTrackMC(groupedMCTracks, mctrack);

      int isig = 0;
      for (auto sig = fGenMCSignals.begin(); sig != fGenMCSignals.end(); sig++, isig++) {
        if (mctrack.mcReducedFlags() & (static_cast<uint16_t>(1) << isig)) {
          fHistMan->FillHistClass(Form("MCTruthGen_%s", sig->GetName()), VarManager::fgValues);
        }
      }
    }
  }

  PresliceUnsorted<ReducedMCTracks> perReducedMcEvent = aod::reducedtrackMC::reducedMCeventId;

  void processFwdPidMatched(MyEvents::iterator const& event, MyMuonTracks const& muons, MyMftTracks const& mftTracks)
  {
    if (muons.size() > 0 && mftTracks.size() > 0) {
      runFwdTrackPid<false, gkEventFillMap, gkMuonFillMap>(event, muons, mftTracks);
    }
  }

  void processFwdPidMatchedOnly(MyEvents::iterator const& event, MyMuonTracks const& muons, MyMftTracks const& mftTracks)
  {
    if (muons.size() > 0) {
      runFwdTrackPid<true, gkEventFillMap, gkMuonFillMap>(event, muons, mftTracks);
    }
  }

  void processFwdPidMatchedMC(MyEventsMC::iterator const& event, MyMuonTracksMC const& muons, MyMftTracksMC const& mftTracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    if (muons.size() > 0 && mftTracks.size() > 0) {
      runFwdTrackPidMC<false, gkEventFillMap, gkMCEventFillMap, gkMuonFillMap>(event, muons, mftTracks, eventsMC, tracksMC);
    }
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }
  void processFwdPidMatchedOnlyMC(MyEventsMC::iterator const& event, MyMuonTracksMC const& muons, MyMftTracksMC const& mftTracks, ReducedMCEvents const& eventsMC, ReducedMCTracks const& tracksMC)
  {
    if (muons.size() > 0) {
      runFwdTrackPidMC<true, gkEventFillMap, gkMCEventFillMap, gkMuonFillMap>(event, muons, mftTracks, eventsMC, tracksMC);
    }
    auto groupedMCTracks = tracksMC.sliceBy(perReducedMcEvent, event.reducedMCevent().globalIndex());
    groupedMCTracks.bindInternalIndicesTo(&tracksMC);
    runMCGen(groupedMCTracks);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(taskFwdTrackPid, processFwdPidMatched, "Run MFT - muon track pairing filling tree with MFT and global tracks", false);
  PROCESS_SWITCH(taskFwdTrackPid, processFwdPidMatchedOnly, "Run MFT - muon track pairing filling tree with global tracks only", false);
  PROCESS_SWITCH(taskFwdTrackPid, processFwdPidMatchedMC, "Run MFT - muon track pairing filling tree with MFT and global tracks and MC info", false);
  PROCESS_SWITCH(taskFwdTrackPid, processFwdPidMatchedOnlyMC, "Run MFT - muon track pairing filling tree with global tracks only and MC info", false);
  PROCESS_SWITCH(taskFwdTrackPid, processDummy, "Dummy function", false);
}; // End of struct taskFwdTrackPid

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskFwdTrackPid>(cfgc)};
}

void DefineHistograms(HistogramManager* histMan, TString histClasses)
{
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("MCTruthGen")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "mctruth");
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Pt_Rapidity", "MC generator p_{T}, y distribution", false, 120, 0.0, 30.0, VarManager::kMCPt, 150, 2.5, 4.0, VarManager::kMCY);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Eta", "MC generator #eta distribution", false, 200, -5.0, 5.0, VarManager::kMCEta);
      histMan->AddHistogram(objArray->At(iclass)->GetName(), "Phi", "MC generator #varphi distribution", false, 50, 0.0, 2. * TMath::Pi(), VarManager::kMCPhi);
    }

  } // end loop over histogram classes
}
