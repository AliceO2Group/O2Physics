// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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

#include <iostream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TH3F.h>
#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "Common/CCDB/EventSelectionParams.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMftTracks = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

struct taskFwdTrackPid {
  Produces<aod::FwdPidsAll> fwdPidAllList;

  Configurable<float> fConfigMaxDCA{"cfgMaxDCA", 0.5f, "Manually set maximum DCA of the track"};
  Configurable<float> downSampleFactor{"downSampleFactor", 1., "Fraction of candidates to keep for ML"};

  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processDummy")) {
      return;
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
        fwdPidAllList(muon.trackType(), event.posX(), event.posY(), event.posZ(), event.numContrib(), muon.pt(), muon.eta(), muon.phi(), muon.sign(), mftTrack.mftClusterSizesAndTrackFlags(), muon.fwdDcaX(), muon.fwdDcaY(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT());
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
          fwdPidAllList(4, event.posX(), event.posY(), event.posZ(), event.numContrib(), mftTrack.pt(), mftTrack.eta(), mftTrack.phi(), mftTrack.sign(), mftTrack.mftClusterSizesAndTrackFlags(), mftTrack.fwdDcaX(), mftTrack.fwdDcaY(), -999, -999);
        }
      }
    }
  }

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

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(taskFwdTrackPid, processFwdPidMatched, "Run MFT - muon track pairing filling tree with MFT and global tracks", false);
  PROCESS_SWITCH(taskFwdTrackPid, processFwdPidMatchedOnly, "Run MFT - muon track pairing filling tree with global tracks only", false);
  PROCESS_SWITCH(taskFwdTrackPid, processDummy, "Dummy function", false);
}; // End of struct taskFwdTrackPid

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskFwdTrackPid>(cfgc)};
}
