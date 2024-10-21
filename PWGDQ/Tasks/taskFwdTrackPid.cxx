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

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
// TODO: the barrel amd muon selection columns are bit maps so unsigned types should be used, however, for now this is not supported in Filter expressions
// TODO: For now in the tasks we just statically convert from unsigned int to int, which should be fine as long as we do
//      not use a large number of bits (>=30)
// Bcandidate columns for ML analysis of B->Jpsi+K
DECLARE_SOA_COLUMN(MixingHash, mixingHash, int);
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
DECLARE_SOA_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, int);
DECLARE_SOA_COLUMN(IsPrefilterVetoed, isPrefilterVetoed, int);
DECLARE_SOA_COLUMN(massBcandidate, MBcandidate, float);
DECLARE_SOA_COLUMN(pTBcandidate, PtBcandidate, float);
DECLARE_SOA_COLUMN(LxyBcandidate, lxyBcandidate, float);
DECLARE_SOA_COLUMN(LxyzBcandidate, lxyzBcandidate, float);
DECLARE_SOA_COLUMN(LzBcandidate, lzBcandidate, float);
DECLARE_SOA_COLUMN(TauxyBcandidate, tauxyBcandidate, float);
DECLARE_SOA_COLUMN(TauzBcandidate, tauzBcandidate, float);
DECLARE_SOA_COLUMN(CosPBcandidate, cosPBcandidate, float);
DECLARE_SOA_COLUMN(Chi2Bcandidate, chi2Bcandidate, float);

// Xcandidate columns
DECLARE_SOA_COLUMN(massXcandidate, MXcandidate, float);
DECLARE_SOA_COLUMN(pTXcandidate, PtXcandidate, float);
DECLARE_SOA_COLUMN(rapidityXcandidate, YXcandidate, float);
DECLARE_SOA_COLUMN(etaXcandidate, EtaXcandidate, float);
DECLARE_SOA_COLUMN(massJpsicandidate, MJpsicandidate, float);
DECLARE_SOA_COLUMN(massDipioncandidate, MDipioncandidate, float);
DECLARE_SOA_COLUMN(pTJpsicandidate, PtJpsicandidate, float);
DECLARE_SOA_COLUMN(pTDipioncandidate, PtDipioncandidate, float);
DECLARE_SOA_COLUMN(massDiff, Q, float);
DECLARE_SOA_COLUMN(angDistPion1, DeltaR1, float);
DECLARE_SOA_COLUMN(angDistPion2, DeltaR2, float);
DECLARE_SOA_COLUMN(cosDileptonDipion, CosDileptonDipion, float);
DECLARE_SOA_COLUMN(dcaxyPion1, DcaXYPion1, float);
DECLARE_SOA_COLUMN(dcazPion1, DcaZPion1, float);
DECLARE_SOA_COLUMN(dcaxyPion2, DcaXYPion2, float);
DECLARE_SOA_COLUMN(dcazPion2, DcaZPion2, float);
DECLARE_SOA_COLUMN(pTPion1, PtPion1, float);
DECLARE_SOA_COLUMN(pTPion2, PtPion2, float);
DECLARE_SOA_COLUMN(dileptonSign, DileptonSign, int);
DECLARE_SOA_COLUMN(ditrackSign, DitrackSign, int);

} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTS", dqanalysisflags::IsEventSelected);
DECLARE_SOA_TABLE(MixingHashes, "AOD", "DQANAMIXHASH", dqanalysisflags::MixingHash);
DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "DQANATRKCUTS", dqanalysisflags::IsBarrelSelected, dqanalysisflags::IsBarrelSelectedPrefilter);
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);
DECLARE_SOA_TABLE(Prefilter, "AOD", "DQPREFILTER", dqanalysisflags::IsPrefilterVetoed);
DECLARE_SOA_TABLE(BmesonCandidates, "AOD", "DQBMESONS", dqanalysisflags::massBcandidate, dqanalysisflags::pTBcandidate, dqanalysisflags::LxyBcandidate, dqanalysisflags::LxyzBcandidate, dqanalysisflags::LzBcandidate, dqanalysisflags::TauxyBcandidate, dqanalysisflags::TauzBcandidate, dqanalysisflags::CosPBcandidate, dqanalysisflags::Chi2Bcandidate);
DECLARE_SOA_TABLE(XCandidates, "AOD", "DQX3872", dqanalysisflags::massXcandidate, dqanalysisflags::pTXcandidate, dqanalysisflags::rapidityXcandidate, dqanalysisflags::etaXcandidate, dqanalysisflags::massJpsicandidate, dqanalysisflags::massDipioncandidate, dqanalysisflags::pTJpsicandidate, dqanalysisflags::pTDipioncandidate, dqanalysisflags::massDiff, dqanalysisflags::angDistPion1, dqanalysisflags::angDistPion2, dqanalysisflags::cosDileptonDipion, dqanalysisflags::dcaxyPion1, dqanalysisflags::dcazPion1, dqanalysisflags::dcaxyPion2, dqanalysisflags::dcazPion2, dqanalysisflags::pTPion1, dqanalysisflags::pTPion2, dqanalysisflags::dileptonSign, dqanalysisflags::ditrackSign);
} // namespace o2::aod

using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMftTracks = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

struct taskFwdTrackPid {
  Produces<aod::FwdPidsAll> fwdPidAllList;

  Configurable<float> fConfigMaxDCA{"cfgMaxDCA", 0.5f, "Manually set maximum DCA of the track"};

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
