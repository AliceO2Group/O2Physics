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
/// @author  Kaare Endrup Iversen
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StringHelpers.h>
#include <Framework/runDataProcessing.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{

namespace dqanalysisflags
{
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 8);                      //! Event decision
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                       //! Muon track decisions (joinable to ReducedMuonsAssoc)
}

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTSA", dqanalysisflags::IsEventSelected);                                                            //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTSA", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
} 

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;

using MyPairCandidatesSelected = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonAssocsSelected = soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

struct DqJPsiMuonCorrelations {

  // Configurables for the dilepton and dilepton cuts
  Configurable<float> fConfigDileptonLowMass{"cfgDileptonLowMass", 2.8, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"cfgDileptonHighMass", 3.4, "High mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigBackgroundLowMass{"cfgBackgroundLowMass", 2.5, "Low mass cut for the background used in analysis"};
  Configurable<float> fConfigBackgroundHighMass{"cfgBackgroundHighMass", 3.7, "High mass cut for the background used in analysis"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisInvMass{"axisInvMass", {80, 1.0f, 5.0f}, "Invariant Mass (GeV/c^{2})"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {10, -constants::math::PI/2.0f, 3.0f*constants::math::PI/2.0f}, "#Delta#phi (rad)"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {10, -2.0f, 2.0f}, "#Delta#eta"};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define the filter for events
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;

  // Define the filter for the dileptons
  Filter dileptonFilter = aod::reducedpair::sign == 0;

  constexpr static uint32_t fgDimuonsFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::Pair;   // fill map

  // use two values array to avoid mixing up the quantities
  float* fValuesDilepton;
  float* fValuesMuon;

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    fValuesDilepton = new float[VarManager::kNVars];
    fValuesMuon = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();

    registry.add("h2dDimuonPtInvVsInvMass", "h2dDimuonPtInvVsInvMass", kTH2D, {axisInvMass, axisPt});
    registry.add("h2dDimuonMuonDeltaEtaVsMuonPtSignal", "h2dDimuonMuonDeltaEtaVsMuonPtSignal", kTH2D, {axisDeltaEta, axisPt});
    registry.add("h2dDimuonMuonDeltaPhiVsMuonPtSignal", "h2dDimuonMuonDeltaPhiVsMuonPtSignal", kTH2D, {axisDeltaPhi, axisPt});
    registry.add("h2dDimuonMuonDeltaEtaVsMuonPtBackground", "h2dDimuonMuonDeltaEtaVsMuonPtBackground", kTH2D, {axisDeltaEta, axisPt});
    registry.add("h2dDimuonMuonDeltaPhiVsMuonPtBackground", "h2dDimuonMuonDeltaPhiVsMuonPtBackground", kTH2D, {axisDeltaPhi, axisPt});
  }

  // Template function to run pair - muon combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuonAssocs, typename TMuonTracks, typename TDileptons>
  void runDileptonMuon(TEvent const& event, TMuonAssocs const& assocs, TMuonTracks const& /*tracks*/, TDileptons const& dileptons)
  {
    // VarManager::ResetValues(0, VarManager::kNVars, fValuesHadron);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesMuon);
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton);
    VarManager::FillEvent<TEventFillMap>(event, fValuesMuon);
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton);

    if (!event.isEventSelected_bit(0)) {
      return;
    }

    if (dileptons.size() > 0) {

      for (auto& dilepton : dileptons) {
        VarManager::FillTrack<fgDimuonsFillMap>(dilepton, fValuesDilepton);
        registry.fill(HIST("h2dDimuonPtInvVsInvMass"), dilepton.mass(), dilepton.pt());

        for (auto& assoc : assocs) {
          if (!assoc.isMuonSelected_bit(0)) {
            continue;
          }
          auto track = assoc.template reducedmuon_as<TMuonTracks>();

          float deltaEta = track.eta() - dilepton.eta();
          float deltaPhi = track.phi() - dilepton.phi();
          if (deltaPhi < -constants::math::PI/2.0f) {
            deltaPhi += 2.0f * constants::math::PI;
          } else if (deltaPhi > constants::math::PI*3.0f/2.0f) {
            deltaPhi -= 2.0f * constants::math::PI;
          }

          if (dilepton.mass() > fConfigDileptonLowMass && dilepton.mass() < fConfigDileptonHighMass) {
            registry.fill(HIST("h2dDimuonMuonDeltaEtaVsMuonPtSignal"), deltaEta, track.pt());
            registry.fill(HIST("h2dDimuonMuonDeltaPhiVsMuonPtSignal"), deltaPhi, track.pt());
          } else if (dilepton.mass() > fConfigBackgroundLowMass && dilepton.mass() < fConfigBackgroundHighMass) {
            registry.fill(HIST("h2dDimuonMuonDeltaEtaVsMuonPtBackground"), deltaEta, track.pt());
            registry.fill(HIST("h2dDimuonMuonDeltaPhiVsMuonPtBackground"), deltaPhi, track.pt());
          }
        }

      }
    }
  }

  void processSkimmedDimuon(MyEventsSelected::iterator const& event, MyMuonAssocsSelected const& muonassocs, MyMuonTracks const& muontracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  {
    runDileptonMuon<VarManager::kDecayToMuMu, gkEventFillMap, gkMuonFillMap>(event, muonassocs, muontracks, dileptons);
  }
  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(DqJPsiMuonCorrelations, processSkimmedDimuon, "Run dilepton-muon pairing, using skimmed data", false);
  PROCESS_SWITCH(DqJPsiMuonCorrelations, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DqJPsiMuonCorrelations>(cfgc)};
}
