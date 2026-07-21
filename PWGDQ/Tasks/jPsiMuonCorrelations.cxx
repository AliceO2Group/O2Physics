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
/// \file jPsiMuonCorrelations.cxx
/// \brief Task to perform correlations between J/Psi candidates and single muon tracks
/// \author Kaare Endrup Iversen <kaare.endrup.iversen@cern.ch> Lund University
//
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
DECLARE_SOA_BITMAP_COLUMN(IsEventSelected, isEventSelected, 8); //! Event decision
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);  //! Muon track decisions (joinable to ReducedMuonsAssoc)
} // namespace dqanalysisflags

DECLARE_SOA_TABLE(EventCuts, "AOD", "DQANAEVCUTSA", dqanalysisflags::IsEventSelected);      //!  joinable to ReducedEvents
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTSA", dqanalysisflags::IsMuonSelected); //!  joinable to ReducedMuonsAssoc
} // namespace o2::aod

// Declarations of various short names
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;

// using MyPairCandidatesSelected = soa::Join<aod::Dimuons, aod::DimuonsExtra>;
using MyPairCandidatesSelected = soa::Join<aod::Dimuons, aod::DimuonsExtra, aod::DimuonsAll>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonAssocsSelected = soa::Join<aod::ReducedMuonsAssoc, aod::MuonTrackCuts>;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t GkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t GkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

// Declare helper function
double getRapidity(const double pT, const double eta);
double getWeight(const double pT, const std::vector<double>& binsPT, const std::vector<double>& efficiency, const double etaMin, const double etaMax);

struct DqJPsiMuonCorrelations {

  // Configurables for the dilepton signal region
  Configurable<float> fConfigDileptonLowMass{"fConfigDileptonLowMass", 2.8, "Low mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigDileptonHighMass{"fConfigDileptonHighMass", 3.4, "High mass cut for the dileptons used in analysis"};
  Configurable<float> fConfigBackgroundLowMass{"fConfigBackgroundLowMass", 2.5, "Low mass cut for the background used in analysis"};
  Configurable<float> fConfigBackgroundHighMass{"fConfigBackgroundHighMass", 3.7, "High mass cut for the background used in analysis"};

  // Configurables for the dilepton and associated muon cuts
  Configurable<float> fConfigDileptonPtMin{"fConfigDileptonPtMin", 1.0, "Minimum pT cut for the dilepton"};
  Configurable<float> fConfigDileptonPtMax{"fConfigDileptonPtMax", 20.0, "Maximum pT cut for the dilepton"};
  Configurable<float> fConfigDileptonEtaMin{"fConfigDileptonEtaMin", -4.0, "Minimum eta cut for the dileptons"};
  Configurable<float> fConfigDileptonEtaMax{"fConfigDileptonEtaMax", -2.5, "Maximum eta cut for the dileptons"};
  Configurable<float> fConfigMuonEtaMin{"fConfigMuonEtaMin", -4.0, "Minimum eta cut for the associated muons"};
  Configurable<float> fConfigMuonEtaMax{"fConfigMuonEtaMax", -2.5, "Maximum eta cut for the associated muons"};

  // Configurables for histograms
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisInvMass{"axisInvMass", {80, 1.0f, 5.0f}, "Invariant Mass (GeV/c^{2})"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {10, -constants::math::PIHalf, 3.0f * constants::math::PIHalf}, "#Delta#phi (rad)"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {10, -2.0f, 2.0f}, "#Delta#eta"};

  // Configurable for acceptance efficiency correction
  Configurable<std::vector<double>> fConfigBinEffJPsi{"fConfigBinEffJPsi", std::vector<double>{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "acceptance efficiency correction factors for each pT bin"};
  Configurable<std::vector<double>> fConfigBinEffMuon{"fConfigBinEffMuon", std::vector<double>{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}, "acceptance efficiency correction factors for each pT bin"};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb{};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define the filter for events
  Filter eventFilter = aod::dqanalysisflags::isEventSelected == 1;

  // Define the filter for the dileptons
  Filter dileptonFilter = aod::reducedpair::sign == 0;

  constexpr static uint32_t FgDimuonsFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::Pair; // fill map

  // use two values array to avoid mixing up the quantities
  std::vector<float> fValuesDilepton;
  std::vector<float> fValuesMuon;

  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    // Assert correct size of the efficiency correction vector
    int nAxisPtBins = axisPt.value.size() - 2;
    if (nAxisPtBins != fConfigBinEffJPsi.value.size() || nAxisPtBins != fConfigBinEffMuon.value.size()) {
      LOGF(fatal, "Configurables axisPt: %zu must have one more value than fConfigBinEffJPsi: %zu and fConfigBinEffMuon: %zu (excluding 'VARIABLE_WIDTH' entry)",
           axisPt.value.size() - 1, fConfigBinEffJPsi.value.size(), fConfigBinEffMuon.value.size());
    }

    fValuesDilepton.resize(VarManager::kNVars, 0.0f);
    fValuesMuon.resize(VarManager::kNVars, 0.0f);
    VarManager::SetDefaultVarNames();

    // Define trigger histograms
    ConfigurableAxis axisTriggerMass{"axisTriggerMass", {VARIABLE_WIDTH, fConfigBackgroundLowMass, fConfigDileptonLowMass, fConfigDileptonHighMass, fConfigBackgroundHighMass}, "Invariant Mass (GeV/c^{2}) for trigger counting"};
    registry.add("h2dDimuonPtInvVsInvMass", "h2dDimuonPtInvVsInvMass", kTH2D, {axisInvMass, axisPt});
    registry.add("h2dTriggersPtInvVsInvMassRegion", "h2dTriggersPtInvVsInvMassRegion", kTH2D, {axisTriggerMass, axisPt});

    // Define histograms for the dilepton-muon correlations
    registry.add("h2dDimuonMuonDeltaEtaVsMuonPtSignal", "h2dDimuonMuonDeltaEtaVsMuonPtSignal", kTH2D, {axisDeltaEta, axisPt});
    registry.add("h2dDimuonMuonDeltaPhiVsMuonPtSignal", "h2dDimuonMuonDeltaPhiVsMuonPtSignal", kTH2D, {axisDeltaPhi, axisPt});
    registry.add("h2dDimuonMuonDeltaEtaVsMuonPtBackground", "h2dDimuonMuonDeltaEtaVsMuonPtBackground", kTH2D, {axisDeltaEta, axisPt});
    registry.add("h2dDimuonMuonDeltaPhiVsMuonPtBackground", "h2dDimuonMuonDeltaPhiVsMuonPtBackground", kTH2D, {axisDeltaPhi, axisPt});

    // QA histograms
    registry.add("hEventPosZMuon", "hEventPosZMuon", kTH1D, {{50, -25, 25}});
  }

  // Template function to run pair - muon combinations
  template <int TCandidateType, uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuonAssocs, typename TMuonTracks, typename TDileptons>
  void runDileptonMuon(TEvent const& event, TMuonAssocs const& assocs, TMuonTracks const& /*tracks*/, TDileptons const& dileptons)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesMuon.data());
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDilepton.data());
    VarManager::FillEvent<TEventFillMap>(event, fValuesMuon.data());
    VarManager::FillEvent<TEventFillMap>(event, fValuesDilepton.data());

    if (!event.isEventSelected_bit(0)) {
      return;
    }

    if (assocs.size() > 0) {
      registry.fill(HIST("hEventPosZMuon"), event.posZ());
    }

    if (dileptons.size() > 0) {

      for (const auto& dilepton : dileptons) {
        VarManager::FillTrack<FgDimuonsFillMap>(dilepton, fValuesDilepton.data());

        // Dilepton kinematic cuts
        if ((dilepton.eta() < fConfigDileptonEtaMin || dilepton.eta() > fConfigDileptonEtaMax) ||
            (dilepton.pt() < fConfigDileptonPtMin || dilepton.pt() > fConfigDileptonPtMax)) {
          continue;
        }
        // Dilepton leg kinematic cuts
        if ((dilepton.eta1() < fConfigMuonEtaMin || dilepton.eta1() > fConfigMuonEtaMax) ||
            (dilepton.pt1() < axisPt.value[1] || dilepton.pt1() > axisPt.value.back()) ||
            (dilepton.eta2() < fConfigMuonEtaMin || dilepton.eta2() > fConfigMuonEtaMax) ||
            (dilepton.pt2() < axisPt.value[1] || dilepton.pt2() > axisPt.value.back())) {
          continue;
        }

        // Fill invariant mass vs pT histogram for the dileptons and for trigger counting
        double weightDilepton = getWeight(dilepton.pt(), axisPt.value, fConfigBinEffJPsi.value, fConfigDileptonEtaMin, fConfigDileptonEtaMax);

        registry.fill(HIST("h2dDimuonPtInvVsInvMass"), dilepton.mass(), dilepton.pt(), weightDilepton);
        registry.fill(HIST("h2dTriggersPtInvVsInvMassRegion"), dilepton.mass(), dilepton.pt(), weightDilepton);

        for (const auto& assoc : assocs) {
          // Check selection bit
          if (!assoc.isMuonSelected_bit(0)) {
            continue;
          }

          // Skip associated muons that are part of the dilepton candidate
          if (dilepton.index0Id() == assoc.reducedmuonId() || dilepton.index1Id() == assoc.reducedmuonId()) {
            continue;
          }

          // Get muon track information
          auto track = assoc.template reducedmuon_as<TMuonTracks>();

          // Muon kinematic cuts
          if ((track.eta() < fConfigMuonEtaMin || track.eta() > fConfigMuonEtaMax) ||
              (track.pt() < axisPt.value[1] || track.pt() > axisPt.value.back())) {
            continue;
          }

          // Compute deltaEta and deltaPhi between the dilepton and the associated muon
          float deltaEta = dilepton.eta() - track.eta();
          float deltaPhi = RecoDecay::constrainAngle(dilepton.phi() - track.phi(), -constants::math::PIHalf);

          // Fill signal and background histograms based on the dilepton mass
          double weightMuon = getWeight(track.pt(), axisPt.value, fConfigBinEffMuon.value, fConfigMuonEtaMin, fConfigMuonEtaMax);

          if (dilepton.mass() > fConfigDileptonLowMass && dilepton.mass() < fConfigDileptonHighMass) {
            registry.fill(HIST("h2dDimuonMuonDeltaEtaVsMuonPtSignal"), deltaEta, track.pt(), weightDilepton * weightMuon);
            registry.fill(HIST("h2dDimuonMuonDeltaPhiVsMuonPtSignal"), deltaPhi, track.pt(), weightDilepton * weightMuon);
          } else if (dilepton.mass() > fConfigBackgroundLowMass && dilepton.mass() < fConfigBackgroundHighMass) {
            registry.fill(HIST("h2dDimuonMuonDeltaEtaVsMuonPtBackground"), deltaEta, track.pt(), weightDilepton * weightMuon);
            registry.fill(HIST("h2dDimuonMuonDeltaPhiVsMuonPtBackground"), deltaPhi, track.pt(), weightDilepton * weightMuon);
          }
        }
      }
    }
  }

  void processSkimmedDimuon(MyEventsSelected::iterator const& event, MyMuonAssocsSelected const& muonassocs, MyMuonTracks const& muontracks, soa::Filtered<MyPairCandidatesSelected> const& dileptons)
  {
    runDileptonMuon<VarManager::kDecayToMuMu, GkEventFillMap, GkMuonFillMap>(event, muonassocs, muontracks, dileptons);
  }
  void processDummy(MyEvents const&)
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

double getRapidity(const double pT, const double eta)
{
  return std::log(
    (std::sqrt(std::pow(o2::constants::physics::MassJPsi, 2) + (std::pow(pT, 2) * std::pow(std::cosh(eta), 2))) + pT * std::sinh(eta)) /
    (std::sqrt(std::pow(o2::constants::physics::MassJPsi, 2) + std::pow(pT, 2))));
}

double getWeight(const double pT, const std::vector<double>& binsPT, const std::vector<double>& efficiency, const double etaMin, const double etaMax)
{

  int binEff = -1;
  for (size_t b = 0; b < binsPT.size() - 1; ++b) {
    // Shift pT index by one to account for the VARIABLE_WIDTH entry in the axis configuration
    if (pT >= binsPT[b + 1] && pT < binsPT[b + 2]) {
      binEff = b;
      break;
    }
  }

  if (binEff == -1) {
    LOGF(warn, "pT value %f is outside the defined pT bins", pT);
    return 0.0;
  } else {
    return 1.0 / (efficiency[binEff] * (getRapidity(pT, etaMax) - getRapidity(pT, etaMin)));
  }
}
