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

/// \file taskNetCharmFluctuations.cxx
/// \brief Producer of per-candidate and per-event charm-counting tables for net-charm fluctuation studies
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University
/// \author Fan Si <fsi@physi.uni-heidelberg.de>, Heidelberg University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;

namespace o2::aod
{
namespace eyefluc
{
DECLARE_SOA_COLUMN(EventId, eventId, uint64_t);
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t);
DECLARE_SOA_COLUMN(CandUid, candUid, uint64_t);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);
DECLARE_SOA_COLUMN(MassD0, massD0, float);
DECLARE_SOA_COLUMN(MassD0bar, massD0bar, float);
DECLARE_SOA_COLUMN(MassDplus, massDplus, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(FitBinId, fitBinId, int16_t);
DECLARE_SOA_COLUMN(OmegaCharm, omegaCharm, float);
DECLARE_SOA_COLUMN(OmegaAntiCharm, omegaAntiCharm, float);
DECLARE_SOA_COLUMN(OmegaBkg, omegaBkg, float);
DECLARE_SOA_COLUMN(WCharm, wCharm, double);
DECLARE_SOA_COLUMN(WAntiCharm, wAntiCharm, double);
DECLARE_SOA_COLUMN(WBkg, wBkg, double);
} // namespace eyefluc

DECLARE_SOA_TABLE(EyeFlucCharmD0Cands, "AOD", "EYEFCD0CAND",
                  eyefluc::EventId,
                  eyefluc::TimeStamp,
                  eyefluc::CandUid,
                  eyefluc::Sign,
                  eyefluc::Pt,
                  eyefluc::Rapidity,
                  eyefluc::MassD0,
                  eyefluc::MassD0bar,
                  eyefluc::FitBinId,
                  eyefluc::OmegaCharm,
                  eyefluc::OmegaAntiCharm,
                  eyefluc::OmegaBkg);

DECLARE_SOA_TABLE(EyeFlucCharmD0Events, "AOD", "EYEFCD0EVT",
                  eyefluc::EventId,
                  eyefluc::TimeStamp,
                  eyefluc::Centrality,
                  eyefluc::WCharm,
                  eyefluc::WAntiCharm,
                  eyefluc::WBkg);

DECLARE_SOA_TABLE(EyeFlucCharmDplusCands, "AOD", "EYEFCDPCAND",
                  eyefluc::EventId,
                  eyefluc::TimeStamp,
                  eyefluc::CandUid,
                  eyefluc::Sign,
                  eyefluc::Pt,
                  eyefluc::Rapidity,
                  eyefluc::MassDplus,
                  eyefluc::FitBinId,
                  eyefluc::OmegaCharm,
                  eyefluc::OmegaAntiCharm,
                  eyefluc::OmegaBkg);

DECLARE_SOA_TABLE(EyeFlucCharmDplusEvents, "AOD", "EYEFCDPEVT",
                  eyefluc::EventId,
                  eyefluc::TimeStamp,
                  eyefluc::Centrality,
                  eyefluc::WCharm,
                  eyefluc::WAntiCharm,
                  eyefluc::WBkg);
} // namespace o2::aod

enum class CharmFamily : uint8_t {
  D0 = 0,
  Dplus = 1
};

enum EventQa : uint8_t {
  All = 0,
  RejHfEventSelection,
  RejNoCharmCandidate,
  CharmCandidateSelected,
  EventWritten,
  NEventQa
};

using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong,
                                           aod::HfCand2Prong0PidPi,
                                           aod::HfCand2Prong1PidPi,
                                           aod::HfCand2Prong0PidKa,
                                           aod::HfCand2Prong1PidKa,
                                           aod::HfCand2ProngKF,
                                           aod::HfSelD0>>;

using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong,
                                              aod::HfCand3Prong0PidPi,
                                              aod::HfCand3Prong1PidPi,
                                              aod::HfCand3Prong2PidPi,
                                              aod::HfCand3Prong0PidKa,
                                              aod::HfCand3Prong1PidKa,
                                              aod::HfCand3Prong2PidKa,
                                              aod::HfSelDplusToPiKPi>>;

using CollData = soa::Join<aod::Collisions,
                           aod::EvSels,
                           aod::CentFT0Cs,
                           aod::CentFT0Ms,
                           aod::CentFT0As>;

struct HfTaskNetCharmFluctuations {
  Produces<aod::EyeFlucCharmD0Cands> outD0Cand;
  Produces<aod::EyeFlucCharmD0Events> outD0Evt;
  Produces<aod::EyeFlucCharmDplusCands> outDplusCand;
  Produces<aod::EyeFlucCharmDplusEvents> outDplusEvt;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Minimum D0 selection flag"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Minimum D0bar selection flag"};
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Minimum Dplus selection flag"};
  Configurable<int> centralityEstimator{"centralityEstimator", static_cast<int>(CentralityEstimator::FT0C), "Centrality estimator used for output tables"};
  Configurable<std::vector<float>> ptFitBins{"ptFitBins", std::vector<float>{1.f, 2.f, 3.f, 4.f, 6.f, 8.f, 12.f, 24.f}, "pT bins used to assign fitBinId"};
  Configurable<bool> fillOmegaRaw{"fillOmegaRaw", true, "Fill omega sums with raw charm/anti-charm candidate counts"};

  SliceCache cache;
  HfEventSelection hfEvSel;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<CandD0Data> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0;
  Partition<CandD0Data> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  HistogramRegistry registry{"registry"};

  struct HfCandInfo {
    uint64_t uid = 0;
    uint8_t family = 0;
    int8_t sign = 0;
    float massD0 = -1.f;
    float massD0bar = -1.f;
    float massDplus = -1.f;
    float pt = -1.f;
    float rapidity = -999.f;
    float omegaCharm = 0.f;
    float omegaAntiCharm = 0.f;
    float omegaBkg = 1.f;
  };

  void init(InitContext const&)
  {
    std::array<int, 2> processes = {doprocessD0, doprocessDplus};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    static constexpr std::array<std::string_view, EventQa::NEventQa> EventLabels = {
      "All events",
      "rejected by HF event selection",
      "without charm candidates",
      "with charm candidates",
      "written events"};
    static constexpr double EventQaAxisMin = 0.5;
    static constexpr double EventQaAxisMax = static_cast<double>(EventQa::NEventQa) + EventQaAxisMin;
    registry.add("hEventQa", "Event QA;;entries", {HistType::kTH1F, {{static_cast<int>(EventQa::NEventQa), EventQaAxisMin, EventQaAxisMax}}});
    for (int iBin = 0; iBin < EventQa::NEventQa; ++iBin) {
      registry.get<TH1>(HIST("hEventQa"))->GetXaxis()->SetBinLabel(iBin + 1, EventLabels[iBin].data());
    }
    registry.add("hMassVsPtD0", "D0 candidates;#it{M}_{#pi K} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{500, 1.65, 2.15}, {200, 0., 50.}}});
    registry.add("hMassVsPtD0bar", "D0bar candidates;#it{M}_{K#pi} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{500, 1.65, 2.15}, {200, 0., 50.}}});
    registry.add("hMassVsPtDplus", "Dplus candidates;#it{M}_{#pi K#pi} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{500, 1.65, 2.15}, {200, 0., 50.}}});
    registry.add("hCandidateCounter", "Candidate counter;family/sign;entries", {HistType::kTH1F, {{4, 0.5, 4.5}}});
    hfEvSel.init(registry);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  uint64_t makeEventId(CollData::iterator const& collision) const
  {
    return static_cast<uint64_t>(collision.globalIndex());
  }

  int64_t getTimeStamp(CollData::iterator const& collision) const
  {
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    return bc.timestamp();
  }

  float getCentrality(CollData::iterator const& collision) const
  {
    return getCentralityColl(collision, centralityEstimator.value);
  }

  uint64_t makeCandUid(CharmFamily family, int8_t sign, int64_t globalIndex) const
  {
    const uint64_t signBit = sign > 0 ? 0ull : 1ull;
    return (static_cast<uint64_t>(family) << 62) | (signBit << 61) | static_cast<uint64_t>(globalIndex);
  }

  int16_t getFitBin(float pt) const
  {
    auto const& bins = ptFitBins.value;
    static constexpr size_t MinFitBinEdges = 2;
    if (bins.size() < MinFitBinEdges) {
      return -1;
    }
    for (size_t iBin = 0; iBin + 1 < bins.size(); ++iBin) {
      if (pt >= bins[iBin] && pt < bins[iBin + 1]) {
        return static_cast<int16_t>(iBin);
      }
    }
    return -1;
  }

  void setOmegaRaw(HfCandInfo& cand) const
  {
    if (!fillOmegaRaw.value) {
      return;
    }
    cand.omegaBkg = 0.f;
    if (cand.sign > 0) {
      cand.omegaCharm = 1.f;
    } else {
      cand.omegaAntiCharm = 1.f;
    }
  }

  bool passEventSelection(CollData::iterator const& collision)
  {
    registry.fill(HIST("hEventQa"), 1 + EventQa::All);
    float centrality = 0.f;
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    hfEvSel.fillHistograms(collision, rejectionMask, centrality);
    if (rejectionMask != 0) {
      registry.fill(HIST("hEventQa"), 1 + EventQa::RejHfEventSelection);
      return false;
    }
    return true;
  }

  void fillD0OutputTables(CollData::iterator const& collision, std::vector<HfCandInfo>& acceptedCands)
  {
    const uint64_t eventId = makeEventId(collision);
    const int64_t timeStamp = getTimeStamp(collision);
    const float centrality = getCentrality(collision);

    if (acceptedCands.empty()) {
      outD0Evt(eventId, timeStamp, centrality, 0., 0., 0.);
      registry.fill(HIST("hEventQa"), 1 + EventQa::RejNoCharmCandidate);
      registry.fill(HIST("hEventQa"), 1 + EventQa::EventWritten);
      return;
    }
    registry.fill(HIST("hEventQa"), 1 + EventQa::CharmCandidateSelected);

    double wCharm = 0.;
    double wAntiCharm = 0.;
    double wBkg = 0.;

    for (const auto& cand : acceptedCands) {
      wCharm += cand.omegaCharm;
      wAntiCharm += cand.omegaAntiCharm;
      wBkg += cand.omegaBkg;

      outD0Cand(eventId,
                timeStamp,
                cand.uid,
                cand.sign,
                cand.pt,
                cand.rapidity,
                cand.massD0,
                cand.massD0bar,
                getFitBin(cand.pt),
                cand.omegaCharm,
                cand.omegaAntiCharm,
                cand.omegaBkg);
    }

    outD0Evt(eventId, timeStamp, centrality, wCharm, wAntiCharm, wBkg);
    registry.fill(HIST("hEventQa"), 1 + EventQa::EventWritten);
  }

  void fillDplusOutputTables(CollData::iterator const& collision, std::vector<HfCandInfo>& acceptedCands)
  {
    const uint64_t eventId = makeEventId(collision);
    const int64_t timeStamp = getTimeStamp(collision);
    const float centrality = getCentrality(collision);

    if (acceptedCands.empty()) {
      outDplusEvt(eventId, timeStamp, centrality, 0., 0., 0.);
      registry.fill(HIST("hEventQa"), 1 + EventQa::RejNoCharmCandidate);
      registry.fill(HIST("hEventQa"), 1 + EventQa::EventWritten);
      return;
    }
    registry.fill(HIST("hEventQa"), 1 + EventQa::CharmCandidateSelected);

    double wCharm = 0.;
    double wAntiCharm = 0.;
    double wBkg = 0.;

    for (const auto& cand : acceptedCands) {
      wCharm += cand.omegaCharm;
      wAntiCharm += cand.omegaAntiCharm;
      wBkg += cand.omegaBkg;

      outDplusCand(eventId,
                   timeStamp,
                   cand.uid,
                   cand.sign,
                   cand.pt,
                   cand.rapidity,
                   cand.massDplus,
                   getFitBin(cand.pt),
                   cand.omegaCharm,
                   cand.omegaAntiCharm,
                   cand.omegaBkg);
    }

    outDplusEvt(eventId, timeStamp, centrality, wCharm, wAntiCharm, wBkg);
    registry.fill(HIST("hEventQa"), 1 + EventQa::EventWritten);
  }

  template <int8_t Sign, typename TCandidates>
  void addD0Candidates(TCandidates const& candidates, std::vector<HfCandInfo>& acceptedCands)
  {
    for (const auto& cand : candidates) {
      const float massD0 = HfHelper::invMassD0ToPiK(cand);
      const float massD0bar = HfHelper::invMassD0barToKPi(cand);

      HfCandInfo info;
      info.uid = makeCandUid(CharmFamily::D0, Sign, cand.globalIndex());
      info.family = static_cast<uint8_t>(CharmFamily::D0);
      info.sign = Sign;
      info.massD0 = massD0;
      info.massD0bar = massD0bar;
      info.pt = cand.pt();
      info.rapidity = HfHelper::yD0(cand);
      setOmegaRaw(info);
      acceptedCands.push_back(info);

      if constexpr (Sign > 0) {
        registry.fill(HIST("hMassVsPtD0"), massD0, cand.pt());
        registry.fill(HIST("hCandidateCounter"), 1.f);
      } else {
        registry.fill(HIST("hMassVsPtD0bar"), massD0bar, cand.pt());
        registry.fill(HIST("hCandidateCounter"), 2.f);
      }
    }
  }

  void processD0(CollData::iterator const& collision,
                 aod::BCsWithTimestamps const&,
                 CandD0Data const&)
  {
    if (!passEventSelection(collision)) {
      return;
    }

    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);

    std::vector<HfCandInfo> acceptedCands;
    addD0Candidates<+1>(candsD0ToPiK, acceptedCands);
    addD0Candidates<-1>(candsD0ToKPi, acceptedCands);
    fillD0OutputTables(collision, acceptedCands);
  }
  PROCESS_SWITCH(HfTaskNetCharmFluctuations, processD0, "Process D0 and D0bar candidates", false);

  void processDplus(CollData::iterator const& collision,
                    aod::BCsWithTimestamps const&,
                    CandDplusData const& candidatesDplus,
                    aod::Tracks const&)
  {
    if (!passEventSelection(collision)) {
      return;
    }

    std::vector<HfCandInfo> acceptedCands;
    for (const auto& cand : candidatesDplus) {
      const auto trackProng0 = cand.template prong0_as<aod::Tracks>();
      const int8_t sign = trackProng0.sign() > 0 ? +1 : -1;
      const float massDplus = HfHelper::invMassDplusToPiKPi(cand);

      HfCandInfo info;
      info.uid = makeCandUid(CharmFamily::Dplus, sign, cand.globalIndex());
      info.family = static_cast<uint8_t>(CharmFamily::Dplus);
      info.sign = sign;
      info.massDplus = massDplus;
      info.pt = cand.pt();
      info.rapidity = HfHelper::yDplus(cand);
      setOmegaRaw(info);
      acceptedCands.push_back(info);
      registry.fill(HIST("hMassVsPtDplus"), massDplus, cand.pt());
      registry.fill(HIST("hCandidateCounter"), sign > 0 ? 3.f : 4.f);
    }

    fillDplusOutputTables(collision, acceptedCands);
  }
  PROCESS_SWITCH(HfTaskNetCharmFluctuations, processDplus, "Process Dplus and Dminus candidates", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskNetCharmFluctuations>(cfgc)};
}
