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

/// \file f0phiproxy.cxx
/// \brief task for dipion f0 candidates with a half-momentum kaon proxy and bachelor kaon.
/// \author Sushanta Tripathy

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <Math/Vector4Dfwd.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct F0phiproxy {
  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using Tracks =
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
              aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi,
              aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FourVector = ROOT::Math::PxPyPzMVector;
  Preslice<Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> vertexMax{"vertexMax", 10.f,
                                "Maximum absolute vertex z (cm)"};
  Configurable<bool> requireSel8{"requireSel8", true, "Require sel8"};
  Configurable<bool> requireINELgt0{"requireINELgt0", true,
                                    "Require at least one PV track in |eta|<1"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15f,
                                 "Minimum track pT (GeV/c)"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.8f, "Maximum track |eta|"};
  Configurable<bool> requireGlobalTrack{"requireGlobalTrack", true,
                                        "Require global track selection"};
  Configurable<bool> requirePVContributor{"requirePVContributor", true,
                                          "Require PV contributor"};
  Configurable<float> pionTPC{"pionTPC", 3.f, "Pion TPC nSigma cut"};
  Configurable<float> pionTOF{"pionTOF", 3.f, "Pion TOF veto when matched"};
  Configurable<float> kaonTPC{"kaonTPC", 3.f, "Kaon TPC nSigma cut"};
  Configurable<float> kaonTOF{"kaonTOF", 3.f, "Kaon TOF veto when matched"};
  Configurable<bool> requireTOF{"requireTOF", false, "Require TOF for both species instead of veto-only"};
  Configurable<float> dipionYMax{"dipionYMax", 0.5f, "Maximum dipion |y|"};
  Configurable<float> proxyYMax{"proxyYMax", 0.5f,
                                "Maximum proxy+bachelor |y|"};
  Configurable<float> dipionMassMin{"dipionMassMin", 0.6f,
                                    "Lower dipion mass for proxy construction"};
  Configurable<float> dipionMassMax{"dipionMassMax", 1.4f,
                                    "Upper dipion mass for proxy construction"};
  Configurable<int> mixDepth{"mixDepth", 5, "Previous events per mixing bin; zero disables mixing"};
  Configurable<float> mixZWidth{"mixZWidth", 2.f, "Vertex bin width (cm)"};
  Configurable<float> mixMultWidth{"mixMultWidth", 20.f,
                                   "PV multiplicity bin width"};

  struct TrackCandidate {
    int64_t id;
    int sign;
    double px, py, pz;
    bool hasTOF;
    [[nodiscard]] FourVector vector(double mass) const
    {
      return {px, py, pz, mass};
    }
  };
  struct Event {
    std::vector<TrackCandidate> pions, kaons;
  };

  void init(InitContext const&)
  {
    if (mixDepth < 0 || mixZWidth <= 0.f || mixMultWidth <= 0.f ||
        vertexMax <= 0.f || pionTPC <= 0.f || pionTOF <= 0.f ||
        kaonTPC <= 0.f || kaonTOF <= 0.f || dipionMassMin >= dipionMassMax) {
      throw std::runtime_error("Invalid F0PhiProxy configuration");
    }
    AxisSpec massPi{300, 0.2, 1.7, "m_{#pi#pi} (GeV/c^{2})"};
    AxisSpec massProxy{600, 0.9, 1.5, "m_{K,proxy K} (GeV/c^{2})"};
    AxisSpec pt{100, 0., 10., "p_{T} (GeV/c)"};
    AxisSpec mult{100, 0., 200., "N_{PV}(|#eta|<1)"};
    AxisSpec charge{3, -0.5, 2.5, "dipion: 0=unlike, 1=++, 2=--"};
    AxisSpec bachelorCharge{2, -1.5, 1.5, "bachelor charge"};
    AxisSpec tof{2, -0.5, 1.5, "bachelor has TOF"};
    histos.add("hEventSelection", "Received; selected", HistType::kTH1D,
               {{2, 0.5, 2.5}});
    histos.add("hMixPartners", "Mixing partners;N", HistType::kTH1D,
               {{101, -0.5, 100.5}});
    histos.add("hPionTPC", "Quality-selected tracks;pT;nSigmaTPC pion",
               HistType::kTH2F, {pt, {100, -10., 10.}});
    histos.add("hPionTOF", "TOF-matched quality tracks;pT;nSigmaTOF pion",
               HistType::kTH2F, {pt, {100, -10., 10.}});
    histos.add("hKaonTPC", "Quality-selected tracks;pT;nSigmaTPC kaon",
               HistType::kTH2F, {pt, {100, -10., 10.}});
    histos.add("hKaonTOF", "TOF-matched quality tracks;pT;nSigmaTOF kaon",
               HistType::kTH2F, {pt, {100, -10., 10.}});
    histos.add("hDipionSE", "Dipion SE", HistType::kTHnSparseF,
               {massPi, pt, mult, charge});
    histos.add("hDipionME", "Dipion ME (raw)", HistType::kTHnSparseF,
               {massPi, pt, mult, charge});
    // Keep measured dipion mass: no f0 identification by a mass window alone.
    histos.add("hProxySE", "Dipion plus bachelor proxy SE",
               HistType::kTHnSparseF,
               {massProxy, pt, massPi, mult, charge, bachelorCharge, tof});
    histos.add("hProxyME", "Intact same-event dipion plus mixed bachelor (raw)",
               HistType::kTHnSparseF,
               {massProxy, pt, massPi, mult, charge, bachelorCharge, tof});
    AxisSpec massKK{600, 0.9, 1.5, "m_{KK} (GeV/c^{2})"};
    histos.add("hPhiKK", "Unlike-sign kaon control SE", HistType::kTHnSparseF,
               {massKK, pt, mult});
    histos.add("hPhiKKLikePP", "K+K+ control SE", HistType::kTHnSparseF,
               {massKK, pt, mult});
    histos.add("hPhiKKLikeMM", "K-K- control SE", HistType::kTHnSparseF,
               {massKK, pt, mult});
    histos.add("hPhiKKMixed", "Unlike-sign kaon control ME (raw)", HistType::kTHnSparseF,
               {massKK, pt, mult});
    histos.add("hPhiKKMixedLikePP", "K+K+ control ME (raw)", HistType::kTHnSparseF,
               {massKK, pt, mult});
    histos.add("hPhiKKMixedLikeMM", "K-K- control ME (raw)", HistType::kTHnSparseF,
               {massKK, pt, mult});
  }

  bool passPID(float tpc, float tof, bool matched, float tpcCut,
               float tofCut) const
  {
    return std::isfinite(tpc) && std::abs(tpc) < tpcCut &&
           (!requireTOF.value || matched) &&
           (!matched || (std::isfinite(tof) && std::abs(tof) < tofCut));
  }

  static int chargeCategory(int a, int b)
  {
    return a != b ? 0 : (a > 0 ? 1 : 2);
  }

  void fillProxy(FourVector const& dipion, TrackCandidate const& first,
                 TrackCandidate const& second,
                 std::vector<TrackCandidate> const& kaons, float mult,
                 bool mixed)
  {
    const FourVector constituent(dipion.Px() * 0.5, dipion.Py() * 0.5,
                                 dipion.Pz() * 0.5,
                                 constants::physics::MassKaonCharged);
    const int category = chargeCategory(first.sign, second.sign);
    for (auto const& kaon : kaons) {
      if (!mixed && (kaon.id == first.id || kaon.id == second.id)) {
        continue;
      }
      const auto pair =
        constituent + kaon.vector(constants::physics::MassKaonCharged);
      if (std::abs(pair.Rapidity()) >= proxyYMax) {
        continue;
      }
      // The neutral dipion has no measured constituent-kaon charge: fill once
      // per bachelor.
      if (mixed) {
        histos.fill(HIST("hProxyME"), pair.M(), pair.Pt(), dipion.M(), mult,
                    category, kaon.sign, kaon.hasTOF);
      } else {
        histos.fill(HIST("hProxySE"), pair.M(), pair.Pt(), dipion.M(), mult,
                    category, kaon.sign, kaon.hasTOF);
      }
    }
  }

  void fillPhiPair(TrackCandidate const& first, TrackCandidate const& second, float mult, bool mixed)
  {
    const auto pair = first.vector(constants::physics::MassKaonCharged) +
                      second.vector(constants::physics::MassKaonCharged);
    if (std::abs(pair.Rapidity()) >= proxyYMax.value) {
      return;
    }
    if (mixed) {
      if (first.sign != second.sign) {
        histos.fill(HIST("hPhiKKMixed"), pair.M(), pair.Pt(), mult);
      } else if (first.sign > 0) {
        histos.fill(HIST("hPhiKKMixedLikePP"), pair.M(), pair.Pt(), mult);
      } else {
        histos.fill(HIST("hPhiKKMixedLikeMM"), pair.M(), pair.Pt(), mult);
      }
    } else {
      if (first.sign != second.sign) {
        histos.fill(HIST("hPhiKK"), pair.M(), pair.Pt(), mult);
      } else if (first.sign > 0) {
        histos.fill(HIST("hPhiKKLikePP"), pair.M(), pair.Pt(), mult);
      } else {
        histos.fill(HIST("hPhiKKLikeMM"), pair.M(), pair.Pt(), mult);
      }
    }
  }

  void fillPhiControls(Event const& event, std::deque<Event> const& pool, float mult)
  {
    for (size_t i = 0; i < event.kaons.size(); ++i) {
      // Each same-event unordered pair is counted once.
      for (size_t j = i + 1; j < event.kaons.size(); ++j) {
        fillPhiPair(event.kaons[i], event.kaons[j], mult, false);
      }
      // Includes K+current K-previous and K-current K+previous once each.
      // Pool insertion occurs after this call; no event is mixed with itself.
      for (auto const& previous : pool) {
        for (auto const& second : previous.kaons) {
          fillPhiPair(event.kaons[i], second, mult, true);
        }
      }
    }
  }

  void process(Collisions const& collisions, Tracks const& tracks,
               aod::BCsWithTimestamps const&)
  {
    // Local pools cannot retain track IDs across dataframes. Run is part of the
    // key.
    std::map<std::pair<int, std::pair<int, int>>, std::deque<Event>> pools;
    for (auto const& collision : collisions) {
      histos.fill(HIST("hEventSelection"), 1.);
      if (std::abs(collision.posZ()) >= vertexMax ||
          (requireSel8 && !collision.sel8()) ||
          (requireINELgt0 && collision.multNTracksPVeta1() < 1)) {
        continue;
      }
      histos.fill(HIST("hEventSelection"), 2.);
      const float mult = collision.multNTracksPVeta1();
      const int run = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
      auto& pool =
        pools[{run,
               {static_cast<int>(std::floor(collision.posZ() / mixZWidth)),
                static_cast<int>(std::floor(mult / mixMultWidth))}}];
      Event event;
      auto selected = tracks.sliceBy(perCollision, collision.globalIndex());
      for (auto const& track : selected) {
        if (track.sign() == 0 || track.pt() < trackPtMin ||
            std::abs(track.eta()) >= trackEtaMax ||
            (requireGlobalTrack && !track.isGlobalTrack()) ||
            (requirePVContributor && !track.isPVContributor())) {
          continue;
        }
        const bool tof = track.hasTOF();
        histos.fill(HIST("hPionTPC"), track.pt(), track.tpcNSigmaPi());
        histos.fill(HIST("hKaonTPC"), track.pt(), track.tpcNSigmaKa());
        if (tof) {
          histos.fill(HIST("hPionTOF"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("hKaonTOF"), track.pt(), track.tofNSigmaKa());
        }
        TrackCandidate candidate{.id = track.globalIndex(), .sign = track.sign(),
                                 .px = track.px(), .py = track.py(), .pz = track.pz(), .hasTOF = tof};
        if (passPID(track.tpcNSigmaPi(), track.tofNSigmaPi(), tof, pionTPC,
                    pionTOF)) {
          event.pions.push_back(candidate);
        }
        if (passPID(track.tpcNSigmaKa(), track.tofNSigmaKa(), tof, kaonTPC,
                    kaonTOF)) {
          event.kaons.push_back(candidate);
        }
      }
      histos.fill(HIST("hMixPartners"), pool.size());
      for (size_t i = 0; i < event.pions.size(); ++i) {
        auto const& first = event.pions[i];
        for (size_t j = i + 1; j < event.pions.size(); ++j) {
          auto const& second = event.pions[j];
          const auto pair = first.vector(constants::physics::MassPionCharged) +
                            second.vector(constants::physics::MassPionCharged);
          if (std::abs(pair.Rapidity()) >= dipionYMax) {
            continue;
          }
          histos.fill(HIST("hDipionSE"), pair.M(), pair.Pt(), mult,
                      chargeCategory(first.sign, second.sign));
          if (pair.M() < dipionMassMin || pair.M() >= dipionMassMax) {
            continue;
          }
          fillProxy(pair, first, second, event.kaons, mult, false);
          for (auto const& previous : pool) {
            fillProxy(pair, first, second, previous.kaons, mult, true);
          }
        }
        for (auto const& previous : pool) {
          for (auto const& second : previous.pions) {
            const auto pair =
              first.vector(constants::physics::MassPionCharged) +
              second.vector(constants::physics::MassPionCharged);
            if (std::abs(pair.Rapidity()) < dipionYMax) {
              histos.fill(HIST("hDipionME"), pair.M(), pair.Pt(), mult,
                          chargeCategory(first.sign, second.sign));
            }
          }
        }
      }
      fillPhiControls(event, pool, mult);
      if (mixDepth > 0) {
        pool.push_back(std::move(event));
        if (pool.size() > static_cast<size_t>(mixDepth.value)) {
          pool.pop_front();
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<F0phiproxy>(context)};
}
