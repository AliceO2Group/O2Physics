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
// Preliminary QA analysis task for resonances

///
/// \file   lfpidqa.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-12-12
/// \brief  Task to produce the PID QA information for the TPC for the purpose of the Light flavor PWG
///

#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

#include <string_view>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct lfpidqa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hdelta[Np] = {"delta/El", "delta/Mu", "delta/Pi",
                                                  "delta/Ka", "delta/Pr", "delta/De",
                                                  "delta/Tr", "delta/He", "delta/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};

  static constexpr std::string_view hdeltalf[Np] = {"deltalf/El", "deltalf/Mu", "deltalf/Pi",
                                                    "deltalf/Ka", "deltalf/Pr", "deltalf/De",
                                                    "deltalf/Tr", "deltalf/He", "deltalf/Al"};
  static constexpr std::string_view hnsigmalf[Np] = {"nsigmalf/El", "nsigmalf/Mu", "nsigmalf/Pi",
                                                     "nsigmalf/Ka", "nsigmalf/Pr", "nsigmalf/De",
                                                     "nsigmalf/Tr", "nsigmalf/He", "nsigmalf/Al"};
  static constexpr std::string_view hexpectedlf[Np] = {"expectedlf/El", "expectedlf/Mu", "expectedlf/Pi",
                                                       "expectedlf/Ka", "expectedlf/Pr", "expectedlf/De",
                                                       "expectedlf/Tr", "expectedlf/He", "expectedlf/Al"};

  // Positive
  static constexpr std::string_view hdeltapos[Np] = {"delta/pos/El", "delta/pos/Mu", "delta/pos/Pi",
                                                     "delta/pos/Ka", "delta/pos/Pr", "delta/pos/De",
                                                     "delta/pos/Tr", "delta/pos/He", "delta/pos/Al"};
  static constexpr std::string_view hnsigmapos[Np] = {"nsigma/pos/El", "nsigma/pos/Mu", "nsigma/pos/Pi",
                                                      "nsigma/pos/Ka", "nsigma/pos/Pr", "nsigma/pos/De",
                                                      "nsigma/pos/Tr", "nsigma/pos/He", "nsigma/pos/Al"};
  static constexpr std::string_view hexpectedpos[Np] = {"expected/pos/El", "expected/pos/Mu", "expected/pos/Pi",
                                                        "expected/pos/Ka", "expected/pos/Pr", "expected/pos/De",
                                                        "expected/pos/Tr", "expected/pos/He", "expected/pos/Al"};

  static constexpr std::string_view hdeltalfpos[Np] = {"deltalf/pos/El", "deltalf/pos/Mu", "deltalf/pos/Pi",
                                                       "deltalf/pos/Ka", "deltalf/pos/Pr", "deltalf/pos/De",
                                                       "deltalf/pos/Tr", "deltalf/pos/He", "deltalf/pos/Al"};
  static constexpr std::string_view hnsigmalfpos[Np] = {"nsigmalf/pos/El", "nsigmalf/pos/Mu", "nsigmalf/pos/Pi",
                                                        "nsigmalf/pos/Ka", "nsigmalf/pos/Pr", "nsigmalf/pos/De",
                                                        "nsigmalf/pos/Tr", "nsigmalf/pos/He", "nsigmalf/pos/Al"};
  static constexpr std::string_view hexpectedlfpos[Np] = {"expectedlf/pos/El", "expectedlf/pos/Mu", "expectedlf/pos/Pi",
                                                          "expectedlf/pos/Ka", "expectedlf/pos/Pr", "expectedlf/pos/De",
                                                          "expectedlf/pos/Tr", "expectedlf/pos/He", "expectedlf/pos/Al"};

  // Negative
  static constexpr std::string_view hdeltaneg[Np] = {"delta/neg/El", "delta/neg/Mu", "delta/neg/Pi",
                                                     "delta/neg/Ka", "delta/neg/Pr", "delta/neg/De",
                                                     "delta/neg/Tr", "delta/neg/He", "delta/neg/Al"};
  static constexpr std::string_view hnsigmaneg[Np] = {"nsigma/neg/El", "nsigma/neg/Mu", "nsigma/neg/Pi",
                                                      "nsigma/neg/Ka", "nsigma/neg/Pr", "nsigma/neg/De",
                                                      "nsigma/neg/Tr", "nsigma/neg/He", "nsigma/neg/Al"};
  static constexpr std::string_view hexpectedneg[Np] = {"expected/neg/El", "expected/neg/Mu", "expected/neg/Pi",
                                                        "expected/neg/Ka", "expected/neg/Pr", "expected/neg/De",
                                                        "expected/neg/Tr", "expected/neg/He", "expected/neg/Al"};

  static constexpr std::string_view hdeltalfneg[Np] = {"deltalf/neg/El", "deltalf/neg/Mu", "deltalf/neg/Pi",
                                                       "deltalf/neg/Ka", "deltalf/neg/Pr", "deltalf/neg/De",
                                                       "deltalf/neg/Tr", "deltalf/neg/He", "deltalf/neg/Al"};
  static constexpr std::string_view hnsigmalfneg[Np] = {"nsigmalf/neg/El", "nsigmalf/neg/Mu", "nsigmalf/neg/Pi",
                                                        "nsigmalf/neg/Ka", "nsigmalf/neg/Pr", "nsigmalf/neg/De",
                                                        "nsigmalf/neg/Tr", "nsigmalf/neg/He", "nsigmalf/neg/Al"};
  static constexpr std::string_view hexpectedlfneg[Np] = {"expectedlf/neg/El", "expectedlf/neg/Mu", "expectedlf/neg/Pi",
                                                          "expectedlf/neg/Ka", "expectedlf/neg/Pr", "expectedlf/neg/De",
                                                          "expectedlf/neg/Tr", "expectedlf/neg/He", "expectedlf/neg/Al"};

  Configurable<uint16_t> minPVcontrib{"minPVcontrib", 0, "Minimum number of PV contributors"};
  Configurable<uint16_t> maxPVcontrib{"maxPVcontrib", 10000, "Maximum number of PV contributors"};
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum value for |z_vtx|"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> nBinsDelta{"nBinsDelta", 400, "Number of bins for the Delta"};
  Configurable<float> minDelta{"minDelta", -1000.f, "Minimum Delta in range"};
  Configurable<float> maxDelta{"maxDelta", 1000.f, "Maximum Delta in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 401, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.025f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.025f, "Maximum NSigma in range"};
  Configurable<bool> enableTrackCuts{"enableTrackCuts", false, "Flag to enable the custom selections for the QA plots"};
  Configurable<float> minPTPC{"minPTPC", 0.f, "Minimum TPC Momentum for the QA plots"};
  Configurable<float> maxPTPC{"maxPTPC", 0.f, "Maximum TPC Momentum for the QA plots"};
  Configurable<float> minTPCdEdx{"minTPCdEdx", 0.f, "Minimum TPC dEdx for the QA plots"};
  Configurable<float> maxTPCdEdx{"maxTPCdEdx", 0.f, "Maximum TPC dEdx for the QA plots"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    if (doprocessPiKaPr && doprocessAll) {
      LOG(fatal) << "Cannot enable processPiKaPr and processAll at the same time";
    }

    const AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    for (int id = 0; id < Np; id++) {
      if (doprocessPiKaPr) {
        bool keepIt = false;
        switch (id) {
          case o2::track::PID::Pion:
          case o2::track::PID::Kaon:
          case o2::track::PID::Proton:
            keepIt = true;
            break;
        }
        if (!keepIt) {
          continue;
        }
      }
      const char* axisTitle = Form("N_{#sigma}^{TPC}(%s)", pT[id]);
      const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
      histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      axisTitle = Form("#Delta^{TPC}(%s)", pT[id]);
      const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[id])};
      histos.add(hdelta[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});
      const AxisSpec expectedAxis{1000, 0, 10000, Form("Expected(%s)", pT[id])};
      histos.add(hexpected[id].data(), "", kTH2F, {pAxis, expectedAxis});

      histos.addClone(hdelta[id].data(), hdeltapos[id].data());
      histos.addClone(hdelta[id].data(), hdeltaneg[id].data());
      histos.addClone(hnsigma[id].data(), hnsigmapos[id].data());
      histos.addClone(hnsigma[id].data(), hnsigmaneg[id].data());
      histos.addClone(hexpected[id].data(), hexpectedpos[id].data());
      histos.addClone(hexpected[id].data(), hexpectedneg[id].data());

      histos.addClone(hdelta[id].data(), hdeltalf[id].data());
      histos.addClone(hnsigma[id].data(), hnsigmalf[id].data());
      histos.addClone(hexpected[id].data(), hexpectedlf[id].data());

      histos.addClone(hdelta[id].data(), hdeltalfpos[id].data());
      histos.addClone(hdelta[id].data(), hdeltalfneg[id].data());
      histos.addClone(hnsigma[id].data(), hnsigmalfpos[id].data());
      histos.addClone(hnsigma[id].data(), hnsigmalfneg[id].data());
      histos.addClone(hexpected[id].data(), hexpectedlfpos[id].data());
      histos.addClone(hexpected[id].data(), hexpectedlfneg[id].data());
    }
  }

  Filter collisionZVtxFilter = (nabs(o2::aod::collision::posZ) < zVtxMax);
  // Filter collisionEvSelFilter = ((o2::aod::evsel::Sel8 > 0));
  Filter collisionNumContribPV = (minPVcontrib <= o2::aod::collision::numContrib) && (o2::aod::collision::numContrib < maxPVcontrib);
  Filter trackFilter = (requireGlobalTrackInFilter());
  Filter trackFilterCustom = (enableTrackCuts.node() == false) || ((o2::aod::track::tpcInnerParam > minPTPC) && (o2::aod::track::tpcInnerParam < maxPTPC) && (o2::aod::track::tpcSignal > minTPCdEdx) && (o2::aod::track::tpcSignal < maxTPCdEdx));

  template <int id, typename T>
  void fillStd(const T& track)
  {
    histos.fill(HIST(hnsigma[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
    histos.fill(HIST(hdelta[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
    histos.fill(HIST(hexpected[id]), track.tpcInnerParam(), o2::aod::pidutils::tpcExpSignal<id>(track));
    if (track.sign() > 0) {
      histos.fill(HIST(hnsigmapos[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
      histos.fill(HIST(hdeltapos[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
      histos.fill(HIST(hexpectedpos[id]), track.tpcInnerParam(), o2::aod::pidutils::tpcExpSignal<id>(track));
    } else {
      histos.fill(HIST(hnsigmaneg[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
      histos.fill(HIST(hdeltaneg[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
      histos.fill(HIST(hexpectedneg[id]), track.tpcInnerParam(), o2::aod::pidutils::tpcExpSignal<id>(track));
    }
  }

  template <int id, typename T>
  void fillLf(const T& track)
  {
    histos.fill(HIST(hnsigmalf[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
    histos.fill(HIST(hdeltalf[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
    histos.fill(HIST(hexpectedlf[id]), track.tpcInnerParam(), o2::aod::pidutils::tpcExpSignal<id>(track));
    if (track.sign() > 0) {
      histos.fill(HIST(hnsigmalfpos[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
      histos.fill(HIST(hdeltalfpos[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
      histos.fill(HIST(hexpectedlfpos[id]), track.tpcInnerParam(), o2::aod::pidutils::tpcExpSignal<id>(track));
    } else {
      histos.fill(HIST(hnsigmalfneg[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
      histos.fill(HIST(hdeltalfneg[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
      histos.fill(HIST(hexpectedlfneg[id]), track.tpcInnerParam(), o2::aod::pidutils::tpcExpSignal<id>(track));
    }
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  void processAll(soa::Filtered<CollisionCandidate> const&,
                  soa::Filtered<soa::Join<TrackCandidates,
                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                          aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                          aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl>> const& tracks,
                  soa::Filtered<soa::Join<TrackCandidates,
                                          aod::pidTPCLfFullEl, aod::pidTPCLfFullMu, aod::pidTPCLfFullPi,
                                          aod::pidTPCLfFullKa, aod::pidTPCLfFullPr, aod::pidTPCLfFullDe,
                                          aod::pidTPCLfFullTr, aod::pidTPCLfFullHe, aod::pidTPCLfFullAl>> const& lftracks)
  {
    for (const auto& trk : tracks) {
      static_for<0, 8>([&](auto i) {
        fillStd<i>(trk);
      });
    }

    for (const auto& trk : lftracks) {
      static_for<0, 8>([&](auto i) {
        fillLf<i>(trk);
      });
    }
  }
  PROCESS_SWITCH(lfpidqa, processAll, "Process All particles", false);

  void processPiKaPr(soa::Filtered<CollisionCandidate> const&,
                     soa::Filtered<soa::Join<TrackCandidates, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>> const& tracks,
                     soa::Filtered<soa::Join<TrackCandidates, aod::pidTPCLfFullPi, aod::pidTPCLfFullKa, aod::pidTPCLfFullPr>> const& lftracks)
  {
    for (const auto& trk : tracks) {
      fillStd<o2::track::PID::Pion>(trk);
      fillStd<o2::track::PID::Kaon>(trk);
      fillStd<o2::track::PID::Proton>(trk);
    }

    for (const auto& trk : lftracks) {
      fillLf<o2::track::PID::Pion>(trk);
      fillLf<o2::track::PID::Kaon>(trk);
      fillLf<o2::track::PID::Proton>(trk);
    }
  }
  PROCESS_SWITCH(lfpidqa, processPiKaPr, "Process only Pi Ka Pr", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<lfpidqa>(cfgc)}; }
