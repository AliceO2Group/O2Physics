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

#include <string_view>

#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Framework/ASoAHelpers.h"

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
  static constexpr std::string_view hdeltalf[Np] = {"deltalf/El", "deltalf/Mu", "deltalf/Pi",
                                                    "deltalf/Ka", "deltalf/Pr", "deltalf/De",
                                                    "deltalf/Tr", "deltalf/He", "deltalf/Al"};
  static constexpr std::string_view hnsigmalf[Np] = {"nsigmalf/El", "nsigmalf/Mu", "nsigmalf/Pi",
                                                     "nsigmalf/Ka", "nsigmalf/Pr", "nsigmalf/De",
                                                     "nsigmalf/Tr", "nsigmalf/He", "nsigmalf/Al"};
  Configurable<uint16_t> minPVcontrib{"minPVcontrib", 0, "Minimum number of PV contributors"};
  Configurable<uint16_t> maxPVcontrib{"maxPVcontrib", 10000, "Maximum number of PV contributors"};
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum value for |z_vtx|"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> nBinsDelta{"nBinsDelta", 200, "Number of bins for the Delta"};
  Configurable<float> minDelta{"minDelta", -1000.f, "Minimum Delta in range"};
  Configurable<float> maxDelta{"maxDelta", 1000.f, "Maximum Delta in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 401, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.025f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.025f, "Maximum NSigma in range"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    const AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    for (int id = 0; id < Np; id++) {

      const char* axisTitle = Form("N_{#sigma}^{TPC}(%s)", pT[id]);
      const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
      histos.add(hnsigma[id].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      axisTitle = Form("#Delta^{TPC}(%s)", pT[id]);
      const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[id])};
      histos.add(hdelta[id].data(), axisTitle, kTH2F, {pAxis, deltaAxis});

      histos.addClone(hdelta[id].data(), hdeltalf[id].data());
      histos.addClone(hnsigma[id].data(), hnsigmalf[id].data());
    }
  }

  Filter collisionZVtxFilter = (nabs(o2::aod::collision::posZ) < zVtxMax);
  // Filter collisionEvSelFilter = ((o2::aod::evsel::Sel8 > 0));
  Filter collisionNumContribPV = (minPVcontrib <= o2::aod::collision::numContrib) && (o2::aod::collision::numContrib < maxPVcontrib);

  template <int id, typename T>
  void fillStd(const T& track)
  {
    histos.fill(HIST(hnsigma[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
    histos.fill(HIST(hdelta[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
  }

  template <int id, typename T>
  void fillLf(const T& track)
  {
    histos.fill(HIST(hnsigmalf[id]), track.p(), o2::aod::pidutils::tpcNSigma<id>(track));
    histos.fill(HIST(hdeltalf[id]), track.p(), o2::aod::pidutils::tpcExpSignalDiff<id>(track));
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  void process(soa::Filtered<CollisionCandidate> const& collisions,
               soa::Filtered<soa::Join<TrackCandidates,
                                       aod::pidTPCFullEl,
                                       aod::pidTPCFullMu,
                                       aod::pidTPCFullPi,
                                       aod::pidTPCFullKa,
                                       aod::pidTPCFullPr,
                                       aod::pidTPCFullDe,
                                       aod::pidTPCFullTr,
                                       aod::pidTPCFullHe,
                                       aod::pidTPCFullAl>> const& tracks,
               soa::Filtered<soa::Join<TrackCandidates,
                                       aod::pidTPCLfFullEl,
                                       aod::pidTPCLfFullMu,
                                       aod::pidTPCLfFullPi,
                                       aod::pidTPCLfFullKa,
                                       aod::pidTPCLfFullPr,
                                       aod::pidTPCLfFullDe,
                                       aod::pidTPCLfFullTr,
                                       aod::pidTPCLfFullHe,
                                       aod::pidTPCLfFullAl>> const& lftracks)
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<lfpidqa>(cfgc)}; }
