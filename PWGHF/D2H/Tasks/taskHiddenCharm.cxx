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

/// \file taskHiddenCharm.cxx
/// \brief Analysis task for hidden charm hadrons at midrapidity
///
/// \author A. Palasciano, <antonio.palasciano@cern.ch>, INFN Bari
/// \author S. Politanò <stefano.politano@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;

enum TrackType : uint8_t {
  Pion = 0,
  Kaon,
  Proton
};

struct HfHiddenCharmHadrons {
  Configurable<int> centEstimator{"centEstimator", 0, "Centrality estimation (None: 0, FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4, NTracksPV: 5, FT0CVariant2: 6)"};
  Configurable<float> centralityMin{"centralityMin", 0.f, "Minimum accepted centrality"};
  Configurable<float> centralityMax{"centralityMax", 100.f, "Maximum accepted centrality"};
  Configurable<bool> fillOnlyUnlikeSign{"fillOnlyUnlikeSign", true, "Fill only unlike-sign proton pairs"};
  Configurable<bool> fillOnlyLikeSign{"fillOnlyLikeSign", true, "Fill only like-sign proton pairs"};

  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {1400, 2.8, 4.2}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 10.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisSign{"thnConfigAxisSign", {2, -1., 1.}, ""};

  SliceCache cache;

  using SelectedCollisionsPP = aod::HfRedCollisions;
  using SelectedCollisionsPbPb = soa::Join<aod::Collisions,
                                           aod::HfRedCollisions,
                                           aod::CentFT0As,
                                           aod::CentFT0Cs,
                                           aod::CentFT0Ms,
                                           aod::CentFV0As,
                                           aod::CentNTPVs,
                                           aod::CentFT0CVariant2s>;

  Partition<aod::HcSelTracks> selectedProtons = aod::hf_track_vars_reduced::tracktype == static_cast<uint8_t>(TrackType::Proton);
  Partition<aod::HcSelTracks> selectedPions = aod::hf_track_vars_reduced::tracktype == static_cast<uint8_t>(TrackType::Pion);
  Partition<aod::HcSelTracks> selectedKaons = aod::hf_track_vars_reduced::tracktype == static_cast<uint8_t>(TrackType::Kaon);

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec axisInvMass{thnConfigAxisInvMass, "M_{p#bar{p}} (GeV/#it{c}^{2})"};
    const AxisSpec axisPt{thnConfigAxisPt, "#it{p}_{T}^{p#bar{p}} (GeV/#it{c})"};
    const AxisSpec axisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec axisSign{thnConfigAxisSign, "q_{1} #times q_{2}"};

    registry.add("hSparseHiddenCharm", "Hidden-charm proton-pair candidates", HistType::kTHnSparseF, {axisInvMass, axisPt, axisSign, axisCent});
    registry.add("hPtVsInvMassLikeSign", "Hidden-charm LS M_{inv}", HistType::kTH2D, {axisInvMass, axisPt});
    registry.add("hPtVsInvMassUnlikeSign", "Hidden-charm proton-pair ULS", HistType::kTH2D, {axisInvMass, axisPt});
    registry.add("hPtVsInvMassAllSign", "Hidden-charm proton-pair LS+ULS M_{inv}", HistType::kTH2D, {axisInvMass, axisPt});
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename TCollisions, typename TProtonIds>
  void fillEtac(TCollisions const& collision,
                TProtonIds const& protonIds)
  {
    float cent{-1.f};
    if constexpr (CentEstimator != o2::hf_centrality::CentralityEstimator::None) {
      cent = o2::hf_centrality::getCentralityColl(collision, centEstimator);
      if (cent < centralityMin || cent >= centralityMax) {
        return; // skip events outside the centrality range
      }
    }

    for (const auto& proton1 : protonIds) {
      for (const auto& proton2 : protonIds) {
        if (proton1.trackId() >= proton2.trackId()) {
          continue; // avoid double counting and self-pairs
        }
        const int sign = (proton1.sign() * proton2.sign() > 0) ? 1 : -1;
        if ((sign == 1 && !fillOnlyLikeSign) || (sign == -1 && !fillOnlyUnlikeSign)) {
          continue;
        }
        std::array<float, 3> pVec1{proton1.px(), proton1.py(), proton1.pz()};
        std::array<float, 3> pVec2{proton2.px(), proton2.py(), proton2.pz()};
        float invMass = RecoDecay::m(std::array{pVec1, pVec2}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassProton});
        float ptEtac = RecoDecay::pt(RecoDecay::sumOfVec(pVec1, pVec2));
        registry.fill(HIST("hSparseHiddenCharm"), invMass, ptEtac, sign, cent);
        if (sign == 1) {
          registry.fill(HIST("hPtVsInvMassLikeSign"), invMass, ptEtac);
        } else if (sign == -1) {
          registry.fill(HIST("hPtVsInvMassUnlikeSign"), invMass, ptEtac);
          registry.fill(HIST("hPtVsInvMassAllSign"), invMass, ptEtac);
        }
      }
    }
  }

  void processEtacPP(SelectedCollisionsPP::iterator const& collision,
                     aod::HcSelTracks const& /*tracks*/)
  {
    auto candProtons = selectedProtons->sliceByCached(aod::hf_track_index_reduced::hfRedCollisionId, collision.globalIndex(), cache);
    fillEtac<CentralityEstimator::None>(collision, candProtons);
  }
  PROCESS_SWITCH(HfHiddenCharmHadrons, processEtacPP, "Process Etac candidates for pp", true);

  // void processEtacPbPb(SelectedCollisionsPbPb::iterator const& collisions,
  //                      aod::HcSelTracks const& protonIds)
  //{
  //   fillEtac(collisions, protonIds, true);
  // }
  // PROCESS_SWITCH(HfHiddenCharmHadrons, processEtacPbPb, "Process Etac candidates for PbPb", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfHiddenCharmHadrons>(cfgc)};
}
