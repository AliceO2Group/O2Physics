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

/// \file taskMixingDstarCandTreeCreator.cxx
/// \brief Writer of D*+ → D0 ( → π+ K-) π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for Mix-candidate analysis in spin alignment measurement.
///        Serving as a correction for detector acceptance and reconstruction efficiency.
///
/// \author Mingze li <mingze.li@cern.ch>, CCNU/UniTo

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/D2H/Utils/utilsFlow.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::analysis::hf_flow_utils;

namespace o2::aod
{
namespace mixing_dstar
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
// D0 related variables
DECLARE_SOA_COLUMN(MD0, mD0, float);
DECLARE_SOA_COLUMN(PtD0, ptD0, float);
DECLARE_SOA_COLUMN(PD0, pD0, float);
DECLARE_SOA_COLUMN(EtaD0, etaD0, float);
DECLARE_SOA_COLUMN(PhiD0, phiD0, float);
DECLARE_SOA_COLUMN(YD0, yD0, float);
// soft pion related variables
DECLARE_SOA_COLUMN(PtSoftPi, ptSoftPi, float);
DECLARE_SOA_COLUMN(PSoftPi, pSoftPi, float);
DECLARE_SOA_COLUMN(EtaSoftPi, etaSoftPi, float);
DECLARE_SOA_COLUMN(PhiSoftPi, phiSoftPi, float);
DECLARE_SOA_COLUMN(YSoftPi, ySoftPi, float);
// Dstar related variables
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(MlProbDstarToD0PiBkg, mlProbDstarToD0PiBkg, float);
DECLARE_SOA_COLUMN(MlProbDstarToD0PiPrompt, mlProbDstarToD0PiPrompt, float);
DECLARE_SOA_COLUMN(MlProbDstarToD0PiNonPrompt, mlProbDstarToD0PiNonPrompt, float);
// Events
DECLARE_SOA_COLUMN(ZVec, zVec, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, int);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_COLUMN(XQvec, xqVec, float);
DECLARE_SOA_COLUMN(YQvec, yqVec, float);
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t);
// Tracks
DECLARE_SOA_COLUMN(MinAbsEtaTrack, MinabsEtaTrack, float);
DECLARE_SOA_COLUMN(MinNumItsCls, minNumItsCls, int);
DECLARE_SOA_COLUMN(MinNumTpcCls, minNumTpcCls, int); 
} // namespace mixing_dstar

DECLARE_SOA_TABLE(HfCandDstMix, "AOD", "HFCANDDSTMIX",
                  mixing_dstar::MD0,
                  // mixing_dstar::PtD0,
                  // mixing_dstar::EtaD0,
                  // mixing_dstar::PhiD0,
                  // mixing_dstar::YD0,
                  mixing_dstar::PtSoftPi,
                  mixing_dstar::EtaSoftPi,
                  mixing_dstar::PhiSoftPi,
                  // mixing_dstar::YSoftPi,
                  mixing_dstar::M,
                  mixing_dstar::Pt,
                  mixing_dstar::Eta,
                  mixing_dstar::Phi,
                  mixing_dstar::Y,
                  mixing_dstar::MlProbDstarToD0PiBkg,
                  mixing_dstar::MlProbDstarToD0PiPrompt,
                  mixing_dstar::MlProbDstarToD0PiNonPrompt,
                  mixing_dstar::ZVec,
                  mixing_dstar::Centrality,
                  // mixing_dstar::Multiplicity,
                  mixing_dstar::Occupancy,
                  mixing_dstar::XQvec,
                  mixing_dstar::YQvec,
                  mixing_dstar::MinAbsEtaTrack,
                  mixing_dstar::MinNumItsCls,
                  mixing_dstar::MinNumTpcCls,
                  mixing_dstar::GIndexCol,
                  mixing_dstar::TimeStamp);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTaskMixingDstarCandTreeCreator {
  Produces<o2::aod::HfCandDstMix> rowCandidateMix;

  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};
  Configurable<bool> fillMixingCandidateTable{"fillMixingCandidateTable", false, "Switch to fill lite table with candidate properties"};

  Configurable<int> qVecDetector{"qVecDetector", 2, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0C: 2)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimator ((None: 0, FT0C: 2, FT0M: 3))"};
  Configurable<int> occEstimator{"occEstimator", 2, "If enabled, replace number of PV contributors with occupancy estimation (0: don't use, 1: ITS, 2: FT0C)"};


  using CollsWithQVecs = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::QvectorBTots, aod::CentFT0Ms, aod::CentFT0Cs>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;
  using CandDstarWSelFlag = soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
  using FilteredCandDstarWSelFlagAndMl = soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfMlDstarToD0Pi>>;

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;

  Preslice<FilteredCandDstarWSelFlagAndMl> dstarWithMlPerCollision = aod::hf_cand::collisionId;

  void init(InitContext const&)
  {
  }

  /// prongTracks is the vector of daughter tracks
  /// etaMin is the minimum eta
  /// nItsClsMin is the minumum number of clusters in ITS
  /// nTpcClsMin is the minumum number of clusters in TPC
  template <typename Trk>
  void getTrackingInfos(std::vector<Trk> const& prongTracks, float& etaMin, int& nItsClsMin, int& nTpcClsMin)
  {
    etaMin = 10.f;
    nItsClsMin = 10;
    nTpcClsMin = 1000;

    for (const auto& track : prongTracks) {
      if (std::abs(track.eta()) < etaMin) {
        etaMin = std::abs(track.eta());
      }
      if (track.itsNCls() < nItsClsMin) {
        nItsClsMin = track.itsNCls();
      }
      if (track.tpcNClsCrossedRows() < nTpcClsMin) {
        nTpcClsMin = track.tpcNClsCrossedRows();
      }
    }
  }

  template <typename CollType, typename T, typename Trk, typename BcType>
  void fillCandidateTable(CollType const& collision, const T& candidate, Trk const& /*tracks*/, BcType const& /*bcs*/)
  {
    const auto bc = collision.template bc_as<BcType>();
    const int64_t timeStamp = bc.timestamp();

    float massD0{-1.f};
    float massDStar{-1.f};
    float etaSoftPi{-1.f};
    float phiSoftPi{-1.f};
    if (candidate.signSoftPi() > 0) {
      massD0 = candidate.invMassD0();
      massDStar = candidate.invMassDstar();
    } else {
      massD0 = candidate.invMassD0Bar();
      massDStar = candidate.invMassAntiDstar();
    }
    etaSoftPi = RecoDecay::eta(std::array{candidate.pxSoftPi(), candidate.pySoftPi(), candidate.pzSoftPi()});
    phiSoftPi = RecoDecay::phi(std::array{candidate.pxSoftPi(), candidate.pySoftPi(), candidate.pzSoftPi()});

    float absEtaTrackMin{-1.f};
    int numItsClsMin{-1}, numTpcClsMin{-1};

    auto trackProng0 = candidate.template prong0_as<Trk>();
    auto trackProng1 = candidate.template prong1_as<Trk>();
    auto trackProng2 = candidate.template prongPi_as<Trk>();
    getTrackingInfos(std::vector{trackProng0, trackProng1, trackProng2}, absEtaTrackMin, numItsClsMin, numTpcClsMin);
    std::array<float, 3> const Qvector = getQvec(collision, qVecDetector.value);

    if (fillMixingCandidateTable) {
      rowCandidateMix(
        massD0,
        // candidate.ptD0(),
        // candidate.etaD0(),
        // candidate.phiD0(),
        // candidate.yD0(),
        candidate.ptSoftPi(),
        etaSoftPi,
        phiSoftPi,
        massDStar,
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        candidate.y(constants::physics::MassDStar),
        candidate.mlProbDstarToD0Pi()[0],
        candidate.mlProbDstarToD0Pi()[1],
        candidate.mlProbDstarToD0Pi()[2],
        collision.posZ(),
        getCentralityColl(collision, centEstimator),
        getOccupancyColl(collision, occEstimator),
        Qvector[0],
        Qvector[1],
        absEtaTrackMin,
        numItsClsMin,
        numTpcClsMin,
        collision.globalIndex(),
        timeStamp);
    }
  }

  void processData(CollsWithQVecs const& collisions,
                   FilteredCandDstarWSelFlagAndMl const& dstarCandidates,
                   TracksWithExtra const& tracks,
                   aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMlPerCollision, thisCollId);
      for (const auto& dstarCandidate : groupedDstarCandidates) {
        fillCandidateTable(collision, dstarCandidate, tracks, bcWithTimeStamps);
      }  
    }
  }
  PROCESS_SWITCH(HfTaskMixingDstarCandTreeCreator, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMixingDstarCandTreeCreator>(cfgc)};
}
