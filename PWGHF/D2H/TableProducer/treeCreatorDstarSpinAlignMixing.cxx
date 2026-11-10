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

/// \file treeCreatorDstarSpinAlignMixing.cxx
/// \brief Writer of D*+ → D0 ( → π+ K-) π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for Mix-candidate analysis in spin alignment measurement.
///        Serving as a correction for detector acceptance and reconstruction efficiency.
///
/// \author Mingze li <mingze.li@cern.ch>, CCNU/UniTo

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/D2H/Utils/utilsFlow.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/GenVector/RotationZ.h>
#include <Math/Vector3D.h> // IWYU pragma: keep (do not replace with Math/Vector3Dfwd.h)
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <cmath>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::analysis::hf_flow_utils;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace mixing_dstar
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
// D0 related variables
DECLARE_SOA_COLUMN(MD0, mD0, float);
// DECLARE_SOA_COLUMN(PtD0, ptD0, float);
// DECLARE_SOA_COLUMN(PD0, pD0, float);
// DECLARE_SOA_COLUMN(EtaD0, etaD0, float);
// DECLARE_SOA_COLUMN(PhiD0, phiD0, float);
// DECLARE_SOA_COLUMN(YD0, yD0, float);
// soft pion related variables
DECLARE_SOA_COLUMN(PtSoftPi, ptSoftPi, float);
// DECLARE_SOA_COLUMN(PSoftPi, pSoftPi, float);
DECLARE_SOA_COLUMN(EtaSoftPi, etaSoftPi, float);
DECLARE_SOA_COLUMN(PhiSoftPi, phiSoftPi, float);
// DECLARE_SOA_COLUMN(YSoftPi, ySoftPi, float);
// Dstar related variables
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
// DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(CosThetaStarEP, cosThetaStarEP, float);
DECLARE_SOA_COLUMN(MlProbDstarToD0PiBkg, mlProbDstarToD0PiBkg, float);
DECLARE_SOA_COLUMN(MlProbDstarToD0PiPrompt, mlProbDstarToD0PiPrompt, float);
DECLARE_SOA_COLUMN(MlProbDstarToD0PiNonPrompt, mlProbDstarToD0PiNonPrompt, float);
// Events
DECLARE_SOA_COLUMN(ZVtx, zVtx, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(NPvContrib, nPvContrib, uint16_t);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_COLUMN(XQVec, xQVec, float);
DECLARE_SOA_COLUMN(YQVec, yQVec, float);
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t);
// Tracks
DECLARE_SOA_COLUMN(MinAbsEtaTrack, minAbsEtaTrack, float);
DECLARE_SOA_COLUMN(MinNumItsCls, minNumItsCls, uint8_t);
DECLARE_SOA_COLUMN(MinNumTpcCls, minNumTpcCls, uint8_t);
} // namespace mixing_dstar

DECLARE_SOA_TABLE(HfCandDstMix, "AOD", "HFCANDDSTMIX",
                  mixing_dstar::MD0,
                  mixing_dstar::PtSoftPi,
                  mixing_dstar::EtaSoftPi,
                  mixing_dstar::PhiSoftPi,
                  mixing_dstar::M,
                  mixing_dstar::Pt,
                  mixing_dstar::Eta,
                  mixing_dstar::Phi,
                  mixing_dstar::Y,
                  mixing_dstar::MlProbDstarToD0PiBkg,
                  mixing_dstar::MlProbDstarToD0PiPrompt,
                  mixing_dstar::MlProbDstarToD0PiNonPrompt,
                  mixing_dstar::ZVtx,
                  mixing_dstar::Centrality,
                  mixing_dstar::NPvContrib,
                  mixing_dstar::Occupancy,
                  mixing_dstar::XQVec,
                  mixing_dstar::YQVec,
                  mixing_dstar::MinAbsEtaTrack,
                  mixing_dstar::MinNumItsCls,
                  mixing_dstar::MinNumTpcCls,
                  mixing_dstar::GIndexCol,
                  mixing_dstar::TimeStamp);
DECLARE_SOA_TABLE(HfMcParticles, "AOD", "HFMCPARTICLES",
                  mixing_dstar::Pt,
                  mixing_dstar::Eta,
                  mixing_dstar::Phi,
                  mixing_dstar::Y,
                  mixing_dstar::CosThetaStarEP,
                  mixing_dstar::Centrality);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorDstarSpinAlignMixing {
  Produces<o2::aod::HfCandDstMix> rowCandidateMix;
  Produces<o2::aod::HfMcParticles> rowGenMcParticles;

  Configurable<int> qVecDetector{"qVecDetector", 2, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0C: 2)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimator ((None: 0, FT0C: 2, FT0M: 3))"};
  Configurable<int> occEstimator{"occEstimator", 2, "If enabled, replace number of PV contributors with occupancy estimation (0: don't use, 1: ITS, 2: FT0C)"};
  Configurable<bool> activateAxisRotation{"activateAxisRotation", false, "Activate the rotation of the quantisation axis for cos(Theta*) calculation for phi dependence studies"};
  Configurable<float> axisRotationAngle{"axisRotationAngle", 0.f, "Rotation angle of the quantisation axis for cos(Theta*) calculation for phi dependence studies (in rad)"};
  Configurable<bool> activatePartRecoDstar{"activatePartRecoDstar", false, "Activate the study of partly reconstructed D*+ -> D0 (-> KPiPi0) Pi decays"};

  using CollsWithQVecs = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::QvectorBTots, aod::CentFT0Ms, aod::CentFT0Cs>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;
  using CandDstarWSelFlag = soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>;

  using CollisionsWithMcLabels = soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels>>;
  using CollisionsWithMcLabelsAndCent = soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0Cs>>;
  using McParticlesDstarMatched = soa::Join<aod::McParticles, aod::HfCandDstarMcGen>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  using FilteredCandDstarWSelFlagAndMl = soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfMlDstarToD0Pi>>;
  using FilteredCandDstarWSelFlagAndMcAndMl = soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfMlDstarToD0Pi, aod::HfCandDstarMcRec>>;

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == true;

  Preslice<FilteredCandDstarWSelFlagAndMl> dstarWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMcAndMl> dstarWithMcAndMlPerCollision = aod::hf_cand::collisionId;

  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mcparticle::mcCollisionId;

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

  /// \param ptCharmHad is the pt of the candidate
  double sampleDeltaPhi(float ptCharmHad)
  {

    double y0 = 6.49e-2;
    double cons = 8.42e-1;
    double coef = -3.78e-1;
    double x0 = -6.37e-1;
    double norm = -3.64e-1;
    double mu = 4.68e-2;
    double sigma = 1.513;

    double sigmoid = y0 + (cons - y0) / (1 + std::exp(-coef * (ptCharmHad - x0)));
    double peak = norm * std::exp(-(ptCharmHad - mu) * (ptCharmHad - mu) / (2 * sigma * sigma));

    double v2 = sigmoid + peak;
    const double funcMax = 1.0 + 2.0 * std::abs(v2);

    while (true) {
      double deltaPhi = gRandom->Uniform(0, constants::math::TwoPI);
      double funcVal = gRandom->Uniform(0, funcMax);

      double funcDeltaPhi = 1.0 + 2.0 * v2 * std::cos(2.0 * deltaPhi);

      if (funcVal < funcDeltaPhi)
        return deltaPhi;
    }
    return -1.;
  }

  /// \param daughter is the daughter particle
  /// \param ptMin is the minimum pt
  /// \param etaMax is the maximum eta
  /// \return true if daughter is in acceptance
  template <typename Part>
  bool isDaughterInAcceptance(Part const& daughter, float ptMin, float etaMax)
  {
    if (daughter.pt() < ptMin) {
      return false;
    }
    if (std::abs(daughter.eta()) > etaMax) {
      return false;
    }

    return true;
  }

  template <bool DoMc, typename CollType, typename T, typename Trk, typename BcType>
  void fillCandidateTable(CollType const& collision, const T& candidate, Trk const& /*tracks*/, BcType const& /*bcs*/)
  {
    int8_t origin{RecoDecay::OriginType::None};
    bool partRecoDstar{false};
    if constexpr (DoMc) {
      partRecoDstar = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPiPi0 && std::abs(candidate.flagMcMatchRecD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiKPi0;
      bool const signalDstar = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi && std::abs(candidate.flagMcMatchRecD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
      if (!signalDstar && (!partRecoDstar || !activatePartRecoDstar)) { // this candidate is not signal and not partially reconstructed signal, skip
        return;
      }
      origin = candidate.originMcRec();
      int const pdgBhadMother = candidate.pdgBhadMotherPart();
      // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
      if (origin == RecoDecay::OriginType::NonPrompt && (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 557)) { // o2-linter: disable=pdg/explicit-code, magic-number (constants not in the PDG header)
        return;
      }
    }

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
    std::array<float, 3> Qvector = getQvec(collision, qVecDetector.value);
    if constexpr (DoMc) {
      double const deltaPhi = sampleDeltaPhi(candidate.pt());
      double psi = candidate.phi() - deltaPhi;
      Qvector = {static_cast<float>(-std::sin(psi)), static_cast<float>(std::cos(psi)), 0.f};
    }

    rowCandidateMix(
      massD0,
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
      collision.numContrib(),
      getOccupancyColl(collision, occEstimator),
      Qvector[0],
      Qvector[1],
      absEtaTrackMin,
      numItsClsMin,
      numTpcClsMin,
      collision.globalIndex(),
      timeStamp);
  }

  /// \param mcParticle is the Mc particle
  /// \param mcParticles is the table of Mc particles
  /// \param centrality is the centrality of the collision
  template <typename Part, typename Particles>
  void fillGenParticleTable(Part const& mcParticle, Particles const& mcParticles, float const centrality)
  {
    int8_t origin{RecoDecay::OriginType::None};
    bool partRecoDstar{false};
    partRecoDstar = (std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPiPi0) && (std::abs(mcParticle.flagMcMatchGenD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiKPi0);
    bool const signalDstar = (std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) && (std::abs(mcParticle.flagMcMatchGenD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);

    if (!signalDstar && (!activatePartRecoDstar || !partRecoDstar)) { // this particle is not signal and not partially reconstructed signal, skip
      return;
    }
    origin = mcParticle.originMcGen();
    float massDau{0.f}, massCharmHad{0.f};
    if (origin == RecoDecay::OriginType::NonPrompt) {
      auto bHadMother = mcParticles.rawIteratorAt(mcParticle.idxBhadMotherPart() - mcParticles.offset());
      int const pdgBhadMother = std::abs(bHadMother.pdgCode());
      // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
      if (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 557) { // o2-linter: disable=pdg/explicit-code, magic-number (constants not in the PDG header)
        return;
      }
    }
    std::vector<int> listDaughters{};
    std::array<int, 2> const dauPdgs = {kPiPlus, o2::constants::physics::Pdg::kD0};
    RecoDecay::getDaughters(mcParticle, &listDaughters, dauPdgs, 1);
    massDau = MassPiPlus;
    massCharmHad = MassDStar;

    float const rapidity = mcParticle.y();
    if (std::abs(rapidity) > 1.f) { // we do not keep particles with |y| > 1
      return;
    }

    bool areDauInAcc{true};
    float pxDau{-1000.f}, pyDau{-1000.f}, pzDau{-1000.f};
    for (const auto& dauIdx : listDaughters) {
      auto dauPart = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
      if (std::abs(dauPart.pdgCode()) == kPiPlus) {
        pxDau = dauPart.px();
        pyDau = dauPart.py();
        pzDau = dauPart.pz();
        if (areDauInAcc) {
          areDauInAcc = isDaughterInAcceptance(dauPart, 0.1, 0.8);
        }
      } else if (areDauInAcc) { // check also D0 daughters
        std::vector<int> listDaughtersD0{};
        std::array<int, 2> const dauPdgsD0 = {kPiPlus, -kKPlus};
        RecoDecay::getDaughters(mcParticle, &listDaughtersD0, dauPdgsD0, 1);
        for (const auto& dauIdxD0 : listDaughtersD0) {
          auto dauPartD0 = mcParticles.rawIteratorAt(dauIdxD0 - mcParticles.offset());
          areDauInAcc = isDaughterInAcceptance(dauPartD0, 0.3, 0.8);
        }
      }
    }

    float const pxCharmHad = mcParticle.px();
    float const pyCharmHad = mcParticle.py();
    float const pzCharmHad = mcParticle.pz();
    ROOT::Math::PxPyPzMVector const fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
    ROOT::Math::PxPyPzMVector const fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, massCharmHad);
    ROOT::Math::Boost const boost{fourVecMother.BoostToCM()};
    ROOT::Math::PxPyPzMVector const fourVecDauCM = boost(fourVecDau);
    ROOT::Math::XYZVector const threeVecDauCM = fourVecDauCM.Vect();

    // EP analysis for MC
    double const deltaPhi = sampleDeltaPhi(mcParticle.pt());
    double psi = mcParticle.phi() - deltaPhi;
    ROOT::Math::XYZVector qVecNorm = ROOT::Math::XYZVector(-std::sin(psi), std::cos(psi), 0.f);
    if (activateAxisRotation) {
      qVecNorm = ROOT::Math::RotationZ(axisRotationAngle)(qVecNorm);
    }
    float const cosThetaStarEP = qVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

    rowGenMcParticles(
      mcParticle.pt(),
      mcParticle.eta(),
      mcParticle.phi(),
      mcParticle.y(),
      cosThetaStarEP,
      centrality);
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
        fillCandidateTable<false>(collision, dstarCandidate, tracks, bcWithTimeStamps);
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorDstarSpinAlignMixing, processData, "Process data", true);

  void processDstarMcWithMlInPbPb(aod::McCollisions::iterator const&,
                                  McParticlesDstarMatched const& mcParticles,
                                  CollisionsWithMcLabelsAndCent const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                                  FilteredCandDstarWSelFlagAndMcAndMl const& dstarCandidates,
                                  TracksWithExtra const& tracks,
                                  BCsInfo const& bcWithTimeStamps)
  {
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMcAndMlPerCollision, thisCollId);
      for (const auto& dstarCandidate : groupedDstarCandidates) {
        fillCandidateTable<true>(collision, dstarCandidate, tracks, bcWithTimeStamps);
      }
    }
    rowGenMcParticles.reserve(mcParticles.size());

    for (const auto& mcParticle : mcParticles) {
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, mcParticle.mcCollision().globalIndex());
      const auto centrality = getCentralityGenColl(recoCollsPerMcColl, centEstimator);
      fillGenParticleTable(mcParticle, mcParticles, centrality);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorDstarSpinAlignMixing, processDstarMcWithMlInPbPb, "Process Dstar candidates in PbPb MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorDstarSpinAlignMixing>(cfgc)};
}
