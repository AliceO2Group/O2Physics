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

/// \file doubleOmegaBuilderHelper.h
/// \brief Building blocks of the double-Omega (Omega + Lambda + kaon) reconstruction, shared by the
///        candidate builder (derived data producer) and the tree creator; and the input features of the
///        double-Omega track-selection BDT, computed from a track as in the training: the columns of the
///        doubleOmegaTrackTreeCreator tree and the derived features of DoubleOmega/track_selection/common.py
///        (add_derived_features), with the same names. Missing values (no TOF, failed DCA propagation) are NaN,
///        as in the training (XGBoost missing values).

#ifndef PWGLF_UTILS_DOUBLEOMEGABUILDERHELPER_H_
#define PWGLF_UTILS_DOUBLEOMEGABUILDERHELPER_H_

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <GPUROOTCartesianFwd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace o2::pwglf::doubleomega
{

struct DoubleOmegaCandidate {
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float x = -999.f;
  float y = -999.f;
  float z = -999.f;
  float pvX = -999.f;
  float pvY = -999.f;
  float pvZ = -999.f;
  float cosPAOmega = -999.f;
  float cosPADirectLambda = -999.f;
  float cosPADoubleOmega = -999.f;
  float cosPAOmegaToPV = -999.f;
  float cosPADirectLambdaToPV = -999.f;
  float dcaXYOmegaToPV = -999.f;
  float dcaXYDirectLambdaToPV = -999.f;
  float dcaXYDirectKaonToPV = -999.f;
  float dcaXYOmegaKaonToPV = -999.f;
  float dcaXYOmegaProtonToPV = -999.f;
  float dcaXYOmegaPionToPV = -999.f;
  float dcaXYDirectProtonToPV = -999.f;
  float dcaXYDirectPionToPV = -999.f;
  float tpcNSigmaDirectKaon = -999.f;
  float tofNSigmaDirectKaon = -999.f;
  uint32_t daughterDetectorMap = 0;
  float mass = -999.f;
  float massOmega = -999.f;
  float massXi = -999.f;
  int8_t sign = 0;
};

struct LambdaCandidate {
  float px = 0.f;
  float py = 0.f;
  float pz = 0.f;
  float mass = -999.f;
};

struct BuiltLambdaCandidate {
  LambdaCandidate candidate;
  o2::track::TrackParCov parentTrack;
  std::array<float, 3> decayVertex{};
  int64_t posTrackId = -1;
  int64_t negTrackId = -1;
  uint32_t daughterDetectorMap = 0; // Proton, pion; three detector bits each.
  float dcaXYProtonToPV = -999.f;
  float dcaXYPionToPV = -999.f;
  float dcaDaughters = -999.f;
};

struct BuiltOmegaCandidate {
  float px = 0.f;
  float py = 0.f;
  float pz = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;
  float massOmega = -999.f;
  float massXi = -999.f;
  int8_t sign = 0;
  int64_t posTrackId = -1;
  int64_t negTrackId = -1;
  int64_t bachelorId = -1;
  o2::track::TrackParCov parentTrack;
  uint32_t daughterDetectorMap = 0; // Kaon, proton, pion; three detector bits each.
  float dcaXYKaonToPV = -999.f;
  float dcaXYProtonToPV = -999.f;
  float dcaXYPionToPV = -999.f;
  float dcaDaughters = -999.f;
};

struct BuiltKaonCandidate {
  o2::track::TrackParCov track;
  int64_t trackId = -1;
  int8_t sign = 0;
  float tpcNSigma = -999.f; // -999 when TPC is unavailable.
  float tofNSigma = -999.f; // -999 when TOF is unavailable.
  uint32_t detectorMap = 0;
};

// Omega, Lambda and kaon candidates of one collision, built once and paired afterwards.
// Lambdas and kaons are stored per charge of the Omega they can be paired with (see chargeIndex).
struct CollisionCandidates {
  std::array<float, 3> primaryVertex{};
  int runNumber = 0;
  float bz = 0.f; // Nominal magnetic field (kG) the collision was recorded with.
  std::vector<BuiltOmegaCandidate> omegas;
  std::array<std::vector<BuiltLambdaCandidate>, 2> lambdas;
  std::array<std::vector<BuiltKaonCandidate>, 2> kaons;
};

struct DoubleOmegaSelections {
  float etaMax = 0.9f;
  float minNITSClus = 0.f;
  float minNTPCClus = 80.f;
  float minNCrossedRows = 100.f;
  float maxNSharedTPCClus = 5.f;
  float minCrossedRowsOverFindable = 0.8f;
  float nSigmaTPCCut = 3.f;
  float dcaBachToPV = 0.05f;
  float dcaKaonToPV = 0.05f;
  float dcaOmegaToPV = 0.f;
  float dcaDirectLambdaToPV = 0.f;
  float dcaV0DauToPV = 0.05f;
  float dcaV0Bach = 1.f;
  float dcaLambdaDaughters = 1.f;
  float mOmegaWindow = 0.01f;
  float mLambdaWindow = 0.01f;
  float maxDoubleOmegaMass = 3.6f;
  float minCosPAOmega = -1.f;
  float minCosPADirectLambda = -1.f;
  float minCosPADoubleOmega = -1.f;
  float minDoubleOmegaDecayRadius = 1.f;
};

class DoubleOmegaBuilderHelper
{
 public:
  DoubleOmegaSelections selections;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3Body;

  void initFitters()
  {
    configureFitter(fitter);
    configureFitter(fitter3Body);
  }

  void setBz(float bz)
  {
    fitter.setBz(bz);
    fitter3Body.setBz(bz);
  }

  void setMatCorrType(o2::base::Propagator::MatCorrType matCorrType)
  {
    matCorr = matCorrType;
    fitter.setMatCorrType(matCorrType);
    fitter3Body.setMatCorrType(matCorrType);
  }

  static int chargeIndex(int8_t charge)
  {
    return charge < 0 ? 0 : 1;
  }

  template <class Track>
  static uint32_t detectorMap(Track const& track)
  {
    return static_cast<uint32_t>(track.hasITS()) |
           (static_cast<uint32_t>(track.hasTPC()) << 1) |
           (static_cast<uint32_t>(track.hasTOF()) << 2);
  }

  template <class Track>
  static BuiltKaonCandidate buildKaon(Track const& track)
  {
    return {getTrackParCov(track),
            track.globalIndex(),
            static_cast<int8_t>(track.sign()),
            track.hasTPC() ? track.tpcNSigmaKa() : -999.f,
            track.hasTOF() ? track.tofNSigmaKa() : -999.f,
            detectorMap(track)};
  }

  // Charge-conjugated kaon with the same trajectory in the opposite magnetic field: a K+ recorded with
  // field -B follows the trajectory of a K- in the field +B. Only q/pt and its correlations change sign.
  static BuiltKaonCandidate mirrorCharge(BuiltKaonCandidate kaon)
  {
    kaon.track.setQ2Pt(-kaon.track.getQ2Pt());
    for (const int element : {o2::track::kSigQ2PtY, o2::track::kSigQ2PtZ, o2::track::kSigQ2PtSnp, o2::track::kSigQ2PtTgl}) {
      kaon.track.setCov(-kaon.track.getCov()[element], element);
    }
    kaon.sign = -kaon.sign;
    return kaon;
  }

  // Kaons of a collision mirrored to the opposite magnetic field (see mirrorCharge), stored per charge
  // of the Omega they can be paired with. Omegas and Lambdas are not mirrored.
  static CollisionCandidates mirroredKaons(CollisionCandidates const& candidates)
  {
    CollisionCandidates mirrored;
    mirrored.primaryVertex = candidates.primaryVertex;
    mirrored.runNumber = candidates.runNumber;
    mirrored.bz = -candidates.bz;
    for (const auto& kaons : candidates.kaons) {
      for (const auto& kaon : kaons) {
        const auto mirroredKaon = mirrorCharge(kaon);
        mirrored.kaons[chargeIndex(mirroredKaon.sign)].push_back(mirroredKaon);
      }
    }
    return mirrored;
  }

  static bool sharesTracks(BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& lambda)
  {
    const std::array<int64_t, 3> omegaTrackIds{omega.posTrackId, omega.negTrackId, omega.bachelorId};
    return std::find(omegaTrackIds.begin(), omegaTrackIds.end(), lambda.posTrackId) != omegaTrackIds.end() ||
           std::find(omegaTrackIds.begin(), omegaTrackIds.end(), lambda.negTrackId) != omegaTrackIds.end();
  }

  static bool sharesTracks(BuiltOmegaCandidate const& omega, BuiltKaonCandidate const& kaon)
  {
    return kaon.trackId == omega.posTrackId || kaon.trackId == omega.negTrackId || kaon.trackId == omega.bachelorId;
  }

  static bool sharesTracks(BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& lambda, BuiltKaonCandidate const& kaon)
  {
    return sharesTracks(omega, lambda) || sharesTracks(omega, kaon) ||
           kaon.trackId == lambda.posTrackId || kaon.trackId == lambda.negTrackId;
  }

  template <class T>
  bool selectTrack(T const& track) const
  {
    if (std::abs(track.eta()) > selections.etaMax) {
      return false;
    }
    if (track.itsNCls() < selections.minNITSClus ||
        track.tpcNClsFound() < selections.minNTPCClus ||
        track.tpcNClsCrossedRows() < selections.minNCrossedRows ||
        track.tpcNClsCrossedRows() < selections.minCrossedRowsOverFindable * track.tpcNClsFindable() ||
        track.tpcNClsShared() > selections.maxNSharedTPCClus) {
      LOG(debug) << "Track failed selection: itsNCls=" << static_cast<int>(track.itsNCls()) << ", tpcNClsFound=" << track.tpcNClsFound() << ", tpcNClsCrossedRows=" << static_cast<int>(track.tpcNClsCrossedRows()) << ", tpcNClsFindable=" << static_cast<int>(track.tpcNClsFindable()) << ", tpcNClsShared=" << static_cast<int>(track.tpcNClsShared());
      return false;
    }
    return true;
  }

  template <class Track>
  bool isSelectedKaon(Track const& track, int8_t charge) const
  {
    return track.sign() == charge && isSelectedKaon(track);
  }

  template <class Track>
  bool isSelectedKaon(Track const& track) const
  {
    return selectTrack(track) &&
           std::abs(track.tpcNSigmaKa()) <= selections.nSigmaTPCCut;
  }

  // Transverse DCA of a track to the collision primary vertex, computed by propagating
  // a copy of the track parameters. Returns false if the propagation fails.
  template <class C, class Track>
  bool dcaXYToPV(C const& collision, Track const& track, float& dcaXY) const
  {
    auto trackParCov = getTrackParCov(track);
    const o2::math_utils::Point3D<float> primaryVertex{collision.posX(), collision.posY(), collision.posZ()};
    std::array<float, 2> dca{};
    if (!o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCov, 2.f, matCorr, &dca)) {
      return false;
    }
    dcaXY = dca[0];
    return true;
  }

  // isLoose: no PID, DCA or mass-window selection (left to the caller)
  template <bool isLoose = false, class C, class Track>
  bool buildLambda(C const& collision, Track const& posTrack, Track const& negTrack, int8_t charge, BuiltLambdaCandidate& builtLambda)
  {
    if (posTrack.sign() <= 0 || negTrack.sign() >= 0 ||
        !selectTrack(posTrack) || !selectTrack(negTrack)) {
      return false;
    }

    const bool isMatter = charge < 0;
    const float protonNSigma = isMatter ? posTrack.tpcNSigmaPr() : negTrack.tpcNSigmaPr();
    const float pionNSigma = isMatter ? negTrack.tpcNSigmaPi() : posTrack.tpcNSigmaPi();
    if ((std::abs(protonNSigma) > selections.nSigmaTPCCut || std::abs(pionNSigma) > selections.nSigmaTPCCut) && !isLoose) {
      return false;
    }
    float dcaXYPos = -999.f;
    float dcaXYNeg = -999.f;
    if (!dcaXYToPV(collision, posTrack, dcaXYPos) || !dcaXYToPV(collision, negTrack, dcaXYNeg)) {
      LOG(debug) << "Propagation to the PV failed for Lambda daughters " << posTrack.globalIndex() << ", " << negTrack.globalIndex();
      return false;
    }
    if ((std::abs(dcaXYPos) < selections.dcaV0DauToPV || std::abs(dcaXYNeg) < selections.dcaV0DauToPV) && !isLoose) {
      LOG(debug) << "Lambda dau don't pass selections: " << "nsigmaPR: " << protonNSigma << ", pionNSigma: " << pionNSigma << ", dcaXYPos: " << dcaXYPos << ", dcaXYNeg: " << dcaXYNeg;
      return false;
    }

    auto posTrackParCov = getTrackParCov(posTrack);
    auto negTrackParCov = getTrackParCov(negTrack);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(posTrackParCov, negTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting Lambda daughter tracks " << posTrack.globalIndex() << ", " << negTrack.globalIndex();
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }

    std::array<float, 3> posMomentum{};
    std::array<float, 3> negMomentum{};
    fitter.getTrack(0).getPxPyPzGlo(posMomentum);
    fitter.getTrack(1).getPxPyPzGlo(negMomentum);
    const std::array<float, 3> lambdaMomentum{
      posMomentum[0] + negMomentum[0],
      posMomentum[1] + negMomentum[1],
      posMomentum[2] + negMomentum[2]};
    const auto& fittedVertex = fitter.getPCACandidate();
    const std::array<float, 3> decayVertex{
      static_cast<float>(fittedVertex[0]), static_cast<float>(fittedVertex[1]), static_cast<float>(fittedVertex[2])};
    const float mass = isMatter ? RecoDecay::m(std::array{posMomentum, negMomentum},
                                               std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged})
                                : RecoDecay::m(std::array{posMomentum, negMomentum},
                                               std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
    const float dcaDaughters = std::sqrt(std::abs(fitter.getChi2AtPCACandidate()));
    const float eta = RecoDecay::eta(lambdaMomentum);
    if (std::abs(eta) > selections.etaMax ||
        ((dcaDaughters > selections.dcaLambdaDaughters ||
          std::abs(mass - o2::constants::physics::MassLambda0) > selections.mLambdaWindow) &&
         !isLoose)) {
      return false;
    }

    builtLambda.candidate = {
      .px = lambdaMomentum[0],
      .py = lambdaMomentum[1],
      .pz = lambdaMomentum[2],
      .mass = mass};
    builtLambda.parentTrack = fitter.createParentTrackParCov(0);
    builtLambda.decayVertex = decayVertex;
    builtLambda.posTrackId = posTrack.globalIndex();
    builtLambda.negTrackId = negTrack.globalIndex();
    // Keep particle-role ordering for matter and antimatter alike.
    builtLambda.daughterDetectorMap = isMatter ? detectorMap(posTrack) | (detectorMap(negTrack) << 3)
                                               : detectorMap(negTrack) | (detectorMap(posTrack) << 3);
    builtLambda.dcaXYProtonToPV = isMatter ? dcaXYPos : dcaXYNeg;
    builtLambda.dcaXYPionToPV = isMatter ? dcaXYNeg : dcaXYPos;
    builtLambda.dcaDaughters = dcaDaughters;
    return true;
  }

  // isLoose: no PID, DCA or mass-window selection (left to the caller)
  template <bool isLoose = false, class C, class Track>
  bool buildOmega(C const& collision, Track const& posTrack, Track const& negTrack, Track const& bachelor, BuiltOmegaCandidate& builtOmega)
  {
    if (bachelor.sign() == 0 ||
        !selectTrack(bachelor) ||
        (std::abs(bachelor.tpcNSigmaKa()) > selections.nSigmaTPCCut && !isLoose)) {
      return false;
    }
    float dcaXYBachelor = -999.f;
    if (!dcaXYToPV(collision, bachelor, dcaXYBachelor) || (!isLoose && std::abs(dcaXYBachelor) < selections.dcaBachToPV)) {
      return false;
    }

    BuiltLambdaCandidate lambda;
    if (!buildLambda<isLoose>(collision, posTrack, negTrack, bachelor.sign(), lambda)) {
      return false;
    }

    auto v0TrackParCov = lambda.parentTrack;
    auto bachelorTrackParCov = getTrackParCov(bachelor);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(v0TrackParCov, bachelorTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting Omega daughter tracks " << posTrack.globalIndex() << ", " << negTrack.globalIndex() << ", " << bachelor.globalIndex();
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }

    std::array<float, 3> lambdaMomentum{};
    std::array<float, 3> bachelorMomentum{};
    fitter.getTrack(0).getPxPyPzGlo(lambdaMomentum);
    fitter.getTrack(1).getPxPyPzGlo(bachelorMomentum);
    const std::array<float, 3> omegaMomentum{
      lambdaMomentum[0] + bachelorMomentum[0],
      lambdaMomentum[1] + bachelorMomentum[1],
      lambdaMomentum[2] + bachelorMomentum[2]};
    const auto& fittedVertex = fitter.getPCACandidate();
    const std::array<float, 3> decayVertex{
      static_cast<float>(fittedVertex[0]), static_cast<float>(fittedVertex[1]), static_cast<float>(fittedVertex[2])};
    const float massOmega = RecoDecay::m(std::array{lambdaMomentum, bachelorMomentum},
                                         std::array{o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged});
    const float massXi = RecoDecay::m(std::array{lambdaMomentum, bachelorMomentum},
                                      std::array{o2::constants::physics::MassLambda0, o2::constants::physics::MassPionCharged});
    const float dcaDaughters = std::sqrt(std::abs(fitter.getChi2AtPCACandidate()));
    const float eta = RecoDecay::eta(omegaMomentum);
    if (std::abs(eta) > selections.etaMax ||
        (!isLoose && dcaDaughters > selections.dcaV0Bach) ||
        (!isLoose && std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > selections.mOmegaWindow)) {
      return false;
    }

    builtOmega = {
      omegaMomentum[0],
      omegaMomentum[1],
      omegaMomentum[2],
      decayVertex[0],
      decayVertex[1],
      decayVertex[2],
      massOmega,
      massXi,
      static_cast<int8_t>(bachelor.sign()),
      posTrack.globalIndex(),
      negTrack.globalIndex(),
      bachelor.globalIndex(),
      fitter.createParentTrackParCov(),
      detectorMap(bachelor) | (lambda.daughterDetectorMap << 3),
      dcaXYBachelor,
      lambda.dcaXYProtonToPV,
      lambda.dcaXYPionToPV,
      dcaDaughters};
    return true;
  }

  // Selections of the Omega, Lambda and kaon candidates that can be applied again on their stored
  // properties (e.g. when reading them back from derived data), to tighten the production selections.
  bool isSelectedOmega(BuiltOmegaCandidate const& omega) const
  {
    return std::abs(RecoDecay::eta(std::array{omega.px, omega.py, omega.pz})) <= selections.etaMax &&
           std::abs(omega.massOmega - o2::constants::physics::MassOmegaMinus) <= selections.mOmegaWindow &&
           omega.dcaDaughters <= selections.dcaV0Bach &&
           std::abs(omega.dcaXYKaonToPV) >= selections.dcaBachToPV &&
           std::abs(omega.dcaXYProtonToPV) >= selections.dcaV0DauToPV &&
           std::abs(omega.dcaXYPionToPV) >= selections.dcaV0DauToPV;
  }

  bool isSelectedLambda(BuiltLambdaCandidate const& lambda) const
  {
    return std::abs(RecoDecay::eta(std::array{lambda.candidate.px, lambda.candidate.py, lambda.candidate.pz})) <= selections.etaMax &&
           std::abs(lambda.candidate.mass - o2::constants::physics::MassLambda0) <= selections.mLambdaWindow &&
           lambda.dcaDaughters <= selections.dcaLambdaDaughters &&
           std::abs(lambda.dcaXYProtonToPV) >= selections.dcaV0DauToPV &&
           std::abs(lambda.dcaXYPionToPV) >= selections.dcaV0DauToPV;
  }

  bool isSelectedKaon(BuiltKaonCandidate const& kaon) const
  {
    return kaon.tpcNSigma > -999.f && std::abs(kaon.tpcNSigma) <= selections.nSigmaTPCCut;
  }

  bool buildDoubleOmega(std::array<float, 3> const& primaryVertex, BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& directLambda, BuiltKaonCandidate const& kaon, DoubleOmegaCandidate& cand)
  {
    auto omegaTrackParCov = omega.parentTrack;
    auto lambdaTrackParCov = directLambda.parentTrack;
    auto kaonTrackParCov = kaon.track;
    int nCandidates = 0;
    try {
      nCandidates = fitter3Body.process(omegaTrackParCov, lambdaTrackParCov, kaonTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting the double-Omega candidate";
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }

    auto omegaTrackAtVertex = fitter3Body.getTrack(0);
    auto lambdaTrackAtVertex = fitter3Body.getTrack(1);
    auto kaonTrackAtVertex = fitter3Body.getTrack(2);
    std::array<float, 3> omegaMomentum{};
    std::array<float, 3> lambdaMomentum{};
    std::array<float, 3> kaonMomentum{};
    omegaTrackAtVertex.getPxPyPzGlo(omegaMomentum);
    lambdaTrackAtVertex.getPxPyPzGlo(lambdaMomentum);
    kaonTrackAtVertex.getPxPyPzGlo(kaonMomentum);

    const std::array<float, 3> totalMomentum{
      omegaMomentum[0] + lambdaMomentum[0] + kaonMomentum[0],
      omegaMomentum[1] + lambdaMomentum[1] + kaonMomentum[1],
      omegaMomentum[2] + lambdaMomentum[2] + kaonMomentum[2]};
    const auto& fittedVertex = fitter3Body.getPCACandidate();
    const std::array<float, 3> decayVertex{
      static_cast<float>(fittedVertex[0]), static_cast<float>(fittedVertex[1]), static_cast<float>(fittedVertex[2])};
    const o2::math_utils::Point3D<float> primaryVertexPoint{primaryVertex[0], primaryVertex[1], primaryVertex[2]};
    const std::array<float, 3> omegaDecayVertex{omega.x, omega.y, omega.z};

    std::array<float, 2> dcaOmega{};
    std::array<float, 2> dcaDirectLambda{};
    std::array<float, 2> dcaDirectKaon{};
    if (!o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertexPoint, omegaTrackAtVertex, 2.f, matCorr, &dcaOmega) ||
        !o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertexPoint, lambdaTrackAtVertex, 2.f, matCorr, &dcaDirectLambda) ||
        !o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertexPoint, kaonTrackAtVertex, 2.f, matCorr, &dcaDirectKaon)) {
      return false;
    }

    const float massSquared = RecoDecay::m2(std::array{omegaMomentum, lambdaMomentum, kaonMomentum},
                                            std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged});

    cand.pt = std::hypot(totalMomentum[0], totalMomentum[1]);
    cand.eta = RecoDecay::eta(totalMomentum);
    cand.phi = std::atan2(totalMomentum[1], totalMomentum[0]);
    cand.x = decayVertex[0];
    cand.y = decayVertex[1];
    cand.z = decayVertex[2];
    cand.pvX = primaryVertex[0];
    cand.pvY = primaryVertex[1];
    cand.pvZ = primaryVertex[2];
    cand.cosPAOmega = RecoDecay::cpa(decayVertex, omegaDecayVertex, omegaMomentum);
    cand.cosPADirectLambda = RecoDecay::cpa(decayVertex, directLambda.decayVertex, lambdaMomentum);
    cand.cosPADoubleOmega = RecoDecay::cpa(primaryVertex, decayVertex, totalMomentum);
    // PV pointing angles use the daughter momenta at their own decay vertices.
    cand.cosPAOmegaToPV = RecoDecay::cpa(primaryVertex, omegaDecayVertex, std::array<float, 3>{omega.px, omega.py, omega.pz});
    cand.cosPADirectLambdaToPV = RecoDecay::cpa(primaryVertex, directLambda.decayVertex, std::array<float, 3>{directLambda.candidate.px, directLambda.candidate.py, directLambda.candidate.pz});
    cand.dcaXYOmegaToPV = dcaOmega[0];
    cand.dcaXYDirectLambdaToPV = dcaDirectLambda[0];
    cand.dcaXYDirectKaonToPV = dcaDirectKaon[0];
    cand.dcaXYOmegaKaonToPV = omega.dcaXYKaonToPV;
    cand.dcaXYOmegaProtonToPV = omega.dcaXYProtonToPV;
    cand.dcaXYOmegaPionToPV = omega.dcaXYPionToPV;
    cand.dcaXYDirectProtonToPV = directLambda.dcaXYProtonToPV;
    cand.dcaXYDirectPionToPV = directLambda.dcaXYPionToPV;
    cand.tpcNSigmaDirectKaon = kaon.tpcNSigma;
    cand.tofNSigmaDirectKaon = kaon.tofNSigma;
    cand.daughterDetectorMap = omega.daughterDetectorMap |
                               (directLambda.daughterDetectorMap << 9) |
                               (kaon.detectorMap << 15);
    cand.mass = std::sqrt(std::max(0.f, massSquared));
    if (cand.cosPAOmega < selections.minCosPAOmega ||
        cand.cosPADirectLambda < selections.minCosPADirectLambda ||
        cand.cosPADoubleOmega < selections.minCosPADoubleOmega ||
        std::hypot(decayVertex[0], decayVertex[1]) < selections.minDoubleOmegaDecayRadius ||
        std::abs(cand.dcaXYOmegaToPV) < selections.dcaOmegaToPV ||
        std::abs(cand.dcaXYDirectLambdaToPV) < selections.dcaDirectLambdaToPV ||
        std::abs(cand.dcaXYDirectKaonToPV) < selections.dcaKaonToPV ||
        cand.mass > selections.maxDoubleOmegaMass) {
      return false;
    }
    cand.massOmega = omega.massOmega;
    cand.massXi = omega.massXi;
    cand.sign = omega.sign;
    return true;
  }

  // Builds the Omega, Lambda and kaon candidates of a collision. If requireOmega is set, Lambdas and
  // kaons are built only if the collision has at least one Omega.
  template <class C, class T, class Cascades, class V0s>
  CollisionCandidates buildCollisionCandidates(C const& collision, T const& tracks, Cascades const& cascades, V0s const& v0s, bool requireOmega)
  {
    CollisionCandidates candidates;
    candidates.primaryVertex = {collision.posX(), collision.posY(), collision.posZ()};
    for (const auto& rawCascade : cascades) {
      const auto omegaV0 = rawCascade.template v0_as<o2::aod::V0s>();
      BuiltOmegaCandidate omega;
      if (buildOmega(collision, omegaV0.template posTrack_as<T>(), omegaV0.template negTrack_as<T>(), rawCascade.template bachelor_as<T>(), omega)) {
        candidates.omegas.push_back(omega);
      }
    }
    if (requireOmega && candidates.omegas.empty()) {
      return candidates;
    }

    for (const int8_t charge : {int8_t{-1}, int8_t{1}}) {
      for (const auto& rawV0 : v0s) {
        BuiltLambdaCandidate lambda;
        if (buildLambda(collision, rawV0.template posTrack_as<T>(), rawV0.template negTrack_as<T>(), charge, lambda)) {
          candidates.lambdas[chargeIndex(charge)].push_back(lambda);
        }
      }
      for (const auto& track : tracks) {
        if (isSelectedKaon(track, charge)) {
          candidates.kaons[chargeIndex(charge)].push_back(buildKaon(track));
        }
      }
    }
    return candidates;
  }

  // Pairs the Omegas, Lambdas and kaons of (possibly different) collisions. The double-Omega candidate
  // is built with respect to the primary vertex of the reference collision.
  template <class F>
  void combine(CollisionCandidates const& reference, CollisionCandidates const& omegaSource,
               CollisionCandidates const& lambdaSource, CollisionCandidates const& kaonSource, F const& onCandidate)
  {
    for (const auto& omega : omegaSource.omegas) {
      for (const auto& lambda : lambdaSource.lambdas[chargeIndex(omega.sign)]) {
        for (const auto& kaon : kaonSource.kaons[chargeIndex(omega.sign)]) {
          // Tracks can also be shared by V0s and cascades associated to different collisions.
          if (sharesTracks(omega, lambda, kaon)) {
            continue;
          }
          DoubleOmegaCandidate cand;
          if (!buildDoubleOmega(reference.primaryVertex, omega, lambda, kaon, cand)) {
            continue;
          }
          onCandidate(cand, omega, lambda, kaon);
        }
      }
    }
  }

  // True if at least one Omega has a Lambda and a kaon of compatible charge that do not share its tracks.
  static bool hasOmegaLambdaKaon(CollisionCandidates const& candidates)
  {
    for (const auto& omega : candidates.omegas) {
      const auto& lambdas = candidates.lambdas[chargeIndex(omega.sign)];
      const auto& kaons = candidates.kaons[chargeIndex(omega.sign)];
      const bool hasLambda = std::any_of(lambdas.begin(), lambdas.end(), [&](BuiltLambdaCandidate const& lambda) { return !sharesTracks(omega, lambda); });
      const bool hasKaon = std::any_of(kaons.begin(), kaons.end(), [&](BuiltKaonCandidate const& kaon) { return !sharesTracks(omega, kaon); });
      if (hasLambda && hasKaon) {
        return true;
      }
    }
    return false;
  }

 private:
  template <class Fitter>
  static void configureFitter(Fitter& dcaFitter)
  {
    dcaFitter.setPropagateToPCA(true);
    dcaFitter.setMaxR(200.);
    dcaFitter.setMaxDZIni(4.);
    dcaFitter.setMinParamChange(1.e-3);
    dcaFitter.setMinRelChi2Change(0.9);
    dcaFitter.setUseAbsDCA(true);
  }
};

// ---------------------------------------------------------------------------------------------
// Input features of the track-selection BDT

enum TrackFeature : uint8_t {
  kFeatPt = 0,
  kFeatEta,
  kFeatPhi,
  kFeatSign,
  kFeatTpcInnerParam,
  kFeatTpcSignal,
  kFeatTpcNSigmaPr,
  kFeatTpcNSigmaPi,
  kFeatTpcNSigmaKa,
  kFeatTofNSigmaPr,
  kFeatTofNSigmaPi,
  kFeatTofNSigmaKa,
  kFeatDcaXY,
  kFeatDcaZ,
  kFeatSigmaDcaXY,
  kFeatSigmaDcaZ,
  kFeatAbsDcaXY,
  kFeatAbsDcaZ,
  kFeatDcaXYSig,
  kFeatDcaZSig,
  kFeatItsClusterMap,
  kFeatItsNCls,
  kFeatItsFirstLayer,
  kFeatItsChi2NCl,
  kFeatTpcNClsFound,
  kFeatTpcNClsCrossedRows,
  kFeatTpcNClsShared,
  kFeatTpcChi2NCl,
  kFeatHasIts,
  kFeatHasTof,
  kFeatIsPVContributor,
  kFeatNumContrib,
  kFeatPosZ,
  kNTrackFeatures
};

// Names of the features, as the columns of the training sample
static constexpr std::array<const char*, kNTrackFeatures> TrackFeatureNames{
  "fPt", "fEta", "fPhi", "fSign", "fTpcInnerParam", "fTpcSignal",
  "fTpcNSigmaPr", "fTpcNSigmaPi", "fTpcNSigmaKa", "fTofNSigmaPr", "fTofNSigmaPi", "fTofNSigmaKa",
  "fDcaXY", "fDcaZ", "fSigmaDcaXY", "fSigmaDcaZ", "fAbsDcaXY", "fAbsDcaZ", "fDcaXYSig", "fDcaZSig",
  "fItsClusterMap", "fItsNCls", "fItsFirstLayer", "fItsChi2NCl",
  "fTpcNClsFound", "fTpcNClsCrossedRows", "fTpcNClsShared", "fTpcChi2NCl",
  "fHasIts", "fHasTof", "fIsPVContributor", "fNumContrib", "fPosZ"};

// Indices of the model inputs (in the order of the model) in the feature array; fatal for unknown names
inline std::vector<int> trackFeatureIndices(std::vector<std::string> const& names)
{
  std::vector<int> indices;
  for (const auto& name : names) {
    int index = -1;
    for (int i = 0; i < kNTrackFeatures; ++i) {
      if (name == TrackFeatureNames[i]) {
        index = i;
        break;
      }
    }
    if (index < 0) {
      LOG(fatal) << "Unknown track-selection feature " << name << " (see PWGLF/Utils/doubleOmegaBuilderHelper.h)";
    }
    indices.push_back(index);
  }
  return indices;
}

// All the features of a track. The DCA to the primary vertex is computed with the uncertainties of the track
// and of the vertex, as in doubleOmegaTrackTreeCreator. Returns false if the propagation to the PV fails (the
// DCA features are then NaN).
template <class C, class T>
bool computeTrackFeatures(C const& collision, T const& track, o2::base::Propagator::MatCorrType matCorr, std::array<float, kNTrackFeatures>& features)
{
  constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
  auto trackParCov = getTrackParCov(track);
  const o2::dataformats::VertexBase primaryVertex{o2::math_utils::Point3D<float>{collision.posX(), collision.posY(), collision.posZ()},
                                                  std::array<float, 6>{collision.covXX(), collision.covXY(), collision.covYY(),
                                                                       collision.covXZ(), collision.covYZ(), collision.covZZ()}};
  o2::dataformats::DCA dca;
  const bool propagated = o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCov, 2.f, matCorr, &dca);
  const float dcaXY = propagated ? dca.getY() : NaN;
  const float dcaZ = propagated ? dca.getZ() : NaN;
  const float sigmaDcaXY = propagated ? std::sqrt(std::max(0.f, static_cast<float>(dca.getSigmaY2()))) : NaN;
  const float sigmaDcaZ = propagated ? std::sqrt(std::max(0.f, static_cast<float>(dca.getSigmaZ2()))) : NaN;

  const uint8_t clusterMap = track.itsClusterMap();
  int firstLayer = 7; // no ITS cluster
  for (int layer = 0; layer < 7; ++layer) {
    if ((clusterMap >> layer) & 1) {
      firstLayer = layer;
      break;
    }
  }

  features[kFeatPt] = track.pt();
  features[kFeatEta] = track.eta();
  features[kFeatPhi] = track.phi();
  features[kFeatSign] = track.sign();
  features[kFeatTpcInnerParam] = track.tpcInnerParam();
  features[kFeatTpcSignal] = track.tpcSignal();
  features[kFeatTpcNSigmaPr] = track.tpcNSigmaPr();
  features[kFeatTpcNSigmaPi] = track.tpcNSigmaPi();
  features[kFeatTpcNSigmaKa] = track.tpcNSigmaKa();
  features[kFeatTofNSigmaPr] = track.hasTOF() ? track.tofNSigmaPr() : NaN;
  features[kFeatTofNSigmaPi] = track.hasTOF() ? track.tofNSigmaPi() : NaN;
  features[kFeatTofNSigmaKa] = track.hasTOF() ? track.tofNSigmaKa() : NaN;
  features[kFeatDcaXY] = dcaXY;
  features[kFeatDcaZ] = dcaZ;
  features[kFeatSigmaDcaXY] = sigmaDcaXY;
  features[kFeatSigmaDcaZ] = sigmaDcaZ;
  features[kFeatAbsDcaXY] = std::abs(dcaXY);
  features[kFeatAbsDcaZ] = std::abs(dcaZ);
  features[kFeatDcaXYSig] = sigmaDcaXY > 0.f ? std::abs(dcaXY) / sigmaDcaXY : NaN;
  features[kFeatDcaZSig] = sigmaDcaZ > 0.f ? std::abs(dcaZ) / sigmaDcaZ : NaN;
  features[kFeatItsClusterMap] = clusterMap;
  features[kFeatItsNCls] = __builtin_popcount(clusterMap);
  features[kFeatItsFirstLayer] = firstLayer;
  features[kFeatItsChi2NCl] = track.itsChi2NCl();
  features[kFeatTpcNClsFound] = track.tpcNClsFound();
  features[kFeatTpcNClsCrossedRows] = track.tpcNClsCrossedRows();
  features[kFeatTpcNClsShared] = track.tpcNClsShared();
  features[kFeatTpcChi2NCl] = track.tpcChi2NCl();
  features[kFeatHasIts] = track.hasITS();
  features[kFeatHasTof] = track.hasTOF();
  features[kFeatIsPVContributor] = track.isPVContributor();
  features[kFeatNumContrib] = collision.numContrib();
  features[kFeatPosZ] = collision.posZ();
  return propagated;
}

} // namespace o2::pwglf::doubleomega

#endif // PWGLF_UTILS_DOUBLEOMEGABUILDERHELPER_H_
