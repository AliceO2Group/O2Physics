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

/// \file fwdtrackUtilities.h
/// \brief Utilities for manipulating parameters of fwdtracks
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>
/// \author Luca Micheletti <luca.micheletti@cern.ch>
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

#ifndef COMMON_CORE_FWDTRACKUTILITIES_H_
#define COMMON_CORE_FWDTRACKUTILITIES_H_

#include "DetectorsBase/GeometryManager.h"
#include "Field/MagneticField.h"
#include "Framework/AnalysisDataModel.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Math/MatrixRepresentationsStatic.h"
#include "Math/SMatrix.h"
#include "TGeoGlobalMagField.h"

#include <type_traits>
#include <vector>

namespace o2::aod
{
namespace fwdtrackutils
{
// Index used to set different options for muon propagation
enum class propagationPoint : int {
  kToVertex = 0,
  kToDCA = 1,
  kToRabs = 2,
  kToMatchingPlane = 3,
};
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

template <typename TFwdTrack, typename TFwdTrackCov>
o2::track::TrackParCovFwd getTrackParCovFwd(TFwdTrack const& track, TFwdTrackCov const& cov)
{
  // This function works for (glMuon, glMuon), (saMuon, saMuon) and (MFTTrack, MFTTrackCov).

  double chi2 = track.chi2();
  if constexpr (std::is_same_v<std::decay_t<TFwdTrackCov>, aod::MFTTracksCov::iterator>) {
    chi2 = track.chi2();
  } else {
    if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      chi2 = track.chi2();
    } else if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      chi2 = track.chi2() * (2.f * track.nClusters() - 5.f);
    }
  }

  SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
  std::vector<double> v1{cov.cXX(), cov.cXY(), cov.cYY(), cov.cPhiX(), cov.cPhiY(),
                         cov.cPhiPhi(), cov.cTglX(), cov.cTglY(), cov.cTglPhi(), cov.cTglTgl(),
                         cov.c1PtX(), cov.c1PtY(), cov.c1PtPhi(), cov.c1PtTgl(), cov.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd trackparCov{track.z(), tpars, tcovs, chi2}; // this is chi2! Not chi2/ndf.
  v1.clear();
  v1.shrink_to_fit();
  return trackparCov;
}

/// propagate fwdtrack to a certain point.
template <typename TFwdTrack, typename TFwdTrackCov, typename TCollision>
o2::dataformats::GlobalFwdTrack propagateMuon(TFwdTrack const& muon, TFwdTrackCov const& cov, TCollision const& collision, const propagationPoint endPoint, const float matchingZ, const float bzkG)
{
  o2::track::TrackParCovFwd trackParCovFwd;
  if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
    trackParCovFwd = getTrackParCovFwd(muon, cov);
  } else if (muon.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
    trackParCovFwd = getTrackParCovFwd(muon, muon);
  } else {
    trackParCovFwd = getTrackParCovFwd(muon, muon);
  }

  o2::dataformats::GlobalFwdTrack propmuon = propagateTrackParCovFwd(trackParCovFwd, muon.trackType(), collision, endPoint, matchingZ, bzkG);
  return propmuon;
}

template <typename TFwdTrackParCov, typename TCollision>
o2::dataformats::GlobalFwdTrack propagateTrackParCovFwd(TFwdTrackParCov const& fwdtrackORG, uint8_t trackType, TCollision const& collision, const propagationPoint endPoint, const float matchingZ, const float bzkG)
{
  // TFwdTrackParCov is o2::track::TrackParCovFwd

  o2::track::TrackParCovFwd fwdtrack(fwdtrackORG);
  o2::dataformats::GlobalFwdTrack propmuon;
  o2::globaltracking::MatchGlobalFwd mMatching;

  if (trackType > 2) { // MCH-MID or MCH standalone
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(fwdtrack.getParameters());
    track.setZ(fwdtrack.getZ());
    track.setCovariances(fwdtrack.getCovariances());
    auto mchTrack = mMatching.FwdtoMCH(track);

    if (endPoint == propagationPoint::kToVertex) {
      o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    } else if (endPoint == propagationPoint::kToDCA) {
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
    } else if (endPoint == propagationPoint::kToRabs) {
      o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
    } else if (endPoint == propagationPoint::kToMatchingPlane) {
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, matchingZ);
    }

    auto proptrack = mMatching.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
  } else if (trackType < 2) { // MFT-MCH-MID
    // const double centerMFT[3] = {0, 0, -61.4};
    // o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    // auto Bz = field->getBz(centerMFT); // Get field at centre of MFT in kG.

    if (endPoint == propagationPoint::kToVertex) {
      // auto geoMan = o2::base::GeometryManager::meanMaterialBudget(fwdtrack.getX(), fwdtrack.getY(), fwdtrack.getZ(), collision.posX(), collision.posY(), collision.posZ());
      // auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      // fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, bzkG, x2x0);
      std::array<double, 3> dcaInfOrig{999.f, 999.f, 999.f};
      fwdtrack.propagateToDCAhelix(bzkG, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);
    } else if (endPoint == propagationPoint::kToDCA) {
      fwdtrack.propagateToZhelix(collision.posZ(), bzkG);
    } else if (endPoint == propagationPoint::kToMatchingPlane) {
      fwdtrack.propagateToZhelix(matchingZ, bzkG);
    }
    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());
  }

  return propmuon;
}

template <typename TFwdTrack, typename TMFTTrack>
o2::dataformats::GlobalFwdTrack refitGlobalMuonCov(TFwdTrack const& muon, TMFTTrack const& mft)
{
  // TFwdTrack and TMFTTrack are o2::track::TrackParCovFwd.

  auto muonCov = muon.getCovariances();
  auto mftCov = mft.getCovariances();

  SMatrix55Std jacob = ROOT::Math::SMatrixIdentity();
  auto tl = muon.getTgl();
  auto invQPt = muon.getInvQPt();
  jacob(4, 3) = tl / (invQPt * std::sqrt(1 + tl * tl));
  jacob(4, 4) = -std::sqrt(1 + tl * tl) / (invQPt * invQPt);

  auto covQP = ROOT::Math::Similarity(jacob, muonCov);
  mftCov(4, 0) = 0;
  mftCov(4, 1) = 0;
  mftCov(4, 2) = 0;
  mftCov(4, 3) = 0;

  mftCov(0, 4) = 0;
  mftCov(1, 4) = 0;
  mftCov(2, 4) = 0;
  mftCov(3, 4) = 0;
  mftCov(4, 4) = covQP(4, 4);

  SMatrix55Std jacobInv = ROOT::Math::SMatrixIdentity();
  auto qp = std::sqrt(1 + tl * tl) / invQPt;
  auto tlMFT = mft.getTgl();
  jacobInv(4, 3) = tlMFT / (qp * std::sqrt(1 + tlMFT * tlMFT));
  jacobInv(4, 4) = -std::sqrt(1 + tlMFT * tlMFT) / (qp * qp);
  auto globalCov = ROOT::Math::Similarity(jacobInv, mftCov);

  auto invQPtGlob = std::sqrt(1 + tlMFT * tlMFT) / qp;

  o2::dataformats::GlobalFwdTrack globalTrack;
  globalTrack.setParameters(mft.getParameters());
  globalTrack.setZ(mft.getZ());
  globalTrack.setInvQPt(invQPtGlob);
  globalTrack.setCovariances(globalCov);

  return globalTrack;
}

} // namespace fwdtrackutils
} // namespace o2::aod

#endif // COMMON_CORE_FWDTRACKUTILITIES_H_
