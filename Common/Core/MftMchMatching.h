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

/// \file MftMchMatching.h
/// \brief Shared MFT MCH propagation matching and inference workflow.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

#ifndef COMMON_CORE_MFTMCHMATCHING_H_
#define COMMON_CORE_MFTMCHMATCHING_H_

#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGUD/Core/UDHelpers.h"

#include "Common/Core/CollisionAssociation.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/MatchMFTFT0.h"
#include "Common/DataModel/MatchMFTMuonData.h"
#include "Common/DataModel/Multiplicity.h"
#include "Tools/ML/model.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DCAFitter/FwdDCAFitterN.h>
#include <DataFormatsGlobalTracking/RecoContainer.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/TableBuilder.h>
#include <Framework/TimingInfo.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackParam.h>
#include <MFTTracking/Constants.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixFunctions.h>
#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TDatabasePDG.h> // o2-linter: disable=pdg/database (Keep the existing PDG table and exact reconstruction masses.)
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TGeoShape.h>
#include <TGeoVolume.h>
#include <TGrid.h>
#include <THnSparse.h>
#include <TLorentzVector.h> // o2-linter: disable=root/lorentz-vector (Preserve the exact existing spectrum arithmetic and rounding.)
#include <TMatrixDSymEigen.h>
#include <TPDGCode.h>
#include <TUUID.h>
#include <TVector2.h>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// Shared data/MC matching implementation. Keep MC annotations independent of
// reconstructed-track selection, propagation and inference.

// -----------------------------------------------------------------------------
// MftMchMatchingWorkflow
// -----------------------------------------------------------------------------
// Single implementation shared by the data and MC entry points.
// Keep truth labels outside this workflow: they must not affect matching.

// -----------------------------------------------------------------------------
// GlobalMuonMatchingPlane
// -----------------------------------------------------------------------------
// Propagation, conversions, timing and chi2 definitions based on
// PWGDQ/Tasks/global-muon-matcher.cxx (0dfb35775). No realignment is applied here.
// Phi residuals in angular chi2 are wrapped to [-pi, pi) for periodicity.
namespace o2::aod::fwdtrackutils
{
struct GlobalMuonMatchingPlane {
  using SMatrix55Sym = o2::track::SMatrix55Sym;
  using SMatrix55Std = o2::track::SMatrix55Std;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;
  static constexpr int MchMidTrackType = 3;
  static constexpr int ThetaAbsBoundaryDeg = 3;
  static constexpr double SlopeResolutionZ = 535.;
  double mBzAtMftCenter{0.};
  std::map<std::string, std::function<std::tuple<double, int>(const o2::track::TrackParCovFwd&, const o2::track::TrackParCovFwd&)>> mMatchingFunctionMap;
  GlobalMuonMatchingPlane() { initMatchingFunctions(); }
  template <class T>
  o2::mch::TrackParam fwdToMch(const T& fwdtrack)
  {
    // Convert Forward Track parameters and covariances matrix to the
    // MCH track format.

    // Parameter conversion
    const double x2 = fwdtrack.getPhi();
    const double x3 = fwdtrack.getTanl();
    const double x4 = fwdtrack.getInvQPt();

    const auto sinX2 = std::sin(x2);
    const auto cosX2 = std::cos(x2);

    const double alpha1 = cosX2 / x3;
    const double alpha3 = sinX2 / x3;
    const double alpha4 = x4 / std::sqrt(x3 * x3 + sinX2 * sinX2);

    const auto kNorm = std::sqrt(x3 * x3 + sinX2 * sinX2);
    const auto kNorm3 = kNorm * kNorm * kNorm;

    // Covariances matrix conversion
    SMatrix55Std jacobian;
    SMatrix55Sym covariances;

    covariances(0, 0) = fwdtrack.getCovariances()(0, 0);
    covariances(0, 1) = fwdtrack.getCovariances()(0, 1);
    covariances(0, 2) = fwdtrack.getCovariances()(0, 2);
    covariances(0, 3) = fwdtrack.getCovariances()(0, 3);
    covariances(0, 4) = fwdtrack.getCovariances()(0, 4);

    covariances(1, 1) = fwdtrack.getCovariances()(1, 1);
    covariances(1, 2) = fwdtrack.getCovariances()(1, 2);
    covariances(1, 3) = fwdtrack.getCovariances()(1, 3);
    covariances(1, 4) = fwdtrack.getCovariances()(1, 4);

    covariances(2, 2) = fwdtrack.getCovariances()(2, 2);
    covariances(2, 3) = fwdtrack.getCovariances()(2, 3);
    covariances(2, 4) = fwdtrack.getCovariances()(2, 4);

    covariances(3, 3) = fwdtrack.getCovariances()(3, 3);
    covariances(3, 4) = fwdtrack.getCovariances()(3, 4);

    covariances(4, 4) = fwdtrack.getCovariances()(4, 4);

    jacobian(0, 0) = 1;

    jacobian(1, 2) = -sinX2 / x3;
    jacobian(1, 3) = -cosX2 / (x3 * x3);

    jacobian(2, 1) = 1;

    jacobian(3, 2) = cosX2 / x3;
    jacobian(3, 3) = -sinX2 / (x3 * x3);

    jacobian(4, 2) = -x4 * sinX2 * cosX2 / kNorm3;
    jacobian(4, 3) = -x3 * x4 / kNorm3;
    jacobian(4, 4) = 1 / kNorm;
    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    std::array<double, 15> cov = {covariances(0, 0), covariances(1, 0), covariances(1, 1), covariances(2, 0), covariances(2, 1), covariances(2, 2), covariances(3, 0), covariances(3, 1), covariances(3, 2), covariances(3, 3), covariances(4, 0), covariances(4, 1), covariances(4, 2), covariances(4, 3), covariances(4, 4)};
    std::array<double, 5> param = {fwdtrack.getX(), alpha1, fwdtrack.getY(), alpha3, alpha4};

    o2::mch::TrackParam convertedTrack(fwdtrack.getZ(), param.data(), cov.data());
    return {convertedTrack};
  }

  o2::track::TrackParCovFwd mchToFwd(const o2::mch::TrackParam& mchParam)
  {
    // Convert a MCH Track parameters and covariances matrix to the
    // Forward track format. Must be called after propagation though the absorber

    o2::track::TrackParCovFwd convertedTrack;

    // Parameter conversion
    const double alpha1 = mchParam.getNonBendingSlope();
    const double alpha3 = mchParam.getBendingSlope();
    const double alpha4 = mchParam.getInverseBendingMomentum();

    const double x2 = std::atan2(-alpha3, -alpha1);
    const double x3 = -1. / std::sqrt(alpha3 * alpha3 + alpha1 * alpha1);
    const double x4 = alpha4 * -x3 * std::sqrt(1 + alpha3 * alpha3);

    const auto kNorm = alpha1 * alpha1 + alpha3 * alpha3;
    const auto kNorm32 = kNorm * std::sqrt(kNorm);
    const auto slopeLen = std::sqrt(alpha3 * alpha3 + 1);

    // Covariances matrix conversion
    SMatrix55Std jacobian;
    SMatrix55Sym covariances;

    covariances(0, 0) = mchParam.getCovariances()(0, 0);
    covariances(0, 1) = mchParam.getCovariances()(0, 1);
    covariances(0, 2) = mchParam.getCovariances()(0, 2);
    covariances(0, 3) = mchParam.getCovariances()(0, 3);
    covariances(0, 4) = mchParam.getCovariances()(0, 4);

    covariances(1, 1) = mchParam.getCovariances()(1, 1);
    covariances(1, 2) = mchParam.getCovariances()(1, 2);
    covariances(1, 3) = mchParam.getCovariances()(1, 3);
    covariances(1, 4) = mchParam.getCovariances()(1, 4);

    covariances(2, 2) = mchParam.getCovariances()(2, 2);
    covariances(2, 3) = mchParam.getCovariances()(2, 3);
    covariances(2, 4) = mchParam.getCovariances()(2, 4);

    covariances(3, 3) = mchParam.getCovariances()(3, 3);
    covariances(3, 4) = mchParam.getCovariances()(3, 4);

    covariances(4, 4) = mchParam.getCovariances()(4, 4);

    jacobian(0, 0) = 1;

    jacobian(1, 2) = 1;

    jacobian(2, 1) = -alpha3 / kNorm;
    jacobian(2, 3) = alpha1 / kNorm;

    jacobian(3, 1) = alpha1 / kNorm32;
    jacobian(3, 3) = alpha3 / kNorm32;

    jacobian(4, 1) = -alpha1 * alpha4 * slopeLen / kNorm32;
    jacobian(4, 3) = alpha3 * alpha4 * (1 / (std::sqrt(kNorm) * slopeLen) - slopeLen / kNorm32);
    jacobian(4, 4) = slopeLen / std::sqrt(kNorm);

    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    // Set output
    convertedTrack.setX(mchParam.getNonBendingCoor());
    convertedTrack.setY(mchParam.getBendingCoor());
    convertedTrack.setZ(mchParam.getZ());
    convertedTrack.setPhi(x2);
    convertedTrack.setTanl(x3);
    convertedTrack.setInvQPt(x4);
    convertedTrack.setCharge(mchParam.getCharge());
    convertedTrack.setCovariances(covariances);

    return convertedTrack;
  }

  void initMatchingFunctions()
  {
    using SVector2 = ROOT::Math::SVector<double, 2>;
    using SVector4 = ROOT::Math::SVector<double, 4>;
    using SVector5 = ROOT::Math::SVector<double, 5>;

    using SMatrix44 = ROOT::Math::SMatrix<double, 4>;
    using SMatrix45 = ROOT::Math::SMatrix<double, 4, 5>;
    using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
    using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;

    // Define built-in matching functions
    //________________________________________________________________________________
    mMatchingFunctionMap["matchALL"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating all parameters: X,Y, phi, tanl & q/pt

      SMatrix55Sym hK, vK;
      SVector5 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl(), mftTrack.getInvQPt()),
        rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      vK(2, 2) = mftTrack.getCovariances()(2, 2);
      vK(3, 3) = mftTrack.getCovariances()(3, 3);
      vK(4, 4) = mftTrack.getCovariances()(4, 4);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;
      hK(2, 2) = 1.0;
      hK(3, 3) = 1.0;
      hK(4, 4) = 1.0;

      // Covariance of residuals
      SMatrix55Std invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Update Parameters
      rKKminus1 = mK - hK * globalMuonTrackParameters; // Residuals of prediction
      rKKminus1[2] = TVector2::Phi_mpi_pi(rKKminus1[2]);

      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 5};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating positions & angles

      SMatrix45 hK;
      SMatrix44 vK;
      SVector4 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl()),
        rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      vK(2, 2) = mftTrack.getCovariances()(2, 2);
      vK(3, 3) = mftTrack.getCovariances()(3, 3);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;
      hK(2, 2) = 1.0;
      hK(3, 3) = 1.0;

      // Covariance of residuals
      SMatrix44 invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      rKKminus1 = mK - hK * globalMuonTrackParameters;
      rKKminus1[2] = TVector2::Phi_mpi_pi(rKKminus1[2]);

      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 4};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Calculate Matching Chi2 - X and Y positions

      SMatrix25 hK;
      SMatrix22 vK;
      SVector2 mK(mftTrack.getX(), mftTrack.getY()), rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;

      // Covariance of residuals
      SMatrix22 invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      rKKminus1 = mK - hK * globalMuonTrackParameters;
      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return reduced chi2
      return {matchChi2Track, 2};
    };
  }

  o2::track::TrackParCovFwd propagateToZMch(const o2::track::TrackParCovFwd& muon, const double z)
  {
    auto mchTrack = fwdToMch(muon);

    float absFront = -90.f;
    float absBack = -505.f;

    if (muon.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else {
      // all other cases
      o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
    }

    return mchToFwd(mchTrack);
  }

  o2::track::TrackParCovFwd propagateToZMft(const o2::track::TrackParCovFwd& mftTrack, const double z)
  {
    o2::track::TrackParCovFwd trackExtrap{mftTrack};
    trackExtrap.propagateToZ(z, mBzAtMftCenter);
    return trackExtrap;
  }

  template <class EVT, class BC, class TMUON, class TMFT>
  bool isMftMchTimeCompatible(EVT const& collisions,
                              BC const& bcs,
                              TMUON const& mchTrack,
                              TMFT const& mftTrack)
  {
    if (!mchTrack.has_collision() || !mftTrack.has_collision()) {
      return false;
    }

    const auto& collMch = collisions.rawIteratorAt(mchTrack.collisionId());
    const auto& bcMch = bcs.rawIteratorAt(collMch.bcId());
    const auto& collMft = collisions.rawIteratorAt(mftTrack.collisionId());
    const auto& bcMft = bcs.rawIteratorAt(collMft.bcId());

    int64_t deltaBc = static_cast<int64_t>(bcMft.globalBC()) - static_cast<int64_t>(bcMch.globalBC());
    double deltaBcNS = o2::constants::lhc::LHCBunchSpacingNS * deltaBc;
    double deltaTrackTime = mftTrack.trackTime() - mchTrack.trackTime() + deltaBcNS;
    double trackTimeResTot = mftTrack.trackTimeRes() + mchTrack.trackTimeRes();

    return std::fabs(deltaTrackTime) <= trackTimeResTot;
  }

  // Equivalent to VarManager::PropagateMuon(..., kToVertex) with its default
  // zero coordinate shifts. Kept here to avoid a Common -> PWGDQ dependency.
  template <class T, class C>
  o2::dataformats::GlobalFwdTrack propagateToVertexMchNative(const T& track, const C& collision)
  {
    o2::globaltracking::MatchGlobalFwd matching;
    auto mch = matching.FwdtoMCH(getTrackParCovFwd3DShift(track, 0.f, 0.f, 0.f, track));
    o2::mch::TrackExtrap::extrapToVertex(mch, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    return matching.MCHtoFwd(mch);
  }

  template <class T, class C>
  bool pDcaCut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    constexpr double AbsorberEndZ = 505.;
    constexpr double RadToDeg = 180. / o2::constants::math::PI;
    double thetaAbs = std::atan(mchTrack.rAtAbsorberEnd() / AbsorberEndZ) * RadToDeg;

    // propagate muon track to vertex
    auto mchTrackAtVertex = propagateToVertexMchNative(mchTrack, collision);

    // double pUncorr = mchTrack.p();
    double p = mchTrackAtVertex.getP();

    double pDCA = mchTrack.pDca();
    double sigmaPDCA = (thetaAbs < ThetaAbsBoundaryDeg) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = SlopeResolutionZ * slopeRes * p;
    double sigmaPDCAWithRes = std::sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
    return pDCA <= nSigmaPDCA * sigmaPDCAWithRes;
  }
  template <class T, class C>
  bool isGoodMuon(const T& mchTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    // chi2 cut
    if (mchTrack.chi2() > chi2Cut) {
      return false;
    }

    // momentum cut
    if (mchTrack.p() < pCut) {
      return false; // skip low-momentum tracks
    }

    // transverse momentum cut
    if (mchTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mchTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    // RAbs cut
    double rAbs = mchTrack.rAtAbsorberEnd();
    if ((rAbs < rAbsCut[0] || rAbs > rAbsCut[1])) {
      return false;
    }

    // pDCA cut
    return pDcaCut(mchTrack, collision, nSigmaPdcaCut);
  }
  template <class TMCH, class C>
  o2::track::TrackParCovFwd propagateToVertexMch(const TMCH& muon,
                                                 const C& collision)
  {
    auto mchTrack = fwdToMch(fwdtrackutils::getTrackParCovFwd(muon, muon));
    o2::mch::TrackExtrap::extrapToVertex(mchTrack,
                                         collision.posX(),
                                         collision.posY(),
                                         collision.posZ(),
                                         collision.covXX(),
                                         collision.covYY());
    return mchToFwd(mchTrack);
  }
  template <class TMUON, class C, class TagConfig>
  void getTaggedMuons(C const& collisions,
                      TMUON const& muonTracks,
                      std::vector<int64_t>& taggedMuons, TagConfig& configMuonTagging)
  {
    taggedMuons.clear();
    for (const auto& muonTrack : muonTracks) {

      // only consider MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) != MchMidTrackType) {
        continue;
      }

      // only select MCH-MID tracks associated to a collision
      if (!muonTrack.has_collision()) {
        continue;
      }

      const auto& collision = collisions.rawIteratorAt(muonTrack.collisionId());

      // select MCH tracks with strict quality cuts
      if (!isGoodMuon(muonTrack, collision,
                      configMuonTagging.cfgMuonTaggingTrackChi2MchUp,
                      configMuonTagging.cfgMuonTaggingPMchLow,
                      configMuonTagging.cfgMuonTaggingPtMchLow,
                      {configMuonTagging.cfgMuonTaggingEtaMchLow, configMuonTagging.cfgMuonTaggingEtaMchUp},
                      {configMuonTagging.cfgMuonTaggingRabsLow, configMuonTagging.cfgMuonTaggingRabsUp},
                      configMuonTagging.cfgMuonTaggingPdcaUp)) {
        continue;
      }

      // propagate MCH track to the vertex
      auto mchTrackAtVertex = propagateToVertexMch(muonTrack, collision);

      // propagate the track from the vertex to the first MFT plane
      const auto& extrapToMFTfirst = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[0]);
      double rFront = std::sqrt(extrapToMFTfirst.getX() * extrapToMFTfirst.getX() + extrapToMFTfirst.getY() * extrapToMFTfirst.getY());
      if (rFront < configMuonTagging.cfgMuonTaggingRadiusAtMftFrontLow.value || rFront > configMuonTagging.cfgMuonTaggingRadiusAtMftFrontUp.value) {
        continue;
      }

      // propagate the track from the vertex to the last MFT plane
      const auto& extrapToMFTlast = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[9]);
      double rBack = std::sqrt(extrapToMFTlast.getX() * extrapToMFTlast.getX() + extrapToMFTlast.getY() * extrapToMFTlast.getY());
      if (rBack < configMuonTagging.cfgMuonTaggingRadiusAtMftBackLow.value || rBack > configMuonTagging.cfgMuonTaggingRadiusAtMftBackUp.value) {
        continue;
      }

      int64_t muonTrackIndex = muonTrack.globalIndex();
      taggedMuons.emplace_back(muonTrackIndex);
    }
  }
};
} // namespace o2::aod::fwdtrackutils

// -----------------------------------------------------------------------------
// MftMchEventMixing
// -----------------------------------------------------------------------------
namespace o2::analysis::mftmch
{
inline constexpr int TrackParameterCount = 5;
inline constexpr int JointParameterCount = 2 * TrackParameterCount;
inline constexpr int PhiParameter = 2;
inline constexpr int InversePtParameter = 4;
inline constexpr float UpstreamPlaneZ = -77.5f;
inline constexpr float DownstreamPlaneZ = -505.f;
inline constexpr double AbsorberEntranceZ = -90.;
inline constexpr double MaximumCentrality = 100.;
inline constexpr int ExcitationEnergySwitchZ = 13;
inline constexpr double MinimumHelixField = 1.e-12;
inline constexpr int AxisConfigurationSize = 3;
inline constexpr int MinimumBinEdgeCount = 2;
inline constexpr int TensorRank = 2;
inline constexpr int FeatureMaskFactor = 2;
inline constexpr int ReconstructionMethodCount = 3;
inline constexpr int PairChargeCategoryCount = 6;
inline constexpr int MlMethodIndex = 2;
inline constexpr int GeometryWarmupPasses = 2;
inline constexpr int MchChamberCount = 10;
inline constexpr int UpcSelectionMode = 2;
inline constexpr int MatchingPlaneCount = 2;
inline constexpr int MinimumDimuonTracks = 2;
} // namespace o2::analysis::mftmch

namespace o2::analysis::mftmch::mixing
{
struct Config {
  int depth = 5;
  int64_t minDeltaBc = 3564;
  double maxDeltaPhi = 3.14159265358979323846 / 10.; // negative disables position-phi selection.
  double maxDeltaR = 10.;
  double maxDeltaAttempts = .1;
  double maxDeltaZ = 1.;
  bool symmetricCuts = true;              // false preserves signed cuts in global-muon-matcher.
  double maxDeltaMftMultiplicityRel = .1; // -1 disables; symmetric relative to target.
  int maxDeltaMchMultiplicity = 0;        // -1 disables; 0 requires equal counts.
  double maxDeltaCentrality = -1.;        // percentile points; -1 disables.
  double maxDeltaMchPhiPV = -1.;          // wrapped momentum-azimuth difference (rad); -1 disables.
  double maxDeltaMchMomentumRel = -1.;    // |donor pPV - target pPV| / target pPV; negative disables.
  double maxDeltaMchEtaPV = -1.;          // momentum-pseudorapidity difference; -1 disables.
  double maxDeltaCbMultiplicityRel = -1.; // PVMults::multNTracksPV, relative to target.
  bool mftMultiplicityStrict = false;     // true excludes the configured boundary.
  bool mchRadiusStrict = false;
  bool mchEtaStrict = false;
  double maxDeltaMchPtRel = -1.; // |donor pT(PV) - target pT(PV)| / target pT(PV); negative disables.
  bool cbMultiplicityStrict = false;
  bool vertexZStrict = false;
};
struct Event {
  int64_t collision = -1, mch = -1, bc = 0;
  int run = -1, attempts = 0;
  double phi = 0., radius = 0., z = 0.;
  int mftMultiplicity = -1, mchMultiplicity = -1;
  int cbMultiplicity = -1;
  double centrality = std::numeric_limits<double>::quiet_NaN();
  double mchPhiPV = std::numeric_limits<double>::quiet_NaN();
  double mchEtaPV = std::numeric_limits<double>::quiet_NaN();
  double mchMomentumPV = std::numeric_limits<double>::quiet_NaN();
  double mchPtPV = std::numeric_limits<double>::quiet_NaN();
  uint64_t frame = 0; // DPL timeslice; row IDs are scoped to this frame and job.
};
struct Deltas {
  int64_t bc = 0;
  double phi = 0., radius = 0., attempts = 0., z = 0.;
};
inline bool compatible(const Event& target, const Event& donor, const Config& c, Deltas& d)
{
  if (target.collision < 0 || donor.collision < 0 || (target.frame == donor.frame && target.collision == donor.collision) ||
      target.run != donor.run || donor.attempts <= 0) {
    return false;
  }
  if (c.maxDeltaMftMultiplicityRel >= 0.) {
    if (target.mftMultiplicity < 0 || donor.mftMultiplicity < 0) {
      return false;
    }
    if (target.mftMultiplicity == 0) {
      if (donor.mftMultiplicity != 0) {
        return false;
      }
    } else {
      const double relative = std::abs(static_cast<double>(donor.mftMultiplicity) - target.mftMultiplicity) / target.mftMultiplicity;
      if (c.mftMultiplicityStrict ? relative >= c.maxDeltaMftMultiplicityRel : relative > c.maxDeltaMftMultiplicityRel) {
        return false;
      }
    }
  }
  if (c.maxDeltaCbMultiplicityRel >= 0.) {
    if (target.cbMultiplicity < 0 || donor.cbMultiplicity < 0 ||
        (target.cbMultiplicity == 0 ? donor.cbMultiplicity != 0 : (c.cbMultiplicityStrict ? std::abs(static_cast<double>(donor.cbMultiplicity) - target.cbMultiplicity) >= c.maxDeltaCbMultiplicityRel * target.cbMultiplicity : std::abs(static_cast<double>(donor.cbMultiplicity) - target.cbMultiplicity) > c.maxDeltaCbMultiplicityRel * target.cbMultiplicity))) {
      return false;
    }
  }
  if (c.maxDeltaMchMultiplicity >= 0 &&
      (target.mchMultiplicity < 0 || donor.mchMultiplicity < 0 ||
       std::abs(donor.mchMultiplicity - target.mchMultiplicity) > c.maxDeltaMchMultiplicity)) {
    return false;
  }
  if (c.maxDeltaCentrality >= 0. &&
      (!std::isfinite(target.centrality) || !std::isfinite(donor.centrality) ||
       target.centrality < 0. || target.centrality > MaximumCentrality || donor.centrality < 0. || donor.centrality > MaximumCentrality ||
       std::abs(donor.centrality - target.centrality) > c.maxDeltaCentrality)) {
    return false;
  }
  if (c.maxDeltaMchMomentumRel >= 0. &&
      (!std::isfinite(target.mchMomentumPV) || !std::isfinite(donor.mchMomentumPV) ||
       target.mchMomentumPV <= 0. || donor.mchMomentumPV <= 0. ||
       std::abs(donor.mchMomentumPV - target.mchMomentumPV) > c.maxDeltaMchMomentumRel * target.mchMomentumPV)) {
    return false;
  }
  if (c.maxDeltaMchPtRel >= 0. &&
      (!std::isfinite(target.mchPtPV) || !std::isfinite(donor.mchPtPV) ||
       target.mchPtPV <= 0. || donor.mchPtPV <= 0. ||
       std::abs(donor.mchPtPV - target.mchPtPV) > c.maxDeltaMchPtRel * target.mchPtPV)) {
    return false;
  }
  // Select donor MCH directions, never the borrowed MFT or the mixed-pair residual.
  if (c.maxDeltaMchPhiPV >= 0. &&
      (!std::isfinite(target.mchPhiPV) || !std::isfinite(donor.mchPhiPV) ||
       std::abs(std::remainder(donor.mchPhiPV - target.mchPhiPV, 2. * 3.14159265358979323846)) > c.maxDeltaMchPhiPV)) {
    return false;
  }
  if (c.maxDeltaMchEtaPV >= 0. &&
      (!std::isfinite(target.mchEtaPV) || !std::isfinite(donor.mchEtaPV) ||
       (c.mchEtaStrict ? std::abs(donor.mchEtaPV - target.mchEtaPV) >= c.maxDeltaMchEtaPV : std::abs(donor.mchEtaPV - target.mchEtaPV) > c.maxDeltaMchEtaPV))) {
    return false;
  }
  d.bc = target.bc >= donor.bc ? target.bc - donor.bc : donor.bc - target.bc;
  d.phi = std::abs(std::remainder(donor.phi - target.phi, 2. * 3.14159265358979323846));
  d.radius = donor.radius - target.radius;
  d.attempts = target.attempts > 0 ? static_cast<double>(donor.attempts - target.attempts) / target.attempts : 0.;
  d.z = std::abs(donor.z - target.z);
  return d.bc >= c.minDeltaBc &&
         (c.maxDeltaPhi < 0. || (std::isfinite(d.phi) && d.phi <= c.maxDeltaPhi)) &&
         (c.maxDeltaZ < 0. || (std::isfinite(d.z) && (c.vertexZStrict ? d.z < c.maxDeltaZ : d.z <= c.maxDeltaZ))) &&
         (c.maxDeltaR < 0. || (std::isfinite(d.radius) && (c.mchRadiusStrict ? (c.symmetricCuts ? std::abs(d.radius) : d.radius) < c.maxDeltaR : (c.symmetricCuts ? std::abs(d.radius) : d.radius) <= c.maxDeltaR))) &&
         (c.maxDeltaAttempts < 0. || (std::isfinite(d.attempts) && (c.symmetricCuts ? std::abs(d.attempts) : d.attempts) <= c.maxDeltaAttempts));
}
} // namespace o2::analysis::mftmch::mixing

// -----------------------------------------------------------------------------
// MftMchMixingHistory
// -----------------------------------------------------------------------------
namespace o2::analysis::mftmch::mixing
{
// Item must own its state and expose stateCount. Never store AO2D iterators.
template <class Item>
struct History {
  std::deque<Item> events;
  size_t states = 0, evicted = 0, oversized = 0;
  int run = -1;
  double solenoid = NAN, dipole = NAN;
  void clear()
  {
    events.clear();
    states = 0;
  }
  bool setRun(int value)
  {
    if (run == value) {
      return false;
    }
    clear();
    run = value;
    solenoid = dipole = NAN;
    return true;
  }
  bool setField(double sol, double dip)
  {
    if (solenoid == sol && dipole == dip) {
      return false;
    }
    const bool reset = std::isfinite(solenoid) || !events.empty();
    clear();
    solenoid = sol;
    dipole = dip;
    return reset;
  }
  bool push(Item item, size_t maxEvents, size_t maxStates)
  {
    if (item.stateCount > maxStates || maxEvents == 0) {
      ++oversized;
      return false;
    }
    while (!events.empty() && (events.size() >= maxEvents || states + item.stateCount > maxStates)) {
      states -= events.front().stateCount;
      events.pop_front();
      ++evicted;
    }
    states += item.stateCount;
    events.push_back(std::move(item));
    return true;
  }
};
} // namespace o2::analysis::mftmch::mixing

namespace o2::analysis::mftmch
{

using MyCollisions = aod::Collisions;
using MyBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedToFT0>;
using MyMUONs = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMFTs = aod::MFTTracks;

using MyCollision = MyCollisions::iterator;
using MyBC = MyBCs::iterator;
using MyMUON = MyMUONs::iterator;
using MyMFT = MyMFTs::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

[[nodiscard]] inline double particleMass(int pdgCode)
{
  const auto* particle = TDatabasePDG::Instance()->GetParticle(pdgCode); // o2-linter: disable=pdg/database (Preserve existing PDG values without initializing extra particle tables.)
  if (!particle) {
    LOGF(fatal, "Particle mass is missing from the PDG database");
  }
  return particle->Mass();
}
[[nodiscard]] inline double muonMass()
{
  return particleMass(kMuonMinus);
}
inline const float muonMassFloat = static_cast<float>(muonMass());
using MuonFourVector = TLorentzVector; // o2-linter: disable=root/lorentz-vector (Changing this type changes floating-point spectrum values.)
[[nodiscard]] inline MuonFourVector muonMomentum(double x, double y, double z, double mass)
{
  MuonFourVector p;
  p.SetXYZM(x, y, z, mass);
  return p;
}
[[nodiscard]] inline MuonFourVector muonMomentumFromPtEtaPhi(double pt, double eta, double phi, double mass)
{
  MuonFourVector p;
  p.SetPtEtaPhiM(pt, eta, phi, mass);
  return p;
}
[[nodiscard]] inline double muonInvariantMass(MuonFourVector const& p) { return p.M(); }
[[nodiscard]] inline double muonPseudorapidity(MuonFourVector const& p) { return p.Eta(); }
[[nodiscard]] inline double muonRapidity(MuonFourVector const& p) { return p.Rapidity(); }

// Adapt raw FIT records to the found-detector interface used by UDHelpers.
// Missing detectors pass, exactly as in udhelpers::cleanFIT.
struct MatchingFITRecord {
  MatchingFITRecord(aod::FV0As const& fv0Records, aod::FT0s const& ft0Records, aod::FDDs const& fddRecords,
                    int64_t fv0Index, int64_t ft0Index, int64_t fddIndex)
    : fv0s(fv0Records), ft0s(ft0Records), fdds(fddRecords), fv0(fv0Index), ft0(ft0Index), fdd(fddIndex) {}
  std::reference_wrapper<aod::FV0As const> fv0s;
  std::reference_wrapper<aod::FT0s const> ft0s;
  std::reference_wrapper<aod::FDDs const> fdds;
  int64_t fv0 = -1, ft0 = -1, fdd = -1;
  [[nodiscard]] bool has_foundFV0() const { return fv0 >= 0; }
  [[nodiscard]] bool has_foundFT0() const { return ft0 >= 0; }
  [[nodiscard]] bool has_foundFDD() const { return fdd >= 0; }
  [[nodiscard]] auto foundFV0() const { return fv0s.get().rawIteratorAt(fv0); }
  [[nodiscard]] auto foundFT0() const { return ft0s.get().rawIteratorAt(ft0); }
  [[nodiscard]] auto foundFDD() const { return fdds.get().rawIteratorAt(fdd); }
};

// Downstream counterpart of O2 TrackExtrap::extrapToVertexWithoutBranson.
// Material/ELoss/MCS helpers below are adapted from O2 TrackExtrap.cxx,
// revision 024f05624a, retaining its linear absorber-path / fixed-p budget
// approximation, two half-loss corrections, and covariance conventions.
// Changes: upstream -> downstream, E + loss -> E - loss; add MCS at the EXIT.
// This is a study extension, not an official O2 downstream transport API.
struct MftAbsorberExtrap {
  using TrackParam = o2::mch::TrackParam;
  static constexpr double SMuMass = 0.105658; // Same mass as O2 TrackExtrap // o2-linter: disable=pdg/explicit-mass (Preserve the exact O2 TrackExtrap mass and transport results.)
  static constexpr double AbsorberBegin = -90.;
  static constexpr double AbsorberEnd = -505.;
  enum class Status { Success,
                      InvalidInput,
                      FieldPropagation,
                      Geometry,
                      RangeOut,
                      Nonfinite,
                      StepLimit,
                      Count };
  struct NoiseOptions {
    double lossPull = 0.;
    bool includeLossVariance = true;
    bool includeMCS = true;
  };
  struct TransportInfo {
    double pathLength = 0., energyLoss = 0., momentumIn = 0., momentumOut = 0.;
    double energyLossSigma = 0.;
    int materialSteps = 0;
    Status status = Status::Success;
  };

  static bool finite(const TrackParam& track)
  {
    if (!std::isfinite(track.getZ()) || !track.hasCovariances() || !std::isfinite(track.p()) || track.p() <= 0.) {
      return false;
    }
    for (int i = 0; i < TrackParameterCount; ++i) {
      if (!std::isfinite(track.getParameters()(i, 0))) {
        return false;
      }
      for (int j = 0; j < TrackParameterCount; ++j) {
        if (!std::isfinite(track.getCovariances()(i, j))) {
          return false;
        }
      }
    }
    return true;
  }

  // One material interval. With the full -90 -> -505 interval, this retains
  // the previous whole-absorber operation order and representative momentum.
  static bool transportMaterialSlice(TrackParam& track, double zTarget, TransportInfo& result, const NoiseOptions* noise = nullptr)
  {
    TrackParam transported(track);
    // 2. Field-only pilot determines the exit of the straight material chord,
    //    analogous to the unconstrained upstream pilot in official TrackExtrap.
    const std::array<double, 3> xyzIn = {transported.getNonBendingCoor(), transported.getBendingCoor(), transported.getZ()};
    TrackParam pilot(transported);
    if (!o2::mch::TrackExtrap::extrapToZ(pilot, zTarget)) {
      result.status = Status::FieldPropagation;
      return false;
    }
    std::array<double, 3> in = {xyzIn[0], xyzIn[1], xyzIn[2]};
    std::array<double, 3> out = {pilot.getNonBendingCoor(), pilot.getBendingCoor(), zTarget};
    double pathLength = 0., f0 = 0., f1 = 0., f2 = 0., meanRho = 0., loss = 0., lossVariance = 0.;
    const double momentumIn = transported.p();
    if (!getAbsorberCorrectionParam(in, out, momentumIn, pathLength, f0, f1, f2, meanRho, loss, lossVariance) ||
        !std::isfinite(pathLength) || !std::isfinite(loss) || loss < 0. ||
        !std::isfinite(lossVariance) || lossVariance < 0. ||
        !std::isfinite(f0) || !std::isfinite(f1) || !std::isfinite(f2) || f0 < 0.) {
      result.status = Status::Geometry;
      return false;
    }
    // O2's Landau-width approximation adds sigma (not variance) along the
    // material path. Preserve that convention across user intervals: adding
    // independent sigma_step^2 would make straggling vanish as step -> 0.
    const double sigmaStep = std::sqrt(lossVariance);
    if (result.materialSteps > 0) {
      lossVariance = sigmaStep * (2. * result.energyLossSigma + sigmaStep);
    }
    if (noise) {
      loss += noise->lossPull * sigmaStep;
      if (!noise->includeLossVariance) {
        lossVariance = 0.;
      }
    }
    // 3. Range-out guard. Unlike upstream reconstruction, E must DECREASE.
    if (std::hypot(momentumIn, SMuMass) - loss <= SMuMass) {
      result.status = Status::RangeOut;
      return false;
    }
    // 4. Half loss -> field-map/covariance propagation -> MCS -> half loss.
    //    correctELossEffectInAbsorber accepts signed ELoss: negative subtracts E.
    //    The variance is POSITIVE in both propagation directions.
    correctELossEffectInAbsorber(transported, -0.5 * loss, 0.5 * lossVariance);
    if (!finite(transported) || !o2::mch::TrackExtrap::extrapToZCov(transported, zTarget)) {
      result.status = Status::FieldPropagation;
      return false;
    }
    // The O2 MCS helper assumes parameters at zOut. zOut-zIn is negative
    // for a forward muon: position/slope cross terms must therefore be negative.
    if (!noise || noise->includeMCS) {
      addMCSEffectInAbsorber(transported, -pathLength, f0, f1, f2);
    }
    correctELossEffectInAbsorber(transported, -0.5 * loss, 0.5 * lossVariance);
    if (!finite(transported)) {
      result.status = Status::Nonfinite;
      return false;
    }
    result.pathLength += pathLength;
    result.energyLoss += loss;
    result.energyLossSigma += sigmaStep;
    ++result.materialSteps;
    track = transported;
    return true;
  }

  // Shorten a positive-step interval at its first material boundary. A small
  // navigation margin crosses the boundary, avoiding repeated zero steps.
  // The material budget still integrates every medium crossed by the chord.
  static bool limitAtMaterialBoundary(const TrackParam& track, double& zTarget, TransportInfo& result)
  {
    TrackParam pilot(track);
    if (!o2::mch::TrackExtrap::extrapToZ(pilot, zTarget)) {
      result.status = Status::FieldPropagation;
      return false;
    }
    std::array<double, 3> point = {track.getNonBendingCoor(), track.getBendingCoor(), track.getZ()};
    std::array<double, 3> direction = {pilot.getNonBendingCoor() - point[0], pilot.getBendingCoor() - point[1], zTarget - point[2]};
    double length = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
    if (!std::isfinite(length) || length <= 0.) {
      result.status = Status::Nonfinite;
      return false;
    }
    for (double& component : direction) { // o2-linter: disable=const-ref-in-for-loop (Normalizes the navigation direction in place.)
      component /= length;
    }
    if (!gGeoManager->InitTrack(point.data(), direction.data())) {
      result.status = Status::Geometry;
      return false;
    }
    gGeoManager->FindNextBoundary(length);
    double boundary = gGeoManager->GetStep();
    if (!std::isfinite(boundary) || boundary < 0.) {
      result.status = Status::Geometry;
      return false;
    }
    constexpr double NavigationMargin = 1.e-6; // cm, same scale as official material integration
    if (boundary + NavigationMargin < length) {
      double fraction = std::min(1., (boundary + NavigationMargin) / length);
      zTarget = point[2] + direction[2] * length * fraction;
    }
    return true;
  }

  // Input: original MFT converted to TrackParam upstream of absorber.
  // stepCm==0: previous whole-absorber model, unchanged at zEnd=-505.
  // stepCm>0: maximum |delta z| per material update, shortened at boundaries.
  // RK magnetic integration remains adaptive inside every interval.
  // Additional material downstream of -505 is intentionally NOT modelled here.
  // Transactional: leave input unchanged on any failure.
  static bool extrapToAbsorberEndWithoutBranson(TrackParam& track, double zEnd = AbsorberEnd,
                                                TransportInfo* info = nullptr, double stepCm = 0., const NoiseOptions* noise = nullptr)
  {
    TransportInfo result;
    auto fail = [&](Status status) {
      result.status = status;
      if (info) {
        *info = result;
      }
      return false;
    };
    if (!finite(track) || track.getZ() < AbsorberBegin || !std::isfinite(zEnd) || zEnd > AbsorberEnd ||
        !std::isfinite(stepCm) || stepCm < 0.) {
      return fail(Status::InvalidInput);
    }
    if (!gGeoManager) {
      return fail(Status::Geometry);
    }
    TrackParam transported(track);
    if (!o2::mch::TrackExtrap::extrapToZCov(transported, AbsorberBegin) || !finite(transported)) {
      return fail(Status::FieldPropagation);
    }
    result.momentumIn = transported.p();
    constexpr int MaxMaterialSteps = 1000000;
    do {
      if (result.materialSteps >= MaxMaterialSteps) {
        return fail(Status::StepLimit);
      }
      double nextZ = stepCm == 0. ? AbsorberEnd : std::max(AbsorberEnd, transported.getZ() - stepCm);
      if (stepCm > 0. && !limitAtMaterialBoundary(transported, nextZ, result)) {
        return fail(result.status);
      }
      if (!(nextZ < transported.getZ())) {
        return fail(Status::StepLimit);
      }
      if (!transportMaterialSlice(transported, nextZ, result, noise)) {
        return fail(result.status);
      }
    } while (transported.getZ() > AbsorberEnd);
    if (zEnd < AbsorberEnd && !o2::mch::TrackExtrap::extrapToZCov(transported, zEnd)) {
      return fail(Status::FieldPropagation);
    }
    if (!finite(transported)) {
      return fail(Status::Nonfinite);
    }
    result.momentumOut = transported.p();
    track = transported;
    if (info) {
      *info = result;
    }
    return true;
  }

  // Supporting calculations copied from the official private helpers.
  static bool getAbsorberCorrectionParam(std::array<double, 3> const& trackXYZIn, std::array<double, 3> const& trackXYZOut, double pTotal,
                                         double& pathLength, double& f0, double& f1, double& f2,
                                         double& meanRho, double& totalELoss, double& sigmaELoss2)
  {
    /// Parameters used to correct for Multiple Coulomb Scattering and energy loss in absorber
    /// Calculated assuming a linear propagation from trackXYZIn to trackXYZOut (order is important)
    /// pathLength:  path length between trackXYZIn and trackXYZOut (cm)
    /// f0:          0th moment of z calculated with the inverse radiation-length distribution
    /// f1:          1st moment of z calculated with the inverse radiation-length distribution
    /// f2:          2nd moment of z calculated with the inverse radiation-length distribution
    /// meanRho:     average density of crossed material (g/cm3)
    /// totalELoss:  total energy loss in absorber
    /// sigmaELoss2: square of energy loss fluctuation in absorber

    // Check whether the geometry is available
    if (!gGeoManager) {
      LOG(warning) << "geometry is missing";
      return false;
    }

    // Initialize starting point and direction
    pathLength = std::sqrt((trackXYZOut[0] - trackXYZIn[0]) * (trackXYZOut[0] - trackXYZIn[0]) +
                           (trackXYZOut[1] - trackXYZIn[1]) * (trackXYZOut[1] - trackXYZIn[1]) +
                           (trackXYZOut[2] - trackXYZIn[2]) * (trackXYZOut[2] - trackXYZIn[2]));
    if (pathLength < TGeoShape::Tolerance()) {
      LOG(warning) << "path length is too small";
      return false;
    }
    std::array<double, 3> b = {(trackXYZOut[0] - trackXYZIn[0]) / pathLength, (trackXYZOut[1] - trackXYZIn[1]) / pathLength, (trackXYZOut[2] - trackXYZIn[2]) / pathLength};
    TGeoNode* currentnode = gGeoManager->InitTrack(trackXYZIn.data(), b.data());
    if (!currentnode) {
      LOG(warning) << "starting point out of geometry";
      return false;
    }

    // loop over absorber slices and calculate absorber's parameters
    f0 = f1 = f2 = meanRho = totalELoss = 0.;
    double sigmaELoss(0.);
    double zB = trackXYZIn[2];
    double remainingPathLength = pathLength;
    do {

      // Get material properties
      TGeoMaterial* material = currentnode->GetVolume()->GetMedium()->GetMaterial();
      double rho = material->GetDensity(); // material density (g/cm3)
      double x0 = material->GetRadLen();   // radiation-length (cm-1)
      double atomicZ = material->GetZ();   // Z of material
      double atomicZoverA(0.);             // Z/A of material
      if (material->IsMixture()) {
        auto* mixture = dynamic_cast<TGeoMixture*>(material);
        if (!mixture) {
          return false;
        }
        double sum(0.);
        for (int iel = 0; iel < mixture->GetNelements(); ++iel) {
          sum += mixture->GetWmixt()[iel];
          atomicZoverA += mixture->GetWmixt()[iel] * mixture->GetZmixt()[iel] / mixture->GetAmixt()[iel];
        }
        atomicZoverA /= sum;
      } else {
        atomicZoverA = atomicZ / material->GetA();
      }

      // Get path length within this material
      gGeoManager->FindNextBoundary(remainingPathLength);
      double localPathLength = gGeoManager->GetStep() + 1.e-6;
      // Check if boundary within remaining path length. If so, make sure to cross the boundary to prepare the next step
      if (localPathLength >= remainingPathLength) {
        localPathLength = remainingPathLength;
      } else {
        currentnode = gGeoManager->Step();
        if (!currentnode) {
          LOG(warning) << "navigation failed";
          return false;
        }
        if (!gGeoManager->IsEntering()) {
          // make another small step to try to enter in new absorber slice
          gGeoManager->SetStep(0.001);
          currentnode = gGeoManager->Step();
          if (!gGeoManager->IsEntering() || !currentnode) {
            LOG(warning) << "navigation failed";
            return false;
          }
          localPathLength += 0.001;
        }
      }

      // calculate absorber's parameters
      double zE = b[2] * localPathLength + zB;
      double dzB = zB - trackXYZIn[2];
      double dzE = zE - trackXYZIn[2];
      f0 += localPathLength / x0;
      f1 += (dzE * dzE - dzB * dzB) / b[2] / b[2] / x0 / 2.;
      f2 += (dzE * dzE * dzE - dzB * dzB * dzB) / b[2] / b[2] / b[2] / x0 / 3.;
      meanRho += localPathLength * rho;
      // Vacuum has no energy loss; avoid log(0) in the official formula.
      if (rho > 0.) {
        totalELoss += betheBloch(pTotal, localPathLength, rho, atomicZ, atomicZoverA);
        sigmaELoss += energyLossFluctuation(pTotal, localPathLength, rho, atomicZoverA);
      }

      // prepare next step
      zB = zE;
      remainingPathLength -= localPathLength;
    } while (remainingPathLength > TGeoShape::Tolerance());

    meanRho /= pathLength;
    sigmaELoss2 = sigmaELoss * sigmaELoss;

    return true;
  }

  static void addMCSEffectInAbsorber(TrackParam& param, double signedPathLength, double f0, double f1, double f2)
  {
    /// Add to the track parameter covariances the effects of multiple Coulomb scattering
    /// signedPathLength must have the sign of (zOut - zIn) where all other parameters are assumed to be given at zOut.

    // absorber related covariance parameters
    double bendingSlope = param.getBendingSlope();
    double nonBendingSlope = param.getNonBendingSlope();
    double inverseBendingMomentum = param.getInverseBendingMomentum();
    double alpha2 = 0.0136 * 0.0136 * inverseBendingMomentum * inverseBendingMomentum * (1.0 + bendingSlope * bendingSlope) /
                    (1.0 + bendingSlope * bendingSlope + nonBendingSlope * nonBendingSlope); // velocity = 1
    double pathLength = std::abs(signedPathLength);
    double varCoor = alpha2 * (pathLength * pathLength * f0 - 2. * pathLength * f1 + f2);
    double covCorrSlope = std::copysign(1., signedPathLength) * alpha2 * (pathLength * f0 - f1);
    double varSlop = alpha2 * f0;

    // Set MCS covariance matrix
    TMatrixD newParamCov(param.getCovariances());
    // Non bending plane
    newParamCov(0, 0) += varCoor;
    newParamCov(0, 1) += covCorrSlope;
    newParamCov(1, 0) += covCorrSlope;
    newParamCov(1, 1) += varSlop;
    // Bending plane
    newParamCov(2, 2) += varCoor;
    newParamCov(2, 3) += covCorrSlope;
    newParamCov(3, 2) += covCorrSlope;
    newParamCov(3, 3) += varSlop;

    // Set momentum related covariances if B!=0
    if (o2::mch::TrackExtrap::isFieldON()) {
      // compute derivative d(q/Pxy) / dSlopeX and d(q/Pxy) / dSlopeY
      double dqPxydSlopeX = inverseBendingMomentum * nonBendingSlope / (1. + nonBendingSlope * nonBendingSlope + bendingSlope * bendingSlope);
      double dqPxydSlopeY = -inverseBendingMomentum * nonBendingSlope * nonBendingSlope * bendingSlope /
                            (1. + bendingSlope * bendingSlope) / (1. + nonBendingSlope * nonBendingSlope + bendingSlope * bendingSlope);
      // Inverse bending momentum (due to dependences with bending and non bending slopes)
      newParamCov(4, 0) += dqPxydSlopeX * covCorrSlope;
      newParamCov(0, 4) += dqPxydSlopeX * covCorrSlope;
      newParamCov(4, 1) += dqPxydSlopeX * varSlop;
      newParamCov(1, 4) += dqPxydSlopeX * varSlop;
      newParamCov(4, 2) += dqPxydSlopeY * covCorrSlope;
      newParamCov(2, 4) += dqPxydSlopeY * covCorrSlope;
      newParamCov(4, 3) += dqPxydSlopeY * varSlop;
      newParamCov(3, 4) += dqPxydSlopeY * varSlop;
      newParamCov(4, 4) += (dqPxydSlopeX * dqPxydSlopeX + dqPxydSlopeY * dqPxydSlopeY) * varSlop;
    }

    // Set new covariances
    param.setCovariances(newParamCov);
  }

  static double betheBloch(double pTotal, double pathLength, double rho, double atomicZ, double atomicZoverA)
  {
    /// Returns the mean total momentum energy loss of muon with total momentum='pTotal'
    /// in the absorber layer of lenght='pathLength', density='rho', Z='atomicZ' and mean Z/A='atomicZoverA'
    /// This is the parameterization of the Bethe-Bloch formula inspired by Geant.

    static constexpr double BetheBlochCoefficient = 0.307075e-3; // [GeV*cm^2/g]
    static constexpr double ElectronMass = 0.511e-3;             // [GeV/c^2]
    static constexpr double DensityJunctionLow = 0.2 * 2.303;    // density effect first junction point * 2.303
    static constexpr double DensityJunctionHigh = 3. * 2.303;    // density effect second junction point * 2.303

    double bg = pTotal / SMuMass; // beta*gamma
    double bg2 = bg * bg;
    double maxT = 2 * ElectronMass * bg2; // neglecting the electron mass

    // mean exitation energy (GeV)
    double mI = (atomicZ < ExcitationEnergySwitchZ) ? (12. * atomicZ + 7.) * 1.e-9 : (9.76 * atomicZ + 58.8 * std::pow(atomicZ, -0.19)) * 1.e-9;

    // density effect
    double x = std::log(bg);
    double lhwI = std::log(28.816 * 1e-9 * std::sqrt(rho * atomicZoverA) / mI);
    double d2(0.);
    if (x > DensityJunctionHigh) {
      d2 = lhwI + x - 0.5;
    } else if (x > DensityJunctionLow) {
      double r = (DensityJunctionHigh - x) / (DensityJunctionHigh - DensityJunctionLow);
      d2 = lhwI + x - 0.5 + (0.5 - lhwI - DensityJunctionLow) * r * r * r;
    }

    return pathLength * rho * (BetheBlochCoefficient * atomicZoverA * (1 + bg2) / bg2 * (0.5 * std::log(2 * ElectronMass * bg2 * maxT / (mI * mI)) - bg2 / (1 + bg2) - d2));
  }

  static double energyLossFluctuation(double pTotal, double pathLength, double rho, double atomicZoverA)
  {
    /// Returns the total momentum energy loss fluctuation of muon with total momentum='pTotal'
    /// in the absorber layer of lenght='pathLength', density='rho', A='atomicA' and Z='atomicZ'

    static constexpr double BetheBlochCoefficient = 0.307075e-3; // GeV.g^-1.cm^2

    double p2 = pTotal * pTotal;
    double beta2 = p2 / (p2 + SMuMass * SMuMass);

    double fwhm = 2. * BetheBlochCoefficient * rho * pathLength * atomicZoverA / beta2; // FWHM of the energy loss Landau distribution

    return fwhm / std::sqrt(8. * std::log(2.)); // gaussian: fwmh = 2 * srqt(2*ln(2)) * sigma (i.e. fwmh = 2.35 * sigma)
  }

  static void correctELossEffectInAbsorber(TrackParam& param, double eLoss, double sigmaELoss2)
  {
    /// Correct parameters for energy loss and add energy loss fluctuation effect to covariances

    // Get parameter covariances in (X, SlopeX, Y, SlopeY, q*PTot) coordinate system
    TMatrixD newParamCov(param.getCovariances());
    cov2CovP(param.getParameters(), newParamCov);

    // Compute new parameters corrected for energy loss
    double p = param.p();
    double e = std::sqrt(p * p + SMuMass * SMuMass);
    double eCorr = e + eLoss;
    double pCorr = std::sqrt(eCorr * eCorr - SMuMass * SMuMass);
    double nonBendingSlope = param.getNonBendingSlope();
    double bendingSlope = param.getBendingSlope();
    param.setInverseBendingMomentum(param.getCharge() / pCorr *
                                    std::sqrt(1.0 + nonBendingSlope * nonBendingSlope + bendingSlope * bendingSlope) /
                                    std::sqrt(1.0 + bendingSlope * bendingSlope));

    // Add effects of energy loss fluctuation to covariances
    newParamCov(4, 4) += eCorr * eCorr / pCorr / pCorr * sigmaELoss2;

    // Get new parameter covariances in (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
    covP2Cov(param.getParameters(), newParamCov);

    // Set new parameter covariances
    param.setCovariances(newParamCov);
  }

  static void cov2CovP(const TMatrixD& param, TMatrixD& cov)
  {
    /// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, SlopeX, Y, SlopeY, q*PTot)
    /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system

    // charge * total momentum
    double qPTot = std::sqrt(1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0)) /
                   std::sqrt(1. + param(3, 0) * param(3, 0)) / param(4, 0);

    // Jacobian of the opposite transformation
    TMatrixD jacob(5, 5);
    jacob.UnitMatrix();
    jacob(4, 1) = qPTot * param(1, 0) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 3) = -qPTot * param(1, 0) * param(1, 0) * param(3, 0) /
                  (1. + param(3, 0) * param(3, 0)) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 4) = -qPTot / param(4, 0);

    // compute covariances in new coordinate system
    TMatrixD tmp(cov, TMatrixD::kMultTranspose, jacob);
    cov.Mult(jacob, tmp);
  }

  static void covP2Cov(const TMatrixD& param, TMatrixD& covP)
  {
    /// change coordinate system: (X, SlopeX, Y, SlopeY, q*PTot) -> (X, SlopeX, Y, SlopeY, q/Pyz)
    /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system

    // charge * total momentum
    double qPTot = std::sqrt(1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0)) /
                   std::sqrt(1. + param(3, 0) * param(3, 0)) / param(4, 0);

    // Jacobian of the transformation
    TMatrixD jacob(5, 5);
    jacob.UnitMatrix();
    jacob(4, 1) = param(4, 0) * param(1, 0) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 3) = -param(4, 0) * param(1, 0) * param(1, 0) * param(3, 0) /
                  (1. + param(3, 0) * param(3, 0)) / (1. + param(1, 0) * param(1, 0) + param(3, 0) * param(3, 0));
    jacob(4, 4) = -param(4, 0) / qPTot;

    // compute covariances in new coordinate system
    TMatrixD tmp(covP, TMatrixD::kMultTranspose, jacob);
    covP.Mult(jacob, tmp);
  }
};

// Conditional MFT prediction using this candidate MCH's momentum. Input detector
// errors are independent; endpoint errors are NOT. Shared absorber energy loss
// uses one common unit-normal pull; MCS is introduced once on the forward path.
struct MftMchJointTransport {
  using Fwd = o2::track::TrackParCovFwd;
  using Conversion = o2::aod::fwdtrackutils::GlobalMuonMatchingPlane;
  struct State {
    Fwd seed, mft, mch;
    double sigmaLossUp = 0.;
  };
  struct Result {
    State state;
    TMatrixD joint{10, 10}, residual{5, 5}, seedCov{5, 5}, seedMchCross{5, 5}, jacobian{10, 10};
    bool hasNominalState = false;
    double chiAll = std::numeric_limits<double>::quiet_NaN(), chiFour = std::numeric_limits<double>::quiet_NaN(), chiXY = std::numeric_limits<double>::quiet_NaN();
    int rankAll = 0, rankFour = 0, rankXY = 0, status = 0;
  };

  static Fwd perturb(const Fwd& track, int index, double amount)
  {
    Fwd out(track);
    auto parameters = out.getParameters();
    parameters[index] += amount;
    out.setParameters(parameters);
    return out;
  }
  static double difference(double a, double b, int index)
  {
    return index % TrackParameterCount == PhiParameter ? TVector2::Phi_mpi_pi(a - b) : a - b;
  }
  static bool evaluate(Conversion& conversion, const Fwd& mft, const Fwd& mch,
                       double zEnd, double step, bool mchCharge, double pull, State& result,
                       double downstreamPull = std::numeric_limits<double>::quiet_NaN())
  {
    if (mft.getZ() < AbsorberEntranceZ || mft.getTanl() >= 0. || !std::isfinite(mch.getInvQPt()) || mch.getInvQPt() == 0.) {
      return false;
    }
    auto up = conversion.fwdToMch(mch);
    if (!o2::mch::TrackExtrap::extrapToVertexWithoutBranson(up, mft.getZ())) {
      return false;
    }
    // Reproduce the official unconstrained upstream material chord for sigmaE.
    auto exit = conversion.fwdToMch(mch);
    if (!o2::mch::TrackExtrap::extrapToZ(exit, -505.)) {
      return false;
    }
    auto entrance = exit;
    if (!o2::mch::TrackExtrap::extrapToZ(entrance, -90.)) {
      return false;
    }
    std::array<double, 3> in = {entrance.getNonBendingCoor(), entrance.getBendingCoor(), -90.};
    std::array<double, 3> out = {exit.getNonBendingCoor(), exit.getBendingCoor(), -505.};
    double length = 0., f0 = 0., f1 = 0., f2 = 0., rho = 0., loss = 0., variance = 0.;
    if (!MftAbsorberExtrap::getAbsorberCorrectionParam(in, out, exit.p(), length, f0, f1, f2, rho, loss, variance) ||
        !std::isfinite(variance) || variance < 0.) {
      return false;
    }
    result.sigmaLossUp = std::sqrt(variance);
    double energy = std::hypot(up.p(), MftAbsorberExtrap::SMuMass) + pull * result.sigmaLossUp;
    if (!std::isfinite(energy) || energy <= MftAbsorberExtrap::SMuMass) {
      return false;
    }
    double p = std::sqrt((energy - MftAbsorberExtrap::SMuMass) * (energy + MftAbsorberExtrap::SMuMass));
    result.seed = mft;
    double charge = std::copysign(1., mchCharge ? mch.getInvQPt() : mft.getInvQPt());
    result.seed.setInvQPt(charge * std::hypot(1., mft.getTanl()) / p);
    // This covariance is a PROCESS-only forward covariance. Measurement errors
    // are propagated below through the joint Jacobian, never added twice here.
    o2::track::SMatrix55Sym zero;
    for (int i = 0; i < TrackParameterCount; ++i) {
      for (int j = 0; j <= i; ++j) {
        zero(i, j) = 0.;
      }
    }
    result.seed.setCovariances(zero);
    auto forward = conversion.fwdToMch(result.seed);
    MftAbsorberExtrap::NoiseOptions noise;
    noise.lossPull = std::isfinite(downstreamPull) ? downstreamPull : pull;
    noise.includeLossVariance = false;
    noise.includeMCS = true;
    if (!MftAbsorberExtrap::extrapToAbsorberEndWithoutBranson(forward, zEnd, nullptr, step, &noise)) {
      return false;
    }
    result.mft = conversion.mchToFwd(forward);
    auto comparison = conversion.fwdToMch(mch);
    if (!o2::mch::TrackExtrap::extrapToZCov(comparison, zEnd)) {
      return false;
    }
    result.mch = conversion.mchToFwd(comparison);
    for (int i = 0; i < TrackParameterCount; ++i) {
      if (!std::isfinite(result.mft.getParameters()[i]) || !std::isfinite(result.mch.getParameters()[i])) {
        return false;
      }
    }
    return true;
  }

  // Use the official matching prescription: directly invert the residual
  // covariance, without a PSD gate or eigenvalue clipping. Preserve the
  // pair correlations already included in the downstream residual covariance.
  template <unsigned int Dimension>
  static bool chiDirect(const TMatrixD& covariance, const std::vector<double>& delta,
                        double& value, int& rank)
  {
    ROOT::Math::SMatrix<double, Dimension> inverse;
    ROOT::Math::SVector<double, Dimension> residual;
    for (unsigned int i = 0; i < Dimension; ++i) {
      residual[i] = delta[i];
      if (!std::isfinite(residual[i])) {
        return false;
      }
      for (unsigned int j = 0; j < Dimension; ++j) {
        inverse(i, j) = covariance(i, j);
        if (!std::isfinite(inverse(i, j))) {
          return false;
        }
      }
    }
    inverse.Invert(); // As in official matchALL: no PSD/inversion-status gate.
    value = ROOT::Math::Similarity(residual, inverse);
    rank = Dimension;
    return std::isfinite(value);
  }

  static TMatrixD inputCovariance(const Fwd& mft, const Fwd& mch)
  {
    TMatrixD input(10, 10);
    input.Zero();
    for (int i = 0; i < TrackParameterCount; ++i) {
      for (int j = 0; j < TrackParameterCount; ++j) {
        // MFT momentum is replaced, so its original uncertainty/correlations
        // are deliberately absent. Its measured angular/position marginal stays.
        if (i < InversePtParameter && j < InversePtParameter) {
          input(i, j) = mft.getCovariances()(i, j);
        }
        input(i + 5, j + 5) = mch.getCovariances()(i, j);
      }
    }
    return input;
  }
  static bool compute(Conversion& conversion, const Fwd& mft, const Fwd& mch,
                      double zEnd, double step, bool mchCharge, double derivativeScale, Result& result)
  {
    result.status = 1;
    if (!evaluate(conversion, mft, mch, zEnd, step, mchCharge, 0., result.state)) {
      return false;
    }
    result.hasNominalState = true;
    TMatrixD seedJac(5, 10);
    seedJac.Zero();
    result.jacobian.Zero();
    auto input = inputCovariance(mft, mch);
    const std::array<double, 5> floors = {1.e-3, 1.e-3, 1.e-5, 1.e-4, 1.e-6};
    for (int column = 0; column < JointParameterCount; ++column) {
      if (column == InversePtParameter) {
        continue; // Never use MFT's own momentum magnitude.
      }
      int index = column % 5;
      const auto& original = column < TrackParameterCount ? mft : mch;
      double h = derivativeScale * std::max({floors[index], 1.e-4 * std::abs(original.getParameters()[index]),
                                             .01 * std::sqrt(std::max(0., input(column, column)))});
      if (index == InversePtParameter) {
        h = std::min(h, .1 * std::abs(original.getInvQPt()));
      }
      State plus, minus;
      auto p = perturb(original, index, h), m = perturb(original, index, -h);
      bool ok = column < TrackParameterCount ? evaluate(conversion, p, mch, zEnd, step, mchCharge, 0., plus) && evaluate(conversion, m, mch, zEnd, step, mchCharge, 0., minus)
                                             : evaluate(conversion, mft, p, zEnd, step, mchCharge, 0., plus) && evaluate(conversion, mft, m, zEnd, step, mchCharge, 0., minus);
      if (!ok) {
        result.status = 2;
        return false;
      }
      for (int row = 0; row < TrackParameterCount; ++row) {
        result.jacobian(row, column) = difference(plus.mft.getParameters()[row], minus.mft.getParameters()[row], row) / (2 * h);
        result.jacobian(row + 5, column) = difference(plus.mch.getParameters()[row], minus.mch.getParameters()[row], row) / (2 * h);
        seedJac(row, column) = difference(plus.seed.getParameters()[row], minus.seed.getParameters()[row], row) / (2 * h);
      }
    }
    TMatrixD temporary(result.jacobian, TMatrixD::kMult, input);
    result.joint = TMatrixD(temporary, TMatrixD::kMultTranspose, result.jacobian);
    TMatrixD seedTmp(seedJac, TMatrixD::kMult, input);
    result.seedCov = TMatrixD(seedTmp, TMatrixD::kMultTranspose, seedJac);
    for (int i = 0; i < TrackParameterCount; ++i) {
      for (int j = 0; j < TrackParameterCount; ++j) {
        result.seedMchCross(i, j) = seedTmp(i, j + 5);
      }
    }
    // One absorber-loss nuisance, shared with opposite correction directions.
    // Perfect upstream/downstream correlation is the explicit same-particle
    // model; widths can differ because the two estimated paths differ.
    State plus, minus;
    constexpr double PullStep = 1.e-3;
    if (!evaluate(conversion, mft, mch, zEnd, step, mchCharge, PullStep, plus) ||
        !evaluate(conversion, mft, mch, zEnd, step, mchCharge, -PullStep, minus)) {
      result.status = 3;
      return false;
    }
    std::array<double, 5> lossDerivative{}, seedDerivative{};
    for (int i = 0; i < TrackParameterCount; ++i) {
      lossDerivative[i] = difference(plus.mft.getParameters()[i], minus.mft.getParameters()[i], i) / (2 * PullStep);
      seedDerivative[i] = difference(plus.seed.getParameters()[i], minus.seed.getParameters()[i], i) / (2 * PullStep);
    }
    for (int i = 0; i < TrackParameterCount; ++i) {
      for (int j = 0; j < TrackParameterCount; ++j) {
        result.joint(i, j) += result.state.mft.getCovariances()(i, j) + lossDerivative[i] * lossDerivative[j];
        result.seedCov(i, j) += seedDerivative[i] * seedDerivative[j];
        result.residual(i, j) = result.joint(i, j) + result.joint(i + 5, j + 5) - result.joint(i, j + 5) - result.joint(i + 5, j);
      }
    }
    std::vector<double> delta(5);
    for (int i = 0; i < TrackParameterCount; ++i) {
      delta[i] = difference(result.state.mft.getParameters()[i], result.state.mch.getParameters()[i], i);
    }
    if (!chiDirect<5>(result.residual, delta, result.chiAll, result.rankAll) ||
        !chiDirect<4>(result.residual, delta, result.chiFour, result.rankFour) ||
        !chiDirect<2>(result.residual, delta, result.chiXY, result.rankXY)) {
      result.status = 4;
      return false;
    }
    o2::track::SMatrix55Sym mftCov, mchCov, seedCov;
    for (int i = 0; i < TrackParameterCount; ++i) {
      for (int j = 0; j <= i; ++j) {
        mftCov(i, j) = .5 * (result.joint(i, j) + result.joint(j, i));
        mchCov(i, j) = .5 * (result.joint(i + 5, j + 5) + result.joint(j + 5, i + 5));
        seedCov(i, j) = .5 * (result.seedCov(i, j) + result.seedCov(j, i));
      }
    }
    result.state.mft.setCovariances(mftCov);
    result.state.mch.setCovariances(mchCov);
    result.state.seed.setCovariances(seedCov);
    result.status = 0;
    return true;
  }
  static std::vector<double> flatten(const TMatrixD& matrix)
  {
    std::vector<double> out;
    out.reserve(matrix.GetNrows() * matrix.GetNcols());
    for (int i = 0; i < matrix.GetNrows(); ++i) {
      for (int j = 0; j < matrix.GetNcols(); ++j) {
        out.push_back(matrix(i, j));
      }
    }
    return out;
  }
};

// Bernoulli thinning is applied to final published training pairs, AFTER MC
// eligibility and donor-group normalization. Both matching planes share a draw.
inline double mixedPairSaveFraction = 1.;
inline uint64_t mixedPairSamplingSeed = 0;
struct MixedPairSampler {
  using Key = std::tuple<uint64_t, uint64_t, uint64_t, int64_t, int64_t, int64_t, int>;
  std::mt19937_64 random;
  std::map<Key, bool> decisions;
  uint64_t attempted = 0, saved = 0;
  MixedPairSampler()
  {
    std::random_device entropy;
    const auto seed = mixedPairSamplingSeed != 0 ? mixedPairSamplingSeed : (static_cast<uint64_t>(entropy()) << 32) ^ static_cast<uint64_t>(entropy());
    random.seed(seed);
  }
  template <typename Row>
  bool keep(Row const& row)
  {
    const Key key{row.mixingStreamId(), row.targetFrameId(), row.donorFrameId(), row.mchCollisionId(), row.mchTrackId(), row.mftTrackId(), row.mixedGroupIndex()};
    const auto found = decisions.find(key);
    if (found != decisions.end()) {
      return found->second;
    }
    const bool selected = mixedPairSaveFraction >= 1. || (mixedPairSaveFraction > 0. && std::bernoulli_distribution(mixedPairSaveFraction)(random));
    decisions.emplace(key, selected);
    ++attempted;
    saved += selected;
    return selected;
  }
  std::string metadata() const
  {
    std::ostringstream text;
    text << "{\"save_fraction\":" << mixedPairSaveFraction << ",\"sampling_seed\":" << mixedPairSamplingSeed
         << ",\"attempted_pairs\":" << attempted << ",\"saved_pairs\":" << saved
         << ",\"mixing_weight\":\"original donor normalization; divide by save_fraction for unbiased total counts\"}";
    return text.str();
  }
};

template <typename Row, typename Output>
inline void copyMixedResidual(const Row& row, Output& output, int groups = -1)
{
  output(row.mchTrackId(), row.mftTrackId(), row.mchCollisionId(), row.mftCollisionId(), row.matchingPlaneZ(), row.bzAtMftCenter(), row.mchX(), row.mchY(), row.mchPhi(), row.mchTanl(), row.mchInvQPt(), row.mftX(), row.mftY(), row.mftPhi(), row.mftTanl(), row.mftInvQPt(), row.deltaX(), row.deltaY(), row.deltaPhi(), row.deltaPhiWrapped(), row.deltaTanl(), row.deltaInvQPt(), row.deltaEta(), row.deltaPt(), row.chi2All(), row.chi2XYPhiTanl(), row.chi2XY(), row.isFinite(), row.absorberStepCm(), row.transportStatus(), row.donorMchTrackId(), row.donorCollisionId(), row.mixedGroupIndex(), groups >= 0 ? groups : row.numMixedGroups(), row.mixDeltaBc(), row.mixDeltaPhi(), row.mixDeltaR(), row.mixDeltaAttemptsRel(), row.mixDeltaZ(), row.targetMatchAttempts(), row.donorMatchAttempts(), row.targetMftMultiplicity(), row.donorMftMultiplicity(), row.targetMchMultiplicity(), row.donorMchMultiplicity(), row.targetCentrality(), row.donorCentrality(), row.targetCbMultiplicity(), row.donorCbMultiplicity(), (groups > 0 ? 1. / groups : row.mixingWeight()), row.targetFrameId(), row.donorFrameId(), row.mixingStreamId());
}

} // namespace o2::analysis::mftmch
// -----------------------------------------------------------------------------
// MftMchInference
// -----------------------------------------------------------------------------

namespace o2::analysis::mftmch
{
struct MatchingOutcome {
  int8_t stage = 0;
  int64_t mft = -1;
  double score = NAN, raw = NAN;
};
using MatchingOutcomes = std::map<int64_t, std::array<MatchingOutcome, 2>>;

// Public audit cursors stay empty in sparse mode. Internal transport carries
// the same rows to MC/sparse devices, independently of a pipeline's AOD writer.
inline bool sparseAuditStorage = false;
struct MatchingTreeStorage {
  bool mcTruth = true;
  bool candidateTruth = true;
  bool candidateRelations = true;
  bool mchTruth = true;
  bool mixedTruth = true;
  bool mchPopulation = true;
};
inline MatchingTreeStorage matchingTreeStorage;

template <typename T>
inline bool useInternalAuditStorage()
{
  if constexpr (std::is_same_v<T, o2::aod::MixedMftMchMatchTruths>) {
    // Mixed labels are training rows, retained in every output mode by default.
    return !(matchingTreeStorage.mcTruth && matchingTreeStorage.mixedTruth);
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingMchPopulation>) {
    return sparseAuditStorage || !matchingTreeStorage.mchPopulation;
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingCandidateTruths>) {
    return sparseAuditStorage || !(matchingTreeStorage.mcTruth && matchingTreeStorage.candidateTruth);
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingCandidateRelations>) {
    return sparseAuditStorage || !(matchingTreeStorage.mcTruth && matchingTreeStorage.candidateRelations);
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingMchTruths>) {
    return sparseAuditStorage || !(matchingTreeStorage.mcTruth && matchingTreeStorage.mchTruth);
  }
  return sparseAuditStorage;
}
template <typename T>
struct MatchingInternalTable;
template <>
struct MatchingInternalTable<o2::aod::MatchingMchPopulation> {
  using type = o2::aod::InternalMatchingMchPopulation;
};
template <>
struct MatchingInternalTable<o2::aod::MatchingCandidates> {
  using type = o2::aod::InternalMatchingCandidates;
};
template <>
struct MatchingInternalTable<o2::aod::MatchingMchOutcomes> {
  using type = o2::aod::InternalMatchingMchOutcomes;
};
template <>
struct MatchingInternalTable<o2::aod::SelectedMatchingMuons> {
  using type = o2::aod::InternalSelectedMatchingMuons;
};
template <>
struct MatchingInternalTable<o2::aod::MatchingDimuonVertices> {
  using type = o2::aod::InternalMatchingDimuonVertices;
};
template <>
struct MatchingInternalTable<o2::aod::MchDimuonsPV> {
  using type = o2::aod::InternalMchDimuonsPV;
};
template <>
struct MatchingInternalTable<o2::aod::MlMatchedMuonsPV> {
  using type = o2::aod::InternalMlMatchedMuonsPV;
};
template <>
struct MatchingInternalTable<o2::aod::MlMatchedDimuonsPV> {
  using type = o2::aod::InternalMlMatchedDimuonsPV;
};
template <>
struct MatchingInternalTable<o2::aod::MatchingCandidateTruths> {
  using type = o2::aod::InternalMatchingCandidateTruths;
};
template <>
struct MatchingInternalTable<o2::aod::MatchingCandidateRelations> {
  using type = o2::aod::InternalMatchingCandidateRelations;
};
template <>
struct MatchingInternalTable<o2::aod::MatchingMchTruths> {
  using type = o2::aod::InternalMatchingMchTruths;
};
template <>
struct MatchingInternalTable<o2::aod::McDimuonTruths> {
  using type = o2::aod::InternalMcDimuonTruths;
};
template <>
struct MatchingInternalTable<o2::aod::MixedMftMchMatchTruths> {
  using type = o2::aod::InternalMixedMftMchMatchTruths;
};
template <typename T>
constexpr uint8_t matchingAuditStorageBit()
{
  if constexpr (std::is_same_v<T, o2::aod::MatchingMchPopulation>) {
    return 1;
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingCandidateTruths>) {
    return 2;
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingCandidateRelations>) {
    return 4;
  } else if constexpr (std::is_same_v<T, o2::aod::MatchingMchTruths>) {
    return 8;
  }
  return 0;
}
inline uint8_t matchingAuditStorageMask()
{
  return (!matchingTreeStorage.mchPopulation ? 1 : 0) |
         (!(matchingTreeStorage.mcTruth && matchingTreeStorage.candidateTruth) ? 2 : 0) |
         (!(matchingTreeStorage.mcTruth && matchingTreeStorage.candidateRelations) ? 4 : 0) |
         (!(matchingTreeStorage.mcTruth && matchingTreeStorage.mchTruth) ? 8 : 0);
}
// Choose typed staged inputs before DPL builds its graph. Merely changing an
// InputSpec matcher does not change the table reference used by InputRecord.
template <uint8_t Mask = 0, typename Factory>
auto withMatchingAuditStorageMask(uint8_t mask, Factory&& factory)
{
  if (mask == Mask) {
    return factory.template operator()<Mask>();
  }
  if constexpr (Mask < 15) {
    return withMatchingAuditStorageMask<Mask + 1>(mask, std::forward<Factory>(factory));
  } else {
    LOGF(fatal, "Invalid matching audit storage mask");
    return factory.template operator()<0>();
  }
}
template <typename T, bool Sparse, uint8_t Mask = 0>
using MatchingAuditInput = std::conditional_t<Sparse || (Mask & matchingAuditStorageBit<T>()) != 0, typename MatchingInternalTable<T>::type, T>;
template <typename T>
struct MatchingAuditOutput : o2::framework::ProducesGroup {
  o2::framework::Produces<T> tree;
  o2::framework::Produces<typename MatchingInternalTable<T>::type> internal;
  template <typename... Args>
  void operator()(Args&&... args)
  {
    if (useInternalAuditStorage<T>()) {
      internal(std::forward<Args>(args)...);
    } else {
      tree(std::forward<Args>(args)...);
    }
  }
  int64_t lastIndex() { return useInternalAuditStorage<T>() ? internal.lastIndex() : tree.lastIndex(); }
};

struct MatchingMuonOutput : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::SelectedMatchingMuons> tree;
  o2::framework::Produces<o2::aod::InternalSelectedMatchingMuons> internal;
  o2::framework::Produces<o2::aod::AssignedMatchingMuons> assigned;
  o2::framework::Produces<o2::aod::RawMchDimuonPairs> mchPairs;
  o2::framework::Produces<o2::aod::RawGlobalDimuonPairs> globalPairs;
  template <typename... Args>
  void operator()(Args&&... args)
  {
    if (useInternalAuditStorage<o2::aod::SelectedMatchingMuons>()) {
      internal(std::forward<Args>(args)...);
    } else {
      tree(std::forward<Args>(args)...);
    }
  }
};

struct MatchingInferenceConfig : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> cfgTrainingMode{"cfgTrainingMode", true, "Training-sample tables; does NOT train a model"};
  o2::framework::Configurable<bool> cfgEvaluationMode{"cfgEvaluationMode", false, "Run frozen model, matching, single-muon and same/mixed dimuon sparse histograms"};
  o2::framework::Configurable<bool> cfgSaveSelectedMuonTable{"cfgSaveSelectedMuonTable", false, "Optional selected single-muon table; no candidate table in inference-only mode"};
  o2::framework::Configurable<bool> cfgEnableMlMatching{"cfgEnableMlMatching", true, "Enable ML matching in inference mode"};
  o2::framework::Configurable<bool> cfgEnableChi2Matching{"cfgEnableChi2Matching", true, "Enable chi2 matching in inference mode"};
  o2::framework::Configurable<bool> cfgEnableMchPV{"cfgEnableMchPV", true, "Enable MCH-only PV reference"};
  o2::framework::Configurable<std::string> cfgMlModel{"cfgMlModel", "NN", "NN or XGBoost; supervised/Same-Mixed selection follows model files and metadata"};
  o2::framework::Configurable<std::string> cfgMlModelDirectory{"cfgMlModelDirectory", "", "Local directory or alien:// directory containing NN/XGBoost .onnx and .json contracts"};
  o2::framework::Configurable<std::string> cfgMlModelFile{"cfgMlModelFile", "", "Local path, file:// URI or alien:// URI of ONNX; overrides directory; metadata must describe the same model"};
  o2::framework::Configurable<std::string> cfgMlMetadataFile{"cfgMlMetadataFile", "", "Local path, file:// URI or alien:// URI of feature/preprocessing JSON; overrides directory"};
  o2::framework::Configurable<std::string> cfgMlModelCacheDirectory{"cfgMlModelCacheDirectory", "", "Parent directory for a unique temporary Grid-model cache; empty uses the OS temporary directory"};
  o2::framework::Configurable<std::vector<std::string>> cfgMlInputFeatures{"cfgMlInputFeatures", std::vector<std::string>{}, "Ordered feature names; empty uses model contract; mismatches are fatal"};
  o2::framework::Configurable<double> cfgMlScoreMin{"cfgMlScoreMin", 0., "Require sigmoid(good logit)>=cut AFTER fixed one-to-one assignment; source score is not a truth probability"};
  o2::framework::Configurable<std::vector<double>> cfgMlScorePtEdges{"cfgMlScorePtEdges", std::vector<double>{}, "Optional MCH PV pT bin edges for ML score cuts; empty uses cfgMlScoreMin"};
  o2::framework::Configurable<std::vector<double>> cfgMlScorePtCuts{"cfgMlScorePtCuts", std::vector<double>{}, "ML score minima for pT bins; size=edges.size()-1; outside edges uses cfgMlScoreMin"};
  o2::framework::Configurable<std::vector<std::string>> cfgMlWorkingPointNames{"cfgMlWorkingPointNames", std::vector<std::string>{}, "Optional named ML dimuon working points from one common uncut assignment"};
  o2::framework::Configurable<std::vector<double>> cfgMlWorkingPointPtEdges{"cfgMlWorkingPointPtEdges", std::vector<double>{}, "Common MCH PV pT bin edges for named working points"};
  o2::framework::Configurable<std::vector<double>> cfgMlWorkingPointCuts{"cfgMlWorkingPointCuts", std::vector<double>{}, "Row-major score minima: working-point index, then pT bin"};
  o2::framework::Configurable<double> cfgSingleDcaXYMax{"cfgSingleDcaXYMax", -1., "Maximum sqrt(DCAx^2+DCAy^2) in cm after matching; inclusive; negative disables"};
  o2::framework::Configurable<double> cfgSingleDcaZMax{"cfgSingleDcaZMax", -1., "Maximum absolute DCAz in cm after matching; inclusive; negative disables"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaPtEdges{"cfgSingleDcaPtEdges", std::vector<double>{}, "Optional MCH PV pT edges for single-muon DCA limits; independent of ML score bins"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaXYMaxPtCuts{"cfgSingleDcaXYMaxPtCuts", std::vector<double>{}, "DCAxy maxima in cm per pT bin; negative disables in that bin; empty/outside edges uses cfgSingleDcaXYMax"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaZMaxPtCuts{"cfgSingleDcaZMaxPtCuts", std::vector<double>{}, "Absolute DCAz maxima in cm per pT bin; negative disables in that bin; empty/outside edges uses cfgSingleDcaZMax"};
  o2::framework::Configurable<bool> cfgApplySingleDcaCutsToMchPV{"cfgApplySingleDcaCutsToMchPV", false, "Also cut MCH-only no-Branson DCA at the PV z plane; differs from the matched MFT 3D DCA and DCAz is zero by construction"};
  o2::framework::Configurable<int> cfgMlBatchSize{"cfgMlBatchSize", 4096, "Maximum candidate batch for CPU ONNX inference"};
  o2::framework::Configurable<int> cfgOnnxThreads{"cfgOnnxThreads", 1, "ONNX CPU threads"};
  o2::framework::Configurable<std::string> cfgChi2Method{"cfgChi2Method", "combined", "Recomputed chi2: 77_all,77_four,77_xy,505_all,505_four,505_xy,combined; cfgComputeMatchingChi2=false uses AO2D official saved scalar, not parameter decompositions"};
  o2::framework::Configurable<double> cfgMatchingChi2Max{"cfgMatchingChi2Max", 40., "Maximum MFT-MCH matching chi2 selected by cfgChi2Method AFTER one-to-one assignment; combined applies to the summed chi2"};
  o2::framework::Configurable<double> cfgChi2Max{"cfgChi2Max", -1., "Deprecated alias of cfgMatchingChi2Max; -1 means unused; a nonnegative legacy value takes precedence"};
  o2::framework::Configurable<bool> cfgRequireMid{"cfgRequireMid", true, "Require MCH-MID standalone track (trackType=3) for every reference/matched method"};
  o2::framework::Configurable<bool> cfgEnableDimuonMixing{"cfgEnableDimuonMixing", true, "Mix already-reconstructed single muons, separately from MFT-MCH training mixing"};
  o2::framework::Configurable<bool> cfgEnableDimuonVertexing{"cfgEnableDimuonVertexing", false, "Fit two-track PCA after matching; failed fits reject dimuons only when enabled"};
  o2::framework::Configurable<bool> cfgRequireDimuonVertexSuccess{"cfgRequireDimuonVertexSuccess", true, "When vertexing is enabled, require fit success even without vertex cuts; false retains failed fits in the uncut mass spectrum; active vertex cuts always require a valid fit"};
  o2::framework::Configurable<bool> cfgDimuonVertexUseCovariance{"cfgDimuonVertexUseCovariance", false, "Covariance-weighted official forward fitter; requires measured MFT covariance for matched muons; false uses geometric distance"};
  o2::framework::Configurable<double> cfgDimuonPairDcaMin{"cfgDimuonPairDcaMin", -1., "Minimum 3D separation of the two fitted track points; inclusive; negative disables; cm"};
  o2::framework::Configurable<double> cfgDimuonPairDcaMax{"cfgDimuonPairDcaMax", -1., "Maximum 3D separation of the two fitted track points in cm; negative disables"};
  o2::framework::Configurable<double> cfgDimuonPcaXYMin{"cfgDimuonPcaXYMin", -1., "Minimum transverse PCA vertex displacement from PV; inclusive; negative disables; cm"};
  o2::framework::Configurable<double> cfgDimuonPcaXYMax{"cfgDimuonPcaXYMax", -1., "Maximum transverse PCA vertex displacement from PV in cm; negative disables"};
  o2::framework::Configurable<double> cfgDimuonPcaZMin{"cfgDimuonPcaZMin", -1., "Minimum absolute PCA vertex z minus PV z; inclusive; negative disables; cm"};
  o2::framework::Configurable<double> cfgDimuonPcaZMax{"cfgDimuonPcaZMax", -1., "Maximum absolute PCA vertex z minus PV z in cm; negative disables"};
  o2::framework::Configurable<double> cfgDimuonPca3DMin{"cfgDimuonPca3DMin", -1., "Minimum 3D PCA vertex displacement from PV; inclusive; negative disables; cm"};
  o2::framework::Configurable<double> cfgDimuonPca3DMax{"cfgDimuonPca3DMax", -1., "Maximum 3D PCA vertex displacement from PV in cm; negative disables"};
  o2::framework::Configurable<double> cfgDimuonVertexChi2Min{"cfgDimuonVertexChi2Min", -1., "Minimum official fitter objective per prong, covariance mode only; inclusive; negative disables"};
  o2::framework::Configurable<double> cfgDimuonVertexChi2Max{"cfgDimuonVertexChi2Max", -1., "Maximum official fitter objective per prong, covariance mode only; not a reduced chi2; negative disables"};
  o2::framework::Configurable<double> cfgDimuonLxyMin{"cfgDimuonLxyMin", -1.e30, "Inclusive signed Lxy lower bound in cm; default unbounded; negative is a physical bound"};
  o2::framework::Configurable<double> cfgDimuonLxyMax{"cfgDimuonLxyMax", 1.e30, "Inclusive signed Lxy upper bound in cm; default unbounded; negative is a physical bound"};
  o2::framework::Configurable<double> cfgDimuonPcaSignedZMin{"cfgDimuonPcaSignedZMin", -1.e30, "Inclusive signed PCA z minus PV z lower bound in cm; default unbounded"};
  o2::framework::Configurable<double> cfgDimuonPcaSignedZMax{"cfgDimuonPcaSignedZMax", 1.e30, "Inclusive signed PCA z minus PV z upper bound in cm; default unbounded"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonDcaPtEdges{"cfgDimuonDcaPtEdges", std::vector<double>{}, "Reconstructed dimuon pT edges for vertex/DCA/Lxy cuts in GeV/c; [low,high); outside edges uses scalar limits"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPairDcaMinPtCuts{"cfgDimuonPairDcaMinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for PairDca; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPairDcaMaxPtCuts{"cfgDimuonPairDcaMaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for PairDca; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaXYMinPtCuts{"cfgDimuonPcaXYMinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for PcaXY; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaXYMaxPtCuts{"cfgDimuonPcaXYMaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for PcaXY; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaZMinPtCuts{"cfgDimuonPcaZMinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for PcaZ; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaZMaxPtCuts{"cfgDimuonPcaZMaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for PcaZ; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPca3DMinPtCuts{"cfgDimuonPca3DMinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for Pca3D; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPca3DMaxPtCuts{"cfgDimuonPca3DMaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for Pca3D; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonVertexChi2MinPtCuts{"cfgDimuonVertexChi2MinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for VertexChi2; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonVertexChi2MaxPtCuts{"cfgDimuonVertexChi2MaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for VertexChi2; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonLxyMinPtCuts{"cfgDimuonLxyMinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for Lxy; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonLxyMaxPtCuts{"cfgDimuonLxyMaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for Lxy; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaSignedZMinPtCuts{"cfgDimuonPcaSignedZMinPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin min for PcaSignedZ; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaSignedZMaxPtCuts{"cfgDimuonPcaSignedZMaxPtCuts", std::vector<double>{}, "Per-dimuon-pT-bin max for PcaSignedZ; empty uses scalar limit; one value per bin; same units and disable convention as scalar"};
  o2::framework::Configurable<double> cfgDimuonVertexMaxR{"cfgDimuonVertexMaxR", 200., "Official fitter maximum seed radius from beam axis in cm"};
  o2::framework::Configurable<int> cfgDimuonVertexMaxIterations{"cfgDimuonVertexMaxIterations", 60, "Maximum PCA minimization iterations"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPairDcaAxis{"cfgDimuonPairDcaAxis", std::vector<double>{400, 0., 2.}, "Two-track DCA axis: nbins,min,max cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaXYAxis{"cfgDimuonPcaXYAxis", std::vector<double>{400, 0., 2.}, "PCA transverse displacement axis: nbins,min,max cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaZAxis{"cfgDimuonPcaZAxis", std::vector<double>{400, -20., 20.}, "Signed PCA z minus PV z axis: nbins,min,max cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonLxyAxis{"cfgDimuonLxyAxis", std::vector<double>{400, -2., 2.}, "Signed Lxy axis: nbins,min,max cm"};
  o2::framework::Configurable<int> cfgDimuonMixingDepth{"cfgDimuonMixingDepth", 5, "Maximum prior donor events per event; 0 disables"};
  o2::framework::Configurable<int> cfgDimuonPoolMaxEvents{"cfgDimuonPoolMaxEvents", 2000, "Maximum retained prior events across DFs"};
  o2::framework::Configurable<int> cfgDimuonPoolMaxMuons{"cfgDimuonPoolMaxMuons", 200000, "Maximum retained single-muon states across all methods"};
  o2::framework::Configurable<bool> cfgDimuonMixAcrossDF{"cfgDimuonMixAcrossDF", true, "Retain owned single-muon states across DFs; pool resets on run/field changes"};
  o2::framework::Configurable<double> cfgDimuonMixDeltaZ{"cfgDimuonMixDeltaZ", 1., "Maximum absolute PV z difference in cm; negative disables"};
  o2::framework::Configurable<double> cfgDimuonMixMftMultiplicityRel{"cfgDimuonMixMftMultiplicityRel", -1., "Relative MFT multiplicity difference / target; negative disables"};
  o2::framework::Configurable<double> cfgDimuonMixPvMultiplicityRel{"cfgDimuonMixPvMultiplicityRel", -1., "PVMults relative difference / target; requires centrality process; negative disables"};
  o2::framework::Configurable<double> cfgDimuonMixDeltaCentrality{"cfgDimuonMixDeltaCentrality", -1., "FT0C centrality difference; requires centrality process; negative disables"};
  o2::framework::Configurable<int64_t> cfgDimuonMixMinDeltaBC{"cfgDimuonMixMinDeltaBC", 0, "Minimum absolute BC separation between distinct mixed events"};
  o2::framework::Configurable<bool> cfgDimuonMixNormalize{"cfgDimuonMixNormalize", true, "Mixed-pair weight=1/number of compatible donor events (empty selected-method events included)"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonMassAxis{"cfgDimuonMassAxis", std::vector<double>{600, 0., 12.}, "Sparse mass: nbins,min,max GeV/c^2"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPtAxis{"cfgDimuonPtAxis", std::vector<double>{200, 0., 20.}, "Sparse dimuon pT: nbins,min,max GeV/c"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonRapidityAxis{"cfgDimuonRapidityAxis", std::vector<double>{80, -5., -1.}, "Sparse rapidity: nbins,min,max"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonDcaXYAxis{"cfgDimuonDcaXYAxis", std::vector<double>{400, 0., 2.}, "Sparse dimuon sqrt(DCAxy1^2+DCAxy2^2): nbins,min,max cm"};
  o2::framework::Configurable<std::vector<double>> cfgSinglePtAxis{"cfgSinglePtAxis", std::vector<double>{200, 0., 20.}, "Single pT axis"};
  o2::framework::Configurable<std::vector<double>> cfgSingleEtaAxis{"cfgSingleEtaAxis", std::vector<double>{100, -5., -1.}, "Single eta axis"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaAxis{"cfgSingleDcaAxis", std::vector<double>{400, -.2, .2}, "Single signed DCA x/y/z in cm"};
  o2::framework::Configurable<std::vector<double>> cfgMlScoreAxis{"cfgMlScoreAxis", std::vector<double>{100, 0., 1.}, "ML selection score axis"};
  o2::framework::Configurable<std::vector<double>> cfgChi2ScoreAxis{"cfgChi2ScoreAxis", std::vector<double>{200, 0., 100.}, "Chi2 score axis (overflow retained)"};
};
struct MatchingInference {
  static constexpr size_t BaseFeatureCount = 38;
  static constexpr size_t Chi2FeaturesPerPlane = 3;
  static constexpr size_t InputFeatureCount = BaseFeatureCount + 2 * Chi2FeaturesPerPlane;
  using Fwd = o2::track::TrackParCovFwd;
  struct PairVertex {
    bool success = false;
    std::array<double, 3> position{NAN, NAN, NAN};
    double dca = NAN, xy = NAN, dz = NAN, distance = NAN, objective = NAN, lxy = NAN;
  };
  struct Muon {
    int64_t mch = -1, mft = -1;
    int charge = 0;
    double score = NAN, raw = NAN, mchPtPV = NAN;
    bool passesMatchingCuts = true; // Baseline dimuons retain assigned tracks; selected singles/working points use this flag.
    MuonFourVector p;
    Fwd vertexTrack;
    std::array<double, 3> vertexPV{NAN, NAN, NAN};
    bool hasVertexTrack = false, hasVertexCovariance = false;
    std::array<double, 3> dca{NAN, NAN, NAN};
    double pdca = NAN, mftChi = NAN, mftClusters = NAN, mchClusters = NAN, midChi = NAN;
  };
  struct Event {
    uint64_t frame = 0, bc = 0;
    int64_t collision = -1;
    int run = -1, mftMult = 0, pvMult = -1;
    double z = 0., cent = NAN;
    std::array<std::vector<Muon>, 3> muons;
  };
  // Each analysis device owns its download directory; no shared cache names.
  // Declare the cache before the ONNX model so the session is destroyed first.
  struct ModelCache {
    std::filesystem::path directory;
    ~ModelCache()
    {
      if (!directory.empty()) {
        std::error_code error;
        std::filesystem::remove_all(directory, error);
      }
    }
  };
  std::function<void(Muon const&, Muon const&, int, PairVertex const&, bool)> savePairVertex;
  std::function<void(Muon const&, Muon const&, double, bool)> saveMchPair;
  std::function<void(Muon const&, Muon const&, double, int, int, PairVertex const&, double, std::string const&)> saveGlobalPair;
  std::unique_ptr<ModelCache> modelCache;
  std::unique_ptr<o2::ml::OnnxModel> model;
  std::vector<std::string> featureNames;
  std::vector<float> mean, stddev;
  double logitSign = -1., requiredAbsorberStep = 0.;
  bool requiredMchCharge = true;
  bool useStoredMatchingChi2 = false;
  std::string requestedChi2Method;
  std::string contract, contractContent;
  std::array<std::array<std::shared_ptr<TH1>, 6>, 3> vertexStatus;
  std::shared_ptr<TH1> status, poolQA;
  std::deque<Event> pool;
  size_t poolMuons = 0;
  int run = -1;
  double field = NAN;
  static bool validGoodLogitSign(double sign) { return std::isfinite(sign) && (sign == 1. || sign == -1.); }
  static float toGoodLogit(float rawLogit, double sign) { return static_cast<float>(sign * rawLogit); }
  static bool isCombinedChi2(const std::string& method)
  {
    return method == "77_all_plus_505_xy" || method == "combined";
  }
  static bool needsChi2Plane(const std::string& method, float planeZ)
  {
    return isCombinedChi2(method) || ((method.starts_with("505")) == (planeZ == DownstreamPlaneZ));
  }
  template <typename Row>
  static double matchingChi2(const std::string& method, const Row& upstream, const Row& downstream)
  {
    if (isCombinedChi2(method)) {
      const double upstreamChi2 = std::get<24>(upstream);
      const double downstreamChi2 = std::get<26>(downstream);
      if (!std::isfinite(upstreamChi2) || upstreamChi2 < 0. || !std::isfinite(downstreamChi2) || downstreamChi2 < 0.) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      return upstreamChi2 + downstreamChi2;
    }
    const auto& row = method.starts_with("505") ? downstream : upstream;
    return method.find("all") != std::string::npos    ? std::get<24>(row)
           : method.find("four") != std::string::npos ? std::get<25>(row)
                                                      : std::get<26>(row);
  }
  static std::string modelLocation(const std::string& file, const std::string& directory, const std::string& name)
  {
    if (!file.empty()) {
      return file;
    }
    if (directory.empty()) {
      LOGF(fatal, "Set cfgMlModelDirectory or both cfgMlModelFile and cfgMlMetadataFile when ML evaluation is enabled");
      return {};
    }
    return directory + (directory.back() == '/' ? "" : "/") + name;
  }
  std::string materializeModelFile(const std::string& location, const MatchingInferenceConfig& config, const char* name)
  {
    std::string path = location;
    if (path.starts_with("alien://")) {
      if ((!gGrid || !gGrid->IsConnected()) && !TGrid::Connect("alien://")) {
        LOGF(fatal, "Cannot connect to AliEn to read %s; check the Grid credentials", location.c_str());
        return {};
      }
      if (!modelCache) {
        auto cache = std::make_unique<ModelCache>();
        const auto parent = config.cfgMlModelCacheDirectory.value.empty()
                              ? std::filesystem::temp_directory_path()
                              : std::filesystem::path(config.cfgMlModelCacheDirectory.value);
        std::filesystem::create_directories(parent);
        cache->directory = parent / (std::string("o2-mft-mch-model-") + TUUID().AsString());
        if (!std::filesystem::create_directory(cache->directory)) {
          LOGF(fatal, "Cannot create unique model cache %s", cache->directory.string().c_str());
          return {};
        }
        modelCache = std::move(cache);
      }
      path = (modelCache->directory / name).string();
      if (!TFile::Cp(location.c_str(), path.c_str(), false)) {
        LOGF(fatal, "Cannot download model asset %s", location.c_str());
        return {};
      }
      LOGF(info, "Downloaded ML asset %s to %s", location.c_str(), path.c_str());
    } else if (path.starts_with("file://")) {
      path.erase(0, 7);
    } else if (path.find("://") != std::string::npos) {
      LOGF(fatal, "Unsupported ML asset URI %s; use a local path, file:// or alien://", location.c_str());
      return {};
    }
    std::ifstream input(path, std::ios::binary);
    if (!input || input.peek() == std::ifstream::traits_type::eof()) {
      LOGF(fatal, "ML asset is missing, empty or unreadable: %s", path.c_str());
      return {};
    }
    return path;
  }
  static std::vector<std::string> supportedFeatures() { return {"MFT77_X", "MFT77_Y", "MFT77_phi", "MFT77_tanLambda", "MFT77_qOverPt", "MCH_PV_pt", "MCH_PV_eta", "MCH_PV_p", "MCH77_X", "MCH77_Y", "MCH77_phi", "MCH77_tanLambda", "MCH77_qOverPt", "MCH505_X", "MCH505_Y", "MCH505_phi", "MCH505_tanLambda", "MCH505_qOverPt", "delta77_X", "delta77_Y", "delta77_phiWrapped", "delta77_tanLambda", "delta77_qOverPt", "delta77_eta", "delta77_pt", "delta505_X", "delta505_Y", "delta505_phiWrapped", "delta505_tanLambda", "delta505_qOverPt", "delta505_eta", "delta505_pt", "MFT77_charge", "MCH77_charge", "MFT505_X", "MFT505_Y", "R77", "R505", "chi2_5_77", "chi2_4_77", "chi2_XY_77", "chi2_5_505", "chi2_4_505", "chi2_XY_505"}; }
  std::vector<int> featureIndices;
  bool mlChi2Upstream = false;
  bool mlChi2Downstream = false;
  [[nodiscard]] bool usesMlChi2Plane(float planeZ) const { return planeZ == UpstreamPlaneZ ? mlChi2Upstream : mlChi2Downstream; }
  static o2::framework::AxisSpec axis(const std::vector<double>& v, const char* title)
  {
    if (v.size() != AxisConfigurationSize || !std::isfinite(v[0]) || v[0] < 1 || v[0] != std::floor(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]) || v[2] <= v[1]) {
      LOGF(fatal, "Invalid sparse axis %s", title);
    }
    return {static_cast<int>(v[0]), v[1], v[2], title};
  }
  static double matchingChi2Limit(MatchingInferenceConfig const& c)
  {
    return c.cfgChi2Max.value >= 0. ? c.cfgChi2Max.value : c.cfgMatchingChi2Max.value;
  }
  static bool validSingleDcaCuts(MatchingInferenceConfig const& c)
  {
    if (!std::isfinite(c.cfgSingleDcaXYMax.value) || !std::isfinite(c.cfgSingleDcaZMax.value)) {
      return false;
    }
    auto const& edges = c.cfgSingleDcaPtEdges.value;
    auto const& xy = c.cfgSingleDcaXYMaxPtCuts.value;
    auto const& z = c.cfgSingleDcaZMaxPtCuts.value;
    if (edges.empty()) {
      return xy.empty() && z.empty();
    }
    if (edges.size() < MinimumBinEdgeCount || (xy.empty() && z.empty()) ||
        (!xy.empty() && xy.size() + 1 != edges.size()) || (!z.empty() && z.size() + 1 != edges.size())) {
      return false;
    }
    for (size_t i = 0; i < edges.size(); ++i) {
      if (!std::isfinite(edges[i]) || ((i != 0u) && edges[i] <= edges[i - 1])) {
        return false;
      }
    }
    const auto finite = [](double value) { return std::isfinite(value); };
    return std::all_of(xy.begin(), xy.end(), finite) && std::all_of(z.begin(), z.end(), finite);
  }
  static std::array<double, 2> singleDcaLimits(MatchingInferenceConfig const& c, double mchPtPV)
  {
    std::array<double, 2> limits{c.cfgSingleDcaXYMax.value, c.cfgSingleDcaZMax.value};
    auto const& edges = c.cfgSingleDcaPtEdges.value;
    if (edges.size() >= MinimumBinEdgeCount && mchPtPV >= edges.front() && mchPtPV < edges.back()) {
      const size_t bin = std::upper_bound(edges.begin(), edges.end(), mchPtPV) - edges.begin() - 1;
      if (!c.cfgSingleDcaXYMaxPtCuts.value.empty()) {
        limits[0] = c.cfgSingleDcaXYMaxPtCuts.value[bin];
      }
      if (!c.cfgSingleDcaZMaxPtCuts.value.empty()) {
        limits[1] = c.cfgSingleDcaZMaxPtCuts.value[bin];
      }
    }
    return limits;
  }
  static bool passesSingleDcaCuts(MatchingInferenceConfig const& c, Muon const& muon, int method)
  {
    if (method == 0 && !c.cfgApplySingleDcaCutsToMchPV.value) {
      return true;
    }
    if (!c.cfgSingleDcaPtEdges.value.empty() && !std::isfinite(muon.mchPtPV)) {
      return false;
    }
    const auto limits = singleDcaLimits(c, muon.mchPtPV);
    if (limits[0] >= 0. && (!std::isfinite(muon.dca[0]) || !std::isfinite(muon.dca[1]) || std::hypot(muon.dca[0], muon.dca[1]) > limits[0])) {
      return false;
    }
    return limits[1] < 0. || (std::isfinite(muon.dca[2]) && std::abs(muon.dca[2]) <= limits[1]);
  }
  void init(MatchingInferenceConfig& c, o2::framework::HistogramRegistry& h)
  {
    if (!c.cfgEvaluationMode) {
      return;
    }
    const std::array<std::array<double, 2>, 5> vertexRanges{{{c.cfgDimuonPairDcaMin.value, c.cfgDimuonPairDcaMax.value},
                                                             {c.cfgDimuonPcaXYMin.value, c.cfgDimuonPcaXYMax.value},
                                                             {c.cfgDimuonPcaZMin.value, c.cfgDimuonPcaZMax.value},
                                                             {c.cfgDimuonPca3DMin.value, c.cfgDimuonPca3DMax.value},
                                                             {c.cfgDimuonVertexChi2Min.value, c.cfgDimuonVertexChi2Max.value}}};
    for (auto const& limits : vertexRanges) {
      if (!std::isfinite(limits[0]) || !std::isfinite(limits[1]) || (limits[0] >= 0. && limits[1] >= 0. && limits[0] > limits[1])) {
        LOGF(fatal, "Invalid dimuon vertex min/max; require finite ordered active bounds");
      }
      if ((limits[0] >= 0. || limits[1] >= 0.) && !c.cfgEnableDimuonVertexing.value) {
        LOGF(fatal, "Enable dimuon vertexing before enabling vertex cuts");
      }
    }
    const std::array<std::array<double, 2>, 2> signedRanges{{{c.cfgDimuonLxyMin.value, c.cfgDimuonLxyMax.value},
                                                             {c.cfgDimuonPcaSignedZMin.value, c.cfgDimuonPcaSignedZMax.value}}};
    for (auto const& limits : signedRanges) {
      if (!std::isfinite(limits[0]) || !std::isfinite(limits[1]) || limits[0] > limits[1]) {
        LOGF(fatal, "Invalid signed vertex min/max; negative bounds are allowed, min must not exceed max");
      }
      if ((limits[0] > -1.e30 || limits[1] < 1.e30) && !c.cfgEnableDimuonVertexing.value) {
        LOGF(fatal, "Enable vertexing before enabling signed Lxy/PCA-z cuts");
      }
    }
    if (!std::isfinite(c.cfgDimuonVertexMaxR.value) || c.cfgDimuonVertexMaxR.value <= 0. || c.cfgDimuonVertexMaxIterations.value < 3 || ((c.cfgDimuonVertexChi2Min.value >= 0. || c.cfgDimuonVertexChi2Max.value >= 0.) && !c.cfgDimuonVertexUseCovariance.value)) {
      LOGF(fatal, "Invalid PCA fitter settings; a chi2 cut requires covariance mode");
    }
    if (!validDimuonDcaPtCuts(c)) {
      LOGF(fatal, "Invalid dimuon pT-dependent vertex limits: require finite increasing nonnegative edges, one value per bin, ordered effective bounds, vertexing enabled for active cuts, and covariance mode for chi2 cuts");
    }
    if (!validSingleDcaCuts(c)) {
      LOGF(fatal, "Invalid single-muon DCA limits: finite limits and increasing pT edges with one value per bin are required");
    }
    if (c.cfgMlScoreMin < 0. || c.cfgMlScoreMin > 1. || !std::isfinite(c.cfgMlScoreMin.value) || c.cfgMatchingChi2Max < 0. || !std::isfinite(c.cfgMatchingChi2Max.value) || (c.cfgChi2Max.value < 0. && c.cfgChi2Max.value != -1.) || !std::isfinite(c.cfgChi2Max.value) || c.cfgMlBatchSize < 1 || c.cfgDimuonMixingDepth < 0 || c.cfgDimuonPoolMaxEvents < 1 || c.cfgDimuonPoolMaxMuons < 1) {
      LOGF(fatal, "Invalid inference cuts/pool limits");
    }
    if (c.cfgChi2Max.value >= 0.) {
      LOGF(warning, "Deprecated cfgChi2Max overrides cfgMatchingChi2Max; use cfgMatchingChi2Max in new JSON configurations");
    }
    for (double const& v : {c.cfgDimuonMixDeltaZ.value, c.cfgDimuonMixMftMultiplicityRel.value, c.cfgDimuonMixPvMultiplicityRel.value, c.cfgDimuonMixDeltaCentrality.value}) {
      if (!std::isfinite(v)) {
        LOGF(fatal, "Nonfinite dimuon pool setting");
      }
    }
    auto const& scoreEdges = c.cfgMlScorePtEdges.value;
    auto const& scoreCuts = c.cfgMlScorePtCuts.value;
    if (!scoreEdges.empty() || !scoreCuts.empty()) {
      if (scoreEdges.size() < MinimumBinEdgeCount || scoreCuts.size() + 1 != scoreEdges.size()) {
        LOGF(fatal, "ML pT score cuts require N+1 edges for N cuts");
      }
      for (size_t i = 0; i < scoreEdges.size(); ++i) {
        if (!std::isfinite(scoreEdges[i]) || ((i != 0u) && scoreEdges[i] <= scoreEdges[i - 1])) {
          LOGF(fatal, "ML score pT edges must be finite and strictly increasing");
        }
      }
      for (double const& v : scoreCuts) {
        if (!std::isfinite(v) || v < 0. || v > 1.) {
          LOGF(fatal, "ML pT score minima must lie in [0,1]");
        }
      }
    }
    auto const& wpNames = c.cfgMlWorkingPointNames.value;
    auto const& wpEdges = c.cfgMlWorkingPointPtEdges.value;
    auto const& wpCuts = c.cfgMlWorkingPointCuts.value;
    if (!wpNames.empty()) {
      if (!c.cfgEnableMlMatching || c.cfgMlScoreMin.value != 0. || !c.cfgMlScorePtEdges.value.empty()) {
        LOGF(fatal, "Named working points require enabled ML with no primary score cut");
      }
      if (wpEdges.size() < MinimumBinEdgeCount || wpCuts.size() != wpNames.size() * (wpEdges.size() - 1)) {
        LOGF(fatal, "Working point cut matrix has wrong dimensions");
      }
      std::set<std::string> uniqueNames;
      for (auto const& name : wpNames) {
        if (name.empty() || name == "ML" || name == "MCHPV" || name == "Chi2" || !uniqueNames.insert(name).second ||
            name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos) {
          LOGF(fatal, "Invalid or duplicate working point name");
        }
      }
      for (size_t i = 0; i < wpEdges.size(); ++i) {
        if (!std::isfinite(wpEdges[i]) || ((i != 0u) && wpEdges[i] <= wpEdges[i - 1])) {
          LOGF(fatal, "Invalid working point pT edges");
        }
      }
      for (double const& cut : wpCuts) {
        if (!std::isfinite(cut) || cut < 0. || cut > 1.) {
          LOGF(fatal, "Invalid working point score cut");
        }
      }
    } else if (!wpEdges.empty() || !wpCuts.empty()) {
      LOGF(fatal, "Working point cuts require names");
    }
    if (c.cfgDimuonMixMinDeltaBC < 0 || c.cfgOnnxThreads < 1) {
      LOGF(fatal, "Invalid BC gap/ONNX thread count");
    }
    const std::set<std::string> methods{"77_all", "77_four", "77_xy", "505_all", "505_four", "505_xy", "77_all_plus_505_xy", "combined"};
    if (!methods.contains(c.cfgChi2Method.value)) {
      LOGF(fatal, "Unknown chi2 method");
    }
    if (c.cfgEnableMlMatching) {
      if (c.cfgMlModel.value != "NN" && c.cfgMlModel.value != "XGBoost") {
        LOGF(fatal, "Choose NN or XGBoost");
      }
      const auto metadataSource = modelLocation(c.cfgMlMetadataFile.value, c.cfgMlModelDirectory.value, c.cfgMlModel.value + ".json");
      const auto modelSource = modelLocation(c.cfgMlModelFile.value, c.cfgMlModelDirectory.value, c.cfgMlModel.value + ".onnx");
      const auto metadata = materializeModelFile(metadataSource, c, "metadata.json");
      const auto modelFile = materializeModelFile(modelSource, c, "model.onnx");
      boost::property_tree::ptree t;
      boost::property_tree::read_json(metadata, t);
      if (t.get<int>("schema") != 1 || t.get<std::string>("transform") != "signed_log1p_finite_then_standardize" || t.get<std::string>("model") != c.cfgMlModel.value) {
        LOGF(fatal, "Unsupported/wrong model contract");
      }
      contract = metadataSource;
      requiredAbsorberStep = t.get<double>("absorber_step_cm");
      requiredMchCharge = t.get<bool>("mch_charge_seed");
      for (auto const& v : t.get_child("features")) {
        featureNames.push_back(v.second.get_value<std::string>());
      }
      if (!c.cfgMlInputFeatures.value.empty() && c.cfgMlInputFeatures.value != featureNames) {
        LOGF(fatal, "JSON feature order differs from trained model contract");
      }
      const auto available = supportedFeatures();
      for (auto const& n : featureNames) {
        auto it = std::find(available.begin(), available.end(), n);
        if (it == available.end()) {
          LOGF(fatal, "Unsupported feature %s", n.c_str());
        }
        featureIndices.push_back(it - available.begin());
        if (n.starts_with("chi2_")) {
          mlChi2Upstream = mlChi2Upstream || n.ends_with("_77");
          mlChi2Downstream = mlChi2Downstream || n.ends_with("_505");
        }
      }
      for (auto const& v : t.get_child("mean")) {
        mean.push_back(v.second.get_value<float>());
      }
      for (auto const& v : t.get_child("std")) {
        stddev.push_back(v.second.get_value<float>());
      }
      if (mean.size() != FeatureMaskFactor * featureNames.size() || stddev.size() != mean.size()) {
        LOGF(fatal, "Wrong feature preprocessing dimensions");
      }
      for (size_t i = 0; i < mean.size(); ++i) {
        if (!std::isfinite(mean[i]) || !std::isfinite(stddev[i]) || stddev[i] <= 0.) {
          LOGF(fatal, "Invalid model normalization");
        }
      }
      logitSign = t.get<double>("good_logit_sign");
      if (!validGoodLogitSign(logitSign)) {
        LOGF(fatal, "Model metadata good_logit_sign must be +1 (correct logit) or -1 (mixed logit)");
      }
      model = std::make_unique<o2::ml::OnnxModel>();
      model->initModel(modelFile, false, c.cfgOnnxThreads);
      const auto shapes = model->getInputShapes();
      if (shapes.size() != 1 || shapes[0].size() != TensorRank || shapes[0][1] != FeatureMaskFactor * static_cast<int64_t>(featureNames.size())) {
        LOGF(fatal, "ONNX input shape disagrees with model contract");
      }
      std::ostringstream contractStream;
      boost::property_tree::write_json(contractStream, t);
      contractContent = contractStream.str();
      LOGF(info, "Inference model %s contract=%s; %zu features; sigmoid(good logit), NOT correct probability", c.cfgMlModel.value.c_str(), metadataSource.c_str(), featureNames.size());
    }
    using namespace o2::framework;
    const std::array<std::string, 3> labels{"MCHPV", "Chi2", "ML"};
    const std::array<std::string, 6> charges{"ULS", "LSpp", "LSmm", "MixULS", "MixLSpp", "MixLSmm"};
    for (int m = 0; m < ReconstructionMethodCount; ++m) {
      if (!enabled(c, m)) {
        continue;
      }
      if (c.cfgEnableDimuonVertexing.value) {
        for (int q = 0; q < PairChargeCategoryCount; ++q) {
          vertexStatus[m][q] = h.add<TH1>(labels[m] + "/VertexStatus/" + charges[q], "0=attempt,1=success,2=accepted,3=missing state/covariance,4=fit failure", o2::framework::HistType::kTH1D, {{5, -.5, 4.5}});
        }
      }
    }
    status = h.add<TH1>("Status", "0=MCH PV failure,1=matched DCA failure,2=ML input pair,3=ML assigned,4=ML accepted,5=chi2 assigned,6=chi2 accepted,7=event,8=pool reset,9=pool eviction", o2::framework::HistType::kTH1D, {{10, -.5, 9.5}});
    poolQA = h.add<TH1>("DimuonMixDonors", "Compatible donor events;N donors;Target events", o2::framework::HistType::kTH1D, {{201, -.5, 200.5}});
  }
  static bool enabled(MatchingInferenceConfig const& c, int m) { return m == 0 ? c.cfgEnableMchPV.value : m == 1 ? c.cfgEnableChi2Matching.value
                                                                                                                 : c.cfgEnableMlMatching.value; }
  static double sigmoid(double x) { return x >= 0. ? 1. / (1. + std::exp(-x)) : std::exp(x) / (1. + std::exp(x)); }
  // Rectangular Hungarian maximum-weight matching; dummy columns allow unpaired MCH.
  static std::vector<int> assign(const std::vector<std::pair<int64_t, int64_t>>& ids, const std::vector<double>& w)
  {
    std::map<int64_t, int> cm, fm;
    for (auto const& p : ids) {
      cm.emplace(p.first, 0);
      fm.emplace(p.second, 0);
    }
    int n = 0, k = 0;
    for (auto& p : cm) { // o2-linter: disable=const-ref-in-for-loop (Assigns dense row indices in place.)
      p.second = n++;
    }
    for (auto& p : fm) { // o2-linter: disable=const-ref-in-for-loop (Assigns dense column indices in place.)
      p.second = k++;
    }
    int cols = k + n;
    std::vector<std::vector<double>> cost(n, std::vector<double>(cols, 1e6));
    std::vector<std::vector<int>> lookup(n, std::vector<int>(k, -1));
    for (int i = 0; i < n; ++i) {
      for (int j = k; j < cols; ++j) {
        cost[i][j] = 0.;
      }
    }
    for (size_t t = 0; t < ids.size(); ++t) {
      int i = cm.at(ids[t].first), j = fm.at(ids[t].second);
      if (std::isfinite(w[t]) && -w[t] < cost[i][j]) {
        cost[i][j] = -w[t];
        lookup[i][j] = t;
      }
    }
    std::vector<double> u(n + 1), v(cols + 1);
    std::vector<int> p(cols + 1), way(cols + 1);
    for (int i = 1; i <= n; ++i) {
      p[0] = i;
      int j0 = 0;
      std::vector<double> minv(cols + 1, 1e20);
      std::vector<bool> used(cols + 1, false);
      do {
        used[j0] = true;
        int i0 = p[j0], j1 = 0;
        double delta = 1e20;
        for (int j = 1; j <= cols; ++j) {
          if (!used[j]) {
            double cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
            if (cur < minv[j]) {
              minv[j] = cur;
              way[j] = j0;
            }
            if (minv[j] < delta) {
              delta = minv[j];
              j1 = j;
            }
          }
        }
        for (int j = 0; j <= cols; ++j) {
          if (used[j]) {
            u[p[j]] += delta;
            v[j] -= delta;
          } else {
            minv[j] -= delta;
          }
        }
        j0 = j1;
      } while (p[j0] != 0);
      do {
        int j1 = way[j0];
        p[j0] = p[j1];
        j0 = j1;
      } while (j0 != 0);
    }
    std::vector<int> out;
    for (int j = 1; j <= k; ++j) {
      if ((p[j] != 0) && lookup[p[j] - 1][j - 1] >= 0) {
        out.push_back(lookup[p[j] - 1][j - 1]);
      }
    }
    std::sort(out.begin(), out.end());
    return out;
  }
  template <typename Row>
  static std::array<float, InputFeatureCount> features(Row const& a, Row const& b, Fwd const& pv)
  {
    std::array<float, InputFeatureCount> f{};
    f[0] = std::get<11>(a);
    f[1] = std::get<12>(a);
    f[2] = std::get<13>(a);
    f[3] = std::get<14>(a);
    f[4] = std::get<15>(a);
    f[5] = pv.getPt();
    f[6] = pv.getEta();
    f[7] = pv.getP();
    f[8] = std::get<6>(a);
    f[9] = std::get<7>(a);
    f[10] = std::get<8>(a);
    f[11] = std::get<9>(a);
    f[12] = std::get<10>(a);
    f[13] = std::get<6>(b);
    f[14] = std::get<7>(b);
    f[15] = std::get<8>(b);
    f[16] = std::get<9>(b);
    f[17] = std::get<10>(b);
    f[18] = std::get<16>(a);
    f[19] = std::get<17>(a);
    f[20] = std::get<19>(a);
    f[21] = std::get<20>(a);
    f[22] = std::get<21>(a);
    f[23] = std::get<22>(a);
    f[24] = std::get<23>(a);
    f[25] = std::get<16>(b);
    f[26] = std::get<17>(b);
    f[27] = std::get<19>(b);
    f[28] = std::get<20>(b);
    f[29] = std::get<21>(b);
    f[30] = std::get<22>(b);
    f[31] = std::get<23>(b);
    f[32] = std::get<15>(a) > 0. ? 1. : -1.;
    f[33] = std::get<10>(a) > 0. ? 1. : -1.;
    f[34] = std::get<11>(b);
    f[35] = std::get<12>(b);
    f[36] = std::hypot(f[18], f[19]);
    f[37] = std::hypot(f[25], f[26]);
    // Match training preprocessing: unavailable, negative and float-overflowing
    // chi2 values remain missing inputs, represented by NaN and a finite mask.
    const auto usableChi2 = [](double value) {
      const auto chi2 = static_cast<float>(value);
      return std::isfinite(chi2) && chi2 >= 0.f ? chi2 : std::numeric_limits<float>::quiet_NaN();
    };
    f[BaseFeatureCount] = usableChi2(std::get<24>(a));
    f[BaseFeatureCount + 1] = usableChi2(std::get<25>(a));
    f[BaseFeatureCount + 2] = usableChi2(std::get<26>(a));
    f[BaseFeatureCount + Chi2FeaturesPerPlane] = usableChi2(std::get<24>(b));
    f[BaseFeatureCount + Chi2FeaturesPerPlane + 1] = usableChi2(std::get<25>(b));
    f[BaseFeatureCount + Chi2FeaturesPerPlane + 2] = usableChi2(std::get<26>(b));
    return f;
  }
  template <size_t N>
  std::vector<float> predict(std::vector<std::array<float, N>> const& x, MatchingInferenceConfig const& c)
  {
    if (std::any_of(featureIndices.begin(), featureIndices.end(), [](int index) { return index < 0 || static_cast<size_t>(index) >= N; })) {
      LOGF(fatal, "Input feature buffer is shorter than the model contract");
    }
    std::vector<float> out;
    out.reserve(x.size());
    const size_t n = featureNames.size();
    for (size_t lo = 0; lo < x.size(); lo += c.cfgMlBatchSize.value) {
      const size_t hi = std::min(x.size(), lo + c.cfgMlBatchSize.value);
      std::vector<float> input((hi - lo) * 2 * n);
      for (size_t i = lo; i < hi; ++i) {
        for (size_t j = 0; j < n; ++j) {
          float v = x[i][featureIndices[j]];
          bool finite = std::isfinite(v);
          float z = finite ? std::copysign(std::log1p(std::abs(v)), v) : 0.f;
          input[(i - lo) * 2 * n + j] = (z - mean[j]) / stddev[j];
          input[(i - lo) * 2 * n + n + j] = ((finite ? 1.f : 0.f) - mean[n + j]) / stddev[n + j];
        }
      }
      auto z = model->evalModel<float>(input);
      if (z.size() != hi - lo) {
        LOGF(fatal, "ONNX output must be one scalar logit per pair");
      }
      for (auto const& v : z) {
        if (!std::isfinite(v)) {
          LOGF(fatal, "Nonfinite ONNX logit");
        }
        out.push_back(toGoodLogit(v, logitSign));
      }
    }
    return out;
  }
  static bool validVertexCovariance(Fwd const& track)
  {
    // This check is local to covariance-weighted vertexing, never matching chi2.
    double lower[TrackParameterCount][TrackParameterCount]{};
    auto const& cov = track.getCovariances();
    for (int i = 0; i < TrackParameterCount; ++i) {
      for (int j = 0; j <= i; ++j) {
        double value = cov(i, j);
        for (int k = 0; k < j; ++k) {
          value -= lower[i][k] * lower[j][k];
        }
        if (!std::isfinite(value) || (i == j && value <= 0.)) {
          return false;
        }
        lower[i][j] = i == j ? std::sqrt(value) : value / lower[j][j];
      }
    }
    return true;
  }
  static PairVertex fitPairVertex(MatchingInferenceConfig const& c, Muon const& a, Muon const& b, double bz, bool mixed)
  {
    PairVertex result;
    if (!a.hasVertexTrack || !b.hasVertexTrack || (c.cfgDimuonVertexUseCovariance.value && (!a.hasVertexCovariance || !b.hasVertexCovariance))) {
      return result;
    }
    auto second = b.vertexTrack;
    if (mixed) {
      // Match PV-relative coordinates, without rotating either track.
      second.setX(second.getX() + a.vertexPV[0] - b.vertexPV[0]);
      second.setY(second.getY() + a.vertexPV[1] - b.vertexPV[1]);
      second.setZ(second.getZ() + a.vertexPV[2] - b.vertexPV[2]);
    }
    if (!c.cfgDimuonVertexUseCovariance.value && std::abs(bz) < MinimumHelixField) {
      // The official helix seed is singular at zero field: solve two straight
      // trajectories analytically, retaining the same three-dimensional metric.
      const auto& first = a.vertexTrack;
      std::array<double, 3> u{std::cos(first.getPhi()) / first.getTgl(), std::sin(first.getPhi()) / first.getTgl(), 1.};
      std::array<double, 3> v{std::cos(second.getPhi()) / second.getTgl(), std::sin(second.getPhi()) / second.getTgl(), 1.};
      std::array<double, 3> r{first.getX() - second.getX(), first.getY() - second.getY(), first.getZ() - second.getZ()};
      double aa = 0., ab = 0., bb = 0., ar = 0., br = 0.;
      for (int i = 0; i < 3; ++i) {
        aa += u[i] * u[i];
        ab += u[i] * v[i];
        bb += v[i] * v[i];
        ar += u[i] * r[i];
        br += v[i] * r[i];
      }
      const double det = aa * bb - ab * ab;
      if (!std::isfinite(det) || det <= 1.e-12 * aa * bb) {
        return result; // No uniquely determined PCA for parallel trajectories.
      }
      const double t = (ab * br - bb * ar) / det, z = (aa * br - ab * ar) / det;
      std::array<double, 3> firstPoint{first.getX(), first.getY(), first.getZ()}, secondPoint{second.getX(), second.getY(), second.getZ()};
      double squared = 0.;
      for (int i = 0; i < 3; ++i) {
        firstPoint[i] += t * u[i];
        secondPoint[i] += z * v[i];
        result.position[i] = .5 * (firstPoint[i] + secondPoint[i]);
        squared += (firstPoint[i] - secondPoint[i]) * (firstPoint[i] - secondPoint[i]);
      }
      result.dca = std::sqrt(squared);
      result.xy = std::hypot(result.position[0] - a.vertexPV[0], result.position[1] - a.vertexPV[1]);
      result.dz = result.position[2] - a.vertexPV[2];
      result.distance = std::hypot(result.xy, result.dz);
      result.objective = squared / 4.; // Official absolute-distance objective per prong for N=2.
      result.success = std::isfinite(result.distance) && std::isfinite(result.dca) && std::hypot(result.position[0], result.position[1]) <= c.cfgDimuonVertexMaxR.value;
      return result;
    }
    // The forward fitter uses a positive cos(phi) in its derivatives. Rotate
    // both tracks and their measured covariance into their common azimuthal
    // bisector frame; this is a coordinate change, not a physical rotation.
    const double angle = std::atan2(std::sin(a.vertexTrack.getPhi()) + std::sin(second.getPhi()), std::cos(a.vertexTrack.getPhi()) + std::cos(second.getPhi()));
    const double cosine = std::cos(angle), sine = std::sin(angle);
    auto rotate = [&](Fwd track) {
      const double x = track.getX(), y = track.getY();
      track.setX(cosine * x + sine * y);
      track.setY(-sine * x + cosine * y);
      track.setPhi(std::atan2(std::sin(track.getPhi() - angle), std::cos(track.getPhi() - angle)));
      o2::track::SMatrix55Std jacobian = ROOT::Math::SMatrixIdentity();
      jacobian(0, 0) = jacobian(1, 1) = cosine;
      jacobian(0, 1) = sine;
      jacobian(1, 0) = -sine;
      track.setCovariances(ROOT::Math::Similarity(jacobian, track.getCovariances()));
      return track;
    };
    const auto first = rotate(a.vertexTrack);
    second = rotate(second);
    struct Fitter : o2::vertexing::FwdDCAFitterN<2> {
      using o2::vertexing::FwdDCAFitterN<2>::getTrackPos;
    };
    Fitter fitter;
    fitter.setBz(bz);
    fitter.setUseAbsDCA(!c.cfgDimuonVertexUseCovariance.value);
    fitter.setPropagateToPCA(false); // Do not request a covariance-based vertex update.
    fitter.setMaxR(c.cfgDimuonVertexMaxR.value);
    fitter.setMaxIter(c.cfgDimuonVertexMaxIterations.value);
    fitter.setMaxDXIni(-1.); // No hidden four-centimetre preselection.
    fitter.setMaxChi2(std::numeric_limits<float>::max());
    try {
      if (!fitter.process(first, second)) {
        return result;
      }
      auto const& p = fitter.getPCACandidate();
      auto const& firstPoint = fitter.getTrackPos(0);
      auto const& secondPoint = fitter.getTrackPos(1);
      for (int i = 0; i < 3; ++i) {
        result.position[i] = p[i];
      }
      // Use the fitted closest points, not tracks propagated to a common z plane.
      const double localX = result.position[0], localY = result.position[1];
      result.position[0] = cosine * localX - sine * localY;
      result.position[1] = sine * localX + cosine * localY;
      result.dca = std::hypot(std::hypot(firstPoint[0] - secondPoint[0], firstPoint[1] - secondPoint[1]), firstPoint[2] - secondPoint[2]);
      result.xy = std::hypot(result.position[0] - a.vertexPV[0], result.position[1] - a.vertexPV[1]);
      result.dz = result.position[2] - a.vertexPV[2];
      result.distance = std::hypot(result.xy, result.dz);
      result.objective = fitter.getChi2AtPCACandidate();
      result.success = std::isfinite(result.dca) && std::isfinite(result.distance) && std::isfinite(result.objective) && result.objective >= 0.;
    } catch (std::exception const&) {
      result.success = false;
    }
    return result;
  }
  static double signedLxy(double dx, double dy, double px, double py)
  {
    const double pt = std::hypot(px, py);
    return pt > 0. && std::isfinite(pt) ? (dx * px + dy * py) / pt : NAN;
  }
  using VertexLimits = std::array<std::array<double, 2>, 7>;
  static std::array<std::vector<double> const*, 14> dimuonPtCutArrays(MatchingInferenceConfig const& c)
  {
    return {&c.cfgDimuonPairDcaMinPtCuts.value,
            &c.cfgDimuonPairDcaMaxPtCuts.value,
            &c.cfgDimuonPcaXYMinPtCuts.value,
            &c.cfgDimuonPcaXYMaxPtCuts.value,
            &c.cfgDimuonPcaZMinPtCuts.value,
            &c.cfgDimuonPcaZMaxPtCuts.value,
            &c.cfgDimuonPca3DMinPtCuts.value,
            &c.cfgDimuonPca3DMaxPtCuts.value,
            &c.cfgDimuonVertexChi2MinPtCuts.value,
            &c.cfgDimuonVertexChi2MaxPtCuts.value,
            &c.cfgDimuonLxyMinPtCuts.value,
            &c.cfgDimuonLxyMaxPtCuts.value,
            &c.cfgDimuonPcaSignedZMinPtCuts.value,
            &c.cfgDimuonPcaSignedZMaxPtCuts.value};
  }
  static VertexLimits dimuonVertexLimits(MatchingInferenceConfig const& c, double dimuonPt)
  {
    VertexLimits limits{{{c.cfgDimuonPairDcaMin.value, c.cfgDimuonPairDcaMax.value},
                         {c.cfgDimuonPcaXYMin.value, c.cfgDimuonPcaXYMax.value},
                         {c.cfgDimuonPcaZMin.value, c.cfgDimuonPcaZMax.value},
                         {c.cfgDimuonPca3DMin.value, c.cfgDimuonPca3DMax.value},
                         {c.cfgDimuonVertexChi2Min.value, c.cfgDimuonVertexChi2Max.value},
                         {c.cfgDimuonLxyMin.value, c.cfgDimuonLxyMax.value},
                         {c.cfgDimuonPcaSignedZMin.value, c.cfgDimuonPcaSignedZMax.value}}};
    auto const& edges = c.cfgDimuonDcaPtEdges.value;
    if (edges.size() >= MinimumBinEdgeCount && std::isfinite(dimuonPt) && dimuonPt >= edges.front() && dimuonPt < edges.back()) {
      const size_t bin = std::upper_bound(edges.begin(), edges.end(), dimuonPt) - edges.begin() - 1;
      const auto arrays = dimuonPtCutArrays(c);
      for (size_t i = 0; i < arrays.size(); ++i) {
        if (!arrays[i]->empty()) {
          limits[i / 2][i % 2] = (*arrays[i])[bin];
        }
      }
    }
    return limits;
  }
  static bool validDimuonDcaPtCuts(MatchingInferenceConfig const& c)
  {
    auto const& edges = c.cfgDimuonDcaPtEdges.value;
    const auto arrays = dimuonPtCutArrays(c);
    if (edges.empty()) {
      return std::all_of(arrays.begin(), arrays.end(), [](auto values) { return values->empty(); });
    }
    if (edges.size() < MinimumBinEdgeCount || edges.front() < 0.) {
      return false;
    }
    for (size_t i = 0; i < edges.size(); ++i) {
      if (!std::isfinite(edges[i]) || (i > 0 && edges[i] <= edges[i - 1])) {
        return false;
      }
    }
    for (auto values : arrays) {
      if (!values->empty() && (values->size() + 1 != edges.size() || !std::all_of(values->begin(), values->end(), [](double value) { return std::isfinite(value); }))) {
        return false;
      }
    }
    for (size_t bin = 0; bin + 1 < edges.size(); ++bin) {
      const auto limits = dimuonVertexLimits(c, edges[bin]);
      for (size_t i = 0; i < limits.size(); ++i) {
        const double low = limits[i][0], high = limits[i][1];
        if (i < 5) {
          if ((low >= 0. && high >= 0. && low > high) || ((low >= 0. || high >= 0.) && !c.cfgEnableDimuonVertexing.value)) {
            return false;
          }
        } else if (low > high || ((low > -1.e30 || high < 1.e30) && !c.cfgEnableDimuonVertexing.value)) {
          return false;
        }
      }
      if ((limits[4][0] >= 0. || limits[4][1] >= 0.) && !c.cfgDimuonVertexUseCovariance.value) {
        return false;
      }
    }
    return true;
  }
  static bool passesPairVertex(MatchingInferenceConfig const& c, PairVertex const& v, double dimuonPt = NAN)
  {
    if (!c.cfgDimuonDcaPtEdges.value.empty() && !std::isfinite(dimuonPt)) {
      return false;
    }
    const auto limits = dimuonVertexLimits(c, dimuonPt);
    auto nonnegativeRange = [](double value, std::array<double, 2> const& bounds) {
      return std::isfinite(value) && (bounds[0] < 0. || value >= bounds[0]) && (bounds[1] < 0. || value <= bounds[1]);
    };
    if (!v.success) {
      bool active = false;
      for (size_t i = 0; i < limits.size(); ++i) {
        active = active || (i < 5 ? limits[i][0] >= 0. || limits[i][1] >= 0. : limits[i][0] > -1.e30 || limits[i][1] < 1.e30);
      }
      return !c.cfgRequireDimuonVertexSuccess.value && !active;
    }
    const bool lxyCut = limits[5][0] > -1.e30 || limits[5][1] < 1.e30;
    return nonnegativeRange(v.dca, limits[0]) && nonnegativeRange(v.xy, limits[1]) &&
           nonnegativeRange(std::abs(v.dz), limits[2]) && nonnegativeRange(v.distance, limits[3]) &&
           nonnegativeRange(v.objective, limits[4]) && v.dz >= limits[6][0] && v.dz <= limits[6][1] &&
           (!lxyCut || (std::isfinite(v.lxy) && v.lxy >= limits[5][0] && v.lxy <= limits[5][1]));
  }
  void fillPair(MatchingInferenceConfig const& c, Muon const& a, Muon const& b, double w, int method, int category, bool recordVertex = true, std::string const& workingPoint = "")
  {
    const auto p = a.p + b.p;
    if (!std::isfinite(muonInvariantMass(p)) || !std::isfinite(p.Pt()) || !std::isfinite(muonRapidity(p))) {
      return;
    }
    const double pairDcaXY = std::hypot(std::hypot(a.dca[0], a.dca[1]), std::hypot(b.dca[0], b.dca[1]));
    if (!std::isfinite(pairDcaXY)) {
      return;
    }
    PairVertex vertex;
    bool accepted = true;
    if (c.cfgEnableDimuonVertexing.value) {
      auto const& counter = vertexStatus[method][category];
      if (recordVertex) {
        counter->Fill(0., w);
      }
      const bool missing = !a.hasVertexTrack || !b.hasVertexTrack || (c.cfgDimuonVertexUseCovariance.value && (!a.hasVertexCovariance || !b.hasVertexCovariance));
      vertex = fitPairVertex(c, a, b, field, category >= 3);
      if (vertex.success) {
        vertex.lxy = signedLxy(vertex.position[0] - a.vertexPV[0], vertex.position[1] - a.vertexPV[1], p.Px(), p.Py());
      }
      accepted = passesPairVertex(c, vertex, p.Pt());
      if (recordVertex && category < 3) {
        savePairVertex(a, b, method, vertex, accepted);
      }
      if (recordVertex) {
        if (vertex.success) {
          counter->Fill(1., w);
        } else {
          counter->Fill(missing ? 3 : 4, w);
        }
        if (accepted) {
          counter->Fill(2., w);
        }
      }
    }
    if (method != 0) {
      // Baseline pools contain assignments BEFORE score/single-DCA cuts; save before vertex cuts too.
      saveGlobalPair(a, b, w, method, category, vertex, pairDcaXY, workingPoint);
    } else if (accepted) {
      // MCH-only reference semantics are unchanged by the GlobalMuon redesign.
      saveMchPair(a, b, w, category >= 3);
    }
  }
  static int signIndex(Muon const& a, Muon const& b) { return a.charge * b.charge < 0 ? 0 : a.charge > 0 ? 1
                                                                                                         : 2; }
  static bool passesWorkingPoint(MatchingInferenceConfig const& c, size_t workingPoint, Muon const& u)
  {
    if (!u.passesMatchingCuts) {
      return false;
    }
    auto const& edges = c.cfgMlWorkingPointPtEdges.value;
    if (!std::isfinite(u.mchPtPV) || u.mchPtPV < edges.front() || u.mchPtPV >= edges.back()) {
      return false;
    }
    size_t bin = std::upper_bound(edges.begin(), edges.end(), u.mchPtPV) - edges.begin() - 1;
    return u.score >= c.cfgMlWorkingPointCuts.value[workingPoint * (edges.size() - 1) + bin];
  }
  static bool rel(double a, double b, double limit) { return limit < 0. || (a == 0. ? b == 0. : std::abs(a - b) / std::abs(a) <= limit); }
  template <typename Collisions, typename BCs, typename Mchs, typename Mfts, typename Rows, typename Output, typename AssignedOutput, typename VertexOutput, typename Conversion>
  void process(MatchingInferenceConfig& c, uint64_t frame, Collisions const& collisions, BCs const& bcs, Mchs const& tracks, Mfts const& mfts, Rows const& upstream, Rows const& downstream, std::unordered_map<int64_t, bool> const& selected, std::unordered_map<int64_t, double> const& cents, std::unordered_map<int64_t, int> const& pvMults, Conversion& conversion, double bz, Output& table, AssignedOutput& assignedTable, VertexOutput& vertexTable, MatchingOutcomes& outcomes, std::map<int64_t, Fwd> const& canonicalPV, std::unordered_map<int64_t, Fwd> const& mftStates, std::unordered_map<int64_t, int64_t> const& mftCovarianceIds)
  {
    if (!c.cfgEvaluationMode || bcs.size() == 0) {
      return;
    }
    if (c.cfgEnableDimuonVertexing.value && c.cfgDimuonVertexUseCovariance.value && (c.cfgEnableMlMatching || c.cfgEnableChi2Matching) && mftCovarianceIds.empty()) {
      LOGF(fatal, "Covariance-weighted dimuon vertexing requires measured MFT covariance; select geometric mode for this input");
    }
    int currentRun = bcs.begin().runNumber();
    if (run != currentRun || field != bz) {
      pool.clear();
      poolMuons = 0;
      run = currentRun;
      field = bz;
      if (c.cfgEnableDimuonVertexing.value) {
        LOGF(info, "Dimuon vertex field: run=%d Bz=%.17g kG", currentRun, bz);
      }
      status->Fill(8);
    }
    if (!c.cfgDimuonMixAcrossDF) {
      pool.clear();
      poolMuons = 0;
    }
    savePairVertex = [&](Muon const& a, Muon const& b, int method, PairVertex const& vertex, bool accepted) {
      auto tr = tracks.rawIteratorAt(a.mch - tracks.offset());
      vertexTable(frame, tr.collisionId(), static_cast<int8_t>(method), a.mch, b.mch,
                  vertex.success, c.cfgDimuonVertexUseCovariance.value, accepted,
                  vertex.position[0], vertex.position[1], vertex.position[2], vertex.xy, vertex.dz, vertex.dca, vertex.objective, vertex.lxy);
    };
    saveMchPair = [&](Muon const& a, Muon const& b, double weight, bool mixed) {
      const auto p = a.p + b.p;
      table.mchPairs(a.mch, b.mch, mixed, weight,
                     muonInvariantMass(p), p.Pt(), muonRapidity(p), TVector2::Phi_0_2pi(p.Phi()), static_cast<int8_t>(a.charge + b.charge));
    };
    saveGlobalPair = [&](Muon const& a, Muon const& b, double weight, int method, int category, PairVertex const& vertex, double dcaXY, std::string const& workingPoint) {
      const auto p = a.p + b.p;
      table.globalPairs(a.mch, b.mch, a.mft, b.mft, static_cast<int8_t>(method), category >= 3,
                        static_cast<int8_t>(a.charge + b.charge), weight, std::vector<uint8_t>(workingPoint.begin(), workingPoint.end()),
                        muonInvariantMass(p), p.Pt(), muonRapidity(p),
                        vertex.success ? vertex.dca : NAN, vertex.success ? vertex.xy : NAN,
                        vertex.success ? vertex.dz : NAN, dcaXY, vertex.success ? vertex.lxy : NAN, vertex.success,
                        std::isfinite(a.score) && std::isfinite(b.score) ? (method == 1 ? std::max(a.score, b.score) : std::min(a.score, b.score)) : NAN);
    };
    std::map<int64_t, Fwd> pvCache;
    std::map<int64_t, Event> events;
    std::map<int64_t, int> mftMult;
    for (auto const& m : mfts) {
      if (m.has_collision()) {
        ++mftMult[m.collisionId()];
      }
    }
    auto baseMuon = [&](auto const& tr) {
      Muon u{};
      u.mch = tr.globalIndex();
      u.charge = tr.signed1Pt() > 0 ? 1 : -1;
      u.pdca = tr.pDca();
      u.midChi = tr.chi2MatchMCHMID();
      u.mchClusters = tr.nClusters();
      return u;
    };
    for (auto const& tr : tracks) {
      auto sel = selected.find(tr.globalIndex());
      if (sel == selected.end() || !sel->second || !tr.has_collision()) {
        continue;
      }
      if (c.cfgRequireMid && tr.trackType() != o2::aod::fwdtrack::MuonStandaloneTrack) {
        for (auto& outcome : outcomes[tr.globalIndex()]) { // o2-linter: disable=const-ref-in-for-loop (Update the audit outcome in place.)
          outcome.stage = 8;
        }
        continue;
      }
      auto& audit = outcomes[tr.globalIndex()];
      for (auto& outcome : audit) { // o2-linter: disable=const-ref-in-for-loop (Update the audit outcome in place.)
        outcome.stage = 1;
      }
      auto col = collisions.rawIteratorAt(tr.collisionId() - collisions.offset());
      auto state = canonicalPV.find(tr.globalIndex());
      if (state == canonicalPV.end()) {
        status->Fill(0);
        continue;
      }
      auto const& pv = state->second;
      if (!std::isfinite(pv.getP()) || pv.getP() <= 0.) {
        status->Fill(0);
        continue;
      }
      pvCache.emplace(tr.globalIndex(), pv);
      for (auto& outcome : audit) { // o2-linter: disable=const-ref-in-for-loop (Update the audit outcome in place.)
        outcome.stage = 0;
      }
      auto& e = events[tr.collisionId()];
      e.frame = frame;
      e.collision = tr.collisionId();
      e.run = currentRun;
      e.z = col.posZ();
      e.bc = bcs.rawIteratorAt(col.bcId() - bcs.offset()).globalBC();
      e.mftMult = mftMult[tr.collisionId()];
      auto ct = cents.find(tr.collisionId());
      if (ct != cents.end()) {
        e.cent = ct->second;
      }
      auto mt = pvMults.find(tr.collisionId());
      if (mt != pvMults.end()) {
        e.pvMult = mt->second;
      }
      if (c.cfgEnableMchPV) {
        auto u = baseMuon(tr);
        u.p = muonMomentum(pv.getPx(), pv.getPy(), pv.getPz(), muonMassFloat);
        auto d = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(tr, tr));
        if (o2::mch::TrackExtrap::extrapToVertexWithoutBranson(d, col.posZ())) {
          auto dd = conversion.mchToFwd(d);
          u.dca = {dd.getX() - col.posX(), dd.getY() - col.posY(), dd.getZ() - col.posZ()};
          if (c.cfgEnableDimuonVertexing.value) {
            u.vertexTrack = dd;
            u.vertexPV = {col.posX(), col.posY(), col.posZ()};
            u.hasVertexTrack = true;
            u.hasVertexCovariance = validVertexCovariance(dd);
          }
        }
        u.mchPtPV = pv.getPt();
        if (passesSingleDcaCuts(c, u, 0)) {
          e.muons[0].push_back(u);
        }
      }
    }
    std::map<int64_t, std::vector<int>> candidates;
    for (size_t i = 0; i < upstream.size(); ++i) {
      if (pvCache.contains(std::get<0>(upstream[i]))) {
        // A reassociated MFT may compete in several collisions: assign once per DF.
        candidates[0].push_back(i);
      }
    }
    auto matched = [&](int row, double score, double raw) -> std::optional<Muon> {
      auto const& r = upstream[row];
      auto tr = tracks.rawIteratorAt(std::get<0>(r) - tracks.offset());
      auto mf = mfts.rawIteratorAt(std::get<1>(r) - mfts.offset());
      auto col = collisions.rawIteratorAt(tr.collisionId() - collisions.offset());
      Muon u = baseMuon(tr);
      u.mft = mf.globalIndex();
      u.score = score;
      u.raw = raw;
      u.mchPtPV = pvCache.at(tr.globalIndex()).getPt();
      u.mftChi = mf.chi2();
      uint64_t packed = mf.mftClusterSizesAndTrackFlags();
      u.mftClusters = 0;
      for (int j = 0; j < JointParameterCount; ++j) {
        u.mftClusters += ((packed >> (6 * j)) & 63) != 0;
      }
      o2::track::SMatrix55Sym zero;
      for (int i = 0; i < TrackParameterCount; ++i) {
        for (int j = 0; j <= i; ++j) {
          zero(i, j) = 0.;
        }
      }
      Fwd native(mf.z(), SMatrix5(mf.x(), mf.y(), mf.phi(), mf.tgl(), mf.signed1Pt()), zero, mf.chi2());
      auto hybrid = o2::aod::fwdtrackutils::refitGlobalMuonCov(pvCache.at(tr.globalIndex()), native);
      std::array<double, 3> dca{NAN, NAN, NAN};
      if (std::abs(bz) > MinimumHelixField) {
        hybrid.propagateToDCAhelix(bz, {col.posX(), col.posY(), col.posZ()}, dca);
      } else {
        double inv = 1. / hybrid.getTgl(), tx = std::cos(hybrid.getPhi()) * inv, ty = std::sin(hybrid.getPhi()) * inv;
        double dz = -((hybrid.getX() - col.posX()) * tx + (hybrid.getY() - col.posY()) * ty + (hybrid.getZ() - col.posZ())) / (1. + tx * tx + ty * ty);
        hybrid.propagateToZlinear(hybrid.getZ() + dz);
        dca = {hybrid.getX() - col.posX(), hybrid.getY() - col.posY(), hybrid.getZ() - col.posZ()};
      }
      if (!std::isfinite(dca[0]) || !std::isfinite(dca[1]) || !std::isfinite(dca[2]) || !std::isfinite(hybrid.getPhi())) {
        status->Fill(1);
        return std::nullopt;
      }
      if (c.cfgEnableDimuonVertexing.value) {
        u.vertexPV = {col.posX(), col.posY(), col.posZ()};
        u.vertexTrack = hybrid;
        u.hasVertexTrack = true;
        if (c.cfgDimuonVertexUseCovariance.value && mftCovarianceIds.contains(mf.globalIndex())) {
          auto covarianceState = o2::aod::fwdtrackutils::refitGlobalMuonCov(pvCache.at(tr.globalIndex()), mftStates.at(mf.globalIndex()));
          // Same central parameters as the hybrid trajectory; propagate the measured errors.
          if (std::abs(bz) > MinimumHelixField) {
            covarianceState.propagateToZhelix(hybrid.getZ(), bz);
          } else {
            covarianceState.propagateToZlinear(hybrid.getZ());
          }
          covarianceState.setParameters(hybrid.getParameters());
          covarianceState.setZ(hybrid.getZ());
          u.vertexTrack = covarianceState;
          u.hasVertexCovariance = validVertexCovariance(u.vertexTrack);
        }
      }
      u.dca = dca;
      const double p = pvCache.at(tr.globalIndex()).getP(), pt = p / std::sqrt(1. + hybrid.getTgl() * hybrid.getTgl());
      u.p = muonMomentum(pt * std::cos(hybrid.getPhi()), pt * std::sin(hybrid.getPhi()), pt * hybrid.getTgl(), muonMassFloat);
      return u;
    };
    std::map<std::pair<int64_t, int64_t>, double> officialChi2;
    if (useStoredMatchingChi2 && c.cfgEnableChi2Matching) {
      for (auto const& tr : tracks) {
        const bool globalMuon = tr.trackType() == o2::aod::fwdtrack::GlobalMuonTrack;
        const bool globalForward = !c.cfgRequireMid && tr.trackType() == o2::aod::fwdtrack::GlobalForwardTrack;
        if (!globalMuon && !globalForward) {
          continue;
        }
        const int64_t mchId = tr.matchMCHTrackId(), mftId = tr.matchMFTTrackId();
        const double chi2 = tr.chi2MatchMCHMFT();
        if (!pvCache.contains(mchId) || mftId < static_cast<int64_t>(mfts.offset()) || mftId - static_cast<int64_t>(mfts.offset()) >= mfts.size() || !std::isfinite(chi2) || chi2 < 0.) {
          continue;
        }
        auto [it, inserted] = officialChi2.emplace(std::make_pair(mchId, mftId), chi2);
        if (!inserted) {
          it->second = std::min(it->second, chi2);
        }
      }
    }
    size_t usableOfficialPairs = 0;
    for (auto const& entry : candidates) {
      auto const& rows = entry.second;
      std::vector<std::pair<int64_t, int64_t>> ids;
      ids.reserve(rows.size());
      for (int const& i : rows) {
        ids.emplace_back(std::get<0>(upstream[i]), std::get<1>(upstream[i]));
      }
      for (int method = 1; method <= MlMethodIndex; ++method) {
        if (!enabled(c, method)) {
          continue;
        }
        std::vector<double> w(rows.size(), NAN), score(rows.size(), NAN), raw(rows.size(), NAN);
        if (method == MlMethodIndex) {
          std::vector<std::array<float, InputFeatureCount>> x;
          x.reserve(rows.size());
          for (int const& i : rows) {
            x.push_back(features(upstream[i], downstream[i], pvCache.at(std::get<0>(upstream[i]))));
          }
          auto z = predict(x, c);
          for (size_t j = 0; j < rows.size(); ++j) {
            raw[j] = z[j];
            score[j] = w[j] = sigmoid(z[j]);
          }
          status->Fill(2, rows.size());
        } else {
          for (size_t j = 0; j < rows.size(); ++j) {
            double chi = NAN;
            if (useStoredMatchingChi2) {
              auto found = officialChi2.find(ids[j]);
              if (found != officialChi2.end()) {
                chi = found->second;
                ++usableOfficialPairs;
              }
            } else {
              chi = matchingChi2(c.cfgChi2Method.value, upstream[rows[j]], downstream[rows[j]]);
            }
            if (std::isfinite(chi) && chi >= 0.) {
              raw[j] = score[j] = chi;
              w[j] = std::exp(-chi / 2.);
            }
          }
        }
        for (size_t j = 0; j < rows.size(); ++j) {
          auto& audit = outcomes[std::get<0>(upstream[rows[j]])][method - 1];
          if (audit.stage == 0) {
            audit.stage = 2;
          }
          if (std::isfinite(w[j])) {
            audit.stage = 3;
          }
        }
        auto selectedRows = assign(ids, w);
        status->Fill(method == 1 ? 5 : 3, selectedRows.size());
        for (int const& j : selectedRows) {
          auto& audit = outcomes[std::get<0>(upstream[rows[j]])][method - 1];
          audit.mft = std::get<1>(upstream[rows[j]]);
          audit.score = score[j];
          audit.raw = raw[j];
          audit.stage = 4;
          // Save the assigned track before score and DCA cuts. No reassignment.
          auto u = matched(rows[j], score[j], raw[j]);
          if (u) {
            assignedTable(frame, std::get<2>(upstream[rows[j]]), u->mch, u->mft, static_cast<int8_t>(method), u->score, u->raw,
                          static_cast<int8_t>(u->charge), u->p.Px(), u->p.Py(), u->p.Pz(), u->dca[0], u->dca[1], u->dca[2]);
          }
          double mlCut = c.cfgMlScoreMin.value;
          const double mchPtPV = pvCache.at(std::get<0>(upstream[rows[j]])).getPt();
          auto const& edges = c.cfgMlScorePtEdges.value;
          if (edges.size() >= MinimumBinEdgeCount && mchPtPV >= edges.front() && mchPtPV < edges.back()) {
            const size_t bin = std::upper_bound(edges.begin(), edges.end(), mchPtPV) - edges.begin() - 1;
            mlCut = c.cfgMlScorePtCuts.value[bin];
          }
          const bool passesScore = method == 1 ? score[j] <= matchingChi2Limit(c) : score[j] >= mlCut;
          if (u) {
            u->passesMatchingCuts = passesScore && passesSingleDcaCuts(c, *u, method);
            // Keep every reconstructable one-to-one assignment for baseline dimuons.
            events.at(std::get<2>(upstream[rows[j]])).muons[method].push_back(*u);
          }
          if (!passesScore) {
            continue;
          }
          audit.stage = u ? 6 : 5;
          if (u && u->passesMatchingCuts) {
            audit.stage = 7;
            status->Fill(method == 1 ? 6 : 4);
          }
        }
      }
    }
    if (useStoredMatchingChi2 && c.cfgEnableChi2Matching) {
      LOGF(info, "Official saved chi2: %zu unique pairs linked to eligible MCHs, %zu present in the task candidate set", officialChi2.size(), usableOfficialPairs);
      if (officialChi2.empty()) {
        LOGF(warning, "No valid AO2D official matching pairs linked to eligible MCH tracks in this DF; no chi2 matches are fabricated");
      } else if (usableOfficialPairs < officialChi2.size()) {
        LOGF(warning, "Official saved chi2: %zu pairs fail current MFT/collision/timing/propagation candidate requirements and are not selected", officialChi2.size() - usableOfficialPairs);
      }
    }
    for (auto const& [col, e] : events) {
      status->Fill(7);
      for (int m = 0; m < ReconstructionMethodCount; ++m) {
        if (!enabled(c, m)) {
          continue;
        }
        for (auto const& u : e.muons[m]) {
          if (c.cfgSaveSelectedMuonTable && u.passesMatchingCuts) {
            table(frame, col, u.mch, u.mft, static_cast<int8_t>(m), u.score, u.raw, static_cast<int8_t>(u.charge), u.p.Px(), u.p.Py(), u.p.Pz(), u.dca[0], u.dca[1], u.dca[2]);
          }
        }
        for (size_t i = 0; i < e.muons[m].size(); ++i) {
          for (size_t j = i + 1; j < e.muons[m].size(); ++j) {
            fillPair(c, e.muons[m][i], e.muons[m][j], 1., m, signIndex(e.muons[m][i], e.muons[m][j]));
          }
        }
      }
      for (size_t w = 0; w < c.cfgMlWorkingPointNames.value.size(); ++w) {
        auto const& muons = e.muons[2];
        for (size_t i = 0; i < muons.size(); ++i) {
          if (!passesWorkingPoint(c, w, muons[i])) {
            continue;
          }
          for (size_t j = i + 1; j < muons.size(); ++j) {
            if (passesWorkingPoint(c, w, muons[j])) {
              fillPair(c, muons[i], muons[j], 1., 2, signIndex(muons[i], muons[j]), false, c.cfgMlWorkingPointNames.value[w]);
            }
          }
        }
      }
      if (!c.cfgEnableDimuonMixing || c.cfgDimuonMixingDepth == 0) {
        continue;
      }
      std::vector<Event const*> donors;
      for (auto it = pool.rbegin(); it != pool.rend() && donors.size() < static_cast<size_t>(c.cfgDimuonMixingDepth.value); ++it) {
        auto const& d = *it;
        if (e.frame == d.frame && e.collision == d.collision) {
          continue;
        }
        auto deltaBC = e.bc > d.bc ? e.bc - d.bc : d.bc - e.bc;
        if (deltaBC < static_cast<uint64_t>(c.cfgDimuonMixMinDeltaBC.value)) {
          continue;
        }
        if (c.cfgDimuonMixDeltaZ >= 0. && std::abs(e.z - d.z) > c.cfgDimuonMixDeltaZ) {
          continue;
        }
        if (!rel(e.mftMult, d.mftMult, c.cfgDimuonMixMftMultiplicityRel)) {
          continue;
        }
        if (c.cfgDimuonMixPvMultiplicityRel >= 0. && (e.pvMult < 0 || d.pvMult < 0 || !rel(e.pvMult, d.pvMult, c.cfgDimuonMixPvMultiplicityRel))) {
          continue;
        }
        if (c.cfgDimuonMixDeltaCentrality >= 0. && (!std::isfinite(e.cent) || !std::isfinite(d.cent) || std::abs(e.cent - d.cent) > c.cfgDimuonMixDeltaCentrality)) {
          continue;
        }
        donors.push_back(&d);
      }
      poolQA->Fill(donors.size());
      double weight = c.cfgDimuonMixNormalize && !donors.empty() ? 1. / donors.size() : 1.;
      for (auto const& d : donors) {
        for (int m = 0; m < ReconstructionMethodCount; ++m) {
          if (enabled(c, m)) {
            for (auto const& u : e.muons[m]) {
              for (auto const& v : d->muons[m]) {
                fillPair(c, u, v, weight, m, 3 + signIndex(u, v));
              }
            }
          }
        }
      }
      for (auto const& donor : donors) {
        for (size_t w = 0; w < c.cfgMlWorkingPointNames.value.size(); ++w) {
          for (auto const& u : e.muons[2]) {
            if (passesWorkingPoint(c, w, u)) {
              for (auto const& v : donor->muons[2]) {
                if (passesWorkingPoint(c, w, v)) {
                  fillPair(c, u, v, weight, 2, 3 + signIndex(u, v), false, c.cfgMlWorkingPointNames.value[w]);
                }
              }
            }
          }
        }
      }
      size_t n = 0;
      for (auto const& v : e.muons) {
        n += v.size();
      }
      if (n > static_cast<size_t>(c.cfgDimuonPoolMaxMuons.value)) {
        continue;
      }
      poolMuons += n;
      pool.push_back(e);
      while (pool.size() > static_cast<size_t>(c.cfgDimuonPoolMaxEvents.value) || poolMuons > static_cast<size_t>(c.cfgDimuonPoolMaxMuons.value)) {
        for (auto const& v : pool.front().muons) {
          poolMuons -= v.size();
        }
        pool.pop_front();
        status->Fill(9);
      }
    }
    savePairVertex = {};
    saveMchPair = {};
    saveGlobalPair = {};
  }
};
} // namespace o2::analysis::mftmch

namespace o2::analysis::mftmch
{
template <typename Crossing>
[[nodiscard]] bool passesMinimumDiskCount(std::vector<double> const& diskZ, int minimum, Crossing const& crosses)
{
  if (minimum < 1 || static_cast<size_t>(minimum) > diskZ.size()) {
    return false;
  }
  const auto required = static_cast<size_t>(minimum);
  size_t crossed = 0;
  for (size_t i = 0; i < diskZ.size(); ++i) {
    crossed += crosses(diskZ[i]) ? 1 : 0;
    if (crossed + diskZ.size() - i - 1 < required) {
      return false;
    }
  }
  return crossed >= required;
}

struct LegacyState {
  int mRunNumber = 0;
  std::unordered_map<int, std::vector<int64_t>> mftTracksByCollision;
  std::unordered_map<int, std::vector<int64_t>> muonTracksByCollision;
  std::unordered_map<int, bool> matchingCollisions;
  std::unordered_map<int, bool> hasMftTracksByCollision;
  std::unordered_map<int, bool> hasMuonTracksByCollision;
  std::unordered_map<int, float> vertexZByCollision;
  std::unordered_map<int, int> mftMultiplicityByCollision;
};
// Local guard around the official algorithm. No fork of the association logic.
struct LocalMftAssociationConfig : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<float> cfgMftAssocNumSigma{"cfgMftAssocNumSigma", 4.f, "Official MFT time-compatibility sigma multiplier for internal association"};
  o2::framework::Configurable<float> cfgMftAssocTimeMargin{"cfgMftAssocTimeMargin", 0.f, "Official MFT time margin in ns for internal association"};
  o2::framework::Configurable<bool> cfgMftAssocIncludeUnassigned{"cfgMftAssocIncludeUnassigned", false, "Include ambiguous/unassigned MFT tracks in internal association"};
  o2::framework::Configurable<int> cfgMftAssocBcWindow{"cfgMftAssocBcWindow", 115, "Official one-sigma BC window for internal association"};
};
struct ProtectedMftAssociation {
  template <typename Collisions, typename BCs>
  static std::vector<std::pair<int64_t, int64_t>> run(Collisions const& collisions, BCs const& bcs,
                                                      aod::MFTTracks const& tracks, aod::AmbiguousMFTTracks const& ambiguous,
                                                      LocalMftAssociationConfig const& config)
  {
    std::vector<std::pair<int64_t, int64_t>> pairs;
    // Some AO2Ds have stale native collision IDs in a zero-collision DF.
    // Never enter the official track-navigation loop without its target tables.
    if (collisions.size() == 0 || tracks.size() == 0 || bcs.size() == 0) {
      return pairs;
    }
    // Bind independent views so the official algorithm sees exactly its native
    // table types, without changing the matching task's navigation bindings.
    aod::Collisions nativeCollisions{collisions.asArrowTableRef()};
    aod::BCs nativeBCs{bcs.asArrowTableRef()};
    aod::MFTTracks nativeTracks{{tracks.asArrowTableRef()}};
    aod::AmbiguousMFTTracks nativeAmbiguous{ambiguous.asArrowTableRef()};
    nativeCollisions.bindExternalIndices(&nativeBCs);
    nativeTracks.bindExternalIndices(&nativeCollisions);
    nativeAmbiguous.bindExternalIndices(&nativeBCs);
    // Some AO2Ds contain ambiguous slices pointing one row past the BC table.
    // The official forward branch dereferences them without the barrel guard.
    // Keep valid rows in their original order; unavailable BCs cannot associate.
    o2::soa::SelectionVector validAmbiguousRows;
    for (const auto& row : nativeAmbiguous) {
      const auto slice = row.bcIds();
      if (slice[0] >= 0 && slice[1] >= slice[0] && slice[1] < static_cast<int64_t>(nativeBCs.size())) {
        validAmbiguousRows.push_back(row.globalIndex());
      }
    }
    if (validAmbiguousRows.size() != static_cast<size_t>(nativeAmbiguous.size())) {
      LOGF(warning, "Ignoring %zu ambiguous MFT rows with unavailable BC slices", static_cast<size_t>(nativeAmbiguous.size()) - validAmbiguousRows.size());
    }
    o2::soa::Filtered<aod::AmbiguousMFTTracks> validAmbiguous{{nativeAmbiguous.asArrowTableRef()}, std::move(validAmbiguousRows)};
    validAmbiguous.bindExternalIndices(&nativeBCs);
    CollisionAssociation<false> associator;
    associator.setNumSigmaForTimeCompat(config.cfgMftAssocNumSigma.value);
    associator.setTimeMargin(config.cfgMftAssocTimeMargin.value);
    associator.setTrackSelectionOptionForStdAssoc(aod::track_association::TrackSelection::None);
    associator.setUsePvAssociation(aod::track_association::PVContrReassocOpt::Disabled);
    associator.setIncludeUnassigned(config.cfgMftAssocIncludeUnassigned.value);
    associator.setFillTableOfCollIdsPerTrack(false);
    associator.setBcWindow(config.cfgMftAssocBcWindow.value);
    auto collect = [&](int64_t collision, int64_t mft) { pairs.emplace_back(collision, mft); };
    auto discardReverseIndices = [](std::vector<int> const&) {};
    associator.runAssocWithTime(nativeCollisions, nativeTracks, nativeTracks, validAmbiguous, nativeBCs,
                                collect, discardReverseIndices);
    return pairs;
  }
};
struct MftMchMatchingTask {
  LocalMftAssociationConfig localAssociation;
  int matchingPlaneProcessCount() const
  {
    return static_cast<int>(doprocessMatchingPlane.value) + static_cast<int>(doprocessMatchingPlaneCentrality.value) +
           static_cast<int>(doprocessMatchingPlaneWithAssociation.value) + static_cast<int>(doprocessMatchingPlaneCentralityWithAssociation.value) +
           static_cast<int>(doprocessMatchingPlaneWithAssociationNoCov.value) + static_cast<int>(doprocessMatchingPlaneCentralityWithAssociationNoCov.value) +
           static_cast<int>(doprocessMatchingPlaneWithInternalAssociation.value) + static_cast<int>(doprocessMatchingPlaneCentralityWithInternalAssociation.value) +
           static_cast<int>(doprocessMatchingPlaneWithInternalAssociationNoCov.value) + static_cast<int>(doprocessMatchingPlaneCentralityWithInternalAssociationNoCov.value);
  }
  bool hasCentralityProcess() const
  {
    return doprocessMatchingPlaneCentrality.value || doprocessMatchingPlaneCentralityWithAssociation.value || doprocessMatchingPlaneCentralityWithAssociationNoCov.value ||
           doprocessMatchingPlaneCentralityWithInternalAssociation.value || doprocessMatchingPlaneCentralityWithInternalAssociationNoCov.value;
  }
  bool hasNoCovarianceProcess() const
  {
    return doprocessMatchingPlaneWithAssociationNoCov.value || doprocessMatchingPlaneCentralityWithAssociationNoCov.value ||
           doprocessMatchingPlaneWithInternalAssociationNoCov.value || doprocessMatchingPlaneCentralityWithInternalAssociationNoCov.value;
  }
  bool hasInternalAssociationProcess() const
  {
    return doprocessMatchingPlaneWithInternalAssociation.value || doprocessMatchingPlaneCentralityWithInternalAssociation.value ||
           doprocessMatchingPlaneWithInternalAssociationNoCov.value || doprocessMatchingPlaneCentralityWithInternalAssociationNoCov.value;
  }

  MatchingAuditOutput<aod::MatchingMchPopulation> mchPopulation;
  MatchingAuditOutput<aod::MatchingCandidates> candidateAudit;
  MatchingAuditOutput<aod::MatchingMchOutcomes> mchOutcomes;
  o2::framework::Configurable<bool> cfgSaveMatchingAudit{"cfgSaveMatchingAudit", true, "Save pre-assignment MCH population, all candidate identities and per-method failures; required for unchanged efficiency denominators"};
  o2::framework::Configurable<bool> cfgIncludeNativeMftCandidates{"cfgIncludeNativeMftCandidates", true, "In association processes, union official MFTTrackAssoc with native same-collision MFTs; deduplicate before timing/acceptance"};
  LegacyState legacy;
  MatchingInferenceConfig inferenceConfig;
  MatchingInference inference;
  o2::framework::HistogramRegistry inferenceHistos{"MatchingInference", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  MatchingMuonOutput selectedMuonTable;
  MatchingAuditOutput<aod::MatchingDimuonVertices> dimuonVertexTable;
  o2::framework::OutputObj<TNamed> inferenceMetadata{TNamed("InferenceMetadata", "Inference disabled")};

  o2::framework::Produces<o2::aod::RawMixedMatchPlaneResiduals> tableMixedResiduals;
  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<bool> cfgEnableEventMixing{"cfgEnableEventMixing", false, "Save separate mixed-event residuals at both planes"};
    o2::framework::Configurable<int> cfgMixingDepth{"cfgMixingDepth", 5, "Maximum distinct donor collisions per MCH; -1=no limit"};
    o2::framework::Configurable<int64_t> cfgMinDeltaBc{"cfgMinDeltaBc", 3564, "Minimum BC separation of mixed collisions"};
    o2::framework::Configurable<float> cfgMaxDeltaPhi{"cfgMaxDeltaPhi", static_cast<float>(o2::constants::math::PI / 10.), "Maximum difference of MCH position azimuths in radians; negative disables"};
    o2::framework::Configurable<float> cfgMaxDeltaR{"cfgMaxDeltaR", 10.f, "Maximum donor-target MCH radius difference in cm; negative disables"};
    o2::framework::Configurable<float> cfgMaxDeltaAttempts{"cfgMaxDeltaAttempts", .1f, "Maximum donor-target relative MFT candidate-count difference; negative disables"};
    o2::framework::Configurable<float> cfgMaxDeltaZ{"cfgMaxDeltaZ", 1.f, "Maximum absolute vertex-z difference in cm; negative disables"};
    o2::framework::Configurable<bool> cfgMixingSymmetricCuts{"cfgMixingSymmetricCuts", true, "Use absolute radius/relative-count differences; false preserves official signed cuts"};
    o2::framework::Configurable<float> cfgMixingMaxDeltaCentrality{"cfgMixingMaxDeltaCentrality", -1.f, "Absolute FT0C centrality difference in percentile points; -1 disables"};
    o2::framework::Configurable<double> cfgMixingMaxDeltaMchEtaPV{"cfgMixingMaxDeltaMchEtaPV", -1., "Maximum absolute donor-target MCH momentum-eta difference at PV; negative disables"};
  } mixingCuts;
  // Group multiplicity controls to stay within the DPL aggregate-field limit.
  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<double> cfgMixingMaxDeltaMftMultiplicityRel{"cfgMixingMaxDeltaMftMultiplicityRel", .1, "Absolute MFT count difference / target count; -1 disables"};
    o2::framework::Configurable<int> cfgMixingMaxDeltaMchMultiplicity{"cfgMixingMaxDeltaMchMultiplicity", 0, "Absolute MCH count difference; 0 exact, -1 disables"};
    o2::framework::Configurable<double> cfgMixingMaxDeltaCbMultiplicityRel{"cfgMixingMaxDeltaCbMultiplicityRel", -1., "Absolute PVMults multNTracksPV difference / target count; negative disables, zero target requires zero donor"};
    o2::framework::Configurable<bool> cfgMixingCbMultiplicityStrict{"cfgMixingCbMultiplicityStrict", false, "Require PV multiplicity relative difference strictly below the limit"};
    o2::framework::Configurable<bool> cfgMixingVertexZStrict{"cfgMixingVertexZStrict", false, "Require vertex-z difference strictly below the limit"};
    o2::framework::Configurable<bool> cfgMixingMftMultiplicityStrict{"cfgMixingMftMultiplicityStrict", false, "Require MFT multiplicity relative difference strictly below the limit"};
    o2::framework::Configurable<bool> cfgMixingMchRadiusStrict{"cfgMixingMchRadiusStrict", false, "Require MCH position-radius difference strictly below its limit"};
    o2::framework::Configurable<bool> cfgMixingMchEtaStrict{"cfgMixingMchEtaStrict", false, "Require MCH PV eta difference strictly below its limit"};
    o2::framework::Configurable<double> cfgMixingMaxDeltaMchPtRel{"cfgMixingMaxDeltaMchPtRel", -1., "Absolute donor-target MCH pT difference at PV / target pT; negative disables"};
    o2::framework::Configurable<float> cfgMixingMaxDeltaMchMomentumRel{"cfgMixingMaxDeltaMchMomentumRel", -1.f, "Absolute donor-target MCH total momentum difference at PV / target momentum; negative disables"};
    o2::framework::Configurable<float> cfgMixingMaxDeltaMchPhiPV{"cfgMixingMaxDeltaMchPhiPV", -1.f, "Maximum wrapped donor-target MCH momentum-phi difference at PV (rad); negative disables"};
  } mixingMultiplicityCuts;
  struct OwnedMixMft {
    int64_t id = -1, collision = -1;
    MftMchJointTransport::Fwd state;
    int64_t referenceBC = 0;
    double time = 0., timeResolution = 0.;
    [[nodiscard]] int64_t globalIndex() const { return id; }
    [[nodiscard]] int64_t collisionId() const { return collision; }
  };
  struct MixingPoolEntry {
    mixing::Event event;
    std::vector<OwnedMixMft> mfts;
  };
  struct OwnedMixEvent {
    uint64_t frame = 0;
    int64_t collision = -1;
    size_t stateCount = 0;
    std::vector<MixingPoolEntry> anchors;
  };
  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<bool> cfgMixingAcrossDF{"cfgMixingAcrossDF", false, "Also use owned donor events from previous DFs in this task instance"};
    o2::framework::Configurable<int> cfgMixingPoolMaxEvents{"cfgMixingPoolMaxEvents", 2000, "Maximum retained donor collisions across DFs (FIFO)"};
    o2::framework::Configurable<int> cfgMixingPoolMaxMftStates{"cfgMixingPoolMaxMftStates", 200000, "Maximum retained MFT state copies; oversized events are skipped intact"};
    o2::framework::Configurable<uint64_t> cfgMixingStreamId{"cfgMixingStreamId", 0, "Source stream ID for merged outputs; 0 generates a random nonzero 64-bit ID per task instance"};
    struct State {
      uint64_t stream = 0;
      mixing::History<OwnedMixEvent> history;
      uint64_t frame = 0;
    } state;
  } mixingPersistence;
  o2::framework::Produces<o2::aod::MftMchMixingPoolStates> tableMixingPoolStates;
  std::unordered_map<int64_t, int> mixingCbMultiplicities;
  std::unordered_map<int64_t, double> mixingCentralities;

  o2::framework::Produces<o2::aod::RawMatchPlaneResiduals> tableMatchPlaneResiduals;
  o2::framework::Produces<o2::aod::DownstreamMatchPlaneResiduals> tableDownstreamResiduals;
  o2::framework::Produces<o2::aod::MatchPlaneEventSels> tableMatchPlaneEventSels;
  o2::framework::Produces<o2::aod::MftMchPairCovariances> tablePairCovariances;
  o2::framework::Produces<o2::aod::MftMchMixingDiagnostics> tableMixingDiagnostics;
  o2::framework::Configurable<bool> cfgSaveMixingDiagnostics{"cfgSaveMixingDiagnostics", false, "Save pair covariance matrices for MCHs with mixing candidates"};
  o2::framework::Configurable<std::string> cfgEventSelection{"cfgEventSelection", "all", "Matching-plane output: all, int8 (Run 3 sel8), upc (MCH-BC FV0-veto candidates)"};
  o2::framework::Configurable<bool> cfgFilterBeforeMatching{"cfgFilterBeforeMatching", true, "For int8/upc: reject MCH tracks before propagation; false keeps the previous late-filter path for numerical comparisons"};
  o2::framework::Configurable<int> cfgUpcMinMchMid{"cfgUpcMinMchMid", 1, "Minimum collision-associated MCH-MID tracks in the same corrected BC"};
  o2::framework::Configurable<int> cfgUpcMaxMchMid{"cfgUpcMaxMchMid", -1, "Maximum MCH-MID tracks in the same corrected BC; -1 disables upper cut"};
  o2::framework::Configurable<int> cfgUpcMinMchTracks{"cfgUpcMinMchTracks", 2, "Minimum standalone MCH tracks in the MCH collision"};
  o2::framework::Configurable<int> cfgUpcMinMftTracks{"cfgUpcMinMftTracks", 2, "Minimum MFT tracks in the MCH collision, including tracks without covariance"};
  o2::framework::Configurable<int> cfgUpcMaxMchTracks{"cfgUpcMaxMchTracks", -1, "Maximum standalone MCH tracks in the MCH collision; -1 disables the upper cut"};
  o2::framework::Configurable<bool> cfgUpcRequireEqualTrackCounts{"cfgUpcRequireEqualTrackCounts", false, "Require equal MFT and standalone MCH counts in the MCH collision; MFT count includes tracks without covariance"};
  o2::framework::Configurable<int> cfgUpcFV0BcWindow{"cfgUpcFV0BcWindow", 0, "FIT veto within +/- this many BCs of corrected MCH BC (legacy option name)"};
  o2::framework::Configurable<float> cfgUpcMaxFV0Amplitude{"cfgUpcMaxFV0Amplitude", 100.f, "UDHelpers FV0A amplitude limit; negative disables FV0A"};
  o2::framework::Configurable<float> cfgUpcFV0ValidTime{"cfgUpcFV0ValidTime", 15.f, "UDHelpers maximum absolute FIT time in ns (legacy option name); out-of-time records fail"};
  o2::framework::Configurable<float> cfgUpcMaxFT0AAmplitude{"cfgUpcMaxFT0AAmplitude", -1.f, "UDHelpers FT0A amplitude limit; negative disables"};
  o2::framework::Configurable<float> cfgUpcMaxFT0CAmplitude{"cfgUpcMaxFT0CAmplitude", -1.f, "UDHelpers FT0C amplitude limit; negative disables"};
  o2::framework::Configurable<float> cfgUpcMaxFDDAAmplitude{"cfgUpcMaxFDDAAmplitude", -1.f, "UDHelpers FDDA charge limit; negative disables"};
  o2::framework::Configurable<float> cfgUpcMaxFDDCAmplitude{"cfgUpcMaxFDDCAmplitude", -1.f, "UDHelpers FDDC charge limit; negative disables"};
  o2::framework::Configurable<bool> cfgUpcRequireUnlikeSign{"cfgUpcRequireUnlikeSign", true, "Require at least one positive and one negative MCH"};
  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<bool> cfgRequireMchMftAcceptance{"cfgRequireMchMftAcceptance", false, "Require PV-propagated MCH to cross the configured minimum number of nominal MFT disk envelopes"};
    o2::framework::Configurable<int> cfgMchMftAcceptanceMinDisks{"cfgMchMftAcceptanceMinDisks", 5, "Minimum crossed nominal MFT disks when acceptance selection is enabled; 1..number of configured distinct disk z positions"};
    o2::framework::Configurable<bool> cfgRequireMchPVEta{"cfgRequireMchPVEta", true, "Require propagated MCH eta at its collision PV within the configured strict range"};
    o2::framework::Configurable<double> cfgMchPVEtaMin{"cfgMchPVEtaMin", -3.6, "Strict lower limit on propagated MCH eta at its collision PV"};
    o2::framework::Configurable<double> cfgMchPVEtaMax{"cfgMchPVEtaMax", -2.5, "Strict upper limit on propagated MCH eta at its collision PV"};
    o2::framework::Configurable<double> cfgMchMftAcceptanceEtaMin{"cfgMchMftAcceptanceEtaMin", -3.6, "Nominal MFT envelope eta minimum for a vertex at (0,0,0)"};
    o2::framework::Configurable<double> cfgMchMftAcceptanceEtaMax{"cfgMchMftAcceptanceEtaMax", -2.5, "Nominal MFT envelope eta maximum for a vertex at (0,0,0)"};
    o2::framework::Configurable<std::vector<double>> cfgMchMftAcceptanceDiskZ{"cfgMchMftAcceptanceDiskZ", {-46., -49.3, -53.1, -68.4, -76.8}, "Nominal disk z positions (cm), from MFT Geometry.xml"};
    o2::framework::Configurable<bool> cfgPropMuon{"cfgPropMuon", true, "Propagate MCH to its collision PV before evaluating quality-cut eta and pt"};
    o2::framework::Configurable<std::string> cfgMuonCuts{"cfgMuonCuts", "muonQualityCuts", "PWGDQ cut: muonQualityCuts, muonQualityCuts5SigmaPDCA_Run3, muonQualityCuts10SigmaPDCA, or none"};
  } mchAcceptance;
  o2::framework::HistogramRegistry mchAcceptanceQA{"MchMftAcceptance", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // This is the nominal eta-defined geometric envelope, not a chip hit-efficiency
  // prediction. Propagate from the ACTUAL collision PV, including its x/y/z;
  // a fixed eta(PV) cut is a separate configurable selection, not this envelope.
  bool crossesNominalMftAcceptance(const o2::mch::TrackParam& mchAtPV)
  {
    return passesMinimumDiskCount(mchAcceptance.cfgMchMftAcceptanceDiskZ.value, mchAcceptance.cfgMchMftAcceptanceMinDisks.value, [&](double z) {
      if (!(mchAtPV.getZ() > z)) {
        return false;
      }
      auto atDisk = mchAtPV;
      // Official field integration; no absorber lies between the PV and MFT.
      // Only central parameters are needed for this geometric selection.
      if (!o2::mch::TrackExtrap::extrapToZ(atDisk, z)) {
        return false;
      }
      const double r = std::hypot(atDisk.getNonBendingCoor(), atDisk.getBendingCoor());
      const double rMin = std::abs(z / std::sinh(mchAcceptance.cfgMchMftAcceptanceEtaMin.value));
      const double rMax = std::abs(z / std::sinh(mchAcceptance.cfgMchMftAcceptanceEtaMax.value));
      return std::isfinite(r) && r > rMin && r < rMax;
    });
  }
  std::unique_ptr<AnalysisCut> mchQualityCut;
  int8_t matchingEventSelectionMode = 0;
  o2::framework::Configurable<double> cfgMftAbsorberStepCm{"cfgMftAbsorberStepCm", 0., "MFT absorber material-update maximum delta-z (cm): 0=whole absorber, positive=stepwise; ignored at z>=-90"};
  o2::framework::Configurable<bool> cfgMftUseMchMomentum{"cfgMftUseMchMomentum", true, "At z<=-505, replace MFT momentum magnitude using candidate MCH; propagate full joint errors"};
  o2::framework::Configurable<bool> cfgMftUseMchCharge{"cfgMftUseMchCharge", true, "Use reliable MCH charge for momentum-updated MFT; original MFT q/pt is saved"};
  o2::framework::Configurable<bool> cfgRequireMftCovariance{"cfgRequireMftCovariance", true, "Require stored MFT covariance; false admits all accepted MFT tracks and requires cfgComputeMatchingChi2=false"};
  o2::framework::Configurable<bool> cfgComputeMatchingChi2{"cfgComputeMatchingChi2", true, "Compute candidate matching chi2 and joint errors; false saves NaN training chi2 and uses official AO2D chi2MatchMCHMFT in chi2 evaluation"};
  o2::framework::Configurable<double> cfgPairJacobianScale{"cfgPairJacobianScale", 1., "Finite-difference scale for joint detector-error propagation; positive"};
  o2::aod::fwdtrackutils::GlobalMuonMatchingPlane globalMatchingPlane;
  // Keep buckets across DFs, exactly as global-muon-matcher::mMchTrackInfos.
  std::unordered_map<int64_t, int> globalMchIterationOrder;

  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<int> cfgMuonTaggingNCrossedMftPlanesLow{"cfgMuonTaggingNCrossedMftPlanesLow", 5, ""};
    o2::framework::Configurable<float> cfgMuonTaggingTrackChi2MchUp{"cfgMuonTaggingTrackChi2MchUp", 5.f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingPMchLow{"cfgMuonTaggingPMchLow", 0.0f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingPtMchLow{"cfgMuonTaggingPtMchLow", 0.7f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingEtaMchLow{"cfgMuonTaggingEtaMchLow", -3.6f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingEtaMchUp{"cfgMuonTaggingEtaMchUp", -2.5f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingRabsLow{"cfgMuonTaggingRabsLow", 17.6f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingRabsUp{"cfgMuonTaggingRabsUp", 89.5f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingPdcaUp{"cfgMuonTaggingPdcaUp", 4.f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingRadiusAtMftFrontLow{"cfgMuonTaggingRadiusAtMftFrontLow", 3.f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingRadiusAtMftFrontUp{"cfgMuonTaggingRadiusAtMftFrontUp", 9.f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingRadiusAtMftBackLow{"cfgMuonTaggingRadiusAtMftBackLow", 5.f, ""};
    o2::framework::Configurable<float> cfgMuonTaggingRadiusAtMftBackUp{"cfgMuonTaggingRadiusAtMftBackUp", 12.f, ""};
  } configMuonTagging;

  ////  Variables for matching method
  o2::framework::Configurable<int> cfgMatchingMethod{"cfgMatchingMethod", 0, ""};

  ////   Variables for selecting muon tracks
  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<float> cfgEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
    o2::framework::Configurable<float> cfgEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
    o2::framework::Configurable<float> cfgRabsLow1{"cfgRabsLow1", 17.6f, ""};
    o2::framework::Configurable<float> cfgRabsUp1{"cfgRabsUp1", 26.5f, ""};
    o2::framework::Configurable<float> cfgRabsLow2{"cfgRabsLow2", 26.5f, ""};
    o2::framework::Configurable<float> cfgRabsUp2{"cfgRabsUp2", 89.5f, ""};
    o2::framework::Configurable<float> cfgPdcaUp1{"cfgPdcaUp1", 594.f, ""};
    o2::framework::Configurable<float> cfgPdcaUp2{"cfgPdcaUp2", 324.f, ""};
    o2::framework::Configurable<float> cfgTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
    o2::framework::Configurable<float> cfgMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};
  } legacySelection;

  ////   Variables for selecting mft tracks
  o2::framework::Configurable<float> cfgEtaMftlow{"cfgEtaMftlow", -3.6f, ""};
  o2::framework::Configurable<float> cfgEtaMftup{"cfgEtaMftup", -2.5f, ""};
  o2::framework::Configurable<int> cfgTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  o2::framework::Configurable<float> cfgTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  ///    Variables to add preselection for the matching table
  o2::framework::Configurable<float> cfgPreselectMatchingX{"cfgPreselectMatchingX", 15.f, ""};
  o2::framework::Configurable<float> cfgPreselectMatchingY{"cfgPreselectMatchingY", 15.f, ""};

  ///    Variables to event mixing criteria
  o2::framework::Configurable<float> cfgSaveMixedMatchingParamsRate{"cfgSaveMixedMatchingParamsRate", 0.002f, ""};
  o2::framework::Configurable<int> cfgEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
  o2::framework::Configurable<float> cfgEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};

  ////   Variables for selecting tag muon
  o2::framework::Configurable<float> cfgTagMassWindowMin{"cfgTagMassWindowMin", 2.8f, ""};
  o2::framework::Configurable<float> cfgTagMassWindowMax{"cfgTagMassWindowMax", 3.3f, ""};

  ////   Variables for ccdb
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"}; // o2-linter: disable=name/configurable (Preserve the public ccdb-url JSON key.)
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  ////   Variables for Tag matching criteria
  o2::framework::Configurable<float> cfgSigmaXTagMuonCut{"cfgSigmaXTagMuonCut", 1.f, ""};
  o2::framework::Configurable<float> cfgMeanXTagMuonCut{"cfgMeanXTagMuonCut", 0.f, ""};
  o2::framework::Configurable<float> cfgSigmaYTagMuonCut{"cfgSigmaYTagMuonCut", 1.f, ""};
  o2::framework::Configurable<float> cfgMeanYTagMuonCut{"cfgMeanYTagMuonCut", 1.f, ""};

  o2::framework::Configurable<float> cfgSigmaEtaTagMuonCut{"cfgSigmaEtaTagMuonCut", 0.2f, ""};
  o2::framework::Configurable<float> cfgMeanEtaTagMuonCut{"cfgMeanEtaTagMuonCut", 0.f, ""};
  o2::framework::Configurable<float> cfgSigmaPhiTagMuonCut{"cfgSigmaPhiTagMuonCut", 0.2f, ""};
  o2::framework::Configurable<float> cfgMeanPhiTagMuonCut{"cfgMeanPhiTagMuonCut", 0.f, ""};

  template <typename MUON, typename Collision>
  class FindTagAndProbe
  {
   private:
    std::array<o2::dataformats::GlobalFwdTrack, 2> muonTrackAtPV{};

    MuonFourVector mDimuon;
    MUON muontrack1;
    MUON muontrack2;
    Collision collision;
    int tagIdx{-1}, probeIdx{-1};

    int16_t mQ{0};

    inline void fillCovarianceArray(MUON const& muontrack, std::array<float, 15>& cov) const
    {
      cov[0] = muontrack.cXX();
      cov[1] = muontrack.cXY();
      cov[2] = muontrack.cYY();
      cov[3] = muontrack.cPhiX();
      cov[4] = muontrack.cPhiY();
      cov[5] = muontrack.cPhiPhi();
      cov[6] = muontrack.cTglX();
      cov[7] = muontrack.cTglY();
      cov[8] = muontrack.cTglPhi();
      cov[9] = muontrack.cTglTgl();
      cov[10] = muontrack.c1PtX();
      cov[11] = muontrack.c1PtY();
      cov[12] = muontrack.c1PtPhi();
      cov[13] = muontrack.c1PtTgl();
      cov[14] = muontrack.c1Pt21Pt2();
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoPV(MUON const& muontrack) const
    {
      const double mz = muontrack.z();
      const double mchi2 = muontrack.chi2();
      const float mx = muontrack.x();
      const float my = muontrack.y();
      const float mphi = muontrack.phi();
      const float mtgl = muontrack.tgl();
      const float m1pt = muontrack.signed1Pt();

      std::array<float, 15> cov{};
      fillCovarianceArray(muontrack, cov);
      SMatrix5 tpars(mx, my, mphi, mtgl, m1pt);
      SMatrix55 tcovs(cov.begin(), cov.end());

      o2::track::TrackParCovFwd parcovmuontrack{mz, tpars, tcovs, mchi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      o2::globaltracking::MatchGlobalFwd mMatching;
      auto mchtrack = mMatching.FwdtoMCH(gtrack);

      o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
      o2::dataformats::GlobalFwdTrack extrapolatedMuonTrack;
      extrapolatedMuonTrack.setParameters(fwdtrack.getParameters());
      extrapolatedMuonTrack.setZ(fwdtrack.getZ());
      extrapolatedMuonTrack.setCovariances(fwdtrack.getCovariances());

      return extrapolatedMuonTrack;
    }

    inline void setTagAndProbe()
    {
      if (muontrack1.pt() > muontrack2.pt()) {
        tagIdx = 0;
        probeIdx = 1;
      } else {
        tagIdx = 1;
        probeIdx = 0;
      }
    }

   public:
    inline FindTagAndProbe(const MUON& muon1, const MUON& muon2, const Collision& coll) // o2-linter: disable=name/function-variable (Constructor spelling must match the class name.)
      : muonTrackAtPV(), muontrack1(muon1), muontrack2(muon2), collision(coll), mQ(muontrack1.sign() + muontrack2.sign())
    {
      setTagAndProbe();
    }

    void calcMuonPairAtPV()
    {
      muonTrackAtPV[0] = propagateMUONtoPV(muontrack1);
      muonTrackAtPV[1] = propagateMUONtoPV(muontrack2);
      MuonFourVector vMuon1, vMuon2;
      vMuon1 = muonMomentumFromPtEtaPhi(muonTrackAtPV[0].getPt(), muonTrackAtPV[0].getEta(), muonTrackAtPV[0].getPhi(), muonMassFloat);
      vMuon2 = muonMomentumFromPtEtaPhi(muonTrackAtPV[1].getPt(), muonTrackAtPV[1].getEta(), muonTrackAtPV[1].getPhi(), muonMassFloat);
      mDimuon = vMuon1 + vMuon2;
    }
    [[nodiscard]] inline int getTagMuonIndex() const { return tagIdx; }
    [[nodiscard]] inline int getProbeMuonIndex() const { return probeIdx; }
    [[nodiscard]] inline float getMass() const { return muonInvariantMass(mDimuon); }
    [[nodiscard]] inline float getPt() const { return mDimuon.Pt(); }
    [[nodiscard]] inline float getRap() const { return muonRapidity(mDimuon); }
    [[nodiscard]] inline int16_t getCharge() const { return mQ; }
    [[nodiscard]] inline const o2::dataformats::GlobalFwdTrack& getMuonAtPV(int idx) const { return muonTrackAtPV[idx]; }
  }; // end of class FindTagAndProbe

  template <typename MUON, typename MFT, typename Collision>
  class MatchingParamsML
  {
   private:
    MUON muontrack;
    MFT mfttrack;
    Collision collision;

    float mDX{0.f}, mDY{0.f}, mDPt{0.f}, mDPhi{0.f}, mDEta{0.f};
    float mGlobalMuonPtAtDCA{0.f}, mGlobalMuonEtaAtDCA{0.f}, mGlobalMuonPhiAtDCA{0.f}, mGlobalMuonDCAx{0.f}, mGlobalMuonDCAy{0.f}, mGlobalMuonQ{0.f};
    int mMatchingType;

    o2::field::MagneticField* fieldB = nullptr;
    o2::globaltracking::MatchGlobalFwd mMatching;

    inline o2::track::TrackParCovFwd propagateMFTtoMatchingPlane()
    {
      std::array<double, 15> covArr{0.0};
      SMatrix55 tmftcovs(covArr.begin(), covArr.end());

      SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      o2::track::TrackParCovFwd extrapolatedMftTrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};

      std::array<double, 3> propVec = {0.};
      float zPlane = 0.f;
      if (mMatchingType == MchFirstCluster) {
        propVec[0] = muontrack.x() - mfttrack.x();
        propVec[1] = muontrack.y() - mfttrack.y();
        propVec[2] = muontrack.z() - mfttrack.z();
        zPlane = muontrack.z();
      } else if (mMatchingType == EndOfAbsorber || mMatchingType == BEGINING_OF_ABSORBER) {
        auto extrapolatedMuonTrack = propagateMUONtoMatchingPlane();
        propVec[0] = extrapolatedMuonTrack.getX() - mfttrack.x();
        propVec[1] = extrapolatedMuonTrack.getY() - mfttrack.y();
        propVec[2] = extrapolatedMuonTrack.getZ() - mfttrack.z();
        zPlane = (mMatchingType == EndOfAbsorber) ? -505.f : -90.f;
      } else {
        zPlane = mfttrack.z();
      }

      std::array<double, 3> centerZ = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
      float bz = fieldB->getBz(centerZ.data());
      extrapolatedMftTrack.propagateToZ(zPlane, bz); // z in cm
      return extrapolatedMftTrack;
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoMatchingPlane()
    {
      std::array<float, 15> cov = {
        muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
        muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
        muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
        muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
        muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};

      SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
      SMatrix55 tcovs(cov.begin(), cov.end());
      double chi2 = muontrack.chi2();

      o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      auto mchtrack = mMatching.FwdtoMCH(gtrack);

      if (mMatchingType == MFT_LAST_CLUSTR) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, mfttrack.z());
      } else if (mMatchingType == EndOfAbsorber) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, -505.);
      } else if (mMatchingType == BEGINING_OF_ABSORBER) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack, -90.);
      }

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);

      o2::dataformats::GlobalFwdTrack extrapolatedMuonTrack;
      extrapolatedMuonTrack.setParameters(fwdtrack.getParameters());
      extrapolatedMuonTrack.setZ(fwdtrack.getZ());
      extrapolatedMuonTrack.setCovariances(fwdtrack.getCovariances());
      return extrapolatedMuonTrack;
    }

    inline o2::track::TrackParCovFwd propagateMFTtoDCA()
    {
      std::array<double, 15> covArr{0.0};
      SMatrix55 tmftcovs(covArr.begin(), covArr.end());

      SMatrix5 tmftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      o2::track::TrackParCovFwd extrapolatedMftTrack{mfttrack.z(), tmftpars, tmftcovs, mfttrack.chi2()};

      std::array<double, 3> propVec = {};
      propVec[0] = collision.posX() - mfttrack.x();
      propVec[1] = collision.posY() - mfttrack.y();
      propVec[2] = collision.posZ() - mfttrack.z();

      std::array<double, 3> centerZ = {mfttrack.x() + propVec[0] / 2., mfttrack.y() + propVec[1] / 2., mfttrack.z() + propVec[2] / 2.};
      float bz = fieldB->getBz(centerZ.data());
      extrapolatedMftTrack.propagateToZ(collision.posZ(), bz); // z in cm
      return extrapolatedMftTrack;
    }

    inline o2::dataformats::GlobalFwdTrack propagateMUONtoPV()
    {
      std::array<float, 15> cov = {
        muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
        muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
        muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
        muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
        muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};

      SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
      SMatrix55 tcovs(cov.begin(), cov.end());
      double chi2 = muontrack.chi2();

      o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

      o2::dataformats::GlobalFwdTrack gtrack;
      gtrack.setParameters(tpars);
      gtrack.setZ(parcovmuontrack.getZ());
      gtrack.setCovariances(tcovs);

      auto mchtrack = mMatching.FwdtoMCH(gtrack);
      o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

      auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
      o2::dataformats::GlobalFwdTrack extrapolatedMuonTrack;
      extrapolatedMuonTrack.setParameters(fwdtrack.getParameters());
      extrapolatedMuonTrack.setZ(fwdtrack.getZ());
      extrapolatedMuonTrack.setCovariances(fwdtrack.getCovariances());

      return extrapolatedMuonTrack;
    }

   public:
    enum MatchingType { MchFirstCluster,
                        MFT_LAST_CLUSTR,
                        EndOfAbsorber,
                        BEGINING_OF_ABSORBER };

    MatchingParamsML(MUON const& muon, MFT const& mft, Collision const& coll, int MType, o2::field::MagneticField* field) : muontrack(muon), mfttrack(mft), collision(coll), mMatchingType(MType), fieldB(field) {}
    void calcMatchingParams()
    {
      auto mftTrackOnMatchingPlane = propagateMFTtoMatchingPlane();
      auto muonTrackOnMatchingPlane = propagateMUONtoMatchingPlane();

      float dphiRaw = mftTrackOnMatchingPlane.getPhi() - muonTrackOnMatchingPlane.getPhi();
      float dphi = TVector2::Phi_mpi_pi(dphiRaw);
      float deta = mftTrackOnMatchingPlane.getEta() - muonTrackOnMatchingPlane.getEta();

      mDX = mftTrackOnMatchingPlane.getX() - muonTrackOnMatchingPlane.getX();
      mDY = mftTrackOnMatchingPlane.getY() - muonTrackOnMatchingPlane.getY();
      mDPt = mftTrackOnMatchingPlane.getPt() - muonTrackOnMatchingPlane.getPt();
      mDPhi = dphi;
      mDEta = deta;
    }

    void calcGlobalMuonParams()
    {
      auto mftTrackAtDca = propagateMFTtoDCA();
      auto muonTrackAtPV = propagateMUONtoPV();

      float momentum = muonTrackAtPV.getP();
      float theta = mftTrackAtDca.getTheta();
      float phiTrack = mftTrackAtDca.getPhi();
      float px = momentum * std::sin(theta) * std::cos(phiTrack);
      float py = momentum * std::sin(theta) * std::sin(phiTrack);

      mGlobalMuonQ = muontrack.sign() + mfttrack.sign();
      mGlobalMuonPtAtDCA = std::sqrt(px * px + py * py);
      mGlobalMuonEtaAtDCA = mftTrackAtDca.getEta();
      mGlobalMuonPhiAtDCA = mftTrackAtDca.getPhi();
      mGlobalMuonDCAx = mftTrackAtDca.getX() - collision.posX();
      mGlobalMuonDCAy = mftTrackAtDca.getY() - collision.posY();
    }

    [[nodiscard]] inline float getDx() const { return mDX; }
    [[nodiscard]] inline float getDy() const { return mDY; }
    [[nodiscard]] inline float getDphi() const { return mDPhi; }
    [[nodiscard]] inline float getDeta() const { return mDEta; }
    [[nodiscard]] inline float getDpt() const { return mDPt; }
    [[nodiscard]] inline float getGMPtAtDCA() const { return mGlobalMuonPtAtDCA; }
    [[nodiscard]] inline float getGMEtaAtDCA() const { return mGlobalMuonEtaAtDCA; }
    [[nodiscard]] inline float getGMPhiAtDCA() const { return mGlobalMuonPhiAtDCA; }
    [[nodiscard]] inline float getGMDcaX() const { return mGlobalMuonDCAx; }
    [[nodiscard]] inline float getGMDcaY() const { return mGlobalMuonDCAy; }
    [[nodiscard]] inline float getGMDcaXY() const { return std::sqrt(mGlobalMuonDCAx * mGlobalMuonDCAx + mGlobalMuonDCAy * mGlobalMuonDCAy); }
    [[nodiscard]] inline int16_t getGMQ() const { return static_cast<int16_t>(mGlobalMuonQ); }

  }; // end of class MatchingParamsML

  template <typename MUON, typename Collisions>
  o2::dataformats::GlobalFwdTrack propagateMUONtoPV(MUON const& muontrack, Collisions const& collisions)
  {
    auto collision = collisions.rawIteratorAt(muontrack.collisionId());
    o2::globaltracking::MatchGlobalFwd mMatching;
    o2::dataformats::GlobalFwdTrack extrapolatedMuonTrack;

    SMatrix5 tpars(muontrack.x(), muontrack.y(), muontrack.phi(), muontrack.tgl(), muontrack.signed1Pt());
    std::vector<float> v1{muontrack.cXX(), muontrack.cXY(), muontrack.cYY(),
                          muontrack.cPhiX(), muontrack.cPhiY(), muontrack.cPhiPhi(),
                          muontrack.cTglX(), muontrack.cTglY(), muontrack.cTglPhi(),
                          muontrack.cTglTgl(), muontrack.c1PtX(), muontrack.c1PtY(),
                          muontrack.c1PtPhi(), muontrack.c1PtTgl(), muontrack.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    double chi2 = muontrack.chi2();
    o2::track::TrackParCovFwd parcovmuontrack{muontrack.z(), tpars, tcovs, chi2};

    o2::dataformats::GlobalFwdTrack gtrack;
    gtrack.setParameters(tpars);
    gtrack.setZ(parcovmuontrack.getZ());
    gtrack.setCovariances(tcovs);
    auto mchtrack = mMatching.FwdtoMCH(gtrack);

    o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

    auto fwdtrack = mMatching.MCHtoFwd(mchtrack);
    extrapolatedMuonTrack.setParameters(fwdtrack.getParameters());
    extrapolatedMuonTrack.setZ(fwdtrack.getZ());
    extrapolatedMuonTrack.setCovariances(fwdtrack.getCovariances());

    return extrapolatedMuonTrack;
  }

  inline bool isGoodTagDimuon(float M)
  {
    return M >= cfgTagMassWindowMin && M <= cfgTagMassWindowMax;
  }

  inline bool isGoodTagMatching(float mDX, float mDY, float mDEta, float mDPhi)
  {
    float dxNorm = (mDX - cfgMeanXTagMuonCut) / (cfgSigmaXTagMuonCut * 3);
    float dyNorm = (mDY - cfgMeanYTagMuonCut) / (cfgSigmaYTagMuonCut * 3);
    float detaNorm = (mDEta - cfgMeanEtaTagMuonCut) / (cfgSigmaEtaTagMuonCut * 3);
    float dphiNorm = (mDPhi - cfgMeanPhiTagMuonCut) / (cfgSigmaPhiTagMuonCut * 3);

    float rTagXY = dxNorm * dxNorm + dyNorm * dyNorm;
    float rTagEtaPhi = detaNorm * detaNorm + dphiNorm * dphiNorm;

    return (rTagXY < 1.f && rTagEtaPhi > 0.f);
  }

  template <typename MUON, typename MFT>
  bool isCorrectMatching(MUON const& muontrack, MFT const& mfttrack)
  {

    int idmuon = muontrack.mcParticleId();
    int idmft = mfttrack.mcParticleId();

    if (idmuon == -1 || idmft == -1) {
      return false;
    }
    return idmuon == idmft;
  };

  template <typename MUON>
  bool isGoodMuonQuality(const MUON& muontrack)
  {
    if (!muontrack.has_collision()) {
      return false;
    }
    if (muontrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      return false;
    }
    if (muontrack.chi2() > legacySelection.cfgTrackChi2MchUp) {
      return false;
    }
    if (legacySelection.cfgRabsLow1 > muontrack.rAtAbsorberEnd() || muontrack.rAtAbsorberEnd() > legacySelection.cfgRabsUp2) {
      return false;
    }
    if (muontrack.rAtAbsorberEnd() < legacySelection.cfgRabsUp1 && legacySelection.cfgPdcaUp1 < muontrack.pDca()) {
      return false;
    }
    if (muontrack.rAtAbsorberEnd() > legacySelection.cfgRabsLow2 && legacySelection.cfgPdcaUp2 < muontrack.pDca()) {
      return false;
    }
    return true;
  }

  template <typename MUON>
  bool isGoodMuonKine(const MUON& muontrack)
  {
    return !static_cast<bool>(legacySelection.cfgEtaMchLow > muontrack.getEta() || muontrack.getEta() > legacySelection.cfgEtaMchUp);
  }

  template <typename MFT>
  bool isGoodMFTQuality(const MFT& mfttrack)
  {
    if (!mfttrack.has_collision()) {
      return false;
    }
    if (mfttrack.chi2() > cfgTrackChi2MftUp) {
      return false;
    }
    if (mfttrack.nClusters() < cfgTrackNClustMftLow) {
      return false;
    }
    return true;
  }

  template <typename MFT>
  bool isGoodMFTKine(const MFT& mfttrack)
  {
    return !static_cast<bool>(cfgEtaMftlow > mfttrack.getEta() || mfttrack.getEta() > cfgEtaMftup);
  }

  inline bool isPassMatchingPreselection(float Dx, float Dy)
  {
    return std::abs(Dx) <= cfgPreselectMatchingX && std::abs(Dy) <= cfgPreselectMatchingY;
  }

  template <typename MUONs, typename Collisions>
  void setMUONs(MUONs const& muontracks, Collisions const& collisions)
  {
    for (const auto& muontrack : muontracks) {
      if (!isGoodMuonQuality(muontrack)) {
        continue;
      }
      o2::dataformats::GlobalFwdTrack muonTrackAtPV = propagateMUONtoPV(muontrack, collisions);
      if (!isGoodMuonKine(muonTrackAtPV)) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(muontrack.collisionId());

      bool& has = legacy.hasMuonTracksByCollision[muontrack.collisionId()];
      has = true;

      std::vector<int64_t>& muonTrackIds = legacy.muonTracksByCollision[collision.globalIndex()];
      muonTrackIds.push_back(muontrack.globalIndex());
    }
  }

  template <typename MFT, typename Collisions>
  o2::track::TrackParCovFwd propagateMFTtoDCA(MFT const& mfttrack, Collisions const& collisions, o2::field::MagneticField* field)
  {
    auto collision = collisions.rawIteratorAt(mfttrack.collisionId());
    std::vector<double> mftv1;
    SMatrix55 mftcovs{mftv1.begin(), mftv1.end()};
    SMatrix5 mftpars = {mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt()};
    o2::track::TrackParCovFwd mftpartrack = {mfttrack.z(), mftpars, mftcovs, mfttrack.chi2()};
    std::array<double, 3> propVec = {std::fabs(mfttrack.x() - collision.posX()),
                                     std::fabs(mfttrack.y() - collision.posY()),
                                     std::fabs(mfttrack.z() - collision.posZ())};
    std::array<double, 3> centerZ = {mfttrack.x() - propVec[0] / 2.,
                                     mfttrack.y() - propVec[1] / 2.,
                                     mfttrack.z() - propVec[2] / 2.};
    float bz = field->getBz(centerZ.data());
    mftpartrack.propagateToZ(collision.posZ(), bz);
    return mftpartrack;
  }

  template <typename MFTs, typename Collisions>
  void setMFTs(MFTs const& mfttracks, Collisions const& collisions, o2::field::MagneticField* field)
  {
    for (const auto& mfttrack : mfttracks) {
      if (!isGoodMFTQuality(mfttrack)) {
        continue;
      }

      o2::track::TrackParCovFwd mftTrackAtDca = propagateMFTtoDCA(mfttrack, collisions, field);
      if (!isGoodMFTKine(mftTrackAtDca)) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(mfttrack.collisionId());

      legacy.vertexZByCollision[mfttrack.collisionId()] = collision.posZ();
      legacy.mftMultiplicityByCollision[mfttrack.collisionId()] += 1;

      bool& has = legacy.hasMftTracksByCollision[mfttrack.collisionId()];
      has = true;

      std::vector<int64_t>& mftTrackIds = legacy.mftTracksByCollision[collision.globalIndex()];
      mftTrackIds.push_back(mfttrack.globalIndex());
    }
  }

  o2::framework::Produces<o2::aod::MatchParams> tableMatchingParams;
  o2::framework::Produces<o2::aod::TagMatchParams> tableTagMatchingParams;
  o2::framework::Produces<o2::aod::ProbeMatchParams> tableProbeMatchingParams;
  o2::framework::Produces<o2::aod::MixMatchParams> tableMixMatchingParams;
  o2::framework::Produces<o2::aod::MuonPair> tableMuonPair;

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdbManager{};

  o2::field::MagneticField* fieldB = nullptr;
  o2::ccdb::CcdbApi ccdbApi;

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (legacy.mRunNumber == bc.runNumber()) {
      return;
    }

    legacy.mRunNumber = bc.runNumber();
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, legacy.mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(geoPath);
    }
    if (matchingPlaneProcessCount() > 0) {
      // Match both geometry loads and detector-element navigation performed by
      // global-muon-matcher, even with its realignment disabled.
      const std::array<int, 10> elementsPerChamber{4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
      for (int geometryCopy = 0; geometryCopy < GeometryWarmupPasses; ++geometryCopy) {
        auto* geometry = ccdbManager->getForTimeStamp<TGeoManager>(geoPath, bc.timestamp());
        if (!geometry) {
          LOGF(fatal, "Aligned geometry is unavailable at the BC timestamp");
        }
        ccdbManager->clearCache(geoPath);
        auto transformation = o2::mch::geo::transformationFromTGeoManager(*geometry);
        for (int chamber = 1; chamber <= MchChamberCount; ++chamber) {
          for (int element = 0; element < elementsPerChamber[chamber - 1]; ++element) {
            (void)transformation(100 * chamber + element);
          }
        }
      }
    }
    o2::mch::TrackExtrap::setField();
    fieldB = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (!fieldB) {
      LOGF(fatal, "Expected an O2 magnetic field after CCDB initialization");
    }
  }

  void run(o2::framework::ProcessingContext& context)
  {
    mixingPersistence.state.frame = context.services().get<o2::framework::TimingInfo>().timeslice;
  }
  void init(o2::framework::InitContext&)
  {
    if (mchAcceptance.cfgRequireMchPVEta &&
        (!std::isfinite(mchAcceptance.cfgMchPVEtaMin.value) || !std::isfinite(mchAcceptance.cfgMchPVEtaMax.value) ||
         !(mchAcceptance.cfgMchPVEtaMin < mchAcceptance.cfgMchPVEtaMax))) {
      LOGF(fatal, "Invalid MCH PV eta selection range");
    }
    if (mchAcceptance.cfgRequireMchMftAcceptance) {
      if (!std::isfinite(mchAcceptance.cfgMchMftAcceptanceEtaMin.value) || !std::isfinite(mchAcceptance.cfgMchMftAcceptanceEtaMax.value) ||
          mchAcceptance.cfgMchMftAcceptanceEtaMin >= mchAcceptance.cfgMchMftAcceptanceEtaMax || mchAcceptance.cfgMchMftAcceptanceEtaMax >= 0. || mchAcceptance.cfgMchMftAcceptanceDiskZ.value.empty()) {
        LOGF(fatal, "Invalid nominal MFT acceptance envelope");
      }
      for (double const& z : mchAcceptance.cfgMchMftAcceptanceDiskZ.value) {
        if (!std::isfinite(z) || z >= 0.) {
          LOGF(fatal, "MFT acceptance disk z must be finite and negative");
        }
      }
      if (mchAcceptance.cfgMchMftAcceptanceMinDisks < 1 ||
          static_cast<size_t>(mchAcceptance.cfgMchMftAcceptanceMinDisks.value) > mchAcceptance.cfgMchMftAcceptanceDiskZ.value.size() ||
          std::set<double>(mchAcceptance.cfgMchMftAcceptanceDiskZ.value.begin(), mchAcceptance.cfgMchMftAcceptanceDiskZ.value.end()).size() != mchAcceptance.cfgMchMftAcceptanceDiskZ.value.size()) {
        LOGF(fatal, "Nominal MFT acceptance requires distinct disk z positions and a minimum between 1 and their count");
      }
    }
    mchAcceptanceQA.add("status", "MCH acceptance: 0=tested, 1=PV failure, 2=envelope failure, 3=pass;status;MCH tracks", o2::framework::HistType::kTH1F, {{4, -0.5, 3.5}});
    mchAcceptanceQA.add("pvEtaStatus", "MCH PV eta: 0=tested, 1=PV failure, 2=eta failure, 3=pass;status;MCH tracks", o2::framework::HistType::kTH1F, {{4, -0.5, 3.5}});
    LOGF(info, "MCH acceptance: disk envelope=%d, PV eta cut=%d, strict eta range=(%g,%g)",
         static_cast<int>(mchAcceptance.cfgRequireMchMftAcceptance.value), static_cast<int>(mchAcceptance.cfgRequireMchPVEta.value),
         mchAcceptance.cfgMchPVEtaMin.value, mchAcceptance.cfgMchPVEtaMax.value);
    mixingPersistence.state.stream = mixingPersistence.cfgMixingStreamId.value;
    if (mixingPersistence.state.stream == 0) {
      std::random_device random;
      do {
        mixingPersistence.state.stream = (static_cast<uint64_t>(random()) << 32) ^ static_cast<uint64_t>(random());
      } while (mixingPersistence.state.stream == 0);
    }
    LOGF(info, "Mixing source stream ID: %" PRIu64 "", static_cast<uint64_t>(mixingPersistence.state.stream));
    if (mixingPersistence.cfgMixingPoolMaxEvents <= 0 || mixingPersistence.cfgMixingPoolMaxMftStates <= 0) {
      LOGF(fatal, "Cross-DF pool limits must be positive");
    }
    if (hasInternalAssociationProcess() &&
        (!std::isfinite(localAssociation.cfgMftAssocNumSigma.value) || localAssociation.cfgMftAssocNumSigma < 0.f ||
         !std::isfinite(localAssociation.cfgMftAssocTimeMargin.value) || localAssociation.cfgMftAssocTimeMargin < 0.f || localAssociation.cfgMftAssocBcWindow < 0)) {
      LOGF(fatal, "Internal MFT association requires finite nonnegative time-compatibility parameters");
    }
    if (hasNoCovarianceProcess() && (cfgRequireMftCovariance || cfgComputeMatchingChi2)) {
      LOGF(fatal, "Association NoCov input requires cfgRequireMftCovariance=false and cfgComputeMatchingChi2=false; chi2 evaluation may use only saved official pairs");
    }
    if (!cfgRequireMftCovariance && cfgComputeMatchingChi2 && !inferenceConfig.cfgEvaluationMode) {
      LOGF(fatal, "All-MFT training without inference requires cfgComputeMatchingChi2=false");
    }
    if (inferenceConfig.cfgEvaluationMode) {
      if (doprocessLegacy || (matchingPlaneProcessCount() == 0)) {
        LOGF(fatal, "Inference requires matching-plane process, no legacy process");
      }
      inference.requestedChi2Method = inferenceConfig.cfgChi2Method.value;
      inference.useStoredMatchingChi2 = !cfgComputeMatchingChi2.value;
      if (inferenceConfig.cfgEnableChi2Matching && inference.useStoredMatchingChi2) {
        if (MatchingInference::isCombinedChi2(inferenceConfig.cfgChi2Method.value)) {
          LOGF(fatal, "cfgChi2Method=%s requires cfgComputeMatchingChi2=true and measured MFT covariance; the combined score cannot be obtained from the AO2D saved matching chi2", inference.requestedChi2Method.c_str());
        }
        if (inferenceConfig.cfgChi2Method.value != "77_all") {
          LOGF(warning, "cfgComputeMatchingChi2=false: requested %s cannot be obtained from AO2D. Switching cfgChi2Method to 77_all using only the official saved chi2MatchMCHMFT scalar; no 4/2-parameter decomposition or -505 cm recomputation. Requested and effective methods are recorded in metadata", inference.requestedChi2Method.c_str());
          inferenceConfig.cfgChi2Method.value = "77_all";
        }
        LOGF(info, "Chi2 evaluation uses AO2D official GlobalMuon/GlobalForward pair IDs and chi2MatchMCHMFT; only pairs passing this task's candidate preselection enter the one-to-one assignment");
      }
      if (!inferenceConfig.cfgTrainingMode && mixingCuts.cfgEnableEventMixing) {
        LOGF(fatal, "Inference-only: disable MFT-MCH training event mixing; use cfgEnableDimuonMixing");
      }
      if ((inferenceConfig.cfgDimuonMixPvMultiplicityRel >= 0. || inferenceConfig.cfgDimuonMixDeltaCentrality >= 0.) && !hasCentralityProcess()) {
        LOGF(fatal, "Requested dimuon centrality/PV-mult pool requires centrality process");
      }
    }
    inference.init(inferenceConfig, inferenceHistos);
    if (inferenceConfig.cfgEvaluationMode) {
      if ((inference.mlChi2Upstream || inference.mlChi2Downstream) && !cfgComputeMatchingChi2) {
        LOGF(fatal, "Selected ML model requires chi2 input features; enable cfgComputeMatchingChi2");
      }
      if (inferenceConfig.cfgEnableMlMatching && (!cfgMftUseMchMomentum || cfgMftUseMchCharge.value != inference.requiredMchCharge || cfgMftAbsorberStepCm.value != inference.requiredAbsorberStep)) {
        LOGF(fatal, "Propagation settings disagree with trained ML feature contract");
      }
      std::ostringstream meta;
      meta << "ML=" << inferenceConfig.cfgMlModel.value << ";contract=" << inference.contract
           << ";MLcut=" << inferenceConfig.cfgMlScoreMin.value << ";chi2Requested=" << inference.requestedChi2Method << ";chi2=" << inferenceConfig.cfgChi2Method.value
           << ";chi2Source=" << (inference.useStoredMatchingChi2 ? "AO2D.chi2MatchMCHMFT (official production-defined scalar)" : "recomputed candidate chi2 selected by cfgChi2Method")
           << ";matchingChi2Max=" << MatchingInference::matchingChi2Limit(inferenceConfig) << ";MIDrequired=" << inferenceConfig.cfgRequireMid.value
           << ";kinematics=MCH total P at PV + MFT position/direction, MCH charge, then official 3D helix DCA"
           << ";globalDimuonAxes=mass,pt,rapidity,pairDCA,PCAxy,PCAz,DCAxy,Lxy,fitSuccess;vertexCutsAppliedToGlobalSparse=false;dimuonDCAxy=hypot(singleDCAxy1,singleDCAxy2) cm, no sqrt(2) divisor"
           << ";mixDepth=" << inferenceConfig.cfgDimuonMixingDepth.value << ";mixDeltaZ=" << inferenceConfig.cfgDimuonMixDeltaZ.value
           << ";mixNormalize=" << inferenceConfig.cfgDimuonMixNormalize.value << ";source score is NOT correct-match probability\n"
           << inference.contractContent;
      meta << "\nSingle DCA limits (cm), applied after assignment and score cuts, no reassignment: xy=" << inferenceConfig.cfgSingleDcaXYMax.value
           << ";absZ=" << inferenceConfig.cfgSingleDcaZMax.value << ";applyToMchPV=" << inferenceConfig.cfgApplySingleDcaCutsToMchPV.value;
      const auto appendValues = [&meta](const char* name, auto const& values) {
        meta << ';' << name << "=[";
        for (size_t i = 0; i < values.size(); ++i) {
          if (i != 0u) {
            meta << ',';
          }
          meta << values[i];
        }
        meta << ']';
      };
      appendValues("mchPtPVEdges", inferenceConfig.cfgSingleDcaPtEdges.value);
      appendValues("xyMaxCm", inferenceConfig.cfgSingleDcaXYMaxPtCuts.value);
      appendValues("absZMaxCm", inferenceConfig.cfgSingleDcaZMaxPtCuts.value);
      inferenceMetadata->SetTitle(meta.str().c_str());
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel.value) || mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel < -1.) {
      LOGF(fatal, "Invalid mixing Central barrel multiplicity limit");
    }
    if (mixingCuts.cfgEnableEventMixing && mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel >= 0. && !hasCentralityProcess()) {
      LOGF(fatal, "Central barrel mixing cut requires processMatchingPlaneCentrality=true and the official PVMults table");
    }
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    ccdbManager->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    legacy.mRunNumber = 0;
    if (mixingCuts.cfgEnableEventMixing && matchingPlaneProcessCount() == 0) {
      LOGF(fatal, "Event mixing requires processMatchingPlane");
    }
    if (matchingPlaneProcessCount() > 1) {
      LOGF(fatal, "Enable only one matching-plane process");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMftMultiplicityRel.value) || mixingMultiplicityCuts.cfgMixingMaxDeltaMftMultiplicityRel < -1.f ||
        mixingMultiplicityCuts.cfgMixingMaxDeltaMchMultiplicity < -1 || !std::isfinite(mixingCuts.cfgMixingMaxDeltaCentrality.value) || mixingCuts.cfgMixingMaxDeltaCentrality < -1.f) {
      LOGF(fatal, "Invalid mixing multiplicity/centrality limits");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMchPtRel.value)) {
      LOGF(fatal, "Invalid mixing MCH pT limit");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMchMomentumRel.value)) {
      LOGF(fatal, "Invalid mixing MCH momentum limit");
    }
    if (!std::isfinite(mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV.value) || !std::isfinite(mixingCuts.cfgMixingMaxDeltaMchEtaPV.value)) {
      LOGF(fatal, "Invalid mixing MCH direction limits");
    }
    if (mixingCuts.cfgEnableEventMixing && mixingCuts.cfgMixingMaxDeltaCentrality >= 0.f && !hasCentralityProcess()) {
      LOGF(fatal, "Centrality mixing cut requires processMatchingPlaneCentrality=true and processMatchingPlane=false");
    }
    if (mixingCuts.cfgMixingDepth < -1 || mixingCuts.cfgMinDeltaBc < 0 || !std::isfinite(mixingCuts.cfgMaxDeltaPhi.value) ||
        !std::isfinite(mixingCuts.cfgMaxDeltaR.value) || !std::isfinite(mixingCuts.cfgMaxDeltaAttempts.value) ||
        !std::isfinite(mixingCuts.cfgMaxDeltaZ.value)) {
      LOGF(fatal, "Invalid event-mixing configuration");
    }
    if (matchingPlaneProcessCount() > 0) {
      if (doprocessLegacy) {
        LOGF(fatal, "Enable only one of processLegacy and processMatchingPlane");
      }
      const auto& cutName = mchAcceptance.cfgMuonCuts.value;
      if (cutName != "none" && cutName != "muonQualityCuts" && cutName != "muonQualityCuts5SigmaPDCA_Run3" && cutName != "muonQualityCuts10SigmaPDCA") {
        LOGF(fatal, "Unsupported mchAcceptance.cfgMuonCuts: %s", cutName.c_str());
      }
      if (cutName != "none") {
        mchQualityCut.reset(o2::aod::dqcuts::GetAnalysisCut(cutName.c_str()));
        if (!mchQualityCut) {
          LOGF(fatal, "Cannot load PWGDQ muon cut");
        }
      }
      if (cfgEventSelection.value == "all") {
        matchingEventSelectionMode = 0;
      } else if (cfgEventSelection.value == "int8" || cfgEventSelection.value == "sel8") {
        matchingEventSelectionMode = 1;
      } else if (cfgEventSelection.value == "upc") {
        matchingEventSelectionMode = 2;
      } else {
        LOGF(fatal, "cfgEventSelection must be all, int8/sel8 or upc");
      }
      if (cfgUpcMinMchTracks < 0 || cfgUpcMinMftTracks < 0 || (cfgUpcMaxMchTracks != -1 && cfgUpcMaxMchTracks < cfgUpcMinMchTracks) || cfgUpcMinMchMid < 0 || (cfgUpcMaxMchMid != -1 && cfgUpcMaxMchMid < cfgUpcMinMchMid) ||
          cfgUpcFV0BcWindow < 0 || !std::isfinite(cfgUpcMaxFT0AAmplitude.value) || !std::isfinite(cfgUpcMaxFT0CAmplitude.value) || !std::isfinite(cfgUpcMaxFDDAAmplitude.value) || !std::isfinite(cfgUpcMaxFDDCAmplitude.value) || !std::isfinite(cfgUpcMaxFV0Amplitude.value) ||
          !std::isfinite(cfgUpcFV0ValidTime.value) || cfgUpcFV0ValidTime < 0.f) {
        LOGF(fatal, "Invalid UPC multiplicity, FV0 threshold, BC window or valid-time configuration");
      }
      LOGF(info, "Matching-plane event selection: %s (mode %d); only MCH output is selected, MFT candidates retain the original time window", cfgEventSelection.value.c_str(), matchingEventSelectionMode);
      if (!std::isfinite(cfgMftAbsorberStepCm.value) || cfgMftAbsorberStepCm < 0.) {
        LOGF(fatal, "cfgMftAbsorberStepCm must be finite and nonnegative");
      }
      LOGF(info, "Always save both planes: -77.5 cm official transport; -505 cm material step=%g cm (0=whole absorber); magnetic integration remains official/adaptive", cfgMftAbsorberStepCm.value);
      if (!std::isfinite(cfgPairJacobianScale.value) || cfgPairJacobianScale <= 0.) {
        LOGF(fatal, "cfgPairJacobianScale must be finite and positive");
      }
      o2::conf::ConfigurableParam::setValue("FieldOriginBias.z", "0");
    }
  }

  void processLegacy(MyCollisions const& collisions,
                     MyBCs const& bcs,
                     MyMUONs const& muontracks,
                     MyMFTs const& mfttracks)
  {
    LOG(info) << "Process()  ";
    legacy.muonTracksByCollision.clear();
    legacy.mftTracksByCollision.clear();
    legacy.matchingCollisions.clear();
    legacy.hasMuonTracksByCollision.clear();
    legacy.hasMftTracksByCollision.clear();

    initCCDB(bcs.begin());
    setMUONs(muontracks, collisions);
    setMFTs(mfttracks, collisions, fieldB);

    for (auto const& mapHasMuonTracksCollision : legacy.hasMuonTracksByCollision) {
      auto muonCollisionId = mapHasMuonTracksCollision.first;
      for (auto const& mapHasMftTracksCollision : legacy.hasMftTracksByCollision) {
        auto mftCollisionId = mapHasMftTracksCollision.first;
        if (muonCollisionId != mftCollisionId) {
          continue;
        }
        legacy.matchingCollisions[mftCollisionId] = true;
      }
    }

    for (auto const& map_collision : legacy.matchingCollisions) {
      auto const& collision = collisions.rawIteratorAt(map_collision.first);

      for (auto const& imuontrack1 : legacy.muonTracksByCollision[map_collision.first]) {
        auto const& muontrack1 = muontracks.rawIteratorAt(imuontrack1);

        for (auto const& imfttrack1 : legacy.mftTracksByCollision[map_collision.first]) {
          auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);

          MatchingParamsML<MyMUON, MyMFT, MyCollision> matching(muontrack1, mfttrack1, collision, cfgMatchingMethod, fieldB);
          matching.calcMatchingParams();

          if (!isPassMatchingPreselection(matching.getDx(), matching.getDy())) {
            continue;
          }

          matching.calcGlobalMuonParams();

          bool isTrue = false;

          tableMatchingParams(matching.getGMPtAtDCA(),
                              matching.getGMEtaAtDCA(),
                              matching.getGMQ(),
                              matching.getDpt(),
                              matching.getDx(),
                              matching.getDy(),
                              matching.getDeta(),
                              matching.getDphi(),
                              isTrue);
        }

        for (auto const& map_mfttrack : legacy.mftTracksByCollision) {
          if (map_mfttrack.first == map_collision.first) {
            continue;
          }
          if (std::fabs(legacy.vertexZByCollision[map_mfttrack.first] - legacy.vertexZByCollision[map_collision.first]) > cfgEventMaxDeltaVtxZ) {
            continue;
          }
          if (std::fabs(legacy.mftMultiplicityByCollision[map_mfttrack.first] - legacy.mftMultiplicityByCollision[map_collision.first]) > cfgEventMaxDeltaNMFT) {
            continue;
          }

          for (auto const& imfttrack1 : map_mfttrack.second) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matching(muontrack1, mfttrack1, collision, cfgMatchingMethod, fieldB);
            matching.calcMatchingParams();
            if (!isPassMatchingPreselection(matching.getDx(), matching.getDy())) {
              continue;
            }
            matching.calcGlobalMuonParams();

            bool isTrue = false;

            tableMixMatchingParams(matching.getGMPtAtDCA(),
                                   matching.getGMEtaAtDCA(),
                                   matching.getGMQ(),
                                   matching.getDpt(),
                                   matching.getDx(),
                                   matching.getDy(),
                                   matching.getDeta(),
                                   matching.getDphi(),
                                   isTrue);
          }
        }

        for (auto const& imuontrack2 : legacy.muonTracksByCollision[map_collision.first]) {

          if (imuontrack1 >= imuontrack2) {
            continue;
          }

          auto const& muontrack2 = muontracks.rawIteratorAt(imuontrack2);

          FindTagAndProbe<MyMUON, MyCollision> tagdimuon(muontrack1, muontrack2, collision);
          tagdimuon.calcMuonPairAtPV();
          tableMuonPair(tagdimuon.getCharge(), tagdimuon.getMass(), tagdimuon.getPt(), tagdimuon.getRap());

          if (!isGoodTagDimuon(tagdimuon.getMass())) {
            continue;
          }

          auto tagmuontrack = muontrack1;
          auto probemuontrack = muontrack2;

          if (tagdimuon.getTagMuonIndex() == 1) {
            tagmuontrack = muontrack2;
            probemuontrack = muontrack1;
          }

          int nTagMFTCand = 0;
          int nProbeMFTCand = 0;

          int indexTagMftCandidate = -1;
          float tagGMPtAtDCA = 0;
          // float tagGMEtaAtDCA = 0;

          float minimumR = 9999.;
          int minimumIndexProbeMFTCand = -1;

          std::unordered_map<int, std::vector<float>> probeMatchingParameters;

          for (auto const& imfttrack1 : legacy.mftTracksByCollision[map_collision.first]) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matchingTag(tagmuontrack, mfttrack1, collision, cfgMatchingMethod, fieldB);
            matchingTag.calcMatchingParams();
            matchingTag.calcGlobalMuonParams();
            if (isGoodTagMatching(matchingTag.getDx(), matchingTag.getDy(), matchingTag.getDeta(), matchingTag.getDphi()) &&
                isPassMatchingPreselection(matchingTag.getDx(), matchingTag.getDy())) {
              bool isTrue = false;
              tableTagMatchingParams(matchingTag.getGMPtAtDCA(),
                                     matchingTag.getGMEtaAtDCA(),
                                     matchingTag.getGMQ(),
                                     matchingTag.getDpt(),
                                     matchingTag.getDx(),
                                     matchingTag.getDy(),
                                     matchingTag.getDeta(),
                                     matchingTag.getDphi(),
                                     isTrue);
              indexTagMftCandidate = mfttrack1.globalIndex();
              tagGMPtAtDCA = matchingTag.getGMPtAtDCA();
              // tagGMEtaAtDCA = matchingTag.getGMEtaAtDCA();
              ++nTagMFTCand;
            }
          } // end of loop imfttrack1

          if (nTagMFTCand != 1) {
            continue;
          }
          for (auto const& imfttrack1 : legacy.mftTracksByCollision[map_collision.first]) {
            auto const& mfttrack1 = mfttracks.rawIteratorAt(imfttrack1);
            if (mfttrack1.globalIndex() == indexTagMftCandidate) {
              continue;
            }
            MatchingParamsML<MyMUON, MyMFT, MyCollision> matchingProbe(probemuontrack, mfttrack1, collision, cfgMatchingMethod, fieldB);
            matchingProbe.calcMatchingParams();
            if (isPassMatchingPreselection(matchingProbe.getDx(), matchingProbe.getDy())) {
              float radius = std::sqrt(matchingProbe.getDx() * matchingProbe.getDx() + matchingProbe.getDy() * matchingProbe.getDy());
              bool isTrue = false;
              matchingProbe.calcGlobalMuonParams();
              std::vector<float>& probeMatchingParams = probeMatchingParameters[nProbeMFTCand];
              probeMatchingParams.push_back(tagGMPtAtDCA);
              probeMatchingParams.push_back(matchingProbe.getGMPtAtDCA());
              probeMatchingParams.push_back(matchingProbe.getGMEtaAtDCA());
              probeMatchingParams.push_back(matchingProbe.getGMQ());
              probeMatchingParams.push_back(matchingProbe.getDpt());
              probeMatchingParams.push_back(matchingProbe.getDx());
              probeMatchingParams.push_back(matchingProbe.getDy());
              probeMatchingParams.push_back(matchingProbe.getDeta());
              probeMatchingParams.push_back(matchingProbe.getDphi());
              probeMatchingParams.push_back(static_cast<float>(isTrue));
              if (radius < minimumR) {
                minimumIndexProbeMFTCand = nProbeMFTCand;
                minimumR = radius;
              }
              ++nProbeMFTCand;
            }
          } // end of loop imfttrack1

          if (nProbeMFTCand < 1) {
            continue;
          }

          if (minimumIndexProbeMFTCand > -1) {
            std::vector<float>& probeMatchingParams = probeMatchingParameters[minimumIndexProbeMFTCand];
            tableProbeMatchingParams(probeMatchingParams[0],
                                     probeMatchingParams[1],
                                     probeMatchingParams[2],
                                     static_cast<int16_t>(probeMatchingParams[3]),
                                     probeMatchingParams[4],
                                     probeMatchingParams[5],
                                     probeMatchingParams[6],
                                     probeMatchingParams[7],
                                     probeMatchingParams[8],
                                     static_cast<bool>(probeMatchingParams[9]));
          }

        } // end of loop imuontrack2
      } // end of loop imuontrack1
    } // end of loop map_collision

  } // end of processMC

  void processMatchingPlaneImpl(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                aod::BCsWithTimestamps const& bcs,
                                MyMUONs const& muontracks,
                                aod::MFTTracks const& mfttracks,
                                aod::MFTTracksCov const& mftCovs,
                                aod::FV0As const& fv0s,
                                aod::FT0s const& ft0s,
                                aod::FDDs const& fdds, aod::Origins const& origins,
                                std::vector<std::pair<int64_t, int64_t>> const& associations = {},
                                std::unordered_map<int64_t, int64_t> const& ambiguousBCs = {}, bool useAssociation = false)
  {
    using o2::aod::fwdtrackutils::getTrackParCovFwd;
    if (!inferenceConfig.cfgTrainingMode && !inferenceConfig.cfgEvaluationMode) {
      return;
    }
    auto& history = mixingPersistence.state.history;
    const uint64_t frame = mixingPersistence.state.frame;
    bool poolReset = false;
    // Publish even on early return, so MC caches see every source frame.
    struct PublishOnExit {
      std::function<void()> f;
      ~PublishOnExit() noexcept
      {
        try {
          f();
        } catch (...) {
          std::terminate();
        }
      }
    } publishPool{[&]() {
      std::set<uint64_t> frames;
      for (const auto& event : history.events) {
        frames.insert(event.frame);
      }
      tableMixingPoolStates(frame, std::vector<uint64_t>(frames.begin(), frames.end()),
                            history.events.size(), history.states, history.evicted, history.oversized, poolReset, mixingPersistence.state.stream);
    }};
    if (!mixingPersistence.cfgMixingAcrossDF) {
      history.clear();
    }
    if (bcs.size() > 0) {
      poolReset = history.setRun(bcs.begin().runNumber());
    }
    if (bcs.size() == 0) {
      return;
    }
    if (collisions.size() == 0) {
      // There is no reconstructed vertex to match to. Retain the MCH rows as
      // ineligible, without following potentially stale collision references.
      if (cfgSaveMatchingAudit) {
        const uint64_t sourceFrame = (origins.size() != 0) ? origins.begin().dataframeID() : frame;
        for (auto const& mch : muontracks) {
          if (mch.trackType() > o2::aod::fwdtrack::GlobalForwardTrack) {
            mchPopulation(sourceFrame, mch.collisionId(), mch.globalIndex(), static_cast<int8_t>(mch.signed1Pt() > 0. ? 1 : -1),
                          false, false, NAN, NAN, NAN, NAN, 0, 0, 0, useAssociation, 0, mch.trackType() == o2::aod::fwdtrack::MuonStandaloneTrack);
          }
        }
      }
      return;
    }
    const bool filterBeforeMatching = cfgFilterBeforeMatching && (matchingEventSelectionMode != 0 || mchQualityCut != nullptr || mchAcceptance.cfgRequireMchMftAcceptance || mchAcceptance.cfgRequireMchPVEta);

    globalMchIterationOrder.clear();
    for (const auto& mch : muontracks) {
      if (mch.trackType() > o2::aod::fwdtrack::GlobalForwardTrack) {
        globalMchIterationOrder.try_emplace(mch.globalIndex(), 0);
      }
    }

    // Event selection is anchored on the MCH side. UPC multiplicity counts use
    // the original collision ID and all associated MFT tracks before covariance
    // selection. Candidate matching retains the original timing definition.
    // Event/pool multiplicities retain their native collision definition.
    std::unordered_map<int64_t, int> numMchPerCollision;
    std::unordered_map<int64_t, int> numMftPerCollision;
    std::unordered_map<int64_t, int> numPositiveMchPerCollision;
    std::unordered_map<int64_t, int> numNegativeMchPerCollision;
    for (const auto& mch : muontracks) {
      // Global tracks (types 0--2) duplicate standalone MCH tracks: do not count.
      if (mch.trackType() > o2::aod::fwdtrack::GlobalForwardTrack && mch.has_collision()) {
        ++numMchPerCollision[mch.collisionId()];
        if (std::isfinite(mch.signed1Pt())) {
          if (mch.signed1Pt() > 0.f) {
            ++numPositiveMchPerCollision[mch.collisionId()];
          } else if (mch.signed1Pt() < 0.f) {
            ++numNegativeMchPerCollision[mch.collisionId()];
          }
        }
      }
    }
    for (const auto& mft : mfttracks) {
      if (mft.has_collision()) {
        ++numMftPerCollision[mft.collisionId()];
      }
    }
    std::map<int64_t, std::set<int64_t>> candidateMfts, associatedMfts;
    std::unordered_map<int64_t, int64_t> referenceBCs = ambiguousBCs;
    for (auto const& mft : mfttracks) {
      if (mft.has_collision()) {
        auto col = collisions.rawIteratorAt(mft.collisionId() - collisions.offset());
        referenceBCs[mft.globalIndex()] = bcs.rawIteratorAt(col.bcId() - bcs.offset()).globalBC();
        if (!useAssociation || cfgIncludeNativeMftCandidates) {
          candidateMfts[mft.collisionId()].insert(mft.globalIndex());
        }
      }
    }
    for (auto const& [collisionId, mftId] : associations) {
      if (collisionId < static_cast<int64_t>(collisions.offset()) || collisionId - static_cast<int64_t>(collisions.offset()) >= collisions.size() ||
          mftId < static_cast<int64_t>(mfttracks.offset()) || mftId - static_cast<int64_t>(mfttracks.offset()) >= mfttracks.size()) {
        LOGF(fatal, "MFTTrackAssoc contains an out-of-frame collision or MFT ID");
      }
      candidateMfts[collisionId].insert(mftId);
      associatedMfts[collisionId].insert(mftId);
    }
    auto timeCompatible = [&](auto const& mch, auto const& mft) {
      if (!mch.has_collision()) {
        return false;
      }
      auto ref = referenceBCs.find(mft.globalIndex());
      if (ref == referenceBCs.end()) {
        return false;
      }
      auto col = collisions.rawIteratorAt(mch.collisionId() - collisions.offset());
      const auto bc = bcs.rawIteratorAt(col.bcId() - bcs.offset()).globalBC();
      const double dt = mft.trackTime() - mch.trackTime() +
                        (ref->second - static_cast<int64_t>(bc)) * o2::constants::lhc::LHCBunchSpacingNS;
      const double resolution = mft.trackTimeRes() + mch.trackTimeRes();
      return std::isfinite(dt) && std::isfinite(resolution) && resolution >= 0. && std::abs(dt) <= resolution;
    };
    std::unordered_map<int64_t, int64_t> mchGlobalBCs;
    std::map<int64_t, int> numMchMidPerBC;
    for (const auto& mch : muontracks) {
      if (mch.trackType() <= o2::aod::fwdtrack::GlobalForwardTrack || !mch.has_collision() || !std::isfinite(mch.trackTime())) {
        continue;
      }
      auto collision = collisions.rawIteratorAt(mch.collisionId());
      auto bc = bcs.rawIteratorAt(collision.bcId());
      int64_t correctedBC = static_cast<int64_t>(bc.globalBC()) + TMath::FloorNint(mch.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + 1.);
      mchGlobalBCs.emplace(mch.globalIndex(), correctedBC);
      if (mch.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        ++numMchMidPerBC[correctedBC];
      }
    }
    // Apply the official helper to every FIT record; duplicate BCs must all pass.
    // Keep amplitudes separately for the existing FV0 diagnostic columns.
    std::map<int64_t, float> fv0Amplitudes;
    std::map<int64_t, bool> cleanFITPerBC;
    const std::vector<float> fitLimits{cfgUpcMaxFV0Amplitude, cfgUpcMaxFT0AAmplitude,
                                       cfgUpcMaxFT0CAmplitude, cfgUpcMaxFDDAAmplitude, cfgUpcMaxFDDCAmplitude};
    auto recordFIT = [&](int64_t bcId, MatchingFITRecord record) {
      int64_t gbc = bcs.rawIteratorAt(bcId).globalBC();
      bool clean = udhelpers::cleanFIT(record, cfgUpcFV0ValidTime, fitLimits);
      auto [it, inserted] = cleanFITPerBC.emplace(gbc, clean);
      if (!inserted) {
        it->second = it->second && clean;
      }
    };
    if (matchingEventSelectionMode == UpcSelectionMode) {
      for (const auto& fv0 : fv0s) {
        recordFIT(fv0.bcId(), MatchingFITRecord{fv0s, ft0s, fdds, fv0.index(), -1, -1});
        auto gbc = static_cast<int64_t>(bcs.rawIteratorAt(fv0.bcId()).globalBC());
        float sum = udhelpers::FV0AmplitudeA(fv0);
        auto [it, inserted] = fv0Amplitudes.emplace(gbc, sum);
        if (!inserted) {
          it->second = std::max(it->second, sum);
        }
      }
      for (const auto& ft0 : ft0s) {
        recordFIT(ft0.bcId(), MatchingFITRecord{fv0s, ft0s, fdds, -1, ft0.index(), -1});
      }
      for (const auto& fdd : fdds) {
        recordFIT(fdd.bcId(), MatchingFITRecord{fv0s, ft0s, fdds, -1, -1, fdd.index()});
      }
    }
    std::unordered_map<int64_t, std::array<double, 4>> mixingMchDirectionsPV;
    const bool mixingNeedsDirectionPV = mixingCuts.cfgEnableEventMixing && mixingCuts.cfgMixingDepth != 0 &&
                                        (mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV >= 0.f || mixingCuts.cfgMixingMaxDeltaMchEtaPV >= 0.f || mixingMultiplicityCuts.cfgMixingMaxDeltaMchMomentumRel >= 0.f || mixingMultiplicityCuts.cfgMixingMaxDeltaMchPtRel >= 0.);
    std::unordered_map<int64_t, bool> selectedMch;
    int numSelectedMch = 0;
    bool qualityFieldReady = false;
    for (const auto& mch : muontracks) {
      if (mch.trackType() <= o2::aod::fwdtrack::GlobalForwardTrack) {
        continue;
      }
      bool sel8 = mch.has_collision() && collisions.rawIteratorAt(mch.collisionId()).sel8();
      auto bcIt = mchGlobalBCs.find(mch.globalIndex());
      bool hasBC = bcIt != mchGlobalBCs.end();
      int64_t globalBC = hasBC ? bcIt->second : 0;
      int numMchMid = hasBC ? numMchMidPerBC[globalBC] : 0;
      bool passFIT = true;
      bool hasFV0 = false;
      float maxFV0Amplitude = 0.f;
      if (hasBC) {
        for (auto it = cleanFITPerBC.lower_bound(globalBC - cfgUpcFV0BcWindow);
             it != cleanFITPerBC.end() && it->first <= globalBC + cfgUpcFV0BcWindow; ++it) {
          passFIT = passFIT && it->second;
        }
        auto first = fv0Amplitudes.lower_bound(globalBC - cfgUpcFV0BcWindow);
        for (auto it = first; it != fv0Amplitudes.end() && it->first <= globalBC + cfgUpcFV0BcWindow; ++it) {
          hasFV0 = true;
          maxFV0Amplitude = std::max(maxFV0Amplitude, it->second);
        }
      }
      int numMch = mch.has_collision() ? numMchPerCollision[mch.collisionId()] : 0;
      int numMft = mch.has_collision() ? numMftPerCollision[mch.collisionId()] : 0;
      bool passTrackMultiplicity = mch.has_collision() && numMch >= cfgUpcMinMchTracks && numMft >= cfgUpcMinMftTracks &&
                                   (cfgUpcMaxMchTracks < 0 || numMch <= cfgUpcMaxMchTracks) &&
                                   (!cfgUpcRequireEqualTrackCounts || numMft == numMch);
      bool passMchUnlikeSign = mch.has_collision() && numPositiveMchPerCollision[mch.collisionId()] > 0 &&
                               numNegativeMchPerCollision[mch.collisionId()] > 0;
      // Inclusive multiplicities: any number of MCH opposite-sign pairs is allowed.
      // Individual J/psi mass tags are attached downstream per participating MCH.
      bool upc = passTrackMultiplicity && (!cfgUpcRequireUnlikeSign || passMchUnlikeSign) && hasBC && numMchMid >= cfgUpcMinMchMid &&
                 (cfgUpcMaxMchMid < 0 || numMchMid <= cfgUpcMaxMchMid) && passFIT;
      bool selected = matchingEventSelectionMode == 0 || (matchingEventSelectionMode == 1 ? sel8 : upc);
      bool passMuonQuality = !mchQualityCut;
      bool propagatedMuon = false;
      float cutEta = std::numeric_limits<float>::quiet_NaN();
      float cutPt = std::numeric_limits<float>::quiet_NaN();
      // Event rejection precedes PV propagation and geometric acceptance.
      if (selected && (mchQualityCut || mchAcceptance.cfgRequireMchMftAcceptance || mchAcceptance.cfgRequireMchPVEta)) {
        cutEta = mch.eta();
        cutPt = mch.pt();
        bool valid = mch.has_collision();
        bool passMftAcceptance = !mchAcceptance.cfgRequireMchMftAcceptance;
        bool passPVEta = !mchAcceptance.cfgRequireMchPVEta;
        if (mchAcceptance.cfgRequireMchPVEta) {
          mchAcceptanceQA.fill(HIST("pvEtaStatus"), 0.);
        }
        if (mchAcceptance.cfgRequireMchMftAcceptance) {
          mchAcceptanceQA.fill(HIST("status"), 0.);
        }
        if (valid && (mchAcceptance.cfgPropMuon || mchAcceptance.cfgRequireMchMftAcceptance || mchAcceptance.cfgRequireMchPVEta)) {
          if (!qualityFieldReady) {
            initCCDB(bcs.begin());
            o2::mch::TrackExtrap::useExtrapV2();
            qualityFieldReady = true;
          }
          const auto collision = collisions.rawIteratorAt(mch.collisionId());
          auto mchAtPV = globalMatchingPlane.fwdToMch(getTrackParCovFwd(mch, mch));
          valid = std::isfinite(collision.posX()) && std::isfinite(collision.posY()) && std::isfinite(collision.posZ()) &&
                  o2::mch::TrackExtrap::extrapToVertex(mchAtPV, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
          if (valid) {
            const auto atPV = globalMatchingPlane.mchToFwd(mchAtPV);
            // Preserve mchAcceptance.cfgPropMuon semantics for quality cuts independently of
            // the mandatory PV propagation used by geometric acceptance.
            if (mchAcceptance.cfgPropMuon) {
              cutEta = atPV.getEta();
              cutPt = atPV.getPt();
            }
            if (mchAcceptance.cfgRequireMchPVEta) {
              const float etaAtPV = atPV.getEta();
              passPVEta = std::isfinite(etaAtPV) &&
                          etaAtPV > static_cast<float>(mchAcceptance.cfgMchPVEtaMin.value) &&
                          etaAtPV < static_cast<float>(mchAcceptance.cfgMchPVEtaMax.value);
            }
            propagatedMuon = true;
            if (mixingNeedsDirectionPV) {
              mixingMchDirectionsPV[mch.globalIndex()] = {atPV.getPhi(), atPV.getEta(), atPV.getP(), atPV.getPt()};
            }
            if (mchAcceptance.cfgRequireMchMftAcceptance) {
              passMftAcceptance = crossesNominalMftAcceptance(mchAtPV);
            }
          }
        }
        if (mchAcceptance.cfgRequireMchMftAcceptance) {
          mchAcceptanceQA.fill(HIST("status"), !valid ? 1. : passMftAcceptance ? 3.
                                                                               : 2.);
        }
        if (mchAcceptance.cfgRequireMchPVEta) {
          mchAcceptanceQA.fill(HIST("pvEtaStatus"), !valid ? 1. : passPVEta ? 3.
                                                                            : 2.);
        }
        if (mchQualityCut) {
          std::array<float, VarManager::kNVars> values{};
          values[VarManager::kEta] = cutEta;
          values[VarManager::kPt] = cutPt;
          values[VarManager::kMuonRAtAbsorberEnd] = mch.rAtAbsorberEnd();
          values[VarManager::kMuonPDca] = mch.pDca();
          values[VarManager::kMuonChi2] = mch.chi2();
          values[VarManager::kMuonChi2MatchMCHMID] = mch.chi2MatchMCHMID();
          const bool validQuality = valid && std::isfinite(cutEta) && std::isfinite(cutPt) &&
                                    std::isfinite(mch.rAtAbsorberEnd()) && std::isfinite(mch.pDca()) &&
                                    std::isfinite(mch.chi2()) && std::isfinite(mch.chi2MatchMCHMID());
          passMuonQuality = validQuality && mchQualityCut->IsSelected(values.data());
        }
        selected = valid && passMuonQuality && passMftAcceptance && passPVEta;
      }
      selectedMch.emplace(mch.globalIndex(), selected);
      numSelectedMch += static_cast<int>(selected);
      if (inferenceConfig.cfgTrainingMode) {
        tableMatchPlaneEventSels(mch.globalIndex(), mch.collisionId(), matchingEventSelectionMode, selected,
                                 hasBC, static_cast<uint64_t>(globalBC), numMchMid, sel8, hasFV0, maxFV0Amplitude,
                                 numMch, numMft, passTrackMultiplicity, passMchUnlikeSign, passFIT, passMuonQuality, propagatedMuon, cutEta, cutPt, cfgMftAbsorberStepCm.value);
      }
    }
    LOGF(info, "Matching-plane event selection mode=%d: %d/%zu MCH tracks selected", matchingEventSelectionMode, numSelectedMch, selectedMch.size());

    // Rejected DFs still publish selection QA, but perform no geometry/field
    // setup, track conversion or propagation in the early-filter path.
    if (filterBeforeMatching && numSelectedMch == 0) {
      LOGF(info, "Matching-plane early selection: no accepted MCH tracks; skip all pair propagation");
      if (cfgSaveMatchingAudit) {
        const uint64_t sourceFrame = (origins.size() != 0) ? origins.begin().dataframeID() : frame;
        for (auto const& mch : muontracks) {
          if (selectedMch.contains(mch.globalIndex())) {
            mchPopulation(sourceFrame, mch.collisionId(), mch.globalIndex(), static_cast<int8_t>(mch.signed1Pt() > 0. ? 1 : -1), false, false,
                          NAN, NAN, NAN, NAN, 0, 0, 0, useAssociation, 0, mch.trackType() == o2::aod::fwdtrack::MuonStandaloneTrack);
          }
        }
      }
      return;
    }
    initCCDB(bcs.begin());
    o2::mch::TrackExtrap::useExtrapV2();
    const std::array<double, 3> centerMft{0., 0., -61.4};
    globalMatchingPlane.mBzAtMftCenter = fieldB->getBz(centerMft.data());
    poolReset = history.setField(fieldB->getFactorSolenoid(), fieldB->getFactorDipole()) || poolReset;

    // Same last-row-wins lookup as global-muon-matcher, including duplicate cov rows.
    std::unordered_map<int64_t, int64_t> covIndexes;
    for (const auto& cov : mftCovs) {
      covIndexes[cov.matchMFTTrackId()] = cov.globalIndex();
    }
    std::unordered_map<int64_t, o2::track::TrackParCovFwd> mftTrackPars;
    for (const auto& mft : mfttracks) {
      // Apply the MFT acceptance to the original track, before either transport.
      // Both matching planes and the data/MC workflows use this candidate set.
      const auto eta = mft.eta();
      if (!std::isfinite(eta) || eta <= cfgEtaMftlow || eta >= cfgEtaMftup) {
        continue;
      }
      auto covIt = covIndexes.find(mft.globalIndex());
      if (covIt == covIndexes.end() && cfgRequireMftCovariance) {
        continue;
      }
      if (covIt != covIndexes.end()) {
        auto cov = mftCovs.rawIteratorAt(covIt->second);
        mftTrackPars.emplace(mft.globalIndex(), getTrackParCovFwd(mft, cov));
      } else {
        // Only the fitted central parameters are needed. This zero matrix is
        // container storage, never a measured uncertainty or an input to chi2.
        SMatrix5 parameters(mft.x(), mft.y(), mft.phi(), mft.tgl(), mft.signed1Pt());
        o2::track::SMatrix55Sym zero;
        for (int i = 0; i < TrackParameterCount; ++i) {
          for (int j = 0; j <= i; ++j) {
            zero(i, j) = 0.;
          }
        }
        mftTrackPars.emplace(mft.globalIndex(), o2::track::TrackParCovFwd(mft.z(), parameters, zero, mft.chi2()));
      }
    }
    // Tagging results are unused here. Only the late-filter/reference path
    // replays these calls to preserve the original navigation history.
    if (!filterBeforeMatching) {
      std::vector<int64_t> taggedMuons;
      globalMatchingPlane.getTaggedMuons(collisions, muontracks, taggedMuons, configMuonTagging);
    }
    using PendingRow = std::tuple<int64_t, int64_t, int64_t, int64_t, float, float, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, bool>;
    std::vector<PendingRow> pendingRows;
    std::vector<PendingRow> downstreamRows;
    std::vector<int> downstreamStatuses;
    size_t numFailedMftAbsorberTransport = 0;
    size_t numFailedJointTransport = 0;
    std::array<size_t, static_cast<size_t>(MftAbsorberExtrap::Status::Count)> mftTransportStatus{};
    std::map<int64_t, std::vector<std::pair<float, int64_t>>> candidatesForOutputPropagation;
    std::vector<PendingRow> mixedPendingRows, mixedDownstreamRows;
    std::vector<int> mixedDownstreamStatuses;
    struct MixingDiagnostic {
      int64_t mch{}, mft{};
      float z{};
      bool mixed{};
      std::vector<double> mchCov, mftCov, residualCov;
      uint64_t donorFrame{};
    };
    std::vector<MixingDiagnostic> mixingDiagnostics;
    auto calculatePair = [&](float planeZ, const auto& mch, const auto& mft, bool mixed,
                             const MftMchJointTransport::Fwd* ownedMft = nullptr, uint64_t donorSourceFrame = 0) {
      const bool chiPlaneRequested = inferenceConfig.cfgTrainingMode || !inferenceConfig.cfgEvaluationMode ||
                                     MatchingInference::needsChi2Plane(inferenceConfig.cfgChi2Method.value, planeZ) ||
                                     (inferenceConfig.cfgEnableMlMatching && inference.usesMlChi2Plane(planeZ));
      const bool calculateChi2 = cfgComputeMatchingChi2 && chiPlaneRequested && (mixed ? cfgRequireMftCovariance.value : covIndexes.count(mft.globalIndex()) > 0);
      const auto mchTrackPar = getTrackParCovFwd(mch, mch);
      const auto it = mftTrackPars.find(mft.globalIndex());
      if (!ownedMft && it == mftTrackPars.end()) {
        LOGF(fatal, "Missing MFT state in pair calculation");
      }
      const auto& mftInput = ownedMft ? *ownedMft : it->second;
      auto& pairPendingRows = mixed ? mixedPendingRows : pendingRows;
      auto& pairDownstreamRows = mixed ? mixedDownstreamRows : downstreamRows;
      auto& pairDownstreamStatuses = mixed ? mixedDownstreamStatuses : downstreamStatuses;
      // Preserve the original per-pair call order. The geometry navigator
      // and field extrapolator have state; do not cache propagated MCH tracks.
      int transportStatus = 0;
      const auto appendFailedDownstream = [&](int status) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        pairDownstreamRows.emplace_back(mch.globalIndex(), mft.globalIndex(), mch.collisionId(), mft.collisionId(),
                                        -505.f, static_cast<float>(globalMatchingPlane.mBzAtMftCenter),
                                        nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
                                        nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, false);
        pairDownstreamStatuses.push_back(status);
      };
      auto mftProp = mftInput;
      auto mchProp = mchTrackPar;
      std::unique_ptr<MftMchJointTransport::Result> joint;
      if (planeZ <= DownstreamPlaneZ && cfgMftUseMchMomentum) {
        joint = std::make_unique<MftMchJointTransport::Result>();
        bool valid = false;
        if (calculateChi2) {
          valid = MftMchJointTransport::compute(globalMatchingPlane, mftInput, mchTrackPar,
                                                planeZ, cfgMftAbsorberStepCm, cfgMftUseMchCharge, cfgPairJacobianScale, *joint);
        } else {
          // Same nominal MCH-momentum-seeded trajectory, without the numeric
          // joint Jacobian, uncertainty construction or chi2 inversion.
          valid = MftMchJointTransport::evaluate(globalMatchingPlane, mftInput, mchTrackPar,
                                                 planeZ, cfgMftAbsorberStepCm, cfgMftUseMchCharge, 0., joint->state);
          joint->hasNominalState = valid;
          joint->status = valid ? 0 : 1;
        }
        if (!valid) {
          joint->chiAll = joint->chiFour = joint->chiXY = std::numeric_limits<double>::quiet_NaN();
          joint->rankAll = joint->rankFour = joint->rankXY = 0;
        }
        if (calculateChi2 && inferenceConfig.cfgTrainingMode && !mixed && selectedMch.at(mch.globalIndex())) {
          tablePairCovariances(mch.globalIndex(), mft.globalIndex(), valid, joint->status,
                               joint->hasNominalState ? joint->state.seed.getP() : std::numeric_limits<double>::quiet_NaN(),
                               mftInput.getInvQPt(), cfgMftUseMchCharge.value,
                               joint->rankAll, joint->rankFour, joint->rankXY,
                               valid ? MftMchJointTransport::flatten(joint->joint) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->residual) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->seedCov) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->seedMchCross) : std::vector<double>{},
                               valid ? MftMchJointTransport::flatten(joint->jacobian) : std::vector<double>{},
                               joint->state.sigmaLossUp, cfgPairJacobianScale.value);
        }
        transportStatus = joint->status;
        if (!valid) {
          ++numFailedJointTransport;
          // Covariance failure must not discard a successful nominal trajectory.
          // All chi2 values/ranks are invalid together, including partial chi2 success.
          joint->chiAll = joint->chiFour = joint->chiXY = std::numeric_limits<double>::quiet_NaN();
          joint->rankAll = joint->rankFour = joint->rankXY = 0;
          if (!joint->hasNominalState) {
            appendFailedDownstream(transportStatus);
            return;
          }
        }
        mftProp = joint->state.mft;
        mchProp = joint->state.mch;
      } else {

        if (planeZ <= DownstreamPlaneZ) {
          // Absorber-exit study: use mapped B and downstream material transport.
          auto mftAsMch = globalMatchingPlane.fwdToMch(mftInput);
          MftAbsorberExtrap::TransportInfo transportInfo;
          if (!MftAbsorberExtrap::extrapToAbsorberEndWithoutBranson(mftAsMch, planeZ, &transportInfo, cfgMftAbsorberStepCm)) {
            ++numFailedMftAbsorberTransport;
            ++mftTransportStatus[static_cast<size_t>(transportInfo.status)];
            appendFailedDownstream(100 + static_cast<int>(transportInfo.status));
            return;
          }
          mftProp = globalMatchingPlane.mchToFwd(mftAsMch);
        } else {
          // In particular z=-77.5: preserve the official fixed-Bz propagation.
          mftProp = globalMatchingPlane.propagateToZMft(mftInput, planeZ);
        }
        mchProp = globalMatchingPlane.propagateToZMch(mchTrackPar, planeZ);
      }
      SMatrix5 delta = mftProp.getParameters() - mchProp.getParameters();
      // Save both raw and wrapped phi residuals; all angular chi2 use the wrapped difference.
      const double wrappedPhi = TVector2::Phi_mpi_pi(delta[2]);
      if (calculateChi2 && cfgSaveMixingDiagnostics && selectedMch.at(mch.globalIndex())) {
        MixingDiagnostic diag{.mch = mch.globalIndex(), .mft = mft.globalIndex(), .z = planeZ, .mixed = mixed, .mchCov = {}, .mftCov = {}, .residualCov = {}, .donorFrame = mixed ? donorSourceFrame : frame};
        for (int r = 0; r < TrackParameterCount; ++r) {
          for (int c = 0; c < TrackParameterCount; ++c) {
            diag.mchCov.push_back(mchProp.getCovariances()(r, c));
            diag.mftCov.push_back(mftProp.getCovariances()(r, c));
            diag.residualCov.push_back(joint ? joint->residual(r, c) : mchProp.getCovariances()(r, c) + (r == c ? mftProp.getCovariances()(r, c) : 0.));
          }
        }
        mixingDiagnostics.push_back(std::move(diag));
      }
      const double nanChi2 = std::numeric_limits<double>::quiet_NaN();
      const double chi2All = !calculateChi2 ? nanChi2 : joint ? joint->chiAll
                                                              : std::get<0>(globalMatchingPlane.mMatchingFunctionMap.at("matchALL")(mchProp, mftProp));
      const double chi2Angles = !calculateChi2 ? nanChi2 : joint ? joint->chiFour
                                                                 : std::get<0>(globalMatchingPlane.mMatchingFunctionMap.at("matchXYPhiTanl")(mchProp, mftProp));
      const double chi2XY = !calculateChi2 ? nanChi2 : joint ? joint->chiXY
                                                             : std::get<0>(globalMatchingPlane.mMatchingFunctionMap.at("matchXY")(mchProp, mftProp));
      if (calculateChi2 && !mixed && !filterBeforeMatching && planeZ == UpstreamPlaneZ) {
        candidatesForOutputPropagation[mch.globalIndex()].emplace_back(static_cast<float>(chi2All), mft.globalIndex());
      }
      bool finite = !calculateChi2 || (std::isfinite(chi2All) && std::isfinite(chi2Angles) && std::isfinite(chi2XY));
      for (int i = 0; i < TrackParameterCount; ++i) {
        finite = finite && std::isfinite(delta[i]);
      }
      auto& outputRows = planeZ == UpstreamPlaneZ ? pairPendingRows : pairDownstreamRows;
      if (planeZ == DownstreamPlaneZ) {
        pairDownstreamStatuses.push_back(transportStatus);
      }
      outputRows.emplace_back(mch.globalIndex(), mft.globalIndex(), mch.collisionId(), mft.collisionId(),
                              planeZ, static_cast<float>(globalMatchingPlane.mBzAtMftCenter),
                              mchProp.getX(), mchProp.getY(), mchProp.getPhi(), mchProp.getTanl(), mchProp.getInvQPt(),
                              mftProp.getX(), mftProp.getY(), mftProp.getPhi(), mftProp.getTanl(), mftProp.getInvQPt(),
                              delta[0], delta[1], delta[2], wrappedPhi, delta[3], delta[4],
                              mftProp.getEta() - mchProp.getEta(), mftProp.getPt() - mchProp.getPt(),
                              chi2All, chi2Angles, chi2XY, finite);
    };
    // Two passes preserve the upstream navigation/call order. Both start from
    // the original AO2D states, never from the other plane's propagated state.
    for (float const& planeZ : {-77.5f, -505.f}) {
      for (const auto& entry : globalMchIterationOrder) {
        if (filterBeforeMatching && !selectedMch.at(entry.first)) {
          continue;
        }
        // Retain the relative unordered traversal order of the accepted MCH tracks.
        const auto mch = muontracks.rawIteratorAt(entry.first);
        for (int64_t const& mftId : candidateMfts[mch.collisionId()]) {
          const auto mft = mfttracks.rawIteratorAt(mftId - mfttracks.offset());
          auto it = mftTrackPars.find(mftId);
          if (it == mftTrackPars.end() || !timeCompatible(mch, mft)) {
            continue;
          }
          calculatePair(planeZ, mch, mft, false);
        }
      }
    } // plane passes
    {
      LOGF(info, "MFT/MCH joint transport: failed candidate pairs=%zu, enabled=%d", numFailedJointTransport, static_cast<int>(cfgMftUseMchMomentum.value));
      LOGF(info, "MFT absorber transport: stepCm=%g, failedPairs=%zu [invalid=%zu field=%zu geometry=%zu rangeOut=%zu nonfinite=%zu stepLimit=%zu]", cfgMftAbsorberStepCm.value, numFailedMftAbsorberTransport, mftTransportStatus[1], mftTransportStatus[2], mftTransportStatus[3], mftTransportStatus[4], mftTransportStatus[5], mftTransportStatus[6]);
    }
    // The reference path replays propagation to the MFT point for compatibility.
    // The early-filter path leaves this map empty: those states are never saved.
    for (auto& [mchId, candidates] : candidatesForOutputPropagation) { // o2-linter: disable=const-ref-in-for-loop (Sorts each candidate vector in place.)
      std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
      const auto mch = muontracks.rawIteratorAt(mchId);
      const auto mchPar = getTrackParCovFwd(mch, mch);
      for (const auto& candidate : candidates) {
        const auto& mftPar = mftTrackPars.at(candidate.second);
        auto mchAtMft = globalMatchingPlane.fwdToMch(mchPar);
        o2::mch::TrackExtrap::extrapToVertex(mchAtMft, mftPar.getX(), mftPar.getY(), mftPar.getZ(), mftPar.getSigma2X(), mftPar.getSigma2Y());
      }
    }
    if (mixingCuts.cfgEnableEventMixing && mixingCuts.cfgMixingDepth != 0) {
      using PoolEntry = MixingPoolEntry;
      struct MixedPair {
        int64_t mch{}, mft{}, donorMch{}, donorCollision{};
        int group{}, groups{}, attemptsTarget{}, attemptsDonor{};
        mixing::Deltas delta;
        int mftTarget{}, mftDonor{}, mchTarget{}, mchDonor{};
        double centTarget{}, centDonor{};
        int cbTarget{}, cbDonor{};
        uint64_t donorFrame{};
        const OwnedMixMft* ownedMft{};
      };
      std::vector<PoolEntry> pool;
      // Current-DF donors are followed by explicitly owned historical events.
      // Use only selected donor events and MFTs associated with that donor collision.
      for (const auto& mch : muontracks) {
        if (mch.trackType() <= o2::aod::fwdtrack::GlobalForwardTrack || !mch.has_collision() || !selectedMch.at(mch.globalIndex())) {
          continue;
        }
        const auto col = collisions.rawIteratorAt(mch.collisionId());
        const auto bc = bcs.rawIteratorAt(col.bcId());
        PoolEntry entry;
        entry.event = {.collision = mch.collisionId(), .mch = mch.globalIndex(), .bc = static_cast<int64_t>(bc.globalBC()), .run = bc.runNumber(), .attempts = 0, .phi = std::atan2(mch.y(), mch.x()), .radius = std::hypot(mch.x(), mch.y()), .z = col.posZ()};
        entry.event.frame = frame;
        if (mixingNeedsDirectionPV) {
          auto direction = mixingMchDirectionsPV.find(mch.globalIndex());
          if (direction == mixingMchDirectionsPV.end()) {
            // Also use PV directions when quality cuts are disabled or mchAcceptance.cfgPropMuon=false.
            // This propagation is exclusively for the mixing pool; same-event selection is unchanged.
            auto mchAtPV = globalMatchingPlane.fwdToMch(getTrackParCovFwd(mch, mch));
            if (!o2::mch::TrackExtrap::extrapToVertex(mchAtPV, col.posX(), col.posY(), col.posZ(), col.covXX(), col.covYY())) {
              continue;
            }
            const auto atPV = globalMatchingPlane.mchToFwd(mchAtPV);
            direction = mixingMchDirectionsPV.emplace(mch.globalIndex(), std::array<double, 4>{atPV.getPhi(), atPV.getEta(), atPV.getP(), atPV.getPt()}).first;
          }
          entry.event.mchPhiPV = direction->second[0];
          entry.event.mchEtaPV = direction->second[1];
          entry.event.mchMomentumPV = direction->second[2];
          entry.event.mchPtPV = direction->second[3];
          if ((mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV >= 0.f && !std::isfinite(entry.event.mchPhiPV)) ||
              (mixingCuts.cfgMixingMaxDeltaMchEtaPV >= 0.f && !std::isfinite(entry.event.mchEtaPV))) {
            continue;
          }
        }
        if (auto count = mixingCbMultiplicities.find(mch.collisionId()); count != mixingCbMultiplicities.end()) {
          entry.event.cbMultiplicity = count->second;
        }
        entry.event.mftMultiplicity = numMftPerCollision[mch.collisionId()];
        entry.event.mchMultiplicity = numMchPerCollision[mch.collisionId()];
        if (auto cent = mixingCentralities.find(mch.collisionId()); cent != mixingCentralities.end()) {
          entry.event.centrality = cent->second;
        }
        for (const auto& mft : mfttracks) {
          if (!timeCompatible(mch, mft)) {
            continue;
          }
          // Official cfgMatchAllTracks=false definition: all time-compatible MFTs,
          // including those without saved covariance, enter the attempts estimator.
          ++entry.event.attempts;
          if (candidateMfts[mch.collisionId()].contains(mft.globalIndex()) && mftTrackPars.contains(mft.globalIndex())) {
            entry.mfts.push_back({mft.globalIndex(), mft.collisionId(), mftTrackPars.at(mft.globalIndex()), referenceBCs.at(mft.globalIndex()), mft.trackTime(), mft.trackTimeRes()});
          }
        }
        pool.push_back(std::move(entry));
      }
      mixing::Config config{.depth = mixingCuts.cfgMixingDepth.value, .minDeltaBc = mixingCuts.cfgMinDeltaBc.value, .maxDeltaPhi = mixingCuts.cfgMaxDeltaPhi.value, .maxDeltaR = mixingCuts.cfgMaxDeltaR.value, .maxDeltaAttempts = mixingCuts.cfgMaxDeltaAttempts.value, .maxDeltaZ = mixingCuts.cfgMaxDeltaZ.value, .symmetricCuts = mixingCuts.cfgMixingSymmetricCuts.value, .maxDeltaMftMultiplicityRel = mixingMultiplicityCuts.cfgMixingMaxDeltaMftMultiplicityRel.value, .maxDeltaMchMultiplicity = mixingMultiplicityCuts.cfgMixingMaxDeltaMchMultiplicity.value, .maxDeltaCentrality = mixingCuts.cfgMixingMaxDeltaCentrality.value, .maxDeltaMchPhiPV = mixingMultiplicityCuts.cfgMixingMaxDeltaMchPhiPV.value, .maxDeltaMchMomentumRel = mixingMultiplicityCuts.cfgMixingMaxDeltaMchMomentumRel.value, .maxDeltaMchEtaPV = mixingCuts.cfgMixingMaxDeltaMchEtaPV.value, .maxDeltaCbMultiplicityRel = mixingMultiplicityCuts.cfgMixingMaxDeltaCbMultiplicityRel.value, .mftMultiplicityStrict = mixingMultiplicityCuts.cfgMixingMftMultiplicityStrict.value, .mchRadiusStrict = mixingMultiplicityCuts.cfgMixingMchRadiusStrict.value, .mchEtaStrict = mixingMultiplicityCuts.cfgMixingMchEtaStrict.value, .maxDeltaMchPtRel = mixingMultiplicityCuts.cfgMixingMaxDeltaMchPtRel.value, .cbMultiplicityStrict = mixingMultiplicityCuts.cfgMixingCbMultiplicityStrict.value, .vertexZStrict = mixingMultiplicityCuts.cfgMixingVertexZStrict.value};
      std::vector<const PoolEntry*> donors;
      donors.reserve(pool.size() + history.states);
      for (const auto& entry : pool) {
        donors.push_back(&entry);
      }
      if (mixingPersistence.cfgMixingAcrossDF) {
        for (const auto& event : history.events) {
          for (const auto& entry : event.anchors) {
            donors.push_back(&entry);
          }
        }
      }
      std::vector<MixedPair> mixedPairs;
      for (const auto& target : pool) {
        std::set<std::pair<uint64_t, int64_t>> usedCollisions;
        std::set<std::pair<uint64_t, int64_t>> usedMfts;
        const auto targetMch = muontracks.rawIteratorAt(target.event.mch);
        const size_t begin = mixedPairs.size();
        int group = 0;
        for (auto const& donorPtr : donors) {
          const auto& donor = *donorPtr;
          if (config.depth >= 0 && group >= config.depth) {
            break;
          }
          if (donor.mfts.empty() || (usedCollisions.contains({donor.event.frame, donor.event.collision}))) {
            continue;
          }
          mixing::Deltas delta;
          if (!mixing::compatible(target.event, donor.event, config, delta)) {
            continue;
          }
          // Reassociation can place one physical MFT in several collision pools.
          // Reject the whole donor group if any borrowed track can belong to the
          // target, or was already borrowed from another accepted donor anchor.
          const bool overlapsTarget = std::any_of(donor.mfts.begin(), donor.mfts.end(), [&](const auto& mft) {
            const double dt = mft.time - targetMch.trackTime() +
                              (mft.referenceBC - target.event.bc) * o2::constants::lhc::LHCBunchSpacingNS;
            const double resolution = mft.timeResolution + targetMch.trackTimeRes();
            return (donor.event.frame == target.event.frame &&
                    (mft.collision == target.event.collision || candidateMfts[target.event.collision].contains(mft.id))) ||
                   usedMfts.contains({donor.event.frame, mft.id}) ||
                   (std::isfinite(dt) && std::isfinite(resolution) && resolution >= 0. && std::abs(dt) <= resolution);
          });
          if (overlapsTarget) {
            continue;
          }
          usedCollisions.insert({donor.event.frame, donor.event.collision});
          for (const auto& mft : donor.mfts) {
            usedMfts.insert({donor.event.frame, mft.id});
            mixedPairs.push_back({target.event.mch, mft.id, donor.event.mch, donor.event.collision,
                                  group, 0, target.event.attempts, donor.event.attempts, delta,
                                  target.event.mftMultiplicity, donor.event.mftMultiplicity, target.event.mchMultiplicity, donor.event.mchMultiplicity,
                                  target.event.centrality, donor.event.centrality, target.event.cbMultiplicity, donor.event.cbMultiplicity, donor.event.frame, &mft});
          }
          ++group;
        }
        for (size_t i = begin; i < mixedPairs.size(); ++i) {
          mixedPairs[i].groups = group;
        }
      }
      // Separate passes preserve same-event output exactly, including navigator history.
      // Cross-event timing is deliberately NOT reapplied: donor membership was checked above.
      for (float const& planeZ : {-77.5f, -505.f}) {
        for (const auto& pair : mixedPairs) {
          const auto mch = muontracks.rawIteratorAt(pair.mch);
          calculatePair(planeZ, mch, *pair.ownedMft, true, &pair.ownedMft->state, pair.donorFrame);
        }
      }
      if (mixedPendingRows.size() != mixedPairs.size() || mixedDownstreamRows.size() != mixedPairs.size() || mixedDownstreamStatuses.size() != mixedPairs.size()) {
        LOGF(fatal, "Mixed dual-plane row alignment failure");
      }
      for (size_t i = 0; i < mixedPairs.size(); ++i) {
        const auto& pair = mixedPairs[i];
        for (int plane = 0; plane < MatchingPlaneCount; ++plane) {
          const auto& row = plane == 0 ? mixedPendingRows[i] : mixedDownstreamRows[i];
          if (std::get<0>(row) != pair.mch || std::get<1>(row) != pair.mft || (frame == pair.donorFrame && std::get<2>(row) == std::get<3>(row))) {
            LOGF(fatal, "Mixed pair identity failure");
          }
          std::apply([&](const auto&... values) { tableMixedResiduals(values..., cfgMftAbsorberStepCm.value,
                                                                      plane == 0 ? 0 : mixedDownstreamStatuses[i], pair.donorMch, pair.donorCollision, pair.group, pair.groups,
                                                                      pair.delta.bc, pair.delta.phi, pair.delta.radius, pair.delta.attempts, pair.delta.z,
                                                                      pair.attemptsTarget, pair.attemptsDonor, pair.mftTarget, pair.mftDonor, pair.mchTarget, pair.mchDonor, pair.centTarget, pair.centDonor, pair.cbTarget, pair.cbDonor, 1. / pair.groups, frame, pair.donorFrame, mixingPersistence.state.stream); }, row);
        }
      }
      if (cfgSaveMixingDiagnostics) {
        std::set<int64_t> targets;
        for (const auto& pair : mixedPairs) {
          targets.insert(pair.mch);
        }
        for (const auto& diag : mixingDiagnostics) {
          if (targets.contains(diag.mch)) {
            tableMixingDiagnostics(diag.mch, diag.mft, diag.z, diag.mixed, diag.mchCov, diag.mftCov, diag.residualCov, frame, diag.donorFrame, mixingPersistence.state.stream);
          }
        }
      }
      const size_t crossPairs = std::count_if(mixedPairs.begin(), mixedPairs.end(), [&](const auto& p) { return p.donorFrame != frame; });
      // Copy/retain only after all pointers used above have gone out of use.
      // One FIFO item per donor collision; all its MCH anchors remain available.
      if (mixingPersistence.cfgMixingAcrossDF) {
        std::map<int64_t, OwnedMixEvent> currentEvents;
        for (auto& entry : pool) { // o2-linter: disable=const-ref-in-for-loop (Moves each donor anchor into the persistent event.)
          if (entry.mfts.empty()) {
            continue;
          }
          auto& event = currentEvents[entry.event.collision];
          event.frame = frame;
          event.collision = entry.event.collision;
          event.stateCount += entry.mfts.size();
          event.anchors.push_back(std::move(entry));
        }
        for (auto& [collision, event] : currentEvents) { // o2-linter: disable=const-ref-in-for-loop (Moves the event into the persistent pool.)
          history.push(std::move(event), mixingPersistence.cfgMixingPoolMaxEvents.value, mixingPersistence.cfgMixingPoolMaxMftStates.value);
        }
      }
      LOGF(info, "Cross-DF mixing: frame=%" PRIu64 " crossPairs=%zu retainedEvents=%zu retainedMftStates=%zu evicted=%zu oversized=%zu", static_cast<uint64_t>(frame), crossPairs, history.events.size(), history.states, history.evicted, history.oversized);
      LOGF(info, "Event mixing: selected MCH=%zu; mixed pairs=%zu; saved rows=%zu; depth=%d", pool.size(), mixedPairs.size(), 2 * mixedPairs.size(), config.depth);
    }
    // calculatePair captures and fills both vectors; Cppcheck does not model its calls.
    // cppcheck-suppress knownConditionTrueFalse
    if (pendingRows.size() != downstreamRows.size() || downstreamRows.size() != downstreamStatuses.size()) {
      LOGF(fatal, "Dual-plane row alignment failure");
    }
    // Calculate the PV state once, after pair transport, in the audit's original
    // track order. Inference and saved features must use the same values: a
    // repeated material extrapolation can differ enough to cross a tree split.
    std::map<int64_t, MatchingInference::Fwd> mchPVStates;
    if (cfgSaveMatchingAudit || inferenceConfig.cfgEvaluationMode) {
      for (auto const& mch : muontracks) {
        auto sel = selectedMch.find(mch.globalIndex());
        if (sel == selectedMch.end() || !mch.has_collision()) {
          continue;
        }
        if (!cfgSaveMatchingAudit && !sel->second) {
          continue;
        }
        auto col = collisions.rawIteratorAt(mch.collisionId() - collisions.offset());
        auto native = globalMatchingPlane.fwdToMch(getTrackParCovFwd(mch, mch));
        if (o2::mch::TrackExtrap::extrapToVertex(native, col.posX(), col.posY(), col.posZ(), col.covXX(), col.covYY())) {
          mchPVStates.emplace(mch.globalIndex(), globalMatchingPlane.mchToFwd(native));
        }
      }
    }
    std::map<int64_t, int> preselectedMftCounts;
    for (auto const& [collisionId, ids] : candidateMfts) {
      for (auto id : ids) {
        preselectedMftCounts[collisionId] += mftTrackPars.contains(id);
      }
    }
    MatchingOutcomes outcomes;
    const uint64_t sourceFrame = (origins.size() != 0) ? origins.begin().dataframeID() : frame;
    if (inferenceConfig.cfgEvaluationMode) {
      inference.process(inferenceConfig, (origins.size() != 0) ? origins.begin().dataframeID() : frame, collisions, bcs, muontracks, mfttracks, pendingRows, downstreamRows, selectedMch,
                        mixingCentralities, mixingCbMultiplicities, globalMatchingPlane, globalMatchingPlane.mBzAtMftCenter, selectedMuonTable, selectedMuonTable.assigned, dimuonVertexTable, outcomes, mchPVStates, mftTrackPars, covIndexes);
    }
    if (cfgSaveMatchingAudit) {
      std::map<int64_t, std::array<int, 3>> counts;
      auto finiteState = [](auto const& row) {
        const std::array<double, 10> values{std::get<6>(row), std::get<7>(row), std::get<8>(row), std::get<9>(row), std::get<10>(row),
                                            std::get<11>(row), std::get<12>(row), std::get<13>(row), std::get<14>(row), std::get<15>(row)};
        return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
      };
      for (size_t i = 0; i < pendingRows.size(); ++i) {
        auto const& row = pendingRows[i];
        const int64_t mchId = std::get<0>(row), mftId = std::get<1>(row), collisionId = std::get<2>(row);
        if (!selectedMch.at(mchId)) {
          continue;
        }
        const bool up = finiteState(row), down = finiteState(downstreamRows[i]);
        auto& count = counts[mchId];
        ++count[0];
        count[1] += static_cast<int>(up);
        count[2] += static_cast<int>(down);
        candidateAudit(sourceFrame, collisionId, mchId, mftId, std::get<3>(row),
                       collisionId == std::get<3>(row), associatedMfts[collisionId].contains(mftId),
                       covIndexes.contains(mftId), up, down, downstreamStatuses[i]);
      }
      // Independent of assignment/cut success, including zero-candidate MCHs.
      // Run after matching, so this QA propagation cannot affect pair transport.
      for (auto const& mch : muontracks) {
        auto sel = selectedMch.find(mch.globalIndex());
        if (sel == selectedMch.end()) {
          continue;
        }
        double pt = NAN, eta = NAN, phi = NAN, momentum = NAN;
        bool hasPV = false;
        if (auto state = mchPVStates.find(mch.globalIndex()); state != mchPVStates.end()) {
          auto const& pv = state->second;
          pt = pv.getPt();
          eta = pv.getEta();
          phi = pv.getPhi();
          momentum = pv.getP();
          hasPV = std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi) && std::isfinite(momentum);
        }
        auto const& count = counts[mch.globalIndex()];
        mchPopulation(sourceFrame, mch.collisionId(), mch.globalIndex(), static_cast<int8_t>(mch.signed1Pt() > 0. ? 1 : -1), sel->second, hasPV,
                      pt, eta, phi, momentum, count[0], count[1], count[2], useAssociation, preselectedMftCounts[mch.collisionId()],
                      mch.trackType() == o2::aod::fwdtrack::MuonStandaloneTrack);
        if (sel->second && inferenceConfig.cfgEvaluationMode) {
          for (int method = 1; method <= MlMethodIndex; ++method) {
            if (MatchingInference::enabled(inferenceConfig, method)) {
              auto const& outcome = outcomes[mch.globalIndex()][method - 1];
              mchOutcomes(sourceFrame, mch.collisionId(), mch.globalIndex(), static_cast<int8_t>(method),
                          outcome.stage, outcome.mft, outcome.score, outcome.raw);
            }
          }
        }
      }
    }
    if (!inferenceConfig.cfgTrainingMode) {
      return;
    }
    // Filled by calculatePair above, including successful downstream transports.
    // cppcheck-suppress unsignedLessThanZero
    for (size_t i = 0; i < downstreamRows.size(); ++i) {
      const auto& row = downstreamRows[i];
      if (std::get<0>(row) != std::get<0>(pendingRows[i]) || std::get<1>(row) != std::get<1>(pendingRows[i])) {
        LOGF(fatal, "Dual-plane pair identity/order mismatch");
      }
      if (selectedMch.at(std::get<0>(row))) {
        std::apply([&](const auto&... values) { tableDownstreamResiduals(values..., cfgMftAbsorberStepCm.value, downstreamStatuses[i]); }, row);
      }
    }
    // Publish only after all propagation, as in global-muon-matcher.
    // Filled by calculatePair above, not an always-empty container.
    // cppcheck-suppress knownEmptyContainer
    for (const auto& row : pendingRows) {
      // The early path already rejected unselected MCH tracks before propagation.
      // Keep this check for the late-filter/reference path as well.
      if (selectedMch.at(std::get<0>(row))) {
        std::apply([this](const auto&... values) { tableMatchPlaneResiduals(values...); }, row);
      }
    }
  }
  void processMatchingPlane(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                            aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                            aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    processMatchingPlaneImpl(collisions, bcs, muontracks, mfttracks, mftCovs, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneCentrality(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                      aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                                      aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds,
                                      aod::CentFT0Cs const& centralities, aod::PVMults const& pvMultiplicities, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    if (pvMultiplicities.size() != collisions.size()) {
      LOGF(fatal, "PVMults/collision row count mismatch");
    }
    if (centralities.size() != collisions.size()) {
      LOGF(fatal, "FT0C centrality/collision row count mismatch");
    }
    for (const auto& col : collisions) {
      mixingCentralities.emplace(col.globalIndex(), centralities.rawIteratorAt(col.globalIndex() - collisions.offset()).centFT0C());
      mixingCbMultiplicities.emplace(col.globalIndex(), pvMultiplicities.rawIteratorAt(col.globalIndex() - collisions.offset()).multNTracksPV());
    }
    processMatchingPlaneImpl(collisions, bcs, muontracks, mfttracks, mftCovs, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneWithAssociation(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                           aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                           aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                                           aod::MFTTrackAssoc const& associations, aod::AmbiguousMFTTracks const& ambiguous,
                                           aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    processAssociated(collisions, bcs, muontracks, mfttracks, mftCovs, associations, ambiguous, fv0s, ft0s, fdds, origins);
  }
  template <typename Collisions>
  void processAssociated(Collisions const& collisions, aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                         aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                         aod::MFTTrackAssoc const& associations, aod::AmbiguousMFTTracks const& ambiguous,
                         aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    std::vector<std::pair<int64_t, int64_t>> pairs;
    for (auto const& association : associations) {
      pairs.emplace_back(association.collisionId(), association.mfttrackId());
    }
    processAssociationPairs(collisions, bcs, muontracks, mfttracks, mftCovs, pairs, ambiguous, fv0s, ft0s, fdds, origins);
  }
  template <typename Collisions>
  void processAssociationPairs(Collisions const& collisions, aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                               aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                               std::vector<std::pair<int64_t, int64_t>> const& pairs, aod::AmbiguousMFTTracks const& ambiguous,
                               aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    std::unordered_map<int64_t, int64_t> refs;
    for (auto const& track : ambiguous) {
      auto const& range = track.bcIds();
      if (range[0] >= static_cast<int64_t>(bcs.offset()) && range[1] >= range[0] && range[1] - static_cast<int64_t>(bcs.offset()) < bcs.size()) {
        refs.try_emplace(track.mfttrackId(), bcs.rawIteratorAt(range[0] - bcs.offset()).globalBC());
      }
    }
    processMatchingPlaneImpl(collisions, bcs, muontracks, mfttracks, mftCovs, fv0s, ft0s, fdds, origins, pairs, refs, true);
  }
  void processMatchingPlaneCentralityWithAssociation(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                     aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                     aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                                                     aod::MFTTrackAssoc const& associations, aod::AmbiguousMFTTracks const& ambiguous,
                                                     aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds,
                                                     aod::CentFT0Cs const& centralities, aod::PVMults const& pvMultiplicities, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    if (pvMultiplicities.size() != collisions.size() || centralities.size() != collisions.size()) {
      LOGF(fatal, "Centrality/PVMults collision row count mismatch");
    }
    for (auto const& col : collisions) {
      mixingCentralities.emplace(col.globalIndex(), centralities.rawIteratorAt(col.globalIndex() - collisions.offset()).centFT0C());
      mixingCbMultiplicities.emplace(col.globalIndex(), pvMultiplicities.rawIteratorAt(col.globalIndex() - collisions.offset()).multNTracksPV());
    }
    processAssociated(collisions, bcs, muontracks, mfttracks, mftCovs, associations, ambiguous, fv0s, ft0s, fdds, origins);
  }
  // Some AO2D productions have no MFT covariance table at all. Use a typed,
  // zero-row table, never synthetic measured covariance rows, for that input.
  static aod::MFTTracksCov emptyMftCovariances()
  {
    o2::framework::TableBuilder builder;
    static_cast<void>(builder.cursor<aod::MFTTracksCov>());
    return aod::MFTTracksCov(std::vector<o2::soa::ArrowTableRef>{o2::soa::ArrowTableRef{builder.finalize()}});
  }
  void processMatchingPlaneWithAssociationNoCov(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                aod::MFTTracks const& mfttracks,
                                                aod::MFTTrackAssoc const& associations, aod::AmbiguousMFTTracks const& ambiguous,
                                                aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    auto empty = emptyMftCovariances();
    processMatchingPlaneWithAssociation(collisions, bcs, muontracks, mfttracks, empty, associations, ambiguous, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneCentralityWithAssociationNoCov(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                          aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                          aod::MFTTracks const& mfttracks,
                                                          aod::MFTTrackAssoc const& associations, aod::AmbiguousMFTTracks const& ambiguous,
                                                          aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds,
                                                          aod::CentFT0Cs const& centralities, aod::PVMults const& pvMultiplicities, aod::Origins const& origins)
  {
    auto empty = emptyMftCovariances();
    processMatchingPlaneCentralityWithAssociation(collisions, bcs, muontracks, mfttracks, empty, associations, ambiguous, fv0s, ft0s, fdds, centralities, pvMultiplicities, origins);
  }
  void processMatchingPlaneWithInternalAssociation(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                   aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                   aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                                                   aod::AmbiguousMFTTracks const& ambiguous,
                                                   aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    auto pairs = ProtectedMftAssociation::run(collisions, bcs, mfttracks, ambiguous, localAssociation);
    processAssociationPairs(collisions, bcs, muontracks, mfttracks, mftCovs, pairs, ambiguous, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneCentralityWithInternalAssociation(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                             aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                             aod::MFTTracks const& mfttracks, aod::MFTTracksCov const& mftCovs,
                                                             aod::AmbiguousMFTTracks const& ambiguous,
                                                             aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds,
                                                             aod::CentFT0Cs const& centralities, aod::PVMults const& pvMultiplicities, aod::Origins const& origins)
  {
    mixingCentralities.clear();
    mixingCbMultiplicities.clear();
    if (pvMultiplicities.size() != collisions.size() || centralities.size() != collisions.size()) {
      LOGF(fatal, "Centrality/PVMults collision row count mismatch");
    }
    for (auto const& col : collisions) {
      mixingCentralities.emplace(col.globalIndex(), centralities.rawIteratorAt(col.globalIndex() - collisions.offset()).centFT0C());
      mixingCbMultiplicities.emplace(col.globalIndex(), pvMultiplicities.rawIteratorAt(col.globalIndex() - collisions.offset()).multNTracksPV());
    }
    auto pairs = ProtectedMftAssociation::run(collisions, bcs, mfttracks, ambiguous, localAssociation);
    processAssociationPairs(collisions, bcs, muontracks, mfttracks, mftCovs, pairs, ambiguous, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneWithInternalAssociationNoCov(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                        aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                        aod::MFTTracks const& mfttracks, aod::AmbiguousMFTTracks const& ambiguous,
                                                        aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds, aod::Origins const& origins)
  {
    auto empty = emptyMftCovariances();
    processMatchingPlaneWithInternalAssociation(collisions, bcs, muontracks, mfttracks, empty, ambiguous, fv0s, ft0s, fdds, origins);
  }
  void processMatchingPlaneCentralityWithInternalAssociationNoCov(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                                  aod::BCsWithTimestamps const& bcs, MyMUONs const& muontracks,
                                                                  aod::MFTTracks const& mfttracks, aod::AmbiguousMFTTracks const& ambiguous,
                                                                  aod::FV0As const& fv0s, aod::FT0s const& ft0s, aod::FDDs const& fdds,
                                                                  aod::CentFT0Cs const& centralities, aod::PVMults const& pvMultiplicities, aod::Origins const& origins)
  {
    auto empty = emptyMftCovariances();
    processMatchingPlaneCentralityWithInternalAssociation(collisions, bcs, muontracks, mfttracks, empty, ambiguous, fv0s, ft0s, fdds, centralities, pvMultiplicities, origins);
  }
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneWithInternalAssociation, "Locally guarded official MFT association; omit external associator executable", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneCentralityWithInternalAssociation, "Locally guarded official MFT association with centrality", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneWithInternalAssociationNoCov, "Locally guarded official MFT association without covariance input", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneCentralityWithInternalAssociationNoCov, "Locally guarded official MFT association with centrality and no covariance input", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneWithAssociationNoCov, "Official MFTTrackAssoc; input without an MFT covariance table; no chi2", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneCentralityWithAssociationNoCov, "Official MFTTrackAssoc with centrality; no MFT covariance input or chi2", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneWithAssociation, "Use official MFTTrackAssoc candidate collisions; run fwdtrack-to-collision-associator upstream", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneCentralityWithAssociation, "Official MFTTrackAssoc with FT0C/PV multiplicities", false);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlaneCentrality, "Matching with FT0C centrality for event mixing", false);
  PROCESS_SWITCH(MftMchMatchingTask, processLegacy, "Legacy residual tables and matching-plane conventions", true);
  PROCESS_SWITCH(MftMchMatchingTask, processMatchingPlane, "Global-muon-matcher propagation, timing and chi2; no realignment", false);
};

// Separate device: PV propagation must not change matching-plane navigation caches.
struct MchDimuonMassAtPv {
  MatchingAuditOutput<aod::MchDimuonsPV> tableMchDimuonsPV;
  o2::framework::Configurable<bool> cfgSaveMchDimuonMass{"cfgSaveMchDimuonMass", true, "Save selected same-collision MCH dimuons using momenta propagated to the primary vertex"};
  o2::framework::Configurable<int> cfgMassBins{"cfgMassBins", 600, "Dimuon mass histogram bins"};
  o2::framework::Configurable<float> cfgMassMax{"cfgMassMax", 12.f, "Dimuon mass display upper edge in GeV/c^2; higher masses kept in overflow and table"};
  o2::framework::Configurable<std::string> ccdbURL{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"}; // o2-linter: disable=name/configurable (Preserve the public ccdb-url JSON key.)
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::framework::HistogramRegistry histos{"MchDimuonMassPV", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  o2::aod::fwdtrackutils::GlobalMuonMatchingPlane conversion;
  int runNumber = -1;

  void init(o2::framework::InitContext&)
  {
    if (!cfgSaveMchDimuonMass) {
      return;
    }
    ccdb->setURL(ccdbURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (cfgMassBins < 1 || !std::isfinite(cfgMassMax.value) || cfgMassMax <= 0.f) {
      LOGF(fatal, "Invalid MCH dimuon mass histogram range or bin count");
    }
    const o2::framework::AxisSpec massAxis{cfgMassBins, 0., cfgMassMax, "m_{#mu#mu} at PV (GeV/c^{2})"};
    histos.add("hMassAll", "Selected MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassUnlikeSign", "Unlike-sign MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassPlusPlus", "Positive like-sign MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassMinusMinus", "Negative like-sign MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); MCH pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassVsPt", "Selected MCH dimuons at PV; m_{#mu#mu} (GeV/c^{2}); p_{T,#mu#mu} (GeV/c)", o2::framework::HistType::kTH2D, {massAxis, {200, 0., 20.}});
    histos.add("hPVPropagation", "MCH PV propagation; status; Tracks", o2::framework::HistType::kTH1D, {{3, 0., 3.}});
    histos.get<TH1>(HIST("hPVPropagation"))->GetXaxis()->SetBinLabel(1, "Attempted");
    histos.get<TH1>(HIST("hPVPropagation"))->GetXaxis()->SetBinLabel(2, "Succeeded");
    histos.get<TH1>(HIST("hPVPropagation"))->GetXaxis()->SetBinLabel(3, "Failed/nonfinite");
    o2::conf::ConfigurableParam::setValue("FieldOriginBias.z", "0");
  }

  void process(aod::Collisions const& collisions,
               aod::BCsWithTimestamps const& bcs,
               MyMUONs const& tracks,
               aod::MatchPlaneEventSels const& selections)
  {
    if (!cfgSaveMchDimuonMass || bcs.size() == 0) {
      return;
    }
    std::map<int64_t, std::vector<int64_t>> selectedByCollision;
    for (const auto& sel : selections) {
      if (sel.isSelected() && sel.mchCollisionId() >= 0) {
        selectedByCollision[sel.mchCollisionId()].push_back(sel.mchTrackId());
      }
    }
    bool hasPair = std::any_of(selectedByCollision.begin(), selectedByCollision.end(), [](const auto& entry) { return entry.second.size() >= MinimumDimuonTracks; });
    if (!hasPair) {
      return;
    }
    const auto& bc = bcs.begin();
    if (runNumber != bc.runNumber()) {
      auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (!grp) {
        LOGF(fatal, "No magnetic field for MCH PV propagation at BC timestamp");
      }
      o2::base::Propagator::initFieldFromGRP(grp);
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->getForTimeStamp<TGeoManager>("GLO/Config/GeometryAligned", bc.timestamp());
      }
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        LOGF(fatal, "No geometry for MCH PV propagation");
      }
      o2::mch::TrackExtrap::setField();
      o2::mch::TrackExtrap::useExtrapV2();
      runNumber = bc.runNumber();
    }
    const double muonRestMass = muonMass();
    struct MuonAtPV {
      int64_t id{};
      int charge{};
      MuonFourVector momentum;
    };
    for (auto& [collisionId, ids] : selectedByCollision) { // o2-linter: disable=const-ref-in-for-loop (Sorts and deduplicates track IDs in place.)
      if (ids.size() < MinimumDimuonTracks) {
        continue;
      }
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      const auto collision = collisions.rawIteratorAt(collisionId);
      std::vector<MuonAtPV> muons;
      for (auto const& id : ids) {
        const auto track = tracks.rawIteratorAt(id);
        if (track.trackType() <= o2::aod::fwdtrack::GlobalForwardTrack || !track.has_collision() || track.collisionId() != collisionId) {
          LOGF(fatal, "MCH selection contains an invalid track/collision association");
        }
        histos.fill(HIST("hPVPropagation"), 0.5);
        // Start from the ORIGINAL AO2D MCH state, never from matching-plane states.
        auto mch = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(track, track));
        bool success = o2::mch::TrackExtrap::extrapToVertex(mch, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
        if (!success) {
          histos.fill(HIST("hPVPropagation"), 2.5);
          continue;
        }
        auto pv = conversion.mchToFwd(mch);
        double px = pv.getPx(), py = pv.getPy(), pz = pv.getPz();
        if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(pv.getInvQPt()) || pv.getInvQPt() == 0.) {
          histos.fill(HIST("hPVPropagation"), 2.5);
          continue;
        }
        histos.fill(HIST("hPVPropagation"), 1.5);
        MuonFourVector momentum;
        momentum = muonMomentum(px, py, pz, muonRestMass);
        muons.push_back({id, pv.getInvQPt() > 0. ? 1 : -1, momentum});
      }
      for (size_t i = 0; i < muons.size(); ++i) {
        for (size_t j = i + 1; j < muons.size(); ++j) {
          const auto& a = muons[i];
          const auto& b = muons[j];
          MuonFourVector pair = a.momentum + b.momentum;
          const double mass = muonInvariantMass(pair), pt = pair.Pt(), rapidity = muonRapidity(pair);
          if (!std::isfinite(mass) || !std::isfinite(pt) || !std::isfinite(rapidity)) {
            continue;
          }
          histos.fill(HIST("hMassAll"), mass);
          histos.fill(HIST("hMassVsPt"), mass, pt);
          if (a.charge * b.charge < 0) {
            histos.fill(HIST("hMassUnlikeSign"), mass);
          } else if (a.charge > 0) {
            histos.fill(HIST("hMassPlusPlus"), mass);
          } else {
            histos.fill(HIST("hMassMinusMinus"), mass);
          }
          tableMchDimuonsPV(collisionId, a.id, b.id, static_cast<int8_t>(a.charge + b.charge),
                            collision.posX(), collision.posY(), collision.posZ(),
                            a.momentum.Px(), a.momentum.Py(), a.momentum.Pz(),
                            b.momentum.Px(), b.momentum.Py(), b.momentum.Pz(), mass, pt, rapidity);
        }
      }
    }
  }
};

// ML inference/one-to-one assignment is performed externally. A per-AO2D CSV
// selects exact candidates; no MC labels are consumed by this device.
struct MlMatchedDimuonMassAtPv {
  MatchingAuditOutput<aod::MlMatchedMuonsPV> muonTable;
  MatchingAuditOutput<aod::MlMatchedDimuonsPV> pairTable;
  o2::framework::Configurable<bool> cfgSaveMlDimuons{"cfgSaveMlDimuons", false, "Build PV dimuons from externally ML-selected one-to-one MFT-MCH pairs"};
  o2::framework::Configurable<std::string> cfgMlSelectionFile{"cfgMlSelectionFile", "", "Per-input CSV: frameId,collisionId,mchTrackId,mftTrackId,goodScore (no MC truth)"};
  o2::framework::Configurable<double> cfgMinMlGoodScore{"cfgMinMlGoodScore", 0., "Minimum supplied good-match score; source score is not a calibrated correct-match probability"};
  o2::framework::Configurable<int> cfgMassBins{"cfgMassBins", 600, "Dimuon mass bins"};
  o2::framework::Configurable<float> cfgMassMax{"cfgMassMax", 12.f, "Mass histogram upper edge; all finite masses retained in table"};
  o2::framework::Configurable<std::string> ccdbURL{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"}; // o2-linter: disable=name/configurable (Preserve the public ccdb-url JSON key.)
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::framework::HistogramRegistry histos{"MlMatchedDimuonPV", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  struct Choice {
    int64_t collision, mch, mft;
    double score;
  };
  std::map<uint64_t, std::vector<Choice>> choices;
  std::set<uint64_t> processedFrames;
  size_t expected = 0, consumed = 0;
  int runNumber = -1;
  double bz = 0.;
  o2::aod::fwdtrackutils::GlobalMuonMatchingPlane conversion;

  void init(o2::framework::InitContext&)
  {
    if (static_cast<int>(doprocessWithCovariance) + static_cast<int>(doprocessNoCovariance) != 1) {
      LOGF(fatal, "Enable exactly one external-CSV dimuon input process");
    }
    if (!cfgSaveMlDimuons) {
      return;
    }
    if (cfgMlSelectionFile.value.empty() || cfgMassBins < 1 || !std::isfinite(cfgMassMax.value) || cfgMassMax <= 0.f ||
        !std::isfinite(cfgMinMlGoodScore.value) || cfgMinMlGoodScore < 0. || cfgMinMlGoodScore > 1.) {
      LOGF(fatal, "Invalid ML dimuon configuration: a per-input selection CSV is required");
    }
    std::ifstream input(cfgMlSelectionFile.value);
    if (!input) {
      LOGF(fatal, "Cannot read ML selection file %s", cfgMlSelectionFile.value.c_str());
    }
    std::string line;
    size_t lineNumber = 0;
    std::map<std::pair<uint64_t, int64_t>, std::set<int64_t>> usedMch, usedMft;
    while (std::getline(input, line)) {
      ++lineNumber;
      if (line.empty() || line[0] == '#') {
        continue;
      }
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream row(line);
      uint64_t frame{};
      Choice c{};
      std::string extra;
      if (!(row >> frame >> c.collision >> c.mch >> c.mft >> c.score) || (row >> extra) ||
          c.collision < 0 || c.mch < 0 || c.mft < 0 || !std::isfinite(c.score) || c.score < 0. || c.score > 1.) {
        LOGF(fatal, "Malformed ML selection CSV at line %zu", lineNumber);
      }
      if (c.score < cfgMinMlGoodScore) {
        continue;
      }
      const auto event = std::make_pair(frame, c.collision);
      if (!usedMch[event].insert(c.mch).second || !usedMft[event].insert(c.mft).second) {
        LOGF(fatal, "ML dimuons require one-to-one selection within each reconstructed event");
      }
      choices[frame].push_back(c);
      ++expected;
    }
    ccdb->setURL(ccdbURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    const o2::framework::AxisSpec massAxis{cfgMassBins, 0., cfgMassMax, "m_{#mu#mu} at PV (GeV/c^{2})"};
    histos.add("hMassAll", "ML-selected matched dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassUnlikeSign", "ML-selected unlike-sign dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassPlusPlus", "ML-selected positive like-sign dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassMinusMinus", "ML-selected negative like-sign dimuons; m_{#mu#mu} (GeV/c^{2}); Pairs", o2::framework::HistType::kTH1D, {massAxis});
    histos.add("hMassVsPt", "ML-selected dimuons; m_{#mu#mu} (GeV/c^{2}); p_{T,#mu#mu} (GeV/c)", o2::framework::HistType::kTH2D, {massAxis, {200, 0., 20.}});
    histos.add("hMuonStatus", "ML-selected PV reconstruction; status; Tracks", o2::framework::HistType::kTH1D, {{4, 0., 4.}});
  }

  void processWithCovariance(aod::Origins const& origins, aod::Collisions const& collisions,
                             aod::BCsWithTimestamps const&, MyMUONs const& mchTracks,
                             aod::MFTTracks const& mftTracks, aod::MFTTracksCov const& mftCovs,
                             aod::MatchPlaneResiduals const& candidates)
  {
    if (!cfgSaveMlDimuons) {
      return;
    }
    if (origins.size() == 0) {
      LOGF(fatal, "ML dimuon selection needs the original AO2D frame identity");
    }
    const uint64_t frame = origins.begin().dataframeID();
    auto found = choices.find(frame);
    if (found == choices.end()) {
      return;
    }
    if (!processedFrames.insert(frame).second) {
      LOGF(fatal, "Repeated AO2D frame in ML dimuon input");
    }
    std::set<std::tuple<int64_t, int64_t, int64_t>> available;
    for (const auto& row : candidates) {
      if (row.mchCollisionId() == row.mftCollisionId()) {
        available.emplace(row.mchCollisionId(), row.mchTrackId(), row.mftTrackId());
      }
    }
    std::unordered_map<int64_t, int64_t> covIndex;
    for (const auto& cov : mftCovs) {
      covIndex[cov.matchMFTTrackId()] = cov.globalIndex();
    }
    struct Muon {
      Choice choice{};
      int charge{};
      MuonFourVector momentum;
    };
    std::map<int64_t, std::vector<Muon>> muons;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double mass = muonMass();
    for (const auto& c : found->second) {
      if (!available.contains({c.collision, c.mch, c.mft})) {
        LOGF(fatal, "ML-selected pair is absent from same-collision matching candidates: frame=%" PRIu64 " MCH=%" PRId64 " MFT=%" PRId64 "", static_cast<uint64_t>(frame), static_cast<int64_t>(c.mch), static_cast<int64_t>(c.mft));
      }
      const auto collision = collisions.rawIteratorAt(c.collision - collisions.offset());
      const auto mch = mchTracks.rawIteratorAt(c.mch - mchTracks.offset());
      const auto mft = mftTracks.rawIteratorAt(c.mft - mftTracks.offset());
      if (mch.collisionId() != c.collision || mft.collisionId() != c.collision || mch.trackType() <= o2::aod::fwdtrack::GlobalForwardTrack) {
        LOGF(fatal, "Invalid selected MFT-MCH collision or track type");
      }
      const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
        if (!grp) {
          LOGF(fatal, "Missing magnetic field for ML dimuon reconstruction");
        }
        o2::base::Propagator::initFieldFromGRP(grp);
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
          ccdb->getForTimeStamp<TGeoManager>("GLO/Config/GeometryAligned", bc.timestamp());
        }
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
          LOGF(fatal, "Missing geometry for ML dimuons");
        }
        o2::mch::TrackExtrap::setField();
        o2::mch::TrackExtrap::useExtrapV2();
        const std::array<double, 3> center{0., 0., -61.4};
        auto* field = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
        if (!field) {
          LOGF(fatal, "Expected an O2 magnetic field for ML dimuons");
        }
        bz = field->getBz(center.data());
        runNumber = bc.runNumber();
      }
      const bool hasCov = covIndex.contains(c.mft);
      auto publishFailure = [&](int status) {
        muonTable(frame, c.collision, c.mch, c.mft, c.score, hasCov, status, static_cast<int8_t>(0), nan, nan, nan);
        histos.fill(HIST("hMuonStatus"), status + 0.5);
      };
      ++consumed;
      auto mchAtPV = conversion.fwdToMch(o2::aod::fwdtrackutils::getTrackParCovFwd(mch, mch));
      if (!o2::mch::TrackExtrap::extrapToVertex(mchAtPV, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY())) {
        publishFailure(1);
        continue;
      }
      const auto momentumAtPV = conversion.mchToFwd(mchAtPV);
      if (!std::isfinite(momentumAtPV.getInvQPt()) || momentumAtPV.getInvQPt() == 0. || !std::isfinite(momentumAtPV.getTgl())) {
        publishFailure(2);
        continue;
      }
      o2::track::SMatrix55Sym zero;
      for (int i = 0; i < TrackParameterCount; ++i) {
        for (int j = 0; j <= i; ++j) {
          zero(i, j) = 0.;
        }
      }
      o2::track::SMatrix5 parameters(mft.x(), mft.y(), mft.phi(), mft.tgl(), mft.signed1Pt());
      o2::track::TrackParCovFwd mftState(mft.z(), parameters, zero, mft.chi2());
      if (hasCov) {
        auto cov = mftCovs.rawIteratorAt(covIndex.at(c.mft) - mftCovs.offset());
        mftState = o2::aod::fwdtrackutils::getTrackParCovFwd(mft, cov);
      }
      // PWGEM construction: MCH momentum at PV, MFT position/direction at IU,
      // MCH charge. No absorber-end hybrid state is used. Without MFT covariance,
      // zero is only container storage: only central kinematics are published.
      auto refit = o2::aod::fwdtrackutils::refitGlobalMuonCov(momentumAtPV, mftState);
      const auto pv = o2::aod::fwdtrackutils::propagateTrackParCovFwd(refit, static_cast<uint8_t>(0), collision,
                                                                      o2::aod::fwdtrackutils::propagationPoint::kToVertex, -77.5f, bz);
      const double px = pv.getPx(), py = pv.getPy(), pz = pv.getPz();
      if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(pv.getInvQPt()) || pv.getInvQPt() == 0.) {
        publishFailure(3);
        continue;
      }
      const int charge = pv.getInvQPt() > 0. ? 1 : -1;
      MuonFourVector momentum;
      momentum = muonMomentum(px, py, pz, mass);
      muonTable(frame, c.collision, c.mch, c.mft, c.score, hasCov, 0, static_cast<int8_t>(charge), px, py, pz);
      histos.fill(HIST("hMuonStatus"), 0.5);
      muons[c.collision].push_back({c, charge, momentum});
    }
    for (const auto& [collisionId, tracks] : muons) {
      for (size_t i = 0; i < tracks.size(); ++i) {
        for (size_t j = i + 1; j < tracks.size(); ++j) {
          const auto& a = tracks[i];
          const auto& b = tracks[j];
          const auto pair = a.momentum + b.momentum;
          if (!std::isfinite(muonInvariantMass(pair)) || !std::isfinite(pair.Pt()) || !std::isfinite(muonRapidity(pair))) {
            continue;
          }
          pairTable(frame, collisionId, a.choice.mch, a.choice.mft, b.choice.mch, b.choice.mft,
                    a.choice.score, b.choice.score, static_cast<int8_t>(a.charge + b.charge),
                    a.momentum.Px(), a.momentum.Py(), a.momentum.Pz(), b.momentum.Px(), b.momentum.Py(), b.momentum.Pz(),
                    muonInvariantMass(pair), pair.Pt(), muonRapidity(pair));
          histos.fill(HIST("hMassAll"), muonInvariantMass(pair));
          histos.fill(HIST("hMassVsPt"), muonInvariantMass(pair), pair.Pt());
          if (a.charge * b.charge < 0) {
            histos.fill(HIST("hMassUnlikeSign"), muonInvariantMass(pair));
          } else if (a.charge > 0) {
            histos.fill(HIST("hMassPlusPlus"), muonInvariantMass(pair));
          } else {
            histos.fill(HIST("hMassMinusMinus"), muonInvariantMass(pair));
          }
        }
      }
    }
  }
  void processNoCovariance(aod::Origins const& origins, aod::Collisions const& collisions,
                           aod::BCsWithTimestamps const& bcs, MyMUONs const& mchTracks,
                           aod::MFTTracks const& mftTracks, aod::MatchPlaneResiduals const& candidates)
  {
    auto empty = MftMchMatchingTask::emptyMftCovariances();
    processWithCovariance(origins, collisions, bcs, mchTracks, mftTracks, empty, candidates);
  }
  PROCESS_SWITCH(MlMatchedDimuonMassAtPv, processWithCovariance, "Legacy external-CSV dimuons with optional stored MFT covariance rows", true);
  PROCESS_SWITCH(MlMatchedDimuonMassAtPv, processNoCovariance, "Legacy external-CSV device without MFT covariance input; also use when this optional device is disabled on no-covariance AO2D", false);
  void endOfStream(o2::framework::EndOfStreamContext&)
  {
    if (cfgSaveMlDimuons && consumed != expected) {
      LOGF(fatal, "ML selection/input mismatch: consumed %zu of %zu selected muons", consumed, expected);
    }
  }
};

// Join PV mass tags to residual rows downstream, without feeding PV propagation
// back into the matching device. Each candidate inherits the tag of its MCH pair.
template <bool Sparse>
struct MchJpsiMatchFlags {
  o2::framework::Produces<aod::MatchPlaneResiduals> residuals;
  o2::framework::Produces<aod::JpsiMatchPlaneResiduals> jpsiResiduals;
  o2::framework::Configurable<double> cfgJpsiMassMin{"cfgJpsiMassMin", 2.8, "Inclusive lower bound on OS MCH-pair PV invariant mass in GeV/c^2"};
  o2::framework::Configurable<double> cfgJpsiMassMax{"cfgJpsiMassMax", 3.3, "Inclusive upper bound on OS MCH-pair PV invariant mass in GeV/c^2"};
  void init(o2::framework::InitContext&)
  {
    if (!std::isfinite(cfgJpsiMassMin.value) || !std::isfinite(cfgJpsiMassMax.value) || cfgJpsiMassMin < 0. || cfgJpsiMassMax <= cfgJpsiMassMin) {
      LOGF(fatal, "Invalid J/psi PV mass window");
    }
  }
  void process(aod::RawMatchPlaneResiduals const& raw,
               aod::MatchPlaneEventSels const& selections,
               MatchingAuditInput<aod::MchDimuonsPV, Sparse> const& pairs)
  {
    std::map<int64_t, std::set<int64_t>> eligible;
    for (const auto& sel : selections) {
      if (sel.isSelected() && sel.mchCollisionId() >= 0) {
        eligible[sel.mchCollisionId()].insert(sel.mchTrackId());
      }
    }
    struct TrackTag {
      double mass = std::numeric_limits<double>::quiet_NaN();
      int64_t partner = -1;
      bool representativeInWindow = false;
      std::set<int64_t> jpsiPartners;
    };
    // The key contains the MCH ID: a signal elsewhere in the collision must not
    // tag unrelated MCH tracks. Multiple qualifying partners do not duplicate rows.
    std::map<std::pair<int64_t, int64_t>, TrackTag> tags;
    const double targetMass = particleMass(o2::constants::physics::kJPsi);
    auto addPartner = [&](int64_t collision, int64_t mch, int64_t partner, double mass) {
      auto& tag = tags[{collision, mch}];
      bool inWindow = mass >= cfgJpsiMassMin && mass <= cfgJpsiMassMax;
      if (inWindow) {
        tag.jpsiPartners.insert(partner);
      }
      // One scalar mass per MFT-MCH row: prioritize an in-window partner, then
      // the mass closest to PDG J/psi. Use partner ID to break exact ties.
      if (tag.partner < 0 || (inWindow && !tag.representativeInWindow) ||
          (inWindow == tag.representativeInWindow &&
           (std::abs(mass - targetMass) < std::abs(tag.mass - targetMass) ||
            (std::abs(mass - targetMass) == std::abs(tag.mass - targetMass) && partner < tag.partner)))) {
        tag.mass = mass;
        tag.partner = partner;
        tag.representativeInWindow = inWindow;
      }
    };
    for (const auto& pair : pairs) {
      auto event = eligible.find(pair.collisionId());
      if (event == eligible.end() || pair.pairCharge() != 0 || !std::isfinite(pair.mass()) ||
          (!event->second.contains(pair.mchTrackId1())) || (!event->second.contains(pair.mchTrackId2()))) {
        continue;
      }
      addPartner(pair.collisionId(), pair.mchTrackId1(), pair.mchTrackId2(), pair.mass());
      addPartner(pair.collisionId(), pair.mchTrackId2(), pair.mchTrackId1(), pair.mass());
    }
    for (const auto& row : raw) {
      auto found = tags.find({row.mchCollisionId(), row.mchTrackId()});
      bool valid = found != tags.end();
      double mass = valid ? found->second.mass : std::numeric_limits<double>::quiet_NaN();
      bool isJpsi = valid && !found->second.jpsiPartners.empty();
      int64_t partner = valid ? found->second.partner : -1;
      int numPartners = valid ? static_cast<int>(found->second.jpsiPartners.size()) : 0;
      residuals(row.mchTrackId(), row.mftTrackId(), row.mchCollisionId(), row.mftCollisionId(), row.matchingPlaneZ(), row.bzAtMftCenter(), row.mchX(), row.mchY(), row.mchPhi(), row.mchTanl(), row.mchInvQPt(), row.mftX(), row.mftY(), row.mftPhi(), row.mftTanl(), row.mftInvQPt(), row.deltaX(), row.deltaY(), row.deltaPhi(), row.deltaPhiWrapped(), row.deltaTanl(), row.deltaInvQPt(), row.deltaEta(), row.deltaPt(), row.chi2All(), row.chi2XYPhiTanl(), row.chi2XY(), row.isFinite(), mass, valid, isJpsi, partner, numPartners);
      if (isJpsi) {
        jpsiResiduals(row.mchTrackId(), row.mftTrackId(), row.mchCollisionId(), row.mftCollisionId(), row.matchingPlaneZ(), row.bzAtMftCenter(), row.mchX(), row.mchY(), row.mchPhi(), row.mchTanl(), row.mchInvQPt(), row.mftX(), row.mftY(), row.mftPhi(), row.mftTanl(), row.mftInvQPt(), row.deltaX(), row.deltaY(), row.deltaPhi(), row.deltaPhiWrapped(), row.deltaTanl(), row.deltaInvQPt(), row.deltaEta(), row.deltaPt(), row.chi2All(), row.chi2XYPhiTanl(), row.chi2XY(), row.isFinite(), mass, valid, isJpsi, partner, numPartners);
      }
    }
  }
};

struct MixedMatchPlanePublisher {
  o2::framework::Produces<aod::MixedMatchPlaneResiduals> residuals;
  MixedPairSampler sampler;
  o2::framework::OutputObj<TNamed> samplingMetadata{TNamed("MixedPairSampling", "")};
  void init(o2::framework::InitContext&) { samplingMetadata->SetTitle(sampler.metadata().c_str()); }
  void process(aod::RawMixedMatchPlaneResiduals const& rows)
  {
    sampler.decisions.clear();
    for (const auto& row : rows) {
      if (sampler.keep(row)) {
        copyMixedResidual(row, residuals);
      }
    }
    samplingMetadata->SetTitle(sampler.metadata().c_str());
  }
};

// These tables remain available to downstream MC devices in every mode. Only
// disk persistence changes; reconstruction, populations and failure rows do not.
inline std::string matchingOutputMode(o2::framework::ConfigContext const& context)
{
  std::string mode = context.options().get<std::string>("matching-output-mode");
  std::string configuredMode = "tree";
  mixedPairSaveFraction = 1.;
  mixedPairSamplingSeed = 0;
  matchingTreeStorage = {};
  for (int i = 1; i < context.argc(); ++i) {
    std::string arg = context.argv()[i], source;
    if (arg == "--configuration" && i + 1 < context.argc()) {
      source = context.argv()[++i];
    } else if (arg.starts_with("--configuration=")) {
      source = arg.substr(std::string("--configuration=").size());
    }
    if (source.starts_with("json://")) {
      boost::property_tree::ptree config;
      boost::property_tree::read_json(source.substr(7), config);
      configuredMode = config.get<std::string>("matching-output.cfgOutputMode", "tree");
      mixedPairSaveFraction = config.get<double>("matching-output.cfgMixedPairSaveFraction", 1.);
      mixedPairSamplingSeed = config.get<uint64_t>("matching-output.cfgMixedPairSamplingSeed", 0);
      matchingTreeStorage.mcTruth = config.get<bool>("matching-output.cfgSaveMcTruthTables", true);
      matchingTreeStorage.candidateTruth = config.get<bool>("matching-output.cfgSaveCandidateTruthTable", true);
      matchingTreeStorage.candidateRelations = config.get<bool>("matching-output.cfgSaveCandidateRelationsTable", true);
      matchingTreeStorage.mchTruth = config.get<bool>("matching-output.cfgSaveMchTruthTable", true);
      matchingTreeStorage.mixedTruth = config.get<bool>("matching-output.cfgSaveMixedTruthTable", true);
      matchingTreeStorage.mchPopulation = config.get<bool>("matching-output.cfgSaveMchPopulationTable", true);
    } else if (!source.empty() && mode == "auto") {
      LOGF(fatal, "Storage auto mode requires local json:// configuration; use --matching-output-mode for other providers");
    }
  }
  if (mode == "auto") {
    mode = configuredMode;
  }
  if (mode != "tree" && mode != "sparse" && mode != "both") {
    LOGF(fatal, "Matching output mode must be tree, sparse or both");
  }
  if (!std::isfinite(mixedPairSaveFraction) || mixedPairSaveFraction < 0. || mixedPairSaveFraction > 1.) {
    LOGF(fatal, "cfgMixedPairSaveFraction must be finite and between 0 and 1");
  }
  return mode;
}

inline void configureMatchingWriter(o2::framework::ConfigContext const&, o2::framework::WorkflowSpec const&, std::string const& mode)
{
  if (mode == "sparse") {
    LOGF(info, "Sparse matching storage: public audit tables are empty; internal staged transport is not written as AOD");
  }
}

inline std::vector<double> matchingUniformEdges(int bins, double low, double high)
{
  std::vector<double> edges(bins + 1);
  for (int i = 0; i <= bins; ++i) {
    edges[i] = low + (high - low) * i / bins;
  }
  return edges;
}

// MC-only categories. Source codes v2 put eta prime immediately after eta.
enum class MuonSignalOrigin : int8_t {
  Other = 0,
  Eta = 1,
  EtaPrime = 2,
  Rho = 3,
  Omega = 4,
  Phi = 5,
  Jpsi = 6,
  Psi2S = 7,
  Upsilon1S = 8,
  Upsilon2S = 9,
  Upsilon3S = 10,
  Charm = 11,
  Beauty = 12,
  BeautyCharm = 13,
  Unknown = 14
};
inline constexpr std::array<const char*, 15> SignalOriginLabels{
  "Other", "eta", "eta prime", "rho", "omega", "phi", "J/psi", "psi(2S)",
  "Upsilon(1S)", "Upsilon(2S)", "Upsilon(3S)", "Charm ancestry",
  "Beauty ancestry", "Beauty -> Charm ancestry", "Unknown"};
enum class MuonHfOrigin : int8_t { None = 0,
                                   Charm = 1,
                                   Beauty = 2,
                                   BeautyCharm = 3,
                                   Unknown = 4 };
inline constexpr std::array<const char*, 5> HfOriginLabels{
  "N: no HF ancestry", "C: Charm, no Beauty", "B: Beauty, no open-Charm cascade",
  "BC: Beauty -> open Charm -> mu", "U: Unknown"};
inline constexpr std::array<const char*, 15> HfPairTypeLabels{
  "N+N (no HF + no HF)", "N+C", "C+C", "N+B", "C+B", "B+B",
  "N+BC", "C+BC", "B+BC", "BC+BC", "N+U", "C+U", "B+U", "BC+U", "U+U"};
inline int hfPairType(int a, int b)
{
  if (a < 0 || a > 4) {
    a = 4;
  }
  if (b < 0 || b > 4) {
    b = 4;
  }
  const int low = std::min(a, b), high = std::max(a, b);
  return high * (high + 1) / 2 + low;
}
inline o2::framework::AxisSpec hfOriginAxis()
{
  return {5, -.5, 4.5, "hfOrigin: 0=N, 1=C, 2=B, 3=BC, 4=U"};
}
inline void labelHfOriginAxis(TAxis* axis)
{
  for (size_t i = 0; i < HfOriginLabels.size(); ++i) {
    axis->SetBinLabel(i + 1, HfOriginLabels[i]);
  }
}
inline o2::framework::AxisSpec signalOriginAxis()
{
  return {static_cast<int>(SignalOriginLabels.size()), -.5, static_cast<double>(SignalOriginLabels.size()) - .5,
          "signalOrigin: nearest target ancestor of MCH muon"};
}
inline void labelSignalOriginAxis(TAxis* axis)
{
  for (size_t i = 0; i < SignalOriginLabels.size(); ++i) {
    axis->SetBinLabel(i + 1, SignalOriginLabels[i]);
  }
}

struct MatchingSparseAxes : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::vector<double>> cfgSinglePtEdges{"cfgSinglePtEdges", {0., .2, .4, .6, .8, 1., 1.5, 2., 4., 6., 10., 100.}, "Reconstructed single pT and legacy matching-diagnostic pT edges, GeV/c"};
  o2::framework::Configurable<std::vector<double>> cfgMchPtPVEdges{"cfgMchPtPVEdges", {0., .2, .4, .6, .8, 1., 1.5, 2., 4., 6., 10., 100.}, "Common MCH pT@PV edges for MCHTrack and GlobalMuonTrack, GeV/c"};
  o2::framework::Configurable<std::vector<double>> cfgSingleEtaEdges{"cfgSingleEtaEdges", matchingUniformEdges(44, -4.5, -2.3), "Single eta edges"};
  o2::framework::Configurable<std::vector<double>> cfgSinglePhiEdges{"cfgSinglePhiEdges", matchingUniformEdges(128, 0., o2::constants::math::TwoPI), "MCHTrack/GlobalMuonTrack momentum phi edges in rad, [0,2pi)"};
  o2::framework::Configurable<std::vector<double>> cfgScoreEdges{"cfgScoreEdges", matchingUniformEdges(200, 0., 1.), "ML good/Same score edges, not a calibrated correct probability"};
  o2::framework::Configurable<std::vector<double>> cfgChi2Edges{"cfgChi2Edges", {0., .1, .2, .5, 1., 2., 3., 5., 10., 20., 40., 60., 100., 200., 500., 1000., 1.e6}, "Matching chi2 score edges"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaXYEdges{"cfgSingleDcaXYEdges", matchingUniformEdges(100, 0., 1.), "Single DCAxy magnitude edges, cm"};
  o2::framework::Configurable<std::vector<double>> cfgSingleDcaZEdges{"cfgSingleDcaZEdges", matchingUniformEdges(200, -2., 2.), "Single signed DCAz edges, cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonMassEdges{"cfgDimuonMassEdges", matchingUniformEdges(2400, 0., 12.), "Dimuon mass edges, GeV/c^2"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPtEdges{"cfgDimuonPtEdges", {0., .5, 1., 2., 3., 4., 6., 10., 100.}, "Reconstructed dimuon pT edges, GeV/c"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonRapidityEdges{"cfgDimuonRapidityEdges", matchingUniformEdges(80, -5., -1.), "Dimuon rapidity edges"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPhiEdges{"cfgDimuonPhiEdges", matchingUniformEdges(128, 0., o2::constants::math::TwoPI), "Dimuon phi edges in rad, [0,2pi)"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonDcaXYEdges{"cfgDimuonDcaXYEdges", matchingUniformEdges(200, 0., 2.), "Quadrature sum of single DCAxy edges, cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPairDcaEdges{"cfgDimuonPairDcaEdges", matchingUniformEdges(100, 0., 2.), "Two-track DCA edges, cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaXYEdges{"cfgDimuonPcaXYEdges", matchingUniformEdges(100, 0., 2.), "PCA-PV transverse distance edges, cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonPcaZEdges{"cfgDimuonPcaZEdges", matchingUniformEdges(200, -20., 20.), "PCA-PV signed z edges, cm"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonLxyEdges{"cfgDimuonLxyEdges", matchingUniformEdges(200, -2., 2.), "Signed Lxy edges, cm"};

  static o2::framework::AxisSpec axis(std::vector<double> const& edges, char const* title)
  {
    if (edges.size() < 2 || edges.size() > 1000001 || !std::all_of(edges.begin(), edges.end(), [](double x) { return std::isfinite(x); }) ||
        std::adjacent_find(edges.begin(), edges.end(), std::greater_equal<double>{}) != edges.end()) {
      LOGF(fatal, "Sparse axis edges must be finite, strictly increasing, with at least two edges: %s", title);
    }
    return {edges, title};
  }
  void validate() const
  {
    axis(cfgSinglePtEdges.value, "single pT");
    axis(cfgMchPtPVEdges.value, "MCH pT at PV");
    axis(cfgSingleEtaEdges.value, "eta");
    axis(cfgSinglePhiEdges.value, "phi");
    axis(cfgScoreEdges.value, "score");
    axis(cfgChi2Edges.value, "chi2");
    axis(cfgSingleDcaXYEdges.value, "single DCAxy");
    axis(cfgSingleDcaZEdges.value, "single DCAz");
    axis(cfgDimuonMassEdges.value, "mass");
    axis(cfgDimuonPtEdges.value, "dimuon pT");
    axis(cfgDimuonRapidityEdges.value, "rapidity");
    axis(cfgDimuonPhiEdges.value, "dimuon phi");
    axis(cfgDimuonDcaXYEdges.value, "dimuon DCAxy");
    axis(cfgDimuonPairDcaEdges.value, "pair DCA");
    axis(cfgDimuonPcaXYEdges.value, "PCA xy");
    axis(cfgDimuonPcaZEdges.value, "PCA z");
    axis(cfgDimuonLxyEdges.value, "Lxy");
  }
  // NaNs must not silently become overflow. Append a labelled missing-value bin.
  static o2::framework::AxisSpec nullable(std::vector<double> edges, char const* title)
  {
    axis(edges, title);
    edges.push_back(edges.back() + (edges.back() - edges[edges.size() - 2]));
    auto a = axis(edges, title);
    return a;
  }
  static double value(double x, std::vector<double> const& edges)
  {
    // Real overflows use the ROOT overflow bin, not the missing bin.
    if (!std::isfinite(x)) {
      return edges.back() + .5 * (edges.back() - edges[edges.size() - 2]);
    }
    return x >= edges.back() ? edges.back() + 2. * (edges.back() - edges[edges.size() - 2]) : x;
  }
};

template <bool Sparse, bool InternalPopulation = false, bool Mc = false, uint8_t TruthMask = 0>
struct MatchingOutputSparse {
  template <typename T>
  using Input = MatchingAuditInput<T, Sparse, (InternalPopulation ? 1 : 0) | TruthMask>;
  o2::framework::Configurable<std::string> cfgOutputMode{"cfgOutputMode", "tree", "tree, sparse or both; selected at workflow topology construction"};
  o2::framework::Configurable<double> cfgMixedPairSaveFraction{"cfgMixedPairSaveFraction", 1., "Random final Mixed pair retention probability, [0,1]; both planes and MC labels share a draw"};
  o2::framework::Configurable<uint64_t> cfgMixedPairSamplingSeed{"cfgMixedPairSamplingSeed", 0, "0 uses fresh random entropy; nonzero makes thinning reproducible"};
  o2::framework::Configurable<bool> cfgSaveMcTruthTables{"cfgSaveMcTruthTables", true, "Master persistence switch for candidate truth/relations, MCH truth and Mixed truth; calculations are unchanged"};
  o2::framework::Configurable<bool> cfgSaveCandidateTruthTable{"cfgSaveCandidateTruthTable", true, "Persist O2fwdcandtruth rows; requires cfgSaveMcTruthTables"};
  o2::framework::Configurable<bool> cfgSaveCandidateRelationsTable{"cfgSaveCandidateRelationsTable", true, "Persist O2fwdcandrel rows; requires cfgSaveMcTruthTables"};
  o2::framework::Configurable<bool> cfgSaveMchTruthTable{"cfgSaveMchTruthTable", true, "Persist O2fwdmchtruth rows; requires cfgSaveMcTruthTables"};
  o2::framework::Configurable<bool> cfgSaveMixedTruthTable{"cfgSaveMixedTruthTable", true, "Persist O2fwdmpmixlbl rows; requires cfgSaveMcTruthTables; Mixed MC eligibility is still checked"};
  o2::framework::Configurable<bool> cfgSaveMchPopulationTable{"cfgSaveMchPopulationTable", true, "Persist reconstructed O2fwdmchpop rows; independent of MC truth persistence"};
  MatchingSparseAxes axes;
  o2::framework::HistogramRegistry histos{"MatchingOutput", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  std::array<std::shared_ptr<THnSparse>, 3> singles{};
  std::shared_ptr<THnSparse> population;
  o2::framework::Configurable<std::vector<double>> cfgMftCandidateCountEdges{"cfgMftCandidateCountEdges", matchingUniformEdges(1001, -.5, 1000.5), "Event preselected MFT multiplicity edges"};
  void init(o2::framework::InitContext&)
  {
    using namespace o2::framework;
    axes.validate();
    std::vector<AxisSpec> mchAxes{axes.nullable(axes.cfgMchPtPVEdges, "MCH pT@PV (GeV/c)"),
                                  axes.nullable(axes.cfgSingleEtaEdges, "MCH eta@PV"),
                                  axes.nullable(axes.cfgSinglePhiEdges, "MCH phi@PV (rad)"),
                                  {3, -1.5, 1.5, "charge (e)"},
                                  axes.axis(cfgMftCandidateCountEdges, "Number of preselected MFT tracks in event")};
    if constexpr (Mc) {
      mchAxes.push_back({2, -.5, 1.5, "pairable: 0=no correct MFT candidate, 1=correct MFT candidate exists"});
      mchAxes.push_back({2, -.5, 1.5, "isDecayMuonMch: 0=no pi/K decay tag, 1=pi/K decay muon"});
      mchAxes.push_back(signalOriginAxis());
      mchAxes.push_back(hfOriginAxis());
    }
    population = histos.add<THnSparse>("MCHTrack", "Pre-matching MCH+MID tracks; one fill per eligible MCH, including no candidates", HistType::kTHnSparseD, mchAxes, true);
    if constexpr (Mc) {
      population->GetAxis(5)->SetBinLabel(1, "Not pairable");
      population->GetAxis(5)->SetBinLabel(2, "Pairable");
      population->GetAxis(6)->SetBinLabel(1, "No pi/K decay tag");
      population->GetAxis(6)->SetBinLabel(2, "Decay muon MCH");
      labelSignalOriginAxis(population->GetAxis(7));
      labelHfOriginAxis(population->GetAxis(8));
    }
    for (int m = 1; m <= 2; ++m) {
      auto const& score = m == 1 ? axes.cfgChi2Edges.value : axes.cfgScoreEdges.value;
      std::vector<AxisSpec> singleAxes{axes.nullable(axes.cfgSinglePtEdges, "GlobalMuon pT (GeV/c)"),
                                       axes.nullable(axes.cfgSingleEtaEdges, "GlobalMuon eta (MFT@DCA direction)"),
                                       axes.nullable(axes.cfgSinglePhiEdges, "GlobalMuon phi (rad; MFT@DCA direction)"),
                                       {3, -1.5, 1.5, "charge (e)"},
                                       axes.nullable(score, m == 1 ? "Matching chi2 (smaller is better)" : "ML score (Same/correct-like; larger is better)"),
                                       axes.nullable(axes.cfgSingleDcaXYEdges, "GlobalMuon DCAxy magnitude (cm)"),
                                       axes.nullable(axes.cfgSingleDcaZEdges, "GlobalMuon signed DCAz (cm)"),
                                       axes.nullable(axes.cfgMchPtPVEdges, "Original MCH pT@PV (GeV/c)")};
      if constexpr (Mc) {
        singleAxes.push_back({2, -.5, 1.5, "matchingLabel: 0=wrong/unknown, 1=correct"});
        singleAxes.push_back({2, -.5, 1.5, "isDecayMuonMch: 0=no pi/K decay tag, 1=pi/K decay muon"});
        singleAxes.push_back(signalOriginAxis());
        singleAxes.push_back(hfOriginAxis());
      }
      singles[m] = histos.add<THnSparse>(m == 1 ? "GlobalMuonTrackChi2" : "GlobalMuonTrackML", "MFT+MCH+MID; one-to-one assignment before score and DCA cuts", HistType::kTHnSparseD, singleAxes, true);
      if constexpr (Mc) {
        singles[m]->GetAxis(8)->SetBinLabel(1, "Wrong/unknown");
        singles[m]->GetAxis(8)->SetBinLabel(2, "Correct");
        singles[m]->GetAxis(9)->SetBinLabel(1, "No pi/K decay tag");
        singles[m]->GetAxis(9)->SetBinLabel(2, "Decay muon MCH");
        labelSignalOriginAxis(singles[m]->GetAxis(10));
        labelHfOriginAxis(singles[m]->GetAxis(11));
      }
    }
  }
  void processData(Input<o2::aod::MatchingMchPopulation> const& pops, Input<o2::aod::MatchingMchOutcomes> const& outs,
                   o2::aod::AssignedMatchingMuons const& assigned, Input<o2::aod::SelectedMatchingMuons> const& muons, Input<o2::aod::MatchingDimuonVertices> const& vtx)
  {
    fill(pops, outs, assigned, muons, vtx, {}, {});
  }
  PROCESS_SWITCH(MatchingOutputSparse, processData, "Reconstructed sparse output without MC truth axes", !Mc);
  void processMc(Input<o2::aod::MatchingMchPopulation> const& pops, Input<o2::aod::MatchingMchOutcomes> const& outs,
                 o2::aod::AssignedMatchingMuons const& assigned, Input<o2::aod::SelectedMatchingMuons> const& muons, Input<o2::aod::MatchingDimuonVertices> const& vertices,
                 Input<o2::aod::MatchingCandidateTruths> const& candidates, Input<o2::aod::MatchingMchTruths> const& truths)
  {
    std::map<std::pair<int64_t, int64_t>, int> labels;
    std::map<int64_t, std::array<int, 4>> origins;
    for (auto const& t : candidates) {
      labels.emplace(std::pair{t.mchTrackId(), t.mftTrackId()}, t.matchTruthStatus() == 1 ? 1 : 0);
    }
    for (auto const& t : truths) {
      origins.emplace(t.mchTrackId(), std::array<int, 4>{t.numTruePairs() > 0 ? 1 : 0, t.isDecayMuonMch() ? 1 : 0, t.signalOrigin(), t.hfOrigin()});
    }
    fill(pops, outs, assigned, muons, vertices, labels, origins);
  }
  PROCESS_SWITCH(MatchingOutputSparse, processMc, "MC MCHTrack/GlobalMuonTrack output with pairable, matchingLabel and decay flags", Mc);
  void fill(Input<o2::aod::MatchingMchPopulation> const& pops, Input<o2::aod::MatchingMchOutcomes> const& outs,
            o2::aod::AssignedMatchingMuons const& assigned, Input<o2::aod::SelectedMatchingMuons> const& muons, Input<o2::aod::MatchingDimuonVertices> const& vtx,
            std::map<std::pair<int64_t, int64_t>, int> const& matchingLabels, std::map<int64_t, std::array<int, 4>> const& muonOrigins)
  {
    (void)outs;
    std::map<int64_t, double> mchPt;
    std::set<int64_t> eligible;
    for (auto const& r : pops) {
      if (!r.isEligible() || !r.isMchMid()) {
        continue;
      }
      mchPt[r.mchTrackId()] = r.mchPtPV();
      eligible.insert(r.mchTrackId());
      std::array<double, 9> x{axes.value(r.mchPtPV(), axes.cfgMchPtPVEdges), axes.value(r.mchEtaPV(), axes.cfgSingleEtaEdges),
                              axes.value(TVector2::Phi_0_2pi(r.mchPhiPV()), axes.cfgSinglePhiEdges), static_cast<double>(r.charge()), static_cast<double>(r.numPreselectedMfts())};
      if constexpr (Mc) {
        const auto it = muonOrigins.find(r.mchTrackId());
        if (it == muonOrigins.end()) {
          LOGF(fatal, "Missing MCH truth for denominator");
        }
        x[5] = it->second[0];
        x[6] = it->second[1];
        x[7] = it->second[2];
        x[8] = it->second[3];
      }
      population->Fill(x.data());
    }
    for (auto const& r : assigned) {
      if (r.method() < 1 || r.method() > 2 || !eligible.contains(r.mchTrackId())) {
        continue;
      }
      auto const& score = r.method() == 1 ? axes.cfgChi2Edges.value : axes.cfgScoreEdges.value;
      const double pt = std::hypot(r.px(), r.py());
      std::array<double, 12> x{axes.value(pt, axes.cfgSinglePtEdges), axes.value(pt > 0. ? std::asinh(r.pz() / pt) : NAN, axes.cfgSingleEtaEdges),
                               axes.value(TVector2::Phi_0_2pi(std::atan2(r.py(), r.px())), axes.cfgSinglePhiEdges), static_cast<double>(r.charge()), axes.value(r.score(), score),
                               axes.value(std::hypot(r.dcaX(), r.dcaY()), axes.cfgSingleDcaXYEdges), axes.value(r.dcaZ(), axes.cfgSingleDcaZEdges), axes.value(mchPt.at(r.mchTrackId()), axes.cfgMchPtPVEdges)};
      if constexpr (Mc) {
        auto const it = matchingLabels.find({r.mchTrackId(), r.mftTrackId()});
        if (it == matchingLabels.end()) {
          LOGF(fatal, "Missing assigned candidate truth");
        }
        x[8] = it->second;
        x[9] = muonOrigins.at(r.mchTrackId())[1];
        x[10] = muonOrigins.at(r.mchTrackId())[2];
        x[11] = muonOrigins.at(r.mchTrackId())[3];
      }
      singles[r.method()]->Fill(x.data());
    }
    (void)muons;
    (void)vtx;
  }
};

template <bool Mc>
void addMchDimuonOutput(o2::framework::WorkflowSpec&, o2::framework::ConfigContext const&, std::string const&);
template <bool Mc>
void addGlobalDimuonOutput(o2::framework::WorkflowSpec&, o2::framework::ConfigContext const&, std::string const&);

template <bool Mc = false>
inline o2::framework::WorkflowSpec makeMatchingWorkflow(o2::framework::ConfigContext const& cfgc, bool modernDefaults = false)
{
  const auto mode = matchingOutputMode(cfgc);
  sparseAuditStorage = mode == "sparse";
  o2::framework::WorkflowSpec workflow{o2::framework::adaptAnalysisTask<MftMchMatchingTask>(cfgc, o2::framework::TaskName{"match_mft_mch_data_mc"},
                                                                                            modernDefaults ? o2::framework::SetDefaultProcesses{{{"processLegacy", false}, {"processMatchingPlane", true}}} : o2::framework::SetDefaultProcesses{}),
                                       o2::framework::adaptAnalysisTask<MchDimuonMassAtPv>(cfgc, o2::framework::TaskName{"mch-dimuon-mass-at-pv"}),
                                       o2::framework::adaptAnalysisTask<MlMatchedDimuonMassAtPv>(cfgc, o2::framework::TaskName{"ml-matched-dimuon-mass-at-pv"}),
                                       (sparseAuditStorage ? o2::framework::adaptAnalysisTask<MchJpsiMatchFlags<true>>(cfgc, o2::framework::TaskName{"mch-jpsi-match-flags"}) : o2::framework::adaptAnalysisTask<MchJpsiMatchFlags<false>>(cfgc, o2::framework::TaskName{"mch-jpsi-match-flags"}))};
  if (!modernDefaults) {
    workflow.push_back(o2::framework::adaptAnalysisTask<MixedMatchPlanePublisher>(cfgc, o2::framework::TaskName{"mixed-match-plane-publisher"}));
  }
  if (mode == "sparse") {
    workflow.push_back(o2::framework::adaptAnalysisTask<MatchingOutputSparse<true, false, Mc>>(cfgc, o2::framework::TaskName{"matching-output"}));
  } else if (mode == "both") {
    workflow.push_back(withMatchingAuditStorageMask(matchingAuditStorageMask(), [&]<uint8_t Mask>() {
      return o2::framework::adaptAnalysisTask<MatchingOutputSparse<false, false, Mc, Mask>>(cfgc, o2::framework::TaskName{"matching-output"});
    }));
  }
  addMchDimuonOutput<Mc>(workflow, cfgc, mode);
  addGlobalDimuonOutput<Mc>(workflow, cfgc, mode);
  if (!modernDefaults) {
    configureMatchingWriter(cfgc, workflow, mode);
  }
  return workflow;
}

} // namespace o2::analysis::mftmch

// -----------------------------------------------------------------------------
// MftMchMatchAncestry
// -----------------------------------------------------------------------------
namespace o2::analysis::mftmch::truth
{
// AO2D MFT and forward-track labels both use bit 7 for a fake label.
inline bool hasFakeLabel(uint8_t mchMask, uint8_t mftMask)
{
  constexpr uint8_t FakeLabelBit = uint8_t{1} << 7;
  return ((mchMask | mftMask) & FakeLabelBit) != 0;
}
enum class MatchStatus : int8_t { Unknown = -1,
                                  False = 0,
                                  Correct = 1 };
inline MatchStatus matchingTruthStatus(bool mchValid, bool mftValid, int64_t mchId, int64_t mftId,
                                       uint8_t mchMask, uint8_t mftMask)
{
  if (hasFakeLabel(mchMask, mftMask)) {
    return MatchStatus::False;
  }
  if (!mchValid || !mftValid) {
    return MatchStatus::Unknown;
  }
  return mchId == mftId ? MatchStatus::Correct : MatchStatus::False;
}

struct Node {
  int64_t id = -1;
  int pdg = 0, process = -1;
  bool valid = false, generator = false, incoming = false;
  float vx = 0.f, vy = 0.f, vz = 0.f;
  std::vector<int64_t> mothers;
};
// Completeness refers to the genealogy stored in this AO2D, not omitted MC tracks.
enum SearchStatus : int { Complete = 0,
                          MissingMother = 1,
                          DepthLimit = 2,
                          Cycle = 4,
                          NodeLimit = 8,
                          MissingLabel = 16 };
struct Path {
  std::vector<Node> childToAncestor;
  bool decayOnly = false;
};
struct Search {
  std::map<int64_t, Path> ancestors;
  int status = Complete;
};
// Search all mother branches. A generator genealogy is accepted as decay-like,
// while every transport-produced child must explicitly have the decay process.
// At least one PDG change is required: same-PDG bookkeeping copies are not decays.
template <typename Getter>
Search findAncestors(const Node& descendant, Getter getNode, int decayProcess, int maxDepth, int maxNodes)
{
  Search out;
  struct State {
    std::vector<Node> path;
    bool clean = true, changed = false;
  };
  std::deque<State> queue;
  queue.push_back({{descendant}, true, false});
  std::set<std::tuple<int64_t, bool, bool>> visited;
  visited.emplace(descendant.id, true, false);
  int expanded = 0;
  while (!queue.empty()) {
    State state = std::move(queue.front());
    queue.pop_front();
    const Node child = state.path.back();
    bool hasMother = std::any_of(child.mothers.begin(), child.mothers.end(), [](int64_t id) { return id != -1; });
    if (static_cast<int>(state.path.size()) - 1 >= maxDepth) {
      if (hasMother) {
        out.status |= DepthLimit;
      }
      continue;
    }
    if (++expanded > maxNodes) {
      out.status |= NodeLimit;
      break;
    }
    auto ids = child.mothers;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    for (int64_t const& id : ids) {
      if (id == -1) {
        continue;
      }
      if (id < 0) {
        out.status |= MissingMother;
        continue;
      }
      if (std::any_of(state.path.begin(), state.path.end(), [&](const Node& n) { return n.id == id; })) {
        out.status |= Cycle;
        continue;
      }
      Node parent = getNode(id);
      if (!parent.valid) {
        out.status |= MissingMother;
        continue;
      }
      const bool clean = state.clean && (child.generator || child.process == decayProcess);
      const bool changed = state.changed || parent.pdg != child.pdg;
      State next = state;
      next.path.push_back(parent);
      next.clean = clean;
      next.changed = changed;
      const bool decayOnly = clean && changed;
      auto found = out.ancestors.find(id);
      if (found == out.ancestors.end() || (decayOnly && !found->second.decayOnly) ||
          (decayOnly == found->second.decayOnly && next.path.size() < found->second.childToAncestor.size())) {
        out.ancestors[id] = Path{next.path, decayOnly};
      }
      if (visited.emplace(id, clean, changed).second) {
        queue.push_back(std::move(next));
      }
    }
  }
  return out;
}
// Open-heavy-flavour mesons/baryons, excluding quarkonia, diquarks and nuclei.
inline bool isOpenHeavyHadron(int pdg, int flavour)
{
  const int code = std::abs(pdg);
  if (code < 100 || code >= 1000000000) {
    return false;
  }
  const int q1 = (code / 1000) % 10, q2 = (code / 100) % 10, q3 = (code / 10) % 10;
  if (q2 < 1 || q2 > 6 || q3 < 1 || q3 > 6) {
    return false;
  }
  if (q1 == 0) {
    return (q2 == flavour) != (q3 == flavour);
  }
  return q1 <= 6 && (q1 == flavour || q2 == flavour || q3 == flavour);
}
inline MuonSignalOrigin resonanceOrigin(int pdg)
{
  switch (std::abs(pdg)) {
    case 221:
      return MuonSignalOrigin::Eta;
    case 331:
      return MuonSignalOrigin::EtaPrime;
    case 113:
      return MuonSignalOrigin::Rho;
    case 223:
      return MuonSignalOrigin::Omega;
    case 333:
      return MuonSignalOrigin::Phi;
    case 443:
      return MuonSignalOrigin::Jpsi;
    case 100443:
      return MuonSignalOrigin::Psi2S;
    case 553:
      return MuonSignalOrigin::Upsilon1S;
    case 100553:
      return MuonSignalOrigin::Upsilon2S;
    case 200553:
      return MuonSignalOrigin::Upsilon3S;
    default:
      return MuonSignalOrigin::Other;
  }
}
// Breadth-first mother search: nearest target wins, equal-depth ties use MC ID.
// All target ancestors, including HF hadrons, may precede pi/K decay.
template <typename Getter>
MuonSignalOrigin classifySignalOrigin(const Node& muon, Getter getNode, int maxDepth, int maxNodes)
{
  if (!muon.valid || std::abs(muon.pdg) != 13 || maxDepth < 1 || maxNodes < 1) {
    return MuonSignalOrigin::Unknown;
  }
  struct State {
    Node node;
    std::vector<int64_t> path;
  };
  std::vector<State> frontier{{muon, {muon.id}}};
  std::set<int64_t> seen{muon.id};
  int expanded = 0;
  bool incomplete = false;
  for (int depth = 0; !frontier.empty(); ++depth) {
    std::sort(frontier.begin(), frontier.end(), [](const auto& a, const auto& b) { return a.node.id < b.node.id; });
    for (auto const& state : frontier) {
      if (++expanded > maxNodes) {
        return MuonSignalOrigin::Unknown;
      }
      if (depth == 0) {
        continue;
      }
      auto origin = resonanceOrigin(state.node.pdg);
      const bool beauty = isOpenHeavyHadron(state.node.pdg, 5);
      const bool charm = isOpenHeavyHadron(state.node.pdg, 4);
      if (origin == MuonSignalOrigin::Other && !beauty && !charm) {
        continue;
      }
      if (incomplete) {
        return MuonSignalOrigin::Unknown;
      }
      if (origin != MuonSignalOrigin::Other) {
        return origin;
      }
      if (beauty) {
        return MuonSignalOrigin::Beauty;
      }
      // Distinguish charm without beauty from b->c on the selected charm branch.
      auto ancestors = findAncestors(state.node, getNode, -1, maxDepth - depth, maxNodes - expanded);
      for (auto const& [id, path] : ancestors.ancestors) {
        if (isOpenHeavyHadron(getNode(id).pdg, 5)) {
          return MuonSignalOrigin::BeautyCharm;
        }
      }
      return ancestors.status == Complete ? MuonSignalOrigin::Charm : MuonSignalOrigin::Unknown;
    }
    std::vector<State> next;
    for (auto const& state : frontier) {
      for (auto id : state.node.mothers) {
        if (id == -1) {
          continue;
        }
        if (depth >= maxDepth || id < 0 || std::find(state.path.begin(), state.path.end(), id) != state.path.end()) {
          incomplete = true;
          continue;
        }
        auto parent = getNode(id);
        if (!parent.valid) {
          incomplete = true;
          continue;
        }
        if (!seen.emplace(id).second) {
          continue;
        }
        auto path = state.path;
        path.push_back(id);
        next.push_back({parent, std::move(path)});
      }
    }
    frontier = std::move(next);
  }
  return incomplete ? MuonSignalOrigin::Unknown : MuonSignalOrigin::Other;
}

// Full muon ancestry, independent of signalOrigin. BC requires an ordered
// open-Charm -> Beauty path on one branch, also through pi/K or resonances.
inline MuonHfOrigin classifyHfOrigin(const Node& muon, const Search& ancestors)
{
  if (!muon.valid || std::abs(muon.pdg) != 13) {
    return MuonHfOrigin::Unknown;
  }
  bool charm = false, beauty = false;
  std::map<int64_t, Node> nodes;
  nodes.emplace(muon.id, muon);
  for (auto const& [id, path] : ancestors.ancestors) {
    for (auto const& node : path.childToAncestor) {
      nodes.emplace(node.id, node);
    }
  }
  for (auto const& [id, node] : nodes) {
    beauty = beauty || isOpenHeavyHadron(node.pdg, 5);
    if (!isOpenHeavyHadron(node.pdg, 4) || isOpenHeavyHadron(node.pdg, 5)) {
      continue;
    }
    charm = true;
    std::vector<int64_t> todo = node.mothers;
    std::set<int64_t> seen;
    for (size_t i = 0; i < todo.size(); ++i) {
      if (!seen.insert(todo[i]).second) {
        continue;
      }
      auto it = nodes.find(todo[i]);
      if (it == nodes.end()) {
        continue;
      }
      if (isOpenHeavyHadron(it->second.pdg, 5)) {
        return MuonHfOrigin::BeautyCharm;
      }
      for (auto const mother : it->second.mothers) {
        todo.push_back(mother);
      }
    }
  }
  // Incomplete genealogy cannot distinguish B from an unseen BC cascade.
  if (ancestors.status != Complete) {
    return MuonHfOrigin::Unknown;
  }
  if (beauty) {
    return MuonHfOrigin::Beauty;
  }
  return charm ? MuonHfOrigin::Charm : MuonHfOrigin::None;
}

struct Relation {
  // -1=unknown, 0=no MFT->MCH-muon ancestry, 1=same particle,
  // 2=MFT ancestor -> MCH muon via decay-only path, 3=other ancestor path.
  int origin = -1, status = Complete, depth = -1;
  bool ancestor = false, decay = false;
  std::vector<Node> parentToChild;
};
inline Relation classify(const Node& mch, const Node& mft, const Search* search)
{
  Relation out;
  if (!mch.valid || !mft.valid) {
    out.status = MissingLabel;
    return out;
  }
  if (mch.id == mft.id) {
    out.origin = 1;
    out.depth = 0;
    out.parentToChild = {mch};
    return out;
  }
  if (mch.pdg != kMuonMinus && mch.pdg != -kMuonMinus) {
    out.origin = 0;
    return out;
  }
  if (!search) {
    out.status = MissingMother;
    return out;
  }
  out.status = search->status;
  auto found = search->ancestors.find(mft.id);
  if (found == search->ancestors.end()) {
    out.origin = out.status == Complete ? 0 : -1;
    return out;
  }
  out.ancestor = true;
  out.decay = found->second.decayOnly;
  out.origin = out.decay ? 2 : 3;
  out.parentToChild = found->second.childToAncestor;
  std::reverse(out.parentToChild.begin(), out.parentToChild.end());
  out.depth = static_cast<int>(out.parentToChild.size()) - 1;
  return out;
}
} // namespace o2::analysis::mftmch::truth

namespace o2::analysis::mftmch
{
struct DimuonMchTag {
  bool pairable = false, decay = false, valid = false;
  int hf = static_cast<int>(MuonHfOrigin::Unknown);
  int64_t particle = -1, collision = -1;
};

template <typename Getter>
std::pair<int, int> commonDimuonMother(DimuonMchTag const& a, DimuonMchTag const& b, Getter getNode, int maxDepth, int maxNodes)
{
  using namespace truth;
  if (!a.valid || !b.valid || a.collision < 0 || b.collision < 0) {
    return {-1, static_cast<int>(MuonSignalOrigin::Unknown)};
  }
  if (a.collision != b.collision) {
    return {0, static_cast<int>(MuonSignalOrigin::Other)};
  }
  auto ancestors = [&](int64_t child) {
    std::map<int64_t, int> hadrons;
    bool complete = true;
    std::vector<std::pair<int64_t, int>> todo;
    for (auto id : getNode(child).mothers) {
      todo.emplace_back(id, 1);
    }
    std::set<int64_t> seen;
    for (size_t i = 0; i < todo.size(); ++i) {
      const auto [id, depth] = todo[i];
      if (id == -1) {
        continue;
      }
      if (!seen.insert(id).second) {
        continue;
      }
      if (depth > maxDepth || seen.size() > static_cast<size_t>(maxNodes)) {
        complete = false;
        break;
      }
      auto node = getNode(id);
      if (!node.valid) {
        complete = false;
        continue;
      }
      const auto pdg = std::abs(node.pdg);
      // Incoming beams, nuclei and partons are not a decay mother.
      const int q1 = (pdg / 1000) % 10, q2 = (pdg / 100) % 10, q3 = (pdg / 10) % 10;
      const bool hadron = pdg >= 100 && pdg < 1000000000 && q1 <= 6 && q2 >= 1 && q2 <= 6 && q3 >= 1 && q3 <= 6;
      if (!node.incoming && hadron) {
        hadrons.emplace(id, depth);
      }
      for (auto mother : node.mothers) {
        todo.emplace_back(mother, depth + 1);
      }
    }
    return std::make_pair(hadrons, complete);
  };
  auto const [aa, completeA] = ancestors(a.particle);
  auto const [bb, completeB] = ancestors(b.particle);
  if (!completeA || !completeB) {
    return {-1, static_cast<int>(MuonSignalOrigin::Unknown)};
  }
  std::tuple<int, int, int64_t> best{std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), -1};
  for (auto const& [id, depth] : aa) {
    auto it = bb.find(id);
    if (it != bb.end()) {
      best = std::min(best, std::make_tuple(std::max(depth, it->second), depth + it->second, id));
    }
  }
  const int64_t id = std::get<2>(best);
  if (id < 0) {
    return {0, static_cast<int>(MuonSignalOrigin::Other)};
  }
  Node proxy;
  proxy.id = -1;
  proxy.pdg = 13;
  proxy.valid = true;
  proxy.mothers = {id};
  // The classifier starts at mothers, so this reports exactly the chosen
  // common hadron (and distinguishes a beauty-fed common charm ancestor).
  const auto origin = classifySignalOrigin(proxy, getNode, maxDepth, maxNodes);
  const auto node = getNode(id);
  const bool target = resonanceOrigin(node.pdg) != MuonSignalOrigin::Other || isOpenHeavyHadron(node.pdg, 4) || isOpenHeavyHadron(node.pdg, 5);
  return {1, static_cast<int>(target ? origin : MuonSignalOrigin::Other)};
}

// The pair stream is produced by exactly the canonical inference event pool.
// Mixed needs no MC labels. Same is annotated solely within its source DF.
template <bool Mc, bool Sparse, uint8_t Mask = 0>
struct MchDimuonOutput {
  template <typename T>
  using Input = MatchingAuditInput<T, Sparse, Mask>;
  MatchingSparseAxes axes;
  o2::framework::Configurable<int> cfgMcAncestryMaxDepth{"cfgMcAncestryMaxDepth", 128, "Maximum common-mother ancestor depth"};
  o2::framework::Configurable<int> cfgMcAncestryMaxNodes{"cfgMcAncestryMaxNodes", 4096, "Maximum common-mother ancestor nodes per muon"};
  o2::framework::HistogramRegistry histos{"MchDimuonOutput", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  std::array<std::shared_ptr<THnSparse>, 2> spectra{};
  using Tag = DimuonMchTag;
  void init(o2::framework::InitContext&)
  {
    axes.validate();
    if (cfgMcAncestryMaxDepth < 1 || cfgMcAncestryMaxNodes < 1) {
      LOGF(fatal, "Common-mother ancestry limits must be positive");
    }
  }
  void createSpectra()
  {
    if (spectra[0]) {
      return;
    }
    using namespace o2::framework;
    for (int mixed = 0; mixed <= 1; ++mixed) {
      std::vector<AxisSpec> aa{axes.axis(axes.cfgDimuonMassEdges, "mass (GeV/c^{2})"),
                               axes.axis(axes.cfgDimuonPtEdges, "dimuon pT (GeV/c)"),
                               axes.axis(axes.cfgDimuonRapidityEdges, "dimuon rapidity"),
                               axes.axis(axes.cfgDimuonPhiEdges, "dimuon phi (rad; [0,2pi))"),
                               {5, -2.5, 2.5, "charge sum: -2=LSmm, 0=ULS, +2=LSpp"}};
      if constexpr (Mc) {
        if (!mixed) {
          aa.push_back({3, -.5, 2.5, "nPairableMuons: correct MFT exists in candidate set"});
          aa.push_back({3, -1.5, 1.5, "isSameMother: -1=unknown, 0=no, 1=shared hadron MC ID"});
          aa.push_back(signalOriginAxis());
          aa.back().title = "motherLabel: nearest common hadron ancestor";
          aa.push_back({3, -.5, 2.5, "nDecayMuonMch: number of MCH muons with pi/K ancestry"});
          aa.push_back({15, -.5, 14.5, "hfPairType: N=no HF, C=Charm, B=Beauty, BC=Beauty-Charm, U=Unknown"});
        }
      }
      spectra[mixed] = histos.add<THnSparse>(mixed ? "MCHDimuonMixed" : "MCHDimuonSame", "MCH+MID PV dimuons; canonical event pool and weights", HistType::kTHnSparseD, aa, true);
      auto charge = spectra[mixed]->GetAxis(4);
      charge->SetBinLabel(1, "LSmm (-2)");
      charge->SetBinLabel(3, "ULS (0)");
      charge->SetBinLabel(5, "LSpp (+2)");
      if constexpr (Mc) {
        if (!mixed) {
          auto flag = spectra[mixed]->GetAxis(6);
          flag->SetBinLabel(1, "Unknown ancestry/label");
          flag->SetBinLabel(2, "No common hadron");
          flag->SetBinLabel(3, "Same hadron MC ID");
          labelSignalOriginAxis(spectra[mixed]->GetAxis(7));
          const int decayAxis = 8, hfAxis = 9;
          for (int n = 0; n <= 2; ++n) {
            spectra[mixed]->GetAxis(5)->SetBinLabel(n + 1, (std::to_string(n) + " pairable muons").c_str());
            spectra[mixed]->GetAxis(decayAxis)->SetBinLabel(n + 1, (std::to_string(n) + " pi/K decay muons").c_str());
          }
          for (size_t k = 0; k < HfPairTypeLabels.size(); ++k) {
            spectra[mixed]->GetAxis(hfAxis)->SetBinLabel(k + 1, HfPairTypeLabels[k]);
          }
        }
      }
    }
  }
  void processData(o2::aod::RawMchDimuonPairs const& pairs)
  {
    if (pairs.size() == 0) {
      return;
    }
    createSpectra();
    for (auto const& row : pairs) {
      const double x[]{row.mass(), row.pairPt(), row.pairRapidity(), row.phi(), static_cast<double>(row.pairCharge())};
      spectra[row.isMixed()]->Fill(x, row.weight());
    }
  }
  PROCESS_SWITCH(MchDimuonOutput, processData, "Data MCH dimuon sparse", !Mc);
  void processMC(o2::aod::RawMchDimuonPairs const& pairs,
                 Input<o2::aod::MatchingMchTruths> const& labels, o2::aod::McParticles const& particles)
  {
    std::map<int64_t, Tag> currentTruth;
    for (auto const& row : labels) {
      Tag tag;
      tag.pairable = row.numTruePairs() > 0;
      tag.decay = row.isDecayMuonMch();
      tag.hf = row.hfOrigin();
      tag.valid = row.hasMchMcParticle() && (row.mchMcMask() & 128) == 0 && std::abs(row.mchPdgCode()) == 13;
      tag.particle = row.mchMcParticleId();
      tag.collision = row.mchMcCollisionId();
      currentTruth.emplace(row.mchTrackId(), tag);
    }
    auto getNode = [&](int64_t id) {
      truth::Node n;
      n.id = id;
      if (id < static_cast<int64_t>(particles.offset()) || id - static_cast<int64_t>(particles.offset()) >= particles.size()) {
        return n;
      }
      auto row = particles.rawIteratorAt(id - particles.offset());
      n.valid = true;
      n.pdg = row.pdgCode();
      n.generator = row.producedByGenerator();
      n.incoming = row.getHepMCStatusCode() == 4;
      for (auto mother : row.mothersIds()) {
        n.mothers.push_back(mother);
      }
      return n;
    };
    if (pairs.size() != 0) {
      createSpectra();
    }
    std::map<std::pair<int64_t, int64_t>, std::pair<int, int>> mothers;
    for (auto const& row : pairs) {
      if (row.isMixed()) {
        const double x[]{row.mass(), row.pairPt(), row.pairRapidity(), row.phi(), static_cast<double>(row.pairCharge())};
        spectra[1]->Fill(x, row.weight());
        continue;
      }
      auto a = currentTruth.find(row.mchId1()), b = currentTruth.find(row.mchId2());
      if (a == currentTruth.end() || b == currentTruth.end()) {
        LOGF(fatal, "Missing Same MCH dimuon truth");
      }
      std::array<double, 10> x{row.mass(), row.pairPt(), row.pairRapidity(), row.phi(), static_cast<double>(row.pairCharge()), static_cast<double>(a->second.pairable + b->second.pairable)};
      const std::pair<int64_t, int64_t> key{std::min(row.mchId1(), row.mchId2()), std::max(row.mchId1(), row.mchId2())};
      auto [it, inserted] = mothers.try_emplace(key);
      if (inserted) {
        it->second = commonDimuonMother(a->second, b->second, getNode, cfgMcAncestryMaxDepth.value, cfgMcAncestryMaxNodes.value);
      }
      x[6] = it->second.first;
      x[7] = it->second.second;
      x[8] = a->second.decay + b->second.decay;
      x[9] = hfPairType(a->second.hf, b->second.hf);
      spectra[0]->Fill(x.data(), row.weight());
    }
  }
  PROCESS_SWITCH(MchDimuonOutput, processMC, "MC Same truth and inclusive Mixed MCH dimuon sparse", Mc);
};
template <bool Mc>
void addMchDimuonOutput(o2::framework::WorkflowSpec& workflow, o2::framework::ConfigContext const& context, std::string const& mode)
{
  using namespace o2::framework;
  if constexpr (!Mc) {
    workflow.push_back(adaptAnalysisTask<MchDimuonOutput<false, false>>(context, TaskName{"mch-dimuon-output"}));
  } else if (mode == "sparse") {
    workflow.push_back(adaptAnalysisTask<MchDimuonOutput<true, true>>(context, TaskName{"mch-dimuon-output"}));
  } else if (useInternalAuditStorage<o2::aod::MatchingMchTruths>()) {
    workflow.push_back(adaptAnalysisTask<MchDimuonOutput<true, false, 8>>(context, TaskName{"mch-dimuon-output"}));
  } else {
    workflow.push_back(adaptAnalysisTask<MchDimuonOutput<true, false>>(context, TaskName{"mch-dimuon-output"}));
  }
}
} // namespace o2::analysis::mftmch

namespace o2::analysis::mftmch
{
template <bool Mc, bool Sparse, uint8_t Mask = 0>
struct GlobalDimuonOutput {
  template <typename T>
  using Input = MatchingAuditInput<T, Sparse, Mask>;
  MatchingSparseAxes axes;
  o2::framework::Configurable<std::vector<double>> cfgDimuonMlScoreEdges{"cfgDimuonMlScoreEdges", [] {
                                                                           auto edges = matchingUniformEdges(200, 0., 1.);
                                                                           edges.back() = 1.000001; // Keep an exactly saturated sigmoid score of 1 in range.
                                                                           return edges;
                                                                         }(),
                                                                         "Pair min ML good/Same score edges; higher is better; not a correct probability"};
  o2::framework::Configurable<std::vector<double>> cfgDimuonChi2ScoreEdges{"cfgDimuonChi2ScoreEdges", {0., .1, .2, .5, 1., 2., 3., 5., 10., 20., 40., 60., 100., 200., 500., 1000., 1.e6}, "Pair max matching chi2 edges; lower is better"};
  o2::framework::Configurable<int> cfgMcAncestryMaxDepth{"cfgMcAncestryMaxDepth", 128, "Maximum common-mother ancestor depth"};
  o2::framework::Configurable<int> cfgMcAncestryMaxNodes{"cfgMcAncestryMaxNodes", 4096, "Maximum common-mother ancestor nodes per muon"};
  o2::framework::HistogramRegistry histos{"GlobalDimuonOutput", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, std::shared_ptr<THnSparse>> spectra;
  void init(o2::framework::InitContext&)
  {
    axes.validate();
    MatchingSparseAxes::axis(cfgDimuonMlScoreEdges.value, "dimuon ML score");
    MatchingSparseAxes::axis(cfgDimuonChi2ScoreEdges.value, "dimuon matching chi2");
    if (cfgMcAncestryMaxDepth < 1 || cfgMcAncestryMaxNodes < 1) {
      LOGF(fatal, "Common-mother ancestry limits must be positive");
    }
  }
  static std::string prefix(int method, std::string const& workingPoint)
  {
    if (method != 1 && method != 2) {
      LOGF(fatal, "Unexpected GlobalMuon matching method");
    }
    return workingPoint.empty() ? (method == 1 ? "Chi2" : "ML") : "WorkingPoints/" + workingPoint;
  }
  static std::string chargeName(int charge)
  {
    if (charge != -2 && charge != 0 && charge != 2) {
      LOGF(fatal, "Unexpected dimuon charge sum");
    }
    return charge == 0 ? "ULS" : charge > 0 ? "LSpp"
                                            : "LSmm";
  }
  void createSpectra(std::string const& method, int matchingMethod)
  {
    using namespace o2::framework;
    if (spectra.contains(method + "/Same/ULS")) {
      return;
    }
    for (bool mixed : {false, true}) {
      for (auto const& charge : {"ULS", "LSpp", "LSmm"}) {
        std::vector<AxisSpec> aa{axes.axis(axes.cfgDimuonMassEdges, "mass (GeV/c^{2})"),
                                 axes.axis(axes.cfgDimuonPtEdges, "dimuon pT (GeV/c)"),
                                 axes.axis(axes.cfgDimuonRapidityEdges, "dimuon rapidity"),
                                 axes.nullable(axes.cfgDimuonPairDcaEdges, "two-track DCA (cm)"),
                                 axes.nullable(axes.cfgDimuonPcaXYEdges, "PCA-PV transverse distance (cm)"),
                                 axes.nullable(axes.cfgDimuonPcaZEdges, "PCA z - PV z (cm)"),
                                 axes.axis(axes.cfgDimuonDcaXYEdges, "dimuon DCAxy: sqrt(DCAxy1^2+DCAxy2^2) (cm)"),
                                 axes.nullable(axes.cfgDimuonLxyEdges, "signed Lxy (cm)"),
                                 {2, -.5, 1.5, "fitSuccess: 0=failed/not attempted, 1=success"}};
        if constexpr (Mc) {
          if (!mixed) {
            aa.push_back({3, -.5, 2.5, "nPairableMuons: correct MFT exists in candidate set"});
            aa.push_back({3, -1.5, 1.5, "isSameMother: -1=unknown, 0=no, 1=shared hadron MC ID"});
            aa.push_back(signalOriginAxis());
            aa.back().title = "motherLabel: nearest common hadron ancestor";
            aa.push_back({3, -.5, 2.5, "nDecayMuonMch: number of MCH muons with pi/K ancestry"});
            aa.push_back({15, -.5, 14.5, "hfPairType: N=no HF, C=Charm, B=Beauty, BC=Beauty-Charm, U=Unknown"});
            aa.push_back({3, -.5, 2.5, "nCorrectMatches: 0=both Wrong, 1=one Correct, 2=both Correct"});
          }
        }
        aa.push_back(axes.nullable(matchingMethod == 1 ? cfgDimuonChi2ScoreEdges.value : cfgDimuonMlScoreEdges.value,
                                   matchingMethod == 1 ? "pairMatchingScore: max matching chi2" : "pairMatchingScore: min ML good score"));
        const std::string path = method + (mixed ? "/Mixed/" : "/Same/") + charge;
        auto h = histos.add<THnSparse>(path, "GlobalMuon dimuons BEFORE score, single-DCA and vertex cuts; named working points keep their explicit selection; failed fits retained", HistType::kTHnSparseD, aa, true);
        auto* scoreAxis = h->GetAxis(h->GetNdimensions() - 1);
        scoreAxis->SetBinLabel(scoreAxis->GetNbins(), "Missing pair score");
        for (int i : {3, 4, 5, 7}) {
          h->GetAxis(i)->SetBinLabel(h->GetAxis(i)->GetNbins(), "Missing/not fitted");
        }
        h->GetAxis(8)->SetBinLabel(1, "Failed/not attempted");
        h->GetAxis(8)->SetBinLabel(2, "Success");
        if constexpr (Mc) {
          if (!mixed) {
            for (int i = 0; i <= 2; ++i) {
              h->GetAxis(9)->SetBinLabel(i + 1, (std::to_string(i) + " pairable muons").c_str());
              h->GetAxis(12)->SetBinLabel(i + 1, (std::to_string(i) + " pi/K decay muons").c_str());
            }
            h->GetAxis(10)->SetBinLabel(1, "Unknown ancestry/label");
            h->GetAxis(10)->SetBinLabel(2, "No common hadron");
            h->GetAxis(10)->SetBinLabel(3, "Same hadron MC ID");
            labelSignalOriginAxis(h->GetAxis(11));
            for (size_t k = 0; k < HfPairTypeLabels.size(); ++k) {
              h->GetAxis(13)->SetBinLabel(k + 1, HfPairTypeLabels[k]);
            }
            h->GetAxis(14)->SetBinLabel(1, "Both Wrong/unknown");
            h->GetAxis(14)->SetBinLabel(2, "One Correct");
            h->GetAxis(14)->SetBinLabel(3, "Both Correct");
          }
        }
        spectra.emplace(path, std::move(h));
      }
    }
  }
  template <typename Row>
  void fill(Row const& row, std::array<int, 6> const& truth = {})
  {
    const auto nameBytes = row.workingPoint();
    const auto method = prefix(row.method(), std::string(nameBytes.begin(), nameBytes.end()));
    createSpectra(method, row.method());
    std::array<double, 16> x{row.mass(), row.pairPt(), row.pairRapidity(),
                             axes.value(row.pairDca(), axes.cfgDimuonPairDcaEdges),
                             axes.value(row.pcaPVXY(), axes.cfgDimuonPcaXYEdges),
                             axes.value(row.pcaPVZ(), axes.cfgDimuonPcaZEdges), row.dimuonDcaXY(),
                             axes.value(row.lxy(), axes.cfgDimuonLxyEdges), static_cast<double>(row.fitSuccess())};
    if constexpr (Mc) {
      if (!row.isMixed()) {
        std::copy(truth.begin(), truth.end(), x.begin() + 9);
      }
    }
    const int scoreAxis = Mc && !row.isMixed() ? 15 : 9;
    x[scoreAxis] = axes.value(row.pairMatchingScore(), row.method() == 1 ? cfgDimuonChi2ScoreEdges.value : cfgDimuonMlScoreEdges.value);
    spectra.at(method + (row.isMixed() ? "/Mixed/" : "/Same/") + chargeName(row.pairCharge()))->Fill(x.data(), row.weight());
  }
  void processData(o2::aod::RawGlobalDimuonPairs const& pairs)
  {
    for (auto const& row : pairs) {
      fill(row);
    }
  }
  PROCESS_SWITCH(GlobalDimuonOutput, processData, "Data unified GlobalMuon dimuon sparse", !Mc);
  void processMC(o2::aod::RawGlobalDimuonPairs const& pairs,
                 Input<o2::aod::MatchingMchTruths> const& labels,
                 Input<o2::aod::MatchingCandidateTruths> const& candidates, o2::aod::McParticles const& particles)
  {
    std::map<int64_t, DimuonMchTag> tags;
    for (auto const& row : labels) {
      DimuonMchTag tag;
      tag.pairable = row.numTruePairs() > 0;
      tag.decay = row.isDecayMuonMch();
      tag.hf = row.hfOrigin();
      tag.valid = row.hasMchMcParticle() && (row.mchMcMask() & 128) == 0 && std::abs(row.mchPdgCode()) == 13;
      tag.particle = row.mchMcParticleId();
      tag.collision = row.mchMcCollisionId();
      tags.emplace(row.mchTrackId(), tag);
    }
    std::map<std::pair<int64_t, int64_t>, int> correct;
    for (auto const& row : candidates) {
      correct.emplace(std::make_pair(row.mchTrackId(), row.mftTrackId()), row.matchTruthStatus() == 1 ? 1 : 0);
    }
    auto getNode = [&](int64_t id) {
      truth::Node n;
      n.id = id;
      if (id < static_cast<int64_t>(particles.offset()) || id - static_cast<int64_t>(particles.offset()) >= particles.size()) {
        return n;
      }
      auto row = particles.rawIteratorAt(id - particles.offset());
      n.valid = true;
      n.pdg = row.pdgCode();
      n.generator = row.producedByGenerator();
      n.incoming = row.getHepMCStatusCode() == 4;
      for (auto mother : row.mothersIds()) {
        n.mothers.push_back(mother);
      }
      return n;
    };
    std::map<std::pair<int64_t, int64_t>, std::pair<int, int>> mothers;
    for (auto const& row : pairs) {
      // No donor MC cache, no truth axes, and no MC lookup for Mixed.
      if (row.isMixed()) {
        fill(row);
        continue;
      }
      auto a = tags.find(row.mchId1()), b = tags.find(row.mchId2());
      auto ca = correct.find({row.mchId1(), row.mftId1()}), cb = correct.find({row.mchId2(), row.mftId2()});
      if (a == tags.end() || b == tags.end() || ca == correct.end() || cb == correct.end()) {
        LOGF(fatal, "Missing Same GlobalMuon dimuon MC annotation");
      }
      const std::pair<int64_t, int64_t> key{std::min(row.mchId1(), row.mchId2()), std::max(row.mchId1(), row.mchId2())};
      auto [it, inserted] = mothers.try_emplace(key);
      if (inserted) {
        it->second = commonDimuonMother(a->second, b->second, getNode, cfgMcAncestryMaxDepth.value, cfgMcAncestryMaxNodes.value);
      }
      fill(row, {a->second.pairable + b->second.pairable, it->second.first, it->second.second,
                 a->second.decay + b->second.decay, hfPairType(a->second.hf, b->second.hf), ca->second + cb->second});
    }
  }
  PROCESS_SWITCH(GlobalDimuonOutput, processMC, "MC Same with 6 truth axes; inclusive Mixed without MC axes", Mc);
};
template <bool Mc>
void addGlobalDimuonOutput(o2::framework::WorkflowSpec& workflow, o2::framework::ConfigContext const& context, std::string const& mode)
{
  using namespace o2::framework;
  if constexpr (!Mc) {
    workflow.push_back(adaptAnalysisTask<GlobalDimuonOutput<false, false>>(context, TaskName{"global-dimuon-output"}));
  } else if (mode == "sparse") {
    workflow.push_back(adaptAnalysisTask<GlobalDimuonOutput<true, true>>(context, TaskName{"global-dimuon-output"}));
  } else {
    const bool internalMch = useInternalAuditStorage<o2::aod::MatchingMchTruths>();
    const bool internalCandidate = useInternalAuditStorage<o2::aod::MatchingCandidateTruths>();
    if (internalMch && internalCandidate) {
      workflow.push_back(adaptAnalysisTask<GlobalDimuonOutput<true, false, 10>>(context, TaskName{"global-dimuon-output"}));
    } else if (internalMch) {
      workflow.push_back(adaptAnalysisTask<GlobalDimuonOutput<true, false, 8>>(context, TaskName{"global-dimuon-output"}));
    } else if (internalCandidate) {
      workflow.push_back(adaptAnalysisTask<GlobalDimuonOutput<true, false, 2>>(context, TaskName{"global-dimuon-output"}));
    } else {
      workflow.push_back(adaptAnalysisTask<GlobalDimuonOutput<true, false>>(context, TaskName{"global-dimuon-output"}));
    }
  }
}
} // namespace o2::analysis::mftmch

#endif // COMMON_CORE_MFTMCHMATCHING_H_
