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

#ifndef PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_
#define PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_

#include <cstdlib>
#include <cmath>
#include <array>
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#ifndef HomogeneousField
#define HomogeneousField
#endif

/// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

namespace o2
{
namespace pwglf
{
//__________________________________________
// V0 group: abstraction to deal with duplicates
// in an intuitive manner
struct V0group {
  std::vector<int> V0Ids;        // index list to original aod::V0s
  std::vector<int> collisionIds; // coll indices
  int posTrackId;
  int negTrackId;
  uint8_t v0Type;
};

//_______________________________________________________________________
template <typename T>
std::vector<std::size_t> sort_indices_posTrack(const std::vector<T>& v)
{
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](std::size_t i1, std::size_t i2) { return v[i1].posTrackId < v[i2].posTrackId; });
  return idx;
}

//_______________________________________________________________________
template <typename T>
std::vector<std::size_t> sort_indices_negTrack(const std::vector<T>& v)
{
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](std::size_t i1, std::size_t i2) { return v[i1].negTrackId < v[i2].negTrackId; });
  return idx;
}

//_______________________________________________________________________
// this function deals with the fact that V0s provided in AO2Ds may
// be duplicated in several collisions and groups them into entries
// of type pwglf::V0group, each entry having the same neg/pos tracks
// but an array of compatible collisions. The original V0 indices
// are preserved in the resulting structure to allow for easy referencing
// back afterwards. Algorithmically, full N^2 loops and/or multiple
// find calls are avoided via sorting.
template <typename T>
std::vector<V0group> groupDuplicates(const T& V0s)
{
  std::vector<V0group> v0table;
  if (V0s.size() == 0) {
    return v0table;
  }
  V0group thisV0;
  thisV0.V0Ids.push_back(-1);        // create one single element
  thisV0.collisionIds.push_back(-1); // create one single element
  for (auto const& V0 : V0s) {
    thisV0.V0Ids[0] = V0.globalIndex();
    thisV0.collisionIds[0] = V0.collisionId();
    thisV0.posTrackId = V0.posTrackId();
    thisV0.negTrackId = V0.negTrackId();
    thisV0.v0Type = V0.v0Type();
    v0table.push_back(thisV0);
  }

  // sort tracks according to positive track index to avoid excessive N^2 searches
  auto posTrackSort = sort_indices_posTrack(v0table);

  // create a proper list of V0s including duplicates: collisionIds is now a vector
  int atPosTrackId = v0table[posTrackSort[0]].posTrackId;
  std::vector<V0group> v0tableFixedPositive; // small list with fixed positive id
  std::vector<V0group> v0tableGrouped;       // final list with proper grouping
  for (size_t iV0 = 0; iV0 < posTrackSort.size(); iV0++) {
    if (atPosTrackId != v0table[posTrackSort[iV0]].posTrackId) {
      // switched pos track id. Process chunk of V0s
      auto negTrackSort = sort_indices_negTrack(v0tableFixedPositive);
      thisV0.collisionIds.clear();
      thisV0.V0Ids.clear();
      thisV0.negTrackId = v0tableFixedPositive[negTrackSort[0]].negTrackId;
      for (size_t iPV0 = 0; iPV0 < v0tableFixedPositive.size(); iPV0++) {
        if (thisV0.negTrackId != v0tableFixedPositive[negTrackSort[iPV0]].negTrackId) {
          v0tableGrouped.push_back(thisV0);
          thisV0.collisionIds.clear(); // clean collision Ids
          thisV0.V0Ids.clear();        // clean aod::V0s Ids
        }
        thisV0.V0Ids.push_back(v0tableFixedPositive[negTrackSort[iPV0]].V0Ids[0]);
        thisV0.collisionIds.push_back(v0tableFixedPositive[negTrackSort[iPV0]].collisionIds[0]);
        thisV0.posTrackId = v0tableFixedPositive[negTrackSort[iPV0]].posTrackId;
        thisV0.negTrackId = v0tableFixedPositive[negTrackSort[iPV0]].negTrackId;
        thisV0.v0Type = v0tableFixedPositive[negTrackSort[iPV0]].v0Type;
      }
      v0tableGrouped.push_back(thisV0); // publish last
      v0tableFixedPositive.clear();
      atPosTrackId = v0table[posTrackSort[iV0]].posTrackId; // move to the next pos index
    }
    v0tableFixedPositive.push_back(v0table[posTrackSort[iV0]]);
  }
  v0tableGrouped.push_back(thisV0); // publish last

  LOGF(debug, "Duplicate V0s grouped. aod::V0s counted: %i, unique index pairs: %i", V0s.size(), v0tableGrouped.size());
  return v0tableGrouped;
}

//__________________________________________
// V0 information storage
struct v0candidate {
  // indexing
  int collisionId = -1;
  int negativeTrack = -1;
  int positiveTrack = -1;

  // daughter properties
  std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> positivePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativePosition = {0.0f, 0.0f, 0.0f};
  float positiveTrackX = 0.0f;
  float negativeTrackX = 0.0f;
  float positiveDCAxy = 0.0f;
  float negativeDCAxy = 0.0f;

  // V0 properties
  std::array<float, 3> momentum = {0.0f, 0.0f, 0.0f}; // necessary for KF
  std::array<float, 3> position = {0.0f, 0.0f, 0.0f};
  float daughterDCA = 1000.0f;
  float pointingAngle = 1.0f;
  float dcaToPV = 0.0f;
  float v0DCAToPVxy = 0.0f;
  float v0DCAToPVz = 0.0f;

  // calculated masses for convenience
  float massGamma;
  float massK0Short;
  float massLambda;
  float massAntiLambda;

  // stored for decay chains
  float positionCovariance[6];
  float momentumCovariance[6];
};

//__________________________________________
// Cascade information storage
struct cascadeCandidate {
  // indexing
  int collisionId = -1;
  int v0Id = -1;
  int negativeTrack = -1;
  int positiveTrack = -1;
  int bachelorTrack = -1;

  // daughter properties
  std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> bachelorMomentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> positivePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> negativePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> bachelorPosition = {0.0f, 0.0f, 0.0f};
  float positiveDCAxy = 0.0f;
  float negativeDCAxy = 0.0f;
  float bachelorDCAxy = 0.0f;
  float positiveTrackX = 0.0f;
  float negativeTrackX = 0.0f;
  float bachelorTrackX = 0.0f;

  // cascade properties
  int charge = -1; // default: []Minus
  float cascadeDaughterDCA = 1000.0f;
  float v0DaughterDCA = 1000.0f;

  float pointingAngle = 0.0f;
  float cascadeDCAxy = 0.0f;
  float cascadeDCAz = 0.0f;
  std::array<float, 3> v0Position = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> v0Momentum = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> cascadePosition = {0.0f, 0.0f, 0.0f};
  std::array<float, 3> cascadeMomentum = {0.0f, 0.0f, 0.0f};

  float bachBaryonCosPA = 0.0f;
  float bachBaryonDCAxyToPV = 1e+3;
  float massXi = 0.0f;
  float massOmega = 0.0f;

  // KF-specific variables
  float kfV0Chi2 = 0.0f;
  float kfCascadeChi2 = 0.0f;
  float kfMLambda = 0.0f;
  float kfTrackCovarianceV0[21];
  float kfTrackCovariancePos[21];
  float kfTrackCovarianceNeg[21];

  // stored for decay chains
  float covariance[21];
};

//__________________________________________
// builder helper class
class strangenessBuilderHelper
{
 public:
  strangenessBuilderHelper()
  {
    // standards hardcoded in builder ...
    // ...but can be changed easily since fitter is public
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxDXYIni(4.0f);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);

    v0selections.minCrossedRows = -1;
    v0selections.dcanegtopv = -1.0f;
    v0selections.dcapostopv = -1.0f;
    v0selections.v0cospa = -2;
    v0selections.dcav0dau = 1e+6;
    v0selections.v0radius = 0.0f;
    v0selections.maxDaughterEta = 2.0;

    // LUT has to be loaded later
    lut = nullptr;
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

    // mag field has to be set later
    fitter.setBz(-999.9f); // will NOT make sense if not changed
  };

  //_______________________________________________________________________
  // standard build V0 function. Populates ::v0 object
  // ::v0 will be initialized to defaults if build fails
  // --- useSelections: meant to maximize recovery, but beware high cost in CPU
  // --- calculateProngDCAtoPV: optionally don't propagate prongs to PV, saves
  //     CPU, of interest when dealing with de-duplication (variable not checked)
  template <bool useSelections = true, bool calculateProngDCAtoPV = true, typename TTrack, typename TTrackParametrization>
  bool buildV0Candidate(int collisionIndex,
                        float pvX, float pvY, float pvZ,
                        TTrack const& positiveTrack,
                        TTrack const& negativeTrack,
                        TTrackParametrization& positiveTrackParam,
                        TTrackParametrization& negativeTrackParam,
                        bool useCollinearFit = false,
                        bool calculateCovariance = false,
                        bool acceptTPCOnly = false)
  {
    v0 = {}; // safe initialization: start new

    if constexpr (useSelections) {
      // verify track quality
      if (positiveTrack.tpcNClsCrossedRows() < v0selections.minCrossedRows) {
        v0 = {};
        return false;
      }
      if (negativeTrack.tpcNClsCrossedRows() < v0selections.minCrossedRows) {
        v0 = {};
        return false;
      }
      // verify eta
      if (std::fabs(positiveTrack.eta()) > v0selections.maxDaughterEta) {
        v0 = {};
        return false;
      }
      if (std::fabs(negativeTrack.eta()) > v0selections.maxDaughterEta) {
        v0 = {};
        return false;
      }
      if (!acceptTPCOnly && !positiveTrack.hasITS() && !positiveTrack.hasTRD() && !positiveTrack.hasTOF()) {
        v0 = {};
        return false;
      }
      if (!acceptTPCOnly && !negativeTrack.hasITS() && !negativeTrack.hasTRD() && !negativeTrack.hasTOF()) {
        v0 = {};
        return false;
      }
    }

    if constexpr (calculateProngDCAtoPV) {
      // Calculate DCA with respect to the collision associated to the V0
      std::array<float, 2> dcaInfo;

      // do DCA to PV on TrackPar copies and not TrackParCov
      // TrackPar preferred: don't calculate multiple scattering / CovMat changes
      // Spares CPU since variables not checked
      o2::track::TrackPar positiveTrackParamCopy(positiveTrackParam);
      o2::track::TrackPar negativeTrackParamCopy(negativeTrackParam);

      dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted
      o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, positiveTrackParamCopy, 2.f, fitter.getMatCorrType(), &dcaInfo);
      v0.positiveDCAxy = dcaInfo[0];

      if constexpr (useSelections) {
        if (std::fabs(v0.positiveDCAxy) < v0selections.dcapostopv) {
          v0 = {};
          return false;
        }
      }

      dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted
      o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, negativeTrackParamCopy, 2.f, fitter.getMatCorrType(), &dcaInfo);
      v0.negativeDCAxy = dcaInfo[0];

      if constexpr (useSelections) {
        if (std::fabs(v0.negativeDCAxy) < v0selections.dcanegtopv) {
          v0 = {};
          return false;
        }
      }
    } else {
      v0.positiveDCAxy = 0.0f; // default invalid
      v0.negativeDCAxy = 0.0f; // default invalid
    }

    // Perform DCA fit
    int nCand = 0;
    fitter.setCollinear(useCollinearFit);
    try {
      nCand = fitter.process(positiveTrackParam, negativeTrackParam);
    } catch (...) {
      v0 = {};
      fitter.setCollinear(false); // even if returned, reset
      return false;
    }
    if (nCand == 0) {
      v0 = {};
      fitter.setCollinear(false); // even if returned, reset
      return false;
    }
    fitter.setCollinear(false); // proper cleaning: when exiting this loop, always reset to not collinear

    // Calculate DCAToPV of the V0
    o2::track::TrackPar V0Temp = fitter.createParentTrackPar();
    V0Temp.setAbsCharge(0); // charge zero
    std::array<float, 2> dcaV0Info;

    // propagate to collision vertex
    dcaV0Info[0] = dcaV0Info[1] = 999.0f; // default DCA: large, use with care if propagation fails
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, V0Temp, 2.f, fitter.getMatCorrType(), &dcaV0Info);
    v0.v0DCAToPVxy = dcaV0Info[0];
    v0.v0DCAToPVz = dcaV0Info[1];

    v0.positiveTrackX = fitter.getTrack(0).getX();
    v0.negativeTrackX = fitter.getTrack(1).getX();
    positiveTrackParam = fitter.getTrack(0);
    negativeTrackParam = fitter.getTrack(1);
    positiveTrackParam.getPxPyPzGlo(v0.positiveMomentum);
    negativeTrackParam.getPxPyPzGlo(v0.negativeMomentum);
    positiveTrackParam.getXYZGlo(v0.positivePosition);
    negativeTrackParam.getXYZGlo(v0.negativePosition);
    for (int i = 0; i < 3; i++) {
      // avoids misuse if mixed with KF particle use
      v0.momentum[i] = v0.positiveMomentum[i] + v0.negativeMomentum[i];
    }

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      v0.position[i] = vtx[i];
    }

    if constexpr (useSelections) {
      if (std::hypot(v0.position[0], v0.position[1]) < v0selections.v0radius) {
        v0 = {};
        return false;
      }
    }

    v0.daughterDCA = TMath::Sqrt(fitter.getChi2AtPCACandidate());

    if constexpr (useSelections) {
      if (v0.daughterDCA > v0selections.dcav0dau) {
        v0 = {};
        return false;
      }
    }

    double cosPA = RecoDecay::cpa(
      std::array{pvX, pvY, pvZ},
      std::array{v0.position[0], v0.position[1], v0.position[2]},
      std::array{v0.positiveMomentum[0] + v0.negativeMomentum[0], v0.positiveMomentum[1] + v0.negativeMomentum[1], v0.positiveMomentum[2] + v0.negativeMomentum[2]});
    if constexpr (useSelections) {
      if (cosPA < v0selections.v0cospa) {
        v0 = {};
        return false;
      }
    }

    v0.pointingAngle = TMath::ACos(cosPA);
    v0.dcaToPV = CalculateDCAStraightToPV(
      v0.position[0], v0.position[1], v0.position[2],
      v0.positiveMomentum[0] + v0.negativeMomentum[0],
      v0.positiveMomentum[1] + v0.negativeMomentum[1],
      v0.positiveMomentum[2] + v0.negativeMomentum[2],
      pvX, pvY, pvZ);

    // Calculate masses
    v0.massGamma = RecoDecay::m(std::array{
                                  std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                  std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
    v0.massK0Short = RecoDecay::m(std::array{
                                    std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                    std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                  std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
    v0.massLambda = RecoDecay::m(std::array{
                                   std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                   std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                 std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
    v0.massAntiLambda = RecoDecay::m(std::array{
                                       std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                       std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                     std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});

    // calculate covariance if requested
    if (calculateCovariance) {
      // Calculate position covariance matrix
      auto covVtxV = fitter.calcPCACovMatrix(0);
      // std::array<float, 6> positionCovariance;
      v0.positionCovariance[0] = covVtxV(0, 0);
      v0.positionCovariance[1] = covVtxV(1, 0);
      v0.positionCovariance[2] = covVtxV(1, 1);
      v0.positionCovariance[3] = covVtxV(2, 0);
      v0.positionCovariance[4] = covVtxV(2, 1);
      v0.positionCovariance[5] = covVtxV(2, 2);
      std::array<float, 21> covTpositive = {0.};
      std::array<float, 21> covTnegative = {0.};
      positiveTrackParam.getCovXYZPxPyPzGlo(covTpositive);
      negativeTrackParam.getCovXYZPxPyPzGlo(covTnegative);
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        v0.momentumCovariance[i] = covTpositive[MomInd[i]] + covTnegative[MomInd[i]];
      }
    }

    // set collision Id correctly
    v0.collisionId = collisionIndex;

    // information validated, V0 built successfully. Signal OK
    return true;
  }

  //_______________________________________________________________________
  // build V0 with KF function. Populates ::v0 object
  // ::v0 will be initialized to defaults if build fails
  template <typename TCollision, typename TTrack, typename TTrackParametrization>
  bool buildV0CandidateWithKF(TCollision const& collision,
                              TTrack const& positiveTrack,
                              TTrack const& negativeTrack,
                              TTrackParametrization& positiveTrackParam,
                              TTrackParametrization& negativeTrackParam,
                              int kfConstructMethod = 2,           // the typical used
                              float kfConstrainedMassValue = 0.0f, // negative: no constraint
                              bool kfConstrainToPrimaryVertex = true)
  {
    int collisionIndex = collision.globalIndex();
    float pvX = collision.posX();
    float pvY = collision.posY();
    float pvZ = collision.posZ();

    // verify track quality
    if (positiveTrack.tpcNClsCrossedRows() < v0selections.minCrossedRows) {
      v0 = {};
      return false;
    }
    if (negativeTrack.tpcNClsCrossedRows() < v0selections.minCrossedRows) {
      v0 = {};
      return false;
    }

    // verify eta
    if (std::fabs(positiveTrack.eta()) > v0selections.maxDaughterEta) {
      v0 = {};
      return false;
    }
    if (std::fabs(negativeTrack.eta()) > v0selections.maxDaughterEta) {
      v0 = {};
      return false;
    }

    // Calculate DCA with respect to the collision associated to the V0
    std::array<float, 2> dcaInfo;

    // do DCA to PV on TrackPar copies and not TrackParCov
    // TrackPar preferred: don't calculate multiple scattering / CovMat changes
    // Spares CPU since variables not checked
    o2::track::TrackPar positiveTrackParamCopy(positiveTrackParam);
    o2::track::TrackPar negativeTrackParamCopy(negativeTrackParam);

    dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, positiveTrackParamCopy, 2.f, fitter.getMatCorrType(), &dcaInfo);
    v0.positiveDCAxy = dcaInfo[0];

    if (std::fabs(v0.positiveDCAxy) < v0selections.dcanegtopv) {
      v0 = {};
      return false;
    }

    dcaInfo[0] = dcaInfo[1] = 999.0f; // reset to default value
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, negativeTrackParamCopy, 2.f, fitter.getMatCorrType(), &dcaInfo);
    v0.negativeDCAxy = dcaInfo[0];

    if (std::fabs(v0.negativeDCAxy) < v0selections.dcanegtopv) {
      v0 = {};
      return false;
    }

    //__________________________________________
    //*>~<* do V0 with KF
    // create KFParticle objects from trackParCovs
    KFParticle kfpPos = createKFParticleFromTrackParCov(positiveTrackParam, positiveTrackParam.getCharge(), o2::constants::physics::MassElectron);
    KFParticle kfpNeg = createKFParticleFromTrackParCov(negativeTrackParam, negativeTrackParam.getCharge(), o2::constants::physics::MassElectron);

    KFParticle kfpPos_DecayVtx = kfpPos;
    KFParticle kfpNeg_DecayVtx = kfpNeg;
    const KFParticle* V0Daughters[2] = {&kfpPos, &kfpNeg};

    // construct V0
    KFParticle KFV0;
    KFV0.SetConstructMethod(kfConstructMethod);
    try {
      KFV0.Construct(V0Daughters, 2);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
      v0 = {};
      return false;
    }

    if (kfConstrainedMassValue > -1e-4) {
      // photon constraint: this one's got no mass
      KFV0.SetNonlinearMassConstraint(kfConstrainedMassValue);
    }

    // V0 constructed, now recovering TrackParCov for dca fitter minimization (with material correction)
    KFV0.TransportToDecayVertex();
    o2::track::TrackParCov v0TrackParCov = getTrackParCovFromKFP(KFV0, o2::track::PID::Lambda, 0);
    v0TrackParCov.setAbsCharge(0); // to be sure

    // estimate momentum of daughters (since KFParticle does not allow us to retrieve it directly...)
    float xyz_decay[3] = {KFV0.GetX(), KFV0.GetY(), KFV0.GetZ()};
    kfpPos_DecayVtx.TransportToPoint(xyz_decay);
    kfpNeg_DecayVtx.TransportToPoint(xyz_decay);

    v0.positiveMomentum = {kfpPos_DecayVtx.GetPx(), kfpPos_DecayVtx.GetPy(), kfpPos_DecayVtx.GetPz()};
    v0.negativeMomentum = {kfpNeg_DecayVtx.GetPx(), kfpNeg_DecayVtx.GetPy(), kfpNeg_DecayVtx.GetPz()};

    v0.daughterDCA = std::hypot(
      kfpPos_DecayVtx.GetX() - kfpNeg_DecayVtx.GetX(),
      kfpPos_DecayVtx.GetY() - kfpNeg_DecayVtx.GetY(),
      kfpPos_DecayVtx.GetZ() - kfpNeg_DecayVtx.GetZ());

    if (v0.daughterDCA > v0selections.dcav0dau) {
      v0 = {};
      return false;
    }

    // check radius
    for (int i = 0; i < 3; i++) {
      v0.position[i] = xyz_decay[i];
    }
    if (std::hypot(v0.position[0], v0.position[1]) < v0selections.v0radius) {
      v0 = {};
      return false;
    }

    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    // deal with pointing angle
    float cosPA = cpaFromKF(KFV0, KFPV);
    if (cosPA < v0selections.v0cospa) {
      v0 = {};
      return false;
    }
    v0.pointingAngle = TMath::ACos(cosPA);

    v0.dcaToPV = CalculateDCAStraightToPV(
      v0.position[0], v0.position[1], v0.position[2],
      v0.momentum[0], v0.momentum[1], v0.momentum[2],
      pvX, pvY, pvZ);

    // apply topological constraint to PV if requested
    // might adjust px py pz
    KFParticle KFV0_PV = KFV0;
    if (kfConstrainToPrimaryVertex) {
      KFV0_PV.SetProductionVertex(KFPV);
    }
    v0.momentum = {KFV0_PV.GetPx(), KFV0_PV.GetPy(), KFV0_PV.GetPz()};

    // set collision Id correctly
    v0.collisionId = collisionIndex;

    // Calculate masses
    v0.massGamma = RecoDecay::m(std::array{
                                  std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                  std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
    v0.massK0Short = RecoDecay::m(std::array{
                                    std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                    std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                  std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
    v0.massLambda = RecoDecay::m(std::array{
                                   std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                   std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                 std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
    v0.massAntiLambda = RecoDecay::m(std::array{
                                       std::array{v0.positiveMomentum[0], v0.positiveMomentum[1], v0.positiveMomentum[2]},
                                       std::array{v0.negativeMomentum[0], v0.negativeMomentum[1], v0.negativeMomentum[2]}},
                                     std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});

    // information validated, V0 built successfully. Signal OK
    return true;
  }

  //_______________________________________________________________________
  // build Cascade from three tracks, including V0 building.
  // Populates ::cascade object.
  // ::cascade will be initialized to defaults if build fails
  // cascade builder creating a cascade from plain tracks
  template <typename TTrack>
  bool buildCascadeCandidate(int collisionIndex,
                             float pvX, float pvY, float pvZ,
                             TTrack const& positiveTrack,
                             TTrack const& negativeTrack,
                             TTrack const& bachelorTrack,
                             bool calculateBachelorBaryonVariables = false,
                             bool useCascadeMomentumAtPV = false,
                             bool processCovariances = false)
  {
    // no special treatment of positive and negative tracks when building V0s for cascades
    auto posTrackPar = getTrackParCov(positiveTrack);
    auto negTrackPar = getTrackParCov(negativeTrack);

    if (!buildV0Candidate(collisionIndex, pvX, pvY, pvZ, positiveTrack, negativeTrack, posTrackPar, negTrackPar, false, processCovariances, false)) {
      return false;
    }
    if (!buildCascadeCandidate(collisionIndex, pvX, pvY, pvZ, v0, positiveTrack, negativeTrack, bachelorTrack, calculateBachelorBaryonVariables, useCascadeMomentumAtPV, processCovariances)) {
      return false;
    }
    return true;
  }

  //_______________________________________________________________________
  // cascade builder using pre-fabricated information, thus not calling
  // the DCAfitter again for the V0 contained in the cascade
  // If building from scratch, prefer previous version!
  // Populates ::cascade object.
  // ::cascade will be initialized to defaults if build fails
  // cascade builder creating a cascade from plain tracks
  template <typename TTrack>
  bool buildCascadeCandidate(int collisionIndex,
                             float pvX, float pvY, float pvZ,
                             v0candidate const& v0input,
                             TTrack const& positiveTrack,
                             TTrack const& negativeTrack,
                             TTrack const& bachelorTrack,
                             bool calculateBachelorBaryonVariables = false,
                             bool useCascadeMomentumAtPV = false,
                             bool processCovariances = false)
  {
    cascade = {}; // initialize / empty (extra safety)

    // verify track quality
    if (positiveTrack.tpcNClsCrossedRows() < cascadeselections.minCrossedRows) {
      cascade = {};
      return false;
    }
    if (negativeTrack.tpcNClsCrossedRows() < cascadeselections.minCrossedRows) {
      cascade = {};
      return false;
    }
    if (bachelorTrack.tpcNClsCrossedRows() < cascadeselections.minCrossedRows) {
      cascade = {};
      return false;
    }

    // verify eta
    if (std::fabs(positiveTrack.eta()) > cascadeselections.maxDaughterEta) {
      cascade = {};
      return false;
    }
    if (std::fabs(negativeTrack.eta()) > cascadeselections.maxDaughterEta) {
      cascade = {};
      return false;
    }
    if (std::fabs(bachelorTrack.eta()) > cascadeselections.maxDaughterEta) {
      cascade = {};
      return false;
    }

    // verify lambda mass
    if (bachelorTrack.sign() < 0 && std::fabs(v0input.massLambda - 1.116) > cascadeselections.lambdaMassWindow) {
      cascade = {};
      return false;
    }
    if (bachelorTrack.sign() > 0 && std::fabs(v0input.massAntiLambda - 1.116) > cascadeselections.lambdaMassWindow) {
      cascade = {};
      return false;
    }

    if (calculateBachelorBaryonVariables) {
      // Calculates properties of the V0 comprised of bachelor and baryon in the cascade
      // baryon: distinguished via bachelor charge
      if (bachelorTrack.sign() < 0) {
        processBachBaryonVariables(pvX, pvY, pvZ, bachelorTrack, positiveTrack);
      } else {
        processBachBaryonVariables(pvX, pvY, pvZ, bachelorTrack, negativeTrack);
      }
    }

    // Overall cascade charge
    cascade.charge = bachelorTrack.signed1Pt() > 0 ? +1 : -1;

    // bachelor DCA track to PV
    // Calculate DCA with respect to the collision associated to the V0, not individual tracks
    std::array<float, 2> dcaInfo;
    dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted

    auto bachTrackPar = getTrackPar(bachelorTrack);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, bachTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.bachelorDCAxy = dcaInfo[0];

    if (std::fabs(cascade.bachelorDCAxy) < cascadeselections.dcabachtopv) {
      cascade = {};
      return false;
    }

    // Do actual minimization
    auto lBachelorTrack = getTrackParCov(bachelorTrack);

    // Set up covariance matrices (should in fact be optional)
    std::array<float, 21> covV = {0.};
    constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      covV[MomInd[i]] = v0input.momentumCovariance[i];
      covV[i] = v0input.positionCovariance[i];
    }
    auto lV0Track = o2::track::TrackParCov(
      {v0input.position[0], v0input.position[1], v0input.position[2]},
      {v0input.positiveMomentum[0] + v0input.negativeMomentum[0], v0input.positiveMomentum[1] + v0input.negativeMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2]},
      covV, 0, true);
    lV0Track.setAbsCharge(0);
    lV0Track.setPID(o2::track::PID::Lambda);

    //---/---/---/
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(lV0Track, lBachelorTrack);
    } catch (...) {
      cascade = {};
      return false;
    }
    if (nCand == 0) {
      cascade = {};
      return false;
    }

    lV0Track = fitter.getTrack(0);
    lBachelorTrack = fitter.getTrack(1);

    // DCA between cascade daughters
    cascade.cascadeDaughterDCA = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    if (cascade.cascadeDaughterDCA > cascadeselections.dcacascdau) {
      cascade = {};
      return false;
    }

    lBachelorTrack.getPxPyPzGlo(cascade.bachelorMomentum);
    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      cascade.cascadePosition[i] = vtx[i];
    }
    if (std::hypot(cascade.cascadePosition[0], cascade.cascadePosition[1]) < cascadeselections.cascradius) {
      cascade = {};
      return false;
    }

    double cosPA = RecoDecay::cpa(
      std::array{pvX, pvY, pvZ},
      std::array{cascade.cascadePosition[0], cascade.cascadePosition[1], cascade.cascadePosition[2]},
      std::array{v0input.positiveMomentum[0] + v0input.negativeMomentum[0] + cascade.bachelorMomentum[0],
                 v0input.positiveMomentum[1] + v0input.negativeMomentum[1] + cascade.bachelorMomentum[1],
                 v0input.positiveMomentum[2] + v0input.negativeMomentum[2] + cascade.bachelorMomentum[2]});
    if (cosPA < cascadeselections.casccospa) {
      cascade = {};
      return false;
    }
    cascade.pointingAngle = TMath::ACos(cosPA);

    // Calculate DCAxy of the cascade (with bending)
    auto lCascadeTrack = fitter.createParentTrackParCov();
    lCascadeTrack.setAbsCharge(cascade.charge);    // to be sure
    lCascadeTrack.setPID(o2::track::PID::XiMinus); // FIXME: not OK for omegas
    dcaInfo[0] = 999;
    dcaInfo[1] = 999;

    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, lCascadeTrack, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.cascadeDCAxy = dcaInfo[0];
    cascade.cascadeDCAz = dcaInfo[1];

    // Calculate masses a priori
    cascade.massXi = RecoDecay::m(std::array{std::array{cascade.bachelorMomentum[0], cascade.bachelorMomentum[1], cascade.bachelorMomentum[2]},
                                             std::array{v0input.positiveMomentum[0] + v0input.negativeMomentum[0], v0input.positiveMomentum[1] + v0input.negativeMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2]}},
                                  std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
    cascade.massOmega = RecoDecay::m(std::array{std::array{cascade.bachelorMomentum[0], cascade.bachelorMomentum[1], cascade.bachelorMomentum[2]},
                                                std::array{v0input.positiveMomentum[0] + v0input.negativeMomentum[0], v0input.positiveMomentum[1] + v0input.negativeMomentum[1], v0input.positiveMomentum[2] + v0input.negativeMomentum[2]}},
                                     std::array{o2::constants::physics::MassKaonCharged, o2::constants::physics::MassLambda});

    // Populate information
    // cascadecandidate.v0Id = v0index.globalIndex();
    cascade.collisionId = collisionIndex;
    cascade.positiveTrack = positiveTrack.globalIndex();
    cascade.negativeTrack = negativeTrack.globalIndex();
    cascade.bachelorTrack = bachelorTrack.globalIndex();
    cascade.positiveTrackX = v0input.positiveTrackX;
    cascade.negativeTrackX = v0input.negativeTrackX;
    cascade.bachelorTrackX = lBachelorTrack.getX(); // from this minimization
    cascade.v0Position[0] = v0input.position[0];
    cascade.v0Position[1] = v0input.position[1];
    cascade.v0Position[2] = v0input.position[2];
    cascade.v0Momentum[0] = v0input.positiveMomentum[0] + v0input.negativeMomentum[0];
    cascade.v0Momentum[1] = v0input.positiveMomentum[1] + v0input.negativeMomentum[1];
    cascade.v0Momentum[2] = v0input.positiveMomentum[2] + v0input.negativeMomentum[2];
    cascade.positiveMomentum[0] = v0input.positiveMomentum[0];
    cascade.positiveMomentum[1] = v0input.positiveMomentum[1];
    cascade.positiveMomentum[2] = v0input.positiveMomentum[2];
    cascade.negativeMomentum[0] = v0input.negativeMomentum[0];
    cascade.negativeMomentum[1] = v0input.negativeMomentum[1];
    cascade.negativeMomentum[2] = v0input.negativeMomentum[2];
    cascade.v0DaughterDCA = v0input.daughterDCA;
    cascade.positiveDCAxy = v0input.positiveDCAxy;
    cascade.negativeDCAxy = v0input.negativeDCAxy;

    if (useCascadeMomentumAtPV) {
      lCascadeTrack.getPxPyPzGlo(cascade.cascadeMomentum);
    } else {
      cascade.cascadeMomentum[0] = cascade.bachelorMomentum[0] + v0input.positiveMomentum[0] + v0input.negativeMomentum[0];
      cascade.cascadeMomentum[1] = cascade.bachelorMomentum[1] + v0input.positiveMomentum[1] + v0input.negativeMomentum[1];
      cascade.cascadeMomentum[2] = cascade.bachelorMomentum[2] + v0input.positiveMomentum[2] + v0input.negativeMomentum[2];
    }

    if (processCovariances) {
      // Calculate position covariance matrix
      auto covVtxV = fitter.calcPCACovMatrix(0);
      // std::array<float, 6> positionCovariance;
      float positionCovariance[6];
      positionCovariance[0] = covVtxV(0, 0);
      positionCovariance[1] = covVtxV(1, 0);
      positionCovariance[2] = covVtxV(1, 1);
      positionCovariance[3] = covVtxV(2, 0);
      positionCovariance[4] = covVtxV(2, 1);
      positionCovariance[5] = covVtxV(2, 2);
      // store momentum covariance matrix
      std::array<float, 21> covTv0 = {0.};
      std::array<float, 21> covTbachelor = {0.};
      // std::array<float, 6> momentumCovariance;
      lV0Track.getCovXYZPxPyPzGlo(covTv0);
      lBachelorTrack.getCovXYZPxPyPzGlo(covTbachelor);
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 21; i++) {
        cascade.covariance[i] = 0.0f;
      }
      for (int i = 0; i < 6; i++) {
        cascade.covariance[i] = positionCovariance[i];
        cascade.covariance[MomInd[i]] = covTv0[MomInd[i]] + covTbachelor[MomInd[i]];
      }
    }

    // Final outcome is YES if I got here!
    return true;
  }

  //_______________________________________________________________________
  // build KF Cascade from three tracks, including V0 building.
  // Populates ::cascade object.
  // ::cascade will be initialized to defaults if build fails
  // cascade builder creating a cascade from plain tracks
  template <typename TTrack>
  bool buildCascadeCandidateWithKF(int collisionIndex,
                                   float pvX, float pvY, float pvZ,
                                   TTrack const& positiveTrack,
                                   TTrack const& negativeTrack,
                                   TTrack const& bachelorTrack,
                                   bool calculateBachelorBaryonVariables = false,
                                   int kfConstructMethod = 2,
                                   bool kfTuneForOmega = false,
                                   bool kfUseV0MassConstraint = true,
                                   bool kfUseCascadeMassConstraint = false,
                                   bool kfDoDCAFitterPreMinimV0 = false,
                                   bool kfDoDCAFitterPreMinimCasc = false)
  {
    cascade = {}; // initialize / empty (extra safety)

    //*>~<*>~<*>~<*>~<*>~<*>~<*>~<*>~<*>~<*
    // KF particle based rebuilding
    // dispenses prior V0 generation, uses constrained (re-)fit based on bachelor charge
    //*>~<*>~<*>~<*>~<*>~<*>~<*>~<*>~<*>~<*

    if (positiveTrack.tpcNClsCrossedRows() < cascadeselections.minCrossedRows) {
      cascade = {};
      return false;
    }
    if (negativeTrack.tpcNClsCrossedRows() < cascadeselections.minCrossedRows) {
      cascade = {};
      return false;
    }
    if (bachelorTrack.tpcNClsCrossedRows() < cascadeselections.minCrossedRows) {
      cascade = {};
      return false;
    }

    // verify eta
    if (std::fabs(positiveTrack.eta()) > cascadeselections.maxDaughterEta) {
      cascade = {};
      return false;
    }
    if (std::fabs(negativeTrack.eta()) > cascadeselections.maxDaughterEta) {
      cascade = {};
      return false;
    }
    if (std::fabs(bachelorTrack.eta()) > cascadeselections.maxDaughterEta) {
      cascade = {};
      return false;
    }

    if (calculateBachelorBaryonVariables) {
      // Calculates properties of the V0 comprised of bachelor and baryon in the cascade
      // baryon: distinguished via bachelor charge
      if (bachelorTrack.sign() < 0) {
        processBachBaryonVariables(pvX, pvY, pvZ, bachelorTrack, positiveTrack);
      } else {
        processBachBaryonVariables(pvX, pvY, pvZ, bachelorTrack, negativeTrack);
      }
    }

    // Overall cascade charge
    cascade.charge = bachelorTrack.signed1Pt() > 0 ? +1 : -1;

    // bachelor DCA track to PV
    // Calculate DCA with respect to the collision associated to the V0, not individual tracks
    std::array<float, 2> dcaInfo;
    dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted

    auto bachTrackPar = getTrackPar(bachelorTrack);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, bachTrackPar, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.bachelorDCAxy = dcaInfo[0];

    dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted
    o2::track::TrackParCov posTrackParCovForDCA = getTrackParCov(positiveTrack);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, posTrackParCovForDCA, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.positiveDCAxy = dcaInfo[0];

    dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value to make sure candidate accepted
    o2::track::TrackParCov negTrackParCovForDCA = getTrackParCov(negativeTrack);
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, negTrackParCovForDCA, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.negativeDCAxy = dcaInfo[0];

    if (std::fabs(cascade.bachelorDCAxy) < cascadeselections.dcabachtopv) {
      cascade = {};
      return false;
    }

    o2::track::TrackParCov lBachelorTrack = getTrackParCov(bachelorTrack);
    o2::track::TrackParCov posTrackParCov = getTrackParCov(positiveTrack);
    o2::track::TrackParCov negTrackParCov = getTrackParCov(negativeTrack);

    float massPosTrack, massNegTrack;
    if (cascade.charge < 0) {
      massPosTrack = o2::constants::physics::MassProton;
      massNegTrack = o2::constants::physics::MassPionCharged;
    } else {
      massPosTrack = o2::constants::physics::MassPionCharged;
      massNegTrack = o2::constants::physics::MassProton;
    }

    //__________________________________________
    //*>~<* step 1 : V0 with dca fitter, uses material corrections implicitly
    // This is optional - move close to minima and therefore take material
    if (kfDoDCAFitterPreMinimV0) {
      int nCand = 0;
      try {
        nCand = fitter.process(posTrackParCov, negTrackParCov);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        cascade = {};
        return false;
      }
      if (nCand == 0) {
        cascade = {};
        return false;
      }
      // save classical DCA daughters
      cascade.v0DaughterDCA = TMath::Sqrt(fitter.getChi2AtPCACandidate());

      // re-acquire from DCA fitter
      posTrackParCov = fitter.getTrack(0);
      negTrackParCov = fitter.getTrack(1);
    }

    //__________________________________________
    //*>~<* step 2 : V0 with KF
    // create KFParticle objects from trackParCovs
    KFParticle kfpPos = createKFParticleFromTrackParCov(posTrackParCov, posTrackParCov.getCharge(), massPosTrack);
    KFParticle kfpNeg = createKFParticleFromTrackParCov(negTrackParCov, negTrackParCov.getCharge(), massNegTrack);
    const KFParticle* V0Daughters[2] = {&kfpPos, &kfpNeg};

    // construct V0
    KFParticle KFV0;
    KFV0.SetConstructMethod(kfConstructMethod);
    try {
      KFV0.Construct(V0Daughters, 2);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
      cascade = {};
      return false;
    }

    if (kfUseV0MassConstraint) {
      KFV0.SetNonlinearMassConstraint(o2::constants::physics::MassLambda);
    }

    // V0 constructed, now recovering TrackParCov for dca fitter minimization (with material correction)
    KFV0.TransportToDecayVertex();
    o2::track::TrackParCov v0TrackParCov = getTrackParCovFromKFP(KFV0, o2::track::PID::Lambda, 0);
    v0TrackParCov.setAbsCharge(0); // to be sure

    //__________________________________________
    //*>~<* step 3 : Cascade with dca fitter (with material corrections)
    if (kfDoDCAFitterPreMinimCasc) {
      int nCandCascade = 0;
      try {
        nCandCascade = fitter.process(v0TrackParCov, lBachelorTrack);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        cascade = {};
        return false;
      }
      if (nCandCascade == 0) {
        cascade = {};
        return false;
      }

      v0TrackParCov = fitter.getTrack(0);
      lBachelorTrack = fitter.getTrack(1);
    }

    //__________________________________________
    //*>~<* step 4 : Cascade with KF particle (potentially mass-constrained if asked)
    float massBachelorPion = o2::constants::physics::MassPionCharged;
    float massBachelorKaon = o2::constants::physics::MassKaonCharged;

    KFParticle kfpV0 = createKFParticleFromTrackParCov(v0TrackParCov, 0, o2::constants::physics::MassLambda);
    KFParticle kfpBachPion = createKFParticleFromTrackParCov(lBachelorTrack, cascade.charge, massBachelorPion);
    KFParticle kfpBachKaon = createKFParticleFromTrackParCov(lBachelorTrack, cascade.charge, massBachelorKaon);
    const KFParticle* XiDaugthers[2] = {&kfpBachPion, &kfpV0};
    const KFParticle* OmegaDaugthers[2] = {&kfpBachKaon, &kfpV0};

    // construct mother
    KFParticle KFXi, KFOmega;
    KFXi.SetConstructMethod(kfConstructMethod);
    KFOmega.SetConstructMethod(kfConstructMethod);
    try {
      KFXi.Construct(XiDaugthers, 2);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to construct xi from V0 and bachelor track: " << e.what();
      cascade = {};
      return false;
    }
    try {
      KFOmega.Construct(OmegaDaugthers, 2);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to construct omega from V0 and bachelor track: " << e.what();
      cascade = {};
      return false;
    }
    if (kfUseCascadeMassConstraint) {
      // set mass constraint if requested
      // WARNING: this is only adequate for decay chains, i.e. XiC -> Xi or OmegaC -> Omega
      KFXi.SetNonlinearMassConstraint(o2::constants::physics::MassXiMinus);
      KFOmega.SetNonlinearMassConstraint(o2::constants::physics::MassOmegaMinus);
    }
    KFXi.TransportToDecayVertex();
    KFOmega.TransportToDecayVertex();

    // get DCA of daughters at vertex
    cascade.cascadeDaughterDCA = kfpBachPion.GetDistanceFromParticle(kfpV0);
    if (cascade.cascadeDaughterDCA > cascadeselections.dcacascdau) {
      cascade = {};
      return false;
    }

    //__________________________________________
    //*>~<* step 5 : propagate cascade to primary vertex with material corrections if asked

    o2::track::TrackParCov lCascadeTrack;
    if (!kfTuneForOmega) {
      lCascadeTrack = getTrackParCovFromKFP(KFXi, o2::track::PID::XiMinus, cascade.charge);
    } else {
      lCascadeTrack = getTrackParCovFromKFP(KFOmega, o2::track::PID::OmegaMinus, cascade.charge);
    }
    dcaInfo[0] = 999;
    dcaInfo[1] = 999;
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, lCascadeTrack, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.cascadeDCAxy = dcaInfo[0];
    cascade.cascadeDCAz = dcaInfo[1];

    //__________________________________________
    //*>~<* step 6 : acquire all parameters for analysis

    // basic indices
    cascade.collisionId = collisionIndex;
    cascade.v0Id = -1;
    cascade.positiveTrack = positiveTrack.globalIndex();
    cascade.negativeTrack = negativeTrack.globalIndex();
    cascade.bachelorTrack = bachelorTrack.globalIndex();

    // KF chi2
    cascade.kfV0Chi2 = KFV0.GetChi2();
    cascade.kfCascadeChi2 = KFXi.GetChi2();
    if (kfTuneForOmega)
      cascade.kfCascadeChi2 = KFOmega.GetChi2();

    // Daughter momentum not KF-updated FIXME --> but DCA fitter updated if pre-minimisation is used
    lBachelorTrack.getPxPyPzGlo(cascade.bachelorMomentum);
    posTrackParCov.getPxPyPzGlo(cascade.positiveMomentum);
    negTrackParCov.getPxPyPzGlo(cascade.negativeMomentum);

    // Daughter track position at vertex not KF-updated FIXME --> but DCA fitter updated if pre-minimisation is used
    posTrackParCov.getXYZGlo(cascade.positivePosition);
    negTrackParCov.getXYZGlo(cascade.negativePosition);

    // mother position information from KF
    cascade.v0Position[0] = KFV0.GetX();
    cascade.v0Position[1] = KFV0.GetY();
    cascade.v0Position[2] = KFV0.GetZ();

    // mother momentumm information from KF
    cascade.v0Momentum[0] = KFV0.GetPx();
    cascade.v0Momentum[1] = KFV0.GetPy();
    cascade.v0Momentum[2] = KFV0.GetPz();

    // Mother position + momentum is KF updated
    if (!kfTuneForOmega) {
      cascade.cascadePosition[0] = KFXi.GetX();
      cascade.cascadePosition[1] = KFXi.GetY();
      cascade.cascadePosition[2] = KFXi.GetZ();
      cascade.cascadeMomentum[0] = KFXi.GetPx();
      cascade.cascadeMomentum[1] = KFXi.GetPy();
      cascade.cascadeMomentum[2] = KFXi.GetPz();
    } else {
      cascade.cascadePosition[0] = KFOmega.GetX();
      cascade.cascadePosition[1] = KFOmega.GetY();
      cascade.cascadePosition[2] = KFOmega.GetZ();
      cascade.cascadeMomentum[0] = KFOmega.GetPx();
      cascade.cascadeMomentum[1] = KFOmega.GetPy();
      cascade.cascadeMomentum[2] = KFOmega.GetPz();
    }
    if (std::hypot(cascade.cascadePosition[0], cascade.cascadePosition[1]) < cascadeselections.cascradius) {
      cascade = {};
      return false;
    }

    // KF-aware cosPA
    double cosPA = RecoDecay::cpa(
      std::array{pvX, pvY, pvZ},
      std::array{cascade.cascadePosition[0], cascade.cascadePosition[1], cascade.cascadePosition[2]},
      std::array{cascade.cascadeMomentum[0], cascade.cascadeMomentum[1], cascade.cascadeMomentum[2]});
    if (cosPA < cascadeselections.casccospa) {
      cascade = {};
      return false;
    }
    cascade.pointingAngle = TMath::ACos(cosPA);

    // Calculate masses a priori
    float MLambda, SigmaLambda, MXi, SigmaXi, MOmega, SigmaOmega;
    KFV0.GetMass(MLambda, SigmaLambda);
    KFXi.GetMass(MXi, SigmaXi);
    KFOmega.GetMass(MOmega, SigmaOmega);
    cascade.kfMLambda = MLambda;
    cascade.massXi = MXi;
    cascade.massOmega = MOmega;

    // KF Cascade covariance matrix
    std::array<float, 21> covCascKF;
    for (int i = 0; i < 21; i++) { // get covariance matrix elements (lower triangle)
      covCascKF[i] = KFXi.GetCovariance(i);
      cascade.covariance[i] = covCascKF[i];
    }

    // KF V0 covariance matrix
    std::array<float, 21> covV0KF;
    for (int i = 0; i < 21; i++) { // get covariance matrix elements (lower triangle)
      covV0KF[i] = KFV0.GetCovariance(i);
      cascade.kfTrackCovarianceV0[i] = covV0KF[i];
    }

    // V0 daughter covariance matrices
    std::array<float, 21> cvPosKF, cvNegKF;
    posTrackParCov.getCovXYZPxPyPzGlo(cvPosKF);
    negTrackParCov.getCovXYZPxPyPzGlo(cvNegKF);
    for (int i = 0; i < 21; i++) {
      cascade.kfTrackCovariancePos[i] = cvPosKF[i];
      cascade.kfTrackCovarianceNeg[i] = cvNegKF[i];
    }

    return true;
  }

  o2::base::MatLayerCylSet* lut;       // material LUT for DCA fitter
  o2::vertexing::DCAFitterN<2> fitter; // 2-prong o2 dca fitter

  v0candidate v0;           // storage for V0 candidate properties
  cascadeCandidate cascade; // storage for cascade candidate properties

  // v0 candidate criteria
  struct {
    int minCrossedRows;
    float dcanegtopv;
    float dcapostopv;
    double v0cospa;
    float dcav0dau;
    float v0radius;
    float maxDaughterEta;
  } v0selections;

  // cascade candidate criteria
  struct {
    int minCrossedRows;
    float dcabachtopv;
    float cascradius;
    float casccospa;
    float dcacascdau;
    float lambdaMassWindow;
    float maxDaughterEta;
  } cascadeselections;

 private:
  // internal helper to calculate DCA (3D) of a straight line to a given PV analytically
  float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }

  template <typename TTrack>
  void processBachBaryonVariables(float pvX, float pvY, float pvZ, TTrack const& track1, TTrack const& track2)
  {
    cascade.bachBaryonCosPA = 0;        // would ordinarily accept all
    cascade.bachBaryonDCAxyToPV = 1e+3; // would ordinarily accept all

    // create tracks from table rows
    o2::track::TrackParCov tr1 = getTrackParCov(track1);
    o2::track::TrackParCov tr2 = getTrackParCov(track2);

    //---/---/---/
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(tr1, tr2);
    } catch (...) {
      return;
    }
    if (nCand == 0)
      return; // variables are such that candidate is accepted (not obvious...)

    // Calculate DCAxy of the cascade (with bending)
    o2::track::TrackPar wrongV0 = fitter.createParentTrackPar();
    wrongV0.setAbsCharge(0); // charge zero
    std::array<float, 2> dcaInfo;
    dcaInfo[0] = dcaInfo[1] = 999.0f; // by default, take large value

    // bachelor-baryon DCAxy to PV
    o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, wrongV0, 2.f, fitter.getMatCorrType(), &dcaInfo);
    cascade.bachBaryonDCAxyToPV = dcaInfo[0];

    const auto& vtx = fitter.getPCACandidate();
    if (!fitter.isPropagateTracksToVertexDone())
      return;

    std::array<float, 3> tr1p;
    std::array<float, 3> tr2p;

    fitter.getTrack(1).getPxPyPzGlo(tr1p);
    fitter.getTrack(2).getPxPyPzGlo(tr2p);

    // bachelor-baryon CosPA
    cascade.bachBaryonCosPA = RecoDecay::cpa(
      std::array{pvX, pvY, pvZ},
      std::array{vtx[0], vtx[1], vtx[2]},
      std::array{tr1p[0] + tr2p[0], tr1p[1] + tr2p[1], tr1p[2] + tr2p[2]});

    // Potentially also to be considered: bachelor-baryon DCA (between the two tracks)
    // to be added here as complementary information in the future
  }
};

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_STRANGENESSBUILDERHELPER_H_
