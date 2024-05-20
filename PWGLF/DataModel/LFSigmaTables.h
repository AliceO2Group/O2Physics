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

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

#ifndef PWGLF_DATAMODEL_LFSIGMATABLES_H_
#define PWGLF_DATAMODEL_LFSIGMATABLES_H_

using namespace o2;
using namespace o2::framework;

// Creating output TTree for sigma analysis
namespace o2::aod
{

// for real data
namespace v0SigmaCandidate
{
DECLARE_SOA_COLUMN(SigmapT, sigmapT, float);
DECLARE_SOA_COLUMN(SigmaMass, sigmaMass, float);
// DECLARE_SOA_COLUMN(SigmaDCAz, sigmaDCAz, float);
// DECLARE_SOA_COLUMN(SigmaDCAxy, sigmaDCAxy, float);
// DECLARE_SOA_COLUMN(SigmaDCADau, sigmaDCADau, float);
DECLARE_SOA_COLUMN(PhotonPt, photonPt, float);
DECLARE_SOA_COLUMN(PhotonMass, photonMass, float);
DECLARE_SOA_COLUMN(PhotonQt, photonQt, float);
DECLARE_SOA_COLUMN(PhotonAlpha, photonAlpha, float);
DECLARE_SOA_COLUMN(PhotonRadius, photonRadius, float);
DECLARE_SOA_COLUMN(PhotonCosPA, photonCosPA, float);
DECLARE_SOA_COLUMN(PhotonDCADau, photonDCADau, float);
DECLARE_SOA_COLUMN(PhotonDCANegPV, photonDCANegPV, float);
DECLARE_SOA_COLUMN(PhotonDCAPosPV, photonDCAPosPV, float);
DECLARE_SOA_COLUMN(PhotonZconv, photonZconv, float);
DECLARE_SOA_COLUMN(LambdaPt, lambdaPt, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(LambdaQt, lambdaQt, float);
DECLARE_SOA_COLUMN(LambdaAlpha, lambdaAlpha, float);
DECLARE_SOA_COLUMN(LambdaRadius, lambdaRadius, float);
DECLARE_SOA_COLUMN(LambdaCosPA, lambdaCosPA, float);
DECLARE_SOA_COLUMN(LambdaDCADau, lambdaDCADau, float);
DECLARE_SOA_COLUMN(LambdaDCANegPV, lambdaDCANegPV, float);
DECLARE_SOA_COLUMN(LambdaDCAPosPV, lambdaDCAPosPV, float);
DECLARE_SOA_COLUMN(GammaProbability, gammaProbability, float);
DECLARE_SOA_COLUMN(LambdaProbability, lambdaProbability, float);

} // namespace v0SigmaCandidate

DECLARE_SOA_TABLE(V0SigmaCandidates, "AOD", "V0SIGMAS",
                  v0SigmaCandidate::SigmapT,
                  v0SigmaCandidate::SigmaMass,
                  // v0SigmaCandidate::SigmaDCAz,
                  // v0SigmaCandidate::SigmaDCAxy,
                  // v0SigmaCandidate::SigmaDCADau,
                  v0SigmaCandidate::PhotonPt,
                  v0SigmaCandidate::PhotonMass,
                  v0SigmaCandidate::PhotonQt,
                  v0SigmaCandidate::PhotonAlpha,
                  v0SigmaCandidate::PhotonRadius,
                  v0SigmaCandidate::PhotonCosPA,
                  v0SigmaCandidate::PhotonDCADau,
                  v0SigmaCandidate::PhotonDCANegPV,
                  v0SigmaCandidate::PhotonDCAPosPV,
                  v0SigmaCandidate::PhotonZconv,
                  v0SigmaCandidate::LambdaPt,
                  v0SigmaCandidate::LambdaMass,
                  v0SigmaCandidate::LambdaQt,
                  v0SigmaCandidate::LambdaAlpha,
                  v0SigmaCandidate::LambdaRadius,
                  v0SigmaCandidate::LambdaCosPA,
                  v0SigmaCandidate::LambdaDCADau,
                  v0SigmaCandidate::LambdaDCANegPV,
                  v0SigmaCandidate::LambdaDCAPosPV,
                  v0SigmaCandidate::GammaProbability,
                  v0SigmaCandidate::LambdaProbability);

// for MC data
namespace v0SigmaMCCandidate
{
DECLARE_SOA_COLUMN(IsSigma, isSigma, bool);

} // namespace v0SigmaMCCandidate

DECLARE_SOA_TABLE(V0SigmaMCCandidates, "AOD", "V0MCSIGMAS",
                  v0SigmaMCCandidate::IsSigma);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMATABLES_H_
