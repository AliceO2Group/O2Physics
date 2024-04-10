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


// Creating output TTree for ML analysis
namespace o2::aod 
{

// for real data
namespace V0MLSigmaCandidate
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

} // namespace V0MLSigmaCandidate

DECLARE_SOA_TABLE(V0MLSigmaCandidates, "AOD", "V0MLSIGMAS", 
				V0MLSigmaCandidate::SigmapT,
				V0MLSigmaCandidate::SigmaMass,
				// V0MLSigmaCandidate::SigmaDCAz,
				// V0MLSigmaCandidate::SigmaDCAxy,
				// V0MLSigmaCandidate::SigmaDCADau,
				V0MLSigmaCandidate::PhotonPt,
				V0MLSigmaCandidate::PhotonMass,
				V0MLSigmaCandidate::PhotonQt,
				V0MLSigmaCandidate::PhotonAlpha,
				V0MLSigmaCandidate::PhotonRadius,
				V0MLSigmaCandidate::PhotonCosPA,
				V0MLSigmaCandidate::PhotonDCADau,
				V0MLSigmaCandidate::PhotonDCANegPV,
				V0MLSigmaCandidate::PhotonDCAPosPV,
				V0MLSigmaCandidate::PhotonZconv,
				V0MLSigmaCandidate::LambdaPt,
				V0MLSigmaCandidate::LambdaMass,
				V0MLSigmaCandidate::LambdaQt,
				V0MLSigmaCandidate::LambdaAlpha,
				V0MLSigmaCandidate::LambdaRadius,
				V0MLSigmaCandidate::LambdaCosPA,
				V0MLSigmaCandidate::LambdaDCADau,
				V0MLSigmaCandidate::LambdaDCANegPV,
				V0MLSigmaCandidate::LambdaDCAPosPV,
				V0MLSigmaCandidate::GammaProbability,
				V0MLSigmaCandidate::LambdaProbability);

// for MC data
namespace V0MLSigmaMCCandidate 
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
DECLARE_SOA_COLUMN(IsSigma, isSigma, bool);

} // namespace V0MLSigmaMCCandidate

DECLARE_SOA_TABLE(V0MLSigmaMCCandidates, "AOD", "V0MLMCSIGMAS",
				V0MLSigmaMCCandidate::SigmapT,
				V0MLSigmaMCCandidate::SigmaMass,
				// V0MLSigmaMCCandidate::SigmaDCAz,
				// V0MLSigmaMCCandidate::SigmaDCAxy,
				// V0MLSigmaMCCandidate::SigmaDCADau,
				V0MLSigmaMCCandidate::PhotonPt,
				V0MLSigmaMCCandidate::PhotonMass,
				V0MLSigmaMCCandidate::PhotonQt,
				V0MLSigmaMCCandidate::PhotonAlpha,
				V0MLSigmaMCCandidate::PhotonRadius,
				V0MLSigmaMCCandidate::PhotonCosPA,
				V0MLSigmaMCCandidate::PhotonDCADau,
				V0MLSigmaMCCandidate::PhotonDCANegPV,
				V0MLSigmaMCCandidate::PhotonDCAPosPV,
				V0MLSigmaMCCandidate::PhotonZconv,
				V0MLSigmaMCCandidate::LambdaPt,
				V0MLSigmaMCCandidate::LambdaMass,
				V0MLSigmaMCCandidate::LambdaQt,
				V0MLSigmaMCCandidate::LambdaAlpha,
				V0MLSigmaMCCandidate::LambdaRadius,
				V0MLSigmaMCCandidate::LambdaCosPA,
				V0MLSigmaMCCandidate::LambdaDCADau,
				V0MLSigmaMCCandidate::LambdaDCANegPV,
				V0MLSigmaMCCandidate::LambdaDCAPosPV,
				V0MLSigmaMCCandidate::GammaProbability,
				V0MLSigmaMCCandidate::LambdaProbability,
				V0MLSigmaMCCandidate::IsSigma);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMATABLES_H_