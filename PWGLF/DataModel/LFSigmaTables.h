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
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"
#include "CommonConstants/PhysicsConstants.h"

#ifndef PWGLF_DATAMODEL_LFSIGMATABLES_H_
#define PWGLF_DATAMODEL_LFSIGMATABLES_H_

// Creating output TTree for sigma analysis
namespace o2::aod
{
DECLARE_SOA_TABLE(Sigma0Collisions, "AOD", "SIGMA0COLLISION", //! basic collision properties: position
                  o2::soa::Index<>, collision::PosX, collision::PosY, collision::PosZ,
                  cent::CentFT0M, cent::CentFT0A, cent::CentFT0C, cent::CentFV0A);

using Sigma0Collision = Sigma0Collisions::iterator;

namespace v0SigmaCandidate
{
//______________________________________________________
// REGULAR COLUMNS FOR INDEXING
// FOR DERIVED
DECLARE_SOA_INDEX_COLUMN(Sigma0Collision, sigma0Collision); //!
} // namespace v0SigmaCandidate

// for real data
namespace v0SigmaCandidate
{
DECLARE_SOA_COLUMN(SigmapT, sigmapT, float);
DECLARE_SOA_COLUMN(SigmaMass, sigmaMass, float);
DECLARE_SOA_COLUMN(SigmaRapidity, sigmaRapidity, float);
// DECLARE_SOA_COLUMN(SigmaDCAz, sigmaDCAz, float);
// DECLARE_SOA_COLUMN(SigmaDCAxy, sigmaDCAxy, float);
// DECLARE_SOA_COLUMN(SigmaDCADau, sigmaDCADau, float);

} // namespace v0SigmaCandidate

DECLARE_SOA_TABLE(V0SigmaCandidates, "AOD", "V0SIGMAS",
                  v0SigmaCandidate::SigmapT,
                  v0SigmaCandidate::SigmaMass,
                  v0SigmaCandidate::SigmaRapidity);

DECLARE_SOA_TABLE(V0Sigma0CollRefs, "AOD", "V0SIGMA0COLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0SigmaCandidate::Sigma0CollisionId);

// For Photon extra info
namespace v0SigmaPhotonExtras
{
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
DECLARE_SOA_COLUMN(PhotonEta, photonEta, float);
DECLARE_SOA_COLUMN(PhotonY, photonY, float);
DECLARE_SOA_COLUMN(PhotonPosTPCNSigma, photonPosTPCNSigma, float);
DECLARE_SOA_COLUMN(PhotonNegTPCNSigma, photonNegTPCNSigma, float);
DECLARE_SOA_COLUMN(PhotonPosTPCCrossedRows, photonPosTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(PhotonNegTPCCrossedRows, photonNegTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(PhotonPosPt, photonPosPt, float);
DECLARE_SOA_COLUMN(PhotonNegPt, photonNegPt, float);
DECLARE_SOA_COLUMN(PhotonPosEta, photonPosEta, float);
DECLARE_SOA_COLUMN(PhotonNegEta, photonNegEta, float);
DECLARE_SOA_COLUMN(PhotonPosY, photonPosY, float);
DECLARE_SOA_COLUMN(PhotonNegY, photonNegY, float);
DECLARE_SOA_COLUMN(PhotonPsiPair, photonPsiPair, float);
DECLARE_SOA_COLUMN(PhotonPosITSCls, photonPosITSCls, int);
DECLARE_SOA_COLUMN(PhotonNegITSCls, photonNegITSCls, int);
DECLARE_SOA_COLUMN(PhotonPosITSClSize, photonPosITSClSize, uint32_t);
DECLARE_SOA_COLUMN(PhotonNegITSClSize, photonNegITSClSize, uint32_t);
DECLARE_SOA_COLUMN(PhotonV0Type, photonV0Type, uint8_t);
DECLARE_SOA_COLUMN(GammaBDTScore, gammaBDTScore, float);

} // namespace v0SigmaPhotonExtras

DECLARE_SOA_TABLE(V0SigmaPhotonExtras, "AOD", "V0SIGMAPHOTON",
                  v0SigmaPhotonExtras::PhotonMass,
                  v0SigmaPhotonExtras::PhotonPt,
                  v0SigmaPhotonExtras::PhotonQt,
                  v0SigmaPhotonExtras::PhotonAlpha,
                  v0SigmaPhotonExtras::PhotonRadius,
                  v0SigmaPhotonExtras::PhotonCosPA,
                  v0SigmaPhotonExtras::PhotonDCADau,
                  v0SigmaPhotonExtras::PhotonDCANegPV,
                  v0SigmaPhotonExtras::PhotonDCAPosPV,
                  v0SigmaPhotonExtras::PhotonZconv,
                  v0SigmaPhotonExtras::PhotonEta,
                  v0SigmaPhotonExtras::PhotonY,
                  v0SigmaPhotonExtras::PhotonPosTPCNSigma,
                  v0SigmaPhotonExtras::PhotonNegTPCNSigma,
                  v0SigmaPhotonExtras::PhotonPosTPCCrossedRows,
                  v0SigmaPhotonExtras::PhotonNegTPCCrossedRows,
                  v0SigmaPhotonExtras::PhotonPosPt,
                  v0SigmaPhotonExtras::PhotonNegPt,
                  v0SigmaPhotonExtras::PhotonPosEta,
                  v0SigmaPhotonExtras::PhotonNegEta,
                  v0SigmaPhotonExtras::PhotonPosY,
                  v0SigmaPhotonExtras::PhotonNegY,
                  v0SigmaPhotonExtras::PhotonPsiPair,
                  v0SigmaPhotonExtras::PhotonPosITSCls,
                  v0SigmaPhotonExtras::PhotonNegITSCls,
                  v0SigmaPhotonExtras::PhotonPosITSClSize,
                  v0SigmaPhotonExtras::PhotonNegITSClSize,
                  v0SigmaPhotonExtras::PhotonV0Type,
                  v0SigmaPhotonExtras::GammaBDTScore);

// For Lambda extra info
namespace v0SigmaLambdaExtras
{
DECLARE_SOA_COLUMN(LambdaPt, lambdaPt, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(LambdaQt, lambdaQt, float);
DECLARE_SOA_COLUMN(LambdaAlpha, lambdaAlpha, float);
DECLARE_SOA_COLUMN(LambdaRadius, lambdaRadius, float);
DECLARE_SOA_COLUMN(LambdaCosPA, lambdaCosPA, float);
DECLARE_SOA_COLUMN(LambdaDCADau, lambdaDCADau, float);
DECLARE_SOA_COLUMN(LambdaDCANegPV, lambdaDCANegPV, float);
DECLARE_SOA_COLUMN(LambdaDCAPosPV, lambdaDCAPosPV, float);
DECLARE_SOA_COLUMN(LambdaEta, lambdaEta, float);
DECLARE_SOA_COLUMN(LambdaY, lambdaY, float);
DECLARE_SOA_COLUMN(LambdaPosPrTPCNSigma, lambdaPosPrTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaPosPiTPCNSigma, lambdaPosPiTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaNegPrTPCNSigma, lambdaNegPrTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaNegPiTPCNSigma, lambdaNegPiTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaPosTPCCrossedRows, lambdaPosTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(LambdaNegTPCCrossedRows, lambdaNegTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(LambdaPosPt, lambdaPosPt, float);
DECLARE_SOA_COLUMN(LambdaNegPt, lambdaNegPt, float);
DECLARE_SOA_COLUMN(LambdaPosEta, lambdaPosEta, float);
DECLARE_SOA_COLUMN(LambdaNegEta, lambdaNegEta, float);
DECLARE_SOA_COLUMN(LambdaPosPrY, lambdaPosPrY, float);
DECLARE_SOA_COLUMN(LambdaPosPiY, lambdaPosPiY, float);
DECLARE_SOA_COLUMN(LambdaNegPrY, lambdaNegPrY, float);
DECLARE_SOA_COLUMN(LambdaNegPiY, lambdaNegPiY, float);
DECLARE_SOA_COLUMN(LambdaPosITSCls, lambdaPosITSCls, int);
DECLARE_SOA_COLUMN(LambdaNegITSCls, lambdaNegITSCls, int);
DECLARE_SOA_COLUMN(LambdaPosITSClSize, lambdaPosITSClSize, uint32_t);
DECLARE_SOA_COLUMN(LambdaNegITSClSize, lambdaNegITSClSize, uint32_t);
DECLARE_SOA_COLUMN(LambdaV0Type, lambdaV0Type, uint8_t);
DECLARE_SOA_COLUMN(LambdaBDTScore, lambdaBDTScore, float);
DECLARE_SOA_COLUMN(AntiLambdaBDTScore, antilambdaBDTScore, float);

} // namespace v0SigmaLambdaExtras

DECLARE_SOA_TABLE(V0SigmaLambdaExtras, "AOD", "V0SIGMALAMBDA",
                  v0SigmaLambdaExtras::LambdaPt,
                  v0SigmaLambdaExtras::LambdaMass,
                  v0SigmaLambdaExtras::LambdaQt,
                  v0SigmaLambdaExtras::LambdaAlpha,
                  v0SigmaLambdaExtras::LambdaRadius,
                  v0SigmaLambdaExtras::LambdaCosPA,
                  v0SigmaLambdaExtras::LambdaDCADau,
                  v0SigmaLambdaExtras::LambdaDCANegPV,
                  v0SigmaLambdaExtras::LambdaDCAPosPV,
                  v0SigmaLambdaExtras::LambdaEta,
                  v0SigmaLambdaExtras::LambdaY,
                  v0SigmaLambdaExtras::LambdaPosPrTPCNSigma,
                  v0SigmaLambdaExtras::LambdaPosPiTPCNSigma,
                  v0SigmaLambdaExtras::LambdaNegPrTPCNSigma,
                  v0SigmaLambdaExtras::LambdaNegPiTPCNSigma,
                  v0SigmaLambdaExtras::LambdaPosTPCCrossedRows,
                  v0SigmaLambdaExtras::LambdaNegTPCCrossedRows,
                  v0SigmaLambdaExtras::LambdaPosPt,
                  v0SigmaLambdaExtras::LambdaNegPt,
                  v0SigmaLambdaExtras::LambdaPosEta,
                  v0SigmaLambdaExtras::LambdaNegEta,
                  v0SigmaLambdaExtras::LambdaPosPrY,
                  v0SigmaLambdaExtras::LambdaPosPiY,
                  v0SigmaLambdaExtras::LambdaNegPrY,
                  v0SigmaLambdaExtras::LambdaNegPiY,
                  v0SigmaLambdaExtras::LambdaPosITSCls,
                  v0SigmaLambdaExtras::LambdaNegITSCls,
                  v0SigmaLambdaExtras::LambdaPosITSClSize,
                  v0SigmaLambdaExtras::LambdaNegITSClSize,
                  v0SigmaLambdaExtras::LambdaV0Type,
                  v0SigmaLambdaExtras::LambdaBDTScore,
                  v0SigmaLambdaExtras::AntiLambdaBDTScore);

// for MC data
namespace v0SigmaMCCandidate
{
DECLARE_SOA_COLUMN(IsSigma, isSigma, bool);

} // namespace v0SigmaMCCandidate

DECLARE_SOA_TABLE(V0SigmaMCCandidates, "AOD", "V0MCSIGMAS",
                  v0SigmaMCCandidate::IsSigma);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMATABLES_H_
