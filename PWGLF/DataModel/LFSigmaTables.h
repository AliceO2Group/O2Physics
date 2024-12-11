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

namespace sigma0Core
{
//______________________________________________________
// REGULAR COLUMNS FOR INDEXING
// FOR DERIVED
DECLARE_SOA_INDEX_COLUMN(Sigma0Collision, sigma0Collision); //!
} // namespace sigma0Core

// for real data
namespace sigma0Core
{
DECLARE_SOA_COLUMN(SigmapT, sigmapT, float);
DECLARE_SOA_COLUMN(SigmaMass, sigmaMass, float);
DECLARE_SOA_COLUMN(SigmaRapidity, sigmaRapidity, float);
DECLARE_SOA_COLUMN(SigmaOPAngle, sigmaOPAngle, float);
DECLARE_SOA_COLUMN(SigmaDeltaEta, sigmaDeltaEta, float);
DECLARE_SOA_COLUMN(SigmaDeltaPhi, sigmaDeltaPhi, float);
} // namespace sigma0Core

DECLARE_SOA_TABLE(Sigma0Cores, "AOD", "SIGMA0CORES",
                  sigma0Core::SigmapT,
                  sigma0Core::SigmaMass,
                  sigma0Core::SigmaRapidity,
                  sigma0Core::SigmaOPAngle,
                  sigma0Core::SigmaDeltaEta,
                  sigma0Core::SigmaDeltaPhi);

DECLARE_SOA_TABLE(Sigma0CollRefs, "AOD", "SIGMA0COLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, sigma0Core::Sigma0CollisionId);

// For Photon extra info
namespace sigmaPhotonExtra
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
DECLARE_SOA_COLUMN(PhotonPhi, photonPhi, float);
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

} // namespace sigmaPhotonExtra

DECLARE_SOA_TABLE(SigmaPhotonExtras, "AOD", "SIGMA0PHOTON",
                  sigmaPhotonExtra::PhotonPt,
                  sigmaPhotonExtra::PhotonMass,
                  sigmaPhotonExtra::PhotonQt,
                  sigmaPhotonExtra::PhotonAlpha,
                  sigmaPhotonExtra::PhotonRadius,
                  sigmaPhotonExtra::PhotonCosPA,
                  sigmaPhotonExtra::PhotonDCADau,
                  sigmaPhotonExtra::PhotonDCANegPV,
                  sigmaPhotonExtra::PhotonDCAPosPV,
                  sigmaPhotonExtra::PhotonZconv,
                  sigmaPhotonExtra::PhotonEta,
                  sigmaPhotonExtra::PhotonY,
                  sigmaPhotonExtra::PhotonPhi,
                  sigmaPhotonExtra::PhotonPosTPCNSigma,
                  sigmaPhotonExtra::PhotonNegTPCNSigma,
                  sigmaPhotonExtra::PhotonPosTPCCrossedRows,
                  sigmaPhotonExtra::PhotonNegTPCCrossedRows,
                  sigmaPhotonExtra::PhotonPosPt,
                  sigmaPhotonExtra::PhotonNegPt,
                  sigmaPhotonExtra::PhotonPosEta,
                  sigmaPhotonExtra::PhotonNegEta,
                  sigmaPhotonExtra::PhotonPosY,
                  sigmaPhotonExtra::PhotonNegY,
                  sigmaPhotonExtra::PhotonPsiPair,
                  sigmaPhotonExtra::PhotonPosITSCls,
                  sigmaPhotonExtra::PhotonNegITSCls,
                  sigmaPhotonExtra::PhotonPosITSClSize,
                  sigmaPhotonExtra::PhotonNegITSClSize,
                  sigmaPhotonExtra::PhotonV0Type,
                  sigmaPhotonExtra::GammaBDTScore);

// For Lambda extra info
namespace sigmaLambdaExtra
{
DECLARE_SOA_COLUMN(LambdaPt, lambdaPt, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(AntiLambdaMass, antilambdaMass, float);
DECLARE_SOA_COLUMN(LambdaQt, lambdaQt, float);
DECLARE_SOA_COLUMN(LambdaAlpha, lambdaAlpha, float);
DECLARE_SOA_COLUMN(LambdaRadius, lambdaRadius, float);
DECLARE_SOA_COLUMN(LambdaCosPA, lambdaCosPA, float);
DECLARE_SOA_COLUMN(LambdaDCADau, lambdaDCADau, float);
DECLARE_SOA_COLUMN(LambdaDCANegPV, lambdaDCANegPV, float);
DECLARE_SOA_COLUMN(LambdaDCAPosPV, lambdaDCAPosPV, float);
DECLARE_SOA_COLUMN(LambdaEta, lambdaEta, float);
DECLARE_SOA_COLUMN(LambdaY, lambdaY, float);
DECLARE_SOA_COLUMN(LambdaPhi, lambdaPhi, float);
DECLARE_SOA_COLUMN(LambdaPosPrTPCNSigma, lambdaPosPrTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaPosPiTPCNSigma, lambdaPosPiTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaNegPrTPCNSigma, lambdaNegPrTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaNegPiTPCNSigma, lambdaNegPiTPCNSigma, float);
DECLARE_SOA_COLUMN(LambdaPrTOFNSigma, lambdaPrTOFNSigma, float);
DECLARE_SOA_COLUMN(LambdaPiTOFNSigma, lambdaPiTOFNSigma, float);
DECLARE_SOA_COLUMN(ALambdaPrTOFNSigma, aLambdaPrTOFNSigma, float);
DECLARE_SOA_COLUMN(ALambdaPiTOFNSigma, aLambdaPiTOFNSigma, float);
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

} // namespace sigmaLambdaExtra

DECLARE_SOA_TABLE(SigmaLambdaExtras, "AOD", "SIGMA0LAMBDA",
                  sigmaLambdaExtra::LambdaPt,
                  sigmaLambdaExtra::LambdaMass,
                  sigmaLambdaExtra::AntiLambdaMass,
                  sigmaLambdaExtra::LambdaQt,
                  sigmaLambdaExtra::LambdaAlpha,
                  sigmaLambdaExtra::LambdaRadius,
                  sigmaLambdaExtra::LambdaCosPA,
                  sigmaLambdaExtra::LambdaDCADau,
                  sigmaLambdaExtra::LambdaDCANegPV,
                  sigmaLambdaExtra::LambdaDCAPosPV,
                  sigmaLambdaExtra::LambdaEta,
                  sigmaLambdaExtra::LambdaY,
                  sigmaLambdaExtra::LambdaPhi,
                  sigmaLambdaExtra::LambdaPosPrTPCNSigma,
                  sigmaLambdaExtra::LambdaPosPiTPCNSigma,
                  sigmaLambdaExtra::LambdaNegPrTPCNSigma,
                  sigmaLambdaExtra::LambdaNegPiTPCNSigma,
                  sigmaLambdaExtra::LambdaPrTOFNSigma,
                  sigmaLambdaExtra::LambdaPiTOFNSigma,
                  sigmaLambdaExtra::ALambdaPrTOFNSigma,
                  sigmaLambdaExtra::ALambdaPiTOFNSigma,
                  sigmaLambdaExtra::LambdaPosTPCCrossedRows,
                  sigmaLambdaExtra::LambdaNegTPCCrossedRows,
                  sigmaLambdaExtra::LambdaPosPt,
                  sigmaLambdaExtra::LambdaNegPt,
                  sigmaLambdaExtra::LambdaPosEta,
                  sigmaLambdaExtra::LambdaNegEta,
                  sigmaLambdaExtra::LambdaPosPrY,
                  sigmaLambdaExtra::LambdaPosPiY,
                  sigmaLambdaExtra::LambdaNegPrY,
                  sigmaLambdaExtra::LambdaNegPiY,
                  sigmaLambdaExtra::LambdaPosITSCls,
                  sigmaLambdaExtra::LambdaNegITSCls,
                  sigmaLambdaExtra::LambdaPosITSClSize,
                  sigmaLambdaExtra::LambdaNegITSClSize,
                  sigmaLambdaExtra::LambdaV0Type,
                  sigmaLambdaExtra::LambdaBDTScore,
                  sigmaLambdaExtra::AntiLambdaBDTScore);

// for MC data
namespace sigmaMCCore
{
DECLARE_SOA_COLUMN(IsSigma, isSigma, bool); // TODO: include PDG + IsPhysicalPrimary
DECLARE_SOA_COLUMN(IsAntiSigma, isAntiSigma, bool);
DECLARE_SOA_COLUMN(PhotonCandPDGCode, photonCandPDGCode, int);
DECLARE_SOA_COLUMN(PhotonCandPDGCodeMother, photonCandPDGCodeMother, int);
DECLARE_SOA_COLUMN(IsPhotonCandPrimary, isPhotonCandPrimary, bool);
DECLARE_SOA_COLUMN(LambdaCandPDGCode, lambdaCandPDGCode, int);
DECLARE_SOA_COLUMN(LambdaCandPDGCodeMother, lambdaCandPDGCodeMother, int);
DECLARE_SOA_COLUMN(IsLambdaCandPrimary, isLambdaCandPrimary, bool);

} // namespace sigmaMCCore

DECLARE_SOA_TABLE(SigmaMCCores, "AOD", "SIGMA0MCCORES",
                  sigmaMCCore::IsSigma,
                  sigmaMCCore::IsAntiSigma,
                  sigmaMCCore::PhotonCandPDGCode,
                  sigmaMCCore::PhotonCandPDGCodeMother,
                  sigmaMCCore::IsPhotonCandPrimary,
                  sigmaMCCore::LambdaCandPDGCode,
                  sigmaMCCore::LambdaCandPDGCodeMother,
                  sigmaMCCore::IsLambdaCandPrimary);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMATABLES_H_
