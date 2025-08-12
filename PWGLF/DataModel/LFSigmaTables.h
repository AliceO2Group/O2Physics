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

// for real data
namespace sigma0Core
{
DECLARE_SOA_COLUMN(SigmapT, sigmapT, float);
DECLARE_SOA_COLUMN(SigmaMass, sigmaMass, float);
DECLARE_SOA_COLUMN(SigmaRapidity, sigmaRapidity, float);
DECLARE_SOA_COLUMN(SigmaOPAngle, sigmaOPAngle, float);

} // namespace sigma0Core

DECLARE_SOA_TABLE(Sigma0Cores, "AOD", "SIGMA0CORES",
                  sigma0Core::SigmapT,
                  sigma0Core::SigmaMass,
                  sigma0Core::SigmaRapidity,
                  sigma0Core::SigmaOPAngle);

// For Photon extra info
namespace sigma0PhotonExtra
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
DECLARE_SOA_COLUMN(PhotonPosTPCNSigmaEl, photonPosTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(PhotonNegTPCNSigmaEl, photonNegTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(PhotonPosTPCNSigmaPi, photonPosTPCNSigmaPi, float);
DECLARE_SOA_COLUMN(PhotonNegTPCNSigmaPi, photonNegTPCNSigmaPi, float);
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
DECLARE_SOA_COLUMN(PhotonPosITSChi2PerNcl, photonPosITSChi2PerNcl, float);
DECLARE_SOA_COLUMN(PhotonNegITSChi2PerNcl, photonNegITSChi2PerNcl, float);
DECLARE_SOA_COLUMN(PhotonPosTrackCode, photonPosTrackCode, uint8_t);
DECLARE_SOA_COLUMN(PhotonNegTrackCode, photonNegTrackCode, uint8_t);
DECLARE_SOA_COLUMN(PhotonV0Type, photonV0Type, uint8_t);
DECLARE_SOA_COLUMN(GammaBDTScore, gammaBDTScore, float);

} // namespace sigma0PhotonExtra

DECLARE_SOA_TABLE(Sigma0PhotonExtras, "AOD", "SIGMA0PHOTON",
                  sigma0PhotonExtra::PhotonPt,
                  sigma0PhotonExtra::PhotonMass,
                  sigma0PhotonExtra::PhotonQt,
                  sigma0PhotonExtra::PhotonAlpha,
                  sigma0PhotonExtra::PhotonRadius,
                  sigma0PhotonExtra::PhotonCosPA,
                  sigma0PhotonExtra::PhotonDCADau,
                  sigma0PhotonExtra::PhotonDCANegPV,
                  sigma0PhotonExtra::PhotonDCAPosPV,
                  sigma0PhotonExtra::PhotonZconv,
                  sigma0PhotonExtra::PhotonEta,
                  sigma0PhotonExtra::PhotonY,
                  sigma0PhotonExtra::PhotonPhi,
                  sigma0PhotonExtra::PhotonPosTPCNSigmaEl,
                  sigma0PhotonExtra::PhotonNegTPCNSigmaEl,
                  sigma0PhotonExtra::PhotonPosTPCNSigmaPi,
                  sigma0PhotonExtra::PhotonNegTPCNSigmaPi,
                  sigma0PhotonExtra::PhotonPosTPCCrossedRows,
                  sigma0PhotonExtra::PhotonNegTPCCrossedRows,
                  sigma0PhotonExtra::PhotonPosPt,
                  sigma0PhotonExtra::PhotonNegPt,
                  sigma0PhotonExtra::PhotonPosEta,
                  sigma0PhotonExtra::PhotonNegEta,
                  sigma0PhotonExtra::PhotonPosY,
                  sigma0PhotonExtra::PhotonNegY,
                  sigma0PhotonExtra::PhotonPsiPair,
                  sigma0PhotonExtra::PhotonPosITSCls,
                  sigma0PhotonExtra::PhotonNegITSCls,
                  sigma0PhotonExtra::PhotonPosITSChi2PerNcl,
                  sigma0PhotonExtra::PhotonNegITSChi2PerNcl,
                  sigma0PhotonExtra::PhotonPosTrackCode,
                  sigma0PhotonExtra::PhotonNegTrackCode,
                  sigma0PhotonExtra::PhotonV0Type,
                  sigma0PhotonExtra::GammaBDTScore);

// For Lambda extra info
namespace sigma0LambdaExtra
{
DECLARE_SOA_COLUMN(LambdaPt, lambdaPt, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(AntiLambdaMass, antilambdaMass, float);
DECLARE_SOA_COLUMN(LambdaQt, lambdaQt, float);
DECLARE_SOA_COLUMN(LambdaAlpha, lambdaAlpha, float);
DECLARE_SOA_COLUMN(LambdaLifeTime, lambdaLifeTime, float);
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
DECLARE_SOA_COLUMN(LambdaPosITSChi2PerNcl, lambdaPosChi2PerNcl, float);
DECLARE_SOA_COLUMN(LambdaNegITSChi2PerNcl, lambdaNegChi2PerNcl, float);
DECLARE_SOA_COLUMN(LambdaPosTrackCode, lambdaPosTrackCode, uint8_t);
DECLARE_SOA_COLUMN(LambdaNegTrackCode, lambdaNegTrackCode, uint8_t);
DECLARE_SOA_COLUMN(LambdaV0Type, lambdaV0Type, uint8_t);
DECLARE_SOA_COLUMN(LambdaBDTScore, lambdaBDTScore, float);
DECLARE_SOA_COLUMN(AntiLambdaBDTScore, antilambdaBDTScore, float);

} // namespace sigma0LambdaExtra

DECLARE_SOA_TABLE(Sigma0LambdaExtras, "AOD", "SIGMA0LAMBDA",
                  sigma0LambdaExtra::LambdaPt,
                  sigma0LambdaExtra::LambdaMass,
                  sigma0LambdaExtra::AntiLambdaMass,
                  sigma0LambdaExtra::LambdaQt,
                  sigma0LambdaExtra::LambdaAlpha,
                  sigma0LambdaExtra::LambdaLifeTime,
                  sigma0LambdaExtra::LambdaRadius,
                  sigma0LambdaExtra::LambdaCosPA,
                  sigma0LambdaExtra::LambdaDCADau,
                  sigma0LambdaExtra::LambdaDCANegPV,
                  sigma0LambdaExtra::LambdaDCAPosPV,
                  sigma0LambdaExtra::LambdaEta,
                  sigma0LambdaExtra::LambdaY,
                  sigma0LambdaExtra::LambdaPhi,
                  sigma0LambdaExtra::LambdaPosPrTPCNSigma,
                  sigma0LambdaExtra::LambdaPosPiTPCNSigma,
                  sigma0LambdaExtra::LambdaNegPrTPCNSigma,
                  sigma0LambdaExtra::LambdaNegPiTPCNSigma,
                  sigma0LambdaExtra::LambdaPrTOFNSigma,
                  sigma0LambdaExtra::LambdaPiTOFNSigma,
                  sigma0LambdaExtra::ALambdaPrTOFNSigma,
                  sigma0LambdaExtra::ALambdaPiTOFNSigma,
                  sigma0LambdaExtra::LambdaPosTPCCrossedRows,
                  sigma0LambdaExtra::LambdaNegTPCCrossedRows,
                  sigma0LambdaExtra::LambdaPosPt,
                  sigma0LambdaExtra::LambdaNegPt,
                  sigma0LambdaExtra::LambdaPosEta,
                  sigma0LambdaExtra::LambdaNegEta,
                  sigma0LambdaExtra::LambdaPosPrY,
                  sigma0LambdaExtra::LambdaPosPiY,
                  sigma0LambdaExtra::LambdaNegPrY,
                  sigma0LambdaExtra::LambdaNegPiY,
                  sigma0LambdaExtra::LambdaPosITSCls,
                  sigma0LambdaExtra::LambdaNegITSCls,
                  sigma0LambdaExtra::LambdaPosITSChi2PerNcl,
                  sigma0LambdaExtra::LambdaNegITSChi2PerNcl,
                  sigma0LambdaExtra::LambdaPosTrackCode,
                  sigma0LambdaExtra::LambdaNegTrackCode,
                  sigma0LambdaExtra::LambdaV0Type,
                  sigma0LambdaExtra::LambdaBDTScore,
                  sigma0LambdaExtra::AntiLambdaBDTScore);

// for MC 
namespace sigma0MCCore
{
DECLARE_SOA_COLUMN(IsSigma, isSigma, bool); 
DECLARE_SOA_COLUMN(IsAntiSigma, isAntiSigma, bool);
DECLARE_SOA_COLUMN(SigmaMCPt, sigmaMCPt, float);
DECLARE_SOA_COLUMN(PhotonCandPDGCode, photonCandPDGCode, int);
DECLARE_SOA_COLUMN(PhotonCandPDGCodeMother, photonCandPDGCodeMother, int);
DECLARE_SOA_COLUMN(IsPhotonCandPrimary, isPhotonCandPrimary, bool);
DECLARE_SOA_COLUMN(PhotonMCPt, photonMCPt, float);
DECLARE_SOA_COLUMN(PhotonIsCorrectlyAssoc, photonIsCorrectlyAssoc, bool);
DECLARE_SOA_COLUMN(LambdaCandPDGCode, lambdaCandPDGCode, int);
DECLARE_SOA_COLUMN(LambdaCandPDGCodeMother, lambdaCandPDGCodeMother, int);
DECLARE_SOA_COLUMN(IsLambdaCandPrimary, isLambdaCandPrimary, bool);
DECLARE_SOA_COLUMN(LambdaMCPt, lambdaMCPt, float);
DECLARE_SOA_COLUMN(LambdaIsCorrectlyAssoc, lambdaIsCorrectlyAssoc, bool);

} // namespace sigma0MCCore

DECLARE_SOA_TABLE(Sigma0MCCores, "AOD", "SIGMA0MCCORES",
                  sigma0MCCore::IsSigma,
                  sigma0MCCore::IsAntiSigma,
                  sigma0MCCore::SigmaMCPt,
                  sigma0MCCore::PhotonCandPDGCode,
                  sigma0MCCore::PhotonCandPDGCodeMother,
                  sigma0MCCore::IsPhotonCandPrimary,
                  sigma0MCCore::PhotonMCPt,
                  sigma0MCCore::PhotonIsCorrectlyAssoc,
                  sigma0MCCore::LambdaCandPDGCode,
                  sigma0MCCore::LambdaCandPDGCodeMother,
                  sigma0MCCore::IsLambdaCandPrimary,
                  sigma0MCCore::LambdaMCPt,
                  sigma0MCCore::LambdaIsCorrectlyAssoc);

namespace sigma0Gen
{
DECLARE_SOA_COLUMN(IsSigma0, isSigma0, bool);           // true: sigma0, false: antisigma0                   
DECLARE_SOA_COLUMN(Sigma0MCPt, sigma0MCPt, float);      // MC pT
DECLARE_SOA_COLUMN(Sigma0Type, sigma0Type, int);        // Decay channel
} // namespace sigma0Gen

DECLARE_SOA_TABLE(Sigma0Gens, "AOD", "SIGMA0GENS",
                  sigma0Gen::IsSigma0,    
                  sigma0Gen::Sigma0MCPt,  
                  sigma0Gen::Sigma0Type); 

DECLARE_SOA_TABLE(SigmaCollRef, "AOD", "SIGMACOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraCollisionId);

DECLARE_SOA_TABLE(SigmaGenCollRef, "AOD", "SIGMAGENCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraMCCollisionId);

// ___________________________________________________________________________
// pi0 QA
namespace Pi0Core
{
DECLARE_SOA_COLUMN(Pi0Pt, pi0Pt, float);
DECLARE_SOA_COLUMN(Pi0Mass, pi0Mass, float);
DECLARE_SOA_COLUMN(Pi0Y, pi0Y, float);

DECLARE_SOA_COLUMN(Photon1Mass, photon1Mass, float);
DECLARE_SOA_COLUMN(Photon1Pt, photon1Pt, float);
DECLARE_SOA_COLUMN(Photon1Qt, photon1Qt, float);
DECLARE_SOA_COLUMN(Photon1Alpha, photon1Alpha, float);
DECLARE_SOA_COLUMN(Photon1DCAPosPV, photon1DCAPosPV, float);
DECLARE_SOA_COLUMN(Photon1DCANegPV, photon1DCANegPV, float);
DECLARE_SOA_COLUMN(Photon1DCADau, photon1DCADau, float);
DECLARE_SOA_COLUMN(Photon1NegEta, photon1NegEta, float);
DECLARE_SOA_COLUMN(Photon1PosEta, photon1PosEta, float);
DECLARE_SOA_COLUMN(Photon1CosPA, photon1CosPA, float);
DECLARE_SOA_COLUMN(Photon1Radius, photon1Radius, float);
DECLARE_SOA_COLUMN(Photon1Zconv, photon1Zconv, float);
DECLARE_SOA_COLUMN(Photon1PosTPCCrossedRows, photon1PosTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(Photon1NegTPCCrossedRows, photon1NegTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(Photon1PosTPCNSigmaEl, photon1PosTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(Photon1NegTPCNSigmaEl, photon1NegTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(Photon1Y, photon1Y, float);
DECLARE_SOA_COLUMN(Photon1V0Type, photon1V0Type, uint8_t);

DECLARE_SOA_COLUMN(Photon2Mass, photon2Mass, float);
DECLARE_SOA_COLUMN(Photon2Pt, photon2Pt, float);
DECLARE_SOA_COLUMN(Photon2Qt, photon2Qt, float);
DECLARE_SOA_COLUMN(Photon2Alpha, photon2Alpha, float);
DECLARE_SOA_COLUMN(Photon2DCAPosPV, photon2DCAPosPV, float);
DECLARE_SOA_COLUMN(Photon2DCANegPV, photon2DCANegPV, float);
DECLARE_SOA_COLUMN(Photon2DCADau, photon2DCADau, float);
DECLARE_SOA_COLUMN(Photon2NegEta, photon2NegEta, float);
DECLARE_SOA_COLUMN(Photon2PosEta, photon2PosEta, float);
DECLARE_SOA_COLUMN(Photon2CosPA, photon2CosPA, float);
DECLARE_SOA_COLUMN(Photon2Radius, photon2Radius, float);
DECLARE_SOA_COLUMN(Photon2Zconv, photon2Zconv, float);
DECLARE_SOA_COLUMN(Photon2PosTPCCrossedRows, photon2PosTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(Photon2NegTPCCrossedRows, photon2NegTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(Photon2PosTPCNSigmaEl, photon2PosTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(Photon2NegTPCNSigmaEl, photon2NegTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(Photon2Y, photon2Y, float);
DECLARE_SOA_COLUMN(Photon2V0Type, photon2V0Type, uint8_t);


} // namespace Pi0Core

DECLARE_SOA_TABLE(Pi0Cores, "AOD", "PI0CORES",
                  Pi0Core::Pi0Pt,
                  Pi0Core::Pi0Mass,
                  Pi0Core::Pi0Y,
                  Pi0Core::Photon1Mass,
                  Pi0Core::Photon1Pt,                  
                  Pi0Core::Photon1Qt,
                  Pi0Core::Photon1Alpha,
                  Pi0Core::Photon1DCAPosPV,
                  Pi0Core::Photon1DCANegPV,
                  Pi0Core::Photon1DCADau,
                  Pi0Core::Photon1NegEta,
                  Pi0Core::Photon1PosEta,
                  Pi0Core::Photon1CosPA,                  
                  Pi0Core::Photon1Radius,                                          
                  Pi0Core::Photon1Zconv,
                  Pi0Core::Photon1PosTPCCrossedRows,
                  Pi0Core::Photon1NegTPCCrossedRows,
                  Pi0Core::Photon1PosTPCNSigmaEl,
                  Pi0Core::Photon1NegTPCNSigmaEl,   
                  Pi0Core::Photon1Y,                                                                        
                  Pi0Core::Photon1V0Type,  
                  Pi0Core::Photon2Mass,
                  Pi0Core::Photon2Pt,                  
                  Pi0Core::Photon2Qt,
                  Pi0Core::Photon2Alpha,
                  Pi0Core::Photon2DCAPosPV,
                  Pi0Core::Photon2DCANegPV,
                  Pi0Core::Photon2DCADau,
                  Pi0Core::Photon2NegEta,
                  Pi0Core::Photon2PosEta,
                  Pi0Core::Photon2CosPA,                  
                  Pi0Core::Photon2Radius,                                          
                  Pi0Core::Photon2Zconv,
                  Pi0Core::Photon2PosTPCCrossedRows,
                  Pi0Core::Photon2NegTPCCrossedRows,
                  Pi0Core::Photon2PosTPCNSigmaEl,
                  Pi0Core::Photon2NegTPCNSigmaEl,   
                  Pi0Core::Photon2Y,                                                                        
                  Pi0Core::Photon2V0Type);

// for MC 
namespace Pi0CoreMC
{
DECLARE_SOA_COLUMN(IsPi0, isPi0, bool); 
DECLARE_SOA_COLUMN(Pi0MCPt, pi0MCPt, float);
DECLARE_SOA_COLUMN(Photon1MCPt, photon1MCPt, float);
DECLARE_SOA_COLUMN(Photon1PDGCode, photon1PDGCode, int);
DECLARE_SOA_COLUMN(Photon1PDGCodeMother, photon1PDGCodeMother, int);
DECLARE_SOA_COLUMN(IsPhoton1Primary, isPhoton1Primary, bool);
DECLARE_SOA_COLUMN(Photon1IsCorrectlyAssoc, photon1IsCorrectlyAssoc, bool);
DECLARE_SOA_COLUMN(Photon2MCPt, photon2MCPt, float);
DECLARE_SOA_COLUMN(Photon2PDGCode, photon2PDGCode, int);
DECLARE_SOA_COLUMN(Photon2PDGCodeMother, photon2PDGCodeMother, int);
DECLARE_SOA_COLUMN(IsPhoton2Primary, isPhoton2Primary, bool);
DECLARE_SOA_COLUMN(Photon2IsCorrectlyAssoc, photon2IsCorrectlyAssoc, bool);

} // namespace Pi0CoreMC

DECLARE_SOA_TABLE(Pi0CoresMC, "AOD", "PI0CORESMC",
                  Pi0CoreMC::IsPi0,                  
                  Pi0CoreMC::Pi0MCPt,
                  Pi0CoreMC::Photon1MCPt,
                  Pi0CoreMC::Photon1PDGCode,
                  Pi0CoreMC::Photon1PDGCodeMother,
                  Pi0CoreMC::IsPhoton1Primary,                  
                  Pi0CoreMC::Photon1IsCorrectlyAssoc,
                  Pi0CoreMC::Photon2MCPt,
                  Pi0CoreMC::Photon2PDGCode,
                  Pi0CoreMC::Photon2PDGCodeMother,
                  Pi0CoreMC::IsPhoton2Primary,                  
                  Pi0CoreMC::Photon2IsCorrectlyAssoc);


DECLARE_SOA_TABLE(Pi0CollRef, "AOD", "PI0COLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraCollisionId);

namespace pi0Gen
{
DECLARE_SOA_COLUMN(Pi0MCPt, pi0MCPt, float);      // MC pT
} // namespace pi0Gen

DECLARE_SOA_TABLE(Pi0Gens, "AOD", "PI0GENS",                     
                  pi0Gen::Pi0MCPt); 

DECLARE_SOA_TABLE(Pi0GenCollRef, "AOD", "PI0GENCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraMCCollisionId);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMATABLES_H_
