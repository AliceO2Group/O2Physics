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

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#include "Math/Vector3D.h"
#include "TVector3.h"

#include <cmath>
#include <vector>

#ifndef PWGLF_DATAMODEL_LFSIGMATABLES_H_
#define PWGLF_DATAMODEL_LFSIGMATABLES_H_

using std::array;

// Creating output TTree for sigma analysis
namespace o2::aod
{

// Indexing
namespace sigma0Core
{
DECLARE_SOA_INDEX_COLUMN_FULL(PhotonV0, photonV0, int, V0Cores, "_PhotonV0"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(LambdaV0, lambdaV0, int, V0Cores, "_LambdaV0"); //!
} // namespace sigma0Core

// for real data
namespace sigma0Core
{
DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(DCADaughters, dcadaughters, float);

DECLARE_SOA_COLUMN(PhotonPx, photonPx, float);
DECLARE_SOA_COLUMN(PhotonPy, photonPy, float);
DECLARE_SOA_COLUMN(PhotonPz, photonPz, float);
DECLARE_SOA_COLUMN(PhotonMass, photonMass, float);

DECLARE_SOA_COLUMN(LambdaPx, lambdaPx, float);
DECLARE_SOA_COLUMN(LambdaPy, lambdaPy, float);
DECLARE_SOA_COLUMN(LambdaPz, lambdaPz, float);
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);
DECLARE_SOA_COLUMN(AntiLambdaMass, antilambdaMass, float);

//______________________________________________________
// DYNAMIC COLUMNS
// Sigma0
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Sigma0 px
                           [](float photonPx, float lambdaPx) -> float { return photonPx + lambdaPx; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Sigma0 py
                           [](float photonPy, float lambdaPy) -> float { return photonPy + lambdaPy; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Sigma0 pz
                           [](float photonPz, float lambdaPz) -> float { return photonPz + lambdaPz; });

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](float photonPx, float photonPy, float lambdaPx, float lambdaPy) -> float {
                             return RecoDecay::pt(array{photonPx + lambdaPx, photonPy + lambdaPy});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Total momentum in GeV/c
                           [](float photonPx, float photonPy, float photonPz, float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photonPx + lambdaPx, photonPy + lambdaPy, photonPz + lambdaPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Sigma0Mass, sigma0Mass,
                           [](float photonPx, float photonPy, float photonPz, float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             std::array<float, 3> pVecPhotons{photonPx, photonPy, photonPz};
                             std::array<float, 3> pVecLambda{lambdaPx, lambdaPy, lambdaPz};
                             auto arrMom = std::array{pVecPhotons, pVecLambda};
                             return RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Sigma0Y, sigma0Y,
                           [](float photonPx, float photonPy, float photonPz, float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             return RecoDecay::y(std::array{photonPx + lambdaPx, photonPy + lambdaPy, photonPz + lambdaPz}, o2::constants::physics::MassSigma0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! Phi in the range [0, 2pi)
                           [](float photonPx, float photonPy, float lambdaPx, float lambdaPy) -> float { return RecoDecay::phi(photonPx + lambdaPx, photonPy + lambdaPy); });

DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! Pseudorapidity
                           [](float photonPx, float photonPy, float photonPz, float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             return RecoDecay::eta(std::array{photonPx + lambdaPx, photonPy + lambdaPy, photonPz + lambdaPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Radius, radius, //! Sigma0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

DECLARE_SOA_DYNAMIC_COLUMN(OPAngle, opAngle,
                           [](float photonPx, float photonPy, float photonPz, float lambdaPx, float lambdaPy, float lambdaPz) {
                             TVector3 v1(photonPx, photonPy, photonPz);
                             TVector3 v2(lambdaPx, lambdaPy, lambdaPz);
                             return v1.Angle(v2);
                           });

// Photon
DECLARE_SOA_DYNAMIC_COLUMN(PhotonPt, photonPt, //! Transverse momentum in GeV/c
                           [](float photonPx, float photonPy) -> float {
                             return RecoDecay::sqrtSumOfSquares(photonPx, photonPy);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonP, photonp, //! Total momentum in GeV/c
                           [](float photonPx, float photonPy, float photonPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photonPx, photonPy, photonPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonEta, photonEta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float photonPx, float photonPy, float photonPz) -> float {
                             return RecoDecay::eta(std::array{photonPx, photonPy, photonPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonY, photonY, //! Rapidity
                           [](float photonPx, float photonPy, float photonPz) -> float {
                             return RecoDecay::y(std::array{photonPx, photonPy, photonPz}, o2::constants::physics::MassGamma);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonPhi, photonPhi, //! Phi in the range [0, 2pi)
                           [](float photonPx, float photonPy) -> float { return RecoDecay::phi(photonPx, photonPy); });

// Lambda/ALambda
DECLARE_SOA_DYNAMIC_COLUMN(LambdaPt, lambdaPt, //! Transverse momentum in GeV/c
                           [](float lambdaPx, float lambdaPy) -> float {
                             return RecoDecay::sqrtSumOfSquares(lambdaPx, lambdaPy);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaP, lambdap, //! Total momentum in GeV/c
                           [](float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(lambdaPx, lambdaPy, lambdaPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaEta, lambdaEta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             return RecoDecay::eta(std::array{lambdaPx, lambdaPy, lambdaPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaY, lambdaY, //! Rapidity
                           [](float lambdaPx, float lambdaPy, float lambdaPz) -> float {
                             return RecoDecay::y(std::array{lambdaPx, lambdaPy, lambdaPz}, o2::constants::physics::MassLambda);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaPhi, lambdaPhi, //! Phi in the range [0, 2pi)
                           [](float lambdaPx, float lambdaPy) -> float { return RecoDecay::phi(lambdaPx, lambdaPy); });

} // namespace sigma0Core

DECLARE_SOA_TABLE(Sigma0Cores, "AOD", "SIGMA0CORES",
                  // Basic properties
                  sigma0Core::X, sigma0Core::Y, sigma0Core::Z, sigma0Core::DCADaughters,
                  sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz, sigma0Core::PhotonMass,
                  sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz, sigma0Core::LambdaMass, sigma0Core::AntiLambdaMass,

                  // Dynamic columns
                  sigma0Core::Px<sigma0Core::PhotonPx, sigma0Core::LambdaPx>,
                  sigma0Core::Py<sigma0Core::PhotonPy, sigma0Core::LambdaPy>,
                  sigma0Core::Pz<sigma0Core::PhotonPz, sigma0Core::LambdaPz>,
                  sigma0Core::Pt<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::LambdaPx, sigma0Core::LambdaPy>,
                  sigma0Core::P<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz, sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::Sigma0Mass<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz, sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::Sigma0Y<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz, sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::Phi<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::LambdaPx, sigma0Core::LambdaPy>,
                  sigma0Core::Eta<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz, sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::Radius<sigma0Core::X, sigma0Core::Y>,
                  sigma0Core::OPAngle<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz, sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,

                  sigma0Core::PhotonPt<sigma0Core::PhotonPx, sigma0Core::PhotonPy>,
                  sigma0Core::PhotonP<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz>,
                  sigma0Core::PhotonEta<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz>,
                  sigma0Core::PhotonY<sigma0Core::PhotonPx, sigma0Core::PhotonPy, sigma0Core::PhotonPz>,
                  sigma0Core::PhotonPhi<sigma0Core::PhotonPx, sigma0Core::PhotonPy>,

                  sigma0Core::LambdaPt<sigma0Core::LambdaPx, sigma0Core::LambdaPy>,
                  sigma0Core::LambdaP<sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::LambdaEta<sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::LambdaY<sigma0Core::LambdaPx, sigma0Core::LambdaPy, sigma0Core::LambdaPz>,
                  sigma0Core::LambdaPhi<sigma0Core::LambdaPx, sigma0Core::LambdaPy>);

// For Photon extra info
namespace sigma0PhotonExtra
{
//______________________________________________________
// REGULAR COLUMNS FOR SIGMA0PHOTON
DECLARE_SOA_COLUMN(PhotonQt, photonQt, float);
DECLARE_SOA_COLUMN(PhotonAlpha, photonAlpha, float);
DECLARE_SOA_COLUMN(PhotonCosPA, photonCosPA, float);
DECLARE_SOA_COLUMN(PhotonDCADau, photonDCADau, float);
DECLARE_SOA_COLUMN(PhotonDCANegPV, photonDCANegPV, float);
DECLARE_SOA_COLUMN(PhotonDCAPosPV, photonDCAPosPV, float);
DECLARE_SOA_COLUMN(PhotonRadius, photonRadius, float);
DECLARE_SOA_COLUMN(PhotonZconv, photonZconv, float);
DECLARE_SOA_COLUMN(PhotonPosTPCNSigmaEl, photonPosTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(PhotonNegTPCNSigmaEl, photonNegTPCNSigmaEl, float);
DECLARE_SOA_COLUMN(PhotonPosTPCCrossedRows, photonPosTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(PhotonNegTPCCrossedRows, photonNegTPCCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(PhotonPosEta, photonPosEta, float);
DECLARE_SOA_COLUMN(PhotonNegEta, photonNegEta, float);
DECLARE_SOA_COLUMN(PhotonPsiPair, photonPsiPair, float);
DECLARE_SOA_COLUMN(PhotonPosITSCls, photonPosITSCls, int);
DECLARE_SOA_COLUMN(PhotonNegITSCls, photonNegITSCls, int);
DECLARE_SOA_COLUMN(PhotonPosITSChi2PerNcl, photonPosITSChi2PerNcl, float);
DECLARE_SOA_COLUMN(PhotonNegITSChi2PerNcl, photonNegITSChi2PerNcl, float);
DECLARE_SOA_COLUMN(PhotonPosTrackCode, photonPosTrackCode, uint8_t);
DECLARE_SOA_COLUMN(PhotonNegTrackCode, photonNegTrackCode, uint8_t);
DECLARE_SOA_COLUMN(PhotonV0Type, photonV0Type, uint8_t);

} // namespace sigma0PhotonExtra

DECLARE_SOA_TABLE(Sigma0PhotonExtras, "AOD", "SIGMA0PHOTON",
                  sigma0PhotonExtra::PhotonQt,
                  sigma0PhotonExtra::PhotonAlpha,
                  sigma0PhotonExtra::PhotonCosPA,
                  sigma0PhotonExtra::PhotonDCADau,
                  sigma0PhotonExtra::PhotonDCANegPV,
                  sigma0PhotonExtra::PhotonDCAPosPV,
                  sigma0PhotonExtra::PhotonRadius,
                  sigma0PhotonExtra::PhotonZconv,
                  sigma0PhotonExtra::PhotonPosTPCNSigmaEl,
                  sigma0PhotonExtra::PhotonNegTPCNSigmaEl,
                  sigma0PhotonExtra::PhotonPosTPCCrossedRows,
                  sigma0PhotonExtra::PhotonNegTPCCrossedRows,
                  sigma0PhotonExtra::PhotonPosEta,
                  sigma0PhotonExtra::PhotonNegEta,
                  sigma0PhotonExtra::PhotonPsiPair,
                  sigma0PhotonExtra::PhotonPosITSCls,
                  sigma0PhotonExtra::PhotonNegITSCls,
                  sigma0PhotonExtra::PhotonPosITSChi2PerNcl,
                  sigma0PhotonExtra::PhotonNegITSChi2PerNcl,
                  sigma0PhotonExtra::PhotonPosTrackCode,
                  sigma0PhotonExtra::PhotonNegTrackCode,
                  sigma0PhotonExtra::PhotonV0Type);

// For Lambda extra info
namespace sigma0LambdaExtra
{
DECLARE_SOA_COLUMN(LambdaQt, lambdaQt, float);
DECLARE_SOA_COLUMN(LambdaAlpha, lambdaAlpha, float);
DECLARE_SOA_COLUMN(LambdaLifeTime, lambdaLifeTime, float);
DECLARE_SOA_COLUMN(LambdaRadius, lambdaRadius, float);
DECLARE_SOA_COLUMN(LambdaCosPA, lambdaCosPA, float);
DECLARE_SOA_COLUMN(LambdaDCADau, lambdaDCADau, float);
DECLARE_SOA_COLUMN(LambdaDCANegPV, lambdaDCANegPV, float);
DECLARE_SOA_COLUMN(LambdaDCAPosPV, lambdaDCAPosPV, float);
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
DECLARE_SOA_COLUMN(LambdaPosEta, lambdaPosEta, float);
DECLARE_SOA_COLUMN(LambdaNegEta, lambdaNegEta, float);
DECLARE_SOA_COLUMN(LambdaPosITSCls, lambdaPosITSCls, int);
DECLARE_SOA_COLUMN(LambdaNegITSCls, lambdaNegITSCls, int);
DECLARE_SOA_COLUMN(LambdaPosITSChi2PerNcl, lambdaPosChi2PerNcl, float);
DECLARE_SOA_COLUMN(LambdaNegITSChi2PerNcl, lambdaNegChi2PerNcl, float);
DECLARE_SOA_COLUMN(LambdaPosTrackCode, lambdaPosTrackCode, uint8_t);
DECLARE_SOA_COLUMN(LambdaNegTrackCode, lambdaNegTrackCode, uint8_t);
DECLARE_SOA_COLUMN(LambdaV0Type, lambdaV0Type, uint8_t);

} // namespace sigma0LambdaExtra

DECLARE_SOA_TABLE(Sigma0LambdaExtras, "AOD", "SIGMA0LAMBDA",
                  sigma0LambdaExtra::LambdaQt,
                  sigma0LambdaExtra::LambdaAlpha,
                  sigma0LambdaExtra::LambdaLifeTime,
                  sigma0LambdaExtra::LambdaRadius,
                  sigma0LambdaExtra::LambdaCosPA,
                  sigma0LambdaExtra::LambdaDCADau,
                  sigma0LambdaExtra::LambdaDCANegPV,
                  sigma0LambdaExtra::LambdaDCAPosPV,
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
                  sigma0LambdaExtra::LambdaPosEta,
                  sigma0LambdaExtra::LambdaNegEta,
                  sigma0LambdaExtra::LambdaPosITSCls,
                  sigma0LambdaExtra::LambdaNegITSCls,
                  sigma0LambdaExtra::LambdaPosITSChi2PerNcl,
                  sigma0LambdaExtra::LambdaNegITSChi2PerNcl,
                  sigma0LambdaExtra::LambdaPosTrackCode,
                  sigma0LambdaExtra::LambdaNegTrackCode,
                  sigma0LambdaExtra::LambdaV0Type);

// for MC
namespace sigma0MCCore
{
DECLARE_SOA_COLUMN(MCradius, mcradius, float);
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);
DECLARE_SOA_COLUMN(MCprocess, mcprocess, int);
DECLARE_SOA_COLUMN(IsProducedByGenerator, isProducedByGenerator, bool);

DECLARE_SOA_COLUMN(PhotonMCPx, photonmcpx, float);
DECLARE_SOA_COLUMN(PhotonMCPy, photonmcpy, float);
DECLARE_SOA_COLUMN(PhotonMCPz, photonmcpz, float);
DECLARE_SOA_COLUMN(IsPhotonPrimary, isPhotonPrimary, bool);
DECLARE_SOA_COLUMN(PhotonPDGCode, photonPDGCode, int);
DECLARE_SOA_COLUMN(PhotonPDGCodeMother, photonPDGCodeMother, int);
DECLARE_SOA_COLUMN(PhotonIsCorrectlyAssoc, photonIsCorrectlyAssoc, bool);

DECLARE_SOA_COLUMN(LambdaMCPx, lambdamcpx, float);
DECLARE_SOA_COLUMN(LambdaMCPy, lambdamcpy, float);
DECLARE_SOA_COLUMN(LambdaMCPz, lambdamcpz, float);
DECLARE_SOA_COLUMN(IsLambdaPrimary, isLambdaPrimary, bool);
DECLARE_SOA_COLUMN(LambdaPDGCode, lambdaPDGCode, int);
DECLARE_SOA_COLUMN(LambdaPDGCodeMother, lambdaPDGCodeMother, int);
DECLARE_SOA_COLUMN(LambdaIsCorrectlyAssoc, lambdaIsCorrectlyAssoc, bool);

DECLARE_SOA_DYNAMIC_COLUMN(IsSigma0, isSigma0, //! IsSigma0
                           [](int pdgCode) -> bool { return pdgCode == 3212; });

DECLARE_SOA_DYNAMIC_COLUMN(IsAntiSigma0, isAntiSigma0, //! IsASigma0
                           [](int pdgCode) -> bool { return pdgCode == -3212; });

DECLARE_SOA_DYNAMIC_COLUMN(MCPx, mcpx, //! Sigma0 px
                           [](float photonMCPx, float lambdaMCPx) -> float { return photonMCPx + lambdaMCPx; });
DECLARE_SOA_DYNAMIC_COLUMN(MCPy, mcpy, //! Sigma0 py
                           [](float photonMCPy, float lambdaMCPy) -> float { return photonMCPy + lambdaMCPy; });
DECLARE_SOA_DYNAMIC_COLUMN(MCPz, mcpz, //! Sigma0 pz
                           [](float photonMCPz, float lambdaMCPz) -> float { return photonMCPz + lambdaMCPz; });

DECLARE_SOA_DYNAMIC_COLUMN(MCPt, mcpt,
                           [](float photonMCPx, float photonMCPy, float lambdaMCPx, float lambdaMCPy) -> float {
                             return RecoDecay::pt(array{photonMCPx + lambdaMCPx, photonMCPy + lambdaMCPy});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(MCP, mcp, //! Total momentum in GeV/c
                           [](float photonMCPx, float photonMCPy, float photonMCPz, float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photonMCPx + lambdaMCPx, photonMCPy + lambdaMCPy, photonMCPz + lambdaMCPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Sigma0MCMass, sigma0MCMass,
                           [](float photonMCPx, float photonMCPy, float photonMCPz, float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             std::array<float, 3> pVecPhotons{photonMCPx, photonMCPy, photonMCPz};
                             std::array<float, 3> pVecLambda{lambdaMCPx, lambdaMCPy, lambdaMCPz};
                             auto arrMom = std::array{pVecPhotons, pVecLambda};
                             return RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Sigma0MCY, sigma0MCY,
                           [](float photonMCPx, float photonMCPy, float photonMCPz, float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             return RecoDecay::y(std::array{photonMCPx + lambdaMCPx, photonMCPy + lambdaMCPy, photonMCPz + lambdaMCPz}, o2::constants::physics::MassSigma0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(MCPhi, mcphi, //! Phi in the range [0, 2pi)
                           [](float photonMCPx, float photonMCPy, float lambdaMCPx, float lambdaMCPy) -> float { return RecoDecay::phi(photonMCPx + lambdaMCPx, photonMCPy + lambdaMCPy); });

DECLARE_SOA_DYNAMIC_COLUMN(MCEta, mceta, //! Pseudorapidity
                           [](float photonMCPx, float photonMCPy, float photonMCPz, float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             return RecoDecay::eta(std::array{photonMCPx + lambdaMCPx, photonMCPy + lambdaMCPy, photonMCPz + lambdaMCPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(MCOPAngle, mcopAngle,
                           [](float photonMCPx, float photonMCPy, float photonMCPz, float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) {
                             TVector3 v1(photonMCPx, photonMCPy, photonMCPz);
                             TVector3 v2(lambdaMCPx, lambdaMCPy, lambdaMCPz);
                             return v1.Angle(v2);
                           });

// Photon
DECLARE_SOA_DYNAMIC_COLUMN(PhotonMCPt, photonmcpt, //! Transverse momentum in GeV/c
                           [](float photonMCPx, float photonMCPy) -> float {
                             return RecoDecay::sqrtSumOfSquares(photonMCPx, photonMCPy);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonMCP, photonmcp, //! Total momentum in GeV/c
                           [](float photonMCPx, float photonMCPy, float photonMCPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photonMCPx, photonMCPy, photonMCPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonMCEta, photonMCEta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float photonMCPx, float photonMCPy, float photonMCPz) -> float {
                             return RecoDecay::eta(std::array{photonMCPx, photonMCPy, photonMCPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonMCY, photonMCY, //! Rapidity
                           [](float photonMCPx, float photonMCPy, float photonMCPz) -> float {
                             return RecoDecay::y(std::array{photonMCPx, photonMCPy, photonMCPz}, o2::constants::physics::MassGamma);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PhotonMCPhi, photonMCPhi, //! Phi in the range [0, 2pi)
                           [](float photonMCPx, float photonMCPy) -> float { return RecoDecay::phi(photonMCPx, photonMCPy); });

// Lambda/ALambda
DECLARE_SOA_DYNAMIC_COLUMN(LambdaMCPt, lambdamcpt, //! Transverse momentum in GeV/c
                           [](float lambdaMCPx, float lambdaMCPy) -> float {
                             return RecoDecay::sqrtSumOfSquares(lambdaMCPx, lambdaMCPy);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaMCP, lambdamcp, //! Total momentum in GeV/c
                           [](float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(lambdaMCPx, lambdaMCPy, lambdaMCPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaMCEta, lambdaMCEta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             return RecoDecay::eta(std::array{lambdaMCPx, lambdaMCPy, lambdaMCPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaMCY, lambdaMCY, //! Rapidity
                           [](float lambdaMCPx, float lambdaMCPy, float lambdaMCPz) -> float {
                             return RecoDecay::y(std::array{lambdaMCPx, lambdaMCPy, lambdaMCPz}, o2::constants::physics::MassLambda);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(LambdaMCPhi, lambdaMCPhi, //! Phi in the range [0, 2pi)
                           [](float lambdaMCPx, float lambdaMCPy) -> float { return RecoDecay::phi(lambdaMCPx, lambdaMCPy); });

} // namespace sigma0MCCore

DECLARE_SOA_TABLE(Sigma0MCCores, "AOD", "SIGMA0MCCORES",

                  // Basic properties
                  sigma0MCCore::MCradius, sigma0MCCore::PDGCode, sigma0MCCore::PDGCodeMother, sigma0MCCore::MCprocess, sigma0MCCore::IsProducedByGenerator,

                  sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz,
                  sigma0MCCore::IsPhotonPrimary, sigma0MCCore::PhotonPDGCode, sigma0MCCore::PhotonPDGCodeMother, sigma0MCCore::PhotonIsCorrectlyAssoc,

                  sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz,
                  sigma0MCCore::IsLambdaPrimary, sigma0MCCore::LambdaPDGCode, sigma0MCCore::LambdaPDGCodeMother, sigma0MCCore::LambdaIsCorrectlyAssoc,

                  // Dynamic columns
                  sigma0MCCore::IsSigma0<sigma0MCCore::PDGCode>,
                  sigma0MCCore::IsAntiSigma0<sigma0MCCore::PDGCode>,

                  sigma0MCCore::MCPx<sigma0MCCore::PhotonMCPx, sigma0MCCore::LambdaMCPx>,
                  sigma0MCCore::MCPy<sigma0MCCore::PhotonMCPy, sigma0MCCore::LambdaMCPy>,
                  sigma0MCCore::MCPz<sigma0MCCore::PhotonMCPz, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::MCPt<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy>,
                  sigma0MCCore::MCP<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::Sigma0MCMass<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::Sigma0MCY<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::MCPhi<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy>,
                  sigma0MCCore::MCEta<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::MCOPAngle<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz, sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,

                  sigma0MCCore::PhotonMCPt<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy>,
                  sigma0MCCore::PhotonMCP<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz>,
                  sigma0MCCore::PhotonMCEta<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz>,
                  sigma0MCCore::PhotonMCY<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy, sigma0MCCore::PhotonMCPz>,
                  sigma0MCCore::PhotonMCPhi<sigma0MCCore::PhotonMCPx, sigma0MCCore::PhotonMCPy>,

                  sigma0MCCore::LambdaMCPt<sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy>,
                  sigma0MCCore::LambdaMCP<sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::LambdaMCEta<sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::LambdaMCY<sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy, sigma0MCCore::LambdaMCPz>,
                  sigma0MCCore::LambdaMCPhi<sigma0MCCore::LambdaMCPx, sigma0MCCore::LambdaMCPy>);

namespace sigma0MCCore
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for Sigma0
}
namespace sigma0Gen
{
DECLARE_SOA_COLUMN(IsSigma0, isSigma0, bool); // true: sigma0, false: antisigma0
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(MCPt, mcpt, float); // MC pT
DECLARE_SOA_COLUMN(MCY, mcy, float);   // MC Y

} // namespace sigma0Gen

DECLARE_SOA_TABLE(Sigma0Gens, "AOD", "SIGMA0GENS",
                  sigma0Gen::IsSigma0,
                  sigma0Gen::ProducedByGenerator,
                  sigma0Gen::MCPt,
                  sigma0Gen::MCY);

DECLARE_SOA_TABLE(SigmaCollRef, "AOD", "SIGMACOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraCollisionId);

DECLARE_SOA_TABLE(SigmaIndices, "AOD", "SIGMAINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, sigma0Core::PhotonV0Id, sigma0Core::LambdaV0Id, o2::soa::Marker<1>);

DECLARE_SOA_TABLE(SigmaMCLabels, "AOD", "SIGMAMCLABEL", //! optional table to refer to mcparticles
                  o2::soa::Index<>, sigma0MCCore::McParticleId);

DECLARE_SOA_TABLE(SigmaGenCollRef, "AOD", "SIGMAGENCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraMCCollisionId);

// ___________________________________________________________________________
// pi0 QA
namespace Pi0Core
{

DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(DCADaughters, dcadaughters, float);
DECLARE_SOA_COLUMN(CosPA, cospa, float);

DECLARE_SOA_COLUMN(Photon1Px, photon1Px, float);
DECLARE_SOA_COLUMN(Photon1Py, photon1Py, float);
DECLARE_SOA_COLUMN(Photon1Pz, photon1Pz, float);
DECLARE_SOA_COLUMN(Photon1Mass, photon1Mass, float);
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
DECLARE_SOA_COLUMN(Photon1V0Type, photon1V0Type, uint8_t);

DECLARE_SOA_COLUMN(Photon2Px, photon2Px, float);
DECLARE_SOA_COLUMN(Photon2Py, photon2Py, float);
DECLARE_SOA_COLUMN(Photon2Pz, photon2Pz, float);
DECLARE_SOA_COLUMN(Photon2Mass, photon2Mass, float);
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
DECLARE_SOA_COLUMN(Photon2V0Type, photon2V0Type, uint8_t);

//______________________________________________________
// DYNAMIC COLUMNS
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Pi0 px
                           [](float photon1Px, float photon2Px) -> float { return photon1Px + photon2Px; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Pi0 py
                           [](float photon1Py, float photon2Py) -> float { return photon1Py + photon2Py; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Pi0 pz
                           [](float photon1Pz, float photon2Pz) -> float { return photon1Pz + photon2Pz; });

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](float photon1Px, float photon1Py, float photon2Px, float photon2Py) -> float {
                             return RecoDecay::pt(array{photon1Px + photon2Px, photon1Py + photon2Py});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Total momentum in GeV/c
                           [](float photon1Px, float photon1Py, float photon1Pz, float photon2Px, float photon2Py, float photon2Pz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon1Px + photon2Px, photon1Py + photon2Py, photon1Pz + photon2Pz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Pi0Mass, pi0Mass,
                           [](float photon1Px, float photon1Py, float photon1Pz, float photon2Px, float photon2Py, float photon2Pz) -> float {
                             std::array<float, 3> pVecPhoton1{photon1Px, photon1Py, photon1Pz};
                             std::array<float, 3> pVecPhoton2{photon2Px, photon2Py, photon2Pz};
                             auto arrMom = std::array{pVecPhoton1, pVecPhoton2};
                             return RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Pi0Y, pi0Y,
                           [](float photon1Px, float photon1Py, float photon1Pz, float photon2Px, float photon2Py, float photon2Pz) -> float {
                             return RecoDecay::y(std::array{photon1Px + photon2Px, photon1Py + photon2Py, photon1Pz + photon2Pz}, o2::constants::physics::MassPi0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! Phi in the range [0, 2pi)
                           [](float photon1Px, float photon1Py, float photon2Px, float photon2Py) -> float { return RecoDecay::phi(photon1Px + photon2Px, photon1Py + photon2Py); });

DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! Pseudorapidity
                           [](float photon1Px, float photon1Py, float photon1Pz, float photon2Px, float photon2Py, float photon2Pz) -> float {
                             return RecoDecay::eta(std::array{photon1Px + photon2Px, photon1Py + photon2Py, photon1Pz + photon2Pz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Radius, radius, //! Pi0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

DECLARE_SOA_DYNAMIC_COLUMN(OPAngle, opAngle,
                           [](float photon1Px, float photon1Py, float photon1Pz, float photon2Px, float photon2Py, float photon2Pz) {
                             TVector3 v1(photon1Px, photon1Py, photon1Pz);
                             TVector3 v2(photon2Px, photon2Py, photon2Pz);
                             return v1.Angle(v2);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1Pt, photon1Pt, //! Transverse momentum in GeV/c
                           [](float photon1Px, float photon1Py) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon1Px, photon1Py);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1P, photon1p, //! Total momentum in GeV/c
                           [](float photon1Px, float photon1Py, float photon1Pz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon1Px, photon1Py, photon1Pz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1Eta, photon1Eta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float photon1Px, float photon1Py, float photon1Pz) -> float {
                             return RecoDecay::eta(std::array{photon1Px, photon1Py, photon1Pz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1Y, photon1Y, //! Rapidity
                           [](float photon1Px, float photon1Py, float photon1Pz) -> float {
                             return RecoDecay::y(std::array{photon1Px, photon1Py, photon1Pz}, o2::constants::physics::MassGamma);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1Phi, photon1Phi, //! Phi in the range [0, 2pi)
                           [](float photon1Px, float photon1Py) -> float { return RecoDecay::phi(photon1Px, photon1Py); });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2Pt, photon2Pt, //! Transverse momentum in GeV/c
                           [](float photon2Px, float photon2Py) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon2Px, photon2Py);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2P, photon2p, //! Total momentum in GeV/c
                           [](float photon2Px, float photon2Py, float photon2Pz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon2Px, photon2Py, photon2Pz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2Eta, photon2Eta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float photon2Px, float photon2Py, float photon2Pz) -> float {
                             return RecoDecay::eta(std::array{photon2Px, photon2Py, photon2Pz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2Y, photon2Y, //! Rapidity
                           [](float photon2Px, float photon2Py, float photon2Pz) -> float {
                             return RecoDecay::y(std::array{photon2Px, photon2Py, photon2Pz}, o2::constants::physics::MassGamma);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2Phi, photon2Phi, //! Phi in the range [0, 2pi)
                           [](float photon2Px, float photon2Py) -> float { return RecoDecay::phi(photon2Px, photon2Py); });

} // namespace Pi0Core

DECLARE_SOA_TABLE(Pi0Cores, "AOD", "PI0CORES",
                  Pi0Core::X, Pi0Core::Y, Pi0Core::Z, Pi0Core::DCADaughters, Pi0Core::CosPA,

                  // Photon 1 base properties
                  Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz,
                  Pi0Core::Photon1Mass, Pi0Core::Photon1Qt, Pi0Core::Photon1Alpha, Pi0Core::Photon1DCAPosPV, Pi0Core::Photon1DCANegPV, Pi0Core::Photon1DCADau,
                  Pi0Core::Photon1NegEta, Pi0Core::Photon1PosEta, Pi0Core::Photon1CosPA, Pi0Core::Photon1Radius, Pi0Core::Photon1Zconv,
                  Pi0Core::Photon1PosTPCCrossedRows, Pi0Core::Photon1NegTPCCrossedRows, Pi0Core::Photon1PosTPCNSigmaEl, Pi0Core::Photon1NegTPCNSigmaEl, Pi0Core::Photon1V0Type,

                  // Photon 2 base properties
                  Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz,
                  Pi0Core::Photon2Mass, Pi0Core::Photon2Qt, Pi0Core::Photon2Alpha, Pi0Core::Photon2DCAPosPV, Pi0Core::Photon2DCANegPV, Pi0Core::Photon2DCADau,
                  Pi0Core::Photon2NegEta, Pi0Core::Photon2PosEta, Pi0Core::Photon2CosPA, Pi0Core::Photon2Radius, Pi0Core::Photon2Zconv,
                  Pi0Core::Photon2PosTPCCrossedRows, Pi0Core::Photon2NegTPCCrossedRows, Pi0Core::Photon2PosTPCNSigmaEl, Pi0Core::Photon2NegTPCNSigmaEl, Pi0Core::Photon2V0Type,

                  // Dynamic columns
                  Pi0Core::Px<Pi0Core::Photon1Px, Pi0Core::Photon2Px>,
                  Pi0Core::Py<Pi0Core::Photon1Py, Pi0Core::Photon2Py>,
                  Pi0Core::Pz<Pi0Core::Photon1Pz, Pi0Core::Photon2Pz>,
                  Pi0Core::Pt<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon2Px, Pi0Core::Photon2Py>,
                  Pi0Core::P<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz, Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Pi0Mass<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz, Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Pi0Y<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz, Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Phi<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon2Px, Pi0Core::Photon2Py>,
                  Pi0Core::Eta<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz, Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Radius<Pi0Core::X, Pi0Core::Y>,
                  Pi0Core::OPAngle<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz, Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,

                  Pi0Core::Photon1Pt<Pi0Core::Photon1Px, Pi0Core::Photon1Py>,
                  Pi0Core::Photon1P<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz>,
                  Pi0Core::Photon1Eta<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz>,
                  Pi0Core::Photon1Y<Pi0Core::Photon1Px, Pi0Core::Photon1Py, Pi0Core::Photon1Pz>,
                  Pi0Core::Photon1Phi<Pi0Core::Photon1Px, Pi0Core::Photon1Py>,

                  Pi0Core::Photon2Pt<Pi0Core::Photon2Px, Pi0Core::Photon2Py>,
                  Pi0Core::Photon2P<Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Photon2Eta<Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Photon2Y<Pi0Core::Photon2Px, Pi0Core::Photon2Py, Pi0Core::Photon2Pz>,
                  Pi0Core::Photon2Phi<Pi0Core::Photon2Px, Pi0Core::Photon2Py>);

// for MC
namespace Pi0CoreMC
{

DECLARE_SOA_COLUMN(MCradius, mcradius, float);
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);
DECLARE_SOA_COLUMN(MCprocess, mcprocess, int);
DECLARE_SOA_COLUMN(IsProducedByGenerator, isProducedByGenerator, bool);

DECLARE_SOA_COLUMN(Photon1MCPx, photon1mcpx, float);
DECLARE_SOA_COLUMN(Photon1MCPy, photon1mcpy, float);
DECLARE_SOA_COLUMN(Photon1MCPz, photon1mcpz, float);
DECLARE_SOA_COLUMN(IsPhoton1Primary, isPhoton1Primary, bool);
DECLARE_SOA_COLUMN(Photon1PDGCode, photon1PDGCode, int);
DECLARE_SOA_COLUMN(Photon1PDGCodeMother, photon1PDGCodeMother, int);
DECLARE_SOA_COLUMN(Photon1IsCorrectlyAssoc, photon1IsCorrectlyAssoc, bool);

DECLARE_SOA_COLUMN(Photon2MCPx, photon2mcpx, float);
DECLARE_SOA_COLUMN(Photon2MCPy, photon2mcpy, float);
DECLARE_SOA_COLUMN(Photon2MCPz, photon2mcpz, float);
DECLARE_SOA_COLUMN(IsPhoton2Primary, isPhoton2Primary, bool);
DECLARE_SOA_COLUMN(Photon2PDGCode, photon2PDGCode, int);
DECLARE_SOA_COLUMN(Photon2PDGCodeMother, photon2PDGCodeMother, int);
DECLARE_SOA_COLUMN(Photon2IsCorrectlyAssoc, photon2IsCorrectlyAssoc, bool);

DECLARE_SOA_DYNAMIC_COLUMN(IsPi0, isPi0, //! IsPi0
                           [](int pdgCode) -> bool { return pdgCode == 111; });

DECLARE_SOA_DYNAMIC_COLUMN(IsFromXi0, isFromXi0, //! Pi0 from Xi0
                           [](int pdgCodeMother) -> bool { return pdgCodeMother == 3322; });

DECLARE_SOA_DYNAMIC_COLUMN(MCPx, mcpx, //! Pi0 MC px
                           [](float photon1MCPx, float photon2MCPx) -> float { return photon1MCPx + photon2MCPx; });
DECLARE_SOA_DYNAMIC_COLUMN(MCPy, mcpy, //! Pi0 MC py
                           [](float photon1MCPy, float photon2MCPy) -> float { return photon1MCPy + photon2MCPy; });
DECLARE_SOA_DYNAMIC_COLUMN(MCPz, mcpz, //! Pi0 MC pz
                           [](float photon1MCPz, float photon2MCPz) -> float { return photon1MCPz + photon2MCPz; });

DECLARE_SOA_DYNAMIC_COLUMN(MCPt, mcpt,
                           [](float photon1MCPx, float photon1MCPy, float photon2MCPx, float photon2MCPy) -> float {
                             return RecoDecay::pt(array{photon1MCPx + photon2MCPx, photon1MCPy + photon2MCPy});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(MCP, mcp, //! Total momentum in GeV/c
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz, float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon1MCPx + photon2MCPx, photon1MCPy + photon2MCPy, photon1MCPz + photon2MCPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Pi0MCMass, pi0MCMass,
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz, float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             std::array<float, 3> pVecPhoton1{photon1MCPx, photon1MCPy, photon1MCPz};
                             std::array<float, 3> pVecPhoton2{photon2MCPx, photon2MCPy, photon2MCPz};
                             auto arrMom = std::array{pVecPhoton1, pVecPhoton2};
                             return RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Pi0MCY, pi0MCY,
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz, float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             return RecoDecay::y(std::array{photon1MCPx + photon2MCPx, photon1MCPy + photon2MCPy, photon1MCPz + photon2MCPz}, o2::constants::physics::MassPi0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(MCPhi, mcphi, //! Phi in the range [0, 2pi)
                           [](float photon1MCPx, float photon1MCPy, float photon2MCPx, float photon2MCPy) -> float { return RecoDecay::phi(photon1MCPx + photon2MCPx, photon1MCPy + photon2MCPy); });

DECLARE_SOA_DYNAMIC_COLUMN(MCEta, mceta, //! Pseudorapidity
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz, float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             return RecoDecay::eta(std::array{photon1MCPx + photon2MCPx, photon1MCPy + photon2MCPy, photon1MCPz + photon2MCPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(MCOPAngle, mcopAngle,
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz, float photon2MCPx, float photon2MCPy, float photon2MCPz) {
                             TVector3 v1(photon1MCPx, photon1MCPy, photon1MCPz);
                             TVector3 v2(photon2MCPx, photon2MCPy, photon2MCPz);
                             return v1.Angle(v2);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1MCPt, photon1MCPt, //! Transverse momentum in GeV/c
                           [](float photon1MCPx, float photon1MCPy) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon1MCPx, photon1MCPy);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1MCP, photon1MCp, //! Total momentum in GeV/c
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon1MCPx, photon1MCPy, photon1MCPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1MCEta, photon1MCEta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz) -> float {
                             return RecoDecay::eta(std::array{photon1MCPx, photon1MCPy, photon1MCPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1MCY, photon1MCY, //! Rapidity
                           [](float photon1MCPx, float photon1MCPy, float photon1MCPz) -> float {
                             return RecoDecay::y(std::array{photon1MCPx, photon1MCPy, photon1MCPz}, o2::constants::physics::MassGamma);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon1MCPhi, photon1MCPhi, //! Phi in the range [0, 2pi)
                           [](float photon1MCPx, float photon1MCPy) -> float { return RecoDecay::phi(photon1MCPx, photon1MCPy); });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2MCPt, photon2MCPt, //! Transverse momentum in GeV/c
                           [](float photon2MCPx, float photon2MCPy) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon2MCPx, photon2MCPy);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2MCP, photon2MCp, //! Total momentum in GeV/c
                           [](float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             return RecoDecay::sqrtSumOfSquares(photon2MCPx, photon2MCPy, photon2MCPz);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2MCEta, photon2MCEta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             return RecoDecay::eta(std::array{photon2MCPx, photon2MCPy, photon2MCPz});
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2MCY, photon2MCY, //! Rapidity
                           [](float photon2MCPx, float photon2MCPy, float photon2MCPz) -> float {
                             return RecoDecay::y(std::array{photon2MCPx, photon2MCPy, photon2MCPz}, o2::constants::physics::MassGamma);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Photon2MCPhi, photon2MCPhi, //! Phi in the range [0, 2pi)
                           [](float photon2MCPx, float photon2MCPy) -> float { return RecoDecay::phi(photon2MCPx, photon2MCPy); });

} // namespace Pi0CoreMC

DECLARE_SOA_TABLE(Pi0CoresMC, "AOD", "PI0CORESMC",
                  // Basic properties
                  Pi0CoreMC::MCradius, Pi0CoreMC::PDGCode, Pi0CoreMC::PDGCodeMother, Pi0CoreMC::MCprocess, Pi0CoreMC::IsProducedByGenerator,

                  Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz,
                  Pi0CoreMC::IsPhoton1Primary, Pi0CoreMC::Photon1PDGCode, Pi0CoreMC::Photon1PDGCodeMother, Pi0CoreMC::Photon1IsCorrectlyAssoc,

                  Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz,
                  Pi0CoreMC::IsPhoton2Primary, Pi0CoreMC::Photon2PDGCode, Pi0CoreMC::Photon2PDGCodeMother, Pi0CoreMC::Photon2IsCorrectlyAssoc,

                  // Dynamic columns
                  Pi0CoreMC::IsPi0<Pi0CoreMC::PDGCode>,
                  Pi0CoreMC::IsFromXi0<Pi0CoreMC::PDGCodeMother>,

                  Pi0CoreMC::MCPx<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon2MCPx>,
                  Pi0CoreMC::MCPy<Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon2MCPy>,
                  Pi0CoreMC::MCPz<Pi0CoreMC::Photon1MCPz, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::MCPt<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy>,
                  Pi0CoreMC::MCP<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::Pi0MCMass<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::Pi0MCY<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::MCPhi<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy>,
                  Pi0CoreMC::MCEta<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::MCOPAngle<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz, Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,

                  Pi0CoreMC::Photon1MCPt<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy>,
                  Pi0CoreMC::Photon1MCP<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz>,
                  Pi0CoreMC::Photon1MCEta<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz>,
                  Pi0CoreMC::Photon1MCY<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy, Pi0CoreMC::Photon1MCPz>,
                  Pi0CoreMC::Photon1MCPhi<Pi0CoreMC::Photon1MCPx, Pi0CoreMC::Photon1MCPy>,

                  Pi0CoreMC::Photon2MCPt<Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy>,
                  Pi0CoreMC::Photon2MCP<Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::Photon2MCEta<Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::Photon2MCY<Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy, Pi0CoreMC::Photon2MCPz>,
                  Pi0CoreMC::Photon2MCPhi<Pi0CoreMC::Photon2MCPx, Pi0CoreMC::Photon2MCPy>);

DECLARE_SOA_TABLE(Pi0CollRef, "AOD", "PI0COLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraCollisionId);

namespace pi0Gen
{
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(MCPt, mcpt, float); // MC pT
DECLARE_SOA_COLUMN(MCY, mcy, float);   // MC Y
} // namespace pi0Gen

DECLARE_SOA_TABLE(Pi0Gens, "AOD", "PI0GENS",
                  pi0Gen::ProducedByGenerator,
                  pi0Gen::MCPt,
                  pi0Gen::MCY);

DECLARE_SOA_TABLE(Pi0GenCollRef, "AOD", "PI0GENCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraMCCollisionId);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMATABLES_H_
