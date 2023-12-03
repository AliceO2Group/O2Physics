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
#ifndef PWGLF_DATAMODEL_LFSTRANGENESSTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"

namespace o2::aod
{
//______________________________________________________
// Collision declarations for derived data analysis
// this is optional but will ensure full flexibility
// if required (for 2pc, etc)
DECLARE_SOA_TABLE(StraCollisions, "AOD", "STRACOLLISION", //! basic collision properties: position
                  collision::PosX, collision::PosY, collision::PosZ);
DECLARE_SOA_TABLE(StraCents, "AOD", "STRACENTS", //! centrality percentiles
                  cent::CentFT0M, cent::CentFT0A,
                  cent::CentFT0C, cent::CentFV0A);
DECLARE_SOA_TABLE(StraEPs, "AOD", "STRAEPS", //! centrality percentiles
                  qvec::QvecFT0ARe, qvec::QvecFT0AIm, qvec::SumAmplFT0A,
                  qvec::QvecFT0CRe, qvec::QvecFT0CIm, qvec::SumAmplFT0C,
                  qvec::QvecFT0MRe, qvec::QvecFT0MIm, qvec::SumAmplFT0M,
                  qvec::QvecFV0ARe, qvec::QvecFV0AIm, qvec::SumAmplFV0A);

using StraCollision = StraCollisions::iterator;
using StraCent = StraCents::iterator;
using StraEP = StraEPs::iterator;

namespace dautrack
{
//______________________________________________________
// Daughter track declarations for derived data analysis
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t);          //! detector map for reference
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t); //! ITS cluster sizes per layer
DECLARE_SOA_COLUMN(TPCClusters, tpcClusters, uint8_t);          //! N TPC clusters
DECLARE_SOA_COLUMN(TPCCrossedRows, tpcCrossedRows, uint8_t);    //! N TPC clusters

//______________________________________________________
// for extras: replicated here to ensure ease of manipulating the ITS information
// directly from the V0 extras table in simple ways for derived data as well
DECLARE_SOA_DYNAMIC_COLUMN(ITSClusterMap, itsClusterMap, //! ITS cluster map, one bit per layer, starting from the innermost
                           [](uint32_t itsClusterSizes) -> uint8_t {
                             uint8_t clmap = 0;
                             for (unsigned int layer = 0; layer < 7; layer++) {
                               if ((itsClusterSizes >> (layer * 4)) & 0xf) {
                                 clmap |= (1 << layer);
                               }
                             }
                             return clmap;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ITSNCls, itsNCls, //! Number of ITS clusters
                           [](uint32_t itsClusterSizes) -> uint8_t {
                             uint8_t itsNcls = 0;
                             for (int layer = 0; layer < 7; layer++) {
                               if ((itsClusterSizes >> (layer * 4)) & 0xf)
                                 itsNcls++;
                             }
                             return itsNcls;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(HasITS, hasITS, //! Flag to check if track has a ITS match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::ITS; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTPC, hasTPC, //! Flag to check if track has a TPC match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TPC; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTRD, hasTRD, //! Flag to check if track has a TRD match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TRD; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF, //! Flag to check if track has a TOF measurement
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TOF; });
} // namespace dautrack

DECLARE_SOA_TABLE(DauTrackExtras, "AOD", "DAUTRACKEXTRA", //! detector properties of decay daughters
                  dautrack::DetectorMap, dautrack::ITSClusterSizes,
                  dautrack::TPCClusters, dautrack::TPCCrossedRows,

                  // Dynamic columns for manipulating information
                  dautrack::ITSClusterMap<dautrack::ITSClusterSizes>,
                  dautrack::ITSNCls<dautrack::ITSClusterSizes>,
                  dautrack::HasITS<dautrack::DetectorMap>,
                  dautrack::HasTPC<dautrack::DetectorMap>,
                  dautrack::HasTRD<dautrack::DetectorMap>,
                  dautrack::HasTOF<dautrack::DetectorMap>);

using DauTrackExtra = DauTrackExtras::iterator;

namespace motherParticle
{
DECLARE_SOA_COLUMN(Px, px, float);                              //! px
DECLARE_SOA_COLUMN(Py, py, float);                              //! py
DECLARE_SOA_COLUMN(Pz, pz, float);                              //! pz
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                      //! pdg code
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! primary criterion
} // namespace motherParticle

DECLARE_SOA_TABLE(MotherMCParticles, "AOD", "MOTHERMCPART", //! mother MC information
                  motherParticle::Px, motherParticle::Py, motherParticle::Pz,
                  motherParticle::PDGCode, motherParticle::IsPhysicalPrimary);

using MotherMCParticle = MotherMCParticles::iterator;

namespace v0data
{
//______________________________________________________
// REGULAR COLUMNS FOR INDEXING
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!
// FOR DERIVED
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrackExtra, posTrackExtra, int, DauTrackExtras, "_PosExtra"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrackExtra, negTrackExtra, int, DauTrackExtras, "_NegExtra"); //!
DECLARE_SOA_INDEX_COLUMN(StraCollision, straCollision);                                        //!
DECLARE_SOA_INDEX_COLUMN(MotherMCParticle, motherMCParticle);                                  //!

//______________________________________________________
// REGULAR COLUMNS FOR V0CORES
// General V0 properties: position, momentum
DECLARE_SOA_COLUMN(PosX, posX, float);   //! positive track X at min
DECLARE_SOA_COLUMN(NegX, negX, float);   //! negative track X at min
DECLARE_SOA_COLUMN(PxPos, pxpos, float); //! positive track px at min
DECLARE_SOA_COLUMN(PyPos, pypos, float); //! positive track py at min
DECLARE_SOA_COLUMN(PzPos, pzpos, float); //! positive track pz at min
DECLARE_SOA_COLUMN(PxNeg, pxneg, float); //! negative track px at min
DECLARE_SOA_COLUMN(PyNeg, pyneg, float); //! negative track py at min
DECLARE_SOA_COLUMN(PzNeg, pzneg, float); //! negative track pz at min
DECLARE_SOA_COLUMN(X, x, float);         //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);         //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);         //! decay position Z

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float); //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);         //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);         //! DCA negative prong to PV
DECLARE_SOA_COLUMN(V0CosPA, v0cosPA, float);               //! V0 CosPA
DECLARE_SOA_COLUMN(DCAV0ToPV, dcav0topv, float);           //! DCA V0 to PV

// Saved from finding: covariance matrix of parent track (on request)
DECLARE_SOA_COLUMN(PositionCovMat, positionCovMat, float[6]); //! covariance matrix elements
DECLARE_SOA_COLUMN(MomentumCovMat, momentumCovMat, float[6]); //! covariance matrix elements

// Saved from KF particle fit for specic table
DECLARE_SOA_COLUMN(KFV0Chi2, kfV0Chi2, float); //!

//______________________________________________________
// REGULAR COLUMNS FOR V0MCCORES
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                      //! V0 PDG Code
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);          //! V0 mother PDG code (for feeddown)
DECLARE_SOA_COLUMN(PDGCodePositive, pdgCodePositive, int);      //! V0 positive prong PDG code
DECLARE_SOA_COLUMN(PDGCodeNegative, pdgCodeNegative, int);      //! V0 negative prong PDG code
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! is V0 physical primary
DECLARE_SOA_COLUMN(XMC, xMC, float);                            //! V0 decay position X (cm)
DECLARE_SOA_COLUMN(YMC, yMC, float);                            //! V0 decay position Y (cm)
DECLARE_SOA_COLUMN(ZMC, zMC, float);                            //! V0 decay position Z (cm)
DECLARE_SOA_COLUMN(PxPosMC, pxPosMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyPosMC, pyPosMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzPosMC, pzPosMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxNegMC, pxNegMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyNegMC, pyNegMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzNegMC, pzNegMC, float);                    //! V0 positive daughter pz (GeV/c)

//______________________________________________________
// EXPRESSION COLUMNS
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! V0 px
                              float, 1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! V0 py
                              float, 1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! V0 pz
                              float, 1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum in GeV/c
                              nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) *
                                      (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                    (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg)));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Total momentum in GeV/c
                              nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) *
                                      (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                    (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                    (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)));
DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float, //! Phi in the range [0, 2pi)
                              o2::constants::math::PI + natan2(-1.0f * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg), -1.0f * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg)));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity, conditionally defined to avoid FPEs
                              ifnode((nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                            (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                            (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) -
                                      (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) < static_cast<float>(1e-7),
                                     ifnode((1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) < 0.f, -100.f, 100.f),
                                     0.5f * nlog((nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                                        (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                                        (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) +
                                                  (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) /
                                                 (nsqrt((1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) * (1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg) +
                                                        (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) * (1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg) +
                                                        (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg) * (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)) -
                                                  (1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg)))));

//______________________________________________________
// DYNAMIC COLUMNS
// Account for rigidity in case of hypertriton
DECLARE_SOA_DYNAMIC_COLUMN(PtHypertriton, ptHypertriton, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(2.0f * pxpos + pxneg, 2.0f * pypos + pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PtAntiHypertriton, ptAntiHypertriton, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + 2.0f * pxneg, pypos + 2.0f * pyneg); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //! V0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distovertotmom, //! PV to V0decay distance over total momentum
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) {
                             float P = RecoDecay::sqrtSumOfSquares(Px, Py, Pz);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (P + 1E-10);
                           });

// Armenteros-Podolanski variables
DECLARE_SOA_DYNAMIC_COLUMN(Alpha, alpha, //! Armenteros Alpha
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0
                           });

DECLARE_SOA_DYNAMIC_COLUMN(QtArm, qtarm, //! Armenteros Qt
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
                             return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qtarm
                           });

// Psi pair angle: angle between the plane defined by the electron and positron momenta and the xy plane
DECLARE_SOA_DYNAMIC_COLUMN(PsiPair, psipair, //! psi pair angle
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
                             float ptot2 = RecoDecay::p2(pxpos, pypos, pzpos) * RecoDecay::p2(pxneg, pyneg, pzneg);
                             float argcos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
                             float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxpos, pypos), pzpos);
                             float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxneg, pyneg), pzneg);
                             float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
                             return std::asin(clipToPM1(argsin));
                           });

// calculate the fraction of the pos/neg momentum of the V0 momentum
DECLARE_SOA_DYNAMIC_COLUMN(PFracPos, pfracpos,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float ppos = RecoDecay::sqrtSumOfSquares(pxpos, pypos, pzpos);
                             float PV0 = RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             return (ppos / PV0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PFracNeg, pfracneg,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float pneg = RecoDecay::sqrtSumOfSquares(pxneg, pyneg, pzneg);
                             float PV0 = RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             return (pneg / PV0);
                           });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiLambda, mAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
DECLARE_SOA_DYNAMIC_COLUMN(MK0Short, mK0Short, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MGamma, mGamma, //! mass under gamma hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron}); });
// Account for rigidity in case of hypertriton
DECLARE_SOA_DYNAMIC_COLUMN(MHypertriton, mHypertriton, //! mass under hypertriton  hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{2.0f * pxpos, 2.0f * pypos, 2.0f * pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHypertriton, mAntiHypertriton, //! mass under antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{2.0f * pxneg, 2.0f * pyneg, 2.0f * pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassHelium3}); });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! mass under a certain hypothesis (0:K0, 1:L, 2:Lbar, 3:gamma, 4:hyp, 5:ahyp)
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg, int value) -> float {
                             if (value == 0)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                             if (value == 1)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
                             if (value == 2)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                             if (value == 3)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
                             if (value == 4)
                               return RecoDecay::m(std::array{std::array{2.0f * pxpos, 2.0f * pypos, 2.0f * pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPionCharged});
                             if (value == 5)
                               return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{2.0f * pxneg, 2.0f * pyneg, 2.0f * pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassHelium3});
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(YK0Short, yK0Short, //! V0 y with K0short hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassKaonNeutral); });
DECLARE_SOA_DYNAMIC_COLUMN(YLambda, yLambda, //! V0 y with lambda or antilambda hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassLambda); });
DECLARE_SOA_DYNAMIC_COLUMN(YHypertriton, yHypertriton, //! V0 y with hypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::y(std::array{2.0f * pxpos + pxneg, 2.0f * pypos + pyneg, 2.0f * pzpos + pzneg}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(YAntiHypertriton, yAntiHypertriton, //! V0 y with antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::y(std::array{pxpos + 2.0f * pxneg, pypos + 2.0f * pyneg, pzpos + 2.0f * pzneg}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! rapidity (0:K0, 1:L, 2:Lbar)
                           [](float Px, float Py, float Pz, int value) -> float {
                             if (value == 0)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassKaonNeutral);
                             if (value == 1 || value == 2)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassLambda);
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativept, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivept, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeeta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::eta(std::array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativephi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveeta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(std::array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivephi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });
} // namespace v0data

DECLARE_SOA_TABLE(V0Indices, "AOD", "V0INDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, v0data::PosTrackId, v0data::NegTrackId, v0data::CollisionId, v0data::V0Id);

DECLARE_SOA_TABLE(V0CollRefs, "AOD", "V0COLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, v0data::StraCollisionId);

DECLARE_SOA_TABLE(V0Extras, "AOD", "V0EXTRA", //! optional table to refer to custom track extras
                  o2::soa::Index<>, v0data::PosTrackExtraId, v0data::NegTrackExtraId);

DECLARE_SOA_TABLE(V0TrackXs, "AOD", "V0TRACKX", //! track X positions at minima when using AO2Ds
                  v0data::PosX, v0data::NegX);

DECLARE_SOA_TABLE_FULL(StoredV0Cores, "V0Cores", "AOD", "V0CORE", //! core information about decay, viable with AO2Ds or derived
                       v0data::X, v0data::Y, v0data::Z,
                       v0data::PxPos, v0data::PyPos, v0data::PzPos,
                       v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                       v0data::DCAV0Daughters, v0data::DCAPosToPV, v0data::DCANegToPV,
                       v0data::V0CosPA, v0data::DCAV0ToPV,

                       // Dynamic columns
                       v0data::PtHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::PtAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::V0Radius<v0data::X, v0data::Y>,
                       v0data::DistOverTotMom<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::Alpha<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::QtArm<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PsiPair<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PFracPos<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PFracNeg<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                       // Invariant masses
                       v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::M<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                       // Longitudinal
                       v0data::YK0Short<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::YLambda<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::YHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::YAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::Rapidity<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::NegativePt<v0data::PxNeg, v0data::PyNeg>,
                       v0data::PositivePt<v0data::PxPos, v0data::PyPos>,
                       v0data::NegativeEta<v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::NegativePhi<v0data::PxNeg, v0data::PyNeg>,
                       v0data::PositiveEta<v0data::PxPos, v0data::PyPos, v0data::PzPos>,
                       v0data::PositivePhi<v0data::PxPos, v0data::PyPos>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(V0Cores, StoredV0Cores, "V0COREEXT",                                                  //!
                                v0data::Px, v0data::Py, v0data::Pz, v0data::Pt, v0data::P, v0data::Phi, v0data::Eta); // the table name has here to be the one with EXT which is not nice and under study

DECLARE_SOA_TABLE_FULL(V0Covs, "V0Covs", "AOD", "V0COVS", //! V0 covariance matrices
                       v0data::PositionCovMat, v0data::MomentumCovMat);

DECLARE_SOA_TABLE(V0MCCores, "AOD", "V0MCCORE", //! MC properties of the V0 for posterior analysis
                  v0data::PDGCode, v0data::PDGCodeMother,
                  v0data::PDGCodePositive, v0data::PDGCodeNegative,
                  v0data::IsPhysicalPrimary, v0data::XMC, v0data::YMC, v0data::ZMC,
                  v0data::PxPosMC, v0data::PyPosMC, v0data::PzPosMC,
                  v0data::PxNegMC, v0data::PyNegMC, v0data::PzNegMC);

DECLARE_SOA_TABLE(V0MCMothers, "AOD", "V0MCMOTHER", //! optional table for MC mothers
                  o2::soa::Index<>, v0data::MotherMCParticleId);

using V0Index = V0Indices::iterator;
using V0Core = V0Cores::iterator;
using V0TrackX = V0TrackXs::iterator;
using V0Datas = soa::Join<V0Indices, V0TrackXs, V0Cores>;
using V0Data = V0Datas::iterator;
using V0MCDatas = soa::Join<V0MCCores, V0MCMothers>;
using V0MCData = V0MCDatas::iterator;

// definitions of indices for interlink tables
namespace v0data
{
DECLARE_SOA_INDEX_COLUMN(V0Data, v0Data); //! Index to V0Data entry
DECLARE_SOA_INDEX_COLUMN_FULL(V0MC, v0MC, int, V0MCCores, "_MC"); //!
}

DECLARE_SOA_TABLE(V0DataLink, "AOD", "V0DATALINK", //! Joinable table with V0s which links to V0Data which is not produced for all entries
                  v0data::V0DataId);
DECLARE_SOA_TABLE(V0MCRefs, "AOD", "V0MCREF", //! index table when using AO2Ds
                  o2::soa::Index<>, v0data::V0MCId);

using V0sLinked = soa::Join<V0s, V0DataLink>;
using V0Linked = V0sLinked::iterator;

// helper for building
namespace v0tag
{
// Global bool
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueGamma, isTrueGamma, bool);                     //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueK0Short, isTrueK0Short, bool);                 //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueLambda, isTrueLambda, bool);                   //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiLambda, isTrueAntiLambda, bool);           //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueHypertriton, isTrueHypertriton, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiHypertriton, isTrueAntiHypertriton, bool); //! PDG checked correctly in MC

// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsdEdxGamma, isdEdxGamma, bool);                     //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxK0Short, isdEdxK0Short, bool);                 //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxLambda, isdEdxLambda, bool);                   //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxAntiLambda, isdEdxAntiLambda, bool);           //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxHypertriton, isdEdxHypertriton, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxAntiHypertriton, isdEdxAntiHypertriton, bool); //! compatible with dE/dx hypotheses

// used in cascades (potentially useful in general, make available as tags)
DECLARE_SOA_COLUMN(IsFromCascade, isFromCascade, bool);               //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsFromTrackedCascade, isFromTrackedCascade, bool); //! compatible with dE/dx hypotheses
} // namespace v0tag
DECLARE_SOA_TABLE(V0Tags, "AOD", "V0TAGS",
                  v0tag::IsInteresting,
                  v0tag::IsTrueGamma,
                  v0tag::IsTrueK0Short,
                  v0tag::IsTrueLambda,
                  v0tag::IsTrueAntiLambda,
                  v0tag::IsTrueHypertriton,
                  v0tag::IsTrueAntiHypertriton,
                  v0tag::IsdEdxGamma,
                  v0tag::IsdEdxK0Short,
                  v0tag::IsdEdxLambda,
                  v0tag::IsdEdxAntiLambda,
                  v0tag::IsdEdxHypertriton,
                  v0tag::IsdEdxAntiHypertriton,
                  v0tag::IsFromCascade,
                  v0tag::IsFromTrackedCascade);

namespace kfcascdata
{
// declare in different namespace to 'overload' operator
DECLARE_SOA_COLUMN(MLambda, mLambda, float); //!
} // namespace kfcascdata

namespace cascdata
{
//______________________________________________________
// REGULAR COLUMNS FOR CASCINDICES
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                           //!
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);                                 //!
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, Tracks, "");         //!
DECLARE_SOA_INDEX_COLUMN_FULL(StrangeTrack, strangeTrack, int, Tracks, ""); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                             //!
// FOR DERIVED
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrackExtra, posTrackExtra, int, DauTrackExtras, "_PosExtra");   //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrackExtra, negTrackExtra, int, DauTrackExtras, "_NegExtra");   //!
DECLARE_SOA_INDEX_COLUMN_FULL(BachTrackExtra, bachTrackExtra, int, DauTrackExtras, "_BachExtra"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(StrangeTrackExtra, strangeTrackExtra, int, DauTrackExtras, "_StrangeExtra"); //!
DECLARE_SOA_INDEX_COLUMN(StraCollision, straCollision);                                          //!
DECLARE_SOA_INDEX_COLUMN(MotherMCParticle, motherMCParticle);                                    //!

//______________________________________________________
// REGULAR COLUMNS FOR CASCCORES
// General cascade properties: position, momentum
DECLARE_SOA_COLUMN(Sign, sign, int);         //!
DECLARE_SOA_COLUMN(MXi, mXi, float);         //!
DECLARE_SOA_COLUMN(MOmega, mOmega, float);   //!
DECLARE_SOA_COLUMN(PxPos, pxpos, float);     //!
DECLARE_SOA_COLUMN(PyPos, pypos, float);     //!
DECLARE_SOA_COLUMN(PzPos, pzpos, float);     //!
DECLARE_SOA_COLUMN(PxNeg, pxneg, float);     //!
DECLARE_SOA_COLUMN(PyNeg, pyneg, float);     //!
DECLARE_SOA_COLUMN(PzNeg, pzneg, float);     //!
DECLARE_SOA_COLUMN(PxBach, pxbach, float);   //!
DECLARE_SOA_COLUMN(PyBach, pybach, float);   //!
DECLARE_SOA_COLUMN(PzBach, pzbach, float);   //!
DECLARE_SOA_COLUMN(Px, px, float);           //! cascade momentum X
DECLARE_SOA_COLUMN(Py, py, float);           //! cascade momentum Y
DECLARE_SOA_COLUMN(Pz, pz, float);           //! cascade momentum Z
DECLARE_SOA_COLUMN(X, x, float);             //!
DECLARE_SOA_COLUMN(Y, y, float);             //!
DECLARE_SOA_COLUMN(Z, z, float);             //!
DECLARE_SOA_COLUMN(Xlambda, xlambda, float); //!
DECLARE_SOA_COLUMN(Ylambda, ylambda, float); //!
DECLARE_SOA_COLUMN(Zlambda, zlambda, float); //!

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float);     //!
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float); //!
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);             //!
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);             //!
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);           //!
DECLARE_SOA_COLUMN(DCAXYCascToPV, dcaXYCascToPV, float);       //!
DECLARE_SOA_COLUMN(DCAZCascToPV, dcaZCascToPV, float);         //!

//______________________________________________________
// REGULAR COLUMNS FOR CASCCOVS
// Saved from finding: covariance matrix of parent track (on request)
DECLARE_SOA_COLUMN(PositionCovMat, positionCovMat, float[6]); //! covariance matrix elements
DECLARE_SOA_COLUMN(MomentumCovMat, momentumCovMat, float[6]); //! covariance matrix elements
DECLARE_SOA_COLUMN(KFTrackMat, kfTrackCovMat, float[15]);     //! covariance matrix elements for KF method

//______________________________________________________
// REGULAR COLUMNS FOR CASCBBS
// General cascade properties: position, momentum
// Selection to avoid spurious invariant mass correlation
// bachelor-baryon cosine of pointing angle / DCA to PV
DECLARE_SOA_COLUMN(BachBaryonCosPA, bachBaryonCosPA, float);         //! avoid bach-baryon correlated inv mass structure in analysis
DECLARE_SOA_COLUMN(BachBaryonDCAxyToPV, bachBaryonDCAxyToPV, float); //! avoid bach-baryon correlated inv mass structure in analysis

//______________________________________________________
// REGULAR COLUMNS FOR KFCASCCORES
// General cascade properties: position, momentum
// Saved from KF particle fit for specic table
// note: separate chi2 is a consequence of fit -> conversion -> propagation -> fit logic
//       which, in turn, is necessary to do material corrections at the moment
//       this could be improved in the future!
DECLARE_SOA_COLUMN(KFV0Chi2, kfV0Chi2, float);           //!
DECLARE_SOA_COLUMN(KFCascadeChi2, kfCascadeChi2, float); //!

//______________________________________________________
// REGULAR COLUMNS FOR TRACASCCORES
// Saved from strangeness tracking
DECLARE_SOA_COLUMN(MatchingChi2, matchingChi2, float); //!
DECLARE_SOA_COLUMN(TopologyChi2, topologyChi2, float); //!
DECLARE_SOA_COLUMN(ItsClsSize, itsCluSize, float);     //!

//______________________________________________________
// REGULAR COLUMNS FOR CASCMCCORES
DECLARE_SOA_COLUMN(PDGCode, pdgCode, int);                      //! cascade PDG Code
DECLARE_SOA_COLUMN(PDGCodeMother, pdgCodeMother, int);          //! cascade mother PDG code (for feeddown)
DECLARE_SOA_COLUMN(PDGCodeV0, pdgCodeV0, int);                  //! cascade PDG Code
DECLARE_SOA_COLUMN(PDGCodePositive, pdgCodePositive, int);      //! V0 positive prong PDG code
DECLARE_SOA_COLUMN(PDGCodeNegative, pdgCodeNegative, int);      //! V0 negative prong PDG code
DECLARE_SOA_COLUMN(PDGCodeBachelor, pdgCodeBachelor, int);      //! cascade bachelor prong PDG code
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool); //! is cascade physical primary
DECLARE_SOA_COLUMN(XMC, xMC, float);                            //! cascade decay position X (cm)
DECLARE_SOA_COLUMN(YMC, yMC, float);                            //! cascade decay position Y (cm)
DECLARE_SOA_COLUMN(ZMC, zMC, float);                            //! cascade decay position Z (cm)
DECLARE_SOA_COLUMN(XlambdaMC, xlambdaMC, float);                //! V0 decay position X (cm)
DECLARE_SOA_COLUMN(YlambdaMC, ylambdaMC, float);                //! V0 decay position Y (cm)
DECLARE_SOA_COLUMN(ZlambdaMC, zlambdaMC, float);                //! V0 decay position Z (cm)
DECLARE_SOA_COLUMN(PxPosMC, pxPosMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyPosMC, pyPosMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzPosMC, pzPosMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxNegMC, pxNegMC, float);                    //! V0 positive daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyNegMC, pyNegMC, float);                    //! V0 positive daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzNegMC, pzNegMC, float);                    //! V0 positive daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxBachMC, pxBachMC, float);                  //! cascade bachelor daughter px (GeV/c)
DECLARE_SOA_COLUMN(PyBachMC, pyBachMC, float);                  //! cascade bachelor daughter py (GeV/c)
DECLARE_SOA_COLUMN(PzBachMC, pzBachMC, float);                  //! cascade bachelor daughter pz (GeV/c)
DECLARE_SOA_COLUMN(PxMC, pxMC, float);                          //! cascade px (GeV/c)
DECLARE_SOA_COLUMN(PyMC, pyMC, float);                          //! cascade py (GeV/c)
DECLARE_SOA_COLUMN(PzMC, pzMC, float);                          //! cascade pz (GeV/c)

//______________________________________________________
// DERIVED
// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //!
                           [](float xlambda, float ylambda) -> float { return RecoDecay::sqrtSumOfSquares(xlambda, ylambda); });
DECLARE_SOA_DYNAMIC_COLUMN(CascRadius, cascradius, //!
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// CosPAs
DECLARE_SOA_DYNAMIC_COLUMN(V0CosPA, v0cosPA, //!
                           [](float Xlambda, float Ylambda, float Zlambda, float PxLambda, float PyLambda, float PzLambda, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{Xlambda, Ylambda, Zlambda}, std::array{PxLambda, PyLambda, PzLambda}); });
DECLARE_SOA_DYNAMIC_COLUMN(CascCosPA, casccosPA, //!
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{X, Y, Z}, std::array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAV0ToPV, dcav0topv, //!
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //!
                           [](int charge, float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, charge < 0 ? std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged} : std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! mass under a certain hypothesis (0:K0, 1:L, 2:Lbar, 3:gamma, 4:hyp, 5:ahyp)
                           [](float mXi, float mOmega, int value) -> float {
                             if (value == 0 || value == 1)
                               return mXi;
                             if (value == 2 || value == 3)
                               return mOmega;
                             return 0.0f;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(YXi, yXi, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassXiMinus); });
DECLARE_SOA_DYNAMIC_COLUMN(YOmega, yOmega, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassOmegaMinus); });
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! rapidity (0, 1: Xi; 2, 3: Omega)
                           [](float Px, float Py, float Pz, int value) -> float {
                             if (value == 0 || value == 1)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassXiMinus);
                             if (value == 2 || value == 3)
                               return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassOmegaMinus);
                             return 0.0f;
                           });
} // namespace cascdata

//______________________________________________________
// EXPRESSION COLUMNS FOR TRACASCCORES
namespace cascdataext
{
DECLARE_SOA_EXPRESSION_COLUMN(PxLambda, pxlambda, //!
                              float, 1.f * aod::cascdata::pxpos + 1.f * aod::cascdata::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(PyLambda, pylambda, //!
                              float, 1.f * aod::cascdata::pypos + 1.f * aod::cascdata::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(PzLambda, pzlambda, //!
                              float, 1.f * aod::cascdata::pzpos + 1.f * aod::cascdata::pzneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum in GeV/c
                              nsqrt(aod::cascdata::px* aod::cascdata::px +
                                    aod::cascdata::py * aod::cascdata::py));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Total momentum in GeV/c
                              nsqrt(aod::cascdata::px* aod::cascdata::px +
                                    aod::cascdata::py * aod::cascdata::py +
                                    aod::cascdata::pz * aod::cascdata::pz));
DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float, //! Phi in the range [0, 2pi)
                              o2::constants::math::PI + natan2(-1.0f * aod::cascdata::py, -1.0f * aod::cascdata::px));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity, conditionally defined to avoid FPEs
                              ifnode((nsqrt(aod::cascdata::px * aod::cascdata::px +
                                            aod::cascdata::py * aod::cascdata::py +
                                            aod::cascdata::pz * aod::cascdata::pz) -
                                      aod::cascdata::pz) < static_cast<float>(1e-7),
                                     ifnode(aod::cascdata::pz < 0.f, -100.f, 100.f),
                                     0.5f * nlog((nsqrt(aod::cascdata::px * aod::cascdata::px +
                                                        aod::cascdata::py * aod::cascdata::py +
                                                        aod::cascdata::pz * aod::cascdata::pz) +
                                                  aod::cascdata::pz) /
                                                 (nsqrt(aod::cascdata::px * aod::cascdata::px +
                                                        aod::cascdata::py * aod::cascdata::py +
                                                        aod::cascdata::pz * aod::cascdata::pz) -
                                                  aod::cascdata::pz))));
} // namespace cascdataext

//______________________________________________________
// Cascade data model:
// --- standard, KF version and tracked version
// includes three variants to be able to account for
// cases in which multiple analyses are executed
// at the same time. Tables are kept as similar
// as possible for ease of use and comparison

DECLARE_SOA_TABLE(CascIndices, "AOD", "CascINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, cascdata::V0Id, cascdata::CascadeId, cascdata::BachelorId, cascdata::CollisionId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(KFCascIndices, "AOD", "KFCascINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, cascdata::V0Id, cascdata::CascadeId, cascdata::BachelorId, cascdata::CollisionId, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(TraCascIndices, "AOD", "TraCascINDEX", //! index table when using AO2Ds
                  o2::soa::Index<>, cascdata::V0Id, cascdata::CascadeId, cascdata::BachelorId, cascdata::StrangeTrackId, cascdata::CollisionId);

DECLARE_SOA_TABLE(CascCollRefs, "AOD", "CASCCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, cascdata::StraCollisionId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(KFCascCollRefs, "AOD", "KFCASCCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, cascdata::StraCollisionId, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(TraCascCollRefs, "AOD", "TRACASCCOLLREF", //! optional table to refer back to a collision
                  o2::soa::Index<>, cascdata::StraCollisionId, o2::soa::Marker<3>);

DECLARE_SOA_TABLE(CascExtras, "AOD", "CASCEXTRA", //! optional table to refer to custom track extras
                  o2::soa::Index<>, cascdata::PosTrackExtraId, cascdata::NegTrackExtraId,
                  cascdata::BachTrackExtraId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(StraTrackExtras, "AOD", "STRATRACKEXTRAS", //! optional table to refer to custom track extras
                  o2::soa::Index<>, cascdata::StrangeTrackExtraId);

DECLARE_SOA_TABLE(StoredCascCores, "AOD", "CASCCORE", //! core information about decay, viable with AO2Ds or derived
                  cascdata::Sign, cascdata::MXi, cascdata::MOmega,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::Px, cascdata::Py, cascdata::Pz,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCAXYCascToPV, cascdata::DCAZCascToPV,

                  // Dynamic columns
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::MLambda<cascdata::Sign, cascdata::PxPos, cascdata::PyPos, cascdata::PzPos, cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,
                  cascdata::M<cascdata::MXi, cascdata::MOmega>,

                  // Longitudinal
                  cascdata::YXi<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::YOmega<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::Rapidity<cascdata::Px, cascdata::Py, cascdata::Pz>);

DECLARE_SOA_TABLE(StoredKFCascCores, "AOD", "KFCASCCORE", //!
                  cascdata::Sign, cascdata::MXi, cascdata::MOmega,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::Px, cascdata::Py, cascdata::Pz,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCAXYCascToPV, cascdata::DCAZCascToPV,

                  // KF particle fit specific
                  kfcascdata::MLambda, cascdata::KFV0Chi2, cascdata::KFCascadeChi2,

                  // Dynamic columns
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::M<cascdata::MXi, cascdata::MOmega>,

                  // Longitudinal
                  cascdata::YXi<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::YOmega<cascdata::Px, cascdata::Py, cascdata::Pz>);

DECLARE_SOA_TABLE(StoredTraCascCores, "AOD", "TRACASCCORE", //!
                  cascdata::Sign, cascdata::MXi, cascdata::MOmega,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::Px, cascdata::Py, cascdata::Pz,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCAXYCascToPV, cascdata::DCAZCascToPV,

                  // Strangeness tracking specific
                  cascdata::MatchingChi2, cascdata::TopologyChi2, cascdata::ItsClsSize,

                  // Dynamic columns
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::MLambda<cascdata::Sign, cascdata::PxPos, cascdata::PyPos, cascdata::PzPos, cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,

                  // Longitudinal
                  cascdata::YXi<cascdata::Px, cascdata::Py, cascdata::Pz>,
                  cascdata::YOmega<cascdata::Px, cascdata::Py, cascdata::Pz>);

DECLARE_SOA_TABLE(CascMCCores, "AOD", "CASCMCCORE", //! bachelor-baryon correlation variables
                  cascdata::PDGCode, cascdata::PDGCodeMother, cascdata::PDGCodeV0, cascdata::IsPhysicalPrimary,
                  cascdata::PDGCodePositive, cascdata::PDGCodeNegative, cascdata::PDGCodeBachelor,
                  cascdata::XMC, cascdata::YMC, cascdata::ZMC,
                  cascdata::XlambdaMC, cascdata::YlambdaMC, cascdata::ZlambdaMC,
                  cascdata::PxPosMC, cascdata::PyPosMC, cascdata::PzPosMC,
                  cascdata::PxNegMC, cascdata::PyNegMC, cascdata::PzNegMC,
                  cascdata::PxBachMC, cascdata::PyBachMC, cascdata::PzBachMC,
                  cascdata::PxMC, cascdata::PyMC, cascdata::PzMC);

DECLARE_SOA_TABLE(CascMCMothers, "AOD", "CASCMCMOTHER", //! optional table for MC mothers
                  o2::soa::Index<>, cascdata::MotherMCParticleId);

DECLARE_SOA_TABLE(CascBBs, "AOD", "CASCBB", //! bachelor-baryon correlation variables
                  cascdata::BachBaryonCosPA, cascdata::BachBaryonDCAxyToPV)

DECLARE_SOA_TABLE_FULL(CascCovs, "CascCovs", "AOD", "CASCCOVS", //!
                       cascdata::PositionCovMat, cascdata::MomentumCovMat);

DECLARE_SOA_TABLE_FULL(KFCascCovs, "KFCascCovs", "AOD", "KFCASCCOVS", //!
                       cascdata::KFTrackMat);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(CascCores, StoredCascCores, "CascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda, cascdataext::Pt, cascdataext::P, cascdataext::Eta, cascdataext::Phi);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(KFCascCores, StoredKFCascCores, "KFCascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda,
                                cascdataext::Pt, cascdataext::P, cascdataext::Eta, cascdataext::Phi);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(TraCascCores, StoredTraCascCores, "TraCascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda,
                                cascdataext::Pt, cascdataext::P, cascdataext::Eta, cascdataext::Phi);

namespace cascdata
{
// For cross-linking all cascade kinds
DECLARE_SOA_INDEX_COLUMN_FULL(TrackedCascade, trackedCascade, int, TraCascCores, "_Refs"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(KFCascade, kfCascade, int, KFCascCores, "_Refs");            //!
DECLARE_SOA_INDEX_COLUMN_FULL(StandardCascade, standardCascade, int, CascCores, "_Refs");  //!
} // namespace cascdata

// interlink different cascade types
DECLARE_SOA_TABLE(CascToTraRefs, "AOD", "CASCTOTRAREFS", //! standard -> tracked
                  o2::soa::Index<>, cascdata::TrackedCascadeId);
DECLARE_SOA_TABLE(CascToKFRefs, "AOD", "CASCTOKFREFS", //! standard -> KF
                  o2::soa::Index<>, cascdata::KFCascadeId);
DECLARE_SOA_TABLE(TraToCascRefs, "AOD", "TRATOCASCREFS", //! standard -> KF
                  o2::soa::Index<>, cascdata::StandardCascadeId, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(KFToCascRefs, "AOD", "KFTOCASCREFS", //! standard -> KF
                  o2::soa::Index<>, cascdata::StandardCascadeId, o2::soa::Marker<2>);

using CascIndex = CascIndices::iterator;
using CascCore = CascCores::iterator;
using KFCascIndex = KFCascIndices::iterator;
using KFCascCore = KFCascCores::iterator;
using TraCascIndex = TraCascIndices::iterator;
using TraCascCore = TraCascCores::iterator;

using CascDatas = soa::Join<CascIndices, CascBBs, CascCores>;
using KFCascDatas = soa::Join<KFCascIndices, KFCascCores>;
using TraCascDatas = soa::Join<TraCascIndices, TraCascCores>;

using CascData = CascDatas::iterator;
using KFCascData = KFCascDatas::iterator;
using TraCascData = TraCascDatas::iterator;

using CascMCCore = CascMCCores::iterator;
using CascMCMother = CascMCMothers::iterator;
using CascMCDatas = soa::Join<CascMCCores, CascMCMothers>;
using CascMCData = CascMCDatas::iterator;

// For compatibility with previous table declarations
using CascDataFull = CascDatas;
using CascDataExt = CascDatas;

namespace cascdata
{
DECLARE_SOA_INDEX_COLUMN(CascData, cascData); //! Index to CascData entry
DECLARE_SOA_INDEX_COLUMN(KFCascData, kfCascData); //! Index to CascData entry
DECLARE_SOA_INDEX_COLUMN(TraCascData, traCascData); //! Index to CascData entry
}

DECLARE_SOA_TABLE(CascDataLink, "AOD", "CASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::CascDataId);
DECLARE_SOA_TABLE(KFCascDataLink, "AOD", "KFCASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::KFCascDataId);
DECLARE_SOA_TABLE(TraCascDataLink, "AOD", "TRACASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::TraCascDataId);

using CascadesLinked = soa::Join<Cascades, CascDataLink>;
using CascadeLinked = CascadesLinked::iterator;
using KFCascadesLinked = soa::Join<Cascades, KFCascDataLink>;
using KFCascadeLinked = KFCascadesLinked::iterator;

namespace casctag
{
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueXiMinus, isTrueXiMinus, bool);       //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueXiPlus, isTrueXiPlus, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueOmegaMinus, isTrueOmegaMinus, bool); //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueOmegaPlus, isTrueOmegaPlus, bool);   //! PDG checked correctly in MC

// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsdEdxXiMinus, isdEdxXiMinus, bool);       //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxXiPlus, isdEdxXiPlus, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxOmegaMinus, isdEdxOmegaMinus, bool); //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsdEdxOmegaPlus, isdEdxOmegaPlus, bool);   //! compatible with dE/dx hypotheses
} // namespace casctag
DECLARE_SOA_TABLE(CascTags, "AOD", "CASCTAGS",
                  casctag::IsInteresting,
                  casctag::IsTrueXiMinus,
                  casctag::IsTrueXiPlus,
                  casctag::IsTrueOmegaMinus,
                  casctag::IsTrueOmegaPlus,
                  casctag::IsdEdxXiMinus,
                  casctag::IsdEdxXiPlus,
                  casctag::IsdEdxOmegaMinus,
                  casctag::IsdEdxOmegaPlus);

// Definition of labels for V0s
namespace mcv0label
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mcv0label

DECLARE_SOA_TABLE(McV0Labels, "AOD", "MCV0LABEL", //! Table joinable with V0Data containing the MC labels
                  mcv0label::McParticleId);
using McV0Label = McV0Labels::iterator;

// Definition of labels for V0s // Full table, joinable with V0 (CAUTION: NOT WITH V0DATA)
namespace mcfullv0label
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mcfullv0label

DECLARE_SOA_TABLE(McFullV0Labels, "AOD", "MCFULLV0LABEL", //! Table joinable with V0
                  mcfullv0label::McParticleId);
using McFullV0Label = McFullV0Labels::iterator;

// Definition of labels for cascades
namespace mccasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for Cascade
DECLARE_SOA_COLUMN(IsBachBaryonCandidate, isBachBaryonCandidate, bool); //! will this be built or not?
} // namespace mccasclabel

DECLARE_SOA_TABLE(McCascLabels, "AOD", "MCCASCLABEL", //! Table joinable with CascData containing the MC labels
                  mccasclabel::McParticleId);
DECLARE_SOA_TABLE(McCascBBTags, "AOD", "MCCASCBBTAG", //! Table joinable with CascData containing yes / no for BB correlation
                  mccasclabel::IsBachBaryonCandidate);
using McCascLabel = McCascLabels::iterator;
using McCascBBTag = McCascBBTags::iterator;

// Definition of labels for kf cascades
namespace mckfcasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mckfcasclabel

DECLARE_SOA_TABLE(McKFCascLabels, "AOD", "MCKFCASCLABEL", //! Table joinable to cascdata containing the MC labels
                  mckfcasclabel::McParticleId);
using McKFCascLabel = McKFCascLabels::iterator;

// Definition of labels for tracked cascades
namespace mctracasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mctracasclabel

DECLARE_SOA_TABLE(McTraCascLabels, "AOD", "MCTRACASCLABEL", //! Table joinable to cascdata containing the MC labels
                  mctracasclabel::McParticleId);
using McTraCascLabel = McTraCascLabels::iterator;

DECLARE_SOA_TABLE(TrackedCascadeColls, "AOD", "TRACASCCOLL", //! Table joinable with TrackedCascades containing collision ids
                  track::CollisionId, o2::soa::Marker<1>);
using TrackedCascadeColl = TrackedCascadeColls::iterator;
using AssignedTrackedCascades = soa::Join<aod::TrackedCascades, aod::TrackedCascadeColls>;
using AssignedTrackedCascade = AssignedTrackedCascades::iterator;

DECLARE_SOA_TABLE(TrackedV0Colls, "AOD", "TRAV0COLL", //! Table joinable with TrackedV0s containing collision ids
                  track::CollisionId, o2::soa::Marker<2>);
using TrackedV0Coll = TrackedV0Colls::iterator;
using AssignedTrackedV0s = soa::Join<aod::TrackedV0s, aod::TrackedV0Colls>;
using AssignedTrackedV0 = AssignedTrackedV0s::iterator;

DECLARE_SOA_TABLE(Tracked3BodyColls, "AOD", "TRA3BODYCOLL", //! Table joinable with Tracked3Bodys containing collision ids
                  track::CollisionId, o2::soa::Marker<3>);
using Tracked3BodyColl = Tracked3BodyColls::iterator;
using AssignedTracked3Bodys = soa::Join<aod::Tracked3Bodys, aod::Tracked3BodyColls>;
using AssignedTracked3Body = AssignedTracked3Bodys::iterator;
} // namespace o2::aod

//______________________________________________________
// Equivalency declarations
namespace o2::soa
{
DECLARE_EQUIVALENT_FOR_INDEX(aod::V0Indices, aod::V0Cores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::V0TrackXs, aod::V0Cores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::V0TrackXs, aod::V0Indices);
DECLARE_EQUIVALENT_FOR_INDEX(aod::CascIndices, aod::CascCores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::CascIndices, aod::CascBBs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::CascCores, aod::CascBBs);
DECLARE_EQUIVALENT_FOR_INDEX(aod::KFCascIndices, aod::KFCascCores);
DECLARE_EQUIVALENT_FOR_INDEX(aod::TraCascIndices, aod::TraCascCores);
} // namespace o2::soa

#endif // PWGLF_DATAMODEL_LFSTRANGENESSTABLES_H_
