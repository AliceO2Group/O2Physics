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

///
/// \file LFLambda1405Tables.h
/// \brief Slim tables for Lambda(1405) candidates
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>
///

#ifndef PWGLF_DATAMODEL_LFLAMBDA1405TABLE_H_
#define PWGLF_DATAMODEL_LFLAMBDA1405TABLE_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{

namespace lambda1405
{

DECLARE_SOA_COLUMN(Px, px, float);                           //! Px of the candidate
DECLARE_SOA_COLUMN(Py, py, float);                           //! Py of the candidate
DECLARE_SOA_COLUMN(Pz, pz, float);                           //! Pz of the candidate
DECLARE_SOA_COLUMN(Pt, pt, float);                           //! Pt of the candidate
DECLARE_SOA_COLUMN(Mass, mass, float);                       //! Invariant mass of the candidate
DECLARE_SOA_COLUMN(MassXi1530, massXi1530, float);           //! Invariant mass of the Xi(1530) candidate
DECLARE_SOA_COLUMN(SigmaMinusMass, sigmaMinusMass, float);   //! Invariant mass of the Sigma- candidate
DECLARE_SOA_COLUMN(SigmaPlusMass, sigmaPlusMass, float);     //! Invariant mass of the Sigma+ candidate
DECLARE_SOA_COLUMN(XiMinusMass, xiMinusMass, float);         //! Invariant mass of the Xi- candidate
DECLARE_SOA_COLUMN(PtSigma, ptSigma, float);                 //! Signed pT of the Sigma daughter
DECLARE_SOA_COLUMN(AlphaAPSigma, alphaAPSigma, float);       //! Alpha of the Sigma
DECLARE_SOA_COLUMN(QtAPSigma, qtAPSigma, float);             //! qT of the Sigma
DECLARE_SOA_COLUMN(RadiusSigma, radiusSigma, float);         //! Radius of the Sigma decay vertex
DECLARE_SOA_COLUMN(PtKink, ptKink, float);                   //! pT of the kink daughter
DECLARE_SOA_COLUMN(NSigmaTPCPiKink, nSigmaTPCPiKink, float); //! Number of sigmas for the pion candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTOFPiKink, nSigmaTOFPiKink, float); //! Number of sigmas for the pion candidate from Sigma kink in TOF
DECLARE_SOA_COLUMN(NSigmaTPCPrKink, nSigmaTPCPrKink, float); //! Number of sigmas for the proton candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTOFPrKink, nSigmaTOFPrKink, float); //! Number of sigmas for the proton candidate from Sigma kink in TOF
DECLARE_SOA_COLUMN(DCAKinkDauToPV, dcaKinkDauToPV, float);   //! DCA of the kink daughter to the primary vertex
DECLARE_SOA_COLUMN(DCASigmaToPV, dcaSigmaToPV, float);       //! DCA of the sigma to the primary vertex
DECLARE_SOA_COLUMN(NSigmaTPCPiDau, nSigmaTPCPiDau, float);   //! Number of sigmas for the lambda1405 pion daughter in TPC
DECLARE_SOA_COLUMN(NSigmaTOFPiDau, nSigmaTOFPiDau, float);   //! Number of sigmas for the lambda1405 pion daughter in TOF

// Event properties
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! Centrality of the candidate
DECLARE_SOA_COLUMN(Occupancy, occupancy, float);   //! Occupancy of the candidate
DECLARE_SOA_COLUMN(PvContrib, pvContrib, float);   //! Number of primary vertex contributors

// Flow columns
DECLARE_SOA_COLUMN(ScalarProd, scalarProd, float); //! Scalar product of the candidate

// MC Columns
DECLARE_SOA_COLUMN(PtMC, ptMC, float);                   //! pT of the candidate in MC
DECLARE_SOA_COLUMN(MassMC, massMC, float);               //! Invariant mass of the candidate in MC
DECLARE_SOA_COLUMN(SigmaPdgCode, sigmaPdgCode, int);     //! PDG code of the Sigma daughter
DECLARE_SOA_COLUMN(KinkDauPdgCode, kinkDauPdgCode, int); //! PDG code of the kink daughter

// Sigma efficiency MC columns
DECLARE_SOA_COLUMN(PxSigma, pxSigma, float);                       //! Px of the sigma candidate
DECLARE_SOA_COLUMN(PySigma, pySigma, float);                       //! Py of the sigma candidate
DECLARE_SOA_COLUMN(PzSigma, pzSigma, float);                       //! Pz of the sigma candidate
DECLARE_SOA_COLUMN(MassSigma, massSigma, float);                   //! Mass of the sigma candidate
DECLARE_SOA_COLUMN(DeltaPxSigma, deltaPxSigma, float);             //! Gen-reco diff of sigma p_x
DECLARE_SOA_COLUMN(DeltaPySigma, deltaPySigma, float);             //! Gen-reco diff of sigma p_y
DECLARE_SOA_COLUMN(DeltaPzSigma, deltaPzSigma, float);             //! Gen-reco diff of sigma p_z
DECLARE_SOA_COLUMN(DeltaPtSigma, deltaPtSigma, float);             //! Gen-reco diff of sigma p_z
DECLARE_SOA_COLUMN(DeltaRadiusSigma, deltaRadiusSigma, float);     //! Gen-reco diff of sigma radius
DECLARE_SOA_COLUMN(DeltaMassSigma, deltaMassSigma, float);         //! Gen-reco diff of sigma mass
DECLARE_SOA_COLUMN(GenPhiSigma, genPhiSigma, float);               //! Gen-reco diff of sigma phi
DECLARE_SOA_COLUMN(GenEtaSigma, genEtaSigma, float);               //! Gen-reco diff of sigma eta
DECLARE_SOA_COLUMN(DeltaPxSigmaRecalc, deltaPxSigmaRecalc, float); //! reco recalc-original diff of sigma p_x
DECLARE_SOA_COLUMN(DeltaPySigmaRecalc, deltaPySigmaRecalc, float); //! reco recalc-original diff of sigma p_y
DECLARE_SOA_COLUMN(DeltaPzSigmaRecalc, deltaPzSigmaRecalc, float); //! reco recalc-original diff of sigma p_z
DECLARE_SOA_COLUMN(AlphaAPSigmaRecalc, alphaAPSigmaRecalc, float); //! Alpha of the Sigma
DECLARE_SOA_COLUMN(QtAPSigmaRecalc, qtAPSigmaRecalc, float);       //! qT of the Sigma
DECLARE_SOA_COLUMN(PxKinkDaug, pxKinkDaug, float);                 //! Px of the sigma candidate
DECLARE_SOA_COLUMN(PyKinkDaug, pyKinkDaug, float);                 //! Py of the sigma candidate
DECLARE_SOA_COLUMN(PzKinkDaug, pzKinkDaug, float);                 //! Pz of the sigma candidate
DECLARE_SOA_COLUMN(PtKinkDaug, ptKinkDaug, float);                 //! Pt of the sigma candidate
DECLARE_SOA_COLUMN(DeltaPxKinkDaug, deltaPxKinkDaug, float);       //! Gen-reco diff of kink daughter p_x
DECLARE_SOA_COLUMN(DeltaPyKinkDaug, deltaPyKinkDaug, float);       //! Gen-reco diff of kink daughter p_y
DECLARE_SOA_COLUMN(DeltaPzKinkDaug, deltaPzKinkDaug, float);       //! Gen-reco diff of kink daughter p_z
DECLARE_SOA_COLUMN(DeltaPtKinkDaug, deltaPtKinkDaug, float);       //! Gen-reco diff of kink daughter p_t
DECLARE_SOA_COLUMN(GenPhiKinkDaug, genPhiKinkDaug, float);         //! Gen-reco diff of kink daughter phi
DECLARE_SOA_COLUMN(GenEtaKinkDaug, genEtaKinkDaug, float);         //! Gen-reco diff of kink daughter eta
DECLARE_SOA_COLUMN(XKinkVtx, xKinkVtx, float);                     //! X of kink vertex
DECLARE_SOA_COLUMN(YKinkVtx, yKinkVtx, float);                     //! Y of kink vertex
DECLARE_SOA_COLUMN(ZKinkVtx, zKinkVtx, float);                     //! Z of kink vertex
DECLARE_SOA_COLUMN(DeltaXKinkVtx, deltaXKinkVtx, float);           //! Gen-reco diff of X of kink vertex
DECLARE_SOA_COLUMN(DeltaYKinkVtx, deltaYKinkVtx, float);           //! Gen-reco diff of Y of kink vertex
DECLARE_SOA_COLUMN(DeltaZKinkVtx, deltaZKinkVtx, float);           //! Gen-reco diff of Z of kink vertex

} // namespace lambda1405

DECLARE_SOA_TABLE(Lambda1405Cands, "AOD", "LAMBDA1405",
                  o2::soa::Index<>,
                  lambda1405::Px, lambda1405::Py, lambda1405::Pz,
                  lambda1405::Mass, lambda1405::MassXi1530,
                  lambda1405::SigmaMinusMass, lambda1405::SigmaPlusMass, lambda1405::XiMinusMass,
                  lambda1405::PtSigma, lambda1405::AlphaAPSigma, lambda1405::QtAPSigma, lambda1405::RadiusSigma,
                  lambda1405::PtKink,
                  lambda1405::NSigmaTPCPiKink, lambda1405::NSigmaTOFPiKink,
                  lambda1405::NSigmaTPCPrKink, lambda1405::NSigmaTOFPrKink,
                  lambda1405::DCAKinkDauToPV,
                  lambda1405::NSigmaTPCPiDau, lambda1405::NSigmaTOFPiDau,
                  lambda1405::Centrality, lambda1405::Occupancy);

DECLARE_SOA_TABLE(Lambda1405Flow, "AOD", "LAMBDA1405FLOW",
                  o2::soa::Index<>,
                  lambda1405::Pt, lambda1405::Mass,
                  lambda1405::PtSigma,
                  lambda1405::SigmaMinusMass, lambda1405::SigmaPlusMass,
                  lambda1405::AlphaAPSigma, lambda1405::QtAPSigma,
                  lambda1405::NSigmaTPCPiKink, lambda1405::NSigmaTOFPiKink,
                  lambda1405::NSigmaTPCPrKink, lambda1405::NSigmaTOFPrKink,
                  lambda1405::DCAKinkDauToPV,
                  lambda1405::NSigmaTPCPiDau, lambda1405::NSigmaTOFPiDau,
                  lambda1405::ScalarProd, lambda1405::Centrality);

DECLARE_SOA_TABLE(Lambda1405CandsMC, "AOD", "MCLAMBDA1405",
                  o2::soa::Index<>,
                  lambda1405::Px, lambda1405::Py, lambda1405::Pz,
                  lambda1405::Mass, lambda1405::MassXi1530,
                  lambda1405::SigmaMinusMass, lambda1405::SigmaPlusMass, lambda1405::XiMinusMass,
                  lambda1405::PtSigma, lambda1405::AlphaAPSigma, lambda1405::QtAPSigma, lambda1405::RadiusSigma,
                  lambda1405::PtKink,
                  lambda1405::NSigmaTPCPiKink, lambda1405::NSigmaTOFPiKink,
                  lambda1405::NSigmaTPCPrKink, lambda1405::NSigmaTOFPrKink,
                  lambda1405::DCAKinkDauToPV,
                  lambda1405::NSigmaTPCPiDau, lambda1405::NSigmaTOFPiDau,
                  lambda1405::PtMC, lambda1405::MassMC, lambda1405::SigmaPdgCode, lambda1405::KinkDauPdgCode,
                  lambda1405::Centrality, lambda1405::Occupancy);

DECLARE_SOA_TABLE(Lambda1405SigmaEffMC, "AOD", "MCL1405SIGEFF",
                  o2::soa::Index<>,
                  lambda1405::PxSigma, lambda1405::DeltaPxSigma,
                  lambda1405::PySigma, lambda1405::DeltaPySigma,
                  lambda1405::PzSigma, lambda1405::DeltaPzSigma,
                  lambda1405::PtSigma, lambda1405::DeltaPtSigma,
                  lambda1405::RadiusSigma, lambda1405::DeltaRadiusSigma,
                  lambda1405::MassSigma, lambda1405::DeltaMassSigma,
                  lambda1405::DeltaPxSigmaRecalc,
                  lambda1405::DeltaPySigmaRecalc,
                  lambda1405::DeltaPzSigmaRecalc,
                  lambda1405::GenPhiSigma,
                  lambda1405::GenEtaSigma,
                  lambda1405::PxKinkDaug, lambda1405::DeltaPxKinkDaug,
                  lambda1405::PyKinkDaug, lambda1405::DeltaPyKinkDaug,
                  lambda1405::PzKinkDaug, lambda1405::DeltaPzKinkDaug,
                  lambda1405::PtKinkDaug, lambda1405::DeltaPtKinkDaug,
                  lambda1405::GenPhiKinkDaug,
                  lambda1405::GenEtaKinkDaug,
                  lambda1405::XKinkVtx, lambda1405::DeltaXKinkVtx,
                  lambda1405::YKinkVtx, lambda1405::DeltaYKinkVtx,
                  lambda1405::ZKinkVtx, lambda1405::DeltaZKinkVtx,
                  lambda1405::AlphaAPSigma, lambda1405::QtAPSigma,
                  lambda1405::AlphaAPSigmaRecalc, lambda1405::QtAPSigmaRecalc,
                  lambda1405::DCAKinkDauToPV, lambda1405::DCASigmaToPV,
                  lambda1405::SigmaPdgCode, lambda1405::KinkDauPdgCode,
                  lambda1405::Centrality, lambda1405::Occupancy, lambda1405::PvContrib);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFLAMBDA1405TABLE_H_
