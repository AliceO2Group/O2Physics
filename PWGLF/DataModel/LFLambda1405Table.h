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

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFLAMBDA1405TABLES_H_
#define PWGLF_DATAMODEL_LFLAMBDA1405TABLES_H_

namespace o2::aod
{

namespace lambda1405
{

DECLARE_SOA_COLUMN(Px, px, float);                           //! Px of the candidate
DECLARE_SOA_COLUMN(Py, py, float);                           //! Py of the candidate
DECLARE_SOA_COLUMN(Pz, pz, float);                           //! Pz of the candidate
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
DECLARE_SOA_COLUMN(NSigmaTPCPiDau, nSigmaTPCPiDau, float);   //! Number of sigmas for the lambda1405 pion daughter in TPC
DECLARE_SOA_COLUMN(NSigmaTOFPiDau, nSigmaTOFPiDau, float);   //! Number of sigmas for the lambda1405 pion daughter in TOF

// MC Columns
DECLARE_SOA_COLUMN(PtMC, ptMC, float);                   //! pT of the candidate in MC
DECLARE_SOA_COLUMN(MassMC, massMC, float);               //! Invariant mass of the candidate in MC
DECLARE_SOA_COLUMN(SigmaPdgCode, sigmaPdgCode, int);     //! PDG code of the Sigma daughter
DECLARE_SOA_COLUMN(KinkDauPdgCode, kinkDauPdgCode, int); //! PDG code of the kink daughter

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
                  lambda1405::NSigmaTPCPiDau, lambda1405::NSigmaTOFPiDau);

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
                  lambda1405::PtMC, lambda1405::MassMC, lambda1405::SigmaPdgCode, lambda1405::KinkDauPdgCode);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFLAMBDA1405TABLES_H_
