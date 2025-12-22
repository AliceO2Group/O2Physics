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
/// \file LFKinkDecayTables.h
/// \brief Slim tables for kinks
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>
///

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFSIGMAPROTONTABLES_H_
#define PWGLF_DATAMODEL_LFSIGMAPROTONTABLES_H_

namespace o2::aod
{

namespace sigmaproton
{
DECLARE_SOA_COLUMN(ChargeSigma, chargeSigma, int);   //! Charge of the sigma candidate
DECLARE_SOA_COLUMN(SigmaDecRad, sigmaDecRad, float); //! Decay radius of the Sigma candidate
DECLARE_SOA_COLUMN(SigmaCosPA, sigmaCosPA, float);   //! Cosine of pointing angle of the Sigma candidate
DECLARE_SOA_COLUMN(ChargePr, chargePr, int);         //! Charge of the proton candidate
DECLARE_SOA_COLUMN(PxPr, pxPr, float);               //! Px of the proton candidate
DECLARE_SOA_COLUMN(PyPr, pyPr, float);               //! Py of the proton candidate
DECLARE_SOA_COLUMN(PzPr, pzPr, float);               //! Pz of the proton candidate
DECLARE_SOA_COLUMN(NSigmaTPCPr, nSigmaTPCPr, float); //! Number of sigmas for the proton candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTOFPr, nSigmaTOFPr, float); //! Number of sigmas for the proton candidate from Sigma kink in TOF

// MC Columns
DECLARE_SOA_COLUMN(SigmaPDG, sigmaPDG, int);       //! PDG code of the Sigma daughter
DECLARE_SOA_COLUMN(DaughterPDG, daughterPDG, int); //! PDG code of the kink daughter
DECLARE_SOA_COLUMN(PrPDG, prPDG, int);             //! PDG code of the proton candidate

} // namespace sigmaproton

DECLARE_SOA_TABLE(SigmaProtonCands, "AOD", "SIGMAPROTONCANDS",
                  o2::soa::Index<>,
                  sigmaproton::ChargeSigma, kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug, sigmaproton::SigmaDecRad, sigmaproton::SigmaCosPA,
                  sigmaproton::ChargePr, sigmaproton::PxPr, sigmaproton::PyPr, sigmaproton::PzPr,
                  sigmaproton::NSigmaTPCPr, sigmaproton::NSigmaTOFPr);

DECLARE_SOA_TABLE(SigmaProtonMCCands, "AOD", "SIGMAPROTONMCCANDS",
                  o2::soa::Index<>,
                  sigmaproton::ChargeSigma, kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug, sigmaproton::SigmaDecRad, sigmaproton::SigmaCosPA,
                  sigmaproton::ChargePr, sigmaproton::PxPr, sigmaproton::PyPr, sigmaproton::PzPr,
                  sigmaproton::NSigmaTPCPr, sigmaproton::NSigmaTOFPr,
                  sigmaproton::SigmaPDG, sigmaproton::DaughterPDG, sigmaproton::PrPDG);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMAPROTONTABLES_H_
