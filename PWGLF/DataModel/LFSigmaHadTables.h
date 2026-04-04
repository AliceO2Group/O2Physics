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
/// \file LFSigmaHadTables.h
/// \brief Slim tables for Sigma-hadron pairs
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>
///

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include <Framework/AnalysisDataModel.h>

#ifndef PWGLF_DATAMODEL_LFSIGMAHADTABLES_H_
#define PWGLF_DATAMODEL_LFSIGMAHADTABLES_H_

namespace o2::aod
{

namespace sigmaproton
{
DECLARE_SOA_COLUMN(ChargeSigma, chargeSigma, int);     //! Charge of the sigma candidate
DECLARE_SOA_COLUMN(SigmaDecRad, sigmaDecRad, float);   //! Decay radius of the Sigma candidate
DECLARE_SOA_COLUMN(SigmaCosPA, sigmaCosPA, float);     //! Cosine of pointing angle of the Sigma candidate
DECLARE_SOA_COLUMN(ChargeHad, chargeHad, int);         //! Charge of the hadron candidate
DECLARE_SOA_COLUMN(PxHad, pxHad, float);               //! Px of the hadron candidate
DECLARE_SOA_COLUMN(PyHad, pyHad, float);               //! Py of the hadron candidate
DECLARE_SOA_COLUMN(PzHad, pzHad, float);               //! Pz of the hadron candidate
DECLARE_SOA_COLUMN(NSigmaTPCHad, nSigmaTPCHad, float); //! Number of sigmas for the hadron candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTOFHad, nSigmaTOFHad, float); //! Number of sigmas for the hadron candidate from Sigma kink in TOF

// MC Columns
DECLARE_SOA_COLUMN(SigmaPDG, sigmaPDG, int);       //! PDG code of the Sigma daughter
DECLARE_SOA_COLUMN(DaughterPDG, daughterPDG, int); //! PDG code of the kink daughter
DECLARE_SOA_COLUMN(HadPDG, hadPDG, int);           //! PDG code of the hadron candidate
DECLARE_SOA_COLUMN(SigmaGenPt, sigmaGenPt, float); //! Generated pT of the Sigma candidate
DECLARE_SOA_COLUMN(HadGenPt, hadGenPt, float);     //! Generated pT of the hadron candidate
DECLARE_SOA_COLUMN(GenKStar, genKStar, float);     //! Generated k* of the Sigma-hadron pair

} // namespace sigmaproton

DECLARE_SOA_TABLE(SigmaProtonCands, "AOD", "SIGMAPROTONCANDS",
                  o2::soa::Index<>,
                  sigmaproton::ChargeSigma, kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug, sigmaproton::SigmaDecRad, sigmaproton::SigmaCosPA,
                  sigmaproton::ChargeHad, sigmaproton::PxHad, sigmaproton::PyHad, sigmaproton::PzHad,
                  sigmaproton::NSigmaTPCHad, sigmaproton::NSigmaTOFHad);

DECLARE_SOA_TABLE(SigmaProtonMCCands, "AOD", "SIGMAPROTONMCCANDS",
                  o2::soa::Index<>,
                  sigmaproton::ChargeSigma, kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug, sigmaproton::SigmaDecRad, sigmaproton::SigmaCosPA,
                  sigmaproton::ChargeHad, sigmaproton::PxHad, sigmaproton::PyHad, sigmaproton::PzHad,
                  sigmaproton::NSigmaTPCHad, sigmaproton::NSigmaTOFHad,
                  sigmaproton::SigmaPDG, sigmaproton::DaughterPDG, sigmaproton::HadPDG,
                  sigmaproton::SigmaGenPt, sigmaproton::HadGenPt, sigmaproton::GenKStar);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSIGMAHADTABLES_H_
