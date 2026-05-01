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
/// \file LFPhotonDeuteronTables.h
/// \brief Tables for photon-deuteron correlation analysis
/// \author Arvind Khuntia <arvind.khuntia@cern.ch> and Francesco Noferini <francesco.noferini@cern.ch>

#ifndef PWGLF_DATAMODEL_LFPHOTONDEUTERONTABLES_H_
#define PWGLF_DATAMODEL_LFPHOTONDEUTERONTABLES_H_

#include <Framework/ASoA.h>

namespace o2::aod
{
namespace photondeuteron
{
DECLARE_SOA_COLUMN(PhotonPt, photonPt, float);                         //! Photon transverse momentum
DECLARE_SOA_COLUMN(PhotonEta, photonEta, float);                       //! Photon pseudorapidity
DECLARE_SOA_COLUMN(PhotonPhi, photonPhi, float);                       //! Photon azimuthal angle
DECLARE_SOA_COLUMN(PhotonMass, photonMass, float);                     //! Photon invariant mass
DECLARE_SOA_COLUMN(PhotonPosPt, photonPosPt, float);                   //! Positive daughter (e+) transverse momentum
DECLARE_SOA_COLUMN(PhotonPosEta, photonPosEta, float);                 //! Positive daughter (e+) pseudorapidity
DECLARE_SOA_COLUMN(PhotonPosPhi, photonPosPhi, float);                 //! Positive daughter (e+) azimuthal angle
DECLARE_SOA_COLUMN(PhotonPosNSigmaElTPC, photonPosNSigmaElTPC, float); //! Positive daughter (e+) NSigma electron TPC
DECLARE_SOA_COLUMN(PhotonNegPt, photonNegPt, float);                   //! Negative daughter (e-) transverse momentum
DECLARE_SOA_COLUMN(PhotonNegEta, photonNegEta, float);                 //! Negative daughter (e-) pseudorapidity
DECLARE_SOA_COLUMN(PhotonNegPhi, photonNegPhi, float);                 //! Negative daughter (e-) azimuthal angle
DECLARE_SOA_COLUMN(PhotonNegNSigmaElTPC, photonNegNSigmaElTPC, float); //! Negative daughter (e-) NSigma electron TPC
DECLARE_SOA_COLUMN(DeuteronPt, deuteronPt, float);                     //! Deuteron transverse momentum
DECLARE_SOA_COLUMN(DeuteronEta, deuteronEta, float);                   //! Deuteron pseudorapidity
DECLARE_SOA_COLUMN(DeuteronPhi, deuteronPhi, float);                   //! Deuteron azimuthal angle
DECLARE_SOA_COLUMN(DeuteronNSigmaTPC, deuteronNSigmaTPC, float);       //! Deuteron NSigma TPC
DECLARE_SOA_COLUMN(DeuteronNSigmaTOF, deuteronNSigmaTOF, float);       //! Deuteron NSigma TOF
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);                         //! Delta phi between photon and deuteron
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);                         //! Delta eta between photon and deuteron
DECLARE_SOA_COLUMN(InvMass, invMass, float);                           //! Photon-deuteron invariant mass
DECLARE_SOA_COLUMN(RelativeMomentum, relativeMomentum, float);         //! Relative momentum k*_pn
} // namespace photondeuteron

DECLARE_SOA_TABLE(PhotonDeuteronPairs, "AOD", "PHOTONDPAIRS", //! Table for photon-deuteron pairs
                  photondeuteron::PhotonPt,
                  photondeuteron::PhotonEta,
                  photondeuteron::PhotonPhi,
                  photondeuteron::PhotonMass,
                  photondeuteron::PhotonPosPt,
                  photondeuteron::PhotonPosEta,
                  photondeuteron::PhotonPosPhi,
                  photondeuteron::PhotonPosNSigmaElTPC,
                  photondeuteron::PhotonNegPt,
                  photondeuteron::PhotonNegEta,
                  photondeuteron::PhotonNegPhi,
                  photondeuteron::PhotonNegNSigmaElTPC,
                  photondeuteron::DeuteronPt,
                  photondeuteron::DeuteronEta,
                  photondeuteron::DeuteronPhi,
                  photondeuteron::DeuteronNSigmaTPC,
                  photondeuteron::DeuteronNSigmaTOF,
                  photondeuteron::DeltaPhi,
                  photondeuteron::DeltaEta,
                  photondeuteron::InvMass,
                  photondeuteron::RelativeMomentum);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPHOTONDEUTERONTABLES_H_
