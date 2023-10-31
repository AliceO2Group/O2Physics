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

/// \file MchTrkEffTables.h
/// \brief MCH tracking efficiency table definition.
/// \author  Zaida Conesa del Valle <zaida.conesa.del.valle@cern.ch>
///

#ifndef PWGDQ_DATAMODEL_MCHTRKEFFTABLES_H_
#define PWGDQ_DATAMODEL_MCHTRKEFFTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{

namespace mch_trk_eff
{
// Define the selection columns
DECLARE_SOA_COLUMN(Eta, eta, float);                //! reconstructed eta
DECLARE_SOA_COLUMN(Pt, pt, float);                  //! reconstructed pt
DECLARE_SOA_COLUMN(Phi, phi, float);                //! reconstructed phi
DECLARE_SOA_COLUMN(MchBitMap, mchBitMap, uint16_t); //! mch bit map
DECLARE_SOA_COLUMN(EtaGen, etaGen, float);          //! simulated eta
DECLARE_SOA_COLUMN(PtGen, ptGen, float);            //! simulated pt
DECLARE_SOA_COLUMN(PhiGen, phiGen, float);          //! simulated phi
} // namespace mch_trk_eff

// Table
DECLARE_SOA_TABLE(MchTrkEffBase, "AOD", "MCHTRKEFFBASE", //! table with muon track properties and mch bit map
                  mch_trk_eff::Eta, mch_trk_eff::Pt, mch_trk_eff::Phi, mch_trk_eff::MchBitMap);

DECLARE_SOA_TABLE(MchTrkEffGen, "AOD", "MCHTRKEFFGEN", //! table with simulated muon track properties
                  mch_trk_eff::EtaGen, mch_trk_eff::PtGen, mch_trk_eff::PhiGen);

} // namespace o2::aod

#endif // PWGDQ_DATAMODEL_MCHTRKEFFTABLES_H_
