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
/// \file   mcCentrality.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \since  2024-06-05
/// \brief  Set of tables for the MC centrality
///

#ifndef PWGLF_DATAMODEL_MCCENTRALITY_H_
#define PWGLF_DATAMODEL_MCCENTRALITY_H_

// O2 includes
#include "Common/DataModel/Centrality.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Logger.h"

namespace o2::aod
{
// DECLARE_SOA_TABLE(CentFV0As, "AOD", "CENTFV0A", cent::CentFV0A); //! Run3 FV0A estimated centrality table
// DECLARE_SOA_TABLE(CentFT0Ms, "AOD", "CENTFT0M", cent::CentFT0M); //! Run3 FT0M estimated centrality table
// DECLARE_SOA_TABLE(CentFT0As, "AOD", "CENTFT0A", cent::CentFT0A); //! Run3 FT0A estimated centrality table
// DECLARE_SOA_TABLE(CentFT0Cs, "AOD", "CENTFT0C", cent::CentFT0C); //! Run3 FT0C estimated centrality table
// DECLARE_SOA_TABLE(CentFDDMs, "AOD", "CENTFDDM", cent::CentFDDM); //! Run3 FDDM estimated centrality table
// DECLARE_SOA_TABLE(CentNTPVs, "AOD", "CENTNTPV", cent::CentNTPV); //! Run3 NTPV estimated centrality table

DECLARE_SOA_TABLE(McCentFV0As, "AOD", "MCCENTFV0A", o2::soa::Marker<1>, cent::CentFV0A);
DECLARE_SOA_TABLE(McCentFT0Ms, "AOD", "MCCENTFT0M", o2::soa::Marker<2>, cent::CentFT0M);
DECLARE_SOA_TABLE(McCentFT0As, "AOD", "MCCENTFT0A", o2::soa::Marker<3>, cent::CentFT0A);
DECLARE_SOA_TABLE(McCentFT0Cs, "AOD", "MCCENTFT0C", o2::soa::Marker<4>, cent::CentFT0C);
DECLARE_SOA_TABLE(McCentFDDMs, "AOD", "MCCENTFDDM", o2::soa::Marker<5>, cent::CentFDDM);
DECLARE_SOA_TABLE(McCentNTPVs, "AOD", "MCCENTNTPV", o2::soa::Marker<6>, cent::CentNTPV);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_MCCENTRALITY_H_
