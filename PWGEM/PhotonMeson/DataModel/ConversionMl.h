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
/// \file ConversionMl.h
/// \brief This header provides the table definitions to store cluster pairs to analyze and potentially tag late conversions
/// \author Marvin Hemmer (marvin.hemmer@cern.ch) - Goethe University Frankfurt

#ifndef PWGEM_PHOTONMESON_DATAMODEL_CONVERSIONML_H_
#define PWGEM_PHOTONMESON_DATAMODEL_CONVERSIONML_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{

namespace convtag
{
DECLARE_SOA_COLUMN(PMEvent, pmevent, int);          //! index column to the pmevent
DECLARE_SOA_COLUMN(Minv, minv, float);              //! invariant mass of the pair
DECLARE_SOA_COLUMN(HarmonicET, harmonicET, float);  //! ET1 * ET2 / (ET1 + ET2) in GeV
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);      //! difference in eta of the two clusters
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);      //! difference in phi of the two clusters in rad
DECLARE_SOA_COLUMN(Phiv, phiv, float);              //! angle between plane spanned by the two clusters and the plane between the B-field direction the their combined momentum
DECLARE_SOA_COLUMN(E1, e1, float);                  //! energy of first cluster in GeV
DECLARE_SOA_COLUMN(E2, e2, float);                  //! energy of second cluster in GeV
DECLARE_SOA_COLUMN(M021, m021, float);              //! M02 of first cluster
DECLARE_SOA_COLUMN(M022, m022, float);              //! M02 of second cluster
DECLARE_SOA_COLUMN(Time1, time1, float);            //! time of first cluster in ns
DECLARE_SOA_COLUMN(Time2, time2, float);            //! time of second cluster in ns
DECLARE_SOA_COLUMN(NCells1, nCells1, uint8_t);      //! NCells of first cluster
DECLARE_SOA_COLUMN(NCells2, nCells2, uint8_t);      //! NCells of second cluster
DECLARE_SOA_COLUMN(TruthLabel, truthLabel, int8_t); //! truth label for ML training -- see TruthClass enum for the mapping
DECLARE_SOA_COLUMN(CentOrMult, centOrMult, float);  //! centrality or multiplicity value of the collision
DECLARE_SOA_COLUMN(Purity1, purity1, float);        // leading-contributor amplitude fraction, cluster 1
DECLARE_SOA_COLUMN(Purity2, purity2, float);        // leading-contributor amplitude fraction, cluster 2

} // namespace convtag
DECLARE_SOA_TABLE(ConvTagCandidates, "AOD", "CONVTAGCAND",
                  convtag::PMEvent,
                  convtag::Minv, convtag::HarmonicET, convtag::DeltaEta, convtag::DeltaPhi, convtag::Phiv,
                  convtag::E1, convtag::E2, convtag::M021, convtag::M022,
                  convtag::Time1, convtag::Time2, convtag::NCells1, convtag::NCells2,
                  convtag::TruthLabel, convtag::CentOrMult);

DECLARE_SOA_TABLE_VERSIONED(ConvTagCandidates_001, "AOD", "CONVTAGCAND", 1,
                            convtag::PMEvent,
                            convtag::Minv, convtag::HarmonicET, convtag::DeltaEta, convtag::DeltaPhi, convtag::Phiv,
                            convtag::E1, convtag::E2, convtag::M021, convtag::M022,
                            convtag::Time1, convtag::Time2, convtag::NCells1, convtag::NCells2, convtag::Purity1, convtag::Purity2,
                            convtag::TruthLabel, convtag::CentOrMult);

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_CONVERSIONML_H_
