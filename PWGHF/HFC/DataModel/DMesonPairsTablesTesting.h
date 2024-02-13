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

/// \file DMesonPairsTablesTesting.h
/// \brief D meson pair table definition.
/// \note Temporary code for tests of D0-D0 correlations
///
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#ifndef PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLESTESTING_H_
#define PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLESTESTING_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables for D(bar) Meson correlation pair studies
namespace hf_correlation_d_meson_pair_testing
{
// Kinematic info
DECLARE_SOA_COLUMN(PtCand1, ptCand1, float); //! Transverse momentum of first candidate
DECLARE_SOA_COLUMN(PtCand2, ptCand2, float); //! Transverse momentum of second candidate
DECLARE_SOA_COLUMN(YCand1, yCand1, float);   //! Rapidity of first candidate
DECLARE_SOA_COLUMN(YCand2, yCand2, float);   //! Rapidity of second candidate
// Invariant mass
DECLARE_SOA_COLUMN(MDCand1, mDCand1, float);       //! Invariant mass of first candidate as D
DECLARE_SOA_COLUMN(MDbarCand1, mDbarCand1, float); //! Invariant mass of first candidate as Dbar
DECLARE_SOA_COLUMN(MDCand2, mDCand2, float);       //! Invariant mass of second candidate as D
DECLARE_SOA_COLUMN(MDbarCand2, mDbarCand2, float); //! Invariant mass of second candidate as Dbar
DECLARE_SOA_COLUMN(PairType, pairType, uint8_t);   //! Bitmap with all pair types (DD, DDbar, etc.) a pair of candidates has passed
// MC info
DECLARE_SOA_COLUMN(Origin1, origin1, uint8_t);       //! candidate 1 origin
DECLARE_SOA_COLUMN(Origin2, origin2, uint8_t);       //! candidate 2 origin
DECLARE_SOA_COLUMN(MatchedMc1, matchedMc1, uint8_t); //! MC matching of candidate 1
DECLARE_SOA_COLUMN(MatchedMc2, matchedMc2, uint8_t); //! MC matching of candidate 2
} // namespace hf_correlation_d_meson_pair_testing

// Definition of the D meson pair table. Contains the info needed at Data level.
#define DECLARE_DMESON_PAIR_TABLE_TESTING(_pair_type_, _marker_value_, _description_)   \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_, o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair_testing::PtCand1,                       \
                    hf_correlation_d_meson_pair_testing::PtCand2,                       \
                    hf_correlation_d_meson_pair_testing::YCand1,                        \
                    hf_correlation_d_meson_pair_testing::YCand2,                        \
                    hf_correlation_d_meson_pair_testing::MDCand1,                       \
                    hf_correlation_d_meson_pair_testing::MDbarCand1,                    \
                    hf_correlation_d_meson_pair_testing::MDCand2,                       \
                    hf_correlation_d_meson_pair_testing::MDbarCand2,                    \
                    hf_correlation_d_meson_pair_testing::PairType);
// Definition of the D meson pair table with info at MC level.
#define DECLARE_DMESON_PAIR_MCINFO_TABLE_TESTING(_pair_type_, _marker_value_, _description_)     \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_ "MCINFO", o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair_testing::Origin1,                                \
                    hf_correlation_d_meson_pair_testing::Origin2,                                \
                    hf_correlation_d_meson_pair_testing::MatchedMc1,                             \
                    hf_correlation_d_meson_pair_testing::MatchedMc2);

// Creation of tables with D Meson Pairs info
DECLARE_DMESON_PAIR_TABLE_TESTING(D0PairTesting, 1, "D0PAIR");              //! D0 pairs Info
DECLARE_DMESON_PAIR_MCINFO_TABLE_TESTING(D0PairMcInfoTesting, 1, "D0PAIR"); //! D0 pairs MC Rec Info

DECLARE_DMESON_PAIR_TABLE_TESTING(D0PairMcGenTesting, 2, "D0PAIRGEN");            //! D0 pairs MC Gen Kinematic Info
DECLARE_DMESON_PAIR_MCINFO_TABLE_TESTING(D0PairMcGenInfoTesting, 2, "D0PAIRGEN"); //! D0 pairs MC Gen Info

} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLESTESTING_H_
