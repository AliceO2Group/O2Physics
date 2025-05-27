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

/// \file DMesonPairsTables.h
/// \brief D meson pair table definition.
///
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#ifndef PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLES_H_
#define PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLES_H_

#include <vector>

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables for D(bar) Meson correlation pair studies
namespace hf_correlation_d_meson_pair
{
// Kinematic info
DECLARE_SOA_COLUMN(PtCand1, ptCand1, float);   //! Transverse momentum of first candidate
DECLARE_SOA_COLUMN(PtCand2, ptCand2, float);   //! Transverse momentum of second candidate
DECLARE_SOA_COLUMN(YCand1, yCand1, float);     //! Rapidity of first candidate
DECLARE_SOA_COLUMN(YCand2, yCand2, float);     //! Rapidity of second candidate
DECLARE_SOA_COLUMN(PhiCand1, phiCand1, float); //! Azimuthal angle of first candidate
DECLARE_SOA_COLUMN(PhiCand2, phiCand2, float); //! Azimuthal angle of second candidate
// Invariant mass
DECLARE_SOA_COLUMN(MDCand1, mDCand1, float);                 //! Invariant mass of first candidate as D
DECLARE_SOA_COLUMN(MDbarCand1, mDbarCand1, float);           //! Invariant mass of first candidate as Dbar
DECLARE_SOA_COLUMN(MDCand2, mDCand2, float);                 //! Invariant mass of second candidate as D
DECLARE_SOA_COLUMN(MDbarCand2, mDbarCand2, float);           //! Invariant mass of second candidate as Dbar
DECLARE_SOA_COLUMN(PairType, pairType, uint8_t);             //! Bitmap with all pair types (DD, DDbar, etc.) a pair of candidates has passed
DECLARE_SOA_COLUMN(CandidateType1, candidateType1, uint8_t); //! Bitmap with Selected and True info of candidate 1
DECLARE_SOA_COLUMN(CandidateType2, candidateType2, uint8_t); //! Bitmap with Selected and True info of candidate 2
// MC info
DECLARE_SOA_COLUMN(Origin1, origin1, uint8_t);       //! candidate 1 origin
DECLARE_SOA_COLUMN(Origin2, origin2, uint8_t);       //! candidate 2 origin
DECLARE_SOA_COLUMN(MatchedMc1, matchedMc1, uint8_t); //! MC matching of candidate 1
DECLARE_SOA_COLUMN(MatchedMc2, matchedMc2, uint8_t); //! MC matching of candidate 2
// ML info
DECLARE_SOA_COLUMN(MlProbD0Cand1, mlProbD0Cand1, std::vector<float>);       //!
DECLARE_SOA_COLUMN(MlProbD0barCand1, mlProbD0barCand1, std::vector<float>); //!
DECLARE_SOA_COLUMN(MlProbD0Cand2, mlProbD0Cand2, std::vector<float>);       //!
DECLARE_SOA_COLUMN(MlProbD0barCand2, mlProbD0barCand2, std::vector<float>); //!
} // namespace hf_correlation_d_meson_pair

// Definition of the D meson pair table. Contains the info needed at Data level.
#define DECLARE_DMESON_PAIR_TABLE(_pair_type_, _marker_value_, _description_)           \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_, o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair::PtCand1,                               \
                    hf_correlation_d_meson_pair::PtCand2,                               \
                    hf_correlation_d_meson_pair::YCand1,                                \
                    hf_correlation_d_meson_pair::YCand2,                                \
                    hf_correlation_d_meson_pair::PhiCand1,                              \
                    hf_correlation_d_meson_pair::PhiCand2,                              \
                    hf_correlation_d_meson_pair::MDCand1,                               \
                    hf_correlation_d_meson_pair::MDbarCand1,                            \
                    hf_correlation_d_meson_pair::MDCand2,                               \
                    hf_correlation_d_meson_pair::MDbarCand2,                            \
                    hf_correlation_d_meson_pair::PairType,                              \
                    hf_correlation_d_meson_pair::CandidateType1,                        \
                    hf_correlation_d_meson_pair::CandidateType2);
// Definition of the D meson pair table with info at MC level.
#define DECLARE_DMESON_PAIR_MCINFO_TABLE(_pair_type_, _marker_value_, _description_)             \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_ "MCINFO", o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair::Origin1,                                        \
                    hf_correlation_d_meson_pair::Origin2,                                        \
                    hf_correlation_d_meson_pair::MatchedMc1,                                     \
                    hf_correlation_d_meson_pair::MatchedMc2);
// Definition of the table with the ML info of the D meson pair.
#define DECLARE_DMESON_PAIR_MLINFO_TABLE(_pair_type_, _marker_value_, _description_)         \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_ "ML", o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair::MlProbD0Cand1,                              \
                    hf_correlation_d_meson_pair::MlProbD0barCand1,                           \
                    hf_correlation_d_meson_pair::MlProbD0Cand2,                              \
                    hf_correlation_d_meson_pair::MlProbD0barCand2);

// Creation of tables with D Meson Pairs info
DECLARE_DMESON_PAIR_TABLE(D0Pair, 1, "D0PAIR");              //! D0 pairs Info
DECLARE_DMESON_PAIR_MCINFO_TABLE(D0PairMcInfo, 1, "D0PAIR"); //! D0 pairs MC Rec Info
DECLARE_DMESON_PAIR_MLINFO_TABLE(D0PairMl, 1, "D0PAIR");     //! D0 pairs ML Info

DECLARE_DMESON_PAIR_TABLE(D0PairMcGen, 2, "D0PAIRGEN");            //! D0 pairs MC Gen Kinematic Info
DECLARE_DMESON_PAIR_MCINFO_TABLE(D0PairMcGenInfo, 2, "D0PAIRGEN"); //! D0 pairs MC Gen Info

} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLES_H_
