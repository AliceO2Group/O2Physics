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
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#ifndef PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLES_H_
#define PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables for D(bar) Meson correlation pair studies
namespace hf_correlation_d_meson_pair
{
// Kinematic info
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float); //! Delta phi of the pair
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float); //! Delta eta of the pair
DECLARE_SOA_COLUMN(PtCand1, ptCand1, float);   //! Transverse momentum of first candidate
DECLARE_SOA_COLUMN(PtCand2, ptCand2, float);   //! Transverse momentum of second candidate
DECLARE_SOA_COLUMN(YCand1, yCand1, float);     //! Rapidity of first candidate
DECLARE_SOA_COLUMN(YCand2, yCand2, float);     //! Rapidity of second candidate
// Invariant mass
DECLARE_SOA_COLUMN(MCand1, mCand1, float); //! Invariant mass of first candidate
DECLARE_SOA_COLUMN(MCand2, mCand2, float); //! Invariant mass of second candidate
// Type of candidate: candidate charge and whether it is signal, reflected, or bkg
DECLARE_SOA_COLUMN(CandidateType1, candidateType1, uint8_t); //! Type of first candidate
DECLARE_SOA_COLUMN(CandidateType2, candidateType2, uint8_t); //! Type of second candidate
DECLARE_SOA_COLUMN(DataType, dataType, uint8_t);             //! 0: data, 1: MC reco, 2: MC gen
// MC info
DECLARE_SOA_COLUMN(Origin1, origin1, uint8_t);       //! candidate 1 origin
DECLARE_SOA_COLUMN(Origin2, origin2, uint8_t);       //! candidate 2 origin
DECLARE_SOA_COLUMN(MatchedMc1, matchedMc1, uint8_t); //! MC matching of candidate 1
DECLARE_SOA_COLUMN(MatchedMc2, matchedMc2, uint8_t); //! MC matching of candidate 2
} // namespace hf_correlation_d_meson_pair

// Definition of the D meson pair table. Contains the info needed at Data level.
#define DECLARE_DMESON_PAIR_TABLE(_pair_type_, _marker_value_, _description_)           \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_, o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair::DeltaPhi,                              \
                    hf_correlation_d_meson_pair::DeltaEta,                              \
                    hf_correlation_d_meson_pair::PtCand1,                               \
                    hf_correlation_d_meson_pair::PtCand2,                               \
                    hf_correlation_d_meson_pair::YCand1,                                \
                    hf_correlation_d_meson_pair::YCand2,                                \
                    hf_correlation_d_meson_pair::MCand1,                                \
                    hf_correlation_d_meson_pair::MCand2,                                \
                    hf_correlation_d_meson_pair::CandidateType1,                        \
                    hf_correlation_d_meson_pair::CandidateType2,                        \
                    hf_correlation_d_meson_pair::DataType);

// definition of the D meson pair table with info at MC level.
#define DECLARE_DMESON_PAIR_RECOINFO_TABLE(_pair_type_, _marker_value_, _description_)             \
  DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_ "RECOINFO", o2::soa::Marker<_marker_value_>, \
                    hf_correlation_d_meson_pair::Origin1,                                          \
                    hf_correlation_d_meson_pair::Origin2,                                          \
                    hf_correlation_d_meson_pair::MatchedMc1,                                       \
                    hf_correlation_d_meson_pair::MatchedMc2);

// Creation of tables with D Meson Pairs info
DECLARE_DMESON_PAIR_TABLE(D0Pair, 1, "D0PAIR");                    //! D0 pairs Info
DECLARE_DMESON_PAIR_RECOINFO_TABLE(D0PairRecoInfo, 1, "D0PAIR");   //! D0 pairs Reconstructed Info
DECLARE_DMESON_PAIR_TABLE(DplusPair, 2, "DPLUSPAIR");              //! D+ pairs Info
DECLARE_DMESON_PAIR_RECOINFO_TABLE(DplusPairRecoInfo, 2, "DPLUS"); //! D+ pairs Reconstructed Info

} // namespace o2::aod

#endif // PWGHF_HFC_DATAMODEL_DMESONPAIRSTABLES_H_
