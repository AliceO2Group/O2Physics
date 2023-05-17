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

#ifndef PWGHF_DATAMODEL_DMESONPAIRTABLES_H_
#define PWGHF_DATAMODEL_DMESONPAIRTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
// definition of columns and tables for D(bar) Meson correlation pair studies
namespace hf_correlation_d_meson_pair
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(PtCand1, ptCand1, float);
DECLARE_SOA_COLUMN(PtCand2, ptCand2, float);
DECLARE_SOA_COLUMN(YCand1, yCand1, float);
DECLARE_SOA_COLUMN(YCand2, yCand2, float);
DECLARE_SOA_COLUMN(MCand1, mCand1, float);
DECLARE_SOA_COLUMN(MCand2, mCand2, float);
DECLARE_SOA_COLUMN(CandidateType1, candidateType1, int);
DECLARE_SOA_COLUMN(CandidateType2, candidateType2, int);
DECLARE_SOA_COLUMN(Origin1, origin1, uint8_t);
DECLARE_SOA_COLUMN(Origin2, origin2, uint8_t);
DECLARE_SOA_COLUMN(MatchedMc1, matchedMc1, uint8_t);
DECLARE_SOA_COLUMN(MatchedMc2, matchedMc2, uint8_t);
DECLARE_SOA_COLUMN(DataType, dataType, uint8_t); // 0: data, 1: MC reco, 2: MC gen
} // namespace hf_correlation_d_meson_pair

// Definition of the D meson pair table. Contains the info needed at Data level.
#define DECLARE_DMESON_PAIR_TABLE(_pair_type_, _marker_value_, _description_)                                      \
    DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_, o2::soa::Marker<_marker_value_>,                          \
                      hf_correlation_d_meson_pair::DeltaPhi,                                                       \
                      hf_correlation_d_meson_pair::DeltaEta,                                                       \
                      hf_correlation_d_meson_pair::PtCand1,                                                        \
                      hf_correlation_d_meson_pair::PtCand2,                                                        \
                      hf_correlation_d_meson_pair::YCand1,                                                         \
                      hf_correlation_d_meson_pair::YCand2,                                                         \
                      hf_correlation_d_meson_pair::MCand1,                                                         \
                      hf_correlation_d_meson_pair::MCand2,                                                         \
                      hf_correlation_d_meson_pair::CandidateType1,                                                 \
                      hf_correlation_d_meson_pair::CandidateType2,                                                 \
                      hf_correlation_d_meson_pair::DataType);

// definition of the D meson pair table with info at MC level.
#define DECLARE_DMESON_PAIR_RECOINFO_TABLE(_pair_type_, _marker_value_, _description_)                             \
    DECLARE_SOA_TABLE(_pair_type_, "AOD", _description_ "RECOINFO", o2::soa::Marker<_marker_value_>,               \
                      hf_correlation_d_meson_pair::Origin1,                                                        \
                      hf_correlation_d_meson_pair::Origin2,                                                        \
                      hf_correlation_d_meson_pair::MatchedMc1,                                                     \
                      hf_correlation_d_meson_pair::MatchedMc2);
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DMESONPAIRTABLES_H_
