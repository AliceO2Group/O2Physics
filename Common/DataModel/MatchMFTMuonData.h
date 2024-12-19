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
#ifndef COMMON_DATAMODEL_MATCHMFTMUONDATA_H_
#define COMMON_DATAMODEL_MATCHMFTMUONDATA_H_
#include "Framework/AnalysisDataModel.h"
#endif // COMMON_DATAMODEL_MATCHMFTMUONDATA_H_

namespace o2::aod
{
namespace matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, mDeltaPt, float);
DECLARE_SOA_COLUMN(DeltaEta, mDeltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, mDeltaPhi, float);
DECLARE_SOA_COLUMN(DeltaX, mDeltaX, float);
DECLARE_SOA_COLUMN(DeltaY, mDeltaY, float);
DECLARE_SOA_COLUMN(GMuonPt, mGMuonPt, float);
DECLARE_SOA_COLUMN(GMuonEta, mGMuonEta, float);
DECLARE_SOA_COLUMN(PairQ, mPairQ, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, mIsCorrectMatch, bool);
} // namespace matching_params

DECLARE_SOA_TABLE(MatchParams, "AOD", "MATCHING",
                  matching_params::GMuonPt,
                  matching_params::GMuonEta,
                  matching_params::PairQ,
                  matching_params::DeltaPt,
                  matching_params::DeltaX,
                  matching_params::DeltaY,
                  matching_params::DeltaEta,
                  matching_params::DeltaPhi,
                  matching_params::IsCorrectMatch);

namespace tag_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, mDeltaPt, float);
DECLARE_SOA_COLUMN(DeltaEta, mDeltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, mDeltaPhi, float);
DECLARE_SOA_COLUMN(DeltaX, mDeltaX, float);
DECLARE_SOA_COLUMN(DeltaY, mDeltaY, float);
DECLARE_SOA_COLUMN(GMuonPt, mGMuonPt, float);
DECLARE_SOA_COLUMN(GMuonEta, mGMuonEta, float);
DECLARE_SOA_COLUMN(PairQ, mPairQ, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, mIsCorrectMatch, bool);
} // namespace tag_matching_params

DECLARE_SOA_TABLE(TagMatchParams, "AOD", "TAGMATCHING",
                  tag_matching_params::GMuonPt,
                  tag_matching_params::GMuonEta,
                  tag_matching_params::PairQ,
                  tag_matching_params::DeltaPt,
                  tag_matching_params::DeltaX,
                  tag_matching_params::DeltaY,
                  tag_matching_params::DeltaEta,
                  tag_matching_params::DeltaPhi,
                  tag_matching_params::IsCorrectMatch);

namespace probe_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, mDeltaPt, float);
DECLARE_SOA_COLUMN(DeltaEta, mDeltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, mDeltaPhi, float);
DECLARE_SOA_COLUMN(DeltaX, mDeltaX, float);
DECLARE_SOA_COLUMN(DeltaY, mDeltaY, float);
DECLARE_SOA_COLUMN(TagGMuonPt, mTagGMuonPt, float);
DECLARE_SOA_COLUMN(GMuonPt, mGMuonPt, float);
DECLARE_SOA_COLUMN(GMuonEta, mGMuonEta, float);
DECLARE_SOA_COLUMN(PairQ, mPairQ, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, mIsCorrectMatch, bool);
} // namespace probe_matching_params

DECLARE_SOA_TABLE(ProbeMatchParams, "AOD", "PROBEMATCHING",
                  probe_matching_params::TagGMuonPt,
                  probe_matching_params::GMuonPt,
                  probe_matching_params::GMuonEta,
                  probe_matching_params::PairQ,
                  probe_matching_params::DeltaPt,
                  probe_matching_params::DeltaX,
                  probe_matching_params::DeltaY,
                  probe_matching_params::DeltaEta,
                  probe_matching_params::DeltaPhi,
                  probe_matching_params::IsCorrectMatch);

namespace mix_matching_params
{
// matching parameters
DECLARE_SOA_COLUMN(DeltaPt, mDeltaPt, float);
DECLARE_SOA_COLUMN(DeltaEta, mDeltaEta, float);
DECLARE_SOA_COLUMN(DeltaPhi, mDeltaPhi, float);
DECLARE_SOA_COLUMN(DeltaX, mDeltaX, float);
DECLARE_SOA_COLUMN(DeltaY, mDeltaY, float);
DECLARE_SOA_COLUMN(GMuonPt, mGMuonPt, float);
DECLARE_SOA_COLUMN(GMuonEta, mGMuonEta, float);
DECLARE_SOA_COLUMN(PairQ, mPairQ, int16_t);
DECLARE_SOA_COLUMN(IsCorrectMatch, mIsCorrectMatch, bool);
} // namespace mix_matching_params

DECLARE_SOA_TABLE(MixMatchParams, "AOD", "MIXMATCHING",
                  mix_matching_params::GMuonPt,
                  mix_matching_params::GMuonEta,
                  mix_matching_params::PairQ,
                  mix_matching_params::DeltaPt,
                  mix_matching_params::DeltaX,
                  mix_matching_params::DeltaY,
                  mix_matching_params::DeltaEta,
                  mix_matching_params::DeltaPhi,
                  mix_matching_params::IsCorrectMatch);

namespace muon_pair
{
// matching parameters
DECLARE_SOA_COLUMN(Mass, mMass, float);
DECLARE_SOA_COLUMN(Pt, mPt, float);
DECLARE_SOA_COLUMN(Rap, mRap, float);
DECLARE_SOA_COLUMN(PairQ, mPairQ, int16_t);
} // namespace muon_pair

DECLARE_SOA_TABLE(MuonPair, "AOD", "MUONPAIR",
                  muon_pair::PairQ,
                  muon_pair::Mass,
                  muon_pair::Pt,
                  muon_pair::Rap);

} // namespace o2::aod
