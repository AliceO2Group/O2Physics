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

#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>

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

namespace matchplaneresidual
{
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(MchCollisionId, mchCollisionId, int64_t);
DECLARE_SOA_COLUMN(MftCollisionId, mftCollisionId, int64_t);
DECLARE_SOA_COLUMN(MatchingPlaneZ, matchingPlaneZ, float);
DECLARE_SOA_COLUMN(BzAtMftCenter, bzAtMftCenter, float);
DECLARE_SOA_COLUMN(MchX, mchX, double);
DECLARE_SOA_COLUMN(MchY, mchY, double);
DECLARE_SOA_COLUMN(MchPhi, mchPhi, double);
DECLARE_SOA_COLUMN(MchTanl, mchTanl, double);
DECLARE_SOA_COLUMN(MchInvQPt, mchInvQPt, double);
DECLARE_SOA_COLUMN(MftX, mftX, double);
DECLARE_SOA_COLUMN(MftY, mftY, double);
DECLARE_SOA_COLUMN(MftPhi, mftPhi, double);
DECLARE_SOA_COLUMN(MftTanl, mftTanl, double);
DECLARE_SOA_COLUMN(MftInvQPt, mftInvQPt, double);
DECLARE_SOA_COLUMN(DeltaX, deltaX, double);
DECLARE_SOA_COLUMN(DeltaY, deltaY, double);
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, double);
DECLARE_SOA_COLUMN(DeltaPhiWrapped, deltaPhiWrapped, double);
DECLARE_SOA_COLUMN(DeltaTanl, deltaTanl, double);
DECLARE_SOA_COLUMN(DeltaInvQPt, deltaInvQPt, double);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, double);
DECLARE_SOA_COLUMN(DeltaPt, deltaPt, double);
DECLARE_SOA_COLUMN(Chi2All, chi2All, double);
DECLARE_SOA_COLUMN(Chi2XYPhiTanl, chi2XYPhiTanl, double);
DECLARE_SOA_COLUMN(Chi2XY, chi2XY, double);
DECLARE_SOA_COLUMN(IsFinite, isFinite, bool);
DECLARE_SOA_COLUMN(MchPairMassPV, mchPairMassPV, double);
DECLARE_SOA_COLUMN(HasMchPairMassPV, hasMchPairMassPV, bool);
DECLARE_SOA_COLUMN(IsMchPairJpsi, isMchPairJpsi, bool);
DECLARE_SOA_COLUMN(MchPairPartnerId, mchPairPartnerId, int64_t);
DECLARE_SOA_COLUMN(NumMchJpsiPartners, numMchJpsiPartners, int);
DECLARE_SOA_COLUMN(AbsorberStepCm, absorberStepCm, double);
// 0=success; 1..5=joint transport/covariance status; 101..106=standalone MFT transport failure.
DECLARE_SOA_COLUMN(TransportStatus, transportStatus, int);
} // namespace matchplaneresidual
DECLARE_SOA_TABLE(DownstreamMatchPlaneResiduals, "AOD", "FWDMPDOWN",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::AbsorberStepCm,
                  matchplaneresidual::TransportStatus);

namespace matchplanemixing
{
DECLARE_SOA_COLUMN(MixingStreamId, mixingStreamId, uint64_t);
DECLARE_SOA_COLUMN(TargetFrameId, targetFrameId, uint64_t);
DECLARE_SOA_COLUMN(DonorFrameId, donorFrameId, uint64_t);
} // namespace matchplanemixing
// Optional diagnostic matrices; row-major 5x5. No change to matching.
namespace matchplanemixdiag
{
DECLARE_SOA_COLUMN(IsMixed, isMixed, bool);
DECLARE_SOA_COLUMN(MchCov, mchCov, std::vector<double>);
DECLARE_SOA_COLUMN(MftCov, mftCov, std::vector<double>);
DECLARE_SOA_COLUMN(ResidualCov, residualCov, std::vector<double>);
} // namespace matchplanemixdiag
DECLARE_SOA_TABLE(MftMchMixingDiagnostics, "AOD", "FWDMPMIXDG",
                  matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  matchplaneresidual::MatchingPlaneZ, matchplanemixdiag::IsMixed,
                  matchplanemixdiag::MchCov, matchplanemixdiag::MftCov, matchplanemixdiag::ResidualCov,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
// Mixed rows use independent tables; both planes carry source IDs and pool diagnostics.
namespace matchplanemixing
{
DECLARE_SOA_COLUMN(DonorMchTrackId, donorMchTrackId, int64_t);
DECLARE_SOA_COLUMN(DonorCollisionId, donorCollisionId, int64_t);
DECLARE_SOA_COLUMN(MixedGroupIndex, mixedGroupIndex, int);
DECLARE_SOA_COLUMN(NumMixedGroups, numMixedGroups, int);
DECLARE_SOA_COLUMN(MixDeltaBc, mixDeltaBc, int64_t);
DECLARE_SOA_COLUMN(MixDeltaPhi, mixDeltaPhi, double);
DECLARE_SOA_COLUMN(MixDeltaR, mixDeltaR, double);
DECLARE_SOA_COLUMN(MixDeltaAttemptsRel, mixDeltaAttemptsRel, double);
DECLARE_SOA_COLUMN(MixDeltaZ, mixDeltaZ, double);
DECLARE_SOA_COLUMN(TargetMatchAttempts, targetMatchAttempts, int);
DECLARE_SOA_COLUMN(DonorMatchAttempts, donorMatchAttempts, int);
DECLARE_SOA_COLUMN(TargetMftMultiplicity, targetMftMultiplicity, int);
DECLARE_SOA_COLUMN(DonorMftMultiplicity, donorMftMultiplicity, int);
DECLARE_SOA_COLUMN(TargetMchMultiplicity, targetMchMultiplicity, int);
DECLARE_SOA_COLUMN(DonorMchMultiplicity, donorMchMultiplicity, int);
DECLARE_SOA_COLUMN(TargetCentrality, targetCentrality, double);
DECLARE_SOA_COLUMN(DonorCentrality, donorCentrality, double);
DECLARE_SOA_COLUMN(TargetCbMultiplicity, targetCbMultiplicity, int);
DECLARE_SOA_COLUMN(DonorCbMultiplicity, donorCbMultiplicity, int);
DECLARE_SOA_COLUMN(MixingWeight, mixingWeight, double);
// Combine frame ID with the output job/stream identity when merging jobs.

DECLARE_SOA_COLUMN(RetainedFrames, retainedFrames, std::vector<uint64_t>);
DECLARE_SOA_COLUMN(RetainedEvents, retainedEvents, uint64_t);
DECLARE_SOA_COLUMN(RetainedMftStates, retainedMftStates, uint64_t);
DECLARE_SOA_COLUMN(EvictedEvents, evictedEvents, uint64_t);
DECLARE_SOA_COLUMN(OversizedEvents, oversizedEvents, uint64_t);
DECLARE_SOA_COLUMN(PoolWasReset, poolWasReset, bool);
} // namespace matchplanemixing
// One state row per input DF, including empty/rejected DFs. Also synchronizes
// MC label-cache lifetime with the exact donor-history frames owned upstream.
DECLARE_SOA_TABLE(MftMchMixingPoolStates, "AOD", "FWDMPMIXPS",
                  matchplanemixing::TargetFrameId, matchplanemixing::RetainedFrames,
                  matchplanemixing::RetainedEvents, matchplanemixing::RetainedMftStates,
                  matchplanemixing::EvictedEvents, matchplanemixing::OversizedEvents,
                  matchplanemixing::PoolWasReset, matchplanemixing::MixingStreamId);
DECLARE_SOA_TABLE(RawMixedMatchPlaneResiduals, "AOD", "FWDMPMIXRAW",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::AbsorberStepCm,
                  matchplaneresidual::TransportStatus,
                  matchplanemixing::DonorMchTrackId,
                  matchplanemixing::DonorCollisionId,
                  matchplanemixing::MixedGroupIndex,
                  matchplanemixing::NumMixedGroups,
                  matchplanemixing::MixDeltaBc,
                  matchplanemixing::MixDeltaPhi,
                  matchplanemixing::MixDeltaR,
                  matchplanemixing::MixDeltaAttemptsRel,
                  matchplanemixing::MixDeltaZ,
                  matchplanemixing::TargetMatchAttempts,
                  matchplanemixing::DonorMatchAttempts,
                  matchplanemixing::TargetMftMultiplicity,
                  matchplanemixing::DonorMftMultiplicity,
                  matchplanemixing::TargetMchMultiplicity,
                  matchplanemixing::DonorMchMultiplicity,
                  matchplanemixing::TargetCentrality,
                  matchplanemixing::DonorCentrality,
                  matchplanemixing::TargetCbMultiplicity,
                  matchplanemixing::DonorCbMultiplicity,
                  matchplanemixing::MixingWeight,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
DECLARE_SOA_TABLE(MixedMatchPlaneResiduals, "AOD", "FWDMPMIX",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::AbsorberStepCm,
                  matchplaneresidual::TransportStatus,
                  matchplanemixing::DonorMchTrackId,
                  matchplanemixing::DonorCollisionId,
                  matchplanemixing::MixedGroupIndex,
                  matchplanemixing::NumMixedGroups,
                  matchplanemixing::MixDeltaBc,
                  matchplanemixing::MixDeltaPhi,
                  matchplanemixing::MixDeltaR,
                  matchplanemixing::MixDeltaAttemptsRel,
                  matchplanemixing::MixDeltaZ,
                  matchplanemixing::TargetMatchAttempts,
                  matchplanemixing::DonorMatchAttempts,
                  matchplanemixing::TargetMftMultiplicity,
                  matchplanemixing::DonorMftMultiplicity,
                  matchplanemixing::TargetMchMultiplicity,
                  matchplanemixing::DonorMchMultiplicity,
                  matchplanemixing::TargetCentrality,
                  matchplanemixing::DonorCentrality,
                  matchplanemixing::TargetCbMultiplicity,
                  matchplanemixing::DonorCbMultiplicity,
                  matchplanemixing::MixingWeight,
                  matchplanemixing::TargetFrameId, matchplanemixing::DonorFrameId, matchplanemixing::MixingStreamId);
DECLARE_SOA_TABLE(RawMatchPlaneResiduals, "AOD", "FWDMPRERAW",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite);
DECLARE_SOA_TABLE(MatchPlaneResiduals, "AOD", "FWDMPRESID",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::MchPairMassPV,
                  matchplaneresidual::HasMchPairMassPV,
                  matchplaneresidual::IsMchPairJpsi,
                  matchplaneresidual::MchPairPartnerId,
                  matchplaneresidual::NumMchJpsiPartners);
DECLARE_SOA_TABLE(JpsiMatchPlaneResiduals, "AOD", "FWDMPJPSI",
                  matchplaneresidual::MchTrackId,
                  matchplaneresidual::MftTrackId,
                  matchplaneresidual::MchCollisionId,
                  matchplaneresidual::MftCollisionId,
                  matchplaneresidual::MatchingPlaneZ,
                  matchplaneresidual::BzAtMftCenter,
                  matchplaneresidual::MchX,
                  matchplaneresidual::MchY,
                  matchplaneresidual::MchPhi,
                  matchplaneresidual::MchTanl,
                  matchplaneresidual::MchInvQPt,
                  matchplaneresidual::MftX,
                  matchplaneresidual::MftY,
                  matchplaneresidual::MftPhi,
                  matchplaneresidual::MftTanl,
                  matchplaneresidual::MftInvQPt,
                  matchplaneresidual::DeltaX,
                  matchplaneresidual::DeltaY,
                  matchplaneresidual::DeltaPhi,
                  matchplaneresidual::DeltaPhiWrapped,
                  matchplaneresidual::DeltaTanl,
                  matchplaneresidual::DeltaInvQPt,
                  matchplaneresidual::DeltaEta,
                  matchplaneresidual::DeltaPt,
                  matchplaneresidual::Chi2All,
                  matchplaneresidual::Chi2XYPhiTanl,
                  matchplaneresidual::Chi2XY,
                  matchplaneresidual::IsFinite,
                  matchplaneresidual::MchPairMassPV,
                  matchplaneresidual::HasMchPairMassPV,
                  matchplaneresidual::IsMchPairJpsi,
                  matchplaneresidual::MchPairPartnerId,
                  matchplaneresidual::NumMchJpsiPartners);
// One row per input MCH(-MID) track, including rejected tracks for selection QA.
// IDs must be combined with the input file and DF identifier.
namespace matchplaneeventsel
{
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MchCollisionId, mchCollisionId, int64_t);
DECLARE_SOA_COLUMN(SelectionMode, selectionMode, int8_t); // 0=all, 1=int8/sel8, 2=UPC candidate
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);
DECLARE_SOA_COLUMN(HasMchBC, hasMchBC, bool);
DECLARE_SOA_COLUMN(MchGlobalBC, mchGlobalBC, uint64_t);
DECLARE_SOA_COLUMN(NumMchMidTracks, numMchMidTracks, int);
DECLARE_SOA_COLUMN(MchSel8, mchSel8, bool);
DECLARE_SOA_COLUMN(HasFV0, hasFV0, bool);
DECLARE_SOA_COLUMN(MaxFV0Amplitude, maxFV0Amplitude, float);
DECLARE_SOA_COLUMN(NumMchTracks, numMchTracks, int); // MCH-MID + MCH-only at original collisionId
DECLARE_SOA_COLUMN(NumMftTracks, numMftTracks, int); // All MFT at same collisionId, before covariance selection
DECLARE_SOA_COLUMN(PassTrackMultiplicity, passTrackMultiplicity, bool);
DECLARE_SOA_COLUMN(PassMchUnlikeSign, passMchUnlikeSign, bool);
DECLARE_SOA_COLUMN(PassFIT, passFIT, bool); // UDHelpers cleanFIT in configured BC window; true outside UPC mode
DECLARE_SOA_COLUMN(PassMuonQuality, passMuonQuality, bool);
DECLARE_SOA_COLUMN(IsMuonPropagated, isMuonPropagated, bool);
DECLARE_SOA_COLUMN(MuonCutEta, muonCutEta, float);
DECLARE_SOA_COLUMN(MuonCutPt, muonCutPt, float);
DECLARE_SOA_COLUMN(MftAbsorberStepCm, mftAbsorberStepCm, double); // effective material delta-z; 0=whole absorber/inactive
} // namespace matchplaneeventsel
DECLARE_SOA_TABLE(MatchPlaneEventSels, "AOD", "FWDMPEVSEL",
                  matchplaneeventsel::MchTrackId,
                  matchplaneeventsel::MchCollisionId,
                  matchplaneeventsel::SelectionMode,
                  matchplaneeventsel::IsSelected,
                  matchplaneeventsel::HasMchBC,
                  matchplaneeventsel::MchGlobalBC,
                  matchplaneeventsel::NumMchMidTracks,
                  matchplaneeventsel::MchSel8,
                  matchplaneeventsel::HasFV0,
                  matchplaneeventsel::MaxFV0Amplitude,
                  matchplaneeventsel::NumMchTracks,
                  matchplaneeventsel::NumMftTracks,
                  matchplaneeventsel::PassTrackMultiplicity,
                  matchplaneeventsel::PassMchUnlikeSign,
                  matchplaneeventsel::PassFIT,
                  matchplaneeventsel::PassMuonQuality,
                  matchplaneeventsel::IsMuonPropagated,
                  matchplaneeventsel::MuonCutEta,
                  matchplaneeventsel::MuonCutPt,
                  matchplaneeventsel::MftAbsorberStepCm);
// One row per selected same-collision unordered MCH pair, independent of MFT candidates.
namespace mchdimuonpv
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId1, mchTrackId1, int64_t);
DECLARE_SOA_COLUMN(MchTrackId2, mchTrackId2, int64_t);
DECLARE_SOA_COLUMN(PairCharge, pairCharge, int8_t);
DECLARE_SOA_COLUMN(VertexX, vertexX, float);
DECLARE_SOA_COLUMN(VertexY, vertexY, float);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
DECLARE_SOA_COLUMN(Px1, px1, double);
DECLARE_SOA_COLUMN(Py1, py1, double);
DECLARE_SOA_COLUMN(Pz1, pz1, double);
DECLARE_SOA_COLUMN(Px2, px2, double);
DECLARE_SOA_COLUMN(Py2, py2, double);
DECLARE_SOA_COLUMN(Pz2, pz2, double);
DECLARE_SOA_COLUMN(Mass, mass, double);
DECLARE_SOA_COLUMN(PairPt, pairPt, double);
DECLARE_SOA_COLUMN(PairRapidity, pairRapidity, double);
} // namespace mchdimuonpv
DECLARE_SOA_TABLE(MchDimuonsPV, "AOD", "FWDMCHPAIR",
                  mchdimuonpv::CollisionId, mchdimuonpv::MchTrackId1, mchdimuonpv::MchTrackId2,
                  mchdimuonpv::PairCharge, mchdimuonpv::VertexX, mchdimuonpv::VertexY, mchdimuonpv::VertexZ,
                  mchdimuonpv::Px1, mchdimuonpv::Py1, mchdimuonpv::Pz1,
                  mchdimuonpv::Px2, mchdimuonpv::Py2, mchdimuonpv::Pz2,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity);
// ML selections are reconstructed from real-data observables only.
namespace mlmatchedpv
{
DECLARE_SOA_COLUMN(FrameId, frameId, uint64_t);
DECLARE_SOA_COLUMN(MlGoodScore, mlGoodScore, double);
DECLARE_SOA_COLUMN(HasMftCovariance, hasMftCovariance, bool);
DECLARE_SOA_COLUMN(Status, status, int); // 0=success, 1=MCH PV failure, 2=invalid momentum, 3=nonfinite matched PV state
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(Px, px, double);
DECLARE_SOA_COLUMN(Py, py, double);
DECLARE_SOA_COLUMN(Pz, pz, double);
DECLARE_SOA_COLUMN(MftTrackId1, mftTrackId1, int64_t);
DECLARE_SOA_COLUMN(MftTrackId2, mftTrackId2, int64_t);
DECLARE_SOA_COLUMN(MlGoodScore1, mlGoodScore1, double);
DECLARE_SOA_COLUMN(MlGoodScore2, mlGoodScore2, double);
} // namespace mlmatchedpv
DECLARE_SOA_TABLE(MlMatchedMuonsPV, "AOD", "FWDMLMUONPV",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  mlmatchedpv::MlGoodScore, mlmatchedpv::HasMftCovariance, mlmatchedpv::Status,
                  mlmatchedpv::Charge, mlmatchedpv::Px, mlmatchedpv::Py, mlmatchedpv::Pz);
DECLARE_SOA_TABLE(MlMatchedDimuonsPV, "AOD", "FWDMLPAIRPV",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId,
                  mchdimuonpv::MchTrackId1, mlmatchedpv::MftTrackId1, mchdimuonpv::MchTrackId2, mlmatchedpv::MftTrackId2,
                  mlmatchedpv::MlGoodScore1, mlmatchedpv::MlGoodScore2, mchdimuonpv::PairCharge,
                  mchdimuonpv::Px1, mchdimuonpv::Py1, mchdimuonpv::Pz1, mchdimuonpv::Px2, mchdimuonpv::Py2, mchdimuonpv::Pz2,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity);

// Pair-specific covariance diagnostics. Matrices are row-major. Joint order:
// [MFT(x,y,phi,tanl,q/pt), MCH(x,y,phi,tanl,q/pt)] at matching plane.
namespace mftmchpaircov
{
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(PairCovValid, pairCovValid, bool);
DECLARE_SOA_COLUMN(PairCovStatus, pairCovStatus, int);
DECLARE_SOA_COLUMN(SeedMomentum, seedMomentum, double);
DECLARE_SOA_COLUMN(OriginalMftInvQPt, originalMftInvQPt, double);
DECLARE_SOA_COLUMN(UsesMchCharge, usesMchCharge, bool);
DECLARE_SOA_COLUMN(RankAll, rankAll, int);
DECLARE_SOA_COLUMN(RankFour, rankFour, int);
DECLARE_SOA_COLUMN(RankXY, rankXY, int);
DECLARE_SOA_COLUMN(JointCov, jointCov, std::vector<double>);
DECLARE_SOA_COLUMN(ResidualCov, residualCov, std::vector<double>);
DECLARE_SOA_COLUMN(SeedCov, seedCov, std::vector<double>);
DECLARE_SOA_COLUMN(SeedMchCrossCov, seedMchCrossCov, std::vector<double>);
DECLARE_SOA_COLUMN(PairJacobian, pairJacobian, std::vector<double>);
DECLARE_SOA_COLUMN(UpstreamSigmaELoss, upstreamSigmaELoss, double);
DECLARE_SOA_COLUMN(JacobianScale, jacobianScale, double);
} // namespace mftmchpaircov
DECLARE_SOA_TABLE(MftMchPairCovariances, "AOD", "FWDMPPAIRCV",
                  mftmchpaircov::MchTrackId,
                  mftmchpaircov::MftTrackId,
                  mftmchpaircov::PairCovValid,
                  mftmchpaircov::PairCovStatus,
                  mftmchpaircov::SeedMomentum,
                  mftmchpaircov::OriginalMftInvQPt,
                  mftmchpaircov::UsesMchCharge,
                  mftmchpaircov::RankAll,
                  mftmchpaircov::RankFour,
                  mftmchpaircov::RankXY,
                  mftmchpaircov::JointCov,
                  mftmchpaircov::ResidualCov,
                  mftmchpaircov::SeedCov,
                  mftmchpaircov::SeedMchCrossCov,
                  mftmchpaircov::PairJacobian,
                  mftmchpaircov::UpstreamSigmaELoss,
                  mftmchpaircov::JacobianScale);
namespace selectedmatchingmuon
{
DECLARE_SOA_COLUMN(Method, method, int8_t);     // 0=MCHPV,1=chi2,2=ML
DECLARE_SOA_COLUMN(Score, score, double);       // chi2 or good-source score; MCHPV NaN
DECLARE_SOA_COLUMN(RawScore, rawScore, double); // chi2 or good logit; MCHPV NaN
DECLARE_SOA_COLUMN(DcaX, dcaX, double);
DECLARE_SOA_COLUMN(DcaY, dcaY, double);
DECLARE_SOA_COLUMN(DcaZ, dcaZ, double);
} // namespace selectedmatchingmuon
DECLARE_SOA_TABLE(SelectedMatchingMuons, "AOD", "FWDEVALMUON",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  selectedmatchingmuon::Method, selectedmatchingmuon::Score, selectedmatchingmuon::RawScore,
                  mlmatchedpv::Charge, mlmatchedpv::Px, mlmatchedpv::Py, mlmatchedpv::Pz,
                  selectedmatchingmuon::DcaX, selectedmatchingmuon::DcaY, selectedmatchingmuon::DcaZ);
// Internal pre-cut one-to-one assignments; never persisted as an AOD training table.
DECLARE_SOA_TABLE(AssignedMatchingMuons, "AOD2", "FWDASSIGNMUON",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId,
                  selectedmatchingmuon::Method, selectedmatchingmuon::Score, selectedmatchingmuon::RawScore,
                  mlmatchedpv::Charge, mlmatchedpv::Px, mlmatchedpv::Py, mlmatchedpv::Pz,
                  selectedmatchingmuon::DcaX, selectedmatchingmuon::DcaY, selectedmatchingmuon::DcaZ);

// Internal accepted MCH-only pair transport to the data/MC sparse writer.
// Same IDs are local to the current DF; Mixed is not annotated with MC truth.
namespace mchdimuontransport
{
DECLARE_SOA_COLUMN(MchId1, mchId1, int64_t);
DECLARE_SOA_COLUMN(MchId2, mchId2, int64_t);
DECLARE_SOA_COLUMN(IsMixed, isMixed, bool);
DECLARE_SOA_COLUMN(Weight, weight, double);
DECLARE_SOA_COLUMN(Phi, phi, double);
} // namespace mchdimuontransport
DECLARE_SOA_TABLE(RawMchDimuonPairs, "AOD2", "FWDMCHDIMURAW",
                  mchdimuontransport::MchId1, mchdimuontransport::MchId2,
                  mchdimuontransport::IsMixed, mchdimuontransport::Weight,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity,
                  mchdimuontransport::Phi, mchdimuonpv::PairCharge);
// Same-event PCA audit: emitted only when dimuon vertexing is enabled.
// Failed fits are retained; no matching or single-muon denominator is changed.
namespace matchingdimuonvertex
{
DECLARE_SOA_COLUMN(MchTrackId1, mchTrackId1, int64_t);
DECLARE_SOA_COLUMN(MchTrackId2, mchTrackId2, int64_t);
DECLARE_SOA_COLUMN(FitSuccess, fitSuccess, bool);
DECLARE_SOA_COLUMN(UsesCovariance, usesCovariance, bool);
DECLARE_SOA_COLUMN(PassesVertexCuts, passesVertexCuts, bool);
DECLARE_SOA_COLUMN(PcaX, pcaX, double);
DECLARE_SOA_COLUMN(PcaY, pcaY, double);
DECLARE_SOA_COLUMN(PcaZ, pcaZ, double);
DECLARE_SOA_COLUMN(PcaPVXY, pcaPVXY, double);
DECLARE_SOA_COLUMN(PcaPVZ, pcaPVZ, double);
DECLARE_SOA_COLUMN(PairDca, pairDca, double);
DECLARE_SOA_COLUMN(VertexObjective, vertexObjective, double);
DECLARE_SOA_COLUMN(Lxy, lxy, double);
} // namespace matchingdimuonvertex
DECLARE_SOA_TABLE(MatchingDimuonVertices, "AOD", "FWDDIMUVTX",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, selectedmatchingmuon::Method,
                  matchingdimuonvertex::MchTrackId1, matchingdimuonvertex::MchTrackId2,
                  matchingdimuonvertex::FitSuccess, matchingdimuonvertex::UsesCovariance, matchingdimuonvertex::PassesVertexCuts,
                  matchingdimuonvertex::PcaX, matchingdimuonvertex::PcaY, matchingdimuonvertex::PcaZ,
                  matchingdimuonvertex::PcaPVXY, matchingdimuonvertex::PcaPVZ, matchingdimuonvertex::PairDca,
                  matchingdimuonvertex::VertexObjective, matchingdimuonvertex::Lxy);

// Pre-score/DCA/vertex-cut baseline GlobalMuon dimuons. Same MC IDs are DF-local; Mixed has
// no truth annotation. This AOD2 transport is never a persisted training table.
namespace globaldimuontransport
{
DECLARE_SOA_COLUMN(MftId1, mftId1, int64_t);
DECLARE_SOA_COLUMN(MftId2, mftId2, int64_t);
DECLARE_SOA_COLUMN(WorkingPoint, workingPoint, std::vector<uint8_t>); // UTF-8 name; ASoA does not support string columns.
DECLARE_SOA_COLUMN(DimuonDcaXY, dimuonDcaXY, double);
DECLARE_SOA_COLUMN(PairMatchingScore, pairMatchingScore, double); // ML min score; chi2 max score of the two assigned muons.
} // namespace globaldimuontransport
DECLARE_SOA_TABLE(RawGlobalDimuonPairs, "AOD2", "FWDGLOBDIMU",
                  mchdimuontransport::MchId1, mchdimuontransport::MchId2,
                  globaldimuontransport::MftId1, globaldimuontransport::MftId2,
                  selectedmatchingmuon::Method, mchdimuontransport::IsMixed,
                  mchdimuonpv::PairCharge, mchdimuontransport::Weight, globaldimuontransport::WorkingPoint,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity,
                  matchingdimuonvertex::PairDca, matchingdimuonvertex::PcaPVXY, matchingdimuonvertex::PcaPVZ,
                  globaldimuontransport::DimuonDcaXY, matchingdimuonvertex::Lxy, matchingdimuonvertex::FitSuccess,
                  globaldimuontransport::PairMatchingScore);

// Audit tables are independent of propagation success, assignment and score cuts.
// FrameId + collision/MCH/MFT IDs are the persistent keys; never join across DFs
// using track IDs alone. Original MFT collision IDs are NOT overwritten.
namespace matchingaudit
{
DECLARE_SOA_COLUMN(IsEligible, isEligible, bool);
DECLARE_SOA_COLUMN(HasPVState, hasPVState, bool);
DECLARE_SOA_COLUMN(MchPtPV, mchPtPV, double);
DECLARE_SOA_COLUMN(MchEtaPV, mchEtaPV, double);
DECLARE_SOA_COLUMN(MchPhiPV, mchPhiPV, double);
DECLARE_SOA_COLUMN(MchMomentumPV, mchMomentumPV, double);
DECLARE_SOA_COLUMN(NumPreselectedMfts, numPreselectedMfts, int);
DECLARE_SOA_COLUMN(IsMchMid, isMchMid, bool);
DECLARE_SOA_COLUMN(NumCandidates, numCandidates, int);
DECLARE_SOA_COLUMN(NumFiniteUpstream, numFiniteUpstream, int);
DECLARE_SOA_COLUMN(NumFiniteDownstream, numFiniteDownstream, int);
DECLARE_SOA_COLUMN(UsesAssociation, usesAssociation, bool);
DECLARE_SOA_COLUMN(IsNativeCandidate, isNativeCandidate, bool);
DECLARE_SOA_COLUMN(IsAssociatedCandidate, isAssociatedCandidate, bool);
DECLARE_SOA_COLUMN(HasMftCovariance, hasMftCovariance, bool);
DECLARE_SOA_COLUMN(UpstreamFinite, upstreamFinite, bool);
DECLARE_SOA_COLUMN(DownstreamFinite, downstreamFinite, bool);
DECLARE_SOA_COLUMN(DownstreamStatus, downstreamStatus, int);
// 0=no candidates; 1=PV failure; 2=no finite score; 3=one-to-one conflict;
// 4=score rejected; 5=matched DCA failure; 6=DCA rejected; 7=accepted;
// 8=MID requirement failed. Disabled methods have no outcome row.
DECLARE_SOA_COLUMN(Outcome, outcome, int8_t);
} // namespace matchingaudit
DECLARE_SOA_TABLE(MatchingMchPopulation, "AOD", "FWDMCHPOP",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId, matchplaneresidual::MchTrackId,
                  mlmatchedpv::Charge, matchingaudit::IsEligible, matchingaudit::HasPVState,
                  matchingaudit::MchPtPV, matchingaudit::MchEtaPV, matchingaudit::MchPhiPV, matchingaudit::MchMomentumPV,
                  matchingaudit::NumCandidates, matchingaudit::NumFiniteUpstream, matchingaudit::NumFiniteDownstream,
                  matchingaudit::UsesAssociation, matchingaudit::NumPreselectedMfts, matchingaudit::IsMchMid);
// Exactly one row per eligible same-event pair BEFORE any propagation/score cuts.
// Failed transports remain here and in the aligned downstream residual table.
DECLARE_SOA_TABLE(MatchingCandidates, "AOD", "FWDMPCAND",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId,
                  matchplaneresidual::MchTrackId, matchplaneresidual::MftTrackId, matchplaneresidual::MftCollisionId,
                  matchingaudit::IsNativeCandidate, matchingaudit::IsAssociatedCandidate,
                  matchingaudit::HasMftCovariance, matchingaudit::UpstreamFinite,
                  matchingaudit::DownstreamFinite, matchingaudit::DownstreamStatus);
// One row per eligible MCH and enabled matching method, including failures.
// Assigned ID/score are saved BEFORE score/DCA cuts, allowing offline scans.
DECLARE_SOA_TABLE(MatchingMchOutcomes, "AOD", "FWDMCHOUT",
                  mlmatchedpv::FrameId, mchdimuonpv::CollisionId,
                  matchplaneresidual::MchTrackId, selectedmatchingmuon::Method,
                  matchingaudit::Outcome, matchplaneresidual::MftTrackId,
                  selectedmatchingmuon::Score, selectedmatchingmuon::RawScore);

// MC annotations share the matching data model; reconstruction does not use these labels.
namespace mftmchmatchtruth
{
DECLARE_SOA_COLUMN(ResidualRowId, residualRowId, int64_t);
DECLARE_SOA_COLUMN(MchTrackId, mchTrackId, int64_t);
DECLARE_SOA_COLUMN(MftTrackId, mftTrackId, int64_t);
DECLARE_SOA_COLUMN(MchMcParticleId, mchMcParticleId, int64_t);
DECLARE_SOA_COLUMN(MftMcParticleId, mftMcParticleId, int64_t);
DECLARE_SOA_COLUMN(HasMchMcParticle, hasMchMcParticle, bool);
DECLARE_SOA_COLUMN(HasMftMcParticle, hasMftMcParticle, bool);
// Strict correct: both labels valid, same MC particle, neither fake bit set.
DECLARE_SOA_COLUMN(IsCorrectMatch, isCorrectMatch, bool);
// False includes either fake label, even for equal IDs. Unknown is separate.
DECLARE_SOA_COLUMN(IsFalseMatch, isFalseMatch, bool);
// -1 = missing/invalid non-fake label; 0 = false (including fake); 1 = strict correct.
DECLARE_SOA_COLUMN(MatchTruthStatus, matchTruthStatus, int8_t);
DECLARE_SOA_COLUMN(MchMcMask, mchMcMask, uint8_t);
DECLARE_SOA_COLUMN(MftMcMask, mftMcMask, uint8_t);
DECLARE_SOA_COLUMN(MchMcCollisionId, mchMcCollisionId, int64_t);
DECLARE_SOA_COLUMN(MftMcCollisionId, mftMcCollisionId, int64_t);
DECLARE_SOA_COLUMN(MchPdgCode, mchPdgCode, int);
DECLARE_SOA_COLUMN(MftPdgCode, mftPdgCode, int);
DECLARE_SOA_COLUMN(MchMcPx, mchMcPx, float);
DECLARE_SOA_COLUMN(MchMcPy, mchMcPy, float);
DECLARE_SOA_COLUMN(MchMcPz, mchMcPz, float);
DECLARE_SOA_COLUMN(MftMcPx, mftMcPx, float);
DECLARE_SOA_COLUMN(MftMcPy, mftMcPy, float);
DECLARE_SOA_COLUMN(MftMcPz, mftMcPz, float);
DECLARE_SOA_COLUMN(MchPhysicalPrimary, mchPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(MftPhysicalPrimary, mftPhysicalPrimary, bool);
} // namespace mftmchmatchtruth
DECLARE_SOA_TABLE(MftMchMatchTruths, "AOD", "FWDMPMCLBL",
                  mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
                  mftmchmatchtruth::MchMcParticleId, mftmchmatchtruth::MftMcParticleId,
                  mftmchmatchtruth::HasMchMcParticle, mftmchmatchtruth::HasMftMcParticle, mftmchmatchtruth::IsCorrectMatch,
                  mftmchmatchtruth::MchMcMask, mftmchmatchtruth::MftMcMask,
                  mftmchmatchtruth::MchMcCollisionId, mftmchmatchtruth::MftMcCollisionId,
                  mftmchmatchtruth::MchPdgCode, mftmchmatchtruth::MftPdgCode,
                  mftmchmatchtruth::MchMcPx, mftmchmatchtruth::MchMcPy, mftmchmatchtruth::MchMcPz,
                  mftmchmatchtruth::MftMcPx, mftmchmatchtruth::MftMcPy, mftmchmatchtruth::MftMcPz,
                  mftmchmatchtruth::MchPhysicalPrimary, mftmchmatchtruth::MftPhysicalPrimary,
                  mftmchmatchtruth::IsFalseMatch, mftmchmatchtruth::MatchTruthStatus);
DECLARE_SOA_TABLE(MixedMftMchMatchTruths, "AOD", "FWDMPMIXLBL",
                  mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
                  mftmchmatchtruth::MchMcParticleId, mftmchmatchtruth::MftMcParticleId,
                  mftmchmatchtruth::HasMchMcParticle, mftmchmatchtruth::HasMftMcParticle, mftmchmatchtruth::IsCorrectMatch,
                  mftmchmatchtruth::MchMcMask, mftmchmatchtruth::MftMcMask,
                  mftmchmatchtruth::MchMcCollisionId, mftmchmatchtruth::MftMcCollisionId,
                  mftmchmatchtruth::MchPdgCode, mftmchmatchtruth::MftPdgCode,
                  mftmchmatchtruth::MchMcPx, mftmchmatchtruth::MchMcPy, mftmchmatchtruth::MchMcPz,
                  mftmchmatchtruth::MftMcPx, mftmchmatchtruth::MftMcPy, mftmchmatchtruth::MftMcPz,
                  mftmchmatchtruth::MchPhysicalPrimary, mftmchmatchtruth::MftPhysicalPrimary,
                  mftmchmatchtruth::IsFalseMatch, mftmchmatchtruth::MatchTruthStatus);
namespace mftmchmatchrelation
{
// -1=unknown, 0=other, 1=strict same particle, 2=ancestor decay,
// 3=ancestor non-decay, 4=either track fake; detailed genealogy is retained.
DECLARE_SOA_COLUMN(PairOrigin, pairOrigin, int8_t);
DECLARE_SOA_COLUMN(AncestryStatus, ancestryStatus, int);
DECLARE_SOA_COLUMN(AncestryComplete, ancestryComplete, bool);
DECLARE_SOA_COLUMN(IsMchMuon, isMchMuon, bool);
DECLARE_SOA_COLUMN(IsMftAncestorOfMch, isMftAncestorOfMch, bool);
DECLARE_SOA_COLUMN(IsMftToMuonDecay, isMftToMuonDecay, bool);
DECLARE_SOA_COLUMN(IsDirectMftToMuonDecay, isDirectMftToMuonDecay, bool);
DECLARE_SOA_COLUMN(AncestryDepth, ancestryDepth, int);
DECLARE_SOA_COLUMN(MchMcProcess, mchMcProcess, int);
DECLARE_SOA_COLUMN(MftMcProcess, mftMcProcess, int);
DECLARE_SOA_COLUMN(MchProducedByGenerator, mchProducedByGenerator, bool);
DECLARE_SOA_COLUMN(MftProducedByGenerator, mftProducedByGenerator, bool);
DECLARE_SOA_COLUMN(MchMcVertexX, mchMcVertexX, float);
DECLARE_SOA_COLUMN(MchMcVertexY, mchMcVertexY, float);
DECLARE_SOA_COLUMN(MchMcVertexZ, mchMcVertexZ, float);
DECLARE_SOA_COLUMN(MftMcVertexX, mftMcVertexX, float);
DECLARE_SOA_COLUMN(MftMcVertexY, mftMcVertexY, float);
DECLARE_SOA_COLUMN(MftMcVertexZ, mftMcVertexZ, float);
DECLARE_SOA_COLUMN(MchMotherIds, mchMotherIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MchMotherPdgs, mchMotherPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(MftMotherIds, mftMotherIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MftMotherPdgs, mftMotherPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(AncestryIds, ancestryIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(AncestryPdgs, ancestryPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(AncestryProcesses, ancestryProcesses, std::vector<int>);
DECLARE_SOA_COLUMN(AncestryGeneratorFlags, ancestryGeneratorFlags, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(AncestryVertexX, ancestryVertexX, std::vector<float>);
DECLARE_SOA_COLUMN(AncestryVertexY, ancestryVertexY, std::vector<float>);
DECLARE_SOA_COLUMN(AncestryVertexZ, ancestryVertexZ, std::vector<float>);
} // namespace mftmchmatchrelation
// Row-aligned with MatchPlaneResiduals, DownstreamMatchPlaneResiduals and MftMchMatchTruths.
DECLARE_SOA_TABLE(MftMchMatchRelations, "AOD", "FWDMPMCREL",
                  mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
                  mftmchmatchrelation::PairOrigin,
                  mftmchmatchrelation::AncestryStatus,
                  mftmchmatchrelation::AncestryComplete,
                  mftmchmatchrelation::IsMchMuon,
                  mftmchmatchrelation::IsMftAncestorOfMch,
                  mftmchmatchrelation::IsMftToMuonDecay,
                  mftmchmatchrelation::IsDirectMftToMuonDecay,
                  mftmchmatchrelation::AncestryDepth,
                  mftmchmatchrelation::MchMcProcess,
                  mftmchmatchrelation::MftMcProcess,
                  mftmchmatchrelation::MchProducedByGenerator,
                  mftmchmatchrelation::MftProducedByGenerator,
                  mftmchmatchrelation::MchMcVertexX,
                  mftmchmatchrelation::MchMcVertexY,
                  mftmchmatchrelation::MchMcVertexZ,
                  mftmchmatchrelation::MftMcVertexX,
                  mftmchmatchrelation::MftMcVertexY,
                  mftmchmatchrelation::MftMcVertexZ,
                  mftmchmatchrelation::MchMotherIds,
                  mftmchmatchrelation::MchMotherPdgs,
                  mftmchmatchrelation::MftMotherIds,
                  mftmchmatchrelation::MftMotherPdgs,
                  mftmchmatchrelation::AncestryIds,
                  mftmchmatchrelation::AncestryPdgs,
                  mftmchmatchrelation::AncestryProcesses,
                  mftmchmatchrelation::AncestryGeneratorFlags,
                  mftmchmatchrelation::AncestryVertexX,
                  mftmchmatchrelation::AncestryVertexY,
                  mftmchmatchrelation::AncestryVertexZ);
// Row-aligned to MatchingCandidates, independent of propagation and assignment.
DECLARE_SOA_TABLE(MatchingCandidateTruths, "AOD", "FWDCANDTRUTH",
                  mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
                  mftmchmatchtruth::MchMcParticleId, mftmchmatchtruth::MftMcParticleId,
                  mftmchmatchtruth::HasMchMcParticle, mftmchmatchtruth::HasMftMcParticle, mftmchmatchtruth::IsCorrectMatch,
                  mftmchmatchtruth::MchMcMask, mftmchmatchtruth::MftMcMask,
                  mftmchmatchtruth::MchMcCollisionId, mftmchmatchtruth::MftMcCollisionId,
                  mftmchmatchtruth::MchPdgCode, mftmchmatchtruth::MftPdgCode,
                  mftmchmatchtruth::MchMcPx, mftmchmatchtruth::MchMcPy, mftmchmatchtruth::MchMcPz,
                  mftmchmatchtruth::MftMcPx, mftmchmatchtruth::MftMcPy, mftmchmatchtruth::MftMcPz,
                  mftmchmatchtruth::MchPhysicalPrimary, mftmchmatchtruth::MftPhysicalPrimary,
                  mftmchmatchtruth::IsFalseMatch, mftmchmatchtruth::MatchTruthStatus);
DECLARE_SOA_TABLE(MatchingCandidateRelations, "AOD", "FWDCANDREL",
                  mftmchmatchtruth::ResidualRowId, mftmchmatchtruth::MchTrackId, mftmchmatchtruth::MftTrackId,
                  mftmchmatchrelation::PairOrigin,
                  mftmchmatchrelation::AncestryStatus,
                  mftmchmatchrelation::AncestryComplete,
                  mftmchmatchrelation::IsMchMuon,
                  mftmchmatchrelation::IsMftAncestorOfMch,
                  mftmchmatchrelation::IsMftToMuonDecay,
                  mftmchmatchrelation::IsDirectMftToMuonDecay,
                  mftmchmatchrelation::AncestryDepth,
                  mftmchmatchrelation::MchMcProcess,
                  mftmchmatchrelation::MftMcProcess,
                  mftmchmatchrelation::MchProducedByGenerator,
                  mftmchmatchrelation::MftProducedByGenerator,
                  mftmchmatchrelation::MchMcVertexX,
                  mftmchmatchrelation::MchMcVertexY,
                  mftmchmatchrelation::MchMcVertexZ,
                  mftmchmatchrelation::MftMcVertexX,
                  mftmchmatchrelation::MftMcVertexY,
                  mftmchmatchrelation::MftMcVertexZ,
                  mftmchmatchrelation::MchMotherIds,
                  mftmchmatchrelation::MchMotherPdgs,
                  mftmchmatchrelation::MftMotherIds,
                  mftmchmatchrelation::MftMotherPdgs,
                  mftmchmatchrelation::AncestryIds,
                  mftmchmatchrelation::AncestryPdgs,
                  mftmchmatchrelation::AncestryProcesses,
                  mftmchmatchrelation::AncestryGeneratorFlags,
                  mftmchmatchrelation::AncestryVertexX,
                  mftmchmatchrelation::AncestryVertexY,
                  mftmchmatchrelation::AncestryVertexZ);
namespace matchingmctruth
{
DECLARE_SOA_COLUMN(PopulationRowId, populationRowId, int64_t);
DECLARE_SOA_COLUMN(NumTruePairs, numTruePairs, int);
DECLARE_SOA_COLUMN(IsDecayMuonMch, isDecayMuonMch, bool);
DECLARE_SOA_COLUMN(AncestryStatus, ancestryStatus, int);
DECLARE_SOA_COLUMN(ChargedAncestorIds, chargedAncestorIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(ChargedAncestorPdgs, chargedAncestorPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(SignalOrigin, signalOrigin, int8_t);
DECLARE_SOA_COLUMN(HfOrigin, hfOrigin, int8_t);
} // namespace matchingmctruth
// Includes MCHs with no candidates/invalid labels. No MC selection is applied.
// Sum NumTruePairs over eligible MCHs is the PRE-selection N_pairable.
DECLARE_SOA_TABLE(MatchingMchTruths, "AOD", "FWDMCHTRUTH",
                  matchingmctruth::PopulationRowId, mftmchmatchtruth::MchTrackId,
                  mftmchmatchtruth::MchMcParticleId, mftmchmatchtruth::MchMcMask,
                  mftmchmatchtruth::HasMchMcParticle, mftmchmatchtruth::MchPdgCode,
                  mftmchmatchtruth::MchMcCollisionId, matchingmctruth::NumTruePairs,
                  matchingmctruth::IsDecayMuonMch, matchingmctruth::AncestryStatus,
                  matchingmctruth::ChargedAncestorIds, matchingmctruth::ChargedAncestorPdgs, matchingmctruth::SignalOrigin, matchingmctruth::HfOrigin);

// MC-only inclusive dimuon annotations. No mother/ancestry selection is applied.
// Reconstruction inputs and matching denominators are never filtered.
namespace mcdimuontruth
{
DECLARE_SOA_COLUMN(MchMcParticleId1, mchMcParticleId1, int64_t);
DECLARE_SOA_COLUMN(MchMcParticleId2, mchMcParticleId2, int64_t);
DECLARE_SOA_COLUMN(CommonJpsiId, commonJpsiId, int64_t);
DECLARE_SOA_COLUMN(IsOppositeSignMuonPair, isOppositeSignMuonPair, bool);
DECLARE_SOA_COLUMN(IsCommonJpsi, isCommonJpsi, bool);
DECLARE_SOA_COLUMN(AncestryComplete, ancestryComplete, bool);
DECLARE_SOA_COLUMN(NumCorrectMatches, numCorrectMatches, int8_t); // -1 for MCH-only
// Direct mothers are distinct from the recursively searched CommonJpsiId.
// These vectors are filled only for reconstructed ULS (PairCharge == 0).
DECLARE_SOA_COLUMN(MchMotherIds1, mchMotherIds1, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MchMotherPdgs1, mchMotherPdgs1, std::vector<int>);
DECLARE_SOA_COLUMN(MchMotherIds2, mchMotherIds2, std::vector<int64_t>);
DECLARE_SOA_COLUMN(MchMotherPdgs2, mchMotherPdgs2, std::vector<int>);
DECLARE_SOA_COLUMN(CommonMotherIds, commonMotherIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(CommonMotherPdgs, commonMotherPdgs, std::vector<int>);
DECLARE_SOA_COLUMN(MotherLabelsComplete, motherLabelsComplete, bool);

// Bit 0..9: eta, etaPrime, rho0, omega, phi, Jpsi, psi2S,
// Upsilon1S, Upsilon2S, Upsilon3S. Zero means no validated direct-mother tag.
constexpr uint16_t resonanceMotherBit(int pdg)
{
  switch (pdg) {
    case 221:
      return 1u << 0;
    case 331:
      return 1u << 1;
    case 113:
      return 1u << 2;
    case 223:
      return 1u << 3;
    case 333:
      return 1u << 4;
    case 443:
      return 1u << 5;
    case 100443:
      return 1u << 6;
    case 553:
      return 1u << 7;
    case 100553:
      return 1u << 8;
    case 200553:
      return 1u << 9;
    default:
      return 0;
  }
}
DECLARE_SOA_COLUMN(ResonanceMotherMask, resonanceMotherMask, uint16_t);
DECLARE_SOA_COLUMN(ResonanceMotherIds, resonanceMotherIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(ResonanceMotherPdgs, resonanceMotherPdgs, std::vector<int>);

} // namespace mcdimuontruth
DECLARE_SOA_TABLE(McDimuonTruths, "AOD", "FWDDIMUTRUTH",
                  mchdimuonpv::CollisionId, selectedmatchingmuon::Method,
                  mchdimuonpv::MchTrackId1, mlmatchedpv::MftTrackId1,
                  mchdimuonpv::MchTrackId2, mlmatchedpv::MftTrackId2,
                  mcdimuontruth::MchMcParticleId1, mcdimuontruth::MchMcParticleId2,
                  mcdimuontruth::CommonJpsiId, mcdimuontruth::IsOppositeSignMuonPair,
                  mcdimuontruth::IsCommonJpsi, mcdimuontruth::AncestryComplete,
                  mcdimuontruth::NumCorrectMatches, mchdimuonpv::PairCharge,
                  mlmatchedpv::MlGoodScore1, mlmatchedpv::MlGoodScore2,
                  mchdimuonpv::Mass, mchdimuonpv::PairPt, mchdimuonpv::PairRapidity,
                  mcdimuontruth::MchMotherIds1, mcdimuontruth::MchMotherPdgs1,
                  mcdimuontruth::MchMotherIds2, mcdimuontruth::MchMotherPdgs2,
                  mcdimuontruth::CommonMotherIds, mcdimuontruth::CommonMotherPdgs,
                  mcdimuontruth::MotherLabelsComplete, mcdimuontruth::ResonanceMotherMask,
                  mcdimuontruth::ResonanceMotherIds, mcdimuontruth::ResonanceMotherPdgs);

// Internal audit transport shares the public schema, but is never an AOD output.
O2HASH("InternalMatchingMchPopulation");
DECLARE_SOA_STAGE(InternalMatchingMchPopulation, "AOD2", "FWDMCHPOP", 0);
O2HASH("InternalMatchingCandidates");
DECLARE_SOA_STAGE(InternalMatchingCandidates, "AOD2", "FWDMPCAND", 0);
O2HASH("InternalMatchingMchOutcomes");
DECLARE_SOA_STAGE(InternalMatchingMchOutcomes, "AOD2", "FWDMCHOUT", 0);
O2HASH("InternalSelectedMatchingMuons");
DECLARE_SOA_STAGE(InternalSelectedMatchingMuons, "AOD2", "FWDEVALMUON", 0);
O2HASH("InternalMatchingDimuonVertices");
DECLARE_SOA_STAGE(InternalMatchingDimuonVertices, "AOD2", "FWDDIMUVTX", 0);
O2HASH("InternalMchDimuonsPV");
DECLARE_SOA_STAGE(InternalMchDimuonsPV, "AOD2", "FWDMCHPAIR", 0);
O2HASH("InternalMlMatchedMuonsPV");
DECLARE_SOA_STAGE(InternalMlMatchedMuonsPV, "AOD2", "FWDMLMUONPV", 0);
O2HASH("InternalMlMatchedDimuonsPV");
DECLARE_SOA_STAGE(InternalMlMatchedDimuonsPV, "AOD2", "FWDMLPAIRPV", 0);
O2HASH("InternalMatchingCandidateTruths");
DECLARE_SOA_STAGE(InternalMatchingCandidateTruths, "AOD2", "FWDCANDTRUTH", 0);
O2HASH("InternalMatchingCandidateRelations");
DECLARE_SOA_STAGE(InternalMatchingCandidateRelations, "AOD2", "FWDCANDREL", 0);
O2HASH("InternalMatchingMchTruths");
DECLARE_SOA_STAGE(InternalMatchingMchTruths, "AOD2", "FWDMCHTRUTH", 0);
O2HASH("InternalMixedMftMchMatchTruths");
DECLARE_SOA_STAGE(InternalMixedMftMchMatchTruths, "AOD2", "FWDMPMIXLBL", 0);
O2HASH("InternalMcDimuonTruths");
DECLARE_SOA_STAGE(InternalMcDimuonTruths, "AOD2", "FWDDIMUTRUTH", 0);

} // namespace o2::aod
#endif // COMMON_DATAMODEL_MATCHMFTMUONDATA_H_
