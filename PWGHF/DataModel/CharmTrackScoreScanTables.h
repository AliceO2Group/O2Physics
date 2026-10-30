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

/// \file CharmTrackScoreScanTables.h
/// \brief Per-candidate ML scores, so a track-level working-point scan can be redone offline
///
/// Written by dplusCombinatorics.cxx and consumed by dplusCombinatoricsQa.cxx. It lives in a
/// header rather than being declared twice because the two declarations must not drift: a
/// mismatch would be read back as silently reinterpreted columns, not as an error.
///
/// The point of the table is that one *loose* run gives the whole scan. A candidate is built iff
/// every prong clears its role's threshold, and thresholds only rise as the working point
/// tightens, so the candidate set at any tighter point is a strict subset of the loose one. Given
/// the per-role minimum score each candidate was built from, membership at a tighter point is a
/// comparison, not another run.
///
/// Scores are per *role*, not one global minimum: the pion and kaon thresholds sit at different
/// absolute values (0.027 vs 0.032 at the 90% working point), so a single minimum over all three
/// prongs could not be compared against either of them.
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Universita and INFN Torino

#ifndef PWGHF_DATAMODEL_CHARMTRACKSCORESCANTABLES_H_
#define PWGHF_DATAMODEL_CHARMTRACKSCORESCANTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
/// Score written for a hypothesis a candidate does not have. Below every real score (which lie in
/// [0, 1]) and below every threshold, so such a hypothesis never survives a scan point.
constexpr float kScoreScanNone = -1.f;

namespace hf_cand_score_scan
{
DECLARE_SOA_COLUMN(PtCand, ptCand, float);                           //! pT of the candidate
DECLARE_SOA_COLUMN(MinScorePionDplus, minScorePionDplus, float);     //! D+: min pion score over the two pion prongs
DECLARE_SOA_COLUMN(ScoreKaonDplus, scoreKaonDplus, float);           //! D+: kaon score of the opposite-sign prong
DECLARE_SOA_COLUMN(MinScoreKaonDsHypo0, minScoreKaonDsHypo0, float); //! Ds (K K pi): min kaon score over the two kaon prongs
DECLARE_SOA_COLUMN(ScorePionDsHypo0, scorePionDsHypo0, float);       //! Ds (K K pi): pion score of the pion prong
DECLARE_SOA_COLUMN(MinScoreKaonDsHypo1, minScoreKaonDsHypo1, float); //! Ds (pi K K): min kaon score over the two kaon prongs
DECLARE_SOA_COLUMN(ScorePionDsHypo1, scorePionDsHypo1, float);       //! Ds (pi K K): pion score of the pion prong
} // namespace hf_cand_score_scan

DECLARE_SOA_TABLE(HfCandScoreScans, "AOD", "HFCANDSCORE", //! One row per aod::Hf3Prongs row, in the same order
                  hf_cand_score_scan::PtCand,
                  hf_cand_score_scan::MinScorePionDplus,
                  hf_cand_score_scan::ScoreKaonDplus,
                  hf_cand_score_scan::MinScoreKaonDsHypo0,
                  hf_cand_score_scan::ScorePionDsHypo0,
                  hf_cand_score_scan::MinScoreKaonDsHypo1,
                  hf_cand_score_scan::ScorePionDsHypo1);
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_CHARMTRACKSCORESCANTABLES_H_
