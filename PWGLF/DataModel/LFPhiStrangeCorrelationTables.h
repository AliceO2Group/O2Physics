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
/// \file LFPhiStrangeCorrelationTables.h
/// \brief Data model Phi-Strangeness correlation analysis
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#ifndef PWGLF_DATAMODEL_LFPHISTRANGECORRELATIONTABLES_H_
#define PWGLF_DATAMODEL_LFPHISTRANGECORRELATIONTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdlib>

namespace o2::aod
{
namespace lf_selection_event
{
DECLARE_SOA_COLUMN(DefaultSel, defaultSel, bool);
DECLARE_SOA_COLUMN(PhimesonSel, phimesonSel, bool);
} // namespace lf_selection_event

DECLARE_SOA_TABLE(PhiStrangeEvtSelDataLike, "AOD", "EVTSELDATA",
                  lf_selection_event::DefaultSel,
                  lf_selection_event::PhimesonSel);

DECLARE_SOA_TABLE(PhiStrangeEvtSelMcGen, "AOD", "EVTSELMCGEN",
                  lf_selection_event::DefaultSel,
                  lf_selection_event::PhimesonSel);

namespace lf_selection_phi_candidate
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);

DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Phi, phi, float);

DECLARE_SOA_DYNAMIC_COLUMN(InMassRegion, inMassRegion,
                           [](float m, float minM, float maxM) -> bool {
                             return (m >= minM && m <= maxM);
                           });
} // namespace lf_selection_phi_candidate

DECLARE_SOA_TABLE(PhimesonCandidatesData, "AOD", "PHICANDDATA",
                  soa::Index<>,
                  lf_selection_phi_candidate::CollisionId,
                  lf_selection_phi_candidate::M,
                  lf_selection_phi_candidate::Pt,
                  lf_selection_phi_candidate::Y,
                  lf_selection_phi_candidate::Phi,
                  lf_selection_phi_candidate::InMassRegion<lf_selection_phi_candidate::M>);

DECLARE_SOA_TABLE(PhimesonCandidatesMcReco, "AOD", "PHICANDMCREC",
                  soa::Index<>,
                  lf_selection_phi_candidate::CollisionId,
                  lf_selection_phi_candidate::M,
                  lf_selection_phi_candidate::Pt,
                  lf_selection_phi_candidate::Y,
                  lf_selection_phi_candidate::Phi,
                  lf_selection_phi_candidate::InMassRegion<lf_selection_phi_candidate::M>);

DECLARE_SOA_TABLE(PhimesonCandidatesMcGen, "AOD", "PHICANDMCGEN",
                  soa::Index<>,
                  lf_selection_phi_candidate::CollisionId,
                  lf_selection_phi_candidate::M,
                  lf_selection_phi_candidate::Pt,
                  lf_selection_phi_candidate::Y,
                  lf_selection_phi_candidate::Phi,
                  lf_selection_phi_candidate::InMassRegion<lf_selection_phi_candidate::M>);

namespace lf_selection_strange_reduced
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);

DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Phi, phi, float);

DECLARE_SOA_DYNAMIC_COLUMN(InMassRegion, inMassRegion,
                           [](float m, float minM, float maxM) -> bool {
                             return (m >= minM && m <= maxM);
                           });
} // namespace lf_selection_strange_reduced

DECLARE_SOA_TABLE(K0sReducedCandidatesData, "AOD", "K0SCANDDATA",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(K0sReducedCandidatesMcReco, "AOD", "K0SCANDMCREC",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(LambdaReducedCandidatesData, "AOD", "LAMBDACANDDATA",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(LambdaReducedCandidatesMcReco, "AOD", "LAMBDACANDMCREC",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(AntiLambdaReducedCandidatesData, "AOD", "ALAMBCANDDATA",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(AntiLambdaReducedCandidatesMcReco, "AOD", "ALAMBCANDMCREC",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(XiReducedCandidatesData, "AOD", "XICANDDATA",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(XiReducedCandidatesMcReco, "AOD", "XICANDMCREC",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(OmegaReducedCandidatesData, "AOD", "OMEGACANDDATA",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

DECLARE_SOA_TABLE(OmegaReducedCandidatesMcReco, "AOD", "OMEGACANDMCREC",
                  soa::Index<>,
                  lf_selection_strange_reduced::CollisionId,
                  lf_selection_strange_reduced::M,
                  lf_selection_strange_reduced::Pt,
                  lf_selection_strange_reduced::Y,
                  lf_selection_strange_reduced::Phi,
                  lf_selection_strange_reduced::InMassRegion<lf_selection_strange_reduced::M>);

namespace lf_selection_pion_track
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);

DECLARE_SOA_COLUMN(NSigmaTPC, nSigmaTPC, float);
DECLARE_SOA_COLUMN(NSigmaTOF, nSigmaTOF, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);

DECLARE_SOA_DYNAMIC_COLUMN(InNSigmaRegion, inNSigmaRegion,
                           [](float nSigmaTPC, bool hasTOF, float nSigmaTOF, float pidTPCMax, float pidTOFMax) -> bool {
                             if (std::abs(nSigmaTPC) >= pidTPCMax) {
                               return false; // TPC check failed
                             }
                             if (hasTOF && std::abs(nSigmaTOF) >= pidTOFMax) {
                               return false; // TOF check failed
                             }
                             return true;
                           });
} // namespace lf_selection_pion_track

DECLARE_SOA_TABLE(PionTracksData, "AOD", "PITRACKSDATA",
                  soa::Index<>,
                  lf_selection_pion_track::CollisionId,
                  lf_selection_pion_track::NSigmaTPC,
                  lf_selection_pion_track::NSigmaTOF,
                  lf_selection_pion_track::Pt,
                  lf_selection_pion_track::Y,
                  lf_selection_pion_track::Phi,
                  lf_selection_pion_track::HasTOF,
                  lf_selection_pion_track::InNSigmaRegion<lf_selection_pion_track::NSigmaTPC,
                                                          lf_selection_pion_track::HasTOF,
                                                          lf_selection_pion_track::NSigmaTOF>);

DECLARE_SOA_TABLE(PionTracksMcReco, "AOD", "PITRACKSMCREC",
                  soa::Index<>,
                  lf_selection_pion_track::CollisionId,
                  lf_selection_pion_track::NSigmaTPC,
                  lf_selection_pion_track::NSigmaTOF,
                  lf_selection_pion_track::Pt,
                  lf_selection_pion_track::Y,
                  lf_selection_pion_track::Phi,
                  lf_selection_pion_track::HasTOF,
                  lf_selection_pion_track::InNSigmaRegion<lf_selection_pion_track::NSigmaTPC,
                                                          lf_selection_pion_track::HasTOF,
                                                          lf_selection_pion_track::NSigmaTOF>);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPHISTRANGECORRELATIONTABLES_H_
