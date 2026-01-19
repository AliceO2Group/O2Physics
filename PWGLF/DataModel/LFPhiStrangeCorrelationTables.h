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

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include <Framework/ASoA.h>

namespace o2::aod
{
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

DECLARE_SOA_TABLE(PhimesonCandidatesData, "AOD", "PHICANDIDATESDATA",
                  lf_selection_phi_candidate::CollisionId,
                  lf_selection_phi_candidate::M,
                  lf_selection_phi_candidate::Pt,
                  lf_selection_phi_candidate::Y,
                  lf_selection_phi_candidate::Phi,
                  lf_selection_phi_candidate::InMassRegion<lf_selection_phi_candidate::M>);

DECLARE_SOA_TABLE(PhimesonCandidatesMcReco, "AOD", "PHICANDIDATESMCRECO",
                  lf_selection_phi_candidate::CollisionId,
                  lf_selection_phi_candidate::M,
                  lf_selection_phi_candidate::Pt,
                  lf_selection_phi_candidate::Y,
                  lf_selection_phi_candidate::Phi,
                  lf_selection_phi_candidate::InMassRegion<lf_selection_phi_candidate::M>);

DECLARE_SOA_TABLE(PhimesonCandidatesMcGen, "AOD", "PHICANDIDATESMCGEN",
                  lf_selection_phi_candidate::CollisionId,
                  lf_selection_phi_candidate::M,
                  lf_selection_phi_candidate::Pt,
                  lf_selection_phi_candidate::Y,
                  lf_selection_phi_candidate::Phi,
                  lf_selection_phi_candidate::InMassRegion<lf_selection_phi_candidate::M>);

namespace lf_selection_phi_collision
{
DECLARE_SOA_COLUMN(PhimesonSel, phimesonSel, bool);
} // namespace lf_selection_phi_collision

DECLARE_SOA_TABLE(PhimesonSelectionData, "AOD", "PHIINCOLLDATA",
                  lf_selection_phi_collision::PhimesonSel);

DECLARE_SOA_TABLE(PhimesonSelectionMcGen, "AOD", "PHIINCOLLMCGEN",
                  lf_selection_phi_collision::PhimesonSel);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFPHISTRANGECORRELATIONTABLES_H_
