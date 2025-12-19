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

/// \file SelectionStudyTables
/// \brief  tables meant to do event selection studies for O-O / light systems
///
/// \author ALICE

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <vector>

#ifndef COMMON_DATAMODEL_SELECTIONSTUDYTABLES_H_
#define COMMON_DATAMODEL_SELECTIONSTUDYTABLES_H_

namespace o2::aod
{
namespace selectionstudy
{
DECLARE_SOA_COLUMN(PtPions, ptPions, std::vector<float>);
DECLARE_SOA_COLUMN(PtKaons, ptKaons, std::vector<float>);
DECLARE_SOA_COLUMN(PtProtons, ptProtons, std::vector<float>);
DECLARE_SOA_COLUMN(PtK0s, ptPK0s, std::vector<float>);
DECLARE_SOA_COLUMN(PtLambdas, ptLambdas, std::vector<float>);
DECLARE_SOA_COLUMN(PtXis, ptXis, std::vector<float>);
DECLARE_SOA_COLUMN(PtOmegas, ptOmegas, std::vector<float>);
DECLARE_SOA_COLUMN(PtPhis, ptPhis, std::vector<float>);
DECLARE_SOA_COLUMN(PtKStars, ptKStars, std::vector<float>);
DECLARE_SOA_COLUMN(PtDs, ptDs, std::vector<float>);
DECLARE_SOA_COLUMN(PtLambdaCs, ptLambdaCs, std::vector<float>);
DECLARE_SOA_COLUMN(PtJPsis, ptJPsis, std::vector<float>);
} // namespace selectionstudy

DECLARE_SOA_TABLE(PIDPts, "AOD", "PIDPTS", o2::soa::Index<>,
                  o2::aod::selectionstudy::PtPions,
                  o2::aod::selectionstudy::PtKaons,
                  o2::aod::selectionstudy::PtProtons,
                  o2::aod::selectionstudy::PtK0s,
                  o2::aod::selectionstudy::PtLambdas,
                  o2::aod::selectionstudy::PtXis,
                  o2::aod::selectionstudy::PtOmegas,
                  o2::aod::selectionstudy::PtPhis,
                  o2::aod::selectionstudy::PtKStars,
                  o2::aod::selectionstudy::PtDs,
                  o2::aod::selectionstudy::PtLambdaCs,
                  o2::aod::selectionstudy::PtJPsis);
} // namespace o2::aod
#endif // COMMON_DATAMODEL_SELECTIONSTUDYTABLES_H_
