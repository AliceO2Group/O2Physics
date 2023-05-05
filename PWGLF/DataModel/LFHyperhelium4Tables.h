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
#ifndef PWGLF_DATAMODEL_LFHYHEFOURTABLES_H_
#define PWGLF_DATAMODEL_LFHYHEFOURTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

namespace o2::aod
{
// helper for building
namespace hyhe4tag
{
// Global bool
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueHyHe4, isTrueHyHe4, bool);                     //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiHyHe4, isTrueAntiHyHe4, bool);                 //! PDG checked correctly in MC
// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsHyHe4Candidate, isHyHe4Candidate, bool);                     //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsAntiHyHe4Candidate, isAntiHyHe4Candidate, bool);                 //! compatible with dE/dx hypotheses
}
DECLARE_SOA_TABLE(HyHe4Tags, "AOD", "HYHE4TAGS",
                  hyhe4tag::IsInteresting,
                  hyhe4tag::IsTrueHyHe4,
                  hyhe4tag::IsTrueAntiHyHe4,
                  hyhe4tag::IsHyHe4Candidate,
                  hyhe4tag::IsAntiHyHe4Candidate);
} // namespace o2::aod


#endif // PWGLF_DATAMODEL_LFHYHEFOURTABLES_H_
