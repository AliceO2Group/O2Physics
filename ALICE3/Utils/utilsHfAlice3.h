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

/// \file utilsHfAlice3.h
/// \brief Class with helper functions for HF analyses
///
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN Turin

#ifndef ALICE3_UTILS_UTILSHFALICE3_H_
#define ALICE3_UTILS_UTILSHFALICE3_H_

#include "PWGHF/Core/HfHelper.h"

namespace o2::analysis
{

enum CharmHadAlice3 { Lc = 1 };

} // namespace o2::analysis

struct HfHelperAlice3 {
  /// Get candidate mass (ALICE3 HF data model)
  /// \tparam TCand candidate type
  /// \param cand candidate
  /// \return candidate mass
  template <o2::analysis::CharmHadAlice3 CharmHad, bool SwapHypo, typename TCand>
  static double getCandMass(const TCand& cand)
  {
    switch (CharmHad) {
      case o2::analysis::CharmHadAlice3::Lc:
        return SwapHypo ? HfHelper::invMassLcToPiKP(cand) : HfHelper::invMassLcToPKPi(cand);
      default:
        LOG(fatal) << "Unsupported charm hadron type";
        return -1.;
    }
  }

  /// Get candidate energy (ALICE3 HF data model)
  /// \tparam TCand candidate type
  /// \param cand candidate
  /// \return candidate energy
  template <o2::analysis::CharmHadAlice3 CharmHad, typename TCand>
  static double getCandEnergy(const TCand& cand)
  {
    switch (CharmHad) {
      case o2::analysis::CharmHadAlice3::Lc:
        return HfHelper::eLc(cand);
      default:
        LOG(fatal) << "Unsupported charm hadron type";
        return -1.;
    }
  }

  /// Get candidate rapidity (ALICE3 HF data model)
  /// \tparam TCand candidate type
  /// \param cand candidate
  /// \return candidate rapidity
  template <o2::analysis::CharmHadAlice3 CharmHad, typename TCand>
  static double getCandY(const TCand& cand)
  {
    if constexpr (requires { cand.flagMcRec(); }) {
      switch (CharmHad) {
        case o2::analysis::CharmHadAlice3::Lc:
          return HfHelper::yLc(cand);
        default:
          LOG(fatal) << "Unsupported charm hadron type";
          return -1.;
      }
    } else {
      switch (CharmHad) {
        case o2::analysis::CharmHadAlice3::Lc:
          return RecoDecay::y(cand.pVector(), o2::constants::physics::MassLambdaCPlus);
        default:
          LOG(fatal) << "Unsupported charm hadron type";
          return -1.;
      }
    }
  }
};

#endif // ALICE3_UTILS_UTILSHFALICE3_H_
