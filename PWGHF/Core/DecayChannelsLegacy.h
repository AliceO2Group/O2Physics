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

/// \file DecayChannelsLegacy.h
/// \brief Legacy definitions of constants for MC flagging of HF decay channels.
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University
/// \note Should be replaced with corresponding constants in DecayChannels.h.

#ifndef PWGHF_CORE_DECAYCHANNELSLEGACY_H_
#define PWGHF_CORE_DECAYCHANNELSLEGACY_H_

namespace o2::aod
{
namespace hf_cand_xic0_omegac0
{
enum DecayType {
  XiczeroToXiPi = 0,
  OmegaczeroToXiPi,
  OmegaczeroToOmegaPi,
  OmegaczeroToOmegaK,
  OmegaczeroToOmegaPiOneMu
};
} // namespace hf_cand_xic0_omegac0

namespace hf_cand_xic_to_xi_pi_pi
{
enum DecayType {
  XicToXiPiPi = 0,      // Ξc± → Ξ∓ π± π±
  XicToXiResPiToXiPiPi, // Ξc± → Ξ(1530) π± → Ξ∓ π± π±
  NDecayType
};
} // namespace hf_cand_xic_to_xi_pi_pi

namespace hf_cand_dstar
{
enum DecayType {
  DstarToD0Pi = 0,
  D0ToPiK,
  DstarToD0PiPi0,
  D0ToPiKPi0,
  NDstarDecayType
};
} // namespace hf_cand_dstar

namespace hf_cand_sigmac
{
enum DecayType {
  Sc0ToPKPiPi = 0,
  ScplusplusToPKPiPi,
  ScStar0ToPKPiPi,
  ScStarPlusPlusToPKPiPi
};
} // namespace hf_cand_sigmac

namespace hf_cand_b0
{
enum DecayType {
  B0ToDPi = 0,
  B0ToDstarPi
};
} // namespace hf_cand_b0

namespace hf_cand_bplus
{
enum DecayType {
  BplusToD0Pi = 0
};
} // namespace hf_cand_bplus

namespace hf_cand_bs
{
enum DecayType {
  BsToDsPi = 0
};
} // namespace hf_cand_bs

namespace hf_cand_lb
{
enum DecayType {
  LbToLcPi
};
} // namespace hf_cand_lb

} // namespace o2::aod

#endif // PWGHF_CORE_DECAYCHANNELSLEGACY_H_
