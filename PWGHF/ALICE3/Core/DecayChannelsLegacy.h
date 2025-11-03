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

#ifndef PWGHF_ALICE3_CORE_DECAYCHANNELSLEGACY_H_
#define PWGHF_ALICE3_CORE_DECAYCHANNELSLEGACY_H_

namespace o2::aod
{
namespace hf_cand_x
{
enum DecayType {
  XToJpsiToEEPiPi = 0,
  XToJpsiToMuMuPiPi
};
} // namespace hf_cand_x

namespace hf_cand_xicc
{
enum DecayType {
  XiccToXicPi = 0
};
} // namespace hf_cand_xicc

namespace hf_cand_chic
{
enum DecayType {
  ChicToJpsiToEEGamma = 0,
  ChicToJpsiToMuMuGamma
};
} // namespace hf_cand_chic

} // namespace o2::aod

#endif // PWGHF_ALICE3_CORE_DECAYCHANNELSLEGACY_H_
