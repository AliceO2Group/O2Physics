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

/// \file EMCConversionCandidate.h
/// \brief Header file that defines EMCConversionCandidate, a struct to be used with the EMCConversion ML model
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>

#ifndef PWGEM_PHOTONMESON_CORE_EMCCONVERSIONCANDIDATE_H_
#define PWGEM_PHOTONMESON_CORE_EMCCONVERSIONCANDIDATE_H_

#include <concepts>

namespace o2::analysis::em
{

// Requires T to expose the specific named getters getInputFeatures() needs,
// each returning something convertible to float.
template <typename T>
concept IsEmcConversionCandidate = requires(T const& c) {
  { c.minv() } -> std::convertible_to<float>;
  { c.deltaEta() } -> std::convertible_to<float>;
  { c.deltaR() } -> std::convertible_to<float>;
  { c.phiv() } -> std::convertible_to<float>;
  { c.rConv() } -> std::convertible_to<float>;
  { c.totE() } -> std::convertible_to<float>;
  { c.e2() } -> std::convertible_to<float>;
  { c.e1() } -> std::convertible_to<float>;
  { c.deltaPhi() } -> std::convertible_to<float>;
};

struct EMCConversionCandidate {
  float mMinv;
  float mDeltaEta;
  float mDeltaR;
  float mPhiv;
  float mRConv;
  float mTotE;
  float mE2;
  float mE1;
  float mDeltaPhi;

  [[nodiscard]] float minv() const { return mMinv; }
  [[nodiscard]] float deltaEta() const { return mDeltaEta; }
  [[nodiscard]] float deltaR() const { return mDeltaR; }
  [[nodiscard]] float phiv() const { return mPhiv; }
  [[nodiscard]] float rConv() const { return mRConv; }
  [[nodiscard]] float totE() const { return mTotE; }
  [[nodiscard]] float e2() const { return mE2; }
  [[nodiscard]] float e1() const { return mE1; }
  [[nodiscard]] float deltaPhi() const { return mDeltaPhi; }
};

} // namespace o2::analysis::em

#endif // PWGEM_PHOTONMESON_CORE_EMCCONVERSIONCANDIDATE_H_
