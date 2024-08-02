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

/// \file CentralityEstimation.h
/// \brief Definitions of centrality estimators
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_CORE_CENTRALITYESTIMATION_H_
#define PWGHF_CORE_CENTRALITYESTIMATION_H_

namespace o2::hf_centrality
{
// centrality selection estimators
enum CentralityEstimator {
  None = 0,
  FT0A,
  FT0C,
  FT0M,
  FV0A,
  NTracksPV,
  NCentralityEstimators
};
} // namespace o2::hf_centrality

#endif // PWGHF_CORE_CENTRALITYESTIMATION_H_
