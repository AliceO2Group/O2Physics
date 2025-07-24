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

/// \file HfMlResponse.h
/// \brief Class to compute the ML response for HF-analysis selections
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_CORE_HFMLRESPONSE_H_
#define PWGHF_CORE_HFMLRESPONSE_H_

#include "Tools/ML/MlResponse.h"

namespace o2::analysis
{

template <typename TypeOutputScore = float>
class HfMlResponse : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  HfMlResponse() = default;
  /// Default destructor
  virtual ~HfMlResponse() = default;
};

} // namespace o2::analysis

#endif // PWGHF_CORE_HFMLRESPONSE_H_
