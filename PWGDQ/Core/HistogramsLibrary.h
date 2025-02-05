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
//
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//

#ifndef PWGDQ_CORE_HISTOGRAMSLIBRARY_H_
#define PWGDQ_CORE_HISTOGRAMSLIBRARY_H_

#include <TString.h>
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/VarManager.h"
#include "CommonConstants/MathConstants.h"
#include "rapidjson/document.h"

namespace o2::aod
{
namespace dqhistograms
{
void DefineHistograms(HistogramManager* hm, const char* histClass, const char* groupName, const char* subGroupName = "");
template <typename T>
bool ValidateJSONHistogram(T hist);
void AddHistogramsFromJSON(HistogramManager* hm, const char* json);
}
} // namespace o2::aod

#endif // PWGDQ_CORE_HISTOGRAMSLIBRARY_H_
