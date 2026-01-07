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

#ifndef PWGDQ_CORE_MCSIGNALLIBRARY_H_
#define PWGDQ_CORE_MCSIGNALLIBRARY_H_

#include "PWGDQ/Core/MCProng.h"
#include "PWGDQ/Core/MCSignal.h"

#include "rapidjson/document.h"

#include <string>
#include <vector>

namespace o2::aod
{
namespace dqmcsignals
{
MCSignal* GetMCSignal(const char* signalName);

std::vector<MCSignal*> GetMCSignalsFromJSON(const char* json);

template <typename T>
bool ValidateJSONMCSignal(T sigJSON, const char* sigName);

template <typename T>
MCProng* ParseJSONMCProng(T prongJSON, const char* prongName);

template <typename T>
bool ValidateJSONMCProng(T prongJSON, const char* prongName);
} // namespace dqmcsignals
} // namespace o2::aod

#endif // PWGDQ_CORE_MCSIGNALLIBRARY_H_
