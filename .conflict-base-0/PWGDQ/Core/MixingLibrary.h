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
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/VarManager.h"

#include <TString.h>

namespace o2::aod
{
namespace dqmixing
{
void SetUpMixing(MixingHandler* mh, const char* mixingVarible);
} // namespace dqmixing
} // namespace o2::aod

void o2::aod::dqmixing::SetUpMixing(MixingHandler* mh, const char* mixingVarible);
