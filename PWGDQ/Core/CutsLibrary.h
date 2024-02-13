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

#ifndef PWGDQ_CORE_CUTSLIBRARY_H_
#define PWGDQ_CORE_CUTSLIBRARY_H_

#include <string>
#include <vector>
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/VarManager.h"

namespace o2::aod
{
namespace dqcuts
{
AnalysisCompositeCut* GetCompositeCut(const char* cutName);
AnalysisCut* GetAnalysisCut(const char* cutName);
} // namespace dqcuts
} // namespace o2::aod

AnalysisCompositeCut* o2::aod::dqcuts::GetCompositeCut(const char* cutName);
#endif // PWGDQ_CORE_CUTSLIBRARY_H_
