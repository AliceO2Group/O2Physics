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
// Contact: daiki.sekihata@cern.ch
//

#ifndef PWGEM_PHOTONMESON_CORE_CUTSLIBRARY_H_
#define PWGEM_PHOTONMESON_CORE_CUTSLIBRARY_H_

#include <string>
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"

namespace o2::aod
{
namespace pcmcuts
{
V0PhotonCut* GetCut(const char* cutName);
} // namespace pcmcuts

namespace phoscuts
{
PHOSPhotonCut* GetCut(const char* cutName);
} // namespace phoscuts

namespace emccuts
{
EMCPhotonCut* GetCut(const char* cutName);
} // namespace emccuts

} // namespace o2::aod
#endif // PWGEM_PHOTONMESON_CORE_CUTSLIBRARY_H_
