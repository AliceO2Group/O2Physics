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

///
/// \file   PIDResponseITS.h
/// \since  2024-11-12
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesco Mazzaschi francesco.mazzaschi@cern.ch
/// \brief  Set of tables, tasks and utilities to provide the interface between
///         the analysis data model and the PID response of the ITS
///

#ifndef COMMON_DATAMODEL_PIDRESPONSEITS_H_
#define COMMON_DATAMODEL_PIDRESPONSEITS_H_

// O2 includes
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/Logger.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

namespace o2::aod
{

struct ITSParams {
  std::array<float, 5> mBetheBlochParams = {0.03209809958934784, 19.9768009185791, 2.5266601063857674e-16, 2.7212300300598145, 6.080920219421387};
  float mChargeFactor = 2.299999952316284f;
} mITSParams;

float averageClusterSize(uint32_t itsClusterSizes)
{
  float average = 0;
  int nclusters = 0;

  for (int i = 0; i < 7; i++) {
    if ((itsClusterSizes >> (layer * 4)) & 0xf) {
      nclusters++;
      average += (itsClusterSizes >> (layer * 4)) & 0xf;
    }
  }
  if (nclusters == 0) {
    return 0;
  }
  return average / nclusters;
};

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaEl, itsNSigmaEl, //! E Nsigma separation with the ITS detector for electron
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             const o2cp::PID::ID id = o2cp::PID::Electron;
                             const float bethe = o2::tpc::BetheBlochAleph(momentum / o2::track::pid_constants::sMasses[id],
                                                                          mITSParams.mBetheBlochParams[0], mITSParams.mBetheBlochParams[1], mITSParams.mBetheBlochParams[2], mITSParams.mBetheBlochParams[3], mITSParams.mBetheBlochParams[4]) *
                                                 std::pow(static_cast<float>(o2::track::pid_constants::sCharges[id]), mITSParams.mChargeFactor);
                             const float resolution = 0.07 * bethe;
                             return (average - bethe) / resolution;
                           });
DECLARE_SOA_TABLE(pidITSEl, "AOD", "pidITSEl", //! Table of the ITS response with expected signal, expected resolution and Nsigma for electron
                  pidtof::ITSNSigmaEl<track::ITSClusterSizes, track::P, track::Eta>);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSEITS_H_
