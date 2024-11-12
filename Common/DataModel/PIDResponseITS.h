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
  float mResolution = 0.07f;

  template <o2::track::PID::ID id>
  float expSignal(const float momentum)
  {
    static constexpr float invmass = o2::track::pid_constants::sMasses[id];
    static constexpr float charge = static_cast<float>(o2::track::pid_constants::sCharges[id]);
    return o2::tpc::BetheBlochAleph(momentum * invmass,
                                    mBetheBlochParams[0],
                                    mBetheBlochParams[1],
                                    mBetheBlochParams[2],
                                    mBetheBlochParams[3],
                                    mBetheBlochParams[4]) *
           std::pow(charge, mChargeFactor);
  }
} mITSParams;

float averageClusterSize(uint32_t itsClusterSizes)
{
  float average = 0;
  int nclusters = 0;

  for (int layer = 0; layer < 7; layer++) {
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

namespace pidits
{
DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaElImp, itsNSigmaEl, //! Nsigma separation with the ITS detector for electrons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Electron>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaMuImp, itsNSigmaMu, //! Nsigma separation with the ITS detector for muons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Muon>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaPiImp, itsNSigmaPi, //! Nsigma separation with the ITS detector for pions
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Pion>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaKaImp, itsNSigmaKa, //! Nsigma separation with the ITS detector for kaons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Kaon>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaPrImp, itsNSigmaPr, //! Nsigma separation with the ITS detector for protons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Proton>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaDeImp, itsNSigmaDe, //! Nsigma separation with the ITS detector for deuterons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Deuteron>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaTrImp, itsNSigmaTr, //! Nsigma separation with the ITS detector for tritons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Triton>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaHeImp, itsNSigmaHe, //! Nsigma separation with the ITS detector for helium3
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Helium3>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaAlImp, itsNSigmaAl, //! Nsigma separation with the ITS detector for alphas
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             const float bethe = mITSParams.expSignal<o2::track::PID::Alpha>(momentum);
                             const float average = averageClusterSize(itsClusterSizes);
                             const float resolution = mITSParams.mResolution * bethe;
                             return (average - bethe) / resolution;
                           });

#define ITSNSigmaEl ITSNSigmaElImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaMu ITSNSigmaMuImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaPi ITSNSigmaPiImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaKa ITSNSigmaKaImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaPr ITSNSigmaPrImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaDe ITSNSigmaDeImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaTr ITSNSigmaTrImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaHe ITSNSigmaHeImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>
#define ITSNSigmaAl ITSNSigmaAlImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>

} // namespace pidits
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSEITS_H_
