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

namespace o2::aod
{

struct ITSResponse {
  static float averageClusterSize(uint32_t itsClusterSizes)
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

  template <o2::track::PID::ID id>
  static float expSignal(const float momentum)
  {
    static constexpr float inverseMass = 1. / o2::track::pid_constants::sMasses[id];
    static constexpr float charge = static_cast<float>(o2::track::pid_constants::sCharges[id]);
    const float bg = momentum * inverseMass;
    return (mITSRespParams[0] / (std::pow(bg, mITSRespParams[1])) + mITSRespParams[2]) * std::pow(charge, mChargeFactor);
  }

  template <o2::track::PID::ID id>
  static float nSigmaITS(uint32_t itsClusterSizes, float momentum)
  {
    const float exp = expSignal<id>(momentum);
    const float average = averageClusterSize(itsClusterSizes);
    const float resolution = mResolution * exp;
    return (average - exp) / resolution;
  };

  static void setParameters(float p0, float p1, float p2, float chargeFactor, float resolution)
  {
    if (mIsInitialized) {
      LOG(fatal) << "ITSResponse parameters already initialized";
    }
    mIsInitialized = true;
    mITSRespParams[0] = p0;
    mITSRespParams[1] = p1;
    mITSRespParams[2] = p2;
    mChargeFactor = chargeFactor;
    mResolution = resolution;
  }

 private:
  static std::array<float, 3> mITSRespParams;
  static float mChargeFactor;
  static float mResolution;
  static bool mIsInitialized;
};

std::array<float, 3> ITSResponse::mITSRespParams = {0.903, 2.014, 2.440};
float ITSResponse::mChargeFactor = 2.299999952316284f;
float ITSResponse::mResolution = 0.15f;
bool ITSResponse::mIsInitialized = false;

namespace pidits
{
DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaElImp, itsNSigmaEl, //! Nsigma separation with the ITS detector for electrons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Electron>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaMuImp, itsNSigmaMu, //! Nsigma separation with the ITS detector for muons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Muon>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaPiImp, itsNSigmaPi, //! Nsigma separation with the ITS detector for pions
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Pion>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaKaImp, itsNSigmaKa, //! Nsigma separation with the ITS detector for kaons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Kaon>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaPrImp, itsNSigmaPr, //! Nsigma separation with the ITS detector for protons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Proton>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaDeImp, itsNSigmaDe, //! Nsigma separation with the ITS detector for deuterons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Deuteron>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaTrImp, itsNSigmaTr, //! Nsigma separation with the ITS detector for tritons
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Triton>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaHeImp, itsNSigmaHe, //! Nsigma separation with the ITS detector for helium3
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Helium3>(itsClusterSizes, momentum);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaAlImp, itsNSigmaAl, //! Nsigma separation with the ITS detector for alphas
                           [](uint32_t itsClusterSizes, float momentum) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Alpha>(itsClusterSizes, momentum);
                           });

// Define user friendly names for the columns to join with the tracks
using ITSNSigmaEl = ITSNSigmaElImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaMu = ITSNSigmaMuImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaPi = ITSNSigmaPiImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaKa = ITSNSigmaKaImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaPr = ITSNSigmaPrImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaDe = ITSNSigmaDeImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaTr = ITSNSigmaTrImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaHe = ITSNSigmaHeImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;
using ITSNSigmaAl = ITSNSigmaAlImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P>;

} // namespace pidits
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSEITS_H_
