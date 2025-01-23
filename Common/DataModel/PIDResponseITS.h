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
/// \author Giorgio Alberto Lucia giorgio.alberto.lucia@cern.ch
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
    float sum = 0;
    int nclusters = 0;
    int max = 0;
    for (int layer = 0; layer < 7; layer++) {
      int clsize = (itsClusterSizes >> (layer * 4)) & 0xf;
      if (clsize > 0) {
        nclusters++;
        sum += clsize;
        if (clsize > max) {
          max = clsize;
        }
      }
    }
    if (nclusters == 0) {
      return 0;
    }
    // truncated mean
    return (sum - max) / (nclusters - 1);
  };

  template <o2::track::PID::ID id>
  static float expSignal(const float momentum)
  {
    static constexpr float inverseMass = 1. / o2::track::pid_constants::sMasses[id];
    // static constexpr float charge = static_cast<float>(o2::track::pid_constants::sCharges[id]);
    const float bg = momentum * inverseMass;
    if (id == o2::track::PID::Helium3 || id == o2::track::PID::Alpha) {
      return (mITSRespParamsZ2[0] / (std::pow(bg, mITSRespParamsZ2[1])) + mITSRespParamsZ2[2]);
    }
    return (mITSRespParams[0] / (std::pow(bg, mITSRespParams[1])) + mITSRespParams[2]);
  }

  template <o2::track::PID::ID id>
  static float expResolution(const float momentum)
  {
    static constexpr float inverseMass = 1. / o2::track::pid_constants::sMasses[id];
    // static constexpr float charge = static_cast<float>(o2::track::pid_constants::sCharges[id]);
    const float bg = momentum * inverseMass;
    if (id == o2::track::PID::Helium3 || id == o2::track::PID::Alpha) {
      return mResolutionParamsZ2[0] * std::erf((bg - mResolutionParamsZ2[1]) / mResolutionParamsZ2[2]);
    }
    return mResolutionParams[0] * std::erf((bg - mResolutionParams[1]) / mResolutionParams[2]);
  }

  template <o2::track::PID::ID id>
  static float nSigmaITS(uint32_t itsClusterSizes, float momentum, float eta)
  {
    const float exp = expSignal<id>(momentum);
    const float average = averageClusterSize(itsClusterSizes);
    const float coslInv = 1. / std::cosh(eta);
    const float resolution = expResolution<id>(momentum) * exp;
    return (average * coslInv - exp) / resolution;
  };

  template <o2::track::PID::ID id, typename T>
  static float nSigmaITS(const T& track)
  {
    return nSigmaITS<id>(track.itsClusterSizes(), track.p(), track.eta());
  }

  static void setParameters(float p0, float p1, float p2,
                            float p0_Z2, float p1_Z2, float p2_Z2,
                            float p0_res, float p1_res, float p2_res,
                            float p0_res_Z2, float p1_res_Z2, float p2_res_Z2)
  {
    if (mIsInitialized) {
      LOG(fatal) << "ITSResponse parameters already initialized";
    }
    mIsInitialized = true;
    mITSRespParams[0] = p0;
    mITSRespParams[1] = p1;
    mITSRespParams[2] = p2;
    mITSRespParamsZ2[0] = p0_Z2;
    mITSRespParamsZ2[1] = p1_Z2;
    mITSRespParamsZ2[2] = p2_Z2;
    mResolutionParams[0] = p0_res;
    mResolutionParams[1] = p1_res;
    mResolutionParams[2] = p2_res;
    mResolutionParamsZ2[0] = p0_res_Z2;
    mResolutionParamsZ2[1] = p1_res_Z2;
    mResolutionParamsZ2[2] = p2_res_Z2;
  }

 private:
  static std::array<float, 3> mITSRespParams;
  static std::array<float, 3> mITSRespParamsZ2;
  static std::array<float, 3> mResolutionParams;
  static std::array<float, 3> mResolutionParamsZ2;
  static bool mIsInitialized;
};

std::array<float, 3> ITSResponse::mITSRespParams = {1.18941, 1.53792, 1.69961};
std::array<float, 3> ITSResponse::mITSRespParamsZ2 = {2.35117, 1.80347, 5.14355};
// relative resolution is modelled with an erf function: [0]*TMath::Erf((x-[1])/[2])
std::array<float, 3> ITSResponse::mResolutionParams = {1.94669e-01, -2.08616e-01, 1.30753};
std::array<float, 3> ITSResponse::mResolutionParamsZ2 = {8.74371e-02, -1.82804, 5.06449e-01};
bool ITSResponse::mIsInitialized = false;

namespace pidits
{
DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaElImp, itsNSigmaEl, //! Nsigma separation with the ITS detector for electrons
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Electron>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaMuImp, itsNSigmaMu, //! Nsigma separation with the ITS detector for muons
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Muon>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaPiImp, itsNSigmaPi, //! Nsigma separation with the ITS detector for pions
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Pion>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaKaImp, itsNSigmaKa, //! Nsigma separation with the ITS detector for kaons
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Kaon>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaPrImp, itsNSigmaPr, //! Nsigma separation with the ITS detector for protons
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Proton>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaDeImp, itsNSigmaDe, //! Nsigma separation with the ITS detector for deuterons
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Deuteron>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaTrImp, itsNSigmaTr, //! Nsigma separation with the ITS detector for tritons
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Triton>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaHeImp, itsNSigmaHe, //! Nsigma separation with the ITS detector for helium3
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Helium3>(itsClusterSizes, momentum, eta);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(ITSNSigmaAlImp, itsNSigmaAl, //! Nsigma separation with the ITS detector for alphas
                           [](uint32_t itsClusterSizes, float momentum, float eta) -> float {
                             return ITSResponse::nSigmaITS<o2::track::PID::Alpha>(itsClusterSizes, momentum, eta);
                           });

// Define user friendly names for the columns to join with the tracks
using ITSNSigmaEl = ITSNSigmaElImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaMu = ITSNSigmaMuImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaPi = ITSNSigmaPiImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaKa = ITSNSigmaKaImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaPr = ITSNSigmaPrImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaDe = ITSNSigmaDeImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaTr = ITSNSigmaTrImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaHe = ITSNSigmaHeImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;
using ITSNSigmaAl = ITSNSigmaAlImp<o2::aod::track::ITSClusterSizes, o2::aod::track::P, o2::aod::track::Eta>;

} // namespace pidits
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSEITS_H_
