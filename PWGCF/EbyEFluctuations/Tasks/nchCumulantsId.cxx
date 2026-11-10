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

/// \file NchCumulantsId.cxx
/// \brief Event by Event conserved charges fluctuations
/// \author Pravata Panigrahi <pravata.panigrahi@cern.ch> :: Sadhana Dash(sadhana@phy.iitb.ac.in)

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
// #include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
// #include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
// #include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TPDGCode.h>
// #include <TProfile.h>
// #include <TProfile2D.h>
#include <TRandom3.h>

// #include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants
using namespace std;

constexpr int IdBitPI = 0; // Identificationi bits for PID checks
constexpr int IdBitKA = 1;
constexpr int IdBitPR = 2;
// constexpr int IdBitEL = 3;
// constexpr int IdBitDE = 4;

constexpr int McBitPI = 0; // MC particle identification bits
constexpr int McBitKA = 1;
constexpr int McBitPR = 2;
constexpr int McBitEL = 3;
constexpr int McBitDE = 4;

constexpr void BITSET(int& mask, int ithBit)
{
  mask |= (1 << ithBit);
}

constexpr bool BITCHECK(int mask, int ithBit)
{
  return (mask & (1 << ithBit)) != 0;
}

enum PidEnum {
  kCh = 0,
  kPi, // dont use kPion, kKaon, as these enumeration
  kKa, // are already defined in $ROOTSYS/root/include/TPDGCode.h
  kPr,
  kEl,
  kDe
};

enum ChargeEnum {
  kPos = 0,
  kNeg = 1
};

enum class FCPrefixEnum {
  Pr = 0,
  APr = 1,
  PiPos = 2,
  PiNeg = 3,
  KaPos = 4,
  KaNeg = 5,
  Pos = 6,
  Neg = 7,
  NetPi = 8,
  NetKa = 9,
  NetPr = 10,
  NetCh = 11
};

static constexpr std::array<std::string_view, 12> FCRecoDir = {
  "Reco/Pr/",
  "Reco/APr/",
  "Reco/PiPos/",
  "Reco/PiNeg/",
  "Reco/KaPos/",
  "Reco/KaNeg/",
  "Reco/Pos/",
  "Reco/Neg/",
  "Reco/NetPi/",
  "Reco/NetKa/",
  "Reco/NetPr/",
  "Reco/NetCh/"};

static constexpr std::array<std::string_view, 12> FCGenDir = {
  "Gen/Pr/",
  "Gen/APr/",
  "Gen/PiPos/",
  "Gen/PiNeg/",
  "Gen/KaPos/",
  "Gen/KaNeg/",
  "Gen/Pos/",
  "Gen/Neg/",
  "Gen/NetPi/",
  "Gen/NetKa/",
  "Gen/NetPr/",
  "Gen/NetCh/"};

static constexpr std::array<std::string_view, 8> CountDir = {
  "CountCheck/Pr/",
  "CountCheck/APr/",
  "CountCheck/PiPos/",
  "CountCheck/PiNeg/",
  "CountCheck/KaPos/",
  "CountCheck/KaNeg/",
  "CountCheck/Pos/",
  "CountCheck/Neg/"};

static constexpr std::array<std::string_view, 6> PidDire = {
  "Ch/",
  "Pi/",
  "Ka/",
  "Pr/",
  "El/",
  "De/"};

static constexpr std::array<std::string_view, 2> ChargeDire = {
  "Pos/",
  "Neg/"};

std::string getModifiedStr(const std::string& myString)
{
  size_t pos = myString.rfind('/');
  if (pos != std::string::npos) {
    std::string subString = myString.substr(0, pos); // remove "/" from end of the string
    return subString;
  }
  return myString;
}

enum PidSelectionEnum {
  kPidPi = 0,
  kPidKa,
  kPidPr,
  kPidEl,
  kPidMu,
  kPidDe,
  kNPidSelection
};

enum CutSettingEnum {
  kThrPforTOF = 0,
  kIdCutTypeLowP,
  kNSigmaTPCLowP,
  kNSigmaTOFLowP,
  kNSigmaRadLowP,
  kIdCutTypeHighP,
  kNSigmaTPCHighP,
  kNSigmaTOFHighP,
  kNSigmaRadHighP,
  kDoVetoOthers,
  kDoRelativeTPCcheck,
  kDoRelativeTOFcheck,
  kDoRelativeTPCTOFcheck,
  kNCutSettings
};

enum VetoSettingEnum {
  kDoVetoTPC = 0,
  kDoVetoTOF,
  kVetoTPC,
  kVetoTOF,
  kNVetoSettings
};

static constexpr std::array<std::array<double, kNCutSettings>, kNPidSelection> DefaultPIDcheckValues = {{
  {0.7, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Pi
  {0.5, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Ka
  {0.8, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Pr
  {0.4, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // El
  {0.4, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Mu
  {0.4, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}  // De
}};

static constexpr std::array<std::array<double, kNVetoSettings>, kNPidSelection> DefaultPidVetoValues = {{
  {1, 1, 3.0, 3.0}, // Pi
  {1, 1, 3.0, 3.0}, // Ka
  {1, 1, 3.0, 3.0}, // Pr
  {0, 0, 3.0, 3.0}, // El
  {0, 0, 3.0, 3.0}, // Mu
  {0, 0, 3.0, 3.0}  // De
}};

template <typename T>
inline T getCfg(const auto& cfgLabelledArray, int32_t row, int32_t col)
{
  double val = cfgLabelledArray->get(row, col); // Always stored as double

  if constexpr (std::is_same_v<T, bool>) {
    return val != 0.0; // will return anything nonzero as true
  } else if constexpr (std::is_same_v<T, float>) {
    // cppcheck-suppress suspiciousFloatingPointCast
    return static_cast<float>(val);
  } else if constexpr (std::is_same_v<T, double>) {
    return val;
  } else if constexpr (std::is_same_v<T, int>) {
    return static_cast<int>(val);
  } else {
    static_assert(
      std::is_same_v<T, bool> || std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, int>,
      "getCfg<T>() only supports T = bool, float, double or int");
  }
}

template <int pidMode, typename T>
bool relativeIdOthersTPC(const T& track)
{
  float distSpecies = 1.e6f;

  if constexpr (pidMode == kPidPi) {
    distSpecies = std::fabs(track.tpcNSigmaPi());
  } else if constexpr (pidMode == kPidKa) {
    distSpecies = std::fabs(track.tpcNSigmaKa());
  } else if constexpr (pidMode == kPidPr) {
    distSpecies = std::fabs(track.tpcNSigmaPr());
  } else {
    // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
    // else if constexpr (pidMode == kPidEl) {
    //   distSpecies = std::fabs(track.tpcNSigmaEl());
    // } else if constexpr (pidMode == kPidMu) {
    //   distSpecies = std::fabs(track.tpcNSigmaMu());
    // } else if constexpr (pidMode == kPidDe) {
    //   distSpecies = std::fabs(track.tpcNSigmaDe());
    // }

    return true;
  }

  if constexpr (pidMode != kPidPi) {
    if (std::fabs(track.tpcNSigmaPi()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kPidKa) {
    if (std::fabs(track.tpcNSigmaKa()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kPidPr) {
    if (std::fabs(track.tpcNSigmaPr()) <= distSpecies) {
      return false;
    }
  }

  // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
  // if constexpr (pidMode != kPidEl)
  //   if (std::fabs(track.tpcNSigmaEl()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kPidMu)
  //   if (std::fabs(track.tpcNSigmaMu()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kPidDe)
  //   if (std::fabs(track.tpcNSigmaDe()) <= distSpecies)
  //     return false;

  return true;
}

template <int pidMode, typename T>
bool relativeIdOthersTOF(const T& track)
{
  float distSpecies = 1.e6f;

  if constexpr (pidMode == kPidPi) {
    distSpecies = std::fabs(track.tofNSigmaPi());
  } else if constexpr (pidMode == kPidKa) {
    distSpecies = std::fabs(track.tofNSigmaKa());
  } else if constexpr (pidMode == kPidPr) {
    distSpecies = std::fabs(track.tofNSigmaPr());

    // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
    // } else if constexpr (pidMode == kPidEl) {
    //   distSpecies = std::fabs(track.tofNSigmaEl());
    // } else if constexpr (pidMode == kPidMu) {
    //   distSpecies = std::fabs(track.tofNSigmaMu());
    // } else if constexpr (pidMode == kPidDe) {
    //   distSpecies = std::fabs(track.tofNSigmaDe());
  } else {
    return true;
  }

  if constexpr (pidMode != kPidPi) {
    if (std::fabs(track.tofNSigmaPi()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kPidKa) {
    if (std::fabs(track.tofNSigmaKa()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kPidPr) {
    if (std::fabs(track.tofNSigmaPr()) <= distSpecies) {
      return false;
    }
  }
  // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
  // if constexpr (pidMode != kPidEl)
  //   if (std::fabs(track.tofNSigmaEl()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kPidMu)
  //   if (std::fabs(track.tofNSigmaMu()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kPidDe)
  //   if (std::fabs(track.tofNSigmaDe()) <= distSpecies)
  //     return false;

  return true;
}

template <int pidMode, typename T>
bool relativeIdOthersTPCTOF(const T& track)
{
  float distSpeciesSq = 1.e6f;

  if constexpr (pidMode == kPidPi) {
    distSpeciesSq = track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi();
  } else if constexpr (pidMode == kPidKa) {
    distSpeciesSq = track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa();
  } else if constexpr (pidMode == kPidPr) {
    distSpeciesSq = track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr();
    // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
    // } else if constexpr (pidMode == kPidEl) {
    //   distSpeciesSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
    // } else if constexpr (pidMode == kPidMu) {
    //   distSpeciesSq = track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu();
    // } else if constexpr (pidMode == kPidDe) {
    //   distSpeciesSq = track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe();
  } else {
    return true;
  }

  if constexpr (pidMode != kPidPi) {
    const float distSq = track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi();
    if (distSq <= distSpeciesSq) {
      return false;
    }
  }

  if constexpr (pidMode != kPidKa) {
    const float distSq = track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa();
    if (distSq <= distSpeciesSq) {
      return false;
    }
  }

  if constexpr (pidMode != kPidPr) {
    const float distSq = track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr();
    if (distSq <= distSpeciesSq) {
      return false;
    }
  }
  // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
  // if constexpr (pidMode != kPidEl) {
  //   const float distSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
  //   if (distSq <= distSpeciesSq)
  //     return false;
  // }

  // if constexpr (pidMode != kPidMu) {
  //   const float distSq = track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu();
  //   if (distSq <= distSpeciesSq)
  //     return false;
  // }

  // if constexpr (pidMode != kPidDe) {
  //   const float distSq = track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe();
  //   if (distSq <= distSpeciesSq)
  //     return false;
  // }

  return true;
}

template <int pidMode, typename T>
bool selIdRectangularCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
{
  if constexpr (pidMode == kPidPi) {
    return std::fabs(track.tpcNSigmaPi()) < nSigmaTPC && std::fabs(track.tofNSigmaPi()) < nSigmaTOF;
  } else if constexpr (pidMode == kPidKa) {
    return std::fabs(track.tpcNSigmaKa()) < nSigmaTPC && std::fabs(track.tofNSigmaKa()) < nSigmaTOF;
  } else if constexpr (pidMode == kPidPr) {
    return std::fabs(track.tpcNSigmaPr()) < nSigmaTPC && std::fabs(track.tofNSigmaPr()) < nSigmaTOF;
  } else if constexpr (pidMode == kPidEl) {
    return std::fabs(track.tpcNSigmaEl()) < nSigmaTPC && std::fabs(track.tofNSigmaEl()) < nSigmaTOF;
  } else if constexpr (pidMode == kPidMu) {
    return std::fabs(track.tpcNSigmaMu()) < nSigmaTPC && std::fabs(track.tofNSigmaMu()) < nSigmaTOF;
  } else if constexpr (pidMode == kPidDe) {
    return std::fabs(track.tpcNSigmaDe()) < nSigmaTPC && std::fabs(track.tofNSigmaDe()) < nSigmaTOF;
  } else {
    return false;
  }
}

template <int pidMode, typename T>
bool selIdEllipsoidalCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
{
  if constexpr (pidMode == kPidPi) {
    float tpc = track.tpcNSigmaPi() / nSigmaTPC;
    float tof = track.tofNSigmaPi() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kPidKa) {
    float tpc = track.tpcNSigmaKa() / nSigmaTPC;
    float tof = track.tofNSigmaKa() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kPidPr) {
    float tpc = track.tpcNSigmaPr() / nSigmaTPC;
    float tof = track.tofNSigmaPr() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kPidEl) {
    float tpc = track.tpcNSigmaEl() / nSigmaTPC;
    float tof = track.tofNSigmaEl() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kPidMu) {
    float tpc = track.tpcNSigmaMu() / nSigmaTPC;
    float tof = track.tofNSigmaMu() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kPidDe) {
    float tpc = track.tpcNSigmaDe() / nSigmaTPC;
    float tof = track.tofNSigmaDe() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else {
    return false;
  }
}

// Circular cut : (tpcNSigma^2 + tofNSigma^2) < nSigmaSquaredRad ; for 3 sigma, nSigmaSquaredRad = 9.0
template <int pidMode, typename T>
constexpr bool selIdCircularCut(const T& track, const float& nSigmaSquaredRad)
{
  float tpc = 0.f;
  float tof = 0.f;

  if constexpr (pidMode == kPidPi) {
    tpc = track.tpcNSigmaPi();
    tof = track.tofNSigmaPi();
  } else if constexpr (pidMode == kPidKa) {
    tpc = track.tpcNSigmaKa();
    tof = track.tofNSigmaKa();
  } else if constexpr (pidMode == kPidPr) {
    tpc = track.tpcNSigmaPr();
    tof = track.tofNSigmaPr();
  } else if constexpr (pidMode == kPidEl) {
    tpc = track.tpcNSigmaEl();
    tof = track.tofNSigmaEl();
  } else if constexpr (pidMode == kPidMu) {
    tpc = track.tpcNSigmaMu();
    tof = track.tofNSigmaMu();
  } else if constexpr (pidMode == kPidDe) {
    tpc = track.tpcNSigmaDe();
    tof = track.tofNSigmaDe();
  } else {
    return false; // unknown pidMode
  }
  return (tpc * tpc + tof * tof) < nSigmaSquaredRad;
}

template <int pidMode, typename T>
bool applyVetoOthersTPC(const T& track, const std::array<bool, kNPidSelection>& doVetoTPC, const std::array<float, kNPidSelection>& vetoTPC)
{
  if constexpr (pidMode != kPidPi) {
    if (doVetoTPC[kPidPi] && std::fabs(track.tpcNSigmaPi()) < vetoTPC[kPidPi]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidKa) {
    if (doVetoTPC[kPidKa] && std::fabs(track.tpcNSigmaKa()) < vetoTPC[kPidKa]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidPr) {
    if (doVetoTPC[kPidPr] && std::fabs(track.tpcNSigmaPr()) < vetoTPC[kPidPr]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidEl) {
    if (doVetoTPC[kPidEl] && std::fabs(track.tpcNSigmaEl()) < vetoTPC[kPidEl]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidMu) {
    if (doVetoTPC[kPidMu] && std::fabs(track.tpcNSigmaMu()) < vetoTPC[kPidMu]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidDe) {
    if (doVetoTPC[kPidDe] && std::fabs(track.tpcNSigmaDe()) < vetoTPC[kPidDe]) {
      return false;
    }
  }

  return true;
}

template <int pidMode, typename T>
bool applyVetoOthersTOF(const T& track, const std::array<bool, kNPidSelection>& doVetoTOF, const std::array<float, kNPidSelection>& vetoTOF)
{
  if constexpr (pidMode != kPidPi) {
    if (doVetoTOF[kPidPi] && std::fabs(track.tofNSigmaPi()) < vetoTOF[kPidPi]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidKa) {
    if (doVetoTOF[kPidKa] && std::fabs(track.tofNSigmaKa()) < vetoTOF[kPidKa]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidPr) {
    if (doVetoTOF[kPidPr] && std::fabs(track.tofNSigmaPr()) < vetoTOF[kPidPr]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidEl) {
    if (doVetoTOF[kPidEl] && std::fabs(track.tofNSigmaEl()) < vetoTOF[kPidEl]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidMu) {
    if (doVetoTOF[kPidMu] && std::fabs(track.tofNSigmaMu()) < vetoTOF[kPidMu]) {
      return false;
    }
  }

  if constexpr (pidMode != kPidDe) {
    if (doVetoTOF[kPidDe] && std::fabs(track.tofNSigmaDe()) < vetoTOF[kPidDe]) {
      return false;
    }
  }

  return true;
}

struct NchCumulantsId {

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry genAnalysis{"genAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoAnalysis{"recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry purityAnalysis{"purityAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PDG data base
  // Service<o2::framework::O2DatabasePDG> pdgDB;

  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 10.0, "cut for vertex Z"};
  Configurable<float> cfgCutDcaXY{"cfgCutDcaXY", 0.12, "cut for dcaXY"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.3, "cut for dcaZ"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "cut for eta"};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 3.0, "max cut for pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.15, "min cut for pT"};

  ConfigurableAxis axisNch{"axisNch", {1200, -60.0, 60.0}, "Net_charge_dN"};
  ConfigurableAxis axisPosCh{"axisPosCh", {3010, -0.5, 300.5}, "Pos_charge"};
  ConfigurableAxis axisNegCh{"axisNegCh", {3010, -0.5, 300.5}, "Neg_charge"};
  ConfigurableAxis axisNt{"axisNt", {8010, -0.5, 800.5}, "Mult_midRap_Nch"};
  ConfigurableAxis axisPrCh{"axisPrCh", {3010, -0.5, 300.5}, "Pr_charge"};
  ConfigurableAxis axisAPrCh{"axisAPrCh", {3010, -0.5, 300.5}, "APr_charge"};
  ConfigurableAxis axisKaCh{"axisKaCh", {3010, -0.5, 300.5}, "Ka_charge"};
  ConfigurableAxis axisAKaCh{"axisAKaCh", {3010, -0.5, 300.5}, "AKa_charge"};
  ConfigurableAxis axisPiCh{"axisPiCh", {3010, -0.5, 300.5}, "Pion_Positive"};
  ConfigurableAxis axisAPiCh{"axisAPiCh", {3010, -0.5, 300.5}, "Pion_Negative"};
  // QA axes as ConfigurableAxis
  ConfigurableAxis axisEvents{"axisEvents", {1, 0, 1}, "Counts"};
  ConfigurableAxis axisEta{"axisEta", {100, -1., +1.}, "#eta"};
  ConfigurableAxis axisRapidity{"axisRapidity", {200, -5, 5}, "Rapidity (y)"};
  ConfigurableAxis axisPt{"axisPt", {100, 0., 5.}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisP{"axisP", {100, 0., 5.}, "p (GeV/c)"};
  ConfigurableAxis axisTPCInnerParam{"axisTPCInnerParam", {100, 0, 3}, "P_innerParam_Gev"};
  ConfigurableAxis axisdEdx{"axisdEdx", {100, 20, 500}, "#frac{dE}{dx}"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {80, -20., 20.}, "V_{Z} (cm)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -3., 3.}, "DCA_{Z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -3., 3.}, "DCA_{XY} (cm)"};
  ConfigurableAxis axisMultFT0{"axisMultFT0", {150, 0, 1500}, "MultFT0"};
  ConfigurableAxis axisCent{"axisCent", {103, -1., 102.}, "FT0C(%)"};
  ConfigurableAxis axisPhi{"axisPhi", {80, -1, 7}, "phi"};
  ConfigurableAxis axisTOFBeta{"axisTOFBeta", {40, -2.0, 2.0}, "tofBeta"};
  ConfigurableAxis axisTPCSignal{"axisTPCSignal", {100, -1, 1000}, "tpcSignal"};
  ConfigurableAxis axisTPCNSigma{"axisTPCNSigma", {200, -10.0, 10.0}, "n#sigma_{TPC}"};
  ConfigurableAxis axisTOFNSigma{"axisTOFNSigma", {200, -10.0, 10.0}, "n#sigma_{TOF}"};
  ConfigurableAxis axisTOFExpMom{"axisTOFExpMom", {200, 0.0f, 10.0f}, "#it{p}_{tofExpMom} (GeV/#it{c})"};

  struct : ConfigurableGroup {
    Configurable<bool> checkCollPosZMc{"checkCollPosZMc", false, "checkCollPosZMc"};
    Configurable<bool> flagUnusedVariableError{"flagUnusedVariableError", false, "flagUnusedVariableError"};
    Configurable<bool> cfgDoRejectionForId{"cfgDoRejectionForId", false, "Apply rejection cut before PID selection (selTrackForId)"};
    Configurable<bool> fillSparseForReco{"fillSparseForReco", false, "Fill sparse for reconstructed tracks"};
    Configurable<bool> fillSparseForGen{"fillSparseForGen", false, "Fill sparse for generated tracks"};
    Configurable<bool> fillSparseForPurity{"fillSparseForPurity", false, "Fill sparse for purity tracks"};
    Configurable<bool> cfgEvSel01doNoSameBunchPileup{"cfgEvSel01doNoSameBunchPileup", true, "apply kNoSameBunchPileup"};
    Configurable<bool> cfgEvSel02doIsGoodZvtxFT0vsPV{"cfgEvSel02doIsGoodZvtxFT0vsPV", true, "apply kIsGoodZvtxFT0vsPV"};
    Configurable<bool> cfgEvSel03doIsGoodITSLayersAll{"cfgEvSel03doIsGoodITSLayersAll", true, "apply kIsGoodITSLayersAll"};
    Configurable<bool> cfgDoSubsampling{"cfgDoSubsampling", true, "do subsampling for error estimation"};
    ConfigurableAxis subSampleAxis{"subSampleAxis", {10, 0., 10.}, "Subsample"};
    TRandom3* fRandom = new TRandom3(0); // Random number generator for subsampling
    int currentSubsample = 0;
    Configurable<float> cfgCutPtMinKa{"cfgCutPtMinKa", 0.30, "min pT for kaon eff-correction inclusion"};
    Configurable<float> cfgCutPtMaxKa{"cfgCutPtMaxKa", 1.20, "max pT for kaon eff-correction inclusion (purity cutoff)"};
    Configurable<float> cfgCutPtMinPr{"cfgCutPtMinPr", 0.40, "min pT for proton eff-correction inclusion"};
    Configurable<float> cfgCutPtMaxPr{"cfgCutPtMaxPr", 1.50, "max pT for proton eff-correction inclusion (purity cutoff)"};
    Configurable<float> cfgMaxWeight{"cfgMaxWeight", 20.0f, "max reco weights, use high values to disable the cap"};
    Configurable<bool> cfgGateGenByWeight{"cfgGateGenByWeight", true, "count a gen particle iff that is counted by reco looking at pt eta"};
    Configurable<bool> cfgRequireSameMcColl{"cfgRequireSameMcColl", true, "ignore reco track that belongs to a diff collision"};
  } cfgEventSelection;

  // PID selection: same defaults and configurables as hParticleCorrelation.cxx.
  struct : ConfigurableGroup {
    Configurable<LabeledArray<double>> c1142pidConfigSetting{"c1142pidConfigSetting", {DefaultPIDcheckValues[0].data(), kNPidSelection, kNCutSettings, {"Pi", "Ka", "Pr", "El", "Mu", "De"}, {"ThrPforTOF", "IdCutTypeLowP", "NSigmaTPCLowP", "NSigmaTOFLowP", "NSigmaRadLowP", "IdCutTypeHighP", "NSigmaTPCHighP", "NSigmaTOFHighP", "NSigmaRadHighP", "doVetoOthers", "doRelativeTPCCheck", "doRelativeTOFcheck", "doRelativeTPCTOFcheck"}}, "Cut values for particle identification"};
    Configurable<LabeledArray<double>> c1143pidVetoSetting{"c1143pidVetoSetting", {DefaultPidVetoValues[0].data(), kNPidSelection, kNVetoSettings, {"Pi", "Ka", "Pr", "El", "Mu", "De"}, {"doVetoTPC", "doVetoTOF", "vetoTPC", "vetoTOF"}}, "Veto cuts for particle identification"};
    Configurable<bool> c1144cfgId07CheckTofBeta{"c1144cfgId07CheckTofBeta", false, "Require beta > 0 for reliable TOF"};
  } cfgIdCut;

  // Optional pre-PID electron rejection already present in NchCumulantsId.
  Configurable<float> cfgIdElRejLowNSigma{"cfgIdElRejLowNSigma", -3.0, "cfgIdElRejLowNSigma"};
  Configurable<float> cfgIdElRejHighNSigma{"cfgIdElRejHighNSigma", 5.0, "cfgIdElRejHighNSigma"};
  Configurable<float> cfgIdPiRejNSigma{"cfgIdPiRejNSigma", 3.0, "cfgIdPiRejNSigma"};
  Configurable<float> cfgIdKaRejNSigma{"cfgIdKaRejNSigma", 3.0, "cfgIdKaRejNSigma"};
  Configurable<float> cfgIdPrRejNSigma{"cfgIdPrRejNSigma", 3.0, "cfgIdPrRejNSigma"};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDB01URL{"cfgCCDB01URL", "http://ccdb-test.cern.ch:8080", "cfgCCDB01URL"};
    Configurable<std::string> cfgCCDB02Path{"cfgCCDB02Path", "Users/p/ppanigra/NchCumulants/correctionWeights/", "cfgCCDB02Path"};
    Configurable<int> cfgCCDB03SOR{"cfgCCDB03SOR", 1, "cfgCCDB03SOR"};
    Configurable<std::string> cfgCCDB04StrHistP{"cfgCCDB04StrHistP", "h01_p", "cfgCCDB04StrHistP"};
    Configurable<std::string> cfgCCDB05StrHistPt{"cfgCCDB05StrHistPt", "h02_pt", "cfgCCDB05StrHistPt"};
    Configurable<std::string> cfgCCDB06StrHistPtEta{"cfgCCDB06StrHistPtEta", "h20_pt_eta", "cfgCCDB06StrHistPtEta"};
  } cfgCCDB;

  struct : ConfigurableGroup {
    Configurable<bool> printDebugMessages{"printDebugMessages", false, "printDebugMessages"};
    Configurable<bool> resetHistograms{"resetHistograms", false, "resetHistograms"};
  } cfgDebug;

  // Efficieny containing Histograms.
  TH2F* hPtEtaForBinSearch = nullptr;
  std::vector<std::array<TH2F*, kNeg + 1>> hPtEtaForEffCorrection{kDe + 1, std::array<TH2F*, kNeg + 1>{}};

  // AxisSpec subSampleAxis = {static_cast<int>(cfgNSubsamples), 0., static_cast<double>(cfgNSubsamples), "Subsample"};

  struct EffPowerSums {
    double q1 = 0.;
    double q2 = 0.;
    double q3 = 0.;
    double q4 = 0.;
    // bookeeping for counts check not to be used for cumulants calculation
    double n = 0.;     // track that entered this species (raw, before correction)
    double nDrop = 0.; // track that was dropped due to weight > cfgMaxWeight
    double sumAll = 0; // sum of all weights (including dropped)
  };

  // static constexpr float KmaxWeight = 10.0f; // secondary safeguard for pr and kaon

  inline bool usableWeight(float weight) const
  {
    return weight > 0.f && weight <= cfgEventSelection.cfgMaxWeight.value;
  }

  inline void fillEffPower(EffPowerSums& p, float weight)
  {
    p.n += 1;
    if (weight > 0.f) {
      p.sumAll += weight;
    }
    if (!usableWeight(weight)) {
      p.nDrop += 1;
      return;
    }

    p.q1 += weight;
    p.q2 += weight * weight;
    p.q3 += weight * weight * weight;
    p.q4 += weight * weight * weight * weight;
  }

  // NEW TEMPLATE FUNCTION FOR RECO PROFILES for individual species
  template <FCPrefixEnum Prefix>
  void addFCRecoProfiles()
  {
    constexpr std::string_view Dir = FCRecoDir[static_cast<int>(Prefix)];

    registry.add(std::string(Dir) + "Q1", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q1Sq", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q2", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q1Cube", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q1Q2", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q3", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q1Pow4", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q1SqQ2", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q2Sq", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q1Q3", "", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q4", "", kTProfile, {axisCent});

    if (cfgEventSelection.cfgDoSubsampling) {
      registry.add(std::string(Dir) + "Q1_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q1Sq_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q2_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q1Cube_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q1Q2_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q3_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q1Pow4_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q1SqQ2_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q2Sq_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q1Q3_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q4_subsample", "", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
    }
  }

  template <FCPrefixEnum Prefix>
  void addNetQVectorProfileHistograms()
  {
    constexpr std::string_view Dir = FCRecoDir[static_cast<int>(Prefix)];

    registry.add(std::string(Dir) + "Q_net_1", "Net Q1", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_1Sq", "Net Q1²", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_2", "Net Q2", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_1Cube", "Net Q1³", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_1Q_net_2", "Net Q1·Q2", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_3", "Net Q3", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_1Pow4", "Net Q1⁴", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_1SqQ_net_2", "Net Q1²·Q2", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_2Sq", "Net Q2²", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_1Q_net_3", "Net Q1·Q3", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "Q_net_4", "Net Q4", kTProfile, {axisCent});

    // joint pos-neg correction profiles, they should be stored here. F11 needed for F3 and F21 and F12 for F4 correction
    registry.add(std::string(Dir) + "JointF11", "Joint <N+N->", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "JointF12", "Joint <N+N-(N--1)>", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "JointF21", "Joint <N+(N+-1)N->", kTProfile, {axisCent});

    // Factorial moments
    if (cfgEventSelection.cfgDoSubsampling) {
      registry.add(std::string(Dir) + "Q_net_1_subsample", "Net Q1 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_1Sq_subsample", "Net Q1² Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_2_subsample", "Net Q2 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_1Cube_subsample", "Net Q1³ Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_1Q_net_2_subsample", "Net Q1·Q2 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_3_subsample", "Net Q3 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_1Pow4_subsample", "Net Q1⁴ Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_1SqQ_net_2_subsample", "Net Q1²·Q2 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_2Sq_subsample", "Net Q2² Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_1Q_net_3_subsample", "Net Q1·Q3 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "Q_net_4_subsample", "Net Q4 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});

      registry.add(std::string(Dir) + "JointF11_subsample", "Joint <N+N-> Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "JointF12_subsample", "Joint <N+N-(N--1)> Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "JointF21_subsample", "Joint <N+(N+-1)N-> Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
    }
  }

  // NEW TEMPLATE FUNCTION FOR GEN PROFILES
  template <FCPrefixEnum Prefix>
  void addFCGenProfiles()
  {
    constexpr std::string_view Dir = FCGenDir[static_cast<int>(Prefix)];

    registry.add(std::string(Dir) + "F1", "Gen FactorialMoment1", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "F2", "Gen FactorialMoment2", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "F3", "Gen FactorialMoment3", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "F4", "Gen FactorialMoment4", kTProfile, {axisCent});

    if (cfgEventSelection.cfgDoSubsampling) {
      registry.add(std::string(Dir) + "F1_subsample", "Gen FactorialMoment1 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "F2_subsample", "Gen FactorialMoment2 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "F3_subsample", "Gen FactorialMoment3 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
      registry.add(std::string(Dir) + "F4_subsample", "Gen FactorialMoment4 Subsample", kTProfile2D, {axisCent, cfgEventSelection.subSampleAxis});
    }
  }

  template <FCPrefixEnum Prefix>
  void addCountCheckProfiles()
  {
    constexpr std::string_view Dir = CountDir[static_cast<int>(Prefix)];
    registry.add(std::string(Dir) + "nGenAll", "Generated counted Before Weight Gate", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "nGenUsed", "Generated counted After Weight Gate", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "nRecRaw", "Reco tracks counted Before Weight Gate", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "nRecDropped", "Reco tracks dropped after correction(weight <0 or weight > maxWeight)", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "sumWUsable", "reco corrected: sum of usable weights (= Q1)", kTProfile, {axisCent});
    registry.add(std::string(Dir) + "sumWNoCap", "reco corrected: sum of all positive weights, no cap", kTProfile, {axisCent});
  }

  template <FCPrefixEnum Prefix, typename Registry>
  void fillCountCheck(float genAll, float genUsed, const EffPowerSums& p, float cent, Registry& hReg)
  {
    auto base = HIST(CountDir[static_cast<int>(Prefix)]);
    hReg.fill(base + HIST("nGenAll"), cent, genAll);
    hReg.fill(base + HIST("nGenUsed"), cent, genUsed);
    hReg.fill(base + HIST("nRecRaw"), cent, p.n);
    hReg.fill(base + HIST("nRecDropped"), cent, p.nDrop);
    hReg.fill(base + HIST("sumWUsable"), cent, p.q1);
    hReg.fill(base + HIST("sumWNoCap"), cent, p.sumAll);
  }

  template <FCPrefixEnum Prefix, typename Registry>
  void fillFCBasis(const EffPowerSums& p, float cent, Registry& hReg)
  {
    auto base = HIST(FCRecoDir[static_cast<int>(Prefix)]);

    // Fill Factorial Moments for Reco directly to ensure Reco stores the same type of data as Gen
    // Compute Efficiency Corrected Factorial Moments from power sums
    // These formulas come from elementary symmetric polynomials and their relation to power sums
    hReg.fill(base + HIST("Q1"), cent, p.q1);
    hReg.fill(base + HIST("Q1Sq"), cent, p.q1 * p.q1);
    hReg.fill(base + HIST("Q2"), cent, p.q2);
    hReg.fill(base + HIST("Q1Cube"), cent, p.q1 * p.q1 * p.q1);
    hReg.fill(base + HIST("Q1Q2"), cent, p.q1 * p.q2);
    hReg.fill(base + HIST("Q3"), cent, p.q3);
    hReg.fill(base + HIST("Q1Pow4"), cent, p.q1 * p.q1 * p.q1 * p.q1);
    hReg.fill(base + HIST("Q1SqQ2"), cent, p.q1 * p.q1 * p.q2);
    hReg.fill(base + HIST("Q2Sq"), cent, p.q2 * p.q2);
    hReg.fill(base + HIST("Q1Q3"), cent, p.q1 * p.q3);
    hReg.fill(base + HIST("Q4"), cent, p.q4);
    if (cfgEventSelection.cfgDoSubsampling) {
      hReg.fill(base + HIST("Q1_subsample"), cent, cfgEventSelection.currentSubsample, p.q1);
      hReg.fill(base + HIST("Q1Sq_subsample"), cent, cfgEventSelection.currentSubsample, p.q1 * p.q1);
      hReg.fill(base + HIST("Q2_subsample"), cent, cfgEventSelection.currentSubsample, p.q2);
      hReg.fill(base + HIST("Q1Cube_subsample"), cent, cfgEventSelection.currentSubsample, p.q1 * p.q1 * p.q1);
      hReg.fill(base + HIST("Q1Q2_subsample"), cent, cfgEventSelection.currentSubsample, p.q1 * p.q2);
      hReg.fill(base + HIST("Q3_subsample"), cent, cfgEventSelection.currentSubsample, p.q3);
      hReg.fill(base + HIST("Q1Pow4_subsample"), cent, cfgEventSelection.currentSubsample, p.q1 * p.q1 * p.q1 * p.q1);
      hReg.fill(base + HIST("Q1SqQ2_subsample"), cent, cfgEventSelection.currentSubsample, p.q1 * p.q1 * p.q2);
      hReg.fill(base + HIST("Q2Sq_subsample"), cent, cfgEventSelection.currentSubsample, p.q2 * p.q2);
      hReg.fill(base + HIST("Q1Q3_subsample"), cent, cfgEventSelection.currentSubsample, p.q1 * p.q3);
      hReg.fill(base + HIST("Q4_subsample"), cent, cfgEventSelection.currentSubsample, p.q4);
    }
  }

  // fill function for reco net species
  template <FCPrefixEnum Prefix, typename Registry>
  void fillNetQVectorProfileHistograms(const EffPowerSums& posPow, const EffPowerSums& negPow, float cent, Registry& hReg)
  {
    // Construct base directory path
    auto base = HIST(FCRecoDir[static_cast<int>(Prefix)]);

    // ─── Compute NET q-vectors (sign alternation by power) ───
    // Odd powers: alternate sign (Q_pos - Q_neg)
    // Even powers: same sign (Q_pos + Q_neg)

    double qNet1 = posPow.q1 - negPow.q1; // Odd
    double qNet2 = posPow.q2 + negPow.q2; // Even
    double qNet3 = posPow.q3 - negPow.q3; // Odd
    double qNet4 = posPow.q4 + negPow.q4; // Even

    // Fill joint pos-neg correction profiles
    // capture inline from pos and neg
    double jointF11 = posPow.q1 * negPow.q1;
    double cf2pos = posPow.q1 * posPow.q1 - posPow.q2;
    double cf2neg = negPow.q1 * negPow.q1 - negPow.q2;
    double jointF12 = posPow.q1 * cf2neg;
    double jointF21 = negPow.q1 * cf2pos;
    // ─── Fill raw q-vectors (for cross-verification)
    hReg.fill(base + HIST("Q_net_1"), cent, qNet1);
    hReg.fill(base + HIST("Q_net_1Sq"), cent, qNet1 * qNet1);
    hReg.fill(base + HIST("Q_net_2"), cent, qNet2);
    hReg.fill(base + HIST("Q_net_1Cube"), cent, qNet1 * qNet1 * qNet1);
    hReg.fill(base + HIST("Q_net_1Q_net_2"), cent, qNet1 * qNet2);
    hReg.fill(base + HIST("Q_net_3"), cent, qNet3);
    hReg.fill(base + HIST("Q_net_1Pow4"), cent, qNet1 * qNet1 * qNet1 * qNet1);
    hReg.fill(base + HIST("Q_net_1SqQ_net_2"), cent, qNet1 * qNet1 * qNet2);
    hReg.fill(base + HIST("Q_net_2Sq"), cent, qNet2 * qNet2);
    hReg.fill(base + HIST("Q_net_1Q_net_3"), cent, qNet1 * qNet3);
    hReg.fill(base + HIST("Q_net_4"), cent, qNet4);

    // ─── Fill factorial moments ───
    hReg.fill(base + HIST("JointF11"), cent, jointF11);
    hReg.fill(base + HIST("JointF12"), cent, jointF12);
    hReg.fill(base + HIST("JointF21"), cent, jointF21);

    if (cfgEventSelection.cfgDoSubsampling) {
      hReg.fill(base + HIST("Q_net_1_subsample"), cent, cfgEventSelection.currentSubsample, qNet1);
      hReg.fill(base + HIST("Q_net_1Sq_subsample"), cent, cfgEventSelection.currentSubsample, qNet1 * qNet1);
      hReg.fill(base + HIST("Q_net_2_subsample"), cent, cfgEventSelection.currentSubsample, qNet2);
      hReg.fill(base + HIST("Q_net_1Cube_subsample"), cent, cfgEventSelection.currentSubsample, qNet1 * qNet1 * qNet1);
      hReg.fill(base + HIST("Q_net_1Q_net_2_subsample"), cent, cfgEventSelection.currentSubsample, qNet1 * qNet2);
      hReg.fill(base + HIST("Q_net_3_subsample"), cent, cfgEventSelection.currentSubsample, qNet3);
      hReg.fill(base + HIST("Q_net_1Pow4_subsample"), cent, cfgEventSelection.currentSubsample, qNet1 * qNet1 * qNet1 * qNet1);
      hReg.fill(base + HIST("Q_net_1SqQ_net_2_subsample"), cent, cfgEventSelection.currentSubsample, qNet1 * qNet1 * qNet2);
      hReg.fill(base + HIST("Q_net_2Sq_subsample"), cent, cfgEventSelection.currentSubsample, qNet2 * qNet2);
      hReg.fill(base + HIST("Q_net_1Q_net_3_subsample"), cent, cfgEventSelection.currentSubsample, qNet1 * qNet3);
      hReg.fill(base + HIST("Q_net_4_subsample"), cent, cfgEventSelection.currentSubsample, qNet4);

      hReg.fill(base + HIST("JointF11_subsample"), cent, cfgEventSelection.currentSubsample, jointF11);
      hReg.fill(base + HIST("JointF12_subsample"), cent, cfgEventSelection.currentSubsample, jointF12);
      hReg.fill(base + HIST("JointF21_subsample"), cent, cfgEventSelection.currentSubsample, jointF21);
    }
  }

  template <FCPrefixEnum Prefix, typename Registry>
  void fillGenFactorialMoments(float n, float cent, Registry& hReg)
  {
    auto base = HIST(FCGenDir[static_cast<int>(Prefix)]);

    double f1 = n;
    double f2 = n * (n - 1.);
    double f3 = n * (n - 1.) * (n - 2.);
    double f4 = n * (n - 1.) * (n - 2.) * (n - 3.);

    hReg.fill(base + HIST("F1"), cent, f1);
    hReg.fill(base + HIST("F2"), cent, f2);
    hReg.fill(base + HIST("F3"), cent, f3);
    hReg.fill(base + HIST("F4"), cent, f4);

    if (cfgEventSelection.cfgDoSubsampling) {
      hReg.fill(base + HIST("F1_subsample"), cent, cfgEventSelection.currentSubsample, f1);
      hReg.fill(base + HIST("F2_subsample"), cent, cfgEventSelection.currentSubsample, f2);
      hReg.fill(base + HIST("F3_subsample"), cent, cfgEventSelection.currentSubsample, f3);
      hReg.fill(base + HIST("F4_subsample"), cent, cfgEventSelection.currentSubsample, f4);
    }
  }
  // ─── GEN level (net speciestruth) ───
  template <FCPrefixEnum Prefix, typename Registry>
  void fillGenNetFactorialMoments(float nPos, float nNeg, float cent, Registry& hReg)
  {
    auto base = HIST(FCGenDir[static_cast<int>(Prefix)]);

    double nNet = nPos - nNeg;

    double f1 = nNet;
    double f2 = nNet * (nNet - 1.0);
    double f3 = nNet * (nNet - 1.0) * (nNet - 2.0);
    double f4 = nNet * (nNet - 1.0) * (nNet - 2.0) * (nNet - 3.0);

    hReg.fill(base + HIST("F1"), cent, f1);
    hReg.fill(base + HIST("F2"), cent, f2);
    hReg.fill(base + HIST("F3"), cent, f3);
    hReg.fill(base + HIST("F4"), cent, f4);

    if (cfgEventSelection.cfgDoSubsampling) {
      hReg.fill(base + HIST("F1_subsample"), cent, cfgEventSelection.currentSubsample, f1);
      hReg.fill(base + HIST("F2_subsample"), cent, cfgEventSelection.currentSubsample, f2);
      hReg.fill(base + HIST("F3_subsample"), cent, cfgEventSelection.currentSubsample, f3);
      hReg.fill(base + HIST("F4_subsample"), cent, cfgEventSelection.currentSubsample, f4);
    }
  }

  void init(InitContext const&)
  {
    auto& mgr = o2::ccdb::BasicCCDBManager::instance();
    mgr.setURL(cfgCCDB.cfgCCDB01URL);
    mgr.setCaching(true);
    auto ccdbObj = mgr.getForTimeStamp<TList>(cfgCCDB.cfgCCDB02Path, cfgCCDB.cfgCCDB03SOR);
    if (!ccdbObj) {
      if (cfgDebug.printDebugMessages) {
        LOG(info) << "DEBUG :: CCDB OBJECT NOT FOUND";
      }
    } else {
      if (cfgDebug.printDebugMessages) {
        LOG(info) << "DEBUG :: CCDB OBJECT FOUND";
      }
    }

    ccdbObj->Print();

    hPtEtaForBinSearch = dynamic_cast<TH2F*>(ccdbObj->FindObject("hPtEta"));
    if (cfgDebug.printDebugMessages) {
      LOG(info) << "DEBUG :: Obj Name = " << hPtEtaForBinSearch->GetName() << " :: entries = " << hPtEtaForBinSearch->GetEntries();
    }
    std::string name;
    for (int i = 0; i <= kDe; i++) {
      for (int j = 0; j < (kNeg + 1); j++) {
        name = "hPtEta" + getModifiedStr(static_cast<std::string>(PidDire[i])) + getModifiedStr(static_cast<std::string>(ChargeDire[j]));
        hPtEtaForEffCorrection[i][j] = dynamic_cast<TH2F*>(ccdbObj->FindObject(name.c_str()));
        if (cfgDebug.printDebugMessages) {
          LOG(info) << "DEBUG :: Obj Name = " << hPtEtaForEffCorrection[i][j]->GetName() << " :: entries = " << hPtEtaForBinSearch->GetEntries();
        }
        if (hPtEtaForEffCorrection[i][j]) {
          LOG(info) << "MAPCHECK " << name << " integral=" << hPtEtaForEffCorrection[i][j]->Integral();
        }
      }
    }
    LOG(info) << "MAPCHECK cfgMaxWeight=" << cfgEventSelection.cfgMaxWeight;

    mgr.setURL("http://alice-ccdb.cern.ch"); // RESET the URL otherwise the other process functions which contains ccdb lookups will fail

    const AxisSpec axisIdTag = {32, -0.5f, 31.5f, "idTag"};
    const AxisSpec axisMcTag = {32, -0.5f, 31.5f, "mcTag"};

    HistogramConfigSpec qnHist1({HistType::kTHnSparseD, {axisNch, axisPosCh, axisNegCh, axisPrCh, axisAPrCh, axisKaCh, axisAKaCh, axisNt, axisCent}});
    HistogramConfigSpec qnHist2({HistType::kTHnSparseD, {axisNch, axisPosCh, axisNegCh, axisPiCh, axisAPiCh, axisKaCh, axisAKaCh, axisNt, axisCent}});
    HistogramConfigSpec histTPCPIDSparse({HistType::kTHnSparseD, {axisP, axisTPCNSigma, axisIdTag, axisMcTag}});
    HistogramConfigSpec histTOFPIDSparse({HistType::kTHnSparseD, {axisP, axisTOFNSigma, axisIdTag, axisMcTag}});

    HistogramConfigSpec histPPt({HistType::kTH2F, {axisP, axisPt}});
    HistogramConfigSpec histPTpcInnerParam({HistType::kTH2F, {axisP, axisTPCInnerParam}});
    HistogramConfigSpec histPTpcSignal({HistType::kTH2F, {axisP, axisTPCSignal}});
    HistogramConfigSpec histTpcInnerParamTpcSignal({HistType::kTH2F, {axisTPCInnerParam, axisTPCSignal}});
    HistogramConfigSpec histPBeta({HistType::kTH2F, {axisP, axisTOFBeta}});
    HistogramConfigSpec histTpcInnerParamBeta({HistType::kTH2F, {axisTPCInnerParam, axisTOFBeta}});
    HistogramConfigSpec histPTpcNSigma({HistType::kTH2F, {axisP, axisTPCNSigma}});
    HistogramConfigSpec histPtTpcNSigma({HistType::kTH2F, {axisPt, axisTPCNSigma}});
    HistogramConfigSpec histTpcInnerParamTpcNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTPCNSigma}});
    HistogramConfigSpec histTofExpMomTpcNSigma({HistType::kTH2F, {axisTOFExpMom, axisTPCNSigma}});
    HistogramConfigSpec histPTofNSigma({HistType::kTH2F, {axisP, axisTOFNSigma}});
    HistogramConfigSpec histPtTofNSigma({HistType::kTH2F, {axisPt, axisTOFNSigma}});
    HistogramConfigSpec histTpcInnerParamTofNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTOFNSigma}});
    HistogramConfigSpec histTofExpMomTofNSigma({HistType::kTH2F, {axisTOFExpMom, axisTOFNSigma}});
    HistogramConfigSpec histTpcNSigmaTofNSigma({HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});

    HistogramConfigSpec histPtMc({HistType::kTH1F, {axisPt}});
    // Register histograms for PID validation
    // Register TPC spares per species
    hist.add("PIDValidation/tpcSparse_Pi", "p vs tpcNSigmaPi vs idTag vs mcTag (Pion)", histTPCPIDSparse);
    hist.add("PIDValidation/tpcSparse_Ka", "p vs tpcNSigmaKa vs idTag vs mcTag (Kaon)", histTPCPIDSparse);
    hist.add("PIDValidation/tpcSparse_Pr", "p vs tpcNSigmaPr vs idTag vs mcTag (Proton)", histTPCPIDSparse);

    // Register TOF spares per species
    hist.add("PIDValidation/tofSparse_Pi", "p vs tofNSigmaPi vs idTag vs mcTag (Pion)", histTOFPIDSparse);
    hist.add("PIDValidation/tofSparse_Ka", "p vs tofNSigmaKa vs idTag vs mcTag (Kaon)", histTOFPIDSparse);
    hist.add("PIDValidation/tofSparse_Pr", "p vs tofNSigmaPr vs idTag vs mcTag (Proton)", histTOFPIDSparse);

    // QA check histos

    hist.add("QA/events/preSel/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/events/preSel/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/events/preSel/multFT0", "multFT0", kTH1F, {axisMultFT0});
    hist.add("QA/events/preSel/centFT0", "centFT0", kTH1F, {axisCent});
    hist.addClone("QA/events/preSel/", "QA/events/postSel/");
    hist.add("QA/events/postSel/net_charge", "net_charge", kTH1F, {axisNch});
    hist.add("QA/events/postSel/Nt_centFT", "Mid_rap_Mult_VS_Cent", kTH2D, {{axisCent}, {axisNt}});

    hist.add("QA/tracks/preSel/h_P", "p (Gev/c)", kTH1D, {axisP});
    hist.add("QA/tracks/preSel/h_P_InnerParameter", "p_InnerParameter (Gev/c)", kTH1D, {axisTPCInnerParam});
    hist.add("QA/tracks/preSel/h_Pt", "p_{T} (TPC & TPC+TOF)", kTH1D, {axisPt});
    hist.add("QA/tracks/preSel/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/tracks/preSel/h_phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/tracks/preSel/h2_Pt_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/tracks/preSel/h2_Pt_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/tracks/preSel/h2_p_DcaZ", "DCA_{z}", kTH2D, {{axisP}, {axisDCAz}});
    hist.add("QA/tracks/preSel/h2_p_DcaXY", "DCA_{xy}", kTH2D, {{axisP}, {axisDCAxy}});
    hist.add("QA/tracks/preSel/dE_dx1", "dE/dx vs p", kTH2F, {axisP, axisdEdx});
    hist.add("QA/tracks/preSel/dE_dx2", "dE/dx vs innerparam", kTH2F, {axisTPCInnerParam, axisdEdx});
    hist.add("QA/tracks/preSel/p_pt", "p_vs_pT", kTH2F, {axisP, axisPt});
    hist.add("QA/tracks/preSel/p_pInnerParameter", "p_vs_innerparameter", kTH2F, {axisP, axisTPCInnerParam});

    // tofBeta
    hist.add("QA/tracks/preSel/p_beta", "p_beta", histPBeta);
    hist.add("QA/tracks/preSel/tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);

    // Look at Pion
    hist.add("QA/tracks/preSel/Pi/NoId/p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    hist.addClone("QA/tracks/preSel/Pi/", "QA/tracks/preSel/Ka/");
    hist.addClone("QA/tracks/preSel/Pi/", "QA/tracks/preSel/Pr/");

    hist.addClone("QA/tracks/preSel/", "QA/tracks/postSel/");

    hist.add("QA/tracks/Idfd/Pi/tpcId/p_pt", "p_pt", histPPt);
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    // tpcSignal
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    // tofBeta
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_beta", "p_beta", histPBeta);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    // Look at Pion
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    hist.addClone("QA/tracks/Idfd/Pi/tpcId/", "QA/tracks/Idfd/Pi/tpctofId/");
    hist.addClone("QA/tracks/Idfd/Pi/", "QA/tracks/Idfd/Ka/");
    hist.addClone("QA/tracks/Idfd/Pi/", "QA/tracks/Idfd/Pr/");

    hist.add("sparse1", "sparse1", qnHist1);
    hist.add("sparse2", "sparse2", qnHist2);

    hist.add("sim/gen/sparse1", "sparse1", qnHist1);
    hist.add("sim/gen/sparse2", "sparse2", qnHist2);

    hist.add("sim/reco/sparse1", "sparse1", qnHist1);
    hist.add("sim/reco/sparse2", "sparse2", qnHist2);

    hist.add("sim/purity/sparse1", "sparse1", qnHist1);
    hist.add("sim/purity/sparse2", "sparse2", qnHist2);

    recoAnalysis.add("recoAnalysis/Pi/Pos/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Ka/");
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Pr/");
    // recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/El/");
    // recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/De/");
    recoAnalysis.add("recoAnalysis/Charge/Pos/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    genAnalysis.add("genAnalysis/Pi/Pos/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Pi/Pos/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Pi/Pos/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Pi/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Pi/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Pi/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    genAnalysis.add("genAnalysis/Pi/Neg/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Pi/Neg/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Pi/Neg/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Pi/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Pi/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Pi/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/Ka/");
    genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/Pr/");
    // genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/El/");
    // genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/De/");

    genAnalysis.add("genAnalysis/Charge/Pos/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Charge/Pos/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Charge/Pos/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Charge/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Charge/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Charge/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    genAnalysis.add("genAnalysis/Charge/Neg/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Charge/Neg/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Charge/Neg/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Charge/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Charge/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Charge/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    purityAnalysis.add("purityAnalysis/Charge/Pos/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    purityAnalysis.add("purityAnalysis/Pi/Pos/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/Ka/");
    purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/Pr/");
    // purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/El/");
    // purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/De/");

    // profiles adding for Reco
    addFCRecoProfiles<FCPrefixEnum::Pr>();
    addFCRecoProfiles<FCPrefixEnum::APr>();
    addFCRecoProfiles<FCPrefixEnum::PiPos>();
    addFCRecoProfiles<FCPrefixEnum::PiNeg>();
    addFCRecoProfiles<FCPrefixEnum::KaPos>();
    addFCRecoProfiles<FCPrefixEnum::KaNeg>();
    addFCRecoProfiles<FCPrefixEnum::Pos>();
    addFCRecoProfiles<FCPrefixEnum::Neg>();

    addNetQVectorProfileHistograms<FCPrefixEnum::NetPi>();
    addNetQVectorProfileHistograms<FCPrefixEnum::NetKa>();
    addNetQVectorProfileHistograms<FCPrefixEnum::NetPr>();
    addNetQVectorProfileHistograms<FCPrefixEnum::NetCh>();

    // profiles adding for Gen
    addFCGenProfiles<FCPrefixEnum::Pr>();
    addFCGenProfiles<FCPrefixEnum::APr>();
    addFCGenProfiles<FCPrefixEnum::PiPos>();
    addFCGenProfiles<FCPrefixEnum::PiNeg>();
    addFCGenProfiles<FCPrefixEnum::KaPos>();
    addFCGenProfiles<FCPrefixEnum::KaNeg>();
    addFCGenProfiles<FCPrefixEnum::Pos>();
    addFCGenProfiles<FCPrefixEnum::Neg>();

    addFCGenProfiles<FCPrefixEnum::NetPi>();
    addFCGenProfiles<FCPrefixEnum::NetKa>();
    addFCGenProfiles<FCPrefixEnum::NetPr>();
    addFCGenProfiles<FCPrefixEnum::NetCh>();

    // profiles for counts
    addCountCheckProfiles<FCPrefixEnum::Pr>();
    addCountCheckProfiles<FCPrefixEnum::APr>();
    addCountCheckProfiles<FCPrefixEnum::PiPos>();
    addCountCheckProfiles<FCPrefixEnum::PiNeg>();
    addCountCheckProfiles<FCPrefixEnum::KaPos>();
    addCountCheckProfiles<FCPrefixEnum::KaNeg>();
    addCountCheckProfiles<FCPrefixEnum::Pos>();
    addCountCheckProfiles<FCPrefixEnum::Neg>();

  } // init ends

  static constexpr std::array<std::string_view, 11> HistRegDire2 = {
    "v0Table/Full/",
    "v0Table/postK0sCheck/",
    "v0Table/postMassCut/",
    "v0Table/postSelectionCut/",
    "recoK0s/PreSel/",
    "recoK0s/PostSel/",
    "recoTracks/PreSel/",
    "recoTracks/PostSel/",
    "recoAnalysis/",
    "genAnalysis/",
    "purityAnalysis/"};

  enum HistRegEnum2 {
    v0Full = 0,
    v0PostK0sCheck,
    v0PostMassCut,
    v0PostSelectionCut,
    recoK0sPreSel,
    recoK0sPostSel,
    recoTracksPreSel,
    recoTracksPostSel,
    recoAnalysisDir,
    genAnalysisDir,
    purityAnalysisDir
  };

  enum IdentificationType {
    kTPCidentified = 0,
    kTOFidentified,
    kTPCTOFidentified,
    kUnidentified
  };

  enum TpcTofCutType {
    kCircularCut = 0,
    kRectangularCut,
    kEllipsoidalCut
  };

  enum DetEnum {
    tpcId = 0,
    tofId,
    tpctofId,
    NoId
  };

  static constexpr std::array<std::string_view, 4> DetDire = {
    "tpcId/",
    "tofId/",
    "tpctofId/",
    "NoId/"};

  enum HistRegEnum {
    qaEventPreSel = 0,
    qaEventPostSel,
    qaTracksPreSel,
    qaTracksPostSel,
    qaTracksIdfd,
  };

  static constexpr std::array<std::string_view, 5> HistRegDire = {
    "QA/events/preSel/",
    "QA/events/postSel/",
    "QA/tracks/preSel/",
    "QA/tracks/postSel/",
    "QA/tracks/Idfd/"};

  enum ProcessTypeEnum {
    doDataProcessing = 0,
    doRecoProcessing,
    doPurityProcessing,
    doGenProcessing,
    doSimProcessing
  };

  // Particle identification: ported from hParticleCorrelation.cxx.
  // 0 = circular, 1 = rectangular, 2 = ellipsoidal for combined TPC-TOF cuts.
  template <int pidMode, typename T>
  bool vetoIdOthersTPC(const T& track)
  {
    // Static is only run once, ever.
    static const std::array<bool, kNPidSelection> doVetoTPC = {
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidPi, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidKa, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidPr, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidEl, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidMu, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidDe, kDoVetoTPC)};

    static const std::array<float, kNPidSelection> vetoTPC = {
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidPi, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidKa, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidPr, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidEl, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidMu, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidDe, kVetoTPC)};

    return applyVetoOthersTPC<pidMode>(track, doVetoTPC, vetoTPC);
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTOF(const T& track)
  {
    // Only computed once
    static const std::array<bool, kNPidSelection> doVetoTOF = {
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidPi, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidKa, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidPr, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidEl, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidMu, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPidDe, kDoVetoTOF)};

    static const std::array<float, kNPidSelection> vetoTOF = {
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidPi, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidKa, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidPr, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidEl, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidMu, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPidDe, kVetoTOF)};

    return applyVetoOthersTOF<pidMode>(track, doVetoTOF, vetoTOF);
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTPCTOF(const T& track)
  {
    // If either veto fails, reject
    return vetoIdOthersTPC<pidMode>(track) && vetoIdOthersTOF<pidMode>(track);
  }

  // Check if TOF is reliable
  template <typename T>
  inline bool checkReliableTOF(const T& track)
  {
    // which check makes the information of TOF relaiable? should track.beta() be checked? e.g.:
    if (cfgIdCut.c1144cfgId07CheckTofBeta) {
      return (track.hasTOF() && track.beta() > 0.0f);
    }
    return track.hasTOF();
  }

  template <int pidMode, typename T>
  bool idTPC(const T& track, const float& nSigmaTPC, float& nSigmaIdDistSq)
  {
    if constexpr (pidMode == kPidPi) {
      nSigmaIdDistSq = track.tpcNSigmaPi() * track.tpcNSigmaPi();
    } else if constexpr (pidMode == kPidKa) {
      nSigmaIdDistSq = track.tpcNSigmaKa() * track.tpcNSigmaKa();
    } else if constexpr (pidMode == kPidPr) {
      nSigmaIdDistSq = track.tpcNSigmaPr() * track.tpcNSigmaPr();
    } else if constexpr (pidMode == kPidEl) {
      nSigmaIdDistSq = track.tpcNSigmaEl() * track.tpcNSigmaEl();
    } else if constexpr (pidMode == kPidMu) {
      nSigmaIdDistSq = track.tpcNSigmaMu() * track.tpcNSigmaMu();
    } else if constexpr (pidMode == kPidDe) {
      nSigmaIdDistSq = track.tpcNSigmaDe() * track.tpcNSigmaDe();
    } else {
      nSigmaIdDistSq = 1000000;
    }

    static const bool doVetoOthers = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoVetoOthers);
    if (doVetoOthers) {
      if (!vetoIdOthersTPC<pidMode>(track)) {
        // If vetoIdOthers = true; it passed all veto checks
        // If vetoIdOthers = false; it failed veto check with some other particle
        return false;
      }
    }

    static const bool doRelativeTPCcheck = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoRelativeTPCcheck);
    if (doRelativeTPCcheck) {
      if (!relativeIdOthersTPC<pidMode>(track)) {
        // If relativeIdOthersTPC = true; particle has stronger nSigma compared to others
        // If relativeIdOthersTPC = false; some particle has stronger nSigma compared to it
        return false;
      }
    }

    if constexpr (pidMode == kPidPi) {
      return std::fabs(track.tpcNSigmaPi()) < nSigmaTPC;
    } else if constexpr (pidMode == kPidKa) {
      return std::fabs(track.tpcNSigmaKa()) < nSigmaTPC;
    } else if constexpr (pidMode == kPidPr) {
      return std::fabs(track.tpcNSigmaPr()) < nSigmaTPC;
    } else if constexpr (pidMode == kPidEl) {
      return std::fabs(track.tpcNSigmaEl()) < nSigmaTPC;
    } else if constexpr (pidMode == kPidMu) {
      return std::fabs(track.tpcNSigmaMu()) < nSigmaTPC;
    } else if constexpr (pidMode == kPidDe) {
      return std::fabs(track.tpcNSigmaDe()) < nSigmaTPC;
    } else {
      return false;
    }
  }

  template <int pidMode, typename T>
  bool idTPCTOF(const T& track, const int& pidCutType, const float& nSigmaTPC, const float& nSigmaTOF, const float& nSigmaSquaredRad, float& nSigmaIdDistSq)
  {
    if constexpr (pidMode == kPidPi) {
      nSigmaIdDistSq = track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi();
    } else if constexpr (pidMode == kPidKa) {
      nSigmaIdDistSq = track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa();
    } else if constexpr (pidMode == kPidPr) {
      nSigmaIdDistSq = track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr();
    } else if constexpr (pidMode == kPidEl) {
      nSigmaIdDistSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
    } else if constexpr (pidMode == kPidMu) {
      nSigmaIdDistSq = track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu();
    } else if constexpr (pidMode == kPidDe) {
      nSigmaIdDistSq = track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe();
    } else {
      nSigmaIdDistSq = 1000000;
    }
    static const bool doVetoOthers = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoVetoOthers);
    if (doVetoOthers) {
      if (!vetoIdOthersTPCTOF<pidMode>(track)) {
        // If vetoIdOthers = true; it passed all veto checks
        // If vetoIdOthers = false; it failed veto check with some other particle
        return false;
      }
    }

    static const bool doRelativeTOFcheck = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoRelativeTOFcheck);
    if (doRelativeTOFcheck) {
      if (!relativeIdOthersTOF<pidMode>(track)) {
        // If relativeIdOthersTOF = true; particle has stronger nSigma compared to others
        // If relativeIdOthersTOF = false; some particle has stronger nSigma compared to it
        return false;
      }
    }

    static const bool doRelativeTPCTOFcheck = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoRelativeTPCTOFcheck);
    if (doRelativeTPCTOFcheck) {
      if (!relativeIdOthersTPCTOF<pidMode>(track)) {
        // If relativeIdOthersTPCTOF = true; particle has stronger nSigma compared to others
        // If relativeIdOthersTPCTOF = false; some particle has stronger nSigma compared to it
        return false;
      }
    }

    switch (pidCutType) {
      case kCircularCut:
        return selIdCircularCut<pidMode>(track, nSigmaSquaredRad);
      case kRectangularCut:
        return selIdRectangularCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
      case kEllipsoidalCut:
        return selIdEllipsoidalCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
      default:
        return false;
    }
  }

  template <int pidMode, typename T>
  bool selPdependent(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    // Static cache inside function - initialized once on first call
    static const auto thrPforTOF = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kThrPforTOF);
    static const auto idCutTypeLowP = getCfg<int>(cfgIdCut.c1142pidConfigSetting, pidMode, kIdCutTypeLowP);
    static const auto nSigmaTPCLowP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTPCLowP);
    static const auto nSigmaTOFLowP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTOFLowP);
    static const auto nSigmaRadLowP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaRadLowP);
    static const auto idCutTypeHighP = getCfg<int>(cfgIdCut.c1142pidConfigSetting, pidMode, kIdCutTypeHighP);
    static const auto nSigmaTPCHighP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTPCHighP);
    static const auto nSigmaTOFHighP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTOFHighP);
    static const auto nSigmaRadHighP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaRadHighP);

    if (track.p() < thrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeLowP, nSigmaTPCLowP, nSigmaTOFLowP, nSigmaRadLowP, nSigmaIdDistSq)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
      } else {
        if (idTPC<pidMode>(track, nSigmaTPCLowP, nSigmaIdDistSq)) {
          IdMethod = kTPCidentified;
          return true;
        }
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeHighP, nSigmaTPCHighP, nSigmaTOFHighP, nSigmaRadHighP, nSigmaIdDistSq)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  int getMCTag(const T& track)
  {
    int mcTag = 0;
    if (!track.has_mcParticle()) {
      return mcTag;
    }
    auto mcPart = track.mcParticle();
    int pdgCode = std::abs(mcPart.pdgCode());

    if (pdgCode == kPiPlus || pdgCode == kPiMinus) {
      BITSET(mcTag, McBitPI);
    } else if (pdgCode == kKPlus || pdgCode == kKMinus) {
      BITSET(mcTag, McBitKA);
    } else if (pdgCode == kProton || pdgCode == kProtonBar) {
      BITSET(mcTag, McBitPR);
    } else if (pdgCode == kElectron || pdgCode == kPositron) {
      BITSET(mcTag, McBitEL);
    } else if (pdgCode == kDeuteron || pdgCode == -kDeuteron) {
      BITSET(mcTag, McBitDE);
    }

    return mcTag;
  }
  //
  //______________________________Identification Functions________________________________________________________________
  // Victor slection cuts for electron and other rejections
  template <typename T>
  bool selTrackForId(const T& track)
  {
    return !(cfgIdElRejLowNSigma < track.tpcNSigmaEl() && track.tpcNSigmaEl() < cfgIdElRejHighNSigma &&
             std::fabs(track.tpcNSigmaPi()) > cfgIdPiRejNSigma &&
             std::fabs(track.tpcNSigmaKa()) > cfgIdKaRejNSigma &&
             std::fabs(track.tpcNSigmaPr()) > cfgIdPrRejNSigma);
  }

  // Existing NchCumulantsId two-argument PID interface (data / MC / purity).
  template <typename T>
  bool selPion(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selPdependent<kPidPi>(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selKaon(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selPdependent<kPidKa>(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selProton(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selPdependent<kPidPr>(track, IdMethod, nSigmaIdDistSq);
  }

  template <int Mode, int pidMode, int detMode, typename H, typename T>
  void fillIdentificationQA(H histReg, const T& track)
  {
    float tpcNSigmaVal = -999, tofNSigmaVal = -999;
    switch (pidMode) {
      case kPi:
        tpcNSigmaVal = track.tpcNSigmaPi();
        tofNSigmaVal = track.tofNSigmaPi();
        break;
      case kKa:
        tpcNSigmaVal = track.tpcNSigmaKa();
        tofNSigmaVal = track.tofNSigmaKa();
        break;
      case kPr:
        tpcNSigmaVal = track.tpcNSigmaPr();
        tofNSigmaVal = track.tofNSigmaPr();
        break;
      case kEl:
        tpcNSigmaVal = track.tpcNSigmaEl();
        tofNSigmaVal = track.tofNSigmaEl();
        break;
      case kDe:
        tpcNSigmaVal = track.tpcNSigmaDe();
        tofNSigmaVal = track.tofNSigmaDe();
        break;
      default:
        tpcNSigmaVal = -999, tofNSigmaVal = -999;
        break;
    }

    // NSigma
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("p_tpcNSigma"), track.p(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("pt_tpcNSigma"), track.pt(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tpcInnerParam_tpcNSigma"), track.tpcInnerParam(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tofExpMom_tpcNSigma"), track.tofExpMom(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("p_tofNSigma"), track.p(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("pt_tofNSigma"), track.pt(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tpcInnerParam_tofNSigma"), track.tpcInnerParam(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tofExpMom_tofNSigma"), track.tofExpMom(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tpcNSigma_tofNSigma"), tpcNSigmaVal, tofNSigmaVal);
  }

  template <int mode, typename T>
  void fillCollQA(const T& coll, const float& nCh, const float& nT)
  {
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_VtxZ"), coll.posZ());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_Counts"), 0.5);
    hist.fill(HIST(HistRegDire[mode]) + HIST("multFT0"), coll.multFT0C());
    hist.fill(HIST(HistRegDire[mode]) + HIST("centFT0"), coll.centFT0M());
    if (mode == qaEventPostSel) {
      hist.fill(HIST(HistRegDire[mode]) + HIST("net_charge"), nCh);
      hist.fill(HIST(HistRegDire[mode]) + HIST("Nt_centFT"), coll.centFT0M(), nT);
    }
  }

  template <int mode, typename T>
  void fillTrackQA(const T& track)
  {
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_P"), track.p());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_P_InnerParameter"), track.tpcInnerParam());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_Pt"), track.pt());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_Eta"), track.eta());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_phi"), track.phi());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_Pt_DcaZ"), track.pt(), track.dcaZ());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_Pt_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_p_DcaZ"), track.p(), track.dcaZ());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_p_DcaXY"), track.p(), track.dcaXY());
    hist.fill(HIST(HistRegDire[mode]) + HIST("dE_dx1"), track.p(), track.tpcSignal());
    hist.fill(HIST(HistRegDire[mode]) + HIST("dE_dx2"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST(HistRegDire[mode]) + HIST("p_pt"), track.p(), track.pt());
    hist.fill(HIST(HistRegDire[mode]) + HIST("p_pInnerParameter"), track.p(), track.tpcInnerParam());

    // tofBeta
    hist.fill(HIST(HistRegDire[mode]) + HIST("p_beta"), track.p(), track.beta());
    hist.fill(HIST(HistRegDire[mode]) + HIST("tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());

    // Look at Pion
    fillIdentificationQA<mode, kPi, NoId>(hist, track); // Look at Pion
    fillIdentificationQA<mode, kKa, NoId>(hist, track); // Look at Kaon
    fillIdentificationQA<mode, kPr, NoId>(hist, track); // Look at Proton
  }

  //  template <int Mode, int pidMode, typename H, typename T>
  //  void fillGenTrackQA(H& histReg, const T& mcTrack)
  //  {
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h12_p"), mcTrack.p());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h13_pt"), mcTrack.pt());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h14_eta"), mcTrack.eta());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h15_phi"), mcTrack.phi());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h16_rapidity"), mcTrack.y());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
  //  }
  //
  //  template <int Mode, int pidMode, typename H, typename T>
  //  void fillRecoTrackQA(H& histReg, const T& track)
  //  {
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h12_p"), track.p());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h13_pt"), track.pt());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h14_eta"), track.eta());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h15_phi"), track.phi());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h16_rapidity"), track.y());
  //    histReg.fill(HIST(HistRegDire2[Mode]) + HIST(PidDire[pidMode]) + HIST("h20_pt_eta"), track.pt(), track.eta());
  //  }

  template <HistRegEnum2 DIR, PidEnum P, ChargeEnum C, typename TrackType>
  void fillGenTrackQA(HistogramRegistry& r, const TrackType& t)
  {
    auto base =
      HIST(HistRegDire2[DIR]) +
      HIST(PidDire[P]) +
      HIST(ChargeDire[C]);
    r.fill(base + HIST("h12_p"), t.p());
    r.fill(base + HIST("h13_pt"), t.pt());
    r.fill(base + HIST("h14_eta"), t.eta());
    r.fill(base + HIST("h15_phi"), t.phi());
    r.fill(base + HIST("h16_rapidity"), t.y());
    r.fill(base + HIST("h20_pt_eta"), t.pt(), t.eta());
  }

  template <HistRegEnum2 DIR, PidEnum P, ChargeEnum C, typename TrackType>
  void fillRecoTrackQA(HistogramRegistry& r, const TrackType& t)
  {
    auto base =
      HIST(HistRegDire2[DIR]) +
      HIST(PidDire[P]) +
      HIST(ChargeDire[C]);
    r.fill(base + HIST("h12_p"), t.p());
    r.fill(base + HIST("h13_pt"), t.pt());
    r.fill(base + HIST("h14_eta"), t.eta());
    r.fill(base + HIST("h15_phi"), t.phi());
    r.fill(base + HIST("h16_rapidity"), t.y());
    r.fill(base + HIST("h20_pt_eta"), t.pt(), t.eta());
  }

  template <HistRegEnum2 DIR, PidEnum P, ChargeEnum C, typename TrackType>
  void fillPurityTrackQA(HistogramRegistry& r, const TrackType& t)
  {
    auto base =
      HIST(HistRegDire2[DIR]) +
      HIST(PidDire[P]) +
      HIST(ChargeDire[C]);
    r.fill(base + HIST("h12_p"), t.p());
    r.fill(base + HIST("h13_pt"), t.pt());
    r.fill(base + HIST("h14_eta"), t.eta());
    r.fill(base + HIST("h15_phi"), t.phi());
    r.fill(base + HIST("h16_rapidity"), t.y());
    r.fill(base + HIST("h20_pt_eta"), t.pt(), t.eta());
  }

  template <typename T, typename H>
  void executeTrackAnalysisPart(const T& track, const int& trackIdTag, float& nP, float& nM,
                                const int& idMethodPi, const bool& trackIsPion, float& nAPi, float& nPi,
                                const int& idMethodKa, const bool& trackIsKaon, float& nAKa, float& nKa,
                                const int& idMethodPr, const bool& trackIsProton, float& nPr, float& nAPr, H& hReg)
  {
    if (cfgEventSelection.flagUnusedVariableError) {
      LOG(info) << trackIdTag << idMethodPi << ":" << idMethodKa << ":" << idMethodPr;
    }
    if (track.sign() > 0) {
      // fillRecoTrackQA<kPos>(hReg, track);
      nP++;
      hReg.fill(HIST("recoAnalysis/Charge/Pos/h12_p"), track.p());
      hReg.fill(HIST("recoAnalysis/Charge/Pos/h13_pt"), track.pt());
      hReg.fill(HIST("recoAnalysis/Charge/Pos/h14_eta"), track.eta());
      hReg.fill(HIST("recoAnalysis/Charge/Pos/h15_phi"), track.phi());
      hReg.fill(HIST("recoAnalysis/Charge/Pos/h16_rapidity"), track.y());
      hReg.fill(HIST("recoAnalysis/Charge/Pos/h20_pt_eta"), track.pt(), track.eta());
    }
    if (track.sign() < 0) {
      // fillRecoTrackQA<kNeg>(hReg, track);
      nM++;
      hReg.fill(HIST("recoAnalysis/Charge/Neg/h12_p"), track.p());
      hReg.fill(HIST("recoAnalysis/Charge/Neg/h13_pt"), track.pt());
      hReg.fill(HIST("recoAnalysis/Charge/Neg/h14_eta"), track.eta());
      hReg.fill(HIST("recoAnalysis/Charge/Neg/h15_phi"), track.phi());
      hReg.fill(HIST("recoAnalysis/Charge/Neg/h16_rapidity"), track.y());
      hReg.fill(HIST("recoAnalysis/Charge/Neg/h20_pt_eta"), track.pt(), track.eta());
    }

    if (trackIsPion) {
      // if (idMethodPi == kTPCidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPi, tpcId, true>(hist, track); // set hist as hReg after tpcId etc add true
      // } else if (idMethodPi == kTPCTOFidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPi, tpctofId, true>(hist, track);
      // } else if (idMethodPi == kUnidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPi, NoId, true>(hist, track);
      // }
      if (track.sign() > 0) {
        nPi++;
        fillRecoTrackQA<recoAnalysisDir, kPi, kPos>(hReg, track);
      } else if (track.sign() < 0) {
        nAPi++;
        fillRecoTrackQA<recoAnalysisDir, kPi, kNeg>(hReg, track);
      }
      // fillRecoTrackQA<recoAnalysisDir, kPi>(hReg, track);
    }
    if (trackIsKaon) {
      // if (idMethodKa == kTPCidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kKa, tpcId, true>(hist, track);
      // } else if (idMethodKa == kTPCTOFidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kKa, tpctofId, true>(hist, track);
      // } else if (idMethodKa == kUnidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kKa, NoId, true>(hist, track);
      // }
      if (track.sign() > 0) {
        nKa++;
        fillRecoTrackQA<recoAnalysisDir, kKa, kPos>(hReg, track);
      } else if (track.sign() < 0) {
        nAKa++;
        fillRecoTrackQA<recoAnalysisDir, kKa, kNeg>(hReg, track);
      }
      // fillRecoTrackQA<recoAnalysisDir, kKa>(hReg, track);
    }
    if (trackIsProton) {
      // if (idMethodPr == kTPCidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPr, tpcId, true>(hist, track);
      // } else if (idMethodPr == kTPCTOFidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPr, tpctofId, true>(hist, track);
      // } else if (idMethodPr == kUnidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPr, NoId, true>(hist, track);
      // }
      if (track.sign() > 0) {
        nPr++;
        fillRecoTrackQA<recoAnalysisDir, kPr, kPos>(hReg, track);
      } else if (track.sign() < 0) {
        nAPr++;
        fillRecoTrackQA<recoAnalysisDir, kPr, kNeg>(hReg, track);
      }
      // fillRecoTrackQA<recoAnalysisDir, kPr>(hReg, track);
    }

    // hReg.fill(HIST("recoAnalysis/SelectedTrack_IdentificationTag"), trackIdTag);
  }
  // fill the basis for factorial cumulants

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                                aod::TracksDCA, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl,
                                aod::pidTPCFullMu, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullDe>;
  using MyCollisions = soa::Join<aod::Collisions, aod::Mults, aod::CentFT0Ms, aod::EvSels>;

  using MyTracksWithMclabels = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTOFmass,
                                         aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullDe,
                                         aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullDe, aod::McTrackLabels>;
  using MyCollisionsWithMcLabels = soa::Join<aod::Collisions, aod::Mults, aod::CentFT0Ms, aod::EvSels, aod::McCollisionLabels>;

  template <typename CollisionType>
  bool isEventSelected(const CollisionType& coll)
  {
    if (cfgEventSelection.cfgEvSel01doNoSameBunchPileup &&
        !coll.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (cfgEventSelection.cfgEvSel02doIsGoodZvtxFT0vsPV &&
        !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (cfgEventSelection.cfgEvSel03doIsGoodITSLayersAll &&
        !coll.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    return true;
  }

  // tracks and collision filters
  Filter colSel8 = aod::evsel::sel8 == true;
  Filter colFilter = nabs(aod::collision::posZ) < cfgCutPosZ;
  Filter trackFilter = requireGlobalTrackInFilter();
  Filter trackPt = (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax);
  Filter trackDCAxy = nabs(aod::track::dcaXY) < cfgCutDcaXY;
  Filter trackDCAz = nabs(aod::track::dcaZ) < cfgCutDcaZ;
  Filter tracketa = nabs(aod::track::eta) < cfgCutEta;

  using MyFilteredCol = soa::Filtered<MyCollisions>;
  using MyFilteredTracks = soa::Filtered<MyAllTracks>;

  using MyFilteredColsWithMcLabels = soa::Filtered<MyCollisionsWithMcLabels>;
  using MyFilteredTracksWithMcLabels = soa::Filtered<MyTracksWithMclabels>;

  // manual sliceby
  SliceCache cache;

  Preslice<MyFilteredTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<MyFilteredTracksWithMcLabels> mctracksPerCollisionPreslice = o2::aod::track::collisionId;

  template <int analysisType, typename C, typename T>
  void executeAnalysis(const C& collisions, const T& tracks)
  {
    float nP = 0;
    float nM = 0;
    float nCh = 0;
    float nT = 0;
    float nPr = 0;
    float nAPr = 0;
    float nKa = 0;
    float nAKa = 0;
    float nPi = 0;
    float nAPi = 0;

    bool trackIsPion = false;
    bool trackIsKaon = false;
    bool trackIsProton = false;

    int trackIdTag = 0;
    int idMethodPi = kUnidentified;
    int idMethodKa = kUnidentified;
    int idMethodPr = kUnidentified;

    int ptEtaBin = -1;

    if constexpr (analysisType == doDataProcessing) {
      for (const auto& col : collisions) {
        if (!isEventSelected(col)) {
          continue;
        }
        nP = 0;
        nM = 0;
        nCh = 0;
        nT = 0;
        nPr = 0;
        nAPr = 0;
        nKa = 0;
        nAKa = 0;
        nPi = 0;
        nAPi = 0;

        // group tracks manually with corresponding collision using col id;
        const uint64_t collIdx = col.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);

        for (const auto& track : tracksTablePerColl) {

          ptEtaBin = hPtEtaForBinSearch->FindBin(track.pt(), track.eta()); // Find Track Bin for efficiency correction

          fillTrackQA<qaTracksPostSel>(track);

          if (track.sign() == 1) {
            // nP++;
            nP += hPtEtaForEffCorrection[kCh][kPos]->GetBinContent(ptEtaBin);
          }
          if (track.sign() == -1) {
            nM += hPtEtaForEffCorrection[kCh][kNeg]->GetBinContent(ptEtaBin);
          }

          // Reject electrons first
          if (cfgEventSelection.cfgDoRejectionForId && !selTrackForId(track)) {
            continue;
          }

          int idMethod = -1;
          // pion
          if (selPion(track, idMethod)) {
            if (track.sign() == 1) {
              nPi += hPtEtaForEffCorrection[kPi][kPos]->GetBinContent(ptEtaBin);
            }
            if (track.sign() == -1) {
              nAPi += hPtEtaForEffCorrection[kPi][kNeg]->GetBinContent(ptEtaBin);
            }

            if (idMethod == kTPCidentified) {
              fillIdentificationQA<qaTracksIdfd, kPi, tpcId>(hist, track);
            }
            if (idMethod == kTPCTOFidentified) {
              fillIdentificationQA<qaTracksIdfd, kPi, tpctofId>(hist, track);
            }
          }
          // kaon
          if (selKaon(track, idMethod)) {
            if (track.sign() == 1) {
              nKa += hPtEtaForEffCorrection[kKa][kPos]->GetBinContent(ptEtaBin);
            }
            if (track.sign() == -1) {
              nAKa += hPtEtaForEffCorrection[kKa][kNeg]->GetBinContent(ptEtaBin);
            }

            if (idMethod == kTPCidentified) {
              fillIdentificationQA<qaTracksIdfd, kKa, tpcId>(hist, track);
            }
            if (idMethod == kTPCTOFidentified) {
              fillIdentificationQA<qaTracksIdfd, kKa, tpctofId>(hist, track);
            }
          }
          // proton
          if (selProton(track, idMethod)) {
            if (track.sign() == 1) {
              nPr += hPtEtaForEffCorrection[kPr][kPos]->GetBinContent(ptEtaBin);
            }
            if (track.sign() == -1) {
              nAPr += hPtEtaForEffCorrection[kPr][kNeg]->GetBinContent(ptEtaBin);
            }

            if (idMethod == kTPCidentified) {
              fillIdentificationQA<qaTracksIdfd, kPr, tpcId>(hist, track);
            }
            if (idMethod == kTPCTOFidentified) {
              fillIdentificationQA<qaTracksIdfd, kPr, tpctofId>(hist, track);
            }
          }
        } // track loop ends
        nCh = nP - nM;
        nT = nP + nM;

        fillCollQA<qaEventPostSel>(col, nCh, nT);

        hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
        hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());

      } // collision loop ends
    } else if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

      for (const auto& col : collisions) {

        if (!col.has_mcCollision()) {
          LOG(warning) << "No MC collision for this collision, skip...";
          continue;
        }
        if (!isEventSelected(col)) {
          continue;
        }
        nP = 0;
        nM = 0;
        nCh = 0;
        nT = 0;
        nPr = 0;
        nAPr = 0;
        nKa = 0;
        nAKa = 0;
        nPi = 0;
        nAPi = 0;
        // group tracks manually with corresponding collision using col id;
        const uint64_t collIdx = col.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(mctracksPerCollisionPreslice, collIdx);

        for (const auto& track : tracksTablePerColl) {

          if (!track.has_mcParticle()) {
            LOG(warning) << "No MC Particle for this track, skip...";
            continue;
          }

          auto mcPart = track.mcParticle();

          fillTrackQA<qaTracksPostSel>(track);

          // Do Proper Track Identification
          trackIsPion = false;
          trackIsKaon = false;
          trackIsProton = false;

          trackIdTag = 0;
          idMethodPi = kUnidentified;
          idMethodKa = kUnidentified;
          idMethodPr = kUnidentified;

          // Get MC tag
          int mcTag = getMCTag(track);

          // Reject electrons first
          if (cfgEventSelection.cfgDoRejectionForId && !selTrackForId(track)) {
            continue;
          }

          if (selPion(track, idMethodPi)) {
            trackIsPion = true;
            BITSET(trackIdTag, IdBitPI);
            hist.fill(HIST("PIDValidation/tpcSparse_Pi"), track.p(), track.tpcNSigmaPi(), trackIdTag, mcTag);
            if (track.hasTOF()) {
              hist.fill(HIST("PIDValidation/tofSparse_Pi"), track.p(), track.tofNSigmaPi(), trackIdTag, mcTag);
            }
          }
          if (selKaon(track, idMethodKa)) {
            trackIsKaon = true;
            BITSET(trackIdTag, IdBitKA);
            hist.fill(HIST("PIDValidation/tpcSparse_Ka"), track.p(), track.tpcNSigmaKa(), trackIdTag, mcTag);
            if (track.hasTOF()) {
              hist.fill(HIST("PIDValidation/tofSparse_Ka"), track.p(), track.tofNSigmaKa(), trackIdTag, mcTag);
            }
          }
          if (selProton(track, idMethodPr)) {
            trackIsProton = true;
            BITSET(trackIdTag, IdBitPR);
            hist.fill(HIST("PIDValidation/tpcSparse_Pr"), track.p(), track.tpcNSigmaPr(), trackIdTag, mcTag);
            if (track.hasTOF()) {
              hist.fill(HIST("PIDValidation/tofSparse_Pr"), track.p(), track.tofNSigmaPr(), trackIdTag, mcTag);
            }
          }

          if constexpr (analysisType == doPurityProcessing) {
            if (trackIsPion) {
              if (track.sign() > 0 && mcPart.pdgCode() != kPiPlus) {
                trackIsPion = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kPiMinus) {
                trackIsPion = false;
              }
            }
            if (trackIsKaon) {
              if (track.sign() > 0 && mcPart.pdgCode() != kKPlus) {
                trackIsKaon = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kKMinus) {
                trackIsKaon = false;
              }
            }
            if (trackIsProton) {
              if (track.sign() > 0 && mcPart.pdgCode() != kProton) {
                trackIsProton = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kProtonBar) {
                trackIsProton = false;
              }
            }
          }

          executeTrackAnalysisPart(track, trackIdTag, nP, nM,
                                   idMethodPi, trackIsPion, nAPi, nPi,
                                   idMethodKa, trackIsKaon, nAKa, nKa,
                                   idMethodPr, trackIsProton, nPr, nAPr, recoAnalysis);

        } // track itteration ends

        nCh = nP - nM;
        nT = nP + nM;

        fillCollQA<qaEventPostSel>(col, nCh, nT);

        hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
        hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());

      } // collision ends
    }
  }

  void processRun3(MyFilteredCol const& collisions, MyFilteredTracks const& tracks)
  {
    executeAnalysis<doDataProcessing>(collisions, tracks);
  }
  PROCESS_SWITCH(NchCumulantsId, processRun3, "Process for Run-3", false);

  void processMcRecco(MyFilteredColsWithMcLabels const& collisions, MyFilteredTracksWithMcLabels const& tracks, aod::McParticles const&)
  {
    executeAnalysis<doRecoProcessing>(collisions, tracks);
  }
  PROCESS_SWITCH(NchCumulantsId, processMcRecco, "Process for MC Recco", false);

  void processPurity(MyFilteredColsWithMcLabels const& collisions, MyFilteredTracksWithMcLabels const& tracks, aod::McParticles const&)
  {
    executeAnalysis<doPurityProcessing>(collisions, tracks);
  }
  PROCESS_SWITCH(NchCumulantsId, processPurity, "Process for MC generated purity", false);
  // process fun is over

  // process for truth and response matrix starts now

  Preslice<aod::McParticles> mcTracksPerMcCollisionPreslice = o2::aod::mcparticle::mcCollisionId;

  using MyMcCollisions = aod::McCollisions;

  void processGen(MyMcCollisions const&, MyFilteredColsWithMcLabels const& collisions, aod::McParticles const& mcParticles)
  {
    for (const auto& col : collisions) {
      if (!col.has_mcCollision()) {
        LOG(warning) << "No MC collision for this event, skip...";
        continue;
      }
      if (!isEventSelected(col)) {
        continue;
      }
      const auto& mcColl = col.mcCollision();

      // ---- apply same Vz cut as data/reco ----
      if (std::abs(mcColl.posZ()) > cfgCutPosZ) {
        continue; // reject GEN events outside acceptance
      }
      // slice MCParticles for this MC collision
      const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcColl.globalIndex());

      float nP = 0;
      float nM = 0;
      float nCh = 0;
      float nT = 0;
      float nPr = 0;
      float nAPr = 0;
      float nKa = 0;
      float nAKa = 0;
      float nPi = 0;
      float nAPi = 0;
      for (const auto& mcTrack : mcTracksTablePerMcColl) {
        if (!mcTrack.isPhysicalPrimary()) {
          continue;
        }
        // pt / eta acceptance (same as RECO acceptance)
        if (mcTrack.pt() <= cfgCutPtMin ||
            mcTrack.pt() >= cfgCutPtMax ||
            std::abs(mcTrack.eta()) >= cfgCutEta) {
          continue;
        }

        int pdg = mcTrack.pdgCode();

        if (pdg == kPiPlus || pdg == kKPlus || pdg == kProton || pdg == kPositron || pdg == kMuonPlus || pdg == kDeuteron) {
          // fillGenTrackQA<genAnalysisDir, kPos>(genAnalysis, mcTrack);
          nP++;

          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
        } else if (pdg == kPiMinus || pdg == kKMinus || pdg == kProtonBar || pdg == kElectron || pdg == kMuonMinus || pdg == -kDeuteron) {
          // fillGenTrackQA<genAnalysisDir, kNeg>(genAnalysis, mcTrack);
          nM++;

          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
        }
        // ----- Pions -----
        if (pdg == kPiPlus) {
          // fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          nPi++;
          fillGenTrackQA<genAnalysisDir, kPi, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kPiMinus) {
          // fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          nAPi++;
          fillGenTrackQA<genAnalysisDir, kPi, kNeg>(genAnalysis, mcTrack);
        } else if (pdg == kKPlus) {
          // fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          nKa++;
          fillGenTrackQA<genAnalysisDir, kKa, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kKMinus) {
          // fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          nAKa++;
          fillGenTrackQA<genAnalysisDir, kKa, kNeg>(genAnalysis, mcTrack);
        } else if (pdg == kProton) {
          // fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          nPr++;
          fillGenTrackQA<genAnalysisDir, kPr, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kProtonBar) {
          // fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          nAPr++;
          fillGenTrackQA<genAnalysisDir, kPr, kNeg>(genAnalysis, mcTrack);
        }
      }
      nCh = nP - nM;
      nT = nP + nM;

      hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
      hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());
    } // collision ends
  }
  PROCESS_SWITCH(NchCumulantsId, processGen, "Process for MC generated or truth", false);

  void processSim(MyFilteredColsWithMcLabels const& collisions, MyFilteredTracksWithMcLabels const& tracks, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {

    if (cfgEventSelection.flagUnusedVariableError) {
      LOG(info) << mcCollisions.size();
    }
    bool trackIsPion = false;
    bool trackIsKaon = false;
    bool trackIsProton = false;

    int trackIdTag = 0;
    int idMethodPi = kUnidentified;
    int idMethodKa = kUnidentified;
    int idMethodPr = kUnidentified;

    int ptEtaBin = -1;

    for (auto const& col : collisions) {

      if (!col.has_mcCollision()) {
        LOG(warning) << "No MC collision for this collision, skip...";
        continue;
      }
      if (!isEventSelected(col)) {
        continue;
      }
      auto mcCollision = col.mcCollision();

      if (cfgEventSelection.checkCollPosZMc && std::abs(mcCollision.posZ()) > cfgCutPosZ) {
        continue;
      }

      // slice reco tracks to this collision
      const uint64_t collIdx = col.globalIndex();
      const auto tracksTablePerColl = tracks.sliceBy(mctracksPerCollisionPreslice, collIdx);

      // slice mc particles to mc collisions
      const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcCollision.globalIndex());
      float cent = col.centFT0M();

      EffPowerSums prPow;
      EffPowerSums aprPow;
      EffPowerSums pipPow;
      EffPowerSums pimPow;
      EffPowerSums kapPow;
      EffPowerSums kamPow;
      EffPowerSums posPow;
      EffPowerSums negPow;

      cfgEventSelection.currentSubsample = static_cast<int>(cfgEventSelection.fRandom->Uniform(0, cfgEventSelection.subSampleAxis.value[0]));

      // Denominator -- Generator level(truth)

      float nPGen = 0, nMGen = 0, nChGen = 0, nTGen = 0;
      float nPrGen = 0, nAPrGen = 0;
      float nKaGen = 0, nAKaGen = 0;
      float nPiGen = 0, nAPiGen = 0;

      std::array<float, 8> genAll{}; // candidates per species (index = FCPrefixEnum 0..7), before the weight gate
      auto okCell = [&](FCPrefixEnum sp, int pid, int sgn, double pt, double eta) {
        genAll[static_cast<int>(sp)] += 1.f;
        if (!cfgEventSelection.cfgGateGenByWeight) {
          return true;
        }
        const int bin = hPtEtaForBinSearch->FindBin(pt, eta);
        return usableWeight(hPtEtaForEffCorrection[pid][sgn]->GetBinContent(bin));
      };

      for (const auto& mcTrack : mcTracksTablePerMcColl) {
        if (!mcTrack.isPhysicalPrimary()) {
          continue;
        }
        if (mcTrack.pt() <= cfgCutPtMin ||
            mcTrack.pt() >= cfgCutPtMax ||
            std::abs(mcTrack.eta()) >= cfgCutEta) {
          continue;
        }

        int pdg = mcTrack.pdgCode();

        if (pdg == kPiPlus || pdg == kKPlus || pdg == kProton || pdg == kPositron || pdg == kMuonPlus || pdg == kDeuteron) {
          // fillGenTrackQA<genAnalysisDir, kPos>(genAnalysis, mcTrack);
          // nPGen++;
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());

          if (okCell(FCPrefixEnum::Pos, kCh, kPos, mcTrack.pt(), mcTrack.eta())) {
            nPGen++;
          }
        } else if (pdg == kPiMinus || pdg == kKMinus || pdg == kProtonBar || pdg == kElectron || pdg == kMuonMinus || pdg == -kDeuteron) {
          // fillGenTrackQA<genAnalysisDir, kNeg>(genAnalysis, mcTrack);
          // nMGen++;
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());

          if (okCell(FCPrefixEnum::Neg, kCh, kNeg, mcTrack.pt(), mcTrack.eta())) {
            nMGen++;
          }
        }
        // ----- Pions -----
        if (pdg == kPiPlus) {
          // fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          // nPiGen++;
          fillGenTrackQA<genAnalysisDir, kPi, kPos>(genAnalysis, mcTrack);
          if (okCell(FCPrefixEnum::PiPos, kPi, kPos, mcTrack.pt(), mcTrack.eta())) {
            nPiGen++;
          }
        } else if (pdg == kPiMinus) {
          // fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          // nAPiGen++;
          fillGenTrackQA<genAnalysisDir, kPi, kNeg>(genAnalysis, mcTrack);
          if (okCell(FCPrefixEnum::PiNeg, kPi, kNeg, mcTrack.pt(), mcTrack.eta())) {
            nAPiGen++;
          }
        } else if (pdg == kKPlus) {
          // fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          fillGenTrackQA<genAnalysisDir, kKa, kPos>(genAnalysis, mcTrack);
          if (mcTrack.pt() >= cfgEventSelection.cfgCutPtMinKa && mcTrack.pt() <= cfgEventSelection.cfgCutPtMaxKa && okCell(FCPrefixEnum::KaPos, kKa, kPos, mcTrack.pt(), mcTrack.eta())) {
            nKaGen++;
          }

        } else if (pdg == kKMinus) {
          // fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          fillGenTrackQA<genAnalysisDir, kKa, kNeg>(genAnalysis, mcTrack);
          if (mcTrack.pt() >= cfgEventSelection.cfgCutPtMinKa && mcTrack.pt() <= cfgEventSelection.cfgCutPtMaxKa && okCell(FCPrefixEnum::KaNeg, kKa, kNeg, mcTrack.pt(), mcTrack.eta())) {
            nAKaGen++;
          }
        } else if (pdg == kProton) {
          // fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          fillGenTrackQA<genAnalysisDir, kPr, kPos>(genAnalysis, mcTrack);
          if (mcTrack.pt() >= cfgEventSelection.cfgCutPtMinPr && mcTrack.pt() <= cfgEventSelection.cfgCutPtMaxPr && okCell(FCPrefixEnum::Pr, kPr, kPos, mcTrack.pt(), mcTrack.eta())) {
            nPrGen++;
          }
        } else if (pdg == kProtonBar) {
          // fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          fillGenTrackQA<genAnalysisDir, kPr, kNeg>(genAnalysis, mcTrack);
          if (mcTrack.pt() >= cfgEventSelection.cfgCutPtMinPr && mcTrack.pt() <= cfgEventSelection.cfgCutPtMaxPr && okCell(FCPrefixEnum::APr, kPr, kNeg, mcTrack.pt(), mcTrack.eta())) {
            nAPrGen++;
          }
        }

      } // particle/track loop for gen ends
      nChGen = nPGen - nMGen;
      nTGen = nPGen + nMGen;

      // ── Fill GEN sparse (denominator) ────────────────────────
      if (cfgEventSelection.fillSparseForGen) {
        hist.fill(HIST("sim/gen/sparse1"), nChGen, nPGen, nMGen,
                  nPrGen, nAPrGen, nKaGen, nAKaGen, nTGen,
                  cent);
        hist.fill(HIST("sim/gen/sparse2"), nChGen, nPGen, nMGen,
                  nPiGen, nAPiGen, nKaGen, nAKaGen, nTGen,
                  cent);
      }
      fillGenFactorialMoments<FCPrefixEnum::Pr>(nPrGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::APr>(nAPrGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::PiPos>(nPiGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::PiNeg>(nAPiGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::KaPos>(nKaGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::KaNeg>(nAKaGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::Pos>(nPGen, cent, registry);
      fillGenFactorialMoments<FCPrefixEnum::Neg>(nMGen, cent, registry);

      fillGenNetFactorialMoments<FCPrefixEnum::NetCh>(nPGen, nMGen, cent, registry);
      fillGenNetFactorialMoments<FCPrefixEnum::NetPr>(nPrGen, nAPrGen, cent, registry);
      fillGenNetFactorialMoments<FCPrefixEnum::NetPi>(nPiGen, nAPiGen, cent, registry);
      fillGenNetFactorialMoments<FCPrefixEnum::NetKa>(nKaGen, nAKaGen, cent, registry);
      // ── End of GEN level filling
      //
      // Numerator - Reconstructed + truth matched
      // reco->selFunc passed, no pdg
      // purity -> selFunc matched with pdg matched
      //
      float nPRec = 0, nMRec = 0, nChRec = 0, nTRec = 0;
      float nPrRec = 0, nAPrRec = 0;
      float nKaRec = 0, nAKaRec = 0;
      float nPiRec = 0, nAPiRec = 0;

      // purity counters — separate from reco
      float nPPur = 0, nMPur = 0, nChPur = 0, nTPur = 0;
      float nPrPur = 0, nAPrPur = 0;
      float nKaPur = 0, nAKaPur = 0;
      float nPiPur = 0, nAPiPur = 0;

      for (const auto& track : tracksTablePerColl) {

        if (!track.has_mcParticle()) {
          LOG(warning) << "No MC Particle for this track, skip...";
          continue;
        }
        auto mcPart = track.mcParticle();
        if (!mcPart.isPhysicalPrimary()) {
          continue;
        }

        if (cfgEventSelection.cfgRequireSameMcColl && mcPart.mcCollisionId() != mcCollision.globalIndex()) {
          continue;
        }

        int pdg = mcPart.pdgCode();

        int mcTag = getMCTag(track);

        fillTrackQA<qaTracksPostSel>(track);

        trackIsPion = false;
        trackIsKaon = false;
        trackIsProton = false;
        trackIdTag = 0;
        idMethodPi = kUnidentified;
        idMethodKa = kUnidentified;
        idMethodPr = kUnidentified;

        // Reject electrons first
        if (cfgEventSelection.cfgDoRejectionForId && !selTrackForId(track)) {
          continue;
        }

        // Fill separate spares for each species if it passes the cut
        if (selPion(track, idMethodPi)) {
          trackIsPion = true;
          BITSET(trackIdTag, IdBitPI);
          // Fill TPC sparse for pion
          hist.fill(HIST("PIDValidation/tpcSparse_Pi"), track.p(), track.tpcNSigmaPi(), trackIdTag, mcTag);
          // Fill TOF sparse for pion if has TOF
          if (track.hasTOF()) {
            hist.fill(HIST("PIDValidation/tofSparse_Pi"), track.p(), track.tofNSigmaPi(), trackIdTag, mcTag);
          }
        }
        if (selKaon(track, idMethodKa)) {
          trackIsKaon = true;
          BITSET(trackIdTag, IdBitKA);
          // Fill TPC sparse for kaon
          hist.fill(HIST("PIDValidation/tpcSparse_Ka"), track.p(), track.tpcNSigmaKa(), trackIdTag, mcTag);
          // Fill TOF sparse for kaon if has TOF
          if (track.hasTOF()) {
            hist.fill(HIST("PIDValidation/tofSparse_Ka"), track.p(), track.tofNSigmaKa(), trackIdTag, mcTag);
          }
        }
        if (selProton(track, idMethodPr)) {
          trackIsProton = true;
          BITSET(trackIdTag, IdBitPR);
          // Fill TPC sparse for proton
          hist.fill(HIST("PIDValidation/tpcSparse_Pr"), track.p(), track.tpcNSigmaPr(), trackIdTag, mcTag);
          // Fill TOF sparse for proton if has TOF
          if (track.hasTOF()) {
            hist.fill(HIST("PIDValidation/tofSparse_Pr"), track.p(), track.tofNSigmaPr(), trackIdTag, mcTag);
          }
        }

        // bool isKnownCharged = (std::abs(pdg) == kPiPlus  ||
        //                          std::abs(pdg) == kKPlus   ||
        //                          std::abs(pdg) == kProton  ||
        //                          std::abs(pdg) == kElectron ||
        //                          std::abs(pdg) == kMuonMinus ||
        //                          std::abs(pdg) == kDeuteron);

        ptEtaBin = hPtEtaForBinSearch->FindBin(track.pt(), track.eta()); // Find Track Bin for efficiency correction

        if (track.sign() > 0) {
          nPRec += hPtEtaForEffCorrection[kCh][kPos]->GetBinContent(ptEtaBin);
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h12_p"), track.p());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h13_pt"), track.pt());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h14_eta"), track.eta());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h15_phi"), track.phi());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h16_rapidity"), track.y());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h20_pt_eta"), track.pt(), track.eta());

          float weight = hPtEtaForEffCorrection[kCh][kPos]->GetBinContent(ptEtaBin);
          fillEffPower(posPow, weight);
        } else if (track.sign() < 0) {
          nMRec += hPtEtaForEffCorrection[kCh][kNeg]->GetBinContent(ptEtaBin);
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h12_p"), track.p());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h13_pt"), track.pt());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h14_eta"), track.eta());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h15_phi"), track.phi());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h16_rapidity"), track.y());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h20_pt_eta"), track.pt(), track.eta());

          float weight = hPtEtaForEffCorrection[kCh][kNeg]->GetBinContent(ptEtaBin);
          fillEffPower(negPow, weight);
        }

        // species reco — sel passes, PDG not checked (raw reco)
        if (trackIsPion) {
          if (track.sign() > 0) {
            nPiRec += hPtEtaForEffCorrection[kPi][kPos]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kPi, kPos>(recoAnalysis, track);

            float weight = hPtEtaForEffCorrection[kPi][kPos]->GetBinContent(ptEtaBin);
            fillEffPower(pipPow, weight);
          } else if (track.sign() < 0) {
            nAPiRec += hPtEtaForEffCorrection[kPi][kNeg]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kPi, kNeg>(recoAnalysis, track);

            float weight = hPtEtaForEffCorrection[kPi][kNeg]->GetBinContent(ptEtaBin);
            fillEffPower(pimPow, weight);
          }
          // PID band QA for pions
          if (idMethodPi == kTPCidentified) {
            fillIdentificationQA<qaTracksIdfd, kPi, tpcId>(hist, track);
          }
          if (idMethodPi == kTPCTOFidentified) {
            fillIdentificationQA<qaTracksIdfd, kPi, tpctofId>(hist, track);
          }
        } else if (trackIsKaon) {
          if (track.sign() > 0) {
            fillRecoTrackQA<recoAnalysisDir, kKa, kPos>(recoAnalysis, track);
            if (track.pt() >= cfgEventSelection.cfgCutPtMinKa && track.pt() <= cfgEventSelection.cfgCutPtMaxKa) {
              nKaRec += hPtEtaForEffCorrection[kKa][kPos]->GetBinContent(ptEtaBin);
              float weight = hPtEtaForEffCorrection[kKa][kPos]->GetBinContent(ptEtaBin);
              fillEffPower(kapPow, weight);
            }
          } else if (track.sign() < 0) {
            fillRecoTrackQA<recoAnalysisDir, kKa, kNeg>(recoAnalysis, track);
            if (track.pt() >= cfgEventSelection.cfgCutPtMinKa && track.pt() <= cfgEventSelection.cfgCutPtMaxKa) {
              nAKaRec += hPtEtaForEffCorrection[kKa][kNeg]->GetBinContent(ptEtaBin);
              float weight = hPtEtaForEffCorrection[kKa][kNeg]->GetBinContent(ptEtaBin);
              fillEffPower(kamPow, weight);
            }
          }
          // PID band QA for kaons
          if (idMethodKa == kTPCidentified) {
            fillIdentificationQA<qaTracksIdfd, kKa, tpcId>(hist, track);
          }
          if (idMethodKa == kTPCTOFidentified) {
            fillIdentificationQA<qaTracksIdfd, kKa, tpctofId>(hist, track);
          }
        } else if (trackIsProton) {
          if (track.sign() > 0) {
            fillRecoTrackQA<recoAnalysisDir, kPr, kPos>(recoAnalysis, track);
            if (track.pt() >= cfgEventSelection.cfgCutPtMinPr && track.pt() <= cfgEventSelection.cfgCutPtMaxPr) {
              nPrRec += hPtEtaForEffCorrection[kPr][kPos]->GetBinContent(ptEtaBin);
              float weight = hPtEtaForEffCorrection[kPr][kPos]->GetBinContent(ptEtaBin);
              fillEffPower(prPow, weight);
            }
          } else if (track.sign() < 0) {
            fillRecoTrackQA<recoAnalysisDir, kPr, kNeg>(recoAnalysis, track);
            if (track.pt() >= cfgEventSelection.cfgCutPtMinPr && track.pt() <= cfgEventSelection.cfgCutPtMaxPr) {
              nAPrRec += hPtEtaForEffCorrection[kPr][kNeg]->GetBinContent(ptEtaBin);
              float weight = hPtEtaForEffCorrection[kPr][kNeg]->GetBinContent(ptEtaBin);
              fillEffPower(aprPow, weight);
            }
          }
          // PID band QA for protons
          if (idMethodPr == kTPCidentified) {
            fillIdentificationQA<qaTracksIdfd, kPr, tpcId>(hist, track);
          }
          if (idMethodPr == kTPCTOFidentified) {
            fillIdentificationQA<qaTracksIdfd, kPr, tpctofId>(hist, track);
          }
        }
        // purity check - check pdg aginst sign
        bool purityPion = false;
        bool purityKaon = false;
        bool purityProton = false;

        if (trackIsPion) {
          if (track.sign() > 0 && pdg == kPiPlus) {
            purityPion = true;
          }
          if (track.sign() < 0 && pdg == kPiMinus) {
            purityPion = true;
          }
        } else if (trackIsKaon) {
          if (track.sign() > 0 && pdg == kKPlus) {
            purityKaon = true;
          }
          if (track.sign() < 0 && pdg == kKMinus) {
            purityKaon = true;
          }
        } else if (trackIsProton) {
          if (track.sign() > 0 && pdg == kProton) {
            purityProton = true;
          }
          if (track.sign() < 0 && pdg == kProtonBar) {
            purityProton = true;
          }
        }

        // charge purity — track.sign() + isKnownCharged + PDG sign consistency
        bool pdgPositive = (pdg == kPiPlus || pdg == kKPlus || pdg == kProton ||
                            pdg == kPositron || pdg == kMuonPlus || pdg == kDeuteron);
        bool pdgNegative = (pdg == kPiMinus || pdg == kKMinus || pdg == kProtonBar ||
                            pdg == kElectron || pdg == kMuonMinus || pdg == -kDeuteron);

        if (track.sign() > 0 && pdgPositive) {
          nPPur++;
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h12_p"), track.p());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h13_pt"), track.pt());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h14_eta"), track.eta());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h15_phi"), track.phi());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h16_rapidity"), track.y());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h20_pt_eta"), track.pt(), track.eta());
        } else if (track.sign() < 0 && pdgNegative) {
          nMPur++;
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h12_p"), track.p());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h13_pt"), track.pt());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h14_eta"), track.eta());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h15_phi"), track.phi());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h16_rapidity"), track.y());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h20_pt_eta"), track.pt(), track.eta());
        }

        // species purity fill
        if (purityPion) {
          if (track.sign() > 0) {
            nPiPur++;
            fillPurityTrackQA<purityAnalysisDir, kPi, kPos>(purityAnalysis, track);
          } else {
            nAPiPur++;
            fillPurityTrackQA<purityAnalysisDir, kPi, kNeg>(purityAnalysis, track);
          }
        } else if (purityKaon) {
          if (track.sign() > 0) {
            nKaPur++;
            fillPurityTrackQA<purityAnalysisDir, kKa, kPos>(purityAnalysis, track);
          } else {
            nAKaPur++;
            fillPurityTrackQA<purityAnalysisDir, kKa, kNeg>(purityAnalysis, track);
          }
        } else if (purityProton) {
          if (track.sign() > 0) {
            nPrPur++;
            fillPurityTrackQA<purityAnalysisDir, kPr, kPos>(purityAnalysis, track);
          } else {
            nAPrPur++;
            fillPurityTrackQA<purityAnalysisDir, kPr, kNeg>(purityAnalysis, track);
          }
        }
      } // reconstructed track loop ends
      nChRec = nPRec - nMRec;
      nTRec = nPRec + nMRec;
      nChPur = nPPur - nMPur;
      nTPur = nPPur + nMPur;

      // ── fill reco histos ─────────────────────────────────────
      if (cfgEventSelection.fillSparseForReco) {
        hist.fill(HIST("sim/reco/sparse1"), nChRec, nPRec, nMRec,
                  nPrRec, nAPrRec, nKaRec, nAKaRec, nTRec, cent);
        hist.fill(HIST("sim/reco/sparse2"), nChRec, nPRec, nMRec,
                  nPiRec, nAPiRec, nKaRec, nAKaRec, nTRec, cent);
      }

      // ── fill purity histos ───────────────────────────────────
      if (cfgEventSelection.fillSparseForPurity) {
        hist.fill(HIST("sim/purity/sparse1"), nChPur, nPPur, nMPur,
                  nPrPur, nAPrPur, nKaPur, nAKaPur, nTPur, cent);
        hist.fill(HIST("sim/purity/sparse2"), nChPur, nPPur, nMPur,
                  nPiPur, nAPiPur, nKaPur, nAKaPur, nTPur, cent);
      }
      fillCollQA<qaEventPostSel>(col, nChRec, nTRec);

      fillCountCheck<FCPrefixEnum::Pr>(genAll[0], nPrGen, prPow, cent, registry);
      fillCountCheck<FCPrefixEnum::APr>(genAll[1], nAPrGen, aprPow, cent, registry);
      fillCountCheck<FCPrefixEnum::PiPos>(genAll[2], nPiGen, pipPow, cent, registry);
      fillCountCheck<FCPrefixEnum::PiNeg>(genAll[3], nAPiGen, pimPow, cent, registry);
      fillCountCheck<FCPrefixEnum::KaPos>(genAll[4], nKaGen, kapPow, cent, registry);
      fillCountCheck<FCPrefixEnum::KaNeg>(genAll[5], nAKaGen, kamPow, cent, registry);
      fillCountCheck<FCPrefixEnum::Pos>(genAll[6], nPGen, posPow, cent, registry);
      fillCountCheck<FCPrefixEnum::Neg>(genAll[7], nMGen, negPow, cent, registry);

      fillFCBasis<FCPrefixEnum::Pr>(prPow, cent, registry);
      fillFCBasis<FCPrefixEnum::APr>(aprPow, cent, registry);
      fillFCBasis<FCPrefixEnum::PiPos>(pipPow, cent, registry);
      fillFCBasis<FCPrefixEnum::PiNeg>(pimPow, cent, registry);
      fillFCBasis<FCPrefixEnum::KaPos>(kapPow, cent, registry);
      fillFCBasis<FCPrefixEnum::KaNeg>(kamPow, cent, registry);
      fillFCBasis<FCPrefixEnum::Pos>(posPow, cent, registry);
      fillFCBasis<FCPrefixEnum::Neg>(negPow, cent, registry);

      fillNetQVectorProfileHistograms<FCPrefixEnum::NetCh>(posPow, negPow, cent, registry);
      fillNetQVectorProfileHistograms<FCPrefixEnum::NetPr>(prPow, aprPow, cent, registry);
      fillNetQVectorProfileHistograms<FCPrefixEnum::NetPi>(pipPow, pimPow, cent, registry);
      fillNetQVectorProfileHistograms<FCPrefixEnum::NetKa>(kapPow, kamPow, cent, registry);
    } // common collision loop ends
  } // process sim ends
  PROCESS_SWITCH(NchCumulantsId, processSim, "Process Sim: Gen + Reco + Purity", true);
}; // structure ends

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NchCumulantsId>(cfgc)};
}
