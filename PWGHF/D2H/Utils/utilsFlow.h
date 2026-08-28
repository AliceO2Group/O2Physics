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

/// \file utilsFlow.h
/// \brief Utilities for flow analyses
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#ifndef PWGHF_D2H_UTILS_UTILSFLOW_H_
#define PWGHF_D2H_UTILS_UTILSFLOW_H_

#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>


namespace o2::analysis
{
namespace hf_flow_utils
{
enum QvecEstimator { FV0A = 0,
                     FT0M,
                     FT0A,
                     FT0C,
                     TPCPos,
                     TPCNeg,
                     TPCTot };

class QVectorResoHelper {
 public:

  QVectorResoHelper() = default;

  /// Call once per run / timestamp update
  void update(o2::framework::Service<o2::ccdb::BasicCCDBManager>& ccdb,
              std::string const& basePath,
              std::string const& tripletPath,
              int64_t timestamp,
              int runNumber)
  {
    if (mCurrentRun == runNumber && mResoValues != nullptr) {
      return;
    }

    std::string fullPath = basePath + "/" + tripletPath + "/";
    LOG(info) << "Fetching resolution array from CCDB path: " << fullPath << " for timestamp " << timestamp << " and run number " << runNumber;

    mResoValues = ccdb->getForTimeStamp<std::array<float, 100>>(fullPath, timestamp);
    if (!mResoValues) {
      LOGF(fatal, "Resolution array not found in CCDB at path: %s", fullPath.c_str());
    }

    mCurrentRun = runNumber;
  }

  /// Get pointer to the full array
  const std::array<float, 100>* getResoValues() const { return mResoValues; }

  /// Direct array access with simple float return type
  float operator[](size_t index) const
  {
    assert(mResoValues != nullptr && "CCDB resolution array accessed before initialization!");
    return (*mResoValues)[index];
  }

 private:
  const std::array<float, 100>* mResoValues{nullptr};
  int64_t mCurrentRun{-1};
};

/// Compute the delta psi in the range [0, pi/harmonic]
/// \param psi1 is the first angle
/// \param psi2 is the second angle
/// \param harmonic is the harmonic
/// \note Ported from AliAnalysisTaskSECharmHadronvn::GetDeltaPsiSubInRange
float getDeltaPsiInRange(float psi1, float psi2, int harmonic)
{
  float deltaPsi = psi1 - psi2;
  deltaPsi = RecoDecay::constrainAngle(deltaPsi, -o2::constants::math::PI / harmonic, harmonic);
  return deltaPsi;
}

/// Get the Q vector
template <typename T>
concept HasQvecFT0A = requires(T collision) {
  collision.qvecFT0ARe();
  collision.qvecFT0AIm();
};

template <typename T>
concept HasQvecFT0C = requires(T collision) {
  collision.qvecFT0CRe();
  collision.qvecFT0CIm();
};

template <typename T>
concept HasQvecFT0M = requires(T collision) {
  collision.qvecFT0MRe();
  collision.qvecFT0MIm();
};

template <typename T>
concept HasQvecFV0A = requires(T collision) {
  collision.qvecFV0ARe();
  collision.qvecFV0AIm();
};

template <typename T>
concept HasQvecTPCpos = requires(T collision) {
  collision.qvecBPosRe();
  collision.qvecBPosIm();
};

template <typename T>
concept HasQvecTPCneg = requires(T collision) {
  collision.qvecBNegRe();
  collision.qvecBNegIm();
};

template <typename T>
concept HasQvecTPCtot = requires(T collision) {
  collision.qvecBTotRe();
  collision.qvecBTotIm();
};

/// Get the Q vector using FT0A estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecFT0A TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A()};
}

/// Get the Q vector using FT0C estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecFT0C TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C()};
}

/// Get the Q vector using FT0C estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecFT0M TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M()};
}

/// Get the Q vector using FV0A estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecFV0A TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecFV0ARe(), collision.qvecFV0AIm(), collision.sumAmplFV0A()};
}

/// Get the Q vector using TPCpos estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecTPCpos TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecBPosRe(), collision.qvecBPosIm(), collision.nTrkBPos()};
}

/// Get the Q vector using TPCneg estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecTPCneg TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecBNegRe(), collision.qvecBNegIm(), collision.nTrkBNeg()};
}

/// Get the Q vector using TPCtot estimator
/// \param collision is the collision
/// \return Q vector of the collision
template <HasQvecTPCtot TCollision>
std::array<float, 3> getQvec(const TCollision& collision)
{
  return std::array<float, 3>{collision.qvecBTotRe(), collision.qvecBTotIm(), collision.nTrkBTot()};
}

/// Get the Q vector choosing your favourite estimator
/// \param collision is the collision with the Q vector information
/// \param qvecEst is the chosen Q-vector estimator
template <typename TCollision>
std::array<float, 3> getQvec(TCollision const& collision, const int qvecEst)
{
  switch (qvecEst) {
    case QvecEstimator::FV0A:
      if constexpr (HasQvecFV0A<TCollision>) {
        return std::array<float, 3>{collision.qvecFV0ARe(), collision.qvecFV0AIm(), collision.sumAmplFV0A()};
      }
      break;
    case QvecEstimator::FT0A:
      if constexpr (HasQvecFT0A<TCollision>) {
        return std::array<float, 3>{collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.sumAmplFT0A()};
      }
      break;
    case QvecEstimator::FT0C:
      if constexpr (HasQvecFT0C<TCollision>) {
        return std::array<float, 3>{collision.qvecFT0CRe(), collision.qvecFT0CIm(), collision.sumAmplFT0C()};
      }
      break;
    case QvecEstimator::FT0M:
      if constexpr (HasQvecFT0M<TCollision>) {
        return std::array<float, 3>{collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.sumAmplFT0M()};
      }
      break;
    case QvecEstimator::TPCPos:
      if constexpr (HasQvecTPCpos<TCollision>) {
        return std::array<float, 3>{collision.qvecBPosRe(), collision.qvecBPosIm(), static_cast<float>(collision.nTrkBPos())};
      }
      break;
    case QvecEstimator::TPCNeg:
      if constexpr (HasQvecTPCneg<TCollision>) {
        return std::array<float, 3>{collision.qvecBNegRe(), collision.qvecBNegIm(), static_cast<float>(collision.nTrkBNeg())};
      }
      break;
    case QvecEstimator::TPCTot:
      if constexpr (HasQvecTPCtot<TCollision>) {
        return std::array<float, 3>{collision.qvecBTotRe(), collision.qvecBTotIm(), static_cast<float>(collision.nTrkBTot())};
      }
      break;
    default:
      LOGP(fatal, "Q-vector estimator not valid. Please choose between FV0A, FT0M, FT0A, FT0C, TPCPos, TPCNeg, TPCTot");
      break;
  }
  return std::array<float, 3>{-999.f, -999.f, -999.f};
}

/// Get the EsE Q vector
template <typename T>
concept HasEsEQvecFT0A = requires(T collision) {
  collision.eseQvecFT0ARe();
  collision.eseQvecFT0AIm();
};

template <typename T>
concept HasEsEQvecFT0C = requires(T collision) {
  collision.eseQvecFT0CRe();
  collision.eseQvecFT0CIm();
};

template <typename T>
concept HasEsEQvecFT0M = requires(T collision) {
  collision.eseQvecFT0MRe();
  collision.eseQvecFT0MIm();
};

template <typename T>
concept HasEsEQvecFV0A = requires(T collision) {
  collision.eseQvecFV0ARe();
  collision.eseQvecFV0AIm();
};

template <typename T>
concept HasEsEQvecTPCpos = requires(T collision) {
  collision.eseQvecTPCposRe();
  collision.eseQvecTPCposIm();
};

template <typename T>
concept HasEsEQvecTPCneg = requires(T collision) {
  collision.eseQvecTPCnegRe();
  collision.eseQvecTPCnegIm();
};

template <typename T>
concept HasEsEQvecTPCtot = requires(T collision) {
  collision.eseQvecTPCallRe();
  collision.eseQvecTPCallIm();
};

/// Get the Ese Q vector choosing your favourite estimator
/// \param collision is the collision with the Q vector information
/// \param qvecEst is the chosen Q-vector estimator
template <typename TCollision>
std::array<float, 3> getEseQvec(TCollision const& collision, const int qvecEst)
{
  switch (qvecEst) {
    case QvecEstimator::FV0A:
      if constexpr (HasEsEQvecFV0A<TCollision>) {
        return std::array<float, 3>{collision.eseQvecFV0ARe(), collision.eseQvecFV0AIm(), collision.sumAmplFV0A()};
      }
      break;
    case QvecEstimator::FT0A:
      if constexpr (HasEsEQvecFT0A<TCollision>) {
        return std::array<float, 3>{collision.eseQvecFT0ARe(), collision.eseQvecFT0AIm(), collision.sumAmplFT0A()};
      }
      break;
    case QvecEstimator::FT0C:
      if constexpr (HasEsEQvecFT0C<TCollision>) {
        return std::array<float, 3>{collision.eseQvecFT0CRe(), collision.eseQvecFT0CIm(), collision.sumAmplFT0C()};
      }
      break;
    case QvecEstimator::FT0M:
      if constexpr (HasEsEQvecFT0M<TCollision>) {
        return std::array<float, 3>{collision.eseQvecFT0MRe(), collision.eseQvecFT0MIm(), collision.sumAmplFT0M()};
      }
      break;
    case QvecEstimator::TPCPos:
      if constexpr (HasEsEQvecTPCpos<TCollision>) {
        return std::array<float, 3>{collision.eseQvecTPCposRe(), collision.eseQvecTPCposIm(), static_cast<float>(collision.nTrkTPCpos())};
      }
      break;
    case QvecEstimator::TPCNeg:
      if constexpr (HasEsEQvecTPCneg<TCollision>) {
        return std::array<float, 3>{collision.eseQvecTPCnegRe(), collision.eseQvecTPCnegIm(), static_cast<float>(collision.nTrkTPCneg())};
      }
      break;
    case QvecEstimator::TPCTot:
      if constexpr (HasEsEQvecTPCtot<TCollision>) {
        return std::array<float, 3>{collision.eseQvecTPCallRe(), collision.eseQvecTPCallIm(), static_cast<float>(collision.nTrkTPCall())};
      }
      break;
    default:
      LOGP(fatal, "Q-vector estimator not valid. Please choose between FV0A, FT0M, FT0A, FT0C, TPCPos, TPCNeg, TPCTot");
      break;
  }
  return std::array<float, 3>{-999.f, -999.f, -999.f};
}

} // namespace hf_flow_utils
} // namespace o2::analysis

#endif // PWGHF_D2H_UTILS_UTILSFLOW_H_
