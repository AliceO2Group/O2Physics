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

/// \file CentralityEstimation.h
/// \brief Definitions of centrality estimators
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_CORE_CENTRALITYESTIMATION_H_
#define PWGHF_CORE_CENTRALITYESTIMATION_H_

namespace o2::hf_centrality
{
// centrality selection estimators
enum CentralityEstimator {
  None = 0,
  FT0A,
  FT0C,
  FT0M,
  FV0A,
  NTracksPV,
  NCentralityEstimators
};

template <typename T>
concept hasFT0ACent = requires(T collision)
{
  collision.centFT0A();
};

template <typename T>
concept hasFT0CCent = requires(T collision)
{
  collision.centFT0C();
};

template <typename T>
concept hasFT0MCent = requires(T collision)
{
  collision.centFT0M();
};

template <typename T>
concept hasFV0ACent = requires(T collision)
{
  collision.centFV0A();
};

template <typename T>
concept hasNTracksPVCent = requires(T collision)
{
  collision.centNTPV();
};

/// Evaluate centrality/multiplicity percentile using FT0A estimator
/// \param candidate is candidate
/// \return centrality/multiplicity percentile of the collision
template <hasFT0ACent Coll>
float getCentralityColl(const Coll& collision)
{
  return collision.centFT0A();
}

/// Evaluate centrality/multiplicity percentile using FT0C estimator
/// \param candidate is candidate
/// \return centrality/multiplicity percentile of the collision
template <hasFT0CCent Coll>
float getCentralityColl(const Coll& collision)
{
  return collision.centFT0C();
}

/// Evaluate centrality/multiplicity percentile using FT0M estimator
/// \param candidate is candidate
/// \return centrality/multiplicity percentile of the collision
template <hasFT0MCent Coll>
float getCentralityColl(const Coll& collision)
{
  return collision.centFT0M();
}

/// Evaluate centrality/multiplicity percentile using FV0A estimator
/// \param candidate is candidate
/// \return centrality/multiplicity percentile of the collision
template <hasFV0ACent Coll>
float getCentralityColl(const Coll& collision)
{
  return collision.centFV0A();
}

/// Evaluate centrality/multiplicity percentile using NTracksPV estimator
/// \param candidate is candidate
/// \return centrality/multiplicity percentile of the collision
template <hasNTracksPVCent Coll>
float getCentralityColl(const Coll& collision)
{
  return collision.centNTPV();
}

/// Default case if no centrality/multiplicity estimator is provided
/// \param candidate is candidate
/// \return dummy value for centrality/multiplicity percentile of the collision
template <typename Coll>
float getCentralityColl(const Coll&)
{
  return 105.0f;
}

/// Get the centrality
/// \param collision is the collision with the centrality information
/// \return collision centrality
template <typename Coll>
float getCentralityColl(const Coll& collision, int centEstimator)
{
  float cent = -999.;
  switch (centEstimator) {
    case CentralityEstimator::FT0A:
      if constexpr (hasFT0ACent<Coll>) {
        cent = collision.centFT0A();
      } else {
        LOG(warning) << "Warning: Collision does not have centFT0A().\n";
      }
      break;
    case CentralityEstimator::FT0C:
      if constexpr (hasFT0CCent<Coll>) {
        cent = collision.centFT0C();
      } else {
        LOG(warning) << "Warning: Collision does not have centFT0C().\n";
      }
      break;
    case CentralityEstimator::FT0M:
      if constexpr (hasFT0MCent<Coll>) {
        cent = collision.centFT0M();
      } else {
        LOG(warning) << "Warning: Collision does not have centFT0M().\n";
      }
      break;
    case CentralityEstimator::FV0A:
      if constexpr (hasFV0ACent<Coll>) {
        cent = collision.centFV0A();
      } else {
        LOG(warning) << "Warning: Collision does not have centFV0A().\n";
      }
      break;
    default:
      LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to V0A";
      cent = collision.centFV0A();
      break;
  }
  return cent;
}

/// \brief Function to get MC collision centrality
/// \param collSlice collection of reconstructed collisions associated to a generated one
/// \return generated MC collision centrality
template <typename CCs>
float getCentralityGenColl(CCs const& collSlice)
{
  float centrality{-1};
  float multiplicity{0.f};
  for (const auto& collision : collSlice) {
    float collMult = collision.numContrib();
    if (collMult > multiplicity) {
      centrality = getCentralityColl(collision);
      multiplicity = collMult;
    }
  }
  return centrality;
}

/// \brief Function to get MC collision centrality
/// \param collSlice collection of reconstructed collisions associated to a generated one
/// \return generated MC collision centrality
template <typename CCs>
float getCentralityGenColl(CCs const& collSlice, int centEstimator)
{
  float centrality{-1};
  float multiplicity{0.f};
  for (const auto& collision : collSlice) {
    float collMult = collision.numContrib();
    if (collMult > multiplicity) {
      centrality = getCentralityColl(collision, centEstimator);
      multiplicity = collMult;
    }
  }
  return centrality;
}

/// Get the centrality
/// \param collision is the collision with the centrality information
/// \param centEstimator is the centrality estimator from hf_centrality::CentralityEstimator
template <o2::hf_centrality::CentralityEstimator centEstimator, typename Coll>
float evalCentralityColl(Coll const& collision)
{
  if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0A) {
    return collision.centFT0A();
  } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
    return collision.centFT0C();
  } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
    return collision.centFT0M();
  } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FV0A) {
    return collision.centFV0A();
  } else {
    LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to FT0c";
    return collision.centFT0C();
  }
}

/// \brief Function to get MC collision centrality
/// \param collSlice collection of reconstructed collisions associated to a generated one
/// \return generated MC collision centrality
template <o2::hf_centrality::CentralityEstimator centEstimator, typename CCs>
float evalCentralityGenColl(CCs const& collSlice)
{
  float centrality{-1};
  float multiplicity{0.f};
  for (const auto& collision : collSlice) {
    float collMult = collision.numContrib();
    if (collMult > multiplicity) {
      centrality = evalCentralityColl<centEstimator>(collision);
      multiplicity = collMult;
    }
  }
  return centrality;
}

} // namespace o2::hf_centrality

#endif // PWGHF_CORE_CENTRALITYESTIMATION_H_
