// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file sigma0Builder.h
/// \brief Sigma0 <- Lambda + gamma builder for femto
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#ifndef PWGCF_FEMTO_CORE_SIGMA0BUILDER_H_
#define PWGCF_FEMTO_CORE_SIGMA0BUILDER_H_

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/photonBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <string>

namespace o2::analysis::femto::sigma0builder
{

// Selection of the Lambda and photon daughters. Applied with partitions in the producer, before the builder.
constexpr const char PrefixSigma0LambdaSelection[] = "Sigma0LambdaSelection";
using ConfSigma0LambdaSelection = v0builder::ConfLambdaSelection<PrefixSigma0LambdaSelection, 0>;
constexpr const char PrefixSigma0PhotonSelection[] = "Sigma0PhotonSelection";
using ConfSigma0PhotonSelection = photonbuilder::ConfPhotonSelection<PrefixSigma0PhotonSelection>;

// Loose kinematic filters on the reconstructed Lambda-photon pair
struct ConfSigma0Filters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Filters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 1.1f, "Minimum invariant mass"};
  o2::framework::Configurable<float> massMax{"massMax", 1.3f, "Maximum invariant mass"};
};

// Base selection configuration for downstream analysis tasks.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                             \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "PDG code. Set sign to -1 for antiparticle"};   \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                       \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                     \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                                  \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                                   \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                    \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                       \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass"};                          \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass"};                          \
  o2::framework::Configurable<datatypes::V0MaskType> lambdaMask{"lambdaMask", 0, "Bitmask required for the Lambda daughter"}; \
  o2::framework::Configurable<datatypes::PhotonMaskType> photonMask{"photonMask", 0, "Bitmask required for the photon daughter"};

template <auto& Prefix, int DefaultSign = 1>
struct ConfSigma0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  SIGMA0_DEFAULT_SELECTIONS(1.1, 1.3, 3212)
  o2::framework::Configurable<int> sign{"sign", DefaultSign, "Sign of the Sigma0 (+1: Sigma0; -1: AntiSigma0; 0: both)"};
};

#undef SIGMA0_DEFAULT_SELECTIONS

constexpr const char PrefixSigma0Selection[] = "Sigma0Selection";
using ConfSigma0SelectionDefault = ConfSigma0Selection<PrefixSigma0Selection>;
constexpr const char PrefixSigma0SelectionStrangeness[] = "Sigma0SelectionStrangeness";
using ConfSigma0SelectionStrangeness = ConfSigma0Selection<PrefixSigma0SelectionStrangeness, 0>;

struct Sigma0BuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FSigma0s> producedSigma0s;
  o2::framework::Produces<o2::aod::FLiteSigma0s> producedLiteSigma0s;
  o2::framework::Produces<o2::aod::FSigma0Masks> producedSigma0Masks;
  o2::framework::Produces<o2::aod::FSigma0Extras> producedSigma0Extras;
};

struct ConfSigma0Tables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Tables");
  o2::framework::Configurable<int> produceSigma0s{"produceSigma0s", -1, "Produce Sigma0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteSigma0s{"produceLiteSigma0s", -1, "Produce LiteSigma0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigma0Masks{"produceSigma0Masks", -1, "Produce Sigma0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigma0Extras{"produceSigma0Extras", -1, "Produce Sigma0Extras (-1: auto; 0 off; 1 on). Needs FLambdaExtras and FPhotonExtras as input"};
};

/// \brief Builder for Sigma0 <- Lambda + PCM gamma candidates
/// Pairs daughters within the collision and applies loose kinematic cuts
class Sigma0Builder
{
 public:
  Sigma0Builder() = default;
  ~Sigma0Builder() = default;

  template <typename T1, typename T2, typename T3>
  void init(T1& confFilter, T2& confTable, T3& initContext)
  {
    LOG(info) << "Initialize femto Sigma0 builder...";
    mMassMin = confFilter.massMin.value;
    mMassMax = confFilter.massMax.value;
    mPtMin = confFilter.ptMin.value;
    mPtMax = confFilter.ptMax.value;
    mEtaMin = confFilter.etaMin.value;
    mEtaMax = confFilter.etaMax.value;
    mPhiMin = confFilter.phiMin.value;
    mPhiMax = confFilter.phiMax.value;

    mProduceSigma0s = utils::enableTable("FSigma0s_001", confTable.produceSigma0s.value, initContext);
    mProduceLiteSigma0s = utils::enableTable("FLiteSigma0s_001", confTable.produceLiteSigma0s.value, initContext);
    mProduceSigma0Masks = utils::enableTable("FSigma0Masks_001", confTable.produceSigma0Masks.value, initContext);
    mProduceSigma0Extras = utils::enableTable("FSigma0Extras_001", confTable.produceSigma0Extras.value, initContext);

    if (mProduceSigma0s && mProduceLiteSigma0s) {
      LOG(fatal) << "FSigma0s and FLiteSigma0s are mutually exclusive -- enable only one. "
                 << "FLiteSigma0s is meant to replace FSigma0s at the producer stage for better compression in derived data.";
    }

    mFillAnyTable = mProduceSigma0s || mProduceLiteSigma0s || mProduceSigma0Masks || mProduceSigma0Extras;
    if (!mFillAnyTable) {
      LOG(info) << "No tables configured...";
    }
    LOG(info) << "Initialization done...";
  }

  /// Pair all Lambdas with all photons of one collision
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillSigma0s(T1 const& col, T2& products, T3& lambdaPartition, T4& photonPartition, T5& cache)
  {
    if (!mFillAnyTable) {
      return;
    }
    auto lambdaSlice = lambdaPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto photonSlice = photonPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& [lambda, photon] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(lambdaSlice, photonSlice))) {
      this->fillSigma0(col, lambda, photon, products);
    }
  }

 private:
  template <typename T1, typename T2>
  void reconstructSigma0(T1 const& lambda, T2 const& photon)
  {
    // use the PDG Lambda mass to improve the Sigma0 mass resolution
    mVecLambda = ROOT::Math::PtEtaPhiMVector{lambda.pt(), lambda.eta(), lambda.phi(), o2::constants::physics::MassLambda0};
    mVecPhoton = ROOT::Math::PtEtaPhiMVector{photon.pt(), photon.eta(), photon.phi(), o2::constants::physics::MassPhoton};
    ROOT::Math::PtEtaPhiMVector vecSigma0 = mVecLambda + mVecPhoton;

    mPt = vecSigma0.Pt();
    mEta = vecSigma0.Eta();
    mPhi = RecoDecay::constrainAngle(vecSigma0.Phi());
    mMass = vecSigma0.M();
  }

  [[nodiscard]] bool checkFilters() const
  {
    return ((mMass > mMassMin && mMass < mMassMax) &&
            (mPt > mPtMin && mPt < mPtMax) &&
            (mEta > mEtaMin && mEta < mEtaMax) &&
            (mPhi >= mPhiMin && mPhi < mPhiMax));
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillSigma0(T1 const& col, T2 const& lambda, T3 const& photon, T4& products)
  {
    reconstructSigma0(lambda, photon);
    if (!checkFilters()) {
      return;
    }

    // the Sigma0 inherits the sign of the Lambda (Sigma0 -> Lambda gamma, AntiSigma0 -> AntiLambda gamma)
    const float signedPt = lambda.sign() * mPt;

    if (mProduceSigma0s) {
      products.producedSigma0s(col.globalIndex(),
                               signedPt,
                               mEta,
                               mPhi,
                               mMass,
                               lambda.posDauId(),
                               lambda.negDauId(),
                               photon.globalIndex());
    }
    if (mProduceLiteSigma0s) {
      products.producedLiteSigma0s(col.globalIndex(),
                                   o2::aod::femtobase::lite::binSignedPt(signedPt),
                                   o2::aod::femtobase::lite::binEta(mEta),
                                   o2::aod::femtobase::lite::binPhi(mPhi),
                                   o2::aod::femtosigma0s::lite::binSigma0Mass(mMass),
                                   lambda.posDauId(),
                                   lambda.negDauId(),
                                   photon.globalIndex());
    }
    if (mProduceSigma0Masks) {
      products.producedSigma0Masks(lambda.mask(), photon.mask());
    }
    if (mProduceSigma0Extras) {
      if constexpr (requires { lambda.cosPa(); lambda.dauDca(); lambda.transRadius(); photon.cosPa(); photon.v0Radius(); photon.chi2Ndf(); }) {
        products.producedSigma0Extras(lambda.mass(),
                                      lambda.pt(),
                                      lambda.cosPa(),
                                      lambda.dauDca(),
                                      lambda.transRadius(),
                                      photon.pt(),
                                      photon.cosPa(),
                                      photon.v0Radius(),
                                      photon.chi2Ndf(),
                                      ROOT::Math::VectorUtil::Angle(mVecLambda, mVecPhoton));
      } else {
        LOG(fatal) << "FSigma0Extras requested, but FLambdaExtras/FPhotonExtras are not joined to the input Lambda/photon tables";
      }
    }
  }

  // cached kinematics of the daughters and the Sigma0
  ROOT::Math::PtEtaPhiMVector mVecLambda;
  ROOT::Math::PtEtaPhiMVector mVecPhoton;
  float mPt = 0.f;
  float mEta = 0.f;
  float mPhi = 0.f;
  float mMass = 0.f;

  float mMassMin = 0.f;
  float mMassMax = 0.f;
  float mPtMin = 0.f;
  float mPtMax = 0.f;
  float mEtaMin = 0.f;
  float mEtaMax = 0.f;
  float mPhiMin = 0.f;
  float mPhiMax = 0.f;

  bool mFillAnyTable = false;
  bool mProduceSigma0s = false;
  bool mProduceLiteSigma0s = false;
  bool mProduceSigma0Masks = false;
  bool mProduceSigma0Extras = false;
};

} // namespace o2::analysis::femto::sigma0builder
#endif // PWGCF_FEMTO_CORE_SIGMA0BUILDER_H_
