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

/// \file twoPhotonResonanceBuilder.h
/// \brief two photon resonance builder (pi0, eta, ... -> gamma gamma)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TWOPHOTONRESONANCEBUILDER_H_
#define PWGCF_FEMTO_CORE_TWOPHOTONRESONANCEBUILDER_H_

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <string>

namespace o2::analysis::femto::twophotonresonancebuilder
{
template <auto& Prefix>
struct ConfTwoPhotonResonanceFilters : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 10.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.f, "Minimum invariant mass for the resonance"};
  o2::framework::Configurable<float> massMax{"massMax", 1.f, "Maximum invariant mass for the resonance"};
};
constexpr const char PrefixPi0Filters[] = "Pi0Filters1";
constexpr const char PrefixEtaFilters[] = "EtaFilters1";

using ConfPi0Filters = ConfTwoPhotonResonanceFilters<PrefixPi0Filters>;
using ConfEtaFilters = ConfTwoPhotonResonanceFilters<PrefixEtaFilters>;

// selection used downstream (QA/pairing) to re-select already produced resonances by mass window + daughter mask.
// unlike TWOTRACKRESONANCE_DEFAULT_SELECTION there is no pos/neg split and no momentum-threshold PID
// switch: photon daughters are unordered and have no momentum-dependent PID regime.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOPHOTONRESONANCE_DEFAULT_SELECTION(defaultMassMin, defaultMassMax)                                               \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                    \
  o2::framework::Configurable<float> ptMax{"ptMax", 10.f, "Maximum pT"};                                                   \
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};                                               \
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};                                                \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                 \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                    \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass for the resonance"};     \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass for the resonance"};     \
  o2::framework::Configurable<datatypes::PhotonMaskType> dau1Mask{"dau1Mask", 0, "Bitmask required for first photon daughter"};  \
  o2::framework::Configurable<datatypes::PhotonMaskType> dau2Mask{"dau2Mask", 0, "Bitmask required for second photon daughter"};

struct ConfPi0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Pi0Selection1");
  TWOPHOTONRESONANCE_DEFAULT_SELECTION(0.1f, 0.17f)
};

struct ConfEtaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("EtaSelection1");
  TWOPHOTONRESONANCE_DEFAULT_SELECTION(0.4f, 0.7f)
};

#undef TWOPHOTONRESONANCE_DEFAULT_SELECTION

struct TwoPhotonResonanceBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FPi0s> producedPi0s;
  o2::framework::Produces<o2::aod::FPi0Masks> producedPi0Masks;
  o2::framework::Produces<o2::aod::FEtas> producedEtas;
  o2::framework::Produces<o2::aod::FEtaMasks> producedEtaMasks;
};

struct ConfTwoPhotonResonanceTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TwoPhotonResonanceTables");
  o2::framework::Configurable<int> producePi0s{"producePi0s", -1, "Produce Pi0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePi0Masks{"producePi0Masks", -1, "Produce Pi0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceEtas{"produceEtas", -1, "Produce Etas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceEtaMasks{"produceEtaMasks", -1, "Produce EtaMasks (-1: auto; 0 off; 1 on)"};
};

/// \brief Builder for resonances reconstructed from two PCM photons (pi0, eta, ... -> gamma gamma)
///
/// Unlike TwoTrackResonanceBuilder, there is no pos/neg daughter distinction (photons are their own
/// antiparticle) and no momentum-threshold PID switch (photon selection has no momentum-dependent PID
/// regime the way track PID does), so daughter combinatorics run over a single photon partition using
/// CombinationsStrictlyUpperIndexPolicy to avoid self-pairing/double-counting.
template <modes::TwoPhotonResonance resoType>
class TwoPhotonResonanceBuilder
{
 public:
  TwoPhotonResonanceBuilder() = default;
  ~TwoPhotonResonanceBuilder() = default;

  template <typename T1, typename T2, typename T3>
  void init(T1& confFilter, T2& confTable, T3& initContext)
  {
    mMassMin = confFilter.massMin.value;
    mMassMax = confFilter.massMax.value;
    mPtMin = confFilter.ptMin.value;
    mPtMax = confFilter.ptMax.value;
    mEtaMin = confFilter.etaMin.value;
    mEtaMax = confFilter.etaMax.value;
    mPhiMin = confFilter.phiMin.value;
    mPhiMax = confFilter.phiMax.value;

    if constexpr (modes::isEqual(resoType, modes::TwoPhotonResonance::kPi0)) {
      LOG(info) << "Initialize femto Pi0 builder...";
      mProducePi0s = utils::enableTable("FPi0s_001", confTable.producePi0s.value, initContext);
      mProducePi0Masks = utils::enableTable("FPi0Masks_001", confTable.producePi0Masks.value, initContext);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoPhotonResonance::kEta)) {
      LOG(info) << "Initialize femto Eta builder...";
      mProduceEtas = utils::enableTable("FEtas_001", confTable.produceEtas.value, initContext);
      mProduceEtaMasks = utils::enableTable("FEtaMasks_001", confTable.produceEtaMasks.value, initContext);
    }

    if (mProducePi0s || mProducePi0Masks || mProduceEtas || mProduceEtaMasks) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillResonances(T1 const& col, T2& resonanceProducts, T3& photonPartition, T4& cache)
  {
    if (!mFillAnyTable) {
      return;
    }
    auto photonSlice = photonPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& [dau1, dau2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(photonSlice, photonSlice))) {
      this->fillResonance(col, dau1, dau2, resonanceProducts);
    }
  }

 private:
  template <typename T1, typename T2>
  void reconstructResonance(T1 const& dau1, T2 const& dau2)
  {
    ROOT::Math::PtEtaPhiMVector vecDau1{dau1.pt(), dau1.eta(), dau1.phi(), 0.f};
    ROOT::Math::PtEtaPhiMVector vecDau2{dau2.pt(), dau2.eta(), dau2.phi(), 0.f};
    ROOT::Math::PtEtaPhiMVector vecResonance = vecDau1 + vecDau2;

    mPt = vecResonance.Pt();
    mEta = vecResonance.Eta();
    mPhi = RecoDecay::constrainAngle(vecResonance.Phi());
    mMass = vecResonance.M();
  }

  [[nodiscard]] bool checkFilters() const
  {
    return ((mMass > mMassMin && mMass < mMassMax) &&
            (mPt > mPtMin && mPt < mPtMax) &&
            (mEta > mEtaMin && mEta < mEtaMax) &&
            (mPhi >= mPhiMin && mPhi < mPhiMax));
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillResonance(T1 const& col, T2 const& dau1, T3 const& dau2, T4& resonanceProducts)
  {
    reconstructResonance(dau1, dau2);
    if (!checkFilters()) {
      return;
    }

    if constexpr (modes::isEqual(resoType, modes::TwoPhotonResonance::kPi0)) {
      if (mProducePi0s) {
        resonanceProducts.producedPi0s(col.globalIndex(),
                                       mPt,
                                       mEta,
                                       mPhi,
                                       mMass,
                                       dau1.globalIndex(),
                                       dau2.globalIndex());
      }
      if (mProducePi0Masks) {
        resonanceProducts.producedPi0Masks(dau1.mask(), dau2.mask());
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoPhotonResonance::kEta)) {
      if (mProduceEtas) {
        resonanceProducts.producedEtas(col.globalIndex(),
                                       mPt,
                                       mEta,
                                       mPhi,
                                       mMass,
                                       dau1.globalIndex(),
                                       dau2.globalIndex());
      }
      if (mProduceEtaMasks) {
        resonanceProducts.producedEtaMasks(dau1.mask(), dau2.mask());
      }
    }
  }

  // cached kinematics of the resonance
  float mPt = 0;
  float mEta = 0;
  float mPhi = 0;
  float mMass = 0;

  float mMassMin = 0.f;
  float mMassMax = 0.f;
  float mPtMin = 0.f;
  float mPtMax = 0.f;
  float mEtaMin = 0.f;
  float mEtaMax = 0.f;
  float mPhiMin = 0.f;
  float mPhiMax = 0.f;

  bool mFillAnyTable = false;
  bool mProducePi0s = false;
  bool mProducePi0Masks = false;
  bool mProduceEtas = false;
  bool mProduceEtaMasks = false;
};

} // namespace o2::analysis::femto::twophotonresonancebuilder
#endif // PWGCF_FEMTO_CORE_TWOPHOTONRESONANCEBUILDER_H_
