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

/// \file photonBuilder.h
/// \brief photon (PCM) builder
///
/// Default selection values for the mandatory cuts (dauTpcNSigmaElAbsMax, keepDaughtersWithoutTpc)
/// are taken from the PCM branch of EventFiltering/PWGEM/EMPhotonFilter.cxx (isSelectedSecondary()
/// and the min_pt_pcm_photon cut in runFilter()): both V0 daughters must have |TPC nSigma_e| < 3.5
/// whenever a TPC signal is available, and a daughter without TPC is not rejected. Note that
/// EMPhotonFilter.cxx also declares minpt_v0, maxeta_v0, max_dcatopv_xy_v0 and max_dcatopv_z_v0,
/// but never applies them -- they are dead configurables in that task. The additional geometric
/// cuts here (cosPA, DCA to PV, chi2/NDF, radius) are standard PCM conversion-quality cuts, not
/// used by that filter task, so they default to permissive (no-op) values.
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PHOTONBUILDER_H_
#define PWGCF_FEMTO_CORE_PHOTONBUILDER_H_

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto::photonbuilder
{

// filters applied in the producer task (loose pre-selection, applied before quality/PID cuts)
struct ConfPhotonFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhotonFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT (matches EMPhotonFilter min_pt_pcm_photon default)"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> mGammaAbsMax{"mGammaAbsMax", 1e10f, "Maximum |m_ee| of the conversion (KF) at the secondary vertex"};
};

// selection bits for photons (PCM)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PHOTON_DEFAULT_BITS                                                                                                                                                       \
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all photons are passed through. Bits for all selections are stored."};                          \
  o2::framework::Configurable<std::vector<float>> cosPaMin{"cosPaMin", {-1.f}, "Minimum cosine of pointing angle (not used by EMPhotonFilter, default is permissive)"};           \
  o2::framework::Configurable<std::vector<float>> dcaToPvXYAbsMax{"dcaToPvXYAbsMax", {1e10f}, "Maximum |DCAxy| of the photon (V0) to the PV (cm)"};                               \
  o2::framework::Configurable<std::vector<float>> dcaToPvZAbsMax{"dcaToPvZAbsMax", {1e10f}, "Maximum |DCAz| of the photon (V0) to the PV (cm)"};                                  \
  o2::framework::Configurable<std::vector<float>> chi2NdfMax{"chi2NdfMax", {1e10f}, "Maximum chi2/NDF of the KF conversion vertex"};                                              \
  o2::framework::Configurable<std::vector<float>> v0RadiusMin{"v0RadiusMin", {0.f}, "Minimum transverse radius of the conversion point (cm)"};                                   \
  o2::framework::Configurable<std::vector<float>> v0RadiusMax{"v0RadiusMax", {1e10f}, "Maximum transverse radius of the conversion point (cm)"};                                 \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.9f}, "Maximum |eta| for daughter tracks (matches EMPhotonFilter maxeta_v0, though unused there)"}; \
  o2::framework::Configurable<std::vector<float>> dauTpcNSigmaElAbsMax{"dauTpcNSigmaElAbsMax", {3.5f}, "Maximum |TPC nSigma_e| for V0 daughters (from EMPhotonFilter isSelectedSecondary)"}; \
  o2::framework::Configurable<bool> keepDaughtersWithoutTpc{"keepDaughtersWithoutTpc", true, "If true, a daughter without a TPC signal is not rejected by the TPC electron PID cut (matches EMPhotonFilter behaviour)"};

struct ConfPhotonBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhotonBits");
  PHOTON_DEFAULT_BITS
};

#undef PHOTON_DEFAULT_BITS

// base selection for QA/analysis tasks, used to re-select already produced photons by mask
// (analogous to ConfLambdaSelection/ConfK0shortSelection in spirit, but photons have neither
// a PDG hypothesis nor a mass window)
template <auto& Prefix>
struct ConfPhotonSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<datatypes::PhotonMaskType> mask{"mask", 0, "Bitmask for photon selection"};
};

constexpr const char PrefixPhotonSelection1[] = "PhotonSelection1";
using ConfPhotonSelection1 = ConfPhotonSelection<PrefixPhotonSelection1>;

/// The different selections for photons (PCM)
enum PhotonSels {
  kCosPaMin,         ///< Min. CPA (cosine pointing angle) of the conversion to the PV
  kDcaToPvXYAbsMax,  ///< Max. |DCAxy| of the photon to the PV
  kDcaToPvZAbsMax,   ///< Max. |DCAz| of the photon to the PV
  kChi2NdfMax,       ///< Max. chi2/NDF of the KF conversion vertex
  kV0RadiusMin,      ///< Min. transverse radius of the conversion point
  kV0RadiusMax,      ///< Max. transverse radius of the conversion point
  kDauAbsEtaMax,     ///< Max. absolute pseudorapidity of the daughters
  kPosDauTpcNSigmaEl, ///< TPC electron PID for positive daughter
  kNegDauTpcNSigmaEl, ///< TPC electron PID for negative daughter
  kPhotonSelsMax
};

constexpr char PhotonSelHistName[] = "hPhotonSelection";
constexpr char PhotonSelsName[] = "Photon selection object";
const std::unordered_map<PhotonSels, std::string> photonSelectionNames = {
  {kCosPaMin, "Min. CPA (cosine pointing angle)"},
  {kDcaToPvXYAbsMax, "Max. |DCAxy| to PV"},
  {kDcaToPvZAbsMax, "Max. |DCAz| to PV"},
  {kChi2NdfMax, "Max. chi2/NDF of conversion vertex"},
  {kV0RadiusMin, "Min. transverse radius"},
  {kV0RadiusMax, "Max. transverse radius"},
  {kDauAbsEtaMax, "Max. absolute pseudorapidity of daughters"},
  {kPosDauTpcNSigmaEl, "TPC electron PID for positive daughter"},
  {kNegDauTpcNSigmaEl, "TPC electron PID for negative daughter"},
};

// enum for all photon (pre)filters (loose pre-selection, applied before quality/PID cuts)
enum PhotonFilters {
  kPtMin,
  kPtMax,
  kEtaMin,
  kEtaMax,
  kPhiMin,
  kPhiMax,
  kMGammaAbsMax,
  kPhotonFiltersMax
};

constexpr char PhotonFilterHistName[] = "hPhotonFilters";
const std::unordered_map<PhotonFilters, std::string> photonFilterNames = {
  {kPtMin, "ptMin"},
  {kPtMax, "ptMax"},
  {kEtaMin, "etaMin"},
  {kEtaMax, "etaMax"},
  {kPhiMin, "phiMin"},
  {kPhiMax, "phiMax"},
  {kMGammaAbsMax, "mGammaAbsMax"},
};

/// \brief Cut class to contain and execute all cuts applied to PCM photons (V0PhotonsKF)
template <auto& SelectionHistName, auto& FilterHistName>
class PhotonSelection : public baseselection::BaseSelection<float, datatypes::PhotonMaskType, kPhotonSelsMax>
{
 public:
  PhotonSelection() = default;
  ~PhotonSelection() override = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1& config, T2& filter)
  {
    this->init(config.passThrough.value);

    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;
    mMGammaAbsMax = filter.mGammaAbsMax.value;
    mKeepDaughtersWithoutTpc = config.keepDaughtersWithoutTpc.value;

    this->addSelection(kCosPaMin, photonSelectionNames.at(kCosPaMin), config.cosPaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kDcaToPvXYAbsMax, photonSelectionNames.at(kDcaToPvXYAbsMax), config.dcaToPvXYAbsMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kDcaToPvZAbsMax, photonSelectionNames.at(kDcaToPvZAbsMax), config.dcaToPvZAbsMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kChi2NdfMax, photonSelectionNames.at(kChi2NdfMax), config.chi2NdfMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kV0RadiusMin, photonSelectionNames.at(kV0RadiusMin), config.v0RadiusMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kV0RadiusMax, photonSelectionNames.at(kV0RadiusMax), config.v0RadiusMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kDauAbsEtaMax, photonSelectionNames.at(kDauAbsEtaMax), config.dauAbsEtaMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kPosDauTpcNSigmaEl, photonSelectionNames.at(kPosDauTpcNSigmaEl), config.dauTpcNSigmaElAbsMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kNegDauTpcNSigmaEl, photonSelectionNames.at(kNegDauTpcNSigmaEl), config.dauTpcNSigmaElAbsMax.value, limits::kAbsUpperLimit, true, true, false);

    this->setupSelectionHistogram<SelectionHistName>(registry);
    this->template setupFilterHistogram<FilterHistName>(
      registry,
      {
        {photonFilterNames.at(kPtMin), mPtMin},
        {photonFilterNames.at(kPtMax), mPtMax},
        {photonFilterNames.at(kEtaMin), mEtaMin},
        {photonFilterNames.at(kEtaMax), mEtaMax},
        {photonFilterNames.at(kPhiMin), mPhiMin},
        {photonFilterNames.at(kPhiMax), mPhiMax},
        {photonFilterNames.at(kMGammaAbsMax), mMGammaAbsMax},
      });
  }

  template <typename T1, typename T2>
  void applySelections(T1 const& photon, T2 const& /*v0legs*/)
  {
    this->reset();
    this->evaluateObservable(kCosPaMin, photon.cospa());
    this->evaluateObservable(kDcaToPvXYAbsMax, photon.dcaXYtopv());
    this->evaluateObservable(kDcaToPvZAbsMax, photon.dcaZtopv());
    this->evaluateObservable(kChi2NdfMax, photon.chiSquareNDF());
    this->evaluateObservable(kV0RadiusMin, photon.v0radius());
    this->evaluateObservable(kV0RadiusMax, photon.v0radius());

    auto posDaughter = photon.template posTrack_as<T2>();
    auto negDaughter = photon.template negTrack_as<T2>();

    std::array<float, 2> etaAbsDaughters = {std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaAbsDaughters.begin(), etaAbsDaughters.end()));

    // TPC electron PID: evaluate only when a TPC signal is available, matching
    // EMPhotonFilter::isSelectedSecondary(), which never rejects a daughter without TPC
    auto evaluateDaughterTpcEl = [this](PhotonSels bit, bool hasTpc, float tpcNSigmaEl) {
      if (hasTpc) {
        this->evaluateObservable(bit, tpcNSigmaEl);
      } else if (mKeepDaughtersWithoutTpc) {
        this->evaluateObservable(bit, 0.f);
      }
    };
    evaluateDaughterTpcEl(kPosDauTpcNSigmaEl, posDaughter.hasTPC(), posDaughter.tpcNSigmaEl());
    evaluateDaughterTpcEl(kNegDauTpcNSigmaEl, negDaughter.hasTPC(), negDaughter.tpcNSigmaEl());

    this->assembleBitmask<SelectionHistName>();
  }

  template <typename T>
  bool checkFilters(const T& photon) const
  {
    bool pass = true;
    bool p = false;

    p = photon.pt() > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = photon.pt() < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = photon.eta() > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = photon.eta() < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = photon.phi() > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = photon.phi() < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    p = std::fabs(photon.mGamma()) < mMGammaAbsMax;
    this->template fillFilter<FilterHistName>(kMGammaAbsMax, p);
    pass &= p;

    this->template fillFilterSummary<FilterHistName>(pass);

    return this->isPassThrough() || pass;
  }

 protected:
  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 99.f;
  float mEtaMin = -10.f;
  float mEtaMax = 10.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
  float mMGammaAbsMax = 1e10f;

  bool mKeepDaughtersWithoutTpc = true;
};

struct PhotonBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FPhotons> producedPhotons;
  o2::framework::Produces<o2::aod::FLitePhotons> producedLitePhotons;
  o2::framework::Produces<o2::aod::FPhotonMasks> producedPhotonMasks;
  o2::framework::Produces<o2::aod::FPhotonExtras> producedPhotonExtras;
};

struct ConfPhotonTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhotonTables");
  o2::framework::Configurable<int> producePhotons{"producePhotons", -1, "Produce Photons (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLitePhotons{"produceLitePhotons", -1, "Produce LitePhotons (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePhotonMasks{"producePhotonMasks", -1, "Produce PhotonMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePhotonExtras{"producePhotonExtras", -1, "Produce PhotonExtras (-1: auto; 0 off; 1 on)"};
};

/// \brief Builder for PCM photons (V0PhotonsKF + V0Legs, from PWGEM/PhotonMeson)
///
/// Unlike the strangeness V0 builder, daughter tracks are not indexed into FTracks: the
/// conversion legs (V0Legs) carry a reduced set of columns and are not meant to be paired
/// individually in femtoscopy, so only photon-level kinematics plus daughter PID/quality
/// diagnostics (FPhotonExtras) are stored.
class PhotonBuilder
{
 public:
  PhotonBuilder() = default;
  ~PhotonBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext)
  {
    LOG(info) << "Initialize femto photon (PCM) builder...";
    mProducePhotons = utils::enableTable("FPhotons_001", table.producePhotons.value, initContext);
    mProduceLitePhotons = utils::enableTable("FLitePhotons_001", table.produceLitePhotons.value, initContext);
    mProducePhotonMasks = utils::enableTable("FPhotonMasks_001", table.producePhotonMasks.value, initContext);
    mProducePhotonExtras = utils::enableTable("FPhotonExtras_001", table.producePhotonExtras.value, initContext);

    if (mProducePhotons && mProduceLitePhotons) {
      LOG(fatal) << "FPhotons and FLitePhotons are mutually exclusive -- enable only one. "
                 << "FLitePhotons is meant to replace FPhotons at the producer stage (for better compression in derived data).";
    }
    if (mProducePhotons || mProduceLitePhotons || mProducePhotonMasks || mProducePhotonExtras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mPhotonSelection.configure(registry, config, filter);
    mPhotonSelection.printSelections(PhotonSelsName);
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fillPhotons(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& photonProducts, T5 const& photons, T6 const& v0legs)
  {
    if (!mFillAnyTable) {
      return;
    }
    for (const auto& photon : photons) {
      if (!mPhotonSelection.checkFilters(photon)) {
        continue;
      }
      mPhotonSelection.applySelections(photon, v0legs);
      if (!mPhotonSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);

      if (mProducePhotons) {
        photonProducts.producedPhotons(collisionBuilder.collisionIndex(),
                                       photon.pt(),
                                       photon.eta(),
                                       photon.phi());
      }
      if (mProduceLitePhotons) {
        photonProducts.producedLitePhotons(collisionBuilder.collisionIndex(),
                                           o2::aod::femtobase::lite::binUnsignedPt(photon.pt()),
                                           o2::aod::femtobase::lite::binEta(photon.eta()),
                                           o2::aod::femtobase::lite::binPhi(photon.phi()));
      }
      if (mProducePhotonMasks) {
        photonProducts.producedPhotonMasks(mPhotonSelection.getBitmask());
      }
      if (mProducePhotonExtras) {
        auto posDaughter = photon.template posTrack_as<T6>();
        auto negDaughter = photon.template negTrack_as<T6>();
        photonProducts.producedPhotonExtras(
          photon.cospa(),
          photon.dcaXYtopv(),
          photon.dcaZtopv(),
          photon.chiSquareNDF(),
          photon.v0radius(),
          posDaughter.pt(),
          negDaughter.pt(),
          posDaughter.tpcInnerParam(),
          negDaughter.tpcInnerParam(),
          posDaughter.tpcSignal(),
          negDaughter.tpcSignal(),
          posDaughter.tpcNSigmaEl(),
          negDaughter.tpcNSigmaEl(),
          photon.vx(),
          photon.vy(),
          photon.vz());
      }
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mPhotonSelection.isPassThrough(); }

 private:
  PhotonSelection<PhotonSelHistName, PhotonFilterHistName> mPhotonSelection;
  bool mFillAnyTable = false;
  bool mProducePhotons = false;
  bool mProduceLitePhotons = false;
  bool mProducePhotonMasks = false;
  bool mProducePhotonExtras = false;
};

} // namespace o2::analysis::femto::photonbuilder
#endif // PWGCF_FEMTO_CORE_PHOTONBUILDER_H_
