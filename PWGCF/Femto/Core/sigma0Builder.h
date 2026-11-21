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

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto::sigma0builder
{
// Filters applied before the Sigma0 selection bits are evaluated.
struct ConfSigma0Filters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Filters");

  // Photon preselection.
  o2::framework::Configurable<datatypes::PhotonMaskType> photonMask{"photonMask", 0, "Photon selection bits that must all be set. 0: no requirement"};
  o2::framework::Configurable<float> photonPtMin{"photonPtMin", 0.f, "Minimum photon pT"};
  o2::framework::Configurable<float> photonRapMax{"photonRapMax", 0.5f, "Maximum photon |y| (equal to |eta| for a massless photon)"};

  // Lambda preselection.
  o2::framework::Configurable<datatypes::V0MaskType> lambdaMask{"lambdaMask", 0, "Lambda selection bits that must all be set. 0: no requirement"};
  o2::framework::Configurable<int> lambdaSign{"lambdaSign", 0, "Lambda sign (+1: Lambda; -1: AntiLambda; 0: both)"};
  o2::framework::Configurable<float> lambdaPtMin{"lambdaPtMin", 0.f, "Minimum Lambda pT"};
  o2::framework::Configurable<float> lambdaRapMax{"lambdaRapMax", 0.5f, "Maximum Lambda |y|, computed with the PDG Lambda mass"};
  o2::framework::Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.015f, "Maximum |m - m_Lambda(PDG)| in GeV/c^2. Set <= 0 to disable"};

  // Loose filters on the Lambda-photon pair.
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum Sigma0 pT"};
  o2::framework::Configurable<float> rapMax{"rapMax", 1.f, "Maximum Sigma0 |y|"};
  o2::framework::Configurable<float> massMin{"massMin", 1.1f, "Minimum Sigma0 invariant mass"};
  o2::framework::Configurable<float> massMax{"massMax", 1.3f, "Maximum Sigma0 invariant mass"};
};

// Working points stored in the Sigma0 selection mask.
struct ConfSigma0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Bits");
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, store all Sigma0 candidates passing the loose filters and record every selection bit"};
  o2::framework::Configurable<std::vector<float>> massWindow{"massWindow", {0.1f}, "Maximum |m(Lambda gamma) - m_Sigma0(PDG)| (GeV/c^2)"};
  o2::framework::Configurable<std::vector<float>> absRapMax{"absRapMax", {0.8f}, "Maximum Sigma0 |y|"};
  o2::framework::Configurable<std::vector<float>> ptMin{"ptMin", {0.f}, "Minimum Sigma0 pT (GeV/c)"};
};

// Base selection configuration for downstream analysis tasks.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                           \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "PDG code. Set sign to -1 for antiparticle"}; \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                     \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                   \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                                \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                                 \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                  \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                     \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass"};                        \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass"};                        \
  o2::framework::Configurable<datatypes::Sigma0MaskType> mask{"mask", 0, "Bitmask for Sigma0 selection"};

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

/// The selection bits stored for Sigma0 candidates.
enum Sigma0Sels {
  kSigmaMassWindow, ///< Maximum |m(Lambda gamma) - m_Sigma0(PDG)|
  kSigmaAbsRapMax,  ///< Maximum absolute Sigma0 rapidity
  kSigmaPtMin,      ///< Minimum Sigma0 pT
  kSigma0SelsMax
};

constexpr char Sigma0SelHistName[] = "hSigma0Selection";
constexpr char Sigma0SelsName[] = "Sigma0 selection object";
const std::unordered_map<Sigma0Sels, std::string> sigma0SelectionNames = {
  {kSigmaMassWindow, "Maximum Sigma0 mass deviation"},
  {kSigmaAbsRapMax, "Maximum absolute Sigma0 rapidity"},
  {kSigmaPtMin, "Minimum Sigma0 pT"}};

// Loose preselection of each Lambda-photon pair.
enum Sigma0Filters {
  kPtMin,
  kRapMax,
  kMassMin,
  kMassMax,
  kSigma0FiltersMax
};

constexpr char Sigma0FilterHistName[] = "hSigma0Filters";
const std::unordered_map<Sigma0Filters, std::string> sigma0FilterNames = {
  {kPtMin, "ptMin"},
  {kRapMax, "rapMax"},
  {kMassMin, "massMin"},
  {kMassMax, "massMax"}};

/// Cut class containing and executing all Sigma0 selection bits and loose filters.
class Sigma0Selection : public baseselection::BaseSelection<float, datatypes::Sigma0MaskType, kSigma0SelsMax>
{
 public:
  Sigma0Selection() = default;
  ~Sigma0Selection() override = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1& config, T2& filter)
  {
    this->init(config.passThrough.value);

    mPtMin = filter.ptMin.value;
    mRapMax = filter.rapMax.value;
    mMassMin = filter.massMin.value;
    mMassMax = filter.massMax.value;

    this->addSelection(kSigmaMassWindow, sigma0SelectionNames.at(kSigmaMassWindow), config.massWindow.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kSigmaAbsRapMax, sigma0SelectionNames.at(kSigmaAbsRapMax), config.absRapMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kSigmaPtMin, sigma0SelectionNames.at(kSigmaPtMin), config.ptMin.value, limits::kLowerLimit, true, true, false);

    this->setupSelectionHistogram<Sigma0SelHistName>(registry);
    this->template setupFilterHistogram<Sigma0FilterHistName>(
      registry,
      {
        {sigma0FilterNames.at(kPtMin), mPtMin},
        {sigma0FilterNames.at(kRapMax), mRapMax},
        {sigma0FilterNames.at(kMassMin), mMassMin},
        {sigma0FilterNames.at(kMassMax), mMassMax},
      });
  }

  /// Evaluate the selection bits for a Lambda-photon pair.
  void applySelections(float mass, float rapidity, float pt)
  {
    this->reset();
    this->evaluateObservable(kSigmaMassWindow, mass - o2::constants::physics::MassSigma0);
    this->evaluateObservable(kSigmaAbsRapMax, rapidity);
    this->evaluateObservable(kSigmaPtMin, pt);
    this->assembleBitmask<Sigma0SelHistName>();
  }

  /// Apply the loose preselection to a Lambda-photon pair.
  bool checkFilters(float mass, float rapidity, float pt) const
  {
    bool pass = true;
    bool p = false;

    p = pt > mPtMin;
    this->template fillFilter<Sigma0FilterHistName>(kPtMin, p);
    pass &= p;

    p = std::abs(rapidity) < mRapMax;
    this->template fillFilter<Sigma0FilterHistName>(kRapMax, p);
    pass &= p;

    p = mass > mMassMin;
    this->template fillFilter<Sigma0FilterHistName>(kMassMin, p);
    pass &= p;

    p = mass < mMassMax;
    this->template fillFilter<Sigma0FilterHistName>(kMassMax, p);
    pass &= p;

    this->template fillFilterSummary<Sigma0FilterHistName>(pass);
    return this->isPassThrough() || pass;
  }

 protected:
  float mPtMin = 0.f;
  float mRapMax = 1.f;
  float mMassMin = 1.1f;
  float mMassMax = 1.3f;
};

struct Sigma0BuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FSigma0s> producedSigma0s;
  o2::framework::Produces<o2::aod::FLiteSigma0s> producedLiteSigma0s;
  o2::framework::Produces<o2::aod::FSigma0Masks> producedSigma0Masks;
};

struct ConfSigma0Tables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Tables");
  o2::framework::Configurable<int> produceSigma0s{"produceSigma0s", -1, "Produce Sigma0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteSigma0s{"produceLiteSigma0s", -1, "Produce LiteSigma0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigma0Masks{"produceSigma0Masks", -1, "Produce Sigma0Masks (-1: auto; 0 off; 1 on)"};
};

class Sigma0Builder
{
 public:
  Sigma0Builder() = default;
  ~Sigma0Builder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry, T1& filter, T2& bits, T3& table, T4& initContext)
  {
    LOG(info) << "Initialize femto Sigma0 builder...";
    mRegistry = registry;

    mProduceSigma0s = utils::enableTable("FSigma0s_001", table.produceSigma0s.value, initContext);
    mProduceLiteSigma0s = utils::enableTable("FLiteSigma0s_001", table.produceLiteSigma0s.value, initContext);
    mProduceSigma0Masks = utils::enableTable("FSigma0Masks_002", table.produceSigma0Masks.value, initContext);

    if (mProduceSigma0s && mProduceLiteSigma0s) {
      LOG(fatal) << "FSigma0s and FLiteSigma0s are mutually exclusive -- enable only one. "
                 << "FLiteSigma0s is meant to replace FSigma0s at the producer stage for better compression in derived data.";
    }

    mFillAnyTable = mProduceSigma0s || mProduceLiteSigma0s || mProduceSigma0Masks;
    if (!mFillAnyTable) {
      LOG(info) << "No tables configured, Sigma0 selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }

    mPhotonMask = filter.photonMask.value;
    mPhotonPtMin = filter.photonPtMin.value;
    mPhotonRapMax = filter.photonRapMax.value;

    mLambdaMask = filter.lambdaMask.value;
    mLambdaSign = filter.lambdaSign.value;
    mLambdaPtMin = filter.lambdaPtMin.value;
    mLambdaRapMax = filter.lambdaRapMax.value;
    mLambdaMassWindow = filter.lambdaMassWindow.value;

    mSigma0Selection.configure(registry, bits, filter);
    mSigma0Selection.printSelections(Sigma0SelsName);

    // Daughter cut-flow histograms. Bin N is filled when the candidate passes cut N.
    mRegistry->add("hPhotonSelection", "Photon selection;;Entries", o2::framework::HistType::kTH1D, {{4, 0.5, 4.5}});
    setBinLabels(mRegistry->get<TH1>(HIST("hPhotonSelection")), {"No Sel", "Mask", "pT", "Y"});

    mRegistry->add("hLambdaSelection", "Lambda selection;;Entries", o2::framework::HistType::kTH1D, {{6, 0.5, 6.5}});
    setBinLabels(mRegistry->get<TH1>(HIST("hLambdaSelection")), {"No Sel", "Mask", "Sign", "pT", "Y", "Mass"});
    mRegistry->add("hLambdaMassBeforeSel", "Lambda mass before Lambda selections;#it{M}_{#Lambda} (GeV/#it{c}^{2});Entries", o2::framework::HistType::kTH1D, {{500, 1.10f, 1.30f}});
    mRegistry->add("hLambdaMassSelected", "Lambda mass after Lambda selections;#it{M}_{#Lambda} (GeV/#it{c}^{2});Entries", o2::framework::HistType::kTH1D, {{500, 1.10f, 1.30f}});
    mRegistry->add("hSigmaMassBeforeSel", "Sigma0 mass before selection;#it{M}_{#Lambda#gamma} (GeV/#it{c}^{2});Entries", o2::framework::HistType::kTH1D, {{500, 1.10f, 1.30f}});
    mRegistry->add("hSigmaMassSelected", "Sigma0 mass after selection;#it{M}_{#Lambda#gamma} (GeV/#it{c}^{2});Entries", o2::framework::HistType::kTH1D, {{500, 1.10f, 1.30f}});

    LOG(info) << "Initialization done...";
  }

  /// Build Sigma0 candidates from Lambda and photon tables belonging to the same collision.
  template <typename TProducts, typename TLambdas, typename TPhotons, typename TTracks>
  void fillSigma0s(TProducts& products, TLambdas const& lambdas, TPhotons const& photons, TTracks const& /*tracks*/)
  {
    if (!mFillAnyTable) {
      return;
    }

    // Preselect photons once and group their row indices by collision.
    std::unordered_map<int64_t, std::vector<int64_t>> photonRowsPerCollision;
    for (const auto& photon : photons) {
      if (selectPhoton(photon)) {
        photonRowsPerCollision[photon.fColId()].push_back(photon.globalIndex());
      }
    }
    if (photonRowsPerCollision.empty()) {
      return;
    }

    for (const auto& lambda : lambdas) {
      if (!selectLambda(lambda)) {
        continue;
      }
      // const auto posDau = lambda.posDau();
      // const auto negDau = lambda.negDau();
      const auto& posDau = lambda.template posDau_as<TTracks>();
      const auto& negDau = lambda.template negDau_as<TTracks>();

      const auto photonRows = photonRowsPerCollision.find(lambda.fColId());
      if (photonRows == photonRowsPerCollision.end()) {
        continue;
      }
      for (const auto photonRow : photonRows->second) {
        const auto photon = photons.rawIteratorAt(photonRow);
        buildSigma0(products, lambda, photon, posDau, negDau);
      }
    }
  }

 private:
  template <typename H>
  static void setBinLabels(H const& hist, std::initializer_list<const char*> labels)
  {
    int bin = 1;
    for (const char* label : labels) {
      hist->GetXaxis()->SetBinLabel(bin++, label);
    }
  }

  static std::array<float, 3> momentum(float pt, float eta, float phi)
  {
    return {pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta)};
  }

  template <typename T>
  bool selectPhoton(T const& photon)
  {
    mRegistry->fill(HIST("hPhotonSelection"), 1.);
    if ((photon.mask() & mPhotonMask) != mPhotonMask) {
      return false;
    }
    mRegistry->fill(HIST("hPhotonSelection"), 2.);
    if (photon.pt() < mPhotonPtMin) {
      return false;
    }
    mRegistry->fill(HIST("hPhotonSelection"), 3.);
    // For a massless photon, rapidity equals pseudorapidity.
    if (std::abs(photon.eta()) > mPhotonRapMax) {
      return false;
    }
    mRegistry->fill(HIST("hPhotonSelection"), 4.);
    return true;
  }

  template <typename T>
  bool selectLambda(T const& lambda)
  {
    mRegistry->fill(HIST("hLambdaMassBeforeSel"), lambda.mass());
    mRegistry->fill(HIST("hLambdaSelection"), 1.);
    if ((lambda.mask() & mLambdaMask) != mLambdaMask) {
      return false;
    }
    mRegistry->fill(HIST("hLambdaSelection"), 2.);

    const bool isAntiLambda = lambda.signedPt() < 0.f;
    if ((mLambdaSign > 0 && isAntiLambda) || (mLambdaSign < 0 && !isAntiLambda)) {
      return false;
    }
    mRegistry->fill(HIST("hLambdaSelection"), 3.);

    const float pt = std::abs(lambda.signedPt());
    if (pt < mLambdaPtMin) {
      return false;
    }
    mRegistry->fill(HIST("hLambdaSelection"), 4.);

    const auto pLambda = momentum(pt, lambda.eta(), lambda.phi());
    if (std::abs(RecoDecay::y(pLambda, o2::constants::physics::MassLambda0)) > mLambdaRapMax) {
      return false;
    }
    mRegistry->fill(HIST("hLambdaSelection"), 5.);

    if (mLambdaMassWindow > 0.f && std::abs(lambda.mass() - o2::constants::physics::MassLambda0) > mLambdaMassWindow) {
      return false;
    }
    mRegistry->fill(HIST("hLambdaSelection"), 6.);
    mRegistry->fill(HIST("hLambdaMassSelected"), lambda.mass());
    return true;
  }

  template <typename TProducts, typename TLambda, typename TPhoton, typename TPosDau, typename TNegDau>
  bool buildSigma0(TProducts& products, TLambda const& lambda, TPhoton const& photon, TPosDau const& posDau, TNegDau const& negDau)
  {
    const float lambdaPt = std::abs(lambda.signedPt());
    const std::array<float, 3> pLambda = momentum(lambdaPt, lambda.eta(), lambda.phi());
    const std::array<float, 3> pPhoton = momentum(photon.pt(), photon.eta(), photon.phi());
    const std::array<float, 3> pSigma = {pLambda[0] + pPhoton[0], pLambda[1] + pPhoton[1], pLambda[2] + pPhoton[2]};

    const float sigmaPt = std::hypot(pSigma[0], pSigma[1]);
    if (sigmaPt <= 0.f) {
      return false;
    }

    // Compute the pair mass with the PDG masses of Lambda and photon.
    const float sigmaMass = RecoDecay::m(std::array{pPhoton, pLambda}, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    const float sigmaY = RecoDecay::y(pSigma, o2::constants::physics::MassSigma0);
    mRegistry->fill(HIST("hSigmaMassBeforeSel"), sigmaMass);

    if (!mSigma0Selection.checkFilters(sigmaMass, sigmaY, sigmaPt)) {
      return false;
    }
    mSigma0Selection.applySelections(sigmaMass, sigmaY, sigmaPt);
    if (!mSigma0Selection.passesAllRequiredSelections()) {
      return false;
    }
    mRegistry->fill(HIST("hSigmaMassSelected"), sigmaMass);

    const float sigmaEta = std::asinh(pSigma[2] / sigmaPt);
    const float sigmaPhi = RecoDecay::phi(pSigma[0], pSigma[1]);
    const float sign = lambda.signedPt() < 0.f ? -1.f : 1.f;

    if (mProduceSigma0s) {
      products.producedSigma0s(lambda.fColId(),
                               sign * sigmaPt,
                               sigmaEta,
                               sigmaPhi,
                               sigmaMass,
                               static_cast<int32_t>(lambda.globalIndex()),
                               static_cast<int32_t>(photon.globalIndex()),
                               static_cast<int32_t>(posDau.globalIndex()),
                               static_cast<int32_t>(negDau.globalIndex()));
    }
    if (mProduceLiteSigma0s) {
      products.producedLiteSigma0s(lambda.fColId(),
                                   o2::aod::femtobase::lite::binSignedPt(sign * sigmaPt),
                                   o2::aod::femtobase::lite::binEta(sigmaEta),
                                   o2::aod::femtobase::lite::binPhi(sigmaPhi),
                                   o2::aod::femtosigma0s::lite::binSigma0Mass(sigmaMass),
                                   static_cast<int32_t>(lambda.globalIndex()),
                                   static_cast<int32_t>(photon.globalIndex()));
    }
    if (mProduceSigma0Masks) {
      products.producedSigma0Masks(mSigma0Selection.getBitmask());
    }
    return true;
  }

  o2::framework::HistogramRegistry* mRegistry = nullptr;
  Sigma0Selection mSigma0Selection;
  bool mFillAnyTable = false;
  bool mProduceSigma0s = false;
  bool mProduceLiteSigma0s = false;
  bool mProduceSigma0Masks = false;

  // Photon preselection.
  datatypes::PhotonMaskType mPhotonMask = 0;
  float mPhotonPtMin = 0.f;
  float mPhotonRapMax = 0.5f;

  // Lambda preselection.
  datatypes::V0MaskType mLambdaMask = 0;
  int mLambdaSign = 0;
  float mLambdaPtMin = 0.f;
  float mLambdaRapMax = 0.5f;
  float mLambdaMassWindow = 0.015f;
};

} // namespace o2::analysis::femto::sigma0builder
#endif // PWGCF_FEMTO_CORE_SIGMA0BUILDER_H_
