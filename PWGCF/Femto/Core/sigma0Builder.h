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
#include "PWGCF/Femto/Core/photonBuilder.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace o2::analysis::femto::sigma0builder
{

// Selection of the input Lambda and photon daughters. Applied with partitions in the producer, before the builder.
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

// Bits of the output daughter masks stored in FSigma0Masks.
// Each entry is a (decimal) bitmask on the daughter mask of the femto producer (see its hLambdaSelection/hPhotonSelection) and defines one
// output bit, which is set if the daughter has all bits of the entry set. The input mask of the partition (Sigma0LambdaSelection.mask,
// Sigma0PhotonSelection.mask) is always added as the loosest selection, without a bit.
struct ConfSigma0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Bits");
  o2::framework::Configurable<std::vector<int>> lambdaMask{"lambdaMask", {}, "Bitmasks (decimal) on the Lambda mask, one bit of FSigma0Masks.maskLambda per entry"};
  o2::framework::Configurable<std::vector<int>> photonMask{"photonMask", {}, "Bitmasks (decimal) on the photon mask, one bit of FSigma0Masks.maskPhoton per entry"};
};

// Minimal QA histograms of the producer (independent of FSigma0Extras)
struct ConfSigma0Qa : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Qa");
  o2::framework::Configurable<bool> fillQa{"fillQa", true, "Fill QA histograms of the accepted Sigma0 candidates and their daughters"};
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "pT (GeV/c) of Sigma0, Lambda and photon"};
  o2::framework::ConfigurableAxis sigma0Mass{"sigma0Mass", {{400, 1.1, 1.3}}, "Invariant mass of the Sigma0 (GeV/c^2)"};
  o2::framework::ConfigurableAxis lambdaMass{"lambdaMass", {{400, 1.0, 1.2}}, "Invariant mass of the Lambda daughter (GeV/c^2)"};
};

/// \brief Histogram documenting a filter, i.e. the configured bounds and how many candidates pass each of them.
/// Same bin labels as BaseSelection::setupFilterHistogram, but the values are passed as strings so that bitmasks are printed exactly.
class FilterHistogram
{
 public:
  void init(o2::framework::HistogramRegistry* registry, std::string_view name, std::vector<std::pair<std::string, std::string>> const& filters)
  {
    mNFilters = static_cast<int>(filters.size());
    mHist = registry->add<TH1>(std::string(name), "; Filter; Entries", o2::framework::HistType::kTH1F, {{mNFilters + 2, -0.5, mNFilters + 1.5}});
    auto* axis = mHist->GetXaxis();
    for (int i = 0; i < mNFilters; ++i) {
      std::ostringstream label;
      label << "FilterName" << selectioncontainer::ValueDelimiter << filters.at(i).first << selectioncontainer::SectionDelimiter
            << "FilterValue" << selectioncontainer::ValueDelimiter << filters.at(i).second;
      axis->SetBinLabel(i + 1, label.str().c_str());
    }
    axis->SetBinLabel(mNFilters + 1, "All analyzed");
    axis->SetBinLabel(mNFilters + 2, "All passed");
  }

  /// fill the bin of one filter bound if it is passed, returns the decision
  bool fill(int filterIndex, bool passed)
  {
    if (passed) {
      mHist->Fill(filterIndex);
    }
    return passed;
  }

  /// fill the summary bins, called once per candidate
  void fillSummary(bool passedAll)
  {
    mHist->Fill(mNFilters);
    if (passedAll) {
      mHist->Fill(mNFilters + 1);
    }
  }

  template <typename T>
  static std::string toString(T value)
  {
    std::ostringstream oss;
    oss << value;
    return oss.str();
  }

  template <typename T>
  static std::string maskToString(T mask)
  {
    return std::to_string(static_cast<uint64_t>(mask));
  }

 private:
  std::shared_ptr<TH1> mHist;
  int mNFilters = 0;
};

// filters applied to the input collisions (Filter), Lambdas and photons (Partitions) of the producer
enum CollisionFilters {
  kColVtxZMin,
  kColVtxZMax,
  kColMultMin,
  kColMultMax,
  kColCentMin,
  kColCentMax,
  kColMagFieldMin,
  kColMagFieldMax,
  kColMask,
  kCollisionFiltersMax
};

enum LambdaFilters {
  kLambdaSign,
  kLambdaPtMin,
  kLambdaPtMax,
  kLambdaEtaMin,
  kLambdaEtaMax,
  kLambdaPhiMin,
  kLambdaPhiMax,
  kLambdaMassMin,
  kLambdaMassMax,
  kLambdaMask,
  kLambdaFiltersMax
};

enum PhotonFilters {
  kPhotonPtMin,
  kPhotonPtMax,
  kPhotonEtaMin,
  kPhotonEtaMax,
  kPhotonPhiMin,
  kPhotonPhiMax,
  kPhotonMask,
  kPhotonFiltersMax
};

// filters on the reconstructed Lambda-photon pair
enum Sigma0Filters {
  kSigma0PtMin,
  kSigma0PtMax,
  kSigma0EtaMin,
  kSigma0EtaMax,
  kSigma0PhiMin,
  kSigma0PhiMax,
  kSigma0MassMin,
  kSigma0MassMax,
  kSigma0FiltersMax
};

constexpr std::string_view CollisionFilterHistName = "hCollisionFilters";
constexpr std::string_view LambdaFilterHistName = "hLambdaFilters";
constexpr std::string_view PhotonFilterHistName = "hPhotonFilters";
constexpr std::string_view Sigma0FilterHistName = "hSigma0Filters";

// minimal QA histograms of the producer
constexpr char QaSigma0Mass[] = "Sigma0Qa/hSigma0Mass";
constexpr char QaLambdaMass[] = "Sigma0Qa/hLambdaMass";
constexpr char QaSigma0Pt[] = "Sigma0Qa/hSigma0Pt";
constexpr char QaLambdaPt[] = "Sigma0Qa/hLambdaPt";
constexpr char QaPhotonPt[] = "Sigma0Qa/hPhotonPt";

/// \brief Documents the filters applied to the input tables of the (derived-to-derived) Sigma0 producer.
/// The filters themselves are applied by the Filter/Partitions of the producer; here they are re-evaluated on the full input
/// to fill the filter histograms. Lambdas and photons are only counted in collisions passing the collision filter.
class Sigma0InputFilters
{
 public:
  template <typename T1, typename T2, typename T3>
  void init(o2::framework::HistogramRegistry* registry, T1 const& confCollision, T2 const& confLambda, T3 const& confPhoton)
  {
    mVtxZMin = confCollision.vtxZMin.value;
    mVtxZMax = confCollision.vtxZMax.value;
    mMultMin = confCollision.multMin.value;
    mMultMax = confCollision.multMax.value;
    mCentMin = confCollision.centMin.value;
    mCentMax = confCollision.centMax.value;
    mMagFieldMin = static_cast<int8_t>(confCollision.magFieldMin.value);
    mMagFieldMax = static_cast<int8_t>(confCollision.magFieldMax.value);
    mCollisionMask = confCollision.collisionMask.value;

    mLambdaSign = confLambda.sign.value;
    mLambdaPtMin = confLambda.ptMin.value;
    mLambdaPtMax = confLambda.ptMax.value;
    mLambdaEtaMin = confLambda.etaMin.value;
    mLambdaEtaMax = confLambda.etaMax.value;
    mLambdaPhiMin = confLambda.phiMin.value;
    mLambdaPhiMax = confLambda.phiMax.value;
    mLambdaMassMin = confLambda.massMin.value;
    mLambdaMassMax = confLambda.massMax.value;
    mLambdaMask = confLambda.mask.value;

    mPhotonPtMin = confPhoton.ptMin.value;
    mPhotonPtMax = confPhoton.ptMax.value;
    mPhotonEtaMin = confPhoton.etaMin.value;
    mPhotonEtaMax = confPhoton.etaMax.value;
    mPhotonPhiMin = confPhoton.phiMin.value;
    mPhotonPhiMax = confPhoton.phiMax.value;
    mPhotonMask = confPhoton.mask.value;

    mCollisionHist.init(registry, CollisionFilterHistName,
                        {
                          {"vtxZMin", FilterHistogram::toString(mVtxZMin)},
                          {"vtxZMax", FilterHistogram::toString(mVtxZMax)},
                          {"multMin", FilterHistogram::toString(mMultMin)},
                          {"multMax", FilterHistogram::toString(mMultMax)},
                          {"centMin", FilterHistogram::toString(mCentMin)},
                          {"centMax", FilterHistogram::toString(mCentMax)},
                          {"magFieldMin", FilterHistogram::toString(static_cast<int>(mMagFieldMin))},
                          {"magFieldMax", FilterHistogram::toString(static_cast<int>(mMagFieldMax))},
                          {"collisionMask", FilterHistogram::maskToString(mCollisionMask)},
                        });
    mLambdaHist.init(registry, LambdaFilterHistName,
                     {
                       {"sign", FilterHistogram::toString(mLambdaSign)},
                       {"ptMin", FilterHistogram::toString(mLambdaPtMin)},
                       {"ptMax", FilterHistogram::toString(mLambdaPtMax)},
                       {"etaMin", FilterHistogram::toString(mLambdaEtaMin)},
                       {"etaMax", FilterHistogram::toString(mLambdaEtaMax)},
                       {"phiMin", FilterHistogram::toString(mLambdaPhiMin)},
                       {"phiMax", FilterHistogram::toString(mLambdaPhiMax)},
                       {"massMin", FilterHistogram::toString(mLambdaMassMin)},
                       {"massMax", FilterHistogram::toString(mLambdaMassMax)},
                       {"mask", FilterHistogram::maskToString(mLambdaMask)},
                     });
    mPhotonHist.init(registry, PhotonFilterHistName,
                     {
                       {"ptMin", FilterHistogram::toString(mPhotonPtMin)},
                       {"ptMax", FilterHistogram::toString(mPhotonPtMax)},
                       {"etaMin", FilterHistogram::toString(mPhotonEtaMin)},
                       {"etaMax", FilterHistogram::toString(mPhotonEtaMax)},
                       {"phiMin", FilterHistogram::toString(mPhotonPhiMin)},
                       {"phiMax", FilterHistogram::toString(mPhotonPhiMax)},
                       {"mask", FilterHistogram::maskToString(mPhotonMask)},
                     });
  }

  template <typename T1, typename T2, typename T3>
  void fill(T1 const& cols, T2 const& lambdas, T3 const& photons)
  {
    mCollisionPassed.assign(cols.size(), false);
    for (auto const& col : cols) {
      mCollisionPassed.at(col.globalIndex()) = checkCollision(col);
    }
    for (auto const& lambda : lambdas) {
      if (mCollisionPassed.at(lambda.fColId())) {
        checkLambda(lambda);
      }
    }
    for (auto const& photon : photons) {
      if (mCollisionPassed.at(photon.fColId())) {
        checkPhoton(photon);
      }
    }
  }

 private:
  // same conditions as MAKE_COLLISION_FILTER
  template <typename T>
  bool checkCollision(T const& col)
  {
    bool pass = true;
    pass &= mCollisionHist.fill(kColVtxZMin, col.posZ() >= mVtxZMin);
    pass &= mCollisionHist.fill(kColVtxZMax, col.posZ() <= mVtxZMax);
    pass &= mCollisionHist.fill(kColMultMin, col.mult() >= mMultMin);
    pass &= mCollisionHist.fill(kColMultMax, col.mult() <= mMultMax);
    pass &= mCollisionHist.fill(kColCentMin, col.cent() >= mCentMin);
    pass &= mCollisionHist.fill(kColCentMax, col.cent() <= mCentMax);
    pass &= mCollisionHist.fill(kColMagFieldMin, col.magField() >= mMagFieldMin);
    pass &= mCollisionHist.fill(kColMagFieldMax, col.magField() <= mMagFieldMax);
    pass &= mCollisionHist.fill(kColMask, (col.mask() & mCollisionMask) == mCollisionMask);
    mCollisionHist.fillSummary(pass);
    return pass;
  }

  // same conditions as MAKE_LAMBDA_PARTITION
  template <typename T>
  void checkLambda(T const& lambda)
  {
    const float signedPt = lambda.signedPt();
    const float pt = std::abs(signedPt);
    bool pass = true;
    pass &= mLambdaHist.fill(kLambdaSign, mLambdaSign == 0 || (mLambdaSign > 0 ? signedPt > 0.f : signedPt < 0.f));
    pass &= mLambdaHist.fill(kLambdaPtMin, pt > mLambdaPtMin);
    pass &= mLambdaHist.fill(kLambdaPtMax, pt < mLambdaPtMax);
    pass &= mLambdaHist.fill(kLambdaEtaMin, lambda.eta() > mLambdaEtaMin);
    pass &= mLambdaHist.fill(kLambdaEtaMax, lambda.eta() < mLambdaEtaMax);
    pass &= mLambdaHist.fill(kLambdaPhiMin, lambda.phi() > mLambdaPhiMin);
    pass &= mLambdaHist.fill(kLambdaPhiMax, lambda.phi() < mLambdaPhiMax);
    pass &= mLambdaHist.fill(kLambdaMassMin, lambda.mass() > mLambdaMassMin);
    pass &= mLambdaHist.fill(kLambdaMassMax, lambda.mass() < mLambdaMassMax);
    pass &= mLambdaHist.fill(kLambdaMask, (lambda.mask() & mLambdaMask) == mLambdaMask);
    mLambdaHist.fillSummary(pass);
  }

  // same conditions as MAKE_PHOTON_PARTITION
  template <typename T>
  void checkPhoton(T const& photon)
  {
    bool pass = true;
    pass &= mPhotonHist.fill(kPhotonPtMin, photon.pt() > mPhotonPtMin);
    pass &= mPhotonHist.fill(kPhotonPtMax, photon.pt() < mPhotonPtMax);
    pass &= mPhotonHist.fill(kPhotonEtaMin, photon.eta() > mPhotonEtaMin);
    pass &= mPhotonHist.fill(kPhotonEtaMax, photon.eta() < mPhotonEtaMax);
    pass &= mPhotonHist.fill(kPhotonPhiMin, photon.phi() > mPhotonPhiMin);
    pass &= mPhotonHist.fill(kPhotonPhiMax, photon.phi() < mPhotonPhiMax);
    pass &= mPhotonHist.fill(kPhotonMask, (photon.mask() & mPhotonMask) == mPhotonMask);
    mPhotonHist.fillSummary(pass);
  }

  FilterHistogram mCollisionHist;
  FilterHistogram mLambdaHist;
  FilterHistogram mPhotonHist;
  std::vector<bool> mCollisionPassed;

  float mVtxZMin = 0.f;
  float mVtxZMax = 0.f;
  float mMultMin = 0.f;
  float mMultMax = 0.f;
  float mCentMin = 0.f;
  float mCentMax = 0.f;
  int8_t mMagFieldMin = 0;
  int8_t mMagFieldMax = 0;
  datatypes::CollisionMaskType mCollisionMask = 0;

  int mLambdaSign = 0;
  float mLambdaPtMin = 0.f;
  float mLambdaPtMax = 0.f;
  float mLambdaEtaMin = 0.f;
  float mLambdaEtaMax = 0.f;
  float mLambdaPhiMin = 0.f;
  float mLambdaPhiMax = 0.f;
  float mLambdaMassMin = 0.f;
  float mLambdaMassMax = 0.f;
  datatypes::V0MaskType mLambdaMask = 0;

  float mPhotonPtMin = 0.f;
  float mPhotonPtMax = 0.f;
  float mPhotonEtaMin = 0.f;
  float mPhotonEtaMax = 0.f;
  float mPhotonPhiMin = 0.f;
  float mPhotonPhiMax = 0.f;
  datatypes::PhotonMaskType mPhotonMask = 0;
};

// selection histograms of the output daughter masks
constexpr char LambdaSelHistName[] = "hLambdaSelection";
constexpr char PhotonSelHistName[] = "hPhotonSelection";
constexpr std::string_view LambdaSelsName = "Sigma0 Lambda daughter selection object";
constexpr std::string_view PhotonSelsName = "Sigma0 photon daughter selection object";

/// \brief Builds an output daughter mask from the daughter mask of the femto producer.
/// The selection values are bitmasks on the daughter mask: the first (loosest) one is the input mask used in the partition of the producer,
/// which is always fulfilled and gets no bit; each further bitmask gets one bit, set if the daughter mask contains all of its bits.
/// Downstream, one of the bits can be chosen as a tighter daughter selection (minimal cut).
/// Values are stored as double to represent 32-bit masks exactly.
template <auto& HistName, typename MaskType>
class DaughterMaskSelection : public baseselection::BaseSelection<double, MaskType, 1>
{
 public:
  DaughterMaskSelection() = default;
  ~DaughterMaskSelection() override = default;

  void configure(o2::framework::HistogramRegistry* registry, std::string const& selectionName, MaskType inputMask, std::vector<int> const& masks)
  {
    mMasks.clear();
    mMasks.push_back(inputMask);
    for (auto const& mask : masks) {
      if (mask < 0 || static_cast<uint64_t>(mask) > std::numeric_limits<MaskType>::max()) {
        LOG(fatal) << "Selection '" << selectionName << "': bitmask " << mask << " is not a valid daughter mask";
      }
      mMasks.push_back(static_cast<MaskType>(mask));
    }
    std::vector<double> values;
    std::vector<std::string> comments;
    for (auto const& mask : mMasks) {
      values.push_back(static_cast<double>(mask));
      comments.emplace_back(values.size() == 1 ? "input filter mask" : "");
    }
    this->addSelection(0, selectionName, values, limits::kEqualArray, true, true, false);
    this->addComments(0, comments);
    this->template setupSelectionHistogram<HistName>(registry);
  }

  MaskType applySelections(MaskType daughterMask)
  {
    this->reset();
    // kEqualArray: the i-th bit is set if the i-th value equals the i-th bitmask
    std::vector<double> values(mMasks.size(), -1.);
    for (std::size_t i = 0; i < mMasks.size(); i++) {
      if ((daughterMask & mMasks.at(i)) == mMasks.at(i)) {
        values.at(i) = static_cast<double>(mMasks.at(i));
      }
    }
    this->evaluateObservable(0, values);
    this->template assembleBitmask<HistName>();
    return this->getBitmask();
  }

 private:
  std::vector<MaskType> mMasks;
};

// Base selection configuration for downstream analysis tasks.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                                     \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "PDG code. Set sign to -1 for antiparticle"};           \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                               \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                             \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                                          \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                                           \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                            \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                               \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass"};                                  \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass"};                                  \
  o2::framework::Configurable<datatypes::V0MaskType> lambdaMask{"lambdaMask", 0, "Bitmask for the Lambda daughter (see Sigma0Bits)"}; \
  o2::framework::Configurable<datatypes::PhotonMaskType> photonMask{"photonMask", 0, "Bitmask for the photon daughter (see Sigma0Bits)"};

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
  o2::framework::Produces<o2::aod::FSigma0Labels> producedSigma0Labels;
};

struct ConfSigma0Tables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Sigma0Tables");
  o2::framework::Configurable<int> produceSigma0s{"produceSigma0s", -1, "Produce Sigma0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteSigma0s{"produceLiteSigma0s", -1, "Produce LiteSigma0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigma0Masks{"produceSigma0Masks", -1, "Produce Sigma0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigma0Extras{"produceSigma0Extras", -1, "Produce Sigma0Extras (-1: auto; 0 off; 1 on). Needs FLambdaExtras and FPhotonExtras as input"};
  o2::framework::Configurable<int> produceSigma0Labels{"produceSigma0Labels", -1, "Produce Sigma0Labels (-1: auto; 0 off; 1 on). Needs FLambdaLabels as input"};
};

/// \brief Builder for Sigma0 <- Lambda + PCM gamma candidates
/// Pairs daughters within the collision and applies loose kinematic cuts
class Sigma0Builder
{
 public:
  Sigma0Builder() = default;
  ~Sigma0Builder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void init(o2::framework::HistogramRegistry* registry, T1 const& confLambdaSelection, T2 const& confPhotonSelection, T3& confFilter, T4& confBits, T5& confQa, T6& confTable, T7& initContext)
  {
    LOG(info) << "Initialize femto Sigma0 builder...";
    mProduceSigma0s = utils::enableTable("FSigma0s_001", confTable.produceSigma0s.value, initContext);
    mProduceLiteSigma0s = utils::enableTable("FLiteSigma0s_001", confTable.produceLiteSigma0s.value, initContext);
    mProduceSigma0Masks = utils::enableTable("FSigma0Masks_001", confTable.produceSigma0Masks.value, initContext);
    mProduceSigma0Extras = utils::enableTable("FSigma0Extras_001", confTable.produceSigma0Extras.value, initContext);
    mProduceSigma0Labels = utils::enableTable("FSigma0Labels_001", confTable.produceSigma0Labels.value, initContext);

    if (mProduceSigma0s && mProduceLiteSigma0s) {
      LOG(fatal) << "FSigma0s and FLiteSigma0s are mutually exclusive -- enable only one. "
                 << "FLiteSigma0s is meant to replace FSigma0s at the producer stage for better compression in derived data.";
    }

    mFillAnyTable = mProduceSigma0s || mProduceLiteSigma0s || mProduceSigma0Masks || mProduceSigma0Extras || mProduceSigma0Labels;
    if (!mFillAnyTable) {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mPtMin = confFilter.ptMin.value;
    mPtMax = confFilter.ptMax.value;
    mEtaMin = confFilter.etaMin.value;
    mEtaMax = confFilter.etaMax.value;
    mPhiMin = confFilter.phiMin.value;
    mPhiMax = confFilter.phiMax.value;
    mMassMin = confFilter.massMin.value;
    mMassMax = confFilter.massMax.value;
    mSigma0FilterHist.init(registry, Sigma0FilterHistName,
                           {
                             {"ptMin", FilterHistogram::toString(mPtMin)},
                             {"ptMax", FilterHistogram::toString(mPtMax)},
                             {"etaMin", FilterHistogram::toString(mEtaMin)},
                             {"etaMax", FilterHistogram::toString(mEtaMax)},
                             {"phiMin", FilterHistogram::toString(mPhiMin)},
                             {"phiMax", FilterHistogram::toString(mPhiMax)},
                             {"massMin", FilterHistogram::toString(mMassMin)},
                             {"massMax", FilterHistogram::toString(mMassMax)},
                           });

    mLambdaSelection.configure(registry, "Lambda daughter mask", confLambdaSelection.mask.value, confBits.lambdaMask.value);
    mLambdaSelection.printSelections(LambdaSelsName.data());
    mPhotonSelection.configure(registry, "Photon daughter mask", confPhotonSelection.mask.value, confBits.photonMask.value);
    mPhotonSelection.printSelections(PhotonSelsName.data());

    mFillQa = confQa.fillQa.value;
    if (mFillQa) {
      mRegistry = registry;
      mRegistry->add(QaSigma0Mass, "#Sigma^{0} invariant mass;#it{M}_{#Lambda#gamma} (GeV/#it{c}^{2});Entries", o2::framework::HistType::kTH1F, {confQa.sigma0Mass});
      mRegistry->add(QaLambdaMass, "#Lambda daughter invariant mass;#it{M}_{#Lambda} (GeV/#it{c}^{2});Entries", o2::framework::HistType::kTH1F, {confQa.lambdaMass});
      mRegistry->add(QaSigma0Pt, "#Sigma^{0} transverse momentum;#it{p}_{T} (GeV/#it{c});Entries", o2::framework::HistType::kTH1F, {confQa.pt});
      mRegistry->add(QaLambdaPt, "#Lambda daughter transverse momentum;#it{p}_{T} (GeV/#it{c});Entries", o2::framework::HistType::kTH1F, {confQa.pt});
      mRegistry->add(QaPhotonPt, "#gamma daughter transverse momentum;#it{p}_{T} (GeV/#it{c});Entries", o2::framework::HistType::kTH1F, {confQa.pt});
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

  bool checkFilters()
  {
    bool pass = true;
    pass &= mSigma0FilterHist.fill(kSigma0PtMin, mPt > mPtMin);
    pass &= mSigma0FilterHist.fill(kSigma0PtMax, mPt < mPtMax);
    pass &= mSigma0FilterHist.fill(kSigma0EtaMin, mEta > mEtaMin);
    pass &= mSigma0FilterHist.fill(kSigma0EtaMax, mEta < mEtaMax);
    pass &= mSigma0FilterHist.fill(kSigma0PhiMin, mPhi >= mPhiMin);
    pass &= mSigma0FilterHist.fill(kSigma0PhiMax, mPhi < mPhiMax);
    pass &= mSigma0FilterHist.fill(kSigma0MassMin, mMass > mMassMin);
    pass &= mSigma0FilterHist.fill(kSigma0MassMax, mMass < mMassMax);
    mSigma0FilterHist.fillSummary(pass);
    return pass;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillSigma0(T1 const& col, T2 const& lambda, T3 const& photon, T4& products)
  {
    reconstructSigma0(lambda, photon);
    if (!checkFilters()) {
      return;
    }
    if (mFillQa) {
      mRegistry->fill(HIST(QaSigma0Mass), mMass);
      mRegistry->fill(HIST(QaLambdaMass), lambda.mass());
      mRegistry->fill(HIST(QaSigma0Pt), mPt);
      mRegistry->fill(HIST(QaLambdaPt), lambda.pt());
      mRegistry->fill(HIST(QaPhotonPt), photon.pt());
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
      products.producedSigma0Masks(mLambdaSelection.applySelections(lambda.mask()), mPhotonSelection.applySelections(photon.mask()));
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
    if (mProduceSigma0Labels) {
      // PCM photons carry no MC labels and the Sigma0 itself is not stored as a femto mc particle, so the
      // Sigma0 is labelled with the mc particle of its Lambda daughter. A true Sigma0 is a Lambda whose
      // mc mother is a Sigma0 (see Sigma0HistManager); the photon is not truth-matched.
      if constexpr (requires { lambda.fMcParticleId(); }) {
        products.producedSigma0Labels(lambda.fMcParticleId());
      } else {
        LOG(fatal) << "FSigma0Labels requested, but FLambdaLabels are not joined to the input Lambda table";
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

  // filters on the Lambda-photon pair
  FilterHistogram mSigma0FilterHist;
  float mPtMin = 0.f;
  float mPtMax = 0.f;
  float mEtaMin = 0.f;
  float mEtaMax = 0.f;
  float mPhiMin = 0.f;
  float mPhiMax = 0.f;
  float mMassMin = 0.f;
  float mMassMax = 0.f;

  // output daughter masks
  DaughterMaskSelection<LambdaSelHistName, datatypes::V0MaskType> mLambdaSelection;
  DaughterMaskSelection<PhotonSelHistName, datatypes::PhotonMaskType> mPhotonSelection;

  // minimal QA
  o2::framework::HistogramRegistry* mRegistry = nullptr;
  bool mFillQa = false;

  bool mFillAnyTable = false;
  bool mProduceSigma0s = false;
  bool mProduceLiteSigma0s = false;
  bool mProduceSigma0Masks = false;
  bool mProduceSigma0Extras = false;
  bool mProduceSigma0Labels = false;
};

} // namespace o2::analysis::femto::sigma0builder
#endif // PWGCF_FEMTO_CORE_SIGMA0BUILDER_H_
