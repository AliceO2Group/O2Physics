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

/// \file xi1530kaoncorrelation.cxx
/// \brief Xi(1530)^0--K femtoscopic correlation task on compact reduced tables.
/// \author Prottay Das, prottay.das@cern.ch
///
/// Some information: Same-event and mixed-event processing are intentionally implemented as
/// separate PROCESS_SWITCH functions.  The task reconstructs all downstream
/// selections from the producer bitmaps/raw kaon PID, keeps M(Xi pi) as an
/// explicit axis for signal extraction in each k* bin, and uses the
/// reconstructed Xi-pi mass in the k* calculation.
///
/// Default PID choices:
///   * prompt Xi* pion: TPC-only, 3 sigma;
///   * Xi bachelor pion: TPC-only, 3 sigma;
///   * Lambda pion/proton: TPC-only, 4 sigma;
///   * external K: TPC-only below pT=0.5 GeV/c and circular TPC+TOF above,
///     requiring TOF above threshold.
///
/// The producer also stores alternative prompt/bachelor PID decisions and raw
/// kaon TPC/TOF PID, so those choices I kept configurable downstream (will change later).

#include "PWGLF/DataModel/LFReducedXi1530KaonTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::soa;

namespace
{
constexpr float kKaonMass = o2::constants::physics::MassKPlus;
constexpr float kPi = 3.14159265358979323846f;
constexpr float kTwoPi = 2.f * kPi;

inline bool hasTrackBit(uint64_t word, redxistark::XiStarTrackSelBit bit)
{
  return (word & redxistark::xiStarTrackSelMask(bit)) != 0u;
}

inline bool hasTopoBit(uint64_t word, redxistark::XiStarTopoSelBit bit)
{
  return (word & redxistark::xiStarTopoSelMask(bit)) != 0u;
}

inline bool hasKaonBit(uint8_t word, redxistarkaon::KaonSelBit bit)
{
  return (word & redxistarkaon::kaonSelMask(bit)) != 0u;
}

inline float wrapToPi(float x)
{
  while (x > kPi) {
    x -= kTwoPi;
  }
  while (x <= -kPi) {
    x += kTwoPi;
  }
  return x;
}
} // namespace

struct xi1530kaoncorrelation {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ---------------------------------------------------------------------------
  // Xi(1530) selection reconstructed from the producer bitmaps.
  //
  // Default daughter PID is TPC-only; producer-stored alternatives remain selectable.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "xiStarSelection";

    Configurable<int> promptPionPIDMode{"promptPionPIDMode", 2, "Prompt pion PID mode: 0=availability hybrid, 1=pT-threshold hybrid, 2=TPC-only"};
    Configurable<int> promptPionPIDWP{"promptPionPIDWP", 0, "Prompt pion TPC PID: 0=3sigma, 1=4sigma, 2=5sigma, 3=6sigma"};
    Configurable<int> bachelorPionPIDMode{"bachelorPionPIDMode", 1, "Xi bachelor pion PID mode: 0=availability hybrid, 1=TPC-only"};
    Configurable<int> bachelorPionPIDWP{"bachelorPionPIDWP", 0, "Xi bachelor pion PID: 0=3sigma, 1=4sigma, 2=4.8sigma, 3=5sigma, 4=6sigma; 4.8 is TPC-only"};
    Configurable<int> lambdaPionPIDWP{"lambdaPionPIDWP", 1, "Lambda pion TPC PID: 0=3sigma, 1=4sigma, 2=4.8sigma, 3=5sigma, 4=6sigma"};
    Configurable<int> lambdaProtonPIDWP{"lambdaProtonPIDWP", 1, "Lambda proton TPC PID: 0=3sigma, 1=4sigma, 2=5sigma, 3=6sigma"};

    Configurable<int> promptPionRowsWP{"promptPionRowsWP", 1, "Prompt-pion rows: 0=>70, 1=>80, 2=>90"};
    Configurable<int> promptPionDCAzWP{"promptPionDCAzWP", 1, "Prompt-pion DCAz: 0=loose 1 cm bit, 1=producer default bit, 2=tight 0.1 cm bit"};
    Configurable<int> v0DaughterRowsWP{"v0DaughterRowsWP", 1, "Cascade/V0 daughter rows: 0=all cascade daughters broad bit, 1=V0 rows nominal bit, 2=Var1, 3=Var2"};

    Configurable<int> v0PionDcaPVWP{"v0PionDcaPVWP", 3, "V0 pion DCA-to-PV bit: 0=0.05, 1=0.06, 2=0.10, 3=0.20 cm"};
    Configurable<int> v0ProtonDcaPVWP{"v0ProtonDcaPVWP", 2, "V0 proton DCA-to-PV bit: 0=0.05, 1=0.06, 2=0.07, 3=0.10 cm"};
    Configurable<int> v0DcaPVWP{"v0DcaPVWP", 1, "Minimum Lambda/V0 DCA-to-PV bit: 0=0, 1=0.03, 2=0.10 cm"};
    Configurable<int> v0DcaDaughtersWP{"v0DcaDaughtersWP", 0, "V0 daughter DCA/native metric bit: 0=<1.0, 1=<0.5, 2=<0.1"};
    Configurable<int> v0CosPAWP{"v0CosPAWP", 4, "V0 cosPA bit: 0=0.97, 1=0.98, 2=0.9876, 3=0.99, 4=0.995"};
    Configurable<int> v0RadiusWP{"v0RadiusWP", 4, "V0 minimum-radius bit: 0=0.9, 1=1.01, 2=1.2, 3=2.5, 4=3.0 cm"};
    Configurable<int> lambdaMassWP{"lambdaMassWP", 1, "Lambda mass-window bit: 0=10, 1=8, 2=6, 3=11.6 MeV/c2"};
    Configurable<bool> applyV0LifetimeCut{"applyV0LifetimeCut", true, "Apply producer Lambda-lifetime bit"};
    Configurable<int> v0LifetimeWP{"v0LifetimeWP", 0, "Lambda lifetime bit: 0=default, 1=Var1, 2=Var2; thresholds are those used in producer JSON"};

    Configurable<int> cascBachelorDcaPVWP{"cascBachelorDcaPVWP", 1, "Xi bachelor DCA-to-PV bit: 0=0.05, 1=0.06, 2=0.10 cm"};
    Configurable<int> cascDcaDaughtersWP{"cascDcaDaughtersWP", 1, "Cascade daughter DCA/native metric bit: 0=<1.0, 1=<0.25, 2=<0.20"};
    Configurable<int> cascCosPAWP{"cascCosPAWP", 1, "Cascade cosPA bit: 0=0.97, 1=0.98, 2=0.9947, 3=0.995"};
    Configurable<int> cascRadiusWP{"cascRadiusWP", 2, "Cascade minimum-radius bit: 0=0.9, 1=1.0, 2=1.01, 3=1.3 cm"};
    Configurable<int> xiMassWP{"xiMassWP", 1, "Xi mass-window bit: 0=10, 1=8, 2=6 MeV/c2"};
    Configurable<bool> applyXiLifetimeCut{"applyXiLifetimeCut", true, "Apply producer Xi-lifetime bit"};
    Configurable<int> xiLifetimeWP{"xiLifetimeWP", 0, "Xi lifetime bit: 0=default, 1=Var1, 2=Var2; thresholds are those used in producer JSON"};

    Configurable<bool> requirePPBachBaryonDCAxy{"requirePPBachBaryonDCAxy", false, "Require pp-Xi bachelor-baryon DCAxy bit"};
    Configurable<bool> requirePPXiRapidity{"requirePPXiRapidity", false, "Require pp-Xi rapidity bit"};
  } xiCuts;

  // ---------------------------------------------------------------------------
  // External K selection.
  // Raw TPC/TOF n-sigma values are stored by the producer.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "kaonSelection";

    Configurable<float> ptMin{"ptMin", 0.15f, "Minimum external-K pT (GeV/c)"};
    Configurable<float> ptMax{"ptMax", 100.f, "Maximum external-K pT (GeV/c)"};
    Configurable<int> tpcRowsMin{"tpcRowsMin", 80, "Minimum TPC crossed rows, strict >"};

    Configurable<int> dcaXYWP{"dcaXYWP", 0, "K DCAxy bit: 0=default, 1=Var1, 2=Var2"};
    Configurable<int> dcaZWP{"dcaZWP", 0, "K DCAz bit: 0=default, 1=Var1, 2=Var2"};

    Configurable<int> pidMode{"pidMode", 2, "K PID mode: 0=TPC-only, 1=availability hybrid, 2=pT-threshold hybrid"};
    Configurable<float> pidPtThreshold{"pidPtThreshold", 0.5f, "pT threshold between TPC-only and circular TPC+TOF PID (GeV/c)"};
    Configurable<float> tpcNSigmaMax{"tpcNSigmaMax", 3.f, "Maximum |TPC nSigma_K| below threshold"};
    Configurable<float> circularNSigmaMax{"circularNSigmaMax", 3.f, "Maximum sqrt(TPC^2+TOF^2) above threshold"};
    Configurable<bool> requireTOFAboveThreshold{"requireTOFAboveThreshold", true, "Require TOF for K at/above PID threshold"};
  } kaonCuts;

  // ---------------------------------------------------------------------------
  // Pair/channel selections.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "pairSelection";

    Configurable<bool> fillMatter{"fillMatter", true, "Process Xi*-K- matter channels"};
    Configurable<bool> fillAntimatter{"fillAntimatter", true, "Process anti-Xi*-K+ antimatter channels"};
    Configurable<bool> fillLikeSignControl{"fillLikeSignControl", true, "Keep like-sign Xi-pi control channels (+/-2)"};
    Configurable<bool> applySharedTrackCleaning{"applySharedTrackCleaning", true, "Reject SE kaons sharing any Xi* daughter track"};
  } pairCuts;

  // ---------------------------------------------------------------------------
  // Event mixing.
  // ---------------------------------------------------------------------------
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of previous events to mix"};
  ConfigurableAxis cfgMixVtxBins{"cfgMixVtxBins", {VARIABLE_WIDTH, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0}, "Mixing bins in z vertex (cm)"};
  ConfigurableAxis cfgMixFT0MBins{"cfgMixFT0MBins", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 101.0, 110.0}, "Mixing bins in FT0M percentile"};
  ConfigurableAxis cfgMixBzBins{"cfgMixBzBins", {VARIABLE_WIDTH, -1.0, -0.1, 0.1, 1.0}, "Mixing bins in Bz (T)"};

  using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ,
                                            aod::redxistarkevent::FT0MPercentile,
                                            aod::redxistarkevent::Bz>;
  MixingBinning mixingBinning{{cfgMixVtxBins, cfgMixFT0MBins, cfgMixBzBins}, true};

  Preslice<aod::XiStarCandidates> xiStarsPerEvent = aod::redxistarcandidate::redXiStarKEventId;
  Preslice<aod::KaonCandidates> kaonsPerEvent = aod::redxistarkaon::redXiStarKEventId;

  // ---------------------------------------------------------------------------
  // Histogram axes.
  // ---------------------------------------------------------------------------
  ConfigurableAxis cfgAxisKstar{"cfgAxisKstar", {200, 0.0, 1.0}, "k* axis (GeV/c)"};
  ConfigurableAxis cfgAxisMass{"cfgAxisMass", {500, 1.45, 1.70}, "M(Xi pi) axis (GeV/c2)"};
  ConfigurableAxis cfgAxisMt{"cfgAxisMt", {100, 0.0, 5.0}, "pair mT axis (GeV/c2)"};
  ConfigurableAxis cfgAxisXiPt{"cfgAxisXiPt", {100, 0.0, 10.0}, "Xi(1530) pT axis (GeV/c)"};
  ConfigurableAxis cfgAxisKaonPt{"cfgAxisKaonPt", {100, 0.0, 5.0}, "external K pT axis (GeV/c)"};
  ConfigurableAxis cfgAxisPIDP{"cfgAxisPIDP", {200, 0.0, 10.0}, "momentum axis for PID QA (GeV/c)"};
  ConfigurableAxis cfgAxisPIDNSigma{"cfgAxisPIDNSigma", {160, -8.0, 8.0}, "n-sigma axis for PID QA"};

  Configurable<bool> fillCandidateQA{"fillCandidateQA", true, "Fill selected Xi* and kaon QA"};
  Configurable<bool> fillKinematicQA{"fillKinematicQA", true, "Fill k*-mT and k*-pT QA histograms"};
  Configurable<bool> fillCPRQA{"fillCPRQA", true, "Fill K--prompt-pion close-pair QA; no CPR cut is applied"};
  Configurable<float> cprKstarMax{"cprKstarMax", 0.5f, "Maximum k* filled in CPR QA (GeV/c)"};
  ConfigurableAxis cfgAxisCPRDEta{"cfgAxisCPRDEta", {120, -0.3, 0.3}, "Delta-eta axis for CPR QA"};
  ConfigurableAxis cfgAxisCPRDPhiStar{"cfgAxisCPRDPhiStar", {120, -0.3, 0.3}, "Delta-phi* axis for CPR QA"};
  ConfigurableAxis cfgAxisCPRKstar{"cfgAxisCPRKstar", {100, 0.0, 0.5}, "k* axis for CPR QA (GeV/c)"};
  Configurable<float> cprRadiusMinCm{"cprRadiusMinCm", 85.f, "Minimum TPC radius used for phi* averaging (cm)"};
  Configurable<float> cprRadiusMaxCm{"cprRadiusMaxCm", 245.f, "Maximum TPC radius used for phi* averaging (cm)"};
  Configurable<float> cprRadiusStepCm{"cprRadiusStepCm", 20.f, "TPC-radius step used for phi* averaging (cm)"};

  // ---------------------------------------------------------------------------
  // Xi(1530) selection helpers.
  // ---------------------------------------------------------------------------
  bool passPromptPionPID(uint64_t bits) const
  {
    const int mode = xiCuts.promptPionPIDMode.value;
    const int wp = xiCuts.promptPionPIDWP.value;

    if (mode == 0) {
      switch (wp) {
        case 0:
          return hasTrackBit(bits, redxistark::kPiFirstPID3Sigma);
        case 1:
          return hasTrackBit(bits, redxistark::kPiFirstPID4Sigma);
        case 2:
          return hasTrackBit(bits, redxistark::kPiFirstPID5Sigma);
        case 3:
          return hasTrackBit(bits, redxistark::kPiFirstPID6Sigma);
        default:
          return false;
      }
    }

    if (mode == 1) {
      switch (wp) {
        case 0:
          return hasTrackBit(bits, redxistark::kPiFirstPtThresholdPID3Sigma);
        case 1:
          return hasTrackBit(bits, redxistark::kPiFirstPtThresholdPID4Sigma);
        case 2:
          return hasTrackBit(bits, redxistark::kPiFirstPtThresholdPID5Sigma);
        case 3:
          return hasTrackBit(bits, redxistark::kPiFirstPtThresholdPID6Sigma);
        default:
          return false;
      }
    }

    switch (wp) {
      case 0:
        return hasTrackBit(bits, redxistark::kPiFirstTPC3Sigma);
      case 1:
        return hasTrackBit(bits, redxistark::kPiFirstTPC4Sigma);
      case 2:
        return hasTrackBit(bits, redxistark::kPiFirstTPC5Sigma);
      case 3:
        return hasTrackBit(bits, redxistark::kPiFirstTPC6Sigma);
      default:
        return false;
    }
  }

  bool passBachelorPionPID(uint64_t bits) const
  {
    const int mode = xiCuts.bachelorPionPIDMode.value;
    const int wp = xiCuts.bachelorPionPIDWP.value;

    if (mode == 0) {
      switch (wp) {
        case 0:
          return hasTrackBit(bits, redxistark::kXiBachelorPiHybridPID3Sigma);
        case 1:
          return hasTrackBit(bits, redxistark::kXiBachelorPiHybridPID4Sigma);
        case 3:
          return hasTrackBit(bits, redxistark::kXiBachelorPiHybridPID5Sigma);
        case 4:
          return hasTrackBit(bits, redxistark::kXiBachelorPiHybridPID6Sigma);
        default:
          return false; // 4.8 sigma exists only for TPC-only.
      }
    }

    switch (wp) {
      case 0:
        return hasTrackBit(bits, redxistark::kXiBachelorPiTPC3Sigma);
      case 1:
        return hasTrackBit(bits, redxistark::kXiBachelorPiTPC4Sigma);
      case 2:
        return hasTrackBit(bits, redxistark::kXiBachelorPiTPC4p8Sigma);
      case 3:
        return hasTrackBit(bits, redxistark::kXiBachelorPiTPC5Sigma);
      case 4:
        return hasTrackBit(bits, redxistark::kXiBachelorPiTPC6Sigma);
      default:
        return false;
    }
  }

  bool passLambdaPionTPCPID(uint64_t bits) const
  {
    switch (xiCuts.lambdaPionPIDWP.value) {
      case 0:
        return hasTrackBit(bits, redxistark::kLambdaPiTPC3Sigma);
      case 1:
        return hasTrackBit(bits, redxistark::kLambdaPiTPC4Sigma);
      case 2:
        return hasTrackBit(bits, redxistark::kLambdaPiTPC4p8Sigma);
      case 3:
        return hasTrackBit(bits, redxistark::kLambdaPiTPC5Sigma);
      case 4:
        return hasTrackBit(bits, redxistark::kLambdaPiTPC6Sigma);
      default:
        return false;
    }
  }

  bool passLambdaProtonTPCPID(uint64_t bits) const
  {
    switch (xiCuts.lambdaProtonPIDWP.value) {
      case 0:
        return hasTrackBit(bits, redxistark::kLambdaPrTPC3Sigma);
      case 1:
        return hasTrackBit(bits, redxistark::kLambdaPrTPC4Sigma);
      case 2:
        return hasTrackBit(bits, redxistark::kLambdaPrTPC5Sigma);
      case 3:
        return hasTrackBit(bits, redxistark::kLambdaPrTPC6Sigma);
      default:
        return false;
    }
  }

  template <typename TXiStar>
  bool passXiStarSelection(TXiStar const& candidate) const
  {
    const uint64_t trackBits = candidate.trackSelectionBits();
    const uint64_t topoBits = candidate.topologySelectionBits();

    if (!passPromptPionPID(trackBits) ||
        !passBachelorPionPID(trackBits) ||
        !passLambdaPionTPCPID(trackBits) ||
        !passLambdaProtonTPCPID(trackBits)) {
      return false;
    }

    const std::array<redxistark::XiStarTrackSelBit, 3> promptRowsBits{
      redxistark::kPiFirstTPCRows70,
      redxistark::kPiFirstTPCRows80,
      redxistark::kPiFirstTPCRows90};
    if (!hasTrackBit(trackBits, promptRowsBits[xiCuts.promptPionRowsWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTrackSelBit, 3> promptDCAzBits{
      redxistark::kPiFirstDCAz1cm,
      redxistark::kPiFirstDCAzDefault,
      redxistark::kPiFirstDCAz0p1cm};
    if (!hasTrackBit(trackBits, promptDCAzBits[xiCuts.promptPionDCAzWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTrackSelBit, 4> v0RowsBits{
      redxistark::kAllCascDaughtersTPCRows50,
      redxistark::kV0DaughtersTPCRows70,
      redxistark::kV0DaughtersTPCRowsVar1,
      redxistark::kV0DaughtersTPCRowsVar2};
    if (!hasTrackBit(trackBits, v0RowsBits[xiCuts.v0DaughterRowsWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 4> pionDcaBits{
      redxistark::kV0PionDcaPV005,
      redxistark::kV0PionDcaPV006,
      redxistark::kV0PionDcaPV010,
      redxistark::kV0PionDcaPV020};
    if (!hasTopoBit(topoBits, pionDcaBits[xiCuts.v0PionDcaPVWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 4> protonDcaBits{
      redxistark::kV0ProtonDcaPV005,
      redxistark::kV0ProtonDcaPV006,
      redxistark::kV0ProtonDcaPV007,
      redxistark::kV0ProtonDcaPV010};
    if (!hasTopoBit(topoBits, protonDcaBits[xiCuts.v0ProtonDcaPVWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 3> v0DcaPVBits{
      redxistark::kV0DcaPV000,
      redxistark::kV0DcaPV003,
      redxistark::kV0DcaPV010};
    if (!hasTopoBit(topoBits, v0DcaPVBits[xiCuts.v0DcaPVWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 3> v0DaughterDCABits{
      redxistark::kV0DcaDaughters1p0,
      redxistark::kV0DcaDaughters0p5,
      redxistark::kV0DcaDaughters0p1};
    if (!hasTopoBit(topoBits, v0DaughterDCABits[xiCuts.v0DcaDaughtersWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 5> v0CosPABits{
      redxistark::kV0CosPA097,
      redxistark::kV0CosPA098,
      redxistark::kV0CosPA09876,
      redxistark::kV0CosPA099,
      redxistark::kV0CosPA0995};
    if (!hasTopoBit(topoBits, v0CosPABits[xiCuts.v0CosPAWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 5> v0RadiusBits{
      redxistark::kV0Radius0p9,
      redxistark::kV0Radius1p01,
      redxistark::kV0Radius1p2,
      redxistark::kV0Radius2p5,
      redxistark::kV0Radius3p0};
    if (!hasTopoBit(topoBits, v0RadiusBits[xiCuts.v0RadiusWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 4> lambdaMassBits{
      redxistark::kLambdaMass10MeV,
      redxistark::kLambdaMass8MeV,
      redxistark::kLambdaMass6MeV,
      redxistark::kLambdaMass11p6MeV};
    if (!hasTopoBit(topoBits, lambdaMassBits[xiCuts.lambdaMassWP.value])) {
      return false;
    }

    if (xiCuts.applyV0LifetimeCut.value) {
      const std::array<redxistark::XiStarTopoSelBit, 3> v0LifetimeBits{
        redxistark::kV0LifetimeDefault,
        redxistark::kV0LifetimeVar1,
        redxistark::kV0LifetimeVar2};
      if (!hasTopoBit(topoBits, v0LifetimeBits[xiCuts.v0LifetimeWP.value])) {
        return false;
      }
    }

    const std::array<redxistark::XiStarTopoSelBit, 3> bachelorDcaBits{
      redxistark::kCascBachelorDcaPV005,
      redxistark::kCascBachelorDcaPV006,
      redxistark::kCascBachelorDcaPV010};
    if (!hasTopoBit(topoBits, bachelorDcaBits[xiCuts.cascBachelorDcaPVWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 3> cascDaughterDCABits{
      redxistark::kCascDcaDaughters1p0,
      redxistark::kCascDcaDaughters0p25,
      redxistark::kCascDcaDaughters0p20};
    if (!hasTopoBit(topoBits, cascDaughterDCABits[xiCuts.cascDcaDaughtersWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 4> cascCosPABits{
      redxistark::kCascCosPA097,
      redxistark::kCascCosPA098,
      redxistark::kCascCosPA09947,
      redxistark::kCascCosPA0995};
    if (!hasTopoBit(topoBits, cascCosPABits[xiCuts.cascCosPAWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 4> cascRadiusBits{
      redxistark::kCascRadius0p9,
      redxistark::kCascRadius1p0,
      redxistark::kCascRadius1p01,
      redxistark::kCascRadius1p3};
    if (!hasTopoBit(topoBits, cascRadiusBits[xiCuts.cascRadiusWP.value])) {
      return false;
    }

    const std::array<redxistark::XiStarTopoSelBit, 3> xiMassBits{
      redxistark::kXiMass10MeV,
      redxistark::kXiMass8MeV,
      redxistark::kXiMass6MeV};
    if (!hasTopoBit(topoBits, xiMassBits[xiCuts.xiMassWP.value])) {
      return false;
    }

    if (xiCuts.applyXiLifetimeCut.value) {
      const std::array<redxistark::XiStarTopoSelBit, 3> xiLifetimeBits{
        redxistark::kXiLifetimeDefault,
        redxistark::kXiLifetimeVar1,
        redxistark::kXiLifetimeVar2};
      if (!hasTopoBit(topoBits, xiLifetimeBits[xiCuts.xiLifetimeWP.value])) {
        return false;
      }
    }

    if (xiCuts.requirePPBachBaryonDCAxy.value &&
        !hasTopoBit(topoBits, redxistark::kBachBaryonDCAxy0020)) {
      return false;
    }
    if (xiCuts.requirePPXiRapidity.value &&
        !hasTopoBit(topoBits, redxistark::kXiRapidity05)) {
      return false;
    }

    return true;
  }

  // ---------------------------------------------------------------------------
  // Kaon selection helpers.
  // IMPORTANT: getter spelling follows the attached data model exactly:
  //   tPCNClsCrossedRows(), tPCNSigmaKa(), tOFNSigmaKa().
  // ---------------------------------------------------------------------------
  bool passKaonDCA(uint8_t bits, int wp, bool xy) const
  {
    if (xy) {
      switch (wp) {
        case 0:
          return hasKaonBit(bits, redxistarkaon::kKaonDCAxyDefault);
        case 1:
          return hasKaonBit(bits, redxistarkaon::kKaonDCAxyVar1);
        case 2:
          return hasKaonBit(bits, redxistarkaon::kKaonDCAxyVar2);
        default:
          return false;
      }
    }

    switch (wp) {
      case 0:
        return hasKaonBit(bits, redxistarkaon::kKaonDCAzDefault);
      case 1:
        return hasKaonBit(bits, redxistarkaon::kKaonDCAzVar1);
      case 2:
        return hasKaonBit(bits, redxistarkaon::kKaonDCAzVar2);
      default:
        return false;
    }
  }

  template <typename TKaon>
  bool passKaonPID(TKaon const& kaon) const
  {
    const float tpc = kaon.tPCNSigmaKa();
    if (!std::isfinite(tpc)) {
      return false;
    }

    if (kaonCuts.pidMode.value == 0) {
      return std::abs(tpc) < kaonCuts.tpcNSigmaMax.value;
    }

    if (kaonCuts.pidMode.value == 1) {
      if (!kaon.hasTOF()) {
        return std::abs(tpc) < kaonCuts.tpcNSigmaMax.value;
      }
      const float tof = kaon.tOFNSigmaKa();
      return std::isfinite(tof) && std::hypot(tpc, tof) < kaonCuts.circularNSigmaMax.value;
    }

    // Default: pT-threshold hybrid PID.
    if (kaon.pt() < kaonCuts.pidPtThreshold.value) {
      return std::abs(tpc) < kaonCuts.tpcNSigmaMax.value;
    }

    if (!kaon.hasTOF()) {
      if (kaonCuts.requireTOFAboveThreshold.value) {
        return false;
      }
      return std::abs(tpc) < kaonCuts.tpcNSigmaMax.value;
    }

    const float tof = kaon.tOFNSigmaKa();
    return std::isfinite(tof) && std::hypot(tpc, tof) < kaonCuts.circularNSigmaMax.value;
  }

  template <typename TKaon>
  bool passKaonSelection(TKaon const& kaon) const
  {
    if (kaon.pt() <= kaonCuts.ptMin.value || kaon.pt() >= kaonCuts.ptMax.value) {
      return false;
    }
    if (kaon.tPCNClsCrossedRows() <= kaonCuts.tpcRowsMin.value) {
      return false;
    }
    if (!passKaonDCA(kaon.selectionBits(), kaonCuts.dcaXYWP.value, true) ||
        !passKaonDCA(kaon.selectionBits(), kaonCuts.dcaZWP.value, false)) {
      return false;
    }
    return passKaonPID(kaon);
  }

  // ---------------------------------------------------------------------------
  // Pair helpers.
  // ---------------------------------------------------------------------------
  template <typename TXiStar, typename TKaon>
  bool sharesTrack(TXiStar const& xiStar, TKaon const& kaon) const
  {
    const int64_t kIndex = kaon.trackIndex();
    return kIndex == xiStar.xiStarPionIndex() ||
           kIndex == xiStar.xiBachelorIndex() ||
           kIndex == xiStar.v0PositiveIndex() ||
           kIndex == xiStar.v0NegativeIndex();
  }

  bool passChannelSelection(int8_t channel) const
  {
    if (channel == 0) {
      return false;
    }
    if (channel > 0 && !pairCuts.fillMatter.value) {
      return false;
    }
    if (channel < 0 && !pairCuts.fillAntimatter.value) {
      return false;
    }
    if (std::abs(static_cast<int>(channel)) == 2 && !pairCuts.fillLikeSignControl.value) {
      return false;
    }
    return true;
  }

  bool isRequestedChargeCombination(int8_t channel, int8_t kaonCharge) const
  {
    // +1/+2 are Xi- based -> K- ; -1/-2 are anti-Xi+ based -> K+.
    if (channel > 0) {
      return kaonCharge < 0;
    }
    if (channel < 0) {
      return kaonCharge > 0;
    }
    return false;
  }

  int promptPionCharge(int8_t channel) const
  {
    // +1 = Xi- pi+ US, +2 = Xi- pi- LS,
    // -1 = antiXi+ pi- US, -2 = antiXi+ pi+ LS.
    if (channel == redxistark::kXiStarUS || channel == redxistark::kAntiXiStarLS) {
      return +1;
    }
    return -1;
  }

  template <typename TXiStar, typename TKaon>
  float getKstar(TXiStar const& xiStar, TKaon const& kaon) const
  {
    const auto m1 = static_cast<double>(xiStar.mass());
    const auto m2 = static_cast<double>(kKaonMass);

    const auto p1x = static_cast<double>(xiStar.px());
    const auto p1y = static_cast<double>(xiStar.py());
    const auto p1z = static_cast<double>(xiStar.pz());
    const auto p2x = static_cast<double>(kaon.px());
    const auto p2y = static_cast<double>(kaon.py());
    const auto p2z = static_cast<double>(kaon.pz());

    const double e1 = std::sqrt(p1x * p1x + p1y * p1y + p1z * p1z + m1 * m1);
    const double e2 = std::sqrt(p2x * p2x + p2y * p2y + p2z * p2z + m2 * m2);

    const double px = p1x + p2x;
    const double py = p1y + p2y;
    const double pz = p1z + p2z;
    const double e = e1 + e2;
    const double s = e * e - px * px - py * py - pz * pz;

    if (!std::isfinite(s) || s <= 0.0) {
      return -1.f;
    }

    const double mPlus = m1 + m2;
    const double mMinus = m1 - m2;
    const double lambda = (s - mPlus * mPlus) * (s - mMinus * mMinus);
    if (!std::isfinite(lambda) || lambda <= 0.0) {
      return -1.f;
    }

    const double kstar = std::sqrt(lambda) / (2.0 * std::sqrt(s));
    return std::isfinite(kstar) ? static_cast<float>(kstar) : -1.f;
  }

  template <typename TXiStar, typename TKaon>
  float getPairMt(TXiStar const& xiStar, TKaon const& kaon) const
  {
    const float pairPx = xiStar.px() + kaon.px();
    const float pairPy = xiStar.py() + kaon.py();
    const float kT = 0.5f * std::hypot(pairPx, pairPy);
    const float meanMass = 0.5f * (xiStar.mass() + kKaonMass);
    return std::sqrt(kT * kT + meanMass * meanMass);
  }

  float phiAtTPCAverage(float pt, float phi, int charge, float bz) const
  {
    if (pt <= 0.f || !std::isfinite(pt) || !std::isfinite(phi) || !std::isfinite(bz)) {
      return phi;
    }

    float sum = 0.f;
    int n = 0;

    const auto nRadiusSteps = static_cast<int>(
      (cprRadiusMaxCm.value - cprRadiusMinCm.value) / cprRadiusStepCm.value + 0.5f);

    for (int iRadius = 0; iRadius <= nRadiusSteps; ++iRadius) {
      const auto radiusCm = cprRadiusMinCm.value +
                            static_cast<float>(iRadius) * cprRadiusStepCm.value;

      const float arg = 0.3f * static_cast<float>(charge) * bz * radiusCm * 0.01f / (2.f * pt);
      if (std::abs(arg) < 1.f) {
        sum += phi - std::asin(arg);
        ++n;
      }
    }
    return n > 0 ? sum / static_cast<float>(n) : phi;
  }

  template <typename TXiStar, typename TKaon>
  void fillCPR(bool mixed,
               TXiStar const& xiStar,
               TKaon const& kaon,
               float kstar,
               float bzXiEvent,
               float bzKaonEvent)
  {
    if (!fillCPRQA.value || kstar < 0.f || kstar >= cprKstarMax.value) {
      return;
    }

    const int pionCharge = promptPionCharge(xiStar.channel());
    const float pionPhiTPC = phiAtTPCAverage(xiStar.xiStarPionPt(), xiStar.xiStarPionPhi(), pionCharge, bzXiEvent);
    const float kaonPhiTPC = phiAtTPCAverage(kaon.pt(), kaon.phi(), kaon.charge(), bzKaonEvent);

    const float dEta = kaon.eta() - xiStar.xiStarPionEta();
    const float dPhiStar = wrapToPi(kaonPhiTPC - pionPhiTPC);

    if (mixed) {
      histos.fill(HIST("hCPRME"), dEta, dPhiStar, kstar, static_cast<float>(xiStar.channel()));
    } else {
      histos.fill(HIST("hCPRSE"), dEta, dPhiStar, kstar, static_cast<float>(xiStar.channel()));
    }
  }

  // ---------------------------------------------------------------------------
  // Histogram filling helpers.
  // ---------------------------------------------------------------------------
  template <typename TCollision, typename TXiStar, typename TKaon>
  void fillSameEventPair(TCollision const& collision, TXiStar const& xiStar, TKaon const& kaon)
  {
    histos.fill(HIST("hPairCounterSE"), 1.f);

    if (!passChannelSelection(xiStar.channel()) ||
        !isRequestedChargeCombination(xiStar.channel(), kaon.charge())) {
      return;
    }
    histos.fill(HIST("hPairCounterSE"), 2.f);

    if (pairCuts.applySharedTrackCleaning.value && sharesTrack(xiStar, kaon)) {
      return;
    }
    histos.fill(HIST("hPairCounterSE"), 3.f);

    const float kstar = getKstar(xiStar, kaon);
    if (kstar < 0.f) {
      histos.fill(HIST("hInvalidKstar"), 0.f);
      return;
    }
    histos.fill(HIST("hPairCounterSE"), 4.f);

    const auto channel = static_cast<float>(xiStar.channel());
    const float mt = getPairMt(xiStar, kaon);

    // Getter spelling is fT0MPercentile(), exactly as declared in the data model.
    histos.fill(HIST("SEPairs"), kstar, xiStar.mass(), xiStar.pt(), collision.fT0MPercentile(), channel);
    histos.fill(HIST("hSEKstarVsChannel"), kstar, channel);
    if (fillKinematicQA.value) {
      histos.fill(HIST("hSEKstarVsMtVsChannel"), kstar, mt, channel);
      histos.fill(HIST("hSEKstarVsXiPtVsChannel"), kstar, xiStar.pt(), channel);
      histos.fill(HIST("hSEKstarVsKaonPtVsChannel"), kstar, kaon.pt(), channel);
    }

    if (xiStar.channel() == redxistark::kXiStarUS) {
      histos.fill(HIST("hSEKstarXiStarKminus"), kstar);
    } else if (xiStar.channel() == redxistark::kAntiXiStarUS) {
      histos.fill(HIST("hSEKstarAntiXiStarKplus"), kstar);
    }

    fillCPR(false, xiStar, kaon, kstar, collision.bz(), collision.bz());
  }

  template <typename TCollision1, typename TCollision2, typename TXiStar, typename TKaon>
  void fillMixedEventPair(TCollision1 const& collision1,
                          TCollision2 const& collision2,
                          TXiStar const& xiStar,
                          TKaon const& kaon)
  {
    histos.fill(HIST("hPairCounterME"), 1.f);

    if (!passChannelSelection(xiStar.channel()) ||
        !isRequestedChargeCombination(xiStar.channel(), kaon.charge())) {
      return;
    }
    histos.fill(HIST("hPairCounterME"), 2.f);

    const float kstar = getKstar(xiStar, kaon);
    if (kstar < 0.f) {
      histos.fill(HIST("hInvalidKstar"), 1.f);
      return;
    }
    histos.fill(HIST("hPairCounterME"), 3.f);

    const auto channel = static_cast<float>(xiStar.channel());
    const float mt = getPairMt(xiStar, kaon);

    histos.fill(HIST("MEPairs"), kstar, xiStar.mass(), xiStar.pt(), collision1.fT0MPercentile(), channel);
    histos.fill(HIST("hMEKstarVsChannel"), kstar, channel);
    if (fillKinematicQA.value) {
      histos.fill(HIST("hMEKstarVsMtVsChannel"), kstar, mt, channel);
    }

    if (xiStar.channel() == redxistark::kXiStarUS) {
      histos.fill(HIST("hMEKstarXiStarKminus"), kstar);
    } else if (xiStar.channel() == redxistark::kAntiXiStarUS) {
      histos.fill(HIST("hMEKstarAntiXiStarKplus"), kstar);
    }

    fillCPR(true, xiStar, kaon, kstar, collision1.bz(), collision2.bz());
  }

  // ---------------------------------------------------------------------------
  // Initialization.
  // ---------------------------------------------------------------------------
  void init(InitContext const&)
  {
    auto inRange = [](int value, int lo, int hi) { return value >= lo && value <= hi; };

    if (!inRange(xiCuts.promptPionPIDMode.value, 0, 2) ||
        !inRange(xiCuts.bachelorPionPIDMode.value, 0, 1) ||
        !inRange(xiCuts.promptPionPIDWP.value, 0, 3) ||
        !inRange(xiCuts.bachelorPionPIDWP.value, 0, 4) ||
        !inRange(xiCuts.lambdaPionPIDWP.value, 0, 4) ||
        !inRange(xiCuts.lambdaProtonPIDWP.value, 0, 3)) {
      LOGF(fatal, "Invalid Xi(1530) daughter PID mode/working point");
    }
    if (xiCuts.bachelorPionPIDMode.value == 0 && xiCuts.bachelorPionPIDWP.value == 2) {
      LOGF(fatal, "Xi-bachelor 4.8-sigma PID working point is available only in TPC-only mode");
    }

    if (!inRange(xiCuts.promptPionRowsWP.value, 0, 2) ||
        !inRange(xiCuts.promptPionDCAzWP.value, 0, 2) ||
        !inRange(xiCuts.v0DaughterRowsWP.value, 0, 3) ||
        !inRange(xiCuts.v0PionDcaPVWP.value, 0, 3) ||
        !inRange(xiCuts.v0ProtonDcaPVWP.value, 0, 3) ||
        !inRange(xiCuts.v0DcaPVWP.value, 0, 2) ||
        !inRange(xiCuts.v0DcaDaughtersWP.value, 0, 2) ||
        !inRange(xiCuts.v0CosPAWP.value, 0, 4) ||
        !inRange(xiCuts.v0RadiusWP.value, 0, 4) ||
        !inRange(xiCuts.lambdaMassWP.value, 0, 3) ||
        !inRange(xiCuts.v0LifetimeWP.value, 0, 2) ||
        !inRange(xiCuts.cascBachelorDcaPVWP.value, 0, 2) ||
        !inRange(xiCuts.cascDcaDaughtersWP.value, 0, 2) ||
        !inRange(xiCuts.cascCosPAWP.value, 0, 3) ||
        !inRange(xiCuts.cascRadiusWP.value, 0, 3) ||
        !inRange(xiCuts.xiMassWP.value, 0, 2) ||
        !inRange(xiCuts.xiLifetimeWP.value, 0, 2)) {
      LOGF(fatal, "Invalid Xi(1530) topology/quality working-point index");
    }

    if (!inRange(kaonCuts.dcaXYWP.value, 0, 2) ||
        !inRange(kaonCuts.dcaZWP.value, 0, 2) ||
        !inRange(kaonCuts.pidMode.value, 0, 2)) {
      LOGF(fatal, "Kaon DCA working points and PID mode must be 0, 1 or 2");
    }
    if (kaonCuts.pidPtThreshold.value < 0.f ||
        kaonCuts.tpcNSigmaMax.value <= 0.f ||
        kaonCuts.circularNSigmaMax.value <= 0.f) {
      LOGF(fatal, "Invalid kaon PID configuration");
    }
    if (nEvtMixing.value < 1) {
      LOGF(fatal, "nEvtMixing must be >= 1");
    }
    if (cprRadiusStepCm.value <= 0.f || cprRadiusMaxCm.value < cprRadiusMinCm.value ||
        cprKstarMax.value <= 0.f) {
      LOGF(fatal, "Invalid CPR configuration");
    }

    const AxisSpec axisKstar{cfgAxisKstar, "#it{k}^{*} (GeV/#it{c})"};
    const AxisSpec axisMass{cfgAxisMass, "#it{M}_{#Xi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec axisFT0M{cfgMixFT0MBins, "FT0M percentile (%)"};
    const AxisSpec axisChannel{5, -2.5, 2.5, "#Xi#pi channel"};
    const AxisSpec axisMt{cfgAxisMt, "#it{m}_{T} (GeV/#it{c}^{2})"};
    const AxisSpec axisXiPt{cfgAxisXiPt, "#it{p}_{T}^{#Xi^{*}} (GeV/#it{c})"};
    const AxisSpec axisKaonPt{cfgAxisKaonPt, "#it{p}_{T}^{K} (GeV/#it{c})"};
    const AxisSpec axisPIDP{cfgAxisPIDP, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPIDNSigma{cfgAxisPIDNSigma, "n#sigma"};

    histos.add("hEventCounterSE", "SE event flow;;events", HistType::kTH1F, {{4, 0.5, 4.5}});
    auto hEventSE = histos.get<TH1>(HIST("hEventCounterSE"));
    hEventSE->GetXaxis()->SetBinLabel(1, "reduced event");
    hEventSE->GetXaxis()->SetBinLabel(2, "selected XiPi");
    hEventSE->GetXaxis()->SetBinLabel(3, "selected K");
    hEventSE->GetXaxis()->SetBinLabel(4, "XiPi and K");

    histos.add("hPairCounterSE", "SE pair flow;;pairs", HistType::kTH1F, {{4, 0.5, 4.5}});
    auto hPairSE = histos.get<TH1>(HIST("hPairCounterSE"));
    hPairSE->GetXaxis()->SetBinLabel(1, "selected XiPi x K");
    hPairSE->GetXaxis()->SetBinLabel(2, "charge matched");
    hPairSE->GetXaxis()->SetBinLabel(3, "pair cleaned");
    hPairSE->GetXaxis()->SetBinLabel(4, "valid kstar");

    histos.add("hEventCounterME", "ME event-pair flow;;event pairs", HistType::kTH1F, {{2, 0.5, 2.5}});
    auto hEventME = histos.get<TH1>(HIST("hEventCounterME"));
    hEventME->GetXaxis()->SetBinLabel(1, "mixed event pair");
    hEventME->GetXaxis()->SetBinLabel(2, "nonempty selected groups");

    histos.add("hPairCounterME", "ME pair flow;;pairs", HistType::kTH1F, {{3, 0.5, 3.5}});
    auto hPairME = histos.get<TH1>(HIST("hPairCounterME"));
    hPairME->GetXaxis()->SetBinLabel(1, "selected XiPi x K");
    hPairME->GetXaxis()->SetBinLabel(2, "charge matched");
    hPairME->GetXaxis()->SetBinLabel(3, "valid kstar");

    histos.add("hInvalidKstar", "Invalid k*;sample;counts", HistType::kTH1F, {{2, -0.5, 1.5}});
    auto hInvalid = histos.get<TH1>(HIST("hInvalidKstar"));
    hInvalid->GetXaxis()->SetBinLabel(1, "SE");
    hInvalid->GetXaxis()->SetBinLabel(2, "ME");

    histos.add("hSelectedXiStarMassVsChannel", "Selected #Xi#pi candidates;#it{M}_{#Xi#pi} (GeV/#it{c}^{2});channel", HistType::kTH2F, {axisMass, axisChannel});
    histos.add("hSelectedKaonTPCnSigmaVsP", "Selected external K;#it{p} (GeV/#it{c});n#sigma_{TPC}^{K}", HistType::kTH2F, {axisPIDP, axisPIDNSigma});
    histos.add("hSelectedKaonTOFnSigmaVsP", "Selected external K with TOF;#it{p} (GeV/#it{c});n#sigma_{TOF}^{K}", HistType::kTH2F, {axisPIDP, axisPIDNSigma});

    histos.add("SEPairs", "SE Xi(1530)-K candidates", HistType::kTHnSparseF, {axisKstar, axisMass, axisXiPt, axisFT0M, axisChannel});
    histos.add("MEPairs", "ME Xi(1530)-K candidates", HistType::kTHnSparseF, {axisKstar, axisMass, axisXiPt, axisFT0M, axisChannel});

    histos.add("hSEKstarVsChannel", "SE;#it{k}^{*} (GeV/#it{c});channel", HistType::kTH2F, {axisKstar, axisChannel});
    histos.add("hMEKstarVsChannel", "ME;#it{k}^{*} (GeV/#it{c});channel", HistType::kTH2F, {axisKstar, axisChannel});

    histos.add("hSEKstarXiStarKminus", "SE #Xi^{*0}K^{-};#it{k}^{*} (GeV/#it{c});pairs", HistType::kTH1F, {axisKstar});
    histos.add("hSEKstarAntiXiStarKplus", "SE #bar{#Xi}^{*0}K^{+};#it{k}^{*} (GeV/#it{c});pairs", HistType::kTH1F, {axisKstar});
    histos.add("hMEKstarXiStarKminus", "ME #Xi^{*0}K^{-};#it{k}^{*} (GeV/#it{c});pairs", HistType::kTH1F, {axisKstar});
    histos.add("hMEKstarAntiXiStarKplus", "ME #bar{#Xi}^{*0}K^{+};#it{k}^{*} (GeV/#it{c});pairs", HistType::kTH1F, {axisKstar});

    histos.add("hSEKstarVsMtVsChannel", "SE pair kinematics;#it{k}^{*};#it{m}_{T};channel", HistType::kTH3F, {axisKstar, axisMt, axisChannel});
    histos.add("hMEKstarVsMtVsChannel", "ME pair kinematics;#it{k}^{*};#it{m}_{T};channel", HistType::kTH3F, {axisKstar, axisMt, axisChannel});
    histos.add("hSEKstarVsXiPtVsChannel", "SE pair kinematics;#it{k}^{*};#it{p}_{T}^{#Xi^{*}};channel", HistType::kTH3F, {axisKstar, axisXiPt, axisChannel});
    histos.add("hSEKstarVsKaonPtVsChannel", "SE pair kinematics;#it{k}^{*};#it{p}_{T}^{K};channel", HistType::kTH3F, {axisKstar, axisKaonPt, axisChannel});

    const AxisSpec axisDEta{cfgAxisCPRDEta, "#Delta#eta(K,#pi_{#Xi^{*}})"};
    const AxisSpec axisDPhiStar{cfgAxisCPRDPhiStar, "#Delta#varphi^{*}(K,#pi_{#Xi^{*}})"};
    const AxisSpec axisCPRKstar{cfgAxisCPRKstar, "#it{k}^{*} (GeV/#it{c})"};
    histos.add("hCPRSE", "SE close-pair QA", HistType::kTHnSparseF, {axisDEta, axisDPhiStar, axisCPRKstar, axisChannel});
    histos.add("hCPRME", "ME close-pair QA", HistType::kTHnSparseF, {axisDEta, axisDPhiStar, axisCPRKstar, axisChannel});
  }

  // ---------------------------------------------------------------------------
  // Same-event numerator.
  // The framework groups the two daughter tables by RedXiStarKEventId for the
  // current reduced-event iterator.
  // ---------------------------------------------------------------------------
  void processSameEvent(aod::RedXiStarKEvents::iterator const& collision,
                        aod::XiStarCandidates const& xiStars,
                        aod::KaonCandidates const& kaons)
  {
    histos.fill(HIST("hEventCounterSE"), 1.f);

    bool hasSelectedXi = false;
    for (auto const& xiStar : xiStars) {
      if (!passXiStarSelection(xiStar) || !passChannelSelection(xiStar.channel())) {
        continue;
      }
      hasSelectedXi = true;
      if (fillCandidateQA.value) {
        histos.fill(HIST("hSelectedXiStarMassVsChannel"), xiStar.mass(), static_cast<float>(xiStar.channel()));
      }
    }
    if (hasSelectedXi) {
      histos.fill(HIST("hEventCounterSE"), 2.f);
    }

    bool hasSelectedK = false;
    for (auto const& kaon : kaons) {
      if (!passKaonSelection(kaon)) {
        continue;
      }
      hasSelectedK = true;
      if (fillCandidateQA.value) {
        histos.fill(HIST("hSelectedKaonTPCnSigmaVsP"), kaon.p(), kaon.tPCNSigmaKa());
        if (kaon.hasTOF()) {
          histos.fill(HIST("hSelectedKaonTOFnSigmaVsP"), kaon.p(), kaon.tOFNSigmaKa());
        }
      }
    }
    if (hasSelectedK) {
      histos.fill(HIST("hEventCounterSE"), 3.f);
    }
    if (!(hasSelectedXi && hasSelectedK)) {
      return;
    }
    histos.fill(HIST("hEventCounterSE"), 4.f);

    for (auto const& xiStar : xiStars) {
      if (!passXiStarSelection(xiStar) || !passChannelSelection(xiStar.channel())) {
        continue;
      }
      for (auto const& kaon : kaons) {
        if (!passKaonSelection(kaon)) {
          continue;
        }
        fillSameEventPair(collision, xiStar, kaon);
      }
    }
  }
  PROCESS_SWITCH(xi1530kaoncorrelation, processSameEvent, "Process same-event Xi(1530)-K pairs", true);

  // ---------------------------------------------------------------------------
  // Mixed-event denominator.
  // Xi*(event A) is paired with K(event B) using z-vertex, FT0M percentile and
  // Bz mixing bins.  No cross-event shared-track cleaning is applied for now.
  // ---------------------------------------------------------------------------
  void processMixedEvent(aod::RedXiStarKEvents const& collisions,
                         aod::XiStarCandidates const& xiStars,
                         aod::KaonCandidates const& kaons)

  {
    for (auto const& [collision1, collision2] : selfCombinations(mixingBinning, nEvtMixing, -1, collisions, collisions)) {
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      histos.fill(HIST("hEventCounterME"), 1.f);

      auto groupXi = xiStars.sliceBy(xiStarsPerEvent, collision1.globalIndex());
      auto groupK = kaons.sliceBy(kaonsPerEvent, collision2.globalIndex());

      bool hasSelectedXi = false;
      for (auto const& xiStar : groupXi) {
        if (passXiStarSelection(xiStar) && passChannelSelection(xiStar.channel())) {
          hasSelectedXi = true;
          break;
        }
      }

      bool hasSelectedK = false;
      for (auto const& kaon : groupK) {
        if (passKaonSelection(kaon)) {
          hasSelectedK = true;
          break;
        }
      }

      if (!(hasSelectedXi && hasSelectedK)) {
        continue;
      }
      histos.fill(HIST("hEventCounterME"), 2.f);

      for (auto const& xiStar : groupXi) {
        if (!passXiStarSelection(xiStar) || !passChannelSelection(xiStar.channel())) {
          continue;
        }
        for (auto const& kaon : groupK) {
          if (!passKaonSelection(kaon)) {
            continue;
          }
          fillMixedEventPair(collision1, collision2, xiStar, kaon);
        }
      }
    }
  }
  PROCESS_SWITCH(xi1530kaoncorrelation, processMixedEvent, "Process mixed-event Xi(1530)-K pairs", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<xi1530kaoncorrelation>(cfgc)};
}
