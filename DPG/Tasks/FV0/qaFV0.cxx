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

/// \file   qaFV0.cxx
/// \author Andreas Molander andreas.molander@cern.ch
/// \brief  FV0 QA

#include <bitset>
#include <cstddef>
#include <cstdint>
#include <string>
#include <map>
#include <vector>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/InitContext.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/runDataProcessing.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;

struct fv0Qa {
  // Constants
  constexpr static std::size_t sMaxBC = o2::constants::lhc::LHCMaxBunches;

  // Event selection conditions
  // Caution required if order is changed. See FillConditionSumHistograms macro below
  enum eConditions {
    kAll,
    kSel8,
    kHasFV0,
    kOrATrgFV0,
    kNChanTrgFV0,
    kChargeTrgFV0,
    kAInTrgFV0,
    kAOutTrgFV0,
    kLaserFV0,
    kOutputsAreBlocked,
    kDataIsValid,
    kDataIsNotValid,
    kNConditions
  };

  // Event selection condition names
  // Caution required if names are changed. See FillConditionHistograms macro below (+ post processing scripts)
  inline static std::map<eConditions, const char*> mapConditionNames = {
    {kAll, "All"},
    {kSel8, "sel8"},
    {kHasFV0, "HasFV0"},
    {kOrATrgFV0, "FV0OrA"},
    {kNChanTrgFV0, "FV0NChan"},
    {kChargeTrgFV0, "FV0Charge"},
    {kAInTrgFV0, "FV0AIn"},
    {kAOutTrgFV0, "FV0AOut"},
    {kLaserFV0, "FV0Laser"},
    {kOutputsAreBlocked, "FV0OutputsAreBlocked"},
    {kDataIsValid, "FV0DataIsValid"},
    {kDataIsNotValid, "FV0DataIsNotValid"}};

  // Observables
  enum eObservables {
    kBc,
    kBcFV0,
    kChannelAmplitudeFV0,
    kChannelAmplitudeSumFV0,
    kAmplitudePerChannelFV0,
    kTimeFV0,
    kNumChannelsFV0,
    kMultiplicityFV0,
    kSumAmpVsNumChannelsFV0,
    kMultiplicityVsNumChannelsFV0,
    kChannelStatsFV0,
    kTriggersFV0,
    kTriggersCorrelationFV0,
    kContributors,
    kNObservables
  };

  // Observable names
  inline static std::map<eObservables, const char*> mapObservableNames = {
    {kBc, "CollisionBC"},
    {kBcFV0, "FV0BC"},
    {kChannelAmplitudeFV0, "FV0ChannelAmplitude"},       // FV0 channel amplitude (ADC channels)
    {kChannelAmplitudeSumFV0, "FV0ChannelAmplitudeSum"}, // Sum of FV0 channel amplitudes per event (ADC channels)
    {kAmplitudePerChannelFV0, "FV0AmplitudePerChannel"},
    {kTimeFV0, "FV0Time"},
    {kNumChannelsFV0, "FV0NumChannels"},
    {kMultiplicityFV0, "FV0Multiplicity"},
    {kSumAmpVsNumChannelsFV0, "Fv0SumAmpVsNumChannels"},
    {kMultiplicityVsNumChannelsFV0, "FV0MultiplicityVsNumChannels"},
    {kChannelStatsFV0, "FV0ChannelStats"},
    {kTriggersFV0, "FV0Triggers"},
    {kTriggersCorrelationFV0, "FV0TriggersCorrelation"},
    {kContributors, "Contributors"}};

  // Histogram registry
  HistogramRegistry histograms{"Histograms"};

  void init(InitContext&)
  {
    if (kNConditions != mapConditionNames.size()) {
      LOG(fatal) << "Number of conditions does not match number of condition names";
    }
    if (kNObservables != mapObservableNames.size()) {
      LOG(fatal) << "Number of observables does not match number of observable names";
    }

    const AxisSpec axisEvSelStats{kNConditions, 0, kNConditions};
    const AxisSpec axisBC{sMaxBC, 0, sMaxBC, "BC"};
    const AxisSpec axisChannelsFV0{49, 0, 49, "Channel ID"};
    const AxisSpec axisTimeNS{500, -5, 5, "Time (ns)"};
    const AxisSpec axisADC{4096, 0, 4096, "Amplitude (ADC channels)"};
    const AxisSpec axisADCSum{4096, 0, 4096 * 49, "Amplitude sum (ADC channels)"};
    const AxisSpec axisMultiplicityFV0{4096, 0, 4096 * 49, "Multiplicity"};
    const AxisSpec axisTriggersFV0{9, 0, 9, "FV0 triggers"};
    const AxisSpec axisContributors{5000, 0, 5000, "# contributors"};

    // Histogram for storing event selection statistics
    auto h = histograms.add<TH1>("EventSelectionStats", "Event selection statistics", kTH1F, {axisEvSelStats});
    for (int iCondition = 0; iCondition < kNConditions; iCondition++) {
      h->GetXaxis()->SetBinLabel(iCondition + 1, mapConditionNames[static_cast<eConditions>(iCondition)]);
    }

    // Lambda that creates one histogram for 'observable' per condition
    auto const makeConditionHistos = [&](eObservables observable, HistType histType, const std::vector<AxisSpec>& axes) {
      for (int iCondition = 0; iCondition < kNConditions; iCondition++) {
        const std::string condition = mapConditionNames[static_cast<eConditions>(iCondition)];
        const std::string histoName = mapObservableNames[observable] + std::string("/") + condition;

        // TODO: make nicer
        if (observable == kTriggersCorrelationFV0) {
          auto h = histograms.add<TH2>(histoName.c_str(), histoName.c_str(), histType, axes);
          h->GetXaxis()->SetBinLabel(1, "orA");
          h->GetXaxis()->SetBinLabel(2, "aOut");
          h->GetXaxis()->SetBinLabel(3, "nChan");
          h->GetXaxis()->SetBinLabel(4, "charge");
          h->GetXaxis()->SetBinLabel(5, "aIn");
          h->GetXaxis()->SetBinLabel(6, "laser");
          h->GetXaxis()->SetBinLabel(7, "outputsBlocked");
          h->GetXaxis()->SetBinLabel(8, "dataIsValid");
          h->GetYaxis()->SetBinLabel(9, "dataIsNotValid");

          h->GetYaxis()->SetBinLabel(1, "orA");
          h->GetYaxis()->SetBinLabel(2, "aOut");
          h->GetYaxis()->SetBinLabel(3, "nChan");
          h->GetYaxis()->SetBinLabel(4, "charge");
          h->GetYaxis()->SetBinLabel(5, "aIn");
          h->GetYaxis()->SetBinLabel(6, "laser");
          h->GetYaxis()->SetBinLabel(7, "outputsBlocked");
          h->GetYaxis()->SetBinLabel(8, "dataIsValid");
          h->GetYaxis()->SetBinLabel(9, "dataIsNotValid");
        } else if (observable == kTriggersFV0) {
          auto h = histograms.add<TH1>(histoName.c_str(), histoName.c_str(), histType, axes);
          h->GetXaxis()->SetBinLabel(1, "orA");
          h->GetXaxis()->SetBinLabel(2, "aOut");
          h->GetXaxis()->SetBinLabel(3, "nChan");
          h->GetXaxis()->SetBinLabel(4, "charge");
          h->GetXaxis()->SetBinLabel(5, "aIn");
          h->GetXaxis()->SetBinLabel(6, "laser");
          h->GetXaxis()->SetBinLabel(7, "outputsBlocked");
          h->GetXaxis()->SetBinLabel(8, "dataIsValid");
          h->GetYaxis()->SetBinLabel(9, "dataIsNotValid");
        } else {
          histograms.add(histoName, histoName.c_str(), histType, axes);
        }
      }
    };

    makeConditionHistos(kBc, kTH1I, {axisBC});
    makeConditionHistos(kBcFV0, kTH1I, {axisBC});
    makeConditionHistos(kChannelAmplitudeFV0, kTH1F, {axisADC});
    makeConditionHistos(kChannelAmplitudeSumFV0, kTH1F, {axisADCSum});
    makeConditionHistos(kAmplitudePerChannelFV0, kTH2F, {axisChannelsFV0, axisADC});
    makeConditionHistos(kTimeFV0, kTH1F, {axisTimeNS});
    makeConditionHistos(kNumChannelsFV0, kTH1I, {axisChannelsFV0});
    makeConditionHistos(kMultiplicityFV0, kTH1F, {axisMultiplicityFV0});
    makeConditionHistos(kMultiplicityVsNumChannelsFV0, kTH2F, {axisChannelsFV0, axisMultiplicityFV0});
    makeConditionHistos(kChannelStatsFV0, kTH1F, {axisChannelsFV0});
    makeConditionHistos(kTriggersFV0, kTH1I, {axisTriggersFV0});
    makeConditionHistos(kTriggersCorrelationFV0, kTH2I, {axisTriggersFV0, axisTriggersFV0});
    makeConditionHistos(kContributors, kTH1I, {axisContributors});

  } // init

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, aod::FV0As const&, aod::BCs const&)
  {
    histograms.fill(HIST("EventSelectionStats"), kAll);

    const int nContributors = collision.numContrib();

    const auto collisionBC = collision.bc_as<aod::BCs>();
    const int localCollisionBC = collisionBC.globalBC() % sMaxBC;

    const bool sel8 = collision.sel8();
    const bool hasFV0 = collision.has_foundFV0();

    bool trgOrAFV0 = false;
    bool trgNChanFV0 = false;
    bool trgChargeFV0 = false;
    bool trgAInFV0 = false;
    bool trgAOutFV0 = false;
    bool trgLaserFV0 = false;
    bool trgOutputsAreBlockedFV0 = false;
    bool trgDataIsValidFV0 = false;
    bool trgDataIsNotValidFV0 = false;

    if (sel8) {
      histograms.fill(HIST("EventSelectionStats"), kSel8);
    }

    if (hasFV0) {
      histograms.fill(HIST("EventSelectionStats"), kHasFV0);

      const auto fv0 = collision.foundFV0();
      const auto fv0BC = fv0.bc_as<aod::BCs>();
      const int localCollisionBCFV0 = fv0BC.globalBC() % sMaxBC;

      const int nFiredChannelsFV0 = fv0.channel().size();
      const float multiplicityFV0 = collision.multFV0A();
      const float timeFV0 = fv0.time();
      const std::bitset<8> triggersFV0 = fv0.triggerMask();

      trgOrAFV0 = triggersFV0[o2::fit::Triggers::bitA];
      trgNChanFV0 = triggersFV0[o2::fit::Triggers::bitTrgNchan];
      trgChargeFV0 = triggersFV0[o2::fit::Triggers::bitTrgCharge];
      trgAInFV0 = triggersFV0[o2::fit::Triggers::bitAIn];
      trgAOutFV0 = triggersFV0[o2::fit::Triggers::bitAOut];
      trgLaserFV0 = triggersFV0[o2::fit::Triggers::bitLaser];
      trgOutputsAreBlockedFV0 = triggersFV0[o2::fit::Triggers::bitOutputsAreBlocked];
      trgDataIsValidFV0 = triggersFV0[o2::fit::Triggers::bitDataIsValid];
      trgDataIsNotValidFV0 = !trgDataIsValidFV0;

      if (trgOrAFV0) {
        histograms.fill(HIST("EventSelectionStats"), kOrATrgFV0);
      }
      if (trgNChanFV0) {
        histograms.fill(HIST("EventSelectionStats"), kNChanTrgFV0);
      }
      if (trgChargeFV0) {
        histograms.fill(HIST("EventSelectionStats"), kChargeTrgFV0);
      }
      if (trgAInFV0) {
        histograms.fill(HIST("EventSelectionStats"), kAInTrgFV0);
      }
      if (trgAOutFV0) {
        histograms.fill(HIST("EventSelectionStats"), kAOutTrgFV0);
      }
      if (trgLaserFV0) {
        histograms.fill(HIST("EventSelectionStats"), kLaserFV0);
      }
      if (trgOutputsAreBlockedFV0) {
        histograms.fill(HIST("EventSelectionStats"), kOutputsAreBlocked);
      }
      if (trgDataIsValidFV0) {
        histograms.fill(HIST("EventSelectionStats"), kDataIsValid);
      }
      if (trgDataIsNotValidFV0) {
        histograms.fill(HIST("EventSelectionStats"), kDataIsNotValid);
      }

      // Sum of amplitudes per condition
      float totalAmplitudes[kNConditions] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Macro for filling histograms for 'observable' based on conditions
#define FillConditionHistograms(observable, ...)                            \
  histograms.fill(HIST(observable "/All"), __VA_ARGS__);                    \
  if (sel8) {                                                               \
    histograms.fill(HIST(observable "/sel8"), __VA_ARGS__);                 \
  }                                                                         \
  if (hasFV0) {                                                             \
    histograms.fill(HIST(observable "/HasFV0"), __VA_ARGS__);               \
  }                                                                         \
  if (trgOrAFV0) {                                                          \
    histograms.fill(HIST(observable "/FV0OrA"), __VA_ARGS__);               \
  }                                                                         \
  if (trgNChanFV0) {                                                        \
    histograms.fill(HIST(observable "/FV0NChan"), __VA_ARGS__);             \
  }                                                                         \
  if (trgChargeFV0) {                                                       \
    histograms.fill(HIST(observable "/FV0Charge"), __VA_ARGS__);            \
  }                                                                         \
  if (trgAInFV0) {                                                          \
    histograms.fill(HIST(observable "/FV0AIn"), __VA_ARGS__);               \
  }                                                                         \
  if (trgAOutFV0) {                                                         \
    histograms.fill(HIST(observable "/FV0AOut"), __VA_ARGS__);              \
  }                                                                         \
  if (trgLaserFV0) {                                                        \
    histograms.fill(HIST(observable "/LaserFV0"), __VA_ARGS__);             \
  }                                                                         \
  if (trgOutputsAreBlockedFV0) {                                            \
    histograms.fill(HIST(observable "/FV0OutputsAreBlocked"), __VA_ARGS__); \
  }                                                                         \
  if (trgDataIsValidFV0) {                                                  \
    histograms.fill(HIST(observable "/FV0DataIsValid"), __VA_ARGS__);       \
  }                                                                         \
  if (trgDataIsNotValidFV0) {                                               \
    histograms.fill(HIST(observable "/FV0DataIsNotValid"), __VA_ARGS__);    \
  }

// Macro for filling histograms for 'observable' with sums (__VA_ARGS__) calculated based on conditions
#define FillConditionSumHistograms(observable, ...)                            \
  histograms.fill(HIST(observable "/All"), __VA_ARGS__[0]);                    \
  if (sel8) {                                                                  \
    histograms.fill(HIST(observable "/sel8"), __VA_ARGS__[1]);                 \
  }                                                                            \
  if (hasFV0) {                                                                \
    histograms.fill(HIST(observable "/HasFV0"), __VA_ARGS__[2]);               \
  }                                                                            \
  if (trgOrAFV0) {                                                             \
    histograms.fill(HIST(observable "/FV0OrA"), __VA_ARGS__[3]);               \
  }                                                                            \
  if (trgNChanFV0) {                                                           \
    histograms.fill(HIST(observable "/FV0NChan"), __VA_ARGS__[4]);             \
  }                                                                            \
  if (trgChargeFV0) {                                                          \
    histograms.fill(HIST(observable "/FV0Charge"), __VA_ARGS__[5]);            \
  }                                                                            \
  if (trgAInFV0) {                                                             \
    histograms.fill(HIST(observable "/FV0AIn"), __VA_ARGS__[6]);               \
  }                                                                            \
  if (trgAOutFV0) {                                                            \
    histograms.fill(HIST(observable "/FV0AOut"), __VA_ARGS__[7]);              \
  }                                                                            \
  if (trgLaserFV0) {                                                           \
    histograms.fill(HIST(observable "/LaserFV0"), __VA_ARGS__[8]);             \
  }                                                                            \
  if (trgOutputsAreBlockedFV0) {                                               \
    histograms.fill(HIST(observable "/FV0OutputsAreBlocked"), __VA_ARGS__[9]); \
  }                                                                            \
  if (trgDataIsValidFV0) {                                                     \
    histograms.fill(HIST(observable "/FV0DataIsValid"), __VA_ARGS__[10]);      \
  }                                                                            \
  if (trgDataIsNotValidFV0) {                                                  \
    histograms.fill(HIST(observable "/FV0DataIsNotValid"), __VA_ARGS__[11]);   \
  }

      // Lambda for filling array of sums based on conditions
      auto const sum = [&](float* sums, float value) {
        sums[kAll] += value;
        if (sel8) {
          sums[kSel8] += value;
        }
        if (hasFV0) {
          sums[kHasFV0] += value;
        }
        if (trgOrAFV0) {
          sums[kOrATrgFV0] += value;
        }
        if (trgNChanFV0) {
          sums[kNChanTrgFV0] += value;
        }
        if (trgChargeFV0) {
          sums[kChargeTrgFV0] += value;
        }
        if (trgAInFV0) {
          sums[kAInTrgFV0] += value;
        }
        if (trgAOutFV0) {
          sums[kAOutTrgFV0] += value;
        }
        if (trgLaserFV0) {
          sums[kLaserFV0] += value;
        }
        if (trgOutputsAreBlockedFV0) {
          sums[kOutputsAreBlocked] += value;
        }
        if (trgDataIsValidFV0) {
          sums[kDataIsValid] += value;
        }
        if (trgDataIsNotValidFV0) {
          sums[kDataIsNotValid] += value;
        }
      };

      FillConditionHistograms("FV0BC", localCollisionBCFV0);

      for (int i = 0; i < fv0.amplitude().size(); i++) {
        FillConditionHistograms("FV0ChannelAmplitude", fv0.amplitude()[i]);
        FillConditionHistograms("FV0AmplitudePerChannel", fv0.channel()[i], fv0.amplitude()[i]);
        sum(totalAmplitudes, fv0.amplitude()[i]);
      }

      FillConditionSumHistograms("FV0ChannelAmplitudeSum", totalAmplitudes);

      FillConditionHistograms("FV0Time", timeFV0);
      FillConditionHistograms("FV0NumChannels", nFiredChannelsFV0);
      FillConditionHistograms("FV0Multiplicity", multiplicityFV0);
      FillConditionHistograms("FV0MultiplicityVsNumChannels", static_cast<float>(nFiredChannelsFV0), multiplicityFV0);
    } // if (collision.has_foundFV0())

    FillConditionHistograms("CollisionBC", localCollisionBC);
    FillConditionHistograms("Contributors", nContributors);

#undef FillConditionHistograms
#undef FillConditionSumHistograms
  } // process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<fv0Qa>(cfgc, TaskName{"fv0-qa"})};
}
