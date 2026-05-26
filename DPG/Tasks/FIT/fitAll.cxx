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

/// \file   fitAll.cxx
/// \brief  FITAll table definition and producer. Standalone derived data used for FIT studies.
///
/// \author Andreas Molander andreas.molander@cern.ch

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"

#include <FV0Base/Constants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <FDDBase/Constants.h>

#include <cstdint>
#include <limits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

namespace o2::aod
{
namespace fit
{

// Constants
static constexpr int NchFt0 = 208;                              ///< Number of FT0 channels
static constexpr int NchFt0A = 96;                              ///< Number of FT0A channels
static constexpr int NchFv0 = o2::fv0::Constants::nFv0Channels; ///< Number of FV0 channels
static constexpr int NchFdd = o2::fdd::Nchannels;               ///< Number of FDD channels

// BCs
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(GlobalBc, globalBc, uint64_t);
DECLARE_SOA_COLUMN(CtpTriggerMask, ctpTriggerMask, uint64_t);
DECLARE_SOA_COLUMN(CtpInputMask, ctpInputMask, uint64_t);

// Timestamps
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);

// EvSels
DECLARE_SOA_BITMAP_COLUMN(TriggerAlias, triggerAlias, 32);
DECLARE_SOA_BITMAP_COLUMN(SelectionFlags, selectionFlags, 64);
DECLARE_SOA_BITMAP_COLUMN(RctFlags, rctFlags, 32);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(HasFt0, hasFt0, bool);
DECLARE_SOA_COLUMN(HasFv0, hasFv0, bool);
DECLARE_SOA_COLUMN(HasFdd, hasFdd, bool);
DECLARE_SOA_COLUMN(HasZdc, hasZdc, bool);

// Collisions
DECLARE_SOA_COLUMN(BcId, bcId, int32_t);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t);
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float);

// Mults
DECLARE_SOA_COLUMN(MultFt0A, multFt0A, float);
DECLARE_SOA_COLUMN(MultFt0C, multFt0C, float);
DECLARE_SOA_COLUMN(MultFv0A, multFv0A, float);
DECLARE_SOA_COLUMN(MultFddA, multFddA, float);
DECLARE_SOA_COLUMN(MultFddC, multFddC, float);
DECLARE_SOA_COLUMN(MultZnA, multZnA, float);
DECLARE_SOA_COLUMN(MultZnC, multZnC, float);
DECLARE_SOA_COLUMN(MultZem1, multZem1, float);
DECLARE_SOA_COLUMN(MultZem2, multZem2, float);
DECLARE_SOA_COLUMN(MultZpA, multZpA, float);
DECLARE_SOA_COLUMN(MultZpC, multZpC, float);

// FT0s
DECLARE_SOA_COLUMN(Ft0BcId, ft0BcId, int32_t);
DECLARE_SOA_COLUMN(Ft0AmplitudeA, ft0AmplitudeA, std::vector<float>);
DECLARE_SOA_COLUMN(Ft0ChannelA, ft0ChannelA, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Ft0AmplitudeC, ft0AmplitudeC, std::vector<float>);
DECLARE_SOA_COLUMN(Ft0ChannelC, ft0ChannelC, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Ft0TimeA, ft0TimeA, float);
DECLARE_SOA_COLUMN(Ft0TimeC, ft0TimeC, float);
DECLARE_SOA_COLUMN(Ft0TriggerMask, ft0TriggerMask, uint8_t);
DECLARE_SOA_COLUMN(Ft0PosZ, ft0PosZ, float);
DECLARE_SOA_COLUMN(Ft0CollTime, ft0CollTime, float);
DECLARE_SOA_COLUMN(Ft0SumAmpA, ft0SumAmpA, float);
DECLARE_SOA_COLUMN(Ft0SumAmpC, ft0SumAmpC, float);

// FT0sCorrected
DECLARE_SOA_COLUMN(T0aCorrected, t0aCorrected, float);
DECLARE_SOA_COLUMN(T0cCorrected, t0cCorrected, float);
DECLARE_SOA_COLUMN(T0ac, t0ac, float);
DECLARE_SOA_COLUMN(T0resolution, t0resolution, float);

// FT0 derived quantities
DECLARE_SOA_COLUMN(Ft0ChAmpl, ft0ChAmpl, std::vector<float>); ///< FT0 channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(Ft0TotAmplA, ft0TotAmplA, float);          ///< FT0-A total amplitude computed from channel amplitudes (for cross check)
DECLARE_SOA_COLUMN(Ft0TotAmplC, ft0TotAmplC, float);          ///< FT0-C total amplitude computed from channel amplitudes (for cross check)

// FV0As
DECLARE_SOA_COLUMN(Fv0BcId, fv0BcId, int32_t);
DECLARE_SOA_COLUMN(Fv0Amplitude, fv0Amplitude, std::vector<float>);
DECLARE_SOA_COLUMN(Fv0Channel, fv0Channel, std::vector<uint8_t>);
DECLARE_SOA_COLUMN(Fv0Time, fv0Time, float);
DECLARE_SOA_COLUMN(Fv0TriggerMask, fv0TriggerMask, uint8_t);

// FV0 derived quantities
DECLARE_SOA_COLUMN(Fv0ChAmpl, fv0ChAmpl, std::vector<float>); ///< FV0 channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(Fv0TotAmpl, fv0TotAmpl, float);            ///< FV0 total amplitude computed from channel amplitudes (for cross check)

// FDDs
DECLARE_SOA_COLUMN(FddBcId, fddBcId, int32_t);
DECLARE_SOA_COLUMN(FddChargeA, fddChargeA, int16_t[8]);
DECLARE_SOA_COLUMN(FddChargeC, fddChargeC, int16_t[8]);
DECLARE_SOA_COLUMN(FddTimeA, fddTimeA, float);
DECLARE_SOA_COLUMN(FddTimeC, fddTimeC, float);
DECLARE_SOA_COLUMN(FddTriggerMask, fddTriggerMask, uint8_t);

// FDD derived quantities
DECLARE_SOA_COLUMN(FddChAmpl, fddChAmpl, std::vector<float>); ///< FDD channel amplitudes, vector idx = ch ID
DECLARE_SOA_COLUMN(FddTotAmplA, fddTotAmplA, float);          ///< FDD-A total amplitude computed from channel amplitudes (for cross check)
DECLARE_SOA_COLUMN(FddTotAmplC, fddTotAmplC, float);          ///< FDD-C total amplitude computed from channel amplitudes (for cross check)

// ZDCs
DECLARE_SOA_COLUMN(EnergyCommonZnA, energyCommonZnA, float);
DECLARE_SOA_COLUMN(EnergyCommonZnC, energyCommonZnC, float);
DECLARE_SOA_COLUMN(TimeZem1, timeZem1, float);
DECLARE_SOA_COLUMN(TimeZem2, timeZem2, float);

} // namespace fit

DECLARE_SOA_TABLE(FITsAll, "AOD", "FITALL", ///< Standalone derived data used for FIT studies.
                  fit::RunNumber, fit::GlobalBc, fit::CtpTriggerMask, fit::CtpInputMask,
                  fit::Timestamp,
                  fit::TriggerAlias, fit::SelectionFlags, fit::RctFlags, fit::Sel8, fit::HasFt0, fit::HasFv0, fit::HasFdd, fit::HasZdc,
                  fit::BcId, fit::PosX, fit::PosY, fit::PosZ, fit::Flags, fit::NumContrib, fit::CollisionTime, fit::CollisionTimeRes,
                  fit::MultFt0A, fit::MultFt0C, fit::MultFv0A, fit::MultFddA, fit::MultFddC, fit::MultZnA, fit::MultZnC, fit::MultZem1, fit::MultZem2, fit::MultZpA, fit::MultZpC,
                  fit::Ft0BcId, fit::Ft0AmplitudeA, fit::Ft0ChannelA, fit::Ft0AmplitudeC, fit::Ft0ChannelC,
                  fit::Ft0TimeA, fit::Ft0TimeC, fit::Ft0TriggerMask, fit::Ft0PosZ, fit::Ft0CollTime, fit::Ft0SumAmpA, fit::Ft0SumAmpC,
                  fit::T0aCorrected, fit::T0cCorrected, fit::T0ac, fit::T0resolution,
                  fit::Ft0ChAmpl, fit::Ft0TotAmplA, fit::Ft0TotAmplC,
                  fit::Fv0BcId, fit::Fv0Amplitude, fit::Fv0Channel, fit::Fv0Time, fit::Fv0TriggerMask,
                  fit::Fv0ChAmpl, fit::Fv0TotAmpl,
                  fit::FddBcId, fit::FddChargeA, fit::FddChargeC, fit::FddTimeA, fit::FddTimeC, fit::FddTriggerMask,
                  fit::FddChAmpl, fit::FddTotAmplA, fit::FddTotAmplC,
                  fit::EnergyCommonZnA, fit::EnergyCommonZnC, fit::TimeZem1, fit::TimeZem2);

using FITAll = FITsAll::iterator;

} // namespace o2::aod

struct FitAll {
  // Producer
  Produces<o2::aod::FITsAll> table;

  void init(InitContext const&)
  {
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::FT0sCorrected> const& collisions,
               aod::BCsWithTimestamps const&,
               aod::FT0s const&, aod::FV0As const&, aod::FDDs const&, aod::Zdcs const&)
  {
    table.reserve(collisions.size());

    // BC
    int runNumber = -1;
    uint64_t globalBC = 0;
    uint64_t ctpTriggerMask = 0;
    uint64_t ctpInputMask = 0;

    // Timestamp
    uint64_t timestamp = 0;

    // EvSel
    uint32_t triggerAlias = 0;
    uint64_t selectionFlags = 0;
    uint32_t rctFlags = 0;
    bool sel8 = false;
    bool hasFT0 = false;
    bool hasFV0 = false;
    bool hasFDD = false;
    bool hasZDC = false;

    // Collision
    int32_t bcId = -1;
    float posX = -200;
    float posY = -200;
    float posZ = -200;
    uint16_t flags = 0;
    uint16_t numContrib = 0;
    float collisionTime = -200;
    float collisionTimeRes = -200;

    // Mult
    float multFT0A = 0;
    float multFT0C = 0;
    float multFV0A = 0;
    float multFDDA = 0;
    float multFDDC = 0;
    float multZNA = 0;
    float multZNC = 0;
    float multZEM1 = 0;
    float multZEM2 = 0;
    float multZPA = 0;
    float multZPC = 0;

    // FT0
    int32_t ft0BCId = -1;
    std::vector<float> ft0AmplitudeA;
    std::vector<uint8_t> ft0ChannelA;
    std::vector<float> ft0AmplitudeC;
    std::vector<uint8_t> ft0ChannelC;
    float ft0TimeA = -200;
    float ft0TimeC = -200;
    uint8_t ft0TriggerMask = 0;
    float ft0PosZ = -200;
    float ft0CollTime = -200;
    float ft0SumAmpA = 0;
    float ft0SumAmpC = 0;

    // FT0Corrected
    float t0ACorrected = -200;
    float t0CCorrected = -200;
    float t0AC = -200;
    float t0resolution = -200;

    // FT0 derived quantities
    std::vector<float> ft0ChAmpl(o2::aod::fit::NchFt0, 0);
    float ft0TotAmplA = 0;
    float ft0TotAmplC = 0;

    // FV0A
    int32_t fv0BCId = -1;
    std::vector<float> fv0Amplitude;
    std::vector<uint8_t> fv0Channel;
    float fv0Time = -200;
    uint8_t fv0TriggerMask = 0;

    // FV0 derived quantities
    std::vector<float> fv0ChAmpl(o2::aod::fit::NchFv0, 0);
    float fv0TotAmpl = 0;

    // FDD
    int32_t fddBCId = -1;
    int16_t fddChargeA[8] = {0};
    int16_t fddChargeC[8] = {0};
    float fddTimeA = -200;
    float fddTimeC = -200;
    uint8_t fddTriggerMask = 0;

    // FDD derived quantities
    std::vector<float> fddChAmpl(o2::aod::fit::NchFdd, 0);
    float fddTotAmplA = 0;
    float fddTotAmplC = 0;

    // ZDC
    float energyCommonZNA = -std::numeric_limits<float>::infinity();
    float energyCommonZNC = -std::numeric_limits<float>::infinity();
    float timeZEM1 = -std::numeric_limits<float>::infinity();
    float timeZEM2 = -std::numeric_limits<float>::infinity();

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();

      // BC
      runNumber = bc.runNumber();
      globalBC = bc.globalBC();
      ctpTriggerMask = bc.triggerMask();
      ctpInputMask = bc.inputMask();

      // Timestamp
      timestamp = bc.timestamp();

      // EvSel
      triggerAlias = collision.alias_raw();
      selectionFlags = collision.selection_raw();
      rctFlags = collision.rct_raw();
      sel8 = collision.sel8();
      hasFT0 = collision.has_foundFT0();
      hasFV0 = collision.has_foundFV0();
      hasFDD = collision.has_foundFDD();
      hasZDC = collision.has_foundZDC();

      // Collision
      bcId = collision.bcId();
      posX = collision.posX();
      posY = collision.posY();
      posZ = collision.posZ();
      flags = collision.flags();
      numContrib = collision.numContrib();
      collisionTime = collision.collisionTime();
      collisionTimeRes = collision.collisionTimeRes();

      // Multiplicities
      multFT0A = collision.multFT0A();
      multFT0C = collision.multFT0C();
      multFV0A = collision.multFV0A();
      multFDDA = collision.multFDDA();
      multFDDC = collision.multFDDC();
      multZNA = collision.multZNA();
      multZNC = collision.multZNC();
      multZEM1 = collision.multZEM1();
      multZEM2 = collision.multZEM2();
      multZPA = collision.multZPA();
      multZPC = collision.multZPC();

      // FT0
      ft0BCId = -1;
      ft0AmplitudeA.clear();
      ft0ChannelA.clear();
      ft0AmplitudeC.clear();
      ft0ChannelC.clear();
      ft0TimeA = -200;
      ft0TimeC = -200;
      ft0TriggerMask = 0;
      ft0PosZ = -200;
      ft0CollTime = -200;
      ft0SumAmpA = 0;
      ft0SumAmpC = 0;

      // FT0Corrected
      t0ACorrected = collision.t0ACorrected();
      t0CCorrected = collision.t0CCorrected();
      t0AC = collision.t0AC();
      t0resolution = collision.t0resolution();

      // FT0 derived quantities
      ft0ChAmpl.assign(o2::aod::fit::NchFt0, 0);
      ft0TotAmplA = 0;
      ft0TotAmplC = 0;

      // FV0A
      fv0BCId = -1;
      fv0Amplitude.clear();
      fv0Channel.clear();
      fv0Time = -200;
      fv0TriggerMask = 0;

      // FV0 derived quantities
      fv0ChAmpl.assign(o2::aod::fit::NchFv0, 0);
      fv0TotAmpl = 0;

      // FDD
      fddBCId = -1;
      std::fill(std::begin(fddChargeA), std::end(fddChargeA), 0);
      std::fill(std::begin(fddChargeC), std::end(fddChargeC), 0);
      fddTimeA = -200;
      fddTimeC = -200;
      fddTriggerMask = 0;

      // FDD derived quantities
      fddChAmpl.assign(o2::aod::fit::NchFdd, 0);
      fddTotAmplA = 0;
      fddTotAmplC = 0;

      // ZDC
      energyCommonZNA = -std::numeric_limits<float>::infinity();
      energyCommonZNC = -std::numeric_limits<float>::infinity();
      timeZEM1 = -std::numeric_limits<float>::infinity();
      timeZEM2 = -std::numeric_limits<float>::infinity();

      if (hasFT0) {
        // FT0
        auto ft0 = collision.foundFT0();
        ft0BCId = ft0.bcId();
        ft0TimeA = ft0.timeA();
        ft0TimeC = ft0.timeC();
        ft0TriggerMask = ft0.triggerMask();
        ft0PosZ = ft0.posZ();
        ft0CollTime = ft0.collTime();
        ft0SumAmpA = ft0.sumAmpA();
        ft0SumAmpC = ft0.sumAmpC();

        for (size_t i = 0; i < ft0.amplitudeA().size(); i++) {
          // FT0
          ft0AmplitudeA.push_back(ft0.amplitudeA()[i]);
          ft0ChannelA.push_back(ft0.channelA()[i]);

          // FT0 derived quantities
          ft0ChAmpl[ft0.channelA()[i]] = ft0.amplitudeA()[i];
          ft0TotAmplA += ft0.amplitudeA()[i];
        }
        for (size_t i = 0; i < ft0.amplitudeC().size(); i++) {
          // FT0
          ft0AmplitudeC.push_back(ft0.amplitudeC()[i]);
          ft0ChannelC.push_back(ft0.channelC()[i]);

          // FT0 derived quantities
          ft0ChAmpl[ft0.channelC()[i] + o2::aod::fit::NchFt0A] = ft0.amplitudeC()[i]; // Channel IDs in the C-side array start from zero in AO2D (JIRA AFIT-129)
          ft0TotAmplC += ft0.amplitudeC()[i];
        }
      }

      if (hasFV0) {
        // FV0A
        auto fv0 = collision.foundFV0();
        fv0BCId = fv0.bcId();
        fv0Time = fv0.time();
        fv0TriggerMask = fv0.triggerMask();

        for (size_t i = 0; i < fv0.amplitude().size(); i++) {
          // FV0A
          fv0Amplitude.push_back(fv0.amplitude()[i]);
          fv0Channel.push_back(fv0.channel()[i]);

          // FV0 derived quantities
          fv0ChAmpl[fv0.channel()[i]] = fv0.amplitude()[i];
          fv0TotAmpl += fv0.amplitude()[i];
        }
      }

      if (hasFDD) {
        // FDD
        auto fdd = collision.foundFDD();
        fddBCId = fdd.bcId();
        fddTimeA = fdd.timeA();
        fddTimeC = fdd.timeC();
        fddTriggerMask = fdd.triggerMask();

        for (size_t i = 0; i < o2::aod::fit::NchFdd / 2; i++) {
          // FDD
          fddChargeA[i] = fdd.chargeA()[i];
          fddChargeC[i] = fdd.chargeC()[i];

          // FDD derived quantities
          fddChAmpl[i + o2::aod::fit::NchFdd / 2] = fdd.chargeA()[i];
          fddChAmpl[i] = fdd.chargeC()[i];
          fddTotAmplA += fdd.chargeA()[i];
          fddTotAmplC += fdd.chargeC()[i];
        }
      }

      if (hasZDC) {
        // ZDC
        auto zdc = collision.foundZDC();
        energyCommonZNA = zdc.energyCommonZNA();
        energyCommonZNC = zdc.energyCommonZNC();
        timeZEM1 = zdc.timeZEM1();
        timeZEM2 = zdc.timeZEM2();
      }

      table(runNumber, globalBC, ctpTriggerMask, ctpInputMask,
            timestamp,
            triggerAlias, selectionFlags, rctFlags, sel8, hasFT0, hasFV0, hasFDD, hasZDC,
            bcId, posX, posY, posZ, flags, numContrib, collisionTime, collisionTimeRes,
            multFT0A, multFT0C, multFV0A, multFDDA, multFDDC, multZNA, multZNC, multZEM1, multZEM2, multZPA, multZPC,
            ft0BCId, ft0AmplitudeA, ft0ChannelA, ft0AmplitudeC, ft0ChannelC,
            ft0TimeA, ft0TimeC, ft0TriggerMask, ft0PosZ, ft0CollTime, ft0SumAmpA, ft0SumAmpC,
            t0ACorrected, t0CCorrected, t0AC, t0resolution,
            ft0ChAmpl, ft0TotAmplA, ft0TotAmplC,
            fv0BCId, fv0Amplitude, fv0Channel, fv0Time, fv0TriggerMask,
            fv0ChAmpl, fv0TotAmpl,
            fddBCId, fddChargeA, fddChargeC, fddTimeA, fddTimeC, fddTriggerMask,
            fddChAmpl, fddTotAmplA, fddTotAmplC,
            energyCommonZNA, energyCommonZNC, timeZEM1, timeZEM2);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FitAll>(cfgc)};
  return workflow;
}
