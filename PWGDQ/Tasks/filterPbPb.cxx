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
//
#include <Framework/AnalysisDataModel.h>
#include <fairlogger/Logger.h>
#include <cstdint>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "CommonConstants/LHCConstants.h"

using namespace std;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct DQFilterPbPbTask {
  Produces<aod::DQEventFilter> eventFilter;
  OutputObj<TH1I> fStats{"Statistics"};
  OutputObj<TH1I> fFilterOutcome{"Filter outcome"};

  Configurable<std::string> fConfigEventTypes{"cfgEventTypes", "doublegap,singlegap", "Which event types to select. doublegap, singlegap or both, comma separated"};
  Configurable<int> fConfigNDtColl{"cfgNDtColl", 4, "Number of standard deviations to consider in BC range"};
  Configurable<int> fConfigMinNBCs{"cfgMinNBCs", 7, "Minimum number of BCs to consider in BC range"};
  Configurable<int> fConfigMinNPVCs{"cfgMinNPVCs", 2, "Minimum number of PV contributors"};
  Configurable<int> fConfigMaxNPVCs{"cfgMaxNPVCs", 5, "Maximum number of PV contributors"};
  Configurable<float> fConfigMaxFITTime{"cfgMaxFITTime", 4, "Maximum time in FIT"};
  Configurable<float> fConfigFV0AmpLimit{"cfgFV0AmpLimit", 0, "FV0 amplitude limit for event selection"};
  Configurable<float> fConfigFT0AAmpLimit{"cfgFT0AAmpLimit", 0, "FT0A amplitude limit for event selection"};
  Configurable<float> fConfigFT0CAmpLimit{"cfgFT0CAmpLimit", 0, "FT0C amplitude limit for event selection"};
  Configurable<float> fConfigFDDAAmpLimit{"cfgFDDAAmpLimit", 0, "FDDA amplitude limit for event selection"};
  Configurable<float> fConfigFDDCAmpLimit{"cfgFDDCAmpLimit", 0, "FDDC amplitude limit for event selection"};
  Configurable<bool> fConfigUseFV0{"cfgUseFV0", true, "Whether to use FV0 for veto"};
  Configurable<bool> fConfigUseFT0{"cfgUseFT0", true, "Whether to use FT0 for veto"};
  Configurable<bool> fConfigUseFDD{"cfgUseFDD", true, "Whether to use FDD for veto"};

  int eventTypeMap = 0;
  std::vector<std::string> eventTypeOptions = {"doublegap", "singlegap"}; // Map for which types of event to select

  void init(o2::framework::InitContext&)
  {
    // setup the Stats histogram
    fStats.setObject(new TH1I("Statistics", "Stats for DQ triggers", 2, -2.5, -0.5));
    fStats->GetXaxis()->SetBinLabel(1, "Events inspected");
    fStats->GetXaxis()->SetBinLabel(2, "Events selected");
    // setup the FilterOutcome histogram
    fFilterOutcome.setObject(new TH1I("Filter outcome", "Filter outcome", 5, -0.5, 4.5));
    fFilterOutcome->GetXaxis()->SetBinLabel(1, "Events inspected");
    fFilterOutcome->GetXaxis()->SetBinLabel(2, "Events selected");
    fFilterOutcome->GetXaxis()->SetBinLabel(3, "Failed FIT veto");
    fFilterOutcome->GetXaxis()->SetBinLabel(4, Form("numContrib not in [%d, %d]", fConfigMinNPVCs.value, fConfigMaxNPVCs.value));
    fFilterOutcome->GetXaxis()->SetBinLabel(5, "BC not found");

    TString eventTypesString = fConfigEventTypes.value;
    for (std::vector<std::string>::size_type i = 0; i < eventTypeOptions.size(); i++) {
      if (eventTypesString.Contains(eventTypeOptions[i])) {
        eventTypeMap |= (uint32_t(1) << i);
        LOGF(info, "filterPbPb will select '%s' events", eventTypeOptions[i]);
      }
    }
    if (eventTypeMap == 0) {
      LOGF(fatal, "No valid choice of event types to select. Use 'doublegap', 'singlegap' or both");
    }
  }

  // Helper function for selecting double gap and single gap events
  template <typename TEvent, typename TBCs>
  uint64_t rapidityGapFilter(TEvent const& collision, TBCs const& bcs,
                             aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
                             int eventTypes, std::vector<float> FITAmpLimits, int nDtColl, int minNBCs, int minNPVCs, int maxNPVCs, float maxFITTime,
                             bool useFV0, bool useFT0, bool useFDD)
  {
    fFilterOutcome->Fill(0., 1.);
    // Find BC associated with collision
    if (!collision.has_foundBC()) {
      fFilterOutcome->Fill(4., 1);
      return 0;
    }
    // foundBCId is stored in EvSels
    auto bc = collision.template foundBC_as<TBCs>();

    // Obtain slice of compatible BCs
    uint64_t mostProbableBC = bc.globalBC();
    uint64_t meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);
    int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * nDtColl);
    if (deltaBC < minNBCs) {
      deltaBC = minNBCs;
    }

    // Range of BCs to consider
    uint64_t minBC = (uint64_t)deltaBC < meanBC ? meanBC - (uint64_t)deltaBC : 0;
    uint64_t maxBC = meanBC + (uint64_t)deltaBC;

    int slicemin = 0;
    int slicemax = 0;
    // Check if there is overlap between acceptable and possible BCs
    if (maxBC > bcs.iteratorAt(0).globalBC() && minBC < bcs.iteratorAt(bcs.size() - 1).globalBC()) {
      // find slice of BCs table with BC in [minBC, maxBC]
      int moveCount = 0;
      int64_t minBCId = bc.globalIndex();
      int64_t maxBCId = bc.globalIndex();
      // lower limit
      if (bc.globalBC() < minBC) {
        while (bc != bcs.end() && bc.globalBC() < minBC) {
          ++bc;
          ++moveCount;
          minBCId = bc.globalIndex();
        }
      } else {
        while (bc.globalIndex() > 0 && bc.globalBC() >= minBC) {
          minBCId = bc.globalIndex();
          --bc;
          --moveCount;
        }
      }
      // upper limit
      if (bc.globalBC() < maxBC) {
        while (bc != bcs.end() && bc.globalBC() <= maxBC) {
          maxBCId = bc.globalIndex();
          ++bc;
          ++moveCount;
        }
      } else {
        while (bc.globalIndex() > 0 && bc.globalBC() > maxBC) {
          --bc;
          --moveCount;
          maxBCId = bc.globalIndex();
        }
      }
      // reset bc
      bc.moveByIndex(-moveCount);
      // Create BC slice
      slicemin = minBCId;
      slicemax = maxBCId - minBCId + 1;
    }
    MyBCs bcrange{{bcs.asArrowTable()->Slice(slicemin, slicemax)}, (uint16_t)slicemin};
    // Rapidity gap condition: Check FIT activity in BC range
    bool isSideAClean = true;
    bool isSideCClean = true;

    for (auto const& bc : bcrange) {
      if (useFV0) {
        if (bc.has_foundFV0()) {
          auto fv0a = fv0as.iteratorAt(bc.foundFV0Id());
          float FV0Amplitude = 0;
          for (auto amp : fv0a.amplitude()) {
            FV0Amplitude += amp;
          }
          float FV0Time = std::abs(fv0a.time());
          if (FV0Amplitude > FITAmpLimits[0] || FV0Time > maxFITTime) {
            isSideAClean = false;
          }
        }
      }
      if (useFT0) {
        if (bc.has_foundFT0()) {
          auto ft0 = ft0s.iteratorAt(bc.foundFT0Id());
          float FT0AAmplitude = 0;
          float FT0CAmplitude = 0;
          for (auto amp : ft0.amplitudeA()) {
            FT0AAmplitude += amp;
          }
          for (auto amp : ft0.amplitudeC()) {
            FT0CAmplitude += amp;
          }
          float FT0ATime = std::abs(ft0.timeA());
          float FT0CTime = std::abs(ft0.timeC());
          if (FT0AAmplitude > FITAmpLimits[1] || FT0ATime > maxFITTime) {
            isSideAClean = false;
          }
          if (FT0CAmplitude > FITAmpLimits[2] || FT0CTime > maxFITTime) {
            isSideCClean = false;
          }
        }
      }
      if (useFDD) {
        if (bc.has_foundFDD()) {
          auto fdd = fdds.iteratorAt(bc.foundFDDId());
          float FDDAAmplitude = 0;
          float FDDCAmplitude = 0;
          for (auto amp : fdd.chargeA()) {
            FDDAAmplitude += amp;
          }
          for (auto amp : fdd.chargeC()) {
            FDDCAmplitude += amp;
          }
          float FDDATime = std::abs(fdd.timeA());
          float FDDCTime = std::abs(fdd.timeC());
          if (FDDAAmplitude > FITAmpLimits[3] || FDDATime > maxFITTime) {
            isSideAClean = false;
          }
          if (FDDCAmplitude > FITAmpLimits[4] || FDDCTime > maxFITTime) {
            isSideCClean = false;
          }
        }
      }
    }

    // Compute FIT decision
    uint64_t FITDecision = 0;
    if (eventTypes & (uint32_t(1) << 0)) {
      if (isSideAClean && isSideCClean) {
        FITDecision |= (uint64_t(1) << VarManager::kDoubleGap);
      }
    }
    if (eventTypes & (uint32_t(1) << 1)) {
      if (isSideAClean && !isSideCClean) {
        FITDecision |= (uint64_t(1) << VarManager::kSingleGapA);
      } else if (!isSideAClean && isSideCClean) {
        FITDecision |= (uint64_t(1) << VarManager::kSingleGapC);
      }
    }
    if (!FITDecision) {
      fFilterOutcome->Fill(2, 1);
      return 0;
    }

    // Number of primary vertex contributors
    if (collision.numContrib() < minNPVCs || collision.numContrib() > maxNPVCs) {
      fFilterOutcome->Fill(3, 1);
      return 0;
    }

    // If we made it here, the event passed
    fFilterOutcome->Fill(1, 1);
    // Return filter bitmap corresponding to FIT decision
    return FITDecision;
  }

  void processFilterPbPb(MyEvents::iterator const& collision, MyBCs const& bcs,
                         aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    fStats->Fill(-2.0);

    std::vector<float> FITAmpLimits = {fConfigFV0AmpLimit, fConfigFT0AAmpLimit, fConfigFT0CAmpLimit, fConfigFDDAAmpLimit, fConfigFDDCAmpLimit};

    uint64_t filter = rapidityGapFilter(collision, bcs, ft0s, fv0as, fdds,
                                        eventTypeMap, FITAmpLimits, fConfigNDtColl, fConfigMinNBCs, fConfigMinNPVCs, fConfigMaxNPVCs, fConfigMaxFITTime,
                                        fConfigUseFV0, fConfigUseFT0, fConfigUseFDD);

    bool isSelected = filter;
    fStats->Fill(-1.0, isSelected);

    eventFilter(filter);
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPbPbTask, processFilterPbPb, "Run filter task", true);
  PROCESS_SWITCH(DQFilterPbPbTask, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQFilterPbPbTask>(cfgc)};
}
