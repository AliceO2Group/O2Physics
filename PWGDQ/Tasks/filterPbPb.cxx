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
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGUD/Core/SGSelector.h"

#include "CommonConstants/LHCConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Vertex.h"
#include <Framework/AnalysisDataModel.h>

#include <fairlogger/Logger.h>

#include <cstdint>
#include <string>
#include <vector>

using namespace std;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct DQFilterPbPbTask {
  Produces<aod::DQRapidityGapFilter> eventRapidityGapFilter;
  OutputObj<TH1D> fFilterOutcome{"Filter outcome"};

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

  SGSelector sgSelector;
  SGCutParHolder sgCuts = SGCutParHolder();

  void init(o2::framework::InitContext&)
  {
    // setup the FilterOutcome histogram
    fFilterOutcome.setObject(new TH1D("Filter outcome", "Filter outcome", 8, 0.5, 8.5));
    fFilterOutcome->GetXaxis()->SetBinLabel(1, "Events inspected");
    fFilterOutcome->GetXaxis()->SetBinLabel(2, "!A && !C");
    fFilterOutcome->GetXaxis()->SetBinLabel(3, "!A && C");
    fFilterOutcome->GetXaxis()->SetBinLabel(4, "A && !C");
    fFilterOutcome->GetXaxis()->SetBinLabel(5, "A && C");
    fFilterOutcome->GetXaxis()->SetBinLabel(6, Form("numContrib not in [%d, %d]", fConfigMinNPVCs.value, fConfigMaxNPVCs.value));
    fFilterOutcome->GetXaxis()->SetBinLabel(7, "BC not found");
    fFilterOutcome->GetXaxis()->SetBinLabel(8, "ITS UPC settings used");

    sgCuts.SetNDtcoll(fConfigNDtColl);
    sgCuts.SetMinNBCs(fConfigMinNBCs);
    sgCuts.SetNTracks(fConfigMinNPVCs, fConfigMaxNPVCs);
    sgCuts.SetMaxFITtime(fConfigMaxFITTime);
    sgCuts.SetFITAmpLimits({static_cast<float>((fConfigUseFV0) ? fConfigFV0AmpLimit : -1.),
                            static_cast<float>((fConfigUseFT0) ? fConfigFT0AAmpLimit : -1.),
                            static_cast<float>((fConfigUseFT0) ? fConfigFT0CAmpLimit : -1.),
                            static_cast<float>((fConfigUseFDD) ? fConfigFDDAAmpLimit : -1.),
                            static_cast<float>((fConfigUseFDD) ? fConfigFDDCAmpLimit : -1.)});
  }

  void processUDSGSelector(MyEvents::iterator const& collision, MyBCs const& bcs,
                           aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds, aod::Zdcs& /*zdcs*/)
  {
    uint64_t filter = 0;
    fFilterOutcome->Fill(1, 1.);
    if (!collision.has_foundBC()) {
      fFilterOutcome->Fill(7, 1);
      return;
    }

    auto bc = collision.foundBC_as<MyBCs>();
    auto newbc = bc;
    auto bcRange = udhelpers::compatibleBCs(collision, fConfigNDtColl, bcs, fConfigMinNBCs);
    auto isSGEvent = sgSelector.IsSelected(sgCuts, collision, bcRange, bc);
    int issgevent = isSGEvent.value;
    // Translate SGSelector values to DQEventFilter values
    if (issgevent == sgselector::SingleGapA) {
      filter |= (static_cast<uint64_t>(1) << VarManager::kSingleGapA);
      fFilterOutcome->Fill(3, 1);
    } else if (issgevent == sgselector::SingleGapC) {
      filter |= (static_cast<uint64_t>(1) << VarManager::kSingleGapC);
      fFilterOutcome->Fill(4, 1);
    } else if (issgevent == sgselector::DoubleGap) {
      filter |= (static_cast<uint64_t>(1) << VarManager::kDoubleGap);
      fFilterOutcome->Fill(2, 1);
    } else if (issgevent == sgselector::NoUpc) {
      fFilterOutcome->Fill(5, 1);
    } else if (issgevent == sgselector::TrkOutOfRange) {
      fFilterOutcome->Fill(6, 1);
    }

    // Get closest bc with FIT activity above threshold
    if (isSGEvent.bc) {
      newbc = *(isSGEvent.bc);
    }
    upchelpers::FITInfo fitInfo{};
    udhelpers::getFITinfo(fitInfo, newbc, bcs, ft0s, fv0as, fdds);

    // Record whether UPC settings were used for this event
    if (collision.flags() & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode) {
      filter |= (static_cast<uint64_t>(1) << VarManager::kITSUPCMode);
      fFilterOutcome->Fill(8, 1);
    }

    eventRapidityGapFilter(filter, newbc.globalIndex());
  }

  void processDummy(MyEvents&)
  {
    // do nothing
  }

  PROCESS_SWITCH(DQFilterPbPbTask, processUDSGSelector, "Run filter task with SG selector from UD", false);
  PROCESS_SWITCH(DQFilterPbPbTask, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DQFilterPbPbTask>(cfgc)};
}
