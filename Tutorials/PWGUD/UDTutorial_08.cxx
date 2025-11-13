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
// \brief UD tutorial demonstrating event selection using SGSelector. Processes raw AO2Ds. Dependency: o2-analysis-event-selection-service
// \author Sigurd Nese
// \since  November 2025

#include "PWGUD/Core/SGSelector.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2::framework;

// EvSels table contains connection between collision and BC
using MyEvents = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
// BcSels table contains connection between BC and FIT info. Run3MatchedToBCSparse table contains index into ZDC table.
using MyBCs = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels, o2::aod::Run3MatchedToBCSparse>;

struct UDTutorial08 {

  // Histogram setup
  OutputObj<TH1I> hSelectionResult{TH1I("hSelectionResult", "hSelectionResult", 5, -0.5, 4.5)};
  OutputObj<TH1I> hSelectionResultAfterCut{TH1I("hSelectionResultAfterCut", "hSelectionResultAfterCut", 5, -0.5, 4.5)};
  OutputObj<TH1F> hZNAEnergy{TH1F("hZNAEnergy", "hZNAEnergy", 200, 0, 20)};
  OutputObj<TH1F> hZNAEnergyAfterCut{TH1F("hZNAEnergyAfterCut", "hZNAEnergyAfterCut", 200, 0, 20)};
  OutputObj<TH1F> hZNCEnergy{TH1F("hZNCEnergy", "hZNCEnergy", 200, 0, 20)};
  OutputObj<TH1F> hZNCEnergyAfterCut{TH1F("hZNCEnergyAfterCut", "hZNCEnergyAfterCut", 200, 0, 20)};
  // Create instance of the selector class which runs the gap selection algorithm
  SGSelector sgSelector;
  // Create instance of cut holder class to contain the user defined cuts
  SGCutParHolder sgCuts = SGCutParHolder();

  void init(o2::framework::InitContext&)
  {
    // Configure the gap selection criteria. Rest of the values are kept default
    sgCuts.SetNDtcoll(1);         // Time range to consider, in units of collision time resolution
    sgCuts.SetMinNBCs(2);         // Minimum number of BCs to check
    sgCuts.SetNTracks(2, 100);    // Reject collisions with < 2 tracks and > 100 tracks
    sgCuts.SetMaxFITtime(34);     // Reject collisions with FIT time > 34 ns in a compatible BC
    sgCuts.SetFITAmpLimits({-1,   // Don't use the FV0A for selection
                            150,  // Require FT0A amplitude to be below 150 in all compatible BCs to classify as gap
                            50,   // Require FT0C amplitude to be below 50 in all compatible BCs to classify as gap
                            -1,   // Don't use the FDDA for selection
                            -1}); // Don't use the FDDC for selection
  }

  void process(MyEvents::iterator const& collision, // Process collision by collision
               MyBCs const& bcs,                    // We will check a range of bunch crossings
               o2::aod::FT0s const&,                // Must subscribe to the FIT tables for the SGSelector to access them
               o2::aod::FDDs const&,
               o2::aod::FV0As const&,
               o2::aod::Zdcs const&) // Want to plot ZDC energies, so we need to subscribe to the ZDC table
  {
    // Find the bunch crossing assigned to this collision
    auto bc = collision.foundBC_as<MyBCs>();
    // Find the range of bunch crossings compatible with this collision
    auto bcRange = udhelpers::compatibleBCs(collision, sgCuts.NDtcoll(), bcs, sgCuts.minNBCs());
    // Determine whether this event is single gap (A or C), double gap, or no gap
    auto selectorResult = sgSelector.IsSelected(sgCuts, collision, bcRange, bc);
    auto newbc = *(selectorResult.bc);

    // --- Process the event here: Apply cuts, save to derived tables, fill histograms... ---

    // Plot the outcome of the gap selection algorithm
    hSelectionResult->Fill(selectorResult.value);
    // Plot ZDC energies
    if (newbc.has_zdc()) {
      auto zdc = newbc.zdc();
      hZNAEnergy->Fill(zdc.energyCommonZNA());
      hZNCEnergy->Fill(zdc.energyCommonZNC());
    } else {
      hZNAEnergy->Fill(-999);
      hZNCEnergy->Fill(-999);
    }

    // Apply a selection -- as an example, keep only events classified as single gap on C side
    if (selectorResult.value != o2::aod::sgselector::SingleGapC) {
      return;
    }
    hSelectionResultAfterCut->Fill(selectorResult.value);
    if (newbc.has_zdc()) {
      auto zdc = newbc.zdc();
      hZNAEnergyAfterCut->Fill(zdc.energyCommonZNA());
      hZNCEnergyAfterCut->Fill(zdc.energyCommonZNC());
    } else {
      hZNAEnergyAfterCut->Fill(-999);
      hZNCEnergyAfterCut->Fill(-999);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial08>(cfgc, TaskName{"udtutorial08"})};
}
