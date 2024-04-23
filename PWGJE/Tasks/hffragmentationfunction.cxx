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

// charm lambda/d0 ratio hadronization task
//
/// \author Christian Reckziegel <christian.reckziegel@cern.ch>
//
#include "TVector3.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// creating table for storing distance data
namespace o2::aod
{
namespace DistanceSpace
{
DECLARE_SOA_COLUMN(JetHfDist, jethfdist, float);
DECLARE_SOA_COLUMN(JetPt, jetpt, float);
DECLARE_SOA_COLUMN(JetEta, jeteta, float);
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);
DECLARE_SOA_COLUMN(HfPt, hfpt, float);
DECLARE_SOA_COLUMN(HfEta, hfeta, float);
DECLARE_SOA_COLUMN(HfPhi, hfphi, float);
DECLARE_SOA_COLUMN(HfMass, hfmass, float);
DECLARE_SOA_COLUMN(HfY, hfy, float);
} // namespace DistanceSpace
DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE",
                  DistanceSpace::JetHfDist,
                  DistanceSpace::JetPt,
                  DistanceSpace::JetEta,
                  DistanceSpace::JetPhi,
                  DistanceSpace::HfPt,
                  DistanceSpace::HfEta,
                  DistanceSpace::HfPhi,
                  DistanceSpace::HfMass,
                  DistanceSpace::HfY);
} // namespace o2::aod

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

struct HfFragmentationFunctionTask {
  // producing new table
  Produces<aod::JetDistanceTable> distJetTable;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // create histograms
    // D0 candidate histograms
    registry.add("h_jet_counter", ";# jets;", {HistType::kTH1F, {{2, 0., 1.}}});
    registry.add("h_d0_jet_projection", ";z^{D^{0},jet}_{||};dN/dz^{D^{0},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_jet_distance_vs_projection", ";#DeltaR_{D^{0},jet};z^{D^{0},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
    registry.add("h_d0_jet_distance", ";#DeltaR_{D^{0},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_jet_pt", ";p_{T,D^{0} jet};dN/dp_{T,D^{0} jet}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add("h_d0_jet_eta", ";#eta_{T,D^{0} jet};dN/d#eta_{D^{0} jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_d0_jet_phi", ";#phi_{T,D^{0} jet};dN/d#phi_{D^{0} jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add("h_d0_mass", ";m_{D^{0}} (GeV/c^{2});dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_eta", ";#eta_{D^{0}} (GeV/c^{2});dN/d#eta_{D^{0}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_d0_phi", ";#phi_{D^{0}} (GeV/c^{2});dN/d#phi_{D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
  }

  void processDummy(aod::TracksIU const&) {}
  PROCESS_SWITCH(HfFragmentationFunctionTask, processDummy, "Dummy process function turned on by default", false);

  void processDataChargedSubstructure(JetCollision const&,
                                      soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                                      JetTracks const&,
                                      CandidatesD0Data const&)
  {
    // CandidatesLcData const& lccands) {

    double axisDistance = 0;

    for (auto& jet : jets) {
      // fill jet counter histogram
      registry.fill(HIST("h_jet_counter"), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (auto& d0Candidate : jet.hfcandidates_as<CandidatesD0Data>()) { // for jet constituents use -> auto& jetConstituent : jet.tracks_as<JetTracks>()

        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());
        // calculating fraction of the jet momentum carried by the D0 along the direction of the jet axis
        double z_parallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        axisDistance = sqrt(pow(jet.eta() - d0Candidate.eta(), 2) + pow(jet.phi() - d0Candidate.phi(), 2));

        // filling histograms
        registry.fill(HIST("h_d0_jet_projection"), z_parallel);
        registry.fill(HIST("h_d0_jet_distance_vs_projection"), axisDistance, z_parallel);
        registry.fill(HIST("h_d0_jet_distance"), axisDistance);
        registry.fill(HIST("h_d0_jet_pt"), jet.pt());
        registry.fill(HIST("h_d0_jet_eta"), jet.eta());
        registry.fill(HIST("h_d0_jet_phi"), jet.phi());
        registry.fill(HIST("h_d0_mass"), d0Candidate.m());
        registry.fill(HIST("h_d0_eta"), d0Candidate.eta());
        registry.fill(HIST("h_d0_phi"), d0Candidate.phi());
        // filling table
        distJetTable(axisDistance, jet.pt(), jet.eta(), jet.phi(), d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.m(), d0Candidate.y());
        break; // get out of candidates' loop after first HF particle is found in jet
      }        // end of D0 candidates loop

    } // end of jets loop

  } // end of process function
  PROCESS_SWITCH(HfFragmentationFunctionTask, processDataChargedSubstructure, "charged HF jet substructure", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfFragmentationFunctionTask>(cfgc, TaskName{"jet-charm-hadronization"})};
}
