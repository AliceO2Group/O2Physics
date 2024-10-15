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
#include "iostream"
#include <TString.h>
#include <TTree.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"
#include "TLorentzVector.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"

using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// \brief This is an example to fill outputs in a tree, useful for post processing analysis
/// \author Amrit Gautam
/// \author Anisa Khatun
/// \date 10.10.2024

namespace o2::aod
{
namespace tree
{
// track tables
DECLARE_SOA_COLUMN(SIGMAPI, sigmapi, std::vector<float>);
DECLARE_SOA_COLUMN(PT, Pt, float);
DECLARE_SOA_COLUMN(RAP, rap, float);
DECLARE_SOA_COLUMN(PHI, Phi, float);
DECLARE_SOA_COLUMN(TOTSIGN, totsign, int);
DECLARE_SOA_COLUMN(MASS, mass, float);
DECLARE_SOA_COLUMN(NPVTRACK, npvtrack, int);
DECLARE_SOA_COLUMN(PTS, Pts, std::vector<float>);
DECLARE_SOA_COLUMN(ETAS, etas, std::vector<float>);
DECLARE_SOA_COLUMN(PHIS, Phis, std::vector<float>);
DECLARE_SOA_COLUMN(SIGNS, Signs, std::vector<float>);
} // namespace tree

DECLARE_SOA_TABLE(TREE, "AOD", "Tree",
                  tree::PT,
                  tree::RAP,
                  tree::PHI,
                  tree::MASS,
                  tree::TOTSIGN,
                  tree::NPVTRACK,
                  tree::SIGMAPI,
                  tree::PTS,
                  tree::ETAS,
                  tree::PHIS,
                  tree::SIGNS);

} // namespace o2::aod

struct UDTutorial07 {
  SGSelector sgSelector;
  Produces<aod::TREE> tree;

  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 200., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};

  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  Configurable<float> TPC_cluster{"TPC_cluster", 50, "No.of TPC cluster"};

  // Kinmatic cuts
  Configurable<float> PID_cut{"PID_cut", 5, "TPC PID"};
  Configurable<float> Rap_cut{"Rap_cut", 0.9, "Track rapidity"};
  Configurable<float> Mass_Max{"Mass_Max", 10, "Invariant Mass range high"};
  Configurable<float> Mass_Min{"Mass_Min", 0, "Invariant Mass range low"};
  Configurable<float> Pt_coherent{"Pt_coherent", 0.15, "Coherent selection"};

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  //-----------------------------------------------------------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {

    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});

    // Fill counter to see effect of each selection criteria
    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});
    TString SelectionCuts[7] = {"NoSelection", "gapside", "goodtracks", "truegap", "2collcontrib", "2goodtrk"};
    for (int i = 0; i < 7; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    // Fill histograms to cross-check the selected numbers
    registry.add("hTracks", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
    registry.add("hTracksPions", "N_{tracks}", kTH1F, {{100, -0.5, 99.5}});
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>;
  //__________________________________________________________________________
  // Main process
  void process(UDCollisionsFull::iterator const& collision, udtracksfull const& tracks)
  {
    registry.fill(HIST("hSelectionCounter"), 0);

    // Accessing gap sides
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;

    registry.fill(HIST("hSelectionCounter"), 1);

    // Accessing FIT information for further exclusivity and/or inclusivity
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);

    // Intiating track parameters to select good tracks, values to be optimized in the configurables, parameters will be taken from SGtrackselector.h task included in the header
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);

    // Gap side to be selected in the configuables
    gapSide = truegapSide;

    registry.fill(HIST("hSelectionCounter"), 2);
    //_____________________________________
    // Create LorentzVector to store track information
    std::vector<TLorentzVector> allTracks;
    std::vector<TLorentzVector> onlyPionTracks;
    std::vector<float> onlyPionSigma;
    std::vector<decltype(tracks.begin())> rawPionTracks;
    std::vector<float> trackpt;
    std::vector<float> tracketa;
    std::vector<float> trackphi;
    std::vector<float> tracksign;
    std::vector<float> pitpcpid;

    TLorentzVector p;

    if (gapSide == gap_Side) {

      registry.fill(HIST("hSelectionCounter"), 3);

      for (auto t : tracks) {

        // Apply good track selection criteria
        if (!trackselector(t, parameters))
          continue;

        // Creating Lorenz vector to store raw tracks and piontracks
        TLorentzVector a;
        a.SetXYZM(t.px(), t.py(), t.pz(), o2::constants::physics::MassPionCharged);
        allTracks.push_back(a);

        // Apply TPC pion sigma
        auto nSigmaPi = t.tpcNSigmaPi();
        if (fabs(nSigmaPi) < PID_cut) {
          onlyPionTracks.push_back(a);
          onlyPionSigma.push_back(nSigmaPi);
          rawPionTracks.push_back(t);
        }
      }
      registry.fill(HIST("hTracksPions"), onlyPionTracks.size());

      //_____________________________________
      // Selecting collisions with Two PV contributors
      if (collision.numContrib() == 2) {

        registry.fill(HIST("hSelectionCounter"), 4);

        // Selecting only Two good tracks
        if ((rawPionTracks.size() == 2) && (allTracks.size() == 2)) {

          registry.fill(HIST("hSelectionCounter"), 4);

          // Creating rhos
          for (auto pion : onlyPionTracks) {
            p += pion;
          }

          for (auto rtrk : rawPionTracks) {
            TLorentzVector itrk;
            itrk.SetXYZM(rtrk.px(), rtrk.py(), rtrk.pz(), o2::constants::physics::MassPionCharged);
            trackpt.push_back(itrk.Pt());
            tracketa.push_back(itrk.Eta());
            trackphi.push_back(itrk.Phi());
            tracksign.push_back(rtrk.sign());
            pitpcpid.push_back(rtrk.tpcNSigmaPi());
          }

          int sign = 0;
          for (auto rawPion : rawPionTracks) {
            sign += rawPion.sign();
          }
          // Filling tree, make to be consistent with the declared tables
          tree(p.Pt(), p.Y(), p.Phi(), p.M(), sign, collision.numContrib(), pitpcpid, trackpt, tracketa, trackphi, tracksign);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial07>(cfgc, TaskName{"udtutorial07"})};
}
