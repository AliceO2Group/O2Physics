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
///
//  \author
//  Johanna Lömker
//  \since Dec 2022
// The goal would be to set up the jet, track and event/collision QA's for the validation framework
// Staring with the hybrid tracks

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/Jet.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksJE = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

struct jetTrackCollisionQa {

  HistogramRegistry mHistManager{"JetCollisionQAHistograms"};
  Configurable<int> nBins{"nBins", 200, "N bins in histos"};
  Configurable<int> nBinsPt{"nBinsPt", 200, "N bins in pT histos"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "N bins in Eta histos"};
  Configurable<int> nBinsPhi{"nBinsPhi", 200, "N bins in Phi histos"};

  TrackSelection validationTracks;

  void init(InitContext const&)
  {
    /*set trackselections - initial idea for now on hold
    //we should not touch the getJEGlobalTrackSelectionRun2 anymore, but here we are free to modify the defaults
    validationTracks = getJEGlobalTrackSelectionRun2();
    validationTracks.SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
    validationTracks.SetMinNCrossedRowsTPC(70);
    validationTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    validationTracks.SetMaxChi2PerClusterTPC(4.f);
    validationTracks.SetRequireTPCRefit(true);
    validationTracks.SetRequireITSRefit(true);
    validationTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
    validationTracks.SetMaxChi2PerClusterITS(36.f);
    validationTracks.SetPtRange(0.15f, 1e10f); // the 0.15 is the original hybrid cuts
    validationTracks.SetEtaRange(-0.9f, 0.9f);
    validationTracks.SetMaxDcaXY(2.4f);
    validationTracks.SetMaxDcaZ(3.2f);
    validationTracks.print();*/
    validationTracks = getGlobalTrackSelection(); // using global tracks instead

    // histograms
    // 1)Jetvalidation from data
    mHistManager.add("collisionVtxZ", "control collsion VtxZ ; z [cm]", HistType::kTH1F, {{nBins, -15, 15}});
    // process jet qa
    mHistManager.add("jetPt", "inclusive jetPt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("jetPhi", "inclusive jet #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, -3.2, 6.4}});
    mHistManager.add("jetEta", "inclusive jet #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
    // process jet constituent qa - constituents as tracks
    mHistManager.add("jetConstTrackPt", "inclusive jet constituent Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("jetConstTrackPhi", "inclusive jet constituent #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("jetConstTrackEta", "inclusive jet constituent #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
    // cross check the cuts from Run2Hybrid selection
    mHistManager.add("selectedTrackPt", "hybrid track Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("selectedTrackPhi", "hybrid track #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("selectedTrackEta", "hybrid track #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});

    // leading jets per collision
    mHistManager.add("leadJetPt", "track Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("leadJetPhi", "track constituent #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("leadJetEta", "track constituent #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
    // leading constituents per jet in collision
    mHistManager.add("leadJetConstPt", "leading jet constituent Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("leadJetConstPhi", "leading jet constituent #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("leadJetConstEta", "leading jet constituent #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
    // leading selected tracks per collision
    mHistManager.add("leadTrackPt", "leading selected track Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("leadTrackPhi", "leading selected track #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("leadTrackEta", "leading selected track #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});

    // 2)Jetvalidation from MC
    mHistManager.add("MCcollisionVtxZ", "MC control collsion VtxZ ; z [cm]", HistType::kTH1F, {{nBins, -15, 15}});
    // process jet qa
    mHistManager.add("MCjetPt", "MC inclusive jetPt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("MCjetPhi", "MC inclusive jet #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, -3.2, 6.4}});
    mHistManager.add("MCjetEta", "MC inclusive jet #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
    // process jet constituent qa - constituents as tracks
    mHistManager.add("MCjetConstTrackPt", "MC inclusive jet constituent Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("MCjetConstTrackPhi", "MC inclusive jet constituent #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("MCjetConstTrackEta", "MC inclusive jet constituent #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
    // cross check the cuts from Run2Hybrid selection
    mHistManager.add("MCselectedTrackPt", "MC hybrid track Pt ; p_{T} (GeV/#it{c})", HistType::kTH1F, {{nBinsPt, 0, 100}});
    mHistManager.add("MCselectedTrackPhi", "MC hybrid track #phi ; #phi ", HistType::kTH1F, {{nBinsPhi, 0, 6.4}});
    mHistManager.add("MCselectedTrackEta", "MC hybrid track #eta ; #eta ", HistType::kTH1F, {{nBinsEta, -0.9, 0.9}});
  }
  // add another process for MC studies !
  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::Jets, aod::JetConstituents> const& jets, TracksJE const& tracks)
  {
    mHistManager.fill(HIST("collisionVtxZ"), collision.posZ());

    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    // jet QA hists per jet in this collision
    for (const auto& j : jets) {
      mHistManager.fill(HIST("jetPt"), j.pt());
      mHistManager.fill(HIST("jetPhi"), j.phi());
      mHistManager.fill(HIST("jetEta"), j.eta());
      if (j.pt() > leadingJetPt) {
        leadingJetPt = j.pt();
        leadingJetPhi = j.phi();
        leadingJetEta = j.eta();
      }

      double leadingConstTrackPt = -1;
      double leadingConstTrackPhi = -1;
      double leadingConstTrackEta = -1;
      // access jet constituents as tracks
      for (auto& jct : j.tracks_as<TracksJE>()) {
        mHistManager.fill(HIST("jetConstTrackPt"), jct.pt());
        mHistManager.fill(HIST("jetConstTrackPhi"), jct.phi());
        mHistManager.fill(HIST("jetConstTrackEta"), jct.eta());
        if (jct.pt() > leadingConstTrackPt) {
          leadingConstTrackPt = jct.pt();
          leadingConstTrackPhi = jct.phi();
          leadingConstTrackEta = jct.eta();
        }
      } // end of jet constituent loop

      // fill leading jet constituent qa
      mHistManager.fill(HIST("leadJetConstPt"), leadingConstTrackPt);
      mHistManager.fill(HIST("leadJetConstPhi"), leadingConstTrackPhi);
      mHistManager.fill(HIST("leadJetConstEta"), leadingConstTrackEta);
    } // end of jet loop

    // fill leading jet qa
    mHistManager.fill(HIST("leadJetPt"), leadingJetPt);
    mHistManager.fill(HIST("leadJetPhi"), leadingJetPhi);
    mHistManager.fill(HIST("leadJetEta"), leadingJetEta);

    double leadingTrackPt = -1;
    double leadingTrackPhi = -1;
    double leadingTrackEta = -1;
    // qa histograms for selected tracks in collision
    for (const auto& t : tracks) {
      if (!validationTracks.IsSelected(t)) {
        continue;
      } // check if this is really a global track - maybe adding hists for rejected tracks ?
      mHistManager.fill(HIST("selectedTrackPt"), t.pt());
      mHistManager.fill(HIST("selectedTrackPhi"), t.phi());
      mHistManager.fill(HIST("selectedTrackEta"), t.eta());
      if (t.pt() > leadingTrackPt) {
        leadingTrackPt = t.pt();
        leadingTrackPhi = t.phi();
        leadingTrackEta = t.eta();
      }
    } // end of tracks loop

    // fill leading selected track qa
    mHistManager.fill(HIST("leadTrackPt"), leadingTrackPt);
    mHistManager.fill(HIST("leadTrackPhi"), leadingTrackPhi);
    mHistManager.fill(HIST("leadTrackEta"), leadingTrackEta);
  } // end process
  PROCESS_SWITCH(jetTrackCollisionQa, processData, "process jets from jet-finder output", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jetTrackCollisionQa>(cfgc, TaskName{"jet-validation-track-collision-qa"})};
}
