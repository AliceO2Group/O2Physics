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
/// \author Johanna Lömker <johanna.lomker@cern.ch>
//  \since Dec 2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// struct for jetfinder validation on run2 ESD's and run3 data
struct jetTrackCollisionQa {
  // Track filter configs
  Configurable<float> ptLow{"ptLow", 0.15f, "lowest pt"};
  Configurable<float> ptUp{"ptUp", 10e10f, "highest pt"};
  Configurable<float> etalow{"etaLow", -0.9f, "lowest eta"};
  Configurable<float> etaup{"etaUp", 0.9f, "highest eta"};
  Configurable<bool> evSel{"evSel", false, "to use event selection 7 for run2 (or 8 for run3)"};

  HistogramRegistry mHistManager{"JetCollisionQAHistograms"};
  Configurable<int> nBins{"nBins", 200, "N bins in histos"}; // keep nBins for vertex and special 2D's
  // change the binning for pT in config file depending on AliPhysics status
  ConfigurableAxis BinsPhi{"BinsPhi", {200, -3.2, 6.4}, "Binning of the phi axis"};
  ConfigurableAxis BinsEta{"BinsEta", {200, -0.9, 0.9}, "Binning of the eta axis"};
  ConfigurableAxis BinsPt{"BinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};

  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int trackSelection = -1;

  void init(InitContext const&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // histograms
    const AxisSpec vtxZAxis{nBins, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec phiAxis{BinsPhi, "#phi "};
    const AxisSpec etaAxis{BinsEta, "#eta "};
    const AxisSpec ptAxis{BinsPt, "#it{p}_{T} (GeV/#it{c})"};

    // set trackselections
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    // histograms
    // 1)Jetvalidation on data
    mHistManager.add("collisionVtxZ", "selected collsion VtxZ ", HistType::kTH1D, {vtxZAxis});
    mHistManager.add("controlCollisionVtxZ", "control collsion VtxZ ", HistType::kTH1D, {vtxZAxis});
    // process jet qa
    mHistManager.add("jetPt", "inclusive jetPt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("jetPhi", "inclusive jet #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("jetEta", "inclusive jet #eta ", HistType::kTH1F, {etaAxis});
    // process jet constituent qa - constituents as tracks
    mHistManager.add("jetConstTrackPt", "inclusive jet constituent Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("jetConstTrackPhi", "inclusive jet constituent #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("jetConstTrackEta", "inclusive jet constituent #eta ", HistType::kTH1F, {etaAxis});
    // cross check the cuts from Run2Hybrid selection
    mHistManager.add("selectedTrackPt", "hybrid track Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("selectedTrackPhi", "hybrid track #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("selectedTrackEta", "hybrid track #eta ", HistType::kTH1F, {etaAxis});

    // leading jets per collision
    mHistManager.add("leadJetPt", "track Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("leadJetPhi", "track constituent #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("leadJetEta", "track constituent #eta ", HistType::kTH1F, {etaAxis});
    // leading constituents per jet in collision
    mHistManager.add("leadJetConstPt", "leading jet constituent Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("leadJetConstPhi", "leading jet constituent #phi", HistType::kTH1F, {phiAxis});
    mHistManager.add("leadJetConstEta", "leading jet constituent #eta", HistType::kTH1F, {etaAxis});
    // leading selected tracks per collision
    mHistManager.add("leadTrackPt", "leading selected track Pt", HistType::kTH1F, {ptAxis});
    mHistManager.add("leadTrackPhi", "leading selected track #phi", HistType::kTH1F, {phiAxis});
    mHistManager.add("leadTrackEta", "leading selected track #eta", HistType::kTH1F, {etaAxis});

    // Some more Track QA histos to test the implemented track cuts
    mHistManager.add("TPC/nClsCrossedRows", "nCls crossed rows TPC", HistType::kTH1F, {{165, -0.5, 164.5}});
    mHistManager.add("TPC/chi2PerCluster", "chi2 per cluster TPC", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TPC/crossedRowsOverFindableCls", "ratio of crossed rows over findable cluster TPC", HistType::kTH1F, {{120, 0.0, 1.2}});
    mHistManager.add("TPC/NClsFindable", "number of findable TPC clusters;# findable clusters TPC", HistType::kTH1F, {{165, -0.5, 164.5}});
    mHistManager.add("TPC/NClsFound", "number of found TPC clusters;# clusters TPC", HistType::kTH1F, {{165, -0.5, 164.5}});
    mHistManager.add("TPC/NClsShared", "number of shared TPC clusters;# shared clusters TPC", HistType::kTH1F, {{165, -0.5, 164.5}});
    mHistManager.add("TPC/FractionSharedCls", "fraction of shared TPC clusters;# fraction shared clusters TPC", HistType::kTH1F, {{100, 0., 1.}});

    mHistManager.add("ITS/chi2PerCluster", "chi2 per ITS cluster", HistType::kTH1F, {{100, 0, 40}});
    mHistManager.add("ITS/itsHits", "hitmap ITS;#it{p}_{T} [GeV/c];layer ITS", HistType::kTH1F, {{7, -0.5, 6.5}});
    mHistManager.add("ITS/itsNCls", "number of found ITS clusters;# clusters ITS", HistType::kTH1F, {{8, -0.5, 7.5}});
  }

  template <typename validationTracks>
  void fillTrackQA(validationTracks const& track)
  {
    mHistManager.fill(HIST("selectedTrackPt"), track.pt());
    mHistManager.fill(HIST("selectedTrackPhi"), track.phi());
    mHistManager.fill(HIST("selectedTrackEta"), track.eta());
    // START OF FILLING TRACK CUTS HISTOS
    mHistManager.fill(HIST("TPC/chi2PerCluster"), track.tpcChi2NCl());
    mHistManager.fill(HIST("TPC/nClsCrossedRows"), track.tpcNClsCrossedRows());
    mHistManager.fill(HIST("TPC/crossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
    mHistManager.fill(HIST("TPC/NClsFindable"), track.tpcNClsFindable());
    mHistManager.fill(HIST("TPC/NClsFound"), track.tpcNClsFound());
    mHistManager.fill(HIST("TPC/NClsShared"), track.tpcNClsShared());
    mHistManager.fill(HIST("TPC/FractionSharedCls"), track.tpcFractionSharedCls());

    mHistManager.fill(HIST("ITS/chi2PerCluster"), track.itsChi2NCl());
    mHistManager.fill(HIST("ITS/itsNCls"), track.itsNCls());
    for (unsigned int i = 0; i < 7; i++) {
      if (track.itsClusterMap() & (1 << i)) {
        mHistManager.fill(HIST("ITS/itsHits"), i);
      }
    }
  } // end of fillTrackQA template

  void fillLeadingTrackQA(double leadingTrackPt, double leadingTrackPhi, double leadingTrackEta)
  {
    mHistManager.fill(HIST("leadTrackPt"), leadingTrackPt);
    mHistManager.fill(HIST("leadTrackPhi"), leadingTrackPhi);
    mHistManager.fill(HIST("leadTrackEta"), leadingTrackEta);
  } // end of fillLeadingTrackQA template

  template <typename dataJet>
  void fillJetQA(dataJet const& jet)
  {
    mHistManager.fill(HIST("jetPt"), jet.pt());
    mHistManager.fill(HIST("jetPhi"), jet.phi());
    mHistManager.fill(HIST("jetEta"), jet.eta());
  } // end of fillJetQA template

  void fillLeadingJetQA(double leadingJetPt, double leadingJetPhi, double leadingJetEta)
  {
    mHistManager.fill(HIST("leadJetPt"), leadingJetPt);
    mHistManager.fill(HIST("leadJetPhi"), leadingJetPhi);
    mHistManager.fill(HIST("leadJetEta"), leadingJetEta);
  } // end of fillLeadingJetQA template

  template <typename dataJetConstituent>
  void fillJetConstituentQA(dataJetConstituent const& jct)
  {
    mHistManager.fill(HIST("jetConstTrackPt"), jct.pt());
    mHistManager.fill(HIST("jetConstTrackPhi"), jct.phi());
    mHistManager.fill(HIST("jetConstTrackEta"), jct.eta());
  } // end of mcDetJetConstituent template

  void fillLeadingJetConstQA(double leadingConstTrackPt, double leadingConstTrackPhi, double leadingConstTrackEta)
  {
    mHistManager.fill(HIST("leadJetConstPt"), leadingConstTrackPt);
    mHistManager.fill(HIST("leadJetConstPhi"), leadingConstTrackPhi);
    mHistManager.fill(HIST("leadJetConstEta"), leadingConstTrackEta);
  } // end of fillLeadingJetConstQA template

  Filter etafilter = (aod::jtrack::eta <= etaup) && (aod::jtrack::eta >= etalow);
  Filter ptfilter = (aod::jtrack::pt <= ptUp) && (aod::jtrack::pt >= ptLow);
  using Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using TracksJE = soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackPIs>>;

  void processESD(aod::JetCollision const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, TracksJE const& tracks, Tracks const&)
  {
    mHistManager.fill(HIST("controlCollisionVtxZ"), collision.posZ());
    if (evSel == true) {
      if (!jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel7")) || fabs(collision.posZ()) > 10) {
        return;
      }
    } else {
      if (fabs(collision.posZ()) > 10) {
        return;
      }
    }

    mHistManager.fill(HIST("collisionVtxZ"), collision.posZ());

    double leadingTrackPt = -1;
    double leadingTrackPhi = -1;
    double leadingTrackEta = -1;
    // qa histograms for selected tracks in collision
    for (const auto& t : tracks) {
      if (t.collisionId() == collision.globalIndex() && jetderiveddatautilities::selectTrack(t, trackSelection)) {
        auto track = t.track_as<Tracks>();
        fillTrackQA(track);
        if (track.pt() > leadingTrackPt) {
          leadingTrackPt = track.pt();
          leadingTrackPhi = track.phi();
          leadingTrackEta = track.eta();
        }
      }
    } // end of tracks loop
    // fill leading track
    fillLeadingTrackQA(leadingTrackPt, leadingTrackPhi, leadingTrackEta);

    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    // jet QA hists per jet in this collision
    for (const auto& j : jets) {
      if (j.collisionId() != collision.globalIndex()) {
        continue;
      }
      fillJetQA(j);
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
        fillJetConstituentQA(jct);
        if (jct.pt() > leadingConstTrackPt) {
          leadingConstTrackPt = jct.pt();
          leadingConstTrackPhi = jct.phi();
          leadingConstTrackEta = jct.eta();
        }
      } // end of jet constituent loop
      // fill leading jet constituent qa
      fillLeadingJetConstQA(leadingConstTrackPt, leadingConstTrackPhi, leadingConstTrackEta);
    } // end of jet loop
    // fill leading jet qa
    fillLeadingJetQA(leadingJetPt, leadingJetPhi, leadingJetEta);
  } // end process
  PROCESS_SWITCH(jetTrackCollisionQa, processESD, "validate jet-finder output on run2 ESD", true);

  // process for run3 AOD's
  void processRun3AOD(aod::JetCollision const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, TracksJE const& tracks, Tracks const&)
  {
    if (evSel == true) {
      if (!jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8")) || fabs(collision.posZ()) > 10) {
        return;
      }
    } else {
      if (fabs(collision.posZ()) > 10) {
        return;
      }
    }
    mHistManager.fill(HIST("collisionVtxZ"), collision.posZ());
    double leadingTrackPt = -1;
    double leadingTrackPhi = -1;
    double leadingTrackEta = -1;
    // qa histograms for selected tracks in collision
    for (const auto& t : tracks) {
      if (t.collisionId() == collision.globalIndex() && jetderiveddatautilities::selectTrack(t, trackSelection)) {
        auto track = t.track_as<Tracks>();
        fillTrackQA(track);
        if (track.pt() > leadingTrackPt) {
          leadingTrackPt = track.pt();
          leadingTrackPhi = track.phi();
          leadingTrackEta = track.eta();
        }
      }
    } // end of tracks loop
    // fill leading track
    fillLeadingTrackQA(leadingTrackPt, leadingTrackPhi, leadingTrackEta);

    double leadingJetPt = -1;
    double leadingJetPhi = -1;
    double leadingJetEta = -1;
    // jet QA hists per jet in this collision
    for (const auto& j : jets) {
      if (j.collisionId() != collision.globalIndex()) {
        continue;
      }
      fillJetQA(j);
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
        fillJetConstituentQA(jct);
        if (jct.pt() > leadingConstTrackPt) {
          leadingConstTrackPt = jct.pt();
          leadingConstTrackPhi = jct.phi();
          leadingConstTrackEta = jct.eta();
        }
      } // end of jet constituent loop
      // fill leading jet constituent qa
      fillLeadingJetConstQA(leadingConstTrackPt, leadingConstTrackPhi, leadingConstTrackEta);
    } // end of jet loop
    // fill leading jet qa
    fillLeadingJetQA(leadingJetPt, leadingJetPhi, leadingJetEta);
  }
  PROCESS_SWITCH(jetTrackCollisionQa, processRun3AOD, "validate jet-finder output on run3 AOD", false);

  // dummy process to run jetfinder validation code on ESD, but MC validation for run3 on hyperloop
  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(jetTrackCollisionQa, processDummy, "Dummy process function turned on by default", false);
};

// MC validation for run2 and run3 on AO2D's
struct mcJetTrackCollisionQa {
  // Track filter configs
  Configurable<float> ptLow{"ptLow", 0.15f, "lowest pt"};
  Configurable<float> ptUp{"ptUp", 10e15f, "highest pt"};
  Configurable<float> etalow{"etaLow", -0.9f, "lowest eta"};
  Configurable<float> etaup{"etaUp", 0.9f, "highest eta"};

  HistogramRegistry mHistManager{"JetCollisionQAHistograms"};
  Configurable<int> nBins{"nBins", 200, "N bins in histos"}; // keep nBins for vertex and special 2D's

  ConfigurableAxis BinsPhi{"BinsPhi", {200, -3.2, 6.4}, "Binning of the phi axis"};
  ConfigurableAxis BinsEta{"BinsEta", {200, -0.9, 0.9}, "Binning of the eta axis"};
  ConfigurableAxis BinsPt{"BinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};

  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int trackSelection;

  void init(InitContext const&)
  {
    // set trackselection
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    // histograms
    const AxisSpec vtxZAxis{nBins, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec ptAxis{BinsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec phiAxis{BinsPhi, "#phi "};
    const AxisSpec etaAxis{BinsEta, "#eta "};

    // 2)Jetvalidation on MC: generator = particle = mcTruth and reconstruction = detector (=tracks)
    mHistManager.add("collisionVtxZ", "Control collsion VtxZ ", HistType::kTH1F, {vtxZAxis});
    mHistManager.add("genMCcollisionVtxZ", "MC control gen.collsion VtxZ ", HistType::kTH1F, {vtxZAxis});
    mHistManager.add("recMCcollisionVtxZ", "MC control rec.collsion VtxZ ", HistType::kTH1F, {vtxZAxis});
    // 2D for reco vs. truth level
    mHistManager.add("collMatchPosZ", "MC reco vs truth; MC truth posZ (cm); MC reco posZ (cm)", {HistType::kTH2F, {{nBins, -15, 15}, {nBins, -15, 15}}});
    // 2D for 'relative resolution figure'
    mHistManager.add("collResolutionPt", "Collision reso #Delta posZ = (MC reco - MC truth)/ MC truth; MC truth posZ (cm); #Delta posZ (cm)", {HistType::kTH2F, {{nBins, -15, 15}, {nBins, -5, 5}}});

    // process jet qa
    mHistManager.add("genMCjetPt", "MC inclusive gen jetPt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("genMCjetPhi", "MC inclusive gen jet #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("genMCjetEta", "MC inclusive gen jet #eta ", HistType::kTH1F, {etaAxis});
    mHistManager.add("recMCjetPt", "MC inclusive rec jetPt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("recMCjetPhi", "MC inclusive rec jet #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("recMCjetEta", "MC inclusive rec jet #eta ", HistType::kTH1F, {etaAxis});
    // process mc matching from particle to det
    mHistManager.add("genRecMCjetPt", "MC rec to part jetPt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("genRecMCjetPhi", "MC reco to part jet #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("genRecMCjetEta", "MC rec to part jet #eta ", HistType::kTH1F, {etaAxis});
    // process jet constituent qa - constituents as tracks
    mHistManager.add("genMCjetConstTrackPt", "MC inclusive part jet constituent Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("genMCjetConstTrackPhi", "MC inclusive part jet constituent #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("genMCjetConstTrackEta", "MC inclusive part jet constituent #eta ", HistType::kTH1F, {etaAxis});
    mHistManager.add("recMCjetConstTrackPt", "MC inclusive reco jet constituent Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("recMCjetConstTrackPhi", "MC inclusive reco jet constituent #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("recMCjetConstTrackEta", "MC inclusive reco jet constituent #eta ", HistType::kTH1F, {etaAxis});
    // process mc matching from partice to detector - needs matching from nime / aimeric has something for it
    mHistManager.add("genRecMCjetConstTrackPt", "MC rec to part jet constituent Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("genRecMCjetConstTrackPhi", "MC rec to part inclusive part jet constituent #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("genRecMCjetConstTrackEta", "MC rec to part part jet constituent #eta ", HistType::kTH1F, {etaAxis});

    // cross check the cuts from Run2Hybrid selection
    mHistManager.add("genMCselectedTrackPt", "MC track Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("genMCselectedTrackPhi", "MC track #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("genMCselectedTrackEta", "MC track #eta ", HistType::kTH1F, {etaAxis});
    mHistManager.add("recMCselectedTrackPt", "reconstructed MC track Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("recMCselectedTrackPhi", "reconstructed MC track #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("recMCselectedTrackEta", "reconstructed MC track #eta ", HistType::kTH1F, {etaAxis});
    // tracks from mc data not mc associated
    mHistManager.add("selectedTrackPt", "selected collission tracks Pt ", HistType::kTH1F, {ptAxis});
    mHistManager.add("selectedTrackPhi", "selected collission tracks #phi ", HistType::kTH1F, {phiAxis});
    mHistManager.add("selectedTrackEta", "selected collission tracks #eta ", HistType::kTH1F, {etaAxis});
    // 2D for reco vs. truth level - we want this for jets too, but first we need proper matching there !
    mHistManager.add("trackMatchPt", "MC reco vs truth; MC truth p_{T} (GeV/#it{c}); MC reco p_{T} (GeV/#it{c})", {HistType::kTH2F, {{nBins, 0, 20}, {nBins, 0, 20}}});
    mHistManager.add("trackMatchEta", "MC reco vs truth; MC truth #eta; MC reco  #eta", {HistType::kTH2F, {{nBins, -0.9, 0.9}, {nBins, -0.9, 0.9}}});
    mHistManager.add("trackMatchPhi", "MC reco vs truth; MC truth #phi; MC reco #phi", {HistType::kTH2F, {{nBins, 0, 6.32}, {nBins, 0, 6.32}}});
    // 2D for 'relative resolution figure'
    mHistManager.add("trackResolutionPt", "Track reso #Delta p_{T} = (MC reco - MC truth)/ MC truth; MC truth p_{T} (GeV/#it{c}); #Delta p_{T}", {HistType::kTH2F, {{nBins, 0, 20}, {nBins, -3, 3}}});
    mHistManager.add("trackResolutionEta", "Track reso #Delta #eta = (MC reco - MC truth)/ MC truth; MC truth #Delta #eta", {HistType::kTH2F, {{nBins, -0.9, 0.9}, {nBins, -1, 1}}});
    mHistManager.add("trackResolutionPhi", "Track reso #Delta #phi = (MC reco - MC truth)/ MC truth; MC truth #phi (GeV/#it{c}); #Delta #phi", {HistType::kTH2F, {{nBins, 0, 6.32}, {nBins, -5, 5}}});
  }

  // fill collision qa histograms
  template <class mcCollisions>
  void fillMcCollisionHistos(mcCollisions const& collision)
  {
    mHistManager.fill(HIST("genMCcollisionVtxZ"), collision.mcCollision().posZ());
    mHistManager.fill(HIST("recMCcollisionVtxZ"), collision.posZ());
    mHistManager.fill(HIST("collMatchPosZ"), collision.mcCollision().posZ(), collision.posZ());
    // 2D for 'relative resolution figure'
    mHistManager.fill(HIST("collResolutionPt"), collision.mcCollision().posZ(), (collision.posZ() - collision.mcCollision().posZ()) / collision.mcCollision().posZ());
  } // end of collision template

  // fill qa histograms for selected tracks in collision
  template <class ValidationTracks, typename coll>
  void fillMcTrackHistos(ValidationTracks const& mct, coll collision, bool mc) // could give collision as argument for additional association
  {
    for (const auto& track : mct) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection) || !(track.collisionId() == collision.globalIndex())) {
        continue;
      }
      if (mc == true) {
        if (track.has_mcParticle()) {
          auto mcParticle = track.mcParticle();
          mHistManager.fill(HIST("genMCselectedTrackPt"), mcParticle.pt());
          mHistManager.fill(HIST("genMCselectedTrackPhi"), mcParticle.phi());
          mHistManager.fill(HIST("genMCselectedTrackEta"), mcParticle.eta());
          // 2D reco vs. gen
          mHistManager.fill(HIST("trackMatchPt"), mcParticle.pt(), track.pt());
          mHistManager.fill(HIST("trackMatchEta"), mcParticle.eta(), track.eta());
          mHistManager.fill(HIST("trackMatchPhi"), mcParticle.phi(), track.phi());
          // 2D for 'relative resolution figure'
          mHistManager.fill(HIST("trackResolutionPt"), mcParticle.pt(), (track.pt() - mcParticle.pt()) / mcParticle.pt());
          mHistManager.fill(HIST("trackResolutionEta"), mcParticle.eta(), (track.eta() - mcParticle.eta()) / mcParticle.eta());
          mHistManager.fill(HIST("trackResolutionPhi"), mcParticle.phi(), (track.phi() - mcParticle.phi()) / mcParticle.phi());
        } // end of track.has_mcParticle
        // fill histos for all reconstructed particles from mc associated collision (includes mcParticles)
        mHistManager.fill(HIST("recMCselectedTrackPt"), track.pt());
        mHistManager.fill(HIST("recMCselectedTrackPhi"), track.phi());
        mHistManager.fill(HIST("recMCselectedTrackEta"), track.eta());
      } // end of if mc
      if (mc == false) {
        mHistManager.fill(HIST("selectedTrackPt"), track.pt());
        mHistManager.fill(HIST("selectedTrackPhi"), track.phi());
        mHistManager.fill(HIST("selectedTrackEta"), track.eta());
      }
    } // end of tracks loop
  } // end of mcTrack template

  template <typename detectorJet>
  void fillMcDetJets(detectorJet const& mcdJet)
  {
    mHistManager.fill(HIST("recMCjetPt"), mcdJet.pt());
    mHistManager.fill(HIST("recMCjetPhi"), mcdJet.phi());
    mHistManager.fill(HIST("recMCjetEta"), mcdJet.eta());
  } // end of mcDetJet template

  template <typename particleJet>
  void fillMcPartJets(particleJet const& mcpJet)
  {
    mHistManager.fill(HIST("genMCjetPt"), mcpJet.pt());
    mHistManager.fill(HIST("genMCjetPhi"), mcpJet.phi());
    mHistManager.fill(HIST("genMCjetEta"), mcpJet.eta());
  } // end of mcPartJet template

  template <typename detectorJetConstituent>
  void fillMcDetJetConstituents(detectorJetConstituent const& mcdJConst)
  {
    mHistManager.fill(HIST("recMCjetConstTrackPt"), mcdJConst.pt());
    mHistManager.fill(HIST("recMCjetConstTrackPhi"), mcdJConst.phi());
    mHistManager.fill(HIST("recMCjetConstTrackEta"), mcdJConst.eta());
  } // end of mcDetJetConstituent template

  template <typename particleJetConstituent>
  void fillMcPartJetConstituents(particleJetConstituent const& mcpJConst)
  {
    mHistManager.fill(HIST("genMCjetConstTrackPt"), mcpJConst.pt());
    mHistManager.fill(HIST("genMCjetConstTrackPhi"), mcpJConst.phi());
    mHistManager.fill(HIST("genMCjetConstTrackEta"), mcpJConst.eta());
  } // end of mcPartJetConstituent template

  Filter etafilter = (aod::jtrack::eta < etaup) && (aod::jtrack::eta > etalow);
  Filter ptfilter = (aod::jtrack::pt < ptUp) && (aod::jtrack::pt > ptLow);
  using MCTracksJE = soa::Filtered<aod::JetTracksMCD>;

  void processMcRun2(aod::JetCollisionsMCD::iterator const& collision,
                     soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& mcPartJets,
                     soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcDetJets,
                     aod::JetParticles const&, aod::JetMcCollisions const&,
                     MCTracksJE const& tracks)
  {
    if (fabs(collision.posZ()) > 10) {
      return;
    }
    mHistManager.fill(HIST("collisionVtxZ"), collision.posZ());
    if (collision.has_mcCollision()) {
      fillMcCollisionHistos(collision);
      fillMcTrackHistos(tracks, collision, true);
      for (const auto& genJet : mcPartJets) {
        if (genJet.mcCollisionId() == collision.globalIndex()) {
          fillMcPartJets(genJet);
          for (auto& mcParticle : genJet.tracks_as<aod::JetParticles>()) {
            fillMcPartJetConstituents(mcParticle);
          }
        }
      } // end of loop particle level jets
    } // end if has mc collision
    fillMcTrackHistos(tracks, collision, false);
    for (const auto& detJet : mcDetJets) {
      if (detJet.collisionId() == collision.globalIndex()) {
        fillMcDetJets(detJet);
        for (auto& detConst : detJet.tracks_as<MCTracksJE>()) {
          fillMcDetJetConstituents(detConst);
        }
      }
    } // end of loop detector level jets
  } // end processMcRun2
  PROCESS_SWITCH(mcJetTrackCollisionQa, processMcRun2, "validate jet-finder output on converted run2 mc AOD's", false);

  void processMcRun3(aod::JetCollisionsMCD::iterator const& collision,
                     soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& mcPartJets,
                     soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& mcDetJets,
                     aod::JetParticles const&, aod::JetMcCollisions const&,
                     MCTracksJE const& tracks)
  {
    if (fabs(collision.posZ()) > 10) {
      return;
    }
    mHistManager.fill(HIST("collisionVtxZ"), collision.posZ());
    if (collision.has_mcCollision()) {
      fillMcCollisionHistos(collision);
      fillMcTrackHistos(tracks, collision, true);
      for (const auto& genJet : mcPartJets) {
        if (genJet.mcCollisionId() == collision.globalIndex()) {
          fillMcPartJets(genJet);
          for (auto& mcParticle : genJet.tracks_as<aod::JetParticles>()) {
            fillMcPartJetConstituents(mcParticle);
          }
        }
      } // end of loop particle level jets
    } //   end of loop if mc collision
    fillMcTrackHistos(tracks, collision, false);
    for (const auto& detJet : mcDetJets) {
      if (detJet.collisionId() == collision.globalIndex()) {
        fillMcDetJets(detJet);
        for (auto& detConst : detJet.tracks_as<MCTracksJE>()) {
          fillMcDetJetConstituents(detConst);
        }
      }
    } // end of loop detector level jets
  } // end processMcRun3
  PROCESS_SWITCH(mcJetTrackCollisionQa, processMcRun3, "validate jet-finder output on run3 mc AOD's", false);

  // dummy process to run jetfinder validation code on AO2D's, but MC validation for run3 on hyperloop
  void processDummy(aod::JetMcCollisions const&)
  {
  }
  PROCESS_SWITCH(mcJetTrackCollisionQa, processDummy, "Dummy process function turned off by default", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jetTrackCollisionQa>(cfgc, TaskName{"jet-validation-track-collision-qa"}),
    adaptAnalysisTask<mcJetTrackCollisionQa>(cfgc, TaskName{"mc-jet-validation-track-collision-qa"})};
}
