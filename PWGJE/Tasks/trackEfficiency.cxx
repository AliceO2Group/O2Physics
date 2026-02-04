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

/// \file   trackEfficiency.cxx
/// \author Aimeric Landou <aimeric.landou@cern.ch>
/// \brief task that creates the histograms necessary for computation of efficiency and purity functions in offline postprocess macros; also can make mcparticle and track QC histograms

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackEfficiency {
  Service<o2::framework::O2DatabasePDG> pdg;

  using JetParticlesWithOriginal = soa::Join<aod::JetParticles, aod::JMcParticlePIs>;

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections; other option: uniformTracks"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events"};

  // Tracking efficiency process function configurables:
  Configurable<bool> checkPrimaryPart{"checkPrimaryPart", true, "0: doesn't check mcparticle.isPhysicalPrimary() - 1: checks particle.isPhysicalPrimary()"};
  Configurable<bool> cutCentrality{"cutCentrality", false, ""};
  Configurable<bool> checkCentFT0M{"checkCentFT0M", false, "0: centFT0C as default, 1: use centFT0M estimator"};
  Configurable<bool> checkOccupancy{"checkOccupancy", false, "check occupancy only in general purpose Pb-Pb MC, default as false"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<float> trackEtaAcceptanceCountQA{"trackEtaAcceptanceCountQA", 0.9, "eta acceptance"}; // removed from actual cuts for now because all the histograms have an eta axis
  Configurable<float> centralityMin{"centralityMin", -999, ""};
  Configurable<float> centralityMax{"centralityMax", 999, ""};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackDcaZmax{"trackDcaZmax", 99, "additional cut on dcaZ to PV for tracks; uniformTracks in particular don't cut on this at all"};
  Configurable<int> nBinsLowPt{"nBinsLowPt", 200, "number of pt bins for low pt (below 10GeV) efficiency histograms"};

  // Track QA process function configurables:
  Configurable<float> trackQAEtaMin{"trackQAEtaMin", -0.9, "minimum eta acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAEtaMax{"trackQAEtaMax", 0.9, "maximum eta acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAPtMin{"trackQAPtMin", 0.15, "minimum pT acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAPtMax{"trackQAPtMax", 100.0, "maximum pT acceptance for tracks in the processTracks QA"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied for reconstructed tracks, not mc particles"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied for reconstructed tracks, not mc particles"};

  Configurable<std::vector<double>> centralityBinning{"centralityBinning", {0., 10., 50., 70., 100}, "binning of centrality histograms"};
  Configurable<int> intRateNBins{"intRateNBins", 50, "number of bins for interaction rate axis"};
  Configurable<float> intRateMax{"intRateMax", 50000.0, "maximum value of interaction rate axis"};
  Configurable<int> phiEffNBins{"phiEffNBins", 200, "number of bins for phi axis in efficiency plots"};
  Configurable<int> etaEffNBins{"etaEffNBins", 200, "number of bins for eta axis in efficiency plots"};

  Configurable<float> ptHatMin{"ptHatMin", 5, "min pT hat of collisions"};
  Configurable<float> ptHatMax{"ptHatMax", 300, "max pT hat of collisions"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMaxFractionMCD{"pTHatMaxFractionMCD", 999.0, "maximum fraction of hard scattering for reconstructed track acceptance in MC"};

  Configurable<bool> getPtHatFromHepMCXSection{"getPtHatFromHepMCXSection", true, "test configurable, configurable should be removed once well tested"};
  Configurable<bool> useTrueTrackWeight{"useTrueTrackWeight", true, "test configurable, should be set to 1 then config removed once well tested"};

  // systematics variation - Run 2 guidelines: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertainty
  TrackSelection customTrackSelection;
  Configurable<bool> useCustomTrackSelection{"useCustomTrackSelection", false, "whether to use the custom cuts (used for cut variation for tracking efficiency systematics)"};
  Configurable<int> effSystMinNCrossedRowsTPC{"effSystMinNCrossedRowsTPC", 70, "min number of crossed rows TPC"};
  Configurable<float> effSystMinNCrossedRowsOverFindableClustersTPC{"effSystMinNCrossedRowsOverFindableClustersTPC", 0.8, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> effSystMaxChi2PerClusterTPC{"effSystMaxChi2PerClusterTPC", 4.0, "max chi2 per cluster TPC"};
  Configurable<float> effSystMaxChi2PerClusterITS{"effSystMaxChi2PerClusterITS", 36.0, "max chi2 per cluster ITS"};
  // Configurable<float> effSystMaxDcaXY{"effSystMaxDcaXY", 0.0105 * 0.035 / pT^1.1 ????, "max DCA to vertex xy"}; not including this for now as it's a function with 3 parameters
  Configurable<float> effSystMaxDcaZ{"effSystMaxDcaZ", 2.0, "max DCA to vertex z"};
  Configurable<int> effSystMinNrequiredHits{"effSystMinNrequiredHits", 1, "minimum number of hits among the 3 innermost layers of the ITS"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  float simPtRef = 10.;

  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  template <typename TJetTrack>
  bool isAcceptedTrack(TJetTrack const& jetTrack)
  {
    if (!useCustomTrackSelection) {
      if (jetderiveddatautilities::selectTrack(jetTrack, trackSelection) && jetderiveddatautilities::selectTrackDcaZ(jetTrack, trackDcaZmax)) { // if track selection is uniformTrack, dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
        return true;
      }
    } else {
      const auto& aodTrack = jetTrack.template track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>>();
      if (customTrackSelection.IsSelected(aodTrack)) {
        return true;
      }
    }
    return false;
  }

  bool isChargedParticle(int code)
  {
    const float chargeUnit = 3.;
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= chargeUnit;
  }

  template <typename TCollision, typename TJetTracks>
  void fillTrackHistograms(TCollision const& collision, TJetTracks const& jetTracks, float weight = 1.0)
  {
    for (auto const& track : jetTracks) {
      if (!isAcceptedTrack(track)) {
        continue;
      }

      float pTHat = simPtRef / (std::pow(weight, 1.0 / pTHatExponent));
      if (track.pt() > pTHatMaxFractionMCD * pTHat) {
        continue;
      }

      float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();

      registry.fill(HIST("h2_centrality_track_pt"), centrality, track.pt(), weight);
      registry.fill(HIST("h2_centrality_track_eta"), centrality, track.eta(), weight);
      registry.fill(HIST("h2_centrality_track_phi"), centrality, track.phi(), weight);
      registry.fill(HIST("h2_centrality_track_energy"), centrality, track.energy(), weight);
      registry.fill(HIST("h2_track_pt_track_sigma1overpt"), track.pt(), track.sigma1Pt(), weight);
      registry.fill(HIST("h2_track_pt_track_sigmapt"), track.pt(), track.sigma1Pt() * track.pt(), weight);
      registry.fill(HIST("h2_track_pt_high_track_sigma1overpt"), track.pt(), track.sigma1Pt(), weight);
      registry.fill(HIST("h2_track_pt_high_track_sigmapt"), track.pt(), track.sigma1Pt() * track.pt(), weight);
      registry.fill(HIST("h3_intrate_centrality_track_pt"), collision.hadronicRate(), centrality, track.pt(), weight);
    }
  }

  template <typename TMCCollision, typename TCollisions, typename TParticles, typename TTracks>
  void fillParticlesHistograms(TMCCollision const& /*mcCollision*/, TCollisions const& collisions, TParticles const& mcparticles, TTracks tracks, float weight = 1.0)
  {
    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C(); mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    float centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();

    for (auto const& mcparticle : mcparticles) {
      registry.fill(HIST("h2_centrality_particle_pt"), centrality, mcparticle.pt(), weight);
      registry.fill(HIST("h2_centrality_particle_eta"), centrality, mcparticle.eta(), weight);
      registry.fill(HIST("h2_centrality_particle_phi"), centrality, mcparticle.phi(), weight);
      registry.fill(HIST("h2_centrality_particle_energy"), centrality, mcparticle.energy(), weight);
      registry.fill(HIST("h3_intrate_centrality_particle_pt"), collisions.begin().hadronicRate(), centrality, mcparticle.pt(), weight);
      auto partTracks = tracks.sliceBy(tracksPerJParticles, mcparticle.globalIndex());
      for (auto const& track : partTracks) {
        registry.fill(HIST("h2_particle_pt_track_pt_deltapt"), mcparticle.pt(), mcparticle.pt() - track.pt(), weight);
        registry.fill(HIST("h2_particle_pt_track_pt_deltaptoverparticlept"), mcparticle.pt(), (mcparticle.pt() - track.pt()) / mcparticle.pt(), weight);
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    if (!(acceptSplitCollisions == NonSplitOnly || acceptSplitCollisions == SplitOkCheckAnyAssocColl || acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly)) {
      LOGP(fatal, "Configurable acceptSplitCollisions has wrong input value; stopping workflow");
    }

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (useCustomTrackSelection) {
      // Custom track cuts
      LOGP(info, "Using custom track selection from values:");
      LOGP(info, "\tminNCrossedRowsTPC= %f", effSystMinNCrossedRowsTPC.value);
      LOGP(info, "\tminNCrossedRowsOverFindableClustersTPC= %f", effSystMinNCrossedRowsOverFindableClustersTPC.value);
      LOGP(info, "\tmaxChi2PerClusterTPC= %f", effSystMaxChi2PerClusterTPC.value);
      LOGP(info, "\tmaxChi2PerClusterITS= %f", effSystMaxChi2PerClusterITS.value);
      // LOGP(info, "\tmaxDcaXY= %f", effSystMaxDcaXY.value);
      LOGP(info, "\tmaxDcaZ= %f", effSystMaxDcaZ.value);
      LOGP(info, "\tRequireHitsInITSLayers= %i", effSystMinNrequiredHits.value);

      LOGP(info, "\trequireITS= true");
      LOGP(info, "\trequireTPC= true");

      LOGP(info, "Customizing track selection:");
      int dcaSetup = 0;                                                                                                               // default dca setup
      customTrackSelection = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, dcaSetup); // takes global tracks configuration, then some of the cuts are edited in the lines below
      customTrackSelection.SetEtaRange(-999, 999);
      customTrackSelection.SetPtRange(0, 1e10f);

      customTrackSelection.SetMinNCrossedRowsTPC(effSystMinNCrossedRowsTPC.value);
      customTrackSelection.SetMinNCrossedRowsOverFindableClustersTPC(effSystMinNCrossedRowsOverFindableClustersTPC.value);
      customTrackSelection.SetMaxChi2PerClusterTPC(effSystMaxChi2PerClusterTPC.value);
      customTrackSelection.SetMaxChi2PerClusterITS(effSystMaxChi2PerClusterITS.value);
      // customTrackSelection.SetMaxDcaXY(effSystMaxDcaXY.value);
      customTrackSelection.SetMaxDcaZ(effSystMaxDcaZ.value);
      customTrackSelection.SetRequireHitsInITSLayers(effSystMinNrequiredHits.value, {0, 1, 2}); // one hit in any SPD layer (#hits, {layer0, layer1,...})

      // customTrackSelection.SetRequireITSRefit(true); already set by default
      // customTrackSelection.SetRequireTPCRefit(true); already set by default
      // customTrackSelection.SetRequireGoldenChi2(requireGoldenChi2.value); already set by default

      customTrackSelection.print();
    } else {
      LOGP(info, "Using standard track selection: %s", trackSelections.value);
    }

    if (doprocessEFficiencyPurity || doprocessEFficiencyPurityWeighted) {

      registry.add("hMcCollCutsCounts", "McColl cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(2, "vertexZ");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(4, "splitColl");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(5, "recoCollEvtSel");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(6, "centralityCut");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(7, "ptHatCut");
      if (checkOccupancy) {
        registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(8, "occupancyCut");
      }

      registry.add("hMcPartCutsCounts", "McPart cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(1, "allPartsInSelMcColl");
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(2, "isCharged");
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(3, "isPrimary");
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(4, "etaAccept"); // not actually applied here but it will give an idea of what will be done in the post processing

      registry.add("hTrackCutsCounts", "Track cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(1, "allTracksInSelColl");
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(2, "trackSel");
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(3, "hasMcParticle");

      if (doprocessEFficiencyPurity) {
        registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(4, "mcPartIsPrimary");
        registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(5, "etaAcc"); // not actually applied here but it will give an idea of what will be done in the post processing
      }
      if (doprocessEFficiencyPurityWeighted) {
        registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(4, "ptHatMaxFraction");
        registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(5, "mcPartIsPrimary");
        registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(6, "etaAcc"); // not actually applied here but it will give an idea of what will be done in the post processing
      }

      AxisSpec ptAxisEff = {nBinsLowPt, 0., 10., "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec ptAxisHighEff = {18, 10., 100., "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec etaAxisEff = {etaEffNBins, -1.0, 1.0, "#eta"};
      AxisSpec phiAxisEff = {phiEffNBins, -1.0, 7., "#phi"};

      // ptAxisLow
      registry.add("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});

      registry.add("h3_track_pt_track_eta_track_phi_nonassociatedtrack", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_split_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_split_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});

      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisEff, etaAxisEff, phiAxisEff}});

      registry.add("h2_particle_pt_track_pt_residual_associatedtrack_primary", "(#it{p}_{T, mcpart} - #it{p}_{T, track}) / #it{p}_{T, mcpart}; #it{p}_{T, mcpart} (GeV/#it{c})", {HistType::kTH2F, {ptAxisEff, {200, -1., 1.}}});

      // ptAxisHigh
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});

      registry.add("h3_track_pt_high_track_eta_track_phi_nonassociatedtrack", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});

      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHighEff, etaAxisEff, phiAxisEff}});

      registry.add("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary", "(#it{p}_{T, mcpart} - #it{p}_{T, track}) / #it{p}_{T, mcpart}; #it{p}_{T, mcpart} (GeV/#it{c})", {HistType::kTH2F, {ptAxisHighEff, {200, -1., 1.}}});
    }

    if (doprocessTracksFromData || doprocessTracksFromMc || doprocessTracksFromMcWeighted) {
      AxisSpec centAxis = {centralityBinning, "centrality (%)"};
      AxisSpec intRateAxis = {intRateNBins, 0., intRateMax, "int. rate (kHz)"};
      registry.add("h2_centrality_track_pt", "centrality vs track pT; centrality; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {centAxis, {200, 0., 200.}}});
      registry.add("h2_centrality_track_eta", "centrality vs track #eta; centrality; #eta_{track}", {HistType::kTH2F, {centAxis, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_track_phi", "centrality vs track #varphi; centrality; #varphi_{track}", {HistType::kTH2F, {centAxis, {160, -1.0, 7.}}});
      registry.add("h2_centrality_track_energy", "centrality vs track energy; centrality; Energy GeV", {HistType::kTH2F, {centAxis, {100, 0.0, 100.0}}});
      registry.add("h2_track_pt_track_sigmapt", "#sigma(#it{p}_{T})/#it{p}_{T}; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 10.}, {100000, 0.0, 100.0}}});
      registry.add("h2_track_pt_high_track_sigmapt", "#sigma(#it{p}_{T})/#it{p}_{T}; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{90, 10., 100.}, {100000, 0.0, 100.0}}});
      registry.add("h2_track_pt_track_sigma1overpt", "#sigma(1/#it{p}_{T}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 10.}, {1000, 0.0, 10.0}}});
      registry.add("h2_track_pt_high_track_sigma1overpt", "#sigma(1/#it{p}_{T}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{90, 10., 100.}, {1000, 0.0, 10.0}}});
      registry.add("h3_intrate_centrality_track_pt", "interaction rate vs centrality vs track pT; int. rate; centrality; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH3F, {intRateAxis, centAxis, {200, 0., 200.}}});
    }

    if (doprocessParticles || doprocessParticlesWeighted) {
      AxisSpec centAxis = {centralityBinning, "centrality (%)"};
      AxisSpec intRateAxis = {intRateNBins, 0., intRateMax, "int. rate (kHz)"};
      registry.add("h2_centrality_particle_pt", "centrality vs particle pT; centrality; #it{p}_{T,part} (GeV/#it{c})", {HistType::kTH2F, {centAxis, {200, 0., 200.}}});
      registry.add("h2_centrality_particle_eta", "centrality vs particle #eta; centrality; #eta_{part}", {HistType::kTH2F, {centAxis, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_particle_phi", "centrality vs particle #varphi; centrality; #varphi_{part}", {HistType::kTH2F, {centAxis, {160, -1.0, 7.}}});
      registry.add("h2_centrality_particle_energy", "centrality vs particle energy; centrality; Energy GeV", {HistType::kTH2F, {centAxis, {100, 0.0, 100.0}}});
      registry.add("h3_intrate_centrality_particle_pt", "interaction rate vs centrality vs particle pT; int. rate; centrality; #it{p}_{T,part} (GeV/#it{c})", {HistType::kTH3F, {intRateAxis, centAxis, {200, 0., 200.}}});

      registry.add("h2_particle_pt_track_pt_deltapt", "track pt vs delta pT; pT; #it{p}_{T, part} - #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, -1., 1.}}});
      registry.add("h2_particle_pt_track_pt_deltaptoverparticlept", "track vs delta pT / MC pT ; pT; #frac{#it{p}_{T, part} - #it{p}_{T,track}}{#it{p}_{T,part}}", {HistType::kTH2F, {{200, 0., 200.}, {200, -1., 1.}}});
    }

    if (doprocessCollisionsFromData || doprocessCollisionsFromMc || doprocessCollisionsFromMcWeighted) {
      AxisSpec centAxis = {centralityBinning, "centrality (%)"};
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_collisions", "centrality vs collisions; centrality; collisions", {HistType::kTH2F, {centAxis, {4, 0.0, 4.0}}});
    }
    if (doprocessMcCollisions || doprocessMcCollisionsWeighted) {
      AxisSpec centAxis = {centralityBinning, "centrality (%)"};
      registry.add("h_mccollisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_mccollisions", "centrality vs mccollisions; centrality; collisions", {HistType::kTH2F, {centAxis, {4, 0.0, 4.0}}});
      registry.add("h2_mccollision_pthardfromweight_pthardfromhepmcxsection", "ptHard from weight vs ptHard from HepMCXSections; ptHard_weight; ptHard_hepmcxsections", {HistType::kTH2F, {{200, 0.0, 200.0}, {200, 0.0, 200.0}}});
    }

    if (doprocessCollisionsFromMc || doprocessCollisionsFromMcWeighted) {
      registry.add("h_fakecollisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    }
    if (doprocessCollisionsFromMcWeighted) {
      AxisSpec centAxis = {centralityBinning, "centrality (%)"};
      registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_collisions_weighted", "centrality vs mccollisions; centrality; collisions", {HistType::kTH2F, {centAxis, {4, 0.0, 4.0}}});
    }
    if (doprocessMcCollisionsWeighted) {
      AxisSpec centAxis = {centralityBinning, "centrality (%)"};
      registry.add("h_mccollisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_mccollisions_weighted", "centrality vs mccollisions; centrality; collisions", {HistType::kTH2F, {centAxis, {4, 0.0, 4.0}}});
      registry.add("h2_mccollision_pthardfromweight_pthardfromhepmcxsection_weighted", "ptHard from weight vs ptHard from HepMCXSections; ptHard_weight; ptHard_hepmcxsections", {HistType::kTH2F, {{200, 0.0, 200.0}, {200, 0.0, 200.0}}});
    }

    if (doprocessTrackSelectionHistograms) {
      registry.add("h_trackselplot_tpccrossedrows", "track selection variable: number of tpc crossed rows", {HistType::kTH1F, {{165, -0.5, 164.5}}});
      registry.add("h_trackselplot_tpccrossedrowsoverfindable", "track selection variable: ratio of of tpc crossed rows over number of findable clusters", {HistType::kTH1F, {{120, 0.0, 1.2}}});
      registry.add("h_trackselplot_chi2ncls_tpc", "track selection variable: Chi2 / cluster for the TPC track segment", {HistType::kTH1F, {{100, 0.0, 10.0}}});
      registry.add("h_trackselplot_chi2ncls_its", "track selection variable: Chi2 / cluster for the ITS track segment", {HistType::kTH1F, {{200, 0.0, 40.0}}});
      registry.add("h_trackselplot_dcaxy", "track selection variable: dca XY", {HistType::kTH1F, {{1000, -1.0, 1.0}}});
      registry.add("h_trackselplot_dcaz", "track selection variable: dca Z", {HistType::kTH1F, {{4000, -4.0, 4.0}}});

      registry.add("h2_trackselplot_pt_tpccrossedrows", "track selection variable: pt vs number of tpc crossed rows", {HistType::kTH2F, {{200, 0., 200.}, {165, -0.5, 164.5}}});
      registry.add("h2_trackselplot_pt_tpccrossedrowsoverfindable", "track selection variable: pt vs ratio of of tpc crossed rows over number of findable clusters", {HistType::kTH2F, {{200, 0., 200.}, {120, 0.0, 1.2}}});
      registry.add("h2_trackselplot_pt_chi2ncls_tpc", "track selection variable: pt vs Chi2 / cluster for the TPC track segment", {HistType::kTH2F, {{200, 0., 200.}, {100, 0.0, 10.0}}});
      registry.add("h2_trackselplot_pt_chi2ncls_its", "track selection variable: pt vs Chi2 / cluster for the ITS track segment", {HistType::kTH2F, {{200, 0., 200.}, {200, 0.0, 40.0}}});
      registry.add("h2_trackselplot_pt_dcaxy", "track selection variable: pt vs dca XY", {HistType::kTH2F, {{200, 0., 200.}, {1000, -1.0, 1.0}}});
      registry.add("h2_trackselplot_pt_dcaz", "track selection variable: pt vs dca Z", {HistType::kTH2F, {{200, 0., 200.}, {4000, -4.0, 4.0}}});
    }
  }

  Preslice<aod::JetTracksMCD> tracksPerJCollision = o2::aod::jtrack::collisionId;
  PresliceUnsorted<aod::JetTracksMCD> tracksPerJParticles = o2::aod::jmctracklb::mcParticleId;

  // filters for processTracks QA functions only:
  Filter trackCuts = (aod::jtrack::pt >= trackQAPtMin && aod::jtrack::pt < trackQAPtMax && aod::jtrack::eta > trackQAEtaMin && aod::jtrack::eta < trackQAEtaMax);
  Filter particleCuts = (aod::jmcparticle::pt >= trackQAPtMin && aod::jmcparticle::pt < trackQAPtMax && aod::jmcparticle::eta > trackQAEtaMin && aod::jmcparticle::eta < trackQAEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);

  void processEFficiencyPurity(soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>::iterator const& mcCollision,
                               soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                               soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, // smallgroups gives only the collisions associated to the current mccollision, thanks to the mccollisionlabel pre-integrated in jetcollisionsmcd
                               soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs> const& jetTracks,
                               soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const&,
                               JetParticlesWithOriginal const& jMcParticles)
  {
    // missing:
    //   * constexpr auto hasCentrality = CollisionMCRecTableCentFT0C::template contains<aod::CentFT0Cs>();
    //           if constexpr (hasCentrality) {
    // At the moment, are only counted mc particles from mc collisions that have at least one reconstructed collision that passes the chosen event selection. Thus, the reconstruction efficiency of mccollision is not counted in this tracking efficiency.

    registry.fill(HIST("hMcCollCutsCounts"), 0.5); // all mcCollisions

    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 1.5); // mcCollision.posZ() condition

    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 2.5); // mcCollisions with at least one reconstructed collision

    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 3.5); // split mcCollisions condition

    float centrality = -1;
    bool hasSel8Coll = false;
    bool centralityCheck = false;
    bool occupancyCheck = false;
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly || acceptSplitCollisions == NonSplitOnly) {    // check only that the first reconstructed collision passes the check (for the NonSplitOnly case, there's only one associated collision)
      if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
        hasSel8Coll = true;
      }
      if (!checkOccupancy || ((trackOccupancyInTimeRangeMin < collisions.begin().trackOccupancyInTimeRange()) && (collisions.begin().trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax))) { // check occupancy only in GP Pb-Pb MC
        occupancyCheck = true;
      }
      centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();
      if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
        centralityCheck = true;
      }
    } else if (acceptSplitCollisions == SplitOkCheckAnyAssocColl) { // check that at least one of the reconstructed collisions passes the checks
      for (auto const& collision : collisions) {
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
          hasSel8Coll = true;
        }
        if (!checkOccupancy || ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax))) { // check occupancy only in GP Pb-Pb MC
          occupancyCheck = true;
        }
        centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
        if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // effect unclear if mcColl is split
          centralityCheck = true;
        }
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 4.5); // at least one of the reconstructed collisions associated with this mcCollision is selected

    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C(); mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    // if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
    //   return;
    // }
    if (!centralityCheck) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 5.5); // at least one of the reconstructed collisions associated with this mcCollision is selected with regard to centrality

    float pTHat = getPtHatFromHepMCXSection ? mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(mcCollision.weight(), 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 6.5); // ptHat condition

    if (checkOccupancy) {
      if (!occupancyCheck) {
        return;
      }
      registry.fill(HIST("hMcCollCutsCounts"), 7.5);
    }

    for (auto const& jMcParticle : jMcParticles) {
      registry.fill(HIST("hMcPartCutsCounts"), 0.5); // allPartsInSelMcColl

      if (!isChargedParticle(jMcParticle.pdgCode())) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 1.5); // isCharged

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());

      if (checkPrimaryPart && !jMcParticle.isPhysicalPrimary()) { // global tracks should be mostly primaries
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 2.5); // isPrimary

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());

      registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());

      if ((std::abs(jMcParticle.eta()) < trackEtaAcceptanceCountQA)) { // removed from actual cuts for now because all the histograms have an eta axis
        registry.fill(HIST("hMcPartCutsCounts"), 3.5);                 // etaAccept // not actually applied here but it will give an idea of what will be done in the post processing
      }
    }

    std::vector<int> seenMcParticlesVector; // is reset every mc collision

    int splitCollCounter = 0;
    for (auto const& collision : collisions) {
      splitCollCounter++;
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly && splitCollCounter > 1) {
        return;
      }

      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents) || !(std::abs(collision.posZ()) < vertexZCut)) {
        continue;
      }

      auto collTracks = jetTracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      for (auto const& track : collTracks) {
        registry.fill(HIST("hTrackCutsCounts"), 0.5);

        if (!isAcceptedTrack(track)) {
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 1.5);

        if (!track.has_mcParticle()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi());

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi());
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 2.5);

        auto jMcParticleFromTrack = track.mcParticle_as<JetParticlesWithOriginal>();
        if (!jMcParticleFromTrack.isPhysicalPrimary()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

          if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
            registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split_nonprimary"), track.pt(), track.eta(), track.phi());
            registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

            registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_nonprimary"), track.pt(), track.eta(), track.phi());
            registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
          } else {
            seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
          }

          continue;
        }

        registry.fill(HIST("hTrackCutsCounts"), 3.5);

        registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        registry.fill(HIST("h2_particle_pt_track_pt_residual_associatedtrack_primary"), jMcParticleFromTrack.pt(), (jMcParticleFromTrack.pt() - track.pt()) / jMcParticleFromTrack.pt());

        registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        registry.fill(HIST("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary"), jMcParticleFromTrack.pt(), (jMcParticleFromTrack.pt() - track.pt()) / jMcParticleFromTrack.pt());

        if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_primary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        } else {
          seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
        }

        if (std::abs(jMcParticleFromTrack.eta()) < trackEtaAcceptanceCountQA) { // not actually applied here but it will give an idea of what will be done in the post processing
          registry.fill(HIST("hTrackCutsCounts"), 4.5);
        }
      }
    }
  }
  PROCESS_SWITCH(TrackEfficiency, processEFficiencyPurity, "Histograms for efficiency and purity quantities", true);

  void processEFficiencyPurityWeighted(soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>::iterator const& mcCollision,
                                       soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                                       soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, // smallgroups gives only the collisions associated to the current mccollision, thanks to the mccollisionlabel pre-integrated in jetcollisionsmcd
                                       soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs> const& jetTracks,
                                       soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const&,
                                       JetParticlesWithOriginal const& jMcParticles)
  {
    // missing:
    //   * constexpr auto hasCentrality = CollisionMCRecTableCentFT0C::template contains<aod::CentFT0Cs>();
    //           if constexpr (hasCentrality) {
    // At the moment, are only counted mc particles from mc collisions that have at least one reconstructed collision that passes the chosen event selection. Thus, the reconstruction efficiency of mccollision is not counted in this tracking efficiency.

    registry.fill(HIST("hMcCollCutsCounts"), 0.5, mcCollision.weight()); // all mcCollisions

    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 1.5, mcCollision.weight()); // mcCollision.posZ() condition

    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 2.5, mcCollision.weight()); // mcCollisions with at least one reconstructed collision

    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 3.5, mcCollision.weight()); // split mcCollisions condition

    float centrality = -1;
    bool hasSel8Coll = false;
    bool centralityCheck = false;
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly || acceptSplitCollisions == NonSplitOnly) {    // check only that the first reconstructed collision passes the check (for the NonSplitOnly case, there's only one associated collision)
      if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
        hasSel8Coll = true;
      }
      centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();
      if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
        centralityCheck = true;
      }
    } else if (acceptSplitCollisions == SplitOkCheckAnyAssocColl) { // check that at least one of the reconstructed collisions passes the checks
      for (auto const& collision : collisions) {
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
          hasSel8Coll = true;
        }
        centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
        if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
          centralityCheck = true;
        }
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 4.5, mcCollision.weight()); // at least one of the reconstructed collisions associated with this mcCollision is selected

    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C();  mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    // if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
    //   return;
    // }
    if (!centralityCheck) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 5.5, mcCollision.weight()); // centrality condition

    float mcCollEventWeight = mcCollision.weight();
    float pTHat = simPtRef / (std::pow(mcCollEventWeight, 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 6.5, mcCollision.weight()); // ptHat condition

    for (auto const& jMcParticle : jMcParticles) {
      registry.fill(HIST("hMcPartCutsCounts"), 0.5, mcCollision.weight()); // allPartsInSelMcColl

      if (!isChargedParticle(jMcParticle.pdgCode())) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 1.5, mcCollision.weight()); // isCharged

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi(), mcCollEventWeight);

      if (checkPrimaryPart && !jMcParticle.isPhysicalPrimary()) { // global tracks should be mostly primaries
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 2.5, mcCollision.weight()); // isPrimary

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi(), mcCollEventWeight);

      registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi(), mcCollEventWeight);

      if ((std::abs(jMcParticle.eta()) < trackEtaAcceptanceCountQA)) { // removed from actual cuts for now because all the histograms have an eta axis
        registry.fill(HIST("hMcPartCutsCounts"), 3.5, mcCollision.weight()); // etaAccept // not actually applied here but it will give an idea of what will be done in the post processing
      }
    }

    std::vector<int> seenMcParticlesVector; // is reset every mc collision

    int splitCollCounter = 0;
    for (auto const& collision : collisions) {
      splitCollCounter++;
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly && splitCollCounter > 1) {
        return;
      }

      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents) || !(std::abs(collision.posZ()) < vertexZCut)) {
        continue;
      }

      auto collTracks = jetTracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      for (auto const& track : collTracks) {
        registry.fill(HIST("hTrackCutsCounts"), 0.5, mcCollision.weight());

        if (!isAcceptedTrack(track)) {
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 1.5, mcCollision.weight());

        if (!track.has_mcParticle()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi(), mcCollEventWeight); // weight attribution here not trivial; I use the one of the current mcCollision, but track belongs to no collision; what should be its weight? could be a moot point but algo has complained about invalid index for mcParticle if I put th etrueTrackCollEventWeight before this cut

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi(), mcCollEventWeight);
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 2.5, mcCollision.weight());

        if (track.pt() > pTHatMaxFractionMCD * pTHat) {
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 3.5, mcCollision.weight());

        auto mcParticle = track.mcParticle_as<JetParticlesWithOriginal>();
        auto trueTrackMcCollision = mcParticle.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>();
        float trueTrackCollEventWeight = useTrueTrackWeight ? trueTrackMcCollision.weight() : mcCollEventWeight;

        auto jMcParticleFromTrack = track.mcParticle_as<JetParticlesWithOriginal>();
        if (!jMcParticleFromTrack.isPhysicalPrimary()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
          registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);

          if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
            registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split_nonprimary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
            registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);

            registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_nonprimary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
            registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);
          } else {
            seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
          }

          continue;
        }

        registry.fill(HIST("hTrackCutsCounts"), 4.5, mcCollision.weight());

        registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
        registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);
        registry.fill(HIST("h2_particle_pt_track_pt_residual_associatedtrack_primary"), jMcParticleFromTrack.pt(), (jMcParticleFromTrack.pt() - track.pt()) / jMcParticleFromTrack.pt(), trueTrackCollEventWeight);

        registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
        registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);
        registry.fill(HIST("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary"), jMcParticleFromTrack.pt(), (jMcParticleFromTrack.pt() - track.pt()) / jMcParticleFromTrack.pt(), trueTrackCollEventWeight);

        if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_primary"), track.pt(), track.eta(), track.phi(), trueTrackCollEventWeight);
          registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi(), trueTrackCollEventWeight);
        } else {
          seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
        }

        if (std::abs(jMcParticleFromTrack.eta()) < trackEtaAcceptanceCountQA) { // not actually applied here but it will give an idea of what will be done in the post processing
          registry.fill(HIST("hTrackCutsCounts"), 5.5, mcCollision.weight());
        }
      }
    }
  }
  PROCESS_SWITCH(TrackEfficiency, processEFficiencyPurityWeighted, "Histograms for efficiency and purity quantities for weighted simulations", false);

  void processTracksFromData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                             soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>> const& jetTracks,
                             soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }

    fillTrackHistograms(collision, jetTracks);
  }
  PROCESS_SWITCH(TrackEfficiency, processTracksFromData, "QA for charged tracks in data", false);

  void processTracksFromMc(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                           soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const&,
                           soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                           soa::Filtered<soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>> const& jetTracks,
                           soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const&)
  {
    if (!collision.has_mcCollision()) { // the collision is fake and has no associated mc coll; skip as .mccollision() cannot be called
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }

    float pTHat = getPtHatFromHepMCXSection ? collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(collision.mcCollision().weight(), 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }

    fillTrackHistograms(collision, jetTracks);
  }
  PROCESS_SWITCH(TrackEfficiency, processTracksFromMc, "QA for charged tracks in MC without weights", false);

  void processTracksFromMcWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                                   soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const&,
                                   soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                                   soa::Filtered<soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>> const& jetTracks,
                                   soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const&)
  {
    if (!collision.has_mcCollision()) { // the collision is fake and has no associated mc coll; skip as .mccollision() cannot be called
      return;
    }
    float eventWeight = collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().weight();
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }

    float pTHat = getPtHatFromHepMCXSection ? collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }

    fillTrackHistograms(collision, jetTracks, eventWeight);
  }
  PROCESS_SWITCH(TrackEfficiency, processTracksFromMcWeighted, "QA for charged tracks in weighted MC", false);

  void processParticles(soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>::iterator const& mcCollision,
                        soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                        soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                        soa::Filtered<aod::JetParticles> const& mcparticles,
                        soa::Filtered<aod::JetTracksMCD> const& tracks)
  {

    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return;
    }

    float pTHat = getPtHatFromHepMCXSection ? mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(mcCollision.weight(), 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }

    float centrality = -1;
    bool hasSel8Coll = false;
    bool centralityCheck = false;
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly || acceptSplitCollisions == NonSplitOnly) {    // check only that the first reconstructed collision passes the check (for the NonSplitOnly case, there's only one associated collision)
      if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
        hasSel8Coll = true;
      }
      centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();
      if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
        centralityCheck = true;
      }
    } else if (acceptSplitCollisions == SplitOkCheckAnyAssocColl) { // check that at least one of the reconstructed collisions passes the checks
      for (auto const& collision : collisions) {
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
          hasSel8Coll = true;
        }
        centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
        if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
          centralityCheck = true;
        }
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C();  mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    // if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
    //   return;
    // }
    if (!centralityCheck) {
      return;
    }

    fillParticlesHistograms(mcCollision, collisions, mcparticles, tracks);
  }
  PROCESS_SWITCH(TrackEfficiency, processParticles, "QA for charged particles", false);

  void processParticlesWeighted(soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>::iterator const& mcCollision,
                                soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                                soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                soa::Filtered<aod::JetParticles> const& mcparticles,
                                soa::Filtered<aod::JetTracksMCD> const& tracks)
  {
    if (skipMBGapEvents && mcCollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }

    float eventWeight = mcCollision.weight();

    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return;
    }

    float pTHat = getPtHatFromHepMCXSection ? mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }

    float centrality = -1;
    bool hasSel8Coll = false;
    bool centralityCheck = false;
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly || acceptSplitCollisions == NonSplitOnly) {    // check only that the first reconstructed collision passes the check (for the NonSplitOnly case, there's only one associated collision)
      if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
        hasSel8Coll = true;
      }
      centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();
      if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
        centralityCheck = true;
      }
    } else if (acceptSplitCollisions == SplitOkCheckAnyAssocColl) { // check that at least one of the reconstructed collisions passes the checks
      for (auto const& collision : collisions) {
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
          hasSel8Coll = true;
        }
        centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
        if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
          centralityCheck = true;
        }
      }
    }
    if (!hasSel8Coll) {
      return;
    }

    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C();  mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    // if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
    //   return;
    // }
    if (!centralityCheck) {
      return;
    }

    fillParticlesHistograms(mcCollision, collisions, mcparticles, tracks, eventWeight);
  }
  PROCESS_SWITCH(TrackEfficiency, processParticlesWeighted, "QA for charged particles weighted", false);

  void processCollisionsFromData(soa::Filtered<aod::JetCollisions>::iterator const& collision)
  {
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();

    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 1.5);
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 2.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 3.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 3.5);
  }
  PROCESS_SWITCH(TrackEfficiency, processCollisionsFromData, "QA for reconstructed collisions in data", false);

  void processCollisionsFromMc(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                               soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const&,
                               soa::Join<aod::McCollisions, aod::HepMCXSections> const&)
  {
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();

    if (!collision.has_mcCollision()) { // the collision is fake and has no associated mc coll; skip as .mccollision() cannot be called
      registry.fill(HIST("h_fakecollisions"), 0.5);
      return;
    }
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 1.5);
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 2.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 3.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 3.5);

    float pTHat = getPtHatFromHepMCXSection ? collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(collision.mcCollision().weight(), 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }
    registry.fill(HIST("h_collisions"), 4.5);
    registry.fill(HIST("h2_centrality_collisions"), centrality, 4.5);
  }
  PROCESS_SWITCH(TrackEfficiency, processCollisionsFromMc, "QA for reconstructed collisions in MC without weights", false);

  void processCollisionsFromMcWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                                       soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const&,
                                       soa::Join<aod::McCollisions, aod::HepMCXSections> const&)
  {
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();

    if (!collision.has_mcCollision()) { // the collision is fake and has no associated mc coll; skip as .mccollision() cannot be called
      registry.fill(HIST("h_fakecollisions"), 0.5);
      return;
    }
    float eventWeight = collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h_collisions_weighted"), 2.5, eventWeight);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 3.5);
    registry.fill(HIST("h_collisions_weighted"), 3.5, eventWeight);

    float pTHat = getPtHatFromHepMCXSection ? collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>().mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }
    registry.fill(HIST("h_collisions"), 4.5);
    registry.fill(HIST("h_collisions_weighted"), 4.5, eventWeight);
  }
  PROCESS_SWITCH(TrackEfficiency, processCollisionsFromMcWeighted, "QA for reconstructed collisions in weighted MC", false);

  void processMcCollisions(soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>::iterator const& mcCollision,
                           soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                           soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
  {
    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C(); mcCollision.centFT0C() isn't filled at the moment; can be added back when it is

    float eventWeight = mcCollision.weight();
    float pTHat = getPtHatFromHepMCXSection ? mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    registry.fill(HIST("h2_mccollision_pthardfromweight_pthardfromhepmcxsection"), 10. / (std::pow(eventWeight, 1.0 / pTHatExponent)), mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard());

    float centrality = -1;
    bool hasSel8Coll = false;
    bool centralityCheck = false;
    if (collisions.size() > 1) {                                                                                 // remove and move the if block below under if (collisions.size() < 1) { when mccoll.centFt0C has been fixed
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly || acceptSplitCollisions == NonSplitOnly) {    // check only that the first reconstructed collision passes the check (for the NonSplitOnly case, there's only one associated collision)
        if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
          hasSel8Coll = true;
        }
        centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();
        if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
          centralityCheck = true;
        }
      } else if (acceptSplitCollisions == SplitOkCheckAnyAssocColl) { // check that at least one of the reconstructed collisions passes the checks
        for (auto const& collision : collisions) {
          if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
            hasSel8Coll = true;
          }
          centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
          if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
            centralityCheck = true;
          }
        }
      }
    }

    registry.fill(HIST("h_mccollisions"), 0.5);
    registry.fill(HIST("h2_centrality_mccollisions"), centrality, 0.5);

    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return;
    }

    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }
    registry.fill(HIST("h_mccollisions"), 1.5);
    registry.fill(HIST("h2_centrality_mccollisions"), centrality, 1.5);

    if (!hasSel8Coll) {
      return;
    }
    // if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) { mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    //   return;
    // }
    if (!centralityCheck) {
      return;
    }

    registry.fill(HIST("h_mccollisions"), 2.5);
    registry.fill(HIST("h2_centrality_mccollisions"), centrality, 2.5);
  }
  PROCESS_SWITCH(TrackEfficiency, processMcCollisions, "QA for McCollisions in MC without weights", false);

  void processMcCollisionsWeighted(soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>::iterator const& mcCollision,
                                   soa::Join<aod::McCollisions, aod::HepMCXSections> const&,
                                   soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
  {
    if (skipMBGapEvents && mcCollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }

    // float centrality = checkCentFT0M ? mcCollision.centFT0M() : mcCollision.centFT0C();  mcCollision.centFT0C() isn't filled at the moment; can be added back when it is

    float eventWeight = mcCollision.weight();
    float pTHat = getPtHatFromHepMCXSection ? mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard() : 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    registry.fill(HIST("h2_mccollision_pthardfromweight_pthardfromhepmcxsection"), 10. / (std::pow(eventWeight, 1.0 / pTHatExponent)), mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard());
    registry.fill(HIST("h2_mccollision_pthardfromweight_pthardfromhepmcxsection_weighted"), 10. / (std::pow(eventWeight, 1.0 / pTHatExponent)), mcCollision.mcCollision_as<soa::Join<aod::McCollisions, aod::HepMCXSections>>().ptHard(), eventWeight);

    float centrality = -1;
    bool hasSel8Coll = false;
    bool centralityCheck = false;
    if (collisions.size() > 1) {                                                                                 // remove and move the if block below under if (collisions.size() < 1) { when mccoll.centFt0C has been fixed
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly || acceptSplitCollisions == NonSplitOnly) {    // check only that the first reconstructed collision passes the check (for the NonSplitOnly case, there's only one associated collision)
        if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have their first associated collision not reconstructed
          hasSel8Coll = true;
        }
        centrality = checkCentFT0M ? collisions.begin().centFT0M() : collisions.begin().centFT0C();
        if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
          centralityCheck = true;
        }
      } else if (acceptSplitCollisions == SplitOkCheckAnyAssocColl) { // check that at least one of the reconstructed collisions passes the checks
        for (auto const& collision : collisions) {
          if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
            hasSel8Coll = true;
          }
          centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
          if (!cutCentrality || ((centralityMin < centrality) && (centrality < centralityMax))) { // mcCollision.centFT0C() isn't filled at the moment; can use it instead when it is added to O2Physics
            centralityCheck = true;
          }
        }
      }
    }

    registry.fill(HIST("h_mccollisions"), 0.5);
    registry.fill(HIST("h_mccollisions_weighted"), 0.5, eventWeight);
    registry.fill(HIST("h2_centrality_mccollisions"), centrality, 0.5);
    registry.fill(HIST("h2_centrality_mccollisions_weighted"), centrality, 0.5, eventWeight);

    if (!(std::abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    if (collisions.size() < 1) {
      return;
    }
    if (acceptSplitCollisions == NonSplitOnly && collisions.size() > 1) {
      return;
    }

    if (pTHat < ptHatMin || pTHat > ptHatMax) { // only allows mcCollisions with weight in between min and max
      return;
    }
    registry.fill(HIST("h_mccollisions"), 1.5);
    registry.fill(HIST("h_mccollisions_weighted"), 1.5, eventWeight);
    registry.fill(HIST("h2_centrality_mccollisions"), centrality, 1.5);
    registry.fill(HIST("h2_centrality_mccollisions_weighted"), centrality, 1.5, eventWeight);

    if (!hasSel8Coll) {
      return;
    }
    // if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) { mcCollision.centFT0C() isn't filled at the moment; can be added back when it is
    //   return;
    // }
    if (!centralityCheck) {
      return;
    }
    registry.fill(HIST("h_mccollisions"), 2.5);
    registry.fill(HIST("h_mccollisions_weighted"), 2.5, eventWeight);
    registry.fill(HIST("h2_centrality_mccollisions"), centrality, 2.5);
    registry.fill(HIST("h2_centrality_mccollisions_weighted"), centrality, 2.5, eventWeight);
  }
  PROCESS_SWITCH(TrackEfficiency, processMcCollisionsWeighted, "QA for McCollisions in weighted MC", false);

  void processTrackSelectionHistograms(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<aod::JetTracks, aod::JTrackPIs> const& jetTracks, soa::Join<aod::Tracks, aod::TracksExtra, o2::aod::TracksDCA> const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    if (cutCentrality && (centrality < centralityMin || centralityMax < centrality)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }

    for (auto const& jetTrack : jetTracks) {
      const auto& aodTrack = jetTrack.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>>();

      registry.fill(HIST("h_trackselplot_tpccrossedrows"), aodTrack.tpcNClsCrossedRows());
      registry.fill(HIST("h_trackselplot_tpccrossedrowsoverfindable"), aodTrack.tpcCrossedRowsOverFindableCls());
      registry.fill(HIST("h_trackselplot_chi2ncls_tpc"), aodTrack.tpcChi2NCl());
      registry.fill(HIST("h_trackselplot_chi2ncls_its"), aodTrack.itsChi2NCl());
      registry.fill(HIST("h_trackselplot_dcaxy"), aodTrack.dcaXY());
      registry.fill(HIST("h_trackselplot_dcaz"), aodTrack.dcaZ());

      registry.fill(HIST("h2_trackselplot_pt_tpccrossedrows"), aodTrack.pt(), aodTrack.tpcNClsCrossedRows());
      registry.fill(HIST("h2_trackselplot_pt_tpccrossedrowsoverfindable"), aodTrack.pt(), aodTrack.tpcCrossedRowsOverFindableCls());
      registry.fill(HIST("h2_trackselplot_pt_chi2ncls_tpc"), aodTrack.pt(), aodTrack.tpcChi2NCl());
      registry.fill(HIST("h2_trackselplot_pt_chi2ncls_its"), aodTrack.pt(), aodTrack.itsChi2NCl());
      registry.fill(HIST("h2_trackselplot_pt_dcaxy"), aodTrack.pt(), aodTrack.dcaXY());
      registry.fill(HIST("h2_trackselplot_pt_dcaz"), aodTrack.pt(), aodTrack.dcaZ());
    }
  }
  PROCESS_SWITCH(TrackEfficiency, processTrackSelectionHistograms, "plots distributions of variables that are cut on during track selection", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackEfficiency>(cfgc)};
}
