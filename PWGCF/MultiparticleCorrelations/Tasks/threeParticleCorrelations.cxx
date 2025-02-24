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
/// \file threeParticleCorrelations.cxx
/// \brief Task for producing particle correlations
/// \author Joey Staa <joey.staa@fysik.lu.se>

#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "RecoDecay.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ThreeParticleCorrelations {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Histogram registry
  HistogramRegistry rMECorrRegistry{"MECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rSECorrRegistry{"SECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rMCRegistry{"MCRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rPhiStarRegistry{"PhiStarRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rQARegistry{"QARegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Collision & Event filters
  Filter collCent = aod::cent::centFT0M > 0.0f && aod::cent::centFT0M < 90.0f;
  Filter collZvtx = nabs(aod::collision::posZ) < 7.0f;
  Filter mcCollZvtx = nabs(aod::mccollision::posZ) < 7.0f;
  Filter evSelect = aod::evsel::sel8 == true;

  // V0 filters
  Filter v0Pt = aod::v0data::pt > 0.6f && aod::v0data::pt < 12.0f;
  Filter v0Eta = nabs(aod::v0data::eta) < 0.72f;

  // Track filters
  Filter trackPt = aod::track::pt > 0.2f && aod::track::pt < 3.0f;
  Filter trackEta = nabs(aod::track::eta) < 0.8f;
  Filter globalTracks = requireGlobalTrackInFilter();

  // Particle filters
  Filter particleEta = nabs(aod::mcparticle::eta) < 0.8f;

  // Table aliases - Data
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;
  using MyFilteredV0s = soa::Filtered<aod::V0Datas>;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                   aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;

  // Table aliases - MC
  using MyFilteredMCGenCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::CentFT0Ms>>;
  using MyFilteredMCGenCollision = MyFilteredMCGenCollisions::iterator;
  using MyFilteredMCParticles = soa::Filtered<aod::McParticles>;
  using MyFilteredMCRecCollision = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>>::iterator;
  using MyFilteredMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels,
                                                     aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                     aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;

  // Mixed-events binning policy
  SliceCache cache;
  ConfigurableAxis confCentBins{"confCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "ME Centrality binning"};
  ConfigurableAxis confZvtxBins{"confZvtxBins", {VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 0.0f, 1.0f, 3.0f, 5.0f, 7.0f}, "ME Zvtx binning"};
  using BinningType = ColumnBinningPolicy<aod::cent::CentFT0M, aod::collision::PosZ>;

  BinningType collBinning{{confCentBins, confZvtxBins}, true};
  Pair<MyFilteredCollisions, MyFilteredV0s, MyFilteredTracks, BinningType> pairData{collBinning, 5, -1, &cache};
  SameKindPair<MyFilteredMCGenCollisions, MyFilteredMCParticles, BinningType> pairMC{collBinning, 5, -1, &cache};

  // Process configurables
  Configurable<bool> confFilterSwitch{"confFilterSwitch", false, "Switch for the fakeV0Filter function"};

  // Particle masses
  double massLambda = constants::physics::MassLambda0;
  double dGaussSigma = 0.0021;

  // Efficiency histograms
  TH1D** hEffPions = new TH1D*[2];
  TH1D** hEffKaons = new TH1D*[2];
  TH1D** hEffProtons = new TH1D*[2];

  // Correlation variables
  int triggSign, assocSign;
  double candMass;
  double* assocPID;

  double deltaPhi, deltaEta;

  //==========================================================================================================================================================================================================================================================================

  void init(InitContext const&)
  {

    // Histograms axes
    const AxisSpec centralityAxis{confCentBins};
    const AxisSpec zvtxAxis{confZvtxBins};
    const AxisSpec phiAxis{36, (-1. / 2) * constants::math::PI, (3. / 2) * constants::math::PI};
    const AxisSpec etaAxis{32, -1.52, 1.52};
    const AxisSpec v0PtAxis{114, 0.6, 12};
    const AxisSpec trackPtAxis{28, 0.2, 3};
    const AxisSpec lambdaInvMassAxis{100, 1.08, 1.16};

    // QA & PID
    rQARegistry.add("hTrackPt", "hTrackPt", {HistType::kTH1D, {{100, 0, 4}}});
    rQARegistry.add("hTrackEta", "hTrackEta", {HistType::kTH1D, {{100, -1, 1}}});
    rQARegistry.add("hTrackPhi", "hTrackPhi", {HistType::kTH1D, {{100, (-1. / 2) * constants::math::PI, (5. / 2) * constants::math::PI}}});
    rQARegistry.add("hEventCentrality", "hEventCentrality", {HistType::kTH1D, {{centralityAxis}}});
    rQARegistry.add("hEventZvtx", "hEventZvtx", {HistType::kTH1D, {{zvtxAxis}}});

    rQARegistry.add("hdEdx", "hdEdx", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxPion", "hdEdxPion", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxKaon", "hdEdxKaon", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxProton", "hdEdxProton", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hBeta", "hBeta", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaPion", "hBetaPion", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaKaon", "hBetaKaon", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaProton", "hBetaProton", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hNSigmaPion", "hNSigmaPion", {HistType::kTH2D, {{201, -5.025, 5.025}, {201, -5.025, 5.025}}});
    rQARegistry.add("hNSigmaKaon", "hNSigmaKaon", {HistType::kTH2D, {{201, -5.025, 5.025}, {201, -5.025, 5.025}}});
    rQARegistry.add("hNSigmaProton", "hNSigmaProton", {HistType::kTH2D, {{201, -5.025, 5.025}, {201, -5.025, 5.025}}});

    rQARegistry.add("hTPCPion", "hTPCPion", {HistType::kTH2D, {{trackPtAxis}, {241, -6, 6}}});
    rQARegistry.add("hTPCKaon", "hTPCKaon", {HistType::kTH2D, {{trackPtAxis}, {241, -6, 6}}});
    rQARegistry.add("hTPCProton", "hTPCProton", {HistType::kTH2D, {{trackPtAxis}, {241, -6, 6}}});
    rQARegistry.add("hTOFPion", "hTOFPion", {HistType::kTH2D, {{trackPtAxis}, {1000, -50, 50}}});
    rQARegistry.add("hTOFKaon", "hTOFKaon", {HistType::kTH2D, {{trackPtAxis}, {1000, -50, 50}}});
    rQARegistry.add("hTOFProton", "hTOFProton", {HistType::kTH2D, {{trackPtAxis}, {1000, -50, 50}}});

    rQARegistry.add("hInvMassLambda", "hInvMassLambda", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hInvMassAntiLambda", "hInvMassAntiLambda", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hNLambdas_MC", "hNLambdas_MC", {HistType::kTH3D, {{2, -2, 2}, {v0PtAxis}, {centralityAxis}}});

    // PhiStar
    rPhiStarRegistry.add("hSEProtonPreCut", "hSEProtonPreCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});
    rPhiStarRegistry.add("hSEProtonPostCut", "hSEProtonPostCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});
    rPhiStarRegistry.add("hMEProtonPreCut", "hMEProtonPreCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});
    rPhiStarRegistry.add("hMEProtonPostCut", "hMEProtonPostCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});

    // Efficiency
    rMCRegistry.add("hGenerated", "hGenerated", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hGenPionP", "hGenPionP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hGenPionN", "hGenPionN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hGenKaonP", "hGenKaonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hGenKaonN", "hGenKaonN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hGenProtonP", "hGenProtonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hGenProtonN", "hGenProtonN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hReconstructed", "hReconstructed", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hRecPionP", "hRecPionP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hRecPionN", "hRecPionN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hRecKaonP", "hRecKaonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hRecKaonN", "hRecKaonN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hRecProtonP", "hRecProtonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hRecProtonN", "hRecProtonN", {HistType::kTH1D, {trackPtAxis}});

    // Purity
    rMCRegistry.add("hSelectPionP", "hSelectPionP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hSelectPionN", "hSelectPionN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hSelectKaonP", "hSelectKaonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hSelectKaonN", "hSelectKaonN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hSelectProtonP", "hSelectProtonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hSelectProtonN", "hSelectProtonN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hTrueSelectPionP", "hTrueSelectPionP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hTrueSelectPionN", "hTrueSelectPionN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hTrueSelectKaonP", "hTrueSelectKaonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hTrueSelectKaonN", "hTrueSelectKaonN", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hTrueSelectProtonP", "hTrueSelectProtonP", {HistType::kTH1D, {trackPtAxis}});
    rMCRegistry.add("hTrueSelectProtonN", "hTrueSelectProtonN", {HistType::kTH1D, {trackPtAxis}});

    // Correlations
    rSECorrRegistry.add("hSameLambdaPion_SGNL", "Same-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_SB", "Same-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_SGNL", "Same-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_SB", "Same-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_SGNL", "Same-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_SB", "Same-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_MC", "Same-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_MC", "Same-event #Lambda - K correlator (MC)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_MC", "Same-event #Lambda - p correlator (MC)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    rMECorrRegistry.add("hMixLambdaPion_SGNL", "Mixed-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_SB", "Mixed-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_SGNL", "Mixed-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_SB", "Mixed-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_SGNL", "Mixed-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_SB", "Mixed-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_MC", "Mixed-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_MC", "Mixed-event #Lambda - K correlator (MC)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_MC", "Mixed-event #Lambda - p correlator (MC)", {HistType::kTHnSparseD, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    TList* efficiencyList = ccdb->getForTimeStamp<TList>("Users/j/jstaa/Efficiency/ChargedParticles", 1);
    hEffPions[0] = static_cast<TH1D*>(efficiencyList->FindObject("hEfficiencyPionP"));
    hEffPions[1] = static_cast<TH1D*>(efficiencyList->FindObject("hEfficiencyPionN"));
    hEffKaons[0] = static_cast<TH1D*>(efficiencyList->FindObject("hEfficiencyKaonP"));
    hEffKaons[1] = static_cast<TH1D*>(efficiencyList->FindObject("hEfficiencyKaonN"));
    hEffProtons[0] = static_cast<TH1D*>(efficiencyList->FindObject("hEfficiencyProtonP"));
    hEffProtons[1] = static_cast<TH1D*>(efficiencyList->FindObject("hEfficiencyProtonN"));
  }

  //==========================================================================================================================================================================================================================================================================

  void processSame(MyFilteredCollision const& collision, MyFilteredV0s const& v0s, MyFilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());
    rQARegistry.fill(HIST("hEventCentrality"), collision.centFT0M());
    rQARegistry.fill(HIST("hEventZvtx"), collision.posZ());

    // Start of the Track QA
    for (const auto& track : tracks) {
      rQARegistry.fill(HIST("hTPCPion"), track.pt(), track.tpcNSigmaPi());
      rQARegistry.fill(HIST("hTPCKaon"), track.pt(), track.tpcNSigmaKa());
      rQARegistry.fill(HIST("hTPCProton"), track.pt(), track.tpcNSigmaPr());
      if (track.hasTOF()) {
        rQARegistry.fill(HIST("hTOFPion"), track.pt(), track.tofNSigmaPi());
        rQARegistry.fill(HIST("hTOFKaon"), track.pt(), track.tofNSigmaKa());
        rQARegistry.fill(HIST("hTOFProton"), track.pt(), track.tofNSigmaPr());
      }

      if (trackFilters(track)) {
        assocPID = trackPID(track);
        rQARegistry.fill(HIST("hTrackPt"), track.pt());
        rQARegistry.fill(HIST("hTrackEta"), track.eta());
        rQARegistry.fill(HIST("hTrackPhi"), track.phi());
        rQARegistry.fill(HIST("hdEdx"), track.pt(), track.tpcSignal());
        rQARegistry.fill(HIST("hBeta"), track.pt(), track.beta());
        if (assocPID[0] == 0.0) { // Pions
          rQARegistry.fill(HIST("hdEdxPion"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaPion"), track.pt(), track.beta());
          rQARegistry.fill(HIST("hNSigmaPion"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        } else if (assocPID[0] == 1.0) { // Kaons
          rQARegistry.fill(HIST("hdEdxKaon"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaKaon"), track.pt(), track.beta());
          rQARegistry.fill(HIST("hNSigmaKaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        } else if (assocPID[0] == 2.0) { // Protons
          rQARegistry.fill(HIST("hdEdxProton"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaProton"), track.pt(), track.beta());
          rQARegistry.fill(HIST("hNSigmaProton"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
      }
    }
    // End of the Track QA

    // Start of the Same-Event Correlations
    for (const auto& trigger : v0s) {
      if (v0Filters(trigger)) {

        triggSign = v0Sign(trigger);
        if (triggSign == 1) {
          candMass = trigger.mLambda();
          rQARegistry.fill(HIST("hInvMassLambda"), trigger.mLambda(), trigger.pt(), collision.centFT0M());
        } else if (triggSign == -1) {
          candMass = trigger.mAntiLambda();
          rQARegistry.fill(HIST("hInvMassAntiLambda"), trigger.mAntiLambda(), trigger.pt(), collision.centFT0M());
        }

        for (const auto& associate : tracks) {
          if (trackFilters(associate)) {
            if (correlationFilters(trigger, associate) && radialDistanceFilter(trigger, associate, bField, false) && fakeV0Filter(trigger, associate)) {

              assocPID = trackPID(associate);
              deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
              deltaEta = trigger.eta() - associate.eta();

              if (candMass >= massLambda - 4 * dGaussSigma && candMass <= massLambda + 4 * dGaussSigma) {
                if (assocPID[0] == 0.0) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SGNL"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
                } else if (assocPID[0] == 1.0) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SGNL"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
                } else if (assocPID[0] == 2.0) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SGNL"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
                }
              } else if (candMass >= massLambda - 8 * dGaussSigma && candMass <= massLambda + 8 * dGaussSigma) {
                if (assocPID[0] == 0.0) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SB"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
                } else if (assocPID[0] == 1.0) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SB"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
                } else if (assocPID[0] == 2.0) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SB"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
                }
              }
            }
          }
        }
      }
    }
    // End of the Same-Event Correlations
  }

  void processMixed(MyFilteredCollisions const&, MyFilteredV0s const&, MyFilteredTracks const&, aod::BCsWithTimestamps const&)
  {

    // Start of the Mixed-events Correlations
    for (const auto& [coll_1, v0_1, coll_2, track_2] : pairData) {

      auto bc = coll_1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      for (const auto& [trigger, associate] : soa::combinations(soa::CombinationsFullIndexPolicy(v0_1, track_2))) {
        if (v0Filters(trigger) && trackFilters(associate)) {
          if (radialDistanceFilter(trigger, associate, bField, true) && fakeV0Filter(trigger, associate)) {

            triggSign = v0Sign(trigger);
            if (triggSign == 1) {
              candMass = trigger.mLambda();
            } else if (triggSign == -1) {
              candMass = trigger.mAntiLambda();
            }

            assocPID = trackPID(associate);
            deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
            deltaEta = trigger.eta() - associate.eta();

            if (candMass >= massLambda - 4 * dGaussSigma && candMass <= massLambda + 4 * dGaussSigma) {
              if (assocPID[0] == 0.0) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SGNL"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
              } else if (assocPID[0] == 1.0) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SGNL"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
              } else if (assocPID[0] == 2.0) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SGNL"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
              }
            } else if (candMass >= massLambda - 8 * dGaussSigma && candMass <= massLambda + 8 * dGaussSigma) {
              if (assocPID[0] == 0.0) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SB"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
              } else if (assocPID[0] == 1.0) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SB"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
              } else if (assocPID[0] == 2.0) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SB"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
              }
            }
          }
        }
      }
    }
    // End of the Mixed-events Correlations
  }

  void processMCSame(MyFilteredMCGenCollision const& collision, MyFilteredMCParticles const& particles)
  {

    Partition<MyFilteredMCParticles> mcTriggers = (aod::mcparticle::pdgCode == static_cast<int>(kLambda0) || aod::mcparticle::pdgCode == static_cast<int>(kLambda0Bar)) && aod::mcparticle::pt > 0.6f && aod::mcparticle::pt < 12.0f && nabs(aod::mcparticle::eta) < 0.72f;
    Partition<MyFilteredMCParticles> mcAssociates = (((aod::mcparticle::pdgCode == static_cast<int>(kPiPlus) || aod::mcparticle::pdgCode == static_cast<int>(kPiMinus)) && aod::mcparticle::pt > 0.3f && aod::mcparticle::pt < 2.3f) ||
                                                     ((aod::mcparticle::pdgCode == static_cast<int>(kKPlus) || aod::mcparticle::pdgCode == static_cast<int>(kKMinus)) && aod::mcparticle::pt > 0.5f && aod::mcparticle::pt < 2.5f) ||
                                                     ((aod::mcparticle::pdgCode == static_cast<int>(kProton) || aod::mcparticle::pdgCode == static_cast<int>(kProtonBar)) && aod::mcparticle::pt > 0.5f));
    mcTriggers.bindTable(particles);
    mcAssociates.bindTable(particles);

    // Start of the MC Same-Event Correlations
    for (const auto& trigger : mcTriggers) {
      if (trigger.isPhysicalPrimary()) {

        if (trigger.pdgCode() > 0) {
          triggSign = 1;
          rQARegistry.fill(HIST("hNLambdas_MC"), 1, trigger.pt(), collision.centFT0M());
        } else if (trigger.pdgCode() < 0) {
          triggSign = -1;
          rQARegistry.fill(HIST("hNLambdas_MC"), -1, trigger.pt(), collision.centFT0M());
        }

        for (const auto& associate : mcAssociates) {
          if (associate.isPhysicalPrimary()) {

            if (associate.pdgCode() > 0) {
              assocSign = 1;
            } else if (associate.pdgCode() < 0) {
              assocSign = -1;
            }

            deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
            deltaEta = trigger.eta() - associate.eta();

            if (std::abs(associate.pdgCode()) == kPiPlus) {
              rSECorrRegistry.fill(HIST("hSameLambdaPion_MC"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, assocSign);
            } else if (std::abs(associate.pdgCode()) == kKPlus) {
              rSECorrRegistry.fill(HIST("hSameLambdaKaon_MC"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, assocSign);
            } else if (std::abs(associate.pdgCode()) == kProton) {
              rSECorrRegistry.fill(HIST("hSameLambdaProton_MC"), deltaPhi, deltaEta, collision.centFT0M(), collision.posZ(), triggSign, assocSign);
            }
          }
        }
      }
    }
    // End of the MC Same-Event Correlations
  }

  void processMCMixed(MyFilteredMCGenCollisions const&, MyFilteredMCParticles const&)
  {

    // Start of the MC Mixed-events Correlations
    for (const auto& [coll_1, v0_1, coll_2, track_2] : pairMC) {
      Partition<MyFilteredMCParticles> mcTriggers = (aod::mcparticle::pdgCode == static_cast<int>(kLambda0) || aod::mcparticle::pdgCode == static_cast<int>(kLambda0Bar)) && aod::mcparticle::pt > 0.6f && aod::mcparticle::pt < 12.0f && nabs(aod::mcparticle::eta) < 0.72f;
      Partition<MyFilteredMCParticles> mcAssociates = (((aod::mcparticle::pdgCode == static_cast<int>(kPiPlus) || aod::mcparticle::pdgCode == static_cast<int>(kPiMinus)) && aod::mcparticle::pt > 0.3f && aod::mcparticle::pt < 2.3f) ||
                                                       ((aod::mcparticle::pdgCode == static_cast<int>(kKPlus) || aod::mcparticle::pdgCode == static_cast<int>(kKMinus)) && aod::mcparticle::pt > 0.5f && aod::mcparticle::pt < 2.5f) ||
                                                       ((aod::mcparticle::pdgCode == static_cast<int>(kProton) || aod::mcparticle::pdgCode == static_cast<int>(kProtonBar)) && aod::mcparticle::pt > 0.5f));
      mcTriggers.bindTable(v0_1);
      mcAssociates.bindTable(track_2);

      for (const auto& [trigger, associate] : soa::combinations(soa::CombinationsFullIndexPolicy(mcTriggers, mcAssociates))) {
        if (trigger.isPhysicalPrimary() && associate.isPhysicalPrimary()) {

          if (trigger.pdgCode() > 0) {
            triggSign = 1;
          } else if (trigger.pdgCode() < 0) {
            triggSign = -1;
          }
          if (associate.pdgCode() > 0) {
            assocSign = 1;
          } else if (associate.pdgCode() < 0) {
            assocSign = -1;
          }

          deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
          deltaEta = trigger.eta() - associate.eta();

          if (std::abs(associate.pdgCode()) == kPiPlus) {
            rMECorrRegistry.fill(HIST("hMixLambdaPion_MC"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, assocSign);
          } else if (std::abs(associate.pdgCode()) == kKPlus) {
            rMECorrRegistry.fill(HIST("hMixLambdaKaon_MC"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, assocSign);
          } else if (std::abs(associate.pdgCode()) == kProton) {
            rMECorrRegistry.fill(HIST("hMixLambdaProton_MC"), deltaPhi, deltaEta, coll_1.centFT0M(), coll_1.posZ(), triggSign, assocSign);
          }
        }
      }
    }
    // End of the MC Mixed-events Correlations
  }

  void processMCGen(MyFilteredMCGenCollision const&, MyFilteredMCParticles const& particles)
  {

    Partition<MyFilteredMCParticles> mcParticles = aod::mcparticle::pt > 0.2f && aod::mcparticle::pt < 3.0f;
    mcParticles.bindTable(particles);

    // Start of the Monte-Carlo generated QA
    for (const auto& particle : mcParticles) {
      if (particle.isPhysicalPrimary()) {

        // Efficiency - Generated
        rMCRegistry.fill(HIST("hGenerated"), particle.pt());
        if (particle.pdgCode() == kPiPlus) { // Pos pions
          rMCRegistry.fill(HIST("hGenPionP"), particle.pt());
        } else if (particle.pdgCode() == kPiMinus) { // Neg pions
          rMCRegistry.fill(HIST("hGenPionN"), particle.pt());
        } else if (particle.pdgCode() == kKPlus) { // Pos kaons
          rMCRegistry.fill(HIST("hGenKaonP"), particle.pt());
        } else if (particle.pdgCode() == kKMinus) { // Neg kaons
          rMCRegistry.fill(HIST("hGenKaonN"), particle.pt());
        } else if (particle.pdgCode() == kProton) { // Pos protons
          rMCRegistry.fill(HIST("hGenProtonP"), particle.pt());
        } else if (particle.pdgCode() == kProtonBar) { // Neg protons
          rMCRegistry.fill(HIST("hGenProtonN"), particle.pt());
        }
      }
    }
    // End of the Monte-Carlo generated QA
  }

  void processMCRec(MyFilteredMCRecCollision const& collision, MyFilteredMCTracks const& tracks, aod::McCollisions const&, aod::McParticles const&)
  {

    if (!collision.has_mcCollision()) {
      return;
    }

    // Start of the Monte-Carlo reconstructed QA
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }

      if (trackFilters(track)) {
        auto particle = track.mcParticle();
        if (particle.isPhysicalPrimary()) {

          // Efficiency - Reconstructed
          rMCRegistry.fill(HIST("hReconstructed"), track.pt());
          if (particle.pdgCode() == kPiPlus) { // Pos pions
            rMCRegistry.fill(HIST("hRecPionP"), track.pt());
          } else if (particle.pdgCode() == kPiMinus) { // Neg pions
            rMCRegistry.fill(HIST("hRecPionN"), track.pt());
          } else if (particle.pdgCode() == kKPlus) { // Pos kaons
            rMCRegistry.fill(HIST("hRecKaonP"), track.pt());
          } else if (particle.pdgCode() == kKMinus) { // Neg kaons
            rMCRegistry.fill(HIST("hRecKaonN"), track.pt());
          } else if (particle.pdgCode() == kProton) { // Pos protons
            rMCRegistry.fill(HIST("hRecProtonP"), track.pt());
          } else if (particle.pdgCode() == kProtonBar) { // Neg protons
            rMCRegistry.fill(HIST("hRecProtonN"), track.pt());
          }

          // Purity (PID)
          assocPID = trackPID(track);

          if (track.sign() > 0) {     // Positive tracks
            if (assocPID[0] == 0.0) { // Pions
              rMCRegistry.fill(HIST("hSelectPionP"), track.pt());
              if (particle.pdgCode() == kPiPlus) {
                rMCRegistry.fill(HIST("hTrueSelectPionP"), track.pt());
              }
            } else if (assocPID[0] == 1.0) { // Kaons
              rMCRegistry.fill(HIST("hSelectKaonP"), track.pt());
              if (particle.pdgCode() == kKPlus) {
                rMCRegistry.fill(HIST("hTrueSelectKaonP"), track.pt());
              }
            } else if (assocPID[0] == 2.0) { // Protons
              rMCRegistry.fill(HIST("hSelectProtonP"), track.pt());
              if (particle.pdgCode() == kProton) {
                rMCRegistry.fill(HIST("hTrueSelectProtonP"), track.pt());
              }
            }
          } else if (track.sign() < 0) { // Negative tracks
            if (assocPID[0] == 0.0) {    // Pions
              rMCRegistry.fill(HIST("hSelectPionN"), track.pt());
              if (particle.pdgCode() == kPiMinus) {
                rMCRegistry.fill(HIST("hTrueSelectPionN"), track.pt());
              }
            } else if (assocPID[0] == 1.0) { // Kaons
              rMCRegistry.fill(HIST("hSelectKaonN"), track.pt());
              if (particle.pdgCode() == kKMinus) {
                rMCRegistry.fill(HIST("hTrueSelectKaonN"), track.pt());
              }
            } else if (assocPID[0] == 2.0) { // Protons
              rMCRegistry.fill(HIST("hSelectProtonN"), track.pt());
              if (particle.pdgCode() == kProtonBar) {
                rMCRegistry.fill(HIST("hTrueSelectProtonN"), track.pt());
              }
            }
          }
        }
      }
    }
    // End of the Monte-Carlo reconstructed QA
  }

  PROCESS_SWITCH(ThreeParticleCorrelations, processSame, "Process same-event correlations", true);
  PROCESS_SWITCH(ThreeParticleCorrelations, processMixed, "Process mixed-event correlations", true);
  PROCESS_SWITCH(ThreeParticleCorrelations, processMCSame, "Process MC same-event correlations", true);
  PROCESS_SWITCH(ThreeParticleCorrelations, processMCMixed, "Process MC mixed-event correlations", true);
  PROCESS_SWITCH(ThreeParticleCorrelations, processMCGen, "Process Monte-Carlo, generator level", false);
  PROCESS_SWITCH(ThreeParticleCorrelations, processMCRec, "Process Monte-Carlo, reconstructed level", false);

  //==========================================================================================================================================================================================================================================================================

  double getMagneticField(uint64_t timestamp)
  {
    static parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }

    return 0.1 * (grpo->getNominalL3Field()); // 1 T = 10 kG
  }

  double trackEff(TH1D** efficiencies, int sign, double pT)
  {

    int index = -999;
    if (sign > 0) {
      index = 0;
    } else if (sign < 0) {
      index = 1;
    }

    double efficiency = efficiencies[index]->GetBinContent(efficiencies[index]->FindBin(pT));
    if (efficiency > 0) {
      return efficiency;
    } else {
      return 1.0;
    }
  }

  template <class V0Cand>
  int v0Sign(const V0Cand& v0)
  {

    if (std::abs(v0.mLambda() - massLambda) <= std::abs(v0.mAntiLambda() - massLambda)) {
      return 1;
    } else if (std::abs(v0.mLambda() - massLambda) > std::abs(v0.mAntiLambda() - massLambda)) {
      return -1;
    }

    return 0;
  }

  template <class TrackCand>
  double* trackPID(const TrackCand& track)
  {

    static double pid[2]; // {PID, NSigma}

    double nSigma[3];
    double nSigmaTOF[3];
    nSigmaTOF[0] = track.tofNSigmaPi();
    nSigmaTOF[1] = track.tofNSigmaKa();
    nSigmaTOF[2] = track.tofNSigmaPr();

    nSigma[0] = std::abs(nSigmaTOF[0]);
    nSigma[1] = std::abs(nSigmaTOF[1]);
    nSigma[2] = std::abs(nSigmaTOF[2]);

    if (nSigma[0] <= std::min(nSigma[1], nSigma[2])) { // Pions
      pid[0] = 0.0;
      pid[1] = nSigmaTOF[0];
    } else if (nSigma[1] <= std::min(nSigma[0], nSigma[2])) { // Kaons
      pid[0] = 1.0;
      pid[1] = nSigmaTOF[1];
    } else if (nSigma[2] < std::min(nSigma[0], nSigma[1])) { // Protons
      pid[0] = 2.0;
      pid[1] = nSigmaTOF[2];
    }

    return pid;
  }

  //==========================================================================================================================================================================================================================================================================

  template <class V0Cand>
  bool v0Filters(const V0Cand& v0)
  {

    if (v0Sign(v0) == 1) {
      const auto& posDaughter = v0.template posTrack_as<MyFilteredTracks>();
      if (std::abs(posDaughter.tpcNSigmaPr()) > 4.0) {
        return kFALSE;
      }
    } else if (v0Sign(v0) == -1) {
      const auto& negDaughter = v0.template negTrack_as<MyFilteredTracks>();
      if (std::abs(negDaughter.tpcNSigmaPr()) > 4.0) {
        return kFALSE;
      }
    }

    return kTRUE;
  }

  template <class TrackCand>
  bool trackFilters(const TrackCand& track)
  {

    if (!track.hasTOF()) {
      return kFALSE;
    }

    if (trackPID(track)[0] == 0.0) { // Pions
      if (std::abs(track.tpcNSigmaPi()) > 4.0) {
        return kFALSE;
      }
      if (track.pt() < 0.3) {
        return kFALSE;
      } else if (track.pt() > 0.3 && track.pt() < 1.5) {
        if (std::abs(track.tofNSigmaPi()) > 4.0) {
          return kFALSE;
        }
      } else if (track.pt() > 1.5 && track.pt() < 2.3) {
        if (track.tofNSigmaPi() < -4.0 || track.tofNSigmaPi() > 0.0) {
          return kFALSE;
        }
      } else if (track.pt() > 2.3) {
        return kFALSE;
      }

    } else if (trackPID(track)[0] == 1.0) { // Kaons
      if (std::abs(track.tpcNSigmaKa()) > 4.0) {
        return kFALSE;
      }
      if (track.pt() < 0.5) {
        return kFALSE;
      } else if (track.pt() > 0.5 && track.pt() < 1.5) {
        if (std::abs(track.tofNSigmaKa()) > 4.0) {
          return kFALSE;
        }
      } else if (track.pt() > 1.5 && track.pt() < 2.0) {
        if (track.tofNSigmaKa() < -2.0 || track.tofNSigmaKa() > 4.0) {
          return kFALSE;
        }
      } else if (track.pt() > 2.0 && track.pt() < 2.5) {
        if (track.tofNSigmaKa() < 0.0 || track.tofNSigmaKa() > 4.0) {
          return kFALSE;
        }
      } else if (track.pt() > 2.5) {
        return kFALSE;
      }

    } else if (trackPID(track)[0] == 2.0) { // Protons
      if (std::abs(track.tpcNSigmaPr()) > 4.0) {
        return kFALSE;
      }
      if (track.pt() < 0.5) {
        return kFALSE;
      } else if (track.pt() > 0.5 && track.pt() < 0.7) {
        if (track.tofNSigmaPr() < -2.0 || track.tofNSigmaPr() > 4.0) {
          return kFALSE;
        }
      } else if (track.pt() > 0.7 && track.pt() < 2.5) {
        if (std::abs(track.tofNSigmaPr()) > 4.0) {
          return kFALSE;
        }
      } else if (track.pt() > 2.5) {
        if (track.tofNSigmaPr() < -2.0 || track.tofNSigmaPr() > 4.0) {
          return kFALSE;
        }
      }
    }

    return kTRUE;
  }

  template <class V0Cand, class TrackCand>
  bool correlationFilters(const V0Cand& v0, const TrackCand& track)
  {

    if (track.globalIndex() == v0.posTrackId() || track.globalIndex() == v0.negTrackId()) {
      return kFALSE;
    }

    return kTRUE;
  }

  template <class V0Cand, class TrackCand>
  bool fakeV0Filter(const V0Cand& v0, const TrackCand& track)
  {

    if (confFilterSwitch) {

      if (trackPID(track)[0] == 1.0) { // Kaons
        return kTRUE;
      }

      std::array<float, 2> massArray;
      std::array<float, 3> dMomArray;
      std::array<float, 3> aMomArray = track.pVector();
      if (trackPID(track)[0] == 0.0) {
        massArray = {constants::physics::MassProton, constants::physics::MassPionCharged};

        if (v0Sign(v0) == 1 && track.sign() == -1) { // Lambda - Pi_min
          const auto& dTrack = v0.template posTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        } else if (v0Sign(v0) == -1 && track.sign() == 1) { // Antilambda - Pi_plus
          const auto& dTrack = v0.template negTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        }
      } else if (trackPID(track)[0] == 2.0) {
        massArray = {constants::physics::MassPionCharged, constants::physics::MassProton};

        if (v0Sign(v0) == 1 && track.sign() == 1) { // Lambda - Proton
          const auto& dTrack = v0.template negTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        } else if (v0Sign(v0) == -1 && track.sign() == -1) { // Antilambda - Antiproton
          const auto& dTrack = v0.template posTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        }
      }

      double invMass = RecoDecay::m(std::array{dMomArray, aMomArray}, massArray);
      if (invMass >= massLambda - 4 * dGaussSigma && invMass <= massLambda + 4 * dGaussSigma) {
        return kFALSE;
      }
    }

    return kTRUE;
  }

  template <class V0Cand, class TrackCand>
  bool radialDistanceFilter(const V0Cand& v0, const TrackCand& track, double B, bool Mix)
  {

    auto proton = v0.template posTrack_as<MyFilteredTracks>();
    if (v0Sign(v0) == -1) {
      proton = v0.template negTrack_as<MyFilteredTracks>();
    }

    double dEta = proton.eta() - track.eta();
    if (std::abs(dEta) > 0.02) {
      return kTRUE;
    }

    double dPhiStar;
    double dPhi = proton.phi() - track.phi();
    double phaseProton = (-0.3 * B * proton.sign()) / (2 * proton.pt());
    double phaseTrack = (-0.3 * B * track.sign()) / (2 * track.pt());

    for (double r = 0.8; r <= 2.5; r += 0.01) {
      dPhiStar = RecoDecay::constrainAngle(dPhi + std::asin(phaseProton * r) - std::asin(phaseTrack * r), -constants::math::PIHalf);

      if (r == 0.8) {
        if (!Mix) {
          rPhiStarRegistry.fill(HIST("hSEProtonPreCut"), dPhiStar, dEta);
        } else {
          rPhiStarRegistry.fill(HIST("hMEProtonPreCut"), dPhiStar, dEta);
        }
      }

      if (std::abs(dPhiStar) < 0.1) {
        return kFALSE;
      }

      if (r == 0.8) {
        if (!Mix) {
          rPhiStarRegistry.fill(HIST("hSEProtonPostCut"), dPhiStar, dEta);
        } else {
          rPhiStarRegistry.fill(HIST("hMEProtonPostCut"), dPhiStar, dEta);
        }
      }
    }

    return kTRUE;
  }
};

//============================================================================================================================================================================================================================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<ThreeParticleCorrelations>(cfgc)};
  return workflow;
}

//============================================================================================================================================================================================================================================================================
