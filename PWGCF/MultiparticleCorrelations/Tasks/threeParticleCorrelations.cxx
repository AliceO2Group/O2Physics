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
#include "Common/DataModel/McCollisionExtra.h"
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

  // Analysis parameters
  float centMin = 0.0, centMax = 90.0;
  float zvtxMax = 7.0;
  float v0PtMin = 0.6, v0PtMax = 12.0;
  float v0EtaMax = 0.72;
  float trackPtMin = 0.2, trackPtMax = 3.0;
  float trackEtaMax = 0.8;

  double pionID = 0.0, kaonID = 1.0, protonID = 2.0;
  float nSigma0 = 0.0, nSigma2 = 2.0, nSigma4 = 4.0;

  float pionPtMin = 0.3, pionPtMax = 2.3, kaonPtMin = 0.5, kaonPtMax = 2.5, protonPtMin = 0.5, protonPtMax = 2.5;
  float pionPtMid = 1.5, kaonPtMid1 = 1.5, kaonPtMid2 = 2.0, protonPtMid = 0.7;

  float dEtaMin = 0.02, dPhiStarMin = 0.1;
  float rMin = 0.8, rMax = 2.5;

  // Particle masses
  double massLambda = constants::physics::MassLambda0;
  double dGaussSigma = 0.0021;

  // Histogram registry
  HistogramRegistry rMECorrRegistry{"MECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rSECorrRegistry{"SECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rMCRegistry{"MCRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rPhiStarRegistry{"PhiStarRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rQARegistry{"QARegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Collision & Event filters
  Filter collCent = aod::cent::centFT0C > centMin && aod::cent::centFT0C < centMax;
  Filter collZvtx = nabs(aod::collision::posZ) < zvtxMax;
  Filter mcCollZvtx = nabs(aod::mccollision::posZ) < zvtxMax;
  Filter evSelect = aod::evsel::sel8 == true;

  // V0 filters
  Filter v0Pt = aod::v0data::pt > v0PtMin && aod::v0data::pt < v0PtMax;
  Filter v0Eta = nabs(aod::v0data::eta) < v0EtaMax;

  // Track filters
  Filter trackPt = aod::track::pt > trackPtMin && aod::track::pt < trackPtMax;
  Filter trackEta = nabs(aod::track::eta) < trackEtaMax;
  Filter globalTracks = requireGlobalTrackInFilter();

  // Particle filters
  Filter particleEta = nabs(aod::mcparticle::eta) < trackEtaMax;

  // Table aliases - Data
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;
  using MyFilteredV0s = soa::Filtered<aod::V0Datas>;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                   aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;

  // Table aliases - MC
  using MyFilteredMCGenCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::McCollsExtra>>;
  using MyFilteredMCGenCollision = MyFilteredMCGenCollisions::iterator;
  using MyFilteredMCParticles = soa::Filtered<aod::McParticles>;
  using MyFilteredMCRecCollision = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>>::iterator;
  using MyFilteredMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels,
                                                     aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                     aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;

  // Partitions
  Partition<MyFilteredMCParticles> mcParticles = aod::mcparticle::pt > trackPtMin && aod::mcparticle::pt < trackPtMax;
  Partition<MyFilteredMCParticles> mcTriggers = ((aod::mcparticle::pdgCode == static_cast<int>(kLambda0) || aod::mcparticle::pdgCode == static_cast<int>(kLambda0Bar)) &&
						 aod::mcparticle::pt > v0PtMin && aod::mcparticle::pt < v0PtMax && nabs(aod::mcparticle::eta) < v0EtaMax);
  Partition<MyFilteredMCParticles> mcAssociates = (((aod::mcparticle::pdgCode == static_cast<int>(kPiPlus) || aod::mcparticle::pdgCode == static_cast<int>(kPiMinus)) && aod::mcparticle::pt > pionPtMin && aod::mcparticle::pt < pionPtMax) ||
						   ((aod::mcparticle::pdgCode == static_cast<int>(kKPlus) || aod::mcparticle::pdgCode == static_cast<int>(kKMinus)) && aod::mcparticle::pt > kaonPtMin && aod::mcparticle::pt < kaonPtMax) ||
						   ((aod::mcparticle::pdgCode == static_cast<int>(kProton) || aod::mcparticle::pdgCode == static_cast<int>(kProtonBar)) && aod::mcparticle::pt > protonPtMin));
  
  // Mixed-events binning policy
  SliceCache cache;
  Preslice<MyFilteredMCParticles> perCol = aod::mcparticle::mcCollisionId;
  
  ConfigurableAxis confCentBins{"confCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "ME Centrality binning"};
  ConfigurableAxis confZvtxBins{"confZvtxBins", {VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 0.0f, 1.0f, 3.0f, 5.0f, 7.0f}, "ME Zvtx binning"};
  using BinningType = ColumnBinningPolicy<aod::cent::CentFT0C, aod::collision::PosZ>;
  using BinningTypeMC = ColumnBinningPolicy<aod::mccollisionprop::BestCollisionCentFT0C, aod::mccollision::PosZ>;

  BinningType collBinning{{confCentBins, confZvtxBins}, true};
  BinningTypeMC collBinningMC{{confCentBins, confZvtxBins}, true};
  Pair<MyFilteredCollisions, MyFilteredV0s, MyFilteredTracks, BinningType> pairData{collBinning, 5, -1, &cache};
  SameKindPair<MyFilteredMCGenCollisions, MyFilteredMCParticles, BinningTypeMC> pairMC{collBinningMC, 5, -1, &cache};

  // Process configurables
  Configurable<bool> confFilterSwitch{"confFilterSwitch", false, "Switch for the fakeV0Filter function"};

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
    const AxisSpec trackEtaAxis{32, -0.8, 0.8};
    const AxisSpec lambdaInvMassAxis{100, 1.08, 1.16};

    // QA & PID
    rQARegistry.add("hTrackPt", "hTrackPt", {HistType::kTH1D, {{100, 0, 4}}});
    rQARegistry.add("hTrackEta", "hTrackEta", {HistType::kTH1D, {{100, -1, 1}}});
    rQARegistry.add("hTrackPhi", "hTrackPhi", {HistType::kTH1D, {{100, (-1. / 2) * constants::math::PI, (5. / 2) * constants::math::PI}}});
    rQARegistry.add("hEventCentrality", "hEventCentrality", {HistType::kTH1D, {{centralityAxis}}});
    rQARegistry.add("hEventCentrality_MC", "hEventCentrality_MC", {HistType::kTH1D, {{centralityAxis}}});
    rQARegistry.add("hEventZvtx", "hEventZvtx", {HistType::kTH1D, {{zvtxAxis}}});

    rQARegistry.add("hPtPion", "hPtPion", {HistType::kTH1D, {{trackPtAxis}}});
    rQARegistry.add("hPtKaon", "hPtKaon", {HistType::kTH1D, {{trackPtAxis}}});
    rQARegistry.add("hPtProton", "hPtProton", {HistType::kTH1D, {{trackPtAxis}}});
    rQARegistry.add("hPtV0", "hPtV0", {HistType::kTH1D, {{v0PtAxis}}});      
    rQARegistry.add("hdEdx", "hdEdx", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxPion", "hdEdxPion", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxKaon", "hdEdxKaon", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxProton", "hdEdxProton", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hBeta", "hBeta", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaPion", "hBetaPion", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaKaon", "hBetaKaon", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaProton", "hBetaProton", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});

    rQARegistry.add("hTPCPion", "hTPCPion", {HistType::kTH2D, {{trackPtAxis}, {241, -6, 6}}});
    rQARegistry.add("hTPCKaon", "hTPCKaon", {HistType::kTH2D, {{trackPtAxis}, {241, -6, 6}}});
    rQARegistry.add("hTPCProton", "hTPCProton", {HistType::kTH2D, {{trackPtAxis}, {241, -6, 6}}});
    rQARegistry.add("hTOFPion", "hTOFPion", {HistType::kTH2D, {{trackPtAxis}, {1000, -50, 50}}});
    rQARegistry.add("hTOFKaon", "hTOFKaon", {HistType::kTH2D, {{trackPtAxis}, {1000, -50, 50}}});
    rQARegistry.add("hTOFProton", "hTOFProton", {HistType::kTH2D, {{trackPtAxis}, {1000, -50, 50}}});

    rQARegistry.add("hInvMassLambda", "hInvMassLambda", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hInvMassAntiLambda", "hInvMassAntiLambda", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hNLambdas", "hNLambdas", {HistType::kTH3D, {{2, -2, 2}, {v0PtAxis}, {centralityAxis}}});

    // PhiStar
    rPhiStarRegistry.add("hSEProtonPreCut", "hSEProtonPreCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});
    rPhiStarRegistry.add("hSEProtonPostCut", "hSEProtonPostCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});
    rPhiStarRegistry.add("hMEProtonPreCut", "hMEProtonPreCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});
    rPhiStarRegistry.add("hMEProtonPostCut", "hMEProtonPostCut", {HistType::kTH2D, {{80, -0.2, 0.2}, {40, -0.1, 0.1}}});

    // Efficiency
    rMCRegistry.add("hGenerated", "hGenerated", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hGenPionP", "hGenPionP", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hGenPionN", "hGenPionN", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hGenKaonP", "hGenKaonP", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hGenKaonN", "hGenKaonN", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hGenProtonP", "hGenProtonP", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hGenProtonN", "hGenProtonN", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hReconstructed", "hReconstructed", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hRecPionP", "hRecPionP", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hRecPionN", "hRecPionN", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hRecKaonP", "hRecKaonP", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hRecKaonN", "hRecKaonN", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hRecProtonP", "hRecProtonP", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});
    rMCRegistry.add("hRecProtonN", "hRecProtonN", {HistType::kTH2D, {{trackPtAxis}, {trackEtaAxis}}});

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
    rSECorrRegistry.add("hSameLambdaPion_SGNL", "Same-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_SB", "Same-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_SGNL", "Same-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_SB", "Same-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_SGNL", "Same-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_SB", "Same-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_MC", "Same-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_MC", "Same-event #Lambda - K correlator (MC)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_MC", "Same-event #Lambda - p correlator (MC)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    rMECorrRegistry.add("hMixLambdaPion_SGNL", "Mixed-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_SB", "Mixed-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_SGNL", "Mixed-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_SB", "Mixed-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_SGNL", "Mixed-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_SB", "Mixed-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_MC", "Mixed-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_MC", "Mixed-event #Lambda - K correlator (MC)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_MC", "Mixed-event #Lambda - p correlator (MC)", {HistType::kTHnSparseF, {{phiAxis}, {etaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

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
    rQARegistry.fill(HIST("hEventCentrality"), collision.centFT0C());
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
        if (assocPID[0] == pionID) { // Pions
	  rQARegistry.fill(HIST("hPtPion"), track.pt(), 1. / trackEff(hEffPions, track.sign(), track.pt()));
          rQARegistry.fill(HIST("hdEdxPion"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaPion"), track.pt(), track.beta());
        } else if (assocPID[0] == kaonID) { // Kaons
	  rQARegistry.fill(HIST("hPtKaon"), track.pt(), 1. / trackEff(hEffKaons, track.sign(), track.pt()));
          rQARegistry.fill(HIST("hdEdxKaon"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaKaon"), track.pt(), track.beta());
        } else if (assocPID[0] == protonID) { // Protons
	  rQARegistry.fill(HIST("hPtProton"), track.pt(), 1. / trackEff(hEffProtons, track.sign(), track.pt()));
          rQARegistry.fill(HIST("hdEdxProton"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaProton"), track.pt(), track.beta());
        }
      }
    }
    // End of the Track QA

    // Start of the Same-Event correlations
    for (const auto& trigger : v0s) {
      if (v0Filters(trigger)) {

	rQARegistry.fill(HIST("hPtV0"), trigger.pt());
        triggSign = v0Sign(trigger);
        if (triggSign == 1) {
          candMass = trigger.mLambda();
          rQARegistry.fill(HIST("hInvMassLambda"), trigger.mLambda(), trigger.pt(), collision.centFT0C());
        } else if (triggSign == -1) {
          candMass = trigger.mAntiLambda();
          rQARegistry.fill(HIST("hInvMassAntiLambda"), trigger.mAntiLambda(), trigger.pt(), collision.centFT0C());
        }

        for (const auto& associate : tracks) {
          if (trackFilters(associate)) {
            if (correlationFilters(trigger, associate) && radialDistanceFilter(trigger, associate, bField, false) && fakeV0Filter(trigger, associate)) {

              assocPID = trackPID(associate);
              deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
              deltaEta = trigger.eta() - associate.eta();

              if (candMass >= massLambda - 4 * dGaussSigma && candMass <= massLambda + 4 * dGaussSigma) {
                if (assocPID[0] == pionID) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
                } else if (assocPID[0] == kaonID) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
                } else if (assocPID[0] == protonID) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
                }
              } else if (candMass >= massLambda - 8 * dGaussSigma && candMass <= massLambda + 8 * dGaussSigma) {
                if (assocPID[0] == pionID) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
                } else if (assocPID[0] == kaonID) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
                } else if (assocPID[0] == protonID) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
                }
              }
            }
          }
        }
      }
    }
    // End of the Same-Event correlations
  }

  void processMixed(MyFilteredCollisions const&, MyFilteredV0s const&, MyFilteredTracks const&, aod::BCsWithTimestamps const&)
  {    

    // Start of the Mixed-Event correlations
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
              if (assocPID[0] == pionID) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
              } else if (assocPID[0] == kaonID) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
              } else if (assocPID[0] == protonID) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
              }
            } else if (candMass >= massLambda - 8 * dGaussSigma && candMass <= massLambda + 8 * dGaussSigma) {
              if (assocPID[0] == pionID) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate.sign(), associate.pt()));
              } else if (assocPID[0] == kaonID) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate.sign(), associate.pt()));
              } else if (assocPID[0] == protonID) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate.sign(), associate.pt()));
              }
            }
          }
        }
      }
    }
    // End of the Mixed-Event Correlations
  }

  void processMCSame(MyFilteredMCGenCollision const& collision, MyFilteredMCParticles const&)
  {

    rQARegistry.fill(HIST("hEventCentrality_MC"), collision.bestCollisionCentFT0C());
    auto groupMCTriggers = mcTriggers->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto groupMCAssociates = mcAssociates->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);

    // Start of the MC Same-Event correlations
    for (const auto& trigger : groupMCTriggers) {
      if (trigger.isPhysicalPrimary()) {

        if (trigger.pdgCode() > 0) {
          triggSign = 1;
        } else if (trigger.pdgCode() < 0) {
          triggSign = -1;
        }
        rQARegistry.fill(HIST("hNLambdas"), triggSign, trigger.pt(), collision.bestCollisionCentFT0C());

        for (const auto& associate : groupMCAssociates) {
          if (associate.isPhysicalPrimary()) {

            if (associate.pdgCode() > 0) {
              assocSign = 1;
            } else if (associate.pdgCode() < 0) {
              assocSign = -1;
            }

            deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
            deltaEta = trigger.eta() - associate.eta();

            if (std::abs(associate.pdgCode()) == kPiPlus) {
              rSECorrRegistry.fill(HIST("hSameLambdaPion_MC"), deltaPhi, deltaEta, collision.bestCollisionCentFT0C(), collision.posZ(), triggSign, assocSign);
            } else if (std::abs(associate.pdgCode()) == kKPlus) {
              rSECorrRegistry.fill(HIST("hSameLambdaKaon_MC"), deltaPhi, deltaEta, collision.bestCollisionCentFT0C(), collision.posZ(), triggSign, assocSign);
            } else if (std::abs(associate.pdgCode()) == kProton) {
              rSECorrRegistry.fill(HIST("hSameLambdaProton_MC"), deltaPhi, deltaEta, collision.bestCollisionCentFT0C(), collision.posZ(), triggSign, assocSign);
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
    for (const auto& [coll_1, v0_1, coll_2, particle_2] : pairMC) {
      auto groupMCTriggers = mcTriggers->sliceByCached(aod::mcparticle::mcCollisionId, coll_1.globalIndex(), cache);
      auto groupMCAssociates = mcAssociates->sliceByCached(aod::mcparticle::mcCollisionId, coll_2.globalIndex(), cache);

      for (const auto& [trigger, associate] : soa::combinations(soa::CombinationsFullIndexPolicy(groupMCTriggers, groupMCAssociates))) {
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
            rMECorrRegistry.fill(HIST("hMixLambdaPion_MC"), deltaPhi, deltaEta, coll_1.bestCollisionCentFT0C(), coll_1.posZ(), triggSign, assocSign);
          } else if (std::abs(associate.pdgCode()) == kKPlus) {
            rMECorrRegistry.fill(HIST("hMixLambdaKaon_MC"), deltaPhi, deltaEta, coll_1.bestCollisionCentFT0C(), coll_1.posZ(), triggSign, assocSign);
          } else if (std::abs(associate.pdgCode()) == kProton) {
            rMECorrRegistry.fill(HIST("hMixLambdaProton_MC"), deltaPhi, deltaEta, coll_1.bestCollisionCentFT0C(), coll_1.posZ(), triggSign, assocSign);
          }
        }
      }
    }
    // End of the MC Mixed-events Correlations
  }

  void processMCGen(MyFilteredMCGenCollision const&, MyFilteredMCParticles const&)
  {

    // Start of the Monte-Carlo generated QA
    for (const auto& particle : mcParticles) {
      if (particle.isPhysicalPrimary()) {

        // Efficiency - Generated
        rMCRegistry.fill(HIST("hGenerated"), particle.pt(), particle.eta());
        if (particle.pdgCode() == kPiPlus) { // Pos pions
          rMCRegistry.fill(HIST("hGenPionP"), particle.pt(), particle.eta());
        } else if (particle.pdgCode() == kPiMinus) { // Neg pions
          rMCRegistry.fill(HIST("hGenPionN"), particle.pt(), particle.eta());
        } else if (particle.pdgCode() == kKPlus) { // Pos kaons
          rMCRegistry.fill(HIST("hGenKaonP"), particle.pt(), particle.eta());
        } else if (particle.pdgCode() == kKMinus) { // Neg kaons
          rMCRegistry.fill(HIST("hGenKaonN"), particle.pt(), particle.eta());
        } else if (particle.pdgCode() == kProton) { // Pos protons
          rMCRegistry.fill(HIST("hGenProtonP"), particle.pt(), particle.eta());
        } else if (particle.pdgCode() == kProtonBar) { // Neg protons
          rMCRegistry.fill(HIST("hGenProtonN"), particle.pt(), particle.eta());
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
          rMCRegistry.fill(HIST("hReconstructed"), track.pt(), track.eta());
          if (particle.pdgCode() == kPiPlus) { // Pos pions
            rMCRegistry.fill(HIST("hRecPionP"), track.pt(), track.eta());
          } else if (particle.pdgCode() == kPiMinus) { // Neg pions
            rMCRegistry.fill(HIST("hRecPionN"), track.pt(), track.eta());
          } else if (particle.pdgCode() == kKPlus) { // Pos kaons
            rMCRegistry.fill(HIST("hRecKaonP"), track.pt(), track.eta());
          } else if (particle.pdgCode() == kKMinus) { // Neg kaons
            rMCRegistry.fill(HIST("hRecKaonN"), track.pt(), track.eta());
          } else if (particle.pdgCode() == kProton) { // Pos protons
            rMCRegistry.fill(HIST("hRecProtonP"), track.pt(), track.eta());
          } else if (particle.pdgCode() == kProtonBar) { // Neg protons
            rMCRegistry.fill(HIST("hRecProtonN"), track.pt(), track.eta());
          }

          // Purity (PID)
          assocPID = trackPID(track);

          if (track.sign() > 0) {     // Positive tracks
            if (assocPID[0] == pionID) { // Pions
              rMCRegistry.fill(HIST("hSelectPionP"), track.pt());
              if (particle.pdgCode() == kPiPlus) {
                rMCRegistry.fill(HIST("hTrueSelectPionP"), track.pt());
              }
            } else if (assocPID[0] == kaonID) { // Kaons
              rMCRegistry.fill(HIST("hSelectKaonP"), track.pt());
              if (particle.pdgCode() == kKPlus) {
                rMCRegistry.fill(HIST("hTrueSelectKaonP"), track.pt());
              }
            } else if (assocPID[0] == protonID) { // Protons
              rMCRegistry.fill(HIST("hSelectProtonP"), track.pt());
              if (particle.pdgCode() == kProton) {
                rMCRegistry.fill(HIST("hTrueSelectProtonP"), track.pt());
              }
            }
          } else if (track.sign() < 0) { // Negative tracks
            if (assocPID[0] == pionID) {    // Pions
              rMCRegistry.fill(HIST("hSelectPionN"), track.pt());
              if (particle.pdgCode() == kPiMinus) {
                rMCRegistry.fill(HIST("hTrueSelectPionN"), track.pt());
              }
            } else if (assocPID[0] == kaonID) { // Kaons
              rMCRegistry.fill(HIST("hSelectKaonN"), track.pt());
              if (particle.pdgCode() == kKMinus) {
                rMCRegistry.fill(HIST("hTrueSelectKaonN"), track.pt());
              }
            } else if (assocPID[0] == protonID) { // Protons
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
  PROCESS_SWITCH(ThreeParticleCorrelations, processMCSame, "Process MC same-event correlations", false);
  PROCESS_SWITCH(ThreeParticleCorrelations, processMCMixed, "Process MC mixed-event correlations", false);
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
      pid[0] = pionID;
      pid[1] = nSigmaTOF[0];
    } else if (nSigma[1] <= std::min(nSigma[0], nSigma[2])) { // Kaons
      pid[0] = kaonID;
      pid[1] = nSigmaTOF[1];
    } else if (nSigma[2] < std::min(nSigma[0], nSigma[1])) { // Protons
      pid[0] = protonID;
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
      if (std::abs(posDaughter.tpcNSigmaPr()) > nSigma4) {
        return kFALSE;
      }
    } else if (v0Sign(v0) == -1) {
      const auto& negDaughter = v0.template negTrack_as<MyFilteredTracks>();
      if (std::abs(negDaughter.tpcNSigmaPr()) > nSigma4) {
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

    if (trackPID(track)[0] == pionID) { // Pions
      if (std::abs(track.tpcNSigmaPi()) > nSigma4) {
        return kFALSE;
      }
      if (track.pt() < pionPtMin) {
        return kFALSE;
      } else if (track.pt() > pionPtMin && track.pt() < pionPtMid) {
        if (std::abs(track.tofNSigmaPi()) > nSigma4) {
          return kFALSE;
        }
      } else if (track.pt() > pionPtMid && track.pt() < pionPtMax) {
        if (track.tofNSigmaPi() < -nSigma4 || track.tofNSigmaPi() > nSigma0) {
          return kFALSE;
        }
      } else if (track.pt() > pionPtMax) {
        return kFALSE;
      }

    } else if (trackPID(track)[0] == kaonID) { // Kaons
      if (std::abs(track.tpcNSigmaKa()) > nSigma4) {
        return kFALSE;
      }
      if (track.pt() < kaonPtMin) {
        return kFALSE;
      } else if (track.pt() > kaonPtMin && track.pt() < kaonPtMid1) {
        if (std::abs(track.tofNSigmaKa()) > nSigma4) {
          return kFALSE;
        }
      } else if (track.pt() > kaonPtMid1 && track.pt() < kaonPtMid2) {
        if (track.tofNSigmaKa() < -nSigma2 || track.tofNSigmaKa() > nSigma4) {
          return kFALSE;
        }
      } else if (track.pt() > kaonPtMid2 && track.pt() < kaonPtMax) {
        if (track.tofNSigmaKa() < nSigma0 || track.tofNSigmaKa() > nSigma4) {
          return kFALSE;
        }
      } else if (track.pt() > kaonPtMax) {
        return kFALSE;
      }

    } else if (trackPID(track)[0] == protonID) { // Protons
      if (std::abs(track.tpcNSigmaPr()) > nSigma4) {
        return kFALSE;
      }
      if (track.pt() < protonPtMin) {
        return kFALSE;
      } else if (track.pt() > protonPtMin && track.pt() < protonPtMid) {
        if (track.tofNSigmaPr() < -nSigma2 || track.tofNSigmaPr() > nSigma4) {
          return kFALSE;
        }
      } else if (track.pt() > protonPtMid && track.pt() < protonPtMax) {
        if (std::abs(track.tofNSigmaPr()) > nSigma4) {
          return kFALSE;
        }
      } else if (track.pt() > protonPtMax) {
        if (track.tofNSigmaPr() < -nSigma2 || track.tofNSigmaPr() > nSigma4) {
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

      if (trackPID(track)[0] == kaonID) { // Kaons
        return kTRUE;
      }

      std::array<float, 2> massArray;
      std::array<float, 3> dMomArray;
      std::array<float, 3> aMomArray = track.pVector();
      if (trackPID(track)[0] == pionID) {
        massArray = {constants::physics::MassProton, constants::physics::MassPionCharged};

        if (v0Sign(v0) == 1 && track.sign() == -1) { // Lambda - Pi_min
          const auto& dTrack = v0.template posTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        } else if (v0Sign(v0) == -1 && track.sign() == 1) { // Antilambda - Pi_plus
          const auto& dTrack = v0.template negTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        }
      } else if (trackPID(track)[0] == protonID) {
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
    if (std::abs(dEta) > dEtaMin) {
      return kTRUE;
    }

    double dPhiStar;
    double dPhi = proton.phi() - track.phi();
    double phaseProton = (-0.3 * B * proton.sign()) / (2 * proton.pt());
    double phaseTrack = (-0.3 * B * track.sign()) / (2 * track.pt());

    for (double r = rMin; r <= rMax; r += 0.01) {
      dPhiStar = RecoDecay::constrainAngle(dPhi + std::asin(phaseProton * r) - std::asin(phaseTrack * r), -constants::math::PIHalf);

      if (r == rMin) {
        if (!Mix) {
          rPhiStarRegistry.fill(HIST("hSEProtonPreCut"), dPhiStar, dEta);
        } else {
          rPhiStarRegistry.fill(HIST("hMEProtonPreCut"), dPhiStar, dEta);
        }
      }

      if (std::abs(dPhiStar) < dPhiStarMin) {
        return kFALSE;
      }

      if (r == rMin) {
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
