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
using namespace constants::physics;

struct ThreeParticleCorrelations {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Analysis parameters
  float centMin = 0.0, centMax = 90.0;
  float zvtxMax = 7.0;
  float v0PtMin = 0.6, v0PtMax = 12.0;
  float v0EtaMax = 0.72;
  float trackPtMin = 0.2, trackPtMax = 3.0;
  float trackEtaMax = 0.8;

  // Track PID parameters
  double pionID = 0.0, kaonID = 1.0, protonID = 2.0;
  float nSigma0 = 0.0, nSigma2 = 2.0, nSigma4 = 4.0, nSigma5 = 5.0;

  // V0 filter parameters
  float tpcNCrossedRowsMin = 70.0;
  float decayRMin = 1.2, ctauMax = 30.0;
  float cosPAMin = 0.995;
  float dcaProtonMin = 0.05, dcaPionMin = 0.2;
  int dcaV0DauMax = 1;
  
  // Track filter parameters
  float pionPtMin = 0.3, pionPtMax = 2.3, kaonPtMin = 0.5, kaonPtMax = 2.5, protonPtMin = 0.5, protonPtMax = 2.5;
  float pionPtMid = 1.5, kaonPtMid1 = 1.5, kaonPtMid2 = 2.0, protonPtMid = 0.7;

  // RD filter parameters
  float dEtaMax = 0.05, dEtaMin = 0.022;
  float dPhiStarMinOS = 0.075, dPhiStarMinSS = 0.12;
  float rMin = 0.8, rMax = 2.5;

  // Lambda invariant mass fit
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

  // Track filters
  Filter trackPt = aod::track::pt > trackPtMin && aod::track::pt < trackPtMax;
  Filter trackEta = nabs(aod::track::eta) < trackEtaMax;
  Filter globalTracks = requireGlobalTrackInFilter();

  // Particle filters
  Filter particleEta = nabs(aod::mcparticle::eta) < trackEtaMax;

  // Table aliases - Data
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                   aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;

  // Table aliases - MC Gen
  using MyFilteredMCGenCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::McCollsExtra>>;
  using MyFilteredMCGenCollision = MyFilteredMCGenCollisions::iterator;
  using MyFilteredMCParticles = soa::Filtered<aod::McParticles>;

  // Table aliases - MC Rec
  using MCRecCollisions = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::McCollisionLabels>;
  using MyFilteredMCRecCollisions = soa::Filtered<MCRecCollisions>;
  using MyMCV0s = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using MyFilteredMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels,
                                                     aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                     aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>>;

  // Partitions
  Partition<MyFilteredMCParticles> mcTracks = aod::mcparticle::pt > trackPtMin && aod::mcparticle::pt < trackPtMax;
  Partition<MyFilteredMCParticles> mcV0s = aod::mcparticle::pt > v0PtMin && aod::mcparticle::pt < v0PtMax && nabs(aod::mcparticle::eta) < v0EtaMax;
  Partition<MyFilteredMCParticles> mcTriggers = ((aod::mcparticle::pdgCode == static_cast<int>(kLambda0) || aod::mcparticle::pdgCode == static_cast<int>(kLambda0Bar)) &&
                                                 aod::mcparticle::pt > v0PtMin && aod::mcparticle::pt < v0PtMax && nabs(aod::mcparticle::eta) < v0EtaMax);
  Partition<MyFilteredMCParticles> mcAssociates = (((aod::mcparticle::pdgCode == static_cast<int>(kPiPlus) || aod::mcparticle::pdgCode == static_cast<int>(kPiMinus)) && aod::mcparticle::pt > pionPtMin && aod::mcparticle::pt < pionPtMax) ||
                                                   ((aod::mcparticle::pdgCode == static_cast<int>(kKPlus) || aod::mcparticle::pdgCode == static_cast<int>(kKMinus)) && aod::mcparticle::pt > kaonPtMin && aod::mcparticle::pt < kaonPtMax) ||
                                                   ((aod::mcparticle::pdgCode == static_cast<int>(kProton) || aod::mcparticle::pdgCode == static_cast<int>(kProtonBar)) && aod::mcparticle::pt > protonPtMin));

  // Mixed-events binning policy
  SliceCache cache;
  Preslice<MyFilteredMCParticles> perCol = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<aod::McCollisionLabels> perMCCol = aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis confCentBins{"confCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "ME Centrality binning"};
  ConfigurableAxis confZvtxBins{"confZvtxBins", {VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 0.0f, 1.0f, 3.0f, 5.0f, 7.0f}, "ME Zvtx binning"};
  using BinningType = ColumnBinningPolicy<aod::cent::CentFT0C, aod::collision::PosZ>;
  using BinningTypeMC = ColumnBinningPolicy<aod::mccollisionprop::BestCollisionCentFT0C, aod::mccollision::PosZ>;

  BinningType collBinning{{confCentBins, confZvtxBins}, true};
  BinningTypeMC collBinningMC{{confCentBins, confZvtxBins}, true};
  Pair<MyFilteredCollisions, aod::V0Datas, MyFilteredTracks, BinningType> pairData{collBinning, 5, -1, &cache};
  SameKindPair<MyFilteredMCGenCollisions, MyFilteredMCParticles, BinningTypeMC> pairMC{collBinningMC, 5, -1, &cache};

  // Process configurables
  Configurable<bool> confFakeV0Switch{"confFakeV0Switch", false, "Switch for the fakeV0Filter function"};
  Configurable<bool> confRDSwitch{"confRDSwitch", true, "Switch for the radialDistanceFilter function"};

  // Efficiency histograms
  TH3D** hEffPions = new TH3D*[2];
  TH3D** hEffKaons = new TH3D*[2];
  TH3D** hEffProtons = new TH3D*[2];

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
    const AxisSpec dPhiAxis{36, (-1. / 2) * constants::math::PI, (3. / 2) * constants::math::PI};
    const AxisSpec dEtaAxis{32, -1.52, 1.52};
    const AxisSpec v0PtAxis{114, 0.6, 12};
    const AxisSpec v0EtaAxis{36, -0.72, 0.72};
    const AxisSpec trackPtAxis{28, 0.2, 3};
    const AxisSpec trackEtaAxis{32, -0.8, 0.8};
    const AxisSpec lambdaInvMassAxis{100, 1.08, 1.16};

    // QA & PID
    rQARegistry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{3, 0, 3}}});
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "All");
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "kIsGoodZvtxFT0vsPV");
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "kNoSameBunchPileup");

    rQARegistry.add("hEventCentrality", "hEventCentrality", {HistType::kTH1D, {{centralityAxis}}});
    rQARegistry.add("hEventCentrality_MC", "hEventCentrality_MC", {HistType::kTH1D, {{centralityAxis}}});
    rQARegistry.add("hEventZvtx", "hEventZvtx", {HistType::kTH1D, {{zvtxAxis}}});
    rQARegistry.add("hTrackPt", "hTrackPt", {HistType::kTH1D, {{100, 0, 4}}});
    rQARegistry.add("hTrackEta", "hTrackEta", {HistType::kTH1D, {{100, -1, 1}}});
    rQARegistry.add("hTrackPhi", "hTrackPhi", {HistType::kTH1D, {{100, (-1. / 2) * constants::math::PI, (5. / 2) * constants::math::PI}}});

    rQARegistry.add("hPtPion", "hPtPion", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtKaon", "hPtKaon", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtProton", "hPtProton", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtV0", "hPtV0", {HistType::kTH2D, {{v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtPion_MC", "hPtPion_MC", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtKaon_MC", "hPtKaon_MC", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtProton_MC", "hPtProton_MC", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rQARegistry.add("hPtV0_MC", "hPtV0_MC", {HistType::kTH2D, {{v0PtAxis}, {centralityAxis}}});

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
    rQARegistry.add("hInvMassLambda_MC", "hInvMassLambda_MC", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hInvMassAntiLambda_MC", "hInvMassAntiLambda_MC", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});

    // PhiStar
    rPhiStarRegistry.add("hSEProtonPreCut_OS", "hSEProtonPreCut_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPreCut_SS", "hSEProtonPreCut_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPreCut_SSP", "hSEProtonPreCut_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPreCut_SSN", "hSEProtonPreCut_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPostCut_OS", "hSEProtonPostCut_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPostCut_SS", "hSEProtonPostCut_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPostCut_SSP", "hSEProtonPostCut_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEProtonPostCut_SSN", "hSEProtonPostCut_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_OS", "hSEPhiStarMean_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_SS", "hSEPhiStarMean_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_SSP", "hSEPhiStarMean_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_SSN", "hSEPhiStarMean_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});

    rPhiStarRegistry.add("hMEProtonPreCut_OS", "hMEProtonPreCut_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPreCut_SS", "hMEProtonPreCut_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPreCut_SSP", "hMEProtonPreCut_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPreCut_SSN", "hMEProtonPreCut_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPostCut_OS", "hMEProtonPostCut_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPostCut_SS", "hMEProtonPostCut_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPostCut_SSP", "hMEProtonPostCut_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEProtonPostCut_SSN", "hMEProtonPostCut_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_OS", "hMEPhiStarMean_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_SS", "hMEPhiStarMean_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_SSP", "hMEPhiStarMean_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_SSN", "hMEPhiStarMean_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});

    // Efficiency
    rMCRegistry.add("hGenerated", "hGenerated", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenPionP", "hGenPionP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenPionN", "hGenPionN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenKaonP", "hGenKaonP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenKaonN", "hGenKaonN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenProtonP", "hGenProtonP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenProtonN", "hGenProtonN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenLambdaP", "hGenLambdaP", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hGenLambdaN", "hGenLambdaN", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hReconstructed", "hReconstructed", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecPionP", "hRecPionP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecPionN", "hRecPionN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecKaonP", "hRecKaonP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecKaonN", "hRecKaonN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecProtonP", "hRecProtonP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecProtonN", "hRecProtonN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecLambdaP", "hRecLambdaP", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hRecLambdaN", "hRecLambdaN", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hIdentified", "hIdentified", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDPionP", "hPIDPionP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDPionN", "hPIDPionN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDKaonP", "hPIDKaonP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDKaonN", "hPIDKaonN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDProtonP", "hPIDProtonP", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDProtonN", "hPIDProtonN", {HistType::kTH3D, {{trackPtAxis}, {trackEtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDLambdaP_SGNL", "hPIDLambdaP_SGNL", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDLambdaP_SB", "hPIDLambdaP_SB", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDLambdaN_SGNL", "hPIDLambdaN_SGNL", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});
    rMCRegistry.add("hPIDLambdaN_SB", "hPIDLambdaN_SB", {HistType::kTH3D, {{v0PtAxis}, {v0EtaAxis}, {centralityAxis}}});

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
    rSECorrRegistry.add("hSameLambdaPion_SGNL", "Same-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_SB", "Same-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_SGNL", "Same-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_SB", "Same-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_SGNL", "Same-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_SB", "Same-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_MC", "Same-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_MC", "Same-event #Lambda - K correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_MC", "Same-event #Lambda - p correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    rSECorrRegistry.add("hSameLambdaPion_leftSB", "Same-event #Lambda - #pi correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaPion_rightSB", "Same-event #Lambda - #pi correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_leftSB", "Same-event #Lambda - K correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaKaon_rightSB", "Same-event #Lambda - K correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_leftSB", "Same-event #Lambda - p correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rSECorrRegistry.add("hSameLambdaProton_rightSB", "Same-event #Lambda - p correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    rMECorrRegistry.add("hMixLambdaPion_SGNL", "Mixed-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_SB", "Mixed-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_SGNL", "Mixed-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_SB", "Mixed-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_SGNL", "Mixed-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_SB", "Mixed-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_MC", "Mixed-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_MC", "Mixed-event #Lambda - K correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_MC", "Mixed-event #Lambda - p correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    rMECorrRegistry.add("hMixLambdaPion_leftSB", "Mixed-event #Lambda - #pi correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaPion_rightSB", "Mixed-event #Lambda - #pi correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_leftSB", "Mixed-event #Lambda - K correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaKaon_rightSB", "Mixed-event #Lambda - K correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_leftSB", "Mixed-event #Lambda - p correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    rMECorrRegistry.add("hMixLambdaProton_rightSB", "Mixed-event #Lambda - p correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    TList* efficiencyList = ccdb->getForTimeStamp<TList>("Users/j/jstaa/Efficiency/ChargedParticles", 1);
    hEffPions[0] = static_cast<TH3D*>(efficiencyList->FindObject("hEfficiencyPionP"));
    hEffPions[1] = static_cast<TH3D*>(efficiencyList->FindObject("hEfficiencyPionN"));
    hEffKaons[0] = static_cast<TH3D*>(efficiencyList->FindObject("hEfficiencyKaonP"));
    hEffKaons[1] = static_cast<TH3D*>(efficiencyList->FindObject("hEfficiencyKaonN"));
    hEffProtons[0] = static_cast<TH3D*>(efficiencyList->FindObject("hEfficiencyProtonP"));
    hEffProtons[1] = static_cast<TH3D*>(efficiencyList->FindObject("hEfficiencyProtonN"));
  }

  //==========================================================================================================================================================================================================================================================================

  void processSame(MyFilteredCollision const& collision, aod::V0Datas const& v0s, MyFilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    if (!acceptEvent(collision, true)) {
      return;
    }

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
          rQARegistry.fill(HIST("hPtPion"), track.pt(), collision.centFT0C(), 1. / trackEff(hEffPions, track, collision.centFT0C()));
          rQARegistry.fill(HIST("hdEdxPion"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaPion"), track.pt(), track.beta());
        } else if (assocPID[0] == kaonID) { // Kaons
          rQARegistry.fill(HIST("hPtKaon"), track.pt(), collision.centFT0C(), 1. / trackEff(hEffKaons, track, collision.centFT0C()));
          rQARegistry.fill(HIST("hdEdxKaon"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaKaon"), track.pt(), track.beta());
        } else if (assocPID[0] == protonID) { // Protons
          rQARegistry.fill(HIST("hPtProton"), track.pt(), collision.centFT0C(), 1. / trackEff(hEffProtons, track, collision.centFT0C()));
          rQARegistry.fill(HIST("hdEdxProton"), track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaProton"), track.pt(), track.beta());
        }
      }
    }
    // End of the Track QA

    // Start of the Same-Event correlations
    for (const auto& trigger : v0s) {
      if (v0Filters(collision, trigger, tracks)) {

        rQARegistry.fill(HIST("hPtV0"), trigger.pt(), collision.centFT0C());
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

              if (candMass >= MassLambda0 - 4 * dGaussSigma && candMass <= MassLambda0 + 4 * dGaussSigma) {
                if (assocPID[0] == pionID) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, collision.centFT0C()));
                } else if (assocPID[0] == kaonID) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, collision.centFT0C()));
                } else if (assocPID[0] == protonID) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, collision.centFT0C()));
                }

              } else if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                if (assocPID[0] == pionID) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, collision.centFT0C()));
                  if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass < MassLambda0 - 4 * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaPion_leftSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, collision.centFT0C()));
                  } else if (candMass > MassLambda0 + 4 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaPion_rightSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, collision.centFT0C()));
                  }

                } else if (assocPID[0] == kaonID) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, collision.centFT0C()));
                  if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass < MassLambda0 - 4 * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaKaon_leftSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, collision.centFT0C()));
                  } else if (candMass > MassLambda0 + 4 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaKaon_rightSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, collision.centFT0C()));
                  }

                } else if (assocPID[0] == protonID) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, collision.centFT0C()));
                  if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass < MassLambda0 - 4 * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaProton_leftSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, collision.centFT0C()));
                  } else if (candMass > MassLambda0 + 4 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaProton_rightSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, collision.centFT0C()));
                  }
                }
              }
            }
          }
        }
      }
    }
    // End of the Same-Event correlations
  }

  void processMixed(MyFilteredCollisions const&, aod::V0Datas const&, MyFilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    // Start of the Mixed-Event correlations
    for (const auto& [coll_1, v0_1, coll_2, track_2] : pairData) {
      if (!acceptEvent(coll_1, false) || !acceptEvent(coll_2, false)) {
        return;
      }

      auto bc = coll_1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      for (const auto& [trigger, associate] : soa::combinations(soa::CombinationsFullIndexPolicy(v0_1, track_2))) {
        if (v0Filters(coll_1, trigger, tracks) && trackFilters(associate)) {
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

            if (candMass >= MassLambda0 - 4 * dGaussSigma && candMass <= MassLambda0 + 4 * dGaussSigma) {
              if (assocPID[0] == pionID) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, coll_1.centFT0C()));
              } else if (assocPID[0] == kaonID) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, coll_1.centFT0C()));
              } else if (assocPID[0] == protonID) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, coll_1.centFT0C()));
              }

            } else if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
              if (assocPID[0] == pionID) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, coll_1.centFT0C()));
                if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass < MassLambda0 - 4 * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaPion_leftSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, coll_1.centFT0C()));
                } else if (candMass > MassLambda0 + 4 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaPion_rightSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffPions, associate, coll_1.centFT0C()));
                }

              } else if (assocPID[0] == kaonID) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, coll_1.centFT0C()));
                if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass < MassLambda0 - 4 * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaKaon_leftSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, coll_1.centFT0C()));
                } else if (candMass > MassLambda0 + 4 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaKaon_rightSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffKaons, associate, coll_1.centFT0C()));
                }

              } else if (assocPID[0] == protonID) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, coll_1.centFT0C()));
                if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass < MassLambda0 - 4 * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaProton_leftSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, coll_1.centFT0C()));
                } else if (candMass > MassLambda0 + 4 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaProton_rightSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), 1. / trackEff(hEffProtons, associate, coll_1.centFT0C()));
                }
              }
            }
          }
        }
      }
    }
    // End of the Mixed-Event Correlations
  }

  void processMCSame(MyFilteredMCGenCollision const& collision, MyFilteredMCParticles const&, soa::SmallGroups<MCRecCollisions> const& recCollisions)
  {

    if (recCollisions.size() == 1) {
      for (const auto& recCollision : recCollisions) {
        if (!acceptEvent(recCollision, false)) {
          return;
        }
      }
    }

    rQARegistry.fill(HIST("hEventCentrality_MC"), collision.bestCollisionCentFT0C());
    auto groupMCTriggers = mcTriggers->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto groupMCAssociates = mcAssociates->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);

    // Start of the MC Track QA
    for (const auto& track : groupMCAssociates) {
      if (track.isPhysicalPrimary()) {

        if (std::abs(track.pdgCode()) == kPiPlus) { // Pions
          rQARegistry.fill(HIST("hPtPion_MC"), track.pt(), collision.bestCollisionCentFT0C());
        } else if (std::abs(track.pdgCode()) == kKPlus) { // Kaons
          rQARegistry.fill(HIST("hPtKaon_MC"), track.pt(), collision.bestCollisionCentFT0C());
        } else if (std::abs(track.pdgCode()) == kProton) { // Protons
          rQARegistry.fill(HIST("hPtProton_MC"), track.pt(), collision.bestCollisionCentFT0C());
        }
      }
    }
    // End of the MC Track QA

    // Start of the MC Same-Event correlations
    for (const auto& trigger : groupMCTriggers) {
      if (trigger.isPhysicalPrimary()) {

        rQARegistry.fill(HIST("hPtV0_MC"), trigger.pt(), collision.bestCollisionCentFT0C());
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

  void processMCMixed(MyFilteredMCGenCollisions const&, MyFilteredMCParticles const&, MyFilteredMCRecCollisions const& recCollisions)
  {

    // Start of the MC Mixed-events Correlations
    for (const auto& [coll_1, v0_1, coll_2, particle_2] : pairMC) {
      auto recCollsA1 = recCollisions.sliceBy(perMCCol, coll_1.globalIndex());
      auto recCollsA2 = recCollisions.sliceBy(perMCCol, coll_2.globalIndex());
      if (recCollsA1.size() == 1 && recCollsA2.size() == 1) {
        for (const auto& recColl_1 : recCollsA1) {
          if (!acceptEvent(recColl_1, false)) {
            return;
          }
        }
        for (const auto& recColl_2 : recCollsA2) {
          if (!acceptEvent(recColl_2, false)) {
            return;
          }
        }
      }

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

  void processMCGen(MyFilteredMCGenCollision const& collision, MyFilteredMCParticles const&, soa::SmallGroups<MCRecCollisions> const& recCollisions)
  {

    if (recCollisions.size() == 1) {
      for (const auto& recCollision : recCollisions) {
        if (!acceptEvent(recCollision, false)) {
          return;
        }
      }
    }

    auto groupMCTracks = mcTracks->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto groupMCV0s = mcV0s->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);

    // Start of the Monte-Carlo generated QA
    for (const auto& particle : groupMCTracks) {
      if (particle.isPhysicalPrimary()) {

        // Track efficiency - Generated
        rMCRegistry.fill(HIST("hGenerated"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        if (particle.pdgCode() == kPiPlus) { // Pos pions
          rMCRegistry.fill(HIST("hGenPionP"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        } else if (particle.pdgCode() == kPiMinus) { // Neg pions
          rMCRegistry.fill(HIST("hGenPionN"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        } else if (particle.pdgCode() == kKPlus) { // Pos kaons
          rMCRegistry.fill(HIST("hGenKaonP"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        } else if (particle.pdgCode() == kKMinus) { // Neg kaons
          rMCRegistry.fill(HIST("hGenKaonN"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        } else if (particle.pdgCode() == kProton) { // Pos protons
          rMCRegistry.fill(HIST("hGenProtonP"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        } else if (particle.pdgCode() == kProtonBar) { // Neg protons
          rMCRegistry.fill(HIST("hGenProtonN"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        }
      }
    }

    for (const auto& particle : groupMCV0s) {
      if (particle.isPhysicalPrimary()) {

        // V0 efficiency - Generated
        if (particle.pdgCode() == kLambda0) { // Lambdas
          rMCRegistry.fill(HIST("hGenLambdaP"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        } else if (particle.pdgCode() == kLambda0Bar) { // AntiLambdas
          rMCRegistry.fill(HIST("hGenLambdaN"), particle.pt(), particle.eta(), collision.bestCollisionCentFT0C());
        }
      }
    }
    // End of the Monte-Carlo generated QA
  }

  void processMCRec(MyFilteredMCRecCollisions::iterator const& collision, MyMCV0s const& v0s, MyFilteredMCTracks const& tracks, aod::McCollisions const&, aod::McParticles const&)
  {

    if (!acceptEvent(collision, false) || !collision.has_mcCollision()) {
      return;
    }

    // Start of the Monte-Carlo reconstructed QA
    for (const auto& track : tracks) {

      if (!track.has_mcParticle()) {
        continue;
      }
      auto particle = track.mcParticle();

      // Track efficiency - Reconstructed
      rMCRegistry.fill(HIST("hReconstructed"), track.pt(), track.eta(), collision.centFT0C());
      if (particle.pdgCode() == kPiPlus) { // Pos pions
        rMCRegistry.fill(HIST("hRecPionP"), track.pt(), track.eta(), collision.centFT0C());
      } else if (particle.pdgCode() == kPiMinus) { // Neg pions
        rMCRegistry.fill(HIST("hRecPionN"), track.pt(), track.eta(), collision.centFT0C());
      } else if (particle.pdgCode() == kKPlus) { // Pos kaons
        rMCRegistry.fill(HIST("hRecKaonP"), track.pt(), track.eta(), collision.centFT0C());
      } else if (particle.pdgCode() == kKMinus) { // Neg kaons
        rMCRegistry.fill(HIST("hRecKaonN"), track.pt(), track.eta(), collision.centFT0C());
      } else if (particle.pdgCode() == kProton) { // Pos protons
        rMCRegistry.fill(HIST("hRecProtonP"), track.pt(), track.eta(), collision.centFT0C());
      } else if (particle.pdgCode() == kProtonBar) { // Neg protons
        rMCRegistry.fill(HIST("hRecProtonN"), track.pt(), track.eta(), collision.centFT0C());
      }

      if (trackFilters(track)) {

        // Track efficiency - Reconstructed & PID filters applied
        assocPID = trackPID(track);
        rMCRegistry.fill(HIST("hIdentified"), track.pt(), track.eta(), collision.centFT0C());
        if (assocPID[0] == pionID && track.sign() > 0) { // Pos pions
          rMCRegistry.fill(HIST("hPIDPionP"), track.pt(), track.eta(), collision.centFT0C());
        } else if (assocPID[0] == pionID && track.sign() < 0) { // Neg pions
          rMCRegistry.fill(HIST("hPIDPionN"), track.pt(), track.eta(), collision.centFT0C());
        } else if (assocPID[0] == kaonID && track.sign() > 0) { // Pos kaons
          rMCRegistry.fill(HIST("hPIDKaonP"), track.pt(), track.eta(), collision.centFT0C());
        } else if (assocPID[0] == kaonID && track.sign() < 0) { // Neg kaons
          rMCRegistry.fill(HIST("hPIDKaonN"), track.pt(), track.eta(), collision.centFT0C());
        } else if (assocPID[0] == protonID && track.sign() > 0) { // Pos protons
          rMCRegistry.fill(HIST("hPIDProtonP"), track.pt(), track.eta(), collision.centFT0C());
        } else if (assocPID[0] == protonID && track.sign() < 0) { // Neg protons
          rMCRegistry.fill(HIST("hPIDProtonN"), track.pt(), track.eta(), collision.centFT0C());
        }

        // Purity (PID)
        if (track.sign() > 0) {        // Positive tracks
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
          if (assocPID[0] == pionID) { // Pions
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

    for (const auto& v0 : v0s) {

      if (!v0.has_mcParticle() || v0.pt() <= v0PtMin || v0.pt() >= v0PtMax || std::abs(v0.eta()) >= v0EtaMax) {
        continue;
      }
      auto particle = v0.mcParticle();

      if (particle.isPhysicalPrimary()) {

        // V0 efficiency - Reconstructed
        if (particle.pdgCode() == kLambda0) { // Lambdas
          rMCRegistry.fill(HIST("hRecLambdaP"), v0.pt(), v0.eta(), collision.centFT0C());
        } else if (particle.pdgCode() == kLambda0Bar) { // AntiLambdas
          rMCRegistry.fill(HIST("hRecLambdaN"), v0.pt(), v0.eta(), collision.centFT0C());
        }

        if (v0Filters(collision, v0, tracks)) {

          // V0 efficiency - Reconstructed
          if (v0Sign(v0) == 1) { // Lambdas
            candMass = v0.mLambda();
            rQARegistry.fill(HIST("hInvMassLambda_MC"), v0.mLambda(), v0.pt(), collision.centFT0C());
            if (candMass >= MassLambda0 - 4 * dGaussSigma && candMass <= MassLambda0 + 4 * dGaussSigma) {
              rMCRegistry.fill(HIST("hPIDLambdaP_SGNL"), v0.pt(), v0.eta(), collision.centFT0C());
            } else if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
              rMCRegistry.fill(HIST("hPIDLambdaP_SB"), v0.pt(), v0.eta(), collision.centFT0C());
            }

          } else if (v0Sign(v0) == -1) { // AntiLambdas
            candMass = v0.mAntiLambda();
            rQARegistry.fill(HIST("hInvMassAntiLambda_MC"), v0.mAntiLambda(), v0.pt(), collision.centFT0C());
            if (candMass >= MassLambda0 - 4 * dGaussSigma && candMass <= MassLambda0 + 4 * dGaussSigma) {
              rMCRegistry.fill(HIST("hPIDLambdaN_SGNL"), v0.pt(), v0.eta(), collision.centFT0C());
            } else if (candMass >= MassLambda0 - 8 * dGaussSigma && candMass <= MassLambda0 + 8 * dGaussSigma) {
              rMCRegistry.fill(HIST("hPIDLambdaN_SB"), v0.pt(), v0.eta(), collision.centFT0C());
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

  template <class TrackCand>
  double trackEff(TH3D** efficiencies, const TrackCand& track, double centrality)
  {

    int index = -999;
    if (track.sign() > 0) {
      index = 0;
    } else if (track.sign() < 0) {
      index = 1;
    }

    double efficiency = efficiencies[index]->GetBinContent(efficiencies[index]->FindBin(track.pt(), track.eta(), centrality));
    if (efficiency > 0) {
      return efficiency;
    } else {
      return 1.0;
    }
  }

  template <class V0Cand>
  int v0Sign(const V0Cand& v0)
  {

    if (std::abs(v0.mLambda() - MassLambda0) <= std::abs(v0.mAntiLambda() - MassLambda0)) {
      return 1;
    } else if (std::abs(v0.mLambda() - MassLambda0) > std::abs(v0.mAntiLambda() - MassLambda0)) {
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

  template <class Col>
  bool acceptEvent(const Col& col, bool FillHist) // Event filter
  {

    if (FillHist) {
      rQARegistry.fill(HIST("hNEvents"), 0.5);
    }

    if (!col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) { // kIsGoodZvtxFT0vsPV
      return false;
    }
    if (FillHist) {
      rQARegistry.fill(HIST("hNEvents"), 1.5);
    }

    if (!col.selection_bit(aod::evsel::kNoSameBunchPileup)) { // kNoSameBunchPileup
      return false;
    }
    if (FillHist) {
      rQARegistry.fill(HIST("hNEvents"), 2.5);
    }
    
    return true;
  }

  template <class Col, class V0Cand, typename T>
  bool v0Filters(const Col& col, const V0Cand& v0, T const&) // V0 filter
  {
    
    // Kinematic cuts
    if (v0.pt() <= v0PtMin || v0.pt() >= v0PtMax || std::abs(v0.eta()) >= v0EtaMax) {
      return false;
    }

    // Daughter cuts
    auto posDaughter = v0.template posTrack_as<T>();
    auto negDaughter = v0.template negTrack_as<T>();
    if (std::abs(posDaughter.eta()) >= trackEtaMax || std::abs(negDaughter.eta()) >= trackEtaMax) {
      return false;
    }
    if (posDaughter.tpcNClsCrossedRows() <= tpcNCrossedRowsMin || negDaughter.tpcNClsCrossedRows() <= tpcNCrossedRowsMin) {
      return false;
    }
    if (v0Sign(v0) == 1) {
      if (std::abs(posDaughter.tpcNSigmaPr()) >= nSigma5 || std::abs(negDaughter.tpcNSigmaPi()) >= nSigma5) {
	return false;
      }
      if (std::abs(v0.dcapostopv()) <= dcaProtonMin || std::abs(v0.dcanegtopv()) <= dcaPionMin) {
	return false;
      }
    } else if (v0Sign(v0) == -1) {
      if (std::abs(posDaughter.tpcNSigmaPi()) >= nSigma5 || std::abs(negDaughter.tpcNSigmaPr()) >= nSigma5) {
	return false;
      }
      if (std::abs(v0.dcapostopv()) <= dcaPionMin || std::abs(v0.dcanegtopv()) <= dcaProtonMin) {
	return false;
      }
    }

    // Topological cuts
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (v0.v0radius() <= decayRMin) {
      return false;
    }
    if (ctau >= ctauMax) {
      return false;
    }
    if (v0.v0cosPA() <= cosPAMin) {
      return false;
    }
    if (v0.dcaV0daughters() >= dcaV0DauMax) {
      return false;
    }

    return true;
  }

  template <class TrackCand>
  bool trackFilters(const TrackCand& track) // Track filter
  {

    if (!track.hasTOF()) {
      return false;
    }

    if (trackPID(track)[0] == pionID) { // Pions
      if (std::abs(track.tpcNSigmaPi()) >= nSigma4) {
        return false;
      }
      if (track.pt() < pionPtMin) {
        return false;
      } else if (track.pt() > pionPtMin && track.pt() < pionPtMid) {
        if (std::abs(track.tofNSigmaPi()) >= nSigma4) {
          return false;
        }
      } else if (track.pt() > pionPtMid && track.pt() < pionPtMax) {
        if (track.tofNSigmaPi() <= -nSigma4 || track.tofNSigmaPi() >= nSigma0) {
          return false;
        }
      } else if (track.pt() > pionPtMax) {
        return false;
      }

    } else if (trackPID(track)[0] == kaonID) { // Kaons
      if (std::abs(track.tpcNSigmaKa()) >= nSigma4) {
        return false;
      }
      if (track.pt() < kaonPtMin) {
        return false;
      } else if (track.pt() > kaonPtMin && track.pt() < kaonPtMid1) {
        if (std::abs(track.tofNSigmaKa()) >= nSigma4) {
          return false;
        }
      } else if (track.pt() > kaonPtMid1 && track.pt() < kaonPtMid2) {
        if (track.tofNSigmaKa() <= -nSigma2 || track.tofNSigmaKa() >= nSigma4) {
          return false;
        }
      } else if (track.pt() > kaonPtMid2 && track.pt() < kaonPtMax) {
        if (track.tofNSigmaKa() <= nSigma0 || track.tofNSigmaKa() >= nSigma4) {
          return false;
        }
      } else if (track.pt() > kaonPtMax) {
        return false;
      }

    } else if (trackPID(track)[0] == protonID) { // Protons
      if (std::abs(track.tpcNSigmaPr()) >= nSigma4) {
        return false;
      }
      if (track.pt() < protonPtMin) {
        return false;
      } else if (track.pt() > protonPtMin && track.pt() < protonPtMid) {
        if (track.tofNSigmaPr() <= -nSigma2 || track.tofNSigmaPr() >= nSigma4) {
          return false;
        }
      } else if (track.pt() > protonPtMid && track.pt() < protonPtMax) {
        if (std::abs(track.tofNSigmaPr()) >= nSigma4) {
          return false;
        }
      } else if (track.pt() > protonPtMax) {
        if (track.tofNSigmaPr() <= -nSigma2 || track.tofNSigmaPr() >= nSigma4) {
          return false;
        }
      }
    }

    return true;
  }

  template <class V0Cand, class TrackCand>
  bool correlationFilters(const V0Cand& v0, const TrackCand& track) // Correlation filter
  {

    if (track.globalIndex() == v0.posTrackId() || track.globalIndex() == v0.negTrackId()) {
      return false;
    }

    return true;
  }

  template <class V0Cand, class TrackCand>
  bool fakeV0Filter(const V0Cand& v0, const TrackCand& track)
  {

    if (confFakeV0Switch) {

      if (trackPID(track)[0] == kaonID) { // Kaons
        return true;
      }

      std::array<float, 2> massArray;
      std::array<float, 3> dMomArray;
      std::array<float, 3> aMomArray = track.pVector();
      if (trackPID(track)[0] == pionID) {
        massArray = {MassProton, MassPionCharged};

        if (v0Sign(v0) == 1 && track.sign() == -1) { // Lambda - Pi_min
          const auto& dTrack = v0.template posTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        } else if (v0Sign(v0) == -1 && track.sign() == 1) { // Antilambda - Pi_plus
          const auto& dTrack = v0.template negTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        }
      } else if (trackPID(track)[0] == protonID) {
        massArray = {MassPionCharged, MassProton};

        if (v0Sign(v0) == 1 && track.sign() == 1) { // Lambda - Proton
          const auto& dTrack = v0.template negTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        } else if (v0Sign(v0) == -1 && track.sign() == -1) { // Antilambda - Antiproton
          const auto& dTrack = v0.template posTrack_as<MyFilteredTracks>();
          dMomArray = dTrack.pVector();
        }
      }

      double invMass = RecoDecay::m(std::array{dMomArray, aMomArray}, massArray);
      if (invMass >= MassLambda0 - 4 * dGaussSigma && invMass <= MassLambda0 + 4 * dGaussSigma) {
        return false;
      }
    }

    return true;
  }

  template <class V0Cand, class TrackCand>
  bool radialDistanceFilter(const V0Cand& v0, const TrackCand& track, double B, bool Mix)
  {

    bool pass = true;
    if (confRDSwitch) {

      auto proton = v0.template posTrack_as<MyFilteredTracks>();
      if (v0Sign(v0) == -1) {
        proton = v0.template negTrack_as<MyFilteredTracks>();
      }

      double dEta = proton.eta() - track.eta();
      if (std::abs(dEta) > dEtaMax) {
        return pass;
      }

      double dPhiStar;
      double dPhi = proton.phi() - track.phi();
      double phaseProton = (-0.3 * B * proton.sign()) / (2 * proton.pt());
      double phaseTrack = (-0.3 * B * track.sign()) / (2 * track.pt());

      double dPhiStarMean = 0;

      // Start of the TPC radius loop
      for (double r = rMin; r <= rMax; r += 0.01) {
        dPhiStar = RecoDecay::constrainAngle(dPhi + std::asin(phaseProton * r) - std::asin(phaseTrack * r), -constants::math::PIHalf);

        if (r == rMin) {
          if (!Mix) {                                 // Same-event
            if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
              rPhiStarRegistry.fill(HIST("hSEProtonPreCut_OS"), dPhiStar, dEta);
            } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
              rPhiStarRegistry.fill(HIST("hSEProtonPreCut_SS"), dPhiStar, dEta);
              if (proton.sign() == 1) { // Positive
                rPhiStarRegistry.fill(HIST("hSEProtonPreCut_SSP"), dPhiStar, dEta);
              } else if (proton.sign() == -1) { // Negative
                rPhiStarRegistry.fill(HIST("hSEProtonPreCut_SSN"), dPhiStar, dEta);
              }
            }

          } else {                                    // Mixed-event
            if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
              rPhiStarRegistry.fill(HIST("hMEProtonPreCut_OS"), dPhiStar, dEta);
            } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
              rPhiStarRegistry.fill(HIST("hMEProtonPreCut_SS"), dPhiStar, dEta);
              if (proton.sign() == 1) { // Positive
                rPhiStarRegistry.fill(HIST("hMEProtonPreCut_SSP"), dPhiStar, dEta);
              } else if (proton.sign() == -1) { // Negative
                rPhiStarRegistry.fill(HIST("hMEProtonPreCut_SSN"), dPhiStar, dEta);
              }
            }
          }
        }

        if (std::abs(dEta) <= dEtaMin) {
          if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
            if (std::abs(dPhiStar) <= dPhiStarMinOS) {
              pass = false;
            }
          } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
            if (std::abs(dPhiStar) <= dPhiStarMinSS) {
              pass = false;
            }
          }
        }

        if (r == rMin && pass) {
          if (!Mix) {                                 // Same-event
            if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
              rPhiStarRegistry.fill(HIST("hSEProtonPostCut_OS"), dPhiStar, dEta);
            } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
              rPhiStarRegistry.fill(HIST("hSEProtonPostCut_SS"), dPhiStar, dEta);
              if (proton.sign() == 1) { // Positive
                rPhiStarRegistry.fill(HIST("hSEProtonPostCut_SSP"), dPhiStar, dEta);
              } else if (proton.sign() == -1) { // Negative
                rPhiStarRegistry.fill(HIST("hSEProtonPostCut_SSN"), dPhiStar, dEta);
              }
            }

          } else {                                    // Mixed-event
            if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
              rPhiStarRegistry.fill(HIST("hMEProtonPostCut_OS"), dPhiStar, dEta);
            } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
              rPhiStarRegistry.fill(HIST("hMEProtonPostCut_SS"), dPhiStar, dEta);
              if (proton.sign() == 1) { // Positive
                rPhiStarRegistry.fill(HIST("hMEProtonPostCut_SSP"), dPhiStar, dEta);
              } else if (proton.sign() == -1) { // Negative
                rPhiStarRegistry.fill(HIST("hMEProtonPostCut_SSN"), dPhiStar, dEta);
              }
            }
          }
        }

        dPhiStarMean += (dPhiStar / 170);
      }
      // End of the TPC radius loop

      if (!Mix) {                                 // Same-event
        if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
          rPhiStarRegistry.fill(HIST("hSEPhiStarMean_OS"), dPhiStarMean, dEta);
        } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
          rPhiStarRegistry.fill(HIST("hSEPhiStarMean_SS"), dPhiStarMean, dEta);
          if (proton.sign() == 1) { // Positive
            rPhiStarRegistry.fill(HIST("hSEPhiStarMean_SSP"), dPhiStarMean, dEta);
          } else if (proton.sign() == -1) { // Negative
            rPhiStarRegistry.fill(HIST("hSEPhiStarMean_SSN"), dPhiStarMean, dEta);
          }
        }

      } else {                                    // Mixed-event
        if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
          rPhiStarRegistry.fill(HIST("hMEPhiStarMean_OS"), dPhiStarMean, dEta);
        } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
          rPhiStarRegistry.fill(HIST("hMEPhiStarMean_SS"), dPhiStarMean, dEta);
          if (proton.sign() == 1) { // Positive
            rPhiStarRegistry.fill(HIST("hMEPhiStarMean_SSP"), dPhiStarMean, dEta);
          } else if (proton.sign() == -1) { // Negative
            rPhiStarRegistry.fill(HIST("hMEPhiStarMean_SSN"), dPhiStarMean, dEta);
          }
        }
      }
    }

    return pass;
  }
};

//============================================================================================================================================================================================================================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<ThreeParticleCorrelations>(cfgc)};
  return workflow;
}

//============================================================================================================================================================================================================================================================================
