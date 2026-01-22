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

#include "RecoDecay.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TPDGCode.h"

#include <algorithm>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::physics;

struct ThreeParticleCorrelations {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Analysis parameters
  float centMin = 0.0, centMax = 90.0;
  float v0PtMin = 0.6, v0PtMax = 12.0;
  float v0EtaMax = 0.72;
  float trackPtMin = 0.2, trackPtMax = 3.0;
  float trackEtaMax = 0.8;

  // Track PID parameters
  double pionID = 0.0, kaonID = 1.0, protonID = 2.0;
  float nSigma0 = 0.0, nSigma1 = 1.0, nSigma2 = 2.0, nSigma4 = 4.0, nSigma5 = 5.0;

  // Event selection parameters
  struct : ConfigurableGroup {
    std::string prefix = "EventSelection";
    Configurable<float> zvtxMax{"zvtxMax", 10.0, "Maximum collision Z-vertex position (cm)"};
    Configurable<int> occupMin{"occupMin", 0, "Minimum collision occupancy"};
    Configurable<int> occupMax{"occupMax", 15000, "Maximum collision occupancy"};
    Configurable<bool> useOccupCut{"useOccupCut", true, "Use the kNoCollInTimeRangeStandard cut"};
  } evSelGroup;

  // V0 filter parameters
  struct : ConfigurableGroup {
    std::string prefix = "V0Selection";
    Configurable<float> tpcNCrossedRows{"tpcNCrossedRows", 70.0, "Minimum number of TPC crossed rows"};
    Configurable<float> decayR{"decayR", 1.2, "Minimum V0 decay radius (cm)"};
    Configurable<float> ctau{"ctau", 30.0, "Maximum V0 proper lifetime (cm)"};
    Configurable<float> cosPA{"cosPA", 0.995, "Minimum V0 cosine of pointing angle"};
    Configurable<float> dcaProton{"dcaProton", 0.05, "Minimum DCA of proton daughter (cm)"};
    Configurable<float> dcaPion{"dcaPion", 0.2, "Minimum DCA of pion daughter (cm)"};
    Configurable<float> dcaV0Dau{"dcaV0Dau", 1.0, "Maximum DCA between V0 daughters"};
  } v0SelGroup;

  // Track filter parameters
  float pionPtMin = 0.3, pionPtMax = 2.3, kaonPtMin = 0.5, kaonPtMax = 2.3, protonPtMin = 0.6;
  float pionPtMid1 = 1.6, pionPtMid2 = 2.0, kaonPtMid1 = 1.5, kaonPtMid2 = 2.0, protonPtMid = 2.3;

  // RD filter parameters
  float dEtaMax = 0.05, dEtaMin = 0.023;
  float dPhiStarMinOS = 0.09, dPhiStarMinSS = 0.095;
  float rMin = 0.8, rMax = 2.5;

  // Lambda invariant mass fit
  Configurable<float> invMassNSigma{"invMassNSigma", 4.0, "Number of standard deviations from the mean of the Lambda invariant mass peak"};
  double dGaussSigma = 0.002;

  // Histogram registry
  HistogramRegistry rMECorrRegistry{"MECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rSECorrRegistry{"SECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rMCRegistry{"MCRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rPhiStarRegistry{"PhiStarRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry rQARegistry{"QARegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Collision & Event filters
  Filter collCent = aod::cent::centFT0C > centMin&& aod::cent::centFT0C < centMax;
  Filter collZvtx = nabs(aod::collision::posZ) < evSelGroup.zvtxMax;
  Filter mcCollZvtx = nabs(aod::mccollision::posZ) < evSelGroup.zvtxMax;
  Filter evSel8 = aod::evsel::sel8 == true;
  Filter evSelOccup = o2::aod::evsel::trackOccupancyInTimeRange >= evSelGroup.occupMin && o2::aod::evsel::trackOccupancyInTimeRange < evSelGroup.occupMax;

  // Track filters
  Filter trackPt = aod::track::pt > trackPtMin&& aod::track::pt < trackPtMax;
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
  Partition<MyFilteredMCParticles> mcTracks = aod::mcparticle::pt > trackPtMin&& aod::mcparticle::pt < trackPtMax;
  Partition<MyFilteredMCParticles> mcV0s = aod::mcparticle::pt > v0PtMin&& aod::mcparticle::pt < v0PtMax&& nabs(aod::mcparticle::eta) < v0EtaMax;
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
  ConfigurableAxis confZvtxBins{"confZvtxBins", {VARIABLE_WIDTH, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f}, "ME Zvtx binning"};
  using BinningType = ColumnBinningPolicy<aod::cent::CentFT0C, aod::collision::PosZ>;
  using BinningTypeMC = ColumnBinningPolicy<aod::mccollisionprop::BestCollisionCentFT0C, aod::mccollision::PosZ>;

  BinningType collBinning{{confCentBins, confZvtxBins}, true};
  BinningTypeMC collBinningMC{{confCentBins, confZvtxBins}, true};
  Pair<MyFilteredCollisions, aod::V0Datas, MyFilteredTracks, BinningType> pairData{collBinning, 5, -1, &cache};
  SameKindPair<MyFilteredMCGenCollisions, MyFilteredMCParticles, BinningTypeMC> pairMC{collBinningMC, 5, -1, &cache};

  // Process configurables
  struct : ConfigurableGroup {
    std::string prefix = "processSwitchBoard";
    Configurable<int> confBfieldSwitch{"confBfieldSwitch", 0, "Switch for the detector magnetic field (1 if Pos, -1 if Neg, 0 if both)"};
    Configurable<bool> confRatioCorrectionSwitch{"confRatioCorrectionSwitch", false, "Switch for correcting the negative spectra back to the positive spectra"};
    Configurable<bool> confFakeV0Switch{"confFakeV0Switch", false, "Switch for the fakeV0Filter function"};
    Configurable<bool> confRDSwitch{"confRDSwitch", true, "Switch for the radialDistanceFilter function"};
  } switchGroup;

  // Efficiency histograms
  TH3D** hEffPions = new TH3D*[2];
  TH3D** hEffKaons = new TH3D*[2];
  TH3D** hEffProtons = new TH3D*[2];
  TH3D** hEffLambdas = new TH3D*[2];

  // Spectra correction histograms
  TH2D** hCorrectionPions = new TH2D*[2];
  TH2D** hCorrectionKaons = new TH2D*[2];
  TH2D** hCorrectionProtons = new TH2D*[2];

  // Correlation variables
  int triggSign, assocSign;
  double v0Efficiency;
  double candMass;
  double* assocPID;

  double deltaPhi, deltaEta;

  //==========================================================================================================================================================================================================================================================================

  void init(InitContext const&)
  {

    TH1::SetDefaultSumw2(true);
    
    // Bins of variable width
    std::vector<double> fineCentBins = {0.0, 2.0, 4.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0};

    // Histograms axes
    const AxisSpec centralityAxis{confCentBins};
    const AxisSpec fineCentralityAxis{fineCentBins};
    const AxisSpec zvtxAxis{confZvtxBins};
    const AxisSpec occupancyAxis{200, 0, 20000};
    const AxisSpec dPhiAxis{36, (-1. / 2) * constants::math::PI, (3. / 2) * constants::math::PI};
    const AxisSpec dEtaAxis{32, -1.52, 1.52};
    const AxisSpec v0PtAxis{114, 0.6, 12};
    const AxisSpec v0EtaAxis{36, -0.72, 0.72};
    const AxisSpec trackPtAxis{28, 0.2, 3};
    const AxisSpec trackEtaAxis{32, -0.8, 0.8};
    const AxisSpec lambdaInvMassAxis{100, 1.08, 1.16};

    // QA & PID
    rQARegistry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{5, 0, 5}}});
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "All");
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "kIsGoodZvtxFT0vsPV");
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "kNoSameBunchPileup");
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, Form("[%i < Occupancy < %i)", static_cast<int>(evSelGroup.occupMin), static_cast<int>(evSelGroup.occupMax)));
    rQARegistry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "kNoCollInTimeRangeStandard");

    rQARegistry.add("hEventCentrality", "hEventCentrality", {HistType::kTH1D, {{fineCentralityAxis}}});
    rQARegistry.add("hEventCentrality_MC", "hEventCentrality_MC", {HistType::kTH1D, {{fineCentralityAxis}}});
    rQARegistry.add("hEventZvtx", "hEventZvtx", {HistType::kTH1D, {{zvtxAxis}}});
    rQARegistry.add("hEventOccupancy", "hEventOccupancy", {HistType::kTH1D, {{occupancyAxis}}});
    rQARegistry.add("hEventBfield", "hEventBfield", {HistType::kTH1D, {{2, -1, 1}}});
    rQARegistry.add("hTrackPt", "hTrackPt", {HistType::kTH1D, {{100, 0, 4}}});
    rQARegistry.add("hTrackEta", "hTrackEta", {HistType::kTH1D, {{100, -1, 1}}});
    rQARegistry.add("hTrackPhi", "hTrackPhi", {HistType::kTH1D, {{100, (-1. / 2) * constants::math::PI, (5. / 2) * constants::math::PI}}});
    rQARegistry.add("hTrackNSharedClusters", "hTrackNSharedClusters", {HistType::kTH1D, {{200, 0, 200}}});

    rQARegistry.add("hPtPion_Uncorrected", "hPtPion_Uncorrected", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtKaon_Uncorrected", "hPtKaon_Uncorrected", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtProton_Uncorrected", "hPtProton_Uncorrected", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtV0_Uncorrected", "hPtV0_Uncorrected", {HistType::kTH3D, {{v0PtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtPion_Corrected", "hPtPion_Corrected", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtKaon_Corrected", "hPtKaon_Corrected", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtProton_Corrected", "hPtProton_Corrected", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtV0_Corrected", "hPtV0_Corrected", {HistType::kTH3D, {{v0PtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtPion_Looped", "hPtPion_Looped", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtKaon_Looped", "hPtKaon_Looped", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtProton_Looped", "hPtProton_Looped", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtPion_MC", "hPtPion_MC", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtKaon_MC", "hPtKaon_MC", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtProton_MC", "hPtProton_MC", {HistType::kTH3D, {{trackPtAxis}, {fineCentralityAxis}, {2, -2, 2}}});
    rQARegistry.add("hPtV0_MC", "hPtV0_MC", {HistType::kTH3D, {{v0PtAxis}, {fineCentralityAxis}, {2, -2, 2}}});

    rQARegistry.add("hdEdx", "hdEdx", {HistType::kTH2D, {{120, -3.0, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxPion", "hdEdxPion", {HistType::kTH2D, {{120, -3.0, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxKaon", "hdEdxKaon", {HistType::kTH2D, {{120, -3.0, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hdEdxProton", "hdEdxProton", {HistType::kTH2D, {{120, -3.0, 3.0}, {180, 20, 200}}});
    rQARegistry.add("hBeta", "hBeta", {HistType::kTH2D, {{120, -3.0, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaPion", "hBetaPion", {HistType::kTH2D, {{120, -3.0, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaKaon", "hBetaKaon", {HistType::kTH2D, {{120, -3.0, 3.0}, {70, 0.4, 1.1}}});
    rQARegistry.add("hBetaProton", "hBetaProton", {HistType::kTH2D, {{120, -3.0, 3.0}, {70, 0.4, 1.1}}});

    rQARegistry.add("hTPCPion", "hTPCPion", {HistType::kTH3D, {{trackPtAxis}, {1001, -50.05, 50.05}, {2, -2, 2}}});
    rQARegistry.add("hTPCKaon", "hTPCKaon", {HistType::kTH3D, {{trackPtAxis}, {1001, -50.05, 50.05}, {2, -2, 2}}});
    rQARegistry.add("hTPCProton", "hTPCProton", {HistType::kTH3D, {{trackPtAxis}, {1001, -50.05, 50.05}, {2, -2, 2}}});
    rQARegistry.add("hTOFPion", "hTOFPion", {HistType::kTH3D, {{trackPtAxis}, {1001, -50.05, 50.05}, {2, -2, 2}}});
    rQARegistry.add("hTOFKaon", "hTOFKaon", {HistType::kTH3D, {{trackPtAxis}, {1001, -50.05, 50.05}, {2, -2, 2}}});
    rQARegistry.add("hTOFProton", "hTOFProton", {HistType::kTH3D, {{trackPtAxis}, {1001, -50.05, 50.05}, {2, -2, 2}}});

    rQARegistry.add("hInvMassLambda", "hInvMassLambda", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hInvMassAntiLambda", "hInvMassAntiLambda", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hNLambdas", "hNLambdas", {HistType::kTH3D, {{2, -2, 2}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hInvMassLambda_MC", "hInvMassLambda_MC", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});
    rQARegistry.add("hInvMassAntiLambda_MC", "hInvMassAntiLambda_MC", {HistType::kTH3D, {{lambdaInvMassAxis}, {v0PtAxis}, {centralityAxis}}});

    // PhiStar
    rPhiStarRegistry.add("hSEPhiStarIR_OS", "hSEPhiStarIR_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarIR_SS", "hSEPhiStarIR_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarIR_SSP", "hSEPhiStarIR_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarIR_SSN", "hSEPhiStarIR_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_OS", "hSEPhiStarMean_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_SS", "hSEPhiStarMean_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_SSP", "hSEPhiStarMean_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarMean_SSN", "hSEPhiStarMean_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hSEPhiStarRadial_OS", "hSEPhiStarRadial_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {170, 0.8, 2.5}}});
    rPhiStarRegistry.add("hSEPhiStarRadial_SS", "hSEPhiStarRadial_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {170, 0.8, 2.5}}});

    rPhiStarRegistry.add("hMEPhiStarIR_OS", "hMEPhiStarIR_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarIR_SS", "hMEPhiStarIR_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarIR_SSP", "hMEPhiStarIR_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarIR_SSN", "hMEPhiStarIR_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_OS", "hMEPhiStarMean_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_SS", "hMEPhiStarMean_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_SSP", "hMEPhiStarMean_SSP", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarMean_SSN", "hMEPhiStarMean_SSN", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {101, -0.0505, 0.0505}}});
    rPhiStarRegistry.add("hMEPhiStarRadial_OS", "hMEPhiStarRadial_OS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {170, 0.8, 2.5}}});
    rPhiStarRegistry.add("hMEPhiStarRadial_SS", "hMEPhiStarRadial_SS", {HistType::kTH2D, {{121, -0.3025, 0.3025}, {170, 0.8, 2.5}}});

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

    // Purity
    rMCRegistry.add("hSelectPionP", "hSelectPionP", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hSelectPionN", "hSelectPionN", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hSelectKaonP", "hSelectKaonP", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hSelectKaonN", "hSelectKaonN", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hSelectProtonP", "hSelectProtonP", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hSelectProtonN", "hSelectProtonN", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hTrueSelectPionP", "hTrueSelectPionP", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hTrueSelectPionN", "hTrueSelectPionN", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hTrueSelectKaonP", "hTrueSelectKaonP", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hTrueSelectKaonN", "hTrueSelectKaonN", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hTrueSelectProtonP", "hTrueSelectProtonP", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});
    rMCRegistry.add("hTrueSelectProtonN", "hTrueSelectProtonN", {HistType::kTH2D, {{trackPtAxis}, {centralityAxis}}});

    // Correlations
    rSECorrRegistry.add("hSameLambdaPion_SGNL", "Same-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaPion_SB", "Same-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaKaon_SGNL", "Same-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaKaon_SB", "Same-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaProton_SGNL", "Same-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaProton_SB", "Same-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaPion_MC", "Same-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaKaon_MC", "Same-event #Lambda - K correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaProton_MC", "Same-event #Lambda - p correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);

    rSECorrRegistry.add("hSameLambdaPion_leftSB", "Same-event #Lambda - #pi correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaPion_rightSB", "Same-event #Lambda - #pi correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaKaon_leftSB", "Same-event #Lambda - K correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaKaon_rightSB", "Same-event #Lambda - K correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaProton_leftSB", "Same-event #Lambda - p correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rSECorrRegistry.add("hSameLambdaProton_rightSB", "Same-event #Lambda - p correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);

    rMECorrRegistry.add("hMixLambdaPion_SGNL", "Mixed-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaPion_SB", "Mixed-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaKaon_SGNL", "Mixed-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaKaon_SB", "Mixed-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaProton_SGNL", "Mixed-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaProton_SB", "Mixed-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaPion_MC", "Mixed-event #Lambda - #pi correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaKaon_MC", "Mixed-event #Lambda - K correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaProton_MC", "Mixed-event #Lambda - p correlator (MC)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);

    rMECorrRegistry.add("hMixLambdaPion_leftSB", "Mixed-event #Lambda - #pi correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaPion_rightSB", "Mixed-event #Lambda - #pi correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaKaon_leftSB", "Mixed-event #Lambda - K correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaKaon_rightSB", "Mixed-event #Lambda - K correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaProton_leftSB", "Mixed-event #Lambda - p correlator (Left SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);
    rMECorrRegistry.add("hMixLambdaProton_rightSB", "Mixed-event #Lambda - p correlator (Right SB region)", {HistType::kTHnSparseF, {{dPhiAxis}, {dEtaAxis}, {centralityAxis}, {zvtxAxis}, {2, -2, 2}, {2, -2, 2}}}, true);

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    TList* effListChargedParticles = ccdb->getForTimeStamp<TList>("Users/j/jstaa/Efficiency/ChargedParticles", 1);
    TList* effListLambdas = ccdb->getForTimeStamp<TList>("Users/j/jstaa/Efficiency/Lambdas", 1);
    hEffPions[0] = static_cast<TH3D*>(effListChargedParticles->FindObject("hEfficiencyPionP"));
    hEffPions[1] = static_cast<TH3D*>(effListChargedParticles->FindObject("hEfficiencyPionN"));
    hEffKaons[0] = static_cast<TH3D*>(effListChargedParticles->FindObject("hEfficiencyKaonP"));
    hEffKaons[1] = static_cast<TH3D*>(effListChargedParticles->FindObject("hEfficiencyKaonN"));
    hEffProtons[0] = static_cast<TH3D*>(effListChargedParticles->FindObject("hEfficiencyProtonP"));
    hEffProtons[1] = static_cast<TH3D*>(effListChargedParticles->FindObject("hEfficiencyProtonN"));
    hEffLambdas[0] = static_cast<TH3D*>(effListLambdas->FindObject("hEfficiencyLambdaP"));
    hEffLambdas[1] = static_cast<TH3D*>(effListLambdas->FindObject("hEfficiencyLambdaN"));

    TList* correctionListChargedParticles = ccdb->getForTimeStamp<TList>("Users/j/jstaa/SpectraRatios/ChargedParticles", 1);
    hCorrectionPions[0] = static_cast<TH2D*>(correctionListChargedParticles->FindObject("h2DRatioPionP"));
    hCorrectionPions[1] = static_cast<TH2D*>(correctionListChargedParticles->FindObject("h2DRatioPionN"));
    hCorrectionKaons[0] = static_cast<TH2D*>(correctionListChargedParticles->FindObject("h2DRatioKaonP"));
    hCorrectionKaons[1] = static_cast<TH2D*>(correctionListChargedParticles->FindObject("h2DRatioKaonN"));
    hCorrectionProtons[0] = static_cast<TH2D*>(correctionListChargedParticles->FindObject("h2DRatioProtonP"));
    hCorrectionProtons[1] = static_cast<TH2D*>(correctionListChargedParticles->FindObject("h2DRatioProtonN"));
  }

  //==========================================================================================================================================================================================================================================================================

  void processSame(MyFilteredCollision const& collision, aod::V0Datas const& v0s, MyFilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    if (!acceptEvent(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());
    if (switchGroup.confBfieldSwitch != 0) {
      if (std::signbit(static_cast<double>(switchGroup.confBfieldSwitch)) != std::signbit(bField)) {
        return;
      }
    }

    rQARegistry.fill(HIST("hEventCentrality"), collision.centFT0C());
    rQARegistry.fill(HIST("hEventZvtx"), collision.posZ());
    rQARegistry.fill(HIST("hEventOccupancy"), collision.trackOccupancyInTimeRange());
    rQARegistry.fill(HIST("hEventBfield"), bField);

    // Start of the Track QA
    for (const auto& track : tracks) {
      rQARegistry.fill(HIST("hTPCPion"), track.pt(), track.tpcNSigmaPi(), track.sign());
      rQARegistry.fill(HIST("hTPCKaon"), track.pt(), track.tpcNSigmaKa(), track.sign());
      rQARegistry.fill(HIST("hTPCProton"), track.pt(), track.tpcNSigmaPr(), track.sign());
      if (track.hasTOF()) {
        rQARegistry.fill(HIST("hTOFPion"), track.pt(), track.tofNSigmaPi(), track.sign());
        rQARegistry.fill(HIST("hTOFKaon"), track.pt(), track.tofNSigmaKa(), track.sign());
        rQARegistry.fill(HIST("hTOFProton"), track.pt(), track.tofNSigmaPr(), track.sign());
      }

      if (trackFilters(track)) {
        assocPID = trackPID(track);
        rQARegistry.fill(HIST("hTrackPt"), track.pt());
        rQARegistry.fill(HIST("hTrackEta"), track.eta());
        rQARegistry.fill(HIST("hTrackPhi"), track.phi());
        rQARegistry.fill(HIST("hTrackNSharedClusters"), track.tpcNClsShared());
        rQARegistry.fill(HIST("hdEdx"), track.sign() * track.pt(), track.tpcSignal());
        rQARegistry.fill(HIST("hBeta"), track.sign() * track.pt(), track.beta());
        if (assocPID[0] == pionID) { // Pions
          rQARegistry.fill(HIST("hPtPion_Uncorrected"), track.pt(), collision.centFT0C(), track.sign());
          rQARegistry.fill(HIST("hPtPion_Corrected"), track.pt(), collision.centFT0C(), track.sign(), ratioCorrection(hCorrectionPions, track, collision.centFT0C()) / trackEff(hEffPions, track, collision.centFT0C()));
          rQARegistry.fill(HIST("hdEdxPion"), track.sign() * track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaPion"), track.sign() * track.pt(), track.beta());
        } else if (assocPID[0] == kaonID) { // Kaons
          rQARegistry.fill(HIST("hPtKaon_Uncorrected"), track.pt(), collision.centFT0C(), track.sign());
          rQARegistry.fill(HIST("hPtKaon_Corrected"), track.pt(), collision.centFT0C(), track.sign(), ratioCorrection(hCorrectionKaons, track, collision.centFT0C()) / trackEff(hEffKaons, track, collision.centFT0C()));
          rQARegistry.fill(HIST("hdEdxKaon"), track.sign() * track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaKaon"), track.sign() * track.pt(), track.beta());
        } else if (assocPID[0] == protonID) { // Protons
          rQARegistry.fill(HIST("hPtProton_Uncorrected"), track.pt(), collision.centFT0C(), track.sign());
          rQARegistry.fill(HIST("hPtProton_Corrected"), track.pt(), collision.centFT0C(), track.sign(), ratioCorrection(hCorrectionProtons, track, collision.centFT0C()) / trackEff(hEffProtons, track, collision.centFT0C()));
          rQARegistry.fill(HIST("hdEdxProton"), track.sign() * track.pt(), track.tpcSignal());
          rQARegistry.fill(HIST("hBetaProton"), track.sign() * track.pt(), track.beta());
        }
      }
    }
    // End of the Track QA

    // Start of the Same-Event correlations
    for (const auto& trigger : v0s) {
      if (v0Filters(collision, trigger, tracks)) {

        triggSign = v0Sign(trigger);
        v0Efficiency = v0Eff(hEffLambdas, trigger, collision.centFT0C());

        rQARegistry.fill(HIST("hPtV0_Uncorrected"), trigger.pt(), collision.centFT0C(), triggSign);
        rQARegistry.fill(HIST("hPtV0_Corrected"), trigger.pt(), collision.centFT0C(), triggSign, 1. / v0Efficiency);
        if (triggSign == 1) {
          candMass = trigger.mLambda();
          rQARegistry.fill(HIST("hInvMassLambda"), trigger.mLambda(), trigger.pt(), collision.centFT0C(), 1. / v0Efficiency);
        } else if (triggSign == -1) {
          candMass = trigger.mAntiLambda();
          rQARegistry.fill(HIST("hInvMassAntiLambda"), trigger.mAntiLambda(), trigger.pt(), collision.centFT0C(), 1. / v0Efficiency);
        }

        for (const auto& associate : tracks) {
          if (trackFilters(associate)) {
            if (correlationFilters(trigger, associate) && radialDistanceFilter(trigger, associate, bField, false) && fakeV0Filter(trigger, associate)) {

              assocPID = trackPID(associate);
              deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
              deltaEta = trigger.eta() - associate.eta();

              if (assocPID[0] == pionID) { // Pions
                rQARegistry.fill(HIST("hPtPion_Looped"), associate.pt(), collision.centFT0C(), associate.sign(), ratioCorrection(hCorrectionPions, associate, collision.centFT0C()) / trackEff(hEffPions, associate, collision.centFT0C()));
              } else if (assocPID[0] == kaonID) { // Kaons
                rQARegistry.fill(HIST("hPtKaon_Looped"), associate.pt(), collision.centFT0C(), associate.sign(), ratioCorrection(hCorrectionKaons, associate, collision.centFT0C()) / trackEff(hEffKaons, associate, collision.centFT0C()));
              } else if (assocPID[0] == protonID) { // Protons
                rQARegistry.fill(HIST("hPtProton_Looped"), associate.pt(), collision.centFT0C(), associate.sign(), ratioCorrection(hCorrectionProtons, associate, collision.centFT0C()) / trackEff(hEffProtons, associate, collision.centFT0C()));
              }

              if (candMass >= MassLambda0 - invMassNSigma * dGaussSigma && candMass <= MassLambda0 + invMassNSigma * dGaussSigma) {
                if (assocPID[0] == pionID) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, collision.centFT0C()) / (trackEff(hEffPions, associate, collision.centFT0C()) * v0Efficiency));
                } else if (assocPID[0] == kaonID) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, collision.centFT0C()) / (trackEff(hEffKaons, associate, collision.centFT0C()) * v0Efficiency));
                } else if (assocPID[0] == protonID) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SGNL"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, collision.centFT0C()) / (trackEff(hEffProtons, associate, collision.centFT0C()) * v0Efficiency));
                }

              } else if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                if (assocPID[0] == pionID) { // Pions
                  rSECorrRegistry.fill(HIST("hSameLambdaPion_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, collision.centFT0C()) / (trackEff(hEffPions, associate, collision.centFT0C()) * v0Efficiency));
                  if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass < MassLambda0 - invMassNSigma * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaPion_leftSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, collision.centFT0C()) / (trackEff(hEffPions, associate, collision.centFT0C()) * v0Efficiency));
                  } else if (candMass > MassLambda0 + invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaPion_rightSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, collision.centFT0C()) / (trackEff(hEffPions, associate, collision.centFT0C()) * v0Efficiency));
                  }

                } else if (assocPID[0] == kaonID) { // Kaons
                  rSECorrRegistry.fill(HIST("hSameLambdaKaon_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, collision.centFT0C()) / (trackEff(hEffKaons, associate, collision.centFT0C()) * v0Efficiency));
                  if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass < MassLambda0 - invMassNSigma * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaKaon_leftSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, collision.centFT0C()) / (trackEff(hEffKaons, associate, collision.centFT0C()) * v0Efficiency));
                  } else if (candMass > MassLambda0 + invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaKaon_rightSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, collision.centFT0C()) / (trackEff(hEffKaons, associate, collision.centFT0C()) * v0Efficiency));
                  }

                } else if (assocPID[0] == protonID) { // Protons
                  rSECorrRegistry.fill(HIST("hSameLambdaProton_SB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, collision.centFT0C()) / (trackEff(hEffProtons, associate, collision.centFT0C()) * v0Efficiency));
                  if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass < MassLambda0 - invMassNSigma * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaProton_leftSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, collision.centFT0C()) / (trackEff(hEffProtons, associate, collision.centFT0C()) * v0Efficiency));
                  } else if (candMass > MassLambda0 + invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                    rSECorrRegistry.fill(HIST("hSameLambdaProton_rightSB"), deltaPhi, deltaEta, collision.centFT0C(), collision.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, collision.centFT0C()) / (trackEff(hEffProtons, associate, collision.centFT0C()) * v0Efficiency));
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
      if (!acceptEvent(coll_1, false)) {
        continue;
      }
      if (!acceptEvent(coll_2, false)) {
        continue;
      }

      auto bc = coll_1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      if (switchGroup.confBfieldSwitch != 0) {
        if (std::signbit(static_cast<double>(switchGroup.confBfieldSwitch)) != std::signbit(bField)) {
          continue;
        }
      }

      for (const auto& [trigger, associate] : soa::combinations(soa::CombinationsFullIndexPolicy(v0_1, track_2))) {
        if (v0Filters(coll_1, trigger, tracks) && trackFilters(associate)) {
          if (radialDistanceFilter(trigger, associate, bField, true) && fakeV0Filter(trigger, associate)) {

            triggSign = v0Sign(trigger);
            v0Efficiency = v0Eff(hEffLambdas, trigger, coll_1.centFT0C());

            if (triggSign == 1) {
              candMass = trigger.mLambda();
            } else if (triggSign == -1) {
              candMass = trigger.mAntiLambda();
            }

            assocPID = trackPID(associate);
            deltaPhi = RecoDecay::constrainAngle(trigger.phi() - associate.phi(), -constants::math::PIHalf);
            deltaEta = trigger.eta() - associate.eta();

            if (candMass >= MassLambda0 - invMassNSigma * dGaussSigma && candMass <= MassLambda0 + invMassNSigma * dGaussSigma) {
              if (assocPID[0] == pionID) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, coll_1.centFT0C()) / (trackEff(hEffPions, associate, coll_1.centFT0C()) * v0Efficiency));
              } else if (assocPID[0] == kaonID) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, coll_1.centFT0C()) / (trackEff(hEffKaons, associate, coll_1.centFT0C()) * v0Efficiency));
              } else if (assocPID[0] == protonID) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SGNL"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, coll_1.centFT0C()) / (trackEff(hEffProtons, associate, coll_1.centFT0C()) * v0Efficiency));
              }

            } else if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
              if (assocPID[0] == pionID) { // Pions
                rMECorrRegistry.fill(HIST("hMixLambdaPion_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, coll_1.centFT0C()) / (trackEff(hEffPions, associate, coll_1.centFT0C()) * v0Efficiency));
                if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass < MassLambda0 - invMassNSigma * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaPion_leftSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, coll_1.centFT0C()) / (trackEff(hEffPions, associate, coll_1.centFT0C()) * v0Efficiency));
                } else if (candMass > MassLambda0 + invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaPion_rightSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionPions, associate, coll_1.centFT0C()) / (trackEff(hEffPions, associate, coll_1.centFT0C()) * v0Efficiency));
                }

              } else if (assocPID[0] == kaonID) { // Kaons
                rMECorrRegistry.fill(HIST("hMixLambdaKaon_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, coll_1.centFT0C()) / (trackEff(hEffKaons, associate, coll_1.centFT0C()) * v0Efficiency));
                if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass < MassLambda0 - invMassNSigma * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaKaon_leftSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, coll_1.centFT0C()) / (trackEff(hEffKaons, associate, coll_1.centFT0C()) * v0Efficiency));
                } else if (candMass > MassLambda0 + invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaKaon_rightSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionKaons, associate, coll_1.centFT0C()) / (trackEff(hEffKaons, associate, coll_1.centFT0C()) * v0Efficiency));
                }

              } else if (assocPID[0] == protonID) { // Protons
                rMECorrRegistry.fill(HIST("hMixLambdaProton_SB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, coll_1.centFT0C()) / (trackEff(hEffProtons, associate, coll_1.centFT0C()) * v0Efficiency));
                if (candMass >= MassLambda0 - 2 * invMassNSigma * dGaussSigma && candMass < MassLambda0 - invMassNSigma * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaProton_leftSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, coll_1.centFT0C()) / (trackEff(hEffProtons, associate, coll_1.centFT0C()) * v0Efficiency));
                } else if (candMass > MassLambda0 + invMassNSigma * dGaussSigma && candMass <= MassLambda0 + 2 * invMassNSigma * dGaussSigma) {
                  rMECorrRegistry.fill(HIST("hMixLambdaProton_rightSB"), deltaPhi, deltaEta, coll_1.centFT0C(), coll_1.posZ(), triggSign, associate.sign(), ratioCorrection(hCorrectionProtons, associate, coll_1.centFT0C()) / (trackEff(hEffProtons, associate, coll_1.centFT0C()) * v0Efficiency));
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

    if (recCollisions.size() == 1 && acceptEvent(recCollisions.begin(), false)) {

      rQARegistry.fill(HIST("hEventCentrality_MC"), collision.bestCollisionCentFT0C());
      auto groupMCTriggers = mcTriggers->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
      auto groupMCAssociates = mcAssociates->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);

      // Start of the MC Track QA
      for (const auto& track : groupMCAssociates) {
        if (track.isPhysicalPrimary()) {

          if (track.pdgCode() > 0) {
            assocSign = 1;
          } else if (track.pdgCode() < 0) {
            assocSign = -1;
          }

          if (std::abs(track.pdgCode()) == kPiPlus) { // Pions
            rQARegistry.fill(HIST("hPtPion_MC"), track.pt(), collision.bestCollisionCentFT0C(), assocSign);
          } else if (std::abs(track.pdgCode()) == kKPlus) { // Kaons
            rQARegistry.fill(HIST("hPtKaon_MC"), track.pt(), collision.bestCollisionCentFT0C(), assocSign);
          } else if (std::abs(track.pdgCode()) == kProton) { // Protons
            rQARegistry.fill(HIST("hPtProton_MC"), track.pt(), collision.bestCollisionCentFT0C(), assocSign);
          }
        }
      }
      // End of the MC Track QA

      // Start of the MC Same-Event correlations
      for (const auto& trigger : groupMCTriggers) {
        if (trigger.isPhysicalPrimary()) {

          if (trigger.pdgCode() > 0) {
            triggSign = 1;
          } else if (trigger.pdgCode() < 0) {
            triggSign = -1;
          }
          rQARegistry.fill(HIST("hPtV0_MC"), trigger.pt(), collision.bestCollisionCentFT0C(), triggSign);
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
  }

  void processMCMixed(MyFilteredMCGenCollisions const&, MyFilteredMCParticles const&, MyFilteredMCRecCollisions const& recCollisions)
  {

    // Start of the MC Mixed-events Correlations
    for (const auto& [coll_1, v0_1, coll_2, particle_2] : pairMC) {
      auto recCollsA1 = recCollisions.sliceBy(perMCCol, coll_1.globalIndex());
      auto recCollsA2 = recCollisions.sliceBy(perMCCol, coll_2.globalIndex());
      if (recCollsA1.size() == 1 && recCollsA2.size() == 1 && acceptEvent(recCollsA1.begin(), false) && acceptEvent(recCollsA2.begin(), false)) {

        LOGF(info, "Size_1 = %i, Size_2 = %i", recCollsA1.size(), recCollsA1.size());
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
    }
    // End of the MC Mixed-events Correlations
  }

  void processMCGen(MyFilteredMCGenCollision const& collision, MyFilteredMCParticles const&, soa::SmallGroups<MCRecCollisions> const& recCollisions)
  {

    if (recCollisions.size() == 1 && acceptEvent(recCollisions.begin(), false)) {

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
            rMCRegistry.fill(HIST("hSelectPionP"), track.pt(), collision.centFT0C());
            if (particle.pdgCode() == kPiPlus) {
              rMCRegistry.fill(HIST("hTrueSelectPionP"), track.pt(), collision.centFT0C());
            }
          } else if (assocPID[0] == kaonID) { // Kaons
            rMCRegistry.fill(HIST("hSelectKaonP"), track.pt(), collision.centFT0C());
            if (particle.pdgCode() == kKPlus) {
              rMCRegistry.fill(HIST("hTrueSelectKaonP"), track.pt(), collision.centFT0C());
            }
          } else if (assocPID[0] == protonID) { // Protons
            rMCRegistry.fill(HIST("hSelectProtonP"), track.pt(), collision.centFT0C());
            if (particle.pdgCode() == kProton) {
              rMCRegistry.fill(HIST("hTrueSelectProtonP"), track.pt(), collision.centFT0C());
            }
          }
        } else if (track.sign() < 0) { // Negative tracks
          if (assocPID[0] == pionID) { // Pions
            rMCRegistry.fill(HIST("hSelectPionN"), track.pt(), collision.centFT0C());
            if (particle.pdgCode() == kPiMinus) {
              rMCRegistry.fill(HIST("hTrueSelectPionN"), track.pt(), collision.centFT0C());
            }
          } else if (assocPID[0] == kaonID) { // Kaons
            rMCRegistry.fill(HIST("hSelectKaonN"), track.pt(), collision.centFT0C());
            if (particle.pdgCode() == kKMinus) {
              rMCRegistry.fill(HIST("hTrueSelectKaonN"), track.pt(), collision.centFT0C());
            }
          } else if (assocPID[0] == protonID) { // Protons
            rMCRegistry.fill(HIST("hSelectProtonN"), track.pt(), collision.centFT0C());
            if (particle.pdgCode() == kProtonBar) {
              rMCRegistry.fill(HIST("hTrueSelectProtonN"), track.pt(), collision.centFT0C());
            }
          }
        }
      }
    }

    for (const auto& v0 : v0s) {

      if (!v0.has_mcParticle()) {
        continue;
      }

      if (v0Filters(collision, v0, tracks)) {

        v0Efficiency = v0Eff(hEffLambdas, v0, collision.centFT0C());

        // V0 efficiency - Reconstructed
        if (v0Sign(v0) == 1) { // Lambdas
          candMass = v0.mLambda();
          rQARegistry.fill(HIST("hInvMassLambda_MC"), v0.mLambda(), v0.pt(), collision.centFT0C(), 1. / v0Efficiency);
          rMCRegistry.fill(HIST("hRecLambdaP"), v0.pt(), v0.eta(), collision.centFT0C());
        } else if (v0Sign(v0) == -1) { // AntiLambdas
          candMass = v0.mAntiLambda();
          rQARegistry.fill(HIST("hInvMassAntiLambda_MC"), v0.mAntiLambda(), v0.pt(), collision.centFT0C(), 1. / v0Efficiency);
          rMCRegistry.fill(HIST("hRecLambdaN"), v0.pt(), v0.eta(), collision.centFT0C());
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

  template <class V0Cand>
  double v0Eff(TH3D** efficiencies, const V0Cand& v0, double centrality)
  {

    int index = -999;
    if (v0Sign(v0) > 0) {
      index = 0;
    } else if (v0Sign(v0) < 0) {
      index = 1;
    }

    double efficiency = efficiencies[index]->GetBinContent(efficiencies[index]->FindBin(v0.pt(), v0.eta(), centrality));
    if (efficiency > 0) {
      return efficiency;
    } else {
      return 1.0;
    }
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

  template <class TrackCand>
  double ratioCorrection(TH2D** ratios, const TrackCand& track, double centrality)
  {

    double ratioCorrection = 1.0;
    if (switchGroup.confRatioCorrectionSwitch) {

      int index = -999;
      if (track.sign() > 0) {
        index = 0;
      } else if (track.sign() < 0) {
        index = 1;
      }

      ratioCorrection = ratios[index]->GetBinContent(ratios[index]->FindBin(track.pt(), centrality));
    }

    return ratioCorrection;
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

    if (evSelGroup.useOccupCut) {
      if (!col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) { // kNoCollInTimeRangeStandard
        return false;
      }
      if (FillHist) {
        rQARegistry.fill(HIST("hNEvents"), 4.5);
      }
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
    if (posDaughter.tpcNClsCrossedRows() <= v0SelGroup.tpcNCrossedRows || negDaughter.tpcNClsCrossedRows() <= v0SelGroup.tpcNCrossedRows) {
      return false;
    }
    if (v0Sign(v0) == 1) {
      if (std::abs(posDaughter.tpcNSigmaPr()) >= nSigma5 || std::abs(negDaughter.tpcNSigmaPi()) >= nSigma5) {
        return false;
      }
      if (std::abs(v0.dcapostopv()) <= v0SelGroup.dcaProton || std::abs(v0.dcanegtopv()) <= v0SelGroup.dcaPion) {
        return false;
      }
    } else if (v0Sign(v0) == -1) {
      if (std::abs(posDaughter.tpcNSigmaPi()) >= nSigma5 || std::abs(negDaughter.tpcNSigmaPr()) >= nSigma5) {
        return false;
      }
      if (std::abs(v0.dcapostopv()) <= v0SelGroup.dcaPion || std::abs(v0.dcanegtopv()) <= v0SelGroup.dcaProton) {
        return false;
      }
    }

    // Topological cuts
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (v0.v0radius() <= v0SelGroup.decayR) {
      return false;
    }
    if (ctau >= v0SelGroup.ctau) {
      return false;
    }
    if (v0.v0cosPA() <= v0SelGroup.cosPA) {
      return false;
    }
    if (v0.dcaV0daughters() >= v0SelGroup.dcaV0Dau) {
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
      } else if (track.pt() > pionPtMin && track.pt() < pionPtMid1) {
        if (std::abs(track.tofNSigmaPi()) >= nSigma4) {
          return false;
        }
      } else if (track.pt() > pionPtMid1 && track.pt() < pionPtMid2) {
        if (track.tofNSigmaPi() <= -nSigma4 || track.tofNSigmaPi() >= nSigma2) {
          return false;
        }
      } else if (track.pt() > pionPtMid2 && track.pt() < pionPtMax) {
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
        if (track.tofNSigmaKa() <= -nSigma1 || track.tofNSigmaKa() >= nSigma4) {
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
        if (std::abs(track.tofNSigmaPr()) >= nSigma4) {
          return false;
        }
      } else if (track.pt() > protonPtMid) {
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

    if (switchGroup.confFakeV0Switch) {

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
      if (invMass >= MassLambda0 - invMassNSigma * dGaussSigma && invMass <= MassLambda0 + invMassNSigma * dGaussSigma) {
        return false;
      }
    }

    return true;
  }

  template <class V0Cand, class TrackCand>
  bool radialDistanceFilter(const V0Cand& v0, const TrackCand& track, double B, bool Mix)
  {

    bool pass = true;
    if (switchGroup.confRDSwitch) {

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

        if (r == rMin) {                              // TPC inner radius
          if (!Mix) {                                 // Same-event
            if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
              rPhiStarRegistry.fill(HIST("hSEPhiStarIR_OS"), dPhiStar, dEta);
            } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
              rPhiStarRegistry.fill(HIST("hSEPhiStarIR_SS"), dPhiStar, dEta);
              if (proton.sign() == 1) { // Positive
                rPhiStarRegistry.fill(HIST("hSEPhiStarIR_SSP"), dPhiStar, dEta);
              } else if (proton.sign() == -1) { // Negative
                rPhiStarRegistry.fill(HIST("hSEPhiStarIR_SSN"), dPhiStar, dEta);
              }
            }

          } else {                                    // Mixed-event
            if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
              rPhiStarRegistry.fill(HIST("hMEPhiStarIR_OS"), dPhiStar, dEta);
            } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
              rPhiStarRegistry.fill(HIST("hMEPhiStarIR_SS"), dPhiStar, dEta);
              if (proton.sign() == 1) { // Positive
                rPhiStarRegistry.fill(HIST("hMEPhiStarIR_SSP"), dPhiStar, dEta);
              } else if (proton.sign() == -1) { // Negative
                rPhiStarRegistry.fill(HIST("hMEPhiStarIR_SSN"), dPhiStar, dEta);
              }
            }
          }

          if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
            if (std::abs(dEta) < dEtaMin && std::abs(dPhiStar) < dPhiStarMinOS) {
              pass = false;
            }
          }
        }

        if (!Mix) {                                 // Same-event
          if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
            rPhiStarRegistry.fill(HIST("hSEPhiStarRadial_OS"), dPhiStar, r);
          } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
            rPhiStarRegistry.fill(HIST("hSEPhiStarRadial_SS"), dPhiStar, r);
          }

        } else {                                    // Mixed-event
          if (proton.sign() * track.sign() == -1) { // OS (Electric charge)
            rPhiStarRegistry.fill(HIST("hMEPhiStarRadial_OS"), dPhiStar, r);
          } else if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
            rPhiStarRegistry.fill(HIST("hMEPhiStarRadial_SS"), dPhiStar, r);
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

      if (proton.sign() * track.sign() == 1) { // SS (Electric charge)
        if (std::abs(dEta) < dEtaMin && std::abs(dPhiStarMean) < dPhiStarMinSS) {
          pass = false;
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
