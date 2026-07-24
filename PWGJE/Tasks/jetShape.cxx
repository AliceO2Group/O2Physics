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

/// \file   jetShape.cxx
/// \author Yuto Nishida <yuto.nishida@cern.ch>
/// \brief Task for measuring the dependence of the jet shape function rho(r) on the distance r from the jet axis.

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetShapeTask {

  static constexpr float PionMass = static_cast<float>(o2::constants::physics::MassPionCharged);
  static constexpr float ProtonMass = static_cast<float>(o2::constants::physics::MassProton);
  static constexpr float RandomConeNegativeSideProbability = 0.5f;

  Configurable<int> nBinsNSigma{"nBinsNSigma", 100, "Number of nsigma bins"};
  Configurable<float> nSigmaMin{"nSigmaMin", -10.0f, "Min value of nsigma"};
  Configurable<float> nSigmaMax{"nSigmaMax", 10.0f, "Max value of nsigma"};
  Configurable<int> nBinsPForDedx{"nBinsPForDedx", 700, "Number of p bins"};
  Configurable<int> nBinsPForBeta{"nBinsPForBeta", 500, "Number of pT bins"};
  Configurable<int> nBinsTpcDedx{"nBinsTpcDedx", 500, "Number of DEdx bins"};
  Configurable<int> nBinsTofBeta{"nBinsTofBeta", 350, "Number of Beta bins"};
  Configurable<float> dcaxyMin{"dcaxyMin", -1.0f, "Min value of dcaXY"};
  Configurable<float> dcaxyMax{"dcaxyMax", 1.0f, "Max value of dcaXY"};
  Configurable<float> pMax{"pMax", 8.0f, "Max value of p"};
  Configurable<float> ptMax{"ptMax", 6.0f, "Max value of pT"};
  Configurable<float> jetPtMinForCut{"jetPtMinForCut", 0.0f, "Minimum value of jet pT cut"};
  Configurable<float> jetPtMaxForCut{"jetPtMaxForCut", 200.0f, "Maximum value of the jet pT cut"};
  Configurable<float> centralityMinForCut{"centralityMinForCut", 0.0f, "Minimum value of ce cut"};
  Configurable<float> centralityMaxForCut{"centralityMaxForCut", 100.0f, "Maximum value of the jet pT cut"};
  Configurable<float> jetShapeFuncMax{"jetShapeFuncMax", 300, "Maximum value of JetShapeFunction"};
  Configurable<int> nBinsJetShapeFunc{"nBinsJetShapeFunc", 900, "Number of JetShapeFunction bins"};
  Configurable<int> nBinsDcaxyForData{"nBinsDcaxyForData", 400, "Number of DcaXY bins for data"};
  Configurable<int> nBinsDcaxyForMc{"nBinsDcaxyForMc", 400, "Number of DcaXY bins for mc data"};
  Configurable<int> nBinsP{"nBinsP", 40, "Number of p bins"};
  Configurable<int> nBinsPt{"nBinsPt", 30, "Number of pT bins"};
  Configurable<int> nBinsPtForDca{"nBinsPtForDca", 15, "Number of pT bins for dcaXY"};
  Configurable<int> nBinsJetPt{"nBinsJetPt", 10, "Number of jet pT bins"};
  Configurable<int> nBinsPForCut{"nBinsPForCut", 30, "Number of p track bins"};
  Configurable<int> nBinsCentrality{"nBinsCentrality", 10, "Number of centrality bins"};
  Configurable<int> nBinsDistance{"nBinsDistance", 7, "Number of distance bins"};
  Configurable<float> distanceMax{"distanceMax", 0.7f, "Max value of distance"};
  Configurable<float> nSigmaTofCut{"nSigmaTofCut", 2.0f, "Number of sigma cut for TOF PID"};
  Configurable<float> tpcNSigmaPrMin{"tpcNSigmaPrMin", -2.0f, "Min value of tpcNsigmaProton"};
  Configurable<float> tpcNSigmaPrMax{"tpcNSigmaPrMax", 2.0f, "Max value of tpcNsigmaProton"};
  Configurable<float> tpcNSigmaPiMin{"tpcNSigmaPiMin", -2.0f, "Min value of tpcNsigmaPion"};
  Configurable<float> tpcNSigmaPiMax{"tpcNSigmaPiMax", 2.0f, "Max value of tpcNsigmaPion"};
  Configurable<float> randomConeDeltaPhiMin{"randomConeDeltaPhiMin", static_cast<float>(o2::constants::math::PIThird), "Minimum delta phi used only when randomConeFullAzimuth=false"};
  Configurable<float> randomConeDeltaPhiMax{"randomConeDeltaPhiMax", static_cast<float>(2.0f * o2::constants::math::PIThird), "Maximum delta phi used only when randomConeFullAzimuth=false"};
  Configurable<bool> enableRandomCone{"enableRandomCone", true, "Fill one accepted random cone per selected jet"};
  Configurable<bool> randomConeFullAzimuth{"randomConeFullAzimuth", true, "Generate random cones over the full azimuth instead of the old perpendicular-side window"};
  Configurable<int> randomConeMaxAttempts{"randomConeMaxAttempts", 100, "Maximum attempts for each random cone"};
  Configurable<unsigned int> randomConeSeed{"randomConeSeed", 0u, "TRandom3 seed; 0 lets ROOT choose a unique seed"};
  Configurable<bool> fillDetailedQA{"fillDetailedQA", false, "Fill high-volume dE/dx, beta and distance QA histograms"};
  Configurable<bool> applyCentralityCut{"applyCentralityCut", true, "Apply centralityMinForCut <= centrality < centralityMaxForCut inside this task"};
  Configurable<bool> useRapiditySelection{"useRapiditySelection", false, "Use species-dependent |y| cut in addition to the detector |eta| acceptance (cross-check only)"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8",
                                            "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "nominal track selection for yields and tracking efficiency"};
  Configurable<std::string> trackSelectionsForDca{"trackSelectionsForDca", "QualityTracksWDCA", "wide-DCA track selection used only for data DCA distributions and MC DCA templates"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", 5.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};

  // for jet shape
  Configurable<std::vector<float>> distanceCategory{"distanceCategory", {0.00f, 0.05f, 0.10f, 0.15f, 0.20f, 0.25f, 0.30f, 0.35f, 0.40f, 0.45f, 0.50f, 0.55f, 0.60f, 0.65f, 0.70f}, "distance of category"};

  // for ppi production
  Configurable<float> etaTrUp{"etaTrUp", 0.9f, "maximum track eta"};
  Configurable<float> trackPtMinForAnalysis{"trackPtMinForAnalysis", 0.15f, "minimum reconstructed and generated particle pT used for yields and efficiencies"};
  Configurable<float> dcaxyCutMax{"dcaxyCutMax", 2.0f, "maximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 4.0f, "tpc chi2 cut"};
  Configurable<float> nclItsMin{"nclItsMin", 2.0f, "its # of cluster cut"};
  Configurable<float> nclTpcMin{"nclTpcMin", 70.0f, "tpc # if cluster cut"};
  Configurable<float> nclcrossTpcMin{"nclcrossTpcMin", 70.0f, "tpc # of crossedRows cut"};
  Configurable<float> mcRapidityMax{"mcRapidityMax", 0.5f, "maximum |y| used only when useRapiditySelection=true"};
  Configurable<double> epsilon{"epsilon", 1e-6, "standard for aboid division of zero"};
  Configurable<float> maxDeltaEtaSafe{"maxDeltaEtaSafe", 0.9f, "maximum track eta for cut"};
  Configurable<float> nSigmaMaxForDcaxy{"nSigmaMaxForDcaxy", 4.0f, "maximum nSigma for DCAxy"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};

  HistogramRegistry registry{"registry"};
  TRandom3 randomNumber;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  int trackSelectionForDca = -1;
  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    trackSelectionForDca = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelectionsForDca));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    randomNumber.SetSeed(randomConeSeed);

    // Histograms definition
    registry.add("tpcTofPi", "tpcTofPi", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcTofPr", "tpcTofPr", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcPi", "tpcPi", HistType::kTH2F, {{nBinsP.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("tofPi", "tofPi", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("tpcPr", "tpcPr", HistType::kTH2F, {{nBinsP.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("tofPr", "tofPr", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("tpcDedx", "tpcDedx", HistType::kTHnSparseD, {{nBinsPForDedx.value, 0, pMax.value}, {nBinsTpcDedx.value, 0, 1000}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tofBeta", "tofBeta", HistType::kTHnSparseD, {{nBinsPForBeta.value, 0, pMax.value}, {nBinsTofBeta.value, 0.4, 1.1}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPr", "pVsPtForPr", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPi", "pVsPtPi", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPiInc", "inclusive pion-candidate DCAxy from wide-DCA tracks", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrInc", "inclusive proton-candidate DCAxy from wide-DCA tracks", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("trackPhi", "trackPhi", HistType::kTH1F, {{80, -1, 7}});
    registry.add("trackEta", "trackEta", HistType::kTH1F, {{100, -1, 1}});
    registry.add("trackTpcNClsCrossedRows", "trackTpcNClsCrossedRows", HistType::kTH1F, {{50, 0, 200}});
    registry.add("trackDcaXY", "trackDcaXY", HistType::kTH1F, {{40, -10, 10}});
    registry.add("trackItsChi2NCl", "trackItsChi2NCl", HistType::kTH1F, {{60, 0, 30}});
    registry.add("trackTpcChi2NCl", "trackTpcChi2NCl", HistType::kTH1F, {{100, 0, 50}});
    registry.add("trackTpcNClsFound", "trackTpcNClsFound", HistType::kTH1F, {{100, 0, 200}});
    registry.add("trackItsNCls", "trackItsNCls", HistType::kTH1F, {{10, 0, 10}});

    registry.add("jetTpcPi", "jetTpcPi", HistType::kTH2F, {{nBinsP.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("jetTofPi", "jetTofPi", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("jetTpcPr", "jetTpcPr", HistType::kTH2F, {{nBinsP.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("jetTofPr", "jetTofPr", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});

    registry.add("jetTpcDedx", "jetTpcDedx", HistType::kTHnSparseD, {{nBinsPForDedx.value, 0, pMax.value}, {nBinsTpcDedx.value, 0, 1000}, {nBinsDistance.value, 0, distanceMax.value}});
    registry.add("tpcDedxPerpJet", "tpcDedxPerpJet", HistType::kTH2F, {{nBinsPForDedx.value, 0, pMax.value}, {nBinsTpcDedx.value, 0, 1000}});
    registry.add("jetTofBeta", "jetTofBeta", HistType::kTH2F, {{nBinsPForBeta.value, 0, pMax.value}, {nBinsTofBeta.value, 0.4, 1.1}});

    // Jet (in Cone)
    registry.add("jetTpcTofPi", "jetTpcTofPi", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetTpcTofPr", "jetTpcTofPr", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetpVsPtForPr", "jetpVsPtForPr", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetpVsPtForPi", "jetpVsPtPi", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetDcaPr", "jetDcaPr", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetDcaPi", "jetDcaPi", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Perp Cone
    registry.add("tpcTofPiPerpJet", "tpcTofPiPerpJet", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcTofPrPerpJet", "tpcTofPrPerpJet", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPrPerpJet", "pVsPtForPrPerpJet", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPiPerpJet", "pVsPtPionPerpJet", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrPerpJet", "dcaPrPerpJet", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPiPerpJet", "dcaPiPerpJet", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Random Cone
    registry.add("tpcTofPiRandCone", "tpcTofPiRandCone", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcTofPrRandCone", "tpcTofPrRandCone", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPrRandCone", "pVsPtForPrRandCone", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPiRandCone", "pVsPtPionRandCone", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrRandCone", "dcaPrRandCone", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPiRandCone", "dcaPiRandCone", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("jetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", HistType::kTH2F, {{200, 0., 200.}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetPtMc", "MC jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", HistType::kTH2F, {{200, 0., 200.}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetEta", "jet #eta;#eta_{jet};entries", HistType::kTH1F, {{100, -1.0, 1.0}});
    registry.add("jetPhi", "jet #phi;#phi_{jet};entries", HistType::kTH1F, {{80, -1.0, 7.}});
    registry.add("area", "area", HistType::kTH2F, {{100, 0, 2}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("rho", "rho", HistType::kTH2F, {{120, 0, 300}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptCorr", "Corrected jet pT; p_{T}^{corr} (GeV/c); Counts", HistType::kTH1F, {{200, 0, 200}});
    registry.add("jetPtCorrCent", "Selected jets;p_{T,jet}^{corr} (GeV/c);Centrality (%)", HistType::kTH2F, {{200, 0, 200}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("areaExposureJet", "Sum of accepted annulus area for signal;r;p_{T,jet}^{corr};Centrality", HistType::kTHnSparseD, {{nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("areaExposurePerp", "Sum of accepted annulus area for two perpendicular cones;r;p_{T,jet}^{corr};Centrality", HistType::kTHnSparseD, {{nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("areaExposureRandCone", "Sum of accepted annulus area for random cones;r;p_{T,jet}^{corr};Centrality", HistType::kTHnSparseD, {{nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptCorrVsDistance", "ptcorr_vs_distance", HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}});
    registry.add("jetDistanceVsTrackpt", "trackpt_vs_distance_injet", HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}});
    registry.add("ptSum", "ptSum", HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc.value, 0, jetShapeFuncMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptSumBg1", "ptSumBg1", HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc.value, 0, jetShapeFuncMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptSumBg2", "ptSumBg2", HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc.value, 0, jetShapeFuncMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("event/vertexz", ";Vtx_{z} (cm);Entries", HistType::kTH1F, {{100, -20, 20}});
    registry.add("eventCounterJetShape", "eventCounterJetShape", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterJet", "eventCounterJet", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterInc", "eventCounterInc", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterRandCone", "Number of Random Cones;Centrality (%);Count", HistType::kTH2F, {{nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterMc", "eventCounterMc", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("ptVsCentrality", "ptvscentrality", HistType::kTH2F, {{100, 0, 100}, {300, 0, 300}});
    registry.add("ptResolution", "ptResolution", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {100, -1.0, +1.0}});
    registry.add("mcCentralityReco", "mcCentralityReco", HistType::kTH1F, {{100, 0, 100}});

    // MC Efficiency Denominators
    // Jet (In-cone)
    registry.add("ptGenPi", "ptGenPi", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenKa", "ptGenKa", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenPr", "ptGenPr", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    // Perp Cone
    registry.add("ptGenPiPerp", "ptGenPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenKaPerp", "ptGenKaPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenPrPerp", "ptGenPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    // Random cone
    registry.add("ptGenPiRand", "ptGenPiRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenKaRand", "ptGenKaRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenPrRand", "ptGenPrRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    // Inclusive (No JetPt axis)
    registry.add("ptGenPiInc", "ptGenPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenKaInc", "ptGenKaInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptGenPrInc", "ptGenPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // MC Efficiency Numerators
    // Inclusive
    registry.add("effNumTpcPiInc", "effNumTpcPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcKaInc", "effNumTpcKaInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcPrInc", "effNumTpcPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPiInc", "effNumTofPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofKaInc", "effNumTofKaInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPrInc", "effNumTofPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPiInc", "True primary pion reconstructed with the full analysis PID selection", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPrInc", "True primary proton reconstructed with the full analysis PID selection", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Jet (In-cone)
    registry.add("effNumTpcPiJet", "effNumTpcPiJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcKaJet", "effNumTpcKaJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcPrJet", "effNumTpcPrJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPiJet", "effNumTofPiJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofKaJet", "effNumTofKaJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPrJet", "effNumTofPrJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPiJet", "Full pion PID efficiency numerator in signal region", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPrJet", "Full proton PID efficiency numerator in signal region", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Perp Cone
    registry.add("effNumTpcPiPerp", "effNumTpcPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcKaPerp", "effNumTpcKaPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcPrPerp", "effNumTpcPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPiPerp", "effNumTofPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofKaPerp", "effNumTofKaPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPrPerp", "effNumTofPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPiPerp", "Full pion PID efficiency numerator in perpendicular cones", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPrPerp", "Full proton PID efficiency numerator in perpendicular cones", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Random cone efficiency numerators
    registry.add("effNumTpcPiRand", "effNumTpcPiRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcKaRand", "effNumTpcKaRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcPrRand", "effNumTpcPrRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPiRand", "effNumTofPiRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofKaRand", "effNumTofKaRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPrRand", "effNumTofPrRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPiRand", "Full pion PID efficiency numerator in random cones", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumPidPrRand", "Full proton PID efficiency numerator in random cones", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // DCA Templates for Primary Fraction
    // Jet (In-cone)
    registry.add("dcaPrimPi", "dcaPrimPi", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPi", "dcaDecayPi", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPi", "dcaMatPi", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrimPr", "dcaPrimPr", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPr", "dcaDecayPr", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPr", "dcaMatPr", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Perp Cone
    registry.add("dcaPrimPiPerp", "dcaPrimPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPiPerp", "dcaDecayPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPiPerp", "dcaMatPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrimPrPerp", "dcaPrimPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPrPerp", "dcaDecayPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPrPerp", "dcaMatPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Random cone
    registry.add("dcaPrimPiRand", "dcaPrimPiRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPiRand", "dcaDecayPiRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPiRand", "dcaMatPiRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrimPrRand", "dcaPrimPrRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPrRand", "dcaDecayPrRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPrRand", "dcaMatPrRand", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Inclusive
    registry.add("dcaPrimPiInc", "dcaPrimPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPiInc", "dcaDecayPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPiInc", "dcaMatPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrimPrInc", "dcaPrimPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPrInc", "dcaDecayPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPrInc", "dcaMatPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Transfer from the wide-DCA fit sample to the nominal global-track sample.
    // Axes:
    //   Jet/Perp: reco pT, reco distance, detector-level corrected jet pT,
    //             centrality, species (0=pi,1=p), origin (0=primary,1=decay,2=material)
    //   Inclusive: reco pT, centrality, species, origin
    registry.add("originCountWideJet", "wide-DCA origin counts around matched detector-level jets", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountGlobalJet", "nominal-global origin counts around matched detector-level jets", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountWidePerp", "wide-DCA origin counts in perpendicular cones", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountGlobalPerp", "nominal-global origin counts in perpendicular cones", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountWideRand", "wide-DCA origin counts in random cones", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountGlobalRand", "nominal-global origin counts in random cones", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountWideInc", "inclusive wide-DCA origin counts", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
    registry.add("originCountGlobalInc", "inclusive nominal-global origin counts", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}, {2, 0, 2}, {3, 0, 3}});
  }

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {
    static constexpr double JetAreaFractionMinValue = -98.0;
    if (jetAreaFractionMin > JetAreaFractionMinValue) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI *
                         (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
      if (jet.area() <
          o2::constants::math::PIHalf * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    static constexpr double LeadingConstituentPtMaxValue = 9998.0;
    bool checkConstituentPt = true;
    // A non-negative value activates the minimum-leading-constituent cut.
    // Therefore the default value 5 GeV/c is actually applied.
    bool checkConstituentMinPt = (leadingConstituentPtMin >= 0.0f);
    bool checkConstituentMaxPt =
      (leadingConstituentPtMax < LeadingConstituentPtMaxValue);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<T>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }

    return true;
  }

  template <typename T>
  bool passCommonTrackQuality(T const& track)
  {
    return track.pt() >= trackPtMinForAnalysis &&
           std::abs(track.eta()) <= etaTrUp &&
           track.tpcNClsCrossedRows() >= nclcrossTpcMin &&
           track.itsChi2NCl() <= chi2ItsMax &&
           track.tpcChi2NCl() <= chi2TpcMax &&
           track.tpcNClsFound() >= nclTpcMin &&
           track.itsNCls() >= nclItsMin;
  }

  template <typename T>
  bool passNominalTrackDca(T const& track)
  {
    return std::abs(track.dcaXY()) <= dcaxyCutMax;
  }

  bool passCentrality(float centrality) const
  {
    return !applyCentralityCut.value ||
           (centrality >= centralityMinForCut.value &&
            centrality < centralityMaxForCut.value);
  }

  static float absDeltaPhi(float phi1, float phi2)
  {
    float dPhi = std::abs(phi1 - phi2);
    if (dPhi > o2::constants::math::PI) {
      dPhi = o2::constants::math::TwoPI - dPhi;
    }
    return dPhi;
  }

  static float rapidityFromPtEtaMass(float pt, float eta, float mass)
  {
    const float pz = pt * std::sinh(eta);
    const float mt = std::sqrt(pt * pt + mass * mass);
    return std::asinh(pz / mt);
  }

  bool passSpeciesKinematics(float pt, float eta, float mass) const
  {
    if (std::abs(eta) >= etaTrUp.value) {
      return false;
    }

    if (!useRapiditySelection.value) {
      return true;
    }

    return std::abs(rapidityFromPtEtaMass(pt, eta, mass)) <
           mcRapidityMax.value;
  }

  template <typename T>
  bool passTruthKinematics(T const& particle) const
  {
    if (particle.pt() < trackPtMinForAnalysis.value ||
        std::abs(particle.eta()) >= etaTrUp.value) {
      return false;
    }

    return !useRapiditySelection.value ||
           std::abs(particle.y()) < mcRapidityMax.value;
  }

  static double diskAreaPrimitive(double radius, double x)
  {
    if (radius <= 0.0) {
      return 0.0;
    }
    x = std::clamp(x, -radius, radius);
    const double root = std::sqrt(std::max(0.0, radius * radius - x * x));
    return x * root + radius * radius * std::asin(x / radius);
  }

  static double effectiveDiskArea(double radius, double axisEta, double etaMin, double etaMax)
  {
    const double lower = std::max(-radius, etaMin - axisEta);
    const double upper = std::min(radius, etaMax - axisEta);
    if (upper <= lower) {
      return 0.0;
    }
    return diskAreaPrimitive(radius, upper) - diskAreaPrimitive(radius, lower);
  }

  static double effectiveAnnulusArea(double rMin, double rMax, double axisEta, double etaMin, double etaMax)
  {
    return effectiveDiskArea(rMax, axisEta, etaMin, etaMax) -
           effectiveDiskArea(rMin, axisEta, etaMin, etaMax);
  }

  struct ConeVetoJet {
    float eta;
    float phi;
    float ptCorr;
  };

  struct RandomConeAxis {
    float eta;
    float phi;
    float referenceJetPtCorr;
  };

  std::vector<RandomConeAxis> makeRandomCones(std::vector<ConeVetoJet> const& jets)
  {
    std::vector<RandomConeAxis> cones;
    if (!enableRandomCone || jets.empty() || etaTrUp <= distanceMax) {
      return cones;
    }

    cones.reserve(jets.size());
    const float etaMin = -etaTrUp + distanceMax;
    const float etaMax = etaTrUp - distanceMax;
    const float vetoDistance = jetR + distanceMax;

    for (const auto& referenceJet : jets) {
      bool accepted = false;
      float rcEta = 0.0f;
      float rcPhi = 0.0f;

      for (int attempt = 0; attempt < randomConeMaxAttempts; ++attempt) {
        rcEta = randomNumber.Uniform(etaMin, etaMax);
        if (randomConeFullAzimuth) {
          rcPhi = randomNumber.Uniform(0.0f, static_cast<float>(o2::constants::math::TwoPI));
        } else {
          float dPhi = randomNumber.Uniform(randomConeDeltaPhiMin, randomConeDeltaPhiMax);
          if (randomNumber.Uniform() < RandomConeNegativeSideProbability) {
            dPhi = -dPhi;
          }
          rcPhi = RecoDecay::constrainAngle(referenceJet.phi + dPhi);
        }

        accepted = true;
        for (const auto& vetoJet : jets) {
          const float dEta = rcEta - vetoJet.eta;
          const float dPhi = absDeltaPhi(rcPhi, vetoJet.phi);
          if (dEta * dEta + dPhi * dPhi <= vetoDistance * vetoDistance) {
            accepted = false;
            break;
          }
        }
        if (accepted) {
          break;
        }
      }

      if (accepted) {
        cones.push_back({rcEta, rcPhi, referenceJet.ptCorr});
      }
    }
    return cones;
  }

  enum class ExposureRegion { Signal,
                              Perpendicular };

  void fillAreaExposure(float axisEta, float jetPtCorr, float centrality, ExposureRegion region, double multiplicity = 1.0)
  {
    const double dr = distanceMax / static_cast<double>(nBinsDistance);
    for (int i = 0; i < nBinsDistance; ++i) {
      const double rMin = i * dr;
      const double rMax = (i + 1) * dr;
      const double rCenter = 0.5 * (rMin + rMax);
      const double area = multiplicity * effectiveAnnulusArea(rMin, rMax, axisEta, -etaTrUp, etaTrUp);
      if (region == ExposureRegion::Signal) {
        registry.fill(HIST("areaExposureJet"), rCenter, jetPtCorr, centrality, area);
      } else {
        registry.fill(HIST("areaExposurePerp"), rCenter, jetPtCorr, centrality, area);
      }
    }
  }

  void fillFullRandomConeArea(float jetPtCorr, float centrality)
  {
    const double dr = distanceMax / static_cast<double>(nBinsDistance);
    for (int i = 0; i < nBinsDistance; ++i) {
      const double rMin = i * dr;
      const double rMax = (i + 1) * dr;
      const double rCenter = 0.5 * (rMin + rMax);
      const double area = o2::constants::math::PI * (rMax * rMax - rMin * rMin);
      registry.fill(HIST("areaExposureRandCone"), rCenter, jetPtCorr, centrality, area);
    }
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  using JetTracksWithPI = soa::Join<aod::JetTracks, aod::JTrackPIs>;
  using JetTracksWithPIMC = soa::Join<aod::JetTracks, aod::JTrackPIs, aod::JMcTrackLbs>;

  using FullTrackInfo = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta>;
  using FullTrackInfoMC = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::pidTOFbeta>;

  void processJetShape(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, aod::JetTracks const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    const float cent = collision.centFT0M();
    if (!passCentrality(cent)) {
      return;
    }
    registry.fill(HIST("eventCounterJetShape"), cent);

    size_t nBins = distanceCategory->size() - 1;

    float maxDistance = distanceCategory->at(nBins);

    struct JetAccumulator {
      float pt, eta, phi, ptCorr;
      float phiBg1, phiBg2;
      float area;
    };

    std::vector<JetAccumulator> cachedJets;
    std::vector<std::vector<float>> accumSumSig;
    std::vector<std::vector<float>> accumSumBg1;
    std::vector<std::vector<float>> accumSumBg2;

    cachedJets.reserve(jets.size());
    accumSumSig.reserve(jets.size());
    accumSumBg1.reserve(jets.size());
    accumSumBg2.reserve(jets.size());

    float rho = collision.rho();

    for (const auto& jet : jets) {
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }

      float ptCorr = jet.pt() - rho * jet.area();

      float phiBg1 =
        RecoDecay::constrainAngle(jet.phi() + o2::constants::math::PIHalf);
      float phiBg2 =
        RecoDecay::constrainAngle(jet.phi() - o2::constants::math::PIHalf);

      // make accumulator
      cachedJets.push_back({jet.pt(), jet.eta(), jet.phi(), ptCorr, phiBg1, phiBg2, jet.area()});

      accumSumSig.push_back(std::vector<float>(nBins, 0.f));
      accumSumBg1.push_back(std::vector<float>(nBins, 0.f));
      accumSumBg2.push_back(std::vector<float>(nBins, 0.f));
    }

    if (cachedJets.empty()) {
      return;
    }

    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      float trkPt = track.pt();
      float trkPhi = track.phi();
      float trkEta = track.eta();

      for (size_t i = 0; i < cachedJets.size(); ++i) {

        const auto& jet = cachedJets[i];

        float dEta = trkEta - jet.eta;

        if (std::abs(dEta) > maxDistance) {
          continue;
        }

        float dPhi = RecoDecay::constrainAngle(trkPhi - jet.phi);
        float distance = std::sqrt(dEta * dEta + dPhi * dPhi);

        registry.fill(HIST("ptCorrVsDistance"), distance, jet.ptCorr);

        if (distance < maxDistance) {
          registry.fill(HIST("ptVsCentrality"), cent, trkPt);
        }

        // summmention in Signal resion
        for (size_t k = 0; k < nBins; k++) {
          if (distanceCategory->at(k) <= distance &&
              distance < distanceCategory->at(k + 1)) {
            accumSumSig[i][k] += trkPt;
            break;
          }
        }

        // calculation in Background region
        // 1. Bg1
        float deltaPhiBg1 = RecoDecay::constrainAngle(trkPhi - jet.phiBg1);
        float distBg1 = std::sqrt(dEta * dEta + deltaPhiBg1 * deltaPhiBg1);

        for (size_t k = 0; k < nBins; k++) {
          if (distanceCategory->at(k) <= distBg1 &&
              distBg1 < distanceCategory->at(k + 1)) {
            accumSumBg1[i][k] += trkPt;
            break;
          }
        }

        // 2. Bg2
        float deltaPhiBg2 = RecoDecay::constrainAngle(trkPhi - jet.phiBg2);
        float distBg2 = std::sqrt(dEta * dEta + deltaPhiBg2 * deltaPhiBg2);

        for (size_t k = 0; k < nBins; k++) {
          if (distanceCategory->at(k) <= distBg2 &&
              distBg2 < distanceCategory->at(k + 1)) {
            accumSumBg2[i][k] += trkPt;
            break;
          }
        }
      }
    }

    for (size_t i = 0; i < cachedJets.size(); ++i) {
      const auto& jet = cachedJets[i];

      registry.fill(HIST("ptCorr"), jet.ptCorr);

      for (size_t k = 0; k < nBins; k++) {
        float binWidth = distanceCategory->at(k + 1) - distanceCategory->at(k);
        double jetX = distanceCategory->at(k) + binWidth / 2.0;

        // Sum(pT) / (binWidth * ptCorr
        float normFactor = (binWidth * jet.ptCorr);

        // to avoid division of zero
        double jetShapeFunction = 0.0;
        double jetShapeFunctionBg1 = 0.0;
        double jetShapeFunctionBg2 = 0.0;

        if (std::abs(normFactor) > epsilon) {
          jetShapeFunction = accumSumSig[i][k] / normFactor;
          jetShapeFunctionBg1 = accumSumBg1[i][k] / normFactor;
          jetShapeFunctionBg2 = accumSumBg2[i][k] / normFactor;
        }

        registry.fill(HIST("ptSum"), jetX, jetShapeFunction, jet.ptCorr, cent);
        registry.fill(HIST("ptSumBg1"), jetX, jetShapeFunctionBg1, jet.ptCorr, cent);
        registry.fill(HIST("ptSumBg2"), jetX, jetShapeFunctionBg2, jet.ptCorr, cent);
      }
    }
  }

  PROCESS_SWITCH(JetShapeTask, processJetShape, "JetShape", false);

  template <typename CollisionT>
  void runJetProductionRatio(CollisionT const& collision,
                             JetTracksWithPI const& tracks,
                             soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                             float rho)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    const float centrality = collision.centFT0M();
    if (!passCentrality(centrality)) {
      return;
    }

    registry.fill(HIST("eventCounterJet"), centrality);
    registry.fill(HIST("event/vertexz"), collision.posZ());
    registry.fill(HIST("rho"), rho, centrality);

    struct CachedJet {
      float pt, eta, phi, ptCorr;
      float phiBg1, phiBg2;
    };
    std::vector<CachedJet> cachedJets;
    std::vector<ConeVetoJet> vetoJets;
    cachedJets.reserve(jets.size());
    vetoJets.reserve(jets.size());

    for (const auto& jet : jets) {
      if (!isAcceptedJet<JetTracksWithPI>(jet)) {
        continue;
      }

      const float ptCorr = jet.pt() - rho * jet.area();
      const float phiBg1 = RecoDecay::constrainAngle(jet.phi() + o2::constants::math::PIHalf);
      const float phiBg2 = RecoDecay::constrainAngle(jet.phi() - o2::constants::math::PIHalf);

      registry.fill(HIST("jetEta"), jet.eta());
      registry.fill(HIST("jetPhi"), jet.phi());
      registry.fill(HIST("area"), jet.area(), centrality);
      registry.fill(HIST("jetPt"), jet.pt(), centrality);
      registry.fill(HIST("ptCorr"), ptCorr);
      registry.fill(HIST("jetPtCorrCent"), ptCorr, centrality);

      fillAreaExposure(jet.eta(), ptCorr, centrality, ExposureRegion::Signal);
      fillAreaExposure(jet.eta(), ptCorr, centrality, ExposureRegion::Perpendicular, 2.0);

      cachedJets.push_back({jet.pt(), jet.eta(), jet.phi(), ptCorr, phiBg1, phiBg2});
      vetoJets.push_back({jet.eta(), jet.phi(), ptCorr});
    }

    if (cachedJets.empty()) {
      return;
    }

    const auto randomCones = makeRandomCones(vetoJets);
    for (const auto& cone : randomCones) {
      registry.fill(HIST("eventCounterRandCone"), cone.referenceJetPtCorr, centrality);
      fillFullRandomConeArea(cone.referenceJetPtCorr, centrality);
    }

    const float maxR2 = distanceMax * distanceMax;

    // The jet finder remains configured with globalTracks. The nominal-yield
    // sample and the wide-DCA sample are selected independently relative to
    // the same reconstructed jet axes.
    for (const auto& jetTrack : tracks) {
      const bool passesNominalBit = jetderiveddatautilities::selectTrack(jetTrack, trackSelection);
      const bool passesWideDcaBit = jetderiveddatautilities::selectTrack(jetTrack, trackSelectionForDca);
      if (!passesNominalBit && !passesWideDcaBit) {
        continue;
      }

      const auto track = jetTrack.track_as<FullTrackInfo>();
      if (!passCommonTrackQuality(track)) {
        continue;
      }

      const bool passesNominal = passesNominalBit && passNominalTrackDca(track);
      const bool passesWideDca = passesWideDcaBit;
      if (!passesNominal && !passesWideDca) {
        continue;
      }

      const float trkPt = track.pt();
      const float trkP = track.p();
      const float trkEta = track.eta();
      const float trkPhi = track.phi();

      const float tpcSig = track.tpcSignal();
      const float tpcPi = track.tpcNSigmaPi();
      const float tofPi = track.tofNSigmaPi();
      const float tpcPr = track.tpcNSigmaPr();
      const float tofPr = track.tofNSigmaPr();
      const float beta = track.beta();

      const bool passPiKinematics = passSpeciesKinematics(trkPt, trkEta, PionMass);
      const bool passPrKinematics = passSpeciesKinematics(trkPt, trkEta, ProtonMass);
      const bool hasTofPi = passPiKinematics && track.hasTOF() && (std::abs(tofPi) < nSigmaTofCut);
      const bool hasTofPr = passPrKinematics && track.hasTOF() && (std::abs(tofPr) < nSigmaTofCut);
      const bool isTpcPiRange = (tpcPi > tpcNSigmaPiMin && tpcPi < tpcNSigmaPiMax);
      const bool isTpcPrRange = (tpcPr > tpcNSigmaPrMin && tpcPr < tpcNSigmaPrMax);

      if (passesNominal && fillDetailedQA) {
        registry.fill(HIST("jetTpcPi"), trkP, tpcPi);
        registry.fill(HIST("jetTofPi"), trkPt, tofPi);
        registry.fill(HIST("jetTpcPr"), trkP, tpcPr);
        registry.fill(HIST("jetTofPr"), trkPt, tofPr);
      }

      // Signal and perpendicular cones: one loop over accepted jets.
      for (const auto& jet : cachedJets) {
        const float dEta = trkEta - jet.eta;
        if (std::abs(dEta) >= distanceMax) {
          continue;
        }

        const float dPhi = absDeltaPhi(trkPhi, jet.phi);
        const float distanceSq = dEta * dEta + dPhi * dPhi;

        auto fillPerpendicular = [&](float conePhi) {
          const float dPhiBg = absDeltaPhi(trkPhi, conePhi);
          const float distBgSq = dEta * dEta + dPhiBg * dPhiBg;
          if (distBgSq >= maxR2) {
            return;
          }
          const float distBg = std::sqrt(distBgSq);

          if (passesWideDca) {
            if (hasTofPr && isTpcPrRange) {
              registry.fill(HIST("dcaPrPerpJet"), trkPt, track.dcaXY(), distBg, jet.ptCorr, centrality);
            }
            if (hasTofPi && isTpcPiRange) {
              registry.fill(HIST("dcaPiPerpJet"), trkPt, track.dcaXY(), distBg, jet.ptCorr, centrality);
            }
          }

          if (passesNominal) {
            if (fillDetailedQA) {
              registry.fill(HIST("tpcDedxPerpJet"), trkP, tpcSig);
            }
            if (hasTofPi) {
              registry.fill(HIST("tpcTofPiPerpJet"), trkP, tpcPi, distBg, jet.ptCorr, centrality);
              if (isTpcPiRange) {
                registry.fill(HIST("pVsPtForPiPerpJet"), trkP, trkPt, distBg, jet.ptCorr, centrality);
              }
            }
            if (hasTofPr) {
              registry.fill(HIST("tpcTofPrPerpJet"), trkP, tpcPr, distBg, jet.ptCorr, centrality);
              if (isTpcPrRange) {
                registry.fill(HIST("pVsPtForPrPerpJet"), trkP, trkPt, distBg, jet.ptCorr, centrality);
              }
            }
          }
        };

        fillPerpendicular(jet.phiBg1);
        fillPerpendicular(jet.phiBg2);

        if (distanceSq >= maxR2) {
          continue;
        }
        const float distance = std::sqrt(distanceSq);

        if (passesWideDca) {
          if (hasTofPr && isTpcPrRange) {
            registry.fill(HIST("jetDcaPr"), trkPt, track.dcaXY(), distance, jet.ptCorr, centrality);
          }
          if (hasTofPi && isTpcPiRange) {
            registry.fill(HIST("jetDcaPi"), trkPt, track.dcaXY(), distance, jet.ptCorr, centrality);
          }
        }

        if (passesNominal) {
          if (fillDetailedQA) {
            registry.fill(HIST("jetDistanceVsTrackpt"), distance, trkPt);
            registry.fill(HIST("jetTpcDedx"), trkP, tpcSig, distance);
            registry.fill(HIST("jetTofBeta"), trkP, beta);
          }

          if (hasTofPr) {
            registry.fill(HIST("jetTpcTofPr"), trkP, tpcPr, distance, jet.ptCorr, centrality);
            if (isTpcPrRange) {
              registry.fill(HIST("jetpVsPtForPr"), trkP, trkPt, distance, jet.ptCorr, centrality);
            }
          }
          if (hasTofPi) {
            registry.fill(HIST("jetTpcTofPi"), trkP, tpcPi, distance, jet.ptCorr, centrality);
            if (isTpcPiRange) {
              registry.fill(HIST("jetpVsPtForPi"), trkP, trkPt, distance, jet.ptCorr, centrality);
            }
          }
        }
      }

      // Random cones are outside the jet loop, so every track-cone pair is
      // filled exactly once. Each cone is associated with one selected jet.
      for (const auto& cone : randomCones) {
        const float dEtaRC = trkEta - cone.eta;
        if (std::abs(dEtaRC) >= distanceMax) {
          continue;
        }
        const float dPhiRC = absDeltaPhi(trkPhi, cone.phi);
        const float distRCSq = dEtaRC * dEtaRC + dPhiRC * dPhiRC;
        if (distRCSq >= maxR2) {
          continue;
        }
        const float distRC = std::sqrt(distRCSq);

        if (passesWideDca) {
          if (hasTofPr && isTpcPrRange) {
            registry.fill(HIST("dcaPrRandCone"), trkPt, track.dcaXY(), distRC, cone.referenceJetPtCorr, centrality);
          }
          if (hasTofPi && isTpcPiRange) {
            registry.fill(HIST("dcaPiRandCone"), trkPt, track.dcaXY(), distRC, cone.referenceJetPtCorr, centrality);
          }
        }

        if (passesNominal) {
          if (hasTofPi) {
            registry.fill(HIST("tpcTofPiRandCone"), trkP, tpcPi, distRC, cone.referenceJetPtCorr, centrality);
            if (isTpcPiRange) {
              registry.fill(HIST("pVsPtForPiRandCone"), trkP, trkPt, distRC, cone.referenceJetPtCorr, centrality);
            }
          }
          if (hasTofPr) {
            registry.fill(HIST("tpcTofPrRandCone"), trkP, tpcPr, distRC, cone.referenceJetPtCorr, centrality);
            if (isTpcPrRange) {
              registry.fill(HIST("pVsPtForPrRandCone"), trkP, trkPt, distRC, cone.referenceJetPtCorr, centrality);
            }
          }
        }
      }
    }
  }

  void processJetProductionRatioWithRho(
    soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
    JetTracksWithPI const& tracks,
    FullTrackInfo const&,
    soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    runJetProductionRatio(collision, tracks, jets, collision.rho());
  }
  PROCESS_SWITCH(JetShapeTask, processJetProductionRatioWithRho, "production ratio around jets with rho subtraction", false);

  void processJetProductionRatioNoRho(
    soa::Filtered<aod::JetCollisions>::iterator const& collision,
    JetTracksWithPI const& tracks,
    FullTrackInfo const&,
    soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    runJetProductionRatio(collision, tracks, jets, 0.0f);
  }
  PROCESS_SWITCH(JetShapeTask, processJetProductionRatioNoRho, "production ratio around jets without rho table (recommended for pp)", false);

  void processInclusiveProductionRatio(soa::Filtered<aod::JetCollisions>::iterator const& collision, JetTracksWithPI const& tracks, FullTrackInfo const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    const float centrality = collision.centFT0M();
    if (!passCentrality(centrality)) {
      return;
    }
    registry.fill(HIST("eventCounterInc"), centrality);

    for (const auto& jetTrack : tracks) {
      const bool passesNominalBit = jetderiveddatautilities::selectTrack(jetTrack, trackSelection);
      const bool passesWideDcaBit = jetderiveddatautilities::selectTrack(jetTrack, trackSelectionForDca);
      if (!passesNominalBit && !passesWideDcaBit) {
        continue;
      }

      const auto track = jetTrack.track_as<FullTrackInfo>();

      // Keep the original QA tied to the nominal global-track selection.
      if (passesNominalBit && fillDetailedQA) {
        registry.fill(HIST("trackTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("trackDcaXY"), track.dcaXY());
        registry.fill(HIST("trackItsChi2NCl"), track.itsChi2NCl());
        registry.fill(HIST("trackTpcChi2NCl"), track.tpcChi2NCl());
        registry.fill(HIST("trackTpcNClsFound"), track.tpcNClsFound());
        registry.fill(HIST("trackItsNCls"), track.itsNCls());
        registry.fill(HIST("trackEta"), track.eta());
        registry.fill(HIST("trackPhi"), track.phi());
      }

      if (!passCommonTrackQuality(track)) {
        continue;
      }

      const bool passesNominal = passesNominalBit && passNominalTrackDca(track);
      const bool passesWideDca = passesWideDcaBit;

      const bool passPiKinematics = passSpeciesKinematics(track.pt(), track.eta(), PionMass);
      const bool passPrKinematics = passSpeciesKinematics(track.pt(), track.eta(), ProtonMass);
      const bool hasTofPi = passPiKinematics && track.hasTOF() && std::abs(track.tofNSigmaPi()) < nSigmaTofCut;
      const bool hasTofPr = passPrKinematics && track.hasTOF() && std::abs(track.tofNSigmaPr()) < nSigmaTofCut;
      const bool isTpcPiRange = track.tpcNSigmaPi() > tpcNSigmaPiMin && track.tpcNSigmaPi() < tpcNSigmaPiMax;
      const bool isTpcPrRange = track.tpcNSigmaPr() > tpcNSigmaPrMin && track.tpcNSigmaPr() < tpcNSigmaPrMax;

      if (passesWideDca) {
        if (hasTofPi && isTpcPiRange) {
          registry.fill(HIST("dcaPiInc"), track.pt(), track.dcaXY(), centrality);
        }
        if (hasTofPr && isTpcPrRange) {
          registry.fill(HIST("dcaPrInc"), track.pt(), track.dcaXY(), centrality);
        }
      }

      if (!passesNominal) {
        continue;
      }

      if (fillDetailedQA) {
        registry.fill(HIST("tpcPi"), track.p(), track.tpcNSigmaPi());
        registry.fill(HIST("tofPi"), track.pt(), track.tofNSigmaPi());
        registry.fill(HIST("tpcPr"), track.p(), track.tpcNSigmaPr());
        registry.fill(HIST("tofPr"), track.pt(), track.tofNSigmaPr());
        registry.fill(HIST("tpcDedx"), track.p(), track.tpcSignal(), centrality);
        registry.fill(HIST("tofBeta"), track.p(), track.beta(), centrality);
      }

      if (hasTofPr) {
        registry.fill(HIST("tpcTofPr"), track.p(), track.tpcNSigmaPr(), centrality);
        if (isTpcPrRange) {
          registry.fill(HIST("pVsPtForPr"), track.p(), track.pt(), centrality);
        }
      }

      if (hasTofPi) {
        registry.fill(HIST("tpcTofPi"), track.p(), track.tpcNSigmaPi(), centrality);
        if (isTpcPiRange) {
          registry.fill(HIST("pVsPtForPi"), track.p(), track.pt(), centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processInclusiveProductionRatio, "inclusive Production ratio", false);

  void processEfficiencyAndPurity(
    soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision,
    JetTracksWithPIMC const& jetTracks,
    FullTrackInfoMC const&,
    soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdJets,
    aod::ChargedMCParticleLevelJets const& mcpJets,
    aod::JetParticles const& mcParticles)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }

    const auto reducedMcCollisionId = collision.mcCollisionId();
    (void)mcpJets;

    const float centrality = collision.centFT0M();
    if (!passCentrality(centrality)) {
      return;
    }
    const float rho = collision.rho();
    const float maxR2 = distanceMax * distanceMax;

    registry.fill(HIST("eventCounterMc"), centrality);

    struct MatchedJet {
      float detPtCorr;
      float detEta;
      float detPhi;
      float detPhiBg1;
      float detPhiBg2;
      float partEta;
      float partPhi;
      float partPhiBg1;
      float partPhiBg2;
    };

    std::vector<MatchedJet> validJets;
    std::vector<ConeVetoJet> vetoJets;
    validJets.reserve(mcdJets.size());
    vetoJets.reserve(mcdJets.size());

    for (const auto& detJet : mcdJets) {
      if (!isAcceptedJet<JetTracksWithPIMC>(detJet)) {
        continue;
      }

      const auto matchedJets = detJet.matchedJetGeo();
      if (matchedJets.empty()) {
        continue;
      }

      const auto& partJet = matchedJets[0];
      if (partJet.mcCollisionId() != reducedMcCollisionId) {
        continue;
      }

      const float detPtCorr = detJet.pt() - rho * detJet.area();
      validJets.push_back({detPtCorr,
                           detJet.eta(),
                           detJet.phi(),
                           RecoDecay::constrainAngle(detJet.phi() + o2::constants::math::PIHalf),
                           RecoDecay::constrainAngle(detJet.phi() - o2::constants::math::PIHalf),
                           partJet.eta(),
                           partJet.phi(),
                           RecoDecay::constrainAngle(partJet.phi() + o2::constants::math::PIHalf),
                           RecoDecay::constrainAngle(partJet.phi() - o2::constants::math::PIHalf)});

      vetoJets.push_back({detJet.eta(), detJet.phi(), detPtCorr});
      registry.fill(HIST("jetPtMc"), detJet.pt(), centrality);
    }

    const auto randomCones = makeRandomCones(vetoJets);

    // Generated denominator for tracking efficiency
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.mcCollisionId() != reducedMcCollisionId) {
        continue;
      }
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (!passTruthKinematics(mcParticle)) {
        continue;
      }

      const int absPdg = std::abs(mcParticle.pdgCode());
      const bool isPi = absPdg == PDG_t::kPiPlus;
      const bool isKa = absPdg == PDG_t::kKPlus;
      const bool isPr = absPdg == PDG_t::kProton;
      if (!isPi && !isKa && !isPr) {
        continue;
      }

      const float mcPt = mcParticle.pt();

      if (isPi) {
        registry.fill(HIST("ptGenPiInc"), mcPt, centrality);
      } else if (isKa) {
        registry.fill(HIST("ptGenKaInc"), mcPt, centrality);
      } else {
        registry.fill(HIST("ptGenPrInc"), mcPt, centrality);
      }

      for (const auto& jet : validJets) {
        const float dEtaPart = mcParticle.eta() - jet.partEta;

        const float dPhiJetPart = absDeltaPhi(mcParticle.phi(), jet.partPhi);
        const float distJetSq = dEtaPart * dEtaPart + dPhiJetPart * dPhiJetPart;
        if (distJetSq < maxR2) {
          const float distJet = std::sqrt(distJetSq);
          if (isPi) {
            registry.fill(HIST("ptGenPi"), mcPt, distJet, jet.detPtCorr, centrality);
          } else if (isKa) {
            registry.fill(HIST("ptGenKa"), mcPt, distJet, jet.detPtCorr, centrality);
          } else {
            registry.fill(HIST("ptGenPr"), mcPt, distJet, jet.detPtCorr, centrality);
          }
        }

        const float dPhiBg1Part = absDeltaPhi(mcParticle.phi(), jet.partPhiBg1);
        const float dPhiBg2Part = absDeltaPhi(mcParticle.phi(), jet.partPhiBg2);
        const float distBg1Sq = dEtaPart * dEtaPart + dPhiBg1Part * dPhiBg1Part;
        const float distBg2Sq = dEtaPart * dEtaPart + dPhiBg2Part * dPhiBg2Part;

        float distBg = -1.0f;
        if (distBg1Sq < maxR2) {
          distBg = std::sqrt(distBg1Sq);
        } else if (distBg2Sq < maxR2) {
          distBg = std::sqrt(distBg2Sq);
        }

        if (distBg >= 0.0f) {
          if (isPi) {
            registry.fill(HIST("ptGenPiPerp"), mcPt, distBg, jet.detPtCorr, centrality);
          } else if (isKa) {
            registry.fill(HIST("ptGenKaPerp"), mcPt, distBg, jet.detPtCorr, centrality);
          } else {
            registry.fill(HIST("ptGenPrPerp"), mcPt, distBg, jet.detPtCorr, centrality);
          }
        }
      }

      for (const auto& cone : randomCones) {
        const float dEta = mcParticle.eta() - cone.eta;
        if (std::abs(dEta) >= distanceMax) {
          continue;
        }
        const float dPhi = absDeltaPhi(mcParticle.phi(), cone.phi);
        const float distSq = dEta * dEta + dPhi * dPhi;
        if (distSq >= maxR2) {
          continue;
        }
        const float distance = std::sqrt(distSq);
        if (isPi) {
          registry.fill(HIST("ptGenPiRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
        } else if (isKa) {
          registry.fill(HIST("ptGenKaRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
        } else {
          registry.fill(HIST("ptGenPrRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
        }
      }
    }

    // Reconstructed tracks
    for (const auto& jetTrack : jetTracks) {
      const bool passesNominalBit = jetderiveddatautilities::selectTrack(jetTrack, trackSelection);
      const bool passesWideDcaBit = jetderiveddatautilities::selectTrack(jetTrack, trackSelectionForDca);
      if (!passesNominalBit && !passesWideDcaBit) {
        continue;
      }

      const auto track = jetTrack.track_as<FullTrackInfoMC>();
      if (!passCommonTrackQuality(track)) {
        continue;
      }

      const bool passesNominal = passesNominalBit && passNominalTrackDca(track);
      const bool passesWideDca = passesWideDcaBit;
      if (!passesNominal && !passesWideDca) {
        continue;
      }

      if (!jetTrack.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = jetTrack.mcParticle_as<aod::JetParticles>();
      if (mcParticle.mcCollisionId() != reducedMcCollisionId) {
        continue;
      }

      const int absPdg = std::abs(mcParticle.pdgCode());
      const bool isTruePi = absPdg == PDG_t::kPiPlus;
      const bool isTrueKa = absPdg == PDG_t::kKPlus;
      const bool isTruePr = absPdg == PDG_t::kProton;
      if (!isTruePi && !isTrueKa && !isTruePr) {
        continue;
      }

      constexpr int ProducedByDecay = 4; // ROOT TMCProcess::kPDecay
      const bool isPrimary = mcParticle.isPhysicalPrimary();
      const bool isDecay = !isPrimary && mcParticle.getProcess() == ProducedByDecay;

      const float recoPt = track.pt();
      const float mcPt = mcParticle.pt();
      const float dcaXy = track.dcaXY();
      const bool hasTof = track.hasTOF();

      const bool passTpcPi = track.tpcNSigmaPi() > tpcNSigmaPiMin &&
                             track.tpcNSigmaPi() < tpcNSigmaPiMax;
      const bool passTpcPr = track.tpcNSigmaPr() > tpcNSigmaPrMin &&
                             track.tpcNSigmaPr() < tpcNSigmaPrMax;
      const bool isRecoPi = passSpeciesKinematics(recoPt, track.eta(), PionMass) &&
                            passTpcPi && hasTof &&
                            std::abs(track.tofNSigmaPi()) < nSigmaTofCut;
      const bool isRecoPr = passSpeciesKinematics(recoPt, track.eta(), ProtonMass) &&
                            passTpcPr && hasTof &&
                            std::abs(track.tofNSigmaPr()) < nSigmaTofCut;

      float originBin = 2.5f;
      if (isPrimary) {
        originBin = 0.5f;
      } else if (isDecay) {
        originBin = 1.5f;
      }

      // DCA templates and wide->global selection transfer
      auto fillInclusiveDcaAndTransfer = [&]() {
        if (isRecoPi && isTruePi) {
          if (passesWideDca) {
            if (isPrimary) {
              registry.fill(HIST("dcaPrimPiInc"), recoPt, dcaXy, centrality);
            } else if (isDecay) {
              registry.fill(HIST("dcaDecayPiInc"), recoPt, dcaXy, centrality);
            } else {
              registry.fill(HIST("dcaMatPiInc"), recoPt, dcaXy, centrality);
            }
            registry.fill(HIST("originCountWideInc"), recoPt, centrality, 0.5f, originBin);
          }
          if (passesNominal) {
            registry.fill(HIST("originCountGlobalInc"), recoPt, centrality, 0.5f, originBin);
          }
        }

        if (isRecoPr && isTruePr) {
          if (passesWideDca) {
            if (isPrimary) {
              registry.fill(HIST("dcaPrimPrInc"), recoPt, dcaXy, centrality);
            } else if (isDecay) {
              registry.fill(HIST("dcaDecayPrInc"), recoPt, dcaXy, centrality);
            } else {
              registry.fill(HIST("dcaMatPrInc"), recoPt, dcaXy, centrality);
            }
            registry.fill(HIST("originCountWideInc"), recoPt, centrality, 1.5f, originBin);
          }
          if (passesNominal) {
            registry.fill(HIST("originCountGlobalInc"), recoPt, centrality, 1.5f, originBin);
          }
        }
      };
      fillInclusiveDcaAndTransfer();

      for (const auto& jet : validJets) {
        const float dEtaDet = track.eta() - jet.detEta;

        const float dPhiJetDet = absDeltaPhi(track.phi(), jet.detPhi);
        const float distJetDetSq = dEtaDet * dEtaDet + dPhiJetDet * dPhiJetDet;
        const bool inJetDet = distJetDetSq < maxR2;
        const float distJetDet = inJetDet ? std::sqrt(distJetDetSq) : -1.0f;

        const float dPhiBg1Det = absDeltaPhi(track.phi(), jet.detPhiBg1);
        const float dPhiBg2Det = absDeltaPhi(track.phi(), jet.detPhiBg2);
        const float distBg1DetSq = dEtaDet * dEtaDet + dPhiBg1Det * dPhiBg1Det;
        const float distBg2DetSq = dEtaDet * dEtaDet + dPhiBg2Det * dPhiBg2Det;

        float distBgDet = -1.0f;
        if (distBg1DetSq < maxR2) {
          distBgDet = std::sqrt(distBg1DetSq);
        } else if (distBg2DetSq < maxR2) {
          distBgDet = std::sqrt(distBg2DetSq);
        }

        if (inJetDet) {
          if (isRecoPi && isTruePi) {
            if (passesWideDca) {
              if (isPrimary) {
                registry.fill(HIST("dcaPrimPi"), recoPt, dcaXy, distJetDet, jet.detPtCorr, centrality);
              } else if (isDecay) {
                registry.fill(HIST("dcaDecayPi"), recoPt, dcaXy, distJetDet, jet.detPtCorr, centrality);
              } else {
                registry.fill(HIST("dcaMatPi"), recoPt, dcaXy, distJetDet, jet.detPtCorr, centrality);
              }
              registry.fill(HIST("originCountWideJet"), recoPt, distJetDet, jet.detPtCorr, centrality, 0.5f, originBin);
            }
            if (passesNominal) {
              registry.fill(HIST("originCountGlobalJet"), recoPt, distJetDet, jet.detPtCorr, centrality, 0.5f, originBin);
            }
          }

          if (isRecoPr && isTruePr) {
            if (passesWideDca) {
              if (isPrimary) {
                registry.fill(HIST("dcaPrimPr"), recoPt, dcaXy, distJetDet, jet.detPtCorr, centrality);
              } else if (isDecay) {
                registry.fill(HIST("dcaDecayPr"), recoPt, dcaXy, distJetDet, jet.detPtCorr, centrality);
              } else {
                registry.fill(HIST("dcaMatPr"), recoPt, dcaXy, distJetDet, jet.detPtCorr, centrality);
              }
              registry.fill(HIST("originCountWideJet"), recoPt, distJetDet, jet.detPtCorr, centrality, 1.5f, originBin);
            }
            if (passesNominal) {
              registry.fill(HIST("originCountGlobalJet"), recoPt, distJetDet, jet.detPtCorr, centrality, 1.5f, originBin);
            }
          }
        }

        if (distBgDet >= 0.0f) {
          if (isRecoPi && isTruePi) {
            if (passesWideDca) {
              if (isPrimary) {
                registry.fill(HIST("dcaPrimPiPerp"), recoPt, dcaXy, distBgDet, jet.detPtCorr, centrality);
              } else if (isDecay) {
                registry.fill(HIST("dcaDecayPiPerp"), recoPt, dcaXy, distBgDet, jet.detPtCorr, centrality);
              } else {
                registry.fill(HIST("dcaMatPiPerp"), recoPt, dcaXy, distBgDet, jet.detPtCorr, centrality);
              }
              registry.fill(HIST("originCountWidePerp"), recoPt, distBgDet, jet.detPtCorr, centrality, 0.5f, originBin);
            }
            if (passesNominal) {
              registry.fill(HIST("originCountGlobalPerp"), recoPt, distBgDet, jet.detPtCorr, centrality, 0.5f, originBin);
            }
          }

          if (isRecoPr && isTruePr) {
            if (passesWideDca) {
              if (isPrimary) {
                registry.fill(HIST("dcaPrimPrPerp"), recoPt, dcaXy, distBgDet, jet.detPtCorr, centrality);
              } else if (isDecay) {
                registry.fill(HIST("dcaDecayPrPerp"), recoPt, dcaXy, distBgDet, jet.detPtCorr, centrality);
              } else {
                registry.fill(HIST("dcaMatPrPerp"), recoPt, dcaXy, distBgDet, jet.detPtCorr, centrality);
              }
              registry.fill(HIST("originCountWidePerp"), recoPt, distBgDet, jet.detPtCorr, centrality, 1.5f, originBin);
            }
            if (passesNominal) {
              registry.fill(HIST("originCountGlobalPerp"), recoPt, distBgDet, jet.detPtCorr, centrality, 1.5f, originBin);
            }
          }
        }
      }

      for (const auto& cone : randomCones) {
        const float dEta = track.eta() - cone.eta;
        if (std::abs(dEta) >= distanceMax) {
          continue;
        }
        const float dPhi = absDeltaPhi(track.phi(), cone.phi);
        const float distSq = dEta * dEta + dPhi * dPhi;
        if (distSq >= maxR2) {
          continue;
        }
        const float distance = std::sqrt(distSq);

        if (isRecoPi && isTruePi) {
          if (passesWideDca) {
            if (isPrimary) {
              registry.fill(HIST("dcaPrimPiRand"), recoPt, dcaXy, distance, cone.referenceJetPtCorr, centrality);
            } else if (isDecay) {
              registry.fill(HIST("dcaDecayPiRand"), recoPt, dcaXy, distance, cone.referenceJetPtCorr, centrality);
            } else {
              registry.fill(HIST("dcaMatPiRand"), recoPt, dcaXy, distance, cone.referenceJetPtCorr, centrality);
            }
            registry.fill(HIST("originCountWideRand"), recoPt, distance, cone.referenceJetPtCorr, centrality, 0.5f, originBin);
          }
          if (passesNominal) {
            registry.fill(HIST("originCountGlobalRand"), recoPt, distance, cone.referenceJetPtCorr, centrality, 0.5f, originBin);
          }
        }

        if (isRecoPr && isTruePr) {
          if (passesWideDca) {
            if (isPrimary) {
              registry.fill(HIST("dcaPrimPrRand"), recoPt, dcaXy, distance, cone.referenceJetPtCorr, centrality);
            } else if (isDecay) {
              registry.fill(HIST("dcaDecayPrRand"), recoPt, dcaXy, distance, cone.referenceJetPtCorr, centrality);
            } else {
              registry.fill(HIST("dcaMatPrRand"), recoPt, dcaXy, distance, cone.referenceJetPtCorr, centrality);
            }
            registry.fill(HIST("originCountWideRand"), recoPt, distance, cone.referenceJetPtCorr, centrality, 1.5f, originBin);
          }
          if (passesNominal) {
            registry.fill(HIST("originCountGlobalRand"), recoPt, distance, cone.referenceJetPtCorr, centrality, 1.5f, originBin);
          }
        }
      }

      // Tracking and TOF-matching efficiency numerators
      if (!passesNominal || !isPrimary) {
        continue;
      }
      if (!passTruthKinematics(mcParticle)) {
        continue;
      }

      registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcPt);

      if (isTruePi) {
        registry.fill(HIST("effNumTpcPiInc"), mcPt, centrality);
        if (hasTof) {
          registry.fill(HIST("effNumTofPiInc"), mcPt, centrality);
        }
        if (isRecoPi) {
          registry.fill(HIST("effNumPidPiInc"), mcPt, centrality);
        }
      } else if (isTrueKa) {
        registry.fill(HIST("effNumTpcKaInc"), mcPt, centrality);
        if (hasTof) {
          registry.fill(HIST("effNumTofKaInc"), mcPt, centrality);
        }
      } else {
        registry.fill(HIST("effNumTpcPrInc"), mcPt, centrality);
        if (hasTof) {
          registry.fill(HIST("effNumTofPrInc"), mcPt, centrality);
        }
        if (isRecoPr) {
          registry.fill(HIST("effNumPidPrInc"), mcPt, centrality);
        }
      }

      for (const auto& jet : validJets) {
        const float dEtaPart = mcParticle.eta() - jet.partEta;

        const float dPhiJetPart = absDeltaPhi(mcParticle.phi(), jet.partPhi);
        const float distJetPartSq = dEtaPart * dEtaPart + dPhiJetPart * dPhiJetPart;
        if (distJetPartSq < maxR2) {
          const float distJetPart = std::sqrt(distJetPartSq);
          if (isTruePi) {
            registry.fill(HIST("effNumTpcPiJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            if (hasTof) {
              registry.fill(HIST("effNumTofPiJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            }
            if (isRecoPi) {
              registry.fill(HIST("effNumPidPiJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            }
          } else if (isTrueKa) {
            registry.fill(HIST("effNumTpcKaJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            if (hasTof) {
              registry.fill(HIST("effNumTofKaJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            }
          } else {
            registry.fill(HIST("effNumTpcPrJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            if (hasTof) {
              registry.fill(HIST("effNumTofPrJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            }
            if (isRecoPr) {
              registry.fill(HIST("effNumPidPrJet"), mcPt, distJetPart, jet.detPtCorr, centrality);
            }
          }
        }

        const float dPhiBg1Part = absDeltaPhi(mcParticle.phi(), jet.partPhiBg1);
        const float dPhiBg2Part = absDeltaPhi(mcParticle.phi(), jet.partPhiBg2);
        const float distBg1PartSq = dEtaPart * dEtaPart + dPhiBg1Part * dPhiBg1Part;
        const float distBg2PartSq = dEtaPart * dEtaPart + dPhiBg2Part * dPhiBg2Part;

        float distBgPart = -1.0f;
        if (distBg1PartSq < maxR2) {
          distBgPart = std::sqrt(distBg1PartSq);
        } else if (distBg2PartSq < maxR2) {
          distBgPart = std::sqrt(distBg2PartSq);
        }

        if (distBgPart >= 0.0f) {
          if (isTruePi) {
            registry.fill(HIST("effNumTpcPiPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            if (hasTof) {
              registry.fill(HIST("effNumTofPiPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            }
            if (isRecoPi) {
              registry.fill(HIST("effNumPidPiPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            }
          } else if (isTrueKa) {
            registry.fill(HIST("effNumTpcKaPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            if (hasTof) {
              registry.fill(HIST("effNumTofKaPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            }
          } else {
            registry.fill(HIST("effNumTpcPrPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            if (hasTof) {
              registry.fill(HIST("effNumTofPrPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            }
            if (isRecoPr) {
              registry.fill(HIST("effNumPidPrPerp"), mcPt, distBgPart, jet.detPtCorr, centrality);
            }
          }
        }
      }

      for (const auto& cone : randomCones) {
        const float dEta = mcParticle.eta() - cone.eta;
        if (std::abs(dEta) >= distanceMax) {
          continue;
        }
        const float dPhi = absDeltaPhi(mcParticle.phi(), cone.phi);
        const float distSq = dEta * dEta + dPhi * dPhi;
        if (distSq >= maxR2) {
          continue;
        }
        const float distance = std::sqrt(distSq);

        if (isTruePi) {
          registry.fill(HIST("effNumTpcPiRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          if (hasTof) {
            registry.fill(HIST("effNumTofPiRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          }
          if (isRecoPi) {
            registry.fill(HIST("effNumPidPiRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          }
        } else if (isTrueKa) {
          registry.fill(HIST("effNumTpcKaRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          if (hasTof) {
            registry.fill(HIST("effNumTofKaRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          }
        } else {
          registry.fill(HIST("effNumTpcPrRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          if (hasTof) {
            registry.fill(HIST("effNumTofPrRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          }
          if (isRecoPr) {
            registry.fill(HIST("effNumPidPrRand"), mcPt, distance, cone.referenceJetPtCorr, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processEfficiencyAndPurity, "process MC information for Efficiency and Purity", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetShapeTask>(cfgc)}; }
