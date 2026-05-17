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

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetShapeTask {

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
  Configurable<float> tpcNSigmaPrMin{"tpcNSigmaPrMin", -3.5f, "Min value of tpcNsigmaProton"};
  Configurable<float> tpcNSigmaPrMax{"tpcNSigmaPrMax", 0.5f, "Max value of tpcNsigmaProton"};
  Configurable<float> tpcNSigmaPiMin{"tpcNSigmaPiMin", -0.5f, "Min value of tpcNsigmaPion"};
  Configurable<float> tpcNSigmaPiMax{"tpcNSigmaPiMax", 3.5f, "Max value of tpcNsigmaPion"};
  Configurable<float> randomConeDeltaPhiMin{"randomConeDeltaPhiMin", static_cast<float>(o2::constants::math::PIThird), "Minimum delta phi for random cone"};
  Configurable<float> randomConeDeltaPhiMax{"randomConeDeltaPhiMax", static_cast<float>(2.0f * o2::constants::math::PIThird), "Maximum delta phi for random cone"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8",
                                            "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", 5.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};

  // for jet shape
  Configurable<std::vector<float>> distanceCategory{"distanceCategory", {0.00f, 0.05f, 0.10f, 0.15f, 0.20f, 0.25f, 0.30f, 0.35f, 0.40f, 0.45f, 0.50f, 0.55f, 0.60f, 0.65f, 0.70f}, "distance of category"};

  // for ppi production
  Configurable<float> etaTrUp{"etaTrUp", 0.9f, "maximum track eta"};
  Configurable<float> dcaxyCutMax{"dcaxyCutMax", 2.0f, "maximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 4.0f, "tpc chi2 cut"};
  Configurable<float> nclItsMin{"nclItsMin", 2.0f, "its # of cluster cut"};
  Configurable<float> nclTpcMin{"nclTpcMin", 70.0f, "tpc # if cluster cut"};
  Configurable<float> nclcrossTpcMin{"nclcrossTpcMin", 70.0f, "tpc # of crossedRows cut"};
  Configurable<float> mcRapidityMax{"mcRapidityMax", 0.5f, "maximum mctrack y"};
  Configurable<double> epsilon{"epsilon", 1e-6, "standard for aboid division of zero"};
  Configurable<float> maxDeltaEtaSafe{"maxDeltaEtaSafe", 0.9f, "maximum track eta for cut"};
  Configurable<float> nSigmaMaxForDcaxy{"nSigmaMaxForDcaxy", 4.0f, "maximum nSigma for DCAxy"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};

  HistogramRegistry registry{"registry"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);

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
    registry.add("trackPhi", "trackPhi", HistType::kTH1F, {{80, -1, 7}});
    registry.add("trackEta", "trackEta", HistType::kTH1F, {{100, -1, 1}});
    registry.add("trackTpcNClsCrossedRows", "trackTpcNClsCrossedRows", HistType::kTH1F, {{50, 0, 200}});
    registry.add("trackDcaXY", "trackDcaXY", HistType::kTH1F, {{40, -10, 10}});
    registry.add("trackItsChi2NCl", "trackItsChi2NCl", HistType::kTH1F, {{60, 0, 30}});
    registry.add("trackTpcChi2NCl", "trackTpcChi2NCl", HistType::kTH1F, {{100, 0, 50}});
    registry.add("trackTpcNClsFound", "trackTpcNClsFound", HistType::kTH1F, {{100, 0, 200}});
    registry.add("trackItsNCls", "trackItsNCls", HistType::kTH1F, {{10, 0, 10}});

    registry.add("jetTpcTofPi", "jetTpcTofPi", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetTpcTofPr", "jetTpcTofPr", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcTofPiPerpJet", "tpcTofPiPerpJet", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcTofPrPerpJet", "tpcTofPrPerpJet", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("jetTpcPi", "jetTpcPi", HistType::kTH2F, {{nBinsP.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("jetTofPi", "jetTofPi", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("jetTpcPr", "jetTpcPr", HistType::kTH2F, {{nBinsP.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});
    registry.add("jetTofPr", "jetTofPr", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}});

    registry.add("jetTpcDedx", "jetTpcDedx", HistType::kTHnSparseD, {{nBinsPForDedx.value, 0, pMax.value}, {nBinsTpcDedx.value, 0, 1000}, {nBinsDistance.value, 0, distanceMax.value}});
    registry.add("tpcDedxPerpJet", "tpcDedxPerpJet", HistType::kTH2F, {{nBinsPForDedx.value, 0, pMax.value}, {nBinsTpcDedx.value, 0, 1000}});
    registry.add("jetTofBeta", "jetTofBeta", HistType::kTH2F, {{nBinsPForBeta.value, 0, pMax.value}, {nBinsTofBeta.value, 0.4, 1.1}});

    registry.add("jetpVsPtForPr", "jetpVsPtForPr", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetpVsPtForPi", "jetpVsPtPi", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPrPerpJet", "pVsPtForPrPerpJet", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPiPerpJet", "pVsPtPionPerpJet", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("jetDcaPr", "jetDcaPr", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetDcaPi", "jetDcaPi", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrPerpJet", "dcaPrPerpJet", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPiPerpJet", "dcaPiPerpJet", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("tpcTofPiRandCone", "tpcTofPiRandCone", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("tpcTofPrRandCone", "tpcTofPrRandCone", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPrRandCone", "pVsPtForPrRandCone", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("pVsPtForPiRandCone", "pVsPtPionRandCone", HistType::kTHnSparseD, {{nBinsP.value, 0, pMax.value}, {nBinsPt.value, 0, ptMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrRandCone", "dcaPrRandCone", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPiRandCone", "dcaPiRandCone", HistType::kTHnSparseD, {{nBinsPtForDca.value, 0, ptMax.value}, {nBinsDcaxyForData.value, dcaxyMin.value, dcaxyMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("jetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", HistType::kTH2F, {{200, 0., 200.}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetPtMc", "MC jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", HistType::kTH2F, {{200, 0., 200.}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("jetEta", "jet #eta;#eta_{jet};entries", HistType::kTH1F, {{100, -1.0, 1.0}});
    registry.add("jetPhi", "jet #phi;#phi_{jet};entries", HistType::kTH1F, {{80, -1.0, 7.}});
    registry.add("area", "area", HistType::kTH2F, {{100, 0, 2}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("rho", "rho", HistType::kTH2F, {{120, 0, 300}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptCorr", "Corrected jet pT; p_{T}^{corr} (GeV/c); Counts", HistType::kTH1F, {{200, 0, 200}});
    registry.add("ptCorrVsDistance", "ptcorr_vs_distance", HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}});
    registry.add("jetDistanceVsTrackpt", "trackpt_vs_distance_injet", HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}});
    registry.add("ptSum", "ptSum", HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc.value, 0, jetShapeFuncMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptSumBg1", "ptSumBg1", HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc.value, 0, jetShapeFuncMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("ptSumBg2", "ptSumBg2", HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc.value, 0, jetShapeFuncMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    registry.add("event/vertexz", ";Vtx_{z} (cm);Entries", HistType::kTH1F, {{100, -20, 20}});
    registry.add("eventCounterJetShape", "eventCounterJetShape", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterJet", "eventCounterJet", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterInc", "eventCounterInc", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("eventCounterRandCone", "Number of Random Cones;Centrality (%);Count", HistType::kTH1F, {{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
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

    // Jet (In-cone)
    registry.add("effNumTpcPiJet", "effNumTpcPiJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcKaJet", "effNumTpcKaJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcPrJet", "effNumTpcPrJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPiJet", "effNumTofPiJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofKaJet", "effNumTofKaJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPrJet", "effNumTofPrJet", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // Perp Cone
    registry.add("effNumTpcPiPerp", "effNumTpcPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcKaPerp", "effNumTpcKaPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcPrPerp", "effNumTpcPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPiPerp", "effNumTofPiPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofKaPerp", "effNumTofKaPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofPrPerp", "effNumTofPrPerp", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDistance.value, 0, distanceMax.value}, {nBinsJetPt.value, jetPtMinForCut.value, jetPtMaxForCut.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

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

    // Inclusive
    registry.add("dcaPrimPiInc", "dcaPrimPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPiInc", "dcaDecayPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPiInc", "dcaMatPiInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaPrimPrInc", "dcaPrimPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaDecayPrInc", "dcaDecayPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("dcaMatPrInc", "dcaMatPrInc", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsDcaxyForMc.value, dcaxyMin.value, dcaxyMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
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
    static constexpr double LeadingConstituentPtMinValue = 5.0;
    static constexpr double LeadingConstituentPtMaxValue = 9998.0;
    bool checkConstituentPt = true;
    bool checkConstituentMinPt =
      (leadingConstituentPtMin > LeadingConstituentPtMinValue);
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

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  using FullTrackInfo = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta>;
  using FullTrackInfoMC = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::pidTOFbeta>;

  void processJetShape(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, aod::JetTracks const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("eventCounterJetShape"), collision.centFT0M());

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
    float cent = collision.centFT0M();

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

        if (std::abs(dEta) > maxDistance)
          continue;

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

  void processJetProductionRatio(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, FullTrackInfo const&, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("eventCounterJet"), collision.centFT0M());
    registry.fill(HIST("event/vertexz"), collision.posZ());

    float rho = collision.rho();
    float centrality = collision.centFT0M();

    struct CachedJet {
      float pt, eta, phi, ptCorr;
      float phiBg1, phiBg2;
    };
    std::vector<CachedJet> cachedJets;
    cachedJets.reserve(jets.size());

    for (const auto& jet : jets) {
      if (!isAcceptedJet<aod::JetTracks>(jet))
        continue;

      float ptCorr = jet.pt() - rho * jet.area();
      float phiBg1 =
        RecoDecay::constrainAngle(jet.phi() + o2::constants::math::PIHalf);
      float phiBg2 =
        RecoDecay::constrainAngle(jet.phi() - o2::constants::math::PIHalf);

      registry.fill(HIST("jetEta"), jet.eta());
      registry.fill(HIST("jetPhi"), jet.phi());
      registry.fill(HIST("area"), jet.area(), centrality);
      registry.fill(HIST("rho"), rho, centrality);
      registry.fill(HIST("jetPt"), jet.pt(), centrality);

      cachedJets.push_back(
        {jet.pt(), jet.eta(), jet.phi(), ptCorr, phiBg1, phiBg2});
    }

    bool isValidRC = false;
    float rcEta = 0.0f;
    float rcPhi = 0.0f;

    if (!cachedJets.empty()) {
      constexpr unsigned int RandomSeed = 0;
      TRandom3 randomNumber(RandomSeed);

      const auto& leadJet = cachedJets[0];

      constexpr float MaxTrackEta = 0.9f;
      constexpr int MaxAttempts = 100; // for RandomCone
      constexpr float FlipProbability = 0.5f;

      // Range to generate
      float rcEtaMin = -MaxTrackEta + distanceMax;
      float rcEtaMax = MaxTrackEta - distanceMax;

      int attempts = 0;
      while (!isValidRC && attempts < MaxAttempts) {
        rcEta = randomNumber.Uniform(rcEtaMin, rcEtaMax);

        float dPhi = randomNumber.Uniform(randomConeDeltaPhiMin, randomConeDeltaPhiMax);

        // flipProbability (0.5)
        if (randomNumber.Uniform() < FlipProbability) {
          dPhi = -dPhi;
        }

        rcPhi = RecoDecay::constrainAngle(leadJet.phi + dPhi);

        float dPhiLead = std::abs(rcPhi - leadJet.phi);
        if (dPhiLead > o2::constants::math::PI)
          dPhiLead = o2::constants::math::TwoPI - dPhiLead;
        float dEtaLead = rcEta - leadJet.eta;
        float distLead = std::sqrt(dEtaLead * dEtaLead + dPhiLead * dPhiLead);

        if (distLead > (jetR + distanceMax)) {
          isValidRC = true;
        }
        attempts++;
      }

      if (isValidRC) {
        registry.fill(HIST("eventCounterRandCone"), centrality);
      }
    }

    for (const auto& jetTrack : tracks) {
      if (!jetderiveddatautilities::selectTrack(jetTrack, trackSelection)) {
        continue;
      }

      auto track = jetTrack.track_as<FullTrackInfo>();

      if (std::abs(track.eta()) > etaTrUp)
        continue;
      if (track.tpcNClsCrossedRows() < nclcrossTpcMin)
        continue;
      if (std::abs(track.dcaXY()) > dcaxyCutMax)
        continue;
      if (track.itsChi2NCl() > chi2ItsMax)
        continue;
      if (track.tpcChi2NCl() > chi2TpcMax)
        continue;
      if (track.tpcNClsFound() < nclTpcMin)
        continue;
      if (track.itsNCls() < nclItsMin)
        continue;

      float trkPt = track.pt();
      float trkP = track.p();
      float trkEta = track.eta();
      float trkPhi = track.phi();

      float tpcSig = track.tpcSignal();
      float tpcPi = track.tpcNSigmaPi();
      float tofPi = track.tofNSigmaPi();
      float tpcPr = track.tpcNSigmaPr();
      float tofPr = track.tofNSigmaPr();
      float beta = track.beta();

      // PID QA
      registry.fill(HIST("jetTpcPi"), trkP, tpcPi);
      registry.fill(HIST("jetTofPi"), trkPt, tofPi);
      registry.fill(HIST("jetTpcPr"), trkP, tpcPr);
      registry.fill(HIST("jetTofPr"), trkPt, tofPr);

      bool hasTofPi = (std::abs(tofPi) < nSigmaTofCut);
      bool hasTofPr = (std::abs(tofPr) < nSigmaTofCut);
      bool isTpcPiRange = (tpcPi > tpcNSigmaPiMin && tpcPi < tpcNSigmaPiMax);
      bool isTpcPrRange = (tpcPr > tpcNSigmaPrMin && tpcPr < tpcNSigmaPrMax);

      for (const auto& jet : cachedJets) {

        float dEta = trkEta - jet.eta;
        if (std::abs(dEta) > maxDeltaEtaSafe)
          continue;

        float dPhi = RecoDecay::constrainAngle(trkPhi - jet.phi);
        float distance = std::sqrt(dEta * dEta + dPhi * dPhi);

        // Background judge
        float deltaPhiBg1 = RecoDecay::constrainAngle(trkPhi - jet.phiBg1);
        float deltaPhiBg2 = RecoDecay::constrainAngle(trkPhi - jet.phiBg2);
        float distBg1 = std::sqrt(dEta * dEta + deltaPhiBg1 * deltaPhiBg1);
        float distBg2 = std::sqrt(dEta * dEta + deltaPhiBg2 * deltaPhiBg2);

        float distBg = -1.0f;
        if (distBg1 < distanceMax) {
          distBg = distBg1;
        } else if (distBg2 < distanceMax) {
          distBg = distBg2;
        }

        // --- Background Fill ---
        if (distBg >= 0.0f) {
          registry.fill(HIST("tpcDedxPerpJet"), trkP, tpcSig);

          // dcaXY
          if (track.hasTOF()) {
            if (hasTofPr && isTpcPrRange) {
              registry.fill(HIST("dcaPrPerpJet"), trkPt, track.dcaXY(), distBg, jet.ptCorr, centrality);
            }

            if (hasTofPi && isTpcPiRange) {
              registry.fill(HIST("dcaPiPerpJet"), trkPt, track.dcaXY(), distBg, jet.ptCorr, centrality);
            }
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

        if (isValidRC) {
          const auto& leadJet = cachedJets[0];

          float dEtaRC = trkEta - rcEta;
          float dPhiRC = std::abs(trkPhi - rcPhi);
          if (dPhiRC > o2::constants::math::PI)
            dPhiRC = o2::constants::math::TwoPI - dPhiRC;
          float distRC = std::sqrt(dEtaRC * dEtaRC + dPhiRC * dPhiRC);

          if (distRC < distanceMax) {

            // dcaXY
            if (track.hasTOF()) {
              if (hasTofPr && isTpcPrRange) {
                registry.fill(HIST("dcaPrRandCone"), trkPt, track.dcaXY(), leadJet.ptCorr, centrality);
              }
              if (hasTofPi && isTpcPiRange) {
                registry.fill(HIST("dcaPiRandCone"), trkPt, track.dcaXY(), leadJet.ptCorr, centrality);
              }
            }

            if (hasTofPi) {
              registry.fill(HIST("tpcTofPiRandCone"), trkP, tpcPi, leadJet.ptCorr, centrality);
              if (isTpcPiRange) {
                registry.fill(HIST("pVsPtForPiRandCone"), trkP, trkPt, leadJet.ptCorr, centrality);
              }
            }
            if (hasTofPr) {
              registry.fill(HIST("tpcTofPrRandCone"), trkP, tpcPr, leadJet.ptCorr, centrality);
              if (isTpcPrRange) {
                registry.fill(HIST("pVsPtForPrRandCone"), trkP, trkPt, leadJet.ptCorr, centrality);
              }
            }
          }
        }

        // --- Signal Fill ---
        registry.fill(HIST("jetDistanceVsTrackpt"), distance, trkPt);
        registry.fill(HIST("jetTpcDedx"), trkP, tpcSig, distance);
        registry.fill(HIST("jetTofBeta"), trkP, beta);

        // dcaXY
        if (track.hasTOF()) {
          if (hasTofPr && isTpcPrRange) {
            registry.fill(HIST("jetDcaPr"), trkPt, track.dcaXY(), distance, jet.ptCorr, centrality);
          }

          if (hasTofPi && isTpcPiRange) {
            registry.fill(HIST("jetDcaPi"), trkPt, track.dcaXY(), distance, jet.ptCorr, centrality);
          }
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
  }
  PROCESS_SWITCH(JetShapeTask, processJetProductionRatio, "production ratio around jets", false);

  void processInclusiveProductionRatio(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, FullTrackInfo const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("eventCounterInc"), collision.centFT0M());
    // tracks conditions
    for (const auto& jetTrack : tracks) {

      if (!jetderiveddatautilities::selectTrack(jetTrack, trackSelection)) {
        continue;
      }

      auto track = jetTrack.track_as<FullTrackInfo>();

      registry.fill(HIST("trackTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("trackDcaXY"), track.dcaXY());
      registry.fill(HIST("trackItsChi2NCl"), track.itsChi2NCl());
      registry.fill(HIST("trackTpcChi2NCl"), track.tpcChi2NCl());
      registry.fill(HIST("trackTpcNClsFound"), track.tpcNClsFound());
      registry.fill(HIST("trackItsNCls"), track.itsNCls());
      registry.fill(HIST("trackEta"), track.eta());
      registry.fill(HIST("trackPhi"), track.phi());

      if (std::abs(track.eta()) > etaTrUp)
        continue;
      if (track.tpcNClsCrossedRows() < nclcrossTpcMin)
        continue;
      if (std::abs(track.dcaXY()) > dcaxyCutMax)
        continue;
      if (track.itsChi2NCl() > chi2ItsMax)
        continue;
      if (track.tpcChi2NCl() > chi2TpcMax)
        continue;
      if (track.tpcNClsFound() < nclTpcMin)
        continue;
      if (track.itsNCls() < nclItsMin)
        continue;

      // PID check
      registry.fill(HIST("tpcPi"), track.p(), track.tpcNSigmaPi());
      registry.fill(HIST("tofPi"), track.pt(), track.tofNSigmaPi());
      registry.fill(HIST("tpcPr"), track.p(), track.tpcNSigmaPr());
      registry.fill(HIST("tofPr"), track.pt(), track.tofNSigmaPr());
      registry.fill(HIST("tpcDedx"), track.p(), track.tpcSignal(), collision.centFT0M());
      registry.fill(HIST("tofBeta"), track.p(), track.beta(), collision.centFT0M());

      if (std::abs(track.tofNSigmaPr()) < nSigmaTofCut) {
        registry.fill(HIST("tpcTofPr"), track.p(), track.tpcNSigmaPr(), collision.centFT0M());

        if (track.tpcNSigmaPr() > tpcNSigmaPrMin && track.tpcNSigmaPr() < tpcNSigmaPrMax) {
          registry.fill(HIST("pVsPtForPr"), track.p(), track.pt(), collision.centFT0M());
        }
      }

      if (std::abs(track.tofNSigmaPi()) < nSigmaTofCut) {
        registry.fill(HIST("tpcTofPi"), track.p(), track.tpcNSigmaPi(), collision.centFT0M());
        if (track.tpcNSigmaPi() > tpcNSigmaPiMin && track.tpcNSigmaPi() < tpcNSigmaPiMax) {
          registry.fill(HIST("pVsPtForPi"), track.p(), track.pt(), collision.centFT0M());
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processInclusiveProductionRatio, "inclusive Production ratio", false);

  void processEfficiencyAndPurity(
    soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision,
    soa::Join<aod::JetTracks, aod::JTrackPIs> const& jetTracks,
    FullTrackInfoMC const&,
    soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdJets,
    aod::ChargedMCParticleLevelJets const& mcpJets,
    aod::McParticles const& mcParticles)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;

    (void)mcpJets;

    float centrality = collision.centFT0M();
    float rho = collision.rho();
    float distMax = distanceMax;
    float maxR2 = distMax * distMax;

    registry.fill(HIST("eventCounterMc"), centrality);

    // jet matching
    struct MatchedJet {
      float detPtCorr, partEta, partPhi, partPhiBg1, partPhiBg2;
    };
    std::vector<MatchedJet> validJets;
    validJets.reserve(mcdJets.size());

    for (auto const& detJet : mcdJets) {
      auto matchedIndices = detJet.matchedJetGeo();

      if (matchedIndices.size() > 0) {
        auto const& partJet = matchedIndices[0];

        float detPtCorr = detJet.pt() - rho * detJet.area();
        float partPhi = partJet.phi();
        float phiBg1 = RecoDecay::constrainAngle(partPhi + o2::constants::math::PIHalf);
        float phiBg2 = RecoDecay::constrainAngle(partPhi - o2::constants::math::PIHalf);

        validJets.push_back({detPtCorr, partJet.eta(), partPhi, phiBg1, phiBg2});
        registry.fill(HIST("jetPtMc"), detJet.pt(), centrality);
      }
    }

    // Denominator: True Primary Particles
    for (const auto& mcParticle : mcParticles) {

      if (mcParticle.mcCollisionId() != collision.globalIndex())
        continue;

      if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.y()) > mcRapidityMax)
        continue;

      int absPdg = std::abs(mcParticle.pdgCode());
      bool isPi = (absPdg == PDG_t::kPiPlus);
      bool isKa = (absPdg == PDG_t::kKPlus);
      bool isPr = (absPdg == PDG_t::kProton);
      if (!isPi && !isKa && !isPr)
        continue;

      float mcPt = mcParticle.pt();

      if (isPi)
        registry.fill(HIST("ptGenPiInc"), mcPt, centrality);
      else if (isKa)
        registry.fill(HIST("ptGenKaInc"), mcPt, centrality);
      else if (isPr)
        registry.fill(HIST("ptGenPrInc"), mcPt, centrality);

      for (const auto& jet : validJets) {
        float dEta = mcParticle.eta() - jet.partEta;

        // Jet In-cone
        float dPhiJet = std::abs(mcParticle.phi() - jet.partPhi);
        if (dPhiJet > o2::constants::math::PI)
          dPhiJet = o2::constants::math::TwoPI - dPhiJet;

        float distJetGen = std::sqrt(dEta * dEta + dPhiJet * dPhiJet);

        if ((dEta * dEta + dPhiJet * dPhiJet) < maxR2) {
          if (isPi)
            registry.fill(HIST("ptGenPi"), mcPt, distJetGen, jet.detPtCorr, centrality);
          else if (isKa)
            registry.fill(HIST("ptGenKa"), mcPt, distJetGen, jet.detPtCorr, centrality);
          else if (isPr)
            registry.fill(HIST("ptGenPr"), mcPt, distJetGen, jet.detPtCorr, centrality);
        }

        // Perp Cone
        float dPhiBg1 = std::abs(mcParticle.phi() - jet.partPhiBg1);
        if (dPhiBg1 > o2::constants::math::PI)
          dPhiBg1 = o2::constants::math::TwoPI - dPhiBg1;
        float dPhiBg2 = std::abs(mcParticle.phi() - jet.partPhiBg2);
        if (dPhiBg2 > o2::constants::math::PI)
          dPhiBg2 = o2::constants::math::TwoPI - dPhiBg2;

        float distBg1Sq = dEta * dEta + dPhiBg1 * dPhiBg1;
        float distBg2Sq = dEta * dEta + dPhiBg2 * dPhiBg2;

        float distBgGen = -1.0f;
        if (distBg1Sq < maxR2)
          distBgGen = std::sqrt(distBg1Sq);
        else if (distBg2Sq < maxR2)
          distBgGen = std::sqrt(distBg2Sq);

        if (distBgGen >= 0.0f) {
          if (isPi)
            registry.fill(HIST("ptGenPiPerp"), mcPt, distBgGen, jet.detPtCorr, centrality);
          else if (isKa)
            registry.fill(HIST("ptGenKaPerp"), mcPt, distBgGen, jet.detPtCorr, centrality);
          else if (isPr)
            registry.fill(HIST("ptGenPrPerp"), mcPt, distBgGen, jet.detPtCorr, centrality);
        }
      }
    }

    // Numerator
    for (const auto& jetTrack : jetTracks) {
      if (!jetderiveddatautilities::selectTrack(jetTrack, trackSelection))
        continue;

      auto track = jetTrack.track_as<FullTrackInfoMC>();

      if (std::abs(track.eta()) > etaTrUp)
        continue;
      if (track.tpcNClsCrossedRows() < nclcrossTpcMin)
        continue;
      if (std::abs(track.dcaXY()) > dcaxyCutMax)
        continue;
      if (track.itsChi2NCl() > chi2ItsMax)
        continue;
      if (track.tpcChi2NCl() > chi2TpcMax)
        continue;
      if (track.tpcNClsFound() < nclTpcMin)
        continue;
      if (track.itsNCls() < nclItsMin)
        continue;

      if (!track.has_mcParticle())
        continue;
      auto mcParticle = track.mcParticle();

      if (mcParticle.mcCollisionId() != collision.globalIndex())
        continue;

      registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());

      if (std::abs(mcParticle.y()) >= mcRapidityMax)
        continue;

      int absPdg = std::abs(mcParticle.pdgCode());
      bool isTruePi = (absPdg == PDG_t::kPiPlus);
      bool isTrueKa = (absPdg == PDG_t::kKPlus);
      bool isTruePr = (absPdg == PDG_t::kProton);
      if (!isTruePi && !isTrueKa && !isTruePr)
        continue;

      const int producedByDecay = 4;
      bool isPrimary = mcParticle.isPhysicalPrimary();
      bool isDecay = !isPrimary && (mcParticle.getProcess() == producedByDecay);
      bool isMaterial = !isPrimary && !isDecay;

      float mcPt = mcParticle.pt();
      float recoPt = track.pt();
      float dcaXy = track.dcaXY();
      bool hasTof = track.hasTOF();

      // PID flag
      // TPC nsigma
      bool passTpcPiOnly = (track.tpcNSigmaPi() > tpcNSigmaPiMin && track.tpcNSigmaPi() < tpcNSigmaPiMax);
      bool passTpcPrOnly = (track.tpcNSigmaPr() > tpcNSigmaPrMin && track.tpcNSigmaPr() < tpcNSigmaPrMax);

      // TPC + TOF
      bool isRecoPi = passTpcPiOnly && hasTof && (std::abs(track.tofNSigmaPi()) < nSigmaTofCut);
      bool isRecoPr = passTpcPrOnly && hasTof && (std::abs(track.tofNSigmaPr()) < nSigmaTofCut);

      if (isPrimary) {
        // Inclusive
        // TPC Efficiency
        if (isTruePi)
          registry.fill(HIST("effNumTpcPiInc"), mcPt, centrality);
        if (isTrueKa)
          registry.fill(HIST("effNumTpcKaInc"), mcPt, centrality);
        if (isTruePr)
          registry.fill(HIST("effNumTpcPrInc"), mcPt, centrality);

        // TOF Matching Efficiency
        if (isTruePi && hasTof)
          registry.fill(HIST("effNumTofPiInc"), mcPt, centrality);
        if (isTrueKa && hasTof)
          registry.fill(HIST("effNumTofKaInc"), mcPt, centrality);
        if (isTruePr && hasTof)
          registry.fill(HIST("effNumTofPrInc"), mcPt, centrality);

        // Purity DCA
        if (isRecoPi && isTruePi)
          registry.fill(HIST("dcaPrimPiInc"), recoPt, dcaXy, centrality);
        if (isRecoPr && isTruePr)
          registry.fill(HIST("dcaPrimPrInc"), recoPt, dcaXy, centrality);
      } else if (isDecay) {
        if (isRecoPi && isTruePi)
          registry.fill(HIST("dcaDecayPiInc"), recoPt, dcaXy, centrality);
        if (isRecoPr && isTruePr)
          registry.fill(HIST("dcaDecayPrInc"), recoPt, dcaXy, centrality);
      } else if (isMaterial) {
        if (isRecoPi && isTruePi)
          registry.fill(HIST("dcaMatPiInc"), recoPt, dcaXy, centrality);
        if (isRecoPr && isTruePr)
          registry.fill(HIST("dcaMatPrInc"), recoPt, dcaXy, centrality);
      }

      for (const auto& jet : validJets) {
        float dEta = mcParticle.eta() - jet.partEta;

        // --- Jet In-cone distance check ---
        float dPhiJet = std::abs(mcParticle.phi() - jet.partPhi);
        if (dPhiJet > o2::constants::math::PI)
          dPhiJet = o2::constants::math::TwoPI - dPhiJet;

        float distJetReco = std::sqrt(dEta * dEta + dPhiJet * dPhiJet);
        bool inJet = (distJetReco < distanceMax);

        // --- Perp Cone distance check ---
        float dPhiBg1 = std::abs(mcParticle.phi() - jet.partPhiBg1);
        if (dPhiBg1 > o2::constants::math::PI)
          dPhiBg1 = o2::constants::math::TwoPI - dPhiBg1;
        float dPhiBg2 = std::abs(mcParticle.phi() - jet.partPhiBg2);
        if (dPhiBg2 > o2::constants::math::PI)
          dPhiBg2 = o2::constants::math::TwoPI - dPhiBg2;

        float distBg1Sq = dEta * dEta + dPhiBg1 * dPhiBg1;
        float distBg2Sq = dEta * dEta + dPhiBg2 * dPhiBg2;

        float distBgReco = -1.0f;
        if (distBg1Sq < maxR2)
          distBgReco = std::sqrt(distBg1Sq);
        else if (distBg2Sq < maxR2)
          distBgReco = std::sqrt(distBg2Sq);

        bool inPerp = (distBgReco >= 0.0f);

        if (inJet) {
          if (isPrimary) {
            // --- Jet In-cone 分子 ---
            // TPC Efficiency
            if (isTruePi)
              registry.fill(HIST("effNumTpcPiJet"), mcPt, distJetReco, jet.detPtCorr, centrality);
            if (isTrueKa)
              registry.fill(HIST("effNumTpcKaJet"), mcPt, distJetReco, jet.detPtCorr, centrality);
            if (isTruePr)
              registry.fill(HIST("effNumTpcPrJet"), mcPt, distJetReco, jet.detPtCorr, centrality);

            // TOF Matching Efficiency
            if (isTruePi && hasTof)
              registry.fill(HIST("effNumTofPiJet"), mcPt, distJetReco, jet.detPtCorr, centrality);
            if (isTrueKa && hasTof)
              registry.fill(HIST("effNumTofKaJet"), mcPt, distJetReco, jet.detPtCorr, centrality);
            if (isTruePr && hasTof)
              registry.fill(HIST("effNumTofPrJet"), mcPt, distJetReco, jet.detPtCorr, centrality);

            // Purity DCA
            if (isRecoPi && isTruePi)
              registry.fill(HIST("dcaPrimPi"), recoPt, dcaXy, distJetReco, jet.detPtCorr, centrality);
            if (isRecoPr && isTruePr)
              registry.fill(HIST("dcaPrimPr"), recoPt, dcaXy, distJetReco, jet.detPtCorr, centrality);
          } else if (isDecay) {
            if (isRecoPi && isTruePi)
              registry.fill(HIST("dcaDecayPi"), recoPt, dcaXy, distJetReco, jet.detPtCorr, centrality);
            if (isRecoPr && isTruePr)
              registry.fill(HIST("dcaDecayPr"), recoPt, dcaXy, distJetReco, jet.detPtCorr, centrality);
          } else if (isMaterial) {
            if (isRecoPi && isTruePi)
              registry.fill(HIST("dcaMatPi"), recoPt, dcaXy, distJetReco, jet.detPtCorr, centrality);
            if (isRecoPr && isTruePr)
              registry.fill(HIST("dcaMatPr"), recoPt, dcaXy, distJetReco, jet.detPtCorr, centrality);
          }
        }

        if (inPerp) {
          if (isPrimary) {
            // Perp Cone
            // TPC Efficiency
            if (isTruePi)
              registry.fill(HIST("effNumTpcPiPerp"), mcPt, distBgReco, jet.detPtCorr, centrality);
            if (isTrueKa)
              registry.fill(HIST("effNumTpcKaPerp"), mcPt, distBgReco, jet.detPtCorr, centrality);
            if (isTruePr)
              registry.fill(HIST("effNumTpcPrPerp"), mcPt, distBgReco, jet.detPtCorr, centrality);

            // TOF Matching Efficiency
            if (isTruePi && hasTof)
              registry.fill(HIST("effNumTofPiPerp"), mcPt, distBgReco, jet.detPtCorr, centrality);
            if (isTrueKa && hasTof)
              registry.fill(HIST("effNumTofKaPerp"), mcPt, distBgReco, jet.detPtCorr, centrality);
            if (isTruePr && hasTof)
              registry.fill(HIST("effNumTofPrPerp"), mcPt, distBgReco, jet.detPtCorr, centrality);

            // Purity DCA
            if (isRecoPi && isTruePi)
              registry.fill(HIST("dcaPrimPiPerp"), recoPt, dcaXy, distBgReco, jet.detPtCorr, centrality);
            if (isRecoPr && isTruePr)
              registry.fill(HIST("dcaPrimPrPerp"), recoPt, dcaXy, distBgReco, jet.detPtCorr, centrality);
          } else if (isDecay) {
            if (isRecoPi && isTruePi)
              registry.fill(HIST("dcaDecayPiPerp"), recoPt, dcaXy, distBgReco, jet.detPtCorr, centrality);
            if (isRecoPr && isTruePr)
              registry.fill(HIST("dcaDecayPrPerp"), recoPt, dcaXy, distBgReco, jet.detPtCorr, centrality);
          } else if (isMaterial) {
            if (isRecoPi && isTruePi)
              registry.fill(HIST("dcaMatPiPerp"), recoPt, dcaXy, distBgReco, jet.detPtCorr, centrality);
            if (isRecoPr && isTruePr)
              registry.fill(HIST("dcaMatPrPerp"), recoPt, dcaXy, distBgReco, jet.detPtCorr, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processEfficiencyAndPurity, "process MC information for Efficiency and Purity", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetShapeTask>(cfgc)}; }
