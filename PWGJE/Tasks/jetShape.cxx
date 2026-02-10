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

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <cmath>
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

  HistogramRegistry registry{
    "registry",
    {{"tpcTofPi", "tpcTofPi", {HistType::kTHnSparseD, {{nBinsPForCut, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"tpcTofPr", "tpcTofPr", {HistType::kTHnSparseD, {{nBinsPForCut, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"tpcPi", "tpcPi", {HistType::kTH2F, {{nBinsP, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"tofPi", "tofPi", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"tpcPr", "tpcPr", {HistType::kTH2F, {{nBinsP, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"tofPr", "tofPr", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"tpcDedx", "tpcDedx", {HistType::kTHnSparseD, {{nBinsPForDedx, 0, pMax}, {nBinsTpcDedx, 0, 1000}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"tofBeta", "tofBeta", {HistType::kTHnSparseD, {{nBinsPForBeta, 0, pMax}, {nBinsTofBeta, 0.4, 1.1}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"pVsPtForPr", "pVsPtForPr", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"pVsPtForPi", "pVsPtPi", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"trackPhi", "trackPhi", {HistType::kTH1F, {{80, -1, 7}}}},
     {"trackEta", "trackEta", {HistType::kTH1F, {{100, -1, 1}}}},
     {"trackTpcNClsCrossedRows", "trackTpcNClsCrossedRows", {HistType::kTH1F, {{50, 0, 200}}}},
     {"trackDcaXY", "trackDcaXY", {HistType::kTH1F, {{40, -10, 10}}}},
     {"trackItsChi2NCl", "trackItsChi2NCl", {HistType::kTH1F, {{60, 0, 30}}}},
     {"trackTpcChi2NCl", "trackTpcChi2NCl", {HistType::kTH1F, {{100, 0, 50}}}},
     {"trackTpcNClsFound", "trackTpcNClsFound", {HistType::kTH1F, {{100, 0, 200}}}},
     {"trackItsNCls", "trackItsNCls", {HistType::kTH1F, {{10, 0, 10}}}},
     {"jetTpcTofPi", "jetTpcTofPi", {HistType::kTHnSparseD, {{nBinsPForCut, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"jetTpcTofPr", "jetTpcTofPr", {HistType::kTHnSparseD, {{nBinsPForCut, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"tpcTofPiOutOfJet", "tpcTofPiOutOfJet", {HistType::kTHnSparseD, {{nBinsPForCut, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"tpcTofPrOutOfJet", "tpcTofPrOutOfJet", {HistType::kTHnSparseD, {{nBinsPForCut, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"jetTpcPi", "jetTpcPi", {HistType::kTH2F, {{nBinsP, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"jetTofPi", "jetTofPi", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"jetTpcPr", "jetTpcPr", {HistType::kTH2F, {{nBinsP, 0, pMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"jetTofPr", "jetTofPr", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {nBinsNSigma, nSigmaMin, nSigmaMax}}}},
     {"jetTpcDedx", "jetTpcDedx", {HistType::kTHnSparseD, {{nBinsPForDedx, 0, pMax}, {nBinsTpcDedx, 0, 1000}, {nBinsDistance, 0, distanceMax}}}},
     {"tpcDedxOutOfJet", "tpcDedxOutOfJet", {HistType::kTH2F, {{nBinsPForDedx, 0, pMax}, {nBinsTpcDedx, 0, 1000}}}},
     {"jetTofBeta", "jetTofBeta", {HistType::kTH2F, {{nBinsPForBeta, 0, pMax}, {nBinsTofBeta, 0.4, 1.1}}}},
     {"jetpVsPtForPr", "jetpVsPtForPr", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"jetpVsPtForPi", "jetpVsPtPi", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"pVsPtForPrOutOfJet", "pVsPtForPrOutOfJet", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"pVsPtForPiOutOfJet", "pVsPtPionOutOfJet", {HistType::kTHnSparseD, {{nBinsP, 0, pMax}, {nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"jetDcaPr", "jetDcaPr", {HistType::kTHnSparseD, {{nBinsPtForDca, 0, ptMax}, {nBinsDcaxyForData, dcaxyMin, dcaxyMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"jetDcaPi", "jetDcaPi", {HistType::kTHnSparseD, {{nBinsPtForDca, 0, ptMax}, {nBinsDcaxyForData, dcaxyMin, dcaxyMax}, {nBinsDistance, 0, distanceMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"dcaPrOutOfJet", "dcaPrOutOfJet", {HistType::kTHnSparseD, {{nBinsPtForDca, 0, ptMax}, {nBinsDcaxyForData, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"dcaPiOutOfJet", "dcaPiOutOfJet", {HistType::kTHnSparseD, {{nBinsPtForDca, 0, ptMax}, {nBinsDcaxyForData, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"jetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
     {"jetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
     {"jetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
     {"area", "area", {HistType::kTH1F, {{100, 0, 4}}}},
     {"rho", "rho", {HistType::kTH1F, {{120, 0, 300}}}},
     {"ptCorr", "Corrected jet pT; p_{T}^{corr} (GeV/c); Counts", {HistType::kTH1F, {{200, 0, 200}}}},
     {"ptCorrVsDistance", "ptcorr_vs_distance", {HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}}}},
     {"jetDistanceVsTrackpt", "trackpt_vs_distance_injet", {HistType::kTH2F, {{70, 0, 0.7}, {100, 0, 100}}}},
     {"ptSum", "ptSum", {HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc, 0, jetShapeFuncMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptSumBg1", "ptSumBg1", {HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc, 0, jetShapeFuncMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptSumBg2", "ptSumBg2", {HistType::kTHnSparseD, {{14, 0, 0.7}, {nBinsJetShapeFunc, 0, jetShapeFuncMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"event/vertexz", ";Vtx_{z} (cm);Entries", {HistType::kTH1F, {{100, -20, 20}}}},
     {"eventCounterJetShape", "eventCounterJetShape", {HistType::kTH1F, {{1, 0, +1, ""}}}},
     {"eventCounterJet", "eventCounterJet", {HistType::kTH1F, {{1, 0, +1, ""}}}},
     {"eventCounterInc", "eventCounterInc", {HistType::kTH1F, {{1, 0, +1, ""}}}},
     {"eventCounterMc", "eventCounterMc", {HistType::kTH1F, {{1, 0, +1, ""}}}},
     {"ptVsCentrality", "ptvscentrality", {HistType::kTH2F, {{100, 0, 100}, {300, 0, 300}}}},
     {"ptResolution", "ptResolution", {HistType::kTH2F, {{nBinsPt, 0, ptMax}, {100, -1.0, +1.0}}}},
     {"mcCentralityReco", "mcCentralityReco", {HistType::kTH1F, {{100, 0, 100}}}},
     {"mcCentralitySim", "mcCentralitySim", {HistType::kTH1F, {{100, 0, 100}}}},
     {"ptHistogramPion", "ptHistogramPion", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptHistogramKaon", "ptHistogramKaon", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptHistogramProton", "ptHistogramProton", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptHistogramPionTof", "ptHistogramPionTof", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptHistogramKaonTof", "ptHistogramKaonTof", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptHistogramProtonTof", "ptHistogramProtonTof", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"dcaDecayPion", "dcaDecayPion", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"dcaDecayProton", "dcaDecayProton", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"dcaMaterialPion", "dcaMaterialPion", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"dcaMaterialProton", "dcaMaterialProton", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsDcaxyForMc, dcaxyMin, dcaxyMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptGeneratedPion", "ptGeneratedPion", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptGeneratedKaon", "ptGeneratedKaon", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}},
     {"ptGeneratedProton", "ptGeneratedProton", {HistType::kTHnSparseD, {{nBinsPt, 0, ptMax}, {nBinsJetPt, jetPtMinForCut, jetPtMaxForCut}, {nBinsCentrality, centralityMinForCut, centralityMaxForCut}}}}}};

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
  Configurable<float> etaTrUp{"etaTrUp", 0.7f, "maximum track eta"};
  Configurable<float> dcaxyCutMax{"dcaxyCutMax", 2.0f, "maximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 4.0f, "tpc chi2 cut"};
  Configurable<float> nclItsMin{"nclItsMin", 2.0f, "its # of cluster cut"};
  Configurable<float> nclTpcMin{"nclTpcMin", 100.0f, "tpc # if cluster cut"};
  Configurable<float> nclcrossTpcMin{"nclcrossTpcMin", 70.0f, "tpc # of crossedRows cut"};
  Configurable<float> mcRapidityMax{"mcRapidityMax", 0.5f, "maximum mctrack y"};
  Configurable<double> epsilon{"epsilon", 1e-6, "standard for aboid division of zero"};
  Configurable<float> maxDeltaEtaSafe{"maxDeltaEtaSafe", 0.9f, "maximum track eta for cut"};
  Configurable<float> nSigmaMaxForDcaxy{"nSigmaMaxForDcaxy", 4.0f, "maximum nSigma for DCAxy"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
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

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r ==
                   nround(jetR.node() * 100.0f);
  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  using FullTrackInfo = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta>;

  void processJetShape(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, aod::JetTracks const& tracks, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("eventCounterJetShape"), 0.5);

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

      registry.fill(HIST("area"), jet.area);
      registry.fill(HIST("rho"), rho);
      registry.fill(HIST("jetPt"), jet.pt);
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
    registry.fill(HIST("eventCounterJet"), 0.5);
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

      cachedJets.push_back(
        {jet.pt(), jet.eta(), jet.phi(), ptCorr, phiBg1, phiBg2});
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

      float nSigmaSqPr = tpcPr * tpcPr + tofPr * tofPr;
      float nSigmaSqPi = tpcPi * tpcPi + tofPi * tofPi;

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

        // --- Background Fill ---
        if (distBg1 < distanceMax || distBg2 < distanceMax) {
          registry.fill(HIST("tpcDedxOutOfJet"), trkP, tpcSig);

          // dcaXY
          if (track.hasTOF()) {
            if (nSigmaSqPr < nSigmaMaxForDcaxy) {
              registry.fill(HIST("dcaPrOutOfJet"), trkPt, track.dcaXY(), jet.ptCorr, centrality);
            }

            if (nSigmaSqPi < nSigmaMaxForDcaxy) {
              registry.fill(HIST("dcaPiOutOfJet"), trkPt, track.dcaXY(), jet.ptCorr, centrality);
            }
          }

          if (hasTofPi) {
            registry.fill(HIST("tpcTofPiOutOfJet"), trkP, tpcPi, jet.ptCorr, centrality);
            if (isTpcPiRange) {
              registry.fill(HIST("pVsPtForPiOutOfJet"), trkP, trkPt, jet.ptCorr, centrality);
            }
          }
          if (hasTofPr) {
            registry.fill(HIST("tpcTofPrOutOfJet"), trkP, tpcPr, jet.ptCorr, centrality);
            if (isTpcPrRange) {
              registry.fill(HIST("pVsPtForPrOutOfJet"), trkP, trkPt, jet.ptCorr, centrality);
            }
          }
        }

        // --- Signal Fill ---
        registry.fill(HIST("jetDistanceVsTrackpt"), distance, trkPt);
        registry.fill(HIST("jetTpcDedx"), trkP, tpcSig, distance);
        registry.fill(HIST("jetTofBeta"), trkP, beta);

        // dcaXY
        if (track.hasTOF()) {
          if (nSigmaSqPr < nSigmaMaxForDcaxy) {
            registry.fill(HIST("jetDcaPr"), trkPt, track.dcaXY(), distance, jet.ptCorr, centrality);
          }

          if (nSigmaSqPi < nSigmaMaxForDcaxy) {
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
    registry.fill(HIST("eventCounterJet"), 0.5);
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

  void processReco(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::JetTracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::ChargedMCDetectorLevelJets const& jets, aod::McParticles const& mcParticles)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    (void)mcParticles;

    float centrality = collision.centFT0M();
    float rho = collision.rho();

    registry.fill(HIST("eventCounterMc"), 0.5);
    registry.fill(HIST("mcCentralityReco"), centrality);

    struct CachedJet {
      float pt;
      float eta;
      float phi;
      float ptCorr;
    };
    std::vector<CachedJet> cachedJets;
    cachedJets.reserve(jets.size());

    for (const auto& jet : jets) {
      float mcdPtCorr = jet.pt() - rho * jet.area();
      cachedJets.push_back({jet.pt(), jet.eta(), jet.phi(), mcdPtCorr});
      registry.fill(HIST("jetPt"), jet.pt());
    }

    // reco track loop
    for (const auto& track : tracks) {

      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      if (!track.has_mcParticle())
        continue;

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

      auto mcParticle = track.mcParticle();
      registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcParticle.pt());

      if (std::fabs(mcParticle.y()) >= mcRapidityMax)
        continue;

      const int producedByDecay = 4;

      bool isPrimary = mcParticle.isPhysicalPrimary();
      bool isSecondDecay = !isPrimary && (mcParticle.getProcess() == producedByDecay);
      bool isSecondMaterial = !isPrimary && !isSecondDecay;

      int pdg = std::abs(mcParticle.pdgCode());
      bool isPion = (pdg == PDG_t::kPiPlus);
      bool isKaon = (pdg == PDG_t::kKPlus);
      bool isProton = (pdg == PDG_t::kProton);

      if (!isPion && !isKaon && !isProton)
        continue;

      bool hasTof = track.hasTOF();

      for (const auto& jet : cachedJets) {

        float dEta = std::abs(track.eta() - jet.eta);
        if (dEta > distanceMax)
          continue;

        float dPhi = std::abs(track.phi() - jet.phi);
        if (dPhi > o2::constants::math::PI) {
          dPhi = o2::constants::math::TwoPI - dPhi;
        }

        if (dPhi > distanceMax)
          continue;

        float deltaR = std::sqrt(dEta * dEta + dPhi * dPhi);
        if (deltaR > distanceMax)
          continue;

        if (isPrimary) {
          // Tracking
          if (isPion)
            registry.fill(HIST("ptHistogramPion"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);
          else if (isKaon)
            registry.fill(HIST("ptHistogramKaon"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);
          else if (isProton)
            registry.fill(HIST("ptHistogramProton"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);

          // TOF matched
          if (hasTof) {
            if (isPion)
              registry.fill(HIST("ptHistogramPionTof"), mcParticle.pt(), jet.ptCorr, centrality);
            else if (isKaon)
              registry.fill(HIST("ptHistogramKaonTof"), mcParticle.pt(), jet.ptCorr, centrality);
            else if (isProton)
              registry.fill(HIST("ptHistogramProtonTof"), mcParticle.pt(), jet.ptCorr, centrality);
          }
        } else { // Secondary
          if (isSecondDecay) {
            // from Decay
            if (isPion)
              registry.fill(HIST("dcaDecayPion"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);
            else if (isProton)
              registry.fill(HIST("dcaDecayProton"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);
          } else if (isSecondMaterial) {
            // from Material
            if (isPion)
              registry.fill(HIST("dcaMaterialPion"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);
            else if (isProton)
              registry.fill(HIST("dcaMaterialProton"), mcParticle.pt(), track.dcaXY(), jet.ptCorr, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processReco, "process reconstructed simulation information", true);

  void processSim(aod::JetMcCollisions::iterator const& mcCollision, soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, aod::ChargedMCParticleLevelJets const& mcpjets, aod::JetParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) > vertexZCut) {
      return;
    }

    if (collisions.size() == 0) {
      return;
    }

    // --- centrality ---
    float centrality = collisions.begin().centFT0M();
    registry.fill(HIST("mcCentralitySim"), centrality);
    const float maxR2 = distanceMax * distanceMax;

    // --- loop over MC particles only once ---
    for (const auto& mcParticle : mcParticles) {

      // --- early cuts on particle ---
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      if (std::abs(mcParticle.y()) > mcRapidityMax) {
        continue;
      }

      int absPdg = std::abs(mcParticle.pdgCode());
      if (absPdg != PDG_t::kPiPlus &&
          absPdg != PDG_t::kKPlus &&
          absPdg != PDG_t::kProton) {
        continue;
      }

      const float partPt = mcParticle.pt();
      const float partEta = mcParticle.eta();
      const float partPhi = mcParticle.phi();

      // --- loop over jets ---
      for (const auto& mcpjet : mcpjets) {

        // --- delta eta cut first ---
        float dEta = partEta - mcpjet.eta();
        if (std::abs(dEta) > distanceMax) {
          continue;
        }

        // --- delta phi ---
        float dPhi = std::abs(partPhi - mcpjet.phi());
        if (dPhi > o2::constants::math::PI) {
          dPhi = o2::constants::math::TwoPI - dPhi;
        }

        // --- delta R^2 ---
        float dR2 = dEta * dEta + dPhi * dPhi;
        if (dR2 > maxR2) {
          continue;
        }

        const float jetPt = mcpjet.pt();

        // --- histogram fill ---
        if (absPdg == PDG_t::kPiPlus) {
          registry.fill(HIST("ptGeneratedPion"), partPt, jetPt, centrality);
        } else if (absPdg == PDG_t::kKPlus) {
          registry.fill(HIST("ptGeneratedKaon"), partPt, jetPt, centrality);
        } else if (absPdg == PDG_t::kProton) {
          registry.fill(HIST("ptGeneratedProton"), partPt, jetPt, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetShapeTask>(cfgc)}; }
