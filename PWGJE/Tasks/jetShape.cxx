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
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
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
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetShapeTask {

  static constexpr float PionMass = static_cast<float>(o2::constants::physics::MassPionCharged);
  static constexpr float KaonMass = static_cast<float>(o2::constants::physics::MassKaonCharged);
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
  Configurable<float> randomConeDeltaPhiMin{"randomConeDeltaPhiMin", +o2::constants::math::PIThird, "Minimum delta phi used only when randomConeFullAzimuth=false"};
  Configurable<float> randomConeDeltaPhiMax{"randomConeDeltaPhiMax", 2.0f * o2::constants::math::PIThird, "Maximum delta phi used only when randomConeFullAzimuth=false"};
  Configurable<bool> enableRandomCone{"enableRandomCone", true, "Fill one accepted random cone per selected jet"};
  Configurable<bool> randomConeFullAzimuth{"randomConeFullAzimuth", true, "Generate random cones over the full azimuth instead of the old perpendicular-side window"};
  Configurable<int> randomConeMaxAttempts{"randomConeMaxAttempts", 100, "Maximum attempts for each random cone"};
  Configurable<unsigned int> randomConeSeed{"randomConeSeed", 0u, "TRandom3 seed; 0 lets ROOT choose a unique seed"};
  Configurable<bool> fillDetailedQA{"fillDetailedQA", false, "Fill high-volume dE/dx, beta and distance QA histograms"};
  Configurable<bool> fillPidNsigmaQA{"fillPidNsigmaQA", true, "Fill the light centrality / TOF-matching / nSigma_TOF histograms used to validate the normalisation and the PID efficiency on data"};
  struct : ConfigurableGroup {
    Configurable<bool> fillPid2D{"fillPid2D", true, "Fill the dense (p, nSigma_TPC^p, nSigma_TOF^p) map of inclusive tracks with TOF of ALL selected events (no jet bias), without any PID cut"};
    Configurable<bool> fillPid2DJet{"fillPid2DJet", true, "Fill dense (p, nSigma_TPC^p, nSigma_TOF^p, r x jet pT) maps for signal and perpendicular cones, plus the inclusive map of the same events, without any PID cut (shared-shape 2D fit)"};
    Configurable<float> pid2DCentMin{"pid2DCentMin", 0.f, "Lower centrality edge for the pid2D jet maps"};
    Configurable<float> pid2DCentMax{"pid2DCentMax", 10.f, "Upper centrality edge for the pid2D jet maps"};
    Configurable<float> pid2DPMin{"pid2DPMin", 0.8f, "Lower momentum edge of the proton-hypothesis pid2D maps (GeV/c)"};
    Configurable<float> pid2DPMax{"pid2DPMax", 7.4f, "Upper momentum edge of the proton-hypothesis pid2D maps (GeV/c); p_T <= 5 at |eta| < 0.9 needs p <= 7.2"};
    Configurable<bool> fillPid2DPi{"fillPid2DPi", true, "Also fill the pion-hypothesis pid2D maps (nSigma^pi: pi at the origin with unit width, no p smearing of the pion) in [pid2DPiPMin, pid2DPiPMax)"};
    Configurable<float> pid2DPiPMin{"pid2DPiPMin", 0.8f, "Lower momentum edge of the pion-hypothesis pid2D maps (GeV/c)"};
    Configurable<float> pid2DPiPMax{"pid2DPiPMax", 4.0f, "Upper momentum edge of the pion-hypothesis pid2D maps (GeV/c); p_T < 2.8 at |eta| < 0.9 needs p <= 4.0, above that pi is taken from the proton map"};
    Configurable<bool> fillPid2DKa{"fillPid2DKa", true, "Also fill the kaon-hypothesis pid2D maps (nSigma^K) in [pid2DKaPMin, pid2DKaPMax)"};
    Configurable<float> pid2DKaPMin{"pid2DKaPMin", 0.8f, "Lower momentum edge of the kaon-hypothesis pid2D maps (GeV/c)"};
    Configurable<float> pid2DKaPMax{"pid2DKaPMax", 7.4f, "Upper momentum edge of the kaon-hypothesis pid2D maps (GeV/c); 7.4 = K to p_T 5 at |eta| < 0.9 from its own map (K at the origin, no extrapolated K position)"};
    Configurable<bool> fillPid2DMismatch{"fillPid2DMismatch", true, "Fill TOF-mismatch templates: nSigma_TOF^p recomputed with the TOF time of another (random / nearest) TOF track of the same collision"};
    Configurable<bool> fillPid2DTpcOnly{"fillPid2DTpcOnly", true, "Fill the TPC-only (p, nSigma_TPC^p, hasTOF, region) map without TOF requirement (TOF matching of p+K relative to pi on data)"};
    Configurable<bool> fillPid2DMc{"fillPid2DMc", true, "MC: fill the pid2D maps split by true species, origin and TOF-match quality (closure of the 2D fit)"};
    Configurable<float> pid2DMcTofBadNSigma{"pid2DMcTofBadNSigma", 5.f, "MC: flag the TOF hit as mismatched when |nSigma_TOF| under the true mass hypothesis exceeds this"};
    Configurable<float> pid2DTpcOnlyPMax{"pid2DTpcOnlyPMax", 1.6f, "Upper momentum edge of the TPC-only pion map and of the TPC-only cone maps (GeV/c); the TPC separates protons only below ~1.4 GeV/c"};
    Configurable<float> pid2DV0PMin{"pid2DV0PMin", 0.6f, "V0 daughters (processPid2DV0Daughters): lower momentum edge (GeV/c)"};
    Configurable<float> pid2DV0PMax{"pid2DV0PMax", 4.0f, "V0 daughters: upper momentum edge (GeV/c)"};
    Configurable<float> pid2DV0CosPAMin{"pid2DV0CosPAMin", 0.995f, "V0 daughters: minimum cosine of the V0 pointing angle"};
    Configurable<float> pid2DV0DcaDauMax{"pid2DV0DcaDauMax", 1.0f, "V0 daughters: maximum DCA between the two daughters (cm)"};
    Configurable<float> pid2DV0DcaDauToPvMin{"pid2DV0DcaDauToPvMin", 0.1f, "V0 daughters: minimum |DCA| of each daughter to the primary vertex (cm)"};
    Configurable<float> pid2DV0RadiusMin{"pid2DV0RadiusMin", 0.5f, "V0 daughters: minimum decay radius (cm)"};
    Configurable<float> pid2DV0RadiusMax{"pid2DV0RadiusMax", 20.f, "V0 daughters: maximum decay radius (cm); the TOF expected time assumes a primary track, so the maps keep the radius as an axis"};
    Configurable<float> pid2DV0LambdaWindow{"pid2DV0LambdaWindow", 0.006f, "V0 daughters: Lambda mass window half width (GeV/c^2); side bands at 2-3 half widths on both sides"};
    Configurable<float> pid2DV0K0SWindow{"pid2DV0K0SWindow", 0.015f, "V0 daughters: K0S mass window half width (GeV/c^2); side bands at 2-3 half widths on both sides"};
    Configurable<float> pid2DV0CompetingVeto{"pid2DV0CompetingVeto", 0.01f, "V0 daughters: reject Lambda candidates within this of the K0S mass, and K0S candidates within half of it of the Lambda mass (GeV/c^2)"};
  } pid2DOpts;
  Configurable<bool> fillTofMismatchQA{"fillTofMismatchQA", true, "Fill the wide-range nSigma_TOF histograms used to measure the TOF-mismatch pedestal; the data side also needs fillPidNsigmaQA"};
  Configurable<float> nSigmaWideMid{"nSigmaWideMid", 100.0f, "Boundary between the medium and the coarse part of the wide nSigma_TOF axis"};
  Configurable<float> nSigmaWideMax{"nSigmaWideMax", 1000.0f, "Outer edge of the wide nSigma_TOF axis; one bunch crossing is about 625 sigma"};
  Configurable<int> nBinsNSigmaWideMid{"nBinsNSigmaWideMid", 18, "Bins per side between nSigmaMax and nSigmaWideMid"};
  Configurable<int> nBinsNSigmaWideOuter{"nBinsNSigmaWideOuter", 18, "Bins per side between nSigmaWideMid and nSigmaWideMax"};
  Configurable<int> nBinsPForWide{"nBinsPForWide", 80, "Momentum bins of the wide nSigma_TOF histograms; nSigma depends on p, not on p_T, so these have to be fine"};
  Configurable<bool> applyCentralityCut{"applyCentralityCut", true, "Apply centralityMinForCut <= centrality < centralityMaxForCut inside this task"};
  Configurable<bool> useRapiditySelection{"useRapiditySelection", false, "Use species-dependent |y| cut in addition to the detector |eta| acceptance (cross-check only)"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8",
                                            "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "nominal track selection for yields and tracking efficiency"};
  Configurable<bool> skipMbGapEvents{"skipMbGapEvents", true, "reject events whose sub-generator id is the MB gap one"};
  Configurable<bool> applyRctSelection{"applyRctSelection", true, "apply the run-condition-table quality selection (has no meaning in MC)"};
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
  o2::aod::rctsel::RCTFlagsChecker pid2DV0RctChecker{"CBT_hadronPID"}; // same RCT label as jetderiveddatautilities::selectCollision
  std::vector<int> emptyEventSelectionBits;
  int trackSelection = -1;
  int trackSelectionForDca = -1;
  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections);
    trackSelectionForDca = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelectionsForDca));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    randomNumber.SetSeed(randomConeSeed);

    std::vector<double> nSigmaWideEdges;
    {
      const double fineEdge = nSigmaMax.value;
      const double midEdge = nSigmaWideMid.value;
      const double outerEdge = nSigmaWideMax.value;
      const int nFine = nBinsNSigma.value;
      const int nMid = nBinsNSigmaWideMid.value;
      const int nOuter = nBinsNSigmaWideOuter.value;
      nSigmaWideEdges.reserve(2 * (nOuter + nMid) + nFine + 1);
      for (int i = 0; i <= nOuter; ++i) {
        nSigmaWideEdges.push_back(-outerEdge + (outerEdge - midEdge) * i / nOuter);
      }
      for (int i = 1; i <= nMid; ++i) {
        nSigmaWideEdges.push_back(-midEdge + (midEdge - fineEdge) * i / nMid);
      }
      for (int i = 1; i <= nFine; ++i) {
        nSigmaWideEdges.push_back(-fineEdge + 2.0 * fineEdge * i / nFine);
      }
      for (int i = 1; i <= nMid; ++i) {
        nSigmaWideEdges.push_back(fineEdge + (midEdge - fineEdge) * i / nMid);
      }
      for (int i = 1; i <= nOuter; ++i) {
        nSigmaWideEdges.push_back(midEdge + (outerEdge - midEdge) * i / nOuter);
      }
    }
    const AxisSpec axisNSigmaWide{nSigmaWideEdges, "n#sigma_{TOF}"};
    const AxisSpec axisPWide{nBinsPForWide.value, 0., pMax.value, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPtWide{nBinsPt.value, 0., ptMax.value, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisCentWide{nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value, "centrality (%)"};

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

    if (fillPidNsigmaQA) {
      // (1) Normalisation QA: every collision reaching processInclusiveProductionRatio,
      //     BEFORE passCentrality().  A percentile estimator must be flat here.
      registry.add("centralityAllInc", "centFT0M of every collision entering the inclusive process (before the centrality cut);centrality (%);events", HistType::kTH1F, {{110, -5., 105.}});
      registry.add("eventCounterIncAll", "collisions entering the inclusive process, before the centrality cut", HistType::kTH1F, {{1, 0., 1.}});

      // (2) PID-free charged-track yield and the data-driven TOF-matching efficiency.
      registry.add("trackPtAllInc", "all nominal tracks (no PID, no TOF requirement);p_{T} (GeV/c);centrality (%)", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
      registry.add("trackPtTofInc", "nominal tracks with a TOF hit;p_{T} (GeV/c);centrality (%)", HistType::kTH2F, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

      // (3) Species-tagged TOF status: 0 = no TOF, 1 = finite nSigma_TOF, 2 = sentinel (no event time).
      //     Same definition as mcTofNSigmaStatus* on the MC side.
      registry.add("tofStatusPi", "TPC-tagged pion candidates;p_{T} (GeV/c);TOF status;centrality (%)", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {3, 0., 3.}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
      registry.add("tofStatusPr", "TPC-tagged proton candidates;p_{T} (GeV/c);TOF status;centrality (%)", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {3, 0., 3.}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

      // (4) nSigma_TOF of TPC-tagged candidates, WITHOUT the TOF window -> eps_TOFwin on data.
      registry.add("tofNSigmaTaggedPi", "TPC-tagged pion candidates with TOF;p_{T} (GeV/c);n#sigma_{TOF}^{#pi};centrality (%)", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
      registry.add("tofNSigmaTaggedPr", "TPC-tagged proton candidates with TOF;p_{T} (GeV/c);n#sigma_{TOF}^{p};centrality (%)", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    }

    if (fillTofMismatchQA) {
      // (5) Same selection as (4) but on the wide axis: the mismatch pedestal is measured
      //     where there is no signal and interpolated under the core, and the species
      //     contamination peaks (whose positions follow from p_T and the path length)
      //     become visible instead of piling up in the underflow.
      registry.add("tofNSigmaWidePi", "TPC-tagged pion candidates with TOF, wide axis", HistType::kTH3F, {axisPWide, axisNSigmaWide, axisCentWide});
      registry.add("tofNSigmaWidePr", "TPC-tagged proton candidates with TOF, wide axis", HistType::kTH3F, {axisPWide, axisNSigmaWide, axisCentWide});
      // Response matrix p -> p_T for exactly the same tracks, needed to bring the yields
      // extracted in p bins back onto the p_T axis. pVsPtFor* cannot be used for this:
      // it is filled only for tracks that passed both nSigma windows.
      registry.add("pVsPtWidePi", "TPC-tagged pion candidates with TOF;#it{p};#it{p}_{T};centrality (%)", HistType::kTH3F, {axisPWide, axisPtWide, axisCentWide});
      registry.add("pVsPtWidePr", "TPC-tagged proton candidates with TOF;#it{p};#it{p}_{T};centrality (%)", HistType::kTH3F, {axisPWide, axisPtWide, axisCentWide});
    }

    // ---- 2D PID maps (nSigma_TPC^p x nSigma_TOF^p), dense THnF, no PID cut ----
    // Binned in p, not p_T, because both nSigma depend on p. One proton-hypothesis map holds p, K and pi
    // (pi lies at nSigma_TOF^p > -12 for p > 2.8 GeV/c). Dense on purpose: at full statistics a THnSparse
    // of this size is filled to >50 % and costs far more per bin. Centrality is a window
    // (pid2DCentMin/Max), not an axis: a THn keeps under/overflow bins on every axis.
    const AxisSpec axisP2D{pid2DPEdges(), "#it{p} (GeV/#it{c})"};
    const AxisSpec axisTpc2D{22, -4., 7., "n#sigma_{TPC}^{p}"}; // 0.5: < 0.5 % information loss (Asimov)
    const AxisSpec axisTof2D{72, -12., 6., "n#sigma_{TOF}^{p}"};
    const AxisSpec axisRJet2D{28, 0., 28., "r bin + 7 #times jet p_{T} bin"};
    const bool anyPid2D = pid2DOpts.fillPid2D || pid2DOpts.fillPid2DJet;
    // Pion / kaon hypothesis maps: the hypothesis species sits at the origin with unit width at every p,
    // so its yield needs no p-bin smearing model (pi at p < 4, where it is far off in nSigma_TOF^p; K at low p).
    // TOF: pi map [-6, 12] holds K (+27 at p 0.8, +2.4 at 2.8); K map [-10, 12] holds pi (< 0) and p (> 0).
    const AxisSpec axisP2DPi{pid2DPEdges(pid2DOpts.pid2DPiPMin.value, pid2DOpts.pid2DPiPMax.value), "#it{p} (GeV/#it{c})"};
    const AxisSpec axisTpc2DPi{28, -6., 8., "n#sigma_{TPC}^{#pi}"};
    const AxisSpec axisTof2DPi{72, -6., 12., "n#sigma_{TOF}^{#pi}"};
    const AxisSpec axisP2DKa{pid2DPEdges(pid2DOpts.pid2DKaPMin.value, pid2DOpts.pid2DKaPMax.value), "#it{p} (GeV/#it{c})"};
    const AxisSpec axisTpc2DKa{32, -8., 8., "n#sigma_{TPC}^{K}"};
    const AxisSpec axisTof2DKa{88, -10., 12., "n#sigma_{TOF}^{K}"};
    const bool pid2DPiOn = anyPid2D && pid2DOpts.fillPid2DPi;
    const bool pid2DKaOn = anyPid2D && pid2DOpts.fillPid2DKa;

    if (pid2DOpts.fillPid2D) {
      registry.add("pid2DPr", "inclusive tracks with TOF, all events, no PID cut;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p}", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D});
      if (pid2DOpts.fillPid2DPi) {
        registry.add("pid2DPi", "inclusive tracks with TOF, all events, no PID cut;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi}", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi});
      }
      if (pid2DOpts.fillPid2DKa) {
        registry.add("pid2DKa", "inclusive tracks with TOF, all events, no PID cut;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K}", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa});
      }
    }

    if (pid2DOpts.fillPid2DJet) {
      registry.add("pid2DJetPr", "signal cone, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};r bin + 7 jet bin", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, axisRJet2D});
      registry.add("pid2DPerpPr", "perpendicular cones, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};r bin + 7 jet bin", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, axisRJet2D});
      registry.add("pid2DRandPr", "random cones, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};r bin + 7 reference-jet bin", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, axisRJet2D});
      registry.add("pid2DIncPr", "all tracks with TOF in events with an accepted jet;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p}", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D});
      if (pid2DOpts.fillPid2DPi) {
        registry.add("pid2DJetPi", "signal cone, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};r bin + 7 jet bin", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, axisRJet2D});
        registry.add("pid2DPerpPi", "perpendicular cones, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};r bin + 7 jet bin", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, axisRJet2D});
        registry.add("pid2DRandPi", "random cones, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};r bin + 7 reference-jet bin", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, axisRJet2D});
        registry.add("pid2DIncPi", "all tracks with TOF in events with an accepted jet;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi}", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi});
      }
      if (pid2DOpts.fillPid2DKa) {
        registry.add("pid2DJetKa", "signal cone, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K};r bin + 7 jet bin", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, axisRJet2D});
        registry.add("pid2DPerpKa", "perpendicular cones, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K};r bin + 7 jet bin", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, axisRJet2D});
        registry.add("pid2DRandKa", "random cones, tracks with TOF, no PID cut;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K};r bin + 7 reference-jet bin", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, axisRJet2D});
        registry.add("pid2DIncKa", "all tracks with TOF in events with an accepted jet;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K}", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa});
      }
    }

    if (pid2DOpts.fillPid2DMismatch && anyPid2D) {
      registry.add("pid2DMisPr", "TOF-mismatch template;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p} (borrowed time);source + 2 region", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, {10, 0., 10., "source + 2 #times region"}});
      if (pid2DOpts.fillPid2DPi) {
        registry.add("pid2DMisPi", "TOF-mismatch template;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi} (borrowed time);source + 2 region", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, {10, 0., 10., "source + 2 #times region"}});
      }
      if (pid2DOpts.fillPid2DKa) {
        registry.add("pid2DMisKa", "TOF-mismatch template;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K} (borrowed time);source + 2 region", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, {10, 0., 10., "source + 2 #times region"}});
      }
    }

    if (pid2DOpts.fillPid2DTpcOnly && anyPid2D) {
      registry.add("pid2DTpcOnlyPr", "TPC only, no TOF requirement;#it{p};n#sigma_{TPC}^{p};hasTOF;region", HistType::kTHnF, {axisP2D, axisTpc2D, {2, 0., 2., "hasTOF"}, {3, 0., 3., "region"}});
      // Low-p TPC-only yields (p and pi with the same method) and the data TOF matching probability.
      const AxisSpec axisP2DTpcOnly{pid2DPEdges(pid2DOpts.pid2DPMin.value, pid2DOpts.pid2DTpcOnlyPMax.value), "#it{p} (GeV/#it{c})"};
      registry.add("pid2DTpcOnlyPi", "TPC only, no TOF requirement;#it{p};n#sigma_{TPC}^{#pi};hasTOF;region", HistType::kTHnF, {axisP2DTpcOnly, axisTpc2DPi, {2, 0., 2., "hasTOF"}, {3, 0., 3., "region"}});
      registry.add("pid2DTpcOnlyKa", "TPC only, no TOF requirement;#it{p};n#sigma_{TPC}^{K};hasTOF;region", HistType::kTHnF, {axisP2DTpcOnly, axisTpc2DKa, {2, 0., 2., "hasTOF"}, {3, 0., 3., "region"}});
      if (pid2DOpts.fillPid2DJet) {
        const AxisSpec axisTpcOnlyCone{56, 0., 56., "r bin + 7 #times jet p_{T} bin + 28 #times (0 signal, 1 perp)"};
        registry.add("pid2DTpcOnlyConePr", "TPC only, signal / perp cones;#it{p};n#sigma_{TPC}^{p};hasTOF;cone index", HistType::kTHnF, {axisP2DTpcOnly, axisTpc2D, {2, 0., 2., "hasTOF"}, axisTpcOnlyCone});
        registry.add("pid2DTpcOnlyConePi", "TPC only, signal / perp cones;#it{p};n#sigma_{TPC}^{#pi};hasTOF;cone index", HistType::kTHnF, {axisP2DTpcOnly, axisTpc2DPi, {2, 0., 2., "hasTOF"}, axisTpcOnlyCone});
        registry.add("pid2DTpcOnlyConeKa", "TPC only, signal / perp cones;#it{p};n#sigma_{TPC}^{K};hasTOF;cone index", HistType::kTHnF, {axisP2DTpcOnly, axisTpc2DKa, {2, 0., 2., "hasTOF"}, axisTpcOnlyCone});
      }
    }

    if (doprocessPid2DV0Daughters) {
      // Identified daughters (Lambda -> p pi, K0S -> pi pi): TOF matching probability of p / pi on data (note-style bookkeeping).
      // No TOF hit -> nSigma_TOF filled as underflow, and hasTOF = 0 in the status axis.
      const AxisSpec axisP2DV0{pid2DPEdges(pid2DOpts.pid2DV0PMin.value, pid2DOpts.pid2DV0PMax.value), "#it{p} (GeV/#it{c})"};
      const AxisSpec axisV0Status{4, 0., 4., "hasTOF + 2 #times (0 mass window, 1 side bands)"};
      const AxisSpec axisV0Radius{std::vector<double>{0., 3., 8., 100.}, "V0 decay radius (cm)"};
      registry.add("pid2DV0Pr", "protons from Lambda / anti-Lambda;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};status;radius", HistType::kTHnF, {axisP2DV0, axisTpc2D, axisTof2D, axisV0Status, axisV0Radius});
      registry.add("pid2DV0Pi", "pions from K0S;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};status;radius", HistType::kTHnF, {axisP2DV0, axisTpc2DPi, axisTof2DPi, axisV0Status, axisV0Radius});
      registry.add("pid2DV0MassLambda", "Lambda and anti-Lambda candidates;#it{p}_{T} (GeV/#it{c});#it{m} (GeV/#it{c}^{2})", HistType::kTH2F, {{25, 0., 5.}, {140, 1.08, 1.15}});
      registry.add("pid2DV0MassK0S", "K0S candidates;#it{p}_{T} (GeV/#it{c});#it{m} (GeV/#it{c}^{2})", HistType::kTH2F, {{25, 0., 5.}, {160, 0.42, 0.58}});
      registry.add("pid2DV0Events", "events entering processPid2DV0Daughters;centrality (%);events", HistType::kTH1F, {{20, 0., 100.}});
    }

    if (anyPid2D) {
      constexpr double PtMin2D = 0.5; // p_T 0.8-1.0 needs tracks down to p_T 0.8 / cosh(0) and the 0.1 grid
      constexpr double PtBinWidth2D = 0.1;
      auto axisPt2D = [&](float pMax) {
        const int nPt2D = static_cast<int>(std::lround((pMax - PtMin2D) / PtBinWidth2D));
        return AxisSpec{nPt2D, PtMin2D, PtMin2D + nPt2D * PtBinWidth2D, "#it{p}_{T} (GeV/#it{c})"};
      };
      registry.add("pid2DPVsPt", "tracks entering the pid2D maps;#it{p};#it{p}_{T};region", HistType::kTHnF, {axisP2D, axisPt2D(pid2DOpts.pid2DPMax.value), {5, 0., 5., "region"}});
      // Separate p -> p_T matrices: with useRapiditySelection the pi / K maps hold a different track set.
      if (pid2DPiOn) {
        registry.add("pid2DPVsPtPi", "tracks entering the pion-hypothesis pid2D maps;#it{p};#it{p}_{T};region", HistType::kTHnF, {axisP2DPi, axisPt2D(pid2DOpts.pid2DPiPMax.value), {5, 0., 5., "region"}});
      }
      if (pid2DKaOn) {
        registry.add("pid2DPVsPtKa", "tracks entering the kaon-hypothesis pid2D maps;#it{p};#it{p}_{T};region", HistType::kTHnF, {axisP2DKa, axisPt2D(pid2DOpts.pid2DKaPMax.value), {5, 0., 5., "region"}});
      }
    }

    if (pid2DOpts.fillPid2DMc && pid2DOpts.fillPid2DTpcOnly && doprocessEfficiencyAndPurity) {
      // MC closure of the TPC-only counting and the true P_hit / P_correct per species (TOF status 0 none, 1 good, 2 bad)
      const AxisSpec axisP2DTpcOnlyMc{pid2DPEdges(pid2DOpts.pid2DPMin.value, pid2DOpts.pid2DTpcOnlyPMax.value), "#it{p} (GeV/#it{c})"};
      const AxisSpec axisTruthTpcOnlyMc{36, 0., 36., "species + 6 #times (origin + 2 #times TOF status)"};
      registry.add("mcPid2DTpcOnlyPr", "MC: TPC only, all tracks in events with a matched jet;#it{p};n#sigma_{TPC}^{p};truth index", HistType::kTHnF, {axisP2DTpcOnlyMc, axisTpc2D, axisTruthTpcOnlyMc});
      registry.add("mcPid2DTpcOnlyPi", "MC: TPC only, all tracks in events with a matched jet;#it{p};n#sigma_{TPC}^{#pi};truth index", HistType::kTHnF, {axisP2DTpcOnlyMc, axisTpc2DPi, axisTruthTpcOnlyMc});
      registry.add("mcPid2DTpcOnlyKa", "MC: TPC only, all tracks in events with a matched jet;#it{p};n#sigma_{TPC}^{K};truth index", HistType::kTHnF, {axisP2DTpcOnlyMc, axisTpc2DKa, axisTruthTpcOnlyMc});
    }

    if (pid2DOpts.fillPid2DMc && doprocessEfficiencyAndPurity) {
      const AxisSpec axisTruthInc2DMc{24, 0., 24., "species + 6 #times (origin + 2 #times tofBad)"};
      const AxisSpec axisTruthCone2DMc{24, 0., 24., "species4 + 4 #times (tofBad + 2 #times r bin)"};
      registry.add("mcPid2DIncPr", "MC: all tracks with TOF in events with a matched jet;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};truth index", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, axisTruthInc2DMc});
      registry.add("mcPid2DJetPr", "MC: signal cone, tracks with TOF;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};truth index", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, axisTruthCone2DMc});
      registry.add("mcPid2DPerpPr", "MC: perpendicular cones, tracks with TOF;#it{p};n#sigma_{TPC}^{p};n#sigma_{TOF}^{p};truth index", HistType::kTHnF, {axisP2D, axisTpc2D, axisTof2D, axisTruthCone2DMc});
      if (pid2DOpts.fillPid2DPi) {
        registry.add("mcPid2DIncPi", "MC: all tracks with TOF in events with a matched jet;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};truth index", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, axisTruthInc2DMc});
        registry.add("mcPid2DJetPi", "MC: signal cone, tracks with TOF;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};truth index", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, axisTruthCone2DMc});
        registry.add("mcPid2DPerpPi", "MC: perpendicular cones, tracks with TOF;#it{p};n#sigma_{TPC}^{#pi};n#sigma_{TOF}^{#pi};truth index", HistType::kTHnF, {axisP2DPi, axisTpc2DPi, axisTof2DPi, axisTruthCone2DMc});
      }
      if (pid2DOpts.fillPid2DKa) {
        registry.add("mcPid2DIncKa", "MC: all tracks with TOF in events with a matched jet;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K};truth index", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, axisTruthInc2DMc});
        registry.add("mcPid2DJetKa", "MC: signal cone, tracks with TOF;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K};truth index", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, axisTruthCone2DMc});
        registry.add("mcPid2DPerpKa", "MC: perpendicular cones, tracks with TOF;#it{p};n#sigma_{TPC}^{K};n#sigma_{TOF}^{K};truth index", HistType::kTHnF, {axisP2DKa, axisTpc2DKa, axisTof2DKa, axisTruthCone2DMc});
      }
    }

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
    registry.add("mcCutFlow", "processEfficiencyAndPurity cut flow: 0 = called, 1 = selectCollision, 2 = has_mcCollision, 3 = passCentrality", HistType::kTH1F, {{4, 0., 4.}});
    registry.add("mcSelDetail", "selectCollision breakdown: 0 = called, 1 = not MB gap, 2 = RCT, 3 = event selection bits", HistType::kTH1F, {{4, 0., 4.}});
    registry.add("mcSubGeneratorId", "sub-generator id of the associated MC collision", HistType::kTH1F, {{6, -2.5, 3.5}});
    registry.add("mcCentralityRaw", "centFT0M seen by processEfficiencyAndPurity, before the cut", HistType::kTH1F, {{124, -12., 112.}});

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

    // nSigma window efficiency breakdown (inclusive only, ~0.3 MB each)
    registry.add("effNumKinePiInc", "True primary pion, TOF matched, passing the species kinematics (denominator of the nSigma windows)", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcWinPiInc", "effNumKinePiInc and the |nSigma_TPC^pi| window ONLY", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofWinPiInc", "effNumKinePiInc and the |nSigma_TOF^pi| window ONLY", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumKinePrInc", "True primary proton, TOF matched, passing the species kinematics (denominator of the nSigma windows)", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTpcWinPrInc", "effNumKinePrInc and the |nSigma_TPC^p| window ONLY", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("effNumTofWinPrInc", "effNumKinePrInc and the |nSigma_TOF^p| window ONLY", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});

    // nSigma response of true primary pions and protons in MC (~1.5 MB each)
    registry.add("mcTpcNSigmaPi", "True primary pion: p vs nSigma_TPC^pi", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("mcTofNSigmaPi", "True primary pion: pT vs nSigma_TOF^pi", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("mcTpcNSigmaPr", "True primary proton: p vs nSigma_TPC^p", HistType::kTHnSparseD, {{nBinsPForCut.value, 0, pMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("mcTofNSigmaPr", "True primary proton: pT vs nSigma_TOF^p", HistType::kTHnSparseD, {{nBinsPt.value, 0, ptMax.value}, {nBinsNSigma.value, nSigmaMin.value, nSigmaMax.value}, {nBinsCentrality.value, centralityMinForCut.value, centralityMaxForCut.value}});
    registry.add("mcTofNSigmaStatusPi", "True primary pion: 0 = no TOF, 1 = finite nSigma_TOF, 2 = unphysical nSigma_TOF", HistType::kTH1F, {{3, 0., 3.}});
    registry.add("mcTofNSigmaStatusPr", "True primary proton: 0 = no TOF, 1 = finite nSigma_TOF, 2 = unphysical nSigma_TOF", HistType::kTH1F, {{3, 0., 3.}});

    if (fillTofMismatchQA) {
      registry.add("mcTofNSigmaWidePi", "True primary pion with TOF, wide axis", HistType::kTH3F, {axisPWide, axisNSigmaWide, axisCentWide});
      registry.add("mcTofNSigmaWidePr", "True primary proton with TOF, wide axis", HistType::kTH3F, {axisPWide, axisNSigmaWide, axisCentWide});
      registry.add("mcTofNSigmaWideKa", "True primary kaon with TOF, wide axis", HistType::kTH3F, {axisPWide, axisNSigmaWide, axisCentWide});
    }

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

  std::vector<double> pid2DPEdges() const
  {
    return pid2DPEdges(pid2DOpts.pid2DPMin.value, pid2DOpts.pid2DPMax.value);
  }

  // 0.2 GeV/c up to 5 GeV/c, 0.4 GeV/c above (common to all hypotheses, so the maps share bin edges).
  static std::vector<double> pid2DPEdges(double pMin, double pMax)
  {
    constexpr double FineWidth = 0.2;
    constexpr double CoarseWidth = 0.4;
    constexpr double FineUpTo = 5.0;
    constexpr double Tolerance = 1e-6;
    constexpr double Rounding = 1000.;
    pMin = std::round(pMin * Rounding) / Rounding;
    pMax = std::round(pMax * Rounding) / Rounding;
    std::vector<double> edges{pMin};
    double p = pMin;
    while (p < pMax - Tolerance) {
      p += (p < FineUpTo - Tolerance) ? FineWidth : CoarseWidth;
      p = std::min(std::round(p * Rounding) / Rounding, pMax);
      edges.push_back(p);
    }
    return edges;
  }

  // TOF hits of one collision, used to build the mismatch templates (8).
  struct TofHit2D {
    int64_t index;
    float eta;
    float phi;
    float time; // t - t0 = L / (beta c), ps
  };
  static constexpr float SpeedOfLightCmPerPs = 0.0299792458f;

  template <typename TrackT>
  static bool tofTimeOfFlight(TrackT const& track, float& time)
  {
    if (!track.hasTOF() || track.beta() <= 0.f || track.length() <= 0.f) {
      return false;
    }
    time = track.length() / (track.beta() * SpeedOfLightCmPerPs);
    return true;
  }

  template <typename FullTrackT, typename JetTracksT>
  void collectTofHits(JetTracksT const& tracks, std::vector<TofHit2D>& hits)
  {
    hits.clear();
    for (const auto& jetTrack : tracks) {
      const auto track = jetTrack.template track_as<FullTrackT>();
      float time = 0.f;
      if (tofTimeOfFlight(track, time)) {
        hits.push_back({track.globalIndex(), track.eta(), track.phi(), time});
      }
    }
  }

  // Time difference (ps) between a borrowed TOF hit of the same collision and the track's own hit:
  // [0] a random other hit, [1] the nearest other hit in (eta, phi) at the vertex. NaN if not available.
  // The same borrowed hit is used for all mass hypotheses, so the p, K and pi templates are consistent.
  template <typename TrackT>
  std::array<float, 2> mismatchDeltaT(TrackT const& track, std::vector<TofHit2D> const& hits)
  {
    constexpr int MaxAttempts = 10;
    constexpr size_t MinHits = 2;
    std::array<float, 2> out{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    float time = 0.f;
    if (hits.size() < MinHits || !tofTimeOfFlight(track, time)) {
      return out;
    }
    const int64_t self = track.globalIndex();
    for (int attempt = 0; attempt < MaxAttempts; ++attempt) {
      const auto& hit = hits[randomNumber.Integer(static_cast<uint32_t>(hits.size()))];
      if (hit.index != self) {
        out[0] = hit.time - time;
        break;
      }
    }
    float best = std::numeric_limits<float>::max();
    const TofHit2D* nearest = nullptr;
    for (const auto& hit : hits) {
      if (hit.index == self) {
        continue;
      }
      const float dEta = hit.eta - track.eta();
      const float dPhi = absDeltaPhi(hit.phi, track.phi());
      const float d2 = dEta * dEta + dPhi * dPhi;
      if (d2 < best) {
        best = d2;
        nearest = &hit;
      }
    }
    if (nearest != nullptr) {
      out[1] = nearest->time - time;
    }
    return out;
  }

  // nSigma_TOF^h with the borrowed time: nSigma_h + dt / sigma_h (sigma_h = tofExpSigma of hypothesis h).
  static std::array<float, 2> mismatchNSigma(std::array<float, 2> const& dt, float nSigma, float sigma)
  {
    if (!(sigma > 0.f)) {
      return {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    }
    return {nSigma + dt[0] / sigma, nSigma + dt[1] / sigma};
  }

  // Regions of the pid2D maps (int constants, not an enum: they are combined with float bin centres).
  static constexpr int RegIncAll = 0; // inclusive, all events
  static constexpr int RegIncJet = 1; // inclusive, events with an accepted jet
  static constexpr int RegSignal = 2;
  static constexpr int RegPerp = 3;
  static constexpr int RegRand = 4;

  // Mass hypotheses of the pid2D maps.
  static constexpr int Hyp2DPr = 0;
  static constexpr int Hyp2DPi = 1;
  static constexpr int Hyp2DKa = 2;

  template <int Hyp = Hyp2DPr>
  void fillMismatchTemplate(float p, float tpc, std::array<float, 2> const& mis, int region)
  {
    constexpr int NSources = 2;
    for (int src = 0; src < NSources; ++src) {
      if (!std::isfinite(mis[src])) {
        continue;
      }
      const float index = src + NSources * region + 0.5f;
      if constexpr (Hyp == Hyp2DPi) {
        registry.fill(HIST("pid2DMisPi"), p, tpc, mis[src], index);
      } else if constexpr (Hyp == Hyp2DKa) {
        registry.fill(HIST("pid2DMisKa"), p, tpc, mis[src], index);
      } else {
        registry.fill(HIST("pid2DMisPr"), p, tpc, mis[src], index);
      }
    }
  }

  // Folded (r, jet pT) index of the pid2D jet maps; -1 if outside.
  static int rJetIndex2D(float r, float jetPt)
  {
    constexpr int NRBins2D = 7;
    constexpr int NJetPtBins2D = 4;
    constexpr float RBinWidth2D = 0.1f;
    constexpr float RMax2D = NRBins2D * RBinWidth2D;
    constexpr int Outside = -1;
    static constexpr std::array<float, NJetPtBins2D + 1> JetPtEdges2D = {20.f, 40.f, 60.f, 80.f, 140.f};
    if (r < 0.f || r >= RMax2D) {
      return Outside;
    }
    const int ir = std::min(NRBins2D - 1, static_cast<int>(r / RBinWidth2D));
    for (int k = 0; k < NJetPtBins2D; ++k) {
      if (jetPt >= JetPtEdges2D[k] && jetPt < JetPtEdges2D[k + 1]) {
        return ir + NRBins2D * k;
      }
    }
    return Outside;
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
          rcPhi = randomNumber.Uniform(0.0f, o2::constants::math::TwoPI);
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

  using FullTrackInfo = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta>;
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

    std::vector<TofHit2D> tofHits;
    if (pid2DOpts.fillPid2DJet && pid2DOpts.fillPid2DMismatch) {
      collectTofHits<FullTrackInfo>(tracks, tofHits);
    }

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

      const bool inPid2DWindow = pid2DOpts.fillPid2DJet && passesNominal && passPrKinematics &&
                                 centrality >= pid2DOpts.pid2DCentMin && centrality < pid2DOpts.pid2DCentMax &&
                                 trkP >= pid2DOpts.pid2DPMin && trkP < pid2DOpts.pid2DPMax;
      const bool fillPid2DTrack = inPid2DWindow && track.hasTOF();
      const bool fillPid2DTpcOnlyTrack = inPid2DWindow && pid2DOpts.fillPid2DTpcOnly;
      const float hasTofBin = track.hasTOF() ? 1.5f : 0.5f;

      // Pion / kaon hypothesis maps: same track and event selection, own p window and kinematics.
      const bool inPid2DOther = pid2DOpts.fillPid2DJet && passesNominal && track.hasTOF() &&
                                centrality >= pid2DOpts.pid2DCentMin && centrality < pid2DOpts.pid2DCentMax;
      const bool fillPid2DPiTrack = inPid2DOther && pid2DOpts.fillPid2DPi && passPiKinematics &&
                                    trkP >= pid2DOpts.pid2DPiPMin && trkP < pid2DOpts.pid2DPiPMax;
      const bool fillPid2DKaTrack = inPid2DOther && pid2DOpts.fillPid2DKa && passSpeciesKinematics(trkPt, trkEta, KaonMass) &&
                                    trkP >= pid2DOpts.pid2DKaPMin && trkP < pid2DOpts.pid2DKaPMax;
      const float tpcKa = track.tpcNSigmaKa();
      const float tofKa = track.tofNSigmaKa();

      std::array<float, 2> misNSigma{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
      std::array<float, 2> misNSigmaPi = misNSigma;
      std::array<float, 2> misNSigmaKa = misNSigma;
      if (pid2DOpts.fillPid2DMismatch && (fillPid2DTrack || fillPid2DPiTrack || fillPid2DKaTrack)) {
        const auto misDt = mismatchDeltaT(track, tofHits);
        misNSigma = mismatchNSigma(misDt, tofPr, track.tofExpSigmaPr());
        misNSigmaPi = mismatchNSigma(misDt, tofPi, track.tofExpSigmaPi());
        misNSigmaKa = mismatchNSigma(misDt, tofKa, track.tofExpSigmaKa());
      }
      if (fillPid2DTrack) {
        registry.fill(HIST("pid2DIncPr"), trkP, tpcPr, tofPr);
        registry.fill(HIST("pid2DPVsPt"), trkP, trkPt, RegIncJet + 0.5f);
        if (pid2DOpts.fillPid2DMismatch) {
          fillMismatchTemplate<Hyp2DPr>(trkP, tpcPr, misNSigma, RegIncJet);
        }
      }
      if (fillPid2DPiTrack) {
        registry.fill(HIST("pid2DIncPi"), trkP, tpcPi, tofPi);
        registry.fill(HIST("pid2DPVsPtPi"), trkP, trkPt, RegIncJet + 0.5f);
        if (pid2DOpts.fillPid2DMismatch) {
          fillMismatchTemplate<Hyp2DPi>(trkP, tpcPi, misNSigmaPi, RegIncJet);
        }
      }
      if (fillPid2DKaTrack) {
        registry.fill(HIST("pid2DIncKa"), trkP, tpcKa, tofKa);
        registry.fill(HIST("pid2DPVsPtKa"), trkP, trkPt, RegIncJet + 0.5f);
        if (pid2DOpts.fillPid2DMismatch) {
          fillMismatchTemplate<Hyp2DKa>(trkP, tpcKa, misNSigmaKa, RegIncJet);
        }
      }
      // Cone maps (signal / perp / random cone) of the pion and kaon hypotheses.
      auto fillPid2DPiKaCone = [&](int region, float r, float jetPtForIndex) {
        if (!fillPid2DPiTrack && !fillPid2DKaTrack) {
          return;
        }
        const int idx = rJetIndex2D(r, jetPtForIndex);
        if (idx < 0) {
          return;
        }
        const float bin = idx + 0.5f;
        if (fillPid2DPiTrack) {
          if (region == RegSignal) {
            registry.fill(HIST("pid2DJetPi"), trkP, tpcPi, tofPi, bin);
          } else if (region == RegPerp) {
            registry.fill(HIST("pid2DPerpPi"), trkP, tpcPi, tofPi, bin);
          } else {
            registry.fill(HIST("pid2DRandPi"), trkP, tpcPi, tofPi, bin);
          }
          registry.fill(HIST("pid2DPVsPtPi"), trkP, trkPt, region + 0.5f);
          if (pid2DOpts.fillPid2DMismatch) {
            fillMismatchTemplate<Hyp2DPi>(trkP, tpcPi, misNSigmaPi, region);
          }
        }
        if (fillPid2DKaTrack) {
          if (region == RegSignal) {
            registry.fill(HIST("pid2DJetKa"), trkP, tpcKa, tofKa, bin);
          } else if (region == RegPerp) {
            registry.fill(HIST("pid2DPerpKa"), trkP, tpcKa, tofKa, bin);
          } else {
            registry.fill(HIST("pid2DRandKa"), trkP, tpcKa, tofKa, bin);
          }
          registry.fill(HIST("pid2DPVsPtKa"), trkP, trkPt, region + 0.5f);
          if (pid2DOpts.fillPid2DMismatch) {
            fillMismatchTemplate<Hyp2DKa>(trkP, tpcKa, misNSigmaKa, region);
          }
        }
      };

      // TPC-only cone maps (p and pi hypotheses, with and without TOF hit) for the low-p_T TPC-only yields.
      const bool inTpcOnlyCone = pid2DOpts.fillPid2DJet && pid2DOpts.fillPid2DTpcOnly && passesNominal &&
                                 centrality >= pid2DOpts.pid2DCentMin && centrality < pid2DOpts.pid2DCentMax &&
                                 trkP >= pid2DOpts.pid2DPMin && trkP < pid2DOpts.pid2DTpcOnlyPMax;
      const bool fillTpcOnlyConePr = inTpcOnlyCone && passPrKinematics;
      const bool fillTpcOnlyConePi = inTpcOnlyCone && passPiKinematics;
      const bool fillTpcOnlyConeKa = inTpcOnlyCone && passSpeciesKinematics(trkPt, trkEta, KaonMass);
      auto fillTpcOnlyCone = [&](int region, float r, float jetPtForIndex) {
        if (!fillTpcOnlyConePr && !fillTpcOnlyConePi && !fillTpcOnlyConeKa) {
          return;
        }
        const int idx = rJetIndex2D(r, jetPtForIndex);
        if (idx < 0) {
          return;
        }
        constexpr int NCells2D = 28;
        const float bin = idx + NCells2D * (region == RegSignal ? 0 : 1) + 0.5f;
        if (fillTpcOnlyConePr) {
          registry.fill(HIST("pid2DTpcOnlyConePr"), trkP, tpcPr, hasTofBin, bin);
        }
        if (fillTpcOnlyConePi) {
          registry.fill(HIST("pid2DTpcOnlyConePi"), trkP, tpcPi, hasTofBin, bin);
          registry.fill(HIST("pid2DTpcOnlyPi"), trkP, tpcPi, hasTofBin, region == RegSignal ? 1.5f : 2.5f);
        }
        if (fillTpcOnlyConeKa) {
          registry.fill(HIST("pid2DTpcOnlyConeKa"), trkP, tpcKa, hasTofBin, bin);
          registry.fill(HIST("pid2DTpcOnlyKa"), trkP, tpcKa, hasTofBin, region == RegSignal ? 1.5f : 2.5f);
        }
      };

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

          if (fillPid2DTpcOnlyTrack || fillPid2DTrack) {
            const int idx2D = rJetIndex2D(distBg, jet.ptCorr);
            if (idx2D >= 0) {
              if (fillPid2DTpcOnlyTrack) {
                registry.fill(HIST("pid2DTpcOnlyPr"), trkP, tpcPr, hasTofBin, 2.5f);
              }
              if (fillPid2DTrack) {
                registry.fill(HIST("pid2DPerpPr"), trkP, tpcPr, tofPr, idx2D + 0.5f);
                registry.fill(HIST("pid2DPVsPt"), trkP, trkPt, RegPerp + 0.5f);
                if (pid2DOpts.fillPid2DMismatch) {
                  fillMismatchTemplate<Hyp2DPr>(trkP, tpcPr, misNSigma, RegPerp);
                }
              }
            }
          }
          fillPid2DPiKaCone(RegPerp, distBg, jet.ptCorr);
          fillTpcOnlyCone(RegPerp, distBg, jet.ptCorr);

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

        if (fillPid2DTpcOnlyTrack || fillPid2DTrack) {
          const int idx2D = rJetIndex2D(distance, jet.ptCorr);
          if (idx2D >= 0) {
            if (fillPid2DTpcOnlyTrack) {
              registry.fill(HIST("pid2DTpcOnlyPr"), trkP, tpcPr, hasTofBin, 1.5f);
            }
            if (fillPid2DTrack) {
              registry.fill(HIST("pid2DJetPr"), trkP, tpcPr, tofPr, idx2D + 0.5f);
              registry.fill(HIST("pid2DPVsPt"), trkP, trkPt, RegSignal + 0.5f);
              if (pid2DOpts.fillPid2DMismatch) {
                fillMismatchTemplate<Hyp2DPr>(trkP, tpcPr, misNSigma, RegSignal);
              }
            }
          }
        }
        fillPid2DPiKaCone(RegSignal, distance, jet.ptCorr);
        fillTpcOnlyCone(RegSignal, distance, jet.ptCorr);

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

        if (fillPid2DTrack) {
          const int idx2D = rJetIndex2D(distRC, cone.referenceJetPtCorr);
          if (idx2D >= 0) {
            registry.fill(HIST("pid2DRandPr"), trkP, tpcPr, tofPr, idx2D + 0.5f);
            registry.fill(HIST("pid2DPVsPt"), trkP, trkPt, RegRand + 0.5f);
            if (pid2DOpts.fillPid2DMismatch) {
              fillMismatchTemplate<Hyp2DPr>(trkP, tpcPr, misNSigma, RegRand);
            }
          }
        }
        fillPid2DPiKaCone(RegRand, distRC, cone.referenceJetPtCorr);

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
    if (fillPidNsigmaQA) {
      registry.fill(HIST("centralityAllInc"), centrality);
      registry.fill(HIST("eventCounterIncAll"), 0.5);
    }
    if (!passCentrality(centrality)) {
      return;
    }
    registry.fill(HIST("eventCounterInc"), centrality);

    std::vector<TofHit2D> tofHits;
    if (pid2DOpts.fillPid2D && pid2DOpts.fillPid2DMismatch) {
      collectTofHits<FullTrackInfo>(tracks, tofHits);
    }

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

      if (fillPidNsigmaQA) {
        constexpr float MaxSaneNSigma = 10.f;
        constexpr float NoTofBin = 0.5f;
        constexpr float FiniteNSigmaBin = 1.5f;
        constexpr float UnphysicalNSigmaBin = 2.5f;
        const bool trackHasTof = track.hasTOF();

        registry.fill(HIST("trackPtAllInc"), track.pt(), centrality);
        if (trackHasTof) {
          registry.fill(HIST("trackPtTofInc"), track.pt(), centrality);
        }
        if (passPiKinematics && isTpcPiRange) {
          registry.fill(HIST("tofStatusPi"), track.pt(),
                        !trackHasTof ? NoTofBin
                                     : (std::abs(track.tofNSigmaPi()) < MaxSaneNSigma ? FiniteNSigmaBin : UnphysicalNSigmaBin),
                        centrality);
          if (trackHasTof) {
            registry.fill(HIST("tofNSigmaTaggedPi"), track.pt(), track.tofNSigmaPi(), centrality);
            if (fillTofMismatchQA) {
              registry.fill(HIST("tofNSigmaWidePi"), track.p(), track.tofNSigmaPi(), centrality);
              registry.fill(HIST("pVsPtWidePi"), track.p(), track.pt(), centrality);
            }
          }
        }
        if (passPrKinematics && isTpcPrRange) {
          registry.fill(HIST("tofStatusPr"), track.pt(),
                        !trackHasTof ? NoTofBin
                                     : (std::abs(track.tofNSigmaPr()) < MaxSaneNSigma ? FiniteNSigmaBin : UnphysicalNSigmaBin),
                        centrality);
          if (trackHasTof) {
            registry.fill(HIST("tofNSigmaTaggedPr"), track.pt(), track.tofNSigmaPr(), centrality);
            if (fillTofMismatchQA) {
              registry.fill(HIST("tofNSigmaWidePr"), track.p(), track.tofNSigmaPr(), centrality);
              registry.fill(HIST("pVsPtWidePr"), track.p(), track.pt(), centrality);
            }
          }
        }
      }

      // pid2D inclusive maps (all events, no jet bias), proton / pion / kaon hypotheses.
      // One borrowed TOF hit per track serves the mismatch templates of all hypotheses.
      if (pid2DOpts.fillPid2D && centrality >= pid2DOpts.pid2DCentMin && centrality < pid2DOpts.pid2DCentMax) {
        const float trkP = track.p();
        const bool inPr = passPrKinematics && trkP >= pid2DOpts.pid2DPMin && trkP < pid2DOpts.pid2DPMax;
        const bool inPi = pid2DOpts.fillPid2DPi && track.hasTOF() && passPiKinematics && trkP >= pid2DOpts.pid2DPiPMin && trkP < pid2DOpts.pid2DPiPMax;
        const bool inKa = pid2DOpts.fillPid2DKa && track.hasTOF() && passSpeciesKinematics(track.pt(), track.eta(), KaonMass) &&
                          trkP >= pid2DOpts.pid2DKaPMin && trkP < pid2DOpts.pid2DKaPMax;
        std::array<float, 2> misDt{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
        if (pid2DOpts.fillPid2DMismatch && ((inPr && track.hasTOF()) || inPi || inKa)) {
          misDt = mismatchDeltaT(track, tofHits);
        }
        if (inPr) {
          if (pid2DOpts.fillPid2DTpcOnly) {
            registry.fill(HIST("pid2DTpcOnlyPr"), trkP, track.tpcNSigmaPr(), track.hasTOF() ? 1.5f : 0.5f, 0.5f);
          }
          if (track.hasTOF()) {
            registry.fill(HIST("pid2DPr"), trkP, track.tpcNSigmaPr(), track.tofNSigmaPr());
            registry.fill(HIST("pid2DPVsPt"), trkP, track.pt(), RegIncAll + 0.5f);
            if (pid2DOpts.fillPid2DMismatch) {
              fillMismatchTemplate<Hyp2DPr>(trkP, track.tpcNSigmaPr(), mismatchNSigma(misDt, track.tofNSigmaPr(), track.tofExpSigmaPr()), RegIncAll);
            }
          }
        }
        if (pid2DOpts.fillPid2DTpcOnly && trkP >= pid2DOpts.pid2DPMin && trkP < pid2DOpts.pid2DTpcOnlyPMax) {
          const float hasTofBinInc = track.hasTOF() ? 1.5f : 0.5f;
          if (passPiKinematics) {
            registry.fill(HIST("pid2DTpcOnlyPi"), trkP, track.tpcNSigmaPi(), hasTofBinInc, 0.5f);
          }
          if (passSpeciesKinematics(track.pt(), track.eta(), KaonMass)) {
            registry.fill(HIST("pid2DTpcOnlyKa"), trkP, track.tpcNSigmaKa(), hasTofBinInc, 0.5f);
          }
        }
        if (inPi) {
          registry.fill(HIST("pid2DPi"), trkP, track.tpcNSigmaPi(), track.tofNSigmaPi());
          registry.fill(HIST("pid2DPVsPtPi"), trkP, track.pt(), RegIncAll + 0.5f);
          if (pid2DOpts.fillPid2DMismatch) {
            fillMismatchTemplate<Hyp2DPi>(trkP, track.tpcNSigmaPi(), mismatchNSigma(misDt, track.tofNSigmaPi(), track.tofExpSigmaPi()), RegIncAll);
          }
        }
        if (inKa) {
          registry.fill(HIST("pid2DKa"), trkP, track.tpcNSigmaKa(), track.tofNSigmaKa());
          registry.fill(HIST("pid2DPVsPtKa"), trkP, track.pt(), RegIncAll + 0.5f);
          if (pid2DOpts.fillPid2DMismatch) {
            fillMismatchTemplate<Hyp2DKa>(trkP, track.tpcNSigmaKa(), mismatchNSigma(misDt, track.tofNSigmaKa(), track.tofExpSigmaKa()), RegIncAll);
          }
        }
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

  // ---- Lambda -> p pi and K0S -> pi pi daughters: TOF matching probability of identified p / pi on data.
  using V0DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr,
                                aod::pidTOFbeta, aod::TrackSelection>;
  using V0Collisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>>;

  // 0 = mass window (|m - m0| < w), 1 = side bands (2w <= |m - m0| < 3w; same total width), -1 = neither
  static int v0MassRegion(float m, float m0, float halfWidth)
  {
    constexpr float SideLo = 2.f;
    constexpr float SideHi = 3.f;
    const float d = std::abs(m - m0);
    if (d < halfWidth) {
      return 0;
    }
    if (d >= SideLo * halfWidth && d < SideHi * halfWidth) {
      return 1;
    }
    return -1;
  }

  // Same track quality as the analysis tracks, without the DCA cut (the daughters are secondaries).
  template <bool IsProton, typename TrackT>
  void fillPid2DV0Daughter(TrackT const& dau, int massRegion, float radius)
  {
    if (!passCommonTrackQuality(dau) || !dau.isQualityTrackITS() || !dau.isQualityTrackTPC()) {
      return;
    }
    const float p = dau.p();
    if (p < pid2DOpts.pid2DV0PMin || p >= pid2DOpts.pid2DV0PMax) {
      return;
    }
    constexpr float NoTofValue = -999.f;
    const bool hasTof = dau.hasTOF();
    const float status = (hasTof ? 1.f : 0.f) + 2.f * massRegion + 0.5f;
    if constexpr (IsProton) {
      registry.fill(HIST("pid2DV0Pr"), p, dau.tpcNSigmaPr(), hasTof ? dau.tofNSigmaPr() : NoTofValue, status, radius);
    } else {
      registry.fill(HIST("pid2DV0Pi"), p, dau.tpcNSigmaPi(), hasTof ? dau.tofNSigmaPi() : NoTofValue, status, radius);
    }
  }

  // Same event selection as the JE collisions of the other processes (eventSelections bits and the RCT flags),
  // rebuilt from the EvSels columns: the TOF matching probability depends on the pileup / occupancy selection.
  template <typename CollisionT>
  bool passPid2DV0EventSelection(CollisionT const& collision)
  {
    if (applyRctSelection && !pid2DV0RctChecker.checkTable(collision)) {
      return false;
    }
    const uint16_t bits = jetderiveddatautilities::setEventSelectionBit(collision);
    bool isOrCondition = false;
    for (const auto& bit : eventSelectionBits) {
      if (bit < 0) {
        isOrCondition = true;
        continue;
      }
      const bool isSet = (bits & (1U << bit)) != 0;
      if (!isOrCondition && !isSet) {
        return false;
      }
      if (isOrCondition && isSet) {
        return true;
      }
    }
    return !isOrCondition;
  }

  void processPid2DV0Daughters(V0Collisions::iterator const& collision, V0DauTracks const&, aod::V0Datas const& v0s)
  {
    if (!passPid2DV0EventSelection(collision)) {
      return;
    }
    const float centrality = collision.centFT0M();
    if (centrality < pid2DOpts.pid2DCentMin || centrality >= pid2DOpts.pid2DCentMax) {
      return;
    }
    registry.fill(HIST("pid2DV0Events"), centrality);
    constexpr float MassLambda = o2::constants::physics::MassLambda0;
    constexpr float MassK0S = o2::constants::physics::MassK0Short;
    constexpr float PionDauTpcNSigmaMax = 5.f; // the pion of a Lambda is slow: clean TPC selection, no bias on the proton
    constexpr float LambdaVetoFraction = 0.5f; // K0S candidates are vetoed within half of the window around the Lambda mass
    const float veto = pid2DOpts.pid2DV0CompetingVeto;
    for (const auto& v0 : v0s) {
      if (v0.v0cosPA() < pid2DOpts.pid2DV0CosPAMin || v0.dcaV0daughters() > pid2DOpts.pid2DV0DcaDauMax) {
        continue;
      }
      const float radius = v0.v0radius();
      if (radius < pid2DOpts.pid2DV0RadiusMin || radius > pid2DOpts.pid2DV0RadiusMax) {
        continue;
      }
      if (std::abs(v0.dcapostopv()) < pid2DOpts.pid2DV0DcaDauToPvMin || std::abs(v0.dcanegtopv()) < pid2DOpts.pid2DV0DcaDauToPvMin) {
        continue;
      }
      const auto pos = v0.template posTrack_as<V0DauTracks>();
      const auto neg = v0.template negTrack_as<V0DauTracks>();
      const bool k0sLike = std::abs(v0.mK0Short() - MassK0S) < veto;
      const bool lambdaLike = std::abs(v0.mLambda() - MassLambda) < LambdaVetoFraction * veto || std::abs(v0.mAntiLambda() - MassLambda) < LambdaVetoFraction * veto;
      if (!k0sLike) {
        // Lambda: proton = positive daughter, pi- = negative. Anti-Lambda: antiproton = negative daughter.
        if (std::abs(neg.tpcNSigmaPi()) < PionDauTpcNSigmaMax) {
          registry.fill(HIST("pid2DV0MassLambda"), v0.pt(), v0.mLambda());
          const int region = v0MassRegion(v0.mLambda(), MassLambda, pid2DOpts.pid2DV0LambdaWindow);
          if (region >= 0) {
            fillPid2DV0Daughter<true>(pos, region, radius);
          }
        }
        if (std::abs(pos.tpcNSigmaPi()) < PionDauTpcNSigmaMax) {
          registry.fill(HIST("pid2DV0MassLambda"), v0.pt(), v0.mAntiLambda());
          const int region = v0MassRegion(v0.mAntiLambda(), MassLambda, pid2DOpts.pid2DV0LambdaWindow);
          if (region >= 0) {
            fillPid2DV0Daughter<true>(neg, region, radius);
          }
        }
      }
      if (!lambdaLike) {
        registry.fill(HIST("pid2DV0MassK0S"), v0.pt(), v0.mK0Short());
        const int region = v0MassRegion(v0.mK0Short(), MassK0S, pid2DOpts.pid2DV0K0SWindow);
        if (region >= 0) {
          fillPid2DV0Daughter<false>(pos, region, radius);
          fillPid2DV0Daughter<false>(neg, region, radius);
        }
      }
    }
  }
  PROCESS_SWITCH(JetShapeTask, processPid2DV0Daughters, "Lambda / K0S daughters: data TOF matching probability of identified p / pi (needs the V0 tables, e.g. from the propagation service)", false);

  // MC: TPC-only maps of all reconstructed tracks (with and without TOF) with the true species, origin and TOF status.
  template <typename JetTrackT, typename TrackT>
  void fillPid2DMcTpcOnlyTrack(JetTrackT const& jetTrack, TrackT const& track)
  {
    constexpr int SpPi = 0;
    constexpr int SpKa = 1;
    constexpr int SpPr = 2;
    constexpr int SpEl = 3;
    constexpr int SpMu = 4;
    constexpr int SpOther = 5;
    constexpr int NSpecies = 6;
    constexpr int TofNone = 0;
    constexpr int TofGood = 1;
    constexpr int TofBad = 2;
    int species = SpOther;
    int origin = 1;
    float tofTrue = 0.f;
    if (jetTrack.has_mcParticle()) {
      const auto mcPart = jetTrack.template mcParticle_as<aod::JetParticles>();
      origin = mcPart.isPhysicalPrimary() ? 0 : 1;
      switch (std::abs(mcPart.pdgCode())) {
        case PDG_t::kPiPlus:
          species = SpPi;
          tofTrue = track.tofNSigmaPi();
          break;
        case PDG_t::kKPlus:
          species = SpKa;
          tofTrue = track.tofNSigmaKa();
          break;
        case PDG_t::kProton:
          species = SpPr;
          tofTrue = track.tofNSigmaPr();
          break;
        case PDG_t::kElectron:
          species = SpEl;
          tofTrue = track.tofNSigmaPi();
          break;
        case PDG_t::kMuonMinus:
          species = SpMu;
          tofTrue = track.tofNSigmaPi();
          break;
        default:
          break;
      }
    }
    int tofStatus = TofNone;
    if (track.hasTOF()) {
      tofStatus = (species != SpOther && std::abs(tofTrue) > pid2DOpts.pid2DMcTofBadNSigma) ? TofBad : TofGood;
    }
    const float index = species + NSpecies * (origin + 2 * tofStatus) + 0.5f;
    const float trkP = track.p();
    if (passSpeciesKinematics(track.pt(), track.eta(), ProtonMass)) {
      registry.fill(HIST("mcPid2DTpcOnlyPr"), trkP, track.tpcNSigmaPr(), index);
    }
    if (passSpeciesKinematics(track.pt(), track.eta(), PionMass)) {
      registry.fill(HIST("mcPid2DTpcOnlyPi"), trkP, track.tpcNSigmaPi(), index);
    }
    if (passSpeciesKinematics(track.pt(), track.eta(), KaonMass)) {
      registry.fill(HIST("mcPid2DTpcOnlyKa"), trkP, track.tpcNSigmaKa(), index);
    }
  }

  template <typename JetTrackT, typename TrackT, typename JetT>
  void fillPid2DMcTrack(JetTrackT const& jetTrack, TrackT const& track, std::vector<JetT> const& jets, float maxR2,
                        bool inPr, bool inPi, bool inKa)
  {
    constexpr int SpPi = 0;
    constexpr int SpKa = 1;
    constexpr int SpPr = 2;
    constexpr int SpEl = 3;
    constexpr int SpMu = 4;
    constexpr int SpOther = 5;
    constexpr int NSpecies = 6;
    constexpr int NSpeciesCone = 4;
    constexpr float JetPtMin2D = 20.f;
    constexpr float JetPtMax2D = 140.f;
    constexpr float RMax2D = 0.7f;
    constexpr float REdge1 = 0.2f;
    constexpr float REdge2 = 0.4f;

    const float trkP = track.p();
    const float tpcPr = track.tpcNSigmaPr();
    const float tofPr = track.tofNSigmaPr();

    int species = SpOther;
    int origin = 1;
    float tofTrue = 0.f;
    if (jetTrack.has_mcParticle()) {
      const auto mcPart = jetTrack.template mcParticle_as<aod::JetParticles>();
      origin = mcPart.isPhysicalPrimary() ? 0 : 1;
      switch (std::abs(mcPart.pdgCode())) {
        case PDG_t::kPiPlus:
          species = SpPi;
          tofTrue = track.tofNSigmaPi();
          break;
        case PDG_t::kKPlus:
          species = SpKa;
          tofTrue = track.tofNSigmaKa();
          break;
        case PDG_t::kProton:
          species = SpPr;
          tofTrue = tofPr;
          break;
        case PDG_t::kElectron: // beta ~ 1 at p > 2 GeV/c: the pion hypothesis is within 0.1 sigma
          species = SpEl;
          tofTrue = track.tofNSigmaPi();
          break;
        case PDG_t::kMuonMinus:
          species = SpMu;
          tofTrue = track.tofNSigmaPi();
          break;
        default:
          break;
      }
    }
    const int tofBad = (species != SpOther && std::abs(tofTrue) > pid2DOpts.pid2DMcTofBadNSigma) ? 1 : 0;
    const float incIndex = species + NSpecies * (origin + 2 * tofBad) + 0.5f;
    const float tpcPi = track.tpcNSigmaPi();
    const float tofPi = track.tofNSigmaPi();
    const float tpcKa = track.tpcNSigmaKa();
    const float tofKa = track.tofNSigmaKa();
    if (inPr) {
      registry.fill(HIST("mcPid2DIncPr"), trkP, tpcPr, tofPr, incIndex);
    }
    if (inPi) {
      registry.fill(HIST("mcPid2DIncPi"), trkP, tpcPi, tofPi, incIndex);
    }
    if (inKa) {
      registry.fill(HIST("mcPid2DIncKa"), trkP, tpcKa, tofKa, incIndex);
    }

    const int speciesCone = std::min(species, NSpeciesCone - 1);
    auto coneIndex = [&](float r) {
      const int ir = r < REdge1 ? 0 : (r < REdge2 ? 1 : 2);
      return speciesCone + NSpeciesCone * (tofBad + 2 * ir) + 0.5f;
    };

    const float trkEta = track.eta();
    const float trkPhi = track.phi();
    for (const auto& jet : jets) {
      if (jet.detPtCorr < JetPtMin2D || jet.detPtCorr >= JetPtMax2D) {
        continue;
      }
      const float dEta = trkEta - jet.detEta;
      if (std::abs(dEta) >= RMax2D) {
        continue;
      }
      const float dPhi = absDeltaPhi(trkPhi, jet.detPhi);
      const float distSq = dEta * dEta + dPhi * dPhi;
      if (distSq < maxR2) {
        const float r = std::sqrt(distSq);
        if (r < RMax2D) {
          const float index = coneIndex(r);
          if (inPr) {
            registry.fill(HIST("mcPid2DJetPr"), trkP, tpcPr, tofPr, index);
          }
          if (inPi) {
            registry.fill(HIST("mcPid2DJetPi"), trkP, tpcPi, tofPi, index);
          }
          if (inKa) {
            registry.fill(HIST("mcPid2DJetKa"), trkP, tpcKa, tofKa, index);
          }
        }
      }
      // Both perpendicular cones, as for pid2DPerpPr in data.
      for (const auto& conePhi : {jet.detPhiBg1, jet.detPhiBg2}) {
        const float dPhiBg = absDeltaPhi(trkPhi, conePhi);
        const float distBgSq = dEta * dEta + dPhiBg * dPhiBg;
        if (distBgSq < maxR2) {
          const float rBg = std::sqrt(distBgSq);
          if (rBg < RMax2D) {
            const float index = coneIndex(rBg);
            if (inPr) {
              registry.fill(HIST("mcPid2DPerpPr"), trkP, tpcPr, tofPr, index);
            }
            if (inPi) {
              registry.fill(HIST("mcPid2DPerpPi"), trkP, tpcPi, tofPi, index);
            }
            if (inKa) {
              registry.fill(HIST("mcPid2DPerpKa"), trkP, tpcKa, tofKa, index);
            }
          }
        }
      }
    }
  }

  void processEfficiencyAndPurity(
    soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator const& collision,
    JetTracksWithPIMC const& jetTracks,
    FullTrackInfoMC const&,
    soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdJets,
    aod::ChargedMCParticleLevelJets const& mcpJets,
    aod::JetParticles const& mcParticles)
  {
    constexpr float CutFlowCalled = 0.5f;
    constexpr float CutFlowEventSelection = 1.5f;
    constexpr float CutFlowHasMcCollision = 2.5f;
    constexpr float CutFlowCentrality = 3.5f;

    constexpr float SelDetailCalled = 0.5f;
    constexpr float SelDetailNotMbGap = 1.5f;
    constexpr float SelDetailRct = 2.5f;
    constexpr float SelDetailEventBits = 3.5f;

    registry.fill(HIST("mcCutFlow"), CutFlowCalled);
    registry.fill(HIST("mcSelDetail"), SelDetailCalled);
    registry.fill(HIST("mcSubGeneratorId"), collision.getSubGeneratorId());
    if (collision.getSubGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("mcSelDetail"), SelDetailNotMbGap);
    }
    if (jetderiveddatautilities::selectCollision(collision, emptyEventSelectionBits, false, true)) {
      registry.fill(HIST("mcSelDetail"), SelDetailRct);
    }
    if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false, false)) {
      registry.fill(HIST("mcSelDetail"), SelDetailEventBits);
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMbGapEvents, applyRctSelection)) {
      return;
    }
    registry.fill(HIST("mcCutFlow"), CutFlowEventSelection);
    if (!collision.has_mcCollision()) {
      return;
    }
    registry.fill(HIST("mcCutFlow"), CutFlowHasMcCollision);

    const auto reducedMcCollisionId = collision.mcCollisionId();
    (void)mcpJets;

    const float centrality = collision.centFT0M();
    registry.fill(HIST("mcCentralityRaw"), centrality);
    if (!passCentrality(centrality)) {
      return;
    }
    registry.fill(HIST("mcCutFlow"), CutFlowCentrality);
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

      if (pid2DOpts.fillPid2DMc && pid2DOpts.fillPid2DTpcOnly && passesNominal && !validJets.empty() &&
          centrality >= pid2DOpts.pid2DCentMin && centrality < pid2DOpts.pid2DCentMax &&
          track.p() >= pid2DOpts.pid2DPMin && track.p() < pid2DOpts.pid2DTpcOnlyPMax) {
        fillPid2DMcTpcOnlyTrack(jetTrack, track);
      }
      if (pid2DOpts.fillPid2DMc && passesNominal && track.hasTOF() && !validJets.empty() &&
          centrality >= pid2DOpts.pid2DCentMin && centrality < pid2DOpts.pid2DCentMax) {
        const float trkPMc = track.p();
        const bool inPr = trkPMc >= pid2DOpts.pid2DPMin && trkPMc < pid2DOpts.pid2DPMax && passSpeciesKinematics(track.pt(), track.eta(), ProtonMass);
        const bool inPi = pid2DOpts.fillPid2DPi && trkPMc >= pid2DOpts.pid2DPiPMin && trkPMc < pid2DOpts.pid2DPiPMax && passSpeciesKinematics(track.pt(), track.eta(), PionMass);
        const bool inKa = pid2DOpts.fillPid2DKa && trkPMc >= pid2DOpts.pid2DKaPMin && trkPMc < pid2DOpts.pid2DKaPMax && passSpeciesKinematics(track.pt(), track.eta(), KaonMass);
        if (inPr || inPi || inKa) {
          fillPid2DMcTrack(jetTrack, track, validJets, maxR2, inPr, inPi, inKa);
        }
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
      const bool passKineRecoPi = passSpeciesKinematics(recoPt, track.eta(), PionMass);
      const bool passKineRecoPr = passSpeciesKinematics(recoPt, track.eta(), ProtonMass);
      const bool passTofPi = hasTof && std::abs(track.tofNSigmaPi()) < nSigmaTofCut;
      const bool passTofPr = hasTof && std::abs(track.tofNSigmaPr()) < nSigmaTofCut;
      const bool isRecoPi = passKineRecoPi && passTpcPi && passTofPi;
      const bool isRecoPr = passKineRecoPr && passTpcPr && passTofPr;

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

      constexpr float MaxSaneNSigma = 10.f; // nSigma_TOF is filled with a sentinel value when the event time is missing
      constexpr float NoTofBin = 0.5f;
      constexpr float FiniteNSigmaBin = 1.5f;
      constexpr float UnphysicalNSigmaBin = 2.5f;

      if (isTruePi) {
        registry.fill(HIST("effNumTpcPiInc"), mcPt, centrality);
        registry.fill(HIST("mcTofNSigmaStatusPi"), !hasTof ? NoTofBin : (std::abs(track.tofNSigmaPi()) < MaxSaneNSigma ? FiniteNSigmaBin : UnphysicalNSigmaBin));
        if (hasTof) {
          registry.fill(HIST("effNumTofPiInc"), mcPt, centrality);
          registry.fill(HIST("mcTpcNSigmaPi"), track.p(), track.tpcNSigmaPi(), centrality);
          registry.fill(HIST("mcTofNSigmaPi"), mcPt, track.tofNSigmaPi(), centrality);
          if (fillTofMismatchQA) {
            registry.fill(HIST("mcTofNSigmaWidePi"), track.p(), track.tofNSigmaPi(), centrality);
          }
          if (passKineRecoPi) {
            registry.fill(HIST("effNumKinePiInc"), mcPt, centrality);
            if (passTpcPi) {
              registry.fill(HIST("effNumTpcWinPiInc"), mcPt, centrality);
            }
            if (passTofPi) {
              registry.fill(HIST("effNumTofWinPiInc"), mcPt, centrality);
            }
          }
        }
        if (isRecoPi) {
          registry.fill(HIST("effNumPidPiInc"), mcPt, centrality);
        }
      } else if (isTrueKa) {
        registry.fill(HIST("effNumTpcKaInc"), mcPt, centrality);
        if (hasTof) {
          registry.fill(HIST("effNumTofKaInc"), mcPt, centrality);
          if (fillTofMismatchQA) {
            registry.fill(HIST("mcTofNSigmaWideKa"), track.p(), track.tofNSigmaKa(), centrality);
          }
        }
      } else {
        registry.fill(HIST("effNumTpcPrInc"), mcPt, centrality);
        registry.fill(HIST("mcTofNSigmaStatusPr"), !hasTof ? NoTofBin : (std::abs(track.tofNSigmaPr()) < MaxSaneNSigma ? FiniteNSigmaBin : UnphysicalNSigmaBin));
        if (hasTof) {
          registry.fill(HIST("effNumTofPrInc"), mcPt, centrality);
          registry.fill(HIST("mcTpcNSigmaPr"), track.p(), track.tpcNSigmaPr(), centrality);
          registry.fill(HIST("mcTofNSigmaPr"), mcPt, track.tofNSigmaPr(), centrality);
          if (fillTofMismatchQA) {
            registry.fill(HIST("mcTofNSigmaWidePr"), track.p(), track.tofNSigmaPr(), centrality);
          }
          if (passKineRecoPr) {
            registry.fill(HIST("effNumKinePrInc"), mcPt, centrality);
            if (passTpcPr) {
              registry.fill(HIST("effNumTpcWinPrInc"), mcPt, centrality);
            }
            if (passTofPr) {
              registry.fill(HIST("effNumTofWinPrInc"), mcPt, centrality);
            }
          }
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
