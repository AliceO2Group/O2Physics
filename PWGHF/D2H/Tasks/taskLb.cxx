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

/// \file taskLb.cxx
/// \brief Λb0 analysis task
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef
/// \author Martin Voelkl <martin.andreas.volkl@cern.ch>, University of Birmingham

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_decay::hf_cand_beauty;

/// Λb0 analysis task
struct HfTaskLb {
  Configurable<int> selectionFlagLb{"selectionFlagLb", 0, "Selection Flag for Lb"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> lengthDCAParameter{"lengthDCAParameter", 0.02, "decay length for DCA"};
  Configurable<float> minLikelihoodRatio{"minLikelihoodRatio", 10., "min. likelihood ratio for combined DCAs"};
  Configurable<float> minLikelihoodRatioLc{"minLikelihoodRatioLc", 10., "min. likelihood ratio for Lc cross check"};
  Configurable<float> mDiffKStar892Max{"mDiffKStar892Max", 0.0473, "Accepted range around KStar mass peak"};
  Configurable<float> mDiffDelta1232Max{"mDiffDelta1232Max", 0.117, "Accepted range around Delta mass peak"};
  Configurable<float> mDiffLambda1520Max{"mDiffLambda1520Max", 0.016 * 2., "Accepted range around Lambda 1520 mass peak"};
  Configurable<float> mDiffLcMax{"mDiffLcMax", 0.1, "Accepted range around LambdaC mass peak for filling two body mass histograms"};
  Configurable<float> maximumImpactParameterForLambdaCCrossChecks{"maximumImpactParameterForLambdaCCrossChecks", 0.2, "maximum d0 for LambdaC checks"};
  Configurable<float> resoCorrectionFactor{"resoCorrectionFactor", 1.1, "Resolution correction compared to reconstruction estimate"};
  Configurable<float> largeLifetimeBG{"largeLifetimeBG", 0.01, "fraction of strange contribution within 2mm"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lb_to_lc_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;
  Service<o2::framework::O2DatabasePDG> pdg;

  using TracksWExt = soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  using TracksWExtMc = soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, McTrackLabels>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lb::isSelLbToLcPi >= selectionFlagLb);

  PresliceUnsorted<aod::TracksWMc> mcPartID = aod::mctracklabel::mcParticleId;

  bool passesImpactParameterResolution(float pT, float d0Resolution)
  {
    float const expectedResolution(0.001 + 0.0052 * std::exp(-0.655 * pT));
    return (d0Resolution <= expectedResolution * 1.5);
  } // Compares to pT dependent cut on impact parameter resolution

  float logLikelihoodRatioSingleTrackDCA(float dca, float reso, float lengthParameter)
  {
    reso *= resoCorrectionFactor; // In case real resolution is worse
    float const numerator = 1. / lengthParameter * std::exp(-dca / lengthParameter);
    float const denominator = (1. - largeLifetimeBG) * TMath::Gaus(dca, 0., reso, true) + largeLifetimeBG / 0.2; // flat distribution to 2 mm
    return std::log(numerator / denominator);
  } // Creates the single track log likelihood assuming an exonential law for the secondaries

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "Lb candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", "Lb candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", "Lb candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hIPs", "Impact parameters;p_{T} (GeV/#it{c});d_{0} (cm)", {HistType::kTH2F, {{20, 0., 20.}, {400, 0., 0.2}}}},
     {"hIPsAfterCut", "Impact parameters;p_{T} (GeV/#it{c});d_{0} (cm)", {HistType::kTH2F, {{20, 0., 20.}, {400, 0., 0.2}}}},
     {"hIPResolution", "Impact parameter resolution;p_{T} (GeV/#it{c});#sigma_{d_{0}} (cm)", {HistType::kTH2F, {{20, 0., 10.}, {400, 0., 0.02}}}},
     {"hPtlogLikelihood", "log Likelihood;p_{T} (GeV/#it{c});log L", {HistType::kTH2F, {{20, 0., 20.}, {400, -10., 70.}}}},
     {"hPtinvMassKStar", "K^{*}(892) invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 0.5, 1.5}}}},
     {"hPtinvMassDelta", "#Delta(1232) invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.0, 2.0}}}},
     {"hPtinvMassLambda1520", "#Lambda(1520) invariant maas;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.0, 2.0}}}},
     {"hPtinvMassLcKStar", "#Lambda_{c} invariant mass from K^{*};p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"hPtinvMassLcDelta", "#Lambda_{c} invariant mass from #Delta;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"hPtinvMassLcLambda1520", "#Lambda_{c} invariant mass from #Lambda;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"hPtinvMassLc", "#Lambda_{c} invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"hPtinvMassLcReso", "#Lambda_{c} from resonances invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"hPtinvMassLb", "#Lambda_{b} invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 3.5, 7.5}}}},
     {"hZVertex", "z Vertex;z_{vtx};counts", {HistType::kTH1F, {{100, -20., 20.}}}},
     {"MC/hPtRecSig", "Lb candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"MC/hPtRecBg", "Lb candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"MC/hPtGenSig", "Lb candidates (matched);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 10.}}}},
     {"MC/hPtGen", "MC particles (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}}}};

  void init(InitContext&)
  {
    registry.add("hMass", "#Lambda_{b}^{0} candidates;inv. mass #Lambda_{c}^{#plus}#pi^{#minus} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{500, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "#Lambda_{b}^{0} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXY", "#Lambda_{b}^{0} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "#Lambda_{b}^{0} candidates;prong 0 (#Lambda_{c}^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "#Lambda_{b}^{0} candidates;prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidity", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hIPProd", "#Lambda_{b}^{0} candidates;#Lambda_{b}^{0} candidate impact parameter product;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hInvMassLc", "#Lambda_{b}^{0} candidates;prong0, #Lambda_{c}^{+} inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0, 5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    // Now add MC histograms
    registry.add("MC/hEtaGen", "MC particles (matched);#Lambda_{b}^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hYGen", "MC particles (matched);#Lambda_{b}^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hPtProng0Gen", "MC particles (matched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hPtProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hYProng0Gen", "MC particles (matched);prong 0 (#Lambda_{c}^{+}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hYProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hEtaProng0Gen", "MC particles (matched);prong 0 (#Lambda_{b}^{0}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hEtaProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hCPARecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hCPARecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hCPAxyRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hCPAxyRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hCPALcRecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hCPALcRecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hEtaRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hEtaRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hRapidityRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hRapidityRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate #it{#y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/hPtProng0RecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hPtProng1RecSig", "#Lambda_{b}^{0} candidates (matched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hPtProng0RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hPtProng1RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hMassRecSig", "#Lambda_{b}^{0} candidates (matched);inv. mass #Lambda_{c}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hMassRecBg", "#Lambda_{b}^{0} candidates (unmatched);inv. mass #Lambda_{c}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hd0Prong0RecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hd0Prong1RecSig", "#Lambda_{b}^{0} candidates (matched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hd0Prong0RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hd0Prong1RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthXYRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthXYRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthLcRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthNormRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hDecLengthNormRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hImpParProdLbRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hImpParProdLbRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/hChi2PCARecSig", "#Lambda_{b}^{0} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hChi2PCARecBg", "#Lambda_{b}^{0} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hThetaStarRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/hThetaStarRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void processData(aod::Collisions::iterator const& collision,
                   soa::Filtered<soa::Join<aod::HfCandLb, aod::HfSelLbToLcPi>> const& candidates,
                   soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidatesLc,
                   TracksWExt const&)
  {
    float const massKStar892 = 0.892;
    float const massDelta1232 = 1.232;
    std::array<float, 3> dca = {0.f, 0.f, 0.f};
    std::array<float, 3> dcaResolution = {0.f, 0.f, 0.f};

    for (const auto& candidateLc : candidatesLc) {
      if ((candidateLc.isSelLcToPKPi() == 0) && (candidateLc.isSelLcToPiKP() == 0)) {
        continue;
      }
      auto track0 = candidateLc.prong0_as<TracksWExt>();
      auto track1 = candidateLc.prong1_as<TracksWExt>();
      auto track2 = candidateLc.prong2_as<TracksWExt>();
      registry.get<TH2>(HIST("hIPs"))->Fill(candidateLc.pt(), candidateLc.impactParameter0());
      registry.get<TH2>(HIST("hIPs"))->Fill(candidateLc.pt(), candidateLc.impactParameter1());
      registry.get<TH2>(HIST("hIPs"))->Fill(candidateLc.pt(), candidateLc.impactParameter2());
      float const reso0 = candidateLc.errorImpactParameter0(); // 0.0023166 *pow(track0.pt(), -0.788);
      float const reso1 = candidateLc.errorImpactParameter1();
      float const reso2 = candidateLc.errorImpactParameter2();
      registry.get<TH2>(HIST("hIPResolution"))->Fill(track0.pt(), reso0);
      registry.get<TH2>(HIST("hIPResolution"))->Fill(track1.pt(), reso1);
      registry.get<TH2>(HIST("hIPResolution"))->Fill(track2.pt(), reso2);
      if (!passesImpactParameterResolution(track0.pt(), reso0)) {
        continue;
      }
      if (!passesImpactParameterResolution(track1.pt(), reso1)) {
        continue;
      }
      if (!passesImpactParameterResolution(track2.pt(), reso2)) {
        continue;
      }

      dca = {
        candidateLc.impactParameter0(),
        candidateLc.impactParameter1(),
        candidateLc.impactParameter2()};

      bool const exceedsMaxDca = std::any_of(dca.begin(), dca.end(), [&](float val) {
        return val > maximumImpactParameterForLambdaCCrossChecks;
      });

      if (exceedsMaxDca) {
        continue;
      }
      dcaResolution = {reso0, reso1, reso2};

      float likelihoodRatio = 0.0f;
      for (size_t i = 0; i < dca.size(); ++i) {
        likelihoodRatio += logLikelihoodRatioSingleTrackDCA(dca[i], dcaResolution[i], lengthDCAParameter);
      }

      registry.get<TH2>(HIST("hPtlogLikelihood"))->Fill(candidateLc.pt(), likelihoodRatio);
      if (likelihoodRatio < minLikelihoodRatioLc) {
        continue;
      }
      registry.get<TH2>(HIST("hIPsAfterCut"))->Fill(candidateLc.pt(), candidateLc.impactParameter0());
      registry.get<TH2>(HIST("hIPsAfterCut"))->Fill(candidateLc.pt(), candidateLc.impactParameter1());
      registry.get<TH2>(HIST("hIPsAfterCut"))->Fill(candidateLc.pt(), candidateLc.impactParameter2());
      if (candidateLc.isSelLcToPKPi() != 0) {
        registry.get<TH2>(HIST("hPtinvMassLc"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPKPi(candidateLc));
        float const mRecoKstar = RecoDecay::m(std::array{track1.pVector(), track2.pVector()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
        float const mRecoDelta1232 = RecoDecay::m(std::array{track0.pVector(), track2.pVector()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
        float const mRecoLambda1520 = RecoDecay::m(std::array{track0.pVector(), track1.pVector()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus});
        float const mDiffKStar892 = std::abs(mRecoKstar - massKStar892);
        float const mDiffDelta1232 = std::abs(mRecoDelta1232 - massDelta1232);
        float const mDiffLambda1520 = std::abs(mRecoLambda1520 - o2::constants::physics::MassLambda1520);
        if (mDiffKStar892 < mDiffKStar892Max || mDiffDelta1232 < mDiffDelta1232Max || mDiffLambda1520 < mDiffLambda1520Max) {
          registry.get<TH2>(HIST("hPtinvMassLcReso"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPKPi(candidateLc));
        }
        if (mDiffKStar892 < mDiffKStar892Max) {
          registry.get<TH2>(HIST("hPtinvMassLcKStar"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPKPi(candidateLc));
        }
        if (mDiffDelta1232 < mDiffDelta1232Max) {
          registry.get<TH2>(HIST("hPtinvMassLcDelta"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPKPi(candidateLc));
        }
        if (mDiffLambda1520 < mDiffLambda1520Max) {
          registry.get<TH2>(HIST("hPtinvMassLcLambda1520"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPKPi(candidateLc));
        }

        if (std::abs(hfHelper.invMassLcToPKPi(candidateLc) - o2::constants::physics::MassLambdaCPlus) < mDiffLcMax) {
          registry.get<TH2>(HIST("hPtinvMassKStar"))->Fill(candidateLc.pt(), mRecoKstar);
          registry.get<TH2>(HIST("hPtinvMassDelta"))->Fill(candidateLc.pt(), mRecoDelta1232);
          registry.get<TH2>(HIST("hPtinvMassLambda1520"))->Fill(candidateLc.pt(), mRecoLambda1520);
        }
      }
      if (candidateLc.isSelLcToPiKP() != 0) {
        registry.get<TH2>(HIST("hPtinvMassLc"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPiKP(candidateLc));
        float const mRecoKstar = RecoDecay::m(std::array{track1.pVector(), track0.pVector()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
        float const mRecoDelta1232 = RecoDecay::m(std::array{track2.pVector(), track0.pVector()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
        float const mRecoLambda1520 = RecoDecay::m(std::array{track2.pVector(), track1.pVector()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus});
        float const mDiffKStar892 = std::abs(mRecoKstar - massKStar892);
        float const mDiffDelta1232 = std::abs(mRecoDelta1232 - massDelta1232);
        float const mDiffLambda1520 = std::abs(mRecoLambda1520 - o2::constants::physics::MassLambda1520);
        if (mDiffKStar892 < mDiffKStar892Max || mDiffDelta1232 < mDiffDelta1232Max || mDiffLambda1520 < mDiffLambda1520Max) {
          registry.get<TH2>(HIST("hPtinvMassLcReso"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPiKP(candidateLc));
        }
        if (mDiffKStar892 < mDiffKStar892Max) {
          registry.get<TH2>(HIST("hPtinvMassLcKStar"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPiKP(candidateLc));
        }
        if (mDiffDelta1232 < mDiffDelta1232Max) {
          registry.get<TH2>(HIST("hPtinvMassLcDelta"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPiKP(candidateLc));
        }
        if (mDiffLambda1520 < mDiffLambda1520Max) {
          registry.get<TH2>(HIST("hPtinvMassLcLambda1520"))->Fill(candidateLc.pt(), hfHelper.invMassLcToPiKP(candidateLc));
        }

        if (std::abs(hfHelper.invMassLcToPiKP(candidateLc) - o2::constants::physics::MassLambdaCPlus) < mDiffLcMax) {
          registry.get<TH2>(HIST("hPtinvMassKStar"))->Fill(candidateLc.pt(), mRecoKstar);
          registry.get<TH2>(HIST("hPtinvMassDelta"))->Fill(candidateLc.pt(), mRecoDelta1232);
          registry.get<TH2>(HIST("hPtinvMassLambda1520"))->Fill(candidateLc.pt(), mRecoLambda1520);
        }
      }
    } // Lambda_c candidates loop for cross checks

    for (const auto& candidate : candidates) {

      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      registry.get<TH1>(HIST("hZVertex"))->Fill(collision.posZ());

      auto candLc = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
      dca = {
        candLc.impactParameter0(),
        candLc.impactParameter1(),
        candLc.impactParameter2()};

      dcaResolution = {
        candLc.errorImpactParameter0(),
        candLc.errorImpactParameter1(),
        candLc.errorImpactParameter2()};

      float likelihoodRatio = 0.0f;
      for (size_t i = 0; i < dca.size(); ++i) {
        likelihoodRatio += logLikelihoodRatioSingleTrackDCA(dca[i], dcaResolution[i], lengthDCAParameter);
      }

      if (likelihoodRatio < minLikelihoodRatio) {
        continue; // Larger likelihood means more likely to be signal
      }
      float const lbMass = hfHelper.invMassLbToLcPi(candidate);
      registry.get<TH2>(HIST("hPtinvMassLb"))->Fill(candidate.pt(), lbMass);

      registry.fill(HIST("hMass"), hfHelper.invMassLbToLcPi(candidate), candidate.pt());
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hIPProd"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hRapidity"), hfHelper.yLb(candidate), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
      registry.fill(HIST("hInvMassLc"), lbMass, candidate.pt());
    } // candidate loop
  }
  PROCESS_SWITCH(HfTaskLb, processData, "Process Data", true);

  void processMc(soa::Filtered<soa::Join<aod::HfCandLb, aod::HfSelLbToLcPi, aod::HfCandLbMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandLbMcGen> const& mcParticles,
                 TracksWExtMc const&,
                 aod::HfCand3Prong const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {

      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      auto candLc = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec>>();
      auto flagMcMatchRecLb = std::abs(candidate.flagMcMatchRec());

      if (flagMcMatchRecLb == DecayChannelMain::LbToLcPi) {

        auto indexMother = RecoDecay::getMother(mcParticles, candidate.prong1_as<TracksWExtMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandLbMcGen>>(), o2::constants::physics::Pdg::kLambdaB0, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        registry.fill(HIST("MC/hPtGenSig"), particleMother.pt());
        registry.fill(HIST("MC/hPtRecSig"), candidate.pt());
        registry.fill(HIST("MC/hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/hCPAxyRecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/hEtaRecSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("MC/hRapidityRecSig"), hfHelper.yLb(candidate), candidate.pt());
        registry.fill(HIST("MC/hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/hDecLengthXYRecSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("MC/hMassRecSig"), hfHelper.invMassLbToLcPi(candidate), candidate.pt());
        registry.fill(HIST("MC/hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("MC/hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("MC/hImpParProdLbRecSig"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("MC/hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("MC/hCPALcRecSig"), candLc.cpa(), candidate.pt());
        registry.fill(HIST("MC/hDecLengthLcRecSig"), candLc.decayLength(), candidate.pt());
        registry.fill(HIST("MC/hChi2PCARecSig"), candidate.chi2PCA(), candidate.pt());
        // registry.fill(HIST("MC/hThetaStarRecSig"), candidate.cosThetaStar(), candidate.pt());
      } else {
        registry.fill(HIST("MC/hPtRecBg"), candidate.pt());
        registry.fill(HIST("MC/hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/hCPAxyRecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("MC/hRapidityRecBg"), hfHelper.yLb(candidate), candidate.pt());
        registry.fill(HIST("MC/hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/hDecLengthXYRecBg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("MC/hMassRecBg"), hfHelper.invMassLbToLcPi(candidate), candidate.pt());
        registry.fill(HIST("MC/hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("MC/hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("MC/hImpParProdLbRecBg"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("MC/hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("MC/hCPALcRecBg"), candLc.cpa(), candidate.pt());
        registry.fill(HIST("MC/hDecLengthLcRecBg"), candLc.decayLength(), candidate.pt());
        registry.fill(HIST("MC/hChi2PCARecBg"), candidate.chi2PCA(), candidate.pt());
        // registry.fill(HIST("MC/hThetaStarRecBg"), candidate.cosThetaStar(), candidate.pt());
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == DecayChannelMain::LbToLcPi) {

        auto yParticle = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaB0);
        if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
          continue;
        }

        float ptProngs[2], yProngs[2], etaProngs[2];
        int counter = 0;
        for (const auto& daught : particle.daughters_as<soa::Join<aod::McParticles, aod::HfCandLbMcGen>>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(daught.pVector(), pdg->Mass(daught.pdgCode()));
          counter++;
        }

        registry.fill(HIST("MC/hPtProng0Gen"), ptProngs[0], particle.pt());
        registry.fill(HIST("MC/hPtProng1Gen"), ptProngs[1], particle.pt());
        registry.fill(HIST("MC/hYProng0Gen"), yProngs[0], particle.pt());
        registry.fill(HIST("MC/hYProng1Gen"), yProngs[1], particle.pt());
        registry.fill(HIST("MC/hEtaProng0Gen"), etaProngs[0], particle.pt());
        registry.fill(HIST("MC/hEtaProng1Gen"), etaProngs[1], particle.pt());

        //  if (yCandMax >= 0. && (std::abs(yProngs[0]) > yCandMax || std::abs(yProngs[1]) > yCandMax))
        //    continue;

        registry.fill(HIST("MC/hPtGen"), particle.pt());
        registry.fill(HIST("MC/hYGen"), yParticle, particle.pt());
        registry.fill(HIST("MC/hEtaGen"), particle.eta(), particle.pt());
      }
    } // gen
  }
  PROCESS_SWITCH(HfTaskLb, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskLb>(cfgc)};
}
