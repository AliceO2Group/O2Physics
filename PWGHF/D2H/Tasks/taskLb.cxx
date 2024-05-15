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

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/Centrality.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Λb0 analysis task
struct HfTaskLb {
  Configurable<int> selectionFlagLb{"selectionFlagLb", 0, "Selection Flag for Lb"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<double> DCALengthParameter{"DCALengthParameter", 0.02, "decay length for DCA"};
  Configurable<double> minLikelihoodRatio{"minLikelihoodRatio", 10., "min. likelihood ratio for combined DCAs"};
  Configurable<double> minLikelihoodRatioLc{"minLikelihoodRatioLc", 10., "min. likelihood ratio for Lc cross check"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lb_to_lc_pi::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lb::isSelLbToLcPi >= selectionFlagLb);
  float vtxCut = 10.;
  Filter posZfilter = nabs(aod::collision::posZ) < vtxCut;

  using TracksWExt = soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;

  double InvariantMass(auto track0, auto track1, double m0, double m1)
  {
    double px0 = track0.px();
    double px1 = track1.px();
    double pxSum = px0 + px1;
    double py0 = track0.py();
    double py1 = track1.py();
    double pySum = py0 + py1;
    double pz0 = track0.pz();
    double pz1 = track1.pz();
    double pzSum = pz0 + pz1;
    double pSq0 = px0 * px0 + py0 * py0 + pz0 * pz0;
    double pSq1 = px1 * px1 + py1 * py1 + pz1 * pz1;
    double pSumSq = pxSum * pxSum + pySum * pySum + pzSum * pzSum;
    double eSq0 = pSq0 + m0 * m0;
    double eSq1 = pSq1 + m1 * m1;
    double eTot = sqrt(eSq0) + sqrt(eSq1);
    double eSq = eTot * eTot;
    return sqrt(eSq - pSumSq);
  }

  bool PassesImpactParameterResolution(double pT, double d0Resolution)
  {
    double expectedResolution(0.001 + 0.0052 * exp(-0.655 * pT));
    if (d0Resolution > expectedResolution * 1.5)
      return false;
    else
      return true;
  } // Compares to pT dependent cut on impact parameter resolution

  double LogLikelihoodRatioSingleTrackDCA(double DCA, double reso, double lengthParameter)
  {
    reso *= 1.1;                         // In case real resolution is worse
    double largeLifetimeFraction = 0.01; // contribution from strange decays or similar
    double numerator = 1 / lengthParameter * exp(-DCA / lengthParameter);
    double denominator = (1. - largeLifetimeFraction) * TMath::Gaus(DCA, 0., reso) + largeLifetimeFraction / 0.2; // flat distribution to 2 mm
    return log(numerator / denominator);
  } // Creates the single track log likelihood assuming an exonential law for the secondaries

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "Lb candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"hPtProng1", "Lb candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hPtCand", "Lb candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 50.}}}},
     {"IPs", "Impact parameters;p_{T} (GeV/#it{c});d_{0} (cm)", {HistType::kTH2F, {{20, 0., 20.}, {400, 0., 0.2}}}},
     {"IPsAfterCut", "Impact parameters;p_{T} (GeV/#it{c});d_{0} (cm)", {HistType::kTH2F, {{20, 0., 20.}, {400, 0., 0.2}}}},
     {"IPResolution", "Impact parameter resolution;p_{T} (GeV/#it{c});#sigma_{d_{0}} (cm)", {HistType::kTH2F, {{20, 0., 10.}, {400, 0., 0.02}}}},
     {"pTlogLikelihood", "log Likelihood;p_{T} (GeV/#it{c});log L", {HistType::kTH2F, {{20, 0., 20.}, {400, -10., 70.}}}},
     {"pTinvMassKStar", "K^{*}(892) invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 0.5, 1.5}}}},
     {"pTinvMassDelta", "#Delta(1232) invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.0, 2.0}}}},
     {"pTinvMassLambda1520", "#Lambda(1520) invariant maas;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.0, 2.0}}}},
     {"pTinvMassLcKStar", "#Lambda_{c} invariant mass from K^{*};p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"pTinvMassLcDelta", "#Lambda_{c} invariant mass from #Delta;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"pTinvMassLcLambda1520", "#Lambda_{c} invariant mass from #Lambda;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"pTinvMassLc", "#Lambda_{c} invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"pTinvMassLcReso", "#Lambda_{c} from resonances invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 1.5, 3.5}}}},
     {"pTinvMassLb", "#Lambda_{b} invariant mass;p_{T} (GeV/#it{c});m_{inv}", {HistType::kTH2F, {{20, 0., 20.}, {400, 3.5, 7.5}}}}}};

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
  }

  void process(soa::Filtered<aod::Collisions>::iterator const& collision,
               soa::Filtered<soa::Join<aod::HfCandLb, aod::HfSelLbToLcPi>> const& candidates,
               soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& lcCandidates,
               TracksWExt const&)
  {

    for (const auto& lccandidate : lcCandidates) {
      if (!lccandidate.isSelLcToPKPi() && !lccandidate.isSelLcToPiKP())
        continue;
      if (abs(lccandidate.eta() > 0.6))
        continue;
      auto track0 = lccandidate.prong0_as<TracksWExt>();
      auto track1 = lccandidate.prong1_as<TracksWExt>();
      auto track2 = lccandidate.prong2_as<TracksWExt>();
      registry.get<TH2>(HIST("IPs"))->Fill(lccandidate.pt(), lccandidate.impactParameter0());
      registry.get<TH2>(HIST("IPs"))->Fill(lccandidate.pt(), lccandidate.impactParameter1());
      registry.get<TH2>(HIST("IPs"))->Fill(lccandidate.pt(), lccandidate.impactParameter2());
      double reso0 = lccandidate.errorImpactParameter0(); // 0.0023166 *pow(track0.pt(), -0.788);
      double reso1 = lccandidate.errorImpactParameter1();
      double reso2 = lccandidate.errorImpactParameter2();
      if (!PassesImpactParameterResolution(track0.pt(), reso0))
        continue;
      if (!PassesImpactParameterResolution(track1.pt(), reso1))
        continue;
      if (!PassesImpactParameterResolution(track2.pt(), reso2))
        continue;
      registry.get<TH2>(HIST("IPResolution"))->Fill(track0.pt(), reso0);
      registry.get<TH2>(HIST("IPResolution"))->Fill(track1.pt(), reso1);
      registry.get<TH2>(HIST("IPResolution"))->Fill(track2.pt(), reso2);
      double DCA0 = lccandidate.impactParameter0();
      double DCA1 = lccandidate.impactParameter1();
      double DCA2 = lccandidate.impactParameter2();
      if (DCA0 > 0.2 || DCA1 > 0.2 || DCA2 > 0.2)
        continue;
      double likelihoodRatio = LogLikelihoodRatioSingleTrackDCA(DCA0, reso0, DCALengthParameter) + LogLikelihoodRatioSingleTrackDCA(DCA1, reso1, DCALengthParameter) + LogLikelihoodRatioSingleTrackDCA(DCA2, reso2, DCALengthParameter);
      registry.get<TH2>(HIST("pTlogLikelihood"))->Fill(lccandidate.pt(), likelihoodRatio);
      // if(lccandidate.impactParameter0()<0.003 || lccandidate.impactParameter1()<0.003 || lccandidate.impactParameter2()<0.003) continue;
      if (likelihoodRatio < minLikelihoodRatioLc)
        continue; // Likelihood ratio; 3.9 corresponds to 50; 3 corresponds to 20; 1.6 corresponds to 5; 5 corresponds to 100
      registry.get<TH2>(HIST("IPsAfterCut"))->Fill(lccandidate.pt(), lccandidate.impactParameter0());
      registry.get<TH2>(HIST("IPsAfterCut"))->Fill(lccandidate.pt(), lccandidate.impactParameter1());
      registry.get<TH2>(HIST("IPsAfterCut"))->Fill(lccandidate.pt(), lccandidate.impactParameter2());
      if (lccandidate.isSelLcToPKPi()) {
        registry.get<TH2>(HIST("pTinvMassLc"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPKPi(lccandidate));
        double mDiffKStar892 = abs(InvariantMass(track1, track2, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus) - 0.892);
        double mDiffDelta1232 = abs(InvariantMass(track0, track2, o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus) - 1.232);
        double mDiffLambda1520 = abs(InvariantMass(track0, track1, o2::constants::physics::MassProton, o2::constants::physics::MassKPlus) - 1.520);
        if (mDiffKStar892 < 0.07 || mDiffDelta1232 < 0.117 || mDiffLambda1520 < 0.05)
          registry.get<TH2>(HIST("pTinvMassLcReso"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPKPi(lccandidate));
        if (mDiffKStar892 < 0.046 * 2.)
          registry.get<TH2>(HIST("pTinvMassLcKStar"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPKPi(lccandidate));
        if (mDiffDelta1232 < 0.117)
          registry.get<TH2>(HIST("pTinvMassLcDelta"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPKPi(lccandidate));
        if (mDiffLambda1520 < 0.032)
          registry.get<TH2>(HIST("pTinvMassLcLambda1520"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPKPi(lccandidate));

        if (abs(hfHelper.invMassLcToPKPi(lccandidate) - o2::constants::physics::MassLambdaCPlus) < 0.05) {
          registry.get<TH2>(HIST("pTinvMassKStar"))->Fill(lccandidate.pt(), InvariantMass(track1, track2, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus));
          registry.get<TH2>(HIST("pTinvMassDelta"))->Fill(lccandidate.pt(), InvariantMass(track0, track2, o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus));
          registry.get<TH2>(HIST("pTinvMassLambda1520"))->Fill(lccandidate.pt(), InvariantMass(track0, track1, o2::constants::physics::MassProton, o2::constants::physics::MassKPlus));
        }
      }
      if (lccandidate.isSelLcToPiKP()) {
        registry.get<TH2>(HIST("pTinvMassLc"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPiKP(lccandidate));
        double mDiffKStar892 = abs(InvariantMass(track1, track0, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus) - 0.892);
        double mDiffDelta1232 = abs(InvariantMass(track2, track0, o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus) - 1.232);
        double mDiffLambda1520 = abs(InvariantMass(track2, track1, o2::constants::physics::MassProton, o2::constants::physics::MassKPlus) - 1.520);
        if (mDiffKStar892 < 0.07 || mDiffDelta1232 < 0.117 || mDiffLambda1520 < 0.05)
          registry.get<TH2>(HIST("pTinvMassLcReso"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPiKP(lccandidate));
        if (mDiffKStar892 < 0.046 * 2.)
          registry.get<TH2>(HIST("pTinvMassLcKStar"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPiKP(lccandidate));
        if (mDiffDelta1232 < 0.117)
          registry.get<TH2>(HIST("pTinvMassLcDelta"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPiKP(lccandidate));
        if (mDiffLambda1520 < 0.032)
          registry.get<TH2>(HIST("pTinvMassLcLambda1520"))->Fill(lccandidate.pt(), hfHelper.invMassLcToPiKP(lccandidate));

        if (abs(hfHelper.invMassLcToPiKP(lccandidate) - o2::constants::physics::MassLambdaCPlus) < 0.05) {
          registry.get<TH2>(HIST("pTinvMassKStar"))->Fill(lccandidate.pt(), InvariantMass(track1, track0, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus));
          registry.get<TH2>(HIST("pTinvMassDelta"))->Fill(lccandidate.pt(), InvariantMass(track2, track0, o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus));
          registry.get<TH2>(HIST("pTinvMassLambda1520"))->Fill(lccandidate.pt(), InvariantMass(track2, track1, o2::constants::physics::MassProton, o2::constants::physics::MassKPlus));
        }
      }
    } // Lambda_c candidates loop for cross checks

    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << hf_cand_lb::DecayType::LbToLcPi)) { // This should never be true as the loop is over Lb candidates
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLb(candidate)) > yCandRecoMax) {
        continue;
      }

      auto candLc = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
      auto candPi = candidate.prong1_as<TracksWExt>();

      auto track0Lc = candLc.prong0_as<TracksWExt>();
      auto track1Lc = candLc.prong1_as<TracksWExt>();
      auto track2Lc = candLc.prong2_as<TracksWExt>();
      double d0resolution0 = candLc.errorImpactParameter0();
      double d0resolution1 = candLc.errorImpactParameter1();
      double d0resolution2 = candLc.errorImpactParameter2();
      if (!PassesImpactParameterResolution(track0Lc.pt(), d0resolution0))
        continue; // Track quality selection here - move to Lb selector later
      if (!PassesImpactParameterResolution(track1Lc.pt(), d0resolution1))
        continue;
      if (!PassesImpactParameterResolution(track2Lc.pt(), d0resolution2))
        continue;
      if (!PassesImpactParameterResolution(candPi.pt(), candidate.errorImpactParameter1()))
        continue;
      double DCA0 = candLc.impactParameter0();
      double DCA1 = candLc.impactParameter1();
      double DCA2 = candLc.impactParameter2();
      if (DCA0 > 0.2 || DCA1 > 0.2 || DCA2 > 0.2)
        continue; // reject clear strangeness feed down - might also be done with the likelihood
      double likelihoodRatio = LogLikelihoodRatioSingleTrackDCA(DCA0, d0resolution0, DCALengthParameter) + LogLikelihoodRatioSingleTrackDCA(DCA1, d0resolution1, DCALengthParameter) + LogLikelihoodRatioSingleTrackDCA(DCA2, d0resolution2, DCALengthParameter);
      if (likelihoodRatio < minLikelihoodRatio)
        continue; // Larger likelihood means more likely to be signal
      if (candidate.impactParameter1() < 0.005)
        continue;
      double LcMass = 0.;
      if (candLc.isSelLcToPKPi())
        LcMass = hfHelper.invMassLcToPKPi(candLc);
      if (candLc.isSelLcToPiKP())
        LcMass = hfHelper.invMassLcToPiKP(candLc);
      if (abs(LcMass - o2::constants::physics::MassLambdaCPlus) > 0.1)
        continue;
      if (candidate.cpa() < 0.7)
        continue;
      if (candidate.errorDecayLengthXY() > 0.01)
        continue;
      if (candidate.decayLengthXY() < 0.03)
        continue;
      if (candidate.errorDecayLength() > 0.015)
        continue;
      if (candidate.decayLength() < 0.04)
        continue;
      registry.get<TH2>(HIST("pTinvMassLb"))->Fill(candidate.pt(), InvariantMass(candLc, candPi, o2::constants::physics::MassLambdaCPlus, o2::constants::physics::MassPiPlus));

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
      if (candPi.sign() < 0) {
        registry.fill(HIST("hInvMassLc"), hfHelper.invMassLcToPKPi(candLc), candidate.pt());
      }
    } // candidate loop
  }
};

/// Lb MC analysis and fill histograms
struct HfTaskLbMc {
  Configurable<int> selectionFlagLb{"selectionFlagLb", 1, "Selection Flag for Lb"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lb_to_lc_pi::vecBinsPt}, "pT bin limits"};

  Service<o2::framework::O2DatabasePDG> pdg;
  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lb::isSelLbToLcPi >= selectionFlagLb);

  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "Lb candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtRecBg", "Lb candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtGenSig", "Lb candidates (matched);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 10.}}}},
     {"hPtGen", "MC particles (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}}}};

  void init(InitContext&)
  {
    registry.add("hEtaGen", "MC particles (matched);#Lambda_{b}^{0} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYGen", "MC particles (matched);#Lambda_{b}^{0} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0Gen", "MC particles (matched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng0Gen", "MC particles (matched);prong 0 (#Lambda_{c}^{+}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng0Gen", "MC particles (matched);prong 0 (#Lambda_{b}^{0}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaProng1Gen", "MC particles (matched);prong 1 (#pi^{-}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAxyRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate CPAxy;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPALcRecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPALcRecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) cosine of pointing angle;entries", {HistType::kTH2F, {{220, 0., 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hRapidityRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate #it{#y};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "#Lambda_{b}^{0} candidates (matched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 1 (#pi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecSig", "#Lambda_{b}^{0} candidates (matched);inv. mass #Lambda_{c}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.00}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "#Lambda_{b}^{0} candidates (unmatched);inv. mass #Lambda_{c}^{+}#pi^{+} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 4.0, 7.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "#Lambda_{b}^{0} candidates (matched);prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "#Lambda_{b}^{0} candidates (matched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "#Lambda_{b}^{0} candidates (unmatched);prong 1 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthXYRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthLcRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthNormRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdLbRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParProdLbRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} candidate impact parameter product ;entries", {HistType::kTH2F, {{100, -0.5, 0.5}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCARecSig", "#Lambda_{b}^{0} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCARecBg", "#Lambda_{b}^{0} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hThetaStarRecSig", "#Lambda_{b}^{0} candidates (matched);#Lambda_{b}^{0} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hThetaStarRecBg", "#Lambda_{b}^{0} candidates (unmatched);#Lambda_{b}^{0} #cos(#theta^{*});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCandLb, aod::HfSelLbToLcPi, aod::HfCandLbMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandLbMcGen> const& mcParticles,
               aod::TracksWMc const& tracks,
               aod::HfCand3Prong const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << hf_cand_lb::DecayType::LbToLcPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      auto candLc = candidate.prong0_as<aod::HfCand3Prong>();
      if (std::abs(candidate.flagMcMatchRec()) == 1 << hf_cand_lb::DecayType::LbToLcPi) {

        auto indexMother = RecoDecay::getMother(mcParticles, candidate.prong1_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandLbMcGen>>(), o2::constants::physics::Pdg::kLambdaB0, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecSig"), hfHelper.yLb(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), hfHelper.invMassLbToLcPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hImpParProdLbRecSig"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecSig"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("hCPALcRecSig"), candLc.cpa(), candidate.pt());
        registry.fill(HIST("hDecLengthLcRecSig"), candLc.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), candidate.pt());
        // registry.fill(HIST("hThetaStarRecSig"), candidate.cosThetaStar(), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());
        registry.fill(HIST("hRapidityRecBg"), hfHelper.yLb(candidate), candidate.pt());
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), hfHelper.invMassLbToLcPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hImpParProdLbRecBg"), candidate.impactParameterProduct(), candidate.pt());
        registry.fill(HIST("hDecLengthNormRecBg"), candidate.decayLengthXYNormalised(), candidate.pt());
        registry.fill(HIST("hCPALcRecBg"), candLc.cpa(), candidate.pt());
        registry.fill(HIST("hDecLengthLcRecBg"), candLc.decayLength(), candidate.pt());
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), candidate.pt());
        // registry.fill(HIST("hThetaStarRecBg"), candidate.cosThetaStar(), candidate.pt());
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << hf_cand_lb::DecayType::LbToLcPi) {

        auto yParticle = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassLambdaB0);
        if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
          continue;
        }

        float ptProngs[2], yProngs[2], etaProngs[2];
        int counter = 0;
        for (const auto& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(std::array{daught.px(), daught.py(), daught.pz()}, pdg->Mass(daught.pdgCode()));
          counter++;
        }

        registry.fill(HIST("hPtProng0Gen"), ptProngs[0], particle.pt());
        registry.fill(HIST("hPtProng1Gen"), ptProngs[1], particle.pt());
        registry.fill(HIST("hYProng0Gen"), yProngs[0], particle.pt());
        registry.fill(HIST("hYProng1Gen"), yProngs[1], particle.pt());
        registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], particle.pt());
        registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], particle.pt());

        //  if (yCandMax >= 0. && (std::abs(yProngs[0]) > yCandMax || std::abs(yProngs[1]) > yCandMax))
        //    continue;

        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hYGen"), yParticle, particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());
      }
    } // gen
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  const bool doMC = cfgc.options().get<bool>("doMC");
  workflow.push_back(adaptAnalysisTask<HfTaskLb>(cfgc));
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskLbMc>(cfgc));
  }
  return workflow;
}
