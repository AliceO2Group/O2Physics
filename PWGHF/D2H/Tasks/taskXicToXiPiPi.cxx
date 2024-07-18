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

/// \file taskXicToXiPiPi.cxx
/// \brief Ξc± → (Ξ∓ → (Λ → p π∓) π∓) π± π± analysis task
/// \note adapted from taskBs.cxx
///
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Xic analysis task
struct HfTaskXicToXiPiPi {
  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi_pi::vecBinsPt}, "pT bin limits"};
  // MC checks
  Configurable<bool> checkDecayTypeMc{"checkDecayTypeMc", false, "Flag to enable DecayType histogram"};

  Service<o2::framework::O2DatabasePDG> pdg;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic);

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "#Xi^{#plus}_{c} candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{400, 0., 40.}}}},
     {"hPtProng0", "#Xi^{#plus}_{c} candidates;prong 0 (#Xi^{#minus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}}},
     {"hPtProng1", "#Xi^{#plus}_{c} candidates;prong 1 (#pi^{#plus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 20.}}}},
     {"hPtProng2", "#Xi^{#plus}_{c} candidates;prong 2 (#pi^{#plus}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 20.}}}}}};

  void init(InitContext const&)
  {
    static const AxisSpec axisMassXic = {300, 1.8, 3.0, "inv. mass (GeV/#it{c}^{2})"};
    static const AxisSpec axisMassXiRes = {300, 1.0, 2.0, "inv. mass (GeV/#it{c}^{2})"};
    static const AxisSpec axisPt = {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"};

    registry.add("hEta", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hRapidity", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    registry.add("hCPA", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hCPAxy", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hMass", "#Xi^{#plus}_{c} candidates;inv. mass #Xi^{#mp} #pi^{#pm} #pi^{#pm} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisMassXic, axisPt}});
    registry.add("hDecLength", "#Xi^{#plus}_{c} candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, axisPt}});
    registry.add("hErrDecLength", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, axisPt}});
    registry.add("hDecLengthXY", "#Xi^{#plus}_{c} candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 0.4}, axisPt}});
    registry.add("hErrDecLengthXY", "#Xi^{#plus}_{c} candidates;#Xi^{#plus}_{c} candidate decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, axisPt}});
    registry.add("hd0Prong0", "#Xi^{#plus}_{c} candidates;prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong1", "#Xi^{#plus}_{c} candidates;prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hd0Prong2", "#Xi^{#plus}_{c} candidates;prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -0.05, 0.05}, axisPt}});
    registry.add("hImpParErr", "#Xi^{#plus}_{c} candidates;prongs impact parameter error (cm);entries", {HistType::kTH2F, {{200, -0.1, 0.1}, axisPt}});
    registry.add("hChi2PCA", "#Xi^{#plus}_{c} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.5}, axisPt}});
    registry.add("hCPAXi", "#Xi^{#plus}_{c} candidates;#Xi^{#minus} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hCPAxyXi", "#Xi^{#plus}_{c} candidates;#Xi^{#minus} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hCPALambda", "#Xi^{#plus}_{c} candidates;#Lambda candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hCPAxyLambda", "#Xi^{#plus}_{c} candidates;#Lambda candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
    registry.add("hMassXiPi1", "#Xi^{#plus}_{c} candidates;inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisMassXiRes, axisPt}});
    registry.add("hMassXiPi2", "#Xi^{#plus}_{c} candidates;inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisMassXiRes, axisPt}});

    if (doprocessMc) {
      // MC reconstructed
      registry.add("hPtGenSig", "#Xi^{#plus}_{c} candidates (gen+rec);candidate #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hPtRecSig", "#Xi^{#plus}_{c} candidates (matched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hPtRecBg", "#Xi^{#plus}_{c} candidates (unmatched);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hPtProng0RecSig", "#Xi^{#plus}_{c} candidates (matched);prong 0 (#Xi^{#mp}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hPtProng0RecBg", "#Xi^{#plus}_{c} candidates (unmatched);prong 0 (#Xi^{#mp}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hPtProng1RecSig", "#Xi^{#plus}_{c} candidates (matched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      registry.add("hPtProng1RecBg", "#Xi^{#plus}_{c} candidates (unmatched);prong 1 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      registry.add("hPtProng2RecSig", "#Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      registry.add("hPtProng2RecBg", "#Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 20.}}});
      registry.add("hEtaRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hEtaRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hRapidityRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hRapidityRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hCPARecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      registry.add("hCPARecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      registry.add("hCPAxyRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate CPAxy;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      registry.add("hCPAxyRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate CPAxy;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      registry.add("hMassRecSig", "#Xi^{#plus}_{c} candidates (matched);inv. mass  #Xi^{#mp} #pi^{#pm} #pi^{#pm} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.8, 3.0}, axisPt}});
      registry.add("hMassRecBg", "#Xi^{#plus}_{c} candidates (unmatched);inv. mass  #Xi^{#mp} #pi^{#pm} #pi^{#pm} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.8, 3.0}, axisPt}});
      registry.add("hDecLengthRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hDecLengthRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hErrDecLengthRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hErrDecLengthRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hDecLengthXYRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hDecLengthXYRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hErrDecLengthXYRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#plus}_{c} candidate decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hErrDecLengthXYRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#plus}_{c} candidate decay length xy(cm);entries", {HistType::kTH2F, {{100, 0., 0.5}, axisPt}});
      registry.add("hd0Prong0RecSig", "#Xi^{#plus}_{c} candidates (matched);prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
      registry.add("hd0Prong0RecBg", "#Xi^{#plus}_{c} candidates (unmatched);prong 0 (#Xi^{#mp}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
      registry.add("hd0Prong1RecSig", "#Xi^{#plus}_{c} candidates (matched);prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
      registry.add("hd0Prong1RecBg", "#Xi^{#plus}_{c} candidates (unmatched);prong 1 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
      registry.add("hd0Prong2RecSig", "#Xi^{#plus}_{c} candidates (matched);prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
      registry.add("hd0Prong2RecBg", "#Xi^{#plus}_{c} candidates (unmatched);prong 2 (#pi^{#pm}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{200, -0.05, 0.05}, axisPt}});
      registry.add("hImpParErrRecSig", "#Xi^{#plus}_{c} candidates (matched);prongs impact parameter error (cm);entries", {HistType::kTH2F, {{100, -0.1, 0.1}, axisPt}});
      registry.add("hImpParErrRecBg", "#Xi^{#plus}_{c} candidates (unmatched);prongs impact parameter error (cm);entries", {HistType::kTH2F, {{100, -0.1, 0.1}, axisPt}});
      registry.add("hChi2PCARecSig", "#Xi^{#plus}_{c} candidates (matched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});
      registry.add("hChi2PCARecBg", "#Xi^{#plus}_{c} candidates (unmatched);sum of distances of the secondary vertex to its prongs;entries", {HistType::kTH2F, {{240, -0.01, 0.1}, axisPt}});
      registry.add("hCPAXiRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#minus} cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      registry.add("hCPAXiRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#minus} cosine of pointing angle;entries", {HistType::kTH2F, {{220, -1.1, 1.1}, axisPt}});
      registry.add("hCPAxyXiRecSig", "#Xi^{#plus}_{c} candidates (matched);#Xi^{#minus} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
      registry.add("hCPAxyXiRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Xi^{#minus} candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
      registry.add("hCPALambdaRecSig", "#Xi^{#plus}_{c} candidates (matched);#Lambda candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
      registry.add("hCPALambdaRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Lambda candidate cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
      registry.add("hCPAxyLambdaRecSig", "#Xi^{#plus}_{c} candidates (matched);#Lambda candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
      registry.add("hCPAxyLambdaRecBg", "#Xi^{#plus}_{c} candidates (unmatched);#Lambda candidate cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, axisPt}});
      registry.add("hMassXiPi1RecSig", "#Xi^{#plus}_{c} candidates (matched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
      registry.add("hMassXiPi1RecBg", "#Xi^{#plus}_{c} candidates (unmatched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 1) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
      registry.add("hMassXiPi2RecSig", "#Xi^{#plus}_{c} candidates (matched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
      registry.add("hMassXiPi2RecBg", "#Xi^{#plus}_{c} candidates (unmatched);inv. mass #Xi^{#mp} #pi^{#pm} (prong 2) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{300, 1.0, 2.0}, axisPt}});
      // MC generated
      registry.add("hPtProng0Gen", "MC particles (generated);prong 0 (#Xi^{#mp}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{300, 0., 30.}, axisPt}});
      registry.add("hPtProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      registry.add("hPtProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{p}_{T}^{gen} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 20.}, axisPt}});
      registry.add("hEtaProng0Gen", "MC particles (generated);prong 0 (#Xi^{#mp}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
      registry.add("hEtaProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
      registry.add("hEtaProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
      registry.add("hYProng0Gen", "MC particles (generated);prong 0 (#Xi^{#mp}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
      registry.add("hYProng1Gen", "MC particles (generated);prong 1 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
      registry.add("hYProng2Gen", "MC particles (generated);prong 2 (#pi^{#pm}) #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2, 2}, axisPt}});
      registry.add("hPtGen", "MC particles (generated);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hEtaGen", "MC particles (generated);#Xi^{#plus}_{c} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hYGen", "MC particles (generated);#Xi^{#plus}_{c} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 30.}}});
      registry.add("hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#Xi^{#plus}_{c} candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
      registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#Xi^{#plus}_{c} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, axisPt}});
    }

    if (checkDecayTypeMc) {
      constexpr uint8_t kNBinsDecayTypeMc = hf_cand_xic_to_xi_pi_pi::DecayType::NDecayType + 1;
      TString labels[kNBinsDecayTypeMc];
      labels[hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi] = "#Xi^{+}_{c} #rightarrow #Xi^{#minus} #pi^{#plus}) #pi^{#plus}";
      labels[hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi] = "#Xi^{+}_{c} #rightarrow #Xi(1530)^{0} #pi^{#plus} #rightarrow #Xi^{#minus} #pi^{#plus}) #pi^{#plus}";
      labels[hf_cand_xic_to_xi_pi_pi::DecayType::NDecayType] = "Other decays";
      static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
      registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassXic, axisPt}});
      for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
        registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
      }
    }
  }

  /// Selection of Xic daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of Xic prong
  /// \param ptProng is the pT of Xic prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  void process(soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>> const& candidates)
  {
    for (const auto& candidate : candidates) {
      auto yCandXic = candidate.y(o2::constants::physics::MassXiCPlus);
      if (yCandRecoMax >= 0. && std::abs(yCandXic) > yCandRecoMax) {
        continue;
      }

      auto ptCandXic = candidate.pt();

      registry.fill(HIST("hPtCand"), ptCandXic);
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("hEta"), candidate.eta(), ptCandXic);
      registry.fill(HIST("hRapidity"), yCandXic, ptCandXic);
      registry.fill(HIST("hCPA"), candidate.cpa(), ptCandXic);
      registry.fill(HIST("hCPAxy"), candidate.cpaXY(), ptCandXic);
      registry.fill(HIST("hMass"), candidate.invMassXic(), ptCandXic);
      registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandXic);
      registry.fill(HIST("hErrDecLength"), candidate.errorDecayLength(), ptCandXic);
      registry.fill(HIST("hDecLengthXY"), candidate.decayLengthXY(), ptCandXic);
      registry.fill(HIST("hErrDecLengthXY"), candidate.errorDecayLengthXY(), ptCandXic);
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandXic);
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandXic);
      registry.fill(HIST("hd0Prong2"), candidate.impactParameter2(), ptCandXic);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), ptCandXic);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), ptCandXic);
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter2(), ptCandXic);
      registry.fill(HIST("hChi2PCA"), candidate.chi2PCA(), ptCandXic);
      registry.fill(HIST("hCPAXi"), candidate.cosPaXi(), ptCandXic);
      registry.fill(HIST("hCPAxyXi"), candidate.cosPaXYXi(), ptCandXic);
      registry.fill(HIST("hCPALambda"), candidate.cosPaLambda(), ptCandXic);
      registry.fill(HIST("hCPAxyLambda"), candidate.cosPaLambda(), ptCandXic);
      registry.fill(HIST("hMassXiPi1"), candidate.invMassXiPi0(), ptCandXic);
      registry.fill(HIST("hMassXiPi2"), candidate.invMassXiPi1(), ptCandXic);
    } // candidate loop
  }   // process

  /// MC analysis and fill histograms
  void processMc(soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi, aod::HfCandXicMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandXicMcGen> const& mcParticles,
                 aod::TracksWMc const&)
  {
    std::vector<int> arrDaughIndex;

    // MC rec
    for (const auto& candidate : candidates) {
      auto yCandXic = candidate.y(o2::constants::physics::MassXiCPlus);
      if (yCandRecoMax >= 0. && std::abs(yCandXic) > yCandRecoMax) {
        continue;
      }

      auto ptCandXic = candidate.pt();
      int flagMcMatchRecXic = std::abs(candidate.flagMcMatchRec());

      if (TESTBIT(flagMcMatchRecXic, hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) || TESTBIT(flagMcMatchRecXic, hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
        auto indexMother = RecoDecay::getMother(mcParticles, candidate.pi0_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandXicMcGen>>(), o2::constants::physics::Pdg::kXiCPlus, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);

        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), ptCandXic);
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2RecSig"), candidate.ptProng2());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandXic);
        registry.fill(HIST("hRapidityRecSig"), yCandXic, ptCandXic);
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), ptCandXic);
        registry.fill(HIST("hCPAxyRecSig"), candidate.cpaXY(), ptCandXic);
        registry.fill(HIST("hMassRecSig"), candidate.invMassXic(), ptCandXic);
        registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandXic);
        registry.fill(HIST("hErrDecLengthRecSig"), candidate.errorDecayLength(), ptCandXic);
        registry.fill(HIST("hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandXic);
        registry.fill(HIST("hErrDecLengthXYRecSig"), candidate.errorDecayLengthXY(), ptCandXic);
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandXic);
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandXic);
        registry.fill(HIST("hd0Prong2RecSig"), candidate.impactParameter2(), ptCandXic);
        registry.fill(HIST("hImpParErrRecSig"), candidate.errorImpactParameter0(), ptCandXic);
        registry.fill(HIST("hImpParErrRecSig"), candidate.errorImpactParameter1(), ptCandXic);
        registry.fill(HIST("hImpParErrRecSig"), candidate.errorImpactParameter2(), ptCandXic);
        registry.fill(HIST("hChi2PCARecSig"), candidate.chi2PCA(), ptCandXic);
        registry.fill(HIST("hCPAXiRecSig"), candidate.cosPaXi(), ptCandXic);
        registry.fill(HIST("hCPAxyXiRecSig"), candidate.cosPaXYXi(), ptCandXic);
        registry.fill(HIST("hCPALambdaRecSig"), candidate.cosPaLambda(), ptCandXic);
        registry.fill(HIST("hCPAxyLambdaRecSig"), candidate.cosPaLambda(), ptCandXic);
      } else {
        registry.fill(HIST("hPtRecBg"), ptCandXic);
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1());
        registry.fill(HIST("hPtProng2RecBg"), candidate.ptProng2());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandXic);
        registry.fill(HIST("hRapidityRecBg"), yCandXic, ptCandXic);
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), ptCandXic);
        registry.fill(HIST("hCPAxyRecBg"), candidate.cpaXY(), ptCandXic);
        registry.fill(HIST("hMassRecBg"), candidate.invMassXic(), ptCandXic);
        registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandXic);
        registry.fill(HIST("hErrDecLengthRecBg"), candidate.errorDecayLength(), ptCandXic);
        registry.fill(HIST("hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandXic);
        registry.fill(HIST("hErrDecLengthXYRecBg"), candidate.errorDecayLengthXY(), ptCandXic);
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandXic);
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandXic);
        registry.fill(HIST("hd0Prong2RecBg"), candidate.impactParameter2(), ptCandXic);
        registry.fill(HIST("hImpParErrRecBg"), candidate.errorImpactParameter0(), ptCandXic);
        registry.fill(HIST("hImpParErrRecBg"), candidate.errorImpactParameter1(), ptCandXic);
        registry.fill(HIST("hImpParErrRecBg"), candidate.errorImpactParameter2(), ptCandXic);
        registry.fill(HIST("hChi2PCARecBg"), candidate.chi2PCA(), ptCandXic);
        registry.fill(HIST("hCPAXiRecBg"), candidate.cosPaXi(), ptCandXic);
        registry.fill(HIST("hCPAxyXiRecBg"), candidate.cosPaXYXi(), ptCandXic);
        registry.fill(HIST("hCPALambdaRecBg"), candidate.cosPaLambda(), ptCandXic);
        registry.fill(HIST("hCPAxyLambdaRecBg"), candidate.cosPaLambda(), ptCandXic);
      }

      if (checkDecayTypeMc) {
        if (TESTBIT(flagMcMatchRecXic, hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi)) {
          registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi, candidate.invMassXic(), ptCandXic);
        } else if (TESTBIT(flagMcMatchRecXic, hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
          registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi, candidate.invMassXic(), ptCandXic);
        } else {
          registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_xic_to_xi_pi_pi::DecayType::NDecayType, candidate.invMassXic(), ptCandXic);
        }
      }
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      if (TESTBIT(std::abs(particle.flagMcMatchGen()), hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) || TESTBIT(std::abs(particle.flagMcMatchGen()), hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi)) {
        arrDaughIndex.clear();

        auto ptParticle = particle.pt();
        auto yParticle = RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus);
        if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
          continue;
        }

        // get kinematic variables of Ξ π π
        std::array<float, 3> ptProngs;
        std::array<float, 3> yProngs;
        std::array<float, 3> etaProngs;
        int counter = 0;
        RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, 2);
        for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
          auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
          ptProngs[counter] = daughI.pt();
          etaProngs[counter] = daughI.eta();
          yProngs[counter] = RecoDecay::y(daughI.pVector(), pdg->Mass(daughI.pdgCode()));
          counter++;
        }

        registry.fill(HIST("hPtProng0Gen"), ptProngs[0], ptParticle);
        registry.fill(HIST("hPtProng1Gen"), ptProngs[1], ptParticle);
        registry.fill(HIST("hPtProng2Gen"), ptProngs[2], ptParticle);
        registry.fill(HIST("hEtaProng0Gen"), etaProngs[0], ptParticle);
        registry.fill(HIST("hEtaProng1Gen"), etaProngs[1], ptParticle);
        registry.fill(HIST("hEtaProng2Gen"), etaProngs[2], ptParticle);
        registry.fill(HIST("hYProng0Gen"), yProngs[0], ptParticle);
        registry.fill(HIST("hYProng1Gen"), yProngs[1], ptParticle);
        registry.fill(HIST("hYProng2Gen"), yProngs[2], ptParticle);
        registry.fill(HIST("hPtGen"), ptParticle);
        registry.fill(HIST("hYGen"), yParticle, ptParticle);
        registry.fill(HIST("hEtaGen"), particle.eta(), ptParticle);

        // reject Xic daughters that are not in geometrical acceptance
        if (!isProngInAcceptance(etaProngs[0], ptProngs[0]) || !isProngInAcceptance(etaProngs[1], ptProngs[1]) || !isProngInAcceptance(etaProngs[2], ptProngs[2])) {
          continue;
        }
        registry.fill(HIST("hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("hEtaGenWithProngsInAcceptance"), particle.eta(), ptParticle);
        registry.fill(HIST("hYGenWithProngsInAcceptance"), yParticle, ptParticle);
      }
    } // gen
  }   // process
  PROCESS_SWITCH(HfTaskXicToXiPiPi, processMc, "Process MC", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXicToXiPiPi>(cfgc)};
}
