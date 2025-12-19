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

/// \file taskCd.cxx
/// \brief Cd± → d± K∓ π±  analysis task
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg Universiity

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <vector> // std::vector

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::hf_evsel;

struct HfTaskCd {
  Configurable<int> selectionFlagCd{"selectionFlagCd", 1, "Selection Flag for Cd"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_cd_to_de_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<bool> fillTHn{"fillTHn", false, "fill THn"};

  SliceCache cache;

  using CollisionsWEvSel = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithEvSelFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithEvSelFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

  using CdCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelCd>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_cd::isSelCdToDeKPi >= selectionFlagCd;
  Preslice<aod::HfCand3Prong> candCdPerCollision = aod::hf_cand::collisionId;

  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {72, 0, 36}, ""};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {400, 2.4, 4.4}, ""};
  ConfigurableAxis thnConfigAxisPtProng{"thnConfigAxisPtProng", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisChi2PCA{"thnConfigAxisChi2PCA", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisCentrality{"thnConfigAxisCentrality", {100, 0, 100}, ""};

  HistogramRegistry registry{
    "registry",
    {/// mass candidate
     {"Data/hMass", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.4, 4.4}}}},
     /// pT
     {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     /// DCAxy to prim. vertex prongs
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     /// decay length candidate
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// decay length xy candidate
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// cosine of pointing angle
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// cosine of pointing angle xy
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// Chi 2 PCA to sec. vertex
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{400, 0., 20.}}}},
     /// eta
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     /// phi
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}}}};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    std::array<bool, 3> doprocess{doprocessDataStd, doprocessDataStdWithFT0C, doprocessDataStdWithFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }
    /// mass candidate
    registry.add("Data/hMassVsPtVsNPvContributors", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2}); p_{T}; Number of PV contributors", {HistType::kTH3F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}, {500, 0., 5000.}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// impact parameter error
    registry.add("Data/hImpParErrProng0", "3-prong candidates;prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng1", "3-prong candidates;prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng2", "3-prong candidates;prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    if (fillTHn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (de K #pi) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T}(C_{d}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng0{thnConfigAxisPtProng, "#it{p}_{T}(prong0) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng1{thnConfigAxisPtProng, "#it{p}_{T}(prong1) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng2{thnConfigAxisPtProng, "#it{p}_{T}(prong2) (GeV/#it{c})"};
      const AxisSpec thnAxisChi2PCA{thnConfigAxisChi2PCA, "Chi2PCA to sec. vertex (cm)"};
      const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length (cm)"};
      const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle"};
      const AxisSpec thnAxisCentrality{thnConfigAxisCentrality, "centrality (FT0C)"};

      std::vector axesStd{thnAxisMass, thnAxisPt, thnAxisPtProng0, thnAxisPtProng1, thnAxisPtProng2, thnAxisChi2PCA, thnAxisDecLength, thnAxisCPA, thnAxisCentrality};
      registry.add("hnCdVars", "THn for Reconstructed Cd candidates for data", HistType::kTHnSparseF, axesStd);
    }
  }

  /// Fill histograms for real data
  template <typename CollType, typename CandType>
  void fillHistosData(CollType const& collision, CandType const& candidates)
  {
    auto thisCollId = collision.globalIndex();
    auto groupedCdCandidates = candidates.sliceBy(candCdPerCollision, thisCollId);
    auto numPvContributors = collision.numContrib();

    for (const auto& candidate : groupedCdCandidates) {
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::CdToDeKPi)) {
        continue;
      }

      const auto pt = candidate.pt();
      const auto ptProng0 = candidate.ptProng0();
      const auto ptProng1 = candidate.ptProng1();
      const auto ptProng2 = candidate.ptProng2();
      const auto decayLength = candidate.decayLength();
      const auto decayLengthXY = candidate.decayLengthXY();
      const auto chi2PCA = candidate.chi2PCA();
      const auto cpa = candidate.cpa();
      const auto cpaXY = candidate.cpaXY();

      if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
        registry.fill(HIST("Data/hMass"), HfHelper::invMassCdToDeKPi(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), HfHelper::invMassCdToDeKPi(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), HfHelper::invMassCdToDeKPi(candidate), pt);
      }
      if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
        registry.fill(HIST("Data/hMass"), HfHelper::invMassCdToPiKDe(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), HfHelper::invMassCdToPiKDe(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), HfHelper::invMassCdToPiKDe(candidate), pt);
      }
      registry.fill(HIST("Data/hPt"), pt);
      registry.fill(HIST("Data/hPtProng0"), ptProng0);
      registry.fill(HIST("Data/hPtProng1"), ptProng1);
      registry.fill(HIST("Data/hPtProng2"), ptProng2);
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0());
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1());
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2());
      registry.fill(HIST("Data/hd0VsPtProng0"), candidate.impactParameter0(), pt);
      registry.fill(HIST("Data/hd0VsPtProng1"), candidate.impactParameter1(), pt);
      registry.fill(HIST("Data/hd0VsPtProng2"), candidate.impactParameter2(), pt);
      registry.fill(HIST("Data/hDecLength"), decayLength);
      registry.fill(HIST("Data/hDecLengthVsPt"), decayLength, pt);
      registry.fill(HIST("Data/hDecLengthxy"), decayLengthXY);
      registry.fill(HIST("Data/hDecLengthxyVsPt"), decayLengthXY, pt);
      registry.fill(HIST("Data/hCPA"), cpa);
      registry.fill(HIST("Data/hCPAVsPt"), cpa, pt);
      registry.fill(HIST("Data/hCPAxy"), cpaXY);
      registry.fill(HIST("Data/hCPAxyVsPt"), cpaXY, pt);
      registry.fill(HIST("Data/hDca2"), chi2PCA);
      registry.fill(HIST("Data/hDca2VsPt"), chi2PCA, pt);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), pt);
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPhiVsPt"), candidate.phi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelCdToDeKPi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelCdToPiKDe(), pt);
      registry.fill(HIST("Data/hImpParErrProng0"), candidate.errorImpactParameter0(), pt);
      registry.fill(HIST("Data/hImpParErrProng1"), candidate.errorImpactParameter1(), pt);
      registry.fill(HIST("Data/hImpParErrProng2"), candidate.errorImpactParameter2(), pt);

      if (fillTHn) {
        float const cent = o2::hf_centrality::getCentralityColl(collision);
        double massCd(-1);
        if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
          massCd = HfHelper::invMassCdToDeKPi(candidate);
          std::vector<double> valuesToFill{massCd, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        }
        if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
          massCd = HfHelper::invMassCdToPiKDe(candidate);
          std::vector<double> valuesToFill{massCd, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        }
      }
    }
  }
  /// Run the analysis on real data
  template <typename CollType, typename CandType>
  void runAnalysisPerCollisionData(CollType const& collisions,
                                   CandType const& candidates)
  {

    for (const auto& collision : collisions) {
      fillHistosData(collision, candidates);
    }
  }

  void processDataStd(CollisionsWEvSel const& collisions,
                      CdCandidates const& selectedCdCandidates,
                      aod::Tracks const&)
  {
    runAnalysisPerCollisionData(collisions, selectedCdCandidates);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStd, "Process Data with the standard method", true);

  void processDataStdWithFT0C(CollisionsWithEvSelFT0C const& collisions,
                              CdCandidates const& selectedCdCandidates,
                              aod::Tracks const&)
  {
    runAnalysisPerCollisionData(collisions, selectedCdCandidates);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStdWithFT0C, "Process real data with the standard method and with FT0C centrality", false);

  void processDataStdWithFT0M(CollisionsWithEvSelFT0M const& collisions,
                              CdCandidates const& selectedCdCandidates,
                              aod::Tracks const&)
  {
    runAnalysisPerCollisionData(collisions, selectedCdCandidates);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStdWithFT0M, "Process real data with the standard method and with FT0M centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCd>(cfgc)};
}
