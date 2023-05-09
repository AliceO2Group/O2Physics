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

/// \file correlatorDMesonPairs.cxx
/// \brief D0(bar) and DPlus(Minus) correlator task - data-like, MC-reco and MC-kine analyses.
///
/// \author Fabio Colamaria <fabio.colamaria@ba.infn.it>, INFN Bari
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::constants::math;

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

// histogram binning definition
const int massAxisBins = 120;
const double massAxisMin = 1.5848;
const double massAxisMax = 2.1848;
const int phiAxisBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = o2::constants::math::TwoPI;
const int yAxisBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

using McParticlesPlus2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
using McParticlesPlus3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

struct HfCorrelatorDMesonPairs {
  SliceCache cache;
  Preslice<aod::HfCand2Prong> perCol2Prong = aod::hf_cand::collisionId;
  Preslice<aod::HfCand3Prong> perCol3Prong = aod::hf_cand::collisionId;

  Produces<aod::D0Pair> entryD0Pair;
  Produces<aod::D0PairRecoInfo> entryD0PairRecoInfo;
  Produces<aod::DPlusPair> entryDPlusPair;
  Produces<aod::DPlusPairRecoInfo> entryDPlusPairRecoInfo;

  Configurable<bool> verbose{"verbose", false, "Enable debug mode"};
  Configurable<int> particlePdgCode{"particlePdgCode", 421, "PDG code of particle to analyse"};
  // Selection Flags
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagDplus{"selectionFlagDPlus", 1, "Selection Flag for DPlus"};
  // Kinematics
  Configurable<double> etaCutMax{"etaCutMax", 4.0, "max. eta accepted"};
  Configurable<double> etaCutMin{"etaCutMin", -4.0, "min. eta accepted"};
  Configurable<double> dcaXYCut{"dcaXYCut", 0.0025, "cut in dcaXY"};
  Configurable<double> dcaZCut{"dcaZCut", 0.0025, "cut in dcaZ"};

  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<double> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots"};

  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> selectedD0candidatesMc = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> selectedDPlusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>> selectedDPlusCandidatesMc = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "D Meson pair candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "D Meson pair candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "D Meson pair candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "D Meson pair candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "D Meson pair candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "D Meson pair candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "D Meson pair candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     // Mc Reco
     {"hOriginMcRec", "D Meson pair candidates - MC reco;prompt vs. non-prompt;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hMatchedMcRec", "D Meson pair candidates - MC reco;MC Matched;entries", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hPtCandMcRec", "D Meson pair candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0McRec", "D Meson pair candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1McRec", "D Meson pair candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMcRec", "D Meson pair candidates - MC reco;selection status;entries", {HistType::kTH1F, {{301, -0.5, 300.5}}}},
     {"hEtaMcRec", "D Meson pair candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMcRec", "D Meson pair candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMcRec", "D Meson pair candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     // Mc Gen
     {"hMcEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMcGen", "D Meson pair particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMcGen", "D Meson pair particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMcGen", "D Meson pair particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMcGen", "D Meson pair candidates - MC gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMatchedMcGen", "D Meson pair candidates - MC gen;MC Matched;entries", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
     {"hOriginMcGen", "D Meson pair candidates - MC gen;prompt vs. non-prompt;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hCountD0D0barPerEvent", "D Meson pair particles - MC gen;Number per event;entries", {HistType::kTH1F, {{20, 0., 20.}}}}}};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMass", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0", "D0,D0bar candidates;inv. mass D0 only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0bar", "D0,D0bar candidates;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0McRecSig", "D0 signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0McRecRefl", "D0 reflection candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0McRecBkg", "D0 background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMcRecSig", "D0bar signal candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMcRecRefl", "D0bar reflection candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMcRecBkg", "D0bar background candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // DPlus plots
    registry.add("hMassDPlus", "Dplus Pair candidates;inv. mass DPlus only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDMinus", "Dplus Pair candidates;inv. mass DMinus only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDPlusMcRecSig", "DPlus signal candidates - MC reco;inv. mass DPlus only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDPlusMcRecBkg", "DPlus background candidates - MC reco;inv. mass DPlus only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDMinusMcRecSig", "DMinus signal candidates - MC reco;inv. mass DMinus only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDMinusMcRecBkg", "DMinus background candidates - MC reco;inv. mass DMinus only (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  template <typename T, typename U>
  void analyseMultiplicity(const T& collision, const U& tracks)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (track.eta() < etaCutMin || track.eta() > etaCutMax) {
          continue;
        }
        if (std::abs(track.dcaXY()) > dcaXYCut || std::abs(track.dcaZ()) > dcaZCut) {
          continue;
        }
        nTracks++;
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);
  }

  template <typename T>
  bool kinematicCuts(const T& candidate)
  {
    // check decay channel flag for candidate
    bool cuts = true;
    if (particlePdgCode == pdg::Code::kD0) {
      if (!(candidate.hfflag() & 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        cuts = false;
      }
    } else if (particlePdgCode == pdg::Code::kDPlus) {
      if (!(candidate.hfflag() & 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi)) {
        cuts = false;
      }
    }
    if (yCandMax >= 0. && std::abs(candidate.y(RecoDecay::getMassPDG(particlePdgCode))) > yCandMax) {
      cuts = false;
    }
    if (ptCandMin >= 0. && candidate.pt() < ptCandMin) {
      cuts = false;
    }
    return cuts;
  }

  template <typename T>
  bool kinematicCutsGen(const T& particle)
  {
    bool cuts = true;
    // check if the particle is D or Dbar (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
    if (std::abs(particle.pdgCode()) != particlePdgCode) {
      cuts = false;
    }
    if (yCandMax >= 0. && std::abs(particle.y()) > yCandMax) {
      cuts = false;
    }
    if (ptCandMin >= 0. && particle.pt() < ptCandMin) {
      cuts = false;
    }
    return cuts;
  }

  template <typename T>
  void fillInfoHists(const T& candidate, bool const& isReco)
  {
    if (isReco) {
      registry.fill(HIST("hPtCandMcRec"), candidate.pt());
      registry.fill(HIST("hPtProng0McRec"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1McRec"), candidate.ptProng1());
      registry.fill(HIST("hEtaMcRec"), candidate.eta());
      registry.fill(HIST("hPhiMcRec"), candidate.phi());
      registry.fill(HIST("hYMcRec"), candidate.y(RecoDecay::getMassPDG(particlePdgCode)));
    } else {
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hEta"), candidate.eta());
      registry.fill(HIST("hPhi"), candidate.phi());
      registry.fill(HIST("hY"), candidate.y(RecoDecay::getMassPDG(particlePdgCode)));
    }
  }

  // candidateType assignations. Possible values of candidateType:
  // |1: Sig D0              |10: Ref D0              |21: Sig D0 + Ref D0bar     |102: Bkg D0 + Sig D0bar      |201: Sig D0 + Bkg D0bar
  // |2: Sig D0bar           |12: Ref D0 + Sig D0bar  |30: Ref D0 + Ref D0bar     |120: Bkg D0 + Ref D0bar      |210: Ref D0 + Bkg D0bar
  // |3: Sig D0 + Sig D0bar  |20: Ref D0bar           |100: Bkg D0                |200: Bkg D0bar               |200: Bkg D0 + Bkg D0bar

  template <typename T, bool isRecoMc>
  uint assignCandidateTypeD0(const T& candidate, bool const& fillMassHists)
  {
    uint candidateType = 0;
    if (candidate.isSelD0() >= selectionFlagD0) {
      if constexpr (isRecoMc) {
        if (candidate.flagMcMatchRec() == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) { // also matched as D0
          if (fillMassHists) {
            registry.fill(HIST("hMassD0McRecSig"), invMassD0ToPiK(candidate), candidate.pt());
          }
          candidateType += 1;
        } else if (candidate.flagMcMatchRec() == -(1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          if (fillMassHists) {
            registry.fill(HIST("hMassD0McRecRefl"), invMassD0ToPiK(candidate), candidate.pt());
          }
          candidateType += 10;
        } else {
          if (fillMassHists) {
            registry.fill(HIST("hMassD0McRecBkg"), invMassD0ToPiK(candidate), candidate.pt());
          }
          candidateType += 100;
        }
      } else {
        if (fillMassHists) {
          registry.fill(HIST("hMass"), invMassD0ToPiK(candidate), candidate.pt());
          registry.fill(HIST("hMassD0"), invMassD0ToPiK(candidate), candidate.pt()); // only reco as D0
        }
        candidateType += 1;
      }
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) { // only reco as D0bar
      if constexpr (isRecoMc) {
        if (candidate.flagMcMatchRec() == -(1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) { // also matched as D0bar
          if (fillMassHists) {
            registry.fill(HIST("hMassD0barMcRecSig"), invMassD0barToKPi(candidate), candidate.pt());
          }
          candidateType += 2;
        } else if (candidate.flagMcMatchRec() == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
          if (fillMassHists) {
            registry.fill(HIST("hMassD0barMcRecRefl"), invMassD0barToKPi(candidate), candidate.pt());
          }
          candidateType += 20;
        } else {
          if (fillMassHists) {
            registry.fill(HIST("hMassD0barMcRecBkg"), invMassD0barToKPi(candidate), candidate.pt());
          }
          candidateType += 200;
        }
      } else {
        if (fillMassHists) {
          registry.fill(HIST("hMass"), invMassD0barToKPi(candidate), candidate.pt());
          registry.fill(HIST("hMassD0bar"), invMassD0barToKPi(candidate), candidate.pt());
        }
        candidateType += 2;
      }
    }
    if (verbose) {
      // candidateType = 0 -> no candidate, candidateType = 1 -> D0, candidateType1 = 2 -> D0bar, candidateType = 3 -> both selections passed
      std::cout << "candidateType: " << candidateType << std::endl;
      std::cout << " ----------------------------- " << std::endl;
    }
    return candidateType;
  }

  template <typename T, bool isRecoMc>
  uint assignCandidateTypeDPlus(const T& candidate, double const& mDPlus, int const& particleSign, bool const& fillMassHists)
  {
    uint candidateType = 0;
    if constexpr (isRecoMc) {
      if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi) {
        if (particleSign == 1) { // reco and matched as Dplus
          if (fillMassHists) {
            registry.fill(HIST("hMassDplusMcRecSig"), mDPlus, candidate.pt());
          }
          candidateType += 1;
        } else { // reco and matched as Dminus
          if (fillMassHists) {
            registry.fill(HIST("hMassDminusMcRecSig"), mDPlus, candidate.pt());
          }
          candidateType += 2;
        }
      } else {
        // fill invariant mass plots from Dplus/Dminus background candidates
        if (particleSign == 1) { // reco as Dplus
          if (fillMassHists) {
            registry.fill(HIST("hMassDplusMcRecBkg"), mDPlus, candidate.pt());
          }
          candidateType += 100;
        } else { // matched as Dminus
          if (fillMassHists) {
            registry.fill(HIST("hMassDminusMcRecBkg"), mDPlus, candidate.pt());
          }
          candidateType += 200;
        }
      }
    } else {
      if (particleSign == 1) {
        if (fillMassHists) {
          registry.fill(HIST("hMass"), mDPlus, candidate.pt());
          registry.fill(HIST("hMassDPlus"), mDPlus, candidate.pt());
        }
        candidateType += 1;
      } else {
        if (fillMassHists) {
          registry.fill(HIST("hMass"), mDPlus, candidate.pt());
          registry.fill(HIST("hMassDMinus"), mDPlus, candidate.pt());
        }
        candidateType += 2;
      }
    }
    if (verbose) {
      std::cout << "candidateType: " << candidateType << std::endl;
      std::cout << " ----------------------------- " << std::endl;
    }
    return candidateType;
  }

  template <typename T>
  uint assignCandidateTypeGen(const T& candidate)
  {
    if (candidate.pdgCode() == particlePdgCode) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
      return 1;
    } else if (candidate.pdgCode() == -particlePdgCode) { // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
      return 2;
    } else {
      return 0;
    }
  }

  template <typename T>
  void analyseMcGen(const T& particlesMc)
  {
    registry.fill(HIST("hMcEvtCount"), 0);
    for (const auto& particle1 : particlesMc) {
      // check if the particle is D0, D0bar, DPlus or DMinus (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) == pdg::Code::kD0 || std::abs(particle1.pdgCode()) == pdg::Code::kDPlus) {
        continue;
      }
      double yD = RecoDecay::y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (!kinematicCutsGen(particle1)) {
        continue;
      }

      registry.fill(HIST("hPtCandMcGen"), particle1.pt());
      registry.fill(HIST("hEtaMcGen"), particle1.eta());
      registry.fill(HIST("hPhiMcGen"), particle1.phi());
      registry.fill(HIST("hYMcGen"), yD);

      auto candidateType1 = assignCandidateTypeGen(particle1); // Candidate sign attribution.

      // check if it's prompt or non-prompt
      int8_t originGen1 = particle1.originMcGen();
      registry.fill(HIST("hOriginMcGen"), originGen1);
      // check if it's MC matched
      int8_t matchedGen1 = particle1.flagMcMatchGen();
      registry.fill(HIST("hMatchedMcGen"), matchedGen1);

      for (const auto& particle2 : particlesMc) {
        // Candidate sign attribution.
        auto candidateType2 = assignCandidateTypeGen(particle2);
        if (!kinematicCutsGen(particle2)) {
          continue;
        }

        // check if it's prompt or non-prompt
        int8_t originGen2 = particle2.originMcGen();
        // check if it's MC matched
        int8_t matchedGen2 = particle2.flagMcMatchGen();

        // If both particles are D0's, fill D0Pair table
        if (std::abs(particle1.pdgCode()) == pdg::Code::kD0 && std::abs(particle2.pdgCode()) == pdg::Code::kD0) {
          entryD0Pair(getDeltaPhi(particle2.phi(), particle1.phi()),
                      particle2.eta() - particle1.eta(),
                      particle1.pt(),
                      particle2.pt(),
                      particle1.y(),
                      particle2.y(),
                      RecoDecay::getMassPDG(pdg::Code::kD0),
                      RecoDecay::getMassPDG(pdg::Code::kD0),
                      candidateType1,
                      candidateType2,
                      0);
          entryD0PairRecoInfo(originGen1,
                              originGen2,
                              matchedGen1,
                              matchedGen2,
                              0);
        // If both particles are DPlus', fill DPlusPair table
        } else if (std::abs(particle1.pdgCode()) == pdg::Code::kDPlus && std::abs(particle2.pdgCode()) == pdg::Code::kDPlus) {
          entryDPlusPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                         particle2.eta() - particle1.eta(),
                         particle1.pt(),
                         particle2.pt(),
                         particle1.y(),
                         particle2.y(),
                         RecoDecay::getMassPDG(pdg::Code::kDPlus),
                         RecoDecay::getMassPDG(pdg::Code::kDPlus),
                         candidateType1,
                         candidateType2,
                         0);
          entryDPlusPairRecoInfo(originGen1,
                                 originGen2,
                                 matchedGen1,
                                 matchedGen2,
                                 0);
        }
      } // end inner loop
    }   // end outer loop
  }

  /// D0(bar)-D0(bar) correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processDataD0(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
  {
    if (particlePdgCode != pdg::Code::kD0) {
      LOGF(fatal, "Wrong PDG code. Change it to 421 (D0)."); // Assure that we have the right PDG code
    }
    analyseMultiplicity(collision, tracks);
    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedD0CandidatesGrouped) {
      if (!kinematicCuts(candidate1)) {
        continue;
      }
      fillInfoHists(candidate1, false);

      auto candidateType1 = assignCandidateTypeD0<decltype(candidate1), false>(candidate1, true); // Candidate type attribution.
      registry.fill(HIST("hSelectionStatus"), candidateType1);

      for (const auto& candidate2 : selectedD0CandidatesGrouped) {
        if (!kinematicCuts(candidate2)) {
          continue;
        }
        // excluding trigger self-correlations (possible in case of both mass hypotheses accepted)
        if (candidate1.mRowIndex == candidate2.mRowIndex) {
          continue;
        }
        auto candidateType2 = assignCandidateTypeD0<decltype(candidate2), false>(candidate2, false); // Candidate type attribution

        // Register the pair only if we have two candidates that have passed all cuts
        if (candidateType1 != 0 && candidateType2 != 0) {
          entryD0Pair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                      candidate2.eta() - candidate1.eta(),
                      candidate1.pt(),
                      candidate2.pt(),
                      yD0(candidate1),
                      yD0(candidate2),
                      invMassD0ToPiK(candidate1),
                      invMassD0barToKPi(candidate2),
                      candidateType1,
                      candidateType2,
                      1);
        } // end if
      }   // end inner loop (Cand2)
    }     // end outer loop (Cand1)
  }
  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processDataD0, "Process data D0", true);

  /// D0(bar)-D0(bar) correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRecD0(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const&)
  {
    if (particlePdgCode != pdg::Code::kD0) {
      LOGF(fatal, "Wrong PDG code. Change it to 421 (D0)."); // Assure that we have the right PDG code
    }
    analyseMultiplicity(collision, tracks);
    auto selectedD0CandidatesGroupedMc = selectedD0candidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedD0CandidatesGroupedMc) {
      if (!kinematicCuts(candidate1)) {
        continue;
      }
      if (std::abs(candidate1.flagMcMatchRec()) == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
        fillInfoHists(candidate1, true);
      }

      auto candidateType1 = assignCandidateTypeD0<decltype(candidate1), true>(candidate1, true); // Candidate type attribution
      registry.fill(HIST("hSelectionStatusMcRec"), candidateType1);
      int8_t origin1 = 0, matchedRec1 = 0;
      if (candidateType1 < 100) { // if our event is not bkg
        // check if it's prompt or non-prompt
        origin1 = candidate1.originMcRec();
        registry.fill(HIST("hOriginMcRec"), origin1);
        // check if it's MC matched
        matchedRec1 = candidate1.flagMcMatchRec();
        registry.fill(HIST("hMatchedMcRec"), matchedRec1);
      }

      for (const auto& candidate2 : selectedD0CandidatesGroupedMc) {
        if (!kinematicCuts(candidate2)) {
          continue;
        }
        // Excluding trigger self-correlations (possible in case of both mass hypotheses accepted)
        if (candidate1.mRowIndex == candidate2.mRowIndex) {
          continue;
        }

        auto candidateType2 = assignCandidateTypeD0<decltype(candidate2), true>(candidate2, false); // Candidate type attribution
        int8_t origin2 = 0, matchedRec2 = 0;
        if (candidateType1 < 100) { // if our event is not bkg
          // check if it's prompt or non-prompt
          origin2 = candidate2.originMcRec();
          // check if it's MC matched
          matchedRec2 = candidate2.flagMcMatchRec();
        }

        // Register the pair only if we have two candidates that have passed all cuts
        if (candidateType1 != 0 && candidateType2 != 0) {
          entryD0Pair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                      candidate2.eta() - candidate1.eta(),
                      candidate1.pt(),
                      candidate2.pt(),
                      yD0(candidate1),
                      yD0(candidate2),
                      invMassD0ToPiK(candidate1),
                      invMassD0barToKPi(candidate2),
                      candidateType1,
                      candidateType2,
                      0);
          entryD0PairRecoInfo(origin1,
                              origin2,
                              matchedRec1,
                              matchedRec2,
                              1);
        }
      } // end inner loop
    }   // end outer loop
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcRecD0, "Process D0 Mc Reco mode", false);

  /// D0(bar)-D0(bar) correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGenD0(aod::McCollision const&, McParticlesPlus2Prong const& particlesMc)
  {
    if (particlePdgCode != pdg::Code::kD0) {
      LOGF(fatal, "Wrong PDG code. Change it to 421 (D0)."); // Assure that we have the right PDG code
    }
    analyseMcGen(particlesMc);
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcGenD0, "Process D0 Mc Gen mode", false);

  /// Dplus(minus)-Dplus(minus) correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processDataDPlus(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const&, aod::BigTracks const&)
  {
    if (particlePdgCode != pdg::Code::kDPlus) {
      LOGF(fatal, "Wrong PDG code. Change it to 431 (DPlus)."); // Assure that we have the right PDG code
    }
    analyseMultiplicity(collision, tracks);
    auto selectedDPlusCandidatesGrouped = selectedDPlusCandidates->sliceByCached(o2::aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedDPlusCandidatesGrouped) {
      if (!kinematicCuts(candidate1)) {
        continue;
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.prong1_as<aod::BigTracks>();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }
      double mDPlus = invMassDplusToPiKPi(candidate1);
      uint candidateType1 = assignCandidateTypeDPlus<decltype(candidate1), false>(candidate1, mDPlus, outerParticleSign, true);
      fillInfoHists(candidate1, false);

      registry.fill(HIST("hSelectionStatus"), candidateType1);

      for (const auto& candidate2 : selectedDPlusCandidatesGrouped) {
        if (!kinematicCuts(candidate2)) {
          continue;
        }
        if (candidate1.mRowIndex == candidate2.mRowIndex) {
          continue;
        }

        int innerParticleSign = 1; // Dplus
        auto innerSecondTrack = candidate2.prong1_as<aod::BigTracks>();
        if (innerSecondTrack.sign() == 1) {
          innerParticleSign = -1; // Dminus (second daughter track is positive)
        }
        uint candidateType2 = assignCandidateTypeDPlus<decltype(candidate2), false>(candidate2, mDPlus, innerParticleSign, false);

        // Register the pair only if we have two candidates that have passed all cuts
        if (candidateType1 != 0 && candidateType2 != 0) {
          entryDPlusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                         candidate2.eta() - candidate1.eta(),
                         candidate1.pt(),
                         candidate2.pt(),
                         yDplus(candidate1),
                         yDplus(candidate2),
                         invMassDplusToPiKPi(candidate1),
                         invMassDplusToPiKPi(candidate2),
                         candidateType1,
                         candidateType2,
                         1);
        } // end if
      }   // end inner loop (cand2)
    }     // end outer loop (cand1)
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processDataDPlus, "Process Data DPlus", false);

  /// Dplus(minus)-Dplus(minus) correlation pair builder - for MC reco-level analysis (candidates matched to true signal only, but also the various bkg sources are studied)
  void processMcRecDPlus(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec> const&, aod::BigTracks const&)
  {
    if (particlePdgCode != pdg::Code::kDPlus) {
      LOGF(fatal, "Wrong PDG code. Change it to 431 (DPlus)."); // Assure that we have the right PDG code
    }
    analyseMultiplicity(collision, tracks);
    auto selectedDPlusCandidatesGroupedMc = selectedDPlusCandidatesMc->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    for (const auto& candidate1 : selectedDPlusCandidatesGroupedMc) {
      if (!kinematicCuts(candidate1)) {
        continue;
      }

      int outerParticleSign = 1; // Dplus
      auto outerSecondTrack = candidate1.prong1_as<aod::BigTracks>();
      if (outerSecondTrack.sign() == 1) {
        outerParticleSign = -1; // Dminus (second daughter track is positive)
      }
      double mDPlus = invMassDplusToPiKPi(candidate1);
      uint candidateType1 = assignCandidateTypeDPlus<decltype(candidate1), true>(candidate1, mDPlus, outerParticleSign, true);

      registry.fill(HIST("hSelectionStatusMcRec"), candidateType1);
      int8_t origin1 = 0, matchedRec1 = 0;
      // if our event is not bkg
      if (candidateType1 < 100) {
        // check if it's prompt or non-prompt
        origin1 = candidate1.originMcRec();
        registry.fill(HIST("hOriginMcRec"), origin1);
        // check if it's MC matched
        matchedRec1 = candidate1.flagMcMatchRec();
        registry.fill(HIST("hMatchedMcRec"), matchedRec1);
      }

      for (const auto& candidate2 : selectedDPlusCandidatesGroupedMc) {
        if (!kinematicCuts(candidate2)) {
          continue;
        }
        if (candidate1.mRowIndex == candidate2.mRowIndex) {
          continue;
        }

        int innerParticleSign = 1; // Dplus
        auto innerSecondTrack = candidate2.prong1_as<aod::BigTracks>();
        if (innerSecondTrack.sign() == 1) {
          innerParticleSign = -1; // Dminus (second daughter track is positive)
        }

        uint candidateType2 = assignCandidateTypeDPlus<decltype(candidate2), true>(candidate2, mDPlus, innerParticleSign, false);
        int8_t origin2 = 0, matchedRec2 = 0;
        // if our event is not bkg
        if (candidateType2 < 100) {
          // check if it's prompt or non-prompt
          origin2 = candidate1.originMcRec();
          // check if it's MC matched
          matchedRec2 = candidate1.flagMcMatchRec();
        }

        // Register the pair only if we have two candidates that have passed all cuts
        if (candidateType1 != 0 && candidateType2 != 0) {
          entryDPlusPair(getDeltaPhi(candidate2.phi(), candidate1.phi()),
                         candidate2.eta() - candidate1.eta(),
                         candidate1.pt(),
                         candidate2.pt(),
                         yDplus(candidate1),
                         yDplus(candidate2),
                         invMassDplusToPiKPi(candidate1),
                         invMassDplusToPiKPi(candidate2),
                         candidateType1,
                         candidateType2,
                         0);
          entryDPlusPairRecoInfo(origin1,
                                 origin2,
                                 matchedRec1,
                                 matchedRec2,
                                 1);
        } // end if
      }   // end inner loop (cand2)
    }     // end outer loop (cand1)
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcRecDPlus, "Process DPlus Mc Reco", false);

  /// Dplus(minus)-Dplus(minus) correlation pair builder - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGenDPlus(aod::McCollision const&, McParticlesPlus3Prong const& particlesMc)
  {
    if (particlePdgCode != pdg::Code::kDPlus) {
      LOGF(fatal, "Wrong PDG code. Change it to 431 (DPlus)."); // Assure that we have the right PDG code
    }
    analyseMcGen(particlesMc);
  }

  PROCESS_SWITCH(HfCorrelatorDMesonPairs, processMcGenDPlus, "Process DPlus Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorDMesonPairs>(cfgc)};
}
