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

/// \file jetHFAngularity.cxx
/// \brief Measure heavy flavour angularity(D0 & Lambda_c)

/// \author Rajdeep Nandi
/// \author Preeti Dhankher

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/JetCandidateUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>

#include <fmt/format.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr int OriginMcPrompt = 1;    // RecoDecay::OriginType::Prompt
constexpr int OriginMcNonPrompt = 2; // RecoDecay::OriginType::NonPrompt
constexpr int NOrigins = 2;          // prompt, non-prompt (index 0, 1 of the per-origin objects)

constexpr int CandSelBitD0 = 0;
constexpr int CandSelBitD0bar = 1;

constexpr int ChanD0ToPiK = static_cast<int>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
constexpr int ChanD0ToPiKPi0 = static_cast<int>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiKPi0);
constexpr int ChanD0ToKK = static_cast<int>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToKK);

/// correlated-background channels: 2 = K pi pi0, 3 = pi pi, 4 = pi pi pi0, 5 = K K
inline bool isCorrBkgChannel(int absFlag)
{
  return absFlag >= ChanD0ToPiKPi0 && absFlag <= ChanD0ToKK;
}

struct McCandClass {
  bool isTrueSignal = false; // D0 -> K pi, selected with the correct K/pi mass hypothesis
  bool isReflection = false; // D0 -> K pi, selected with the swapped mass hypothesis
  bool isCorrBkg = false;    // other D0 decay channels, either hypothesis
  int channel = 0;           // |flagMcMatchRec|
};

template <typename CAND>
McCandClass classifyMcCandidate(CAND const& cand)
{
  McCandClass c;
  const int flag = cand.flagMcMatchRec();
  c.channel = std::abs(flag);
  if (c.channel == ChanD0ToPiK) {
    int hypothesis = 0; // +1: row selected as D0, -1: row selected as D0bar
    const int sel = cand.candidateSelFlag();
    if ((sel & BIT(CandSelBitD0)) != 0) {
      hypothesis = 1;
    } else if ((sel & BIT(CandSelBitD0bar)) != 0) {
      hypothesis = -1;
    }
    // reflection: true species differs from the hypothesis used for this row
    c.isReflection = (hypothesis != 0) && ((flag > 0) != (hypothesis > 0));
    c.isTrueSignal = !c.isReflection;
  } else if (isCorrBkgChannel(c.channel)) {
    c.isCorrBkg = true;
  }
  return c;
}

/// Returns 1 (input = pseudo-data sample) for a fraction 'fraction` of MC collisions, else 0 ("correction" sample).
inline int closureSample(int64_t mcCollisionIndex, float fraction)
{
  if (fraction <= 0.f) {
    return 0;
  }
  uint64_t h = static_cast<uint64_t>(mcCollisionIndex) * 0x9E3779B97F4A7C15ULL;
  h ^= (h >> 29);
  h *= 0xBF58476D1CE4E5B9ULL;
  h ^= (h >> 32);
  return (static_cast<float>(h % 100000ULL) < fraction * 100000.f) ? 1 : 0;
}
} // namespace

namespace o2::aod
{
// Jet-side columns (shared by all tables)
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(JetLambda, jetLambda, float);
DECLARE_SOA_COLUMN(JetMass, jetMass, float);

// D0-candidate columns
DECLARE_SOA_COLUMN(HfZParallel, hfZParallel, float);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfSelFlag, hfSelFlag, int8_t);

// ML scores (one entry per BDT class), -1 if the model has fewer classes
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);

// MC-truth columns
DECLARE_SOA_COLUMN(HfFlagMcMatch, hfFlagMcMatch, int8_t);
DECLARE_SOA_COLUMN(HfOriginMc, hfOriginMc, int8_t);

DECLARE_SOA_COLUMN(JetMatchGeoStatus, jetMatchGeoStatus, int8_t);
DECLARE_SOA_COLUMN(JetMatchCandStatus, jetMatchCandStatus, int8_t);
DECLARE_SOA_COLUMN(JetMatchClean, jetMatchClean, int8_t);
DECLARE_SOA_COLUMN(JetMatchedPt, jetMatchedPt, float);
DECLARE_SOA_COLUMN(JetMatchedEta, jetMatchedEta, float);
DECLARE_SOA_COLUMN(JetMatchedPhi, jetMatchedPhi, float);
DECLARE_SOA_COLUMN(JetMatchedDR, jetMatchedDR, float);
DECLARE_SOA_COLUMN(JetMatchedD0Pt, jetMatchedD0Pt, float);
DECLARE_SOA_COLUMN(JetMatchedLambda, jetMatchedLambda, float);
DECLARE_SOA_COLUMN(ClosureSample, closureSample, int8_t);

DECLARE_SOA_TABLE(D0JetObjTable, "AOD", "D0JETOBJTABLE",
                  JetHfDist, JetPt, JetEta, JetPhi, JetNConst,
                  JetLambda, JetMass,
                  HfZParallel, HfPt, HfEta, HfPhi, HfMass, HfY, HfSelFlag,
                  HfMlScore0, HfMlScore1, HfMlScore2);

DECLARE_SOA_TABLE(D0JetMCDObjTable, "AOD", "D0JETMCDOBJ",
                  JetHfDist, JetPt, JetEta, JetPhi, JetNConst,
                  JetLambda, JetMass,
                  HfZParallel, HfPt, HfEta, HfPhi, HfMass, HfY, HfSelFlag,
                  HfMlScore0, HfMlScore1, HfMlScore2,
                  HfFlagMcMatch, HfOriginMc);

DECLARE_SOA_TABLE(D0JetMCPObjTable, "AOD", "D0JETMCPOBJ",
                  JetHfDist, JetPt, JetEta, JetPhi, JetNConst,
                  JetLambda, JetMass,
                  HfZParallel, HfPt, HfEta, HfPhi, HfY,
                  HfFlagMcMatch, HfOriginMc);

DECLARE_SOA_TABLE(D0JetMCDObjMatchedTable, "AOD", "D0JMCDOBJMATCH",
                  JetHfDist, JetPt, JetEta, JetPhi, JetNConst,
                  JetLambda, JetMass,
                  HfZParallel, HfPt, HfEta, HfPhi, HfMass, HfY, HfSelFlag,
                  HfMlScore0, HfMlScore1, HfMlScore2,
                  HfFlagMcMatch, HfOriginMc,
                  JetMatchGeoStatus, JetMatchCandStatus, JetMatchClean,
                  JetMatchedPt, JetMatchedEta, JetMatchedPhi, JetMatchedDR, JetMatchedD0Pt,
                  JetMatchedLambda,
                  ClosureSample);
} // namespace o2::aod

consteval float getValFromBin(int bin)
{
  return static_cast<float>(bin) - 0.5f;
}

/// Collision-counter bins (data)
enum BinCollCntr {
  AllCollisions = 1,
  SelectedCollisions = 2,
  CollisionsWithD0Jet = 3,
  RejectedCollisions = 4
};

/// Jet-counter bins
enum BinJetCntr {
  D0JetsInFiducialRange = 1,
  D0JetsWithSelectedCandidate = 2
};

/// MC collision-counter bins
enum BinMCColCntr {
  AllMCCollisions = 1,
  SelectedMCCollisions = 2,
  AssociatedRecoCollisions = 3,
  SelectedAssociatedRecoCollisions = 4
};

/// MC jet-counter bins
enum BinMCJetCntr {
  DetectorLevelJetInMCCollision = 1,
  ParticleLevelJetInMCCollision = 2,
  DetectorLevelJetWithMatchedCandidate = 3,
  ParticleLevelJetWithMatchedCandidate = 4
};

struct JetHFAngularityTask {

  using D0CandidatesData = aod::CandidatesD0Data;
  using D0DataJets = soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>;

  using D0CandidatesMCD = aod::CandidatesD0MCD;
  using D0MCDJets = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>;
  using D0MCDJetsMatched = soa::Join<aod::D0ChargedMCDetectorLevelJets,
                                     aod::D0ChargedMCDetectorLevelJetConstituents,
                                     aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;

  using D0CandidatesMCP = aod::CandidatesD0MCP;
  using D0MCPJets = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>;
  using D0MCPJetsMatched = soa::Join<aod::D0ChargedMCParticleLevelJets,
                                     aod::D0ChargedMCParticleLevelJetConstituents,
                                     aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;

  // ---------------------------------------------------------------- configurables
  // Event selection
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range (cm), same as jet finder"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "Event selection string, same as jet finder"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "Skip MB-gap events (same as jet finder skipMBGapEvents)"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "Apply RCT selection (same as jet finder applyRCTSelections)"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "Track selection string (track QA only)"};

  // Jets
  Configurable<float> jetPtMin{"jetPtMin", 5.0f, "Minimum jet pT (GeV/c). Use 0 in MC so the response covers migrations from below the analysis range"};
  Configurable<float> jetR{"jetR", 0.4f, "Jet resolution parameter R"};

  // Observables
  Configurable<float> alpha{"alpha", 1.f, "Angularity exponent alpha of DeltaR/R (one value per run; e.g. 1, 1.5, 2, 3)"};
  Configurable<float> kappa{"kappa", 1.f, "Angularity exponent kappa of the momentum fraction z (one value per run; 1 = IRC safe)"};
  Configurable<float> lambdaMax{"lambdaMax", 1.f, "Upper edge of the lambda axis (1 for kappa = 1; adjust e.g. 0.8/0.6/0.5 for alpha = 1.5/2/3)"};
  Configurable<int> nBinsLambda{"nBinsLambda", 100, "Number of lambda bins in the data / pseudo-data sparses"};
  Configurable<int> nBinsLambdaResp{"nBinsLambdaResp", 50, "Number of lambda bins in response / efficiency objects, mass-template and legacy sparses"};

  // Optional BDT cut on top of the HF selector (for the BDT systematic). Off by default.
  Configurable<bool> applyMlCut{"applyMlCut", false, "Apply an additional pT-dependent cut on ML score 0 (background)"};
  Configurable<std::vector<float>> mlCutPtEdges{"mlCutPtEdges", std::vector<float>{1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36., 50.}, "pT_D0 bin edges of the ML cut"};
  Configurable<std::vector<float>> mlCutBkgMax{"mlCutBkgMax", std::vector<float>{0.02, 0.02, 0.02, 0.05, 0.06, 0.08, 0.08, 0.1, 0.1, 0.2, 0.25, 0.3, 0.4}, "Max ML background score per pT_D0 bin (n_edges - 1 values)"};

  // MC
  Configurable<bool> rejectReflections{"rejectReflections", true, "Exclude reflections from response matrices and efficiency numerators"};
  Configurable<bool> useMcWeights{"useMcWeights", false, "Weight MC entries with mcCollision.weight() (pT-hard weighted productions)"};
  Configurable<float> closureInputFraction{"closureInputFraction", 0.2f, "Fraction of MC collisions flagged as closure 'input' (sample axis = 1)"};
  Configurable<int> matchingMode{"matchingMode", 0, "MCD-MCP jet matching for response/efficiency: 0 = clean (geo AND HF point to the same MCP jet, analysis note), 1 = HF (candidate) matching only"};

  Configurable<bool> fillLegacyOutputs{"fillLegacyOutputs", true, "Also fill the histograms of the previous task version"};
  Configurable<bool> fillTables{"fillTables", true, "Fill the per-candidate output tables"};

  // Axes
  ConfigurableAxis binsMassD0{"binsMassD0", {200, 1.6f, 2.2f}, "m_{K#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis binsPtD0{"binsPtD0", {100, 0.f, 50.f}, "#it{p}_{T,D^{0}} (GeV/#it{c})"};
  ConfigurableAxis binsZParallel{"binsZParallel", {100, 0.f, 1.f}, "z_{||}"};
  ConfigurableAxis binsDeltaR{"binsDeltaR", {50, 0.f, 0.5f}, "#DeltaR"};
  ConfigurableAxis binsPtJetData{"binsPtJetData", {60, 0.f, 60.f}, "#it{p}_{T,jet} (GeV/#it{c})"};
  ConfigurableAxis binsPtJetMC{"binsPtJetMC", {60, 0.f, 60.f}, "#it{p}_{T,jet} (GeV/#it{c})"};
  ConfigurableAxis binsPtD0Resp{"binsPtD0Resp", {VARIABLE_WIDTH, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36., 50.}, "#it{p}_{T,D^{0}} (GeV/#it{c})"};
  ConfigurableAxis binsPtJetResp{"binsPtJetResp", {VARIABLE_WIDTH, 0., 3., 5., 7., 10., 15., 20., 30., 50., 80.}, "#it{p}_{T,jet} (GeV/#it{c})"};

  Produces<aod::D0JetObjTable> objJetTable;
  Produces<aod::D0JetMCDObjTable> objJetMCDTable;
  Produces<aod::D0JetMCPObjTable> objJetMCPTable;
  Produces<aod::D0JetMCDObjMatchedTable> objJetMCDMatchedTable;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  float massD0Pdg = -1.f;
  int jetRBin = 40;

  std::shared_ptr<THnSparse> hMassPtLambda;   // data
  std::shared_ptr<THnSparse> hMassPtLambdaMC; // MC pseudo-data (closure)
  std::array<std::shared_ptr<THnSparse>, NOrigins> hResponse;
  std::array<std::shared_ptr<THnSparse>, NOrigins> hGen;
  std::array<std::shared_ptr<THnSparse>, NOrigins> hGenMatched;
  std::array<std::shared_ptr<THnSparse>, NOrigins> hDetTrue;
  std::array<std::shared_ptr<THnSparse>, NOrigins> hDetMatched;
  std::shared_ptr<THnSparse> hGenD0Jet; // MCP  particle level

  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<D0MCDJetsMatched> d0MCDJetsPerCollisionPreslice = aod::jet::collisionId;
  Preslice<D0MCPJetsMatched> d0MCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;

  Filter jetCutsPt = aod::jet::pt > jetPtMin;
  Filter jetCutsR = aod::jet::r == nround(jetR.node() * 100.0f);

  HistogramRegistry registry{"registry"};

  // ================================================================ helpers

  struct JetObservables {
    float lambda = -1.f;
    float mass = 0.f;
    int nConst = 0;
  };

  struct D0InJet {
    float zParallel = -1.f;
    float deltaR = -1.f;
  };

  template <typename JET, typename TRACKS, typename CANDIDATES>
  static float computeLambda(JET const& jet, TRACKS const& tracks, CANDIDATES const& candidates, float expAlpha, float expKappa)
  {
    if (jet.pt() <= 0.f) {
      return -1.f;
    }
    double sum = 0.;
    for (auto const& trk : tracks) {
      sum += std::pow(trk.pt(), expKappa) * std::pow(jetutilities::deltaR(jet, trk), expAlpha);
    }
    for (auto const& cand : candidates) {
      sum += std::pow(cand.pt(), expKappa) * std::pow(jetutilities::deltaR(jet, cand), expAlpha);
    }
    const double jetRadius = jet.r() / 100.; // stored as R x 100
    const double denom = std::pow(jet.pt(), expKappa) * std::pow(jetRadius, expAlpha);
    return denom > 0. ? static_cast<float>(sum / denom) : -1.f;
  }

  template <typename TRACKS, typename CANDIDATES>
  float computeJetMass(TRACKS const& tracks, CANDIDATES const& candidates) const
  {
    std::array<double, 3> momTotal{0., 0., 0.};
    double energyTot = 0.;
    for (auto const& trk : tracks) {
      momTotal[0] += trk.px();
      momTotal[1] += trk.py();
      momTotal[2] += trk.pz();
      energyTot += RecoDecay::p(trk.px(), trk.py(), trk.pz());
    }
    for (auto const& cand : candidates) {
      momTotal[0] += cand.px();
      momTotal[1] += cand.py();
      momTotal[2] += cand.pz();
      energyTot += RecoDecay::e(cand.px(), cand.py(), cand.pz(), static_cast<double>(massD0Pdg));
    }
    const double mass2 = RecoDecay::m2(momTotal, energyTot);
    return (mass2 > 0.) ? static_cast<float>(std::sqrt(mass2)) : 0.f;
  }

  template <typename JET, typename TRACKS, typename CANDIDATES>
  JetObservables jetObservables(JET const& jet, TRACKS const& tracks, CANDIDATES const& candidates) const
  {
    JetObservables obs;
    obs.lambda = computeLambda(jet, tracks, candidates, alpha.value, kappa.value);
    obs.mass = computeJetMass(tracks, candidates);
    obs.nConst = static_cast<int>(tracks.size()) + static_cast<int>(candidates.size());
    return obs;
  }

  template <typename JET, typename CAND>
  static D0InJet d0JetKinematics(JET const& jet, CAND const& cand)
  {
    const std::array<double, 3> jetMom{jet.px(), jet.py(), jet.pz()};
    const std::array<double, 3> d0Mom{cand.px(), cand.py(), cand.pz()};
    D0InJet k;
    k.zParallel = static_cast<float>(RecoDecay::dotProd(d0Mom, jetMom) / RecoDecay::mag2(jetMom));
    k.deltaR = jetutilities::deltaR(jet, cand);
    return k;
  }

  template <typename CAND>
  static float mlScore(CAND const& cand, size_t i)
  {
    const auto scores = cand.mlScores();
    return (scores.size() > i) ? scores[i] : -1.f;
  }

  /// optional pT-dependent cut on the background score (BDT systematic)
  template <typename CAND>
  bool passMlCut(CAND const& cand) const
  {
    if (!applyMlCut.value) {
      return true;
    }
    const auto& edges = mlCutPtEdges.value;
    const auto& cuts = mlCutBkgMax.value;
    const float pt = cand.pt();
    for (size_t i = 0; i + 1 < edges.size() && i < cuts.size(); ++i) {
      if (pt >= edges[i] && pt < edges[i + 1]) {
        const float bkg = mlScore(cand, 0);
        return bkg >= 0.f && bkg < cuts[i];
      }
    }
    return false;
  }

  template <typename COLL>
  bool selectRecoCollision(COLL const& collision) const
  {
    if (std::abs(collision.posZ()) >= vertexZCut.value) {
      return false;
    }
    return jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents.value, applyRCTSelections.value);
  }

  template <typename MCCOLL>
  bool selectMcCollision(MCCOLL const& mcCollision) const
  {
    if (std::abs(mcCollision.posZ()) >= vertexZCut.value) {
      return false;
    }
    if (skipMBGapEvents.value && mcCollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return false;
    }
    return true;
  }

  template <typename MCCOLL>
  float mcWeight(MCCOLL const& mcCollision) const
  {
    return useMcWeights.value ? mcCollision.weight() : 1.f;
  }

  static int originIndex(int origin)
  {
    if (origin == OriginMcPrompt) {
      return 0;
    }
    if (origin == OriginMcNonPrompt) {
      return 1;
    }
    return -1;
  }

  template <typename CAND>
  void fillMcQa(CAND const& cand)
  {
    const int flag = cand.flagMcMatchRec();
    const int sel = cand.candidateSelFlag();
    registry.fill(HIST("h_d0_flag_mcd"), static_cast<float>(flag));
    registry.fill(HIST("h_d0_selflag"), static_cast<float>(sel));
    // mass of true D0 -> K pi rows per (sign of flag, selection bit):
    //   1: flag>0 & D0 hyp (signal)       2: flag>0 & D0bar hyp (reflection)
    //   3: flag<0 & D0 hyp (reflection)   4: flag<0 & D0bar hyp (signal)
    //   5: both bits                      6: no bit
    if (std::abs(flag) == ChanD0ToPiK) {
      const bool bit0 = (sel & BIT(CandSelBitD0)) != 0;
      const bool bit1 = (sel & BIT(CandSelBitD0bar)) != 0;
      float category = 6.f;
      if (bit0 && bit1) {
        category = 5.f;
      } else if (bit0) {
        category = (flag > 0) ? 1.f : 3.f;
      } else if (bit1) {
        category = (flag > 0) ? 2.f : 4.f;
      }
      registry.fill(HIST("h_d0_mass_hypothesis_check"), cand.m(), category);
    }
  }

  // ================================================================ init

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections.value);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections.value);
    massD0Pdg = jetcandidateutilities::getTablePDGMass<D0CandidatesMCP>();
    jetRBin = static_cast<int>(std::lround(static_cast<float>(jetR) * 100.0f));
    if (alpha.value <= 0.f || kappa.value <= 0.f || lambdaMax.value <= 0.f) {
      LOGF(fatal, "alpha, kappa and lambdaMax must be > 0, got %g, %g, %g", static_cast<float>(alpha), static_cast<float>(kappa), static_cast<float>(lambdaMax));
    }
    LOGF(info, "angularity lambda_alpha^kappa with alpha = %g, kappa = %g, lambda axis [0, %g]", static_cast<float>(alpha), static_cast<float>(kappa), static_cast<float>(lambdaMax));
    if (matchingMode != 0 && matchingMode != 1) {
      LOGF(fatal, "matchingMode must be 0 (clean: geo AND HF) or 1 (HF only), got %d", static_cast<int>(matchingMode));
    }
    if (applyMlCut && mlCutBkgMax.value.size() + 1 != mlCutPtEdges.value.size()) {
      LOGF(fatal, "mlCutBkgMax needs exactly (mlCutPtEdges - 1) entries");
    }

    const AxisSpec axMass{binsMassD0, "m_{K#pi} (GeV/#it{c}^{2})"};
    const AxisSpec axPtD0{binsPtD0, "#it{p}_{T,D^{0}} (GeV/#it{c})"};
    const AxisSpec axZ{binsZParallel, "z_{||}"};
    const AxisSpec axDR{binsDeltaR, "#DeltaR"};

    const std::string alphaStr = fmt::format("{:g}", static_cast<float>(alpha));
    const std::string kappaStr = fmt::format("{:g}", static_cast<float>(kappa));
    const std::string lamLabel = fmt::format("#lambda_{{{}}}^{{{}}}", alphaStr, kappaStr);
    const AxisSpec axLam{nBinsLambda.value, 0.f, lambdaMax.value, lamLabel}; // data / pseudo-data sparses
    const AxisSpec axLamT{nBinsLambdaResp.value, 0.f, lambdaMax.value, lamLabel};
    const AxisSpec axPtJetData{binsPtJetData, "#it{p}_{T,jet} (GeV/#it{c})"};
    const AxisSpec axPtJetDet{binsPtJetMC, "#it{p}_{T,jet}^{det} (GeV/#it{c})"};
    const AxisSpec axPtJetPart{binsPtJetMC, "#it{p}_{T,jet}^{part} (GeV/#it{c})"};
    const AxisSpec axPtD0Det{binsPtD0Resp, "#it{p}_{T,D^{0}}^{det} (GeV/#it{c})"};
    const AxisSpec axPtD0Part{binsPtD0Resp, "#it{p}_{T,D^{0}}^{part} (GeV/#it{c})"};
    const AxisSpec axPtJetRespDet{binsPtJetResp, "#it{p}_{T,jet}^{det} (GeV/#it{c})"};
    const AxisSpec axPtJetRespPart{binsPtJetResp, "#it{p}_{T,jet}^{part} (GeV/#it{c})"};
    const AxisSpec axLamDet{nBinsLambdaResp.value, 0.f, lambdaMax.value, fmt::format("#lambda_{{{}}}^{{{},det}}", alphaStr, kappaStr)};
    const AxisSpec axLamPart{nBinsLambdaResp.value, 0.f, lambdaMax.value, fmt::format("#lambda_{{{}}}^{{{},part}}", alphaStr, kappaStr)};
    const AxisSpec axOrigin{3, -0.5, 2.5, "Origin(D^{0})"};
    const AxisSpec axMatch{2, -0.5, 1.5, "Matching status"};
    const AxisSpec axRefl{2, -0.5, 1.5, "IsReflection"};
    const AxisSpec axSample{2, -0.5, 1.5, "closure sample (0 = correction, 1 = input)"};
    const AxisSpec axChannelCorr{4, 1.5, 5.5, "|flagMcMatch| (correlated bkg channel)"};
    const AxisSpec axChannelGen{5, 0.5, 5.5, "|flagMcMatchGen| (D^{0} decay channel)"};
    const AxisSpec axHypCat{6, 0.5, 6.5, "category"};

    // ---------------- QA (all levels)
    registry.add("h_collisions", "Event status;;entries", {HistType::kTH1F, {{4, 0., 4.}}});
    registry.add("h_track_pt", ";#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_track_eta", ";#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_track_phi", ";#varphi_{track};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    registry.add("h_jet_counter", ";;counts", {HistType::kTH1F, {{2, 0., 2.}}});
    registry.add("h_jet_pt", ";#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_jet_eta", ";#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi", ";#phi_{jet};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    registry.add("h_d0_mass", ";m_{K#pi} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axMass}});
    registry.add("h_d0_pt", ";#it{p}_{T,D^{0}} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}});
    registry.add("h_d0_eta", ";#eta_{D^{0}};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_d0_phi", ";#phi_{D^{0}};entries", {HistType::kTH1F, {{100, -1., 7.}}});
    registry.add("h_d0_y", ";y_{D^{0}};entries", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("h_d0_jet_projection", ";z_{||}^{D^{0},jet};entries", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_d0_jet_distance", ";#DeltaR_{D^{0},jet};entries", {HistType::kTH1F, {{100, 0., 0.5}}});
    registry.add("h_d0_jet_lambda", (";" + lamLabel + ";entries").c_str(), {HistType::kTH1F, {{200, 0., lambdaMax.value}}});
    registry.add("h_d0_jet_mass", ";m_{jet} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}});
    registry.add("h_jet_nconst", ";N_{const};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}});
    registry.add("h_d0_mlscore0", ";ML score 0 (bkg);entries", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_d0_mlscore1", ";ML score 1 (prompt);entries", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_d0_mlscore2", ";ML score 2 (non-prompt);entries", {HistType::kTH1F, {{100, 0., 1.}}});

    auto hColl = registry.get<TH1>(HIST("h_collisions"));
    hColl->GetXaxis()->SetBinLabel(AllCollisions, "all");
    hColl->GetXaxis()->SetBinLabel(SelectedCollisions, "selected (ev. sel + RCT + |z|)");
    hColl->GetXaxis()->SetBinLabel(CollisionsWithD0Jet, "has D0 jet");
    hColl->GetXaxis()->SetBinLabel(RejectedCollisions, "rejected");
    auto hJet = registry.get<TH1>(HIST("h_jet_counter"));
    hJet->GetXaxis()->SetBinLabel(D0JetsInFiducialRange, "D0 jets (pT, R)");
    hJet->GetXaxis()->SetBinLabel(D0JetsWithSelectedCandidate, "D0 jets with selected cand.");

    const std::string title = fmt::format("alpha = {}, kappa = {}", alphaStr, kappaStr);
    if (doprocessDataChargedSubstructure) {
      hMassPtLambda = registry.add<THnSparse>("ang/hMassPtLambda_a0", title.c_str(), HistType::kTHnSparseF,
                                              {axMass, axPtD0, axPtJetData, axLam});
    }
    if (doprocessMCDChargedSubstructureMatched) {
      hMassPtLambdaMC = registry.add<THnSparse>("ang/hMassPtLambdaMC_a0", title.c_str(), HistType::kTHnSparseF,
                                                {axMass, axPtD0, axPtJetData, axLam, axSample});
      const std::array<std::string, NOrigins> originTag{"pr", "np"};
      for (int o = 0; o < NOrigins; ++o) {
        const std::string sfx = "a0_" + originTag[o];

        hResponse[o] = registry.add<THnSparse>(("ang/hResponse_" + sfx).c_str(), title.c_str(), HistType::kTHnSparseF,
                                               {axPtJetRespDet, axLamDet, axPtD0Det, axPtJetRespPart, axLamPart, axPtD0Part, axSample});
        hGen[o] = registry.add<THnSparse>(("ang/hGen_" + sfx).c_str(), title.c_str(), HistType::kTHnSparseF,
                                          {axPtJetRespPart, axLamPart, axPtD0Part, axSample});
        hGenMatched[o] = registry.add<THnSparse>(("ang/hGenMatched_" + sfx).c_str(), title.c_str(), HistType::kTHnSparseF,
                                                 {axPtJetRespPart, axLamPart, axPtD0Part, axSample});
        hDetTrue[o] = registry.add<THnSparse>(("ang/hDetTrue_" + sfx).c_str(), title.c_str(), HistType::kTHnSparseF,
                                              {axPtJetRespDet, axLamDet, axPtD0Det, axSample});
        hDetMatched[o] = registry.add<THnSparse>(("ang/hDetMatched_" + sfx).c_str(), title.c_str(), HistType::kTHnSparseF,
                                                 {axPtJetRespDet, axLamDet, axPtD0Det, axSample});
      }
    }
    if (doprocessMCPChargedSubstructure) {
      hGenD0Jet = registry.add<THnSparse>("ang/hGenD0Jet_a0", title.c_str(), HistType::kTHnSparseF,
                                          {axPtD0Part, axPtJetRespPart, axLamPart, axOrigin});
    }

    if (doprocessDataChargedSubstructure && fillLegacyOutputs) {
      registry.add("hSparse_d0", "hSparse_d0", HistType::kTHnSparseF, {axMass, axPtD0, axPtJetData, axZ, axDR, axLamT});
    }

    // ---------------- MC detector level (templates, QA)
    if (doprocessMCDChargedSubstructure || doprocessMCDChargedSubstructureMatched) {
      registry.add("h_d0_mass_mcd", ";m_{K#pi}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axMass}});
      registry.add("h_d0_mass_mcd_signal", ";m_{K#pi}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axMass}});
      registry.add("h_d0_mass_mcd_reflection", ";m_{K#pi}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axMass}});
      registry.add("h_d0_origin_mcd", ";origin;entries", {HistType::kTH1F, {{3, 0., 3.}}});
      registry.add("h_d0_flag_mcd", ";flagMcMatchRec;entries", HistType::kTH1F, {{11, -5.5, 5.5}});
      registry.add("h_d0_selflag", ";candidateSelFlag;entries", HistType::kTH1F, {{5, -0.5, 4.5}});
      registry.add("h_d0_mass_hypothesis_check", ";m_{K#pi} (GeV/#it{c}^{2});category", HistType::kTH2F, {axMass, axHypCat});
      registry.add("hSparse_d0_mcd", "hSparse_d0_mcd", HistType::kTHnSparseF, {axMass, axPtD0, axPtJetDet, axZ, axDR, axLamT});
      registry.add("hSparse_d0_mcd_calib", "hSparse_d0_mcd_calib", HistType::kTHnSparseF, {axMass, axPtD0, axPtJetDet, axZ, axOrigin, axLamT});
      registry.add("hSparse_d0_mcd_reflection", "hSparse_d0_mcd_reflection", HistType::kTHnSparseF, {axMass, axPtD0, axPtJetDet, axZ, axOrigin, axLamT});
      registry.add("hSparse_d0_mcd_corrbkg", "hSparse_d0_mcd_corrbkg", HistType::kTHnSparseF, {axMass, axPtD0, axPtJetDet, axZ, axChannelCorr, axOrigin, axLamT});
      registry.add("h_d0_mass_mcd_corrbkg", ";m_{K#pi}^{det} (GeV/#it{c}^{2});channel", HistType::kTH2F, {axMass, axChannelCorr});
      auto hCorr = registry.get<THnSparse>(HIST("hSparse_d0_mcd_corrbkg"));
      hCorr->GetAxis(4)->SetBinLabel(1, "K#pi#pi^{0}");
      hCorr->GetAxis(4)->SetBinLabel(2, "#pi#pi");
      hCorr->GetAxis(4)->SetBinLabel(3, "#pi#pi#pi^{0}");
      hCorr->GetAxis(4)->SetBinLabel(4, "KK");
    }

    // ---------------- MC particle level / matched
    if (doprocessMCPChargedSubstructure || doprocessMCDChargedSubstructureMatched) {
      registry.add("MC_CollisionCounter", "N_{collisions};", {HistType::kTH1F, {{4, 0., 4.}}});
      registry.add("MC_JetCounter", "N_{jet};", {HistType::kTH1F, {{4, 0., 4.}}});
      registry.add("h_d0_flag_mcp", ";flagMcMatchGen;entries", HistType::kTH1F, {{11, -5.5, 5.5}});
      registry.add("hSparse_d0_mcp_channel", "hSparse_d0_mcp_channel", HistType::kTHnSparseF, {axPtD0, axPtJetPart, axZ, axChannelGen, axOrigin, axLamT});
      auto c = registry.get<TH1>(HIST("MC_CollisionCounter"));
      c->GetXaxis()->SetBinLabel(AllMCCollisions, "All MC collisions");
      c->GetXaxis()->SetBinLabel(SelectedMCCollisions, "Selected MC collisions");
      c->GetXaxis()->SetBinLabel(AssociatedRecoCollisions, "Associated reco collisions");
      c->GetXaxis()->SetBinLabel(SelectedAssociatedRecoCollisions, "Selected associated reco collisions");
      auto j = registry.get<TH1>(HIST("MC_JetCounter"));
      j->GetXaxis()->SetBinLabel(DetectorLevelJetInMCCollision, "Detector-level D0 jets");
      j->GetXaxis()->SetBinLabel(ParticleLevelJetInMCCollision, "Particle-level D0 jets");
      j->GetXaxis()->SetBinLabel(DetectorLevelJetWithMatchedCandidate, "Det. jets matched (matchingMode)");
      j->GetXaxis()->SetBinLabel(ParticleLevelJetWithMatchedCandidate, "Part. jets matched (matchingMode)");
      auto hChan = registry.get<THnSparse>(HIST("hSparse_d0_mcp_channel"));
      hChan->GetAxis(3)->SetBinLabel(1, "K#pi (signal)");
      hChan->GetAxis(3)->SetBinLabel(2, "K#pi#pi^{0}");
      hChan->GetAxis(3)->SetBinLabel(3, "#pi#pi");
      hChan->GetAxis(3)->SetBinLabel(4, "#pi#pi#pi^{0}");
      hChan->GetAxis(3)->SetBinLabel(5, "KK");
    }

    if (doprocessMCDChargedSubstructureMatched) {
      registry.add("h_jet_matching_geo_status_mcd", ";status (0=unmatched,1=matched);entries", {HistType::kTH1F, {{2, 0., 2.}}});
      registry.add("h_jet_matching_cand_status_mcd", ";status (0=unmatched,1=matched);entries", {HistType::kTH1F, {{2, 0., 2.}}});
      registry.add("h_jet_matching_clean_mcd", ";status (0=not clean,1=clean);entries", {HistType::kTH1F, {{2, 0., 2.}}});
      registry.add("h_jet_matching_dr_mcd", ";#DeltaR(jet^{det},jet^{part});entries", {HistType::kTH1F, {{100, 0., 0.5}}});
      registry.add("h_jet_matching_ngeo_mcd", ";number of geo-matched jets;entries", {HistType::kTH1F, {{6, 0., 6.}}});
      registry.add("h_jet_matching_ncand_mcd", ";number of cand-matched jets;entries", {HistType::kTH1F, {{6, 0., 6.}}});
      registry.add("h_closure_sample", ";closure sample;MC collisions", {HistType::kTH1F, {{2, -0.5, 1.5}}});

      registry.add("h_eff_run2_num_pr", ";#it{p}_{T,D^{0}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}});
      registry.add("h_eff_run2_num_np", ";#it{p}_{T,D^{0}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}});
      registry.add("h_eff_run2_den_pr", ";#it{p}_{T,D^{0}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}});
      registry.add("h_eff_run2_den_np", ";#it{p}_{T,D^{0}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}});

      if (fillLegacyOutputs) {
        registry.add("h_jet_pt_response_matrix", "", HistType::kTH2F, {axPtJetDet, axPtJetPart});
        registry.add("h_jet_pt_response_matrix_pr", "prompt", HistType::kTH2F, {axPtJetDet, axPtJetPart});
        registry.add("h_jet_pt_response_matrix_np", "non-prompt", HistType::kTH2F, {axPtJetDet, axPtJetPart});
        registry.add("h_response_d0pt_pr", "", HistType::kTHnSparseF, {axPtJetDet, axPtD0Det, axPtJetPart, axPtD0Part});
        registry.add("h_response_d0pt_np", "", HistType::kTHnSparseF, {axPtJetDet, axPtD0Det, axPtJetPart, axPtD0Part});
        registry.add("h_response_angularity_pr", "", HistType::kTHnSparseF, {axPtJetDet, axLamDet, axPtJetPart, axLamPart});
        registry.add("h_response_angularity_np", "", HistType::kTHnSparseF, {axPtJetDet, axLamDet, axPtJetPart, axLamPart});
        registry.add("hSparse_d0_mcd_eff", "hSparse_d0_mcd_eff", HistType::kTHnSparseF, {axMass, axPtD0, axPtJetDet, axZ, axOrigin, axMatch, axLamT, axRefl});
        registry.add("hSparse_d0_mcp_eff", "hSparse_d0_mcp_eff", HistType::kTHnSparseF, {axPtD0, axPtJetPart, axZ, axOrigin, axMatch, axLamT});
      }
    }

    if (doprocessMCPChargedSubstructure) {
      registry.add("h_jet_pt_mcp", ";#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta_mcp", ";#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_d0_pt_mcp", ";#it{p}_{T,D^{0}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}});
      registry.add("h_d0_y_mcp", ";y_{D^{0}}^{part};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("h_d0_origin_mcp", ";origin;entries", {HistType::kTH1F, {{3, 0., 3.}}});
    }
  }

  // ================================================================ data

  void processCollisions(aod::JetCollision const& collision, aod::JetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), getValFromBin(AllCollisions));
    if (!selectRecoCollision(collision)) {
      registry.fill(HIST("h_collisions"), getValFromBin(RejectedCollisions));
      return;
    }
    registry.fill(HIST("h_collisions"), getValFromBin(SelectedCollisions)); // N_ev for the feed-down normalisation
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }
  PROCESS_SWITCH(JetHFAngularityTask, processCollisions, "DATA: event counter (N_ev) and track QA", true);

  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      soa::Filtered<D0DataJets> const& jets,
                                      D0CandidatesData const&,
                                      aod::JetTracks const&)
  {
    if (!selectRecoCollision(collision)) {
      return;
    }
    bool collisionHasD0Jet = false;
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_counter"), getValFromBin(D0JetsInFiducialRange));
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());

      auto jetTracks = jet.template tracks_as<aod::JetTracks>();
      auto jetCandidates = jet.template candidates_as<D0CandidatesData>();
      const auto obs = jetObservables(jet, jetTracks, jetCandidates);

      bool hasSelectedD0 = false;
      for (const auto& cand : jetCandidates) {
        if (!passMlCut(cand)) {
          continue;
        }
        hasSelectedD0 = true;
        const auto kin = d0JetKinematics(jet, cand);

        registry.fill(HIST("h_d0_mass"), cand.m());
        registry.fill(HIST("h_d0_pt"), cand.pt());
        registry.fill(HIST("h_d0_eta"), cand.eta());
        registry.fill(HIST("h_d0_phi"), cand.phi());
        registry.fill(HIST("h_d0_y"), cand.y());
        registry.fill(HIST("h_d0_jet_projection"), kin.zParallel);
        registry.fill(HIST("h_d0_jet_distance"), kin.deltaR);
        registry.fill(HIST("h_jet_nconst"), obs.nConst);
        registry.fill(HIST("h_d0_mlscore0"), mlScore(cand, 0));
        registry.fill(HIST("h_d0_mlscore1"), mlScore(cand, 1));
        registry.fill(HIST("h_d0_mlscore2"), mlScore(cand, 2));

        hMassPtLambda->Fill(cand.m(), cand.pt(), jet.pt(), obs.lambda);
        if (fillLegacyOutputs) {
          registry.fill(HIST("hSparse_d0"), cand.m(), cand.pt(), jet.pt(), kin.zParallel, kin.deltaR, obs.lambda);
        }
        if (fillTables) {
          objJetTable(kin.deltaR, jet.pt(), jet.eta(), jet.phi(), obs.nConst,
                      obs.lambda, obs.mass,
                      kin.zParallel, cand.pt(), cand.eta(), cand.phi(), cand.m(), cand.y(),
                      static_cast<int8_t>(cand.candidateSelFlag()),
                      mlScore(cand, 0), mlScore(cand, 1), mlScore(cand, 2));
        }
      }
      if (hasSelectedD0) {
        registry.fill(HIST("h_jet_counter"), getValFromBin(D0JetsWithSelectedCandidate));
        registry.fill(HIST("h_d0_jet_lambda"), obs.lambda);
        registry.fill(HIST("h_d0_jet_mass"), obs.mass);
        collisionHasD0Jet = true;
      }
    }
    if (collisionHasD0Jet) {
      registry.fill(HIST("h_collisions"), getValFromBin(CollisionsWithD0Jet));
    }
  }
  PROCESS_SWITCH(JetHFAngularityTask, processDataChargedSubstructure, "DATA: D0-tagged charged-jet angularities", true);

  // ================================================================ MC detector level (QA only)

  void processMCDChargedSubstructure(aod::JetCollision const& collision,
                                     soa::Filtered<D0MCDJets> const& mcdjets,
                                     D0CandidatesMCD const&,
                                     aod::JetTracks const&)
  {
    if (!selectRecoCollision(collision)) {
      return;
    }
    for (const auto& jet : mcdjets) {
      auto jetTracks = jet.template tracks_as<aod::JetTracks>();
      auto jetCandidates = jet.template candidates_as<D0CandidatesMCD>();
      const auto obs = jetObservables(jet, jetTracks, jetCandidates);
      for (const auto& cand : jetCandidates) {
        if (!passMlCut(cand)) {
          continue;
        }
        const auto kin = d0JetKinematics(jet, cand);
        fillMcQa(cand);
        registry.fill(HIST("h_d0_mass_mcd"), cand.m());
        registry.fill(HIST("h_d0_origin_mcd"), static_cast<float>(cand.originMcRec()));
        registry.fill(HIST("hSparse_d0_mcd"), cand.m(), cand.pt(), jet.pt(), kin.zParallel, kin.deltaR, obs.lambda);
        if (fillTables) {
          objJetMCDTable(kin.deltaR, jet.pt(), jet.eta(), jet.phi(), obs.nConst,
                         obs.lambda, obs.mass,
                         kin.zParallel, cand.pt(), cand.eta(), cand.phi(), cand.m(), cand.y(),
                         static_cast<int8_t>(cand.candidateSelFlag()),
                         mlScore(cand, 0), mlScore(cand, 1), mlScore(cand, 2),
                         static_cast<int8_t>(cand.flagMcMatchRec()), static_cast<int8_t>(cand.originMcRec()));
        }
      }
    }
  }
  PROCESS_SWITCH(JetHFAngularityTask, processMCDChargedSubstructure, "MC DET LEVEL QA (no efficiency objects)", false);

  // ================================================================ MC particle level

  void processMCPChargedSubstructure(aod::JetMcCollision const& mcCollision,
                                     soa::Filtered<D0MCPJets> const& mcpjets,
                                     D0CandidatesMCP const&,
                                     aod::JetParticles const&)
  {
    registry.fill(HIST("MC_CollisionCounter"), getValFromBin(AllMCCollisions));
    if (!selectMcCollision(mcCollision)) {
      return;
    }
    registry.fill(HIST("MC_CollisionCounter"), getValFromBin(SelectedMCCollisions)); // N_ev,gen for sigma/N_ev scaling
    const float w = mcWeight(mcCollision);

    for (const auto& jet : mcpjets) {
      registry.fill(HIST("MC_JetCounter"), getValFromBin(ParticleLevelJetInMCCollision));
      registry.fill(HIST("h_jet_pt_mcp"), jet.pt(), w);
      registry.fill(HIST("h_jet_eta_mcp"), jet.eta(), w);

      auto jetParticles = jet.template tracks_as<aod::JetParticles>();
      auto jetCandidates = jet.template candidates_as<D0CandidatesMCP>();
      const auto obs = jetObservables(jet, jetParticles, jetCandidates);

      for (const auto& d0 : jetCandidates) {
        const auto kin = d0JetKinematics(jet, d0);
        const int flag = d0.flagMcMatchGen();
        const int absFlag = std::abs(flag);
        const int origin = d0.originMcGen();
        registry.fill(HIST("h_d0_pt_mcp"), d0.pt(), w);
        registry.fill(HIST("h_d0_y_mcp"), d0.y(), w);
        registry.fill(HIST("h_d0_origin_mcp"), static_cast<float>(origin), w);
        registry.fill(HIST("h_d0_flag_mcp"), static_cast<float>(flag), w);
        if (absFlag >= ChanD0ToPiK && absFlag <= ChanD0ToKK) {
          registry.fill(HIST("hSparse_d0_mcp_channel"), d0.pt(), jet.pt(), kin.zParallel, static_cast<float>(absFlag),
                        static_cast<float>(origin), obs.lambda, w);
        }

        if (isCorrBkgChannel(absFlag)) {
          continue;
        }
        const std::array<double, 4> x{d0.pt(), jet.pt(), obs.lambda, static_cast<double>(origin)};
        hGenD0Jet->Fill(x.data(), w);
        if (fillTables) {
          objJetMCPTable(kin.deltaR, jet.pt(), jet.eta(), jet.phi(), obs.nConst,
                         obs.lambda, obs.mass,
                         kin.zParallel, d0.pt(), d0.eta(), d0.phi(), d0.y(),
                         static_cast<int8_t>(flag), static_cast<int8_t>(origin));
        }
      }
    }
  }
  PROCESS_SWITCH(JetHFAngularityTask, processMCPChargedSubstructure, "MC PART LEVEL (POWHEG feed-down input)", false);

  // ================================================================ MC matched: templates, response, efficiency

  struct MatchInfo {
    bool geo = false;
    bool cand = false;
    bool clean = false;
    bool used = false;
    int nGeo = 0;
    int nCand = 0;
    float pt = -1.f;
    float eta = -999.f;
    float phi = -999.f;
    float dR = -1.f;
    float d0Pt = -1.f;
    float lambda = -1.f;
  };

  template <typename MCDJET>
  MatchInfo findMatch(MCDJET const& mcdjet)
  {
    MatchInfo m;
    int geoIndex = -1;
    if (mcdjet.has_matchedJetGeo()) {
      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<D0MCPJetsMatched>()) {
        if (m.nGeo++ == 0) {
          geoIndex = mcpjet.globalIndex();
          m.geo = true;
        }
      }
    }
    int candIndex = -1;
    if (mcdjet.has_matchedJetCand()) {
      for (const auto& mcpjet : mcdjet.template matchedJetCand_as<D0MCPJetsMatched>()) {
        if (m.nCand++ == 0) {
          candIndex = mcpjet.globalIndex();
          m.cand = true;
        }
      }
    }
    m.clean = m.geo && m.cand && (geoIndex == candIndex);
    m.used = (matchingMode == 1) ? m.cand : m.clean;
    if (!m.used) {
      return m;
    }
    for (const auto& mcpjet : mcdjet.template matchedJetCand_as<D0MCPJetsMatched>()) {
      m.pt = mcpjet.pt();
      m.eta = mcpjet.eta();
      m.phi = mcpjet.phi();
      m.dR = jetutilities::deltaR(mcdjet, mcpjet);
      auto parts = mcpjet.template tracks_as<aod::JetParticles>();
      auto cands = mcpjet.template candidates_as<D0CandidatesMCP>();
      const auto obs = jetObservables(mcpjet, parts, cands);
      m.lambda = obs.lambda;
      if (cands.size() > 0) {
        m.d0Pt = mcpjet.template candidates_first_as<D0CandidatesMCP>().pt();
      }
      break; // first candidate match (== first geo match)
    }
    return m;
  }

  void processMCDChargedSubstructureMatched(aod::JetMcCollisions const& mccollisions,
                                            aod::JetCollisionsMCD const& collisions,
                                            D0MCDJetsMatched const& mcdjets,
                                            D0MCPJetsMatched const& mcpjets,
                                            D0CandidatesMCD const&,
                                            D0CandidatesMCP const&,
                                            aod::JetTracks const&,
                                            aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {
      registry.fill(HIST("MC_CollisionCounter"), getValFromBin(AllMCCollisions));
      if (!selectMcCollision(mccollision)) {
        continue;
      }
      registry.fill(HIST("MC_CollisionCounter"), getValFromBin(SelectedMCCollisions));
      const float w = mcWeight(mccollision);
      const int sample = closureSample(mccollision.globalIndex(), closureInputFraction);
      registry.fill(HIST("h_closure_sample"), sample);

      // ------------------------------------------------------------ detector level
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {
        registry.fill(HIST("MC_CollisionCounter"), getValFromBin(AssociatedRecoCollisions));
        if (!selectRecoCollision(collision)) {
          continue;
        }
        registry.fill(HIST("MC_CollisionCounter"), getValFromBin(SelectedAssociatedRecoCollisions));

        const auto jetsThisColl = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : jetsThisColl) {
          if (mcdjet.pt() <= jetPtMin || mcdjet.r() != jetRBin) {
            continue;
          }
          registry.fill(HIST("MC_JetCounter"), getValFromBin(DetectorLevelJetInMCCollision));
          registry.fill(HIST("h_jet_eta"), mcdjet.eta());

          auto jetTracks = mcdjet.template tracks_as<aod::JetTracks>();
          auto jetCandidates = mcdjet.template candidates_as<D0CandidatesMCD>();
          const auto obs = jetObservables(mcdjet, jetTracks, jetCandidates);
          const auto match = findMatch(mcdjet);

          registry.fill(HIST("h_jet_matching_ngeo_mcd"), match.nGeo);
          registry.fill(HIST("h_jet_matching_ncand_mcd"), match.nCand);
          registry.fill(HIST("h_jet_matching_geo_status_mcd"), match.geo ? 1.5f : 0.5f);
          registry.fill(HIST("h_jet_matching_cand_status_mcd"), match.cand ? 1.5f : 0.5f);
          registry.fill(HIST("h_jet_matching_clean_mcd"), match.clean ? 1.5f : 0.5f);
          if (match.used) {
            registry.fill(HIST("h_jet_matching_dr_mcd"), match.dR);
            registry.fill(HIST("MC_JetCounter"), getValFromBin(DetectorLevelJetWithMatchedCandidate));
          }

          for (const auto& cand : jetCandidates) {
            if (!passMlCut(cand)) {
              continue;
            }
            const auto kin = d0JetKinematics(mcdjet, cand);
            const int origin = cand.originMcRec();
            const auto mcClass = classifyMcCandidate(cand);
            fillMcQa(cand);
            registry.fill(HIST("h_d0_y"), cand.y());
            registry.fill(HIST("h_d0_mass_mcd"), cand.m(), w);
            registry.fill(HIST("h_d0_origin_mcd"), static_cast<float>(origin), w);

            // pseudo-data for the full-chain closure: every candidate, like in data
            const std::array<double, 5> xMc{cand.m(), cand.pt(), mcdjet.pt(), obs.lambda, static_cast<double>(sample)};
            hMassPtLambdaMC->Fill(xMc.data(), w);

            // mass templates
            if (mcClass.isReflection) {
              registry.fill(HIST("h_d0_mass_mcd_reflection"), cand.m(), w);
              registry.fill(HIST("hSparse_d0_mcd_reflection"), cand.m(), cand.pt(), mcdjet.pt(), kin.zParallel, origin, obs.lambda, w);
            } else if (mcClass.isTrueSignal) {
              registry.fill(HIST("h_d0_mass_mcd_signal"), cand.m(), w);
              registry.fill(HIST("hSparse_d0_mcd_calib"), cand.m(), cand.pt(), mcdjet.pt(), kin.zParallel, origin, obs.lambda, w);
            } else if (mcClass.isCorrBkg) {
              registry.fill(HIST("h_d0_mass_mcd_corrbkg"), cand.m(), static_cast<float>(mcClass.channel), w);
              registry.fill(HIST("hSparse_d0_mcd_corrbkg"), cand.m(), cand.pt(), mcdjet.pt(), kin.zParallel,
                            static_cast<float>(mcClass.channel), static_cast<float>(origin), obs.lambda, w);
            }
            if (fillLegacyOutputs && !mcClass.isCorrBkg) {
              registry.fill(HIST("hSparse_d0_mcd_eff"), cand.m(), cand.pt(), mcdjet.pt(), kin.zParallel, origin,
                            match.clean ? 1 : 0, obs.lambda, mcClass.isReflection ? 1 : 0, w);
            }

            // true D0 -> K pi with origin: efficiency numerators and response
            const bool skip = !(mcClass.isTrueSignal || mcClass.isReflection) || (rejectReflections && mcClass.isReflection);
            const int o = originIndex(origin);
            if (skip || o < 0) {
              continue;
            }
            if (o == 0) {
              registry.fill(HIST("h_eff_run2_num_pr"), cand.pt(), w);
            } else {
              registry.fill(HIST("h_eff_run2_num_np"), cand.pt(), w);
            }
            const std::array<double, 4> xd{mcdjet.pt(), obs.lambda, cand.pt(), static_cast<double>(sample)};
            hDetTrue[o]->Fill(xd.data(), w);
            if (!match.used) {
              continue;
            }
            const std::array<double, 4> xg{match.pt, match.lambda, match.d0Pt, static_cast<double>(sample)};
            const std::array<double, 7> xr{mcdjet.pt(), obs.lambda, cand.pt(), match.pt, match.lambda, match.d0Pt, static_cast<double>(sample)};
            hDetMatched[o]->Fill(xd.data(), w);
            hGenMatched[o]->Fill(xg.data(), w);
            hResponse[o]->Fill(xr.data(), w);
            if (fillLegacyOutputs) {
              registry.fill(HIST("h_jet_pt_response_matrix"), mcdjet.pt(), match.pt, w);
              if (o == 0) {
                registry.fill(HIST("h_jet_pt_response_matrix_pr"), mcdjet.pt(), match.pt, w);
                registry.fill(HIST("h_response_angularity_pr"), mcdjet.pt(), obs.lambda, match.pt, match.lambda, w);
                registry.fill(HIST("h_response_d0pt_pr"), mcdjet.pt(), cand.pt(), match.pt, match.d0Pt, w);
              } else {
                registry.fill(HIST("h_jet_pt_response_matrix_np"), mcdjet.pt(), match.pt, w);
                registry.fill(HIST("h_response_angularity_np"), mcdjet.pt(), obs.lambda, match.pt, match.lambda, w);
                registry.fill(HIST("h_response_d0pt_np"), mcdjet.pt(), cand.pt(), match.pt, match.d0Pt, w);
              }
            }
            if (fillTables) {
              objJetMCDMatchedTable(kin.deltaR, mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), obs.nConst,
                                    obs.lambda, obs.mass,
                                    kin.zParallel, cand.pt(), cand.eta(), cand.phi(), cand.m(), cand.y(),
                                    static_cast<int8_t>(cand.candidateSelFlag()),
                                    mlScore(cand, 0), mlScore(cand, 1), mlScore(cand, 2),
                                    static_cast<int8_t>(cand.flagMcMatchRec()), static_cast<int8_t>(origin),
                                    static_cast<int8_t>(match.geo), static_cast<int8_t>(match.cand), static_cast<int8_t>(match.clean),
                                    match.pt, match.eta, match.phi, match.dR, match.d0Pt,
                                    match.lambda,
                                    static_cast<int8_t>(sample));
            }
          } // candidates
        } // MCD jets
      } // reco collisions

      // ------------------------------------------------------------ particle level
      const auto mcpJetsThisColl = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : mcpJetsThisColl) {
        if (mcpjet.pt() <= jetPtMin || mcpjet.r() != jetRBin) {
          continue;
        }
        registry.fill(HIST("MC_JetCounter"), getValFromBin(ParticleLevelJetInMCCollision));
        auto parts = mcpjet.template tracks_as<aod::JetParticles>();
        auto cands = mcpjet.template candidates_as<D0CandidatesMCP>();
        if (cands.size() == 0) {
          continue;
        }
        const auto obs = jetObservables(mcpjet, parts, cands);
        const auto d0 = mcpjet.template candidates_first_as<D0CandidatesMCP>();
        const auto kin = d0JetKinematics(mcpjet, d0);
        const int flag = d0.flagMcMatchGen();
        const int absFlag = std::abs(flag);
        const int origin = d0.originMcGen();

        registry.fill(HIST("h_d0_flag_mcp"), static_cast<float>(flag), w);
        if (absFlag >= ChanD0ToPiK && absFlag <= ChanD0ToKK) {
          registry.fill(HIST("hSparse_d0_mcp_channel"), d0.pt(), mcpjet.pt(), kin.zParallel, static_cast<float>(absFlag),
                        static_cast<float>(origin), obs.lambda, w);
        }
        if (absFlag != ChanD0ToPiK) {
          continue; // only D0 -> K pi enters the efficiency denominators
        }

        bool cleanMatched = false;
        if (mcpjet.has_matchedJetGeo() && mcpjet.has_matchedJetCand()) {
          int geoTarget = -1;
          for (const auto& j : mcpjet.template matchedJetGeo_as<D0MCDJetsMatched>()) {
            geoTarget = j.globalIndex();
            break;
          }
          int candTarget = -1;
          for (const auto& j : mcpjet.template matchedJetCand_as<D0MCDJetsMatched>()) {
            candTarget = j.globalIndex();
            break;
          }
          cleanMatched = (geoTarget != -1) && (geoTarget == candTarget);
        }
        if (matchingMode == 1) {
          cleanMatched = mcpjet.has_matchedJetCand(); // HF matching only
        }
        if (cleanMatched) {
          registry.fill(HIST("MC_JetCounter"), getValFromBin(ParticleLevelJetWithMatchedCandidate));
        }

        if (fillLegacyOutputs) {
          registry.fill(HIST("hSparse_d0_mcp_eff"), d0.pt(), mcpjet.pt(), kin.zParallel, origin, cleanMatched ? 1 : 0,
                        obs.lambda, w);
        }
        const int o = originIndex(origin);
        if (o < 0) {
          continue;
        }
        if (o == 0) {
          registry.fill(HIST("h_eff_run2_den_pr"), d0.pt(), w);
        } else {
          registry.fill(HIST("h_eff_run2_den_np"), d0.pt(), w);
        }
        const std::array<double, 4> xg{mcpjet.pt(), obs.lambda, d0.pt(), static_cast<double>(sample)};
        hGen[o]->Fill(xg.data(), w);
      } // MCP jets
    } // MC collisions
  }
  PROCESS_SWITCH(JetHFAngularityTask, processMCDChargedSubstructureMatched,
                 "MC: templates + 6D responses + efficiencies + closure split (needs jet-matching-mc-d0-ch)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetHFAngularityTask>(cfgc)};
}
