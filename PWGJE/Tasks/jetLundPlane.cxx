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

/// \file jetLundPlane.cxx
/// \brief Task for jet Lund plane. Creates histograms for offline unfolding (including QA histos), and optionally tables.
/// \author Zoltan Varga <zoltan.varga@cern.ch>

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <THnSparse.h>

#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/PseudoJet.hh>

#include "Common/Core/RecoDecay.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Mini-AOD (tables)
namespace o2::aod
{
// Parent table: one row per collision stored in the MiniAOD
DECLARE_SOA_TABLE(MiniCollisions, "AOD", "MINICOLL");

// MiniJets -> MiniCollisions
DECLARE_SOA_INDEX_COLUMN(MiniCollision, miniCollision);

// Jet payload
DECLARE_SOA_COLUMN(Level, level, uint8_t);     // JetLevel::Det=reco(det), JetLevel::Part=truth(part)
DECLARE_SOA_COLUMN(JetRint, jetRint, int32_t); // jet.r() as stored (int R*100)
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);

DECLARE_SOA_TABLE(MiniJets, "AOD", "MINIJET",
                  MiniCollisionId, Level, JetRint, JetPt, JetEta, JetPhi);

// MiniSplittings -> MiniJets
DECLARE_SOA_INDEX_COLUMN(MiniJet, miniJet);

// Per-splitting observables (primary branch)
DECLARE_SOA_COLUMN(SplitId, splitId, uint16_t);
DECLARE_SOA_COLUMN(DeltaR, deltaR, float);
DECLARE_SOA_COLUMN(PtSoft, ptSoft, float);
DECLARE_SOA_COLUMN(PtHard, ptHard, float);
DECLARE_SOA_COLUMN(SoftEta, softEta, float);
DECLARE_SOA_COLUMN(SoftPhi, softPhi, float);

DECLARE_SOA_TABLE(MiniSplittings, "AOD", "MINISPL",
                  MiniJetId, SplitId, DeltaR, PtSoft, PtHard, SoftEta, SoftPhi, JetPt);

// Jet-jet matching (MC)
DECLARE_SOA_COLUMN(MatchDR, matchDR, float);
DECLARE_SOA_COLUMN(MatchRelPt, matchRelPt, float);

DECLARE_SOA_INDEX_COLUMN_FULL(DetMiniJet, detMiniJet, int, MiniJets, "_Det");
DECLARE_SOA_INDEX_COLUMN_FULL(PartMiniJet, partMiniJet, int, MiniJets, "_Part");

DECLARE_SOA_TABLE(MiniJetMatches, "AOD", "MINIMCH",
                  DetMiniJetId,
                  PartMiniJetId,
                  MatchDR,
                  MatchRelPt);
} // namespace o2::aod

namespace
{
constexpr float kTiny = 1e-12f;

struct JetLevel {
  enum Type : uint8_t {
    Det = 0,
    Part = 1
  };
};

template <typename JetT>
float jetRfromTable(const JetT& jet)
{
  return static_cast<float>(jet.r()) * 0.01f;
}

struct SplittingObs {
  float deltaR{};
  float ptSoft{};
  float ptHard{};
  float kT{};
  float lnRoverDR{};
  float lnkt{};
  float z{};
  float softEta{};
  float softPhi{};
};

struct SplitMatchPair {
  size_t recoIdx{};
  size_t truthIdx{};
  float dR{};
};

struct JetMatchInfo {
  uint64_t otherKey{};
  float dR{1e9f};
  float relPt{-1.f};
  float otherPt{};
};

inline float deltaPhi(float phi1, float phi2)
{
  return RecoDecay::constrainAngle(phi1 - phi2, -o2::constants::math::PI);
}

inline float splittingMatchDistance(const SplittingObs& a, const SplittingObs& b)
{
  return std::hypot(a.softEta - b.softEta, deltaPhi(a.softPhi, b.softPhi));
}

std::vector<SplitMatchPair> buildUniqueSplittingMatches(const std::vector<SplittingObs>& recoSpl,
                                                        const std::vector<SplittingObs>& truthSpl,
                                                        float maxDR)
{
  std::vector<SplitMatchPair> matches;
  matches.reserve(std::min(recoSpl.size(), truthSpl.size()));

  if (recoSpl.empty() || truthSpl.empty()) {
    return matches;
  }

  // Mutual-nearest-neighbor matching in the soft-prong (eta,phi) coordinates:
  // 1) for each truth splitting, find the closest reco splitting within maxDR
  // 2) for that reco splitting, find the closest truth splitting within maxDR
  // 3) keep the pair only if the reco splitting points back to the original truth splitting
  for (size_t truthIdx = 0; truthIdx < truthSpl.size(); ++truthIdx) {
    float bestRecoDR = maxDR;
    size_t bestRecoIdx = recoSpl.size();

    for (size_t recoIdx = 0; recoIdx < recoSpl.size(); ++recoIdx) {
      const float dR = splittingMatchDistance(recoSpl[recoIdx], truthSpl[truthIdx]);
      if (dR < bestRecoDR) {
        bestRecoDR = dR;
        bestRecoIdx = recoIdx;
      }
    }

    if (bestRecoIdx == recoSpl.size()) {
      continue;
    }

    float bestTruthDR = maxDR;
    size_t reverseTruthIdx = truthSpl.size();
    for (size_t candTruthIdx = 0; candTruthIdx < truthSpl.size(); ++candTruthIdx) {
      const float dR = splittingMatchDistance(recoSpl[bestRecoIdx], truthSpl[candTruthIdx]);
      if (dR < bestTruthDR) {
        bestTruthDR = dR;
        reverseTruthIdx = candTruthIdx;
      }
    }

    if (reverseTruthIdx == truthIdx) {
      matches.push_back({bestRecoIdx, truthIdx, bestRecoDR});
    }
  }

  return matches;
}

template <typename ConstituentRangeT>
std::vector<fastjet::PseudoJet> buildFastJetInputs(ConstituentRangeT&& constituents, float trackPtMin)
{
  std::vector<fastjet::PseudoJet> fjInputs;
  fjInputs.reserve(64);
  for (auto const& c : constituents) {
    if (c.pt() < trackPtMin) {
      continue;
    }
    const float mPi = o2::constants::physics::MassPiPlus;
    const float e = std::sqrt(c.p() * c.p() + mPi * mPi);
    fjInputs.emplace_back(c.px(), c.py(), c.pz(), e);
  }
  return fjInputs;
}

std::vector<SplittingObs> primaryDeclusteringSplittings(fastjet::PseudoJet jetCA, float jetR)
{
  std::vector<SplittingObs> out;
  out.reserve(32);

  fastjet::PseudoJet p1, p2;
  while (jetCA.has_parents(p1, p2)) {
    if (p2.perp() > p1.perp()) {
      std::swap(p1, p2);
    }
    const float pt1 = p1.perp();
    const float pt2 = p2.perp();
    const float pt = pt1 + pt2;
    if (pt <= kTiny) {
      break;
    }
    const float z = pt2 / pt;
    const float dR = p1.delta_R(p2);
    const float kT = pt2 * std::sin(dR);

    if (dR > kTiny && kT > kTiny && z > kTiny) {
      SplittingObs s;
      s.deltaR = dR;
      s.ptSoft = pt2;
      s.ptHard = pt1;
      s.kT = kT;
      s.lnRoverDR = std::log(jetR / dR);
      s.lnkt = std::log(kT);
      s.z = z;
      s.softEta = p2.eta();
      s.softPhi = p2.phi_std();
      out.push_back(s);
    }
    jetCA = p1; // follow hard branch (primary Lund)
  }
  return out;
}

} // namespace

struct JetLundPlaneUnfolding {
  // Config
  Configurable<float> vertexZCut{"vertexZCut", 10.f, "|z_vtx| cut"};
  Configurable<float> jetPtMin{"jetPtMin", 20.f, "min reco jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.5f, "min jet eta"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.5f, "max jet eta"};
  Configurable<float> jetR{"jetR", 0.4f, "jet radius (must match derived tables)"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15f, "min constituent pT"};

  Configurable<int> nBinsJetPt{"nBinsJetPt", 200, "jet pT bins"};
  Configurable<float> jetPtMax{"jetPtMax", 200.f, "jet pT max"};

  Configurable<int> nBinsLund{"nBinsLund", 120, "lund bins (lnR/DR and lnkt)"};
  Configurable<float> lnRoverDRMin{"lnRoverDRMin", -1.f, "min ln(R/DR)"};
  Configurable<float> lnRoverDRMax{"lnRoverDRMax", 6.f, "max ln(R/DR)"};
  Configurable<float> lnKtMin{"lnKtMin", -6.f, "min ln(kT)"};
  Configurable<float> lnKtMax{"lnKtMax", 6.f, "max ln(kT)"};

  // switches (runtime)
  Configurable<bool> writeMiniAOD{"writeMiniAOD", true, "write mini-AOD tables (jets, splittings, matching)"};

  // matching knobs (applied on top of matching relations)
  Configurable<float> matchMaxDR{"matchMaxDR", 0.2f, "max ΔR between det and part jet"};
  Configurable<float> splitMatchMaxDR{"splitMatchMaxDR", 0.1f, "max ΔR between reco and truth splittings (soft prong)"};
  Configurable<bool> matchUseRelPt{"matchUseRelPt", true, "apply relative pT compatibility cut"};
  Configurable<float> matchMaxRelPtDiff{"matchMaxRelPtDiff", 0.5f, "max |pTdet-pTpart|/pTpart"};

  // Registry
  HistogramRegistry registry{"registry"};

  // Mini-AOD outputs (optional)
  Produces<aod::MiniCollisions> outMiniCollisions;
  Produces<aod::MiniJets> outMiniJets;
  Produces<aod::MiniSplittings> outMiniSplittings;
  Produces<aod::MiniJetMatches> outMiniJetMatches;

  // FastJet reclustering setup (C/A)
  JetFinder reclusterer;
  std::vector<fastjet::PseudoJet> jetReclustered;

  void init(InitContext const&)
  {
    const int nbPt = nBinsJetPt.value;
    const float ptMax = jetPtMax.value;
    const int nbL = nBinsLund.value;

    registry.add("hEventCount", "Event counter;step;counts", HistType::kTH1F, {{10, 0.f, 10.f}});
    registry.add("hJetCountSummary", "Jet count summary;category;counts", HistType::kTH1F, {{6, 0.5f, 6.5f}});

    // Jet spectra / matching QA
    registry.add("hJetPtRecoAll", "Reco jets;p_{T}^{jet} (GeV/c);counts", HistType::kTH1F, {{nbPt, 0.f, ptMax}});
    registry.add("hJetPtRecoMatched", "Reco matched jets;p_{T}^{jet} (GeV/c);counts", HistType::kTH1F, {{nbPt, 0.f, ptMax}});
    registry.add("hJetPtRecoFake", "Reco unmatched (fake) jets;p_{T}^{jet} (GeV/c);counts", HistType::kTH1F, {{nbPt, 0.f, ptMax}});

    registry.add("hJetPtTruthAll", "Truth jets;p_{T}^{jet} (GeV/c);counts", HistType::kTH1F, {{nbPt, 0.f, ptMax}});
    registry.add("hJetPtTruthMatched", "Truth matched jets;p_{T}^{jet} (GeV/c);counts", HistType::kTH1F, {{nbPt, 0.f, ptMax}});
    registry.add("hJetPtTruthMiss", "Truth unmatched (miss) jets;p_{T}^{jet} (GeV/c);counts", HistType::kTH1F, {{nbPt, 0.f, ptMax}});

    // Jet pT response (1D unfolding)
    registry.add("hJetPtResponse", "Jet pT response;p_{T}^{det};p_{T}^{part}", HistType::kTH2F,
                 {{nbPt, 0.f, ptMax}, {nbPt, 0.f, ptMax}});

    // Lund reco / truth (3D)
    registry.add("hLundReco3D", "Primary Lund (reco);ln(R/#DeltaR);ln(k_{T});p_{T}^{jet}",
                 HistType::kTH3F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value},
                  {nbL, lnKtMin.value, lnKtMax.value},
                  {nbPt, 0.f, ptMax}});

    registry.add("hLundTruth3D", "Primary Lund (truth);ln(R/#DeltaR);ln(k_{T});p_{T}^{jet}",
                 HistType::kTH3F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value},
                  {nbL, lnKtMin.value, lnKtMax.value},
                  {nbPt, 0.f, ptMax}});

    // 6D response for 3D unfolding:
    // (det lnR/DR, det lnkt, det jetpt, part lnR/DR, part lnkt, part jetpt)
    registry.add("hLundResponse6D",
                 "Lund response 6D;"
                 "ln(R/#DeltaR)_{det};ln(k_{T})_{det};p_{T}^{det};"
                 "ln(R/#DeltaR)_{part};ln(k_{T})_{part};p_{T}^{part}",
                 HistType::kTHnSparseF,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value},
                  {nbL, lnKtMin.value, lnKtMax.value},
                  {nbPt, 0.f, ptMax},
                  {nbL, lnRoverDRMin.value, lnRoverDRMax.value},
                  {nbL, lnKtMin.value, lnKtMax.value},
                  {nbPt, 0.f, ptMax}});

    // Early QA histograms for declustering and matching
    registry.add("hNSplittingsReco", "Number of primary splittings (reco);N_{split};counts", HistType::kTH1F, {{40, -0.5f, 39.5f}});
    registry.add("hNSplittingsTruth", "Number of primary splittings (truth);N_{split};counts", HistType::kTH1F, {{40, -0.5f, 39.5f}});
    registry.add("hNSplittingsVsJetPtReco", "Reco primary splittings;#it{p}_{T}^{jet} (GeV/#it{c});N_{split}", HistType::kTH2F, {{nbPt, 0.f, ptMax}, {40, -0.5f, 39.5f}});
    registry.add("hNSplittingsVsJetPtTruth", "Truth primary splittings;#it{p}_{T}^{jet} (GeV/#it{c});N_{split}", HistType::kTH2F, {{nbPt, 0.f, ptMax}, {40, -0.5f, 39.5f}});

    registry.add("hDeltaRReco", "Reco splitting #DeltaR;#DeltaR;counts", HistType::kTH1F, {{120, 0.f, 1.2f}});
    registry.add("hDeltaRTruth", "Truth splitting #DeltaR;#DeltaR;counts", HistType::kTH1F, {{120, 0.f, 1.2f}});
    registry.add("hLnRoverDRReco", "Reco splitting ln(R/#DeltaR);ln(R/#DeltaR);counts", HistType::kTH1F, {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}});
    registry.add("hLnRoverDRTruth", "Truth splitting ln(R/#DeltaR);ln(R/#DeltaR);counts", HistType::kTH1F, {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}});
    registry.add("hLnKtReco", "Reco splitting ln(k_{T});ln(k_{T});counts", HistType::kTH1F, {{nbL, lnKtMin.value, lnKtMax.value}});
    registry.add("hLnKtTruth", "Truth splitting ln(k_{T});ln(k_{T});counts", HistType::kTH1F, {{nbL, lnKtMin.value, lnKtMax.value}});
    registry.add("hKtReco", "Reco splitting k_{T};k_{T};counts", HistType::kTH1F, {{200, 0.f, 20.f}});
    registry.add("hKtTruth", "Truth splitting k_{T};k_{T};counts", HistType::kTH1F, {{200, 0.f, 20.f}});

    registry.add("hJetMatchDR", "Matched jet #DeltaR;#DeltaR(det,part);counts", HistType::kTH1F, {{100, 0.f, 1.f}});
    registry.add("hJetMatchRelPt", "Matched jet relative #it{p}_{T} difference;|#it{p}_{T}^{det}-#it{p}_{T}^{part}|/#it{p}_{T}^{part};counts", HistType::kTH1F, {{100, 0.f, 2.f}});
    registry.add("hJetPtResidual", "Jet #it{p}_{T} residual;(#it{p}_{T}^{det}-#it{p}_{T}^{part})/#it{p}_{T}^{part};counts", HistType::kTH1F, {{160, -2.f, 2.f}});
    registry.add("hLnRoverDRResidual", "Splitting ln(R/#DeltaR) residual;ln(R/#DeltaR)_{det}-ln(R/#DeltaR)_{part};counts", HistType::kTH1F, {{160, -4.f, 4.f}});
    registry.add("hLnKtResidual", "Splitting ln(k_{T}) residual;ln(k_{T})_{det}-ln(k_{T})_{part};counts", HistType::kTH1F, {{160, -6.f, 6.f}});
    registry.add("hSplitMatchDR", "Matched splitting #DeltaR (soft prong);#DeltaR_{match};counts", HistType::kTH1F, {{120, 0.f, 1.2f}});

    registry.add("hSplitTruthDen", "Truth splitting denominator;ln(R/#DeltaR);ln(k_{T})", HistType::kTH2F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}, {nbL, lnKtMin.value, lnKtMax.value}});
    registry.add("hSplitTruthMatched", "Truth matched splitting numerator;ln(R/#DeltaR);ln(k_{T})", HistType::kTH2F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}, {nbL, lnKtMin.value, lnKtMax.value}});
    registry.add("hSplitTruthMiss", "Truth missed splittings;ln(R/#DeltaR);ln(k_{T})", HistType::kTH2F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}, {nbL, lnKtMin.value, lnKtMax.value}});

    registry.add("hSplitRecoDen", "Reco splitting denominator;ln(R/#DeltaR);ln(k_{T})", HistType::kTH2F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}, {nbL, lnKtMin.value, lnKtMax.value}});
    registry.add("hSplitRecoMatched", "Reco matched splitting numerator;ln(R/#DeltaR);ln(k_{T})", HistType::kTH2F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}, {nbL, lnKtMin.value, lnKtMax.value}});
    registry.add("hSplitRecoFake", "Reco fake splittings;ln(R/#DeltaR);ln(k_{T})", HistType::kTH2F,
                 {{nbL, lnRoverDRMin.value, lnRoverDRMax.value}, {nbL, lnKtMin.value, lnKtMax.value}});

    // reclusterer config
    reclusterer.isReclustering = true;
    reclusterer.algorithm = fastjet::cambridge_algorithm;
    reclusterer.jetR = 1.0; // recluster radius for declustering tree
    reclusterer.recombScheme = fastjet::E_scheme;
    reclusterer.strategy = fastjet::Best;
    reclusterer.areaType = fastjet::active_area;
  }

  // Filters for reco jets (data or MC reco)
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut.node());
  Filter jetFilter = (aod::jet::pt > jetPtMin.node()) &&
                     (aod::jet::eta > jetEtaMin.node()) &&
                     (aod::jet::eta < jetEtaMax.node()) &&
                     (aod::jet::r == nround(jetR.node() * 100.f));

  // Type aliases
  using RecoJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
  using DetJetsMatched = soa::Join<aod::ChargedMCDetectorLevelJets,
                                   aod::ChargedMCDetectorLevelJetConstituents,
                                   aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
  using PartJetsMatched = soa::Join<aod::ChargedMCParticleLevelJets,
                                    aod::ChargedMCParticleLevelJetConstituents,
                                    aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

  template <typename JetT>
  bool passJetFiducial(JetT const& jet, int rWanted) const
  {
    return (jet.r() == rWanted) && (jet.pt() > jetPtMin.value) && (jet.eta() > jetEtaMin.value) && (jet.eta() < jetEtaMax.value);
  }

  void fillJetCountSummary(int bin)
  {
    registry.fill(HIST("hJetCountSummary"), static_cast<float>(bin));
  }

  // Helpers to fill Lund from a jet row
  template <typename JetRowT, typename ConstituentTableT>
  std::vector<SplittingObs> getPrimarySplittings(JetRowT const& jet, ConstituentTableT const&)
  {
    auto fjInputs = buildFastJetInputs(jet.template tracks_as<ConstituentTableT>(), trackPtMin.value);
    if (fjInputs.size() < 2) {
      return {};
    }

    jetReclustered.clear();
    fastjet::ClusterSequenceArea cs(reclusterer.findJets(fjInputs, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    if (jetReclustered.empty()) {
      return {};
    }

    const float rjet = jetRfromTable(jet);
    return primaryDeclusteringSplittings(jetReclustered[0], rjet);
  }

  void fillPrimaryLund3DFromSplittings(std::vector<SplittingObs> const& spl, float jetPt, bool isTruth, float weight = 1.f)
  {
    for (auto const& s : spl) {
      if (isTruth) {
        registry.fill(HIST("hLundTruth3D"), s.lnRoverDR, s.lnkt, jetPt, weight);
      } else {
        registry.fill(HIST("hLundReco3D"), s.lnRoverDR, s.lnkt, jetPt, weight);
      }
    }
  }

  template <typename JetRowT, typename ConstituentTableT>
  void fillPrimaryLund3D(JetRowT const& jet, ConstituentTableT const& constituents, bool isTruth, float weight = 1.f)
  {
    fillPrimaryLund3DFromSplittings(getPrimarySplittings(jet, constituents), jet.pt(), isTruth, weight);
  }

  void fillSplittingQAHists(std::vector<SplittingObs> const& splittings, bool isTruth, float jetPt)
  {
    if (isTruth) {
      registry.fill(HIST("hNSplittingsTruth"), static_cast<float>(splittings.size()));
      registry.fill(HIST("hNSplittingsVsJetPtTruth"), jetPt, static_cast<float>(splittings.size()));
    } else {
      registry.fill(HIST("hNSplittingsReco"), static_cast<float>(splittings.size()));
      registry.fill(HIST("hNSplittingsVsJetPtReco"), jetPt, static_cast<float>(splittings.size()));
    }

    for (auto const& s : splittings) {
      if (isTruth) {
        registry.fill(HIST("hDeltaRTruth"), s.deltaR);
        registry.fill(HIST("hLnRoverDRTruth"), s.lnRoverDR);
        registry.fill(HIST("hLnKtTruth"), s.lnkt);
        registry.fill(HIST("hKtTruth"), s.kT);
      } else {
        registry.fill(HIST("hDeltaRReco"), s.deltaR);
        registry.fill(HIST("hLnRoverDRReco"), s.lnRoverDR);
        registry.fill(HIST("hLnKtReco"), s.lnkt);
        registry.fill(HIST("hKtReco"), s.kT);
      }
    }
  }

  std::vector<SplitMatchPair> fillSplittingCorrectionHists(std::vector<SplittingObs> const& recoSpl, std::vector<SplittingObs> const& truthSpl)
  {
    for (auto const& s : truthSpl) {
      registry.fill(HIST("hSplitTruthDen"), s.lnRoverDR, s.lnkt);
    }
    for (auto const& s : recoSpl) {
      registry.fill(HIST("hSplitRecoDen"), s.lnRoverDR, s.lnkt);
    }

    const auto matches = buildUniqueSplittingMatches(recoSpl, truthSpl, splitMatchMaxDR.value);
    std::vector<bool> recoUsed(recoSpl.size(), false);
    std::vector<bool> truthUsed(truthSpl.size(), false);

    for (const auto& m : matches) {
      recoUsed[m.recoIdx] = true;
      truthUsed[m.truthIdx] = true;
      registry.fill(HIST("hSplitRecoMatched"), recoSpl[m.recoIdx].lnRoverDR, recoSpl[m.recoIdx].lnkt);
      registry.fill(HIST("hSplitTruthMatched"), truthSpl[m.truthIdx].lnRoverDR, truthSpl[m.truthIdx].lnkt);
      registry.fill(HIST("hSplitMatchDR"), m.dR);
    }

    for (size_t i = 0; i < truthSpl.size(); ++i) {
      if (!truthUsed[i]) {
        registry.fill(HIST("hSplitTruthMiss"), truthSpl[i].lnRoverDR, truthSpl[i].lnkt);
      }
    }
    for (size_t i = 0; i < recoSpl.size(); ++i) {
      if (!recoUsed[i]) {
        registry.fill(HIST("hSplitRecoFake"), recoSpl[i].lnRoverDR, recoSpl[i].lnkt);
      }
    }
    return matches;
  }

  // DATA / RECO PROCESSING
  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   soa::Filtered<RecoJets> const& jets,
                   aod::JetTracks const& tracks)
  {
    registry.fill(HIST("hEventCount"), 0.5);

    (void)collisions; // collision ids are used only transiently for grouping MiniJets by source event

    int miniCollIdx = -1;
    if (writeMiniAOD.value) {
      outMiniCollisions();
      miniCollIdx = outMiniCollisions.lastIndex();
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("hJetPtRecoAll"), jet.pt());
      auto spl = getPrimarySplittings(jet, tracks);
      fillPrimaryLund3DFromSplittings(spl, jet.pt(), /*isTruth*/ false);
      fillSplittingQAHists(spl, /*isTruth*/ false, jet.pt());

      if (writeMiniAOD.value) {
        outMiniJets(miniCollIdx, /*level*/ JetLevel::Det, jet.r(), jet.pt(), jet.eta(), jet.phi());
        const int miniJetIdx = outMiniJets.lastIndex();

        uint16_t sid = 0;
        for (auto const& s : spl) {
          outMiniSplittings(miniJetIdx, sid++, s.deltaR, s.ptSoft, s.ptHard, s.softEta, s.softPhi, jet.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(JetLundPlaneUnfolding, processData, "Reco/data Lund + jet spectra", true);

  // MC PROCESSING (det + part + response)

  void processMC(DetJetsMatched const& detJets,
                 PartJetsMatched const& partJets,
                 aod::JetCollisions const& collisions,
                 aod::JetTracks const& tracks,
                 aod::JetParticles const& particles)
  {
    registry.fill(HIST("hEventCount"), 1.5);

    (void)collisions; // collision ids are used only transiently for grouping MiniJets by source event

    const int rWanted = static_cast<int>(std::lround(jetR.value * 100.f));
    std::unordered_map<uint64_t, bool> truthMatchedById;
    std::unordered_map<uint64_t, std::vector<SplittingObs>> truthSplittingsById;
    std::unordered_map<uint64_t, JetMatchInfo> truthBestDet;
    std::unordered_set<uint64_t> detSplittingsWritten;
    std::unordered_set<uint64_t> truthSplittingsWritten;
    auto h6 = registry.get<THnSparse>(HIST("hLundResponse6D"));

    // Transient maps used only during this task execution to avoid duplicating
    // MiniCollision / MiniJet rows and to translate framework-local indices into
    // merge-safe MiniAOD row indices.
    std::unordered_map<uint64_t, int> detMiniCollByKey;
    std::unordered_map<uint64_t, int> partMiniCollByKey;
    std::unordered_map<uint64_t, int> detJetToMiniJetIdx;
    std::unordered_map<uint64_t, int> partJetToMiniJetIdx;

    // --- Truth pass ---
    // Fill inclusive truth histograms, cache truth splittings, write all accepted
    // truth jets to MiniJets, and determine the best detector candidate for each truth jet.
    for (auto const& partJet : partJets) {
      if (!passJetFiducial(partJet, rWanted)) {
        continue;
      }

      registry.fill(HIST("hJetPtTruthAll"), partJet.pt());
      fillJetCountSummary(4);

      const uint64_t truthJetKey = partJet.globalIndex();
      auto spl = getPrimarySplittings(partJet, particles);
      truthSplittingsById[truthJetKey] = spl;

      fillPrimaryLund3DFromSplittings(spl, partJet.pt(), /*isTruth*/ true);
      fillSplittingQAHists(spl, /*isTruth*/ true, partJet.pt());

      if (writeMiniAOD.value) {
        const uint64_t partCollKey = (static_cast<uint64_t>(partJet.mcCollisionId()) << 1U) | 1ULL;
        int partMiniCollIdx = -1;
        auto collIt = partMiniCollByKey.find(partCollKey);
        if (collIt == partMiniCollByKey.end()) {
          outMiniCollisions();
          partMiniCollIdx = outMiniCollisions.lastIndex();
          partMiniCollByKey.emplace(partCollKey, partMiniCollIdx);
        } else {
          partMiniCollIdx = collIt->second;
        }

        outMiniJets(partMiniCollIdx, /*level*/ JetLevel::Part, partJet.r(), partJet.pt(), partJet.eta(), partJet.phi());
        partJetToMiniJetIdx[truthJetKey] = outMiniJets.lastIndex();
      }

      if (!partJet.has_matchedJetGeo()) {
        continue;
      }

      JetMatchInfo bestDet{};
      bool foundDet = false;
      for (auto const& candDetJet : partJet.template matchedJetGeo_as<DetJetsMatched>()) {
        if (!passJetFiducial(candDetJet, rWanted)) {
          continue;
        }

        const float dR = std::hypot(candDetJet.eta() - partJet.eta(),
                                    deltaPhi(candDetJet.phi(), partJet.phi()));
        if (dR > matchMaxDR.value) {
          continue;
        }

        const float rel = std::abs(candDetJet.pt() - partJet.pt()) / std::max(partJet.pt(), kTiny);
        if (matchUseRelPt.value && rel > matchMaxRelPtDiff.value) {
          continue;
        }

        if (!foundDet || dR < bestDet.dR) {
          bestDet.otherKey = candDetJet.globalIndex();
          bestDet.dR = dR;
          bestDet.relPt = rel;
          bestDet.otherPt = candDetJet.pt();
          foundDet = true;
        }
      }
      if (foundDet) {
        truthBestDet[truthJetKey] = bestDet;
      }
    }

    // --- Detector loop ---
    // Write all accepted detector jets to MiniJets. A final matched pair is accepted only
    // if the detector jet and truth jet are mutual best matches under the same cuts.
    for (auto const& detJet : detJets) {
      if (!passJetFiducial(detJet, rWanted)) {
        continue;
      }

      registry.fill(HIST("hJetPtRecoAll"), detJet.pt());
      fillJetCountSummary(1);

      const uint64_t detJetKey = detJet.globalIndex();
      auto detSpl = getPrimarySplittings(detJet, tracks);
      fillPrimaryLund3DFromSplittings(detSpl, detJet.pt(), /*isTruth*/ false);
      fillSplittingQAHists(detSpl, /*isTruth*/ false, detJet.pt());

      if (writeMiniAOD.value) {
        const uint64_t detCollKey = (static_cast<uint64_t>(detJet.collisionId()) << 1U);
        int detMiniCollIdx = -1;
        auto collIt = detMiniCollByKey.find(detCollKey);
        if (collIt == detMiniCollByKey.end()) {
          outMiniCollisions();
          detMiniCollIdx = outMiniCollisions.lastIndex();
          detMiniCollByKey.emplace(detCollKey, detMiniCollIdx);
        } else {
          detMiniCollIdx = collIt->second;
        }

        outMiniJets(detMiniCollIdx, /*level*/ JetLevel::Det, detJet.r(), detJet.pt(), detJet.eta(), detJet.phi());
        detJetToMiniJetIdx[detJetKey] = outMiniJets.lastIndex();
      }

      if (!detJet.has_matchedJetGeo()) {
        registry.fill(HIST("hJetPtRecoFake"), detJet.pt());
        fillJetCountSummary(3);
        continue;
      }

      JetMatchInfo bestTruth{};
      bool foundMatch = false;
      std::vector<SplittingObs> bestPartSpl;

      for (auto const& candPartJet : detJet.template matchedJetGeo_as<PartJetsMatched>()) {
        if (!passJetFiducial(candPartJet, rWanted)) {
          continue;
        }

        const float dR = std::hypot(detJet.eta() - candPartJet.eta(),
                                    deltaPhi(detJet.phi(), candPartJet.phi()));
        if (dR > matchMaxDR.value) {
          continue;
        }

        const float rel = std::abs(detJet.pt() - candPartJet.pt()) / std::max(candPartJet.pt(), kTiny);
        if (matchUseRelPt.value && rel > matchMaxRelPtDiff.value) {
          continue;
        }

        if (!foundMatch || dR < bestTruth.dR) {
          const uint64_t candTruthKey = candPartJet.globalIndex();
          bestTruth.otherKey = candTruthKey;
          bestTruth.dR = dR;
          bestTruth.relPt = rel;
          bestTruth.otherPt = candPartJet.pt();

          auto splIt = truthSplittingsById.find(candTruthKey);
          if (splIt != truthSplittingsById.end()) {
            bestPartSpl = splIt->second;
          } else {
            bestPartSpl = getPrimarySplittings(candPartJet, particles);
            truthSplittingsById[candTruthKey] = bestPartSpl;
          }
          foundMatch = true;
        }
      }

      if (!foundMatch) {
        registry.fill(HIST("hJetPtRecoFake"), detJet.pt());
        fillJetCountSummary(3);
        continue;
      }

      const auto reverseIt = truthBestDet.find(bestTruth.otherKey);
      if (reverseIt == truthBestDet.end() || reverseIt->second.otherKey != detJetKey) {
        registry.fill(HIST("hJetPtRecoFake"), detJet.pt());
        fillJetCountSummary(3);
        continue;
      }

      const float bestDR = bestTruth.dR;
      const float bestRelPt = bestTruth.relPt;
      const uint64_t bestPartId = bestTruth.otherKey;
      const float bestPartPt = bestTruth.otherPt;

      registry.fill(HIST("hJetPtRecoMatched"), detJet.pt());
      fillJetCountSummary(2);
      registry.fill(HIST("hJetPtTruthMatched"), bestPartPt);
      fillJetCountSummary(5);
      registry.fill(HIST("hJetPtResponse"), detJet.pt(), bestPartPt);
      registry.fill(HIST("hJetMatchDR"), bestDR);
      registry.fill(HIST("hJetMatchRelPt"), bestRelPt);
      registry.fill(HIST("hJetPtResidual"), (detJet.pt() - bestPartPt) / std::max(bestPartPt, kTiny));
      truthMatchedById[bestPartId] = true;

      if (writeMiniAOD.value) {
        auto detIdxIt = detJetToMiniJetIdx.find(detJetKey);
        auto partIdxIt = partJetToMiniJetIdx.find(bestPartId);
        if (detIdxIt != detJetToMiniJetIdx.end() && partIdxIt != partJetToMiniJetIdx.end()) {
          const int detMiniJetIdx = detIdxIt->second;
          const int partMiniJetIdx = partIdxIt->second;

          outMiniJetMatches(detMiniJetIdx, partMiniJetIdx, bestDR, bestRelPt);

          if (detSplittingsWritten.insert(detJetKey).second) {
            uint16_t sidTmp = 0;
            for (auto const& s : detSpl) {
              outMiniSplittings(detMiniJetIdx, sidTmp++, s.deltaR, s.ptSoft, s.ptHard, s.softEta, s.softPhi, detJet.pt());
            }
          }

          if (truthSplittingsWritten.insert(bestPartId).second) {
            uint16_t sid = 0;
            for (auto const& s : bestPartSpl) {
              outMiniSplittings(partMiniJetIdx, sid++, s.deltaR, s.ptSoft, s.ptHard, s.softEta, s.softPhi, bestPartPt);
            }
          }
        }
      }

      const auto splitMatches = fillSplittingCorrectionHists(detSpl, bestPartSpl);

      // Fill the 6D response and residual QA using the unique geometric splitting matches
      for (const auto& m : splitMatches) {
        const auto& detS = detSpl[m.recoIdx];
        const auto& partS = bestPartSpl[m.truthIdx];
        double x[6] = {detS.lnRoverDR, detS.lnkt, detJet.pt(),
                       partS.lnRoverDR, partS.lnkt, bestPartPt};
        h6->Fill(x);
        registry.fill(HIST("hLnRoverDRResidual"), detS.lnRoverDR - partS.lnRoverDR);
        registry.fill(HIST("hLnKtResidual"), detS.lnkt - partS.lnkt);
      }
    }

    // --- Final truth pass for misses ---
    for (auto const& partJet : partJets) {
      if (!passJetFiducial(partJet, rWanted)) {
        continue;
      }

      const uint64_t truthJetKey = partJet.globalIndex();
      if (!truthMatchedById[truthJetKey]) {
        registry.fill(HIST("hJetPtTruthMiss"), partJet.pt());
        fillJetCountSummary(6);
      }
    }
  }
  PROCESS_SWITCH(JetLundPlaneUnfolding, processMC, "MC det+part + responses", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetLundPlaneUnfolding>(cfgc, TaskName{"jet-lund-plane"})};
}
