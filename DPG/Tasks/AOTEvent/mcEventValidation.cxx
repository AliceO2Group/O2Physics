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

/// \file mcEventValidation.cxx
/// \brief MC validation of the event reconstruction: generated multiplicity vs FT0C, vertex reconstruction efficiency and purity
///
/// Investigates the off-diagonal populations in Nch_gen vs FT0C and in Ntracks vs FT0C (O2-6499):
/// - MC truth: charged physical primaries (as in MultModule) split into those produced at the MC collision vertex and the displaced ones;
/// - reconstructed collisions: FT0C from the found BC, categorised by the BC association (correct BC, BC shared with other MC collisions, wrong BC);
/// - number of MC collisions per BC (same-bunch pile-up in MC), PV contributors from other MC collisions, number of reconstructed collisions per MC collision;
/// - vertex reconstruction efficiency and purity vs the generated multiplicity (|eta|<0.8) and vs the distance (in BCs) to the nearest other MC collision.
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TProfile2D.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using ColsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

struct McEventValidationTask {
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TracksMC> perCollision = aod::track::collisionId;

  // event selection (as in the Nch_gen vs FT0C studies: sel8, |vZ| < 10 cm, kNoSameBunchPileup)
  Configurable<bool> requireSel8{"requireSel8", true, "require sel8 for reconstructed collisions"};
  Configurable<bool> requireNoSameBunchPileup{"requireNoSameBunchPileup", true, "require kNoSameBunchPileup for reconstructed collisions"};
  Configurable<float> maxAbsVtxZ{"maxAbsVtxZ", 10.f, "max |vZ| of reconstructed collisions (cm)"};
  // MC particles
  Configurable<float> maxDistFromVertex{"maxDistFromVertex", 1.f, "particles produced within this distance (cm, in r and in z) from the MC collision vertex are counted as 'from the vertex'"};
  // tracks: PV contributors with ITS, and "global-like" tracks with TPC
  Configurable<float> minPt{"minPt", 0.2f, "min pT of tracks (GeV/c)"};
  Configurable<float> maxAbsEta{"maxAbsEta", 0.8f, "max |eta| of tracks"};
  Configurable<int> minITSNCls{"minITSNCls", 5, "min number of ITS clusters of tracks"};
  Configurable<int> minTPCNCls{"minTPCNCls", 50, "min number of TPC clusters of global-like tracks"};
  // vertex reconstruction efficiency: an MC collision is reconstructed if a reconstructed collision has more than minFractionMatch of its PV contributors from it
  // and the reconstructed vZ agrees with the generated one within maxDeltaVzMatch
  Configurable<float> minFractionMatch{"minFractionMatch", 0.5f, "min fraction of PV contributors from the MC collision for the matching"};
  Configurable<float> maxDeltaVzMatch{"maxDeltaVzMatch", 0.3f, "max |vZ(reco) - vZ(gen)| for the matching (cm)"};
  Configurable<int> maxDistBC{"maxDistBC", 800, "upper edge of the axis for the distance (in BCs) to the nearest other MC collision"};
  Configurable<int> nBinsDistBC{"nBinsDistBC", 160, "number of bins for the distance (in BCs) to the nearest other MC collision"};
  // axes
  Configurable<int> nBinsNch{"nBinsNch", 250, "number of bins for Nch and number of tracks"};
  Configurable<float> maxNch{"maxNch", 5000.f, "upper edge of the axes for Nch and number of tracks"};
  Configurable<int> nBinsFT0C{"nBinsFT0C", 300, "number of bins for the FT0C amplitude"};
  Configurable<float> maxFT0C{"maxFT0C", 60000.f, "upper edge of the axis for the FT0C amplitude"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // categories of reconstructed collisions by the association of the found BC to the true BC of the MC collision
  enum BcCategory { kCorrectBC = 0, // found BC = BC of the MC collision, no other MC collisions in this BC
                    kSharedBC,      // found BC = BC of the MC collision, other MC collisions in this BC (same-bunch pile-up in MC)
                    kWrongBC,       // found BC != BC of the MC collision
                    kNBcCategories };

  static constexpr float MaxEtaNch = 0.5f;       // |eta| range for Nch, as in MultModule
  static constexpr float MinEtaFT0C = -3.3f;     // FT0C acceptance, as in MultModule
  static constexpr float MaxEtaFT0C = -2.1f;     // FT0C acceptance, as in MultModule
  static constexpr float MinCharge = 1e-3f;      // charged particles
  static constexpr int MaxDeltaBC = 50;          // axis range for found BC - true BC
  static constexpr int MaxNumMcCollsInBC = 10;   // axis range for the number of MC collisions in a BC
  static constexpr int MaxNumRecoPerMcColl = 10; // axis range for the number of reconstructed collisions per MC collision
  static constexpr float MaxEtaNch08 = 0.8f;     // |eta| range for Nch in the efficiency and purity studies
  static constexpr int MaxNchLow = 200;          // axis range for Nch (|eta|<0.8) at low multiplicity, one bin per particle
  static constexpr int NNchClasses = 4;          // classes of MC collisions by Nch (|eta|<0.8): [0, 10), [10, 20), [20, 30), >= 30
  static constexpr std::array<int, NNchClasses - 1> NchClassLimits = {10, 20, 30};

  // class of an MC collision by Nch (|eta|<0.8): index of the first class limit above nch
  static int nchClass(int nch)
  {
    return static_cast<int>(std::upper_bound(NchClassLimits.begin(), NchClassLimits.end(), nch) - NchClassLimits.begin());
  }

  void init(InitContext&)
  {
    const AxisSpec axisNch{nBinsNch.value, 0., maxNch.value, "N_{ch} (|#eta|<0.5)"};
    const AxisSpec axisNchFT0C{nBinsNch.value, 0., maxNch.value, "N_{ch} in FT0C acceptance (-3.3<#eta<-2.1)"};
    const AxisSpec axisNtracks{nBinsNch.value, 0., maxNch.value, "number of tracks"};
    const AxisSpec axisFT0C{nBinsFT0C.value, 0., maxFT0C.value, "FT0C amplitude (found BC)"};
    const AxisSpec axisB{200, 0., 20., "impact parameter b (fm)"};
    const AxisSpec axisCategory{kNBcCategories, -0.5, static_cast<double>(kNBcCategories) - 0.5, ""};

    // MC truth
    histos.add("MC/hNch", "MC collisions: charged physical primaries in |#eta|<0.5", kTH2F, {axisB, axisNch});
    histos.add("MC/hNchFromVertex", "MC collisions: charged physical primaries in |#eta|<0.5 produced at the vertex", kTH2F, {axisB, axisNch});
    histos.add("MC/pFractionDisplacedVsB", "MC collisions: fraction of N_{ch} (|#eta|<0.5) not produced at the vertex;impact parameter b (fm);fraction", kTProfile, {axisB});
    histos.add("MC/hNchVsNchFT0C", "MC collisions: all charged physical primaries", kTH2F, {axisNchFT0C, axisNch});
    histos.add("MC/hNchVsNchFT0CFromVertex", "MC collisions: charged physical primaries produced at the vertex", kTH2F, {axisNchFT0C, axisNch});
    // displaced charged physical primaries (all eta)
    histos.add("MC/Displaced/hPt", "displaced charged physical primaries;p_{T} (GeV/c)", kTH1F, {{500, 0., 0.05}});
    histos.add("MC/Displaced/hDeltaZ", "displaced charged physical primaries;z_{prod} - z_{vertex} (cm)", kTH1F, {{500, -500., 500.}});
    histos.add("MC/Displaced/hR", "displaced charged physical primaries;r_{prod} (cm)", kTH1F, {{500, 0., 500.}});
    histos.add("MC/Displaced/hSpecies", "displaced charged physical primaries", kTH1F, {{3, -0.5, 2.5}});
    histos.add("MC/Displaced/hGenStatus", "displaced charged physical primaries;generator status code", kTH1F, {{201, -1.5, 199.5}});
    histos.add("MC/FromVertex/hGenStatus", "charged physical primaries produced at the vertex;generator status code", kTH1F, {{201, -1.5, 199.5}});
    // same-bunch pile-up in MC and splitting of MC collisions into several reconstructed ones
    histos.add("MC/hNumMcCollsInBC", "number of MC collisions in the BC of a given MC collision;number of MC collisions in the BC", kTH1F, {{MaxNumMcCollsInBC, 0.5, MaxNumMcCollsInBC + 0.5}});
    histos.add("MC/hNumRecoPerMcColl", "MC collisions;number of reconstructed collisions", kTH1F, {{MaxNumRecoPerMcColl + 1, -0.5, MaxNumRecoPerMcColl + 0.5}});
    histos.add("MC/hNumRecoPerMcCollVsB", "MC collisions", kTH2F, {axisB, {MaxNumRecoPerMcColl + 1, -0.5, MaxNumRecoPerMcColl + 0.5, "number of reconstructed collisions"}});

    // reconstructed collisions
    histos.add("Reco/hCounter", "reconstructed collisions", kTH1F, {{6, -0.5, 5.5}});
    histos.add("Reco/hCategory", "selected reconstructed collisions", kTH1F, {axisCategory});
    histos.add("Reco/hDeltaBC", "selected reconstructed collisions;found BC - BC of the MC collision", kTH1F, {{2 * MaxDeltaBC + 1, -MaxDeltaBC - 0.5, MaxDeltaBC + 0.5}});
    histos.add("Reco/hNumMcCollsInFoundBC", "selected reconstructed collisions;number of MC collisions in the found BC", kTH1F, {{MaxNumMcCollsInBC + 1, -0.5, MaxNumMcCollsInBC + 0.5}});
    histos.add("Reco/hPurity", "selected reconstructed collisions;fraction of PV contributors from the labelled MC collision", kTH2F, {axisCategory, {110, 0., 1.1, "purity"}});
    histos.add("Reco/hNchVsFT0C", "selected reconstructed collisions: all charged physical primaries", kTH3F, {axisFT0C, axisNch, axisCategory});
    histos.add("Reco/hNchFromVertexVsFT0C", "selected reconstructed collisions: charged physical primaries produced at the vertex", kTH3F, {axisFT0C, axisNch, axisCategory});
    histos.add("Reco/hNchFT0CFromVertexVsFT0C", "selected reconstructed collisions: N_{ch} in FT0C acceptance produced at the vertex", kTH3F, {axisFT0C, axisNchFT0C, axisCategory});
    histos.add("Reco/hNPVVsFT0C", "selected reconstructed collisions: PV contributors with ITS", kTH3F, {axisFT0C, axisNtracks, axisCategory});
    histos.add("Reco/hNGlobalVsFT0C", "selected reconstructed collisions: global-like tracks", kTH3F, {axisFT0C, axisNtracks, axisCategory});
    histos.add("Reco/hNGlobalVsNchFromVertex", "selected reconstructed collisions: global-like tracks", kTH3F, {axisNch, axisNtracks, axisCategory});

    // vertex reconstruction efficiency and purity vs Nch (|eta|<0.8, produced at the vertex) and vs the distance to the nearest other MC collision
    const AxisSpec axisNch08Low{MaxNchLow, -0.5, MaxNchLow - 0.5, "N_{ch} (|#eta|<0.8)"};
    const AxisSpec axisNch08{nBinsNch.value, 0., maxNch.value, "N_{ch} (|#eta|<0.8)"};
    const AxisSpec axisDist{nBinsDistBC.value, 0., static_cast<double>(maxDistBC.value), "distance to the nearest other MC collision (BCs)"};
    const AxisSpec axisNchClass{NNchClasses, -0.5, NNchClasses - 0.5, ""};
    const AxisSpec axisPurity{110, 0., 1.1, "purity"};
    histos.add("Efficiency/hNchGen", "MC collisions", kTH1F, {axisNch08Low});
    histos.add("Efficiency/hNchGenReco", "MC collisions with a matched reconstructed collision", kTH1F, {axisNch08Low});
    histos.add("Efficiency/hNchGenWide", "MC collisions", kTH1F, {axisNch08});
    histos.add("Efficiency/hNchGenRecoWide", "MC collisions with a matched reconstructed collision", kTH1F, {axisNch08});
    histos.add("Efficiency/hDistGen", "MC collisions", kTH2F, {axisDist, axisNchClass});
    histos.add("Efficiency/hDistGenReco", "MC collisions with a matched reconstructed collision", kTH2F, {axisDist, axisNchClass});
    histos.add("Purity/pPurityPVVsNch", "reconstructed collisions;N_{ch} (|#eta|<0.8) of the labelled MC collision;purity of PV contributors", kTProfile, {axisNch08Low});
    histos.add("Purity/pPurityAssocVsNch", "reconstructed collisions;N_{ch} (|#eta|<0.8) of the labelled MC collision;purity of associated tracks", kTProfile, {axisNch08Low});
    histos.add("Purity/hPurityPVVsNch", "reconstructed collisions: PV contributors", kTH2F, {axisNch08, axisPurity});
    histos.add("Purity/hPurityAssocVsNch", "reconstructed collisions: associated tracks", kTH2F, {axisNch08, axisPurity});
    histos.add("Purity/pPurityPVVsDist", "reconstructed collisions: PV contributors;distance to the nearest other MC collision (BCs);N_{ch} class", kTProfile2D, {axisDist, axisNchClass});
    histos.add("Purity/pPurityAssocVsDist", "reconstructed collisions: associated tracks;distance to the nearest other MC collision (BCs);N_{ch} class", kTProfile2D, {axisDist, axisNchClass});
    for (int i = 0; i < NNchClasses; i++) {
      const char* label = nullptr;
      if (i == 0) {
        label = Form("N_{ch} < %d", NchClassLimits[i]);
      } else if (i == NNchClasses - 1) {
        label = Form("N_{ch} #geq %d", NchClassLimits[i - 1]);
      } else {
        label = Form("%d #leq N_{ch} < %d", NchClassLimits[i - 1], NchClassLimits[i]);
      }
      histos.get<TH2>(HIST("Efficiency/hDistGen"))->GetYaxis()->SetBinLabel(i + 1, label);
      histos.get<TH2>(HIST("Efficiency/hDistGenReco"))->GetYaxis()->SetBinLabel(i + 1, label);
      histos.get<TProfile2D>(HIST("Purity/pPurityPVVsDist"))->GetYaxis()->SetBinLabel(i + 1, label);
      histos.get<TProfile2D>(HIST("Purity/pPurityAssocVsDist"))->GetYaxis()->SetBinLabel(i + 1, label);
    }

    // bin labels
    const std::array<const char*, kNBcCategories> categoryNames = {"correct BC", "BC shared with other MC colls", "wrong BC"};
    for (int i = 0; i < kNBcCategories; i++) {
      histos.get<TH1>(HIST("Reco/hCategory"))->GetXaxis()->SetBinLabel(i + 1, categoryNames[i]);
      histos.get<TH2>(HIST("Reco/hPurity"))->GetXaxis()->SetBinLabel(i + 1, categoryNames[i]);
      for (const auto& h : {histos.get<TH3>(HIST("Reco/hNchVsFT0C")).get(), histos.get<TH3>(HIST("Reco/hNchFromVertexVsFT0C")).get(), histos.get<TH3>(HIST("Reco/hNchFT0CFromVertexVsFT0C")).get(),
                            histos.get<TH3>(HIST("Reco/hNPVVsFT0C")).get(), histos.get<TH3>(HIST("Reco/hNGlobalVsFT0C")).get(), histos.get<TH3>(HIST("Reco/hNGlobalVsNchFromVertex")).get()}) {
        h->GetZaxis()->SetBinLabel(i + 1, categoryNames[i]);
      }
    }
    auto hSpecies = histos.get<TH1>(HIST("MC/Displaced/hSpecies"));
    hSpecies->GetXaxis()->SetBinLabel(1, "e^{-}");
    hSpecies->GetXaxis()->SetBinLabel(2, "e^{+}");
    hSpecies->GetXaxis()->SetBinLabel(3, "other");
    auto hCounter = histos.get<TH1>(HIST("Reco/hCounter"));
    const std::array<const char*, 6> counterNames = {"all", "with MC label", "with found FT0", "sel8", "|vZ| cut", "kNoSameBunchPileup"};
    for (size_t i = 0; i < counterNames.size(); i++) {
      hCounter->GetXaxis()->SetBinLabel(i + 1, counterNames[i]);
    }
  }

  void process(ColsMC const& cols, aod::McCollisions const& mcCols, aod::McParticles const& mcParts, TracksMC const& tracks, aod::BCs const& bcs, aod::FT0s const&)
  {
    // ### MC truth: charged physical primaries (as in MultModule), split into produced at the MC collision vertex and displaced
    const size_t nMcColls = mcCols.size();
    std::vector<int> nch(nMcColls, 0), nchFromVertex(nMcColls, 0), nchFT0C(nMcColls, 0), nchFT0CFromVertex(nMcColls, 0), nch08FromVertex(nMcColls, 0);
    for (const auto& particle : mcParts) {
      if (!particle.isPhysicalPrimary() || !particle.has_mcCollision()) {
        continue;
      }
      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle == nullptr || std::abs(pdgParticle->Charge()) < MinCharge) {
        continue;
      }
      const auto mcCol = particle.mcCollision();
      const float dz = particle.vz() - mcCol.posZ();
      const float dr = std::hypot(particle.vx() - mcCol.posX(), particle.vy() - mcCol.posY());
      const bool isFromVertex = std::abs(dz) < maxDistFromVertex && dr < maxDistFromVertex;
      const float eta = particle.eta();
      const bool isInNchAcc = std::abs(eta) < MaxEtaNch;
      const bool isInFT0CAcc = eta > MinEtaFT0C && eta < MaxEtaFT0C;
      const auto iMc = particle.mcCollisionId();
      nch[iMc] += isInNchAcc;
      nchFromVertex[iMc] += isInNchAcc && isFromVertex;
      nchFT0C[iMc] += isInFT0CAcc;
      nchFT0CFromVertex[iMc] += isInFT0CAcc && isFromVertex;
      nch08FromVertex[iMc] += std::abs(eta) < MaxEtaNch08 && isFromVertex;
      if (isFromVertex) {
        histos.fill(HIST("MC/FromVertex/hGenStatus"), particle.getGenStatusCode());
        continue;
      }
      histos.fill(HIST("MC/Displaced/hPt"), particle.pt());
      histos.fill(HIST("MC/Displaced/hDeltaZ"), dz);
      histos.fill(HIST("MC/Displaced/hR"), std::hypot(particle.vx(), particle.vy()));
      histos.fill(HIST("MC/Displaced/hSpecies"), particle.pdgCode() == PDG_t::kElectron ? 0 : (particle.pdgCode() == PDG_t::kPositron ? 1 : 2));
      histos.fill(HIST("MC/Displaced/hGenStatus"), particle.getGenStatusCode());
    }

    // number of MC collisions per BC (same-bunch pile-up in MC)
    std::unordered_map<uint64_t, int> numMcCollsInBC;
    for (const auto& mcCol : mcCols) {
      numMcCollsInBC[mcCol.bc_as<aod::BCs>().globalBC()]++;
    }
    // number of reconstructed collisions per MC collision
    std::vector<int> numRecoPerMcColl(nMcColls, 0);
    for (const auto& col : cols) {
      if (col.has_mcCollision()) {
        numRecoPerMcColl[col.mcCollisionId()]++;
      }
    }

    for (const auto& mcCol : mcCols) {
      const auto iMc = mcCol.globalIndex();
      const float b = mcCol.impactParameter();
      histos.fill(HIST("MC/hNch"), b, nch[iMc]);
      histos.fill(HIST("MC/hNchFromVertex"), b, nchFromVertex[iMc]);
      if (nch[iMc] > 0) {
        histos.fill(HIST("MC/pFractionDisplacedVsB"), b, 1. - static_cast<double>(nchFromVertex[iMc]) / nch[iMc]);
      }
      histos.fill(HIST("MC/hNchVsNchFT0C"), nchFT0C[iMc], nch[iMc]);
      histos.fill(HIST("MC/hNchVsNchFT0CFromVertex"), nchFT0CFromVertex[iMc], nchFromVertex[iMc]);
      histos.fill(HIST("MC/hNumMcCollsInBC"), numMcCollsInBC[mcCol.bc_as<aod::BCs>().globalBC()]);
      histos.fill(HIST("MC/hNumRecoPerMcColl"), numRecoPerMcColl[iMc]);
      histos.fill(HIST("MC/hNumRecoPerMcCollVsB"), b, numRecoPerMcColl[iMc]);
    }

    // ### vertex reconstruction efficiency and purity
    // distance (in BCs) from each MC collision to the nearest other MC collision in this DF; -1 if it can be outside the DF (closer to the DF edges)
    int64_t firstBCinDF = std::numeric_limits<int64_t>::max(), lastBCinDF = 0;
    for (const auto& bc : bcs) {
      firstBCinDF = std::min(firstBCinDF, static_cast<int64_t>(bc.globalBC()));
      lastBCinDF = std::max(lastBCinDF, static_cast<int64_t>(bc.globalBC()));
    }
    std::vector<std::pair<int64_t, int64_t>> mcCollBCs; // (global BC, MC collision index)
    for (const auto& mcCol : mcCols) {
      mcCollBCs.emplace_back(static_cast<int64_t>(mcCol.bc_as<aod::BCs>().globalBC()), mcCol.globalIndex());
    }
    std::sort(mcCollBCs.begin(), mcCollBCs.end());
    std::vector<int64_t> distToNearest(nMcColls, -1);
    for (size_t k = 0; k < mcCollBCs.size(); k++) {
      const int64_t bc = mcCollBCs[k].first;
      int64_t dist = std::numeric_limits<int64_t>::max();
      if (k > 0) {
        dist = std::min(dist, bc - mcCollBCs[k - 1].first);
      }
      if (k + 1 < mcCollBCs.size()) {
        dist = std::min(dist, mcCollBCs[k + 1].first - bc);
      }
      if (dist <= std::min(bc - firstBCinDF, lastBCinDF - bc)) { // the nearest one cannot be outside the DF
        distToNearest[mcCollBCs[k].second] = dist;
      }
    }

    // matching of reconstructed collisions to MC collisions (PV contributors and vZ), and purity of the reconstructed collisions
    std::vector<bool> isReconstructed(nMcColls, false);
    std::unordered_map<int64_t, int> nPVPerMcColl, nAssocPerMcColl;
    for (const auto& col : cols) {
      nPVPerMcColl.clear();
      nAssocPerMcColl.clear();
      int nPVAll = 0, nAssocAll = 0;
      for (const auto& track : tracks.sliceBy(perCollision, col.globalIndex())) {
        const int64_t iMcOfTrack = track.has_mcParticle() ? mcParts.iteratorAt(track.mcParticleId()).mcCollisionId() : -1;
        nAssocAll++;
        nAssocPerMcColl[iMcOfTrack]++;
        if (track.isPVContributor()) {
          nPVAll++;
          nPVPerMcColl[iMcOfTrack]++;
        }
      }
      // matching: MC collision with most PV contributors
      int64_t iMcDominant = -1;
      int nPVDominant = 0;
      for (const auto& [iMc, n] : nPVPerMcColl) {
        if (iMc >= 0 && n > nPVDominant) {
          iMcDominant = iMc;
          nPVDominant = n;
        }
      }
      if (iMcDominant >= 0 && nPVDominant > minFractionMatch.value * nPVAll && std::abs(col.posZ() - mcCols.iteratorAt(iMcDominant).posZ()) < maxDeltaVzMatch.value) {
        isReconstructed[iMcDominant] = true;
      }
      // purity wrt the labelled MC collision
      if (!col.has_mcCollision()) {
        continue;
      }
      const auto iMc = col.mcCollisionId();
      const int nch08 = nch08FromVertex[iMc];
      const int64_t dist = distToNearest[iMc];
      const int iClass = nchClass(nch08);
      if (nPVAll > 0) {
        const double purityPV = static_cast<double>(nPVPerMcColl[iMc]) / nPVAll;
        histos.fill(HIST("Purity/pPurityPVVsNch"), nch08, purityPV);
        histos.fill(HIST("Purity/hPurityPVVsNch"), nch08, purityPV);
        if (dist >= 0) {
          histos.fill(HIST("Purity/pPurityPVVsDist"), dist, iClass, purityPV);
        }
      }
      if (nAssocAll > 0) {
        const double purityAssoc = static_cast<double>(nAssocPerMcColl[iMc]) / nAssocAll;
        histos.fill(HIST("Purity/pPurityAssocVsNch"), nch08, purityAssoc);
        histos.fill(HIST("Purity/hPurityAssocVsNch"), nch08, purityAssoc);
        if (dist >= 0) {
          histos.fill(HIST("Purity/pPurityAssocVsDist"), dist, iClass, purityAssoc);
        }
      }
    }
    // efficiency
    for (size_t iMc = 0; iMc < nMcColls; iMc++) {
      const int nch08 = nch08FromVertex[iMc];
      histos.fill(HIST("Efficiency/hNchGen"), nch08);
      histos.fill(HIST("Efficiency/hNchGenWide"), nch08);
      if (isReconstructed[iMc]) {
        histos.fill(HIST("Efficiency/hNchGenReco"), nch08);
        histos.fill(HIST("Efficiency/hNchGenRecoWide"), nch08);
      }
      if (distToNearest[iMc] < 0) {
        continue;
      }
      histos.fill(HIST("Efficiency/hDistGen"), distToNearest[iMc], nchClass(nch08));
      if (isReconstructed[iMc]) {
        histos.fill(HIST("Efficiency/hDistGenReco"), distToNearest[iMc], nchClass(nch08));
      }
    }

    // ### reconstructed collisions
    for (const auto& col : cols) {
      histos.fill(HIST("Reco/hCounter"), 0);
      if (!col.has_mcCollision()) {
        continue;
      }
      histos.fill(HIST("Reco/hCounter"), 1);
      if (!col.has_foundFT0()) {
        continue;
      }
      histos.fill(HIST("Reco/hCounter"), 2);
      if (requireSel8 && !col.sel8()) {
        continue;
      }
      histos.fill(HIST("Reco/hCounter"), 3);
      if (std::abs(col.posZ()) > maxAbsVtxZ) {
        continue;
      }
      histos.fill(HIST("Reco/hCounter"), 4);
      if (requireNoSameBunchPileup && !col.selection_bit(kNoSameBunchPileup)) {
        continue;
      }
      histos.fill(HIST("Reco/hCounter"), 5);

      const float ampFT0C = col.foundFT0().sumAmpC(); // as in MultModule: FT0 of the found BC
      const auto mcCol = col.mcCollision();
      const auto iMc = col.mcCollisionId();
      const auto foundBC = col.foundBC_as<aod::BCs>().globalBC();
      const auto trueBC = mcCol.bc_as<aod::BCs>().globalBC();
      const int numMcInFoundBC = numMcCollsInBC.contains(foundBC) ? numMcCollsInBC[foundBC] : 0;
      int category = kCorrectBC;
      if (foundBC != trueBC) {
        category = kWrongBC;
      } else if (numMcInFoundBC > 1) {
        category = kSharedBC;
      }

      // tracks: PV contributors with ITS, global-like tracks, and PV contributors from the labelled MC collision
      int nPV = 0, nGlobal = 0, nPVAll = 0, nPVFromLabelledColl = 0;
      for (const auto& track : tracks.sliceBy(perCollision, col.globalIndex())) {
        if (!track.isPVContributor()) {
          continue;
        }
        nPVAll++;
        if (track.has_mcParticle() && mcParts.iteratorAt(track.mcParticleId()).mcCollisionId() == iMc) {
          nPVFromLabelledColl++;
        }
        if (track.pt() < minPt || std::abs(track.eta()) > maxAbsEta || track.itsNCls() < minITSNCls) {
          continue;
        }
        nPV++;
        if (track.hasTPC() && track.tpcNClsFound() >= minTPCNCls) {
          nGlobal++;
        }
      }

      histos.fill(HIST("Reco/hCategory"), category);
      histos.fill(HIST("Reco/hDeltaBC"), static_cast<int64_t>(foundBC) - static_cast<int64_t>(trueBC));
      histos.fill(HIST("Reco/hNumMcCollsInFoundBC"), numMcInFoundBC);
      if (nPVAll > 0) {
        histos.fill(HIST("Reco/hPurity"), category, static_cast<double>(nPVFromLabelledColl) / nPVAll);
      }
      histos.fill(HIST("Reco/hNchVsFT0C"), ampFT0C, nch[iMc], category);
      histos.fill(HIST("Reco/hNchFromVertexVsFT0C"), ampFT0C, nchFromVertex[iMc], category);
      histos.fill(HIST("Reco/hNchFT0CFromVertexVsFT0C"), ampFT0C, nchFT0CFromVertex[iMc], category);
      histos.fill(HIST("Reco/hNPVVsFT0C"), ampFT0C, nPV, category);
      histos.fill(HIST("Reco/hNGlobalVsFT0C"), ampFT0C, nGlobal, category);
      histos.fill(HIST("Reco/hNGlobalVsNchFromVertex"), nchFromVertex[iMc], nGlobal, category);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<McEventValidationTask>(cfgc)};
}
