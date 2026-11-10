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

/// \file match-mft-mch-data-mc.cxx
/// \brief MC entry point and truth annotation for shared MFT MCH matching.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

// o2-linter: disable=name/workflow-file (Thin entry point delegates to the shared workflow factory.)
// o2-linter: disable=name/file-cpp (Preserve the established matching task source path.)

#include "Common/Core/MftMchMatching.h"
#include "Common/DataModel/MatchMFTMuonData.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
void customize(std::vector<o2::framework::ConfigParamSpec>& options)
{
  options.push_back({"matching-output-mode", o2::framework::VariantType::String, "auto", {"auto reads matching-output.cfgOutputMode from local JSON; tree, sparse or both overrides it"}});
}

#include <Framework/runDataProcessing.h>

#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

// This device only annotates final residual rows. No matching/propagation here.
template <bool Sparse, bool InternalPopulation = false>
struct MftMchMcTruth {
  template <typename T>
  using Input = o2::analysis::mftmch::MatchingAuditInput<T, Sparse, InternalPopulation ? 1 : 0>;
  Produces<aod::MftMchMatchTruths> truths;
  o2::analysis::mftmch::MatchingAuditOutput<aod::MatchingCandidateTruths> candidateTruths;
  o2::analysis::mftmch::MatchingAuditOutput<aod::MatchingCandidateRelations> candidateRelations;
  o2::analysis::mftmch::MatchingAuditOutput<aod::MatchingMchTruths> mchTruths;
  Produces<aod::MixedMatchPlaneResiduals> mixedResiduals;
  o2::analysis::mftmch::MatchingAuditOutput<aod::MixedMftMchMatchTruths> mixedTruths;
  Produces<aod::MftMchMatchRelations> relations;
  o2::analysis::mftmch::MixedPairSampler mixedSampler;
  OutputObj<TNamed> samplingMetadata{TNamed("MixedPairSampling", "")};
  Configurable<int> cfgMcAncestryMaxDepth{"cfgMcAncestryMaxDepth", 128, "Maximum stored MC mother-graph depth; truncation is explicitly flagged"};
  Configurable<int> cfgMcAncestryMaxNodes{"cfgMcAncestryMaxNodes", 4096, "Maximum mother-graph states expanded per distinct MCH MC particle"};
  void init(InitContext&)
  {
    samplingMetadata->SetTitle(mixedSampler.metadata().c_str());
    if (cfgMcAncestryMaxDepth < 1 || cfgMcAncestryMaxNodes < 1) {
      LOGF(fatal, "MC ancestry limits must be positive");
    }
  }
  using MchLabels = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
  using MftLabels = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  struct ParticleInfo {
    int64_t id = -1, collision = -1;
    int pdg = 0;
    bool valid = false, primary = false;
    float px = std::numeric_limits<float>::quiet_NaN();
    float py = std::numeric_limits<float>::quiet_NaN();
    float pz = std::numeric_limits<float>::quiet_NaN();
  };
  static ParticleInfo particle(int64_t id, aod::McParticles const& particles)
  {
    ParticleInfo out;
    out.id = id;
    const auto offset = static_cast<int64_t>(particles.offset());
    if (id < offset || id - offset >= particles.size()) {
      return out;
    }
    const auto p = particles.rawIteratorAt(id - offset);
    out.valid = true;
    out.collision = p.mcCollisionId();
    out.pdg = p.pdgCode();
    out.px = p.px();
    out.py = p.py();
    out.pz = p.pz();
    out.primary = p.isPhysicalPrimary();
    return out;
  }
  template <typename Rows, typename TruthOutput, typename RelationOutput>
  void annotate(Rows const& residuals, MchLabels const& mchTracks,
                MftLabels const& mftTracks, aod::McParticles const& particles,
                TruthOutput& truthOutput, RelationOutput& relationOutput)
  {
    using namespace o2::analysis::mftmch::truth;
    // Local caches must not cross DF boundaries: MC IDs are DF scoped.
    std::unordered_map<int64_t, Node> nodes;
    std::unordered_map<int64_t, Search> searches;
    auto node = [&](int64_t id) -> Node {
      auto cached = nodes.find(id);
      if (cached != nodes.end()) {
        return cached->second;
      }
      Node n;
      n.id = id;
      const auto offset = static_cast<int64_t>(particles.offset());
      if (id >= offset && id - offset < particles.size()) {
        const auto p = particles.rawIteratorAt(id - offset);
        n.valid = true;
        n.pdg = p.pdgCode();
        n.process = p.getProcess();
        n.generator = p.producedByGenerator();
        n.vx = p.vx();
        n.vy = p.vy();
        n.vz = p.vz();
        for (auto const& mother : p.mothersIds()) {
          n.mothers.push_back(mother);
        }
      }
      nodes.emplace(id, n);
      return n;
    };
    auto motherPdgs = [&](const Node& n) {
      std::vector<int> pdgs;
      for (auto const& id : n.mothers) {
        auto p = node(id);
        pdgs.push_back(p.valid ? p.pdg : 0);
      }
      return pdgs;
    };
    int64_t residualRow = residuals.offset();
    for (const auto& row : residuals) {
      const auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      const auto mft = mftTracks.rawIteratorAt(row.mftTrackId() - mftTracks.offset());
      const auto a = particle(mch.mcParticleId(), particles);
      const auto b = particle(mft.mcParticleId(), particles);
      // Strict truth policy is annotation-only; it never changes matching.
      const bool labelsValid = a.valid && b.valid;
      const bool fake = hasFakeLabel(mch.mcMask(), mft.mcMask());
      const auto matchStatus = matchingTruthStatus(a.valid, b.valid, a.id, b.id, mch.mcMask(), mft.mcMask());
      const bool correct = matchStatus == MatchStatus::Correct;
      const bool falseMatch = matchStatus == MatchStatus::False;
      const auto status = static_cast<int8_t>(matchStatus);
      truthOutput(residualRow, row.mchTrackId(), row.mftTrackId(), a.id, b.id,
                  a.valid, b.valid, correct, mch.mcMask(), mft.mcMask(),
                  a.collision, b.collision, a.pdg, b.pdg,
                  a.px, a.py, a.pz, b.px, b.py, b.pz, a.primary, b.primary, falseMatch, status);
      const Node mchNode = node(a.id), mftNode = node(b.id);
      const Search* search = nullptr;
      if (labelsValid && !correct && (a.pdg == kMuonMinus || a.pdg == -kMuonMinus)) {
        auto cached = searches.find(a.id);
        if (cached == searches.end()) {
          cached = searches.emplace(a.id, findAncestors(mchNode, node, static_cast<int>(TMCProcess::kPDecay), cfgMcAncestryMaxDepth, cfgMcAncestryMaxNodes)).first;
        }
        search = &cached->second;
      }
      auto relation = classify(mchNode, mftNode, search);
      if (fake) {
        // Keep genealogy/quality columns, but do not call a fake pair correct.
        relation.origin = 4;
      }
      std::vector<int64_t> pathIds;
      std::vector<int> pathPdgs, pathProcesses;
      std::vector<uint8_t> pathGenerator;
      std::vector<float> pathX, pathY, pathZ;
      for (const auto& n : relation.parentToChild) {
        pathIds.push_back(n.id);
        pathPdgs.push_back(n.pdg);
        pathProcesses.push_back(n.process);
        pathGenerator.push_back(static_cast<uint8_t>(n.generator));
        pathX.push_back(n.vx);
        pathY.push_back(n.vy);
        pathZ.push_back(n.vz);
      }
      const float nan = std::numeric_limits<float>::quiet_NaN();
      relationOutput(residualRow++, row.mchTrackId(), row.mftTrackId(),
                     static_cast<int8_t>(relation.origin), relation.status, relation.status == Complete,
                     a.valid && (a.pdg == kMuonMinus || a.pdg == -kMuonMinus), relation.ancestor, relation.decay,
                     relation.decay && relation.depth == 1, relation.depth,
                     mchNode.process, mftNode.process, mchNode.generator, mftNode.generator,
                     mchNode.valid ? mchNode.vx : nan, mchNode.valid ? mchNode.vy : nan, mchNode.valid ? mchNode.vz : nan,
                     mftNode.valid ? mftNode.vx : nan, mftNode.valid ? mftNode.vy : nan, mftNode.valid ? mftNode.vz : nan,
                     mchNode.mothers, motherPdgs(mchNode), mftNode.mothers, motherPdgs(mftNode),
                     pathIds, pathPdgs, pathProcesses, pathGenerator, pathX, pathY, pathZ);
    }
  }
  void process(aod::MatchPlaneResiduals const& residuals, MchLabels const& mchTracks,
               MftLabels const& mftTracks, aod::McParticles const& particles)
  {
    annotate(residuals, mchTracks, mftTracks, particles, truths, relations);
  }
  void processCandidates(Input<aod::MatchingCandidates> const& candidates, Input<aod::MatchingMchPopulation> const& population,
                         MchLabels const& mchTracks, MftLabels const& mftTracks, aod::McParticles const& particles)
  {
    using namespace o2::analysis::mftmch::truth;
    annotate(candidates, mchTracks, mftTracks, particles, candidateTruths, candidateRelations);
    std::unordered_map<int64_t, int> trueCounts;
    std::set<std::pair<int64_t, int64_t>> seen;
    for (auto const& row : candidates) {
      if (!seen.emplace(row.mchTrackId(), row.mftTrackId()).second) {
        LOGF(fatal, "Duplicate candidate pair would corrupt N_pairable");
      }
      auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      auto mft = mftTracks.rawIteratorAt(row.mftTrackId() - mftTracks.offset());
      auto a = particle(mch.mcParticleId(), particles), b = particle(mft.mcParticleId(), particles);
      if (matchingTruthStatus(a.valid, b.valid, a.id, b.id, mch.mcMask(), mft.mcMask()) == MatchStatus::Correct) {
        ++trueCounts[row.mchTrackId()];
      }
    }
    auto getNode = [&](int64_t id) {
      Node n;
      const auto info = particle(id, particles);
      n.id = id;
      n.pdg = info.pdg;
      n.valid = info.valid;
      if (info.valid) {
        auto p = particles.rawIteratorAt(id - particles.offset());
        n.process = p.getProcess();
        n.generator = p.producedByGenerator();
        n.vx = p.vx();
        n.vy = p.vy();
        n.vz = p.vz();
        for (auto const& mother : p.mothersIds()) {
          n.mothers.push_back(mother);
        }
      }
      return n;
    };
    int64_t populationRow = population.offset();
    for (auto const& row : population) {
      auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      auto info = particle(mch.mcParticleId(), particles);
      std::vector<int64_t> ancestorIds;
      std::vector<int> ancestorPdgs;
      int status = info.valid ? Complete : MissingLabel;
      auto hfOrigin = o2::analysis::mftmch::MuonHfOrigin::Unknown;
      if (info.valid && std::abs(info.pdg) == kMuonMinus) {
        auto search = findAncestors(getNode(info.id), getNode, static_cast<int>(TMCProcess::kPDecay), cfgMcAncestryMaxDepth, cfgMcAncestryMaxNodes);
        status = search.status;
        hfOrigin = classifyHfOrigin(getNode(info.id), search);
        for (auto const& [id, path] : search.ancestors) {
          auto ancestor = getNode(id);
          if (std::abs(ancestor.pdg) == kPiPlus || std::abs(ancestor.pdg) == kKPlus) {
            ancestorIds.push_back(id);
            ancestorPdgs.push_back(ancestor.pdg);
          }
        }
      }
      mchTruths(populationRow++, row.mchTrackId(), info.id, mch.mcMask(), info.valid, info.pdg, info.collision,
                trueCounts[row.mchTrackId()], !ancestorIds.empty(), status, ancestorIds, ancestorPdgs,
                static_cast<int8_t>(classifySignalOrigin(getNode(info.id), getNode, cfgMcAncestryMaxDepth, cfgMcAncestryMaxNodes)), static_cast<int8_t>(hfOrigin));
    }
  }
  PROCESS_SWITCH(MftMchMcTruth, processCandidates, "Truth for ALL candidate identities and all MCHs before matching selection", true);
  // Only MC-independent matching is done upstream. This publication guard makes
  // the MC mixed sample strictly different generator collisions and particles.
  // Labels are never forced false to hide a true/unknown pair.
  struct CachedMftTruth {
    ParticleInfo info;
    uint16_t mask = 0;
  };
  // MC particle/collision IDs are only comparable within the same source frame.
  std::map<uint64_t, std::unordered_map<int64_t, CachedMftTruth>> mixedLabelFrames;
  void processMixed(aod::RawMixedMatchPlaneResiduals const& rows, MchLabels const& mchTracks,
                    MftLabels const& mftTracks, aod::McParticles const& particles,
                    aod::MftMchMixingPoolStates const& poolStates)
  {
    if (poolStates.size() == 0) {
      if (rows.size() != 0) {
        LOGF(fatal, "Missing mixed-event source-frame metadata");
      }
      mixedLabelFrames.clear();
      return;
    }
    if (poolStates.size() != 1) {
      LOGF(fatal, "Expected one pool state per DF");
    }
    const auto state = poolStates.begin();
    const auto frame = state.targetFrameId();
    auto& current = mixedLabelFrames[frame];
    current.clear();
    for (const auto& mft : mftTracks) {
      current.emplace(mft.globalIndex(), CachedMftTruth{.info = particle(mft.mcParticleId(), particles), .mask = mft.mcMask()});
    }
    auto donorTruth = [&](const auto& row) -> const CachedMftTruth& {
      if (row.targetFrameId() != frame || row.mixingStreamId() != state.mixingStreamId()) {
        LOGF(fatal, "Mixed source frame/stream mismatch");
      }
      const auto source = mixedLabelFrames.find(row.donorFrameId());
      if (source == mixedLabelFrames.end()) {
        LOGF(fatal, "Historical MFT label frame is unavailable");
      }
      const auto track = source->second.find(row.mftTrackId());
      if (track == source->second.end()) {
        LOGF(fatal, "Historical MFT label ID is unavailable");
      }
      return track->second;
    };
    std::set<int64_t> eligible;
    std::map<int64_t, std::set<int>> groups;
    size_t rejectedUnknown = 0, rejectedSameMcEvent = 0;
    int64_t sourceIndex = 0;
    for (const auto& row : rows) {
      const auto sourceRow = sourceIndex++;
      const auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      const auto a = particle(mch.mcParticleId(), particles);
      const auto& b = donorTruth(row).info;
      if (!a.valid || !b.valid || a.collision < 0 || b.collision < 0) {
        ++rejectedUnknown;
        continue;
      }
      if (row.targetFrameId() == row.donorFrameId() && (a.id == b.id || a.collision == b.collision)) {
        ++rejectedSameMcEvent;
        continue;
      }
      eligible.insert(sourceRow);
      groups[row.mchTrackId()].insert(row.mixedGroupIndex());
    }
    int64_t index = 0;
    sourceIndex = 0;
    mixedSampler.decisions.clear();
    for (const auto& row : rows) {
      if (!eligible.contains(sourceIndex++)) {
        continue;
      }
      if (!mixedSampler.keep(row)) {
        continue;
      }
      const auto mch = mchTracks.rawIteratorAt(row.mchTrackId() - mchTracks.offset());
      const auto a = particle(mch.mcParticleId(), particles);
      const auto& donor = donorTruth(row);
      const auto& b = donor.info;
      o2::analysis::mftmch::copyMixedResidual(row, mixedResiduals, groups.at(row.mchTrackId()).size());
      mixedTruths(index++, row.mchTrackId(), row.mftTrackId(), a.id, b.id,
                  true, true, false, mch.mcMask(), donor.mask, a.collision, b.collision, a.pdg, b.pdg,
                  a.px, a.py, a.pz, b.px, b.py, b.pz, a.primary, b.primary, true, static_cast<int8_t>(0));
    }
    samplingMetadata->SetTitle(mixedSampler.metadata().c_str());
    // Evicted upstream frames may still have been used in the current rows.
    // Prune only after publication, against the exact retained-frame set.
    const auto retained = state.retainedFrames();
    std::set<uint64_t> keep(retained.begin(), retained.end());
    for (auto it = mixedLabelFrames.begin(); it != mixedLabelFrames.end();) {
      if (!keep.contains(it->first)) {
        it = mixedLabelFrames.erase(it);
      } else {
        ++it;
      }
    }
    LOGF(info, "Mixed MC publication: saved=%" PRId64 ", unknown rejected=%zu, same-MC-event/particle rejected=%zu, retained label frames=%zu",
         static_cast<int64_t>(index), rejectedUnknown, rejectedSameMcEvent, mixedLabelFrames.size());
  }
  PROCESS_SWITCH(MftMchMcTruth, processMixed, "Publish mixed rows only for valid, distinct MC collisions/particles", true);
};

// This device adds a separate MC-annotated output. The shared reconstruction,
// inclusive dimuons and pre-selection MCH denominator remain independent of MC.
template <bool Sparse>
struct McDimuonTruth {
  template <typename T>
  using Input = o2::analysis::mftmch::MatchingAuditInput<T, Sparse>;
  o2::analysis::mftmch::MatchingAuditOutput<aod::McDimuonTruths> pairs;
  Configurable<bool> cfgEnableMcDimuonTruth{"cfgEnableMcDimuonTruth", false, "Save MC-annotated MCH-PV and selected chi2/ML dimuons; matched methods require cfgSaveSelectedMuonTable"};
  Configurable<int> cfgMcAncestryMaxDepth{"cfgMcAncestryMaxDepth", 128, "Maximum MC ancestor depth; incomplete ancestry is flagged without rejecting the dimuon"};
  Configurable<int> cfgMcAncestryMaxNodes{"cfgMcAncestryMaxNodes", 4096, "Maximum distinct MC ancestor nodes per muon"};
  Configurable<std::vector<double>> cfgDimuonMassAxis{"cfgDimuonMassAxis", {600, 0., 12.}, "Sparse mass: nbins,min,max GeV/c^2"};
  Configurable<std::vector<double>> cfgDimuonPtAxis{"cfgDimuonPtAxis", {200, 0., 20.}, "Sparse dimuon pT: nbins,min,max GeV/c"};
  Configurable<std::vector<double>> cfgDimuonRapidityAxis{"cfgDimuonRapidityAxis", {80, -5., -1.}, "Sparse dimuon rapidity: nbins,min,max"};
  static AxisSpec axis(std::vector<double> const& v, char const* title)
  {
    if (v.size() != 3 || !std::isfinite(v[0]) || v[0] < 1 || v[0] > 1000000 || v[0] != std::floor(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]) || v[1] >= v[2]) {
      LOGF(fatal, "Invalid MC dimuon axis: require positive integer nbins,min,max");
    }
    return {static_cast<int>(v[0]), v[1], v[2], title};
  }
  void init(InitContext&)
  {
    if (!cfgEnableMcDimuonTruth) {
      return;
    }
    if (cfgMcAncestryMaxDepth < 1 || cfgMcAncestryMaxNodes < 1) {
      LOGF(fatal, "MC dimuon ancestry limits must be positive");
    }
  }
  struct Ancestry {
    std::set<int64_t> jpsis;
    bool complete = true;
  };
  void process(Input<aod::MchDimuonsPV> const& mchPairs, Input<aod::SelectedMatchingMuons> const& muons, Input<aod::MatchingDimuonVertices> const& vertices,
               typename MftMchMcTruth<Sparse>::MchLabels const& mchTracks, typename MftMchMcTruth<Sparse>::MftLabels const& mftTracks,
               aod::McParticles const& particles)
  {
    if (!cfgEnableMcDimuonTruth) {
      return;
    }
    std::map<std::tuple<int64_t, int8_t, int64_t, int64_t>, bool> vertexCuts;
    for (auto const& vertex : vertices) {
      vertexCuts[{vertex.collisionId(), vertex.method(), std::min(vertex.mchTrackId1(), vertex.mchTrackId2()), std::max(vertex.mchTrackId1(), vertex.mchTrackId2())}] = vertex.passesVertexCuts();
    }
    auto inRange = [](int64_t id, auto const& table) { return id >= static_cast<int64_t>(table.offset()) && id - static_cast<int64_t>(table.offset()) < table.size(); };
    auto particle = [&](int64_t id) {
      typename MftMchMcTruth<Sparse>::ParticleInfo info;
      info.id = id;
      if (inRange(id, particles)) {
        auto p = particles.rawIteratorAt(id - particles.offset());
        info.valid = true;
        info.pdg = p.pdgCode();
        info.collision = p.mcCollisionId();
      }
      return info;
    };
    auto label = [&](int64_t id, auto const& tracks) -> int64_t {
      if (!inRange(id, tracks)) {
        return -1;
      }
      auto t = tracks.rawIteratorAt(id - tracks.offset());
      if ((t.mcMask() & 128) != 0 || !particle(t.mcParticleId()).valid) {
        return -1;
      }
      return t.mcParticleId();
    };
    // MC IDs and caches are local to one process call / DF.
    std::map<int64_t, Ancestry> cache;
    auto ancestors = [&](int64_t id) -> Ancestry const& {
      auto found = cache.find(id);
      if (found != cache.end()) {
        return found->second;
      }
      Ancestry result;
      std::vector<std::pair<int64_t, int>> todo{{id, 0}};
      std::set<int64_t> seen;
      for (size_t i = 0; i < todo.size(); ++i) {
        auto [current, depth] = todo[i];
        if (!seen.insert(current).second) {
          continue;
        }
        if (seen.size() > static_cast<size_t>(cfgMcAncestryMaxNodes.value) || !inRange(current, particles)) {
          result.complete = false;
          break;
        }
        auto p = particles.rawIteratorAt(current - particles.offset());
        if (p.pdgCode() == 443) {
          result.jpsis.insert(current);
        }
        for (auto mother : p.mothersIds()) {
          if (depth >= cfgMcAncestryMaxDepth) {
            result.complete = false;
          } else {
            todo.emplace_back(mother, depth + 1);
          }
        }
      }
      return cache.emplace(id, std::move(result)).first->second;
    };
    auto save = [&](int64_t collision, int8_t method, int64_t mch1, int64_t mft1, int64_t mch2, int64_t mft2,
                    int8_t charge, double score1, double score2, double mass, double pt, double rapidity) {
      if (!std::isfinite(mass) || !std::isfinite(pt) || !std::isfinite(rapidity)) {
        return;
      }
      if (vertices.size() != 0) {
        auto vertex = vertexCuts.find({collision, method, std::min(mch1, mch2), std::max(mch1, mch2)});
        if (vertex == vertexCuts.end() || !vertex->second) {
          return;
        }
      }
      const auto id1 = label(mch1, mchTracks), id2 = label(mch2, mchTracks);
      const auto p1 = particle(id1), p2 = particle(id2);
      const bool oppositeMuons = p1.valid && p2.valid && std::abs(p1.pdg) == 13 && p1.pdg == -p2.pdg;
      bool complete = false;
      int64_t jpsi = -1;
      if (oppositeMuons && p1.collision >= 0 && p1.collision == p2.collision) {
        auto const& a = ancestors(id1);
        auto const& b = ancestors(id2);
        complete = a.complete && b.complete;
        if (complete) {
          for (auto ancestor : a.jpsis) {
            if (b.jpsis.contains(ancestor)) {
              jpsi = ancestor;
              break;
            }
          }
        }
      }
      const bool commonJpsi = jpsi >= 0;
      const int8_t correct = method == 0 ? -1 : static_cast<int8_t>((id1 >= 0 && id1 == label(mft1, mftTracks)) + (id2 >= 0 && id2 == label(mft2, mftTracks)));
      std::array<std::vector<int64_t>, 2> motherIds;
      std::array<std::vector<int>, 2> motherPdgs;
      std::vector<int64_t> commonMotherIds;
      std::vector<int> commonMotherPdgs;
      bool mothersComplete = charge == 0 && p1.valid && p2.valid;
      if (charge == 0) {
        const std::array<int64_t, 2> childIds{id1, id2};
        for (size_t leg = 0; leg < childIds.size(); ++leg) {
          if (!inRange(childIds[leg], particles)) {
            continue;
          }
          const auto child = particles.rawIteratorAt(childIds[leg] - particles.offset());
          for (auto mother : child.mothersIds()) {
            motherIds[leg].push_back(mother);
            const auto info = particle(mother);
            motherPdgs[leg].push_back(info.valid ? info.pdg : 0);
            mothersComplete = mothersComplete && info.valid;
          }
        }
        // Compare particle IDs, never just PDG codes. Preserve multiple mothers.
        std::set<int64_t> shared;
        for (auto mother : motherIds[0]) {
          if (inRange(mother, particles) && std::find(motherIds[1].begin(), motherIds[1].end(), mother) != motherIds[1].end()) {
            shared.insert(mother);
          }
        }
        for (auto mother : shared) {
          commonMotherIds.push_back(mother);
          commonMotherPdgs.push_back(particle(mother).pdg);
        }
      }

      uint16_t resonanceMask = 0;
      std::vector<int64_t> resonanceIds;
      std::vector<int> resonancePdgs;
      // The tag depends ONLY on the MCH particles, not on selected MFT labels.
      // A pion/kaon at either MCH leg or a merely common remote ancestor fails.
      if (charge == 0 && mothersComplete && std::abs(p1.pdg) == 13 && std::abs(p2.pdg) == 13) {
        for (size_t i = 0; i < commonMotherIds.size(); ++i) {
          const auto bit = aod::mcdimuontruth::resonanceMotherBit(commonMotherPdgs[i]);
          if (bit != 0) {
            resonanceMask |= bit;
            resonanceIds.push_back(commonMotherIds[i]);
            resonancePdgs.push_back(commonMotherPdgs[i]);
          }
        }
      }
      pairs(collision, method, mch1, mft1, mch2, mft2, id1, id2, jpsi, oppositeMuons, commonJpsi, complete,
            correct, charge, score1, score2, mass, pt, rapidity,
            motherIds[0], motherPdgs[0], motherIds[1], motherPdgs[1], commonMotherIds, commonMotherPdgs, mothersComplete, resonanceMask, resonanceIds, resonancePdgs);
    };
    // MCH-only Same/Mixed dimuons are written once by mch-dimuon-output.
    (void)mchPairs;
    using MuonRow = typename Input<aod::SelectedMatchingMuons>::iterator;
    std::map<std::tuple<uint64_t, int64_t, int8_t>, std::vector<MuonRow>> groups;
    for (auto const& muon : muons) {
      if (muon.method() == 1 || muon.method() == 2) {
        groups[{muon.frameId(), muon.collisionId(), muon.method()}].push_back(muon);
      }
    }
    for (auto const& [key, rows] : groups) {
      for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = i + 1; j < rows.size(); ++j) {
          auto const& a = rows[i];
          auto const& b = rows[j];
          if (a.mchTrackId() == b.mchTrackId() || a.mftTrackId() == b.mftTrackId()) {
            continue;
          }
          auto p = o2::analysis::mftmch::muonMomentum(a.px(), a.py(), a.pz(), o2::constants::physics::MassMuon) + o2::analysis::mftmch::muonMomentum(b.px(), b.py(), b.pz(), o2::constants::physics::MassMuon);
          save(std::get<1>(key), std::get<2>(key), a.mchTrackId(), a.mftTrackId(), b.mchTrackId(), b.mftTrackId(),
               static_cast<int8_t>(a.charge() + b.charge()), a.score(), b.score(), p.M(), p.Pt(), p.Rapidity());
        }
      }
    }
  }
};

// Consume optional legacy MC audit transport without duplicating sparse output.
// The unified GlobalMuon dimuon sparse is filled independently of these tables.
template <bool Sparse, uint8_t Mask = 0>
struct MatchingMcOutputSparse {
  template <typename T>
  using Input = o2::analysis::mftmch::MatchingAuditInput<T, Sparse, Mask>;
  void process(Input<aod::McDimuonTruths> const&, Input<aod::MatchingDimuonVertices> const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  auto workflow = o2::analysis::mftmch::makeMatchingWorkflow<true>(context, true);
  const auto mode = o2::analysis::mftmch::matchingOutputMode(context);
  if (mode == "sparse") {
    workflow.push_back(adaptAnalysisTask<MftMchMcTruth<true>>(context, TaskName{"mft-mch-mc-truth"}));
    workflow.push_back(adaptAnalysisTask<McDimuonTruth<true>>(context, TaskName{"mc-dimuon-truth"}));
    workflow.push_back(adaptAnalysisTask<MatchingMcOutputSparse<true>>(context, TaskName{"matching-mc-output"}));
  } else {
    if (o2::analysis::mftmch::matchingTreeStorage.mchPopulation) {
      workflow.push_back(adaptAnalysisTask<MftMchMcTruth<false>>(context, TaskName{"mft-mch-mc-truth"}));
    } else {
      workflow.push_back(adaptAnalysisTask<MftMchMcTruth<false, true>>(context, TaskName{"mft-mch-mc-truth"}));
    }
    workflow.push_back(adaptAnalysisTask<McDimuonTruth<false>>(context, TaskName{"mc-dimuon-truth"}));
    if (mode == "both") {
      workflow.push_back(o2::analysis::mftmch::withMatchingAuditStorageMask(o2::analysis::mftmch::matchingAuditStorageMask(), [&]<uint8_t Mask>() {
        return adaptAnalysisTask<MatchingMcOutputSparse<false, Mask>>(context, TaskName{"matching-mc-output"});
      }));
    }
  }
  o2::analysis::mftmch::configureMatchingWriter(context, workflow, mode);
  return workflow;
}
