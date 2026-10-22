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
///
/// \file omega2012TrainingTable.cxx
/// \brief Derived Omega(2012) training tables, separately for the XiK0s and Xi1530K decay modes
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// Mode A (XiK0s): processXiK0s, processXiK0sMC. Mode B (Xi1530K, Xi- pi+ K- with micro tracks):
/// processXi1530KMicro, processXi1530KMCMicro. Generated parents: processMCTrue. MC truth of reconstructed MC: processMcTruth.
/// Each candidate is written only to the tables of its own mode (PWGLF/DataModel/LFOmega2012MlTables.h);
/// the selection is the one of omega2012Analysis.cxx (PWGLF/Core/Omega2012AnalysisCore.h, same JSON keys).
/// Only collisions with written candidates get an event row. In reconstructed MC the daughter labels of the written
/// objects and the pruned McParticles rows (Omega(2012) decay trees, ancestry of the labels) carry the truth, so the
/// output needs no source table.

#include "PWGLF/Core/Omega2012AnalysisCore.h"
#include "PWGLF/Core/Omega2012MlFeatures.h"
#include "PWGLF/Core/ResoAnalysisSelectionCore.h"
#include "PWGLF/DataModel/LFOmega2012MlTables.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis::omega2012;

namespace ml = o2::analysis::omega2012ml;

static_assert(std::extent_v<aod::omega2012ml::OmXiK0sFeatures::type> == ml::XiK0sNMasterFeatures,
              "OmXiK0sFeatures column size must match the XiK0s feature contract");
static_assert(std::extent_v<aod::omega2012ml::OmXi1530KFeatures::type> == ml::Xi1530KNMasterFeatures,
              "OmXi1530KFeatures column size must match the Xi1530K feature contract");

struct Omega2012TrainingTable {
  using ResoCollisions = aod::ResoCollisions_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001, aod::ResoMCCollisionIds>;
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCMicroTracks = soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001>;
  using ResoMCCascades = soa::Join<aod::ResoCascades, aod::ResoMCCascades>;
  using ResoMCV0s = soa::Join<aod::ResoV0s, aod::ResoMCV0s>;
  using ResoMCParents = aod::ResoMCParents_001;

  // FNV-1a over the bit patterns of the master features, for parity logs
  static constexpr uint64_t FnvOffsetBasis = 14695981039346656037ULL;
  static constexpr uint64_t FnvPrime = 1099511628211ULL;
  static constexpr unsigned int BitsPerByte = 8;
  static constexpr unsigned int BitsPerFloat = 32;
  static constexpr uint32_t ByteMask = 0xffU;
  static constexpr int NBuildStatus = 6;
  // Pruned McParticles: generations of Omega(2012) descendants and mother steps above a daughter label
  static constexpr int McTreeDepth = 3;
  static constexpr int MaxAncestorSteps = 8;
  static constexpr uint8_t McSelectionTree = 1;
  static constexpr uint8_t McSelectionAncestry = 2;
  static constexpr std::size_t NStoredMothers = 2;

  Produces<aod::Omega2012MlEvents> mlEvents;
  Produces<aod::Omega2012MlCascades> mlCascades;
  Produces<aod::Omega2012MlV0s> mlV0s;
  Produces<aod::Omega2012MlTracks> mlTracks;
  Produces<aod::Omega2012MlXiK0sCandidates> mlXiK0sCandidates;
  Produces<aod::Omega2012MlXiK0sInputs> mlXiK0sInputs;
  Produces<aod::Omega2012MlXiK0sTruth> mlXiK0sTruth;
  Produces<aod::Omega2012MlXi1530KCandidates> mlXi1530KCandidates;
  Produces<aod::Omega2012MlXi1530KInputs> mlXi1530KInputs;
  Produces<aod::Omega2012MlXi1530KTruth> mlXi1530KTruth;
  Produces<aod::Omega2012MlGenAudit> mlGenAudit;
  Produces<aod::Omega2012MlEventMcs> mlEventMcs;
  Produces<aod::Omega2012MlCascadeMcs> mlCascadeMcs;
  Produces<aod::Omega2012MlV0Mcs> mlV0Mcs;
  Produces<aod::Omega2012MlTrackMcs> mlTrackMcs;
  Produces<aod::Omega2012MlMcParticles> mlMcParticles;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Selection shared with the Omega(2012) histogram task (same JSON keys)
  o2::analysis::resonance::EventCuts eventCuts;
  XiCuts xiCuts;
  K0sCuts k0sCuts;
  Xi1530KCuts xi1530KCuts;
  o2::analysis::resonance::TrackCuts trackCuts = makeXi1530KTrackCuts();
  PionPidCuts pionPID;
  KaonPidCuts kaonPID;
  CandidateCuts candidateCuts;

  Configurable<std::string> omega2012MlExportStage{"omega2012MlExportStage", "loose", "Candidate export stage: loose (pass bits >= 1) or selected (all pass bits of the mode)"};
  Configurable<bool> omega2012MlLooseAudit{"omega2012MlLooseAudit", true, "Record the loose cut flow and mass/pT/activity spectrum per mode"};
  Configurable<int> omega2012MlParityRows{"omega2012MlParityRows", 0, "Log feature hashes of the first N exported candidates of each mode"};
  Configurable<bool> omega2012MlExportWrongSign{"omega2012MlExportWrongSign", false, "Also export the wrong-sign charge-pattern controls of the Xi1530K mode"};

  // Event row written before the first candidate of the collision
  struct PendingEvent {
    DecayMode mode = DecayMode::XiK0s;
    int64_t recoCollision = -1;
    float posZ = 0.f;
    float bField = 0.f;
    float centrality = 0.f;
    float multiplicity = 0.f;
    bool recINELgt0 = false;
    bool isMC = false;
    int64_t mcCollision = -1;
  };

  Omega2012AnalysisCore core;
  PendingEvent pendingEvent;
  int64_t mlEventRow = -1;
  // MC tables of the current DF, set only by the reconstructed-MC process functions
  aod::McParticles const* mcParticlesDf = nullptr;
  aod::McTrackLabels const* mcTrackLabelsDf = nullptr;
  std::unordered_set<int64_t> mlAncestryRows;
  std::unordered_map<int64_t, int64_t> mlCascadeRows;
  std::unordered_map<int64_t, int64_t> mlV0Rows;
  std::unordered_map<int64_t, int64_t> mlTrackRows;
  int parityLoggedXiK0s = 0;
  int parityLoggedXi1530K = 0;

  void init(InitContext&)
  {
    const bool xiK0s = doprocessXiK0s || doprocessXiK0sMC;
    const bool xi1530K = doprocessXi1530KMicro || doprocessXi1530KMCMicro;
    if ((doprocessXiK0s && doprocessXiK0sMC) || (doprocessXi1530KMicro && doprocessXi1530KMCMicro)) {
      LOG(fatal) << "Enable at most one of processXiK0s/processXiK0sMC and one of processXi1530KMicro/processXi1530KMCMicro";
    }
    if (!xiK0s && !xi1530K) {
      LOG(fatal) << "Enable at least one candidate process function";
    }
    if (doprocessMCTrue && !doprocessXiK0sMC && !doprocessXi1530KMCMicro) {
      LOG(fatal) << "processMCTrue requires processXiK0sMC or processXi1530KMCMicro";
    }
    if (doprocessMcTruth != (doprocessXiK0sMC || doprocessXi1530KMCMicro)) {
      LOG(fatal) << "processMcTruth must be enabled exactly when processXiK0sMC or processXi1530KMCMicro is enabled";
    }
    if (omega2012MlParityRows < 0) {
      LOG(fatal) << "omega2012MlParityRows must not be negative";
    }
    if (omega2012MlExportStage.value != "loose" && omega2012MlExportStage.value != "selected") {
      LOG(fatal) << "omega2012MlExportStage must be loose or selected";
    }

    ProcessModes modes;
    modes.xiK0s = xiK0s;
    modes.xi1530K = xi1530K;
    modes.microTracks = xi1530K;
    modes.mcReco = doprocessXiK0sMC || doprocessXi1530KMCMicro;
    modes.mcGen = doprocessMCTrue;
    LooseStageOptions looseOptions;
    looseOptions.audit = omega2012MlLooseAudit;
    looseOptions.exportSelected = omega2012MlExportStage.value == "selected";
    core.init(histos, eventCuts, xiCuts, k0sCuts, xi1530KCuts, trackCuts, pionPID, kaonPID, candidateCuts, modes, looseOptions);

    // Candidates that violate the canonical or feature contract are skipped, never written.
    const std::array<const char*, NBuildStatus> statusLabels{"Ok", "InvalidChargePattern", "ReusedDaughter", "InvalidMomentum", "InvalidKinematics", "InvalidContract"};
    if (xiK0s) {
      auto skipped = histos.add<TH1>("ML/xiK0s/exportSkipped", "XiK0s mode skipped candidates;build status;candidates", HistType::kTH1D, {{NBuildStatus, -0.5, NBuildStatus - 0.5}});
      for (std::size_t i = 0; i < statusLabels.size(); ++i) {
        skipped->GetXaxis()->SetBinLabel(i + 1, statusLabels[i]);
      }
    }
    if (xi1530K) {
      auto skipped = histos.add<TH1>("ML/xi1530K/exportSkipped", "Xi1530K mode skipped candidates;build status;candidates", HistType::kTH1D, {{NBuildStatus, -0.5, NBuildStatus - 0.5}});
      for (std::size_t i = 0; i < statusLabels.size(); ++i) {
        skipped->GetXaxis()->SetBinLabel(i + 1, statusLabels[i]);
      }
      histos.add("ML/xi1530K/chargePattern", "Xi1530K mode candidates handed to the export;charge pattern (0 signal);candidates", HistType::kTH1D, {{4, -0.5, 3.5}});
    }

    LOG(info) << "Size of the histograms in Omega(2012) training table task";
    histos.print();
  }

  template <std::size_t N>
  static uint64_t featureHash(std::array<float, N> const& features)
  {
    uint64_t hash = FnvOffsetBasis;
    for (const auto& value : features) {
      const auto bits = std::bit_cast<uint32_t>(value);
      for (unsigned int shift = 0; shift < BitsPerFloat; shift += BitsPerByte) {
        hash = (hash ^ ((bits >> shift) & ByteMask)) * FnvPrime;
      }
    }
    return hash;
  }

  template <typename Collision>
  void beginEvent(Collision const& collision, DecayMode mode)
  {
    mlCascadeRows.clear();
    mlV0Rows.clear();
    mlTrackRows.clear();
    mlEventRow = -1;
    // ResoCollisions_001 carries no run number or BC; the reduced collision row identifies the event within its DF.
    pendingEvent = PendingEvent{mode, static_cast<int64_t>(collision.globalIndex()), collision.posZ(), collision.bMagField(),
                                collision.cent(), collision.multiplicity(), collision.isRecINELgt0(), false, -1};
    if constexpr (requires { collision.originalMcCollisionId(); }) {
      pendingEvent.isMC = true;
      pendingEvent.mcCollision = static_cast<int64_t>(collision.originalMcCollisionId());
    }
  }

  void ensureEventRow()
  {
    if (mlEventRow >= 0) {
      return;
    }
    mlEvents(static_cast<uint8_t>(pendingEvent.mode), pendingEvent.recoCollision, pendingEvent.posZ, pendingEvent.bField,
             pendingEvent.centrality, pendingEvent.multiplicity, pendingEvent.recINELgt0);
    mlEventRow = mlEvents.lastIndex();
    if (pendingEvent.isMC) {
      mlEventMcs(mlEventRow, pendingEvent.mcCollision);
    }
  }

  template <typename Particle>
  void writeMcParticle(Particle const& particle, uint8_t selection)
  {
    const auto mothers = particle.mothersIds();
    std::array<int, NStoredMothers> motherRows{-1, -1};
    for (std::size_t i = 0; i < NStoredMothers && i < mothers.size(); ++i) {
      motherRows[i] = mothers[i];
    }
    const auto daughters = particle.daughtersIds();
    std::array<int, 2> daughterRows{daughters[0], daughters[1]};
    const auto nMothers = static_cast<uint8_t>(std::min<std::size_t>(mothers.size(), std::numeric_limits<uint8_t>::max()));
    mlMcParticles(static_cast<int64_t>(particle.globalIndex()), selection, particle.pdgCode(), static_cast<int64_t>(particle.mcCollisionId()),
                  nMothers, motherRows.data(), daughterRows.data(), particle.px(), particle.py(), particle.pz(), particle.e(), particle.weight());
  }

  // McParticles row of a source track, -1 without MC tables or label
  int mcLabel(int64_t trackId) const
  {
    if (mcTrackLabelsDf == nullptr) {
      return -1;
    }
    const auto row = trackId - static_cast<int64_t>(mcTrackLabelsDf->offset());
    if (row < 0 || row >= mcTrackLabelsDf->size()) {
      return -1;
    }
    return mcTrackLabelsDf->rawIteratorAt(row).mcParticleId();
  }

  // Writes the label particle and its first-mother chain, once per DF
  void writeAncestry(int label)
  {
    if (mcParticlesDf == nullptr) {
      return;
    }
    int64_t row = label;
    for (int step = 0; step <= MaxAncestorSteps; ++step) {
      const auto local = row - static_cast<int64_t>(mcParticlesDf->offset());
      if (row < 0 || local < 0 || local >= mcParticlesDf->size()) {
        return;
      }
      const auto particle = mcParticlesDf->rawIteratorAt(local);
      if (mlAncestryRows.insert(row).second) {
        writeMcParticle(particle, McSelectionAncestry);
      }
      const auto mothers = particle.mothersIds();
      if (mothers.empty()) {
        return;
      }
      row = mothers[0];
    }
  }

  template <typename Cascade>
  int64_t writeCascade(Cascade const& xi)
  {
    const auto id = static_cast<int64_t>(xi.globalIndex());
    if (auto it = mlCascadeRows.find(id); it != mlCascadeRows.end()) {
      return it->second;
    }
    const auto indices = xi.cascadeIndices();
    const std::array<int, 3> daughterIds{indices[0], indices[1], indices[2]};
    mlCascades(mlEventRow, id, xi.px(), xi.py(), xi.pz(), static_cast<int8_t>(xi.sign()), xi.mXi(), xi.mLambda(),
               xi.v0CosPA(), xi.cascCosPA(), xi.daughDCA(), xi.cascDaughDCA(),
               xi.dcapostopv(), xi.dcanegtopv(), xi.dcabachtopv(), xi.dcav0topv(), xi.dcaXYCascToPV(), xi.dcaZCascToPV(),
               xi.transRadius(), xi.cascTransRadius(), xi.decayVtxX(), xi.decayVtxY(), xi.decayVtxZ(),
               xi.daughterTPCNSigmaPosPi10(), xi.daughterTPCNSigmaPosKa10(), xi.daughterTPCNSigmaPosPr10(),
               xi.daughterTPCNSigmaNegPi10(), xi.daughterTPCNSigmaNegKa10(), xi.daughterTPCNSigmaNegPr10(),
               xi.daughterTPCNSigmaBachPi10(), xi.daughterTPCNSigmaBachKa10(), xi.daughterTPCNSigmaBachPr10(),
               xi.daughterTOFNSigmaPosPi10(), xi.daughterTOFNSigmaPosKa10(), xi.daughterTOFNSigmaPosPr10(),
               xi.daughterTOFNSigmaNegPi10(), xi.daughterTOFNSigmaNegKa10(), xi.daughterTOFNSigmaNegPr10(),
               xi.daughterTOFNSigmaBachPi10(), xi.daughterTOFNSigmaBachKa10(), xi.daughterTOFNSigmaBachPr10(),
               xi.nCrossedRowsPos(), xi.nCrossedRowsNeg(), xi.nCrossedRowsBach(), daughterIds.data());
    const auto row = mlCascades.lastIndex();
    mlCascadeRows.emplace(id, row);
    if (mcTrackLabelsDf != nullptr) {
      std::array<int, 3> labels{};
      for (std::size_t i = 0; i < labels.size(); ++i) {
        labels[i] = mcLabel(daughterIds[i]);
        writeAncestry(labels[i]);
      }
      mlCascadeMcs(row, labels.data());
    }
    return row;
  }

  template <typename V0>
  int64_t writeV0(V0 const& v0)
  {
    const auto id = static_cast<int64_t>(v0.globalIndex());
    if (auto it = mlV0Rows.find(id); it != mlV0Rows.end()) {
      return it->second;
    }
    const auto indices = v0.indices();
    const std::array<int, 2> daughterIds{indices[0], indices[1]};
    mlV0s(mlEventRow, id, v0.px(), v0.py(), v0.pz(), v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(),
          v0.v0CosPA(), v0.daughDCA(), v0.dcapostopv(), v0.dcanegtopv(), v0.dcav0topv(), v0.transRadius(),
          v0.decayVtxX(), v0.decayVtxY(), v0.decayVtxZ(), v0.alpha(), v0.qtarm(),
          v0.daughterTPCNSigmaPosPi10(), v0.daughterTPCNSigmaPosKa10(), v0.daughterTPCNSigmaPosPr10(),
          v0.daughterTPCNSigmaNegPi10(), v0.daughterTPCNSigmaNegKa10(), v0.daughterTPCNSigmaNegPr10(),
          v0.daughterTOFNSigmaPosPi10(), v0.daughterTOFNSigmaPosKa10(), v0.daughterTOFNSigmaPosPr10(),
          v0.daughterTOFNSigmaNegPi10(), v0.daughterTOFNSigmaNegKa10(), v0.daughterTOFNSigmaNegPr10(),
          v0.nCrossedRowsPos(), v0.nCrossedRowsNeg(), daughterIds.data());
    const auto row = mlV0s.lastIndex();
    mlV0Rows.emplace(id, row);
    if (mcTrackLabelsDf != nullptr) {
      std::array<int, 2> labels{};
      for (std::size_t i = 0; i < labels.size(); ++i) {
        labels[i] = mcLabel(daughterIds[i]);
        writeAncestry(labels[i]);
      }
      mlV0Mcs(row, labels.data());
    }
    return row;
  }

  template <typename Track>
  int64_t writeTrack(Track const& track)
  {
    const auto id = static_cast<int64_t>(track.globalIndex());
    if (auto it = mlTrackRows.find(id); it != mlTrackRows.end()) {
      return it->second;
    }
    mlTracks(mlEventRow, static_cast<int64_t>(track.trackId()), track.px(), track.py(), track.pz(),
             track.pidNSigmaPiFlag(), track.pidNSigmaKaFlag(), track.pidNSigmaPrFlag(),
             track.trackSelectionFlags(), track.trackFlags(), track.tpcNClsCrossedRows(), track.itsClusterMap());
    const auto row = mlTracks.lastIndex();
    mlTrackRows.emplace(id, row);
    if (mcTrackLabelsDf != nullptr) {
      const auto label = mcLabel(static_cast<int64_t>(track.trackId()));
      writeAncestry(label);
      mlTrackMcs(row, label);
    }
    return row;
  }

  template <bool IsMC, typename Collision, typename Cascade, typename V0>
  void writeXiK0sCandidate(Collision const& collision, Cascade const& xi, V0 const& v0, XiK0sCandidateValues const& values, uint16_t passBits)
  {
    const auto canonical = ml::canonicalizeXiK0s(ml::makeCascadeSnapshot(collision, xi), ml::makeV0Snapshot(collision, v0));
    if (canonical.status != ml::BuildStatus::Ok) {
      histos.fill(HIST("ML/xiK0s/exportSkipped"), static_cast<int>(canonical.status));
      return;
    }
    const auto pack = ml::buildXiK0sFeatures(canonical.candidate);
    if (pack.status != ml::BuildStatus::Ok) {
      histos.fill(HIST("ML/xiK0s/exportSkipped"), static_cast<int>(pack.status));
      return;
    }
    if (parityLoggedXiK0s < omega2012MlParityRows) {
      LOGP(info, "OMEGA2012MLPARITY mode=XiK0s collision={} cascade={} v0={} passBits={} featureHash={}",
           collision.globalIndex(), xi.globalIndex(), v0.globalIndex(), passBits, featureHash(pack.master));
      ++parityLoggedXiK0s;
    }
    ensureEventRow();
    const auto cascadeRow = writeCascade(xi);
    const auto v0Row = writeV0(v0);
    mlXiK0sCandidates(mlEventRow, cascadeRow, v0Row,
                      static_cast<float>(values.omega.M()), static_cast<float>(values.omega.Pt()),
                      static_cast<float>(values.omega.Rapidity()), static_cast<float>(values.omega.Eta()),
                      static_cast<float>(values.omega.Phi()), values.alpha, static_cast<int8_t>(xi.sign()), passBits);
    const auto row = mlXiK0sCandidates.lastIndex();
    mlXiK0sInputs(row, pack.master.data(), static_cast<uint8_t>(pack.status));
    if constexpr (IsMC) {
      const bool matched = classifyXiK0sTruth(xi, v0) == XiK0sTruth::Matched;
      mlXiK0sTruth(row, matched ? uint8_t{1} : uint8_t{2}, matched ? xi.motherPDG() : 0,
                   matched ? static_cast<int64_t>(xi.motherId()) : int64_t{-1}, v0.motherPDG(), xi.motherPDG());
    } else {
      mlXiK0sTruth(row, uint8_t{0}, 0, int64_t{-1}, 0, 0);
    }
  }

  template <bool IsMC, typename Collision, typename Cascade, typename Track>
  void writeXi1530KCandidate(Collision const& collision, Cascade const& xi, Track const& pion, Track const& kaon,
                             Xi1530KCandidateValues const& values, uint16_t passBits)
  {
    histos.fill(HIST("ML/xi1530K/chargePattern"), values.chargePattern);
    if (values.chargePattern != ml::kSignalPattern && !omega2012MlExportWrongSign) {
      return; // charge-pattern controls are exported only on request
    }
    const auto canonical = ml::canonicalizeXi1530K(ml::makeCascadeSnapshot(collision, xi), ml::makeTrackSnapshot(pion), ml::makeTrackSnapshot(kaon));
    if (canonical.status != ml::BuildStatus::Ok) {
      histos.fill(HIST("ML/xi1530K/exportSkipped"), static_cast<int>(canonical.status));
      return;
    }
    const auto pack = ml::buildXi1530KFeatures(canonical.candidate);
    if (pack.status != ml::BuildStatus::Ok) {
      histos.fill(HIST("ML/xi1530K/exportSkipped"), static_cast<int>(pack.status));
      return;
    }
    if (parityLoggedXi1530K < omega2012MlParityRows) {
      LOGP(info, "OMEGA2012MLPARITY mode=Xi1530K collision={} cascade={} tracks={},{} chargePattern={} passBits={} featureHash={}",
           collision.globalIndex(), xi.globalIndex(), canonical.candidate.pion.sourceTrackId, canonical.candidate.kaon.sourceTrackId,
           values.chargePattern, passBits, featureHash(pack.master));
      ++parityLoggedXi1530K;
    }
    ensureEventRow();
    const auto cascadeRow = writeCascade(xi);
    const auto pionRow = writeTrack(pion);
    const auto kaonRow = writeTrack(kaon);
    mlXi1530KCandidates(mlEventRow, cascadeRow, pionRow, kaonRow,
                        static_cast<float>(values.omega.M()), values.massXiPi, values.massXiK, values.massPiK,
                        static_cast<float>(values.omega.Pt()), static_cast<float>(values.omega.Rapidity()),
                        static_cast<float>(values.omega.Eta()), static_cast<float>(values.omega.Phi()),
                        values.chargePattern, passBits);
    const auto row = mlXi1530KCandidates.lastIndex();
    mlXi1530KInputs(row, pack.master.data(), static_cast<uint8_t>(pack.status));
    if constexpr (IsMC) {
      const bool matched = classifyXi1530KTruth(xi, pion, kaon) == Xi1530KTruth::Matched;
      mlXi1530KTruth(row, matched ? uint8_t{1} : uint8_t{2}, matched ? kaon.motherPDG() : 0,
                     matched ? static_cast<int64_t>(kaon.motherId()) : int64_t{-1}, static_cast<int64_t>(xi.motherId()));
    } else {
      mlXi1530KTruth(row, uint8_t{0}, 0, int64_t{-1}, int64_t{-1});
    }
  }

  template <bool IsMC, typename Collision, typename Cascades, typename V0s>
  void exportXiK0s(Collision const& collision, Cascades const& cascades, V0s const& v0s)
  {
    beginEvent(collision, DecayMode::XiK0s);
    // Selection only: the analysis histograms belong to the Omega(2012) analysis task
    core.forEachXiK0sCandidate<IsMC, false>(histos, collision, collision, cascades, v0s, true, nullptr, nullptr, nullptr,
                                            [&](auto const& coll, auto const& xi, auto const& v0, XiK0sCandidateValues const& values, uint16_t passBits) {
                                              writeXiK0sCandidate<IsMC>(coll, xi, v0, values, passBits);
                                            });
  }

  template <bool IsMC, typename Collision, typename Cascades, typename Tracks>
  void exportXi1530K(Collision const& collision, Cascades const& cascades, Tracks const& tracks)
  {
    beginEvent(collision, DecayMode::Xi1530K);
    core.forEachXi1530KCandidate<IsMC, false, true>(histos, collision, cascades, tracks, nullptr, true, nullptr, nullptr, nullptr,
                                                    [&](auto const& coll, auto const& xi, auto const& pion, auto const& kaon,
                                                        Xi1530KCandidateValues const& values, uint16_t passBits) {
                                                      writeXi1530KCandidate<IsMC>(coll, xi, pion, kaon, values, passBits);
                                                    });
  }

  void processXiK0s(ResoCollisions::iterator const& collision, aod::ResoCascades const& cascades, aod::ResoV0s const& v0s)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    exportXiK0s<false>(collision, cascades, v0s);
  }
  PROCESS_SWITCH(Omega2012TrainingTable, processXiK0s, "Write XiK0s-mode candidates from data", true);

  void processXiK0sMC(ResoMCCols::iterator const& collision, ResoMCCascades const& cascades, ResoMCV0s const& v0s,
                      aod::McParticles const& mcParticles, aod::McTrackLabels const& mcTrackLabels)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    mcParticlesDf = &mcParticles;
    mcTrackLabelsDf = &mcTrackLabels;
    exportXiK0s<true>(collision, cascades, v0s);
    // cppcheck-suppress redundantAssignment
    mcParticlesDf = nullptr;
    // cppcheck-suppress redundantAssignment
    mcTrackLabelsDf = nullptr;
  }
  PROCESS_SWITCH(Omega2012TrainingTable, processXiK0sMC, "Write XiK0s-mode candidates with truth from reconstructed MC", false);

  void processXi1530KMicro(ResoCollisions::iterator const& collision, aod::ResoCascades const& cascades, ResoMicroTracks const& tracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    exportXi1530K<false>(collision, cascades, tracks);
  }
  PROCESS_SWITCH(Omega2012TrainingTable, processXi1530KMicro, "Write Xi1530K-mode candidates from data micro v001 tables", false);

  void processXi1530KMCMicro(ResoMCCols::iterator const& collision, ResoMCCascades const& cascades, ResoMCMicroTracks const& tracks,
                             aod::McParticles const& mcParticles, aod::McTrackLabels const& mcTrackLabels)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    mcParticlesDf = &mcParticles;
    mcTrackLabelsDf = &mcTrackLabels;
    exportXi1530K<true>(collision, cascades, tracks);
    // cppcheck-suppress redundantAssignment
    mcParticlesDf = nullptr;
    // cppcheck-suppress redundantAssignment
    mcTrackLabelsDf = nullptr;
  }
  PROCESS_SWITCH(Omega2012TrainingTable, processXi1530KMCMicro, "Write Xi1530K-mode candidates with truth from reconstructed MC micro v001 tables", false);

  void processMCTrue(ResoMCCols::iterator const& collision, ResoMCParents const& resoParents)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    core.forEachGeneratedOmega2012(histos, resoParents, [&](auto const& part, GeneratedChannel channel) {
      mlGenAudit(collision.globalIndex(), static_cast<int64_t>(part.originalMcParticleId()),
                 part.pdgCode(), part.daughterPDG1(), part.daughterPDG2(), static_cast<uint8_t>(channel), part.pt(), part.y());
    });
  }
  PROCESS_SWITCH(Omega2012TrainingTable, processMCTrue, "Write generated Omega(2012) parents of selected reconstructed MC events", false);

  // Once per DF: every Omega(2012) with McTreeDepth generations of descendants. The ancestry set is reset here, so
  // the result does not depend on the order of this call relative to the per-collision process functions.
  void processMcTruth(aod::McParticles const& mcParticles)
  {
    mlAncestryRows.clear();
    const auto offset = static_cast<int64_t>(mcParticles.offset());
    const auto nRows = mcParticles.size();
    std::unordered_set<int64_t> treeRows;
    std::vector<std::pair<int64_t, int>> stack;
    for (const auto& parent : mcParticles) {
      if (std::abs(parent.pdgCode()) != o2::constants::physics::Pdg::kOmega2012Minus) {
        continue;
      }
      stack.assign(1, {parent.globalIndex(), 0});
      while (!stack.empty()) {
        const auto [row, depth] = stack.back();
        stack.pop_back();
        const auto local = row - offset;
        if (local < 0 || local >= nRows || !treeRows.insert(row).second) {
          continue;
        }
        const auto particle = mcParticles.rawIteratorAt(local);
        writeMcParticle(particle, McSelectionTree);
        const auto daughters = particle.daughtersIds();
        if (depth >= McTreeDepth || daughters[0] < 0 || daughters[1] < daughters[0]) {
          continue;
        }
        for (int64_t daughter = daughters[0]; daughter <= daughters[1]; ++daughter) {
          stack.emplace_back(daughter, depth + 1);
        }
      }
    }
  }
  PROCESS_SWITCH(Omega2012TrainingTable, processMcTruth, "Write the pruned McParticles rows of the Omega(2012) decay trees; required with the reconstructed-MC process functions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Omega2012TrainingTable>(cfgc)};
}
