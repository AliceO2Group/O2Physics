// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License version 3, copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file k1TrainingTable.cxx
/// \brief Derived K1 training table workflow
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#include "PWGLF/Core/K1AnalysisMicroCore.h"
#include "PWGLF/Core/K1MlFeatures.h"
#include "PWGLF/Core/ResoAnalysisSelectionCore.h"
#include "PWGLF/DataModel/LFK1MlTables.h"
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

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <TH1.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::physics;
using namespace o2::analysis::resonance;
using namespace o2::analysis::k1micro;

static_assert(std::extent_v<aod::k1ml::K1MasterFeatures::type> == o2::analysis::k1ml::NMasterFeatures,
              "K1MasterFeatures column size must match the K1 ML feature contract");

/// Writes the unlike-sign K1 micro candidates of the shared K1 selection, with the
/// canonical (kaon, same-sign pion, opposite-sign pion) tracks and the master features.
struct K1TrainingTable {
  using ResoCollisions = aod::ResoCollisions_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCMicroTracks = soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001>;
  using ResoMCParents = aod::ResoMCParents_001;

  // FNV-1a over the bit patterns of the master features, for parity logs
  static constexpr uint64_t FnvOffsetBasis = 14695981039346656037ULL;
  static constexpr uint64_t FnvPrime = 1099511628211ULL;
  static constexpr unsigned int BitsPerByte = 8;
  static constexpr unsigned int BitsPerFloat = 32;
  static constexpr uint32_t ByteMask = 0xffU;

  Produces<aod::K1MlEvents> k1MlEvents;
  Produces<aod::K1MlTracks> k1MlTracks;
  Produces<aod::K1MlCandidates> k1MlCandidates;
  Produces<aod::K1MlInputs> k1MlInputs;
  Produces<aod::K1MlTruth> k1MlTruth;
  Produces<aod::K1MlGenAudit> k1MlGenAudit;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Selection shared with the K1 histogram task (plain JSON keys, no group prefix)
  EventCuts eventCuts;
  TrackCuts trackCuts;
  PionPidCuts pionPID;
  KaonPidCuts kaonPID;
  SecondaryCuts secondaryCuts;
  CandidateCuts candidateCuts;

  Configurable<std::string> k1MlExportStage{"k1MlExportStage", "loose", "Candidate export stage: loose (pass bits >= 1) or selected (pass bits = 31)"};
  Configurable<bool> k1MlLooseAudit{"k1MlLooseAudit", true, "Record loose US cutflow and mass/pT/activity spectrum"};
  Configurable<int> k1MlParityRows{"k1MlParityRows", 0, "Log feature hashes of the first N exported candidates"};

  K1AnalysisMicroCore core;
  int64_t k1MlEventRow = -1;
  std::unordered_map<int64_t, int64_t> k1MlTrackRows;
  int k1MlParityLogged = 0;

  void init(InitContext&)
  {
    if (static_cast<int>(doprocessResoMicroTracks) + static_cast<int>(doprocessMCMicro) != 1 ||
        (doprocessMCTrue && !doprocessMCMicro)) {
      LOG(fatal) << "K1 training table requires exactly one of processResoMicroTracks and processMCMicro; processMCTrue requires processMCMicro";
    }
    if (k1MlParityRows < 0) {
      LOG(fatal) << "k1MlParityRows must not be negative";
    }
    if (k1MlExportStage.value != "loose" && k1MlExportStage.value != "selected") {
      LOG(fatal) << "k1MlExportStage must be loose or selected";
    }

    ProcessModes modes;
    modes.microTracks = true;
    modes.mcReco = doprocessMCMicro;
    modes.mcRecoMicro = doprocessMCMicro;
    modes.mcGen = doprocessMCTrue;
    LooseStageOptions looseOptions;
    looseOptions.audit = k1MlLooseAudit;
    looseOptions.exportSelected = k1MlExportStage.value == "selected";
    core.init(histos, eventCuts, trackCuts, pionPID, kaonPID, secondaryCuts, candidateCuts, modes, looseOptions);
    if (doprocessMCMicro) {
      histos.add("MCReco/collisions", "Selected reconstructed MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/microTracks", "Input micro tracks in selected MC collisions", HistType::kTH1D, {{1, 0, 1}});
    }

    // Candidates that violate the canonical or feature contract are skipped, never written.
    auto skipped = histos.add<TH1>("ML/exportSkipped", "Skipped candidates;K1 ML build status;candidates", HistType::kTH1D, {{6, -0.5, 5.5}});
    const std::array<const char*, 6> statusLabels{"Ok", "InvalidChargePattern", "ReusedTrack", "InvalidMomentum", "InvalidKinematics", "InvalidContract"};
    for (std::size_t i = 0; i < statusLabels.size(); ++i) {
      skipped->GetXaxis()->SetBinLabel(i + 1, statusLabels[i]);
    }

    LOG(info) << "Size of the histograms in K1 training table task";
    histos.print();
  }

  template <typename Track>
  int64_t writeK1MlTrack(Track const& track)
  {
    const auto id = static_cast<int64_t>(track.globalIndex());
    if (auto it = k1MlTrackRows.find(id); it != k1MlTrackRows.end()) {
      return it->second;
    }
    k1MlTracks(k1MlEventRow, static_cast<int64_t>(track.trackId()), track.px(), track.py(), track.pz(),
               track.pidNSigmaPiFlag(), track.pidNSigmaKaFlag(), track.pidNSigmaPrFlag(),
               track.trackSelectionFlags(), track.trackFlags(), track.tpcNClsCrossedRows(), track.itsClusterMap());
    const auto row = static_cast<int64_t>(k1MlTracks.lastIndex());
    k1MlTrackRows.emplace(id, row);
    return row;
  }

  template <typename Collision>
  void writeK1MlEvent(Collision const& collision)
  {
    k1MlTrackRows.clear();
    // ResoCollisions_001 carries no run number or BC; the reduced collision row identifies the event within its DF.
    k1MlEvents(static_cast<int64_t>(collision.globalIndex()),
               collision.posZ(), collision.bMagField(), collision.cent(), collision.multiplicity(), collision.isRecINELgt0());
    k1MlEventRow = static_cast<int64_t>(k1MlEvents.lastIndex());
  }

  template <typename Collision>
  void logParity(Collision const& collision, o2::analysis::k1ml::CandidateSnapshot const& candidate, o2::analysis::k1ml::FeaturePack const& pack)
  {
    if (k1MlParityLogged >= k1MlParityRows) {
      return;
    }
    uint64_t hash = FnvOffsetBasis;
    for (const auto& value : pack.master) {
      const uint32_t bits = std::bit_cast<uint32_t>(value);
      for (unsigned int shift = 0; shift < BitsPerFloat; shift += BitsPerByte) {
        hash = (hash ^ ((bits >> shift) & ByteMask)) * FnvPrime;
      }
    }
    LOGP(info, "K1MLPARITY collision={} tracks={},{},{} featureHash={}",
         collision.globalIndex(), candidate.tracks[0].sourceTrackId, candidate.tracks[1].sourceTrackId,
         candidate.tracks[2].sourceTrackId, hash);
    ++k1MlParityLogged;
  }

  template <bool IsMC, typename Collision, typename Kaon, typename Pion>
  void writeK1MlCandidate(Collision const& collision, Kaon const& kaon, Pion const& samePion, Pion const& oppPion,
                          K1TruthChannel channel, uint16_t passBits)
  {
    using o2::analysis::k1ml::BuildStatus;
    const auto canonical = o2::analysis::k1ml::canonicalizeUS(o2::analysis::k1ml::makeTrackSnapshot(kaon), o2::analysis::k1ml::makeTrackSnapshot(samePion), o2::analysis::k1ml::makeTrackSnapshot(oppPion));
    if (canonical.status != BuildStatus::Ok) {
      histos.fill(HIST("ML/exportSkipped"), static_cast<int>(canonical.status));
      return;
    }
    const auto pack = o2::analysis::k1ml::buildMasterFeatures(canonical.candidate);
    if (pack.status != BuildStatus::Ok) {
      histos.fill(HIST("ML/exportSkipped"), static_cast<int>(pack.status));
      return;
    }
    logParity(collision, canonical.candidate, pack);
    const auto kaonRow = writeK1MlTrack(kaon);
    const auto sameRow = writeK1MlTrack(samePion);
    const auto oppRow = writeK1MlTrack(oppPion);
    ROOT::Math::PxPyPzMVector k{kaon.px(), kaon.py(), kaon.pz(), MassKaonCharged};
    ROOT::Math::PxPyPzMVector s{samePion.px(), samePion.py(), samePion.pz(), MassPionCharged};
    ROOT::Math::PxPyPzMVector o{oppPion.px(), oppPion.py(), oppPion.pz(), MassPionCharged};
    const auto mother = k + s + o;
    k1MlCandidates(k1MlEventRow, kaonRow, sameRow, oppRow,
                   static_cast<float>(mother.M()), static_cast<float>((s + o).M()),
                   static_cast<float>((k + s).M()), static_cast<float>((k + o).M()),
                   pack.kinematics.scalarSumPt, static_cast<float>((s + o).Pt()),
                   static_cast<float>(mother.Pt()), static_cast<float>(mother.Rapidity()),
                   static_cast<float>(mother.Eta()), static_cast<float>(mother.Phi()),
                   static_cast<int8_t>(kaon.sign()), passBits);
    const auto row = static_cast<int64_t>(k1MlCandidates.lastIndex());
    k1MlInputs(row, pack.master.data(), static_cast<uint8_t>(pack.status));
    if constexpr (IsMC) {
      const bool matched = channel != K1TruthChannel::None;
      int motherId = -1;
      if (matched) {
        motherId = channel == K1TruthChannel::RhoK ? kaon.motherId() : std::abs(samePion.motherPDG()) == Pdg::kK1_1270Plus ? samePion.motherId()
                                                                                                                           : oppPion.motherId();
      }
      // Immediate-mother information alone does not prove a physical UID.
      k1MlTruth(row, matched ? uint8_t{1} : uint8_t{2}, static_cast<uint8_t>(channel),
                matched ? kaon.sign() * Pdg::kK1_1270Plus : 0, static_cast<int64_t>(motherId),
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
    } else {
      k1MlTruth(row, uint8_t{0}, uint8_t{0}, 0, int64_t{-1},
                std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
    }
  }

  void processResoMicroTracks(ResoCollisions::iterator const& collision, ResoMicroTracks const& tracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    writeK1MlEvent(collision);
    // Selection only: the K1 analysis histograms belong to the K1 analysis task
    core.forEachCandidate<false, false, true>(histos, collision, tracks, tracks, false, nullptr, nullptr,
                                              [this](auto const& coll, auto const& kaon, auto const& samePion, auto const& oppPion,
                                                     K1TruthChannel channel, uint16_t passBits) {
                                                writeK1MlCandidate<false>(coll, kaon, samePion, oppPion, channel, passBits);
                                              });
  }
  PROCESS_SWITCH(K1TrainingTable, processResoMicroTracks, "Write K1 candidates from data micro v001 tables", true);

  void processMCMicro(ResoMCCols::iterator const& collision, ResoMCMicroTracks const& tracks)
  {
    // The modular producer already selected these reconstructed collisions.
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    histos.fill(HIST("MCReco/microTracks"), 0.5, tracks.size());
    writeK1MlEvent(collision);
    core.forEachCandidate<true, false, true>(histos, collision, tracks, tracks, false, nullptr, nullptr,
                                             [this](auto const& coll, auto const& kaon, auto const& samePion, auto const& oppPion,
                                                    K1TruthChannel channel, uint16_t passBits) {
                                               writeK1MlCandidate<true>(coll, kaon, samePion, oppPion, channel, passBits);
                                             });
  }
  PROCESS_SWITCH(K1TrainingTable, processMCMicro, "Write K1 candidates with truth from reconstructed MC micro v001 tables", false);

  void processMCTrue(ResoMCCols::iterator const& collision, ResoMCParents const& resoParents)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    core.forEachGeneratedK1(histos, resoParents, [&](auto const& part, K1TruthChannel channel) {
      k1MlGenAudit(static_cast<int64_t>(collision.globalIndex()), static_cast<int64_t>(part.originalMcParticleId()),
                   part.pdgCode(), part.daughterPDG1(), part.daughterPDG2(), static_cast<uint8_t>(channel),
                   part.pt(), part.y(), true);
    });
  }
  PROCESS_SWITCH(K1TrainingTable, processMCTrue, "Write generated K1 parents of selected reconstructed MC events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K1TrainingTable>(cfgc)};
}
