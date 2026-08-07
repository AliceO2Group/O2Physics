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

/// \file treeCreatorPidTpcQa.cxx
/// \brief Creates trees with PID QA variables along with variables used for NN training
///
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>

#include "treeCreatorPidTpcQa.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DPG/Tasks/TPC/tpcSkimsTableCreator.h"
#include "DPG/Tasks/TPC/utilsTpcSkimsTableCreator.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TString.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::track;
using namespace o2::dpg_tpcskimstablecreator;

#define DO_FOR_ALL_PARTICLES(MACRO) \
  MACRO(El, Electron)               \
  MACRO(Mu, Muon)                   \
  MACRO(Pi, Pion)                   \
  MACRO(Ka, Kaon)                   \
  MACRO(Pr, Proton)                 \
  MACRO(De, Deuteron)               \
  MACRO(Tr, Triton)                 \
  MACRO(He, Helium3)                \
  MACRO(Al, Alpha)

struct TreeCreatorPidTpcQa {
  Produces<o2::aod::QaPidTpc> rowPidTpcQa;

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<float> cutVtxZ{"cutVtxZ", 10.f, "Cut on vertex Z position [cm]"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> requireGlobalTrack{"requireGlobalTrack", true, "Skip non-global tracks"};
  Configurable<bool> requireIts{"requireIts", true, "Skip tracks without ITS"};
  Configurable<int16_t> cutMinTPCNcls{"cutMinTPCNcls", 0, "Minimum number or TPC Clusters for tracks"};
  Configurable<float> cutRapidity{"cutRapidity", 999.f, "Rapidity cut"};
  Configurable<float> nClNorm{"nClNorm", 152.f, "Number of cluster normalization. Run 2: 159, Run 3 152"};
  // Configurable for the path of CCDB General Run Parameters LHC Interface information
  Configurable<std::string> ccdbPathGrpLhcIf{"ccdbPathGrpLhcIf", "GLO/Config/GRPLHCIF", "Path on the CCDB for the GRPLHCIF object"};
  // Configurables for output tables reservation size
  Configurable<float> reserveRatio{"reserveRatio", 1.f, "Ratio of how many rows expected in the output table to the input Tracks table size"};
  Configurable<bool> saveReserveQaHisto{"saveReserveQaHisto", true, "Flag to save the DF-wise ratio of output table size to that of input table"};
  // Configurables for run condtion table
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> checkZdc{"checkZdc", false, "set ZDC flag for PbPb"};
  Configurable<bool> treatLimitedAcceptanceAsBad{"treatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  Configurable<bool> requireGoodRct{"requireGoodRct", false, "require good detector flag in run condtion table"};

#define DECLARE_PARTICLE_WISE_CONFIGURABLES(ParticleNameShort, ParticleNameLong)                                                                                                                                                                          \
  Configurable<float> cutTpcInnerParameterMin##ParticleNameLong{"cutTpcInnerParameterMin" #ParticleNameLong, 0.f, "Lower-value cut on tpcInnerParam for " #ParticleNameLong};   /* o2-linter: disable=name/configurable (Configurable defined in macro)*/ \
  Configurable<float> cutTpcInnerParameterMax##ParticleNameLong{"cutTpcInnerParameterMax" #ParticleNameLong, 999.f, "Upper-value cut on tpcInnerParam for " #ParticleNameLong}; /* o2-linter: disable=name/configurable (Configurable defined in macro)*/ \
  Configurable<float> cutNSigmaTpcAbs##ParticleNameLong{"cutNSigmaTpcAbs" #ParticleNameLong, 999.f, "Cut on absolute value of nSigmaTpc for " #ParticleNameLong};               // o2-linter: disable=name/configurable (Configurable defined in macro)

  DO_FOR_ALL_PARTICLES(DECLARE_PARTICLE_WISE_CONFIGURABLES)
#undef DECLARE_PARTICLE_WISE_CONFIGURABLES

#define PACK_CONFIGURABLES_TO_ARRAY(ParticleNameShort, ParticleNameLong) &cutTpcInnerParameterMin##ParticleNameLong,
  std::array<Configurable<float>*, PID::Alpha + 1> cutTpcInnerParameterMin{
    DO_FOR_ALL_PARTICLES(PACK_CONFIGURABLES_TO_ARRAY)};
#undef PACK_CONFIGURABLES_TO_ARRAY

#define PACK_CONFIGURABLES_TO_ARRAY(ParticleNameShort, ParticleNameLong) &cutTpcInnerParameterMax##ParticleNameLong,
  std::array<Configurable<float>*, PID::Alpha + 1> cutTpcInnerParameterMax{
    DO_FOR_ALL_PARTICLES(PACK_CONFIGURABLES_TO_ARRAY)};
#undef PACK_CONFIGURABLES_TO_ARRAY

#define PACK_CONFIGURABLES_TO_ARRAY(ParticleNameShort, ParticleNameLong) &cutNSigmaTpcAbs##ParticleNameLong,
  std::array<Configurable<float>*, PID::Alpha + 1> cutNSigmaTpcAbs{
    DO_FOR_ALL_PARTICLES(PACK_CONFIGURABLES_TO_ARRAY)};
#undef PACK_CONFIGURABLES_TO_ARRAY

  HistogramRegistry registry{"registry", {}};

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  ctpRateFetcher mRateFetcher;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  using CollisionsExtra = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;

  Preslice<TrackCandidates> perCollisionTracks = aod::track::collisionId;

  int mEnabledParticles{0};
  int mProcessedParticles{0};

  template <o2::track::PID::ID ParticleId>
  bool initPerParticle()
  {
    static_assert(ParticleId >= 0 && ParticleId <= PID::Alpha && "Particle index outside limits");
    int enabledProcesses{0};

    switch (ParticleId) {
#define INIT_PARTICLE(ParticleNameShort, ParticleNameLong)                                                             \
  case PID::ParticleNameLong:                                                                                          \
    if (!doprocess##ParticleNameLong && !doprocessFull##ParticleNameLong && !doprocessFullWithTOF##ParticleNameLong) { \
      return false;                                                                                                    \
    }                                                                                                                  \
    if (doprocess##ParticleNameLong) {                                                                                 \
      ++enabledProcesses;                                                                                              \
    }                                                                                                                  \
    if (doprocessFull##ParticleNameLong) {                                                                             \
      ++enabledProcesses;                                                                                              \
    }                                                                                                                  \
    if (doprocessFullWithTOF##ParticleNameLong) {                                                                      \
      ++enabledProcesses;                                                                                              \
    }                                                                                                                  \
    LOG(info) << "Enabled TPC QA for " << #ParticleNameLong;                                                           \
    break;

      DO_FOR_ALL_PARTICLES(INIT_PARTICLE)
#undef INIT_PARTICLE
    }
    if (enabledProcesses != 1) {
      LOG(fatal) << "Cannot enable more than one process function per particle, check and retry!";
    }
    return true;
  }

  void init(o2::framework::InitContext&)
  {
    static_for<0, PID::Alpha>([&](auto ParticleId) {
      mEnabledParticles += static_cast<int>(initPerParticle<ParticleId>());
    });

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    rctChecker.init(rctLabel, checkZdc, treatLimitedAcceptanceAsBad);

    if (saveReserveQaHisto) {
      registry.add("hOutputRatio", "Table out/in ratio;Table out/in ratio;Entries", {HistType::kTH1F, {{100, 0, reserveRatio}}});
    }
  }

  template <o2::track::PID::ID ParticleId, bool IsFullTable, bool IsTofTable, typename TrackType>
  void processSingleParticle(CollisionsExtra const& collisions,
                             TrackType const& tracks)
  {
    if (mProcessedParticles == 0) {
      rowPidTpcQa.reserve(tracks.size() * reserveRatio);
    }

    std::string irSource{};
    float sqrtSNN{}; // placeholder to satisfy evaluateIrSourceAndSqrtSnn's signature
    bool isFirstCollision{true};
    for (const auto& collision : collisions) {
      if (!isEventSelected(collision, applyEvSel) || (std::abs(collision.posZ()) > cutVtxZ)) {
        continue;
      }

      const bool isGoodRctEvent = rctChecker.checkTable(collision);
      if (requireGoodRct && !isGoodRctEvent) {
        continue;
      }

      const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (isFirstCollision) {
        evaluateIrSourceAndSqrtSnn(ccdb, ccdbPathGrpLhcIf, bc.timestamp(), irSource, sqrtSNN);
      }
      isFirstCollision = false;
      const float ft0Occ = collision.ft0cOccupancyInTimeRange();
      const float multTPC = collision.multTPC() / MultiplicityNorm;
      const auto hadronicRate = !irSource.empty() ? mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * OneToKilo : 0.;

      const auto tracksFromCollision = tracks.sliceBy(perCollisionTracks, static_cast<int>(collision.globalIndex()));

      for (const auto& track : tracksFromCollision) {
        bool isGoodTrack = isTrackSelected(track, trackSelection);
        isGoodTrack &= (!requireGlobalTrack || track.isGlobalTrack());
        isGoodTrack &= (!requireIts || track.hasITS());
        isGoodTrack &= track.hasTPC();
        isGoodTrack &= (track.tpcNClsFound() >= cutMinTPCNcls);

        const float rapidity = track.rapidity(PID::getMass(ParticleId));
        const float momentum = track.p();
        const float nClNormalized = std::sqrt(nClNorm / track.tpcNClsFound());
        const auto nclPID = static_cast<float>(track.tpcNClsPID());
        const float phi = track.phi();
        const float tgl = track.tgl();
        const float tpcInnerParam = track.tpcInnerParam();
        const float signed1Pt = track.signed1Pt();
        const float nSigmaTpc = o2::aod::pidutils::tpcNSigma<ParticleId>(track);

        isGoodTrack &= (std::fabs(rapidity) <= cutRapidity);
        isGoodTrack &= (tpcInnerParam >= *cutTpcInnerParameterMin.at(ParticleId));
        isGoodTrack &= (tpcInnerParam <= *cutTpcInnerParameterMax.at(ParticleId));
        isGoodTrack &= (std::fabs(nSigmaTpc) <= *cutNSigmaTpcAbs.at(ParticleId));

        if (!isGoodTrack) {
          continue;
        }

        float dedxDiff{UndefValueFloat};
        float dedxExpected{UndefValueFloat};
        float expSigma{UndefValueFloat};

        if constexpr (IsFullTable) {
          dedxDiff = o2::aod::pidutils::tpcExpSignalDiff<ParticleId>(track);
          dedxExpected = track.tpcSignal() - dedxDiff;
          expSigma = o2::aod::pidutils::tpcExpSigma<ParticleId>(track);
        }

        float nSigmaTof{UndefValueFloat};

        if constexpr (IsTofTable) {
          nSigmaTof = o2::aod::pidutils::tofNSigma<ParticleId>(track);
        }

        rowPidTpcQa(isGoodRctEvent, ParticleId, ft0Occ, hadronicRate, multTPC, nClNormalized, nclPID, phi, tgl, tpcInnerParam, rapidity, momentum, signed1Pt, nSigmaTpc, dedxExpected, dedxDiff, expSigma, nSigmaTof);
      } // tracksFromCollision
    } // collisions
    ++mProcessedParticles;
    if (mProcessedParticles == mEnabledParticles) {
      mProcessedParticles = 0;
      if (saveReserveQaHisto) {
        registry.fill(HIST("hOutputRatio"), static_cast<double>((rowPidTpcQa.lastIndex() + 1)) / tracks.size());
      }
    }
  }

#define MAKE_PROCESS_FUNCTIONS(ParticleNameShort, ParticleNameLong)                                                                                           \
  void process##ParticleNameLong(CollisionsExtra const& collisions,                                                                                           \
                                 soa::Join<TrackCandidates, aod::pidTPC##ParticleNameShort> const& tracks,                                                    \
                                 aod::BCsWithTimestamps const&)                                                                                               \
  {                                                                                                                                                           \
    processSingleParticle<PID::ParticleNameLong, false, false>(collisions, tracks);                                                                           \
  }                                                                                                                                                           \
  PROCESS_SWITCH(TreeCreatorPidTpcQa, process##ParticleNameLong, Form("Process for the %s hypothesis for TPC NSigma QA", #ParticleNameLong), false);          \
                                                                                                                                                              \
  void processFull##ParticleNameLong(CollisionsExtra const& collisions,                                                                                       \
                                     soa::Join<TrackCandidates, aod::pidTPCFull##ParticleNameShort> const& tracks,                                            \
                                     aod::BCsWithTimestamps const&)                                                                                           \
  {                                                                                                                                                           \
    processSingleParticle<PID::ParticleNameLong, true, false>(collisions, tracks);                                                                            \
  }                                                                                                                                                           \
  PROCESS_SWITCH(TreeCreatorPidTpcQa, processFull##ParticleNameLong, Form("Process for the %s hypothesis for full TPC PID QA", #ParticleNameLong), false);    \
                                                                                                                                                              \
  void processFullWithTOF##ParticleNameLong(CollisionsExtra const& collisions,                                                                                \
                                            soa::Join<TrackCandidates, aod::pidTPCFull##ParticleNameShort, aod::pidTOFFull##ParticleNameShort> const& tracks, \
                                            aod::BCsWithTimestamps const&)                                                                                    \
  {                                                                                                                                                           \
    processSingleParticle<PID::ParticleNameLong, true, true>(collisions, tracks);                                                                             \
  }                                                                                                                                                           \
  PROCESS_SWITCH(TreeCreatorPidTpcQa, processFullWithTOF##ParticleNameLong, Form("Process for the %s hypothesis for full TPC PID QA with the TOF info added", #ParticleNameLong), false);

  DO_FOR_ALL_PARTICLES(MAKE_PROCESS_FUNCTIONS)
#undef MAKE_PROCESS_FUNCTIONS
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TreeCreatorPidTpcQa>(cfgc)};
}
#undef DO_FOR_ALL_PARTICLES
