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

/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"
#include "FemtoDreamPairCleaner.h"
#include "FemtoDreamContainer.h"
#include "FemtoDreamDetaDphiStar.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 5.f, 5.f, 100.f}};

static const std::vector<float> kNsigma = {3.5f, 3.f, 2.5f};

} // namespace

struct femtoDreamPairTaskTrackV0 {

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};

  /// Particle 1 (track)
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FemtoDreamParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2 (V0)
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 3122, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 338, "Particle 2 - Selection bit"};

  /// Partition for particle 2
  Partition<aod::FemtoDreamParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) && ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry);
    trackHistoPartTwo.init(&qaRegistry);

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;

    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processSameEvent(o2::aod::FemtoDreamCollision& col,
                        o2::aod::FemtoDreamParticles& parts)
  {
    const auto& tmstamp = col.timestamp();
    const auto& magFieldTesla = getMagneticFieldTesla(tmstamp, ccdb);

    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    const int multCol = col.multV0M();
    eventHisto.fillQA(col);
    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if (part.p() > cfgCutTable->get("PartOne", "MaxP") || part.pt() > cfgCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgNspecies, kNsigma, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }
      trackHistoPartOne.fillQA(part);
    }
    for (auto& part : groupPartsTwo) {
      trackHistoPartTwo.fillQA(part);
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
      if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgNspecies, kNsigma, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      sameEventCont.setPair(p1, p2, multCol);
    }
  }

  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processMixedEvent(o2::aod::FemtoDreamCollisions& cols,
                         o2::aod::FemtoDreamParticles& parts)
  {
    BinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision1.globalIndex());
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision2.globalIndex());

      const auto& tmstamp1 = collision1.timestamp();
      const auto& magFieldTesla1 = getMagneticFieldTesla(tmstamp1, ccdb);

      const auto& tmstamp2 = collision2.timestamp();
      const auto& magFieldTesla2 = getMagneticFieldTesla(tmstamp2, ccdb);

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgNspecies, kNsigma, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
          continue;
        }
        if (ConfIsCPR) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1)) {
            continue;
          }
        }
        mixedEventCont.setPair(p1, p2, collision1.multV0M());
      }
    }
  }

  PROCESS_SWITCH(femtoDreamPairTaskTrackV0, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackV0>(cfgc),
  };
  return workflow;
}
