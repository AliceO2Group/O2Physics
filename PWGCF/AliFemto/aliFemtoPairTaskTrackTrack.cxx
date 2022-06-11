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

/// \file aliFemtoPairTaskTrackTrack.cxx
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

//#include "AliFemtoDerived.h"
#include "PWGCF/DataModel/FemtoDerived.h"
#include "AliFemtoParticleHisto.h" //new AliFemto name
#include "AliFemtoEventHisto.h" //new AliFemto name
#include "AliFemtoPairCleaner.h" //new AliFemto name
#include "AliFemtoContainer.h" //new AliFemto name
#include "AliFemtoDetaDphiStar.h" //new AliFemto name
#include "AliFemtoUtils.h" //new AliFemto name

using namespace o2;
using namespace o2::analysis::aliFemto; //new AliFemto name
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
  {4.05f, 1.f, 3.f, 3.f, 100.f}};

static const std::vector<float> kNsigma = {3.5f, 3.f, 2.5f};

} // namespace

struct aliFemtoPairTaskTrackTrack { //new AliFemto name

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};

  /// Particle 1
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 84035877, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1 //new AliFemto name
  Partition<aod::FemtoDreamParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for particle 1 //new AliFemto name
  AliFemtoParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 84035877, "Particle 2 - Selection bit"};
  Configurable<std::vector<int>> ConfPIDPartTwo{"ConfPIDPartTwo", std::vector<int>{2}, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2 //new AliFemto name
  Partition<aod::FemtoDreamParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 //  (aod::femtodreamparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                                 ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);

  /// Histogramming for particle 2 //new AliFemto name
  AliFemtoParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event //new AliFemto name
  AliFemtoEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne, vPIDPartTwo;

  /// Correlation part
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  //new AliFemto name
  AliFemtoContainer<aliFemtoContainer::EventType::same, aliFemtoContainer::Observable::kstar> sameEventCont;
  //new AliFemto name
  AliFemtoContainer<aliFemtoContainer::EventType::mixed, aliFemtoContainer::Observable::kstar> mixedEventCont;
  AliFemtoPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  AliFemtoDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry);
    if (!ConfIsSame) {
      trackHistoPartTwo.init(&qaRegistry);
    }

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;
    vPIDPartTwo = ConfPIDPartTwo;

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
    //new AliFemto name
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    //new AliFemto name
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
    if (!ConfIsSame) {
      for (auto& part : groupPartsTwo) {
        if (part.p() > cfgCutTable->get("PartTwo", "MaxP") || part.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgNspecies, kNsigma, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartTwo.fillQA(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
      if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.p() > cfgCutTable->get("PartTwo", "MaxP") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgNspecies, kNsigma, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF")) || !isFullPIDSelected(p2.pidcut(), p2.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgNspecies, kNsigma, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
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

  PROCESS_SWITCH(aliFemtoPairTaskTrackTrack, processSameEvent, "Enable processing same event", true); //new AliFemto name

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processMixedEvent(o2::aod::FemtoDreamCollisions& cols,
                         o2::aod::FemtoDreamParticles& parts)
  {

    BinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      //new AliFemto name
      auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision1.globalIndex());
      //new AliFemto name
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision2.globalIndex());

      const auto& tmstamp1 = collision1.timestamp();
      const auto& magFieldTesla1 = getMagneticFieldTesla(tmstamp1, ccdb);

      const auto& tmstamp2 = collision2.timestamp();
      const auto& magFieldTesla2 = getMagneticFieldTesla(tmstamp2, ccdb);

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.p() > cfgCutTable->get("PartTwo", "MaxP") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgNspecies, kNsigma, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF")) || !isFullPIDSelected(p2.pidcut(), p2.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgNspecies, kNsigma, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
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

  PROCESS_SWITCH(aliFemtoPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true); //new AliFemto name
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<aliFemtoPairTaskTrackTrack>(cfgc), //new AliFemto name
  };
  return workflow;
}
