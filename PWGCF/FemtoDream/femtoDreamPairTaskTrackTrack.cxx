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

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 4;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 0.75f, 3.f, 3.f},
  {4.05f, 0.75f, 3.f, 3.f}};

static constexpr int nNsigma = 3;
static constexpr float kNsigma[nNsigma] = {3.5f, 3.f, 2.5f};

enum kPIDselection {
  k3d5sigma = 0,
  k3sigma = 1,
  k2d5sigma = 2
};

enum kDetector {
  kTPC = 0,
  kTPCTOF = 1,
  kNdetectors = 2
};
} // namespace

struct femtoDreamPairTaskTrackTrack {

  /// Particle selection part
  // uint trackTypeSel = aod::femtodreamparticle::ParticleType::kTrack; // \todo at some point filters will be able to cope with enums
  // \todo checky why is the enum not working in the partition

  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};

  /// Particle 1
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FemtoDreamParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (aod::femtodreamparticle::pt < cfgCutTable->get("PartOne", "MaxPt")) &&
                                                 ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 5542474, "Particle 2 - Selection bit"};
  Configurable<std::vector<int>> ConfPIDPartTwo{"ConfPIDPartTwo", std::vector<int>{2}, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2
  Partition<aod::FemtoDreamParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 (aod::femtodreamparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                                 ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne, vPIDPartTwo;

  /// Correlation part
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
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
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, false); /// \todo add config for Δη and ΔΦ cut values
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

  /// Function to retrieve the nominal mgnetic field in kG (0.1T) and convert it directly to T
  float getMagneticFieldTesla(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = 0.1 * (grpo->getNominalL3Field());
    return output;
  }

  /// internal function that returns the kPIDselection element corresponding to a specifica n-sigma value
  /// \param number of sigmas for PIF
  /// \return kPIDselection corresponing to n-sigma
  kPIDselection getPIDselection(float nSigma)
  {
    for (int i = 0; i < nNsigma; i++) {
      if (abs(nSigma - kNsigma[i]) < 1e-3) {
        return static_cast<kPIDselection>(i);
      }
    }
    LOG(info) << "Invalid value of nSigma: " << nSigma << ". Standard 3 sigma returned." << std::endl;
    return kPIDselection::k3sigma;
  }

  /// function that checks whether the PID selection specified in the vectors is fulfilled
  /// \param pidcut Bit-wise container for the PID
  /// \param vSpecies Vector with the different selections
  /// \return Whether the PID selection specified in the vectors is fulfilled
  bool isPIDSelected(aod::femtodreamparticle::cutContainerType const& pidcut, std::vector<int> const& vSpecies, float nSigma, kDetector iDet = kDetector::kTPC)
  {
    bool pidSelection = true;
    kPIDselection kNsigma = getPIDselection(nSigma);
    for (auto iSpecies : vSpecies) {
      //\todo we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
      // if (!((pidcut >> it.first) & it.second)) {
      int bit_to_check = cfgNspecies * kDetector::kNdetectors * kNsigma + iSpecies * kDetector::kNdetectors + iDet;
      if (!(pidcut & (1UL << bit_to_check))) {
        pidSelection = false;
      }
    }
    return pidSelection;
  };

  /// function that checks whether the PID selection specified in the vectors is fulfilled, depending on the momentum TPC or TPC+TOF PID is conducted
  /// \param pidcut Bit-wise container for the PID
  /// \param mom Momentum of the track
  /// \param pidThresh Momentum threshold that separates between TPC and TPC+TOF PID
  /// \param vecTPC Vector with the different selections for the TPC PID
  /// \param vecComb Vector with the different selections for the TPC+TOF PID
  /// \return Whether the PID selection is fulfilled
  bool isFullPIDSelected(aod::femtodreamparticle::cutContainerType const& pidCut, float const& momentum, float const& pidThresh, std::vector<int> const& vSpecies, float nSigmaTPC = 3.5, float nSigmaTPCTOF = 3.5)
  {
    bool pidSelection = true;
    if (momentum < pidThresh) {
      /// TPC PID only
      pidSelection = isPIDSelected(pidCut, vSpecies, nSigmaTPC, kDetector::kTPC);
    } else {
      /// TPC + TOF PID
      pidSelection = isPIDSelected(pidCut, vSpecies, nSigmaTPCTOF, kDetector::kTPCTOF);
    }
    return pidSelection;
  };

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processSameEvent(o2::aod::FemtoDreamCollision& col,
                        o2::aod::FemtoDreamParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    const int multCol = col.multV0M();
    eventHisto.fillQA(col);
    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if (part.pt() > cfgCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }
      trackHistoPartOne.fillQA(part);
    }
    if (!ConfIsSame) {
      for (auto& part : groupPartsTwo) {
        if (part.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartTwo.fillQA(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
      if (p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF")) || !isFullPIDSelected(p2.pidcut(), p2.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      auto tmstamp = col.timestamp();
      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, getMagneticFieldTesla(tmstamp))) {
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

  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processMixedEvent(o2::aod::FemtoDreamCollisions& cols,
                         o2::aod::Hashes& hashes,
                         o2::aod::FemtoDreamParticles& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations("fBin", ConfNEventsMix, -1, soa::join(hashes, cols), soa::join(hashes, cols))) {

      auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision1.globalIndex());
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision2.globalIndex());

      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      /// could work like that, but only if PID is contained within the partitioning!
      // auto particlesEvent1 = std::get<aod::FemtoDreamParticles>(it1.associatedTables());
      // particlesEvent1.bindExternalIndices(&cols);
      // auto particlesEvent2 = std::get<aod::FemtoDreamParticles>(it2.associatedTables());
      // particlesEvent2.bindExternalIndices(&cols);
      /// for the x-check
      // partsOne.bindTable(particlesEvent2);
      // auto nPart1Evt2 = partsOne.size();
      // partsTwo.bindTable(particlesEvent1);
      // auto nPart2Evt1 = partsTwo.size();
      /// for actual event mixing
      // partsOne.bindTable(particlesEvent1);
      // partsTwo.bindTable(particlesEvent2);
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF")) || !isFullPIDSelected(p2.pidcut(), p2.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }

        if (ConfIsCPR) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, getMagneticFieldTesla(collision1.timestamp()))) {
            continue;
          }
        }
        mixedEventCont.setPair(p1, p2, collision1.multV0M()); // < \todo dirty trick, the multiplicity will be of course within the bin width used for the hashes
      }
    }
  }

  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
