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

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamPairCleaner.h"
#include "FemtoDreamContainer.h"
#include "FemtoDreamDetaDphiStar.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 7;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MinPt", "MaxPt", "MaxEta", "MaxDCAxy", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF"};
static const float cutsTable[nPart][nCuts]{
  {0.5f, 4.05f, 0.8f, 0.1f, 0.75f, 3.f, 3.f},
  {0.5f, 4.05f, 0.8f, 0.1f, 0.75f, 3.f, 3.f}};

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
  Configurable<int> ConfCutPartOne{"ConfCutPartOne", 693318, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FemtoDreamParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 ((aod::femtodreamparticle::cut & aod::femtodreamparticle::cutContainerType(ConfCutPartOne)) == aod::femtodreamparticle::cutContainerType(ConfCutPartOne)) &&
                                                 (aod::femtodreamparticle::pt > cfgCutTable->get("PartOne", "MinPt")) &&
                                                 (aod::femtodreamparticle::pt < cfgCutTable->get("PartOne", "MaxPt")) &&
                                                 (nabs(aod::femtodreamparticle::tempFitVar) < cfgCutTable->get("PartOne", "MaxDCAxy")) &&
                                                 (nabs(aod::femtodreamparticle::eta) < cfgCutTable->get("PartOne", "MaxEta"));

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<int> ConfCutPartTwo{"ConfCutPartTwo", 693318, "Particle 2 - Selection bit"};
  Configurable<std::vector<int>> ConfPIDPartTwo{"ConfPIDPartTwo", std::vector<int>{2}, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2
  Partition<aod::FemtoDreamParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 ((aod::femtodreamparticle::cut & aod::femtodreamparticle::cutContainerType(ConfCutPartTwo)) == aod::femtodreamparticle::cutContainerType(ConfCutPartTwo)) &&
                                                 (aod::femtodreamparticle::pt > cfgCutTable->get("PartTwo", "MinPt")) &&
                                                 (aod::femtodreamparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                                 (nabs(aod::femtodreamparticle::tempFitVar) < cfgCutTable->get("PartTwo", "MaxDCAxy")) &&
                                                 (nabs(aod::femtodreamparticle::eta) < cfgCutTable->get("PartTwo", "MaxEta"));

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne, vPIDPartTwo;

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {26, 1, 27}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<float> ConfBField{"ConfBField", +0.5, "Magnetic Field"};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
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
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfBField, false); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;
    vPIDPartTwo = ConfPIDPartTwo;
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
    const int multCol = col.multV0M();
    /// Histogramming same event
    for (auto& part : partsOne) {
      if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }
      trackHistoPartOne.fillQA(part);
    }
    if (!ConfIsSame) {
      for (auto& part : partsTwo) {
        if (!isFullPIDSelected(part.pidcut(), part.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartTwo.fillQA(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(partsOne, partsTwo)) {
      if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF")) || !isFullPIDSelected(p2.pidcut(), p2.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      /// close pair rejection
      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts)) {
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
    cols.bindExternalIndices(&parts);
    auto particlesTuple = std::make_tuple(parts);
    AnalysisDataProcessorBuilder::GroupSlicer slicer(cols, particlesTuple);

    for (auto& [collision1, collision2] : soa::selfCombinations("fBin", ConfNEventsMix, -1, soa::join(hashes, cols), soa::join(hashes, cols))) {
      auto it1 = slicer.begin();
      auto it2 = slicer.begin();
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision2.index()) {
          it2 = slice;
          break;
        }
      }

      auto particles1 = std::get<aod::FemtoDreamParticles>(it1.associatedTables());
      particles1.bindExternalIndices(&cols);
      auto particles2 = std::get<aod::FemtoDreamParticles>(it2.associatedTables());
      particles2.bindExternalIndices(&cols);

      partsOne.bindTable(particles1);
      partsTwo.bindTable(particles2);

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

      for (auto& [p1, p2] : combinations(partsOne, partsTwo)) {
        if (!isFullPIDSelected(p1.pidcut(), p1.p(), cfgCutTable->get("PartOne", "PIDthr"), vPIDPartOne, cfgCutTable->get("PartOne", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF")) || !isFullPIDSelected(p2.pidcut(), p2.p(), cfgCutTable->get("PartTwo", "PIDthr"), vPIDPartTwo, cfgCutTable->get("PartTwo", "nSigmaTPC"), cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
          continue;
        }
        if (ConfIsCPR) {
          if (pairCloseRejection.isClosePair(p1, p2, parts)) {
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
