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

/// \file CFFilterTwoN.cxx
/// \brief Selection of events with different kind of pairs for femtoscopic studies
///
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "../filterTables.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/FemtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/FemtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/FemtoDreamContainer.h"
#include "PWGCF/FemtoDream/FemtoDreamMath.h"
#include "PWGCF/FemtoDream/FemtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/FemtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/FemtoDreamContainer.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/PhysicsConstants.h"

#include <vector>

namespace
{

enum kCFTwoBodyTriggers {
  kPD, //=0
  kLD, //=1
  kLAST_CFTwoBodyTriggers
};

enum kDetector {
  kTPC,
  kTPCTOF,
  kNdetectors,
  kLAST_Detector
};

static const std::vector<std::string> CfTriggerNames{"kPD", "kLD"};
static constexpr uint8_t Track = 0;
static constexpr uint8_t V0 = 1; // V0
// static constexpr uint8_t V0Daughter = 2; // V0  daughters

static constexpr uint32_t kSignMinusMask = 1;
static constexpr uint32_t kSignPlusMask = 2;

// static constexpr uint32_t knSigmaProton = 48;
static constexpr uint32_t kValue0 = 0;

} // namespace

namespace o2::aod
{
using FullCollision = soa::Join<aod::Collisions,
                                aod::EvSels,
                                aod::Mults>::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct CFFilterTwoN {

  Produces<aod::CFFiltersTwoN> tags;

  Configurable<float> confKstarTriggerLimit{"KstarTriggerLimitUpper", 1.4f, "Kstar limit for selection"};
  Configurable<int> KstarTrigger{"KstarTrigger", 0, "Choice which trigger to run"};
  Configurable<float> confPIDThreshold{"PThreshold", 0.75f, "P threshold for TPC/TPC&TOF selection"};
  Configurable<float> ldeltaPhiMax{"ldeltaPhiMax", 0.010, "Max limit of delta phi"};
  Configurable<float> ldeltaEtaMax{"ldeltaEtaMax", 0.010, "Max limit of delta eta"};

  // obtain particle candidates of protons, deuterons as well as antiprotons and antideuterons
  Partition<o2::aod::FemtoDreamParticles> partPD = (o2::aod::femtodreamparticle::partType == Track) &&
                                                   ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partAntiPD = (o2::aod::femtodreamparticle::partType == Track) &&
                                                       ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);
  // obtain lambdas and antilambdas
  Partition<o2::aod::FemtoDreamParticles> partL = (o2::aod::femtodreamparticle::partType == V0) &&
                                                  ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partAntiL = (o2::aod::femtodreamparticle::partType == V0) &&
                                                      ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  // containers for close pair rejection
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> closePairRejectionTT;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> closePairRejectionTV0;

  // helper function for checking paritcle pid
  bool PIDHelper(aod::femtodreamparticle::cutContainerType const& pidcut, int vSpecies, kDetector iDet = kDetector::kTPC)
  {
    bool pidSelection = true;
    int bit_to_check = vSpecies * kDetector::kNdetectors + iDet;
    if (!(pidcut & (1UL << bit_to_check))) {
      pidSelection = false;
    }
    return pidSelection;
  };

  // function for checking paritcle pid
  bool SelectParticlePID(aod::femtodreamparticle::cutContainerType const& pidCut, int vSpecies, float momentum, float const PIDThreshold)
  {
    bool pidSelection = false;
    if (momentum < PIDThreshold) {
      /// TPC PID only
      pidSelection = PIDHelper(pidCut, vSpecies, kDetector::kTPC);
    } else {
      /// TPC + TOF PID
      pidSelection = PIDHelper(pidCut, vSpecies, kDetector::kTPCTOF);
    }
    return pidSelection;
  };

  void init(o2::framework::InitContext&)
  {
    registry.add("fProcessedEvents", "CF Two Body - event filtered;;events}", HistType::kTH1F, {{2 + kLAST_CFTwoBodyTriggers, 0, 2 + kLAST_CFTwoBodyTriggers}});

    std::array<std::string, 2 + kLAST_CFTwoBodyTriggers> eventTitles = {"all", "rejected", "p-d", "l-d"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    registry.add("fMultiplicityBefore", "Multiplicity before trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fMultiplicityAfter", "Multiplicity after trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fZvtxBefore", "Zvtx before trigger", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("fZvtxAfter", "Zvtx after trigger", HistType::kTH1F, {{1000, -15, 15}});

    registry.add("fPtBeforeSel", "Transverse momentum of positive tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtAntiBeforeSel", "Transverse momentum of negative tracks", HistType::kTH1F, {{1000, 0, 10}});

    bool plotPerRadii = true;

    if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
      registry.add("fKstarPD", "CF - same event pd distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiPD", "CF - same event pd distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fPtAntiProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});

      closePairRejectionTT.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      registry.add("fKstarLD", "CF - same event ld distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiLD", "CF - same event ld distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtLambdaAfterSel", "Transverse momentum of Lambdas which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fPtAntiLambdaAfterSel", "Transverse momentum of Antidlambdas which passed selections", HistType::kTH1F, {{1000, 0, 10}});

      closePairRejectionTV0.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    }

    registry.add("fPtDeuteronAfterSel", "Transverse momentum of Deuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtAntiDeuteronAfterSel", "Transverse momentum of Antideuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
  }

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassDeuteron = o2::constants::physics::MassDeuteron;
  float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(o2::aod::FemtoDreamCollision& col, o2::aod::FemtoDreamParticles& partsFemto)
  {
    // get partitions of all paritcles and antiparticles
    auto partsPD = partPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsAntiPD = partAntiPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    // get partions of V0s
    auto partsL = partL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsAntiL = partAntiL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    // magnetic field is need for close pair rejection
    auto magneticField = col.magField();

    registry.fill(HIST("fProcessedEvents"), 0);
    registry.fill(HIST("fMultiplicityBefore"), col.multV0M());
    registry.fill(HIST("fZvtxBefore"), col.posZ());

    // pass through the particles once to check if there are any particles of interest in the first place
    int Nproton = 0;
    int Nantiproton = 0;
    int Nlambda = 0;
    int Nantilambda = 0;
    int Ndeuteron = 0;
    int Nantideuteron = 0;

    for (auto pd : partsPD) {
      registry.fill(HIST("fPtBeforeSel"), pd.pt());
      // select protons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
        if (SelectParticlePID(pd.pidcut(), o2::track::PID::Proton, pd.p(), confPIDThreshold)) {
          registry.fill(HIST("fPtProtonAfterSel"), pd.pt());
          Nproton++;
        }
      }
      // select deuterons
      if (SelectParticlePID(pd.pidcut(), o2::track::PID::Deuteron, pd.p(), confPIDThreshold)) {
        registry.fill(HIST("fPtDeuteronAfterSel"), pd.pt());
        Ndeuteron++;
      }
    }
    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      // select lambdas
      for (auto lambda : partsL) {
        registry.fill(HIST("fPtLambdaAfterSel"), lambda.pt());
        Nlambda++;
      }
    }

    for (auto antipd : partsAntiPD) {
      registry.fill(HIST("fPtAntiBeforeSel"), antipd.pt());
      // select antiprotons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
        if (SelectParticlePID(antipd.pidcut(), o2::track::PID::Proton, antipd.p(), confPIDThreshold)) {
          registry.fill(HIST("fPtAntiProtonAfterSel"), antipd.pt());
          Nantiproton++;
        }
      }
      // select antideuterons
      if (SelectParticlePID(antipd.pidcut(), o2::track::PID::Deuteron, antipd.p(), confPIDThreshold)) {
        registry.fill(HIST("fPtAntiDeuteronAfterSel"), antipd.pt());
        Nantideuteron++;
      }
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      for (auto antilambda : partsAntiL) {
        // select antilambdas
        registry.fill(HIST("fPtAntiLambdaAfterSel"), antilambda.pt());
        Nlambda++;
      }
    }

    bool keepEvent[kLAST_CFTwoBodyTriggers] = {false, false};
    int lowKstarPairs[kLAST_CFTwoBodyTriggers] = {0, 0};

    bool pdPair = false;
    bool dpPair = false;
    double kStar = 0.;

    // trigger for pd pairs
    if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
      if (Ndeuteron > 0 && Nproton > 0) {
        // loop over all unique combinations of particles, excluding the self combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsPD, partsPD))) {

          // check if it is a pd pair
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Proton, p1.p(), confPIDThreshold.value) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Deuteron, p2.p(), confPIDThreshold.value)) {
            pdPair = true;
          } else {
            pdPair = false;
          }

          // check if it is dp pair
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), confPIDThreshold.value) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Proton, p2.p(), confPIDThreshold.value)) {
            dpPair = true;
          } else {
            dpPair = false;
          }

          // if neither is the case, skip
          if (!(pdPair || dpPair)) {
            continue;
          }

          // reject close pairs
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }

          // compute kstar depending on the pairing
          if (pdPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
          } else if (dpPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassDeuteron, p2, mMassProton);
          } else {
            kStar = confKstarTriggerLimit;
          }
          // check if the kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[kPD]++;
            registry.fill(HIST("fKstarPD"), kStar);
          }
        }
      }

      if (Nantideuteron > 0 && Nantiproton > 0) {
        // loop over all unique combinations of antiparticles, excluding the self combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsAntiPD, partsAntiPD))) {

          // check if it is a (anti)pd pair
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Proton, p1.p(), confPIDThreshold.value) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Deuteron, p2.p(), confPIDThreshold.value)) {
            pdPair = true;
          } else {
            pdPair = false;
          }

          // check if it is (anti)dp pair
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), confPIDThreshold.value) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Proton, p2.p(), confPIDThreshold.value)) {
            dpPair = true;
          } else {
            dpPair = false;
          }

          // if neither is the case, skip
          if (!(pdPair || dpPair)) {
            continue;
          }

          // reject close pairs
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }

          // compute kstar depending on the pairing
          if (pdPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
          } else if (dpPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassDeuteron, p2, mMassProton);
          } else {
            kStar = confKstarTriggerLimit;
          }

          // check if the kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[kPD]++;
            registry.fill(HIST("fKstarAntiPD"), kStar);
          }
        }
      }
    }

    // trigger for ld pairs
    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      if (Ndeuteron > 0 && Nlambda > 0) {
        // loop over all unique combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsUpperIndexPolicy(partsPD, partsL))) {
          // check if the particle is a deuteron
          // we do not need to check the V0s
          if (!SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), confPIDThreshold.value)) {
            continue;
          }
          if (closePairRejectionTV0.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          kStar = FemtoDreamMath::getkstar(p1, mMassLambda, p2, mMassDeuteron);
          // check if kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[1]++;
            registry.fill(HIST("fKstarLD"), kStar);
          }
        }
      }
      if (Nantideuteron > 0 && Nantilambda > 0) {
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsAntiPD, partsAntiL))) {
          // check if the particle is a antideuteron
          if (!SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), confPIDThreshold.value)) {
            continue;
          }
          if (closePairRejectionTV0.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassLambda, p2, mMassDeuteron);
          // check if kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[1]++;
            registry.fill(HIST("fKstarAntiLD"), kstar);
          }
        }
      }
    }

    // if we found any pair below the kstar limit, keep the event
    if (lowKstarPairs[kPD] > 0) {
      keepEvent[kPD] = true;
      registry.fill(HIST("fProcessedEvents"), 2);
    }
    if (lowKstarPairs[kLD] > 0) {
      keepEvent[kLD] = true;
      registry.fill(HIST("fProcessedEvents"), 3);
    }

    // fill table for the trigger
    tags(keepEvent[kPD], keepEvent[kLD]);

    if (keepEvent[kPD] > 0 || keepEvent[kLD] > 0) {
      registry.fill(HIST("fMultiplicityAfter"), col.multV0M());
      registry.fill(HIST("fZvtxAfter"), col.posZ());
    } else {
      registry.fill(HIST("fProcessedEvents"), 1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilterTwoN>(cfg)};
}
