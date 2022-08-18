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

// #include <CCDB/BasicCCDBManager.h>
#include <fairlogger/Logger.h>
#include <math.h>
#include "DataFormatsParameters/GRPObject.h"

#include <cmath>
#include <string>

#include <bitset>
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

  Configurable<float> confNSigma{"nSigma", 3.5f, "nSigma for TPC/TPC&TOF selection"};
  Configurable<float> confPIDThreshold{"PThreshold", 0.75f, "P threshold for TPC/TPC&TOF selection"};

  Configurable<float> ldeltaPhiMax{"ldeltaPhiMax", 0.010, "Max limit of delta phi"};
  Configurable<float> ldeltaEtaMax{"ldeltaEtaMax", 0.010, "Max limit of delta eta"};

  // Obtain particle candidates of protons, deuterons as well as antiprotons and antideuterons
  Partition<o2::aod::FemtoDreamParticles> partPD = (o2::aod::femtodreamparticle::partType == Track) &&
                                                   ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partAntiPD = (o2::aod::femtodreamparticle::partType == Track) &&
                                                       ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partL = (o2::aod::femtodreamparticle::partType == V0) &&
                                                  ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partAntiL = (o2::aod::femtodreamparticle::partType == V0) &&
                                                      ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> closePairRejectionTT;
  // FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> closePairRejectionTV0;

  bool PIDHelper(aod::femtodreamparticle::cutContainerType const& pidcut, std::vector<int> const& vSpecies, float nSigma, kDetector iDet = kDetector::kTPC)
  {
    bool pidSelection = true;
    for (auto iSpecies : vSpecies) {
      int bit_to_check = iSpecies * kDetector::kNdetectors + iDet;
      if (!(pidcut & (1UL << bit_to_check))) {
        pidSelection = false;
      }
    }
    return pidSelection;
  };

  bool SelectParticlePID(aod::femtodreamparticle::cutContainerType const& pidCut, std::vector<int> const& vSpecies, float nSigma, float momentum, float const PIDThreshold)
  {
    bool pidSelection = false;
    if (momentum < PIDThreshold) {
      /// TPC PID only
      pidSelection = PIDHelper(pidCut, vSpecies, nSigma, kDetector::kTPC);
    } else {
      /// TPC + TOF PID
      pidSelection = PIDHelper(pidCut, vSpecies, nSigma, kDetector::kTPCTOF);
    }
    return pidSelection;
  };

  void init(o2::framework::InitContext&)
  {

    bool plotPerRadii = true;

    closePairRejectionTT.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    // closePairRejectionTV0.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);

    registry.add("fProcessedEvents", "CF Two Body - event filtered;;events}", HistType::kTH1F, {{2 + kLAST_CFTwoBodyTriggers, 0, 2 + kLAST_CFTwoBodyTriggers}});

    std::array<std::string, 2 + kLAST_CFTwoBodyTriggers> eventTitles = {"all", "rejected", "p-d", "L-d"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    registry.add("fMultiplicityBefore", "Multiplicity before trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fMultiplicityAfter", "Multiplicity after trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fZvtxBefore", "Zvtx before trigger", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("fZvtxAfter", "Zvtx after trigger", HistType::kTH1F, {{1000, -15, 15}});

    registry.add("fPtBeforeSel", "Transverse momentum of positive tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtAntiBeforeSel", "Transverse momentum of negative tracks", HistType::kTH1F, {{1000, 0, 10}});

    if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
      LOGF(info, "PD Trigger selected");
      registry.add("fKstarPD", "CF - same event pd distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiPD", "CF - same event pd distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fPtAntiProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      LOGF(info, "LD Trigger selected");
      registry.add("fKstarLD", "CF - same event ld distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiLD", "CF - same event ld distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtLambdaAfterSel", "Transverse momentum of Lambdas which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fPtAntiLambdaAfterSel", "Transverse momentum of Antidlambdas which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    }

    registry.add("fPtDeuteronAfterSel", "Transverse momentum of Deuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtAntiDeuteronAfterSel", "Transverse momentum of Antideuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
  }

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassDeuteron = 1.87561f; // get mass from pdg? TDatabasePDG::Instance()->GetParticle(1000010020)->Mass()
  float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(o2::aod::FemtoDreamCollision& col, o2::aod::FemtoDreamParticles& partsFemto)
  {
    auto partsPD = partPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsAntiPD = partAntiPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    auto partsL = partL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsAntiL = partAntiL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    auto magneticField = col.magField();
    registry.fill(HIST("fProcessedEvents"), 0);
    registry.fill(HIST("fMultiplicityBefore"), col.multV0M());
    registry.fill(HIST("fZvtxBefore"), col.posZ());

    LOGF(info, "Collision: %d/%d", col.index(), col.size());

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
        if (SelectParticlePID(pd.pidcut(), std::vector<int>{o2::track::PID::Proton}, confNSigma, pd.p(), confPIDThreshold)) {
          registry.fill(HIST("fPtProtonAfterSel"), pd.pt());
          Nproton++;
        }
      }
      // select deuterons
      if (SelectParticlePID(pd.pidcut(), std::vector<int>{o2::track::PID::Deuteron}, confNSigma, pd.p(), confPIDThreshold)) {
        registry.fill(HIST("fPtDeuteronAfterSel"), pd.pt());
        Ndeuteron++;
      }
    }
    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      for (auto lambda : partsL) {
        // select lambdas
        if (SelectParticlePID(lambda.pidcut(), std::vector<int>{o2::track::PID::Lambda}, confNSigma, lambda.p(), confPIDThreshold)) {
          registry.fill(HIST("fPtLambdaAfterSel"), lambda.pt());
          Nlambda++;
        }
      }
    }

    for (auto antipd : partsAntiPD) {
      registry.fill(HIST("fPtAntiBeforeSel"), antipd.pt());
      // select antiprotons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
        if (SelectParticlePID(antipd.pidcut(), std::vector<int>{o2::track::PID::Proton}, confNSigma, antipd.p(), confPIDThreshold)) {
          registry.fill(HIST("fPtAntiProtonAfterSel"), antipd.pt());
          Nantiproton++;
        }
      }
      // select antideuterons
      if (SelectParticlePID(antipd.pidcut(), std::vector<int>{o2::track::PID::Deuteron}, confNSigma, antipd.p(), confPIDThreshold)) {
        registry.fill(HIST("fPtAntiDeuteronAfterSel"), antipd.pt());
        Nantideuteron++;
      }
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      for (auto antilambda : partsAntiL) {
        // select antilambdas
        if (SelectParticlePID(antilambda.pidcut(), std::vector<int>{o2::track::PID::Lambda}, confNSigma, antilambda.p(), confPIDThreshold)) {
          registry.fill(HIST("fPtAntiLambdaAfterSel"), antilambda.pt());
          Nlambda++;
        }
      }
    }

    bool keepEvent[kLAST_CFTwoBodyTriggers] = {false, false};
    int lowKstarPairs[kLAST_CFTwoBodyTriggers] = {0, 0};

    // trigger for pd pairs
    if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
      if (Ndeuteron > 0 && Nproton > 0) {
        for (auto& [p1, p2] : combinations(partsPD, partsPD)) {
          if (!SelectParticlePID(p1.pidcut(), std::vector<int>{o2::track::PID::Proton, o2::track::PID::Deuteron}, confNSigma, p1.p(), confPIDThreshold)) {
            continue;
          }
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
          if (kstar < confKstarTriggerLimit.value) {
            lowKstarPairs[kPD]++;
            registry.fill(HIST("fKstarPD"), kstar);
          }
        }
      }
      if (Nantideuteron > 0 && Nantiproton > 0) {
        for (auto& [p1, p2] : combinations(partsAntiPD, partsAntiPD)) {
          if (!SelectParticlePID(p1.pidcut(), std::vector<int>{o2::track::PID::Proton, o2::track::PID::Deuteron}, confNSigma, p1.p(), confPIDThreshold)) {
            continue;
          }
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
          if (kstar < confKstarTriggerLimit.value) {
            lowKstarPairs[kPD]++;
            registry.fill(HIST("fKstarAntiPD"), kstar);
          }
        }
      }
    }

    // trigger for ld pairs
    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      if (Ndeuteron > 0 && Nlambda > 0) {
        for (auto& [p1, p2] : combinations(partsL, partsPD)) {
          if (!SelectParticlePID(p1.pidcut(), std::vector<int>{o2::track::PID::Lambda, o2::track::PID::Deuteron}, confNSigma, p1.p(), confPIDThreshold)) {
            continue;
          }
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassLambda, p2, mMassDeuteron);
          if (kstar < confKstarTriggerLimit.value) {
            lowKstarPairs[1]++;
            registry.fill(HIST("fKstarLD"), kstar);
          }
        }
      }
      if (Nantideuteron > 0 && Nantilambda > 0) {
        for (auto& [p1, p2] : combinations(partsAntiL, partsAntiPD)) {
          if (!SelectParticlePID(p1.pidcut(), std::vector<int>{o2::track::PID::Lambda, o2::track::PID::Deuteron}, confNSigma, p1.p(), confPIDThreshold)) {
            continue;
          }
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassLambda, p2, mMassDeuteron);
          if (kstar < confKstarTriggerLimit.value) {
            lowKstarPairs[1]++;
            registry.fill(HIST("fKstarAntiLD"), kstar);
          }
        }
      }
    }

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
