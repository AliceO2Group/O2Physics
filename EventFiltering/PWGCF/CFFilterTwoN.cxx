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

/// \file CFFilter.cxx
/// \brief Selection of events with triplets for femtoscopic studies
///
/// \author Bhawani Singh, TU München, bhawani.singh@cern.ch

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

#include <CCDB/BasicCCDBManager.h>
#include <fairlogger/Logger.h>
#include "DataFormatsParameters/GRPObject.h"

#include <cmath>
#include <string>

#include <bitset>
#include <vector>

namespace
{
static constexpr int nPairs{2};

enum CFPdLdTriggers {
  kPD = 0,
  kLD = 1
};

enum kDetector {
  kTPC = 0,
  kTPCTOF = 1,
  kNdetectors = 2
};

static const std::vector<std::string> CfTriggerNames{"kPD", "kLD"};
static constexpr uint8_t Track = 0; // Track
static constexpr uint8_t V0 = 1;    // V0
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

  Configurable<std::vector<float>> confKstarTriggerLimit{"KstarTriggerLimitUpper", std::vector<float>{2.0f}, "Kstar limit for selection"};
  Configurable<int> KstarTrigger{"KstarTrigger", 0, "Choice which trigger to run"};

  Configurable<float> confNSigma{"nSigma", 3.5f, "nSigma for TPC/TPC&TOF selection"};
  Configurable<float> confPIDThreshold{"PThreshold", 0.75f, "P threshold for TPC/TPC&TOF selection"};

  Configurable<float> ldeltaPhiMax{"ldeltaPhiMax", 0.010, "Max limit of delta phi"};
  Configurable<float> ldeltaEtaMax{"ldeltaEtaMax", 0.010, "Max limit of delta eta"};

  // Obtain particle candidates of protons, deuterons as well as antiprotons and antideuterons
  Partition<o2::aod::FemtoDreamParticles> partsProtonDeuteron0Part = (o2::aod::femtodreamparticle::partType == Track) && ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partsProtonDeuteron1Part = (o2::aod::femtodreamparticle::partType == Track) && ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  // FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleanerTT; Currently not used, will be needed later
  // FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleanerTV;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> closePairRejectionTT;
  // FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> closePairRejectionTV0;

  bool isPIDSelected(aod::femtodreamparticle::cutContainerType const& pidcut, std::vector<int> const& vSpecies, float nSigma, kDetector iDet = kDetector::kTPC)
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

  bool isFullPIDSelected(aod::femtodreamparticle::cutContainerType const& pidCut, std::vector<int> const& vSpecies, float nSigma, float momentum, float const PIDThreshold)
  {
    bool pidSelection = false;
    if (momentum < PIDThreshold) {
      /// TPC PID only
      pidSelection = isPIDSelected(pidCut, vSpecies, nSigma, kDetector::kTPC);
    } else {
      /// TPC + TOF PID
      pidSelection = isPIDSelected(pidCut, vSpecies, nSigma, kDetector::kTPCTOF);
    }
    return pidSelection;
  };

  void init(o2::framework::InitContext&)
  {
    bool plotPerRadii = true;

    closePairRejectionTT.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    // closePairRejectionTV0.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    registry.add("fProcessedEvents", "CF - event filtered;;events", HistType::kTH1F, {{6, -0.5, 5.5}});
    std::array<std::string, 3> eventTitles = {"all", "rejected", "p-d"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }
    registry.add("fMultiplicityBefore", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fMultiplicityAfter", "Multiplicity of events which passed pd trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fZvtxBefore", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("fZvtxAfter", "Zvtx of events which passed pd trigger", HistType::kTH1F, {{1000, -15, 15}});

    if (KstarTrigger == 0 || KstarTrigger == 11) {
      registry.add("fSameEventPartPD", "CF - same event pd distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fSameEventAntiPartPD", "CF - same event pd distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtBeforePD", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
      // registry.add("fPtAfterPD", "Transverse momentum of processed tracks which passed  selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fProtonPtAfterPD", "Transverse momentum of processed Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fDeuteronPtAfterPD", "Transverse momentum of processed Deuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});

      registry.add("fPtBeforeAntiPD", "Transverse momentum of all processed antitracks", HistType::kTH1F, {{1000, 0, 10}});
      // registry.add("fPtAfterAntiPD", "Transverse momentum  of processed antitracks passed selection", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fAntiProtonPtAfterPD", "Transverse momentum of processed Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fAntiDeuteronPtAfterPD", "Transverse momentum of processed Deuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    }
    // if (KstarTrigger == 1 || KstarTrigger == 11) {
    //   registry.add("fSameEventPartLD", "CF - same event ld distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
    //   registry.add("fSameEventAntiPartLD", "CF - same event ld distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});
    //
    //   registry.add("fPtBeforeLd", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    //   registry.add("fPtAfterLd", "Transverse momentum  of processed tracks which passed  selections", HistType::kTH1F, {{1000, 0, 10}});
    //   registry.add("fPtBeforeAntiLd", "Transverse momentum of all processed antitracks", HistType::kTH1F, {{1000, 0, 10}});
    //   registry.add("fPtAfterAntiLd", "Transverse momentum  of processed antitracks passed selection", HistType::kTH1F, {{1000, 0, 10}});
    // }

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

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassDeuteron = 1.87561f;
  // float mMassDeuteron = TDatabasePDG::Instance()->GetParticle(1000010020)->Mass();
  // float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(o2::aod::FemtoDreamCollision& col, o2::aod::FemtoDreamParticles& partsFemto)
  {
    auto partsProtonDeuteron0 = partsProtonDeuteron0Part->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsProtonDeuteron1 = partsProtonDeuteron1Part->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    auto magneticField = col.magField();
    registry.get<TH1>(HIST("fProcessedEvents"))->Fill(0);
    registry.get<TH1>(HIST("fMultiplicityBefore"))->Fill(col.multV0M());
    registry.get<TH1>(HIST("fZvtxBefore"))->Fill(col.posZ());

    for (auto p1pt : partsProtonDeuteron0) {
      registry.get<TH1>(HIST("fPtBeforePD"))->Fill(p1pt.pt());
    }

    LOGF(info, "Collision: %d", col.size());

    int prot = 0;
    int antiprot = 0;

    int deut = 0;
    int antideut = 0;

    for (auto p1pt : partsProtonDeuteron0) {
      // select protons
      if (isFullPIDSelected(p1pt.pidcut(), std::vector<int>{o2::track::PID::Proton}, confNSigma, p1pt.p(), confPIDThreshold)) {
        registry.get<TH1>(HIST("fProtonPtAfterPD"))->Fill(p1pt.pt());
        prot++;
      }
      // select deutrons
      if (isFullPIDSelected(p1pt.pidcut(), std::vector<int>{o2::track::PID::Deuteron}, confNSigma, p1pt.p(), confPIDThreshold)) {
        registry.get<TH1>(HIST("fDeuteronPtAfterPD"))->Fill(p1pt.pt());
        deut++;
      }
    }

    for (auto p1pt : partsProtonDeuteron1) {
      // select antiprotons
      if (isFullPIDSelected(p1pt.pidcut(), std::vector<int>{o2::track::PID::Proton}, confNSigma, p1pt.p(), confPIDThreshold)) {
        registry.get<TH1>(HIST("fAntiProtonPtAfterPD"))->Fill(p1pt.pt());
        antiprot++;
      }
      // select antideutrons
      if (isFullPIDSelected(p1pt.pidcut(), std::vector<int>{o2::track::PID::Deuteron}, confNSigma, p1pt.p(), confPIDThreshold)) {
        registry.get<TH1>(HIST("fAntiDeuteronPtAfterPD"))->Fill(p1pt.pt());
        antideut++;
      }
    }

    LOGF(info, "Protons: %d", prot);
    LOGF(info, "AntiProtons: %d", antideut);
    LOGF(info, "Deuterons: %d", deut);
    LOGF(info, "AntiDeuterons: %d", antideut);

    bool keepEvent[nPairs]{false};
    int lowKstarPairs[2] = {0, 0};

    if (deut > 0 && prot > 0) {

      registry.get<TH1>(HIST("fMultiplicityAfter"))->Fill(col.multV0M());
      registry.get<TH1>(HIST("fZvtxAfter"))->Fill(col.posZ());

      // TRIGGER FOR PD PAIRS
      if (KstarTrigger == 0 || KstarTrigger == 11) {

        if (partsProtonDeuteron0.size() >= 2) {
          for (auto& [p1, p2] : combinations(partsProtonDeuteron0Part, partsProtonDeuteron0Part)) {
            // would that be necessary when considering pd pairs ?
            if (!isFullPIDSelected(p1.pidcut(), std::vector<int>{o2::track::PID::Proton}, confNSigma, p1.p(), confPIDThreshold) ||
                !isFullPIDSelected(p2.pidcut(), std::vector<int>{o2::track::PID::Deuteron}, confNSigma, p2.p(), confPIDThreshold)) {
              continue;
            }
            // Think if pair cleaning is needed in current framework
            // Run close pair rejection
            if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
              continue;
            }

            auto kstar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
            registry.get<TH1>(HIST("fSameEventPartPD"))->Fill(kstar);
            if (kstar < confKstarTriggerLimit.value.at(0)) {
              lowKstarPairs[0]++;
            }
          }
        }
      }
    }
    // __________________________________________________________________________________________________________
    // TRIGGER FOR LD PAIRS
    // if (KstarTrigger == 1 || KstarTrigger == 11) {
    //   if (partsLambda0.size() >= 1 && partsDeuteron0.size() >= 1) {
    //     for (auto& partLambda : partsLambda0) {
    //       registry.get<TH1>(HIST("fPtLD"))->Fill(partLambda.pt());
    //       if (!pairCleanerTV.isCleanPair(partLambda, partLambda, partsFemto)) {
    //         continue;
    //       }
    //       for (auto& [p1, p2] : combinations(partsLambda0, partsDeuteron0)) {
    //         if (!isFullPIDSelectedProton(p1.pidcut(), p1.p()) || !isFullPIDSelectedProton(p2.pidcut(), p2.p())) {
    //           continue;
    //         }
    //         if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, getMagneticFieldTesla(tmstamp))) {
    //           continue;
    //         }
    //         auto kstar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
    //         registry.get<TH1>(HIST("fSameEventPartLd"))->Fill(kstar);
    //         if (kstar < KstarTriggerLimit.at(0)) {
    //         }
    //       }
    //     }
    //     }
    //     }

    if (lowKstarPairs[1] == 0) { // if at least one triplet found in particles, no need to check antiparticles
      if (antideut > 0 && antiprot > 0) {
        if (KstarTrigger == 0 || KstarTrigger == 11) {

          if (partsProtonDeuteron1.size() >= 2) {
            for (auto& [p1, p2] : combinations(partsProtonDeuteron1Part, partsProtonDeuteron1Part)) {
              // would that be necessary when considering pd pairs ?
              if (!isFullPIDSelected(p1.pidcut(), std::vector<int>{o2::track::PID::Proton}, confNSigma, p1.p(), confPIDThreshold) ||
                  !isFullPIDSelected(p2.pidcut(), std::vector<int>{o2::track::PID::Deuteron}, confNSigma, p2.p(), confPIDThreshold)) {
                continue;
              }
              // Think if pair cleaning is needed in current framework
              // Run close pair rejection
              if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
                continue;
              }

              auto kstar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
              registry.get<TH1>(HIST("fSameEventAntiPartPD"))->Fill(kstar);
              if (kstar < confKstarTriggerLimit.value.at(0)) {
                lowKstarPairs[1]++;
              }
            }
          }
        }
      }
    }

    if (lowKstarPairs[0] > 0) {
      keepEvent[kPD] = true;
    }

    if (lowKstarPairs[1] > 0) {
      keepEvent[kLD] = true;
    }

    tags(keepEvent[kPD], keepEvent[kLD]);

    if (keepEvent[kPD] || keepEvent[kLD]) {
      registry.get<TH1>(HIST("fProcessedEvents"))->Fill(2);
    } else {
      registry.get<TH1>(HIST("fProcessedEvents"))->Fill(1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilterTwoN>(cfg)};
}
