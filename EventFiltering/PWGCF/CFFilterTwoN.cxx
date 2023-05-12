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

#include <Framework/Configurable.h>
#include <TMath.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "../filterTables.h"
#include "../../PWGCF/FemtoDream/FemtoUtils.h"

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

static const std::vector<std::string> CfTriggerNames{"kPD", "kLD"};
static constexpr uint8_t Track = 0;
static constexpr uint8_t V0 = 1;
// static constexpr uint8_t V0Daughter = 2; // V0  daughters

static constexpr uint32_t kSignMinusMask = 1;
static constexpr uint32_t kSignPlusMask = 2;

// static constexpr uint32_t knSigmaProton = 48;
static constexpr uint32_t kValue0 = 0;

} // namespace

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct CFFilterTwoN {
  SliceCache cache;

  Produces<aod::CFFiltersTwoN> tags;

  Configurable<float> confKstarTriggerLimit{"KstarTriggerLimitUpper", 1.0f, "Kstar limit for selection"};
  Configurable<int> KstarTrigger{"KstarTrigger", 0, "Choice which trigger to run"};
  // 0 for pd
  // 1 for ld
  // 10 for both

  Configurable<float> confProtonPtMin{"ProtonPtMin", 0.5, "Minimal Pt for Protons"};
  Configurable<float> confProtonPtMax{"ProtonPtMax", 4.05, "Maximal Pt for Protons"};
  Configurable<float> confPIDThreshold{"PThreshold", 0.75f, "P threshold for TPC/TPC&TOF selection (Protons only)"};

  Configurable<float> confDeuteronPtMin{"DeuteronPtMin", 0.5, "Minimal Pt for Deuterons"};
  Configurable<float> confDeuteronPtMax{"DeuteronPtMax", 1.4, "Maximal Pt for Deuterons"};

  Configurable<float> confPIDnSigmaTPCAcceptance{"PIDnSigmaTPCPIDAcceptance",
                                                 3.,
                                                 "nSigmaTPC for accepting Protons and Deuterons (this value needs to be listed in PIDnSgimaTPCMax)"};
  // the value used in this configurable needs to one of the values listed in ConfPIDnSigmaTPCMax!
  Configurable<float> confPIDnSigmTPCTOFAcceptance{"PIDnSigmaTPCTOFPIDAcceptance",
                                                   3.,
                                                   "nSigmaTPCTOF for accepting Protons (this values needs to be listed in PIDnSigmaTPCMax)"};
  // the value used in this configurable needs to one of the values listed in ConfPIDnSigmaTPCMax!

  Configurable<std::vector<float>> ConfPIDnSigmaTPCMax{"PIDnSigmaTPCMax",
                                                       std::vector<float>{3.5f, 3.f, 2.5f},
                                                       "Vector of all possible nSigma values for Acceptance and Rejection (this needs to be in sync with FemtoDreamProducerTask.ConfTrkPIDnSigmaMax)"};
  // this configurable needs to be in sync with FemtoDreamProducerTask.ConfTrkPIDnSigmaMax
  // do not use more than 3 values, otherwise the FemtoDreamProducerTask will break!

  Configurable<float> confPIDRejection{"PIDRejection", 3.5, "nSigma for rejection bogus Deuterons (set it to a negative value to disable the rejection)"};
  // the value used in this configurable needs to one of the values listed in ConfPIDnSigmaTPCMax
  // set it to a negative value to disable the rejection

  // suggestion for setting ConfTrkTPIDspecies of the FemtoDreamProducerTask
  // ConfTrkTPIDspecies = {0, <- Electron PID at Index 0
  //                       2, <- Pion PID at Index 1
  //                       4, <- Proton PID at Index 2
  //                       5} <- Deuteron PID at Index 3
  //                       total of 4 indices
  // will become clear in the following

  Configurable<int> confPIDProtonIndex{"PIDProtonIndex", 2, "Index of Proton PID in ConfTrkTPIDspecies of the FemtoDreamProducerTask"};
  Configurable<int> confPIDDeuteronIndex{"PIDDeuteronIndex", 3, "Index of Deuteron PID in ConfTrkTPIDspecies of the FemtoDreamProducerTask"};
  Configurable<int> confPIDIndexMax{"PIDIndexMax", 4, "Number of Indices in ConfTrkTPIDspecies of the FemtoDreamProducerTask"};

  Configurable<std::vector<int>> ConfPIDRejectionSpeciesIndex{"PIDRejectionSpeciesIndex", std::vector<int>{0, 1, 2}, "Indices of the particles we want to reject from the Deuteron signal"};
  // when configuring the FemtoDreamProducerTask, select for ConfTrkTPIDspecies at least Protons (=4) and Deuterons (=5)
  // for the rejection to work properly also select electron (=0) and pions (=2)
  // suppose the FemtoDreamProducerTask is configured as above
  // and we want to rejection electrons (at index 0), pions (at index 1) and proton (at index 2)
  // set this configurable to PIDRejectionSpeciesIndex = {1,2,3}

  Configurable<bool> ConfClosePairRejection{"ClosePairRejection", true, "Reject close pairs or not"};
  Configurable<float>
    ldeltaPhiMax{"ldeltaPhiMax", 0.010, "Max limit of delta phi"};
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

  Preslice<o2::aod::FemtoDreamParticles> perCol = aod::femtodreamparticle::femtoDreamCollisionId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  // containers for close pair rejection
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> closePairRejectionTT;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> closePairRejectionTV0;

  bool SelectParticlePID(aod::femtodreamparticle::cutContainerType const& pidCut, int vSpecies, float momentum)
  {
    bool pidSelection = false;
    bool rejectDeuteron = false;
    if (vSpecies == o2::track::PID::Proton) {
      // use momentum dependend (TPC or TPC&TOF) pid selection for protons
      pidSelection = isFullPIDSelected(pidCut,
                                       momentum,
                                       confPIDThreshold.value,
                                       std::vector<int>{confPIDProtonIndex.value},
                                       confPIDIndexMax.value,
                                       ConfPIDnSigmaTPCMax.value,
                                       confPIDnSigmaTPCAcceptance.value,
                                       confPIDnSigmTPCTOFAcceptance.value);
    } else if (vSpecies == o2::track::PID::Deuteron) {
      // use additional rejection for deuterons
      if (confPIDRejection.value > 0.) {
        for (auto rejectSpecies : ConfPIDRejectionSpeciesIndex.value) {
          // if the PID is selected for another particle species, reject this candidate
          rejectDeuteron = isPIDSelected(pidCut, std::vector<int>{rejectSpecies},
                                         confPIDIndexMax.value,
                                         confPIDRejection.value,
                                         ConfPIDnSigmaTPCMax.value,
                                         kDetector::kTPC);
          // if a deuteron candidate is found which could also be another particle we want to reject, break out of the loop
          if (rejectDeuteron) {
            break;
          }
        }
      }
      // and reject the candidate, otherwise check if it is fullfills the deuteron hypothesis
      if (!rejectDeuteron) {
        pidSelection = isPIDSelected(pidCut,
                                     std::vector<int>{confPIDDeuteronIndex.value},
                                     confPIDIndexMax.value,
                                     confPIDnSigmaTPCAcceptance.value,
                                     ConfPIDnSigmaTPCMax.value,
                                     kDetector::kTPC);
      }
    } else {
      LOG(fatal) << "Other PID selections are not supported by this trigger" << std::endl;
    }
    return pidSelection;
  }

  void
    init(o2::framework::InitContext&)
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

    registry.add("fPtBeforeSel", "Transverse momentum of positive tracks", HistType::kTH1F, {{6000, 0, 6}});
    registry.add("fEtaBeforeSel", "Pseudorapidity of positive tracks", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiBeforeSel", "Azimuthal angle of positive tracks", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

    registry.add("fPtAntiBeforeSel", "Transverse momentum of negative tracks", HistType::kTH1F, {{6000, 0, 6}});
    registry.add("fEtaAntiBeforeSel", "Pseudorapidity of negative tracks", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiAntiBeforeSel", "Azimuthal angle of negative tracks", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

    bool plotPerRadii = true;

    if (KstarTrigger.value == 0 || KstarTrigger.value == 10) {
      registry.add("fKstarPD", "CF - same event pd distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiPD", "CF - same event pd distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{6000, 0, 6}});
      registry.add("fEtaProtonAfterSel", "Pseudorapidity of Protons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
      registry.add("fPhiProtonAfterSel", "Azimuthal angle of Protons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

      registry.add("fPtAntiProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{6000, 0, 6}});
      registry.add("fEtaAntiProtonAfterSel", "Pseudorapidity of AntiProtons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
      registry.add("fPhiAntiProtonAfterSel", "Azimuthal angle of AntiProtons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

      closePairRejectionTT.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 10) {
      registry.add("fKstarLD", "CF - same event ld distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiLD", "CF - same event ld distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtLambdaAfterSel", "Transverse momentum of Lambdas which passed selections", HistType::kTH1F, {{6000, 0, 6}});
      registry.add("fPtAntiLambdaAfterSel", "Transverse momentum of AntidLambdas which passed selections", HistType::kTH1F, {{6000, 0, 6}});

      registry.add("fMinvLambda", "Invariant mass of lambdas ", HistType::kTH1F, {{1000, 0.7, 1.5}});
      registry.add("fMinvAntiLambda", "Invariant mass of antilambdas ", HistType::kTH1F, {{1000, 0.7, 1.5}});

      closePairRejectionTV0.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    }

    registry.add("fPtDeuteronAfterSel", "Transverse momentum of Deuterons which passed selections", HistType::kTH1F, {{6000, 0, 6}});
    registry.add("fEtaDeuteronAfterSel", "Pseudorapidity of Deuterons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiDeuteronAfterSel", "Azimuthal angle of Deuterons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

    registry.add("fPtAntiDeuteronAfterSel", "Transverse momentum of Antideuterons which passed selections", HistType::kTH1F, {{6000, 0, 6}});
    registry.add("fEtaAntiDeuteronAfterSel", "Pseudorapidity of AntiDeuterons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiAntiDeuteronAfterSel", "Azimuthal angle of AntiDeuterons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});
  }

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassDeuteron = o2::constants::physics::MassDeuteron;
  float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(o2::aod::FemtoDreamCollision& col, o2::aod::FemtoDreamParticles& partsFemto)
  {
    // get partitions of all paritcles and antiparticles
    auto partsPD = partPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto partsAntiPD = partAntiPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);

    // get partions of V0s
    auto partsL = partL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto partsAntiL = partAntiL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);

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
      registry.fill(HIST("fEtaBeforeSel"), pd.eta());
      registry.fill(HIST("fPhiBeforeSel"), pd.phi());

      // select deuterons
      if (SelectParticlePID(pd.pidcut(), o2::track::PID::Deuteron, pd.p()) &&
          pd.pt() < confDeuteronPtMax.value &&
          pd.pt() > confDeuteronPtMin.value) {
        registry.fill(HIST("fPtDeuteronAfterSel"), pd.pt());
        registry.fill(HIST("fEtaDeuteronAfterSel"), pd.eta());
        registry.fill(HIST("fPhiDeuteronAfterSel"), pd.phi());
        Ndeuteron++;
      }

      // select protons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 10) {
        if (SelectParticlePID(pd.pidcut(), o2::track::PID::Proton, pd.p()) &&
            pd.pt() < confProtonPtMax.value &&
            pd.pt() > confProtonPtMin.value) {
          registry.fill(HIST("fPtProtonAfterSel"), pd.pt());
          registry.fill(HIST("fEtaProtonAfterSel"), pd.eta());
          registry.fill(HIST("fPhiProtonAfterSel"), pd.phi());
          Nproton++;
        }
      }
    }

    for (auto antipd : partsAntiPD) {
      registry.fill(HIST("fPtAntiBeforeSel"), antipd.pt());
      registry.fill(HIST("fEtaAntiBeforeSel"), antipd.eta());
      registry.fill(HIST("fPhiAntiBeforeSel"), antipd.phi());

      // select antideuterons
      if (SelectParticlePID(antipd.pidcut(), o2::track::PID::Deuteron, antipd.p()) &&
          antipd.pt() < confDeuteronPtMax.value &&
          antipd.pt() > confDeuteronPtMin.value) {
        registry.fill(HIST("fPtAntiDeuteronAfterSel"), antipd.pt());
        registry.fill(HIST("fEtaAntiDeuteronAfterSel"), antipd.eta());
        registry.fill(HIST("fPhiAntiDeuteronAfterSel"), antipd.phi());
        Nantideuteron++;
      }

      // select antiprotons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 10) {
        if (SelectParticlePID(antipd.pidcut(), o2::track::PID::Proton, antipd.p()) &&
            antipd.pt() < confProtonPtMax.value &&
            antipd.pt() > confProtonPtMin.value) {
          registry.fill(HIST("fPtAntiProtonAfterSel"), antipd.pt());
          registry.fill(HIST("fEtaAntiProtonAfterSel"), antipd.eta());
          registry.fill(HIST("fPhiAntiProtonAfterSel"), antipd.phi());
          Nantiproton++;
        }
      }
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 10) {
      // select lambdas
      for (auto lambda : partsL) {
        registry.fill(HIST("fPtLambdaAfterSel"), lambda.pt());
        registry.fill(HIST("fMinvLambda"), lambda.mLambda());
        Nlambda++;
      }
      for (auto antilambda : partsAntiL) {
        // select antilambdas
        registry.fill(HIST("fPtAntiLambdaAfterSel"), antilambda.pt());
        registry.fill(HIST("fMinvAntiLambda"), antilambda.mAntiLambda());
        Nantilambda++;
      }
    }

    bool keepEvent[kLAST_CFTwoBodyTriggers] = {false, false};
    int lowKstarPairs[kLAST_CFTwoBodyTriggers] = {0, 0};

    bool pdPair = false;
    bool dpPair = false;
    double kStar = 0.;

    // trigger for pd pairs
    if (KstarTrigger.value == 0 || KstarTrigger.value == 10) {
      if (Ndeuteron > 0 && Nproton > 0) {
        // loop over all unique combinations of particles, excluding the self combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsPD, partsPD))) {

          // check if it is a pd pair
          // p1 => proton
          // p2 => deuteron
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Proton, p1.p()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Deuteron, p2.p()) &&
              p1.pt() < confProtonPtMax.value &&
              p1.pt() > confProtonPtMin.value &&
              p2.pt() < confDeuteronPtMax.value &&
              p2.pt() > confDeuteronPtMin.value) {
            pdPair = true;
          } else {
            pdPair = false;
          }

          // check if it is dp pair
          // p1 => deuteron
          // p2 => proton
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Proton, p2.p()) &&
              p1.pt() < confDeuteronPtMax.value &&
              p1.pt() > confDeuteronPtMin.value &&
              p2.pt() < confProtonPtMax.value &&
              p2.pt() > confProtonPtMin.value) {
            dpPair = true;
          } else {
            dpPair = false;
          }

          // if neither is the case, skip
          if (!(pdPair || dpPair)) {
            continue;
          }

          // reject close pairs
          if (ConfClosePairRejection.value && closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
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
          // p1 => antiproton
          // p2 => antideuteron
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Proton, p1.p()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Deuteron, p2.p()) &&
              p1.pt() < confProtonPtMax.value &&
              p1.pt() > confProtonPtMin.value &&
              p2.pt() < confDeuteronPtMax.value &&
              p2.pt() > confDeuteronPtMin.value) {
            pdPair = true;
          } else {
            pdPair = false;
          }

          // check if it is (anti)dp pair
          // p1 => antideuteron
          // p2 => antiproton
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Proton, p2.p()) &&
              p1.pt() < confDeuteronPtMax.value &&
              p1.pt() > confDeuteronPtMin.value &&
              p2.pt() < confProtonPtMax.value &&
              p2.pt() > confProtonPtMin.value) {
            dpPair = true;
          } else {
            dpPair = false;
          }

          // if neither is the case, skip
          if (!(pdPair || dpPair)) {
            continue;
          }

          // reject close pairs
          if (ConfClosePairRejection.value && closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
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
    if (KstarTrigger.value == 1 || KstarTrigger.value == 10) {
      if (Ndeuteron > 0 && Nlambda > 0) {
        // loop over all unique combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsUpperIndexPolicy(partsPD, partsL))) {
          // check if the particle is a deuteron
          // we do not need to check the V0s
          if (!SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p())) {
            continue;
          }
          if (ConfClosePairRejection.value && closePairRejectionTV0.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          kStar = FemtoDreamMath::getkstar(p1, mMassDeuteron, p2, mMassLambda);
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
          if (!SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p())) {
            continue;
          }
          if (ConfClosePairRejection.value && closePairRejectionTV0.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassDeuteron, p2, mMassLambda);
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
      registry.fill(HIST("fProcessedEvents"), 2 + kPD); // first two bins are counter of all and rejected events
    }
    if (lowKstarPairs[kLD] > 0) {
      keepEvent[kLD] = true;
      registry.fill(HIST("fProcessedEvents"), 2 + kLD); // first two bins are counter of all and rejected events
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
