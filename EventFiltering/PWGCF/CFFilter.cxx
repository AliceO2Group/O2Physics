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
/// \author Laura Serksnyte, TU München, laura.serksnyte@cern.ch

#include <string>

#include "../filterTables.h"
#include "../../PWGCF/FemtoDream/FemtoUtils.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/FemtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/FemtoDreamMath.h"
#include "PWGCF/FemtoDream/FemtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/FemtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/FemtoDreamContainer.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

namespace
{

enum CFTriggers {
  kPPP = 0,
  kPPL,
  kPLL,
  kLLL,
  kNTriggers
};

static const std::vector<std::string> CfTriggerNames{"ppp", "ppL", "pLL", "LLL"};
// uint8_t trackTypeSel = o2::aod::femtodreamparticle::ParticleType::kTrack; Fix this to work instead of below hardcoded lines
// uint V0TypeSel = o2::aod::femtodreamparticle::ParticleType::kV0; Fix this to work instead of below hardcoded lines
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

struct CFFilter {
  SliceCache cache;

  Produces<aod::CFFilters> tags;

  Configurable<std::vector<float>> confQ3TriggerLimit{"Q3TriggerLimitC", std::vector<float>{0.6f, 0.6f, 0.6f, 0.6f}, "Q3 limit for selection"};
  Configurable<int> Q3Trigger{"Q3Trigger", 0, "Choice which trigger to run"};
  Configurable<bool> performCPR{"performCPR", true, "Perform or not the close pair rejection"};
  Configurable<float> ldeltaPhiMax{"ldeltaPhiMax", 0.010, "Max limit of delta phi"};
  Configurable<float> ldeltaEtaMax{"ldeltaEtaMax", 0.010, "Max limit of delta eta"};
  Configurable<std::vector<float>> ConfPIDnSigmaMax{"PIDnSigmaMax",
                                                    std::vector<float>{4.f, 3.f},
                                                    "Vector of all possible nSigma values for Acceptance and Rejection (this needs to be in sync with FemtoDreamProducerTask.ConfTrkPIDnSigmaMax)"};
  Configurable<int> ConfPIDProtonIndex{"PIDProtonIndex", 2, "Index of Proton PID in ConfTrkTPIDspecies of the FemtoDreamProducerTask"};
  Configurable<int> ConfPIDIndexMax{"PIDIndexMax", 4, "Number of Indices in ConfTrkTPIDspecies of the FemtoDreamProducerTask"};
  Configurable<float> ConfPIDThreshold{"PThreshold", 0.75f, "P threshold for TPC/TPC&TOF selection"};
  Configurable<float> ConfPIDnSigma{"PIDnSigma", 4.f, "nSigma value for Proton PID"};

  // Obtain particle and antiparticle candidates of protons and lambda hyperons for current femto collision
  Partition<o2::aod::FemtoDreamParticles> partsProton0Part = (o2::aod::femtodreamparticle::partType == Track) && ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0); // Consider later: && ((o2::aod::femtodreamparticle::pidcut & knSigmaProton) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partsLambda0Part = (o2::aod::femtodreamparticle::partType == V0) && ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partsProton1Part = (o2::aod::femtodreamparticle::partType == Track) && ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0); // Consider later: && ((o2::aod::femtodreamparticle::pidcut & knSigmaProton) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partsLambda1Part = (o2::aod::femtodreamparticle::partType == V0) && ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);

  Preslice<o2::aod::FemtoDreamParticles> perCol = aod::femtodreamparticle::femtoDreamCollisionId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  // FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleanerTT; Currently not used, will be needed later
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> pairCleanerTV;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> closePairRejectionTT;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> closePairRejectionTV0;

  void init(o2::framework::InitContext&)
  {
    bool plotPerRadii = true;

    closePairRejectionTT.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    closePairRejectionTV0.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    registry.add("fProcessedEvents", "CF - event filtered;;events", HistType::kTH1F, {{6, -0.5, 5.5}});
    // dont use hardcoded 6
    std::array<std::string, kNTriggers + 2> eventTitles = {"all", "rejected", "p-p-p", "p-p-L", "p-L-L", "L-L-L"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }
    registry.add("fMultiplicityBefore", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fMultiplicityAfter", "Multiplicity of events which passed ppp trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fZvtxBefore", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("fZvtxAfter", "Zvtx of events which passed ppp trigger", HistType::kTH1F, {{1000, -15, 15}});

    registry.add("fPtBeforePPP", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtAfterPPP", "Transverse momentum  of processed tracks which passed  selections", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtBeforeAntiPPP", "Transverse momentum of all processed antitracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fPtAfterAntiPPP", "Transverse momentum  of processed antitracks passed selection", HistType::kTH1F, {{1000, 0, 10}});

    if (Q3Trigger == 0 || Q3Trigger == 11 || Q3Trigger == 1111) {
      registry.add("fSameEventPartPPP", "CF - same event ppp distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fSameEventAntiPartPPP", "CF - same event ppp distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});
    }
    if (Q3Trigger == 1 || Q3Trigger == 11 || Q3Trigger == 1111) {
      registry.add("fSameEventPartPPL", "CF - same event ppL distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fSameEventAntiPartPPL", "CF - same event ppL distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtPPL", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fPtAntiPPL", "Transverse momentum of all processed antitracks", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fMinvLambda", "Invariant mass of lambdas ", HistType::kTH1F, {{1000, 0.7, 1.5}});
      registry.add("fMinvAntiLambda", "Invariant mass of antilambdas ", HistType::kTH1F, {{1000, 0.7, 1.5}});
      registry.add("fMinvLambdaBefore", "Invariant mass of lambdas before ", HistType::kTH1F, {{1000, 0.7, 1.5}});
      registry.add("fMinvAntiLambdaBefore", "Invariant mass of antilambdas before", HistType::kTH1F, {{1000, 0.7, 1.5}});
    }
    if (Q3Trigger == 2 || Q3Trigger == 1111) {
      registry.add("fSameEventPartPLL", "CF - same event pLL distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fSameEventAntiPartPLL", "CF - same event pLL distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});
    }
    if (Q3Trigger == 3 || Q3Trigger == 1111) {
      registry.add("fSameEventPartLLL", "CF - same event LLL distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fSameEventAntiPartLLL", "CF - same event LLL distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});
    }
  }

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(o2::aod::FemtoDreamCollision& col, o2::aod::FemtoDreamParticles& partsFemto)
  {
    auto partsProton0 = partsProton0Part->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto partsLambda0 = partsLambda0Part->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto partsProton1 = partsProton1Part->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto partsLambda1 = partsLambda1Part->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);

    auto magneticField = col.magField();

    registry.get<TH1>(HIST("fProcessedEvents"))->Fill(0);
    registry.get<TH1>(HIST("fMultiplicityBefore"))->Fill(col.multV0M());
    registry.get<TH1>(HIST("fZvtxBefore"))->Fill(col.posZ());

    int prot = 0;
    int antiprot = 0;
    for (auto p1pt : partsProton0) {
      registry.get<TH1>(HIST("fPtBeforePPP"))->Fill(p1pt.pt());
      if (isFullPIDSelected(p1pt.pidcut(),
                            p1pt.p(),
                            ConfPIDThreshold.value,
                            std::vector<int>{ConfPIDProtonIndex.value},
                            ConfPIDIndexMax.value,
                            ConfPIDnSigmaMax.value,
                            ConfPIDnSigma.value,
                            ConfPIDnSigma.value)) {
        registry.get<TH1>(HIST("fPtAfterPPP"))->Fill(p1pt.pt());
        prot++;
      }
    }
    for (auto p1pt : partsProton1) {
      registry.get<TH1>(HIST("fPtBeforeAntiPPP"))->Fill(p1pt.pt());
      if (isFullPIDSelected(p1pt.pidcut(),
                            p1pt.p(),
                            ConfPIDThreshold.value,
                            std::vector<int>{ConfPIDProtonIndex.value},
                            ConfPIDIndexMax.value,
                            ConfPIDnSigmaMax.value,
                            ConfPIDnSigma.value,
                            ConfPIDnSigma.value)) {
        registry.get<TH1>(HIST("fPtAfterAntiPPP"))->Fill(p1pt.pt());
        antiprot++;
      }
    }

    for (auto lambda : partsLambda0) {
      registry.get<TH1>(HIST("fMinvLambdaBefore"))->Fill(lambda.mLambda());
    }
    for (auto antilambda : partsLambda1) {
      registry.get<TH1>(HIST("fMinvAntiLambdaBefore"))->Fill(antilambda.mAntiLambda());
    }

    bool keepEvent[kNTriggers] = {false, false, false, false};
    int lowQ3Triplets[kNTriggers] = {0, 0, 0, 0};

    if (partsFemto.size() != 0) {
      registry.get<TH1>(HIST("fMultiplicityAfter"))->Fill(col.multV0M());
      registry.get<TH1>(HIST("fZvtxAfter"))->Fill(col.posZ());
      auto Q3TriggerLimit = (std::vector<float>)confQ3TriggerLimit;
      // __________________________________________________________________________________________________________
      // TRIGGER FOR PPP TRIPLETS
      if (Q3Trigger == 0 || Q3Trigger == 1111 || Q3Trigger == 11) {
        if (prot >= 3) {
          // test default combinations options
          for (auto& [p1, p2, p3] : combinations(partsProton0, partsProton0, partsProton0)) {
            if (!isFullPIDSelected(p1.pidcut(),
                                   p1.p(),
                                   ConfPIDThreshold.value,
                                   std::vector<int>{ConfPIDProtonIndex.value},
                                   ConfPIDIndexMax.value,
                                   ConfPIDnSigmaMax.value,
                                   ConfPIDnSigma.value,
                                   ConfPIDnSigma.value) ||
                !isFullPIDSelected(p2.pidcut(),
                                   p2.p(),
                                   ConfPIDThreshold.value,
                                   std::vector<int>{ConfPIDProtonIndex.value},
                                   ConfPIDIndexMax.value,
                                   ConfPIDnSigmaMax.value,
                                   ConfPIDnSigma.value,
                                   ConfPIDnSigma.value) ||
                !isFullPIDSelected(p3.pidcut(),
                                   p3.p(),
                                   ConfPIDThreshold.value,
                                   std::vector<int>{ConfPIDProtonIndex.value},
                                   ConfPIDIndexMax.value,
                                   ConfPIDnSigmaMax.value,
                                   ConfPIDnSigma.value,
                                   ConfPIDnSigma.value)) {
              continue;
            }
            // Think if pair cleaning is needed in current framework
            // Run close pair rejection
            if (performCPR) {
              if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
                continue;
              }
              if (closePairRejectionTT.isClosePair(p1, p3, partsFemto, magneticField)) {
                continue;
              }
              if (closePairRejectionTT.isClosePair(p2, p3, partsFemto, magneticField)) {
                continue;
              }
            }
            auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, p2, mMassProton, p3, mMassProton);
            registry.get<TH1>(HIST("fSameEventPartPPP"))->Fill(Q3);
            if (Q3 < Q3TriggerLimit.at(0)) {
              lowQ3Triplets[0]++;
            }
          }
        } // end if

        if (lowQ3Triplets[0] == 0) { // Use this in final version only, for testing comment { // if at least one triplet found in particles, no need to check antiparticles
          if (antiprot >= 3) {
            for (auto& [p1, p2, p3] : combinations(partsProton1, partsProton1, partsProton1)) {

              if (!isFullPIDSelected(p1.pidcut(),
                                     p1.p(),
                                     ConfPIDThreshold.value,
                                     std::vector<int>{ConfPIDProtonIndex.value},
                                     ConfPIDIndexMax.value,
                                     ConfPIDnSigmaMax.value,
                                     ConfPIDnSigma.value,
                                     ConfPIDnSigma.value) ||
                  !isFullPIDSelected(p2.pidcut(),
                                     p2.p(),
                                     ConfPIDThreshold.value,
                                     std::vector<int>{ConfPIDProtonIndex.value},
                                     ConfPIDIndexMax.value,
                                     ConfPIDnSigmaMax.value,
                                     ConfPIDnSigma.value,
                                     ConfPIDnSigma.value) ||
                  !isFullPIDSelected(p3.pidcut(),
                                     p3.p(),
                                     ConfPIDThreshold.value,
                                     std::vector<int>{ConfPIDProtonIndex.value},
                                     ConfPIDIndexMax.value,
                                     ConfPIDnSigmaMax.value,
                                     ConfPIDnSigma.value,
                                     ConfPIDnSigma.value)) {
                continue;
              }
              // Think if pair cleaning is needed in current framework
              // Run close pair rejection
              if (performCPR) {
                if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
                  continue;
                }
                if (closePairRejectionTT.isClosePair(p1, p3, partsFemto, magneticField)) {
                  continue;
                }
                if (closePairRejectionTT.isClosePair(p2, p3, partsFemto, magneticField)) {
                  continue;
                }
              }
              auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, p2, mMassProton, p3, mMassProton);
              registry.get<TH1>(HIST("fSameEventAntiPartPPP"))->Fill(Q3);
              if (Q3 < Q3TriggerLimit.at(0)) {
                lowQ3Triplets[0]++;
              }
            }
          }
        } // end if
        //}
      }
      // __________________________________________________________________________________________________________
      // TRIGGER FOR PPL TRIPLETS
      if (Q3Trigger == 1 || Q3Trigger == 1111 || Q3Trigger == 11) {
        if (partsLambda0.size() >= 1 && prot >= 2) {
          for (auto& partLambda : partsLambda0) {
            registry.get<TH1>(HIST("fPtPPL"))->Fill(partLambda.pt());
            registry.get<TH1>(HIST("fMinvLambda"))->Fill(partLambda.mLambda());
            if (!pairCleanerTV.isCleanPair(partLambda, partLambda, partsFemto)) {
              continue;
            }
            for (auto& [p1, p2] : combinations(partsProton0, partsProton0)) {
              if (!isFullPIDSelected(p1.pidcut(),
                                     p1.p(),
                                     ConfPIDThreshold.value,
                                     std::vector<int>{ConfPIDProtonIndex.value},
                                     ConfPIDIndexMax.value,
                                     ConfPIDnSigmaMax.value,
                                     ConfPIDnSigma.value,
                                     ConfPIDnSigma.value) ||
                  !isFullPIDSelected(p2.pidcut(),
                                     p2.p(),
                                     ConfPIDThreshold.value,
                                     std::vector<int>{ConfPIDProtonIndex.value},
                                     ConfPIDIndexMax.value,
                                     ConfPIDnSigmaMax.value,
                                     ConfPIDnSigma.value,
                                     ConfPIDnSigma.value)) {
                continue;
              }
              if (performCPR) {
                if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
                  continue;
                }
                if (closePairRejectionTV0.isClosePair(p1, partLambda, partsFemto, magneticField)) {
                  continue;
                }
                if (closePairRejectionTV0.isClosePair(p2, partLambda, partsFemto, magneticField)) {
                  continue;
                }
              }
              auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, p2, mMassProton, partLambda, mMassLambda);
              registry.get<TH1>(HIST("fSameEventPartPPL"))->Fill(Q3);
              if (Q3 < Q3TriggerLimit.at(1)) {
                lowQ3Triplets[1]++;
              }
            }
          }
        } // end if

        if (lowQ3Triplets[1] == 0) { // if at least one triplet found in particles, no need to check antiparticles
          if (partsLambda1.size() >= 1 && antiprot >= 2) {
            for (auto& partLambda : partsLambda1) {
              registry.get<TH1>(HIST("fPtAntiPPL"))->Fill(partLambda.pt());
              registry.get<TH1>(HIST("fMinvAntiLambda"))->Fill(partLambda.mAntiLambda());
              if (!pairCleanerTV.isCleanPair(partLambda, partLambda, partsFemto)) {
                continue;
              }
              for (auto& [p1, p2] : combinations(partsProton1, partsProton1)) {
                if (!isFullPIDSelected(p1.pidcut(),
                                       p1.p(),
                                       ConfPIDThreshold.value,
                                       std::vector<int>{ConfPIDProtonIndex.value},
                                       ConfPIDIndexMax.value,
                                       ConfPIDnSigmaMax.value,
                                       ConfPIDnSigma.value,
                                       ConfPIDnSigma.value) ||
                    !isFullPIDSelected(p2.pidcut(),
                                       p2.p(),
                                       ConfPIDThreshold.value,
                                       std::vector<int>{ConfPIDProtonIndex.value},
                                       ConfPIDIndexMax.value,
                                       ConfPIDnSigmaMax.value,
                                       ConfPIDnSigma.value,
                                       ConfPIDnSigma.value)) {
                  continue;
                }
                if (performCPR) {
                  if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
                    continue;
                  }
                  if (closePairRejectionTV0.isClosePair(p1, partLambda, partsFemto, magneticField)) {
                    continue;
                  }
                  if (closePairRejectionTV0.isClosePair(p2, partLambda, partsFemto, magneticField)) {
                    continue;
                  }
                }
                auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, p2, mMassProton, partLambda, mMassLambda);
                registry.get<TH1>(HIST("fSameEventAntiPartPPL"))->Fill(Q3);
                if (Q3 < Q3TriggerLimit.at(1)) {
                  lowQ3Triplets[1]++;
                }
              }
            }
          } // end if
        }
      }

      // __________________________________________________________________________________________________________
      // TRIGGER FOR PLL TRIPLETS
      if (Q3Trigger == 2 || Q3Trigger == 1111) {
        if (partsLambda0.size() >= 2 && prot >= 1) {
          for (auto& p1 : partsProton0) {
            if (!isFullPIDSelected(p1.pidcut(),
                                   p1.p(),
                                   ConfPIDThreshold.value,
                                   std::vector<int>{ConfPIDProtonIndex.value},
                                   ConfPIDIndexMax.value,
                                   ConfPIDnSigmaMax.value,
                                   ConfPIDnSigma.value,
                                   ConfPIDnSigma.value)) {
              continue;
            }
            for (auto& [partLambda1, partLambda2] : combinations(partsLambda0, partsLambda0)) {
              // maybe implement L1-L2 no shared tracks
              if (!pairCleanerTV.isCleanPair(partLambda1, partLambda1, partsFemto)) {
                continue;
              }
              if (!pairCleanerTV.isCleanPair(partLambda2, partLambda2, partsFemto)) {
                continue;
              }

              if (performCPR) {
                if (closePairRejectionTV0.isClosePair(p1, partLambda1, partsFemto, magneticField)) {
                  continue;
                }
                if (closePairRejectionTV0.isClosePair(p1, partLambda2, partsFemto, magneticField)) {
                  continue;
                }
                // maybe implement L-L cpr
              }
              auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, partLambda1, mMassLambda, partLambda2, mMassLambda);
              registry.get<TH1>(HIST("fSameEventPartPLL"))->Fill(Q3);
              if (Q3 < Q3TriggerLimit.at(2)) {
                lowQ3Triplets[2]++;
              }
            }
          }
        } // end if

        if (lowQ3Triplets[2] == 0) { // if at least one triplet found in particles, no need to check antiparticles
          if (partsLambda1.size() >= 2 && antiprot >= 1) {
            for (auto& p1 : partsProton1) {
              if (!isFullPIDSelected(p1.pidcut(),
                                     p1.p(),
                                     ConfPIDThreshold.value,
                                     std::vector<int>{ConfPIDProtonIndex.value},
                                     ConfPIDIndexMax.value,
                                     ConfPIDnSigmaMax.value,
                                     ConfPIDnSigma.value,
                                     ConfPIDnSigma.value)) {
                continue;
              }
              for (auto& [partLambda1, partLambda2] : combinations(partsLambda1, partsLambda1)) {
                // maybe implement L1-L2 no shared tracks
                if (!pairCleanerTV.isCleanPair(partLambda1, partLambda1, partsFemto)) {
                  continue;
                }
                if (!pairCleanerTV.isCleanPair(partLambda2, partLambda2, partsFemto)) {
                  continue;
                }

                if (performCPR) {
                  if (closePairRejectionTV0.isClosePair(p1, partLambda1, partsFemto, magneticField)) {
                    continue;
                  }
                  if (closePairRejectionTV0.isClosePair(p1, partLambda2, partsFemto, magneticField)) {
                    continue;
                  }
                  // maybe implement L-L cpr
                }
                auto Q3 = FemtoDreamMath::getQ3(p1, mMassProton, partLambda1, mMassLambda, partLambda2, mMassLambda);
                registry.get<TH1>(HIST("fSameEventAntiPartPLL"))->Fill(Q3);
                if (Q3 < Q3TriggerLimit.at(2)) {
                  lowQ3Triplets[2]++;
                }
              }
            }
          } // end if
        }
      }

      // __________________________________________________________________________________________________________
      // TRIGGER FOR LLL TRIPLETS
      if (Q3Trigger == 3 || Q3Trigger == 1111) {
        if (partsLambda0.size() >= 3) {
          // test default combinations options
          for (auto& [partLambda1, partLambda2, partLambda3] : combinations(partsLambda0, partsLambda0, partsLambda0)) {
            if (!pairCleanerTV.isCleanPair(partLambda1, partLambda1, partsFemto)) {
              continue;
            }
            if (!pairCleanerTV.isCleanPair(partLambda2, partLambda2, partsFemto)) {
              continue;
            }
            if (!pairCleanerTV.isCleanPair(partLambda3, partLambda3, partsFemto)) {
              continue;
            }
            // Run close pair rejection
            if (performCPR) {
              // check close pair rejection for L-L
            }
            auto Q3 = FemtoDreamMath::getQ3(partLambda1, mMassLambda, partLambda2, mMassLambda, partLambda3, mMassLambda);
            registry.get<TH1>(HIST("fSameEventPartLLL"))->Fill(Q3);
            if (Q3 < Q3TriggerLimit.at(3)) {
              lowQ3Triplets[3]++;
            }
          }
        } // end if

        if (lowQ3Triplets[3] == 0) { // Use this in final version only, for testing comment { // if at least one triplet found in particles, no need to check antiparticles

          if (partsLambda1.size() >= 3) {
            // test default combinations options
            for (auto& [partLambda1, partLambda2, partLambda3] : combinations(partsLambda1, partsLambda1, partsLambda1)) {
              if (!pairCleanerTV.isCleanPair(partLambda1, partLambda1, partsFemto)) {
                continue;
              }
              if (!pairCleanerTV.isCleanPair(partLambda2, partLambda2, partsFemto)) {
                continue;
              }
              if (!pairCleanerTV.isCleanPair(partLambda3, partLambda3, partsFemto)) {
                continue;
              }
              // Run close pair rejection
              if (performCPR) {
                // check close pair rejection for L-L
              }
              auto Q3 = FemtoDreamMath::getQ3(partLambda1, mMassLambda, partLambda2, mMassLambda, partLambda3, mMassLambda);
              registry.get<TH1>(HIST("fSameEventAntiPartLLL"))->Fill(Q3);
              if (Q3 < Q3TriggerLimit.at(3)) {
                lowQ3Triplets[3]++;
              }
            }
          }

        } // end if
        //}
      }
    }

    if (lowQ3Triplets[0] > 0) {
      keepEvent[kPPP] = true;
    }

    if (lowQ3Triplets[1] > 0) {
      keepEvent[kPPL] = true;
    }

    if (lowQ3Triplets[2] > 0) {
      keepEvent[kPLL] = true;
    }

    if (lowQ3Triplets[3] > 0) {
      keepEvent[kLLL] = true;
    }

    tags(keepEvent[kPPP], keepEvent[kPPL], keepEvent[kPLL], keepEvent[kLLL]);

    if (!keepEvent[kPPP] && !keepEvent[kPPL] && !keepEvent[kPLL] && !keepEvent[kLLL]) {
      registry.get<TH1>(HIST("fProcessedEvents"))->Fill(1);
    } else {
      for (int iTrigger{0}; iTrigger < kNTriggers; iTrigger++) {
        if (keepEvent[iTrigger]) {
          registry.get<TH1>(HIST("fProcessedEvents"))->Fill(iTrigger + 2);
        }
      }
    } // end else
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilter>(cfg)};
}
