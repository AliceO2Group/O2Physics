// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamProducerTaskForSpecificAnalysis.cxx
/// \brief Tasks that reads the track tables and creates track triplets; only three identical particles can be used
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainerThreeBody.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoDreamProducerTaskForSpecificAnalysis {

  SliceCache cache;

  Produces<aod::StoredFDCollisions> outputCollision;
  Produces<aod::StoredFDParticles> outputParts;

  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  float mMassOne = -999, mMassTwo = -999, mMassThree = -999;
  int collisions = 0;

  // Require bitmask selection for candidates
  Configurable<bool> confRequireBitmask{"confRequireBitmask", false, "Require bitmask selection for candidates"};

  // Number of candidates required
  Configurable<int> confNumberOfTracks{"confNumberOfTracks", 0, "Number of tracks"};
  Configurable<int> confNumberOfV0{"confNumberOfV0", 1, "Number of V0"};
  Configurable<int> confNumberOfCascades{"confNumberOfCascades", 0, "Number of Cascades"};
   Configurable<int> confNumberOfReso{"confNumberOfV0", 1, "Number of Reso"};

  /// Track selection
  Configurable<float> confPIDthrMom{"confPIDthrMom", 1.f, "Momentum threshold from which TPC and TOF are required for PID"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confTPCPIDBit{"confTPCPIDBit", 16, "PID TPC bit from cutCulator "};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confTPCTOFPIDBit{"confTPCTOFPIDBit", 8, "PID TPCTOF bit from cutCulator"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confCutPart{"confCutPart", 0, "Track - Selection bit from cutCulator for part"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confCutPartAntiPart{"confCutPartAntiPart", 0, "Track - Selection bit from cutCulator for antipart"};

  /// Partition for selected particles
  Partition<aod::FDParticles> selectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= confPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, confTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, confTPCTOFPIDBit));

  /// V0 selection
  Configurable<float> confMinInvMassV0{"confMinInvMassV0", 1.08, "Minimum invariant mass of V0 (particle)"};
  Configurable<float> confMaxInvMassV0{"confMaxInvMassV0", 1.15, "Maximum invariant mass of V0 (particle)"};
  Configurable<float> confMinInvMassAntiV0{"confMinInvMassAntiV0", 1.08, "Minimum invariant mass of V0 (antiparticle)"};
  Configurable<float> confMaxInvMassAntiV0{"confMaxInvMassAntiV0", 1.15, "Maximum invariant mass of V0 (antiparticle)"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confCutV0SameForAntipart{"confCutV0SameForAntipart", 0, "V0 - Selection bit from cutCulator for part/antipart"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildPosCutV0{"confChildPosCutV0", 149, "Selection bit for positive child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildPosTPCBitV0{"confChildPosTPCBitV0", 2, "PID TPC bit for positive child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildNegCutV0{"confChildNegCutV0", 149, "Selection bit for negative child of V0"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> confChildNegTPCBitV0{"confChildNegTPCBitV0", 2, "PID TPC bit for negative child of V0"};

  /// Cascade selection
  Configurable<float> confMinInvMassCascade{"confMinInvMassCascade", 1.2, "Minimum invariant mass of Cascade (particle)"};
  Configurable<float> confMaxInvMassCascade{"confMaxInvMassCascade", 1.5, "Maximum invariant mass of Cascade (particle)"};


   struct : ConfigurableGroup { //set loosest cuts
    std::string prefix = std::string("Reso");
    Configurable<int> pdgCode{"pdgCode", 333, "PDG code of particle 2 Reso"};

    Configurable<float> confMinInvMassReso{"confMinInvMassReso", 0.86, "Minimum invariant mass of Reso (particle)"};
    Configurable<float> confMaxInvMassReso{"confMaxInvMassReso", 1.3, "Maximum invariant mass of Reso (particle)"};


    Configurable<femtodreamparticle::cutContainerType> daughPosCutBit{"daughPosCutBit", 2401446, "Selection bit for positive child of V02"}; // K+
    Configurable<femtodreamparticle::cutContainerType> daughPosTPCBit{"daughPosTPCBit", 4096, "PID TPC bit for positive child of V02"};
    Configurable<femtodreamparticle::cutContainerType> daughPosTPCTOFBit{"daughPosTPCTOFBit", 2048, "PID TOF bit for positive child of V02"};
    Configurable<femtodreamparticle::cutContainerType> daughNegCutBit{"daughNegCutBit", 2401445, "Selection bit for negative child of V02"}; // K-
    Configurable<femtodreamparticle::cutContainerType> daughNegMergedTPCBit{"daughNegMergedTPCBit", 16386, "PID TPC bit for negative child of V02"};
    Configurable<femtodreamparticle::cutContainerType> daughNegMergedTPCTOFBit{"daughNegMergedTPCTOFBit", 8194, "PID TOF bit for negative child of V02"};

    Configurable<aod::femtodreamparticle::partType>  partType1{"partType1", kResoPosdaughTPC_NegdaughTPC};
    Configurable<aod::femtodreamparticle::partType>  partType2{"partType2", kResoPosdaughTOF_NegdaughTOF};
    Configurable<aod::femtodreamparticle::partType>  partType3{"partType3", kResoPosdaughTOF_NegdaughTPC};
    Configurable<aod::femtodreamparticle::partType>  partType4{"partType4", kResoPosdaughTPC_NegdaughTOF};

  } Reso;

  // Partition for selected particles
  Partition<aod::FDParticles> selectedV0s = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0));
  Partition<aod::FDParticles> selectedCascades = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade));

  Partition<aod::FDParticles> selectedResos =  (ifnode(aod::femtodreamparticle::partType == uint8_t(Reso.partType1), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCBit), false) ||
                                                ifnode(aod::femtodreamparticle::partType == uint8_t(Reso.partType1), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCTOFBit), false) ||
                                                ifnode(aod::femtodreamparticle::partType == uint8_t(Reso.partType1), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCBit), false) ||
                                                ifnode(aod::femtodreamparticle::partType == uint8_t(Reso.partType1), ncheckbit(aod::femtodreamparticle::pidcut, Reso2.daughPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, Reso2.daughNegMergedTPCTOFBit), false)) &&
                                               (aod::femtodreamparticle::mLambda > Reso.confMinInvMassReso) &&
                                               (aod::femtodreamparticle::mLambda < Reso.confMaxInvMassReso);

  
  HistogramRegistry eventRegistry{"eventRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr uint32_t kSignPlusMask = 1 << 1;

  template <typename T>
  int getRowDaughters(int daughID, T const& vecID)
  {
    int rowInPrimaryTrackTableDaugh = -1;
    for (size_t i = 0; i < vecID.size(); i++) {
      if (vecID.at(i) == daughID) {
        rowInPrimaryTrackTableDaugh = i;
        break;
      }
    }
    return rowInPrimaryTrackTableDaugh;
  }

  void init(InitContext&)
  {
    eventRegistry.add("hStatistiscs", ";bin;Entries", kTH1F, {{3, 0, 3}});
    // Never run V0s and Cascades together as this will DOUBLE the track number and induce self correlations
    if ((doprocessCollisionsWithNTracksAndNCascades && doprocessCollisionsWithNTracksAndNV0)) {
      LOG(fatal) << "Never run V0s and Cascades together as this will DOUBLE the track number and induce self correlations!";
    }
  }

  /// This function stores accepted collisions in derived data
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param groupSelectedV0s partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table
  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedData(const o2::aod::FDCollision& col, PartitionType groupSelectedTracks, PartitionType groupSelectedV0s, PartType parts)
  {
    /// check tracks
    int tracksCount = 0;
    int antitracksCount = 0;
    for (const auto& part : groupSelectedTracks) {
      if (part.cut() & 1) {
        if (!confRequireBitmask || ncheckbit(part.cut(), confCutPartAntiPart)) {
          antitracksCount++;
        }
      } else {
        if (!confRequireBitmask || ncheckbit(part.cut(), confCutPart)) {
          tracksCount++;
        }
      }
    }

    /// check V0s
    int v0Count = 0;
    int antiV0Count = 0;
    for (const auto& V0 : groupSelectedV0s) {
      if ((V0.mLambda() > confMinInvMassV0) && (V0.mLambda() < confMaxInvMassV0)) {
        if (confRequireBitmask) {
          if (ncheckbit(V0.cut(), confCutV0SameForAntipart)) {
            const auto& posChild = parts.iteratorAt(V0.index() - 2);
            const auto& negChild = parts.iteratorAt(V0.index() - 1);
            if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
                 (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
                 (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
                 (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
              v0Count++;
            }
          }
        } else {
          v0Count++;
        }
      } else if ((V0.mAntiLambda() > confMinInvMassAntiV0) && (V0.mAntiLambda() < confMaxInvMassAntiV0)) {
        if (confRequireBitmask) {
          if (ncheckbit(V0.cut(), confCutV0SameForAntipart)) {
            const auto& posChild = parts.iteratorAt(V0.index() - 2);
            const auto& negChild = parts.iteratorAt(V0.index() - 1);
            if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
                 (posChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0 && // exchanged values because checking antiparticle daughters and pid of particles exchange
                 (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
                 (negChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0)) { // exchanged values because checking antiparticle daughters and pid of particles exchange
              antiV0Count++;
            }
          }
        } else {
          antiV0Count++;
        }
      }
    }

    std::vector<int> tmpIDtrack;

    if ((v0Count >= confNumberOfV0 && tracksCount >= confNumberOfTracks) || (antiV0Count >= confNumberOfV0 && antitracksCount >= confNumberOfTracks)) {
      eventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());
      for (const auto& femtoParticle : parts) {
        if (aod::femtodreamparticle::ParticleType::kTrack == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
          tmpIDtrack.push_back(femtoParticle.index());
        }
        if (aod::femtodreamparticle::ParticleType::kV0Child == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          const auto& children = femtoParticle.childrenIds();
          int childId = (children[0] != 0) ? children[0] : children[1];
          if (childId != -1) {
            int rowInPrimaryTrackTable = getRowDaughters(childId, tmpIDtrack);
            childIDs = (children[0] != 0) ? std::vector<int>{rowInPrimaryTrackTable, 0} : std::vector<int>{0, rowInPrimaryTrackTable};
          } else {
            childIDs = (children[0] != 0) ? std::vector<int>{-1, 0} : std::vector<int>{0, -1};
          }
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
        if (aod::femtodreamparticle::ParticleType::kV0 == femtoParticle.partType()) {
          // If the order in primary producer is changed of storing first pos, neg daughters and then V0 - this must be updated
          const int rowOfLastTrack = outputParts.lastIndex();
          std::vector<int> childIDs = {rowOfLastTrack - 1, rowOfLastTrack};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
      }
    } else {
      eventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }

  /// process function to create derived data with only collisions containing n tracks
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processCollisionsWithNTracksAndNV0(const o2::aod::FDCollision& col,
                                          const o2::aod::FDParticles& parts)
  {
    eventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedParts = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedData<false>(col, thegroupSelectedParts, thegroupSelectedV0s, parts);
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNTracksAndNV0, "Enable producing data with ppp collisions for data", true);

  /// This function stores accepted collisions in derived data
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param groupSelectedV0s partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table
  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedDataTrkCascade(const o2::aod::FDCollision& col, PartitionType groupSelectedTracks, PartitionType groupSelectedCascades, PartType parts)
  {

    /// check tracks
    int tracksCount = 0;
    int antitracksCount = 0;
    for (const auto& part : groupSelectedTracks) {
      if (part.cut() & 1) {
        antitracksCount++;
      } else {
        tracksCount++;
      }
    }

    /// check Cascades
    int ascadeCount = 0;
    int antiCascadeCount = 0;
    for (const auto& casc : groupSelectedCascades) {
      if ((casc.cut() & kSignPlusMask) == kSignPlusMask) {
        ascadeCount++;
      } else {
        antiCascadeCount++;
      }
    }

    std::vector<int> tmpIDtrack;

    if ((ascadeCount >= confNumberOfCascades && tracksCount >= confNumberOfTracks) || (antiCascadeCount >= confNumberOfCascades && antitracksCount >= confNumberOfTracks)) {
      eventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());

      for (const auto& femtoParticle : parts) {
        if (aod::femtodreamparticle::ParticleType::kTrack == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
          tmpIDtrack.push_back(femtoParticle.index());
        }
        if (aod::femtodreamparticle::ParticleType::kCascadeV0Child == femtoParticle.partType() || aod::femtodreamparticle::ParticleType::kCascadeBachelor == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0, 0};
          const auto& children = femtoParticle.childrenIds();
          int childId = 0;
          if (children[0] != 0) {
            childId = children[0];
          } else if (children[1] != 0) {
            childId = children[1];
          } else if (children[2] != 0) {
            childId = children[2];
          }

          if (childId != -1) {
            int rowInPrimaryTrackTable = getRowDaughters(childId, tmpIDtrack);
            if (children[0] != 0) {
              childIDs = std::vector<int>{rowInPrimaryTrackTable, 0, 0};
            } else if (children[1] != 0) {
              childIDs = std::vector<int>{0, rowInPrimaryTrackTable, 0};
            } else if (children[2] != 0) {
              childIDs = std::vector<int>{0, 0, rowInPrimaryTrackTable};
            }
          } else {
            if (children[0] != 0) {
              childIDs = std::vector<int>{-1, 0, 0};
            } else if (children[1] != 0) {
              childIDs = std::vector<int>{0, -1, 0};
            } else if (children[2] != 0) {
              childIDs = std::vector<int>{0, 0, -1};
            }
          }
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
        if (aod::femtodreamparticle::ParticleType::kCascade == femtoParticle.partType()) {
          // If the order in primary producer is changed of storing first pos, neg daughters and then V0 - this must be updated
          const int rowOfLastTrack = outputParts.lastIndex();
          std::vector<int> childIDs = {rowOfLastTrack - 2, rowOfLastTrack - 1, rowOfLastTrack};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
      }
    } else {
      eventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }

  /// process function to create derived data with only collisions containing n tracks
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processCollisionsWithNTracksAndNCascades(const o2::aod::FDCollision& col,
                                                const o2::aod::FDParticles& parts)
  {
    eventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedParts = selectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedCascades = selectedCascades->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedDataTrkCascade<false>(col, thegroupSelectedParts, thegroupSelectedCascades, parts);
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNTracksAndNCascades, "Enable producing data with tracks and Cascades collisions for data", true);


  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedDataV0Reso(const o2::aod::FDCollision& col, PartitionType groupSelectedV0s, PartitionType groupSelectedResos,  PartType parts)
  {
    //check v0's
    int v0Count = 0;
    int antiV0Count = 0;
    int ResoCount = 0; // no antiparticles

    for (const auto& V0 : groupSelectedV0s) {
      if ((V0.mLambda() > confMinInvMassV0) && (V0.mLambda() < confMaxInvMassV0)) {
        if (confRequireBitmask) {
          if (ncheckbit(V0.cut(), confCutV0SameForAntipart)) {
            const auto& posChild = parts.iteratorAt(V0.index() - 2);
            const auto& negChild = parts.iteratorAt(V0.index() - 1);
            if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
                 (posChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0 &&
                 (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
                 (negChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0)) {
              v0Count++;
            }
          }
        } else {
          v0Count++;
        }
      } else if ((V0.mAntiLambda() > confMinInvMassAntiV0) && (V0.mAntiLambda() < confMaxInvMassAntiV0)) {
        if (confRequireBitmask) {
          if (ncheckbit(V0.cut(), confCutV0SameForAntipart)) {
            const auto& posChild = parts.iteratorAt(V0.index() - 2);
            const auto& negChild = parts.iteratorAt(V0.index() - 1);
            if (((posChild.cut() & confChildPosCutV0) == confChildPosCutV0 &&
                 (posChild.pidcut() & confChildNegTPCBitV0) == confChildNegTPCBitV0 && // exchanged values because checking antiparticle daughters and pid of particles exchange
                 (negChild.cut() & confChildNegCutV0) == confChildNegCutV0 &&
                 (negChild.pidcut() & confChildPosTPCBitV0) == confChildPosTPCBitV0)) { // exchanged values because checking antiparticle daughters and pid of particles exchange
              antiV0Count++;
            }
          }
        } else {
          antiV0Count++;
        }
      }
    }


    for (const auto& reso : groupSelectedResos) {

      if (confRequireBitmask) {

        const auto& posChild = parts.iteratorAt(reso.index() - 2);
        const auto& negChild = parts.iteratorAt(reso.index() - 1);

        if (((posresoChild.cut() & Reso2.daughPosCutBit) == Reso2.daughPosCutBit) &&
            ((negresoChild.cut() & Reso2.daughNegCutBit) == Reso2.daughNegCutBit)) {

          ResoCount++;   
        }
      } else {
        ResoCount++;  
      }

    }

    std::vector<int> tmpIDtrack;

    if ((v0Count >= confNumberOfV0 && ResoCount >= confNumberOfReso) || (antiV0Count >= confNumberOfV0 && confNumberOfReso >= confNumberOfReso)) {
      eventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());
      for (const auto& femtoParticle : parts) {
        if (aod::femtodreamparticle::ParticleType::kResoChild == femtoParticle.partType()) { //
          std::vector<int> childIDs = {0, 0};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
          tmpIDtrack.push_back(femtoParticle.index());
        }
        if (aod::femtodreamparticle::ParticleType::kV0Child == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          const auto& children = femtoParticle.childrenIds();
          int childId = (children[0] != 0) ? children[0] : children[1]; //??
          if (childId != -1) {
            int rowInPrimaryTrackTable = getRowDaughters(childId, tmpIDtrack);
            childIDs = (children[0] != 0) ? std::vector<int>{rowInPrimaryTrackTable, 0} : std::vector<int>{0, rowInPrimaryTrackTable};
          } else {
            childIDs = (children[0] != 0) ? std::vector<int>{-1, 0} : std::vector<int>{0, -1};
          }
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
        if (aod::femtodreamparticle::ParticleType::kV0 == femtoParticle.partType()) {
          // If the order in primary producer is changed of storing first pos, neg daughters and then V0 - this must be updated
          const int rowOfLastTrack = outputParts.lastIndex();
          std::vector<int> childIDs = {rowOfLastTrack - 1, rowOfLastTrack};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }

        if ( (Reso.partType1 == femtoParticle.partType()) ||
             (Reso.partType2 == femtoParticle.partType()) ||
             (Reso.partType3 == femtoParticle.partType()) ||
             (Reso.partType4 == femtoParticle.partType()) ) {

            const int rowOfLastTrack = outputParts.lastIndex();               //überprüfen
            std::vector<int> childIDs = {rowOfLastTrack - 1, rowOfLastTrack}; //überprüfen
            outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
      }
    } else {
      eventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }


   void processCollisionsWithNV0AndNReso(const o2::aod::FDCollision& col,
                                          const o2::aod::FDParticles& parts)
  {
    eventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedResos = selectedResos->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = selectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedDataV0Reso<false>(col, thegroupSelectedV0s, thegroupSelectedResos, parts);
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNV0AndNReso, "Enable producing data with pp collisions for data v0-reso", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoDreamProducerTaskForSpecificAnalysis>(cfgc),
  };
  return workflow;
}
