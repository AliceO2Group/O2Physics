// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
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

#include <vector>
#include <string>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainerThreeBody.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamProducerTaskForSpecificAnalysis {

  SliceCache cache;

  Produces<aod::StoredFDCollisions> outputCollision;
  Produces<aod::StoredFDParticles> outputParts;

  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  float mMassOne = -999, mMassTwo = -999, mMassThree = -999;
  int collisions = 0;

  Configurable<int> ConfNumberOfTracks{"ConfNumberOfTracks", 3, "Number of tracks"};
  Configurable<int> ConfNumberOfV0{"ConfNumberOfV0", 0, "Number of V0"};
  Configurable<int> ConfNumberOfCascades{"ConfNumberOfCascades", 0, "Number of Cascades"};

  /// Track selection
  Configurable<float> ConfPIDthrMom{"ConfPIDthrMom", 1.f, "Momentum threshold from which TPC and TOF are required for PID"};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCPIDBit{"ConfTPCPIDBit", 16, "PID TPC bit from cutCulator "};
  Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCTOFPIDBit{"ConfTPCTOFPIDBit", 8, "PID TPCTOF bit from cutCulator"};

  /// Partition for selected particles
  Partition<aod::FDParticles> SelectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCTOFPIDBit));

  /// V0 selection
  Configurable<float> Conf_minInvMass_V0{"Conf_minInvMass_V0", 1.08, "Minimum invariant mass of V0 (particle)"};
  Configurable<float> Conf_maxInvMass_V0{"Conf_maxInvMass_V0", 1.15, "Maximum invariant mass of V0 (particle)"};
  Configurable<float> Conf_minInvMassAnti_V0{"Conf_minInvMassAnti_V0", 1.08, "Minimum invariant mass of V0 (antiparticle)"};
  Configurable<float> Conf_maxInvMassAnti_V0{"Conf_maxInvMassAnti_V0", 1.15, "Maximum invariant mass of V0 (antiparticle)"};
  /// Cascade selection
  Configurable<float> Conf_minInvMass_Cascade{"Conf_minInvMass_Cascade", 1.2, "Minimum invariant mass of Cascade (particle)"};
  Configurable<float> Conf_maxInvMass_Cascade{"Conf_maxInvMass_Cascade", 1.5, "Maximum invariant mass of Cascade (particle)"};

  // Partition for selected particles
  Partition<aod::FDParticles> SelectedV0s = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0));
  Partition<aod::FDParticles> SelectedCascades = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade));

  HistogramRegistry EventRegistry{"EventRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

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
    EventRegistry.add("hStatistiscs", ";bin;Entries", kTH1F, {{3, 0, 3}});
    // get bit for the collision mask
  }
  /// This function stores accepted collisions in derived data
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param groupSelectedV0s partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table
  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedData(o2::aod::FDCollision& col, PartitionType groupSelectedTracks, PartitionType groupSelectedV0s, PartType parts)
  {
    /// check tracks
    int tracksCount = 0;
    int antitracksCount = 0;
    for (auto& part : groupSelectedTracks) {
      if (part.cut() & 1) {
        antitracksCount++;
      } else {
        tracksCount++;
      }
    }

    /// check V0s
    int V0Count = 0;
    int antiV0Count = 0;
    for (auto& V0 : groupSelectedV0s) {
      if ((V0.mLambda() > Conf_minInvMass_V0) && (V0.mLambda() < Conf_maxInvMass_V0)) {
        V0Count++;
      } else if ((V0.mAntiLambda() > Conf_minInvMassAnti_V0) && (V0.mAntiLambda() < Conf_maxInvMassAnti_V0)) {
        antiV0Count++;
      }
    }

    std::vector<int> tmpIDtrack;

    if ((V0Count >= ConfNumberOfV0 && tracksCount >= ConfNumberOfTracks) || (antiV0Count >= ConfNumberOfV0 && antitracksCount >= ConfNumberOfTracks)) {
      EventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());
      for (auto& femtoParticle : parts) {
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
      EventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }

  /// process function to create derived data with only collisions containing n tracks
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processCollisionsWithNTracksAndNV0(o2::aod::FDCollision& col,
                                          o2::aod::FDParticles& parts)
  {
    EventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedParts = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedData<false>(col, thegroupSelectedParts, thegroupSelectedV0s, parts);
  }
  PROCESS_SWITCH(femtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNTracksAndNV0, "Enable producing data with ppp collisions for data", true);

  /// This function stores accepted collisions in derived data
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param groupSelectedV0s partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table
  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedData_TrkCascade(o2::aod::FDCollision& col, PartitionType groupSelectedTracks, PartitionType groupSelectedCascades, PartType parts)
  {

    /// check tracks
    int tracksCount = 0;
    int antitracksCount = 0;
    for (auto& part : groupSelectedTracks) {
      if (part.cut() & 1) {
        antitracksCount++;
      } else {
        tracksCount++;
      }
    }

    /// check Cascades
    int CascadeCount = 0;
    int antiCascadeCount = 0;
    for (auto& casc : groupSelectedCascades) {
      if ((casc.cut() & kSignPlusMask) == kSignPlusMask) {
        CascadeCount++;
      } else {
        antiCascadeCount++;
      }
    }

    std::vector<int> tmpIDtrack;

    if ((CascadeCount >= ConfNumberOfCascades && tracksCount >= ConfNumberOfTracks) || (antiCascadeCount >= ConfNumberOfCascades && antitracksCount >= ConfNumberOfTracks)) {
      EventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());

      for (auto& femtoParticle : parts) {
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
      EventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }

  /// process function to create derived data with only collisions containing n tracks
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processCollisionsWithNTracksAndNCascades(o2::aod::FDCollision& col,
                                                o2::aod::FDParticles& parts)
  {
    EventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedParts = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedCascades = SelectedCascades->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedData_TrkCascade<false>(col, thegroupSelectedParts, thegroupSelectedCascades, parts);
  }
  PROCESS_SWITCH(femtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNTracksAndNCascades, "Enable producing data with tracks and Cascades collisions for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamProducerTaskForSpecificAnalysis>(cfgc),
  };
  return workflow;
}
