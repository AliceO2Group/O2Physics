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

// task for planar flow analysis
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include <vector>
#include <MathUtils/Utils.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/Logger.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "fastjet/contrib/AxesDefinition.hh"
#include "fastjet/contrib/MeasureDefinition.hh"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace planarflow
{
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(PlanarFlow, planarFlow, float);
DECLARE_SOA_COLUMN(Nsub2to1CASD, nsub2to1CASD, float);
DECLARE_SOA_COLUMN(DeltaRCASD, deltaRCASD, float);
DECLARE_SOA_COLUMN(TrackPt, trackPt, std::vector<float>);
DECLARE_SOA_COLUMN(TrackPhiRot, trackPhiRot, std::vector<float>);
DECLARE_SOA_COLUMN(TrackThetaRot, trackThetaRot, std::vector<float>);
DECLARE_SOA_COLUMN(TrackDR, trackDr, std::vector<float>);
DECLARE_SOA_COLUMN(TrackInJet, trackInJet, std::vector<uint8_t>);
DECLARE_SOA_DYNAMIC_COLUMN(MC, mc, []() -> int { return 0; });
} // namespace planarflow
DECLARE_SOA_TABLE(JetTable, "AOD", "JETTABLE",
                  planarflow::JetPt,
                  planarflow::PlanarFlow,
                  planarflow::Nsub2to1CASD,
                  planarflow::DeltaRCASD,
                  planarflow::TrackPt,
                  planarflow::TrackPhiRot,
                  planarflow::TrackThetaRot,
                  planarflow::TrackDR,
                  planarflow::TrackInJet);

DECLARE_SOA_TABLE(JetTableMC, "AOD", "JETTABLEMC",
                  planarflow::JetPt,
                  planarflow::PlanarFlow,
                  planarflow::Nsub2to1CASD,
                  planarflow::DeltaRCASD,
                  planarflow::TrackPt,
                  planarflow::TrackPhiRot,
                  planarflow::TrackThetaRot,
                  planarflow::TrackDR,
                  planarflow::TrackInJet,
                  planarflow::MC<>);

} // namespace o2::aod

struct JetPlanarFlowTask {

  Produces<aod::JetTable> jetTable;
  Produces<aod::JetTableMC> jetTableMC;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackDRMax{"trackDRMax", 0.8, "maximum track distance from jet axis"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> betaSD{"betaSD", 0.0, "SoftDrop beta"};
  Configurable<float> zCutSD{"zCutSD", 0.10, "SoftDrop z cut"};

  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  std::string particleSelection;

  uint32_t precisionMask;

  std::vector<float> trackPtVector;
  std::vector<float> trackThetaRotVector;
  std::vector<float> trackPhiRotVector;
  std::vector<float> trackDRVector;
  std::vector<uint8_t> trackIsInJetVector;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);
    precisionMask = 0xFFFFFC00;
  }
  // jet pT, tau2/tau1, jetdR, track pt, phi', eta', dR, isInJet

  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <bool isMc, bool isAreaSubtracted, typename T, typename U, typename V>
  void fillHistograms(T const& collision, U const& jet, V const& tracks)
  {

    std::vector<float> nSubCASDResults = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, tracks, 2, fastjet::contrib::CA_Axes(), true, zCutSD, betaSD);

    float rotationMatrix[3][3];
    float rotationMatrix2D[2][2];

    float jetUnitVector[3] = {static_cast<float>(TMath::Cos(jet.phi()) / TMath::CosH(jet.eta())), static_cast<float>(TMath::Sin(jet.phi()) / TMath::CosH(jet.eta())), static_cast<float>(TMath::SinH(jet.eta()) / TMath::CosH(jet.eta()))};
    float magPt = TMath::Sqrt((jetUnitVector[0] * jetUnitVector[0]) + (jetUnitVector[1] * jetUnitVector[1]));
    float cosTheta = jetUnitVector[2];
    float sinTheta = magPt;
    float cosPhi = TMath::Cos(jet.phi());
    float sinPhi = TMath::Sin(jet.phi());

    rotationMatrix[0][0] = -1.0 * cosTheta * cosPhi;
    rotationMatrix[0][1] = -1.0 * cosTheta * sinPhi;
    rotationMatrix[0][2] = sinTheta;
    rotationMatrix[1][0] = sinPhi;
    rotationMatrix[1][1] = -1.0 * cosPhi;
    rotationMatrix[1][2] = 0.;
    rotationMatrix[2][0] = sinTheta * cosPhi;
    rotationMatrix[2][1] = sinTheta * sinPhi;
    rotationMatrix[2][2] = cosTheta;

    float principleMatrix[2][2];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        principleMatrix[i][j] = 0.0;
      }
    }

    for (auto& jetConstituent : jet.template tracks_as<V>()) {

      double normalisation_factor = (1.0 / (jetConstituent.energy() * jet.mass()));
      float pxRotated = (rotationMatrix[0][0] * jetConstituent.px()) + (rotationMatrix[0][1] * jetConstituent.py()) + (rotationMatrix[0][2] * jetConstituent.pz());
      float pyRotated = (rotationMatrix[1][0] * jetConstituent.px()) + (rotationMatrix[1][1] * jetConstituent.py()) + (rotationMatrix[1][2] * jetConstituent.pz());
      // float pzRotated = (rotationMatrix[2][0] * jetConstituent.px()) + (rotationMatrix[2][1] * jetConstituent.py()) + (rotationMatrix[2][2] * jetConstituent.pz());

      principleMatrix[0][0] += normalisation_factor * pxRotated * pxRotated;
      principleMatrix[0][1] += normalisation_factor * pxRotated * pyRotated;
      principleMatrix[1][0] += normalisation_factor * pyRotated * pxRotated;
      principleMatrix[1][1] += normalisation_factor * pyRotated * pyRotated;
    }

    float principleMatrixTrace = principleMatrix[0][0] + principleMatrix[1][1];
    float PrinciplMatrixDeterminant = (principleMatrix[0][0] * principleMatrix[1][1]) - (principleMatrix[0][1] * principleMatrix[1][0]);
    float eigenValue1 = 0.5 * (principleMatrixTrace + TMath::Sqrt(principleMatrixTrace * principleMatrixTrace - 4 * PrinciplMatrixDeterminant));
    float eigenValue2 = 0.5 * (principleMatrixTrace - TMath::Sqrt(principleMatrixTrace * principleMatrixTrace - 4 * PrinciplMatrixDeterminant));

    auto planarFlow = (4.0 * PrinciplMatrixDeterminant) / (principleMatrixTrace * principleMatrixTrace);

    float eigenVector1[2];
    float eigenVector2[2];
    if (principleMatrix[1][0] == 0.0 || principleMatrix[0][1] == 0.0) {
      eigenVector1[0] = principleMatrix[0][0];
      eigenVector1[1] = principleMatrix[1][1];
      eigenVector2[0] = principleMatrix[0][0];
      eigenVector2[1] = principleMatrix[1][1];
    } else {
      eigenVector1[0] = eigenValue1 - principleMatrix[1][1];
      eigenVector1[1] = principleMatrix[1][0];
      eigenVector2[0] = principleMatrix[0][1];
      eigenVector2[1] = eigenValue2 - principleMatrix[0][0];
    }
    if (eigenValue1 < eigenValue2) {
      eigenVector1[0] = eigenVector2[0];
      eigenVector1[1] = eigenVector2[1];
    }

    float theta = TMath::ATan(eigenVector1[1] / eigenVector1[0]);
    if (theta < 0)
      theta += TMath::Pi();
    rotationMatrix2D[0][0] = TMath::Cos(theta);
    rotationMatrix2D[0][1] = TMath::Sin(theta);
    rotationMatrix2D[1][0] = -TMath::Sin(theta);
    rotationMatrix2D[1][1] = TMath::Cos(theta);

    trackPtVector.clear();
    trackThetaRotVector.clear();
    trackPhiRotVector.clear();
    trackDRVector.clear();
    trackIsInJetVector.clear();

    for (auto const& track : tracks) {

      if (track.pt() < trackPtMin || track.pt() >= trackPtMax || track.eta() < trackEtaMin || track.eta() > trackEtaMax) {
        continue;
      }
      auto trackDR = jetutilities::deltaR(jet, track);
      if (trackDR > trackDRMax) {
        continue;
      }

      if constexpr (isMc) {

        if (particleSelection == "PhysicalPrimary" && !track.isPhysicalPrimary()) {
          continue;
        }
        if (std::isinf(track.eta())) {
          continue;
        }
      } else {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }
      }

      uint8_t isInJet = 0;
      for (auto& jetConstituentId : jet.tracksIds()) {
        if (track.globalIndex() == jetConstituentId) {
          isInJet = 1;
          break;
        }
      }

      float pxRotatedPrincipleAxis = 0.0;
      float pyRotatedPrincipleAxis = 0.0;
      float pxRotated = (rotationMatrix[0][0] * track.px()) + (rotationMatrix[0][1] * track.py()) + (rotationMatrix[0][2] * track.pz());
      float pyRotated = (rotationMatrix[1][0] * track.px()) + (rotationMatrix[1][1] * track.py()) + (rotationMatrix[1][2] * track.pz());
      float pzRotated = (rotationMatrix[2][0] * track.px()) + (rotationMatrix[2][1] * track.py()) + (rotationMatrix[2][2] * track.pz());

      pxRotatedPrincipleAxis = (rotationMatrix2D[0][0] * pxRotated) + (rotationMatrix2D[0][1] * pyRotated);
      pyRotatedPrincipleAxis = (rotationMatrix2D[1][0] * pxRotated) + (rotationMatrix2D[1][1] * pyRotated);
      trackPtVector.push_back(o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMask));
      trackThetaRotVector.push_back(o2::math_utils::detail::truncateFloatFraction(TMath::ACos(pzRotated / TMath::Sqrt((pxRotated * pxRotated) + (pyRotated * pyRotated) + (pzRotated * pzRotated))), precisionMask));
      trackPhiRotVector.push_back(o2::math_utils::detail::truncateFloatFraction(TMath::ATan2(pyRotatedPrincipleAxis, pxRotatedPrincipleAxis), precisionMask));
      trackDRVector.push_back(o2::math_utils::detail::truncateFloatFraction(jetutilities::deltaR(jet, track), precisionMask));
      trackIsInJetVector.push_back(isInJet);
    }

    if (isMc) {
      jetTableMC(o2::math_utils::detail::truncateFloatFraction(jet.pt(), precisionMask), o2::math_utils::detail::truncateFloatFraction(planarFlow, precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[2] / nSubCASDResults[1], precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[0], precisionMask), trackPtVector, trackPhiRotVector, trackThetaRotVector, trackDRVector, trackIsInJetVector);
    } else {
      float jetPt = 0.0;
      if constexpr (isAreaSubtracted) {
        jetPt = jet.pt() - (jet.area() * collision.rho());
      } else {
        jetPt = jet.pt();
      }
      jetTable(o2::math_utils::detail::truncateFloatFraction(jetPt, precisionMask), o2::math_utils::detail::truncateFloatFraction(planarFlow, precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[2] / nSubCASDResults[1], precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[0], precisionMask), trackPtVector, trackPhiRotVector, trackThetaRotVector, trackDRVector, trackIsInJetVector);
    }
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                              soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                              aod::JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<false, false>(collision, jet, tracks);
      break; // only fill for the highest pT jet in the collision
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsData, "charged jet analysis", false);

  void processChargedRhoAreaSubtractedJetsData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                               soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                                               aod::JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() - (jet.area() * collision.rho()) < jetPtMin) {
        continue;
      }
      fillHistograms<false, true>(collision, jet, tracks);
      break;
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedRhoAreaSubtractedJetsData, "charged rho-area subtracted jet analysis", false);

  void processChargedJetsEventWiseSubData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                          soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets,
                                          aod::JetTracksSub const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<false, false>(collision, jet, tracks);
      break;
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsEventWiseSubData, "charged event-wise subtracted jet analysis", false);

  void processChargedJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                             soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                             aod::JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<false, false>(collision, jet, tracks);
      break;
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsMCD, "charged detector level jet analysis", false);

  void processChargedJetsMCP(aod::JetMcCollisions const& collision,
                             soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                             aod::JetParticles const& particles)
  {
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<true, false>(collision, jet, particles);
      break;
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsMCP, "charged particle level jet analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetPlanarFlowTask>(
    cfgc, TaskName{"jet-planarflow"})};
}
