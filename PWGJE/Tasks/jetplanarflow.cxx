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

#include "Framework/runDataProcessing.h"

struct JetPlanarFlowTask {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> betaSD{"betaSD", 0.0, "SoftDrop beta"};
  Configurable<float> zCutSD{"zCutSD", 0.10, "SoftDrop z cut"};
  Configurable<std::vector<double>> jetRadiiPlot{"jetRadiiPlot", std::vector<double>{0.2, 0.4, 0.6}, "jet resolution parameters"};

  int trackSelection = -1;
  int eventSelection = -1;
  std::string particleSelection;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    auto jetRadiiPlotBins = (std::vector<double>)jetRadiiPlot;
    if (jetRadiiPlotBins.size() > 1) {
      jetRadiiPlotBins.push_back(jetRadiiPlotBins[jetRadiiPlotBins.size() - 1] + (TMath::Abs(jetRadiiPlotBins[jetRadiiPlotBins.size() - 1] - jetRadiiPlotBins[jetRadiiPlotBins.size() - 2])));
    } else {
      jetRadiiPlotBins.push_back(jetRadiiPlotBins[jetRadiiPlotBins.size() - 1] + 0.1);
    }
    registry.add("Thn_tracks", "Thn for tracks", {HistType::kTHnC, {{jetRadiiPlotBins, ""}, {150, 0., 150.}, {100, -1.0, 1.0}, {12, 0.0, 1.2}, {12, 0.0, 1.2}, {10, 0.0, 1.0}, {150, 0., 150.}, {80, -1.0, 7.0}, {100, -1.0, 1.0}, {300, 0.0, 3.0}, {2, -0.5, 1.5}}});
  }
  // jet radii, jet pT, tau1, tau2, jetdR pt, phi', eta', dR', isInJet

  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter trackSubCuts = (aod::jtracksub::pt >= trackPtMin && aod::jtracksub::pt < trackPtMax && aod::jtracksub::eta > trackEtaMin && aod::jtracksub::eta < trackEtaMax && aod::jtracksub::phi >= trackPhiMin && aod::jtracksub::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta > trackEtaMin && aod::jmcparticle::eta < trackEtaMax);

  template <bool isMc, typename T, typename U>
  void fillHistograms(T const& jet, U const& tracks)
  {

    std::vector<float> nSubCASDResults = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, tracks, 2, fastjet::contrib::CA_Axes(), true, zCutSD, betaSD);

    float thetaTrack = -1.0;
    float phiTrack = -1.0;
    float rotationMatrix[3][3];
    float rotationMatrix2D[2][2];

    float jetUnitVector[3] = {float(TMath::Cos(jet.phi()) / TMath::CosH(jet.eta())), float(TMath::Sin(jet.phi()) / TMath::CosH(jet.eta())), float(TMath::SinH(jet.eta()) / TMath::CosH(jet.eta()))};
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

    for (auto& jetConstituent : jet.template tracks_as<U>()) {

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

    for (auto const& track : tracks) {

      if constexpr (isMc) {

        if (particleSelection == "PhysicalPrimary" && !track.isPhysicalPrimary()) {
          continue;
        }
        if (isinf(track.eta())) {
          continue;
        }
      } else {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }
      }

      bool isInJet = false;
      for (auto& jetConstituentId : jet.tracksIds()) {
        if (track.globalIndex() == jetConstituentId) {
          isInJet = true;
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
      thetaTrack = TMath::ACos(pzRotated / TMath::Sqrt((pxRotated * pxRotated) + (pyRotated * pyRotated) + (pzRotated * pzRotated)));
      phiTrack = TMath::ATan2(pyRotatedPrincipleAxis, pxRotatedPrincipleAxis);

      registry.fill(HIST("Thn_tracks"), jet.r() / 100.0, jet.pt(), planarFlow, nSubCASDResults[1], nSubCASDResults[2], nSubCASDResults[0], track.pt(), phiTrack, thetaTrack, jetutilities::deltaR(jet, track), isInJet);
    }
  }

  void processDummy(JetTracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<JetCollisions>::iterator const& collision,
                              soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                              soa::Filtered<JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<false>(jet, tracks);
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsData, "charged jet analysis", false);

  void processChargedJetsEventWiseSubData(soa::Filtered<JetCollisions>::iterator const& collision,
                                          soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets,
                                          soa::Filtered<JetTracksSub> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<false>(jet, tracks);
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsEventWiseSubData, "charged event-wise subtracted jet analysis", false);

  void processChargedJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision,
                             soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                             soa::Filtered<JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<false>(jet, tracks);
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsMCD, "charged detector level jet analysis", false);

  void processChargedJetsMCP(JetMcCollisions const& collision,
                             soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                             soa::Filtered<JetParticles> const& particles)
  {
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillHistograms<true>(jet, particles);
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processChargedJetsMCP, "charged particle level jet analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetPlanarFlowTask>(
    cfgc, TaskName{"jet-planarflow"})};
}
