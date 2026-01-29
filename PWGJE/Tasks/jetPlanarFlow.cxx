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
/// \author Wenhui Feng <wenhui.feng@cern.ch>
//

#include "JetDerivedDataUtilities.h"

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/detail/TypeTruncation.h>

#include <TMath.h>

#include "fastjet/contrib/AxesDefinition.hh"

#include <cstdint>
#include <string>
#include <vector>

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
  Configurable<float> trackDRMax{"trackDRMax", 0.8, "maximum track distance from jet axis"};
  Configurable<float> planarTrackPtMin{"planarTrackPtMin", 5.0, "minimum constituent pT in planar flow calculation"};
  Configurable<float> relatedTrackptlow{"relatedTrackptlow", 0.15, "after rotation related track pT low edge"};
  Configurable<float> relatedTrackpthigh{"relatedTrackpthigh", 5.0, "after rotation related track pT high edge"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // jet level configurables
  Configurable<float> jetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
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

    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};
    AxisSpec trackPtAxis = {200, -0.5, 199.5, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec trackEtaAxis = {200, -1.0, 1.0, "#eta"};
    AxisSpec trackphiAxis = {160, -1.0, 7.0, "#varphi"};
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {400, -200., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetEtaAxis = {200, -1.0, 1.0, "#eta"};

    if (doprocessJetsQCData || doprocessJetsQCMCD) {
      registry.add("h_collisions", "number of events;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_collisions_zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h2_centrality_collisions", "event status vs. centrality;entries;centrality", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}});
      registry.add("h_track_pt", "track #it{p}_{T} ; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH1F, {trackPtAxis}});
      registry.add("h2_track_eta_track_phi", "track eta vs. track phi; #eta; #phi; counts", {HistType::kTH2F, {trackEtaAxis, trackphiAxis}});
      registry.add("h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_eta", "jet eta;#eta; counts", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_jet_phi", "jet phi;#phi; counts", {HistType::kTH1F, {trackphiAxis}});
      registry.add("h_leadingjet_pt", "leading jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h2_jet_pt_jet_ntracks", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_leadingjet_pt_ntracks", "leading jet pT vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_jet_pt_track_pt", "jet #it{p}_{T,jet} vs. #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}});
      registry.add("h2_leadingjet_pt_track_pt", "leading jet #it{p}_{T,jet} vs. #it{p}_{T,track}; leading #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}});
    }

    if (doprocessChargedJetsData || doprocessChargedJetsMCD) {
      registry.add("h2_leadingjet_pt_planarflow", "leading jet #it{p}_{T,jet} vs. planar flow; #it{p}_{T,jet} (GeV/#it{c}); planar flow", {HistType::kTH2F, {jetPtAxis, {50, 0, 1}}});
      registry.add("h2_jetpt_Nsub2to1CASD", "leading jet #it{p}_{T,jet} vs. Nsubjettiness tau2/tua1; #it{p}_{T,jet} (GeV/#it{c}); tau2 / tau1", {HistType::kTH2F, {jetPtAxis, {50, 0, 1}}});
      registry.add("h2_trackDR_Nsub2to1CASD", "track Delta R vs. Nsubjettiness tau2/tua1; track to jet axis Delta R ; tau2 / tau1", {HistType::kTH2F, {{100, 0, 1}, {50, 0, 1}}});
      registry.add("h2_jetDR_Nsub2to1CASD", "jet axis Delta R vs. Nsubjettiness tau2/tua1; jet reclustering axis Delta R ; tau2 / tau1", {HistType::kTH2F, {{100, 0, 1}, {50, 0, 1}}});
      registry.add("h3_jet_pt_jetDR_Nsub2to1CASD", "leading jet #it{p}_{T,jet} vs. jet axis Delta R vs. Nsubjettiness tau2/tua1; #it{p}_{T,jet} (GeV/#it{c}); jet reclustering axis Delta R ; tau2 / tau1", {HistType::kTH3F, {jetPtAxis, {100, 0, 1}, {50, 0, 1}}});
      registry.add("h3_track_pt_track_theta_track_phi_incone_beforeRot", "track #it{p}_{T} vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {trackPtAxis, {80, -0.5, 3.5}, {140, -1.0, 7.0}}});
      registry.add("h3_track_pt_track_theta_track_phi_outcone_beforeRot", "track #it{p}_{T} vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {trackPtAxis, {80, -0.5, 3.5}, {140, -1.0, 7.0}}});
      registry.add("h3_track_pt_track_theta_track_phi_incone_afterRot", "track #it{p}_{T} vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {trackPtAxis, {80, -0.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h3_track_pt_track_theta_track_phi_outcone_afterRot", "track #it{p}_{T} vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {trackPtAxis, {80, -0.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h3_track_pt_thetaSinphi_track_thetaCosphi_afterRot", "track #it{p}_{T} vs. track theta * Sin phi vs. track theta * Cos phi; #it{p}_{T,track} (GeV/#it{c}); #theta * Sin #phi; #theta * Cos #phi", {HistType::kTH3F, {trackPtAxis, {140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h2_track_thetaSinphi_track_thetaCosphi_incone_lowpt_afterRot", "track theta * Sin phi vs. track theta * Cos phi; #theta * Sin #phi; #theta * Cos #phi", {HistType::kTH2F, {{140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h2_track_thetaSinphi_track_thetaCosphi_incone_highpt_afterRot", "track theta * Sin phi vs. track theta * Cos phi; #theta * Sin #phi; #theta * Cos #phi", {HistType::kTH2F, {{140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h2_track_thetaSinphi_track_thetaCosphi_outcone_lowpt_afterRot", "track theta * Sin phi vs. track theta * Cos phi; #theta * Sin #phi; #theta * Cos #phi", {HistType::kTH2F, {{140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h2_track_thetaSinphi_track_thetaCosphi_outcone_highpt_afterRot", "track theta * Sin phi vs. track theta * Cos phi; #theta * Sin #phi; #theta * Cos #phi", {HistType::kTH2F, {{140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_2prolong_lowpt", "jet axis Delta R vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {{50, 0, 1}, {140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_2prolong_highpt", "jet axis Delta R vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {{50, 0, 1}, {140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_1prolong_lowpt", "jet axis Delta R vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {{50, 0, 1}, {140, -3.5, 3.5}, {140, -3.5, 3.5}}});
      registry.add("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_1prolong_highpt", "jet axis Delta R vs. track theta vs. track phi; #it{p}_{T,track} (GeV/#it{c}); #theta ; #phi", {HistType::kTH3F, {{50, 0, 1}, {140, -3.5, 3.5}, {140, -3.5, 3.5}}});
    }
  }
  // jet pT, tau2/tau1, jetdR, track pt, phi', eta', dR, isInJet

  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);

  template <bool isMc, bool isAreaSubtracted, typename Tcoll, typename Tjet, typename Ttrk>
  void fillHistograms(Tcoll const& collision, Tjet const& jet, Ttrk const& tracks)
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

    for (auto& jetConstituent : jet.template tracks_as<Ttrk>()) {

      if (jetConstituent.pt() < planarTrackPtMin) {
        continue; // calculate flow eigenvector using high pT constituents
      }
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
    auto NsubTau2to1 = nSubCASDResults[2] / nSubCASDResults[1];

    float jetPt = 0.0;
    if constexpr (isAreaSubtracted) {
      jetPt = jet.pt() - (jet.area() * collision.rho());
    } else {
      jetPt = jet.pt();
    }
    registry.fill(HIST("h2_leadingjet_pt_planarflow"), jetPt, planarFlow);
    registry.fill(HIST("h2_jetpt_Nsub2to1CASD"), jetPt, NsubTau2to1);
    registry.fill(HIST("h2_jetDR_Nsub2to1CASD"), nSubCASDResults[0], NsubTau2to1);
    registry.fill(HIST("h3_jet_pt_jetDR_Nsub2to1CASD"), jetPt, nSubCASDResults[0], NsubTau2to1);

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
    if (theta < 0) {
      theta += TMath::Pi();
    }
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
      for (const auto& jetConstituentId : jet.tracksIds()) {
        if (track.globalIndex() == jetConstituentId) {
          isInJet = 1;
          break;
        }
      }
      float thetaOrigin = TMath::ACos(track.pz() / TMath::Sqrt((track.pt() * track.pt()) + (track.pz() * track.pz())));
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

      float thetaRot = TMath::ACos(pzRotated / TMath::Sqrt((pxRotated * pxRotated) + (pyRotated * pyRotated) + (pzRotated * pzRotated)));
      float phiRot = TMath::ATan2(pyRotatedPrincipleAxis, pxRotatedPrincipleAxis);

      registry.fill(HIST("h2_trackDR_Nsub2to1CASD"), trackDR, NsubTau2to1);
      registry.fill(HIST("h3_track_pt_thetaSinphi_track_thetaCosphi_afterRot"), track.pt(), thetaRot * TMath::Sin(phiRot), thetaRot * TMath::Cos(phiRot));
      if (isInJet) {
        registry.fill(HIST("h3_track_pt_track_theta_track_phi_incone_beforeRot"), track.pt(), thetaOrigin, track.phi());
        registry.fill(HIST("h3_track_pt_track_theta_track_phi_incone_afterRot"), track.pt(), thetaRot, phiRot);
        if (track.pt() > relatedTrackptlow &&  track.pt() < relatedTrackpthigh) { // low pt tracks, no contribute to rotation axis
          registry.fill(HIST("h2_track_thetaSinphi_track_thetaCosphi_incone_lowpt_afterRot"), thetaRot * TMath::Sin(phiRot), thetaRot * TMath::Cos(phiRot));
        } else if (track.pt() > relatedTrackpthigh) {
          registry.fill(HIST("h2_track_thetaSinphi_track_thetaCosphi_incone_highpt_afterRot"), thetaRot * TMath::Sin(phiRot), thetaRot * TMath::Cos(phiRot));
        }
      } else {
        registry.fill(HIST("h3_track_pt_track_theta_track_phi_outcone_beforeRot"), track.pt(), thetaOrigin, track.phi());
        registry.fill(HIST("h3_track_pt_track_theta_track_phi_outcone_afterRot"), track.pt(), thetaRot, phiRot);
        if (track.pt() > relatedTrackptlow &&  track.pt() < relatedTrackpthigh) { // low pt tracks, no contribute to rotation axis
          registry.fill(HIST("h2_track_thetaSinphi_track_thetaCosphi_outcone_lowpt_afterRot"), thetaRot * TMath::Sin(phiRot), thetaRot * TMath::Cos(phiRot));
        } else if (track.pt() > relatedTrackpthigh) {
          registry.fill(HIST("h2_track_thetaSinphi_track_thetaCosphi_outcone_highpt_afterRot"), thetaRot * TMath::Sin(phiRot), thetaRot * TMath::Cos(phiRot));
        }
      }

      if (NsubTau2to1 < 0.3) {
        if (track.pt() > relatedTrackptlow &&  track.pt() < relatedTrackpthigh) {
          registry.fill(HIST("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_2prolong_lowpt"), nSubCASDResults[0], thetaRot * TMath::Sin(track.phi()), thetaRot * TMath::Cos(track.phi()));
        } else if (track.pt() > relatedTrackpthigh) {
          registry.fill(HIST("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_2prolong_highpt"), nSubCASDResults[0], thetaRot * TMath::Sin(track.phi()), thetaRot * TMath::Cos(track.phi()));
        }
      } else if (NsubTau2to1 > 0.6){
        if (track.pt() > relatedTrackptlow &&  track.pt() < relatedTrackpthigh) {
          registry.fill(HIST("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_1prolong_lowpt"), nSubCASDResults[0], thetaRot * TMath::Sin(track.phi()), thetaRot * TMath::Cos(track.phi()));
        } else if (track.pt() > relatedTrackpthigh) {
          registry.fill(HIST("h3_track_thetaSinPhi_track_thetaCosphi_jetDR_1prolong_highpt"), nSubCASDResults[0], thetaRot * TMath::Sin(track.phi()), thetaRot * TMath::Cos(track.phi()));
        }
      }
    }

    if (isMc) {
      jetTableMC(o2::math_utils::detail::truncateFloatFraction(jetPt, precisionMask), o2::math_utils::detail::truncateFloatFraction(planarFlow, precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[2] / nSubCASDResults[1], precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[0], precisionMask), trackPtVector, trackPhiRotVector, trackThetaRotVector, trackDRVector, trackIsInJetVector);
    } else {
      jetTable(o2::math_utils::detail::truncateFloatFraction(jetPt, precisionMask), o2::math_utils::detail::truncateFloatFraction(planarFlow, precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[2] / nSubCASDResults[1], precisionMask), o2::math_utils::detail::truncateFloatFraction(nSubCASDResults[0], precisionMask), trackPtVector, trackPhiRotVector, trackThetaRotVector, trackDRVector, trackIsInJetVector);
    }
  }

  template <typename Tjet>
  void fillJetQCHistograms(Tjet const& jet, bool isjetlead = false, float weight = 1.0)
  {
    if (jet.r() == round(jetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_pt_jet_ntracks"), jet.pt(), jet.tracksIds().size(), weight);
      if (isjetlead) {
        registry.fill(HIST("h_leadingjet_pt"), jet.pt(), weight);
        registry.fill(HIST("h2_leadingjet_pt_ntracks"), jet.pt(), jet.tracksIds().size(), weight);
      }
    }

    for (const auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
      if (isjetlead) {
        registry.fill(HIST("h2_leadingjet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
      }
    }
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processDummy, "Dummy process function turned on by default", true);

  void processJetsQCData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                         aod::JetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centFT0M(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centFT0M(), 1.5);
    registry.fill(HIST("h_collisions_zvertex"), collision.posZ());

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h2_track_eta_track_phi"), track.eta(), track.phi());
    }

    bool isjetlead = true;
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillJetQCHistograms(jet, isjetlead);
      isjetlead = false;
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processJetsQCData, "charged jet QC data ", true);

  void processJetsQCMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                        aod::JetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centFT0M(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centFT0M(), 1.5);
    registry.fill(HIST("h_collisions_zvertex"), collision.posZ());

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h2_track_eta_track_phi"), track.eta(), track.phi());
    }

    bool isjetlead = true;
    for (auto const& jet : jets) {

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax) || jet.pt() < jetPtMin) {
        continue;
      }
      fillJetQCHistograms(jet, isjetlead);
      isjetlead = false;
    }
  }
  PROCESS_SWITCH(JetPlanarFlowTask, processJetsQCMCD, "charged jet QC MCD ", false);

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
