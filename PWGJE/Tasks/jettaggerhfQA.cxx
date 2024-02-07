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

/// \file jettaggerhfQA.cxx
/// \brief Jet tagging general QA
///
/// \author Hanseo Park <hanseo.park@cern.ch>

#include "TF1.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTagTableData, typename JetTagTableMCD, typename JetTagTableMCP>
struct JetTaggerHFQA {

  // Cut configuration
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  // Binning
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};
  ConfigurableAxis binJetPt{"binJetPt", {201, -0.5f, 200.5f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, -0.5, 99.5}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {801, -400.5f, 400.5f}, ""};
  ConfigurableAxis binImpactParameterXYSignificance{"binImpactParameterXYSignificance", {801, -40.5f, 40.5f}, ""}; // test
  ConfigurableAxis binImpactParameterZ{"binImpactParameterZ", {801, -400.5f, 400.5f}, ""};
  ConfigurableAxis binImpactParameterZSignificance{"binImpactParameterZSignificance", {801, -40.5f, 40.5f}, ""}; // test
  ConfigurableAxis binImpactParameterXYZ{"binImpactParameterXYZ", {2001, -1000.5f, 1000.5f}, ""};
  ConfigurableAxis binImpactParameterXYZSignificance{"binImpactParameterXYZSignificance", {2001, -1000.5f, 1000.5f}, ""};

  // Axis
  AxisSpec jetFlavourAxis = {binJetFlavour, "Jet flavour"};
  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
  AxisSpec etaAxis = {binEta, "#eta"};
  AxisSpec phiAxis = {binPhi, "#phi"};
  AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
  AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
  AxisSpec impactParameterXYAxis = {binImpactParameterXY, "IP_{XY} [#mum]"};
  AxisSpec impactParameterXYSignificanceAxis = {binImpactParameterXYSignificance, "IPs_{XY}"};
  AxisSpec impactParameterZAxis = {binImpactParameterZ, "IP_{Z} [#mum]"};
  AxisSpec impactParameterZSignificanceAxis = {binImpactParameterZSignificance, "IPs_{Z}"};
  AxisSpec impactParameterXYZAxis = {binImpactParameterXYZ, "IP_{XYZ} [#mum]"};
  AxisSpec impactParameterXYZSignificanceAxis = {binImpactParameterXYZSignificance, "IPs_{XYZ}"};

  int eventSelection = -1;
  int trackSelection = -1;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (doprocessIPsData) {
      registry.add("h3_jet_pt_track_pt_track_eta", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {etaAxis}}});
      registry.add("h3_jet_pt_track_pt_track_phi", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {phiAxis}}});
      registry.add("h3_jet_pt_track_pt_impact_parameter_xy", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYAxis}}});
      registry.add("h3_jet_pt_track_pt_impact_parameter_xy_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYSignificanceAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYSignificanceAxis}}});

      registry.add("h3_jet_pt_track_pt_impact_parameter_z", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterZAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_z", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterZAxis}}});
      registry.add("h3_jet_pt_track_pt_impact_parameter_z_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterZSignificanceAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_z_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterZSignificanceAxis}}});
      registry.add("h3_jet_pt_track_pt_impact_parameter_xyz", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYZAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xyz", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYZAxis}}});
      registry.add("h3_jet_pt_track_pt_impact_parameter_xyz_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYZSignificanceAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xyz_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYZSignificanceAxis}}});

      // TC
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance_N1", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYSignificanceAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance_N2", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYSignificanceAxis}}});
      registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance_N3", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYSignificanceAxis}}});
    }
    if (doprocessIPsMCD) {
      registry.add("h3_jet_pt_track_pt_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_track_eta_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {etaAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_track_phi_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {phiAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_impact_parameter_z_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_z_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});

      registry.add("h3_track_pt_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_sign_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_sign_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_impact_parameter_z_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_sign_impact_parameter_z_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_sign_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_sign_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
      registry.add("h3_track_pt_sign_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});

      // TC
      registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
    }
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  // Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);

  using JetTagTracksData = soa::Join<aod::JTracks, aod::JTrackPIs>;
  using JetTagTracksMCD = soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>;
  using OriTracksData = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection>;
  using OriTracksMCD = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::McTrackLabels>;

  template <typename T, typename U, typename V, typename W>
  void fillHistogramIPsData(T const& collision, U const& jets)
  {
    for (auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      std::vector<std::vector<float>> TracksImpXY, TracksSignImpXY, TracksImpXYSig, TracksSignImpXYSig;
      std::vector<std::vector<float>> TracksImpZ, TracksSignImpZ, TracksImpZSig, TracksSignImpZSig;
      std::vector<std::vector<float>> TracksImpXYZ, TracksSignImpXYZ, TracksImpXYZSig, TracksSignImpXYZSig;
      for (auto& jtrack : jet.template tracks_as<JetTagTracksData>()) {
        auto track = jtrack.template track_as<OriTracksData>();

        // General parameters
        registry.fill(HIST("h3_jet_pt_track_pt_track_eta"), jet.pt(), track.pt(), track.eta());
        registry.fill(HIST("h3_jet_pt_track_pt_track_phi"), jet.pt(), track.pt(), track.phi());

        float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpZ, varSignImpZ, varImpZSig, varSignImpZSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        int geoSign = jettaggingutilities::getGeoSign(collision, jet, track);
        varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        varSignImpXY = geoSign * TMath::Abs(track.dcaXY()) * jettaggingutilities::cmTomum;
        varImpXYSig = track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2());
        varSignImpXYSig = geoSign * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2());
        varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
        varSignImpZ = geoSign * TMath::Abs(track.dcaZ()) * jettaggingutilities::cmTomum;
        varImpZSig = track.dcaZ() / TMath::Sqrt(track.sigmaDcaZ2());
        varSignImpZSig = geoSign * TMath::Abs(track.dcaZ()) / TMath::Sqrt(track.sigmaDcaZ2());
        float dcaXYZ, sigmaDcaXYZ2;
        jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(), track.sigmaDcaXY2(), track.sigmaDcaZ2());
        varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
        varSignImpXYZ = geoSign * TMath::Abs(dcaXYZ) * jettaggingutilities::cmTomum;
        varImpXYZSig = dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2);
        varSignImpXYZSig = geoSign * TMath::Abs(dcaXYZ) / TMath::Sqrt(sigmaDcaXYZ2);

        registry.fill(HIST("h3_jet_pt_track_pt_impact_parameter_xy"), jet.pt(), track.pt(), varImpXY);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xy"), jet.pt(), track.pt(), varSignImpXY);
        registry.fill(HIST("h3_jet_pt_track_pt_impact_parameter_xy_significance"), jet.pt(), track.pt(), varImpXYSig);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance"), jet.pt(), track.pt(), varSignImpXYSig);

        registry.fill(HIST("h3_jet_pt_track_pt_impact_parameter_z"), jet.pt(), track.pt(), varImpZ);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_z"), jet.pt(), track.pt(), varSignImpZ);
        registry.fill(HIST("h3_jet_pt_track_pt_impact_parameter_z_significance"), jet.pt(), track.pt(), varImpZSig);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_z_significance"), jet.pt(), track.pt(), varSignImpZSig);
        registry.fill(HIST("h3_jet_pt_track_pt_impact_parameter_xyz"), jet.pt(), track.pt(), varImpXYZ);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xyz"), jet.pt(), track.pt(), varSignImpXYZ);
        registry.fill(HIST("h3_jet_pt_track_pt_impact_parameter_xyz_significance"), jet.pt(), track.pt(), varImpXYZSig);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xyz_significance"), jet.pt(), track.pt(), varSignImpXYZSig);

        TracksImpXY.push_back({varImpXY, track.pt()});
        TracksSignImpXY.push_back({varSignImpXY, track.pt()});
        TracksImpXYSig.push_back({varImpXYSig, track.pt()});
        TracksSignImpXYSig.push_back({varSignImpXYSig, track.pt()});
      }
      auto sortImp = [](const std::vector<float>& a, const std::vector<float>& b) {
        return a[0] > b[0];
      };

      std::sort(TracksImpXY.begin(), TracksImpXY.end(), sortImp);
      std::sort(TracksSignImpXY.begin(), TracksSignImpXY.end(), sortImp);
      std::sort(TracksImpXYSig.begin(), TracksImpXYSig.end(), sortImp);
      std::sort(TracksSignImpXYSig.begin(), TracksSignImpXYSig.end(), sortImp);

      if (TracksImpXY.size() > 0) { // N1
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance_N1"), jet.pt(), TracksSignImpXYSig[0][1], TracksSignImpXYSig[0][0]);
      }

      if (TracksImpXY.size() > 1) { // N2
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance_N2"), jet.pt(), TracksSignImpXYSig[1][1], TracksSignImpXYSig[1][0]);
      }

      if (TracksImpXY.size() > 2) { // N3
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance_N3"), jet.pt(), TracksSignImpXYSig[2][1], TracksSignImpXYSig[2][0]);
      }
    }
  }

  template <typename T, typename U, typename V, typename W>
  void fillHistogramIPsMCD(T const& collision, U const& mcdjets)
  {
    int numberOfJetFlavourSpecies = 6;
    for (auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      std::vector<float> TracksImpXY[numberOfJetFlavourSpecies], TracksSignImpXY[numberOfJetFlavourSpecies], TracksImpXYSig[numberOfJetFlavourSpecies], TracksSignImpXYSig[numberOfJetFlavourSpecies];
      std::vector<float> TracksImpZ[numberOfJetFlavourSpecies], TracksSignImpZ[numberOfJetFlavourSpecies], TracksImpZSig[numberOfJetFlavourSpecies], TracksSignImpZSig[numberOfJetFlavourSpecies];
      std::vector<float> TracksImpXYZ[numberOfJetFlavourSpecies], TracksSignImpXYZ[numberOfJetFlavourSpecies], TracksImpXYZSig[numberOfJetFlavourSpecies], TracksSignImpXYZSig[numberOfJetFlavourSpecies];
      int jetflavour = mcdjet.origin();
      if (jetflavour == JetTaggingSpecies::none) {
        LOGF(debug, "NOT DEFINE JET FLAVOR");
      }
      for (auto& jtrack : mcdjet.template tracks_as<JetTagTracksMCD>()) {
        auto track = jtrack.template track_as<OriTracksMCD>();
        // General parameters
        registry.fill(HIST("h3_jet_pt_track_pt_flavour"), mcdjet.pt(), track.pt(), jetflavour);
        registry.fill(HIST("h3_jet_pt_track_eta_flavour"), mcdjet.pt(), track.eta(), jetflavour);
        registry.fill(HIST("h3_jet_pt_track_phi_flavour"), mcdjet.pt(), track.phi(), jetflavour);

        float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig, varImpZ, varSignImpZ, varImpZSig, varSignImpZSig, varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        int geoSign = jettaggingutilities::getGeoSign(collision, mcdjet, track);
        varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        varSignImpXY = geoSign * TMath::Abs(track.dcaXY()) * jettaggingutilities::cmTomum;
        varImpXYSig = track.dcaXY() / TMath::Sqrt(track.sigmaDcaXY2());
        varSignImpXYSig = geoSign * TMath::Abs(track.dcaXY()) / TMath::Sqrt(track.sigmaDcaXY2());
        varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
        varSignImpZ = geoSign * TMath::Abs(track.dcaZ()) * jettaggingutilities::cmTomum;
        varImpZSig = track.dcaZ() / TMath::Sqrt(track.sigmaDcaZ2());
        varSignImpZSig = geoSign * TMath::Abs(track.dcaZ()) / TMath::Sqrt(track.sigmaDcaZ2());
        float dcaXYZ, sigmaDcaXYZ2;
        jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(), track.sigmaDcaXY2(), track.sigmaDcaZ2());
        varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
        varSignImpXYZ = geoSign * TMath::Abs(dcaXYZ) * jettaggingutilities::cmTomum;
        varImpXYZSig = dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2);
        varSignImpXYZSig = geoSign * TMath::Abs(dcaXYZ) / TMath::Sqrt(sigmaDcaXYZ2);

        registry.fill(HIST("h3_jet_pt_impact_parameter_xy_flavour"), mcdjet.pt(), varImpXY, jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_flavour"), mcdjet.pt(), varSignImpXY, jetflavour);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xy_significance_flavour"), mcdjet.pt(), varImpXYSig, jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour"), mcdjet.pt(), varSignImpXYSig, jetflavour);

        registry.fill(HIST("h3_jet_pt_impact_parameter_z_flavour"), mcdjet.pt(), varImpZ, jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_flavour"), mcdjet.pt(), varSignImpZ, jetflavour);
        registry.fill(HIST("h3_jet_pt_impact_parameter_z_significance_flavour"), mcdjet.pt(), varImpZSig, jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour"), mcdjet.pt(), varSignImpZSig, jetflavour);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xyz_flavour"), mcdjet.pt(), varImpXYZ, jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_flavour"), mcdjet.pt(), varSignImpXYZ, jetflavour);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xyz_significance_flavour"), mcdjet.pt(), varImpXYZSig, jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour"), mcdjet.pt(), varSignImpXYZSig, jetflavour);
        registry.fill(HIST("h3_track_pt_impact_parameter_xy_flavour"), track.pt(), varImpXY, jetflavour);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_flavour"), track.pt(), varSignImpXY, jetflavour);
        registry.fill(HIST("h3_track_pt_impact_parameter_xy_significance_flavour"), track.pt(), varImpXYSig, jetflavour);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_significance_flavour"), track.pt(), varSignImpXYSig, jetflavour);
        registry.fill(HIST("h3_track_pt_impact_parameter_z_flavour"), track.pt(), varImpZ, jetflavour);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_flavour"), track.pt(), varSignImpZ, jetflavour);
        registry.fill(HIST("h3_track_pt_impact_parameter_z_significance_flavour"), track.pt(), varImpZSig, jetflavour);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_significance_flavour"), track.pt(), varSignImpZSig, jetflavour);
        registry.fill(HIST("h3_track_pt_impact_parameter_xyz_flavour"), track.pt(), varImpXYZ, jetflavour);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_flavour"), track.pt(), varSignImpXYZ, jetflavour);
        registry.fill(HIST("h3_track_pt_impact_parameter_xyz_significance_flavour"), track.pt(), varImpXYZSig, jetflavour);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_significance_flavour"), track.pt(), varSignImpXYZSig, jetflavour);

        // For TC
        TracksImpXY[jetflavour].push_back(varImpXY);
        TracksSignImpXY[jetflavour].push_back(varSignImpXY);
        TracksImpXYSig[jetflavour].push_back(varImpXYSig);
        TracksSignImpXYSig[jetflavour].push_back(varSignImpXYSig);
        TracksImpZ[jetflavour].push_back(varImpZ);
        TracksSignImpZ[jetflavour].push_back(varSignImpZ);
        TracksImpZSig[jetflavour].push_back(varImpZSig);
        TracksSignImpZSig[jetflavour].push_back(varSignImpZSig);
        TracksImpXYZ[jetflavour].push_back(varImpXYZ);
        TracksSignImpXYZ[jetflavour].push_back(varSignImpXYZ);
        TracksImpXYZSig[jetflavour].push_back(varImpXYZSig);
        TracksSignImpXYZSig[jetflavour].push_back(varSignImpXYZSig);
      }

      // For TC
      sort(TracksImpXY[jetflavour].begin(), TracksImpXY[jetflavour].end(), std::greater<float>());
      sort(TracksSignImpXY[jetflavour].begin(), TracksSignImpXY[jetflavour].end(), std::greater<float>());
      sort(TracksImpXYSig[jetflavour].begin(), TracksImpXYSig[jetflavour].end(), std::greater<float>());
      sort(TracksSignImpXYSig[jetflavour].begin(), TracksSignImpXYSig[jetflavour].end(), std::greater<float>());
      sort(TracksImpZ[jetflavour].begin(), TracksImpZ[jetflavour].end(), std::greater<float>());
      sort(TracksSignImpZ[jetflavour].begin(), TracksSignImpZ[jetflavour].end(), std::greater<float>());
      sort(TracksImpZSig[jetflavour].begin(), TracksImpZSig[jetflavour].end(), std::greater<float>());
      sort(TracksSignImpZSig[jetflavour].begin(), TracksSignImpZSig[jetflavour].end(), std::greater<float>());
      sort(TracksImpXYZ[jetflavour].begin(), TracksImpXYZ[jetflavour].end(), std::greater<float>());
      sort(TracksSignImpXYZ[jetflavour].begin(), TracksSignImpXYZ[jetflavour].end(), std::greater<float>());
      sort(TracksImpXYZSig[jetflavour].begin(), TracksImpXYZSig[jetflavour].end(), std::greater<float>());
      sort(TracksSignImpXYZSig[jetflavour].begin(), TracksSignImpXYZSig[jetflavour].end(), std::greater<float>());

      if (TracksImpXY[jetflavour].size() > 0) { // N1
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1"), mcdjet.pt(), TracksSignImpXYSig[jetflavour][0], jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1"), mcdjet.pt(), TracksSignImpZSig[jetflavour][0], jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1"), mcdjet.pt(), TracksSignImpXYZSig[jetflavour][0], jetflavour);
      }
      if (TracksImpXY[jetflavour].size() > 1) { // N2
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2"), mcdjet.pt(), TracksSignImpXYSig[jetflavour][1], jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2"), mcdjet.pt(), TracksSignImpZSig[jetflavour][1], jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2"), mcdjet.pt(), TracksSignImpXYZSig[jetflavour][1], jetflavour);
      }
      if (TracksImpXY[jetflavour].size() > 2) { // N3
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2"), mcdjet.pt(), TracksSignImpXYSig[jetflavour][2], jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2"), mcdjet.pt(), TracksSignImpZSig[jetflavour][2], jetflavour);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2"), mcdjet.pt(), TracksSignImpXYZSig[jetflavour][2], jetflavour);
      }
    }
  }

  void processDummy(aod::Collision const& collision, aod::Tracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetTaggerHFQA, processDummy, "Dummy process", true);

  void processIPsData(soa::Join<aod::JCollisions, aod::JCollisionPIs>::iterator const& jcollision, aod::Collisions&, JetTagTableData const& jets, JetTagTracksData const& jtracks, OriTracksData const& tracks)
  {
    auto oricoll = jcollision.template collision_as<aod::Collisions>();
    fillHistogramIPsData<aod::Collision, JetTagTableData, JetTagTracksData, OriTracksData>(oricoll, jets);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsData, "Fill impact parameter inpormation for data jets", false);

  void processIPsMCD(soa::Join<aod::JCollisions, aod::JCollisionPIs>::iterator const& jcollision, aod::Collisions&, JetTagTableMCD const& mcdjets, JetTagTracksMCD const& jtracks, OriTracksMCD const& tracks, JetParticles&)
  {
    auto oricoll = jcollision.template collision_as<aod::Collisions>();
    fillHistogramIPsMCD<aod::Collision, JetTagTableMCD, JetTagTracksMCD, OriTracksMCD>(oricoll, mcdjets);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCD, "Fill impact parameter inpormation for mcd jets", false);
};

using JetTaggerQAChargedDataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetTags>;
using JetTaggerQAChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetTags>;
using JetTaggerQAChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

using JetTaggerQACharged = JetTaggerHFQA<JetTaggerQAChargedDataJets, JetTaggerQAChargedMCDJets, JetTaggerQAChargedMCPJets>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerQACharged>(cfgc,
                                          SetDefaultProcesses{}, TaskName{"jet-taggerhf-qa-charged"}));
  /*
  tasks.emplace_back(
    adaptAnalysisTask<JetTaggerQAFull>(cfgc,
                                         SetDefaultProcesses{}, TaskName{"jet-taggerhf-qa-full"}));

    tasks.emplace_back(
      adaptAnalysisTask<JetTaggerQANeutral>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-taggerhf-qa-neutral"}));
  */
  return WorkflowSpec{tasks};
}
