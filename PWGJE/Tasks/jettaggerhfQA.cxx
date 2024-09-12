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

template <typename JetTableData, typename TagTableData, typename JetTableMCD, typename weightMCD, typename TagTableMCD, typename JetTableMCP, typename weightMCP, typename TagTableMCP, typename JetTableMCDMCP, typename JetTableMCPMCD>
struct JetTaggerHFQA {

  // task on/off configuration
  Configurable<bool> fillIPxy{"fillIPxy", true, "process of xy plane of dca"};
  Configurable<bool> fillIPz{"fillIPz", false, "process of z plane of dca"};
  Configurable<bool> fillIPxyz{"fillIPxyz", false, "process of xyz plane of dca"};
  Configurable<bool> fillTrackCounting{"fillTrackCounting", false, "process of track counting method"};

  // Cut configuration
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<float> trackDcaXYMax{"trackDcaXYMax", 1, "minimum DCA xy acceptance for tracks [cm]"};
  Configurable<float> trackDcaZMax{"trackDcaZMax", 2, "minimum DCA z acceptance for tracks [cm]"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> prongChi2PCAMin{"prongChi2PCAMin", 1, "minimum Chi2 PCA of decay length of prongs"};
  Configurable<float> prongChi2PCAMax{"prongChi2PCAMax", 100, "maximum Chi2 PCA of decay length of prongs"};
  Configurable<float> prongsigmaLxyMax{"prongsigmaLxyMax", 100, "maximum sigma of decay length of prongs on xy plane"};
  Configurable<float> prongsigmaLxyzMax{"prongsigmaLxyzMax", 100, "maximum sigma of decay length of prongs on xyz plane"};
  Configurable<float> prongIPxyMin{"prongIPxyMin", 0.008, "maximum impact paramter of prongs on xy plane"};
  Configurable<float> prongIPxyMax{"prongIpxyMax", 1, "minimum impact parmeter of prongs on xy plane"};
  Configurable<int> numFlavourSpecies{"numFlavourSpecies", 6, "number of jet flavour species"};
  Configurable<int> numOrder{"numOrder", 6, "number of ordering"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  // Binning
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, 0., 100.}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {801, -400.5f, 400.5f}, ""};
  ConfigurableAxis binSigmaImpactParameterXY{"binImpactSigmaParameterXY", {800, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXYSignificance{"binImpactParameterXYSignificance", {801, -40.5f, 40.5f}, ""};
  ConfigurableAxis binImpactParameterZ{"binImpactParameterZ", {801, -400.5f, 400.5f}, ""};
  ConfigurableAxis binImpactParameterZSignificance{"binImpactParameterZSignificance", {801, -40.5f, 40.5f}, ""};
  ConfigurableAxis binImpactParameterXYZ{"binImpactParameterXYZ", {2001, -1000.5f, 1000.5f}, ""};
  ConfigurableAxis binImpactParameterXYZSignificance{"binImpactParameterXYZSignificance", {2001, -100.5f, 100.5f}, ""};
  ConfigurableAxis binNumOrder{"binNumOrder", {6, 0.5, 6.5}, ""};
  ConfigurableAxis binJetProbability{"binJetProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbabilityLog{"binJetProbabilityLog", {100, 0.f, 10.f}, ""};
  ConfigurableAxis binNprongs{"binNprongs", {100, 0., 100.}, ""};
  ConfigurableAxis binLxy{"binLxy", {200, 0, 20.f}, ""};
  ConfigurableAxis binSxy{"binSxy", {1000, 0, 1000.f}, ""};
  ConfigurableAxis binLxyz{"binLxyz", {200, 0, 20.f}, ""};
  ConfigurableAxis binSxyz{"binSxyz", {1000, 0, 1000.f}, ""};
  ConfigurableAxis binMass{"binMass", {50, 0, 10.f}, ""};
  ConfigurableAxis binSigmaLxy{"binSigmaLxy", {100, 0., 0.1}, ""};
  ConfigurableAxis binSigmaLxyz{"binSigmaLxyz", {100, 0., 0.1}, ""};

  int numberOfJetFlavourSpecies = 6;
  int eventSelection = -1;
  int trackSelection = -1;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    numberOfJetFlavourSpecies = static_cast<int>(numFlavourSpecies);

    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    // Axis
    AxisSpec jetFlavourAxis = {binJetFlavour, "Jet flavour"};
    AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
    AxisSpec etaAxis = {binEta, "#eta"};
    AxisSpec phiAxis = {binPhi, "#phi"};
    AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
    AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
    AxisSpec impactParameterXYAxis = {binImpactParameterXY, "IP_{XY} [#mum]"};
    AxisSpec sigmaImpactParameterXYAxis = {binSigmaImpactParameterXY, "#sigma_{XY} [#mum]"};
    AxisSpec sigmaImpactParameterXYZAxis = {binSigmaImpactParameterXY, "#sigma_{XYZ} [#mum]"};
    AxisSpec impactParameterXYSignificanceAxis = {binImpactParameterXYSignificance, "IPs_{XY}"};
    AxisSpec impactParameterZAxis = {binImpactParameterZ, "IP_{Z} [#mum]"};
    AxisSpec impactParameterZSignificanceAxis = {binImpactParameterZSignificance, "IPs_{Z}"};
    AxisSpec impactParameterXYZAxis = {binImpactParameterXYZ, "IP_{XYZ} [#mum]"};
    AxisSpec impactParameterXYZSignificanceAxis = {binImpactParameterXYZSignificance, "IPs_{XYZ}"};
    AxisSpec numOrderAxis = {binNumOrder, "N_{order}"};
    AxisSpec JetProbabilityAxis = {binJetProbability, "JP"};
    AxisSpec JetProbabilityLogAxis = {binJetProbabilityLog, "-Log(JP)"};
    AxisSpec nprongsAxis = {binNprongs, "#it{N}_{SV}"};
    AxisSpec LxyAxis = {binLxy, "L_{XY} [cm]"};
    AxisSpec SxyAxis = {binSxy, "S_{XY}"};
    AxisSpec LxyzAxis = {binLxyz, "L_{XYZ} [cm]"};
    AxisSpec SxyzAxis = {binSxyz, "S_{XYZ}"};
    AxisSpec massAxis = {binMass, "#it{m}_{SV}"};
    AxisSpec sigmaLxyAxis = {binSigmaLxy, "#sigma_{L_{XY}} [cm]"};
    AxisSpec sigmaLxyzAxis = {binSigmaLxyz, "#sigma_{L_{XYZ}} [cm]"};

    if (doprocessTracksDca) {
      registry.add("h_impact_parameter_xy", "", {HistType::kTH1F, {{impactParameterXYAxis}}});
      registry.add("h_impact_parameter_xy_significance", "", {HistType::kTH1F, {{impactParameterXYSignificanceAxis}}});
      registry.add("h_impact_parameter_z", "", {HistType::kTH1F, {{impactParameterZAxis}}});
      registry.add("h_impact_parameter_z_significance", "", {HistType::kTH1F, {{impactParameterZSignificanceAxis}}});
      registry.add("h_impact_parameter_xyz", "", {HistType::kTH1F, {{impactParameterXYZAxis}}});
      registry.add("h_impact_parameter_xyz_significance", "", {HistType::kTH1F, {{impactParameterXYZSignificanceAxis}}});
    }
    if (doprocessIPsData) {
      registry.add("h3_jet_pt_track_pt_track_eta", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {etaAxis}}});
      registry.add("h3_jet_pt_track_pt_track_phi", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {phiAxis}}});
      if (fillIPxy) {
        registry.add("h2_jet_pt_impact_parameter_xy", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterXYAxis}}});
        registry.add("h2_jet_pt_sign_impact_parameter_xy", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterXYAxis}}});
        registry.add("h2_jet_pt_impact_parameter_xy_significance", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}}});
        registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYSignificanceAxis}}});
      }
      if (fillIPz) {
        registry.add("h2_jet_pt_impact_parameter_z", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterZAxis}}});
        registry.add("h2_jet_pt_sign_impact_parameter_z", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterZAxis}}});
        registry.add("h2_jet_pt_impact_parameter_z_significance", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterZSignificanceAxis}}});
        registry.add("h3_jet_pt_track_pt_sign_impact_parameter_z_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterZSignificanceAxis}}});
      }
      if (fillIPxyz) {
        registry.add("h2_jet_pt_impact_parameter_xyz", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterXYZAxis}}});
        registry.add("h2_jet_pt_sign_impact_parameter_xyz", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterXYZAxis}}});
        registry.add("h2_jet_pt_impact_parameter_xyz_significance", "", {HistType::kTH2F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}}});
        registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xyz_significance", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {impactParameterXYZSignificanceAxis}}});
      }
      if (fillTrackCounting) {
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_tc", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {numOrderAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_significance_tc", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {numOrderAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_tc", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {numOrderAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_significance_tc", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {numOrderAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_z_significance_tc", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZSignificanceAxis}, {numOrderAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xyz_significance_tc", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZSignificanceAxis}, {numOrderAxis}}});
      }
    }
    if (doprocessIPsMCD || doprocessIPsMCDWeighted) {
      registry.add("h2_jet_pt_flavour", "", {HistType::kTH2F, {{jetPtAxis}, {jetFlavourAxis}}});
      registry.add("h2_jet_eta_flavour", "", {HistType::kTH2F, {{etaAxis}, {jetFlavourAxis}}});
      registry.add("h2_jet_phi_flavour", "", {HistType::kTH2F, {{phiAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_track_pt_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {trackPtAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_track_eta_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {etaAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_track_phi_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {phiAxis}, {jetFlavourAxis}}});
      if (fillIPxy) {
        registry.add("h3_jet_pt_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sigma_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaImpactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      }
      if (fillIPz) {
        registry.add("h3_jet_pt_impact_parameter_z_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_z_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_z_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
      }
      if (fillIPxyz) {
        registry.add("h3_jet_pt_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
      }
      if (fillTrackCounting) {
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYZSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_sign_impact_parameter_xy_significance_tc_flavour", "", {HistType::kTH3F, {{impactParameterXYSignificanceAxis}, {numOrderAxis}, {jetFlavourAxis}}});
        registry.add("h3_sign_impact_parameter_z_significance_tc_flavour", "", {HistType::kTH3F, {{impactParameterZSignificanceAxis}, {numOrderAxis}, {jetFlavourAxis}}});
        registry.add("h3_sign_impact_parameter_xyz_significance_tc_flavour", "", {HistType::kTH3F, {{impactParameterXYZSignificanceAxis}, {numOrderAxis}, {jetFlavourAxis}}});
      }
    }
    if (doprocessIPsMCP || doprocessIPsMCPWeighted) {
      registry.add("h2_jet_pt_part_flavour", "", {HistType::kTH2F, {{jetPtAxis}, {jetFlavourAxis}}});
      registry.add("h2_jet_eta_part_flavour", "", {HistType::kTH2F, {{etaAxis}, {jetFlavourAxis}}});
      registry.add("h2_jet_phi_part_flavour", "", {HistType::kTH2F, {{phiAxis}, {jetFlavourAxis}}});
    }

    if (doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) {
      registry.add("h3_jet_pt_jet_pt_part_matchedgeo_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {jetPtAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_jet_pt_part_matchedgeo_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {jetPtAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_flavour_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {jetFlavourAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_eta_flavour_flavour_run2", "", {HistType::kTH3F, {{etaAxis}, {jetFlavourAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_phi_flavour_flavour_run2", "", {HistType::kTH3F, {{phiAxis}, {jetFlavourAxis}, {jetFlavourAxis}}});
      registry.add("h_compare_flavour_flavour_run2", "", {HistType::kTH1F, {{2, 0, 2}}});
      if (fillIPxy) {
        registry.add("h3_jet_pt_impact_parameter_xy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sigma_impact_parameter_xy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaImpactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_impact_parameter_xy_significance_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_xy_flavour_run2", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_flavour_run2", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_impact_parameter_xy_significance_flavour_run2", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_significance_flavour_run2", "", {HistType::kTH3F, {{trackPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      }
      if (fillTrackCounting) {
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2_N1", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2_N2", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2_N3", "", {HistType::kTH3F, {{jetPtAxis}, {impactParameterXYSignificanceAxis}, {jetFlavourAxis}}});
      }
    }

    if (doprocessJPData) {
      registry.add("h2_jet_pt_JP", "jet pt jet probability untagged", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityAxis}}});
      registry.add("h2_jet_pt_neg_log_JP", "jet pt jet probabilityun tagged", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityLogAxis}}});
      registry.add("h2_jet_pt_JP_N1", "jet pt jet probability N1", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityAxis}}});
      registry.add("h2_jet_pt_neg_log_JP_N1", "jet pt jet probabilityun N1", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityLogAxis}}});
      registry.add("h2_jet_pt_JP_N2", "jet pt jet probability N2", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityAxis}}});
      registry.add("h2_jet_pt_neg_log_JP_N2", "jet pt jet probabilityun N2", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityLogAxis}}});
      registry.add("h2_jet_pt_JP_N3", "jet pt jet probability N3", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityAxis}}});
      registry.add("h2_jet_pt_neg_log_JP_N3", "jet pt jet probabilityun N3", {HistType::kTH2F, {{jetPtAxis}, {JetProbabilityLogAxis}}});
    }
    if (doprocessJPMCD || doprocessJPMCDWeighted) {
      registry.add("h3_jet_pt_JP_flavour", "jet pt jet probability flavour untagged", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_neg_log_JP_flavour", "jet pt log jet probability flavour untagged", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityLogAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_JP_N1_flavour", "jet pt jet probability flavour N1", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_neg_log_JP_N1_flavour", "jet pt log jet probability flavour N1", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityLogAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_JP_N2_flavour", "jet pt jet probability flavour N2", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_neg_log_JP_N2_flavour", "jet pt log jet probability flavour N2", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityLogAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_JP_N3_flavour", "jet pt jet probability flavour N3", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_neg_log_JP_N3_flavour", "jet pt log jet probability flavour N3", {HistType::kTH3F, {{jetPtAxis}, {JetProbabilityLogAxis}, {jetFlavourAxis}}});
    }
    if (doprocessSV2ProngData) {
      registry.add("h_2prong_nprongs", "", {HistType::kTH1F, {{nprongsAxis}}});
      registry.add("h2_jet_pt_2prong_Lxy", "", {HistType::kTH2F, {{jetPtAxis}, {LxyAxis}}});
      registry.add("h2_jet_pt_2prong_sigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyAxis}}});
      registry.add("h2_jet_pt_2prong_Sxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_2prong_Lxyz", "", {HistType::kTH2F, {{jetPtAxis}, {LxyzAxis}}});
      registry.add("h2_jet_pt_2prong_sigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyzAxis}}});
      registry.add("h2_jet_pt_2prong_Sxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_2prong_Sxy_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_2prong_Sxyz_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_2prong_mass_N1", "", {HistType::kTH2F, {{jetPtAxis}, {massAxis}}});
    }
    if (doprocessSV3ProngData) {
      registry.add("h_3prong_nprongs", "", {HistType::kTH1F, {{nprongsAxis}}});
      registry.add("h2_jet_pt_3prong_Lxy", "", {HistType::kTH2F, {{jetPtAxis}, {LxyAxis}}});
      registry.add("h2_jet_pt_3prong_sigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyAxis}}});
      registry.add("h2_jet_pt_3prong_Sxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_3prong_Lxyz", "", {HistType::kTH2F, {{jetPtAxis}, {LxyzAxis}}});
      registry.add("h2_jet_pt_3prong_sigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyzAxis}}});
      registry.add("h2_jet_pt_3prong_Sxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_3prong_Sxy_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_3prong_Sxyz_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_3prong_mass_N1", "", {HistType::kTH2F, {{jetPtAxis}, {massAxis}}});
    }
    if (doprocessSV2ProngMCD || doprocessSV2ProngMCDWeighted) {
      registry.add("h2_2prong_nprongs_flavour", "", {HistType::kTH2F, {{nprongsAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Lxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_sigmaLxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Lxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_sigmaLxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_mass_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_2prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_2prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_2prong_mass_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
    }
    if (doprocessSV3ProngMCD || doprocessSV3ProngMCDWeighted) {
      registry.add("h2_3prong_nprongs_flavour", "", {HistType::kTH2F, {{nprongsAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Lxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_sigmaLxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Lxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_sigmaLxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_mass_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_3prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_3prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_3prong_mass_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
    }
    if (doprocessSV2ProngMCPMCDMatched || doprocessSV2ProngMCPMCDMatchedWeighted) {
      registry.add("h3_jet_pt_2prong_Lxy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {LxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_sigmaLxy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Lxyz_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {LxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_sigmaLxyz_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_mass_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
    }
    if (doprocessSV3ProngMCPMCDMatched || doprocessSV3ProngMCPMCDMatchedWeighted) {
      registry.add("h3_jet_pt_3prong_Lxy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {LxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_sigmaLxy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Lxyz_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {LxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_sigmaLxyz_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_mass_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_3prong_Sxy_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_3prong_Sxyz_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_taggedjet_pt_3prong_mass_N1_flavour_run2", "", {HistType::kTH3F, {{jetPtAxis}, {massAxis}, {jetFlavourAxis}}});
    }
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  PresliceUnsorted<soa::Filtered<JetCollisionsMCD>> CollisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetParticles> particlesPerCollision = aod::jmcparticle::mcCollisionId;

  using JetTagTracksData = soa::Join<JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using JetTagTracksMCD = soa::Join<JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;

  std::function<bool(const std::vector<float>&, const std::vector<float>&)> sortImp =
    [](const std::vector<float>& a, const std::vector<float>& b) {
      return a[0] > b[0];
    };

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {
    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    if (leadingConstituentPtMin > -98.0) {
      bool isMinleadingConstituent = false;
      for (auto& constituent : jet.template tracks_as<T>()) {
        if (constituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }
      if (!isMinleadingConstituent) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool trackAcceptance(T const& track)
  {
    if (track.pt() < trackPtMin || track.pt() > trackPtMax)
      return false;

    return true;
  }

  template <typename T, typename U>
  void fillHistogramIPsData(T const& jet, U const& /*jtracks*/)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    std::vector<std::vector<float>> vecSignImpXYSig, vecSignImpZSig, vecSignImpXYZSig;
    for (auto& track : jet.template tracks_as<U>()) {
      if (!trackAcceptance(track))
        continue;
      if (!jettaggingutilities::trackAcceptanceWithDca(track, trackDcaXYMax, trackDcaZMax))
        continue;
      // General parameters
      registry.fill(HIST("h3_jet_pt_track_pt_track_eta"), jet.pt(), track.pt(), track.eta());
      registry.fill(HIST("h3_jet_pt_track_pt_track_phi"), jet.pt(), track.pt(), track.phi());
      int geoSign = jettaggingutilities::getGeoSign(jet, track);
      if (fillIPxy) {
        float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig;
        varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        varSignImpXY = geoSign * std::abs(track.dcaXY()) * jettaggingutilities::cmTomum;
        varImpXYSig = track.dcaXY() / track.sigmadcaXY();
        varSignImpXYSig = geoSign * std::abs(track.dcaXY()) / track.sigmadcaXY();
        registry.fill(HIST("h2_jet_pt_impact_parameter_xy"), jet.pt(), varImpXY);
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xy"), jet.pt(), varSignImpXY);
        registry.fill(HIST("h2_jet_pt_impact_parameter_xy_significance"), jet.pt(), varImpXYSig);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance"), jet.pt(), track.pt(), varSignImpXYSig);
        vecSignImpXYSig.push_back({varSignImpXYSig, track.pt()});
      }
      if (fillIPz) {
        float varImpZ, varSignImpZ, varImpZSig, varSignImpZSig;
        varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
        varSignImpZ = geoSign * std::abs(track.dcaZ()) * jettaggingutilities::cmTomum;
        varImpZSig = track.dcaZ() / track.sigmadcaZ();
        varSignImpZSig = geoSign * std::abs(track.dcaZ()) / track.sigmadcaZ();
        registry.fill(HIST("h2_jet_pt_impact_parameter_z"), jet.pt(), varImpZ);
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_z"), jet.pt(), varSignImpZ);
        registry.fill(HIST("h2_jet_pt_impact_parameter_z_significance"), jet.pt(), varImpZSig);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_z_significance"), jet.pt(), track.pt(), varSignImpZSig);
        vecSignImpZSig.push_back({varSignImpZSig, track.pt()});
      }
      if (fillIPxyz) {
        float varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        varImpXYZ = track.dcaXYZ() * jettaggingutilities::cmTomum;
        varSignImpXYZ = geoSign * std::abs(track.dcaXYZ()) * jettaggingutilities::cmTomum;
        varImpXYZSig = track.dcaXYZ() / track.sigmadcaXYZ();
        varSignImpXYZSig = geoSign * std::abs(track.dcaXYZ()) / track.sigmadcaXYZ();
        registry.fill(HIST("h2_jet_pt_impact_parameter_xyz"), jet.pt(), varImpXYZ);
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xyz"), jet.pt(), varSignImpXYZ);
        registry.fill(HIST("h2_jet_pt_impact_parameter_xyz_significance"), jet.pt(), varImpXYZSig);
        registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xyz_significance"), jet.pt(), track.pt(), varSignImpXYZSig);
        vecSignImpXYZSig.push_back({varSignImpXYZSig, track.pt()});
      }
    }

    if (!fillTrackCounting)
      return;
    if (fillIPxy)
      std::sort(vecSignImpXYSig.begin(), vecSignImpXYSig.end(), sortImp);
    if (fillIPz)
      std::sort(vecSignImpZSig.begin(), vecSignImpZSig.end(), sortImp);
    if (fillIPxyz)
      std::sort(vecSignImpXYZSig.begin(), vecSignImpXYZSig.end(), sortImp);
    if (fillIPxy && vecSignImpXYSig.empty())
      return;
    if (fillIPz && vecSignImpZSig.empty())
      return;
    if (fillIPxyz && vecSignImpXYZSig.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSig.size()) {
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_tc"), jet.pt(), vecSignImpXYSig[order - 1][0], order);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_significance_tc"), vecSignImpXYSig[order - 1][1], vecSignImpXYSig[order - 1][0], order);
      }
      if (fillIPz && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpZSig.size()) {
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_tc"), jet.pt(), vecSignImpZSig[order - 1][0], order);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_significance_tc"), vecSignImpZSig[order - 1][1], vecSignImpZSig[order - 1][0], order);
      }
      if (fillIPxyz && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYZSig.size()) {
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_tc"), jet.pt(), vecSignImpXYZSig[order - 1][0], order);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_significance_tc"), vecSignImpXYZSig[order - 1][1], vecSignImpXYZSig[order - 1][0], order);
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramIPsMCD(T const& mcdjet, U const& /*jtracks*/, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    std::vector<float> vecImpXY[numberOfJetFlavourSpecies], vecSignImpXY[numberOfJetFlavourSpecies], vecImpXYSig[numberOfJetFlavourSpecies], vecSignImpXYSig[numberOfJetFlavourSpecies];
    std::vector<float> vecImpZ[numberOfJetFlavourSpecies], vecSignImpZ[numberOfJetFlavourSpecies], vecImpZSig[numberOfJetFlavourSpecies], vecSignImpZSig[numberOfJetFlavourSpecies];
    std::vector<float> vecImpXYZ[numberOfJetFlavourSpecies], vecSignImpXYZ[numberOfJetFlavourSpecies], vecImpXYZSig[numberOfJetFlavourSpecies], vecSignImpXYZSig[numberOfJetFlavourSpecies];
    std::vector<std::vector<float>> vecSignImpXYSigTC, vecSignImpZSigTC, vecSignImpXYZSigTC;
    int jetflavour = mcdjet.origin();
    if (jetflavour == JetTaggingSpecies::none) {
      LOGF(debug, "NOT DEFINE JET FLAVOR");
    }
    registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), jetflavour, eventWeight);
    registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), jetflavour, eventWeight);
    registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), jetflavour, eventWeight);
    for (auto& track : mcdjet.template tracks_as<U>()) {
      if (!trackAcceptance(track))
        continue;
      if (!jettaggingutilities::trackAcceptanceWithDca(track, trackDcaXYMax, trackDcaZMax))
        continue;

      // General parameters
      registry.fill(HIST("h3_jet_pt_track_pt_flavour"), mcdjet.pt(), track.pt(), jetflavour, eventWeight);
      registry.fill(HIST("h3_jet_pt_track_eta_flavour"), mcdjet.pt(), track.eta(), jetflavour, eventWeight);
      registry.fill(HIST("h3_jet_pt_track_phi_flavour"), mcdjet.pt(), track.phi(), jetflavour, eventWeight);
      int geoSign = jettaggingutilities::getGeoSign(mcdjet, track);
      if (fillIPxy) {
        float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig;
        varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        float varSigmaImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        varSignImpXY = geoSign * std::abs(track.dcaXY()) * jettaggingutilities::cmTomum;
        varImpXYSig = track.dcaXY() / track.sigmadcaXY();
        varSignImpXYSig = geoSign * std::abs(track.dcaXY()) / track.sigmadcaXY();
        registry.fill(HIST("h3_jet_pt_impact_parameter_xy_flavour"), mcdjet.pt(), varImpXY, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sigma_impact_parameter_xy_flavour"), mcdjet.pt(), varSigmaImpXY, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_flavour"), mcdjet.pt(), varSignImpXY, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xy_significance_flavour"), mcdjet.pt(), varImpXYSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour"), mcdjet.pt(), varSignImpXYSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_xy_flavour"), track.pt(), varImpXY, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_flavour"), track.pt(), varSignImpXY, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_xy_significance_flavour"), track.pt(), varImpXYSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_significance_flavour"), track.pt(), varSignImpXYSig, jetflavour, eventWeight);
        vecImpXY[jetflavour].push_back(varImpXY);
        vecSignImpXY[jetflavour].push_back(varSignImpXY);
        vecImpXYSig[jetflavour].push_back(varImpXYSig);
        vecSignImpXYSig[jetflavour].push_back(varSignImpXYSig);
        vecSignImpXYSigTC.push_back({varSignImpXYSig, track.pt()});
      }
      if (fillIPz) {
        float varImpZ, varSignImpZ, varImpZSig, varSignImpZSig;
        varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
        varSignImpZ = geoSign * std::abs(track.dcaZ()) * jettaggingutilities::cmTomum;
        varImpZSig = track.dcaZ() / track.sigmadcaZ();
        varSignImpZSig = geoSign * std::abs(track.dcaZ()) / track.sigmadcaZ();
        registry.fill(HIST("h3_jet_pt_impact_parameter_z_flavour"), mcdjet.pt(), varImpZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_flavour"), mcdjet.pt(), varSignImpZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_impact_parameter_z_significance_flavour"), mcdjet.pt(), varImpZSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour"), mcdjet.pt(), varSignImpZSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_z_flavour"), track.pt(), varImpZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_flavour"), track.pt(), varSignImpZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_z_significance_flavour"), track.pt(), varImpZSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_significance_flavour"), track.pt(), varSignImpZSig, jetflavour, eventWeight);
        vecImpZ[jetflavour].push_back(varImpZ);
        vecSignImpZ[jetflavour].push_back(varSignImpZ);
        vecImpZSig[jetflavour].push_back(varImpZSig);
        vecSignImpZSig[jetflavour].push_back(varSignImpZSig);
        vecSignImpZSigTC.push_back({varSignImpZSig, track.pt()});
      }
      if (fillIPxyz) {
        float varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
        float dcaXYZ = track.dcaXYZ();
        float sigmadcaXYZ2 = track.sigmadcaXYZ();
        varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
        varSignImpXYZ = geoSign * std::abs(dcaXYZ) * jettaggingutilities::cmTomum;
        varImpXYZSig = dcaXYZ / std::sqrt(sigmadcaXYZ2);
        varSignImpXYZSig = geoSign * std::abs(dcaXYZ) / std::sqrt(sigmadcaXYZ2);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xyz_flavour"), mcdjet.pt(), varImpXYZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_flavour"), mcdjet.pt(), varSignImpXYZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xyz_significance_flavour"), mcdjet.pt(), varImpXYZSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour"), mcdjet.pt(), varSignImpXYZSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_xyz_flavour"), track.pt(), varImpXYZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_flavour"), track.pt(), varSignImpXYZ, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_xyz_significance_flavour"), track.pt(), varImpXYZSig, jetflavour, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_significance_flavour"), track.pt(), varSignImpXYZSig, jetflavour, eventWeight);
        vecImpXYZ[jetflavour].push_back(varImpXYZ);
        vecSignImpXYZ[jetflavour].push_back(varSignImpXYZ);
        vecImpXYZSig[jetflavour].push_back(varImpXYZSig);
        vecSignImpXYZSig[jetflavour].push_back(varSignImpXYZSig);
        vecSignImpXYZSigTC.push_back({varSignImpXYZSig, track.pt()});
      }
    }

    if (!fillTrackCounting)
      return;
    sort(vecImpXY[jetflavour].begin(), vecImpXY[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXY[jetflavour].begin(), vecSignImpXY[jetflavour].end(), std::greater<float>());
    sort(vecImpXYSig[jetflavour].begin(), vecImpXYSig[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXYSig[jetflavour].begin(), vecSignImpXYSig[jetflavour].end(), std::greater<float>());
    sort(vecImpZ[jetflavour].begin(), vecImpZ[jetflavour].end(), std::greater<float>());
    sort(vecSignImpZ[jetflavour].begin(), vecSignImpZ[jetflavour].end(), std::greater<float>());
    sort(vecImpZSig[jetflavour].begin(), vecImpZSig[jetflavour].end(), std::greater<float>());
    sort(vecSignImpZSig[jetflavour].begin(), vecSignImpZSig[jetflavour].end(), std::greater<float>());
    sort(vecImpXYZ[jetflavour].begin(), vecImpXYZ[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXYZ[jetflavour].begin(), vecSignImpXYZ[jetflavour].end(), std::greater<float>());
    sort(vecImpXYZSig[jetflavour].begin(), vecImpXYZSig[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXYZSig[jetflavour].begin(), vecSignImpXYZSig[jetflavour].end(), std::greater<float>());

    if (vecImpXY[jetflavour].size() > 0) { // N1
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1"), mcdjet.pt(), vecSignImpXYSig[jetflavour][0], jetflavour, eventWeight);
      if (fillIPz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1"), mcdjet.pt(), vecSignImpZSig[jetflavour][0], jetflavour, eventWeight);
      if (fillIPxyz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][0], jetflavour, eventWeight);
    }
    if (vecImpXY[jetflavour].size() > 1) { // N2
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2"), mcdjet.pt(), vecSignImpXYSig[jetflavour][1], jetflavour, eventWeight);
      if (fillIPz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2"), mcdjet.pt(), vecSignImpZSig[jetflavour][1], jetflavour, eventWeight);
      if (fillIPxyz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][1], jetflavour, eventWeight);
    }
    if (vecImpXY[jetflavour].size() > 2) { // N3
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N3"), mcdjet.pt(), vecSignImpXYSig[jetflavour][2], jetflavour, eventWeight);
      if (fillIPz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N3"), mcdjet.pt(), vecSignImpZSig[jetflavour][2], jetflavour, eventWeight);
      if (fillIPxyz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N3"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][2], jetflavour, eventWeight);
    }

    std::sort(vecSignImpXYSigTC.begin(), vecSignImpXYSigTC.end(), sortImp);
    std::sort(vecSignImpZSigTC.begin(), vecSignImpZSigTC.end(), sortImp);
    std::sort(vecSignImpXYZSigTC.begin(), vecSignImpXYZSigTC.end(), sortImp);

    if (vecSignImpXYSigTC.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSigTC.size()) {
        registry.fill(HIST("h3_sign_impact_parameter_xy_significance_tc_flavour"), vecSignImpXYSigTC[order - 1][0], order, jetflavour, eventWeight);
      }
    }
    if (vecSignImpZSigTC.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpZSigTC.size()) {
        registry.fill(HIST("h3_sign_impact_parameter_z_significance_tc_flavour"), vecSignImpZSigTC[order - 1][0], order, jetflavour, eventWeight);
      }
    }

    if (vecSignImpXYZSigTC.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYZSigTC.size()) {
        registry.fill(HIST("h3_sign_impact_parameter_xyz_significance_tc_flavour"), vecSignImpXYZSigTC[order - 1][0], order, jetflavour, eventWeight);
      }
    }
  }

  template <typename T>
  void fillHistogramIPsMCP(T const& mcpjet, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcpjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    int jetflavour = mcpjet.origin();
    if (jetflavour == JetTaggingSpecies::none) {
      LOGF(debug, "NOT DEFINE JET FLAVOR");
    }
    registry.fill(HIST("h2_jet_pt_part_flavour"), mcpjet.pt(), jetflavour, eventWeight);
    registry.fill(HIST("h2_jet_eta_part_flavour"), mcpjet.eta(), jetflavour, eventWeight);
    registry.fill(HIST("h2_jet_phi_part_flavour"), mcpjet.phi(), jetflavour, eventWeight);
  }

  template <typename T, typename U, typename V, typename W>
  void fillHistogramIPsMatched(T const& mcdjet, U const&, V const&, W const& particlesPerColl, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    std::vector<float> vecImpXY[numberOfJetFlavourSpecies], vecSignImpXY[numberOfJetFlavourSpecies], vecImpXYSig[numberOfJetFlavourSpecies], vecSignImpXYSig[numberOfJetFlavourSpecies];
    int jetflavour = mcdjet.origin();
    if (jetflavour == JetTaggingSpecies::none)
      jetflavour = JetTaggingSpecies::lightflavour;
    int jetflavourRun2Def = -1;
    // if (!mcdjet.has_matchedJetGeo()) continue;
    for (auto& mcpjet : mcdjet.template matchedJetGeo_as<U>()) {
      jetflavourRun2Def = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
      registry.fill(HIST("h3_jet_pt_jet_pt_part_matchedgeo_flavour"), mcpjet.pt(), mcdjet.pt(), jetflavour, eventWeight);
      registry.fill(HIST("h3_jet_pt_jet_pt_part_matchedgeo_flavour_run2"), mcpjet.pt(), mcdjet.pt(), jetflavourRun2Def, eventWeight);
    }
    if (jetflavourRun2Def < 0)
      return;
    if (jetflavour == jetflavourRun2Def)
      registry.fill(HIST("h_compare_flavour_flavour_run2"), 0.5);
    else
      registry.fill(HIST("h_compare_flavour_flavour_run2"), 1.5);
    registry.fill(HIST("h3_jet_pt_flavour_flavour_run2"), mcdjet.pt(), jetflavour, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_jet_eta_flavour_flavour_run2"), mcdjet.eta(), jetflavour, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_jet_phi_flavour_flavour_run2"), mcdjet.phi(), jetflavour, jetflavourRun2Def, eventWeight);
    for (auto& track : mcdjet.template tracks_as<V>()) {
      int geoSign = jettaggingutilities::getGeoSign(mcdjet, track);
      if (fillIPxy) {
        float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig;
        varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        float varSigmaImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        varSignImpXY = geoSign * std::abs(track.dcaXY()) * jettaggingutilities::cmTomum;
        varImpXYSig = track.dcaXY() / track.sigmadcaXY();
        varSignImpXYSig = geoSign * std::abs(track.dcaXY()) / track.sigmadcaXY();
        registry.fill(HIST("h3_jet_pt_impact_parameter_xy_flavour_run2"), mcdjet.pt(), varImpXY, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_jet_pt_sigma_impact_parameter_xy_flavour_run2"), mcdjet.pt(), varSigmaImpXY, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_flavour_run2"), mcdjet.pt(), varSignImpXY, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_jet_pt_impact_parameter_xy_significance_flavour_run2"), mcdjet.pt(), varImpXYSig, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2"), mcdjet.pt(), varSignImpXYSig, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_xy_flavour_run2"), track.pt(), varImpXY, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_flavour_run2"), track.pt(), varSignImpXY, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_track_pt_impact_parameter_xy_significance_flavour_run2"), track.pt(), varImpXYSig, jetflavourRun2Def, eventWeight);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_significance_flavour_run2"), track.pt(), varSignImpXYSig, jetflavourRun2Def, eventWeight);
        vecImpXY[jetflavour].push_back(varImpXY);
        vecSignImpXY[jetflavour].push_back(varSignImpXY);
        vecImpXYSig[jetflavour].push_back(varImpXYSig);
        vecSignImpXYSig[jetflavour].push_back(varSignImpXYSig);
      }
    }
    if (!fillTrackCounting)
      return;
    sort(vecImpXY[jetflavour].begin(), vecImpXY[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXY[jetflavour].begin(), vecSignImpXY[jetflavour].end(), std::greater<float>());
    sort(vecImpXYSig[jetflavour].begin(), vecImpXYSig[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXYSig[jetflavour].begin(), vecSignImpXYSig[jetflavour].end(), std::greater<float>());
    if (vecImpXY[jetflavour].size() > 0) { // N1
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2_N1"), mcdjet.pt(), vecSignImpXYSig[jetflavour][0], jetflavourRun2Def, eventWeight);
    }
    if (vecImpXY[jetflavour].size() > 1) { // N2
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2_N2"), mcdjet.pt(), vecSignImpXYSig[jetflavour][1], jetflavourRun2Def, eventWeight);
    }
    if (vecImpXY[jetflavour].size() > 2) { // N3
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_run2_N3"), mcdjet.pt(), vecSignImpXYSig[jetflavour][2], jetflavourRun2Def, eventWeight);
    }
  }

  template <typename T>
  void fillHistogramJPData(T const& jet)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    registry.fill(HIST("h2_jet_pt_JP"), jet.pt(), jet.jetProb()[0]);
    registry.fill(HIST("h2_jet_pt_neg_log_JP"), jet.pt(), -1 * std::log(jet.jetProb()[0]));
    registry.fill(HIST("h2_jet_pt_JP_N1"), jet.pt(), jet.jetProb()[1]);
    registry.fill(HIST("h2_jet_pt_neg_log_JP_N1"), jet.pt(), -1 * TMath::Log(jet.jetProb()[1]));
    registry.fill(HIST("h2_jet_pt_JP_N2"), jet.pt(), jet.jetProb()[2]);
    registry.fill(HIST("h2_jet_pt_neg_log_JP_N2"), jet.pt(), -1 * TMath::Log(jet.jetProb()[2]));
    registry.fill(HIST("h2_jet_pt_JP_N3"), jet.pt(), jet.jetProb()[3]);
    registry.fill(HIST("h2_jet_pt_neg_log_JP_N3"), jet.pt(), -1 * TMath::Log(jet.jetProb()[3]));
  }

  template <typename T>
  void fillHistogramJPMCD(T const& mcdjet, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    registry.fill(HIST("h3_jet_pt_JP_flavour"), mcdjet.pt(), mcdjet.jetProb()[0], mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_neg_log_JP_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[0]), mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_JP_N1_flavour"), mcdjet.pt(), mcdjet.jetProb()[1], mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_neg_log_JP_N1_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[1]), mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_JP_N2_flavour"), mcdjet.pt(), mcdjet.jetProb()[2], mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_neg_log_JP_N2_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[2]), mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_JP_N3_flavour"), mcdjet.pt(), mcdjet.jetProb()[3], mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_neg_log_JP_N3_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[3]), mcdjet.origin(), eventWeight);
  }

  template <typename T, typename U>
  void fillHistogramSV2ProngData(T const& jet, U const& /*prongs*/)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    registry.fill(HIST("h_2prong_nprongs"), jet.template secondaryVertices_as<U>().size());
    for (const auto& prong : jet.template secondaryVertices_as<U>()) {
      auto Lxy = prong.decayLengthXY();
      auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      auto Lxyz = prong.decayLength();
      auto Sxyz = prong.decayLength() / prong.errorDecayLength();
      registry.fill(HIST("h2_jet_pt_2prong_Lxy"), jet.pt(), Lxy);
      registry.fill(HIST("h2_jet_pt_2prong_Sxy"), jet.pt(), Sxy);
      registry.fill(HIST("h2_jet_pt_2prong_Lxyz"), jet.pt(), Lxyz);
      registry.fill(HIST("h2_jet_pt_2prong_Sxyz"), jet.pt(), Sxyz);
      registry.fill(HIST("h2_jet_pt_2prong_sigmaLxy"), jet.pt(), prong.errorDecayLengthXY());
      registry.fill(HIST("h2_jet_pt_2prong_sigmaLxyz"), jet.pt(), prong.errorDecayLength());
    }
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax);
    if (!jettaggingutilities::prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false)) {
      return;
    }
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    auto bjetCandForXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, true);
    auto maxSxyz = bjetCandForXYZ.decayLength() / bjetCandForXYZ.errorDecayLength();
    registry.fill(HIST("h2_jet_pt_2prong_Sxy_N1"), jet.pt(), maxSxy);
    registry.fill(HIST("h2_jet_pt_2prong_Sxyz_N1"), jet.pt(), maxSxyz);
  }

  template <typename T, typename U>
  void fillHistogramSV3ProngData(T const& jet, U const& /*prongs*/)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    registry.fill(HIST("h_3prong_nprongs"), jet.template secondaryVertices_as<U>().size());
    for (const auto& prong : jet.template secondaryVertices_as<U>()) {
      auto Lxy = prong.decayLengthXY();
      auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      auto Lxyz = prong.decayLength();
      auto Sxyz = prong.decayLength() / prong.errorDecayLength();
      registry.fill(HIST("h2_jet_pt_3prong_Lxy"), jet.pt(), Lxy);
      registry.fill(HIST("h2_jet_pt_3prong_Sxy"), jet.pt(), Sxy);
      registry.fill(HIST("h2_jet_pt_3prong_Lxyz"), jet.pt(), Lxyz);
      registry.fill(HIST("h2_jet_pt_3prong_Sxyz"), jet.pt(), Sxyz);
      registry.fill(HIST("h2_jet_pt_3prong_sigmaLxy"), jet.pt(), prong.errorDecayLengthXY());
      registry.fill(HIST("h2_jet_pt_3prong_sigmaLxyz"), jet.pt(), prong.errorDecayLength());
    }
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax);
    if (!jettaggingutilities::prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false)) {
      return;
    }
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    auto bjetCandForXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, true);
    auto maxSxyz = bjetCandForXYZ.decayLength() / bjetCandForXYZ.errorDecayLength();
    registry.fill(HIST("h2_jet_pt_3prong_Sxy_N1"), jet.pt(), maxSxy);
    registry.fill(HIST("h2_jet_pt_3prong_Sxyz_N1"), jet.pt(), maxSxyz);
  }

  template <typename T, typename U>
  void fillHistogramSV2ProngMCD(T const& mcdjet, U const& /*prongs*/, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    auto origin = mcdjet.origin();
    registry.fill(HIST("h2_2prong_nprongs_flavour"), mcdjet.template secondaryVertices_as<U>().size(), origin, eventWeight);
    if (mcdjet.template secondaryVertices_as<U>().size() < 1)
      return;
    for (const auto& prong : mcdjet.template secondaryVertices_as<U>()) {
      auto Lxy = prong.decayLengthXY();
      auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      auto Lxyz = prong.decayLength();
      auto Sxyz = prong.decayLength() / prong.errorDecayLength();
      registry.fill(HIST("h3_jet_pt_2prong_Lxy_flavour"), mcdjet.pt(), Lxy, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_Sxy_flavour"), mcdjet.pt(), Sxy, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_Lxyz_flavour"), mcdjet.pt(), Lxyz, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_Sxyz_flavour"), mcdjet.pt(), Sxyz, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_sigmaLxy_flavour"), mcdjet.pt(), prong.errorDecayLengthXY(), origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_sigmaLxyz_flavour"), mcdjet.pt(), prong.errorDecayLength(), origin, eventWeight);
    }
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax);
    if (!jettaggingutilities::prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false)) {
      return;
    }
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    auto massSV = bjetCand.m();
    auto bjetCandForXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, true);
    auto maxSxyz = bjetCandForXYZ.decayLength() / bjetCandForXYZ.errorDecayLength();
    registry.fill(HIST("h3_jet_pt_2prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
    registry.fill(HIST("h3_jet_pt_2prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
    registry.fill(HIST("h3_jet_pt_2prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
    if (!mcdjet.flagtaggedjetSV())
      return;
    registry.fill(HIST("h3_taggedjet_pt_2prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_2prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_2prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
  }

  template <typename T, typename U, typename V, typename W>
  void fillHistogramSV2ProngMCPMCDMatched(T const& mcdjet, U const& /*mcpjets*/, V const& /*prongs*/, W const& particlesPerColl, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    int jetflavourRun2Def = -1;
    for (auto& mcpjet : mcdjet.template matchedJetGeo_as<U>()) {
      jetflavourRun2Def = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
    }
    if (jetflavourRun2Def < 0)
      return;
    if (mcdjet.template secondaryVertices_as<V>().size() < 1)
      return;
    for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
      auto Lxy = prong.decayLengthXY();
      auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      auto Lxyz = prong.decayLength();
      auto Sxyz = prong.decayLength() / prong.errorDecayLength();
      registry.fill(HIST("h3_jet_pt_2prong_Lxy_flavour_run2"), mcdjet.pt(), Lxy, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_Sxy_flavour_run2"), mcdjet.pt(), Sxy, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_Lxyz_flavour_run2"), mcdjet.pt(), Lxyz, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_Sxyz_flavour_run2"), mcdjet.pt(), Sxyz, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_sigmaLxy_flavour_run2"), mcdjet.pt(), prong.errorDecayLengthXY(), jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_sigmaLxyz_flavour_run2"), mcdjet.pt(), prong.errorDecayLength(), jetflavourRun2Def, eventWeight);
    }
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<V>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax);
    if (!jettaggingutilities::prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false)) {
      return;
    }
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    auto massSV = bjetCand.m();
    auto bjetCandForXYZ = jettaggingutilities::jetFromProngMaxDecayLength<V>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, true);
    auto maxSxyz = bjetCandForXYZ.decayLength() / bjetCandForXYZ.errorDecayLength();
    registry.fill(HIST("h3_jet_pt_2prong_Sxy_N1_flavour_run2"), mcdjet.pt(), maxSxy, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_jet_pt_2prong_Sxyz_N1_flavour_run2"), mcdjet.pt(), maxSxyz, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_jet_pt_2prong_mass_N1_flavour_run2"), mcdjet.pt(), massSV, jetflavourRun2Def, eventWeight);
  }

  template <typename T, typename U>
  void fillHistogramSV3ProngMCD(T const& mcdjet, U const& /*prongs*/, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    auto origin = mcdjet.origin();
    registry.fill(HIST("h2_3prong_nprongs_flavour"), mcdjet.template secondaryVertices_as<U>().size(), origin);
    if (mcdjet.template secondaryVertices_as<U>().size() < 1)
      return;
    for (const auto& prong : mcdjet.template secondaryVertices_as<U>()) {
      auto Lxy = prong.decayLengthXY();
      auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      auto Lxyz = prong.decayLength();
      auto Sxyz = prong.decayLength() / prong.errorDecayLength();
      registry.fill(HIST("h3_jet_pt_3prong_Lxy_flavour"), mcdjet.pt(), Lxy, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_Sxy_flavour"), mcdjet.pt(), Sxy, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_Lxyz_flavour"), mcdjet.pt(), Lxyz, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_Sxyz_flavour"), mcdjet.pt(), Sxyz, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_sigmaLxy_flavour"), mcdjet.pt(), prong.errorDecayLengthXY(), origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_sigmaLxyz_flavour"), mcdjet.pt(), prong.errorDecayLength(), origin, eventWeight);
    }
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax);
    if (!jettaggingutilities::prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false)) {
      return;
    }
    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    auto massSV = bjetCand.m();
    auto bjetCandForXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, true);
    auto maxSxyz = bjetCandForXYZ.decayLength() / bjetCandForXYZ.errorDecayLength();
    registry.fill(HIST("h3_jet_pt_3prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
    registry.fill(HIST("h3_jet_pt_3prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
    registry.fill(HIST("h3_jet_pt_3prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
    if (!mcdjet.flagtaggedjetSV())
      return;
    registry.fill(HIST("h3_taggedjet_pt_3prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_3prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_3prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
  }

  template <typename T, typename U, typename V, typename W>
  void fillHistogramSV3ProngMCPMCDMatched(T const& mcdjet, U const& /*mcpjets*/, V const& /*prongs*/, W const& particlesPerColl, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    int jetflavourRun2Def = -1;
    for (auto& mcpjet : mcdjet.template matchedJetGeo_as<U>()) {
      jetflavourRun2Def = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
    }
    if (jetflavourRun2Def < 0)
      return;
    if (mcdjet.template secondaryVertices_as<V>().size() < 1)
      return;
    for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
      auto Lxy = prong.decayLengthXY();
      auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      auto Lxyz = prong.decayLength();
      auto Sxyz = prong.decayLength() / prong.errorDecayLength();
      registry.fill(HIST("h3_jet_pt_3prong_Lxy_flavour_run2"), mcdjet.pt(), Lxy, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_Sxy_flavour_run2"), mcdjet.pt(), Sxy, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_Lxyz_flavour_run2"), mcdjet.pt(), Lxyz, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_Sxyz_flavour_run2"), mcdjet.pt(), Sxyz, jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_sigmaLxy_flavour_run2"), mcdjet.pt(), prong.errorDecayLengthXY(), jetflavourRun2Def, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_sigmaLxyz_flavour_run2"), mcdjet.pt(), prong.errorDecayLength(), jetflavourRun2Def, eventWeight);
    }
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<V>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax);
    if (!jettaggingutilities::prongAcceptance(bjetCand, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false)) {
      return;
    }

    auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
    auto massSV = bjetCand.m();
    auto bjetCandForXYZ = jettaggingutilities::jetFromProngMaxDecayLength<V>(mcdjet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyzMax, prongIPxyMin, prongIPxyMax, true);
    auto maxSxyz = bjetCandForXYZ.decayLength() / bjetCandForXYZ.errorDecayLength();
    registry.fill(HIST("h3_jet_pt_3prong_Sxy_N1_flavour_run2"), mcdjet.pt(), maxSxy, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_jet_pt_3prong_Sxyz_N1_flavour_run2"), mcdjet.pt(), maxSxyz, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_jet_pt_3prong_mass_N1_flavour_run2"), mcdjet.pt(), massSV, jetflavourRun2Def, eventWeight);
    if (!mcdjet.flagtaggedjetSV())
      return;
    registry.fill(HIST("h3_taggedjet_pt_3prong_Sxy_N1_flavour_run2"), mcdjet.pt(), maxSxy, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_3prong_Sxyz_N1_flavour_run2"), mcdjet.pt(), maxSxyz, jetflavourRun2Def, eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_3prong_mass_N1_flavour_run2"), mcdjet.pt(), massSV, jetflavourRun2Def, eventWeight);
  }

  void processDummy(aod::Collision const&, aod::Tracks const&)
  {
  }
  PROCESS_SWITCH(JetTaggerHFQA, processDummy, "Dummy process", true);

  void processTracksDca(JetTagTracksData& jtracks)
  {
    for (auto const& jtrack : jtracks) {
      if (!jetderiveddatautilities::selectTrack(jtrack, trackSelection)) {
        continue;
      }

      float varImpXY, varImpXYSig, varImpZ, varImpZSig, varImpXYZ, varImpXYZSig;
      varImpXY = jtrack.dcaXY() * jettaggingutilities::cmTomum;
      varImpXYSig = jtrack.dcaXY() / jtrack.sigmadcaXY();
      varImpZ = jtrack.dcaZ() * jettaggingutilities::cmTomum;
      varImpZSig = jtrack.dcaZ() / jtrack.sigmadcaZ();
      float dcaXYZ = jtrack.dcaXYZ();
      float sigmadcaXYZ2 = jtrack.sigmadcaXYZ();
      varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
      varImpXYZSig = dcaXYZ / std::sqrt(sigmadcaXYZ2);

      registry.fill(HIST("h_impact_parameter_xy"), varImpXY);
      registry.fill(HIST("h_impact_parameter_xy_significance"), varImpXYSig);
      registry.fill(HIST("h_impact_parameter_z"), varImpZ);
      registry.fill(HIST("h_impact_parameter_z_significance"), varImpZSig);
      registry.fill(HIST("h_impact_parameter_xyz"), varImpXYZ);
      registry.fill(HIST("h_impact_parameter_xyz_significance"), varImpXYZSig);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processTracksDca, "Fill inclusive tracks' imformation for data", false);

  void processIPsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData> const& jets, JetTagTracksData const& jtracks)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillHistogramIPsData(jet, jtracks);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsData, "Fill impact parameter imformation for data jets", false);

  void processIPsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD> const& mcdjets, JetTagTracksMCD const& jtracks, JetParticles&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramIPsMCD(mcdjet, jtracks);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCD, "Fill impact parameter imformation for mcd jets", false);

  void processIPsMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::ChargedMCDetectorLevelJetEventWeights> const& mcdjets, JetTagTracksMCD const& jtracks, JetParticles&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillHistogramIPsMCD(mcdjet, jtracks, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCDWeighted, "Fill impact parameter imformation for mcd jets", false);

  void processIPsMCP(soa::Join<JetTableMCP, TagTableMCP> const& mcpjets, JetParticles&, JetMcCollisions const&, soa::Filtered<JetCollisionsMCD> const& collisions)
  {
    for (auto mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetParticles>(mcpjet)) {
        return;
      }
      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, mcpjet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelection)) {
          fillHistogramIPsMCP(mcpjet);
        }
      } else {
        fillHistogramIPsMCP(mcpjet);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCP, "Fill impact parameter imformation for mcp jets", false);

  void processIPsMCPWeighted(soa::Filtered<JetCollisionsMCD> const& collisions, soa::Join<JetTableMCP, TagTableMCP, aod::ChargedMCParticleLevelJetEventWeights> const& mcpjets, JetParticles&)
  {
    for (auto mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetParticles>(mcpjet)) {
        return;
      }
      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, mcpjet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelection)) {
          fillHistogramIPsMCP(mcpjet, mcpjet.eventWeight());
        } else {
          fillHistogramIPsMCP(mcpjet, mcpjet.eventWeight());
        }
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCPWeighted, "Fill impact parameter imformation for mcp jets weighted", false);

  void processIPsMCPMCDMatched(soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& mcpjets, JetTagTracksMCD const& jtracks, JetParticles& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramIPsMatched(mcdjet, mcpjets, jtracks, particlesPerColl);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCPMCDMatched, "Fill impact parameter imformation for mcp mcd matched jets", false);

  void processIPsMCPMCDMatchedWeighted(soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& mcpjets, JetTagTracksMCD const& jtracks, JetParticles& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramIPsMatched(mcdjet, mcpjets, jtracks, particlesPerColl, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCPMCDMatchedWeighted, "Fill impact parameter imformation for mcp mcd matched jets", false);

  void processJPData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD> const& jets, JetTagTracksData const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillHistogramJPData(jet);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPData, "Fill jet probability imformation for data jets", false);

  void processJPMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD> const& mcdjets, JetTagTracksMCD const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramJPMCD(mcdjet);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCD, "Fill jet probability imformation for mcd jets", false);

  void processJPMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD> const& mcdjets, JetTagTracksMCD const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramJPMCD(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCDWeighted, "Fill jet probability imformation for mcd jets", false);

  void processSV2ProngData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex2ProngIndices> const& jets, aod::DataSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillHistogramSV2ProngData(jet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngData, "Fill 2prong imformation for data jets", false);

  void processSV3ProngData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillHistogramSV3ProngData(jet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngData, "Fill 2prong imformation for data jets", false);

  void processSV2ProngMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV2ProngMCD(mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCD, "Fill 2prong imformation for mcd jets", false);

  void processSV2ProngMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV2ProngMCD(mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCDWeighted, "Fill 2prong imformation for mcd jets", false);

  void processSV2ProngMCPMCDMatched(soa::Filtered<soa::Join<JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& mcpjets, aod::MCDSecondaryVertex2Prongs const& prongs, JetParticles& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV2ProngMCPMCDMatched(mcdjet, mcpjets, prongs, particlesPerColl);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCPMCDMatched, "Fill 2prong imformation for mcd jets", false);

  void processSV2ProngMCPMCDMatchedWeighted(soa::Filtered<soa::Join<JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& mcpjets, aod::MCDSecondaryVertex2Prongs const& prongs, JetParticles& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV2ProngMCPMCDMatched(mcdjet, mcpjets, prongs, particlesPerColl, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCPMCDMatchedWeighted, "Fill 2prong imformation for mcd jets", false);

  void processSV3ProngMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCD(mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCD, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCD(mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCDWeighted, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCPMCDMatched(soa::Filtered<soa::Join<JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& mcpjets, aod::MCDSecondaryVertex3Prongs const& prongs, JetParticles& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCPMCDMatched(mcdjet, mcpjets, prongs, particlesPerColl);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCPMCDMatched, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCPMCDMatchedWeighted(soa::Filtered<soa::Join<JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& mcpjets, aod::MCDSecondaryVertex3Prongs const& prongs, JetParticles& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCPMCDMatched(mcdjet, mcpjets, prongs, particlesPerColl, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCPMCDMatchedWeighted, "Fill 3prong imformation for mcd jets", false);
};

using JetTaggerQAChargedDataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
using JetTaggerQAChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
using JetTaggerQAChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

using JetTaggerQACharged = JetTaggerHFQA<JetTaggerQAChargedDataJets, aod::ChargedJetTags, JetTaggerQAChargedMCDJets, aod::ChargedMCDetectorLevelJetEventWeights, aod::ChargedMCDetectorLevelJetTags, JetTaggerQAChargedMCPJets, aod::ChargedMCParticleLevelJetEventWeights, aod::ChargedMCParticleLevelJetTags, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

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
