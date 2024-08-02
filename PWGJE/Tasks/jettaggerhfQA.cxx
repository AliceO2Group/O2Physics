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
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> prong2sigmaLxyMax{"prong2simgaLxyMax", 0.03, "maximum sigma of decay length of 2-prong on xy plane"};
  Configurable<float> prong2SxyMin{"prong2SxyMin", 7, "minimum decay length significance of 2-prong on xy plane"};
  Configurable<float> prong2sigmaLxyzMax{"prong2sigmaLxyzMax", 0.03, "maximum sigma of decay length of 2-prong on xyz plane"};
  Configurable<float> prong2SxyzMin{"prong2SxyzMin", 7, "minimum decay length significance of 2-prong on xyz plane"};
  Configurable<float> prong3sigmaLxyMax{"prong3sigmaLxyMax", 0.03, "maximum sigma of decay length of 3-prong on xy plane"};
  Configurable<float> prong3SxyMin{"prong3SxyMin", 7, "minimum decay length significance of 3-prong on xy plane"};
  Configurable<float> prong3sigmaLxyzMax{"prong3sigmaLxyzMax", 0.03, "maximum sigma of decay length of 3-prong on xyz plane"};
  Configurable<float> prong3SxyzMin{"prong3SxyzMin", 7, "minimum decay length significance of 3-prong on xyz plane"};
  Configurable<int> numFlavourSpecies{"numFlavourSpecies", 6, "number of jet flavour species"};
  Configurable<int> numOrder{"numOrder", 6, "number of ordering"};

  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  // Binning
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, 0., 100.}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {800, -400.5f, 400.5f}, ""};
  ConfigurableAxis binSigmaImpactParameterXY{"binImpactSigmaParameterXY", {800, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXYSignificance{"binImpactParameterXYSignificance", {800, -40.5f, 40.5f}, ""};
  ConfigurableAxis binImpactParameterZ{"binImpactParameterZ", {800, -400.5f, 400.5f}, ""};
  ConfigurableAxis binImpactParameterZSignificance{"binImpactParameterZSignificance", {800, -40.5f, 40.5f}, ""};
  ConfigurableAxis binImpactParameterXYZ{"binImpactParameterXYZ", {2000, -1000.5f, 1000.5f}, ""};
  ConfigurableAxis binImpactParameterXYZSignificance{"binImpactParameterXYZSignificance", {2000, -100.5f, 100.5f}, ""};
  ConfigurableAxis binNumOrder{"binNumOrder", {6, 0.5, 6.5}, ""};
  ConfigurableAxis binJetProbability{"binJetProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetProbabilityLog{"binJetProbabilityLog", {100, 0.f, 10.f}, ""};
  ConfigurableAxis binNprongs{"binNprongs", {100, 0., 100.}, ""};
  ConfigurableAxis binLxy{"binLxy", {200, 0, 20.f}, ""};
  ConfigurableAxis binSxy{"binSxy", {200, 0, 200.f}, ""};
  ConfigurableAxis binLxyz{"binLxyz", {200, 0, 20.f}, ""};
  ConfigurableAxis binSxyz{"binSxyz", {200, 0, 200.f}, ""};
  ConfigurableAxis binSigmaLxy{"binSigmaLxy", {100, 0., 0.1}, ""};
  ConfigurableAxis binSigmaLxyz{"binSigmaLxyz", {100, 0., 0.1}, ""};

  int numberOfJetFlavourSpecies = 6;
  int trackSelection = -1;
  float maxSigmaLxy2Prong = 0.;
  float minSxy2Prong = 0.;
  float maxSigmaLxyz2Prong = 0.;
  float minSxyz2Prong = 0.;
  float maxSigmaLxy3Prong = 0.;
  float minSxy3Prong = 0.;
  float maxSigmaLxyz3Prong = 0.;
  float minSxyz3Prong = 0.;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // Axis
    AxisSpec jetFlavourAxis = {binJetFlavour, "Jet flavour"};
    AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T, jet}"};
    AxisSpec etaAxis = {binEta, "#eta"};
    AxisSpec phiAxis = {binPhi, "#phi"};
    AxisSpec ntracksAxis = {binNtracks, "#it{N}_{tracks}"};
    AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{track}"};
    AxisSpec impactParameterXYAxis = {binImpactParameterXY, "IP_{XY} [#mum]"};
    AxisSpec sigmaImpactParameterXYAxis = {binSigmaImpactParameterXY, "#sigma_{XY} [#mum]"};
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
    AxisSpec sigmaLxyAxis = {binSigmaLxy, "#sigma_{L_{XY}} [cm]"};
    AxisSpec sigmaLxyzAxis = {binSigmaLxyz, "#sigma_{L_{XYZ}} [cm]"};

    numberOfJetFlavourSpecies = static_cast<int>(numFlavourSpecies);
    maxSigmaLxy2Prong = static_cast<float>(prong2sigmaLxyMax);
    minSxy2Prong = static_cast<float>(prong2SxyMin);
    maxSigmaLxyz2Prong = static_cast<float>(prong2sigmaLxyzMax);
    minSxyz2Prong = static_cast<float>(prong2SxyzMin);
    maxSigmaLxy3Prong = static_cast<float>(prong3sigmaLxyMax);
    minSxy3Prong = static_cast<float>(prong3SxyMin);
    maxSigmaLxyz3Prong = static_cast<float>(prong3sigmaLxyzMax);
    minSxyz3Prong = static_cast<float>(prong3SxyzMin);

    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
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
    if (doprocessIPsMCD) {
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
    if (doprocessJPMCD) {
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
      registry.add("h2_jet_pt_2prong_Sxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_2prong_Lxyz", "", {HistType::kTH2F, {{jetPtAxis}, {LxyzAxis}}});
      registry.add("h2_jet_pt_2prong_Sxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_2prong_Sxy_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_2prong_Sxyz_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_2prong_Sxy_sigmaLxy", "", {HistType::kTH2F, {{SxyAxis}, {sigmaLxyAxis}}});
      registry.add("h2_2prong_Sxyz_sigmaLxyz", "", {HistType::kTH2F, {{SxyzAxis}, {sigmaLxyzAxis}}});
      registry.add("h2_jet_pt_2prong_sigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyAxis}}});
      registry.add("h2_jet_pt_2prong_sigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyzAxis}}});
      registry.add("h2_jet_pt_2prong_Sxy_cutSxyAndsigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_2prong_Sxyz_cutSxyzAndsigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_2prong_Sxy_N1_cutSxyAndsigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_2prong_Sxyz_N1_cutSxyzAndsigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
    }
    if (doprocessSV3ProngData) {
      registry.add("h_3prong_nprongs", "", {HistType::kTH1F, {{nprongsAxis}}});
      registry.add("h2_jet_pt_3prong_Lxy", "", {HistType::kTH2F, {{jetPtAxis}, {LxyAxis}}});
      registry.add("h2_jet_pt_3prong_Sxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_3prong_Lxyz", "", {HistType::kTH2F, {{jetPtAxis}, {LxyzAxis}}});
      registry.add("h2_jet_pt_3prong_Sxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_3prong_Sxy_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_3prong_Sxyz_N1", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_3prong_Sxy_sigmaLxy", "", {HistType::kTH2F, {{SxyAxis}, {sigmaLxyAxis}}});
      registry.add("h2_3prong_Sxyz_sigmaLxyz", "", {HistType::kTH2F, {{SxyzAxis}, {sigmaLxyzAxis}}});
      registry.add("h2_jet_pt_3prong_sigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyAxis}}});
      registry.add("h2_jet_pt_3prong_sigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {sigmaLxyzAxis}}});
      registry.add("h2_jet_pt_3prong_Sxy_cutSxyAndsigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_3prong_Sxyz_cutSxyzAndsigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
      registry.add("h2_jet_pt_3prong_Sxy_N1_cutSxyAndsigmaLxy", "", {HistType::kTH2F, {{jetPtAxis}, {SxyAxis}}});
      registry.add("h2_jet_pt_3prong_Sxyz_N1_cutSxyzAndsigmaLxyz", "", {HistType::kTH2F, {{jetPtAxis}, {SxyzAxis}}});
    }
    if (doprocessSV2ProngMCD) {
      registry.add("h2_2prong_nprongs_flavour", "", {HistType::kTH2F, {{nprongsAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Lxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Lxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_2prong_Sxy_sigmaLxy_flavour", "", {HistType::kTH3F, {{SxyAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_2prong_Sxyz_sigmaLxyz_flavour", "", {HistType::kTH3F, {{SxyzAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_sigmaLxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_sigmaLxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_flavour_cutSxyAndsigmaLxy", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_flavour_cutSxyzAndsigmaLxyz", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxy_N1_flavour_cutSxyAndsigmaLxy", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_2prong_Sxyz_N1_flavour_cutSxyzAndsigmaLxyz", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
    }
    if (doprocessSV3ProngMCD) {
      registry.add("h2_3prong_nprongs_flavour", "", {HistType::kTH2F, {{nprongsAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Lxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Lxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {LxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_3prong_Sxy_sigmaLxy_flavour", "", {HistType::kTH3F, {{SxyAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_3prong_Sxyz_sigmaLxyz_flavour", "", {HistType::kTH3F, {{SxyzAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_sigmaLxy_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_sigmaLxyz_flavour", "", {HistType::kTH3F, {{jetPtAxis}, {sigmaLxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_flavour_cutSxyAndsigmaLxy", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_flavour_cutSxyzAndsigmaLxyz", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxy_N1_flavour_cutSxyAndsigmaLxy", "", {HistType::kTH3F, {{jetPtAxis}, {SxyAxis}, {jetFlavourAxis}}});
      registry.add("h3_jet_pt_3prong_Sxyz_N1_flavour_cutSxyzAndsigmaLxyz", "", {HistType::kTH3F, {{jetPtAxis}, {SxyzAxis}, {jetFlavourAxis}}});
    }
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);

  using JetTagTracksData = soa::Join<JetTracks, aod::JTrackPIs, aod::JTracksTag>;
  using JetTagTracksMCD = soa::Join<JetTracksMCD, aod::JTrackPIs, aod::JTracksTag>;
  using OriTracksData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection>;
  using OriTracksMCD = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::McTrackLabels>;

  std::function<bool(const std::vector<float>&, const std::vector<float>&)> sortImp =
    [](const std::vector<float>& a, const std::vector<float>& b) {
      return a[0] > b[0];
    };

  template <typename T>
  bool trackAcceptance(T const& track)
  {
    if (track.pt() < trackPtMin || track.pt() > trackPtMax)
      return false;

    return true;
  }

  bool prongAcceptance(float sigmaDecayLength, float decayLengthSig, float maxSigmaDecayLength, float minDecayLengthSig)
  {
    if ((sigmaDecayLength < maxSigmaDecayLength) && (decayLengthSig > minDecayLengthSig))
      return true;

    return false;
  }

  template <typename V, typename JetType>
  std::tuple<float, float> getMaxSxyForJet(const JetType& mcdjet)
  {
    float maxSxy = 0;
    float correspondingErrorDecayLengthXY = 0;

    for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
      float Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
      if (maxSxy < Sxy) {
        maxSxy = Sxy;
        correspondingErrorDecayLengthXY = prong.errorDecayLengthXY();
      }
    }
    return std::make_tuple(maxSxy, correspondingErrorDecayLengthXY);
  }

  template <typename V, typename JetType>
  std::tuple<float, float> getMaxSxyzForJet(const JetType& mcdjet)
  {
    float maxSxyz = 0;
    float correspondingErrorDecayLength = 0;

    for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
      float Sxyz = prong.decayLength() / prong.errorDecayLength();
      if (maxSxyz < Sxyz) {
        maxSxyz = Sxyz;
        correspondingErrorDecayLength = prong.errorDecayLength();
      }
    }
    return std::make_tuple(maxSxyz, correspondingErrorDecayLength);
  }

  template <typename T, typename U, typename V, typename W>
  void fillHistogramIPsData(T const& collision, U const& jets, V const& /*jtracks*/, W const& /*tracks*/)
  {
    for (auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      std::vector<std::vector<float>> vecSignImpXYSig, vecSignImpZSig, vecSignImpXYZSig;
      for (auto& jtrack : jet.template tracks_as<V>()) {
        auto track = jtrack.template track_as<W>();
        if (!trackAcceptance(track))
          continue;

        // General parameters
        registry.fill(HIST("h3_jet_pt_track_pt_track_eta"), jet.pt(), track.pt(), track.eta());
        registry.fill(HIST("h3_jet_pt_track_pt_track_phi"), jet.pt(), track.pt(), track.phi());
        int geoSign = jettaggingutilities::getGeoSign(collision, jet, track);
        if (fillIPxy) {
          float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig;
          varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
          varSignImpXY = geoSign * std::abs(track.dcaXY()) * jettaggingutilities::cmTomum;
          varImpXYSig = track.dcaXY() / std::sqrt(track.sigmaDcaXY2());
          varSignImpXYSig = geoSign * std::abs(track.dcaXY()) / std::sqrt(track.sigmaDcaXY2());
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
          varImpZSig = track.dcaZ() / std::sqrt(track.sigmaDcaZ2());
          varSignImpZSig = geoSign * std::abs(track.dcaZ()) / std::sqrt(track.sigmaDcaZ2());
          registry.fill(HIST("h2_jet_pt_impact_parameter_z"), jet.pt(), varImpZ);
          registry.fill(HIST("h2_jet_pt_sign_impact_parameter_z"), jet.pt(), varSignImpZ);
          registry.fill(HIST("h2_jet_pt_impact_parameter_z_significance"), jet.pt(), varImpZSig);
          registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_z_significance"), jet.pt(), track.pt(), varSignImpZSig);
          vecSignImpZSig.push_back({varSignImpZSig, track.pt()});
        }
        if (fillIPxyz) {
          float varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
          float dcaXYZ = jtrack.dcaXYZ();
          float sigmaDcaXYZ2 = jtrack.sigmaDcaXYZ2();
          varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
          varSignImpXYZ = geoSign * std::abs(dcaXYZ) * jettaggingutilities::cmTomum;
          varImpXYZSig = dcaXYZ / std::sqrt(sigmaDcaXYZ2);
          varSignImpXYZSig = geoSign * std::abs(dcaXYZ) / std::sqrt(sigmaDcaXYZ2);
          registry.fill(HIST("h2_jet_pt_impact_parameter_xyz"), jet.pt(), varImpXYZ);
          registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xyz"), jet.pt(), varSignImpXYZ);
          registry.fill(HIST("h2_jet_pt_impact_parameter_xyz_significance"), jet.pt(), varImpXYZSig);
          registry.fill(HIST("h3_jet_pt_track_pt_sign_impact_parameter_xyz_significance"), jet.pt(), track.pt(), varSignImpXYZSig);
          vecSignImpXYZSig.push_back({varSignImpXYZSig, track.pt()});
        }
      }

      if (!fillTrackCounting)
        continue;
      if (fillIPxy)
        std::sort(vecSignImpXYSig.begin(), vecSignImpXYSig.end(), sortImp);
      if (fillIPz)
        std::sort(vecSignImpZSig.begin(), vecSignImpZSig.end(), sortImp);
      if (fillIPxyz)
        std::sort(vecSignImpXYZSig.begin(), vecSignImpXYZSig.end(), sortImp);
      if (fillIPxy && vecSignImpXYSig.empty())
        continue;
      if (fillIPz && vecSignImpZSig.empty())
        continue;
      if (fillIPxyz && vecSignImpXYZSig.empty())
        continue;
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
  }

  template <typename T, typename U, typename V, typename W>
  void fillHistogramIPsMCD(T const& collision, U const& mcdjets, V const& /*jtracks*/, W const& /*tracks*/)
  {
    for (auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      std::vector<float> vecImpXY[numberOfJetFlavourSpecies], vecSignImpXY[numberOfJetFlavourSpecies], vecImpXYSig[numberOfJetFlavourSpecies], vecSignImpXYSig[numberOfJetFlavourSpecies];
      std::vector<float> vecImpZ[numberOfJetFlavourSpecies], vecSignImpZ[numberOfJetFlavourSpecies], vecImpZSig[numberOfJetFlavourSpecies], vecSignImpZSig[numberOfJetFlavourSpecies];
      std::vector<float> vecImpXYZ[numberOfJetFlavourSpecies], vecSignImpXYZ[numberOfJetFlavourSpecies], vecImpXYZSig[numberOfJetFlavourSpecies], vecSignImpXYZSig[numberOfJetFlavourSpecies];
      std::vector<std::vector<float>> vecSignImpXYSigTC, vecSignImpZSigTC, vecSignImpXYZSigTC;
      int jetflavour = mcdjet.origin();
      if (jetflavour == JetTaggingSpecies::none) {
        LOGF(debug, "NOT DEFINE JET FLAVOR");
      }
      registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), jetflavour);
      registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), jetflavour);
      registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), jetflavour);
      for (auto& jtrack : mcdjet.template tracks_as<V>()) {
        auto track = jtrack.template track_as<W>();
        if (!trackAcceptance(track))
          continue;
        // General parameters
        registry.fill(HIST("h3_jet_pt_track_pt_flavour"), mcdjet.pt(), track.pt(), jetflavour);
        registry.fill(HIST("h3_jet_pt_track_eta_flavour"), mcdjet.pt(), track.eta(), jetflavour);
        registry.fill(HIST("h3_jet_pt_track_phi_flavour"), mcdjet.pt(), track.phi(), jetflavour);
        int geoSign = jettaggingutilities::getGeoSign(collision, mcdjet, track);
        if (fillIPxy) {
          float varImpXY, varSignImpXY, varImpXYSig, varSignImpXYSig;
          varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
          float varSigmaImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
          varSignImpXY = geoSign * std::abs(track.dcaXY()) * jettaggingutilities::cmTomum;
          varImpXYSig = std::abs(track.dcaXY()) / std::sqrt(track.sigmaDcaXY2());
          varSignImpXYSig = geoSign * std::abs(track.dcaXY()) / std::sqrt(track.sigmaDcaXY2());
          registry.fill(HIST("h3_jet_pt_impact_parameter_xy_flavour"), mcdjet.pt(), varImpXY, jetflavour);
          registry.fill(HIST("h3_jet_pt_sigma_impact_parameter_xy_flavour"), mcdjet.pt(), varSigmaImpXY, jetflavour);
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_flavour"), mcdjet.pt(), varSignImpXY, jetflavour);
          registry.fill(HIST("h3_jet_pt_impact_parameter_xy_significance_flavour"), mcdjet.pt(), varImpXYSig, jetflavour);
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour"), mcdjet.pt(), varSignImpXYSig, jetflavour);
          registry.fill(HIST("h3_track_pt_impact_parameter_xy_flavour"), jtrack.pt(), varImpXY, jetflavour);
          registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_flavour"), jtrack.pt(), varSignImpXY, jetflavour);
          registry.fill(HIST("h3_track_pt_impact_parameter_xy_significance_flavour"), jtrack.pt(), varImpXYSig, jetflavour);
          registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_significance_flavour"), jtrack.pt(), varSignImpXYSig, jetflavour);
          vecImpXY[jetflavour].push_back(varImpXY);
          vecSignImpXY[jetflavour].push_back(varSignImpXY);
          vecImpXYSig[jetflavour].push_back(varImpXYSig);
          vecSignImpXYSig[jetflavour].push_back(varSignImpXYSig);
          vecSignImpXYSigTC.push_back({varSignImpXYSig, jtrack.pt()});
        }
        if (fillIPz) {
          float varImpZ, varSignImpZ, varImpZSig, varSignImpZSig;
          varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
          varSignImpZ = geoSign * std::abs(track.dcaZ()) * jettaggingutilities::cmTomum;
          varImpZSig = track.dcaZ() / std::sqrt(track.sigmaDcaZ2());
          varSignImpZSig = geoSign * std::abs(track.dcaZ()) / std::sqrt(track.sigmaDcaZ2());
          registry.fill(HIST("h3_jet_pt_impact_parameter_z_flavour"), mcdjet.pt(), varImpZ, jetflavour);
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_flavour"), mcdjet.pt(), varSignImpZ, jetflavour);
          registry.fill(HIST("h3_jet_pt_impact_parameter_z_significance_flavour"), mcdjet.pt(), varImpZSig, jetflavour);
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour"), mcdjet.pt(), varSignImpZSig, jetflavour);
          registry.fill(HIST("h3_track_pt_impact_parameter_z_flavour"), jtrack.pt(), varImpZ, jetflavour);
          registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_flavour"), jtrack.pt(), varSignImpZ, jetflavour);
          registry.fill(HIST("h3_track_pt_impact_parameter_z_significance_flavour"), jtrack.pt(), varImpZSig, jetflavour);
          registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_significance_flavour"), jtrack.pt(), varSignImpZSig, jetflavour);
          vecImpZ[jetflavour].push_back(varImpZ);
          vecSignImpZ[jetflavour].push_back(varSignImpZ);
          vecImpZSig[jetflavour].push_back(varImpZSig);
          vecSignImpZSig[jetflavour].push_back(varSignImpZSig);
          vecSignImpZSigTC.push_back({varSignImpZSig, jtrack.pt()});
        }
        if (fillIPxyz) {
          float varImpXYZ, varSignImpXYZ, varImpXYZSig, varSignImpXYZSig;
          float dcaXYZ = jtrack.dcaXYZ();
          float sigmaDcaXYZ2 = jtrack.sigmaDcaXYZ2();
          varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
          varSignImpXYZ = geoSign * std::abs(dcaXYZ) * jettaggingutilities::cmTomum;
          varImpXYZSig = dcaXYZ / std::sqrt(sigmaDcaXYZ2);
          varSignImpXYZSig = geoSign * std::abs(dcaXYZ) / std::sqrt(sigmaDcaXYZ2);
          registry.fill(HIST("h3_jet_pt_impact_parameter_xyz_flavour"), mcdjet.pt(), varImpXYZ, jetflavour);
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_flavour"), mcdjet.pt(), varSignImpXYZ, jetflavour);
          registry.fill(HIST("h3_jet_pt_impact_parameter_xyz_significance_flavour"), mcdjet.pt(), varImpXYZSig, jetflavour);
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour"), mcdjet.pt(), varSignImpXYZSig, jetflavour);
          registry.fill(HIST("h3_track_pt_impact_parameter_xyz_flavour"), jtrack.pt(), varImpXYZ, jetflavour);
          registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_flavour"), jtrack.pt(), varSignImpXYZ, jetflavour);
          registry.fill(HIST("h3_track_pt_impact_parameter_xyz_significance_flavour"), jtrack.pt(), varImpXYZSig, jetflavour);
          registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_significance_flavour"), jtrack.pt(), varSignImpXYZSig, jetflavour);
          vecImpXYZ[jetflavour].push_back(varImpXYZ);
          vecSignImpXYZ[jetflavour].push_back(varSignImpXYZ);
          vecImpXYZSig[jetflavour].push_back(varImpXYZSig);
          vecSignImpXYZSig[jetflavour].push_back(varSignImpXYZSig);
          vecSignImpXYZSigTC.push_back({varSignImpXYZSig, jtrack.pt()});
        }
      }

      if (!fillTrackCounting)
        continue;
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
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1"), mcdjet.pt(), vecSignImpXYSig[jetflavour][0], jetflavour);
        if (fillIPz)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1"), mcdjet.pt(), vecSignImpZSig[jetflavour][0], jetflavour);
        if (fillIPxyz)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][0], jetflavour);
      }
      if (vecImpXY[jetflavour].size() > 1) { // N2
        if (fillIPxy)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2"), mcdjet.pt(), vecSignImpXYSig[jetflavour][1], jetflavour);
        if (fillIPz)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2"), mcdjet.pt(), vecSignImpZSig[jetflavour][1], jetflavour);
        if (fillIPxyz)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][1], jetflavour);
      }
      if (vecImpXY[jetflavour].size() > 2) { // N3
        if (fillIPxy)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N3"), mcdjet.pt(), vecSignImpXYSig[jetflavour][2], jetflavour);
        if (fillIPz)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N3"), mcdjet.pt(), vecSignImpZSig[jetflavour][2], jetflavour);
        if (fillIPxyz)
          registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N3"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][2], jetflavour);
      }

      std::sort(vecSignImpXYSigTC.begin(), vecSignImpXYSigTC.end(), sortImp);
      std::sort(vecSignImpZSigTC.begin(), vecSignImpZSigTC.end(), sortImp);
      std::sort(vecSignImpXYZSigTC.begin(), vecSignImpXYZSigTC.end(), sortImp);

      if (vecSignImpXYSigTC.empty())
        continue;
      for (int order = 1; order <= numOrder; order++) {
        if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSigTC.size()) {
          registry.fill(HIST("h3_sign_impact_parameter_xy_significance_tc_flavour"), vecSignImpXYSigTC[order - 1][0], order, jetflavour);
        }
      }
      if (vecSignImpZSigTC.empty())
        continue;
      for (int order = 1; order <= numOrder; order++) {
        if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpZSigTC.size()) {
          registry.fill(HIST("h3_sign_impact_parameter_z_significance_tc_flavour"), vecSignImpZSigTC[order - 1][0], order, jetflavour);
        }
      }

      if (vecSignImpXYZSigTC.empty())
        continue;
      for (int order = 1; order <= numOrder; order++) {
        if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYZSigTC.size()) {
          registry.fill(HIST("h3_sign_impact_parameter_xyz_significance_tc_flavour"), vecSignImpXYZSigTC[order - 1][0], order, jetflavour);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramJPData(T const& /*collision*/, U const& jets)
  {
    for (auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
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
  }

  template <typename T, typename U>
  void fillHistogramJPMCD(T const& /*collision*/, U const& mcdjets)
  {
    for (auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      registry.fill(HIST("h3_jet_pt_JP_flavour"), mcdjet.pt(), mcdjet.jetProb()[0], mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_neg_log_JP_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[0]), mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_JP_N1_flavour"), mcdjet.pt(), mcdjet.jetProb()[1], mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_neg_log_JP_N1_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[1]), mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_JP_N2_flavour"), mcdjet.pt(), mcdjet.jetProb()[2], mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_neg_log_JP_N2_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[2]), mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_JP_N3_flavour"), mcdjet.pt(), mcdjet.jetProb()[3], mcdjet.origin());
      registry.fill(HIST("h3_jet_pt_neg_log_JP_N3_flavour"), mcdjet.pt(), -1 * TMath::Log(mcdjet.jetProb()[3]), mcdjet.origin());
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramSV2ProngData(T const& /*collision*/, U const& jets, V const& /*prongs*/)
  {
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      auto [maxSxy, sigmaLxy] = getMaxSxyForJet<V>(jet);
      auto [maxSxyz, sigmaLxyz] = getMaxSxyzForJet<V>(jet);
      registry.fill(HIST("h_2prong_nprongs"), jet.template secondaryVertices_as<V>().size());
      for (const auto& prong : jet.template secondaryVertices_as<V>()) {
        auto Lxy = prong.decayLengthXY();
        auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
        auto Lxyz = prong.decayLength();
        auto Sxyz = prong.decayLength() / prong.errorDecayLength();
        registry.fill(HIST("h2_jet_pt_2prong_Lxy"), jet.pt(), Lxy);
        registry.fill(HIST("h2_jet_pt_2prong_Sxy"), jet.pt(), Sxy);
        registry.fill(HIST("h2_jet_pt_2prong_Lxyz"), jet.pt(), Lxyz);
        registry.fill(HIST("h2_jet_pt_2prong_Sxyz"), jet.pt(), Sxyz);
        registry.fill(HIST("h2_2prong_Sxy_sigmaLxy"), Sxy, prong.errorDecayLengthXY());
        registry.fill(HIST("h2_2prong_Sxyz_sigmaLxyz"), Sxyz, prong.errorDecayLength());
        registry.fill(HIST("h2_jet_pt_2prong_sigmaLxy"), jet.pt(), prong.errorDecayLengthXY());
        registry.fill(HIST("h2_jet_pt_2prong_sigmaLxyz"), jet.pt(), prong.errorDecayLength());
        if (prongAcceptance(prong.errorDecayLengthXY(), Sxy, maxSigmaLxy2Prong, minSxy2Prong)) {
          registry.fill(HIST("h2_jet_pt_2prong_Sxy_cutSxyAndsigmaLxy"), jet.pt(), Sxy);
        }
        if (prongAcceptance(prong.errorDecayLength(), Sxyz, maxSigmaLxyz2Prong, minSxyz2Prong)) {
          registry.fill(HIST("h2_jet_pt_2prong_Sxyz_cutSxyzAndsigmaLxyz"), jet.pt(), Sxyz);
        }
      }
      registry.fill(HIST("h2_jet_pt_2prong_Sxy_N1"), jet.pt(), maxSxy);
      registry.fill(HIST("h2_jet_pt_2prong_Sxyz_N1"), jet.pt(), maxSxyz);
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy2Prong, minSxy2Prong)) {
        registry.fill(HIST("h2_jet_pt_2prong_Sxy_N1_cutSxyAndsigmaLxy"), jet.pt(), maxSxy);
      }
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy2Prong, minSxy2Prong)) {
        registry.fill(HIST("h2_jet_pt_2prong_Sxyz_N1_cutSxyzAndsigmaLxyz"), jet.pt(), maxSxyz);
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramSV3ProngData(T const& /*collision*/, U const& jets, V const& /*prongs*/)
  {
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      auto [maxSxy, sigmaLxy] = getMaxSxyForJet<V>(jet);
      auto [maxSxyz, sigmaLxyz] = getMaxSxyzForJet<V>(jet);
      registry.fill(HIST("h_3prong_nprongs"), jet.template secondaryVertices_as<V>().size());
      for (const auto& prong : jet.template secondaryVertices_as<V>()) {
        auto Lxy = prong.decayLengthXY();
        auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
        auto Lxyz = prong.decayLength();
        auto Sxyz = prong.decayLength() / prong.errorDecayLength();
        registry.fill(HIST("h2_jet_pt_3prong_Lxy"), jet.pt(), Lxy);
        registry.fill(HIST("h2_jet_pt_3prong_Sxy"), jet.pt(), Sxy);
        registry.fill(HIST("h2_jet_pt_3prong_Lxyz"), jet.pt(), Lxyz);
        registry.fill(HIST("h2_jet_pt_3prong_Sxyz"), jet.pt(), Sxyz);
        registry.fill(HIST("h2_3prong_Sxy_sigmaLxy"), Sxy, prong.errorDecayLengthXY());
        registry.fill(HIST("h2_3prong_Sxyz_sigmaLxyz"), Sxyz, prong.errorDecayLength());
        registry.fill(HIST("h2_jet_pt_3prong_sigmaLxy"), jet.pt(), prong.errorDecayLengthXY());
        registry.fill(HIST("h2_jet_pt_3prong_sigmaLxyz"), jet.pt(), prong.errorDecayLength());
        if (prongAcceptance(prong.errorDecayLengthXY(), Sxy, maxSigmaLxy3Prong, minSxy3Prong)) {
          registry.fill(HIST("h2_jet_pt_3prong_Sxy_cutSxyAndsigmaLxy"), jet.pt(), Sxy);
        }
        if (prongAcceptance(prong.errorDecayLength(), Sxyz, maxSigmaLxyz3Prong, minSxyz3Prong)) {
          registry.fill(HIST("h2_jet_pt_3prong_Sxyz_cutSxyzAndsigmaLxyz"), jet.pt(), Sxyz);
        }
      }
      registry.fill(HIST("h2_jet_pt_3prong_Sxy_N1"), jet.pt(), maxSxy);
      registry.fill(HIST("h2_jet_pt_3prong_Sxyz_N1"), jet.pt(), maxSxyz);
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy3Prong, minSxy3Prong)) {
        registry.fill(HIST("h2_jet_pt_3prong_Sxy_N1_cutSxyAndsigmaLxy"), jet.pt(), maxSxy);
      }
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy3Prong, minSxy3Prong)) {
        registry.fill(HIST("h2_jet_pt_3prong_Sxyz_N1_cutSxyzAndsigmaLxyz"), jet.pt(), maxSxyz);
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramSV2ProngMCD(T const& /*collision*/, U const& mcdjets, V const& /*prongs*/)
  {
    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      auto origin = mcdjet.origin();
      auto [maxSxy, sigmaLxy] = getMaxSxyForJet<V>(mcdjet);
      auto [maxSxyz, sigmaLxyz] = getMaxSxyzForJet<V>(mcdjet);
      registry.fill(HIST("h2_2prong_nprongs_flavour"), mcdjet.template secondaryVertices_as<V>().size(), origin);
      for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
        auto Lxy = prong.decayLengthXY();
        auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
        auto Lxyz = prong.decayLength();
        auto Sxyz = prong.decayLength() / prong.errorDecayLength();
        registry.fill(HIST("h3_jet_pt_2prong_Lxy_flavour"), mcdjet.pt(), Lxy, origin);
        registry.fill(HIST("h3_jet_pt_2prong_Sxy_flavour"), mcdjet.pt(), Sxy, origin);
        registry.fill(HIST("h3_jet_pt_2prong_Lxyz_flavour"), mcdjet.pt(), Lxyz, origin);
        registry.fill(HIST("h3_jet_pt_2prong_Sxyz_flavour"), mcdjet.pt(), Sxyz, origin);
        registry.fill(HIST("h3_2prong_Sxy_sigmaLxy_flavour"), Sxy, prong.errorDecayLengthXY(), origin);
        registry.fill(HIST("h3_2prong_Sxyz_sigmaLxyz_flavour"), Sxyz, prong.errorDecayLength(), origin);
        registry.fill(HIST("h3_jet_pt_2prong_sigmaLxy_flavour"), mcdjet.pt(), prong.errorDecayLengthXY(), origin);
        registry.fill(HIST("h3_jet_pt_2prong_sigmaLxyz_flavour"), mcdjet.pt(), prong.errorDecayLength(), origin);
        if (prongAcceptance(prong.errorDecayLengthXY(), Sxy, maxSigmaLxy2Prong, minSxy2Prong)) {
          registry.fill(HIST("h3_jet_pt_2prong_Sxy_flavour_cutSxyAndsigmaLxy"), mcdjet.pt(), Sxy, origin);
        }
        if (prongAcceptance(prong.errorDecayLength(), Sxyz, maxSigmaLxyz2Prong, minSxyz2Prong)) {
          registry.fill(HIST("h3_jet_pt_2prong_Sxyz_flavour_cutSxyzAndsigmaLxyz"), mcdjet.pt(), Sxyz, origin);
        }
      }
      registry.fill(HIST("h3_jet_pt_2prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin);
      registry.fill(HIST("h3_jet_pt_2prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin);
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy2Prong, minSxy2Prong)) {
        registry.fill(HIST("h3_jet_pt_2prong_Sxy_N1_flavour_cutSxyAndsigmaLxy"), mcdjet.pt(), maxSxy, origin);
      }
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy2Prong, minSxy2Prong)) {
        registry.fill(HIST("h3_jet_pt_2prong_Sxyz_N1_flavour_cutSxyzAndsigmaLxyz"), mcdjet.pt(), maxSxyz, origin);
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramSV3ProngMCD(T const& /*collision*/, U const& mcdjets, V const& /*prongs*/)
  {
    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      auto origin = mcdjet.origin();
      auto [maxSxy, sigmaLxy] = getMaxSxyForJet<V>(mcdjet);
      auto [maxSxyz, sigmaLxyz] = getMaxSxyzForJet<V>(mcdjet);
      registry.fill(HIST("h2_3prong_nprongs_flavour"), mcdjet.template secondaryVertices_as<V>().size(), origin);
      for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
        auto Lxy = prong.decayLengthXY();
        auto Sxy = prong.decayLengthXY() / prong.errorDecayLengthXY();
        auto Lxyz = prong.decayLength();
        auto Sxyz = prong.decayLength() / prong.errorDecayLength();
        registry.fill(HIST("h3_jet_pt_3prong_Lxy_flavour"), mcdjet.pt(), Lxy, origin);
        registry.fill(HIST("h3_jet_pt_3prong_Sxy_flavour"), mcdjet.pt(), Sxy, origin);
        registry.fill(HIST("h3_jet_pt_3prong_Lxyz_flavour"), mcdjet.pt(), Lxyz, origin);
        registry.fill(HIST("h3_jet_pt_3prong_Sxyz_flavour"), mcdjet.pt(), Sxyz, origin);
        registry.fill(HIST("h3_3prong_Sxy_sigmaLxy_flavour"), Sxy, prong.errorDecayLengthXY(), origin);
        registry.fill(HIST("h3_3prong_Sxyz_sigmaLxyz_flavour"), Sxyz, prong.errorDecayLength(), origin);
        registry.fill(HIST("h3_jet_pt_3prong_sigmaLxy_flavour"), mcdjet.pt(), prong.errorDecayLengthXY(), origin);
        registry.fill(HIST("h3_jet_pt_3prong_sigmaLxyz_flavour"), mcdjet.pt(), prong.errorDecayLength(), origin);
        if (prongAcceptance(prong.errorDecayLengthXY(), Sxy, maxSigmaLxy3Prong, minSxy3Prong)) {
          registry.fill(HIST("h3_jet_pt_3prong_Sxy_flavour_cutSxyAndsigmaLxy"), mcdjet.pt(), Sxy, origin);
        }
        if (prongAcceptance(prong.errorDecayLength(), Sxyz, maxSigmaLxyz3Prong, minSxyz3Prong)) {
          registry.fill(HIST("h3_jet_pt_3prong_Sxyz_flavour_cutSxyzAndsigmaLxyz"), mcdjet.pt(), Sxyz, origin);
        }
      }
      registry.fill(HIST("h3_jet_pt_3prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin);
      registry.fill(HIST("h3_jet_pt_3prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin);
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy3Prong, minSxy3Prong)) {
        registry.fill(HIST("h3_jet_pt_3prong_Sxy_N1_flavour_cutSxyAndsigmaLxy"), mcdjet.pt(), maxSxy, origin);
      }
      if (prongAcceptance(sigmaLxy, maxSxy, maxSigmaLxy3Prong, minSxy3Prong)) {
        registry.fill(HIST("h3_jet_pt_3prong_Sxyz_N1_flavour_cutSxyzAndsigmaLxyz"), mcdjet.pt(), maxSxyz, origin);
      }
    }
  }

  void processDummy(aod::Collision const&, aod::Tracks const&)
  {
  }
  PROCESS_SWITCH(JetTaggerHFQA, processDummy, "Dummy process", true);

  void processTracksDca(JetTagTracksData& jtracks, OriTracksData const&)
  {
    for (auto const& jtrack : jtracks) {
      if (!jetderiveddatautilities::selectTrack(jtrack, trackSelection)) {
        continue;
      }
      auto track = jtrack.track_as<OriTracksData>();

      float varImpXY, varImpXYSig, varImpZ, varImpZSig, varImpXYZ, varImpXYZSig;
      varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
      varImpXYSig = track.dcaXY() / std::sqrt(track.sigmaDcaXY2());
      varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
      varImpZSig = track.dcaZ() / std::sqrt(track.sigmaDcaZ2());
      float dcaXYZ = jtrack.dcaXYZ();
      float sigmaDcaXYZ2 = jtrack.sigmaDcaXYZ2();
      varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
      varImpXYZSig = dcaXYZ / std::sqrt(sigmaDcaXYZ2);

      registry.fill(HIST("h_impact_parameter_xy"), varImpXY);
      registry.fill(HIST("h_impact_parameter_xy_significance"), varImpXYSig);
      registry.fill(HIST("h_impact_parameter_z"), varImpZ);
      registry.fill(HIST("h_impact_parameter_z_significance"), varImpZSig);
      registry.fill(HIST("h_impact_parameter_xyz"), varImpXYZ);
      registry.fill(HIST("h_impact_parameter_xyz_significance"), varImpXYZSig);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processTracksDca, "Fill inclusive tracks' imformation for data", false);

  void processIPsData(soa::Filtered<JetCollisions>::iterator const& jcollision, JetTagTableData const& jets, JetTagTracksData const& jtracks, OriTracksData const& tracks)
  {
    fillHistogramIPsData(jcollision, jets, jtracks, tracks);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsData, "Fill impact parameter imformation for data jets", false);

  void processIPsMCD(soa::Filtered<JetCollisions>::iterator const& jcollision, JetTagTableMCD const& mcdjets, JetTagTracksMCD const& jtracks, OriTracksMCD const& tracks, JetParticles&)
  {
    fillHistogramIPsMCD(jcollision, mcdjets, jtracks, tracks);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCD, "Fill impact parameter imformation for mcd jets", false);

  void processJPData(soa::Filtered<JetCollisions>::iterator const& jcollision, JetTagTableData const& jets, JetTagTracksData const&)
  {
    fillHistogramJPData(jcollision, jets);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPData, "Fill jet probability imformation for data jets", false);

  void processJPMCD(soa::Filtered<JetCollisions>::iterator const& jcollision, JetTagTableMCD const& mcdjets, JetTagTracksMCD const&)
  {
    fillHistogramJPMCD(jcollision, mcdjets);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCD, "Fill jet probability imformation for mcd jets", false);

  void processSV2ProngData(soa::Filtered<JetCollisions>::iterator const& jcollision, soa::Join<JetTagTableData, aod::DataSecondaryVertex2ProngIndices> const& jets, aod::DataSecondaryVertex2Prongs const& prongs)
  {
    fillHistogramSV2ProngData(jcollision, jets, prongs);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngData, "Fill 2prong imformation for data jets", false);

  void processSV3ProngData(soa::Filtered<JetCollisions>::iterator const& jcollision, soa::Join<JetTagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    fillHistogramSV3ProngData(jcollision, jets, prongs);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngData, "Fill 2prong imformation for data jets", false);

  void processSV2ProngMCD(soa::Filtered<JetCollisions>::iterator const& jcollision, soa::Join<JetTagTableMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    fillHistogramSV2ProngMCD(jcollision, mcdjets, prongs);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCD, "Fill 2prong imformation for mcd jets", false);

  void processSV3ProngMCD(soa::Filtered<JetCollisions>::iterator const& jcollision, soa::Join<JetTagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    fillHistogramSV3ProngMCD(jcollision, mcdjets, prongs);
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCD, "Fill 3prong imformation for mcd jets", false);
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
