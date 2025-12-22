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

/// \file jetTaggerHFQA.cxx
/// \brief Jet tagging general QA
///
/// \author Hanseo Park <hanseo.park@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/DataModel/Multiplicity.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename TagTableData, typename JetTableMCD, typename weightMCD, typename TagTableMCD, typename JetTableMCP, typename weightMCP, typename JetTableMCDMCP, typename JetTableMCPMCD>
struct JetTaggerHFQA {

  // task on/off configuration
  Configurable<bool> fillIPxy{"fillIPxy", true, "process of xy plane of dca"};
  Configurable<bool> fillIPz{"fillIPz", false, "process of z plane of dca"};
  Configurable<bool> fillIPxyz{"fillIPxyz", false, "process of xyz plane of dca"};
  Configurable<bool> fillTrackCounting{"fillTrackCounting", false, "process of track counting method"};
  Configurable<bool> fillGeneralSVQA{"fillGeneralSVQA", true, "process of general QA for sv"};
  Configurable<bool> fillSVxyz{"fillSVxyz", true, "process of decay lenngth of xyz for sv"};

  // Cut configuration
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::vector<float>> trackCuts{"trackCuts", std::vector<float>{0.15, 100.0, -0.9, 0.9}, "Track cuts: ptMin, ptMax, etaMin, etaMax"};
  Configurable<float> trackDcaXYMax{"trackDcaXYMax", 1, "minimum DCA xy acceptance for tracks [cm]"};
  Configurable<float> trackDcaZMax{"trackDcaZMax", 2, "minimum DCA z acceptance for tracks [cm]"};
  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};
  Configurable<std::vector<float>> jetEtaCuts{"jetEtaCuts", std::vector<float>{-99.0, 99.0}, "Jet cuts: etaMin, etaMax"};
  Configurable<std::vector<float>> prongCuts{"prongCuts", std::vector<float>{1, 100, 100, 100, 0.008, 1}, "prong cuts: chi2PCAMin, chi2PCAMax, sigmaLxyMax, sigmaLxyzMax, IPxyMin, IPxyMax"};
  Configurable<float> svDispersionMax{"svDispersionMax", 0.03f, "maximum dispersion of sv"};
  Configurable<int> numFlavourSpecies{"numFlavourSpecies", 6, "number of jet flavour species"};
  Configurable<int> numOrder{"numOrder", 6, "number of ordering"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<float> meanFT0A{"meanFT0A", -1., "Mean value of FT0A signal"};
  Configurable<float> meanFT0C{"meanFT0C", -1., "Mean value of FT0C signal"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // Binning
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, o2::constants::math::TwoPI}, ""};
  ConfigurableAxis binNtracks{"binNtracks", {100, 0., 100.}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binImpactParameterXY{"binImpactParameterXY", {801, -400.5f, 400.5f}, ""};
  ConfigurableAxis binSigmaImpactParameterXY{"binSigmaImpactParameterXY", {800, 0.f, 100.f}, ""};
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
  ConfigurableAxis binMultScaledFT0M{"binMultScaledFT0M", {VARIABLE_WIDTH, 0, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.4, 1.8, 2.4, 3.6, 5., 20.}, "Percentiles of scaled FT0M: 100-90%, 90-80%, 80-70%, 70-60%, 60-50%, 50-40%, 40-30%, 30-20%, 20-10%, 10-1%, 1-0.1%"};

  int numberOfJetFlavourSpecies = 6;
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    numberOfJetFlavourSpecies = static_cast<int>(numFlavourSpecies);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    // Axis
    AxisSpec axisJetFlavour = {binJetFlavour, "Jet flavour"};
    AxisSpec axisJetPt = {binJetPt, "#it{p}_{T, jet}"};
    AxisSpec axisMCDJetPt = {binJetPt, "#it{p}_{T, jet}^{rec}"};
    AxisSpec axisMCPJetPt = {binJetPt, "#it{p}_{T, jet}^{gen}"};
    AxisSpec axisEta = {binEta, "#eta"};
    AxisSpec axisPhi = {binPhi, "#phi"};
    AxisSpec axisNTracks = {binNtracks, "#it{N}_{tracks}"};
    AxisSpec axisTrackPt = {binTrackPt, "#it{p}_{T}^{track}"};
    AxisSpec axisImpactParameterXY = {binImpactParameterXY, "IP_{XY} [#mum]"};
    AxisSpec axisSigmaImpactParameterXY = {binSigmaImpactParameterXY, "#sigma_{XY} [#mum]"};
    AxisSpec axisImpactParameterXYSignificance = {binImpactParameterXYSignificance, "IPs_{XY}"};
    AxisSpec axisImpactParameterZ = {binImpactParameterZ, "IP_{Z} [#mum]"};
    AxisSpec axisImpactParameterZSignificance = {binImpactParameterZSignificance, "IPs_{Z}"};
    AxisSpec axisImpactParameterXYZ = {binImpactParameterXYZ, "IP_{XYZ} [#mum]"};
    AxisSpec axisImpactParameterXYZSignificance = {binImpactParameterXYZSignificance, "IPs_{XYZ}"};
    AxisSpec axisNumOrder = {binNumOrder, "N_{order}"};
    AxisSpec axisJetProbability = {binJetProbability, "JP"};
    AxisSpec axisJetProbabilityLog = {binJetProbabilityLog, "-Log(JP)"};
    AxisSpec axisNprongs = {binNprongs, "#it{N}_{SV}"};
    AxisSpec axisLxy = {binLxy, "L_{XY} [cm]"};
    AxisSpec axisSxy = {binSxy, "S_{XY}"};
    AxisSpec axisLxyz = {binLxyz, "L_{XYZ} [cm]"};
    AxisSpec axisSxyz = {binSxyz, "S_{XYZ}"};
    AxisSpec axisMass = {binMass, "#it{m}_{SV}"};
    AxisSpec axisSigmaLxy = {binSigmaLxy, "#sigma_{L_{XY}} [cm]"};
    AxisSpec axisSigmaLxyz = {binSigmaLxyz, "#sigma_{L_{XYZ}} [cm]"};
    AxisSpec axisFracSecPt = {100, 0, 1, "#frac{#Sigma#it{p}_{T}^{secondary track}}{#it{p}_{T, jet}}"};
    AxisSpec axisMultScaledFT0M = {binMultScaledFT0M, "Multiplicity classes"};

    registry.add("h_collision_events", "data;mcd;mcp evnets", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    if (doprocessTracksDca) {
      if (fillIPxy) {
        registry.add("h_impact_parameter_xy", "", {HistType::kTH1F, {{axisImpactParameterXY}}});
        registry.add("h_impact_parameter_xy_significance", "", {HistType::kTH1F, {{axisImpactParameterXYSignificance}}});
      }
      if (fillIPz) {
        registry.add("h_impact_parameter_z", "", {HistType::kTH1F, {{axisImpactParameterZ}}});
        registry.add("h_impact_parameter_z_significance", "", {HistType::kTH1F, {{axisImpactParameterZSignificance}}});
      }
      if (fillIPxyz) {
        registry.add("h_impact_parameter_xyz", "", {HistType::kTH1F, {{axisImpactParameterXYZ}}});
        registry.add("h_impact_parameter_xyz_significance", "", {HistType::kTH1F, {{axisImpactParameterXYZSignificance}}});
      }
    }
    if (doprocessTracksInJetsData) {
      registry.add("h2_track_pt_impact_parameter_xy", "", {HistType::kTH2F, {{axisTrackPt}, {axisImpactParameterXY}}});
    }
    if (doprocessSecondaryContaminationMCD) {
      registry.add("hn_jet_pt_track_pt_impact_parameter_xy_physical_primary_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisTrackPt}, {axisImpactParameterXY}, {axisJetFlavour}}});
      registry.add("hn_jet_pt_track_pt_impact_parameter_xy_secondary_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisTrackPt}, {axisImpactParameterXY}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_frac_secondary_pt_per_jet_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisFracSecPt}, {axisJetFlavour}}});
    }
    if (doprocessValFlavourDefMCD) {
      registry.add("h3_jet_pt_flavour_dist_quark_flavour_dist_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_flavour_const_quark_flavour_const_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_flavour_const_hadron_flavour_dist_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_flavour_const_quark_flavour_dist_quark", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
    }
    if (doprocessValFlavourDefMCP) {
      registry.add("h3_part_jet_pt_flavour_dist_quark_part_flavour_dist_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_part_jet_pt_flavour_const_quark_part_flavour_const_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_part_jet_pt_flavour_const_hadron_part_flavour_dist_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_part_jet_pt_flavour_const_quark_part_flavour_dist_quark", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
    }
    if (doprocessValFlavourDefMCD) {
      registry.add("h2_flavour_dist_quark_flavour_dist_hadron", "", {HistType::kTH2F, {{axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h2_flavour_const_quark_flavour_const_hadron", "", {HistType::kTH2F, {{axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h2_flavour_const_hadron_flavour_dist_hadron", "", {HistType::kTH2F, {{axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h2_flavour_const_quark_flavour_dist_quark", "", {HistType::kTH2F, {{axisJetFlavour}, {axisJetFlavour}}});
    }
    if (doprocessValFlavourDefMCP) {
      registry.add("h3_part_jet_pt_flavour_dist_quark_part_flavour_dist_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_part_jet_pt_flavour_const_quark_part_flavour_const_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_part_jet_pt_flavour_const_hadron_part_flavour_dist_hadron", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
      registry.add("h3_part_jet_pt_flavour_const_quark_part_flavour_dist_quark", "", {HistType::kTH3F, {{axisJetPt}, {axisJetFlavour}, {axisJetFlavour}}});
    }
    if (doprocessIPsData) {
      registry.add("h_jet_pt", "", {HistType::kTH1F, {{axisJetPt}}});
      registry.add("h_jet_eta", "", {HistType::kTH1F, {{axisEta}}});
      registry.add("h_jet_phi", "", {HistType::kTH1F, {{axisPhi}}});
      registry.add("h3_jet_pt_track_pt_track_eta", "", {HistType::kTH3F, {{axisJetPt}, {axisTrackPt}, {axisEta}}});
      registry.add("h3_jet_pt_track_pt_track_phi", "", {HistType::kTH3F, {{axisJetPt}, {axisTrackPt}, {axisPhi}}});
      if (fillIPxy) {
        registry.add("h2_jet_pt_impact_parameter_xy", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXY}}});
        registry.add("h2_jet_pt_sign_impact_parameter_xy", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXY}}});
        registry.add("h2_jet_pt_impact_parameter_xy_significance", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
        registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xy_significance", "", {HistType::kTH3F, {{axisJetPt}, {axisTrackPt}, {axisImpactParameterXYSignificance}}});
      }
      if (fillIPz) {
        registry.add("h2_jet_pt_impact_parameter_z", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterZ}}});
        registry.add("h2_jet_pt_sign_impact_parameter_z", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterZ}}});
        registry.add("h2_jet_pt_impact_parameter_z_significance", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterZSignificance}}});
        registry.add("h3_jet_pt_track_pt_sign_impact_parameter_z_significance", "", {HistType::kTH3F, {{axisJetPt}, {axisTrackPt}, {axisImpactParameterZSignificance}}});
      }
      if (fillIPxyz) {
        registry.add("h2_jet_pt_impact_parameter_xyz", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYZ}}});
        registry.add("h2_jet_pt_sign_impact_parameter_xyz", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYZ}}});
        registry.add("h2_jet_pt_impact_parameter_xyz_significance", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYZSignificance}}});
        registry.add("h3_jet_pt_track_pt_sign_impact_parameter_xyz_significance", "", {HistType::kTH3F, {{axisJetPt}, {axisTrackPt}, {axisImpactParameterXYZSignificance}}});
      }
      if (fillTrackCounting) {
        if (fillIPxy) {
          registry.add("h2_jet_pt_sign_impact_parameter_xy_significance_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h2_jet_pt_sign_impact_parameter_xy_significance_N2", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h2_jet_pt_sign_impact_parameter_xy_significance_N3", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_tc", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYSignificance}, {axisNumOrder}}});
          registry.add("h3_track_pt_sign_impact_parameter_xy_significance_tc", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYSignificance}, {axisNumOrder}}});
        }
        if (fillIPz) {
          registry.add("h2_jet_pt_sign_impact_parameter_z_significance_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h2_jet_pt_sign_impact_parameter_z_significance_N2", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h2_jet_pt_sign_impact_parameter_z_significance_N3", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h3_jet_pt_sign_impact_parameter_z_significance_tc", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZSignificance}, {axisNumOrder}}});
          registry.add("h3_track_pt_sign_impact_parameter_z_significance_tc", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterZSignificance}, {axisNumOrder}}});
        }
        if (fillIPxyz) {
          registry.add("h2_jet_pt_sign_impact_parameter_xyz_significance_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h2_jet_pt_sign_impact_parameter_xyz_significance_N2", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h2_jet_pt_sign_impact_parameter_xyz_significance_N3", "", {HistType::kTH2F, {{axisJetPt}, {axisImpactParameterXYSignificance}}});
          registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_tc", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZSignificance}, {axisNumOrder}}});
          registry.add("h3_track_pt_sign_impact_parameter_xyz_significance_tc", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYZSignificance}, {axisNumOrder}}});
        }
      }
    }
    if (doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) {
      registry.add("h2_jet_pt_flavour", "", {HistType::kTH2F, {{axisJetPt}, {axisJetFlavour}}});
      registry.add("h2_jet_eta_flavour", "", {HistType::kTH2F, {{axisEta}, {axisJetFlavour}}});
      registry.add("h2_jet_phi_flavour", "", {HistType::kTH2F, {{axisPhi}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_track_pt_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisTrackPt}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_track_eta_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisEta}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_track_phi_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisPhi}, {axisJetFlavour}}});
      if (fillIPxy) {
        registry.add("h3_jet_pt_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXY}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sigma_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSigmaImpactParameterXY}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXY}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
        registry.add("h3_track_pt_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXY}, {axisJetFlavour}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXY}, {axisJetFlavour}}});
        registry.add("h3_track_pt_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
        registry.add("h3_track_pt_sign_impact_parameter_xy_significance_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
      }
      if (fillIPz) {
        registry.add("h3_jet_pt_impact_parameter_z_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZ}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZ}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
        registry.add("h3_track_pt_impact_parameter_z_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterZ}, {axisJetFlavour}}});
        registry.add("h3_track_pt_sign_impact_parameter_z_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterZ}, {axisJetFlavour}}});
        registry.add("h3_track_pt_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
        registry.add("h3_track_pt_sign_impact_parameter_z_significance_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
      }
      if (fillIPxyz) {
        registry.add("h3_jet_pt_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZ}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZ}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
        registry.add("h3_track_pt_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYZ}, {axisJetFlavour}}});
        registry.add("h3_track_pt_sign_impact_parameter_xyz_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYZ}, {axisJetFlavour}}});
        registry.add("h3_track_pt_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
        registry.add("h3_track_pt_sign_impact_parameter_xyz_significance_flavour", "", {HistType::kTH3F, {{axisTrackPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
      }
      if (fillTrackCounting) {
        if (fillIPxy) {
          registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
          registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
          registry.add("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N3", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYSignificance}, {axisJetFlavour}}});
          registry.add("h3_sign_impact_parameter_xy_significance_tc_flavour", "", {HistType::kTH3F, {{axisImpactParameterXYSignificance}, {axisNumOrder}, {axisJetFlavour}}});
        }
        if (fillIPz) {
          registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
          registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
          registry.add("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N3", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterZSignificance}, {axisJetFlavour}}});
          registry.add("h3_sign_impact_parameter_z_significance_tc_flavour", "", {HistType::kTH3F, {{axisImpactParameterZSignificance}, {axisNumOrder}, {axisJetFlavour}}});
        }
        if (fillIPxyz) {
          registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
          registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
          registry.add("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N3", "", {HistType::kTH3F, {{axisJetPt}, {axisImpactParameterXYZSignificance}, {axisJetFlavour}}});
          registry.add("h3_sign_impact_parameter_xyz_significance_tc_flavour", "", {HistType::kTH3F, {{axisImpactParameterXYZSignificance}, {axisNumOrder}, {axisJetFlavour}}});
        }
      }
    }
    if (doprocessIPsMCP || doprocessIPsMCPWeighted) {
      registry.add("h2_jet_pt_part_flavour", "", {HistType::kTH2F, {{axisJetPt}, {axisJetFlavour}}});
      registry.add("h2_jet_eta_part_flavour", "", {HistType::kTH2F, {{axisEta}, {axisJetFlavour}}});
      registry.add("h2_jet_phi_part_flavour", "", {HistType::kTH2F, {{axisPhi}, {axisJetFlavour}}});
    }
    if (doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) {
      registry.add("h3_jet_pt_jet_pt_part_matchedgeo_flavour", "", {HistType::kTH3F, {{axisMCDJetPt}, {axisMCPJetPt}, {axisJetFlavour}}});
    }
    if (doprocessJPData) {
      if (!doprocessIPsData && !doprocessSV2ProngData && !doprocessSV3ProngData) {
        registry.add("h_jet_pt", "", {HistType::kTH1F, {{axisJetPt}}});
        registry.add("h_jet_eta", "", {HistType::kTH1F, {{axisEta}}});
        registry.add("h_jet_phi", "", {HistType::kTH1F, {{axisPhi}}});
      }
      registry.add("h2_jet_pt_JP", "jet pt jet probability untagged", {HistType::kTH2F, {{axisJetPt}, {axisJetProbability}}});
      registry.add("h2_jet_pt_neg_log_JP", "jet pt jet probabilityun tagged", {HistType::kTH2F, {{axisJetPt}, {axisJetProbabilityLog}}});
      registry.add("h2_taggedjet_pt_JP_N1", "jet pt jet probability N1", {HistType::kTH2F, {{axisJetPt}, {axisJetProbability}}});
      registry.add("h2_taggedjet_pt_neg_log_JP_N1", "jet pt jet probabilityun N1", {HistType::kTH2F, {{axisJetPt}, {axisJetProbabilityLog}}});
      registry.add("h2_taggedjet_pt_JP_N2", "jet pt jet probability N2", {HistType::kTH2F, {{axisJetPt}, {axisJetProbability}}});
      registry.add("h2_taggedjet_pt_neg_log_JP_N2", "jet pt jet probabilityun N2", {HistType::kTH2F, {{axisJetPt}, {axisJetProbabilityLog}}});
      registry.add("h2_taggedjet_pt_JP_N3", "jet pt jet probability N3", {HistType::kTH2F, {{axisJetPt}, {axisJetProbability}}});
      registry.add("h2_taggedjet_pt_neg_log_JP_N3", "jet pt jet probabilityun N3", {HistType::kTH2F, {{axisJetPt}, {axisJetProbabilityLog}}});
    }
    if (doprocessJPMCD || doprocessJPMCDWeighted || doprocessJPMCPMCDMatched || doprocessJPMCPMCDMatchedWeighted) {
      if (!(doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) && !(doprocessSV2ProngMCD || doprocessSV2ProngMCDWeighted || doprocessSV2ProngMCPMCDMatched || doprocessSV2ProngMCPMCDMatchedWeighted) && !(doprocessSV3ProngMCD || doprocessSV3ProngMCDWeighted || doprocessSV3ProngMCPMCDMatched || doprocessSV3ProngMCPMCDMatchedWeighted)) {
        registry.add("h2_jet_pt_flavour", "", {HistType::kTH2F, {{axisJetPt}, {axisJetFlavour}}});
        registry.add("h2_jet_eta_flavour", "", {HistType::kTH2F, {{axisEta}, {axisJetFlavour}}});
        registry.add("h2_jet_phi_flavour", "", {HistType::kTH2F, {{axisPhi}, {axisJetFlavour}}});
      }
      registry.add("h3_jet_pt_JP_flavour", "jet pt jet probability flavour untagged", {HistType::kTH3F, {{axisJetPt}, {axisJetProbability}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_neg_log_JP_flavour", "jet pt log jet probability flavour untagged", {HistType::kTH3F, {{axisJetPt}, {axisJetProbabilityLog}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_JP_N1_flavour", "jet pt jet probability flavour N1", {HistType::kTH3F, {{axisJetPt}, {axisJetProbability}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_neg_log_JP_N1_flavour", "jet pt log jet probability flavour N1", {HistType::kTH3F, {{axisJetPt}, {axisJetProbabilityLog}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_JP_N2_flavour", "jet pt jet probability flavour N2", {HistType::kTH3F, {{axisJetPt}, {axisJetProbability}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_neg_log_JP_N2_flavour", "jet pt log jet probability flavour N2", {HistType::kTH3F, {{axisJetPt}, {axisJetProbabilityLog}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_JP_N3_flavour", "jet pt jet probability flavour N3", {HistType::kTH3F, {{axisJetPt}, {axisJetProbability}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_neg_log_JP_N3_flavour", "jet pt log jet probability flavour N3", {HistType::kTH3F, {{axisJetPt}, {axisJetProbabilityLog}, {axisJetFlavour}}});
    }
    if (doprocessSV2ProngData) {
      if (!doprocessIPsData && !doprocessJPData && !doprocessSV3ProngData) {
        registry.add("h_jet_pt", "", {HistType::kTH1F, {{axisJetPt}}});
        registry.add("h_jet_eta", "", {HistType::kTH1F, {{axisEta}}});
        registry.add("h_jet_phi", "", {HistType::kTH1F, {{axisPhi}}});
      }
      if (fillGeneralSVQA) {
        registry.add("h_2prong_nprongs", "", {HistType::kTH1F, {{axisNprongs}}});
        registry.add("h2_jet_pt_2prong_Lxy", "", {HistType::kTH2F, {{axisJetPt}, {axisLxy}}});
        registry.add("h2_jet_pt_2prong_sigmaLxy", "", {HistType::kTH2F, {{axisJetPt}, {axisSigmaLxy}}});
        registry.add("h2_jet_pt_2prong_Sxy", "", {HistType::kTH2F, {{axisJetPt}, {axisSxy}}});
        registry.add("h2_jet_pt_2prong_Lxyz", "", {HistType::kTH2F, {{axisJetPt}, {axisLxyz}}});
        registry.add("h2_jet_pt_2prong_sigmaLxyz", "", {HistType::kTH2F, {{axisJetPt}, {axisSigmaLxyz}}});
        registry.add("h2_jet_pt_2prong_Sxyz", "", {HistType::kTH2F, {{axisJetPt}, {axisSxyz}}});
      }
      registry.add("h2_jet_pt_2prong_Sxy_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxy}}});
      registry.add("h2_jet_pt_2prong_Sxyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxyz}}});
      registry.add("h2_jet_pt_2prong_mass_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
      registry.add("h2_jet_pt_2prong_mass_xyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
      registry.add("h2_taggedjet_pt_2prong_Sxy_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxy}}});
      registry.add("h2_taggedjet_pt_2prong_Sxyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxyz}}});
      registry.add("h2_taggedjet_pt_2prong_mass_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
      registry.add("h2_taggedjet_pt_2prong_mass_xyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
    }
    if (doprocessSV3ProngData) {
      if (!doprocessIPsData && !doprocessJPData && !doprocessSV2ProngData) {
        registry.add("h_jet_pt", "", {HistType::kTH1F, {{axisJetPt}}});
        registry.add("h_jet_eta", "", {HistType::kTH1F, {{axisEta}}});
        registry.add("h_jet_phi", "", {HistType::kTH1F, {{axisPhi}}});
      }
      if (fillGeneralSVQA) {
        registry.add("h_3prong_nprongs", "", {HistType::kTH1F, {{axisNprongs}}});
        registry.add("h2_jet_pt_3prong_Lxy", "", {HistType::kTH2F, {{axisJetPt}, {axisLxy}}});
        registry.add("h2_jet_pt_3prong_sigmaLxy", "", {HistType::kTH2F, {{axisJetPt}, {axisSigmaLxy}}});
        registry.add("h2_jet_pt_3prong_Sxy", "", {HistType::kTH2F, {{axisJetPt}, {axisSxy}}});
        registry.add("h2_jet_pt_3prong_Lxyz", "", {HistType::kTH2F, {{axisJetPt}, {axisLxyz}}});
        registry.add("h2_jet_pt_3prong_sigmaLxyz", "", {HistType::kTH2F, {{axisJetPt}, {axisSigmaLxyz}}});
        registry.add("h2_jet_pt_3prong_Sxyz", "", {HistType::kTH2F, {{axisJetPt}, {axisSxyz}}});
      }
      registry.add("h2_jet_pt_3prong_Sxy_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxy}}});
      registry.add("h2_jet_pt_3prong_Sxyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxyz}}});
      registry.add("h2_jet_pt_3prong_mass_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
      registry.add("h2_jet_pt_3prong_mass_xyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
      registry.add("h2_taggedjet_pt_3prong_Sxy_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxy}}});
      registry.add("h2_taggedjet_pt_3prong_Sxyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisSxyz}}});
      registry.add("h2_taggedjet_pt_3prong_mass_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
      registry.add("h2_taggedjet_pt_3prong_mass_xyz_N1", "", {HistType::kTH2F, {{axisJetPt}, {axisMass}}});
    }
    if (doprocessSV3ProngDataMult) {
      registry.add("h_event_mult", "", {HistType::kTH1F, {{axisMultScaledFT0M}}});
      registry.add("h2_jet_pt_mult", "", {HistType::kTH2F, {{axisJetPt}, {axisMultScaledFT0M}}});
      registry.add("h2_jet_eta_mult", "", {HistType::kTH2F, {{axisEta}, {axisMultScaledFT0M}}});
      registry.add("h2_jet_phi_mult", "", {HistType::kTH2F, {{axisPhi}, {axisMultScaledFT0M}}});
      if (fillGeneralSVQA) {
        registry.add("h2_3prong_nprongs_mult", "", {HistType::kTH2F, {{axisNprongs}, {axisMultScaledFT0M}}});
        registry.add("hn_jet_3prong_Sxy_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisMultScaledFT0M}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisMultScaledFT0M}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisMultScaledFT0M}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisMultScaledFT0M}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisMultScaledFT0M}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisMultScaledFT0M}}});
      }
    }
    if (doprocessSV2ProngMCD || doprocessSV2ProngMCDWeighted || doprocessSV2ProngMCPMCDMatched || doprocessSV2ProngMCPMCDMatchedWeighted) {
      if (!(doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) && !(doprocessJPMCD || doprocessJPMCDWeighted || doprocessJPMCPMCDMatched || doprocessJPMCPMCDMatchedWeighted) && !(doprocessSV3ProngMCD || doprocessSV3ProngMCDWeighted || doprocessSV3ProngMCPMCDMatched || doprocessSV3ProngMCPMCDMatchedWeighted)) {
        registry.add("h2_jet_pt_flavour", "", {HistType::kTH2F, {{axisJetPt}, {axisJetFlavour}}});
        registry.add("h2_jet_eta_flavour", "", {HistType::kTH2F, {{axisEta}, {axisJetFlavour}}});
        registry.add("h2_jet_phi_flavour", "", {HistType::kTH2F, {{axisPhi}, {axisJetFlavour}}});
      }
      if (fillGeneralSVQA) {
        registry.add("h2_2prong_nprongs_flavour", "", {HistType::kTH2F, {{axisNprongs}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_2prong_Lxy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisLxy}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_2prong_sigmaLxy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSigmaLxy}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_2prong_Sxy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxy}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_2prong_Lxyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisLxyz}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_2prong_sigmaLxyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSigmaLxyz}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_2prong_Sxyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxyz}, {axisJetFlavour}}});
      }
      registry.add("h3_jet_pt_2prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxy}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_2prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxyz}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_2prong_mass_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_2prong_mass_xyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_2prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxy}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_2prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxyz}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_2prong_mass_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_2prong_mass_xyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
    }
    if (doprocessSV3ProngMCD || doprocessSV3ProngMCDWeighted || doprocessSV3ProngMCPMCDMatched || doprocessSV3ProngMCPMCDMatchedWeighted) {
      if (!(doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) && !(doprocessJPMCD || doprocessJPMCDWeighted || doprocessJPMCPMCDMatched || doprocessJPMCPMCDMatchedWeighted) && !(doprocessSV2ProngMCD || doprocessSV2ProngMCDWeighted || doprocessSV2ProngMCPMCDMatched || doprocessSV2ProngMCPMCDMatchedWeighted)) {
        registry.add("h2_jet_pt_flavour", "", {HistType::kTH2F, {{axisJetPt}, {axisJetFlavour}}});
        registry.add("h2_jet_eta_flavour", "", {HistType::kTH2F, {{axisEta}, {axisJetFlavour}}});
        registry.add("h2_jet_phi_flavour", "", {HistType::kTH2F, {{axisPhi}, {axisJetFlavour}}});
      }
      if (fillGeneralSVQA) {
        registry.add("h2_3prong_nprongs_flavour", "", {HistType::kTH2F, {{axisNprongs}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_3prong_Lxy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisLxy}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_3prong_sigmaLxy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSigmaLxy}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_3prong_Sxy_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxy}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_3prong_Lxyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisLxyz}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_3prong_sigmaLxyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSigmaLxyz}, {axisJetFlavour}}});
        registry.add("h3_jet_pt_3prong_Sxyz_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxyz}, {axisJetFlavour}}});
      }
      registry.add("h3_jet_pt_3prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxy}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_3prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxyz}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_3prong_mass_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
      registry.add("h3_jet_pt_3prong_mass_xyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_3prong_Sxy_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxy}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_3prong_Sxyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisSxyz}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_3prong_mass_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
      registry.add("h3_taggedjet_pt_3prong_mass_xyz_N1_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMass}, {axisJetFlavour}}});
    }
    if (doprocessSV3ProngMCDMult || doprocessSV3ProngMCDMultWeighted || doprocessSV3ProngMCPMCDMatchedMult || doprocessSV3ProngMCPMCDMatchedMultWeighted) {
      registry.add("h_event_mult", "", {HistType::kTH1F, {{axisMultScaledFT0M}}});
      registry.add("h3_jet_pt_mult_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMultScaledFT0M}, {axisJetFlavour}}});
      registry.add("h3_jet_eta_mult_flavour", "", {HistType::kTH3F, {{axisEta}, {axisMultScaledFT0M}, {axisJetFlavour}}});
      registry.add("h3_jet_phi_mult_flavour", "", {HistType::kTH3F, {{axisPhi}, {axisMultScaledFT0M}, {axisJetFlavour}}});
      if (fillGeneralSVQA) {
        registry.add("h3_3prong_nprongs_mult_flavour", "", {HistType::kTH3F, {{axisNprongs}, {axisMultScaledFT0M}, {axisJetFlavour}}});
        registry.add("hn_jet_3prong_Sxy_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisMultScaledFT0M}, {axisJetFlavour}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisMultScaledFT0M}, {axisJetFlavour}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisMultScaledFT0M}, {axisJetFlavour}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisMultScaledFT0M}, {axisJetFlavour}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisMultScaledFT0M}, {axisJetFlavour}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisMultScaledFT0M}, {axisJetFlavour}}});
      }
    }
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackCuts->at(0) && aod::jtrack::pt < trackCuts->at(1) && aod::jtrack::eta > trackCuts->at(2) && aod::jtrack::eta < trackCuts->at(3));
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<aod::JetParticles> particlesPerCollision = aod::jmcparticle::mcCollisionId;

  using JetTagTracksData = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using JetTagTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;

  std::function<bool(const std::vector<float>&, const std::vector<float>&)> sortImp =
    [](const std::vector<float>& a, const std::vector<float>& b) {
      return a[0] > b[0];
    };

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {
    const float noJetAreaFractionFilter = -98.0;
    const float noConstituentPtMinFilter = -98.0;
    const float noConstituentPtMaxFilter = 9998.0;
    if (jetAreaFractionMin > noJetAreaFractionFilter) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > noConstituentPtMinFilter);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < noConstituentPtMaxFilter);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<T>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }

    return true;
  }

  template <typename T>
  bool trackAcceptance(T const& track)
  {
    if (track.pt() < trackCuts->at(0) || track.pt() > trackCuts->at(1))
      return false;

    return true;
  }

  float getScaledFT0A(const float multFT0A)
  {
    return multFT0A / meanFT0A;
  }

  float getScaledFT0C(const float multFT0C)
  {
    return multFT0C / meanFT0C;
  }

  float getScaledFT0M(const float multFT0A, const float multFT0C)
  {
    return 0.5 * (getScaledFT0A(multFT0A) + getScaledFT0C(multFT0C));
  }

  template <typename U, typename T, typename V, typename W, typename X>
  void fillValidationFlavourDefMCD(T const& mcdjet, V const& tracks, W const& particles, X const& particlesPerColl, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    typename V::iterator hftrack;
    int jetflavourConstQuark = jettaggingutilities::mcdJetFromHFShower(mcdjet, tracks, particles, maxDeltaR, true);
    int jetflavourConstHadron = jettaggingutilities::mcdJetFromHFShower(mcdjet, tracks, particles, maxDeltaR, false);
    int jetflavourDistQuark = -1;
    int jetflavourDistHadron = -1;
    for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<U>()) {
      jetflavourDistQuark = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
      jetflavourDistHadron = jettaggingutilities::getJetFlavorHadron(mcpjet, particlesPerColl);
    }
    if (jetflavourDistQuark < 0 || jetflavourDistHadron < 0)
      return;
    registry.fill(HIST("h3_jet_pt_flavour_dist_quark_flavour_dist_hadron"), mcdjet.pt(), jetflavourDistQuark, jetflavourDistHadron, eventWeight);
    registry.fill(HIST("h3_jet_pt_flavour_const_quark_flavour_const_hadron"), mcdjet.pt(), jetflavourConstQuark, jetflavourConstHadron, eventWeight);
    registry.fill(HIST("h3_jet_pt_flavour_const_hadron_flavour_dist_hadron"), mcdjet.pt(), jetflavourConstHadron, jetflavourDistHadron, eventWeight);
    registry.fill(HIST("h3_jet_pt_flavour_const_quark_flavour_dist_quark"), mcdjet.pt(), jetflavourConstQuark, jetflavourDistQuark, eventWeight);
  }

  template <typename T, typename U, typename V>
  void fillValidationFlavourDefMCP(T const& mcpjet, U const& particles, V const& particlesPerColl, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
      return;
    }
    int jetflavourConstQuark = jettaggingutilities::mcpJetFromHFShower(mcpjet, particles, maxDeltaR, true);
    int jetflavourConstHadron = jettaggingutilities::mcpJetFromHFShower(mcpjet, particles, maxDeltaR, false);
    int jetflavourDistQuark = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
    int jetflavourDistHadron = jettaggingutilities::getJetFlavorHadron(mcpjet, particlesPerColl);
    registry.fill(HIST("h3_part_jet_pt_flavour_dist_quark_part_flavour_dist_hadron"), mcpjet.pt(), jetflavourDistQuark, jetflavourDistHadron, eventWeight);
    registry.fill(HIST("h3_part_jet_pt_flavour_const_quark_part_flavour_const_hadron"), mcpjet.pt(), jetflavourConstQuark, jetflavourConstHadron, eventWeight);
    registry.fill(HIST("h3_part_jet_pt_flavour_const_hadron_part_flavour_dist_hadron"), mcpjet.pt(), jetflavourConstHadron, jetflavourDistHadron, eventWeight);
    registry.fill(HIST("h3_part_jet_pt_flavour_const_quark_part_flavour_dist_quark"), mcpjet.pt(), jetflavourConstQuark, jetflavourDistQuark, eventWeight);
  }

  template <typename T, typename U>
  void fillHistogramIPsData(T const& jet, U const& /*tracks*/)
  {
    std::size_t firstTaggerForTrackCounting = 0;
    std::size_t secondTaggerForTrackCounting = 1;
    std::size_t thirdTaggerForTrackCounting = 2;
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
    std::vector<std::vector<float>> vecSignImpXYSig, vecSignImpZSig, vecSignImpXYZSig;
    for (auto const& track : jet.template tracks_as<U>()) {
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

    if (vecSignImpXYSig.size() > firstTaggerForTrackCounting) { // N1
      if (fillIPxy)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xy_significance_N1"), jet.pt(), vecSignImpXYSig[0][0]);
      if (fillIPz)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_z_significance_N1"), jet.pt(), vecSignImpZSig[0][0]);
      if (fillIPxyz)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xyz_significance_N1"), jet.pt(), vecSignImpXYZSig[0][0]);
    }
    if (vecSignImpXYSig.size() > secondTaggerForTrackCounting) { // N2
      if (fillIPxy)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xy_significance_N2"), jet.pt(), vecSignImpXYSig[1][0]);
      if (fillIPz)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_z_significance_N2"), jet.pt(), vecSignImpZSig[1][0]);
      if (fillIPxyz)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xyz_significance_N2"), jet.pt(), vecSignImpXYZSig[1][0]);
    }
    if (vecSignImpXYSig.size() > thirdTaggerForTrackCounting) { // N3
      if (fillIPxy)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xy_significance_N3"), jet.pt(), vecSignImpXYSig[2][0]);
      if (fillIPz)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_z_significance_N3"), jet.pt(), vecSignImpZSig[2][0]);
      if (fillIPxyz)
        registry.fill(HIST("h2_jet_pt_sign_impact_parameter_xyz_significance_N3"), jet.pt(), vecSignImpXYZSig[2][0]);
    }
    if (fillIPxy && vecSignImpXYSig.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxy && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSig.size()) {
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_tc"), jet.pt(), vecSignImpXYSig[order - 1][0], order);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xy_significance_tc"), vecSignImpXYSig[order - 1][1], vecSignImpXYSig[order - 1][0], order);
      }
    }
    if (fillIPz && vecSignImpZSig.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPz && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSig.size()) {
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_tc"), jet.pt(), vecSignImpZSig[order - 1][0], order);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_z_significance_tc"), vecSignImpZSig[order - 1][1], vecSignImpZSig[order - 1][0], order);
      }
    }
    if (fillIPxyz && vecSignImpXYZSig.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxyz && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSig.size()) {
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_tc"), jet.pt(), vecSignImpXYZSig[order - 1][0], order);
        registry.fill(HIST("h3_track_pt_sign_impact_parameter_xyz_significance_tc"), vecSignImpXYZSig[order - 1][1], vecSignImpXYZSig[order - 1][0], order);
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramIPsMCD(T const& mcdjet, U const& /*tracks*/, float eventWeight = 1.0)
  {
    std::size_t firstTaggerForTrackCounting = 0;
    std::size_t secondTaggerForTrackCounting = 1;
    std::size_t thirdTaggerForTrackCounting = 2;
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
    if (jetflavour == -1) {
      jetflavour = JetTaggingSpecies::none;
    }
    registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), jetflavour, eventWeight);
    registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), jetflavour, eventWeight);
    registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), jetflavour, eventWeight);
    for (auto const& track : mcdjet.template tracks_as<U>()) {
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
        float sigmadcaXYZ = track.sigmadcaXYZ();
        varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
        varSignImpXYZ = geoSign * std::abs(dcaXYZ) * jettaggingutilities::cmTomum;
        varImpXYZSig = dcaXYZ / sigmadcaXYZ;
        varSignImpXYZSig = geoSign * std::abs(dcaXYZ) / sigmadcaXYZ;
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
    sort(vecSignImpXYSig[jetflavour].begin(), vecSignImpXYSig[jetflavour].end(), std::greater<float>());
    sort(vecSignImpZSig[jetflavour].begin(), vecSignImpZSig[jetflavour].end(), std::greater<float>());
    sort(vecSignImpXYZSig[jetflavour].begin(), vecSignImpXYZSig[jetflavour].end(), std::greater<float>());

    if (vecImpXY[jetflavour].size() > firstTaggerForTrackCounting) { // N1
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N1"), mcdjet.pt(), vecSignImpXYSig[jetflavour][0], jetflavour, eventWeight);
      if (fillIPz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N1"), mcdjet.pt(), vecSignImpZSig[jetflavour][0], jetflavour, eventWeight);
      if (fillIPxyz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N1"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][0], jetflavour, eventWeight);
    }
    if (vecImpXY[jetflavour].size() > secondTaggerForTrackCounting) { // N2
      if (fillIPxy)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xy_significance_flavour_N2"), mcdjet.pt(), vecSignImpXYSig[jetflavour][1], jetflavour, eventWeight);
      if (fillIPz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_z_significance_flavour_N2"), mcdjet.pt(), vecSignImpZSig[jetflavour][1], jetflavour, eventWeight);
      if (fillIPxyz)
        registry.fill(HIST("h3_jet_pt_sign_impact_parameter_xyz_significance_flavour_N2"), mcdjet.pt(), vecSignImpXYZSig[jetflavour][1], jetflavour, eventWeight);
    }
    if (vecImpXY[jetflavour].size() > thirdTaggerForTrackCounting) { // N3
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
      if (fillIPz && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSigTC.size()) {
        registry.fill(HIST("h3_sign_impact_parameter_z_significance_tc_flavour"), vecSignImpZSigTC[order - 1][0], order, jetflavour, eventWeight);
      }
    }

    if (vecSignImpXYZSigTC.empty())
      return;
    for (int order = 1; order <= numOrder; order++) {
      if (fillIPxyz && static_cast<std::vector<std::vector<float>>::size_type>(order) < vecSignImpXYSigTC.size()) {
        registry.fill(HIST("h3_sign_impact_parameter_xyz_significance_tc_flavour"), vecSignImpXYZSigTC[order - 1][0], order, jetflavour, eventWeight);
      }
    }
  }

  template <typename T>
  void fillHistogramIPsMCP(T const& mcpjet, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
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

  template <typename T>
  void fillHistogramJPData(T const& jet)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (!doprocessIPsData && !doprocessSV2ProngData && !doprocessSV3ProngData) {
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());
    }
    registry.fill(HIST("h2_jet_pt_JP"), jet.pt(), jet.jetProb());
    registry.fill(HIST("h2_jet_pt_neg_log_JP"), jet.pt(), -1 * std::log(jet.jetProb()));
    registry.fill(HIST("h2_taggedjet_pt_JP_N1"), jet.pt(), jet.isTagged(BJetTaggingMethod::IPsN1) ? jet.jetProb() : -1);
    registry.fill(HIST("h2_taggedjet_pt_neg_log_JP_N1"), jet.pt(), jet.isTagged(BJetTaggingMethod::IPsN1) ? -1 * std::log(jet.jetProb()) : -1);
    registry.fill(HIST("h2_taggedjet_pt_JP_N2"), jet.pt(), jet.isTagged(BJetTaggingMethod::IPsN2) ? jet.jetProb() : -1);
    registry.fill(HIST("h2_taggedjet_pt_neg_log_JP_N2"), jet.pt(), jet.isTagged(BJetTaggingMethod::IPsN2) ? -1 * std::log(jet.jetProb()) : -1);
    registry.fill(HIST("h2_taggedjet_pt_JP_N3"), jet.pt(), jet.isTagged(BJetTaggingMethod::IPsN3) ? jet.jetProb() : -1);
    registry.fill(HIST("h2_taggedjet_pt_neg_log_JP_N3"), jet.pt(), jet.isTagged(BJetTaggingMethod::IPsN3) ? -1 * std::log(jet.jetProb()) : -1);
  }

  template <typename T>
  void fillHistogramJPMCD(T const& mcdjet, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (!((doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) && (doprocessSV2ProngMCD || doprocessSV2ProngMCDWeighted || doprocessSV2ProngMCPMCDMatched || doprocessSV2ProngMCPMCDMatchedWeighted) && (doprocessSV3ProngMCD || doprocessSV3ProngMCDWeighted || doprocessSV3ProngMCPMCDMatched || doprocessSV3ProngMCPMCDMatchedWeighted))) {
      registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), mcdjet.origin(), eventWeight);
      registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), mcdjet.origin(), eventWeight);
      registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), mcdjet.origin(), eventWeight);
    }
    registry.fill(HIST("h3_jet_pt_JP_flavour"), mcdjet.pt(), mcdjet.jetProb(), mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_jet_pt_neg_log_JP_flavour"), mcdjet.pt(), -1 * std::log(mcdjet.jetProb()), mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_JP_N1_flavour"), mcdjet.pt(), mcdjet.isTagged(BJetTaggingMethod::IPsN1) ? mcdjet.jetProb() : -1, mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_neg_log_JP_N1_flavour"), mcdjet.pt(), mcdjet.isTagged(BJetTaggingMethod::IPsN1) ? -1 * std::log(mcdjet.jetProb()) : -1, mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_JP_N2_flavour"), mcdjet.pt(), mcdjet.isTagged(BJetTaggingMethod::IPsN2) ? mcdjet.jetProb() : -1, mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_neg_log_JP_N2_flavour"), mcdjet.pt(), mcdjet.isTagged(BJetTaggingMethod::IPsN2) ? -1 * std::log(mcdjet.jetProb()) : -1, mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_JP_N3_flavour"), mcdjet.pt(), mcdjet.isTagged(BJetTaggingMethod::IPsN3) ? mcdjet.jetProb() : -1, mcdjet.origin(), eventWeight);
    registry.fill(HIST("h3_taggedjet_pt_neg_log_JP_N3_flavour"), mcdjet.pt(), mcdjet.isTagged(BJetTaggingMethod::IPsN3) ? -1 * std::log(mcdjet.jetProb()) : -1, mcdjet.origin(), eventWeight);
  }

  template <typename T, typename U>
  void fillHistogramSV2ProngData(T const& jet, U const& /*prongs*/)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (jet.template secondaryVertices_as<U>().size() < 1)
      return;
    if (!doprocessIPsData && !doprocessJPData && !doprocessSV3ProngData) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), eventWeight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), eventWeight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), eventWeight);
    }
    if (fillGeneralSVQA) {
      registry.fill(HIST("h_2prong_nprongs"), jet.template secondaryVertices_as<U>().size());
      for (const auto& prong : jet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("h2_jet_pt_2prong_Lxy"), jet.pt(), prong.decayLengthXY());
        registry.fill(HIST("h2_jet_pt_2prong_sigmaLxy"), jet.pt(), prong.errorDecayLengthXY());
        registry.fill(HIST("h2_jet_pt_2prong_Sxy"), jet.pt(), prong.decayLengthXY() / prong.errorDecayLengthXY());
        if (fillSVxyz) {
          registry.fill(HIST("h2_jet_pt_2prong_Lxyz"), jet.pt(), prong.decayLength());
          registry.fill(HIST("h2_jet_pt_2prong_sigmaLxyz"), jet.pt(), prong.errorDecayLength());
          registry.fill(HIST("h2_jet_pt_2prong_Sxyz"), jet.pt(), prong.decayLength() / prong.errorDecayLength());
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("h2_jet_pt_2prong_Sxy_N1"), jet.pt(), maxSxy);
      registry.fill(HIST("h2_jet_pt_2prong_mass_N1"), jet.pt(), massSV);
      if (jet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("h2_taggedjet_pt_2prong_Sxy_N1"), jet.pt(), maxSxy);
        registry.fill(HIST("h2_taggedjet_pt_2prong_mass_N1"), jet.pt(), massSV);
      }
    }
    if (fillSVxyz) {
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("h2_jet_pt_2prong_Sxyz_N1"), jet.pt(), maxSxyz);
        registry.fill(HIST("h2_jet_pt_2prong_mass_xyz_N1"), jet.pt(), massSV);
        if (jet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("h2_taggedjet_pt_2prong_Sxyz_N1"), jet.pt(), maxSxyz);
          registry.fill(HIST("h2_taggedjet_pt_2prong_mass_xyz_N1"), jet.pt(), massSV);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramSV3ProngData(T const& jet, U const& /*prongs*/)
  {
    float eventWeight = 1.0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (jet.template secondaryVertices_as<U>().size() < 1)
      return;
    if (!doprocessIPsData && !doprocessJPData && !doprocessSV2ProngData) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), eventWeight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), eventWeight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), eventWeight);
    }
    if (fillGeneralSVQA) {
      registry.fill(HIST("h_3prong_nprongs"), jet.template secondaryVertices_as<U>().size());
      for (const auto& prong : jet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("h2_jet_pt_3prong_Lxy"), jet.pt(), prong.decayLengthXY());
        registry.fill(HIST("h2_jet_pt_3prong_sigmaLxy"), jet.pt(), prong.errorDecayLengthXY());
        registry.fill(HIST("h2_jet_pt_3prong_Sxy"), jet.pt(), prong.decayLengthXY() / prong.errorDecayLengthXY());
        if (fillSVxyz) {
          registry.fill(HIST("h2_jet_pt_3prong_Lxyz"), jet.pt(), prong.decayLength());
          registry.fill(HIST("h2_jet_pt_3prong_sigmaLxyz"), jet.pt(), prong.errorDecayLength());
          registry.fill(HIST("h2_jet_pt_3prong_Sxyz"), jet.pt(), prong.decayLength() / prong.errorDecayLength());
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("h2_jet_pt_3prong_Sxy_N1"), jet.pt(), maxSxy);
      registry.fill(HIST("h2_jet_pt_3prong_mass_N1"), jet.pt(), massSV);
      if (jet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("h2_taggedjet_pt_3prong_Sxy_N1"), jet.pt(), maxSxy);
        registry.fill(HIST("h2_taggedjet_pt_3prong_mass_N1"), jet.pt(), massSV);
      }
    }
    if (fillSVxyz) {
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("h2_jet_pt_3prong_Sxyz_N1"), jet.pt(), maxSxyz);
        registry.fill(HIST("h2_jet_pt_3prong_mass_xyz_N1"), jet.pt(), massSV);
        if (jet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("h2_taggedjet_pt_3prong_Sxyz_N1"), jet.pt(), maxSxyz);
          registry.fill(HIST("h2_taggedjet_pt_3prong_mass_xyz_N1"), jet.pt(), massSV);
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramSV3ProngDataMult(T const& collision, U const& jet, V const& /*prongs*/)
  {
    if (jet.template secondaryVertices_as<V>().size() < 1)
      return;
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h2_jet_pt_mult"), jet.pt(), scaledFT0M);
    registry.fill(HIST("h2_jet_eta_mult"), jet.eta(), scaledFT0M);
    registry.fill(HIST("h2_jet_phi_mult"), jet.phi(), scaledFT0M);
    if (fillGeneralSVQA) {
      registry.fill(HIST("h2_3prong_nprongs_mult"), jet.template secondaryVertices_as<V>().size(), scaledFT0M);
      for (const auto& prong : jet.template secondaryVertices_as<V>()) {
        registry.fill(HIST("hn_jet_3prong_Sxy_mult"), jet.pt(), prong.decayLengthXY(), prong.errorDecayLengthXY(), prong.decayLengthXY() / prong.errorDecayLengthXY(), scaledFT0M);
        if (fillSVxyz) {
          registry.fill(HIST("hn_jet_3prong_Sxyz_mult"), jet.pt(), prong.decayLength(), prong.errorDecayLength(), prong.decayLength() / prong.errorDecayLength(), scaledFT0M);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<V>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("hn_jet_3prong_Sxy_N1_mult"), jet.pt(), maxSxy, massSV, scaledFT0M);
      if (jet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("hn_taggedjet_3prong_Sxy_N1_mult"), jet.pt(), maxSxy, massSV, scaledFT0M);
      }
    }
    if (fillSVxyz) {
      checkSv = false;
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<V>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("hn_jet_3prong_Sxyz_N1_mult"), jet.pt(), maxSxyz, massSV, scaledFT0M);
        if (jet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("hn_taggedjet_3prong_Sxyz_N1_mult"), jet.pt(), maxSxyz, massSV, scaledFT0M);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramSV2ProngMCD(T const& mcdjet, U const& /*prongs*/, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    auto origin = mcdjet.origin();
    if (mcdjet.template secondaryVertices_as<U>().size() < 1)
      return;
    if (!((doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) && (doprocessJPMCD || doprocessJPMCDWeighted || doprocessJPMCPMCDMatched || doprocessJPMCPMCDMatchedWeighted) && (doprocessSV3ProngMCD || doprocessSV3ProngMCDWeighted || doprocessSV3ProngMCPMCDMatched || doprocessSV3ProngMCPMCDMatchedWeighted))) {
      registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), origin, eventWeight);
      registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), origin, eventWeight);
      registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), origin, eventWeight);
    }
    if (fillGeneralSVQA) {
      registry.fill(HIST("h2_2prong_nprongs_flavour"), mcdjet.template secondaryVertices_as<U>().size(), origin, eventWeight);
      for (const auto& prong : mcdjet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("h3_jet_pt_2prong_Lxy_flavour"), mcdjet.pt(), prong.decayLengthXY(), origin, eventWeight);
        registry.fill(HIST("h3_jet_pt_2prong_sigmaLxy_flavour"), mcdjet.pt(), prong.errorDecayLengthXY(), origin, eventWeight);
        registry.fill(HIST("h3_jet_pt_2prong_Sxy_flavour"), mcdjet.pt(), prong.decayLengthXY() / prong.errorDecayLengthXY(), origin, eventWeight);
        if (fillSVxyz) {
          registry.fill(HIST("h3_jet_pt_2prong_Lxyz_flavour"), mcdjet.pt(), prong.decayLength(), origin, eventWeight);
          registry.fill(HIST("h3_jet_pt_2prong_sigmaLxyz_flavour"), mcdjet.pt(), prong.errorDecayLength(), origin, eventWeight);
          registry.fill(HIST("h3_jet_pt_2prong_Sxyz_flavour"), mcdjet.pt(), prong.decayLength() / prong.errorDecayLength(), origin, eventWeight);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("h3_jet_pt_2prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_2prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
      if (mcdjet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("h3_taggedjet_pt_2prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
        registry.fill(HIST("h3_taggedjet_pt_2prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
      }
    }
    if (fillSVxyz) {
      checkSv = false;
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("h3_jet_pt_2prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
        registry.fill(HIST("h3_jet_pt_2prong_mass_xyz_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
        if (mcdjet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("h3_taggedjet_pt_2prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
          registry.fill(HIST("h3_taggedjet_pt_2prong_mass_xyz_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramSV3ProngMCD(T const& mcdjet, U const& /*prongs*/, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    auto origin = mcdjet.origin();
    if (mcdjet.template secondaryVertices_as<U>().size() < 1)
      return;
    if (!((doprocessIPsMCD || doprocessIPsMCDWeighted || doprocessIPsMCPMCDMatched || doprocessIPsMCPMCDMatchedWeighted) && (doprocessJPMCD || doprocessJPMCDWeighted || doprocessJPMCPMCDMatched || doprocessJPMCPMCDMatchedWeighted) && (doprocessSV2ProngMCD || doprocessSV2ProngMCDWeighted || doprocessSV2ProngMCPMCDMatched || doprocessSV2ProngMCPMCDMatchedWeighted))) {
      registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), origin, eventWeight);
      registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), origin, eventWeight);
      registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), origin, eventWeight);
    }
    if (fillGeneralSVQA) {
      registry.fill(HIST("h2_3prong_nprongs_flavour"), mcdjet.template secondaryVertices_as<U>().size(), origin);
      for (const auto& prong : mcdjet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("h3_jet_pt_3prong_Lxy_flavour"), mcdjet.pt(), prong.decayLengthXY(), origin, eventWeight);
        registry.fill(HIST("h3_jet_pt_3prong_sigmaLxy_flavour"), mcdjet.pt(), prong.errorDecayLengthXY(), origin, eventWeight);
        registry.fill(HIST("h3_jet_pt_3prong_Sxy_flavour"), mcdjet.pt(), prong.decayLengthXY() / prong.errorDecayLengthXY(), origin, eventWeight);
        if (fillSVxyz) {
          registry.fill(HIST("h3_jet_pt_3prong_Lxyz_flavour"), mcdjet.pt(), prong.decayLength(), origin, eventWeight);
          registry.fill(HIST("h3_jet_pt_3prong_Sxyz_flavour"), mcdjet.pt(), prong.decayLength() / prong.errorDecayLength(), origin, eventWeight);
          registry.fill(HIST("h3_jet_pt_3prong_sigmaLxyz_flavour"), mcdjet.pt(), prong.errorDecayLength(), origin, eventWeight);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("h3_jet_pt_3prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
      registry.fill(HIST("h3_jet_pt_3prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
      if (mcdjet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("h3_taggedjet_pt_3prong_Sxy_N1_flavour"), mcdjet.pt(), maxSxy, origin, eventWeight);
        registry.fill(HIST("h3_taggedjet_pt_3prong_mass_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
      }
    }
    if (fillSVxyz) {
      checkSv = false;
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("h3_jet_pt_3prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
        registry.fill(HIST("h3_jet_pt_3prong_mass_xyz_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
        if (mcdjet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("h3_taggedjet_pt_3prong_Sxyz_N1_flavour"), mcdjet.pt(), maxSxyz, origin, eventWeight);
          registry.fill(HIST("h3_taggedjet_pt_3prong_mass_xyz_N1_flavour"), mcdjet.pt(), massSV, origin, eventWeight);
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void fillHistogramSV3ProngMCDMult(T const& collision, U const& mcdjet, V const& /*prongs*/, float eventWeight = 1.0)
  {
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    auto origin = mcdjet.origin();
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    if (mcdjet.template secondaryVertices_as<V>().size() < 1)
      return;
    registry.fill(HIST("h3_jet_pt_mult_flavour"), mcdjet.pt(), scaledFT0M, origin, eventWeight);
    registry.fill(HIST("h3_jet_eta_mult_flavour"), mcdjet.eta(), scaledFT0M, origin, eventWeight);
    registry.fill(HIST("h3_jet_phi_mult_flavour"), mcdjet.phi(), scaledFT0M, origin, eventWeight);
    if (fillGeneralSVQA) {
      registry.fill(HIST("h3_3prong_nprongs_mult_flavour"), mcdjet.template secondaryVertices_as<V>().size(), scaledFT0M, origin, eventWeight);
      for (const auto& prong : mcdjet.template secondaryVertices_as<V>()) {
        registry.fill(HIST("hn_jet_3prong_Sxy_mult_flavour"), mcdjet.pt(), prong.decayLengthXY(), prong.errorDecayLengthXY(), prong.decayLengthXY() / prong.errorDecayLengthXY(), scaledFT0M, origin, eventWeight);
        if (fillSVxyz) {
          registry.fill(HIST("hn_jet_3prong_Sxyz_mult_flavour"), mcdjet.pt(), prong.decayLength(), prong.errorDecayLength(), prong.decayLength() / prong.errorDecayLength(), scaledFT0M, origin, eventWeight);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<V>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("hn_jet_3prong_Sxy_N1_mult_flavour"), mcdjet.pt(), maxSxy, massSV, scaledFT0M, origin, eventWeight);
      if (mcdjet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("hn_taggedjet_3prong_Sxy_N1_mult_flavour"), mcdjet.pt(), maxSxy, massSV, scaledFT0M, origin, eventWeight);
      }
    }
    if (fillSVxyz) {
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<V>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("hn_jet_3prong_Sxyz_N1_mult_flavour"), mcdjet.pt(), maxSxyz, massSV, scaledFT0M, origin, eventWeight);
        if (mcdjet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("hn_taggedjet_3prong_Sxyz_N1_mult_flavour"), mcdjet.pt(), maxSxyz, massSV, scaledFT0M, origin, eventWeight);
        }
      }
    }
  }

  void processDummy(aod::Collision const&, aod::Tracks const&)
  {
  }
  PROCESS_SWITCH(JetTaggerHFQA, processDummy, "Dummy process", true);

  void processTracksDca(JetTagTracksData const& tracks)
  {
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      float varImpXY, varImpXYSig, varImpZ, varImpZSig, varImpXYZ, varImpXYZSig;
      varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
      varImpXYSig = track.dcaXY() / track.sigmadcaXY();
      varImpZ = track.dcaZ() * jettaggingutilities::cmTomum;
      varImpZSig = track.dcaZ() / track.sigmadcaZ();
      float dcaXYZ = track.dcaXYZ();
      float sigmadcaXYZ = track.sigmadcaXYZ();
      varImpXYZ = dcaXYZ * jettaggingutilities::cmTomum;
      varImpXYZSig = dcaXYZ / sigmadcaXYZ;

      if (fillIPxy) {
        registry.fill(HIST("h_impact_parameter_xy"), varImpXY);
        registry.fill(HIST("h_impact_parameter_xy_significance"), varImpXYSig);
      }
      if (fillIPz) {
        registry.fill(HIST("h_impact_parameter_z"), varImpZ);
        registry.fill(HIST("h_impact_parameter_z_significance"), varImpZSig);
      }
      if (fillIPxyz) {
        registry.fill(HIST("h_impact_parameter_xyz"), varImpXYZ);
        registry.fill(HIST("h_impact_parameter_xyz_significance"), varImpXYZSig);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processTracksDca, "Fill inclusive tracks' imformation for data", false);

  void processTracksInJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData> const& jets, JetTagTracksData const& /*tracks*/)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      for (auto const& track : jet.template tracks_as<JetTagTracksData>()) {
        float varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        registry.fill(HIST("h2_track_pt_impact_parameter_xy"), track.pt(), varImpXY);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processTracksInJetsData, "Fill QA comtamination of secondary-track inside jets for data jets", false);

  void processSecondaryContaminationMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD> const& mcdjets, JetTagTracksMCD const& /*tracks*/, aod::JetParticles const& /*particles*/)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      float pTHat = 10. / (std::pow(mcdjet.eventWeight(), 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }
      int jetflavour = mcdjet.origin();
      float secondaryPt = 0;
      float totalJetPt = 0;
      for (auto const& track : mcdjet.template tracks_as<JetTagTracksMCD>()) {
        float varImpXY = track.dcaXY() * jettaggingutilities::cmTomum;
        if (!track.has_mcParticle())
          continue;
        auto mcParticle = track.mcParticle();
        totalJetPt += track.pt();
        if (mcParticle.isPhysicalPrimary()) {
          registry.fill(HIST("hn_jet_pt_track_pt_impact_parameter_xy_physical_primary_flavour"), mcdjet.pt(), track.pt(), varImpXY, jetflavour, mcdjet.eventWeight());
        } else {
          registry.fill(HIST("hn_jet_pt_track_pt_impact_parameter_xy_secondary_flavour"), mcdjet.pt(), track.pt(), varImpXY, jetflavour, mcdjet.eventWeight());
          secondaryPt += track.pt();
        }
      }
      if (totalJetPt > 0) {
        float fracSecondary = secondaryPt / totalJetPt;
        registry.fill(HIST("h3_jet_pt_frac_secondary_pt_per_jet_flavour"), mcdjet.pt(), fracSecondary, jetflavour, mcdjet.eventWeight());
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSecondaryContaminationMCD, "Fill QA comtamination of secondary-track inside jets for mcd jets", false);

  void processValFlavourDefMCD(soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, JetTagTracksMCD const& tracks, aod::JetParticles const& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillValidationFlavourDefMCD<soa::Join<JetTableMCP, JetTableMCPMCD>>(mcdjet, tracks, particles, particlesPerColl, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processValFlavourDefMCD, "to check the validation of jet-flavour definition when compared to distance for mcd jets", false);

  void processValFlavourDefMCP(soa::Join<JetTableMCP, weightMCP> const& mcpjets, aod::JetParticles const& particles, aod::JetMcCollisions const&)
  {
    for (auto const& mcpjet : mcpjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, mcpjet.globalIndex());
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(mcpjet)) {
        continue;
      }
      int eventWeight = mcpjet.eventWeight();
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }
      fillValidationFlavourDefMCP(mcpjet, particles, particlesPerColl);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processValFlavourDefMCP, "to check the validation of jet-flavour definition when compared to distance for mcp jets", false);

  void processResponseMatrix(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      auto jetflavour = mcdjet.origin();
      registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), jetflavour);
      registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), jetflavour);
      registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), jetflavour);
      if (!mcdjet.has_matchedJetGeo())
        continue;
      for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<soa::Join<JetTableMCP, JetTableMCPMCD>>()) {
        registry.fill(HIST("h3_jet_pt_jet_pt_part_matchedgeo_flavour"), mcdjet.pt(), mcpjet.pt(), mcdjet.origin());
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processResponseMatrix, "create response matrix with flavour", false);

  void processResponseMatrixWeighted(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs, aod::JCollisionOutliers>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    auto eventWeight = collision.weight();
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }
      auto jetflavour = mcdjet.origin();
      registry.fill(HIST("h2_jet_pt_flavour"), mcdjet.pt(), jetflavour, eventWeight);
      registry.fill(HIST("h2_jet_eta_flavour"), mcdjet.eta(), jetflavour, eventWeight);
      registry.fill(HIST("h2_jet_phi_flavour"), mcdjet.phi(), jetflavour, eventWeight);
      if (!mcdjet.has_matchedJetGeo())
        continue;
      for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<soa::Join<JetTableMCP, JetTableMCPMCD>>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        registry.fill(HIST("h3_jet_pt_jet_pt_part_matchedgeo_flavour"), mcdjet.pt(), mcpjet.pt(), mcdjet.origin(), eventWeight);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processResponseMatrixWeighted, "create response matrix with flavour weighted", false);

  void processMCP(JetTableMCP const& mcpjets, aod::JetParticles const&, soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Filtered<aod::JetCollisionsMCD> const& collisions)
  {
    registry.fill(HIST("h_collision_events"), 2.5); // mcp events
    for (auto const& mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(mcpjet)) {
        continue;
      }
      auto mcCollision = mcCollisions.sliceBy(McCollisionsPerMCPCollision, mcpjet.mcCollisionId());
      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(collisionsPerMCPCollision, mcpjet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          fillHistogramMCP(mcpjet, 1., mcCollision.begin().ptHard());
        }
      } else {
        fillHistogramMCP(mcpjet, 1., mcCollision.begin().ptHard());
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processMCP, "Fill impact parameter imformation for mcp jets", false);

  void processMCPWeighted(JetTableMCP const& mcpjets, aod::JetParticles const&, soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers>> const& collisions)
  {
    registry.fill(HIST("h_collision_events"), 2.5); // mcp events
    for (auto const& mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(mcpjet)) {
        continue;
      }
      auto mcCollision = mcCollisions.sliceBy(McCollisionsPerMCPCollision, mcpjet.mcCollisionId());
      if (mcCollision.size() == 1) {
        if (checkMcCollisionIsMatched) {
          auto collisionspermcpjet = collisions.sliceBy(collisionsPerMCPCollision, mcpjet.mcCollisionId());
          if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits) && !collisionspermcpjet.begin().isOutlier()) {
            fillHistogramMCP(mcpjet, mcCollision.begin().weight(), mcCollision.begin().ptHard());
          }
        } else {
          fillHistogramMCP(mcpjet, mcCollision.begin().weight(), mcCollision.begin().ptHard());
        }
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processMCPWeighted, "Fill impact parameter imformation for mcp jets weighted", false);

  void processIPsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData> const& jets, JetTagTracksData const& tracks)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 0.5); // data events
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramIPsData(jet, tracks);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsData, "Fill impact parameter imformation for data jets", false);

  void processIPsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD> const& mcdjets, JetTagTracksMCD const& tracks, aod::JetParticles const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramIPsMCD(mcdjet, tracks);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCD, "Fill impact parameter imformation for mcd jets", false);

  void processIPsMCDWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::ChargedMCDetectorLevelJetEventWeights> const& mcdjets, JetTagTracksMCD const& tracks, aod::JetParticles const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      fillHistogramIPsMCD(mcdjet, tracks, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCDWeighted, "Fill impact parameter imformation for mcd jets", false);

  void processIPsMCP(JetTableMCP const& mcpjets, aod::JetParticles const&, aod::JetMcCollisions const&, soa::Filtered<aod::JetCollisionsMCD> const& collisions)
  {
    for (auto const& mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(mcpjet)) {
        continue;
      }
      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(collisionsPerMCPCollision, mcpjet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          fillHistogramIPsMCP(mcpjet);
        }
      } else {
        fillHistogramIPsMCP(mcpjet);
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCP, "Fill impact parameter imformation for mcp jets", false);

  void processIPsMCPWeighted(soa::Join<JetTableMCP, aod::ChargedMCParticleLevelJetEventWeights> const& mcpjets, soa::Filtered<aod::JetCollisionsMCD> const& collisions, aod::JetParticles const&)
  {
    for (auto const& mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(mcpjet)) {
        continue;
      }
      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(collisionsPerMCPCollision, mcpjet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          fillHistogramIPsMCP(mcpjet, mcpjet.eventWeight());
        }
      } else {
        fillHistogramIPsMCP(mcpjet, mcpjet.eventWeight());
      }
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCPWeighted, "Fill impact parameter imformation for mcp jets weighted", false);

  void processIPsMCPMCDMatched(soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, JetTagTracksMCD const& tracks, aod::JetParticles const& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo())
        continue;
      for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<soa::Join<JetTableMCP, JetTableMCPMCD>>()) {
        registry.fill(HIST("h3_jet_pt_jet_pt_part_matchedgeo_flavour"), mcdjet.pt(), mcpjet.pt(), mcdjet.origin());
      }
      if (!doprocessIPsMCD)
        fillHistogramIPsMCD(mcdjet, tracks);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCPMCDMatched, "Fill impact parameter imformation for mcp mcd matched jets", false);

  void processIPsMCPMCDMatchedWeighted(soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, JetTagTracksMCD const& tracks, aod::JetParticles const& particles)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo())
        continue;
      float pTHat = 10. / (std::pow(mcdjet.eventWeight(), 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }
      for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<soa::Join<JetTableMCP, JetTableMCPMCD>>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        registry.fill(HIST("h3_jet_pt_jet_pt_part_matchedgeo_flavour"), mcdjet.pt(), mcpjet.pt(), mcdjet.origin(), mcdjet.eventWeight());
      }
      if (!doprocessIPsMCDWeighted)
        fillHistogramIPsMCD(mcdjet, tracks, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processIPsMCPMCDMatchedWeighted, "Fill impact parameter imformation for mcp mcd matched jets", false);

  void processJPData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData> const& jets, JetTagTracksData const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (!doprocessIPsData)
      registry.fill(HIST("h_collision_events"), 0.5); // data events
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramJPData(jet);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPData, "Fill jet probability imformation for data jets", false);

  void processJPMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD> const& mcdjets, JetTagTracksMCD const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramJPMCD(mcdjet);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCD, "Fill jet probability imformation for mcd jets", false);

  void processJPMCDWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD> const& mcdjets, JetTagTracksMCD const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramJPMCD(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCDWeighted, "Fill jet probability imformation for mcd jets", false);

  void processJPMCPMCDMatched(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, JetTagTracksMCD const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo())
        return;
      if (!doprocessJPMCD)
        fillHistogramJPMCD(mcdjet);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCPMCDMatched, "Fill jet probability imformation for mcd jets", false);

  void processJPMCPMCDMatchedWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, JetTagTracksMCD const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo())
        return;
      if (!doprocessJPMCDWeighted)
        fillHistogramJPMCD(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processJPMCPMCDMatchedWeighted, "Fill jet probability imformation for mcd jets", false);

  void processSV2ProngData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex2ProngIndices> const& jets, aod::DataSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 0.5); // mcd events
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramSV2ProngData(jet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngData, "Fill 2prong imformation for data jets", false);

  void processSV3ProngData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 0.5); // mcd events
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramSV3ProngData(jet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngData, "Fill 3prong imformation for data jets", false);

  void processSV3ProngDataMult(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramSV3ProngDataMult(collision, jet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngDataMult, "Fill 3prong imformation for data jets", false);

  void processSV2ProngMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV2ProngMCD(mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCD, "Fill 2prong imformation for mcd jets", false);

  void processSV2ProngMCDWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV2ProngMCD(mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCDWeighted, "Fill 2prong imformation for mcd jets", false);

  void processSV2ProngMCPMCDMatched(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo()) {
        continue;
      }
      if (!doprocessSV2ProngMCD)
        fillHistogramSV2ProngMCD(mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCPMCDMatched, "Fill 2prong imformation for mcd jets", false);

  void processSV2ProngMCPMCDMatchedWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD, aod::MCDSecondaryVertex2ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex2Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo()) {
        continue;
      }
      if (!doprocessSV2ProngMCDWeighted)
        fillHistogramSV2ProngMCD(mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV2ProngMCPMCDMatchedWeighted, "Fill 2prong imformation for mcd jets", false);

  void processSV3ProngMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCD(mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCD, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCDWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCD(mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCDWeighted, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCPMCDMatched(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo()) {
        continue;
      }
      if (!doprocessSV3ProngMCD)
        fillHistogramSV3ProngMCD(mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCPMCDMatched, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCPMCDMatchedWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    registry.fill(HIST("h_collision_events"), 1.5, collision.weight()); // mcd events
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo()) {
        continue;
      }
      if (!doprocessSV3ProngMCDWeighted)
        fillHistogramSV3ProngMCD(mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCPMCDMatchedWeighted, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCDMult(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCDMult(collision, mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCDMult, "Fill 3prong imformation for mcd jets with multiplicity", false);

  void processSV3ProngMCDMultWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCDMult(collision, mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCDMultWeighted, "Fill 3prong imformation for mcd jets with multiplicity weighted", false);

  void processSV3ProngMCPMCDMatchedMult(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo()) {
        continue;
      }
      if (!doprocessSV3ProngMCDMult)
        fillHistogramSV3ProngMCDMult(collision, mcdjet, prongs);
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCPMCDMatchedMult, "Fill 3prong imformation for mcd jets matched with multiplicity", false);

  void processSV3ProngMCPMCDMatchedMultWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, weightMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      if (!mcdjet.has_matchedJetGeo()) {
        continue;
      }
      if (!doprocessSV3ProngMCDMultWeighted)
        fillHistogramSV3ProngMCDMult(collision, mcdjet, prongs, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetTaggerHFQA, processSV3ProngMCPMCDMatchedMultWeighted, "Fill 3prong imformation for mcd jets matched with multiplicity weightd", false);
};

using JetTaggerQAChargedDataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
using JetTaggerQAChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetFlavourDef>;
using JetTaggerQAChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetFlavourDef>;

using JetTaggerhfQaCharged = JetTaggerHFQA<JetTaggerQAChargedDataJets, aod::ChargedJetTags, JetTaggerQAChargedMCDJets, aod::ChargedMCDetectorLevelJetEventWeights, aod::ChargedMCDetectorLevelJetTags, JetTaggerQAChargedMCPJets, aod::ChargedMCParticleLevelJetEventWeights, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetTaggerhfQaCharged>(cfgc, TaskName{"jet-taggerhf-qa-charged"})}; // o2-linter: disable=name/o2-task (templated struct)
}
