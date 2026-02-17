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

/// \file bjetCentMult.cxx
/// \brief bjet multiplicity and centrality analysis
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
struct BjetCentMultTask {

  // task on/off configuration
  Configurable<bool> fillGeneralSVQA{"fillGeneralSVQA", true, "process of general QA for sv"};
  Configurable<bool> fillSVxyz{"fillSVxyz", true, "process of decay lenngth of xyz for sv"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

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
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};
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
  ConfigurableAxis binNprongs{"binNprongs", {100, 0., 100.}, ""};
  ConfigurableAxis binLxy{"binLxy", {200, 0, 20.f}, ""};
  ConfigurableAxis binSxy{"binSxy", {1000, 0, 1000.f}, ""};
  ConfigurableAxis binLxyz{"binLxyz", {200, 0, 20.f}, ""};
  ConfigurableAxis binSxyz{"binSxyz", {1000, 0, 1000.f}, ""};
  ConfigurableAxis binMass{"binMass", {50, 0, 10.f}, ""};
  ConfigurableAxis binSigmaLxy{"binSigmaLxy", {100, 0., 0.1}, ""};
  ConfigurableAxis binSigmaLxyz{"binSigmaLxyz", {100, 0., 0.1}, ""};
  ConfigurableAxis binFT0{"binFT0", {1000, 0, 1000.f}, ""};
  ConfigurableAxis binMultScaledFT0M{"binMultScaledFT0M", {100, 0, 20.f}, "scaled FT0M"};
  ConfigurableAxis binMultScaledFT0MClass{"binMultScaledFT0MClass", {VARIABLE_WIDTH, 0, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.4, 1.8, 2.4, 3.6, 5., 20.}, "Percentiles of scaled FT0M: 100-90%, 90-80%, 80-70%, 70-60%, 60-50%, 50-40%, 40-30%, 30-20%, 20-10%, 10-1%, 1-0.1%"};
  ConfigurableAxis rebinnedCent{"rebinnedCent", {VARIABLE_WIDTH, 0, 0.1, 1.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100., 110.}, "Rebinned percentiles of centality: 100-90%, 90-80%, 80-70%, 70-60%, 60-50%, 50-40%, 40-30%, 30-20%, 20-10%, 10-1%, 1-0.1%"};

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
    AxisSpec axisNprongs = {binNprongs, "#it{N}_{SV}"};
    AxisSpec axisLxy = {binLxy, "L_{XY} [cm]"};
    AxisSpec axisSxy = {binSxy, "S_{XY}"};
    AxisSpec axisLxyz = {binLxyz, "L_{XYZ} [cm]"};
    AxisSpec axisSxyz = {binSxyz, "S_{XYZ}"};
    AxisSpec axisMass = {binMass, "#it{m}_{SV}"};
    AxisSpec axisSigmaLxy = {binSigmaLxy, "#sigma_{L_{XY}} [cm]"};
    AxisSpec axisSigmaLxyz = {binSigmaLxyz, "#sigma_{L_{XYZ}} [cm]"};
    AxisSpec axisFracSecPt = {100, 0, 1, "#frac{#Sigma#it{p}_{T}^{secondary track}}{#it{p}_{T, jet}}"};
    AxisSpec axisFT0 = {binFT0, "amplitude FT0"};
    AxisSpec axisMultScaledFT0M = {binMultScaledFT0M, "scaled FT0M"};
    AxisSpec axisMultScaledFT0MClass = {binMultScaledFT0MClass, "scaled FT0M classes"};
    AxisSpec axisCentrality = {110, -5., 105., "Centrality"};
    AxisSpec axisRebinnedCentrality = {rebinnedCent, "Centrality"};
    AxisSpec axisPercentileMultiplicity = {110, -5., 105., "Percentile multiplicity"};

    if (doprocessCentMultQa) {
      registry.add("h_amplitude_FT0A", "", {HistType::kTH1F, {{axisFT0}}});
      registry.add("h_amplitude_FT0C", "", {HistType::kTH1F, {{axisFT0}}});
      registry.add("h_scaled_FT0M", "", {HistType::kTH1F, {{axisMultScaledFT0M}}});
      registry.add("h_scaled_FT0M_class", "", {HistType::kTH1F, {{axisMultScaledFT0MClass}}});
      registry.add("h2_centrality_percentile_multiplicity", "mcd collision centrality; centrality; counts", {HistType::kTH2F, {{axisRebinnedCentrality}, {axisPercentileMultiplicity}}});
    }
    if (doprocessSV3ProngData) {
      registry.add("h_event_centrality", "", {HistType::kTH1F, {{axisCentrality}}});
      registry.add("h2_jet_pt_centrality", "", {HistType::kTH2F, {{axisJetPt}, {axisCentrality}}});
      registry.add("h2_jet_eta_centrality", "", {HistType::kTH2F, {{axisEta}, {axisCentrality}}});
      registry.add("h2_jet_phi_centrality", "", {HistType::kTH2F, {{axisPhi}, {axisCentrality}}});
      if (fillGeneralSVQA) {
        registry.add("h2_3prong_nprongs_centrality", "", {HistType::kTH2F, {{axisNprongs}, {axisCentrality}}});
        registry.add("hn_jet_3prong_Sxy_centrality", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisCentrality}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_centrality", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisCentrality}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_centrality", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisCentrality}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_centrality", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisCentrality}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_centrality", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisCentrality}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_centrality", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisCentrality}}});
      }
    }
    if (doprocessSV3ProngMCD || doprocessSV3ProngMCPMCDMatched) {
      registry.add("h_event_centrality", "", {HistType::kTH1F, {{axisCentrality}}});
      registry.add("h3_jet_pt_centrality_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisCentrality}, {axisJetFlavour}}});
      registry.add("h3_jet_eta_centrality_flavour", "", {HistType::kTH3F, {{axisEta}, {axisCentrality}, {axisJetFlavour}}});
      registry.add("h3_jet_phi_centrality_flavour", "", {HistType::kTH3F, {{axisPhi}, {axisCentrality}, {axisJetFlavour}}});
      if (fillGeneralSVQA) {
        registry.add("h3_3prong_nprongs_centrality_flavour", "", {HistType::kTH3F, {{axisNprongs}, {axisCentrality}, {axisJetFlavour}}});
        registry.add("hn_jet_3prong_Sxy_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisCentrality}, {axisJetFlavour}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisCentrality}, {axisJetFlavour}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
      }
    }
    if (doprocessRhoAreaSubSV3ProngMCD || doprocessRhoAreaSubSV3ProngMCPMCDMatched) {
      registry.add("h_event_rhoareasubtracted_centrality", "", {HistType::kTH1F, {{axisCentrality}}});
      registry.add("h3_jet_pt_rhoareasubtracted_centrality_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisCentrality}, {axisJetFlavour}}});
      registry.add("h3_jet_eta_rhoareasubtracted_centrality_flavour", "", {HistType::kTH3F, {{axisEta}, {axisCentrality}, {axisJetFlavour}}});
      registry.add("h3_jet_phi_rhoareasubtracted_centrality_flavour", "", {HistType::kTH3F, {{axisPhi}, {axisCentrality}, {axisJetFlavour}}});
      if (fillGeneralSVQA) {
        registry.add("h3_3prong_nprongs_rhoareasubtracted_centrality_flavour", "", {HistType::kTH3F, {{axisNprongs}, {axisCentrality}, {axisJetFlavour}}});
        registry.add("hn_jet_3prong_Sxy_rhoareasubtracted_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisCentrality}, {axisJetFlavour}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_rhoareasubtracted_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisCentrality}, {axisJetFlavour}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_rhoareasubtracted_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_rhoareasubtracted_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_rhoareasubtracted_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_rhoareasubtracted_centrality_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisCentrality}, {axisJetFlavour}}});
      }
    }
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackCuts->at(0) && aod::jtrack::pt < trackCuts->at(1) && aod::jtrack::eta > trackCuts->at(2) && aod::jtrack::eta < trackCuts->at(3));
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>> McCollisionsPerMCPCollision = aod::jmccollision::mcCollisionId;
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

  float getPercentile(const float scaledFT0M)
  {
    static const std::vector<float> x = {0.0, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.4, 3.6, 5.0, 20.0};
    static const std::vector<float> y = {100.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 1.0, 0.1, 0.0};

    if (scaledFT0M <= x.front())
      return y.front();
    if (scaledFT0M >= x.back())
      return y.back();

    for (size_t i = 0; i < x.size() - 1; ++i) {
      if (scaledFT0M >= x[i] && scaledFT0M < x[i + 1]) {
        float slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        return y[i] + slope * (scaledFT0M - x[i]);
      }
    }
    return -1.0f;
  }

  template <typename T>
  void fillHistogramMCP(T const& mcpjet, float eventWeight = 1.0, float pTHat = 999.)
  {
    if (mcpjet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
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

  template <typename T, typename U>
  void fillHistogramSV3ProngData(T const& jet, U const& /*prongs*/, float centrality)
  {
    if (jet.template secondaryVertices_as<U>().size() < 1)
      return;
    registry.fill(HIST("h2_jet_pt_centrality"), jet.pt(), centrality);
    registry.fill(HIST("h2_jet_eta_centrality"), jet.eta(), centrality);
    registry.fill(HIST("h2_jet_phi_centrality"), jet.phi(), centrality);
    if (fillGeneralSVQA) {
      registry.fill(HIST("h2_3prong_nprongs_centrality"), jet.template secondaryVertices_as<U>().size(), centrality);
      for (const auto& prong : jet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("hn_jet_3prong_Sxy_centrality"), jet.pt(), prong.decayLengthXY(), prong.errorDecayLengthXY(), prong.decayLengthXY() / prong.errorDecayLengthXY(), centrality);
        if (fillSVxyz) {
          registry.fill(HIST("hn_jet_3prong_Sxyz_centrality"), jet.pt(), prong.decayLength(), prong.errorDecayLength(), prong.decayLength() / prong.errorDecayLength(), centrality);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("hn_jet_3prong_Sxy_N1_centrality"), jet.pt(), maxSxy, massSV, centrality);
      if (jet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("hn_taggedjet_3prong_Sxy_N1_centrality"), jet.pt(), maxSxy, massSV, centrality);
      }
    }
    if (fillSVxyz) {
      checkSv = false;
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(jet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("hn_jet_3prong_Sxyz_N1_centrality"), jet.pt(), maxSxyz, massSV, centrality);
        if (jet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("hn_taggedjet_3prong_Sxyz_N1_centrality"), jet.pt(), maxSxyz, massSV, centrality);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillHistogramSV3ProngMCD(T const& mcdjet, U const& /*prongs*/, float centrality, float eventWeight = 1.0)
  {
    if (useEventWeight) {
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }
    }
    auto origin = mcdjet.origin();
    if (mcdjet.template secondaryVertices_as<U>().size() < 1)
      return;
    registry.fill(HIST("h3_jet_pt_centrality_flavour"), mcdjet.pt(), centrality, origin, eventWeight);
    registry.fill(HIST("h3_jet_eta_centrality_flavour"), mcdjet.eta(), centrality, origin, eventWeight);
    registry.fill(HIST("h3_jet_phi_centrality_flavour"), mcdjet.phi(), centrality, origin, eventWeight);
    if (fillGeneralSVQA) {
      registry.fill(HIST("h3_3prong_nprongs_centrality_flavour"), mcdjet.template secondaryVertices_as<U>().size(), centrality, origin, eventWeight);
      for (const auto& prong : mcdjet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("hn_jet_3prong_Sxy_centrality_flavour"), mcdjet.pt(), prong.decayLengthXY(), prong.errorDecayLengthXY(), prong.decayLengthXY() / prong.errorDecayLengthXY(), centrality, origin, eventWeight);
        if (fillSVxyz) {
          registry.fill(HIST("hn_jet_3prong_Sxyz_centrality_flavour"), mcdjet.pt(), prong.decayLength(), prong.errorDecayLength(), prong.decayLength() / prong.errorDecayLength(), centrality, origin, eventWeight);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("hn_jet_3prong_Sxy_N1_centrality_flavour"), mcdjet.pt(), maxSxy, massSV, centrality, origin, eventWeight);
      if (mcdjet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("hn_taggedjet_3prong_Sxy_N1_centrality_flavour"), mcdjet.pt(), maxSxy, massSV, centrality, origin, eventWeight);
      }
    }
    if (fillSVxyz) {
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("hn_jet_3prong_Sxyz_N1_centrality_flavour"), mcdjet.pt(), maxSxyz, massSV, centrality, origin, eventWeight);
        if (mcdjet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("hn_taggedjet_3prong_Sxyz_N1_centrality_flavour"), mcdjet.pt(), maxSxyz, massSV, centrality, origin, eventWeight);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillRhoAreaSubtractedHistogramSV3ProngMCD(T const& mcdjet, U const& /*prongs*/, float centrality, float rho, float eventWeight = 1.0)
  {
    if (useEventWeight) {
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }
    }
    auto origin = mcdjet.origin();
    if (mcdjet.template secondaryVertices_as<U>().size() < 1)
      return;
    registry.fill(HIST("h3_jet_pt_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), centrality, origin, eventWeight);
    registry.fill(HIST("h3_jet_eta_rhoareasubtracted_centrality_flavour"), mcdjet.eta(), centrality, origin, eventWeight);
    registry.fill(HIST("h3_jet_phi_rhoareasubtracted_centrality_flavour"), mcdjet.phi(), centrality, origin, eventWeight);
    if (fillGeneralSVQA) {
      registry.fill(HIST("h3_3prong_nprongs_rhoareasubtracted_centrality_flavour"), mcdjet.template secondaryVertices_as<U>().size(), centrality, origin, eventWeight);
      for (const auto& prong : mcdjet.template secondaryVertices_as<U>()) {
        registry.fill(HIST("hn_jet_3prong_Sxy_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), prong.decayLengthXY(), prong.errorDecayLengthXY(), prong.decayLengthXY() / prong.errorDecayLengthXY(), centrality, origin, eventWeight);
        if (fillSVxyz) {
          registry.fill(HIST("hn_jet_3prong_Sxyz_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), prong.decayLength(), prong.errorDecayLength(), prong.decayLength() / prong.errorDecayLength(), centrality, origin, eventWeight);
        }
      }
    }
    bool checkSv = false;
    auto bjetCand = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(2), prongCuts->at(4), prongCuts->at(5), false, &checkSv);
    if (checkSv && jettaggingutilities::svAcceptance(bjetCand, svDispersionMax)) {
      auto maxSxy = bjetCand.decayLengthXY() / bjetCand.errorDecayLengthXY();
      auto massSV = bjetCand.m();
      registry.fill(HIST("hn_jet_3prong_Sxy_N1_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), maxSxy, massSV, centrality, origin, eventWeight);
      if (mcdjet.isTagged(BJetTaggingMethod::SV)) {
        registry.fill(HIST("hn_taggedjet_3prong_Sxy_N1_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), maxSxy, massSV, centrality, origin, eventWeight);
      }
    }
    if (fillSVxyz) {
      auto bjetCandXYZ = jettaggingutilities::jetFromProngMaxDecayLength<U>(mcdjet, prongCuts->at(0), prongCuts->at(1), prongCuts->at(3), prongCuts->at(4), prongCuts->at(5), true, &checkSv);
      if (checkSv && jettaggingutilities::svAcceptance(bjetCandXYZ, svDispersionMax)) {
        auto maxSxyz = bjetCandXYZ.decayLength() / bjetCandXYZ.errorDecayLength();
        auto massSV = bjetCandXYZ.m();
        registry.fill(HIST("hn_jet_3prong_Sxyz_N1_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), maxSxyz, massSV, centrality, origin, eventWeight);
        if (mcdjet.isTagged(BJetTaggingMethod::SV3D)) {
          registry.fill(HIST("hn_taggedjet_3prong_Sxyz_N1_rhoareasubtracted_centrality_flavour"), mcdjet.pt() - (rho * mcdjet.area()), maxSxyz, massSV, centrality, origin, eventWeight);
        }
      }
    }
  }

  void processDummy(aod::Collision const&, aod::Tracks const&)
  {
  }
  PROCESS_SWITCH(BjetCentMultTask, processDummy, "Dummy process", true);

  void processCentMultQa(soa::Filtered<aod::JetCollisions>::iterator const& collision)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    float percentileMult = getPercentile(scaledFT0M);
    registry.fill(HIST("h_amplitude_FT0A"), multFT0A);
    registry.fill(HIST("h_amplitude_FT0C"), multFT0C);
    registry.fill(HIST("h_scaled_FT0M"), scaledFT0M);
    registry.fill(HIST("h_scaled_FT0M_class"), scaledFT0M);
    registry.fill(HIST("h2_centrality_percentile_multiplicity"), collision.centFT0M(), percentileMult);
  }
  PROCESS_SWITCH(BjetCentMultTask, processCentMultQa, "Fill centality and mulitplicity qa", false);

  void processMCP(JetTableMCP const& mcpjets, aod::JetParticles const&, soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers>> const& collisions)
  {
    for (auto const& mcpjet : mcpjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(mcpjet)) {
        continue;
      }
      auto mcCollision = mcCollisions.sliceBy(McCollisionsPerMCPCollision, mcpjet.mcCollisionId());
      if (mcCollision.size() == 1) {
        float weight = useEventWeight ? mcCollision.begin().weight() : 1.f;
        if (checkMcCollisionIsMatched) {
          auto collisionspermcpjet = collisions.sliceBy(collisionsPerMCPCollision, mcpjet.mcCollisionId());
          if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits) && !collisionspermcpjet.begin().isOutlier()) {
            fillHistogramMCP(mcpjet, weight, mcCollision.begin().ptHard());
          }
        } else {
          fillHistogramMCP(mcpjet, weight, mcCollision.begin().ptHard());
        }
      }
    }
  }
  PROCESS_SWITCH(BjetCentMultTask, processMCP, "Fill impact parameter imformation for mcp jets", false);

  void processSV3ProngData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float centrality = collision.centFT0M();
    registry.fill(HIST("h_event_centrality"), centrality);
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramSV3ProngData(jet, prongs, centrality);
    }
  }
  PROCESS_SWITCH(BjetCentMultTask, processSV3ProngData, "Fill 3prong imformation for data jets", false);

  void processSV3ProngMCD(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    float centrality = collision.centFT0M();
    float weight = useEventWeight ? collision.weight() : 1.f;
    registry.fill(HIST("h_event_centrality"), collision.centFT0M());
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCD(mcdjet, prongs, centrality, weight);
    }
  }
  PROCESS_SWITCH(BjetCentMultTask, processSV3ProngMCD, "Fill 3prong imformation for mcd jets", false);

  void processSV3ProngMCPMCDMatched(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs, aod::JCollisionOutliers>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    float centrality = collision.centFT0M();
    float weight = useEventWeight ? collision.weight() : 1.f;
    registry.fill(HIST("h_event_centrality"), collision.centFT0M());
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
      fillHistogramSV3ProngMCD(mcdjet, prongs, centrality, weight);
    }
  }
  PROCESS_SWITCH(BjetCentMultTask, processSV3ProngMCPMCDMatched, "Fill 3prong imformation for mcd jets matched", false);

  void processRhoAreaSubSV3ProngMCD(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    float centrality = collision.centFT0M();
    float rho = collision.rho();
    float weight = useEventWeight ? collision.weight() : 1.f;
    registry.fill(HIST("h_event_rhoareasubtracted_centrality"), collision.centFT0M());
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillRhoAreaSubtractedHistogramSV3ProngMCD(mcdjet, prongs, centrality, rho, weight);
    }
  }
  PROCESS_SWITCH(BjetCentMultTask, processRhoAreaSubSV3ProngMCD, "Fill 3prong imformation for mcd jets with background subraction", false);

  void processRhoAreaSubSV3ProngMCPMCDMatched(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs, aod::JCollisionOutliers, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    float centrality = collision.centFT0M();
    float rho = collision.rho();
    float weight = useEventWeight ? collision.weight() : 1.f;
    registry.fill(HIST("h_event_rhoareasubtracted_centrality"), collision.centFT0M());
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
      fillRhoAreaSubtractedHistogramSV3ProngMCD(mcdjet, prongs, centrality, rho, weight);
    }
  }
  PROCESS_SWITCH(BjetCentMultTask, processRhoAreaSubSV3ProngMCPMCDMatched, "Fill 3prong imformation for mcd jets matched with background subtraction", false);
};

using BjetCentMultChargedDataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
using BjetCentMultChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetFlavourDef>;
using BjetCentMultChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetFlavourDef>;

using BjetCentMultCharged = BjetCentMultTask<BjetCentMultChargedDataJets, aod::ChargedJetTags, BjetCentMultChargedMCDJets, aod::ChargedMCDetectorLevelJetEventWeights, aod::ChargedMCDetectorLevelJetTags, BjetCentMultChargedMCPJets, aod::ChargedMCParticleLevelJetEventWeights, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BjetCentMultCharged>(cfgc, TaskName{"bjet-cent-mult-charged"})}; // o2-linter: disable=name/o2-task (templated struct)
}
