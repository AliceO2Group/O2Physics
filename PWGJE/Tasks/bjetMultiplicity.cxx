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

/// \file bjetMultiplicity.cxx
/// \brief bjet multiplicity analysis
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
struct BjetMultTask {

  // task on/off configuration
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
  ConfigurableAxis binMultScaledFT0M{"binMultScaledFT0M", {100, 0, 20.f}, "scaled FT0M"};
  ConfigurableAxis binMultScaledFT0MClass{"binMultScaledFT0MClass", {VARIABLE_WIDTH, 0, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.4, 1.8, 2.4, 3.6, 5., 20.}, "Percentiles of scaled FT0M: 100-90%, 90-80%, 80-70%, 70-60%, 60-50%, 50-40%, 40-30%, 30-20%, 20-10%, 10-1%, 1-0.1%"};

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
    AxisSpec axisMultScaledFT0M = {binMultScaledFT0M, "Multiplicity classes"};
    AxisSpec axisMultScaledFT0MClass = {binMultScaledFT0MClass, "Multiplicity classes"};
    AxisSpec axisCentrality = {110, -5., 105., "Centrality"};

    if (doprocessSV3ProngDataMult) {
      registry.add("h_event_mult", "", {HistType::kTH1F, {{axisMultScaledFT0M}}});
      registry.add("h_event_mult_class", "", {HistType::kTH1F, {{axisMultScaledFT0MClass}}});
      registry.add("h2_jet_pt_mult", "", {HistType::kTH2F, {{axisJetPt}, {axisMultScaledFT0MClass}}});
      registry.add("h2_jet_eta_mult", "", {HistType::kTH2F, {{axisEta}, {axisMultScaledFT0MClass}}});
      registry.add("h2_jet_phi_mult", "", {HistType::kTH2F, {{axisPhi}, {axisMultScaledFT0MClass}}});
      if (fillGeneralSVQA) {
        registry.add("h2_3prong_nprongs_mult", "", {HistType::kTH2F, {{axisNprongs}, {axisMultScaledFT0MClass}}});
        registry.add("hn_jet_3prong_Sxy_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisMultScaledFT0MClass}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisMultScaledFT0MClass}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisMultScaledFT0MClass}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisMultScaledFT0MClass}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisMultScaledFT0MClass}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_mult", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisMultScaledFT0MClass}}});
      }
    }
    if (doprocessSV3ProngMCDMult || doprocessSV3ProngMCDMultWeighted || doprocessSV3ProngMCPMCDMatchedMult || doprocessSV3ProngMCPMCDMatchedMultWeighted) {
      registry.add("h_event_mult", "", {HistType::kTH1F, {{axisMultScaledFT0M}}});
      registry.add("h_event_mult_class", "", {HistType::kTH1F, {{axisMultScaledFT0MClass}}});
      registry.add("h_cent", "mcd collision centrality; centrality; counts", {HistType::kTH1F, {axisCentrality}});
      registry.add("h3_jet_pt_mult_flavour", "", {HistType::kTH3F, {{axisJetPt}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
      registry.add("h3_jet_eta_mult_flavour", "", {HistType::kTH3F, {{axisEta}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
      registry.add("h3_jet_phi_mult_flavour", "", {HistType::kTH3F, {{axisPhi}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
      if (fillGeneralSVQA) {
        registry.add("h3_3prong_nprongs_mult_flavour", "", {HistType::kTH3F, {{axisNprongs}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
        registry.add("hn_jet_3prong_Sxy_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxy}, {axisSigmaLxy}, {axisSxy}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
        if (fillSVxyz) {
          registry.add("hn_jet_3prong_Sxyz_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisLxyz}, {axisSigmaLxyz}, {axisSxyz}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
        }
      }
      registry.add("hn_jet_3prong_Sxy_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisMass}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
      registry.add("hn_taggedjet_3prong_Sxy_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
      if (fillSVxyz) {
        registry.add("hn_jet_3prong_Sxyz_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxyz}, {axisMass}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
        registry.add("hn_taggedjet_3prong_Sxyz_N1_mult_flavour", "", {HistType::kTHnSparseF, {{axisJetPt}, {axisSxy}, {axisSxyz}, {axisMass}, {axisMultScaledFT0MClass}, {axisJetFlavour}}});
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
  PROCESS_SWITCH(BjetMultTask, processDummy, "Dummy process", true);

  void processMCP(JetTableMCP const& mcpjets, aod::JetParticles const&, soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Filtered<aod::JetCollisionsMCD> const& collisions)
  {
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
  PROCESS_SWITCH(BjetMultTask, processMCP, "Fill impact parameter imformation for mcp jets", false);

  void processMCPWeighted(JetTableMCP const& mcpjets, aod::JetParticles const&, soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers>> const& collisions)
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
  PROCESS_SWITCH(BjetMultTask, processMCPWeighted, "Fill impact parameter imformation for mcp jets weighted", false);

  void processSV3ProngDataMult(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    registry.fill(HIST("h_event_mult_class"), scaledFT0M);
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
  PROCESS_SWITCH(BjetMultTask, processSV3ProngDataMult, "Fill 3prong imformation for data jets", false);

//  void processSV3ProngDataMult(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableData, TagTableData, aod::DataSecondaryVertex3ProngIndices> const& jets, aod::DataSecondaryVertex3Prongs const& prongs)
//  {
//    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
//      return;
//    }
//  }
//  PROCESS_SWITCH(BjetMultTask, processSV3ProngDataMult, "Fill 3prong imformation for data jets", false);

//  void processJetsRhoAreaSubMultData()
//  {
//  }
//  PROCESS_SWITCH(BjetMultTask, processJetsRhoAreaSubMultData, "Fill 3prong imformation for mcd jets with multiplicity", false);

  void processSV3ProngMCDMult(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float centFT0A = collision.centFT0A();
    float centFT0C = collision.centFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    std::cout << "FT0A: " << multFT0A << " FT0C: " << multFT0C << " centFT0A: " << centFT0A << " centFT0C: " << centFT0C << std::endl;
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    registry.fill(HIST("h_event_mult_class"), scaledFT0M);
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
  PROCESS_SWITCH(BjetMultTask, processSV3ProngMCDMult, "Fill 3prong imformation for mcd jets with multiplicity", false);

  void processSV3ProngMCDMultWeighted(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    registry.fill(HIST("h_event_mult_class"), scaledFT0M);
    for (auto const& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaCuts->at(0), jetEtaCuts->at(1), trackCuts->at(2), trackCuts->at(3))) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillHistogramSV3ProngMCDMult(collision, mcdjet, prongs, collision.weight());
    }
  }
  PROCESS_SWITCH(BjetMultTask, processSV3ProngMCDMultWeighted, "Fill 3prong imformation for mcd jets with multiplicity weighted", false);

  void processSV3ProngMCPMCDMatchedMult(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    registry.fill(HIST("h_event_mult_class"), scaledFT0M);
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
  PROCESS_SWITCH(BjetMultTask, processSV3ProngMCPMCDMatchedMult, "Fill 3prong imformation for mcd jets matched with multiplicity", false);

  void processSV3ProngMCPMCDMatchedMultWeighted(soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs, aod::JCollisionOutliers>>::iterator const& collision, soa::Join<JetTableMCD, TagTableMCD, JetTableMCDMCP, aod::MCDSecondaryVertex3ProngIndices> const& mcdjets, soa::Join<JetTableMCP, JetTableMCPMCD> const& /*mcpjets*/, aod::MCDSecondaryVertex3Prongs const& prongs)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);
    float cent = collision.centFT0M();
    registry.fill(HIST("h_event_mult"), scaledFT0M);
    registry.fill(HIST("h_event_mult_class"), scaledFT0M);
    registry.fill(HIST("h_cent"), cent);
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
        fillHistogramSV3ProngMCDMult(collision, mcdjet, prongs, collision.weight());
    }
  }
  PROCESS_SWITCH(BjetMultTask, processSV3ProngMCPMCDMatchedMultWeighted, "Fill 3prong imformation for mcd jets matched with multiplicity weightd", false);
};

using BjetMultChargedDataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
using BjetMultChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetFlavourDef>;
using BjetMultChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetFlavourDef>;

using BjetMultCharged = BjetMultTask<BjetMultChargedDataJets, aod::ChargedJetTags, BjetMultChargedMCDJets, aod::ChargedMCDetectorLevelJetEventWeights, aod::ChargedMCDetectorLevelJetTags, BjetMultChargedMCPJets, aod::ChargedMCParticleLevelJetEventWeights, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BjetMultCharged>(cfgc, TaskName{"bjet-multiplicity-charged"})}; // o2-linter: disable=name/o2-task (templated struct)
}
