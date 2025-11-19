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

/// \file bjetTaggingML.cxx
/// \brief Task for tagging the beauty jets using ML algorithm (in onnx format) loaded from ccdb
///
/// \author Hadi Hassan <hadi.hassan@cern.ch>, University of Jyväskylä

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BJetTaggingML {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.5, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> maxIPxy{"maxIPxy", 10, "maximum track DCA in xy plane"};
  Configurable<float> maxIPz{"maxIPz", 10, "maximum track DCA in z direction"};

  // track level configurables
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> nJetConst{"nJetConst", 10, "maximum number of jet consistuents to be used for ML evaluation"};
  Configurable<bool> useDb{"useDb", false, "Flag whether to use the Db instead of the score for tagging"};
  Configurable<bool> callSumw2{"callSumw2", false, "Flag whether to use sum weight squared for THnSparse error calculation"};

  Configurable<bool> doDataDriven{"doDataDriven", false, "Flag whether to use fill THnSpase for data driven methods"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  std::vector<int> eventSelectionBits;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    // Seed the random number generator using current time
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    const AxisSpec axisDb{300, -10., 20., "#it{D}_{b}"};
    const AxisSpec axisScore{120, -0.1, 1.1, "Score"};
    const AxisSpec axisLogScore{120, 0., 30, "-log(1-score)"};

    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{40, -20.0, 20.0}}});

    registry.add("h2_score_jetpT", "ML scores for inclusive jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, useDb ? axisDb : axisScore}});
    registry.add("h2_logscore_jetpT", "ML scores for inclusive jets;#it{p}_{T,jet} (GeV/#it{c});- log(1 - Score)", {HistType::kTH2F, {{200, 0., 200.}, axisLogScore}});

    registry.add("h2_nTracks_jetpT", "Number of tracks;#it{p}_{T,jet} (GeV/#it{c});nTracks", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 100.0}}});
    registry.add("h2_nSV_jetpT", "Number of secondary vertices;#it{p}_{T,jet} (GeV/#it{c});nSVs", {HistType::kTH2F, {{200, 0., 200.}, {250, 0, 250.0}}});

    registry.add("h2_SIPs2D_jetpT", "2D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_SIPs3D_jetpT", "3D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_LxyS_jetpT", "Decay length in XY;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
    registry.add("h2_Dispersion_jetpT", "SV dispersion;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
    registry.add("h2_jetMass_jetpT", "Jet mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
    registry.add("h2_SVMass_jetpT", "Secondary vertex mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10}}});

    if (doDataDriven) {
      registry.add("hSparse_Incljets", "Inclusive jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;-log(1-score);#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, useDb ? axisDb : axisScore, axisLogScore, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}}, callSumw2);
      if (doprocessMCJets || doprocessMCJetsWeighted) {
        registry.add("hSparse_bjets", "Tagged b-jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;-log(1-score);#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, useDb ? axisDb : axisScore, axisLogScore, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}}, callSumw2);
        registry.add("hSparse_cjets", "Tagged c-jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;-log(1-score);#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, useDb ? axisDb : axisScore, axisLogScore, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}}, callSumw2);
        registry.add("hSparse_lfjets", "Tagged lf-jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;-log(1-score);#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, useDb ? axisDb : axisScore, axisLogScore, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}}, callSumw2);
      }
    }

    if (doprocessMCJets || doprocessMCJetsWeighted || doprocessMCTruthJets || doprocessMCTruthJetsWeighted) {

      registry.add("h2_score_jetpT_bjet", "ML scores for b-jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, useDb ? axisDb : axisScore}});
      registry.add("h2_logscore_jetpT_bjet", "ML scores for b-jets;#it{p}_{T,jet} (GeV/#it{c});- log(1 - Score)", {HistType::kTH2F, {{200, 0., 200.}, axisLogScore}});
      registry.add("h2_SIPs2D_jetpT_bjet", "2D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_bjet", "3D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_bjet", "Decay length in XY b-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_bjet", "SV dispersion b-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_bjet", "Jet mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_bjet", "Secondary vertex mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_score_jetpT_cjet", "ML scores for c-jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, useDb ? axisDb : axisScore}});
      registry.add("h2_logscore_jetpT_cjet", "ML scores for c-jets;#it{p}_{T,jet} (GeV/#it{c});- log(1 - Score)", {HistType::kTH2F, {{200, 0., 200.}, axisLogScore}});
      registry.add("h2_SIPs2D_jetpT_cjet", "2D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_cjet", "3D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_cjet", "Decay length in XY c-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_cjet", "SV dispersion c-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_cjet", "Jet mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_cjet", "Secondary vertex mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_score_jetpT_lfjet", "ML scores for lf-jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, useDb ? axisDb : axisScore}});
      registry.add("h2_logscore_jetpT_lfjet", "ML scores for lf-jets;#it{p}_{T,jet} (GeV/#it{c});- log(1 - Score)", {HistType::kTH2F, {{200, 0., 200.}, axisLogScore}});
      registry.add("h2_SIPs2D_jetpT_lfjet", "2D IP significance lf-jet;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_lfjet", "3D IP significance lf-jet;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_lfjet", "Decay length in XY lf-jet;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_lfjet", "SV dispersion lf-jet;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_lfjet", "Jet mass lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_lfjet", "Secondary vertex mass lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h_jetpT_detector_bjet", "Jet transverse momentum b-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_detector_cjet", "Jet transverse momentum c-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_detector_lfjet", "Jet transverse momentum lf-jet;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h_jetpT_particle_DetColl", "Jet transverse momentum particle level inclusive jets (Detector-level collisions);#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_DetColl_bjet", "Jet transverse momentum particle level b-jets (Detector-level collisions);#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_DetColl_cjet", "Jet transverse momentum particle level c-jets (Detector-level collisions);#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_DetColl_lfjet", "Jet transverse momentum particle level lf-jet (Detector-level collisions);#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h_jetpT_particle_bjet", "Jet transverse momentum particle level b-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_cjet", "Jet transverse momentum particle level c-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_lfjet", "Jet transverse momentum particle level lf-jet;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h2_Response_DetjetpT_PartjetpT_bjet", "Response matrix b-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_cjet", "Response matrix c-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_lfjet", "Response matrix lf-jet;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}, callSumw2);
    }
  }

  // FIXME filtering only works when you loop directly over the list, but if you loop over it as a constituent they will not be filtered
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt > trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt <= jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using FilteredCollision = soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs>>;
  using JetTrackswID = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using JetTracksMCDwID = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  using DataJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::DataSecondaryVertex3ProngIndices, aod::ChargedJetTags>>;

  // Looping over the SV info and writing them to a table
  template <typename AnalysisJet, typename SecondaryVertices>
  void analyzeJetSVInfo(AnalysisJet const& myJet, SecondaryVertices const& /*allSVs*/, std::vector<jettaggingutilities::BJetSVParams>& svsParams, int8_t jetFlavor = 0, double eventweight = 1.0)
  {
    using SVType = typename SecondaryVertices::iterator;

    // Min-heap to store the top 30 SVs by decayLengthXY/errorDecayLengthXY
    auto compare = [](SVType& sv1, SVType& sv2) {
      return (sv1.decayLengthXY() / sv1.errorDecayLengthXY()) > (sv2.decayLengthXY() / sv2.errorDecayLengthXY());
    };

    auto svs = myJet.template secondaryVertices_as<SecondaryVertices>();

    // Sort the SVs based on their decay length significance in descending order
    // This is needed in order to select longest SVs since some jets could have thousands of SVs
    std::sort(svs.begin(), svs.end(), compare);

    for (const auto& candSV : svs) {

      if (candSV.pt() < svPtMin) {
        continue;
      }

      double deltaRJetSV = jetutilities::deltaR(myJet, candSV);
      double massSV = candSV.m();
      double energySV = candSV.e();

      if (static_cast<int>(svsParams.size()) < nJetConst) {
        svsParams.emplace_back(jettaggingutilities::BJetSVParams{candSV.pt(), deltaRJetSV, massSV, energySV / myJet.energy(), candSV.impactParameterXY(), candSV.cpa(), candSV.chi2PCA(), candSV.dispersion(), candSV.decayLengthXY(), candSV.errorDecayLengthXY(), candSV.decayLength(), candSV.errorDecayLength()});
      }

      registry.fill(HIST("h2_LxyS_jetpT"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
      registry.fill(HIST("h2_Dispersion_jetpT"), myJet.pt(), candSV.dispersion(), eventweight);
      registry.fill(HIST("h2_SVMass_jetpT"), myJet.pt(), massSV, eventweight);

      if (doprocessMCJets || doprocessMCJetsWeighted) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_LxyS_jetpT_bjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_bjet"), myJet.pt(), candSV.dispersion(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_bjet"), myJet.pt(), massSV, eventweight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_LxyS_jetpT_cjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_cjet"), myJet.pt(), candSV.dispersion(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_cjet"), myJet.pt(), massSV, eventweight);
        } else {
          registry.fill(HIST("h2_LxyS_jetpT_lfjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_lfjet"), myJet.pt(), candSV.dispersion(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_lfjet"), myJet.pt(), massSV, eventweight);
        }
      }
    }
    svsParams.resize(nJetConst);
  }

  template <typename AnalysisJet, typename AnyTracks>
  void analyzeJetTrackInfo(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, std::vector<jettaggingutilities::BJetTrackParams>& tracksParams, int8_t jetFlavor = 0, double eventweight = 1.0)
  {

    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin || !jettaggingutilities::trackAcceptanceWithDca(constituent, maxIPxy, maxIPz)) {
        continue;
      }

      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

      float rClosestSV = -1.0;

      registry.fill(HIST("h2_SIPs2D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
      registry.fill(HIST("h2_SIPs3D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);

      if (doprocessMCJets || doprocessMCJetsWeighted) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_SIPs2D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_SIPs2D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else {
          registry.fill(HIST("h2_SIPs2D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        }
      }

      tracksParams.emplace_back(jettaggingutilities::BJetTrackParams{constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaXYZ()) * sign, constituent.sigmadcaXYZ(), constituent.p() / analysisJet.p(), rClosestSV});
    }

    auto compare = [](const jettaggingutilities::BJetTrackParams& tr1, const jettaggingutilities::BJetTrackParams& tr2) {
      return (tr1.signedIP2D / tr1.signedIP2DSign) > (tr2.signedIP2D / tr2.signedIP2DSign);
    };

    // Sort the tracks based on their IP significance in descending order
    std::sort(tracksParams.begin(), tracksParams.end(), compare);
    tracksParams.resize(nJetConst);
  }

  template <typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
  void processJetInfo(AnalysisJet const& myJet, AnyTracks const& allTracks, SecondaryVertices const& allSVs, int8_t jetFlavor = 0, double eventWeight = 1.0)
  {
    std::vector<jettaggingutilities::BJetTrackParams> tracksParams;
    std::vector<jettaggingutilities::BJetSVParams> svsParams;

    analyzeJetSVInfo(myJet, allSVs, svsParams, jetFlavor, eventWeight);
    analyzeJetTrackInfo(myJet, allTracks, tracksParams, jetFlavor, eventWeight);

    int nSVs = myJet.template secondaryVertices_as<SecondaryVertices>().size();
    int nTracks = myJet.template tracks_as<AnyTracks>().size();

    registry.fill(HIST("h2_nTracks_jetpT"), myJet.pt(), nTracks);
    registry.fill(HIST("h2_nSV_jetpT"), myJet.pt(), nSVs < 250 ? nSVs : 249);

    registry.fill(HIST("h2_score_jetpT"), myJet.pt(), myJet.scoreML(), eventWeight);
    if (!useDb) {
      registry.fill(HIST("h2_logscore_jetpT"), myJet.pt(), -1 * std::log(1 - myJet.scoreML()), eventWeight);
    }
    registry.fill(HIST("h2_jetMass_jetpT"), myJet.pt(), myJet.mass(), eventWeight);

    if (doDataDriven) {
      registry.fill(HIST("hSparse_Incljets"), myJet.pt(), myJet.scoreML(), useDb ? 0 : -1 * std::log(1 - myJet.scoreML()), myJet.mass(), -1 * std::log(myJet.jetProb()), svsParams[0].svMass, svsParams[0].svfE, eventWeight);
      if (doprocessMCJets || doprocessMCJetsWeighted) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("hSparse_bjets"), myJet.pt(), myJet.scoreML(), useDb ? 0 : -1 * std::log(1 - myJet.scoreML()), myJet.mass(), -1 * std::log(myJet.jetProb()), svsParams[0].svMass, svsParams[0].svfE, eventWeight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("hSparse_cjets"), myJet.pt(), myJet.scoreML(), useDb ? 0 : -1 * std::log(1 - myJet.scoreML()), myJet.mass(), -1 * std::log(myJet.jetProb()), svsParams[0].svMass, svsParams[0].svfE, eventWeight);
        } else {
          registry.fill(HIST("hSparse_lfjets"), myJet.pt(), myJet.scoreML(), useDb ? 0 : -1 * std::log(1 - myJet.scoreML()), myJet.mass(), -1 * std::log(myJet.jetProb()), svsParams[0].svMass, svsParams[0].svfE, eventWeight);
        }
      }
    }

    if (doprocessMCJets || doprocessMCJetsWeighted) {
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_score_jetpT_bjet"), myJet.pt(), myJet.scoreML(), eventWeight);
        if (!useDb) {
          registry.fill(HIST("h2_logscore_jetpT_bjet"), myJet.pt(), -1 * std::log(1 - myJet.scoreML()), eventWeight);
        }
        registry.fill(HIST("h2_jetMass_jetpT_bjet"), myJet.pt(), myJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_bjet"), myJet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_score_jetpT_cjet"), myJet.pt(), myJet.scoreML(), eventWeight);
        if (!useDb) {
          registry.fill(HIST("h2_logscore_jetpT_cjet"), myJet.pt(), -1 * std::log(1 - myJet.scoreML()), eventWeight);
        }
        registry.fill(HIST("h2_jetMass_jetpT_cjet"), myJet.pt(), myJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_cjet"), myJet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h2_score_jetpT_lfjet"), myJet.pt(), myJet.scoreML(), eventWeight);
        if (!useDb) {
          registry.fill(HIST("h2_logscore_jetpT_lfjet"), myJet.pt(), -1 * std::log(1 - myJet.scoreML()), eventWeight);
        }
        registry.fill(HIST("h2_jetMass_jetpT_lfjet"), myJet.pt(), myJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_lfjet"), myJet.pt(), eventWeight);
      }
    }
  }

  template <typename AnalysisJetMCP>
  void fillMCPHistograms(AnalysisJetMCP const& mcpjet, bool detColl = false, double eventWeight = 1.0)
  {
    int8_t jetFlavor = mcpjet.origin();

    if (detColl) {

      registry.fill(HIST("h_jetpT_particle_DetColl"), mcpjet.pt(), eventWeight);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_particle_DetColl_bjet"), mcpjet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_jetpT_particle_DetColl_cjet"), mcpjet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_DetColl_lfjet"), mcpjet.pt(), eventWeight);
      }
    } else {

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_particle_bjet"), mcpjet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_jetpT_particle_cjet"), mcpjet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_lfjet"), mcpjet.pt(), eventWeight);
      }
    }
  }

  void processDummy(FilteredCollision::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BJetTaggingML, processDummy, "Dummy process function turned on by default", true);

  void processDataJets(FilteredCollision::iterator const& collision, DataJets const& alljets, JetTrackswID const& allTracks, aod::DataSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    for (const auto& analysisJet : alljets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      processJetInfo(analysisJet, allTracks, allSVs);
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processDataJets, "jet information in Data", false);

  using MCDJetTableWeighted = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::MCDSecondaryVertex3ProngIndices, aod::ChargedMCDetectorLevelJetTags, aod::ChargedMCDetectorLevelJetEventWeights>>;
  using MCPJetTableWeighted = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef, aod::ChargedMCParticleLevelJetEventWeights>>;
  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs>>;

  Preslice<MCPJetTableWeighted> mcpJetsPerCollisionWeighted = aod::jet::mcCollisionId;

  void processMCJetsWeighted(FilteredCollisionMCD::iterator const& collision, MCDJetTableWeighted const& MCDjets, MCPJetTableWeighted const& MCPjets, JetTracksMCDwID const& allTracks, aod::JetParticles const& /*MCParticles*/, aod::MCDSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    auto const mcPJetsPerColl = MCPjets.sliceBy(mcpJetsPerCollisionWeighted, collision.mcCollisionId());

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float eventWeight = analysisJet.eventWeight();
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int jetFlavor = analysisJet.origin();

      processJetInfo(analysisJet, allTracks, allSVs, jetFlavor, eventWeight);

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTableWeighted>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_bjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_cjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lfjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        }
      }
    }

    // For filling histograms used for the jet matching efficiency
    for (const auto& mcpjet : mcPJetsPerColl) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float eventWeight = mcpjet.eventWeight();
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      fillMCPHistograms(mcpjet, true, eventWeight);
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processMCJetsWeighted, "jet information in MC with event weight", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::MCDSecondaryVertex3ProngIndices, aod::ChargedMCDetectorLevelJetTags>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef>>;

  Preslice<MCPJetTable> mcpJetsPerCollision = aod::jet::mcCollisionId;

  void processMCJets(FilteredCollisionMCD::iterator const& collision, MCDJetTable const& MCDjets, MCPJetTable const& MCPjets, JetTracksMCDwID const& allTracks, aod::JetParticles const& /*MCParticles*/, aod::MCDSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    auto const mcPJetsPerColl = MCPjets.sliceBy(mcpJetsPerCollision, collision.mcCollisionId());

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      int jetFlavor = analysisJet.origin();

      processJetInfo(analysisJet, allTracks, allSVs, jetFlavor);

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_bjet"), analysisJet.pt(), mcpjet.pt());
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_cjet"), analysisJet.pt(), mcpjet.pt());
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lfjet"), analysisJet.pt(), mcpjet.pt());
        }
      }
    }

    // For filling histograms used for the jet matching efficiency
    for (const auto& mcpjet : mcPJetsPerColl) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      fillMCPHistograms(mcpjet, true);
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processMCJets, "jet information in MC without event weight", false);

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionMCP = soa::Filtered<aod::JetMcCollisions>;

  void processMCTruthJetsWeighted(FilteredCollisionMCP::iterator const& /*collision*/, MCPJetTableWeighted const& MCPjets, aod::JetParticles const& /*MCParticles*/)
  {

    for (const auto& mcpjet : MCPjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float eventWeight = mcpjet.eventWeight();
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      fillMCPHistograms(mcpjet, false, eventWeight);
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processMCTruthJetsWeighted, "truth jet information with event weight", false);

  void processMCTruthJets(FilteredCollisionMCP::iterator const& /*collision*/, MCPJetTable const& MCPjets, aod::JetParticles const& /*MCParticles*/)
  {

    for (const auto& mcpjet : MCPjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      fillMCPHistograms(mcpjet);
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processMCTruthJets, "truth jet information without event weight", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<BJetTaggingML>(cfgc, TaskName{"bjet-tagging-ml"})}; // o2-linter: disable=name/o2-task,name/workflow-file
}
