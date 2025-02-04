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

#include <string>
#include <algorithm>
#include <vector>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BJetTaggingML {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.5, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};

  // track level configurables
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> nJetConst{"nJetConst", 10, "maximum number of jet consistuents to be used for ML evaluation"};

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

    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{40, -20.0, 20.0}}});

    registry.add("h2_score_jetpT", "ML scores for inclusive jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, {120, -0.1, 1.1}}});

    registry.add("h2_nTracks_jetpT", "Number of tracks;#it{p}_{T,jet} (GeV/#it{c});nTracks", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 100.0}}});
    registry.add("h2_nSV_jetpT", "Number of secondary vertices;#it{p}_{T,jet} (GeV/#it{c});nSVs", {HistType::kTH2F, {{200, 0., 200.}, {250, 0, 250.0}}});

    registry.add("h2_SIPs2D_jetpT", "2D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_SIPs3D_jetpT", "3D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_LxyS_jetpT", "Decay length in XY;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
    registry.add("h2_Dispersion_jetpT", "SV dispersion;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
    registry.add("h2_jetMass_jetpT", "Jet mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
    registry.add("h2_SVMass_jetpT", "Secondary vertex mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10}}});

    if (doDataDriven) {
      registry.add("hSparse_Incljets", "Inclusive jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, {120, -0.1, 1.1}, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}});
      if (doprocessMCJets) {
        registry.add("hSparse_bjets", "Tagged b-jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, {120, -0.1, 1.1}, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}});
        registry.add("hSparse_cjets", "Tagged c-jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, {120, -0.1, 1.1}, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}});
        registry.add("hSparse_lfjets", "Tagged lf-jets Info;#it{p}_{T,jet} (GeV/#it{c});Score;#it{m}_{jet} (GeV/#it{c}^{2});-log(JP);#it{m}_{SV} (GeV/#it{c}^{2});SVfE;", {HistType::kTHnSparseF, {{200, 0., 200.}, {120, -0.1, 1.1}, {50, 0, 50}, {375, 0, 30}, {50, 0, 10}, {50, 0, 1}}});
      }
    }

    if (doprocessMCJets) {

      registry.add("h2_score_jetpT_bjet", "ML scores for b-jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, {120, -0.1, 1.1}}});
      registry.add("h2_SIPs2D_jetpT_bjet", "2D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_bjet", "3D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_bjet", "Decay length in XY b-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_bjet", "SV dispersion b-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_bjet", "Jet mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_bjet", "Secondary vertex mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_score_jetpT_cjet", "ML scores for c-jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, {120, -0.1, 1.1}}});
      registry.add("h2_SIPs2D_jetpT_cjet", "2D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_cjet", "3D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_cjet", "Decay length in XY c-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_cjet", "SV dispersion c-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_cjet", "Jet mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_cjet", "Secondary vertex mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_score_jetpT_lfjet", "ML scores for lf-jets;#it{p}_{T,jet} (GeV/#it{c});Score", {HistType::kTH2F, {{200, 0., 200.}, {120, -0.1, 1.1}}});
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

      registry.add("h2_Response_DetjetpT_PartjetpT_bjet", "Response matrix b-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h2_Response_DetjetpT_PartjetpT_cjet", "Response matrix c-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h2_Response_DetjetpT_PartjetpT_lfjet", "Response matrix lf-jet;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
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

      if (doprocessMCJets) {
        if (jetFlavor == 2) {
          registry.fill(HIST("h2_LxyS_jetpT_bjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_bjet"), myJet.pt(), candSV.dispersion(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_bjet"), myJet.pt(), massSV, eventweight);
        } else if (jetFlavor == 1) {
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

      if (constituent.pt() < trackPtMin) {
        continue;
      }

      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

      float rClosestSV = -1.0;

      registry.fill(HIST("h2_SIPs2D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
      registry.fill(HIST("h2_SIPs3D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);

      if (doprocessMCJets) {
        if (jetFlavor == 2) {
          registry.fill(HIST("h2_SIPs2D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else if (jetFlavor == 1) {
          registry.fill(HIST("h2_SIPs2D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else {
          registry.fill(HIST("h2_SIPs2D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        }
      }

      tracksParams.emplace_back(jettaggingutilities::BJetTrackParams{constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaXYZ()) * sign, constituent.sigmadcaXYZ(), constituent.p() / analysisJet.p(), rClosestSV});
    }

    auto compare = [](jettaggingutilities::BJetTrackParams& tr1, jettaggingutilities::BJetTrackParams& tr2) {
      return (tr1.mSignedIP2D / tr1.mSignedIP2DSign) > (tr2.mSignedIP2D / tr2.mSignedIP2DSign);
    };

    // Sort the tracks based on their IP significance in descending order
    std::sort(tracksParams.begin(), tracksParams.end(), compare);
    tracksParams.resize(nJetConst);
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

      std::vector<jettaggingutilities::BJetTrackParams> tracksParams;
      std::vector<jettaggingutilities::BJetSVParams> svsParams;

      analyzeJetSVInfo(analysisJet, allSVs, svsParams);
      analyzeJetTrackInfo(analysisJet, allTracks, tracksParams);

      int nSVs = analysisJet.template secondaryVertices_as<aod::DataSecondaryVertex3Prongs>().size();
      int nTracks = analysisJet.template tracks_as<JetTrackswID>().size();

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), nTracks);
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), nSVs < 250 ? nSVs : 249);

      registry.fill(HIST("h2_score_jetpT"), analysisJet.pt(), analysisJet.scoreML());

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass());

      if (doDataDriven) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), analysisJet.mass(), -1 * std::log(analysisJet.jetProb()), svsParams[0].mSVMass, svsParams[0].mSVfE);
      }
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processDataJets, "jet information in Data", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::MCDSecondaryVertex3ProngIndices, aod::ChargedMCDetectorLevelJetTags, aod::ChargedMCDetectorLevelJetEventWeights>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef, aod::ChargedMCParticleLevelJetEventWeights>>;
  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs>>;

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

      float eventWeight = useEventWeight ? analysisJet.eventWeight() : 1.0;
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      std::vector<jettaggingutilities::BJetTrackParams> tracksParams;
      std::vector<jettaggingutilities::BJetSVParams> svsParams;

      int jetFlavor = analysisJet.origin();

      analyzeJetSVInfo(analysisJet, allSVs, svsParams, jetFlavor, eventWeight);
      analyzeJetTrackInfo(analysisJet, allTracks, tracksParams, jetFlavor, eventWeight);

      int nSVs = analysisJet.template secondaryVertices_as<aod::MCDSecondaryVertex3Prongs>().size();
      int nTracks = analysisJet.template tracks_as<JetTracksMCDwID>().size();

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), nTracks);
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), nSVs < 250 ? nSVs : 249);

      registry.fill(HIST("h2_score_jetpT"), analysisJet.pt(), analysisJet.scoreML(), eventWeight);
      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass(), eventWeight);

      if (doDataDriven) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), analysisJet.mass(), -1 * std::log(analysisJet.jetProb()), svsParams[0].mSVMass, svsParams[0].mSVfE);
        if (jetFlavor == 2) {
          registry.fill(HIST("hSparse_bjets"), analysisJet.pt(), analysisJet.scoreML(), analysisJet.mass(), -1 * std::log(analysisJet.jetProb()), svsParams[0].mSVMass, svsParams[0].mSVfE);
        } else if (jetFlavor == 1) {
          registry.fill(HIST("hSparse_cjets"), analysisJet.pt(), analysisJet.scoreML(), analysisJet.mass(), -1 * std::log(analysisJet.jetProb()), svsParams[0].mSVMass, svsParams[0].mSVfE);
        } else {
          registry.fill(HIST("hSparse_lfjets"), analysisJet.pt(), analysisJet.scoreML(), analysisJet.mass(), -1 * std::log(analysisJet.jetProb()), svsParams[0].mSVMass, svsParams[0].mSVfE);
        }
      }

      if (jetFlavor == 2) {
        registry.fill(HIST("h2_score_jetpT_bjet"), analysisJet.pt(), analysisJet.scoreML(), eventWeight);
        registry.fill(HIST("h2_jetMass_jetpT_bjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_bjet"), analysisJet.pt(), eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h2_score_jetpT_cjet"), analysisJet.pt(), analysisJet.scoreML(), eventWeight);
        registry.fill(HIST("h2_jetMass_jetpT_cjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_cjet"), analysisJet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h2_score_jetpT_lfjet"), analysisJet.pt(), analysisJet.scoreML(), eventWeight);
        registry.fill(HIST("h2_jetMass_jetpT_lfjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_lfjet"), analysisJet.pt(), eventWeight);
      }

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        if (jetFlavor == 2) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_bjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else if (jetFlavor == 1) {
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

      float eventWeight = useEventWeight ? mcpjet.eventWeight() : 1.0;
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      int8_t jetFlavor = mcpjet.origin();

      registry.fill(HIST("h_jetpT_particle_DetColl"), mcpjet.pt(), eventWeight);

      if (jetFlavor == 2) {
        registry.fill(HIST("h_jetpT_particle_DetColl_bjet"), mcpjet.pt(), eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h_jetpT_particle_DetColl_cjet"), mcpjet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_DetColl_lfjet"), mcpjet.pt(), eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processMCJets, "jet information in MC", false);

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionMCP = soa::Filtered<aod::JetMcCollisions>;

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

      float eventWeight = useEventWeight ? mcpjet.eventWeight() : 1.0;
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      int8_t jetFlavor = mcpjet.origin();

      if (jetFlavor == 2) {
        registry.fill(HIST("h_jetpT_particle_bjet"), mcpjet.pt(), eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h_jetpT_particle_cjet"), mcpjet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_lfjet"), mcpjet.pt(), eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processMCTruthJets, "truth jet information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<BJetTaggingML>(cfgc, TaskName{"bjet-tagging-ml"})}; // o2-linter: disable=name/o2-task,name/workflow-file
}
