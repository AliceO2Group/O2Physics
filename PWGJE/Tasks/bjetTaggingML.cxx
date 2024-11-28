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
#include "Tools/ML/MlResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BJetTaggingML {

  struct bjetParams {
    float mJetpT = 0.0;
    float mJetEta = 0.0;
    float mJetPhi = 0.0;
    int mNTracks = -1;
    int mNSV = -1;
    float mJetMass = 0.0;
  };

  struct bjetTrackParams {
    double mTrackpT = 0.0;
    double mTrackEta = 0.0;
    double mDotProdTrackJet = 0.0;
    double mDotProdTrackJetOverJet = 0.0;
    double mDeltaRJetTrack = 0.0;
    double mSignedIP2D = 0.0;
    double mSignedIP2DSign = 0.0;
    double mSignedIP3D = 0.0;
    double mSignedIP3DSign = 0.0;
    double mMomFraction = 0.0;
    double mDeltaRTrackVertex = 0.0;
  };

  struct bjetSVParams {
    double mSVpT = 0.0;
    double mDeltaRSVJet = 0.0;
    double mSVMass = 0.0;
    double mSVfE = 0.0;
    double mIPXY = 0.0;
    double mCPA = 0.0;
    double mChi2PCA = 0.0;
    double mDispersion = 0.0;
    double mDecayLength2D = 0.0;
    double mDecayLength2DError = 0.0;
    double mDecayLength3D = 0.0;
    double mDecayLength3DError = 0.0;
  };

  HistogramRegistry registry;

  static constexpr double defaultCutsMl[1][2] = {{0.5, 0.5}};

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

  Configurable<bool> useQuarkDef{"useQuarkDef", true, "Flag whether to use quarks or hadrons for determining the jet flavor"};

  Configurable<float> svReductionFactor{"svReductionFactor", 1.0, "factor for how many SVs to keep"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{5., 1000.}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutSmaller, cuts_ml::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {defaultCutsMl[0], 1, 2, {"pT bin 0"}, {"score for default b-jet tagging", "uncer 1"}}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", 2, "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/h/hahassan"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ML_bjets/01-MVA/Models/LHC23d4_5_20_90Percent/model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::MlResponse<float> bMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  int eventSelection = -1;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    // Seed the random number generator using current time
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));

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

    bMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
    if (loadModelsFromCCDB) {
      ccdbApi.init(ccdbUrl);
      bMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
    } else {
      bMlResponse.setModelPathsLocal(onnxFileNames);
    }
    // bMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
    bMlResponse.init();
  }

  // FIXME filtering only works when you loop directly over the list, but if you loop over it as a constituent they will not be filtered
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt > trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt <= jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using FilteredCollision = soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs>>;
  using JetTrackswID = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using JetTracksMCDwID = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  using DataJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::DataSecondaryVertex3ProngIndices>>;

  std::vector<std::vector<float>> getInputsForML(bjetParams jetparams, std::vector<bjetTrackParams>& tracksParams, std::vector<bjetSVParams>& svsParams)
  {
    std::vector<float> jetInput = {jetparams.mJetpT, jetparams.mJetEta, jetparams.mJetPhi, static_cast<float>(jetparams.mNTracks), static_cast<float>(jetparams.mNSV), jetparams.mJetMass};
    std::vector<float> tracksInputFlat;
    std::vector<float> svsInputFlat;

    for (int iconstit = 0; iconstit < nJetConst; iconstit++) {

      tracksInputFlat.push_back(tracksParams[iconstit].mTrackpT);
      tracksInputFlat.push_back(tracksParams[iconstit].mTrackEta);
      tracksInputFlat.push_back(tracksParams[iconstit].mDotProdTrackJet);
      tracksInputFlat.push_back(tracksParams[iconstit].mDotProdTrackJetOverJet);
      tracksInputFlat.push_back(tracksParams[iconstit].mDeltaRJetTrack);
      tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP2D);
      tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP2DSign);
      tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP3D);
      tracksInputFlat.push_back(tracksParams[iconstit].mSignedIP3DSign);
      tracksInputFlat.push_back(tracksParams[iconstit].mMomFraction);
      tracksInputFlat.push_back(tracksParams[iconstit].mDeltaRTrackVertex);

      svsInputFlat.push_back(svsParams[iconstit].mSVpT);
      svsInputFlat.push_back(svsParams[iconstit].mDeltaRSVJet);
      svsInputFlat.push_back(svsParams[iconstit].mSVMass);
      svsInputFlat.push_back(svsParams[iconstit].mSVfE);
      svsInputFlat.push_back(svsParams[iconstit].mIPXY);
      svsInputFlat.push_back(svsParams[iconstit].mCPA);
      svsInputFlat.push_back(svsParams[iconstit].mChi2PCA);
      svsInputFlat.push_back(svsParams[iconstit].mDispersion);
      svsInputFlat.push_back(svsParams[iconstit].mDecayLength2D);
      svsInputFlat.push_back(svsParams[iconstit].mDecayLength2DError);
      svsInputFlat.push_back(svsParams[iconstit].mDecayLength3D);
      svsInputFlat.push_back(svsParams[iconstit].mDecayLength3DError);
    }

    std::vector<std::vector<float>> totalInput;
    totalInput.push_back(jetInput);
    totalInput.push_back(tracksInputFlat);
    totalInput.push_back(svsInputFlat);

    return totalInput;
  }

  // Looping over the SV info and writing them to a table
  template <typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetSVInfo(AnalysisJet const& myJet, AnyTracks const& /*allTracks*/, SecondaryVertices const& /*allSVs*/, std::vector<bjetSVParams>& svsParams, int jetFlavor = 0, double eventweight = 1.0)
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

      if (svsParams.size() < (svReductionFactor * myJet.template tracks_as<AnyTracks>().size())) {
        svsParams.emplace_back(bjetSVParams{candSV.pt(), deltaRJetSV, massSV, energySV / myJet.energy(), candSV.impactParameterXY(), candSV.cpa(), candSV.chi2PCA(), candSV.dispersion(), candSV.decayLengthXY(), candSV.errorDecayLengthXY(), candSV.decayLength(), candSV.errorDecayLength()});
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
  }

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetTrackInfo(AnyCollision const& /*collision*/, AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, SecondaryVertices const& /*allSVs*/, std::vector<bjetTrackParams>& tracksParams, int jetFlavor = 0, double eventweight = 1.0)
  {

    for (auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin) {
        continue;
      }

      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

      float RClosestSV = 10.;
      for (const auto& candSV : analysisJet.template secondaryVertices_as<SecondaryVertices>()) {
        double deltaRTrackSV = jetutilities::deltaR(constituent, candSV);
        if (deltaRTrackSV < RClosestSV) {
          RClosestSV = deltaRTrackSV;
        }
      }

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

      tracksParams.emplace_back(bjetTrackParams{constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaXYZ()) * sign, constituent.sigmadcaXYZ(), constituent.p() / analysisJet.p(), RClosestSV});
    }

    auto compare = [](bjetTrackParams& tr1, bjetTrackParams& tr2) {
      return (tr1.mSignedIP2D / tr1.mSignedIP2DSign) > (tr2.mSignedIP2D / tr2.mSignedIP2DSign);
    };

    // Sort the tracks based on their IP significance in descending order
    std::sort(tracksParams.begin(), tracksParams.end(), compare);
  }

  void processDummy(FilteredCollision::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BJetTaggingML, processDummy, "Dummy process function turned on by default", true);

  void processDataJets(FilteredCollision::iterator const& collision, DataJets const& alljets, JetTrackswID const& allTracks, aod::DataSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    for (const auto& analysisJet : alljets) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      std::vector<bjetTrackParams> tracksParams;
      std::vector<bjetSVParams> SVsParams;

      analyzeJetSVInfo(analysisJet, allTracks, allSVs, SVsParams);
      analyzeJetTrackInfo(collision, analysisJet, allTracks, allSVs, tracksParams);

      int nSVs = analysisJet.template secondaryVertices_as<aod::DataSecondaryVertex3Prongs>().size();

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), tracksParams.size());
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), nSVs < 250 ? nSVs : 249);

      bjetParams jetparam = {analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), static_cast<int>(tracksParams.size()), static_cast<int>(nSVs), analysisJet.mass()};
      tracksParams.resize(nJetConst); // resize to the number of inputs of the ML
      SVsParams.resize(nJetConst);    // resize to the number of inputs of the ML

      auto inputML = getInputsForML(jetparam, tracksParams, SVsParams);

      std::vector<float> output;
      // bool isSelectedMl = bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);

      registry.fill(HIST("h2_score_jetpT"), analysisJet.pt(), output[0]);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass());
    }
  }
  PROCESS_SWITCH(BJetTaggingML, processDataJets, "jet information in Data", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::MCDSecondaryVertex3ProngIndices, aod::ChargedMCDetectorLevelJetEventWeights>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>>;
  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>;

  Preslice<aod::JMcParticles> McParticlesPerCollision = aod::jmcparticle::mcCollisionId;
  Preslice<MCPJetTable> McPJetsPerCollision = aod::jet::mcCollisionId;

  void processMCJets(FilteredCollisionMCD::iterator const& collision, MCDJetTable const& MCDjets, MCPJetTable const& MCPjets, JetTracksMCDwID const& allTracks, aod::JetParticles const& MCParticles, aod::MCDSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    auto const mcParticlesPerColl = MCParticles.sliceBy(McParticlesPerCollision, collision.mcCollisionId());
    auto const mcPJetsPerColl = MCPjets.sliceBy(McPJetsPerCollision, collision.mcCollisionId());

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
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

      std::vector<bjetTrackParams> tracksParams;
      std::vector<bjetSVParams> SVsParams;

      int jetFlavor = 0;

      for (auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        if (useQuarkDef) {
          jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, mcParticlesPerColl);
        } else {
          jetFlavor = jettaggingutilities::getJetFlavorHadron(mcpjet, mcParticlesPerColl);
        }
      }

      analyzeJetSVInfo(analysisJet, allTracks, allSVs, SVsParams, jetFlavor, eventWeight);
      analyzeJetTrackInfo(collision, analysisJet, allTracks, allSVs, tracksParams, jetFlavor, eventWeight);

      int nSVs = analysisJet.template secondaryVertices_as<aod::MCDSecondaryVertex3Prongs>().size();

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), tracksParams.size());
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), nSVs < 250 ? nSVs : 249);

      bjetParams jetparam = {analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), static_cast<int>(tracksParams.size()), static_cast<int>(nSVs), analysisJet.mass()};
      tracksParams.resize(nJetConst); // resize to the number of inputs of the ML
      SVsParams.resize(nJetConst);    // resize to the number of inputs of the ML

      auto inputML = getInputsForML(jetparam, tracksParams, SVsParams);

      std::vector<float> output;
      // bool isSelectedMl = bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);

      registry.fill(HIST("h2_score_jetpT"), analysisJet.pt(), output[0], eventWeight);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass(), eventWeight);

      if (jetFlavor == 2) {
        registry.fill(HIST("h2_score_jetpT_bjet"), analysisJet.pt(), output[0], eventWeight);
        registry.fill(HIST("h2_jetMass_jetpT_bjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_bjet"), analysisJet.pt(), eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h2_score_jetpT_cjet"), analysisJet.pt(), output[0], eventWeight);
        registry.fill(HIST("h2_jetMass_jetpT_cjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_cjet"), analysisJet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h2_score_jetpT_lfjet"), analysisJet.pt(), output[0], eventWeight);
        registry.fill(HIST("h2_jetMass_jetpT_lfjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_lfjet"), analysisJet.pt(), eventWeight);
      }

      for (auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
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
      for (auto jetR : jetRadiiValues) {
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

      int8_t jetFlavor = 0;

      if (useQuarkDef) {
        jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, mcParticlesPerColl);
      } else {
        jetFlavor = jettaggingutilities::getJetFlavorHadron(mcpjet, mcParticlesPerColl);
      }

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
  using FilteredCollisionMCP = soa::Filtered<aod::JMcCollisions>;

  void processMCTruthJets(FilteredCollisionMCP::iterator const& /*collision*/, MCPJetTable const& MCPjets, aod::JetParticles const& MCParticles)
  {

    for (const auto& mcpjet : MCPjets) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
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

      int8_t jetFlavor = 0;

      if (useQuarkDef) {
        jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, MCParticles);
      } else {
        jetFlavor = jettaggingutilities::getJetFlavorHadron(mcpjet, MCParticles);
      }

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
  return WorkflowSpec{adaptAnalysisTask<BJetTaggingML>(cfgc, TaskName{"bjet-tagging-ml"})};
}
