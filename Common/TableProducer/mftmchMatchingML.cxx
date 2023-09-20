// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <math.h>
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <string>
#include <regex>
#include <TLorentzVector.h>
#include "Common/DataModel/MftmchMatchingML.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CCDB/CcdbApi.h"
#include "Tools/ML/model.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::ml;
using o2::globaltracking::MatchingFunc_t;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct mftmchMatchingML {
  Produces<aod::FwdTracksML> fwdtrackml;

  float etalow = -4;
  float etaup = -2.5;
  float pDCAcutrAtBsorberEndlow1 = 17.6;
  float pDCAcutrAtBsorberEndup1 = 26.5;
  float pDCAcutrAtBsorberEndlow2 = 26.5;
  float pDCAcutrAtBsorberEndup2 = 89.5;
  float pDCAcutdcaup1 = 594;
  float pDCAcutdcaup2 = 324;
  float chi2up = 1000000;
  float chi2MatchMCHMIDup = 1000000;

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (etaup < aod::fwdtrack::eta));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://ccdb-test.cern.ch:8080", "URL of the CCDB repository"};
  Configurable<std::string> cfgModelDir{"ccdb-path", "Users/m/mooya/models", "base path to the ONNX models"};
  Configurable<std::string> cfgModelName{"ccdb-file", "model_LHC22o.onnx", "name of ONNX model file"};
  Configurable<float> cfgThrScore{"threshold-score", 0.5, "Threshold value for matching score"};

  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "model-explorer"};
  Ort::SessionOptions session_options;
  std::shared_ptr<Ort::Experimental::Session> onnx_session = nullptr;
  OnnxModel model;

  template <typename F, typename M>
  std::vector<float> getVariables(F const& fwdtrack, M const& mfttrack)
  {

    static constexpr Double_t MatchingPlaneZ = -77.5;

    // propagate muontrack to matching position
    double muonchi2 = fwdtrack.chi2();
    SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
    std::vector<double> muonv1;
    SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
    o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
    muonpars1.propagateToZlinear(MatchingPlaneZ);

    // propagate mfttrack to matching position
    double mftchi2 = mfttrack.chi2();
    SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
    std::vector<double> mftv1;
    SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
    o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
    mftpars1.propagateToZlinear(MatchingPlaneZ);

    Float_t MFT_X = mftpars1.getX();
    Float_t MFT_Y = mftpars1.getY();
    Float_t MFT_Phi = mftpars1.getPhi();
    Float_t MFT_Tanl = mftpars1.getTanl();

    Float_t MCH_X = muonpars1.getX();
    Float_t MCH_Y = muonpars1.getY();
    Float_t MCH_Phi = muonpars1.getPhi();
    Float_t MCH_Tanl = muonpars1.getTanl();

    Float_t Ratio_X = MFT_X / MCH_X;
    Float_t Ratio_Y = MFT_Y / MCH_Y;
    Float_t Ratio_Phi = MFT_Phi / MCH_Phi;
    Float_t Ratio_Tanl = MFT_Tanl / MCH_Tanl;

    Float_t Delta_X = MFT_X - MCH_X;
    Float_t Delta_Y = MFT_Y - MCH_Y;
    Float_t Delta_Phi = MFT_Phi - MCH_Phi;
    Float_t Delta_Tanl = MFT_Tanl - MCH_Tanl;

    Float_t Delta_XY = sqrt(Delta_X * Delta_X + Delta_Y * Delta_Y);

    std::vector<float> input_tensor_values{
      MFT_X,
      MFT_Y,
      MFT_Phi,
      MFT_Tanl,
      MCH_X,
      MCH_Y,
      MCH_Phi,
      MCH_Tanl,
      Delta_XY,
      Delta_X,
      Delta_Y,
      Delta_Phi,
      Delta_Tanl,
      Ratio_X,
      Ratio_Y,
      Ratio_Phi,
      Ratio_Tanl,
    };
    return input_tensor_values;
  }

  template <typename F, typename M>
  double matchONNX(F const& fwdtrack, M const& mfttrack)
  {
    std::vector<std::string> input_names;
    std::vector<std::vector<int64_t>> input_shapes;
    std::vector<std::string> output_names;
    std::vector<std::vector<int64_t>> output_shapes;

    input_names = onnx_session->GetInputNames();
    input_shapes = onnx_session->GetInputShapes();
    output_names = onnx_session->GetOutputNames();
    output_shapes = onnx_session->GetOutputShapes();

    auto input_shape = input_shapes[0];
    input_shape[0] = 1;

    std::vector<float> input_tensor_values;
    input_tensor_values = getVariables(fwdtrack, mfttrack);

    if (input_tensor_values[8] < 3) {
      std::vector<Ort::Value> input_tensors;
      input_tensors.push_back(Ort::Experimental::Value::CreateTensor<float>(input_tensor_values.data(), input_tensor_values.size(), input_shape));

      std::vector<Ort::Value> output_tensors = onnx_session->Run(input_names, input_tensors, output_names);

      const float* output_value = output_tensors[0].GetTensorData<float>();

      auto score = output_value[0];

      return score;
    } else {
      auto score = 0;
      return score;
    }
  };

  void init(o2::framework::InitContext&)
  {
    o2::ccdb::CcdbApi ccdbApi;
    std::map<std::string, std::string> metadata;

    ccdbApi.init(cfgCCDBURL);
    // retrieving onnx file from ccdb
    std::string modelFile = cfgModelDir.value;
    bool retrieveSuccess = ccdbApi.retrieveBlob(modelFile, ".", metadata, 1642502592629, false, cfgModelName.value);

    // start session
    if (retrieveSuccess) {
      std::map<std::string, std::string> headers = ccdbApi.retrieveHeaders(modelFile, metadata, -1);
      LOG(info) << "Network file downloaded from: " << modelFile << " to: "
                << "."
                << "/" << cfgModelName.value;
      model.initModel(cfgModelName, false, 1, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
      onnx_session = model.getSession();
    } else {
      LOG(info) << "Failed to retrieve Network file";
    }
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        double result = matchONNX(fwdtrack, mfttrack);
        if (result > cfgThrScore) {
          double mftchi2 = mfttrack.chi2();
          SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
          std::vector<double> mftv1;
          SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
          o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
          mftpars1.propagateToZlinear(collision.posZ());

          float dcaX = (mftpars1.getX() - collision.posX());
          float dcaY = (mftpars1.getY() - collision.posY());
          double px = fwdtrack.p() * sin(M_PI / 2 - atan(mfttrack.tgl())) * cos(mfttrack.phi());
          double py = fwdtrack.p() * sin(M_PI / 2 - atan(mfttrack.tgl())) * sin(mfttrack.phi());
          double pz = fwdtrack.p() * cos(M_PI / 2 - atan(mfttrack.tgl()));
          fwdtrackml(fwdtrack.collisionId(), 0, mfttrack.x(), mfttrack.y(), mfttrack.z(), mfttrack.phi(), mfttrack.tgl(), fwdtrack.sign() / std::sqrt(std::pow(px, 2) + std::pow(py, 2)), fwdtrack.nClusters(), -1, -1, -1, -1, -1, result, mfttrack.globalIndex(), fwdtrack.globalIndex(), fwdtrack.mchBitMap(), fwdtrack.midBitMap(), fwdtrack.midBoards(), mfttrack.trackTime(), mfttrack.trackTimeRes(), mfttrack.eta(), std::sqrt(std::pow(px, 2) + std::pow(py, 2)), std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2)), dcaX, dcaY);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftmchMatchingML>(cfgc)};
}
