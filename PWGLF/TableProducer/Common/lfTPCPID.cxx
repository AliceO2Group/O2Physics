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

///
/// \file   lfTPCPID.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-11-20
/// \brief  Task to produce the PID information for the TPC for the purpose of the Light flavor PWG
///

// ROOT includes
#include "TFile.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraph.h"
#include "TList.h"

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/CcdbApi.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "TableHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

static constexpr int nSpecies = 9;
static constexpr int nParameters = 11;
static const std::vector<std::string> particleNames{"El", "Mu", "Pi", "Ka", "Pr", "De", "Tr", "He", "Al"};
static const std::vector<std::string> parameterNames{"Use default tiny",
                                                     "Use default full",
                                                     "Set parameters",
                                                     "bb1", "bb2", "bb3", "bb4", "bb5",
                                                     "MIP value", "Charge exponent", "Resolution"};
static constexpr float defaultParameters[nSpecies][nParameters]{{2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f},
                                                                {2.f, 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
static constexpr int nOptions = 4;
static const std::vector<std::string> optionNames{"Bethe Bloch path",      // If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://
                                                  "Post calibration path", // If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://
                                                  "Simple Bethe Bloch",    // (true/false)
                                                  "RequirePostCalib"};     // (true/false)
const std::string defaultOptions[nSpecies][nOptions]{{"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"},
                                                     {"", "", "false", "false"}};

// Structure to hold the parameters
struct bbParams {
  const std::string name;
  bbParams(const std::string& n) : name(n) {}
  // Parameters for the Bethe-Bloch parametrization
  float bb1 = 0.03209809958934784f;    // Aleph Bethe Bloch parameter 1
  float bb2 = 19.9768009185791f;       // Aleph Bethe Bloch parameter 2
  float bb3 = 2.5266601063857674e-16f; // Aleph Bethe Bloch parameter 3
  float bb4 = 2.7212300300598145f;     // Aleph Bethe Bloch parameter 4
  float bb5 = 6.080920219421387f;      // Aleph Bethe Bloch parameter 5
  float mip = 50.f;                    // MIP value
  float exp = 2.299999952316284f;      // Exponent of the charge factor
  float res = 0.1f;                    // Resolution

  TGraph* postCorrection = nullptr;
  TGraph* postCorrectionSigma = nullptr;

  TF1* postCorrectionFun = nullptr;
  TF1* postCorrectionFunSigma = nullptr;

  // Utility parameters for the usage
  bool betheBlochSet = true;     // Flag to check if the Bethe-Bloch parameters have been set. By default is true as the default values are set. Used to check if the parameters have been set after a CCDB update
  bool requirePostCalib = true;  // Flag to force the post calib. to be required, if not found, it will trigger a fatal error
  bool takeFromCcdb = false;     // Flag to get the parameters from the CCDB
  bool takePostFromCcdb = false; // Flag to get the post calib parameters from the CCDB
  std::string ccdbPath = "";     // Path to the CCDB object
  std::string ccdbPathPost = ""; // Path to the CCDB object for the post calib
  int lastRunNumber = 0;         // Last processed run
  bool isSimple = false;         // Flag to use only the Bethe-Bloch parameters without the charge exponent and the MIP value

  ///
  /// Set the values of the BetheBloch from an array of parameters
  bool setValues(std::vector<float> v)
  {
    if (v.size() != 8) {
      LOG(fatal) << "bbParams `" << name << "` :: The vector of Bethe-Bloch parameters has the wrong size " << v.size() << " while expecting 8";
      return false;
    }
    LOG(info) << "bbParams `" << name << "` :: Before: set of parameters -> bb1: " << bb1 << ", bb2: " << bb2 << ", bb3: " << bb3 << ", bb4: " << bb4 << ", bb5: " << bb5 << ", mip: " << mip << ", exp: " << exp << ", resolution " << res;
    bb1 = v[0];
    bb2 = v[1];
    bb3 = v[2];
    bb4 = v[3];
    bb5 = v[4];
    mip = v[5];
    exp = v[6];
    res = v[7];
    LOG(info) << "bbParams `" << name << "` :: After: set of parameters -> bb1: " << bb1 << ", bb2: " << bb2 << ", bb3: " << bb3 << ", bb4: " << bb4 << ", bb5: " << bb5 << ", mip: " << mip << ", exp: " << exp << ", resolution " << res;
    betheBlochSet = true; // Bethe bloch parameters are set!
    return true;
  }

  ///
  /// Set the values of the BetheBloch from TH1F
  bool setValues(TH1F* h)
  {
    if (!h) {
      LOG(fatal) << "bbParams `" << name << "` :: The input histogram of Bethe-Bloch parameters is not valid";
      return false;
    }
    if (isSimple) {
      return setValuesSimple(h);
    }
    const int n = h->GetNbinsX();
    TAxis* axis = h->GetXaxis();
    std::vector<float> v{static_cast<float>(h->GetBinContent(axis->FindBin("bb1"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb2"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb3"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb4"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb5"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("MIP value"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("Charge exponent"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("Resolution")))};
    if (h->GetNbinsX() != n) {
      LOG(fatal) << "The input histogram of Bethe-Bloch parameters has the wrong size " << n << " while expecting " << h->GetNbinsX();
      return false;
    }
    LOG(info) << "bbParams `" << name << "` :: Setting custom Bethe-Bloch parameters from histogram " << h->GetName();
    return setValues(v);
  }

  ///
  /// Sets only the Bethe-Bloch parameters without the charge exponent and the MIP value
  bool setValuesSimple(TH1F* h)
  {
    const int n = h->GetNbinsX();
    TAxis* axis = h->GetXaxis();
    std::vector<float> v{static_cast<float>(h->GetBinContent(axis->FindBin("bb1"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb2"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb3"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb4"))),
                         static_cast<float>(h->GetBinContent(axis->FindBin("bb5"))),
                         1.f,
                         0.f,
                         static_cast<float>(h->GetBinContent(axis->FindBin("Resolution")))};
    if (h->GetNbinsX() != n) {
      LOG(fatal) << "The input histogram of Bethe-Bloch parameters has the wrong size " << n << " while expecting " << h->GetNbinsX();
      return false;
    }
    // Check that it is indeed simple
    axis->FindBin("MIP value");
    axis->FindBin("Charge exponent");
    if (h->GetNbinsX() == n) {
      LOG(fatal) << "The input histogram of Bethe-Bloch parameters is compatible with the full parametrization";
      return false;
    }
    LOG(info) << "bbParams `" << name << "` :: Setting simple custom Bethe-Bloch parameters from histogram " << h->GetName();
    return setValues(v);
  }

  ///
  /// @brief Set the Bethe Bloch calibration from a TList
  /// @param l List containing the post calibration objects as TF1 or TGraph
  void setValues(TList* l)
  {
    const char* bbname = Form("BBParameter_%s", name.c_str());
    TH1F* h = static_cast<TH1F*>(l->FindObject(bbname));
    if (!h) {
      l->ls();
      LOG(fatal) << "bbParams `" << name << "` :: did not find BetheBloch parametrization " << bbname << " in the input list";
    }
    setValues(h);
  }

  ///
  /// @brief Set the post calibration from a TList
  /// @param l List containing the post calibration objects as TF1 or TGraph
  void setPostCorrection(TList* l)
  {
    if (!l) {
      LOG(fatal) << "bbParams `" << name << "` :: Did not find the post calib list";
      return;
    }

    // Reset post calibration
    postCorrectionFun = nullptr;
    postCorrection = nullptr;
    postCorrectionFunSigma = nullptr;
    postCorrectionSigma = nullptr;

    TString objname = Form("dEdx_MinusBB_%s", name.c_str());
    TObject* obj = l->FindObject(objname);
    if (!obj) {
      if (requirePostCalib) {
        l->ls();
        LOG(fatal) << "Did not find the post calib object " << objname << ", cannot assign post calibration";
      }
    } else {
      TString cn = obj->ClassName();
      if (cn.Contains("TF1")) {
        postCorrectionFun = static_cast<TF1*>(l->FindObject(obj->GetName()));
      } else if (cn.Contains("TGraph")) {
        postCorrection = static_cast<TGraph*>(l->FindObject(obj->GetName()));
      } else {
        LOG(fatal) << "Cannot hanlde class " << cn << " for post calib object " << objname;
      }
    }

    objname = Form("sigmaFitOverSigmaParam_%s", name.c_str());
    obj = l->FindObject(objname);
    if (!obj) {
      if (requirePostCalib) {
        l->ls();
        LOG(fatal) << "Did not find the post calib sigma object " << objname << ", cannot assign post calibration";
      }
    } else {
      TString cn = obj->ClassName();
      if (cn.Contains("TF1")) {
        postCorrectionFunSigma = static_cast<TF1*>(l->FindObject(obj->GetName()));
      } else if (cn.Contains("TGraph")) {
        postCorrectionSigma = static_cast<TGraph*>(l->FindObject(obj->GetName()));
      } else {
        LOG(fatal) << "Cannot hanlde class " << cn << " for post calib sigma object " << objname;
      }
    }
    LOG(info) << "bbParams `" << name << "` :: Setting post calibration";
    if (l->FindObject("Tag")) {
      LOG(info) << "          Tag: " << l->FindObject("Tag")->GetTitle();
    }
    if (postCorrectionFun) {
      LOG(info) << "          postCorrection: " << postCorrectionFun->ClassName() << " " << postCorrectionFun->GetName();
    } else if (postCorrection) {
      LOG(info) << "          postCorrection: " << postCorrection->ClassName() << " " << postCorrection->GetName();
    } else {
      LOG(info) << "          postCorrection: Not assigned";
    }
    if (postCorrectionFunSigma) {
      LOG(info) << "          postCorrectionSigma: " << postCorrectionFunSigma->ClassName() << " " << postCorrectionFunSigma->GetName();
    } else if (postCorrectionSigma) {
      LOG(info) << "          postCorrectionSigma: " << postCorrectionSigma->ClassName() << " " << postCorrectionSigma->GetName();
    } else {
      LOG(info) << "          postCorrectionSigma: Not assigned";
    }
  }

  ///
  /// Set values from a configuration. In this case also the post calibration is checked
  void init(const char* particle,
            const Configurable<LabeledArray<float>>& confParams,
            const Configurable<LabeledArray<std::string>>& cfg,
            o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj)
  {
    LOG(info) << "bbParams `" << name << "` :: inizializing parameters from json with particle " << particle;
    // First check the json configuration
    if (confParams->get(particle, "Set parameters") < 1.5f) { // Keep the default hardcoded values
      LOG(info) << "bbParams `" << name << "` :: Using default parameters for " << particle << " as entry `Set parameters` " << confParams->get(particle, "Set parameters") << " < 1.5";
    } else {
      LOG(info) << "bbParams `" << name << "` :: Setting custom Bethe-Bloch parameters from JSON for mass hypothesis " << particle;
      std::vector<float> v{confParams->get(particle, "bb1"),
                           confParams->get(particle, "bb2"),
                           confParams->get(particle, "bb3"),
                           confParams->get(particle, "bb4"),
                           confParams->get(particle, "bb5"),
                           confParams->get(particle, "MIP value"),
                           confParams->get(particle, "Charge exponent"),
                           confParams->get(particle, "Resolution")};
      setValues(v);
    }

    // Check the TFile/CCDB configuration
    const std::string bb = cfg->get(particle, "Bethe Bloch path");
    const std::string post = cfg->get(particle, "Post calibration path");
    const std::string simple = cfg->get(particle, "Simple Bethe Bloch");
    requirePostCalib = (cfg->get(particle, "RequirePostCalib") == "true");
    LOG(info) << "bbParams `" << name << "` :: initializing from configurable '" << cfg.name << "' bb = " << bb << " post = " << post << " simple = " << simple;
    if (simple == "true" || simple == "yes") { // Check if the parametrization to be used is simple
      LOG(info) << "bbParams `" << name << "` :: using simple BetheBloch parametrization";
      isSimple = true;
    } else {
      LOG(info) << "bbParams `" << name << "` :: using full BetheBloch parametrization";
    }

    // First we check the post calib. In case the post calib is enabled, the BB parameters should be also coming from same source
    if (post.size() > 1) {
      if (bb.size() > 1) {
        LOG(fatal) << "Cannot define post calibration path and BB parameters path at the same time as the BB parameters will be taken from the post. calib. source. Pick one!";
      }
      if (requirePostCalib) {
        LOG(info) << "bbParams `" << name << "` :: fatal errors if post calibrations will not be found!";
      }
      LOG(info) << "bbParams `" << name << "` :: Loading calib. and post calib. from configurable '" << cfg.name << "' with value '" << post << "'";
      std::string s = post;
      if (s.rfind("ccdb://", 0) == 0) { // Getting post calib. from CCDB
        s.replace(0, 7, "");
        ccdbPathPost = s;
        if (ccdbObj->getTimestamp() == 0) { // If the timestamp is 0 we expect to get the timestamp from the run number
          takePostFromCcdb = true;
          LOG(info) << "bbParams `" << name << "` :: For post calib asked to query the parameters from the CCDB and got timestamp " << ccdbObj->getTimestamp() << " -> will take the object corresponding to the run number";
        } else {
          LOG(info) << "bbParams `" << name << "` :: For post calib fetching parameters from CCDB (only once) using timestamp " << ccdbObj->getTimestamp() << " and path " << s;
          TList* l = ccdbObj->get<TList>(s);
          setPostCorrection(l);
          setValues(l);
        }
      } else { // Getting post calib. from file
        TFile* f = TFile::Open(post.c_str(), "READ");
        if (!f) {
          LOG(fatal) << "The input file " << post << " is not valid";
        }
        if (!f->IsOpen()) {
          LOG(fatal) << "The input file " << f->GetName() << " is not open";
        }
        TList* l = static_cast<TList*>(f->Get("ccdb_object"));
        if (!l) {
          f->ls();
          LOG(fatal) << "The input file " << post << " does not contain the TList ccdb_object";
        }
        setPostCorrection(l);
        setValues(l);
      }
    }

    // Check the BetheBloch parameters
    if (bb.size() > 1) {
      LOG(info) << "bbParams `" << name << "` :: Loading parameters from configurable '" << cfg.name << "' with value '" << bb << "'";
      std::string s = bb;
      if (s.rfind("ccdb://", 0) == 0) { // Check if the path is a CCDB path
        s.replace(0, 7, "");
        ccdbPath = s;
        if (ccdbObj->getTimestamp() == 0) { // If the timestamp is 0 we expect to get the timestamp from the run number
          takeFromCcdb = true;
          LOG(info) << "bbParams `" << name << "` :: Asked to query the parameters from the CCDB and got timestamp " << ccdbObj->getTimestamp() << " -> will take the object corresponding to the run number";
        } else {
          LOG(info) << "bbParams `" << name << "` :: Fetching parameters from CCDB (only once) using timestamp " << ccdbObj->getTimestamp() << " and path " << s;
          setValues(ccdbObj->get<TH1F>(s));
        }
      } else { // Getting BetheBloch parameters from file
        TFile* f = TFile::Open(bb.c_str(), "READ");
        if (!f) {
          LOG(fatal) << "The input file " << post << " is not valid";
        }
        if (!f->IsOpen()) {
          LOG(fatal) << "The input file " << f->GetName() << " is not open";
        }
        TH1F* h = nullptr;
        f->GetObject("hpar", h);
        if (!h) {
          // Reattempt with ccdb name
          f->GetObject("ccdb_object", h);
          if (!h) {
            f->ls();
            LOG(fatal) << "The input file does not contain the histograms hpar or ccdb_object";
          }
        }
        LOG(info) << "bbParams `" << name << "` :: Setting parameters from TH1F " << h->GetName() << " in file " << f->GetName();
        setValues(h);
      }
    }
  }

  /// @brief Function to update the values of the parameters from the CCDB
  /// @param ccdbObj Managare of the CCDB
  /// @param timestamp Timestamp to ask the new parameters
  /// @return false if not succesfull
  void updateValues(aod::BCsWithTimestamps::iterator const& bunchCrossing, o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj)
  {
    if (!takePostFromCcdb && !takeFromCcdb) {
      LOG(debug) << "bbParams `" << name << "` :: Requested to not update parameters, ccdb timestamp " << ccdbObj->getTimestamp() << " != 0";
      return;
    }
    // Check that the last updated number is different
    if (lastRunNumber == bunchCrossing.runNumber()) {
      LOG(debug) << "bbParams `" << name << "` :: Not updating parameters of " << name << " from run number " << lastRunNumber << " as they are already up to date";
      return;
    }
    betheBlochSet = false; // Reset the set betheBloch flag as it could also be not found

    // First check the post calib
    if (takePostFromCcdb) {
      LOG(info) << "bbParams `" << name << "` :: Updating post calib from run number " << lastRunNumber << " to " << bunchCrossing.runNumber() << ". Taking them from CCDB path '" << ccdbPathPost << "' with timestamp " << bunchCrossing.timestamp();
      TList* postL = ccdbObj->getForTimeStamp<TList>(ccdbPathPost, bunchCrossing.timestamp());
      if (!postL) { // If the CCDB does not throw a fatal error, we use the centrally provided values
        LOG(info) << "bbParams `" << name << "` :: Post calibration input not found on CCDB, using default values";
        return;
      }
      setPostCorrection(postL);
      setValues(postL);
      lastRunNumber = bunchCrossing.runNumber();
      return;
    }
    // Secondly we check the BB parameters
    LOG(info) << "bbParams `" << name << "` :: Updating parameters of " << name << " from run number " << lastRunNumber << " to " << bunchCrossing.runNumber() << ". Taking them from CCDB path '" << ccdbPath << "' with timestamp " << bunchCrossing.timestamp();
    lastRunNumber = bunchCrossing.runNumber();
    betheBlochSet = false; // Reset the set betheBloch flag as it could also be not found
    TH1F* h = ccdbObj->getForTimeStamp<TH1F>(ccdbPath, bunchCrossing.timestamp());
    if (!h) { // If the CCDB does not throw a fatal error, we use the centrally provided values
      LOG(info) << "bbParams `" << name << "` :: Bethe Bloch calibration input not found on CCDB, using default values";
      return;
    }
    setValues(h);
  }
};

/// Task to produce the response table
struct lfTpcPid {
  using Trks = soa::Join<aod::TracksIU, aod::TracksExtra>;
  using Colls = aod::Collisions;

  // Tables to produce
  Produces<o2::aod::pidTPCLfEl> tablePIDEl;
  Produces<o2::aod::pidTPCLfMu> tablePIDMu;
  Produces<o2::aod::pidTPCLfPi> tablePIDPi;
  Produces<o2::aod::pidTPCLfKa> tablePIDKa;
  Produces<o2::aod::pidTPCLfPr> tablePIDPr;
  Produces<o2::aod::pidTPCLfDe> tablePIDDe;
  Produces<o2::aod::pidTPCLfTr> tablePIDTr;
  Produces<o2::aod::pidTPCLfHe> tablePIDHe;
  Produces<o2::aod::pidTPCLfAl> tablePIDAl;

  Produces<o2::aod::pidTPCLfFullEl> tablePIDFullEl;
  Produces<o2::aod::pidTPCLfFullMu> tablePIDFullMu;
  Produces<o2::aod::pidTPCLfFullPi> tablePIDFullPi;
  Produces<o2::aod::pidTPCLfFullKa> tablePIDFullKa;
  Produces<o2::aod::pidTPCLfFullPr> tablePIDFullPr;
  Produces<o2::aod::pidTPCLfFullDe> tablePIDFullDe;
  Produces<o2::aod::pidTPCLfFullTr> tablePIDFullTr;
  Produces<o2::aod::pidTPCLfFullHe> tablePIDFullHe;
  Produces<o2::aod::pidTPCLfFullAl> tablePIDFullAl;

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<bool> skipTPCOnly{"skipTPCOnly", true, "Flag to skip TPC only tracks (faster but affects the analyses that use TPC only tracks)"};
  Configurable<bool> fatalOnNonExisting{"fatalOnNonExisting", true, "Fatal message if calibrations not found on the CCDB"};

  // Parameters setting from json
  Configurable<LabeledArray<float>> bbParameters{"bbParameters",
                                                 {defaultParameters[0], nSpecies, nParameters, particleNames, parameterNames},
                                                 "Bethe Bloch parameters"};
  // Parameter setting from input file (including the ccdb)
  Configurable<LabeledArray<std::string>> fileParamBbPositive{"fileParamBbPositive",
                                                              {defaultOptions[0], nSpecies, nOptions, particleNames, optionNames},
                                                              "Input for the parametrization for positive particles. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};

  Configurable<LabeledArray<std::string>> fileParamBbNegative{"fileParamBbNegative",
                                                              {defaultOptions[0], nSpecies, nOptions, particleNames, optionNames},
                                                              "Input for the parametrization for negative particles. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};

  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> ccdbTimestamp{"ccdb-timestamp", -1, "timestamp of the object used to query in CCDB the detector response. If 0 the object corresponding to the run number is used, if < 0 the latest object is used"};

  bbParams bbPosEl{"PosEl"};
  bbParams bbPosMu{"PosMu"};
  bbParams bbPosPi{"PosPi"};
  bbParams bbPosKa{"PosKa"};
  bbParams bbPosPr{"PosPr"};
  bbParams bbPosDe{"PosDe"};
  bbParams bbPosTr{"PosTr"};
  bbParams bbPosHe{"PosHe"};
  bbParams bbPosAl{"PosAl"};

  bbParams bbNegEl{"NegEl"};
  bbParams bbNegMu{"NegMu"};
  bbParams bbNegPi{"NegPi"};
  bbParams bbNegKa{"NegKa"};
  bbParams bbNegPr{"NegPr"};
  bbParams bbNegDe{"NegDe"};
  bbParams bbNegTr{"NegTr"};
  bbParams bbNegHe{"NegHe"};
  bbParams bbNegAl{"NegAl"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <o2::track::PID::ID id, typename T>
  float BetheBlochLf(const T& track, const bbParams& params) const
  {
    static constexpr float invmass = 1.f / o2::track::pid_constants::sMasses2Z[id];
    static constexpr float charge = o2::track::pid_constants::sCharges[id];
    float corr = 0.f;
    if (params.postCorrection != nullptr) {
      corr = params.postCorrection->Eval(track.tpcInnerParam());
    }
    if (params.postCorrectionFun != nullptr) {
      corr = params.postCorrectionFun->Eval(track.tpcInnerParam());
    }
    if (params.isSimple) {
      return o2::tpc::BetheBlochAleph(track.tpcInnerParam() * invmass, params.bb1, params.bb2, params.bb3, params.bb4, params.bb5) + corr;
    }
    return params.mip * o2::tpc::BetheBlochAleph(track.tpcInnerParam() * invmass, params.bb1, params.bb2, params.bb3, params.bb4, params.bb5) * std::pow(charge, params.exp) + corr;
  }

  template <o2::track::PID::ID id, typename T>
  float BetheBlochResolutionLf(const T& track, const bbParams& params, const float bb) const
  {
    float corr = 1.f;
    if (params.postCorrectionSigma != nullptr) {
      corr = params.postCorrectionSigma->Eval(track.tpcInnerParam());
    }
    if (params.postCorrectionFunSigma != nullptr) {
      corr = params.postCorrectionFunSigma->Eval(track.tpcInnerParam());
    }
    return params.res * bb * corr;
  }

#define makeBetheBlochPerParticle(id, Particle)                                          \
  template <typename T>                                                                  \
  float BetheBlochPos##Particle(const T& track) const                                    \
  {                                                                                      \
    return BetheBlochLf<o2::track::PID::id>(track, bbPos##Particle);                     \
  }                                                                                      \
                                                                                         \
  template <typename T>                                                                  \
  float BetheBlochNeg##Particle(const T& track) const                                    \
  {                                                                                      \
    return BetheBlochLf<o2::track::PID::id>(track, bbNeg##Particle);                     \
  }                                                                                      \
  template <typename T>                                                                  \
  float BetheBlochResPos##Particle(const T& track, const float bb) const                 \
  {                                                                                      \
    return BetheBlochResolutionLf<o2::track::PID::id>(track, bbPos##Particle, bb);       \
  }                                                                                      \
  template <typename T>                                                                  \
  float BetheBlochResNeg##Particle(const T& track, const float bb) const                 \
  {                                                                                      \
    return BetheBlochResolutionLf<o2::track::PID::Electron>(track, bbNeg##Particle, bb); \
  }

  makeBetheBlochPerParticle(Electron, El);
  makeBetheBlochPerParticle(Muon, Mu);
  makeBetheBlochPerParticle(Pion, Pi);
  makeBetheBlochPerParticle(Kaon, Ka);
  makeBetheBlochPerParticle(Proton, Pr);
  makeBetheBlochPerParticle(Deuteron, De);
  makeBetheBlochPerParticle(Triton, Tr);
  makeBetheBlochPerParticle(Helium3, He);
  makeBetheBlochPerParticle(Alpha, Al);
#undef makeBetheBlochPerParticle

  void init(o2::framework::InitContext& initContext)
  {

    // Set up the CCDB
    const auto ts = ccdbTimestamp.value;
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    if (ts >= 0) {
      LOGP(info, "Initialising LF TPC PID response for fixed timestamp {}:", ts);
      ccdb->setTimestamp(ts);
    }
    if (!fatalOnNonExisting) {
      LOG(warning) << "Setting the CCDB to not fatal when the object is not found. This means that you have no garantuee that the parameters used are correct. Expert user only!!!";
      ccdb->setFatalWhenNull(false);
    }

#define InitPerParticle(Particle)                                                                          \
  if (doprocess##Particle || doprocessFull##Particle) {                                                    \
    LOG(info) << "Enabling " << #Particle;                                                                 \
    bbPos##Particle.init(#Particle, bbParameters, fileParamBbPositive, ccdb);                              \
    auto h = histos.add<TH1>(Form("%s", #Particle), "", kTH1F, {{10, 0, 10}});                             \
    h->SetBit(TH1::kIsAverage);                                                                            \
    h->SetBinContent(1, bbPos##Particle.bb1);                                                              \
    h->SetBinContent(2, bbPos##Particle.bb2);                                                              \
    h->SetBinContent(3, bbPos##Particle.bb3);                                                              \
    h->SetBinContent(4, bbPos##Particle.bb4);                                                              \
    h->SetBinContent(5, bbPos##Particle.bb5);                                                              \
    h->SetBinContent(6, bbPos##Particle.mip);                                                              \
    h->SetBinContent(7, bbPos##Particle.exp);                                                              \
    h->SetBinContent(8, bbPos##Particle.res);                                                              \
    h->SetBinContent(9, 1.f);                                                                              \
    bbNeg##Particle.init(#Particle, bbParameters, fileParamBbNegative, ccdb);                              \
    h = histos.add<TH1>(Form("Neg%s", #Particle), "", kTH1F, {{10, 0, 10}});                               \
    h->SetBit(TH1::kIsAverage);                                                                            \
    h->SetBinContent(1, bbNeg##Particle.bb1);                                                              \
    h->SetBinContent(2, bbNeg##Particle.bb2);                                                              \
    h->SetBinContent(3, bbNeg##Particle.bb3);                                                              \
    h->SetBinContent(4, bbNeg##Particle.bb4);                                                              \
    h->SetBinContent(5, bbNeg##Particle.bb5);                                                              \
    h->SetBinContent(6, bbNeg##Particle.mip);                                                              \
    h->SetBinContent(7, bbNeg##Particle.exp);                                                              \
    h->SetBinContent(8, bbNeg##Particle.res);                                                              \
    h->SetBinContent(9, 1.f);                                                                              \
  } else {                                                                                                 \
    LOG(info) << "Skipping " << #Particle;                                                                 \
    const bool requireTiny = isTableRequiredInWorkflow(initContext, Form("pidTPCLf%s", #Particle));        \
    const bool requireFull = isTableRequiredInWorkflow(initContext, Form("pidTPCLfFull%s", #Particle));    \
    if (requireTiny || requireFull) {                                                                      \
      LOG(fatal) << "Requested "                                                                           \
                 << #Particle << " table but not enabled in configuration: pidTPCLf" << #Particle          \
                 << " -> " << requireTiny << " (" << doprocess##Particle << "), pidTPCLfFull" << #Particle \
                 << " -> " << requireFull << " (" << doprocessFull##Particle << ")";                       \
    }                                                                                                      \
  }

    if (doprocessStandalone) { // If in standalone mode we enable the configuration of tables of interest
      LOG(info) << "Processing in standalone mode";
      doprocessFullPi.value = true;
      doprocessFullKa.value = true;
      doprocessFullPr.value = true;
    }

    InitPerParticle(El);
    InitPerParticle(Mu);
    InitPerParticle(Pi);
    InitPerParticle(Ka);
    InitPerParticle(Pr);
    InitPerParticle(De);
    InitPerParticle(Tr);
    InitPerParticle(He);
    InitPerParticle(Al);

    if (doprocessStandalone) { // If in standalone mode we disable after their configuration the process functions
      doprocessFullPi.value = false;
      doprocessFullKa.value = false;
      doprocessFullPr.value = false;
    }

#undef InitPerParticle
  }

#define makeProcess(Particle)                                                                       \
  void process##Particle(Colls const& collisions,                                                   \
                         soa::Join<Trks, aod::pidTPC##Particle> const& tracks,                      \
                         aod::BCsWithTimestamps const&)                                             \
  {                                                                                                 \
    LOG(debug) << "Filling table for particle: " << #Particle;                                      \
    tablePID##Particle.reserve(tracks.size());                                                      \
    if (bbParameters->get(#Particle, "Use default tiny") >= 1.5f) {                                 \
      for (auto const& trk : tracks) {                                                              \
        tablePID##Particle(trk.tpcNSigmaStore##Particle());                                         \
      }                                                                                             \
    } else {                                                                                        \
      if (!collisions.size()) {                                                                     \
        LOG(warn) << "No collisions in the data frame. Dummy PID table for " << #Particle;          \
        for (unsigned int i{0}; i < tracks.size(); ++i) {                                           \
          tablePID##Particle(aod::pidtpc_tiny::binning::underflowBin);                              \
        }                                                                                           \
        return;                                                                                     \
      }                                                                                             \
      bbPos##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb); \
      bbNeg##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb); \
      float bb = 0.f;                                                                               \
      float expSigma = 1.f;                                                                         \
      for (auto const& trk : tracks) {                                                              \
        if (!trk.hasTPC()) {                                                                        \
          tablePID##Particle(aod::pidtpc_tiny::binning::underflowBin);                              \
          continue;                                                                                 \
        }                                                                                           \
        if (skipTPCOnly) {                                                                          \
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {                                    \
            tablePID##Particle(aod::pidtpc_tiny::binning::underflowBin);                            \
            continue;                                                                               \
          }                                                                                         \
        }                                                                                           \
        if (trk.sign() > 0) {                                                                       \
          if (!bbPos##Particle.betheBlochSet) {                                                     \
            tablePID##Particle(trk.tpcNSigmaStore##Particle());                                     \
            continue;                                                                               \
          }                                                                                         \
          bb = BetheBlochPos##Particle(trk);                                                        \
          expSigma = BetheBlochResPos##Particle(trk, bb);                                           \
        } else {                                                                                    \
          if (!bbNeg##Particle.betheBlochSet) {                                                     \
            tablePID##Particle(trk.tpcNSigmaStore##Particle());                                     \
            continue;                                                                               \
          }                                                                                         \
          bb = BetheBlochNeg##Particle(trk);                                                        \
          expSigma = BetheBlochResNeg##Particle(trk, bb);                                           \
        }                                                                                           \
        aod::pidtpc_tiny::binning::packInTable((trk.tpcSignal() - bb) / expSigma,                   \
                                               tablePID##Particle);                                 \
      }                                                                                             \
    }                                                                                               \
  }                                                                                                 \
  PROCESS_SWITCH(lfTpcPid, process##Particle, "Produce a table for the " #Particle " hypothesis", false);

  makeProcess(El);
  makeProcess(Mu);
  makeProcess(Pi);
  makeProcess(Ka);
  makeProcess(Pr);
  makeProcess(De);
  makeProcess(Tr);
  makeProcess(He);
  makeProcess(Al);

#undef makeProcess

// Full tables
#define makeProcess(Particle)                                                                       \
  void processFull##Particle(Colls const& collisions,                                               \
                             soa::Join<Trks, aod::pidTPCFull##Particle> const& tracks,              \
                             aod::BCsWithTimestamps const&)                                         \
  {                                                                                                 \
    LOG(debug) << "Filling full table for particle: " << #Particle;                                 \
    tablePIDFull##Particle.reserve(tracks.size());                                                  \
    if (bbParameters->get(#Particle, "Use default full") >= 1.5f) {                                 \
      for (auto const& trk : tracks) {                                                              \
        tablePIDFull##Particle(trk.tpcExpSigma##Particle(), trk.tpcNSigma##Particle());             \
      }                                                                                             \
    } else {                                                                                        \
      if (!collisions.size()) {                                                                     \
        LOG(warn) << "No collisions in the data frame. Dummy PID table for " << #Particle;          \
        for (unsigned int i{0}; i < tracks.size(); ++i) {                                           \
          tablePIDFull##Particle(-999.f, -999.f);                                                   \
        }                                                                                           \
        return;                                                                                     \
      }                                                                                             \
      bbPos##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb); \
      bbNeg##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb); \
      float bb = 0.f;                                                                               \
      float expSigma = 1.f;                                                                         \
      for (auto const& trk : tracks) {                                                              \
        if (!trk.hasTPC()) {                                                                        \
          tablePIDFull##Particle(-999.f, -999.f);                                                   \
          continue;                                                                                 \
        }                                                                                           \
        if (skipTPCOnly) {                                                                          \
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {                                    \
            tablePIDFull##Particle(-999.f, -999.f);                                                 \
            continue;                                                                               \
          }                                                                                         \
        }                                                                                           \
        if (trk.sign() > 0) {                                                                       \
          if (!bbPos##Particle.betheBlochSet) {                                                     \
            tablePIDFull##Particle(trk.tpcExpSigma##Particle(), trk.tpcNSigma##Particle());         \
            continue;                                                                               \
          }                                                                                         \
          bb = BetheBlochPos##Particle(trk);                                                        \
          expSigma = BetheBlochResPos##Particle(trk, bb);                                           \
        } else {                                                                                    \
          if (!bbNeg##Particle.betheBlochSet) {                                                     \
            tablePIDFull##Particle(trk.tpcExpSigma##Particle(), trk.tpcNSigma##Particle());         \
            continue;                                                                               \
          }                                                                                         \
          bb = BetheBlochNeg##Particle(trk);                                                        \
          expSigma = BetheBlochResNeg##Particle(trk, bb);                                           \
        }                                                                                           \
        tablePIDFull##Particle(expSigma, (trk.tpcSignal() - bb) / expSigma);                        \
      }                                                                                             \
    }                                                                                               \
  }                                                                                                 \
  PROCESS_SWITCH(lfTpcPid, processFull##Particle, "Produce a full table for the " #Particle " hypothesis", false);

  makeProcess(El);
  makeProcess(Mu);
  makeProcess(Pi);
  makeProcess(Ka);
  makeProcess(Pr);
  makeProcess(De);
  makeProcess(Tr);
  makeProcess(He);
  makeProcess(Al);

#undef makeProcess

  // Full tables (independent on central calibrations)
  void processStandalone(Colls const& collisions,
                         Trks const& tracks,
                         aod::BCsWithTimestamps const&)
  {
    bool dummyPID = false;
    if (!collisions.size()) {
      LOG(warn) << "No collisions in the data frame. Dummy PID table for Pi-Ka-Pr";
      dummyPID = true;
    } else {
      tablePIDFullPi.reserve(tracks.size());
      bbPosPi.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
      bbNegPi.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
      tablePIDFullKa.reserve(tracks.size());
      bbPosKa.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
      bbNegKa.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
      tablePIDFullPr.reserve(tracks.size());
      bbPosPr.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
      bbNegPr.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    }
    float bb = 0.f;
    float expSigma = 1.f;
    for (auto const& trk : tracks) {
      if (!trk.hasTPC() || dummyPID) {
        tablePIDFullPi(-999.f, -999.f);
        tablePIDFullKa(-999.f, -999.f);
        tablePIDFullPr(-999.f, -999.f);
        continue;
      }
      if (skipTPCOnly) {
        if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
          tablePIDFullPi(-999.f, -999.f);
          tablePIDFullKa(-999.f, -999.f);
          tablePIDFullPr(-999.f, -999.f);
          continue;
        }
      }
      if (trk.sign() > 0) {
        bb = BetheBlochPosPi(trk);
        expSigma = BetheBlochResPosPi(trk, bb);
      } else {
        bb = BetheBlochNegPi(trk);
        expSigma = BetheBlochResNegPi(trk, bb);
      }
      tablePIDFullPi(expSigma, (trk.tpcSignal() - bb) / expSigma);

      if (trk.sign() > 0) {
        bb = BetheBlochPosKa(trk);
        expSigma = BetheBlochResPosKa(trk, bb);
      } else {
        bb = BetheBlochNegKa(trk);
        expSigma = BetheBlochResNegKa(trk, bb);
      }
      tablePIDFullKa(expSigma, (trk.tpcSignal() - bb) / expSigma);

      if (trk.sign() > 0) {
        bb = BetheBlochPosPr(trk);
        expSigma = BetheBlochResPosPr(trk, bb);
      } else {
        bb = BetheBlochNegPr(trk);
        expSigma = BetheBlochResNegPr(trk, bb);
      }
      tablePIDFullPr(expSigma, (trk.tpcSignal() - bb) / expSigma);
    }
  }
  PROCESS_SWITCH(lfTpcPid, processStandalone, "Produce full tables in a standalone way for Pi-Ka-Pr", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<lfTpcPid>(cfgc)}; }
