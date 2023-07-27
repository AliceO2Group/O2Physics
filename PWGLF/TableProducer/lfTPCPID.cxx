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

  // Utility parameters for the usage
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
      LOG(error) << "bbParams `" << name << "` :: The vector of Bethe-Bloch parameters has the wrong size";
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
    return true;
  }

  ///
  /// Set the values of the BetheBloch from a configurable array
  bool setValues(const char* particle, const Configurable<LabeledArray<float>>& p)
  {
    if (p->get(particle, "Set parameters") < 1.5f) {
      LOG(info) << "bbParams `" << name << "` :: Using default parameters for " << particle << " as entry `Set parameters` " << p->get(particle, "Set parameters") << " < 1.5";
      return false;
    }
    std::vector<float> v{p->get(particle, "bb1"),
                         p->get(particle, "bb2"),
                         p->get(particle, "bb3"),
                         p->get(particle, "bb4"),
                         p->get(particle, "bb5"),
                         p->get(particle, "MIP value"),
                         p->get(particle, "Charge exponent"),
                         p->get(particle, "Resolution")};
    LOG(info) << "bbParams `" << name << "` :: Setting custom Bethe-Bloch parameters for mass hypothesis " << particle;
    return setValues(v);
  }

  ///
  /// Set the values of the BetheBloch from TH1F
  bool setValues(TH1F* h)
  {
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
      LOG(error) << "The input histogram of Bethe-Bloch parameters has the wrong size " << n << " while expecting " << h->GetNbinsX();
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
      LOG(error) << "The input histogram of Bethe-Bloch parameters has the wrong size " << n << " while expecting " << h->GetNbinsX();
      return false;
    }
    LOG(info) << "bbParams `" << name << "` :: Setting simple custom Bethe-Bloch parameters from histogram " << h->GetName();
    return setValues(v);
  }

  ///
  /// Set values from an input file
  bool setValues(TFile* f)
  {
    if (!f) {
      LOG(fatal) << "The input file is not valid";
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
        return false;
      }
    }
    LOG(info) << "bbParams `" << name << "` :: Setting parameters from TH1F " << h->GetName() << " in file " << f->GetName();
    return setValues(h);
  }

  ///
  /// Set values from a configuration. In this case also the post calibration is checked
  bool init(const Configurable<std::vector<std::string>>& cfg, o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj)
  {
    LOG(info) << "bbParams `" << name << "` :: initializing from configurable '" << cfg.name << "' of size " << cfg.value.size();
    if (cfg.value.size() != 3) {
      LOG(fatal) << "bbParams `" << name << "` :: The input configurable has the wrong size " << cfg.value.size() << " while expecting 3";
    }
    const std::string bb = cfg.value.at(0);
    const std::string post = cfg.value.at(1);
    const std::string simple = cfg.value.at(2);
    // First we check the post calib
    if (post.size() > 1) {
      LOG(info) << "bbParams `" << name << "` :: Loading parameters from configurable '" << cfg.name << "' with value '" << post << "'";
      std::string s = post;
      if (s.rfind("ccdb://", 0) == 0) {
        s.replace(0, 7, "");
        ccdbPathPost = s;
        if (ccdbObj->getTimestamp() == 0) { // If the timestamp is 0 we expect to get the timestamp from the run number
          takePostFromCcdb = true;
          LOG(info) << "bbParams `" << name << "` :: For post calib asked to query the parameters from the CCDB and got timestamp " << ccdbObj->getTimestamp() << " -> will take the object corresponding to the run number";
        } else {
          LOG(info) << "bbParams `" << name << "` :: For post calib fetching parameters from CCDB (only once) using timestamp " << ccdbObj->getTimestamp() << " and path " << s;
          TList* postL = ccdbObj->get<TList>(s);
          postCorrection = static_cast<TGraph*>(postL->FindObject(Form("%s_postCorrection", name.c_str())));
          postCorrectionSigma = static_cast<TGraph*>(postL->FindObject(Form("%s_postCorrectionSigma", name.c_str())));
        }
      }
    }

    // Check the BetheBloch parameters
    if (bb.size() <= 1) {
      return false;
    }
    if (simple == "true") { // Check if the parametrization to be used is simple
      LOG(info) << "bbParams `" << name << "` :: using simple parametrization";
      isSimple = true;
    } else {
      LOG(info) << "bbParams `" << name << "` :: using full parametrization";
    }
    LOG(info) << "bbParams `" << name << "` :: Loading parameters from configurable '" << cfg.name << "' with value '" << bb << "'";
    std::string s = bb;
    if (s.rfind("ccdb://", 0) == 0) {
      s.replace(0, 7, "");
      ccdbPath = s;
      if (ccdbObj->getTimestamp() == 0) { // If the timestamp is 0 we expect to get the timestamp from the run number
        takeFromCcdb = true;
        LOG(info) << "bbParams `" << name << "` :: Asked to query the parameters from the CCDB and got timestamp " << ccdbObj->getTimestamp() << " -> will take the object corresponding to the run number";
        return false;
      }
      LOG(info) << "bbParams `" << name << "` :: Fetching parameters from CCDB (only once) using timestamp " << ccdbObj->getTimestamp() << " and path " << s;
      TH1F* h = ccdbObj->get<TH1F>(s);
      return setValues(h);
    }
    return setValues(TFile::Open(bb.c_str(), "READ"));
  }

  /// @brief Function to update the values of the parameters from the CCDB
  /// @param ccdbObj Managare of the CCDB
  /// @param timestamp Timestamp to ask the new parameters
  /// @return false if not succesfull
  bool updateValues(aod::BCsWithTimestamps::iterator const& bunchCrossing, o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj)
  {
    // Check that the last updated number is different
    if (lastRunNumber == bunchCrossing.runNumber()) {
      LOG(debug) << "bbParams `" << name << "` :: Not updating parameters of " << name << " from run number " << lastRunNumber << " as they are already up to date";
      return false;
    }
    // First check the post calib
    if (takePostFromCcdb) {
      LOG(info) << "bbParams `" << name << "` :: Updating post calib from run number " << lastRunNumber << " to " << bunchCrossing.runNumber() << ". Taking them from CCDB path '" << ccdbPathPost << "' with timestamp " << bunchCrossing.timestamp();
      TList* postL = ccdbObj->getForTimeStamp<TList>(ccdbPathPost, bunchCrossing.timestamp());
      postCorrection = static_cast<TGraph*>(postL->FindObject(Form("%s_postCorrection", name.c_str())));
      postCorrectionSigma = static_cast<TGraph*>(postL->FindObject(Form("%s_postCorrectionSigma", name.c_str())));
    }
    // Secondly we check the BB parameters
    if (!takeFromCcdb) {
      LOG(debug) << "bbParams `" << name << "` :: Not taking parameters from CCDB";
      return false;
    }
    LOG(info) << "bbParams `" << name << "` :: Updating parameters of " << name << " from run number " << lastRunNumber << " to " << bunchCrossing.runNumber() << ". Taking them from CCDB path '" << ccdbPath << "' with timestamp " << bunchCrossing.timestamp();
    lastRunNumber = bunchCrossing.runNumber();
    return setValues(ccdbObj->getForTimeStamp<TH1F>(ccdbPath, bunchCrossing.timestamp()));
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

  // Parameters setting from json
  Configurable<LabeledArray<float>> bbParameters{"bbParameters",
                                                 {defaultParameters[0], nSpecies, nParameters, particleNames, parameterNames},
                                                 "Bethe Bloch parameters"};
  // Parameter setting from input file (including the ccdb)
  Configurable<std::vector<std::string>> fileParamBbEl{"fileParamBbEl",
                                                       {"", "", ""},
                                                       "Input for the parametrization for electrons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbMu{"fileParamBbMu",
                                                       {"", "", ""},
                                                       "Input for the parametrization for muons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbPi{"fileParamBbPi",
                                                       {"", "", ""},
                                                       "Input for the parametrization for pions 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbKa{"fileParamBbKa",
                                                       {"", "", ""},
                                                       "Input for the parametrization for kaons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbPr{"fileParamBbPr",
                                                       {"", "", ""},
                                                       "Input for the parametrization for protons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbDe{"fileParamBbDe",
                                                       {"", "", ""},
                                                       "Input for the parametrization for deuterons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbTr{"fileParamBbTr",
                                                       {"", "", ""},
                                                       "Input for the parametrization for tritons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbHe{"fileParamBbHe",
                                                       {"", "", ""},
                                                       "Input for the parametrization for helium3 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbAl{"fileParamBbAl",
                                                       {"", "", ""},
                                                       "Input for the parametrization for helium4 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegEl{"fileParamBbNegEl",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative electrons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegMu{"fileParamBbNegMu",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative muons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegPi{"fileParamBbNegPi",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative pions 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegKa{"fileParamBbNegKa",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative kaons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegPr{"fileParamBbNegPr",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative protons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegDe{"fileParamBbNegDe",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative deuterons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegTr{"fileParamBbNegTr",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative tritons 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegHe{"fileParamBbNegHe",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative helium3 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};
  Configurable<std::vector<std::string>> fileParamBbNegAl{"fileParamBbNegAl",
                                                          {"", "", ""},
                                                          "Input for the parametrization for negative helium4 1) BB 2) Post 3) Simple. If empty using the default/json values. Can be a CCDB path if the string starts with ccdb://"};

  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> ccdbTimestamp{"ccdb-timestamp", -1, "timestamp of the object used to query in CCDB the detector response. If 0 the object corresponding to the run number is used, if < 0 the latest object is used"};

  bbParams bbEl{"El"};
  bbParams bbMu{"Mu"};
  bbParams bbPi{"Pi"};
  bbParams bbKa{"Ka"};
  bbParams bbPr{"Pr"};
  bbParams bbDe{"De"};
  bbParams bbTr{"Tr"};
  bbParams bbHe{"He"};
  bbParams bbAl{"Al"};

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
  float BetheBlochLf(const T& track, const bbParams& params)
  {
    static constexpr float invmass = 1.f / o2::track::pid_constants::sMasses2Z[id];
    static constexpr float charge = o2::track::pid_constants::sCharges[id];
    float corr = 0.f;
    if (params.postCorrection != nullptr) {
      corr = params.postCorrection->Eval(track.tpcInnerParam());
    }
    if (params.isSimple) {
      return o2::tpc::BetheBlochAleph(track.tpcInnerParam() * invmass, params.bb1, params.bb2, params.bb3, params.bb4, params.bb5) - corr;
    }
    return params.mip * o2::tpc::BetheBlochAleph(track.tpcInnerParam() * invmass, params.bb1, params.bb2, params.bb3, params.bb4, params.bb5) * std::pow(charge, params.exp) - corr;
  }

  template <typename T>
  float BetheBlochEl(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Electron>(track, bbEl);
  }
  template <typename T>
  float BetheBlochMu(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Muon>(track, bbMu);
  }
  template <typename T>
  float BetheBlochPi(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Pion>(track, bbPi);
  }
  template <typename T>
  float BetheBlochKa(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Kaon>(track, bbKa);
  }
  template <typename T>
  float BetheBlochPr(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Proton>(track, bbPr);
  }
  template <typename T>
  float BetheBlochDe(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Deuteron>(track, bbDe);
  }
  template <typename T>
  float BetheBlochTr(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Triton>(track, bbTr);
  }
  template <typename T>
  float BetheBlochHe(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Helium3>(track, bbHe);
  }
  template <typename T>
  float BetheBlochAl(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Alpha>(track, bbAl);
  }

  template <typename T>
  float BetheBlochNegEl(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Electron>(track, bbNegEl);
  }
  template <typename T>
  float BetheBlochNegMu(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Muon>(track, bbNegMu);
  }
  template <typename T>
  float BetheBlochNegPi(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Pion>(track, bbNegPi);
  }
  template <typename T>
  float BetheBlochNegKa(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Kaon>(track, bbNegKa);
  }
  template <typename T>
  float BetheBlochNegPr(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Proton>(track, bbNegPr);
  }
  template <typename T>
  float BetheBlochNegDe(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Deuteron>(track, bbNegDe);
  }
  template <typename T>
  float BetheBlochNegTr(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Triton>(track, bbNegTr);
  }
  template <typename T>
  float BetheBlochNegHe(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Helium3>(track, bbNegHe);
  }
  template <typename T>
  float BetheBlochNegAl(const T& track)
  {
    return BetheBlochLf<o2::track::PID::Alpha>(track, bbNegAl);
  }

  template <o2::track::PID::ID id, typename T>
  float BetheBlochResolutionLf(const T& track, const bbParams& params, const float bb)
  {
    float corr = 1.f;
    if (params.postCorrectionSigma != nullptr) {
      corr = params.postCorrectionSigma->Eval(track.tpcInnerParam());
    }
    return params.res * bb * corr;
  }

  template <typename T>
  float BetheBlochResEl(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Electron>(track, bbEl, bb);
  }
  template <typename T>
  float BetheBlochResMu(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Muon>(track, bbMu, bb);
  }
  template <typename T>
  float BetheBlochResPi(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Pion>(track, bbPi, bb);
  }
  template <typename T>
  float BetheBlochResKa(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Kaon>(track, bbKa, bb);
  }
  template <typename T>
  float BetheBlochResPr(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Proton>(track, bbPr, bb);
  }
  template <typename T>
  float BetheBlochResDe(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Deuteron>(track, bbDe, bb);
  }
  template <typename T>
  float BetheBlochResTr(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Triton>(track, bbTr, bb);
  }
  template <typename T>
  float BetheBlochResHe(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Helium3>(track, bbHe, bb);
  }
  template <typename T>
  float BetheBlochResAl(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Alpha>(track, bbAl, bb);
  }
  template <typename T>
  float BetheBlochResNegEl(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Electron>(track, bbEl, bb);
  }
  template <typename T>
  float BetheBlochResNegMu(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Muon>(track, bbMu, bb);
  }
  template <typename T>
  float BetheBlochResNegPi(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Pion>(track, bbPi, bb);
  }
  template <typename T>
  float BetheBlochResNegKa(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Kaon>(track, bbKa, bb);
  }
  template <typename T>
  float BetheBlochResNegPr(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Proton>(track, bbPr, bb);
  }
  template <typename T>
  float BetheBlochResNegDe(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Deuteron>(track, bbDe, bb);
  }
  template <typename T>
  float BetheBlochResNegTr(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Triton>(track, bbTr, bb);
  }
  template <typename T>
  float BetheBlochResNegHe(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Helium3>(track, bbHe, bb);
  }
  template <typename T>
  float BetheBlochResNegAl(const T& track, const float bb)
  {
    return BetheBlochResolutionLf<o2::track::PID::Alpha>(track, bbAl, bb);
  }

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

#define InitPerParticle(Particle)                                                                          \
  if (doprocess##Particle || doprocessFull##Particle) {                                                    \
    LOG(info) << "Enabling " << #Particle;                                                                 \
    bb##Particle.setValues(#Particle, bbParameters);                                                       \
    bb##Particle.init(fileParamBb##Particle, ccdb);                                                        \
    bbNeg##Particle.setValues(#Particle, bbParameters);                                                    \
    bbNeg##Particle.init(fileParamBbNeg##Particle, ccdb);                                                  \
    auto h = histos.add<TH1>(Form("%s", #Particle), "", kTH1F, {{10, 0, 10}});                             \
    h->SetBit(TH1::kIsAverage);                                                                            \
    h->SetBinContent(1, bb##Particle.bb1);                                                                 \
    h->SetBinContent(2, bb##Particle.bb2);                                                                 \
    h->SetBinContent(3, bb##Particle.bb3);                                                                 \
    h->SetBinContent(4, bb##Particle.bb4);                                                                 \
    h->SetBinContent(5, bb##Particle.bb5);                                                                 \
    h->SetBinContent(6, bb##Particle.mip);                                                                 \
    h->SetBinContent(7, bb##Particle.exp);                                                                 \
    h->SetBinContent(8, bb##Particle.res);                                                                 \
    h->SetBinContent(9, 1.f);                                                                              \
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

#define makeProcess(Particle)                                                                                            \
  void process##Particle(Colls const& collisions,                                                                        \
                         soa::Join<Trks, aod::pidTPC##Particle> const& tracks,                                           \
                         aod::BCsWithTimestamps const&)                                                                  \
  {                                                                                                                      \
    LOG(debug) << "Filling table for particle: " << #Particle;                                                           \
    tablePID##Particle.reserve(tracks.size());                                                                           \
    if (bbParameters->get(#Particle, "Use default tiny") >= 1.5f) {                                                      \
      for (auto const& trk : tracks) {                                                                                   \
        tablePID##Particle(trk.tpcNSigmaStore##Particle());                                                              \
      }                                                                                                                  \
    } else {                                                                                                             \
      bb##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);                         \
      bbNeg##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);                      \
      float bb = 0.f;                                                                                                    \
      for (auto const& trk : tracks) {                                                                                   \
        if (!trk.hasTPC()) {                                                                                             \
          tablePID##Particle(aod::pidtpc_tiny::binning::underflowBin);                                                   \
          continue;                                                                                                      \
        }                                                                                                                \
        if (skipTPCOnly) {                                                                                               \
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {                                                         \
            tablePID##Particle(aod::pidtpc_tiny::binning::underflowBin);                                                 \
            continue;                                                                                                    \
          }                                                                                                              \
        }                                                                                                                \
        if (trk.sign() > 0) {                                                                                            \
          bb = BetheBloch##Particle(trk);                                                                                \
        } else {                                                                                                         \
          bb = BetheBlochNeg##Particle(trk);                                                                             \
        }                                                                                                                \
        aod::pidutils::packInTable<aod::pidtpc_tiny::binning>((trk.tpcSignal() - bb) / BetheBlochRes##Particle(trk, bb), \
                                                              tablePID##Particle);                                       \
      }                                                                                                                  \
    }                                                                                                                    \
  }                                                                                                                      \
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
      bb##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);    \
      bbNeg##Particle.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb); \
      float bb = 0.f;                                                                               \
      float expSigma = 1.f;                                                                         \
      for (auto const& trk : tracks) {                                                              \
        if (skipTPCOnly) {                                                                          \
          if (!trk.hasTPC()) {                                                                      \
            tablePIDFull##Particle(-999.f, -999.f);                                                 \
            continue;                                                                               \
          }                                                                                         \
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {                                    \
            tablePIDFull##Particle(-999.f, -999.f);                                                 \
            continue;                                                                               \
          }                                                                                         \
        }                                                                                           \
        if (trk.sign() > 0) {                                                                       \
          bb = BetheBloch##Particle(trk);                                                           \
          expSigma = BetheBlochRes##Particle(trk, bb);                                              \
        } else {                                                                                    \
          bb = BetheBlochNeg##Particle(trk);                                                        \
          expSigma = BetheBlochResNeg##Particle(trk, bb);                                           \
        }                                                                                           \
        tablePIDFull##Particle(expSigma,                                                            \
                               (trk.tpcSignal() - bb) / expSigma);                                  \
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
    tablePIDFullPi.reserve(tracks.size());
    bbPi.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    bbNegPi.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    tablePIDFullKa.reserve(tracks.size());
    bbKa.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    bbNegKa.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    tablePIDFullPr.reserve(tracks.size());
    bbPr.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    bbNegPr.updateValues(collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>(), ccdb);
    float bb = 0.f;
    float expSigma = 1.f;
    for (auto const& trk : tracks) {
      if (skipTPCOnly) {
        if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
          tablePIDFullPi(-999.f, -999.f);
          tablePIDFullKa(-999.f, -999.f);
          tablePIDFullPr(-999.f, -999.f);
          continue;
        }
        if (!trk.hasTPC()) {
          tablePIDFullPi(-999.f, -999.f);
          tablePIDFullKa(-999.f, -999.f);
          tablePIDFullPr(-999.f, -999.f);
        }
      }
      if (trk.sign() > 0) {
        bb = BetheBlochPi(trk);
        expSigma = BetheBlochResPi(trk, bb);
      } else {
        bb = BetheBlochNegPi(trk);
        expSigma = BetheBlochResNegPi(trk, bb);
      }
      tablePIDFullPi(expSigma, (trk.tpcSignal() - bb) / expSigma);

      if (trk.sign() > 0) {
        bb = BetheBlochKa(trk);
        expSigma = BetheBlochResKa(trk, bb);
      } else {
        bb = BetheBlochNegKa(trk);
        expSigma = BetheBlochResNegKa(trk, bb);
      }
      tablePIDFullKa(expSigma, (trk.tpcSignal() - bb) / expSigma);

      if (trk.sign() > 0) {
        bb = BetheBlochPr(trk);
        expSigma = BetheBlochResPr(trk, bb);
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
