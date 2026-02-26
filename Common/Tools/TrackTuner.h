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

/// \file TrackTuner.h
/// \brief Helper class to "smear" track parameters in MC simulations
/// \author Andrea Rossi (andrea.rossi@cern.ch), INFN Padova, Italy
/// \author Himanshu Sharma (himanshu.sharma@cern.ch), INFN Padova, Italy
/// \author Mattia Faggin (mattia.faggin@cern.ch), CERN

#ifndef COMMON_TOOLS_TRACKTUNER_H_
#define COMMON_TOOLS_TRACKTUNER_H_

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h> // FIXME: remove
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TObject.h>

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace o2::aod
{
namespace track_tuner
{
DECLARE_SOA_COLUMN(TunedQOverPt, tunedQOverPt, float);

/// configuration source
enum configSource : int { InputString = 1,
                          Configurables };
} // namespace track_tuner

DECLARE_SOA_TABLE(TrackTunerTable, "AOD", "TRACKTUNERTABLE", //!
                  track_tuner::TunedQOverPt);
} // namespace o2::aod

struct TrackTuner : o2::framework::ConfigurableGroup {

  std::string prefix = "trackTuner"; // JSON group name
  o2::framework::Configurable<bool> cfgDebugInfo{"debugInfo", false, "Flag to switch on the debug printout"};
  o2::framework::Configurable<bool> cfgUpdateTrackDCAs{"updateTrackDCAs", false, "Flag to enable the DCA smearing"};
  o2::framework::Configurable<bool> cfgUpdateTrackCovMat{"updateTrackCovMat", false, "Flag to enable the DCA covariance-matrix smearing"};
  o2::framework::Configurable<bool> cfgUpdateCurvature{"updateCurvature", false, "Flag to enable the Q/Pt smearing after the propagation to the production point"};
  o2::framework::Configurable<bool> cfgUpdateCurvatureIU{"updateCurvatureIU", false, "Flag to enable the Q/Pt smearing before the propagation to the production point"};
  o2::framework::Configurable<bool> cfgUpdatePulls{"updatePulls", false, "Flag to enable the pulls smearing"};
  o2::framework::Configurable<bool> cfgIsInputFileFromCCDB{"isInputFileFromCCDB", false, "True: files from CCDB; False: fils from local path (debug)"};
  o2::framework::Configurable<std::string> cfgPathInputFile{"pathInputFile", "", "Path to file containing DCAxy, DCAz graphs from data and MC"};
  o2::framework::Configurable<std::string> cfgNameInputFile{"nameInputFile", "", "Name of the file containing DCAxy, DCAz graphs from data and MC"};
  o2::framework::Configurable<std::string> cfgPathFileQoverPt{"pathFileQoverPt", "", "Path to file containing Q/Pt correction graphs from data and MC"};
  o2::framework::Configurable<std::string> cfgNameFileQoverPt{"nameFileQoverPt", "", "Name of file containing Q/Pt correction graphs from data and MC"};
  o2::framework::Configurable<bool> cfgUsePvRefitCorrections{"usePvRefitCorrections", false, "Flag to establish whether to use corrections obtained with or w/o PV refit"};
  o2::framework::Configurable<float> cfgQOverPtMC{"qOverPtMC", -1., "Scaling factor on q/pt of MC"};
  o2::framework::Configurable<float> cfgQOverPtData{"qOverPtData", -1., "Scaling factor on q/pt of data"};
  o2::framework::Configurable<int> cfgNPhiBins{"nPhiBins", 0, "Number of phi bins"};
  o2::framework::Configurable<bool> cfgAutoDetectDcaCalib{"autoDetectDcaCalib", false, "Flag to enable the dca-calibration file autodetect from CCDB (list of predefined cases)"};
  ///////////////////////////////
  /// parameters to be configured
  bool debugInfo = false;
  bool updateTrackDCAs = false; // To update the track DCAs;
  bool updateTrackCovMat = false;
  bool updateCurvature = false;
  bool updateCurvatureIU = false; // To update the track parameter Q/Pt in trackIU table, particularly used for V0 mass width dependence on Q/Pt
  bool updatePulls = false;
  bool isInputFileFromCCDB = false;   // query input file from CCDB or local folder
  std::string pathInputFile = "";     // Path to file containing DCAxy, DCAz graphs from data and MC
  std::string nameInputFile = "";     // Name of the file containing DCAxy, DCAz graphs from data and MC
  std::string pathFileQoverPt = "";   // Path to file containing Q/Pt correction graphs from data and MC (only one proxy provided, i.e. D0 sigma graphs from data and MC)
  std::string nameFileQoverPt = "";   // file name containing Q/Pt correction graphs from data and MC
  bool usePvRefitCorrections = false; // establish whether to use corrections obtained with or w/o PV refit
  float qOverPtMC = -1.;              // 1/pt MC
  float qOverPtData = -1.;            // 1/pt data
  bool autoDetectDcaCalib = false;    // enable automatic pick-up of dca calibration files from CCDB (list of predefined cases)
  ///////////////////////////////
  bool isConfigFromString = false;
  bool isConfigFromConfigurables = false;
  int nPhiBins = 1;
  int runNumber = 0; // first run number considered in analysis (useful only if autoDetectDcaCalib = true)
  bool areGraphsConfigured = false;
  std::string outputString = "";

  std::map<std::string, std::string> metadata;

  std::vector<std::unique_ptr<TGraphErrors>> grDcaXYResVsPtPionMC;
  std::vector<std::unique_ptr<TGraphErrors>> grDcaXYResVsPtPionData;

  std::vector<std::unique_ptr<TGraphErrors>> grDcaZResVsPtPionMC;
  std::vector<std::unique_ptr<TGraphErrors>> grDcaZResVsPtPionData;

  std::vector<std::unique_ptr<TGraphErrors>> grDcaXYMeanVsPtPionMC;
  std::vector<std::unique_ptr<TGraphErrors>> grDcaXYMeanVsPtPionData;

  std::vector<std::unique_ptr<TGraphErrors>> grDcaZMeanVsPtPionMC;
  std::vector<std::unique_ptr<TGraphErrors>> grDcaZMeanVsPtPionData;

  std::unique_ptr<TGraphErrors> grOneOverPtPionMC;   // MC
  std::unique_ptr<TGraphErrors> grOneOverPtPionData; // Data

  std::vector<std::unique_ptr<TGraphErrors>> grDcaXYPullVsPtPionMC;
  std::vector<std::unique_ptr<TGraphErrors>> grDcaXYPullVsPtPionData;

  std::vector<std::unique_ptr<TGraphErrors>> grDcaZPullVsPtPionMC;
  std::vector<std::unique_ptr<TGraphErrors>> grDcaZPullVsPtPionData;

  /// @brief Function to initialize the run number to that of the 1st considered bunch crossing (useful only if autoDetectDcaCalib = true)
  void setRunNumber(int n)
  {
    runNumber = n;
  }

  /// @brief Function doing a few sanity-checks on the configurations
  void checkConfig()
  {
    /// check configuration source
    if (isConfigFromString && isConfigFromConfigurables) {
      LOG(fatal) << " [ isConfigFromString==kTRUE and isConfigFromConfigurables==kTRUE ] Configuration done both via string and via configurables -> Only one of them can be set to kTRUE at once! Please refer to the trackTuner documentation.";
    }
    /// check Q/pt update
    if ((updateCurvatureIU) && (updateCurvature)) {
      LOG(fatal) << " [ updateCurvatureIU==kTRUE and updateCurvature==kTRUE ] -> Only one of them can be set to kTRUE at once! Please refer to the trackTuner documentation.";
    }
  }

  void getPathInputFileAutomaticFromCCDB()
  {

    /// check: no CCDB autodetection if the desired input file is not in CCDB
    if (!isInputFileFromCCDB) {
      LOG(fatal) << "[TrackTuner::getPathInputFileAutomaticFromCCDB] Trying to auto detect the dca calibration file from CCDB, but you ask the input file to not come from CCDB (isInputFileFromCCDB==" << isInputFileFromCCDB << "). Fix it!";
    }
    /// check that the run number has been already properly set
    if (runNumber == 0) {
      LOG(fatal) << "[TrackTuner::getPathInputFileAutomaticFromCCDB] runNumber==" << runNumber << ", automatic detection of dca calibration file from CCDB not possible. Did you call the function TrackTuner::setrunNumber()?";
    }
    /// check than the number of phi bins for the track tuner calibrations is 24
    if (nPhiBins != 24) {
      LOG(fatal) << "[TrackTuner::getPathInputFileAutomaticFromCCDB] nPhiBins==" << nPhiBins << ", but the automatic detection of dca calibration file from CCDB is supported only for nPhiBins==24. Either put nPhiBins=24, or disable the auto-detection (autoDetectDcaCalib=false)";
    }

    pathInputFile = "invalid";

    ///////////////////////////////////////////////////////////////////////////
    ///                                                                     ///
    ///   establish some pre-defined cases based only on the run numbers    ///
    ///                                                                     ///
    ///////////////////////////////////////////////////////////////////////////
    LOG(info) << "";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++   TrackTuner configuration                                                                                                          +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++   Autodetect mode activated for the DCA calibration files                                                                           +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++   The DCA calibration files are picked-up from CCDB based on the analysed run number                                                +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++   NB: only the number of the first analysed run is considered to configure the TrackTuner object                                    +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                                                                                                                                     +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++   Supported cases:                                                                                                                  +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++      [CASE 1]: pp, 13.6 TeV 2022, 2023: CCDB path Users/m/mfaggin/test/inputsTrackTuner/pp2023/pass4/vsPhi                          +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                Run list: (520259 (LHC22f) <= runNumber <= 529691 (LHC22t)) || (534998 (LHC23zc) <= runNumber <= 543113 (LHC23zw))   +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                                                                                                                                     +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++      [CASE 2]: Pb-Pb, 5.34 TeV 2022, 2023, 2024: CCDB path Users/m/mfaggin/test/inputsTrackTuner/PbPb2023/apass4/vsPhi              +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                Run list: (529397 <= runNumber <= 529418 (LHC22o)) || (543437 (LHC23zx) <= runNumber <= 545367 (LHC23zzo))           +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                                                                                                                                     +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++      [CASE 3]: pp, 13.6 TeV 2024: CCDB path Users/m/mfaggin/test/inputsTrackTuner/pp2024/pass1_minBias/vsPhi                        +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                Run list: 549559 (LHC24ac) <= runNumber && runNumber <= 558807 (LHC24ao)                                             +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++      [CASE 4]: OO, 5.36 TeV 2025, period LHC25ae: CCDB path Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25ae                        +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                Run list: 564356 <= runNumber && runNumber <= 564445                                                                 +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++      [CASE 5]: OO, 5.36 TeV 2025, period LHC25af: CCDB path Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25af                        +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++                Run list: 564468 <= runNumber && runNumber <= 564472                                                                 +++";
    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    LOG(info) << "";

    LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]: current run number = " << runNumber;

    if ((520259 <= runNumber && runNumber <= 529691) || (534998 <= runNumber && runNumber <= 543113)) {
      ///
      ///   [CASE 1]: pp, 13.6 TeV 2022, 2023: CCDB path Users/m/mfaggin/test/inputsTrackTuner/pp2023/pass4/vsPhi
      ///             Run list: (520259 (LHC22f) <= runNumber <= 529691 (LHC22t)) || (534998 (LHC23zc) <= runNumber <= 543113 (LHC23zw))
      ///
      pathInputFile = "Users/m/mfaggin/test/inputsTrackTuner/pp2023/pass4/vsPhi";
      LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]:   >>>   pp, 13.6 TeV 2022, 2023: CCDB path " << pathInputFile;
      LOG(info) << "                                                   >>>   Run list: (520259 (LHC22f) <= runNumber <= 529691 (LHC22t)) || (534998 (LHC23zc) <= runNumber <= 543113 (LHC23zw))";
    } else if ((529397 <= runNumber && runNumber <= 529418) || (543437 <= runNumber && runNumber <= 545367)) {
      ///
      ///   [CASE 2]: Pb-Pb, 5.34 TeV 2022, 2023, 2024: CCDB path Users/m/mfaggin/test/inputsTrackTuner/PbPb2023/apass4/vsPhi
      ///             Run list: (529397 <= runNumber <= 529418 (LHC22o)) || (543437 (LHC23zx) <= runNumber <= 545367 (LHC23zzo))
      ///
      pathInputFile = "Users/m/mfaggin/test/inputsTrackTuner/PbPb2023/apass4/vsPhi";
      LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]:   >>>   Pb-Pb, 5.34 TeV 2022, 2023, 2024: CCDB path " << pathInputFile;
      LOG(info) << "                                                   >>>   Run list: (529397 <= runNumber <= 529418 (LHC22o)) || (543437 (LHC23zx) <= runNumber <= 545367 (LHC23zzo))";
    } else if (549559 <= runNumber && runNumber <= 558807) {
      ///
      ///   [CASE 3]: pp, 13.6 TeV 2024: CCDB path Users/m/mfaggin/test/inputsTrackTuner/pp2024/pass1_minBias/vsPhi
      ///             Run list: 549559 (LHC24ac) <= runNumber && runNumber <= 558807 (LHC24ao)
      ///
      pathInputFile = "Users/m/mfaggin/test/inputsTrackTuner/pp2024/pass1_minBias/vsPhi";
      LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]:   >>>   pp, 13.6 TeV 2024: CCDB path " << pathInputFile;
      LOG(info) << "                                                   >>>   Run list: 549559 (LHC24ac) <= runNumber && runNumber <= 558807 (LHC24ao)";
    } else if (564356 <= runNumber && runNumber <= 564445) {
      ///
      ///   [CASE 4]: OO, 5.36 TeV 2025, period LHC25ae: CCDB path Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25ae
      ///             Run list: 564356 <= runNumber && runNumber <= 564445
      ///
      pathInputFile = "Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25ae";
      LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]:   >>>   OO, 5.36 TeV 2025, period LHC25ae: CCDB path " << pathInputFile;
      LOG(info) << "                                                   >>>   Run list: 564356 <= runNumber && runNumber <= 564445";
    } else if (564468 <= runNumber && runNumber <= 564472) {
      ///
      ///   [CASE 5]: OO, 5.36 TeV 2025, period LHC25af: CCDB path Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25af
      ///             Run list: 564468 <= runNumber && runNumber <= 564472
      ///
      pathInputFile = "Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25af";
      LOG(info) << "[TrackTuner::getPathInputFileAutomaticFromCCDB]:   >>>   OO, 5.36 TeV 2025, period LHC25af: CCDB path " << pathInputFile;
      LOG(info) << "                                                   >>>   Run list: 564468 <= runNumber && runNumber <= 564472";
    } else {
      LOG(fatal) << "runNumber " << runNumber << " not supported for the autodetection. Please switch to manual configuration of the TrackTuner object. Aborting...";
    }

    outputString += ", pathInputFile=" + pathInputFile;
  }

  /// @brief Function to configure the TrackTuner parameters with an input string
  /// @param inputString Input string with all parameter configuration. Format: <name>=<value>|<name>=<value>
  /// @return String with the values of all parameters after configurations are listed, to cross check that everything worked well
  std::string configParams(std::string inputString)
  {

    LOG(info) << "[TrackTuner] /*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/";
    LOG(info) << "[TrackTuner] /*/*/                                             /*/*/";
    LOG(info) << "[TrackTuner] /*/*/   Configuring the TrackTuner via a string   /*/*/";
    LOG(info) << "[TrackTuner] /*/*/                                             /*/*/";
    LOG(info) << "[TrackTuner] /*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/";

    std::string delimiter = "|";
    std::string assignmentSymbol = "=";

    LOG(info) << "[TrackTuner] === ";
    LOG(info) << "[TrackTuner] === Parameter configuration via std::string";
    LOG(info) << "[TrackTuner] === Required format: \"<name>" << assignmentSymbol << "<value>" << delimiter << "<name>" << assignmentSymbol << "<value>" << delimiter << "<name>" << assignmentSymbol << "<value>\"";
    LOG(info) << "[TrackTuner] === Delimiter symbol: \"" << delimiter << "\"";
    LOG(info) << "[TrackTuner] === Assignment symbol: \"" << assignmentSymbol << "\"";
    LOG(info) << "[TrackTuner] === ";
    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> Original input string = \"" << inputString << "\"";

    /// Check the format of the input string
    if (inputString.find(delimiter) == std::string::npos) {
      // wrong delimiter symbol used
      LOG(fatal) << "delimiter symbol \"" << delimiter << "\" not found in the configuration string. Fix it!";
    }
    if (inputString.find(assignmentSymbol) == std::string::npos) {
      // wrong assignment symbol used
      LOG(fatal) << "assignment symbol \"" << assignmentSymbol << "\" not found in the configuration string. Fix it!";
    }
    int spaces = std::count(inputString.begin(), inputString.end(), ' ');
    if (spaces > 0) {
      // white spaces to be removed
      LOG(fatal) << spaces << " white spaces found in the configuration string. Remove them!";
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// Parameters to be configured via the string
    /// +++ to be manually updated every time one adds a new parameter to the TrackTuner.h +++
    enum Pars : uint8_t { DebugInfo = 0,
                          UpdateTrackCovMat,
                          UpdateTrackDCAs,
                          UpdateCurvature,
                          UpdateCurvatureIU,
                          UpdatePulls,
                          IsInputFileFromCCDB,
                          PathInputFile,
                          NameInputFile,
                          PathFileQoverPt,
                          NameFileQoverPt,
                          UsePvRefitCorrections,
                          QOverPtMC,
                          QOverPtData,
                          NPhiBins,
                          AutoDetectDcaCalib,
                          NPars };
    std::map<uint8_t, std::string> mapParNames = {
      std::make_pair(DebugInfo, "debugInfo"),
      std::make_pair(UpdateTrackDCAs, "updateTrackDCAs"),
      std::make_pair(UpdateTrackCovMat, "updateTrackCovMat"),
      std::make_pair(UpdateCurvature, "updateCurvature"),
      std::make_pair(UpdateCurvatureIU, "updateCurvatureIU"),
      std::make_pair(UpdatePulls, "updatePulls"),
      std::make_pair(IsInputFileFromCCDB, "isInputFileFromCCDB"),
      std::make_pair(PathInputFile, "pathInputFile"),
      std::make_pair(PathFileQoverPt, "pathFileQoverPt"),
      std::make_pair(NameInputFile, "nameInputFile"),
      std::make_pair(NameFileQoverPt, "nameFileQoverPt"),
      std::make_pair(UsePvRefitCorrections, "usePvRefitCorrections"),
      std::make_pair(QOverPtMC, "qOverPtMC"),
      std::make_pair(QOverPtData, "qOverPtData"),
      std::make_pair(NPhiBins, "nPhiBins"),
      std::make_pair(AutoDetectDcaCalib, "autoDetectDcaCalib")};
    ///////////////////////////////////////////////////////////////////////////////////
    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> Parameters before the custom settings";
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    LOG(info) << "[TrackTuner]     updateTrackDCAs = " << updateTrackDCAs;
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    LOG(info) << "[TrackTuner]     updateCurvatureIU = " << updateCurvatureIU;
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    LOG(info) << "[TrackTuner]     isInputFileFromCCDB = " << isInputFileFromCCDB;
    LOG(info) << "[TrackTuner]     pathInputFile = " << pathInputFile;
    LOG(info) << "[TrackTuner]     nameInputFile = " << nameInputFile;
    LOG(info) << "[TrackTuner]     pathFileQoverPt = " << pathFileQoverPt;
    LOG(info) << "[TrackTuner]     nameFileQoverPt = " << nameFileQoverPt;
    LOG(info) << "[TrackTuner]     usePvRefitCorrections = " << usePvRefitCorrections;
    LOG(info) << "[TrackTuner]     qOverPtMC = " << qOverPtMC;
    LOG(info) << "[TrackTuner]     qOverPtData = " << qOverPtData;
    LOG(info) << "[TrackTuner]     nPhiBins = " << nPhiBins;
    LOG(info) << "[TrackTuner]     autoDetectDcaCalib = " << autoDetectDcaCalib;

    // ##############################################################################################
    // ########   split the original string, separating substrings delimited by "|" symbol   ########

    std::vector<std::string> slices = {};

    while (inputString.find(delimiter) != std::string::npos) {
      /// we are not at the end of our string --> let's find out the next parameter to be configured!
      slices.push_back(inputString.substr(0, inputString.find(delimiter))); // this gives us the substring until the next "|" character
      inputString.erase(0, slices.back().length() + delimiter.length());    // erase the found substring for next iteration
    }
    /// at this point, the input string is erased until the last delimiter (included)
    slices.push_back(inputString); // necessary to read the last parameter, after the last "|" symbol

    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> String slices:";
    for (const std::string& s : slices)
      LOG(info) << "[TrackTuner]     " << s;

    /// check if the number of input parameters is correct
    if (static_cast<uint8_t>(slices.size()) != NPars) {
      LOG(fatal) << "[TrackTuner] " << slices.size() << " parameters provided, while " << static_cast<int>(NPars) << " are expected. Fix it!";
    }

    // ###################################################################################################################
    // ########   each split is now a std::string "<parName>=<value>" --> let's really configure each parameter   ########

    /// lambda expression to search for the parameter value (as string) in the configuration string
    auto getValueString = [&](uint8_t iPar) {
      /// this allows to search the parameter configuration even if they are not written in order
      auto it = std::find_if(slices.begin(), slices.end(), [&](const std::string& s) { return s.find(mapParNames[iPar]) != std::string::npos; });
      if (it == std::end(slices)) {
        // parameter not found
        LOG(fatal) << "\"" << mapParNames[iPar] << "\" not found in the configuration string";
      }
      std::string str = *it;
      if (str.find('=') == std::string::npos || str.back() == '=') {
        // value of the parameter missing in the configuration string
        LOG(fatal) << "Missing value for \"" << mapParNames[iPar] << "\" in the configuration string";
      }
      return str.substr(str.find(assignmentSymbol) + 1, str.length());
    };

    /// further lambda expression to handle bool initialization
    auto setBoolFromString = [=](bool& b, const std::string& str) {
      if (!str.compare("1") || str.find("true") != std::string::npos || str.find("True") != std::string::npos || str.find("TRUE") != std::string::npos) {
        b = true;
      } else if (!str.compare("0") || str.find("false") != std::string::npos || str.find("False") != std::string::npos || str.find("FALSE") != std::string::npos) {
        b = false;
      } else {
        LOG(fatal) << "[TrackTuner] Wrong bool initialization from configuration ";
      }
    };

    LOG(info) << "[TrackTuner] ";
    LOG(info) << "[TrackTuner] >>> Parameters after the custom settings";
    // Configure autoDetectDcaCalib
    setBoolFromString(autoDetectDcaCalib, getValueString(AutoDetectDcaCalib));
    outputString += "autoDetectDcaCalib=" + std::to_string(autoDetectDcaCalib);
    LOG(info) << "[TrackTuner]     autoDetectDcaCalib = " << autoDetectDcaCalib;
    // Configure debugInfo
    setBoolFromString(debugInfo, getValueString(DebugInfo));
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    outputString += ", debugInfo=" + std::to_string(debugInfo);
    // Configure updateTrackDCAs
    setBoolFromString(updateTrackDCAs, getValueString(UpdateTrackDCAs));
    LOG(info) << "[TrackTuner]     updateTrackDCAs = " << updateTrackDCAs;
    outputString += ", updateTrackDCAs=" + std::to_string(updateTrackDCAs);
    // Configure updateTrackCovMat
    setBoolFromString(updateTrackCovMat, getValueString(UpdateTrackCovMat));
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    outputString += ", updateTrackCovMat=" + std::to_string(updateTrackCovMat);
    // Configure updateCurvature
    setBoolFromString(updateCurvature, getValueString(UpdateCurvature));
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    outputString += ", updateCurvature=" + std::to_string(updateCurvature);
    // Configure updateCurvatureIU
    setBoolFromString(updateCurvatureIU, getValueString(UpdateCurvatureIU));
    LOG(info) << "[TrackTuner]     updateCurvatureIU = " << updateCurvatureIU;
    outputString += ", updateCurvatureIU=" + std::to_string(updateCurvatureIU);
    // Configure updatePulls
    setBoolFromString(updatePulls, getValueString(UpdatePulls));
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    outputString += ", updatePulls=" + std::to_string(updatePulls);
    // Configure isInputFileFromCCDB
    setBoolFromString(isInputFileFromCCDB, getValueString(IsInputFileFromCCDB));
    LOG(info) << "[TrackTuner]     isInputFileFromCCDB = " << isInputFileFromCCDB;
    outputString += ", isInputFileFromCCDB=" + std::to_string(isInputFileFromCCDB);
    // Configure pathInputFile
    pathFileQoverPt = getValueString(PathFileQoverPt);
    outputString += ", pathFileQoverPt=" + pathFileQoverPt;
    LOG(info) << "[TrackTuner]     pathFileQoverPt = " << pathFileQoverPt;
    // Configure nameInputFile
    nameInputFile = getValueString(NameInputFile);
    outputString += ", nameInputFile=" + nameInputFile;
    LOG(info) << "[TrackTuner]     nameInputFile = " << nameInputFile;
    // Configure nameFileQoverPt
    nameFileQoverPt = getValueString(NameFileQoverPt);
    outputString += ", nameFileQoverPt=" + nameFileQoverPt;
    LOG(info) << "[TrackTuner]     nameFileQoverPt = " << nameFileQoverPt;
    // Configure usePvRefitCorrections
    setBoolFromString(usePvRefitCorrections, getValueString(UsePvRefitCorrections));
    outputString += ", usePvRefitCorrections=" + std::to_string(usePvRefitCorrections);
    LOG(info) << "[TrackTuner]     usePvRefitCorrections = " << usePvRefitCorrections;
    // Configure qOverPtMC
    qOverPtMC = std::stof(getValueString(QOverPtMC));
    outputString += ", qOverPtMC=" + std::to_string(qOverPtMC);
    LOG(info) << "[TrackTuner]     qOverPtMC = " << qOverPtMC;
    // Configure qOverPtData
    qOverPtData = std::stof(getValueString(QOverPtData));
    outputString += ", qOverPtData=" + std::to_string(qOverPtData);
    LOG(info) << "[TrackTuner]     qOverPtData = " << qOverPtData;
    // Configure nPhiBins
    nPhiBins = std::stoi(getValueString(NPhiBins));
    outputString += ", nPhiBins=" + std::to_string(nPhiBins);
    if (nPhiBins < 0)
      LOG(fatal) << "[TrackTuner]   negative nPhiBins!" << nPhiBins;
    LOG(info) << "[TrackTuner]     nPhiBins = " << nPhiBins;
    // Configure pathInputFile
    if (!autoDetectDcaCalib) {
      // path input file from the input string
      pathInputFile = getValueString(PathInputFile);
      outputString += ", pathInputFile=" + pathInputFile;
      LOG(info) << "[TrackTuner]     pathInputFile = " << pathInputFile;
    } else {
      LOG(info) << "[TrackTuner]     pathInputFile still invalid for the moment --> it will be updated by the \"auto-detect\"";
    }
    /// declare that the configuration is done via an input string
    isConfigFromString = true;

    /// sanity-checks on the configurations
    checkConfig();

    return outputString;
  }

  /// @brief Function to configure the TrackTuner parameters with an input string
  /// @return String with the values of all parameters after configurations are listed, to cross check that everything worked well
  std::string configParams()
  {

    LOG(info) << "[TrackTuner] /=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#";
    LOG(info) << "[TrackTuner] /=/#/                                                               /=/#/";
    LOG(info) << "[TrackTuner] /=/#/   Configuring the TrackTuner using the input Configurables    /=/#/";
    LOG(info) << "[TrackTuner] /=/#/                                                               /=/#/";
    LOG(info) << "[TrackTuner] /=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/=/#/";

    LOG(info) << "[TrackTuner] ";
    LOG(info) << "[TrackTuner] >>> Parameters after the custom settings";
    // Configure autoDetectDcaCalib
    autoDetectDcaCalib = cfgAutoDetectDcaCalib;
    outputString += "autoDetectDcaCalib=" + std::to_string(autoDetectDcaCalib);
    LOG(info) << "[TrackTuner]     autoDetectDcaCalib = " << autoDetectDcaCalib;
    // Configure debugInfo
    debugInfo = cfgDebugInfo;
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    outputString += ", debugInfo=" + std::to_string(debugInfo);
    // Configure updateTrackDCAs
    updateTrackDCAs = cfgUpdateTrackDCAs;
    LOG(info) << "[TrackTuner]     updateTrackDCAs = " << updateTrackDCAs;
    outputString += ", updateTrackDCAs=" + std::to_string(updateTrackDCAs);
    // Configure updateTrackCovMat
    updateTrackCovMat = cfgUpdateTrackCovMat;
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    outputString += ", updateTrackCovMat=" + std::to_string(updateTrackCovMat);
    // Configure updateCurvature
    updateCurvature = cfgUpdateCurvature;
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    outputString += ", updateCurvature=" + std::to_string(updateCurvature);
    // Configure updateCurvatureIU
    updateCurvatureIU = cfgUpdateCurvatureIU;
    LOG(info) << "[TrackTuner]     updateCurvatureIU = " << updateCurvatureIU;
    outputString += ", updateCurvatureIU=" + std::to_string(updateCurvatureIU);
    // Configure updatePulls
    updatePulls = cfgUpdatePulls;
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    outputString += ", updatePulls=" + std::to_string(updatePulls);
    // Configure isInputFileFromCCDB
    isInputFileFromCCDB = cfgIsInputFileFromCCDB;
    LOG(info) << "[TrackTuner]     isInputFileFromCCDB = " << isInputFileFromCCDB;
    outputString += ", isInputFileFromCCDB=" + std::to_string(isInputFileFromCCDB);
    // Configure pathQoverPt
    pathFileQoverPt = cfgPathFileQoverPt;
    outputString += ", pathFileQoverPt=" + pathFileQoverPt;
    LOG(info) << "[TrackTuner]     pathFileQoverPt = " << pathFileQoverPt;
    // Configure nameInputFile
    nameInputFile = cfgNameInputFile;
    outputString += ", nameInputFile=" + nameInputFile;
    LOG(info) << "[TrackTuner]     nameInputFile = " << nameInputFile;
    // Configure nameFileQoverPt
    nameFileQoverPt = cfgNameFileQoverPt;
    outputString += ", nameFileQoverPt=" + nameFileQoverPt;
    LOG(info) << "[TrackTuner]     nameFileQoverPt = " << nameFileQoverPt;
    // Configure usePvRefitCorrections
    usePvRefitCorrections = cfgUsePvRefitCorrections;
    outputString += ", usePvRefitCorrections=" + std::to_string(usePvRefitCorrections);
    LOG(info) << "[TrackTuner]     usePvRefitCorrections = " << usePvRefitCorrections;
    // Configure qOverPtMC
    qOverPtMC = cfgQOverPtMC;
    outputString += ", qOverPtMC=" + std::to_string(qOverPtMC);
    LOG(info) << "[TrackTuner]     qOverPtMC = " << qOverPtMC;
    // Configure qOverPtData
    qOverPtData = cfgQOverPtData;
    outputString += ", qOverPtData=" + std::to_string(qOverPtData);
    LOG(info) << "[TrackTuner]     qOverPtData = " << qOverPtData;
    // Configure nPhiBins
    nPhiBins = cfgNPhiBins;
    outputString += ", nPhiBins=" + std::to_string(nPhiBins);
    if (nPhiBins < 0)
      LOG(fatal) << "[TrackTuner]   negative nPhiBins!" << nPhiBins;
    LOG(info) << "[TrackTuner]     nPhiBins = " << nPhiBins;
    // Configure pathInputFile
    if (!autoDetectDcaCalib) {
      // path input file from configurable
      pathInputFile = cfgPathInputFile;
      outputString += ", pathInputFile=" + pathInputFile;
      LOG(info) << "[TrackTuner]     pathInputFile = " << pathInputFile;
    } else {
      LOG(info) << "[TrackTuner]     pathInputFile still invalid for the moment --> it will be updated by the \"auto-detect\"";
    }

    /// declare that the configuration is done via the Configurables
    isConfigFromConfigurables = true;

    /// sanity-checks on the configurations
    checkConfig();

    return outputString;
  }

  void getDcaGraphs()
  {
    /// abort if the graphs were already loaded
    if (areGraphsConfigured) {
      LOG(fatal) << "[TrackTuner::getDcaGraphs()] Function already called, i.e. the calibrations are already loaded. This further call should never happen. Aborting...";
    }

    std::string fullNameInputFile = pathInputFile + std::string("/") + nameInputFile;
    std::string fullNameFileQoverPt = pathFileQoverPt + std::string("/") + nameFileQoverPt;
    TList* ccdb_object_dca = nullptr;
    TList* ccdb_object_qoverpt = nullptr;

    std::string grOneOverPtPionNameMC = "sigmaVsPtMc";
    std::string grOneOverPtPionNameData = "sigmaVsPtData";

    if (isInputFileFromCCDB) {
      /// use input correction file from CCDB

      // get the TList from the DCA correction file present in CCDB
      ccdb_object_dca = o2::ccdb::BasicCCDBManager::instance().get<TList>(pathInputFile);
      LOG(info) << "[TrackTuner] ccdb_object_dca " << ccdb_object_dca;

      // get the TList from the Q/Pt correction file from CCDB
      if (updateCurvature || updateCurvatureIU) {
        ccdb_object_qoverpt = o2::ccdb::BasicCCDBManager::instance().get<TList>(pathFileQoverPt);
        LOG(info) << "[TrackTuner] ccdb_object_qoverpt " << ccdb_object_qoverpt;
      }
    } else {
      /// use input correction file from local filesystem

      /// open the input correction file - dca correction
      TFile* inputFile = TFile::Open(fullNameInputFile.c_str(), "READ");
      if (!inputFile) {
        LOG(fatal) << "[TrackTuner] Something wrong with the local input file" << fullNameInputFile << " for dca correction. Fix it!";
      }
      ccdb_object_dca = dynamic_cast<TList*>(inputFile->Get("ccdb_object"));

      /// open the input correction file - q/pt correction
      TFile* inputFileQoverPt = TFile::Open(fullNameFileQoverPt.c_str(), "READ");
      if (!inputFileQoverPt && (updateCurvature || updateCurvatureIU)) {
        LOG(fatal) << "Something wrong with the Q/Pt input file" << fullNameFileQoverPt << " for Q/Pt correction. Fix it!";
      }
      ccdb_object_qoverpt = dynamic_cast<TList*>(inputFileQoverPt->Get("ccdb_object"));
    }

    // choose wheter to use corrections w/ PV refit or w/o it, and retrieve the proper TList
    std::string dir = "woPvRefit";
    if (usePvRefitCorrections) {
      dir = "withPvRefit";
    }
    TList* td = dynamic_cast<TList*>(ccdb_object_dca->FindObject(dir.c_str()));
    if (!td) {
      LOG(fatal) << "[TrackTuner] TList " << td << " not found in ccdb_object_dca. Fix it!";
    }

    int inputNphiBins = nPhiBins;
    if (inputNphiBins == 0)
      nPhiBins = 1; // old phi_independent settings

    // reserve  memory and initialize vector for needed number of graphs
    grDcaXYResVsPtPionMC.resize(nPhiBins);
    grDcaXYResVsPtPionData.resize(nPhiBins);

    grDcaZResVsPtPionMC.resize(nPhiBins);
    grDcaZResVsPtPionData.resize(nPhiBins);

    grDcaXYMeanVsPtPionMC.resize(nPhiBins);
    grDcaXYMeanVsPtPionData.resize(nPhiBins);

    grDcaZMeanVsPtPionMC.resize(nPhiBins);
    grDcaZMeanVsPtPionData.resize(nPhiBins);

    grDcaXYPullVsPtPionMC.resize(nPhiBins);
    grDcaXYPullVsPtPionData.resize(nPhiBins);

    grDcaZPullVsPtPionMC.resize(nPhiBins);
    grDcaZPullVsPtPionData.resize(nPhiBins);

    /// Lambda expression to get the TGraphErrors from file
    auto loadGraph = [&](int phiBin, const std::string& strBaseName) -> TGraphErrors* {
      std::string strGraphName = inputNphiBins != 0 ? fmt::format("{}_{}", strBaseName, phiBin) : strBaseName;
      TObject* obj = td->FindObject(strGraphName.c_str());
      if (!obj) {
        LOG(fatal) << "[TrackTuner]     TGraphErrors not found in the Input Root file: " << strGraphName;
        td->ls();
        return nullptr;
      }
      return dynamic_cast<TGraphErrors*>(obj);
    };

    if (inputNphiBins != 0) {
      LOG(info) << "[TrackTuner]    Loading phi-dependent XY TGraphErrors";
    }
    for (int iPhiBin = 0; iPhiBin < nPhiBins; ++iPhiBin) {

      grDcaXYResVsPtPionMC[iPhiBin].reset(loadGraph(iPhiBin, "resCurrentDcaXY"));
      grDcaXYResVsPtPionData[iPhiBin].reset(loadGraph(iPhiBin, "resUpgrDcaXY"));
      grDcaXYMeanVsPtPionMC[iPhiBin].reset(loadGraph(iPhiBin, "meanCurrentDcaXY"));
      grDcaXYMeanVsPtPionData[iPhiBin].reset(loadGraph(iPhiBin, "meanUpgrDcaXY"));
      grDcaXYPullVsPtPionMC[iPhiBin].reset(loadGraph(iPhiBin, "pullsCurrentDcaXY"));
      grDcaXYPullVsPtPionData[iPhiBin].reset(loadGraph(iPhiBin, "pullsUpgrDcaXY"));

      if (!grDcaXYResVsPtPionMC[iPhiBin].get() || !grDcaXYResVsPtPionData[iPhiBin].get() || !grDcaXYMeanVsPtPionMC[iPhiBin].get() || !grDcaXYMeanVsPtPionData[iPhiBin].get() || !grDcaXYPullVsPtPionMC[iPhiBin].get() || !grDcaXYPullVsPtPionData[iPhiBin].get()) {
        LOG(fatal) << "[TrackTuner]     Something wrong with the names of the correction graphs for dcaXY. Fix it! Problematic phi bin is" << iPhiBin;
      }
    }

    if (inputNphiBins != 0) {
      LOG(info) << "[TrackTuner]    Loading phi-dependent Z TGraphErrors";
    }
    for (int iPhiBin = 0; iPhiBin < nPhiBins; ++iPhiBin) {
      grDcaZResVsPtPionMC[iPhiBin].reset(loadGraph(iPhiBin, "resCurrentDcaZ"));
      grDcaZMeanVsPtPionMC[iPhiBin].reset(loadGraph(iPhiBin, "meanCurrentDcaZ"));
      grDcaZPullVsPtPionMC[iPhiBin].reset(loadGraph(iPhiBin, "pullsCurrentDcaZ"));
      grDcaZResVsPtPionData[iPhiBin].reset(loadGraph(iPhiBin, "resUpgrDcaZ"));
      grDcaZMeanVsPtPionData[iPhiBin].reset(loadGraph(iPhiBin, "meanUpgrDcaZ"));
      grDcaZPullVsPtPionData[iPhiBin].reset(loadGraph(iPhiBin, "pullsUpgrDcaZ"));

      if (!grDcaZResVsPtPionMC[iPhiBin].get() || !grDcaZResVsPtPionData[iPhiBin].get() || !grDcaZMeanVsPtPionMC[iPhiBin].get() || !grDcaZMeanVsPtPionData[iPhiBin].get() || !grDcaZPullVsPtPionMC[iPhiBin].get() || !grDcaZPullVsPtPionData[iPhiBin].get()) {
        LOG(fatal) << "[TrackTuner] Something wrong with the names of the correction graphs for dcaZ. Fix it! Problematic phi bin is" << iPhiBin;
      }
    }

    if (updateCurvature || updateCurvatureIU) {
      grOneOverPtPionMC.reset(dynamic_cast<TGraphErrors*>(ccdb_object_qoverpt->FindObject(grOneOverPtPionNameMC.c_str())));
      grOneOverPtPionData.reset(dynamic_cast<TGraphErrors*>(ccdb_object_qoverpt->FindObject(grOneOverPtPionNameData.c_str())));
    }

    /// if we arrive here, it means that the graphs are all set
    areGraphsConfigured = true;

  } // getDcaGraphs() ends here

  template <typename T1, typename T2, typename T3, typename T4, typename H>
  void tuneTrackParams(T1 const& mcparticle, T2& trackParCov, T3 const& matCorr, T4 dcaInfoCov, H hQA)
  {
    /// abort if the calibrations are not loaded
    if (!areGraphsConfigured) {
      LOG(fatal) << "[TrackTuner::tuneTrackParams()] Function called, but calibration graphs not configured. Have you called the function TrackTuner::getDcaGraphs()? Aborting...";
    }

    double ptMC = mcparticle.pt();
    double dcaXYResMC = 0.0; // sd0rpo=0.;
    double dcaZResMC = 0.0;  // sd0zo =0.;

    double dcaXYResData = 0.0; // sd0rpn=0.;
    double dcaZResData = 0.0;  // sd0zn =0.;

    double dcaXYMeanMC = 0.0;   // sd0mrpo=0.;
    double dcaXYMeanData = 0.0; // sd0mrpn=0.;

    double dcaXYPullMC = 1.0;
    double dcaXYPullData = 1.0;

    double dcaZPullMC = 1.0;
    double dcaZPullData = 1.0;

    // get phibin
    double phiMC = mcparticle.phi();
    if (phiMC < 0.)
      phiMC += o2::constants::math::TwoPI;                                    // 2 * std::numbers::pi;//
    int phiBin = phiMC / (o2::constants::math::TwoPI + 0.0000001) * nPhiBins; // 0.0000001 just a numerical protection

    dcaXYResMC = evalGraph(ptMC, grDcaXYResVsPtPionMC[phiBin].get());
    dcaXYResData = evalGraph(ptMC, grDcaXYResVsPtPionData[phiBin].get());

    dcaZResMC = evalGraph(ptMC, grDcaZResVsPtPionMC[phiBin].get());
    dcaZResData = evalGraph(ptMC, grDcaZResVsPtPionData[phiBin].get());

    // For Q/Pt corrections, files on CCDB will be used if both qOverPtMC and qOverPtData are null
    if (updateCurvature || updateCurvatureIU) {
      if ((qOverPtMC < 0) || (qOverPtData < 0)) {
        if (debugInfo) {
          LOG(info) << "### q/pt smearing: qOverPtMC=" << qOverPtMC << ", qOverPtData=" << qOverPtData << ". One of them is negative. Retrieving then values from graphs from input .root file";
        }
        /// check that input graphs for q/pt smearing are correctly retrieved
        if (!grOneOverPtPionData.get() || !grOneOverPtPionMC.get()) {
          LOG(fatal) << "### q/pt smearing: input graphs not correctly retrieved. Aborting.";
        }
        qOverPtMC = std::max(0.0, evalGraph(ptMC, grOneOverPtPionMC.get()));
        qOverPtData = std::max(0.0, evalGraph(ptMC, grOneOverPtPionData.get()));
      } // qOverPtMC, qOverPtData block ends here
    } // updateCurvature, updateCurvatureIU block ends here

    if (updateTrackDCAs) {

      dcaXYMeanMC = evalGraph(ptMC, grDcaXYMeanVsPtPionMC[phiBin].get());
      dcaXYMeanData = evalGraph(ptMC, grDcaXYMeanVsPtPionData[phiBin].get());

      dcaXYPullMC = evalGraph(ptMC, grDcaXYPullVsPtPionMC[phiBin].get());
      dcaXYPullData = evalGraph(ptMC, grDcaXYPullVsPtPionData[phiBin].get());

      dcaZPullMC = evalGraph(ptMC, grDcaZPullVsPtPionMC[phiBin].get());
      dcaZPullData = evalGraph(ptMC, grDcaZPullVsPtPionData[phiBin].get());
    }
    //  Unit conversion, is it required ??
    dcaXYResMC *= 1.e-4;
    dcaZResMC *= 1.e-4;

    dcaXYResData *= 1.e-4;
    dcaZResData *= 1.e-4;

    dcaXYMeanMC *= 1.e-4;
    dcaXYMeanData *= 1.e-4;

    // Apply the smearing
    // ---------------------------------------------
    // double pt1o  =param  [4];
    double trackParQPtMCRec = trackParCov.getQ2Pt();
    int sign = trackParCov.getQ2Pt() / std::abs(trackParCov.getQ2Pt());
    // double pt1mc =parammc[4];
    double trackParQPtMC = sign / mcparticle.pt();
    o2::dataformats::VertexBase vtxMC;
    vtxMC.setPos({mcparticle.vx(), mcparticle.vy(), mcparticle.vz()});
    vtxMC.setCov(0, 0, 0, 0, 0, 0); // ??? or All ZEROs // == 1 cm2? wrt prop point

    if (debugInfo) {
      LOG(info) << mcparticle.pt() << " " << 1 / trackParCov.getPtInv() << " " << (trackParQPtMCRec - trackParQPtMC) / trackParQPtMC * 100.0;
    }

    // for updating the Q/pT from the tracksIU.
    // This is placed before track propagation to production point, so that Q/pT can be updated before propagation
    // Q/Pt is modified and set before track propagation

    double deltaQpt = 0.0;
    double deltaQptTuned = 0.0;
    double trackParQPtTuned = 0.0;

    // variables for track cov matrix elements update
    double sigmaY2 = 0.0;
    double sigmaZY = 0.0;
    double sigmaZ2 = 0.0;
    double sigmaSnpY = 0.0;
    double sigmaSnpZ = 0.0;
    double sigmaTglY = 0.0;
    double sigmaTglZ = 0.0;
    double sigma1PtY = 0.0;
    double sigma1PtZ = 0.0;
    double sigma1PtSnp = 0.0;
    double sigma1PtTgl = 0.0;
    double sigma1Pt2 = 0.0;

    double sigmaY2orig = trackParCov.getSigmaY2();
    double sigmaZYorig = trackParCov.getSigmaZY();
    double sigmaZ2orig = trackParCov.getSigmaZ2();

    if (updateCurvatureIU) {
      // double dpt1o =pt1o-pt1mc;
      deltaQpt = trackParQPtMCRec - trackParQPtMC;
      // double dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
      deltaQptTuned = deltaQpt * (qOverPtMC > 0. ? (qOverPtData / qOverPtMC) : 1.);
      // double pt1n  = pt1mc+dpt1n;
      trackParQPtTuned = trackParQPtMC + deltaQptTuned;
      trackParCov.setQ2Pt(trackParQPtTuned);

      // updating track cov matrix elements for 1/Pt at innermost update point
      //       if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
      sigma1PtY = trackParCov.getSigma1PtY();
      if (dcaXYResMC > 0. && qOverPtMC > 0.) {
        sigma1PtY *= ((dcaXYResData / dcaXYResMC) * (qOverPtData / qOverPtMC));
        trackParCov.setCov(sigma1PtY, 10);
      }

      //       if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
      sigma1PtZ = trackParCov.getSigma1PtZ();
      if (dcaZResMC > 0. && qOverPtMC > 0.) {
        sigma1PtZ *= ((dcaZResData / dcaZResMC) * (qOverPtData / qOverPtMC));
        trackParCov.setCov(sigma1PtZ, 11);
      }

      //       if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
      sigma1PtSnp = trackParCov.getSigma1PtSnp();
      if (qOverPtMC > 0.) {
        sigma1PtSnp *= (qOverPtData / qOverPtMC);
        trackParCov.setCov(sigma1PtSnp, 12);
      }

      //       if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
      sigma1PtTgl = trackParCov.getSigma1PtTgl();
      if (qOverPtMC > 0.) {
        sigma1PtTgl *= (qOverPtData / qOverPtMC);
        trackParCov.setCov(sigma1PtTgl, 13);
      }

      //       if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
      sigma1Pt2 = trackParCov.getSigma1Pt2();
      if (qOverPtMC > 0.) {
        sigma1Pt2 *= (qOverPtData / qOverPtMC) * (qOverPtData / qOverPtMC);
        trackParCov.setCov(sigma1Pt2, 14);
      }
    } // updateCurvatureIU block ends here
    // propagate to DCA with respect to the Production point
    // if (!updateCurvatureIU) {
    //   o2::base::Propagator::Instance()->propagateToDCABxByBz(vtxMC, trackParCov, 2.f, matCorr, dcaInfoCov);
    // }

    double trackParDcaXYoriginal = trackParCov.getY();
    double trackParDcaZoriginal = trackParCov.getZ();

    if (updateTrackDCAs) {
      // propagate to DCA with respect to the Production point
      o2::base::Propagator::Instance()->propagateToDCABxByBz(vtxMC, trackParCov, 2.f, matCorr, dcaInfoCov);
      if (debugInfo) {
        LOG(info) << "phi MC" << mcparticle.phi();
        LOG(info) << "alpha track" << trackParCov.getAlpha();
      }
      // double d0zo  =param  [1];
      double trackParDcaZRec = trackParCov.getZ();
      // double d0rpo =param  [0];
      double trackParDcaXYRec = trackParCov.getY();
      float mcVxRotated = mcparticle.vx() * std::cos(trackParCov.getAlpha()) + mcparticle.vy() * std::sin(trackParCov.getAlpha()); // invert
      float mcVyRotated = mcparticle.vy() * std::cos(trackParCov.getAlpha()) - mcparticle.vx() * std::sin(trackParCov.getAlpha());

      if (debugInfo) {
        // LOG(info) << "mcVy " <<  mcparticle.vy()   <<  std::endl;
        LOG(info) << "mcVxRotated " << mcVxRotated;
        LOG(info) << "mcVyRotated " << mcVyRotated;
      }
      // std::array<float, 3>  arrayXYZ = { mcVxRotated, mcVyRotated , mcparticle.vz()};
      // std::array<float, 3>  arrayPxPyPz = {mcparticle.px(), mcparticle.py(), mcparticle.pz()};

      // const int matSize = 21;
      // std::array<float, matSize> arrayCovMatMC = {0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0. };
      // o2::track::TrackParametrizationWithError<float> trackParCovMC{arrayXYZ, arrayPxPyPz,  arrayCovMatMC, trackParCov.getSign()};

      // double d0rpmc=parammc[0];
      // double trackParDcaXYMC = trackParCovMC.getY(); // here

      double trackParDcaXYMC = mcVyRotated; // here

      // double d0zmc =parammc[1];
      double trackParDcaZMC = mcparticle.vz();

      // double dd0zo =d0zo-d0zmc;
      double deltaDcaZ = trackParDcaZRec - trackParDcaZMC;

      // double dd0zn =dd0zo *(sd0zo >0. ? (sd0zn /sd0zo ) : 1.);
      double deltaDcaZTuned = deltaDcaZ * (dcaZResMC > 0. ? (dcaZResData / dcaZResMC) : 1.);

      // double d0zn  =d0zmc+dd0zn;
      double trackParDcaZTuned = trackParDcaZMC + deltaDcaZTuned;

      // double dd0rpo=d0rpo-d0rpmc;
      double deltaDcaXY = trackParDcaXYRec - trackParDcaXYMC - dcaXYMeanMC;

      // double dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
      double deltaDcaXYTuned = deltaDcaXY * (dcaXYResMC > 0. ? (dcaXYResData / dcaXYResMC) : 1.);

      // double dd0mrpn=std::abs(sd0mrpn)-std::abs(sd0mrpo);
      // double deltaDcaXYmean = std::abs(dcaXYMeanData) - std::abs(dcaXYMeanMC) ;
      double deltaDcaXYmean = dcaXYMeanData - dcaXYMeanMC;

      // double d0rpn =d0rpmc+dd0rpn-dd0mrpn;
      double trackParDcaXYTuned = trackParDcaXYMC + deltaDcaXYTuned + deltaDcaXYmean;

      if (debugInfo) {
        LOG(info) << dcaZResMC << ", " << dcaZResData << ", diff(DcaZ - DcaZMC): " << deltaDcaZ << ", diff upgraded: " << deltaDcaZTuned << ", DcaZ Data : " << trackParDcaZTuned;
        LOG(info) << dcaXYResMC << ", " << dcaXYResData << ", " << dcaXYMeanMC << ", diff(DcaY - DcaYMC - dcaXYMeanMC): " << deltaDcaXY << ", diff upgraded: " << deltaDcaXYTuned << ", DcaY Data :" << trackParDcaXYTuned;
      }
      // option mimic data
      // ----------------------
      // if(fMimicData){
      //     // dd0mrpn=sd0mrpn-sd0mrpo;
      //     deltaDcaXYmean = dcaXYMeanData - dcaXYMeanMC;
      //     // d0rpn = d0rpmc+dd0rpn+dd0mrpn;
      //     trackParDcaXYTuned =  deltaDcaXY + deltaDcaXYTuned + deltaDcaXYmean;
      // }

      // setting updated track parameters
      // --------------------------------
      trackParDcaXYoriginal = trackParCov.getY();
      trackParCov.setY(trackParDcaXYTuned);
      trackParDcaZoriginal = trackParCov.getZ();
      trackParCov.setZ(trackParDcaZTuned);
    } // ----> updateTrackDCAs block ends here

    if ((updateCurvature) && (!updateCurvatureIU)) { // ...block begins here
      if (!updateTrackDCAs) {
        /// propagation to production point not done yet, doing it now
        o2::base::Propagator::Instance()->propagateToDCABxByBz(vtxMC, trackParCov, 2.f, matCorr, dcaInfoCov);
      }
      deltaQpt = trackParQPtMCRec - trackParQPtMC;
      // double dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
      deltaQptTuned = deltaQpt * (qOverPtMC > 0. ? (qOverPtData / qOverPtMC) : 1.);
      // double pt1n  = pt1mc+dpt1n;
      trackParQPtTuned = trackParQPtMC + deltaQptTuned;
      trackParCov.setQ2Pt(trackParQPtTuned);
    } // ...block ends here

    if (updateTrackCovMat) {
      //       if(sd0rpo>0.)            covar[0]*=(sd0rpn/sd0rpo)*(sd0rpn/sd0rpo);//yy
      sigmaY2 = trackParCov.getSigmaY2();
      if (dcaXYResMC > 0.) {
        sigmaY2 *= ((dcaXYResData / dcaXYResMC) * (dcaXYResData / dcaXYResMC));
        trackParCov.setCov(sigmaY2, 0);
      }

      //       if(sd0zo>0. && sd0rpo>0.)covar[1]*=(sd0rpn/sd0rpo)*(sd0zn/sd0zo);//yz
      sigmaZY = trackParCov.getSigmaZY();
      if (dcaZResMC > 0. && dcaXYResMC > 0.) {
        // sigmaZY *= ((dcaXYResData / dcaXYResMC) * (dcaZResData / dcaZResMC));
        sigmaZY *= ((dcaXYResData / dcaXYResMC) * (dcaZResData / dcaZResMC));
        trackParCov.setCov(sigmaZY, 1);
      }

      //       if(sd0zo>0.)             covar[2]*=(sd0zn/sd0zo)*(sd0zn/sd0zo);//zz
      sigmaZ2 = trackParCov.getSigmaZ2();
      if (dcaZResMC > 0.) {
        sigmaZ2 *= ((dcaZResData / dcaZResMC) * (dcaZResData / dcaZResMC));
        trackParCov.setCov(sigmaZ2, 2);
      }

      //       if(sd0rpo>0.)            covar[3]*=(sd0rpn/sd0rpo);//yl
      sigmaSnpY = trackParCov.getSigmaSnpY();
      if (dcaXYResMC > 0.) {
        sigmaSnpY *= ((dcaXYResData / dcaXYResMC));
        trackParCov.setCov(sigmaSnpY, 3);
      }

      //       if(sd0zo>0.)             covar[4]*=(sd0zn/sd0zo);//zl
      sigmaSnpZ = trackParCov.getSigmaSnpZ();
      if (dcaZResMC > 0.) {
        sigmaSnpZ *= ((dcaZResData / dcaZResMC));
        trackParCov.setCov(sigmaSnpZ, 4);
      }

      //       if(sd0rpo>0.)            covar[6]*=(sd0rpn/sd0rpo);//ysenT
      sigmaTglY = trackParCov.getSigmaTglY();
      if (dcaXYResMC > 0.) {
        sigmaTglY *= ((dcaXYResData / dcaXYResMC));
        trackParCov.setCov(sigmaTglY, 6);
      }

      //       if(sd0zo>0.)             covar[7]*=(sd0zn/sd0zo);//zsenT
      sigmaTglZ = trackParCov.getSigmaTglZ();
      if (dcaZResMC > 0.) {
        sigmaTglZ *= ((dcaZResData / dcaZResMC));
        trackParCov.setCov(sigmaTglZ, 7);
      }

      // checking and updating track cov matrix elements for 1/Pt begins
      if ((updateCurvature) && (!updateCurvatureIU)) {
        //       if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
        sigma1PtY = trackParCov.getSigma1PtY();
        if (dcaXYResMC > 0. && qOverPtMC > 0.) {
          sigma1PtY *= ((dcaXYResData / dcaXYResMC) * (qOverPtData / qOverPtMC));
          trackParCov.setCov(sigma1PtY, 10);
        }

        //       if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
        sigma1PtZ = trackParCov.getSigma1PtZ();
        if (dcaZResMC > 0. && qOverPtMC > 0.) {
          sigma1PtZ *= ((dcaZResData / dcaZResMC) * (qOverPtData / qOverPtMC));
          trackParCov.setCov(sigma1PtZ, 11);
        }

        //       if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
        sigma1PtSnp = trackParCov.getSigma1PtSnp();
        if (qOverPtMC > 0.) {
          sigma1PtSnp *= (qOverPtData / qOverPtMC);
          trackParCov.setCov(sigma1PtSnp, 12);
        }

        //       if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
        sigma1PtTgl = trackParCov.getSigma1PtTgl();
        if (qOverPtMC > 0.) {
          sigma1PtTgl *= (qOverPtData / qOverPtMC);
          trackParCov.setCov(sigma1PtTgl, 13);
        }

        //       if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
        sigma1Pt2 = trackParCov.getSigma1Pt2();
        if (qOverPtMC > 0.) {
          sigma1Pt2 *= (qOverPtData / qOverPtMC) * (qOverPtData / qOverPtMC);
          trackParCov.setCov(sigma1Pt2, 14);
        }
      } // ---> track cov matrix elements for 1/Pt ends here
    } // ---> updateTrackCovMat block ends here

    if (updatePulls) {
      double ratioDCAxyPulls = 1.0;
      double ratioDCAzPulls = 1.0;
      if (dcaZPullData > 0.0) {
        ratioDCAzPulls = dcaZPullMC / dcaZPullData;
      }
      if (dcaXYPullData > 0.0) {
        ratioDCAxyPulls = dcaXYPullMC / dcaXYPullData;
      }

      // covar[0]*=pullcorr*pullcorr;//yy
      sigmaY2 *= (ratioDCAxyPulls * ratioDCAxyPulls);
      trackParCov.setCov(sigmaY2, 0);

      // covar[1]*=pullcorr;//yz
      sigmaZY *= (ratioDCAxyPulls * ratioDCAzPulls);
      trackParCov.setCov(sigmaZY, 1);

      sigmaZ2 *= (ratioDCAzPulls * ratioDCAzPulls);
      trackParCov.setCov(sigmaZ2, 2);

      // covar[3]*=pullcorr;//yl
      sigmaSnpY *= ratioDCAxyPulls;
      trackParCov.setCov(sigmaSnpY, 3);

      sigmaSnpZ *= ratioDCAzPulls;
      trackParCov.setCov(sigmaSnpZ, 4);

      // covar[6]*=pullcorr;//ysenT
      sigmaTglY *= ratioDCAxyPulls;
      trackParCov.setCov(sigmaTglY, 6);

      sigmaTglZ *= ratioDCAzPulls;
      trackParCov.setCov(sigmaTglZ, 7);

      // covar[10]*=pullcorr;//ypt
      sigma1PtY *= ratioDCAxyPulls;
      trackParCov.setCov(sigma1PtY, 10);

      sigma1PtZ *= ratioDCAzPulls;
      trackParCov.setCov(sigma1PtZ, 11);
    } // ---> updatePulls block ends here

    /// sanity check for track covariance matrix element
    /// see https://github.com/AliceO2Group/AliceO2/blob/66de30958153cd7badf522150e8554f9fcf975ff/Common/DCAFitter/include/DCAFitter/DCAFitterN.h#L38-L54
    double detYZ = sigmaY2 * sigmaZ2 - sigmaZY * sigmaZY;
    double detYZorig = sigmaY2orig * sigmaZ2orig - sigmaZYorig * sigmaZYorig;
    if (detYZ < 0) {
      /// restore the original values for these elements, since in this case the cov. matrix gets ill-defined
      if (debugInfo) {
        LOG(info) << "\n>>> ILL DEFINED COV. MATRIX <<<";
        LOG(info) << "    sigmaY2 = " << sigmaY2;
        LOG(info) << "    sigmaZY = " << sigmaZY;
        LOG(info) << "    sigmaZ2 = " << sigmaZ2;
        LOG(info) << "    The product (sigmaY2*sigmaZ2 - sigmaZY*sigmaZY) = " << detYZ << " is negative. Restoring the original sigmaY2, sigmaZY, sigmaZ2 cov. matrix elements, as well as Y and Z parameters:";
        LOG(info) << "    sigmaY2orig = " << sigmaY2orig;
        LOG(info) << "    sigmaZYorig = " << sigmaZYorig;
        LOG(info) << "    sigmaZ2orig = " << sigmaZ2orig;
        LOG(info) << "    Original product (sigmaY2orig*sigmaZ2orig - sigmaZYorig*sigmaZYorgi) = " << detYZorig;
        LOG(info) << "    ===> track pt = " << trackParCov.getPt();
      }

      // check if this was pathological already w/o track smearing
      if (detYZorig < 0) {
        hQA->Fill(4);
      }

      // restore original Y and Z parameters
      trackParCov.setY(trackParDcaXYoriginal);
      trackParCov.setZ(trackParDcaZoriginal);

      // restore original cov. matrix elements
      trackParCov.setCov(sigmaY2orig, 0);
      trackParCov.setCov(sigmaZYorig, 1);
      trackParCov.setCov(sigmaZ2orig, 2);

      hQA->Fill(3);

    } else {
      hQA->Fill(2);
    }
  } // tuneTrackParams() ends here

  // to be declared
  // ---------------
  // int getPhiBin(double phi) const
  // {
  //   double pi = TMath::Pi();
  //   if (phi > 2. * pi || phi < 0.)
  //     return -1;
  //   if ((phi <= (pi / 4.)) || (phi > 7. * (pi / 4.)))
  //     return 0;
  //   if ((phi > (pi / 4.)) && (phi <= 3. * (pi / 4.)))
  //     return 1;
  //   if ((phi > 3. * (pi / 4.)) && (phi <= 5. * (pi / 4.)))
  //     return 2;
  //   if ((phi > (5. * pi / 4.)) && (phi <= 7. * (pi / 4.)))
  //     return 3;
  //
  //   return -1;
  // }

  double evalGraph(double x, const TGraphErrors* graph) const
  {

    if (!graph) {
      LOG(fatal) << "\t evalGraph fails !\n";
      return 0.;
    }
    int nPoints = graph->GetN();
    double xMin = graph->GetX()[0];
    double xMax = graph->GetX()[nPoints - 1];
    if (x > xMax)
      return graph->Eval(xMax);
    if (x < xMin)
      return graph->Eval(xMin);
    return graph->Eval(x);
  }
};

#endif // COMMON_TOOLS_TRACKTUNER_H_
