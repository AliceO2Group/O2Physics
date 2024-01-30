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

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/GeomConstants.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

#include <TGraphErrors.h>

struct TrackTuner {

  ///////////////////////////////
  /// parameters to be configured
  bool debugInfo = false;
  bool updateTrackCovMat = false;
  bool updateCurvature = false;
  bool updatePulls = false;
  bool isInputFileFromCCDB = false;   // query input file from CCDB or local folder
  std::string pathInputFile = "";     // Path to file containing DCAxy, DCAz graphs from data (upgr) and MC (current)
  std::string nameInputFile = "";     // Common Name of different files containing graphs, found in the above paths
  bool usePvRefitCorrections = false; // establish whether to use corrections obtained with or w/o PV refit
  float oneOverPtCurrent = 0.;        // 1/pt old
  float oneOverPtUpgr = 0.;           // 1/pt new
  ///////////////////////////////

  o2::ccdb::CcdbApi ccdbApi;
  std::map<std::string, std::string> metadata;

  std::unique_ptr<TGraphErrors> grDcaXYResVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaXYResVsPtPionUpgr;

  std::unique_ptr<TGraphErrors> grDcaZResVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaZResVsPtPionUpgr;

  std::unique_ptr<TGraphErrors> grDcaXYMeanVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaXYMeanVsPtPionUpgr;

  std::unique_ptr<TGraphErrors> grDcaZMeanVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaZMeanVsPtPionUpgr;

  std::unique_ptr<TGraphErrors> grOneOverPtPionCurrent;
  std::unique_ptr<TGraphErrors> grOneOverPtPionUpgr;

  std::unique_ptr<TGraphErrors> grDcaXYPullVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaXYPullVsPtPionUpgr;

  std::unique_ptr<TGraphErrors> grDcaZPullVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaZPullVsPtPionUpgr;

  /// @brief Function to configure the TrackTuner parameters
  /// @param inputString Input string with all parameter configuration. Format: <name>=<value>|<name>=<value>
  /// @return String with the values of all parameters after configurations are listed, to cross check that everything worked well
  std::string configParams(std::string inputString)
  {

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
                          UpdateCurvature,
                          UpdatePulls,
                          PathInputFile,
                          IsInputFileFromCCDB,
                          NameInputFile,
                          UsePvRefitCorrections,
                          OneOverPtCurrent,
                          OneOverPtUpgr,
                          NPars };
    std::map<uint8_t, std::string> mapParNames = {
      std::make_pair(DebugInfo, "debugInfo"),
      std::make_pair(UpdateTrackCovMat, "updateTrackCovMat"),
      std::make_pair(UpdateCurvature, "updateCurvature"),
      std::make_pair(UpdatePulls, "updatePulls"),
      std::make_pair(IsInputFileFromCCDB, "isInputFileFromCCDB"),
      std::make_pair(PathInputFile, "pathInputFile"),
      std::make_pair(NameInputFile, "nameInputFile"),
      std::make_pair(UsePvRefitCorrections, "usePvRefitCorrections"),
      std::make_pair(OneOverPtCurrent, "oneOverPtCurrent"),
      std::make_pair(OneOverPtUpgr, "oneOverPtUpgr")};
    ///////////////////////////////////////////////////////////////////////////////////
    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> Parameters before the custom settings";
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    LOG(info) << "[TrackTuner]     isInputFileFromCCDB = " << isInputFileFromCCDB;
    LOG(info) << "[TrackTuner]     pathInputFile = " << pathInputFile;
    LOG(info) << "[TrackTuner]     nameInputFile = " << nameInputFile;
    LOG(info) << "[TrackTuner]     usePvRefitCorrections = " << usePvRefitCorrections;
    LOG(info) << "[TrackTuner]     oneOverPtCurrent = " << oneOverPtCurrent;
    LOG(info) << "[TrackTuner]     oneOverPtUpgr = " << oneOverPtUpgr;

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
    for (std::string& s : slices)
      LOG(info) << "[TrackTuner]     " << s;

    /// check if the number of input parameters is correct
    if (static_cast<uint8_t>(slices.size()) != NPars) {
      LOG(fatal) << "[TrackTuner] " << slices.size() << " parameters provided, while " << NPars << " are expected. Fix it!";
    }

    // ###################################################################################################################
    // ########   each split is now a std::string "<parName>=<value>" --> let's really configure each parameter   ########

    /// lambda expression to search for the parameter value (as string) in the configuration string
    auto getValueString = [&](uint8_t iPar) {
      /// this allows to search the parameter configuration even if they are not written in order
      auto it = std::find_if(slices.begin(), slices.end(), [&](std::string s) { return s.find(mapParNames[iPar]) != std::string::npos; });
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
    auto setBoolFromString = [=](bool& b, std::string str) {
      if (!str.compare("1") || str.find("true") != std::string::npos || str.find("True") != std::string::npos || str.find("TRUE") != std::string::npos) {
        b = true;
      } else if (!str.compare("0") || str.find("false") != std::string::npos || str.find("False") != std::string::npos || str.find("FALSE") != std::string::npos) {
        b = false;
      } else {
        LOG(fatal) << "[TrackTuner] Wrong bool initialization from configuration ";
      }
    };

    std::string outputString = "";
    LOG(info) << "[TrackTuner] ";
    LOG(info) << "[TrackTuner] >>> Parameters after the custom settings";
    // Configure debugInfo
    setBoolFromString(debugInfo, getValueString(DebugInfo));
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    outputString += "debugInfo=" + std::to_string(debugInfo);
    // Configure updateTrackCovMat
    setBoolFromString(updateTrackCovMat, getValueString(UpdateTrackCovMat));
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    outputString += ", updateTrackCovMat=" + std::to_string(updateTrackCovMat);
    // Configure updateCurvature
    setBoolFromString(updateCurvature, getValueString(UpdateCurvature));
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    outputString += ", updateCurvature=" + std::to_string(updateCurvature);
    // Configure updatePulls
    setBoolFromString(updatePulls, getValueString(UpdatePulls));
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    outputString += ", updatePulls=" + std::to_string(updatePulls);
    // Configure isInputFileFromCCDB
    setBoolFromString(isInputFileFromCCDB, getValueString(IsInputFileFromCCDB));
    LOG(info) << "[TrackTuner]     isInputFileFromCCDB = " << isInputFileFromCCDB;
    outputString += ", isInputFileFromCCDB=" + std::to_string(isInputFileFromCCDB);
    // Configure pathInputFile
    pathInputFile = getValueString(PathInputFile);
    outputString += ", pathInputFile=" + pathInputFile;
    LOG(info) << "[TrackTuner]     pathInputFile = " << pathInputFile;
    // Configure nameInputFile
    nameInputFile = getValueString(NameInputFile);
    outputString += ", nameInputFile=" + nameInputFile;
    LOG(info) << "[TrackTuner]     nameInputFile = " << nameInputFile;
    // Configure usePvRefitCorrections
    setBoolFromString(usePvRefitCorrections, getValueString(UsePvRefitCorrections));
    outputString += ", usePvRefitCorrections=" + usePvRefitCorrections;
    LOG(info) << "[TrackTuner]     usePvRefitCorrections = " << usePvRefitCorrections;
    // Configure oneOverPtCurr
    oneOverPtCurrent = std::stof(getValueString(OneOverPtCurrent));
    outputString += ", oneOverPtCurrent=" + std::to_string(oneOverPtCurrent);
    LOG(info) << "[TrackTuner]     oneOverPtCurrent = " << oneOverPtCurrent;
    // Configure oneOverPtUpgr
    oneOverPtUpgr = std::stof(getValueString(OneOverPtUpgr));
    outputString += ", oneOverPtUpgr=" + std::to_string(oneOverPtUpgr);
    LOG(info) << "[TrackTuner]     oneOverPtUpgr = " << oneOverPtUpgr;

    return outputString;
  }

  void getDcaGraphs()
  {
    std::string fullNameInputFile = "";

    if (isInputFileFromCCDB) {
      /// use input correction file from CCDB

      // properly init the ccdb
      std::string tmpDir = ".";
      ccdbApi.init("http://alice-ccdb.cern.ch");

      // get the file from CCDB
      if (!ccdbApi.retrieveBlob(pathInputFile.data(), tmpDir, metadata, 0, false, nameInputFile.data())) {
        LOG(fatal) << "[TrackTuner] input file not found on CCDB, please check the pathInputFile and nameInputFile!";
      }

      // point to the file in the tmp local folder
      fullNameInputFile = tmpDir + std::string("/") + nameInputFile;
    } else {
      /// use input correction file from local filesystem
      fullNameInputFile = pathInputFile + std::string("/") + nameInputFile;
    }

    /// open the input correction file
    std::unique_ptr<TFile> inputFile(TFile::Open(fullNameInputFile.c_str(), "READ"));
    if (!inputFile.get()) {
      LOG(fatal) << "Something wrong with the input file" << fullNameInputFile << " for dca correction. Fix it!";
    }

    // choose wheter to use corrections w/ PV refit or w/o it, and retrieve the proper TDirectory
    std::string dir = "woPvRefit";
    if (usePvRefitCorrections) {
      dir = "withPvRefit";
    }
    TDirectory* td = dynamic_cast<TDirectory*>(inputFile->Get(dir.c_str()));
    if (!td) {
      LOG(fatal) << "TDirectory " << td << " not found in input file" << inputFile->GetName() << ". Fix it!";
    }

    std::string grDcaXYResNameCurr = "resCurrentDcaXY";
    std::string grDcaXYMeanNameCurr = "meanCurrentDcaXY";
    std::string grDcaXYPullNameCurr = "pullsCurrentDcaXY";
    std::string grDcaXYResNameUpgr = "resUpgrDcaXY";
    std::string grDcaXYMeanNameUpgr = "meanUpgrDcaXY";
    std::string grDcaXYPullNameUpgr = "pullsUpgrDcaXY";

    grDcaXYResVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaXYResNameCurr.c_str())));
    grDcaXYResVsPtPionUpgr.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaXYResNameUpgr.c_str())));
    grDcaXYMeanVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaXYMeanNameCurr.c_str())));
    grDcaXYMeanVsPtPionUpgr.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaXYMeanNameUpgr.c_str())));
    grDcaXYPullVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaXYPullNameCurr.c_str())));
    grDcaXYPullVsPtPionUpgr.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaXYPullNameUpgr.c_str())));
    if (!grDcaXYResVsPtPionCurrent.get() || !grDcaXYResVsPtPionUpgr.get() || !grDcaXYMeanVsPtPionCurrent.get() || !grDcaXYMeanVsPtPionUpgr.get() || !grDcaXYPullVsPtPionCurrent.get() || !grDcaXYPullVsPtPionUpgr.get()) {
      LOG(fatal) << "Something wrong with the names of the correction graphs for dcaXY. Fix it!";
    }

    std::string grDcaZResNameCurr = "resCurrentDcaZ";
    std::string grDcaZMeanNameCurr = "meanCurrentDcaZ";
    std::string grDcaZPullNameCurr = "pullsCurrentDcaZ";
    std::string grDcaZResNameUpgr = "resUpgrDcaZ";
    std::string grDcaZMeanNameUpgr = "meanUpgrDcaZ";
    std::string grDcaZPullNameUpgr = "pullsUpgrDcaZ";

    grDcaZResVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaZResNameCurr.c_str())));
    grDcaZResVsPtPionUpgr.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaZResNameUpgr.c_str())));
    grDcaZMeanVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaZMeanNameCurr.c_str())));
    grDcaZMeanVsPtPionUpgr.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaZMeanNameUpgr.c_str())));
    grDcaZPullVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaZPullNameCurr.c_str())));
    grDcaZPullVsPtPionUpgr.reset(dynamic_cast<TGraphErrors*>(td->Get(grDcaZPullNameUpgr.c_str())));
    if (!grDcaZResVsPtPionCurrent.get() || !grDcaZResVsPtPionUpgr.get() || !grDcaZMeanVsPtPionCurrent.get() || !grDcaZMeanVsPtPionUpgr.get() || !grDcaZPullVsPtPionCurrent.get() || !grDcaZPullVsPtPionUpgr.get()) {
      LOG(fatal) << "Something wrong with the names of the correction graphs for dcaZ. Fix it!";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void tuneTrackParams(T1 const& mcparticle, T2& trackParCov, T3 const& matCorr, T4 dcaInfoCov)
  {

    double ptMC = mcparticle.pt();

    double dcaXYResCurrent = 0.0; // sd0rpo=0.;
    double dcaZResCurrent = 0.0;  // sd0zo =0.;

    double dcaXYResUpgr = 0.0; // sd0rpn=0.;
    double dcaZResUpgr = 0.0;  // sd0zn =0.;

    // double OneOverPtCurrent = 0.0; // spt1o =0.;
    // double OneOverPtUpgr = 0.0; // spt1n =0.;

    double dcaXYMeanCurrent = 0.0; // sd0mrpo=0.;
    double dcaXYMeanUpgr = 0.0;    // sd0mrpn=0.;

    double dcaXYPullCurrent = 1.0;
    double dcaXYPullUpgr = 1.0;

    double dcaZPullCurrent = 1.0;
    double dcaZPullUpgr = 1.0;

    dcaXYResCurrent = evalGraph(ptMC, grDcaXYResVsPtPionCurrent.get());
    dcaXYResUpgr = evalGraph(ptMC, grDcaXYResVsPtPionUpgr.get());

    // dcaXYResCurrent = 1.0;
    // dcaXYResUpgr = 1.5;

    dcaZResCurrent = evalGraph(ptMC, grDcaZResVsPtPionCurrent.get());
    dcaZResUpgr = evalGraph(ptMC, grDcaZResVsPtPionUpgr.get());

    // dcaZResCurrent = 1.0;
    // dcaZResUpgr = 1.0;

    // OneOverPtCurrent = evalGraph(ptMC, grOneOverPtPionCurrent.get() );
    // OneOverPtUpgr = evalGraph(ptMC, grOneOverPtPionUpgr.get() );

    // OneOverPtCurrent = 1.0;
    // OneOverPtUpgr = 2.0;

    dcaXYMeanCurrent = evalGraph(ptMC, grDcaXYMeanVsPtPionCurrent.get());
    dcaXYMeanUpgr = evalGraph(ptMC, grDcaXYMeanVsPtPionUpgr.get());

    // dcaXYMeanCurrent = 0.0;
    // dcaXYMeanUpgr = 0.0;

    dcaXYPullCurrent = evalGraph(ptMC, grDcaXYPullVsPtPionCurrent.get());
    dcaXYPullUpgr = evalGraph(ptMC, grDcaXYPullVsPtPionUpgr.get());

    dcaZPullCurrent = evalGraph(ptMC, grDcaZPullVsPtPionCurrent.get());
    dcaZPullUpgr = evalGraph(ptMC, grDcaZPullVsPtPionUpgr.get());

    //  Unit conversion, is it required ??
    dcaXYResCurrent *= 1.e-4;
    dcaZResCurrent *= 1.e-4;

    dcaXYResUpgr *= 1.e-4;
    dcaZResUpgr *= 1.e-4;

    dcaXYMeanCurrent *= 1.e-4;
    dcaXYMeanUpgr *= 1.e-4;

    // Apply the smearing
    // ---------------------------------------------
    // double pt1o  =param  [4];
    double trackParOneOverPtCurrent = trackParCov.getQ2Pt();
    int sign = trackParCov.getQ2Pt() / std::abs(trackParCov.getQ2Pt());
    // double pt1mc =parammc[4];
    double trackParOneOverPtMC = sign / mcparticle.pt();
    o2::dataformats::VertexBase vtxMC;
    vtxMC.setPos({mcparticle.vx(), mcparticle.vy(), mcparticle.vz()});
    vtxMC.setCov(0, 0, 0, 0, 0, 0); // ??? or All ZEROs // == 1 cm2? wrt prop point

    if (debugInfo) {
      // LOG(info) << " sign " << sign;
      // LOG(info) << " trackParCov.getQ2Pt() " << trackParOneOverPtCurrent << " " << trackParCov.getQ2Pt();
      // LOG(info) << " sign/mcparticle.pt() " << trackParOneOverPtMC;
      // LOG(info) << " (curvReco-curvMC)/curvMC " << (trackParOneOverPtCurrent - trackParOneOverPtMC)/trackParOneOverPtMC * 100.0 << "%";
      // LOG(info) << " trackParCov.getPtInv() " << trackParCov.getPtInv() <<  std::endl;
      // LOG(info) << " 1/trackParCov.getPtInv() " << 1./trackParCov.getPtInv() << " & mcparticle.pt() " << mcparticle.pt();
      LOG(info) << mcparticle.pt() << " " << 1 / trackParCov.getPtInv() << " " << (trackParOneOverPtCurrent - trackParOneOverPtMC) / trackParOneOverPtMC * 100.0;
      // LOG(info) << "Before Propagation to Production Point -> alpha: " << trackParCov.getAlpha() << ", DCAxy: " << trackParCov.getY() << ", DCAz: " << trackParCov.getZ();
    }
    // propagate to DCA with respect to the Production point
    o2::base::Propagator::Instance()->propagateToDCABxByBz(vtxMC, trackParCov, 2.f, matCorr, dcaInfoCov);
    if (debugInfo) {
      // LOG(info) << "After Propagation to Production Point -> alpha: " << trackParCov.getAlpha() << ", DCAxy: " << trackParCov.getY() << ", DCAz: " << trackParCov.getZ();
      LOG(info) << "track.y(): " << trackParCov.getY();
    }

    //////////////////////////////  DCAs modifications Start  /////////////////////////////////

    // double d0zo  =param  [1];
    double trackParDcaZCurrent = trackParCov.getZ();

    // double d0rpo =param  [0];
    double trackParDcaXYCurrent = trackParCov.getY();

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
    // LOG(info) << "trackParCovMC.getY() " <<  trackParCovMC.getY()   <<  std::endl;

    // double d0zmc =parammc[1];
    double trackParDcaZMC = mcparticle.vz();

    // double dd0zo =d0zo-d0zmc;
    double diffDcaZFromMCCurent = trackParDcaZCurrent - trackParDcaZMC;

    // double dd0zn =dd0zo *(sd0zo >0. ? (sd0zn /sd0zo ) : 1.);
    double diffDcaZFromMCUpgr = diffDcaZFromMCCurent * (dcaZResCurrent > 0. ? (dcaZResUpgr / dcaZResCurrent) : 1.);

    // double d0zn  =d0zmc+dd0zn;
    double trackParDcaZUpgr = trackParDcaZMC + diffDcaZFromMCUpgr;

    // double dd0rpo=d0rpo-d0rpmc;
    double diffDcaXYFromMCCurent = trackParDcaXYCurrent - trackParDcaXYMC;

    // double dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
    double diffDcaXYFromMCUpgr = diffDcaXYFromMCCurent * (dcaXYResCurrent > 0. ? (dcaXYResUpgr / dcaXYResCurrent) : 1.);

    // double dd0mrpn=std::abs(sd0mrpn)-std::abs(sd0mrpo);
    // double diffDcaXYMeanUpgMinusCur = std::abs(dcaXYMeanUpgr) - std::abs(dcaXYMeanCurrent) ;
    double diffDcaXYMeanUpgMinusCur = dcaXYMeanUpgr - dcaXYMeanCurrent;

    // double d0rpn =d0rpmc+dd0rpn-dd0mrpn;
    double trackParDcaXYUpgr = trackParDcaXYMC + diffDcaXYFromMCUpgr - diffDcaXYMeanUpgMinusCur;

    if (debugInfo) {
      LOG(info) << dcaZResCurrent << ", " << dcaZResUpgr << ", diff(DcaZ - DcaZMC): " << diffDcaZFromMCCurent << ", diff upgraded: " << diffDcaZFromMCUpgr << ", DcaZ Upgr : " << trackParDcaZUpgr;
      LOG(info) << dcaXYResCurrent << ", " << dcaXYResUpgr << ", diff(DcaY - DcaYMC): " << diffDcaXYFromMCCurent << ", diff upgraded: " << diffDcaXYFromMCUpgr << ", DcaY Upgr :" << trackParDcaXYUpgr;
    }

    // option mimic data
    // ----------------------
    // if(fMimicData){
    //     // dd0mrpn=sd0mrpn-sd0mrpo;
    //     diffDcaXYMeanUpgMinusCur = dcaXYMeanUpgr - dcaXYMeanCurrent;
    //     // d0rpn = d0rpmc+dd0rpn+dd0mrpn;
    //     trackParDcaXYUpgr =  diffDcaXYFromMCCurent + diffDcaXYFromMCUpgr + diffDcaXYMeanUpgMinusCur;
    // }

    // setting updated track parameters
    // --------------------------------
    // param[0]=d0rpn;
    // double oldDCAxyValue = trackParCov.getY();
    trackParCov.setY(trackParDcaXYUpgr);
    // trackParCov.setY(oldDCAxyValue);
    // param[1]=d0zn ;
    trackParCov.setZ(trackParDcaZUpgr);

    if (updateCurvature) {
      // --------------------------------------
      // double dpt1o =pt1o-pt1mc;
      double diffOneOverPtFromMCCurent = trackParOneOverPtCurrent - trackParOneOverPtMC;

      // double dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
      double diffOneOverPtFromMCUpgr = diffOneOverPtFromMCCurent * (oneOverPtCurrent > 0. ? (oneOverPtUpgr / oneOverPtCurrent) : 1.);

      // double pt1n  = pt1mc+dpt1n;
      double trackParOneOverPtUpgr = trackParOneOverPtMC + diffOneOverPtFromMCUpgr;

      // param[4]=pt1n ;
      trackParCov.setQ2Pt(trackParOneOverPtUpgr);
    }
    // if(debugInfo){
    //  LOG(info) << "Inside tuneTrackParams() before modifying trackParCov.getY(): " << trackParCov.getY() << " trackParOneOverPtMC = " << trackParOneOverPtMC << " diffOneOverPtFromMCUpgr = "  << diffOneOverPtFromMCUpgr <<  std::endl;
    //}

    // Updating Single Track Covariance matrices

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

    if (updateTrackCovMat) {
      //       if(sd0rpo>0.)            covar[0]*=(sd0rpn/sd0rpo)*(sd0rpn/sd0rpo);//yy
      sigmaY2 = trackParCov.getSigmaY2();
      if (dcaXYResCurrent > 0.)
        sigmaY2 *= ((dcaXYResUpgr / dcaXYResCurrent) * (dcaXYResUpgr / dcaXYResCurrent));
      trackParCov.setCov(sigmaY2, 0);

      //       if(sd0zo>0. && sd0rpo>0.)covar[1]*=(sd0rpn/sd0rpo)*(sd0zn/sd0zo);//yz
      sigmaZY = trackParCov.getSigmaZY();
      if (dcaZResCurrent > 0. && dcaXYResCurrent > 0.)
        sigmaZY *= ((dcaXYResUpgr / dcaXYResCurrent) * (dcaXYResUpgr / dcaXYResCurrent));
      trackParCov.setCov(sigmaZY, 1);

      //       if(sd0zo>0.)             covar[2]*=(sd0zn/sd0zo)*(sd0zn/sd0zo);//zz
      sigmaZ2 = trackParCov.getSigmaZ2();
      if (dcaZResCurrent > 0.)
        sigmaZ2 *= ((dcaZResUpgr / dcaZResCurrent) * (dcaZResUpgr / dcaZResCurrent));
      trackParCov.setCov(sigmaZ2, 2);

      //       if(sd0rpo>0.)            covar[3]*=(sd0rpn/sd0rpo);//yl
      sigmaSnpY = trackParCov.getSigmaSnpY();
      if (dcaXYResCurrent > 0.)
        sigmaSnpY *= ((dcaXYResUpgr / dcaXYResCurrent));
      trackParCov.setCov(sigmaSnpY, 3);

      //       if(sd0zo>0.)             covar[4]*=(sd0zn/sd0zo);//zl
      sigmaSnpZ = trackParCov.getSigmaSnpZ();
      if (dcaZResCurrent > 0.)
        sigmaSnpZ *= ((dcaZResUpgr / dcaZResCurrent));
      trackParCov.setCov(sigmaSnpZ, 4);

      //       if(sd0rpo>0.)            covar[6]*=(sd0rpn/sd0rpo);//ysenT
      sigmaTglY = trackParCov.getSigmaTglY();
      if (dcaXYResCurrent > 0.)
        sigmaTglY *= ((dcaXYResUpgr / dcaXYResCurrent));
      trackParCov.setCov(sigmaTglY, 6);

      //       if(sd0zo>0.)             covar[7]*=(sd0zn/sd0zo);//zsenT
      sigmaTglZ = trackParCov.getSigmaTglZ();
      if (dcaZResCurrent > 0.)
        sigmaTglZ *= ((dcaZResUpgr / dcaZResCurrent));
      trackParCov.setCov(sigmaTglZ, 7);

      //       if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
      sigma1PtY = trackParCov.getSigma1PtY();
      if (dcaXYResCurrent > 0. && oneOverPtCurrent > 0.)
        sigma1PtY *= ((dcaXYResUpgr / dcaXYResCurrent) * (oneOverPtUpgr / oneOverPtCurrent));
      trackParCov.setCov(sigma1PtY, 10);

      //       if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
      sigma1PtZ = trackParCov.getSigma1PtZ();
      if (dcaZResCurrent > 0. && oneOverPtCurrent > 0.)
        sigma1PtZ *= ((dcaZResUpgr / dcaZResCurrent) * (oneOverPtUpgr / oneOverPtCurrent));
      trackParCov.setCov(sigma1PtZ, 11);

      //       if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
      sigma1PtSnp = trackParCov.getSigma1PtSnp();
      if (oneOverPtCurrent > 0.)
        sigma1PtSnp *= (oneOverPtUpgr / oneOverPtCurrent);
      trackParCov.setCov(sigma1PtSnp, 12);

      //       if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
      sigma1PtTgl = trackParCov.getSigma1PtTgl();
      if (oneOverPtCurrent > 0.)
        sigma1PtTgl *= (oneOverPtUpgr / oneOverPtCurrent);
      trackParCov.setCov(sigma1PtTgl, 13);

      //       if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
      sigma1Pt2 = trackParCov.getSigma1Pt2();
      if (oneOverPtCurrent > 0.)
        sigma1Pt2 *= (oneOverPtUpgr / oneOverPtCurrent);
      trackParCov.setCov(sigma1Pt2, 14);
    }

    if (updatePulls) {
      double ratioDCAxyPulls = dcaXYPullCurrent / dcaXYPullUpgr;
      double ratioDCAzPulls = dcaZPullCurrent / dcaZPullUpgr;

      // covar[0]*=pullcorr*pullcorr;//yy

      sigmaY2 *= (ratioDCAxyPulls * ratioDCAxyPulls);
      trackParCov.setCov(sigmaY2, 0);

      // covar[1]*=pullcorr;//yz
      sigmaZY *= ratioDCAxyPulls;
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
    }
  }

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
      printf("\tevalGraph fails !\n");
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
