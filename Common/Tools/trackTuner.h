
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include <TGraphErrors.h>
#include <string>

using namespace o2;
using namespace o2::framework;

struct TrackTuner {

  ///////////////////////////////
  /// parameters to be configured
  bool debugInfo = false;
  bool updateTrackCovMat = false;
  bool updateCurvature = false;
  bool updatePulls = false;
  std::string pathCurrFileDcaXY = ""; // Path to file containing current DCAxy graphs
  std::string pathUpgrFileDcaXY = ""; // Path to file containing New DCAxy graphs
  std::string pathCurrFileDcaZ = ""; // Path to file containing current DCAz graphs
  std::string pathUpgrFileDcaZ = ""; // sPath to file containing New DCAz graphs
  float oneOverPtCurrent = 0.; // 1/pt old
  float oneOverPtUpgrded = 0.; // 1/pt new
  ///////////////////////////////

  std::unique_ptr<TGraphErrors> grDcaXYResVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaXYResVsPtPionUpgrded;

  std::unique_ptr<TGraphErrors> grDcaZResVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaZResVsPtPionUpgrded;

  std::unique_ptr<TGraphErrors> grDcaXYMeanVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaXYMeanVsPtPionUpgrded;

  std::unique_ptr<TGraphErrors> grDcaZMeanVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaZMeanVsPtPionUpgrded;

  std::unique_ptr<TGraphErrors> grOneOverPtPionCurrent;
  std::unique_ptr<TGraphErrors> grOneOverPtPionUpgrded;

  std::unique_ptr<TGraphErrors> grDcaXYPullVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaXYPullVsPtPionUpgrded;

  std::unique_ptr<TGraphErrors> grDcaZPullVsPtPionCurrent;
  std::unique_ptr<TGraphErrors> grDcaZPullVsPtPionUpgrded;


  /// @brief Function to configure the TrackTuner parameters
  /// @param inputString Input string with all parameter configuration. Format: <name>=<value>|<name>=<value>
  /// @return String with the values of all parameters after configurations are listed, to cross check that everything worked well
  std::string configParams(std::string inputString) {
    
    std::string delimiter = "|";
    std::string assignmentSymbol = "=";
    
    LOG(info) << "[TrackTuner] === ";
    LOG(info) << "[TrackTuner] === Parameter configuration via std::string";
    LOG(info) << "[TrackTuner] === Required format: \"<name>" << assignmentSymbol <<"<value>" << delimiter << "<name>" << assignmentSymbol << "<value>" << delimiter << "<name>" << assignmentSymbol << "<value>\"";
    LOG(info) << "[TrackTuner] === Delimiter symbol: \"" << delimiter << "\"";
    LOG(info) << "[TrackTuner] === Assignment symbol: \"" << assignmentSymbol << "\"";
    LOG(info) << "[TrackTuner] === ";
    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> Original input string = \"" << inputString << "\"";

    /// Check the format of the input string
    if(inputString.find(delimiter) == std::string::npos) {
        // wrong delimiter symbol used
        LOG(fatal) << "ERROR: delimiter symbol \"" << delimiter << "\" not found in the configuration string. Fix it!";
    }
    if(inputString.find(assignmentSymbol) == std::string::npos) {
        // wrong assignment symbol used
        LOG(fatal) << "ERROR: assignment symbol \"" << assignmentSymbol << "\" not found in the configuration string. Fix it!";
    }
    int spaces = std::count(inputString.begin(), inputString.end(), ' ');
    if(spaces > 0) {
        // white spaces to be removed
        LOG(fatal) << "ERROR: " << spaces << " white spaces found in the configuration string. Remove them!";
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// Parameters to be configured via the string
    /// +++ to be manually updated every time one adds a new parameter to the TrackTuner.h +++
    enum ePars {eDebugInfo=0, eUpdateTrackCovMat, eUpdateCurvature, eUpdatePulls, ePathCurrFileDcaXY, ePathUpgrFileDcaXY, ePathCurrFileDcaZ, ePathUpgrFileDcaZ, eOneOverPtCurrent, eOneOverPtUpgrded, eNPars};
    std::map<uint8_t, std::string> mapParNames = {
         std::make_pair(static_cast<uint8_t>(eDebugInfo), "debugInfo")
        ,std::make_pair(static_cast<uint8_t>(eUpdateTrackCovMat), "updateTrackCovMat")
        ,std::make_pair(static_cast<uint8_t>(eUpdateCurvature), "updateCurvature")
        ,std::make_pair(static_cast<uint8_t>(eUpdatePulls), "updatePulls")
        ,std::make_pair(static_cast<uint8_t>(ePathCurrFileDcaXY), "pathCurrFileDcaXY")
        ,std::make_pair(static_cast<uint8_t>(ePathUpgrFileDcaXY), "pathUpgrFileDcaXY")
        ,std::make_pair(static_cast<uint8_t>(ePathCurrFileDcaZ), "pathCurrFileDcaZ")
        ,std::make_pair(static_cast<uint8_t>(ePathUpgrFileDcaZ), "pathUpgrFileDcaZ")
        ,std::make_pair(static_cast<uint8_t>(eOneOverPtCurrent), "oneOverPtCurrent")
        ,std::make_pair(static_cast<uint8_t>(eOneOverPtUpgrded), "oneOverPtUpgrded")
    };
    ///////////////////////////////////////////////////////////////////////////////////
    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> Parameters before the custom settings";
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    LOG(info) << "[TrackTuner]     pathCurrFileDcaXY = " << pathCurrFileDcaXY;
    LOG(info) << "[TrackTuner]     pathUpgrFileDcaXY = " << pathUpgrFileDcaXY;
    LOG(info) << "[TrackTuner]     pathCurrFileDcaZ = " << pathCurrFileDcaZ;
    LOG(info) << "[TrackTuner]     pathUpgrFileDcaZ = " << pathUpgrFileDcaZ;
    LOG(info) << "[TrackTuner]     oneOverPtCurrent = " << oneOverPtCurrent;
    LOG(info) << "[TrackTuner]     oneOverPtUpgrded = " << oneOverPtUpgrded;


    //##############################################################################################
    //########   split the original string, separating substrings delimited by "|" symbol   ########

    std::vector<std::string> slices = {};

    while(inputString.find(delimiter) != std::string::npos) {  
        /// we are not at the end of our string --> let's find out the next parameter to be configured!
        slices.push_back( inputString.substr(0, inputString.find(delimiter)) ); // this gives us the substring until the next "|" character
        inputString.erase(0, slices.back().length() + delimiter.length()); // erase the found substring for next iteration
    }
    /// at this point, the input string is erased until the last delimiter (included)
    slices.push_back(inputString); // necessary to read the last parameter, after the last "|" symbol

    LOG(info) << "[TrackTuner]";
    LOG(info) << "[TrackTuner] >>> String slices:";
    for(std::string& s : slices)    LOG(info) << "[TrackTuner]     " << s;

    /// check if the number of input parameters is correct
    if (static_cast<uint8_t>(slices.size()) != eNPars) {
        LOG(fatal) << "[TrackTuner] " << slices.size() << " parameters provided, while " << eNPars << " are expected. Fix it!";
    }

    //###################################################################################################################
    //########   each split is now a std::string "<parName>=<value>" --> let's really configure each parameter   ########

    /// lambda expression to search for the parameter value (as string) in the configuration string
    auto getValueString = [&](uint8_t iPar) {
        
        /// this allows to search the parameter configuration even if they are not written in order
        auto it = std::find_if(slices.begin(), slices.end(), [&](std::string s) {return s.find(mapParNames[iPar]) != std::string::npos;});
        if(it == std::end(slices)) {
          // parameter not found
          LOG(fatal) << "\"" << mapParNames[iPar] << "\" not found in the configuration string";
        }
        std::string str = *it;
        if(str.find('=') == std::string::npos || str.back() == '=') {
          // value of the parameter missing in the configuration string
          LOG(fatal) << "Missing value for \"" << mapParNames[iPar] << "\" in the configuration string";
        }
        return str.substr(str.find(assignmentSymbol)+1, str.length());
    };

    /// further lambda expression to handle bool initialization
    auto setBoolFromString = [=](bool& b, std::string str) {
        if(!str.compare("1") || str.find("true") != std::string::npos || str.find("True") != std::string::npos || str.find("TRUE") != std::string::npos ) {
            b = true;
        } else if (!str.compare("0") || str.find("false") != std::string::npos || str.find("False") != std::string::npos || str.find("FALSE") != std::string::npos ) {
            b = false;
        } else {
            LOG(fatal) << "[TrackTuner] Wrong bool initialization from configuration ";
        }
    };

    std::string outputString = "";
    LOG(info) << "[TrackTuner] ";
    LOG(info) << "[TrackTuner] >>> Parameters after the custom settings";
    // Configure debugInfo
    setBoolFromString(debugInfo, getValueString(eDebugInfo));
    LOG(info) << "[TrackTuner]     debugInfo = " << debugInfo;
    outputString += "debugInfo=" + std::to_string(debugInfo);
    // Configure updateTrackCovMat
    setBoolFromString(updateTrackCovMat, getValueString(eUpdateTrackCovMat));
    LOG(info) << "[TrackTuner]     updateTrackCovMat = " << updateTrackCovMat;
    outputString += ", updateTrackCovMat=" + std::to_string(updateTrackCovMat);
    // Configure updateCurvature
    setBoolFromString(updateCurvature, getValueString(eUpdateCurvature));
    LOG(info) << "[TrackTuner]     updateCurvature = " << updateCurvature;
    outputString += ", updateCurvature=" + std::to_string(updateCurvature);
    // Configure updatePulls
    setBoolFromString(updatePulls, getValueString(eUpdatePulls));
    LOG(info) << "[TrackTuner]     updatePulls = " << updatePulls;
    outputString += ", updatePulls=" + std::to_string(updatePulls);
    // Configure pathCurrFileDcaXY
    pathCurrFileDcaXY = getValueString(ePathCurrFileDcaXY);
    outputString += ", pathCurrFileDcaXY=" + pathCurrFileDcaXY;
    LOG(info) << "[TrackTuner]     pathCurrFileDcaXY = " << pathCurrFileDcaXY;
    // Configure pathUpgrFileDcaXY
    pathUpgrFileDcaXY = getValueString(ePathUpgrFileDcaXY);
    outputString += ", pathUpgrFileDcaXY=" + pathUpgrFileDcaXY;
    LOG(info) << "[TrackTuner]     pathUpgrFileDcaXY = " << pathUpgrFileDcaXY;
    // Configure pathCurrFileDcaZ
    pathCurrFileDcaZ = getValueString(ePathCurrFileDcaZ);
    outputString += ", pathCurrFileDcaZ=" + pathCurrFileDcaZ;
    LOG(info) << "[TrackTuner]     pathCurrFileDcaZ = " << pathCurrFileDcaZ;
    // Configure pathUpgrFileDcaZ
    pathUpgrFileDcaZ = getValueString(ePathUpgrFileDcaZ);
    outputString += ", pathUpgrFileDcaZ=" + pathUpgrFileDcaZ;
    LOG(info) << "[TrackTuner]     pathUpgrFileDcaZ = " << pathUpgrFileDcaZ;
    // Configure oneOverPtCurr 
    oneOverPtCurrent = std::stof(getValueString(eOneOverPtCurrent));
    outputString += ", oneOverPtCurrent=" + std::to_string(oneOverPtCurrent);
    LOG(info) << "[TrackTuner]     oneOverPtCurrent = " << oneOverPtCurrent;
    // Configure oneOverPtUpgrded
    oneOverPtUpgrded = std::stof(getValueString(eOneOverPtUpgrded));
    outputString += ", oneOverPtUpgrded=" + std::to_string(oneOverPtUpgrded);
    LOG(info) << "[TrackTuner]     oneOverPtUpgrded = " << oneOverPtUpgrded;

    return outputString;    
  }


  void getDcaGraphs()
  {

    const std::string fnameDCAxyFileCurr = pathCurrFileDcaXY.data();
    const std::string fnameDCAxyFileUpgr = pathUpgrFileDcaXY.data();
    const std::string fnameDCAzFileCurr = pathCurrFileDcaZ.data();
    const std::string fnameDCAzFileUpgr = pathUpgrFileDcaZ.data();

    /*
    TODO
    --> add possibility to pick-up the file from CCDB
    */

    if ((fnameDCAxyFileCurr == "") || (fnameDCAzFileCurr == "") || (fnameDCAxyFileUpgr == "") || (fnameDCAzFileUpgr == "")) {
      LOG(fatal) << "[TrackTuner] No Correction DCA files provided!";
      return;
    }

    std::unique_ptr<TFile> fileCurrDcaXY(TFile::Open(fnameDCAxyFileCurr.c_str(), "READ"));
    std::unique_ptr<TFile> fileUpgrDcaXY(TFile::Open(fnameDCAxyFileUpgr.c_str(), "READ"));

    std::unique_ptr<TFile> fileCurrDcaZ(TFile::Open(fnameDCAzFileCurr.c_str(), "READ"));
    std::unique_ptr<TFile> fileUpgrDcaZ(TFile::Open(fnameDCAzFileUpgr.c_str(), "READ"));

    if(!fileCurrDcaXY.get() || !fileUpgrDcaXY.get() || !fileCurrDcaZ.get() || !fileUpgrDcaZ.get()) {
      LOG(fatal) << "Something wrong with the input files for dca correction. Fix it!";
    }

    std::string grDcaResName = "tge_DCA_res_withoutPVrefit_all";
    std::string grDcaMeanName = "tge_DCA_mean_withoutPVrefit_all";
    std::string grDcaPullName = "tge_DCAPulls_res_withoutPVrefit_all";

    grDcaXYResVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(fileCurrDcaXY->Get(grDcaResName.c_str())));
    grDcaXYResVsPtPionUpgrded.reset(dynamic_cast<TGraphErrors*>(fileUpgrDcaXY->Get(grDcaResName.c_str())));
    grDcaXYMeanVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(fileCurrDcaXY->Get(grDcaMeanName.c_str())));
    grDcaXYMeanVsPtPionUpgrded.reset(dynamic_cast<TGraphErrors*>(fileUpgrDcaXY->Get(grDcaMeanName.c_str())));
    grDcaXYPullVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(fileCurrDcaXY->Get(grDcaPullName.c_str())));
    grDcaXYPullVsPtPionUpgrded.reset(dynamic_cast<TGraphErrors*>(fileUpgrDcaXY->Get(grDcaPullName.c_str())));
    if(!grDcaXYResVsPtPionCurrent.get() || !grDcaXYResVsPtPionUpgrded.get() || !grDcaXYMeanVsPtPionCurrent.get() || !grDcaXYMeanVsPtPionUpgrded.get() || !grDcaXYPullVsPtPionCurrent.get() || !grDcaXYPullVsPtPionUpgrded.get()) {
      LOG(fatal) << "Something wrong with the names of the correction graphs for dcaXY. Fix it!";
    }

    grDcaZResVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(fileCurrDcaZ->Get(grDcaResName.c_str())));
    grDcaZResVsPtPionUpgrded.reset(dynamic_cast<TGraphErrors*>(fileUpgrDcaZ->Get(grDcaResName.c_str())));
    grDcaZMeanVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(fileCurrDcaZ->Get(grDcaMeanName.c_str())));
    grDcaZMeanVsPtPionUpgrded.reset(dynamic_cast<TGraphErrors*>(fileUpgrDcaZ->Get(grDcaMeanName.c_str())));
    grDcaZPullVsPtPionCurrent.reset(dynamic_cast<TGraphErrors*>(fileCurrDcaZ->Get(grDcaPullName.c_str())));
    grDcaZPullVsPtPionUpgrded.reset(dynamic_cast<TGraphErrors*>(fileUpgrDcaZ->Get(grDcaPullName.c_str())));
    if(!grDcaZResVsPtPionCurrent.get() || !grDcaZResVsPtPionUpgrded.get() || !grDcaZMeanVsPtPionCurrent.get() || !grDcaZMeanVsPtPionUpgrded.get() || !grDcaZPullVsPtPionCurrent.get() || !grDcaZPullVsPtPionUpgrded.get()) {
      LOG(fatal) << "Something wrong with the names of the correction graphs for dcaZ. Fix it!";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void tuneTrackParams(T1 mcparticle, T2& trackParCov, T3 matCorr, T4 dcaInfoCov)
  {

    double ptMC = TMath::Abs(mcparticle.pt());

    double DcaXYResCurrent = 0.0; // sd0rpo=0.;
    double DcaZResCurrent = 0.0;  // sd0zo =0.;

    double DcaXYResUpgrded = 0.0; // sd0rpn=0.;
    double DcaZResUpgrded = 0.0;  // sd0zn =0.;

    //double OneOverPtCurrent = 0.0; // spt1o =0.;
    //double OneOverPtUpgrded = 0.0; // spt1n =0.;

    double DcaXYMeanCurrent = 0.0; // sd0mrpo=0.;
    double DcaXYMeanUpgrded = 0.0; // sd0mrpn=0.;

    double DcaXYPullCurrent = 1.0;
    double DcaXYPullUpgrded = 1.0;

    double DcaZPullCurrent = 1.0;
    double DcaZPullUpgrded = 1.0;

    DcaXYResCurrent = EvalGraph(ptMC, grDcaXYResVsPtPionCurrent.get());
    DcaXYResUpgrded = EvalGraph(ptMC, grDcaXYResVsPtPionUpgrded.get());

    // DcaXYResCurrent = 1.0;
    // DcaXYResUpgrded = 1.5;

    DcaZResCurrent = EvalGraph(ptMC, grDcaZResVsPtPionCurrent.get());
    DcaZResUpgrded = EvalGraph(ptMC, grDcaZResVsPtPionUpgrded.get());

    // DcaZResCurrent = 1.0;
    // DcaZResUpgrded = 1.0;

    // OneOverPtCurrent = EvalGraph(ptMC, grOneOverPtPionCurrent.get() );
    // OneOverPtUpgrded = EvalGraph(ptMC, grOneOverPtPionUpgrded.get() );

    // OneOverPtCurrent = 1.0;
    // OneOverPtUpgrded = 2.0;

    DcaXYMeanCurrent = EvalGraph(ptMC, grDcaXYMeanVsPtPionCurrent.get());
    DcaXYMeanUpgrded = EvalGraph(ptMC, grDcaXYMeanVsPtPionUpgrded.get());

    // DcaXYMeanCurrent = 0.0;
    // DcaXYMeanUpgrded = 0.0;

    DcaXYPullCurrent = EvalGraph(ptMC, grDcaXYPullVsPtPionCurrent.get());
    DcaXYPullUpgrded = EvalGraph(ptMC, grDcaXYPullVsPtPionUpgrded.get());

    DcaZPullCurrent = EvalGraph(ptMC, grDcaZPullVsPtPionCurrent.get());
    DcaZPullUpgrded = EvalGraph(ptMC, grDcaZPullVsPtPionUpgrded.get());

    //  Unit conversion, is it required ??
    DcaXYResCurrent *= 1.e-4;
    DcaZResCurrent *= 1.e-4;

    DcaXYResUpgrded *= 1.e-4;
    DcaZResUpgrded *= 1.e-4;

    DcaXYMeanCurrent *= 1.e-4;
    DcaXYMeanUpgrded *= 1.e-4;

    // Apply the smearing
    // ---------------------------------------------
    // double pt1o  =param  [4];
    double trackParOneOverPtCurrent = trackParCov.getQ2Pt();
    int sign = trackParCov.getQ2Pt() / TMath::Abs(trackParCov.getQ2Pt());

    // double pt1mc =parammc[4];
    double trackParOneOverPtMC = sign / mcparticle.pt();

    //LOG(info) << std::scientific << std::setprecision(5);
    o2::dataformats::VertexBase vtxMC;
    vtxMC.setPos({mcparticle.vx(), mcparticle.vy(), mcparticle.vz()});
    vtxMC.setCov(0, 0, 0, 0, 0, 0); // ??? or All ZEROs // == 1 cm2? wrt prop point

    if(debugInfo){
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

    if(debugInfo){
      // LOG(info) << "After Propagation to Production Point -> alpha: " << trackParCov.getAlpha() << ", DCAxy: " << trackParCov.getY() << ", DCAz: " << trackParCov.getZ();

      LOG(info) << "track.y(): " << trackParCov.getY();
    }

    //////////////////////////////  DCAs modifications Start  /////////////////////////////////

    // double d0zo  =param  [1];
    double trackParDcaZCurrent = trackParCov.getZ();

    // double d0rpo =param  [0];
    double trackParDcaXYCurrent = trackParCov.getY();

    float mcVxRotated = mcparticle.vx() * TMath::Cos(trackParCov.getAlpha()) + mcparticle.vy() * TMath::Sin(trackParCov.getAlpha()); // invert
    float mcVyRotated = mcparticle.vy() * TMath::Cos(trackParCov.getAlpha()) - mcparticle.vx() * TMath::Sin(trackParCov.getAlpha());

    if(debugInfo){
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
    double diffDcaZFromMCUpgrded = diffDcaZFromMCCurent * (DcaZResCurrent > 0. ? (DcaZResUpgrded / DcaZResCurrent) : 1.);

    // double d0zn  =d0zmc+dd0zn;
    double trackParDcaZUpgrded = trackParDcaZMC + diffDcaZFromMCUpgrded;

    // double dd0rpo=d0rpo-d0rpmc;
    double diffDcaXYFromMCCurent = trackParDcaXYCurrent - trackParDcaXYMC;

    // double dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
    double diffDcaXYFromMCUpgrded = diffDcaXYFromMCCurent * (DcaXYResCurrent > 0. ? (DcaXYResUpgrded / DcaXYResCurrent) : 1.);

    // double dd0mrpn=TMath::Abs(sd0mrpn)-TMath::Abs(sd0mrpo);
    // double diffDcaXYMeanUpgMinusCur = TMath::Abs(DcaXYMeanUpgrded) - TMath::Abs(DcaXYMeanCurrent) ;
    double diffDcaXYMeanUpgMinusCur = DcaXYMeanUpgrded - DcaXYMeanCurrent;

    // double d0rpn =d0rpmc+dd0rpn-dd0mrpn;
    double trackParDcaXYUpgrded = trackParDcaXYMC + diffDcaXYFromMCUpgrded - diffDcaXYMeanUpgMinusCur;

    if(debugInfo){
      LOG(info) << DcaZResCurrent << ", " << DcaZResUpgrded << ", diff(DcaZ - DcaZMC): " << diffDcaZFromMCCurent << ", diff upgraded: " << diffDcaZFromMCUpgrded << ", DcaZ Upgrded : " << trackParDcaZUpgrded;
      LOG(info) << DcaXYResCurrent << ", " << DcaXYResUpgrded << ", diff(DcaY - DcaYMC): " << diffDcaXYFromMCCurent << ", diff upgraded: " << diffDcaXYFromMCUpgrded << ", DcaY Upgrded :" << trackParDcaXYUpgrded;
    }

    // option mimic data
    // ----------------------
    // if(fMimicData){
    //     // dd0mrpn=sd0mrpn-sd0mrpo;
    //     diffDcaXYMeanUpgMinusCur = DcaXYMeanUpgrded - DcaXYMeanCurrent;
    //     // d0rpn = d0rpmc+dd0rpn+dd0mrpn;
    //     trackParDcaXYUpgrded =  diffDcaXYFromMCCurent + diffDcaXYFromMCUpgrded + diffDcaXYMeanUpgMinusCur;
    // }

    // setting updated track parameters
    // --------------------------------
    // param[0]=d0rpn;
    // double oldDCAxyValue = trackParCov.getY();
    trackParCov.setY(trackParDcaXYUpgrded);
    // trackParCov.setY(oldDCAxyValue);
    // param[1]=d0zn ;
    trackParCov.setZ(trackParDcaZUpgrded);

    if (updateCurvature) {
      // --------------------------------------
      // double dpt1o =pt1o-pt1mc;
      double diffOneOverPtFromMCCurent = trackParOneOverPtCurrent - trackParOneOverPtMC;

      // double dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
      double diffOneOverPtFromMCUpgrded = diffOneOverPtFromMCCurent * (oneOverPtCurrent > 0. ? (oneOverPtUpgrded / oneOverPtCurrent) : 1.);

      // double pt1n  = pt1mc+dpt1n;
      double trackParOneOverPtUpgrded = trackParOneOverPtMC + diffOneOverPtFromMCUpgrded;

      // param[4]=pt1n ;
      trackParCov.setQ2Pt(trackParOneOverPtUpgrded);
    }
    //if(debugInfo){
      // LOG(info) << "Inside tuneTrackParams() before modifying trackParCov.getY(): " << trackParCov.getY() << " trackParOneOverPtMC = " << trackParOneOverPtMC << " diffOneOverPtFromMCUpgrded = "  << diffOneOverPtFromMCUpgrded <<  std::endl;
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
      if (DcaXYResCurrent > 0.)
        sigmaY2 *= ((DcaXYResUpgrded / DcaXYResCurrent) * (DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaY2, 0);

      //       if(sd0zo>0. && sd0rpo>0.)covar[1]*=(sd0rpn/sd0rpo)*(sd0zn/sd0zo);//yz
      sigmaZY = trackParCov.getSigmaZY();
      if (DcaZResCurrent > 0. && DcaXYResCurrent > 0.)
        sigmaZY *= ((DcaXYResUpgrded / DcaXYResCurrent) * (DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaZY, 1);

      //       if(sd0zo>0.)             covar[2]*=(sd0zn/sd0zo)*(sd0zn/sd0zo);//zz
      sigmaZ2 = trackParCov.getSigmaZ2();
      if (DcaZResCurrent > 0.)
        sigmaZ2 *= ((DcaZResUpgrded / DcaZResCurrent) * (DcaZResUpgrded / DcaZResCurrent));
      trackParCov.setCov(sigmaZ2, 2);

      //       if(sd0rpo>0.)            covar[3]*=(sd0rpn/sd0rpo);//yl
      sigmaSnpY = trackParCov.getSigmaSnpY();
      if (DcaXYResCurrent > 0.)
        sigmaSnpY *= ((DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaSnpY, 3);

      //       if(sd0zo>0.)             covar[4]*=(sd0zn/sd0zo);//zl
      sigmaSnpZ = trackParCov.getSigmaSnpZ();
      if (DcaZResCurrent > 0.)
        sigmaSnpZ *= ((DcaZResUpgrded / DcaZResCurrent));
      trackParCov.setCov(sigmaSnpZ, 4);

      //       if(sd0rpo>0.)            covar[6]*=(sd0rpn/sd0rpo);//ysenT
      sigmaTglY = trackParCov.getSigmaTglY();
      if (DcaXYResCurrent > 0.)
        sigmaTglY *= ((DcaXYResUpgrded / DcaXYResCurrent));
      trackParCov.setCov(sigmaTglY, 6);

      //       if(sd0zo>0.)             covar[7]*=(sd0zn/sd0zo);//zsenT
      sigmaTglZ = trackParCov.getSigmaTglZ();
      if (DcaZResCurrent > 0.)
        sigmaTglZ *= ((DcaZResUpgrded / DcaZResCurrent));
      trackParCov.setCov(sigmaTglZ, 7);

      //       if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
      sigma1PtY = trackParCov.getSigma1PtY();
      if (DcaXYResCurrent > 0. && oneOverPtCurrent > 0.)
        sigma1PtY *= ((DcaXYResUpgrded / DcaXYResCurrent) * (oneOverPtUpgrded / oneOverPtCurrent));
      trackParCov.setCov(sigma1PtY, 10);

      //       if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
      sigma1PtZ = trackParCov.getSigma1PtZ();
      if (DcaZResCurrent > 0. && oneOverPtCurrent > 0.)
        sigma1PtZ *= ((DcaZResUpgrded / DcaZResCurrent) * (oneOverPtUpgrded / oneOverPtCurrent));
      trackParCov.setCov(sigma1PtZ, 11);

      //       if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
      sigma1PtSnp = trackParCov.getSigma1PtSnp();
      if (oneOverPtCurrent > 0.)
        sigma1PtSnp *= (oneOverPtUpgrded / oneOverPtCurrent);
      trackParCov.setCov(sigma1PtSnp, 12);

      //       if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
      sigma1PtTgl = trackParCov.getSigma1PtTgl();
      if (oneOverPtCurrent > 0.)
        sigma1PtTgl *= (oneOverPtUpgrded / oneOverPtCurrent);
      trackParCov.setCov(sigma1PtTgl, 13);

      //       if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
      sigma1Pt2 = trackParCov.getSigma1Pt2();
      if (oneOverPtCurrent > 0.)
        sigma1Pt2 *= (oneOverPtUpgrded / oneOverPtCurrent);
      trackParCov.setCov(sigma1Pt2, 14);
    }

    if (updatePulls) {
      double ratioDCAxyPulls = DcaXYPullCurrent / DcaXYPullUpgrded;
      double ratioDCAzPulls = DcaZPullCurrent / DcaZPullUpgrded;

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
  int getPhiBin(double phi) const
  {
    double pi = TMath::Pi();
    if (phi > 2. * pi || phi < 0.)
      return -1;
    if ((phi <= (pi / 4.)) || (phi > 7. * (pi / 4.)))
      return 0;
    if ((phi > (pi / 4.)) && (phi <= 3. * (pi / 4.)))
      return 1;
    if ((phi > 3. * (pi / 4.)) && (phi <= 5. * (pi / 4.)))
      return 2;
    if ((phi > (5. * pi / 4.)) && (phi <= 7. * (pi / 4.)))
      return 3;

    return -1;
  }

  double EvalGraph(double x, const TGraphErrors* graph) const
  {

    if (!graph) {
      printf("\tEvalGraph fails !\n");
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
