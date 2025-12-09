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
//
// \author: prottay das 23/12/2024
// \email: prottay.das@cern.ch

#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h>
#include <TComplex.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom3.h>

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::rctsel;

using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

struct spvector {

  Produces<aod::SPCalibrationTables> spcalibrationtable;

  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Centrality cut Max"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Centrality cut Min"};
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> usemem{"usemem", true, "usemem"};

  struct : ConfigurableGroup {
    Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
    Configurable<int> PhiNbins{"PhiNbins", 100, "Number of bins in phi histogram"};
    Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
    Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};
    Configurable<int> VxNbins{"VxNbins", 25, "Number of bins in Vx histograms"};
    Configurable<float> lbinVx{"lbinVx", -0.05, "lower bin value in Vx histograms"};
    Configurable<float> hbinVx{"hbinVx", 0.0, "higher bin value in Vx histograms"};
    Configurable<int> VyNbins{"VyNbins", 25, "Number of bins in Vy histograms"};
    Configurable<float> lbinVy{"lbinVy", -0.02, "lower bin value in Vy histograms"};
    Configurable<float> hbinVy{"hbinVy", 0.02, "higher bin value in Vy histograms"};
    Configurable<int> VzNbins{"VzNbins", 20, "Number of bins in Vz histograms"};
    Configurable<float> lbinVz{"lbinVz", -10.0, "lower bin value in Vz histograms"};
    Configurable<float> hbinVz{"hbinVz", 10.0, "higher bin value in Vz histograms"};
    Configurable<int> CentNbins{"CentNbins", 16, "Number of bins in cent histograms"};
    Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
    Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
    Configurable<int> VxfineNbins{"VxfineNbins", 25, "Number of bins in Vx fine histograms"};
    Configurable<float> lfinebinVx{"lfinebinVx", -0.05, "lower bin value in Vx fine histograms"};
    Configurable<float> hfinebinVx{"hfinebinVx", 0.0, "higher bin value in Vx fine histograms"};
    Configurable<int> VyfineNbins{"VyfineNbins", 25, "Number of bins in Vy fine histograms"};
    Configurable<float> lfinebinVy{"lfinebinVy", -0.02, "lower bin value in Vy fine histograms"};
    Configurable<float> hfinebinVy{"hfinebinVy", 0.02, "higher bin value in Vy fine histograms"};
    Configurable<int> VzfineNbins{"VzfineNbins", 20, "Number of bins in Vz fine histograms"};
    Configurable<float> lfinebinVz{"lfinebinVz", -10.0, "lower bin value in Vz fine histograms"};
    Configurable<float> hfinebinVz{"hfinebinVz", 10.0, "higher bin value in Vz fine histograms"};
    Configurable<int> CentfineNbins{"CentfineNbins", 16, "Number of bins in cent fine histograms"};
    Configurable<float> lfinebinCent{"lfinebinCent", 0.0, "lower bin value in cent fine histograms"};
    Configurable<float> hfinebinCent{"hfinebinCent", 80.0, "higher bin value in cent fine histograms"};
  } configbins;

  Configurable<bool> useShift{"useShift", false, "shift histograms"};
  Configurable<bool> ispolarization{"ispolarization", false, "Flag to check polarization"};
  Configurable<bool> followpub{"followpub", true, "flag to use alphaZDC"};
  Configurable<bool> useGainCallib{"useGainCallib", false, "use gain calibration"};
  Configurable<bool> useCallibvertex{"useCallibvertex", false, "use calibration for vxy"};
  Configurable<bool> coarse1{"coarse1", false, "RE1"};
  Configurable<bool> fine1{"fine1", false, "REfine1"};
  Configurable<bool> coarse2{"coarse2", false, "RE2"};
  Configurable<bool> fine2{"fine2", false, "REfine2"};
  Configurable<bool> coarse3{"coarse3", false, "RE3"};
  Configurable<bool> fine3{"fine3", false, "REfine3"};
  Configurable<bool> coarse4{"coarse4", false, "RE4"};
  Configurable<bool> fine4{"fine4", false, "REfine4"};
  Configurable<bool> coarse5{"coarse5", false, "RE5"};
  Configurable<bool> fine5{"fine5", false, "REfine5"};
  Configurable<bool> coarse6{"coarse6", false, "RE6"};
  Configurable<bool> fine6{"fine6", false, "REfine6"};
  Configurable<bool> useRecentereSp{"useRecentereSp", false, "use Recentering with Sparse or THn"};
  Configurable<bool> useRecenterefineSp{"useRecenterefineSp", false, "use fine Recentering with THn"};
  Configurable<std::string> ConfGainPath{"ConfGainPath", "Users/p/prottay/My/Object/NewPbPbpass4_10092024/gaincallib", "Path to gain calibration"};
  Configurable<std::string> ConfGainPathvxy{"ConfGainPathvxy", "Users/p/prottay/My/Object/swapcoords/PbPbpass4_20112024/recentervert", "Path to gain calibration for vxy"};
  Configurable<std::string> ConfRecentereSp{"ConfRecentereSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere"};
  Configurable<std::string> ConfRecentereSp2{"ConfRecentereSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere2"};
  Configurable<std::string> ConfRecentereSp3{"ConfRecentereSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere3"};
  Configurable<std::string> ConfRecentereSp4{"ConfRecentereSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere4"};
  Configurable<std::string> ConfRecentereSp5{"ConfRecentereSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere5"};
  Configurable<std::string> ConfRecentereSp6{"ConfRecentereSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere6"};
  Configurable<std::string> ConfRecenterecentSp{"ConfRecenterecentSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere"};
  Configurable<std::string> ConfRecenterevxSp{"ConfRecenterevxSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere"};
  Configurable<std::string> ConfRecenterevySp{"ConfRecenterevySp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere"};
  Configurable<std::string> ConfRecenterevzSp{"ConfRecenterevzSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere"};
  Configurable<std::string> ConfRecenterecentSp2{"ConfRecenterecentSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere2"};
  Configurable<std::string> ConfRecenterevxSp2{"ConfRecenterevxSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere2"};
  Configurable<std::string> ConfRecenterevySp2{"ConfRecenterevySp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere2"};
  Configurable<std::string> ConfRecenterevzSp2{"ConfRecenterevzSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere2"};
  Configurable<std::string> ConfRecenterecentSp3{"ConfRecenterecentSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere3"};
  Configurable<std::string> ConfRecenterevxSp3{"ConfRecenterevxSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere3"};
  Configurable<std::string> ConfRecenterevySp3{"ConfRecenterevySp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere3"};
  Configurable<std::string> ConfRecenterevzSp3{"ConfRecenterevzSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere3"};
  Configurable<std::string> ConfRecenterecentSp4{"ConfRecenterecentSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere4"};
  Configurable<std::string> ConfRecenterevxSp4{"ConfRecenterevxSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere4"};
  Configurable<std::string> ConfRecenterevySp4{"ConfRecenterevySp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere4"};
  Configurable<std::string> ConfRecenterevzSp4{"ConfRecenterevzSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere4"};
  Configurable<std::string> ConfRecenterecentSp5{"ConfRecenterecentSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere5"};
  Configurable<std::string> ConfRecenterevxSp5{"ConfRecenterevxSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere5"};
  Configurable<std::string> ConfRecenterevySp5{"ConfRecenterevySp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere5"};
  Configurable<std::string> ConfRecenterevzSp5{"ConfRecenterevzSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere5"};
  Configurable<std::string> ConfRecenterecentSp6{"ConfRecenterecentSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere6"};
  Configurable<std::string> ConfRecenterevxSp6{"ConfRecenterevxSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere6"};
  Configurable<std::string> ConfRecenterevySp6{"ConfRecenterevySp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere6"};
  Configurable<std::string> ConfRecenterevzSp6{"ConfRecenterevzSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere6"};
  Configurable<std::string> ConfShiftC{"ConfShiftC", "Users/p/prottay/My/Object/Testinglocaltree/shiftcallib2", "Path to shift C"};
  Configurable<std::string> ConfShiftA{"ConfShiftA", "Users/p/prottay/My/Object/Testinglocaltree/shiftcallib2", "Path to shift A"};

  // Event selection cuts - Alex
  /*
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  */
  /*
  template <typename TCollision>
  bool eventSelected(TCollision collision, const double& centrality)
  {
     auto multNTracksPV = collision.multNTracksPV();
     if (multNTracksPV < fMultPVCutLow->Eval(centrality))
       return 0;
     if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
       return 0;

    return 1;
  }
  */

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", true, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  RCTFlagsChecker rctChecker;

  void init(o2::framework::InitContext&)
  {

    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    AxisSpec channelZDCAxis = {8, 0.0, 8.0, "ZDC tower"};
    AxisSpec qxZDCAxis = {configbins.QxyNbins, configbins.lbinQxy, configbins.hbinQxy, "Qx"};
    AxisSpec phiAxis = {configbins.PhiNbins, -6.28, 6.28, "phi"};
    AxisSpec vzAxis = {configbins.VzNbins, configbins.lbinVz, configbins.hbinVz, "vz"};
    AxisSpec vxAxis = {configbins.VxNbins, configbins.lbinVx, configbins.hbinVx, "vx"};
    AxisSpec vyAxis = {configbins.VyNbins, configbins.lbinVy, configbins.hbinVy, "vy"};
    AxisSpec centAxis = {configbins.CentNbins, configbins.lbinCent, configbins.hbinCent, "V0M (%)"};
    AxisSpec vzfineAxis = {configbins.VzfineNbins, configbins.lfinebinVz, configbins.hfinebinVz, "vzfine"};
    AxisSpec vxfineAxis = {configbins.VxfineNbins, configbins.lfinebinVx, configbins.hfinebinVx, "vxfine"};
    AxisSpec vyfineAxis = {configbins.VyfineNbins, configbins.lfinebinVy, configbins.hfinebinVy, "vyfine"};
    AxisSpec centfineAxis = {configbins.CentfineNbins, configbins.lfinebinCent, configbins.hfinebinCent, "V0M (%) fine"};
    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {2, 0, 2, "basis"};
    AxisSpec VxyAxis = {2, 0, 2, "Vxy"};

    histos.add("htpcnsigmapi", "htpcnsigmapi", kTH1F, {{50, -10, 10.0}});
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0, 10.0}});
    histos.add("hCentrality", "hCentrality", kTH1F, {{centfineAxis}});
    histos.add("Vz", "Vz", kTH1F, {vzfineAxis});
    histos.add("hpQxZDCAC", "hpQxZDCAC", kTProfile, {centfineAxis});
    histos.add("hpQyZDCAC", "hpQyZDCAC", kTProfile, {centfineAxis});
    histos.add("hpQxZDCAQyZDCC", "hpQxZDCAQyZDCC", kTProfile, {centfineAxis});
    histos.add("hpQxZDCCQyZDCA", "hpQxZDCCQyZDCA", kTProfile, {centfineAxis});

    if (!ispolarization) {
      histos.add("hnQxZDCA", "hnQxZDCA", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
      histos.add("hnQyZDCA", "hnQyZDCA", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
      histos.add("hnQxZDCC", "hnQxZDCC", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
      histos.add("hnQyZDCC", "hnQyZDCC", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});

      histos.add("hcentQxZDCA", "hcentQxZDCA", kTH2F, {{centfineAxis}, {qxZDCAxis}});
      histos.add("hcentQyZDCA", "hcentQyZDCA", kTH2F, {{centfineAxis}, {qxZDCAxis}});
      histos.add("hcentQxZDCC", "hcentQxZDCC", kTH2F, {{centfineAxis}, {qxZDCAxis}});
      histos.add("hcentQyZDCC", "hcentQyZDCC", kTH2F, {{centfineAxis}, {qxZDCAxis}});

      histos.add("hvxQxZDCA", "hvxQxZDCA", kTH2F, {{vxfineAxis}, {qxZDCAxis}});
      histos.add("hvxQyZDCA", "hvxQyZDCA", kTH2F, {{vxfineAxis}, {qxZDCAxis}});
      histos.add("hvxQxZDCC", "hvxQxZDCC", kTH2F, {{vxfineAxis}, {qxZDCAxis}});
      histos.add("hvxQyZDCC", "hvxQyZDCC", kTH2F, {{vxfineAxis}, {qxZDCAxis}});

      histos.add("hvyQxZDCA", "hvyQxZDCA", kTH2F, {{vyfineAxis}, {qxZDCAxis}});
      histos.add("hvyQyZDCA", "hvyQyZDCA", kTH2F, {{vyfineAxis}, {qxZDCAxis}});
      histos.add("hvyQxZDCC", "hvyQxZDCC", kTH2F, {{vyfineAxis}, {qxZDCAxis}});
      histos.add("hvyQyZDCC", "hvyQyZDCC", kTH2F, {{vyfineAxis}, {qxZDCAxis}});

      histos.add("hvzQxZDCA", "hvzQxZDCA", kTH2F, {{vzfineAxis}, {qxZDCAxis}});
      histos.add("hvzQyZDCA", "hvzQyZDCA", kTH2F, {{vzfineAxis}, {qxZDCAxis}});
      histos.add("hvzQxZDCC", "hvzQxZDCC", kTH2F, {{vzfineAxis}, {qxZDCAxis}});
      histos.add("hvzQyZDCC", "hvzQyZDCC", kTH2F, {{vzfineAxis}, {qxZDCAxis}});
    }

    histos.add("PsiZDCC", "PsiZDCC", kTH2F, {centfineAxis, phiAxis});
    histos.add("PsiZDCA", "PsiZDCA", kTH2F, {centfineAxis, phiAxis});
    histos.add("ZDCAmp", "ZDCAmp", kTProfile2D, {channelZDCAxis, vzfineAxis});
    histos.add("ZDCAmpCommon", "ZDCAmpCommon", kTProfile2D, {{2, 0.0, 2.0}, vzfineAxis});
    histos.add("ShiftZDCC", "ShiftZDCC", kTProfile3D, {centfineAxis, basisAxis, shiftAxis});
    histos.add("ShiftZDCA", "ShiftZDCA", kTProfile3D, {centfineAxis, basisAxis, shiftAxis});
    histos.add("hpCosPsiAPsiC", "hpCosPsiAPsiC", kTProfile, {centfineAxis});
    histos.add("hpSinPsiAPsiC", "hpSinPsiAPsiC", kTProfile, {centfineAxis});
    histos.add("AvgVxy", "AvgVxy", kTProfile, {VxyAxis});

    // Event selection cut additional - Alex
    /*
    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
    fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
    fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
    fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
    fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
    fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    */
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TH2D* gainprofile;
  TProfile* gainprofilevxy;
  /*THnF* hrecentereSp;
  TH2F* hrecenterecentSp;
  TH2F* hrecenterevxSp;
  TH2F* hrecenterevySp;
  TH2F* hrecenterevzSp;*/
  std::array<THnF*, 6> hrecentereSpA;     // Array of 6 histograms
  std::array<TH2F*, 6> hrecenterecentSpA; // Array of 5 histograms
  std::array<TH2F*, 6> hrecenterevxSpA;   // Array of 5 histograms
  std::array<TH2F*, 6> hrecenterevySpA;   // Array of 5 histograms
  std::array<TH2F*, 6> hrecenterevzSpA;   // Array of 5 histograms
  TProfile3D* shiftprofileA;
  TProfile3D* shiftprofileC;

  // Bool_t Correctcoarse(int64_t ts, Configurable<std::string>& ConfRecentereSpp, bool useRecentereSp, int currentRunNumber, int lastRunNumber, auto centrality, auto vx, auto vy, auto vz, auto& qxZDCA, auto& qyZDCA, auto& qxZDCC, auto& qyZDCC)
  //{
  Bool_t Correctcoarse(const THnF* hrecentereSp, auto centrality, auto vx, auto vy, auto vz, auto& qxZDCA, auto& qyZDCA, auto& qxZDCC, auto& qyZDCC)
  {

    /*
    if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
      hrecentereSp = ccdb->getForTimeStamp<THnF>(ConfRecentereSpp.value, ts);
      }*/

    int binCoords[5];

    // Get axes of the THnSparse
    TAxis* centralityAxis = hrecentereSp->GetAxis(0); // Axis 0: centrality
    TAxis* vxAxis = hrecentereSp->GetAxis(1);         // Axis 1: vx
    TAxis* vyAxis = hrecentereSp->GetAxis(2);         // Axis 2: vy
    TAxis* vzAxis = hrecentereSp->GetAxis(3);         // Axis 3: vz
    TAxis* channelAxis = hrecentereSp->GetAxis(4);    // Axis 4: channel

    // Find bin indices for centrality, vx, vy, vz, and channel (for meanxA, 0.5)
    binCoords[0] = centralityAxis->FindBin(centrality + 0.00001); // Centrality
    binCoords[1] = vxAxis->FindBin(vx + 0.00001);                 // vx
    binCoords[2] = vyAxis->FindBin(vy + 0.00001);                 // vy
    binCoords[3] = vzAxis->FindBin(vz + 0.00001);                 // vz
    binCoords[4] = channelAxis->FindBin(0.5);                     // Channel for meanxA

    // Get the global bin for meanxA
    int globalBinMeanxA = hrecentereSp->GetBin(binCoords);
    float meanxA = hrecentereSp->GetBinContent(globalBinMeanxA);

    // Repeat for other channels (meanyA, meanxC, meanyC)
    binCoords[4] = channelAxis->FindBin(1.5); // Channel for meanyA
    int globalBinMeanyA = hrecentereSp->GetBin(binCoords);
    float meanyA = hrecentereSp->GetBinContent(globalBinMeanyA);

    binCoords[4] = channelAxis->FindBin(2.5); // Channel for meanxC
    int globalBinMeanxC = hrecentereSp->GetBin(binCoords);
    float meanxC = hrecentereSp->GetBinContent(globalBinMeanxC);

    binCoords[4] = channelAxis->FindBin(3.5); // Channel for meanyC
    int globalBinMeanyC = hrecentereSp->GetBin(binCoords);
    float meanyC = hrecentereSp->GetBinContent(globalBinMeanyC);

    qxZDCA = qxZDCA - meanxA;
    qyZDCA = qyZDCA - meanyA;
    qxZDCC = qxZDCC - meanxC;
    qyZDCC = qyZDCC - meanyC;

    return kTRUE;
  }

  // Bool_t Correctfine(int64_t ts, Configurable<std::string>& ConfRecenterecentSpp, Configurable<std::string>& ConfRecenterevxSpp, Configurable<std::string>& ConfRecenterevySpp, Configurable<std::string>& ConfRecenterevzSpp, bool useRecenterefineSp, int currentRunNumber, int lastRunNumber, auto centrality, auto vx, auto vy, auto vz, auto& qxZDCA, auto& qyZDCA, auto& qxZDCC, auto& qyZDCC)
  //{
  Bool_t Correctfine(TH2F* hrecenterecentSp, TH2F* hrecenterevxSp, TH2F* hrecenterevySp, TH2F* hrecenterevzSp, auto centrality, auto vx, auto vy, auto vz, auto& qxZDCA, auto& qyZDCA, auto& qxZDCC, auto& qyZDCC)
  {

    if (!hrecenterecentSp || !hrecenterevxSp || !hrecenterevySp || !hrecenterevzSp) {
      std::cerr << "Error: One or more histograms are null." << std::endl;
      return false;
    }
    /*
    if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
      hrecenterecentSp = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSpp.value, ts);
      hrecenterevxSp = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSpp.value, ts);
      hrecenterevySp = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySpp.value, ts);
      hrecenterevzSp = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSpp.value, ts);
      }*/

    double meanxAcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality + 0.00001, 0.5));
    double meanyAcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality + 0.00001, 1.5));
    double meanxCcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality + 0.00001, 2.5));
    double meanyCcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality + 0.00001, 3.5));

    double meanxAvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx + 0.00001, 0.5));
    double meanyAvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx + 0.00001, 1.5));
    double meanxCvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx + 0.00001, 2.5));
    double meanyCvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx + 0.00001, 3.5));

    double meanxAvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy + 0.00001, 0.5));
    double meanyAvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy + 0.00001, 1.5));
    double meanxCvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy + 0.00001, 2.5));
    double meanyCvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy + 0.00001, 3.5));

    double meanxAvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz + 0.00001, 0.5));
    double meanyAvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz + 0.00001, 1.5));
    double meanxCvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz + 0.00001, 2.5));
    double meanyCvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz + 0.00001, 3.5));

    qxZDCA = qxZDCA - meanxAcent - meanxAvx - meanxAvy - meanxAvz;
    qyZDCA = qyZDCA - meanyAcent - meanyAvx - meanyAvy - meanyAvz;
    qxZDCC = qxZDCC - meanxCcent - meanxCvx - meanxCvy - meanxCvz;
    qyZDCC = qyZDCC - meanyCcent - meanyCvx - meanyCvy - meanyCvz;

    return kTRUE;
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>;
  Preslice<aod::Zdcs> zdcPerCollision = aod::collision::bcId;

  void process(MyCollisions::iterator const& collision, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0s*/, BCsRun3 const& bcs, aod::Zdcs const&, AllTrackCandidates const& tracks)
  {

    if (usemem) {
      for (const auto& track : tracks) {
        histos.fill(HIST("htpcnsigmapi"), track.tpcNSigmaPi());
      }
    }

    histos.fill(HIST("hEvtSelInfo"), 0.5);
    auto centrality = collision.centFT0C();
    bool triggerevent = false;

    if (bcs.size() != 0) {
      gRandom->SetSeed(bcs.iteratorAt(0).globalBC());
    }

    currentRunNumber = collision.foundBC_as<BCsRun3>().runNumber();
    auto vz = collision.posZ();
    auto vx = collision.posX();
    auto vy = collision.posY();

    double psiZDCC = -99;
    double psiZDCA = -99;
    auto qxZDCA = 0.0;
    auto qxZDCC = 0.0;
    auto qyZDCA = 0.0;
    auto qyZDCC = 0.0;
    auto sumA = 0.0;
    auto sumC = 0.0;

    auto bc = collision.foundBC_as<BCsRun3>();

    if (!bc.has_zdc()) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 1.5);

    auto zdc = bc.zdc();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergy = zdc.energySectorZNA();
    auto zncEnergycommon = zdc.energyCommonZNC();
    auto znaEnergycommon = zdc.energyCommonZNA();

    if (znaEnergycommon <= 0.0 || zncEnergycommon <= 0.0) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 2.5);

    if (znaEnergy[0] <= 0.0 || znaEnergy[1] <= 0.0 || znaEnergy[2] <= 0.0 || znaEnergy[3] <= 0.0) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);

    if (zncEnergy[0] <= 0.0 || zncEnergy[1] <= 0.0 || zncEnergy[2] <= 0.0 || zncEnergy[3] <= 0.0) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 4.5);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 5.5);

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 6.5);

    if (collision.sel8() && centrality > cfgCutCentralityMin && centrality < cfgCutCentralityMax && TMath::Abs(vz) < cfgCutVertex && collision.has_foundFT0() && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TH2D>(ConfGainPath.value, bc.timestamp());
      }

      histos.fill(HIST("hEvtSelInfo"), 7.5);

      auto gainequal = 1.0;
      auto alphaZDC = 0.395;
      constexpr double x[4] = {-1.75, 1.75, -1.75, 1.75};
      constexpr double y[4] = {-1.75, -1.75, 1.75, 1.75};

      histos.fill(HIST("ZDCAmpCommon"), 0.5, vz, znaEnergycommon);
      histos.fill(HIST("ZDCAmpCommon"), 1.5, vz, zncEnergycommon);

      for (std::size_t iChA = 0; iChA < 8; iChA++) {
        auto chanelid = iChA;
        if (useGainCallib && gainprofile) {
          gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz + 0.00001, chanelid + 0.5));
        }

        if (iChA < 4) {

          if (znaEnergy[iChA] <= 0.0) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          } else {
            double ampl = gainequal * znaEnergy[iChA];
            if (followpub) {
              ampl = TMath::Power(ampl, alphaZDC);
            }
            qxZDCA = qxZDCA - ampl * x[iChA];
            qyZDCA = qyZDCA + ampl * y[iChA];
            sumA = sumA + ampl;
            histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
          }
        } else {
          if (zncEnergy[iChA - 4] <= 0.0) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          } else {
            double ampl = gainequal * zncEnergy[iChA - 4];
            if (followpub) {
              ampl = TMath::Power(ampl, alphaZDC);
            }
            qxZDCC = qxZDCC + ampl * x[iChA - 4];
            qyZDCC = qyZDCC + ampl * y[iChA - 4];
            sumC = sumC + ampl;
            histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
          }
        }
      }

      if (sumA > 0) {
        qxZDCA = qxZDCA / sumA;
        qyZDCA = qyZDCA / sumA;
      }
      if (sumC > 0) {
        qxZDCC = qxZDCC / sumC;
        qyZDCC = qyZDCC / sumC;
      }

      if (sumA <= 1e-4 || sumC <= 1e-4) {
        qxZDCA = 0.0;
        qxZDCC = 0.0;
        qyZDCA = 0.0;
        qyZDCC = 0.0;
        triggerevent = false;
        spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
        return;
      }

      histos.fill(HIST("hEvtSelInfo"), 8.5);
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);

      histos.fill(HIST("AvgVxy"), 0.5, vx);
      histos.fill(HIST("AvgVxy"), 1.5, vy);

      if (useCallibvertex && (currentRunNumber != lastRunNumber)) {
        gainprofilevxy = ccdb->getForTimeStamp<TProfile>(ConfGainPathvxy.value, bc.timestamp());
      }

      if (useCallibvertex) {
        vx = vx - gainprofilevxy->GetBinContent(1);
        vy = vy - gainprofilevxy->GetBinContent(2);
      }

      Bool_t res = 0;
      Bool_t resfine = 0;
      Int_t check = 1;

      if (coarse1) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[0] = ccdb->getForTimeStamp<THnF>(ConfRecentereSp.value, bc.timestamp());
        }
        res = Correctcoarse(hrecentereSpA[0], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine1) {

        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[0] = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp.value, bc.timestamp());
          hrecenterevxSpA[0] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp.value, bc.timestamp());
          hrecenterevySpA[0] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp.value, bc.timestamp());
          hrecenterevzSpA[0] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp.value, bc.timestamp());
        }
        resfine = Correctfine(hrecenterecentSpA[0], hrecenterevxSpA[0], hrecenterevySpA[0], hrecenterevzSpA[0], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse2) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[1] = ccdb->getForTimeStamp<THnF>(ConfRecentereSp2.value, bc.timestamp());
        }
        res = Correctcoarse(hrecentereSpA[1], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine2) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[1] = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp2.value, bc.timestamp());
          hrecenterevxSpA[1] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp2.value, bc.timestamp());
          hrecenterevySpA[1] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp2.value, bc.timestamp());
          hrecenterevzSpA[1] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp2.value, bc.timestamp());
        }
        resfine = Correctfine(hrecenterecentSpA[1], hrecenterevxSpA[1], hrecenterevySpA[1], hrecenterevzSpA[1], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse3) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[2] = ccdb->getForTimeStamp<THnF>(ConfRecentereSp3.value, bc.timestamp());
        }
        res = Correctcoarse(hrecentereSpA[2], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine3) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[2] = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp3.value, bc.timestamp());
          hrecenterevxSpA[2] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp3.value, bc.timestamp());
          hrecenterevySpA[2] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp3.value, bc.timestamp());
          hrecenterevzSpA[2] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp3.value, bc.timestamp());
        }
        resfine = Correctfine(hrecenterecentSpA[2], hrecenterevxSpA[2], hrecenterevySpA[2], hrecenterevzSpA[2], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse4) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[3] = ccdb->getForTimeStamp<THnF>(ConfRecentereSp4.value, bc.timestamp());
        }
        res = Correctcoarse(hrecentereSpA[3], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine4) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[3] = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp4.value, bc.timestamp());
          hrecenterevxSpA[3] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp4.value, bc.timestamp());
          hrecenterevySpA[3] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp4.value, bc.timestamp());
          hrecenterevzSpA[3] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp4.value, bc.timestamp());
        }
        resfine = Correctfine(hrecenterecentSpA[3], hrecenterevxSpA[3], hrecenterevySpA[3], hrecenterevzSpA[3], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse5) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[4] = ccdb->getForTimeStamp<THnF>(ConfRecentereSp5.value, bc.timestamp());
        }
        res = Correctcoarse(hrecentereSpA[4], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine5) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[4] = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp5.value, bc.timestamp());
          hrecenterevxSpA[4] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp5.value, bc.timestamp());
          hrecenterevySpA[4] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp5.value, bc.timestamp());
          hrecenterevzSpA[4] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp5.value, bc.timestamp());
        }
        resfine = Correctfine(hrecenterecentSpA[4], hrecenterevxSpA[4], hrecenterevySpA[4], hrecenterevzSpA[4], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse6) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[5] = ccdb->getForTimeStamp<THnF>(ConfRecentereSp6.value, bc.timestamp());
        }
        res = Correctcoarse(hrecentereSpA[5], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine6) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[5] = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp6.value, bc.timestamp());
          hrecenterevxSpA[5] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp6.value, bc.timestamp());
          hrecenterevySpA[5] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp6.value, bc.timestamp());
          hrecenterevzSpA[5] = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp6.value, bc.timestamp());
        }
        resfine = Correctfine(hrecenterecentSpA[5], hrecenterevxSpA[5], hrecenterevySpA[5], hrecenterevzSpA[5], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (res == 0 && resfine == 0 && check == 0) {
        LOG(info) << "Histograms are null";
      }
      psiZDCC = 1.0 * TMath::ATan2(qyZDCC, qxZDCC);
      psiZDCA = 1.0 * TMath::ATan2(qyZDCA, qxZDCA);

      int nshift = 10; // no. of iterations

      if (useShift && (currentRunNumber != lastRunNumber)) {
        shiftprofileC = ccdb->getForTimeStamp<TProfile3D>(ConfShiftC.value, bc.timestamp());
        shiftprofileA = ccdb->getForTimeStamp<TProfile3D>(ConfShiftA.value, bc.timestamp());
      }

      if (useShift) {
        auto deltapsiZDCC = 0.0;
        auto deltapsiZDCA = 0.0;
        for (int ishift = 1; ishift <= nshift; ishift++) {
          auto coeffshiftxZDCC = shiftprofileC->GetBinContent(shiftprofileC->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyZDCC = shiftprofileC->GetBinContent(shiftprofileC->FindBin(centrality, 1.5, ishift - 0.5));
          auto coeffshiftxZDCA = shiftprofileA->GetBinContent(shiftprofileA->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyZDCA = shiftprofileA->GetBinContent(shiftprofileA->FindBin(centrality, 1.5, ishift - 0.5));
          deltapsiZDCC = deltapsiZDCC + ((2 / (1.0 * ishift)) * (-coeffshiftxZDCC * TMath::Cos(ishift * 1.0 * psiZDCC) + coeffshiftyZDCC * TMath::Sin(ishift * 1.0 * psiZDCC)));
          deltapsiZDCA = deltapsiZDCA + ((2 / (1.0 * ishift)) * (-coeffshiftxZDCA * TMath::Cos(ishift * 1.0 * psiZDCA) + coeffshiftyZDCA * TMath::Sin(ishift * 1.0 * psiZDCA)));
        }
        psiZDCC = psiZDCC + deltapsiZDCC;
        psiZDCA = psiZDCA + deltapsiZDCA;
      }

      for (int ishift = 1; ishift <= nshift; ishift++) {
        histos.fill(HIST("ShiftZDCC"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * 1.0 * psiZDCC));
        histos.fill(HIST("ShiftZDCC"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * 1.0 * psiZDCC));
        histos.fill(HIST("ShiftZDCA"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * 1.0 * psiZDCA));
        histos.fill(HIST("ShiftZDCA"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * 1.0 * psiZDCA));
      }

      histos.fill(HIST("hpQxZDCAC"), centrality, (qxZDCA * qxZDCC));
      histos.fill(HIST("hpQyZDCAC"), centrality, (qyZDCA * qyZDCC));
      histos.fill(HIST("hpQxZDCAQyZDCC"), centrality, (qxZDCA * qyZDCC));
      histos.fill(HIST("hpQxZDCCQyZDCA"), centrality, (qxZDCC * qyZDCA));

      if (!ispolarization) {
        histos.fill(HIST("hnQxZDCA"), centrality, vx, vy, vz, qxZDCA);
        histos.fill(HIST("hnQyZDCA"), centrality, vx, vy, vz, qyZDCA);
        histos.fill(HIST("hnQxZDCC"), centrality, vx, vy, vz, qxZDCC);
        histos.fill(HIST("hnQyZDCC"), centrality, vx, vy, vz, qyZDCC);

        histos.fill(HIST("hcentQxZDCA"), centrality, qxZDCA);
        histos.fill(HIST("hcentQyZDCA"), centrality, qyZDCA);
        histos.fill(HIST("hcentQxZDCC"), centrality, qxZDCC);
        histos.fill(HIST("hcentQyZDCC"), centrality, qyZDCC);

        histos.fill(HIST("hvxQxZDCA"), vx, qxZDCA);
        histos.fill(HIST("hvxQyZDCA"), vx, qyZDCA);
        histos.fill(HIST("hvxQxZDCC"), vx, qxZDCC);
        histos.fill(HIST("hvxQyZDCC"), vx, qyZDCC);

        histos.fill(HIST("hvyQxZDCA"), vy, qxZDCA);
        histos.fill(HIST("hvyQyZDCA"), vy, qyZDCA);
        histos.fill(HIST("hvyQxZDCC"), vy, qxZDCC);
        histos.fill(HIST("hvyQyZDCC"), vy, qyZDCC);

        histos.fill(HIST("hvzQxZDCA"), vz, qxZDCA);
        histos.fill(HIST("hvzQyZDCA"), vz, qyZDCA);
        histos.fill(HIST("hvzQxZDCC"), vz, qxZDCC);
        histos.fill(HIST("hvzQyZDCC"), vz, qyZDCC);
      }

      histos.fill(HIST("hpCosPsiAPsiC"), centrality, (TMath::Cos(psiZDCA - psiZDCC)));
      histos.fill(HIST("hpSinPsiAPsiC"), centrality, (TMath::Sin(psiZDCA - psiZDCC)));
      histos.fill(HIST("PsiZDCA"), centrality, psiZDCA);
      histos.fill(HIST("PsiZDCC"), centrality, psiZDCC);

      lastRunNumber = currentRunNumber;
    }
    spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<spvector>(cfgc)};
}
