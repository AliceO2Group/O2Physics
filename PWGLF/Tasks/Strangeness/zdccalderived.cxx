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

#include "PWGLF/DataModel/ZDCCalTables.h"

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

#include "Framework/Logger.h"
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
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::rctsel;

// using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

struct zdccalderived {

  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Centrality cut Max"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Centrality cut Min"};
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};

  struct : ConfigurableGroup {
    Configurable<int> qxyNbins{"qxyNbins", 100, "Number of bins in QxQy histograms"};
    Configurable<int> phiNbins{"phiNbins", 100, "Number of bins in phi histogram"};
    Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
    Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};
    Configurable<int> vxNbins{"vxNbins", 25, "Number of bins in Vx histograms"};
    Configurable<float> lbinVx{"lbinVx", -0.05, "lower bin value in Vx histograms"};
    Configurable<float> hbinVx{"hbinVx", 0.0, "higher bin value in Vx histograms"};
    Configurable<int> vyNbins{"vyNbins", 25, "Number of bins in Vy histograms"};
    Configurable<float> lbinVy{"lbinVy", -0.02, "lower bin value in Vy histograms"};
    Configurable<float> hbinVy{"hbinVy", 0.02, "higher bin value in Vy histograms"};
    Configurable<int> vzNbins{"vzNbins", 20, "Number of bins in Vz histograms"};
    Configurable<float> lbinVz{"lbinVz", -10.0, "lower bin value in Vz histograms"};
    Configurable<float> hbinVz{"hbinVz", 10.0, "higher bin value in Vz histograms"};
    Configurable<int> centNbins{"centNbins", 16, "Number of bins in cent histograms"};
    Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
    Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
    Configurable<int> vxFineNbins{"vxFineNbins", 25, "Number of bins in Vx fine histograms"};
    Configurable<float> lfinebinVx{"lfinebinVx", -0.05, "lower bin value in Vx fine histograms"};
    Configurable<float> hfinebinVx{"hfinebinVx", 0.0, "higher bin value in Vx fine histograms"};
    Configurable<int> vyFineNbins{"vyFineNbins", 25, "Number of bins in Vy fine histograms"};
    Configurable<float> lfinebinVy{"lfinebinVy", -0.02, "lower bin value in Vy fine histograms"};
    Configurable<float> hfinebinVy{"hfinebinVy", 0.02, "higher bin value in Vy fine histograms"};
    Configurable<int> vzFineNbins{"vzFineNbins", 20, "Number of bins in Vz fine histograms"};
    Configurable<float> lfinebinVz{"lfinebinVz", -10.0, "lower bin value in Vz fine histograms"};
    Configurable<float> hfinebinVz{"hfinebinVz", 10.0, "higher bin value in Vz fine histograms"};
    Configurable<int> centFineNbins{"centFineNbins", 16, "Number of bins in cent fine histograms"};
    Configurable<float> lfinebinCent{"lfinebinCent", 0.0, "lower bin value in cent fine histograms"};
    Configurable<float> hfinebinCent{"hfinebinCent", 80.0, "higher bin value in cent fine histograms"};
  } configbins;

  Configurable<bool> useShift{"useShift", false, "shift histograms"};
  Configurable<bool> ispolarization{"ispolarization", false, "Flag to check polarization"};
  Configurable<bool> followpub{"followpub", true, "flag to use alphaZDC"};
  // Configurable<bool> useGainCallib{"useGainCallib", false, "use gain calibration"};
  // Configurable<bool> useCallibvertex{"useCallibvertex", false, "use calibration for vxy"};
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

  Configurable<std::string> confRecentereSp{"confRecentereSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for recentering"};
  Configurable<std::string> confRecentereSp2{"confRecentereSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for recentering 2"};
  Configurable<std::string> confRecentereSp3{"confRecentereSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for recentering 3"};
  Configurable<std::string> confRecentereSp4{"confRecentereSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for recentering 4"};
  Configurable<std::string> confRecentereSp5{"confRecentereSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for recentering 5"};
  Configurable<std::string> confRecentereSp6{"confRecentereSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for recentering 6"};
  Configurable<std::string> confRecentereCentSp{"confRecentereCentSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for cent recentering"};
  Configurable<std::string> confRecentereVxSp{"confRecentereVxSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vx recentering"};
  Configurable<std::string> confRecentereVySp{"confRecentereVySp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vy recentering"};
  Configurable<std::string> confRecentereVzSp{"confRecentereVzSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vz recentering"};
  Configurable<std::string> confRecentereCentSp2{"confRecentereCentSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for cent recentering 2"};
  Configurable<std::string> confRecentereVxSp2{"confRecentereVxSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vx recentering 2"};
  Configurable<std::string> confRecentereVySp2{"confRecentereVySp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vy recentering 2"};
  Configurable<std::string> confRecentereVzSp2{"confRecentereVzSp2", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vz recentering 2"};
  Configurable<std::string> confRecentereCentSp3{"confRecentereCentSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for cent recentering 3"};
  Configurable<std::string> confRecentereVxSp3{"confRecentereVxSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vx recentering 3"};
  Configurable<std::string> confRecentereVySp3{"confRecentereVySp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vy recentering 3"};
  Configurable<std::string> confRecentereVzSp3{"confRecentereVzSp3", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vz recentering 3"};
  Configurable<std::string> confRecentereCentSp4{"confRecentereCentSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for cent recentering 4"};
  Configurable<std::string> confRecentereVxSp4{"confRecentereVxSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vx recentering 4"};
  Configurable<std::string> confRecentereVySp4{"confRecentereVySp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vy recentering 4"};
  Configurable<std::string> confRecentereVzSp4{"confRecentereVzSp4", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vz recentering 4"};
  Configurable<std::string> confRecentereCentSp5{"confRecentereCentSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for cent recentering 5"};
  Configurable<std::string> confRecentereVxSp5{"confRecentereVxSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vx recentering 5"};
  Configurable<std::string> confRecentereVySp5{"confRecentereVySp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vy recentering 5"};
  Configurable<std::string> confRecentereVzSp5{"confRecentereVzSp5", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vz recentering 5"};
  Configurable<std::string> confRecentereCentSp6{"confRecentereCentSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for cent recentering 6"};
  Configurable<std::string> confRecentereVxSp6{"confRecentereVxSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vx recentering 6"};
  Configurable<std::string> confRecentereVySp6{"confRecentereVySp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vy recentering 6"};
  Configurable<std::string> confRecentereVzSp6{"confRecentereVzSp6", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn path for vz recentering 6"};
  Configurable<std::string> confShiftC{"confShiftC", "Users/p/prottay/My/Object/Testinglocaltree/shiftcallib2", "Path to shift C"};
  Configurable<std::string> confShiftA{"confShiftA", "Users/p/prottay/My/Object/Testinglocaltree/shiftcallib2", "Path to shift A"};

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", true, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  // RCTFlagsChecker rctChecker;

  void init(o2::framework::InitContext&)
  {

    // rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    AxisSpec channelZDCAxis = {8, 0.0, 8.0, "ZDC tower"};
    AxisSpec qxZDCAxis = {configbins.qxyNbins, configbins.lbinQxy, configbins.hbinQxy, "Qx"};
    AxisSpec phiAxis = {configbins.phiNbins, -6.28, 6.28, "phi"};
    AxisSpec vzAxis = {configbins.vzNbins, configbins.lbinVz, configbins.hbinVz, "vz"};
    AxisSpec vxAxis = {configbins.vxNbins, configbins.lbinVx, configbins.hbinVx, "vx"};
    AxisSpec vyAxis = {configbins.vyNbins, configbins.lbinVy, configbins.hbinVy, "vy"};
    AxisSpec centAxis = {configbins.centNbins, configbins.lbinCent, configbins.hbinCent, "V0M (%)"};
    AxisSpec vzfineAxis = {configbins.vzFineNbins, configbins.lfinebinVz, configbins.hfinebinVz, "vzfine"};
    AxisSpec vxfineAxis = {configbins.vxFineNbins, configbins.lfinebinVx, configbins.hfinebinVx, "vxfine"};
    AxisSpec vyfineAxis = {configbins.vyFineNbins, configbins.lfinebinVy, configbins.hfinebinVy, "vyfine"};
    AxisSpec centfineAxis = {configbins.centFineNbins, configbins.lfinebinCent, configbins.hfinebinCent, "V0M (%) fine"};
    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {2, 0, 2, "basis"};
    AxisSpec VxyAxis = {2, 0, 2, "Vxy"};

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
  std::array<THnF*, 6> hrecentereSpA;     // Array of 6 histograms
  std::array<TH2F*, 6> hrecenterecentSpA; // Array of 5 histograms
  std::array<TH2F*, 6> hrecenterevxSpA;   // Array of 5 histograms
  std::array<TH2F*, 6> hrecenterevySpA;   // Array of 5 histograms
  std::array<TH2F*, 6> hrecenterevzSpA;   // Array of 5 histograms
  TProfile3D* shiftprofileA;
  TProfile3D* shiftprofileC;

  bool Correctcoarse(const THnF* hrecentereSp, auto centrality, auto vx, auto vy, auto vz, auto& qxZDCA, auto& qyZDCA, auto& qxZDCC, auto& qyZDCC)
  {

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

  bool Correctfine(TH2F* hrecenterecentSp, TH2F* hrecenterevxSp, TH2F* hrecenterevySp, TH2F* hrecenterevzSp, auto centrality, auto vx, auto vy, auto vz, auto& qxZDCA, auto& qyZDCA, auto& qxZDCC, auto& qyZDCC)
  {

    if (!hrecenterecentSp || !hrecenterevxSp || !hrecenterevySp || !hrecenterevzSp) {
      // std::cerr << "Error: One or more histograms are null." << std::endl;
      return false;
    }

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

  using EventCandidates = aod::ZDCCalTables;
  void process(EventCandidates::iterator const& collision)
  {

    if (collision.triggerEventZDC()) {
      auto centrality = collision.cent();
      currentRunNumber = collision.triggerEventRunNo();

      auto vz = collision.vz();
      auto vx = collision.vx();
      auto vy = collision.vy();

      double psiZDCC = -99;
      double psiZDCA = -99;
      auto qxZDCA = collision.qxA();
      auto qxZDCC = collision.qxC();
      auto qyZDCA = collision.qyA();
      auto qyZDCC = collision.qyC();
      // auto sumA = 0.0;
      // auto sumC = 0.0;

      auto timestamps = ccdb->getRunDuration(currentRunNumber, true); /// fatalise if timestamps are not found
      int64_t sorTimestamp = timestamps.first;                        // timestamp of the SOR/SOX/STF in ms
      int64_t eorTimestamp = timestamps.second;                       // timestamp of the EOR/EOX/ETF in ms
      int64_t ts = eorTimestamp / 2 + sorTimestamp / 2;               // timestamp of the middle of the run
      /*
      std::array<float, 4> znaEnergy = {
        collision.znaE0(),
        collision.znaE1(),
        collision.znaE2(),
        collision.znaE3()};

      std::array<float, 4> zncEnergy = {
        collision.zncE0(),
        collision.zncE1(),
        collision.zncE2(),
        collision.zncE3()};

      auto znaEnergycommon = collision.znaC();
      auto zncEnergycommon = collision.zncC();

      histos.fill(HIST("hEvtSelInfo"), 7.5);

      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TH2D>(confGainPath.value, ts);
      }

      auto gainequal = 1.0;
      auto alphaZDC = 0.395;
      constexpr double x[4] = {-1.75, 1.75, -1.75, 1.75};
      constexpr double y[4] = {-1.75, -1.75, 1.75, 1.75};

      histos.fill(HIST("ZDCAmpCommon"), 0.5, vz, znaEnergycommon);
      histos.fill(HIST("ZDCAmpCommon"), 1.5, vz, zncEnergycommon);

      constexpr std::size_t ntow = 8;
      for (std::size_t iChA = 0; iChA < ntow; iChA++) {
        auto chanelid = iChA;
        if (useGainCallib && gainprofile) {
          gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz + 0.00001, chanelid + 0.5));
        }

        if (iChA < 4) {

          double ampl = gainequal * znaEnergy[iChA];
          if (followpub) {
            ampl = std::pow(ampl, alphaZDC);
          }
          qxZDCA = qxZDCA - ampl * x[iChA];
          qyZDCA = qyZDCA + ampl * y[iChA];
          sumA = sumA + ampl;
          histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
        } else {
          double ampl = gainequal * zncEnergy[iChA - 4];
          if (followpub) {
            ampl = std::pow(ampl, alphaZDC);
          }
          qxZDCC = qxZDCC + ampl * x[iChA - 4];
          qyZDCC = qyZDCC + ampl * y[iChA - 4];
          sumC = sumC + ampl;
          histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
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
        return;
      }
      */
      histos.fill(HIST("hEvtSelInfo"), 8.5);
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);

      histos.fill(HIST("AvgVxy"), 0.5, vx);
      histos.fill(HIST("AvgVxy"), 1.5, vy);
      /*
      if (useCallibvertex && (currentRunNumber != lastRunNumber)) {
        gainprofilevxy = ccdb->getForTimeStamp<TProfile>(confGainPathVxy.value, ts);
      }

      if (useCallibvertex) {
        vx = vx - gainprofilevxy->GetBinContent(1);
        vy = vy - gainprofilevxy->GetBinContent(2);
      }
      */
      bool res = 0;
      bool resfine = 0;
      int check = 1;

      if (coarse1) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[0] = ccdb->getForTimeStamp<THnF>(confRecentereSp.value, ts);
        }
        res = Correctcoarse(hrecentereSpA[0], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine1) {

        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[0] = ccdb->getForTimeStamp<TH2F>(confRecentereCentSp.value, ts);
          hrecenterevxSpA[0] = ccdb->getForTimeStamp<TH2F>(confRecentereVxSp.value, ts);
          hrecenterevySpA[0] = ccdb->getForTimeStamp<TH2F>(confRecentereVySp.value, ts);
          hrecenterevzSpA[0] = ccdb->getForTimeStamp<TH2F>(confRecentereVzSp.value, ts);
        }
        resfine = Correctfine(hrecenterecentSpA[0], hrecenterevxSpA[0], hrecenterevySpA[0], hrecenterevzSpA[0], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse2) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[1] = ccdb->getForTimeStamp<THnF>(confRecentereSp2.value, ts);
        }
        res = Correctcoarse(hrecentereSpA[1], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine2) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[1] = ccdb->getForTimeStamp<TH2F>(confRecentereCentSp2.value, ts);
          hrecenterevxSpA[1] = ccdb->getForTimeStamp<TH2F>(confRecentereVxSp2.value, ts);
          hrecenterevySpA[1] = ccdb->getForTimeStamp<TH2F>(confRecentereVySp2.value, ts);
          hrecenterevzSpA[1] = ccdb->getForTimeStamp<TH2F>(confRecentereVzSp2.value, ts);
        }
        resfine = Correctfine(hrecenterecentSpA[1], hrecenterevxSpA[1], hrecenterevySpA[1], hrecenterevzSpA[1], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse3) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[2] = ccdb->getForTimeStamp<THnF>(confRecentereSp3.value, ts);
        }
        res = Correctcoarse(hrecentereSpA[2], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine3) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[2] = ccdb->getForTimeStamp<TH2F>(confRecentereCentSp3.value, ts);
          hrecenterevxSpA[2] = ccdb->getForTimeStamp<TH2F>(confRecentereVxSp3.value, ts);
          hrecenterevySpA[2] = ccdb->getForTimeStamp<TH2F>(confRecentereVySp3.value, ts);
          hrecenterevzSpA[2] = ccdb->getForTimeStamp<TH2F>(confRecentereVzSp3.value, ts);
        }
        resfine = Correctfine(hrecenterecentSpA[2], hrecenterevxSpA[2], hrecenterevySpA[2], hrecenterevzSpA[2], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse4) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[3] = ccdb->getForTimeStamp<THnF>(confRecentereSp4.value, ts);
        }
        res = Correctcoarse(hrecentereSpA[3], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine4) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[3] = ccdb->getForTimeStamp<TH2F>(confRecentereCentSp4.value, ts);
          hrecenterevxSpA[3] = ccdb->getForTimeStamp<TH2F>(confRecentereVxSp4.value, ts);
          hrecenterevySpA[3] = ccdb->getForTimeStamp<TH2F>(confRecentereVySp4.value, ts);
          hrecenterevzSpA[3] = ccdb->getForTimeStamp<TH2F>(confRecentereVzSp4.value, ts);
        }
        resfine = Correctfine(hrecenterecentSpA[3], hrecenterevxSpA[3], hrecenterevySpA[3], hrecenterevzSpA[3], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse5) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[4] = ccdb->getForTimeStamp<THnF>(confRecentereSp5.value, ts);
        }
        res = Correctcoarse(hrecentereSpA[4], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine5) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[4] = ccdb->getForTimeStamp<TH2F>(confRecentereCentSp5.value, ts);
          hrecenterevxSpA[4] = ccdb->getForTimeStamp<TH2F>(confRecentereVxSp5.value, ts);
          hrecenterevySpA[4] = ccdb->getForTimeStamp<TH2F>(confRecentereVySp5.value, ts);
          hrecenterevzSpA[4] = ccdb->getForTimeStamp<TH2F>(confRecentereVzSp5.value, ts);
        }
        resfine = Correctfine(hrecenterecentSpA[4], hrecenterevxSpA[4], hrecenterevySpA[4], hrecenterevzSpA[4], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (coarse6) {
        if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
          hrecentereSpA[5] = ccdb->getForTimeStamp<THnF>(confRecentereSp6.value, ts);
        }
        res = Correctcoarse(hrecentereSpA[5], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (fine6) {
        if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
          hrecenterecentSpA[5] = ccdb->getForTimeStamp<TH2F>(confRecentereCentSp6.value, ts);
          hrecenterevxSpA[5] = ccdb->getForTimeStamp<TH2F>(confRecentereVxSp6.value, ts);
          hrecenterevySpA[5] = ccdb->getForTimeStamp<TH2F>(confRecentereVySp6.value, ts);
          hrecenterevzSpA[5] = ccdb->getForTimeStamp<TH2F>(confRecentereVzSp6.value, ts);
        }
        resfine = Correctfine(hrecenterecentSpA[5], hrecenterevxSpA[5], hrecenterevySpA[5], hrecenterevzSpA[5], centrality, vx, vy, vz, qxZDCA, qyZDCA, qxZDCC, qyZDCC);
      }

      if (res == 0 && resfine == 0 && check == 0) {
        LOG(info) << "Histograms are null";
      }
      psiZDCC = 1.0 * std::atan2(qyZDCC, qxZDCC);
      psiZDCA = 1.0 * std::atan2(qyZDCA, qxZDCA);

      int nshift = 10; // no. of iterations

      if (useShift && (currentRunNumber != lastRunNumber)) {
        shiftprofileC = ccdb->getForTimeStamp<TProfile3D>(confShiftC.value, ts);
        shiftprofileA = ccdb->getForTimeStamp<TProfile3D>(confShiftA.value, ts);
      }

      if (useShift) {
        auto deltapsiZDCC = 0.0;
        auto deltapsiZDCA = 0.0;
        for (int ishift = 1; ishift <= nshift; ishift++) {
          auto coeffshiftxZDCC = shiftprofileC->GetBinContent(shiftprofileC->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyZDCC = shiftprofileC->GetBinContent(shiftprofileC->FindBin(centrality, 1.5, ishift - 0.5));
          auto coeffshiftxZDCA = shiftprofileA->GetBinContent(shiftprofileA->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyZDCA = shiftprofileA->GetBinContent(shiftprofileA->FindBin(centrality, 1.5, ishift - 0.5));
          deltapsiZDCC = deltapsiZDCC + ((2 / (1.0 * ishift)) * (-coeffshiftxZDCC * std::cos(ishift * 1.0 * psiZDCC) + coeffshiftyZDCC * std::sin(ishift * 1.0 * psiZDCC)));
          deltapsiZDCA = deltapsiZDCA + ((2 / (1.0 * ishift)) * (-coeffshiftxZDCA * std::cos(ishift * 1.0 * psiZDCA) + coeffshiftyZDCA * std::sin(ishift * 1.0 * psiZDCA)));
        }
        psiZDCC = psiZDCC + deltapsiZDCC;
        psiZDCA = psiZDCA + deltapsiZDCA;
      }

      for (int ishift = 1; ishift <= nshift; ishift++) {
        histos.fill(HIST("ShiftZDCC"), centrality, 0.5, ishift - 0.5, std::sin(ishift * 1.0 * psiZDCC));
        histos.fill(HIST("ShiftZDCC"), centrality, 1.5, ishift - 0.5, std::cos(ishift * 1.0 * psiZDCC));
        histos.fill(HIST("ShiftZDCA"), centrality, 0.5, ishift - 0.5, std::sin(ishift * 1.0 * psiZDCA));
        histos.fill(HIST("ShiftZDCA"), centrality, 1.5, ishift - 0.5, std::cos(ishift * 1.0 * psiZDCA));
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

      histos.fill(HIST("hpCosPsiAPsiC"), centrality, (std::cos(psiZDCA - psiZDCC)));
      histos.fill(HIST("hpSinPsiAPsiC"), centrality, (std::sin(psiZDCA - psiZDCC)));
      histos.fill(HIST("PsiZDCA"), centrality, psiZDCA);
      histos.fill(HIST("PsiZDCC"), centrality, psiZDCC);

      lastRunNumber = currentRunNumber;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdccalderived>(cfgc)};
}
