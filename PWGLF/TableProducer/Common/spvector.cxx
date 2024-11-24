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
// \author: prottay das 07/09/2024
// \email: prottay.das@cern.ch

// C++/ROOT includes.
#include <TH1F.h>
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TMath.h>
#include <array>
#include <cmath>
#include "Math/Vector4D.h"
#include "TRandom3.h"
#include "TF1.h"

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"
// #include "SPCalibrationTableswrite.h"

// o2 includes.
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

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
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutPTMax{"cfgCutPTMax", 3.0, "Max PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgMinEta{"cfgMinEta", 0.1, "Min Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};

  Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
  Configurable<int> PhiNbins{"PhiNbins", 100, "Number of bins in phi histogram"};
  Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
  Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};
  // Configurable<int> ZDCgainNbins{"ZDCgainNbins", 500, "Number of bins in Gaineq histograms"};
  // Configurable<float> lbinZDCgain{"lbinZDCgain", 0.0, "lower bin value in Gaineq histograms"};
  // Configurable<float> hbinZDCgain{"hbinZDCgain", 1000.0, "higher bin value in Gaineq histograms"};
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
  Configurable<bool> QA{"QA", false, "QA histograms"};
  Configurable<bool> ispolarization{"ispolarization", false, "Flag to check polarization"};
  Configurable<bool> finecorrection{"finecorrection", false, "Flag to check fine correction"};
  Configurable<bool> rejbadevent{"rejbadevent", true, "Flag to check bad events"};
  Configurable<bool> rejbadeventcent{"rejbadeventcent", true, "Flag to check bad events for cent"};
  Configurable<bool> rejbadeventvx{"rejbadeventvx", true, "Flag to check bad events for vx"};
  Configurable<bool> rejbadeventvy{"rejbadeventvy", true, "Flag to check bad events for vy"};
  Configurable<bool> rejbadeventvz{"rejbadeventvz", true, "Flag to check bad events for vz"};
  Configurable<bool> usesparse{"usesparse", false, "flag to use sparse histogram"};
  Configurable<bool> usenormqn{"usenormqn", true, "flag to use normalized qs"};
  Configurable<bool> refsys{"refsys", true, "flag to use own reference system"};
  Configurable<bool> followpub{"followpub", true, "flag to use alphaZDC"};
  Configurable<bool> useGainCallib{"useGainCallib", false, "use gain calibration"};
  Configurable<bool> useRecentereSp{"useRecentereSp", false, "use Recentering with Sparse or THn"};
  Configurable<bool> useRecenterefineSp{"useRecenterefineSp", false, "use fine Recentering with Sparse or THn"};
  Configurable<bool> useRecenteresqSp{"useRecenteresqSp", false, "use Recenteringsq with Sparse or THn"};
  Configurable<bool> recwitherror{"recwitherror", false, "use Recentering with error"};
  Configurable<bool> recfinewitherror{"recfinewitherror", false, "use Recentering fine with error"};
  Configurable<std::string> ConfGainPath{"ConfGainPath", "Users/p/prottay/My/Object/NewPbPbpass4_10092024/gaincallib", "Path to gain calibration"};
  Configurable<std::string> ConfRecentereSp{"ConfRecentereSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere"};
  Configurable<std::string> ConfRecenterecentSp{"ConfRecenterecentSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for cent recentere"};
  Configurable<std::string> ConfRecenterevxSp{"ConfRecenterevxSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vx recentere"};
  Configurable<std::string> ConfRecenterevySp{"ConfRecenterevySp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vy recentere"};
  Configurable<std::string> ConfRecenterevzSp{"ConfRecenterevzSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for vz recentere"};
  Configurable<std::string> ConfRecenteresqSp{"ConfRecenteresqSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recenteresq"};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  int mRunNumber{-1};

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

  void initCCDB(BCsRun3::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
  }

  void init(o2::framework::InitContext&)
  {

    // const AxisSpec centAxis{configAxisCentrality, "V0M (%)"};
    //  AxisSpec amplitudeZDC = {configAxisZDCgain, "ZDC amplitude"};
    // AxisSpec amplitudeZDC = {ZDCgainNbins, lbinZDCgain, hbinZDCgain, "ZDC amplitude"};
    AxisSpec channelZDCAxis = {8, 0.0, 8.0, "ZDC tower"};
    AxisSpec qxZDCAxis = {QxyNbins, lbinQxy, hbinQxy, "Qx"};
    AxisSpec phiAxis = {PhiNbins, -6.28, 6.28, "phi"};
    AxisSpec vzAxis = {VzNbins, lbinVz, hbinVz, "vz"};
    AxisSpec vxAxis = {VxNbins, lbinVx, hbinVx, "vx"};
    AxisSpec vyAxis = {VyNbins, lbinVy, hbinVy, "vy"};
    AxisSpec centAxis = {CentNbins, lbinCent, hbinCent, "V0M (%)"};
    AxisSpec vzfineAxis = {VzfineNbins, lfinebinVz, hfinebinVz, "vzfine"};
    AxisSpec vxfineAxis = {VxfineNbins, lfinebinVx, hfinebinVx, "vxfine"};
    AxisSpec vyfineAxis = {VyfineNbins, lfinebinVy, hfinebinVy, "vyfine"};
    AxisSpec centfineAxis = {CentfineNbins, lfinebinCent, hfinebinCent, "V0M (%) fine"};

    histos.add("hCentrality", "hCentrality", kTH1F, {{centfineAxis}});
    histos.add("hpQxZDCAC", "hpQxZDCAC", kTProfile, {centfineAxis});
    histos.add("hpQyZDCAC", "hpQyZDCAC", kTProfile, {centfineAxis});
    histos.add("hpQxZDCAQyZDCC", "hpQxZDCAQyZDCC", kTProfile, {centfineAxis});
    histos.add("hpQxZDCCQyZDCA", "hpQxZDCCQyZDCA", kTProfile, {centfineAxis});

    if (!ispolarization) {
      if (usesparse == 1) {
        histos.add("hsQxZDCA", "hsQxZDCA", kTHnSparseF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
        histos.add("hsQyZDCA", "hsQyZDCA", kTHnSparseF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
        histos.add("hsQxZDCC", "hsQxZDCC", kTHnSparseF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
        histos.add("hsQyZDCC", "hsQyZDCC", kTHnSparseF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
      } else {
        histos.add("hnQxZDCA", "hnQxZDCA", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
        histos.add("hnQyZDCA", "hnQyZDCA", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
        histos.add("hnQxZDCC", "hnQxZDCC", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});
        histos.add("hnQyZDCC", "hnQyZDCC", kTHnF, {{centAxis}, {vxAxis}, {vyAxis}, {vzAxis}, {qxZDCAxis}});

        if (finecorrection) {
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
      }
    }

    histos.add("PsiZDCC", "PsiZDCC", kTH2F, {centfineAxis, phiAxis});
    histos.add("PsiZDCA", "PsiZDCA", kTH2F, {centfineAxis, phiAxis});
    // histos.add("ZDCAmp", "ZDCAmp", kTProfile3D, {channelZDCAxis, vzfineAxis, centfineAxis});
    histos.add("ZDCAmp", "ZDCAmp", kTProfile2D, {channelZDCAxis, vzfineAxis});
    histos.add("ZDCAmpCommon", "ZDCAmpCommon", kTProfile2D, {{2, 0.0, 2.0}, vzfineAxis});
    // histos.add("ZDCAmpCommon", "ZDCAmpCommon", kTProfile3D, {{2,0.0,2.0}, vzfineAxis, centfineAxis});

    if (QA) {
      histos.add("Vz", "Vz", kTH1F, {vzfineAxis});
      histos.add("hpCosPsiAPsiC", "hpCosPsiAPsiC", kTProfile, {centfineAxis});
      histos.add("hpSinPsiAPsiC", "hpSinPsiAPsiC", kTProfile, {centfineAxis});
    }
    // histos.add("hZDCAmp", "hZDCAmp", kTHnF, {channelZDCAxis, vzAxis, centfineAxis, {1000, 0, 1000}});

    // Event selection cut additional - Alex
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

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  // TH3D* gainprofile;
  TH2D* gainprofile;
  THnF* hrecentereSp;
  TH2F* hrecenterecentSp;
  TH2F* hrecenterevxSp;
  TH2F* hrecenterevySp;
  TH2F* hrecenterevzSp;
  THnF* hrecenteresqSp;

  // Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  // Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  // Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs>;
  // using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>>;

  Preslice<aod::Zdcs> zdcPerCollision = aod::collision::bcId;

  void process(MyCollisions::iterator const& collision, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0s*/, BCsRun3 const& bcs, aod::Zdcs const&)
  {

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

    if (znaEnergy[0] <= 0.0 || znaEnergy[1] <= 0.0 || znaEnergy[2] <= 0.0 || znaEnergy[3] <= 0.0) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }
    if (zncEnergy[0] <= 0.0 || zncEnergy[1] <= 0.0 || zncEnergy[2] <= 0.0 || zncEnergy[3] <= 0.0) {
      triggerevent = false;
      spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
      return;
    }

    if (collision.sel8() && centrality > cfgCutCentralityMin && centrality < cfgCutCentralityMax && TMath::Abs(vz) < cfgCutVertex && collision.has_foundFT0() && eventSelected(collision, centrality) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        // gainprofile = ccdb->getForTimeStamp<TH3D>(ConfGainPath.value, bc.timestamp());
        gainprofile = ccdb->getForTimeStamp<TH2D>(ConfGainPath.value, bc.timestamp());
      }

      // initCCDB(bc);

      auto gainequal = 1.0;
      auto alphaZDC = 0.395;
      constexpr double x[4] = {-1.75, 1.75, -1.75, 1.75};
      constexpr double y[4] = {-1.75, -1.75, 1.75, 1.75};

      // histos.fill(HIST("ZDCAmpCommon"), 0.5, vz, centrality, znaEnergycommon);
      // histos.fill(HIST("ZDCAmpCommon"), 1.5, vz, centrality, zncEnergycommon);
      histos.fill(HIST("ZDCAmpCommon"), 0.5, vz, znaEnergycommon);
      histos.fill(HIST("ZDCAmpCommon"), 1.5, vz, zncEnergycommon);

      // LOG(info) << "**********energy values************" << znaEnergycommon<<" "<<znaEnergy[0]<<" "<<znaEnergy[1]<<" "<<znaEnergy[2]<<" "<<znaEnergy[3]<<" "<<znaEnergy[0]+znaEnergy[1]+znaEnergy[2]+znaEnergy[3];

      for (std::size_t iChA = 0; iChA < 8; iChA++) {
        auto chanelid = iChA;
        if (useGainCallib && gainprofile) {
          // gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz, centrality, chanelid + 0.5));
          gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz, chanelid + 0.5));
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
            // histos.fill(HIST("hZDCAmp"), chanelid + 0.5, vz, centrality, ampl);
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
            // histos.fill(HIST("hZDCAmp"), chanelid + 0.5, vz, centrality, ampl);
          }
        }
      }

      if (usenormqn) {
        if (sumA > 0) {
          qxZDCA = qxZDCA / sumA;
          qyZDCA = qyZDCA / sumA;
        }
        if (sumC > 0) {
          qxZDCC = qxZDCC / sumC;
          qyZDCC = qyZDCC / sumC;
        }
      } else {
        qxZDCA = qxZDCA;
        qxZDCC = qxZDCC;
        qyZDCA = qyZDCA;
        qyZDCC = qyZDCC;
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

      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);

      if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
        hrecentereSp = ccdb->getForTimeStamp<THnF>(ConfRecentereSp.value, bc.timestamp());
      }
      if (useRecenterefineSp && (currentRunNumber != lastRunNumber)) {
        hrecenterecentSp = ccdb->getForTimeStamp<TH2F>(ConfRecenterecentSp.value, bc.timestamp());
        hrecenterevxSp = ccdb->getForTimeStamp<TH2F>(ConfRecenterevxSp.value, bc.timestamp());
        hrecenterevySp = ccdb->getForTimeStamp<TH2F>(ConfRecenterevySp.value, bc.timestamp());
        hrecenterevzSp = ccdb->getForTimeStamp<TH2F>(ConfRecenterevzSp.value, bc.timestamp());
      }
      if (useRecenteresqSp && (currentRunNumber != lastRunNumber)) {
        hrecenteresqSp = ccdb->getForTimeStamp<THnF>(ConfRecenteresqSp.value, bc.timestamp());
      }

      if (useRecentereSp && hrecentereSp) {

        int binCoords[5];

        // Get axes of the THnSparse
        TAxis* centralityAxis = hrecentereSp->GetAxis(0); // Axis 0: centrality
        TAxis* vxAxis = hrecentereSp->GetAxis(1);         // Axis 1: vx
        TAxis* vyAxis = hrecentereSp->GetAxis(2);         // Axis 2: vy
        TAxis* vzAxis = hrecentereSp->GetAxis(3);         // Axis 3: vz
        TAxis* channelAxis = hrecentereSp->GetAxis(4);    // Axis 4: channel

        // Find bin indices for centrality, vx, vy, vz, and channel (for meanxA, 0.5)
        binCoords[0] = centralityAxis->FindBin(centrality); // Centrality
        binCoords[1] = vxAxis->FindBin(vx);                 // vx
        binCoords[2] = vyAxis->FindBin(vy);                 // vy
        binCoords[3] = vzAxis->FindBin(vz);                 // vz
        binCoords[4] = channelAxis->FindBin(0.5);           // Channel for meanxA

        // Get the global bin for meanxA
        int globalBinMeanxA = hrecentereSp->GetBin(binCoords);
        double meanxA = hrecentereSp->GetBinContent(globalBinMeanxA);
        double meanxAerror = hrecentereSp->GetBinError(globalBinMeanxA);

        // Repeat for other channels (meanyA, meanxC, meanyC)
        binCoords[4] = channelAxis->FindBin(1.5); // Channel for meanyA
        int globalBinMeanyA = hrecentereSp->GetBin(binCoords);
        double meanyA = hrecentereSp->GetBinContent(globalBinMeanyA);
        double meanyAerror = hrecentereSp->GetBinError(globalBinMeanyA);

        binCoords[4] = channelAxis->FindBin(2.5); // Channel for meanxC
        int globalBinMeanxC = hrecentereSp->GetBin(binCoords);
        double meanxC = hrecentereSp->GetBinContent(globalBinMeanxC);
        double meanxCerror = hrecentereSp->GetBinError(globalBinMeanxC);

        binCoords[4] = channelAxis->FindBin(3.5); // Channel for meanyC
        int globalBinMeanyC = hrecentereSp->GetBin(binCoords);
        double meanyC = hrecentereSp->GetBinContent(globalBinMeanyC);
        double meanyCerror = hrecentereSp->GetBinError(globalBinMeanyC);

        if (rejbadevent) {
          if ((TMath::Abs(meanxA) > 90000.0 || TMath::Abs(meanxC) > 90000.0 || TMath::Abs(meanyA) > 90000.0 || TMath::Abs(meanyC) > 90000.0) && (TMath::Abs(meanxAerror) > 9000.0 || TMath::Abs(meanxCerror) > 9000.0 || TMath::Abs(meanyAerror) > 9000.0 || TMath::Abs(meanyCerror) > 9000.0)) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          }
        }

        qxZDCA = qxZDCA - meanxA;
        qyZDCA = qyZDCA - meanyA;
        qxZDCC = qxZDCC - meanxC;
        qyZDCC = qyZDCC - meanyC;

        if (recwitherror) {
          if (meanxAerror != 0.0) {
            qxZDCA = qxZDCA / meanxAerror;
          }
          if (meanyAerror != 0.0) {
            qyZDCA = qyZDCA / meanyAerror;
          }
          if (meanxCerror != 0.0) {
            qxZDCC = qxZDCC / meanxCerror;
          }
          if (meanyCerror != 0.0) {
            qyZDCC = qyZDCC / meanyCerror;
          }
        }
      }

      if (useRecenterefineSp && hrecenterecentSp) {

        double meanxAcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 0.5));
        double meanyAcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 1.5));
        double meanxCcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 2.5));
        double meanyCcent = hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 3.5));

        double meanxAcenterror = hrecenterecentSp->GetBinError(hrecenterecentSp->FindBin(centrality, 0.5));
        double meanyAcenterror = hrecenterecentSp->GetBinError(hrecenterecentSp->FindBin(centrality, 1.5));
        double meanxCcenterror = hrecenterecentSp->GetBinError(hrecenterecentSp->FindBin(centrality, 2.5));
        double meanyCcenterror = hrecenterecentSp->GetBinError(hrecenterecentSp->FindBin(centrality, 3.5));

        if (rejbadeventcent) {
          if ((TMath::Abs(meanxAcent) > 90000.0 || TMath::Abs(meanxCcent) > 90000.0 || TMath::Abs(meanyAcent) > 90000.0 || TMath::Abs(meanyCcent) > 90000.0) && (TMath::Abs(meanxAcenterror) > 9000.0 || TMath::Abs(meanxCcenterror) > 9000.0 || TMath::Abs(meanyAcenterror) > 9000.0 || TMath::Abs(meanyCcenterror) > 9000.0)) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          }
        }

        qxZDCA = qxZDCA - hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 0.5));
        qyZDCA = qyZDCA - hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 1.5));
        qxZDCC = qxZDCC - hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 2.5));
        qyZDCC = qyZDCC - hrecenterecentSp->GetBinContent(hrecenterecentSp->FindBin(centrality, 3.5));

        if (recfinewitherror) {
          qxZDCA = qxZDCA / meanxAcenterror;
          qyZDCA = qyZDCA / meanyAcenterror;
          qxZDCC = qxZDCC / meanxCcenterror;
          qyZDCC = qyZDCC / meanyCcenterror;
        }
      }

      if (useRecenterefineSp && hrecenterevxSp) {

        double meanxAvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 0.5));
        double meanyAvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 1.5));
        double meanxCvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 2.5));
        double meanyCvx = hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 3.5));

        double meanxAvxerror = hrecenterevxSp->GetBinError(hrecenterevxSp->FindBin(vx, 0.5));
        double meanyAvxerror = hrecenterevxSp->GetBinError(hrecenterevxSp->FindBin(vx, 1.5));
        double meanxCvxerror = hrecenterevxSp->GetBinError(hrecenterevxSp->FindBin(vx, 2.5));
        double meanyCvxerror = hrecenterevxSp->GetBinError(hrecenterevxSp->FindBin(vx, 3.5));

        if (rejbadeventvx) {
          if ((TMath::Abs(meanxAvx) > 90000.0 || TMath::Abs(meanxCvx) > 90000.0 || TMath::Abs(meanyAvx) > 90000.0 || TMath::Abs(meanyCvx) > 90000.0) && (TMath::Abs(meanxAvxerror) > 9000.0 || TMath::Abs(meanxCvxerror) > 9000.0 || TMath::Abs(meanyAvxerror) > 9000.0 || TMath::Abs(meanyCvxerror) > 9000.0)) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          }
        }

        qxZDCA = qxZDCA - hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 0.5));
        qyZDCA = qyZDCA - hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 1.5));
        qxZDCC = qxZDCC - hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 2.5));
        qyZDCC = qyZDCC - hrecenterevxSp->GetBinContent(hrecenterevxSp->FindBin(vx, 3.5));

        if (recfinewitherror) {
          qxZDCA = qxZDCA / meanxAvxerror;
          qyZDCA = qyZDCA / meanyAvxerror;
          qxZDCC = qxZDCC / meanxCvxerror;
          qyZDCC = qyZDCC / meanyCvxerror;
        }
      }

      if (useRecenterefineSp && hrecenterevySp) {

        double meanxAvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 0.5));
        double meanyAvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 1.5));
        double meanxCvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 2.5));
        double meanyCvy = hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 3.5));

        double meanxAvyerror = hrecenterevySp->GetBinError(hrecenterevySp->FindBin(vy, 0.5));
        double meanyAvyerror = hrecenterevySp->GetBinError(hrecenterevySp->FindBin(vy, 1.5));
        double meanxCvyerror = hrecenterevySp->GetBinError(hrecenterevySp->FindBin(vy, 2.5));
        double meanyCvyerror = hrecenterevySp->GetBinError(hrecenterevySp->FindBin(vy, 3.5));

        if (rejbadeventvy) {
          if ((TMath::Abs(meanxAvy) > 90000.0 || TMath::Abs(meanxCvy) > 90000.0 || TMath::Abs(meanyAvy) > 90000.0 || TMath::Abs(meanyCvy) > 90000.0) && (TMath::Abs(meanxAvyerror) > 9000.0 || TMath::Abs(meanxCvyerror) > 9000.0 || TMath::Abs(meanyAvyerror) > 9000.0 || TMath::Abs(meanyCvyerror) > 9000.0)) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          }
        }

        qxZDCA = qxZDCA - hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 0.5));
        qyZDCA = qyZDCA - hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 1.5));
        qxZDCC = qxZDCC - hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 2.5));
        qyZDCC = qyZDCC - hrecenterevySp->GetBinContent(hrecenterevySp->FindBin(vy, 3.5));

        if (recfinewitherror) {
          qxZDCA = qxZDCA / meanxAvyerror;
          qyZDCA = qyZDCA / meanyAvyerror;
          qxZDCC = qxZDCC / meanxCvyerror;
          qyZDCC = qyZDCC / meanyCvyerror;
        }
      }

      if (useRecenterefineSp && hrecenterevzSp) {

        double meanxAvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 0.5));
        double meanyAvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 1.5));
        double meanxCvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 2.5));
        double meanyCvz = hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 3.5));

        double meanxAvzerror = hrecenterevzSp->GetBinError(hrecenterevzSp->FindBin(vz, 0.5));
        double meanyAvzerror = hrecenterevzSp->GetBinError(hrecenterevzSp->FindBin(vz, 1.5));
        double meanxCvzerror = hrecenterevzSp->GetBinError(hrecenterevzSp->FindBin(vz, 2.5));
        double meanyCvzerror = hrecenterevzSp->GetBinError(hrecenterevzSp->FindBin(vz, 3.5));

        if (rejbadeventvz) {
          if ((TMath::Abs(meanxAvz) > 90000.0 || TMath::Abs(meanxCvz) > 90000.0 || TMath::Abs(meanyAvz) > 90000.0 || TMath::Abs(meanyCvz) > 90000.0) && (TMath::Abs(meanxAvzerror) > 9000.0 || TMath::Abs(meanxCvzerror) > 9000.0 || TMath::Abs(meanyAvzerror) > 9000.0 || TMath::Abs(meanyCvzerror) > 9000.0)) {
            triggerevent = false;
            spcalibrationtable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qxZDCC, qyZDCA, qyZDCC, psiZDCC, psiZDCA);
            return;
          }
        }

        qxZDCA = qxZDCA - hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 0.5));
        qyZDCA = qyZDCA - hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 1.5));
        qxZDCC = qxZDCC - hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 2.5));
        qyZDCC = qyZDCC - hrecenterevzSp->GetBinContent(hrecenterevzSp->FindBin(vz, 3.5));

        if (recfinewitherror) {
          qxZDCA = qxZDCA / meanxAvzerror;
          qyZDCA = qyZDCA / meanyAvzerror;
          qxZDCC = qxZDCC / meanxCvzerror;
          qyZDCC = qyZDCC / meanyCvzerror;
        }
      }

      // LOG(info) << "**********qxa values in spvector************" << qxZDCA<<" "<<centrality<<" "<<vx<<" "<<vy<<" "<<vz;

      psiZDCC = 1.0 * TMath::ATan2(qyZDCC, qxZDCC);
      psiZDCA = 1.0 * TMath::ATan2(qyZDCA, qxZDCA);

      histos.fill(HIST("hpQxZDCAC"), centrality, (qxZDCA * qxZDCC));
      histos.fill(HIST("hpQyZDCAC"), centrality, (qyZDCA * qyZDCC));
      histos.fill(HIST("hpQxZDCAQyZDCC"), centrality, (qxZDCA * qyZDCC));
      histos.fill(HIST("hpQxZDCCQyZDCA"), centrality, (qxZDCC * qyZDCA));

      if (!ispolarization) {
        if (usesparse) {
          histos.fill(HIST("hsQxZDCA"), centrality, vx, vy, vz, qxZDCA);
          histos.fill(HIST("hsQyZDCA"), centrality, vx, vy, vz, qyZDCA);
          histos.fill(HIST("hsQxZDCC"), centrality, vx, vy, vz, qxZDCC);
          histos.fill(HIST("hsQyZDCC"), centrality, vx, vy, vz, qyZDCC);
        } else {
          histos.fill(HIST("hnQxZDCA"), centrality, vx, vy, vz, qxZDCA);
          histos.fill(HIST("hnQyZDCA"), centrality, vx, vy, vz, qyZDCA);
          histos.fill(HIST("hnQxZDCC"), centrality, vx, vy, vz, qxZDCC);
          histos.fill(HIST("hnQyZDCC"), centrality, vx, vy, vz, qyZDCC);

          if (finecorrection) {
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
        }
      }

      if (QA) {
        histos.fill(HIST("hpCosPsiAPsiC"), centrality, (TMath::Cos(psiZDCA - psiZDCC)));
        histos.fill(HIST("hpSinPsiAPsiC"), centrality, (TMath::Sin(psiZDCA - psiZDCC)));
      }
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
