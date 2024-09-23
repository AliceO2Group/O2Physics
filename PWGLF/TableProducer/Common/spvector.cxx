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
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Centrality cut"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutPTMax{"cfgCutPTMax", 3.0, "Max PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgMinEta{"cfgMinEta", 0.1, "Min Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};

  Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
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
  Configurable<int> CentNbins{"CentNbins", 16, "Number of bins in cent histograms"};
  Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
  Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
  Configurable<bool> QA{"QA", false, "QA histograms"};
  Configurable<bool> usesparse{"usesparse", false, "flag to use sparse histogram"};
  Configurable<bool> usenormqn{"usenormqn", true, "flag to use normalized qs"};
  Configurable<bool> tablewrite{"tablewrite", false, "Boolean for writing table"};
  Configurable<bool> useGainCallib{"useGainCallib", false, "use gain calibration"};
  // Configurable<bool> useRecentere{"useRecentere", false, "use Recentering"};
  Configurable<bool> useRecentereSp{"useRecentereSp", false, "use Recentering with Sparse or THn"};
  Configurable<bool> useRecenteresqSp{"useRecenteresqSp", false, "use Recenteringsq with Sparse or THn"};
  // Configurable<bool> useRecentereVxy{"useRecentereVxy", false, "use Recentering for Vxy"};
  Configurable<bool> recwitherror{"recwitherror", false, "use Recentering with error"};
  // Configurable<bool> useShift{"useShift", false, "use Shift"};
  Configurable<std::string> ConfGainPath{"ConfGainPath", "Users/p/prottay/My/Object/NewPbPbpass4_10092024/gaincallib", "Path to gain calibration"};
  // Configurable<std::string> ConfRecentere{"ConfRecentere", "Users/p/prottay/My/Object/NewPbPbpass4_23082024/recenter", "Path for recentere"};
  Configurable<std::string> ConfRecentereSp{"ConfRecentereSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recentere"};
  Configurable<std::string> ConfRecenteresqSp{"ConfRecenteresqSp", "Users/p/prottay/My/Object/Testingwithsparse/NewPbPbpass4_17092024/recenter", "Sparse or THn Path for recenteresq"};
  /*
  Configurable<std::string> ConfRecentereVxyQxA{"ConfRecentereVxyQxA", "Users/p/prottay/My/Object/NewPbPbpass4_23082024/recenter", "Path for recentereVxyQxA"};
  Configurable<std::string> ConfRecentereVxyQyA{"ConfRecentereVxyQyA", "Users/p/prottay/My/Object/NewPbPbpass4_23082024/recenter", "Path for recentereVxyQyA"};
  Configurable<std::string> ConfRecentereVxyQxC{"ConfRecentereVxyQxC", "Users/p/prottay/My/Object/NewPbPbpass4_23082024/recenter", "Path for recentereVxyQxC"};
  Configurable<std::string> ConfRecentereVxyQyC{"ConfRecentereVxyQyC", "Users/p/prottay/My/Object/NewPbPbpass4_23082024/recenter", "Path for recentereVxyQyC"};
  Configurable<std::string> ConfShift{"ConfShift", "Users/p/prottay/My/Object/Finaltest2/recenereall", "Path for Shift"};
  */
  // ConfigurableAxis configAxisCentrality{"configAxisCentrality", {16, 0.0, 80}, "centrality bining"};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  int mRunNumber{-1};

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
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
    AxisSpec phiAxis = {50, -6.28, 6.28, "phi"};
    AxisSpec vzAxis = {20, -10, 10, "vz"};
    AxisSpec vxAxis = {VxNbins, lbinVx, hbinVx, "vx"};
    AxisSpec vyAxis = {VyNbins, lbinVy, hbinVy, "vy"};
    AxisSpec centAxis = {CentNbins, lbinCent, hbinCent, "V0M (%)"};

    histos.add("hCentrality", "hCentrality", kTH1F, {{centAxis}});
    histos.add("hpQxZDCAC", "hpQxZDCAC", kTProfile, {centAxis});
    histos.add("hpQyZDCAC", "hpQyZDCAC", kTProfile, {centAxis});
    histos.add("hpQxZDCAQyZDCC", "hpQxZDCAQyZDCC", kTProfile, {centAxis});
    histos.add("hpQxZDCCQyZDCA", "hpQxZDCCQyZDCA", kTProfile, {centAxis});
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
    }

    /*histos.add("hpQxVZDCC", "hpQxVZDCC", kTProfile3D, {centAxis, vxAxis, vyAxis});
    histos.add("hpQyVZDCC", "hpQyVZDCC", kTProfile3D, {centAxis, vxAxis, vyAxis});
    histos.add("hpQxVZDCA", "hpQxVZDCA", kTProfile3D, {centAxis, vxAxis, vyAxis});
    histos.add("hpQyVZDCA", "hpQyVZDCA", kTProfile3D, {centAxis, vxAxis, vyAxis});*/
    histos.add("PsiZDCC", "PsiZDCC", kTH2F, {centAxis, phiAxis});
    histos.add("PsiZDCA", "PsiZDCA", kTH2F, {centAxis, phiAxis});
    histos.add("ZDCAmp", "ZDCAmp", kTProfile2D, {channelZDCAxis, vzAxis});

    if (QA) {
      histos.add("Vz", "Vz", kTH1F, {vzAxis});
      histos.add("hpCosPsiAPsiC", "hpCosPsiAPsiC", kTProfile, {centAxis});
      histos.add("hpSinPsiAPsiC", "hpSinPsiAPsiC", kTProfile, {centAxis});
    }
    // histos.add("hZDCAmp", "hZDCAmp", kTH3F, {channelZDCAxis, vzAxis, amplitudeZDC});

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
  TH2D* gainprofile;
  // TH3D* hrecentere;
  // THnSparseF* hrecentereSp;
  THnF* hrecentereSp;
  THnF* hrecenteresqSp;

  /*
  TH3D* hrecentereVxyQxA;
  TH3D* hrecentereVxyQyA;
  TH3D* hrecentereVxyQxC;
  TH3D* hrecentereVxyQyC;
  */

  // Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  // Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs>;
  // using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>>;

  Preslice<aod::Zdcs> zdcPerCollision = aod::collision::bcId;

  // void process(MyCollisions::iterator const& collision, aod::FT0s const& /*ft0s*/, BCsRun3 const& bcs, aod::Zdcs const&, MyTracks const&)
  void process(MyCollisions::iterator const& collision, aod::FT0s const& /*ft0s*/, BCsRun3 const& bcs, aod::Zdcs const&)
  {

    auto centrality = collision.centFT0C();

    if (bcs.size() != 0) {
      gRandom->SetSeed(bcs.iteratorAt(0).globalBC());
    }

    auto bc = collision.foundBC_as<BCsRun3>();
    if (!bc.has_zdc()) {
      return;
    }

    currentRunNumber = collision.foundBC_as<BCsRun3>().runNumber();
    auto vz = collision.posZ();
    auto vx = collision.posX();
    auto vy = collision.posY();
    bool triggerevent = false;

    float psiZDCC = -99;
    float psiZDCA = -99;
    auto qxZDCA = 0.0;
    auto qxZDCC = 0.0;
    auto qyZDCA = 0.0;
    auto qyZDCC = 0.0;
    auto sumA = 0.0;
    auto sumC = 0.0;

    auto zdc = bc.zdc();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergy = zdc.energySectorZNA();

    if (znaEnergy[0] < 0.0 || znaEnergy[1] < 0.0 || znaEnergy[2] < 0.0 || znaEnergy[3] < 0.0)
      return;
    if (zncEnergy[0] < 0.0 || zncEnergy[1] < 0.0 || zncEnergy[2] < 0.0 || zncEnergy[3] < 0.0)
      return;

    if (collision.sel8() && centrality < cfgCutCentrality && TMath::Abs(vz) < cfgCutVertex && collision.has_foundFT0() && eventSelected(collision, centrality) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TH2D>(ConfGainPath.value, bc.timestamp());
      }

      histos.fill(HIST("hCentrality"), centrality);
      if (QA) {
        histos.fill(HIST("Vz"), vz);
      }

      initCCDB(bc);

      auto gainequal = 1.0;
      constexpr float x[4] = {-1.75, 1.75, -1.75, 1.75};
      constexpr float y[4] = {-1.75, -1.75, 1.75, 1.75};

      for (std::size_t iChA = 0; iChA < 8; iChA++) {
        auto chanelid = iChA;
        if (useGainCallib && gainprofile) {
          gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz, chanelid + 0.5));
        }

        if (iChA < 4) {

          if (znaEnergy[iChA] <= 0.0) {
            return;
          } else {
            float ampl = gainequal * znaEnergy[iChA];
            qxZDCA = qxZDCA + ampl * x[iChA];
            qyZDCA = qyZDCA + ampl * y[iChA];
            sumA = sumA + ampl;
            histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
            // histos.fill(HIST("hZDCAmp"), chanelid + 0.5, vz, ampl);
          }
        } else {

          if (zncEnergy[iChA - 4] <= 0.0) {
            return;
          } else {
            float ampl = gainequal * zncEnergy[iChA - 4];
            qxZDCC = qxZDCC + ampl * x[iChA - 4];
            qyZDCC = qyZDCC + ampl * y[iChA - 4];
            sumC = sumC + ampl;
            histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
            // histos.fill(HIST("hZDCAmp"), chanelid + 0.5, vz, ampl);
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
        return;
      }

      if (useRecentereSp && (currentRunNumber != lastRunNumber)) {
        hrecentereSp = ccdb->getForTimeStamp<THnF>(ConfRecentereSp.value, bc.timestamp());
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

        if (useRecenteresqSp && hrecenteresqSp) {

          binCoords[4] = channelAxis->FindBin(0.5); // Channel for meanyA
          int globalBinMeansqxA = hrecenteresqSp->GetBin(binCoords);
          float meansqxA = hrecenteresqSp->GetBinContent(globalBinMeansqxA);

          // Repeat for other channels (meanyA, meanxC, meanyC)
          binCoords[4] = channelAxis->FindBin(1.5); // Channel for meanyA
          int globalBinMeansqyA = hrecenteresqSp->GetBin(binCoords);
          float meansqyA = hrecenteresqSp->GetBinContent(globalBinMeansqyA);

          binCoords[4] = channelAxis->FindBin(2.5); // Channel for meanxC
          int globalBinMeansqxC = hrecenteresqSp->GetBin(binCoords);
          float meansqxC = hrecenteresqSp->GetBinContent(globalBinMeansqxC);

          binCoords[4] = channelAxis->FindBin(3.5); // Channel for meanyC
          int globalBinMeansqyC = hrecenteresqSp->GetBin(binCoords);
          float meansqyC = hrecenteresqSp->GetBinContent(globalBinMeansqyC);

          qxZDCA = qxZDCA / meansqxA;
          qyZDCA = qyZDCA / meansqyA;
          qxZDCC = qxZDCC / meansqxC;
          qyZDCC = qyZDCC / meansqyC;

        } else {
          qxZDCA = qxZDCA;
          qyZDCA = qyZDCA;
          qxZDCC = qxZDCC;
          qyZDCC = qyZDCC;
        }
      }

      /*
  if (useRecentere && (currentRunNumber != lastRunNumber)) {
        hrecentere = ccdb->getForTimeStamp<TH3D>(ConfRecentere.value, bc.timestamp());
  }

  if (useRecentere && hrecentere) {
        if (recwitherror) {
  qxZDCA = (qxZDCA - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 0.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, vz, 0.5));
  qyZDCA = (qyZDCA - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 1.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, vz, 1.5));
  qxZDCC = (qxZDCC - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 2.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, vz, 2.5));
  qyZDCC = (qyZDCC - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 3.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, vz, 3.5));
        } else {

  qxZDCA = (qxZDCA - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 0.5)));
  qyZDCA = (qyZDCA - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 1.5)));
  qxZDCC = (qxZDCC - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 2.5)));
  qyZDCC = (qyZDCC - hrecentere->GetBinContent(hrecentere->FindBin(centrality, vz, 3.5)));
        }
  }


  if (useRecentereVxy && (currentRunNumber != lastRunNumber)) {
  hrecentereVxyQxA = ccdb->getForTimeStamp<TH3D>(ConfRecentereVxyQxA.value, bc.timestamp());
  hrecentereVxyQyA = ccdb->getForTimeStamp<TH3D>(ConfRecentereVxyQyA.value, bc.timestamp());
  hrecentereVxyQxC = ccdb->getForTimeStamp<TH3D>(ConfRecentereVxyQxC.value, bc.timestamp());
  hrecentereVxyQyC = ccdb->getForTimeStamp<TH3D>(ConfRecentereVxyQyC.value, bc.timestamp());
  }

  if (useRecentereVxy && hrecentereVxyQxA && hrecentereVxyQyA && hrecentereVxyQxC && hrecentereVxyQyC) {
  if (recwitherror) {
  qxZDCA = (qxZDCA - hrecentereVxyQxA->GetBinContent(hrecentereVxyQxA->FindBin(centrality, vx, vy))) / hrecentereVxyQxA->GetBinError(hrecentereVxyQxA->FindBin(centrality, vx, vy));
  qyZDCA = (qyZDCA - hrecentereVxyQyA->GetBinContent(hrecentereVxyQyA->FindBin(centrality, vx, vy))) / hrecentereVxyQyA->GetBinError(hrecentereVxyQyA->FindBin(centrality, vx, vy));
  qxZDCC = (qxZDCC - hrecentereVxyQxC->GetBinContent(hrecentereVxyQxC->FindBin(centrality, vx, vy))) / hrecentereVxyQxC->GetBinError(hrecentereVxyQxC->FindBin(centrality, vx, vy));
  qyZDCC = (qyZDCC - hrecentereVxyQyC->GetBinContent(hrecentereVxyQyC->FindBin(centrality, vx, vy))) / hrecentereVxyQyC->GetBinError(hrecentereVxyQyC->FindBin(centrality, vx, vy));
  } else {

  qxZDCA = (qxZDCA - hrecentereVxyQxA->GetBinContent(hrecentereVxyQxA->FindBin(centrality, vx, vy)));
  qyZDCA = (qyZDCA - hrecentereVxyQyA->GetBinContent(hrecentereVxyQyA->FindBin(centrality, vx, vy)));
  qxZDCC = (qxZDCC - hrecentereVxyQxC->GetBinContent(hrecentereVxyQxC->FindBin(centrality, vx, vy)));
  qyZDCC = (qyZDCC - hrecentereVxyQyC->GetBinContent(hrecentereVxyQyC->FindBin(centrality, vx, vy)));
  }
  }

      */

      psiZDCC = 1.0 * TMath::ATan2(qyZDCC, qxZDCC);
      psiZDCA = 1.0 * TMath::ATan2(qyZDCA, qxZDCA);

      histos.fill(HIST("hpQxZDCAC"), centrality, (qxZDCA * qxZDCC));
      histos.fill(HIST("hpQyZDCAC"), centrality, (qyZDCA * qyZDCC));
      histos.fill(HIST("hpQxZDCAQyZDCC"), centrality, (qxZDCA * qyZDCC));
      histos.fill(HIST("hpQxZDCCQyZDCA"), centrality, (qxZDCC * qyZDCA));

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
      }
      /*
      histos.fill(HIST("hpQxVZDCC"), centrality, vx, vy, qxZDCC);
      histos.fill(HIST("hpQyVZDCC"), centrality, vx, vy, qyZDCC);
      histos.fill(HIST("hpQxVZDCA"), centrality, vx, vy, qxZDCA);
      histos.fill(HIST("hpQyVZDCA"), centrality, vx, vy, qyZDCA);*/
      if (QA) {
        histos.fill(HIST("hpCosPsiAPsiC"), centrality, (TMath::Cos(psiZDCA - psiZDCC)));
        histos.fill(HIST("hpSinPsiAPsiC"), centrality, (TMath::Sin(psiZDCA - psiZDCC)));
      }
      histos.fill(HIST("PsiZDCA"), centrality, psiZDCA);
      histos.fill(HIST("PsiZDCC"), centrality, psiZDCC);

      lastRunNumber = currentRunNumber;
    }

    if (tablewrite) {
      spcalibrationtable(triggerevent, lastRunNumber, centrality, vz, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3], qxZDCA, qyZDCA, qxZDCC, qyZDCC, psiZDCC, psiZDCA);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<spvector>(cfgc)};
}
