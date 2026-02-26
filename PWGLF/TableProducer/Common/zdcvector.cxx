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
#include "Common/DataModel/PIDResponseTPC.h"
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

using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

struct zdcvector {

  Produces<aod::ZDCCalTables> zdccaltable;

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
  Configurable<bool> usemem{"usemem", true, "usemem"};

  struct : ConfigurableGroup {
    Configurable<int> vzFineNbins{"vzFineNbins", 20, "Number of bins in Vz fine histograms"};
    Configurable<float> lfinebinVz{"lfinebinVz", -10.0, "lower bin value in Vz fine histograms"};
    Configurable<float> hfinebinVz{"hfinebinVz", 10.0, "higher bin value in Vz fine histograms"};
    Configurable<int> centFineNbins{"centFineNbins", 16, "Number of bins in cent fine histograms"};
    Configurable<float> lfinebinCent{"lfinebinCent", 0.0, "lower bin value in cent fine histograms"};
    Configurable<float> hfinebinCent{"hfinebinCent", 80.0, "higher bin value in cent fine histograms"};
  } configbins;

  Configurable<bool> followpub{"followpub", true, "flag to use alphaZDC"};
  Configurable<bool> useGainCallib{"useGainCallib", false, "use gain calibration"};
  Configurable<bool> useCallibvertex{"useCallibvertex", false, "use calibration for vxy"};
  Configurable<std::string> confGainPath{"confGainPath", "Users/p/prottay/My/Object/NewPbPbpass4_10092024/gaincallib", "Path to gain calibration"};
  Configurable<std::string> confGainPathVxy{"confGainPathVxy", "Users/p/prottay/My/Object/swapcoords/PbPbpass4_20112024/recentervert", "Path to gain calibration for vxy"};

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
    AxisSpec vzfineAxis = {configbins.vzFineNbins, configbins.lfinebinVz, configbins.hfinebinVz, "vzfine"};
    AxisSpec centfineAxis = {configbins.centFineNbins, configbins.lfinebinCent, configbins.hfinebinCent, "V0M (%) fine"};
    AxisSpec VxyAxis = {2, 0, 2, "Vxy"};

    histos.add("htpcnsigmapi", "htpcnsigmapi", kTH1F, {{50, -10, 10.0}});
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0, 10.0}});
    histos.add("hCentrality", "hCentrality", kTH1F, {{centfineAxis}});
    histos.add("Vz", "Vz", kTH1F, {vzfineAxis});

    histos.add("ZDCAmp", "ZDCAmp", kTProfile2D, {channelZDCAxis, vzfineAxis});
    histos.add("ZDCAmpCommon", "ZDCAmpCommon", kTProfile2D, {{2, 0.0, 2.0}, vzfineAxis});
    histos.add("AvgVxy", "AvgVxy", kTProfile, {VxyAxis});

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

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>;

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
    auto qxZDCA = 0.0;
    auto qxZDCC = 0.0;
    auto qyZDCA = 0.0;
    auto qyZDCC = 0.0;
    auto sumA = 0.0;
    auto sumC = 0.0;

    auto bc = collision.foundBC_as<BCsRun3>();

    if (!bc.has_zdc()) {
      triggerevent = false;
      zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
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
      zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
      // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 2.5);

    if (znaEnergy[0] <= 0.0 || znaEnergy[1] <= 0.0 || znaEnergy[2] <= 0.0 || znaEnergy[3] <= 0.0) {
      triggerevent = false;
      zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
      // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);

    if (zncEnergy[0] <= 0.0 || zncEnergy[1] <= 0.0 || zncEnergy[2] <= 0.0 || zncEnergy[3] <= 0.0) {
      triggerevent = false;
      zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
      // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 4.5);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      triggerevent = false;
      zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
      // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 5.5);

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      triggerevent = false;
      zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
      // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 6.5);

    if (collision.sel8() && centrality > cfgCutCentralityMin && centrality < cfgCutCentralityMax && std::abs(vz) < cfgCutVertex && collision.has_foundFT0() && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TH2D>(confGainPath.value, bc.timestamp());
      }

      histos.fill(HIST("hEvtSelInfo"), 7.5);

      auto gainequal = 1.0;
      auto alphaZDC = 0.395;
      constexpr double x[4] = {-1.75, 1.75, -1.75, 1.75};
      constexpr double y[4] = {-1.75, -1.75, 1.75, 1.75};

      histos.fill(HIST("ZDCAmpCommon"), 0.5, vz, znaEnergycommon);
      histos.fill(HIST("ZDCAmpCommon"), 1.5, vz, zncEnergycommon);

      // int ntow = 8;
      constexpr std::size_t ntow = 8;
      for (std::size_t iChA = 0; iChA < ntow; iChA++) {
        auto chanelid = iChA;
        if (useGainCallib && gainprofile) {
          gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz + 0.00001, chanelid + 0.5));
        }

        if (iChA < 4) {

          if (znaEnergy[iChA] <= 0.0) {
            triggerevent = false;
            zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
            return;
          } else {
            double ampl = gainequal * znaEnergy[iChA];
            if (followpub) {
              ampl = std::pow(ampl, alphaZDC);
            }
            qxZDCA = qxZDCA - ampl * x[iChA];
            qyZDCA = qyZDCA + ampl * y[iChA];
            sumA = sumA + ampl;
            histos.fill(HIST("ZDCAmp"), chanelid + 0.5, vz, ampl);
          }
        } else {
          if (zncEnergy[iChA - 4] <= 0.0) {
            triggerevent = false;
            zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
            // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
            return;
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
        zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, 0.0, 0.0, 0.0, 0.0);
        // zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, znaEnergycommon, zncEnergycommon, znaEnergy[0], znaEnergy[1], znaEnergy[2], znaEnergy[3], zncEnergy[0], zncEnergy[1], zncEnergy[2], zncEnergy[3]);
        return;
      }

      histos.fill(HIST("hEvtSelInfo"), 8.5);
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);

      histos.fill(HIST("AvgVxy"), 0.5, vx);
      histos.fill(HIST("AvgVxy"), 1.5, vy);

      if (useCallibvertex && (currentRunNumber != lastRunNumber)) {
        gainprofilevxy = ccdb->getForTimeStamp<TProfile>(confGainPathVxy.value, bc.timestamp());
      }

      if (useCallibvertex) {
        vx = vx - gainprofilevxy->GetBinContent(1);
        vy = vy - gainprofilevxy->GetBinContent(2);
      }

      lastRunNumber = currentRunNumber;
    }
    zdccaltable(triggerevent, currentRunNumber, centrality, vx, vy, vz, qxZDCA, qxZDCC, qyZDCA, qyZDCC);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdcvector>(cfgc)};
}
