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
// Calibration task
// prottay.das@cern.ch

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/DataModel/FT0Corrected.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "PWGLF/DataModel/SPCalibrationTablesFT0C.h"
// #include "SPCalibrationTablesFT0C.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct spvectorft0c {

  Produces<aod::SPCalibrationTablesFT0C> spcalibrationtableFT0C;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  std::vector<o2::detectors::AlignParam>* offsetFV0;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};

  Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
  Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
  Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};
  Configurable<int> PhiNbins{"PhiNbins", 100, "Number of bins in phi histogram"};

  Configurable<bool> useGainCallib{"useGainCallib", false, "use gain calibration"};
  Configurable<bool> useRecentere{"useRecentere", false, "use Recentering"};
  Configurable<bool> useShift{"useShift", false, "shift histograms"};
  Configurable<std::string> ConfGainPath{"ConfGainPath", "Users/p/prottay/My/Object/NewPbPbpass4_10092024/gaincallib", "Path to gain calibration"};
  Configurable<std::string> ConfRecentere{"ConfRecentere", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for recentere"};
  Configurable<std::string> ConfShift{"ConfShift", "Users/p/prottay/My/Object/Testinglocaltree/shiftcallib2", "Path to shift"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configchannelFT0CAxis{"configchanelFT0CAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 90.0, 120, 150, 180, 200, 250}, "Channel Axis"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {

    AxisSpec vzfineAxis = {20, -10, 10, "vzfine"};
    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {2, 0, 2, "basis"};
    AxisSpec qxFT0CAxis = {QxyNbins, lbinQxy, hbinQxy, "Qx"};
    AxisSpec phiAxis = {PhiNbins, -6.28, 6.28, "phi"};

    histos.add("hCentrality", "hCentrality", kTH1F, {{configcentAxis}});
    histos.add("Vz", "Vz", kTH1F, {vzfineAxis});
    histos.add("FT0CAmp", "FT0CAmp", kTProfile2D, {configchannelFT0CAxis, vzfineAxis});
    histos.add("ShiftFT0C", "ShiftFT0C", kTProfile3D, {configcentAxis, basisAxis, shiftAxis});
    histos.add("hcentQxFT0C", "hcentQxFT0C", kTH2F, {{configcentAxis}, {qxFT0CAxis}});
    histos.add("hcentQyFT0C", "hcentQyFT0C", kTH2F, {{configcentAxis}, {qxFT0CAxis}});
    histos.add("PsiFT0C", "PsiFT0C", kTH2F, {configcentAxis, phiAxis});

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", cfgCcdbParam.nolaterthan.value);
    offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", cfgCcdbParam.nolaterthan.value);
    printf("Offset for FT0A: x = %.3f y = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY());
    printf("Offset for FT0C: x = %.3f y = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY());
    printf("Offset for FV0-left: x = %.3f y = %.3f\n", (*offsetFV0)[0].getX(), (*offsetFV0)[0].getY());
    printf("Offset for FV0-right: x = %.3f y = %.3f\n", (*offsetFV0)[1].getX(), (*offsetFV0)[1].getY());
  }

  double GetPhiFT0(int chno, double offsetX, double offsetY)
  {
    o2::ft0::Geometry ft0Det;
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return TMath::ATan2(chPos.Y() + offsetY, chPos.X() + offsetX);
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::FT0sCorrected, aod::Mults>>;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TH2D* gainprofile;
  TH2D* hrecentere;
  TProfile3D* shiftprofile;

  void process(EventCandidates::iterator const& collision, aod::FT0s const& /*ft0s*/, aod::BCsWithTimestamps const&)
  {

    auto centrality = collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto vz = collision.posZ();

    bool triggerevent = false;
    double psiFT0C = -99;
    auto qxFT0C = 0.0;
    auto qyFT0C = 0.0;
    auto sum = 0.0;

    if (collision.sel8() && centrality > cfgCutCentralityMin && centrality < cfgCutCentralityMax && TMath::Abs(vz) < cfgCutVertex && collision.has_foundFT0() && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TH2D>(ConfGainPath.value, bc.timestamp());
      }

      auto ft0 = collision.foundFT0();
      auto offsetFT0Cx = (*offsetFT0)[1].getX();
      auto offsetFT0Cy = (*offsetFT0)[1].getY();

      for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
        auto chanelid = ft0.channelC()[iChC] + 96;
        auto gainequal = 1.0;
        if (useGainCallib) {
          gainequal = gainprofile->GetBinContent(gainprofile->FindBin(vz, chanelid + 0.5));
        }
        float ampl = gainequal * ft0.amplitudeC()[iChC];
        auto phiC = GetPhiFT0(chanelid, offsetFT0Cx, offsetFT0Cy);
        qxFT0C = qxFT0C + ampl * TMath::Cos(1.0 * phiC);
        qyFT0C = qyFT0C + ampl * TMath::Sin(1.0 * phiC);
        sum = sum + ampl;
        histos.fill(HIST("FT0CAmp"), chanelid + 0.5, vz, ampl);
      }

      if (sum > 0) {
        qxFT0C = qxFT0C / sum;
        qyFT0C = qyFT0C / sum;
      }

      if (sum <= 1e-4) {
        qxFT0C = 0.0;
        qyFT0C = 0.0;
        triggerevent = false;
        spcalibrationtableFT0C(triggerevent, currentRunNumber, centrality, qxFT0C, qyFT0C, psiFT0C);
        return;
      }

      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);

      if (useRecentere && (currentRunNumber != lastRunNumber)) {
        hrecentere = ccdb->getForTimeStamp<TH2D>(ConfRecentere.value, bc.timestamp());
      }

      if (useRecentere) {
        qxFT0C = (qxFT0C - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 0.5)));
        qyFT0C = (qyFT0C - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 1.5)));
      }

      psiFT0C = 0.5 * TMath::ATan2(qyFT0C, qxFT0C);

      int nshift = 10; // no. of iterations

      if (useShift && (currentRunNumber != lastRunNumber)) {
        shiftprofile = ccdb->getForTimeStamp<TProfile3D>(ConfShift.value, bc.timestamp());
      }

      if (useShift) {
        auto deltapsiFT0C = 0.0;
        for (int ishift = 1; ishift <= nshift; ishift++) {
          auto coeffshiftxFT0C = shiftprofile->GetBinContent(shiftprofile->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyFT0C = shiftprofile->GetBinContent(shiftprofile->FindBin(centrality, 1.5, ishift - 0.5));
          deltapsiFT0C = deltapsiFT0C + ((2 / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * 1.0 * psiFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * 1.0 * psiFT0C)));
        }
        psiFT0C = psiFT0C + deltapsiFT0C;
      }

      for (int ishift = 1; ishift <= nshift; ishift++) {
        histos.fill(HIST("ShiftFT0C"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * 1.0 * psiFT0C));
        histos.fill(HIST("ShiftFT0C"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * 1.0 * psiFT0C));
      }

      histos.fill(HIST("hcentQxFT0C"), centrality, qxFT0C);
      histos.fill(HIST("hcentQyFT0C"), centrality, qyFT0C);
      histos.fill(HIST("PsiFT0C"), centrality, psiFT0C);

      lastRunNumber = currentRunNumber;
    }

    spcalibrationtableFT0C(triggerevent, currentRunNumber, centrality, qxFT0C, qyFT0C, psiFT0C);
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<spvectorft0c>(cfgc, TaskName{"spvectorft0c"})};
}
