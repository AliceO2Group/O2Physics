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
/// \file   epvector.cxx
/// \author Sourav Kundu <sourav.kundu@cern.ch>
///
///
/// \brief  Task calculating the Q-vectors for each collision in a bunch crossing
///         (with or without corrections) and save the results in a dedicated table.
///

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
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

#include <TComplex.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct epvector {
  Produces<aod::EPCalibrationTables> epcalibrationtable;
  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  std::vector<o2::detectors::AlignParam>* offsetFV0;
  std::vector<int> TrkBPosLabel;
  std::vector<int> TrkBNegLabel;
  std::vector<float> qvecRe;
  std::vector<float> qvecIm;
  std::vector<float> qvecAmp;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  // Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Centrality cut"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutPTMax{"cfgCutPTMax", 3.0, "Max PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgMinEta{"cfgMinEta", 0.1, "Min Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 4, "Number of ITS cluster"};
  // Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgHarmonic{"cfgHarmonic", 2, "Harmonic for event plane calculation"};
  Configurable<bool> useGainCallib{"useGainCallib", true, "use gain calibration"};
  Configurable<bool> useRecentere{"useRecentere", true, "use Recentering"};
  Configurable<bool> useShift{"useShift", false, "use Shift"};
  Configurable<bool> useShift2{"useShift2", false, "use Shift for others"};
  Configurable<bool> useEventSelection{"useEventSelection", true, "Apply event selection centrality wise"};
  Configurable<bool> useTimeFrameCut{"useTimeFrameCut", true, "Reject Time Frame border events"};
  Configurable<bool> useITSFrameCut{"useITSFrameCut", true, "Reject ITS RO Frame border events"};
  Configurable<bool> usePileupCut{"usePileupCut", false, "Reject same bunch pileup"};
  Configurable<bool> useITSLayerCut{"useITSLayerCut", false, "Require good ITS layers"};
  Configurable<std::string> ConfGainPath{"ConfGainPath", "Users/s/skundu/My/Object/test100", "Path to gain calibration"};
  Configurable<std::string> ConfRecentere{"ConfRecentere", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for recentere"};
  Configurable<std::string> ConfShift{"ConfShift", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for Shift"};
  Configurable<std::string> ConfShiftFT0A{"ConfShiftFT0A", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for Shift FT0A"};
  Configurable<std::string> ConfShiftTPC{"ConfShiftTPC", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for Shift TPC"};
  Configurable<std::string> ConfShiftTPCL{"ConfShiftTPCL", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for Shift TPCL"};
  Configurable<std::string> ConfShiftTPCR{"ConfShiftTPCR", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for Shift TPCR"};
  ConfigurableAxis configAxisCentrality{"configAxisCentrality", {80, 0.0, 80}, "centrality bining"};
  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(o2::framework::InitContext&)
  {
    std::vector<double> occupancyBinning = {-0.5, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};

    const AxisSpec centAxis{configAxisCentrality, "V0M (%)"};
    // AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec multiplicity = {5000, -500, 500, "TPC Multiplicity"};
    AxisSpec amplitudeFT0 = {5000, 0, 10000, "FT0 amplitude"};
    AxisSpec channelFT0Axis = {220, 0.0, 220.0, "FT0 channel"};
    AxisSpec qxFT0Axis = {80000, -10000.0, 10000.0, "Qx"};
    AxisSpec qyFT0Axis = {80000, -10000.0, 10000.0, "Qy"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec vzAxis = {400, -20, 20, "vz"};
    AxisSpec resAxis = {200, -2, 2, "Resv2"};
    AxisSpec resAxisSP = {800, -80, 80, "ResSP"};
    AxisSpec qAxis = {100, 0, 10, "Q axis"};
    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {2, 0, 2, "basis"};
    AxisSpec occupancyAxis = {occupancyBinning, "occupancy"};

    histos.add("hCentrality", "hCentrality", kTH1F, {{8, 0, 80.0}});
    histos.add("Vz", "Vz", kTH1F, {{400, -20.0, 20.0}});
    histos.add("QxFT0C", "QxFT0C", kTH2F, {centAxis, qxFT0Axis});
    histos.add("QyFT0C", "QyFT0C", kTH2F, {centAxis, qyFT0Axis});
    histos.add("QxFT0A", "QxFT0A", kTH2F, {centAxis, qxFT0Axis});
    histos.add("QyFT0A", "QyFT0A", kTH2F, {centAxis, qyFT0Axis});
    histos.add("PsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
    histos.add("PsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
    histos.add("FT0Amp", "FT0Amp", kTH2F, {channelFT0Axis, amplitudeFT0});
    histos.add("QxTPC", "QxTPC", kTH2F, {centAxis, multiplicity});
    histos.add("QyTPC", "QyTPC", kTH2F, {centAxis, multiplicity});
    histos.add("PsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
    histos.add("QxTPCL", "QxTPCL", kTH2F, {centAxis, multiplicity});
    histos.add("QyTPCL", "QyTPCL", kTH2F, {centAxis, multiplicity});
    histos.add("PsiTPCL", "PsiTPCL", kTH2F, {centAxis, phiAxis});
    histos.add("QxTPCR", "QxTPCR", kTH2F, {centAxis, multiplicity});
    histos.add("QyTPCR", "QyTPCR", kTH2F, {centAxis, multiplicity});
    histos.add("PsiTPCR", "PsiTPCR", kTH2F, {centAxis, phiAxis});

    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH3F, {centAxis, resAxis, occupancyAxis});
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH3F, {centAxis, resAxis, occupancyAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH3F, {centAxis, resAxis, occupancyAxis});
    histos.add("ResFT0CTPCL", "ResFT0CTPCL", kTH3F, {centAxis, resAxis, occupancyAxis});
    histos.add("ResFT0CTPCR", "ResFT0CTPCR", kTH3F, {centAxis, resAxis, occupancyAxis});
    histos.add("ResTPCRTPCL", "ResTPCRTPCL", kTH3F, {centAxis, resAxis, occupancyAxis});

    histos.add("ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH3F, {centAxis, resAxisSP, occupancyAxis});
    histos.add("ResFT0CTPCSP", "ResFT0CTPCSP", kTH3F, {centAxis, resAxisSP, occupancyAxis});
    histos.add("ResFT0ATPCSP", "ResFT0ATPCSP", kTH3F, {centAxis, resAxisSP, occupancyAxis});
    histos.add("ResFT0CTPCLSP", "ResFT0CTPCLSP", kTH3F, {centAxis, resAxisSP, occupancyAxis});
    histos.add("ResFT0CTPCRSP", "ResFT0CTPCRSP", kTH3F, {centAxis, resAxisSP, occupancyAxis});
    histos.add("ResTPCRTPCLSP", "ResTPCRTPCLSP", kTH3F, {centAxis, resAxisSP, occupancyAxis});

    histos.add("QFT0C", "QFT0C", kTH3F, {centAxis, qAxis, occupancyAxis});
    histos.add("QFT0A", "QFT0A", kTH3F, {centAxis, qAxis, occupancyAxis});
    histos.add("QTPCL", "QTPCL", kTH3F, {centAxis, qAxis, occupancyAxis});
    histos.add("QTPCR", "QTPCR", kTH3F, {centAxis, qAxis, occupancyAxis});
    histos.add("QTPC", "QTPC", kTH3F, {centAxis, qAxis, occupancyAxis});

    histos.add("ShiftFT0C", "ShiftFT0C", kTProfile3D, {centAxis, basisAxis, shiftAxis});
    histos.add("ShiftFT0A", "ShiftFT0A", kTProfile3D, {centAxis, basisAxis, shiftAxis});
    histos.add("ShiftTPC", "ShiftTPC", kTProfile3D, {centAxis, basisAxis, shiftAxis});
    histos.add("ShiftTPCL", "ShiftTPCL", kTProfile3D, {centAxis, basisAxis, shiftAxis});
    histos.add("ShiftTPCR", "ShiftTPCR", kTProfile3D, {centAxis, basisAxis, shiftAxis});

    // Event selection cut additional - Alex
    // if (additionalEvsel) {
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
    // }

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

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered                                                                                                                                                                                                                                                                                                               // return 0;
    }
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;

    return 1;
  }

  double GetPhiFT0(int chno, double offsetX, double offsetY)
  {
    o2::ft0::Geometry ft0Det;
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return TMath::ATan2(chPos.Y() + offsetY, chPos.X() + offsetX);
  }

  double GetPhiInRange(double phi, double harmonic = 2)
  {
    double result = phi;
    double period = 2. * TMath::Pi() / harmonic;
    while (result < 0) {
      result = result + period;
    }
    while (result > period) {
      result = result - period;
    }
    return result;
  }

  double GetDeltaPsiSubInRange(double psi1, double psi2, double harmonic = 2)
  {
    double delta = psi1 - psi2;
    double period = TMath::Pi() / harmonic;
    if (TMath::Abs(delta) > period) {
      if (delta > 0.)
        delta -= 2. * period;
      else
        delta += 2. * period;
    }
    return delta;
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }
  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TProfile* gainprofile;
  TH2D* hrecentere;
  TProfile3D* shiftprofile;
  TProfile3D* shiftprofile2;
  TProfile3D* shiftprofile3;
  TProfile3D* shiftprofile4;
  TProfile3D* shiftprofile5;

  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs>;
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>>;

  void process(MyCollisions::iterator const& coll, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0s*/, aod::BCsWithTimestamps const&, MyTracks const& tracks)
  {
    auto centrality = coll.centFT0C();
    auto bc = coll.bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = coll.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto vz = coll.posZ();
    bool triggerevent = false;
    float psiFT0C = -99;
    float psiFT0A = -99;
    float psiTPC = -99;
    float psiTPCL = -99;
    float psiTPCR = -99;
    auto qxFT0A = 0.0;
    auto qxFT0C = 0.0;
    auto qyFT0A = 0.0;
    auto qyFT0C = 0.0;
    auto qxTPC = 0.0;
    auto qyTPC = 0.0;
    auto qxTPCL = 0.0;
    auto qyTPCL = 0.0;
    auto qxTPCR = 0.0;
    auto qyTPCR = 0.0;
    if (coll.sel8() && centrality < cfgCutCentrality && TMath::Abs(vz) < cfgCutVertex && coll.has_foundFT0() && (!useEventSelection || eventSelected(coll, centrality)) && (!useTimeFrameCut || coll.selection_bit(aod::evsel::kNoTimeFrameBorder)) && (!useITSFrameCut || coll.selection_bit(aod::evsel::kNoITSROFrameBorder)) && (!usePileupCut || coll.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!useITSLayerCut || coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TProfile>(ConfGainPath.value, bc.timestamp());
      }
      int occupancy = coll.trackOccupancyInTimeRange();
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);

      auto ft0 = coll.foundFT0();
      auto offsetFT0Ax = (*offsetFT0)[0].getX();
      auto offsetFT0Ay = (*offsetFT0)[0].getY();
      auto offsetFT0Cx = (*offsetFT0)[1].getX();
      auto offsetFT0Cy = (*offsetFT0)[1].getY();
      for (std::size_t iChA = 0; iChA < ft0.channelA().size(); iChA++) {
        auto chanelid = ft0.channelA()[iChA];
        auto gainequal = 1.0;
        if (useGainCallib) {
          gainequal = 1 / gainprofile->GetBinContent(gainprofile->FindBin(chanelid));
        }
        float ampl = gainequal * ft0.amplitudeA()[iChA];
        histos.fill(HIST("FT0Amp"), chanelid, ampl);
        auto phiA = GetPhiFT0(chanelid, offsetFT0Ax, offsetFT0Ay);
        qxFT0A = qxFT0A + ampl * TMath::Cos(cfgHarmonic.value * phiA);
        qyFT0A = qyFT0A + ampl * TMath::Sin(cfgHarmonic.value * phiA);
      }
      for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
        auto chanelid = ft0.channelC()[iChC] + 96;
        auto gainequal = 1.0;
        if (useGainCallib) {
          gainequal = 1 / gainprofile->GetBinContent(gainprofile->FindBin(chanelid));
        }
        float ampl = gainequal * ft0.amplitudeC()[iChC];
        histos.fill(HIST("FT0Amp"), chanelid, ampl);
        auto phiC = GetPhiFT0(chanelid, offsetFT0Cx, offsetFT0Cy);
        qxFT0C = qxFT0C + ampl * TMath::Cos(cfgHarmonic.value * phiC);
        qyFT0C = qyFT0C + ampl * TMath::Sin(cfgHarmonic.value * phiC);
      }

      for (auto& trk : tracks) {
        if (!selectionTrack(trk) || TMath::Abs(trk.eta()) > 0.8 || trk.pt() > cfgCutPTMax || TMath::Abs(trk.eta()) < cfgMinEta) {
          continue;
        }
        qxTPC = qxTPC + trk.pt() * TMath::Cos(cfgHarmonic.value * trk.phi());
        qyTPC = qyTPC + trk.pt() * TMath::Sin(cfgHarmonic.value * trk.phi());
        if (trk.eta() < 0.0) {
          qxTPCL = qxTPCL + trk.pt() * TMath::Cos(cfgHarmonic.value * trk.phi());
          qyTPCL = qyTPCL + trk.pt() * TMath::Sin(cfgHarmonic.value * trk.phi());
        }
        if (trk.eta() > 0.0) {
          qxTPCR = qxTPCR + trk.pt() * TMath::Cos(cfgHarmonic.value * trk.phi());
          qyTPCR = qyTPCR + trk.pt() * TMath::Sin(cfgHarmonic.value * trk.phi());
        }
      }
      if (useRecentere && (currentRunNumber != lastRunNumber)) {
        hrecentere = ccdb->getForTimeStamp<TH2D>(ConfRecentere.value, bc.timestamp());
      }
      if (useRecentere) {
        qxFT0A = (qxFT0A - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 0.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 0.5));
        qyFT0A = (qyFT0A - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 1.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 1.5));
        qxFT0C = (qxFT0C - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 2.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 2.5));
        qyFT0C = (qyFT0C - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 3.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 3.5));
        qxTPC = (qxTPC - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 4.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 4.5));
        qyTPC = (qyTPC - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 5.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 5.5));
        qxTPCL = (qxTPCL - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 6.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 6.5));
        qyTPCL = (qyTPCL - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 7.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 7.5));
        qxTPCR = (qxTPCR - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 8.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 8.5));
        qyTPCR = (qyTPCR - hrecentere->GetBinContent(hrecentere->FindBin(centrality, 9.5))) / hrecentere->GetBinError(hrecentere->FindBin(centrality, 9.5));
      }
      psiFT0C = (1.0 / cfgHarmonic.value) * TMath::ATan2(qyFT0C, qxFT0C);
      psiFT0A = (1.0 / cfgHarmonic.value) * TMath::ATan2(qyFT0A, qxFT0A);
      psiTPC = (1.0 / cfgHarmonic.value) * TMath::ATan2(qyTPC, qxTPC);
      psiTPCL = (1.0 / cfgHarmonic.value) * TMath::ATan2(qyTPCL, qxTPCL);
      psiTPCR = (1.0 / cfgHarmonic.value) * TMath::ATan2(qyTPCR, qxTPCR);

      if (useShift && (currentRunNumber != lastRunNumber)) {
        shiftprofile = ccdb->getForTimeStamp<TProfile3D>(ConfShift.value, bc.timestamp());
        if (useShift2) {
          shiftprofile2 = ccdb->getForTimeStamp<TProfile3D>(ConfShiftFT0A.value, bc.timestamp());
          shiftprofile3 = ccdb->getForTimeStamp<TProfile3D>(ConfShiftTPC.value, bc.timestamp());
          shiftprofile4 = ccdb->getForTimeStamp<TProfile3D>(ConfShiftTPCL.value, bc.timestamp());
          shiftprofile5 = ccdb->getForTimeStamp<TProfile3D>(ConfShiftTPCR.value, bc.timestamp());
        }
      }
      if (useShift) {
        auto deltapsiFT0C = 0.0;
        auto deltapsiFT0A = 0.0;
        auto deltapsiTPC = 0.0;
        auto deltapsiTPCL = 0.0;
        auto deltapsiTPCR = 0.0;
        for (int ishift = 1; ishift <= 10; ishift++) {
          auto coeffshiftxFT0C = shiftprofile->GetBinContent(shiftprofile->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyFT0C = shiftprofile->GetBinContent(shiftprofile->FindBin(centrality, 1.5, ishift - 0.5));
          deltapsiFT0C = deltapsiFT0C + ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * cfgHarmonic.value * psiFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * cfgHarmonic.value * psiFT0C)));
          if (useShift2) {
            auto coeffshiftxFT0A = shiftprofile2->GetBinContent(shiftprofile2->FindBin(centrality, 0.5, ishift - 0.5));
            auto coeffshiftyFT0A = shiftprofile2->GetBinContent(shiftprofile2->FindBin(centrality, 1.5, ishift - 0.5));

            auto coeffshiftxTPC = shiftprofile3->GetBinContent(shiftprofile3->FindBin(centrality, 0.5, ishift - 0.5));
            auto coeffshiftyTPC = shiftprofile3->GetBinContent(shiftprofile3->FindBin(centrality, 1.5, ishift - 0.5));

            auto coeffshiftxTPCL = shiftprofile4->GetBinContent(shiftprofile4->FindBin(centrality, 0.5, ishift - 0.5));
            auto coeffshiftyTPCL = shiftprofile4->GetBinContent(shiftprofile4->FindBin(centrality, 1.5, ishift - 0.5));

            auto coeffshiftxTPCR = shiftprofile5->GetBinContent(shiftprofile5->FindBin(centrality, 0.5, ishift - 0.5));
            auto coeffshiftyTPCR = shiftprofile5->GetBinContent(shiftprofile5->FindBin(centrality, 1.5, ishift - 0.5));
            deltapsiFT0A = deltapsiFT0A + ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * cfgHarmonic.value * psiFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * cfgHarmonic.value * psiFT0A)));
            deltapsiTPC = deltapsiTPC + ((1 / (1.0 * ishift)) * (-coeffshiftxTPC * TMath::Cos(ishift * cfgHarmonic.value * psiTPC) + coeffshiftyTPC * TMath::Sin(ishift * cfgHarmonic.value * psiTPC)));
            deltapsiTPCL = deltapsiTPCL + ((1 / (1.0 * ishift)) * (-coeffshiftxTPCL * TMath::Cos(ishift * cfgHarmonic.value * psiTPCL) + coeffshiftyTPCL * TMath::Sin(ishift * cfgHarmonic.value * psiTPCL)));
            deltapsiTPCR = deltapsiTPCR + ((1 / (1.0 * ishift)) * (-coeffshiftxTPCR * TMath::Cos(ishift * cfgHarmonic.value * psiTPCR) + coeffshiftyTPCR * TMath::Sin(ishift * cfgHarmonic.value * psiTPCR)));
          }
        }
        psiFT0C = psiFT0C + deltapsiFT0C;
        psiFT0A = psiFT0A + deltapsiFT0A;
        psiTPC = psiTPC + deltapsiTPC;
        psiTPCL = psiTPCL + deltapsiTPCL;
        psiTPCR = psiTPCR + deltapsiTPCR;
      }
      histos.fill(HIST("QxFT0C"), centrality, qxFT0C);
      histos.fill(HIST("QyFT0C"), centrality, qyFT0C);
      histos.fill(HIST("QxFT0A"), centrality, qxFT0A);
      histos.fill(HIST("QyFT0A"), centrality, qyFT0A);
      histos.fill(HIST("PsiFT0A"), centrality, psiFT0A);
      histos.fill(HIST("PsiFT0C"), centrality, psiFT0C);
      histos.fill(HIST("QxTPC"), centrality, qxTPC);
      histos.fill(HIST("QyTPC"), centrality, qyTPC);
      histos.fill(HIST("PsiTPC"), centrality, psiTPC);
      histos.fill(HIST("QxTPCL"), centrality, qxTPCL);
      histos.fill(HIST("QyTPCL"), centrality, qyTPCL);
      histos.fill(HIST("PsiTPCL"), centrality, psiTPCL);
      histos.fill(HIST("QxTPCR"), centrality, qxTPCR);
      histos.fill(HIST("QyTPCR"), centrality, qyTPCR);
      histos.fill(HIST("PsiTPCR"), centrality, psiTPCR);

      histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(cfgHarmonic.value * (psiFT0C - psiFT0A)), occupancy);
      histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(cfgHarmonic.value * (psiFT0C - psiTPC)), occupancy);
      histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(cfgHarmonic.value * (psiFT0A - psiTPC)), occupancy);
      histos.fill(HIST("ResFT0CTPCL"), centrality, TMath::Cos(cfgHarmonic.value * (psiFT0C - psiTPCL)), occupancy);
      histos.fill(HIST("ResFT0CTPCR"), centrality, TMath::Cos(cfgHarmonic.value * (psiFT0C - psiTPCR)), occupancy);
      histos.fill(HIST("ResTPCRTPCL"), centrality, TMath::Cos(cfgHarmonic.value * (psiTPCR - psiTPCL)), occupancy);

      double qFT0Cmag = TMath::Sqrt(qxFT0C * qxFT0C + qyFT0C * qyFT0C);
      double qFT0Amag = TMath::Sqrt(qxFT0A * qxFT0A + qyFT0A * qyFT0A);
      double qTPCmag = TMath::Sqrt(qxTPC * qxTPC + qyTPC * qyTPC);
      double qTPCLmag = TMath::Sqrt(qxTPCL * qxTPCL + qyTPCL * qyTPCL);
      double qTPCRmag = TMath::Sqrt(qxTPCR * qxTPCR + qyTPCR * qyTPCR);

      histos.fill(HIST("QTPC"), centrality, qTPCmag, occupancy);
      histos.fill(HIST("QTPCL"), centrality, qTPCLmag, occupancy);
      histos.fill(HIST("QTPCR"), centrality, qTPCRmag, occupancy);
      histos.fill(HIST("QFT0C"), centrality, qFT0Cmag, occupancy);
      histos.fill(HIST("QFT0A"), centrality, qFT0Amag, occupancy);

      histos.fill(HIST("ResFT0CFT0ASP"), centrality, qFT0Cmag * qFT0Amag * TMath::Cos(cfgHarmonic.value * (psiFT0C - psiFT0A)), occupancy);
      histos.fill(HIST("ResFT0CTPCSP"), centrality, qFT0Cmag * qTPCmag * TMath::Cos(cfgHarmonic.value * (psiFT0C - psiTPC)), occupancy);
      histos.fill(HIST("ResFT0ATPCSP"), centrality, qFT0Amag * qTPCmag * TMath::Cos(cfgHarmonic.value * (psiFT0A - psiTPC)), occupancy);
      histos.fill(HIST("ResFT0CTPCLSP"), centrality, qFT0Cmag * qTPCLmag * TMath::Cos(cfgHarmonic.value * (psiFT0C - psiTPCL)), occupancy);
      histos.fill(HIST("ResFT0CTPCRSP"), centrality, qFT0Cmag * qTPCRmag * TMath::Cos(cfgHarmonic.value * (psiFT0C - psiTPCR)), occupancy);
      histos.fill(HIST("ResTPCRTPCLSP"), centrality, qTPCRmag * qTPCLmag * TMath::Cos(cfgHarmonic.value * (psiTPCR - psiTPCL)), occupancy);

      for (int ishift = 1; ishift <= 10; ishift++) {
        histos.fill(HIST("ShiftFT0C"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * cfgHarmonic.value * psiFT0C));
        histos.fill(HIST("ShiftFT0C"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * cfgHarmonic.value * psiFT0C));

        histos.fill(HIST("ShiftFT0A"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * cfgHarmonic.value * psiFT0A));
        histos.fill(HIST("ShiftFT0A"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * cfgHarmonic.value * psiFT0A));

        histos.fill(HIST("ShiftTPC"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * cfgHarmonic.value * psiTPC));
        histos.fill(HIST("ShiftTPC"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * cfgHarmonic.value * psiTPC));

        histos.fill(HIST("ShiftTPCL"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * cfgHarmonic.value * psiTPCL));
        histos.fill(HIST("ShiftTPCL"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * cfgHarmonic.value * psiTPCL));

        histos.fill(HIST("ShiftTPCR"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * cfgHarmonic.value * psiTPCR));
        histos.fill(HIST("ShiftTPCR"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * cfgHarmonic.value * psiTPCR));
      }
      lastRunNumber = currentRunNumber;
    }
    epcalibrationtable(triggerevent, centrality, psiFT0C, psiFT0A, psiTPC, psiTPCL, psiTPCR, TMath::Sqrt(qxFT0C * qxFT0C + qyFT0C * qyFT0C), TMath::Sqrt(qxFT0A * qxFT0A + qyFT0A * qyFT0A), TMath::Sqrt(qxTPC * qxTPC + qyTPC * qyTPC), TMath::Sqrt(qxTPCL * qxTPCL + qyTPCL * qyTPCL), TMath::Sqrt(qxTPCR * qxTPCR + qyTPCR * qyTPCR));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<epvector>(cfgc)};
}
