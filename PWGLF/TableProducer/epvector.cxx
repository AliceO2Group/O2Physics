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

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TMath.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/FT0Corrected.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"

// #include "Common/Core/EventPlaneHelper.h"
// #include "Common/DataModel/Qvectors.h"

// o2 includes.
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct epvector {
  Produces<aod::EPCalibrationTables> epcalibrationtable;
  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
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
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Centrality cut"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutPTMax{"cfgCutPTMax", 3.0, "Max PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 4, "Number of ITS cluster"};
  Configurable<bool> useGainCallib{"useGainCallib", true, "use gain calibration"};
  Configurable<bool> useRecentere{"useRecentere", true, "use Recentering"};
  Configurable<std::string> ConfGainPath{"ConfGainPath", "Users/s/skundu/My/Object/test100", "Path to gain calibration"};
  Configurable<std::string> ConfRecentere{"ConfRecentere", "Users/s/skundu/My/Object/Finaltest2/recenereall", "Path for recentere"};
  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
    AxisSpec multiplicity = {5000, -500, 500, "TPC Multiplicity"};
    AxisSpec amplitudeFT0 = {5000, 0, 10000, "FT0 amplitude"};
    AxisSpec channelFT0Axis = {220, 0.0, 220.0, "FT0 channel"};
    AxisSpec qxFT0Axis = {2000, -100.0, 100.0, "Qx"};
    AxisSpec qyFT0Axis = {2000, -100.0, 100.0, "Qy"};
    AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    AxisSpec vzAxis = {400, -20, 20, "vz"};
    AxisSpec resAxis = {400, -2, 2, "vz"};

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
    histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});

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

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  double GetDeltaPsiSubInRange(double psi1, double psi2)
  {
    double delta = psi1 - psi2;
    if (TMath::Abs(delta) > TMath::Pi() / 2) {
      if (delta > 0.)
        delta -= 2. * TMath::Pi() / 2;
      else
        delta += 2. * TMath::Pi() / 2;
    }
    return delta;
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }
  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TProfile* gainprofile;
  TH2D* hrecentere;

  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0Cs>;
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>>;

  void process(MyCollisions::iterator const& coll, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::BCsWithTimestamps const&, MyTracks const& tracks)
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
    if (coll.sel8() && centrality < cfgCutCentrality && TMath::Abs(vz) < cfgCutVertex && coll.has_foundFT0()) {
      triggerevent = true;
      if (useGainCallib && (currentRunNumber != lastRunNumber)) {
        gainprofile = ccdb->getForTimeStamp<TProfile>(ConfGainPath.value, bc.timestamp());
      }

      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("Vz"), vz);
      auto qxFT0A = 0.0;
      auto qxFT0C = 0.0;
      auto qyFT0A = 0.0;
      auto qyFT0C = 0.0;

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
        qxFT0A = qxFT0A + ampl * TMath::Cos(2.0 * phiA);
        qyFT0A = qyFT0A + ampl * TMath::Sin(2.0 * phiA);
      }
      for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
        auto chanelid = ft0.channelC()[iChC] + 96;
        // printf("Offset for FT0A: x = %d y = %d\n", chanelid, chanelid-96);
        auto gainequal = 1.0;
        if (useGainCallib) {
          gainequal = 1 / gainprofile->GetBinContent(gainprofile->FindBin(chanelid));
        }
        float ampl = gainequal * ft0.amplitudeC()[iChC];
        histos.fill(HIST("FT0Amp"), chanelid, ampl);
        auto phiC = GetPhiFT0(chanelid, offsetFT0Cx, offsetFT0Cy);
        qxFT0C = qxFT0C + ampl * TMath::Cos(2.0 * phiC);
        qyFT0C = qyFT0C + ampl * TMath::Sin(2.0 * phiC);
      }

      auto qxTPC = 0.0;
      auto qyTPC = 0.0;
      auto qxTPCL = 0.0;
      auto qyTPCL = 0.0;
      auto qxTPCR = 0.0;
      auto qyTPCR = 0.0;

      for (auto& trk : tracks) {
        if (!selectionTrack(trk) || abs(trk.eta()) > 0.8 || trk.pt() > cfgCutPTMax) {
          continue;
        }
        qxTPC = qxTPC + trk.pt() * TMath::Cos(2.0 * trk.phi());
        qyTPC = qyTPC + trk.pt() * TMath::Sin(2.0 * trk.phi());
        if (trk.eta() < 0.0) {
          qxTPCL = qxTPCL + trk.pt() * TMath::Cos(2.0 * trk.phi());
          qyTPCL = qyTPCL + trk.pt() * TMath::Sin(2.0 * trk.phi());
        }
        if (trk.eta() > 0.0) {
          qxTPCR = qxTPCR + trk.pt() * TMath::Cos(2.0 * trk.phi());
          qyTPCR = qyTPCR + trk.pt() * TMath::Sin(2.0 * trk.phi());
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
      psiFT0C = 0.5 * TMath::ATan2(qyFT0C, qxFT0C);
      psiFT0A = 0.5 * TMath::ATan2(qyFT0A, qxFT0A);
      psiTPC = 0.5 * TMath::ATan2(qyTPC, qxTPC);
      psiTPCL = 0.5 * TMath::ATan2(qyTPCL, qxTPCL);
      psiTPCR = 0.5 * TMath::ATan2(qyTPCR, qxTPCR);
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

      histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));
      lastRunNumber = currentRunNumber;
    }
    epcalibrationtable(triggerevent, centrality, psiFT0C, psiFT0A, psiTPC, psiTPCL, psiTPCR);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<epvector>(cfgc)};
}
