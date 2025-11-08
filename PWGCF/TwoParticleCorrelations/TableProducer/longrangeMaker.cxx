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
/// \file longrangeMaker.cxx
///
/// \brief task derived table definition for long range correlation
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since October 28, 2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/LongRangeDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TComplex.h>
#include <TH1F.h>
#include <TMath.h>
#include <TPDGCode.h>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::aod::evsel;
using namespace o2::constants::math;

auto static constexpr KminFt0cCell = 96;
auto static constexpr PionTrackN = 1;
auto static constexpr KaonTrackN = 2;
auto static constexpr ProtonTrackN = 3;
AxisSpec axisEvent{15, 0.5, 15.5, "#Event", "EventAxis"};

enum KindOfParticles {
  PIONS,
  KAONS,
  PROTONS
};

enum KindOfV0 {
  kLambda = 0,
  kAntiLambda
};

struct LongrangeMaker {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  struct : ConfigurableGroup {
    Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", false, "Enable SameBunchPileup cut"};
    Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", false, "Enable GoodZvtxFT0vsPV cut"};
    Configurable<bool> isApplyGoodITSLayersAll{"isApplyGoodITSLayersAll", false, "Enable GoodITSLayersAll cut"};
    Configurable<bool> isApplyExtraCorrCut{"isApplyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
    Configurable<float> npvTracksCut{"npvTracksCut", 1.0f, "Apply extra NPVtracks cut"};
    Configurable<float> ft0cCut{"ft0cCut", 1.0f, "Apply extra FT0C cut"};
    Configurable<bool> isApplyNoCollInTimeRangeStandard{"isApplyNoCollInTimeRangeStandard", false, "Enable NoCollInTimeRangeStandard cut"};
    Configurable<bool> isApplyNoCollInRofStandard{"isApplyNoCollInRofStandard", false, "Enable NoCollInRofStandard cut"};
    Configurable<bool> isApplyNoHighMultCollInPrevRof{"isApplyNoHighMultCollInPrevRof", false, "Enable NoHighMultCollInPrevRof cut"};
    Configurable<bool> isApplyCentFT0C{"isApplyCentFT0C", false, "Centrality based on FT0C"};
    Configurable<bool> isApplyCentFV0A{"isApplyCentFV0A", false, "Centrality based on FV0A"};
    Configurable<bool> isApplyCentFT0M{"isApplyCentFT0M", false, "Centrality based on FT0A + FT0C"};
  } cfgevtsel;

  struct : ConfigurableGroup {
    Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "Eta range to consider"};
    Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.2f, "minimum accepted track pT"};
    Configurable<float> cfgPtCutMax{"cfgPtCutMax", 10.0f, "maximum accepted track pT"};
    Configurable<float> cfgPtCutMult{"cfgPtCutMult", 3.0f, "maximum track pT for multiplicity classification"};
    Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "cut on minimum number of TPC crossed rows"};
    Configurable<float> minTPCNClsFound{"minTPCNClsFound", 50.f, "cut on minimum value of TPC found clusters"};
    Configurable<float> maxDcaZ{"maxDcaZ", 2.f, "cut on maximum abs value of DCA z"};
    Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "cut on maximum value of TPC chi2 per cluster"};
  } cfgtrksel;

  struct : ConfigurableGroup {
    Configurable<float> cfigMftEtaMax{"cfigMftEtaMax", -2.5f, "Maximum MFT eta cut"};
    Configurable<float> cfigMftEtaMin{"cfigMftEtaMin", -3.6f, "Minimum MFT eta cut"};
    Configurable<float> cfigMftDcaxy{"cfigMftDcaxy", 2.0f, "cut on DCA xy for MFT tracks"};
    Configurable<int> cfigMftCluster{"cfigMftCluster", 5, "cut on MFT Cluster"};
    Configurable<bool> useMftPtCut{"useMftPtCut", true, "Choose to apply MFT track pT cut"};
    Configurable<float> cfgMftPtCutMin{"cfgMftPtCutMin", 0.f, "minimum accepted MFT track pT"};
    Configurable<float> cfgMftPtCutMax{"cfgMftPtCutMax", 10.f, "maximum accepted MFT track pT"};
  } cfgmfttrksel;

  struct : ConfigurableGroup {
    Configurable<float> minTPCcrossedrows{"minTPCcrossedrows", 70.f, "cut on minimum number of crossed rows in TPC"};
    Configurable<float> minTPCcrossedrowsoverfindcls{"minTPCcrossedrowsoverfindcls", 0.8f, "cut on minimum value of the ratio between crossed rows and findable clusters in the TPC"};
    Configurable<float> v0etaCut{"v0etaCut", 0.8f, "maximum v0 track pseudorapidity"};
    Configurable<float> v0ptCut{"v0ptCut", 0.15f, "minimum v0 track pT"};
    Configurable<float> minK0sMass{"minK0sMass", 0.4f, "minimum mass for K0s"};
    Configurable<float> maxK0sMass{"maxK0sMass", 0.6f, "maximum mass for K0s"};
    Configurable<float> maxK0sRadius{"maxK0sRadius", 30.0f, "Maximum decay radius (cm) for K0s"};
    Configurable<float> minK0sRadius{"minK0sRadius", 1.2f, "Minimum decay radius (cm) for K0s"};
    Configurable<float> minK0sCosPa{"minK0sCosPa", 0.993f, "Minimum cosine of pointing angle for K0s"};
    Configurable<float> maxDcaV0DauK0s{"maxDcaV0DauK0s", 0.8f, "Maximum DCA among the V0 daughters (cm) for K0s"};
    Configurable<float> minqtArmenterosForK0s{"minqtArmenterosForK0s", 0.2f, "Minimum Armenteros' qt for K0s"};
    Configurable<float> maxK0sLifeTime{"maxK0sLifeTime", 20.0f, "Maximum K0s lifetime (in cm)"};
    Configurable<float> daughPIDCuts{"daughPIDCuts", 4.0f, "PID nsigma for V0s"};
    Configurable<float> minV0DcaPiK0s{"minV0DcaPiK0s", 0.1f, "Min V0 pion DCA for K0s"};

    Configurable<float> minLambdaMass{"minLambdaMass", 1.07f, "minimum mass for Lambda"};
    Configurable<float> maxLambdaMass{"maxLambdaMass", 1.17f, "maximum mass for Lambda"};
    Configurable<float> maxLambdaRadius{"maxLambdaRadius", 30.0f, "Maximum decay radius (cm) for Lambda"};
    Configurable<float> minLambdaRadius{"minLambdaRadius", 1.2f, "Minimum decay radius (cm) for Lambda"};
    Configurable<float> minLambdaCosPa{"minLambdaCosPa", 0.993f, "Minimum cosine of pointing angle for Lambda"};
    Configurable<float> maxDcaV0DauLambda{"maxDcaV0DauLambda", 0.8f, "Maximum DCA among the V0 daughters (cm) for K0s"};
    Configurable<float> minV0DcaPiLambda{"minV0DcaPiLambda", 0.2f, "Min V0 pion DCA for Lambda"};
    Configurable<float> minV0DcaPr{"minV0DcaPr", 0.07f, "Min V0 proton DCA for Lambda"};
    Configurable<float> maxLambdaLifeTime{"maxLambdaLifeTime", 30.0f, "Maximum Lambda lifetime (in cm)"};
  } cfgv0trksel;

  Configurable<std::vector<double>> itsNsigmaPidCut{"itsNsigmaPidCut", std::vector<double>{3, 2.5, 2, -3, -2.5, -2}, "ITS n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<double>> tpcNsigmaPidCut{"tpcNsigmaPidCut", std::vector<double>{1.5, 1.5, 1.5, -1.5, -1.5, -1.5}, "TPC n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<double>> tofNsigmaPidCut{"tofNsigmaPidCut", std::vector<double>{1.5, 1.5, 1.5, -1.5, -1.5, -1.5}, "TOF n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<float> cfgTofPidPtCut{"cfgTofPidPtCut", 0.3f, "Minimum pt to use TOF N-sigma"};
  Configurable<bool> isUseItsPid{"isUseItsPid", false, "Use ITS PID for particle identification"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;
  o2::ft0::Geometry ft0Det;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  TrackSelection myTrackFilter;

  std::vector<double> tofNsigmaCut;
  std::vector<double> itsNsigmaCut;
  std::vector<double> tpcNsigmaCut;
  o2::aod::ITSResponse itsResponse;

  void init(InitContext&)
  {
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", cfgCcdbParam.noLaterThan.value);
    LOGF(info, "Offset for FT0A: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY(), (*offsetFT0)[0].getZ());
    LOGF(info, "Offset for FT0C: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY(), (*offsetFT0)[1].getZ());

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "Sel8");
    x->SetBinLabel(3, "ApplySameBunchPileup");
    x->SetBinLabel(4, "ApplyGoodZvtxFT0vsPV");
    x->SetBinLabel(5, "ApplyGoodITSLayersAll");
    x->SetBinLabel(6, "ApplyExtraCorrCut");
    x->SetBinLabel(7, "ApplyNoCollInTimeRangeStandard");
    x->SetBinLabel(8, "ApplyNoCollInRofStandard");
    x->SetBinLabel(9, "ApplyNoHighMultCollInPrevRof");

    myTrackFilter = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    myTrackFilter.SetPtRange(cfgtrksel.cfgPtCutMin, cfgtrksel.cfgPtCutMax);
    myTrackFilter.SetEtaRange(-cfgtrksel.cfgEtaCut, cfgtrksel.cfgEtaCut);
    myTrackFilter.SetMinNCrossedRowsTPC(cfgtrksel.minNCrossedRowsTPC);
    myTrackFilter.SetMinNClustersTPC(cfgtrksel.minTPCNClsFound);
    myTrackFilter.SetMaxDcaZ(cfgtrksel.maxDcaZ);
    myTrackFilter.SetMaxChi2PerClusterTPC(cfgtrksel.maxChi2PerClusterTPC);
    myTrackFilter.print();

    tofNsigmaCut = tofNsigmaPidCut;
    itsNsigmaCut = itsNsigmaPidCut;
    tpcNsigmaCut = tpcNsigmaPidCut;
  }

  Produces<aod::CollLRTables> collisionLRTable;
  Produces<aod::TrkLRTables> tracksLRTable;
  Produces<aod::Ft0aLRTables> ft0aLRTable;
  Produces<aod::Ft0cLRTables> ft0cLRTable;
  Produces<aod::MftTrkLRTables> mftLRTable;
  Produces<aod::MftBestTrkLRTables> mftbestLRTable;
  Produces<aod::V0TrkLRTables> v0LRTable;

  Filter fTracksEta = nabs(aod::track::eta) < cfgtrksel.cfgEtaCut;
  Filter fTracksPt = (aod::track::pt > cfgtrksel.cfgPtCutMin) && (aod::track::pt < cfgtrksel.cfgPtCutMax);

  using CollTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0Ms>;
  using TrksTable = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using MftTrkTable = aod::MFTTracks;

  void process(CollTable::iterator const& col, TrksTable const& tracks, aod::FT0s const&, MftTrkTable const& mfttracks, soa::SmallGroups<aod::BestCollisionsFwd> const& retracks, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    if (!isEventSelected(col)) {
      return;
    }

    auto multiplicity = countNTracks(tracks);
    auto centrality = selColCent(col);
    auto bc = col.bc_as<aod::BCsWithTimestamps>();

    collisionLRTable(bc.runNumber(), col.posZ(), multiplicity, centrality, bc.timestamp());

    // track loop
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack())
        continue;
      if (!myTrackFilter.IsSelected(track))
        continue;
      tracksLRTable(collisionLRTable.lastIndex(), track.pt(), track.eta(), track.phi(), aod::lrcorrtrktable::kSpCharge);
      if (getTrackPID(track) == PionTrackN)
        tracksLRTable(collisionLRTable.lastIndex(), track.pt(), track.eta(), track.phi(), aod::lrcorrtrktable::kSpPion);
      if (getTrackPID(track) == KaonTrackN)
        tracksLRTable(collisionLRTable.lastIndex(), track.pt(), track.eta(), track.phi(), aod::lrcorrtrktable::kSpKaon);
      if (getTrackPID(track) == ProtonTrackN)
        tracksLRTable(collisionLRTable.lastIndex(), track.pt(), track.eta(), track.phi(), aod::lrcorrtrktable::kSpProton);
    }

    // ft0 loop
    if (col.has_foundFT0()) {
      const auto& ft0 = col.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];
        auto phi = getPhiFT0(chanelid, 0);
        auto eta = getEtaFT0(chanelid, 0);
        ft0aLRTable(collisionLRTable.lastIndex(), chanelid, ampl, eta, phi);
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        float ampl = ft0.amplitudeC()[iCh];
        auto phi = getPhiFT0(chanelid, 1);
        auto eta = getEtaFT0(chanelid, 1);
        ft0cLRTable(collisionLRTable.lastIndex(), chanelid, ampl, eta, phi);
      }
    }

    // mft loop
    for (const auto& track : mfttracks) {
      if (!isMftTrackSelected(track))
        continue;
      auto phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      mftLRTable(collisionLRTable.lastIndex(), track.pt(), track.eta(), phi);
    }

    if (retracks.size() > 0) {
      for (const auto& retrack : retracks) {
        if (std::abs(retrack.bestDCAXY()) > cfgmfttrksel.cfigMftDcaxy) {
          continue; // does not point to PV properly
        }
        auto track = retrack.mfttrack();
        if (!isMftTrackSelected(track)) {
          continue;
        }
        auto phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        mftbestLRTable(collisionLRTable.lastIndex(), track.pt(), track.eta(), phi);
      }
    }

    // v0 loop
    for (const auto& v0 : V0s) {
      if (!isSelectV0Track(v0)) { // Quality selection for V0 prongs
        continue;
      }
      const auto& posTrack = v0.template posTrack_as<TrksTable>();
      const auto& negTrack = v0.template negTrack_as<TrksTable>();
      double massV0 = 0.0;

      // K0short
      if (isSelectK0s(col, v0)) { // candidate is K0s
        v0LRTable(collisionLRTable.lastIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                  v0.pt(), v0.eta(), v0.phi(), v0.mK0Short(), aod::lrcorrtrktable::kSpK0short);
      }

      // Lambda and Anti-Lambda
      bool lambdaTag = isSelectLambda<KindOfV0::kLambda>(col, v0);
      bool antilambdaTag = isSelectLambda<KindOfV0::kAntiLambda>(col, v0);

      // Note: candidate compatible with Lambda and Anti-Lambda hypothesis are counted twice (once for each hypothesis)
      if (lambdaTag) { // candidate is Lambda
        massV0 = v0.mLambda();
        v0LRTable(collisionLRTable.lastIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                  v0.pt(), v0.eta(), v0.phi(), massV0, aod::lrcorrtrktable::kSpLambda);
      }
      if (antilambdaTag) { // candidate is Anti-lambda
        massV0 = v0.mAntiLambda();
        v0LRTable(collisionLRTable.lastIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                  v0.pt(), v0.eta(), v0.phi(), massV0, aod::lrcorrtrktable::kSpALambda);
      } // end of Lambda and Anti-Lambda processing
    }
  } // process function

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);
    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);
    if (cfgevtsel.isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);
    if (cfgevtsel.isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);
    if (cfgevtsel.isApplyGoodITSLayersAll && !col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);
    if (cfgevtsel.isApplyExtraCorrCut && col.multNTracksPV() > cfgevtsel.npvTracksCut && col.multFT0C() < (10 * col.multNTracksPV() - cfgevtsel.ft0cCut)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);
    if (cfgevtsel.isApplyNoCollInTimeRangeStandard && !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);
    if (cfgevtsel.isApplyNoCollInRofStandard && !col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);
    if (cfgevtsel.isApplyNoHighMultCollInPrevRof && !col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 9);
    return true;
  }

  template <typename countTrk>
  int countNTracks(countTrk const& tracks)
  {
    auto nTrk = 0;
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack())
        continue;
      if (!myTrackFilter.IsSelected(track))
        continue;
      if (track.pt() < cfgtrksel.cfgPtCutMin || track.pt() > cfgtrksel.cfgPtCutMult) {
        continue;
      }
      nTrk++;
    }
    return nTrk;
  }

  template <typename CheckColCent>
  float selColCent(CheckColCent const& col)
  {
    auto cent = -1;
    if (cfgevtsel.isApplyCentFT0C) {
      cent = col.centFT0C();
    }
    if (cfgevtsel.isApplyCentFV0A) {
      cent = col.centFV0A();
    }
    if (cfgevtsel.isApplyCentFT0M) {
      cent = col.centFT0M();
    }
    return cent;
  }

  template <typename TTrack>
  int getTrackPID(TTrack const& track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    std::array<float, 3> nSigmaToUse = isUseItsPid ? nSigmaITS : nSigmaTPC;            // Choose which nSigma to use: TPC or ITS
    std::vector<double> detectorNsigmaCut = isUseItsPid ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS
    int pid = -1;
    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[0] < detectorNsigmaCut[0] && nSigmaToUse[0] > detectorNsigmaCut[0 + 3];
    bool isDetectedKaon = nSigmaToUse[1] < detectorNsigmaCut[1] && nSigmaToUse[1] > detectorNsigmaCut[1 + 3];
    bool isDetectedProton = nSigmaToUse[2] < detectorNsigmaCut[2] && nSigmaToUse[2] > detectorNsigmaCut[2 + 3];

    bool isTofPion = nSigmaTOF[0] < tofNsigmaCut[0] && nSigmaTOF[0] > tofNsigmaCut[0 + 3];
    bool isTofKaon = nSigmaTOF[1] < tofNsigmaCut[1] && nSigmaTOF[1] > tofNsigmaCut[1 + 3];
    bool isTofProton = nSigmaTOF[2] < tofNsigmaCut[2] && nSigmaTOF[2] > tofNsigmaCut[2 + 3];

    if (track.pt() > cfgTofPidPtCut && !track.hasTOF()) {
      return 0;
    } else if (track.pt() > cfgTofPidPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return 0; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PIONS;
    } else if (isKaon) {
      pid = KAONS;
    } else if (isProton) {
      pid = PROTONS;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  double getPhiFT0(uint chno, int i)
  {
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + (*offsetFT0)[i].getX(), chPos.Y() + (*offsetFT0)[i].getY());
  }

  double getEtaFT0(uint chno, int i)
  {
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    if (chno >= KminFt0cCell)
      z = -z;
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  template <typename CheckMftTrack>
  bool isMftTrackSelected(CheckMftTrack const& track)
  {
    if (track.nClusters() < cfgmfttrksel.cfigMftCluster) {
      return false;
    }
    if (track.eta() > cfgmfttrksel.cfigMftEtaMax || track.eta() < cfgmfttrksel.cfigMftEtaMin) {
      return false;
    }
    if (cfgmfttrksel.useMftPtCut && (track.pt() < cfgmfttrksel.cfgMftPtCutMin || track.pt() > cfgmfttrksel.cfgMftPtCutMax)) {
      return false;
    }
    return true;
  }

  template <typename T1>
  bool isSelectV0Track(const T1& v0)
  {
    const auto& posTrack = v0.template posTrack_as<TrksTable>();
    const auto& negTrack = v0.template negTrack_as<TrksTable>();

    if (!posTrack.hasTPC() || !negTrack.hasTPC()) {
      return false;
    }
    if (posTrack.tpcNClsCrossedRows() < cfgv0trksel.minTPCcrossedrows || negTrack.tpcNClsCrossedRows() < cfgv0trksel.minTPCcrossedrows) {
      return false;
    }
    if (posTrack.tpcCrossedRowsOverFindableCls() < cfgv0trksel.minTPCcrossedrowsoverfindcls || negTrack.tpcCrossedRowsOverFindableCls() < cfgv0trksel.minTPCcrossedrowsoverfindcls) {
      return false;
    }
    if (std::abs(v0.positiveeta()) > cfgv0trksel.v0etaCut || std::abs(v0.negativeeta()) > cfgv0trksel.v0etaCut) {
      return false;
    }
    if (v0.positivept() < cfgv0trksel.v0ptCut || v0.negativept() < cfgv0trksel.v0ptCut) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0candidate>
  bool isSelectK0s(Collision const& col, const V0candidate& v0)
  {
    const auto& posTrack = v0.template posTrack_as<TrksTable>();
    const auto& negTrack = v0.template negTrack_as<TrksTable>();

    float ctauK0s = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * o2::constants::physics::MassK0;

    if (v0.mK0Short() < cfgv0trksel.minK0sMass || v0.mK0Short() > cfgv0trksel.maxK0sMass) {
      return false;
    }
    if ((v0.qtarm() / std::abs(v0.alpha())) < cfgv0trksel.minqtArmenterosForK0s) {
      return false;
    }
    if (v0.v0radius() > cfgv0trksel.maxK0sRadius || v0.v0radius() < cfgv0trksel.minK0sRadius) {
      return false;
    }
    if (v0.v0cosPA() < cfgv0trksel.minK0sCosPa) {
      return false;
    }
    if (v0.dcaV0daughters() > cfgv0trksel.maxDcaV0DauK0s) {
      return false;
    }
    if (std::abs(ctauK0s) > cfgv0trksel.maxK0sLifeTime) {
      return false;
    }
    if (((std::abs(posTrack.tpcNSigmaPi()) > cfgv0trksel.daughPIDCuts) || (std::abs(negTrack.tpcNSigmaPi()) > cfgv0trksel.daughPIDCuts))) {
      return false;
    }
    if ((std::abs(v0.dcapostopv()) < cfgv0trksel.minV0DcaPiK0s || std::abs(v0.dcanegtopv()) < cfgv0trksel.minV0DcaPiK0s)) {
      return false;
    }
    return true;
  }

  template <KindOfV0 pid, typename Collision, typename V0candidate>
  bool isSelectLambda(Collision const& col, const V0candidate& v0)
  {
    const auto& posTrack = v0.template posTrack_as<TrksTable>();
    const auto& negTrack = v0.template negTrack_as<TrksTable>();
    float ctauLambda = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * o2::constants::physics::MassLambda;
    if ((v0.mLambda() < cfgv0trksel.minLambdaMass || v0.mLambda() > cfgv0trksel.maxLambdaMass) &&
        (v0.mAntiLambda() < cfgv0trksel.minLambdaMass || v0.mAntiLambda() > cfgv0trksel.maxLambdaMass)) {
      return false;
    }
    if (v0.v0radius() > cfgv0trksel.maxLambdaRadius || v0.v0radius() < cfgv0trksel.minLambdaRadius) {
      return false;
    }
    if (v0.v0cosPA() < cfgv0trksel.minLambdaCosPa) {
      return false;
    }
    if (v0.dcaV0daughters() > cfgv0trksel.maxDcaV0DauLambda) {
      return false;
    }
    if (pid == KindOfV0::kLambda && (std::abs(v0.dcapostopv()) < cfgv0trksel.minV0DcaPr || std::abs(v0.dcanegtopv()) < cfgv0trksel.minV0DcaPiLambda)) {
      return false;
    }
    if (pid == KindOfV0::kAntiLambda && (std::abs(v0.dcapostopv()) < cfgv0trksel.minV0DcaPiLambda || std::abs(v0.dcanegtopv()) < cfgv0trksel.minV0DcaPr)) {
      return false;
    }
    if (pid == KindOfV0::kLambda && ((std::abs(posTrack.tpcNSigmaPr()) > cfgv0trksel.daughPIDCuts) || (std::abs(negTrack.tpcNSigmaPi()) > cfgv0trksel.daughPIDCuts))) {
      return false;
    }
    if (pid == KindOfV0::kAntiLambda && ((std::abs(posTrack.tpcNSigmaPi()) > cfgv0trksel.daughPIDCuts) || (std::abs(negTrack.tpcNSigmaPr()) > cfgv0trksel.daughPIDCuts))) {
      return false;
    }
    if (std::abs(ctauLambda) > cfgv0trksel.maxLambdaLifeTime) {
      return false;
    }
    return true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LongrangeMaker>(cfgc)};
}
