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
/// \file longrangeCorrelation.cxx
///
/// \brief task for long range correlation analysis
/// \author Abhi Modak (abhi.modak@cern.ch) and Debojit Sarkar (debojit.sarkar@cern.ch)
/// \since April 22, 2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TComplex.h>
#include <TH1F.h>
#include <TMath.h>

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

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDcaxyOnly =
  TrackSelectionFlags::kDCAxy;

AxisSpec axisEvent{10, 0.5, 9.5, "#Event", "EventAxis"};

struct LongrangeCorrelation {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  std::vector<o2::detectors::AlignParam>* offsetFV0;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Vertex Z range to consider"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 1.0f, "Eta range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.2f, "minimum accepted track pT"};
  Configurable<float> cfgPtCutMax{"cfgPtCutMax", 10.0f, "maximum accepted track pT"};
  Configurable<int> mixingParameter{"mixingParameter", 5, "how many events are mixed"};
  Configurable<int> cfgMinMult{"cfgMinMult", 0, "Minimum multiplicity for collision"};
  Configurable<int> cfgMaxMult{"cfgMaxMult", 10, "Maximum multiplicity for collision"};
  Configurable<float> cfigMftEtaMax{"cfigMftEtaMax", -2.5f, "Maximum MFT eta cut"};
  Configurable<float> cfigMftEtaMin{"cfigMftEtaMin", -3.6f, "Minimum MFT eta cut"};
  Configurable<float> cfigMftDcaxy{"cfigMftDcaxy", 2.0f, "cut on DCA xy for MFT tracks"};
  Configurable<int> cfigMftCluster{"cfigMftCluster", 5, "cut on MFT Cluster"};
  Configurable<double> cfgSampleSize{"cfgSampleSize", 10, "Sample size for mixed event"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", false, "Enable SameBunchPileup cut"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -6, -2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultME{"axisMultME", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20, 20}, "vertex axis"};
  ConfigurableAxis axisPhi{"axisPhi", {96, 0, TwoPI}, "#phi axis"};
  ConfigurableAxis axisEtaTrig{"axisEtaTrig", {40, -1., 1.}, "#eta trig axis"};
  ConfigurableAxis axisEtaAssoc{"axisEtaAssoc", {96, 3.5, 4.9}, "#eta assoc axis"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis amplitudeFt0a{"amplitudeFt0a", {5000, 0, 10000}, "FT0A amplitude"};
  ConfigurableAxis channelFt0aAxis{"channelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};

  using CollTable = soa::Join<aod::Collisions, aod::EvSels>;
  using TrksTable = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  using MftTrkTable = soa::Filtered<aod::MFTTracks>;
  Preslice<TrksTable> perColGlobal = aod::track::collisionId;
  Preslice<MftTrkTable> perColMft = aod::fwdtrack::collisionId;

  OutputObj<CorrelationContainer> sameFt0aGlobal{Form("sameEventFt0aGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixedFt0aGlobal{Form("mixedEventFt0aGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> sameFt0cGlobal{Form("sameEventFt0cGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixedFt0cGlobal{Form("mixedEventFt0cGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> sameMftGlobal{Form("sameEventMftGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixedMftGlobal{Form("mixedEventMftGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> sameFv0Global{Form("sameEventFv0Global_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixedFv0Global{Form("mixedEventFv0Global_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> sameFv0Mft{Form("sameEventFv0Mft_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixedFv0Mft{Form("mixedEventFv0Mft_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};

  void init(InitContext const&)
  {
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", cfgCcdbParam.noLaterThan.value);
    offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", cfgCcdbParam.noLaterThan.value);
    LOGF(info, "Offset for FT0A: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY(), (*offsetFT0)[0].getZ());
    LOGF(info, "Offset for FT0C: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY(), (*offsetFT0)[1].getZ());
    LOGF(info, "Offset for FV0-left: x = %.3f y = %.3f\n", (*offsetFV0)[0].getX(), (*offsetFV0)[0].getY());
    LOGF(info, "Offset for FV0-right: x = %.3f y = %.3f\n", (*offsetFV0)[1].getX(), (*offsetFV0)[1].getY());

    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVtxZ, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {{axisVertexEfficiency, "z-vtx (cm)"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisEtaEfficiency, "#eta"}};

    std::vector<AxisSpec> userAxis;

    if (doprocessEventStat) {
      histos.add("QA/EventHist", "events", kTH1F, {axisEvent}, false);
      histos.add("QA/VtxZHist", "v_{z} (cm)", kTH1F, {axisVtxZ}, false);

      auto hstat = histos.get<TH1>(HIST("QA/EventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All events");
      x->SetBinLabel(2, "sel8");
      x->SetBinLabel(3, "kNoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
      x->SetBinLabel(4, "|vz|<10");
    }

    histos.add("Ft0aGlobal/SE/hMult", "", kTH1D, {axisMultiplicity});
    histos.add("Ft0aGlobal/SE/Trig_etavsphi", "", kTH2D, {axisPhi, axisEtaTrig});
    histos.add("Ft0aGlobal/SE/Trig_eta", "", kTH1D, {axisEtaTrig});
    histos.add("Ft0aGlobal/SE/Trig_phi", "", kTH1D, {axisPhi});
    histos.add("Ft0aGlobal/SE/Trig_pt", "", kTH1D, {axisPtTrigger});
    histos.add("Ft0aGlobal/SE/hMult_used", "", kTH1F, {axisMultiplicity});
    histos.add("Ft0aGlobal/SE/Trig_hist", "", kTHnSparseF, {axisSample, axisVtxZ, axisPtTrigger});
    histos.add("Ft0aGlobal/SE/Assoc_amp", "", kTH2D, {channelFt0aAxis, amplitudeFt0a});
    histos.add("Ft0aGlobal/SE/Assoc_eta", "", kTH1D, {axisEtaAssoc});
    histos.add("Ft0aGlobal/SE/Assoc_phi", "", kTH1D, {axisPhi});
    histos.add("Ft0aGlobal/SE/Assoc_etavsphi", "", kTH2D, {axisPhi, axisEtaAssoc});
    histos.add("Ft0aGlobal/SE/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});

    histos.add("Ft0aGlobal/ME/hMult", "", kTH1D, {axisMultiplicity});
    histos.add("Ft0aGlobal/ME/Trig_etavsphi", "", kTH2D, {axisPhi, axisEtaTrig});
    histos.add("Ft0aGlobal/ME/Trig_eta", "", kTH1D, {axisEtaTrig});
    histos.add("Ft0aGlobal/ME/Trig_phi", "", kTH1D, {axisPhi});
    histos.add("Ft0aGlobal/ME/Trig_pt", "", kTH1D, {axisPtTrigger});
    histos.add("Ft0aGlobal/ME/Assoc_amp", "", kTH2D, {channelFt0aAxis, amplitudeFt0a});
    histos.add("Ft0aGlobal/ME/Assoc_eta", "", kTH1D, {axisEtaAssoc});
    histos.add("Ft0aGlobal/ME/Assoc_phi", "", kTH1D, {axisPhi});
    histos.add("Ft0aGlobal/ME/Assoc_etavsphi", "", kTH2D, {axisPhi, axisEtaAssoc});
    histos.add("Ft0aGlobal/ME/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});

    if (doprocessFt0aGlobalSE || doprocessFt0aGlobalME) {
      sameFt0aGlobal.setObject(new CorrelationContainer(Form("sameEventFt0aGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEventFt0aGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
      mixedFt0aGlobal.setObject(new CorrelationContainer(Form("mixedEventFt0aGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEventFt0aGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    }

    if (doprocessFt0cGlobalSE || doprocessFt0cGlobalME) {
      histos.addClone("Ft0aGlobal/SE/", "Ft0cGlobal/SE/");
      histos.addClone("Ft0aGlobal/ME/", "Ft0cGlobal/ME/");
      sameFt0cGlobal.setObject(new CorrelationContainer(Form("sameEventFt0cGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEventFt0cGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
      mixedFt0cGlobal.setObject(new CorrelationContainer(Form("mixedEventFt0cGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEventFt0cGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    }

    if (doprocessMftGlobalSE || doprocessMftGlobalME) {
      histos.addClone("Ft0aGlobal/SE/", "MftGlobal/SE/");
      histos.addClone("Ft0aGlobal/ME/", "MftGlobal/ME/");
      sameMftGlobal.setObject(new CorrelationContainer(Form("sameEventMftGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEventMftGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
      mixedMftGlobal.setObject(new CorrelationContainer(Form("mixedEventMftGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEventMftGlobal_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    }

    if (doprocessFv0GlobalSE || doprocessFv0GlobalME) {
      histos.addClone("Ft0aGlobal/SE/", "Fv0Global/SE/");
      histos.addClone("Ft0aGlobal/ME/", "Fv0Global/ME/");
      sameFv0Global.setObject(new CorrelationContainer(Form("sameEventFv0Global_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEventFv0Global_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
      mixedFv0Global.setObject(new CorrelationContainer(Form("mixedEventFv0Global_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEventFv0Global_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    }

    if (doprocessFv0MftSE || doprocessFv0MftME) {
      histos.addClone("Ft0aGlobal/SE/", "Fv0Mft/SE/");
      histos.addClone("Ft0aGlobal/ME/", "Fv0Mft/ME/");
      sameFv0Mft.setObject(new CorrelationContainer(Form("sameEventFv0Mft_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEventFv0Mft_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
      mixedFv0Mft.setObject(new CorrelationContainer(Form("mixedEventFv0Mft_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEventFv0Mft_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    }
  }

  Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));
  Filter fTracksEta = nabs(aod::track::eta) < cfgEtaCut;
  Filter fTracksPt = (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax);

  Filter fMftTrackEta = (aod::fwdtrack::eta < cfigMftEtaMax) && (aod::fwdtrack::eta > cfigMftEtaMin);
  Filter fMftTrackColID = (aod::fwdtrack::bestCollisionId >= 0);
  Filter fMftTrackDca = (nabs(aod::fwdtrack::bestDCAXY) < cfigMftDcaxy);

  double getPhiFT0(int chno, double offsetX, double offsetY)
  {
    o2::ft0::Geometry ft0Det;
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + offsetX, chPos.Y() + offsetY);
  }

  double getPhiFV0(int chno)
  {
    o2::fv0::Geometry fv0Det;
    int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);
    float offsetX, offsetY;
    if (isChnoInLeft) {
      offsetX = (*offsetFV0)[0].getX();
      offsetY = (*offsetFV0)[0].getY();
    } else {
      offsetX = (*offsetFV0)[1].getX();
      offsetY = (*offsetFV0)[1].getY();
    }

    auto chPos = fv0Det.getReadoutCenter(chno);
    return RecoDecay::phi(chPos.x + offsetX, chPos.y + offsetY);
  }

  double getEtaFT0(int chno, double offsetX, double offsetY, double offsetZ)
  {
    o2::ft0::Geometry ft0Det;
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + offsetX;
    auto y = chPos.Y() + offsetY;
    auto z = chPos.Z() + offsetZ;
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  double getEtaFV0(int chno)
  {
    o2::fv0::Geometry fv0Det;
    int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);
    float offsetX, offsetY, offsetZ;
    if (isChnoInLeft) {
      offsetX = (*offsetFV0)[0].getX();
      offsetY = (*offsetFV0)[0].getY();
      offsetZ = (*offsetFV0)[0].getZ();
    } else {
      offsetX = (*offsetFV0)[1].getX();
      offsetY = (*offsetFV0)[1].getY();
      offsetZ = (*offsetFV0)[1].getZ();
    }

    auto chPos = fv0Det.getReadoutCenter(chno);
    auto x = chPos.x + offsetX;
    auto y = chPos.y + offsetY;
    auto z = chPos.z + offsetZ;
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    if (!col.sel8()) {
      return false;
    }
    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (std::abs(col.posZ()) >= cfgVtxCut) {
      return false;
    }
    return true;
  }

  template <typename CheckMftTrack>
  bool isMftTrackSelected(CheckMftTrack const& track)
  {
    if (track.nClusters() < cfigMftCluster) {
      return false;
    }
    return true;
  }

  template <int mode, typename TTracks>
  void fillYield(TTracks tracks, bool mixing)
  {
    static constexpr std::string_view SubDirSE[] = {"Ft0aGlobal/SE/", "Ft0cGlobal/SE/", "Fv0Global/SE/",
                                                    "MftGlobal/SE/", "Fv0Mft/SE/"};
    static constexpr std::string_view SubDirME[] = {"Ft0aGlobal/ME/", "Ft0cGlobal/ME/", "Fv0Global/ME/",
                                                    "MftGlobal/ME/", "Fv0Mft/ME/"};

    if (mixing) {
      histos.fill(HIST(SubDirME[mode]) + HIST("hMult"), tracks.size());
      for (auto const& triggerTrack : tracks) {
        histos.fill(HIST(SubDirME[mode]) + HIST("Trig_etavsphi"), triggerTrack.phi(), triggerTrack.eta());
        histos.fill(HIST(SubDirME[mode]) + HIST("Trig_eta"), triggerTrack.eta());
        histos.fill(HIST(SubDirME[mode]) + HIST("Trig_phi"), triggerTrack.phi());
        histos.fill(HIST(SubDirME[mode]) + HIST("Trig_pt"), triggerTrack.pt());
      }
    } else {
      histos.fill(HIST(SubDirSE[mode]) + HIST("hMult"), tracks.size());
      for (auto const& triggerTrack : tracks) {
        histos.fill(HIST(SubDirSE[mode]) + HIST("Trig_etavsphi"), triggerTrack.phi(), triggerTrack.eta());
        histos.fill(HIST(SubDirSE[mode]) + HIST("Trig_eta"), triggerTrack.eta());
        histos.fill(HIST(SubDirSE[mode]) + HIST("Trig_phi"), triggerTrack.phi());
        histos.fill(HIST(SubDirSE[mode]) + HIST("Trig_pt"), triggerTrack.pt());
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFT0s>
  void fillCorrFt0aGlobal(TTarget target, TTriggers const& triggers, TFT0s const& ft0, bool mixing, float vz)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    if (!mixing)
      histos.fill(HIST("Ft0aGlobal/SE/hMult_used"), triggers.size());
    for (auto const& triggerTrack : triggers) {
      if (!mixing)
        histos.fill(HIST("Ft0aGlobal/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt());

      auto offsetX = (*offsetFT0)[0].getX();
      auto offsetY = (*offsetFT0)[0].getY();
      auto offsetZ = (*offsetFT0)[0].getZ();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];
        if (ampl <= 0)
          continue;
        if (mixing)
          histos.fill(HIST("Ft0aGlobal/ME/Assoc_amp"), chanelid, ampl);
        else
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_amp"), chanelid, ampl);

        auto phi = getPhiFT0(chanelid, offsetX, offsetY);
        auto eta = getEtaFT0(chanelid, offsetX, offsetY, offsetZ);

        if (mixing) {
          histos.fill(HIST("Ft0aGlobal/ME/Assoc_eta"), eta);
          histos.fill(HIST("Ft0aGlobal/ME/Assoc_phi"), phi);
          histos.fill(HIST("Ft0aGlobal/ME/Assoc_etavsphi"), phi, eta);
        } else {
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_eta"), eta);
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_phi"), phi);
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_etavsphi"), phi, eta);
        }
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Ft0aGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("Ft0aGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrFt0aGlobal

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFT0s>
  void fillCorrFt0cGlobal(TTarget target, TTriggers const& triggers, TFT0s const& ft0, bool mixing, float vz)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    if (!mixing)
      histos.fill(HIST("Ft0cGlobal/SE/hMult_used"), triggers.size());
    for (auto const& triggerTrack : triggers) {
      if (!mixing)
        histos.fill(HIST("Ft0cGlobal/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt());

      auto offsetX = (*offsetFT0)[1].getX();
      auto offsetY = (*offsetFT0)[1].getY();
      auto offsetZ = (*offsetFT0)[1].getZ();
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        float ampl = ft0.amplitudeC()[iCh];
        if (ampl <= 0)
          continue;
        if (mixing)
          histos.fill(HIST("Ft0cGlobal/ME/Assoc_amp"), chanelid, ampl);
        else
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_amp"), chanelid, ampl);

        auto phi = getPhiFT0(chanelid, offsetX, offsetY);
        auto eta = getEtaFT0(chanelid, offsetX, offsetY, offsetZ);

        if (mixing) {
          histos.fill(HIST("Ft0cGlobal/ME/Assoc_eta"), eta);
          histos.fill(HIST("Ft0cGlobal/ME/Assoc_phi"), phi);
          histos.fill(HIST("Ft0cGlobal/ME/Assoc_etavsphi"), phi, eta);
        } else {
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_eta"), eta);
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_phi"), phi);
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_etavsphi"), phi, eta);
        }
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Ft0cGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("Ft0cGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrFt0cGlobal

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TMFTs>
  void fillCorrMftGlobal(TTarget target, TTriggers const& triggers, TMFTs const& mft, bool mixing, float vz)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    if (!mixing)
      histos.fill(HIST("MftGlobal/SE/hMult_used"), triggers.size());
    for (auto const& triggerTrack : triggers) {
      if (!mixing)
        histos.fill(HIST("MftGlobal/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt());

      for (auto const& assoTrack : mft) {
        if (!isMftTrackSelected(assoTrack)) {
          continue;
        }
        auto phi = assoTrack.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (mixing) {
          histos.fill(HIST("MftGlobal/ME/Assoc_eta"), assoTrack.eta());
          histos.fill(HIST("MftGlobal/ME/Assoc_phi"), phi);
          histos.fill(HIST("MftGlobal/ME/Assoc_etavsphi"), phi, assoTrack.eta());
        } else {
          histos.fill(HIST("MftGlobal/SE/Assoc_eta"), assoTrack.eta());
          histos.fill(HIST("MftGlobal/SE/Assoc_phi"), phi);
          histos.fill(HIST("MftGlobal/SE/Assoc_etavsphi"), phi, assoTrack.eta());
        }
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - assoTrack.eta();
        if (mixing)
          histos.fill(HIST("MftGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("MftGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), assoTrack.pt(), deltaPhi, deltaEta);
      } // associated mft tracks
    } // trigger tracks
  } // fillCorrMftGlobal

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFV0s>
  void fillCorrFv0Global(TTarget target, TTriggers const& triggers, TFV0s const& fv0, bool mixing, float vz)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    if (!mixing)
      histos.fill(HIST("Fv0Global/SE/hMult_used"), triggers.size());
    for (auto const& triggerTrack : triggers) {
      if (!mixing)
        histos.fill(HIST("Fv0Global/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt());

      for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
        auto chanelid = fv0.channel()[iCh];
        float ampl = fv0.amplitude()[iCh];
        if (ampl <= 0)
          continue;
        if (mixing)
          histos.fill(HIST("Fv0Global/ME/Assoc_amp"), chanelid, ampl);
        else
          histos.fill(HIST("Fv0Global/SE/Assoc_amp"), chanelid, ampl);

        auto phi = getPhiFV0(chanelid);
        auto eta = getEtaFV0(chanelid);

        if (mixing) {
          histos.fill(HIST("Fv0Global/ME/Assoc_eta"), eta);
          histos.fill(HIST("Fv0Global/ME/Assoc_phi"), phi);
          histos.fill(HIST("Fv0Global/ME/Assoc_etavsphi"), phi, eta);
        } else {
          histos.fill(HIST("Fv0Global/SE/Assoc_eta"), eta);
          histos.fill(HIST("Fv0Global/SE/Assoc_phi"), phi);
          histos.fill(HIST("Fv0Global/SE/Assoc_etavsphi"), phi, eta);
        }
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Fv0Global/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("Fv0Global/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta);
      } // associated fv0 tracks
    } // trigger tracks
  } // fillCorrFv0Global

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks, typename TTriggers, typename TFV0s>
  void fillCorrFv0Mft(TTarget target, TTracks const& tracks, TTriggers const& triggers, TFV0s const& fv0, bool mixing, float vz)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    if (!mixing)
      histos.fill(HIST("Fv0Mft/SE/hMult_used"), tracks.size());
    for (auto const& triggerTrack : triggers) {
      if (!isMftTrackSelected(triggerTrack)) {
        continue;
      }
      if (!mixing)
        histos.fill(HIST("Fv0Mft/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt());

      auto trigphi = triggerTrack.phi();
      o2::math_utils::bringTo02Pi(trigphi);

      for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
        auto chanelid = fv0.channel()[iCh];
        float ampl = fv0.amplitude()[iCh];
        if (ampl <= 0)
          continue;
        if (mixing)
          histos.fill(HIST("Fv0Mft/ME/Assoc_amp"), chanelid, ampl);
        else
          histos.fill(HIST("Fv0Mft/SE/Assoc_amp"), chanelid, ampl);

        auto phi = getPhiFV0(chanelid);
        auto eta = getEtaFV0(chanelid);

        if (mixing) {
          histos.fill(HIST("Fv0Mft/ME/Assoc_eta"), eta);
          histos.fill(HIST("Fv0Mft/ME/Assoc_phi"), phi);
          histos.fill(HIST("Fv0Mft/ME/Assoc_etavsphi"), phi, eta);
        } else {
          histos.fill(HIST("Fv0Mft/SE/Assoc_eta"), eta);
          histos.fill(HIST("Fv0Mft/SE/Assoc_phi"), phi);
          histos.fill(HIST("Fv0Mft/SE/Assoc_etavsphi"), phi, eta);
        }

        float deltaPhi = RecoDecay::constrainAngle(trigphi - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Fv0Mft/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("Fv0Mft/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta);
      } // associated fv0 tracks
    } // trigger tracks
  } // fillCorrFv0Mft

  void processEventStat(CollTable::iterator const& col)
  {
    histos.fill(HIST("QA/EventHist"), 1);
    if (!col.sel8()) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 2);
    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 3);
    if (std::abs(col.posZ()) >= cfgVtxCut) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 4);
    histos.fill(HIST("QA/VtxZHist"), col.posZ());
  }

  void processFt0aGlobalSE(CollTable::iterator const& col, aod::FT0s const&, TrksTable const& tracks)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      fillYield<0>(tracks, false);
      const auto& ft0 = col.foundFT0();
      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }
      fillCorrFt0aGlobal<CorrelationContainer::kCFStepReconstructed>(sameFt0aGlobal, tracks, ft0, false, col.posZ());
    }
  } // same event

  void processFt0cGlobalSE(CollTable::iterator const& col, aod::FT0s const&, TrksTable const& tracks)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      fillYield<1>(tracks, false);
      const auto& ft0 = col.foundFT0();
      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }
      fillCorrFt0cGlobal<CorrelationContainer::kCFStepReconstructed>(sameFt0cGlobal, tracks, ft0, false, col.posZ());
    }
  } // same event

  void processMftGlobalSE(CollTable::iterator const& col, MftTrkTable const& mfttracks, TrksTable const& tracks)
  {
    if (!isEventSelected(col)) {
      return;
    }
    fillYield<3>(tracks, false);
    if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
      return;
    }
    fillCorrMftGlobal<CorrelationContainer::kCFStepReconstructed>(sameMftGlobal, tracks, mfttracks, false, col.posZ());
  } // same event

  void processFv0GlobalSE(CollTable::iterator const& col, aod::FV0As const&, TrksTable const& tracks)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFV0()) {
      fillYield<2>(tracks, false);
      const auto& fv0 = col.foundFV0();
      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }
      fillCorrFv0Global<CorrelationContainer::kCFStepReconstructed>(sameFv0Global, tracks, fv0, false, col.posZ());
    }
  } // same event

  void processFv0MftSE(CollTable::iterator const& col, aod::FV0As const&, TrksTable const& tracks, MftTrkTable const& mfttracks)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFV0()) {
      fillYield<4>(mfttracks, false);
      const auto& fv0 = col.foundFV0();
      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }
      fillCorrFv0Mft<CorrelationContainer::kCFStepReconstructed>(sameFv0Mft, tracks, mfttracks, fv0, false, col.posZ());
    }
  } // same event

  void processFt0aGlobalME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        fillYield<0>(slicedTriggerTracks, true);
        const auto& ft0 = col2.foundFT0();
        if (slicedTriggerTracks.size() < cfgMinMult || slicedTriggerTracks.size() >= cfgMaxMult) {
          continue;
        }
        fillCorrFt0aGlobal<CorrelationContainer::kCFStepReconstructed>(mixedFt0aGlobal, slicedTriggerTracks, ft0, true, col1.posZ());
      }
    }
  } // mixed event

  void processFt0cGlobalME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        fillYield<1>(slicedTriggerTracks, true);
        const auto& ft0 = col2.foundFT0();
        if (slicedTriggerTracks.size() < cfgMinMult || slicedTriggerTracks.size() >= cfgMaxMult) {
          continue;
        }
        fillCorrFt0cGlobal<CorrelationContainer::kCFStepReconstructed>(mixedFt0cGlobal, slicedTriggerTracks, ft0, true, col1.posZ());
      }
    }
  } // mixed event

  void processMftGlobalME(CollTable const& col, MftTrkTable const& mfttracks, TrksTable const& tracks)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    auto tracksTuple = std::make_tuple(tracks, mfttracks);
    Pair<CollTable, TrksTable, MftTrkTable, MixedBinning> pairs{binningOnVtxAndMult, mixingParameter, -1, col, tracksTuple, &cache};
    for (auto const& [col1, tracks1, col2, tracks2] : pairs) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if ((tracks1.size() < cfgMinMult || tracks1.size() >= cfgMaxMult)) {
        continue;
      }
      fillCorrMftGlobal<CorrelationContainer::kCFStepReconstructed>(mixedMftGlobal, tracks1, tracks2, true, col1.posZ());
    }
  } // mixed event

  void processFv0GlobalME(CollTable const& col, aod::FV0As const&, TrksTable const& tracks)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFV0() && col2.has_foundFV0()) {
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        fillYield<2>(slicedTriggerTracks, true);
        const auto& fv0 = col2.foundFV0();
        if (slicedTriggerTracks.size() < cfgMinMult || slicedTriggerTracks.size() >= cfgMaxMult) {
          continue;
        }
        fillCorrFv0Global<CorrelationContainer::kCFStepReconstructed>(mixedFv0Global, slicedTriggerTracks, fv0, true, col1.posZ());
      }
    }
  } // mixed event

  void processFv0MftME(CollTable const& col, aod::FV0As const&, TrksTable const& tracks, MftTrkTable const& mfttracks)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFV0() && col2.has_foundFV0()) {
        auto slicedGlobalTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        auto slicedTriggerMftTracks = mfttracks.sliceBy(perColMft, col1.globalIndex());
        fillYield<4>(slicedTriggerMftTracks, true);
        const auto& fv0 = col2.foundFV0();
        if (slicedGlobalTracks.size() < cfgMinMult || slicedGlobalTracks.size() >= cfgMaxMult) {
          continue;
        }
        fillCorrFv0Mft<CorrelationContainer::kCFStepReconstructed>(mixedFv0Mft, slicedGlobalTracks, slicedTriggerMftTracks, fv0, true, col1.posZ());
      }
    }
  } // mixed event

  PROCESS_SWITCH(LongrangeCorrelation, processEventStat, "event stat", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aGlobalSE, "same event FT0a vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aGlobalME, "mixed event FT0a vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0cGlobalSE, "same event FT0c vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0cGlobalME, "mixed event FT0c vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processMftGlobalSE, "same event MFT vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processMftGlobalME, "mixed event MFT vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFv0GlobalSE, "same event FV0 vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFv0GlobalME, "mixed event FV0 vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFv0MftSE, "same event FV0 vs MFT", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFv0MftME, "mixed event FV0 vs MFT", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LongrangeCorrelation>(cfgc)};
}
