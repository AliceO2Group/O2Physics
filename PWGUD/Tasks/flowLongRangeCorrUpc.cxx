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

/// \file flowLongRangeCorrUpc.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch), Yongxi Du (yongxi.du@cern.ch)
/// \since  Apr/7/2026
/// \brief task for TPC-FT0 correlations in UPC

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGUD/Core/UpcService.h"
#include "PWGUD/Core/UDHelpers.h"     // udhelpers::Bits256, makeBits256, testBit, getPhiEtaFromFitBit
#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <FT0Base/Geometry.h> // o2::ft0::Geometry
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TRandom3.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowLongRangeCorrUpc {
  // event selection
  O2_DEFINE_CONFIGURABLE(cfgGapSide, int, 1, "choose one side 0:A; 1:C")
  O2_DEFINE_CONFIGURABLE(cfgGapSideMerge, bool, true, "merge A and C side")
  // flags that already enabled in sgcandproducer: kNoTimeFrameBorder, kNoITSROFrameBorder, kNoSameBunchPileup, kIsGoodZvtxFT0vsPV, kIsVertexITSTPC, isGoodRCTCollision
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMinMult, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgMaxMult, int, 5, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")

  // track selection
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutRefMin, float, 0.2f, "minimum accepted reference track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutRefMax, float, 3.0f, "minimum accepted reference track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.9f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgGlobalTrack, bool, true, "use global track")
  O2_DEFINE_CONFIGURABLE(cfgDcaXY, bool, true, "enable dcaxy cut")
  O2_DEFINE_CONFIGURABLE(cfgDcaZ, bool, false, "enable dcaz cut")
  O2_DEFINE_CONFIGURABLE(cfgDcaZCut, float, 10.0, "dcaz cut")
  O2_DEFINE_CONFIGURABLE(cfgMaxTPCChi2NCl, int, 4, "max chi2 for tpc fit")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum number of crossed TPC Rows")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum number of found TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum number of ITS clusters")

  // user define
  O2_DEFINE_CONFIGURABLE(cfgUseNchCorrected, bool, true, "use efficiency corrected tracks")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, true, "enable trigger pT < associated pT cut")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, true, "enable trigger pT < associated pT cut in mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgDrawEtaPhiDis, bool, false, "draw eta-phi distribution for detectors in used")

  SliceCache cache;

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {36, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {20, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0a{"axisDeltaEtaTpcFt0a", {32, -5.8, -2.6}, "delta eta axis, -5.8~-2.6 for TPC-FT0A,"};
  ConfigurableAxis axisDeltaEtaTpcFt0c{"axisDeltaEtaTpcFt0c", {32, 1.2, 4.2}, "delta eta axis, 1.2~4.2 for TPC-FT0C"};
  ConfigurableAxis axisDeltaEtaFt0aFt0c{"axisDeltaEtaFt0aFt0c", {32, 5.5, 8.5}, "delta eta axis"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30}, "multiplicity axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, 0, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30}, "multiplicity / centrality axis for mixed event histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisThresholdFt0a{"axisThresholdFt0a", {5, 0, 5}, "FT0A amplitude"};
  AxisSpec axisChID = {220, 0, 220};

  // make the filters and cuts.
  Filter trackFilter = (aod::udtrack::isPVContributor == true);
  Filter collisionFilter = (cfgGapSideMerge == true)
  ? ((aod::udcollision::gapSide == (uint8_t)0 || aod::udcollision::gapSide == (uint8_t)1) &&
     (aod::upcservice::truegapside == 0 || aod::upcservice::truegapside == 1))
  : ((aod::udcollision::gapSide == (uint8_t)cfgGapSide) &&
     (aod::upcservice::truegapside == cfgGapSide));

  using UDTracksFull = soa::Filtered<soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>>;
  using UDCollisionsFull = soa::Filtered<soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::Truegapside, aod::UDCollisionSelExtras>>;

  // FT0 geometry
  o2::ft0::Geometry ft0Det{};
  // Minimal offset container compatible with UDHelpers.h expectations: getX/getY/getZ
  struct OffsetXYZ {
    double x{0}, y{0}, z{0};
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
  };
  std::array<OffsetXYZ, 1> offsetFT0{}; // iRunOffset = 0 for now
  int iRunOffset = 0;
  static constexpr uint64_t Ft0IndexA = 96;

  // Corrections
  TH3D* mEfficiency = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<CorrelationContainer> sameTpcFt0a{"sameEvent_TPC_FT0A"};
  OutputObj<CorrelationContainer> mixedTpcFt0a{"mixedEvent_TPC_FT0A"};
  OutputObj<CorrelationContainer> sameTpcFt0c{"sameEvent_TPC_FT0C"};
  OutputObj<CorrelationContainer> mixedTpcFt0c{"mixedEvent_TPC_FT0C"};
  OutputObj<CorrelationContainer> sameFt0aFt0c{"sameEvent_FT0A_FT0C"};
  OutputObj<CorrelationContainer> mixedFt0aFt0c{"mixedEvent_FT0A_FT0C"};
  HistogramRegistry registry{"registry"};

  // define global variables
  Service<ccdb::BasicCCDBManager> ccdb;
  TRandom3* gRandom = new TRandom3();
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };
  enum FITIndex {
    kFT0A = 0,
    kFT0C = 1
  };
  enum DataType {
    kReco,
    kGen
  };

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    LOGF(info, "Starting init");

    const AxisSpec axisEtaFull{90, -4., 5., "#eta"};
    std::vector<double> binEdges = {-1.5, -0.5, 0.5, 1.5, 2.5, 3.5};
    AxisSpec axisgap = {binEdges, "true gap side"};
    if (doprocessSameTpcFt0a || doprocessSameTpcFt0c || doprocessSameFt0aFt0c) {
      registry.add("truegap_inused", "truegap_inused", {HistType::kTH1D, {axisgap}});
      registry.add("trackSign", "trackSign", {HistType::kTH1D, {{10,-5,5}}});
      registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      registry.add("Phi_FT0A", "Phi_FT0A", {HistType::kTH1D, {axisPhi}});
      registry.add("Phi_FT0C", "Phi_FT0C", {HistType::kTH1D, {axisPhi}});
      registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      registry.add("Eta_FT0A", "Eta_FT0A", {HistType::kTH1D, {{20,3.,5.}}});
      registry.add("Eta_FT0C", "Eta_FT0C", {HistType::kTH1D, {{20,-4.0,-2.0}}});
      registry.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      registry.add("Nch_used_TPCFT0A", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      registry.add("Nch_used_TPCFT0C", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      registry.add("FT0Amp", "", {HistType::kTH2F, {axisChID, axisThresholdFt0a}});
      registry.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, axisThresholdFt0a}});
      registry.add("hDCAz", "DCAz after cuts; DCAz (cm); Pt", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, 0, 10}}});
      registry.add("hDCAxy", "DCAxy after cuts; DCAxy (cm); Pt", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, 0, 10}}});
      if (cfgDrawEtaPhiDis) {
        registry.add("EtaPhi", "", {HistType::kTH2F, {axisEtaFull, axisPhi}});
      }
    }
    

    if (doprocessSameTpcFt0a) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
      registry.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }
    if (doprocessSameTpcFt0c) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
      registry.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }
    if (doprocessSameFt0aFt0c) {
      registry.add("deltaEta_deltaPhi_same_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}});
      registry.add("Trig_hist_FT0A_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }

    registry.add("eventcount_TPCFT0A", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event
    registry.add("eventcount_TPCFT0C", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event
    registry.add("eventcount_FT0AFT0C", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");
    std::vector<AxisSpec> corrAxisTpcFt0a = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisPtAssoc, "p_{T} (GeV/c)"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0a, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;
    std::vector<AxisSpec> corrAxisTpcFt0c = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisPtAssoc, "p_{T} (GeV/c)"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0c, "#Delta#eta"}};
    std::vector<AxisSpec> corrAxisFt0aFt0c = {{axisSample, "Sample"},
                                              {axisVertex, "z-vtx (cm)"},
                                              {axisPtTrigger, "p_{T} (GeV/c)"},
                                              {axisPtAssoc, "p_{T} (GeV/c)"},
                                              {axisDeltaPhi, "#Delta#varphi (rad)"},
                                              {axisDeltaEtaFt0aFt0c, "#Delta#eta"}};

    if (doprocessSameTpcFt0a) {
      sameTpcFt0a.setObject(new CorrelationContainer("sameEvent_TPC_FT0A", "sameEvent_TPC_FT0A", corrAxisTpcFt0a, effAxis, userAxis));
      mixedTpcFt0a.setObject(new CorrelationContainer("mixedEvent_TPC_FT0A", "mixedEvent_TPC_FT0A", corrAxisTpcFt0a, effAxis, userAxis));
    }
    if (doprocessSameTpcFt0c) {
      sameTpcFt0c.setObject(new CorrelationContainer("sameEvent_TPC_FT0C", "sameEvent_TPC_FT0C", corrAxisTpcFt0c, effAxis, userAxis));
      mixedTpcFt0c.setObject(new CorrelationContainer("mixedEvent_TPC_FT0C", "mixedEvent_TPC_FT0C", corrAxisTpcFt0c, effAxis, userAxis));
    }
    if (doprocessSameFt0aFt0c) {
      sameFt0aFt0c.setObject(new CorrelationContainer("sameEvent_FT0A_FT0C", "sameEvent_FT0A_FT0C", corrAxisFt0aFt0c, effAxis, userAxis));
      mixedFt0aFt0c.setObject(new CorrelationContainer("mixedEvent_FT0A_FT0C", "mixedEvent_FT0A_FT0C", corrAxisFt0aFt0c, effAxis, userAxis));
    }

    ft0Det.calculateChannelCenter();

    LOGF(info, "End of init");
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMax) {
      return false;
    }
    if (!track.isPVContributor()) {
      return false;
    }
    if (cfgGlobalTrack && !(track.hasITS() && track.hasTPC())) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < cfgCutTPCCrossedRows) {
      return false;
    }
    auto tpcClu = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
    if (tpcClu < cfgCutTPCclu) {
      return false;
    }
    if (track.tpcChi2NCl() >= cfgMaxTPCChi2NCl) {
      return false;
    }
    if (track.itsNCls() < cfgCutITSclu) {
      return false;
    }
    if (cfgDcaZ && !(std::fabs(track.dcaZ()) < cfgDcaZCut)) {
      return false;
    }
    double dcaLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
    if (cfgDcaXY && !(std::fabs(track.dcaXY()) < dcaLimit)) {
      return false;
    }

    return true;
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta);
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int zBin = mEfficiency->GetZaxis()->FindBin(posZ);
      eff = mEfficiency->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    return true;
  }

  template <typename TTracks>
  double trackCounter(TTracks tracks, float posZ)
  {
    float weff1 = 1;
    double nTracksRaw = 0.;
    double nTracksCorrected = 0.;
    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      nTracksRaw += 1.;
      auto momentum = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
      double pt = track1.pt();
      double phi = RecoDecay::phi(momentum);
      double eta = RecoDecay::eta(momentum);
      if (!getEfficiencyCorrection(weff1, eta, pt, posZ))
        continue;
      nTracksCorrected += weff1;
    }
    if (cfgUseNchCorrected)
     return nTracksCorrected;
    else
     return nTracksRaw;
  }

  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    float weff1 = 1;
    float vtxz = collision.posZ();
    registry.fill(HIST("truegap_inused"), collision.truegapside());
    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), vtxz);
    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      auto momentum = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
      double pt = track1.pt();
      double phi = RecoDecay::phi(momentum);
      double eta = RecoDecay::eta(momentum);
      if (!getEfficiencyCorrection(weff1, eta, pt, vtxz))
        continue;
      registry.fill(HIST("Phi"), RecoDecay::constrainAngle(phi, 0.0));
      registry.fill(HIST("Eta"), eta);
      registry.fill(HIST("EtaCorrected"), eta, weff1);
      registry.fill(HIST("pT"), pt);
      registry.fill(HIST("pTCorrected"), pt, weff1);
      registry.fill(HIST("hDCAz"), track1.dcaZ(), pt);
      registry.fill(HIST("hDCAxy"), track1.dcaXY(), pt);
      registry.fill(HIST("trackSign"), track1.sign());
    }
    
  }

  template <CorrelationContainer::CFStep step, typename TTracks>
  void fillCorrelationsTPCFT0(TTracks tracks1, aod::UDCollisionFITBits const& fitBits, float posZ, int system, int corType) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    if (fitBits.size() > 1) {
      LOGF(fatal, "fillCorrelationsTPCFT0(): fitBits.size() = %d (expected 0 or 1)", fitBits.size());
      return;
    }
    if (fitBits.size() == 0) {
      // no TPC-FT0 correlations for this collision
      return;
    }
    
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {
      if (!trackSelected(track1))
        continue;
      auto momentum = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
      double pt = track1.pt();
      double tr_phi = RecoDecay::phi(momentum);
      double tr_eta = RecoDecay::eta(momentum);
      if (!getEfficiencyCorrection(triggerWeight, tr_eta, pt, posZ))
        continue;
      if (system == SameEvent) {
        if (corType == kFT0C) {
          registry.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, pt, triggerWeight);
        } else if (corType == kFT0A) {
          registry.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, pt, triggerWeight);
        }
        if (cfgDrawEtaPhiDis) {
          registry.fill(HIST("EtaPhi"), tr_eta, tr_phi, triggerWeight);
        }
      }

      // Only one row per collision is expected.
      auto row = fitBits.begin();
      const auto w1 = udhelpers::makeBits256(row.thr1W0(), row.thr1W1(), row.thr1W2(), row.thr1W3());
      const auto w2 = udhelpers::makeBits256(row.thr2W0(), row.thr2W1(), row.thr2W2(), row.thr2W3());
      for (int chanelid = 0; chanelid < udhelpers::kFT0Bits; ++chanelid) {
        // skip A channels if filling TPC-FT0C
        if (corType == kFT0C && chanelid < udhelpers::kFT0AChannels)
          continue;
        // skip C channels if filling TPC-FT0A
        if (corType == kFT0A && chanelid >= udhelpers::kFT0AChannels)
          continue;

        // at least 1 fired channel
        if (!udhelpers::testBit(w1, chanelid)) {
          continue;
        }
        double fitCh_phi = 0., fitCh_eta = 0.;
        const bool ok = udhelpers::getPhiEtaFromFitBit(ft0Det, chanelid, offsetFT0, iRunOffset, fitCh_phi, fitCh_eta);
        if (!ok)
          continue;
        float thr = 1.;
        if (udhelpers::testBit(w2, chanelid))
          thr = 2.;
        
        if (system == SameEvent) {
          if (corType == kFT0C) {
            registry.fill(HIST("Phi_FT0C"), RecoDecay::constrainAngle(fitCh_phi, 0.0));
            registry.fill(HIST("Eta_FT0C"), fitCh_eta);
          } else if (corType == kFT0A) {
            registry.fill(HIST("Phi_FT0A"), RecoDecay::constrainAngle(fitCh_phi, 0.0));
            registry.fill(HIST("Eta_FT0A"), fitCh_eta);
          }
          if (cfgDrawEtaPhiDis) {
            registry.fill(HIST("EtaPhi"), fitCh_eta, fitCh_phi, thr);
          }
          registry.fill(HIST("FT0Amp"), chanelid, thr); 
        }

        float deltaPhi = RecoDecay::constrainAngle(tr_phi - fitCh_phi, -PIHalf);
        float deltaEta = tr_eta - fitCh_eta;
        // fill the right sparse and histograms
        if (system == SameEvent) {
          if (corType == kFT0A) {
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0A"), deltaPhi, deltaEta, thr * triggerWeight);
            sameTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, pt, pt, deltaPhi, deltaEta, thr * triggerWeight);
          } else if (corType == kFT0C) {
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0C"), deltaPhi, deltaEta, thr * triggerWeight);
            sameTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, pt, pt, deltaPhi, deltaEta, thr * triggerWeight);
          }
        } else if (system == MixedEvent) {
          if (corType == kFT0A) {
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0A"), deltaPhi, deltaEta, thr * triggerWeight);
            mixedTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, pt, pt, deltaPhi, deltaEta, thr * triggerWeight);
          } else if (corType == kFT0C) {
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0C"), deltaPhi, deltaEta, thr * triggerWeight);
            mixedTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, pt, pt, deltaPhi, deltaEta, thr * triggerWeight);
          }
        }
      }
    }

  }

  template <CorrelationContainer::CFStep step>
  void fillCorrelationsFT0AFT0C(aod::UDCollisionFITBits const& fitBits1, aod::UDCollisionFITBits const& fitBits2, float posZ, int system) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    if (fitBits1.size() > 1 || fitBits2.size() > 1) {
      LOGF(fatal, "fillCorrelationsFT0AFT0C(): fitBits1.size() = %d, fitBits2.size() = %d (expected 0 or 1)", fitBits1.size(), fitBits2.size());
      return;
    }
    if (fitBits1.size() == 0 || fitBits2.size() == 0) {
      // no TPC-FT0 correlations for this collision
      return;
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    // Only one row per collision is expected.
    auto rowA = fitBits1.begin();
    const auto w1A = udhelpers::makeBits256(rowA.thr1W0(), rowA.thr1W1(), rowA.thr1W2(), rowA.thr1W3());
    const auto w2A = udhelpers::makeBits256(rowA.thr2W0(), rowA.thr2W1(), rowA.thr2W2(), rowA.thr2W3());
    auto rowC = fitBits2.begin();
    const auto w1C = udhelpers::makeBits256(rowC.thr1W0(), rowC.thr1W1(), rowC.thr1W2(), rowC.thr1W3());
    const auto w2C = udhelpers::makeBits256(rowC.thr2W0(), rowC.thr2W1(), rowC.thr2W2(), rowC.thr2W3());
    for (int chanelidA = 0; chanelidA < udhelpers::kFT0AChannels; ++chanelidA) {
      // at least 1 fired channel
      if (!udhelpers::testBit(w1A, chanelidA)) {
        continue;
      }
      double fitCh_phiA = 0., fitCh_etaA = 0.;
      const bool okA = udhelpers::getPhiEtaFromFitBit(ft0Det, chanelidA, offsetFT0, iRunOffset, fitCh_phiA, fitCh_etaA);
      if (!okA)
        continue;
      float thrA = 1.;
      if (udhelpers::testBit(w2A, chanelidA))
        thrA = 2.;
      
      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist_FT0A_FT0C"), fSampleIndex, posZ, 0.5, thrA);
      }
      for (int chanelidC = udhelpers::kFT0AChannels; chanelidC < udhelpers::kFT0Bits; ++chanelidC) {
        // at least 1 fired channel
        if (!udhelpers::testBit(w1C, chanelidC)) {
          continue;
        }
        double fitCh_phiC = 0., fitCh_etaC = 0.;
        const bool okC = udhelpers::getPhiEtaFromFitBit(ft0Det, chanelidC, offsetFT0, iRunOffset, fitCh_phiC, fitCh_etaC);
        if (!okC)
          continue;
        float thrC = 1.;
        if (udhelpers::testBit(w2A, chanelidC))
          thrC = 2.;

        float deltaPhi = RecoDecay::constrainAngle(fitCh_phiA - fitCh_phiC, -PIHalf);
        float deltaEta = fitCh_etaA - fitCh_etaC;

        // fill the right sparse and histograms
        if (system == SameEvent) {
          registry.fill(HIST("deltaEta_deltaPhi_same_FT0A_FT0C"), deltaPhi, deltaEta, thrA * thrC * triggerWeight);
          sameFt0aFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, 0.5, 0.5, deltaPhi, deltaEta, thrA * thrC * triggerWeight);
        } else if (system == MixedEvent) {
          registry.fill(HIST("deltaEta_deltaPhi_mixed_FT0A_FT0C"), deltaPhi, deltaEta, thrA * thrC * triggerWeight);
          mixedFt0aFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, 0.5, 0.5, deltaPhi, deltaEta, thrA * thrC * triggerWeight);
        }
      }
    }

  }

  void processSameTpcFt0a(UDCollisionsFull::iterator const& collision, UDTracksFull const& tracks, aod::UDCollisionFITBits const& fitBits)
  {
    LOGF(info, "processSameTpcFt0a: collisionId:%d, tracks.size() = %d", collision.globalIndex(), tracks.size());
    LOGF(info, "processSameTpcFt0a: collisionId:%d, fitBits.size() = %d", collision.globalIndex(), fitBits.size());
    if (fitBits.size() > 1) {
      LOGF(fatal, "processSameTpcFt0a: fitBits.size() = %d (expected 0 or 1)", fitBits.size());
      return;
    }

    auto currentRunNumber = collision.runNumber();
    auto runDuration = ccdb->getRunDuration(currentRunNumber);
    loadCorrection(runDuration.first);

    double Nch = trackCounter(tracks, collision.posZ());
    if (Nch < cfgMinMult || Nch >= cfgMaxMult) {
      return;
    }
    registry.fill(HIST("eventcount_TPCFT0A"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);
    registry.fill(HIST("Nch_used_TPCFT0A"), Nch);

    sameTpcFt0a->fillEvent(Nch, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks, fitBits, collision.posZ(), SameEvent, kFT0A);
  }
  PROCESS_SWITCH(FlowLongRangeCorrUpc, processSameTpcFt0a, "Process same event for TPC-FT0 correlation", true);

  void processMixedTpcFt0a(UDCollisionsFull const& collisions, UDTracksFull const& tracks, aod::UDCollisionFITBits const& fitBits)
  {
    auto getTracksSize = [&tracks, this](UDCollisionsFull::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::udtrack::udCollisionId, collision.globalIndex(), this->cache);
      auto mult = trackCounter(associatedTracks, collision.posZ());
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<UDCollisionsFull, UDTracksFull, UDTracksFull, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      double Nch1 = trackCounter(tracks1, collision1.posZ());
      double Nch2 = trackCounter(tracks2, collision2.posZ());
      if (Nch1 < cfgMinMult || Nch1 >= cfgMaxMult)
        continue;
      if (Nch2 < cfgMinMult || Nch2 >= cfgMaxMult)
        continue;

      auto fitBits2 = fitBits.sliceByCached(o2::aod::udtrack::udCollisionId, collision2.globalIndex(), this->cache);
      if (fitBits2.size() > 1) {
        LOGF(fatal, "processMixedTpcFt0a: collision2 Id:%d, fitBits2.size() = %d (expected 0 or 1)", collision2.globalIndex(), fitBits2.size());
        return;
      }
      if (fitBits2.size() == 0) {
        continue;
      }

      registry.fill(HIST("eventcount_TPCFT0A"), MixedEvent); // fill the mixed event in the 3 bin
      auto currentRunNumber = collision1.runNumber();
      auto runDuration = ccdb->getRunDuration(currentRunNumber);
      loadCorrection(runDuration.first);

      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks1, fitBits2, collision1.posZ(), MixedEvent, kFT0A);
    }
  }
  PROCESS_SWITCH(FlowLongRangeCorrUpc, processMixedTpcFt0a, "Process mixed events for TPC-FT0A correlation", true);

  void processSameTpcFt0c(UDCollisionsFull::iterator const& collision, UDTracksFull const& tracks, aod::UDCollisionFITBits const& fitBits)
  {
    if (fitBits.size() > 1) {
      LOGF(fatal, "processSameTpcFt0c: fitBits.size() = %d (expected 0 or 1)", fitBits.size());
      return;
    }

    auto currentRunNumber = collision.runNumber();
    auto runDuration = ccdb->getRunDuration(currentRunNumber);
    loadCorrection(runDuration.first);

    double Nch = trackCounter(tracks, collision.posZ());
    if (Nch < cfgMinMult || Nch >= cfgMaxMult) {
      return;
    }
    registry.fill(HIST("eventcount_TPCFT0C"), SameEvent); // because its same event i put it in the 1 bin
    if (!doprocessSameTpcFt0a)
      fillYield(collision, tracks);
    registry.fill(HIST("Nch_used_TPCFT0C"), Nch);

    sameTpcFt0c->fillEvent(Nch, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks, fitBits, collision.posZ(), SameEvent, kFT0C);
  }
  PROCESS_SWITCH(FlowLongRangeCorrUpc, processSameTpcFt0c, "Process same events for TPC-FT0C correlation", false);

  void processMixedTpcFt0c(UDCollisionsFull const& collisions, UDTracksFull const& tracks, aod::UDCollisionFITBits const& fitBits)
  {
    auto getTracksSize = [&tracks, this](UDCollisionsFull::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::udtrack::udCollisionId, collision.globalIndex(), this->cache);
      auto mult = trackCounter(associatedTracks, collision.posZ());
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<UDCollisionsFull, UDTracksFull, UDTracksFull, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      double Nch1 = trackCounter(tracks1, collision1.posZ());
      double Nch2 = trackCounter(tracks2, collision2.posZ());
      if (Nch1 < cfgMinMult || Nch1 >= cfgMaxMult)
        continue;
      if (Nch2 < cfgMinMult || Nch2 >= cfgMaxMult)
        continue;

      auto fitBits2 = fitBits.sliceByCached(o2::aod::udtrack::udCollisionId, collision2.globalIndex(), this->cache);
      if (fitBits2.size() > 1) {
        LOGF(fatal, "processMixedTpcFt0a: collision2 Id:%d, fitBits2.size() = %d (expected 0 or 1)", collision2.globalIndex(), fitBits2.size());
        return;
      }
      if (fitBits2.size() == 0) {
        continue;
      }

      registry.fill(HIST("eventcount_TPCFT0C"), MixedEvent); // fill the mixed event in the 3 bin
      auto currentRunNumber = collision1.runNumber();
      auto runDuration = ccdb->getRunDuration(currentRunNumber);
      loadCorrection(runDuration.first);

      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks1, fitBits2, collision1.posZ(), MixedEvent, kFT0C);
    }
  }
  PROCESS_SWITCH(FlowLongRangeCorrUpc, processMixedTpcFt0c, "Process mixed events for TPC-FT0C correlation", false);

  void processSameFt0aFt0c(UDCollisionsFull::iterator const& collision, UDTracksFull const& tracks, aod::UDCollisionFITBits const& fitBits)
  {
    if (fitBits.size() > 1) {
      LOGF(fatal, "processSameTpcFt0c: fitBits.size() = %d (expected 0 or 1)", fitBits.size());
      return;
    }

    auto currentRunNumber = collision.runNumber();
    auto runDuration = ccdb->getRunDuration(currentRunNumber);
    loadCorrection(runDuration.first);

    double Nch = trackCounter(tracks, collision.posZ());
    if (Nch < cfgMinMult || Nch >= cfgMaxMult) {
      return;
    }
    registry.fill(HIST("eventcount_FT0AFT0C"), SameEvent); // because its same event i put it in the 1 bin
    if (!doprocessSameTpcFt0a)
      fillYield(collision, tracks);

    sameFt0aFt0c->fillEvent(Nch, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(fitBits, fitBits, collision.posZ(), SameEvent);
  }
  PROCESS_SWITCH(FlowLongRangeCorrUpc, processSameFt0aFt0c, "Process same events for FT0A-FT0C correlation", false);

  void processMixedFt0aFt0c(UDCollisionsFull const& collisions, UDTracksFull const& tracks, aod::UDCollisionFITBits const& fitBits)
  {
    auto getTracksSize = [&tracks, this](UDCollisionsFull::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::udtrack::udCollisionId, collision.globalIndex(), this->cache);
      auto mult = trackCounter(associatedTracks, collision.posZ());
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<UDCollisionsFull, UDTracksFull, UDTracksFull, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      double Nch1 = trackCounter(tracks1, collision1.posZ());
      double Nch2 = trackCounter(tracks2, collision2.posZ());
      if (Nch1 < cfgMinMult || Nch1 >= cfgMaxMult)
        continue;
      if (Nch2 < cfgMinMult || Nch2 >= cfgMaxMult)
        continue;

      auto fitBits1 = fitBits.sliceByCached(o2::aod::udtrack::udCollisionId, collision1.globalIndex(), this->cache);
      auto fitBits2 = fitBits.sliceByCached(o2::aod::udtrack::udCollisionId, collision2.globalIndex(), this->cache);
      if (fitBits1.size() > 1 || fitBits2.size() > 1) {
        LOGF(fatal, "processMixedTpcFt0a: fitBits1.size() = %d, fitBits2.size() = %d (expected 0 or 1)", fitBits1.size(), fitBits2.size());
        return;
      }
      if (fitBits1.size() == 0 || fitBits2.size() == 0) {
        continue;
      }

      registry.fill(HIST("eventcount_FT0AFT0C"), MixedEvent); // fill the mixed event in the 3 bin
      auto currentRunNumber = collision1.runNumber();
      auto runDuration = ccdb->getRunDuration(currentRunNumber);
      loadCorrection(runDuration.first);

      fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(fitBits1, fitBits2, collision1.posZ(), MixedEvent);
    }
  }
  PROCESS_SWITCH(FlowLongRangeCorrUpc, processMixedFt0aFt0c, "Process mixed events for FT0A-FT0C correlation", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowLongRangeCorrUpc>(cfgc),
  };
}
