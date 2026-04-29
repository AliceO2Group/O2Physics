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

/// \file flowCorrelationsUpc.cxx
/// \brief Provides a sparse with usefull two particle correlation info
/// \author Yongxi Du (yongxi.du@cern.ch), Mingrui Zhao (mingrui.zhao@cern.ch, mingrui.zhao@mail.labz0.org)
/// copied from Thor Jensen (thor.kjaersgaard.jensen@cern.ch) and Debojit Sarkar (debojit.sarkar@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
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
#include <Framework/StringHelpers.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/PxPyPzE4D.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace o2::aod
{
namespace flowcorrupc
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, int);
DECLARE_SOA_COLUMN(Truegapside, truegapside, int);
} // namespace flowcorrupc
DECLARE_SOA_TABLE(Multiplicity, "AOD", "MULTIPLICITY",
                  flowcorrupc::Multiplicity);
DECLARE_SOA_TABLE(Truegapside, "AOD", "TRUEGAPSIDE", flowcorrupc::Truegapside);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct CalcNchUpc {
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.1f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.9f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")

  // Added UPC Cuts
  SGSelector sgSelector;
  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};

  // Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  using UdTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using UdTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDCollisionSelExtras>;

  Produces<aod::Multiplicity> multiplicityNch;
  Produces<aod::Truegapside> truegapside;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    AxisSpec axisNch = {100, 0, 100};
    AxisSpec axisVrtx = {10, -10, 10};
    // AxisSpec axisgap = {12, -6, 6};
    // std::vector<AxisSpec> trueGapBins = {-2, -1, 0, 1, 2, 3};
    // AxisSpec axisgap = {trueGapBins, "true gap side"};

    std::vector<double> binEdges = {-1.5, -0.5, 0.5, 1.5, 2.5, 3.5};
    AxisSpec axisgap = {binEdges, "true gap side"};
    registry.add("truegap", "truegap", {HistType::kTH1D, {axisgap}});

    registry.add("Ncharge", "N_{charge}", {HistType::kTH1D, {axisNch}});
    registry.add("zVtx_all", "zVtx_all", {HistType::kTH1D, {axisVrtx}});
    registry.add("Nch_vs_zVtx", "Nch vs zVtx", {HistType::kTH2D, {axisVrtx, axisNch}});
    // registry.add("truegap", "truegap", {HistType::kTH1D, {axisgap}});
  }

  void process(UDCollisionsFull::iterator const& collision, UdTracksFull const& tracks)
  {
    multiplicityNch(tracks.size());
    truegapside(sgSelector.trueGap(collision, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC));
    // LOG(info) << "truegapside=" <<  sgSelector.trueGap(collision, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC);
    registry.fill(HIST("Ncharge"), tracks.size());
    registry.fill(HIST("zVtx_all"), collision.posZ());
  }
};

struct FlowCorrelationsUpc {
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgIfVertex, bool, false, "choose vertex or not")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.1f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.9f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMinMult, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgMaxMult, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, true, "enable trigger pT < associated pT cut")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, true, "enable trigger pT < associated pT cut in mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgDcaxy, bool, true, "choose dcaxy")
  O2_DEFINE_CONFIGURABLE(cfgDcaz, bool, false, "choose dcaz")
  O2_DEFINE_CONFIGURABLE(cfgDcazCut, float, 10.0, "dcaz cut")
  O2_DEFINE_CONFIGURABLE(cfgMaxTPCChi2NCl, int, 4, "tpcchi2")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancy, bool, true, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 1000, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum number of crossed TPC Rows")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum number of found TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum number of ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgGlobalTrack, bool, true, "require TPC+ITS track")
  O2_DEFINE_CONFIGURABLE(cfgUseNchCorrected, bool, true, "use corrected Nch for X axis")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis vtxMix{"vtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis multMix{"multMix", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for mixed event histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};

  // Added UPC Cuts
  SGSelector sgSelector;
  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};
  ConfigurableAxis axisIndependent{"axisIndependent", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60}, "X axis for histograms"};
  ConfigurableAxis axisNch{"axisNch", {300, 0, 300}, "N_{ch}"};

  // Corrections
  TH3D* mEfficiency = nullptr;
  bool correctionsLoaded = false;

  // make the filters and cuts.
  Filter trackFilter = (aod::udtrack::isPVContributor == true);
  Filter collisionFilter = ((aod::udcollision::gapSide == (uint8_t)1 || aod::udcollision::gapSide == (uint8_t)0) && (cfgIfVertex == false || aod::collision::posZ < cfgZVtxCut) && (!cfgCutOccupancy || (aod::udcollision::occupancyInTime > 0 && aod::udcollision::occupancyInTime < cfgCutOccupancyHigh)) && (aod::flowcorrupc::truegapside == 1 || aod::flowcorrupc::truegapside == 0));

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};

  TAxis* fPtAxis = nullptr;
  int lastRunNumber = -1;
  std::vector<int> runNumbers; // map of TH3 histograms for all runs
  std::vector<float> efficiencyCache;

  using UdTracks = soa::Filtered<soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>>;
  using UdTracksFull = soa::Filtered<soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>>;

  using UDCollisionsFull = soa::Filtered<soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::Multiplicity, aod::Truegapside, aod::UDCollisionSelExtras>>;

  // Define the outputs
  OutputObj<CorrelationContainer> same{Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixed{Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    LOGF(info, "Starting init");
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    // Make histograms to check the distributions after cuts
    registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
    registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
    registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
    registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
    registry.add("EtaCorrected", "Eta corrected", {HistType::kTH1D, {axisEta}});
    registry.add("pTCorrected", "pT corrected", {HistType::kTH1D, {axisPtTrigger}});

    registry.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisIndependent, axisPtTrigger}}});

    registry.add("eventcont", "bin", {HistType::kTH1F, {{10, 0, 10, "bin"}}});                                     // histogram to see how many events are in the same and mixed event
    registry.add("deltaEta_deltaPhi_same", "deltaeta-deltaphi", {HistType::kTH2D, {axisDeltaEta, axisDeltaPhi}});  // histogram to check the delta eta and delta phi distribution
    registry.add("deltaEta_deltaPhi_mixed", "deltaeta-deltaphi", {HistType::kTH2D, {axisDeltaEta, axisDeltaPhi}}); // histogram to check the delta eta and delta phi distribution
    registry.add("Nch_raw_vs_independent", "Raw vs Independent", {HistType::kTH2D, {axisMultiplicity, axisIndependent}});

    // if (doprocessSim) {
    //   registry.add("eventCounterMC", "Number of MC Events;; Count", {HistType::kTH1D, {{5, 0, 5}}});
    //   registry.add("hVtxZMC", "Vexter Z distribution (MC)", {HistType::kTH1D, {axisVertex}});
    //   registry.add("hMultMC", "Multiplicity distribution (MC)", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    //   registry.add("numberOfTracksMC", "Number of MC tracks;; Count", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    // }

    o2::framework::AxisSpec axis = axisPtTrigger;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisIndependent, "Independent (N_{ch} corrected)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisVertexEfficiency, "z-vtx (cm)"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisEtaEfficiency, "#eta"},
    };
    std::vector<AxisSpec> userAxis;

    same.setObject(new CorrelationContainer(Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer(Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
  }
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  template <typename TTrack>
  float getDPhiStar(TTrack const& track1, TTrack const& track2, float radius, int runnum, float phi1, float phi2)
  {
    float charge1 = track1.sign();
    float charge2 = track2.sign();

    float pt1 = track1.pt();
    float pt2 = track2.pt();

    int fbSign = 1;

    int zzo = 544868;
    if (runnum >= zzo) {
      fbSign = -1;
    }

    float dPhiStar = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);

    if (dPhiStar > constants::math::PI)
      dPhiStar = constants::math::TwoPI - dPhiStar;
    return dPhiStar;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    // registry.fill(HIST("hTrackCount"), 0.5);
    // UPC selection
    if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMax) {
      return false;
    }
    if (cfgGlobalTrack && !(track.hasITS() && track.hasTPC())) {
      return false;
    }
    // registry.fill(HIST("hTrackCount"), 1.5);
    if (cfgDcaz && !(std::fabs(track.dcaZ()) < cfgDcazCut)) {
      return false;
    }
    // registry.fill(HIST("hTrackCount"), 2.5);
    double dcaLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
    if (cfgDcaxy && !(std::fabs(track.dcaXY()) < dcaLimit)) {
      return false;
    }
    // registry.fill(HIST("hTrackCount"), 3.5);
    if (track.itsNCls() <= cfgCutITSclu) {
      return false;
    }
    // registry.fill(HIST("hTrackCount"), 4.5);
    if (track.tpcChi2NCl() >= cfgMaxTPCChi2NCl) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < cfgCutTPCCrossedRows) {
      return false;
    }
    auto tpcClu = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
    if (tpcClu < cfgCutTPCclu) {
      return false;
    }
    // registry.fill(HIST("hTrackCount"), 5.5);
    return true;
  }

  void loadCorrections(uint64_t timestamp)
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
    if (eff <= 0)
      return false;
    weight_nue = 1. / eff;
    return true;
  }

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua)
  {
    float eff = 1.;
    if (mEfficiency)
      eff = mEfficiency->GetBinContent(mEfficiency->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    weight_nua = 1.; // Set to 1 as NUA weight is not being used
    return true;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks, float vtxz) // function to fill the yield and etaphi histograms.
  {
    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      auto momentum = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
      double pt = RecoDecay::pt(momentum);
      double phi = RecoDecay::phi(momentum);
      double eta = RecoDecay::eta(momentum);
      float weff = 1.;
      if (!getEfficiencyCorrection(weff, eta, pt, vtxz))
        continue;

      registry.fill(HIST("Phi"), phi);
      registry.fill(HIST("Eta"), eta);
      registry.fill(HIST("pT"), pt);
      registry.fill(HIST("EtaCorrected"), eta, weff);
      registry.fill(HIST("pTCorrected"), pt, weff);
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks>
  void fillCorrelations(TTracks tracks1, TTracks tracks2, float posZ, int system, int runnum, float vtxz, float eventWeight, double independent) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    if (mEfficiency) {
      efficiencyCache.clear();
      efficiencyCache.reserve(static_cast<int>(tracks2.size()));
      for (const auto& track2 : tracks2) {
        auto momentum = std::array<double, 3>{track2.px(), track2.py(), track2.pz()};
        double pt = RecoDecay::pt(momentum);
        double eta = RecoDecay::eta(momentum);
        float weff = 1.;
        getEfficiencyCorrection(weff, eta, pt, vtxz);
        efficiencyCache.push_back(weff);
      }
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    // loop over all tracks
    for (auto const& track1 : tracks1) {
      if (!trackSelected(track1))
        continue;

      auto momentum = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
      double pt1 = RecoDecay::pt(momentum);
      double phi1 = RecoDecay::phi(momentum);
      double eta1 = RecoDecay::eta(momentum);

      // 计算track1的权重
      float weff1 = 1., wacc1 = 1.;
      if (!setCurrentParticleWeights(weff1, wacc1, phi1, eta1, pt1, vtxz))
        continue;

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, independent, pt1, eventWeight * weff1 * wacc1);
      }

      for (auto const& track2 : tracks2) {
        if (!trackSelected(track2))
          continue;

        if (track1.globalIndex() == track2.globalIndex())
          continue;
        if (system == SameEvent && cfgUsePtOrder && pt1 <= track2.pt())
          continue;
        if (system == MixedEvent && cfgUsePtOrderInMixEvent && pt1 <= track2.pt())
          continue;

        auto momentum = std::array<double, 3>{track2.px(), track2.py(), track2.pz()};
        double pt2 = RecoDecay::pt(momentum);
        double phi2 = RecoDecay::phi(momentum);
        double eta2 = RecoDecay::eta(momentum);

        float weff2 = 1., wacc2 = 1.;
        if (mEfficiency) {
          weff2 = efficiencyCache[track2.filteredIndex()];
        } else {
          getEfficiencyCorrection(weff2, eta2, pt2, vtxz);
        }

        float deltaPhi = RecoDecay::constrainAngle(phi1 - phi2, -PIHalf);
        float deltaEta = eta1 - eta2;

        float weight = eventWeight * weff1 * weff2 * wacc1 * wacc2;

        // Merging cut
        if (std::abs(deltaEta) < cfgCutMerging) {
          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, runnum, phi1, phi2);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, runnum, phi1, phi2);
          const double kLimit = 3.0 * cfgCutMerging;
          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1, track2, rad, runnum, phi1, phi2);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms with weights
        if (system == SameEvent) {
          same->getPairHist()->Fill(step, fSampleIndex, posZ, independent, pt1, pt2, deltaPhi, deltaEta, weight);
          registry.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, weight);
        } else if (system == MixedEvent) {
          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, independent, pt1, pt2, deltaPhi, deltaEta, weight);
          registry.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, weight);
        }
      }
    }
  }

  void processSame(UDCollisionsFull::iterator const& collision, UdTracksFull const& tracks)
  {
    // LOG(info) << "Event passed filter: truegapside=" << collision.truegapside();
    if (tracks.size() < cfgMinMult || tracks.size() > cfgMaxMult) {
      return;
    }

    float vtxz = collision.posZ();
    auto currentRunNumber = collision.runNumber();
    auto runDuration = ccdb->getRunDuration(currentRunNumber);

    loadCorrections(runDuration.first);

    registry.fill(HIST("eventcont"), 3.5);

    //-----------independent---------------
    double nTracksRaw = 0.;
    double nTracksCorrected = 0.;

    for (const auto& track : tracks) {
      if (!trackSelected(track))
        continue;

      auto momentum = std::array<double, 3>{track.px(), track.py(), track.pz()};
      double pt = RecoDecay::pt(momentum);
      double eta = RecoDecay::eta(momentum);

      nTracksRaw += 1.;

      if (cfgUseNchCorrected) {
        float weff = 1.;
        if (getEfficiencyCorrection(weff, eta, pt, vtxz)) {
          nTracksCorrected += weff;
        }
      }
    }
    registry.fill(HIST("Nch_raw_vs_independent"), nTracksRaw, nTracksCorrected);

    double independent = nTracksRaw;
    if (cfgUseNchCorrected) {
      independent = nTracksCorrected;
    }

    fillYield(collision, tracks, currentRunNumber, vtxz);

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(
      tracks, tracks, collision.posZ(), SameEvent,
      currentRunNumber, vtxz, 1.0f, independent);
  }
  PROCESS_SWITCH(FlowCorrelationsUpc, processSame, "Process same event", true);

  // event mixing

  SliceCache cache;
  // using MixedBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::flowcorrupc::Multiplicity>;

  // the process for filling the mixed events
  void processMixed(UDCollisionsFull const& collisions, UdTracksFull const& tracks)
  {
    auto getTracksSize = [&tracks, this](UDCollisionsFull::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::udtrack::udCollisionId, collision.udCollisionId(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {vtxMix, multMix}, true};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<UDCollisionsFull, UdTracksFull, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache};

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (tracks1.size() < cfgMinMult || tracks1.size() > cfgMaxMult ||
          tracks2.size() < cfgMinMult || tracks2.size() > cfgMaxMult) {
        continue;
      }

      auto runDuration1 = ccdb->getRunDuration(collision1.runNumber());
      loadCorrections(runDuration1.first);

      registry.fill(HIST("eventcont"), 4.5);

      double nTracksRaw = 0.;
      double nTracksCorrected = 0.;

      for (const auto& track : tracks1) {
        if (!trackSelected(track))
          continue;

        auto momentum = std::array<double, 3>{track.px(), track.py(), track.pz()};
        double pt = RecoDecay::pt(momentum);
        double eta = RecoDecay::eta(momentum);

        nTracksRaw += 1.;

        if (cfgUseNchCorrected) {
          float weff = 1.;
          if (getEfficiencyCorrection(weff, eta, pt, collision1.posZ())) {
            nTracksCorrected += weff;
          }
        }
      }

      double independent = nTracksRaw;
      if (cfgUseNchCorrected) {
        independent = nTracksCorrected;
      }

      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(
        tracks1, tracks2, collision1.posZ(), MixedEvent,
        collision1.runNumber(), collision1.posZ(), eventWeight, independent);
    }
  }
  PROCESS_SWITCH(FlowCorrelationsUpc, processMixed, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CalcNchUpc>(cfgc),
    adaptAnalysisTask<FlowCorrelationsUpc>(cfgc),
  };
}
