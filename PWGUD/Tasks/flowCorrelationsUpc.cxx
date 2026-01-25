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
/// \author Mingrui Zhao (mingrui.zhao@cern.ch, mingrui.zhao@mail.labz0.org)
/// copied from Thor Jensen (thor.kjaersgaard.jensen@cern.ch) and Debojit Sarkar (debojit.sarkar@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TRandom3.h"

#include <vector>

namespace o2::aod
{
namespace flowcorrupc
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, int);
}
DECLARE_SOA_TABLE(Multiplicity, "AOD", "MULTIPLICITY",
                  flowcorrupc::Multiplicity);

} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct CalcNchUpc {
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")

  // Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  using UdTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using UdTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDCollisionSelExtras>;

  Produces<aod::Multiplicity> multiplicityNch;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    AxisSpec axisNch = {100, 0, 100};
    AxisSpec axisVrtx = {10, -10, 10};

    registry.add("Ncharge", "N_{charge}", {HistType::kTH1D, {axisNch}});
    registry.add("zVtx_all", "zVtx_all", {HistType::kTH1D, {axisVrtx}});
  }

  void process(UDCollisionsFull::iterator const& collision, UdTracksFull const& tracks)
  {
    multiplicityNch(tracks.size());
    registry.fill(HIST("Ncharge"), tracks.size());
    registry.fill(HIST("zVtx_all"), collision.posZ());
  }
};

struct FlowCorrelationsUpc {
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMinMult, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgMaxMult, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, true, "enable trigger pT < associated pT cut")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, true, "enable trigger pT < associated pT cut in mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")

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
  Configurable<float> cfgGapSideSelection{"cfgGapSideSelection", 2, "gap selection"};

  // make the filters and cuts.
  // Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZVtxCut) && (aod::flowcorrupc::multiplicity) > cfgMinMult && (aod::flowcorrupc::multiplicity) < cfgMaxMult && (aod::evsel::sel8) == true;
  // Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  using UdTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using UdTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::Multiplicity, aod::UDCollisionSelExtras>;

  // Define the outputs
  OutputObj<CorrelationContainer> same{Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixed{Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    LOGF(info, "Starting init");
    // Make histograms to check the distributions after cuts
    registry.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
    registry.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
    registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
    registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
    registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});

    registry.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});

    registry.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVertex, "z-vtx (cm)"},
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
    if (dPhiStar < -constants::math::PI)
      dPhiStar = -constants::math::TwoPI - dPhiStar;

    return dPhiStar;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    // UPC selection
    if (!track.isPVContributor()) {
      return false;
    }
    constexpr float kDcazCut = 2.0;
    if (!(std::fabs(track.dcaZ()) < kDcazCut)) {
      return false;
    }
    double dcaLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
    if (!(std::fabs(track.dcaXY()) < dcaLimit)) {
      return false;
    }
    constexpr int kMinITSClusters = 5;
    constexpr int kMaxTPCChi2NCl = 4;

    if (track.itsClusterSizes() <= kMinITSClusters) {
      return false;
    }
    if (track.tpcChi2NCl() >= kMaxTPCChi2NCl) {
      return false;
    }
    if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMax)
      return false;
    return true;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    for (auto const& track1 : tracks) {
      auto momentum1 = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
      registry.fill(HIST("Phi"), RecoDecay::phi(momentum1));
      registry.fill(HIST("Eta"), RecoDecay::eta(momentum1));
      registry.fill(HIST("pT"), track1.pt());
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks>
  void fillCorrelations(TTracks tracks1, TTracks tracks2, float posZ, int system, int runnum) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt());
      }

      for (auto const& track2 : tracks2) {
        if (!trackSelected(track2))
          continue;

        if (track1.globalIndex() == track2.globalIndex())
          continue; // For pt-differential correlations, skip if the trigger and associate are the same track
        if (system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt

        auto momentum1 = std::array<double, 3>{track1.px(), track1.py(), track1.pz()};
        auto momentum2 = std::array<double, 3>{track2.px(), track2.py(), track2.pz()};
        double pt2 = RecoDecay::pt(momentum2);
        double phi1 = RecoDecay::phi(momentum1);
        double phi2 = RecoDecay::phi(momentum2);
        float deltaPhi = RecoDecay::constrainAngle(phi1 - phi2, -PIHalf);
        float deltaEta = RecoDecay::eta(momentum1) - RecoDecay::eta(momentum2);

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

        // fill the right sparse and histograms
        if (system == SameEvent) {
          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta);
        } else if (system == MixedEvent) {
          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta);
        }
      }
    }
  }

  void processSame(UDCollisionsFull::iterator const& collision, UdTracksFull const& tracks)
  {
    if (tracks.size() < cfgMinMult || tracks.size() > cfgMaxMult) {
      return;
    }
    if (collision.trs() == 0) {
      return;
    }

    int gapSide = collision.gapSide();
    const int minGapSide = 0;
    const int maxGapSide = 2;
    if (gapSide > minGapSide && gapSide < maxGapSide) {
      return;
    }

    int trueGapSide = sgSelector.trueGap(collision, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC);
    gapSide = trueGapSide;
    if (gapSide == cfgGapSideSelection) {
      return;
    }

    int runIndex = collision.runNumber();

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, runIndex); // fill the SE histogram and Sparse
  }
  PROCESS_SWITCH(FlowCorrelationsUpc, processSame, "Process same event", true);

  // event mixing

  SliceCache cache;
  using MixedBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::flowcorrupc::Multiplicity>;

  // the process for filling the mixed events
  void processMixed(UDCollisionsFull const& collisions, UdTracksFull const& tracks)
  {
    MixedBinning binningOnVtxAndMult{{vtxMix, multMix}, true}; // true is for 'ignore overflows' (true by default)
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<UDCollisionsFull, UdTracksFull, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto const& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (tracks1.size() < cfgMinMult || tracks1.size() > cfgMaxMult || tracks2.size() < cfgMinMult || tracks2.size() > cfgMaxMult) {
        continue;
      }
      if (collision1.trs() == 0 || collision2.trs() == 0) {
        continue;
      }

      const int minGapSide = 0;
      const int maxGapSide = 2;
      if (collision1.gapSide() > minGapSide && collision1.gapSide() < maxGapSide) {
        continue;
      }
      if (collision2.gapSide() > minGapSide && collision2.gapSide() < maxGapSide) {
        continue;
      }

      int trueGapSide = sgSelector.trueGap(collision1, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC);
      int gapSide = trueGapSide;
      if (gapSide == cfgGapSideSelection) {
        continue;
      }
      trueGapSide = sgSelector.trueGap(collision2, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC);
      gapSide = trueGapSide;
      if (gapSide == cfgGapSideSelection) {
        continue;
      }
      registry.fill(HIST("eventcount"), MixedEvent);                                                                                         // fill the mixed event in the 3 bin
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, collision1.runNumber()); // fill the ME histogram and Sparse
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
