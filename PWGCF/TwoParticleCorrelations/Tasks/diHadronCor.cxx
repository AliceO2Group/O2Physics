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

/// \file diHadronCor.cxx
/// \brief di-hadron correlation for O-O, Pb-Pb collisions
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  May/03/2025

#include <CCDB/BasicCCDBManager.h>
#include "TRandom3.h"
#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StepTHn.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"
#include "Common/Core/RecoDecay.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
namespace o2::aod
{
namespace di_hadron_cor
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, int);
}
DECLARE_SOA_TABLE(Multiplicity, "AOD", "MULTIPLICITY",
                  di_hadron_cor::Multiplicity);

} // namespace o2::aod

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct DiHadronCor {
  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgCutVtxZ, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.0, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMin, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMax, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMin, float, 60.0f, "Minimum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMax, float, 80.0f, "Maximum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgMixEventNumMin, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentFT0CMin, float, 0.0f, "Minimum centrality (FT0C) to cut events in filter")
  O2_DEFINE_CONFIGURABLE(cfgCentFT0CMax, float, 100.0f, "Maximum centrality (FT0C) to cut events in filter")

  SliceCache cache;
  SliceCache cacheNch;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVtxZ) && (aod::evsel::sel8) == true && (aod::cent::centFT0C > cfgCentFT0CMin) && (aod::cent::centFT0C < cfgCentFT0CMax);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  // Define the outputs
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};
  HistogramRegistry registry{"registry"};

  // define global variables
  TRandom3* gRandom = new TRandom3();
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>; // aod::CentFT0Cs
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  void init(InitContext&)
  {
    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    // Make histograms to check the distributions after cuts
    registry.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
    registry.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
    registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
    registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
    registry.add("Nch_used", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    registry.add("Centrality", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}});
    registry.add("Centrality_used", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}}); // histogram to see how many events are in the same and mixed event
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

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));
  }

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        cent = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        cent = collision.centFT0M();
        break;
      case kCentFV0A:
        cent = collision.centFV0A();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    float cent = getCentrality(collision);
    registry.fill(HIST("Centrality"), cent);

    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    for (auto const& track1 : tracks) {
      registry.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("pT"), track1.pt());
    }
  }

  template <typename TTrack, typename TTrackAssoc>
  float getDPhiStar(TTrack const& track1, TTrackAssoc const& track2, float radius, int magField)
  {
    float charge1 = track1.sign();
    float charge2 = track2.sign();

    float phi1 = track1.phi();
    float phi2 = track2.phi();

    float pt1 = track1.pt();
    float pt2 = track2.pt();

    int fbSign = (magField > 0) ? 1 : -1;

    float dPhiStar = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);

    if (dPhiStar > constants::math::PI)
      dPhiStar = constants::math::TwoPI - dPhiStar;
    if (dPhiStar < -constants::math::PI)
      dPhiStar = -constants::math::TwoPI - dPhiStar;

    return dPhiStar;
  }

  //
  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField, float cent) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    if (system == SameEvent) {
      registry.fill(HIST("Centrality_used"), cent);
      registry.fill(HIST("Nch_used"), tracks1.size());
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt());
      }

      for (auto const& track2 : tracks2) {

        if (track1.pt() <= track2.pt())
          continue; // skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (std::abs(deltaEta) < cfgCutMerging) {

          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1, track2, rad, magneticField);
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

  void processSame(AodCollisions::iterator const& collision, AodTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    fillYield(collision, tracks);
    float cent = getCentrality(collision);

    if (cfgSelCollByNch && (tracks.size() < cfgCutMultMin || tracks.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()), cent);
  }
  PROCESS_SWITCH(DiHadronCor, processSame, "Process same event", true);

  // the process for filling the mixed events
  void processMixed(AodCollisions const& collisions, AodTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](AodCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<AodCollisions, AodTracks, AodTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto const& [collision1, tracks1, collision2, tracks2] : pair) {
      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = getCentrality(collision1);
      float cent2 = getCentrality(collision2);

      if (!cfgSelCollByNch && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()), cent1);
    }
  }

  PROCESS_SWITCH(DiHadronCor, processMixed, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiHadronCor>(cfgc),
  };
}
