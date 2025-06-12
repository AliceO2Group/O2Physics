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

/// \file corrSparse.cxx
/// \brief Provides a sparse with usefull two particle correlation info
/// \author Thor Jensen (thor.kjaersgaard.jensen@cern.ch)

#include <CCDB/BasicCCDBManager.h>
#include "TRandom3.h"
#include <vector>

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
namespace corrsparse
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, int);
}
DECLARE_SOA_TABLE(Multiplicity, "AOD", "MULTIPLICITY",
                  corrsparse::Multiplicity);

} // namespace o2::aod

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct CorrSparse {
  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMinMult, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgMaxMult, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgMergingCut, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgApplyTwoTrackEfficiency, bool, true, "Apply two track efficiency for tpc tpc")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(etaMftTrackMin, float, -3.6, "Minimum eta for MFT track")
  O2_DEFINE_CONFIGURABLE(etaMftTrackMax, float, -2.5, "Maximum eta for MFT track")
  O2_DEFINE_CONFIGURABLE(nClustersMftTrack, int, 5, "Minimum number of clusters for MFT track")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")

  Configurable<bool> processMFT{"processMFT", true, "Associate particle from MFT"};

  SliceCache cache;
  SliceCache cacheNch;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis vtxMix{"vtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis multMix{"multMix", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZVtxCut) && (aod::evsel::sel8) == true;
  Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (cfgPtCutMin < aod::track::pt) && (cfgPtCutMax > aod::track::pt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  // Define the outputs
  OutputObj<CorrelationContainer> same{Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixed{Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};

  HistogramRegistry registry{"registry"};

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel>>; // aod::CentFT0Cs
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  void init(InitContext&)
  {
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

  TRandom3* gRandom = new TRandom3();

  template <typename TTrackAssoc>
  bool isAcceptedMftTrack(TTrackAssoc const& mftTrack)
  {
    // cut on the eta of MFT tracks
    if (mftTrack.eta() < etaMftTrackMin || mftTrack.eta() > etaMftTrackMax) {
      return false;
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < nClustersMftTrack) {
      return false;
    }

    return true;
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

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {

    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    for (auto const& track1 : tracks) {

      if (processMFT) {
        if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) {
          if (!isAcceptedMftTrack(track1)) {
            continue;
          }
        }
      }

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
  void fillCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (system == SameEvent) {
        registry.fill(HIST("Nch_used"), tracks1.size());
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt());
      }

      for (auto const& track2 : tracks2) {

        if (processMFT) {
          if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {
            if (!isAcceptedMftTrack(track2)) {
              continue;
            }
          }
        } else {
          if (track1.pt() <= track2.pt())
            continue; // skip if the trigger pt is less than the associate pt
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (cfgApplyTwoTrackEfficiency && std::abs(deltaEta) < cfgMergingCut) {

          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgMergingCut;

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

  void processSame(AodCollisions::iterator const& collision, AodTracks const& tracks, aod::MFTTracks const& mfts, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    if (processMFT) {
      fillYield(collision, mfts);

      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, mfts, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()));

    } else {
      fillYield(collision, tracks);

      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrSparse, processSame, "Process same event", true);

  // the process for filling the mixed events
  void processMixed(AodCollisions const& collisions, AodTracks const& tracks, aod::MFTTracks const& MFTtracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](AodCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {vtxMix, multMix}, true};

    if (processMFT) {

      auto tracksTuple = std::make_tuple(tracks, MFTtracks);
      Pair<AodCollisions, AodTracks, aod::MFTTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto const& [collision1, tracks1, collision2, tracks2] : pair) {
        registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
        auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

        if ((tracks1.size() < cfgMinMult || tracks1.size() >= cfgMaxMult))
          continue;

        fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()));
      }
    } else {
      auto tracksTuple = std::make_tuple(tracks, tracks);
      Pair<AodCollisions, AodTracks, AodTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto const& [collision1, tracks1, collision2, tracks2] : pair) {
        registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
        auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

        if ((tracks1.size() < cfgMinMult || tracks1.size() >= cfgMaxMult))
          continue;

        if ((tracks2.size() < cfgMinMult || tracks2.size() >= cfgMaxMult))
          continue;

        fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()));
      }
    }
  }

  PROCESS_SWITCH(CorrSparse, processMixed, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrSparse>(cfgc),
  };
}
