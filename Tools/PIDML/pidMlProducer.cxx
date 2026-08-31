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

/// \file pidMlProducer.cxx
/// \brief Produce PID ML skimmed data from MC or data files.
///
/// \author Maja Kabus <mkabus@cern.ch>
/// \author Marek Mytkowski <marek.mytkowski@cern.ch>

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/PIDML/pidMl.h"
#include "Tools/PIDML/pidUtils.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <string_view>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace pidml::pidutils;

// Naming convention
//  Data: experimental data without simulation
//  MC: experimental data with Monte Carlo simulation
//  ML: only columns used by Machine Learning network
struct PidMlProducer {
  Produces<aod::PidTracksDataMl> pidTracksTableDataML;
  Produces<aod::PidTracksData> pidTracksTableData;
  Produces<aod::PidTracksMcMl> pidTracksTableMCML;
  Produces<aod::PidTracksMc> pidTracksTableMC;

  // Basic event selection: sel8() + (zVertex < zVertexMax) + (centMin <= FT0C <= centMax)
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "confBasicEventSelection";

    Configurable<bool> enableSel8{"enableSel8", true, "enable sel8 cut."};
    Configurable<float> zVertexMax{"zVertexMax", 10., "Maximum Z Vertex value cm."};
    Configurable<float> centMin{"centMin", 0.f, "Minimum centrality (FT0C)."};
    Configurable<float> centMax{"centMax", 200.f, "Maximum centrality (FT0C)."};
  } confBasicEventSelection;


  // Strict event selection: criteria boolean flags (see their details) + TPC occupancy
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "confStrictEventSelection";
    Configurable<bool> enable{"enable", false, "Enable stricter event selection filtering"};

    struct : o2::framework::ConfigurableGroup {
      std::string prefix = "confStrictEventSelection.criteria";
      Configurable<bool> evIsGoodZvtxFT0vsPV{"evIsGoodZvtxFT0vsPV", true, "Require kIsGoodZvtxFT0vsPV selection on events"};
      Configurable<bool> evNoSameBunchPileup{"evNoSameBunchPileup", true, "Require kNoSameBunchPileup selection on events"};
      Configurable<bool> evIsVertexITSTPC{"evIsVertexITSTPC", true, "Require kIsVertexITSTPC selection on events"};
      Configurable<bool> isGoodITSLayersAll{"isGoodITSLayersAll", true, "Require kIsGoodITSLayersAll selection on events"};
      Configurable<bool> noITSROFrameBorder{"noITSROFrameBorder", true, "Require kNoITSROFrameBorder selection on events"};
      Configurable<bool> noTimeFrameBorder{"noTimeFrameBorder", true, "Require kNoTimeFrameBorder selection on events"};
      Configurable<bool> noCollInRofStandard{"noCollInRofStandard", true, "Require kNoCollInRofStandard selection on events"};
      Configurable<bool> noHighMultCollInPrevRof{"noHighMultCollInPrevRof", true, "Require kNoHighMultCollInPrevRof selection on events"};
      Configurable<bool> noCollInTimeRangeStandard{"noCollInTimeRangeStandard", true, "Require kNoCollInTimeRangeStandard selection on events"};
      Configurable<int> tpcOccupancyMin{"tpcOccupancyMin", 0, "Minimum value for TPC occupancy selection"};
      Configurable<int> tpcOccupancyMax{"tpcOccupancyMax", 500, "Maximum value for TPC occupancy selection"};
    } criteria;
  } confStrictEventSelection;

  Filter trackFilter = requireGlobalTrackInFilter();

  // Data tracks
  using BigTracksDataML = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;
  using BigTracksData = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal>>;

  // MC tracks
  using BigTracksMCML = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;
  using BigTracksMC = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;

  using MyCollisionML = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::Mults, aod::EvSels>::iterator;
  using MyCollision = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::Mults, aod::EvSels>::iterator;

  static constexpr uint32_t NCharges = 2;

  static constexpr std::string_view HistPrefixes[NCharges] = {"minus", "plus"};

  // 2D
  std::array<std::shared_ptr<TH2>, NCharges> hTPCSigvsP;
  std::array<std::shared_ptr<TH2>, NCharges> hTOFBetavsP;
  std::array<std::shared_ptr<TH2>, NCharges> hTOFSigvsP;
  std::array<std::shared_ptr<TH2>, NCharges> hFilteredTOFSigvsP;
  std::array<std::shared_ptr<TH2>, NCharges> hTRDPattvsP;
  std::array<std::shared_ptr<TH2>, NCharges> hTRDSigvsP;

  // 1D
  std::array<std::shared_ptr<TH1>, NCharges> hP;
  std::array<std::shared_ptr<TH1>, NCharges> hPt;
  std::array<std::shared_ptr<TH1>, NCharges> hPx;
  std::array<std::shared_ptr<TH1>, NCharges> hPy;
  std::array<std::shared_ptr<TH1>, NCharges> hPz;
  std::array<std::shared_ptr<TH1>, NCharges> hX;
  std::array<std::shared_ptr<TH1>, NCharges> hY;
  std::array<std::shared_ptr<TH1>, NCharges> hZ;
  std::array<std::shared_ptr<TH1>, NCharges> hAlpha;
  std::array<std::shared_ptr<TH1>, NCharges> hTrackType;
  std::array<std::shared_ptr<TH1>, NCharges> hTPCNClsShared;
  std::array<std::shared_ptr<TH1>, NCharges> hDcaXY;
  std::array<std::shared_ptr<TH1>, NCharges> hDcaZ;
  std::array<std::shared_ptr<TH1>, NCharges> hPdgCode;
  std::array<std::shared_ptr<TH1>, NCharges> hIsPrimary;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  static const char* genTagvsP(const char* name)
  {
    return Form("%s vs  #it{p};#it{p} (GeV/#it{c});%s", name, name);
  }

  template <uint32_t prefixInd>
  void initHistSign()
  {
    hTPCSigvsP[prefixInd] = registry.add<TH2>(Form("%s/hTPCSigvsP", HistPrefixes[prefixInd].data()), genTagvsP("TPC Signal"), HistType::kTH2F, {{500, 0., 10.}, {1000, 0., 600.}});
    hTOFBetavsP[prefixInd] = registry.add<TH2>(Form("%s/hTOFBetavsP", HistPrefixes[prefixInd].data()), genTagvsP("TOF beta"), HistType::kTH2F, {{500, 0., 10.}, {6000, -3., 3.}});
    hTOFSigvsP[prefixInd] = registry.add<TH2>(Form("%s/hTOFSigvsP", HistPrefixes[prefixInd].data()), genTagvsP("TOF signal"), HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}});
    hFilteredTOFSigvsP[prefixInd] = registry.add<TH2>(Form("%s/filtered/hTOFSigvsP", HistPrefixes[prefixInd].data()), genTagvsP("TOF signal (filtered)"), HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}});
    hTRDPattvsP[prefixInd] = registry.add<TH2>(Form("%s/hTRDPattvsP", HistPrefixes[prefixInd].data()), genTagvsP("TRD pattern"), HistType::kTH2F, {{500, 0., 10.}, {110, -10., 100.}});
    hTRDSigvsP[prefixInd] = registry.add<TH2>(Form("%s/hTRDSigvsP", HistPrefixes[prefixInd].data()), genTagvsP("TRD signal"), HistType::kTH2F, {{500, 0., 10.}, {2500, -2., 100.}});
    hP[prefixInd] = registry.add<TH1>(Form("%s/hP", HistPrefixes[prefixInd].data()), "#it{p};#it{p} (GeV/#it{c})", HistType::kTH1F, {{500, 0., 6.}});
    hPt[prefixInd] = registry.add<TH1>(Form("%s/hPt", HistPrefixes[prefixInd].data()), "#it{p}_{t};#it{p}_{t} (GeV/#it{c})", HistType::kTH1F, {{500, 0., 6.}});
    hPx[prefixInd] = registry.add<TH1>(Form("%s/hPx", HistPrefixes[prefixInd].data()), "#it{p}_{x};#it{p}_{x} (GeV/#it{c})", HistType::kTH1F, {{1000, -6., 6.}});
    hPy[prefixInd] = registry.add<TH1>(Form("%s/hPy", HistPrefixes[prefixInd].data()), "#it{p}_{y};#it{p}_{y} (GeV/#it{c})", HistType::kTH1F, {{1000, -6., 6.}});
    hPz[prefixInd] = registry.add<TH1>(Form("%s/hPz", HistPrefixes[prefixInd].data()), "#it{p}_{z};#it{p}_{z} (GeV/#it{c})", HistType::kTH1F, {{1000, -6., 6.}});
    hX[prefixInd] = registry.add<TH1>(Form("%s/hX", HistPrefixes[prefixInd].data()), "#it{x};#it{x}", HistType::kTH1F, {{1000, -2., 2.}});
    hY[prefixInd] = registry.add<TH1>(Form("%s/hY", HistPrefixes[prefixInd].data()), "#it{y};#it{y}", HistType::kTH1F, {{1000, -2., 2.}});
    hZ[prefixInd] = registry.add<TH1>(Form("%s/hZ", HistPrefixes[prefixInd].data()), "#it{z};#it{z}", HistType::kTH1F, {{1000, -10., 10.}});
    hAlpha[prefixInd] = registry.add<TH1>(Form("%s/hAlpha", HistPrefixes[prefixInd].data()), "alpha;alpha", HistType::kTH1F, {{1000, -5., 5.}});
    hTrackType[prefixInd] = registry.add<TH1>(Form("%s/hTrackType", HistPrefixes[prefixInd].data()), "Track Type;Track Type", HistType::kTH1F, {{300, 0., 300.}});
    hTPCNClsShared[prefixInd] = registry.add<TH1>(Form("%s/hTPCNClsShared", HistPrefixes[prefixInd].data()), "hTPCNClsShared;hTPCNClsShared", HistType::kTH1F, {{100, 0., 100.}});
    hDcaXY[prefixInd] = registry.add<TH1>(Form("%s/hDcaXY", HistPrefixes[prefixInd].data()), "#it{DcaXY};#it{DcaXY}", HistType::kTH1F, {{1000, -1., 1.}});
    hDcaZ[prefixInd] = registry.add<TH1>(Form("%s/hDcaZ", HistPrefixes[prefixInd].data()), "#it{DcaZ};#it{DcaZ}", HistType::kTH1F, {{1000, -1., 1.}});
  }

  template <uint32_t prefixInd>
  void initHistSignMC()
  {
    initHistSign<prefixInd>();
    hPdgCode[prefixInd] = registry.add<TH1>(Form("%s/hPdgCode", HistPrefixes[prefixInd].data()), "#it{PdgCode};#it{PdgCode}", HistType::kTH1F, {{2500, 0., 2500.}});
    hIsPrimary[prefixInd] = registry.add<TH1>(Form("%s/hIsPrimary", HistPrefixes[prefixInd].data()), "#it{IsPrimary};#it{IsPrimary}", HistType::kTH1F, {{4, -0.5, 1.5}});
  }

  template <uint32_t prefixInd, typename T>
  void fillHistSign(const T& track)
  {
    hTPCSigvsP[prefixInd]->Fill(track.p(), track.tpcSignal());
    hTOFBetavsP[prefixInd]->Fill(track.p(), track.beta());
    hTOFSigvsP[prefixInd]->Fill(track.p(), track.tofSignal());
    if (tofMissing(track)) {
      hFilteredTOFSigvsP[prefixInd]->Fill(track.p(), kTOFMissingSignal);
    } else {
      hFilteredTOFSigvsP[prefixInd]->Fill(track.p(), track.tofSignal());
    }
    hTRDPattvsP[prefixInd]->Fill(track.p(), track.trdPattern());
    hTRDSigvsP[prefixInd]->Fill(track.p(), track.trdSignal());
    hP[prefixInd]->Fill(track.p());
    hPt[prefixInd]->Fill(track.pt());
    hPx[prefixInd]->Fill(track.px());
    hPy[prefixInd]->Fill(track.py());
    hPz[prefixInd]->Fill(track.pz());
    hX[prefixInd]->Fill(track.x());
    hY[prefixInd]->Fill(track.y());
    hZ[prefixInd]->Fill(track.z());
    hAlpha[prefixInd]->Fill(track.alpha());
    hTrackType[prefixInd]->Fill(track.trackType());
    hTPCNClsShared[prefixInd]->Fill(track.tpcNClsShared());
    hDcaXY[prefixInd]->Fill(track.dcaXY());
    hDcaZ[prefixInd]->Fill(track.dcaZ());
  }

  template <uint32_t prefixInd, typename T>
  void fillHistSignMC(const T& track, uint32_t pdgCode, uint8_t isPrimary)
  {
    fillHistSign<prefixInd>(track);
    hPdgCode[prefixInd]->Fill(pdgCode);
    hIsPrimary[prefixInd]->Fill(isPrimary);
  }

  template <typename T>
  void fillHistMC(const T& track, uint32_t pdgCode, uint8_t isPrimary)
  {
    if (track.sign() < 0) {
      fillHistSignMC<0>(track, pdgCode, isPrimary);
    } else {
      fillHistSignMC<1>(track, pdgCode, isPrimary);
    }
  }

  template <typename T>
  void fillHist(const T& track)
  {
    if (track.sign() < 0) {
      fillHistSign<0>(track);
    } else {
      fillHistSign<1>(track);
    }
  }

  void init(InitContext&)
  {
    if (doprocessMcMl || doprocessMcAll) {
      initHistSignMC<0>();
      initHistSignMC<1>();
    } else {
      initHistSign<0>();
      initHistSign<1>();
    }
  }

  template <typename T>
  float getTOFSignal(T const& track)
  {
    return tofMissing(track) ? std::numeric_limits<float>::quiet_NaN() : track.tofSignal();
  }

  template <typename T>
  float getTOFBeta(T const& track)
  {
    return tofMissing(track) ? std::numeric_limits<float>::quiet_NaN() : track.beta();
  }

  template <typename T>
  float getTRDSignal(T const& track)
  {
    return trdMissing(track) ? std::numeric_limits<float>::quiet_NaN() : track.trdSignal();
  }

  template <typename T>
  uint8_t getTRDPattern(T const& track)
  {
    return trdMissing(track) ? static_cast<uint8_t>(0U) : track.trdPattern();
  }

  // inspired by PWGCF/FemtoUniverse/TableProducer/femtoUniverseProducerTask.cxx
  template <bool isMC, typename CollisionType>
  bool fillCollisionsCentRun3(CollisionType const& col)
  {
    // Basic Event Selection
    {
      const auto& basic_criteria = confBasicEventSelection;

      // Z Vertex filter
      if (std::abs(col.posZ()) > basic_criteria.zVertexMax.value)
        return false;

      // Selection 8 filters
      if (basic_criteria.enableSel8.value && !col.sel8())
        return false;

      // Centrality filter
      const auto cent = col.centFT0C();
      if ((cent < basic_criteria.centMin.value) || (cent > basic_criteria.centMax.value))
        return false;
    }

    // Strict Event Selection
    {
      if (!confStrictEventSelection.enable)
        return true;

      const auto& criteria = confStrictEventSelection.criteria;

      // TPC Occupancy filter
      const auto occupancy = col.trackOccupancyInTimeRange();
      if (occupancy < criteria.tpcOccupancyMin.value || occupancy > criteria.tpcOccupancyMax.value)
        return false;

      // Event Selection boolean filters
      if ((!criteria.evNoSameBunchPileup.value || col.selection_bit(aod::evsel::kNoSameBunchPileup)) &&
          (!criteria.evIsGoodZvtxFT0vsPV.value || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) &&
          (!criteria.isGoodITSLayersAll.value || col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) &&
          (!criteria.noCollInRofStandard.value || col.selection_bit(aod::evsel::kNoCollInRofStandard)) &&
          (!criteria.noHighMultCollInPrevRof.value || col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof)) &&
          (!criteria.evIsVertexITSTPC.value || col.selection_bit(aod::evsel::kIsVertexITSTPC)) &&
          (!criteria.noCollInTimeRangeStandard.value || col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) &&
          (!criteria.noITSROFrameBorder.value || col.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
          (!criteria.noTimeFrameBorder.value || col.selection_bit(aod::evsel::kNoTimeFrameBorder)))
        return true;

      return false;
    }
  }

  void processDataML(MyCollisionML const& collision, BigTracksDataML const& tracks)
  {
    if (!fillCollisionsCentRun3<false>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      pidTracksTableDataML(track.collisionId(), track.tpcSignal(), getTRDSignal(track), getTRDPattern(track),
                           getTOFSignal(track), getTOFBeta(track),
                           track.p(), track.pt(), track.px(), track.py(), track.pz(),
                           track.sign(),
                           track.x(), track.y(), track.z(),
                           track.alpha(),
                           track.trackType(),
                           track.tpcNClsShared(),
                           track.dcaXY(), track.dcaZ());

      fillHist(track);
    }
  }
  PROCESS_SWITCH(PidMlProducer, processDataML, "Produce only ML real data", true);

  void processDataAll(MyCollision const& collision, BigTracksData const& tracks)
  {
    if (!fillCollisionsCentRun3<false>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      pidTracksTableData(track.collisionId(), collision.multFV0A(), collision.multFV0C(), collision.multFV0M(),
                         collision.multFT0A(), collision.multFT0C(), collision.multFT0M(),
                         collision.multZNA(), collision.multZNC(), collision.multTracklets(),
                         collision.multTPC(), collision.multNTracksPV(),
                         collision.posZ(), collision.trackOccupancyInTimeRange(), collision.centFT0C(),
                         collision.selection_bit(aod::evsel::kNoSameBunchPileup),
                         collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV),
                         collision.selection_bit(aod::evsel::kIsGoodITSLayersAll),
                         collision.selection_bit(aod::evsel::kNoCollInRofStandard),
                         collision.selection_bit(aod::evsel::kNoHighMultCollInPrevRof),
                         collision.selection_bit(aod::evsel::kIsVertexITSTPC),
                         collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard),
                         collision.selection_bit(aod::evsel::kNoITSROFrameBorder),
                         collision.selection_bit(aod::evsel::kNoTimeFrameBorder),
                         track.tpcSignal(), getTRDSignal(track), getTRDPattern(track),
                         track.trackEtaEmcal(), track.trackPhiEmcal(),
                         getTOFSignal(track), getTOFBeta(track),
                         track.p(), track.pt(), track.px(), track.py(), track.pz(),
                         track.sign(),
                         track.x(), track.y(), track.z(),
                         track.alpha(),
                         track.trackType(),
                         track.tpcNClsShared(),
                         track.dcaXY(), track.dcaZ(),
                         track.tpcNSigmaEl(), track.tpcExpSigmaEl(), track.tpcExpSignalDiffEl(),
                         track.tofNSigmaEl(), track.tofExpSigmaEl(), track.tofExpSignalDiffEl(),
                         track.tpcNSigmaMu(), track.tpcExpSigmaMu(), track.tpcExpSignalDiffMu(),
                         track.tofNSigmaMu(), track.tofExpSigmaMu(), track.tofExpSignalDiffMu(),
                         track.tpcNSigmaPi(), track.tpcExpSigmaPi(), track.tpcExpSignalDiffPi(),
                         track.tofNSigmaPi(), track.tofExpSigmaPi(), track.tofExpSignalDiffPi(),
                         track.tpcNSigmaKa(), track.tpcExpSigmaKa(), track.tpcExpSignalDiffKa(),
                         track.tofNSigmaKa(), track.tofExpSigmaKa(), track.tofExpSignalDiffKa(),
                         track.tpcNSigmaPr(), track.tpcExpSigmaPr(), track.tpcExpSignalDiffPr(),
                         track.tofNSigmaPr(), track.tofExpSigmaPr(), track.tofExpSignalDiffPr(),
                         track.snp(), track.tgl(), track.sigmaY(), track.sigmaZ(), track.sigmaSnp(),
                         track.sigmaTgl(), track.sigma1Pt(), track.tpcNClsFindable(), track.tpcNClsFound(),
                         track.tpcNClsPID(), track.tpcChi2NCl(), track.trackTime(), track.length(),
                         track.tofChi2(), track.itsClusterMap(), track.trdChi2());

      fillHist(track);
    }
  }
  PROCESS_SWITCH(PidMlProducer, processDataAll, "Produce all real data", false);

  void processMcMl(MyCollisionML const& collision, BigTracksMCML const& tracks, aod::McParticles const& /*mctracks*/)
  {
    if (!fillCollisionsCentRun3<true>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      uint8_t isPrimary = static_cast<uint8_t>(mcParticle.isPhysicalPrimary());
      uint32_t pdgCode = mcParticle.pdgCode();
      pidTracksTableMCML(track.collisionId(), track.tpcSignal(), getTRDSignal(track), getTRDPattern(track),
                         getTOFSignal(track), getTOFBeta(track),
                         track.p(), track.pt(), track.px(), track.py(), track.pz(),
                         track.sign(),
                         track.x(), track.y(), track.z(),
                         track.alpha(),
                         track.trackType(),
                         track.tpcNClsShared(),
                         track.dcaXY(), track.dcaZ(),
                         pdgCode,
                         isPrimary);

      fillHistMC(track, pdgCode, isPrimary);
    }
  }
  PROCESS_SWITCH(PidMlProducer, processMcMl, "Produce only ML MC essential data", false);

  void processMcAll(MyCollision const& collision, BigTracksMC const& tracks, aod::McParticles const& /*mctracks*/)
  {
    if (!fillCollisionsCentRun3<true>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      uint8_t isPrimary = static_cast<uint8_t>(mcParticle.isPhysicalPrimary());
      uint32_t pdgCode = mcParticle.pdgCode();
      pidTracksTableMC(track.collisionId(), collision.multFV0A(), collision.multFV0C(), collision.multFV0M(),
                       collision.multFT0A(), collision.multFT0C(), collision.multFT0M(),
                       collision.multZNA(), collision.multZNC(), collision.multTracklets(),
                       collision.multTPC(), collision.multNTracksPV(),
                       collision.posZ(), collision.trackOccupancyInTimeRange(), collision.centFT0C(),
                       collision.selection_bit(aod::evsel::kNoSameBunchPileup),
                       collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV),
                       collision.selection_bit(aod::evsel::kIsGoodITSLayersAll),
                       collision.selection_bit(aod::evsel::kNoCollInRofStandard),
                       collision.selection_bit(aod::evsel::kNoHighMultCollInPrevRof),
                       collision.selection_bit(aod::evsel::kIsVertexITSTPC),
                       collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard),
                       collision.selection_bit(aod::evsel::kNoITSROFrameBorder),
                       collision.selection_bit(aod::evsel::kNoTimeFrameBorder),
                       track.tpcSignal(), getTRDSignal(track), getTRDPattern(track),
                       track.trackEtaEmcal(), track.trackPhiEmcal(),
                       getTOFSignal(track), getTOFBeta(track),
                       track.p(), track.pt(), track.px(), track.py(), track.pz(),
                       track.sign(),
                       track.x(), track.y(), track.z(),
                       track.alpha(),
                       track.trackType(),
                       track.tpcNClsShared(),
                       track.dcaXY(), track.dcaZ(),
                       track.tpcNSigmaEl(), track.tpcExpSigmaEl(), track.tpcExpSignalDiffEl(),
                       track.tofNSigmaEl(), track.tofExpSigmaEl(), track.tofExpSignalDiffEl(),
                       track.tpcNSigmaMu(), track.tpcExpSigmaMu(), track.tpcExpSignalDiffMu(),
                       track.tofNSigmaMu(), track.tofExpSigmaMu(), track.tofExpSignalDiffMu(),
                       track.tpcNSigmaPi(), track.tpcExpSigmaPi(), track.tpcExpSignalDiffPi(),
                       track.tofNSigmaPi(), track.tofExpSigmaPi(), track.tofExpSignalDiffPi(),
                       track.tpcNSigmaKa(), track.tpcExpSigmaKa(), track.tpcExpSignalDiffKa(),
                       track.tofNSigmaKa(), track.tofExpSigmaKa(), track.tofExpSignalDiffKa(),
                       track.tpcNSigmaPr(), track.tpcExpSigmaPr(), track.tpcExpSignalDiffPr(),
                       track.tofNSigmaPr(), track.tofExpSigmaPr(), track.tofExpSignalDiffPr(),
                       track.snp(), track.tgl(), track.sigmaY(), track.sigmaZ(), track.sigmaSnp(),
                       track.sigmaTgl(), track.sigma1Pt(), track.tpcNClsFindable(), track.tpcNClsFound(),
                       track.tpcNClsPID(), track.tpcChi2NCl(), track.trackTime(), track.length(),
                       track.tofChi2(), track.itsClusterMap(), track.trdChi2(),
                       pdgCode,
                       isPrimary);

      fillHistMC(track, pdgCode, isPrimary);
    }
  }
  PROCESS_SWITCH(PidMlProducer, processMcAll, "Produce all MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidMlProducer>(cfgc)};
}
