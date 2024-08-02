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

/// \file pidMLProducerMc.cxx
/// \brief Produce PID ML skimmed data from MC files.
///
/// \author Maja Kabus <mkabus@cern.ch>

#include <string_view>
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/PIDML/pidML.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

// Naming convention
//  Data: experimental data without simulation
//  MC: experimental data with Monte Carlo simulation
//  ML: only columns used by Machine Learning network
struct PidMlProducer {
  Produces<aod::PidTracksDataMl> pidTracksTableDataML;
  Produces<aod::PidTracksData> pidTracksTableData;
  Produces<aod::PidTracksMcMl> pidTracksTableMCML;
  Produces<aod::PidTracksMc> pidTracksTableMC;

  Filter trackFilter = requireGlobalTrackInFilter();

  // Data tracks
  using BigTracksDataML = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;
  using BigTracksData = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal>>;

  // MC tracks
  using BigTracksMCML = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;
  using BigTracksMC = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels>>;

  using MyCollisionML = aod::Collisions::iterator;
  using MyCollision = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>::iterator;

  static constexpr float kEps = 1e-10f;
  static constexpr float kMissingBeta = -999.0f;
  static constexpr float kMissingTOFSignal = -999.0f;

  static constexpr std::string_view histPrefixes = {"minus", "plus"};

  HistogramRegistry registry{"registry", {}};

  static const char* genTagvsP(const char* name)
  {
    return Form("%s vs  #it{p};#it{p} (GeV/#it{c});%s", name, name);
  }

  template <int32_t prefixInd>
  void initHistSign()
  {
    registry.add<TH2F>(Form("%s/hTPCSigvsP", histPrefixes[prefixInd]), genTagvsP("TPC Signal"), {HistType::kTH2F, {{500, 0., 10.}, {1000, 0., 600.}}});
    registry.add<TH2F>(Form("%s/hTOFBetavsP", histPrefixes[prefixInd]), genTagvsP("TOF beta"), {HistType::kTH2F, {{500, 0., 10.}, {6000, -3., 3.}}});
    registry.add<TH2F>(Form("%s/hTOFSigvsP", histPrefixes[prefixInd]), genTagvsP("TOF signal"), {HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}}});
    registry.add<TH2F>(Form("%s/filtered/hTOFSigvsP", histPrefixes[prefixInd]), genTagvsP("TOF signal (filtered)"), {HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}}});
    registry.add<TH2F>(Form("%s/hTRDPattvsP", histPrefixes[prefixInd]), genTagvsP("TRD pattern"), {HistType::kTH2F, {{500, 0., 10.}, {110, -10., 100.}}});
    registry.add<TH2F>(Form("%s/hTRDSigvsP", histPrefixes[prefixInd]), genTagvsP("TRD signal"), {HistType::kTH2F, {{500, 0., 10.}, {2500, -2., 100.}}});
    registry.add<TH1F>(Form("%s/hP", histPrefixes[prefixInd]), "#it{p};#it{p} (GeV/#it{c})", {HistType::kTH1F, {{500, 0., 6.}}});
    registry.add<TH1F>(Form("%s/hPt", histPrefixes[prefixInd]), "#it{p}_{t};#it{p}_{t} (GeV/#it{c})", {HistType::kTH1F, {{500, 0., 6.}}});
    registry.add<TH1F>(Form("%s/hPx", histPrefixes[prefixInd]), "#it{p}_{x};#it{p}_{x} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}});
    registry.add<TH1F>(Form("%s/hPy", histPrefixes[prefixInd]), "#it{p}_{y};#it{p}_{y} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}});
    registry.add<TH1F>(Form("%s/hPz", histPrefixes[prefixInd]), "#it{p}_{z};#it{p}_{z} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}});
    registry.add<TH1F>(Form("%s/hX", histPrefixes[prefixInd]), "#it{x};#it{x}", {HistType::kTH1F, {{1000, -2., 2.}}});
    registry.add<TH1F>(Form("%s/hY", histPrefixes[prefixInd]), "#it{y};#it{y}", {HistType::kTH1F, {{1000, -2., 2.}}});
    registry.add<TH1F>(Form("%s/hZ", histPrefixes[prefixInd]), "#it{z};#it{z}", {HistType::kTH1F, {{1000, -10., 10.}}});
    registry.add<TH1F>(Form("%s/hAlpha", histPrefixes[prefixInd]), "#{alpha};#{alpha}", {HistType::kTH1F, {{1000, -5., 5.}}});
    registry.add<TH1F>(Form("%s/hTrackType", histPrefixes[prefixInd]), "Track Type;Track Type", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add<TH1F>(Form("%s/hTPCNClsShared", histPrefixes[prefixInd]), "hTPCNClsShared;hTPCNClsShared", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add<TH1F>(Form("%s/hDcaXY", histPrefixes[prefixInd]), "#it{DcaXY};#it{DcaXY}", {HistType::kTH1F, {{1000, -1., 1.}}});
    registry.add<TH1F>(Form("%s/hDcaZ", histPrefixes[prefixInd]), "#it{DcaZ};#it{DcaZ}", {HistType::kTH1F, {{1000, -1., 1.}}});
  }

  template <int32_t prefixInd>
  void initHistSignMC()
  {
    initHistSign<prefixInd>();
    registry.add<TH1F>(Form("%s/hPdgCode", histPrefixes[prefixInd]), "#it{PdgCode};#it{PdgCode}", {HistType::kTH1F, {{2500, 0., 2500.}}});
    registry.add<TH1F>(Form("%s/hIsPrimary", histPrefixes[prefixInd]), "#it{IsPrimary};#it{IsPrimary}", {HistType::kTH1F, {{4, -0.5, 1.5}}});
  }

  template <int32_t prefixInd, typename T>
  void fillHistSign(const T& track)
  {
    registry.fill(HIST(Form("%s/hTPCSigvsP", histPrefixes[prefixInd])), track.p(), track.tpcSignal());
    registry.fill(HIST(Form("%s/hTOFBetavsP", histPrefixes[prefixInd])), track.p(), track.beta());
    registry.fill(HIST(Form("%s/hTOFSigvsP", histPrefixes[prefixInd])), track.p(), track.tofSignal());
    if (TMath::Abs(track.beta() - kMissingBeta) >= kEps) {
      registry.fill(HIST(Form("%s/filtered/hTOFSigvsP", histPrefixes[prefixInd])), track.p(), track.tofSignal());
    } else {
      registry.fill(HIST(Form("%s/filtered/hTOFSigvsP", histPrefixes[prefixInd])), track.p(), kMissingTOFSignal);
    }
    registry.fill(HIST(Form("%s/hTRDPattvsP", histPrefixes[prefixInd])), track.p(), track.trdPattern());
    registry.fill(HIST(Form("%s/hTRDSigvsP", histPrefixes[prefixInd])), track.p(), track.trdSignal());
    registry.fill(HIST(Form("%s/hP", histPrefixes[prefixInd])), track.p());
    registry.fill(HIST(Form("%s/hPt", histPrefixes[prefixInd])), track.pt());
    registry.fill(HIST(Form("%s/hPx", histPrefixes[prefixInd])), track.px());
    registry.fill(HIST(Form("%s/hPy", histPrefixes[prefixInd])), track.py());
    registry.fill(HIST(Form("%s/hPz", histPrefixes[prefixInd])), track.pz());
    registry.fill(HIST(Form("%s/hX", histPrefixes[prefixInd])), track.x());
    registry.fill(HIST(Form("%s/hY", histPrefixes[prefixInd])), track.y());
    registry.fill(HIST(Form("%s/hZ", histPrefixes[prefixInd])), track.z());
    registry.fill(HIST(Form("%s/hAlpha", histPrefixes[prefixInd])), track.alpha());
    registry.fill(HIST(Form("%s/hTrackType", histPrefixes[prefixInd])), track.trackType());
    registry.fill(HIST(Form("%s/hTPCNClsShared", histPrefixes[prefixInd])), track.tpcNClsShared());
    registry.fill(HIST(Form("%s/hDcaXY", histPrefixes[prefixInd])), track.dcaXY());
    registry.fill(HIST(Form("%s/hDcaZ", histPrefixes[prefixInd])), track.dcaZ());
  }

  template <int32_t prefixInd, typename T>
  void fillHistSignMC(const T& track, uint32_t pdgCode, uint8_t isPrimary)
  {
    fillHistSign<prefixInd>(const T& track);
    registry.fill(HIST(Form("%s/hPdgCode", histPrefixes[prefixInd])), pdgCode);
    registry.fill(HIST(Form("%s/hIsPrimary", histPrefixes[prefixInd])), isPrimary);
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

  void initHistosMC()
  {
    static_for<0, 1>([&](auto prefixInd) {
      initHistSignMC<prefixInd>();
    })
  }

  void init(InitContext&)
  {
    if (doProcessMcMl || doProcessMcAll) {
      static_for<0, 1>([&](auto prefixInd) {
        initHistSignMC<prefixInd>();
      });
    } else {
      static_for<0, 1>([&](auto prefixInd) {
        initHistSign<prefixInd>();
      });
    }
  }

  void processDataML(MyCollisionML const& /*collision*/, BigTracksML const& tracks)
  {
    for (const auto& track : tracks) {
      pidTracksTableDataML(track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                           track.tofSignal(), track.beta(),
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
  PROCESS_SWITCH(PidMlProducerData, processDataML, "Produce only ML real data", true);

  void processDataAll(MyCollision const& collision, BigTracks const& tracks)
  {
    for (const auto& track : tracks) {
      pidTracksTableData(collision.centRun2V0M(),
                         collision.multFV0A(), collision.multFV0C(), collision.multFV0M(),
                         collision.multFT0A(), collision.multFT0C(), collision.multFT0M(),
                         collision.multZNA(), collision.multZNC(),
                         collision.multTracklets(), collision.multTPC(),
                         track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                         track.trackEtaEmcal(), track.trackPhiEmcal(),
                         track.tofSignal(), track.beta(),
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
                         track.tofNSigmaPr(), track.tofExpSigmaPr(), track.tofExpSignalDiffPr());

      fillHist(track);
    }
  }
  PROCESS_SWITCH(PidMlProducerData, processDataAll, "Produce all real data", false);

  void processMcMl(MyCollisionML const& /*collision*/, BigTracksMCML const& tracks, aod::McParticles const& /*mctracks*/)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      uint8_t isPrimary = static_cast<uint8_t>(mcParticle.isPhysicalPrimary());
      uint32_t pdgCode = mcParticle.pdgCode();
      pidTracksTableMCML(track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                         track.tofSignal(), track.beta(),
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
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      uint8_t isPrimary = static_cast<uint8_t>(mcParticle.isPhysicalPrimary());
      uint32_t pdgCode = mcParticle.pdgCode();
      pidTracksTableMC(collision.centRun2V0M(),
                       collision.multFV0A(), collision.multFV0C(), collision.multFV0M(),
                       collision.multFT0A(), collision.multFT0C(), collision.multFT0M(),
                       collision.multZNA(), collision.multZNC(),
                       collision.multTracklets(), collision.multTPC(),
                       track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                       track.trackEtaEmcal(), track.trackPhiEmcal(),
                       track.tofSignal(), track.beta(),
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
