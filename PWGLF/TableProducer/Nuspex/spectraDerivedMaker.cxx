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
/// \file   spectraTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///
/// \brief Task for the analysis of the spectra with the TOF detector.
///        Depending on the configuration it can also run on tiny tables.
///

// O2 includes
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/spectraTOF.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct spectraDerivedMaker {
  Configurable<float> cfgNSigmaCut{"cfgNSigmaCut", 10.f, "Value of the Nsigma cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 2.f, "Downsampling factor for the events for derived data"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {100, 0, 100}, "Binning for multiplicity"};
  ConfigurableAxis binsPercentile{"binsPercentile", {100, 0, 100}, "Binning for percentiles"};
  Configurable<int> multiplicityEstimator{"multiplicityEstimator", 0, "Flag to use a multiplicity estimator: 0 no multiplicity, 1 MultFV0M, 2 MultFT0M, 3 MultFDDM, 4 MultTracklets, 5 MultTPC, 6 MultNTracksPV, 7 MultNTracksPVeta1, 8 CentralityFT0C, 9 CentralityFT0M, 10 CentralityFV0A"};
  // Custom track cuts for the cut variation study
  TrackSelection customTrackCuts;
  Configurable<int> itsPattern{"itsPattern", 1, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 60.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.7f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 7.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 1.f, "Additional cut on the maximum value of the DCA xy (multiplicative factor)"};
  Configurable<float> maxDcaZ{"maxDcaZ", 3.f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // Custom track cuts
    LOG(info) << "Using custom track cuts from values:";
    LOG(info) << "\trequireITS=" << requireITS.value;
    LOG(info) << "\trequireTPC=" << requireTPC.value;
    LOG(info) << "\trequireGoldenChi2=" << requireGoldenChi2.value;
    LOG(info) << "\tmaxChi2PerClusterTPC=" << maxChi2PerClusterTPC.value;
    LOG(info) << "\tminNCrossedRowsTPC=" << minNCrossedRowsTPC.value;
    LOG(info) << "\tminTPCNClsFound=" << minTPCNClsFound.value;
    LOG(info) << "\tmaxChi2PerClusterITS=" << maxChi2PerClusterITS.value;
    LOG(info) << "\tmaxDcaZ=" << maxDcaZ.value;

    customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
    LOG(info) << "Customizing track cuts:";
    customTrackCuts.SetRequireITSRefit(requireITS.value);
    customTrackCuts.SetRequireTPCRefit(requireTPC.value);
    customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
    customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
    customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
    customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
    customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
    customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
    customTrackCuts.SetMaxDcaXYPtDep([](float /*pt*/) { return 10.f; }); // No DCAxy cut will be used, this is done via the member function of the task
    customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
    customTrackCuts.print();
    // Histograms
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};

    histos.add("event/vertexz", "", HistType::kTH1D, {vtxZAxis});
    histos.add("event/sampledvertexz", "Sampled collisions", HistType::kTH1D, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1D, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "INEL>0 (fraction)");
    h->GetXaxis()->SetBinLabel(3, "INEL>1 (fraction)");
    h->GetXaxis()->SetBinLabel(4, "Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(5, "INEL>0 (fraction)");
    h->GetXaxis()->SetBinLabel(6, "INEL>1 (fraction)");
    h->GetXaxis()->SetBinLabel(7, "posZ passed");
    h->GetXaxis()->SetBinLabel(8, "INEL>0 (fraction)");
    h->GetXaxis()->SetBinLabel(9, "INEL>1 (fraction)");

    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1D, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, Form("|#eta| < %.2f", cfgCutEta.value));
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    h->GetXaxis()->SetBinLabel(4, "TOF passed (partial)");

    h = histos.add<TH1>("evtime_tof", "event time selections from pidEvTimeFlags", kTH1D, {{10, -0.5, 9.5}});
    h->GetXaxis()->SetBinLabel(1, "AnyEvTime");
    h->GetXaxis()->SetBinLabel(2, "EvTimeDefined");
    h->GetXaxis()->SetBinLabel(3, "EvTimeTOF");
    h->GetXaxis()->SetBinLabel(4, "EvTimeT0AC");
    h->GetXaxis()->SetBinLabel(5, "EvTimeTOFT0AC");
    h->GetXaxis()->SetBinLabel(6, "AnyEvTime (selected)");
    h->GetXaxis()->SetBinLabel(7, "EvTimeDefined (selected)");
    h->GetXaxis()->SetBinLabel(8, "EvTimeTOF (selected)");
    h->GetXaxis()->SetBinLabel(9, "EvTimeT0AC (selected)");
    h->GetXaxis()->SetBinLabel(10, "EvTimeTOFT0AC (selected)");

    histos.add("Centrality/FV0A", "FV0A", HistType::kTH1D, {{binsPercentile, "Centrality FV0A"}});
    histos.add("Centrality/FT0M", "FT0M", HistType::kTH1D, {{binsPercentile, "Centrality FT0M"}});
    histos.add("Centrality/FT0A", "FT0A", HistType::kTH1D, {{binsPercentile, "Centrality FT0A"}});
    histos.add("Centrality/FT0C", "FT0C", HistType::kTH1D, {{binsPercentile, "Centrality FT0C"}});
    histos.add("Centrality/FDDM", "FDDM", HistType::kTH1D, {{binsPercentile, "Centrality FDDM"}});
    histos.add("Centrality/NTPV", "NTPV", HistType::kTH1D, {{binsPercentile, "Centrality NTPV"}});

    histos.add("Mult/FV0M", "MultFV0M", HistType::kTH1D, {{binsMultiplicity, "MultFV0M"}});
    histos.add("Mult/FT0M", "MultFT0M", HistType::kTH1D, {{binsMultiplicity, "MultFT0M"}});
    histos.add("Mult/FDDM", "MultFDDM", HistType::kTH1D, {{binsMultiplicity, "MultFDDM"}});

    // histos.add("Mult/Tracklets", "MultTracklets", HistType::kTH1D, {{binsMultiplicity, "MultTracklets"}});
    histos.add("Mult/TPC", "MultTPC", HistType::kTH1D, {{binsMultiplicity, "MultTPC"}});
    histos.add("Mult/NTracksPV", "MultNTracksPV", HistType::kTH1D, {{binsMultiplicity, "MultNTracksPV"}});
    histos.add("Mult/NTracksPVeta1", "MultNTracksPVeta1", HistType::kTH1D, {{binsMultiplicity, "MultNTracksPVeta1"}});

    if (doprocessMC) {
      auto hh = histos.add<TH1>("MC/GenRecoCollisions", "Generated and Reconstructed MC Collisions", kTH1D, {{10, 0.5, 10.5}});
      hh->GetXaxis()->SetBinLabel(1, "Collisions generated");
      hh->GetXaxis()->SetBinLabel(2, "Collisions reconstructed");
      hh->GetXaxis()->SetBinLabel(3, "INEL>0");
      hh->GetXaxis()->SetBinLabel(4, "INEL>1");
      hh->GetXaxis()->SetBinLabel(5, "hasParticleInFT0C && hasParticleInFT0A");
    }
  }

  template <bool fillHistograms = false, bool fillMultiplicity = false, typename CollisionType, typename TrackType>
  bool isEventSelected(CollisionType const& collision, TrackType const& /*tracks*/)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 1.f);
    }
    if constexpr (fillHistograms) {
      if (collision.multNTracksPVeta1() >= 1) {
        histos.fill(HIST("evsel"), 2.f);
      }
      if (collision.multNTracksPVeta1() >= 2) {
        histos.fill(HIST("evsel"), 3.f);
      }
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 4.f);
      if (collision.multNTracksPVeta1() >= 1) {
        histos.fill(HIST("evsel"), 5.f);
      }
      if (collision.multNTracksPVeta1() >= 2) {
        histos.fill(HIST("evsel"), 6.f);
      }
    }
    if (std::abs(collision.posZ()) > cfgCutVertex) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 7.f);
      if (collision.multNTracksPVeta1() >= 1) {
        histos.fill(HIST("evsel"), 8.f);
      }
      if (collision.multNTracksPVeta1() >= 2) {
        histos.fill(HIST("evsel"), 9.f);
      }
      histos.fill(HIST("event/vertexz"), collision.posZ());

      if constexpr (fillMultiplicity) {
        histos.fill(HIST("Centrality/FV0A"), collision.centFV0A());
        histos.fill(HIST("Centrality/FT0M"), collision.centFT0M());
        histos.fill(HIST("Centrality/FT0A"), collision.centFT0A());
        histos.fill(HIST("Centrality/FT0C"), collision.centFT0C());
        // histos.fill(HIST("Centrality/FDDM"), collision.centFDDM());
        // histos.fill(HIST("Centrality/NTPV"), collision.centNTPV());

        histos.fill(HIST("Mult/FV0M"), collision.multZeqFV0A());
        histos.fill(HIST("Mult/FT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
        histos.fill(HIST("Mult/FDDM"), collision.multZeqFDDA() + collision.multZeqFDDC());

        // histos.fill(HIST("Mult/Tracklets"), collision.multTracklets());
        histos.fill(HIST("Mult/TPC"), collision.multTPC());
        histos.fill(HIST("Mult/NTracksPV"), collision.multZeqNTracksPV());
        histos.fill(HIST("Mult/NTracksPVeta1"), collision.multNTracksPVeta1());
      }
    }

    // Last thing, check the sampling
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/sampledvertexz"), collision.posZ());
    }
    return true;
  }

  template <typename TrackType>
  bool passesCutWoDCA(TrackType const& track) const
  {
    if (customTrackCuts.IsSelected(track)) {
      return true;
    }
    return false;
  }

  template <bool fillHistograms = false, typename TrackType>
  bool isTrackSelected(TrackType const& track)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 1);
    }
    if (std::abs(track.eta()) > cfgCutEta) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 2);
    }

    if (!passesCutWoDCA(track)) {
      return false;
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 3);
      if (track.hasTOF()) {
        histos.fill(HIST("tracksel"), 4);
      }
    }
    return true;
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::pidEvTimeFlags, aod::TrackSelection, aod::TOFSignal>;

  Produces<o2::aod::SpColls> tableColl;
  Produces<o2::aod::SpTracks> tableTrack;
  unsigned int randomSeed = 0;
  void processData(CollisionCandidate::iterator const& collision,
                   soa::Join<TrackCandidates,
                             aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                             aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr> const& tracks,
                   aod::BCs const&)
  {
    if (!isEventSelected<true, true>(collision, tracks)) {
      return;
    }

    tableColl(collision.numContrib(),
              collision.posX(),
              collision.posY(),
              collision.posZ(),
              collision.centFT0M(),
              collision.sel8(),
              collision.multNTracksPVeta1(),
              collision.bc().runNumber());

    tableTrack.reserve(tracks.size());
    for (const auto& trk : tracks) {
      if (!isTrackSelected<false>(trk)) {
        continue;
      }

      tableTrack(tableColl.lastIndex(),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningNSigma>(trk.tpcNSigmaPi()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningNSigma>(trk.tpcNSigmaKa()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningNSigma>(trk.tpcNSigmaPr()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningNSigma>(trk.tofNSigmaPi()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningNSigma>(trk.tofNSigmaKa()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningNSigma>(trk.tofNSigmaPr()),
                 trk.pt() * trk.sign(), trk.eta(), trk.phi(),
                 trk.length(),
                 trk.tpcSignal(),
                 trk.tpcChi2NCl(), trk.itsChi2NCl(), trk.tofChi2(),
                 trk.tpcNClsShared(),
                 trk.tpcNClsFindable(),
                 trk.tpcNClsFindableMinusFound(),
                 trk.tpcNClsFindableMinusCrossedRows(),
                 trk.isPVContributor(),
                 trk.itsClusterSizes(),
                 trk.hasTRD(),
                 trk.tofFlags(),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningDCA>(trk.dcaXY()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningDCA>(trk.dcaZ()),
                 trk.isGlobalTrack(),
                 trk.isGlobalTrackWoDCA());
    }
  }
  PROCESS_SWITCH(spectraDerivedMaker, processData, "Process data for derived dataset production", true);

  using CollisionCandidateMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFT0Ms>;

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;
  void processMC(soa::Join<aod::Tracks, aod::TracksExtra,
                           aod::TracksDCA, aod::McTrackLabels,
                           aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                           aod::TrackSelection> const& /*tracks*/,
                 aod::McParticles const& mcParticles,
                 aod::McCollisions const& mcCollisions,
                 CollisionCandidateMC const& collisions)
  {
    // Fill number of generated and reconstructed collisions for normalization
    histos.fill(HIST("MC/GenRecoCollisions"), 1.f, mcCollisions.size());
    histos.fill(HIST("MC/GenRecoCollisions"), 2.f, collisions.size());

    // Loop on generated particles
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) > cfgCutY) {
        continue;
      }
    }

    // Loop on reconstructed collisions
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache);
      for (const auto& mcParticle : particlesInCollision) {
        if (std::abs(mcParticle.y()) > cfgCutY) {
          continue;
        }
      }
    }

    // Loop on generated collisions
    for (const auto& mcCollision : mcCollisions) {
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      bool hasParticleInFT0C = false;
      bool hasParticleInFT0A = false;

      int nInelPart = 0;
      for (const auto& mcParticle : particlesInCollision) {
        if (mcParticle.isPhysicalPrimary()) {
          if (mcParticle.eta() >= -3.4f && mcParticle.eta() <= -2.3f) { // Acceptance of the FT0C
            hasParticleInFT0C = true;
          }
          if (mcParticle.eta() >= 3.8f && mcParticle.eta() <= 5.0f) { // Acceptance of the FT0A
            hasParticleInFT0A = true;
          }
          if (std::abs(mcParticle.eta()) < 1.f) {
            nInelPart++;
          }
        }

        if (std::abs(mcParticle.y()) > cfgCutY) {
          continue;
        }
      }
      if (nInelPart >= 1) {
        histos.fill(HIST("MC/GenRecoCollisions"), 3.f);
      }
      if (nInelPart >= 2) {
        histos.fill(HIST("MC/GenRecoCollisions"), 4.f);
      }
      if (hasParticleInFT0C && hasParticleInFT0A) {
        histos.fill(HIST("MC/GenRecoCollisions"), 5.f);
      }
    }
  }
  PROCESS_SWITCH(spectraDerivedMaker, processMC, "Process MC for derived dataset production", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<spectraDerivedMaker>(cfgc)}; }
