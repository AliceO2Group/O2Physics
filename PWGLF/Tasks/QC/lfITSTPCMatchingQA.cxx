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

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "MathUtils/BetheBlochAleph.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TDatabasePDG.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPr, aod::pidTPCKa, aod::pidTPCPi>;
using TracksFullMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPr, aod::pidTPCKa, aod::pidTPCPi, aod::McTrackLabels>;
namespace
{
// PID custom parametrisation for d and He3
constexpr double betheBlochDefault[2][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNamesBB{"d", "He3"};
static const std::vector<std::string> pidHypotheses{"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "He3", "Alpha", "Photon", "K0", "Lambda", "HyperTriton", "Hyperhydrog4", "XiMinus", "OmegaMinus"};
static const std::vector<int> pdgCodes{11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040, 22, 311, 3122, 1010010030, 1010010040, 3312, 3334};
} // namespace

struct lfmatchingqa {

  // ConfigurableAxis for the histograms
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 2, 6, particleNamesBB, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for d and He3"};
  ConfigurableAxis tpcNsigmaAxis{"tpcNsigmaAxis", {100, -4.f, 4.f}, "tpc nsigma axis"};
  ConfigurableAxis momAxis{"momAxis", {60., -3.f, 3.f}, "momentum axis binning"};
  ConfigurableAxis momAxisFine{"momAxisFine", {2.e3, -5.f, 5.f}, "momentum axis binning"};
  ConfigurableAxis momResAxis{"momResAxis", {1.e2, -2.f, 2.f}, "momentum resolution binning"};
  ConfigurableAxis tpcAxis{"tpcAxis", {4e2, 0.f, 1.e3f}, "tpc signal axis binning"};
  ConfigurableAxis tpcAxisFine{"tpcAxisFine", {4e3, 0.f, 4.e3f}, "tpc signal axis binning"};

  ConfigurableAxis dcaAxis{"dcaAxis", {50, -0.1f, 0.1f}, "dca axis binning"};
  ConfigurableAxis itsClusSizeAxis{"itsClusSizeAxis", {120, 1, 15}, "its cluster size axis binning"};
  ConfigurableAxis trackingPidAxis{"trackingPidAxis", {static_cast<double>(pidHypotheses.size()), 0, static_cast<double>(pidHypotheses.size())}, "tracking pid hypothesis binning"};

  // Cut values
  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
  Configurable<float> ptMin{"ptMin", 0.05f, "minimum pT (GeV/c)"};
  Configurable<float> nClusITSCut{"nClusITSCut", 7, "Minimum number of ITS clusters"};
  Configurable<float> nClusTPCCut{"nClusTPCCut", 70, "Minimum number of TPC clusters"};
  Configurable<float> dcaCut{"dcaCut", 0.1f, "DCA to PV"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <class T>
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax || track.pt() < ptMin || std::abs(track.dcaXY()) > dcaCut) {
      return false;
    }
    if (track.itsNCls() < nClusITSCut ||
        track.tpcNClsFound() < nClusTPCCut ||
        track.tpcNClsCrossedRows() < 100 ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcChi2NCl() > 4.f ||
        track.itsChi2NCl() > 36.f) {
      return false;
    }
    return true;
  }

  template <class T>
  float getITSClSize(T const& track)
  {
    float sum{0.f};
    for (int iL{0}; iL < 7; ++iL) {
      sum += (track.itsClusterSizes() >> (iL * 4)) & 0xf;
    }
    return sum / track.itsNCls();
  }

  template <class T>
  void fillHistograms(T const& track, int pidMC, float genPMC)
  {

    auto sign = track.sign();
    auto signTPCMom = sign * track.tpcInnerParam();
    auto signGloMom = sign * track.p();
    // calculating cos(L) of the track
    float cosL = 1 / std::sqrt(1.f + track.tgl() * track.tgl());

    // compute custom tpcNsigmaDeu and tpcNsigmaHe3
    double expBetheDeu{common::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / constants::physics::MassDeuteron), cfgBetheBlochParams->get("d", "p0"), cfgBetheBlochParams->get("d", "p1"), cfgBetheBlochParams->get("d", "p2"), cfgBetheBlochParams->get("d", "p3"), cfgBetheBlochParams->get("d", "p4"))};
    double expSigmaDeu{expBetheDeu * cfgBetheBlochParams->get("d", "resolution")};
    double expBetheHe3{common::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / constants::physics::MassHelium3), cfgBetheBlochParams->get("He3", "p0"), cfgBetheBlochParams->get("He3", "p1"), cfgBetheBlochParams->get("He3", "p2"), cfgBetheBlochParams->get("He3", "p3"), cfgBetheBlochParams->get("He3", "p4"))};
    double expSigmaHe3{expBetheHe3 * cfgBetheBlochParams->get("He3", "resolution")};
    auto tpcNSigmaDeu = static_cast<float>((track.tpcSignal() - expBetheDeu) / expSigmaDeu);
    auto tpcNSigmaHe3 = static_cast<float>((track.tpcSignal() - expBetheHe3) / expSigmaHe3);

    // filling the nsigma histograms
    histos.fill(HIST("tpcSignal"), signTPCMom, track.tpcSignal());

    if (abs(track.tpcNSigmaPi()) < 4) {
      histos.fill(HIST("thnPi"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY(), getITSClSize(track) * cosL, track.tpcSignal(), track.tpcNSigmaPi(), track.pidForTracking(), pidMC, genPMC);
    }
    if (abs(track.tpcNSigmaKa()) < 4) {
      histos.fill(HIST("thnKa"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY(), getITSClSize(track) * cosL, track.tpcSignal(), track.tpcNSigmaKa(), track.pidForTracking(), pidMC, genPMC);
    }
    if (abs(track.tpcNSigmaPr()) < 4) {
      histos.fill(HIST("thnPr"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY(), getITSClSize(track) * cosL, track.tpcSignal(), track.tpcNSigmaPr(), track.pidForTracking(), pidMC, genPMC);
    }
    if (abs(tpcNSigmaDeu) < 4) {
      histos.fill(HIST("thnDe"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY(), getITSClSize(track) * cosL, track.tpcSignal(), tpcNSigmaDeu, track.pidForTracking(), pidMC, genPMC);
    }
    if (abs(tpcNSigmaHe3) < 4) {
      histos.fill(HIST("thnHe"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY(), getITSClSize(track) * cosL, track.tpcSignal(), tpcNSigmaHe3, track.pidForTracking(), pidMC, genPMC);
    }
  }

  void init(o2::framework::InitContext&)
  {
    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});TPC signal (a.u.)", HistType::kTH2F, {momAxisFine, tpcAxisFine});

    // use a THnSparse for all the information
    auto thnPi = histos.add("thnPi", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); #it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm); <ITS Cluster size> x cos(#lambda); TPC signal (a.u.); n#sigma_{TPC} (pi); PID hypothesis; MC PDG code; MC Gen P (GeV/c)",
                            {HistType::kTHnSparseF, {momAxis, momAxis, momResAxis, dcaAxis, itsClusSizeAxis, tpcAxis, tpcNsigmaAxis, trackingPidAxis, trackingPidAxis, momAxis}});
    auto thnKa = histos.add("thnKa", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); #it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm); <ITS Cluster size> x cos(#lambda); TPC signal (a.u.); n#sigma_{TPC} (Ka); PID hypothesis; MC PDG code; MC Gen P (GeV/c)",
                            {HistType::kTHnSparseF, {momAxis, momAxis, momResAxis, dcaAxis, itsClusSizeAxis, tpcAxis, tpcNsigmaAxis, trackingPidAxis, trackingPidAxis, momAxis}});
    auto thnPr = histos.add("thnPr", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); #it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm); <ITS Cluster size> x cos(#lambda); TPC signal (a.u.); n#sigma_{TPC} (Pr); PID hypothesis; MC PDG code; MC Gen P (GeV/c)",
                            {HistType::kTHnSparseF, {momAxis, momAxis, momResAxis, dcaAxis, itsClusSizeAxis, tpcAxis, tpcNsigmaAxis, trackingPidAxis, trackingPidAxis, momAxis}});
    auto thnDe = histos.add("thnDe", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); #it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm); <ITS Cluster size> x cos(#lambda); TPC signal (a.u.); n#sigma_{TPC} (De); PID hypothesis; MC PDG code; MC Gen P (GeV/c)",
                            {HistType::kTHnSparseF, {momAxis, momAxis, momResAxis, dcaAxis, itsClusSizeAxis, tpcAxis, tpcNsigmaAxis, trackingPidAxis, trackingPidAxis, momAxis}});
    auto thnHe = histos.add("thnHe", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); #it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm); <ITS Cluster size> x cos(#lambda); TPC signal (a.u.); n#sigma_{TPC} (He); PID hypothesis; MC PDG code; MC Gen P (GeV/c)",
                            {HistType::kTHnSparseF, {momAxis, momAxis, momResAxis, dcaAxis, itsClusSizeAxis, tpcAxis, tpcNsigmaAxis, trackingPidAxis, trackingPidAxis, momAxis}});
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TracksFull const& tracks, aod::BCs const&)
  {

    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    histos.fill(HIST("zVtx"), collision.posZ());

    for (const auto& track : tracks) {
      if (!selectTrack(track)) {
        continue;
      }
      fillHistograms(track, -1, -999);
    }
  }
  PROCESS_SWITCH(lfmatchingqa, processData, "Data analysis", true);

  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TracksFullMC const& tracks, aod::McParticles const&, aod::BCs const&)
  {

    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    for (const auto& track : tracks) {
      if (!selectTrack(track)) {
        continue;
      }

      int pidMC = -1;
      float genPMC = -999;

      if (!track.has_mcParticle()) {
        fillHistograms(track, pidMC, genPMC);
        continue;
      }

      auto mcParticle = track.mcParticle();
      auto pdg = mcParticle.pdgCode();
      bool isPdgFound = false;
      for (size_t iPid = 0; iPid < pdgCodes.size(); ++iPid) {
        if (abs(pdg) == pdgCodes[iPid]) {
          pidMC = iPid;
          int sign = pdg > 0 ? 1 : -1;
          genPMC = mcParticle.p() * sign;
          isPdgFound = true;
          break;
        }
      }
      isPdgFound ? fillHistograms(track, pidMC, genPMC) : fillHistograms(track, -1, -999);
    }
  }
  PROCESS_SWITCH(lfmatchingqa, processMC, "MC analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lfmatchingqa>(cfgc)};
}
