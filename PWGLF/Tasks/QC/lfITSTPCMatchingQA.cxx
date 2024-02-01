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

#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPr, aod::pidTPCKa, aod::pidTPCPi>;

namespace
{
// PID custom parametrisation for d and He3
constexpr double betheBlochDefault[2][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNamesBB{"d", "He3"};
static const std::vector<std::string> pidHypotheses{"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "He3", "Alpha", "Pion0", "Photon", "K0", "Lambda", "HyperTriton", "Hyperhydrog4", "XiMinus", "OmegaMinus"};
} // namespace

struct lfmatchingqa {

  // ConfigurableAxis for the histograms
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 2, 6, particleNamesBB, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for d and He3"};
  ConfigurableAxis tpcNsigmaAxis{"tpcNsigmaAxis", {100, -4.f, 4.f}, "tpc nsigma axis"};
  ConfigurableAxis momAxis{"momAxis", {60., -3.f, 3.f}, "momentum axis binning"};
  ConfigurableAxis momAxisFine{"momAxisFine", {2.e3, -5.f, 5.f}, "momentum axis binning"};
  ConfigurableAxis momResAxis{"momResAxis", {2.e2, -2.f, 2.f}, "momentum resolution binning"};
  ConfigurableAxis tpcAxis{"tpcAxis", {4e3, 0.f, 4.e3f}, "tpc signal axis binning"};
  ConfigurableAxis dcaAxis{"dcaAxis", {100, -0.1f, 0.1f}, "dca axis binning"};
  ConfigurableAxis itsClusSizeAxis{"itsClusSizeAxis", {90, 1, 15}, "its cluster size axis binning"};
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
    for (int iL{0}; iL < 6; ++iL) {
      sum += (track.itsClusterSizes() >> (iL * 4)) & 0xf;
    }
    return sum / track.itsNCls();
  }

  void init(o2::framework::InitContext&)
  {
    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    histos.add<TH2>("tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});TPC signal (a.u.)", {HistType::kTH2F}, {momAxisFine, tpcAxis});
    histos.add<TH3>("tpcSignalPIDHypo", ";#it{p}_{TPC} (GeV/#it{c});TPC signal (a.u.); PID hypothesis", {HistType::kTH3F}, {momAxisFine, tpcAxis, trackingPidAxis});
    histos.add<TH3>("tpcNsigmaPi", ";#it{p}_{TPC}; #it{p}_{GLO}; n#sigma_{TPC} (#pi)", {HistType::kTH3F, {momAxis, momAxis, tpcNsigmaAxis}});
    histos.add<TH3>("tpcNsigmaKa", ";#it{p}_{TPC}; #it{p}_{GLO}; n#sigma_{TPC} (K)", {HistType::kTH3F, {momAxis, momAxis, tpcNsigmaAxis}});
    histos.add<TH3>("tpcNsigmaPr", ";#it{p}_{TPC}; #it{p}_{GLO}; n#sigma_{TPC} (p)", {HistType::kTH3F, {momAxis, momAxis, tpcNsigmaAxis}});
    histos.add<TH3>("tpcNsigmaDe", ";#it{p}_{TPC}; #it{p}_{GLO}; n#sigma_{TPC} (d)", {HistType::kTH3F, {momAxis, momAxis, tpcNsigmaAxis}});
    histos.add<TH3>("tpcNsigmaHe", ";#it{p}_{TPC}; #it{p}_{GLO}; n#sigma_{TPC} (He3)", {HistType::kTH3F, {momAxis, momAxis, tpcNsigmaAxis}});

    auto pidHypoPi = histos.add<TH3>("pidHypoPi", ";#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC} (#pi);", {HistType::kTH3F}, {momAxisFine, tpcNsigmaAxis, trackingPidAxis});
    auto pidHypoKa = histos.add<TH3>("pidHypoKa", ";#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC} (K);", {HistType::kTH3F}, {momAxisFine, tpcNsigmaAxis, trackingPidAxis});
    auto pidHypoPr = histos.add<TH3>("pidHypoPr", ";#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC} (p);", {HistType::kTH3F}, {momAxisFine, tpcNsigmaAxis, trackingPidAxis});
    auto pidHypoDe = histos.add<TH3>("pidHypoDe", ";#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC} (d);", {HistType::kTH3F}, {momAxisFine, tpcNsigmaAxis, trackingPidAxis});
    auto pidHypoHe = histos.add<TH3>("pidHypoHe", ";#it{p}_{TPC} (GeV/#it{c}); n#sigma_{TPC} (He3);", {HistType::kTH3F}, {momAxisFine, tpcNsigmaAxis, trackingPidAxis});
    for (int i{1}; i < pidHypoPi->GetNbinsZ() + 1; ++i) {
      pidHypoPi->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
      pidHypoKa->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
      pidHypoPr->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
      pidHypoDe->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
      pidHypoHe->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
    }

    histos.add<TH3>("momCorrPi", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c})", {HistType::kTH3F, {momAxis, momAxis, momResAxis}});
    histos.add<TH3>("momCorrKa", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c})", {HistType::kTH3F, {momAxis, momAxis, momResAxis}});
    histos.add<TH3>("momCorrPr", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c})", {HistType::kTH3F, {momAxis, momAxis, momResAxis}});
    histos.add<TH3>("momCorrDe", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c})", {HistType::kTH3F, {momAxis, momAxis, momResAxis}});
    histos.add<TH3>("momCorrHe", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c})", {HistType::kTH3F, {momAxis, momAxis, momResAxis}});

    histos.add<TH3>("dcaPi", "; #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm)", {HistType::kTH3F, {momAxis, momResAxis, dcaAxis}});
    histos.add<TH3>("dcaKa", "; #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm)", {HistType::kTH3F, {momAxis, momResAxis, dcaAxis}});
    histos.add<TH3>("dcaPr", "; #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm)", {HistType::kTH3F, {momAxis, momResAxis, dcaAxis}});
    histos.add<TH3>("dcaDe", "; #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm)", {HistType::kTH3F, {momAxis, momResAxis, dcaAxis}});
    histos.add<TH3>("dcaHe", "; #it{p}_{GLO} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c}); DCA_{xy} (cm)", {HistType::kTH3F, {momAxis, momResAxis, dcaAxis}});

    histos.add<TH3>("itsClusSizePi", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); <ITS Cluster size> x cos(#lambda) (#pi)", {HistType::kTH3F, {momAxis, momAxis, itsClusSizeAxis}});
    histos.add<TH3>("itsClusSizeKa", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); <ITS Cluster size> x cos(#lambda) (K)", {HistType::kTH3F, {momAxis, momAxis, itsClusSizeAxis}});
    histos.add<TH3>("itsClusSizePr", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); <ITS Cluster size> x cos(#lambda) (p)", {HistType::kTH3F, {momAxis, momAxis, itsClusSizeAxis}});
    histos.add<TH3>("itsClusSizeDe", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); <ITS Cluster size> x cos(#lambda) (d)", {HistType::kTH3F, {momAxis, momAxis, itsClusSizeAxis}});
    histos.add<TH3>("itsClusSizeHe", ";#it{p}_{TPC} (GeV/#it{c}); #it{p}_{GLO} (GeV/#it{c}); <ITS Cluster size> x cos(#lambda) (He3)", {HistType::kTH3F, {momAxis, momAxis, itsClusSizeAxis}});
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TracksFull const& tracks, aod::BCs const&)
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

      auto sign = track.sign();
      auto signTPCMom = sign * track.tpcInnerParam();
      auto signGloMom = sign * track.p();

      // compute custom tpcNsigmaDeu and tpcNsigmaHe3
      double expBetheDeu{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / constants::physics::MassDeuteron), cfgBetheBlochParams->get("d", "p0"), cfgBetheBlochParams->get("d", "p1"), cfgBetheBlochParams->get("d", "p2"), cfgBetheBlochParams->get("d", "p3"), cfgBetheBlochParams->get("d", "p4"))};
      double expSigmaDeu{expBetheDeu * cfgBetheBlochParams->get("d", "resolution")};
      double expBetheHe3{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / constants::physics::MassHelium3), cfgBetheBlochParams->get("He3", "p0"), cfgBetheBlochParams->get("He3", "p1"), cfgBetheBlochParams->get("He3", "p2"), cfgBetheBlochParams->get("He3", "p3"), cfgBetheBlochParams->get("He3", "p4"))};
      double expSigmaHe3{expBetheHe3 * cfgBetheBlochParams->get("He3", "resolution")};
      auto tpcNSigmaDeu = static_cast<float>((track.tpcSignal() - expBetheDeu) / expSigmaDeu);
      auto tpcNSigmaHe3 = static_cast<float>((track.tpcSignal() - expBetheHe3) / expSigmaHe3);

      // filling the nsigma histograms
      histos.fill(HIST("tpcSignal"), signTPCMom, track.tpcSignal());
      histos.fill(HIST("tpcSignalPIDHypo"), signTPCMom, track.tpcSignal(), track.pidForTracking());
      if (abs(track.tpcNSigmaPi()) < 4) {
        histos.fill(HIST("tpcNsigmaPi"), signTPCMom, signGloMom, track.tpcNSigmaPi());
        histos.fill(HIST("pidHypoPi"), signTPCMom, track.tpcNSigmaPi(), track.pidForTracking());
      }
      if (abs(track.tpcNSigmaKa()) < 4) {
        histos.fill(HIST("tpcNsigmaKa"), signTPCMom, signGloMom, track.tpcNSigmaKa());
        histos.fill(HIST("pidHypoKa"), signTPCMom, track.tpcNSigmaKa(), track.pidForTracking());
      }
      if (abs(track.tpcNSigmaPr()) < 4) {
        histos.fill(HIST("tpcNsigmaPr"), signTPCMom, signGloMom, track.tpcNSigmaPr());
        histos.fill(HIST("pidHypoPr"), signTPCMom, track.tpcNSigmaPr(), track.pidForTracking());
      }
      if (abs(tpcNSigmaDeu) < 4) {
        histos.fill(HIST("tpcNsigmaDe"), signTPCMom, signGloMom, tpcNSigmaDeu);
        histos.fill(HIST("pidHypoDe"), signTPCMom, tpcNSigmaDeu, track.pidForTracking());
      }
      if (abs(tpcNSigmaHe3) < 4) {
        histos.fill(HIST("tpcNsigmaHe"), signTPCMom, signGloMom, tpcNSigmaHe3);
        histos.fill(HIST("pidHypoHe"), signTPCMom, tpcNSigmaHe3, track.pidForTracking());
      }

      // Filling the mom corr and cl sizes histograms (nSigma < 2 required)
      // calculating cos(L) of the track
      float cosL = 1 / std::sqrt(1.f + track.tgl() * track.tgl());

      if (abs(track.tpcNSigmaPi()) < 2) {
        histos.fill(HIST("momCorrPi"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p());
        histos.fill(HIST("dcaPi"), signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY());
        histos.fill(HIST("itsClusSizePi"), signTPCMom, signGloMom, getITSClSize(track) * cosL);
      }
      if (abs(track.tpcNSigmaKa()) < 2) {
        histos.fill(HIST("momCorrKa"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p());
        histos.fill(HIST("dcaKa"), signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY());
        histos.fill(HIST("itsClusSizeKa"), signTPCMom, signGloMom, getITSClSize(track) * cosL);
      }
      if (abs(track.tpcNSigmaPr()) < 2) {
        histos.fill(HIST("momCorrPr"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p());
        histos.fill(HIST("dcaPr"), signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY());
        histos.fill(HIST("itsClusSizePr"), signTPCMom, signGloMom, getITSClSize(track) * cosL);
      }
      if (abs(tpcNSigmaDeu) < 2) {
        histos.fill(HIST("momCorrDe"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p());
        histos.fill(HIST("dcaDe"), signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY());
        histos.fill(HIST("itsClusSizeDe"), signTPCMom, signGloMom, getITSClSize(track) * cosL);
      }
      if (abs(tpcNSigmaHe3) < 2) {
        histos.fill(HIST("momCorrHe"), signTPCMom, signGloMom, track.tpcInnerParam() - track.p());
        histos.fill(HIST("dcaHe"), signGloMom, track.tpcInnerParam() - track.p(), track.dcaXY());
        histos.fill(HIST("itsClusSizeHe"), signTPCMom, signGloMom, getITSClSize(track) * cosL);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lfmatchingqa>(cfgc)};
}
