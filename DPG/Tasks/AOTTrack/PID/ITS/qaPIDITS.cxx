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
/// \file   qaPIDITS.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Implementation for QA tasks of the ITS PID quantities
///

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

static constexpr int nParameters = 1;
static const std::vector<std::string> tableNames{"Electron", // 0
                                                 "Muon",     // 1
                                                 "Pion",     // 2
                                                 "Kaon",     // 3
                                                 "Proton",   // 4
                                                 "Deuteron", // 5
                                                 "Triton",   // 6
                                                 "Helium",   // 7
                                                 "Alpha"};   // 8
static const std::vector<std::string> parameterNames{"enable"};
static const int defaultParameters[9][nParameters]{{0}, {0}, {1}, {1}, {1}, {0}, {0}, {0}, {0}};
static const float defaultPIDSelection[9][nParameters]{{-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}};
static constexpr int Np = 9;
std::array<std::shared_ptr<TH2>, Np> hNsigmaPos;
std::array<std::shared_ptr<TH2>, Np> hNsigmaNeg;

/// Task to produce the ITS QA plots
struct itsPidQa {
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr const char* pN[Np] = {"El", "Mu", "Pi", "Ka", "Pr", "De", "Tr", "He", "Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<LabeledArray<int>> enabledTables{"enabledTables",
                                                {defaultParameters[0], 9, nParameters, tableNames, parameterNames},
                                                "Produce QA for this species: 0 - no, 1 - yes"};
  Configurable<LabeledArray<int>> tofSelection{"tofSelection",
                                               {defaultPIDSelection[0], 9, nParameters, tableNames, parameterNames},
                                               "Selection on the TOF nsigma"};
  Configurable<LabeledArray<int>> tpcSelection{"tpcSelection",
                                               {defaultPIDSelection[0], 9, nParameters, tableNames, parameterNames},
                                               "Selection on the TPC nsigma"};

  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  ConfigurableAxis etaBins{"etaBins", {100, -1.f, 1.f}, "Binning in eta"};
  ConfigurableAxis phiBins{"phiBins", {100, 0, TMath::TwoPi()}, "Binning in phi"};
  ConfigurableAxis trackLengthBins{"trackLengthBins", {100, 0, 1000.f}, "Binning in track length plot"};
  ConfigurableAxis deltaBins{"deltaBins", {200, -1000.f, 1000.f}, "Binning in Delta (dEdx - expected dEdx)"};
  ConfigurableAxis expSigmaBins{"expSigmaBins", {200, 0.f, 200.f}, "Binning in expected Sigma"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {401, -10.025f, 10.025f}, "Binning in NSigma"};
  ConfigurableAxis dEdxBins{"dEdxBins", {5000, 0.f, 5000.f}, "Binning in dE/dx"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply event selection cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};
  Configurable<bool> enableDeDxPlot{"enableDeDxPlot", true, "Enables the dEdx plot (reduces memory footprint if off)"};
  Configurable<int16_t> minTPCNcls{"minTPCNcls", 0, "Minimum number or TPC Clusters for tracks"};
  ConfigurableAxis tpcNclsBins{"tpcNclsBins", {16, 0, 160}, "Binning in number of clusters in TPC"};
  Configurable<bool> fillTHnSparses{"fillTHnSparses", false, "Flag to fill multidimensional histograms for nsigma vs pt, eta, Ncls"};

  float tpcSelValues[9];
  float tofSelValues[9];
  void init(o2::framework::InitContext&)
  {
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{etaBins, "#it{#eta}"};
    const AxisSpec phiAxis{phiBins, "#it{#phi}"};
    const AxisSpec lAxis{trackLengthBins, "Track length (cm)"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T}/|Z| (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p}/|Z| (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
    }
    const AxisSpec dedxAxis{dEdxBins, "d#it{E}/d#it{x} Arb. units"};
    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1D, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed vtx Z");

    h = histos.add<TH1>("event/trackselection", "", kTH1D, {{10, 0.5, 10.5, "Selection passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTPC");
    h->GetXaxis()->SetBinLabel(5, Form("tpcNClsFound > %i", minTPCNcls.value));

    histos.add("event/vertexz", "", kTH1D, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1D, {{10, 0, 10, "PID in tracking"}});
    for (int id = 0; id < 9; id++) {
      h->GetXaxis()->SetBinLabel(id + 1, PID::getName(id));
      tpcSelValue[id] = tpcSelection->get(tableNames[id].c_str(), "enable");
      if (tpcSelValue[id] <= 0.f) {
        tpcSelValue[id] = 999.f;
      }
      tofSelValue[id] = tofSelection->get(tableNames[id].c_str(), "enable");
      if (tofSelValue[id] <= 0.f) {
        tofSelValue[id] = 999.f;
      }
    }
    histos.add("event/eta", "", kTH1D, {etaAxis});
    histos.add("event/phi", "", kTH1D, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1D, {lAxis});
    histos.add("event/pt", "", kTH1D, {ptAxis});
    histos.add("event/p", "", kTH1D, {pAxis});

    for (int id = 0; id < 9; id++) {
      const int f = enabledTables->get(tableNames[id].c_str(), "enable");
      if (f != 1) {
        continue;
      }
      // NSigma
      const char* axisTitle = Form("N_{#sigma}^{ITS}(%s)", pT[id]);
      const AxisSpec nSigmaAxis{nSigmaBins, axisTitle};
      hNsigmaPos[id] = histos.add<TH2>(Form("nsigmaPos/%s", pN[id]), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      hNsigmaNeg[id] = histos.add<TH2>(Form("nsigmaNeg/%s", pN[id]), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    }
    LOG(info) << "QA PID ITS histograms:";
    histos.print();
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& /*tracks*/)
  {

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 1);
    }
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 2);
    }

    if (abs(collision.posZ()) > 10.f) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/evsel"), 3);
      histos.fill(HIST("event/vertexz"), collision.posZ());
    }
    return true;
  }

  template <bool fillHistograms, typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType& /*collision*/, const TrackType& track)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 1.f);
    }
    if (!track.isGlobalTrack()) { // Skipping non global tracks
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 2.f);
    }
    if (!track.hasITS()) { // Skipping tracks without ITS
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 3.f);
    }
    if (!track.hasTPC()) { // Skipping tracks without TPC
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 4.f);
    }
    if (track.tpcNClsFound() < minTPCNcls) { // Skipping tracks without enough TPC clusters
      return false;
    }

    if constexpr (fillHistograms) {
      histos.fill(HIST("event/trackselection"), 5.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/phi"), track.phi());
      histos.fill(HIST("event/etaphi"), track.eta(), track.phi());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
      // histos.fill(HIST("event/ptreso"), track.p(), track.sigma1Pt() * track.pt() * track.pt());
    }
    return true;
  }

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));
  Filter trackFilter = ((trackSelection.node() == 0) ||
                        ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                        ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                        ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                        ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                        ((trackSelection.node() == 5) && requireInAcceptanceTracksInFilter()));
  using CollisionCandidate = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                    aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                    aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                    aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl,
                                    aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                    aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe,
                                    aod::pidTOFTr, aod::pidTOFHe, aod::pidTOFAl>;

  template <typename TrackType>
  float nsigmaITS(const TrackType& track, const o2::track::PID::ID id)
  {
    switch (id) {
      case o2::track::PID::Electron:
        nsigma = track.itsNSigmaEl();
        break;
      case o2::track::PID::Muon:
        nsigma = track.itsNSigmaMu();
        break;
      case o2::track::PID::Pion:
        nsigma = track.itsNSigmaPi();
        break;
      case o2::track::PID::Kaon:
        nsigma = track.itsNSigmaKa();
        break;
      case o2::track::PID::Proton:
        nsigma = track.itsNSigmaPr();
        break;
      case o2::track::PID::Deuteron:
        nsigma = track.itsNSigmaDe();
        break;
      case o2::track::PID::Triton:
        nsigma = track.itsNSigmaTr();
        break;
      case o2::track::PID::Helium3:
        nsigma = track.itsNSigmaHe();
        break;
      case o2::track::PID::Alpha:
        nsigma = track.itsNSigmaAl();
        break;
      default:
        LOG(fatal) << "PID not implemented";
    }
  }
  template <typename TrackType>
  float nsigmaTOF(const TrackType& track, const o2::track::PID::ID id)
  {
    switch (id) {
      case o2::track::PID::Electron:
        nsigma = track.tofNSigmaEl();
        break;
      case o2::track::PID::Muon:
        nsigma = track.tofNSigmaMu();
        break;
      case o2::track::PID::Pion:
        nsigma = track.tofNSigmaPi();
        break;
      case o2::track::PID::Kaon:
        nsigma = track.tofNSigmaKa();
        break;
      case o2::track::PID::Proton:
        nsigma = track.tofNSigmaPr();
        break;
      case o2::track::PID::Deuteron:
        nsigma = track.tofNSigmaDe();
        break;
      case o2::track::PID::Triton:
        nsigma = track.tofNSigmaTr();
        break;
      case o2::track::PID::Helium3:
        nsigma = track.tofNSigmaHe();
        break;
      case o2::track::PID::Alpha:
        nsigma = track.tofNSigmaAl();
        break;
      default:
        LOG(fatal) << "PID not implemented";
    }
  }
  template <typename TrackType>
  float nsigmaTPC(const TrackType& track, const o2::track::PID::ID id)
  {
    switch (id) {
      case o2::track::PID::Electron:
        nsigma = track.tpcNSigmaEl();
        break;
      case o2::track::PID::Muon:
        nsigma = track.tpcNSigmaMu();
        break;
      case o2::track::PID::Pion:
        nsigma = track.tpcNSigmaPi();
        break;
      case o2::track::PID::Kaon:
        nsigma = track.tpcNSigmaKa();
        break;
      case o2::track::PID::Proton:
        nsigma = track.tpcNSigmaPr();
        break;
      case o2::track::PID::Deuteron:
        nsigma = track.tpcNSigmaDe();
        break;
      case o2::track::PID::Triton:
        nsigma = track.tpcNSigmaTr();
        break;
      case o2::track::PID::Helium3:
        nsigma = track.tpcNSigmaHe();
        break;
      case o2::track::PID::Alpha:
        nsigma = track.tpcNSigmaAl();
        break;
      default:
        LOG(fatal) << "PID not implemented";
    }
  }

  void process(CollisionCandidate const& collision,
               soa::Filtered<TrackCandidates> const& tracks)
  {
    auto tracksWithPid = soa::Attach<TrackCandidates,
                                     aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                     aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                     aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

    isEventSelected<true>(collision, tracks);
    for (const auto& track : tracksWithPid) {
      isTrackSelected<true>(collision, track);
      bool discard = false;
      for (int id = 0; id < 9; id++) {
        if (std::abs(nsigmaTPC(track, id)) > tpcSelValue[id]) {
          discard = true;
        }
        if (std::abs(nsigmaTOF(track, id)) > tofSelValue[id]) {
          discard = true;
        }
      }
      if (discard) {
        continue;
      }
      for (o2::track::PID::ID id = 0; id <= o2::track::PID::Last; id++) {
        if (applyRapidityCut) {
          if (std::abs(track.rapidity(PID::getMass(id))) > 0.5) {
            continue;
          }
        }
        const float nsigma = nsigmaITS(track, id);
        switch (id) {
          case o2::track::PID::Electron:
            nsigma = track.itsNSigmaEl();
            break;
          case o2::track::PID::Muon:
            nsigma = track.itsNSigmaMu();
            break;
          case o2::track::PID::Pion:
            nsigma = track.itsNSigmaPi();
            break;
          case o2::track::PID::Kaon:
            nsigma = track.itsNSigmaKa();
            break;
          case o2::track::PID::Proton:
            nsigma = track.itsNSigmaPr();
            break;
          case o2::track::PID::Deuteron:
            nsigma = track.itsNSigmaDe();
            break;
          case o2::track::PID::Triton:
            nsigma = track.itsNSigmaTr();
            break;
          case o2::track::PID::Helium3:
            nsigma = track.itsNSigmaHe();
            break;
          case o2::track::PID::Alpha:
            nsigma = track.itsNSigmaAl();
            break;
          default:
            LOG(fatal) << "PID not implemented";
        }
        if (track.sign() > 0) {
          hNsigmaPos[id]->Fill(track.p(), nsigma);
        } else {
          hNsigmaNeg[id]->Fill(track.p(), nsigma);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<itsPidQa>(cfgc)}; }
