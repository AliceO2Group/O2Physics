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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

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
static const std::vector<std::string> selectionNames{"selection"};
static const int defaultParameters[9][nParameters]{{0}, {0}, {1}, {1}, {1}, {0}, {0}, {0}, {0}};
static const float defaultPIDSelection[9][nParameters]{{-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}, {-1.f}};
static constexpr int Np = 9;
bool enableParticle[Np] = {false, false, false,
                           false, false, false,
                           false, false, false};
std::array<std::shared_ptr<TH2>, Np> hNsigmaPos;
std::array<std::shared_ptr<TH2>, Np> hNsigmaNeg;

template <typename TrackType>
float nsigmaITS(const TrackType& track, const o2::track::PID::ID id)
{
  switch (id) {
    case o2::track::PID::Electron:
      return track.itsNSigmaEl();
    case o2::track::PID::Muon:
      return track.itsNSigmaMu();
    case o2::track::PID::Pion:
      return track.itsNSigmaPi();
    case o2::track::PID::Kaon:
      return track.itsNSigmaKa();
    case o2::track::PID::Proton:
      return track.itsNSigmaPr();
    case o2::track::PID::Deuteron:
      return track.itsNSigmaDe();
    case o2::track::PID::Triton:
      return track.itsNSigmaTr();
    case o2::track::PID::Helium3:
      return track.itsNSigmaHe();
    case o2::track::PID::Alpha:
      return track.itsNSigmaAl();
    default:
      LOG(fatal) << "PID not implemented";
      return 0.f;
  }
}
template <typename TrackType>
float nsigmaTOF(const TrackType& track, const o2::track::PID::ID id)
{
  switch (id) {
    case o2::track::PID::Electron:
      return track.tofNSigmaEl();
    case o2::track::PID::Muon:
      return track.tofNSigmaMu();
    case o2::track::PID::Pion:
      return track.tofNSigmaPi();
    case o2::track::PID::Kaon:
      return track.tofNSigmaKa();
    case o2::track::PID::Proton:
      return track.tofNSigmaPr();
    case o2::track::PID::Deuteron:
      return track.tofNSigmaDe();
    case o2::track::PID::Triton:
      return track.tofNSigmaTr();
    case o2::track::PID::Helium3:
      return track.tofNSigmaHe();
    case o2::track::PID::Alpha:
      return track.tofNSigmaAl();
    default:
      LOG(fatal) << "PID not implemented";
      return 0.f;
  }
}
template <typename TrackType>
float nsigmaTPC(const TrackType& track, const o2::track::PID::ID id)
{
  switch (id) {
    case o2::track::PID::Electron:
      return track.tpcNSigmaEl();
    case o2::track::PID::Muon:
      return track.tpcNSigmaMu();
    case o2::track::PID::Pion:
      return track.tpcNSigmaPi();
    case o2::track::PID::Kaon:
      return track.tpcNSigmaKa();
    case o2::track::PID::Proton:
      return track.tpcNSigmaPr();
    case o2::track::PID::Deuteron:
      return track.tpcNSigmaDe();
    case o2::track::PID::Triton:
      return track.tpcNSigmaTr();
    case o2::track::PID::Helium3:
      return track.tpcNSigmaHe();
    case o2::track::PID::Alpha:
      return track.tpcNSigmaAl();
    default:
      LOG(fatal) << "PID not implemented";
      return 0.f;
  }
}

float tpcSelValues[9];
float tofSelValues[9];

/// Task to produce the ITS QA plots
struct itsPidQa {
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr const char* pN[Np] = {"El", "Mu", "Pi", "Ka", "Pr", "De", "Tr", "He", "Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<LabeledArray<int>> enabledParticle{"enabledParticle",
                                                  {defaultParameters[0], 9, nParameters, tableNames, parameterNames},
                                                  "Produce QA for this species: 0 - no, 1 - yes"};
  Configurable<LabeledArray<float>> tofSelection{"tofSelection",
                                                 {defaultPIDSelection[0], 9, nParameters, tableNames, selectionNames},
                                                 "Selection on the TOF nsigma"};
  Configurable<LabeledArray<float>> tpcSelection{"tpcSelection",
                                                 {defaultPIDSelection[0], 9, nParameters, tableNames, selectionNames},
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
  ConfigurableAxis avClsBins{"avClsBins", {200, 0, 20}, "Binning in average cluster size"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};
  Configurable<int16_t> minTPCNcls{"minTPCNcls", 0, "Minimum number or TPC Clusters for tracks"};
  ConfigurableAxis tpcNclsBins{"tpcNclsBins", {16, 0, 160}, "Binning in number of clusters in TPC"};

  template <typename TrackType>
  float averageClusterSizeTrk(const TrackType& track)
  {
    return o2::aod::ITSResponse::averageClusterSize(track.itsClusterSizes());
  }

  float averageClusterSizePerCoslInv(uint32_t itsClusterSizes, float eta) { return o2::aod::ITSResponse::averageClusterSize(itsClusterSizes) * std::cosh(eta); }

  template <typename TrackType>
  float averageClusterSizePerCoslInv(const TrackType& track)
  {
    return averageClusterSizePerCoslInv(track.itsClusterSizes(), track.eta());
  }

  void init(o2::framework::InitContext& context)
  {
    o2::aod::ITSResponse::setParameters(context);
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
    const AxisSpec avClsAxis{avClsBins, "<ITS Cls. Size>"};
    const AxisSpec avClsEffAxis{avClsBins, "<ITS Cls. Size> / cosh(#eta)"};

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
      tpcSelValues[id] = tpcSelection->get(tableNames[id].c_str(), "selection");
      if (tpcSelValues[id] <= 0.f) {
        tpcSelValues[id] = 999.f;
      }
      tofSelValues[id] = tofSelection->get(tableNames[id].c_str(), "selection");
      if (tofSelValues[id] <= 0.f) {
        tofSelValues[id] = 999.f;
      }
    }
    histos.add("event/eta", "", kTH1D, {etaAxis});
    histos.add("event/phi", "", kTH1D, {phiAxis});
    histos.add("event/etaphi", "", kTH2F, {etaAxis, phiAxis});
    histos.add("event/length", "", kTH1D, {lAxis});
    histos.add("event/pt", "", kTH1D, {ptAxis});
    histos.add("event/p", "", kTH1D, {pAxis});

    for (int id = 0; id < 9; id++) {
      const int f = enabledParticle->get(tableNames[id].c_str(), "enable");
      if (f != 1) {
        continue;
      }
      // NSigma
      const char* axisTitle = Form("N_{#sigma}^{ITS}(%s)", pT[id]);
      const AxisSpec nSigmaAxis{nSigmaBins, axisTitle};
      enableParticle[id] = true;
      hNsigmaPos[id] = histos.add<TH2>(Form("nsigmaPos/%s", pN[id]), axisTitle, kTH2F, {pAxis, nSigmaAxis});
      hNsigmaNeg[id] = histos.add<TH2>(Form("nsigmaNeg/%s", pN[id]), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    }
    histos.add("event/averageClusterSize", "", kTH2D, {pAxis, avClsAxis});
    histos.add("event/averageClusterSizePerCoslInv", "", kTH2D, {pAxis, avClsEffAxis});
    histos.add("event/SelectedAverageClusterSize", "", kTH2D, {pAxis, avClsAxis});
    histos.add("event/SelectedAverageClusterSizePerCoslInv", "", kTH2D, {pAxis, avClsEffAxis});
    LOG(info) << "QA PID ITS histograms:";
    histos.print();
  }

  Filter eventFilter = (o2::aod::evsel::sel8 == true && nabs(o2::aod::collision::posZ) < 10.f);
  // Filter trackFilter = (requireGlobalTrackInFilter());
  using CollisionCandidate = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                    aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                    aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                    aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl,
                                    aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                    aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe,
                                    aod::pidTOFTr, aod::pidTOFHe, aod::pidTOFAl>;
  void process(CollisionCandidate const& collision,
               TrackCandidates const& tracks)
  {
    auto tracksWithPid = soa::Attach<TrackCandidates,
                                     aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                     aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                     aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

    if (tracks.size() != tracksWithPid.size()) {
      LOG(fatal) << "Mismatch in track table size!" << tracks.size() << " vs " << tracksWithPid.size();
    }
    histos.fill(HIST("event/evsel"), 1);
    histos.fill(HIST("event/evsel"), 2);
    histos.fill(HIST("event/evsel"), 3);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& track : tracksWithPid) {
      histos.fill(HIST("event/trackselection"), 1.f);
      if (!track.isGlobalTrack()) { // Skipping non global tracks
        continue;
      }
      histos.fill(HIST("event/trackselection"), 2.f);
      if (!track.hasITS()) { // Skipping tracks without ITS
        continue;
      }
      histos.fill(HIST("event/trackselection"), 3.f);
      if (!track.hasTPC()) { // Skipping tracks without TPC
        continue;
      }
      histos.fill(HIST("event/trackselection"), 4.f);
      if (track.tpcNClsFound() < minTPCNcls) { // Skipping tracks without enough TPC clusters
        continue;
      }

      histos.fill(HIST("event/trackselection"), 5.f);
      histos.fill(HIST("event/particlehypo"), track.pidForTracking());
      histos.fill(HIST("event/eta"), track.eta());
      histos.fill(HIST("event/phi"), track.phi());
      histos.fill(HIST("event/etaphi"), track.eta(), track.phi());
      histos.fill(HIST("event/length"), track.length());
      histos.fill(HIST("event/pt"), track.pt());
      histos.fill(HIST("event/p"), track.p());
      histos.fill(HIST("event/averageClusterSize"), track.p(), averageClusterSizeTrk(track));
      histos.fill(HIST("event/averageClusterSizePerCoslInv"), track.p(), averageClusterSizePerCoslInv(track));
      bool discard = false;
      for (int id = 0; id < 9; id++) {
        if (std::abs(nsigmaTPC(track, id)) > tpcSelValues[id]) {
          LOG(debug) << "Discarding based on TPC hypothesis " << id << " " << std::abs(nsigmaTPC(track, id)) << ">" << tpcSelValues[id];
          discard = true;
          break;
        }
        if (track.hasTOF()) {
          if (std::abs(nsigmaTOF(track, id)) > tofSelValues[id]) {
            LOG(debug) << "Discarding based on TOF hypothesis " << id << " " << std::abs(nsigmaTOF(track, id)) << ">" << tofSelValues[id];
            discard = true;
            break;
          }
        }
      }
      if (discard) {
        continue;
      }
      histos.fill(HIST("event/SelectedAverageClusterSize"), track.p(), averageClusterSizeTrk(track));
      histos.fill(HIST("event/SelectedAverageClusterSizePerCoslInv"), track.p(), averageClusterSizePerCoslInv(track));

      for (o2::track::PID::ID id = 0; id <= o2::track::PID::Last; id++) {
        if (!enableParticle[id]) {
          continue;
        }
        if (applyRapidityCut) {
          if (std::abs(track.rapidity(PID::getMass(id))) > 0.5) {
            continue;
          }
        }
        const float nsigma = nsigmaITS(track, id);
        if (track.sign() > 0) {
          hNsigmaPos[id]->Fill(track.pt(), nsigma);
        } else {
          hNsigmaNeg[id]->Fill(track.pt(), nsigma);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<itsPidQa>(cfgc)}; }
