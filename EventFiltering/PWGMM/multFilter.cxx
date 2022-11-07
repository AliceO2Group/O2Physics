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
#include <cmath>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"

#include "EventFiltering/filterTables.h"

#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> mmObjectsNames{"kHmTrk", "kHmFv0", "kHfFv0", "kHmFt0", "kHfFt0", "kHmFt0cFv0", "kHfFt0cFv0", "kHtPt"};

struct multFilter {
  enum { kHighTrackMult = 0,
         kHighFv0Mult,
         kHighFv0Flat,
         kHighFt0Mult,
         kHighFt0Flat,
         kHighFt0cFv0Mult,
         kHighFt0cFv0Flat,
         kLeadingPtTrack,
         kNtriggersMM };
  // my track selection, discussed with Mesut and Matia
  TrackSelection myTrackSelection();
  // event selection cuts
  Configurable<float> selHTrkMult{"selHTrkMult", 45., "global trk multiplicity threshold"};
  Configurable<float> selHMFv0{"selHMFv0", 33559.5, "FV0-amplitude threshold"};
  Configurable<float> sel1Ffv0{"sel1Ffv0", 0.93, "1-flatenicity FV0  threshold"};
  Configurable<float> sel1Mft0{"sel1Mft0", 116.0, "FT0 mult threshold"};
  Configurable<float> sel1Fft0{"sel1Fft0", 0.85, "1-flatenicity FT0 threshold"};
  Configurable<float> sel1Mft0cFv0{"sel1Mft0cfv0", 202.0, "FT0C+FV0 mult threshold"};
  Configurable<float> sel1Fft0cFv0{"sel1Fft0cfv0", 0.892, "1-flatenicity FT0C+FV0 threshold"};
  Configurable<float> selPtTrig{"selPtTrig", 11., "track pT leading threshold"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 1.5f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  static constexpr std::string_view nhEst_before[8] = {
    "eGlobaltrack", "eFV0", "e1flatencityFV0", "eFT0", "e1flatencityFT0", "eFT0CFV0", "e1flatencityFT0CFV0", "ePtl"};
  static constexpr std::string_view nhEst_after[8] = {
    "eGlobaltrack_selected", "eFV0_selected", "e1flatencityFV0_selected", "eFT0_selected", "e1flatencityFT0_selected", "eFT0CFV0_selected", "e1flatencityFT0CFV0_selected", "ePtl_selected"};
  static constexpr std::string_view npEst[8] = {
    "epGlobaltrack", "epFV0", "ep1flatencityFV0", "epFT0", "ep1flatencityFT0", "epFT0CFV0", "ep1flatencityFT0CFV0", "epPtl"};

  static constexpr std::string_view tEst[8] = {
    "GlobalTrk", "FV0", "1-flatencity_FV0", "FT0", "1-flatencityFT0", "FT0C_FV0", "1-flatencity_FT0C_FV0", "pT^{trig} (GeV/#it{c})"};

  void init(o2::framework::InitContext&)
  {
    int nBinsEst[8] = {100, 500, 102, 500, 102, 500, 102, 150};
    float lowEdgeEst[8] = {-0.5, -0.5, -0.01, -0.5, -0.01, -0.5, -0.01, .0};
    float upEdgeEst[8] = {99.5, 49999.5, 1.01, 499.5, 1.01, 499.5, 1.01, 150.0};

    // QA event level
    multiplicity.add("fCollZpos", "Vtx_z", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    // QA FT0
    multiplicity.add("hAmpT0AVsCh", "", HistType::kTH2F, {{24, -0.5, 23.5, "ch"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    multiplicity.add("hAmpT0CVsCh", "", HistType::kTH2F, {{28, -0.5, 27.5, "ch"}, {600, -0.5, +5999.5, "FT0C amplitude"}});
    multiplicity.add("hFT0C", "FT0C", HistType::kTH1F, {{600, -0.5, 599.5, "FT0C amplitudes"}});
    multiplicity.add("hFT0A", "FT0A", HistType::kTH1F, {{600, -0.5, 599.5, "FT0A amplitudes"}});
    multiplicity.add("hAmpT0AvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    multiplicity.add("hAmpT0CvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

    // QA global tracks
    multiplicity.add("hdNdetaGlobal", "dNdeta", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiGlobal", "Phi", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});
    multiplicity.add("hDCAxyGlobal", "DCA_{xy}", HistType::kTH1F, {{120, -4.0, 4.0, " "}});

    // estimators
    for (int i_e = 0; i_e < 8; ++i_e) {
      multiplicity.add(
        npEst[i_e].data(), "", HistType::kTProfile, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
    }
    for (int i_e = 0; i_e < 8; ++i_e) {
      multiplicity.add(
        nhEst_before[i_e].data(), "", HistType::kTH1F, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
      multiplicity.add(
        nhEst_after[i_e].data(), "", HistType::kTH1F, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
    }

    std::array<std::string, 2> eventTitles = {"all", "rejected"};

    auto scalers{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalers->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalers->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
    // overlap with HtrackMult
    auto scalerso1{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents_overlap1", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalerso1->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalerso1->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
    auto scalerso2{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents_overlap2", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalerso2->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalerso2->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
    auto scalerso3{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents_overlap3", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalerso3->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalerso3->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
  }

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::pt > cfgTrkLowPtCut);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;

  int getT0ASector(int i_ch)
  {
    int i_sec_t0a = -1;
    for (int i_sec = 0; i_sec < 24; ++i_sec) {
      if (i_ch >= 4 * i_sec && i_ch <= 3 + 4 * i_sec) {
        i_sec_t0a = i_sec;
        break;
      }
    }
    return i_sec_t0a;
  }
  int getT0CSector(int i_ch)
  {
    int i_sec_t0c = -1;
    for (int i_sec = 0; i_sec < 28; ++i_sec) {
      if (i_ch >= 4 * i_sec && i_ch <= 3 + 4 * i_sec) {
        i_sec_t0c = i_sec;
        break;
      }
    }
    return i_sec_t0c;
  }
  float GetFlatenicity(float signals[], int entries)
  {
    float flat = 9999;
    float mRho = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      mRho += 1.0 * signals[iCell];
    }
    // average activity per cell
    mRho /= (1.0 * entries);
    if (mRho <= 0) {
      return 9999;
    }
    // get sigma
    float sRho_tmp = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      sRho_tmp += (1.0 * signals[iCell] - mRho) * (1.0 * signals[iCell] - mRho);
    }
    sRho_tmp /= (1.0 * entries * entries);
    float sRho = sqrt(sRho_tmp);
    if (mRho > 0) {
      flat = sRho / mRho;
    }
    return flat;
  }
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TrackCandidates const& tracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {

    bool keepEvent[kNtriggersMM]{false};
    auto vtxZ = collision.posZ();
    multiplicity.fill(HIST("fProcessedEvents"), 0);
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    // global observables
    int multTrack = 0;
    float flPt = 0; // leading pT
    double flatenicity_fv0 = 9999;

    float sumAmpFV0 = 0;
    // float sumAmpFV01to4Ch = 0;
    int innerFV0 = 32;
    const int nCells = 48; // 48 sectors in FV0

    if (collision.has_foundFV0()) {
      float RhoLattice[nCells];
      for (Int_t iCh = 0; iCh < nCells; iCh++) {
        RhoLattice[iCh] = 0;
      }

      auto fv0 = collision.foundFV0();
      // LOGP(info, "amplitude.size()={}", fv0.amplitude().size());
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {

        int channelv0 = fv0.channel()[ich];
        float ampl_ch = fv0.amplitude()[ich];
        sumAmpFV0 += ampl_ch;
        if (channelv0 < innerFV0) {
          RhoLattice[channelv0] = ampl_ch;
        } else {
          RhoLattice[channelv0] = ampl_ch / 2.0; // two channels per bin
        }
      }

      flatenicity_fv0 = GetFlatenicity(RhoLattice, nCells);
    }

    // FT0
    float sumAmpFT0A = 0.f;
    float sumAmpFT0C = 0.f;
    const int nCellsT0A = 24;
    float RhoLatticeT0A[nCellsT0A];
    for (int iCh = 0; iCh < nCellsT0A; iCh++) {
      RhoLatticeT0A[iCh] = 0.0;
    }
    const int nCellsT0C = 28;
    float RhoLatticeT0C[nCellsT0C];
    for (int iCh = 0; iCh < nCellsT0C; iCh++) {
      RhoLatticeT0C[iCh] = 0.0;
    }

    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
        float amplitude = ft0.amplitudeA()[i_a];
        uint8_t channel = ft0.channelA()[i_a];
        int sector = getT0ASector(channel);
        if (sector >= 0 && sector < 24) {
          RhoLatticeT0A[sector] += amplitude;
          multiplicity.fill(HIST("hAmpT0AVsCh"), sector, amplitude);
        }
        sumAmpFT0A += amplitude;
        multiplicity.fill(HIST("hFT0A"), amplitude);
      }
      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        sumAmpFT0C += amplitude;
        uint8_t channel = ft0.channelC()[i_c];
        int sector = getT0CSector(channel);
        if (sector >= 0 && sector < 28) {
          RhoLatticeT0C[sector] += amplitude;
          multiplicity.fill(HIST("hAmpT0CVsCh"), sector, amplitude);
        }
        multiplicity.fill(HIST("hFT0C"), amplitude);
      }
      multiplicity.fill(HIST("hAmpT0AvsVtx"), vtxZ, sumAmpFT0A);
      multiplicity.fill(HIST("hAmpT0CvsVtx"), vtxZ, sumAmpFT0C);
    }
    float flatenicity_t0a = GetFlatenicity(RhoLatticeT0A, nCellsT0A);
    float flatenicity_t0c = GetFlatenicity(RhoLatticeT0C, nCellsT0C);

    // Globaltracks

    for (auto& track : tracks) {
      if (!myTrackSelection().IsSelected(track)) {
        continue;
      }
      float eta_a = track.eta();
      float phi_a = track.phi();
      float pt_a = track.pt();
      multTrack++;
      multiplicity.fill(HIST("hdNdetaGlobal"), eta_a);
      multiplicity.fill(HIST("hPhiGlobal"), phi_a);
      multiplicity.fill(HIST("hDCAxyGlobal"), track.dcaXY());
      if (flPt < pt_a) {
        flPt = pt_a;
      }
    }
    bool isOK_estimator5 = false;
    bool isOK_estimator6 = false;

    float combined_estimator5 = 0;
    float combined_estimator6 = 0;
    float estimator[8];
    float cut[8];
    for (int i_e = 0; i_e < 8; ++i_e) {
      estimator[i_e] = 0;
      cut[i_e] = 0;
    }
    cut[0] = selHTrkMult;
    cut[1] = selHMFv0;
    cut[2] = sel1Ffv0;
    cut[3] = sel1Mft0;
    cut[4] = sel1Fft0;
    cut[5] = sel1Mft0cFv0;
    cut[6] = sel1Fft0cFv0;
    cut[7] = selPtTrig;
    // option 5
    const int nEta5 = 2; // FT0C + FT0A
    // float weigthsEta5[nEta5] = {0.0569502, 0.014552069};// values for pilot run, 900 GeV
    float weigthsEta5[nEta5] = {0.0490638, 0.010958415};
    float ampl5[nEta5] = {0, 0};
    ampl5[0] = sumAmpFT0C;
    ampl5[1] = sumAmpFT0A;
    if (sumAmpFT0C > 0 && sumAmpFT0A > 0) {
      isOK_estimator5 = true;
    }
    if (isOK_estimator5) {
      for (int i_5 = 0; i_5 < nEta5; ++i_5) {
        combined_estimator5 += ampl5[i_5] * weigthsEta5[i_5];
      }
    }
    // option 6
    const int nEta6 = 2; //  FT0C + FV0
    // float weigthsEta6[nEta6] = {0.0569502, 0.00535717};
    float weigthsEta6[nEta6] = {0.0490638, 0.00353962};
    float ampl6[nEta6] = {0, 0};
    ampl6[0] = sumAmpFT0C;
    ampl6[1] = sumAmpFV0;
    if (sumAmpFT0C > 0 && sumAmpFV0 > 0) {
      isOK_estimator6 = true;
    }
    if (isOK_estimator6) {
      for (int i_6 = 0; i_6 < nEta6; ++i_6) {
        combined_estimator6 += ampl6[i_6] * weigthsEta6[i_6];
      }
    }

    estimator[0] = multTrack;
    estimator[1] = sumAmpFV0;
    estimator[2] = 1.0 - flatenicity_fv0;
    estimator[3] = combined_estimator5;
    float flatenicity_ft0 = (flatenicity_t0a + flatenicity_t0c) / 2.0;
    estimator[4] = 1.0 - flatenicity_ft0;
    estimator[5] = combined_estimator6;
    estimator[6] = 1.0 - (flatenicity_fv0 + flatenicity_t0c) / 2.0;
    estimator[7] = flPt;

    static_for<0, 7>([&](auto i) {
      constexpr int index = i.value;
      multiplicity.fill(HIST(npEst[index]), estimator[index], estimator[0]);
      multiplicity.fill(HIST(nhEst_before[index]), estimator[index]);
    });

    static_for<0, 7>([&](auto i) {
      constexpr int index = i.value;
      if (estimator[index] > cut[index]) {
        multiplicity.fill(HIST(nhEst_after[index]), estimator[index]);
        keepEvent[index] = true;
      }
    });

    tags(keepEvent[kHighTrackMult], keepEvent[kHighFv0Mult], keepEvent[kHighFv0Flat], keepEvent[kHighFt0Mult], keepEvent[kHighFt0Flat], keepEvent[kHighFt0cFv0Mult], keepEvent[kHighFt0cFv0Flat], keepEvent[kLeadingPtTrack]);

    if (!keepEvent[kHighTrackMult] && !keepEvent[kHighFv0Mult] && !keepEvent[kHighFv0Flat] && !keepEvent[kHighFt0Mult] && !keepEvent[kHighFt0Flat] && !keepEvent[kHighFt0cFv0Mult] && !keepEvent[kHighFt0cFv0Flat] && !keepEvent[kLeadingPtTrack]) {
      multiplicity.fill(HIST("fProcessedEvents"), 1);
    } else {
      for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
        if (keepEvent[iTrigger]) {
          multiplicity.fill(HIST("fProcessedEvents"), iTrigger + 2);
        }
      }
      if (keepEvent[kHighTrackMult]) {
        for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
          if (keepEvent[iTrigger]) {
            multiplicity.fill(HIST("fProcessedEvents_overlap1"), iTrigger + 2);
          }
        }
      }
      if (keepEvent[kHighFt0cFv0Mult]) {
        for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
          if (keepEvent[iTrigger]) {
            multiplicity.fill(HIST("fProcessedEvents_overlap2"), iTrigger + 2);
          }
        }
      }
      if (keepEvent[kHighFt0Mult]) {
        for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
          if (keepEvent[iTrigger]) {
            multiplicity.fill(HIST("fProcessedEvents_overlap3"), iTrigger + 2);
          }
        }
      }
    }
  }
};
TrackSelection multFilter::myTrackSelection()
{
  TrackSelection selectedTracks;
  selectedTracks.SetPtRange(0.15f, 1e10f);
  selectedTracks.SetEtaRange(-0.8f, 0.8f);
  selectedTracks.SetRequireITSRefit(true);
  selectedTracks.SetRequireTPCRefit(true);
  selectedTracks.SetRequireGoldenChi2(false);
  selectedTracks.SetMinNClustersTPC(60);
  selectedTracks.SetMinNCrossedRowsTPC(70);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
  selectedTracks.SetMaxChi2PerClusterTPC(4.f);
  selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
  selectedTracks.SetMaxChi2PerClusterITS(36.f);
  selectedTracks.SetMaxDcaXY(1.f);
  selectedTracks.SetMaxDcaZ(1.f);
  return selectedTracks;
}

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<multFilter>(cfg)};
}
