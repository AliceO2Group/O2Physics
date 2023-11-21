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

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsFT0/Digit.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> mmObjectsNames{"kHmTrk", "kHmFv0", "kHmFt0", "kHfFt0", "kHmFt0cFv0", "kHfFt0cFv0", "kHtPt"};
float meanMultT0A = 0.f;
float meanMultT0C = 0.f;
float meanMultV0A = 0.f;

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
  TrackSelection mTrackSelector;
  // event selection cuts
  Configurable<float> selHTrkMult{"selHTrkMult", 38., "global trk multiplicity threshold"};
  Configurable<float> selHMFv0{"selHMFv0", 33559.5, "FV0-amplitude threshold"};
  Configurable<float> sel1Mft0{"sel1Mft0", 220.0, "FT0 mult threshold"};
  Configurable<float> sel1Fft0{"sel1Fft0", 0.855, "1-flatenicity FT0 threshold"};
  Configurable<float> sel1Mft0cFv0{"sel1Mft0cfv0", 280.0, "FT0C+FV0 mult threshold"};
  Configurable<float> sel1Fft0cFv0{"sel1Fft0cfv0", 0.895, "1-flatenicity FT0C+FV0 threshold"};
  Configurable<float> selPtTrig{"selPtTrig", 5., "track pT leading threshold"};
  Configurable<bool> sel8{"sel8", 1, "apply sel8 event selection"};
  Configurable<bool> selt0time{"selt0time", 0, "apply 1ns cut T0A and T0C"};
  Configurable<bool> selt0vtx{"selt0vtx", 0, "apply T0 vertext trigger"};

  Configurable<float> avPyT0A{"avPyT0A", 8.16, "nch from pythia T0A"};
  Configurable<float> avPyT0C{"avPyT0C", 8.83, "nch from pythia T0C"};
  Configurable<float> avPyFV0{"avPyFV0", 21.44, "nch from pythia FV0"};

  Configurable<float> maxFV0FT0Cm{"maxFV0FT0Cm", 10000., "upper cut on FV0+FT0C mult"};
  Configurable<float> maxFT0m{"maxFT0m", 10000., "upper cut on FT0 mult"};

  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch",
                                "URL of the CCDB database"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  static constexpr std::string_view nhEst_before[7] = {
    "eGlobaltrack", "eFV0", "eFT0", "e1flatencityFT0", "eFT0CFV0", "e1flatencityFT0CFV0", "ePtl"};
  static constexpr std::string_view nhEst_after[7] = {
    "eGlobaltrack_selected", "eFV0_selected", "eFT0_selected", "e1flatencityFT0_selected", "eFT0CFV0_selected", "e1flatencityFT0CFV0_selected", "ePtl_selected"};
  static constexpr std::string_view npEst[7] = {
    "epGlobaltrack", "epFV0", "epFT0", "ep1flatencityFT0", "epFT0CFV0", "ep1flatencityFT0CFV0", "epPtl"};

  static constexpr std::string_view tEst[7] = {
    "GlobalTrk", "FV0", "FT0", "1-flatencityFT0", "FT0C_FV0", "1-flatencity_FT0C_FV0", "pT^{trig} (GeV/#it{c})"};

  int RunNumber = 0;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  void init(o2::framework::InitContext&)
  {

    ccdbApi.init(o2::base::NameConf::getCCDBServer());
    ccdb->setURL(url.value); //
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (!ccdbApi.isHostReachable()) {
      LOGF(fatal, "CCDB host %s is not reacheable, cannot go forward",
           url.value.data());
    }
    mTrackSelector.SetPtRange(0.15f, 1e10f);
    mTrackSelector.SetEtaRange(-0.8f, 0.8f);
    mTrackSelector.SetRequireITSRefit(true);
    mTrackSelector.SetRequireTPCRefit(true);
    mTrackSelector.SetRequireGoldenChi2(false);
    mTrackSelector.SetMinNClustersTPC(60);
    mTrackSelector.SetMinNCrossedRowsTPC(70);
    mTrackSelector.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    mTrackSelector.SetMaxChi2PerClusterTPC(4.f);
    mTrackSelector.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
    mTrackSelector.SetMaxChi2PerClusterITS(36.f);
    mTrackSelector.SetMaxDcaXY(1.f);
    mTrackSelector.SetMaxDcaZ(1.f);

    int nBinsEst[7] = {100, 1000, 700, 102, 700, 102, 150};
    float lowEdgeEst[7] = {-0.5, -0.5, -0.5, -0.01, -0.5, -0.01, .0};
    float upEdgeEst[7] = {99.5, 99999.5, 699.5, 1.01, 699.5, 1.01, 150.0};

    // QA event level
    multiplicity.add("fCollZpos", "Vtx_z", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    // QA FT0
    multiplicity.add("hAmpT0AVsCh", "", HistType::kTH2F, {{24, -0.5, 23.5, "ch"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    multiplicity.add("hAmpT0CVsCh", "", HistType::kTH2F, {{28, -0.5, 27.5, "ch"}, {600, -0.5, +5999.5, "FT0C amplitude"}});
    multiplicity.add("hFT0C", "FT0C", HistType::kTH1F, {{600, -0.5, 599.5, "FT0C amplitudes"}});
    multiplicity.add("hFT0A", "FT0A", HistType::kTH1F, {{600, -0.5, 599.5, "FT0A amplitudes"}});

    multiplicity.add("hMultFT0C", "hMultFT0C", HistType::kTH1F, {{600, -0.5, 5999.5, "FT0C amplitude"}});
    multiplicity.add("hMultFT0A", "hMultFT0A", HistType::kTH1F, {{600, -0.5, 5999.5, "FT0A amplitude"}});
    multiplicity.add("hMultFV0", "hMultFV0", HistType::kTH1F, {{1000, -0.5, 99999.5, "FV0 amplitude"}});
    multiplicity.add("hMultFV01to4Ring", "hMultFV01to4Ring", HistType::kTH1F, {{1000, -0.5, 99999.5, "FV0 amplitude (rings 1-4)"}});
    multiplicity.add("hMultFV05Ring", "hMultFV05Ring", HistType::kTH1F, {{1000, -0.5, 99999.5, "FV0 amplitude (ring 5)"}});

    multiplicity.add("hMultFV0sel", "hMultFV0sel", HistType::kTH1F, {{1000, -0.5, 99999.5, "FV0 amplitude"}});
    multiplicity.add("hMultFV01to4Ringsel", "hMultFV01to4Ringsel", HistType::kTH1F, {{1000, -0.5, 99999.5, "FV0 amplitude (rings 1-4)"}});
    multiplicity.add("hMultFV05Ringsel", "hMultFV05Ringsel", HistType::kTH1F, {{1000, -0.5, 99999.5, "FV0 amplitude (ring 5)"}});

    multiplicity.add("hT0C_time", "T0C_time", HistType::kTH1F, {{160, -40., 40., "FT0C time"}});
    multiplicity.add("hT0A_time", "T0A_time", HistType::kTH1F, {{160, -40., 40., "FT0C time"}});
    multiplicity.add("hT0Cafter_time", "T0C_time", HistType::kTH1F, {{160, -40., 40., "FT0C time"}});
    multiplicity.add("hT0Aafter_time", "T0A_time", HistType::kTH1F, {{160, -40., 40., "FT0C time"}});

    multiplicity.add("hAmpT0AvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    multiplicity.add("hAmpT0CvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

    // QA global tracks
    multiplicity.add("hdNdetaGlobal", "dNdeta", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiGlobal", "Phi", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});
    multiplicity.add("hDCAxyGlobal", "DCA_{xy}", HistType::kTH1F, {{120, -4.0, 4.0, " "}});

    // estimators
    for (int i_e = 0; i_e < 7; ++i_e) {
      multiplicity.add(
        npEst[i_e].data(), "", HistType::kTProfile, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
    }
    for (int i_e = 0; i_e < 7; ++i_e) {
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
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::BCsWithTimestamps const&, TrackCandidates const& tracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

    meanMultT0C = 0.f;
    auto vMeanMultT0C = ccdb->getForTimeStamp<std::vector<double>>("Users/e/ekryshen/meanT0C", bc.timestamp());
    meanMultT0C = (*vMeanMultT0C)[0];
    meanMultT0A = 0.f;
    auto vMeanMultT0A = ccdb->getForTimeStamp<std::vector<double>>("Users/e/ekryshen/meanT0A", bc.timestamp());
    meanMultT0A = (*vMeanMultT0A)[0];
    meanMultV0A = 0.f;
    auto vMeanMultV0A = ccdb->getForTimeStamp<std::vector<double>>("Users/e/ekryshen/meanV0A", bc.timestamp());
    meanMultV0A = (*vMeanMultV0A)[0];

    float fac_FT0A_ebe = 1.;
    float fac_FT0C_ebe = 1.;
    float fac_FV0_ebe = 1.;
    if (meanMultT0A > 0) {
      fac_FT0A_ebe = avPyT0A / meanMultT0A;
    }
    if (meanMultT0C > 0) {
      fac_FT0C_ebe = avPyT0C / meanMultT0C;
    }
    if (meanMultV0A > 0) {
      fac_FV0_ebe = avPyFV0 / meanMultV0A;
    }

    bool keepEvent[kNtriggersMM]{false};
    auto vtxZ = collision.posZ();
    multiplicity.fill(HIST("fProcessedEvents"), 0);
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());

    float sumAmpFT0A = 0.f;
    float sumAmpFT0C = 0.f;
    bool isOkTimeFT0 = false;
    bool isOkvtxtrig = false;
    bool isOkFV0OrA = false;
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      std::bitset<8> triggers = ft0.triggerMask();
      isOkvtxtrig = triggers[o2::fit::Triggers::bitVertex];
      float t0_a = ft0.timeA();
      float t0_c = ft0.timeC();
      if (abs(t0_a) < 1. && abs(t0_c) < 1.) {
        isOkTimeFT0 = true;
      }
      multiplicity.fill(HIST("hT0C_time"), t0_c);
      multiplicity.fill(HIST("hT0A_time"), t0_a);

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
        float amplitude = ft0.amplitudeA()[i_a];
        sumAmpFT0A += amplitude;
      }
      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        sumAmpFT0C += amplitude;
      }
      multiplicity.fill(HIST("hMultFT0A"), sumAmpFT0A);
      multiplicity.fill(HIST("hMultFT0C"), sumAmpFT0C);
    }

    double flatenicity_fv0 = 9999;
    float sumAmpFV0 = 0;
    float sumAmpFV01to4Ring = 0;
    float sumAmpFV05Ring = 0;
    int innerFV0 = 32;
    const int nCells = 48; // 48 sectors in FV0

    if (collision.has_foundFV0()) {
      float RhoLattice[nCells];
      for (Int_t iCh = 0; iCh < nCells; iCh++) {
        RhoLattice[iCh] = 0;
      }

      auto fv0 = collision.foundFV0();
      std::bitset<8> fV0Triggers = fv0.triggerMask();
      isOkFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
      // LOGP(info, "amplitude.size()={}", fv0.amplitude().size());
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {

        int channelv0 = fv0.channel()[ich];
        float ampl_ch = fv0.amplitude()[ich];
        sumAmpFV0 += ampl_ch;
        if (channelv0 < innerFV0) {
          RhoLattice[channelv0] = ampl_ch;
          sumAmpFV01to4Ring += ampl_ch;
        } else {
          RhoLattice[channelv0] = ampl_ch / 2.0; // two channels per bin
          sumAmpFV05Ring += ampl_ch;
        }
      }

      flatenicity_fv0 = GetFlatenicity(RhoLattice, nCells);
    }
    multiplicity.fill(HIST("hMultFV0"), sumAmpFV0);
    multiplicity.fill(HIST("hMultFV01to4Ring"), sumAmpFV01to4Ring);
    multiplicity.fill(HIST("hMultFV05Ring"), sumAmpFV05Ring);

    if (selt0vtx && !isOkvtxtrig && !isOkFV0OrA) {
      tags(false, false, false, false, false, false, false);
      return;
    }

    if (selt0time && !isOkTimeFT0) { // this cut is expected to reduce the beam-gas bckgnd
      tags(false, false, false, false, false, false, false);
      return;
    }
    if (sel8 && !collision.sel8()) {
      tags(false, false, false, false, false, false, false);
      return;
    }

    multiplicity.fill(HIST("hMultFV0sel"), sumAmpFV0);
    multiplicity.fill(HIST("hMultFV01to4Ringsel"), sumAmpFV01to4Ring);
    multiplicity.fill(HIST("hMultFV05Ringsel"), sumAmpFV05Ring);

    // global observables
    int multTrack = 0;
    float flPt = 0; // leading pT

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
      float t0_a = ft0.timeA();
      float t0_c = ft0.timeC();
      multiplicity.fill(HIST("hT0Cafter_time"), t0_c);
      multiplicity.fill(HIST("hT0Aafter_time"), t0_a);

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
        float amplitude = ft0.amplitudeA()[i_a];
        uint8_t channel = ft0.channelA()[i_a];
        int sector = getT0ASector(channel);
        if (sector >= 0 && sector < 24) {
          RhoLatticeT0A[sector] += amplitude;
          multiplicity.fill(HIST("hAmpT0AVsCh"), sector, amplitude);
        }
        // sumAmpFT0A += amplitude;
        multiplicity.fill(HIST("hFT0A"), amplitude);
      }
      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        // sumAmpFT0C += amplitude;
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
      if (!mTrackSelector.IsSelected(track)) {
        continue;
      }
      // Has this track contributed to the collision vertex fit
      if (!track.isPVContributor()) {
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
    float cut[7];
    for (int i_e = 0; i_e < 7; ++i_e) {
      estimator[i_e] = 0;
      cut[i_e] = 0;
    }
    cut[0] = selHTrkMult;
    cut[1] = selHMFv0;
    cut[2] = sel1Mft0;
    cut[3] = sel1Fft0;
    cut[4] = sel1Mft0cFv0;
    cut[5] = sel1Fft0cFv0;
    cut[6] = selPtTrig;
    // option 5
    const int nEta5 = 2; // FT0C + FT0A
    float weigthsEta5[nEta5] = {0.0490638, 0.010958415};
    if (sumAmpFT0C > 0 && sumAmpFT0A > 0) {
      isOK_estimator5 = true;
    }
    if (isOK_estimator5) {
      if (meanMultT0A > 0 && meanMultT0C > 0) {
        combined_estimator5 = sumAmpFT0C * fac_FT0C_ebe + sumAmpFT0A * fac_FT0A_ebe;
      } else {
        combined_estimator5 = sumAmpFT0C * weigthsEta5[0] + sumAmpFT0A * weigthsEta5[1];
      }
    }
    // option 6
    const int nEta6 = 2; //  FT0C + FV0
    float weigthsEta6[nEta6] = {0.0490638, 0.00353962};
    if (sumAmpFT0C > 0 && sumAmpFV0 > 0) {
      isOK_estimator6 = true;
    }
    if (isOK_estimator6) {
      if (meanMultV0A > 0 && meanMultT0C > 0) {
        combined_estimator6 = sumAmpFT0C * fac_FT0C_ebe + sumAmpFV0 * fac_FV0_ebe;
      } else {
        combined_estimator6 = sumAmpFT0C * weigthsEta6[0] + sumAmpFV0 * weigthsEta6[1];
      }
    }

    estimator[0] = multTrack;
    estimator[1] = sumAmpFV0;
    estimator[2] = combined_estimator5;
    float flatenicity_ft0 = (flatenicity_t0a + flatenicity_t0c) / 2.0;
    estimator[3] = 1.0 - flatenicity_ft0;
    estimator[4] = combined_estimator6;
    estimator[5] = 1.0 - (flatenicity_fv0 + flatenicity_t0c) / 2.0;
    estimator[6] = flPt;

    static_for<0, 6>([&](auto i) {
      constexpr int index = i.value;
      multiplicity.fill(HIST(npEst[index]), estimator[index], estimator[0]);
      multiplicity.fill(HIST(nhEst_before[index]), estimator[index]);
    });

    static_for<0, 6>([&](auto i) {
      constexpr int index = i.value;
      if (estimator[index] > cut[index]) {
        if (index == 2) {
          if (estimator[index] < maxFT0m) {
            multiplicity.fill(HIST(nhEst_after[index]), estimator[index]);
            keepEvent[index] = true;
          }
        } else if (index == 4) {
          if (estimator[index] < maxFV0FT0Cm) {
            multiplicity.fill(HIST(nhEst_after[index]), estimator[index]);
            keepEvent[index] = true;
          }
        } else {
          multiplicity.fill(HIST(nhEst_after[index]), estimator[index]);
          keepEvent[index] = true;
        }
      }
    });

    tags(keepEvent[kHighTrackMult], keepEvent[kHighFv0Mult], keepEvent[kHighFt0Mult], keepEvent[kHighFt0Flat], keepEvent[kHighFt0cFv0Mult], keepEvent[kHighFt0cFv0Flat], keepEvent[kLeadingPtTrack]);

    if (!keepEvent[kHighTrackMult] && !keepEvent[kHighFv0Mult] && !keepEvent[kHighFt0Mult] && !keepEvent[kHighFt0Flat] && !keepEvent[kHighFt0cFv0Mult] && !keepEvent[kHighFt0cFv0Flat] && !keepEvent[kLeadingPtTrack]) {
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
WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<multFilter>(cfg)};
}
