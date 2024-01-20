// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \file flatenicty-chrg.cxx
/// \author Gyula Bencedi (gyula.bencedi@cern.ch), Antonio Ortiz
/// (antonio.ortiz@cern.ch) \brief Task to produce inclusive charged particle pT
/// distributions as a function of charged-particle flattenicity \since 2023

#include <cmath>
#include <string>
#include <string_view>
#include <vector>

#include "EventFiltering/filterTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsFT0/Digit.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

float meanMultT0A = 0.f;
float meanMultT0C = 0.f;
float meanMultV0A = 0.f;

struct FlattenictyCharged {

  HistogramRegistry flatchrg{
    "flatchrg",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  TrackSelection mTrackSelector;

  Configurable<float> cutTrkEta{"cutTrkEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cutTrkPtMin{"cutTrkPtMin", 0.15f, "Minimum pT of tracks"};
  //   Configurable<uint32_t> cutTrkMult{"cutTrkMult", 200, "max measured
  //   multiplicity"};
  Configurable<float> cutVtxzMin{"cutVtxzMin", -10.f,
                                 "Minimum value for z-vertex"};
  Configurable<float> cutVtxzMax{"cutVtxzMax", 10.f,
                                 "Maximum value for z-vertex"};

  ConfigurableAxis multBins{"multBins", {1001, -0.5, 1000.5}, ""};

  Configurable<bool> sel8{"sel8", 1, "apply sel8 event selection"};
  Configurable<bool> selt0vtx{"selt0vtx", 0, "apply T0 vertext trigger"};
  Configurable<bool> selt0time{"selt0time", 0, "apply 1ns cut T0A and T0C"};

  // //   Configurable<float> avPyT0A{"avPyT0A", 8.16, "nch from pythia T0A"};
  Configurable<float> avPyT0C{"avPyT0C", 8.83, "nch from pythia T0C"};
  Configurable<float> avPyFV0{"avPyFV0", 21.44, "nch from pythia FV0"};

  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch",
                                "URL of the CCDB database"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  void init(InitContext&)
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
    mTrackSelector.SetRequireHitsInITSLayers(
      1, {0, 1}); // one hit in any SPD layer
    mTrackSelector.SetMaxChi2PerClusterITS(36.f);
    mTrackSelector.SetMaxDcaXY(1.f);
    mTrackSelector.SetMaxDcaZ(1.f);

    std::vector<double> ptBinEdges = {
      0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
      0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2,
      1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6,
      2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0,
      6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
      18.0, 20.0, 22.0, 24.0, 26.0, 30.0};
    const AxisSpec PtAxis{ptBinEdges, "#it{p}_{T} (GeV/#it{c})", "pt"};
    const AxisSpec ZAxis = {40, -20.0, 20.0};
    const AxisSpec PhiAxis = {600, 0, 2 * M_PI};
    const AxisSpec EtaAxisGlobal = {50, -5.0, 5.0};
    const AxisSpec EtaAxis = {20, 2.2, 5.1}; // FV0
    const AxisSpec CombEstAxis = {500, -0.5, 499.5};
    AxisSpec MultAxis = {multBins, "N_{trk}"};

    flatchrg.add("hMultFV0", "hMultFV0", HistType::kTH1F,
                 {{1000, -0.5, 99999.5, "FV0 amplitude"}});
    flatchrg.add("hMultFV01to4Ring", "hMultFV01to4Ring", HistType::kTH1F,
                 {{1000, -0.5, 99999.5, "FV0 amplitude (rings 1-4)"}});
    flatchrg.add("hMultFV05Ring", "hMultFV05Ring", HistType::kTH1F,
                 {{1000, -0.5, 99999.5, "FV0 amplitude (ring 5)"}});

    flatchrg.add("hMultFV0sel", "hMultFV0sel", HistType::kTH1F,
                 {{1000, -0.5, 99999.5, "FV0 amplitude"}});
    flatchrg.add("hMultFV01to4Ringsel", "hMultFV01to4Ringsel", HistType::kTH1F,
                 {{1000, -0.5, 99999.5, "FV0 amplitude (rings 1-4)"}});
    flatchrg.add("hMultFV05Ringsel", "hMultFV05Ringsel", HistType::kTH1F,
                 {{1000, -0.5, 99999.5, "FV0 amplitude (ring 5)"}});

    flatchrg.add("hT0C_time", "T0C_time", HistType::kTH1F,
                 {{160, -40., 40., "FT0C time"}});
    flatchrg.add("hT0A_time", "T0A_time", HistType::kTH1F,
                 {{160, -40., 40., "FT0C time"}});
    flatchrg.add("hT0C_time_sel", "T0C_time", HistType::kTH1F,
                 {{160, -40., 40., "FT0C time"}});
    flatchrg.add("hT0A_time_sel", "T0A_time", HistType::kTH1F,
                 {{160, -40., 40., "FT0C time"}});

    flatchrg.add("hAmpT0AVsCh", "", HistType::kTH2F,
                 {{24, -0.5, 23.5, "ch"},
                  {600, -0.5, +5999.5, "FT0A amplitude vs channel"}});
    flatchrg.add("hFT0A", "FT0A", HistType::kTH1F,
                 {{600, -0.5, 599.5, "FT0A amplitudes"}});

    flatchrg.add("hAmpT0CVsCh", "", HistType::kTH2F,
                 {{28, -0.5, 27.5, "ch"},
                  {600, -0.5, +5999.5, "FT0C amplitude vs channel"}});
    flatchrg.add("hFT0C", "FT0C", HistType::kTH1F,
                 {{600, -0.5, 599.5, "FT0C amplitudes"}});

    flatchrg.add("hMultFT0C", "hMultFT0C", HistType::kTH1F,
                 {{600, -0.5, 5999.5, "FT0C amplitude"}});
    flatchrg.add("hMultFT0Csel", "hMultFT0C", HistType::kTH1F,
                 {{600, -0.5, 5999.5, "FT0C amplitude"}});
    flatchrg.add("hMultFT0A", "hMultFT0A", HistType::kTH1F,
                 {{600, -0.5, 5999.5, "FT0A amplitude"}});
    flatchrg.add("hMultFT0Asel", "hMultFT0A", HistType::kTH1F,
                 {{600, -0.5, 5999.5, "FT0A amplitude"}});

    flatchrg.add("h1flatencityFV0", "", HistType::kTH1F,
                 {{102, -0.01, 1.01, "1-flatencityFV0"}});
    flatchrg.add("h1flatencityFT0A", "", HistType::kTH1F,
                 {{102, -0.01, 1.01, "1-flatencityFT0A"}});
    flatchrg.add("h1flatencityFT0C", "", HistType::kTH1F,
                 {{102, -0.01, 1.01, "1-flatencityFT0C"}});
    flatchrg.add("h1flatencityFV0FT0C", "", HistType::kTH1F,
                 {{102, -0.01, 1.01, "1-flatencityFV0FT0C"}});
    flatchrg.add("hFV0FT0C", "", HistType::kTH1F,
                 {{102, -0.5, 499.5, "FV0_FT0C"}});

    flatchrg.add("hPtVsFV0FT0C", " ; #it{p}_{T} (GeV/#it{c}); FV0_FT0C",
                 HistType::kTH2F,
                 {{500, -0.5, 499.5, "fv0ft0c"},
                  {100, -0.5, 99.5, "#it{p}_{T} (GeV/#it{c})"}});
    flatchrg.add("hPtVs1flatencityFV0", " ; #it{p}_{T} (GeV/#it{c}); FV0_FT0C",
                 HistType::kTH2F,
                 {{102, -0.01, 1.01, "1flatFV0"},
                  {100, -0.5, 99.5, "#it{p}_{T} (GeV/#it{c})"}});

    // event level histos
    flatchrg.add({"Events/selection",
                  ";status;events",
                  {HistType::kTH1F, {{4, 0.5, 4.5}}}});
    auto hstat = flatchrg.get<TH1>(HIST("Events/selection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected trigger");
    x->SetBinLabel(3, "Selected zvtx");
    x->SetBinLabel(4, "Selected INEL>0");

    // track level histos
    flatchrg.add(
      {"Tracks/VtxZ", " ; #it{z}_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}});
    flatchrg.add({"Tracks/EtaVtxZGlobal",
                  "; #eta; #it{z}_{vtx} (cm); tracks",
                  {HistType::kTH2F, {EtaAxisGlobal, ZAxis}}});
    flatchrg.add({"Tracks/EtaGlobal",
                  "; #eta; #it{z}_{vtx} (cm); tracks",
                  {HistType::kTH1F, {EtaAxisGlobal}}});
    flatchrg.add({"Tracks/PhiEtaGlobal",
                  "; #varphi; #eta; tracks",
                  {HistType::kTH2F, {PhiAxis, EtaAxisGlobal}}});
    flatchrg.add({"Tracks/PtEtaGlobal",
                  " ; #it{p}_{T} (GeV/#it{c}); #eta",
                  {HistType::kTH2F, {PtAxis, EtaAxisGlobal}}});
  }

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

  Filter trackFilter =
    (nabs(aod::track::eta) < cutTrkEta) && (aod::track::pt > cutTrkPtMin);
  using TrackTableData =
    soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                            aod::TrackSelection>>;
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;

  void process(CollisionTableData::iterator const& collision,
               TrackTableData const& tracks,
               soa::Join<aod::BCs, aod::Timestamps> const& bcs,
               /*aod::MFTTracks const& mfttracks,*/ aod::FT0s const& ft0s,
               aod::FV0As const& fv0s)
  {
    LOGF(debug, "<FlattenictyCharged> Collision %d", collision.globalIndex());

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

    meanMultT0C = 0.f;
    auto vMeanMultT0C = ccdb->getForTimeStamp<std::vector<double>>(
      "Users/e/ekryshen/meanT0C", bc.timestamp());
    meanMultT0C = (*vMeanMultT0C)[0];

    //   meanMultT0A = 0.f;
    //   auto vMeanMultT0A =
    //   ccdb->getForTimeStamp<std::vector<double>>("Users/e/ekryshen/meanT0A",
    //   bc.timestamp()); meanMultT0A = (*vMeanMultT0A)[0];

    meanMultV0A = 0.f;
    auto vMeanMultV0A = ccdb->getForTimeStamp<std::vector<double>>(
      "Users/e/ekryshen/meanV0A", bc.timestamp());
    meanMultV0A = (*vMeanMultV0A)[0];

    //   float fac_FT0A_ebe = 1.;
    float fac_FT0C_ebe = 1.;
    float fac_FV0_ebe = 1.;

    //   if (meanMultT0A > 0) {
    //     fac_FT0A_ebe = avPyT0A / meanMultT0A;
    //   }
    if (meanMultT0C > 0) {
      fac_FT0C_ebe = avPyT0C / meanMultT0C;
    }
    if (meanMultV0A > 0) {
      fac_FV0_ebe = avPyFV0 / meanMultV0A;
    }

    bool isAcceptedEvent = false;
    flatchrg.fill(HIST("Events/selection"), 1.);

    if (collision.sel8()) {
      isAcceptedEvent = true;
    }

    if (!isAcceptedEvent) {
      return;
    }

    flatchrg.fill(HIST("Events/selection"), 2.);

    auto vtxZ = collision.posZ();
    flatchrg.fill(HIST("Tracks/VtxZ"), vtxZ);

    bool isGoodEvent = false;
    if (vtxZ > cutVtxzMin && vtxZ < cutVtxzMax) {
      isGoodEvent = true;
    }

    if (!isGoodEvent) {
      return;
    }

    flatchrg.fill(HIST("Events/selection"), 3.);

    //__________  FT0   mult. esimator
    //
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
      flatchrg.fill(HIST("hT0C_time"), t0_c);
      flatchrg.fill(HIST("hT0A_time"), t0_a);

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
        float amplitude = ft0.amplitudeA()[i_a];
        sumAmpFT0A += amplitude;
      }
      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        sumAmpFT0C += amplitude;
      }
      flatchrg.fill(HIST("hMultFT0A"), sumAmpFT0A);
      flatchrg.fill(HIST("hMultFT0C"), sumAmpFT0C);
    }

    //__________  FLATTENICITY - FV0 mult. esimator
    //

    double flatenicity_fv0 = 9999.;
    float sumAmpFV0 = 0;
    float sumAmpFV01to4Ring = 0;
    float sumAmpFV05Ring = 0;
    int innerFV0 = 32;

    const int nCells = 48;
    float RhoLattice[nCells];
    for (Int_t iCh = 0; iCh < nCells; iCh++) {
      RhoLattice[iCh] = 0;
    }

    if (collision.has_foundFV0()) {
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

    flatchrg.fill(HIST("h1flatencityFV0"), 1. - flatenicity_fv0);

    flatchrg.fill(HIST("hMultFV0"), sumAmpFV0);
    flatchrg.fill(HIST("hMultFV01to4Ring"), sumAmpFV01to4Ring);
    flatchrg.fill(HIST("hMultFV05Ring"), sumAmpFV05Ring);

    if (selt0vtx && !isOkvtxtrig && !isOkFV0OrA) {
      return;
    }
    if (selt0time && !isOkTimeFT0) { // to reduce beam-gas background
      return;
    }
    if (sel8 && !collision.sel8()) {
      return;
    }

    flatchrg.fill(HIST("hMultFV0sel"), sumAmpFV0);
    flatchrg.fill(HIST("hMultFV01to4Ringsel"), sumAmpFV01to4Ring);
    flatchrg.fill(HIST("hMultFV05Ringsel"), sumAmpFV05Ring);

    //__________  FLATTENICITY - FT0 mult. esimator
    //

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

    float sumAmpFT0Asel = 0.f;
    float sumAmpFT0Csel = 0.f;
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      float t0_a = ft0.timeA();
      float t0_c = ft0.timeC();
      flatchrg.fill(HIST("hT0C_time_sel"), t0_c);
      flatchrg.fill(HIST("hT0A_time_sel"), t0_a);

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
        float amplitude = ft0.amplitudeA()[i_a];
        uint8_t channel = ft0.channelA()[i_a];
        int sector = getT0ASector(channel);
        if (sector >= 0 && sector < 24) {
          RhoLatticeT0A[sector] += amplitude;
          flatchrg.fill(HIST("hAmpT0AVsCh"), sector, amplitude);
        }
        sumAmpFT0Asel += amplitude;
        flatchrg.fill(HIST("hFT0A"), amplitude);
      }
      flatchrg.fill(HIST("hMultFT0Asel"), sumAmpFT0Asel);

      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        uint8_t channel = ft0.channelC()[i_c];
        int sector = getT0CSector(channel);
        if (sector >= 0 && sector < 28) {
          RhoLatticeT0C[sector] += amplitude;
          flatchrg.fill(HIST("hAmpT0CVsCh"), sector, amplitude);
        }
        sumAmpFT0Csel += amplitude;
        flatchrg.fill(HIST("hFT0C"), amplitude);
      }
      flatchrg.fill(HIST("hMultFT0Csel"), sumAmpFT0Csel);
    }
    float flatenicity_t0a = GetFlatenicity(RhoLatticeT0A, nCellsT0A);
    float flatenicity_t0c = GetFlatenicity(RhoLatticeT0C, nCellsT0C);

    flatchrg.fill(HIST("h1flatencityFT0A"), 1. - flatenicity_t0a);
    flatchrg.fill(HIST("h1flatencityFT0C"), 1. - flatenicity_t0c);

    //__________  FLATTENICITY - combined

    float combest = 0.;

    // option 1
    flatchrg.fill(HIST("h1flatencityFV0FT0C"),
                  (1.0 - (flatenicity_fv0 + flatenicity_t0c) / 2.0));

    // option 2
    const int nEta = 2; // FT0C + FV0
    float weightsEta[nEta] = {0.0490638, 0.00353962};
    float factebye[nEta] = {0., 0.};
    float deltaEeta[nEta] = {1.1, 2.9};
    float ampl[nEta] = {0, 0};

    if (collision.has_foundFV0() && collision.has_foundFT0()) {

      float all_weights = 0;
      ampl[0] = sumAmpFT0C;
      ampl[1] = sumAmpFV0;
      factebye[0] = fac_FT0C_ebe;
      factebye[1] = fac_FV0_ebe;

      if (sumAmpFT0C > 0 && sumAmpFV0 > 0) {
        for (int ie = 0; ie < nEta; ++ie) {
          if (meanMultV0A > 0 && meanMultT0C > 0) {
            combest += ampl[ie] * weightsEta[ie] / deltaEeta[ie];
          } else {
            combest += ampl[ie] * factebye[ie] / deltaEeta[ie];
          }
          all_weights += weightsEta[ie];
        }
        combest /= all_weights;
      }
    }

    flatchrg.fill(HIST("hFV0FT0C"), combest);

    for (auto& track : tracks) {
      // if (!track.isGlobalTrack()) {
      if (!mTrackSelector.IsSelected(track)) {
        continue;
      }

      if (!track.isPVContributor()) {
        continue;
      }

      flatchrg.fill(HIST("Tracks/EtaGlobal"), track.eta());
      flatchrg.fill(HIST("Tracks/EtaVtxZGlobal"), track.eta(), vtxZ);
      flatchrg.fill(HIST("Tracks/PtEtaGlobal"), track.pt(), track.eta());

      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      flatchrg.fill(HIST("Tracks/PhiEtaGlobal"), phi, track.eta());

      flatchrg.fill(HIST("hPtVsFV0FT0C"), combest, track.pt());
      flatchrg.fill(HIST("hPtVs1flatencityFV0"), 1. - flatenicity_fv0,
                    track.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//       adaptAnalysisTask<FlattenictyQA>(cfgc),
                      adaptAnalysisTask<FlattenictyCharged>(cfgc)};
}
