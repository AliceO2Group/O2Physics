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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TRandom.h>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct flatenictyFV0 {

  Configurable<bool> isMC{"isMC", false, "option to flag mc"};
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<bool> applyCalibCh{"applyCalibCh", false, "equalize FV0"};
  Configurable<bool> applyCalibVtx{"applyCalibVtx", false, "equalize FV0 vs vtx"};
  Configurable<bool> applyNorm{"applyNorm", false, "normalization to eta"};
  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 1.5f,
                                   "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum  pT"};

  HistogramRegistry flatenicity{
    "flatenicity",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  static constexpr std::string_view nhEst[16] = {
    "eGlobaltrack", "eFDDAFDDCFT0CFV0MFT", "eFDDAFDDCFV0MFT", "eFV0MFT", "eFV0", "eMFTmult", "e1flatencityFV0", "e1flatencitytrkMFT", "e1flatencitytrkMFTFV0", "e1flatencityMFTFV0", "eMFTmultFT0A", "eFT0", "e1flatencityFT0", "e1flatencityMFTFT0A", "eFV0FT0C", "e1flatencityFV0FT0C"};
  static constexpr std::string_view tEst[16] = {
    "GlobalTrk", "FDDA_FDDC_FT0C_FV0_MFT", "FDDA_FDDC_FV0_MFT", "FV0_MFT", "FV0", "MFTTrk", "1-flatencity_FV0", "1-flatencity_trk_MFT", "1-flatencity_trk_MFT_FV0", "1-flatencity_MFT_FV0", "MFT_FT0A", "FT0", "1-flatencityFT0", "1-flatencity_MFT_FT0A", "FV0_FT0C", "1-flatencity_FV0_FT0C"};
  static constexpr std::string_view nhPtEst[16] = {
    "ptVsGlobaltrack", "ptVsFDDAFDDCFT0CFV0MFT", "ptVsFDDAFDDCFV0MFT", "ptVsFV0MFT", "ptVsFV0", "ptVsMFTmult", "ptVs1flatencityFV0", "ptVs1flatencitytrkMFT", "ptVs1flatencitytrkMFTFV0", "ptVs1flatencityMFTFV0", "ptVsMFTmultFT0A", "ptVsFT0", "ptVs1flatencityFT0", "ptVs1flatencityMFTFT0A", "ptVsFV0FT0C", "ptVs1flatencityFV0FT0C"};

  void init(o2::framework::InitContext&)
  {
    int nBinsEst[16] = {100, 600, 600, 600, 500, 200, 102, 102, 102, 102, 400, 200, 102, 102, 400, 102};
    float lowEdgeEst[16] = {-0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                            -0.01, -0.01, -0.01, -0.01, -0.5, -0.5, -0.01, -0.01, -0.5, -0.01};
    float upEdgeEst[16] = {99.5, 599.5, 599.5, 599.5, 39999.5, 199.5,
                           1.01, 1.01, 1.01, 1.01, 399.5, 199.5, 1.01, 1.01, 399.5, 1.01};

    ConfigurableAxis ptBinning{
      "ptBinning",
      {0, 0.15, 0.3, 0.6, 1.0, 2.0, 4.0, 6.0, 10.0, 20.0, 50.0},
      "pTassoc bin limits"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    flatenicity.add("hEv", "Ev", HistType::kTH1F,
                    {{6, -0.5, 5.5, "index activated detector"}});
    flatenicity.add("hFV0amplRing1to4", "FV01to4", HistType::kTH1F,
                    {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hFT0Aampl", "FTAampl", HistType::kTH1F,
                    {{30000, -0.5, +3999.5, "FT0A amplitude"}});
    flatenicity.add("hFT0Campl", "FTCampl", HistType::kTH1F,
                    {{10000, -0.5, +1999.5, "FT0C amplitude"}});
    flatenicity.add("hFDDAampl", "FDDAampl", HistType::kTH1F,
                    {{10000, -0.5, 39999.5, "FDDA amplitude"}});
    flatenicity.add("hFDDCampl", "FDDCampl", HistType::kTH1F,
                    {{10000, -0.5, 39999.5, "FDDC amplitude"}});
    flatenicity.add("hFT0C", "FT0C", HistType::kTH1F,
                    {{600, -0.5, 599.5, "FT0C amplitudes"}});
    flatenicity.add("hFT0A", "FT0A", HistType::kTH1F,
                    {{600, -0.5, 599.5, "FT0A amplitudes"}});

    flatenicity.add("hMFTmult", "", HistType::kTH1F,
                    {{200, -0.5, +199.5, "MFT mult (-3.6<#eta<-3.4)"}});
    flatenicity.add("hetaMFT", "", HistType::kTH1F,
                    {{100, -5.0, +5.0, "#eta"}});
    flatenicity.add("hphiMFT", "", HistType::kTH1F,
                    {{8, 0.0, 2.0 * M_PI, "#varphi"}});

    flatenicity.add("hMFTmultAll", "", HistType::kTH1F,
                    {{200, -0.5, +199.5, "MFT mult (-3.6<#eta<-2.5)"}});
    // estimators
    for (int i_e = 0; i_e < 16; ++i_e) {
      flatenicity.add(
        nhEst[i_e].data(), "", HistType::kTH2F,
        {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()},
         {100, -0.5, +99.5, "Global track"}});
    }

    // vs pT
    for (int i_e = 0; i_e < 16; ++i_e) {
      flatenicity.add(nhPtEst[i_e].data(), "", HistType::kTProfile,
                      {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
    }

    flatenicity.add("hvtxZ", "vtxZ", HistType::kTH1F, {{40, -20.0, 20.0, " "}});
    // QA flatenicity
    flatenicity.add("fMultFv0", "FV0 amp", HistType::kTH1F,
                    {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add(
      "hAmpV0VsCh", "", HistType::kTH2F,
      {{48, -0.5, 47.5, "channel"}, {500, -0.5, +19999.5, "FV0 amplitude"}});
    flatenicity.add(
      "hAmpV0VsChBeforeCalibration", "", HistType::kTH2F,
      {{48, -0.5, 47.5, "channel"}, {500, -0.5, +19999.5, "FV0 amplitude"}});

    flatenicity.add(
      "hAmpT0AVsChBeforeCalibration", "", HistType::kTH2F,
      {{24, -0.5, 23.5, "channel"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    flatenicity.add(
      "hAmpT0CVsChBeforeCalibration", "", HistType::kTH2F,
      {{28, -0.5, 27.5, "channel"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

    flatenicity.add(
      "hAmpT0AVsCh", "", HistType::kTH2F,
      {{24, -0.5, 23.5, "channel"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    flatenicity.add(
      "hAmpT0CVsCh", "", HistType::kTH2F,
      {{28, -0.5, 27.5, "channel"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

    flatenicity.add(
      "hAmpFDAVsChBeforeCalibration", "", HistType::kTH2F,
      {{8, -0.5, 7.5, "channel"}, {600, -0.5, +599.5, "FDA amplitude"}});
    flatenicity.add(
      "hAmpFDCVsChBeforeCalibration", "", HistType::kTH2F,
      {{8, -0.5, 7.5, "channel"}, {600, -0.5, +599.5, "FDC amplitude"}});
    flatenicity.add(
      "hAmpFDAVsCh", "", HistType::kTH2F,
      {{8, -0.5, 7.5, "channel"}, {600, -0.5, +599.5, "FDA amplitude"}});
    flatenicity.add(
      "hAmpFDCVsCh", "", HistType::kTH2F,
      {{8, -0.5, 7.5, "channel"}, {600, -0.5, +599.5, "FDC amplitude"}});

    flatenicity.add("hFlatMFTvsFlatGlob", "", HistType::kTH2F,
                    {{20, -0.01, +1.01, "flatenicity (Glob)"},
                     {20, -0.01, +1.01, "flatenicity (MFT)"}});
    flatenicity.add("hFlatMFTvsFlatFV0", "", HistType::kTH2F,
                    {{20, -0.01, +1.01, "flatenicity (FV0)"},
                     {20, -0.01, +1.01, "flatenicity (MFT)"}});
    flatenicity.add("hFlatFT0CvsFlatFT0A", "", HistType::kTH2F,
                    {{20, -0.01, +1.01, "flatenicity (FT0C)"},
                     {20, -0.01, +1.01, "flatenicity (FT0A)"}});
    flatenicity.add("fEtaPhiFv0", "eta vs phi", HistType::kTH2F,
                    {{8, 0.0, 2 * M_PI, "#phi (rad)"}, {5, 2.2, 5.1, "#eta"}});

    flatenicity.add("hAmpV0vsVtxBeforeCalibration", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Trk mult"},
                     {1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hAmpT0AvsVtxBeforeCalibration", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {600, -0.5, +5999.5, "FT0A amplitude"}});
    flatenicity.add("hAmpT0CvsVtxBeforeCalibration", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {600, -0.5, +5999.5, "FT0C amplitude"}});
    flatenicity.add("hMFTvsVtxBeforeCalibration", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {200, -0.5, +199.5, "MFT mult (-3.6<#eta<-2.5)"}});
    flatenicity.add("hAmpFDAvsVtxBeforeCalibration", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {6000, -0.5, 7999.5, "Ampl. FDA"}});
    flatenicity.add("hAmpFDCvsVtxBeforeCalibration", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {6000, -0.5, 7999.5, "Ampl. FDC"}});

    flatenicity.add("hAmpV0vsVtx", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Trk mult"},
                     {1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hAmpT0AvsVtx", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {600, -0.5, +5999.5, "FT0A amplitude"}});
    flatenicity.add("hAmpT0CvsVtx", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {600, -0.5, +5999.5, "FT0C amplitude"}});
    flatenicity.add("hMFTvsVtx", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {200, -0.5, +199.5, "MFT mult (-3.6<#eta<-2.5)"}});
    flatenicity.add("hAmpFDAvsVtx", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {6000, -0.5, 7999.5, "Ampl. FDA"}});
    flatenicity.add("hAmpFDCvsVtx", "", HistType::kTH2F,
                    {{30, -15.0, +15.0, "Vtx_z"},
                     {6000, -0.5, 7999.5, "Ampl. FDC"}});
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

  int getFV0Ring(int i_ch)
  {
    int i_ring = -1;
    if (i_ch < 8) {
      i_ring = 0;
    } else if (i_ch >= 8 && i_ch < 16) {
      i_ring = 1;
    } else if (i_ch >= 16 && i_ch < 24) {
      i_ring = 2;
    } else if (i_ch >= 24 && i_ch < 32) {
      i_ring = 3;
    } else {
      i_ring = 4;
    }
    return i_ring;
  }
  int getFV0IndexPhi(int i_ch)
  {
    int i_ring = -1;

    if (i_ch >= 0 && i_ch < 8) {
      if (i_ch < 4) {
        i_ring = i_ch;
      } else {
        if (i_ch == 7) {
          i_ring = 4;
        } else if (i_ch == 6) {
          i_ring = 5;
        } else if (i_ch == 5) {
          i_ring = 6;
        } else if (i_ch == 4) {
          i_ring = 7;
        }
      }
    } else if (i_ch >= 8 && i_ch < 16) {
      if (i_ch < 12) {
        i_ring = i_ch;
      } else {
        if (i_ch == 15) {
          i_ring = 12;
        } else if (i_ch == 14) {
          i_ring = 13;
        } else if (i_ch == 13) {
          i_ring = 14;
        } else if (i_ch == 12) {
          i_ring = 15;
        }
      }
    } else if (i_ch >= 16 && i_ch < 24) {
      if (i_ch < 20) {
        i_ring = i_ch;
      } else {
        if (i_ch == 23) {
          i_ring = 20;
        } else if (i_ch == 22) {
          i_ring = 21;
        } else if (i_ch == 21) {
          i_ring = 22;
        } else if (i_ch == 20) {
          i_ring = 23;
        }
      }
    } else if (i_ch >= 24 && i_ch < 32) {
      if (i_ch < 28) {
        i_ring = i_ch;
      } else {
        if (i_ch == 31) {
          i_ring = 28;
        } else if (i_ch == 30) {
          i_ring = 29;
        } else if (i_ch == 29) {
          i_ring = 30;
        } else if (i_ch == 28) {
          i_ring = 31;
        }
      }
    } else if (i_ch == 32) {
      i_ring = 32;
    } else if (i_ch == 40) {
      i_ring = 33;
    } else if (i_ch == 33) {
      i_ring = 34;
    } else if (i_ch == 41) {
      i_ring = 35;
    } else if (i_ch == 34) {
      i_ring = 36;
    } else if (i_ch == 42) {
      i_ring = 37;
    } else if (i_ch == 35) {
      i_ring = 38;
    } else if (i_ch == 43) {
      i_ring = 39;
    } else if (i_ch == 47) {
      i_ring = 40;
    } else if (i_ch == 39) {
      i_ring = 41;
    } else if (i_ch == 46) {
      i_ring = 42;
    } else if (i_ch == 38) {
      i_ring = 43;
    } else if (i_ch == 45) {
      i_ring = 44;
    } else if (i_ch == 37) {
      i_ring = 45;
    } else if (i_ch == 44) {
      i_ring = 46;
    } else if (i_ch == 36) {
      i_ring = 47;
    }
    return i_ring;
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
    // get sigma
    float sRho_tmp = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      sRho_tmp += TMath::Power(1.0 * signals[iCell] - mRho, 2);
    }
    sRho_tmp /= (1.0 * entries * entries);
    float sRho = TMath::Sqrt(sRho_tmp);
    if (mRho > 0) {
      flat = sRho / mRho;
    }
    return flat;
  }

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) &&
                       (aod::track::pt > cfgTrkLowPtCut);
  using TrackCandidates =
    soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra,
                            aod::TracksDCA, aod::TrackSelection>>;
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;

  void process(CollisionTableData::iterator const& collision,
               TrackCandidates const& tracks,
               soa::Join<aod::BCs, aod::Timestamps> const& bcs,
               aod::MFTTracks const& mfttracks, aod::FT0s const& ft0s,
               aod::FV0As const& fv0s, aod::FDDs_001 const& fdds_001)
  {

    auto vtxZ = collision.posZ();
    flatenicity.fill(HIST("hEv"), 0);
    flatenicity.fill(HIST("hvtxZ"), vtxZ);

    flatenicity.fill(HIST("hEv"), 1);
    bool isGoodEvent = false;
    if (vtxZ > -15.f && vtxZ < 15.0f) {
      flatenicity.fill(HIST("hEv"), 2);
      isGoodEvent = true;
    }

    if (!isGoodEvent) {
      return;
    }

    flatenicity.fill(HIST("hEv"), 3);
    // these values are from equalized pass 4 signals
    const int nEta1 = 5; // FDDC + MFTparc + FT0C + FV0 (rings 1-4) + FDDA
    // float weigthsEta1[nEta1] = {0.0117997, 1.66515, 0.0569502, 0.00548221, 0.0037175}; // after calibration
    float weigthsEta1[nEta1] = {0.000710054, 1.94347, 0.04924, 0.00451969, 0.00215551};
    float deltaEeta1[nEta1] = {2.0, 0.2, 1.1, 2.32, 1.6};
    float ampl1[nEta1] = {0, 0, 0, 0, 0};

    const int nEta2 = 4; // FDDC + MFT + FV0 (rings 1-4) + FDDA
    // float weigthsEta2[nEta2] = {0.0117997, 1.05258, 0.00569502, 0.0037175}; // after calibration
    float weigthsEta2[nEta2] = {0.000710054, 1.27606, 0.00451969, 0.00215551};
    float deltaEeta2[nEta2] = {2.0, 1.1, 2.32, 1.6};
    float ampl2[nEta2] = {0, 0, 0, 0};

    const int nEta3 = 2; // MFT + FV0
    // float weigthsEta3[nEta3] = {1.05258, 0.00535717}; // after calibration
    float weigthsEta3[nEta3] = {1.27606, 0.00437892};
    float deltaEeta3[nEta3] = {1.1, 2.9};
    float ampl3[nEta3] = {0, 0};

    const int nEta4 = 2; // MFT + FT0A
    // float weigthsEta4[nEta4] = {1.05258, 0.014552069}; // after calibration
    float weigthsEta4[nEta4] = {0.00437892, 0.0097735442};
    float deltaEeta4[nEta4] = {1.1, 1.2};
    float ampl4[nEta4] = {0, 0};

    const int nEta5 = 2; // FT0C + FT0A
    // float weigthsEta5[nEta5] = {0.0569502, 0.014552069}; // after calibration
    float weigthsEta5[nEta5] = {0.04924, 0.0097735442};
    float deltaEeta5[nEta5] = {1.1, 1.2};
    float ampl5[nEta5] = {0, 0};

    const int nEta6 = 2; // FT0C + FV0
    // float weigthsEta6[nEta6] = {0.0569502, 0.00535717}; // after calibration
    float weigthsEta6[nEta6] = {0.04924, 0.00437892};
    float deltaEeta6[nEta6] = {1.1, 2.9};
    float ampl6[nEta6] = {0, 0};

    // V0A signal and flatenicity calculation
    float flatenicity_fv0;
    float calib[48] = {1.01697, 1.122, 1.03854, 1.108, 1.11634, 1.14971, 1.19321, 1.06866, 0.954675, 0.952695, 0.969853, 0.957557, 0.989784, 1.01549, 1.02182, 0.976005, 1.01865, 1.06871, 1.06264, 1.02969, 1.07378, 1.06622, 1.15057, 1.0433, 0.83654, 0.847178, 0.890027, 0.920814, 0.888271, 1.04662, 0.8869, 0.856348, 0.863181, 0.906312, 0.902166, 1.00122, 1.03303, 0.887866, 0.892437, 0.906278, 0.884976, 0.864251, 0.917221, 1.10618, 1.04028, 0.893184, 0.915734, 0.892676};
    // calibration T0C
    float calibT0C[28] = {0.949829, 1.05408, 1.00681, 1.00724, 0.990663, 0.973571, 0.9855, 1.03726, 1.02526, 1.00467, 0.983008, 0.979349, 0.952352, 0.985775, 1.013, 1.01721, 0.993948, 0.996421, 0.971871, 1.02921, 0.989641, 1.01885, 1.01259, 0.929502, 1.03969, 1.02496, 1.01385, 1.01711};
    // calibration T0A
    float calibT0A[24] = {0.86041, 1.10607, 1.17724, 0.756397, 1.14954, 1.0879, 0.829438, 1.09014, 1.16515, 0.730077, 1.06722, 0.906344, 0.824167, 1.14716, 1.20692, 0.755034, 1.11734, 1.00556, 0.790522, 1.09138, 1.16225, 0.692458, 1.12428, 1.01127};
    // calibration FDA
    float calibFDA[8] = {0.933485, 1.00743, 0.768484, 0.837354, 1.26397, 1.2159, 0.876259, 1.00434};
    // calibration FDC
    float calibFDC[8] = {0.772909, 1.95841, 0.966258, 0.913508, 0.96176, 0.650286, 0.619638, 0.694932};
    // calibration factor MFT vs vtx
    float biningVtxt[30] = {-14.5, -13.5, -12.5, -11.5, -10.5, -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5};
    float calibMFTvtx[30] = {1.27106, 1.239, 1.23555, 1.21115, 1.18538, 1.16585, 1.15052, 1.11002, 1.09909, 1.07341, 1.05511, 1.04538, 1.02291, 1.01254, 0.994646, 0.981465, 0.962609, 0.955254, 0.942117, 0.932826, 0.925284, 0.906386, 0.904621, 0.887349, 0.884237, 0.860509, 0.853699, 0.836858, 0.819485, 0.802514};

    // calibration factor FV0 vs vtx
    float calibFV0vtx[30] = {0.907962, 0.934607, 0.938929, 0.950987, 0.950817, 0.966362, 0.968509, 0.972741, 0.982412, 0.984872, 0.994543, 0.996003, 0.99435, 1.00266, 0.998245, 1.00584, 1.01078, 1.01003, 1.00726, 1.00872, 1.01726, 1.02015, 1.0193, 1.01106, 1.02229, 1.02104, 1.03435, 1.00822, 1.01921, 1.01736};
    // calibration FT0A vs vtx
    float calibFT0Avtx[30] = {0.924334, 0.950988, 0.959604, 0.965607, 0.970016, 0.979057, 0.978384, 0.982005, 0.992825, 0.990048, 0.998588, 0.997338, 1.00102, 1.00385, 0.99492, 1.01083, 1.00703, 1.00494, 1.00063, 1.0013, 1.00777, 1.01238, 1.01179, 1.00577, 1.01028, 1.017, 1.02975, 1.0085, 1.00856, 1.01662};
    // calibration FT0C vs vtx
    float calibFT0Cvtx[30] = {1.02096, 1.01245, 1.02148, 1.03605, 1.03561, 1.03667, 1.04229, 1.0327, 1.03674, 1.02764, 1.01828, 1.02331, 1.01864, 1.015, 1.01197, 1.00615, 0.996845, 0.993051, 0.985635, 0.982883, 0.981914, 0.964635, 0.967812, 0.95475, 0.956687, 0.932816, 0.92773, 0.914892, 0.891724, 0.872382};
    // calibration FDA vs vtx
    float calibFDAvtx[30] = {1.05852, 1.07943, 1.03542, 1.02851, 1.00617, 1.01377, 0.997411, 1.00899, 0.98556, 0.994506, 0.999267, 0.997915, 0.993361, 0.988509, 0.993386, 0.992661, 0.997844, 1.00428, 0.991939, 0.995139, 0.999882, 1.00976, 1.01239, 0.989125, 1.01432, 1.00652, 1.02114, 1.00582, 0.996546, 1.05708};
    // calibration FDC vs vtx
    float calibFDCvtx[30] = {0.965937, 0.924875, 0.903078, 0.901217, 0.914662, 0.933606, 0.899319, 0.912335, 0.900467, 0.928819, 0.943352, 0.955755, 0.954358, 0.939799, 0.973757, 0.972073, 1.00221, 0.997195, 1.01809, 1.02407, 1.02722, 1.05898, 1.09503, 1.13893, 1.12981, 1.15516, 1.17394, 1.28468, 1.37351, 1.24345};

    const int nDetVtx = 6;
    TGraph* gVtx[nDetVtx];
    const char* nameDet[nDetVtx] = {"AmpV0", "AmpT0A", "AmpT0C", "MFT", "AmpFDA", "AmpFDC"};
    for (int i_d = 0; i_d < nDetVtx; ++i_d) {
      gVtx[i_d] = 0;
      gVtx[i_d] = new TGraph();
    }
    for (int i_v = 0; i_v < 30; ++i_v) {
      gVtx[0]->SetPoint(i_v, biningVtxt[i_v], calibFV0vtx[i_v]);
    }
    for (int i_v = 0; i_v < 30; ++i_v) {
      gVtx[1]->SetPoint(i_v, biningVtxt[i_v], calibFT0Avtx[i_v]);
    }
    for (int i_v = 0; i_v < 30; ++i_v) {
      gVtx[2]->SetPoint(i_v, biningVtxt[i_v], calibFT0Cvtx[i_v]);
    }
    for (int i_v = 0; i_v < 30; ++i_v) {
      gVtx[3]->SetPoint(i_v, biningVtxt[i_v], calibMFTvtx[i_v]);
    }
    for (int i_v = 0; i_v < 30; ++i_v) {
      gVtx[4]->SetPoint(i_v, biningVtxt[i_v], calibFDAvtx[i_v]);
    }
    for (int i_v = 0; i_v < 30; ++i_v) {
      gVtx[5]->SetPoint(i_v, biningVtxt[i_v], calibFDCvtx[i_v]);
    }

    for (int i_d = 0; i_d < nDetVtx; ++i_d) {
      gVtx[i_d]->SetName(Form("g%s", nameDet[i_d]));
    }

    float sumAmpFV0 = 0;
    float sumAmpFV01to4Ch = 0;
    int innerFV0 = 32;
    float maxEtaFV0 = 5.1;
    float minEtaFV0 = 2.2;
    float detaFV0 = (maxEtaFV0 - minEtaFV0) / 5.0;

    const int nCells = 48; // 48 sectors in FV0
    float amp_channel[nCells];
    for (int iCell = 0; iCell < nCells; ++iCell) {
      amp_channel[iCell] = 0.0;
    }
    float amp_channelBefore[nCells];
    for (int iCell = 0; iCell < nCells; ++iCell) {
      amp_channelBefore[iCell] = 0.0;
    }

    if (collision.has_foundFV0()) {

      float RhoLattice[nCells];
      for (Int_t iCh = 0; iCh < nCells; iCh++) {
        RhoLattice[iCh] = 0;
      }
      auto fv0 = collision.foundFV0();
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
        float phiv0 = -999.0;
        float etav0 = -999.0;
        int channelv0 = fv0.channel()[ich];
        float ampl_ch = fv0.amplitude()[ich];
        int ringindex = getFV0Ring(channelv0);
        int channelv0phi = getFV0IndexPhi(channelv0);
        etav0 = maxEtaFV0 - (detaFV0 / 2.0) * (2.0 * ringindex + 1);
        if (channelv0 < innerFV0) {
          phiv0 = (2.0 * (channelv0phi - 8 * ringindex) + 1) * M_PI / (8.0);
        } else {
          phiv0 = ((2.0 * channelv0phi) + 1 - 64.0) * 2.0 * M_PI / (32.0);
        }
        amp_channelBefore[channelv0phi] = ampl_ch;
        if (applyCalibCh) {
          ampl_ch *= calib[channelv0phi];
        }
        sumAmpFV0 += ampl_ch;

        if (channelv0 >= 8) { // exclude the 1st ch, eta 2.2,4.52
          sumAmpFV01to4Ch += ampl_ch;
        }
        flatenicity.fill(HIST("fEtaPhiFv0"), phiv0, etav0, ampl_ch);
        amp_channel[channelv0phi] = ampl_ch;
        if (channelv0 < innerFV0) {
          RhoLattice[channelv0phi] = ampl_ch;
        } else {
          RhoLattice[channelv0phi] = ampl_ch / 2.0; // two channels per bin
        }
      }
      flatenicity_fv0 = GetFlatenicity(RhoLattice, nCells);
      flatenicity.fill(HIST("hAmpV0vsVtxBeforeCalibration"), vtxZ, sumAmpFV0);
      if (applyCalibVtx) {
        sumAmpFV0 *= gVtx[0]->Eval(vtxZ);
        sumAmpFV01to4Ch *= gVtx[0]->Eval(vtxZ);
      }
      flatenicity.fill(HIST("hAmpV0vsVtx"), vtxZ, sumAmpFV0);
    }

    // MFTtracks + flatenicity
    float multMFTTrackParc = 0;
    float multMFTTrack = 0;

    const int nRings1 = 2;
    const int nSectors1 = 8;
    const int nCells1 = nRings1 * nSectors1;
    float maxEta1[nRings1] = {-3.05, -2.50};
    float minEta1[nRings1] = {-3.60, -3.05};
    float maxPhi1[nSectors1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    float minPhi1[nSectors1] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    float RhoLattice1[nCells1];
    for (int iCh = 0; iCh < nCells1; iCh++) {
      RhoLattice1[iCh] = 0.0;
    }

    for (auto& track : mfttracks) {

      float eta_a = track.eta();
      float phi_a = track.phi();
      o2::math_utils::bringTo02Pi(phi_a);
      flatenicity.fill(HIST("hetaMFT"), eta_a);
      flatenicity.fill(HIST("hphiMFT"), phi_a);
      if (eta_a > -2.5 || eta_a < -3.6) {
        continue;
      }

      int i_ch = 0;
      for (int ir = 0; ir < nRings1; ir++) {
        for (int is = 0; is < nSectors1; is++) {
          if (eta_a >= minEta1[ir] && eta_a < maxEta1[ir] &&
              phi_a >= minPhi1[is] * 2.0 * M_PI / (1.0 * nSectors1) &&
              phi_a < maxPhi1[is] * 2.0 * M_PI / (1.0 * nSectors1)) {
            RhoLattice1[i_ch]++;
          }
          i_ch++;
        }
      }

      multMFTTrack++;

      if (eta_a > -3.4 || eta_a < -3.6) {
        continue;
      }
      multMFTTrackParc++;
    }
    flatenicity.fill(HIST("hMFTvsVtxBeforeCalibration"), vtxZ, multMFTTrack);
    if (applyCalibVtx) {
      multMFTTrack *= gVtx[3]->Eval(vtxZ);
      multMFTTrackParc *= gVtx[3]->Eval(vtxZ);
    }
    flatenicity.fill(HIST("hMFTvsVtx"), vtxZ, multMFTTrack);

    float flatenicity_mft = 9999;
    flatenicity_mft = GetFlatenicity(RhoLattice1, nCells1);

    // global tracks + flatenicity
    const int nRings2 = 4;
    const int nSectors2 = 8;
    const int nCells2 = nRings2 * nSectors2;
    float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
    float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
    float maxPhi2[nSectors2] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    float minPhi2[nSectors2] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    float RhoLattice2[nCells2];
    for (int iCh = 0; iCh < nCells2; iCh++) {
      RhoLattice2[iCh] = 0.0;
    }
    int multGlob = 0;
    for (auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
      float eta_a = track.eta();
      float phi_a = track.phi();
      multGlob++;

      int i_ch = 0;
      for (int ir = 0; ir < nRings2; ir++) {
        for (int is = 0; is < nSectors2; is++) {
          if (eta_a >= minEta2[ir] && eta_a < maxEta2[ir] &&
              phi_a >= minPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2) &&
              phi_a < maxPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2)) {
            RhoLattice2[i_ch]++;
          }
          i_ch++;
        }
      }
    }

    float flatenicity_glob = 9999;
    flatenicity_glob = GetFlatenicity(RhoLattice2, nCells2);

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
          flatenicity.fill(HIST("hAmpT0AVsChBeforeCalibration"), sector, amplitude);
          if (applyCalibCh) {
            amplitude *= calibT0A[sector];
          }
          flatenicity.fill(HIST("hAmpT0AVsCh"), sector, amplitude);
        }
        sumAmpFT0A += amplitude;
        flatenicity.fill(HIST("hFT0A"), amplitude);
      }

      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        sumAmpFT0C += amplitude;
        uint8_t channel = ft0.channelC()[i_c];
        int sector = getT0CSector(channel);
        if (sector >= 0 && sector < 28) {
          RhoLatticeT0C[sector] += amplitude;
          flatenicity.fill(HIST("hAmpT0CVsChBeforeCalibration"), sector, amplitude);
          if (applyCalibCh) {
            amplitude *= calibT0C[sector];
          }
          flatenicity.fill(HIST("hAmpT0CVsCh"), sector, amplitude);
        }
        flatenicity.fill(HIST("hFT0C"), amplitude);
      }
      flatenicity.fill(HIST("hAmpT0AvsVtxBeforeCalibration"), vtxZ, sumAmpFT0A);
      flatenicity.fill(HIST("hAmpT0CvsVtxBeforeCalibration"), vtxZ, sumAmpFT0C);
      if (applyCalibVtx) {
        sumAmpFT0A *= gVtx[1]->Eval(vtxZ);
        sumAmpFT0C *= gVtx[2]->Eval(vtxZ);
      }
      flatenicity.fill(HIST("hAmpT0AvsVtx"), vtxZ, sumAmpFT0A);
      flatenicity.fill(HIST("hAmpT0CvsVtx"), vtxZ, sumAmpFT0C);
    }
    float flatenicity_t0a = 9999;
    flatenicity_t0a = GetFlatenicity(RhoLatticeT0A, nCellsT0A);
    float flatenicity_t0c = 9999;
    flatenicity_t0c = GetFlatenicity(RhoLatticeT0C, nCellsT0C);

    // FDD
    float sumAmpFDDA = 0;
    float sumAmpFDDC = 0;
    // LOGP(info, "before sumAmpFDDA={}", sumAmpFDDA);

    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      for (std::size_t ich = 0; ich < 8; ich++) {
        float amplitude = fdd.chargeA()[ich];
        sumAmpFDDA += amplitude;
        flatenicity.fill(HIST("hAmpFDAVsChBeforeCalibration"), ich, amplitude);
        if (applyCalibCh) {
          amplitude *= calibFDA[ich];
        }
        flatenicity.fill(HIST("hAmpFDAVsCh"), ich, amplitude);
      }
      for (std::size_t ich = 0; ich < 8; ich++) {
        float amplitude = fdd.chargeC()[ich];
        sumAmpFDDC += amplitude;
        flatenicity.fill(HIST("hAmpFDCVsChBeforeCalibration"), ich, amplitude);
        if (applyCalibCh) {
          amplitude *= calibFDC[ich];
        }
        flatenicity.fill(HIST("hAmpFDCVsCh"), ich, amplitude);
      }
      flatenicity.fill(HIST("hAmpFDAvsVtxBeforeCalibration"), vtxZ, sumAmpFDDA);
      flatenicity.fill(HIST("hAmpFDCvsVtxBeforeCalibration"), vtxZ, sumAmpFDDC);
      if (applyCalibVtx) {
        sumAmpFDDA *= gVtx[4]->Eval(vtxZ);
        sumAmpFDDC *= gVtx[5]->Eval(vtxZ);
      }
      flatenicity.fill(HIST("hAmpFDAvsVtx"), vtxZ, sumAmpFDDA);
      flatenicity.fill(HIST("hAmpFDCvsVtx"), vtxZ, sumAmpFDDC);
    }

    float combined_estimator1 = 0;
    float combined_estimator2 = 0;
    float combined_estimator3 = 0;
    float combined_estimator4 = 0;
    float combined_estimator5 = 0;
    float combined_estimator6 = 0;

    float estimator[16];
    for (int i_e = 0; i_e < 16; ++i_e) {
      estimator[i_e] = 0;
    }

    if (collision.has_foundFV0() && collision.has_foundFT0() &&
        collision.has_foundFDD()) {
      // option 1
      ampl1[0] = sumAmpFDDC;
      ampl1[1] = multMFTTrackParc;
      ampl1[2] = sumAmpFT0C;
      ampl1[3] = sumAmpFV01to4Ch;
      ampl1[4] = sumAmpFDDA;
      float all_weights = 0;
      if (applyNorm) {
        all_weights = 0;
        for (int i_1 = 0; i_1 < nEta1; ++i_1) {
          combined_estimator1 +=
            ampl1[i_1] * weigthsEta1[i_1] / deltaEeta1[i_1];
          all_weights += weigthsEta1[i_1];
        }
        combined_estimator1 /= all_weights;
      } else {
        for (int i_1 = 0; i_1 < nEta1; ++i_1) {
          combined_estimator1 += ampl1[i_1] * weigthsEta1[i_1];
        }
      }
      // option 2
      ampl2[0] = sumAmpFDDC;
      ampl2[1] = multMFTTrack;
      ampl2[2] = sumAmpFV01to4Ch;
      ampl2[3] = sumAmpFDDA;
      if (applyNorm) {
        all_weights = 0;
        for (int i_2 = 0; i_2 < nEta2; ++i_2) {
          combined_estimator2 +=
            ampl2[i_2] * weigthsEta2[i_2] / deltaEeta2[i_2];
          all_weights += weigthsEta2[i_2];
        }
        combined_estimator2 /= all_weights;
      } else {
        for (int i_2 = 0; i_2 < nEta2; ++i_2) {
          combined_estimator2 += ampl2[i_2] * weigthsEta2[i_2];
        }
      }
      // option 3
      ampl3[0] = multMFTTrack;
      ampl3[1] = sumAmpFV0;
      if (applyNorm) {
        all_weights = 0;
        for (int i_3 = 0; i_3 < nEta3; ++i_3) {
          combined_estimator3 +=
            ampl3[i_3] * weigthsEta3[i_3] / deltaEeta3[i_3];
          all_weights += weigthsEta3[i_3];
        }
        combined_estimator3 /= all_weights;
      } else {
        for (int i_3 = 0; i_3 < nEta3; ++i_3) {
          combined_estimator3 += ampl3[i_3] * weigthsEta3[i_3];
        }
      }
      // option 4
      ampl4[0] = multMFTTrack;
      ampl4[1] = sumAmpFT0A;
      if (applyNorm) {
        all_weights = 0;
        for (int i_4 = 0; i_4 < nEta4; ++i_4) {
          combined_estimator4 +=
            ampl4[i_4] * weigthsEta4[i_4] / deltaEeta4[i_4];
          all_weights += weigthsEta4[i_4];
        }
        combined_estimator4 /= all_weights;
      } else {
        for (int i_4 = 0; i_4 < nEta4; ++i_4) {
          combined_estimator4 += ampl4[i_4] * weigthsEta4[i_4];
        }
      }
      // option 5
      ampl5[0] = sumAmpFT0C;
      ampl5[1] = sumAmpFT0A;
      if (applyNorm) {
        all_weights = 0;
        for (int i_5 = 0; i_5 < nEta5; ++i_5) {
          combined_estimator5 +=
            ampl5[i_5] * weigthsEta5[i_5] / deltaEeta5[i_5];
          all_weights += weigthsEta5[i_5];
        }
        combined_estimator5 /= all_weights;
      } else {
        for (int i_5 = 0; i_5 < nEta5; ++i_5) {
          combined_estimator5 += ampl5[i_5] * weigthsEta5[i_5];
        }
      }

      // option 6: FT0C + FV0
      ampl6[0] = sumAmpFT0C;
      ampl6[1] = sumAmpFV0;
      if (applyNorm) {
        all_weights = 0;
        for (int i_6 = 0; i_6 < nEta6; ++i_6) {
          combined_estimator6 +=
            ampl6[i_6] * weigthsEta6[i_6] / deltaEeta6[i_6];
          all_weights += weigthsEta6[i_6];
        }
        combined_estimator6 /= all_weights;
      } else {
        for (int i_6 = 0; i_6 < nEta6; ++i_6) {
          combined_estimator6 += ampl6[i_6] * weigthsEta6[i_6];
        }
      }

      flatenicity.fill(HIST("hMFTmult"), multMFTTrackParc);
      flatenicity.fill(HIST("hMFTmultAll"), multMFTTrack);
      flatenicity.fill(HIST("hFDDAampl"), sumAmpFDDA);
      flatenicity.fill(HIST("hFDDCampl"), sumAmpFDDC);
      flatenicity.fill(HIST("hFT0Aampl"), sumAmpFT0A);
      flatenicity.fill(HIST("hFT0Campl"), sumAmpFT0C);
      flatenicity.fill(HIST("hFV0amplRing1to4"), sumAmpFV01to4Ch);
      flatenicity.fill(HIST("hEv"), 4);
      float flatenicity_mft_glob = (flatenicity_mft + flatenicity_glob) / 2.0;
      float flatenicity_mft_fv0 = (flatenicity_mft + flatenicity_fv0) / 2.0;
      float flatenicity_mft_glob_fv0 =
        (flatenicity_mft + flatenicity_glob + flatenicity_fv0) / 3.0;
      estimator[0] = multGlob;
      estimator[1] = combined_estimator1;
      estimator[2] = combined_estimator2;
      estimator[3] = combined_estimator3;
      estimator[4] = sumAmpFV0;
      estimator[5] = multMFTTrack;
      estimator[6] = 1.0 - flatenicity_fv0;
      estimator[7] = 1.0 - flatenicity_mft_glob;
      estimator[8] = 1.0 - flatenicity_mft_glob_fv0;
      estimator[9] = 1.0 - flatenicity_mft_fv0;
      estimator[10] = combined_estimator4;
      estimator[11] = combined_estimator5;
      float flatenicity_ft0 = (flatenicity_t0a + flatenicity_t0c) / 2.0;
      estimator[12] = 1.0 - flatenicity_ft0;
      float flatenicity_mft_ft0a = (flatenicity_mft + flatenicity_t0a) / 2.0;
      estimator[13] = 1.0 - flatenicity_mft_ft0a;
      estimator[14] = combined_estimator6;
      estimator[15] = 1.0 - (flatenicity_fv0 + flatenicity_t0c) / 2.0;

      static_for<0, 15>([&](auto i) {
        constexpr int index = i.value;
        flatenicity.fill(HIST(nhEst[index]), estimator[index], estimator[0]);
      });
      flatenicity.fill(HIST("hFlatMFTvsFlatGlob"), flatenicity_mft,
                       flatenicity_glob);
      flatenicity.fill(HIST("hFlatMFTvsFlatFV0"), flatenicity_mft,
                       flatenicity_fv0);
      // plot pt vs estimators
      for (auto& track : tracks) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        float pt = track.pt();
        static_for<0, 15>([&](auto i) {
          constexpr int index = i.value;
          flatenicity.fill(HIST(nhPtEst[index]), estimator[index], pt);
        });
      }

      for (int iCh = 0; iCh < 48; ++iCh) {
        flatenicity.fill(HIST("hAmpV0VsCh"), iCh, amp_channel[iCh]);
        flatenicity.fill(HIST("hAmpV0VsChBeforeCalibration"), iCh,
                         amp_channelBefore[iCh]);
      }
      flatenicity.fill(HIST("fMultFv0"), sumAmpFV0);
      flatenicity.fill(HIST("hFlatFT0CvsFlatFT0A"), flatenicity_t0c, flatenicity_t0a);
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<flatenictyFV0>(cfgc));
  return workflow;
}
