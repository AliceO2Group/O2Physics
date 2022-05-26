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
/// \author Antonio Ortiz (antonio.ortiz.velasquez.@cern.ch)
/// \since May 2022 (flatenicity definition from arXiv:2204.13733)

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom.h>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct flatenictyFV0 {

  // Configurable<int> customDeltaBC{"customDeltaBC", 300, "custom BC delta for FIT-collision matching"};
  Configurable<bool> isMC{"isMC", false, "option to flag mc"};
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<bool> applyCalibFV0{"applyCalibFV0", true, "equalize FV0"};
  Configurable<bool> applyNorm{"applyNorm", false, "normalization to eta"};
  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 1.5f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum  pT"};

  HistogramRegistry flatenicity{"flatenicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  static constexpr std::string_view nhEst[10] = {"eGlobaltrack", "eFDDAFDDCFT0CFV0MFT", "eFDDAFDDCFV0MFT", "eFV0MFT", "eFV0", "eMFTmult", "e1flatencityFV0", "e1flatencitytrkMFT", "e1flatencitytrkMFTFV0", "e1flatencityMFTFV0"};
  static constexpr std::string_view tEst[10] = {"GlobalTrk", "FDDA_FDDC_FT0C_FV0_MFT", "FDDA_FDDC_FV0_MFT", "FV0_MFT", "FV0", "MFTTrk", "1-flatencity_FV0", "1-flatencity_trk_MFT", "1-flatencity_trk_MFT_FV0", "1-flatencity_MFT_FV0"};

  static constexpr std::string_view nhPtEst[10] = {"ptVsGlobaltrack", "ptVsFDDAFDDCFT0CFV0MFT", "ptVsFDDAFDDCFV0MFT", "ptVsFV0MFT", "ptVsFV0", "ptVsMFTmult", "ptVs1flatencityFV0", "ptVs1flatencitytrkMFT", "ptVs1flatencitytrkMFTFV0", "ptVs1flatencityMFTFV0"};

  void init(o2::framework::InitContext&)
  {
    int nBinsEst[10] = {100, 400, 400, 400, 1000, 200, 102, 102, 102, 102};
    float lowEdgeEst[10] = {-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.01, -0.01, -0.01, -0.01};
    float upEdgeEst[10] = {99.5, 399.5, 399.5, 399.5, 39999.5, 199.5, 1.01, 1.01, 1.01, 1.01};

    ConfigurableAxis ptBinning{"ptBinning", {0, 0.0, 0.1, 0.15, 0.3, 0.6, 1.0, 2.0, 4.0, 6.0, 10.0, 20.0, 50.0}, "pTassoc bin limits"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    flatenicity.add("hEv", "Ev", HistType::kTH1F, {{6, -0.5, 5.5, "index activated detector"}});
    flatenicity.add("hFV0amplRing1to4", "FV01to4", HistType::kTH1F, {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hFT0Aampl", "FTAampl", HistType::kTH1F, {{10000, -0.5, +1999.5, "FT0A amplitude"}});
    flatenicity.add("hFT0Campl", "FTCampl", HistType::kTH1F, {{10000, -0.5, +1999.5, "FT0C amplitude"}});
    flatenicity.add("hFDDAampl", "FDDAampl", HistType::kTH1F, {{6000, -0.5, 7999.5, "FDDA amplitude"}});
    flatenicity.add("hFDDCampl", "FDDCampl", HistType::kTH1F, {{6000, -0.5, 7999.5, "FDDC amplitude"}});

    flatenicity.add("hMFTmult", "", HistType::kTH1F, {{200, -0.5, +199.5, "MFT mult (-3.6<#eta<-3.4)"}});
    flatenicity.add("hetaMFT", "", HistType::kTH1F, {{100, -5.0, +5.0, "#eta"}});
    flatenicity.add("hphiMFT", "", HistType::kTH1F, {{8, 0.0, 2.0 * M_PI, "#varphi"}});

    flatenicity.add("hMFTmultAll", "", HistType::kTH1F, {{200, -0.5, +199.5, "MFT mult (-3.6<#eta<-2.5)"}});
    // estimators
    for (int i_e = 0; i_e < 10; ++i_e) {
      flatenicity.add(nhEst[i_e].data(), "", HistType::kTH2F, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}, {100, -0.5, +99.5, "Global track"}});
    }

    // vs pT
    for (int i_e = 0; i_e < 10; ++i_e) {
      flatenicity.add(nhPtEst[i_e].data(), "", HistType::kTH2F, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}, {ptAxis}});
    }

    flatenicity.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
    flatenicity.add("vtxZEta", ";#eta;vtxZ", HistType::kTH2F, {{50, -2.5, 2.5, " "}, {60, -30, 30, " "}});
    flatenicity.add("phiEta", ";#eta;#varphi", HistType::kTH2F, {{50, -2.5, 2.5}, {200, 0., 2 * M_PI, " "}});
    flatenicity.add("hvtxZ", "vtxZ", HistType::kTH1F, {{40, -20.0, 20.0, " "}});
    // QA flatenicity
    flatenicity.add("fMultFv0", "FV0 amp", HistType::kTH1F, {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("fMultFv0Check", "FV0 amp", HistType::kTH1F, {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hAmpVsCh", "", HistType::kTH2F, {{48, -0.5, 47.5, "channel"}, {500, -0.5, +19999.5, "FV0 amplitude"}});
    flatenicity.add("hAmpVsChBeforeCalibration", "", HistType::kTH2F, {{48, -0.5, 47.5, "channel"}, {500, -0.5, +19999.5, "FV0 amplitude"}});

    flatenicity.add("fFlatenicityVsFV0", "", HistType::kTH2F, {{1000, -0.5, +9.5, "flatenicity"}, {1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hFlatMFTvsFlatGlob", "", HistType::kTH2F, {{102, -0.01, +1.01, "flatenicity (Glob)"}, {102, -0.01, +1.01, "flatenicity (MFT)"}});
    flatenicity.add("hFlatMFTvsFlatFV0", "", HistType::kTH2F, {{102, -0.01, +1.01, "flatenicity (FV0)"}, {102, -0.01, +1.01, "flatenicity (MFT)"}});

    flatenicity.add("fMultGlobVsFV0", "", HistType::kTH2F, {{100, -0.5, +99.5, "Trk mult"}, {1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("fEtaPhiFv0", "eta vs phi", HistType::kTH2F, {{8, 0.0, 2 * M_PI, "#phi (rad)"}, {5, 2.2, 5.1, "#eta"}});
    flatenicity.add("fFlatenicityBefore", "flatenicity (before vtx cut)", HistType::kTH1F, {{1000, -0.5, +9.5, "flatenicity"}});
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

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::pt > cfgTrkLowPtCut);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection>>;
  using CollisionTableData = soa::Join<aod::Collisions, aod::FT0sCorrected, aod::EvSels>;

  void process(CollisionTableData::iterator const& collision, TrackCandidates const& tracks, aod::MFTTracks const& mfttracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::FDDs_001 const& fdds_001)
  {

    auto vtxZ = collision.posZ();
    flatenicity.fill(HIST("hEv"), 0);
    flatenicity.fill(HIST("hvtxZ"), vtxZ);

    bool isGoodEvent = false;
    if (vtxZ < -10.f || vtxZ > 10.0f) {
      flatenicity.fill(HIST("hEv"), 1);
      if (collision.has_foundFT0()) {
        flatenicity.fill(HIST("hEv"), 2);
        auto ft0 = collision.foundFT0();
        int triggersignals = ft0.triggerMask();
        if (triggersignals) {
          isGoodEvent = true;
        }
      }
    }

    if (!isGoodEvent) {
      return;
    }

    flatenicity.fill(HIST("hEv"), 3);

    const int nEta1 = 5;
    float weigthsEta1[nEta1] = {0.0126266, 1.90601, 0.101574, 0.00238926, 0.00376294};
    float deltaEeta1[nEta1] = {2.0, 0.2, 1.1, 2.32, 1.6};
    float ampl1[nEta1] = {0, 0, 0, 0, 0};

    const int nEta2 = 4;
    float weigthsEta2[nEta2] = {0.0126266, 2.76128, 0.00238926, 0.00376294};
    float deltaEeta2[nEta2] = {2.0, 1.1, 2.32, 1.6};
    float ampl2[nEta2] = {0, 0, 0, 0};

    const int nEta3 = 2;
    float weigthsEta3[nEta3] = {2.76128, 0.00231781};
    float deltaEeta3[nEta3] = {1.1, 2.9};
    float ampl3[nEta3] = {0, 0};

    // V0A signal and flatenicity calculation
    float flatenicity_fv0;
    float calib[48] = {1.05205, 1.13434, 1.09915, 1.07806, 1.13375, 1.19626, 1.26371, 1.11951, 0.931061, 0.934837, 0.951456, 0.960015, 0.973887, 1.00131, 1.03056, 0.972653, 0.986339, 1.07219, 1.03677, 1.02667, 1.0611, 1.034, 1.10586, 1.08495, 0.822695, 0.876057, 0.869338, 0.92351, 0.854066, 1.14313, 0.889261, 0.83169, 0.841845, 0.868859, 0.884049, 1.06192, 1.03866, 0.925584, 0.885038, 0.932367, 0.918977, 0.827109, 0.887456, 1.13838, 1.03751, 0.872285, 0.926142, 0.867931};

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
        RhoLattice[iCh] = 0.0;
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
        if (applyCalibFV0) {
          ampl_ch *= calib[channelv0phi];
        }
        sumAmpFV0 += ampl_ch;

        if (channelv0 >= 8) { // exclude the first channel, eta coverage 2.2,4.52
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
      flatenicity.fill(HIST("fFlatenicityBefore"), flatenicity_fv0);
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
      flatenicity.fill(HIST("hdNdeta"), eta_a);
      flatenicity.fill(HIST("vtxZEta"), eta_a, vtxZ);
      flatenicity.fill(HIST("phiEta"), eta_a, phi_a);
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
    float sumAmpFT0A = 0;
    float sumAmpFT0C = 0;
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      if (collision.t0ACorrectedValid()) {
        for (auto amplitude : ft0.amplitudeA()) {
          sumAmpFT0A += amplitude;
        }
      }
      if (collision.t0CCorrectedValid()) {
        for (auto amplitude : ft0.amplitudeC()) {
          sumAmpFT0C += amplitude;
        }
      }
    }
    // FDD
    float sumAmpFDDA = 0;
    float sumAmpFDDC = 0;
    // LOGP(info, "before sumAmpFDDA={}", sumAmpFDDA);

    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      for (std::size_t ich = 0; ich < 8; ich++) {
        sumAmpFDDA += fdd.chargeA()[ich];
      }
      for (std::size_t ich = 0; ich < 8; ich++) {
        sumAmpFDDC += fdd.chargeC()[ich];
      }
    }

    float combined_estimator1 = 0;
    float combined_estimator2 = 0;
    float combined_estimator3 = 0;
    float estimator[10];
    for (int i_e = 0; i_e < 10; ++i_e) {
      estimator[i_e] = 0;
    }

    if (collision.has_foundFV0() && collision.has_foundFT0() && collision.has_foundFDD()) {
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
          combined_estimator1 += ampl1[i_1] * weigthsEta1[i_1] / deltaEeta1[i_1];
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
          combined_estimator2 += ampl2[i_2] * weigthsEta2[i_2] / deltaEeta2[i_2];
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
          combined_estimator3 += ampl3[i_3] * weigthsEta3[i_3] / deltaEeta3[i_3];
          all_weights += weigthsEta3[i_3];
        }
        combined_estimator3 /= all_weights;
      } else {
        for (int i_3 = 0; i_3 < nEta3; ++i_3) {
          combined_estimator3 += ampl3[i_3] * weigthsEta3[i_3];
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
      float flatenicity_mft_glob_fv0 = (flatenicity_mft + flatenicity_glob + flatenicity_fv0) / 3.0;
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
      static_for<0, 9>([&](auto i) {
        constexpr int index = i.value;
        flatenicity.fill(HIST(nhEst[index]), estimator[index], estimator[0]);
      });
      flatenicity.fill(HIST("hFlatMFTvsFlatGlob"), flatenicity_mft, flatenicity_glob);
      flatenicity.fill(HIST("hFlatMFTvsFlatFV0"), flatenicity_mft, flatenicity_fv0);
      // plot pt vs estimators
      for (auto& track : tracks) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        float pt = track.pt();
        static_for<0, 9>([&](auto i) {
          constexpr int index = i.value;
          flatenicity.fill(HIST(nhPtEst[index]), estimator[index], pt);
        });
      }
      float multcheck = 0;
      for (int iCh = 0; iCh < 48; ++iCh) {
        flatenicity.fill(HIST("hAmpVsCh"), iCh, amp_channel[iCh]);
        flatenicity.fill(HIST("hAmpVsChBeforeCalibration"), iCh, amp_channelBefore[iCh]);
        multcheck += amp_channel[iCh];
      }
      flatenicity.fill(HIST("fMultFv0Check"), multcheck);
      flatenicity.fill(HIST("fMultFv0"), sumAmpFV0);
      flatenicity.fill(HIST("fFlatenicityVsFV0"), flatenicity_fv0, sumAmpFV0);
      flatenicity.fill(HIST("fMultGlobVsFV0"), multGlob, sumAmpFV0);
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<flatenictyFV0>(cfgc));
  return workflow;
}
