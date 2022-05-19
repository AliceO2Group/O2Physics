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

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 1.5f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum  pT"};
  // flatenicity cuts
  Configurable<float> cfgRhoMinEtaCut{"cfgRhoMinEtaCut", 2.2f, "min eta for fv0 cells"};
  Configurable<float> cfgRhoMaxEtaCut{"cfgRhoMaxEtaCut", 5.1f, "max eta for fv0 cells"};
  Configurable<float> cfgRhoMinMult{"cfgRhoMinMult", 50.0f, "max eta for fv0 cells"};

  HistogramRegistry flatenicity{"flatenicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  static constexpr std::string_view nMultSect[48] = {"hCh0", "hCh1", "hCh2", "hCh3", "hCh4", "hCh5", "hCh6", "hCh7", "hCh8", "hCh9", "hCh10", "hCh11", "hCh12", "hCh13", "hCh14", "hCh15", "hCh16", "hCh17", "hCh18", "hCh19", "hCh20", "hCh21", "hCh22", "hCh23", "hCh24", "hCh25", "hCh26", "hCh27", "hCh28", "hCh29", "hCh30", "hCh31", "hCh32", "hCh33", "hCh34", "hCh35", "hCh36", "hCh37", "hCh38", "hCh39", "hCh40", "hCh41", "hCh42", "hCh43", "hCh44", "hCh45", "hCh46", "hCh47"};

  void init(o2::framework::InitContext&)
  {

    ConfigurableAxis ptBinning{"ptBinning", {0, 0.0, 0.1, 0.15, 0.3, 0.6, 1.0, 2.0, 4.0, 6.0, 10.0, 20.0, 50.0}, "pTassoc bin limits"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    ConfigurableAxis etaBinning{"etaBinning", {0, -7.00, -6.90, -4.90, -3.60, -3.40, -2.30, +2.20, +4.52, +4.70, +6.30, +7.00}, "eta bin limits"};
    AxisSpec etaAxis = {etaBinning, "#eta"};

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
    flatenicity.add("hOpt0", "", HistType::kTH2F, {{100, -0.5, +99.5, "Global track"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt1", "", HistType::kTH2F, {{400, -0.5, +399.5, "Combined FDDA+FDDC+FT0C+FV0+MFT"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt2", "", HistType::kTH2F, {{400, -0.5, +399.5, "Combined FDDA+FDDC+FV0+MFT"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt3", "", HistType::kTH2F, {{400, -0.5, +399.5, "Combined FV0+MFT"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt4", "", HistType::kTH2F, {{1000, -0.5, +39999.5, "FV0"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt5", "", HistType::kTH2F, {{200, -0.5, +199.5, "MFT mult (-3.6<#eta<-3.4)"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt6", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (FV0)"}, {100, -0.5, +99.5, "Global track"}});

    flatenicity.add("hOpt7", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (trk+MFT)"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt8", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (trk+MFT+FV0)"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt9", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (MFT+FV0)"}, {100, -0.5, +99.5, "Global track"}});
    // vs pT
    flatenicity.add("hpTVsOpt0", "", HistType::kTH2F, {{100, -0.5, +99.5, "Global mult"}, {ptAxis}});
    flatenicity.add("hpTVsOpt1", "", HistType::kTH2F, {{400, -0.5, +399.5, "Combined FDDA+FDDC+FT0C+FV0+MFT"}, {ptAxis}});
    flatenicity.add("hpTVsOpt2", "", HistType::kTH2F, {{400, -0.5, +399.5, "Combined FDDA+FDDC+FV0+MFT"}, {ptAxis}});
    flatenicity.add("hpTVsOpt3", "", HistType::kTH2F, {{400, -0.5, +399.5, "Combined FV0+MFT"}, {ptAxis}});
    flatenicity.add("hpTVsOpt4", "", HistType::kTH2F, {{1000, -0.5, +39999.5, "FV0"}, {ptAxis}});
    flatenicity.add("hpTVsOpt5", "", HistType::kTH2F, {{200, -0.5, +199.5, "MFT mult (-3.6<#eta<-3.4)"}, {ptAxis}});
    flatenicity.add("hpTVsOpt6", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (FV0)"}, {ptAxis}});
    flatenicity.add("hpTVsOpt7", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (trk+MFT)"}, {ptAxis}});
    flatenicity.add("hpTVsOpt8", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (trk+MFT+FV0)"}, {ptAxis}});
    flatenicity.add("hpTVsOpt9", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (MFT+FV0)"}, {ptAxis}});

    flatenicity.add("hOpt7b", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (ftrk<0.9&fMFT<0.9)"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt8b", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (ftrk<0.9&fMFT<0.9&&FV0<0.9)"}, {100, -0.5, +99.5, "Global track"}});
    flatenicity.add("hOpt9b", "", HistType::kTH2F, {{102, -0.01, +1.01, "1-flatencity (fMFT<0.9&fFV0<0.9)"}, {100, -0.5, +99.5, "Global track"}});

    flatenicity.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
    flatenicity.add("vtxZEta", ";#eta;vtxZ", HistType::kTH2F, {{50, -2.5, 2.5, " "}, {60, -30, 30, " "}});
    flatenicity.add("phiEta", ";#eta;#varphi", HistType::kTH2F, {{50, -2.5, 2.5}, {200, 0., 2 * M_PI, " "}});
    flatenicity.add("hvtxZ", "vtxZ", HistType::kTH1F, {{40, -20.0, 20.0, " "}});
    flatenicity.add("hCounter", "Counter; sel; Nev", HistType::kTH1D, {{3, 0, 3, " "}});
    // QA flatenicity
    flatenicity.add("fMultFv0", "FV0 amp", HistType::kTH1F, {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("fMultFv0Sect", "FV0 amp sectors", HistType::kTProfile, {{48, -0.5, +47.5, "sectors"}});
    flatenicity.add("fMultFv0Check", "FV0 amp", HistType::kTH1F, {{1000, -0.5, +39999.5, "FV0 amplitude"}});
    for (int iCh = 0; iCh < 48; ++iCh) {
      flatenicity.add(nMultSect[iCh].data(), "", HistType::kTH1F, {{500, -0.5, +19999.5, "FV0 amplitude"}});
    }
    flatenicity.add("fFlatenicityVsFV0", "", HistType::kTH2F, {{1000, -0.5, +9.5, "flatenicity"}, {1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("hFlatMFTvsFlatGlob", "", HistType::kTH2F, {{102, -0.01, +1.01, "flatenicity (Glob)"}, {102, -0.01, +1.01, "flatenicity (MFT)"}});
    flatenicity.add("hFlatMFTvsFlatFV0", "", HistType::kTH2F, {{102, -0.01, +1.01, "flatenicity (FV0)"}, {102, -0.01, +1.01, "flatenicity (MFT)"}});

    flatenicity.add("fMultGlobVsFV0", "", HistType::kTH2F, {{100, -0.5, +99.5, "Trk mult"}, {1000, -0.5, +39999.5, "FV0 amplitude"}});
    flatenicity.add("fEtaPhiFv0", "eta vs phi", HistType::kTH2F, {{8, 0.0, 2 * M_PI, "#phi (rad)"}, {5, 2.2, 5.1, "#eta"}});
    flatenicity.add("fFlatenicityBefore", "flatenicity (before vtx cut)", HistType::kTH1F, {{1000, -0.5, +9.5, "flatenicity"}});
    flatenicity.add("fFlatenicityAfter", "flatenicity (after vtx cut)", HistType::kTH1F, {{1000, -0.5, +9.5, "flatenicity"}});
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

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::pt > cfgTrkLowPtCut);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection>>;
  void process(soa::Join<aod::Collisions, aod::FT0sCorrected, aod::EvSels>::iterator const& collision, TrackCandidates const& tracks, aod::MFTTracks const& mfttracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::FDDs_001 const& fdds_001)
  {
    const int nEta1 = 5;
    float weigthsEta1[nEta1] = {0.0115, 1.4991, 0.1104, 0.0023, 0.0037};
    float deltaEeta1[nEta1] = {2.0, 0.2, 1.1, 2.32, 1.6};
    float ampl1[nEta1] = {0, 0, 0, 0, 0};

    const int nEta2 = 4;
    float weigthsEta2[nEta2] = {0.0115, 1.07908, 0.0023, 0.0037};
    float deltaEeta2[nEta2] = {2.0, 1.1, 2.32, 1.6};
    float ampl2[nEta2] = {0, 0, 0, 0};

    const int nEta3 = 2;
    float weigthsEta3[nEta3] = {1.07908, 0.002322};
    float deltaEeta3[nEta3] = {1.1, 2.9};
    float ampl3[nEta3] = {0, 0};

    // V0A signal and flatenicity calculation
    float flatenicity_fv0;

    float sumAmpFV0 = 0;
    float sumAmpFV01to4Ch = 0;
    int innerFV0 = 32;
    bool hasValidV0 = false;
    float detaFV0 = (cfgRhoMaxEtaCut - cfgRhoMinEtaCut) / 5.0;

    const int nCells = 48; // 48 sectors in FV0
    double amp_channel[nCells];
    for (int iCell = 0; iCell < nCells; ++iCell) {
      amp_channel[iCell] = 0.0;
    }

    if (collision.has_foundFV0()) {

      float RhoLattice[nCells];
      for (Int_t iCh = 0; iCh < nCells; iCh++) {
        RhoLattice[iCh] = 0.0;
      }
      hasValidV0 = true;
      auto fv0 = collision.foundFV0();
      // LOGP(info, "amplitude.size()={}", fv0.amplitude().size());
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
        sumAmpFV0 += fv0.amplitude()[ich];
        double phiv0 = -999.0;
        double etav0 = -999.0;
        int channelv0 = fv0.channel()[ich];
        if (channelv0 >= 8) { // exclude the first channel, eta coverage 2.2,4.52
          sumAmpFV01to4Ch += fv0.amplitude()[ich];
        }
        int ringindex = getFV0Ring(channelv0);
        int channelv0phi = getFV0IndexPhi(channelv0);
        etav0 = cfgRhoMaxEtaCut - (detaFV0 / 2.0) * (2.0 * ringindex + 1);
        if (channelv0 < innerFV0) {
          phiv0 = (2.0 * (channelv0phi - 8 * ringindex) + 1) * M_PI / (8.0);
        } else {
          phiv0 = ((2.0 * channelv0phi) + 1 - 64.0) * 2.0 * M_PI / (32.0);
        }
        flatenicity.fill(HIST("fEtaPhiFv0"), phiv0, etav0, fv0.amplitude()[ich]);
        amp_channel[channelv0phi] = fv0.amplitude()[ich];
        if (channelv0 < innerFV0) {
          RhoLattice[channelv0phi] = fv0.amplitude()[ich];
        } else {
          RhoLattice[channelv0phi] = fv0.amplitude()[ich] / 2.0; // two channels per bin
        }
      }

      double mRho = 0;
      for (int iCell = 0; iCell < nCells; ++iCell) {
        mRho += 1.0 * RhoLattice[iCell];
      }
      // average activity per cell
      mRho /= (1.0 * nCells);
      // get sigma
      double sRho_tmp = 0;
      for (int iCell = 0; iCell < nCells; ++iCell) {
        sRho_tmp += TMath::Power(1.0 * RhoLattice[iCell] - mRho, 2);
      }
      sRho_tmp /= (1.0 * nCells * nCells);
      double sRho = TMath::Sqrt(sRho_tmp);
      flatenicity_fv0 = 9999;
      if (mRho > 0) {
        flatenicity_fv0 = sRho / mRho;
      }
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

    double mRho_mft = 0;
    for (int iCell = 0; iCell < nCells1; ++iCell) {
      mRho_mft += 1.0 * RhoLattice1[iCell];
    }
    // average activity per cell
    mRho_mft /= (1.0 * nCells1);
    // get sigma
    double sRho_mft_tmp = 0;
    for (int iCell = 0; iCell < nCells1; ++iCell) {
      sRho_mft_tmp += TMath::Power(1.0 * RhoLattice1[iCell] - mRho_mft, 2);
    }
    sRho_mft_tmp /= (1.0 * nCells1 * nCells1);
    double sRho_mft = TMath::Sqrt(sRho_mft_tmp);
    float flatenicity_mft = 9999;
    if (mRho_mft > 0) {
      flatenicity_mft = sRho_mft / mRho_mft;
    }
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
    double mRho_glob = 0;
    for (int iCell = 0; iCell < nCells2; ++iCell) {
      mRho_glob += 1.0 * RhoLattice2[iCell];
    }
    // average activity per cell
    mRho_glob /= (1.0 * nCells2);
    // get sigma
    double sRho_glob_tmp = 0;
    for (int iCell = 0; iCell < nCells2; ++iCell) {
      sRho_glob_tmp += TMath::Power(1.0 * RhoLattice2[iCell] - mRho_glob, 2);
    }
    sRho_glob_tmp /= (1.0 * nCells2 * nCells2);
    double sRho_glob = TMath::Sqrt(sRho_glob_tmp);
    float flatenicity_glob = 9999;
    if (mRho_glob > 0) {
      flatenicity_glob = sRho_glob / mRho_glob;
    }
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

    flatenicity.fill(HIST("hCounter"), 0);

    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }

    float combined_estimator1 = -1;
    float combined_estimator2 = -1;
    float combined_estimator3 = -1;
    if (collision.has_foundFDD() && collision.has_foundFT0() && collision.has_foundFV0()) {
      // option 1
      ampl1[0] = sumAmpFDDC;
      ampl1[1] = multMFTTrackParc;
      ampl1[2] = sumAmpFT0C;
      ampl1[3] = sumAmpFV01to4Ch;
      ampl1[4] = sumAmpFDDA;
      float all_weights = 0;
      for (int i_1 = 0; i_1 < nEta1; ++i_1) {
        combined_estimator1 += ampl1[i_1] * weigthsEta1[i_1] / deltaEeta1[i_1];
        all_weights += weigthsEta1[i_1];
      }
      combined_estimator1 /= all_weights;
      // option 2
      ampl2[0] = sumAmpFDDC;
      ampl2[1] = multMFTTrack;
      ampl2[2] = sumAmpFV01to4Ch;
      ampl2[3] = sumAmpFDDA;
      all_weights = 0;
      for (int i_2 = 0; i_2 < nEta2; ++i_2) {
        combined_estimator2 += ampl2[i_2] * weigthsEta2[i_2] / deltaEeta2[i_2];
        all_weights += weigthsEta2[i_2];
      }
      combined_estimator2 /= all_weights;
      // option 3
      ampl3[0] = multMFTTrack;
      ampl3[1] = sumAmpFV0;
      all_weights = 0;
      for (int i_3 = 0; i_3 < nEta3; ++i_3) {
        combined_estimator3 += ampl3[i_3] * weigthsEta3[i_3] / deltaEeta3[i_3];
        all_weights += weigthsEta3[i_3];
      }
      combined_estimator3 /= all_weights;

      flatenicity.fill(HIST("hMFTmult"), multMFTTrackParc);
      flatenicity.fill(HIST("hMFTmultAll"), multMFTTrack);
      flatenicity.fill(HIST("hEv"), 2);
      flatenicity.fill(HIST("hFDDAampl"), sumAmpFDDA);
      flatenicity.fill(HIST("hFDDCampl"), sumAmpFDDC);
      flatenicity.fill(HIST("hEv"), 1);
      flatenicity.fill(HIST("hEv"), 4);
      flatenicity.fill(HIST("hFT0Aampl"), sumAmpFT0A);
      flatenicity.fill(HIST("hFT0Campl"), sumAmpFT0C);
      flatenicity.fill(HIST("hEv"), 5);
      flatenicity.fill(HIST("hFV0amplRing1to4"), sumAmpFV01to4Ch);
      flatenicity.fill(HIST("hEv"), 3);
      flatenicity.fill(HIST("hOpt0"), multGlob, multGlob);
      flatenicity.fill(HIST("hOpt1"), combined_estimator1, multGlob);
      flatenicity.fill(HIST("hOpt2"), combined_estimator2, multGlob);
      flatenicity.fill(HIST("hOpt3"), combined_estimator3, multGlob);
      flatenicity.fill(HIST("hOpt4"), sumAmpFV0, multGlob);
      flatenicity.fill(HIST("hOpt5"), multMFTTrack, multGlob);
      flatenicity.fill(HIST("hOpt6"), 1.0 - flatenicity_fv0, multGlob);
      flatenicity.fill(HIST("hFlatMFTvsFlatGlob"), flatenicity_mft, flatenicity_glob);
      flatenicity.fill(HIST("hFlatMFTvsFlatFV0"), flatenicity_mft, flatenicity_fv0);
      float flatenicity_mft_glob = (flatenicity_mft + flatenicity_glob) / 2.0;
      float flatenicity_mft_fv0 = (flatenicity_mft + flatenicity_fv0) / 2.0;
      float flatenicity_mft_glob_fv0 = (flatenicity_mft + flatenicity_glob + flatenicity_fv0) / 3.0;
      flatenicity.fill(HIST("hOpt7"), 1.0 - flatenicity_mft_glob, multGlob);
      flatenicity.fill(HIST("hOpt8"), 1.0 - flatenicity_mft_glob_fv0, multGlob);
      flatenicity.fill(HIST("hOpt9"), 1.0 - flatenicity_mft_fv0, multGlob);
      if (flatenicity_mft < 0.9 && flatenicity_glob < 0.9) {
        flatenicity.fill(HIST("hOpt7b"), 1.0 - flatenicity_mft_glob, multGlob);
      }
      if (flatenicity_mft < 0.9 && flatenicity_fv0 < 0.9) {
        flatenicity.fill(HIST("hOpt9b"), 1.0 - flatenicity_mft_glob, multGlob);
      }
      if (flatenicity_mft < 0.9 && flatenicity_fv0 < 0.9 && flatenicity_glob < 0.9) {
        flatenicity.fill(HIST("hOpt8b"), 1.0 - flatenicity_mft_glob_fv0, multGlob);
      }
      // plot pt vs estimators
      for (auto& track : tracks) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        float pt = track.pt();
        flatenicity.fill(HIST("hpTVsOpt0"), multGlob, pt);
        flatenicity.fill(HIST("hpTVsOpt1"), combined_estimator1, pt);
        flatenicity.fill(HIST("hpTVsOpt2"), combined_estimator2, pt);
        flatenicity.fill(HIST("hpTVsOpt3"), combined_estimator3, pt);
        flatenicity.fill(HIST("hpTVsOpt4"), sumAmpFV0, pt);
        flatenicity.fill(HIST("hpTVsOpt5"), multMFTTrack, pt);
        flatenicity.fill(HIST("hpTVsOpt6"), 1.0 - flatenicity_fv0, pt);
        flatenicity.fill(HIST("hpTVsOpt7"), 1.0 - flatenicity_mft_glob, pt);
        flatenicity.fill(HIST("hpTVsOpt8"), 1.0 - flatenicity_mft_glob_fv0, pt);
        flatenicity.fill(HIST("hpTVsOpt9"), 1.0 - flatenicity_mft_fv0, pt);
      }
    }

    if (hasValidV0) {
      double multcheck = 0;
      static_for<0, 47>([&](auto i) {
        constexpr int index = i.value;
        flatenicity.fill(HIST(nMultSect[index]), amp_channel[index]);
        flatenicity.fill(HIST("fMultFv0Sect"), index, amp_channel[index]);
        multcheck += amp_channel[index];
      });

      flatenicity.fill(HIST("fMultFv0Check"), multcheck);
      flatenicity.fill(HIST("fMultFv0"), sumAmpFV0);
      flatenicity.fill(HIST("fFlatenicityVsFV0"), flatenicity_fv0, sumAmpFV0);
      flatenicity.fill(HIST("fMultGlobVsFV0"), multGlob, sumAmpFV0);
    }
    if (hasValidV0 && sumAmpFV0 > cfgRhoMinMult) {
      flatenicity.fill(HIST("fFlatenicityAfter"), flatenicity_fv0);
    }

    auto vtxZ = collision.posZ();

    flatenicity.fill(HIST("hCounter"), 1);

    flatenicity.fill(HIST("hvtxZ"), vtxZ);
    // loop over selected tracks
    for (auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }

      flatenicity.fill(HIST("hdNdeta"), track.eta());
      flatenicity.fill(HIST("vtxZEta"), track.eta(), vtxZ);
      flatenicity.fill(HIST("phiEta"), track.eta(), track.phi());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<flatenictyFV0>(cfgc));
  return workflow;
}
