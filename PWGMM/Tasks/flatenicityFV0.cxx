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
  Configurable<float> cfgRhoMinMult{"cfgRhoMinMult", 300.0f, "max eta for fv0 cells"};

  HistogramRegistry flatenicity{"flatenicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  static constexpr std::string_view nMultSect[48] = {"hCh0", "hCh1", "hCh2", "hCh3", "hCh4", "hCh5", "hCh6", "hCh7", "hCh8", "hCh9", "hCh10", "hCh11", "hCh12", "hCh13", "hCh14", "hCh15", "hCh16", "hCh17", "hCh18", "hCh19", "hCh20", "hCh21", "hCh22", "hCh23", "hCh24", "hCh25", "hCh26", "hCh27", "hCh28", "hCh29", "hCh30", "hCh31", "hCh32", "hCh33", "hCh34", "hCh35", "hCh36", "hCh37", "hCh38", "hCh39", "hCh40", "hCh41", "hCh42", "hCh43", "hCh44", "hCh45", "hCh46", "hCh47"};

  void init(o2::framework::InitContext&)
  {

    ConfigurableAxis ptBinning{"ptBinning", {0, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0}, "pTassoc bin limits"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T}^{assoc} (GeV/#it{c})"};

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
  void process(soa::Join<aod::Collisions, aod::FT0sCorrected, aod::EvSels>::iterator const& collision, TrackCandidates const& tracks, aod::FV0As const& fv0s)
  {

    const int nCells = 48; // 48 sectors in FV0
    float RhoLattice[nCells];
    for (Int_t iCh = 0; iCh < nCells; iCh++) {
      RhoLattice[iCh] = 0.0;
    }
    float rho_m;
    // V0A signal and rho calculation
    float sumAmpFV0 = 0;
    int innerFV0 = 32;
    bool hasValidV0 = false;
    float detaFV0 = (cfgRhoMaxEtaCut - cfgRhoMinEtaCut) / 5.0;

    double amp_channel[nCells];
    for (int iCell = 0; iCell < nCells; ++iCell) {
      amp_channel[iCell] = 0.0;
    }

    if (collision.has_foundFV0()) {
      hasValidV0 = true;
      auto fv0 = collision.foundFV0();
      // LOGP(info, "amplitude.size()={}", fv0.amplitude().size());
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
        sumAmpFV0 += fv0.amplitude()[ich];
        double phiv0 = -999.0;
        double etav0 = -999.0;
        int channelv0 = fv0.channel()[ich];
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
      sRho_tmp /= (1.0 * nCells);
      double sRho = TMath::Sqrt(sRho_tmp);
      rho_m = sRho / mRho;

      flatenicity.fill(HIST("fFlatenicityBefore"), rho_m);
    }

    flatenicity.fill(HIST("hCounter"), 0);

    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }
    int multGlob = 0;
    for (auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
      multGlob++;
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
      flatenicity.fill(HIST("fFlatenicityVsFV0"), rho_m, sumAmpFV0);
      flatenicity.fill(HIST("fMultGlobVsFV0"), multGlob, sumAmpFV0);
    }
    if (hasValidV0 && sumAmpFV0 > cfgRhoMinMult) {
      flatenicity.fill(HIST("fFlatenicityAfter"), rho_m);
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
