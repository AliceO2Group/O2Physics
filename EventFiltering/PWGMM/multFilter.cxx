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

// \file multFilter.cxx
// \brief Selection of events with high activity at mid-pseudorapidity, forward (FV0), or with a leading particle
//
// \author Antonio Ortiz, UNAM Mx, antonio.ortiz.velasquez@cern.ch

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "../filterTables.h"
#include "Framework/HistogramRegistry.h"
#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> mmObjectsNames{"kHtPt", "kHmFv0", "kHmTrk"};

struct multFilter {
  enum { kLeadingPtTrack = 0,
         kHighMultFv0,
         kHighTrackMult,
         kNtriggersMM };

  // event selection cuts
  Configurable<float> selPtTrig{"selPtTrig", 5., "Minimum track pT leading threshold"};
  Configurable<float> selHMFv0{"selHMFv0", 15000., "Minimum FV0-amplitude  threshold"};
  Configurable<float> selHTrkMult{"selHTrkMult", 24., "Minimum charged particle multiplicity threshold"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};
  // rho cuts
  Configurable<int> cfgRhoNetaBins{"cfgRhoNetaBins", 5, "number of bins in eta"};
  Configurable<int> cfgRhoNphiBins{"cfgRhoNphiBins", 8, "number of bins in phi"};
  Configurable<float> cfgRhoMinEtaCut{"cfgRhoMinEtaCut", 2.2f, "min eta for fv0 cells"};
  Configurable<float> cfgRhoMaxEtaCut{"cfgRhoMaxEtaCut", 5.1f, "max eta for fv0 cells"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  void init(o2::framework::InitContext&)
  {

    multiplicity.add("fCollZpos", "Vtx_z", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});
    multiplicity.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
    multiplicity.add("hPhi", "Phi", HistType::kTH1F, {{20, -2.0 * M_PI, 2.0 * M_PI, " "}});
    multiplicity.add("fLeadingTrackPt", "Lead trk pT", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fMultFv0", "FV0 amp", HistType::kTH1F, {{800, -0.5, +39999.5, "FV0 amplitude"}});
    multiplicity.add("fTrackMult", "trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}});
    multiplicity.add("fRhoVsTrkMult", "rhoN vs trk mult", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {1000, -0.5, +9.5, "#rho"}});
    multiplicity.add("fTrackMultVsV0A", "trk mult vs FV0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {800, -0.5, +39999.5, "sum AmpFV0"}});
    multiplicity.add("fTrackMultVsT0A", "trk mult vs FT0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {200, -0.5, +999.5, "sum AmpFT0"}});
    multiplicity.add("fLeadingTrackPtSelected", "sel lead trk pT", HistType::kTH1F, {{150, 0., +150., "trk #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fMultFv0Selected", "sel FV0 amp", HistType::kTH1F, {{800, -0.5, +39999.5, "FV0 amplitude"}});
    multiplicity.add("fTrackMultSelected", "sel trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}});
    // control histograms vs trk mult
    multiplicity.add("fMpTVsTrackMult", "Av pT vs trk mult", HistType::kTProfile, {{200, -0.5, +199.5, "trk mult"}});
    multiplicity.add("fMpTVsTrackMultFlat", "Av pT (#eta<0) vs trk mult (flat)", HistType::kTProfile, {{200, -0.5, +199.5, "trk mult"}});
    // detector response FV0
    multiplicity.add("fEtaPhiFv0", "eta vs phi", HistType::kTH2F, {{8, 0.0, 2 * M_PI, "#phi (rad)"}, {5, 2.2, 5.1, "#eta"}});

    std::array<std::string, 2> eventTitles = {"all", "rejected"};

    auto scalers{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalers->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalers->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
  }

  // declare filters on tracks and charged jets
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::pt > cfgTrkLowPtCut);

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;
  float computeDeltaPhi(float phia, float phib,
                        float rangeMin = -M_PI / 2.0, float rangeMax = 3.0 * M_PI / 2.0)
  {
    float dphi = -999;
    if (phia < 0) {
      phia += 2 * M_PI;
    } else if (phia > 2 * M_PI) {
      phia -= 2 * M_PI;
    }
    if (phib < 0) {
      phib += 2 * M_PI;
    } else if (phib > 2 * M_PI) {
      phib -= 2 * M_PI;
    }
    dphi = phib - phia;
    if (dphi < rangeMin) {
      dphi += 2 * M_PI;
    } else if (dphi > rangeMax) {
      dphi -= 2 * M_PI;
    }
    return dphi;
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

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>::iterator const& collision, TrackCandidates const& tracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {

    bool keepEvent[kNtriggersMM]{false};

    multiplicity.fill(HIST("fProcessedEvents"), 0);
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    // global observables
    int multTrack = 0;
    float flPt = 0; // leading pT
    double rho_m = -1;

    const int nEta = cfgRhoNetaBins;
    const int nPhi = cfgRhoNphiBins;
    double EtaBins[nEta + 1];
    double deltaEta = TMath::Abs(cfgRhoMaxEtaCut - cfgRhoMinEtaCut) / (1.0 * cfgRhoNetaBins);
    for (int i_eta = 0; i_eta < cfgRhoNetaBins + 1; ++i_eta) {
      EtaBins[i_eta] = 0;
      if (i_eta < cfgRhoNetaBins) {
        EtaBins[i_eta] = i_eta * deltaEta + 1.0 * cfgRhoMinEtaCut;
      } else {
        EtaBins[i_eta] = 1.0 * cfgRhoMaxEtaCut;
      }
    }
    double PhiBins[nPhi + 1];
    double deltaPhi = (2.0 * M_PI) / (1.0 * cfgRhoNphiBins);
    for (int i_phi = 0; i_phi < cfgRhoNphiBins + 1; ++i_phi) {
      PhiBins[i_phi] = 0;
      if (i_phi < cfgRhoNphiBins) {
        PhiBins[i_phi] = i_phi * deltaPhi;
      } else {
        PhiBins[i_phi] = 2.0 * M_PI;
      }
    }
    int RhoLattice[nEta][nPhi];
    for (int i_eta = 0; i_eta < cfgRhoNetaBins; ++i_eta) {
      for (int i_phi = 0; i_phi < cfgRhoNphiBins; ++i_phi) {
        RhoLattice[i_eta][i_phi] = 0;
      }
    }
    // V0A signal and rho calculation
    float sumAmpFV0 = 0;
    int innerFV0 = 32;
    bool hasValidV0 = false;
    float detaFV0 = (cfgRhoMaxEtaCut - cfgRhoMinEtaCut) / 5.0;
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
        etav0 = cfgRhoMaxEtaCut - (detaFV0 / 2.0) * (2.0 * ringindex + 1);
        if (channelv0 < innerFV0) {
          phiv0 = (2.0 * (channelv0 - 8 * ringindex) + 1) * M_PI / (8.0);
        } else {
          phiv0 = ((2.0 * channelv0) + 1 - 64.0) * 2.0 * M_PI / (32.0);
        }
        multiplicity.fill(HIST("fEtaPhiFv0"), phiv0, etav0);
        for (int i_eta = 0; i_eta < cfgRhoNetaBins; ++i_eta) {
          for (int i_phi = 0; i_phi < cfgRhoNphiBins; ++i_phi) {
            if (etav0 >= EtaBins[i_eta] && etav0 < EtaBins[i_eta + 1] &&
                phiv0 >= PhiBins[i_phi] && phiv0 < PhiBins[i_phi + 1]) {
              if (channelv0 < innerFV0) {
                RhoLattice[i_eta][i_phi] = fv0.amplitude()[ich];
              } else {
                RhoLattice[i_eta][i_phi] = fv0.amplitude()[ich] / 2.0; // two channels per bin
              }
            }
          }
        }
      }

      double mRho = 0;
      for (int i_eta = 0; i_eta < cfgRhoNetaBins; ++i_eta) {
        for (int i_phi = 0; i_phi < cfgRhoNphiBins; ++i_phi) {
          mRho += 1.0 * RhoLattice[i_eta][i_phi];
        }
      }
      // average activity per cell
      mRho /= (1.0 * cfgRhoNetaBins * cfgRhoNphiBins);
      // get sigma
      double sRho_tmp = 0;
      for (int i_eta = 0; i_eta < cfgRhoNetaBins; ++i_eta) {
        for (int i_phi = 0; i_phi < cfgRhoNphiBins; ++i_phi) {
          sRho_tmp += TMath::Power(1.0 * RhoLattice[i_eta][i_phi] - mRho, 2);
        }
      }
      sRho_tmp /= (1.0 * cfgRhoNetaBins * cfgRhoNphiBins);
      double sRho = TMath::Sqrt(sRho_tmp);
      rho_m = sRho / mRho;
    }

    // T0 signal
    float sumAmpFT0 = 0;
    bool hasValidT0 = false;
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      if (collision.t0ACorrectedValid()) {
        hasValidT0 = true;
        for (auto amplitude : ft0.amplitudeA()) {
          sumAmpFT0 += amplitude;
        }
      }
    }

    for (auto& track : tracks) {
      multTrack++;
      multiplicity.fill(HIST("hdNdeta"), track.eta());
      multiplicity.fill(HIST("hPhi"), track.phi());
      if (flPt < track.pt()) {
        flPt = track.pt();
      }
    }
    // filling tprofiles
    for (auto& track : tracks) {
      multiplicity.fill(HIST("fMpTVsTrackMult"), multTrack, track.pt());
      if (rho_m < 0.8 && collision.has_foundFV0()) {
        multiplicity.fill(HIST("fMpTVsTrackMultFlat"), multTrack, track.pt());
      }
    }

    multiplicity.fill(HIST("fRhoVsTrkMult"), multTrack, rho_m);
    multiplicity.fill(HIST("fLeadingTrackPt"), flPt);
    multiplicity.fill(HIST("fTrackMult"), multTrack);
    if (hasValidV0) {
      multiplicity.fill(HIST("fMultFv0"), sumAmpFV0);
      multiplicity.fill(HIST("fTrackMultVsV0A"), multTrack, sumAmpFV0);
    }
    if (hasValidT0) {
      multiplicity.fill(HIST("fTrackMultVsT0A"), multTrack, sumAmpFT0);
    }
    // Check whether this is a high trk multiplicity event
    if (multTrack >= selHTrkMult) {
      keepEvent[kHighTrackMult] = true; // accepted HM events
      multiplicity.fill(HIST("fTrackMultSelected"), multTrack);
    }
    // Check whether this event has a leading track candidate
    if (flPt >= selPtTrig) {
      multiplicity.fill(HIST("fLeadingTrackPtSelected"), flPt); // track pT which passed the cut
      keepEvent[kLeadingPtTrack] = true;
    }
    // Check whether this is a high FV0 multiplicity event
    if (sumAmpFV0 >= selHMFv0) {
      keepEvent[kHighMultFv0] = true; // accepted HM events based on FV0
      multiplicity.fill(HIST("fMultFv0Selected"), sumAmpFV0);
    }

    tags(keepEvent[kLeadingPtTrack], keepEvent[kHighMultFv0], keepEvent[kHighTrackMult]);

    if (!keepEvent[kLeadingPtTrack] && !keepEvent[kHighMultFv0] && !keepEvent[kHighTrackMult]) {
      multiplicity.fill(HIST("fProcessedEvents"), 1);
    } else {
      for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
        if (keepEvent[iTrigger]) {
          multiplicity.fill(HIST("fProcessedEvents"), iTrigger + 2);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<multFilter>(cfg)};
}
