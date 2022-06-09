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
  Configurable<float> selHMMFTMult{"selHMMFTMult", 24., "Minimum MFT multilicity  threshold"};
  Configurable<float> selHTrkMult{"selHTrkMult", 24., "Minimum global trk multiplicity threshold"};
  Configurable<float> selHITSTrkMult{"selHITSTrkMult", 135., "Minimum ITS trk multiplicity threshold"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 1.5f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};
  // rho cuts
  Configurable<int> cfgRhoNetaBins{"cfgRhoNetaBins", 5, "number of bins in eta"};
  Configurable<int> cfgRhoNphiBins{"cfgRhoNphiBins", 8, "number of bins in phi"};
  Configurable<float> cfgRhoMinEtaCut{"cfgRhoMinEtaCut", 2.2f, "min eta for fv0 cells"};
  Configurable<float> cfgRhoMaxEtaCut{"cfgRhoMaxEtaCut", 5.1f, "max eta for fv0 cells"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  void init(o2::framework::InitContext&)
  {

    // QA event level
    multiplicity.add("fCollZpos", "Vtx_z", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    // QA FV0
    multiplicity.add("fEtaPhiFv0", "eta vs phi", HistType::kTH2F, {{8, 0.0, 2 * M_PI, "#phi (rad)"}, {5, 2.2, 5.1, "#eta"}});

    // QA global tracks
    multiplicity.add("hdNdetaGlobal", "dNdeta", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiGlobal", "Phi", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});

    // QA ITS tracks
    multiplicity.add("hdNdetaITS", "dNdeta (ITStracks)", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiITS", "Phi (ITStracks)", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});

    // QA MFT tracks
    multiplicity.add("hdNdetaMFT", "dNdeta (MFTtracks)", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiMFT", "Phi (MFTtracks)", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});

    // QA leading track filter
    multiplicity.add("fLeadingTrackPt", "Lead trk pT", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fLeadingTrackPtSelected", "sel lead trk pT", HistType::kTH1F, {{150, 0., +150., "trk #it{p}_{T} (GeV/#it{c})"}});

    // QA hm FV0 filter
    multiplicity.add("fMultFv0", "FV0 amp", HistType::kTH1F, {{800, -0.5, +39999.5, "FV0 amplitude"}});
    multiplicity.add("fMultFv0Selected", "sel FV0 amp", HistType::kTH1F, {{800, -0.5, +39999.5, "FV0 amplitude"}});
    multiplicity.add("fTrackMultVsV0A", "trk mult vs FV0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {800, -0.5, +39999.5, "sum AmpFV0"}});
    multiplicity.add("fITSTrackMultVsV0A", "its trk mult vs FV0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "its trk mult (|#eta|<1.5)"}, {800, -0.5, +39999.5, "sum AmpFV0"}});

    // QA hm MFT filter
    multiplicity.add("fMFTTrackMult", "MFT trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "MFT trk mult (-3.6 < #eta < -2.5)"}});
    multiplicity.add("fMFTTrackMultSelected", "sel MFT trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "MFT trk mult (-3.6 < #eta < -2.5)"}});
    multiplicity.add("fTrackMultVsMFT", "trk mult vs MFT", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {200, -0.5, +199.5, "MFT trk mult (-3.6 < #eta < -2.5)"}});
    multiplicity.add("fITSTrackMultVsMFT", "its trk mult vs FV0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "its trk mult (|#eta|<1.5)"}, {200, -0.5, +199.5, "MFT trk mult (-3.6 < #eta < -2.5)"}});

    // QA hm global track filter
    multiplicity.add("fTrackMult", "trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}});
    multiplicity.add("fTrackMultSelected", "sel trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}});

    // QA its track filter
    multiplicity.add("fITSTrackMult", "ITS trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<1.5)"}});
    multiplicity.add("fITSTrackMultSelected", "sel ITS trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "itstrk mult (|#eta|<1.5)"}});

    // T0A
    multiplicity.add("fTrackMultVsT0A", "trk mult vs FT0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {200, -0.5, +999.5, "sum AmpFT0"}});

    // QA flatenicity
    multiplicity.add("fRhoVsTrkMult", "rhoN vs trk mult", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {1000, -0.5, +9.5, "#rho"}});
    multiplicity.add("fMpTVsTrackMult", "Av pT vs trk mult", HistType::kTProfile, {{200, -0.5, +199.5, "trk mult"}});
    multiplicity.add("fMpTVsTrackMultFlat", "Av pT vs trk mult (flat)", HistType::kTProfile, {{200, -0.5, +199.5, "trk mult"}});
    multiplicity.add("fMpTVsTrackMultJetty", "Av pT vs trk mult (jetty)", HistType::kTProfile, {{200, -0.5, +199.5, "trk mult"}});

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
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::pt > cfgTrkLowPtCut);

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
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
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>::iterator const& collision, TrackCandidates const& tracks, aod::MFTTracks const& mfttracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {

    bool keepEvent[kNtriggersMM]{false};

    multiplicity.fill(HIST("fProcessedEvents"), 0);
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    // global observables
    int multTrack = 0;
    int multITSTrack = 0;
    int multMFTTrack = 0;
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
        int channelv0phi = getFV0IndexPhi(channelv0);
        etav0 = cfgRhoMaxEtaCut - (detaFV0 / 2.0) * (2.0 * ringindex + 1);
        if (channelv0 < innerFV0) {
          phiv0 = (2.0 * (channelv0phi - 8 * ringindex) + 1) * M_PI / (8.0);
        } else {
          phiv0 = ((2.0 * channelv0phi) + 1 - 64.0) * 2.0 * M_PI / (32.0);
        }
        multiplicity.fill(HIST("fEtaPhiFv0"), phiv0, etav0);
        for (int i_eta = 0; i_eta < cfgRhoNetaBins; ++i_eta) {
          for (int i_phi = 0; i_phi < cfgRhoNphiBins; ++i_phi) {
            if (etav0 >= EtaBins[i_eta] && etav0 < EtaBins[i_eta + 1] &&
                phiv0 >= PhiBins[i_phi] && phiv0 < PhiBins[i_phi + 1]) {
              if (channelv0 < innerFV0) {
                RhoLattice[i_eta][i_phi] = fv0.amplitude()[ich];
              } else {
                RhoLattice[i_eta][i_phi] += fv0.amplitude()[ich] / 2.0; // two channels per bin
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
    // ITStracks
    for (auto& track : tracks) {
      if (track.detectorMap() && track.hasITS() != (uint8_t)0 && track.itsChi2NCl() < 36.0) {
        continue;
      }
      multITSTrack++;
      multiplicity.fill(HIST("hdNdetaITS"), track.eta());
      multiplicity.fill(HIST("hPhiITS"), track.phi());
    }
    // MFTtracks
    for (auto& track : mfttracks) {
      float eta = track.eta();
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (eta > -2.5 || eta < -3.6) { // the MFT eta coverage
        continue;
      }
      multMFTTrack++;
      multiplicity.fill(HIST("hdNdetaMFT"), track.eta());
      multiplicity.fill(HIST("hPhiMFT"), phi);
    }
    // Globaltracks
    for (auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
      multTrack++;
      multiplicity.fill(HIST("hdNdetaGlobal"), track.eta());
      multiplicity.fill(HIST("hPhiGlobal"), track.phi());
      if (flPt < track.pt()) {
        flPt = track.pt();
      }
    }
    // filling tprofiles
    for (auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
      multiplicity.fill(HIST("fMpTVsTrackMult"), multTrack, track.pt());
      if (rho_m < 0.8 && collision.has_foundFV0()) {
        multiplicity.fill(HIST("fMpTVsTrackMultFlat"), multTrack, track.pt());
      }
      if (rho_m > 2.0 && collision.has_foundFV0()) {
        multiplicity.fill(HIST("fMpTVsTrackMultJetty"), multTrack, track.pt());
      }
    }

    multiplicity.fill(HIST("fRhoVsTrkMult"), multTrack, rho_m);
    multiplicity.fill(HIST("fLeadingTrackPt"), flPt);
    multiplicity.fill(HIST("fTrackMult"), multTrack);
    multiplicity.fill(HIST("fITSTrackMult"), multITSTrack);
    multiplicity.fill(HIST("fMFTTrackMult"), multMFTTrack);
    multiplicity.fill(HIST("fTrackMultVsMFT"), multTrack, multMFTTrack);
    multiplicity.fill(HIST("fITSTrackMultVsMFT"), multITSTrack, multMFTTrack);

    if (hasValidV0) {
      multiplicity.fill(HIST("fMultFv0"), sumAmpFV0);
      multiplicity.fill(HIST("fTrackMultVsV0A"), multTrack, sumAmpFV0);
      multiplicity.fill(HIST("fITSTrackMultVsV0A"), multITSTrack, sumAmpFV0);
    }
    if (hasValidT0) {
      multiplicity.fill(HIST("fTrackMultVsT0A"), multTrack, sumAmpFT0);
    }
    // Check whether this is a high trk multiplicity event
    if (multTrack >= selHTrkMult) {
      keepEvent[kHighTrackMult] = true; // accepted HM events
      multiplicity.fill(HIST("fTrackMultSelected"), multTrack);
    }
    if (multITSTrack >= selHITSTrkMult) {
      multiplicity.fill(HIST("fITSTrackMultSelected"), multITSTrack);
    }
    if (multMFTTrack >= selHMMFTMult) {
      multiplicity.fill(HIST("fMFTTrackMultSelected"), multMFTTrack);
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
