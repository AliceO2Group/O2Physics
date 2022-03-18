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
//
// \file multFilter.cxx
// \task for selection of high multiplicity events
//
// \author Antonio Ortiz <antonio.ortiz.velasquez@cern.ch>, ICN-UNAM
//
// O2 includes

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

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

static const std::vector<std::string> mmObjectsNames{"kHtPt", "kHmFv0", "kHmTrk", "kHmTrkTns", "kHmTrk&kHmTrkTns"};

struct multFilter {
  enum { kLeadingPtTrack = 0,
         kHighMultFv0,
         kHighTrackMult,
         kHighTrackMultTrans,
         kHighTrackMultOverlap,
         kNtriggersMM };

  // event selection cuts
  Configurable<float> selPtTrig{"selPtTrig", 5., "Minimum track pT leading threshold"};
  Configurable<float> selHMFv0{"selHMFv0", 24000., "Minimum FV0-amplitude  threshold"};
  Configurable<float> selHTrkMult{"selHTrkMult", 24., "Minimum charged particle multiplicity threshold"};
  Configurable<float> selHTrkMultTrans{"selHTrkMultTrans", 6., "Minimum charged particle multiplicity (in transverse region) threshold"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};
  Configurable<int> cfgNetaBins{"cfgNetaBins", 4, "number of bins in eta"};
  Configurable<int> cfgNphiBins{"cfgNphiBins", 10, "number of bins in phi"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  void init(o2::framework::InitContext&)
  {

    multiplicity.add("fCollZpos", "Vtx_z", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});
    multiplicity.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
    multiplicity.add("fLeadingTrackPt", "Lead trk pT", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fMultFv0", "FV0 amp", HistType::kTH1F, {{800, -0.5, +39999.5, "FV0 amplitude"}});
    multiplicity.add("fTrackMult", "trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}});
    multiplicity.add("fRhoNVsTrkMult", "rhoN vs trk mult", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {1000, -0.5, +9.5, "#rho"}});
    multiplicity.add("fRhoPVsTrkMult", "rhoP vs trk mult", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {1000, -0.5, +9.5, "#rho"}});
    multiplicity.add("fTrackMultVsV0A", "trk mult vs FV0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {800, -0.5, +39999.5, "sum AmpFV0"}});
    multiplicity.add("fTrackMultVsT0A", "trk mult vs FT0 amp", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {200, -0.5, +999.5, "sum AmpFT0"}});
    multiplicity.add("fTrackMultTrans", "trk mult (TS)", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (TS)"}});
    multiplicity.add("fTrackMultVsTrans", "trk mult vs trk mut (TS)", HistType::kTH2F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}, {200, -0.5, +199.5, "trk mult (TS)"}});
    multiplicity.add("fLeadingTrackPtSelected", "sel lead trk pT", HistType::kTH1F, {{150, 0., +150., "trk #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fMultFv0Selected", "sel FV0 amp", HistType::kTH1F, {{800, -0.5, +39999.5, "FV0 amplitude"}});
    multiplicity.add("fTrackMultSelected", "sel trk mult", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (|#eta|<0.8)"}});
    multiplicity.add("fTrackMultTransSelected", "sel trk mult (TS)", HistType::kTH1F, {{200, -0.5, +199.5, "trk mult (TS)"}});
    multiplicity.add("fDeltaPhiTrans", "delta phi distribution in TS", HistType::kTH1F, {{64, -M_PI / 2.0, +3.0 * M_PI / 2.0, "#Delta#varphi (rad)"}});

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
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxCut;
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
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TrackCandidates const& tracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {

    bool keepEvent[kNtriggersMM]{false};

    multiplicity.fill(HIST("fProcessedEvents"), 0);
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    // V0A signal
    float sumAmpFT0 = 0;
    float sumAmpFV0 = 0;
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      if (collision.has_foundFT0()) {
        auto ft0 = collision.foundFT0();
        for (auto amplitude : ft0.amplitudeA()) {
          sumAmpFT0 += amplitude;
        }
        for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
          sumAmpFV0 += fv0.amplitude()[ich];
        }
      }
    }

    // global observables
    int multTrack = 0;
    float flPt = 0; // leading pT
    float flPhi = 0;
    double rho_mMpT = -1;
    double rho_mNch = -1;

    int minMult = 3;
    const int nEta = cfgNetaBins;
    const int nPhi = cfgNphiBins;
    double EtaBins[nEta + 1];
    double deltaEta = (2.0 * cfgTrkEtaCut) / (1.0 * cfgNetaBins);
    for (int i_eta = 0; i_eta < cfgNetaBins + 1; ++i_eta) {
      EtaBins[i_eta] = 0;
      if (i_eta < cfgNetaBins) {
        EtaBins[i_eta] = i_eta * deltaEta - 1.0 * cfgTrkEtaCut;
      } else {
        EtaBins[i_eta] = 1.0 * cfgTrkEtaCut;
      }
    }
    double PhiBins[nPhi + 1];
    double deltaPhi = (2.0 * M_PI) / (1.0 * cfgNphiBins);
    for (int i_phi = 0; i_phi < cfgNphiBins + 1; ++i_phi) {
      PhiBins[i_phi] = 0;
      if (i_phi < cfgNphiBins) {
        PhiBins[i_phi] = i_phi * deltaPhi - 1.0 * M_PI;
      } else {
        PhiBins[i_phi] = 1.0 * M_PI;
      }
    }
    int NchLattice[nEta][nPhi];
    for (int i_eta = 0; i_eta < cfgNetaBins; ++i_eta) {
      for (int i_phi = 0; i_phi < cfgNphiBins; ++i_phi) {
        NchLattice[i_eta][i_phi] = 0;
      }
    }
    double MpTLattice[nEta][nPhi];
    for (int i_eta = 0; i_eta < cfgNetaBins; ++i_eta) {
      for (int i_phi = 0; i_phi < cfgNphiBins; ++i_phi) {
        MpTLattice[i_eta][i_phi] = 0;
      }
    }

    int nchtotal = 0;
    for (auto& track : tracks) {
      multTrack++;
      multiplicity.fill(HIST("hdNdeta"), track.eta());
      for (int i_eta = 0; i_eta < cfgNetaBins; ++i_eta) {
        for (int i_phi = 0; i_phi < cfgNphiBins; ++i_phi) {
          if (track.eta() >= EtaBins[i_eta] && track.eta() < EtaBins[i_eta + 1] &&
              track.phi() >= PhiBins[i_phi] && track.phi() < PhiBins[i_phi + 1]) {
            NchLattice[i_eta][i_phi]++;
            MpTLattice[i_eta][i_phi] += track.pt();
          }
        }
      }

      if (flPt < track.pt()) {
        flPt = track.pt();
        flPhi = track.phi();
      }
    }
    if (nchtotal < minMult) {
      rho_mMpT = -1;
      rho_mNch = -1;
    }
    for (int i_eta = 0; i_eta < cfgNetaBins; ++i_eta) {
      for (int i_phi = 0; i_phi < cfgNphiBins; ++i_phi) {
        if (NchLattice[i_eta][i_phi] > 0)
          MpTLattice[i_eta][i_phi] /= (1.0 * NchLattice[i_eta][i_phi]);
        else
          MpTLattice[i_eta][i_phi] = 0.0;
      }
    }

    double mNch = 0;
    double mMpT = 0;
    for (int i_eta = 0; i_eta < cfgNetaBins; ++i_eta) {
      for (int i_phi = 0; i_phi < cfgNphiBins; ++i_phi) {
        mMpT += MpTLattice[i_eta][i_phi];
        mNch += 1.0 * NchLattice[i_eta][i_phi];
      }
    }
    // average activity per cell
    mMpT /= (1.0 * cfgNetaBins * cfgNphiBins);
    mNch /= (1.0 * cfgNetaBins * cfgNphiBins);
    // get sigma
    double sNch_tmp = 0;
    double sMpT_tmp = 0;
    for (int i_eta = 0; i_eta < cfgNetaBins; ++i_eta) {
      for (int i_phi = 0; i_phi < cfgNphiBins; ++i_phi) {
        sMpT_tmp += TMath::Power(MpTLattice[i_eta][i_phi] - mMpT, 2);
        sNch_tmp += TMath::Power(1.0 * NchLattice[i_eta][i_phi] - mNch, 2);
      }
    }
    sMpT_tmp /= (1.0 * cfgNetaBins * cfgNphiBins);
    sNch_tmp /= (1.0 * cfgNetaBins * cfgNphiBins);
    double sNch = TMath::Sqrt(sNch_tmp);
    double sMpT = TMath::Sqrt(sMpT_tmp);
    rho_mMpT = sMpT / mMpT;
    rho_mNch = sNch / mNch;
    multiplicity.fill(HIST("fRhoNVsTrkMult"), multTrack, rho_mNch);
    multiplicity.fill(HIST("fRhoPVsTrkMult"), multTrack, rho_mMpT);
    multiplicity.fill(HIST("fLeadingTrackPt"), flPt);
    multiplicity.fill(HIST("fTrackMult"), multTrack);
    multiplicity.fill(HIST("fMultFv0"), sumAmpFV0);
    multiplicity.fill(HIST("fTrackMultVsV0A"), multTrack, sumAmpFV0);
    multiplicity.fill(HIST("fTrackMultVsT0A"), multTrack, sumAmpFT0);

    // Check whether this is a high multiplicity event
    if (multTrack >= selHTrkMult) {
      keepEvent[kHighTrackMult] = true; // accepted HM events
      multiplicity.fill(HIST("fTrackMultSelected"), multTrack);
    }

    // Check whether this event has a leading track candidate
    int multTrackTrans = 0;
    if (flPt >= selPtTrig) {
      // compute Nch transverse
      for (auto& track : tracks) { // start loop over tracks
        float DPhi = computeDeltaPhi(track.phi(), flPhi);
        if (!(TMath::Abs(DPhi) < M_PI / 3.0 || TMath::Abs(DPhi - M_PI) < M_PI / 3.0)) { // transverse side
          multTrackTrans++;
          multiplicity.fill(HIST("fDeltaPhiTrans"), DPhi);
        }
      }
      multiplicity.fill(HIST("fTrackMultVsTrans"), multTrack, multTrackTrans);
      multiplicity.fill(HIST("fTrackMultTrans"), multTrackTrans);
      if (multTrackTrans >= selHTrkMultTrans) {
        keepEvent[kHighTrackMultTrans] = true; // accepted HM events
        multiplicity.fill(HIST("fTrackMultTransSelected"), multTrackTrans);
        if (keepEvent[kHighTrackMult] && keepEvent[kHighTrackMultTrans]) {
          keepEvent[kHighTrackMultOverlap] = true;
        }
      }
      multiplicity.fill(HIST("fLeadingTrackPtSelected"), flPt); // track pT which passed the cut
      keepEvent[kLeadingPtTrack] = true;
    }

    if (sumAmpFV0 >= selHMFv0) {
      keepEvent[kHighMultFv0] = true; // accepted HM events based on FV0
      multiplicity.fill(HIST("fMultFv0Selected"), sumAmpFV0);
    }

    tags(keepEvent[kLeadingPtTrack], keepEvent[kHighMultFv0], keepEvent[kHighTrackMult], keepEvent[kHighTrackMultTrans], keepEvent[kHighTrackMultOverlap]);

    if (!keepEvent[kLeadingPtTrack] && !keepEvent[kHighMultFv0] && !keepEvent[kHighTrackMult] && !keepEvent[kHighTrackMultTrans]) {
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

  return WorkflowSpec{
    adaptAnalysisTask<multFilter>(cfg)};
}
