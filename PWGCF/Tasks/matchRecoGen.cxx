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

/// \file matchRecoGen.cxx
/// \brief basic check for the matching between generator level and detector level
/// \author victor.gonzalez.sebastian@gmail.com

#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptDptFilter.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TParameter.h>
#include <TProfile3D.h>
#include <TROOT.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define MATCHRECGENLOGCOLLISIONS debug
#define MATCHRECGENLOGTRACKS debug

namespace o2::analysis::recogenmap
{
std::vector<std::vector<int64_t>> mclabelpos[2];
std::vector<std::vector<int64_t>> mclabelneg[2];
} // namespace o2::analysis::recogenmap

/// \brief Checks the correspondence generator level <=> detector level
struct MatchRecoGen {
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"cfgTraceDCAOutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"cfgTraceOutOfSpeciesParticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<bool> cfgTraceCollId0{"cfgTraceCollId0", false, "Trace particles in collisions id 0. Default false"};
  Configurable<bool> cfgTrackMultiRec{"cfgTrackMultiRec", false, "Track muli-reconstructed particles: true, false. Default false"};
  Configurable<bool> cfgTrackCollAssoc{"cfgTrackCollAssoc", false, "Track collision id association, track-mcparticle-mccollision vs. track-collision-mccollision: true, false. Default false"};

  HistogramRegistry histos{"RecoGenHistograms", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Service<o2::framework::O2DatabasePDG> fPDG;
  typedef enum { kBEFORE = 0,
                 kAFTER } beforeafterselection;
  typedef enum { kPOSITIVE = 0,
                 kNEGATIVE,
                 kNOOFCOLLSIGNS } colllabelsign;
  enum { kMATCH = 0,
         kDONTMATCH };

  void init(InitContext const&)
  {
    using namespace o2::analysis::recogenmap;
    using namespace o2::analysis::dptdptfilter;

    /* update with the configurable values */
    traceDCAOutliers = cfgTraceDCAOutliers;
    traceOutOfSpeciesParticles = cfgTraceOutOfSpeciesParticles;
    traceCollId0 = cfgTraceCollId0;

    AxisSpec deltaEta = {100, -2, 2, "#Delta#eta"};
    AxisSpec deltaPhi = {100, 0, constants::math::TwoPI, "#Delta#varphi (rad)"};
    AxisSpec deltaPt = {1000, 0, 4, "#Delta#it{p}_{T} (GeV/#it{c})"};
    AxisSpec mrectimes = {11, -0.5f, 10.5f, "##/particle"};
    AxisSpec detectors = {32, -0.5, 31.5, "Detectors"};
    std::vector<std::string> detectorLabels = {"", "ITS", "TPC", "ITS+TPC", "TRD", "ITS+TRD", "TPC+TRD", "ITS+TPC+TRD",
                                               "TOF", "ITS+TOF", "TPC+TOF", "ITS+TPC+TOF", "TRD+TOF", "ITS+TRD+TOF", "TPC+TRD+TOF", "ITS+TPC+TRD+TOF",
                                               "UNKN", "ITS+UNKN", "TPC+UNKN", "ITS+TPC+UNKN", "TRD+UNKN", "ITS+TRD+UNKN", "TPC+TRD+UNKN", "ITS+TPC+TRD+UNKN",
                                               "TOF+UNKN", "ITS+TOF+UNKN", "TPC+TOF+UNKN", "ITS+TPC+TOF+UNKN", "TRD+TOF+UNKN", "ITS+TRD+TOF+UNKN", "TPC+TRD+TOF+UNKN", "ITS+TPC+TRD+TOF+UNKN"};
    std::vector<std::string> matchLabels = {"match", "don't match"};

    histos.add("before/positivecolid/mrDeltaEta", "#Delta#eta multirec tracks", kTH1F, {deltaEta});
    histos.add("before/positivecolid/mrDeltaPhi", "#Delta#varphi multirec tracks", kTH1F, {deltaPhi});
    histos.add("before/positivecolid/mrDeltaPt", "#Delta#it{p}_{T} multirec tracks", kTH1F, {deltaPt});
    histos.add("before/positivecolid/multirec", "Multiple reconstruction", kTH1F, {mrectimes});
    histos.add("before/positivecolid/genrecoeta", "#eta Generated vs reconstructed", kTH2F, {{100, -1.0, 1.0, "#eta reco"}, {100, -1.0, 1.0, "#eta gen"}});
    histos.add("before/positivecolid/genrecophi", "#varphi Generated vs reconstructed", kTH2F, {{100, 0, constants::math::TwoPI, "#varphi (rad) reco"}, {100, 0, constants::math::TwoPI, "#varphi (rad) gen"}});
    histos.add("before/positivecolid/genrecopt", "#it{p}_{T} Generated vs reconstructed", kTH2F, {{1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) reco"}, {1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) gen"}});
    histos.add("before/positivecolid/detectormap", "Active detectors", kTH1F, {detectors});
    histos.add("before/positivecolid/matchcollid", "particle MC coll Id <=> track coll MC coll Id", kTH1F, {{2, 0.0, 2.0}});
    histos.add("before/positivecolid/genrecomreta", "#eta Generated vs reconstructed (mr)", kTH2F, {{100, -1.0, 1.0, "#eta reco"}, {100, -1.0, 1.0, "#eta gen"}});
    histos.add("before/positivecolid/genrecomrphi", "#varphi Generated vs reconstructed (mr)", kTH2F, {{100, 0, constants::math::TwoPI, "#varphi (rad) reco"}, {100, 0, constants::math::TwoPI, "#varphi (rad) gen"}});
    histos.add("before/positivecolid/genrecomrpt", "#it{p}_{T} Generated vs reconstructed (mr)", kTH2F, {{1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) reco"}, {1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c}) gen"}});
    histos.add("before/positivecolid/recomreta", "#eta Reconstructed (mr)", kTH1F, {{100, -1.0, 1.0, "#eta"}});
    histos.add("before/positivecolid/recomrphi", "#varphi Reconstructed (mr)", kTH1F, {{100, 0, constants::math::TwoPI, "#varphi (rad)"}});
    histos.add("before/positivecolid/recomrpt", "#it{p}_{T} Reconstructed (mr)", kTH1F, {{1000, 0, 10.0, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("before/positivecolid/detectormapmr", "Active detectors (mr)", kTH1F, {detectors});
    histos.add("before/positivecolid/matchcollidmr", "particle MC coll Id <=> track coll MC coll Id (mr)", kTH1F, {{2, 0.0, 2.0}});
    histos.add("before/positivecolid/dcaxy", "DCA_{xy} Reconstructed", kTH1F, {{1000, -4.0, 4.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/dcaz", "DCA_{z} Reconstructed", kTH1F, {{1000, -4.0, 4.0, "DCA_{z} (cm)"}});
    histos.add("before/positivecolid/finedcaxy", "DCA_{xy} Reconstructed", kTH1F, {{2000, -1.0, 1.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/finedcaz", "DCA_{z} Reconstructed", kTH1F, {{2000, -1.0, 1.0, "DCA_{z} (cm)"}});
    histos.add("before/positivecolid/dcaxymr", "DCA_{xy} Reconstructed (mr)", kTH1F, {{1000, -4.0, 4.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/dcazmr", "DCA_{z} Reconstructed (mr)", kTH1F, {{1000, -4.0, 4.0, "DCA_{z} (cm)"}});
    histos.add("before/positivecolid/finedcaxymr", "DCA_{xy} Reconstructed (mr)", kTH1F, {{2000, -1.0, 1.0, "DCA_{xy} (cm)"}});
    histos.add("before/positivecolid/finedcazmr", "DCA_{z} Reconstructed (mr)", kTH1F, {{2000, -1.0, 1.0, "DCA_{z} (cm)"}});
    for (unsigned int i = 0; i < detectorLabels.size(); ++i) {
      histos.get<TH1>(HIST("before/positivecolid/detectormap"))->GetXaxis()->SetBinLabel(i + 1, detectorLabels[i].c_str());
      histos.get<TH1>(HIST("before/positivecolid/detectormapmr"))->GetXaxis()->SetBinLabel(i + 1, detectorLabels[i].c_str());
    }
    for (unsigned int i = 0; i < matchLabels.size(); ++i) {
      histos.get<TH1>(HIST("before/positivecolid/matchcollid"))->GetXaxis()->SetBinLabel(i + 1, matchLabels[i].c_str());
      histos.get<TH1>(HIST("before/positivecolid/matchcollidmr"))->GetXaxis()->SetBinLabel(i + 1, matchLabels[i].c_str());
    }

    /* clone the set for the other cases */
    histos.addClone("before/positivecolid/", "after/positivecolid/");
    histos.addClone("before/positivecolid/", "before/negativecolid/");
    histos.addClone("before/positivecolid/", "after/negativecolid/");
    histos.add("after/positivecolid/pdgcodemr", "PDG code x-collision multi-reconstructed", kTH1F, {{100, 0.5, 100.5, "PDG code"}});
  }

  template <beforeafterselection ba, colllabelsign collsign, typename TracskListObject, typename ParticlesListObject, typename CollisionsListObject>
  void collectData(TracskListObject const& tracks, ParticlesListObject const& mcParticles, CollisionsListObject const& colls)
  {
    using namespace o2::analysis::recogenmap;

    static constexpr std::string_view Dir[] = {"before/", "after/"};
    static constexpr std::string_view Colldir[] = {"positivecolid/", "negativecolid/"};

    int nRecPosLabel = 0;
    int nRecNegLabel = 0;
    int nRecPosLabelCrossColl = 0;

    for (int ixpart = 0; ixpart < mcParticles.size(); ++ixpart) {
      auto particle = mcParticles.iteratorAt(ixpart);
      /* multireconstructed tracks only for positive labels */
      int nrec = mclabelpos[collsign][ixpart].size();
      nRecPosLabel += mclabelpos[collsign][ixpart].size();
      nRecNegLabel += mclabelneg[collsign][ixpart].size();

      if (nrec > 1) {
        /* multireconstruction only from positive labels */
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("multirec"), nrec);

        if (collsign == kPOSITIVE) {
          /* check the cross collision reconstruction */
          bool crosscollfound = false;
          for (unsigned int i = 0; (i < mclabelpos[collsign][ixpart].size()) && !crosscollfound;
               ++i) {
            for (unsigned int j = i + 1;
                 (j < mclabelpos[collsign][ixpart].size()) && !crosscollfound;
                 ++j) {
              auto track1 = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
              auto track2 = tracks.iteratorAt(mclabelpos[collsign][ixpart][j]);

              if (track1.collisionId() != track2.collisionId()) {
                nRecPosLabelCrossColl++;
                crosscollfound = true;
              }
            }
          }
          if (crosscollfound && (ba == kAFTER)) {
            if (cfgTrackMultiRec) {
              LOGF(info,
                   "BEGIN multi-reconstructed: "
                   "==================================================================");
              LOGF(info,
                   "Particle with index %d and pdg code %d assigned to MC collision %d, pT: "
                   "%f, phi: %f, eta: %f",
                   particle.globalIndex(),
                   particle.pdgCode(),
                   particle.mcCollisionId(),
                   particle.pt(),
                   particle.phi(),
                   particle.eta());
              LOGF(info,
                   "With status %d and flags %0X and multi-reconstructed as: "
                   "==================================",
                   particle.statusCode(),
                   particle.flags());
              for (unsigned int i = 0; i < mclabelpos[collsign][ixpart].size(); ++i) {
                auto track = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
                auto coll = colls.iteratorAt(track.collisionId());
                LOGF(info,
                     "Track with index %d and label %d assigned to collision %d, with "
                     "associated MC collision %d",
                     track.globalIndex(),
                     ixpart,
                     track.collisionId(),
                     coll.mcCollisionId());
              }
              LOGF(info,
                   "END multi-reconstructed:   "
                   "==================================================================");
            }
            histos.get<TH1>(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("pdgcodemr"))->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
          }
        }

        for (unsigned int i = 0; i < mclabelpos[collsign][ixpart].size(); ++i) {
          auto track1 = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
          for (unsigned int j = i + 1; j < mclabelpos[collsign][ixpart].size(); ++j) {
            auto track2 = tracks.iteratorAt(mclabelpos[collsign][ixpart][j]);

            float deltaeta = track1.eta() - track2.eta();
            float deltaphi = track1.phi() - track2.phi();
            deltaphi = RecoDecay::constrainAngle(deltaphi, 0.0f);
            float deltapt = (track1.pt() > track2.pt()) ? track1.pt() - track2.pt() : track2.pt() - track1.pt();

            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("mrDeltaEta"), deltaeta);
            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("mrDeltaPhi"), deltaphi);
            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("mrDeltaPt"), deltapt);
          }
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("recomreta"), track1.eta());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("recomrphi"), track1.phi());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("recomrpt"), track1.pt());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("detectormapmr"), track1.detectorMap());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("dcaxymr"), track1.dcaXY());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("dcazmr"), track1.dcaZ());
          if (std::fabs(track1.dcaXY()) < 1.0) {
            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("finedcaxymr"), track1.dcaXY());
          }
          if (std::fabs(track1.dcaZ()) < 1.0) {
            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("finedcazmr"), track1.dcaZ());
          }
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("genrecomreta"), track1.eta(), particle.eta());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("genrecomrphi"), track1.phi(), particle.phi());
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("genrecomrpt"), track1.pt(), particle.pt());
          if (particle.mcCollisionId() != colls.iteratorAt(track1.collisionId()).mcCollisionId()) {
            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("matchcollidmr"), static_cast<float>(kDONTMATCH) + 0.5f);
          } else {
            histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("matchcollidmr"), static_cast<float>(kMATCH) + 0.5f);
          }
        }
      } else if (nrec > 0) {
        auto track = tracks.iteratorAt(mclabelpos[collsign][ixpart][0]);
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("genrecoeta"), track.eta(), particle.eta());
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("genrecophi"), track.phi(), particle.phi());
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("genrecopt"), track.pt(), particle.pt());
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("detectormap"), track.detectorMap());
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("dcaxy"), track.dcaXY());
        histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("dcaz"), track.dcaZ());
        if (std::fabs(track.dcaXY()) < 1.0) {
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("finedcaxy"), track.dcaXY());
        }
        if (std::fabs(track.dcaZ()) < 1.0) {
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("finedcaz"), track.dcaZ());
        }
        if (particle.mcCollisionId() != colls.iteratorAt(track.collisionId()).mcCollisionId()) {
          if ((ba == kAFTER) && (collsign == kPOSITIVE) && cfgTrackCollAssoc) {
            LOGF(info, "Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
                 particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
            LOGF(info, "        with status %d and flags %0X and", particle.statusCode(), particle.flags());
            LOGF(info, "        associated to track with index %d and label %d assigned to collision %d, with associated MC collision %d",
                 track.globalIndex(), ixpart, track.collisionId(), colls.iteratorAt(track.collisionId()).mcCollisionId());
          }
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("matchcollid"), static_cast<float>(kDONTMATCH) + 0.5f);
        } else {
          histos.fill(HIST(Dir[ba]) + HIST(Colldir[collsign]) + HIST("matchcollid"), static_cast<float>(kMATCH) + 0.5f);
        }
      }
    }

    if (collsign == kPOSITIVE) {
      LOGF(info, "Reconstructed tracks (%s) with positive collision ID: %d with positive label, %d with negative label, %d with cross collision",
           ba == kAFTER ? "after" : "before", nRecPosLabel, nRecNegLabel, nRecPosLabelCrossColl);
    } else {
      LOGF(info, "Reconstructed tracks (%s) with negative collision ID: %d with positive label, %d with negative label",
           ba == kAFTER ? "after" : "before", nRecPosLabel, nRecNegLabel);
    }
  }

  template <typename TracksObject, typename CollisionsObject>
  void processMapChecksBeforeCuts(TracksObject const& tracks, CollisionsObject const& collisions, aod::McParticles const& mcParticles)
  {
    using namespace o2::analysis::recogenmap;
    using namespace o2::analysis::dptdptfilter;

    for (int i = 0; i < kNOOFCOLLSIGNS; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = tracks.size();
    size_t ngen = 0;

    for (auto const& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = getCharge(pdgpart->Charge());
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    // For the time being we are only interested in the information based on the reconstructed tracks
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed tracks", ngen, nreco);

    for (auto const& track : tracks) {
      int64_t recix = track.globalIndex();
      int32_t label = track.mcParticleId();

      LOGF(MATCHRECGENLOGTRACKS, "Track with global Id %d and collision Id %d has label %d associated to MC collision %d", recix, track.collisionId(), label, track.template mcParticle_as<aod::McParticles>().mcCollisionId());
      if (track.collisionId() < 0) {
        if (label >= 0) {
          mclabelpos[kNEGATIVE][label].push_back(recix);
        } else {
          mclabelneg[kNEGATIVE][-label].push_back(recix);
        }
      } else {
        if (label >= 0) {
          mclabelpos[kPOSITIVE][label].push_back(recix);
        } else {
          mclabelneg[kPOSITIVE][-label].push_back(recix);
        }
      }
    }

    collectData<kBEFORE, kPOSITIVE>(tracks, mcParticles, collisions);
    collectData<kBEFORE, kNEGATIVE>(tracks, mcParticles, collisions);
  }

  template <typename TracksObject, typename CollisionsObject>
  void processMapChecksAfterCuts(TracksObject const& tracks, CollisionsObject const& collisions, aod::McParticles const& mcParticles)
  {
    using namespace o2::analysis::recogenmap;
    using namespace o2::analysis::dptdptfilter;

    for (int i = 0; i < kNOOFCOLLSIGNS; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = 0;
    size_t ngen = 0;

    for (auto const& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = getCharge(pdgpart->Charge());
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    for (auto const& track : tracks) {
      int64_t recix = track.globalIndex();
      int32_t label = track.mcParticleId();
      if (!(label < 0)) {
        if (!(track.collisionId() < 0)) {
          typename CollisionsObject::iterator coll = collisions.iteratorAt(track.collisionId());
          if (coll.collisionaccepted() == uint8_t(true)) {
            /* TODO: AcceptTrack does not consider PID */
            if (!(track.trackacceptedid() < 0)) {
              /* the track has been accepted */
              nreco++;
              LOGF(MATCHRECGENLOGTRACKS, "Accepted track with global Id %d and collision Id %d has label %d associated to MC collision %d", recix, track.collisionId(), label, track.template mcParticle_as<aod::McParticles>().mcCollisionId());
              mclabelpos[kPOSITIVE][label].push_back(recix);
            }
          }
        }
      }
    }
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed accepted tracks", ngen, nreco);

    collectData<kAFTER, kPOSITIVE>(tracks, mcParticles, collisions);
  }

  void processMapChecks(soa::Join<aod::FullTracks, aod::TracksDCA, aod::DptDptCFTracksInfo, aod::McTrackLabels> const& tracks,
                        soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo, aod::McCollisionLabels> const& collisions,
                        aod::McParticles const& mcParticles)
  {
    processMapChecksBeforeCuts(tracks, collisions, mcParticles);
    processMapChecksAfterCuts(tracks, collisions, mcParticles);
  }
  PROCESS_SWITCH(MatchRecoGen, processMapChecks, "Process detector <=> generator levels", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<MatchRecoGen>(cfgc)};
  return workflow;
}
