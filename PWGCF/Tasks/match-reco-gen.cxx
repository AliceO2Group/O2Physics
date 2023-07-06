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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TParameter.h>
#include <TProfile3D.h>
#include <TROOT.h>

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
struct CheckGeneratorLevelVsDetectorLevel {
  Configurable<int> cfgTrackType{"trktype", 1, "Type of selected tracks: 0 = no selection, 1 = global tracks FB96"};
  Configurable<std::string> cfgCentMultEstimator{"centmultestimator", "V0M", "Centrality/multiplicity estimator detector:  V0M, NOCM: none. Default V0M"};
  Configurable<std::string> cfgSystem{"syst", "PbPb", "System: pp, PbPb, Pbp, pPb, XeXe, ppRun3. Default PbPb"};
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, datanoevsel, MC, FastMC, OnTheFlyMC. Default data"};
  Configurable<std::string> cfgTriggSel{"triggsel", "MB", "Trigger selection: MB, None. Default MB"};
  Configurable<o2::analysis::DptDptBinningCuts> cfgBinning{"binning",
                                                           {28, -7.0, 7.0, 18, 0.2, 2.0, 16, -0.8, 0.8, 72, 0.5},
                                                           "triplets - nbins, min, max - for z_vtx, pT, eta and phi, binning plus bin fraction of phi origin shift"};
  Configurable<o2::analysis::CheckRangeCfg> cfgTraceDCAOutliers{"trackdcaoutliers", {false, 0.0, 0.0}, "Track the generator level DCAxy outliers: false/true, low dcaxy, up dcaxy. Default {false,0.0,0.0}"};
  Configurable<float> cfgTraceOutOfSpeciesParticles{"trackoutparticles", false, "Track the particles which are not e,mu,pi,K,p: false/true. Default false"};
  Configurable<int> cfgRecoIdMethod{"recoidmethod", 0, "Method for identifying reconstructed tracks: 0 PID, 1 mcparticle. Default 0"};
  Configurable<o2::analysis::TrackSelectionCfg> cfgTrackSelection{"tracksel", {false, false, 0, 70, 0.8, 2.4, 3.2}, "Track selection: {useit: true/false, ongen: true/false, tpccls, tpcxrws, tpcxrfc, dcaxy, dcaz}. Default {false,0.70.0.8,2.4,3.2}"};
  Configurable<bool> cfgTraceCollId0{"tracecollid0", false, "Trace particles in collisions id 0. Default false"};
  Configurable<bool> cfgTrackMultiRec{"trackmultirec", false, "Track muli-reconstructed particles: true, false. Default false"};
  Configurable<bool> cfgTrackCollAssoc{"trackcollassoc", false, "Track collision id association, track-mcparticle-mccollision vs. track-collision-mccollision: true, false. Default false"};

  HistogramRegistry histos{"RecoGenHistograms", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  typedef enum { kBEFORE = 0,
                 kAFTER } beforeafterselection;
  typedef enum { kPOSITIVE = 0,
                 kNEGATIVE } colllabelsign;
  enum { kMATCH = 0,
         kDONTMATCH };

  void init(InitContext const&)
  {
    using namespace o2::analysis::recogenmap;
    using namespace o2::analysis::dptdptfilter;

    /* update with the configurable values */
    /* the binning */
    ptbins = cfgBinning->mPTbins;
    ptlow = cfgBinning->mPTmin;
    ptup = cfgBinning->mPTmax;
    etabins = cfgBinning->mEtabins;
    etalow = cfgBinning->mEtamin;
    etaup = cfgBinning->mEtamax;
    zvtxbins = cfgBinning->mZVtxbins;
    zvtxlow = cfgBinning->mZVtxmin;
    zvtxup = cfgBinning->mZVtxmax;
    /* the track types and combinations */
    tracktype = cfgTrackType.value;
    initializeTrackSelection();
    /* the centrality/multiplicity estimation */
    fCentMultEstimator = getCentMultEstimator(cfgCentMultEstimator);
    /* the trigger selection */
    fTriggerSelection = getTriggerSelection(cfgTriggSel);
    traceDCAOutliers = cfgTraceDCAOutliers;
    traceOutOfSpeciesParticles = cfgTraceOutOfSpeciesParticles;
    recoIdMethod = cfgRecoIdMethod;
    if (cfgTrackSelection->mUseIt) {
      useOwnTrackSelection = true;
      if (cfgTrackSelection->mOnGen) {
        useOwnParticleSelection = true;
        particleMaxDCAxy = cfgTrackSelection->mDCAxy;
        particleMaxDCAZ = cfgTrackSelection->mDCAz;
      }
      ownTrackSelection.SetMinNClustersTPC(cfgTrackSelection->mTPCclusters);
      ownTrackSelection.SetMinNCrossedRowsTPC(cfgTrackSelection->mTPCxRows);
      ownTrackSelection.SetMinNCrossedRowsOverFindableClustersTPC(cfgTrackSelection->mTPCXRoFClusters);
      ownTrackSelection.SetMaxDcaXYPtDep(std::function<float(float)>{});
      ownTrackSelection.SetMaxDcaXY(cfgTrackSelection->mDCAxy);
      ownTrackSelection.SetMaxDcaZ(cfgTrackSelection->mDCAz);
      o2::aod::track::TrackTypeEnum ttype;
      switch (tracktype) {
        case 1:
          ttype = o2::aod::track::Run2Track;
          break;
        case 3:
          ttype = o2::aod::track::Track;
          break;
        default:
          ttype = o2::aod::track::Track;
          break;
      }
      ownTrackSelection.SetTrackType(ttype);
    } else {
      useOwnTrackSelection = false;
    }
    traceCollId0 = cfgTraceCollId0;

    /* if the system type is not known at this time, we have to put the initialization somewhere else */
    fSystem = getSystemType(cfgSystem);
    fDataType = getDataType(cfgDataType);
    fPDG = TDatabasePDG::Instance();

    AxisSpec deltaEta = {100, -2, 2, "#Delta#eta"};
    AxisSpec deltaPhi = {100, 0, constants::math::TwoPI, "#Delta#varphi (rad)"};
    AxisSpec deltaPt = {1000, 0, 4, "#Delta#it{p}_{T} (GeV/#it{c})"};
    AxisSpec mrectimes = {11, -0.5f, 10.5f, "##/particle"};
    AxisSpec detectors = {32, -0.5, 31.5, "Detectors"};
    std::vector<std::string> detectorlbls = {"", "ITS", "TPC", "ITS+TPC", "TRD", "ITS+TRD", "TPC+TRD", "ITS+TPC+TRD",
                                             "TOF", "ITS+TOF", "TPC+TOF", "ITS+TPC+TOF", "TRD+TOF", "ITS+TRD+TOF", "TPC+TRD+TOF", "ITS+TPC+TRD+TOF",
                                             "UNKN", "ITS+UNKN", "TPC+UNKN", "ITS+TPC+UNKN", "TRD+UNKN", "ITS+TRD+UNKN", "TPC+TRD+UNKN", "ITS+TPC+TRD+UNKN",
                                             "TOF+UNKN", "ITS+TOF+UNKN", "TPC+TOF+UNKN", "ITS+TPC+TOF+UNKN", "TRD+TOF+UNKN", "ITS+TRD+TOF+UNKN", "TPC+TRD+TOF+UNKN", "ITS+TPC+TRD+TOF+UNKN"};
    std::vector<std::string> matchlbs = {"match", "don't match"};

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
    for (unsigned int i = 0; i < detectorlbls.size(); ++i) {
      histos.get<TH1>(HIST("before/positivecolid/detectormap"))->GetXaxis()->SetBinLabel(i + 1, detectorlbls[i].c_str());
      histos.get<TH1>(HIST("before/positivecolid/detectormapmr"))->GetXaxis()->SetBinLabel(i + 1, detectorlbls[i].c_str());
    }
    for (unsigned int i = 0; i < matchlbs.size(); ++i) {
      histos.get<TH1>(HIST("before/positivecolid/matchcollid"))->GetXaxis()->SetBinLabel(i + 1, matchlbs[i].c_str());
      histos.get<TH1>(HIST("before/positivecolid/matchcollidmr"))->GetXaxis()->SetBinLabel(i + 1, matchlbs[i].c_str());
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

    static constexpr std::string_view dir[] = {"before/", "after/"};
    static constexpr std::string_view colldir[] = {"positivecolid/", "negativecolid/"};

    int nrec_poslabel = 0;
    int nrec_neglabel = 0;
    int nrec_poslabel_crosscoll = 0;

    for (int ixpart = 0; ixpart < mcParticles.size(); ++ixpart) {
      auto particle = mcParticles.iteratorAt(ixpart);
      /* multireconstructed tracks only for positive labels */
      int nrec = mclabelpos[collsign][ixpart].size();
      nrec_poslabel += mclabelpos[collsign][ixpart].size();
      nrec_neglabel += mclabelneg[collsign][ixpart].size();

      if (nrec > 1) {
        /* multireconstruction only from positive labels */
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("multirec"), nrec);

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
                nrec_poslabel_crosscoll++;
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
            histos.get<TH1>(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("pdgcodemr"))
              ->Fill(TString::Format("%d", particle.pdgCode()).Data(), 1.0);
          }
        }

        for (unsigned int i = 0; i < mclabelpos[collsign][ixpart].size(); ++i) {
          auto track1 = tracks.iteratorAt(mclabelpos[collsign][ixpart][i]);
          for (unsigned int j = i + 1; j < mclabelpos[collsign][ixpart].size(); ++j) {
            auto track2 = tracks.iteratorAt(mclabelpos[collsign][ixpart][j]);

            float deltaeta = track1.eta() - track2.eta();
            float deltaphi = track1.phi() - track2.phi();
            if (deltaphi < 0) {
              deltaphi += constants::math::TwoPI;
            }
            if (deltaphi > constants::math::TwoPI) {
              deltaphi -= constants::math::TwoPI;
            }
            float deltapt = (track1.pt() > track2.pt()) ? track1.pt() - track2.pt() : track2.pt() - track1.pt();

            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("mrDeltaEta"), deltaeta);
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("mrDeltaPhi"), deltaphi);
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("mrDeltaPt"), deltapt);
          }
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("recomreta"), track1.eta());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("recomrphi"), track1.phi());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("recomrpt"), track1.pt());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("detectormapmr"), track1.detectorMap());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcaxymr"), track1.dcaXY());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcazmr"), track1.dcaZ());
          if (track1.dcaXY() < 1.0) {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcaxymr"), track1.dcaXY());
          }
          if (track1.dcaZ() < 1.0) {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcazmr"), track1.dcaZ());
          }
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecomreta"), track1.eta(), particle.eta());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecomrphi"), track1.phi(), particle.phi());
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecomrpt"), track1.pt(), particle.pt());
          if (particle.mcCollisionId() != colls.iteratorAt(track1.collisionId()).mcCollisionId()) {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollidmr"), kDONTMATCH + 0.5f);
          } else {
            histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollidmr"), kMATCH + 0.5f);
          }
        }
      } else if (nrec > 0) {
        auto track = tracks.iteratorAt(mclabelpos[collsign][ixpart][0]);
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecoeta"), track.eta(), particle.eta());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecophi"), track.phi(), particle.phi());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("genrecopt"), track.pt(), particle.pt());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("detectormap"), track.detectorMap());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcaxy"), track.dcaXY());
        histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("dcaz"), track.dcaZ());
        if (track.dcaXY() < 1.0) {
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcaxy"), track.dcaXY());
        }
        if (track.dcaZ() < 1.0) {
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("finedcaz"), track.dcaZ());
        }
        if (particle.mcCollisionId() != colls.iteratorAt(track.collisionId()).mcCollisionId()) {
          if ((ba == kAFTER) && (collsign == kPOSITIVE) && cfgTrackCollAssoc) {
            LOGF(info, "Particle with index %d and pdg code %d assigned to MC collision %d, pT: %f, phi: %f, eta: %f",
                 particle.globalIndex(), particle.pdgCode(), particle.mcCollisionId(), particle.pt(), particle.phi(), particle.eta());
            LOGF(info, "        with status %d and flags %0X and", particle.statusCode(), particle.flags());
            LOGF(info, "        associated to track with index %d and label %d assigned to collision %d, with associated MC collision %d",
                 track.globalIndex(), ixpart, track.collisionId(), colls.iteratorAt(track.collisionId()).mcCollisionId());
          }
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollid"), kDONTMATCH + 0.5f);
        } else {
          histos.fill(HIST(dir[ba]) + HIST(colldir[collsign]) + HIST("matchcollid"), kMATCH + 0.5f);
        }
      }
    }

    if (collsign == kPOSITIVE) {
      LOGF(info, "Reconstructed tracks (%s) with positive collision ID: %d with positive label, %d with negative label, %d with cross collision",
           ba == kAFTER ? "after" : "before", nrec_poslabel, nrec_neglabel, nrec_poslabel_crosscoll);
    } else {
      LOGF(info, "Reconstructed tracks (%s) with negative collision ID: %d with positive label, %d with negative label",
           ba == kAFTER ? "after" : "before", nrec_poslabel, nrec_neglabel);
    }
  }

  template <typename TracksObject, typename CollisionsObject>
  void processMapChecksBeforeCuts(TracksObject const& tracks, CollisionsObject const& collisions, aod::McParticles const& mcParticles)
  {
    using namespace o2::analysis::recogenmap;
    using namespace o2::analysis::dptdptfilter;

    for (int i = 0; i < 2; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = tracks.size();
    size_t ngen = 0;

    for (auto& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = (pdgpart->Charge() >= 3) ? 1.0 : ((pdgpart->Charge() <= -3) ? -1.0 : 0.0);
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    // For the time being we are only interested in the information based on the reconstructed tracks
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed tracks", ngen, nreco);

    for (auto& track : tracks) {
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

    for (int i = 0; i < 2; ++i) {
      mclabelpos[i].clear();
      mclabelneg[i].clear();
      mclabelpos[i].resize(mcParticles.size());
      mclabelneg[i].resize(mcParticles.size());
    }

    size_t nreco = 0;
    size_t ngen = 0;

    for (auto& part : mcParticles) {
      auto pdgpart = fPDG->GetParticle(part.pdgCode());
      if (pdgpart != nullptr) {
        float charge = (pdgpart->Charge() >= 3) ? 1.0 : ((pdgpart->Charge() <= -3) ? -1.0 : 0.0);
        if (charge != 0.0) {
          ngen++;
        }
      }
    }

    // Let's go through the reco-gen mapping to detect multi-reconstructed particles
    for (auto& track : tracks) {
      int64_t recix = track.globalIndex();
      int32_t label = track.mcParticleId();
      if (!(label < 0)) {
        if (!(track.collisionId() < 0)) {
          typename CollisionsObject::iterator coll = collisions.iteratorAt(track.collisionId());
          float centormult = -100.0f;
          if (IsEvtSelected(coll, centormult)) {
            /* TODO: AcceptTrack does not consider PID */
            int pid = AcceptTrack(track);
            if ((pid == 0) || (pid == 1)) {
              /* the track has been accepted */
              nreco++;
              LOGF(MATCHRECGENLOGTRACKS, "Accepted track with global Id %d and collision Id %d has label %d associated to MC collision %d", recix, track.collisionId(), label, track.template mcParticle_as<aod::McParticles>().mcCollisionId());
              mclabelpos[kPOSITIVE][label].push_back(recix);
            } else {
              if (pid > 1) {
                LOGF(fatal, "Task not prepared for PID");
              }
            }
          }
        }
      }
    }
    LOGF(info, "New dataframe (DF) with %d generated charged particles and %d reconstructed accepted tracks", ngen, nreco);

    collectData<kAFTER, kPOSITIVE>(tracks, mcParticles, collisions);
  }

  void processMapChecksWithCent(soa::Join<aod::FullTracks, aod::TracksDCA, aod::McTrackLabels> const& tracks,
                                soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels> const& collisions,
                                aod::McParticles const& mcParticles)
  {
    processMapChecksBeforeCuts(tracks, collisions, mcParticles);
    processMapChecksAfterCuts(tracks, collisions, mcParticles);
  }
  PROCESS_SWITCH(CheckGeneratorLevelVsDetectorLevel, processMapChecksWithCent, "Process detector <=> generator levels with centrality/multiplicity information", false);

  void processMapChecksWithoutCent(soa::Join<aod::FullTracks, aod::TracksDCA, aod::McTrackLabels> const& tracks,
                                   soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels> const& collisions,
                                   aod::McParticles const& mcParticles)
  {
    processMapChecksBeforeCuts(tracks, collisions, mcParticles);
    processMapChecksAfterCuts(tracks, collisions, mcParticles);
  }
  PROCESS_SWITCH(CheckGeneratorLevelVsDetectorLevel, processMapChecksWithoutCent, "Process detector <=> generator levels without centrality/multiplicity information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<CheckGeneratorLevelVsDetectorLevel>(cfgc)};
  return workflow;
}
