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
// V0 Finder task
// ==============
//
// This code loops over positive and negative tracks and finds
// valid V0 candidates from scratch using a certain set of
// minimum (configurable) selection criteria.
//
// It is different than the producer: the producer merely
// loops over an *existing* list of V0s (pos+neg track
// indices) and calculates the corresponding full V0 information
//
// In both cases, any analysis should loop over the "V0Data"
// table as that table contains all information.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessFinderTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;

// For dE/dx association in pre-selection
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA>;
using TracksExtraWithDCAnPID = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::TracksDCA>;

struct lambdakzeroprefilter {
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> dcaSingleTrackToPV{"dcaSingleTrackToPV", .1, "DCA single track to PV"};

  // dEdx pre-selection compatibility
  Configurable<float> ddEdxPreSelectionWindow{"ddEdxPreSelectionWindow", 7, "Nsigma window for dE/dx preselection"};

  Produces<aod::VFinderTracks> VFinderTracks;

  void processAll(aod::Collision const& /*collision*/,
                  soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA> const& tracks)
  {
    for (auto& t0 : tracks) {
      if (TMath::Abs(t0.dcaXY()) < dcaSingleTrackToPV)
        continue; // single track DCA to PV
      if (t0.tpcNClsCrossedRows() < mincrossedrows)
        continue; // crossed rows
      bool isPositive = false;
      if (t0.signed1Pt() > 0.0f)
        isPositive = true;
      VFinderTracks(t0.globalIndex(), isPositive, true, true, true);
    }
  }
  PROCESS_SWITCH(lambdakzeroprefilter, processAll, "Take all tracks, select only on crossed rows + TPC refit", false);

  void processWithdEdx(aod::Collision const& /*collision*/,
                       TracksExtraWithDCAnPID const& tracks)
  {
    for (auto& t0 : tracks) {
      if (TMath::Abs(t0.dcaXY()) < dcaSingleTrackToPV)
        continue; // single track DCA to PV
      if (t0.tpcNClsCrossedRows() < mincrossedrows)
        continue; // crossed rows
      bool compatiblePion = false;
      bool compatibleKaon = false;
      bool compatibleProton = false;
      if (TMath::Abs(t0.tpcNSigmaPi()) < ddEdxPreSelectionWindow)
        compatiblePion = true;
      if (TMath::Abs(t0.tpcNSigmaKa()) < ddEdxPreSelectionWindow)
        compatibleKaon = true;
      if (TMath::Abs(t0.tpcNSigmaPr()) < ddEdxPreSelectionWindow)
        compatibleProton = true;
      bool isPositive = false;
      if (t0.signed1Pt() > 0.0f)
        isPositive = true;
      VFinderTracks(t0.globalIndex(), isPositive, compatiblePion, compatibleKaon, compatibleProton);
    }
  }
  PROCESS_SWITCH(lambdakzeroprefilter, processWithdEdx, "Utilise TPC dE/dx pre-selections", true);
};

struct lambdakzerofinder {
  Produces<aod::V0Indices> v0indices;
  Produces<aod::V0CoresBase> v0cores;
  Produces<aod::V0TrackXs> v0trackXs;
  Produces<aod::V0s_001> v0;
  Produces<aod::V0DataLink> v0datalink;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry",
    {
      {"hCandPerEvent", "hCandPerEvent", {HistType::kTH1F, {{1000, 0.0f, 1000.0f}}}},
    },
  };

  // Configurables
  Configurable<double> d_UseAbsDCA{"d_UseAbsDCA", kTRUE, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> maxV0DCAtoPV{"maxV0DCAtoPV", 0.5, "maximum V0 DCA to PV"};

  // Configurables for selecting which particles to generate
  Configurable<bool> findK0Short{"findK0Short", true, "findK0Short"};
  Configurable<bool> findLambda{"findLambda", true, "findLambda"};
  Configurable<bool> findAntiLambda{"findAntiLambda", true, "findAntiLambda"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Partition<aod::VFinderTracks> pTracks = o2::aod::vFinderTrack::isPositive == true;
  Partition<aod::VFinderTracks> nTracks = o2::aod::vFinderTrack::isPositive == false;

  // Define o2 fitter, 2-prong
  o2::vertexing::DCAFitterN<2> fitter;
  int mRunNumber;
  float d_bz;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    ccdb->setURL("https://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);
  }

  float getDCAtoPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }

  template <class TTrack, class TCollisions>
  int buildV0Candidate(TTrack const& t1, TTrack const& t2, TCollisions const& collisions)
  {
    auto Track1 = getTrackParCov(t1);
    auto Track2 = getTrackParCov(t2);

    // Try to progate to dca
    int nCand = fitter.process(Track1, Track2);
    if (nCand == 0) {
      return 0;
    }
    const auto& vtx = fitter.getPCACandidate();

    // Fiducial: min radius
    auto thisv0radius = TMath::Sqrt(TMath::Power(vtx[0], 2) + TMath::Power(vtx[1], 2));
    if (thisv0radius < v0radius) {
      return 0;
    }

    // DCA V0 daughters
    auto thisdcav0dau = fitter.getChi2AtPCACandidate();
    if (thisdcav0dau > dcav0dau) {
      return 0;
    }

    std::array<float, 3> pos = {0.};
    std::array<float, 3> pvec0;
    std::array<float, 3> pvec1;
    for (int i = 0; i < 3; i++) {
      pos[i] = vtx[i];
    }
    fitter.getTrack(0).getPxPyPzGlo(pvec0);
    fitter.getTrack(1).getPxPyPzGlo(pvec1);

    // Attempt collision association on pure geometrical basis
    // FIXME this can of course be far better
    float smallestDCA = 1e+3;
    float cosPA = -1;
    int collisionIndex = -1;
    //   float getDCAtoPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ){
    for (auto const& collision : collisions) {
      float thisDCA = TMath::Abs(getDCAtoPV(vtx[0], vtx[1], vtx[2], pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2], collision.posX(), collision.posY(), collision.posZ()));
      if (thisDCA < smallestDCA) {
        collisionIndex = collision.globalIndex();
        smallestDCA = thisDCA;
        cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, array{vtx[0], vtx[1], vtx[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
      }
    }
    if (smallestDCA > maxV0DCAtoPV)
      return 0; // unassociated

    v0(collisionIndex, t1.globalIndex(), t2.globalIndex());

    // populates the various tables for analysis
    v0indices(t1.globalIndex(), t2.globalIndex(), collisionIndex, 0);
    v0trackXs(fitter.getTrack(0).getX(), fitter.getTrack(1).getX());
    v0cores(pos[0], pos[1], pos[2],
            pvec0[0], pvec0[1], pvec0[2],
            pvec1[0], pvec1[1], pvec1[2],
            TMath::Sqrt(fitter.getChi2AtPCACandidate()),
            t1.dcaXY(), t2.dcaXY(), cosPA, smallestDCA, 1);
    v0datalink(v0cores.lastIndex(), -1);
    return 1;
  }

  void process(aod::Collisions const& collisions, FullTracksExtIU const& /*tracks*/,
               aod::VFinderTracks const& /*v0findertracks*/, aod::BCsWithTimestamps const&)
  {
    auto firstcollision = collisions.begin();
    auto bc = firstcollision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    Long_t lNCand = 0;

    for (auto& pTrack : pTracks) { // FIXME: turn into combination(...)
      for (auto& nTrack : nTracks) {
        // Check compatibility with certain hypotheses and desired building
        bool keepCandidate = false;
        if (pTrack.compatiblePi() && nTrack.compatiblePi() && findK0Short)
          keepCandidate = true;
        if (pTrack.compatiblePr() && nTrack.compatiblePi() && findLambda)
          keepCandidate = true;
        if (pTrack.compatiblePi() && nTrack.compatiblePr() && findAntiLambda)
          keepCandidate = true;
        if (!keepCandidate)
          continue;

        auto t1 = pTrack.track_as<FullTracksExtIU>();
        auto t2 = nTrack.track_as<FullTracksExtIU>();

        lNCand += buildV0Candidate(t1, t2, collisions);
      }
    }
    registry.fill(HIST("hCandPerEvent"), lNCand);
  }
};

struct lambdakzerofinderQa {
  // Basic checks
  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.998, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", .6, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};

  HistogramRegistry registry{
    "registry",
    {
      {"hCandPerEvent", "hCandPerEvent", {HistType::kTH1F, {{1000, 0.0f, 1000.0f}}}},

      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},

      {"h3dMassK0Short", "h3dMassK0Short", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 0.450f, 0.550f}}}},
      {"h3dMassLambda", "h3dMassLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
      {"h3dMassAntiLambda", "h3dMassAntiLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
    },
  };

  void init(InitContext const&)
  {
    if (doprocessRun3 && doprocessRun2) {
      LOGF(fatal, "processRun3 and processRun2 are both set to true; try again with only one of them set to true");
    }
    if (!doprocessRun3 && !doprocessRun2) {
      LOGF(fatal, "processRun3 nor processRun2 are both set to false; try again with only one of them set to false");
    }
  }

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  /// Connect to V0Data: newly indexed, note: V0Datas table incompatible with standard V0 table!
  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As>::iterator const& collision,
                   soa::Filtered<aod::V0Datas> const& fullV0s)
  {
    if (!collision.sel8()) {
      return;
    }

    Long_t lNCand = 0;
    for (auto& v0 : fullV0s) {
      if (v0.v0radius() > v0radius && v0.v0cosPA() > v0cospa) {
        registry.fill(HIST("hV0Radius"), v0.v0radius());
        registry.fill(HIST("hV0CosPA"), v0.v0cosPA());
        registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
        registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
        registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());

        if (TMath::Abs(v0.yLambda()) < 0.5) {
          registry.fill(HIST("h3dMassLambda"), collision.centFV0A(), v0.pt(), v0.mLambda());
          registry.fill(HIST("h3dMassAntiLambda"), collision.centFV0A(), v0.pt(), v0.mAntiLambda());
        }
        if (TMath::Abs(v0.yK0Short()) < 0.5) {
          registry.fill(HIST("h3dMassK0Short"), collision.centFV0A(), v0.pt(), v0.mK0Short());
        }
        lNCand++;
      }
    }
    registry.fill(HIST("hCandPerEvent"), lNCand);
  }
  PROCESS_SWITCH(lambdakzerofinderQa, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision,
                   soa::Filtered<aod::V0Datas> const& fullV0s)
  {
    if (!collision.alias_bit(kINT7)) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }

    Long_t lNCand = 0;
    for (auto& v0 : fullV0s) {
      if (v0.v0radius() > v0radius && v0.v0cosPA() > v0cospa) {
        registry.fill(HIST("hV0Radius"), v0.v0radius());
        registry.fill(HIST("hV0CosPA"), v0.v0cosPA());
        registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
        registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
        registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());

        if (TMath::Abs(v0.yLambda()) < 0.5) {
          registry.fill(HIST("h3dMassLambda"), collision.centRun2V0M(), v0.pt(), v0.mLambda());
          registry.fill(HIST("h3dMassAntiLambda"), collision.centRun2V0M(), v0.pt(), v0.mAntiLambda());
        }
        if (TMath::Abs(v0.yK0Short()) < 0.5) {
          registry.fill(HIST("h3dMassK0Short"), collision.centRun2V0M(), v0.pt(), v0.mK0Short());
        }
        lNCand++;
      }
    }
    registry.fill(HIST("hCandPerEvent"), lNCand);
  }
  PROCESS_SWITCH(lambdakzerofinderQa, processRun2, "Process Run 2 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeroprefilter>(cfgc, TaskName{"lf-lambdakzeroprefilter"}),
    adaptAnalysisTask<lambdakzerofinder>(cfgc, TaskName{"lf-lambdakzerofinder"}),
    adaptAnalysisTask<lambdakzerofinderQa>(cfgc, TaskName{"lf-lambdakzerofinderQA"})};
}
