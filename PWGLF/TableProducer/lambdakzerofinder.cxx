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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "DataFormatsParameters/GRPObject.h"
#include <CCDB/BasicCCDBManager.h>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;

namespace o2::aod
{
namespace v0goodpostracks
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace v0goodpostracks
DECLARE_SOA_TABLE(V0GoodPosTracks, "AOD", "V0GOODPOSTRACKS", o2::soa::Index<>, v0goodpostracks::GoodTrackId, v0goodpostracks::CollisionId, v0goodpostracks::DCAXY);
namespace v0goodnegtracks
{
DECLARE_SOA_INDEX_COLUMN_FULL(GoodTrack, goodTrack, int, Tracks, "_GoodTrack");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
} // namespace v0goodnegtracks
DECLARE_SOA_TABLE(V0GoodNegTracks, "AOD", "V0GOODNEGTRACKS", o2::soa::Index<>, v0goodnegtracks::GoodTrackId, v0goodnegtracks::CollisionId, v0goodnegtracks::DCAXY);
} // namespace o2::aod

struct lambdakzeroprefilter {
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> tpcrefit{"tpcrefit", 1, "demand TPC refit"};

  Produces<aod::V0GoodPosTracks> v0GoodPosTracks;
  Produces<aod::V0GoodNegTracks> v0GoodNegTracks;

  // still exhibiting issues? To be checked
  // Partition<soa::Join<aod::FullTracks, aod::TracksDCA>> goodPosTracks = aod::track::signed1Pt > 0.0f && aod::track::dcaXY > dcapostopv;
  // Partition<soa::Join<aod::FullTracks, aod::TracksDCA>> goodNegTracks = aod::track::signed1Pt < 0.0f && aod::track::dcaXY < -dcanegtopv;

  void process(aod::Collision const& collision,
               soa::Join<aod::FullTracks, aod::TracksDCA> const& tracks)
  {
    for (auto& t0 : tracks) {
      if (tpcrefit) {
        if (!(t0.trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
      }
      if (t0.tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      if (t0.signed1Pt() > 0.0f) {
        if (fabs(t0.dcaXY()) < dcapostopv) {
          continue;
        }
        v0GoodPosTracks(t0.globalIndex(), t0.collisionId(), t0.dcaXY());
      }
      if (t0.signed1Pt() < 0.0f) {
        if (fabs(t0.dcaXY()) < dcanegtopv) {
          continue;
        }
        v0GoodNegTracks(t0.globalIndex(), t0.collisionId(), -t0.dcaXY());
      }
    }
  }
};

struct lambdakzerofinder {
  Produces<aod::StoredV0Datas> v0data;
  Produces<aod::V0s> v0;
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

  void init(InitContext& context)
  {
    // using namespace analysis::lambdakzerofinder;

    ccdb->setURL("https://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  float getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = grpo->getNominalL3Field();
    return output;
  }

  void process(aod::Collision const& collision, soa::Join<aod::FullTracks, aod::TracksCov> const& tracks,
               aod::V0GoodPosTracks const& ptracks, aod::V0GoodNegTracks const& ntracks, aod::BCsWithTimestamps const&)
  {

    float d_bz;
    if (d_bz_input < -990) {
      // Fetch magnetic field from ccdb for current collision
      d_bz = getMagneticField(collision.bc_as<aod::BCsWithTimestamps>().timestamp());
    } else {
      d_bz = d_bz_input;
    }

    // Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz);
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);

    Long_t lNCand = 0;

    for (auto& t0id : ptracks) { // FIXME: turn into combination(...)
      auto t0 = t0id.goodTrack_as<soa::Join<aod::FullTracks, aod::TracksCov>>();
      auto Track1 = getTrackParCov(t0);
      for (auto& t1id : ntracks) {
        auto t1 = t1id.goodTrack_as<soa::Join<aod::FullTracks, aod::TracksCov>>();
        auto Track2 = getTrackParCov(t1);

        // Try to progate to dca
        int nCand = fitter.process(Track1, Track2);
        if (nCand == 0) {
          continue;
        }
        const auto& vtx = fitter.getPCACandidate();

        // Fiducial: min radius
        auto thisv0radius = TMath::Sqrt(TMath::Power(vtx[0], 2) + TMath::Power(vtx[1], 2));
        if (thisv0radius < v0radius) {
          continue;
        }

        // DCA V0 daughters
        auto thisdcav0dau = fitter.getChi2AtPCACandidate();
        if (thisdcav0dau > dcav0dau) {
          continue;
        }

        std::array<float, 3> pos = {0.};
        std::array<float, 3> pvec0;
        std::array<float, 3> pvec1;
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
        }
        fitter.getTrack(0).getPxPyPzGlo(pvec0);
        fitter.getTrack(1).getPxPyPzGlo(pvec1);

        auto thisv0cospa = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()},
                                          array{vtx[0], vtx[1], vtx[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
        if (thisv0cospa < v0cospa) {
          continue;
        }

        lNCand++;
        v0(t0.collisionId(), t0.globalIndex(), t1.globalIndex());
        v0data(t0.globalIndex(), t1.globalIndex(), t0.collisionId(), 0,
               fitter.getTrack(0).getX(), fitter.getTrack(1).getX(),
               pos[0], pos[1], pos[2],
               pvec0[0], pvec0[1], pvec0[2],
               pvec1[0], pvec1[1], pvec1[2],
               fitter.getChi2AtPCACandidate(),
               t0id.dcaXY(), t1id.dcaXY());
        v0datalink(v0data.lastIndex());
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
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("hV0Radius"), v0.v0radius());
        registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
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
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }

    Long_t lNCand = 0;
    for (auto& v0 : fullV0s) {
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        registry.fill(HIST("hV0Radius"), v0.v0radius());
        registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
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

/// Extends the v0data table with expression columns
struct lambdakzeroinitializer {
  Spawns<aod::V0Datas> v0datas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeroprefilter>(cfgc, TaskName{"lf-lambdakzeroprefilter"}),
    adaptAnalysisTask<lambdakzerofinder>(cfgc, TaskName{"lf-lambdakzerofinder"}),
    adaptAnalysisTask<lambdakzerofinderQa>(cfgc, TaskName{"lf-lambdakzerofinderQA"}),
    adaptAnalysisTask<lambdakzeroinitializer>(cfgc, TaskName{"lf-lambdakzeroinitializer"})};
}
