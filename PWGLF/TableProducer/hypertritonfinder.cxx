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
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
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

//use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullHe>;
using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels, aod::pidTPCFullPi, aod::pidTPCFullHe>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullHe>;
using FullTracksExtMCIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

using MyTracks = FullTracksExt;
using MyTracksIU = FullTracksExtIU;

inline float GetTPCNSigmaHe3(float p, float TPCSignal)
{
  float bg = p/2.80839;
  return  (TPCSignal - o2::tpc::BetheBlochAleph(bg, -9.973f, -18.5543f, 29.5704f, 2.02064f, -3.85076f)) / (TPCSignal*0.0812);
}

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

struct hypertritonprefilter {
  HistogramRegistry registry{
    "registry",
      {
        {"hCrossedRows", "hCrossedRows", {HistType::kTH1F, {{50, 0.0f, 200.0f}}}},
        {"hGoodTrackCount", "hGoodTrackCount",{HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
        {"hGoodPosTrackCount", "hGoodPosTrackCount",{HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
        {"hGoodNegTrackCount", "hGoodNegTrackCount",{HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      },
  };

  //change the dca cut for helium3
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> tpcrefit{"tpcrefit", 0, "demand TPC refit"};

  Produces<aod::V0GoodPosTracks> v0GoodPosTracks;
  Produces<aod::V0GoodNegTracks> v0GoodNegTracks;

  void process(aod::Collision const& collision,
      MyTracksIU const& tracks)
  {
    for (auto& t0 : tracks) {
      registry.fill(HIST("hGoodTrackCount"), 0.5);
      registry.fill(HIST("hCrossedRows"), t0.tpcNClsCrossedRows());
      if (tpcrefit) {
        if (!(t0.trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
      }
      registry.fill(HIST("hGoodTrackCount"), 1.5);
      if (t0.tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      registry.fill(HIST("hGoodTrackCount"), 2.5);
      if (t0.signed1Pt() > 0.0f) {
        /*if (fabs(t0.dcaXY()) < dcapostopv) {
          continue;
          }*/
        v0GoodPosTracks(t0.globalIndex(), t0.collisionId(), t0.dcaXY());
        registry.fill(HIST("hGoodPosTrackCount"), 0.5);
        registry.fill(HIST("hGoodTrackCount"), 3.5);
      }
      if (t0.signed1Pt() < 0.0f) {
        /*if (fabs(t0.dcaXY()) < dcanegtopv) {
          continue;
          }*/
        v0GoodNegTracks(t0.globalIndex(), t0.collisionId(), -t0.dcaXY());
        registry.fill(HIST("hGoodNegTrackCount"), 0.5);
        registry.fill(HIST("hGoodTrackCount"), 3.5);
      }
    }
  }
};

struct hypertritonfinder {

  Produces<aod::StoredV0Datas> v0data;
  Produces<aod::V0s> v0;
  Produces<aod::V0DataLink> v0datalink;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  Configurable<double> d_UseAbsDCA{"d_UseAbsDCA", kTRUE, "Use Abs DCAs"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  //Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};

  Configurable<int> useMatCorrType{"useMatCorrType", 2, "0: none, 1: TGeo, 2: LUT"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  HistogramRegistry registry{
    "registry",
      {
        {"hCandPerEvent", "hCandPerEvent", {HistType::kTH1F, {{1000, 0.0f, 1000.0f}}}},
        {"hV0CutCounter", "hV0CutCounter", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
      },
  };

  int mRunNumber;
  float d_bz;
  float maxSnp;  //max sine phi for propagation
  float maxStep; //max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  //could be changed later
    maxStep = 2.00f; //could be changed later

    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(1, "DiffCol");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(2, "hasSV");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(3, "hasSV2");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(4, "Dcav0Dau");
    registry.get<TH1>(HIST("hV0CutCounter"))->GetXaxis()->SetBinLabel(5, "CosPA");

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);

    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }
    fitter.setMatCorrType(matCorr);
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

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }
  //------------------------------------------------------------------

  void process(aod::Collision const& collision, MyTracksIU const& tracks,
      aod::V0GoodPosTracks const& ptracks, aod::V0GoodNegTracks const& ntracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    Long_t lNCand = 0;

    for (auto& t0id : ptracks) { // FIXME: turn into combination(...)
      for (auto& t1id : ntracks) {

        if (t0id.collisionId() != t1id.collisionId()) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 0.5);

        auto t0 = t0id.goodTrack_as<MyTracksIU>();
        auto t1 = t1id.goodTrack_as<MyTracksIU>();
        auto Track1 = getTrackParCov(t0);
        auto Track2 = getTrackParCov(t1);
        auto pTrack = getTrackParCov(t0);
        auto nTrack = getTrackParCov(t1);

        // Try to progate to dca
        int nCand = fitter.process(Track1, Track2);
        if (nCand == 0) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 1.5);

        //------------------copy from lamdakzerobuilder---------------------
        double finalXpos = fitter.getTrack(0).getX();
        double finalXneg = fitter.getTrack(1).getX();

        // Rotate to desired alpha
        pTrack.rotateParam(fitter.getTrack(0).getAlpha());
        nTrack.rotateParam(fitter.getTrack(1).getAlpha());

        // Retry closer to minimum with material corrections
        o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
        if (useMatCorrType == 1)
          matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
        if (useMatCorrType == 2)
          matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

        o2::base::Propagator::Instance()->propagateToX(pTrack, finalXpos, d_bz, maxSnp, maxStep, matCorr);
        o2::base::Propagator::Instance()->propagateToX(nTrack, finalXneg, d_bz, maxSnp, maxStep, matCorr);

        nCand = fitter.process(pTrack, nTrack);
        if (nCand == 0) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 2.5);

        //------------------------------------------------------------------

        const auto& vtx = fitter.getPCACandidate();
        // Fiducial: min radius
        /*auto thisv0radius = TMath::Sqrt(TMath::Power(vtx[0], 2) + TMath::Power(vtx[1], 2));
          if (thisv0radius < v0radius) {
          continue;
          }*/

        // DCA V0 daughters
        auto thisdcav0dau = fitter.getChi2AtPCACandidate();
        if (thisdcav0dau > dcav0dau) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 3.5);

        std::array<float, 3> pos = {0.};
        std::array<float, 3> pvec0;
        std::array<float, 3> pvec1;
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
        }
        //fitter.getTrack(0).getPxPyPzGlo(pvec0);
        //fitter.getTrack(1).getPxPyPzGlo(pvec1);

        //------------------copy from lamdakzerobuilder---------------------
        pTrack.getPxPyPzGlo(pvec0);
        nTrack.getPxPyPzGlo(pvec1);
        //------------------------------------------------------------------
        /*uint32_t pTrackPID = t0.pidForTracking();
          uint32_t nTrackPID = t1.pidForTracking();
          int pTrackCharge = o2::track::pid_constants::sCharges[pTrackPID];
          int nTrackCharge = o2::track::pid_constants::sCharges[nTrackPID];
          for (int i=0; i<3; i++){
          pvec0[i] = pvec0[i] * pTrackCharge;
          pvec1[i] = pvec1[i] * nTrackCharge;
          }*/
        int pTrackCharge = 1, nTrackCharge = 1;
        if (TMath::Abs( GetTPCNSigmaHe3( 2*t0.p(), t0.tpcSignal()) ) < 5){
          pTrackCharge = 2;
        } 
        if (TMath::Abs( GetTPCNSigmaHe3( 2*t1.p(), t1.tpcSignal()) ) < 5){
          nTrackCharge = 2;
        } 
        for (int i=0; i<3; i++){
          pvec0[i] = pvec0[i] * pTrackCharge;
          pvec1[i] = pvec1[i] * nTrackCharge;
        }


        auto thisv0cospa = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()},
            array{vtx[0], vtx[1], vtx[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
        if (thisv0cospa < v0cospa) {
          continue;
        }
        registry.fill(HIST("hV0CutCounter"), 4.5);

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

/// Extends the v0data table with expression columns
struct hypertritoninitializer {
  Spawns<aod::V0Datas> v0datas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonprefilter>(cfgc),
    adaptAnalysisTask<hypertritonfinder>(cfgc),
    adaptAnalysisTask<hypertritoninitializer>(cfgc)};
}
