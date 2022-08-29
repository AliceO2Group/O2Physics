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
// V0 builder task
// ================
//
// This task loops over an *existing* list of V0s (neg/pos track
// indices) and calculates the corresponding full V0 information
//
// Any analysis should loop over the "V0Data"
// table as that table contains all information
//
// WARNING: adding filters to the builder IS NOT
// equivalent to re-running the finders. This will only
// ever produce *tighter* selection sections. It is your
// responsibility to check if, by setting a loose filter
// setting, you are going into a region in which no
// candidates exist because the original indices were generated
// using tigher selections.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <CCDB/BasicCCDBManager.h>

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"
#include "PWGHF/Utils/UtilsDebugLcK0Sp.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using FullTracksExtMCIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

// in case requested
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

//#define MY_DEBUG
#ifdef MY_DEBUG
using MyTracks = FullTracksExtMC;
using MyTracksIU = FullTracksExtMCIU;
#define MY_DEBUG_MSG(condition, cmd) \
  if (condition) {                   \
    cmd;                             \
  }
#else
using MyTracks = FullTracksExt;
using MyTracksIU = FullTracksExtIU;
#define MY_DEBUG_MSG(condition, cmd)
#endif

// Builder task: rebuilds V0s
// The prefilter part skims the list of good V0s to re-reconstruct so that
// CPU is saved in case there are specific selections that are to be done
// Note: more configurables, more options to be added as needed
struct lambdakzeroBuilder {
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<int> isRun2{"isRun2", 0, "if Run2: demand TPC refit"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;

  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> v_labelK0Spos{"v_labelK0Spos", {729, 2866, 4754}, "labels of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> v_labelK0Sneg{"v_labelK0Sneg", {730, 2867, 4755}, "labels of K0S negative daughters, for debug"};
#endif

  Produces<aod::StoredV0Datas> v0data;
  Produces<aod::V0DataLink> v0dataLink;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      {"hV0Criteria", "hV0Criteria", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
    },
  };

  // Configurables
  // Configurable<int> d_UseAbsDCA{"d_UseAbsDCA", 1, "Use Abs DCAs"}; uncomment this once we want to use the weighted DCA
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<int> rejDiffCollTracks{"rejDiffCollTracks", 0, "rejDiffCollTracks"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> v_labelK0Spos{"v_labelK0Spos", {729, 2866, 4754}, "labels of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> v_labelK0Sneg{"v_labelK0Sneg", {730, 2867, 4755}, "labels of K0S positive daughters, for debug"};
#endif

  void init(InitContext& context)
  {
    // using namespace analysis::lambdakzerobuilder;
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    if (doprocessRun3 && doprocessRun2) {
      LOGF(fatal, "processRun3 and processRun2 are both set to true; try again with only one of them set to true");
    }
    if (!doprocessRun3 && !doprocessRun2) {
      LOGF(fatal, "processRun3 nor processRun2 are both set to false; try again with only one of them set to false");
    }
  }
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!grpo) {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
    }
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", run3grp_timestamp, d_bz);
      } else {
        d_bz = d_bz_input;
      }
    } else {
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    mRunNumber = bc.runNumber();
  }

  template <class TCascTracksTo>
  void buildLambdaKZeroTable(aod::Collision const& collision, aod::V0s const& V0s, Bool_t lRun3 = kTRUE)
  {
    // Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitter;
    fitter.setBz(d_bz);
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true); // use d_UseAbsDCA once we want to use the weighted DCA

    registry.fill(HIST("hEventCounter"), 0.5);

    for (auto& V0 : V0s) {

      // Track preselection part
      auto posTrackCast = V0.posTrack_as<TCascTracksTo>();
      auto negTrackCast = V0.negTrack_as<TCascTracksTo>();

#ifdef MY_DEBUG
      auto labelPos = posTrackCast.mcParticleId();
      auto labelNeg = negTrackCast.mcParticleId();
      bool isK0SfromLc = isK0SfromLcFunc(labelPos, labelNeg, v_labelK0Spos, v_labelK0Sneg);
#endif
      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "V0 builder: found K0S from Lc, posTrack --> " << labelPos << ", negTrack --> " << labelNeg);

      // value 0.5: any considered V0
      registry.fill(HIST("hV0Criteria"), 0.5);
      if (isRun2) {
        if (!(posTrackCast.trackType() & o2::aod::track::TPCrefit) && !lRun3) {
          MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack " << labelPos << " has no TPC refit");
          v0dataLink(-1);
          continue; // TPC refit
        }
        if (!(negTrackCast.trackType() & o2::aod::track::TPCrefit) && !lRun3) {
          MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "negTrack " << labelNeg << " has no TPC refit");
          v0dataLink(-1);
          continue; // TPC refit
        }
      }
      // Passes TPC refit
      registry.fill(HIST("hV0Criteria"), 1.5);
      if (posTrackCast.tpcNClsCrossedRows() < mincrossedrows) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack " << labelPos << " has " << posTrackCast.tpcNClsCrossedRows() << " crossed rows, cut at " << mincrossedrows);
        v0dataLink(-1);
        continue;
      }
      if (negTrackCast.tpcNClsCrossedRows() < mincrossedrows) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "negTrack " << labelNeg << " has " << negTrackCast.tpcNClsCrossedRows() << " crossed rows, cut at " << mincrossedrows);
        v0dataLink(-1);
        continue;
      }
      // passes crossed rows
      registry.fill(HIST("hV0Criteria"), 2.5);
      if (fabs(posTrackCast.dcaXY()) < dcapostopv) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack " << labelPos << " has dcaXY " << posTrackCast.dcaXY() << " , cut at " << dcanegtopv);
        v0dataLink(-1);
        continue;
      }
      if (fabs(negTrackCast.dcaXY()) < dcanegtopv) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "negTrack " << labelNeg << " has dcaXY " << negTrackCast.dcaXY() << " , cut at " << dcanegtopv);
        v0dataLink(-1);
        continue;
      }
      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "Filling good indices: posTrack --> " << labelPos << ", negTrack --> " << labelNeg);
      // passes DCAxy
      registry.fill(HIST("hV0Criteria"), 3.5);

      // Candidate building part
      std::array<float, 3> pos = {0.};
      std::array<float, 3> pvec0 = {0.};
      std::array<float, 3> pvec1 = {0.};

#ifdef MY_DEBUG
      auto labelPos = posTrackCast.mcParticleId();
      auto labelNeg = negTrackCast.mcParticleId();
      bool isK0SfromLc = isK0SfromLcFunc(labelPos, labelNeg, v_labelK0Spos, v_labelK0Sneg);
#endif

      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "labelPos = " << labelPos << ", labelNeg = " << labelNeg);

      auto pTrack = getTrackParCov(posTrackCast);
      auto nTrack = getTrackParCov(negTrackCast);

      // Require collision-ID
      if (posTrackCast.collisionId() != negTrackCast.collisionId() && rejDiffCollTracks) {
        v0dataLink(-1);
        continue;
      }

      // passes diff coll check
      registry.fill(HIST("hV0Criteria"), 4.5);

      // Act on copies for minimization
      auto pTrackCopy = o2::track::TrackParCov(pTrack);
      auto nTrackCopy = o2::track::TrackParCov(nTrack);

      //---/---/---/
      // Move close to minima
      int nCand = fitter.process(pTrackCopy, nTrackCopy);
      if (nCand == 0) {
        v0dataLink(-1);
        continue;
      }

      // passes V0 fitter minimization successfully
      registry.fill(HIST("hV0Criteria"), 5.5);

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
        v0dataLink(-1);
        continue;
      }

      // Passes step 2 of V0 fitter
      registry.fill(HIST("hV0Criteria"), 6.5);

      pTrack.getPxPyPzGlo(pvec0);
      nTrack.getPxPyPzGlo(pvec1);

      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        pos[i] = vtx[i];
      }

      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "in builder 0: posTrack --> " << labelPos << ", negTrack --> " << labelNeg);

      // Apply selections so a skimmed table is created only
      if (fitter.getChi2AtPCACandidate() > dcav0dau) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack --> " << labelPos << ", negTrack --> " << labelNeg << " will be skipped due to dca cut");
        v0dataLink(-1);
        continue;
      }

      // Passes DCA between daughters check
      registry.fill(HIST("hV0Criteria"), 7.5);

      auto V0CosinePA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
      if (V0CosinePA < v0cospa) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack --> " << labelPos << ", negTrack --> " << labelNeg << " will be skipped due to CPA cut");
        v0dataLink(-1);
        continue;
      }

      // Passes CosPA check
      registry.fill(HIST("hV0Criteria"), 8.5);

      auto V0radius = RecoDecay::sqrtSumOfSquares(pos[0], pos[1]); // probably find better name to differentiate the cut from the variable
      if (V0radius < v0radius) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack --> " << labelPos << ", negTrack --> " << labelNeg << " will be skipped due to radius cut");
        v0dataLink(-1);
        continue;
      }

      // Passes radius check
      registry.fill(HIST("hV0Criteria"), 9.5);

      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "in builder 1, keeping K0S candidate: posTrack --> " << labelPos << ", negTrack --> " << labelNeg);

      v0data(
        V0.posTrackId(),
        V0.negTrackId(),
        V0.collisionId(),
        V0.globalIndex(),
        fitter.getTrack(0).getX(), fitter.getTrack(1).getX(),
        pos[0], pos[1], pos[2],
        pvec0[0], pvec0[1], pvec0[2],
        pvec1[0], pvec1[1], pvec1[2],
        fitter.getChi2AtPCACandidate(),
        posTrackCast.dcaXY(),
        negTrackCast.dcaXY());
      v0dataLink(v0data.lastIndex());
    }
  }

  void processRun2(aod::Collision const& collision, aod::V0s const& V0s, MyTracks const& tracks, aod::BCsWithTimestamps const&
#ifdef MY_DEBUG
                   ,
                   aod::McParticles const& particlesMC
#endif
  )
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    // do v0s, typecase correctly into tracks (Run 2 use case)
    buildLambdaKZeroTable<MyTracks>(collision, V0s, kFALSE);
  }
  PROCESS_SWITCH(lambdakzeroBuilder, processRun2, "Produce Run 2 V0 tables", true);

  void processRun3(aod::Collision const& collision, aod::V0s const& V0s, MyTracksIU const& tracks, aod::BCsWithTimestamps const&
#ifdef MY_DEBUG
                   ,
                   aod::McParticles const& particlesMC
#endif
  )
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    // do v0s, typecase correctly into tracksIU (Run 3 use case)
    buildLambdaKZeroTable<MyTracksIU>(collision, V0s, kTRUE);
  }
  PROCESS_SWITCH(lambdakzeroBuilder, processRun3, "Produce Run 3 V0 tables", false);
};

struct lambdakzeroLabelBuilder {

  Produces<aod::McV0Labels> v0labels;

  // for bookkeeping purposes: how many V0s come from same mother etc
  HistogramRegistry registry{
    "registry",
    {
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hK0Short", "hK0Short", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hLambda", "hLambda", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hAntiLambda", "hAntiLambda", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
    },
  };

  void init(InitContext const&) {}

  void processDoNotBuildLabels(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(lambdakzeroLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  void processBuildLabels(aod::Collisions::iterator const& collision, aod::V0Datas const& v0table, LabeledTracks const&, aod::McParticles const& particlesMC)

  {
    for (auto& v0 : v0table) {

      int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      float lFillVal = 0.5f; // all considered V0s

      auto lNegTrack = v0.negTrack_as<LabeledTracks>();
      auto lPosTrack = v0.posTrack_as<LabeledTracks>();

      // Association check
      // There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle()) {
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {

          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother.globalIndex() == lPosMother.globalIndex()) {
                lLabel = lNegMother.globalIndex();
                lPt = lNegMother.pt();
                lPDG = lNegMother.pdgCode();
                lFillVal = 1.5f; // v0s with the same mother
              }
            }
          }
        }
      } // end association check
      registry.fill(HIST("hLabelCounter"), lFillVal);

      // Intended for cross-checks only
      // N.B. no rapidity cut!
      if (lPDG == 310)
        registry.fill(HIST("hK0Short"), lPt);
      if (lPDG == 3122)
        registry.fill(HIST("hLambda"), lPt);
      if (lPDG == -3122)
        registry.fill(HIST("hAntiLambda"), lPt);

      // Construct label table (note: this will be joinable with V0Datas)
      v0labels(
        lLabel);
    }
  }
  PROCESS_SWITCH(lambdakzeroLabelBuilder, processBuildLabels, "Produce MC label tables", false);
};

/// Extends the v0data table with expression columns
struct lambdakzeroInitializer {
  Spawns<aod::V0Datas> v0datas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeroBuilder>(cfgc),
    adaptAnalysisTask<lambdakzeroLabelBuilder>(cfgc),
    adaptAnalysisTask<lambdakzeroInitializer>(cfgc)};
}
