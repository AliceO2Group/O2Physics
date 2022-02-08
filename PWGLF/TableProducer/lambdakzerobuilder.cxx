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
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullTracksExt = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksExtended>;
using FullTracksExtMC = soa::Join<FullTracksExt, aod::McTrackLabels>;

//#define MY_DEBUG
#ifdef MY_DEBUG
using MyTracks = FullTracksExtMC;
#define MY_DEBUG_MSG(condition, cmd) \
  if (condition) {                   \
    cmd;                             \
  }
#else
using MyTracks = FullTracksExt;
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
      {"hV0Candidate", "hV0Candidate", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hGoodIndices", "hGoodIndices", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
    },
  };

  // Configurables
  Configurable<double> d_bz{"d_bz", -5.0, "bz field"};
  // Configurable<int> d_UseAbsDCA{"d_UseAbsDCA", 1, "Use Abs DCAs"}; uncomment this once we want to use the weighted DCA

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<int> rejDiffCollTracks{"rejDiffCollTracks", 0, "rejDiffCollTracks"};

  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> v_labelK0Spos{"v_labelK0Spos", {729, 2866, 4754}, "labels of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> v_labelK0Sneg{"v_labelK0Sneg", {730, 2867, 4755}, "labels of K0S positive daughters, for debug"};
#endif

  void init(InitContext& context)
  {
    // using namespace analysis::lambdakzerobuilder;

    ccdb->setURL("https://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/Geometry");
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  void process(aod::Collision const& collision, aod::V0s const& V0s, MyTracks const& tracks
#ifdef MY_DEBUG
               ,
               aod::McParticles const& particlesMC
#endif
  )
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
#ifdef MY_DEBUG
      auto labelPos = V0.posTrack_as<MyTracks>().mcParticleId();
      auto labelNeg = V0.negTrack_as<MyTracks>().mcParticleId();
      bool isK0SfromLc = isK0SfromLcFunc(labelPos, labelNeg, v_labelK0Spos, v_labelK0Sneg);
#endif
      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "V0 builder: found K0S from Lc, posTrack --> " << labelPos << ", negTrack --> " << labelNeg);

      registry.fill(HIST("hGoodIndices"), 0.5);
      if (isRun2) {
        if (!(V0.posTrack_as<MyTracks>().trackType() & o2::aod::track::TPCrefit)) {
          MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack " << labelPos << " has no TPC refit");
          v0dataLink(-1);
          continue; // TPC refit
        }
        if (!(V0.negTrack_as<MyTracks>().trackType() & o2::aod::track::TPCrefit)) {
          MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "negTrack " << labelNeg << " has no TPC refit");
          v0dataLink(-1);
          continue; // TPC refit
        }
      }
      registry.fill(HIST("hGoodIndices"), 1.5);
      if (V0.posTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack " << labelPos << " has " << V0.posTrack_as<MyTracks>().tpcNClsCrossedRows() << " crossed rows, cut at " << mincrossedrows);
        v0dataLink(-1);
        continue;
      }
      if (V0.negTrack_as<MyTracks>().tpcNClsCrossedRows() < mincrossedrows) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "negTrack " << labelNeg << " has " << V0.negTrack_as<MyTracks>().tpcNClsCrossedRows() << " crossed rows, cut at " << mincrossedrows);
        v0dataLink(-1);
        continue;
      }
      registry.fill(HIST("hGoodIndices"), 2.5);
      if (fabs(V0.posTrack_as<MyTracks>().dcaXY()) < dcapostopv) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack " << labelPos << " has dcaXY " << V0.posTrack_as<MyTracks>().dcaXY() << " , cut at " << dcanegtopv);
        v0dataLink(-1);
        continue;
      }
      if (fabs(V0.negTrack_as<MyTracks>().dcaXY()) < dcanegtopv) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "negTrack " << labelNeg << " has dcaXY " << V0.negTrack_as<MyTracks>().dcaXY() << " , cut at " << dcanegtopv);
        v0dataLink(-1);
        continue;
      }
      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "Filling good indices: posTrack --> " << labelPos << ", negTrack --> " << labelNeg);
      registry.fill(HIST("hGoodIndices"), 3.5);

      // Candidate building part
      std::array<float, 3> pos = {0.};
      std::array<float, 3> pvec0 = {0.};
      std::array<float, 3> pvec1 = {0.};

#ifdef MY_DEBUG
      auto labelPos = V0.posTrack_as<MyTracks>().mcParticleId();
      auto labelNeg = V0.negTrack_as<MyTracks>().mcParticleId();
      bool isK0SfromLc = isK0SfromLcFunc(labelPos, labelNeg, v_labelK0Spos, v_labelK0Sneg);
#endif

      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "labelPos = " << labelPos << ", labelNeg = " << labelNeg);

      registry.fill(HIST("hV0Candidate"), 0.5);

      auto pTrack = getTrackParCov(V0.posTrack_as<MyTracks>());
      auto nTrack = getTrackParCov(V0.negTrack_as<MyTracks>());

      // Require collision-ID
      if (V0.posTrack_as<MyTracks>().collisionId() != V0.negTrack_as<MyTracks>().collisionId() && rejDiffCollTracks) {
        v0dataLink(-1);
        continue;
      }

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

      o2::base::Propagator::Instance()->propagateToX(pTrack, finalXpos, d_bz, 0.85f, 2.0f, matCorr);
      o2::base::Propagator::Instance()->propagateToX(nTrack, finalXneg, d_bz, 0.85f, 2.0f, matCorr);

      nCand = fitter.process(pTrack, nTrack);
      if (nCand == 0) {
        v0dataLink(-1);
        continue;
      }

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

      auto V0CosinePA = RecoDecay::CPA(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
      if (V0CosinePA < v0cospa) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack --> " << labelPos << ", negTrack --> " << labelNeg << " will be skipped due to CPA cut");
        v0dataLink(-1);
        continue;
      }

      auto V0radius = RecoDecay::sqrtSumOfSquares(pos[0], pos[1]); // probably find better name to differentiate the cut from the variable
      if (V0radius < v0radius) {
        MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "posTrack --> " << labelPos << ", negTrack --> " << labelNeg << " will be skipped due to radius cut");
        v0dataLink(-1);
        continue;
      }

      MY_DEBUG_MSG(isK0SfromLc, LOG(info) << "in builder 1, keeping K0S candidate: posTrack --> " << labelPos << ", negTrack --> " << labelNeg);

      registry.fill(HIST("hV0Candidate"), 1.5);
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
        V0.posTrack_as<MyTracks>().dcaXY(),
        V0.negTrack_as<MyTracks>().dcaXY());
      v0dataLink(v0data.lastIndex());
    }
  }
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
    adaptAnalysisTask<lambdakzeroInitializer>(cfgc)};
}
