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
// Cascade builder task
// =====================
//
// This task loops over an *existing* list of cascades (V0+bachelor track
// indices) and calculates the corresponding full cascade information
//
// Any analysis should loop over the "CascData"
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
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTables.h"

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullTracksExt = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksExtended>;

struct cascadedoindexing {
  Builds<aod::MatchedV0Cascades> var;
  void init(InitContext const&) {}
};

/// Cascade builder task: rebuilds cascades
struct cascadebuilder {
  Produces<aod::CascData> cascdata;

  OutputObj<TH1F> hEventCounter{TH1F("hEventCounter", "", 1, 0, 1)};
  OutputObj<TH1F> hCascFiltered{TH1F("hCascFiltered", "", 15, 0, 15)};
  OutputObj<TH1F> hCascCandidate{TH1F("hCascCandidate", "", 10, 0, 10)};

  // Configurables
  Configurable<double> d_bz{"d_bz", -5.0, "bz field"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  Configurable<int> mincrossedrows{"mincrossedrows", -1, "min crossed rows"};
  Configurable<float> dcav0topv{"dcav0topv", .1, "DCA V0 To PV"};
  Configurable<double> cospaV0{"cospaV0", .98, "CosPA V0"};
  Configurable<float> lambdamasswindow{"lambdamasswindow", .006, "Distance from Lambda mass"};
  Configurable<float> dcav0dau{"dcav0dau", .6, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", .1, "DCA Bach To PV"};
  Configurable<bool> tpcrefit{"tpcrefit", true, "demand TPC refit"};
  Configurable<double> v0radius{"v0radius", 0.9, "v0radius"};

  void process(aod::MatchedV0Cascades const& MatchedV0Cascades, aod::V0Datas const&, aod::Cascades const&, aod::Collisions const&, FullTracksExt const&)
  {
    // Define o2 fitter, 2-prong
    o2::vertexing::DCAFitterN<2> fitterV0, fitterCasc;
    fitterV0.setBz(d_bz);
    fitterV0.setPropagateToPCA(true);
    fitterV0.setMaxR(200.);
    fitterV0.setMinParamChange(1e-3);
    fitterV0.setMinRelChi2Change(0.9);
    fitterV0.setMaxDZIni(1e9);
    fitterV0.setMaxChi2(1e9);
    fitterV0.setUseAbsDCA(d_UseAbsDCA);

    fitterCasc.setBz(d_bz);
    fitterCasc.setPropagateToPCA(true);
    fitterCasc.setMaxR(200.);
    fitterCasc.setMinParamChange(1e-3);
    fitterCasc.setMinRelChi2Change(0.9);
    fitterCasc.setMaxDZIni(1e9);
    fitterCasc.setMaxChi2(1e9);
    fitterCasc.setUseAbsDCA(d_UseAbsDCA);

    hEventCounter->Fill(0.5);

    for (auto& cascIndexLUT : MatchedV0Cascades) {
      auto v0 = cascIndexLUT.v0data();
      auto casc = cascIndexLUT.cascade();

      std::array<float, 3> pVtx = {v0.collision().posX(), v0.collision().posY(), v0.collision().posZ()};

      if (tpcrefit) {
        if (!(v0.posTrack_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascFiltered->Fill(1.5);
        if (!(v0.negTrack_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascFiltered->Fill(2.5);
        if (!(casc.bachelor_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascFiltered->Fill(3.5);
      }
      if (v0.posTrack_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      hCascFiltered->Fill(4.5);
      if (v0.negTrack_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      hCascFiltered->Fill(5.5);
      if (casc.bachelor_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      hCascFiltered->Fill(6.5);
      if (v0.posTrack_as<FullTracksExt>().dcaXY() < dcapostopv) {
        continue;
      }
      hCascFiltered->Fill(7.5);
      if (v0.negTrack_as<FullTracksExt>().dcaXY() < dcanegtopv) {
        continue;
      }
      hCascFiltered->Fill(8.5);
      if (casc.bachelor_as<FullTracksExt>().dcaXY() < dcabachtopv) {
        continue;
      }
      hCascFiltered->Fill(9.5);

      // V0 selections
      if (fabs(v0.mLambda() - 1.116) > lambdamasswindow && fabs(v0.mAntiLambda() - 1.116) > lambdamasswindow) {
        continue;
      }
      hCascFiltered->Fill(10.5);
      if (v0.dcaV0daughters() > dcav0dau) {
        continue;
      }
      hCascFiltered->Fill(11.5);
      if (v0.v0radius() < v0radius) {
        continue;
      }
      hCascFiltered->Fill(12.5);
      if (v0.v0cosPA(pVtx[0], pVtx[1], pVtx[2]) < cospaV0) {
        continue;
      }
      hCascFiltered->Fill(13.5);
      if (v0.dcav0topv(pVtx[0], pVtx[1], pVtx[2]) < dcav0topv) {
        continue;
      }
      hCascFiltered->Fill(14.5);

      auto charge = -1;
      std::array<float, 3> pos = {0.};
      std::array<float, 3> posXi = {0.};
      std::array<float, 3> pvecpos = {0.};
      std::array<float, 3> pvecneg = {0.};
      std::array<float, 3> pvecbach = {0.};

      hCascCandidate->Fill(0.5);

      // Acquire basic tracks
      auto pTrack = getTrackParCov(v0.posTrack_as<FullTracksExt>());
      auto nTrack = getTrackParCov(v0.negTrack_as<FullTracksExt>());
      auto bTrack = getTrackParCov(casc.bachelor_as<FullTracksExt>());
      if (casc.bachelor_as<FullTracksExt>().signed1Pt() > 0) {
        charge = +1;
      }

      int nCand = fitterV0.process(pTrack, nTrack);
      if (nCand != 0) {
        fitterV0.propagateTracksToVertex();
        hCascCandidate->Fill(1.5);
        const auto& v0vtx = fitterV0.getPCACandidate();
        for (int i = 0; i < 3; i++) {
          pos[i] = v0vtx[i];
        }

        std::array<float, 21> cov0 = {0};
        std::array<float, 21> cov1 = {0};
        std::array<float, 21> covV0 = {0};

        // Covariance matrix calculation
        const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        fitterV0.getTrack(0).getPxPyPzGlo(pvecpos);
        fitterV0.getTrack(1).getPxPyPzGlo(pvecneg);
        fitterV0.getTrack(0).getCovXYZPxPyPzGlo(cov0);
        fitterV0.getTrack(1).getCovXYZPxPyPzGlo(cov1);
        for (int i = 0; i < 6; i++) {
          int j = momInd[i];
          covV0[j] = cov0[j] + cov1[j];
        }
        auto covVtxV0 = fitterV0.calcPCACovMatrix();
        covV0[0] = covVtxV0(0, 0);
        covV0[1] = covVtxV0(1, 0);
        covV0[2] = covVtxV0(1, 1);
        covV0[3] = covVtxV0(2, 0);
        covV0[4] = covVtxV0(2, 1);
        covV0[5] = covVtxV0(2, 2);

        const std::array<float, 3> vertex = {(float)v0vtx[0], (float)v0vtx[1], (float)v0vtx[2]};
        const std::array<float, 3> momentum = {pvecpos[0] + pvecneg[0], pvecpos[1] + pvecneg[1], pvecpos[2] + pvecneg[2]};

        auto tV0 = o2::track::TrackParCov(vertex, momentum, covV0, 0);
        tV0.setQ2Pt(0); // No bending, please
        int nCand2 = fitterCasc.process(tV0, bTrack);
        if (nCand2 != 0) {
          fitterCasc.propagateTracksToVertex();
          hCascCandidate->Fill(2.5);
          const auto& cascvtx = fitterCasc.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            posXi[i] = cascvtx[i];
          }
          fitterCasc.getTrack(1).getPxPyPzGlo(pvecbach);
        } // end if cascade recoed
      }   // end if v0 recoed
      // Fill table, please
      cascdata(
        v0.globalIndex(),
        casc.bachelor_as<FullTracksExt>().globalIndex(),
        casc.bachelor_as<FullTracksExt>().collisionId(),
        charge, posXi[0], posXi[1], posXi[2], pos[0], pos[1], pos[2],
        pvecpos[0], pvecpos[1], pvecpos[2],
        pvecneg[0], pvecneg[1], pvecneg[2],
        pvecbach[0], pvecbach[1], pvecbach[2],
        fitterV0.getChi2AtPCACandidate(), fitterCasc.getChi2AtPCACandidate(),
        v0.posTrack_as<FullTracksExt>().dcaXY(),
        v0.negTrack_as<FullTracksExt>().dcaXY(),
        casc.bachelor_as<FullTracksExt>().dcaXY());
    }
  }
};

/// Extends the cascdata table with expression columns
struct cascadeinitializer {
  Spawns<aod::CascDataExt> cascdataext;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadedoindexing>(cfgc, TaskName{"lf-cascadedoindexing"}),
    adaptAnalysisTask<cascadebuilder>(cfgc, TaskName{"lf-cascadebuilder"}),
    adaptAnalysisTask<cascadeinitializer>(cfgc, TaskName{"lf-cascadeinitializer"})};
}