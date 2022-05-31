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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

//use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksExtended>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksExtended>;

using MyTracks = FullTracksExt;
using MyTracksIU = FullTracksExtIU;

/// Cascade builder task: rebuilds cascades
struct cascadeBuilder {
  Produces<aod::CascData> cascdata;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  OutputObj<TH1F> hEventCounter{TH1F("hEventCounter", "", 1, 0, 1)};
  OutputObj<TH1F> hCascCandidate{TH1F("hCascCandidate", "", 20, 0, 20)};

  // Configurables
  Configurable<double> d_bz_input{"d_bz", -999.0, "bz field"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  Configurable<int> mincrossedrows{"mincrossedrows", -1, "min crossed rows"};
  Configurable<float> dcav0topv{"dcav0topv", .1, "DCA V0 To PV"};
  Configurable<double> cospaV0{"cospaV0", .95, "CosPA V0"};
  Configurable<float> lambdamasswindow{"lambdamasswindow", .012, "Distance from Lambda mass"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", .1, "DCA Bach To PV"};
  Configurable<bool> tpcrefit{"tpcrefit", false, "demand TPC refit"};
  Configurable<double> v0radius{"v0radius", 0.9, "v0radius"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  //max sine phi for propagation
  float maxStep; //max step size (cm) for propagation

  void init(InitContext& context)
  {
    // using namespace analysis::lambdakzerobuilder;
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  //could be changed later
    maxStep = 2.00f; //could be changed later

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", run3grp_timestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
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

  void processRun2(aod::Collision const& collision, aod::V0sLinked const&, aod::V0Datas const& v0data, aod::Cascades const& cascades, FullTracksExt const&, aod::BCsWithTimestamps const&)
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

    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (bc.runNumber() != mRunNumber) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(collision.bc_as<aod::BCsWithTimestamps>().timestamp());
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = bc.runNumber();
    }

    for (auto& casc : cascades) {
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      hCascCandidate->Fill(0.5); //considered
      if (!(v0index.has_v0Data())) {
        //cascdataLink(-1);
        continue; //skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data(); //de-reference index to correct v0data in case it exists

      std::array<float, 3> pVtx = {v0.collision().posX(), v0.collision().posY(), v0.collision().posZ()};
      hCascCandidate->Fill(1.5); //has matched V0
      if (tpcrefit) {
        if (!(v0.posTrack_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascCandidate->Fill(2.5);
        if (!(v0.negTrack_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascCandidate->Fill(3.5);
        if (!(casc.bachelor_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
          //cascdataLink(-1);
          continue; // TPC refit
        }
        hCascCandidate->Fill(4.5);
      }
      if (v0.posTrack_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(5.5);
      if (v0.negTrack_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(6.5);
      if (casc.bachelor_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(7.5);
      if (fabs(v0.posTrack_as<FullTracksExt>().dcaXY()) < dcapostopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(8.5);
      if (fabs(v0.negTrack_as<FullTracksExt>().dcaXY()) < dcanegtopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(9.5);
      if (fabs(casc.bachelor_as<FullTracksExt>().dcaXY()) < dcabachtopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(10.5);

      // V0 selections
      if (fabs(v0.mLambda() - 1.116) > lambdamasswindow && fabs(v0.mAntiLambda() - 1.116) > lambdamasswindow) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(11.5);
      if (v0.dcaV0daughters() > dcav0dau) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(12.5);
      if (v0.v0radius() < v0radius) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(13.5);
      if (v0.v0cosPA(pVtx[0], pVtx[1], pVtx[2]) < cospaV0) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(14.5);
      if (v0.dcav0topv(pVtx[0], pVtx[1], pVtx[2]) < dcav0topv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(15.5);

      auto charge = -1;
      std::array<float, 3> pos = {0.};
      std::array<float, 3> posXi = {0.};
      std::array<float, 3> pvecpos = {0.};
      std::array<float, 3> pvecneg = {0.};
      std::array<float, 3> pvecbach = {0.};

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
      } else {
        //cascdataLink(-1);
        continue;
      }
      // Fill table, please
      hCascCandidate->Fill(16.5); //this is the master fill: if this is filled, viable candidate
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
  PROCESS_SWITCH(cascadeBuilder, processRun2, "Produce Run 2 multiplicity tables", true);

  void processRun3(aod::Collision const& collision, aod::V0sLinked const&, aod::V0Datas const& v0data, aod::Cascades const& cascades, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
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

    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (bc.runNumber() != mRunNumber) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(collision.bc_as<aod::BCsWithTimestamps>().timestamp());
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = bc.runNumber();
    }

    for (auto& casc : cascades) {
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      hCascCandidate->Fill(0.5); //considered
      if (!(v0index.has_v0Data())) {
        //cascdataLink(-1);
        continue; //skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data(); //de-reference index to correct v0data in case it exists

      std::array<float, 3> pVtx = {v0.collision().posX(), v0.collision().posY(), v0.collision().posZ()};
      hCascCandidate->Fill(1.5); //has matched V0
      if (tpcrefit) {
        if (!(v0.posTrack_as<FullTracksExtIU>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascCandidate->Fill(2.5);
        if (!(v0.negTrack_as<FullTracksExtIU>().trackType() & o2::aod::track::TPCrefit)) {
          continue; // TPC refit
        }
        hCascCandidate->Fill(3.5);
        if (!(casc.bachelor_as<FullTracksExtIU>().trackType() & o2::aod::track::TPCrefit)) {
          //cascdataLink(-1);
          continue; // TPC refit
        }
        hCascCandidate->Fill(4.5);
      }
      if (v0.posTrack_as<FullTracksExtIU>().tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(5.5);
      if (v0.negTrack_as<FullTracksExtIU>().tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(6.5);
      if (casc.bachelor_as<FullTracksExtIU>().tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(7.5);
      if (fabs(v0.posTrack_as<FullTracksExtIU>().dcaXY()) < dcapostopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(8.5);
      if (fabs(v0.negTrack_as<FullTracksExtIU>().dcaXY()) < dcanegtopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(9.5);
      if (fabs(casc.bachelor_as<FullTracksExtIU>().dcaXY()) < dcabachtopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(10.5);

      // V0 selections
      if (fabs(v0.mLambda() - 1.116) > lambdamasswindow && fabs(v0.mAntiLambda() - 1.116) > lambdamasswindow) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(11.5);
      if (v0.dcaV0daughters() > dcav0dau) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(12.5);
      if (v0.v0radius() < v0radius) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(13.5);
      if (v0.v0cosPA(pVtx[0], pVtx[1], pVtx[2]) < cospaV0) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(14.5);
      if (v0.dcav0topv(pVtx[0], pVtx[1], pVtx[2]) < dcav0topv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(15.5);

      auto charge = -1;
      std::array<float, 3> pos = {0.};
      std::array<float, 3> posXi = {0.};
      std::array<float, 3> pvecpos = {0.};
      std::array<float, 3> pvecneg = {0.};
      std::array<float, 3> pvecbach = {0.};

      // Acquire basic tracks
      auto pTrack = getTrackParCov(v0.posTrack_as<FullTracksExtIU>());
      auto nTrack = getTrackParCov(v0.negTrack_as<FullTracksExtIU>());
      auto bTrack = getTrackParCov(casc.bachelor_as<FullTracksExtIU>());
      if (casc.bachelor_as<FullTracksExtIU>().signed1Pt() > 0) {
        charge = +1;
      }

      int nCand = fitterV0.process(pTrack, nTrack);
      if (nCand != 0) {
        fitterV0.propagateTracksToVertex();
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
      } else {
        //cascdataLink(-1);
        continue;
      }
      // Fill table, please
      hCascCandidate->Fill(16.5); //this is the master fill: if this is filled, viable candidate
      cascdata(
        v0.globalIndex(),
        casc.bachelor_as<FullTracksExtIU>().globalIndex(),
        casc.bachelor_as<FullTracksExtIU>().collisionId(),
        charge, posXi[0], posXi[1], posXi[2], pos[0], pos[1], pos[2],
        pvecpos[0], pvecpos[1], pvecpos[2],
        pvecneg[0], pvecneg[1], pvecneg[2],
        pvecbach[0], pvecbach[1], pvecbach[2],
        fitterV0.getChi2AtPCACandidate(), fitterCasc.getChi2AtPCACandidate(),
        v0.posTrack_as<FullTracksExtIU>().dcaXY(),
        v0.negTrack_as<FullTracksExtIU>().dcaXY(),
        casc.bachelor_as<FullTracksExtIU>().dcaXY());
    }
  }
  PROCESS_SWITCH(cascadeBuilder, processRun3, "Produce Run 3 multiplicity tables", false);
};

/// Extends the cascdata table with expression columns
struct cascadeInitializer {
  Spawns<aod::CascDataExt> cascdataext;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeBuilder>(cfgc),
    adaptAnalysisTask<cascadeInitializer>(cfgc)};
}
