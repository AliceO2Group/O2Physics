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
#include "PWGLF/DataModel/LFStrangenessTables.h"
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

//use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;

//in case requested
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

/// Cascade builder task: rebuilds cascades
struct cascadeBuilder {
  Produces<aod::CascData> cascdata;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  OutputObj<TH1F> hEventCounter{TH1F("hEventCounter", "", 1, 0, 1)};
  OutputObj<TH1F> hCascCandidate{TH1F("hCascCandidate", "", 20, 0, 20)};

  // Configurables
  Configurable<double> d_bz_input{"d_bz", 5, "bz field"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  Configurable<int> mincrossedrows{"mincrossedrows", -1, "min crossed rows"};
  Configurable<float> dcav0topv{"dcav0topv", .01, "DCA V0 To PV"};
  Configurable<double> cospaV0{"cospaV0", .9, "CosPA V0"};
  Configurable<float> lambdamasswindow{"lambdamasswindow", .02, "Distance from Lambda mass"};
  Configurable<float> dcav0dau{"dcav0dau", 5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .01, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .01, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", .01, "DCA Bach To PV"};
  Configurable<bool> tpcrefit{"tpcrefit", false, "demand TPC refit"};
  Configurable<double> v0radius{"v0radius", 0.9, "v0radius"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

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
    ccdb->setFatalWhenNull(false);

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
      /* it seems this is needed at this level for the material LUT to work properly */
      /* but what happens if the run changes while doing the processing?             */
      constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;

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
      } else {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
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

  void CheckAndUpdate(Int_t lRunNumber, uint64_t lTimeStamp)
  {
    if (lRunNumber != mRunNumber) {
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = getMagneticField(lTimeStamp);
      } else {
        d_bz = d_bz_input;
      }
      mRunNumber = lRunNumber;
    }
  }

  template <class TCascTracksTo>
  void buildCascadeTable(aod::Collision const& collision, aod::V0Datas const& v0data, aod::Cascades const& cascades, Bool_t lRun3 = kTRUE)
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

    for (auto& casc : cascades) {
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      hCascCandidate->Fill(0.5); //considered
      if (!(v0.has_v0Data())) {
        //cascdataLink(-1);
        continue; //skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); //de-reference index to correct v0data in case it exists

      std::array<float, 3> pVtx = {v0data.collision().posX(), v0data.collision().posY(), v0data.collision().posZ()};

      auto bachTrackCast = casc.bachelor_as<TCascTracksTo>();
      auto posTrackCast = v0data.posTrack_as<TCascTracksTo>();
      auto negTrackCast = v0data.negTrack_as<TCascTracksTo>();

      hCascCandidate->Fill(1.5); //has matched V0
      if (tpcrefit) {
        if (!(posTrackCast.trackType() & o2::aod::track::TPCrefit) && !lRun3) {
          continue; // TPC refit
        }
        hCascCandidate->Fill(2.5);
        if (!(negTrackCast.trackType() & o2::aod::track::TPCrefit) && !lRun3) {
          continue; // TPC refit
        }
        hCascCandidate->Fill(3.5);
        if (!(bachTrackCast.trackType() & o2::aod::track::TPCrefit) && !lRun3) {
          //cascdataLink(-1);
          continue; // TPC refit
        }
        hCascCandidate->Fill(4.5);
      }
      if (posTrackCast.tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(5.5);
      if (negTrackCast.tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(6.5);
      if (bachTrackCast.tpcNClsCrossedRows() < mincrossedrows) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(7.5);
      if (fabs(posTrackCast.dcaXY()) < dcapostopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(8.5);
      if (fabs(negTrackCast.dcaXY()) < dcanegtopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(9.5);
      if (fabs(bachTrackCast.dcaXY()) < dcabachtopv) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(10.5);

      // V0 selections
      if (fabs(v0data.mLambda() - 1.116) > lambdamasswindow && fabs(v0data.mAntiLambda() - 1.116) > lambdamasswindow) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(11.5);
      if (v0data.dcaV0daughters() > dcav0dau) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(12.5);
      if (v0data.v0radius() < v0radius) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(13.5);
      if (v0data.v0cosPA(pVtx[0], pVtx[1], pVtx[2]) < cospaV0) {
        //cascdataLink(-1);
        continue;
      }
      hCascCandidate->Fill(14.5);
      if (v0data.dcav0topv(pVtx[0], pVtx[1], pVtx[2]) < dcav0topv) {
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
      auto pTrack = getTrackParCov(posTrackCast);
      auto nTrack = getTrackParCov(negTrackCast);
      auto bTrack = getTrackParCov(bachTrackCast);
      if (bachTrackCast.signed1Pt() > 0) {
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

      if (casc.collisionId() < 0)
        hCascCandidate->Fill(17.5);
      if (casc.collisionId() >= 0)
        hCascCandidate->Fill(18.5);

      cascdata(
        v0.globalIndex(),
        bachTrackCast.globalIndex(),
        casc.collisionId(),
        charge, posXi[0], posXi[1], posXi[2], pos[0], pos[1], pos[2],
        pvecpos[0], pvecpos[1], pvecpos[2],
        pvecneg[0], pvecneg[1], pvecneg[2],
        pvecbach[0], pvecbach[1], pvecbach[2],
        fitterV0.getChi2AtPCACandidate(), fitterCasc.getChi2AtPCACandidate(),
        posTrackCast.dcaXY(),
        negTrackCast.dcaXY(),
        bachTrackCast.dcaXY());
    }
  }

  void processRun2(aod::Collision const& collision, aod::V0sLinked const&, aod::V0Datas const& v0data, aod::Cascades const& cascades, FullTracksExt const&, aod::BCsWithTimestamps const&)
  {
    hEventCounter->Fill(0.5);

    // check previous run number, update if necessary
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    CheckAndUpdate(bc.runNumber(), bc.timestamp());

    // do cascades, typecase correctly into tracks
    buildCascadeTable<FullTracksExt>(collision, v0data, cascades, kFALSE);
  }
  PROCESS_SWITCH(cascadeBuilder, processRun2, "Produce Run 2 cascade tables", true);

  void processRun3(aod::Collision const& collision, aod::V0sLinked const&, aod::V0Datas const& v0data, aod::Cascades const& cascades, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    hEventCounter->Fill(0.5);

    // check previous run number, update if necessary
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    CheckAndUpdate(bc.runNumber(), bc.timestamp());

    // do cascades, typecase correctly into tracksIU (Run 3 use case)
    buildCascadeTable<FullTracksExtIU>(collision, v0data, cascades, kTRUE);
  }
  PROCESS_SWITCH(cascadeBuilder, processRun3, "Produce Run 3 cascade tables", false);
};

struct cascadeLabelBuilder {
  Produces<aod::McCascLabels> casclabels; //optionally produced

  //Bookkeeping (used for labeler)
  HistogramRegistry registry{
    "registry",
    {
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
      {"hXiMinus", "hXiMinus", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hXiPlus", "hXiPlus", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hOmegaMinus", "hOmegaMinus", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hOmegaPlus", "hOmegaPlus", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
    },
  };

  void init(InitContext const&) {}

  void processDoNotBuildLabels(aod::Collision const& collision)
  {
    //dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(cascadeLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  void processBuildLabels(aod::Collision const& collision, aod::CascDataExt const& casctable, aod::V0sLinked const&, aod::V0Datas const& v0table, LabeledTracks const&, aod::McParticles const&)
  {
    for (auto& casc : casctable) {
      float lFillVal = 0.5f; //all considered V0s
      //Loop over those that actually have the corresponding V0 associated to them
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        registry.fill(HIST("hLabelCounter"), lFillVal);
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists

      int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      lFillVal = 1.5f; //all considered V0s

      //Acquire all three daughter tracks, please
      auto lBachTrack = casc.bachelor_as<LabeledTracks>();
      auto lNegTrack = v0data.negTrack_as<LabeledTracks>();
      auto lPosTrack = v0data.posTrack_as<LabeledTracks>();

      //Association check
      //There might be smarter ways of doing this in the future
      if (lNegTrack.has_mcParticle() && lPosTrack.has_mcParticle() && lBachTrack.has_mcParticle()) {
        auto lMCBachTrack = lBachTrack.mcParticle_as<aod::McParticles>();
        auto lMCNegTrack = lNegTrack.mcParticle_as<aod::McParticles>();
        auto lMCPosTrack = lPosTrack.mcParticle_as<aod::McParticles>();
        lFillVal = 2.5f;

        //Step 1: check if the mother is the same, go up a level
        if (lMCNegTrack.has_mothers() && lMCPosTrack.has_mothers()) {
          lFillVal = 3.5f;
          for (auto& lNegMother : lMCNegTrack.mothers_as<aod::McParticles>()) {
            for (auto& lPosMother : lMCPosTrack.mothers_as<aod::McParticles>()) {
              if (lNegMother == lPosMother) {
                //if we got to this level, it means the mother particle exists and is the same
                //now we have to go one level up and compare to the bachelor mother too
                lFillVal = 4.5f;
                for (auto& lV0Mother : lNegMother.mothers_as<aod::McParticles>()) {
                  for (auto& lBachMother : lMCBachTrack.mothers_as<aod::McParticles>()) {
                    if (lV0Mother == lBachMother) {
                      lLabel = lV0Mother.globalIndex();
                      lPt = lV0Mother.pt();
                      lPDG = lV0Mother.pdgCode();
                      lFillVal = 5.5f; //v0s with the same mother
                    }
                  }
                } //end conditional V0-bach pair
              }   //end neg = pos mother conditional
            }
          } //end loop neg/pos mothers
        }   //end conditional of mothers existing
      }     //end association check

      registry.fill(HIST("hLabelCounter"), lFillVal);

      //Intended for cross-checks only
      //N.B. no rapidity cut!
      if (lPDG == 3312)
        registry.fill(HIST("hXiMinus"), lPt);
      if (lPDG == -3312)
        registry.fill(HIST("hXiPlus"), lPt);
      if (lPDG == 3334)
        registry.fill(HIST("hOmegaMinus"), lPt);
      if (lPDG == -3334)
        registry.fill(HIST("hOmegaPlus"), lPt);

      //Construct label table (note: this will be joinable with CascDatas)
      casclabels(
        lLabel);
    } //end casctable loop
  }
  PROCESS_SWITCH(cascadeLabelBuilder, processBuildLabels, "Produce MC label tables", false);
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
    adaptAnalysisTask<cascadeLabelBuilder>(cfgc),
    adaptAnalysisTask<cascadeInitializer>(cfgc)};
}
