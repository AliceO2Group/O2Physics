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
/// \brief this task allows for the direct one-to-one comparison of
//         cascades computed with standard DCAFitter methods and the KFparticle
//         package. It is meant for the purposes of larger-scale QA of KF reco.

/// \brief cascadebuilder.cxx and lambdakzerobuilder.cxx tasks need to be added to the workflow. Flag createCascCovMats needs to be enabled!

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/kfStrangenessStudy.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include <cmath>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// allows for candidate-by-candidate comparison using Cascade to []CascData link table
using CascadesCrossLinked = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink>;
using CascDataLabeled = soa::Join<aod::CascDatas, aod::CascCovs, aod::McCascLabels>;
using KFCascDataLabeled = soa::Join<aod::KFCascDatas, aod::KFCascCovs, aod::McKFCascLabels>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksCovIU>;
using V0s = aod::V0sLinked;

struct kfStrangenessStudy {

  Produces<aod::CascCand> rowCasc;
  Produces<aod::CascCandMC> rowCascMC;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Filters
  Filter collisionFilter = (aod::evsel::sel8 == true);

  /// Cascade data
  int isDCAfitter = 0, isKF = 0;
  float ptRec = -1.0f, ptRecKF = -1.0f;
  float ptRecV0 = -1.0f, ptRecKFV0 = -1.0f;
  float massXi = -1.0f, massXiKF = -1.0f;
  float massLambda = -1.0f, massLambdaKF = -1.0f;
  float cascRad = -1.0f, cascRadKF = -1.0f;
  float vtxXrec = -1.0f, vtxYrec = -1.0f, vtxZrec = -1.0f, vtxXrecKF = -1.0f, vtxYrecKF = -1.0f, vtxZrecKF = -1.0f;
  float vtxXrecErr = -1.0f, vtxYrecErr = -1.0f, vtxZrecErr = -1.0f, vtxXrecErrKF = -1.0f, vtxYrecErrKF = -1.0f, vtxZrecErrKF = -1.0f;
  float vtxXrecV0 = -1.0f, vtxYrecV0 = -1.0f, vtxZrecV0 = -1.0f, vtxXrecKFV0 = -1.0f, vtxYrecKFV0 = -1.0f, vtxZrecKFV0 = -1.0f;
  float vtxXrecErrV0 = -1.0f, vtxYrecErrV0 = -1.0f, vtxZrecErrV0 = -1.0f, vtxXrecErrKFV0 = -1.0f, vtxYrecErrKFV0 = -1.0f, vtxZrecErrKFV0 = -1.0f;
  float dcaXYCascToPV = -1.0f, dcaXYCascToPVKF = -1.0f;
  float dcaZCascToPV = -1.0f, dcaZCascToPVKF = -1.0f;
  float dcaCascDaughters = -1.0f, dcaCascDaughtersKF = -1.0f;
  float dcaV0Daughters = -1.0f, dcaV0DaughtersKF = -1.0f;
  float dcaPosToPV = -1.0f, dcaPosToPVKF = -1.0f;
  float dcaNegToPV = -1.0f, dcaNegToPVKF = -1.0f;
  float dcaBachToPV = -1.0f, dcaBachToPVKF = -1.0f;
  float cascPointingAngle = -1.0f, cascPointingAngleKF = -1.0f;
  float v0PointingAngle = -1.0f, v0PointingAngleKF = -1.0f;
  float V0Rad = -1.0f, V0RadKF = -1.0f;
  int charge = 0;
  std::array<float, 3> momPosRecIU;
  std::array<float, 3> momNegRecIU;
  std::array<float, 3> momPosRecIUErr;
  std::array<float, 3> momNegRecIUErr;
  std::array<float, 3> momPosRec;
  std::array<float, 3> momNegRec;
  std::array<float, 3> momPosRecErr;
  std::array<float, 3> momNegRecErr;
  std::array<float, 3> posPosRec;
  std::array<float, 3> posNegRec;
  std::array<float, 3> posPosRecErr;
  std::array<float, 3> posNegRecErr;

  /// Additional cascade MC data
  int isTrueCasc = 0;
  float ptGen = -1.0f;
  float ptGenV0 = -1.0f;
  float vtxXgen = -1.0f, vtxYgen = -1.0f, vtxZgen = -1.0f;
  float vtxXgenV0 = -1.0f, vtxYgenV0 = -1.0f, vtxZgenV0 = -1.0f;
  float vtxXgen_firstDau = -1., vtxYgen_firstDau = -1., vtxZgen_firstDau = -1;
  float vtxXgen_firstV0Dau = -1., vtxYgen_firstV0Dau = -1., vtxZgen_firstV0Dau = -1;
  float prodVtxXgen  = -1.0f, prodVtxYgen  = -1.0f, prodVtxZgen  = -1.0f;
  int source = 0;
  std::array<float, 3> momPosGen;
  std::array<float, 3> momNegGen;

  // counter and checks
  int recocase = 0;


  void init(InitContext const&)
  {
    /// QA histos
    histos.add("hChargeCounter", "hChargeCounter", kTH1F, {{3, -1.5f, 1.5f}});
    histos.add("hChargeCounterCascDatas", "hChargeCounterCascDatas", kTH1F, {{3, -1.5f, 1.5f}});
    auto hEventSelectionFlow = histos.add<TH1>("hEventSelectionFlow", "Event selection flow", kTH1F, {{2, 0.5f, 2.5f}});
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(1), "Sel8");
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(2), "|Vtx_{z}|<10cm");

    histos.add("hVertexX", "hVertexX", kTH1F, {{1000, -3.0f, 3.0f}});
    histos.add("hKFVertexX", "hKFVertexX", kTH1F, {{1000, -3.0f, 3.0f}});
    histos.add("hVertexY", "hVertexY", kTH1F, {{1000, -3.0f, 3.0f}});
    histos.add("hKFVertexY", "hKFVertexY", kTH1F, {{1000, -3.0f, 3.0f}});
    histos.add("hVertexZ", "hVertexZ", kTH1F, {{1000, -10.0f, 10.0f}});
    histos.add("hKFVertexZ", "hKFVertexZ", kTH1F, {{1000, -10.0f, 10.0f}});

    histos.add("hCascRadius", "hCascRadius", kTH1F, {{1000, 0.0f, 3.0f}});
    histos.add("hKFCascRadius", "hKFCascRadius", kTH1F, {{1000, 0.0f, 3.0f}});

    histos.add("hDCAxy", "hDCAxy", kTH1F, {{500, -1.0f, 1.0f}});
    histos.add("hKFDCAxy", "hKFDCAxy", kTH1F, {{500, -1.0f, 1.0f}});

    histos.add("hPointingAngle", "hPointingAngle", kTH1F, {{800, 0.0f, 3.5f}});
    histos.add("hKFPointingAngle", "hKFPointingAngle", kTH1F, {{800, 0.0f, 3.5f}});
    histos.add("hCosPointingAngle", "hCosPointingAngle", kTH1F, {{800, -1.0f, 1.0f}});
    histos.add("hKFCosPointingAngle", "hKFCosPointingAngle", kTH1F, {{800, -1.0f, 1.0f}});
    histos.add("hV0PointingAngle", "hV0PointingAngle", kTH1F, {{800, 0.0f, 3.5f}});
    histos.add("hKFV0PointingAngle", "hKFV0PointingAngle", kTH1F, {{800, 0.0f, 3.5f}});
    histos.add("hCosV0PointingAngle", "hCosV0PointingAngle", kTH1F, {{800, -1.0f, 1.0f}});
    histos.add("hKFCosV0PointingAngle", "hKFCosV0PointingAngle", kTH1F, {{800, -1.0f, 1.0f}});

    histos.add("hGenDecayVtxX_firstDau", "hGenDecayVtxX_firstDau", kTH1F, {{1000, -3.0f, 3.0f}});
    histos.add("hGenDecayVtxY_firstDau", "hGenDecayVtxY_firstDau", kTH1F, {{1000, -3.0f, 3.0f}});
    histos.add("hGenDecayVtxZ_firstDau", "hGenDecayVtxZ_firstDau", kTH1F, {{1000, -3.0f, 3.0f}});

    histos.add("hGenSource", "hGenSource", kTH1F, {{5, -2, 3}});
    histos.add("hCase", "hCase", kTH1F, {{5, 0, 5}});
  }

  template <typename TCollision, typename TCascade, typename TCascDatas, typename TKFCascDatas, typename TV0, typename TV0Datas, typename TV0fCDatas>
  void getCascDatas(TCollision const& collision, TCascade const& cascade, TCascDatas const&, TKFCascDatas const&, TV0 const&, TV0Datas const&, TV0fCDatas const&)
  {
    if (cascade.has_cascData()) {
      LOG(info) << "Cascade has CascData!";
      // check aod::Cascades -> aod::CascData link
      // if present: this candidate was accepted by default DCAfitter building
      isDCAfitter = 1;
      auto cascdata = cascade.template cascData_as<TCascDatas>();
      auto v0index = cascade.template v0_as<TV0>();
      if (v0index.has_v0Data()) {
        // V0 passed both standard and cascade V0 selections
        auto v0data = v0index.template v0Data_as<TV0Datas>();
        vtxXrecErrV0 = sqrt(v0data.positionCovMat()[0]);
        vtxYrecErrV0 = sqrt(v0data.positionCovMat()[2]);
        vtxZrecErrV0 = sqrt(v0data.positionCovMat()[5]);
      } else if (v0index.has_v0fCData()) {
        // V0 passed only cascade V0 selections, use this instead
        auto v0data = v0index.template v0fCData_as<TV0fCDatas>();
        vtxXrecErrV0 = sqrt(v0data.positionCovMat()[0]);
        vtxYrecErrV0 = sqrt(v0data.positionCovMat()[2]);
        vtxZrecErrV0 = sqrt(v0data.positionCovMat()[5]);
      }
      ptRec = cascdata.pt();
      vtxXrec = cascdata.x();
      vtxYrec = cascdata.y();
      vtxZrec = cascdata.z();
      vtxXrecErr = sqrt(cascdata.positionCovMat()[0]);
      vtxYrecErr = sqrt(cascdata.positionCovMat()[2]);
      vtxZrecErr = sqrt(cascdata.positionCovMat()[5]);
      ptRecV0 = sqrt((cascdata.pxpos()+cascdata.pxneg())*(cascdata.pxpos()+cascdata.pxneg()) + (cascdata.pypos()+cascdata.pyneg())*(cascdata.pypos()+cascdata.pyneg())); // taken from daughters from lambdakzerobuiler
      vtxXrecV0 = cascdata.xlambda();
      vtxYrecV0 = cascdata.ylambda();
      vtxZrecV0 = cascdata.zlambda();
      massLambda = cascdata.mLambda();
      massXi = cascdata.mXi();
      dcaXYCascToPV = cascdata.dcaXYCascToPV();
      dcaZCascToPV = cascdata.dcaZCascToPV();
      dcaCascDaughters = cascdata.dcacascdaughters();
      dcaV0Daughters = cascdata.dcaV0daughters();
      dcaPosToPV = cascdata.dcapostopv();
      dcaNegToPV = cascdata.dcanegtopv();
      dcaBachToPV = cascdata.dcabachtopv();
      cascPointingAngle = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRad = cascdata.cascradius();
      V0Rad = cascdata.v0radius();
      v0PointingAngle = TMath::ACos(cascdata.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));

      // fill QA histos
      histos.fill(HIST("hVertexX"), vtxXrec);
      histos.fill(HIST("hVertexY"), vtxYrec);
      histos.fill(HIST("hVertexZ"), vtxZrec);
      histos.fill(HIST("hDCAxy"), dcaXYCascToPV);
      histos.fill(HIST("hCascRadius"), cascRad);
      histos.fill(HIST("hPointingAngle"), cascPointingAngle);
      histos.fill(HIST("hCosPointingAngle"), cos(cascPointingAngle));
      histos.fill(HIST("hV0PointingAngle"), v0PointingAngle);
      histos.fill(HIST("hCosV0PointingAngle"), cos(v0PointingAngle));
    }

    if (cascade.has_kfCascData()) {
      LOG(info) << "Cascade has KFcascData!";
      // check aod::Cascades -> aod::KFCascData link
      // if present: this candidate was accepted by KF building
      isKF = 1;
      auto cascdatakf = cascade.template kfCascData_as<TKFCascDatas>();
      ptRecKF = cascdatakf.pt();
      vtxXrecKF = cascdatakf.x();
      vtxYrecKF = cascdatakf.y();
      vtxZrecKF = cascdatakf.z();
      vtxXrecErrKF = sqrt(cascdatakf.kfTrackCovMat()[0]);
      vtxYrecErrKF = sqrt(cascdatakf.kfTrackCovMat()[2]);
      vtxZrecErrKF = sqrt(cascdatakf.kfTrackCovMat()[5]);
      ptRecKFV0 = sqrt(cascdatakf.pxv0()*cascdatakf.pxv0() + cascdatakf.pyv0()*cascdatakf.pyv0()); // taken from KFV0
      vtxXrecKFV0 = cascdatakf.xlambda();
      vtxYrecKFV0 = cascdatakf.ylambda();
      vtxZrecKFV0 = cascdatakf.zlambda();
      vtxXrecErrKFV0 = sqrt(cascdatakf.kfTrackCovMatV0()[0]);
      vtxYrecErrKFV0 = sqrt(cascdatakf.kfTrackCovMatV0()[2]);
      vtxZrecErrKFV0 = sqrt(cascdatakf.kfTrackCovMatV0()[5]);
      massLambdaKF = cascdatakf.mLambda();
      massXiKF = cascdatakf.mXi();
      dcaXYCascToPVKF = cascdatakf.dcaXYCascToPV();
      dcaZCascToPVKF = cascdatakf.dcaZCascToPV();
      dcaCascDaughtersKF = cascdatakf.dcacascdaughters();
      dcaV0DaughtersKF = cascdatakf.dcaV0daughters();
      dcaPosToPVKF = cascdatakf.dcapostopv();
      dcaNegToPVKF = cascdatakf.dcanegtopv();
      dcaBachToPVKF = cascdatakf.dcabachtopv();
      cascPointingAngleKF = TMath::ACos(cascdatakf.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRadKF = cascdatakf.cascradius();
      V0RadKF = cascdatakf.v0radius();
      v0PointingAngleKF = TMath::ACos(cascdatakf.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      // daughter momenta at vertex (from pre-minimisation if enabled)
      momPosRec[0] = cascdatakf.pxpos();
      momPosRec[1] = cascdatakf.pypos();
      momPosRec[2] = cascdatakf.pzpos();
      momNegRec[0] = cascdatakf.pxneg();
      momNegRec[1] = cascdatakf.pyneg();
      momNegRec[2] = cascdatakf.pzneg();
      momPosRecErr[0] = cascdatakf.pxposerr();
      momPosRecErr[1] = cascdatakf.pyposerr();
      momPosRecErr[2] = cascdatakf.pzposerr();
      momNegRecErr[0] = cascdatakf.pxnegerr();
      momNegRecErr[1] = cascdatakf.pynegerr();
      momNegRecErr[2] = cascdatakf.pznegerr();
      // daughter track position at vertex (after pre-minimisation!)
      posPosRec[0] = cascdatakf.xpos();
      posPosRec[1] = cascdatakf.ypos();
      posPosRec[2] = cascdatakf.zpos();
      posNegRec[0] = cascdatakf.xneg();
      posNegRec[1] = cascdatakf.yneg();
      posNegRec[2] = cascdatakf.zneg();
      posPosRecErr[0] = cascdatakf.xposerr();
      posPosRecErr[1] = cascdatakf.yposerr();
      posPosRecErr[2] = cascdatakf.zposerr();
      posNegRecErr[0] = cascdatakf.xnegerr();
      posNegRecErr[1] = cascdatakf.ynegerr();
      posNegRecErr[2] = cascdatakf.znegerr();

      // fill QA histos
      histos.fill(HIST("hKFVertexX"), vtxXrecKF);
      histos.fill(HIST("hKFVertexY"), vtxYrecKF);
      histos.fill(HIST("hKFVertexZ"), vtxZrecKF);
      histos.fill(HIST("hKFDCAxy"), dcaXYCascToPVKF);
      histos.fill(HIST("hKFCascRadius"), cascRadKF);
      histos.fill(HIST("hKFPointingAngle"), cascPointingAngleKF);
      histos.fill(HIST("hKFCosPointingAngle"), cos(cascPointingAngleKF));
      histos.fill(HIST("hKFV0PointingAngle"), v0PointingAngleKF);
      histos.fill(HIST("hKFCosV0PointingAngle"), cos(v0PointingAngleKF));
    }
  }

  template <typename TCollision, typename TCascData, typename TV0s, typename TMCParticle>
  void getCascMCdata(TCollision const& collision, TCascData const& cascdata, TV0s const&, TMCParticle const& mcparticles)
  {
    if (cascdata.has_mcParticle() && cascdata.mcParticleId() > -1 && cascdata.mcParticleId() <= mcparticles.size()) {
      auto MCcascade = cascdata.template mcParticle_as<TMCParticle>();

      if (MCcascade.has_daughters()) {
        LOG(info) << "MC cascade has daughters, getting MC info.";
        // get MC V0
        auto MCv0 = MCcascade.template daughters_as<TMCParticle>().begin();
        if (!(abs(MCv0.pdgCode()) == 3122)) {
          LOG(info) << "V0 is no Lambda";
        }
        // cascade
        ptGen = MCcascade.pt();
        prodVtxXgen = MCcascade.vx();
        prodVtxYgen = MCcascade.vy();
        prodVtxZgen = MCcascade.vz();
        vtxXgen_firstDau = MCv0.vx();
        vtxYgen_firstDau = MCv0.vy();
        vtxZgen_firstDau = MCv0.vz();
        // V0
        ptGenV0 = MCv0.pt();
        vtxXgen_firstV0Dau = MCv0.template daughters_as<TMCParticle>().begin().vx();
        vtxYgen_firstV0Dau = MCv0.template daughters_as<TMCParticle>().begin().vy();
        vtxZgen_firstV0Dau = MCv0.template daughters_as<TMCParticle>().begin().vz();
        // daughters
        momPosGen[0] = MCv0.template daughters_as<TMCParticle>().begin().px();
        momPosGen[1] = MCv0.template daughters_as<TMCParticle>().begin().py();
        momPosGen[2] = MCv0.template daughters_as<TMCParticle>().begin().pz();
        momNegGen[0] = (MCv0.template daughters_as<TMCParticle>().begin() + 1).px();
        momNegGen[1] = (MCv0.template daughters_as<TMCParticle>().begin() + 1).py();
        momNegGen[2] = (MCv0.template daughters_as<TMCParticle>().begin() + 1).pz();

        if (abs(MCcascade.pdgCode()) == 3312) { // Xi
          isTrueCasc = 1;
        } else {
          isTrueCasc = 0;
        }
        if (MCcascade.isPhysicalPrimary()) {
          source = 1;
        } else if (MCcascade.getProcess() == 4) { // from particle decay
          source = 2;
        } else if (MCcascade.fromBackgroundEvent()) {
          source = -1;
        } else {
          source = -2;
        }
        // fill cascade table
        fillCascMCTable(collision);

        // fill QA histos --> vertex position from daughters!
        histos.fill(HIST("hGenDecayVtxX_firstDau"), vtxXgen_firstDau);
        histos.fill(HIST("hGenDecayVtxY_firstDau"), vtxYgen_firstDau);
        histos.fill(HIST("hGenDecayVtxZ_firstDau"), vtxZgen_firstDau);

        histos.fill(HIST("hGenSource"), source);
        histos.fill(HIST("hChargeCounterCascDatas"), charge);
      }
    }
  }

  template <typename TCollision>
  void fillCascDataTable(TCollision const& collision)
  {
    rowCasc(collision.globalIndex(),
            ptRec, ptRecKF,
            massXi, massXiKF,
            cascRad, cascRadKF,
            vtxXrec, vtxYrec, vtxZrec, vtxXrecErr, vtxYrecErr, vtxZrecErr,
            vtxXrecKF, vtxYrecKF, vtxZrecKF, vtxXrecErrKF, vtxYrecErrKF, vtxZrecErrKF,
            dcaXYCascToPV, dcaXYCascToPVKF,
            dcaZCascToPV, dcaZCascToPVKF,
            dcaCascDaughters, dcaCascDaughtersKF,
            cascPointingAngle, cascPointingAngleKF,
            charge,
            ptRecV0, ptRecKFV0,
            massLambda, massLambdaKF,
            V0Rad, V0RadKF,
            vtxXrecV0, vtxYrecV0, vtxZrecV0, vtxXrecErrV0, vtxYrecErrV0, vtxZrecErrV0,
            vtxXrecKFV0, vtxYrecKFV0, vtxZrecKFV0, vtxXrecErrKFV0, vtxYrecErrKFV0, vtxZrecErrKFV0,
            dcaV0Daughters, dcaV0DaughtersKF,
            dcaPosToPV, dcaPosToPVKF,
            dcaNegToPV, dcaNegToPVKF,
            dcaBachToPV, dcaBachToPVKF,
            v0PointingAngle, v0PointingAngleKF,
            momPosRecIU[0], momPosRecIU[1], momPosRecIU[2], momPosRecIUErr[0], momPosRecIUErr[1], momPosRecIUErr[2],
            momNegRecIU[0], momNegRecIU[1], momNegRecIU[2], momNegRecIUErr[0], momNegRecIUErr[1], momNegRecIUErr[2],
            momPosRec[0], momPosRec[1], momPosRec[2], momPosRecErr[0], momPosRecErr[1], momPosRecErr[2],
            momNegRec[0], momNegRec[1], momNegRec[2], momNegRecErr[0], momNegRecErr[1], momNegRecErr[2],
            isDCAfitter, isKF);
  }

  template <typename TCollision>
  void fillCascMCTable(TCollision const& collision)
  {
    rowCascMC(collision.globalIndex(),
              ptRec, ptRecKF, ptGen,
              massXi, massXiKF,
              cascRad, cascRadKF,
              vtxXrec, vtxYrec, vtxZrec, vtxXrecErr, vtxYrecErr, vtxZrecErr,
              vtxXrecKF, vtxYrecKF, vtxZrecKF, vtxXrecErrKF, vtxYrecErrKF, vtxZrecErrKF,
              vtxXgen_firstDau, vtxYgen_firstDau, vtxZgen_firstDau,
              prodVtxXgen, prodVtxYgen, prodVtxZgen,
              dcaXYCascToPV, dcaXYCascToPVKF,
              dcaZCascToPV, dcaZCascToPVKF,
              dcaCascDaughters, dcaCascDaughtersKF,
              cascPointingAngle, cascPointingAngleKF,
              charge,
              ptRecV0, ptRecKFV0, ptGenV0,
              massLambda, massLambdaKF,
              V0Rad, V0RadKF,
              vtxXrecV0, vtxYrecV0, vtxZrecV0, vtxXrecErrV0, vtxYrecErrV0, vtxZrecErrV0,
              vtxXrecKFV0, vtxYrecKFV0, vtxZrecKFV0, vtxXrecErrKFV0, vtxYrecErrKFV0, vtxZrecErrKFV0,
              vtxXgen_firstV0Dau, vtxYgen_firstV0Dau, vtxZgen_firstV0Dau,
              dcaV0Daughters, dcaV0DaughtersKF,
              dcaPosToPV, dcaPosToPVKF,
              dcaNegToPV, dcaNegToPVKF,
              dcaBachToPV, dcaBachToPVKF,
              v0PointingAngle, v0PointingAngleKF,
              momPosRecIU[0], momPosRecIU[1], momPosRecIU[2], momPosRecIUErr[0], momPosRecIUErr[1], momPosRecIUErr[2],
              momNegRecIU[0], momNegRecIU[1], momNegRecIU[2], momNegRecIUErr[0], momNegRecIUErr[1], momNegRecIUErr[2],
              momPosRec[0], momPosRec[1], momPosRec[2], momPosRecErr[0], momPosRecErr[1], momPosRecErr[2],
              momNegRec[0], momNegRec[1], momNegRec[2], momNegRecErr[0], momNegRecErr[1], momNegRecErr[2],
              momPosGen[0], momPosGen[1], momPosGen[2],
              momNegGen[0], momNegGen[1], momNegGen[2],
              posPosRec[0], posPosRec[1], posPosRec[2], posPosRecErr[0], posPosRecErr[1], posPosRecErr[2],
              posNegRec[0], posNegRec[1], posNegRec[2], posNegRecErr[0], posNegRecErr[1], posNegRecErr[2],
              isDCAfitter, isKF,
              isTrueCasc,
              source);
    LOG(info) << "CascMCTable filled!";
    histos.fill(HIST("hCase"), recocase);
  }

  void resetVars()
  {
    /// Cascade data
    isDCAfitter = 0; isKF = 0;
    ptRec = -1.0f; ptRecKF = -1.0f;
    ptRecV0 = -1.0f; ptRecKFV0 = -1.0f;
    massXi = -1.0f; massXiKF = -1.0f;
    massLambda = -1.0f; massLambdaKF = -1.0f;
    cascRad = -1.0f; cascRadKF = -1.0f;
    vtxXrec = -1.0f; vtxYrec = -1.0f; vtxZrec = -1.0f; vtxXrecKF = -1.0f; vtxYrecKF = -1.0f; vtxZrecKF = -1.0f;
    vtxXrecErr = -1.0f; vtxYrecErr = -1.0f; vtxZrecErr = -1.0f; vtxXrecErrKF = -1.0f; vtxYrecErrKF = -1.0f; vtxZrecErrKF = -1.0f;
    vtxXrecV0 = -1.0f; vtxYrecV0 = -1.0f; vtxZrecV0 = -1.0f; vtxXrecKFV0 = -1.0f; vtxYrecKFV0 = -1.0f; vtxZrecKFV0 = -1.0f;
    vtxXrecErrV0 = -1.0f; vtxYrecErrV0 = -1.0f; vtxZrecErrV0 = -1.0f; vtxXrecErrKFV0 = -1.0f; vtxYrecErrKFV0 = -1.0f; vtxZrecErrKFV0 = -1.0f;
    dcaXYCascToPV = -1.0f; dcaXYCascToPVKF = -1.0f;
    dcaZCascToPV = -1.0f; dcaZCascToPVKF = -1.0f;
    dcaCascDaughters = -1.0f; dcaCascDaughtersKF = -1.0f;
    dcaV0Daughters = -1.0f; dcaV0DaughtersKF = -1.0f;
    dcaPosToPV = -1.0f; dcaPosToPVKF = -1.0f;
    dcaNegToPV = -1.0f; dcaNegToPVKF = -1.0f;
    dcaBachToPV = -1.0f; dcaBachToPVKF = -1.0f;
    cascPointingAngle = -1.0f; cascPointingAngleKF = -1.0f;
    v0PointingAngle = -1.0f; v0PointingAngleKF = -1.0f;
    V0Rad = -1.0f; V0RadKF = -1.0f;
    charge = 0;
    for (int i = 0; i < 3; i++) {
      momPosRecIU[i] = -1.0f;
      momNegRecIU[i] = -1.0f;
      momPosRecIUErr[i] = -1.0f;
      momNegRecIUErr[i] = -1.0f;
      momPosRec[i] = -1.0f;
      momNegRec[i] = -1.0f;
      momPosRecErr[i] = -1.0f;
      momNegRecErr[i] = -1.0f;
      posPosRec[i] = -1.0f;
      posNegRec[i] = -1.0f;
      posPosRecErr[i] = -1.0f;
      posNegRecErr[i] = -1.0f;
      
      // Additional cascade MC data
      momPosGen[i] = -1.0f;
      momNegGen[i] = -1.0f;
    }

    /// Additional cascade MC data
    isTrueCasc = 0;
    ptGen = -1.0f;
    ptGenV0 = -1.0f;
    vtxXgen = -1.0f; vtxYgen = -1.0f; vtxZgen = -1.0f;
    vtxXgenV0 = -1.0f; vtxYgenV0 = -1.0f; vtxZgenV0 = -1.0f;
    vtxXgen_firstDau = -1.; vtxYgen_firstDau = -1.; vtxZgen_firstDau = -1;
    vtxXgen_firstV0Dau = -1.; vtxYgen_firstV0Dau = -1.; vtxZgen_firstV0Dau = -1;
    prodVtxXgen  = -1.0f; prodVtxYgen  = -1.0f; prodVtxZgen  = -1.0f;
    source = 0;

  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, aod::V0sLinked const& V0s, soa::Join<aod::V0fCDatas, aod::V0fCCovs> const& V0fCDatas, soa::Join<aod::V0Datas, aod::V0Covs> const& V0Datas, CascadesCrossLinked const& Cascades, soa::Join<aod::CascDatas, aod::CascCovs> const& CascDatas, soa::Join<aod::KFCascDatas, aod::KFCascCovs> const& KFCascDatas, FullTracksIU const&)
  {
    /// Event selection
    histos.fill(HIST("hEventSelectionFlow"), 1.f);
    // select collisions in acceptance
    if (!(abs(collision.posZ()) < 10.)) return;
    histos.fill(HIST("hEventSelectionFlow"), 2.f);

    for (auto& cascade : Cascades) { // allows for cross-referencing everything

      resetVars();

      // get charge from bachelor (unambiguous wrt to building)
      auto bachTrack = cascade.bachelor_as<FullTracksIU>();
      if (bachTrack.sign() < 0) {
        charge = -1;
      } else {
        charge = +1;
      }
      histos.fill(HIST("hChargeCounter"), charge);

      // store daughter momenta and uncertainties at IU
      auto v0 = cascade.v0_as<aod::V0sLinked>();
      auto posTrack = v0.posTrack_as<FullTracksIU>();
      auto negTrack = v0.negTrack_as<FullTracksIU>();
      o2::track::TrackParCov posTrackParCov = getTrackParCov(posTrack);
      o2::track::TrackParCov negTrackParCov = getTrackParCov(negTrack);
      posTrackParCov.getPxPyPzGlo(momPosRecIU);
      negTrackParCov.getPxPyPzGlo(momNegRecIU);
      std::array<float, 21> cvposini, cvnegini;
      posTrackParCov.getCovXYZPxPyPzGlo(cvposini);
      negTrackParCov.getCovXYZPxPyPzGlo(cvnegini);
      momPosRecIUErr[0] = cvposini[9];
      momPosRecIUErr[1] = cvposini[14];
      momPosRecIUErr[2] = cvposini[21];
      momNegRecIUErr[0] = cvnegini[9];
      momNegRecIUErr[1] = cvnegini[14];
      momNegRecIUErr[2] = cvnegini[21];

      // get cascade data and fill table
      getCascDatas(collision, cascade, CascDatas, KFCascDatas, V0s, V0Datas, V0fCDatas);
      if (cascade.has_cascData() || cascade.has_kfCascData()) {
        fillCascDataTable(collision);
      }
    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processData, "process data", true);

  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, aod::V0sLinked const& V0s, soa::Join<aod::V0Datas, aod::V0Covs> const& V0Datas, soa::Join<aod::V0fCDatas, aod::V0fCCovs> const& V0fCDatas, CascadesCrossLinked const& Cascades, CascDataLabeled const& CascDatas, KFCascDataLabeled const& KFCascDatas, FullTracksIU const&, aod::McParticles const& particlesMC)
  {
    /// Event selection
    histos.fill(HIST("hEventSelectionFlow"), 1.f);
    // select collisions in acceptance
    if (!(abs(collision.posZ()) < 10.)) return;
    histos.fill(HIST("hEventSelectionFlow"), 2.f);

    for (auto& cascade : Cascades) {

      resetVars();

      // get charge from bachelor (unambiguous wrt to building)
      auto bachTrack = cascade.bachelor_as<FullTracksIU>();
      if (bachTrack.sign() < 0) {
        charge = -1;
      } else {
        charge = +1;
      }
      histos.fill(HIST("hChargeCounter"), charge);

      // store daughter momenta and uncertainties at IU
      auto v0 = cascade.v0_as<aod::V0sLinked>();
      auto posTrack = v0.posTrack_as<FullTracksIU>();
      auto negTrack = v0.negTrack_as<FullTracksIU>();
      o2::track::TrackParCov posTrackParCov = getTrackParCov(posTrack);
      o2::track::TrackParCov negTrackParCov = getTrackParCov(negTrack);
      posTrackParCov.getPxPyPzGlo(momPosRecIU);
      negTrackParCov.getPxPyPzGlo(momNegRecIU);
      std::array<float, 21> cvpos, cvneg;
      posTrackParCov.getCovXYZPxPyPzGlo(cvpos);
      negTrackParCov.getCovXYZPxPyPzGlo(cvneg);
      momPosRecIUErr[0] = sqrt(cvpos[9]);
      momPosRecIUErr[1] = sqrt(cvpos[14]);
      momPosRecIUErr[2] = sqrt(cvpos[21]);
      momNegRecIUErr[0] = sqrt(cvneg[9]);
      momNegRecIUErr[1] = sqrt(cvneg[14]);
      momNegRecIUErr[2] = sqrt(cvneg[21]);

      // get cascade data
      getCascDatas(collision, cascade, CascDatas, KFCascDatas, V0s, V0Datas, V0fCDatas);

      // ========== get cascade MC information ===========
      if (cascade.has_kfCascData() && cascade.has_cascData()) {
        LOG(info) << "Both fitters were successful!";
        recocase = 1;
        auto cascdata = cascade.cascData_as<CascDataLabeled>();
        getCascMCdata(collision, cascdata, V0s, particlesMC);
      }
      if (cascade.has_kfCascData() && !cascade.has_cascData()) {
        LOG(info) << "Only KF was successful!";
        recocase = 2;
        auto cascdata = cascade.kfCascData_as<KFCascDataLabeled>();
        getCascMCdata(collision, cascdata, V0s, particlesMC);
      }
      if (!cascade.has_kfCascData() && cascade.has_cascData()) {
        LOG(info) << "Only DCA fitter was successful!";
        recocase = 3;
        auto cascdata = cascade.cascData_as<CascDataLabeled>();
        getCascMCdata(collision, cascdata, V0s, particlesMC);
      }

    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processMC, "process MC", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kfStrangenessStudy>(cfgc)};
}
