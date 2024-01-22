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
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
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
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksCovIU, o2::aod::TracksExtra>;

struct kfStrangenessStudy {

  Produces<aod::CascCand> rowCasc;
  Produces<aod::CascCandMC> rowCascMC;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Filters
  Filter collisionFilter = (aod::evsel::sel8 == true);

  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{0};
  double bz = 0.;

  // magnetic field settings for CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  /// Cascade data
  int isDCAfitter = 0, isKF = 0;
  float ptRec = -1.0f, ptRecKF = -1.0f;
  float ptRecV0 = -1.0f, ptRecKFV0 = -1.0f;
  float massXi = -1.0f, massXiKF = -1.0f;
  float massLambda = -1.0f, massLambdaKF = -1.0f;
  float cascRad = -1.0f, cascRadKF = -1.0f;
  std::array<float, 3> vtxRec, vtxRecKF;
  std::array<float, 3> vtxRecErr, vtxRecErrKF;
  std::array<float, 3> vtxRecV0, vtxRecKFV0;
  std::array<float, 3> vtxRecErrV0, vtxRecErrKFV0;
  float dcaXYCascToPV = -1.0f, dcaXYCascToPVKF = -1.0f;
  float dcaZCascToPV = -1.0f, dcaZCascToPVKF = -1.0f;
  float dcaCascDaughters = -1.0f, dcaCascDaughtersKF = -1.0f;
  float dcaV0Daughters = -1.0f, dcaV0DaughtersKF = -1.0f;
  float dcaProtonToPV = -1.0f, dcaProtonToPVKF = -1.0f;
  float dcaPionToPV = -1.0f, dcaPionToPVKF = -1.0f;
  float dcaBachToPV = -1.0f, dcaBachToPVKF = -1.0f;
  float cascPointingAngle = -1.0f, cascPointingAngleKF = -1.0f;
  float v0PointingAngle = -1.0f, v0PointingAngleKF = -1.0f;
  float V0Rad = -1.0f, V0RadKF = -1.0f;
  int charge = 0;
  float etaProton = -1.0, etaPion = -1.0;
  int tpcNClsProton = 0, tpcNClsPion = 0;
  std::array<float, 3> momProtonRecIU;
  std::array<float, 3> momPionRecIU;
  std::array<float, 3> momProtonRecIUErr;
  std::array<float, 3> momPionRecIUErr;
  std::array<float, 3> momProtonRec;
  std::array<float, 3> momPionRec;
  std::array<float, 3> momProtonRecErr;
  std::array<float, 3> momPionRecErr;
  std::array<float, 3> posProtonRec;
  std::array<float, 3> posPionRec;
  std::array<float, 3> posProtonRecErr;
  std::array<float, 3> posPionRecErr;

  /// Additional cascade MC data
  int isTrueCasc = 0;
  float ptGen = -1.0f;
  float ptGenV0 = -1.0f;
  std::array<float, 3> vtxGen;
  std::array<float, 3> vtxGenV0;
  std::array<float, 3> prodVtxGen;
  int source = 0;
  std::array<float, 3> momProtonGen;
  std::array<float, 3> momPionGen;

  // counter and checks
  int recocase = 0;

  void init(InitContext const&)
  {

    /// CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

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

  void initCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
                o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string const& ccdbPathGrp, o2::base::MatLayerCylSet* lut)
  {
    if (mRunNumber != bc.runNumber()) {
      LOGF(info, "====== initCCDB function called");
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
    
      mRunNumber = bc.runNumber();
    }
  } /// end initCCDB

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
        vtxRecErrV0[0] = sqrt(v0data.positionCovMat()[0]);
        vtxRecErrV0[1] = sqrt(v0data.positionCovMat()[2]);
        vtxRecErrV0[2] = sqrt(v0data.positionCovMat()[5]);
      } else if (v0index.has_v0fCData()) {
        // V0 passed only cascade V0 selections, use this instead
        auto v0data = v0index.template v0fCData_as<TV0fCDatas>();
        vtxRecErrV0[0] = sqrt(v0data.positionCovMat()[0]);
        vtxRecErrV0[1] = sqrt(v0data.positionCovMat()[2]);
        vtxRecErrV0[2] = sqrt(v0data.positionCovMat()[5]);
      }
      ptRec = cascdata.pt();
      vtxRec[0] = cascdata.x();
      vtxRec[1] = cascdata.y();
      vtxRec[2] = cascdata.z();
      vtxRecErr[0] = sqrt(cascdata.positionCovMat()[0]);
      vtxRecErr[1] = sqrt(cascdata.positionCovMat()[2]);
      vtxRecErr[2] = sqrt(cascdata.positionCovMat()[5]);
      ptRecV0 = sqrt((cascdata.pxpos()+cascdata.pxneg())*(cascdata.pxpos()+cascdata.pxneg()) + (cascdata.pypos()+cascdata.pyneg())*(cascdata.pypos()+cascdata.pyneg())); // taken from daughters from lambdakzerobuiler
      vtxRecV0[0] = cascdata.xlambda();
      vtxRecV0[1] = cascdata.ylambda();
      vtxRecV0[2] = cascdata.zlambda();
      massLambda = cascdata.mLambda();
      massXi = cascdata.mXi();
      dcaXYCascToPV = cascdata.dcaXYCascToPV();
      dcaZCascToPV = cascdata.dcaZCascToPV();
      dcaCascDaughters = cascdata.dcacascdaughters();
      dcaV0Daughters = cascdata.dcaV0daughters();
      if (charge == -1) {
        dcaProtonToPV = cascdata.dcapostopv();
        dcaPionToPV = cascdata.dcanegtopv();
      } else if (charge == +1) {
        dcaProtonToPV = cascdata.dcanegtopv();
        dcaPionToPV = cascdata.dcapostopv();
      }
      dcaBachToPV = cascdata.dcabachtopv();
      cascPointingAngle = TMath::ACos(cascdata.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRad = cascdata.cascradius();
      V0Rad = cascdata.v0radius();
      v0PointingAngle = TMath::ACos(cascdata.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));

      // fill QA histos
      histos.fill(HIST("hVertexX"), vtxRec[0]);
      histos.fill(HIST("hVertexY"), vtxRec[1]);
      histos.fill(HIST("hVertexZ"), vtxRec[2]);
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
      vtxRecKF[0] = cascdatakf.x();
      vtxRecKF[1] = cascdatakf.y();
      vtxRecKF[2] = cascdatakf.z();
      vtxRecErrKF[0] = sqrt(cascdatakf.kfTrackCovMat()[0]);
      vtxRecErrKF[1] = sqrt(cascdatakf.kfTrackCovMat()[2]);
      vtxRecErrKF[2] = sqrt(cascdatakf.kfTrackCovMat()[5]);
      ptRecKFV0 = sqrt(cascdatakf.pxv0()*cascdatakf.pxv0() + cascdatakf.pyv0()*cascdatakf.pyv0()); // taken from KFV0
      vtxRecKFV0[0] = cascdatakf.xlambda();
      vtxRecKFV0[1] = cascdatakf.ylambda();
      vtxRecKFV0[2] = cascdatakf.zlambda();
      vtxRecErrKFV0[0] = sqrt(cascdatakf.kfTrackCovMatV0()[0]);
      vtxRecErrKFV0[1] = sqrt(cascdatakf.kfTrackCovMatV0()[2]);
      vtxRecErrKFV0[2] = sqrt(cascdatakf.kfTrackCovMatV0()[5]);
      massLambdaKF = cascdatakf.mLambda();
      massXiKF = cascdatakf.mXi();
      dcaXYCascToPVKF = cascdatakf.dcaXYCascToPV();
      dcaZCascToPVKF = cascdatakf.dcaZCascToPV();
      dcaCascDaughtersKF = cascdatakf.dcacascdaughters();
      dcaV0DaughtersKF = cascdatakf.dcaV0daughters();
      if (charge == -1) {
        dcaProtonToPVKF = cascdatakf.dcapostopv();
        dcaPionToPVKF = cascdatakf.dcanegtopv();
      } else if (charge == +1) {
        dcaProtonToPVKF = cascdatakf.dcanegtopv();
        dcaPionToPVKF = cascdatakf.dcapostopv();
      }
      dcaBachToPVKF = cascdatakf.dcabachtopv();
      cascPointingAngleKF = TMath::ACos(cascdatakf.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRadKF = cascdatakf.cascradius();
      V0RadKF = cascdatakf.v0radius();
      v0PointingAngleKF = TMath::ACos(cascdatakf.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
  
      if (charge == -1) {
        // daughter momenta at vertex (from pre-minimisation if enabled)
        momProtonRec[0] = cascdatakf.pxpos();
        momProtonRec[1] = cascdatakf.pypos();
        momProtonRec[2] = cascdatakf.pzpos();
        momPionRec[0] = cascdatakf.pxneg();
        momPionRec[1] = cascdatakf.pyneg();
        momPionRec[2] = cascdatakf.pzneg();
        momProtonRecErr[0] = cascdatakf.pxposerr();
        momProtonRecErr[1] = cascdatakf.pyposerr();
        momProtonRecErr[2] = cascdatakf.pzposerr();
        momPionRecErr[0] = cascdatakf.pxnegerr();
        momPionRecErr[1] = cascdatakf.pynegerr();
        momPionRecErr[2] = cascdatakf.pznegerr();
        // daughter track position at vertex (after pre-minimisation!)
        // posProtonRec[0] = cascdatakf.xpos();
        // posProtonRec[1] = cascdatakf.ypos();
        // posProtonRec[2] = cascdatakf.zpos();
        // posPionRec[0] = cascdatakf.xneg();
        // posPionRec[1] = cascdatakf.yneg();
        // posPionRec[2] = cascdatakf.zneg();
        // posProtonRecErr[0] = cascdatakf.xposerr();
        // posProtonRecErr[1] = cascdatakf.yposerr();
        // posProtonRecErr[2] = cascdatakf.zposerr();
        // posPionRecErr[0] = cascdatakf.xnegerr();
        // posPionRecErr[1] = cascdatakf.ynegerr();
        // posPionRecErr[2] = cascdatakf.znegerr();
      } else if (charge == +1) {
        // daughter momenta at vertex (from pre-minimisation if enabled)
        momProtonRec[0] = cascdatakf.pxneg();
        momProtonRec[1] = cascdatakf.pyneg();
        momProtonRec[2] = cascdatakf.pzneg();
        momPionRec[0] = cascdatakf.pxpos();
        momPionRec[1] = cascdatakf.pypos();
        momPionRec[2] = cascdatakf.pzpos();
        momProtonRecErr[0] = cascdatakf.pxnegerr();
        momProtonRecErr[1] = cascdatakf.pynegerr();
        momProtonRecErr[2] = cascdatakf.pynegerr();
        momPionRecErr[0] = cascdatakf.pxposerr();
        momPionRecErr[1] = cascdatakf.pyposerr();
        momPionRecErr[2] = cascdatakf.pzposerr();
        // daughter track position at vertex (after pre-minimisation!)
        // posProtonRec[0] = cascdatakf.xneg();
        // posProtonRec[1] = cascdatakf.yneg();
        // posProtonRec[2] = cascdatakf.zneg();
        // posPionRec[0] = cascdatakf.xpos();
        // posPionRec[1] = cascdatakf.ypos();
        // posPionRec[2] = cascdatakf.zpos();
        // posProtonRecErr[0] = cascdatakf.xnegerr();
        // posProtonRecErr[1] = cascdatakf.ynegerr();
        // posProtonRecErr[2] = cascdatakf.znegerr();
        // posPionRecErr[0] = cascdatakf.xposerr();
        // posPionRecErr[1] = cascdatakf.yposerr();
        // posPionRecErr[2] = cascdatakf.zposerr();
      }
      
      // fill QA histos
      histos.fill(HIST("hKFVertexX"), vtxRecKF[0]);
      histos.fill(HIST("hKFVertexY"), vtxRecKF[1]);
      histos.fill(HIST("hKFVertexZ"), vtxRecKF[2]);
      histos.fill(HIST("hKFDCAxy"), dcaXYCascToPVKF);
      histos.fill(HIST("hKFCascRadius"), cascRadKF);
      histos.fill(HIST("hKFPointingAngle"), cascPointingAngleKF);
      histos.fill(HIST("hKFCosPointingAngle"), cos(cascPointingAngleKF));
      histos.fill(HIST("hKFV0PointingAngle"), v0PointingAngleKF);
      histos.fill(HIST("hKFCosV0PointingAngle"), cos(v0PointingAngleKF));
    }
  }

  template <typename TCollision, typename TCascData, typename TFullTracksIU, typename TMCParticle>
  void getCascMCdata(TCollision const& collision, TCascData const& cascdata, TFullTracksIU const&, TMCParticle const& mcparticles)
  {
    if (cascdata.has_mcParticle() && cascdata.mcParticleId() > -1 && cascdata.mcParticleId() <= mcparticles.size()) {
      auto MCcascade = cascdata.template mcParticle_as<TMCParticle>();

      if (MCcascade.has_daughters()) {
        LOG(info) << "MC cascade has daughters, getting MC info.";
        // get MC V0
        auto MCv0 = MCcascade.template daughters_as<TMCParticle>().begin();
        if (abs(MCv0.pdgCode()) == 3122 && MCv0.has_daughters()) {
          // cascade
          ptGen = MCcascade.pt();
          prodVtxGen[0] = MCcascade.vx();
          prodVtxGen[1] = MCcascade.vy();
          prodVtxGen[2] = MCcascade.vz();
          vtxGen[0] = MCv0.vx();
          vtxGen[1] = MCv0.vy();
          vtxGen[2] = MCv0.vz();
          // V0
          ptGenV0 = MCv0.pt();
          vtxGenV0[0] = MCv0.template daughters_as<TMCParticle>().begin().vx(); // MC V0 vertex
          vtxGenV0[1] = MCv0.template daughters_as<TMCParticle>().begin().vy();
          vtxGenV0[2] = MCv0.template daughters_as<TMCParticle>().begin().vz();
          // daughters
          for (auto& d : MCv0.template daughters_as<TMCParticle>()) {
            if (abs(d.pdgCode()) == 2212) {
              momProtonGen[0] = d.px();
              momProtonGen[1] = d.py();
              momProtonGen[2] = d.pz();
            } else if (abs(d.pdgCode()) == 211) {
              momPionGen[0] = d.px();
              momPionGen[1] = d.py();
              momPionGen[2] = d.pz();
            }
          }
          
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

          // reconstructed daughter track position at MC truth V0 vertex
          // get daughter tracks from V0
          auto posTrack = cascdata.template posTrack_as<TFullTracksIU>();
          auto negTrack = cascdata.template negTrack_as<TFullTracksIU>();
          o2::track::TrackParCov posTrackParCov = getTrackParCov(posTrack);
          o2::track::TrackParCov negTrackParCov = getTrackParCov(negTrack);
          // propagate to V0 MC vertex
          posTrackParCov.propagateTo(vtxGenV0[0], bz); // propagate the track to the X closest to the V0 vertex
          negTrackParCov.propagateTo(vtxGenV0[0], bz); // propagate the track to the X closest to the V0 vertex
          // get new daughter track position and uncertainties
          std::array<float, 3> protonxyz, pionxyz;
          std::array<float, 21> protoncv, pioncv;
          if (charge == -1) {
            posTrackParCov.getXYZGlo(protonxyz);
            negTrackParCov.getXYZGlo(pionxyz);
            posTrackParCov.getCovXYZPxPyPzGlo(protoncv);
            negTrackParCov.getCovXYZPxPyPzGlo(pioncv);
          } else if (charge == +1) {
            negTrackParCov.getXYZGlo(protonxyz);
            posTrackParCov.getXYZGlo(pionxyz);
            posTrackParCov.getCovXYZPxPyPzGlo(pioncv);
            negTrackParCov.getCovXYZPxPyPzGlo(protoncv);
          }
          for (int i = 0; i < 3; i++) {
            posProtonRec[i] = protonxyz[i];
            posPionRec[i] = pionxyz[i];
          }
          posProtonRecErr[0] = sqrt(protoncv[0]);
          posProtonRecErr[1] = sqrt(protoncv[2]);
          posProtonRecErr[2] = sqrt(protoncv[5]);
          posPionRecErr[0] = sqrt(pioncv[0]);
          posPionRecErr[1] = sqrt(pioncv[2]);
          posPionRecErr[2] = sqrt(pioncv[5]);
          
          // fill cascade table
          fillCascMCTable(collision);

          // fill QA histos --> vertex position from daughters!
          histos.fill(HIST("hGenDecayVtxX_firstDau"), vtxGen[0]);
          histos.fill(HIST("hGenDecayVtxY_firstDau"), vtxGen[1]);
          histos.fill(HIST("hGenDecayVtxZ_firstDau"), vtxGen[2]);

          histos.fill(HIST("hGenSource"), source);
          histos.fill(HIST("hChargeCounterCascDatas"), charge);

        } else {
          LOG(info) << "V0 is no Lambda and/or has no daughters. V0 PDG code: " << MCv0.pdgCode();
        } // end v0 has daughters and is Lambda
      } // end cascade has daughters
    } // end cascade has MC particle
  }

  template <typename TCollision>
  void fillCascDataTable(TCollision const& collision)
  {
    rowCasc(collision.globalIndex(),
            ptRec, ptRecKF,
            massXi, massXiKF,
            cascRad, cascRadKF,
            vtxRec[0], vtxRec[1], vtxRec[2], vtxRecErr[0], vtxRecErr[1], vtxRecErr[2],
            vtxRecKF[0], vtxRecKF[1], vtxRecKF[2], vtxRecErrKF[0], vtxRecErrKF[1], vtxRecErrKF[2],
            dcaXYCascToPV, dcaXYCascToPVKF,
            dcaZCascToPV, dcaZCascToPVKF,
            dcaCascDaughters, dcaCascDaughtersKF,
            cascPointingAngle, cascPointingAngleKF,
            charge,
            ptRecV0, ptRecKFV0,
            massLambda, massLambdaKF,
            V0Rad, V0RadKF,
            vtxRecV0[0], vtxRecV0[1], vtxRecV0[2], vtxRecErrV0[0], vtxRecErrV0[1], vtxRecErrV0[2],
            vtxRecKFV0[0], vtxRecKFV0[1], vtxRecKFV0[2], vtxRecErrKFV0[0], vtxRecErrKFV0[1], vtxRecErrKFV0[2],
            dcaV0Daughters, dcaV0DaughtersKF,
            dcaProtonToPV, dcaProtonToPVKF,
            dcaPionToPV, dcaPionToPVKF,
            dcaBachToPV, dcaBachToPVKF,
            v0PointingAngle, v0PointingAngleKF,
            momProtonRecIU[0], momProtonRecIU[1], momProtonRecIU[2], momProtonRecIUErr[0], momProtonRecIUErr[1], momProtonRecIUErr[2],
            momPionRecIU[0], momPionRecIU[1], momPionRecIU[2], momPionRecIUErr[0], momPionRecIUErr[1], momPionRecIUErr[2],
            momProtonRec[0], momProtonRec[1], momProtonRec[2], momProtonRecErr[0], momProtonRecErr[1], momProtonRecErr[2],
            momPionRec[0], momPionRec[1], momPionRec[2], momPionRecErr[0], momPionRecErr[1], momPionRecErr[2],
            etaProton, etaPion,
            tpcNClsProton, tpcNClsPion,
            isDCAfitter, isKF);
  }

  template <typename TCollision>
  void fillCascMCTable(TCollision const& collision)
  {
    rowCascMC(collision.globalIndex(),
              ptRec, ptRecKF, ptGen,
              massXi, massXiKF,
              cascRad, cascRadKF,
              vtxRec[0], vtxRec[1], vtxRec[2], vtxRecErr[0], vtxRecErr[1], vtxRecErr[2],
              vtxRecKF[0], vtxRecKF[1], vtxRecKF[2], vtxRecErrKF[0], vtxRecErrKF[1], vtxRecErrKF[2],
              vtxGen[0], vtxGen[1], vtxGen[2],
              prodVtxGen[0], prodVtxGen[1], prodVtxGen[2],
              dcaXYCascToPV, dcaXYCascToPVKF,
              dcaZCascToPV, dcaZCascToPVKF,
              dcaCascDaughters, dcaCascDaughtersKF,
              cascPointingAngle, cascPointingAngleKF,
              charge,
              ptRecV0, ptRecKFV0, ptGenV0,
              massLambda, massLambdaKF,
              V0Rad, V0RadKF,
              vtxRecV0[0], vtxRecV0[1], vtxRecV0[2], vtxRecErrV0[0], vtxRecErrV0[1], vtxRecErrV0[2],
              vtxRecKFV0[0], vtxRecKFV0[1], vtxRecKFV0[2], vtxRecErrKFV0[0], vtxRecErrKFV0[1], vtxRecErrKFV0[2],
              vtxGenV0[0], vtxGenV0[1], vtxGenV0[2],
              dcaV0Daughters, dcaV0DaughtersKF,
              dcaProtonToPV, dcaProtonToPVKF,
              dcaPionToPV, dcaPionToPVKF,
              dcaBachToPV, dcaBachToPVKF,
              v0PointingAngle, v0PointingAngleKF,
              momProtonRecIU[0], momProtonRecIU[1], momProtonRecIU[2], momProtonRecIUErr[0], momProtonRecIUErr[1], momProtonRecIUErr[2],
              momPionRecIU[0], momPionRecIU[1], momPionRecIU[2], momPionRecIUErr[0], momPionRecIUErr[1], momPionRecIUErr[2],
              momProtonRec[0], momProtonRec[1], momProtonRec[2], momProtonRecErr[0], momProtonRecErr[1], momProtonRecErr[2],
              momPionRec[0], momPionRec[1], momPionRec[2], momPionRecErr[0], momPionRecErr[1], momPionRecErr[2],
              momProtonGen[0], momProtonGen[1], momProtonGen[2],
              momPionGen[0], momPionGen[1], momPionGen[2],
              posProtonRec[0], posProtonRec[1], posProtonRec[2], posProtonRecErr[0], posProtonRecErr[1], posProtonRecErr[2],
              posPionRec[0], posPionRec[1], posPionRec[2], posPionRecErr[0], posPionRecErr[1], posPionRecErr[2],
              etaProton, etaPion,
              tpcNClsProton, tpcNClsPion,
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
    dcaXYCascToPV = -1.0f; dcaXYCascToPVKF = -1.0f;
    dcaZCascToPV = -1.0f; dcaZCascToPVKF = -1.0f;
    dcaCascDaughters = -1.0f; dcaCascDaughtersKF = -1.0f;
    dcaV0Daughters = -1.0f; dcaV0DaughtersKF = -1.0f;
    dcaProtonToPV = -1.0f; dcaProtonToPVKF = -1.0f;
    dcaPionToPV = -1.0f; dcaPionToPVKF = -1.0f;
    dcaBachToPV = -1.0f; dcaBachToPVKF = -1.0f;
    cascPointingAngle = -1.0f; cascPointingAngleKF = -1.0f;
    v0PointingAngle = -1.0f; v0PointingAngleKF = -1.0f;
    V0Rad = -1.0f; V0RadKF = -1.0f;
    charge = 0;
    etaProton = -1.0, etaPion = -1.0;
    tpcNClsProton = 0, tpcNClsPion = 0;
    for (int i = 0; i < 3; i++) {
      vtxRec[i] = -1.0f;
      vtxRecKF[i] = -1.0f;
      vtxRecErr[i] = -1.0f;
      vtxRecErrKF[i] = -1.0f;
      vtxRecV0[i] = -1.0f;
      vtxRecKFV0[i] = -1.0f;
      vtxRecErrV0[i] = -1.0f;
      vtxRecErrKFV0[i] = -1.0f;

      momProtonRecIU[i] = -1.0f;
      momPionRecIU[i] = -1.0f;
      momProtonRecIUErr[i] = -1.0f;
      momPionRecIUErr[i] = -1.0f;
      momProtonRec[i] = -1.0f;
      momPionRec[i] = -1.0f;
      momProtonRecErr[i] = -1.0f;
      momPionRecErr[i] = -1.0f;
      posProtonRec[i] = -1.0f;
      posPionRec[i] = -1.0f;
      posProtonRecErr[i] = -1.0f;
      posPionRecErr[i] = -1.0f;
      
      // Additional cascade MC data
      vtxGen[i] = -1.0f;
      vtxGenV0[i] = -1.0f;
      prodVtxGen[i]  = -1.0f;
      momProtonGen[i] = -1.0f;
      momPionGen[i] = -1.0f;
    }

    /// Additional cascade MC data
    isTrueCasc = 0;
    ptGen = -1.0f;
    ptGenV0 = -1.0f;
    source = 0;

  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, 
                  aod::V0sLinked const& V0s, 
                  soa::Join<aod::V0fCDatas, aod::V0fCCovs> const& V0fCDatas, 
                  soa::Join<aod::V0Datas, aod::V0Covs> const& V0Datas, 
                  CascadesCrossLinked const& Cascades, 
                  soa::Join<aod::CascDatas, aod::CascCovs> const& CascDatas, 
                  soa::Join<aod::KFCascDatas, aod::KFCascCovs> const& KFCascDatas, 
                  FullTracksIU const&,
                  aod::BCsWithTimestamps const&)
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

      // store daughter momenta and uncertainties at IU (an eta)
      auto v0 = cascade.v0_as<aod::V0sLinked>(); // soa::Join<o2::aod::V0s, o2::aod::V0DataLink>
      auto posTrack = v0.posTrack_as<FullTracksIU>();
      auto negTrack = v0.negTrack_as<FullTracksIU>();
      o2::track::TrackParCov posTrackParCov = getTrackParCov(posTrack);
      o2::track::TrackParCov negTrackParCov = getTrackParCov(negTrack);
      std::array<float, 21> cvproton, cvpion;
      if (charge == -1) {
        posTrackParCov.getPxPyPzGlo(momProtonRecIU);
        negTrackParCov.getPxPyPzGlo(momPionRecIU);
        posTrackParCov.getCovXYZPxPyPzGlo(cvproton);
        negTrackParCov.getCovXYZPxPyPzGlo(cvpion);
        etaProton = posTrack.eta();
        etaPion = negTrack.eta();
        tpcNClsProton = posTrack.tpcNClsFound();
        tpcNClsPion = negTrack.tpcNClsFound();
      } else if (charge == +1) {
        posTrackParCov.getPxPyPzGlo(momPionRecIU);
        negTrackParCov.getPxPyPzGlo(momProtonRecIU);
        posTrackParCov.getCovXYZPxPyPzGlo(cvproton);
        negTrackParCov.getCovXYZPxPyPzGlo(cvpion);
        etaProton = negTrack.eta();
        etaPion = posTrack.eta();
        tpcNClsProton = negTrack.tpcNClsFound();
        tpcNClsPion = posTrack.tpcNClsFound();
      }
      momProtonRecIUErr[0] = sqrt(cvproton[9]);
      momProtonRecIUErr[1] = sqrt(cvproton[14]);
      momProtonRecIUErr[2] = sqrt(cvproton[21]);
      momPionRecIUErr[0] = sqrt(cvpion[9]);
      momPionRecIUErr[1] = sqrt(cvpion[14]);
      momPionRecIUErr[2] = sqrt(cvpion[21]);
      

      // get cascade data and fill table
      getCascDatas(collision, cascade, CascDatas, KFCascDatas, V0s, V0Datas, V0fCDatas);
      if (cascade.has_cascData() || cascade.has_kfCascData()) {
        fillCascDataTable(collision);
      }
    } // end cascade loop
  } // end process
  PROCESS_SWITCH(kfStrangenessStudy, processData, "process data", true);

  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, 
                aod::V0sLinked const& V0s, 
                soa::Join<aod::V0Datas, aod::V0Covs> const& V0Datas, 
                soa::Join<aod::V0fCDatas, aod::V0fCCovs> const& V0fCDatas, 
                CascadesCrossLinked const& Cascades, 
                CascDataLabeled const& CascDatas, 
                KFCascDataLabeled const& KFCascDatas, 
                FullTracksIU const& TracksIU, 
                aod::McParticles const& particlesMC,
                aod::BCsWithTimestamps const&)
  {
    /// magnetic field from CCDB
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initCCDB(bc, runNumber, ccdb, ccdbPathGrpMag, lut);
      bz = o2::base::Propagator::Instance()->getNominalBz();
    }

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

      // store daughter momenta and uncertainties at IU (and eta)
      auto v0 = cascade.v0_as<aod::V0sLinked>();
      auto posTrack = v0.posTrack_as<FullTracksIU>();
      auto negTrack = v0.negTrack_as<FullTracksIU>();
      o2::track::TrackParCov posTrackParCov = getTrackParCov(posTrack);
      o2::track::TrackParCov negTrackParCov = getTrackParCov(negTrack);
      std::array<float, 21> cvproton, cvpion;
      if (charge == -1) {
        posTrackParCov.getPxPyPzGlo(momProtonRecIU);
        negTrackParCov.getPxPyPzGlo(momPionRecIU);
        posTrackParCov.getCovXYZPxPyPzGlo(cvproton);
        negTrackParCov.getCovXYZPxPyPzGlo(cvpion);
        etaProton = posTrack.eta();
        etaPion = negTrack.eta();
        tpcNClsProton = posTrack.tpcNClsFound();
        tpcNClsPion = negTrack.tpcNClsFound();
      } else if (charge == +1) {
        posTrackParCov.getPxPyPzGlo(momPionRecIU);
        negTrackParCov.getPxPyPzGlo(momProtonRecIU);
        posTrackParCov.getCovXYZPxPyPzGlo(cvproton);
        negTrackParCov.getCovXYZPxPyPzGlo(cvpion);
        etaProton = negTrack.eta();
        etaPion = posTrack.eta();
        tpcNClsProton = negTrack.tpcNClsFound();
        tpcNClsPion = posTrack.tpcNClsFound();
      }
      momProtonRecIUErr[0] = sqrt(cvproton[9]);
      momProtonRecIUErr[1] = sqrt(cvproton[14]);
      momProtonRecIUErr[2] = sqrt(cvproton[21]);
      momPionRecIUErr[0] = sqrt(cvpion[9]);
      momPionRecIUErr[1] = sqrt(cvpion[14]);
      momPionRecIUErr[2] = sqrt(cvpion[21]);

      // get cascade data
      getCascDatas(collision, cascade, CascDatas, KFCascDatas, V0s, V0Datas, V0fCDatas);

      // ========== get cascade MC information ===========
      if (cascade.has_kfCascData() && cascade.has_cascData()) {
        LOG(info) << "Both fitters were successful!";
        recocase = 1;
        auto cascdata = cascade.cascData_as<CascDataLabeled>();
        getCascMCdata(collision, cascdata, TracksIU, particlesMC);
      }
      if (cascade.has_kfCascData() && !cascade.has_cascData()) {
        LOG(info) << "Only KF was successful!";
        recocase = 2;
        auto cascdata = cascade.kfCascData_as<KFCascDataLabeled>();
        getCascMCdata(collision, cascdata, TracksIU, particlesMC);
      }
      if (!cascade.has_kfCascData() && cascade.has_cascData()) {
        LOG(info) << "Only DCA fitter was successful!";
        recocase = 3;
        auto cascdata = cascade.cascData_as<CascDataLabeled>();
        getCascMCdata(collision, cascdata, TracksIU, particlesMC);
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
