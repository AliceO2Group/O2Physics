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

///
/// \file   qaLamMomResolution.cxx
/// \author Carolina Reetz c.reetz@cern.ch
/// \brief  QA task to study momentum resolution of Lambda daughter tracks

#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DPG/Tasks/AOTTrack/V0Cascades/qaLamMomResolution.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Tools/TrackTuner.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using V0DatasLabeled = soa::Join<aod::V0Datas, aod::V0Covs, aod::V0DauCovs, aod::McV0Labels>;
using MyTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::McTrackLabels>;

struct qaLamMomResolution {

  Produces<aod::LamDaughters> lamdaughters;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  OutputObj<TH1D> trackTunedTracks{TH1D("trackTunedTracks", "", 4, 0.5, 4.5), OutputObjHandlingPolicy::AnalysisObject};

  TrackTuner trackTuner;

  Configurable<bool> collSelection{"collSelection", true, "Apply collision selection"};
  Configurable<bool> useTrackTuner{"useTrackTuner", false, "Apply pT/DCA corrections to MC"};
  Configurable<std::string> trackTunerParams{"trackTunerParams", "debugInfo=0|updateTrackDCAs=0|updateTrackCovMat=0|updateCurvature=0|updateCurvatureIU=0|updatePulls=0|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/PbPb2022|nameInputFile=trackTuner_DataLHC22sPass5_McLHC22l1b2_run529397.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Configurable<int> itsAfterburnerPos{"itsAfterburnerPos", 0, "Flag for the ITS afterburner tracks on positive daughters: -1 no AB, 0 no selection, 1 AB"};
  Configurable<int> itsAfterburnerNeg{"itsAfterburnerNeg", 0, "Flag for the ITS afterburner tracks on negative daughters: -1 no AB, 0 no selection, 1 AB"};

  HistogramRegistry hist{"Histograms"};

  o2::dataformats::VertexBase mVtx;
  o2::dataformats::VertexBase mPV;
  o2::dataformats::DCA mDcaInfoCovPos;
  o2::dataformats::DCA mDcaInfoCovNeg;
  o2::track::TrackParametrizationWithError<float> mTrackParCovPosVtx;
  o2::track::TrackParametrizationWithError<float> mTrackParCovNegVtx;
  o2::track::TrackParametrizationWithError<float> mTrackParCovPosPV;
  o2::track::TrackParametrizationWithError<float> mTrackParCovNegPV;

  int runNumber = 0;

  int LambdaPDG = 3122;
  int AntiLambdaPDG = -3122;
  int ProtonPDG = 2212;
  int AntiProtonPDG = -2212;
  int PosPionPDG = 211;
  int NegPionPDG = -211;

  float massLambda = -1.0f;
  float radiusLambda = -1.0f;
  float ptLambda = -1.0f;
  float etaProton = -1.0f, etaPion = -1.0f;
  float phiProton = -1.0f, phiPion = -1.0f;
  int tpcNClsProton = 0, tpcNClsPion = 0;
  int chargeProton = 0, chargePion = 0;
  // daughter momenta
  std::array<float, 3> momProtonRecIU;
  std::array<float, 3> momPionRecIU;
  std::array<float, 3> momProtonRecIUErr;
  std::array<float, 3> momPionRecIUErr;
  std::array<float, 3> momProtonRec;
  std::array<float, 3> momPionRec;
  std::array<float, 3> momProtonRecErr;
  std::array<float, 3> momPionRecErr;
  float sigma1PtProtonIU = -1.0f, sigma1PtPionIU = -1.0f;
  // daughter IU position
  std::array<float, 3> posProtonRecIU;
  std::array<float, 3> posPionRecIU;
  std::array<float, 3> posProtonRecIUErr;
  std::array<float, 3> posPionRecIUErr;
  // daughter DCA
  std::array<float, 2> DCAProtonRec;    // 0: xy, 1: z // updated if tuner is used!
  std::array<float, 2> DCAPionRec;      // 0: xy, 1: z // updated if tuner is used!
  std::array<float, 2> DCAProtonRecErr; // 0: xy, 1: z // updated if tuner is used!
  std::array<float, 2> DCAPionRecErr;   // 0: xy, 1: z // updated if tuner is used!
  // MC info
  std::array<float, 3> momProtonGen;
  std::array<float, 3> momPionGen;
  std::array<float, 2> DCAProtonGen; // 0: xy, 1: z
  std::array<float, 2> DCAPionGen;   // 0: xy, 1: z

  void init(InitContext const&)
  {
    auto hEventSelectionFlow = hist.add<TH1>("hEventSelectionFlow", "Event selection flow", kTH1F, {{2, 0.5f, 2.5f}});
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(1), "Sel8");
    hEventSelectionFlow->GetXaxis()->SetBinLabel(hEventSelectionFlow->FindBin(2), "|Vtx_{z}|<10cm");

    if (useTrackTuner) {
      ccdb->setURL(ccdburl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      std::string outputStringParams = trackTuner.configParams(trackTunerParams);
      trackTuner.getDcaGraphs();
      // QA is done in tuneTrackParams method
      trackTunedTracks->SetTitle(outputStringParams.c_str());
      trackTunedTracks->GetXaxis()->SetBinLabel(1, "all tracks");
      trackTunedTracks->GetXaxis()->SetBinLabel(2, "tracks tuned (no negative detXY)");
      trackTunedTracks->GetXaxis()->SetBinLabel(3, "untouched tracks due to negative detXY");
      trackTunedTracks->GetXaxis()->SetBinLabel(4, "original detXY<0");
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath;
    }
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current() << " A for run " << bc.runNumber() << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    runNumber = bc.runNumber();
  }

  template <typename TTrack>
  bool selectAfterburner(TTrack const& posTrack, TTrack const& negTrack)
  {
    switch (itsAfterburnerPos) {
      case -1:
        if (posTrack.itsChi2NCl() >= 0) {
          return false;
        }
        break;
      case 0:
        break;
      case 1:
        if (posTrack.itsChi2NCl() < 0) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid AB selection for positive daughter";
        break;
    }
    switch (itsAfterburnerNeg) {
      case -1:
        if (negTrack.itsChi2NCl() >= 0) {
          return false;
        }
        break;
      case 0:
        break;
      case 1:
        if (negTrack.itsChi2NCl() < 0) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid AB selection for negative daughter";
        break;
    }
    return true;
  }

  template <typename TV0, typename TV0Track, typename TCollision>
  void tuneV0(TV0 const& v0,
              TV0Track const& posTrack,
              TV0Track const& negTrack,
              aod::McParticles const&,
              TCollision const& collision,
              aod::BCsWithTimestamps const& bcs)
  {
    initCCDB(bcs.begin());
    trackTunedTracks->Fill(1, 2); // tune 2 tracks
    o2::track::TrackParametrizationWithError<float> mTrackParCovPos;
    o2::track::TrackParametrizationWithError<float> mTrackParCovNeg;
    setTrackParCov(posTrack, mTrackParCovPos);
    setTrackParCov(negTrack, mTrackParCovNeg);
    mTrackParCovPos.setPID(posTrack.pidForTracking());
    mTrackParCovNeg.setPID(negTrack.pidForTracking());
    mDcaInfoCovPos.set(999, 999, 999, 999, 999);
    mDcaInfoCovNeg.set(999, 999, 999, 999, 999);
    auto mcParticlePos = posTrack.mcParticle();
    auto mcParticleNeg = negTrack.mcParticle();

    // tune parameters at IU
    trackTuner.tuneTrackParams(mcParticlePos, mTrackParCovPos, matCorr, &mDcaInfoCovPos, trackTunedTracks);
    trackTuner.tuneTrackParams(mcParticleNeg, mTrackParCovNeg, matCorr, &mDcaInfoCovNeg, trackTunedTracks);
    // propagate tuned tracks to Lambda vertex
    mTrackParCovPosVtx = mTrackParCovPos;
    mTrackParCovNegVtx = mTrackParCovNeg;
    mVtx.setPos({v0.x(), v0.y(), v0.z()});
    mVtx.setCov(v0.positionCovMat()[0], v0.positionCovMat()[1], v0.positionCovMat()[2], v0.positionCovMat()[3], v0.positionCovMat()[4], v0.positionCovMat()[5]);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCovPosVtx, 2.f, matCorr, &mDcaInfoCovPos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, mTrackParCovNegVtx, 2.f, matCorr, &mDcaInfoCovNeg);

    // get DCAs to PV
    // DCA with respect to collision associated to the V0
    mTrackParCovPosPV = mTrackParCovPos;
    mTrackParCovNegPV = mTrackParCovNeg;
    mPV.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mPV.setCov(collision.covXX(), collision.covXX(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, mTrackParCovPosPV, 2.f, matCorr, &mDcaInfoCovPos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mPV, mTrackParCovNegPV, 2.f, matCorr, &mDcaInfoCovNeg);
  }

  template <typename TCollision>
  void fillTable(TCollision const& collision)
  {
    lamdaughters(collision.globalIndex(),
                 massLambda, radiusLambda, ptLambda,
                 chargeProton, chargePion,
                 etaProton, etaPion,
                 phiProton, phiPion,
                 tpcNClsProton, tpcNClsPion,
                 momProtonRec[0], momProtonRec[1], momProtonRec[2],
                 momProtonRecErr[0], momProtonRecErr[1], momProtonRecErr[2],
                 momPionRec[0], momPionRec[1], momPionRec[2],
                 momPionRecErr[0], momPionRecErr[1], momPionRecErr[2],
                 momProtonRecIU[0], momProtonRecIU[1], momProtonRecIU[2],
                 momProtonRecIUErr[0], momProtonRecIUErr[1], momProtonRecIUErr[2],
                 momPionRecIU[0], momPionRecIU[1], momPionRecIU[2],
                 momPionRecIUErr[0], momPionRecIUErr[1], momPionRecIUErr[2],
                 momProtonGen[0], momProtonGen[1], momProtonGen[2],
                 momPionGen[0], momPionGen[1], momPionGen[2],
                 sigma1PtProtonIU, sigma1PtPionIU,
                 posProtonRecIU[0], posProtonRecIU[1], posProtonRecIU[2],
                 posProtonRecIUErr[0], posProtonRecIUErr[1], posProtonRecIUErr[2],
                 posPionRecIU[0], posPionRecIU[1], posPionRecIU[2],
                 posPionRecIUErr[0], posPionRecIUErr[1], posPionRecIUErr[2],
                 DCAProtonRec[0], DCAProtonRec[1],
                 DCAProtonRecErr[0], DCAProtonRecErr[1],
                 DCAPionRec[0], DCAPionRec[1],
                 DCAPionRecErr[0], DCAPionRecErr[1]);
  }

  template <typename TProton, typename TPion>
  void getTrackInfo(TProton const& protonTrackIU, TPion const& pionTrackIU)
  {
    // daughter momenta at IU
    o2::track::TrackParCov protonTrackParCov, pionTrackParCov;
    std::array<float, 21> protoncv, pioncv;
    protonTrackParCov = getTrackParCov(protonTrackIU);
    pionTrackParCov = getTrackParCov(pionTrackIU);
    protonTrackParCov.getCovXYZPxPyPzGlo(protoncv);
    pionTrackParCov.getCovXYZPxPyPzGlo(pioncv);
    // proton
    momProtonRecIU[0] = protonTrackIU.px();
    momProtonRecIU[1] = protonTrackIU.py();
    momProtonRecIU[2] = protonTrackIU.pz();
    momProtonRecIUErr[0] = sqrt(protoncv[9]);
    momProtonRecIUErr[1] = sqrt(protoncv[14]);
    momProtonRecIUErr[2] = sqrt(protoncv[20]);
    sigma1PtProtonIU = protonTrackIU.sigma1Pt();
    // pion
    momPionRecIU[0] = pionTrackIU.px();
    momPionRecIU[1] = pionTrackIU.py();
    momPionRecIU[2] = pionTrackIU.pz();
    momPionRecIUErr[0] = sqrt(pioncv[9]);
    momPionRecIUErr[1] = sqrt(pioncv[14]);
    momPionRecIUErr[2] = sqrt(pioncv[20]);
    sigma1PtPionIU = pionTrackIU.sigma1Pt();

    // daughter position at IU
    // proton
    posProtonRecIU[0] = protonTrackIU.x();
    posProtonRecIU[1] = protonTrackIU.y();
    posProtonRecIU[2] = protonTrackIU.z();
    posProtonRecIUErr[0] = sqrt(protoncv[0]);
    posProtonRecIUErr[1] = sqrt(protoncv[2]);
    posProtonRecIUErr[2] = sqrt(protoncv[5]);
    // pion
    posPionRecIU[0] = pionTrackIU.px();
    posPionRecIU[1] = pionTrackIU.py();
    posPionRecIU[2] = pionTrackIU.pz();
    posPionRecIUErr[0] = sqrt(pioncv[0]);
    posPionRecIUErr[1] = sqrt(pioncv[2]);
    posPionRecIUErr[2] = sqrt(pioncv[5]);

    // daughter DCA to PV
    DCAProtonRec[0] = protonTrackIU.dcaXY();
    DCAProtonRec[1] = protonTrackIU.dcaZ();
    DCAProtonRecErr[0] = sqrt(protonTrackIU.sigmaDcaXY2());
    DCAProtonRecErr[1] = sqrt(protonTrackIU.sigmaDcaZ2());
    DCAPionRec[0] = pionTrackIU.dcaXY();
    DCAPionRec[1] = pionTrackIU.dcaZ();
    DCAPionRecErr[0] = sqrt(pionTrackIU.sigmaDcaXY2());
    DCAPionRecErr[1] = sqrt(pionTrackIU.sigmaDcaZ2());

    // daughter charges, eta, nTPCclusters
    chargeProton = protonTrackIU.sign();
    chargePion = pionTrackIU.sign();
    etaProton = protonTrackIU.eta();
    etaPion = pionTrackIU.eta();
    phiProton = protonTrackIU.phi();
    phiPion = pionTrackIU.phi();
    tpcNClsProton = protonTrackIU.tpcNClsFound();
    tpcNClsPion = pionTrackIU.tpcNClsFound();
  }

  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                 V0DatasLabeled const& V0Datas,
                 aod::McParticles const& mcparticles,
                 MyTracks const&,
                 aod::BCsWithTimestamps const& bcs)
  {

    // event selection
    if (collSelection && !collision.sel8()) {
      return;
    }
    hist.fill(HIST("hEventSelectionFlow"), 1.f);
    if (collSelection && (std::abs(collision.posZ()) >= 10.)) {
      return;
    }
    hist.fill(HIST("hEventSelectionFlow"), 2.f);

    for (auto& v0data : V0Datas) {

      if (v0data.has_mcParticle() && v0data.mcParticleId() > -1 && v0data.mcParticleId() <= mcparticles.size()) {
        auto MCv0 = v0data.mcParticle_as<aod::McParticles>();

        // Lambda
        if (MCv0.has_daughters() && MCv0.pdgCode() == LambdaPDG) {
          LOG(debug) << "V0 is a Lambda.";
          const auto& protonTrackIU = v0data.posTrack_as<MyTracks>();
          const auto& pionTrackIU = v0data.negTrack_as<MyTracks>();

          // afterburner selection
          if (!selectAfterburner(protonTrackIU, pionTrackIU))
            continue;

          if (protonTrackIU.has_mcParticle() && pionTrackIU.has_mcParticle()) {

            const auto& MCproton = protonTrackIU.mcParticle_as<aod::McParticles>();
            const auto& MCpion = pionTrackIU.mcParticle_as<aod::McParticles>();

            if (MCproton.pdgCode() == ProtonPDG && MCpion.pdgCode() == NegPionPDG) {

              // lambda
              massLambda = v0data.mLambda();
              radiusLambda = v0data.v0radius();
              ptLambda = v0data.pt();
              /// daughter momenta at Lambda vertex
              // proton
              momProtonRec[0] = v0data.pxpos();
              momProtonRec[1] = v0data.pypos();
              momProtonRec[2] = v0data.pzpos();
              momProtonRecErr[0] = sqrt(v0data.covMatPosDau()[9]);
              momProtonRecErr[1] = sqrt(v0data.covMatPosDau()[14]);
              momProtonRecErr[2] = sqrt(v0data.covMatPosDau()[20]);
              momProtonGen[0] = MCproton.px();
              momProtonGen[1] = MCproton.py();
              momProtonGen[2] = MCproton.pz();
              // pion
              momPionRec[0] = v0data.pxneg();
              momPionRec[1] = v0data.pyneg();
              momPionRec[2] = v0data.pzneg();
              momPionRecErr[0] = sqrt(v0data.covMatNegDau()[9]);
              momPionRecErr[1] = sqrt(v0data.covMatNegDau()[14]);
              momPionRecErr[2] = sqrt(v0data.covMatNegDau()[20]);
              momPionGen[0] = MCpion.px();
              momPionGen[1] = MCpion.py();
              momPionGen[2] = MCpion.pz();

              // get daughter momenta at IU, charge, eta, nTPCclusters
              getTrackInfo(protonTrackIU, pionTrackIU);

              // optionally use track tuner to tune daughter tracks' pT at decay vertex
              if (useTrackTuner) {
                // tune V0
                tuneV0(v0data, protonTrackIU, pionTrackIU, mcparticles, collision, bcs);
                // get smeared parameters and cov matrix
                std::array<float, 3> pPos{0., 0., 0.};
                std::array<float, 3> pNeg{0., 0., 0.};
                std::array<float, 21> cPos, cNeg;
                mTrackParCovPosVtx.getPxPyPzGlo(pPos);
                mTrackParCovNegVtx.getPxPyPzGlo(pNeg);
                mTrackParCovPosVtx.getCovXYZPxPyPzGlo(cPos);
                mTrackParCovNegVtx.getCovXYZPxPyPzGlo(cNeg);
                // lambda
                massLambda = RecoDecay::m(std::array{std::array{pPos[0], pPos[1], pPos[2]},
                                                     std::array{pNeg[0], pNeg[1], pNeg[2]}},
                                          std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
                /// smeared daughter momenta at Lambda vertex
                // proton
                momProtonRec[0] = pPos[0];
                momProtonRec[1] = pPos[1];
                momProtonRec[2] = pPos[2];
                momProtonRecErr[0] = sqrt(cPos[9]);
                momProtonRecErr[1] = sqrt(cPos[14]);
                momProtonRecErr[2] = sqrt(cPos[20]);
                // pion
                momPionRec[0] = pNeg[0];
                momPionRec[1] = pNeg[1];
                momPionRec[2] = pNeg[2];
                momPionRecErr[0] = sqrt(cNeg[9]);
                momPionRecErr[1] = sqrt(cNeg[14]);
                momPionRecErr[2] = sqrt(cNeg[20]);
                // smeared DCAs at PV
                // proton
                DCAProtonRec[0] = mDcaInfoCovPos.getY();
                DCAProtonRec[1] = mDcaInfoCovPos.getZ();
                DCAProtonRecErr[0] = sqrt(mDcaInfoCovPos.getSigmaY2());
                DCAProtonRecErr[1] = sqrt(mDcaInfoCovPos.getSigmaZ2());
                // pion
                DCAPionRec[0] = mDcaInfoCovNeg.getY();
                DCAPionRec[1] = mDcaInfoCovNeg.getZ();
                DCAPionRecErr[0] = sqrt(mDcaInfoCovNeg.getSigmaY2());
                DCAPionRecErr[1] = sqrt(mDcaInfoCovNeg.getSigmaZ2());
              }

              // fill table
              fillTable(collision);
            }
          }
        } // end Lambda

        // Anti-Lambda
        if (MCv0.pdgCode() == AntiLambdaPDG) {
          LOG(debug) << "V0 is an Anti-Lambda.";
          const auto& protonTrackIU = v0data.negTrack_as<MyTracks>();
          const auto& pionTrackIU = v0data.posTrack_as<MyTracks>();

          // afterburner selection
          if (!selectAfterburner(pionTrackIU, protonTrackIU))
            continue;

          if (protonTrackIU.has_mcParticle() && pionTrackIU.has_mcParticle()) {

            const auto& MCproton = protonTrackIU.mcParticle_as<aod::McParticles>();
            const auto& MCpion = pionTrackIU.mcParticle_as<aod::McParticles>();

            if (MCproton.pdgCode() == AntiProtonPDG && MCpion.pdgCode() == PosPionPDG) {

              // lambda mass and radius
              massLambda = v0data.mAntiLambda();
              radiusLambda = v0data.v0radius();
              ptLambda = v0data.pt();
              /// daughter momenta at Lambda vertex
              // proton
              momProtonRec[0] = v0data.pxneg();
              momProtonRec[1] = v0data.pyneg();
              momProtonRec[2] = v0data.pzneg();
              momProtonRecErr[0] = sqrt(v0data.covMatNegDau()[9]);
              momProtonRecErr[1] = sqrt(v0data.covMatNegDau()[14]);
              momProtonRecErr[2] = sqrt(v0data.covMatNegDau()[20]);
              momProtonGen[0] = MCproton.px();
              momProtonGen[1] = MCproton.py();
              momProtonGen[2] = MCproton.pz();
              // pion
              momPionRec[0] = v0data.pxpos();
              momPionRec[1] = v0data.pypos();
              momPionRec[2] = v0data.pzpos();
              momPionRecErr[0] = sqrt(v0data.covMatPosDau()[9]);
              momPionRecErr[1] = sqrt(v0data.covMatPosDau()[14]);
              momPionRecErr[2] = sqrt(v0data.covMatPosDau()[20]);
              momPionGen[0] = MCpion.px();
              momPionGen[1] = MCpion.py();
              momPionGen[2] = MCpion.pz();

              // get daughter momenta at IU, charge, eta, nTPCclusters
              getTrackInfo(protonTrackIU, pionTrackIU);

              // optionally use track tuner to tune daughter tracks' pT at decay vertex
              if (useTrackTuner) {
                // tune V0
                tuneV0(v0data, pionTrackIU, protonTrackIU, mcparticles, collision, bcs);
                // get smeared parameters and cov matrix
                std::array<float, 3> pPos{0., 0., 0.};
                std::array<float, 3> pNeg{0., 0., 0.};
                std::array<float, 21> cPos, cNeg;
                mTrackParCovPosVtx.getPxPyPzGlo(pPos);
                mTrackParCovNegVtx.getPxPyPzGlo(pNeg);
                mTrackParCovPosVtx.getCovXYZPxPyPzGlo(cPos);
                mTrackParCovNegVtx.getCovXYZPxPyPzGlo(cNeg);
                // lambda
                massLambda = RecoDecay::m(std::array{std::array{pPos[0], pPos[1], pPos[2]},
                                                     std::array{pNeg[0], pNeg[1], pNeg[2]}},
                                          std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                /// smeared daughter momenta at Lambda vertex
                // proton
                momProtonRec[0] = pNeg[0];
                momProtonRec[1] = pNeg[1];
                momProtonRec[2] = pNeg[2];
                momProtonRecErr[0] = sqrt(cNeg[9]);
                momProtonRecErr[1] = sqrt(cNeg[14]);
                momProtonRecErr[2] = sqrt(cNeg[20]);
                // pion
                momPionRec[0] = pPos[0];
                momPionRec[1] = pPos[1];
                momPionRec[2] = pPos[2];
                momPionRecErr[0] = sqrt(cPos[9]);
                momPionRecErr[1] = sqrt(cPos[14]);
                momPionRecErr[2] = sqrt(cPos[20]);
                // smeared DCAs at PV
                // proton
                DCAProtonRec[0] = mDcaInfoCovNeg.getY();
                DCAProtonRec[1] = mDcaInfoCovNeg.getZ();
                DCAProtonRecErr[0] = sqrt(mDcaInfoCovNeg.getSigmaY2());
                DCAProtonRecErr[1] = sqrt(mDcaInfoCovNeg.getSigmaZ2());
                // pion
                DCAPionRec[0] = mDcaInfoCovPos.getY();
                DCAPionRec[1] = mDcaInfoCovPos.getZ();
                DCAPionRecErr[0] = sqrt(mDcaInfoCovPos.getSigmaY2());
                DCAPionRecErr[1] = sqrt(mDcaInfoCovPos.getSigmaZ2());
              }

              // fill table
              fillTable(collision);
            }
          }
        } // end Anti-Lambda
      }   // end MC
    }     // end V0 loop
  }
  PROCESS_SWITCH(qaLamMomResolution, processMC, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaLamMomResolution>(cfgc)};
}
