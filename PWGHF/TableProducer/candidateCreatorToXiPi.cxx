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

/// \file candidateCreatorToXiPi.cxx
/// \brief Reconstruction of Omegac0 and Xic0 -> xi pi candidates
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::track;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::aod::hf_track_index;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Reconstruction of omegac0 and xic0 candidates
struct HfCandidateCreatorToXiPi {
  Produces<aod::HfCandToXiPi> rowCandidate;

  Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> maxDXYIni{"maxDXYIni", 4., "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};
  Configurable<bool> refitWithMatCorr{"refitWithMatCorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber;

  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>; // to use strangeness tracking, use aod::TraCascDatas instead of aod::CascDatas
  using CascadesLinked = soa::Join<Cascades, CascDataLink>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MySkimIdx = soa::Filtered<aod::HfCascLf2Prongs>;

  Filter filterSelectIndexes = (aod::hf_track_index::hfflag > static_cast<uint8_t>(0));

  OutputObj<TH1F> hInvMassCharmBaryon{TH1F("hInvMassCharmBaryon", "Charm baryon invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hFitterStatus{TH1F("hFitterStatus", "Charm DCAFitter status;status;entries", 3, -0.5, 2.5)};                     // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
  OutputObj<TH1F> hCandidateCounter{TH1F("hCandidateCounter", "Candidate counter wrt derived data;status;entries", 4, -0.5, 3.5)}; // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
  OutputObj<TH1F> hCascadesCounter{TH1F("hCascadesCounter", "Cascades counter wrt derived data;status;entries", 2, -0.5, 1.5)};    // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  void process(aod::Collisions const&,
               aod::BCsWithTimestamps const&,
               TracksWCovDca const&,
               MyCascTable const&, CascadesLinked const&,
               MySkimIdx const& candidates)
  {

    double massPionFromPDG = MassPiPlus;    // pdg code 211
    double massLambdaFromPDG = MassLambda0; // pdg code 3122
    double massXiFromPDG = MassXiMinus;     // pdg code 3312
    double massOmegacFromPDG = MassOmegaC0; // pdg code 4332
    double massXicFromPDG = MassXiC0;       // pdg code 4132

    // 2-prong vertex fitter to build the omegac/xic vertex
    o2::vertexing::DCAFitterN<2> df;
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMaxDXYIni(maxDXYIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setMaxChi2(maxChi2);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    for (const auto& cand : candidates) {

      hCandidateCounter->Fill(0);

      if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
        continue;
      }

      hCandidateCounter->Fill(1);

      auto collision = cand.collision_as<aod::Collisions>();

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component

      df.setBz(magneticField);
      df.setRefitWithMatCorr(refitWithMatCorr);

      auto trackPion = cand.prong0_as<TracksWCovDca>();
      auto cascAodElement = cand.cascade_as<aod::CascadesLinked>();
      hCascadesCounter->Fill(0);
      int v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_cascData())
        continue;
      auto casc = cascAodElement.cascData_as<MyCascTable>();
      hCascadesCounter->Fill(1);
      auto trackXiDauCharged = casc.bachelor_as<TracksWCovDca>(); // pion <- xi track
      auto trackV0Dau0 = casc.posTrack_as<TracksWCovDca>();       // V0 positive daughter track
      auto trackV0Dau1 = casc.negTrack_as<TracksWCovDca>();       // V0 negative daughter track

      //-------------------------- V0 info---------------------------
      // pseudorapidity
      double pseudorapV0Dau0 = trackV0Dau0.eta();
      double pseudorapV0Dau1 = trackV0Dau1.eta();

      // info from LF table
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      double pseudorapPiFromCas = trackXiDauCharged.eta();

      // info from LF table
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
        covCasc[i] = casc.positionCovMat()[i];
      }
      // create cascade track
      o2::track::TrackParCov trackCasc;
      if (trackXiDauCharged.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackXiDauCharged.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }
      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      //------------reconstruct charm baryon decay vtx---------------
      auto trackParVarPi = getTrackParCov(trackPion); // charm bachelor pion track to be processed with DCAFitter

      // reconstruct charm baryon with DCAFitter
      int nVtxFromFitterCharmBaryon = 0;
      try {
        nVtxFromFitterCharmBaryon = df.process(trackCasc, trackParVarPi);
      } catch (...) {
        LOG(error) << "Exception caught in charm DCA fitter process call!";
        hFitterStatus->Fill(1);
        continue;
      }
      if (nVtxFromFitterCharmBaryon == 0) {
        hFitterStatus->Fill(2);
        continue;
      }
      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);
      auto vertexCharmBaryonFromFitter = df.getPCACandidate();
      std::array<float, 3> pVecCascAsD;
      std::array<float, 3> pVecPionFromCharmBaryon;
      df.propagateTracksToVertex();
      if (!df.isPropagateTracksToVertexDone()) {
        continue;
      }
      df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
      df.getTrack(1).getPxPyPzGlo(pVecPionFromCharmBaryon);
      std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecPionFromCharmBaryon[0], pVecCascAsD[1] + pVecPionFromCharmBaryon[1], pVecCascAsD[2] + pVecPionFromCharmBaryon[2]};

      std::array<float, 3> coordVtxCharmBaryon = df.getPCACandidatePos();
      std::array<float, 6> covVtxCharmBaryon = df.calcPCACovMatrixFlat();

      // pseudorapidity
      double pseudorapPiFromCharmBaryon = trackPion.eta();

      // DCAxy (computed with propagateToDCABxByBz method)
      float dcaxyV0Dau0 = trackV0Dau0.dcaXY();
      float dcaxyV0Dau1 = trackV0Dau1.dcaXY();
      float dcaxyPiFromCasc = trackXiDauCharged.dcaXY();

      // DCAz (computed with propagateToDCABxByBz method)
      float dcazV0Dau0 = trackV0Dau0.dcaZ();
      float dcazV0Dau1 = trackV0Dau1.dcaZ();
      float dcazPiFromCasc = trackXiDauCharged.dcaZ();

      // primary vertex of the collision
      auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterPiFromCharmBaryon;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarPi, 2.f, matCorr, &impactParameterPiFromCharmBaryon);
      float impactParPiFromCharmBaryonXY = impactParameterPiFromCharmBaryon.getY();
      float impactParPiFromCharmBaryonZ = impactParameterPiFromCharmBaryon.getZ();

      // invariant mass under the hypothesis of particles ID corresponding to the decay chain
      double mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
      double mCasc = casc.mXi();
      const std::array<double, 2> arrMassCharmBaryon = {massXiFromPDG, massPionFromPDG};
      double mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecPionFromCharmBaryon}, arrMassCharmBaryon);

      // computing cosPA
      double cpaV0 = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      double cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      double cpaCasc = RecoDecay::cpa(coordVtxCharmBaryon, vertexCasc, pVecCasc);
      double cpaxyV0 = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
      double cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      double cpaxyCasc = RecoDecay::cpaXY(coordVtxCharmBaryon, vertexCasc, pVecCasc);

      // computing decay length and ctau
      double decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      double decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      double decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
      auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
      auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

      double ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, massOmegacFromPDG);
      double ctXic = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, massXicFromPDG);
      double ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, massXiFromPDG);
      double ctV0 = RecoDecay::ct(pVecV0, decLenV0, massLambdaFromPDG);

      // computing eta
      double pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
      double pseudorapCascade = RecoDecay::eta(pVecCasc);
      double pseudorapV0 = RecoDecay::eta(pVecV0);

      // DCA between daughters
      float dcaCascDau = casc.dcacascdaughters();
      float dcaV0Dau = casc.dcaV0daughters();
      float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

      // set hfFlag
      int hfFlag = 1 << aod::hf_cand_toxipi::DecayType::DecayToXiPi;

      // fill test histograms
      hInvMassCharmBaryon->Fill(mCharmBaryon);
      hCandidateCounter->Fill(3);

      // fill the table
      rowCandidate(collision.globalIndex(),
                   pvCoord[0], pvCoord[1], pvCoord[2],
                   vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                   vertexCasc[0], vertexCasc[1], vertexCasc[2],
                   vertexV0[0], vertexV0[1], vertexV0[2],
                   trackXiDauCharged.sign(),
                   covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                   pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                   pVecCasc[0], pVecCasc[1], pVecCasc[2],
                   pVecPionFromCharmBaryon[0], pVecPionFromCharmBaryon[1], pVecPionFromCharmBaryon[2],
                   pVecV0[0], pVecV0[1], pVecV0[2],
                   pVecPionFromCasc[0], pVecPionFromCasc[1], pVecPionFromCasc[2],
                   pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                   pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                   impactParameterCasc.getY(), impactParPiFromCharmBaryonXY,
                   impactParameterCasc.getZ(), impactParPiFromCharmBaryonZ,
                   std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPiFromCharmBaryon.getSigmaY2()),
                   v0index, casc.posTrackId(), casc.negTrackId(),
                   casc.cascadeId(), trackPion.globalIndex(), casc.bachelorId(),
                   mLambda, mCasc, mCharmBaryon,
                   cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                   ctOmegac, ctCascade, ctV0, ctXic,
                   pseudorapV0Dau0, pseudorapV0Dau1, pseudorapPiFromCas, pseudorapPiFromCharmBaryon,
                   pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                   dcaxyV0Dau0, dcaxyV0Dau1, dcaxyPiFromCasc,
                   dcazV0Dau0, dcazV0Dau1, dcazPiFromCasc,
                   dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                   decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon,
                   hfFlag);

    } // loop over LF Cascade-bachelor candidates
  }   // end of process
};    // end of struct

/// Performs MC matching.
struct HfCandidateCreatorToXiPiMc {
  Produces<aod::HfToXiPiMCRec> rowMCMatchRec;
  Produces<aod::HfToXiPiMCGen> rowMCMatchGen;

  Configurable<bool> matchOmegacMc{"matchOmegacMc", true, "Do MC matching for Omegac0"};
  Configurable<bool> matchXicMc{"matchXicMc", false, "Do MC matching for Xic0"};

  void init(InitContext const&) {}

  void processDoNoMc(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPiMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandToXiPi const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    float ptCharmBaryonGen = -999.;
    int indexRec = -1;
    int indexRecCharmBaryon = -1;
    int8_t sign = -9;
    int8_t flag = 0;
    int8_t origin = 0; // to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenXi = 0;
    int8_t debugGenLambda = 0;

    int pdgCodeOmegac0 = Pdg::kOmegaC0; // 4332
    int pdgCodeXic0 = Pdg::kXiC0;       // 4132
    int pdgCodeXiMinus = kXiMinus;      // 3312
    int pdgCodeLambda = kLambda0;       // 3122
    int pdgCodePiPlus = kPiPlus;        // 211
    int pdgCodePiMinus = kPiMinus;      // -211
    int pdgCodeProton = kProton;        // 2212

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = 0;
      auto arrayDaughters = std::array{candidate.piFromCharmBaryon_as<aod::TracksWMc>(), // pi <- charm baryon
                                       candidate.bachelor_as<aod::TracksWMc>(),          // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),          // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()};         // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Omegac matching
      if (matchOmegacMc) {
        // Omegac → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi);
            }
          }
        }

        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }
      // Xic matching
      if (matchXicMc) {
        // Xic → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeXic0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi);
            }
          }
        }

        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }

      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug, origin);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      ptCharmBaryonGen = -999.;
      flag = 0;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenXi = 0;
      debugGenLambda = 0;
      origin = RecoDecay::OriginType::None;
      if (matchOmegacMc) {
        //  Omegac → Xi pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          // Xi -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            debugGenXi = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi);
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }
      if (matchXicMc) {
        //  Xic → Xi pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeXic0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          // Xi- -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            debugGenXi = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi);
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }

      rowMCMatchGen(flag, debugGenCharmBar, debugGenXi, debugGenLambda, ptCharmBaryonGen, origin);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorToXiPiMc, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorToXiPi>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorToXiPiMc>(cfgc)};
}
