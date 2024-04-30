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
#include "Framework/RunningWorkflowInfo.h"
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
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::track;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::aod::hf_track_index;
using namespace o2::aod::hf_collision_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_evsel;

// Reconstruction of omegac0 and xic0 candidates
struct HfCandidateCreatorToXiPi {
  Produces<aod::HfCandToXiPi> rowCandidate;

  // centrality
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality"};
  // event selection
  Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "apply the sel8 event selection"};
  Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "max. PV posZ (cm)"};
  Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "apply TF border cut"};

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

  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter to build the omegac/xic vertex
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{-1};
  double magneticField{0.};

  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>; // to use strangeness tracking, use aod::TraCascDatas instead of aod::CascDatas
  using CascadesLinked = soa::Join<Cascades, CascDataLink>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MySkimIdx = soa::Filtered<aod::HfCascLf2Prongs>;

  Filter filterSelectIndexes = (aod::hf_track_index::hfflag > static_cast<uint8_t>(0));

  std::shared_ptr<TH1> hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 3> processes = {doprocessNoCent, doprocessCentFT0C, doprocessCentFT0M};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    std::array<bool, 3> processesCollisions = {doprocessCollisions, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};
    const int nProcessesCollisions = std::accumulate(processesCollisions.begin(), processesCollisions.end(), 0);
    if (nProcessesCollisions > 1) {
      LOGP(fatal, "At most one process function for collision monitoring can be enabled at a time.");
    }
    if (nProcessesCollisions == 1) {
      if (doprocessNoCent && !doprocessCollisions) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisions\"?");
      }
      if (doprocessCentFT0C && !doprocessCollisionsCentFT0C) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0C\"?");
      }
      if (doprocessCentFT0M && !doprocessCollisionsCentFT0M) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0M\"?");
      }
    }

    registry.add("hInvMassCharmBaryon", "Charm baryon invariant mass;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 2.2, 3.1}}});
    registry.add("hFitterStatus", "Charm DCAFitter status;status;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}});                 // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    registry.add("hCandidateCounter", "Candidate counter wrt derived data;status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}); // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    registry.add("hCascadesCounter", "Cascades counter wrt derived data;status;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});   // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table
    hCollisions = registry.add<TH1>("hCollisions", "HF event counter;;entries", {HistType::kTH1D, {axisEvents}});
    hPosZBeforeEvSel = registry.add<TH1>("hPosZBeforeEvSel", "all events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{400, -20., 20.}}});
    hPosZAfterEvSel = registry.add<TH1>("hPosZAfterEvSel", "selected events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{400, -20., 20.}}});
    hPosXAfterEvSel = registry.add<TH1>("hPosXAfterEvSel", "selected events;#it{x}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{200, -0.5, 0.5}}});
    hPosYAfterEvSel = registry.add<TH1>("hPosYAfterEvSel", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{200, -0.5, 0.5}}});
    hNumPvContributorsAfterSel = registry.add<TH1>("hNumPvContributorsAfterSel", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{500, -0.5, 499.5}}});

    /// collision monitoring
    setLabelHistoEvSel(hCollisions);

    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMaxDXYIni(maxDXYIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setMaxChi2(maxChi2);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);
    df.setRefitWithMatCorr(refitWithMatCorr);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  template <o2::aod::hf_collision_centrality::CentralityEstimator centEstimator, typename Coll>
  void runCreatorToXiPi(Coll const&,
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

    for (const auto& cand : candidates) {

      registry.fill(HIST("hCandidateCounter"), 0);

      if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
        continue;
      }

      registry.fill(HIST("hCandidateCounter"), 1);

      auto collision = cand.collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = getHfCollisionRejectionMask<true, centEstimator>(collision, centrality, centralityMin, centralityMax, useSel8Trigger, -1, useTimeFrameBorderCut, -zPvPosMax, zPvPosMax, 0, -1.f);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      // set the magnetic field from CCDB
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << magneticField;
        runNumber = bc.runNumber();
      }
      df.setBz(magneticField);

      auto trackPion = cand.prong0_as<TracksWCovDca>();
      auto cascAodElement = cand.cascade_as<aod::CascadesLinked>();
      registry.fill(HIST("hCascadesCounter"), 0);
      int v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_cascData()) {
        continue;
      }
      auto casc = cascAodElement.cascData_as<MyCascTable>();
      registry.fill(HIST("hCascadesCounter"), 1);
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
        registry.fill(HIST("hFitterStatus"), 1);
        continue;
      }
      if (nVtxFromFitterCharmBaryon == 0) {
        registry.fill(HIST("hFitterStatus"), 2);
        continue;
      }
      registry.fill(HIST("hFitterStatus"), 0);
      registry.fill(HIST("hCandidateCounter"), 2);
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
      double cpaV0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      double cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      double cpaCasc = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      double cpaxyV0 = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      double cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      double cpaxyCasc = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);

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
      registry.fill(HIST("hInvMassCharmBaryon"), mCharmBaryon);
      registry.fill(HIST("hCandidateCounter"), 3);

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
  }   // end of run function

  /// @brief process function w/o centrality selections
  void processNoCent(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                     aod::BCsWithTimestamps const& bcWithTimeStamps,
                     TracksWCovDca const& tracks,
                     MyCascTable const& cascades,
                     CascadesLinked const& cascadeLinks,
                     MySkimIdx const& candidates)
  {
    runCreatorToXiPi<CentralityEstimator::None>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates);
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPi, processNoCent, "Run candidate creator w/o centrality selections", true);

  /// @brief process function w/o centrality selections
  void processCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                       aod::BCsWithTimestamps const& bcWithTimeStamps,
                       TracksWCovDca const& tracks,
                       MyCascTable const& cascades,
                       CascadesLinked const& cascadeLinks,
                       MySkimIdx const& candidates)
  {
    runCreatorToXiPi<CentralityEstimator::FT0C>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates);
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPi, processCentFT0C, "Run candidate creator w/ centrality selection on FT0C", false);

  /// @brief process function w/o centrality selections
  void processCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                       aod::BCsWithTimestamps const& bcWithTimeStamps,
                       TracksWCovDca const& tracks,
                       MyCascTable const& cascades,
                       CascadesLinked const& cascadeLinks,
                       MySkimIdx const& candidates)
  {
    runCreatorToXiPi<CentralityEstimator::FT0M>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates);
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPi, processCentFT0M, "Run candidate creator w/ centrality selection on FT0M", false);

  ///////////////////////////////////////////////////////////
  ///                                                     ///
  ///   Process functions only for collision monitoring   ///
  ///                                                     ///
  ///////////////////////////////////////////////////////////

  /// @brief process function to monitor collisions - no centrality
  void processCollisions(soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = getHfCollisionRejectionMask<true, CentralityEstimator::None>(collision, centrality, centralityMin, centralityMax, useSel8Trigger, -1, useTimeFrameBorderCut, -zPvPosMax, zPvPosMax, 0, -1.f);

      /// monitor the satisfied event selections
      monitorCollision(collision, rejectionMask, hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPi, processCollisions, "Collision monitoring - no centrality", true);

  /// @brief process function to monitor collisions - FT0C centrality
  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = getHfCollisionRejectionMask<true, CentralityEstimator::FT0C>(collision, centrality, centralityMin, centralityMax, useSel8Trigger, -1, useTimeFrameBorderCut, -zPvPosMax, zPvPosMax, 0, -1.f);

      /// monitor the satisfied event selections
      monitorCollision(collision, rejectionMask, hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPi, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

  /// @brief process function to monitor collisions - FT0M centrality
  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = getHfCollisionRejectionMask<true, CentralityEstimator::FT0M>(collision, centrality, centralityMin, centralityMax, useSel8Trigger, -1, useTimeFrameBorderCut, -zPvPosMax, zPvPosMax, 0, -1.f);

      /// monitor the satisfied event selections
      monitorCollision(collision, rejectionMask, hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPi, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);

};    // end of struct

/// Performs MC matching.
struct HfCandidateCreatorToXiPiMc {
  Produces<aod::HfToXiPiMCRec> rowMCMatchRec;
  Produces<aod::HfToXiPiMCGen> rowMCMatchGen;

  Configurable<bool> matchOmegacMc{"matchOmegacMc", true, "Do MC matching for Omegac0"};
  Configurable<bool> matchXicMc{"matchXicMc", false, "Do MC matching for Xic0"};

  float zPvPosMax{1000.f};

  // inspect for which zPvPosMax cut was set for reconstructed
  void init(InitContext& initContext)
  {
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-to-xi-pi") == 0) {
        for (const auto& option : device.options) {
          if (option.name.compare("zPvPosMax") == 0) {
            zPvPosMax = option.defaultValue.get<float>();
            break;
          }
        }
        break;
      }
    }
  }

  void processDoNoMc(aod::Collisions::iterator const&)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPiMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandToXiPi const& candidates,
                 aod::TracksWMc const&,
                 aod::McParticles const& mcParticles,
                 aod::McCollisionLabels const&)
  {
    float ptCharmBaryonGen = -999.;
    float etaCharmBaryonGen = -999.;
    int indexRec = -1;
    int indexRecCharmBaryon = -1;
    int8_t sign = -9;
    int8_t signCasc = -9;
    int8_t signV0 = -9;
    int8_t flag = 0;
    int8_t origin = 0; // to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenXi = 0;
    int8_t debugGenLambda = 0;
    bool collisionMatched = false;

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
      collisionMatched = false;

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
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi);
              collisionMatched = candidate.collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
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
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeXic0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi);
              collisionMatched = candidate.collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
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
      rowMCMatchRec(flag, debug, origin, collisionMatched);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      ptCharmBaryonGen = -999.;
      etaCharmBaryonGen = -999.;
      flag = 0;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenXi = 0;
      debugGenLambda = 0;
      origin = RecoDecay::OriginType::None;

      auto mcCollision = particle.mcCollision();
      float zPv = mcCollision.posZ();
      if (zPv < -zPvPosMax || zPv > zPvPosMax) { // to avoid counting particles in collisions with Zvtx larger than the maximum, we do not match them
        rowMCMatchGen(flag, debugGenCharmBar, debugGenXi, debugGenLambda, ptCharmBaryonGen, etaCharmBaryonGen, origin);
        continue;
      }

      if (matchOmegacMc) {
        //  Omegac → Xi pi
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          etaCharmBaryonGen = particle.eta();
          for (const auto& daughterCharm : particle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgCodeXiMinus) {
              continue;
            }
            // Xi -> Lambda pi
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
              debugGenXi = 1;
              for (const auto& daughterCascade : daughterCharm.daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                  continue;
                }
                // Lambda -> p pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                  debugGenLambda = 1;
                  flag = sign * (1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi);
                }
              }
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
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeXic0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          etaCharmBaryonGen = particle.eta();
          for (const auto& daughterCharm : particle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgCodeXiMinus) {
              continue;
            }
            // Xi -> Lambda pi
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
              debugGenXi = 1;
              for (const auto& daughterCascade : daughterCharm.daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                  continue;
                }
                // Lambda -> p pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                  debugGenLambda = 1;
                  flag = sign * (1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi);
                }
              }
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }

      rowMCMatchGen(flag, debugGenCharmBar, debugGenXi, debugGenLambda, ptCharmBaryonGen, etaCharmBaryonGen, origin);
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
