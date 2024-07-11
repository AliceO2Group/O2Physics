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

/// \file candidateCreatorXic0Omegac0.cxx
/// \brief Reconstruction of Omegac0 and Xic0 decays candidates
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

#include "PWGHF/Core/CentralityEstimation.h"
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
using namespace o2::hf_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_evsel;

// Reconstruction of omegac0 and xic0 candidates
struct HfCandidateCreatorXic0Omegac0 {
  Produces<aod::HfCandToXiPi> rowCandToXiPi;
  Produces<aod::HfCandToOmegaPi> rowCandToOmegaPi;
  Produces<aod::HfCandToOmegaK> rowCandToOmegaK;

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

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter to build the omegac/xic vertex
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{-1};
  double magneticField{0.};

  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>; // to use strangeness tracking, use aod::TraCascDatas instead of aod::CascDatas
  using CascadesLinked = soa::Join<Cascades, CascDataLink>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;

  std::shared_ptr<TH1> hInvMassCharmBaryonToXiPi, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaK, hFitterStatusToXiPi, hFitterStatusToOmegaPi, hFitterStatusToOmegaK, hCandidateCounterToXiPi, hCandidateCounterToOmegaPi, hCandidateCounterToOmegaK, hCascadesCounterToXiPi, hCascadesCounterToOmegaPi, hCascadesCounterToOmegaK;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 9> allProcesses = {doprocessNoCentToXiPi, doprocessCentFT0CToXiPi, doprocessCentFT0MToXiPi, doprocessNoCentToOmegaPi, doprocessCentFT0CToOmegaPi, doprocessCentFT0MToOmegaPi, doprocessNoCentToOmegaK, doprocessCentFT0CToOmegaK, doprocessCentFT0MToOmegaK};
    if (std::accumulate(allProcesses.begin(), allProcesses.end(), 0) == 0) {
      LOGP(fatal, "No process function enabled, please select one for at least one channel.");
    }

    std::array<bool, 3> processesToXiPi = {doprocessNoCentToXiPi, doprocessCentFT0CToXiPi, doprocessCentFT0MToXiPi};
    if (std::accumulate(processesToXiPi.begin(), processesToXiPi.end(), 0) > 1) {
      LOGP(fatal, "One and only one ToXiPi process function must be enabled at a time.");
    }
    std::array<bool, 3> processesToOmegaPi = {doprocessNoCentToOmegaPi, doprocessCentFT0CToOmegaPi, doprocessCentFT0MToOmegaPi};
    if (std::accumulate(processesToOmegaPi.begin(), processesToOmegaPi.end(), 0) > 1) {
      LOGP(fatal, "One and only one process ToOmegaPi function must be enabled at a time.");
    }
    std::array<bool, 3> processesToOmegaK = {doprocessNoCentToOmegaK, doprocessCentFT0CToOmegaK, doprocessCentFT0MToOmegaK};
    if (std::accumulate(processesToOmegaK.begin(), processesToOmegaK.end(), 0) > 1) {
      LOGP(fatal, "One and only one process ToOmegaK function must be enabled at a time.");
    }

    std::array<bool, 3> processesCollisions = {doprocessCollisions, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};
    const int nProcessesCollisions = std::accumulate(processesCollisions.begin(), processesCollisions.end(), 0);
    if (nProcessesCollisions > 1) {
      LOGP(fatal, "At most one process function for collision monitoring can be enabled at a time.");
    }
    if (nProcessesCollisions == 1) {
      if ((doprocessNoCentToXiPi && !doprocessCollisions) || (doprocessNoCentToOmegaPi && !doprocessCollisions) || (doprocessNoCentToOmegaK && !doprocessCollisions)) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisions\"?");
      }
      if ((doprocessCentFT0CToXiPi && !doprocessCollisionsCentFT0C) || (doprocessCentFT0CToOmegaPi && !doprocessCollisionsCentFT0C) || (doprocessCentFT0CToOmegaK && !doprocessCollisionsCentFT0C)) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0C\"?");
      }
      if ((doprocessCentFT0MToXiPi && !doprocessCollisionsCentFT0M) || (doprocessCentFT0MToOmegaPi && !doprocessCollisionsCentFT0M) || (doprocessCentFT0MToOmegaK && !doprocessCollisionsCentFT0M)) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0M\"?");
      }
    }

    hInvMassCharmBaryonToXiPi = registry.add<TH1>("hInvMassCharmBaryonToXiPi", "Charm baryon invariant mass - #Xi #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 2.2, 3.1}}});
    hInvMassCharmBaryonToOmegaPi = registry.add<TH1>("hInvMassCharmBaryonToOmegaPi", "Charm baryon invariant mass - #Omega #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 2.2, 3.1}}});
    hInvMassCharmBaryonToOmegaK = registry.add<TH1>("hInvMassCharmBaryonToOmegaK", "Charm baryon invariant mass - #Omega K decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 2.2, 3.1}}});
    hFitterStatusToXiPi = registry.add<TH1>("hFitterStatusToXiPi", "Charm DCAFitter status - #Xi #pi vtx;status;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}});                                // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    hFitterStatusToOmegaPi = registry.add<TH1>("hFitterStatusToOmegaPi", "Charm DCAFitter status - #Omega #pi vtx ;status;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}});                      // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    hFitterStatusToOmegaK = registry.add<TH1>("hFitterStatusToOmegaK", "Charm DCAFitter status - #Omega K vtx;status;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}});                           // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    hCandidateCounterToXiPi = registry.add<TH1>("hCandidateCounterToXiPi", "Candidate counter wrt derived data - #Xi #pi decay;status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}});          // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    hCandidateCounterToOmegaPi = registry.add<TH1>("hCandidateCounterToOmegaPi", "Candidate counter wrt derived data - #Omega #pi decay;status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}); // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    hCandidateCounterToOmegaK = registry.add<TH1>("hCandidateCounterToOmegaK", "Candidate counter wrt derived data - #Omega K decay;status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}});     // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    hCascadesCounterToXiPi = registry.add<TH1>("hCascadesCounterToXiPi", "Cascades counter wrt derived data - #Xi #pi decay;status;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});             // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table
    hCascadesCounterToOmegaPi = registry.add<TH1>("hCascadesCounterToOmegaPi", "Cascades counter wrt derived data - #Omega #pi decay;status;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});    // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table
    hCascadesCounterToOmegaK = registry.add<TH1>("hCascadesCounterToOmegaK", "Cascades counter wrt derived data - #Omega K decay;status;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});        // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table

    hfEvSel.addHistograms(registry); // collision monitoring

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

  template <o2::hf_centrality::CentralityEstimator centEstimator, int decayChannel, typename Coll, typename Hist>
  void runXic0Omegac0Creator(Coll const&,
                             aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                             TracksWCovDca const&,
                             MyCascTable const&, CascadesLinked const&,
                             aod::HfCascLf2Prongs const& candidates,
                             Hist& hInvMassCharmBaryon,
                             Hist& hFitterStatus,
                             Hist& hCandidateCounter,
                             Hist& hCascadesCounter)
  {

    if constexpr (decayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi && decayChannel != hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi && decayChannel != hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
      LOGP(fatal, "Decay channel not recognized!");
    }

    for (const auto& cand : candidates) {

      hCandidateCounter->Fill(0);

      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
          continue;
        }
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
          continue;
        }
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK)) {
          continue;
        }
      }

      hCandidateCounter->Fill(1);

      auto collision = cand.collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb);
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

      auto trackCharmBachelor = cand.prong0_as<TracksWCovDca>();

      auto cascAodElement = cand.cascade_as<aod::CascadesLinked>();
      hCascadesCounter->Fill(0);
      int v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_cascData()) {
        continue;
      }
      auto casc = cascAodElement.cascData_as<MyCascTable>();
      hCascadesCounter->Fill(1);
      auto trackCascDauCharged = casc.bachelor_as<TracksWCovDca>(); // pion <- xi track
      auto trackV0Dau0 = casc.posTrack_as<TracksWCovDca>();         // V0 positive daughter track
      auto trackV0Dau1 = casc.negTrack_as<TracksWCovDca>();         // V0 negative daughter track

      //-------------------------- V0 info---------------------------
      // pseudorapidity
      float pseudorapV0Dau0 = trackV0Dau0.eta();
      float pseudorapV0Dau1 = trackV0Dau1.eta();

      // info from LF table
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      float pseudorapCascBachelor = trackCascDauCharged.eta();

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
      if (trackCascDauCharged.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackCascDauCharged.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }
      trackCasc.setAbsCharge(1);
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        trackCasc.setPID(o2::track::PID::XiMinus);
      } else {
        trackCasc.setPID(o2::track::PID::OmegaMinus);
      }

      std::array<float, 3> pVecCascBachelor = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      //------------reconstruct charm baryon decay vtx---------------
      auto trackParVarCharmBachelor = getTrackParCov(trackCharmBachelor); // charm bachelor pion track to be processed with DCAFitter

      // reconstruct charm baryon with DCAFitter
      int nVtxFromFitterCharmBaryon = 0;
      try {
        nVtxFromFitterCharmBaryon = df.process(trackCasc, trackParVarCharmBachelor);
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
      std::array<float, 3> pVecCharmBachelorAsD;
      df.propagateTracksToVertex();
      if (!df.isPropagateTracksToVertexDone()) {
        continue;
      }
      df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
      df.getTrack(1).getPxPyPzGlo(pVecCharmBachelorAsD);
      std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecCharmBachelorAsD[0], pVecCascAsD[1] + pVecCharmBachelorAsD[1], pVecCascAsD[2] + pVecCharmBachelorAsD[2]};

      std::array<float, 3> coordVtxCharmBaryon = df.getPCACandidatePos();
      std::array<float, 6> covVtxCharmBaryon = df.calcPCACovMatrixFlat();

      // pseudorapidity
      float pseudorapCharmBachelor = trackCharmBachelor.eta();

      // DCAxy (computed with propagateToDCABxByBz method)
      float dcaxyV0Dau0 = trackV0Dau0.dcaXY();
      float dcaxyV0Dau1 = trackV0Dau1.dcaXY();
      float dcaxyCascBachelor = trackCascDauCharged.dcaXY();

      // DCAz (computed with propagateToDCABxByBz method)
      float dcazV0Dau0 = trackV0Dau0.dcaZ();
      float dcazV0Dau1 = trackV0Dau1.dcaZ();
      float dcazCascBachelor = trackCascDauCharged.dcaZ();

      // primary vertex of the collision
      auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterCharmBachelor;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);
      float impactParBachFromCharmBaryonXY = impactParameterCharmBachelor.getY();
      float impactParBachFromCharmBaryonZ = impactParameterCharmBachelor.getZ();

      // invariant mass under the hypothesis of particles ID corresponding to the decay chain
      float mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
      float mCasc = 0.;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        mCasc = casc.mXi();
      } else {
        mCasc = casc.mOmega();
      }
      auto arrMassCharmBaryon = std::array{0., 0.};
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        arrMassCharmBaryon = {MassXiMinus, MassPiPlus};
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        arrMassCharmBaryon = {MassOmegaMinus, MassPiPlus};
      } else {
        arrMassCharmBaryon = {MassOmegaMinus, MassKPlus};
      }
      float mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecCharmBachelorAsD}, arrMassCharmBaryon);

      // computing cosPA
      float cpaV0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      float cpaCasc = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaxyV0 = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      float cpaxyCasc = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);

      // computing decay length and ctau
      float decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      float decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      float decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
      auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
      auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

      float ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassOmegaC0);
      float ctXic = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassXiC0);
      float ctCascade = 0.;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassXiMinus);
      } else {
        ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassOmegaMinus);
      }
      float ctV0 = RecoDecay::ct(pVecV0, decLenV0, MassLambda0);

      // computing eta
      float pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
      float pseudorapCascade = RecoDecay::eta(pVecCasc);
      float pseudorapV0 = RecoDecay::eta(pVecV0);

      // DCA between daughters
      float dcaCascDau = casc.dcacascdaughters();
      float dcaV0Dau = casc.dcaV0daughters();
      float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

      // fill test histograms
      hInvMassCharmBaryon->Fill(mCharmBaryon);
      hCandidateCounter->Fill(3);

      // fill the table
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        rowCandToXiPi(collision.globalIndex(),
                      pvCoord[0], pvCoord[1], pvCoord[2],
                      vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                      vertexCasc[0], vertexCasc[1], vertexCasc[2],
                      vertexV0[0], vertexV0[1], vertexV0[2],
                      trackCascDauCharged.sign(),
                      covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                      pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                      pVecCasc[0], pVecCasc[1], pVecCasc[2],
                      pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                      pVecV0[0], pVecV0[1], pVecV0[2],
                      pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
                      pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                      pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                      impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
                      impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
                      std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                      v0index, casc.posTrackId(), casc.negTrackId(),
                      casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                      mLambda, mCasc, mCharmBaryon,
                      cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                      ctOmegac, ctCascade, ctV0, ctXic,
                      pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                      pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                      dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                      dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
                      dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                      decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);

      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        rowCandToOmegaPi(collision.globalIndex(),
                         pvCoord[0], pvCoord[1], pvCoord[2],
                         vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                         vertexCasc[0], vertexCasc[1], vertexCasc[2],
                         vertexV0[0], vertexV0[1], vertexV0[2],
                         trackCascDauCharged.sign(),
                         covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                         pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                         pVecCasc[0], pVecCasc[1], pVecCasc[2],
                         pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                         pVecV0[0], pVecV0[1], pVecV0[2],
                         pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
                         pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                         pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                         impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
                         impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
                         std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                         v0index, casc.posTrackId(), casc.negTrackId(),
                         casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                         mLambda, mCasc, mCharmBaryon,
                         cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                         ctOmegac, ctCascade, ctV0,
                         pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                         pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                         dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                         dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
                         dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                         decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);

      } else {
        rowCandToOmegaK(
          collision.globalIndex(), pvCoord[0], pvCoord[1], pvCoord[2],
          vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
          vertexCasc[0], vertexCasc[1], vertexCasc[2],
          vertexV0[0], vertexV0[1], vertexV0[2],
          trackCascDauCharged.sign(),
          covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
          pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
          pVecCasc[0], pVecCasc[1], pVecCasc[2],
          pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
          pVecV0[0], pVecV0[1], pVecV0[2],
          pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
          pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
          pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
          impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
          impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
          std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
          v0index, casc.posTrackId(), casc.negTrackId(),
          casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
          mLambda, mCasc, mCharmBaryon,
          cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
          ctOmegac, ctCascade, ctV0,
          pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
          pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
          dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
          dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
          dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
          decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);
      }

    } // loop over LF Cascade-bachelor candidates
  }   // end of run function

  /// @brief process function w/o centrality selections
  void processNoCentToXiPi(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                           aod::BCsWithTimestamps const& bcWithTimeStamps,
                           TracksWCovDca const& tracks,
                           MyCascTable const& cascades,
                           CascadesLinked const& cascadeLinks,
                           aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToXiPi, "Run candidate creator w/o centrality selections for xi pi decay channel", true);

  void processNoCentToOmegaPi(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              aod::BCsWithTimestamps const& bcWithTimeStamps,
                              TracksWCovDca const& tracks,
                              MyCascTable const& cascades,
                              CascadesLinked const& cascadeLinks,
                              aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToOmegaPi, "Run candidate creator w/o centrality selections for omega pi decay channel", false);

  void processNoCentToOmegaK(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                             aod::BCsWithTimestamps const& bcWithTimeStamps,
                             TracksWCovDca const& tracks,
                             MyCascTable const& cascades,
                             CascadesLinked const& cascadeLinks,
                             aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToOmegaK, "Run candidate creator w/o centrality selections for omega K decay channel", false);

  /// @brief process function w/ FT0C centrality selections
  void processCentFT0CToXiPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                             aod::BCsWithTimestamps const& bcWithTimeStamps,
                             TracksWCovDca const& tracks,
                             MyCascTable const& cascades,
                             CascadesLinked const& cascadeLinks,
                             aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CToXiPi, "Run candidate creator w/ centrality selection on FT0C for xi pi channel", false);

  void processCentFT0CToOmegaPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                TracksWCovDca const& tracks,
                                MyCascTable const& cascades,
                                CascadesLinked const& cascadeLinks,
                                aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CToOmegaPi, "Run candidate creator w/ centrality selection on FT0C for omega pi channel", false);

  void processCentFT0CToOmegaK(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                               aod::BCsWithTimestamps const& bcWithTimeStamps,
                               TracksWCovDca const& tracks,
                               MyCascTable const& cascades,
                               CascadesLinked const& cascadeLinks,
                               aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CToOmegaK, "Run candidate creator w/ centrality selection on FT0C for omega K channel", false);

  /// @brief process function w/ FT0M centrality selections
  void processCentFT0MToXiPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                             aod::BCsWithTimestamps const& bcWithTimeStamps,
                             TracksWCovDca const& tracks,
                             MyCascTable const& cascades,
                             CascadesLinked const& cascadeLinks,
                             aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MToXiPi, "Run candidate creator w/ centrality selection on FT0M for xi pi channel", false);

  void processCentFT0MToOmegaPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                TracksWCovDca const& tracks,
                                MyCascTable const& cascades,
                                CascadesLinked const& cascadeLinks,
                                aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MToOmegaPi, "Run candidate creator w/ centrality selection on FT0M for omega pi channel", false);

  void processCentFT0MToOmegaK(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                               aod::BCsWithTimestamps const& bcWithTimeStamps,
                               TracksWCovDca const& tracks,
                               MyCascTable const& cascades,
                               CascadesLinked const& cascadeLinks,
                               aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MToOmegaK, "Run candidate creator w/ centrality selection on FT0M for omega K channel", false);

  ///////////////////////////////////////////////////////////
  ///                                                     ///
  ///   Process functions only for collision monitoring   ///
  ///                                                     ///
  ///////////////////////////////////////////////////////////

  /// @brief process function to monitor collisions - no centrality
  void processCollisions(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCollisions, "Collision monitoring - no centrality", true);

  /// @brief process function to monitor collisions - FT0C centrality
  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

  /// @brief process function to monitor collisions - FT0M centrality
  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);

}; // end of struct

/// Performs MC matching.
struct HfCandidateCreatorXic0Omegac0Mc {
  Produces<aod::HfXicToXiPiMCRec> rowMCMatchRecXicToXiPi;
  Produces<aod::HfXicToXiPiMCGen> rowMCMatchGenXicToXiPi;
  Produces<aod::HfOmegacToXiPiMCRec> rowMCMatchRecOmegacToXiPi;
  Produces<aod::HfOmegacToXiPiMCGen> rowMCMatchGenOmegacToXiPi;
  Produces<aod::HfToOmegaPiMCRec> rowMCMatchRecToOmegaPi;
  Produces<aod::HfToOmegaPiMCGen> rowMCMatchGenToOmegaPi;
  Produces<aod::HfToOmegaKMCRec> rowMCMatchRecToOmegaK;
  Produces<aod::HfToOmegaKMCGen> rowMCMatchGenToOmegaK;

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  HistogramRegistry registry{"registry"};

  // inspect for which zPvPosMax cut was set for reconstructed
  void init(InitContext& initContext)
  {
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-xic0-omegac0") == 0) {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }
    hfEvSelMc.addHistograms(registry); // particles monitoring
  }

  template <int decayChannel, typename TMyRecoCand>
  void runXic0Omegac0Mc(TMyRecoCand const& candidates,
                        aod::TracksWMc const&,
                        aod::McParticles const& mcParticles,
                        aod::McCollisions const&,
                        aod::McCollisionLabels const&,
                        BCsInfo const&)
  {
    float ptCharmBaryonGen = -999.;
    float rapidityCharmBaryonGen = -999.;
    int indexRec = -1;
    int indexRecCharmBaryon = -1;
    int8_t sign = -9;
    int8_t signCasc = -9;
    int8_t signV0 = -9;
    int8_t flag = 0;
    int8_t origin = 0; // to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenCasc = 0;
    int8_t debugGenLambda = 0;
    bool collisionMatched = false;

    int pdgCodeOmegac0 = Pdg::kOmegaC0;  // 4332
    int pdgCodeXic0 = Pdg::kXiC0;        // 4132
    int pdgCodeXiMinus = kXiMinus;       // 3312
    int pdgCodeOmegaMinus = kOmegaMinus; // 3334
    int pdgCodeLambda = kLambda0;        // 3122
    int pdgCodePiPlus = kPiPlus;         // 211
    int pdgCodePiMinus = kPiMinus;       // -211
    int pdgCodeProton = kProton;         // 2212
    int pdgCodeKaonPlus = kKPlus;        // 321
    int pdgCodeKaonMinus = kKMinus;      // -321

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = 0;
      collisionMatched = false;
      std::vector<int> idxBhadMothers{};

      auto arrayDaughters = std::array{candidate.template bachelorFromCharmBaryon_as<aod::TracksWMc>(), // bachelor <- charm baryon
                                       candidate.template bachelor_as<aod::TracksWMc>(),                // bachelor <- cascade
                                       candidate.template posTrack_as<aod::TracksWMc>(),                // p <- lambda
                                       candidate.template negTrack_as<aod::TracksWMc>()};               // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.template bachelor_as<aod::TracksWMc>(),
                                           candidate.template posTrack_as<aod::TracksWMc>(),
                                           candidate.template negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.template posTrack_as<aod::TracksWMc>(),
                                         candidate.template negTrack_as<aod::TracksWMc>()};
      // Xic0 -> xi pi matching
      if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
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
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi);
              collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecXicToXiPi(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecXicToXiPi(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Xic0ToXiPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) { // Omegac -> xi pi matching
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
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi);
              collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecOmegacToXiPi(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecOmegacToXiPi(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Omegac0ToXiPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) { // Omegac0 -> omega pi matching
        // Omegac → pi K pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Omega- → K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeOmegaMinus, std::array{pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
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
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi);
              collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecToOmegaPi(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecToOmegaPi(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Omegac0ToOmegaPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) { // Omegac0 -> omega K matching
        // Omegac → K K pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodeKaonPlus, pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Omega- → K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeOmegaMinus, std::array{pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
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
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK);
              collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecToOmegaK(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecToOmegaK(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Omegac0ToOmegaK decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      }
    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      ptCharmBaryonGen = -999.;
      rapidityCharmBaryonGen = -999.;
      flag = 0;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenCasc = 0;
      debugGenLambda = 0;
      origin = RecoDecay::OriginType::None;
      std::vector<int> idxBhadMothers{};

      auto mcCollision = particle.mcCollision();
      const auto rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo>(mcCollision);
      hfEvSelMc.fillHistograms(rejectionMask);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the gen particle
        if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
          rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
          rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
          rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
          rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        }
        continue;
      }

      if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
        //  Xic → Xi pi
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeXic0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          rapidityCharmBaryonGen = particle.y();
          for (const auto& daughterCharm : particle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgCodeXiMinus) {
              continue;
            }
            // Xi -> Lambda pi
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
              debugGenCasc = 1;
              for (const auto& daughterCascade : daughterCharm.daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                  continue;
                }
                // Lambda -> p pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                  debugGenLambda = 1;
                  flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi);
                }
              }
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
        } else {
          rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        }

      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
        //  Omegac → Xi pi
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          rapidityCharmBaryonGen = particle.y();
          for (const auto& daughterCharm : particle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgCodeXiMinus) {
              continue;
            }
            // Xi -> Lambda pi
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
              debugGenCasc = 1;
              for (const auto& daughterCascade : daughterCharm.daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                  continue;
                }
                // Lambda -> p pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                  debugGenLambda = 1;
                  flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi);
                }
              }
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
        } else {
          rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        }

      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
        //  Omegac → Omega pi
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeOmegaMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          rapidityCharmBaryonGen = particle.y();
          for (const auto& daughterCharm : particle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgCodeOmegaMinus) {
              continue;
            }
            // Omega -> Lambda K
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeOmegaMinus, std::array{pdgCodeLambda, pdgCodeKaonMinus}, true)) {
              debugGenCasc = 1;
              for (const auto& daughterCascade : daughterCharm.daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                  continue;
                }
                // Lambda -> p pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                  debugGenLambda = 1;
                  flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi);
                }
              }
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
        } else {
          rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        }

      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
        //  Omegac → Omega K
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeOmegaMinus, pdgCodeKaonPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
          rapidityCharmBaryonGen = particle.y();
          for (const auto& daughterCharm : particle.daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgCodeOmegaMinus) {
              continue;
            }
            // Omega -> Lambda K
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeOmegaMinus, std::array{pdgCodeLambda, pdgCodeKaonMinus}, true)) {
              debugGenCasc = 1;
              for (const auto& daughterCascade : daughterCharm.daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                  continue;
                }
                // Lambda -> p pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                  debugGenLambda = 1;
                  flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK);
                }
              }
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
        } else {
          rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
        }
      }
    } // close loop on MCParticles
  }   // close process

  void processDoNoMc(aod::Collisions::iterator const&)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processDoNoMc, "Do not run any MC process function", true);

  void processMcXicToXiPi(aod::HfCandToXiPi const& candidates,
                          aod::TracksWMc const& tracks,
                          aod::McParticles const& mcParticles,
                          aod::McCollisions const& mcColls,
                          aod::McCollisionLabels const& mcLabels,
                          BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, mcColls, mcLabels, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcXicToXiPi, "Run Xic0 to xi pi MC process function", false);

  void processMcOmegacToXiPi(aod::HfCandToXiPi const& candidates,
                             aod::TracksWMc const& tracks,
                             aod::McParticles const& mcParticles,
                             aod::McCollisions const& mcColls,
                             aod::McCollisionLabels const& mcLabels,
                             BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, mcColls, mcLabels, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToXiPi, "Run Omegac0 to xi pi MC process function", false);

  void processMcOmegacToOmegaPi(aod::HfCandToOmegaPi const& candidates,
                                aod::TracksWMc const& tracks,
                                aod::McParticles const& mcParticles,
                                aod::McCollisions const& mcColls,
                                aod::McCollisionLabels const& mcLabels,
                                BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, mcColls, mcLabels, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaPi, "Run Omegac0 to omega pi MC process function", false);

  void processMcOmegacToOmegaK(aod::HfCandToOmegaK const& candidates,
                               aod::TracksWMc const& tracks,
                               aod::McParticles const& mcParticles,
                               aod::McCollisions const& mcColls,
                               aod::McCollisionLabels const& mcLabels,
                               BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, mcColls, mcLabels, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaK, "Run Omegac0 to omega K MC process function", false);

}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0Mc>(cfgc)};
}
