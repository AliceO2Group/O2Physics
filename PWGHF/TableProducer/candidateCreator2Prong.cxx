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

/// \file candidateCreator2Prong.cxx
/// \brief Reconstruction of heavy-flavour 2-prong decay candidates
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Pengzhong Lu <pengzhong.lu@cern.ch>, GSI Darmstadt, USTC

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include <memory>
#include <string>
#include <vector>

#include <KFParticleBase.h>
#include <KFParticle.h>
#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFVertex.h>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGHF/Utils/utilsMcGen.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::aod::pid_tpc_tof_utils;

/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfCandidateCreator2Prong {
  Produces<aod::HfCand2ProngBase> rowCandidateBase;
  Produces<aod::HfProng0PidPi> rowProng0PidPi;
  Produces<aod::HfProng0PidKa> rowProng0PidKa;
  Produces<aod::HfProng1PidPi> rowProng1PidPi;
  Produces<aod::HfProng1PidKa> rowProng1PidKa;
  Produces<aod::HfCand2ProngKF> rowCandidateKF;

  // vertexing
  Configurable<bool> constrainKfToPv{"constrainKfToPv", true, "constraint KFParticle to PV"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using TracksWCovExtraPidPiKa = soa::Join<aod::TracksWCovExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;

  int runNumber{0};
  float toMicrometers = 10000.; // from cm to µm
  double massPi{0.};
  double massK{0.};
  double massPiK{0.};
  double massKPi{0.};
  double bz{0.};

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 6> doprocessDF{doprocessPvRefitWithDCAFitterN, doprocessNoPvRefitWithDCAFitterN,
                                    doprocessPvRefitWithDCAFitterNCentFT0C, doprocessNoPvRefitWithDCAFitterNCentFT0C,
                                    doprocessPvRefitWithDCAFitterNCentFT0M, doprocessNoPvRefitWithDCAFitterNCentFT0M};
    std::array<bool, 6> doprocessKF{doprocessPvRefitWithKFParticle, doprocessNoPvRefitWithKFParticle,
                                    doprocessPvRefitWithKFParticleCentFT0C, doprocessNoPvRefitWithKFParticleCentFT0C,
                                    doprocessPvRefitWithKFParticleCentFT0M, doprocessNoPvRefitWithKFParticleCentFT0M};
    if ((std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) + std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    std::array<bool, 3> processesCollisions = {doprocessCollisions, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};
    const int nProcessesCollisions = std::accumulate(processesCollisions.begin(), processesCollisions.end(), 0);
    if (nProcessesCollisions > 1) {
      LOGP(fatal, "At most one process function for collision monitoring can be enabled at a time.");
    }
    if (nProcessesCollisions == 1) {
      if ((doprocessPvRefitWithDCAFitterN || doprocessNoPvRefitWithDCAFitterN || doprocessPvRefitWithKFParticle || doprocessNoPvRefitWithKFParticle) && !doprocessCollisions) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisions\"?");
      }
      if ((doprocessPvRefitWithDCAFitterNCentFT0C || doprocessNoPvRefitWithDCAFitterNCentFT0C || doprocessPvRefitWithKFParticleCentFT0C || doprocessNoPvRefitWithKFParticleCentFT0C) && !doprocessCollisionsCentFT0C) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0C\"?");
      }
      if ((doprocessPvRefitWithDCAFitterNCentFT0M || doprocessNoPvRefitWithDCAFitterNCentFT0M || doprocessPvRefitWithKFParticleCentFT0M || doprocessNoPvRefitWithKFParticleCentFT0M) && !doprocessCollisionsCentFT0M) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0M\"?");
      }
    }

    // histograms
    registry.add("hMass2", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.6, 2.1}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hCovPVYY", "2-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVYY", "2-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hCovPVXZ", "2-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, -1.e-4, 1.e-4}}});
    registry.add("hCovSVXZ", "2-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, -1.e-4, 0.2}}});
    registry.add("hCovPVZZ", "2-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVZZ", "2-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hDcaXYProngs", "DCAxy of 2-prong candidate daughters;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
    registry.add("hDcaZProngs", "DCAz of 2-prong candidate daughters;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
    registry.add("hVertexerType", "Use KF or DCAFitterN;Vertexer type;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}}); // See o2::aod::hf_cand::VertexerType
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    hfEvSel.addHistograms(registry); // collision monitoring

    massPi = MassPiPlus;
    massK = MassKPlus;

    if (std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) == 1) {
      registry.fill(HIST("hVertexerType"), aod::hf_cand::VertexerType::DCAFitter);
      // Configure DCAFitterN
      // df.setBz(bz);
      df.setPropagateToPCA(propagateToPCA);
      df.setMaxR(maxR);
      df.setMaxDZIni(maxDZIni);
      df.setMinParamChange(minParamChange);
      df.setMinRelChi2Change(minRelChi2Change);
      df.setUseAbsDCA(useAbsDCA);
      df.setWeightedFinalPCA(useWeightedFinalPCA);
    }
    if (std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0) == 1) {
      registry.fill(HIST("hVertexerType"), aod::hf_cand::VertexerType::KfParticle);
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    /// candidate monitoring
    setLabelHistoCands(hCandidates);
  }

  template <bool doPvRefit, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll, typename CandType, typename TTracks>
  void runCreator2ProngWithDCAFitterN(Coll const&,
                                      CandType const& rowsTrackIndexProng2,
                                      TTracks const&,
                                      aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    // loop over pairs of track indices
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {

      /// reject candidates not satisfying the event selections
      auto collision = rowTrackIndexProng2.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      auto track0 = rowTrackIndexProng2.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexProng2.template prong1_as<TTracks>();
      auto trackParVarPos1 = getTrackParCov(track0);
      auto trackParVarNeg1 = getTrackParCov(track1);

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, nullptr, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      df.setBz(bz);

      // reconstruct the 2-prong secondary vertex
      hCandidates->Fill(SVFitting::BeforeFit);
      try {
        if (df.process(trackParVarPos1, trackParVarNeg1) == 0) {
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        hCandidates->Fill(SVFitting::Fail);
        continue;
      }
      hCandidates->Fill(SVFitting::FitOk);

      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      registry.fill(HIST("hCovSVYY"), covMatrixPCA[2]);
      registry.fill(HIST("hCovSVXZ"), covMatrixPCA[3]);
      registry.fill(HIST("hCovSVZZ"), covMatrixPCA[5]);
      auto trackParVar0 = df.getTrack(0);
      auto trackParVar1 = df.getTrack(1);

      // get track momenta
      std::array<float, 3> pvec0;
      std::array<float, 3> pvec1;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexProng2.pvRefitX());
        primaryVertex.setY(rowTrackIndexProng2.pvRefitY());
        primaryVertex.setZ(rowTrackIndexProng2.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexProng2.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexProng2.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexProng2.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexProng2.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexProng2.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexProng2.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov();
      }
      registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
      registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
      registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
      registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
      registry.fill(HIST("hDcaXYProngs"), track0.pt(), impactParameter0.getY() * toMicrometers);
      registry.fill(HIST("hDcaXYProngs"), track1.pt(), impactParameter1.getY() * toMicrometers);
      registry.fill(HIST("hDcaZProngs"), track0.pt(), impactParameter0.getZ() * toMicrometers);
      registry.fill(HIST("hDcaZProngs"), track1.pt(), impactParameter1.getZ() * toMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      auto indexCollision = collision.globalIndex();
      uint8_t bitmapProngsContributorsPV = 0;
      if (indexCollision == track0.collisionId() && track0.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 0);
      }
      if (indexCollision == track1.collisionId() && track1.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 1);
      }
      uint8_t nProngsContributorsPV = hf_trkcandsel::countOnesInBinary(bitmapProngsContributorsPV);

      // fill candidate table rows
      rowCandidateBase(indexCollision,
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2PCA,
                       pvec0[0], pvec0[1], pvec0[2],
                       pvec1[0], pvec1[1], pvec1[2],
                       impactParameter0.getY(), impactParameter1.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                       impactParameter0.getZ(), impactParameter1.getZ(),
                       std::sqrt(impactParameter0.getSigmaZ2()), std::sqrt(impactParameter1.getSigmaZ2()),
                       rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id(), nProngsContributorsPV, bitmapProngsContributorsPV,
                       rowTrackIndexProng2.hfflag());

      // fill candidate prong PID rows
      fillProngPid<HfProngSpecies::Pion>(track0, rowProng0PidPi);
      fillProngPid<HfProngSpecies::Kaon>(track0, rowProng0PidKa);
      fillProngPid<HfProngSpecies::Pion>(track1, rowProng1PidPi);
      fillProngPid<HfProngSpecies::Kaon>(track1, rowProng1PidKa);

      // fill histograms
      if (fillHistograms) {
        // calculate invariant masses
        auto arrayMomenta = std::array{pvec0, pvec1};
        massPiK = RecoDecay::m(arrayMomenta, std::array{massPi, massK});
        massKPi = RecoDecay::m(arrayMomenta, std::array{massK, massPi});
        registry.fill(HIST("hMass2"), massPiK);
        registry.fill(HIST("hMass2"), massKPi);
      }
    }
  }

  template <bool doPvRefit, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll, typename CandType, typename TTracks>
  void runCreator2ProngWithKFParticle(Coll const&,
                                      CandType const& rowsTrackIndexProng2,
                                      TTracks const&,
                                      aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {

    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng2.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      auto track0 = rowTrackIndexProng2.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexProng2.template prong1_as<TTracks>();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, nullptr, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      float covMatrixPV[6];

      KFParticle::SetField(bz);
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);

      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        kfpVertex.SetXYZ(rowTrackIndexProng2.pvRefitX(), rowTrackIndexProng2.pvRefitY(), rowTrackIndexProng2.pvRefitZ());
        // covariance matrix
        kfpVertex.SetCovarianceMatrix(rowTrackIndexProng2.pvRefitSigmaX2(), rowTrackIndexProng2.pvRefitSigmaXY(), rowTrackIndexProng2.pvRefitSigmaY2(), rowTrackIndexProng2.pvRefitSigmaXZ(), rowTrackIndexProng2.pvRefitSigmaYZ(), rowTrackIndexProng2.pvRefitSigmaZ2());
      }
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle kfpV(kfpVertex);
      registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
      registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
      registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
      registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);

      KFPTrack kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(track1);

      KFParticle kfPosPion(kfpTrack0, kPiPlus);
      KFParticle kfNegPion(kfpTrack1, kPiPlus);
      KFParticle kfPosKaon(kfpTrack0, kKPlus);
      KFParticle kfNegKaon(kfpTrack1, kKPlus);

      float impactParameter0XY = 0., errImpactParameter0XY = 0., impactParameter1XY = 0., errImpactParameter1XY = 0.;
      if (!kfPosPion.GetDistanceFromVertexXY(kfpV, impactParameter0XY, errImpactParameter0XY)) {
        registry.fill(HIST("hDcaXYProngs"), track0.pt(), impactParameter0XY * toMicrometers);
        registry.fill(HIST("hDcaZProngs"), track0.pt(), std::sqrt(kfPosPion.GetDistanceFromVertex(kfpV) * kfPosPion.GetDistanceFromVertex(kfpV) - impactParameter0XY * impactParameter0XY) * toMicrometers);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track0.pt(), -999.f);
        registry.fill(HIST("hDcaZProngs"), track0.pt(), -999.f);
      }
      if (!kfNegPion.GetDistanceFromVertexXY(kfpV, impactParameter1XY, errImpactParameter1XY)) {
        registry.fill(HIST("hDcaXYProngs"), track1.pt(), impactParameter1XY * toMicrometers);
        registry.fill(HIST("hDcaZProngs"), track1.pt(), std::sqrt(kfNegPion.GetDistanceFromVertex(kfpV) * kfNegPion.GetDistanceFromVertex(kfpV) - impactParameter1XY * impactParameter1XY) * toMicrometers);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track1.pt(), -999.f);
        registry.fill(HIST("hDcaZProngs"), track1.pt(), -999.f);
      }

      KFParticle kfCandD0;
      const KFParticle* kfDaughtersD0[2] = {&kfPosPion, &kfNegKaon};
      kfCandD0.SetConstructMethod(2);
      kfCandD0.Construct(kfDaughtersD0, 2);
      KFParticle kfCandD0bar;
      const KFParticle* kfDaughtersD0bar[2] = {&kfNegPion, &kfPosKaon};
      kfCandD0bar.SetConstructMethod(2);
      kfCandD0bar.Construct(kfDaughtersD0bar, 2);

      auto massD0 = kfCandD0.GetMass();
      auto massD0bar = kfCandD0bar.GetMass();

      registry.fill(HIST("hCovSVXX"), kfCandD0.Covariance(0, 0));
      registry.fill(HIST("hCovSVYY"), kfCandD0.Covariance(1, 1));
      registry.fill(HIST("hCovSVXZ"), kfCandD0.Covariance(2, 0));
      registry.fill(HIST("hCovSVZZ"), kfCandD0.Covariance(2, 2));

      auto covMatrixSV = kfCandD0.CovarianceMatrix();

      double phi, theta;
      getPointDirection(std::array{kfpV.GetX(), kfpV.GetY(), kfpV.GetZ()}, std::array{kfCandD0.GetX(), kfCandD0.GetY(), kfCandD0.GetZ()}, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      float topolChi2PerNdfD0 = -999.;
      KFParticle kfCandD0Topol2PV;
      if (constrainKfToPv) {
        kfCandD0Topol2PV = kfCandD0;
        kfCandD0Topol2PV.SetProductionVertex(kfpV);
        topolChi2PerNdfD0 = kfCandD0Topol2PV.GetChi2() / kfCandD0Topol2PV.GetNDF();
      }

      auto indexCollision = collision.globalIndex();
      uint8_t bitmapProngsContributorsPV = 0;
      if (indexCollision == track0.collisionId() && track0.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 0);
      }
      if (indexCollision == track1.collisionId() && track1.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 1);
      }
      uint8_t nProngsContributorsPV = hf_trkcandsel::countOnesInBinary(bitmapProngsContributorsPV);

      // fill candidate table rows
      rowCandidateBase(indexCollision,
                       kfpV.GetX(), kfpV.GetY(), kfpV.GetZ(),
                       kfCandD0.GetX(), kfCandD0.GetY(), kfCandD0.GetZ(),
                       errorDecayLength, errorDecayLengthXY,   // TODO: much different from the DCAFitterN one
                       kfCandD0.GetChi2() / kfCandD0.GetNDF(), // TODO: to make sure it should be chi2 only or chi2/ndf, much different from the DCAFitterN one
                       kfPosPion.GetPx(), kfPosPion.GetPy(), kfPosPion.GetPz(),
                       kfNegKaon.GetPx(), kfNegKaon.GetPy(), kfNegKaon.GetPz(),
                       impactParameter0XY, impactParameter1XY,
                       errImpactParameter0XY, errImpactParameter1XY,
                       0.f, 0.f,
                       0.f, 0.f,
                       rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id(), nProngsContributorsPV, bitmapProngsContributorsPV,
                       rowTrackIndexProng2.hfflag());

      // fill candidate prong PID rows
      fillProngPid<HfProngSpecies::Pion>(track0, rowProng0PidPi);
      fillProngPid<HfProngSpecies::Kaon>(track0, rowProng0PidKa);
      fillProngPid<HfProngSpecies::Pion>(track1, rowProng1PidPi);
      fillProngPid<HfProngSpecies::Kaon>(track1, rowProng1PidKa);

      // fill KF info
      rowCandidateKF(topolChi2PerNdfD0,
                     massD0, massD0bar);

      // fill histograms
      if (fillHistograms) {
        registry.fill(HIST("hMass2"), massD0);
        registry.fill(HIST("hMass2"), massD0bar);
      }
    }
  }

  ///////////////////////////////////
  ///                             ///
  ///   No centrality selection   ///
  ///                             ///
  ///////////////////////////////////

  /// @brief process function using DCA fitter w/ PV refit and w/o centrality selections
  void processPvRefitWithDCAFitterN(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                    TracksWCovExtraPidPiKa const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter w/ PV refit and w/o centrality selections", false);

  /// @brief process function using DCA fitter w/o PV refit and w/o centrality selections
  void processNoPvRefitWithDCAFitterN(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      aod::Hf2Prongs const& rowsTrackIndexProng2,
                                      TracksWCovExtraPidPiKa const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter w/o PV refit and w/o centrality selections", true);

  /// @brief process function using KFParticle package w/ PV refit and w/o centrality selections
  void processPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                    TracksWCovExtraPidPiKa const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticle, "Run candidate creator using KFParticle package w/ PV refit and w/o centrality selections", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      aod::Hf2Prongs const& rowsTrackIndexProng2,
                                      TracksWCovExtraPidPiKa const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithKFParticle, "Run candidate creator using KFParticle package w/o PV refit and w/o centrality selections", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0C   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function using DCA fitter w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitWithDCAFitterNCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                            TracksWCovExtraPidPiKa const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter w/ PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using DCA fitter w/o PV refit and w/ centrality selection FT0C
  void processNoPvRefitWithDCAFitterNCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter w/o PV refit and w/ centrality selection FT0C", false);

  /// @brief process function using KFParticle package w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                            TracksWCovExtraPidPiKa const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticleCentFT0C, "Run candidate creator using KFParticle package w/ PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithKFParticleCentFT0C, "Run candidate creator using KFParticle package w/o PV refit and w/ centrality selection on FT0C", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0M   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function using DCA fitter w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitWithDCAFitterNCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                            TracksWCovExtraPidPiKa const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter w/ PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using DCA fitter w/o PV refit and w/ centrality selection FT0M
  void processNoPvRefitWithDCAFitterNCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter w/o PV refit and w/ centrality selection FT0M", false);

  /// @brief process function using KFParticle package w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                            TracksWCovExtraPidPiKa const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package w/ PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package w/o PV refit and w/ centrality selection on FT0M", false);

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
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processCollisions, "Collision monitoring - no centrality", true);

  /// @brief process function to monitor collisions - FT0C centrality
  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

  /// @brief process function to monitor collisions - FT0M centrality
  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
};

/// Extends the base table with expression columns.
struct HfCandidateCreator2ProngExpressions {
  Spawns<aod::HfCand2ProngExt> rowCandidateProng2;
  Produces<aod::HfCand2ProngMcRec> rowMcMatchRec;
  Produces<aod::HfCand2ProngMcGen> rowMcMatchGen;

  // Configuration
  Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"};
  Configurable<bool> matchKinkedDecayTopology{"matchKinkedDecayTopology", false, "Match also candidates with tracks that decay with kinked topology"};
  Configurable<bool> matchInteractionsWithMaterial{"matchInteractionsWithMaterial", false, "Match also candidates with tracks that interact with material"};

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring

  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};

  // inspect for which zPvPosMax cut was set for reconstructed
  void init(InitContext& initContext)
  {
    std::array<bool, 3> procCollisions = {doprocessMc, doprocessMcCentFT0C, doprocessMcCentFT0M};
    if (std::accumulate(procCollisions.begin(), procCollisions.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for collision study can be enabled at a time.");
    }

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-2prong") == 0) {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }
    hfEvSelMc.addHistograms(registry); // particles monitoring
  }

  /// Performs MC matching.
  template <o2::hf_centrality::CentralityEstimator centEstimator, typename CCs, typename McCollisions>
  void runCreator2ProngMc(aod::TracksWMc const& tracks,
                          aod::McParticles const& mcParticles,
                          CCs const& collInfos,
                          McCollisions const& mcCollisions,
                          BCsInfo const&)
  {
    rowCandidateProng2->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t nKinkedTracks = 0;
    int8_t nInteractionsWithMaterial = 0;

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (const auto& candidate : *rowCandidateProng2) {
      flag = 0;
      origin = 0;
      auto arrayDaughters = std::array{candidate.prong0_as<aod::TracksWMc>(), candidate.prong1_as<aod::TracksWMc>()};

      // Check whether the particle is from background events. If so, reject it.
      if (rejectBackground) {
        bool fromBkg{false};
        for (const auto& daughter : arrayDaughters) {
          if (daughter.has_mcParticle()) {
            auto mcParticle = daughter.mcParticle();
            if (mcParticle.fromBackgroundEvent()) {
              fromBkg = true;
              break;
            }
          }
        }
        if (fromBkg) {
          rowMcMatchRec(flag, origin, -1.f, 0, 0, 0);
          continue;
        }
      }
      std::vector<int> idxBhadMothers{};

      // D0(bar) → π± K∓
      if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1, &nKinkedTracks, &nInteractionsWithMaterial);
      } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1, &nKinkedTracks);
      } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
        indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign, 1, nullptr, &nInteractionsWithMaterial);
      } else {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign);
      }
      if (indexRec > -1) {
        flag = sign * (1 << DecayType::D0ToPiK);
      }

      // J/ψ → e+ e−
      if (flag == 0) {
        if (matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kElectron, -kElectron}, true, &sign, 1, nullptr, &nInteractionsWithMaterial);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kElectron, -kElectron}, true);
        }
        if (indexRec > -1) {
          flag = 1 << DecayType::JpsiToEE;
        }
      }

      // J/ψ → μ+ μ−
      if (flag == 0) {
        if (matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true, &sign, 1, nullptr, &nInteractionsWithMaterial);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true);
        }
        if (indexRec > -1) {
          flag = 1 << DecayType::JpsiToMuMu;
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMcMatchRec(flag, origin, bHadMother.pt(), bHadMother.pdgCode(), nKinkedTracks, nInteractionsWithMaterial);
      } else {
        rowMcMatchRec(flag, origin, -1.f, 0, nKinkedTracks, nInteractionsWithMaterial);
      }
    }

    for (const auto& mcCollision : mcCollisions) {

      // Slice the particles table to get the particles for the current MC collision
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      // Slice the collisions table to get the collision info for the current MC collision
      float centrality{-1.f};
      uint16_t rejectionMask{0};
      int nSplitColl = 0;
      if constexpr (centEstimator == CentralityEstimator::FT0C) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == CentralityEstimator::FT0M) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        nSplitColl = collSlice.size();
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == CentralityEstimator::None) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      }
      hfEvSelMc.fillHistograms<centEstimator>(mcCollision, rejectionMask, nSplitColl);
      if (rejectionMask != 0) {
        // at least one event selection not satisfied --> reject all particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          rowMcMatchGen(0, 0, -1);
        }
        continue;
      }
      hf_mc_gen::fillMcMatchGen2Prong(mcParticles, mcParticlesPerMcColl, rowMcMatchGen, rejectBackground);
    }
  }

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisionsNoCents const& collInfos,
                 aod::McCollisions const& mcCollisions,
                 BCsInfo const& BCsInfo)
  {
    runCreator2ProngMc<CentralityEstimator::None>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMc, "Process MC - no centrality", false);

  void processMcCentFT0C(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Cs const& collInfos,
                         aod::McCollisions const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreator2ProngMc<CentralityEstimator::FT0C>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMcCentFT0C, "Process MC - FT0c centrality", false);

  void processMcCentFT0M(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Ms const& collInfos,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreator2ProngMc<CentralityEstimator::FT0M>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMcCentFT0M, "Process MC - FT0m centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator2Prong>(cfgc, TaskName{"hf-candidate-creator-2prong"}),                         // o2-linter: disable=name/o2-task (wrong hyphenation)
    adaptAnalysisTask<HfCandidateCreator2ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-2prong-expressions"})}; // o2-linter: disable=name/o2-task (wrong hyphenation)
}
