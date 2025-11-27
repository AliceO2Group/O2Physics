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

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsMcGen.h"
#include "PWGHF/Utils/utilsMcMatching.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <Rtypes.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::hf_decay;
using namespace o2::hf_decay::hf_cand_2prong;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::aod::pid_tpc_tof_utils;

/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfCandidateCreator2Prong {
  Produces<aod::HfCand2ProngBase> rowCandidateBase;
  Produces<aod::HfCand2ProngKF> rowCandidateKF;
  Produces<aod::HfCand2Prong0PidPi> rowProng0PidPi;
  Produces<aod::HfCand2Prong0PidKa> rowProng0PidKa;
  Produces<aod::HfCand2Prong1PidPi> rowProng1PidPi;
  Produces<aod::HfCand2Prong1PidKa> rowProng1PidKa;

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
  ctpRateFetcher mRateFetcher;

  int runNumber{0};
  double bz{0.};

  const float toMicrometers = 10000.; // from cm to µm

  std::shared_ptr<TH1> hCandidates;

  using TracksWCovExtraPidPiKa = soa::Join<aod::TracksWCovExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;

  ConfigurableAxis axisMass{"axisMass", {500, 1.6, 2.1}, "axis for mass (GeV/c^2)"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 8> doprocessDF{doprocessPvRefitWithDCAFitterN, doprocessNoPvRefitWithDCAFitterN,
                                    doprocessPvRefitWithDCAFitterNCentFT0C, doprocessNoPvRefitWithDCAFitterNCentFT0C,
                                    doprocessPvRefitWithDCAFitterNCentFT0M, doprocessNoPvRefitWithDCAFitterNCentFT0M, doprocessPvRefitWithDCAFitterNUpc, doprocessNoPvRefitWithDCAFitterNUpc};
    std::array<bool, 8> doprocessKF{doprocessPvRefitWithKFParticle, doprocessNoPvRefitWithKFParticle,
                                    doprocessPvRefitWithKFParticleCentFT0C, doprocessNoPvRefitWithKFParticleCentFT0C,
                                    doprocessPvRefitWithKFParticleCentFT0M, doprocessNoPvRefitWithKFParticleCentFT0M, doprocessPvRefitWithKFParticleUpc, doprocessNoPvRefitWithKFParticleUpc};
    if ((std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) + std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    std::array<bool, 4> processesCollisions = {doprocessCollisions, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M, doprocessCollisionsUpc};
    const int nProcessesCollisions = std::accumulate(processesCollisions.begin(), processesCollisions.end(), 0);
    std::array<bool, 5> processesUpc = {doprocessPvRefitWithDCAFitterNUpc, doprocessNoPvRefitWithDCAFitterNUpc, doprocessPvRefitWithKFParticleUpc, doprocessNoPvRefitWithKFParticleUpc, doprocessCollisionsUpc};
    const int nProcessesUpc = std::accumulate(processesUpc.begin(), processesUpc.end(), 0);
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
      if ((doprocessPvRefitWithDCAFitterNUpc || doprocessNoPvRefitWithDCAFitterNUpc || doprocessPvRefitWithKFParticleUpc || doprocessNoPvRefitWithKFParticleUpc) && !doprocessCollisionsUpc) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsUpc\"?");
      }
    }
    if (nProcessesUpc > 0 && isRun2) {
      LOGP(fatal, "Process function for UPC is only available in Run 3!");
    }

    // histograms
    registry.add("hMass2", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}});
    registry.add("hMassEE", "2-prong candidates;inv. mass (e e) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}});
    registry.add("hMassMuMu", "2-prong candidates;inv. mass (#mu #mu) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMass}});
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

    // init HF event selection helper
    hfEvSel.init(registry);

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

  template <bool DoPvRefit, bool ApplyUpcSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename Coll, typename CandType, typename TTracks, typename BCsType>
  void runCreator2ProngWithDCAFitterN(Coll const&,
                                      CandType const& rowsTrackIndexProng2,
                                      TTracks const&,
                                      BCsType const& bcs)
  {
    // loop over pairs of track indices
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {

      /// reject candidates not satisfying the event selections
      auto collision = rowTrackIndexProng2.template collision_as<Coll>();
      float centrality{-1.f};
      o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
      if constexpr (ApplyUpcSel) {
        rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<true, CentEstimator, BCsType>(collision, centrality, ccdb, registry, bcs);
      } else {
        rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, BCsType>(collision, centrality, ccdb, registry);
      }
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
      auto bc = collision.template bc_as<BCsType>();
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
      std::array<float, 3> pvec0{};
      std::array<float, 3> pvec1{};
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (DoPvRefit) {
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
      const auto nProngsContributorsPV = hf_trkcandsel::countOnesInBinary(bitmapProngsContributorsPV);

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
        const auto arrayMomenta = std::array{pvec0, pvec1};
        const auto massPiK = RecoDecay::m(arrayMomenta, std::array{MassPiPlus, MassKPlus});
        const auto massKPi = RecoDecay::m(arrayMomenta, std::array{MassKPlus, MassPiPlus});
        const auto massEE = RecoDecay::m(arrayMomenta, std::array{MassElectron, MassElectron});
        const auto massMuMu = RecoDecay::m(arrayMomenta, std::array{MassMuon, MassMuon});
        registry.fill(HIST("hMass2"), massPiK);
        registry.fill(HIST("hMass2"), massKPi);
        registry.fill(HIST("hMassEE"), massEE);
        registry.fill(HIST("hMassMuMu"), massMuMu);
      }
    }
  }

  template <bool DoPvRefit, bool ApplyUpcSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename Coll, typename CandType, typename TTracks, typename BCsType>
  void runCreator2ProngWithKFParticle(Coll const&,
                                      CandType const& rowsTrackIndexProng2,
                                      TTracks const&,
                                      BCsType const& bcs)
  {

    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng2.template collision_as<Coll>();
      float centrality{-1.f};
      o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
      if constexpr (ApplyUpcSel) {
        rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<true, CentEstimator, BCsType>(collision, centrality, ccdb, registry, bcs);
      } else {
        rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, BCsType>(collision, centrality, ccdb, registry);
      }
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      auto track0 = rowTrackIndexProng2.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexProng2.template prong1_as<TTracks>();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<BCsType>();
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

      if constexpr (DoPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        kfpVertex.SetXYZ(rowTrackIndexProng2.pvRefitX(), rowTrackIndexProng2.pvRefitY(), rowTrackIndexProng2.pvRefitZ());
        // covariance matrix
        kfpVertex.SetCovarianceMatrix(rowTrackIndexProng2.pvRefitSigmaX2(), rowTrackIndexProng2.pvRefitSigmaXY(), rowTrackIndexProng2.pvRefitSigmaY2(), rowTrackIndexProng2.pvRefitSigmaXZ(), rowTrackIndexProng2.pvRefitSigmaYZ(), rowTrackIndexProng2.pvRefitSigmaZ2());
      }
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle const kfpV(kfpVertex);
      registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
      registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
      registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
      registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);

      KFPTrack const kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack const kfpTrack1 = createKFPTrackFromTrack(track1);

      KFParticle const kfPosPion(kfpTrack0, kPiPlus);
      KFParticle const kfNegPion(kfpTrack1, kPiPlus);
      KFParticle const kfPosKaon(kfpTrack0, kKPlus);
      KFParticle const kfNegKaon(kfpTrack1, kKPlus);

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

      auto* covMatrixSV = kfCandD0.CovarianceMatrix();

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
      const auto nProngsContributorsPV = hf_trkcandsel::countOnesInBinary(bitmapProngsContributorsPV);

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
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter w/ PV refit and w/o centrality selections", false);

  /// @brief process function using DCA fitter w/o PV refit and w/o centrality selections
  void processNoPvRefitWithDCAFitterN(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      aod::Hf2Prongs const& rowsTrackIndexProng2,
                                      TracksWCovExtraPidPiKa const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter w/o PV refit and w/o centrality selections", true);

  /// @brief process function using KFParticle package w/ PV refit and w/o centrality selections
  void processPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                    TracksWCovExtraPidPiKa const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticle, "Run candidate creator using KFParticle package w/ PV refit and w/o centrality selections", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      aod::Hf2Prongs const& rowsTrackIndexProng2,
                                      TracksWCovExtraPidPiKa const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
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
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter w/ PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using DCA fitter w/o PV refit and w/ centrality selection FT0C
  void processNoPvRefitWithDCAFitterNCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter w/o PV refit and w/ centrality selection FT0C", false);

  /// @brief process function using KFParticle package w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                            TracksWCovExtraPidPiKa const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticleCentFT0C, "Run candidate creator using KFParticle package w/ PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
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
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter w/ PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using DCA fitter w/o PV refit and w/ centrality selection FT0M
  void processNoPvRefitWithDCAFitterNCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter w/o PV refit and w/ centrality selection FT0M", false);

  /// @brief process function using KFParticle package w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                            TracksWCovExtraPidPiKa const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package w/ PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::Hf2Prongs const& rowsTrackIndexProng2,
                                              TracksWCovExtraPidPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package w/o PV refit and w/ centrality selection on FT0M", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on UPC   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function using DCA fitter w/ PV refit and w/ centrality selection on UPC
  void processPvRefitWithDCAFitterNUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                       soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                       TracksWCovExtraPidPiKa const& tracks,
                                       aod::BcFullInfos const& bcWithTimeStamps,
                                       aod::FT0s const& /*ft0s*/,
                                       aod::FV0As const& /*fv0as*/,
                                       aod::FDDs const& /*fdds*/,
                                       aod::Zdcs const& /*zdcs*/)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ true, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithDCAFitterNUpc, "Run candidate creator using DCA fitter w/ PV refit and w/ centrality selection on UltraPeripheral Collision", false);

  /// @brief process function using DCA fitter w/o PV refit and w/ centrality selection UPC
  void processNoPvRefitWithDCAFitterNUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                         aod::Hf2Prongs const& rowsTrackIndexProng2,
                                         TracksWCovExtraPidPiKa const& tracks,
                                         aod::BcFullInfos const& bcWithTimeStamps,
                                         aod::FT0s const& /*ft0s*/,
                                         aod::FV0As const& /*fv0as*/,
                                         aod::FDDs const& /*fdds*/,
                                         aod::Zdcs const& /*zdcs*/)
  {
    runCreator2ProngWithDCAFitterN</*doPvRefit*/ false, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithDCAFitterNUpc, "Run candidate creator using DCA fitter w/o PV refit and w/ centrality selection UltraPeripheral Collision", false);

  /// @brief process function using KFParticle package w/ PV refit and w/ centrality selection on UPC
  void processPvRefitWithKFParticleUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                       soa::Join<aod::Hf2Prongs, aod::HfPvRefit2Prong> const& rowsTrackIndexProng2,
                                       TracksWCovExtraPidPiKa const& tracks,
                                       aod::BcFullInfos const& bcWithTimeStamps,
                                       aod::FT0s const& /*ft0s*/,
                                       aod::FV0As const& /*fv0as*/,
                                       aod::FDDs const& /*fdds*/,
                                       aod::Zdcs const& /*zdcs*/)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ true, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processPvRefitWithKFParticleUpc, "Run candidate creator using KFParticle package w/ PV refit and w/ centrality selection on UltraPeripheral Collision", false);

  /// @brief process function using KFParticle package w/o PV refit and w/o centrality selections on UPC
  void processNoPvRefitWithKFParticleUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                         aod::Hf2Prongs const& rowsTrackIndexProng2,
                                         TracksWCovExtraPidPiKa const& tracks,
                                         aod::BcFullInfos const& bcWithTimeStamps,
                                         aod::FT0s const& /*ft0s*/,
                                         aod::FV0As const& /*fv0as*/,
                                         aod::FDDs const& /*fdds*/,
                                         aod::Zdcs const& /*zdcs*/)
  {
    runCreator2ProngWithKFParticle</*doPvRefit*/ false, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng2, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processNoPvRefitWithKFParticleUpc, "Run candidate creator using KFParticle package w/o PV refit and w/ centrality selection on UltraPeripheral Collision", false);

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
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), hfEvSel.irSource, true); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

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
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), hfEvSel.irSource, true); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

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
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), hfEvSel.irSource, true); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);

  /// @brief process function to monitor collisions - UPC collision
  void processCollisionsUpc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                            aod::BcFullInfos const& bcs,
                            aod::FT0s const& /*ft0s*/,
                            aod::FV0As const& /*fv0as*/,
                            aod::FDDs const& /*fdds*/,
                            aod::Zdcs const& /*zdcs*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<true, CentralityEstimator::None, aod::BcFullInfos>(collision, centrality, ccdb, registry, bcs);
      const auto bc = collision.template foundBC_as<aod::BcFullInfos>();
      const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), hfEvSel.irSource, true); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator2Prong, processCollisionsUpc, "Collision monitoring - UPC", false);
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
  Configurable<bool> matchCorrelatedBackground{"matchCorrelatedBackground", false, "Match correlated background candidates"};

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
      if (device.name == "hf-candidate-creator-2prong") {
        // init HF event selection helper
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  /// Performs MC matching.
  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename CCs, typename McCollisions>
  void runCreator2ProngMc(aod::TracksWMc const& tracks,
                          aod::McParticles const& mcParticles,
                          CCs const& collInfos,
                          McCollisions const& mcCollisions,
                          BCsInfo const&)
  {
    rowCandidateProng2->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flagChannelMain = 0;
    int8_t flagChannelResonant = 0;
    int8_t origin = 0;
    int8_t nKinkedTracks = 0;
    int8_t nInteractionsWithMaterial = 0;
    constexpr std::size_t NDaughtersResonant{2u};

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (const auto& candidate : *rowCandidateProng2) {
      flagChannelMain = 0;
      flagChannelResonant = 0;
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
          rowMcMatchRec(flagChannelMain, origin, flagChannelResonant, -1.f, 0, 0, 0);
          continue;
        }
      }
      std::vector<int> idxBhadMothers{};

      if (matchCorrelatedBackground) {
        indexRec = -1; // Index of the matched reconstructed candidate
        constexpr int FinalStateDepth = 2;
        constexpr int ResoDepth = 1;

        // D0(bar) → π+ K−, π+ K− π0, π+ π−, π+ π− π0, K+ K−
        for (const auto& [channelMain, finalState] : daughtersD0Main) {
          std::array<int, 2> const arrPdgDaughtersMain2Prongs = std::array{finalState[0], finalState[1]};
          if (finalState.size() == 3) { // o2-linter: disable=magic-number (partially reconstructed 3-prong decays)
            if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, true, true, true>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth, &nKinkedTracks, &nInteractionsWithMaterial);
            } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, true, true, false>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth, &nKinkedTracks);
            } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, true, false, true>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth, nullptr, &nInteractionsWithMaterial);
            } else {
              indexRec = RecoDecay::getMatchedMCRec<false, false, true, false, false>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth);
            }

            if (indexRec > -1) {
              auto motherParticle = mcParticles.rawIteratorAt(indexRec);
              std::array<int, 3> arrPdgDaughtersMain3Prongs = std::array{finalState[0], finalState[1], finalState[2]};
              flipPdgSign(motherParticle.pdgCode(), +kPi0, arrPdgDaughtersMain3Prongs);
              if (!RecoDecay::isMatchedMCGen(mcParticles, motherParticle, Pdg::kD0, arrPdgDaughtersMain3Prongs, true, &sign, FinalStateDepth)) {
                indexRec = -1; // Reset indexRec if the generated decay does not match the reconstructed one
              }
            }
          } else if (finalState.size() == 2) { // o2-linter: disable=magic-number (fully reconstructed 2-prong decays)
            if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth, &nKinkedTracks, &nInteractionsWithMaterial);
            } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth, &nKinkedTracks);
            } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth, nullptr, &nInteractionsWithMaterial);
            } else {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, false>(mcParticles, arrayDaughters, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, FinalStateDepth);
            }
          } else {
            LOG(fatal) << "Final state size not supported: " << finalState.size();
            return;
          }
          if (indexRec > -1) {
            flagChannelMain = sign * channelMain;

            // Flag the resonant decay channel
            std::vector<int> arrResoDaughIndex = {};
            RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrResoDaughIndex, std::array{0}, ResoDepth);
            std::array<int, NDaughtersResonant> arrPdgDaughters = {};
            if (arrResoDaughIndex.size() == NDaughtersResonant) {
              for (auto iProng = 0u; iProng < arrResoDaughIndex.size(); ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrResoDaughIndex[iProng]);
                arrPdgDaughters[iProng] = daughI.pdgCode();
              }
              flagChannelResonant = getDecayChannelResonant(Pdg::kD0, arrPdgDaughters);
            }
            break;
          }
        }
      } else {
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
          flagChannelMain = sign * DecayChannelMain::D0ToPiK;
        }

        // J/ψ → e+ e−
        if (flagChannelMain == 0) {
          if (matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kElectron, +kPositron}, true, &sign, 1, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kElectron, +kPositron}, true);
          }
          if (indexRec > -1) {
            flagChannelMain = DecayChannelMain::JpsiToEE;
          }
        }

        // J/ψ → μ+ μ−
        if (flagChannelMain == 0) {
          if (matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kMuonMinus, +kMuonPlus}, true, &sign, 1, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kJPsi, std::array{+kMuonMinus, +kMuonPlus}, true);
          }
          if (indexRec > -1) {
            flagChannelMain = DecayChannelMain::JpsiToMuMu;
          }
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flagChannelMain != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMcMatchRec(flagChannelMain, origin, flagChannelResonant, bHadMother.pt(), bHadMother.pdgCode(), nKinkedTracks, nInteractionsWithMaterial);
      } else {
        rowMcMatchRec(flagChannelMain, origin, flagChannelResonant, -1.f, 0, nKinkedTracks, nInteractionsWithMaterial);
      }
    }

    for (const auto& mcCollision : mcCollisions) {

      // Slice the particles table to get the particles for the current MC collision
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      // Slice the collisions table to get the collision info for the current MC collision
      float centrality{-1.f};
      o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
      int nSplitColl = 0;
      if constexpr (CentEstimator == CentralityEstimator::FT0C) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (CentEstimator == CentralityEstimator::FT0M) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        nSplitColl = collSlice.size();
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (CentEstimator == CentralityEstimator::None) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
      }
      hfEvSelMc.fillHistograms<CentEstimator>(mcCollision, rejectionMask, nSplitColl);
      if (rejectionMask != 0) {
        // at least one event selection not satisfied --> reject all particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          rowMcMatchGen(0, 0, 0, -1);
        }
        continue;
      }
      hf_mc_gen::fillMcMatchGen2Prong(mcParticles, mcParticlesPerMcColl, rowMcMatchGen, rejectBackground, matchCorrelatedBackground);
    }
  }

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisionsNoCents const& collInfos,
                 aod::McCollisions const& mcCollisions,
                 BCsInfo const& bcsInfo)
  {
    runCreator2ProngMc<CentralityEstimator::None>(tracks, mcParticles, collInfos, mcCollisions, bcsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMc, "Process MC - no centrality", false);

  void processMcCentFT0C(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Cs const& collInfos,
                         aod::McCollisions const& mcCollisions,
                         BCsInfo const& bcsInfo)
  {
    runCreator2ProngMc<CentralityEstimator::FT0C>(tracks, mcParticles, collInfos, mcCollisions, bcsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMcCentFT0C, "Process MC - FT0c centrality", false);

  void processMcCentFT0M(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Ms const& collInfos,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         BCsInfo const& bcsInfo)
  {
    runCreator2ProngMc<CentralityEstimator::FT0M>(tracks, mcParticles, collInfos, mcCollisions, bcsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator2ProngExpressions, processMcCentFT0M, "Process MC - FT0m centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator2Prong>(cfgc, TaskName{"hf-candidate-creator-2prong"}),                         // o2-linter: disable=name/o2-task (wrong hyphenation)
    adaptAnalysisTask<HfCandidateCreator2ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-2prong-expressions"})}; // o2-linter: disable=name/o2-task (wrong hyphenation)
}
