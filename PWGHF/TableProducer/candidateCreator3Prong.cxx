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

/// \file candidateCreator3Prong.cxx
/// \brief Reconstruction of heavy-flavour 3-prong decay candidates
/// \note Extended from candidateCreator2Prong.cxx
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

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

#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
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
#include <cstdlib>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::hf_decay;
using namespace o2::hf_decay::hf_cand_3prong;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pid_tpc_tof_utils;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreator3Prong {
  Produces<aod::HfCand3ProngBase> rowCandidateBase;
  Produces<aod::HfCand3ProngKF> rowCandidateKF;
  Produces<aod::HfCand3Prong0PidPi> rowProng0PidPi;
  Produces<aod::HfCand3Prong0PidKa> rowProng0PidKa;
  Produces<aod::HfCand3Prong0PidPr> rowProng0PidPr;
  Produces<aod::HfCand3Prong0PidDe> rowProng0PidDe;
  Produces<aod::HfCand3Prong1PidPi> rowProng1PidPi;
  Produces<aod::HfCand3Prong1PidKa> rowProng1PidKa;
  Produces<aod::HfCand3Prong1PidPr> rowProng1PidPr;
  Produces<aod::HfCand3Prong1PidDe> rowProng1PidDe;
  Produces<aod::HfCand3Prong2PidPi> rowProng2PidPi;
  Produces<aod::HfCand3Prong2PidKa> rowProng2PidKa;
  Produces<aod::HfCand3Prong2PidPr> rowProng2PidPr;
  Produces<aod::HfCand3Prong2PidDe> rowProng2PidDe;

  // vertexing
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
  // flags to enable creation for different particle species separately
  Configurable<bool> createDplus{"createDplus", false, "enable D+/- candidate creation"};
  Configurable<bool> createDs{"createDs", false, "enable Ds+/- candidate creation"};
  Configurable<bool> createLc{"createLc", false, "enable Lc+/- candidate creation"};
  Configurable<bool> createXic{"createXic", false, "enable Xic+/- candidate creation"};
  Configurable<bool> createCd{"createCd", false, "enable Cd candidate creation"};
  // KF
  Configurable<bool> applyTopoConstraint{"applyTopoConstraint", false, "apply origin from PV hypothesis for created candidate, works only in KF mode"};
  Configurable<bool> applyInvMassConstraint{"applyInvMassConstraint", false, "apply particle type hypothesis to recalculate created candidate's momentum, works only in KF mode"};

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<3> df; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  int runNumber{0};
  double bz{0.};

  constexpr static float CentiToMicro{10000.f}; // from cm to µm
  constexpr static float UndefValueFloat{-999.f};

  using FilteredHf3Prongs = soa::Filtered<aod::Hf3Prongs>;
  using FilteredPvRefitHf3Prongs = soa::Filtered<soa::Join<aod::Hf3Prongs, aod::HfPvRefit3Prong>>;
  using TracksWCovExtraPidPiKaPrDe = soa::Join<aod::TracksWCovExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr, aod::TracksPidDe, aod::PidTpcTofFullDe>;

  // filter candidates
  Filter filterSelected3Prongs = (createDplus && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(DecayType::DplusToPiKPi))) != static_cast<uint8_t>(0)) || (createDs && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(DecayType::DsToKKPi))) != static_cast<uint8_t>(0)) || (createLc && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(DecayType::LcToPKPi))) != static_cast<uint8_t>(0)) || (createXic && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(DecayType::XicToPKPi))) != static_cast<uint8_t>(0)) || (createCd && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(DecayType::CdToDeKPi))) != static_cast<uint8_t>(0));

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

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
    std::array<bool, 5> creationFlags = {createDplus, createDs, createLc, createXic, createCd};
    if (std::accumulate(creationFlags.begin(), creationFlags.end(), 0) == 0) {
      LOGP(fatal, "At least one particle specie should be enabled for the creation.");
    }

    if (createLc && createXic && applyInvMassConstraint) {
      LOGP(fatal, "Unable to apply invariant mass constraint due to ambiguity of mass hypothesis: only one of Lc and Xic and Cd can be reconstructed.");
    }

    // histograms
    registry.add("hMass3PKPi", "3-prong candidates;inv. mass (pK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1200, 1.8, 3.0}}});
    registry.add("hMass3PiKP", "3-prong candidates;inv. mass (#pi Kp) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1200, 1.8, 3.0}}});
    registry.add("hMass3PiKPi", "3-prong candidates;inv. mass (#pi K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{600, 1.6, 2.2}}});
    registry.add("hMass3KKPi", "3-prong candidates;inv. mass (KK #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{600, 1.7, 2.3}}});
    registry.add("hMass3PiKK", "3-prong candidates;inv. mass (#pi KK) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{600, 1.7, 2.3}}});
    registry.add("hMass3DeKPi", "3-prong candidates;inv. mass (deK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{2000, 2.5, 4.5}}});
    registry.add("hMass3PiKDe", "3-prong candidates;inv. mass (#pi Kde) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{2000, 2.5, 4.5}}});
    registry.add("hMass2KPi", "2-prong pairs;inv. mass (K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1200, 0.8, 2.0}}});
    registry.add("hMass2PiK", "2-prong pairs;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1200, 0.8, 2.0}}});
    registry.add("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hCovPVYY", "3-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVYY", "3-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hCovPVXZ", "3-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, -1.e-4, 1.e-4}}});
    registry.add("hCovSVXZ", "3-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, -1.e-4, 0.2}}});
    registry.add("hCovPVZZ", "3-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVZZ", "3-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hDcaXYProngs", "DCAxy of 3-prong candidate daughters;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
    registry.add("hDcaZProngs", "DCAz of 3-prong candidate daughters;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", {HistType::kTH2F, {{100, 0., 20.}, {200, -500., 500.}}});
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});

    // init HF event selection helper
    hfEvSel.init(registry, &zorroSummary);

    // Configure DCAFitterN
    // df.setBz(bz);
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(static_cast<float>(maxR));
    df.setMaxDZIni(static_cast<float>(maxDZIni));
    df.setMinParamChange(static_cast<float>(minParamChange));
    df.setMinRelChi2Change(static_cast<float>(minRelChi2Change));
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    /// candidate monitoring
    setLabelHistoCands(hCandidates);
  }

  template <typename TRK>
  void fillProngsPid(TRK const& track0, TRK const& track1, TRK const& track2)
  {
    fillProngPid<HfProngSpecies::Pion>(track0, rowProng0PidPi);
    fillProngPid<HfProngSpecies::Kaon>(track0, rowProng0PidKa);
    fillProngPid<HfProngSpecies::Pion>(track1, rowProng1PidPi);
    fillProngPid<HfProngSpecies::Kaon>(track1, rowProng1PidKa);
    fillProngPid<HfProngSpecies::Pion>(track2, rowProng2PidPi);
    fillProngPid<HfProngSpecies::Kaon>(track2, rowProng2PidKa);

    /// fill proton PID information only if necessary
    if (createLc || createXic) {
      fillProngPid<HfProngSpecies::Proton>(track0, rowProng0PidPr);
      fillProngPid<HfProngSpecies::Proton>(track1, rowProng1PidPr);
      fillProngPid<HfProngSpecies::Proton>(track2, rowProng2PidPr);
    }
    if (createCd) {
      fillProngPid<HfProngSpecies::Deuteron>(track0, rowProng0PidDe);
      fillProngPid<HfProngSpecies::Deuteron>(track1, rowProng1PidDe);
      fillProngPid<HfProngSpecies::Deuteron>(track2, rowProng2PidDe);
    }
  }

  template <bool DoPvRefit, bool ApplyUpcSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename Coll, typename Cand, typename BCsType>
  void runCreator3ProngWithDCAFitterN(Coll const&,
                                      Cand const& rowsTrackIndexProng3,
                                      TracksWCovExtraPidPiKaPrDe const&,
                                      BCsType const& bcs)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng3.template collision_as<Coll>();
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

      auto track0 = rowTrackIndexProng3.template prong0_as<TracksWCovExtraPidPiKaPrDe>();
      auto track1 = rowTrackIndexProng3.template prong1_as<TracksWCovExtraPidPiKaPrDe>();
      auto track2 = rowTrackIndexProng3.template prong2_as<TracksWCovExtraPidPiKaPrDe>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);

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
      df.setBz(static_cast<float>(bz));

      // reconstruct the 3-prong secondary vertex
      hCandidates->Fill(SVFitting::BeforeFit);
      try {
        if (df.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
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
      trackParVar0 = df.getTrack(0);
      trackParVar1 = df.getTrack(1);
      trackParVar2 = df.getTrack(2);

      // get track momenta
      std::array<float, 3> pvec0{};
      std::array<float, 3> pvec1{};
      std::array<float, 3> pvec2{};
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);
      trackParVar2.getPxPyPzGlo(pvec2);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (DoPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexProng3.pvRefitX());
        primaryVertex.setY(rowTrackIndexProng3.pvRefitY());
        primaryVertex.setZ(rowTrackIndexProng3.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexProng3.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexProng3.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexProng3.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexProng3.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexProng3.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexProng3.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov();
      }
      registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
      registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
      registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
      registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      o2::dataformats::DCA impactParameter2;
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
      trackParVar2.propagateToDCA(primaryVertex, bz, &impactParameter2);
      registry.fill(HIST("hDcaXYProngs"), track0.pt(), impactParameter0.getY() * CentiToMicro);
      registry.fill(HIST("hDcaXYProngs"), track1.pt(), impactParameter1.getY() * CentiToMicro);
      registry.fill(HIST("hDcaXYProngs"), track2.pt(), impactParameter2.getY() * CentiToMicro);
      registry.fill(HIST("hDcaZProngs"), track0.pt(), impactParameter0.getZ() * CentiToMicro);
      registry.fill(HIST("hDcaZProngs"), track1.pt(), impactParameter1.getZ() * CentiToMicro);
      registry.fill(HIST("hDcaZProngs"), track2.pt(), impactParameter2.getZ() * CentiToMicro);

      // get uncertainty of the decay length
      double phi{}, theta{};
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
      if (indexCollision == track2.collisionId() && track2.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 2);
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
                       pvec2[0], pvec2[1], pvec2[2],
                       impactParameter0.getY(), impactParameter1.getY(), impactParameter2.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameter2.getSigmaY2()),
                       impactParameter0.getZ(), impactParameter1.getZ(), impactParameter2.getZ(),
                       std::sqrt(impactParameter0.getSigmaZ2()), std::sqrt(impactParameter1.getSigmaZ2()), std::sqrt(impactParameter2.getSigmaZ2()),
                       rowTrackIndexProng3.prong0Id(), rowTrackIndexProng3.prong1Id(), rowTrackIndexProng3.prong2Id(), nProngsContributorsPV, bitmapProngsContributorsPV,
                       rowTrackIndexProng3.hfflag());

      // fill candidate prong PID rows
      fillProngsPid(track0, track1, track2);

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        const auto arrayMomenta = std::array{pvec0, pvec1, pvec2};
        const auto massPKPi = RecoDecay::m(arrayMomenta, std::array{MassProton, MassKPlus, MassPiPlus});
        const auto massPiKP = RecoDecay::m(arrayMomenta, std::array{MassPiPlus, MassKPlus, MassProton});
        const auto massPiKPi = RecoDecay::m(arrayMomenta, std::array{MassPiPlus, MassKPlus, MassPiPlus});
        const auto massKKPi = RecoDecay::m(arrayMomenta, std::array{MassKPlus, MassKPlus, MassPiPlus});
        const auto massPiKK = RecoDecay::m(arrayMomenta, std::array{MassPiPlus, MassKPlus, MassKPlus});
        const auto massDeKPi = RecoDecay::m(arrayMomenta, std::array{MassDeuteron, MassKPlus, MassPiPlus});
        const auto massPiKDe = RecoDecay::m(arrayMomenta, std::array{MassPiPlus, MassKPlus, MassDeuteron});
        const auto massKPi = RecoDecay::m(std::array{arrayMomenta.at(1), arrayMomenta.at(2)}, std::array{MassKPlus, MassPiPlus});
        const auto massPiK = RecoDecay::m(std::array{arrayMomenta.at(0), arrayMomenta.at(1)}, std::array{MassPiPlus, MassKPlus});
        registry.fill(HIST("hMass3PiKPi"), massPiKPi);
        registry.fill(HIST("hMass3PKPi"), massPKPi);
        registry.fill(HIST("hMass3PiKP"), massPiKP);
        registry.fill(HIST("hMass3KKPi"), massKKPi);
        registry.fill(HIST("hMass3PiKK"), massPiKK);
        registry.fill(HIST("hMass3DeKPi"), massDeKPi);
        registry.fill(HIST("hMass3PiKDe"), massPiKDe);
        registry.fill(HIST("hMass2KPi"), massKPi);
        registry.fill(HIST("hMass2PiK"), massPiK);
      }
    }
  }

  template <bool DoPvRefit, bool ApplyUpcSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename Coll, typename Cand, typename BCsType>
  void runCreator3ProngWithKFParticle(Coll const&,
                                      Cand const& rowsTrackIndexProng3,
                                      TracksWCovExtraPidPiKaPrDe const&,
                                      BCsType const& bcs)
  {
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {
      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng3.template collision_as<Coll>();
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

      auto track0 = rowTrackIndexProng3.template prong0_as<TracksWCovExtraPidPiKaPrDe>();
      auto track1 = rowTrackIndexProng3.template prong1_as<TracksWCovExtraPidPiKaPrDe>();
      auto track2 = rowTrackIndexProng3.template prong2_as<TracksWCovExtraPidPiKaPrDe>();

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

      KFParticle::SetField(static_cast<float>(bz));
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);

      if constexpr (DoPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        kfpVertex.SetXYZ(rowTrackIndexProng3.pvRefitX(), rowTrackIndexProng3.pvRefitY(), rowTrackIndexProng3.pvRefitZ());
        // covariance matrix
        kfpVertex.SetCovarianceMatrix(rowTrackIndexProng3.pvRefitSigmaX2(), rowTrackIndexProng3.pvRefitSigmaXY(), rowTrackIndexProng3.pvRefitSigmaY2(), rowTrackIndexProng3.pvRefitSigmaXZ(), rowTrackIndexProng3.pvRefitSigmaYZ(), rowTrackIndexProng3.pvRefitSigmaZ2());
      }
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle kfpV(kfpVertex);
      registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
      registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
      registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
      registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);

      KFPTrack const kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack const kfpTrack1 = createKFPTrackFromTrack(track1);
      KFPTrack const kfpTrack2 = createKFPTrackFromTrack(track2);

      KFParticle const kfFirstProton(kfpTrack0, kProton);
      KFParticle const kfFirstPion(kfpTrack0, kPiPlus);
      KFParticle const kfFirstKaon(kfpTrack0, kKPlus);
      KFParticle const kfSecondKaon(kfpTrack1, kKPlus);
      KFParticle const kfThirdProton(kfpTrack2, kProton);
      KFParticle const kfThirdPion(kfpTrack2, kPiPlus);
      KFParticle const kfThirdKaon(kfpTrack2, kKPlus);

      float impactParameter0XY = 0., errImpactParameter0XY = 0., impactParameter1XY = 0., errImpactParameter1XY = 0., impactParameter2XY = 0., errImpactParameter2XY = 0.;
      if (!kfFirstProton.GetDistanceFromVertexXY(kfpV, impactParameter0XY, errImpactParameter0XY)) {
        registry.fill(HIST("hDcaXYProngs"), track0.pt(), impactParameter0XY * CentiToMicro);
        registry.fill(HIST("hDcaZProngs"), track0.pt(), std::sqrt(kfFirstProton.GetDistanceFromVertex(kfpV) * kfFirstProton.GetDistanceFromVertex(kfpV) - impactParameter0XY * impactParameter0XY) * CentiToMicro);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track0.pt(), UndefValueFloat);
        registry.fill(HIST("hDcaZProngs"), track0.pt(), UndefValueFloat);
      }
      if (!kfSecondKaon.GetDistanceFromVertexXY(kfpV, impactParameter1XY, errImpactParameter1XY)) {
        registry.fill(HIST("hDcaXYProngs"), track1.pt(), impactParameter1XY * CentiToMicro);
        registry.fill(HIST("hDcaZProngs"), track1.pt(), std::sqrt(kfSecondKaon.GetDistanceFromVertex(kfpV) * kfSecondKaon.GetDistanceFromVertex(kfpV) - impactParameter1XY * impactParameter1XY) * CentiToMicro);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track1.pt(), UndefValueFloat);
        registry.fill(HIST("hDcaZProngs"), track1.pt(), UndefValueFloat);
      }
      if (!kfThirdProton.GetDistanceFromVertexXY(kfpV, impactParameter2XY, errImpactParameter2XY)) {
        registry.fill(HIST("hDcaXYProngs"), track2.pt(), impactParameter2XY * CentiToMicro);
        registry.fill(HIST("hDcaZProngs"), track2.pt(), std::sqrt(kfThirdProton.GetDistanceFromVertex(kfpV) * kfThirdProton.GetDistanceFromVertex(kfpV) - impactParameter2XY * impactParameter2XY) * CentiToMicro);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track2.pt(), UndefValueFloat);
        registry.fill(HIST("hDcaZProngs"), track2.pt(), UndefValueFloat);
      }

      auto [impactParameter0Z, errImpactParameter0Z] = kfCalculateImpactParameterZ(kfFirstProton, kfpV);
      auto [impactParameter1Z, errImpactParameter1Z] = kfCalculateImpactParameterZ(kfSecondKaon, kfpV);
      auto [impactParameter2Z, errImpactParameter2Z] = kfCalculateImpactParameterZ(kfThirdProton, kfpV);

      const float chi2primFirst = kfCalculateChi2ToPrimaryVertex(kfFirstProton, kfpV);
      const float chi2primSecond = kfCalculateChi2ToPrimaryVertex(kfSecondKaon, kfpV);
      const float chi2primThird = kfCalculateChi2ToPrimaryVertex(kfThirdPion, kfpV);

      const float dcaSecondThird = kfCalculateDistanceBetweenParticles(kfSecondKaon, kfThirdPion);
      const float dcaFirstThird = kfCalculateDistanceBetweenParticles(kfFirstProton, kfThirdPion);
      const float dcaFirstSecond = kfCalculateDistanceBetweenParticles(kfFirstProton, kfSecondKaon);

      const float chi2geoSecondThird = kfCalculateChi2geoBetweenParticles(kfSecondKaon, kfThirdPion);
      const float chi2geoFirstThird = kfCalculateChi2geoBetweenParticles(kfFirstProton, kfThirdPion);
      const float chi2geoFirstSecond = kfCalculateChi2geoBetweenParticles(kfFirstProton, kfSecondKaon);

      // Λc± → p± K∓ π±,  Ξc± → p± K∓ π±
      KFParticle kfCandPKPi;
      const KFParticle* kfDaughtersPKPi[3] = {&kfFirstProton, &kfSecondKaon, &kfThirdPion};
      kfCandPKPi.SetConstructMethod(2);
      kfCandPKPi.Construct(kfDaughtersPKPi, 3);
      KFParticle kfCandPiKP;
      const KFParticle* kfDaughtersPiKP[3] = {&kfFirstPion, &kfSecondKaon, &kfThirdProton};
      kfCandPiKP.SetConstructMethod(2);
      kfCandPiKP.Construct(kfDaughtersPiKP, 3);

      // D± → π± K∓ π±
      KFParticle kfCandPiKPi;
      const KFParticle* kfDaughtersPiKPi[3] = {&kfFirstPion, &kfSecondKaon, &kfThirdPion};
      kfCandPiKPi.SetConstructMethod(2);
      kfCandPiKPi.Construct(kfDaughtersPiKPi, 3);

      // Ds± → K± K∓ π±
      KFParticle kfCandKKPi;
      const KFParticle* kfDaughtersKKPi[3] = {&kfFirstKaon, &kfSecondKaon, &kfThirdPion};
      kfCandKKPi.SetConstructMethod(2);
      kfCandKKPi.Construct(kfDaughtersKKPi, 3);
      KFParticle kfCandPiKK;
      const KFParticle* kfDaughtersPiKK[3] = {&kfFirstPion, &kfSecondKaon, &kfThirdKaon};
      kfCandPiKK.SetConstructMethod(2);
      kfCandPiKK.Construct(kfDaughtersPiKK, 3);

      const float chi2topo = kfCalculateChi2ToPrimaryVertex(kfCandPKPi, kfpV);

      if (applyTopoConstraint) { // constraints applied after chi2topo getter - to preserve unbiased value of chi2topo
        for (const auto& kfCand : std::array<KFParticle*, 5>{&kfCandPKPi, &kfCandPiKP, &kfCandPiKPi, &kfCandKKPi, &kfCandPiKK}) {
          kfCand->SetProductionVertex(kfpV);
          kfCand->TransportToDecayVertex();
        }
      }

      KFParticle kfPairKPi;
      const KFParticle* kfDaughtersKPi[3] = {&kfSecondKaon, &kfThirdPion};
      kfPairKPi.SetConstructMethod(2);
      kfPairKPi.Construct(kfDaughtersKPi, 2);

      KFParticle kfPairPiK;
      const KFParticle* kfDaughtersPiK[3] = {&kfFirstPion, &kfSecondKaon};
      kfPairPiK.SetConstructMethod(2);
      kfPairPiK.Construct(kfDaughtersPiK, 2);

      const float massPKPi = kfCandPKPi.GetMass();
      const float massPiKP = kfCandPiKP.GetMass();
      const float massPiKPi = kfCandPiKPi.GetMass();
      const float massKKPi = kfCandKKPi.GetMass();
      const float massPiKK = kfCandPiKK.GetMass();
      const float massKPi = kfPairKPi.GetMass();
      const float massPiK = kfPairPiK.GetMass();

      if (applyInvMassConstraint) { // constraints applied after minv getters - to preserve unbiased values of minv
        kfCandPKPi.SetNonlinearMassConstraint(createLc ? static_cast<float>(MassLambdaCPlus) : static_cast<float>(MassXiCPlus));
        kfCandPiKP.SetNonlinearMassConstraint(createLc ? static_cast<float>(MassLambdaCPlus) : static_cast<float>(MassXiCPlus));
        kfCandPiKPi.SetNonlinearMassConstraint(static_cast<float>(MassDPlus));
        kfCandKKPi.SetNonlinearMassConstraint(static_cast<float>(MassDS));
        kfCandPiKK.SetNonlinearMassConstraint(static_cast<float>(MassDS));
      }

      const float chi2geo = kfCandPKPi.Chi2() / static_cast<float>(kfCandPKPi.NDF());
      const std::pair<float, float> ldl = kfCalculateLdL(kfCandPKPi, kfpV);

      std::array<float, 3> pProng0 = kfCalculateProngMomentumInSecondaryVertex(kfFirstProton, kfCandPiKP);
      std::array<float, 3> pProng1 = kfCalculateProngMomentumInSecondaryVertex(kfSecondKaon, kfCandPiKP);
      std::array<float, 3> pProng2 = kfCalculateProngMomentumInSecondaryVertex(kfThirdPion, kfCandPiKP);

      registry.fill(HIST("hCovSVXX"), kfCandPKPi.Covariance(0, 0));
      registry.fill(HIST("hCovSVYY"), kfCandPKPi.Covariance(1, 1));
      registry.fill(HIST("hCovSVXZ"), kfCandPKPi.Covariance(2, 0));
      registry.fill(HIST("hCovSVZZ"), kfCandPKPi.Covariance(2, 2));

      auto* covMatrixSV = kfCandPKPi.CovarianceMatrix();
      double phi{}, theta{};
      getPointDirection(std::array{kfpV.GetX(), kfpV.GetY(), kfpV.GetZ()}, std::array{kfCandPKPi.GetX(), kfCandPKPi.GetY(), kfCandPKPi.GetZ()}, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      auto indexCollision = collision.globalIndex();
      uint8_t bitmapProngsContributorsPV = 0;
      if (indexCollision == track0.collisionId() && track0.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 0);
      }
      if (indexCollision == track1.collisionId() && track1.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 1);
      }
      if (indexCollision == track2.collisionId() && track2.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 2);
      }
      const auto nProngsContributorsPV = hf_trkcandsel::countOnesInBinary(bitmapProngsContributorsPV);

      // fill candidate table rows
      rowCandidateBase(indexCollision,
                       kfpV.GetX(), kfpV.GetY(), kfpV.GetZ(),
                       kfCandPKPi.GetX(), kfCandPKPi.GetY(), kfCandPKPi.GetZ(),
                       errorDecayLength, errorDecayLengthXY,
                       kfCandPKPi.GetChi2() / static_cast<float>(kfCandPKPi.GetNDF()),
                       pProng0.at(0), pProng0.at(1), pProng0.at(2),
                       pProng1.at(0), pProng1.at(1), pProng1.at(2),
                       pProng2.at(0), pProng2.at(1), pProng2.at(2),
                       impactParameter0XY, impactParameter1XY, impactParameter2XY,
                       errImpactParameter0XY, errImpactParameter1XY, errImpactParameter2XY,
                       impactParameter0Z, impactParameter1Z, impactParameter2Z,
                       errImpactParameter0Z, errImpactParameter1Z, errImpactParameter2Z,
                       rowTrackIndexProng3.prong0Id(), rowTrackIndexProng3.prong1Id(), rowTrackIndexProng3.prong2Id(), nProngsContributorsPV, bitmapProngsContributorsPV,
                       rowTrackIndexProng3.hfflag());

      // fill KF info
      rowCandidateKF(kfCandPKPi.GetErrX(), kfCandPKPi.GetErrY(), kfCandPKPi.GetErrZ(),
                     std::sqrt(kfpV.Covariance(0, 0)), std::sqrt(kfpV.Covariance(1, 1)), std::sqrt(kfpV.Covariance(2, 2)),
                     massPKPi, massPiKP, massPiKPi, massKKPi, massPiKK, massKPi, massPiK,
                     kfCandPKPi.GetPx(), kfCandPKPi.GetPy(), kfCandPKPi.GetPz(),
                     kfCandPKPi.GetErrPx(), kfCandPKPi.GetErrPy(), kfCandPKPi.GetErrPz(),
                     chi2primFirst, chi2primSecond, chi2primThird,
                     dcaSecondThird, dcaFirstThird, dcaFirstSecond,
                     chi2geoSecondThird, chi2geoFirstThird, chi2geoFirstSecond,
                     chi2geo, ldl.first, ldl.second, chi2topo);

      // fill candidate prong PID rows
      fillProngsPid(track0, track1, track2);

      // fill histograms
      if (fillHistograms) {
        registry.fill(HIST("hMass3PiKPi"), massPiKPi);
        registry.fill(HIST("hMass3PKPi"), massPKPi);
        registry.fill(HIST("hMass3PiKP"), massPiKP);
        registry.fill(HIST("hMass3KKPi"), massKKPi);
        registry.fill(HIST("hMass3PiKK"), massPiKK);
        registry.fill(HIST("hMass2KPi"), massKPi);
        registry.fill(HIST("hMass2PiK"), massPiK);
      }
    }
  }

  ///////////////////////////////////
  ///                             ///
  ///   No centrality selection   ///
  ///                             ///
  ///////////////////////////////////

  /// @brief process function using DCA fitter  w/ PV refit and w/o centrality selections
  void processPvRefitWithDCAFitterN(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                    TracksWCovExtraPidPiKaPrDe const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter with PV refit and w/o centrality selections", false);

  /// @brief process function using DCA fitter  w/o PV refit and w/o centrality selections
  void processNoPvRefitWithDCAFitterN(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      FilteredHf3Prongs const& rowsTrackIndexProng3,
                                      TracksWCovExtraPidPiKaPrDe const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter without PV refit and w/o centrality selections", true);

  /// @brief process function using KFParticle package  w/ PV refit and w/o centrality selections
  void processPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                    TracksWCovExtraPidPiKaPrDe const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticle, "Run candidate creator using KFParticle package with PV refit and w/o centrality selections", false);

  /// @brief process function using KFParticle package  w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      FilteredHf3Prongs const& rowsTrackIndexProng3,
                                      TracksWCovExtraPidPiKaPrDe const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithKFParticle, "Run candidate creator using KFParticle package without PV refit and w/o centrality selections", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0C   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function using DCA fitter  w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitWithDCAFitterNCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                            TracksWCovExtraPidPiKaPrDe const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter with PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using DCA fitter  w/o PV refit and  w/ centrality selection on FT0C
  void processNoPvRefitWithDCAFitterNCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              TracksWCovExtraPidPiKaPrDe const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter without PV refit and  w/ centrality selection on FT0C", false);

  /// @brief process function using KFParticle package  w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                            TracksWCovExtraPidPiKaPrDe const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticleCentFT0C, "Run candidate creator using KFParticle package with PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using KFParticle package  w/o PV refit and  w/ centrality selection on FT0C
  void processNoPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              TracksWCovExtraPidPiKaPrDe const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithKFParticleCentFT0C, "Run candidate creator using KFParticle package without PV refit and  w/ centrality selection on FT0C", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0M   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function using DCA fitter  w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitWithDCAFitterNCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                            TracksWCovExtraPidPiKaPrDe const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter with PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using DCA fitter  w/o PV refit and  w/ centrality selection on FT0M
  void processNoPvRefitWithDCAFitterNCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              TracksWCovExtraPidPiKaPrDe const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter without PV refit and  w/ centrality selection on FT0M", false);

  /// @brief process function using KFParticle package  w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                            TracksWCovExtraPidPiKaPrDe const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package with PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using KFParticle package  w/o PV refit and  w/ centrality selection on FT0M
  void processNoPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              TracksWCovExtraPidPiKaPrDe const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package without PV refit and  w/ centrality selection on FT0M", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on UPC   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function using DCA fitter  w/ PV refit and w/ centrality selection on UPC
  void processPvRefitWithDCAFitterNUpc(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                       FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                       TracksWCovExtraPidPiKaPrDe const& tracks,
                                       aod::BcFullInfos const& bcWithTimeStamps,
                                       aod::FT0s const& /*ft0s*/,
                                       aod::FV0As const& /*fv0as*/,
                                       aod::FDDs const& /*fdds*/,
                                       aod::Zdcs const& /*zdcs*/)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterNUpc, "Run candidate creator using DCA fitter with PV refit and w/ centrality selection on UPC", false);

  /// @brief process function using DCA fitter  w/o PV refit and  w/ centrality selection on UPC
  void processNoPvRefitWithDCAFitterNUpc(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                         FilteredHf3Prongs const& rowsTrackIndexProng3,
                                         TracksWCovExtraPidPiKaPrDe const& tracks,
                                         aod::BcFullInfos const& bcWithTimeStamps,
                                         aod::FT0s const& /*ft0s*/,
                                         aod::FV0As const& /*fv0as*/,
                                         aod::FDDs const& /*fdds*/,
                                         aod::Zdcs const& /*zdcs*/)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterNUpc, "Run candidate creator using DCA fitter without PV refit and  w/ centrality selection on UPC", false);

  /// @brief process function using KFParticle package  w/ PV refit and w/ centrality selection on UPC
  void processPvRefitWithKFParticleUpc(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                       FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                       TracksWCovExtraPidPiKaPrDe const& tracks,
                                       aod::BcFullInfos const& bcWithTimeStamps,
                                       aod::FT0s const& /*ft0s*/,
                                       aod::FV0As const& /*fv0as*/,
                                       aod::FDDs const& /*fdds*/,
                                       aod::Zdcs const& /*zdcs*/)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticleUpc, "Run candidate creator using KFParticle package with PV refit and w/ centrality selection on UPC", false);

  /// @brief process function using KFParticle package  w/o PV refit and  w/ centrality selection on UPC
  void processNoPvRefitWithKFParticleUpc(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                         FilteredHf3Prongs const& rowsTrackIndexProng3,
                                         TracksWCovExtraPidPiKaPrDe const& tracks,
                                         aod::BcFullInfos const& bcWithTimeStamps,
                                         aod::FT0s const& /*ft0s*/,
                                         aod::FV0As const& /*fv0as*/,
                                         aod::FDDs const& /*fdds*/,
                                         aod::Zdcs const& /*zdcs*/)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithKFParticleUpc, "Run candidate creator using KFParticle package without PV refit and  w/ centrality selection on UPC", false);

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
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, static_cast<float>(ir));

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisions, "Collision monitoring - no centrality", true);

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
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, static_cast<float>(ir));

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

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
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, static_cast<float>(ir));

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);

  /// @brief process function to monitor collisions - UPC
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
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, static_cast<float>(ir));

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsUpc, "Collision monitoring - UPC", false);
};

/// Extends the base table with expression columns.
struct HfCandidateCreator3ProngExpressions {
  Spawns<aod::HfCand3ProngExt> rowCandidateProng3;
  Produces<aod::HfCand3ProngMcRec> rowMcMatchRec;
  Produces<aod::HfCand3ProngMcGen> rowMcMatchGen;

  // Configuration
  Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"};
  Configurable<bool> matchKinkedDecayTopology{"matchKinkedDecayTopology", false, "Match also candidates with tracks that decay with kinked topology"};
  Configurable<bool> matchInteractionsWithMaterial{"matchInteractionsWithMaterial", false, "Match also candidates with tracks that interact with material"};
  Configurable<bool> matchCorrelatedBackground{"matchCorrelatedBackground", false, "Match correlated background candidates"};
  Configurable<std::vector<int>> pdgMothersCorrelBkg{"pdgMothersCorrelBkg", {Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kLambdaCPlus, Pdg::kXiCPlus}, "PDG codes of the mother particles of correlated background candidates"};

  constexpr static std::size_t NDaughtersResonant{2u};

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring

  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    std::array<bool, 3> procCollisions = {doprocessMc, doprocessMcCentFT0C, doprocessMcCentFT0M};
    if (std::accumulate(procCollisions.begin(), procCollisions.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for collision study can be enabled at a time.");
    }

    // inspect for which particle species the candidates were created and which zPvPosMax cut was set for reconstructed
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-candidate-creator-3prong") {
        // init HF event selection helper
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  /// Performs MC matching.
  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename CCs, typename McCollisions>
  void runCreator3ProngMc(aod::TracksWMc const& tracks,
                          aod::McParticles const& mcParticles,
                          CCs const& collInfos,
                          McCollisions const& mcCollisions,
                          BCsInfo const&)
  {
    rowCandidateProng3->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flagChannelMain = 0;
    int8_t flagChannelResonant = 0;
    int8_t origin = 0;
    int8_t swapping = 0;
    int8_t nKinkedTracks = 0;
    int8_t nInteractionsWithMaterial = 0;
    std::vector<int> arrDaughIndex{};
    std::array<int, NDaughtersResonant> arrPdgDaugResonant{};
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantLcToPKstar0{daughtersLcResonant.at(DecayChannelResonant::LcToPKstar0)};               // Λc± → p± K*
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantLcToDeltaplusplusK{daughtersLcResonant.at(DecayChannelResonant::LcToDeltaplusplusK)}; // Λc± → Δ(1232)±± K∓
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantLcToL1520Pi{daughtersLcResonant.at(DecayChannelResonant::LcToL1520Pi)};               // Λc± → Λ(1520) π±
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantDToPhiPi{daughtersDsResonant.at(DecayChannelResonant::DsToPhiPi)};                    // Ds± → φ π± and D± → φ π±
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantDToKstar0K{daughtersDsResonant.at(DecayChannelResonant::DsToKstar0K)};                // Ds± → anti-K*(892)0 K± and D± → anti-K*(892)0 K±

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (const auto& candidate : *rowCandidateProng3) {
      flagChannelMain = 0;
      flagChannelResonant = 0;
      origin = 0;
      swapping = 0;
      indexRec = -1;
      arrDaughIndex.clear();
      std::vector<int> idxBhadMothers{};
      auto arrayDaughters = std::array{candidate.prong0_as<aod::TracksWMc>(), candidate.prong1_as<aod::TracksWMc>(), candidate.prong2_as<aod::TracksWMc>()};

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
          rowMcMatchRec(flagChannelMain, origin, swapping, flagChannelResonant, -1.f, 0, 0, 0);
          continue;
        }
      }

      if (matchCorrelatedBackground) {
        indexRec = -1;                  // Index of the matched reconstructed candidate
        constexpr int DepthMainMax = 2; // Depth for final state matching
        constexpr int DepthResoMax = 1; // Depth for resonant decay matching

        for (const auto& pdgMother : pdgMothersCorrelBkg.value) {
          int depthMainMax = DepthMainMax;
          if (pdgMother == Pdg::kDStar) {
            depthMainMax = DepthMainMax + 1; // D0 resonant decays are active
          }
          auto finalStates = getDecayChannelsMain(pdgMother);
          for (const auto& [channelMain, finalState] : finalStates) {
            std::array<int, 3> const arrPdgDaughtersMain3Prongs = std::array{finalState[0], finalState[1], finalState[2]};
            if (finalState.size() > 3) { // o2-linter: disable=magic-number (partially reconstructed decays with 4 or 5 final state particles)
              if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                indexRec = RecoDecay::getMatchedMCRec<false, false, true, true, true>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax, &nKinkedTracks, &nInteractionsWithMaterial);
              } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
                indexRec = RecoDecay::getMatchedMCRec<false, false, true, true, false>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax, &nKinkedTracks);
              } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                indexRec = RecoDecay::getMatchedMCRec<false, false, true, false, true>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax, nullptr, &nInteractionsWithMaterial);
              } else {
                indexRec = RecoDecay::getMatchedMCRec<false, false, true, false, false>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax);
              }

              if (indexRec > -1) {
                auto motherParticle = mcParticles.rawIteratorAt(indexRec);
                if (finalState.size() == 4) { // o2-linter: disable=magic-number (Check if the final state has 4 particles)
                  std::array<int, 4> arrPdgDaughtersMain4Prongs = std::array{finalState[0], finalState[1], finalState[2], finalState[3]};
                  flipPdgSign(motherParticle.pdgCode(), +kPi0, arrPdgDaughtersMain4Prongs);
                  if (!RecoDecay::isMatchedMCGen(mcParticles, motherParticle, pdgMother, arrPdgDaughtersMain4Prongs, true, &sign, depthMainMax)) {
                    indexRec = -1; // Reset indexRec if the generated decay does not match the reconstructed one is not matched
                  }
                } else if (finalState.size() == 5) { // o2-linter: disable=magic-number (Check if the final state has 5 particles)
                  std::array<int, 5> arrPdgDaughtersMain5Prongs = std::array{finalState[0], finalState[1], finalState[2], finalState[3], finalState[4]};
                  flipPdgSign(motherParticle.pdgCode(), +kPi0, arrPdgDaughtersMain5Prongs);
                  if (!RecoDecay::isMatchedMCGen(mcParticles, motherParticle, pdgMother, arrPdgDaughtersMain5Prongs, true, &sign, depthMainMax)) {
                    indexRec = -1; // Reset indexRec if the generated decay does not match the reconstructed one is not matched
                  }
                }
              }
            } else if (finalState.size() == 3) { // o2-linter: disable=magic-number (fully reconstructed 3-prong decays)
              if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax, &nKinkedTracks, &nInteractionsWithMaterial);
              } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
                indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax, &nKinkedTracks);
              } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
                indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax, nullptr, &nInteractionsWithMaterial);
              } else {
                indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, false>(mcParticles, arrayDaughters, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax);
              }
            } else {
              LOG(fatal) << "Final state size not supported: " << finalState.size();
              return;
            }
            if (indexRec > -1) {
              flagChannelMain = static_cast<int8_t>(sign * channelMain);

              /// swapping for D+, Ds->Kpipi; Lc, Xic->pKpi
              if (std::abs(flagChannelMain) == DecayChannelMain::DplusToPiKK || std::abs(flagChannelMain) == DecayChannelMain::DsToPiKK || std::abs(flagChannelMain) == DecayChannelMain::LcToPKPi || std::abs(flagChannelMain) == DecayChannelMain::XicToPKPi) {
                if (arrayDaughters[0].has_mcParticle()) {
                  swapping = static_cast<int8_t>(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
                }
              }

              // Flag the resonant decay channel
              std::vector<int> arrResoDaughIndex = {};
              if (pdgMother == Pdg::kDStar) {
                std::vector<int> arrResoDaughIndexDstar = {};
                RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrResoDaughIndexDstar, std::array{0}, DepthResoMax);
                for (const int iDaug : arrResoDaughIndexDstar) { // o2-linter: disable=const-ref-in-for-loop (int elements)
                  auto daughDstar = mcParticles.rawIteratorAt(iDaug);
                  if (std::abs(daughDstar.pdgCode()) == Pdg::kD0 || std::abs(daughDstar.pdgCode()) == Pdg::kDPlus) {
                    RecoDecay::getDaughters(daughDstar, &arrResoDaughIndex, std::array{0}, DepthResoMax);
                    break;
                  }
                }
              } else {
                RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrResoDaughIndex, std::array{0}, DepthResoMax);
              }
              std::array<int, NDaughtersResonant> arrPdgDaughters = {};
              if (arrResoDaughIndex.size() == NDaughtersResonant) {
                for (auto iProng = 0u; iProng < NDaughtersResonant; ++iProng) {
                  auto daughI = mcParticles.rawIteratorAt(arrResoDaughIndex[iProng]);
                  arrPdgDaughters[iProng] = daughI.pdgCode();
                }
                flagChannelResonant = getDecayChannelResonant(pdgMother, arrPdgDaughters);
              }
              break; // Exit loop if a match is found
            }
          }
          if (indexRec > -1) {
            break; // Exit loop if a match is found
          }
        }
      } else {
        // D± → π± K∓ π±
        if (flagChannelMain == 0) {
          auto arrPdgDaughtersDplusToPiKPi{std::array{+kPiPlus, -kKPlus, +kPiPlus}};
          if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDplusToPiKPi, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
          } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDplusToPiKPi, true, &sign, 2, &nKinkedTracks);
          } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDplusToPiKPi, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDplusToPiKPi, true, &sign, 2);
          }
          if (indexRec > -1) {
            flagChannelMain = static_cast<int8_t>(sign * DecayChannelMain::DplusToPiKPi);
          }
        }

        // Ds± → K± K∓ π± and D± → K± K∓ π±
        if (flagChannelMain == 0) {
          auto arrPdgDaughtersDToPiKK{std::array{+kKPlus, -kKPlus, +kPiPlus}};
          bool isDplus = false;
          if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kDS, arrPdgDaughtersDToPiKK, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
          } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kDS, arrPdgDaughtersDToPiKK, true, &sign, 2, &nKinkedTracks);
          } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDS, arrPdgDaughtersDToPiKK, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDS, arrPdgDaughtersDToPiKK, true, &sign, 2);
          }
          if (indexRec == -1) {
            isDplus = true;
            if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDToPiKK, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
            } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDToPiKK, true, &sign, 2, &nKinkedTracks);
            } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
              indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDToPiKK, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
            } else {
              indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDPlus, arrPdgDaughtersDToPiKK, true, &sign, 2);
            }
          }
          if (indexRec > -1) {
            flagChannelMain = sign * (isDplus ? static_cast<int8_t>(DecayChannelMain::DplusToPiKK) : static_cast<int8_t>(DecayChannelMain::DsToPiKK));
            if (arrayDaughters[0].has_mcParticle()) {
              swapping = static_cast<int8_t>(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
            }
            RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrDaughIndex, std::array{0}, 1);
            if (arrDaughIndex.size() == NDaughtersResonant) {
              for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
                arrPdgDaugResonant[iProng] = std::abs(daughI.pdgCode());
              }
              if ((arrPdgDaugResonant[0] == arrPdgDaugResonantDToPhiPi[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToPhiPi[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantDToPhiPi[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToPhiPi[0])) {
                flagChannelResonant = isDplus ? DecayChannelResonant::DplusToPhiPi : DecayChannelResonant::DsToPhiPi;
              } else if ((arrPdgDaugResonant[0] == arrPdgDaugResonantDToKstar0K[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToKstar0K[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantDToKstar0K[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToKstar0K[0])) {
                flagChannelResonant = isDplus ? DecayChannelResonant::DplusToKstar0K : DecayChannelResonant::DsToKstar0K;
              }
            }
          }
        }

        // D* → D0 π → K π π
        if (flagChannelMain == 0) {
          auto arrPdgDaughtersDstarToPiKPi{std::array{+kPiPlus, +kPiPlus, -kKPlus}};
          if (matchKinkedDecayTopology) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDStar, arrPdgDaughtersDstarToPiKPi, true, &sign, 2, &nKinkedTracks);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDStar, arrPdgDaughtersDstarToPiKPi, true, &sign, 2);
          }
          if (indexRec > -1) {
            flagChannelMain = static_cast<int8_t>(sign * DecayChannelMain::DstarToPiKPi);
            flagChannelResonant = 0;
          }
        }

        // Λc± → p± K∓ π±
        if (flagChannelMain == 0) {
          auto arrPdgDaughtersLcToPKPi{std::array{+kProton, -kKPlus, +kPiPlus}};
          if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, arrPdgDaughtersLcToPKPi, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
          } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, arrPdgDaughtersLcToPKPi, true, &sign, 2, &nKinkedTracks);
          } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, arrPdgDaughtersLcToPKPi, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, arrPdgDaughtersLcToPKPi, true, &sign, 2);
          }
          if (indexRec > -1) {
            flagChannelMain = static_cast<int8_t>(sign * DecayChannelMain::LcToPKPi);

            // Flagging the different Λc± → p± K∓ π± decay channels
            if (arrayDaughters[0].has_mcParticle()) {
              swapping = static_cast<int8_t>(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
            }
            RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrDaughIndex, std::array{0}, 1);
            if (arrDaughIndex.size() == NDaughtersResonant) {
              for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
                arrPdgDaugResonant[iProng] = std::abs(daughI.pdgCode());
              }
              if ((arrPdgDaugResonant[0] == std::abs(arrPdgDaugResonantLcToPKstar0[0]) && arrPdgDaugResonant[1] == std::abs(arrPdgDaugResonantLcToPKstar0[1])) || (arrPdgDaugResonant[0] == std::abs(arrPdgDaugResonantLcToPKstar0[1]) && arrPdgDaugResonant[1] == std::abs(arrPdgDaugResonantLcToPKstar0[0]))) {
                flagChannelResonant = DecayChannelResonant::LcToPKstar0;
              } else if ((arrPdgDaugResonant[0] == std::abs(arrPdgDaugResonantLcToDeltaplusplusK[0]) && arrPdgDaugResonant[1] == std::abs(arrPdgDaugResonantLcToDeltaplusplusK[1])) || (arrPdgDaugResonant[0] == std::abs(arrPdgDaugResonantLcToDeltaplusplusK[1]) && arrPdgDaugResonant[1] == std::abs(arrPdgDaugResonantLcToDeltaplusplusK[0]))) {
                flagChannelResonant = DecayChannelResonant::LcToDeltaplusplusK;
              } else if ((arrPdgDaugResonant[0] == std::abs(arrPdgDaugResonantLcToL1520Pi[0]) && arrPdgDaugResonant[1] == std::abs(arrPdgDaugResonantLcToL1520Pi[1])) || (arrPdgDaugResonant[0] == std::abs(arrPdgDaugResonantLcToL1520Pi[1]) && arrPdgDaugResonant[1] == std::abs(arrPdgDaugResonantLcToL1520Pi[0]))) {
                flagChannelResonant = DecayChannelResonant::LcToL1520Pi;
              }
            }
          }
        }

        // Ξc± → p± K∓ π±
        if (flagChannelMain == 0) {
          auto arrPdgDaughtersXicToPKPi{std::array{+kProton, -kKPlus, +kPiPlus}};
          if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, arrPdgDaughtersXicToPKPi, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
          } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kXiCPlus, arrPdgDaughtersXicToPKPi, true, &sign, 2, &nKinkedTracks);
          } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, arrPdgDaughtersXicToPKPi, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kXiCPlus, arrPdgDaughtersXicToPKPi, true, &sign, 2);
          }
          if (indexRec > -1) {
            flagChannelMain = static_cast<int8_t>(sign * DecayChannelMain::XicToPKPi);
            flagChannelResonant = 0; // TODO
            if (arrayDaughters[0].has_mcParticle()) {
              swapping = static_cast<int8_t>(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
            }
          }
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flagChannelMain != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = static_cast<int8_t>(RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers));
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMcMatchRec(flagChannelMain, origin, swapping, flagChannelResonant, bHadMother.pt(), bHadMother.pdgCode(), nKinkedTracks, nInteractionsWithMaterial);
      } else {
        rowMcMatchRec(flagChannelMain, origin, swapping, flagChannelResonant, -1.f, 0, nKinkedTracks, nInteractionsWithMaterial);
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
        // at least one event selection not satisfied --> reject all gen particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          rowMcMatchGen(0, 0, 0, -1);
        }
        continue;
      }
      hf_mc_gen::fillMcMatchGen3Prong(mcParticles, mcParticlesPerMcColl, rowMcMatchGen, rejectBackground, matchCorrelatedBackground ? pdgMothersCorrelBkg : std::vector<int>{});
    }
  }

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisionsNoCents const& collInfos,
                 aod::McCollisions const& mcCollisions,
                 BCsInfo const& bcsInfo)
  {
    runCreator3ProngMc<CentralityEstimator::None>(tracks, mcParticles, collInfos, mcCollisions, bcsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMc, "Process MC - no centrality", false);

  void processMcCentFT0C(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Cs const& collInfos,
                         aod::McCollisions const& mcCollisions,
                         BCsInfo const& bcsInfo)
  {
    runCreator3ProngMc<CentralityEstimator::FT0C>(tracks, mcParticles, collInfos, mcCollisions, bcsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMcCentFT0C, "Process MC - FT0c centrality", false);

  void processMcCentFT0M(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Ms const& collInfos,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         BCsInfo const& bcsInfo)
  {
    runCreator3ProngMc<CentralityEstimator::FT0M>(tracks, mcParticles, collInfos, mcCollisions, bcsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMcCentFT0M, "Process MC - FT0m centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator3Prong>(cfgc, TaskName{"hf-candidate-creator-3prong"}),                         // o2-linter: disable=name/o2-task (wrong hyphenation)
    adaptAnalysisTask<HfCandidateCreator3ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-3prong-expressions"})}; // o2-linter: disable=name/o2-task (wrong hyphenation)
}
