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

#include <memory>
#include <string>
#include <utility>
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
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGHF/Utils/utilsMcGen.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreator3Prong {
  Produces<aod::HfCand3ProngBase> rowCandidateBase;
  Produces<aod::HfCand3ProngKF> rowCandidateKF;

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

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<3> df; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int runNumber{0};
  float toMicrometers = 10000.; // from cm to µm
  double massP{0.};
  double massPi{0.};
  double massK{0.};
  double massPKPi{0.};
  double massPiKP{0.};
  double massPiKPi{0.};
  double massKKPi{0.};
  double massPiKK{0.};
  double massKPi{0.};
  double massPiK{0.};
  double bz{0.};

  constexpr static float UndefValueFloat{-999.f};

  using FilteredHf3Prongs = soa::Filtered<aod::Hf3Prongs>;
  using FilteredPvRefitHf3Prongs = soa::Filtered<soa::Join<aod::Hf3Prongs, aod::HfPvRefit3Prong>>;

  // filter candidates
  Filter filterSelected3Prongs = (createDplus && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi))) != static_cast<uint8_t>(0)) || (createDs && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0)) || (createLc && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi))) != static_cast<uint8_t>(0)) || (createXic && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::XicToPKPi))) != static_cast<uint8_t>(0));

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

    std::array<bool, 4> creationFlags = {createDplus, createDs, createLc, createXic};
    if (std::accumulate(creationFlags.begin(), creationFlags.end(), 0) == 0) {
      LOGP(fatal, "At least one particle specie should be enabled for the creation.");
    }

    // histograms
    registry.add("hMass3PKPi", "3-prong candidates;inv. mass (pK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1200, 1.8, 3.0}}});
    registry.add("hMass3PiKP", "3-prong candidates;inv. mass (#pi Kp) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1200, 1.8, 3.0}}});
    registry.add("hMass3PiKPi", "3-prong candidates;inv. mass (#pi K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{600, 1.6, 2.2}}});
    registry.add("hMass3KKPi", "3-prong candidates;inv. mass (KK #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{600, 1.7, 2.3}}});
    registry.add("hMass3PiKK", "3-prong candidates;inv. mass (#pi KK) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{600, 1.7, 2.3}}});
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
    hfEvSel.addHistograms(registry); // collision monitoring

    massP = MassProton;
    massPi = MassPiPlus;
    massK = MassKPlus;

    // Configure DCAFitterN
    // df.setBz(bz);
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    runNumber = 0;

    /// candidate monitoring
    setLabelHistoCands(hCandidates);
  }

  template <bool doPvRefit = false, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll, typename Cand>
  void runCreator3ProngWithDCAFitterN(Coll const&,
                                      Cand const& rowsTrackIndexProng3,
                                      aod::TracksWCovExtra const&,
                                      aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng3.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      auto track0 = rowTrackIndexProng3.template prong0_as<aod::TracksWCovExtra>();
      auto track1 = rowTrackIndexProng3.template prong1_as<aod::TracksWCovExtra>();
      auto track2 = rowTrackIndexProng3.template prong2_as<aod::TracksWCovExtra>();
      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);

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
      std::array<float, 3> pvec0;
      std::array<float, 3> pvec1;
      std::array<float, 3> pvec2;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);
      trackParVar2.getPxPyPzGlo(pvec2);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (doPvRefit) {
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
      registry.fill(HIST("hDcaXYProngs"), track0.pt(), impactParameter0.getY() * toMicrometers);
      registry.fill(HIST("hDcaXYProngs"), track1.pt(), impactParameter1.getY() * toMicrometers);
      registry.fill(HIST("hDcaXYProngs"), track2.pt(), impactParameter2.getY() * toMicrometers);
      registry.fill(HIST("hDcaZProngs"), track0.pt(), impactParameter0.getZ() * toMicrometers);
      registry.fill(HIST("hDcaZProngs"), track1.pt(), impactParameter1.getZ() * toMicrometers);
      registry.fill(HIST("hDcaZProngs"), track2.pt(), impactParameter2.getZ() * toMicrometers);

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
      if (indexCollision == track2.collisionId() && track2.isPVContributor()) {
        SETBIT(bitmapProngsContributorsPV, 2);
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
                       pvec2[0], pvec2[1], pvec2[2],
                       impactParameter0.getY(), impactParameter1.getY(), impactParameter2.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameter2.getSigmaY2()),
                       impactParameter0.getZ(), impactParameter1.getZ(), impactParameter2.getZ(),
                       std::sqrt(impactParameter0.getSigmaZ2()), std::sqrt(impactParameter1.getSigmaZ2()), std::sqrt(impactParameter2.getSigmaZ2()),
                       rowTrackIndexProng3.prong0Id(), rowTrackIndexProng3.prong1Id(), rowTrackIndexProng3.prong2Id(), nProngsContributorsPV, bitmapProngsContributorsPV,
                       rowTrackIndexProng3.hfflag());

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        auto arrayMomenta = std::array{pvec0, pvec1, pvec2};
        massPKPi = RecoDecay::m(arrayMomenta, std::array{massP, massK, massPi});
        massPiKP = RecoDecay::m(arrayMomenta, std::array{massPi, massK, massP});
        massPiKPi = RecoDecay::m(arrayMomenta, std::array{massPi, massK, massPi});
        massKKPi = RecoDecay::m(arrayMomenta, std::array{massK, massK, massPi});
        massPiKK = RecoDecay::m(arrayMomenta, std::array{massPi, massK, massK});
        massKPi = RecoDecay::m(std::array{arrayMomenta.at(1), arrayMomenta.at(2)}, std::array{massK, massPi});
        massPiK = RecoDecay::m(std::array{arrayMomenta.at(0), arrayMomenta.at(1)}, std::array{massPi, massK});
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

  template <bool doPvRefit = false, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll, typename Cand>
  void runCreator3ProngWithKFParticle(Coll const&,
                                      Cand const& rowsTrackIndexProng3,
                                      aod::TracksWCovExtra const&,
                                      aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {
      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng3.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      auto track0 = rowTrackIndexProng3.template prong0_as<aod::TracksWCovExtra>();
      auto track1 = rowTrackIndexProng3.template prong1_as<aod::TracksWCovExtra>();
      auto track2 = rowTrackIndexProng3.template prong2_as<aod::TracksWCovExtra>();

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

      KFPTrack kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(track1);
      KFPTrack kfpTrack2 = createKFPTrackFromTrack(track2);

      KFParticle kfFirstProton(kfpTrack0, kProton);
      KFParticle kfFirstPion(kfpTrack0, kPiPlus);
      KFParticle kfFirstKaon(kfpTrack0, kKPlus);
      KFParticle kfSecondKaon(kfpTrack1, kKPlus);
      KFParticle kfThirdProton(kfpTrack2, kProton);
      KFParticle kfThirdPion(kfpTrack2, kPiPlus);
      KFParticle kfThirdKaon(kfpTrack2, kKPlus);

      float impactParameter0XY = 0., errImpactParameter0XY = 0., impactParameter1XY = 0., errImpactParameter1XY = 0., impactParameter2XY = 0., errImpactParameter2XY = 0.;
      if (!kfFirstProton.GetDistanceFromVertexXY(kfpV, impactParameter0XY, errImpactParameter0XY)) {
        registry.fill(HIST("hDcaXYProngs"), track0.pt(), impactParameter0XY * toMicrometers);
        registry.fill(HIST("hDcaZProngs"), track0.pt(), std::sqrt(kfFirstProton.GetDistanceFromVertex(kfpV) * kfFirstProton.GetDistanceFromVertex(kfpV) - impactParameter0XY * impactParameter0XY) * toMicrometers);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track0.pt(), UndefValueFloat);
        registry.fill(HIST("hDcaZProngs"), track0.pt(), UndefValueFloat);
      }
      if (!kfSecondKaon.GetDistanceFromVertexXY(kfpV, impactParameter1XY, errImpactParameter1XY)) {
        registry.fill(HIST("hDcaXYProngs"), track1.pt(), impactParameter1XY * toMicrometers);
        registry.fill(HIST("hDcaZProngs"), track1.pt(), std::sqrt(kfSecondKaon.GetDistanceFromVertex(kfpV) * kfSecondKaon.GetDistanceFromVertex(kfpV) - impactParameter1XY * impactParameter1XY) * toMicrometers);
      } else {
        registry.fill(HIST("hDcaXYProngs"), track1.pt(), UndefValueFloat);
        registry.fill(HIST("hDcaZProngs"), track1.pt(), UndefValueFloat);
      }
      if (!kfThirdProton.GetDistanceFromVertexXY(kfpV, impactParameter2XY, errImpactParameter2XY)) {
        registry.fill(HIST("hDcaXYProngs"), track2.pt(), impactParameter2XY * toMicrometers);
        registry.fill(HIST("hDcaZProngs"), track2.pt(), std::sqrt(kfThirdProton.GetDistanceFromVertex(kfpV) * kfThirdProton.GetDistanceFromVertex(kfpV) - impactParameter2XY * impactParameter2XY) * toMicrometers);
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

      const float chi2geo = kfCandPKPi.Chi2() / kfCandPKPi.NDF();
      const float chi2topo = kfCalculateChi2ToPrimaryVertex(kfCandPKPi, kfpV);
      const std::pair<float, float> ldl = kfCalculateLdL(kfCandPKPi, kfpV);

      std::array<float, 3> pProng0 = kfCalculateProngMomentumInSecondaryVertex(kfFirstProton, kfCandPiKP);
      std::array<float, 3> pProng1 = kfCalculateProngMomentumInSecondaryVertex(kfSecondKaon, kfCandPiKP);
      std::array<float, 3> pProng2 = kfCalculateProngMomentumInSecondaryVertex(kfThirdPion, kfCandPiKP);

      registry.fill(HIST("hCovSVXX"), kfCandPKPi.Covariance(0, 0));
      registry.fill(HIST("hCovSVYY"), kfCandPKPi.Covariance(1, 1));
      registry.fill(HIST("hCovSVXZ"), kfCandPKPi.Covariance(2, 0));
      registry.fill(HIST("hCovSVZZ"), kfCandPKPi.Covariance(2, 2));

      auto covMatrixSV = kfCandPKPi.CovarianceMatrix();
      double phi, theta;
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
      uint8_t nProngsContributorsPV = hf_trkcandsel::countOnesInBinary(bitmapProngsContributorsPV);

      // fill candidate table rows
      rowCandidateBase(indexCollision,
                       kfpV.GetX(), kfpV.GetY(), kfpV.GetZ(),
                       kfCandPKPi.GetX(), kfCandPKPi.GetY(), kfCandPKPi.GetZ(),
                       errorDecayLength, errorDecayLengthXY,
                       kfCandPKPi.GetChi2() / kfCandPKPi.GetNDF(),
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
                                    aod::TracksWCovExtra const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter with PV refit and w/o centrality selections", false);

  /// @brief process function using DCA fitter  w/o PV refit and w/o centrality selections
  void processNoPvRefitWithDCAFitterN(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      FilteredHf3Prongs const& rowsTrackIndexProng3,
                                      aod::TracksWCovExtra const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterN, "Run candidate creator using DCA fitter without PV refit and w/o centrality selections", true);

  /// @brief process function using KFParticle package  w/ PV refit and w/o centrality selections
  void processPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                    FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                    aod::TracksWCovExtra const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticle, "Run candidate creator using KFParticle package with PV refit and w/o centrality selections", false);

  /// @brief process function using KFParticle package  w/o PV refit and w/o centrality selections
  void processNoPvRefitWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                      FilteredHf3Prongs const& rowsTrackIndexProng3,
                                      aod::TracksWCovExtra const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
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
                                            aod::TracksWCovExtra const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter with PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using DCA fitter  w/o PV refit and  w/ centrality selection on FT0C
  void processNoPvRefitWithDCAFitterNCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              aod::TracksWCovExtra const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterNCentFT0C, "Run candidate creator using DCA fitter without PV refit and  w/ centrality selection on FT0C", false);

  /// @brief process function using KFParticle package  w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                            aod::TracksWCovExtra const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticleCentFT0C, "Run candidate creator using KFParticle package with PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function using KFParticle package  w/o PV refit and  w/ centrality selection on FT0C
  void processNoPvRefitWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              aod::TracksWCovExtra const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
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
                                            aod::TracksWCovExtra const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ true, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter with PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using DCA fitter  w/o PV refit and  w/ centrality selection on FT0M
  void processNoPvRefitWithDCAFitterNCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              aod::TracksWCovExtra const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithDCAFitterN</*doPvRefit*/ false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithDCAFitterNCentFT0M, "Run candidate creator using DCA fitter without PV refit and  w/ centrality selection on FT0M", false);

  /// @brief process function using KFParticle package  w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                                            aod::TracksWCovExtra const& tracks,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ true, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package with PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function using KFParticle package  w/o PV refit and  w/ centrality selection on FT0M
  void processNoPvRefitWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              FilteredHf3Prongs const& rowsTrackIndexProng3,
                                              aod::TracksWCovExtra const& tracks,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3ProngWithKFParticle</*doPvRefit*/ false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitWithKFParticleCentFT0M, "Run candidate creator using KFParticle package without PV refit and  w/ centrality selection on FT0M", false);

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
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisions, "Collision monitoring - no centrality", true);

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
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

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
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
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
      if (device.name.compare("hf-candidate-creator-3prong") == 0) {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }

    hfEvSelMc.addHistograms(registry); // particles monitoring
  }

  /// Performs MC matching.
  template <o2::hf_centrality::CentralityEstimator centEstimator, typename CCs, typename McCollisions>
  void runCreator3ProngMc(aod::TracksWMc const& tracks,
                          aod::McParticles const& mcParticles,
                          CCs const& collInfos,
                          McCollisions const& mcCollisions,
                          BCsInfo const&)
  {
    rowCandidateProng3->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t swapping = 0;
    int8_t channel = 0;
    int8_t nKinkedTracks = 0;
    int8_t nInteractionsWithMaterial = 0;
    std::vector<int> arrDaughIndex;
    std::array<int, 2> arrPDGDaugh;
    std::array<int, 2> arrPDGResonant1 = {kProton, 313};      // Λc± → p± K*
    std::array<int, 2> arrPDGResonant2 = {2224, kKPlus};      // Λc± → Δ(1232)±± K∓
    std::array<int, 2> arrPDGResonant3 = {102134, kPiPlus};   // Λc± → Λ(1520) π±
    std::array<int, 2> arrPDGResonantDPhiPi = {333, kPiPlus}; // Ds± → Phi π± and D± → Phi π±
    std::array<int, 2> arrPDGResonantDKstarK = {313, kKPlus}; // Ds± → K*(892)0bar K± and D± → K*(892)0bar K±

    // Match reconstructed candidates.
    // Spawned table can be used directly
    for (const auto& candidate : *rowCandidateProng3) {
      flag = 0;
      origin = 0;
      swapping = 0;
      channel = 0;
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
          rowMcMatchRec(flag, origin, swapping, channel, -1.f, 0, 0, 0);
          continue;
        }
      }

      // D± → π± K∓ π±
      if (flag == 0) {
        if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
        } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks);
        } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
        }
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::DplusToPiKPi);
        }
      }

      // Ds± → K± K∓ π± and D± → K± K∓ π±
      if (flag == 0) {
        bool isDplus = false;
        if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
        } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks);
        } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
        }
        if (indexRec == -1) {
          isDplus = true;
          if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
          } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks);
          } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
          } else {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
          }
        }
        if (indexRec > -1) {
          // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
          // TODO: move to different and explicit flags
          flag = sign * (1 << DecayType::DsToKKPi);
          if (arrayDaughters[0].has_mcParticle()) {
            swapping = int8_t(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
          }
          RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
              arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonantDPhiPi[0] && arrPDGDaugh[1] == arrPDGResonantDPhiPi[1]) || (arrPDGDaugh[0] == arrPDGResonantDPhiPi[1] && arrPDGDaugh[1] == arrPDGResonantDPhiPi[0])) {
              channel = isDplus ? DecayChannelDToKKPi::DplusToPhiPi : DecayChannelDToKKPi::DsToPhiPi;
            } else if ((arrPDGDaugh[0] == arrPDGResonantDKstarK[0] && arrPDGDaugh[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaugh[0] == arrPDGResonantDKstarK[1] && arrPDGDaugh[1] == arrPDGResonantDKstarK[0])) {
              channel = isDplus ? DecayChannelDToKKPi::DplusToK0starK : DecayChannelDToKKPi::DsToK0starK;
            }
          }
        }
      }

      // D* → D0π → Kππ
      if (flag == 0) {
        if (matchKinkedDecayTopology) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true>(mcParticles, arrayDaughters, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &sign, 2, &nKinkedTracks);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &sign, 2);
        }
        if (indexRec > -1) {
          flag = sign * (1 << DstarToPiKPiBkg);
          channel = 1;
        }
      }

      // Λc± → p± K∓ π±
      if (flag == 0) {
        if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
        } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks);
        } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        }
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::LcToPKPi);

          // Flagging the different Λc± → p± K∓ π± decay channels
          if (arrayDaughters[0].has_mcParticle()) {
            swapping = int8_t(std::abs(arrayDaughters[0].mcParticle().pdgCode()) == kPiPlus);
          }
          RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRec), &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
              arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonant1[0] && arrPDGDaugh[1] == arrPDGResonant1[1]) || (arrPDGDaugh[0] == arrPDGResonant1[1] && arrPDGDaugh[1] == arrPDGResonant1[0])) {
              channel = 1;
            } else if ((arrPDGDaugh[0] == arrPDGResonant2[0] && arrPDGDaugh[1] == arrPDGResonant2[1]) || (arrPDGDaugh[0] == arrPDGResonant2[1] && arrPDGDaugh[1] == arrPDGResonant2[0])) {
              channel = 2;
            } else if ((arrPDGDaugh[0] == arrPDGResonant3[0] && arrPDGDaugh[1] == arrPDGResonant3[1]) || (arrPDGDaugh[0] == arrPDGResonant3[1] && arrPDGDaugh[1] == arrPDGResonant3[0])) {
              channel = 3;
            }
          }
        }
      }

      // Ξc± → p± K∓ π±
      if (flag == 0) {
        if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks, &nInteractionsWithMaterial);
        } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, &nKinkedTracks);
        } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, nullptr, &nInteractionsWithMaterial);
        } else {
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        }
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::XicToPKPi);
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMcMatchRec(flag, origin, swapping, channel, bHadMother.pt(), bHadMother.pdgCode(), nKinkedTracks, nInteractionsWithMaterial);
      } else {
        rowMcMatchRec(flag, origin, swapping, channel, -1.f, 0, nKinkedTracks, nInteractionsWithMaterial);
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
        // at least one event selection not satisfied --> reject all gen particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          rowMcMatchGen(0, 0, 0, -1);
        }
        continue;
      }
      hf_mc_gen::fillMcMatchGen3Prong(mcParticles, mcParticlesPerMcColl, rowMcMatchGen, rejectBackground);
    }
  }

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisionsNoCents const& collInfos,
                 aod::McCollisions const& mcCollisions,
                 BCsInfo const& BCsInfo)
  {
    runCreator3ProngMc<CentralityEstimator::None>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMc, "Process MC - no centrality", false);

  void processMcCentFT0C(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Cs const& collInfos,
                         aod::McCollisions const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreator3ProngMc<CentralityEstimator::FT0C>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMcCentFT0C, "Process MC - FT0c centrality", false);

  void processMcCentFT0M(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Ms const& collInfos,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreator3ProngMc<CentralityEstimator::FT0M>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMcCentFT0M, "Process MC - FT0m centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator3Prong>(cfgc, TaskName{"hf-candidate-creator-3prong"}),                         // o2-linter: disable=name/o2-task (wrong hyphenation)
    adaptAnalysisTask<HfCandidateCreator3ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-3prong-expressions"})}; // o2-linter: disable=name/o2-task (wrong hyphenation)
}
