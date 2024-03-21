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

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_collision_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreator3Prong {
  Produces<aod::HfCand3ProngBase> rowCandidateBase;

  // centrality
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality"};
  // event selection
  Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "apply the sel8 event selection"};
  Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "max. PV posZ (cm)"};
  Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "apply TF border cut"};
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
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  // flags to enable creation for different particle species separately
  Configurable<bool> createDplus{"createDplus", false, "enable D+/- candidate creation"};
  Configurable<bool> createDs{"createDs", false, "enable Ds+/- candidate creation"};
  Configurable<bool> createLc{"createLc", false, "enable Lc+/- candidate creation"};
  Configurable<bool> createXic{"createXic", false, "enable Xic+/- candidate creation"};

  o2::vertexing::DCAFitterN<3> df; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber{0};
  float toMicrometers = 10000.; // from cm to µm
  double massPi{0.};
  double massK{0.};
  double massPiKPi{0.};
  double bz{0.};

  using FilteredHf3Prongs = soa::Filtered<aod::Hf3Prongs>;
  using FilteredPvRefitHf3Prongs = soa::Filtered<soa::Join<aod::Hf3Prongs, aod::HfPvRefit3Prong>>;

  // filter candidates
  Filter filterSelected3Prongs = (createDplus && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi))) != static_cast<uint8_t>(0)) || (createDs && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0)) || (createLc && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::LcToPKPi))) != static_cast<uint8_t>(0)) || (createXic && (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::XicToPKPi))) != static_cast<uint8_t>(0));

  std::shared_ptr<TH1> hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 6> processes = {doprocessPvRefit, doprocessNoPvRefit,
                                     doprocessPvRefitCentFT0C, doprocessNoPvRefitCentFT0C,
                                     doprocessPvRefitCentFT0M, doprocessNoPvRefitCentFT0M};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    std::array<bool, 3> processesCollisions = {doprocessCollisions, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};
    const int nProcessesCollisions = std::accumulate(processesCollisions.begin(), processesCollisions.end(), 0);
    if (nProcessesCollisions > 1) {
      LOGP(fatal, "At most one process function for collision monitoring can be enabled at a time.");
    }
    if (nProcessesCollisions == 1) {
      if ((doprocessPvRefit || doprocessNoPvRefit) && !doprocessCollisions) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisions\"?");
      }
      if ((doprocessPvRefitCentFT0C || doprocessNoPvRefitCentFT0C) && !doprocessCollisionsCentFT0C) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0C\"?");
      }
      if ((doprocessPvRefitCentFT0M || doprocessNoPvRefitCentFT0M) && !doprocessCollisionsCentFT0M) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0M\"?");
      }
    }

    std::array<bool, 4> creationFlags = {createDplus, createDs, createLc, createXic};
    if (std::accumulate(creationFlags.begin(), creationFlags.end(), 0) == 0) {
      LOGP(fatal, "At least one particle specie should be enabled for the creation.");
    }

    // histograms
    registry.add("hMass3", "3-prong candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.6, 2.1}}});
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
    hCollisions = registry.add<TH1>("hCollisions", "HF event counter;;entries", {HistType::kTH1D, {axisEvents}});
    hPosZBeforeEvSel = registry.add<TH1>("hPosZBeforeEvSel", "all events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{400, -20., 20.}}});
    hPosZAfterEvSel = registry.add<TH1>("hPosZAfterEvSel", "selected events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{400, -20., 20.}}});
    hPosXAfterEvSel = registry.add<TH1>("hPosXAfterEvSel", "selected events;#it{x}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{200, -0.5, 0.5}}});
    hPosYAfterEvSel = registry.add<TH1>("hPosYAfterEvSel", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{200, -0.5, 0.5}}});
    hNumPvContributorsAfterSel = registry.add<TH1>("hNumPvContributorsAfterSel", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{500, -0.5, 499.5}}});

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
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    /// collision monitoring
    setLabelHistoEvSel(hCollisions);
  }

  template <bool doPvRefit = false, o2::aod::hf_collision_centrality::CentralityEstimator centEstimator, typename Coll, typename Cand>
  void runCreator3Prong(Coll const& collisions,
                        Cand const& rowsTrackIndexProng3,
                        aod::TracksWCovExtra const& tracks,
                        aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexProng3.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = getHfCollisionRejectionMask<true, centEstimator>(collision, centrality, centralityMin, centralityMax, useSel8Trigger, -1, useTimeFrameBorderCut, -zPvPosMax, zPvPosMax, 0, -1.f);
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
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      df.setBz(bz);

      // reconstruct the 3-prong secondary vertex
      if (df.process(trackParVar0, trackParVar1, trackParVar2) == 0) {
        continue;
      }
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
      uint8_t nProngsContributorsPV = 0;
      if (indexCollision == track0.collisionId() && track0.isPVContributor()) {
        nProngsContributorsPV += 1;
      }
      if (indexCollision == track1.collisionId() && track1.isPVContributor()) {
        nProngsContributorsPV += 1;
      }
      if (indexCollision == track2.collisionId() && track2.isPVContributor()) {
        nProngsContributorsPV += 1;
      }

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
                       rowTrackIndexProng3.prong0Id(), rowTrackIndexProng3.prong1Id(), rowTrackIndexProng3.prong2Id(), nProngsContributorsPV,
                       rowTrackIndexProng3.hfflag());

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        auto arrayMomenta = std::array{pvec0, pvec1, pvec2};
        massPiKPi = RecoDecay::m(std::move(arrayMomenta), std::array{massPi, massK, massPi});
        registry.fill(HIST("hMass3"), massPiKPi);
      }
    }
  }

  ///////////////////////////////////
  ///                             ///
  ///   No centrality selection   ///
  ///                             ///
  ///////////////////////////////////

  /// @brief process function w/ PV refit and w/o centrality selections
  void processPvRefit(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                      FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                      aod::TracksWCovExtra const& tracks,
                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3Prong</*doPvRefit*/ true, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefit, "Run candidate creator with PV refit and w/o centrality selections", false);

  /// @brief process function w/o PV refit and w/o centrality selections
  void processNoPvRefit(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        FilteredHf3Prongs const& rowsTrackIndexProng3,
                        aod::TracksWCovExtra const& tracks,
                        aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3Prong</*doPvRefit*/ false, CentralityEstimator::None>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefit, "Run candidate creator without PV refit and w/o centrality selections", true);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0C   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                              FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                              aod::TracksWCovExtra const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3Prong</*doPvRefit*/ true, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitCentFT0C, "Run candidate creator with PV refit and w/ centrality selection on FT0C", false);

  /// @brief process function w/o PV refit and  w/ centrality selection on FT0C
  void processNoPvRefitCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                FilteredHf3Prongs const& rowsTrackIndexProng3,
                                aod::TracksWCovExtra const& tracks,
                                aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3Prong</*doPvRefit*/ false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitCentFT0C, "Run candidate creator without PV refit and  w/ centrality selection on FT0C", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0M   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                              FilteredPvRefitHf3Prongs const& rowsTrackIndexProng3,
                              aod::TracksWCovExtra const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3Prong</*doPvRefit*/ true, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processPvRefitCentFT0M, "Run candidate creator with PV refit and w/ centrality selection on FT0M", false);

  /// @brief process function w/o PV refit and  w/ centrality selection on FT0M
  void processNoPvRefitCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                FilteredHf3Prongs const& rowsTrackIndexProng3,
                                aod::TracksWCovExtra const& tracks,
                                aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreator3Prong</*doPvRefit*/ false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexProng3, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreator3Prong, processNoPvRefitCentFT0M, "Run candidate creator without PV refit and  w/ centrality selection on FT0M", false);

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
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisions, "Collision monitoring - no centrality", true);

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
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

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
  PROCESS_SWITCH(HfCandidateCreator3Prong, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
};

/// Extends the base table with expression columns.
struct HfCandidateCreator3ProngExpressions {
  Spawns<aod::HfCand3ProngExt> rowCandidateProng3;
  Produces<aod::HfCand3ProngMcRec> rowMcMatchRec;
  Produces<aod::HfCand3ProngMcGen> rowMcMatchGen;

  bool createDplus{false};
  bool createDs{false};
  bool createLc{false};
  bool createXic{false};

  void init(InitContext& initContext)
  {

    // inspect for which particle species the candidates were created
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-3prong") == 0) {
        for (const auto& option : device.options) {
          if (option.name.compare("createDplus") == 0) {
            createDplus = option.defaultValue.get<bool>();
          } else if (option.name.compare("createDs") == 0) {
            createDs = option.defaultValue.get<bool>();
          } else if (option.name.compare("createLc") == 0) {
            createLc = option.defaultValue.get<bool>();
          } else if (option.name.compare("createXic") == 0) {
            createXic = option.defaultValue.get<bool>();
          }
        }
      }
    }

    LOGP(info, "Flags for candidate creation from the reco workflow:");
    LOGP(info, "    --> createDplus = {}", createDplus);
    LOGP(info, "    --> createDs = {}", createDs);
    LOGP(info, "    --> createLc = {}", createLc);
    LOGP(info, "    --> createXic = {}", createXic);
  }

  /// Performs MC matching.
  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    rowCandidateProng3->bindExternalIndices(&tracks);

    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t swapping = 0;
    int8_t channel = 0;
    std::vector<int> arrDaughIndex;
    std::array<int, 2> arrPDGDaugh;
    std::array<int, 2> arrPDGResonant1 = {kProton, 313};      // Λc± → p± K*
    std::array<int, 2> arrPDGResonant2 = {2224, kKPlus};      // Λc± → Δ(1232)±± K∓
    std::array<int, 2> arrPDGResonant3 = {3124, kPiPlus};     // Λc± → Λ(1520) π±
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
      auto arrayDaughters = std::array{candidate.prong0_as<aod::TracksWMc>(), candidate.prong1_as<aod::TracksWMc>(), candidate.prong2_as<aod::TracksWMc>()};

      // D± → π± K∓ π±
      if (createDplus) {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::DplusToPiKPi);
        }
      }

      // Ds± → K± K∓ π± and D± → K± K∓ π±
      if (flag == 0 && createDs) {
        bool isDplus = false;
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec == -1) {
          isDplus = true;
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
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

      // Λc± → p± K∓ π±
      if (flag == 0 && createLc) {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
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
      if (flag == 0 && createXic) {
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec > -1) {
          flag = sign * (1 << DecayType::XicToPKPi);
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }

      rowMcMatchRec(flag, origin, swapping, channel);
    }

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = 0;
      origin = 0;
      channel = 0;
      arrDaughIndex.clear();

      // D± → π± K∓ π±
      if (createDplus) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flag = sign * (1 << DecayType::DplusToPiKPi);
        }
      }

      // Ds± → K± K∓ π± and D± → K± K∓ π±
      if (flag == 0 && createDs) {
        bool isDplus = false;
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
          // TODO: move to different and explicit flags
          flag = sign * (1 << DecayType::DsToKKPi);
        } else if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
          // TODO: move to different and explicit flags
          flag = sign * (1 << DecayType::DsToKKPi);
          isDplus = true;
        }
        if (flag != 0) {
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) {
              auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
              arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
            }
            if ((arrPDGDaugh[0] == arrPDGResonantDPhiPi[0] && arrPDGDaugh[1] == arrPDGResonantDPhiPi[1]) || (arrPDGDaugh[0] == arrPDGResonantDPhiPi[1] && arrPDGDaugh[1] == arrPDGResonantDPhiPi[0])) {
              channel = isDplus ? DecayChannelDToKKPi::DplusToPhiPi : DecayChannelDToKKPi::DsToPhiPi;
            } else if ((arrPDGDaugh[0] == arrPDGResonantDKstarK[0] && arrPDGDaugh[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaugh[0] == arrPDGResonantDKstarK[1] && arrPDGDaugh[1] == arrPDGResonantDKstarK[0])) {
              channel = isDplus ? DecayChannelDToKKPi::DplusToK0starK : DecayChannelDToKKPi::DsToK0starK;
            }
          }
        }
      }

      // Λc± → p± K∓ π±
      if (flag == 0 && createLc) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flag = sign * (1 << DecayType::LcToPKPi);

          // Flagging the different Λc± → p± K∓ π± decay channels
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == 2) {
            for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) {
              auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
              arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
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
      if (flag == 0 && createXic) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flag = sign * (1 << DecayType::XicToPKPi);
        }
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }

      rowMcMatchGen(flag, origin, channel);
    }
  }

  PROCESS_SWITCH(HfCandidateCreator3ProngExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreator3Prong>(cfgc, TaskName{"hf-candidate-creator-3prong"}),
    adaptAnalysisTask<HfCandidateCreator3ProngExpressions>(cfgc, TaskName{"hf-candidate-creator-3prong-expressions"})};
}
