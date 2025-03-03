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

/// \file candidateCreatorDstar.cxx
/// \brief Reconstruction of D* decay candidates
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
// std
#include <memory>
#include <string>
#include <vector>
// ROOT
#include <TPDGCode.h>
// O2
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
// O2Physics
#include "Common/Core/trackUtilities.h"
// PWGHF
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

using namespace o2;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::hf_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;

namespace o2::aod
{
using HfDstarsWithPvRefitInfo = soa::Join<aod::HfDstars, aod::HfPvRefitDstar>;
} // namespace o2::aod

/// Reconstruction of D* decay candidates
struct HfCandidateCreatorDstar {
  Produces<aod::HfD0FromDstarBase> rowCandD0Base;
  Produces<aod::HfCandDstarBase> rowCandDstarBase;

  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};

  // magnetic field setting from CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};                                   // ........... what is unit of this?
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"}; // ..........What it DZ?
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};

  HfEventSelection hfEvSel;                 // event selection and monitoring
  Service<o2::ccdb::BasicCCDBManager> ccdb; // From utilsBfieldCCDB.h
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  // D0-prong vertex fitter
  o2::vertexing::DCAFitterN<2> df;
  int runNumber;
  double bz;
  static constexpr float CmToMicrometers = 10000.; // from cm to µm
  double massPi, massK, massD0;

  AxisSpec ptAxis = {100, 0., 2.0, "#it{p}_{T} (GeV/#it{c}"};
  AxisSpec dcaAxis = {200, -500., 500., "#it{d}_{xy,z} (#mum)"};

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{
    "registry",
    {{"Refit/hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}}},
     {"Refit/hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"Refit/hCovPVYY", "2-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}}},
     {"Refit/hCovSVYY", "2-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"Refit/hCovPVXZ", "2-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, -1.e-4, 1.e-4}}}},
     {"Refit/hCovSVXZ", "2-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, -1.e-4, 0.2}}}},
     {"Refit/hCovPVZZ", "2-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}}},
     {"Refit/hCovSVZZ", "2-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}}},

     {"QA/hDcaXYProngsD0", "DCAxy of 2-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", {HistType::kTH2F, {ptAxis, dcaAxis}}},
     {"QA/hDcaZProngsD0", "DCAz of 2-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", {HistType::kTH2F, {ptAxis, dcaAxis}}},
     {"QA/hDCAXYPi", "DCAxy of Soft Pi;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", {HistType::kTH2F, {ptAxis, dcaAxis}}},
     {"QA/hDCAZPi", "DCAz of Soft Pi;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", {HistType::kTH2F, {ptAxis, dcaAxis}}},

     {"QA/hPtPi", "#pi candidates", {HistType::kTH1F, {ptAxis}}},
     {"QA/hPtD0Prong0", "D^{0} candidates' prong0", {HistType::kTH1F, {ptAxis}}},
     {"QA/hPtD0Prong1", "D^{0} candidates' prong1", {HistType::kTH1F, {ptAxis}}},
     {"QA/hPtD0", "D^{0} candidates", {HistType::kTH1F, {ptAxis}}},
     {"QA/hPtDstar", "D* candidates", {HistType::kTH1F, {ptAxis}}}}};

  /// @brief This function initializes the ccdb setting, vertex fitter and runs function MatLayerCylSet::rectifyPtrFromFile(..args..)
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

    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    hfEvSel.addHistograms(registry); // collision monitoring

    // LOG(info) << "Init Function Invoked";
    massPi = MassPiPlus;
    massK = MassKPlus;
    massD0 = MassD0;

    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);
    df.setMatCorrType(noMatCorr);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking(); // set the flag to check object validity before CCDB query
    runNumber = 0;
    bz = 0;

    /// candidate monitoring
    setLabelHistoCands(hCandidates);
  }

  /// @brief function for secondary vertex reconstruction and candidate creator
  /// @tparam CandsDstar Table type of Dstar table object
  /// @tparam doPvRefit True/False PV refit option
  /// @param collisions collision object
  /// @param rowsTrackIndexDstar Dstar table object from trackIndexSkimCreator.cxx
  /// @param rowsTrackIndexD0 D0 table object from trackIndexSkimCreator.cxx
  /// @param tracks track table with Cov object
  /// @param bcWithTimeStamps Bunch Crossing with timestamps
  template <bool doPvRefit, o2::hf_centrality::CentralityEstimator centEstimator, typename Coll, typename CandsDstar>
  void runCreatorDstar(Coll const&,
                       CandsDstar const& rowsTrackIndexDstar,
                       aod::Hf2Prongs const&,
                       aod::TracksWCov const&,
                       aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    // LOG(info) << "runCreatorDstar function called";
    // LOG(info) << "candidate loop starts";
    // loop over suspected Dstar Candidate
    for (const auto& rowTrackIndexDstar : rowsTrackIndexDstar) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = rowTrackIndexDstar.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      auto trackPi = rowTrackIndexDstar.template prong0_as<aod::TracksWCov>();
      auto prongD0 = rowTrackIndexDstar.template prongD0_as<aod::Hf2Prongs>();
      auto trackD0Prong0 = prongD0.template prong0_as<aod::TracksWCov>();
      auto trackD0Prong1 = prongD0.template prong1_as<aod::TracksWCov>();

      // Extracts primary vertex position and covariance matrix from a collision
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

      // Extracts track parameters and covariance matrix from a track
      auto trackPiParVar = getTrackParCov(trackPi);
      // These will be used in DCA Fitter to reconstruct secondary vertex
      auto trackD0Prong0ParVarPos1 = getTrackParCov(trackD0Prong0); // from trackUtilities.h
      auto trackD0Prong1ParVarNeg1 = getTrackParCov(trackD0Prong1);

      // auto collisionPiId = trackPi.collisionId();
      // auto collisionD0Id = trackD0Prong0.collisionId();
      // LOGF(info, "Pi collision %ld, D0 collision %ld", collisionPiId, collisionD0Id);
      //..................................................Doubt: Should I apply a condition of (collisionPiId == collisionD0Id)

      /// Set the magnetic field from ccdb.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        // LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        // LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        runNumber = bc.runNumber();
      }
      df.setBz(bz);

      // reconstruct the 2-prong secondary vertex
      hCandidates->Fill(SVFitting::BeforeFit);
      try {
        if (df.process(trackD0Prong0ParVarPos1, trackD0Prong1ParVarNeg1) == 0) {
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

      registry.fill(HIST("Refit/hCovSVXX"), covMatrixPCA[0]);
      registry.fill(HIST("Refit/hCovSVYY"), covMatrixPCA[2]);
      registry.fill(HIST("Refit/hCovSVXZ"), covMatrixPCA[3]);
      registry.fill(HIST("Refit/hCovSVZZ"), covMatrixPCA[5]);

      // Doubt:................Below, track object are at secondary vertex!
      // < track param propagated to V0 candidate (no check for the candidate validity). propagateTracksToVertex must be called in advance
      auto trackD0ProngParVar0 = df.getTrack(0);
      auto trackD0ProngParVar1 = df.getTrack(1);

      std::array<float, 3> pVecD0Prong0;
      std::array<float, 3> pVecD0Prong1;
      trackD0ProngParVar0.getPxPyPzGlo(pVecD0Prong0);
      trackD0ProngParVar1.getPxPyPzGlo(pVecD0Prong1);

      // This modifies track momenta!
      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the *HfCand3ProngBase/HfCand2ProngBase* all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexDstar.pvRefitX());
        primaryVertex.setY(rowTrackIndexDstar.pvRefitY());
        primaryVertex.setZ(rowTrackIndexDstar.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexDstar.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexDstar.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexDstar.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexDstar.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexDstar.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexDstar.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov(); /// Here covMatrixPV Updated!
      }
      registry.fill(HIST("Refit/hCovPVXX"), covMatrixPV[0]);
      registry.fill(HIST("Refit/hCovPVYY"), covMatrixPV[2]);
      registry.fill(HIST("Refit/hCovPVXZ"), covMatrixPV[3]);
      registry.fill(HIST("Refit/hCovPVZZ"), covMatrixPV[5]);

      // get track impact parameters
      o2::dataformats::DCA impactParameter0; // GPUROOTCartesianFwd.h
      o2::dataformats::DCA impactParameter1;
      // Propagating D0 prongs to DCA
      trackD0ProngParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackD0ProngParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);

      // Propagating Soft Pi to DCA
      o2::dataformats::DCA impactParameterPi;
      trackPiParVar.propagateToDCA(primaryVertex, bz, &impactParameterPi);
      registry.fill(HIST("QA/hDcaXYProngsD0"), trackD0Prong0.pt(), impactParameter0.getY() * CmToMicrometers);
      registry.fill(HIST("QA/hDcaXYProngsD0"), trackD0Prong1.pt(), impactParameter1.getY() * CmToMicrometers);
      registry.fill(HIST("QA/hDcaZProngsD0"), trackD0Prong0.pt(), impactParameter0.getZ() * CmToMicrometers);
      registry.fill(HIST("QA/hDcaZProngsD0"), trackD0Prong1.pt(), impactParameter1.getZ() * CmToMicrometers);

      registry.fill(HIST("QA/hDCAXYPi"), trackPi.pt(), impactParameterPi.getY() * CmToMicrometers);
      registry.fill(HIST("QA/hDCAZPi"), trackPi.pt(), impactParameterPi.getZ() * CmToMicrometers);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      // Calculates the XX element of a XYZ covariance matrix after rotation of the coordinate system by phi around the z-axis and by minus theta around the new y-axis.
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // Calculation of kinematics for inv mass
      auto pVecD0 = RecoDecay::pVec(pVecD0Prong0, pVecD0Prong1);

      // D0 pt magnitude
      auto ptD0 = RecoDecay::pt(pVecD0);

      // Soft pi momentum vector and sign
      std::array<float, 3> pVecSoftPi;
      trackPiParVar.getPxPyPzGlo(pVecSoftPi);
      int8_t signSoftPi = static_cast<int8_t>(trackPi.sign());

      // D* pt magnitude
      auto ptDstar = RecoDecay::pt(pVecD0, pVecSoftPi);

      // Fill candidate Table for DStar
      rowCandDstarBase(collision.globalIndex(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       rowTrackIndexDstar.prong0Id(), rowTrackIndexDstar.prongD0Id(),
                       pVecSoftPi[0], pVecSoftPi[1], pVecSoftPi[2],
                       signSoftPi,
                       impactParameterPi.getY(), impactParameterPi.getZ(),
                       std::sqrt(impactParameterPi.getSigmaY2()), std::sqrt(impactParameterPi.getSigmaZ2()),
                       pVecD0Prong0[0], pVecD0Prong0[1], pVecD0Prong0[2],
                       pVecD0Prong1[0], pVecD0Prong1[1], pVecD0Prong1[2],
                       prongD0.prong0Id(), prongD0.prong1Id());
      // Fill candidate Table for D0
      rowCandD0Base(collision.globalIndex(),
                    primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                    secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                    errorDecayLength, errorDecayLengthXY,
                    chi2PCA,
                    pVecD0Prong0[0], pVecD0Prong0[1], pVecD0Prong0[2],
                    pVecD0Prong1[0], pVecD0Prong1[1], pVecD0Prong1[2],
                    impactParameter0.getY(), impactParameter1.getY(),
                    impactParameter0.getZ(), impactParameter1.getZ(),
                    std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                    std::sqrt(impactParameter0.getSigmaZ2()), std::sqrt(impactParameter1.getSigmaZ2()),
                    prongD0.prong0Id(), prongD0.prong1Id(),
                    prongD0.hfflag());

      if (fillHistograms) {
        registry.fill(HIST("QA/hPtD0"), ptD0);
        registry.fill(HIST("QA/hPtPi"), RecoDecay::pt(pVecSoftPi));
        registry.fill(HIST("QA/hPtD0Prong0"), RecoDecay::pt(pVecD0Prong0));
        registry.fill(HIST("QA/hPtD0Prong1"), RecoDecay::pt(pVecD0Prong1));
        registry.fill(HIST("QA/hPtDstar"), ptDstar);
      }
    }
    // LOG(info) << "Candidate for loop ends";
  }

  ///////////////////////////////////
  ///                             ///
  ///   No centrality selection   ///
  ///                             ///
  ///////////////////////////////////

  /// @brief process function w/ PV refit and w/o centrality selections
  void processPvRefit(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                      aod::Hf2Prongs const& rowsTrackIndexD0,
                      aod::HfDstarsWithPvRefitInfo const& rowsTrackIndexDstar,
                      aod::TracksWCov const& tracks,
                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar</*doPvRefit*/ true, CentralityEstimator::None>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processPvRefit, " Run candidate creator with PV refit and w/o centrality selections", false);

  /// @brief process function w/o PV refit and w/o centrality selections
  void processNoPvRefit(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                        aod::Hf2Prongs const& rowsTrackIndexD0,
                        aod::HfDstars const& rowsTrackIndexDstar,
                        aod::TracksWCov const& tracks,
                        aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar</*doPvRefit*/ false, CentralityEstimator::None>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processNoPvRefit, " Run candidate creator without PV refit and w/o centrality selections", true);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0C   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function w/ PV refit and w/ centrality selection on FT0C
  void processPvRefitCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                              aod::Hf2Prongs const& rowsTrackIndexD0,
                              aod::HfDstarsWithPvRefitInfo const& rowsTrackIndexDstar,
                              aod::TracksWCov const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar</*doPvRefit*/ true, CentralityEstimator::FT0C>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processPvRefitCentFT0C, " Run candidate creator with PV refit nad w/ centrality selection on FT0C", false);

  /// @brief process function w/o PV refit and w/ centrality selection on FT0C
  void processNoPvRefitCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                aod::Hf2Prongs const& rowsTrackIndexD0,
                                aod::HfDstars const& rowsTrackIndexDstar,
                                aod::TracksWCov const& tracks,
                                aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar</*doPvRefit*/ false, CentralityEstimator::FT0C>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processNoPvRefitCentFT0C, " Run candidate creator without PV refit and w centrality selection on FT0C", false);

  /////////////////////////////////////////////
  ///                                       ///
  ///   with centrality selection on FT0M   ///
  ///                                       ///
  /////////////////////////////////////////////

  /// @brief process function w/ PV refit and w/ centrality selection on FT0M
  void processPvRefitCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                              aod::Hf2Prongs const& rowsTrackIndexD0,
                              aod::HfDstarsWithPvRefitInfo const& rowsTrackIndexDstar,
                              aod::TracksWCov const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar</*doPvRefit*/ true, CentralityEstimator::FT0M>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processPvRefitCentFT0M, " Run candidate creator with PV refit nad w/ centrality selection on FT0M", false);

  /// @brief process function w/o PV refit and w/ centrality selection on FT0M
  void processNoPvRefitCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                aod::Hf2Prongs const& rowsTrackIndexD0,
                                aod::HfDstars const& rowsTrackIndexDstar,
                                aod::TracksWCov const& tracks,
                                aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar</*doPvRefit*/ false, CentralityEstimator::FT0M>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processNoPvRefitCentFT0M, " Run candidate creator without PV refit and w/ centrality selection on FT0M", false);

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
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processCollisions, "Collision monitoring - no centrality", true);

  /// @brief process function to monitor collisions - FT0C centrality
  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

  /// @brief process function to monitor collisions - FT0M centrality
  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
};

struct HfCandidateCreatorDstarExpressions {
  Spawns<aod::HfD0FromDstarExt> rowsCandidateD0;
  Spawns<aod::HfCandDstarExt> rowsCandidateDstar;
  Produces<aod::HfCand2ProngMcRec> rowsMcMatchRecD0;
  Produces<aod::HfCand2ProngMcGen> rowsMcMatchGenD0;
  Produces<aod::HfCandDstarMcRec> rowsMcMatchRecDstar;
  Produces<aod::HfCandDstarMcGen> rowsMcMatchGenDstar;

  // Configuration
  Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"};
  Configurable<bool> matchKinkedDecayTopology{"matchKinkedDecayTopology", false, "Match also candidates with tracks that decay with kinked topology"};
  Configurable<bool> matchInteractionsWithMaterial{"matchInteractionsWithMaterial", false, "Match also candidates with tracks that interact with material"};

  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring
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
      if (device.name.compare("hf-candidate-creator-dstar") == 0) {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }
    hfEvSelMc.addHistograms(registry); // particles monitoring
  }

  /// Perform MC Matching.
  template <o2::hf_centrality::CentralityEstimator centEstimator, typename CCs, typename McCollisions>
  void runCreatorDstarMc(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         CCs const& collInfos,
                         McCollisions const& mcCollisions,
                         BCsInfo const&)
  {
    rowsCandidateD0->bindExternalIndices(&tracks);
    rowsCandidateDstar->bindExternalIndices(&tracks);

    int indexRecDstar = -1, indexRecD0 = -1;
    int8_t signDstar = 0, signD0 = 0;
    int8_t flagDstar = 0, flagD0 = 0;
    int8_t originDstar = 0, originD0 = 0;
    int8_t nKinkedTracksDstar = 0, nKinkedTracksD0 = 0;
    int8_t nInteractionsWithMaterialDstar = 0, nInteractionsWithMaterialD0 = 0;

    // Match reconstructed candidates.
    for (const auto& rowCandidateDstar : *rowsCandidateDstar) {
      flagDstar = 0;
      flagD0 = 0;
      originDstar = 0;
      originD0 = 0;
      std::vector<int> idxBhadMothers{};

      auto indexDstar = rowCandidateDstar.globalIndex();
      auto candD0 = rowsCandidateD0->iteratorAt(indexDstar);
      auto candSoftPi = rowCandidateDstar.prongPi_as<aod::TracksWMc>();

      auto arrayDaughtersDstar = std::array{candSoftPi, candD0.prong0_as<aod::TracksWMc>(), candD0.prong1_as<aod::TracksWMc>()};
      auto arrayDaughtersofD0 = std::array{candD0.prong0_as<aod::TracksWMc>(), candD0.prong1_as<aod::TracksWMc>()};

      // Check whether the particle is from background events. If so, reject it.
      if (rejectBackground) {
        bool fromBkg{false};
        for (const auto& daughter : arrayDaughtersDstar) {
          if (daughter.has_mcParticle()) {
            auto mcParticle = daughter.mcParticle();
            if (mcParticle.fromBackgroundEvent()) {
              fromBkg = true;
              break;
            }
          }
        }
        if (fromBkg) {
          rowsMcMatchRecDstar(flagDstar, originDstar, -1.f, 0, 0, 0);
          continue;
        }
      }

      if (matchKinkedDecayTopology && matchInteractionsWithMaterial) {
        // D*± → D0(bar) π±
        indexRecDstar = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughtersDstar, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &signDstar, 2, &nKinkedTracksDstar, &nInteractionsWithMaterialDstar);
        // D0(bar) → π± K∓
        indexRecD0 = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughtersofD0, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0, 1, &nKinkedTracksD0, &nInteractionsWithMaterialD0);
      } else if (matchKinkedDecayTopology && !matchInteractionsWithMaterial) {
        // D*± → D0(bar) π±
        indexRecDstar = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughtersDstar, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &signDstar, 2, &nKinkedTracksDstar);
        // D0(bar) → π± K∓
        indexRecD0 = RecoDecay::getMatchedMCRec<false, false, false, true, false>(mcParticles, arrayDaughtersofD0, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0, 1, &nKinkedTracksD0);
      } else if (!matchKinkedDecayTopology && matchInteractionsWithMaterial) {
        // D*± → D0(bar) π±
        indexRecDstar = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughtersDstar, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &signDstar, 2, nullptr, &nInteractionsWithMaterialDstar);
        // D0(bar) → π± K∓
        indexRecD0 = RecoDecay::getMatchedMCRec<false, false, false, false, true>(mcParticles, arrayDaughtersofD0, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0, 1, nullptr, &nInteractionsWithMaterialD0);
      } else {
        // D*± → D0(bar) π±
        indexRecDstar = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersDstar, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &signDstar, 2);
        // D0(bar) → π± K∓
        indexRecD0 = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersofD0, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0);
      }

      if (indexRecDstar > -1) {
        flagDstar = signDstar * (BIT(aod::hf_cand_dstar::DecayType::DstarToD0Pi));
      }
      if (indexRecD0 > -1) {
        flagD0 = signD0 * (BIT(aod::hf_cand_dstar::DecayType::D0ToPiK));
      }

      // check wether the particle is non-promt (from a B0 hadron)
      if (flagDstar != 0) {
        auto particleDstar = mcParticles.iteratorAt(indexRecDstar);
        originDstar = RecoDecay::getCharmHadronOrigin(mcParticles, particleDstar, false, &idxBhadMothers);
      }
      if (flagD0 != 0) {
        auto particleD0 = mcParticles.iteratorAt(indexRecD0);
        originD0 = RecoDecay::getCharmHadronOrigin(mcParticles, particleD0);
      }
      if (originDstar == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowsMcMatchRecDstar(flagDstar, originDstar, bHadMother.pt(), bHadMother.pdgCode(), nKinkedTracksDstar, nInteractionsWithMaterialDstar);
      } else {
        rowsMcMatchRecDstar(flagDstar, originDstar, -1.f, 0, nKinkedTracksDstar, nInteractionsWithMaterialDstar);
      }
      rowsMcMatchRecD0(flagD0, originD0, -1.f, 0, nKinkedTracksD0, nInteractionsWithMaterialD0);
    }

    for (const auto& mcCollision : mcCollisions) {
      // Slice the MC particles table to get the particles for the current MC collision
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
          rowsMcMatchGenDstar(0, 0, -1);
          rowsMcMatchGenD0(0, 0, -1);
        }
        continue;
      }

      // Match generated particles.
      for (const auto& particle : mcParticlesPerMcColl) {
        flagDstar = 0;
        flagD0 = 0;
        originDstar = 0;
        originD0 = 0;
        std::vector<int> idxBhadMothers{};
        // Reject particles from background events
        if (particle.fromBackgroundEvent() && rejectBackground) {
          rowsMcMatchGenDstar(flagDstar, originDstar, -1);
          rowsMcMatchGenD0(flagD0, originD0, -1);
          continue;
        }

        // D*± → D0(bar) π±
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &signDstar, 2)) {
          flagDstar = signDstar * (BIT(aod::hf_cand_dstar::DecayType::DstarToD0Pi));
        }
        // D0(bar) → π± K∓
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0)) {
          flagD0 = signD0 * (BIT(aod::hf_cand_dstar::DecayType::D0ToPiK));
        }

        // check wether the particle is non-promt (from a B0 hadron)
        if (flagDstar != 0) {
          originDstar = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (flagD0 != 0) {
          originD0 = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
        }

        if (originDstar == RecoDecay::OriginType::NonPrompt) {
          rowsMcMatchGenDstar(flagDstar, originDstar, idxBhadMothers[0]);
        } else {
          rowsMcMatchGenDstar(flagDstar, originDstar, -1);
        }
        rowsMcMatchGenD0(flagD0, originD0, -1.);
      }
    }
  }

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisionsNoCents const& collInfos,
                 aod::McCollisions const& mcCollisions,
                 BCsInfo const& BCsInfo)
  {
    runCreatorDstarMc<CentralityEstimator::None>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstarExpressions, processMc, "Process MC - no centrality", false);

  void processMcCentFT0C(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Cs const& collInfos,
                         aod::McCollisions const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreatorDstarMc<CentralityEstimator::FT0C>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstarExpressions, processMcCentFT0C, "Process MC - FT0c centrality", false);

  void processMcCentFT0M(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Ms const& collInfos,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreatorDstarMc<CentralityEstimator::FT0M>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstarExpressions, processMcCentFT0M, "Process MC - FT0m centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorDstarExpressions>(cfgc)};
}
