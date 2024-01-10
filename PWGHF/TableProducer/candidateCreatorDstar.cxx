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

// ROOT
#include <TPDGCode.h>
// O2
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
// O2Physics
#include "Common/Core/trackUtilities.h"
// PWGHF
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
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

  Service<o2::ccdb::BasicCCDBManager> ccdb; // From utilsBfieldCCDB.h
  o2::base::MatLayerCylSet* lut;            // From MatLayercylSet.h
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  // D0-prong vertex fitter
  o2::vertexing::DCAFitterN<2> df;
  int runNumber;
  double bz;
  static constexpr float CmToMicrometers = 10000.; // from cm to µm
  double massPi, massK, massD0;

  AxisSpec ptAxis = {100, 0., 2.0, "#it{p}_{T} (GeV/#it{c}"};
  AxisSpec dcaAxis = {200, -500., 500., "#it{d}_{xy,z} (#mum)"};

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
    if (doprocessPvRefit && doprocessNoPvRefit) { //............Warning! remove this if any of this function is removed
      LOGP(fatal, "Only one process function between processPvRefit and processNoPvRefit can be enabled at a time.");
    }
    // LOG(info) << "Init Function Invoked";
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking(); // set the flag to check object validity before CCDB query
    // LOG(info) << "Retriving ccdb object";
    auto rectification = ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut); // retrieve an object of type T from CCDB as stored under path; will use the timestamp member
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(rectification);
    // LOG(info) << "Successfully Retrived";
    runNumber = 0;
    bz = 0;

    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    massPi = MassPiPlus;
    massK = MassKPlus;
    massD0 = MassD0;
  }

  /// @brief function for secondary vertex reconstruction and candidate creator
  /// @tparam CandsDstar Table type of Dstar table object
  /// @tparam doPvRefit True/False PV refit option
  /// @param collisions collision object
  /// @param rowsTrackIndexDstar Dstar table object from trackIndexSkimCreator.cxx
  /// @param rowsTrackIndexD0 D0 table object from trackIndexSkimCreator.cxx
  /// @param tracks track table with Cov object
  /// @param bcWithTimeStamps Bunch Crossing with timestamps
  template <bool doPvRefit, typename CandsDstar>
  void runCreatorDstar(aod::Collisions const& collisions,
                       CandsDstar const& rowsTrackIndexDstar,
                       aod::Hf2Prongs const& rowsTrackIndexD0,
                       aod::TracksWCov const& tracks,
                       aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // LOG(info) << "runCreatorDstar function called";
    // LOG(info) << "candidate loop starts";
    // loop over suspected Dstar Candidate
    for (const auto& rowTrackIndexDstar : rowsTrackIndexDstar) {

      auto trackPi = rowTrackIndexDstar.template prong0_as<aod::TracksWCov>();
      auto prongD0 = rowTrackIndexDstar.template prongD0_as<aod::Hf2Prongs>();
      auto trackD0Prong0 = prongD0.template prong0_as<aod::TracksWCov>();
      auto trackD0Prong1 = prongD0.template prong1_as<aod::TracksWCov>();

      auto collision = rowTrackIndexDstar.collision();

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
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2); // Sets up the grp object for magnetic field (w/o matCorr for propagation)
        bz = o2::base::Propagator::Instance()->getNominalBz();
        // LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      df.setBz(bz);

      // reconstruct the 2-prong secondary vertex
      if (df.process(trackD0Prong0ParVarPos1, trackD0Prong1ParVarNeg1) == 0) {
        continue;
      }
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

      // Soft pi momentum vector
      std::array<float, 3> pVecSoftPi;
      trackPiParVar.getPxPyPzGlo(pVecSoftPi);

      // D* pt magnitude
      auto ptDstar = RecoDecay::pt(pVecD0, pVecSoftPi);

      // Fill candidate Table for DStar
      rowCandDstarBase(collision.globalIndex(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       rowTrackIndexDstar.prong0Id(), rowTrackIndexDstar.prongD0Id(),
                       pVecSoftPi[0], pVecSoftPi[1], pVecSoftPi[2],
                       impactParameterPi.getY(), std::sqrt(impactParameterPi.getSigmaY2()),
                       pVecD0Prong0[0], pVecD0Prong0[1], pVecD0Prong0[2],
                       pVecD0Prong1[0], pVecD0Prong1[1], pVecD0Prong1[2]);
      // Fill candidate Table for D0
      rowCandD0Base(collision.globalIndex(),
                    primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                    secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                    errorDecayLength, errorDecayLengthXY,
                    chi2PCA,
                    pVecD0Prong0[0], pVecD0Prong0[1], pVecD0Prong0[2],
                    pVecD0Prong1[0], pVecD0Prong1[1], pVecD0Prong1[2],
                    impactParameter0.getY(), impactParameter1.getY(),
                    std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
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

  void processPvRefit(aod::Collisions const& collisions,
                      aod::Hf2Prongs const& rowsTrackIndexD0,
                      aod::HfDstarsWithPvRefitInfo const& rowsTrackIndexDstar,
                      aod::TracksWCov const& tracks,
                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar<true>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processPvRefit, " Run candidate creator with PV refit", false);

  void processNoPvRefit(aod::Collisions const& collisions,
                        aod::Hf2Prongs const& rowsTrackIndexD0,
                        aod::HfDstars const& rowsTrackIndexDstar,
                        aod::TracksWCov const& tracks,
                        aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorDstar<false>(collisions, rowsTrackIndexDstar, rowsTrackIndexD0, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorDstar, processNoPvRefit, " Run candidate creator without PV refit", true);
};

struct HfCandidateCreatorDstarExpressions {
  Spawns<aod::HfD0FromDstarExt> rowsCandidateD0;
  Produces<aod::HfCand2ProngMcRec> rowsMcMatchRecD0;
  Produces<aod::HfCand2ProngMcGen> rowsMcMatchGenD0;

  Spawns<aod::HfCandDstarExt> rowsCandidateDstar;
  Produces<aod::HfCandDstarMcRec> rowsMcMatchRecDstar;
  Produces<aod::HfCandDstarMcGen> rowsMcMatchGenDstar;

  void init(InitContext const&) {}

  /// Perform MC Matching.
  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    rowsCandidateD0->bindExternalIndices(&tracks);
    rowsCandidateDstar->bindExternalIndices(&tracks);

    int indexRecDstar = -1, indexRecD0 = -1;
    int8_t signDstar = 0, signD0 = 0;
    int8_t flagDstar = 0, flagD0 = 0;
    int8_t originDstar = 0, originD0 = 0;

    // Match reconstructed candidates.
    for (const auto& rowCandidateDstar : *rowsCandidateDstar) {
      flagDstar = 0;
      flagD0 = 0;
      originDstar = 0;
      originD0 = 0;

      auto indexDstar = rowCandidateDstar.globalIndex();
      auto candD0 = rowsCandidateD0->iteratorAt(indexDstar);
      auto candSoftPi = rowCandidateDstar.prongPi_as<aod::TracksWMc>();

      auto arrayDaughtersDstar = std::array{candSoftPi, candD0.prong0_as<aod::TracksWMc>(), candD0.prong1_as<aod::TracksWMc>()};
      auto arrayDaughtersofD0 = std::array{candD0.prong0_as<aod::TracksWMc>(), candD0.prong1_as<aod::TracksWMc>()};

      // D*± → D0(bar) π±
      indexRecDstar = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersDstar, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &signDstar, 2);
      // D0(bar) → π± K∓
      indexRecD0 = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersofD0, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &signD0);

      if (indexRecDstar > -1) {
        flagDstar = signDstar * (BIT(aod::hf_cand_dstar::DecayType::DstarToD0Pi));
      }
      if (indexRecD0 > -1) {
        flagD0 = signD0 * (BIT(aod::hf_cand_dstar::DecayType::D0ToPiK));
      }

      // check wether the particle is non-promt (from a B0 hadron)
      if (flagDstar != 0) {
        auto particleDstar = mcParticles.iteratorAt(indexRecDstar);
        originDstar = RecoDecay::getCharmHadronOrigin(mcParticles, particleDstar);
      }
      if (flagD0 != 0) {
        auto particleD0 = mcParticles.iteratorAt(indexRecD0);
        originD0 = RecoDecay::getCharmHadronOrigin(mcParticles, particleD0);
      }
      rowsMcMatchRecDstar(flagDstar, originDstar);
      rowsMcMatchRecD0(flagD0, originD0);
    }

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flagDstar = 0;
      flagD0 = 0;
      originDstar = 0;
      originD0 = 0;

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
        originDstar = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }
      if (flagD0 != 0) {
        originD0 = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }
      rowsMcMatchGenDstar(flagDstar, originDstar);
      rowsMcMatchGenD0(flagD0, originD0);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorDstarExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorDstar>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorDstarExpressions>(cfgc)};
}
