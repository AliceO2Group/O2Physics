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

/// \file candidateCreatorCascade.cxx
/// \brief Reconstruction of heavy-flavour cascade decay candidates
///
/// \author Chiara Zampolli, <Chiara.Zampolli@cern.ch>, CERN
///         Paul Buehler, <paul.buehler@oeaw.ac.at>, Vienna

#include <string>
#include <memory>
#include <vector>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/trackUtilities.h"

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::hf_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;

/// Reconstruction of heavy-flavour cascade decay candidates
struct HfCandidateCreatorCascade {
  Produces<aod::HfCandCascBase> rowCandidateBase;

  // vertexing
  // Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "fill validation histograms"};
  Configurable<bool> silenceV0DataWarning{"silenceV0DataWarning", false, "do not print a warning for not found V0s and silently skip them"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber{0};
  double massP{0.};
  double massK0s{0.};
  double massPi{0.};
  double massLc{0.};
  double mass2K0sP{0.};
  double bz = 0.;

  using V0full = soa::Join<aod::V0Datas, aod::V0Covs>;
  using V0fCfull = soa::Join<aod::V0fCDatas, aod::V0fCCovs>;

  std::shared_ptr<TH1> hCandidates;
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

    // histograms
    registry.add("hMass2", "2-prong candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 2.05, 2.55}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    hfEvSel.addHistograms(registry); // collision monitoring

    massP = MassProton;
    massK0s = MassK0Short;
    massPi = MassPiPlus;
    massLc = MassLambdaCPlus;

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

    /// candidate monitoring
    setLabelHistoCands(hCandidates);
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, typename Coll>
  void runCreatorCascade(Coll const&,
                         aod::HfCascades const& rowsTrackIndexCasc,
                         aod::V0sLinked const&,
                         V0full const&,
                         V0fCfull const&,
                         aod::TracksWCov const&,
                         aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {

    // loop over pairs of track indices
    for (const auto& casc : rowsTrackIndexCasc) {

      auto collision = casc.template collision_as<Coll>();
      /// reject candidates in collisions not satisfying the event selections
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      const auto& bach = casc.prong0_as<aod::TracksWCov>();
      LOGF(debug, "V0 %d in HF cascade %d.", casc.v0Id(), casc.globalIndex());
      if (!casc.has_v0()) {
        LOGF(error, "V0 not there for HF cascade %d. Skipping candidate.", casc.globalIndex());
        continue;
      }

      int posGlobalIndex = -1, negGlobalIndex = -1;
      float v0X, v0Y, v0Z, v0px, v0py, v0pz;
      float v0PosPx, v0PosPy, v0PosPz, v0NegPx, v0NegPy, v0NegPz;
      float dcaV0dau, dcaPosToPV, dcaNegToPV, v0cosPA;
      std::array<float, 21> covV = {0.};

      auto v0index = casc.template v0_as<o2::aod::V0sLinked>();
      if (v0index.has_v0Data()) {
        // this V0 passed both standard V0 and cascade V0 selections
        auto v0row = v0index.template v0Data_as<V0full>();
        const auto& trackV0DaughPos = v0row.posTrack_as<aod::TracksWCov>();
        const auto& trackV0DaughNeg = v0row.negTrack_as<aod::TracksWCov>();
        posGlobalIndex = trackV0DaughPos.globalIndex();
        negGlobalIndex = trackV0DaughNeg.globalIndex();
        v0X = v0row.x();
        v0Y = v0row.y();
        v0Z = v0row.z();
        v0px = v0row.px();
        v0py = v0row.py();
        v0pz = v0row.pz();
        v0PosPx = v0row.pxpos();
        v0PosPy = v0row.pypos();
        v0PosPz = v0row.pzpos();
        v0NegPx = v0row.pxneg();
        v0NegPy = v0row.pyneg();
        v0NegPz = v0row.pzneg();
        dcaV0dau = v0row.dcaV0daughters();
        dcaPosToPV = v0row.dcapostopv();
        dcaNegToPV = v0row.dcanegtopv();
        v0cosPA = v0row.v0cosPA();

        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covV[MomInd[i]] = v0row.momentumCovMat()[i];
          covV[i] = v0row.positionCovMat()[i];
        }
      } else if (v0index.has_v0fCData()) {
        // this V0 passes only V0-for-cascade selections, use that instead
        auto v0row = v0index.template v0fCData_as<V0fCfull>();
        const auto& trackV0DaughPos = v0row.posTrack_as<aod::TracksWCov>();
        const auto& trackV0DaughNeg = v0row.negTrack_as<aod::TracksWCov>();
        posGlobalIndex = trackV0DaughPos.globalIndex();
        negGlobalIndex = trackV0DaughNeg.globalIndex();
        v0X = v0row.x();
        v0Y = v0row.y();
        v0Z = v0row.z();
        v0px = v0row.px();
        v0py = v0row.py();
        v0pz = v0row.pz();
        v0PosPx = v0row.pxpos();
        v0PosPy = v0row.pypos();
        v0PosPz = v0row.pzpos();
        v0NegPx = v0row.pxneg();
        v0NegPy = v0row.pyneg();
        v0NegPz = v0row.pzneg();
        dcaV0dau = v0row.dcaV0daughters();
        dcaPosToPV = v0row.dcapostopv();
        dcaNegToPV = v0row.dcanegtopv();
        v0cosPA = v0row.v0cosPA();

        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covV[MomInd[i]] = v0row.momentumCovMat()[i];
          covV[i] = v0row.positionCovMat()[i];
        }
      } else {
        LOGF(warning, "V0Data/V0fCData not there for V0 %d in HF cascade %d. Skipping candidate.", casc.v0Id(), casc.globalIndex());
        continue; // this was inadequately linked, should not happen
      }

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
      }
      df.setBz(bz);

      auto trackBach = getTrackParCov(bach);
      const std::array<float, 3> vertexV0 = {v0X, v0Y, v0Z};
      const std::array<float, 3> momentumV0 = {v0px, v0py, v0pz};
      // we build the neutral track to then build the cascade
      auto trackV0 = o2::track::TrackParCov(vertexV0, momentumV0, covV, 0, true);
      trackV0.setAbsCharge(0);
      trackV0.setPID(o2::track::PID::K0);

      // reconstruct the cascade secondary vertex
      hCandidates->Fill(SVFitting::BeforeFit);
      try {
        if (df.process(trackV0, trackBach) == 0) {
          continue;
        } else {
          LOG(debug) << "Vertexing succeeded for Lc candidate";
        }
      } catch (const std::runtime_error& error) {
        LOG(debug) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        hCandidates->Fill(SVFitting::Fail);
        continue;
      }
      hCandidates->Fill(SVFitting::FitOk);

      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2PCA = df.getChi2AtPCACandidate();
      auto covMatrixPCA = df.calcPCACovMatrixFlat();
      registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]); // FIXME: Calculation of errorDecayLength(XY) gives wrong values without this line.
      // do I have to call "df.propagateTracksToVertex();"?
      auto trackParVarV0 = df.getTrack(0);
      auto trackParVarBach = df.getTrack(1);

      // get track momenta
      std::array<float, 3> pVecV0;
      std::array<float, 3> pVecBach;
      trackParVarV0.getPxPyPzGlo(pVecV0);
      trackParVarBach.getPxPyPzGlo(pVecBach);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
      o2::dataformats::DCA impactParameterV0;
      o2::dataformats::DCA impactParameterBach;
      trackParVarV0.propagateToDCA(primaryVertex, bz, &impactParameterV0); // we do this wrt the primary vtx
      trackParVarBach.propagateToDCA(primaryVertex, bz, &impactParameterBach);

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2PCA,
                       pVecBach[0], pVecBach[1], pVecBach[2],
                       pVecV0[0], pVecV0[1], pVecV0[2],
                       impactParameterBach.getY(), impactParameterV0.getY(),
                       std::sqrt(impactParameterBach.getSigmaY2()), std::sqrt(impactParameterV0.getSigmaY2()),
                       casc.prong0Id(), casc.v0Id(),
                       v0X, v0Y, v0Z,
                       posGlobalIndex, negGlobalIndex,
                       v0PosPx, v0PosPy, v0PosPz,
                       v0NegPx, v0NegPy, v0NegPz,
                       dcaV0dau,
                       dcaPosToPV,
                       dcaNegToPV,
                       v0cosPA);

      // fill histograms
      if (fillHistograms) {
        // calculate invariant masses
        mass2K0sP = RecoDecay::m(std::array{pVecBach, pVecV0}, std::array{massP, massK0s});
        registry.fill(HIST("hMass2"), mass2K0sP);
      }
    }

    return;
  }

  /// @brief process function w/o centrality selections
  void processNoCent(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                     aod::HfCascades const& rowsTrackIndexCasc,
                     aod::V0sLinked const& v0sLinked,
                     V0full const& v0Full,
                     V0fCfull const& v0fcFull,
                     aod::TracksWCov const& tracks,
                     aod::BCsWithTimestamps const& bcs)
  {
    runCreatorCascade<CentralityEstimator::None>(collisions, rowsTrackIndexCasc, v0sLinked, v0Full, v0fcFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascade, processNoCent, " Run candidate creator w/o centrality selections", true);

  /// @brief process function w/ centrality selection on FT0C
  void processCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                       aod::HfCascades const& rowsTrackIndexCasc,
                       aod::V0sLinked const& v0sLinked,
                       V0full const& v0Full,
                       V0fCfull const& v0fcFull,
                       aod::TracksWCov const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    runCreatorCascade<CentralityEstimator::FT0C>(collisions, rowsTrackIndexCasc, v0sLinked, v0Full, v0fcFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCentFT0C, " Run candidate creator w/ centrality selection on FT0C", false);

  /// @brief process function w/ centrality selection on FT0M
  void processCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                       aod::HfCascades const& rowsTrackIndexCasc,
                       aod::V0sLinked const& v0sLinked,
                       V0full const& v0Full,
                       V0fCfull const& v0fcFull,
                       aod::TracksWCov const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    runCreatorCascade<CentralityEstimator::FT0M>(collisions, rowsTrackIndexCasc, v0sLinked, v0Full, v0fcFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCentFT0M, " Run candidate creator w/ centrality selection on FT0M", false);

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
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCollisions, "Collision monitoring - no centrality", true);

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
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

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
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
};

/// Performs MC matching.
struct HfCandidateCreatorCascadeMc {
  Spawns<aod::HfCandCascExt> rowCandidateCasc;
  Produces<aod::HfCandCascadeMcRec> rowMcMatchRec;
  Produces<aod::HfCandCascadeMcGen> rowMcMatchGen;

  // Configuration
  Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"};

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring

  using MyTracksWMc = soa::Join<aod::TracksWCov, aod::McTrackLabels>;
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

  // inspect for which zPvPosMax cut was set for reconstructed
  void init(InitContext& initContext)
  {
    std::array<bool, 3> procCollisions = {doprocessMc, doprocessMcCentFT0C, doprocessMcCentFT0M};
    if (std::accumulate(procCollisions.begin(), procCollisions.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for collision study can be enabled at a time.");
    }

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-cascade") == 0) {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }
    hfEvSelMc.addHistograms(registry); // particles monitoring
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, typename CCs, typename McCollisions>
  void runCreatorCascMc(MyTracksWMc const& tracks,
                        aod::McParticles const& mcParticles,
                        CCs const& collInfos,
                        McCollisions const& mcCollisions,
                        BCsInfo const&)
  {
    // Match reconstructed candidates.
    rowCandidateCasc->bindExternalIndices(&tracks);
    for (const auto& candidate : *rowCandidateCasc) {

      int8_t sign = 0;
      int8_t origin = 0;
      int indexRec = -1;

      std::vector<int> idxBhadMothers{};

      const auto& bach = candidate.prong0_as<MyTracksWMc>();
      const auto& trackV0DaughPos = candidate.posTrack_as<MyTracksWMc>();
      const auto& trackV0DaughNeg = candidate.negTrack_as<MyTracksWMc>();

      auto arrayDaughtersV0 = std::array{trackV0DaughPos, trackV0DaughNeg};
      auto arrayDaughtersLc = std::array{bach, trackV0DaughPos, trackV0DaughNeg};
      // Check whether the particle is from background events. If so, reject it.
      if (rejectBackground) {
        bool fromBkg{false};
        for (const auto& daughter : arrayDaughtersLc) {
          if (daughter.has_mcParticle()) {
            auto mcParticle = daughter.mcParticle();
            if (mcParticle.fromBackgroundEvent()) {
              fromBkg = true;
              break;
            }
          }
        }
        if (fromBkg) {
          rowMcMatchRec(sign, origin, -1.f, 0);
          continue;
        }
      }

      int indexK0SRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, kK0Short, std::array{+kPiPlus, -kPiPlus}, false, &sign, 1);
      if (indexK0SRec >= 0) { // we have already positively checked the K0s
        // then we check the Lc
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersLc, Pdg::kLambdaCPlus, std::array{+kProton, +kPiPlus, -kPiPlus}, true, &sign, 3); // 3-levels Lc --> p + K0 --> p + K0s --> p + pi+ pi-
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (indexRec >= 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMcMatchRec(sign, origin, bHadMother.pt(), bHadMother.pdgCode());
      } else {
        rowMcMatchRec(sign, origin, -1.f, 0);
      }
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
          rowMcMatchGen(0, 0, -1);
        }
        continue;
      }

      // Match generated particles.
      for (const auto& particle : mcParticlesPerMcColl) {

        int8_t sign = 0;
        int8_t origin = 0;
        int8_t flag = 0;

        std::vector<int> idxBhadMothers{};
        // Reject particles from background events
        if (particle.fromBackgroundEvent() && rejectBackground) {
          rowMcMatchGen(sign, origin, -1);
          continue;
        }
        // checking if I have a Lc --> K0S + p
        RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, +kK0Short}, false, &sign, 2);
        if (sign == 0) { // now check for anti-Lc
          RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, -Pdg::kLambdaCPlus, std::array{-kProton, +kK0Short}, false, &sign, 2);
          sign = -sign;
        }
        if (sign != 0) {
          // we check the K0S
          for (const auto& daughterK0 : particle.template daughters_as<aod::McParticles>()) {
            if (std::abs(daughterK0.pdgCode()) != kK0) {
              continue;
            }
            for (const auto& daughterK0S : daughterK0.template daughters_as<aod::McParticles>()) {
              if (daughterK0S.pdgCode() != kK0Short) {
                continue;
              }
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterK0S, kK0Short, std::array{+kPiPlus, -kPiPlus}, true)) {
                flag = sign;
              }
            }
          }
        }
        // Check whether the particle is non-prompt (from a b quark).
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          rowMcMatchGen(flag, origin, idxBhadMothers[0]);
        } else {
          rowMcMatchGen(flag, origin, -1);
        }
      }
    }
  }

  void processMc(MyTracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisionsNoCents const& collInfos,
                 aod::McCollisions const& mcCollisions,
                 BCsInfo const& BCsInfo)
  {
    runCreatorCascMc<CentralityEstimator::None>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascadeMc, processMc, "Process MC - no centrality", false);

  void processMcCentFT0C(MyTracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Cs const& collInfos,
                         aod::McCollisions const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreatorCascMc<CentralityEstimator::FT0C>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascadeMc, processMcCentFT0C, "Process MC - FT0c centrality", false);

  void processMcCentFT0M(MyTracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsFT0Ms const& collInfos,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         BCsInfo const& BCsInfo)
  {
    runCreatorCascMc<CentralityEstimator::FT0M>(tracks, mcParticles, collInfos, mcCollisions, BCsInfo);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascadeMc, processMcCentFT0M, "Process MC - FT0m centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorCascade>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorCascadeMc>(cfgc)};
}
