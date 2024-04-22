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

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/trackUtilities.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::hf_trkcandsel;
using namespace o2::aod::hf_collision_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;

/// Reconstruction of heavy-flavour cascade decay candidates
struct HfCandidateCreatorCascade {
  Produces<aod::HfCandCascBase> rowCandidateBase;

  // centrality
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality"};
  // event selection
  Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "apply the sel8 event selection"};
  Configurable<float> zPvPosMax{"zPvPosMax", 10.f, "max. PV posZ (cm)"};
  Configurable<bool> useTimeFrameBorderCut{"useTimeFrameBorderCut", true, "apply TF border cut"};
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

  std::shared_ptr<TH1> hCollisions, hPosZBeforeEvSel, hPosZAfterEvSel, hPosXAfterEvSel, hPosYAfterEvSel, hNumPvContributorsAfterSel, hCandidates;
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
    hCollisions = registry.add<TH1>("hCollisions", "HF event counter;;entries", {HistType::kTH1D, {axisEvents}});
    hPosZBeforeEvSel = registry.add<TH1>("hPosZBeforeEvSel", "all events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{400, -20., 20.}}});
    hPosZAfterEvSel = registry.add<TH1>("hPosZAfterEvSel", "selected events;#it{z}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{400, -20., 20.}}});
    hPosXAfterEvSel = registry.add<TH1>("hPosXAfterEvSel", "selected events;#it{x}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{200, -0.5, 0.5}}});
    hPosYAfterEvSel = registry.add<TH1>("hPosYAfterEvSel", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{200, -0.5, 0.5}}});
    hNumPvContributorsAfterSel = registry.add<TH1>("hNumPvContributorsAfterSel", "selected events;#it{y}_{prim. vtx.} (cm);entries", {HistType::kTH1D, {{500, -0.5, 499.5}}});
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});

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

    /// collision monitoring
    setLabelHistoEvSel(hCollisions);

    /// candidate monitoring
    setLabelHistoCands(hCandidates);
  }

  template <o2::aod::hf_collision_centrality::CentralityEstimator centEstimator, typename Coll>
  void runCreatorCascade(Coll const&,
                         aod::HfCascades const& rowsTrackIndexCasc,
                         aod::V0sLinked const&,
                         aod::V0Datas const&,
                         aod::V0fCDatas const&,
                         aod::TracksWCov const&,
                         aod::BCsWithTimestamps const&)
  {
    // loop over pairs of track indices
    for (const auto& casc : rowsTrackIndexCasc) {

      /// reject candidates in collisions not satisfying the event selections
      auto collision = casc.template collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = getHfCollisionRejectionMask<true, centEstimator>(collision, centrality, centralityMin, centralityMax, useSel8Trigger, -1, useTimeFrameBorderCut, -zPvPosMax, zPvPosMax, 0, -1.f);
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
      float v0x, v0y, v0z, v0px, v0py, v0pz;
      float v0PosPx, v0PosPy, v0PosPz, v0NegPx, v0NegPy, v0NegPz;
      float dcaV0dau, dcaPosToPV, dcaNegToPV, v0cosPA;
      float posTrackX, negTrackX;
      o2::track::TrackParCov trackParCovV0DaughPos;
      o2::track::TrackParCov trackParCovV0DaughNeg;

      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (v0index.has_v0Data()) {
        // this V0 passed both standard V0 and cascade V0 selections
        auto v0row = v0index.v0Data();
        const auto& trackV0DaughPos = v0row.posTrack_as<aod::TracksWCov>();
        const auto& trackV0DaughNeg = v0row.negTrack_as<aod::TracksWCov>();
        trackParCovV0DaughPos = getTrackParCov(trackV0DaughPos); // check that aod::TracksWCov does not need TracksDCA!
        trackParCovV0DaughNeg = getTrackParCov(trackV0DaughNeg); // check that aod::TracksWCov does not need TracksDCA!
        posGlobalIndex = trackV0DaughPos.globalIndex();
        negGlobalIndex = trackV0DaughNeg.globalIndex();
        v0x = v0row.x();
        v0y = v0row.y();
        v0z = v0row.z();
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
        posTrackX = v0row.posX();
        negTrackX = v0row.negX();
      } else if (v0index.has_v0fCData()) {
        // this V0 passes only V0-for-cascade selections, use that instead
        auto v0row = v0index.v0fCData();
        const auto& trackV0DaughPos = v0row.posTrack_as<aod::TracksWCov>();
        const auto& trackV0DaughNeg = v0row.negTrack_as<aod::TracksWCov>();
        trackParCovV0DaughPos = getTrackParCov(trackV0DaughPos); // check that aod::TracksWCov does not need TracksDCA!
        trackParCovV0DaughNeg = getTrackParCov(trackV0DaughNeg); // check that aod::TracksWCov does not need TracksDCA!
        posGlobalIndex = trackV0DaughPos.globalIndex();
        negGlobalIndex = trackV0DaughNeg.globalIndex();
        v0x = v0row.x();
        v0y = v0row.y();
        v0z = v0row.z();
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
        posTrackX = v0row.posX();
        negTrackX = v0row.negX();
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

      auto trackParCovBach = getTrackParCov(bach);
      trackParCovV0DaughPos.propagateTo(posTrackX, bz); // propagate the track to the X closest to the V0 vertex
      trackParCovV0DaughNeg.propagateTo(negTrackX, bz); // propagate the track to the X closest to the V0 vertex
      const std::array<float, 3> vertexV0 = {v0x, v0y, v0z};
      const std::array<float, 3> momentumV0 = {v0px, v0py, v0pz};
      // we build the neutral track to then build the cascade
      auto trackV0 = o2::dataformats::V0(vertexV0, momentumV0, {0, 0, 0, 0, 0, 0}, trackParCovV0DaughPos, trackParCovV0DaughNeg); // build the V0 track (indices for v0 daughters set to 0 for now)

      // reconstruct the cascade secondary vertex
      hCandidates->Fill(SVFitting::BeforeFit);
      try {
        if (df.process(trackV0, trackParCovBach) == 0) {
          continue;
        } else {
          // LOG(info) << "Vertexing succeeded for Lc candidate";
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCFitterN cannot work, skipping the candidate.";
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
                       v0x, v0y, v0z,
                       // v0.posTrack(), v0.negTrack(), // why this was not fine?
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
                     aod::V0Datas const& v0Data,
                     aod::V0fCDatas const& v0fCDatas,
                     aod::TracksWCov const& tracks,
                     aod::BCsWithTimestamps const& bcs)
  {
    runCreatorCascade<CentralityEstimator::None>(collisions, rowsTrackIndexCasc, v0sLinked, v0Data, v0fCDatas, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascade, processNoCent, " Run candidate creator w/o centrality selections", true);

  /// @brief process function w/ centrality selection on FT0C
  void processCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                       aod::HfCascades const& rowsTrackIndexCasc,
                       aod::V0sLinked const& v0sLinked,
                       aod::V0Datas const& v0Data,
                       aod::V0fCDatas const& v0fCDatas,
                       aod::TracksWCov const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    runCreatorCascade<CentralityEstimator::FT0C>(collisions, rowsTrackIndexCasc, v0sLinked, v0Data, v0fCDatas, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCentFT0C, " Run candidate creator w/ centrality selection on FT0C", false);

  /// @brief process function w/ centrality selection on FT0M
  void processCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                       aod::HfCascades const& rowsTrackIndexCasc,
                       aod::V0sLinked const& v0sLinked,
                       aod::V0Datas const& v0Data,
                       aod::V0fCDatas const& v0fCDatas,
                       aod::TracksWCov const& tracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    runCreatorCascade<CentralityEstimator::FT0M>(collisions, rowsTrackIndexCasc, v0sLinked, v0Data, v0fCDatas, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCentFT0M, " Run candidate creator w/ centrality selection on FT0M", false);

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
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCollisions, "Collision monitoring - no centrality", true);

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
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

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
  PROCESS_SWITCH(HfCandidateCreatorCascade, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
};

/// Performs MC matching.
struct HfCandidateCreatorCascadeMc {
  Spawns<aod::HfCandCascExt> rowCandidateCasc;
  Produces<aod::HfCandCascadeMcRec> rowMcMatchRec;
  Produces<aod::HfCandCascadeMcGen> rowMcMatchGen;

  using MyTracksWMc = soa::Join<aod::TracksWCov, aod::McTrackLabels>;

  void processMc(MyTracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    int8_t sign = 0;
    int8_t origin = 0;
    int indexRec = -1;
    std::vector<int> arrDaughLcIndex;
    std::array<int, 3> arrDaughLcPDG;
    std::array<int, 3> arrDaughLcPDGRef = {+kProton, +kPiPlus, -kPiPlus};

    // Match reconstructed candidates.
    rowCandidateCasc->bindExternalIndices(&tracks);
    for (const auto& candidate : *rowCandidateCasc) {

      origin = 0;

      const auto& bach = candidate.prong0_as<MyTracksWMc>();
      const auto& trackV0DaughPos = candidate.posTrack_as<MyTracksWMc>();
      const auto& trackV0DaughNeg = candidate.negTrack_as<MyTracksWMc>();

      auto arrayDaughtersV0 = std::array{trackV0DaughPos, trackV0DaughNeg};
      auto arrayDaughtersLc = std::array{bach, trackV0DaughPos, trackV0DaughNeg};

      // First we check the K0s
      LOG(debug) << "\n";
      LOG(debug) << "Checking MC for candidate!";
      LOG(debug) << "Looking for K0s";

      // if (isLc) {
      RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, kK0Short, std::array{+kPiPlus, -kPiPlus}, false, &sign, 1);
      if (sign != 0) { // we have already positively checked the K0s
        // then we check the Lc
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersLc, Pdg::kLambdaCPlus, std::array{+kProton, +kPiPlus, -kPiPlus}, true, &sign, 3); // 3-levels Lc --> p + K0 --> p + K0s --> p + pi+ pi-
      }

      // Check whether the particle is non-prompt (from a b quark).
      if (sign != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }

      rowMcMatchRec(sign, origin);
    }
    //}

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      origin = 0;
      // checking if I have a Lc --> K0S + p
      RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, +kK0Short}, false, &sign, 2);
      if (sign == 0) { // now check for anti-Lc
        RecoDecay::isMatchedMCGen(mcParticles, particle, -Pdg::kLambdaCPlus, std::array{-kProton, +kK0Short}, false, &sign, 2);
        sign = -sign;
      }
      if (sign != 0) {
        arrDaughLcIndex.clear();
        // checking that the final daughters (decay depth = 3) are p, pi+, pi-
        RecoDecay::getDaughters(particle, &arrDaughLcIndex, arrDaughLcPDGRef, 3); // best would be to check the K0S daughters
        if (arrDaughLcIndex.size() == 3) {
          for (std::size_t iProng = 0; iProng < arrDaughLcIndex.size(); ++iProng) {
            auto daughI = mcParticles.rawIteratorAt(arrDaughLcIndex[iProng]);
            arrDaughLcPDG[iProng] = daughI.pdgCode();
          }
          if (!(arrDaughLcPDG[0] == sign * arrDaughLcPDGRef[0] && arrDaughLcPDG[1] == arrDaughLcPDGRef[1] && arrDaughLcPDG[2] == arrDaughLcPDGRef[2])) { // this should be the condition, first bach, then v0
            sign = 0;
          } else {
            LOG(debug) << "Lc --> K0S+p found in MC table";
          }
        }
      }
      // Check whether the particle is non-prompt (from a b quark).
      if (sign != 0) {
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
      }
      rowMcMatchGen(sign, origin);
    }
  }

  PROCESS_SWITCH(HfCandidateCreatorCascadeMc, processMc, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorCascade>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorCascadeMc>(cfgc)};
}
