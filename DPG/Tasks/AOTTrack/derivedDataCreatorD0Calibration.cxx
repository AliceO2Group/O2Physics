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

/// \file derivedDataCreatorD0Calibration.cxx
/// \brief Producer of derived tables of D0 candidates, daughter tracks and collisions for calibration studies
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "D0CalibTables.h"

#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsPid.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"

#include "CommonDataFormat/InteractionRecord.h"
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/AnalysisTask.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/detail/TypeTruncation.h>
#include <ReconstructionDataFormats/DCA.h>

#include <TH1D.h>
#include <TRandom3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_calib;

struct DerivedDataCreatorD0Calibration {

  Produces<aod::D0CalibColls> collTable;
  Produces<aod::D0CalibTracks> trackTable;
  Produces<aod::D0CalibCands> candTable;

  struct : ConfigurableGroup {
    Configurable<float> ptMin{"ptMin", 0.4, "min. track pT"};
    Configurable<float> absEtaMax{"absEtaMax", 1., "max. track absolute eta"};
    Configurable<std::vector<float>> binsPt{"binsPt", std::vector<float>{hf_calib::vecBinsPtTrack}, "track pT bin limits for DCA pT-dependent cut"};
    Configurable<LabeledArray<float>> limitsDca{"limitsDca", {hf_calib::CutsTrack[0], hf_calib::NBinsPtTrack, hf_calib::NCutVarsTrack, hf_calib::labelsPtTrack, hf_calib::labelsCutVarTrack}, "Single-track selections per pT bin"};
    // TPC PID
    Configurable<float> ptPidTpcMin{"ptPidTpcMin", 0., "Lower bound of track pT for TPC PID"};
    Configurable<float> ptPidTpcMax{"ptPidTpcMax", 1000., "Upper bound of track pT for TPC PID"};
    Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
    Configurable<bool> usePidTpcOnly{"usePidTpcOnly", false, "Only use TPC PID"};
    // TOF PID
    Configurable<float> ptPidTofMin{"ptPidTofMin", 0., "Lower bound of track pT for TOF PID"};
    Configurable<float> ptPidTofMax{"ptPidTofMax", 1000., "Upper bound of track pT for TOF PID"};
    Configurable<float> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
    std::string prefix = "trackCuts";
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    Configurable<float> ptMin{"ptMin", 0., "min. D0-candidate pT"};
    Configurable<std::vector<float>> binsPt{"binsPt", std::vector<float>{hf_calib::vecBinsPtCand}, "pT bin limits"};
    Configurable<LabeledArray<float>> topologicalCuts{"topologicalCuts", {hf_calib::CutsCand[0], hf_calib::NBinsPtCand, hf_calib::NCutVarsCand, hf_calib::labelsPtCand, hf_calib::labelsCutVarCand}, "D0 candidate selection per pT bin"};
    std::string prefix = "candidateCuts";
  } cfgCandCuts;

  struct : ConfigurableGroup {
    Configurable<bool> apply{"apply", false, "flag to apply downsampling"};
    Configurable<std::string> pathCcdbWeights{"pathCcdbWeights", "", "CCDB path containing pT-differential weights"};
    std::string prefix = "downsampling";
  } cfgDownsampling;

  struct : ConfigurableGroup {
    Configurable<bool> apply{"apply", false, "flag to apply downsampling"};
    Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_calib::vecBinsPtMl}, "pT bin limits for ML models inference"};
    Configurable<LabeledArray<double>> thresholdMlScores{"thresholdMlScores", {hf_calib::CutsMl[0], hf_calib::NBinsPtMl, 3, hf_calib::labelsPtMl, hf_calib::labelsCutMl}, "Threshold values for Ml output scores of D0 candidates"};
    Configurable<bool> loadMlModelsFromCCDB{"loadMlModelsFromCCDB", true, "Flag to enable or disable the loading of ML models from CCDB"};
    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/f/fgrosa/D0Calib/BDT/Pt0_1"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_pT_0_1.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    std::string prefix = "ml";
  } cfgMl;

  using TracksWCovExtraPid = soa::Join<aod::Tracks, aod::TrackToTmo, aod::TracksCov, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using TracksWCovExtraPidAndQa = soa::Join<aod::Tracks, aod::TrackToTmo, aod::TrackToTracksQA, aod::TracksCov, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using CollisionsWEvSel = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  using TrackMeanOccs = soa::Join<aod::TmoTrackIds, aod::TmoPrim, aod::TmoT0V0, aod::TmoRT0V0Prim, aod::TwmoPrim, aod::TwmoT0V0, aod::TwmoRT0V0Prim>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  o2::analysis::MlResponse<float> mlResponse;

  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;

  int runNumber{0};
  double bz{0.};
  const float zVtxMax{10.f};
  // tolerances for preselections before vertex reconstruction
  const float ptTolerance{0.1f};
  const float invMassTolerance{0.05f};
  uint32_t precisionCovMask{0xFFFFE000}; // 10 bits
  uint32_t precisionDcaMask{0xFFFFFC00}; // 13 bits

  OutputObj<TH1D> histDownSampl{"histDownSampl"};

  void init(InitContext const&)
  {
    // First we set the CCDB manager
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    if (cfgDownsampling.apply) {
      histDownSampl.setObject(reinterpret_cast<TH1D*>(ccdb->getSpecific<TH1D>(cfgDownsampling.pathCcdbWeights)));
    }

    if (cfgMl.apply) {
      std::vector<int> cutDir = {o2::cuts_ml::CutDirection::CutGreater, o2::cuts_ml::CutDirection::CutSmaller, o2::cuts_ml::CutDirection::CutSmaller};
      mlResponse.configure(cfgMl.binsPt, cfgMl.thresholdMlScores, cutDir, 3);
      if (cfgMl.loadMlModelsFromCCDB) {
        ccdbApi.init("http://alice-ccdb.cern.ch");
        mlResponse.setModelPathsCCDB(cfgMl.onnxFileNames, ccdbApi, cfgMl.modelPathsCCDB, -1);
      } else {
        mlResponse.setModelPathsLocal(cfgMl.onnxFileNames);
      }
      mlResponse.init();
    }

    df.setPropagateToPCA(true);
    df.setMaxR(200.f);
    df.setMaxDZIni(4.f);
    df.setMinParamChange(1.e-3f);
    df.setMinRelChi2Change(0.9f);
    df.setUseAbsDCA(false);
    df.setWeightedFinalPCA(false);
    df.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE); // we are always inside the beampipe

    selectorPion.setRangePtTpc(cfgTrackCuts.ptPidTpcMin, cfgTrackCuts.ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-cfgTrackCuts.nSigmaTpcMax, cfgTrackCuts.nSigmaTpcMax);
    selectorPion.setRangePtTof(cfgTrackCuts.ptPidTofMin, cfgTrackCuts.ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-cfgTrackCuts.nSigmaTofMax, cfgTrackCuts.nSigmaTofMax);
    selectorKaon = selectorPion;
  }

  // main function
  template <bool withTrackQa, typename TTrackQa, typename TTracks>
  void runDataCreation(CollisionsWEvSel const& collisions,
                       aod::TrackAssoc const& trackIndices,
                       TTracks const&,
                       aod::BCsWithTimestamps const&,
                       TrackMeanOccs const&,
                       TTrackQa const&)
  {
    std::map<int, int> selectedCollisions; // map with indices of selected collisions (key: original AOD Collision table index, value: D0 collision index)
    std::map<int, int> selectedTracks;     // map with indices of selected tracks (key: original AOD Track table index, value: D0 daughter track index)

    for (auto const& collision : collisions) {

      // minimal event selection
      if (!collision.sel8()) {
        continue;
      }
      auto primaryVertex = getPrimaryVertex(collision);
      if (std::abs(primaryVertex.getZ()) > zVtxMax) {
        continue;
      }

      auto covMatrixPV = primaryVertex.getCov();

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        initCCDB(bc, runNumber, ccdb, "GLO/Config/GRPMagField", nullptr, false);
        bz = o2::base::Propagator::Instance()->getNominalBz();
      }
      o2::InteractionRecord eventIR;
      eventIR.setFromLong(bc.globalBC());

      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      for (auto const& trackIndexPos : groupedTrackIndices) {
        auto trackPos = trackIndexPos.template track_as<TTracks>();
        // track selections
        if (trackPos.sign() < 0) { // first positive track
          continue;
        }
        if (!trackPos.isGlobalTrackWoDCA()) {
          continue;
        }
        if (trackPos.pt() < cfgTrackCuts.ptMin) {
          continue;
        }
        if (std::abs(trackPos.eta()) > cfgTrackCuts.absEtaMax) {
          continue;
        }
        auto trackParCovPos = getTrackParCov(trackPos);
        o2::dataformats::DCA dcaPos;
        trackParCovPos.propagateToDCA(primaryVertex, bz, &dcaPos);
        if (!isSelectedTrackDca(cfgTrackCuts.binsPt, cfgTrackCuts.limitsDca, trackParCovPos.getPt(), dcaPos.getY(), dcaPos.getZ())) {
          continue;
        }

        int pidTrackPosKaon{-1};
        int pidTrackPosPion{-1};
        if (cfgTrackCuts.usePidTpcOnly) {
          /// kaon TPC PID positive daughter
          pidTrackPosKaon = selectorKaon.statusTpc(trackPos);
          /// pion TPC PID positive daughter
          pidTrackPosPion = selectorPion.statusTpc(trackPos);
        } else {
          /// kaon TPC, TOF PID positive daughter
          pidTrackPosKaon = selectorKaon.statusTpcAndTof(trackPos);
          /// pion TPC, TOF PID positive daughter
          pidTrackPosPion = selectorPion.statusTpcAndTof(trackPos);
        }

        for (auto const& trackIndexNeg : groupedTrackIndices) {
          auto trackNeg = trackIndexNeg.template track_as<TTracks>();
          // track selections
          if (trackNeg.sign() > 0) { // second negative track
            continue;
          }
          if (!trackNeg.isGlobalTrackWoDCA()) {
            continue;
          }
          if (trackNeg.pt() < cfgTrackCuts.ptMin) {
            continue;
          }
          if (std::abs(trackNeg.eta()) > cfgTrackCuts.absEtaMax) {
            continue;
          }
          auto trackParCovNeg = getTrackParCov(trackNeg);
          o2::dataformats::DCA dcaNeg;
          trackParCovNeg.propagateToDCA(primaryVertex, bz, &dcaNeg);
          if (!isSelectedTrackDca(cfgTrackCuts.binsPt, cfgTrackCuts.limitsDca, trackParCovNeg.getPt(), dcaNeg.getY(), dcaNeg.getZ())) {
            continue;
          }

          int pidTrackNegKaon{-1};
          int pidTrackNegPion{-1};
          if (cfgTrackCuts.usePidTpcOnly) {
            /// kaon TPC PID negative daughter
            pidTrackNegKaon = selectorKaon.statusTpc(trackNeg);
            /// pion TPC PID negative daughter
            pidTrackNegPion = selectorPion.statusTpc(trackNeg);
          } else {
            /// kaon TPC, TOF PID negative daughter
            pidTrackNegKaon = selectorKaon.statusTpcAndTof(trackNeg);
            /// pion TPC, TOF PID negative daughter
            pidTrackNegPion = selectorPion.statusTpcAndTof(trackNeg);
          }

          // preselections
          // PID
          uint8_t massHypo{D0MassHypo::D0AndD0Bar}; // both mass hypotheses a priori
          if (pidTrackPosPion == TrackSelectorPID::Rejected || pidTrackNegKaon == TrackSelectorPID::Rejected) {
            massHypo -= D0MassHypo::D0; // exclude D0
          }
          if (pidTrackNegPion == TrackSelectorPID::Rejected || pidTrackPosKaon == TrackSelectorPID::Rejected) {
            massHypo -= D0MassHypo::D0Bar; // exclude D0Bar
          }
          if (massHypo == 0) {
            continue;
          }

          // pt
          std::array<float, 3> pVecNoVtxD0 = RecoDecay::pVec(trackPos.pVector(), trackNeg.pVector());
          float ptNoVtxD0 = RecoDecay::pt(pVecNoVtxD0);
          if (ptNoVtxD0 - ptTolerance < cfgCandCuts.ptMin) {
            continue;
          }
          int ptBinNoVtxD0 = findBin(cfgTrackCuts.binsPt, ptNoVtxD0 + ptTolerance); // assuming tighter selections at lower pT
          if (ptBinNoVtxD0 < 0) {
            continue;
          }

          // d0xd0
          if (dcaPos.getY() * dcaNeg.getY() > cfgCandCuts.topologicalCuts->get(ptBinNoVtxD0, "max d0d0")) {
            continue;
          }

          // invariant mass
          if (massHypo == D0MassHypo::D0 || massHypo == D0MassHypo::D0AndD0Bar) {
            float invMassNoVtxD0 = RecoDecay::m(std::array{trackPos.pVector(), trackNeg.pVector()}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
            if (std::abs(invMassNoVtxD0 - o2::constants::physics::MassD0) > cfgCandCuts.topologicalCuts->get(ptBinNoVtxD0, "delta inv. mass") + invMassTolerance) {
              massHypo -= D0MassHypo::D0;
            }
          }
          if (massHypo >= D0MassHypo::D0Bar) {
            float invMassNoVtxD0bar = RecoDecay::m(std::array{trackNeg.pVector(), trackPos.pVector()}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
            if (std::abs(invMassNoVtxD0bar - o2::constants::physics::MassD0) > cfgCandCuts.topologicalCuts->get(ptBinNoVtxD0, "delta inv. mass") + invMassTolerance) {
              massHypo -= D0MassHypo::D0Bar;
            }
          }
          if (massHypo == 0) {
            continue;
          }

          // reconstruct vertex
          if (df.process(trackParCovPos, trackParCovNeg) == 0) {
            continue;
          }
          const auto& secondaryVertex = df.getPCACandidate();
          auto chi2PCA = df.getChi2AtPCACandidate();
          auto covMatrixPCA = df.calcPCACovMatrixFlat();
          auto trackParAtSecVtxPos = df.getTrack(0);
          auto trackParAtSecVtxNeg = df.getTrack(1);

          std::array<float, 3> pVecPos{};
          std::array<float, 3> pVecNeg{};
          trackParAtSecVtxPos.getPxPyPzGlo(pVecPos);
          trackParAtSecVtxNeg.getPxPyPzGlo(pVecNeg);
          std::array<float, 3> pVecD0 = RecoDecay::pVec(pVecPos, pVecNeg);

          // select D0
          // pt
          float ptD0 = RecoDecay::pt(pVecD0);
          if (ptD0 < cfgCandCuts.ptMin) {
            continue;
          }
          int ptBinD0 = findBin(cfgTrackCuts.binsPt, ptD0);
          if (ptBinD0 < 0) {
            continue;
          }

          // random downsampling already here
          if (cfgDownsampling.apply) {
            int ptBinWeights{0};
            if (ptD0 < histDownSampl->GetBinLowEdge(1)) {
              ptBinWeights = 1;
            } else if (ptD0 > histDownSampl->GetXaxis()->GetBinUpEdge(histDownSampl->GetNbinsX())) {
              ptBinWeights = histDownSampl->GetNbinsX();
            } else {
              ptBinWeights = histDownSampl->GetXaxis()->FindBin(ptD0);
            }
            float weight = histDownSampl->GetBinContent(ptBinWeights);
            if (gRandom->Rndm() > weight) {
              continue;
            }
          }

          // d0xd0
          if (dcaPos.getY() * dcaNeg.getY() > cfgCandCuts.topologicalCuts->get(ptBinD0, "max d0d0")) {
            continue;
          }
          // cospa
          float cosPaD0 = RecoDecay::cpa(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, pVecD0);
          if (cosPaD0 < cfgCandCuts.topologicalCuts->get(ptBinD0, "min cos pointing angle")) {
            continue;
          }
          // pointing angle
          float paD0 = std::acos(cosPaD0);
          if (paD0 > cfgCandCuts.topologicalCuts->get(ptBinD0, "max pointing angle")) {
            continue;
          }
          // cospa XY
          float cosPaXYD0 = RecoDecay::cpaXY(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, pVecD0);
          if (cosPaXYD0 < cfgCandCuts.topologicalCuts->get(ptBinD0, "min cos pointing angle XY")) {
            continue;
          }
          // pointing angle XY
          float paXYD0 = std::acos(cosPaXYD0);
          if (paXYD0 > cfgCandCuts.topologicalCuts->get(ptBinD0, "max pointing angle XY")) {
            continue;
          }
          // decay length
          float decLenD0 = RecoDecay::distance(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex);
          if (decLenD0 < cfgCandCuts.topologicalCuts->get(ptBinD0, "min decay length")) {
            continue;
          }
          // decay length XY
          float decLenXYD0 = RecoDecay::distanceXY(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex);
          if (decLenXYD0 < cfgCandCuts.topologicalCuts->get(ptBinD0, "min decay length XY")) {
            continue;
          }
          // normalised decay length
          float phi{0.f}, theta{0.f};
          getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
          float errorDecayLengthD0 = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          if (decLenD0 / errorDecayLengthD0 < cfgCandCuts.topologicalCuts->get(ptBinD0, "min norm decay length")) {
            continue;
          }
          // normalised decay length XY
          float errorDecayLengthXYD0 = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.f) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.f));
          if (decLenXYD0 / errorDecayLengthXYD0 < cfgCandCuts.topologicalCuts->get(ptBinD0, "min norm decay length XY")) {
            continue;
          }

          float invMassD0{0.f}, invMassD0bar{0.f};
          std::vector<float> bdtScoresD0{0.f, 1.f, 1.f}, bdtScoresD0bar{0.f, 1.f, 1.f}; // always selected a priori
          if (massHypo == D0MassHypo::D0 || massHypo == D0MassHypo::D0AndD0Bar) {
            invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
            if (std::abs(invMassD0 - o2::constants::physics::MassD0) > cfgCandCuts.topologicalCuts->get(ptBinD0, "delta inv. mass")) {
              massHypo -= D0MassHypo::D0;
              bdtScoresD0 = std::vector<float>{1.f, 0.f, 0.f};
            } else {
              // apply BDT models
              if (cfgMl.apply) {
                std::vector<float> featuresCandD0 = {dcaPos.getY(), dcaNeg.getY(), chi2PCA, cosPaD0, cosPaXYD0, decLenXYD0, decLenD0, dcaPos.getY() * dcaNeg.getY(), aod::pid_tpc_tof_utils::combineNSigma<false>(trackPos.tpcNSigmaPi(), trackPos.tofNSigmaPi()), aod::pid_tpc_tof_utils::combineNSigma<false>(trackNeg.tpcNSigmaKa(), trackNeg.tofNSigmaKa()), trackPos.tpcNSigmaPi(), trackPos.tpcNSigmaKa(), aod::pid_tpc_tof_utils::combineNSigma<false>(trackPos.tpcNSigmaKa(), trackPos.tofNSigmaKa()), trackNeg.tpcNSigmaPi(), trackNeg.tpcNSigmaKa(), aod::pid_tpc_tof_utils::combineNSigma<false>(trackNeg.tpcNSigmaPi(), trackNeg.tofNSigmaPi())};
                if (!mlResponse.isSelectedMl(featuresCandD0, ptD0, bdtScoresD0)) {
                  massHypo -= D0MassHypo::D0;
                }
              }
            }
          }
          if (massHypo >= D0MassHypo::D0Bar) {
            invMassD0bar = RecoDecay::m(std::array{pVecNeg, pVecPos}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
            if (std::abs(invMassD0bar - o2::constants::physics::MassD0) > cfgCandCuts.topologicalCuts->get(ptBinD0, "delta inv. mass")) {
              massHypo -= D0MassHypo::D0Bar;
              bdtScoresD0bar = std::vector<float>{1.f, 0.f, 0.f};
            } else {
              // apply BDT models
              if (cfgMl.apply) {
                std::vector<float> featuresCandD0bar = {dcaPos.getY(), dcaNeg.getY(), chi2PCA, cosPaD0, cosPaXYD0, decLenXYD0, decLenD0, dcaPos.getY() * dcaNeg.getY(), aod::pid_tpc_tof_utils::combineNSigma<false>(trackNeg.tpcNSigmaPi(), trackNeg.tofNSigmaPi()), aod::pid_tpc_tof_utils::combineNSigma<false>(trackPos.tpcNSigmaKa(), trackPos.tofNSigmaKa()), trackNeg.tpcNSigmaPi(), trackNeg.tpcNSigmaKa(), aod::pid_tpc_tof_utils::combineNSigma<false>(trackNeg.tpcNSigmaKa(), trackNeg.tofNSigmaKa()), trackPos.tpcNSigmaPi(), trackPos.tpcNSigmaKa(), aod::pid_tpc_tof_utils::combineNSigma<false>(trackPos.tpcNSigmaPi(), trackPos.tofNSigmaPi())};
                if (!mlResponse.isSelectedMl(featuresCandD0bar, ptD0, bdtScoresD0bar)) {
                  massHypo -= D0MassHypo::D0Bar;
                }
              }
            }
          }
          if (massHypo == 0) {
            continue;
          }

          float etaD0 = RecoDecay::eta(pVecD0);
          float phiD0 = RecoDecay::phi(pVecD0);

          // fill tables
          // collision
          if (!selectedCollisions.count(collision.globalIndex())) {
            // fill collision table if not yet present
            collTable(collision.posX(),
                      collision.posY(),
                      collision.posZ(),
                      collision.covXX(),
                      collision.covXY(),
                      collision.covXZ(),
                      collision.covYY(),
                      collision.covYZ(),
                      collision.covZZ(),
                      collision.numContrib(),
                      uint8_t(std::round(collision.centFT0C())),
                      getCompressedOccupancy(collision.trackOccupancyInTimeRange()),
                      getCompressedOccupancy(collision.ft0cOccupancyInTimeRange()),
                      eventIR.orbit,
                      runNumber);
            selectedCollisions[collision.globalIndex()] = collTable.lastIndex();
          }
          // tracks
          if (!selectedTracks.count(trackPos.globalIndex())) {
            // fill track table with positive track if not yet present
            uint8_t tmoPrimUnfm80{0u};
            uint8_t tmoFV0AUnfm80{0u};
            uint8_t tmoFT0AUnfm80{0u};
            uint8_t tmoFT0CUnfm80{0u};
            uint8_t twmoPrimUnfm80{0u};
            uint8_t twmoFV0AUnfm80{0u};
            uint8_t twmoFT0AUnfm80{0u};
            uint8_t twmoFT0CUnfm80{0u};
            uint8_t tmoRobustT0V0PrimUnfm80{0u};
            uint8_t twmoRobustT0V0PrimUnfm80{0u};
            if (trackPos.has_tmo()) {
              auto tmoFromTrack = trackPos.template tmo_as<TrackMeanOccs>(); // obtain track mean occupancies
              tmoPrimUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoPrimUnfm80());
              tmoFV0AUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoFV0AUnfm80());
              tmoFT0AUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoFT0AUnfm80());
              tmoFT0CUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoFT0CUnfm80());
              twmoPrimUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoPrimUnfm80());
              twmoFV0AUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoFV0AUnfm80());
              twmoFT0AUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoFT0AUnfm80());
              twmoFT0CUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoFT0CUnfm80());
              tmoRobustT0V0PrimUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoRobustT0V0PrimUnfm80());
              twmoRobustT0V0PrimUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoRobustT0V0PrimUnfm80());
            }
            float tpcTime0{0.f};
            float tpcdEdxNorm{0.f};
            int16_t tpcDcaR{0};
            int16_t tpcDcaZ{0};
            uint8_t tpcClusterByteMask{0u};
            uint8_t tpcdEdxMax0R{0u};
            uint8_t tpcdEdxMax1R{0u};
            uint8_t tpcdEdxMax2R{0u};
            uint8_t tpcdEdxMax3R{0u};
            uint8_t tpcdEdxTot0R{0u};
            uint8_t tpcdEdxTot1R{0u};
            uint8_t tpcdEdxTot2R{0u};
            uint8_t tpcdEdxTot3R{0u};
            int8_t deltaRefContParamY{0};
            int8_t deltaRefITSParamZ{0};
            int8_t deltaRefContParamSnp{0};
            int8_t deltaRefContParamTgl{0};
            int8_t deltaRefContParamQ2Pt{0};
            int8_t deltaRefGloParamY{0};
            int8_t deltaRefGloParamZ{0};
            int8_t deltaRefGloParamSnp{0};
            int8_t deltaRefGloParamTgl{0};
            int8_t deltaRefGloParamQ2Pt{0};
            int8_t deltaTOFdX{0};
            int8_t deltaTOFdZ{0};
            if constexpr (withTrackQa) {
              if (trackPos.has_trackQA()) {
                auto trackQA = trackPos.template trackQA_as<TTrackQa>(); // obtain track QA
                tpcTime0 = trackQA.tpcTime0();
                tpcdEdxNorm = trackQA.tpcdEdxNorm();
                tpcDcaR = trackQA.tpcdcaR();
                tpcDcaZ = trackQA.tpcdcaZ();
                tpcClusterByteMask = trackQA.tpcClusterByteMask();
                tpcdEdxMax0R = trackQA.tpcdEdxMax0R();
                tpcdEdxMax1R = trackQA.tpcdEdxMax1R();
                tpcdEdxMax2R = trackQA.tpcdEdxMax2R();
                tpcdEdxMax3R = trackQA.tpcdEdxMax3R();
                tpcdEdxTot0R = trackQA.tpcdEdxTot0R();
                tpcdEdxTot1R = trackQA.tpcdEdxTot1R();
                tpcdEdxTot2R = trackQA.tpcdEdxTot2R();
                tpcdEdxTot3R = trackQA.tpcdEdxTot3R();
                deltaRefContParamY = trackQA.deltaRefContParamY();
                deltaRefITSParamZ = trackQA.deltaRefITSParamZ();
                deltaRefContParamSnp = trackQA.deltaRefContParamSnp();
                deltaRefContParamTgl = trackQA.deltaRefContParamTgl();
                deltaRefContParamQ2Pt = trackQA.deltaRefContParamQ2Pt();
                deltaRefGloParamY = trackQA.deltaRefGloParamY();
                deltaRefGloParamZ = trackQA.deltaRefGloParamZ();
                deltaRefGloParamSnp = trackQA.deltaRefGloParamSnp();
                deltaRefGloParamTgl = trackQA.deltaRefGloParamTgl();
                deltaRefGloParamQ2Pt = trackQA.deltaRefGloParamQ2Pt();
                deltaTOFdX = trackQA.deltaTOFdX();
                deltaTOFdZ = trackQA.deltaTOFdZ();
              }
            }

            trackTable(selectedCollisions[collision.globalIndex()], // stored at PV
                       trackPos.x(),
                       trackPos.alpha(),
                       trackPos.y(),
                       trackPos.z(),
                       trackPos.snp(),
                       trackPos.tgl(),
                       trackPos.signed1Pt(),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cYY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cZY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cZZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cSnpY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cSnpZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cSnpSnp(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cTglY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cTglZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cTglSnp(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.cTglTgl(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.c1PtY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.c1PtZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.c1PtSnp(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.c1PtTgl(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackPos.c1Pt21Pt2(), precisionCovMask),
                       trackPos.tpcInnerParam(),
                       trackPos.flags(),
                       trackPos.itsClusterSizes(),
                       trackPos.tpcNClsFindable(),
                       trackPos.tpcNClsFindableMinusFound(),
                       trackPos.tpcNClsFindableMinusCrossedRows(),
                       trackPos.tpcNClsShared(),
                       trackPos.trdPattern(),
                       getCompressedChi2(trackPos.itsChi2NCl()),
                       getCompressedChi2(trackPos.tpcChi2NCl()),
                       getCompressedChi2(trackPos.trdChi2()),
                       getCompressedChi2(trackPos.tofChi2()),
                       trackPos.tpcSignal(),
                       trackPos.trdSignal(),
                       trackPos.length(),
                       trackPos.tofExpMom(),
                       trackPos.trackTime(),
                       trackPos.trackTimeRes(),
                       tpcTime0,
                       tpcdEdxNorm,
                       tpcDcaR,
                       tpcDcaZ,
                       tpcClusterByteMask,
                       tpcdEdxMax0R,
                       tpcdEdxMax1R,
                       tpcdEdxMax2R,
                       tpcdEdxMax3R,
                       tpcdEdxTot0R,
                       tpcdEdxTot1R,
                       tpcdEdxTot2R,
                       tpcdEdxTot3R,
                       deltaRefContParamY,
                       deltaRefITSParamZ,
                       deltaRefContParamSnp,
                       deltaRefContParamTgl,
                       deltaRefContParamQ2Pt,
                       deltaRefGloParamY,
                       deltaRefGloParamZ,
                       deltaRefGloParamSnp,
                       deltaRefGloParamTgl,
                       deltaRefGloParamQ2Pt,
                       deltaTOFdX,
                       deltaTOFdZ,
                       o2::math_utils::detail::truncateFloatFraction(dcaPos.getY(), precisionDcaMask),
                       o2::math_utils::detail::truncateFloatFraction(dcaPos.getZ(), precisionDcaMask),
                       getCompressedNumSigmaPid(trackPos.tpcNSigmaPi()),
                       getCompressedNumSigmaPid(trackPos.tpcNSigmaKa()),
                       getCompressedNumSigmaPid(trackPos.tofNSigmaPi()),
                       getCompressedNumSigmaPid(trackPos.tofNSigmaKa()),
                       tmoPrimUnfm80,
                       tmoFV0AUnfm80,
                       tmoFT0AUnfm80,
                       tmoFT0CUnfm80,
                       twmoPrimUnfm80,
                       twmoFV0AUnfm80,
                       twmoFT0AUnfm80,
                       twmoFT0CUnfm80,
                       tmoRobustT0V0PrimUnfm80,
                       twmoRobustT0V0PrimUnfm80);
            selectedTracks[trackPos.globalIndex()] = trackTable.lastIndex();
          }
          if (!selectedTracks.count(trackNeg.globalIndex())) {
            // fill track table with negative track if not yet present
            uint8_t tmoPrimUnfm80{0u};
            uint8_t tmoFV0AUnfm80{0u};
            uint8_t tmoFT0AUnfm80{0u};
            uint8_t tmoFT0CUnfm80{0u};
            uint8_t twmoPrimUnfm80{0u};
            uint8_t twmoFV0AUnfm80{0u};
            uint8_t twmoFT0AUnfm80{0u};
            uint8_t twmoFT0CUnfm80{0u};
            uint8_t tmoRobustT0V0PrimUnfm80{0u};
            uint8_t twmoRobustT0V0PrimUnfm80{0u};
            if (trackNeg.has_tmo()) {
              auto tmoFromTrack = trackNeg.template tmo_as<TrackMeanOccs>(); // obtain track mean occupancies
              tmoPrimUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoPrimUnfm80());
              tmoFV0AUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoFV0AUnfm80());
              tmoFT0AUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoFT0AUnfm80());
              tmoFT0CUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoFT0CUnfm80());
              twmoPrimUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoPrimUnfm80());
              twmoFV0AUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoFV0AUnfm80());
              twmoFT0AUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoFT0AUnfm80());
              twmoFT0CUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoFT0CUnfm80());
              tmoRobustT0V0PrimUnfm80 = getCompressedOccupancy(tmoFromTrack.tmoRobustT0V0PrimUnfm80());
              twmoRobustT0V0PrimUnfm80 = getCompressedOccupancy(tmoFromTrack.twmoRobustT0V0PrimUnfm80());
            }
            float tpcTime0{0.f};
            float tpcdEdxNorm{0.f};
            int16_t tpcDcaR{0};
            int16_t tpcDcaZ{0};
            uint8_t tpcClusterByteMask{0u};
            uint8_t tpcdEdxMax0R{0u};
            uint8_t tpcdEdxMax1R{0u};
            uint8_t tpcdEdxMax2R{0u};
            uint8_t tpcdEdxMax3R{0u};
            uint8_t tpcdEdxTot0R{0u};
            uint8_t tpcdEdxTot1R{0u};
            uint8_t tpcdEdxTot2R{0u};
            uint8_t tpcdEdxTot3R{0u};
            int8_t deltaRefContParamY{0};
            int8_t deltaRefITSParamZ{0};
            int8_t deltaRefContParamSnp{0};
            int8_t deltaRefContParamTgl{0};
            int8_t deltaRefContParamQ2Pt{0};
            int8_t deltaRefGloParamY{0};
            int8_t deltaRefGloParamZ{0};
            int8_t deltaRefGloParamSnp{0};
            int8_t deltaRefGloParamTgl{0};
            int8_t deltaRefGloParamQ2Pt{0};
            int8_t deltaTOFdX{0};
            int8_t deltaTOFdZ{0};
            if constexpr (withTrackQa) {
              if (trackNeg.has_trackQA()) {
                auto trackQA = trackNeg.template trackQA_as<TTrackQa>(); // obtain track QA
                tpcTime0 = trackQA.tpcTime0();
                tpcdEdxNorm = trackQA.tpcdEdxNorm();
                tpcDcaR = trackQA.tpcdcaR();
                tpcDcaZ = trackQA.tpcdcaZ();
                tpcClusterByteMask = trackQA.tpcClusterByteMask();
                tpcdEdxMax0R = trackQA.tpcdEdxMax0R();
                tpcdEdxMax1R = trackQA.tpcdEdxMax1R();
                tpcdEdxMax2R = trackQA.tpcdEdxMax2R();
                tpcdEdxMax3R = trackQA.tpcdEdxMax3R();
                tpcdEdxTot0R = trackQA.tpcdEdxTot0R();
                tpcdEdxTot1R = trackQA.tpcdEdxTot1R();
                tpcdEdxTot2R = trackQA.tpcdEdxTot2R();
                tpcdEdxTot3R = trackQA.tpcdEdxTot3R();
                deltaRefContParamY = trackQA.deltaRefContParamY();
                deltaRefITSParamZ = trackQA.deltaRefITSParamZ();
                deltaRefContParamSnp = trackQA.deltaRefContParamSnp();
                deltaRefContParamTgl = trackQA.deltaRefContParamTgl();
                deltaRefContParamQ2Pt = trackQA.deltaRefContParamQ2Pt();
                deltaRefGloParamY = trackQA.deltaRefGloParamY();
                deltaRefGloParamZ = trackQA.deltaRefGloParamZ();
                deltaRefGloParamSnp = trackQA.deltaRefGloParamSnp();
                deltaRefGloParamTgl = trackQA.deltaRefGloParamTgl();
                deltaRefGloParamQ2Pt = trackQA.deltaRefGloParamQ2Pt();
                deltaTOFdX = trackQA.deltaTOFdX();
                deltaTOFdZ = trackQA.deltaTOFdZ();
              }
            }

            trackTable(selectedCollisions[collision.globalIndex()], // stored at PV
                       trackNeg.x(),
                       trackNeg.alpha(),
                       trackNeg.y(),
                       trackNeg.z(),
                       trackNeg.snp(),
                       trackNeg.tgl(),
                       trackNeg.signed1Pt(),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cYY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cZY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cZZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cSnpY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cSnpZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cSnpSnp(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cTglY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cTglZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cTglSnp(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.cTglTgl(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.c1PtY(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.c1PtZ(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.c1PtSnp(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.c1PtTgl(), precisionCovMask),
                       o2::math_utils::detail::truncateFloatFraction(trackNeg.c1Pt21Pt2(), precisionCovMask),
                       trackNeg.tpcInnerParam(),
                       trackNeg.flags(),
                       trackNeg.itsClusterSizes(),
                       trackNeg.tpcNClsFindable(),
                       trackNeg.tpcNClsFindableMinusFound(),
                       trackNeg.tpcNClsFindableMinusCrossedRows(),
                       trackNeg.tpcNClsShared(),
                       trackNeg.trdPattern(),
                       getCompressedChi2(trackNeg.itsChi2NCl()),
                       getCompressedChi2(trackNeg.tpcChi2NCl()),
                       getCompressedChi2(trackNeg.trdChi2()),
                       getCompressedChi2(trackNeg.tofChi2()),
                       trackNeg.tpcSignal(),
                       trackNeg.trdSignal(),
                       trackNeg.length(),
                       trackNeg.tofExpMom(),
                       trackNeg.trackTime(),
                       trackNeg.trackTimeRes(),
                       tpcTime0,
                       tpcdEdxNorm,
                       tpcDcaR,
                       tpcDcaZ,
                       tpcClusterByteMask,
                       tpcdEdxMax0R,
                       tpcdEdxMax1R,
                       tpcdEdxMax2R,
                       tpcdEdxMax3R,
                       tpcdEdxTot0R,
                       tpcdEdxTot1R,
                       tpcdEdxTot2R,
                       tpcdEdxTot3R,
                       deltaRefContParamY,
                       deltaRefITSParamZ,
                       deltaRefContParamSnp,
                       deltaRefContParamTgl,
                       deltaRefContParamQ2Pt,
                       deltaRefGloParamY,
                       deltaRefGloParamZ,
                       deltaRefGloParamSnp,
                       deltaRefGloParamTgl,
                       deltaRefGloParamQ2Pt,
                       deltaTOFdX,
                       deltaTOFdZ,
                       o2::math_utils::detail::truncateFloatFraction(dcaNeg.getY(), precisionDcaMask),
                       o2::math_utils::detail::truncateFloatFraction(dcaNeg.getZ(), precisionDcaMask),
                       getCompressedNumSigmaPid(trackNeg.tpcNSigmaPi()),
                       getCompressedNumSigmaPid(trackNeg.tpcNSigmaKa()),
                       getCompressedNumSigmaPid(trackNeg.tofNSigmaPi()),
                       getCompressedNumSigmaPid(trackNeg.tofNSigmaKa()),
                       tmoPrimUnfm80,
                       tmoFV0AUnfm80,
                       tmoFT0AUnfm80,
                       tmoFT0CUnfm80,
                       twmoPrimUnfm80,
                       twmoFV0AUnfm80,
                       twmoFT0AUnfm80,
                       twmoFT0CUnfm80,
                       tmoRobustT0V0PrimUnfm80,
                       twmoRobustT0V0PrimUnfm80);
            selectedTracks[trackNeg.globalIndex()] = trackTable.lastIndex();
          }

          // candidate
          candTable(selectedCollisions[collision.globalIndex()],
                    selectedTracks[trackPos.globalIndex()],
                    selectedTracks[trackNeg.globalIndex()],
                    massHypo,
                    ptD0,
                    etaD0,
                    phiD0,
                    invMassD0,
                    invMassD0bar,
                    getCompressedDecayLength(decLenD0),
                    getCompressedDecayLength(decLenXYD0),
                    getCompressedNormDecayLength(decLenD0 / errorDecayLengthD0),
                    getCompressedNormDecayLength(decLenXYD0 / errorDecayLengthXYD0),
                    getCompressedCosPa(cosPaD0),
                    getCompressedCosPa(cosPaXYD0),
                    getCompressedPointingAngle(paD0),
                    getCompressedPointingAngle(paXYD0),
                    getCompressedChi2(chi2PCA),
                    getCompressedBdtScoreBkg(bdtScoresD0[0]),
                    getCompressedBdtScoreSgn(bdtScoresD0[1]),
                    getCompressedBdtScoreSgn(bdtScoresD0[2]),
                    getCompressedBdtScoreBkg(bdtScoresD0bar[0]),
                    getCompressedBdtScoreSgn(bdtScoresD0bar[1]),
                    getCompressedBdtScoreSgn(bdtScoresD0bar[2]));
        } // end loop over negative tracks
      } // end loop over positive tracks
    } // end loop over collisions tracks
  }

  // process functions
  void processWithTrackQa(CollisionsWEvSel const& collisions,
                          aod::TrackAssoc const& trackIndices,
                          TracksWCovExtraPidAndQa const& tracks,
                          aod::BCsWithTimestamps const& bcs,
                          TrackMeanOccs const& occ,
                          aod::TracksQAVersion const& trackQa)
  {
    runDataCreation<true>(collisions, trackIndices, tracks, bcs, occ, trackQa);
  }
  PROCESS_SWITCH(DerivedDataCreatorD0Calibration, processWithTrackQa, "Process with trackQA enabled", false);

  void processNoTrackQa(CollisionsWEvSel const& collisions,
                        aod::TrackAssoc const& trackIndices,
                        TracksWCovExtraPid const& tracks,
                        aod::BCsWithTimestamps const& bcs,
                        TrackMeanOccs const& occ)
  {
    runDataCreation<false>(collisions, trackIndices, tracks, bcs, occ, nullptr);
  }
  PROCESS_SWITCH(DerivedDataCreatorD0Calibration, processNoTrackQa, "Process without trackQA enabled", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DerivedDataCreatorD0Calibration>(cfgc)};
}
