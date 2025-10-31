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

/// \file treeCreatorTccToD0D0Pi.cxx
/// \brief tree creator for studying the charm exotic state Tcc to D0D0pi
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/V0.h>

#include <TH1.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(PxProng0D1, pxProng0D1, float);
DECLARE_SOA_COLUMN(PxProng1D1, pxProng1D1, float);
DECLARE_SOA_COLUMN(PyProng0D1, pyProng0D1, float);
DECLARE_SOA_COLUMN(PyProng1D1, pyProng1D1, float);
DECLARE_SOA_COLUMN(PzProng0D1, pzProng0D1, float);
DECLARE_SOA_COLUMN(PzProng1D1, pzProng1D1, float);
DECLARE_SOA_COLUMN(PxProng0D2, pxProng0D2, float);
DECLARE_SOA_COLUMN(PxProng1D2, pxProng1D2, float);
DECLARE_SOA_COLUMN(PyProng0D2, pyProng0D2, float);
DECLARE_SOA_COLUMN(PyProng1D2, pyProng1D2, float);
DECLARE_SOA_COLUMN(PzProng0D2, pzProng0D2, float);
DECLARE_SOA_COLUMN(PzProng1D2, pzProng1D2, float);
DECLARE_SOA_COLUMN(PxSoftPi, pxSoftPi, float);
DECLARE_SOA_COLUMN(PySoftPi, pySoftPi, float);
DECLARE_SOA_COLUMN(PzSoftPi, pzSoftPi, float);
DECLARE_SOA_COLUMN(SelFlagD1, selFlagD1, int8_t);
DECLARE_SOA_COLUMN(SelFlagD2, selFlagD2, int8_t);
DECLARE_SOA_COLUMN(MD1, mD1, float);
DECLARE_SOA_COLUMN(MD2, mD2, float);
DECLARE_SOA_COLUMN(DeltaMD1, deltaMD1, float);
DECLARE_SOA_COLUMN(DeltaMD2, deltaMD2, float);
DECLARE_SOA_COLUMN(MDD, mDD, float);
DECLARE_SOA_COLUMN(MDPi1, mDPi1, float);
DECLARE_SOA_COLUMN(MDPi2, mDPi2, float);
DECLARE_SOA_COLUMN(MDDPi, mDDPi, float);
DECLARE_SOA_COLUMN(DeltaMDDPi, deltaMDDPi, float);
DECLARE_SOA_COLUMN(EtaD1, etaD1, float);
DECLARE_SOA_COLUMN(EtaD2, etaD2, float);
DECLARE_SOA_COLUMN(EtaSoftPi, etaSoftPi, float);
DECLARE_SOA_COLUMN(PhiD1, phiD1, float);
DECLARE_SOA_COLUMN(PhiD2, phiD2, float);
DECLARE_SOA_COLUMN(PhiSoftPi, phiSoftPi, float);
DECLARE_SOA_COLUMN(YD1, yD1, float);
DECLARE_SOA_COLUMN(YD2, yD2, float);
DECLARE_SOA_COLUMN(YSoftPi, ySoftPi, float);
DECLARE_SOA_COLUMN(NSigTpcSoftPi, nSigTpcSoftPi, float);
DECLARE_SOA_COLUMN(NSigTofSoftPi, nSigTofSoftPi, float);
DECLARE_SOA_COLUMN(MlScoreD1, mlScoreD1, float);
DECLARE_SOA_COLUMN(MlScoreD2, mlScoreD2, float);
DECLARE_SOA_COLUMN(ImpactParameterD1, impactParameterD1, float);
DECLARE_SOA_COLUMN(ImpactParameterD2, impactParameterD2, float);
DECLARE_SOA_COLUMN(ImpactParameterSoftPi, impactParameterSoftPi, float);
DECLARE_SOA_COLUMN(CpaD1, cpaD1, float);
DECLARE_SOA_COLUMN(CpaD2, cpaD2, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);
DECLARE_SOA_COLUMN(SignSoftPi, signSoftPi, float);
DECLARE_SOA_COLUMN(DcaXYSoftPi, dcaXYSoftPi, float);
DECLARE_SOA_COLUMN(DcaZSoftPi, dcaZSoftPi, float);
DECLARE_SOA_COLUMN(NITSClsSoftPi, nITSClsSoftPi, float);
DECLARE_SOA_COLUMN(NTPCClsCrossedRowsSoftPi, nTPCClsCrossedRowsSoftPi, float);
DECLARE_SOA_COLUMN(NTPCChi2NClSoftPi, nTPCChi2NClSoftPi, float);
DECLARE_SOA_COLUMN(CentOfCand, centOfCand, float);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);     //! Global index for the collisionAdd commentMore actions
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t); //! Timestamp for the collision
} // namespace full

DECLARE_SOA_TABLE(HfCandTccLites, "AOD", "HFCANDTCCLITE",
                  full::PxProng0D1,
                  full::PxProng1D1,
                  full::PyProng0D1,
                  full::PyProng1D1,
                  full::PzProng0D1,
                  full::PzProng1D1,
                  full::PxProng0D2,
                  full::PxProng1D2,
                  full::PyProng0D2,
                  full::PyProng1D2,
                  full::PzProng0D2,
                  full::PzProng1D2,
                  full::PxSoftPi,
                  full::PySoftPi,
                  full::PzSoftPi,
                  full::SelFlagD1,
                  full::SelFlagD2,
                  full::MD1,
                  full::MD2,
                  full::DeltaMD1,
                  full::DeltaMD2,
                  full::MDPi1,
                  full::MDPi2,
                  full::MDDPi,
                  full::DeltaMDDPi,
                  full::EtaD1,
                  full::EtaD2,
                  full::EtaSoftPi,
                  full::PhiD1,
                  full::PhiD2,
                  full::PhiSoftPi,
                  full::YD1,
                  full::YD2,
                  full::YSoftPi,
                  full::NSigTpcSoftPi,
                  full::NSigTofSoftPi,
                  full::MlScoreD1,
                  full::MlScoreD2,
                  full::ImpactParameterD1,
                  full::ImpactParameterD2,
                  full::ImpactParameterSoftPi,
                  full::CpaD1,
                  full::CpaD2,
                  full::Chi2PCA,
                  full::SignSoftPi,
                  full::DcaXYSoftPi,
                  full::DcaZSoftPi,
                  full::NITSClsSoftPi,
                  full::NTPCClsCrossedRowsSoftPi,
                  full::NTPCChi2NClSoftPi,
                  full::CentOfCand,
                  full::GIndexCol,
                  full::TimeStamp);

DECLARE_SOA_TABLE(HfCandDDPairs, "AOD", "HFCANDDDPAIR",
                  full::PxProng0D1,
                  full::PxProng1D1,
                  full::PyProng0D1,
                  full::PyProng1D1,
                  full::PzProng0D1,
                  full::PzProng1D1,
                  full::PxProng0D2,
                  full::PxProng1D2,
                  full::PyProng0D2,
                  full::PyProng1D2,
                  full::PzProng0D2,
                  full::PzProng1D2,
                  full::SelFlagD1,
                  full::SelFlagD2,
                  full::EtaD1,
                  full::EtaD2,
                  full::PhiD1,
                  full::PhiD2,
                  full::MlScoreD1,
                  full::MlScoreD2,
                  full::CentOfCand,
                  full::GIndexCol,
                  full::TimeStamp);

DECLARE_SOA_TABLE(HfCandTccFullEvs, "AOD", "HFCANDTCCFULLEV",
                  full::CollisionId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorTccToD0D0Pi {
  Produces<o2::aod::HfCandTccLites> rowCandidateLite;
  Produces<o2::aod::HfCandDDPairs> rowCandidateDDPair;
  Produces<o2::aod::HfCandTccFullEvs> rowCandidateFullEvents;

  Configurable<float> ptMinSoftPion{"ptMinSoftPion", 0.0, "Min pt for the soft pion"};
  Configurable<bool> usePionIsGlobalTrackWoDCA{"usePionIsGlobalTrackWoDCA", true, "check isGlobalTrackWoDCA status for pions"};

  // vertexing
  Configurable<bool> buildVertex{"buildVertex", false, "build vertext for Tcc"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any Lb is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};

  Configurable<float> softPiDcaXYMax{"softPiDcaXYMax", 0.065, "Soft pion max dcaXY (cm)"};
  Configurable<float> softPiDcaZMax{"softPiDcaZMax", 0.065, "Soft pion max dcaZ (cm)"};
  Configurable<float> deltaMassCanMax{"deltaMassCanMax", 2, "delta candidate max mass (DDPi-D0D0) ((GeV/c2)"};
  Configurable<float> massCanMax{"massCanMax", 4.0, "candidate max mass (DDPi) ((GeV/c2)"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  o2::vertexing::DCAFitterN<3> dfTcc; // Tcc vertex fitter
  o2::vertexing::DCAFitterN<2> dfDD;  // DD pair vertex fitter
  o2::vertexing::DCAFitterN<2> dfD1;  // 2-prong vertex fitter (to rebuild D01 vertex)
  o2::vertexing::DCAFitterN<2> dfD2;  // 2-prong vertex fitter (to rebuild D02 vertex)

  HfHelper hfHelper;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  double bz{0.};
  int runNumber{0};

  using TracksPid = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using TracksWPid = soa::Join<aod::TracksWCovDcaExtra, TracksPid, aod::TrackSelection>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

  using SelectedCandidatesMl = soa::Filtered<soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;

  Preslice<SelectedCandidatesMl> candsD0PerCollisionWithMl = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  // Partition<SelectedCandidatesMl> candidatesMlAll = aod::hf_sel_candidate_d0::isSelD0 >= 0;
  std::shared_ptr<TH1> hCandidatesD1, hCandidatesD2, hCandidatesTcc, hCandidatesDD;
  HistogramRegistry registry{"registry"};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "Tcc candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "Tcc candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  void init(InitContext const&)
  {

    std::array<bool, 3> doprocess{doprocessDataWithMl, doprocessDataWithMlWithFT0C, doprocessDataWithMlWithFT0M};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    if (buildVertex) {
      dfD1.setPropagateToPCA(propagateToPCA);
      dfD1.setMaxR(maxR);
      dfD1.setMaxDZIni(maxDZIni);
      dfD1.setMinParamChange(minParamChange);
      dfD1.setMinRelChi2Change(minRelChi2Change);
      dfD1.setUseAbsDCA(useAbsDCA);
      dfD1.setWeightedFinalPCA(useWeightedFinalPCA);

      dfD2.setPropagateToPCA(propagateToPCA);
      dfD2.setMaxR(maxR);
      dfD2.setMaxDZIni(maxDZIni);
      dfD2.setMinParamChange(minParamChange);
      dfD2.setMinRelChi2Change(minRelChi2Change);
      dfD2.setUseAbsDCA(useAbsDCA);
      dfD2.setWeightedFinalPCA(useWeightedFinalPCA);

      dfDD.setPropagateToPCA(propagateToPCA);
      dfDD.setMaxR(maxR);
      dfDD.setMaxDZIni(maxDZIni);
      dfDD.setMinParamChange(minParamChange);
      dfDD.setMinRelChi2Change(minRelChi2Change);
      dfDD.setUseAbsDCA(useAbsDCA);
      dfDD.setWeightedFinalPCA(useWeightedFinalPCA);

      dfTcc.setPropagateToPCA(propagateToPCA);
      dfTcc.setMaxR(maxR);
      dfTcc.setMaxDZIni(maxDZIni);
      dfTcc.setMinParamChange(minParamChange);
      dfTcc.setMinRelChi2Change(minRelChi2Change);
      dfTcc.setUseAbsDCA(useAbsDCA);
      dfTcc.setWeightedFinalPCA(useWeightedFinalPCA);

      // Configure CCDB access
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));

      hCandidatesD1 = registry.add<TH1>("hCandidatesD1", "D1 candidate counter", {HistType::kTH1D, {axisCands}});
      hCandidatesD2 = registry.add<TH1>("hCandidatesD2", "D2 candidate counter", {HistType::kTH1D, {axisCands}});
      hCandidatesTcc = registry.add<TH1>("hCandidatesTcc", "Tcc candidate counter", {HistType::kTH1D, {axisCands}});
      hCandidatesDD = registry.add<TH1>("hCandidatesDD", "DD pair candidate counter", {HistType::kTH1D, {axisCands}});

      setLabelHistoCands(hCandidatesD1);
      setLabelHistoCands(hCandidatesD2);
      setLabelHistoCands(hCandidatesTcc);
      setLabelHistoCands(hCandidatesDD);
    }
  }

  template <typename T>
  void fillEvent(const T& collision, int isEventReject, int runNumber)
  {
    rowCandidateFullEvents(
      collision.globalIndex(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      isEventReject,
      runNumber);
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  template <typename CollType, typename CandType, typename TrkType>
  void runCandCreatorData(CollType const& collisions,
                          CandType const& candidates,
                          aod::TrackAssoc const& trackIndices,
                          TrkType const& tracks)
  {

    for (const auto& collision : collisions) {
      auto primaryVertex = getPrimaryVertex(collision);
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      fillEvent(collision, 0, bc.runNumber());
      int64_t timeStamp = bc.timestamp();
      if (buildVertex) {
        if (runNumber != bc.runNumber()) {
          LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
          initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
          bz = o2::base::Propagator::Instance()->getNominalBz();
          LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        }
        dfTcc.setBz(bz);
        dfDD.setBz(bz);
        dfD1.setBz(bz);
        dfD2.setBz(bz);
      }

      o2::dataformats::V0 trackD1;
      o2::dataformats::V0 trackD2;
      auto thisCollId = collision.globalIndex();
      auto candwD0ThisColl = candidates.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      if (candwD0ThisColl.size() <= 1) {
        continue; // only loop the collision that include at least 2 D candidates
      }
      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);

      for (const auto& candidateD1 : candwD0ThisColl) {

        auto trackD1Prong0 = tracks.rawIteratorAt(candidateD1.prong0Id()); // positive daughter for D1
        auto trackD1Prong1 = tracks.rawIteratorAt(candidateD1.prong1Id()); // negative daughter for D1
        std::array<float, 3> const pVecD1Prong0{trackD1Prong0.pVector()};
        std::array<float, 3> const pVecD1Prong1{trackD1Prong1.pVector()};
        std::array<float, 3> pVecD1 = RecoDecay::pVec(pVecD1Prong0, pVecD1Prong1);

        for (auto candidateD2 = candidateD1 + 1; candidateD2 != candwD0ThisColl.end(); ++candidateD2) {
          // avoid shared tracks
          if (
            (candidateD1.prong0Id() == candidateD2.prong0Id()) ||
            (candidateD1.prong0Id() == candidateD2.prong1Id()) ||
            (candidateD1.prong1Id() == candidateD2.prong0Id()) ||
            (candidateD1.prong1Id() == candidateD2.prong1Id())) {
            continue;
          }

          auto trackD2Prong0 = tracks.rawIteratorAt(candidateD2.prong0Id()); // positive daughter for D2
          auto trackD2Prong1 = tracks.rawIteratorAt(candidateD2.prong1Id()); // negative daughter for D2
          std::array<float, 3> const pVecD2Prong0{trackD2Prong0.pVector()};
          std::array<float, 3> const pVecD2Prong1{trackD2Prong1.pVector()};
          std::array<float, 3> pVecD2 = RecoDecay::pVec(pVecD2Prong0, pVecD2Prong1);

          if (buildVertex) {
            auto trackParVarD1Prong0 = getTrackParCov(trackD1Prong0);
            auto trackParVarD1Prong1 = getTrackParCov(trackD1Prong1);
            auto dca0D1 = o2::dataformats::DCA(trackD1Prong0.dcaXY(), trackD1Prong0.dcaZ(), trackD1Prong0.cYY(), trackD1Prong0.cZY(), trackD1Prong0.cZZ());
            auto dca1D1 = o2::dataformats::DCA(trackD1Prong1.dcaXY(), trackD1Prong1.dcaZ(), trackD1Prong1.cYY(), trackD1Prong1.cZY(), trackD1Prong1.cZZ());

            // repropagate tracks to this collision if needed
            if (trackD1Prong0.collisionId() != thisCollId) {
              trackParVarD1Prong0.propagateToDCA(primaryVertex, bz, &dca0D1);
            }

            if (trackD1Prong1.collisionId() != thisCollId) {
              trackParVarD1Prong1.propagateToDCA(primaryVertex, bz, &dca1D1);
            }
            // reconstruct the 2-prong secondary vertex
            hCandidatesD1->Fill(SVFitting::BeforeFit);
            try {
              if (dfD1.process(trackParVarD1Prong0, trackParVarD1Prong1) == 0) {
                continue;
              }
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for first D0 cannot work, skipping the candidate.";
              hCandidatesD1->Fill(SVFitting::Fail);
              continue;
            }
            hCandidatesD1->Fill(SVFitting::FitOk);
            const auto& vertexD1 = dfD1.getPCACandidatePos();
            trackParVarD1Prong0.propagateTo(vertexD1[0], bz);
            trackParVarD1Prong1.propagateTo(vertexD1[0], bz);

            // build a D1 neutral track
            trackD1 = o2::dataformats::V0(vertexD1, pVecD1, dfD1.calcPCACovMatrixFlat(), trackParVarD1Prong0, trackParVarD1Prong1);

            auto trackParVarD2Prong0 = getTrackParCov(trackD2Prong0);
            auto trackParVarD2Prong1 = getTrackParCov(trackD2Prong1);
            auto dca0D2 = o2::dataformats::DCA(trackD2Prong0.dcaXY(), trackD2Prong0.dcaZ(), trackD2Prong0.cYY(), trackD2Prong0.cZY(), trackD2Prong0.cZZ());
            auto dca1D2 = o2::dataformats::DCA(trackD2Prong1.dcaXY(), trackD2Prong1.dcaZ(), trackD2Prong1.cYY(), trackD2Prong1.cZY(), trackD2Prong1.cZZ());

            // repropagate tracks to this collision if needed
            if (trackD2Prong0.collisionId() != thisCollId) {
              trackParVarD2Prong0.propagateToDCA(primaryVertex, bz, &dca0D2);
            }
            if (trackD2Prong1.collisionId() != thisCollId) {
              trackParVarD2Prong1.propagateToDCA(primaryVertex, bz, &dca1D2);
            }

            // reconstruct the 2-prong secondary vertex
            hCandidatesD2->Fill(SVFitting::BeforeFit);
            try {
              if (dfD2.process(trackParVarD2Prong0, trackParVarD2Prong1) == 0) {
                continue;
              }
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for second D0 cannot work, skipping the candidate.";
              hCandidatesD2->Fill(SVFitting::Fail);
              continue;
            }

            hCandidatesD2->Fill(SVFitting::FitOk);
            const auto& vertexD2 = dfD2.getPCACandidatePos();
            trackParVarD2Prong0.propagateTo(vertexD2[0], bz);
            trackParVarD2Prong1.propagateTo(vertexD2[0], bz);
            // build a D2 neutral track
            trackD2 = o2::dataformats::V0(vertexD2, pVecD2, dfD2.calcPCACovMatrixFlat(), trackParVarD2Prong0, trackParVarD2Prong1);

            hCandidatesDD->Fill(SVFitting::BeforeFit);
            try {
              if (dfDD.process(trackD1, trackD2) == 0) {
                continue;
              }
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for DD cannot work, skipping the candidate.";
              hCandidatesDD->Fill(SVFitting::Fail);
              continue;
            }
            hCandidatesDD->Fill(SVFitting::FitOk);
          }

          int candFlagD1 = -999;
          int candFlagD2 = -999;
          float cent = evaluateCentralityColl(collision);
          float massD01 = -999;
          float massD02 = -999;
          std::vector<float> mlScoresD1;
          std::vector<float> mlScoresD2;

          if (candidateD1.isSelD0()) {
            candFlagD1 = (candidateD1.isSelD0bar()) ? 3 : 1;
            std::copy(candidateD1.mlProbD0().begin(), candidateD1.mlProbD0().end(), std::back_inserter(mlScoresD1));
            massD01 = hfHelper.invMassD0ToPiK(candidateD1);
          }
          if (candidateD1.isSelD0bar() && !candidateD1.isSelD0()) {
            candFlagD1 = 2;
            std::copy(candidateD1.mlProbD0bar().begin(), candidateD1.mlProbD0bar().end(), std::back_inserter(mlScoresD1));
            massD01 = hfHelper.invMassD0barToKPi(candidateD1);
          }

          if (candidateD2.isSelD0()) {
            candFlagD2 = (candidateD2.isSelD0bar()) ? 3 : 1;
            std::copy(candidateD2.mlProbD0().begin(), candidateD2.mlProbD0().end(), std::back_inserter(mlScoresD2));
            massD02 = hfHelper.invMassD0ToPiK(candidateD2);
          }
          if (candidateD2.isSelD0bar() && !candidateD2.isSelD0()) {
            candFlagD2 = 2;
            std::copy(candidateD2.mlProbD0bar().begin(), candidateD2.mlProbD0bar().end(), std::back_inserter(mlScoresD2));
            massD02 = hfHelper.invMassD0barToKPi(candidateD2);
          }

          // const auto massD0D0Pair = RecoDecay::m(std::array{pVecD1, pVecD2}, std::array{MassD0, MassD0});

          rowCandidateDDPair(
            candidateD1.pxProng0(),
            candidateD1.pxProng1(),
            candidateD1.pyProng0(),
            candidateD1.pyProng1(),
            candidateD1.pzProng0(),
            candidateD1.pzProng1(),
            candidateD2.pxProng0(),
            candidateD2.pxProng1(),
            candidateD2.pyProng0(),
            candidateD2.pyProng1(),
            candidateD2.pzProng0(),
            candidateD2.pzProng1(),
            candFlagD1,
            candFlagD2,
            candidateD1.eta(),
            candidateD2.eta(),
            candidateD1.phi(),
            candidateD2.phi(),
            mlScoresD1[0],
            mlScoresD2[0],
            cent,
            collision.globalIndex(),
            timeStamp);

          // start to add the track of softpi to reconstruct Tcc
          for (const auto& trackId : trackIdsThisCollision) {

            auto trackPion = trackId.template track_as<TrkType>();
            if (usePionIsGlobalTrackWoDCA && !trackPion.isGlobalTrackWoDCA()) {
              continue;
            }
            // minimum pT selection
            if (trackPion.pt() < ptMinSoftPion) {
              continue;
            }
            if (std::abs(trackPion.dcaXY()) > softPiDcaXYMax || std::abs(trackPion.dcaZ()) > softPiDcaZMax) {
              continue;
            }
            // avoid shared tracks
            if (
              (candidateD1.prong0Id() == trackPion.globalIndex()) ||
              (candidateD1.prong1Id() == trackPion.globalIndex()) ||
              (candidateD2.prong0Id() == trackPion.globalIndex()) ||
              (candidateD2.prong1Id() == trackPion.globalIndex())) {
              continue;
            }

            std::array<float, 3> pVecSoftPi = {trackPion.pVector()};

            float impactParameterYD1 = -999.f;
            float impactParameterYD2 = -999.f;
            float impactParameterYSoftPi = -999.f;
            float chi2PCA = -999.f;
            if (buildVertex) {
              auto trackParCovPi = getTrackParCov(trackPion);
              // find the DCA between the D01, D02 and the bachelor track, for Tcc
              hCandidatesTcc->Fill(SVFitting::BeforeFit);
              try {
                if (dfTcc.process(trackD1, trackD2, trackParCovPi) == 0) {
                  continue;
                }
              } catch (const std::runtime_error& error) {
                LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for Tcc cannot work, skipping the candidate.";
                hCandidatesTcc->Fill(SVFitting::Fail);
                continue;
              }
              hCandidatesTcc->Fill(SVFitting::FitOk);

              dfTcc.propagateTracksToVertex();        // propagate the softpi and D0 pair to the Tcc vertex
              trackD1.getPxPyPzGlo(pVecD1);           // momentum of D1 at the Tcc vertex
              trackD2.getPxPyPzGlo(pVecD2);           // momentum of D2 at the Tcc vertex
              trackParCovPi.getPxPyPzGlo(pVecSoftPi); // momentum of pi at the Tcc vertex

              chi2PCA = dfTcc.getChi2AtPCACandidate();
              auto covMatrixPCA = dfTcc.calcPCACovMatrixFlat();
              hCovSVXX->Fill(covMatrixPCA[0]);

              // get track impact parameters
              // This modifies track momenta!
              auto covMatrixPV = primaryVertex.getCov();
              hCovPVXX->Fill(covMatrixPV[0]);
              o2::dataformats::DCA impactParameterD1;
              o2::dataformats::DCA impactParameterD2;
              o2::dataformats::DCA impactParameterSoftPi;

              trackD1.propagateToDCA(primaryVertex, bz, &impactParameterD1);
              trackD2.propagateToDCA(primaryVertex, bz, &impactParameterD2);
              trackParCovPi.propagateToDCA(primaryVertex, bz, &impactParameterSoftPi);

              impactParameterYD1 = impactParameterD1.getY();
              impactParameterYD2 = impactParameterD2.getY();
              impactParameterYSoftPi = impactParameterSoftPi.getY();
            }
            // Retrieve properties of the two D0 candidates
            float yD1 = hfHelper.yD0(candidateD1);
            float yD2 = hfHelper.yD0(candidateD2);
            float deltaMassD01 = -999;
            float deltaMassD02 = -999;

            std::array<double, 2> massD1Daus{MassPiPlus, MassKPlus};
            std::array<double, 2> massD2Daus{MassPiPlus, MassKPlus};

            if (candidateD1.isSelD0bar()) {

              massD1Daus[0] = MassKPlus;
              massD1Daus[1] = MassPiPlus;
            }
            if (candidateD2.isSelD0bar()) {
              massD2Daus[0] = MassKPlus;
              massD2Daus[1] = MassPiPlus;
            }

            auto massKpipi1 = RecoDecay::m(std::array{pVecD1Prong0, pVecD1Prong1, pVecSoftPi}, std::array{massD1Daus[0], massD1Daus[1], MassPiPlus});
            auto massKpipi2 = RecoDecay::m(std::array{pVecD2Prong0, pVecD2Prong1, pVecSoftPi}, std::array{massD2Daus[0], massD2Daus[1], MassPiPlus});
            auto arrayMomentaDDpi = std::array{pVecD1, pVecD2, pVecSoftPi};
            const auto massD0D0Pi = RecoDecay::m(arrayMomentaDDpi, std::array{MassD0, MassD0, MassPiPlus});
            const auto deltaMassD0D0Pi = massD0D0Pi - (massD01 + massD02);

            deltaMassD01 = massKpipi1 - massD01;
            deltaMassD02 = massKpipi2 - massD02;

            if (deltaMassD0D0Pi > deltaMassCanMax || massD0D0Pi > massCanMax) {
              continue;
            }

            rowCandidateLite(
              candidateD1.pxProng0(),
              candidateD1.pxProng1(),
              candidateD1.pyProng0(),
              candidateD1.pyProng1(),
              candidateD1.pzProng0(),
              candidateD1.pzProng1(),
              candidateD2.pxProng0(),
              candidateD2.pxProng1(),
              candidateD2.pyProng0(),
              candidateD2.pyProng1(),
              candidateD2.pzProng0(),
              candidateD2.pzProng1(),
              trackPion.px(),
              trackPion.py(),
              trackPion.pz(),
              candFlagD1,
              candFlagD2,
              massD01,
              massD02,
              deltaMassD01,
              deltaMassD02,
              massKpipi1,
              massKpipi2,
              massD0D0Pi,
              deltaMassD0D0Pi,
              candidateD1.eta(),
              candidateD2.eta(),
              trackPion.eta(),
              candidateD1.phi(),
              candidateD2.phi(),
              trackPion.phi(),
              yD1,
              yD2,
              trackPion.y(),
              trackPion.tpcNSigmaPi(),
              trackPion.tofNSigmaPi(),
              mlScoresD1[0],
              mlScoresD2[0],
              impactParameterYD1,
              impactParameterYD2,
              impactParameterYSoftPi,
              candidateD1.cpa(),
              candidateD2.cpa(),
              chi2PCA,
              trackPion.sign(),
              trackPion.dcaXY(),
              trackPion.dcaZ(),
              trackPion.itsNCls(),
              trackPion.tpcNClsCrossedRows(),
              trackPion.tpcChi2NCl(),
              cent,
              collision.globalIndex(),
              timeStamp);
          } // end of loop track
        } // end of loop second D0
      } // end of loop first D0
    } // end of loop collision
  }
  void processDataWithMl(Collisions const& collisions,
                         SelectedCandidatesMl const& candidates,
                         aod::TrackAssoc const& trackIndices,
                         TracksWPid const& tracks,
                         aod::BCsWithTimestamps const&)
  {
    rowCandidateFullEvents.reserve(collisions.size());
    runCandCreatorData(collisions, candidates, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTreeCreatorTccToD0D0Pi, processDataWithMl, "Process data with DCAFitterN with the ML method and without centrality", false);

  void processDataWithMlWithFT0C(CollisionsWithFT0C const& collisions,
                                 SelectedCandidatesMl const& candidates,
                                 aod::TrackAssoc const& trackIndices,
                                 TracksWPid const& tracks,
                                 aod::BCsWithTimestamps const&)
  {
    rowCandidateFullEvents.reserve(collisions.size());
    runCandCreatorData(collisions, candidates, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTreeCreatorTccToD0D0Pi, processDataWithMlWithFT0C, "Process data with DCAFitterN with the ML method and with FT0C centrality", true);

  void processDataWithMlWithFT0M(CollisionsWithFT0M const& collisions,
                                 SelectedCandidatesMl const& candidates,
                                 aod::TrackAssoc const& trackIndices,
                                 TracksWPid const& tracks,
                                 aod::BCsWithTimestamps const&)
  {
    rowCandidateFullEvents.reserve(collisions.size());
    runCandCreatorData(collisions, candidates, trackIndices, tracks);
  }
  PROCESS_SWITCH(HfTreeCreatorTccToD0D0Pi, processDataWithMlWithFT0M, "Process data with DCAFitterN with the ML method and with FT0M centrality", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorTccToD0D0Pi>(cfgc)};
}
