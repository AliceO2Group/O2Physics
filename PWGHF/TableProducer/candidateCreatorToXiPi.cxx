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

/// \file candidateCreatorToXiPi.cxx
/// \brief Reconstruction of Omegac0 and Xic0 -> xi pi candidates
/// \author Federica Zanone <federica.zanone@cern.ch>, HEIDELBERG UNIVERSITY & GSI

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::analysis;
// using namespace o2::analysis::pdg;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Reconstruction of omegac candidates
struct HfCandidateCreatorToXiPi {
  Produces<aod::HfCandToXiPi> rowCandidate;

  Configurable<bool> doPvRefit{"doPvRefit", false, "set to true if you do PV refit in trackIndexSkimCreator.cxx"};

  Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> maxDXYIni{"maxDXYIni", 4., "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};
  Configurable<bool> refitWithMatCorr{"refitWithMatCorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // cascade invariant mass cuts
  Configurable<bool> doCascadeInvMassCut{"doCascadeInvMassCut", false, "Use invariant mass cut to select cascade candidates"};
  Configurable<double> sigmaInvMassCascade{"sigmaInvMassCascade", 0.0025, "Invariant mass cut for cascade (sigma)"};
  Configurable<int> nSigmaInvMassCut{"nSigmaInvMassCut", 4, "Number of sigma for invariant mass cut"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber;

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using MyTracks = soa::Join<aod::TracksWCovDca, aod::HfPvRefitTrack>;
  using FilteredHfTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;
  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>; // to use strangeness tracking, use aod::TraCascDatas instead of aod::CascDatas
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == 0); // filter to use only HF selected collisions
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng >= 4);

  Preslice<MyTracks> tracksPerCollision = aod::track::collisionId;                                  // needed for PV refit
  Preslice<FilteredHfTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId; // aod::hf_track_association::collisionId
  Preslice<MyCascTable> cascadesPerCollision = aod::cascdata::collisionId;

  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  void process(SelectedCollisions const& collisions,
               aod::BCsWithTimestamps const& bcWithTimeStamps,
               MyTracks const& tracks,
               FilteredHfTrackAssocSel const& trackIndices,
               MyCascTable const& cascades,
               MyV0Table const&,
               aod::V0sLinked const&)
  {

    for (const auto& collision : collisions) {

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component

      // 2-prong vertex fitter to build the omegac vertex
      o2::vertexing::DCAFitterN<2> df;
      df.setBz(magneticField);
      df.setPropagateToPCA(propagateToPCA);
      df.setMaxR(maxR);
      df.setMaxDZIni(maxDZIni);
      df.setMaxDXYIni(maxDXYIni);
      df.setMinParamChange(minParamChange);
      df.setMinRelChi2Change(minRelChi2Change);
      df.setMaxChi2(maxChi2);
      df.setUseAbsDCA(useAbsDCA);
      df.setWeightedFinalPCA(useWeightedFinalPCA);
      df.setRefitWithMatCorr(refitWithMatCorr);

      double massPionFromPDG = o2::analysis::pdg::MassPiPlus;    // pdg code 211
      double massLambdaFromPDG = o2::analysis::pdg::MassLambda0; // pdg code 3122
      double massXiFromPDG = o2::analysis::pdg::MassXiMinus;     // pdg code 3312
      double massOmegacFromPDG = o2::analysis::pdg::MassOmegaC0; // pdg code 4332
      double massXicFromPDG = o2::analysis::pdg::MassXiCZero;    // pdg code 4132

      // loop over cascades reconstructed by cascadebuilder.cxx
      auto thisCollId = collision.globalIndex();
      auto groupedCascades = cascades.sliceBy(cascadesPerCollision, thisCollId);

      for (const auto& casc : groupedCascades) {

        //----------------accessing particles in the decay chain-------------
        // cascade daughter - charged particle
        // int indexTrackXiDauCharged = casc.bachelorId();     // pion <- xi index from cascade table (not used)
        auto trackXiDauCharged = casc.bachelor_as<MyTracks>(); // pion <- xi track from MyTracks table
        // cascade daughter - V0
        if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) { // check that V0 data are stored
          continue;
        }
        auto v0 = casc.v0_as<aod::V0sLinked>();
        auto v0Element = v0.v0Data_as<MyV0Table>(); // V0 element from LF table containing V0 info
        // V0 positive daughter
        auto trackV0Dau0 = v0Element.posTrack_as<MyTracks>(); // p <- V0 track (positive track) from MyTracks table
        // V0 negative daughter
        auto trackV0Dau1 = v0Element.negTrack_as<MyTracks>(); // pion <- V0 track (negative track) from MyTracks table

        // check that particles come from the same collision
        if (rejDiffCollTrack) {
          if (trackV0Dau0.collisionId() != trackV0Dau1.collisionId()) {
            continue;
          }
          if (trackXiDauCharged.collisionId() != trackV0Dau0.collisionId()) {
            continue;
          }
        }

        // use invariant mass cut to select cascades candidates
        if (doCascadeInvMassCut) {
          if (std::abs(casc.mXi() - massXiFromPDG) > (nSigmaInvMassCut * sigmaInvMassCascade)) {
            continue;
          }
        }

        //--------------------------reconstruct V0 track---------------------------
        // pseudorapidity
        double pseudorapV0PosDau = trackV0Dau0.eta();
        double pseudorapV0NegDau = trackV0Dau1.eta();

        // pion & p <- V0 tracks
        auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
        auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);

        // info from LF table
        std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()}; // pVec stands for vector containing the 3-momentum components
        std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
        std::array<float, 21> covV0 = {0.};
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covV0[MomInd[i]] = v0Element.momentumCovMat()[i];
          covV0[i] = v0Element.positionCovMat()[i];
        }
        // create V0 track
        auto trackV0 = o2::track::TrackParCov(vertexV0, pVecV0, covV0, 0, true);
        trackV0.setAbsCharge(0);
        trackV0.setPID(o2::track::PID::Lambda);

        std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
        std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

        auto trackV0Copy = trackV0;

        //-----------------------------reconstruct cascade track-----------------------------
        // pseudorapidity
        double pseudorapPiFromCas = trackXiDauCharged.eta();

        // pion <- casc track to be processed with DCAfitter
        auto trackParCovXiDauCharged = getTrackParCov(trackXiDauCharged);

        // info from LF table
        std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
        std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
        std::array<float, 21> covCasc = {0.};
        for (int i = 0; i < 6; i++) {
          covCasc[MomInd[i]] = casc.momentumCovMat()[i];
          covCasc[i] = casc.positionCovMat()[i];
        }
        // create cascade track
        o2::track::TrackParCov trackCasc;
        if (trackXiDauCharged.sign() > 0) {
          trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
        } else if (trackXiDauCharged.sign() < 0) {
          trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
        } else {
          continue;
        }
        trackCasc.setAbsCharge(1);
        trackCasc.setPID(o2::track::PID::XiMinus);

        std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

        auto trackCascCopy = trackCasc;

        //-------------------combining cascade and pion tracks--------------------------
        auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackIndexPion : groupedTrackIndices) {

          auto trackPion = trackIndexPion.track_as<MyTracks>();

          if ((rejDiffCollTrack) && (trackXiDauCharged.collisionId() != trackPion.collisionId())) {
            continue;
          }

          // ask for opposite sign daughters (omegac daughters)
          if (trackPion.sign() * trackXiDauCharged.sign() >= 0) {
            continue;
          }

          // check not to take the same particle twice in the decay chain
          if (trackPion.globalIndex() == trackXiDauCharged.globalIndex() || trackPion.globalIndex() == trackV0Dau0.globalIndex() || trackPion.globalIndex() == trackV0Dau1.globalIndex()) {
            continue;
          }

          // pseudorapidity
          double pseudorapPiFromOme = trackPion.eta();

          // primary pion track to be processed with DCAFitter
          auto trackParVarPi = getTrackParCov(trackPion);
          auto trackParVarPiCopy = trackParVarPi;

          // reconstruct omegac with DCAFitter
          int nVtxFromFitterOmegac = df.process(trackCasc, trackParVarPi);
          if (nVtxFromFitterOmegac == 0) {
            continue;
          }
          auto vertexOmegacFromFitter = df.getPCACandidate();
          auto chi2PCAOmegac = df.getChi2AtPCACandidate();
          std::array<float, 3> pVecCascAsD;
          std::array<float, 3> pVecPionFromOmegac;
          df.propagateTracksToVertex();
          if (!df.isPropagateTracksToVertexDone()) {
            continue;
          }
          df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
          df.getTrack(1).getPxPyPzGlo(pVecPionFromOmegac);
          std::array<float, 3> pVecOmegac = {pVecCascAsD[0] + pVecPionFromOmegac[0], pVecCascAsD[1] + pVecPionFromOmegac[1], pVecCascAsD[2] + pVecPionFromOmegac[2]};

          std::array<float, 3> coordVtxOmegac = df.getPCACandidatePos();
          std::array<float, 6> covVtxOmegac = df.calcPCACovMatrixFlat();

          // create omegac track
          o2::track::TrackParCov trackOmegac = df.createParentTrackParCov();
          trackOmegac.setAbsCharge(0);

          // DCAxy (computed with propagateToDCABxByBz method)
          float dcaxyV0Dau0 = trackV0Dau0.dcaXY();
          float dcaxyV0Dau1 = trackV0Dau1.dcaXY();
          float dcaxyPiFromCasc = trackXiDauCharged.dcaXY();

          // DCAz (computed with propagateToDCABxByBz method)
          float dcazV0Dau0 = trackV0Dau0.dcaZ();
          float dcazV0Dau1 = trackV0Dau1.dcaZ();
          float dcazPiFromCasc = trackXiDauCharged.dcaZ();

          // primary vertex of the collision
          auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
          std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

          if (doPvRefit && ((trackPion.pvRefitSigmaX2() != 1e10f) || (trackPion.pvRefitSigmaY2() != 1e10f) || (trackPion.pvRefitSigmaZ2() != 1e10f))) { // if I asked for PV refit in trackIndexSkimCreator.cxx
            pvCoord[0] = trackPion.pvRefitX();
            pvCoord[1] = trackPion.pvRefitY();
            pvCoord[2] = trackPion.pvRefitZ();

            // o2::dataformats::VertexBase Pvtx;
            primaryVertex.setX(trackPion.pvRefitX());
            primaryVertex.setY(trackPion.pvRefitY());
            primaryVertex.setZ(trackPion.pvRefitZ());
            primaryVertex.setCov(trackPion.pvRefitSigmaX2(), trackPion.pvRefitSigmaXY(), trackPion.pvRefitSigmaY2(), trackPion.pvRefitSigmaXZ(), trackPion.pvRefitSigmaYZ(), trackPion.pvRefitSigmaZ2());

            o2::dataformats::DCA impactParameterV0Dau0;
            o2::dataformats::DCA impactParameterV0Dau1;
            o2::dataformats::DCA impactParameterPiFromCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovXiDauCharged, 2.f, matCorr, &impactParameterPiFromCasc);
            dcaxyV0Dau0 = impactParameterV0Dau0.getY();
            dcaxyV0Dau1 = impactParameterV0Dau1.getY();
            dcaxyPiFromCasc = impactParameterPiFromCasc.getY();
            dcazV0Dau0 = impactParameterV0Dau0.getZ();
            dcazV0Dau1 = impactParameterV0Dau1.getZ();
            dcazPiFromCasc = impactParameterPiFromCasc.getZ();
          }

          // impact parameters
          o2::dataformats::DCA impactParameterCasc;
          o2::dataformats::DCA impactParameterPrimaryPi;
          o2::dataformats::DCA impactParameterV0;
          o2::dataformats::DCA impactParameterOmegac;
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCascCopy, 2.f, matCorr, &impactParameterCasc);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarPiCopy, 2.f, matCorr, &impactParameterPrimaryPi);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackV0Copy, 2.f, matCorr, &impactParameterV0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackOmegac, 2.f, matCorr, &impactParameterOmegac);

          // invariant mass under the hypothesis of particles ID corresponding to the decay chain
          double mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
          double mCasc = casc.mXi();
          const std::array<double, 2> arrMassOmegac = {massXiFromPDG, massPionFromPDG};
          double mOmegac = RecoDecay::m(std::array{pVecCascAsD, pVecPionFromOmegac}, arrMassOmegac);

          // computing cosPA
          double cpaV0 = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
          double cpaOmegac = RecoDecay::cpa(pvCoord, coordVtxOmegac, pVecOmegac);
          double cpaCasc = RecoDecay::cpa(coordVtxOmegac, vertexCasc, pVecCasc);
          double cpaxyV0 = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
          double cpaxyOmegac = RecoDecay::cpaXY(pvCoord, coordVtxOmegac, pVecOmegac);
          double cpaxyCasc = RecoDecay::cpaXY(coordVtxOmegac, vertexCasc, pVecCasc);

          // computing decay length and ctau
          double decLenOmegac = RecoDecay::distance(pvCoord, coordVtxOmegac);
          double decLenCascade = RecoDecay::distance(coordVtxOmegac, vertexCasc);
          double decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);
          double ctOmegac = RecoDecay::ct(pVecOmegac, decLenOmegac, massOmegacFromPDG);
          double ctXic = RecoDecay::ct(pVecOmegac, decLenOmegac, massXicFromPDG);
          double ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, massXiFromPDG);
          double ctV0 = RecoDecay::ct(pVecV0, decLenV0, massLambdaFromPDG);

          // computing eta
          double pseudorapOmegac = RecoDecay::eta(pVecOmegac);
          double pseudorapCascade = RecoDecay::eta(pVecCasc);
          double pseudorapV0 = RecoDecay::eta(pVecV0);

          // DCA between daughters
          float dcaCascDau = casc.dcacascdaughters();
          float dcaV0Dau = casc.dcaV0daughters();
          float dcaOmegacDau = std::sqrt(df.getChi2AtPCACandidate());

          // set hfFlag
          int hfFlag = 1 << aod::hf_cand_toxipi::DecayType::DecayToXiPi;

          // fill test histograms
          hInvMassOmegac->Fill(mOmegac);

          // fill the table
          rowCandidate(collision.globalIndex(),
                       pvCoord[0], pvCoord[1], pvCoord[2],
                       vertexOmegacFromFitter[0], vertexOmegacFromFitter[1], vertexOmegacFromFitter[2],
                       vertexCasc[0], vertexCasc[1], vertexCasc[2],
                       vertexV0[0], vertexV0[1], vertexV0[2],
                       trackXiDauCharged.sign(),
                       chi2PCAOmegac, covVtxOmegac[0], covVtxOmegac[1], covVtxOmegac[2], covVtxOmegac[3], covVtxOmegac[4], covVtxOmegac[5],
                       covV0[0], covV0[1], covV0[2], covV0[3], covV0[4], covV0[5],
                       covCasc[0], covCasc[1], covCasc[2], covCasc[3], covCasc[4], covCasc[5],
                       pVecOmegac[0], pVecOmegac[1], pVecOmegac[2],
                       pVecCasc[0], pVecCasc[1], pVecCasc[2],
                       pVecPionFromOmegac[0], pVecPionFromOmegac[1], pVecPionFromOmegac[2],
                       pVecV0[0], pVecV0[1], pVecV0[2],
                       pVecPionFromCasc[0], pVecPionFromCasc[1], pVecPionFromCasc[2],
                       pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                       pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                       impactParameterCasc.getY(), impactParameterPrimaryPi.getY(),
                       impactParameterCasc.getZ(), impactParameterPrimaryPi.getZ(),
                       impactParameterV0.getY(), impactParameterV0.getZ(),
                       std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPrimaryPi.getSigmaY2()), std::sqrt(impactParameterV0.getSigmaY2()),
                       v0Element.globalIndex(), v0Element.posTrackId(), v0Element.negTrackId(),
                       casc.globalIndex(), trackPion.globalIndex(), trackXiDauCharged.globalIndex(),
                       impactParameterOmegac.getY(), impactParameterOmegac.getZ(),
                       mLambda, mCasc, mOmegac,
                       cpaV0, cpaOmegac, cpaCasc, cpaxyV0, cpaxyOmegac, cpaxyCasc,
                       ctOmegac, ctCascade, ctV0, ctXic,
                       pseudorapV0PosDau, pseudorapV0NegDau, pseudorapPiFromCas, pseudorapPiFromOme,
                       pseudorapOmegac, pseudorapCascade, pseudorapV0,
                       dcaxyV0Dau0, dcaxyV0Dau1, dcaxyPiFromCasc,
                       dcazV0Dau0, dcazV0Dau1, dcazPiFromCasc,
                       dcaCascDau, dcaV0Dau, dcaOmegacDau, hfFlag);

        } // loop over pions
      }   // loop over cascades
    }     // close loop collisions
  }       // end of process
};        // end of struct

/// Performs MC matching.
struct HfCandidateCreatorToXiPiMc {
  Produces<aod::HfToXiPiMCRec> rowMCMatchRec;
  Produces<aod::HfToXiPiMCGen> rowMCMatchGen;

  Configurable<bool> matchOmegacMc{"matchOmegacMc", true, "Do MC matching for Omegac0"};
  Configurable<bool> matchXicMc{"matchXicMc", false, "Do MC matching for Xic0"};

  void init(InitContext const&) {}

  void processDoNoMc(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorToXiPiMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandToXiPi const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    int indexRec = -1;
    int8_t sign = -9;
    int8_t flag = -9;
    // int8_t origin = 0; //to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenXi = 0;
    int8_t debugGenLambda = 0;

    int pdgCodeOmegac0 = pdg::Code::kOmegaC0; // 4332
    int pdgCodeXic0 = pdg::Code::kXiCZero;    // 4132
    int pdgCodeXiMinus = kXiMinus;            // 3312
    int pdgCodeLambda = kLambda0;             // 3122
    int pdgCodePiPlus = kPiPlus;              // 211
    int pdgCodePiMinus = kPiMinus;            // -211
    int pdgCodeProton = kProton;              // 2212

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      // origin = 0;
      debug = 0;
      auto arrayDaughters = std::array{candidate.primaryPi_as<aod::TracksWMc>(), // pi <- omegac
                                       candidate.bachelor_as<aod::TracksWMc>(),  // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),  // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Omegac matching
      if (matchOmegacMc) {
        // Omegac → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi);
            }
          }
        }
      }
      // Xic matching
      if (matchXicMc) {
        // Xic → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeXic0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi);
            }
          }
        }
      }

      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = -9;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenXi = 0;
      debugGenLambda = 0;
      // origin = 0;
      if (matchOmegacMc) {
        //  Omegac → Xi pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          // Xi -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            debugGenXi = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi);
            }
          }
        }
      }
      if (matchXicMc) {
        //  Xic → Xi pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeXic0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          // Xi- -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            debugGenXi = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi);
            }
          }
        }
      }

      // rowMCMatchGen(flag, origin);
      rowMCMatchGen(flag, debugGenCharmBar, debugGenXi, debugGenLambda);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorToXiPiMc, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorToXiPi>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorToXiPiMc>(cfgc)};
}
