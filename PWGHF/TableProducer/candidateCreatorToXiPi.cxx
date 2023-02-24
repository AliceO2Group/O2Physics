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
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Core/SelectorCuts.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::pdg;
using namespace o2::aod::v0data;
using namespace o2::aod::cascdata;
using namespace o2::aod::hf_track_index;
using namespace o2::aod::hf_sel_collision;
using namespace o2::aod::hf_cand_toxipi;

// Reconstruction of omegac candidates
struct HfCandidateCreatorToXiPi {
  Produces<aod::HfCandToXiPi> rowCandidate;

  Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> refitWithMatCorr{"refitWithMatCorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;

  // filter to use only HF selected collisions
  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == 0);

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using MyTracks = soa::Join<aod::BigTracks, aod::TracksDCA, aod::TrackSelection>;
  using MyCascTable = soa::Join<aod::CascDataExt, aod::CascCovs>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;

  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  void process(SelectedCollisions::iterator const& collision,
               aod::BCsWithTimestamps const& bcWithTimeStamps,
               MyCascTable const& cascades,
               MyTracks const& tracks,
               MyV0Table const&,
               aod::V0sLinked const&)
  {

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
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);
    df.setRefitWithMatCorr(refitWithMatCorr);

    double massPionFromPDG = RecoDecay::getMassPDG(kPiPlus);    // pdg code 211
    double massLambdaFromPDG = RecoDecay::getMassPDG(kLambda0); // pdg code 3122
    double massXiFromPDG = RecoDecay::getMassPDG(kXiMinus);     // pdg code 3312
    double massOmegacFromPDG = RecoDecay::getMassPDG(kOmegaC0); // pdg code 4332
    double massXicFromPDG = RecoDecay::getMassPDG(kXiCZero);    // pdg code 4132

    // loop over cascades reconstructed by cascadebuilder.cxx
    for (auto const& casc : cascades) {

      if (collision.globalIndex() != casc.collisionId()) { // check to be further processed when the problem of ambiguous tracks will be solved
        continue;
      }

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
      if (rejDiffCollTrack) { // check to be further processed when the problem of ambiguous tracks will be solved
        if (trackV0Dau0.collisionId() != trackV0Dau1.collisionId()) {
          continue;
        }
        if (trackXiDauCharged.collisionId() != trackV0Dau0.collisionId()) {
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
      std::array<float, 3> pVecV0 = {v0Element.px(), v0Element.py(), v0Element.pz()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {v0Element.x(), v0Element.y(), v0Element.z()};
      const std::array<float, 6> covVtxV0 = {v0Element.positionCovMat()[0], v0Element.positionCovMat()[1], v0Element.positionCovMat()[2], v0Element.positionCovMat()[3], v0Element.positionCovMat()[4], v0Element.positionCovMat()[5]};

      std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      // create V0 track
      auto trackV0 = o2::dataformats::V0(vertexV0, pVecV0, covVtxV0, trackParCovV0Dau0, trackParCovV0Dau1, {0, 0}, {0, 0});
      auto trackV0Copy = trackV0;

      //-----------------------------reconstruct cascade track-----------------------------
      // pseudorapidity
      double pseudorapPiFromCas = trackXiDauCharged.eta();

      // info from LF table
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      const std::array<float, 6> covVtxCasc = {casc.positionCovMat()[0], casc.positionCovMat()[1], casc.positionCovMat()[2], casc.positionCovMat()[3], casc.positionCovMat()[4], casc.positionCovMat()[5]};

      std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      // pion <- casc track to be processed with DCAfitter
      auto trackParVarXiDauCharged = getTrackParCov(trackXiDauCharged);

      // create cascade track
      auto trackCasc = o2::dataformats::V0(vertexCasc, pVecCasc, covVtxCasc, trackV0, trackParVarXiDauCharged, {0, 0}, {0, 0});
      auto trackCascCopy = trackCasc;

      //-------------------combining cascade and pion tracks--------------------------
      for (auto const& trackPion : tracks) {

        if ((rejDiffCollTrack) && (trackXiDauCharged.collisionId() != trackPion.collisionId())) { // check to be further processed when the problem of ambiguous tracks will be solved
          continue;
        }

        // ask for opposite sign daughters (omegac daughters)
        if (trackPion.sign() * trackXiDauCharged.sign() >= 0) {
          continue;
        }

        // check not to take the same particle twice in the decay chain
        if (trackPion.globalIndex() == trackXiDauCharged.globalIndex() || trackPion.globalIndex() == trackV0Dau0.globalIndex() || trackPion.globalIndex() == trackV0Dau1.globalIndex() || trackPion.globalIndex() == casc.globalIndex()) {
          continue;
        }

        // pseudirapidity
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
        auto trackOmegac = o2::dataformats::V0(coordVtxOmegac, pVecOmegac, covVtxOmegac, trackCasc, trackParVarPi, {0, 0}, {0, 0});

        // impact parameter omegac
        o2::dataformats::DCA impactParameterOmegac;
        auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
        trackOmegac.propagateToDCA(primaryVertex, magneticField, &impactParameterOmegac);

        // impact parameter
        o2::dataformats::DCA impactParameterCasc;
        o2::dataformats::DCA impactParameterPrimaryPi;
        o2::dataformats::DCA impactParameterV0;
        trackCascCopy.propagateToDCA(primaryVertex, magneticField, &impactParameterCasc);
        trackParVarPiCopy.propagateToDCA(primaryVertex, magneticField, &impactParameterPrimaryPi);
        trackV0Copy.propagateToDCA(primaryVertex, magneticField, &impactParameterV0);

        // DCAxy
        double dcaxyPrimaryPi = trackPion.dcaXY();
        double dcaxyV0Dau0 = trackV0Dau0.dcaXY();
        double dcaxyV0Dau1 = trackV0Dau1.dcaXY();
        double dcaxyCascDau = trackXiDauCharged.dcaXY();

        // invariant mass under the hypothesis of particles ID corresponding to the decay chain
        double mLambda = v0Element.mLambda();         // from LF table, V0 mass under lambda hypothesis
        double mAntiLambda = v0Element.mAntiLambda(); // from LF table, V0 mass under anti-lambda hypothesis
        double mCasc = casc.mXi();
        const std::array<double, 2> arrMassOmegac = {massXiFromPDG, massPionFromPDG};
        double mOmegac = RecoDecay::m(std::array{pVecCascAsD, pVecPionFromOmegac}, arrMassOmegac);

        // computing cosPA
        double cpaV0 = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
        double cpaOmegac = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac, pVecOmegac);
        double cpaCasc = RecoDecay::cpa(coordVtxOmegac, vertexCasc, pVecCasc);
        double cpaxyV0 = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
        double cpaxyOmegac = RecoDecay::cpaXY(std::array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac, pVecOmegac);
        double cpaxyCasc = RecoDecay::cpaXY(coordVtxOmegac, vertexCasc, pVecCasc);

        // computing decay length and ctau
        double decLenOmegac = RecoDecay::distance(std::array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac);
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
        double dcaCascDau = casc.dcacascdaughters();
        double dcaV0Dau = casc.dcaV0daughters();
        double dcaOmegacDau = std::sqrt(df.getChi2AtPCACandidate());

        // set hfFlag
        int hfFlag = 1 << DecayType::DecayToXiPi;

        // fill test histograms

        hInvMassOmegac->Fill(mOmegac);

        // fill the table
        rowCandidate(collision.globalIndex(),
                     collision.posX(), collision.posY(), collision.posZ(),
                     vertexOmegacFromFitter[0], vertexOmegacFromFitter[1], vertexOmegacFromFitter[2],
                     vertexCasc[0], vertexCasc[1], vertexCasc[2],
                     vertexV0[0], vertexV0[1], vertexV0[2],
                     trackXiDauCharged.sign(),
                     chi2PCAOmegac, covVtxOmegac[0], covVtxOmegac[1], covVtxOmegac[2], covVtxOmegac[3], covVtxOmegac[4], covVtxOmegac[5],
                     covVtxV0[0], covVtxV0[1], covVtxV0[2], covVtxV0[3], covVtxV0[4], covVtxV0[5],
                     covVtxCasc[0], covVtxCasc[1], covVtxCasc[2], covVtxCasc[3], covVtxCasc[4], covVtxCasc[5],
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
                     mLambda, mAntiLambda, mCasc, mOmegac,
                     cpaV0, cpaOmegac, cpaCasc, cpaxyV0, cpaxyOmegac, cpaxyCasc,
                     ctOmegac, ctCascade, ctV0, ctXic,
                     pseudorapV0PosDau, pseudorapV0NegDau, pseudorapPiFromCas, pseudorapPiFromOme,
                     pseudorapOmegac, pseudorapCascade, pseudorapV0,
                     dcaxyPrimaryPi, dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascDau,
                     dcaCascDau, dcaV0Dau, dcaOmegacDau, hfFlag);

      } // loop over pions
    }   // loop over candidates
  }     // end of process
};      // end of struct

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
                 aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    // int8_t origin = 0; //to be used for prompt/non prompt
    int8_t debug = 0;

    int pdgCodeOmegac0 = pdg::Code::kOmegaC0; // 4332
    int pdgCodeXic0 = pdg::Code::kXiCZero;    // 4132
    int pdgCodeXiMinus = kXiMinus;            // 3312
    int pdgCodeLambda = kLambda0;             // 3122
    int pdgCodePiPlus = kPiPlus;              // 211
    int pdgCodePiMinus = kPiMinus;            // -211
    int pdgCodeProton = kProton;              // 2212

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      // Printf("New rec. candidate");
      flag = 0;
      // origin = 0;
      debug = 0;
      auto arrayDaughters = std::array{candidate.primaryPi_as<aod::BigTracksMC>(), // pi <- omegac
                                       candidate.bachelor_as<aod::BigTracksMC>(),  // pi <- cascade
                                       candidate.posTrack_as<aod::BigTracksMC>(),  // p <- lambda
                                       candidate.negTrack_as<aod::BigTracksMC>()}; // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::BigTracksMC>(),
                                           candidate.posTrack_as<aod::BigTracksMC>(),
                                           candidate.negTrack_as<aod::BigTracksMC>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::BigTracksMC>(),
                                         candidate.negTrack_as<aod::BigTracksMC>()};

      // Omegac matching
      if (matchOmegacMc) {
        // Omegac → pi pi pi p
        // Printf("Checking Omegac → pi pi pi p");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // cascade → lambda pi
          // Printf("Checking cascade → pi pi p");
          indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // v0 → p pi
            // Printf("Checking v0 → p pi");
            indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << DecayType::OmegaczeroToXiPi);
            }
          }
        }

        // Xic matching
      }
      if (matchXicMc) {
        // Xic → pi pi pi p
        // Printf("Checking Xic → pi pi pi p");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdgCodeXic0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // cascade → lambda pi
          // Printf("Checking cascade → pi pi p");
          indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // v0 → p pi
            // Printf("Checking v0 → p pi");
            indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << DecayType::XiczeroToXiPi);
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
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      // origin = 0;
      if (matchOmegacMc) {
        //  Omegac → Xi pi
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true)) {
          // Match Xi -> lambda pi
          auto cascMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking cascade → lambda pi");
          if (RecoDecay::isMatchedMCGen(particlesMC, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            // lambda -> p pi
            auto v0MC = particlesMC.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(particlesMC, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign)) {
              flag = sign * (1 << DecayType::OmegaczeroToXiPi);
            }
          }
        }
      }
      if (matchXicMc) {
        //  Xic → Xi pi
        if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdgCodeXic0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true)) {
          // Match Xi -> lambda pi
          auto cascMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
          // Printf("Checking cascade → lambda pi");
          if (RecoDecay::isMatchedMCGen(particlesMC, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            // lambda -> p pi
            auto v0MC = particlesMC.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(particlesMC, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign)) {
              flag = sign * (1 << DecayType::XiczeroToXiPi);
            }
          }
        }
      }

      // rowMCMatchGen(flag, origin);
      rowMCMatchGen(flag);
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
