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

/// \file candidateCreatorOmegac.cxx
/// \brief Reconstruction of Omegac0 candidates
/// \author Federica Zanone <federica.zanone@cern.ch>, HEIDELBERG UNIVERSITY & GSI

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsVertexing/DCAFitterN.h"
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
using namespace o2::aod::hf_cand_omegac;

// Reconstruction of omegac candidates
struct HfCandidateCreatorOmegac {
  Produces<aod::HfCandOmegac> rowCandidate;

  Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "use absolute DCA for vertexing"};
  Configurable<bool> useWeightedPCA{"useWeightedPCA", false, "vertices use cov matrices"};
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

  OutputObj<TH1F> hPtPrimaryPi{TH1F("hPtPrimaryPi", "p_T primary #pi;p_T (GeV/#it{c});entries", 500, 0, 20)};
  OutputObj<TH1F> hxVertexOmegac{TH1F("hxVertexOmegac", "x Omegac vertex;xVtx;entries", 500, -10, 10)};
  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hMassOmegacNotFixed{TH1F("hMassOmegacNotFixed", "hMassOmegacNotFixed;invmass;entries", 500, 2.2, 3.1)};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;
  }

  void process(SelectedCollisions::iterator const& collision,
               aod::BCsWithTimestamps const& bcWithTimeStamps,
               aod::CascDataExt const& cascades,
               MyTracks const& tracks,
               aod::V0Datas const&,
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
    df.setWeightedFinalPCA(useWeightedPCA);
    df.setRefitWithMatCorr(refitWithMatCorr);

    // 2-prong vertex fitter to build the cascade vertex
    o2::vertexing::DCAFitterN<2> dfc;
    dfc.setBz(magneticField);
    dfc.setPropagateToPCA(propagateToPCA);
    dfc.setMaxR(maxR);
    dfc.setMaxDZIni(maxDZIni);
    dfc.setMinParamChange(minParamChange);
    dfc.setMinRelChi2Change(minRelChi2Change);
    dfc.setUseAbsDCA(useAbsDCA);
    dfc.setWeightedFinalPCA(useWeightedPCA);
    dfc.setRefitWithMatCorr(refitWithMatCorr);

    // 2-prong vertex fitter to build the V0 vertex
    o2::vertexing::DCAFitterN<2> dfv;
    dfv.setBz(magneticField);
    dfv.setPropagateToPCA(propagateToPCA);
    dfv.setMaxR(maxR);
    dfv.setMaxDZIni(maxDZIni);
    dfv.setMinParamChange(minParamChange);
    dfv.setMinRelChi2Change(minRelChi2Change);
    dfv.setUseAbsDCA(useAbsDCA);
    dfv.setWeightedFinalPCA(useWeightedPCA);
    dfv.setRefitWithMatCorr(refitWithMatCorr);

    double massPionFromPDG = RecoDecay::getMassPDG(kPiPlus);    // pdg code 211
    double massProtonFromPDG = RecoDecay::getMassPDG(kProton);  // pdg code 2212
    double massLambdaFromPDG = RecoDecay::getMassPDG(kLambda0); // pdg code 3122
    double massXiFromPDG = RecoDecay::getMassPDG(kXiMinus);     // pdg code 3312
    double massOmegacFromPDG = RecoDecay::getMassPDG(kOmegaC0); // pdg code 4332

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
      // int indexV0 = casc.v0Id(); // VO index from cascades table (not used)
      auto v0 = casc.v0_as<aod::V0sLinked>();
      auto v0Element = v0.v0Data(); // V0 element from LF table containing V0 info
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

      //--------------------------reconstruct V0---------------------------
      // pseudorapidity
      double pseudorapV0PosDau = trackV0Dau0.eta();
      double pseudorapV0NegDau = trackV0Dau1.eta();

      // pion & p <- V0 track to be processed with DCAfitter
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);

      // info from LF table (not used in the task - for cross checks)
      // std::array<float, 3> pvecV0LFTable = {v0Element.px(), v0Element.py(), v0Element.pz()}; // pvec stands for vector containing the 3-momentum components
      // std::array<float, 3> vertexV0 = {v0Element.x(), v0Element.y(), v0Element.z()};

      // reconstruct V0 with DCAFitter
      int nVtxFromFitterV0 = dfv.process(trackParCovV0Dau0, trackParCovV0Dau1);
      if (nVtxFromFitterV0 == 0) {
        continue;
      }
      auto vertexV0FromFitter = dfv.getPCACandidate();
      auto chi2PCAV0 = dfv.getChi2AtPCACandidate();
      std::array<float, 3> pvecV0Dau0;
      std::array<float, 3> pvecV0Dau1;
      dfv.propagateTracksToVertex();
      if (!dfv.isPropagateTracksToVertexDone()) {
        continue;
      }
      dfv.getTrack(0).getPxPyPzGlo(pvecV0Dau0);
      dfv.getTrack(1).getPxPyPzGlo(pvecV0Dau1);
      std::array<float, 3> pvecV0AsM = {pvecV0Dau0[0] + pvecV0Dau1[0], pvecV0Dau0[1] + pvecV0Dau1[1], pvecV0Dau0[2] + pvecV0Dau1[2]}; // AsM stands for "as mother"

      std::array<float, 3> coordVtxV0 = dfv.getPCACandidatePos();
      std::array<float, 6> covVtxV0 = dfv.calcPCACovMatrixFlat();

      // create V0 track
      auto trackV0 = o2::dataformats::V0(coordVtxV0, pvecV0AsM, covVtxV0, trackParCovV0Dau0, trackParCovV0Dau1, {0, 0}, {0, 0});
      auto trackV0Copy = trackV0;

      //-----------------------------reconstruct cascade------------------------------
      // pseudorapidity
      double pseudorapPiFromCas = trackXiDauCharged.eta();

      // info from LF table
      std::array<float, 3> vertexCascLFTable = {casc.x(), casc.y(), casc.z()};

      // pion <- casc track to be processed with DCAfitter
      auto trackParVarXiDauCharged = getTrackParCov(trackXiDauCharged);

      // reconstruct cascade with DCAFitter
      int nVtxFromFitterCasc = dfc.process(trackV0, trackParVarXiDauCharged);
      if (nVtxFromFitterCasc == 0) {
        continue;
      }
      auto vertexCascFromFitter = dfc.getPCACandidate();
      auto chi2PCACascade = dfc.getChi2AtPCACandidate();
      std::array<float, 3> pvecV0AsD; // AsD stands for "as daughter"
      std::array<float, 3> pvecPionFromCasc;
      dfc.propagateTracksToVertex();
      if (!dfc.isPropagateTracksToVertexDone()) {
        continue;
      }
      dfc.getTrack(0).getPxPyPzGlo(pvecV0AsD);
      dfc.getTrack(1).getPxPyPzGlo(pvecPionFromCasc);
      std::array<float, 3> pvecCascAsM = {pvecV0AsD[0] + pvecPionFromCasc[0], pvecV0AsD[1] + pvecPionFromCasc[1], pvecV0AsD[2] + pvecPionFromCasc[2]};

      std::array<float, 3> coordVtxCasc = dfc.getPCACandidatePos();
      std::array<float, 6> covVtxCasc = dfc.calcPCACovMatrixFlat();

      // create cascade track
      auto trackCasc = o2::dataformats::V0(coordVtxCasc, pvecCascAsM, covVtxCasc, trackV0, trackParVarXiDauCharged, {0, 0}, {0, 0});
      auto trackCascCopy = trackCasc;

      //-------------------combining cascade and pion tracks--------------------------
      for (auto const& trackPion : tracks) {

        if ((rejDiffCollTrack) && (trackXiDauCharged.collisionId() != trackPion.collisionId())) {
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
        auto vertexOmegacFromFitter = df.getPCACandidate(); // use df.calcPCACovMatrixFlat() to get the covariance matrix
        auto chi2PCAOmegac = df.getChi2AtPCACandidate();
        std::array<float, 3> pvecCascAsD;
        std::array<float, 3> pvecPionFromOmegac;
        df.propagateTracksToVertex();
        if (!df.isPropagateTracksToVertexDone()) {
          continue;
        }
        df.getTrack(0).getPxPyPzGlo(pvecCascAsD);
        df.getTrack(1).getPxPyPzGlo(pvecPionFromOmegac);
        std::array<float, 3> pvecOmegac = {pvecCascAsD[0] + pvecPionFromOmegac[0], pvecCascAsD[1] + pvecPionFromOmegac[1], pvecCascAsD[2] + pvecPionFromOmegac[2]};

        std::array<float, 3> coordVtxOmegac = df.getPCACandidatePos();
        std::array<float, 6> covVtxOmegac = df.calcPCACovMatrixFlat();

        // create omegac track
        auto trackOmegac = o2::dataformats::V0(coordVtxOmegac, pvecOmegac, covVtxOmegac, trackCasc, trackParVarPi, {0, 0}, {0, 0});

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

        hxVertexOmegac->Fill(vertexOmegacFromFitter[0]);

        // primary pi  pT spectrum
        double ptPrimaryPi = std::sqrt((pvecPionFromOmegac[0] * pvecPionFromOmegac[0]) + (pvecPionFromOmegac[1] * pvecPionFromOmegac[1]));
        hPtPrimaryPi->Fill(ptPrimaryPi);

        // computing invariant mass under the hypothesis of particles ID corresponding to the decay chain
        double mLambda = v0Element.mLambda();         // from LF table, V0 mass under lambda hypothesis
        double mAntiLambda = v0Element.mAntiLambda(); // from LF table, V0 mass under anti-lambda hypothesis

        double myMLambda = 0.;
        const std::array<double, 2> arrMassLambda = {massProtonFromPDG, massPionFromPDG};
        const std::array<double, 2> arrMassAntiLambda = {massPionFromPDG, massProtonFromPDG};
        if (trackXiDauCharged.sign() > 0) {
          myMLambda = RecoDecay::m(std::array{pvecV0Dau0, pvecV0Dau1}, arrMassAntiLambda);
        } else if (trackXiDauCharged.sign() < 0) {
          myMLambda = RecoDecay::m(std::array{pvecV0Dau0, pvecV0Dau1}, arrMassLambda);
        }

        const std::array<double, 2> arrMassCascade = {massLambdaFromPDG, massPionFromPDG};
        double mCascade = RecoDecay::m(std::array{pvecV0AsD, pvecPionFromCasc}, arrMassCascade);
        double mCascadeNotFixed = RecoDecay::m(std::array{pvecV0AsD, pvecPionFromCasc}, std::array{myMLambda, massPionFromPDG});
        double mCascLF = casc.mXi();

        const std::array<double, 2> arrMassOmegac = {massXiFromPDG, massPionFromPDG};
        double mOmegac = RecoDecay::m(std::array{pvecCascAsD, pvecPionFromOmegac}, arrMassOmegac);
        double mOmegacNotFixed = RecoDecay::m(std::array{pvecCascAsD, pvecPionFromOmegac}, std::array{mCascadeNotFixed, massPionFromPDG});

        // computing cosPA
        double cpaV0 = RecoDecay::cpa(coordVtxCasc, coordVtxV0, pvecV0AsM);
        double cpaOmegac = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac, pvecOmegac);
        double cpaCasc = RecoDecay::cpa(coordVtxOmegac, coordVtxCasc, pvecCascAsD);
        double cpaxyV0 = RecoDecay::cpaXY(coordVtxCasc, coordVtxV0, pvecV0AsM);
        double cpaxyOmegac = RecoDecay::cpaXY(std::array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac, pvecOmegac);
        double cpaxyCasc = RecoDecay::cpaXY(coordVtxOmegac, coordVtxCasc, pvecCascAsD);

        // computing decay length and ctau
        double decLenOmegac = RecoDecay::distance(std::array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac);
        double decLenCascade = RecoDecay::distance(coordVtxOmegac, coordVtxCasc);
        double decLenV0 = RecoDecay::distance(coordVtxCasc, coordVtxV0);
        double ctOmegac = RecoDecay::ct(pvecOmegac, decLenOmegac, massOmegacFromPDG);
        double ctCascade = RecoDecay::ct(pvecCascAsD, decLenCascade, massXiFromPDG);
        double ctV0 = RecoDecay::ct(pvecV0AsM, decLenV0, massLambdaFromPDG);

        // computing eta
        double pseudorapOmegac = RecoDecay::eta(pvecOmegac);
        double pseudorapCascade = RecoDecay::eta(pvecCascAsD);
        double pseudorapV0 = RecoDecay::eta(pvecV0AsM);

        // DCA between cascade daughters (from LF table)
        double dcaCascDau = dfc.getChi2AtPCACandidate();
        double dcaV0Dau = dfv.getChi2AtPCACandidate();
        double dcaOmegacDau = df.getChi2AtPCACandidate();

        // set hfFlag
        int hfFlag = 1 << DecayType::OmegacToXiPi;

        // fill test histograms

        hInvMassOmegac->Fill(mOmegac);
        hMassOmegacNotFixed->Fill(mOmegacNotFixed);

        // fill the table
        rowCandidate(collision.globalIndex(),
                     collision.posX(), collision.posY(), collision.posZ(),
                     vertexOmegacFromFitter[0], vertexOmegacFromFitter[1], vertexOmegacFromFitter[2],
                     vertexCascFromFitter[0], vertexCascFromFitter[1], vertexCascFromFitter[2],
                     vertexV0FromFitter[0], vertexV0FromFitter[1], vertexV0FromFitter[2],
                     trackXiDauCharged.sign(),
                     chi2PCAOmegac, chi2PCAV0, chi2PCACascade,
                     pvecOmegac[0], pvecOmegac[1], pvecOmegac[2],
                     pvecCascAsD[0], pvecCascAsD[1], pvecCascAsD[2],
                     pvecPionFromOmegac[0], pvecPionFromOmegac[1], pvecPionFromOmegac[2],
                     pvecV0AsD[0], pvecV0AsD[1], pvecV0AsD[2],
                     pvecPionFromCasc[0], pvecPionFromCasc[1], pvecPionFromCasc[2],
                     pvecV0Dau0[0], pvecV0Dau0[1], pvecV0Dau0[2],
                     pvecV0Dau1[0], pvecV0Dau1[1], pvecV0Dau1[2],
                     impactParameterCasc.getY(), impactParameterPrimaryPi.getY(),
                     impactParameterCasc.getZ(), impactParameterPrimaryPi.getZ(),
                     impactParameterV0.getY(), impactParameterV0.getZ(),
                     std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPrimaryPi.getSigmaY2()), std::sqrt(impactParameterV0.getSigmaY2()),
                     v0.globalIndex(),
                     v0Element.posTrackId(), v0Element.negTrackId(),
                     casc.globalIndex(),
                     trackPion.globalIndex(),         // index pi <- omegac
                     trackXiDauCharged.globalIndex(), // index pi <- cascade
                     impactParameterOmegac.getY(), impactParameterOmegac.getZ(),
                     ptPrimaryPi,
                     mLambda, mAntiLambda, mCascade, mOmegac,
                     cpaV0, cpaOmegac, cpaCasc, cpaxyV0, cpaxyOmegac, cpaxyCasc,
                     ctOmegac, ctCascade, ctV0,
                     pseudorapV0PosDau, pseudorapV0NegDau, pseudorapPiFromCas, pseudorapPiFromOme,
                     pseudorapOmegac, pseudorapCascade, pseudorapV0,
                     myMLambda, mCascadeNotFixed, mOmegacNotFixed,
                     vertexCascLFTable[0], vertexCascLFTable[0], vertexCascLFTable[0], mCascLF, dcaxyPrimaryPi, dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascDau,
                     dcaCascDau, dcaV0Dau, dcaOmegacDau, hfFlag);

      } // loop over pions
    }   // loop over candidates
  }     // end of process
};      // end of struct

/// Performs MC matching.
struct HfCandidateCreatorOmegacMc {
  Produces<aod::HfCandOmegacMCRec> rowMCMatchRec;
  Produces<aod::HfCandOmegacMCGen> rowMCMatchGen;

  void init(InitContext const&) {}

  void processDoNoMc(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorOmegacMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandOmegac const& candidates,
                 aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    // int8_t origin = 0; //to be used for prompt/non prompt
    int8_t debug = 0;

    int pdgCodeOmegac0 = pdg::Code::kOmegaC0; // 4332
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
            flag = sign * (1 << DecayType::OmegacToXiPi);
          }
        }
      }
      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Omegac0 decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug);

    } // close loop over candidates

    // Match generated particles.
    for (auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      // origin = 0;
      //  Omegac → Xi pi
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true)) {
        // Match Xi -> lambda pi
        auto cascMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
        // Printf("Checking cascade → lambda pi");
        if (RecoDecay::isMatchedMCGen(particlesMC, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
          // lambda -> p pi
          auto v0MC = particlesMC.rawIteratorAt(cascMC.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(particlesMC, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign)) {
            flag = sign * (1 << DecayType::OmegacToXiPi);
          }
        }
      }
      // rowMCMatchGen(flag, origin);
      rowMCMatchGen(flag);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorOmegacMc, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorOmegac>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorOmegacMc>(cfgc)};
}
