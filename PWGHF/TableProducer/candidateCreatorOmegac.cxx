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
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_sel_collision;
using namespace o2::aod::cascdata;
using namespace o2::aod::cascdataext;
using namespace o2::aod::hf_cand_omegac;
using namespace o2::framework::expressions;

// Reconstruction of omegac candidates
struct HfCandidateCreatorOmegac {
  Produces<aod::HfCandOmegacBase> rowCandidateBase;

  Configurable<bool> PropToPCA{"PropToPCA", false, "create tracks version propagated to PCA"};
  Configurable<double> MaxR{"MaxR", 200., "reject PCA's above this radius"};
  Configurable<double> MaxDzIni{"MaxDzIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> MinParamChange{"MinParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> MinRelChi2Change{"MinRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> UseAbsDCA{"UseAbsDCA", true, "use absolute DCA for vertexing"};
  Configurable<bool> UseWeightedPCA{"UseWeightedPCA", false, "vertices use cov matrices"};
  Configurable<bool> RefitWithMatCorr{"RefitWithMatCorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};

  Configurable<bool> RejDiffCollTrack{"RejDiffCollTrack", true, "Reject tracks coming from different collisions"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
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
    df.setPropagateToPCA(PropToPCA);
    df.setMaxR(MaxR);
    df.setMaxDZIni(MaxDzIni);
    df.setMinParamChange(MinParamChange);
    df.setMinRelChi2Change(MinRelChi2Change);
    df.setUseAbsDCA(UseAbsDCA);
    df.setWeightedFinalPCA(UseWeightedPCA);
    df.setRefitWithMatCorr(RefitWithMatCorr);

    // 2-prong vertex fitter to build the cascade vertex
    o2::vertexing::DCAFitterN<2> dfc;
    dfc.setBz(magneticField);
    dfc.setPropagateToPCA(PropToPCA);
    dfc.setMaxR(MaxR);
    dfc.setMaxDZIni(MaxDzIni);
    dfc.setMinParamChange(MinParamChange);
    dfc.setMinRelChi2Change(MinRelChi2Change);
    dfc.setUseAbsDCA(UseAbsDCA);
    dfc.setWeightedFinalPCA(UseWeightedPCA);
    dfc.setRefitWithMatCorr(RefitWithMatCorr);

    // 2-prong vertex fitter to build the V0 vertex
    o2::vertexing::DCAFitterN<2> dfv;
    dfv.setBz(magneticField);
    dfv.setPropagateToPCA(PropToPCA);
    dfv.setMaxR(MaxR);
    dfv.setMaxDZIni(MaxDzIni);
    dfv.setMinParamChange(MinParamChange);
    dfv.setMinRelChi2Change(MinRelChi2Change);
    dfv.setUseAbsDCA(UseAbsDCA);
    dfv.setWeightedFinalPCA(UseWeightedPCA);
    dfv.setRefitWithMatCorr(RefitWithMatCorr);

    // loop over cascades reconstructed by cascadebuilder.cxx
    for (auto const& casc : cascades) {

      if (collision.globalIndex() != casc.collisionId()) {
        continue;
      }

      //----------------accessing particles in the decay chain-------------
      // cascade daughter - charged particle
      int indexTrackXiDauCharged = casc.bachelorId();        // pion <- xi index from cascade table (not used)
      auto trackXiDauCharged = casc.bachelor_as<MyTracks>(); // pion <- xi track from MyTracks table
      // cascade daughter - V0
      int indexV0 = casc.v0Id();                        // VO index from cascades table
      if (!casc.v0_as<aod::V0sLinked>().has_v0Data()) { // check that V0 data are stored
        continue;
      }
      auto v0 = casc.v0_as<aod::V0sLinked>();
      auto v0Element = v0.v0Data(); // V0 element from LF table containing V0 info
      // V0 positive daughter
      auto trackV0Dau0 = v0Element.posTrack_as<MyTracks>(); // p <- V0 track (positive track) from MyTracks table
      // V0 negative daughter
      auto trackV0Dau1 = v0Element.negTrack_as<MyTracks>(); // pion <- V0 track (negative track) from MyTracks table

      // check that particles come from the same collision
      if (RejDiffCollTrack) {
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

      // info from LF table
      array<float, 3> pvecV0LFTable = {v0Element.px(), v0Element.py(), v0Element.pz()};
      array<float, 3> vertexV0 = {v0Element.x(), v0Element.y(), v0Element.z()};

      // reconstruct V0 with DCAFitter
      int nfv = dfv.process(trackParCovV0Dau0, trackParCovV0Dau1);
      if (nfv == 0) {
        continue;
      }
      auto vertexV0FromFitter = dfv.getPCACandidate();
      auto chi2PCAV0 = dfv.getChi2AtPCACandidate();
      array<float, 3> pvecV0Dau0;
      array<float, 3> pvecV0Dau1;
      dfv.propagateTracksToVertex();
      if (!dfv.isPropagateTracksToVertexDone()) {
        continue;
      }
      dfv.getTrack(0).getPxPyPzGlo(pvecV0Dau0);
      dfv.getTrack(1).getPxPyPzGlo(pvecV0Dau1);
      array<float, 3> pvecV0_m = {pvecV0Dau0[0] + pvecV0Dau1[0], pvecV0Dau0[1] + pvecV0Dau1[1], pvecV0Dau0[2] + pvecV0Dau1[2]}; // m stands for mother

      std::array<float, 3> coordVtxV0 = dfv.getPCACandidatePos();
      std::array<float, 6> covVtxV0 = dfv.calcPCACovMatrixFlat();
      std::array<float, 3> momVtxV0 = pvecV0_m;

      // create V0 track
      auto trackV0 = o2::dataformats::V0(dfv.getPCACandidatePos(), pvecV0_m, dfv.calcPCACovMatrixFlat(), trackParCovV0Dau0, trackParCovV0Dau1, {0, 0}, {0, 0});
      auto trackV0_copy = trackV0;

      //-----------------------------reconstruct cascade------------------------------
      // pseudorapidity
      double pseudorapPiFromCas = trackXiDauCharged.eta();

      // info from LF table
      array<float, 3> vertexCascLFTable = {casc.x(), casc.y(), casc.z()};

      // pion <- casc track to be processed with DCAfitter
      auto trackParVarXiDauCharged = getTrackParCov(trackXiDauCharged);

      // reconstruct cascade with DCAFitter
      int nfc = dfc.process(trackV0, trackParVarXiDauCharged);
      if (nfc == 0) {
        continue;
      }
      auto vertexCascFromFitter = dfc.getPCACandidate();
      auto chi2PCACascade = dfc.getChi2AtPCACandidate();
      array<float, 3> pvecV0_d; // d stands for daughter
      array<float, 3> pvecPionFromCasc;
      dfc.propagateTracksToVertex();
      if (!dfc.isPropagateTracksToVertexDone()) {
        continue;
      }
      dfc.getTrack(0).getPxPyPzGlo(pvecV0_d);
      dfc.getTrack(1).getPxPyPzGlo(pvecPionFromCasc);
      array<float, 3> pvecCasc_m = {pvecV0_d[0] + pvecPionFromCasc[0], pvecV0_d[1] + pvecPionFromCasc[1], pvecV0_d[2] + pvecPionFromCasc[2]};

      std::array<float, 3> coordVtxCasc = dfc.getPCACandidatePos();
      std::array<float, 6> covVtxCasc = dfc.calcPCACovMatrixFlat();
      std::array<float, 3> momVtxCasc = pvecCasc_m;

      // create cascade track
      auto trackCasc = o2::dataformats::V0(dfc.getPCACandidatePos(), pvecCasc_m, dfc.calcPCACovMatrixFlat(), trackV0, trackParVarXiDauCharged, {0, 0}, {0, 0});
      auto trackCasc_copy = trackCasc;

      //-------------------combining cascade and pion tracks--------------------------
      for (auto const& trackPion : tracks) {

        if ((RejDiffCollTrack) && (trackXiDauCharged.collisionId() != trackPion.collisionId())) {
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
        auto trackParVarPi_copy = trackParVarPi;

        // reconstruct omegac with DCAFitter
        int nfo = df.process(trackCasc, trackParVarPi);
        if (nfo == 0) {
          continue;
        }
        auto vertexOmegacFromFitter = df.getPCACandidate();
        auto chi2PCAOmegac = df.getChi2AtPCACandidate();
        auto covMatrixPCA = df.calcPCACovMatrixFlat();
        array<float, 3> pvecCasc_d;
        array<float, 3> pvecPionFromOmegac;
        df.propagateTracksToVertex();
        if (!df.isPropagateTracksToVertexDone()) {
          continue;
        }
        df.getTrack(0).getPxPyPzGlo(pvecCasc_d);
        df.getTrack(1).getPxPyPzGlo(pvecPionFromOmegac);
        array<float, 3> pvecOmegac = {pvecCasc_d[0] + pvecPionFromOmegac[0], pvecCasc_d[1] + pvecPionFromOmegac[1], pvecCasc_d[2] + pvecPionFromOmegac[2]};

        std::array<float, 3> coordVtxOmegac = df.getPCACandidatePos();
        std::array<float, 6> covVtxOmegac = df.calcPCACovMatrixFlat();
        std::array<float, 3> momVtxOmegac = pvecOmegac;

        // create omegac track
        auto trackOmegac = o2::dataformats::V0(df.getPCACandidatePos(), pvecOmegac, df.calcPCACovMatrixFlat(), trackCasc, trackParVarPi, {0, 0}, {0, 0});

        // impact parameter omegac
        o2::dataformats::DCA impactParameterOmegac;
        auto primaryVertex = getPrimaryVertex(collision);
        trackOmegac.propagateToDCA(primaryVertex, magneticField, &impactParameterOmegac);

        // impact parameter
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameterCasc;
        o2::dataformats::DCA impactParameterPrimaryPi;
        o2::dataformats::DCA impactParameterV0;
        trackCasc_copy.propagateToDCA(primaryVertex, magneticField, &impactParameterCasc);
        trackParVarPi_copy.propagateToDCA(primaryVertex, magneticField, &impactParameterPrimaryPi);
        trackV0_copy.propagateToDCA(primaryVertex, magneticField, &impactParameterV0);

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

        double my_mLambda = 0.;
        const std::array<float, 2> arrMassLambda = {0.93827, 0.13957};
        const std::array<float, 2> arrMassAntiLambda = {0.13957, 0.93827};
        if (trackXiDauCharged.sign() > 0) {
          my_mLambda = RecoDecay::m(array{pvecV0Dau0, pvecV0Dau1}, arrMassAntiLambda);
        } else if (trackXiDauCharged.sign() < 0) {
          my_mLambda = RecoDecay::m(array{pvecV0Dau0, pvecV0Dau1}, arrMassLambda);
        }

        const std::array<float, 2> arrMassCascade = {1.11568, 0.13957};
        double mCascade = RecoDecay::m(array{pvecV0_d, pvecPionFromCasc}, arrMassCascade);
        double mCascadeNotFixed = RecoDecay::m(array{pvecV0_d, pvecPionFromCasc}, array{my_mLambda, 0.13957});
        double mCascLF = casc.mXi();

        const std::array<float, 2> arrMassOmegac = {1.32171, 0.13957};
        double mOmegac = RecoDecay::m(array{pvecCasc_d, pvecPionFromOmegac}, arrMassOmegac);
        double mOmegacNotFixed = RecoDecay::m(array{pvecCasc_d, pvecPionFromOmegac}, array{mCascadeNotFixed, 0.13957});

        // computing cosPA
        double cpaV0 = RecoDecay::cpa(coordVtxCasc, coordVtxV0, pvecV0_m);
        double cpaOmegac = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac, pvecOmegac);
        double cpaCasc = RecoDecay::cpa(coordVtxOmegac, coordVtxCasc, pvecCasc_d);
        double cpaxyV0 = RecoDecay::cpaXY(coordVtxCasc, coordVtxV0, pvecV0_m);
        double cpaxyOmegac = RecoDecay::cpaXY(array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac, pvecOmegac);
        double cpaxyCasc = RecoDecay::cpaXY(coordVtxOmegac, coordVtxCasc, pvecCasc_d);

        // computing decay length and ctau
        double declenOmegac = RecoDecay::distance(array{collision.posX(), collision.posY(), collision.posZ()}, coordVtxOmegac);
        double declenCascade = RecoDecay::distance(coordVtxOmegac, coordVtxCasc);
        double declenV0 = RecoDecay::distance(coordVtxCasc, coordVtxV0);
        double ctOmegac = RecoDecay::ct(pvecOmegac, declenOmegac, 2.6952);
        double ctCascade = RecoDecay::ct(pvecCasc_d, declenCascade, 1.32171);
        double ctV0 = RecoDecay::ct(pvecV0_m, declenV0, 1.11568);

        // computing eta
        double pseudorapOmegac = RecoDecay::eta(pvecOmegac);
        double pseudorapCascade = RecoDecay::eta(pvecCasc_d);
        double pseudorapV0 = RecoDecay::eta(pvecV0_m);

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
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         vertexOmegacFromFitter[0], vertexOmegacFromFitter[1], vertexOmegacFromFitter[2],
                         vertexCascFromFitter[0], vertexCascFromFitter[1], vertexCascFromFitter[2],
                         vertexV0FromFitter[0], vertexV0FromFitter[1], vertexV0FromFitter[2],
                         trackXiDauCharged.sign(),
                         chi2PCAOmegac, chi2PCAV0, chi2PCACascade,
                         pvecOmegac[0], pvecOmegac[1], pvecOmegac[2],
                         pvecCasc_d[0], pvecCasc_d[1], pvecCasc_d[2],
                         pvecPionFromOmegac[0], pvecPionFromOmegac[1], pvecPionFromOmegac[2],
                         pvecV0_d[0], pvecV0_d[1], pvecV0_d[2],
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
                         my_mLambda, mCascadeNotFixed, mOmegacNotFixed,
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

  void processMC(aod::HfCandOmegacBase const& candidates,
                 aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    // int8_t origin = 0; //to be used for prompt/non prompt
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      // Printf("New rec. candidate");
      flag = 0;
      // origin = 0;
      debug = 0;
      auto arrayDaughters = array{candidate.primaryPi_as<aod::BigTracksMC>(), // pi <- omegac
                                  candidate.bachelor_as<aod::BigTracksMC>(),  // pi <- cascade
                                  candidate.posTrack_as<aod::BigTracksMC>(),  // p <- lambda
                                  candidate.negTrack_as<aod::BigTracksMC>()}; // pi <- lambda
      auto arrayDaughtersCasc = array{candidate.bachelor_as<aod::BigTracksMC>(),
                                      candidate.posTrack_as<aod::BigTracksMC>(),
                                      candidate.negTrack_as<aod::BigTracksMC>()};
      auto arrayDaughtersV0 = array{candidate.posTrack_as<aod::BigTracksMC>(),
                                    candidate.negTrack_as<aod::BigTracksMC>()};

      // Omegac → p pi pi pi
      // Printf("Checking Omegac → p π π π");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 4332, array{211, -211, 2212, -211}, true, &sign, 3); // pdg::Code::kOmegac0=4332 - proton 2212 - pi+ 211
      if (indexRec == -1) {
        debug = 1;
      }
      if (indexRec > -1) {
        // cascade → lambda pi
        // Printf("Checking cascade → lambda pi");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersCasc, 3312, array{-211, 2212, -211}, true, &sign, 2); // pdg::Code::kXiMinus=3312
        if (indexRec == -1) {
          debug = 2;
        }
        if (indexRec > -1) {
          // v0 → p pi
          // Printf("Checking v0 → p pi");
          indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersV0, 3122, array{2212, -211}, true, &sign, 1); // pdg::Code::kLambda=3122
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
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, 4332, array{3312, 211}, true)) {
        // Match Xi -> lambda pi
        auto cascMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
        // Printf("Checking cascade → lambda pi");
        if (RecoDecay::isMatchedMCGen(particlesMC, cascMC, 3312, array{3122, -211}, true)) {
          // lambda -> p pi
          auto v0MC = particlesMC.rawIteratorAt(cascMC.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(particlesMC, v0MC, 3122, array{2212, -211}, true, &sign)) {
            flag = sign * (1 << DecayType::OmegacToXiPi);
          }
        }
      }
      // rowMCMatchGen(flag, origin);
      rowMCMatchGen(flag);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorOmegacMc, processMC, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorOmegac>(cfgc, TaskName{"hf-candidate-creator-omegac"}),
    adaptAnalysisTask<HfCandidateCreatorOmegacMc>(cfgc, TaskName{"hf-candidate-creator-omegac-mc"})};
}
