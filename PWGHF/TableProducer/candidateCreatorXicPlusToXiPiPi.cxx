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

/// \file candidateCreatorXicToXiPiPi.cxx
/// \brief Reconstruction of Ξc± → Ξ∓ π± π± candidates
///
/// \author Jinjoo Seo <jseo@cern.ch>, Heidelberg University
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, CERN

#include <KFParticleBase.h>
#include <KFParticle.h>
#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFVertex.h>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_xicplustoxipipi;
using namespace o2::constants::physics;
using namespace o2::framework;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreatorXicPlus {
  Produces<aod::HfCandXicPlus> rowCandidateBase;

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // cascade cuts
  Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
  Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascade"};
  Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in xy plane"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  double massXiMinusFromPdg = MassXiMinus;
  double massPionFromPdg = MassPiPlus;

  int runNumber{0};
  float toMicrometers = 10000.; // from cm to µm
  double massPi{0.};
  double massXi{0.};
  double massXiPiPi{0.};
  double bz{0.};

  OutputObj<TH1F> hMass3{TH1F("hMass3", "3-prong candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});entries", 500, 2.3, 2.7)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVYY{TH1F("hCovPVYY", "3-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVYY{TH1F("hCovSVYY", "3-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH1F> hCovPVXZ{TH1F("hCovPVXZ", "3-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, -1.e-4, 1.e-4)};
  OutputObj<TH1F> hCovSVXZ{TH1F("hCovSVXZ", "3-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, -1.e-4, 0.2)};
  OutputObj<TH1F> hCovPVZZ{TH1F("hCovPVZZ", "3-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVZZ{TH1F("hCovSVZZ", "3-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};
  OutputObj<TH2F> hDcaXYProngs{TH2F("hDcaXYProngs", "DCAxy of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH2F> hDcaZProngs{TH2F("hDcaZProngs", "DCAz of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", 100, 0., 20., 200, -500., 500.)};
  OutputObj<TH1F> hVertexerType{TH1F("hVertexerType", "Use KF or DCAFitterN;Vertexer type;entries", 2, -0.5, 1.5)}; // See o2::aod::hf_cand::VertexerType

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    std::array<bool, 2> doprocessDF{doprocessPvRefitWithDCAFitterN, doprocessNoPvRefitWithDCAFitterN};
    std::array<bool, 2> doprocessKF{doprocessPvRefitWithKFParticle, doprocessNoPvRefitWithKFParticle};
    if ((std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) + std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }
    if (std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) == 1) {
      hVertexerType->Fill(aod::hf_cand::VertexerType::DCAFitter);
    }
    if (std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0) == 1) {
      hVertexerType->Fill(aod::hf_cand::VertexerType::KfParticle);
    }
  }

  // Will be modified
  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  using V0Full = soa::Join<aod::V0Datas, aod::V0Covs>;

  template <bool doPvRefit, typename CandType, typename TTracks>
  void runCreatorXicPlusWithDcaFitter(aod::Collisions const& collisions,
                                          CandType const& rowsTrackIndexXicPlus,
                                          aod::V0sLinked const& v0Linked,
                                          V0Full const& v0s,
                                          CascFull const& cascs,
                                          TTracks const& tracks,
                                          aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // 3-prong vertex fitter
    o2::vertexing::DCAFitterN<3> df;
    // df.setBz(bz);
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      auto track0 = rowTrackIndexXicPlus.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexXicPlus.template prong1_as<TTracks>();
      auto casc = rowTrackIndexXicPlus.template cascade_as<CascFull>();

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - massXiMinusFromPdg) > massToleranceCascade) {
          continue;
        }
      }

      //----------------accessing particles in the decay chain-------------
      // cascade daughter - charged particle
      auto trackPionFromXi = casc.template bachelor_as<TTracks>(); // pion <- xi track from TTracks table
      // cascade daughter - V0
      auto v0 = casc.template v0_as<aod::V0sLinked>();
      auto v0Element = v0.template v0Data_as<V0Full>(); // V0 element from LF table containing V0 info
      if (!v0.has_v0Data()) { // check that V0 data are stored
          continue;
        }
      
      // V0 positive daughter
      auto trackV0PosDau = v0Element.template posTrack_as<TTracks>(); // p <- V0 track (positive track) from TTracks table
      // V0 negative daughter
      auto trackV0NegDau = v0Element.template negTrack_as<TTracks>(); // pion <- V0 track (negative track) from TTracks table

      // check that particles come from the same collision
      if (rejDiffCollTrack) {
        if (trackV0PosDau.collisionId() != trackV0NegDau.collisionId()) {
          continue;
        }
        if (trackPionFromXi.collisionId() != trackV0PosDau.collisionId()) {
          continue;
        }
      }

      //--------------------------reconstruct V0 track---------------------------
      // pseudorapidity
      // double pseudorapV0PosDau = trackV0PosDau.eta();
      // double pseudorapV0NegDau = trackV0NegDau.eta();

      // pion & p <- V0 tracks
      auto trackParCovV0PosDau = getTrackParCov(trackV0PosDau);
      auto trackParCovV0NegDau = getTrackParCov(trackV0NegDau);

      // info from LF table
      // std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0PosDau = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0NegDau = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-----------------------------reconstruct cascade track-----------------------------
      // pseudorapidity
      // double pseudorapPionFromXi = trackPionFromXi.eta();

      // pion <- casc track to be processed with DCAfitter
      auto trackParCovPionFromXi = getTrackParCov(trackPionFromXi);

      // info from LF table
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
        covCasc[i] = casc.positionCovMat()[i];
      }

      // create cascade track
      o2::track::TrackParCov trackCasc;
      if (trackPionFromXi.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackPionFromXi.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }

      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      std::array<float, 3> pVecPionFromXi = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto collision = rowTrackIndexXicPlus.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      df.setBz(bz);

      // reconstruct the 3-prong secondary vertex
      if (df.process(trackParVar0, trackParVar1, trackCasc) == 0) {
        continue;
      }

      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2SV = df.getChi2AtPCACandidate();
      auto covMatrixSV = df.calcPCACovMatrixFlat();
      hCovSVXX->Fill(covMatrixSV[0]);
      hCovSVYY->Fill(covMatrixSV[2]);
      hCovSVXZ->Fill(covMatrixSV[3]);
      hCovSVZZ->Fill(covMatrixSV[5]);
      trackParVar0 = df.getTrack(0);
      trackParVar1 = df.getTrack(1);
      trackCasc = df.getTrack(2);

      // get track momenta
      std::array<float, 3> pVec0;
      std::array<float, 3> pVec1;
      std::array<float, 3> pVec2;
      trackParVar0.getPxPyPzGlo(pVec0);
      trackParVar1.getPxPyPzGlo(pVec1);
      trackCasc.getPxPyPzGlo(pVec2);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        //primaryVertex.setX(rowTrackIndexXicPlus.pvRefitX());
        //primaryVertex.setY(rowTrackIndexXicPlus.pvRefitY());
        //primaryVertex.setZ(rowTrackIndexXicPlus.pvRefitZ());
        // covariance matrix
        //primaryVertex.setSigmaX2(rowTrackIndexXicPlus.pvRefitSigmaX2());
        //primaryVertex.setSigmaXY(rowTrackIndexXicPlus.pvRefitSigmaXY());
        //primaryVertex.setSigmaY2(rowTrackIndexXicPlus.pvRefitSigmaY2());
        //primaryVertex.setSigmaXZ(rowTrackIndexXicPlus.pvRefitSigmaXZ());
        //primaryVertex.setSigmaYZ(rowTrackIndexXicPlus.pvRefitSigmaYZ());
        //primaryVertex.setSigmaZ2(rowTrackIndexXicPlus.pvRefitSigmaZ2());
        //covMatrixPV = primaryVertex.getCov();
      }
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);

      // calculate impact parameter
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      o2::dataformats::DCA impactParameterCasc;
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);
      trackCasc.propagateToDCA(primaryVertex, bz, &impactParameterCasc);
      hDcaXYProngs->Fill(track0.pt(), impactParameter0.getY());
      hDcaXYProngs->Fill(track1.pt(), impactParameter1.getY());
      //hDcaXYProngs->Fill(trackCasc.pt(), impactParameterCasc.getY());
      hDcaZProngs->Fill(track0.pt(), impactParameter0.getZ());
      hDcaZProngs->Fill(track1.pt(), impactParameter1.getZ());
      //hDcaZProngs->Fill(trackCasc.pt(), impactParameterCasc.getZ());

      // set hfFlag
      int hfFlag = 1 << aod::hf_cand_xicplustoxipipi::DecayType::XicPlusToXiPiPi;

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       rowTrackIndexXicPlus.prong0Id(), rowTrackIndexXicPlus.prong1Id(), rowTrackIndexXicPlus.cascadeId(),
                       trackPionFromXi.globalIndex(),v0Element.globalIndex(), v0Element.posTrackId(), v0Element.negTrackId(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       errorDecayLength, errorDecayLengthXY,
                       chi2SV,
                       pVec0[0], pVec0[1], pVec0[2],
                       pVec1[0], pVec1[1], pVec1[2],
                       pVec2[0], pVec2[1], pVec2[2],
                       impactParameter0.getY(), impactParameter1.getY(), impactParameterCasc.getY(),
                       std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()), std::sqrt(impactParameterCasc.getSigmaY2()),
                       vertexCasc[0], vertexCasc[1], vertexCasc[2],
                       vertexV0[0], vertexV0[1], vertexV0[2],
                       pVecV0PosDau[0], pVecV0PosDau[1], pVecV0PosDau[2],
                       pVecV0NegDau[0], pVecV0NegDau[1], pVecV0NegDau[2],
                       pVecPionFromXi[0], pVecPionFromXi[1], pVecPionFromXi[2],
                       trackV0PosDau.dcaXY(), trackV0NegDau.dcaXY(), trackPionFromXi.dcaXY(),
                       hfFlag);

      // fill histograms
      if (fillHistograms) {
        // calculate invariant mass
        auto arrayMomenta = std::array{pVec0, pVec1, pVec2};
        massXiPiPi = RecoDecay::m(std::move(arrayMomenta), std::array{massXiMinusFromPdg, massPionFromPdg, massPionFromPdg});
        hMass3->Fill(massXiPiPi);
      }
    }
  }

  template <bool doPvRefit, typename CandType, typename TTracks>
  void runCreatorXicPlusWithKFParticle(aod::Collisions const& collisions,
                                          CandType const& rowsTrackIndexXicPlus,
                                          aod::V0sLinked const&,
                                          V0Full const&,
                                          CascFull const&,
                                          TTracks const& tracks,
                                          aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      auto track0 = rowTrackIndexXicPlus.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexXicPlus.template prong1_as<TTracks>();
      auto casc = rowTrackIndexXicPlus.template cascade_as<CascFull>();

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - massXiMinusFromPdg) > massToleranceCascade) {
          continue;
        }
      }

      //----------------accessing particles in the decay chain-------------
      // cascade daughter - charged particle
      auto trackPionFromXi = casc.template bachelor_as<TTracks>(); // pion <- xi track from TTracks table
      // cascade daughter - V0
      auto v0 = casc.template v0_as<aod::V0sLinked>();
      auto v0Element = v0.template v0Data_as<V0Full>(); // V0 element from LF table containing V0 info
      if (!v0.has_v0Data()) { // check that V0 data are stored
          continue;
        }
      // V0 positive daughter
      auto trackV0PosDau = v0Element.template posTrack_as<TTracks>(); // p <- V0 track (positive track) from TTracks table
      // V0 negative daughter
      auto trackV0NegDau = v0Element.template negTrack_as<TTracks>(); // pion <- V0 track (negative track) from TTracks table

      // check that particles come from the same collision
      if (rejDiffCollTrack) {
        if (trackV0PosDau.collisionId() != trackV0NegDau.collisionId()) {
          continue;
        }
        if (trackPionFromXi.collisionId() != trackV0PosDau.collisionId()) {
          continue;
        }
      }

      //--------------------------reconstruct V0 track---------------------------
      // pseudorapidity
      // double pseudorapV0PosDau = trackV0PosDau.eta();
      // double pseudorapV0NegDau = trackV0NegDau.eta();

      // pion & p <- V0 tracks
      auto trackParCovV0PosDau = getTrackParCov(trackV0PosDau);
      auto trackParCovV0NegDau = getTrackParCov(trackV0NegDau);

      // info from LF table
      // std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()}; // pVec stands for vector containing the 3-momentum components
      // std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      //std::array<float, 3> pVecV0PosDau = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      //std::array<float, 3> pVecV0NegDau = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-----------------------------reconstruct cascade track-----------------------------
      // pseudorapidity
      // double pseudorapPionFromXi = trackPionFromXi.eta();

      // pion <- casc track to be processed with DCAfitter
      auto trackParCovPionFromXi = getTrackParCov(trackPionFromXi);

      // info from LF table
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
        covCasc[i] = casc.positionCovMat()[i];
      }

      // create cascade track
      o2::track::TrackParCov trackCasc;
      if (trackPionFromXi.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackPionFromXi.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }

      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      //std::array<float, 3> pVecPionFromXi = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      auto trackParVar0 = getTrackParCov(track0);
      auto trackParVar1 = getTrackParCov(track1);
      auto collision = rowTrackIndexXicPlus.collision();

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
        // df.setBz(bz); /// put it outside the 'if'! Otherwise we have a difference wrt bz Configurable (< 1 permille) in Run2 conv. data
        // df.print();
      }
      float covMatrixPV[6];

      //KFParticle::SetField(bz);
      //KFPVertex kfpVertex = createKFPVertexFromCollision(collision);

      if constexpr (doPvRefit) {
        // to do
      }
      //kfpVertex.GetCovarianceMatrix(covMatrixPV);
      //KFParticle KFPV(kfpVertex);
      hCovPVXX->Fill(covMatrixPV[0]);
      hCovPVYY->Fill(covMatrixPV[2]);
      hCovPVXZ->Fill(covMatrixPV[3]);
      hCovPVZZ->Fill(covMatrixPV[5]);

      /*KFPTrack kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(track1);
      KFPTrack kfpCasc = createKFPTrackFromTrack(trackCasc);*/

      // fill candidate table rows
      //rowCandidateBase();

      // fill histograms
      // if (fillHistograms) {}
    }
  }

  void processPvRefitWithDCAFitterN(aod::Collisions const& collisions,
                                    aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    aod::TracksWCovDca const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorXicPlusWithDcaFitter<true>(collisions, rowsTrackIndexXicPlus, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicPlus, processPvRefitWithDCAFitterN, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithDCAFitterN(aod::Collisions const& collisions,
                                    aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    aod::TracksWCovDca const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorXicPlusWithDcaFitter<false>(collisions, rowsTrackIndexXicPlus, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicPlus, processNoPvRefitWithDCAFitterN, "Run candidate creator without PV refit", true);
  
  void processPvRefitWithKFParticle(aod::Collisions const& collisions,
                                    aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus, 
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    aod::TracksWCovDca const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // runCreatorXicPlusWithKFParticle<true>(collisions, rowsTrackIndexXicPlus, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicPlus, processPvRefitWithKFParticle, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithKFParticle(aod::Collisions const& collisions,
                                    aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                    aod::V0sLinked const& v0Linked,
                                    V0Full const& v0s,
                                    CascFull const& cascs,
                                    aod::TracksWCovDca const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // runCreatorXicPlusWithKFParticle<false>(collisions, rowsTrackIndexXicPlus, v0Linked, v0s, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicPlus, processNoPvRefitWithKFParticle, "Run candidate creator without PV refit", false);
}; //struct

/// Performs MC matching.
struct HfCandidateCreatorXicPlusMc {
  Produces<aod::HfXicPlusMcRec> rowMcMatchRec;
  Produces<aod::HfXicPlusMcGen> rowMcMatchGen;

  Configurable<bool> matchXicPlusMc{"matchXicPlusMc", true, "Do MC matching for XicPlus"};

  void init(InitContext const&) {}

  void processMc(aod::HfCandXicPlus const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    int indexRec = -1;
    int8_t sign = -9;
    int8_t flag = -9;
    int8_t origin =-9;
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenXi = 0;
    int8_t debugGenLambda = 0;

    int pdgCodeXicPlus = Pdg::kXiCPlus;       // 4232
    int pdgCodeXiMinus = kXiMinus;            // 3312
    int pdgCodeLambda = kLambda0;             // 3122
    int pdgCodePiPlus = kPiPlus;              // 211
    int pdgCodePiMinus = kPiMinus;            // -211
    int pdgCodeProton = kProton;              // 2212

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = 0;
      debug = 0;
      auto arrayDaughters = std::array{candidate.pi1FromXicPlus_as<aod::TracksWMc>(), // pi <- Xic
                                       candidate.pi2FromXicPlus_as<aod::TracksWMc>(), // pi <- Xic
                                       candidate.bachelor_as<aod::TracksWMc>(),  // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),  // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Xic matching
      if (matchXicPlusMc) {
        // Xic → pi pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeXicPlus, std::array{pdgCodePiPlus, pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
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
              flag = sign * (1 << aod::hf_cand_xicplustoxipipi::DecayType::XicPlusToXiPiPi);
            }
          }
        }
      }

      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMcMatchRec(flag, debug, origin);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = -9;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenXi = 0;
      debugGenLambda = 0;
      origin = 0;
      if (matchXicPlusMc) {
        //  Xic → Xi pi pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeXicPlus, std::array{pdgCodeXiMinus, pdgCodePiPlus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          // Xi- -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
            debugGenXi = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_xicplustoxipipi::DecayType::XicPlusToXiPiPi);
            }
          }
        }
      }
      rowMcMatchGen(flag, debugGenCharmBar, debugGenXi, debugGenLambda, origin);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorXicPlusMc, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXicPlus>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXicPlusMc>(cfgc)};
}
