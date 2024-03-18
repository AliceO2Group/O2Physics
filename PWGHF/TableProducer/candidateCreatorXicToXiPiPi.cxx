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
/// \brief Reconstruction of Ξc± → (Ξ∓ → (Λ → p π∓) π∓) π± π± candidates
///
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, Heidelberg University
/// \author Jinjoo Seo <jseo@cern.ch>, Heidelberg University

#ifndef HomogeneousField
#define HomogeneousField
#endif

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
#include "ReconstructionDataFormats/V0.h" // for creating XicPlus track with DCA fitter

#include "Common/Core/trackUtilities.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_xictoxipipi;
using namespace o2::constants::physics;
using namespace o2::framework;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreatorXic {
  Produces<aod::HfCandXicBase> rowCandidateBase;
  Produces<aod::HfCandXicKF> rowCandidateKF;

  // vertexing
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};
  // DCA fitter
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  //  KFParticle
  Configurable<bool> constrainXicPlusToPv{"constrainXicPlusToPv", false, "Constrain XicPlus to PV"};

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
    if ((std::accumulate(doprocessDF.begin(), doprocessDF.end(), 0) == 1) && fillHistograms) {
      hVertexerType->Fill(aod::hf_cand::VertexerType::DCAFitter);
    }
    if ((std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0) == 1) && fillHistograms) {
      hVertexerType->Fill(aod::hf_cand::VertexerType::KfParticle);
    }
  }

  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  using CascFullKF = soa::Join<aod::KFCascDatas, aod::KFCascCovs>;

  template <bool doPvRefit, typename CandType, typename TTracks>
  void runCreatorXicPlusWithDcaFitter(aod::Collisions const& collisions,
                                      CandType const& rowsTrackIndexXicPlus,
                                      CascFull const& cascs,
                                      TTracks const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // 3-prong vertex fitter
    o2::vertexing::DCAFitterN<3> df;
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
      // V0 daughters
      auto trackV0PosDau = casc.template posTrack_as<TTracks>(); // p <- V0 track (positive track) from TTracks table
      auto trackV0NegDau = casc.template negTrack_as<TTracks>(); // pion <- V0 track (negative track) from TTracks table

      // check that particles come from the same collision
      if (rejDiffCollTrack) {
        if (trackV0PosDau.collisionId() != trackV0NegDau.collisionId()) {
          continue;
        }
        if (trackPionFromXi.collisionId() != trackV0PosDau.collisionId()) {
          continue;
        }
      }

      //--------------------------info of V0 track from LF-table---------------------------
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};

      //-----------------------------info of cascade track from Lf-table-------------------
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
      }
      df.setBz(bz);

      // reconstruct the 3-prong secondary vertex
      if (df.process(trackParVar0, trackParVar1, trackCasc) == 0) {
        continue;
      }

      // set hfFlag
      int hfFlag = BIT(aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi);
      int signXicPlus = casc.sign() < 0 ? +1 : -1;

      // get SV properties
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2SV = df.getChi2AtPCACandidate();
      auto covMatrixSV = df.calcPCACovMatrixFlat();

      // get track momenta
      trackParVar0 = df.getTrack(0);
      trackParVar1 = df.getTrack(1);
      trackCasc = df.getTrack(2);
      std::array<float, 3> pVecXi;
      std::array<float, 3> pVecPi0;
      std::array<float, 3> pVecPi1;
      trackCasc.getPxPyPzGlo(pVecXi);
      trackParVar0.getPxPyPzGlo(pVecPi0);
      trackParVar1.getPxPyPzGlo(pVecPi1);

      // get parent track
      std::array<float, 3> pVec2Pi = RecoDecay::pVec(pVecPi0, pVecPi1);
      std::array<float, 3> pVecXicPlus = RecoDecay::pVec(pVecXi, pVecPi0, pVecPi1);
      auto trackParCov2Pi = o2::dataformats::V0(df.getPCACandidatePos(), pVec2Pi, df.calcPCACovMatrixFlat(), trackParVar0, trackParVar1);
      auto trackParCovXicPlus = o2::dataformats::V0(df.getPCACandidatePos(), pVecXicPlus, df.calcPCACovMatrixFlat(), trackParCov2Pi, trackCasc);
      trackParCovXicPlus.getPxPyPzGlo(pVecXicPlus);

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      if constexpr (doPvRefit) {
        /// use PV refit
        /// Using it in the rowCandidateBase all dynamic columns shall take it into account
        // coordinates
        primaryVertex.setX(rowTrackIndexXicPlus.pvRefitX());
        primaryVertex.setY(rowTrackIndexXicPlus.pvRefitY());
        primaryVertex.setZ(rowTrackIndexXicPlus.pvRefitZ());
        // covariance matrix
        primaryVertex.setSigmaX2(rowTrackIndexXicPlus.pvRefitSigmaX2());
        primaryVertex.setSigmaXY(rowTrackIndexXicPlus.pvRefitSigmaXY());
        primaryVertex.setSigmaY2(rowTrackIndexXicPlus.pvRefitSigmaY2());
        primaryVertex.setSigmaXZ(rowTrackIndexXicPlus.pvRefitSigmaXZ());
        primaryVertex.setSigmaYZ(rowTrackIndexXicPlus.pvRefitSigmaYZ());
        primaryVertex.setSigmaZ2(rowTrackIndexXicPlus.pvRefitSigmaZ2());
        covMatrixPV = primaryVertex.getCov();
      }
      // calculate impact parameter
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      o2::dataformats::DCA impactParameterCasc;
      trackCasc.propagateToDCA(primaryVertex, bz, &impactParameterCasc);
      trackParVar0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParVar1.propagateToDCA(primaryVertex, bz, &impactParameter1);

      // calculate cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      double cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      double cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      double cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      double cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);

      // get invariant mass
      auto arrayMomenta = std::array{pVecXi, pVecPi0, pVecPi1};
      massXiPiPi = RecoDecay::m(std::move(arrayMomenta), std::array{massXiMinusFromPdg, massPionFromPdg, massPionFromPdg});

      // get uncertainty of the decay length
      double phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      // fill histograms
      if (fillHistograms) {
        // invariant mass
        hMass3->Fill(massXiPiPi);
        // covariance matrix elements of PV
        hCovPVXX->Fill(covMatrixPV[0]);
        hCovPVYY->Fill(covMatrixPV[2]);
        hCovPVXZ->Fill(covMatrixPV[3]);
        hCovPVZZ->Fill(covMatrixPV[5]);
        // covariance matrix elements of SV
        hCovSVXX->Fill(covMatrixSV[0]);
        hCovSVYY->Fill(covMatrixSV[2]);
        hCovSVXZ->Fill(covMatrixSV[3]);
        hCovSVZZ->Fill(covMatrixSV[5]);
        // DCAs of prongs
        hDcaXYProngs->Fill(track0.pt(), impactParameter0.getY());
        hDcaXYProngs->Fill(track1.pt(), impactParameter1.getY());
        hDcaXYProngs->Fill(trackCasc.getPt(), impactParameterCasc.getY());
        hDcaZProngs->Fill(track0.pt(), impactParameter0.getZ());
        hDcaZProngs->Fill(track1.pt(), impactParameter1.getZ());
        hDcaZProngs->Fill(trackCasc.getPt(), impactParameterCasc.getZ());
      }

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       covMatrixPV[0], covMatrixPV[2], covMatrixPV[5],
                       /*3-prong specific columns*/
                       rowTrackIndexXicPlus.cascadeId(), rowTrackIndexXicPlus.prong0Id(), rowTrackIndexXicPlus.prong1Id(),
                       casc.bachelorId(), casc.posTrackId(), casc.negTrackId(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       covMatrixSV[0], covMatrixSV[2], covMatrixSV[5],
                       errorDecayLength, errorDecayLengthXY,
                       chi2SV, massXiPiPi, signXicPlus,
                       pVecXi[0], pVecXi[1], pVecXi[2],
                       pVecPi0[0], pVecPi0[1], pVecPi0[2],
                       pVecPi1[0], pVecPi1[1], pVecPi1[2],
                       impactParameterCasc.getY(), impactParameter0.getY(), impactParameter1.getY(), 
                       std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                       hfFlag,
                       /*cascade specific columns*/
                       vertexCasc[0], vertexCasc[1], vertexCasc[2],
                       vertexV0[0], vertexV0[1], vertexV0[2],
                       cpaXi, cpaXYXi, cpaLambda, cpaXYLambda);
    } // loop over track triplets
  }   // template

  template <bool doPvRefit, typename CandType, typename TTracks>
  void runCreatorXicPlusWithKFParticle(aod::Collisions const& collisions,
                                       CandType const& rowsTrackIndexXicPlus,
                                       CascFullKF const&,
                                       TTracks const& tracks,
                                       aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      auto casc = rowTrackIndexXicPlus.template cascade_as<CascFullKF>();
      auto track0 = rowTrackIndexXicPlus.template prong0_as<TTracks>();
      auto track1 = rowTrackIndexXicPlus.template prong1_as<TTracks>();
      auto collision = rowTrackIndexXicPlus.collision();

      //--------------------------info of V0 track from LF-table---------------------------
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};

      //-----------------------------info of cascade track from Lf-table-------------------
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - massXiMinusFromPdg) > massToleranceCascade) {
          continue;
        }
      }

      /// Set the magnetic field from ccdb.
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      KFParticle::SetField(bz);

      // initialize primary vertex
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      float covMatrixPV[6];
      if constexpr (doPvRefit) {
        kfpVertex.SetXYZ(rowTrackIndexXicPlus.pvRefitX(), rowTrackIndexXicPlus.pvRefitY(), rowTrackIndexXicPlus.pvRefitZ());
        kfpVertex.SetCovarianceMatrix(rowTrackIndexXicPlus.pvRefitSigmaX2(), rowTrackIndexXicPlus.pvRefitSigmaXY(), rowTrackIndexXicPlus.pvRefitSigmaY2(), rowTrackIndexXicPlus.pvRefitSigmaXZ(), rowTrackIndexXicPlus.pvRefitSigmaYZ(), rowTrackIndexXicPlus.pvRefitSigmaZ2());
      }
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle KFPV(kfpVertex); // for calculation of DCAs to PV

      // convert pion tracks into KFParticle object
      KFPTrack kfpTrack0 = createKFPTrackFromTrack(track0);
      KFPTrack kfpTrack1 = createKFPTrackFromTrack(track1);
      KFParticle kfPion0(kfpTrack0, kPiPlus);
      KFParticle kfPion1(kfpTrack1, kPiPlus);

      // create Xi as KFParticle object
      // read {X,Y,Z,Px,Py,Pz} and corresponding covariance matrix from KF cascade Tables
      std::array<float, 6> xyzpxpypz = {casc.x(), casc.y(), casc.z(), casc.px(), casc.py(), casc.pz()};
      float parPosMom[6];
      for (int i{0}; i < 6; ++i) {
        parPosMom[i] = xyzpxpypz[i];
      }
      // create KFParticle
      KFParticle kfXi;
      kfXi.Create(parPosMom, casc.kfTrackCovMat(), casc.sign(), casc.m(1)); // casc.m(1): 1 for Lambda mass hypothesis

      // create XicPlus as KFParticle object
      KFParticle kfXicPlus;
      const KFParticle* kfDaughtersXicPlus[3] = {&kfPion0, &kfPion1, &kfXi};
      kfXicPlus.SetConstructMethod(2); // 0: no mass constraint, 2: mass constraint
      kfXicPlus.Construct(kfDaughtersXicPlus, 3);
      // topological constraint
      if (constrainXicPlusToPv) {
        kfXicPlus.SetProductionVertex(KFPV);
      }
      auto covMatrixXicPlus = kfXicPlus.CovarianceMatrix();

      // transport daughter particles to XicPlus decay vertex
      kfPion0.TransportToParticle(kfXicPlus);
      kfPion1.TransportToParticle(kfXicPlus);
      kfXi.TransportToParticle(kfXicPlus);

      // set hfFlag
      int hfFlag = BIT(aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi);
      int signXicPlus = casc.sign() < 0 ? +1 : -1;

      // get impact parameters of XicPlus daughters
      float impactParameter0XY = 0., errImpactParameter0XY = 0.;
      float impactParameter1XY = 0., errImpactParameter1XY = 0.;
      float impactParameterXiXY = 0., errImpactParameterXiXY = 0.;
      kfPion0.GetDistanceFromVertexXY(KFPV, impactParameter0XY, errImpactParameter0XY);
      kfPion1.GetDistanceFromVertexXY(KFPV, impactParameter1XY, errImpactParameter1XY);
      kfXi.GetDistanceFromVertexXY(KFPV, impactParameterXiXY, errImpactParameterXiXY);

      // calculate cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      double cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      double cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      double cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      double cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);

      // get DCAs of Pi0-Pi1, Pi0-Xi, Pi1-Xi
      float dcaXYPi0Pi1 = kfPion0.GetDistanceFromParticleXY(kfPion1);
      float dcaXYPi0Xi = kfPion0.GetDistanceFromParticleXY(kfXi);
      float dcaXYPi1Xi = kfPion1.GetDistanceFromParticleXY(kfXi);

      // fill histograms
      if (fillHistograms) {
        // invariant mass
        hMass3->Fill(kfXicPlus.GetMass());
        // covariance matrix elements of PV
        hCovPVXX->Fill(covMatrixPV[0]);
        hCovPVYY->Fill(covMatrixPV[2]);
        hCovPVXZ->Fill(covMatrixPV[3]);
        hCovPVZZ->Fill(covMatrixPV[5]);
        // covariance matrix elements of SV
        hCovSVXX->Fill(covMatrixXicPlus[0]);
        hCovSVYY->Fill(covMatrixXicPlus[2]);
        hCovSVXZ->Fill(covMatrixXicPlus[3]);
        hCovSVZZ->Fill(covMatrixXicPlus[5]);
        // DCAs of prongs
        hDcaXYProngs->Fill(kfPion0.GetPt(), impactParameter0XY);
        hDcaXYProngs->Fill(kfPion1.GetPt(), impactParameter1XY);
        hDcaXYProngs->Fill(kfXi.GetPt(), impactParameterXiXY);
      }
      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       KFPV.GetX(), KFPV.GetY(), KFPV.GetZ(),
                       covMatrixPV[0], covMatrixPV[2],covMatrixPV[5],
                       /*3-prong specific columns*/
                       rowTrackIndexXicPlus.cascadeId(), rowTrackIndexXicPlus.prong0Id(), rowTrackIndexXicPlus.prong1Id(),
                       casc.bachelorId(), casc.posTrackId(), casc.negTrackId(),
                       kfXicPlus.GetX(), kfXicPlus.GetY(), kfXicPlus.GetZ(),
                       kfXicPlus.GetErrX(), kfXicPlus.GetErrY(), kfXicPlus.GetErrZ(),
                       kfXicPlus.GetErrDecayLength(), kfXicPlus.GetErrDecayLengthXY(),
                       kfXicPlus.GetChi2(), kfXicPlus.GetMass(), signXicPlus,
                       kfXi.GetPx(), kfXi.GetPy(), kfXi.GetPz(),
                       kfPion0.GetPx(), kfPion0.GetPy(), kfPion0.GetPz(),
                       kfPion1.GetPx(), kfPion1.GetPy(), kfPion1.GetPz(),
                       impactParameterXiXY, impactParameter0XY, impactParameter1XY,
                       errImpactParameter0XY, errImpactParameter1XY, errImpactParameterXiXY,
                       hfFlag,
                       /*cascade specific columns*/
                       casc.x(), casc.y(), casc.z(),
                       casc.xlambda(), casc.ylambda(), casc.zlambda(),
                       cpaXi, cpaXYXi, cpaLambda, cpaXYLambda);
      rowCandidateKF(casc.kfV0Chi2(), casc.kfCascadeChi2(),
                     dcaXYPi0Pi1, dcaXYPi0Xi, dcaXYPi1Xi);
    } // loop over track triplets
  }   // template

  void processPvRefitWithDCAFitterN(aod::Collisions const& collisions,
                                    soa::Join<aod::HfCascLf3Prongs, aod::HfPvRefit3Prong> const& rowsTrackIndexXicPlus,
                                    CascFull const& cascs,
                                    aod::TracksWCovDca const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorXicPlusWithDcaFitter<true>(collisions, rowsTrackIndexXicPlus, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic, processPvRefitWithDCAFitterN, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithDCAFitterN(aod::Collisions const& collisions,
                                      aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                      CascFull const& cascs,
                                      aod::TracksWCovDca const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorXicPlusWithDcaFitter<false>(collisions, rowsTrackIndexXicPlus, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic, processNoPvRefitWithDCAFitterN, "Run candidate creator without PV refit", true);

  void processPvRefitWithKFParticle(aod::Collisions const& collisions,
                                    soa::Join<aod::HfCascLf3Prongs, aod::HfPvRefit3Prong> const& rowsTrackIndexXicPlus,
                                    CascFullKF const& cascs,
                                    aod::TracksWCovExtra const& tracks,
                                    aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorXicPlusWithKFParticle<true>(collisions, rowsTrackIndexXicPlus, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic, processPvRefitWithKFParticle, "Run candidate creator with PV refit", false);

  void processNoPvRefitWithKFParticle(aod::Collisions const& collisions,
                                      aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                      CascFullKF const& cascs,
                                      aod::TracksWCovExtra const& tracks,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    runCreatorXicPlusWithKFParticle<false>(collisions, rowsTrackIndexXicPlus, cascs, tracks, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic, processNoPvRefitWithKFParticle, "Run candidate creator without PV refit", false);
}; // struct

/// Performs MC matching.
struct HfCandidateCreatorXicMc {
  Spawns<aod::HfCandXicExt> rowCandXicToXiPiPi;
  Produces<aod::HfCandXicMcRec> rowMcMatchRec;
  Produces<aod::HfCandXicMcGen> rowMcMatchGen;

  Configurable<bool> matchXicPlusMc{"matchXicPlusMc", true, "Do MC matching for XicPlus"};

  void init(InitContext const&) {}

  void processMc(aod::HfCandXic const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    int indexRec = -1;
    int indexRecXicPlus = -1;
    int8_t sign = -9;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;
    int8_t debugGen = 0;

    int pdgCodeXicPlus = Pdg::kXiCPlus; // 4232
    int pdgCodeXiMinus = kXiMinus;      // 3312
    int pdgCodeLambda = kLambda0;       // 3122
    int pdgCodePiPlus = kPiPlus;        // 211
    int pdgCodePiMinus = kPiMinus;      // -211
    int pdgCodeProton = kProton;        // 2212

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      sign = -9;
      origin = RecoDecay::OriginType::None;
      debug = 0;

      auto arrayDaughters = std::array{candidate.pi0_as<aod::TracksWMc>(), // pi <- Xic
                                       candidate.pi1_as<aod::TracksWMc>(), // pi <- Xic
                                       candidate.bachelor_as<aod::TracksWMc>(),       // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),       // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()};      // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Xic → pi pi pi pi p
      indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeXicPlus, std::array{pdgCodePiPlus, pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
      indexRecXicPlus = indexRec;
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
            flag = sign * (1 << aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi);
          }
        }
      }

      // Check whether the charm baryon is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRecXicPlus);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
      }

      rowMcMatchRec(flag, debug, origin);
    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = 0;
      sign = -9;
      debugGen = 0;
      origin = RecoDecay::OriginType::None;

      //  Xic → Xi pi pi
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeXicPlus, std::array{pdgCodeXiMinus, pdgCodePiPlus, pdgCodePiPlus}, true, &sign)) {
        debugGen = 1;
        // Xi- -> Lambda pi
        auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
        if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
          debugGen = 2;
          // Lambda -> p pi
          auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
            debugGen = 3;
            flag = sign * (1 << aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi);
          }
        }
      }

      // Check whether the charm baryon is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
      }

      rowMcMatchGen(flag, debugGen, origin);
    } // close loop over generated particles
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorXicMc, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXicMc>(cfgc)};
}
