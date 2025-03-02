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
/// \author Carolina Reetz <c.reetz@cern.ch>, Heidelberg University
/// \author Jaeyoon Cho <jaeyoon.cho@cern.ch>, Inha University
/// \author Jinjoo Seo <jseo@cern.ch>, Heidelberg University

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include <string>
#include <utility>
#include <vector>

#include <KFParticleBase.h>
#include <KFParticle.h>
#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFVertex.h>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_xic_to_xi_pi_pi;
using namespace o2::constants::physics;
using namespace o2::framework;

/// Reconstruction of heavy-flavour 3-prong decay candidates
struct HfCandidateCreatorXicToXiPiPi {
  Produces<aod::HfCandXicBase> rowCandidateBase;
  Produces<aod::HfCandXicKF> rowCandidateKF;

  Configurable<bool> fillHistograms{"fillHistograms", true, "do validation plots"};
  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  // cascade preselections
  Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
  Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascade"};
  Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in xy plane"};
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
  Configurable<bool> constrainXiToXicPlus{"constrainXiToXicPlus", false, "Constrain Xi to XicPlus"};
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "Construct method of XicPlus: 0 fast mathematics without constraint of fixed daughter particle masses, 2 daughter particle masses stay fixed in construction process"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions (effective only for KFParticle w/o derived data)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  o2::vertexing::DCAFitterN<3> df;

  int runNumber{0};
  float massXiPiPi{0.};
  float massXiPi0{0.};
  float massXiPi1{0.};
  double bz{0.};
  enum XicCandCounter { AllIdTriplets = 0,
                        CascPreSel,
                        VertexFit };

  using CascadesLinked = soa::Join<aod::Cascades, aod::CascDataLink>;
  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  using KFCascadesLinked = soa::Join<aod::Cascades, aod::KFCascDataLink>;
  using KFCascFull = soa::Join<aod::KFCascDatas, aod::KFCascCovs>;
  using TracksWCovDcaPidPrPi = soa::Join<aod::TracksWCovDca, aod::TracksPidPr, aod::TracksPidPi>;
  using TracksWCovExtraPidPrPi = soa::Join<aod::TracksWCovExtra, aod::TracksPidPr, aod::TracksPidPi>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if ((doprocessXicplusWithDcaFitter + doprocessXicplusWithKFParticle) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }

    // add histograms to registry
    if (fillHistograms) {
      // counter
      registry.add("hVertexerType", "Use KF or DCAFitterN;Vertexer type;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}}); // See o2::aod::hf_cand::VertexerType
      registry.add("hCandCounter", "hCandCounter", {HistType::kTH1F, {{3, 0.5, 3.5}}});
      registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + AllIdTriplets, "total");
      registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + CascPreSel, "Cascade preselection");
      registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + VertexFit, "Successful vertex fit");
      // physical variables
      registry.add("hMass3", "3-prong candidates;inv. mass (#Xi #pi #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2.3, 2.7}}});
      registry.add("hCovPVXX", "3-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, 0., 1.e-4}}});
      registry.add("hCovSVXX", "3-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, 0., 0.2}}});
      registry.add("hCovPVYY", "3-prong candidates;YY element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, 0., 1.e-4}}});
      registry.add("hCovSVYY", "3-prong candidates;YY element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, 0., 0.2}}});
      registry.add("hCovPVXZ", "3-prong candidates;XZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, -1.e-4, 1.e-4}}});
      registry.add("hCovSVXZ", "3-prong candidates;XZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, -1.e-4, 0.2}}});
      registry.add("hCovPVZZ", "3-prong candidates;ZZ element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, 0., 1.e-4}}});
      registry.add("hCovSVZZ", "3-prong candidates;ZZ element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1D, {{100, 0., 0.2}}});
      registry.add("hDcaXYProngs", "DCAxy of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{xy}) (#mum);entries", {HistType::kTH2D, {{100, 0., 20.}, {200, -500., 500.}}});
      registry.add("hDcaZProngs", "DCAz of 3-prong candidates;#it{p}_{T} (GeV/#it{c};#it{d}_{z}) (#mum);entries", {HistType::kTH2D, {{100, 0., 20.}, {200, -500., 500.}}});
    }

    // fill hVertexerType histogram
    if (doprocessXicplusWithDcaFitter && fillHistograms) {
      registry.fill(HIST("hVertexerType"), aod::hf_cand::VertexerType::DCAFitter);
    }
    if (doprocessXicplusWithKFParticle && fillHistograms) {
      registry.fill(HIST("hVertexerType"), aod::hf_cand::VertexerType::KfParticle);
    }

    // initialize CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    // initialize 3-prong vertex fitter
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);
  }

  void processXicplusWithDcaFitter(aod::Collisions const&,
                                   aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                   CascadesLinked const&,
                                   CascFull const&,
                                   TracksWCovDcaPidPrPi const&,
                                   aod::BCsWithTimestamps const&)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      auto cascAodElement = rowTrackIndexXicPlus.cascade_as<CascadesLinked>();
      if (!cascAodElement.has_cascData()) {
        continue;
      }
      auto casc = cascAodElement.cascData_as<CascFull>();
      auto trackCharmBachelor0 = rowTrackIndexXicPlus.prong0_as<TracksWCovDcaPidPrPi>();
      auto trackCharmBachelor1 = rowTrackIndexXicPlus.prong1_as<TracksWCovDcaPidPrPi>();
      auto collision = rowTrackIndexXicPlus.collision();
      if (fillHistograms) {
        registry.fill(HIST("hCandCounter"), 1 + AllIdTriplets);
      }

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - MassXiMinus) > massToleranceCascade) {
          continue;
        }
      }
      if (fillHistograms) {
        registry.fill(HIST("hCandCounter"), 1 + CascPreSel);
      }

      //----------------------Set the magnetic field from ccdb---------------------------------------
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      df.setBz(bz);

      //--------------------------info of V0 and cascades track from LF-tables---------------------------
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};

      //----------------create cascade track------------------------------------------------------------
      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
        covCasc[i] = casc.positionCovMat()[i];
      }
      // create cascade track
      o2::track::TrackParCov trackCasc;
      if (casc.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (casc.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }
      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      //----------------------------fit SV and create XicPlus track------------------
      auto trackParCovCharmBachelor0 = getTrackParCov(trackCharmBachelor0);
      auto trackParCovCharmBachelor1 = getTrackParCov(trackCharmBachelor1);

      // reconstruct the 3-prong secondary vertex
      try {
        if (df.process(trackCasc, trackParCovCharmBachelor0, trackParCovCharmBachelor1) == 0) {
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        continue;
      }
      if (fillHistograms) {
        registry.fill(HIST("hCandCounter"), 1 + VertexFit);
      }

      //----------------------------calculate physical properties-----------------------
      // Charge of charm baryon
      int signXic = casc.sign() < 0 ? +1 : -1;

      // get SV properties
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2SV = df.getChi2AtPCACandidate();
      auto covMatrixSV = df.calcPCACovMatrixFlat();

      // get track momenta
      trackCasc = df.getTrack(0);
      trackParCovCharmBachelor0 = df.getTrack(1);
      trackParCovCharmBachelor1 = df.getTrack(2);
      std::array<float, 3> pVecXi;
      std::array<float, 3> pVecPi0;
      std::array<float, 3> pVecPi1;
      trackCasc.getPxPyPzGlo(pVecXi);
      trackParCovCharmBachelor0.getPxPyPzGlo(pVecPi0);
      trackParCovCharmBachelor1.getPxPyPzGlo(pVecPi1);

      // get invariant mass of Xic candidate
      auto arrayMomenta = std::array{pVecXi, pVecPi0, pVecPi1};
      massXiPiPi = RecoDecay::m(std::move(arrayMomenta), std::array{MassXiMinus, MassPiPlus, MassPiPlus});

      // get track impact parameters
      // This modifies track momenta!
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      // calculate impact parameter
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameter0;
      o2::dataformats::DCA impactParameter1;
      trackCasc.propagateToDCA(primaryVertex, bz, &impactParameterCasc);
      trackParCovCharmBachelor0.propagateToDCA(primaryVertex, bz, &impactParameter0);
      trackParCovCharmBachelor1.propagateToDCA(primaryVertex, bz, &impactParameter1);

      // calculate cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      float cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float cpaLambdaToXi = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float cpaXYLambdaToXi = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

      // get invariant mass of Xi-pi pairs
      auto arrayMomentaXiPi0 = std::array{pVecXi, pVecPi0};
      massXiPi0 = RecoDecay::m(std::move(arrayMomentaXiPi0), std::array{MassXiMinus, MassPiPlus});
      auto arrayMomentaXiPi1 = std::array{pVecXi, pVecPi1};
      massXiPi1 = RecoDecay::m(std::move(arrayMomentaXiPi1), std::array{MassXiMinus, MassPiPlus});

      // get uncertainty of the decay length
      float phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      //--------------------- get PID information-----------------------
      float nSigTpcPiFromXicPlus0 = trackCharmBachelor0.tpcNSigmaPi();
      float nSigTofPiFromXicPlus0 = trackCharmBachelor0.tofNSigmaPi();
      float nSigTpcPiFromXicPlus1 = trackCharmBachelor1.tpcNSigmaPi();
      float nSigTofPiFromXicPlus1 = trackCharmBachelor1.tofNSigmaPi();
      // Bachelor pion
      auto trackPionFromXi = casc.bachelor_as<TracksWCovDcaPidPrPi>();
      float nSigTpcBachelorPi = trackPionFromXi.tpcNSigmaPi();
      float nSigTofBachelorPi = trackPionFromXi.tofNSigmaPi();
      // Lambda daughters
      auto trackPosLambdaDaughter = casc.posTrack_as<TracksWCovDcaPidPrPi>();
      auto trackNegLambdaDaughter = casc.negTrack_as<TracksWCovDcaPidPrPi>();
      float pPiFromLambda, pPrFromLambda, nSigTpcPiFromLambda, nSigTofPiFromLambda, nSigTpcPrFromLambda, nSigTofPrFromLambda;
      if (signXic == +1) {
        pPiFromLambda = trackNegLambdaDaughter.p();
        nSigTpcPiFromLambda = trackNegLambdaDaughter.tpcNSigmaPi();
        nSigTofPiFromLambda = trackNegLambdaDaughter.tofNSigmaPi();
        pPrFromLambda = trackPosLambdaDaughter.p();
        nSigTpcPrFromLambda = trackPosLambdaDaughter.tpcNSigmaPr();
        nSigTofPrFromLambda = trackPosLambdaDaughter.tofNSigmaPr();
      } else {
        pPiFromLambda = trackPosLambdaDaughter.p();
        nSigTpcPiFromLambda = trackPosLambdaDaughter.tpcNSigmaPi();
        nSigTofPiFromLambda = trackPosLambdaDaughter.tofNSigmaPi();
        pPrFromLambda = trackNegLambdaDaughter.p();
        nSigTpcPrFromLambda = trackNegLambdaDaughter.tpcNSigmaPr();
        nSigTofPrFromLambda = trackNegLambdaDaughter.tofNSigmaPr();
      }

      //--------------------------------------------fill histograms----------------------------------------------------------------
      if (fillHistograms) {
        // invariant mass
        registry.fill(HIST("hMass3"), massXiPiPi);
        // covariance matrix elements of PV
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
        registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
        registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
        registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);
        // covariance matrix elements of SV
        registry.fill(HIST("hCovSVXX"), covMatrixSV[0]);
        registry.fill(HIST("hCovSVYY"), covMatrixSV[2]);
        registry.fill(HIST("hCovSVXZ"), covMatrixSV[3]);
        registry.fill(HIST("hCovSVZZ"), covMatrixSV[5]);
        // DCAs of prongs
        registry.fill(HIST("hDcaXYProngs"), trackCasc.getPt(), impactParameterCasc.getY());
        registry.fill(HIST("hDcaXYProngs"), trackCharmBachelor0.pt(), impactParameter0.getY());
        registry.fill(HIST("hDcaXYProngs"), trackCharmBachelor1.pt(), impactParameter1.getY());
        registry.fill(HIST("hDcaZProngs"), trackCasc.getPt(), impactParameterCasc.getZ());
        registry.fill(HIST("hDcaZProngs"), trackCharmBachelor0.pt(), impactParameter0.getZ());
        registry.fill(HIST("hDcaZProngs"), trackCharmBachelor1.pt(), impactParameter1.getZ());
      }

      //---------------------------------fill candidate table rows-------------------------------------------------------------------------------------------
      rowCandidateBase(collision.globalIndex(),
                       primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       std::sqrt(covMatrixPV[0]), std::sqrt(covMatrixPV[2]), std::sqrt(covMatrixPV[5]),
                       /*3-prong specific columns*/
                       rowTrackIndexXicPlus.cascadeId(), rowTrackIndexXicPlus.prong0Id(), rowTrackIndexXicPlus.prong1Id(),
                       casc.bachelorId(), casc.posTrackId(), casc.negTrackId(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       std::sqrt(covMatrixSV[0]), std::sqrt(covMatrixSV[2]), std::sqrt(covMatrixSV[5]),
                       errorDecayLength, errorDecayLengthXY,
                       chi2SV, massXiPiPi, signXic,
                       pVecXi[0], pVecXi[1], pVecXi[2],
                       pVecPi0[0], pVecPi0[1], pVecPi0[2],
                       pVecPi1[0], pVecPi1[1], pVecPi1[2],
                       impactParameterCasc.getY(), impactParameter0.getY(), impactParameter1.getY(),
                       std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                       /*cascade specific columns*/
                       trackPionFromXi.p(), pPiFromLambda, pPrFromLambda,
                       cpaXi, cpaXYXi, cpaLambda, cpaXYLambda, cpaLambdaToXi, cpaXYLambdaToXi,
                       casc.mXi(), casc.mLambda(), massXiPi0, massXiPi1,
                       /*DCA information*/
                       casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(),
                       casc.dcaXYCascToPV(), casc.dcaZCascToPV(),
                       /*PID information*/
                       nSigTpcPiFromXicPlus0, nSigTpcPiFromXicPlus1, nSigTpcBachelorPi, nSigTpcPiFromLambda, nSigTpcPrFromLambda,
                       nSigTofPiFromXicPlus0, nSigTofPiFromXicPlus1, nSigTofBachelorPi, nSigTofPiFromLambda, nSigTofPrFromLambda);
    } // loop over track triplets
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processXicplusWithDcaFitter, "Run candidate creator with DCAFitter.", true);

  void processXicplusWithKFParticle(aod::Collisions const&,
                                    aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                    KFCascadesLinked const&,
                                    KFCascFull const&,
                                    TracksWCovExtraPidPrPi const&,
                                    aod::BCsWithTimestamps const&)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      auto cascAodElement = rowTrackIndexXicPlus.cascade_as<aod::KFCascadesLinked>();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<KFCascFull>();
      auto trackCharmBachelor0 = rowTrackIndexXicPlus.prong0_as<TracksWCovExtraPidPrPi>();
      auto trackCharmBachelor1 = rowTrackIndexXicPlus.prong1_as<TracksWCovExtraPidPrPi>();
      auto collision = rowTrackIndexXicPlus.collision();
      if (fillHistograms) {
        registry.fill(HIST("hCandCounter"), 1 + AllIdTriplets);
      }

      //-------------------preselect cascade candidates--------------------------------------
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - MassXiMinus) > massToleranceCascade) {
          continue;
        }
      }
      if (fillHistograms) {
        registry.fill(HIST("hCandCounter"), 1 + CascPreSel);
      }

      //----------------------Set the magnetic field from ccdb-----------------------------
      /// The static instance of the propagator was already modified in the HFTrackIndexSkimCreator,
      /// but this is not true when running on Run2 data/MC already converted into AO2Ds.
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      }
      KFParticle::SetField(bz);

      //----------------------info of V0 and cascade tracks from LF-table------------------
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};

      //----------------------Create XicPlus as KFParticle object-------------------------------------------
      // initialize primary vertex
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      float covMatrixPV[6];
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle kfPv(kfpVertex); // for calculation of DCAs to PV

      // convert pion tracks into KFParticle object
      KFPTrack kfpTrackCharmBachelor0 = createKFPTrackFromTrack(trackCharmBachelor0);
      KFPTrack kfpTrackCharmBachelor1 = createKFPTrackFromTrack(trackCharmBachelor1);
      KFParticle kfCharmBachelor0(kfpTrackCharmBachelor0, kPiPlus);
      KFParticle kfCharmBachelor1(kfpTrackCharmBachelor1, kPiPlus);

      // create Xi as KFParticle object
      // read {X,Y,Z,Px,Py,Pz} and corresponding covariance matrix from KF cascade Tables
      std::array<float, 6> xyzpxpypz = {casc.x(), casc.y(), casc.z(), casc.px(), casc.py(), casc.pz()};
      float parPosMom[6];
      for (int i{0}; i < 6; ++i) {
        parPosMom[i] = xyzpxpypz[i];
      }
      // create KFParticle
      KFParticle kfXi;
      kfXi.Create(parPosMom, casc.kfTrackCovMat(), casc.sign(), casc.mXi());

      // create XicPlus as KFParticle object
      KFParticle kfXicPlus;
      const KFParticle* kfDaughtersXicPlus[3] = {&kfCharmBachelor0, &kfCharmBachelor1, &kfXi};
      kfXicPlus.SetConstructMethod(kfConstructMethod);
      try {
        kfXicPlus.Construct(kfDaughtersXicPlus, 3);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct XicPlus : " << e.what();
        continue;
      }
      if (fillHistograms) {
        registry.fill(HIST("hCandCounter"), 1 + VertexFit);
      }

      // get geometrical chi2 of XicPlus
      float chi2GeoXicPlus = kfXicPlus.GetChi2() / kfXicPlus.GetNDF();

      // topological constraint of Xic to PV
      float chi2topoXicPlusToPVBeforeConstraint = kfXicPlus.GetDeviationFromVertex(kfPv);
      KFParticle kfXicPlusToPV = kfXicPlus;
      kfXicPlusToPV.SetProductionVertex(kfPv);
      float chi2topoXicPlusToPV = kfXicPlusToPV.GetChi2() / kfXicPlusToPV.GetNDF();
      if (constrainXicPlusToPv) {
        kfXicPlus = kfXicPlusToPV;
        kfXicPlus.TransportToDecayVertex();
      }

      // topological constraint of Xi to XicPlus
      float chi2topoXiToXicPlusBeforeConstraint = kfXi.GetDeviationFromVertex(kfXicPlus);
      KFParticle kfXiToXicPlus = kfXi;
      kfXiToXicPlus.SetProductionVertex(kfXicPlus);
      float chi2topoXiToXicPlus = kfXiToXicPlus.GetChi2() / kfXiToXicPlus.GetNDF();
      kfXiToXicPlus.TransportToDecayVertex();
      if (constrainXiToXicPlus) {
        KFParticle kfXicPlusWithXiToXicPlus;
        const KFParticle* kfDaughtersXicPlusWithXiToXicPlus[3] = {&kfCharmBachelor0, &kfCharmBachelor1, &kfXiToXicPlus};
        kfXicPlusWithXiToXicPlus.SetConstructMethod(kfConstructMethod);
        try {
          kfXicPlusWithXiToXicPlus.Construct(kfDaughtersXicPlusWithXiToXicPlus, 3);
        } catch (std::runtime_error& e) {
          LOG(debug) << "Failed to construct XicPlus with Xi connstrained to XicPlus: " << e.what();
          continue;
        }
        kfXicPlus = kfXicPlusWithXiToXicPlus;
      }

      // get covariance matrix of XicPlus
      auto covMatrixXicPlus = kfXicPlus.CovarianceMatrix();

      //---------------------calculate physical parameters of XicPlus candidate----------------------
      // sign of charm baryon
      int signXic = casc.sign() < 0 ? +1 : -1;

      // transport XicPlus daughters to XicPlus decay vertex (secondary vertex)
      float secondaryVertex[3] = {0.};
      secondaryVertex[0] = kfXicPlus.GetX();
      secondaryVertex[1] = kfXicPlus.GetY();
      secondaryVertex[2] = kfXicPlus.GetZ();
      kfXi.TransportToPoint(secondaryVertex);
      kfCharmBachelor0.TransportToPoint(secondaryVertex);
      kfCharmBachelor1.TransportToPoint(secondaryVertex);

      // get impact parameters of XicPlus daughters
      float impactParameterPi0XY = 0., errImpactParameterPi0XY = 0.;
      float impactParameterPi1XY = 0., errImpactParameterPi1XY = 0.;
      float impactParameterXiXY = 0., errImpactParameterXiXY = 0.;
      kfCharmBachelor0.GetDistanceFromVertexXY(kfPv, impactParameterPi0XY, errImpactParameterPi0XY);
      kfCharmBachelor1.GetDistanceFromVertexXY(kfPv, impactParameterPi1XY, errImpactParameterPi1XY);
      kfXi.GetDistanceFromVertexXY(kfPv, impactParameterXiXY, errImpactParameterXiXY);

      // calculate cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      float cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float cpaLambdaToXi = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float cpaXYLambdaToXi = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

      // get DCAs of Pi0-Pi1, Pi0-Xi, Pi1-Xi
      float dcaXYPi0Pi1 = kfCharmBachelor0.GetDistanceFromParticleXY(kfCharmBachelor1);
      float dcaXYPi0Xi = kfCharmBachelor0.GetDistanceFromParticleXY(kfXi);
      float dcaXYPi1Xi = kfCharmBachelor1.GetDistanceFromParticleXY(kfXi);
      float dcaPi0Pi1 = kfCharmBachelor0.GetDistanceFromParticle(kfCharmBachelor1);
      float dcaPi0Xi = kfCharmBachelor0.GetDistanceFromParticle(kfXi);
      float dcaPi1Xi = kfCharmBachelor1.GetDistanceFromParticle(kfXi);

      // mass of Xi-Pi0 pair
      KFParticle kfXiPi0;
      float errMassXiPi0;
      const KFParticle* kfXiResonanceDaughtersPi0[2] = {&kfXi, &kfCharmBachelor0};
      kfXiPi0.SetConstructMethod(kfConstructMethod);
      try {
        kfXiPi0.Construct(kfXiResonanceDaughtersPi0, 2);
      } catch (...) {
        LOG(info) << "Failed to construct Xi(1530) with Pi 0";
      }
      kfXiPi0.GetMass(massXiPi0, errMassXiPi0);

      // mass of Xi-Pi1 pair
      KFParticle kfXiPi1;
      float errMassXiPi1;
      const KFParticle* kfXiResonanceDaughtersPi1[2] = {&kfXi, &kfCharmBachelor1};
      kfXiPi1.SetConstructMethod(kfConstructMethod);
      try {
        kfXiPi1.Construct(kfXiResonanceDaughtersPi1, 2);
      } catch (...) {
        LOG(info) << "Failed to construct Xi(1530) with Pi 1";
      }
      kfXiPi1.GetMass(massXiPi1, errMassXiPi1);

      // get invariant mass of Xic candidate
      float errMassXiPiPi;
      kfXicPlus.GetMass(massXiPiPi, errMassXiPiPi);

      //--------------------- get PID information-----------------------
      float nSigTpcPiFromXicPlus0 = trackCharmBachelor0.tpcNSigmaPi();
      float nSigTofPiFromXicPlus0 = trackCharmBachelor0.tofNSigmaPi();
      float nSigTpcPiFromXicPlus1 = trackCharmBachelor1.tpcNSigmaPi();
      float nSigTofPiFromXicPlus1 = trackCharmBachelor1.tofNSigmaPi();
      // Bachelor pion
      auto trackPionFromXi = casc.bachelor_as<TracksWCovExtraPidPrPi>();
      float nSigTpcBachelorPi = trackPionFromXi.tpcNSigmaPi();
      float nSigTofBachelorPi = trackPionFromXi.tofNSigmaPi();
      // Lambda daughters
      auto trackPosLambdaDaughter = casc.posTrack_as<TracksWCovExtraPidPrPi>();
      auto trackNegLambdaDaughter = casc.negTrack_as<TracksWCovExtraPidPrPi>();
      float pPiFromLambda, pPrFromLambda, nSigTpcPiFromLambda, nSigTofPiFromLambda, nSigTpcPrFromLambda, nSigTofPrFromLambda;
      if (signXic == +1) {
        pPiFromLambda = trackNegLambdaDaughter.p();
        nSigTpcPiFromLambda = trackNegLambdaDaughter.tpcNSigmaPi();
        nSigTofPiFromLambda = trackNegLambdaDaughter.tofNSigmaPi();
        pPrFromLambda = trackPosLambdaDaughter.p();
        nSigTpcPrFromLambda = trackPosLambdaDaughter.tpcNSigmaPr();
        nSigTofPrFromLambda = trackPosLambdaDaughter.tofNSigmaPr();
      } else {
        pPiFromLambda = trackPosLambdaDaughter.p();
        nSigTpcPiFromLambda = trackPosLambdaDaughter.tpcNSigmaPi();
        nSigTofPiFromLambda = trackPosLambdaDaughter.tofNSigmaPi();
        pPrFromLambda = trackNegLambdaDaughter.p();
        nSigTpcPrFromLambda = trackNegLambdaDaughter.tpcNSigmaPr();
        nSigTofPrFromLambda = trackNegLambdaDaughter.tofNSigmaPr();
      }

      //-------------------------------fill histograms--------------------------------------------
      if (fillHistograms) {
        // invariant mass
        registry.fill(HIST("hMass3"), massXiPiPi);
        // covariance matrix elements of PV
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
        registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
        registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
        registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);
        // covariance matrix elements of SV
        registry.fill(HIST("hCovSVXX"), covMatrixXicPlus[0]);
        registry.fill(HIST("hCovSVYY"), covMatrixXicPlus[2]);
        registry.fill(HIST("hCovSVXZ"), covMatrixXicPlus[3]);
        registry.fill(HIST("hCovSVZZ"), covMatrixXicPlus[5]);
        // DCAs of prongs
        registry.fill(HIST("hDcaXYProngs"), kfXi.GetPt(), impactParameterXiXY);
        registry.fill(HIST("hDcaXYProngs"), kfCharmBachelor0.GetPt(), impactParameterPi0XY);
        registry.fill(HIST("hDcaXYProngs"), kfCharmBachelor1.GetPt(), impactParameterPi1XY);
      }

      //------------------------------fill candidate table rows--------------------------------------
      rowCandidateBase(collision.globalIndex(),
                       kfPv.GetX(), kfPv.GetY(), kfPv.GetZ(),
                       std::sqrt(covMatrixPV[0]), std::sqrt(covMatrixPV[2]), std::sqrt(covMatrixPV[5]),
                       /*3-prong specific columns*/
                       rowTrackIndexXicPlus.cascadeId(), rowTrackIndexXicPlus.prong0Id(), rowTrackIndexXicPlus.prong1Id(),
                       casc.bachelorId(), casc.posTrackId(), casc.negTrackId(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       kfXicPlus.GetErrX(), kfXicPlus.GetErrY(), kfXicPlus.GetErrZ(),
                       kfXicPlus.GetErrDecayLength(), kfXicPlus.GetErrDecayLengthXY(),
                       chi2GeoXicPlus, massXiPiPi, signXic,
                       kfXi.GetPx(), kfXi.GetPy(), kfXi.GetPz(),
                       kfCharmBachelor0.GetPx(), kfCharmBachelor0.GetPy(), kfCharmBachelor0.GetPz(),
                       kfCharmBachelor1.GetPx(), kfCharmBachelor1.GetPy(), kfCharmBachelor1.GetPz(),
                       impactParameterXiXY, impactParameterPi0XY, impactParameterPi1XY,
                       errImpactParameterXiXY, errImpactParameterPi0XY, errImpactParameterPi1XY,
                       /*cascade specific columns*/
                       trackPionFromXi.p(), pPiFromLambda, pPrFromLambda,
                       cpaXi, cpaXYXi, cpaLambda, cpaXYLambda, cpaLambdaToXi, cpaXYLambdaToXi,
                       casc.mXi(), casc.mLambda(), massXiPi0, massXiPi1,
                       /*DCA information*/
                       casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(),
                       casc.dcaXYCascToPV(), casc.dcaZCascToPV(),
                       /*PID information*/
                       nSigTpcPiFromXicPlus0, nSigTpcPiFromXicPlus1, nSigTpcBachelorPi, nSigTpcPiFromLambda, nSigTpcPrFromLambda,
                       nSigTofPiFromXicPlus0, nSigTofPiFromXicPlus1, nSigTofBachelorPi, nSigTofPiFromLambda, nSigTofPrFromLambda);
      rowCandidateKF(casc.kfCascadeChi2(), casc.kfV0Chi2(),
                     chi2topoXicPlusToPVBeforeConstraint, chi2topoXicPlusToPV, chi2topoXiToXicPlusBeforeConstraint, chi2topoXiToXicPlus,
                     dcaXYPi0Pi1, dcaXYPi0Xi, dcaXYPi1Xi,
                     dcaPi0Pi1, dcaPi0Xi, dcaPi1Xi);
    } // loop over track triplets
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processXicplusWithKFParticle, "Run candidate creator with KFParticle using derived data from HfTrackIndexSkimCreatorLfCascades.", false);
}; // struct

/// Performs MC matching.
struct HfCandidateCreatorXicToXiPiPiExpressions {
  Spawns<aod::HfCandXicExt> rowCandidateXic;
  Produces<aod::HfCandXicMcRec> rowMcMatchRec;
  Produces<aod::HfCandXicMcGen> rowMcMatchGen;

  void init(InitContext const&) {}

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles)
  {
    rowCandidateXic->bindExternalIndices(&tracks);

    int indexRec = -1;
    int indexRecXicPlus = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;
    // for resonance matching:
    std::vector<int> arrDaughIndex;
    std::array<int, 2> arrPDGDaugh;
    std::array<int, 2> arrXiResonance = {3324, kPiPlus}; // 3324: Ξ(1530)

    // Match reconstructed candidates.
    for (const auto& candidate : *rowCandidateXic) {
      flag = 0;
      sign = 0;
      origin = RecoDecay::OriginType::None;
      debug = 0;
      arrDaughIndex.clear();

      auto arrayDaughters = std::array{candidate.pi0_as<aod::TracksWMc>(),       // pi <- Xic
                                       candidate.pi1_as<aod::TracksWMc>(),       // pi <- Xic
                                       candidate.bachelor_as<aod::TracksWMc>(),  // pi <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),  // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Xic → pi pi pi pi p
      indexRec = RecoDecay::getMatchedMCRec<false, true, false, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 4);
      indexRecXicPlus = indexRec;
      if (indexRec == -1) {
        debug = 1;
      }
      if (indexRec > -1) {
        // Xi- → pi pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, &sign, 2);
        if (indexRec == -1) {
          debug = 2;
        }
        if (indexRec > -1) {
          // Lambda → p pi
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &sign, 1);
          if (indexRec == -1) {
            debug = 3;
          }
          if (indexRec > -1) {
            RecoDecay::getDaughters(mcParticles.rawIteratorAt(indexRecXicPlus), &arrDaughIndex, std::array{0}, 1);
            if (arrDaughIndex.size() == 2) {
              for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
                arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
              }
              if ((arrPDGDaugh[0] == arrXiResonance[0] && arrPDGDaugh[1] == arrXiResonance[1]) || (arrPDGDaugh[0] == arrXiResonance[1] && arrPDGDaugh[1] == arrXiResonance[0])) {
                flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi);
              } else {
                debug = 4;
              }
            } else {
              flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi);
            }
          }
        }
      }

      // Check whether the charm baryon is non-prompt (from a b quark).
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRecXicPlus);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false);
      }

      rowMcMatchRec(flag, debug, origin);
    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      flag = 0;
      sign = 0;
      debug = 0;
      origin = RecoDecay::OriginType::None;
      arrDaughIndex.clear();

      //  Xic → Xi pi pi
      if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, Pdg::kXiCPlus, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, true, &sign, 2)) {
        debug = 1;
        // Xi- -> Lambda pi
        auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
        // Find Xi- from Xi(1530) -> Xi pi in case of resonant decay
        RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
        if (arrDaughIndex.size() == 2) {
          auto cascStarMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, cascStarMC, +3324, std::array{+kXiMinus, +kPiPlus}, true)) {
            cascMC = mcParticles.rawIteratorAt(cascStarMC.daughtersIds().front());
          }
        }
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, cascMC, +kXiMinus, std::array{+kLambda0, +kPiMinus}, true)) {
          debug = 2;
          // Lambda -> p pi
          auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, v0MC, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
            debug = 3;
            if (arrDaughIndex.size() == 2) {
              for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
                arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
              }
              if ((arrPDGDaugh[0] == arrXiResonance[0] && arrPDGDaugh[1] == arrXiResonance[1]) || (arrPDGDaugh[0] == arrXiResonance[1] && arrPDGDaugh[1] == arrXiResonance[0])) {
                flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi);
              } else {
                debug = 4;
              }
            } else {
              flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi);
            }
          }
        }
      }

      // Check whether the charm baryon is non-prompt (from a b quark).
      if (flag != 0) {
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false);
      }

      rowMcMatchGen(flag, debug, origin);
    } // close loop over generated particles
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPiExpressions, processMc, "Process MC", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXicToXiPiPi>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXicToXiPiPiExpressions>(cfgc)};
}
