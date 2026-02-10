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

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/KFparticle/KFUtilities.h"

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
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::hf_evsel;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::aod::hf_cand_xic_to_xi_pi_pi;

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
  Configurable<bool> useXiMassConstraint{"useXiMassConstraint", true, "Use mass constraint for Xi"};
  Configurable<bool> constrainXicPlusToPv{"constrainXicPlusToPv", false, "Constrain XicPlus to PV"};
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "Construct method of XicPlus: 0 fast mathematics without constraint of fixed daughter particle masses, 2 daughter particle masses stay fixed in construction process"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions (effective only for KFParticle w/o derived data)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  o2::vertexing::DCAFitterN<3> df;

  HfEventSelection hfEvSel;

  int runNumber{0};
  float massXiPiPi{0.};
  float massXiPi0{0.};
  float massXiPi1{0.};
  double bz{0.};
  enum XicCandCounter { TotalSkimmedTriplets = 0,
                        SelEvent,
                        CascPreSel,
                        VertexFit };

  using CascadesLinked = soa::Join<aod::Cascades, aod::CascDataLink>;
  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  using KFCascadesLinked = soa::Join<aod::Cascades, aod::KFCascDataLink>;
  using KFCascFull = soa::Join<aod::KFCascDatas, aod::KFCascCovs>;
  using TracksWCovDcaPidPrPi = soa::Join<aod::TracksWCovDca, aod::TracksPidPr, aod::TracksPidPi>;
  using TracksWCovExtraPidPrPi = soa::Join<aod::TracksWCovExtra, aod::TracksPidPr, aod::TracksPidPi>;

  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext const&)
  {
    if ((doprocessNoCentXicplusWithDcaFitter + doprocessCentFT0CXicplusWithDcaFitter + doprocessCentFT0MXicplusWithDcaFitter + doprocessNoCentXicplusWithKFParticle + doprocessCentFT0CXicplusWithKFParticle + doprocessCentFT0MXicplusWithKFParticle) != 1) {
      LOGP(fatal, "Only one process function for the Xic reconstruction can be enabled at a time.");
    }

    // add histograms to registry
    registry.add("hVertexerType", "Use DCAFitter or KFParticle;;entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.get<TH1>(HIST("hVertexerType"))->GetXaxis()->SetBinLabel(1 + aod::hf_cand::VertexerType::DCAFitter, "DCAFitter");
    registry.get<TH1>(HIST("hVertexerType"))->GetXaxis()->SetBinLabel(1 + aod::hf_cand::VertexerType::KfParticle, "KFParticle");
    registry.add("hCandCounter", "hCandCounter", {HistType::kTH1D, {{4, -0.5, 3.5}}});
    registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + TotalSkimmedTriplets, "total");
    registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + SelEvent, "Event selected");
    registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + CascPreSel, "Cascade preselection");
    registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + VertexFit, "Successful vertex fit");
    // physical variables
    if (fillHistograms) {
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
    if ((doprocessNoCentXicplusWithDcaFitter || doprocessCentFT0CXicplusWithDcaFitter || doprocessCentFT0MXicplusWithDcaFitter) && fillHistograms) {
      registry.fill(HIST("hVertexerType"), aod::hf_cand::VertexerType::DCAFitter);
    }
    if ((doprocessNoCentXicplusWithKFParticle || doprocessCentFT0CXicplusWithKFParticle || doprocessCentFT0MXicplusWithKFParticle) && fillHistograms) {
      registry.fill(HIST("hVertexerType"), aod::hf_cand::VertexerType::KfParticle);
    }

    // initialize CCDB
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    // initialize HF event selection helper
    hfEvSel.init(registry, &zorroSummary);

    // initialize 3-prong vertex fitter
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename Collision>
  void runXicplusCreatorWithDcaFitter(Collision const&,
                                      aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                      CascadesLinked const&,
                                      CascFull const&,
                                      TracksWCovDcaPidPrPi const&,
                                      aod::BCsWithTimestamps const&)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      registry.fill(HIST("hCandCounter"), TotalSkimmedTriplets);

      // check if the event is selected
      auto collision = rowTrackIndexXicPlus.collision_as<Collision>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }
      registry.fill(HIST("hCandCounter"), SelEvent);

      // Retrieve skimmed cascade and pion tracks
      auto cascAodElement = rowTrackIndexXicPlus.cascade_as<CascadesLinked>();
      if (!cascAodElement.has_cascData()) {
        continue;
      }
      auto casc = cascAodElement.cascData_as<CascFull>();
      auto trackCharmBachelor0 = rowTrackIndexXicPlus.prong0_as<TracksWCovDcaPidPrPi>();
      auto trackCharmBachelor1 = rowTrackIndexXicPlus.prong1_as<TracksWCovDcaPidPrPi>();

      // preselect cascade candidates
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - MassXiMinus) > massToleranceCascade) {
          continue;
        }
      }
      registry.fill(HIST("hCandCounter"), CascPreSel);

      //----------------------Set the magnetic field from ccdb---------------------------------------
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

      //--------------------------info of V0 and cascades track from LF-tables---------------------------
      std::array<float, 3> const vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> const pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> const vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> const pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};

      //----------------create cascade track------------------------------------------------------------
      constexpr std::size_t NElementsCovMatrix{6u};
      constexpr std::array<int, NElementsCovMatrix> MomInd = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (auto i = 0u; i < NElementsCovMatrix; i++) {
        covCasc[i] = casc.positionCovMat()[i];
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
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
      registry.fill(HIST("hCandCounter"), VertexFit);

      //----------------------------calculate physical properties-----------------------
      // Charge of charm baryon
      int8_t const signXic = casc.sign() < 0 ? +1 : -1;

      // get SV properties
      const auto& secondaryVertex = df.getPCACandidate();
      auto chi2SV = df.getChi2AtPCACandidate();
      auto covMatrixSV = df.calcPCACovMatrixFlat();

      // get track momenta
      trackCasc = df.getTrack(0);
      trackParCovCharmBachelor0 = df.getTrack(1);
      trackParCovCharmBachelor1 = df.getTrack(2);
      std::array<float, 3> pVecXi{};
      std::array<float, 3> pVecPi0{};
      std::array<float, 3> pVecPi1{};
      trackCasc.getPxPyPzGlo(pVecXi);
      trackParCovCharmBachelor0.getPxPyPzGlo(pVecPi0);
      trackParCovCharmBachelor1.getPxPyPzGlo(pVecPi1);

      // get invariant mass of Xic candidate
      auto arrayMomenta = std::array{pVecXi, pVecPi0, pVecPi1};
      massXiPiPi = RecoDecay::m(arrayMomenta, std::array{MassXiMinus, MassPiPlus, MassPiPlus});

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
      std::array<float, 3> const pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      float cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float const cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float const cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float const cpaLambdaToXi = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float const cpaXYLambdaToXi = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

      // get invariant mass of Xi-pi pairs
      auto arrayMomentaXiPi0 = std::array{pVecXi, pVecPi0};
      massXiPi0 = RecoDecay::m(arrayMomentaXiPi0, std::array{MassXiMinus, MassPiPlus});
      auto arrayMomentaXiPi1 = std::array{pVecXi, pVecPi1};
      massXiPi1 = RecoDecay::m(arrayMomentaXiPi1, std::array{MassXiMinus, MassPiPlus});

      // get uncertainty of the decay length
      float phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      //--------------------- get PID information-----------------------
      float const nSigTpcPiFromXicPlus0 = trackCharmBachelor0.tpcNSigmaPi();
      float const nSigTofPiFromXicPlus0 = trackCharmBachelor0.tofNSigmaPi();
      float const nSigTpcPiFromXicPlus1 = trackCharmBachelor1.tpcNSigmaPi();
      float const nSigTofPiFromXicPlus1 = trackCharmBachelor1.tofNSigmaPi();
      // Bachelor pion
      auto trackPionFromXi = casc.bachelor_as<TracksWCovDcaPidPrPi>();
      float const nSigTpcBachelorPi = trackPionFromXi.tpcNSigmaPi();
      float const nSigTofBachelorPi = trackPionFromXi.tofNSigmaPi();
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

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename Collision>
  void runXicplusCreatorWithKFParticle(Collision const&,
                                       aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                       KFCascadesLinked const&,
                                       KFCascFull const&,
                                       TracksWCovExtraPidPrPi const&,
                                       aod::BCsWithTimestamps const&)
  {
    // loop over triplets of track indices
    for (const auto& rowTrackIndexXicPlus : rowsTrackIndexXicPlus) {
      registry.fill(HIST("hCandCounter"), TotalSkimmedTriplets);

      // check if the event is selected
      auto collision = rowTrackIndexXicPlus.collision_as<Collision>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }
      registry.fill(HIST("hCandCounter"), SelEvent);

      // Retrieve skimmed cascade and pion tracks
      auto cascAodElement = rowTrackIndexXicPlus.cascade_as<aod::KFCascadesLinked>();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<KFCascFull>();
      auto trackCharmBachelor0 = rowTrackIndexXicPlus.prong0_as<TracksWCovExtraPidPrPi>();
      auto trackCharmBachelor1 = rowTrackIndexXicPlus.prong1_as<TracksWCovExtraPidPrPi>();

      //-------------------preselect cascade candidates--------------------------------------
      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - MassXiMinus) > massToleranceCascade) {
          continue;
        }
      }
      registry.fill(HIST("hCandCounter"), CascPreSel);

      //----------------------Set the magnetic field from ccdb-----------------------------
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

      //----------------------info of V0 and cascade tracks from LF-table------------------
      std::array<float, 3> const vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> const pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> const vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> const pVecCasc = {casc.px(), casc.py(), casc.pz()};

      //----------------------Create XicPlus as KFParticle object-------------------------------------------
      // initialize primary vertex
      KFPVertex const kfpVertex = createKFPVertexFromCollision(collision);
      float covMatrixPV[6];
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle const kfPv(kfpVertex); // for calculation of DCAs to PV

      // convert pion tracks into KFParticle object
      KFPTrack const kfpTrackCharmBachelor0 = createKFPTrackFromTrack(trackCharmBachelor0);
      KFPTrack const kfpTrackCharmBachelor1 = createKFPTrackFromTrack(trackCharmBachelor1);
      KFParticle kfCharmBachelor0(kfpTrackCharmBachelor0, kPiPlus);
      KFParticle kfCharmBachelor1(kfpTrackCharmBachelor1, kPiPlus);

      // create Xi as KFParticle object
      // read {X,Y,Z,Px,Py,Pz} and corresponding covariance matrix from KF cascade Tables
      constexpr std::size_t NElementsStateVector{6};
      std::array<float, NElementsStateVector> xyzpxpypz = {casc.x(), casc.y(), casc.z(), casc.px(), casc.py(), casc.pz()};
      float parPosMom[NElementsStateVector];
      std::copy(xyzpxpypz.begin(), xyzpxpypz.end(), parPosMom);
      // create KFParticle
      KFParticle kfXi;
      float const massXi = casc.mXi();
      kfXi.Create(parPosMom, casc.kfTrackCovMat(), casc.sign(), massXi);
      if (useXiMassConstraint) {
        kfXi.SetNonlinearMassConstraint(MassXiMinus);
      }

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
      registry.fill(HIST("hCandCounter"), VertexFit);

      // get chi2 values
      float const chi2GeoXicPlus = kfXicPlus.GetChi2() / kfXicPlus.GetNDF();
      float chi2PrimXi = kfXi.GetDeviationFromVertex(kfPv);
      float chi2PrimPi0 = kfCharmBachelor0.GetDeviationFromVertex(kfPv);
      float chi2PrimPi1 = kfCharmBachelor1.GetDeviationFromVertex(kfPv);

      // topological constraint of Xic to PV
      float chi2TopoXicPlusToPVBefConst = kfXicPlus.GetDeviationFromVertex(kfPv);
      KFParticle kfXicPlusToPV = kfXicPlus;
      kfXicPlusToPV.SetProductionVertex(kfPv);
      float chi2TopoXicPlusToPV = kfXicPlusToPV.GetChi2() / kfXicPlusToPV.GetNDF();
      if (constrainXicPlusToPv) {
        kfXicPlus = kfXicPlusToPV;
        kfXicPlus.TransportToDecayVertex();
      }

      //---------------------calculate physical parameters of XicPlus candidate----------------------
      // sign of charm baryon
      int8_t const signXic = casc.sign() < 0 ? +1 : -1;

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
      std::array<float, 3> const pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      float cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float const cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float const cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float const cpaLambdaToXi = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float const cpaXYLambdaToXi = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

      // get chi2 deviation of Pi0-Pi1, Pi0-Xi, Pi1-Xi
      float chi2DevPi0Pi1 = kfCharmBachelor0.GetDeviationFromParticle(kfCharmBachelor1);
      float chi2DevPi0Xi = kfCharmBachelor0.GetDeviationFromParticle(kfXi);
      float chi2DevPi1Xi = kfCharmBachelor1.GetDeviationFromParticle(kfXi);

      // get DCAs of Pi0-Pi1, Pi0-Xi, Pi1-Xi
      float dcaPi0Pi1 = kfCharmBachelor0.GetDistanceFromParticle(kfCharmBachelor1);
      float dcaPi0Xi = kfCharmBachelor0.GetDistanceFromParticle(kfXi);
      float dcaPi1Xi = kfCharmBachelor1.GetDistanceFromParticle(kfXi);
      float dcaXYPi0Pi1 = kfCharmBachelor0.GetDistanceFromParticleXY(kfCharmBachelor1);
      float dcaXYPi0Xi = kfCharmBachelor0.GetDistanceFromParticleXY(kfXi);
      float dcaXYPi1Xi = kfCharmBachelor1.GetDistanceFromParticleXY(kfXi);

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

      // decay length of XicPlus
      // use XicPlus constrained to PV (kfXicPlusToPV), since production point must be set before calling GetDecayLength(XY) on KFParticle
      float kfDecayLength = 0., errorKfDecayLength = 0., kfDecayLengthXY = 0., errorKfDecayLengthXY = 0.;
      kfXicPlusToPV.GetDecayLength(kfDecayLength, errorKfDecayLength);
      kfXicPlusToPV.GetDecayLengthXY(kfDecayLengthXY, errorKfDecayLengthXY);
      float kfDecayLengthNormalised = ldlFromKF(kfXicPlus, kfPv);
      float kfDecayLengthXYNormalised = ldlXYFromKF(kfXicPlus, kfPv);

      //--------------------- get PID information-----------------------
      float const nSigTpcPiFromXicPlus0 = trackCharmBachelor0.tpcNSigmaPi();
      float const nSigTofPiFromXicPlus0 = trackCharmBachelor0.tofNSigmaPi();
      float const nSigTpcPiFromXicPlus1 = trackCharmBachelor1.tpcNSigmaPi();
      float const nSigTofPiFromXicPlus1 = trackCharmBachelor1.tofNSigmaPi();
      // Bachelor pion
      auto trackPionFromXi = casc.bachelor_as<TracksWCovExtraPidPrPi>();
      float const nSigTpcBachelorPi = trackPionFromXi.tpcNSigmaPi();
      float const nSigTofBachelorPi = trackPionFromXi.tofNSigmaPi();
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
        auto* covMatrixXicPlus = kfXicPlus.CovarianceMatrix();
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
                       errorKfDecayLength, errorKfDecayLengthXY,
                       chi2GeoXicPlus, massXiPiPi, signXic,
                       kfXi.GetPx(), kfXi.GetPy(), kfXi.GetPz(),
                       kfCharmBachelor0.GetPx(), kfCharmBachelor0.GetPy(), kfCharmBachelor0.GetPz(),
                       kfCharmBachelor1.GetPx(), kfCharmBachelor1.GetPy(), kfCharmBachelor1.GetPz(),
                       impactParameterXiXY, impactParameterPi0XY, impactParameterPi1XY,
                       errImpactParameterXiXY, errImpactParameterPi0XY, errImpactParameterPi1XY,
                       /*cascade specific columns*/
                       trackPionFromXi.p(), pPiFromLambda, pPrFromLambda,
                       cpaXi, cpaXYXi, cpaLambda, cpaXYLambda, cpaLambdaToXi, cpaXYLambdaToXi,
                       massXi, casc.mLambda(), massXiPi0, massXiPi1,
                       /*DCA information*/
                       casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(),
                       casc.dcaXYCascToPV(), casc.dcaZCascToPV(),
                       /*PID information*/
                       nSigTpcPiFromXicPlus0, nSigTpcPiFromXicPlus1, nSigTpcBachelorPi, nSigTpcPiFromLambda, nSigTpcPrFromLambda,
                       nSigTofPiFromXicPlus0, nSigTofPiFromXicPlus1, nSigTofBachelorPi, nSigTofPiFromLambda, nSigTofPrFromLambda);
      rowCandidateKF(kfDecayLength, kfDecayLengthNormalised, kfDecayLengthXY, kfDecayLengthXYNormalised,
                     casc.kfCascadeChi2(), casc.kfV0Chi2(),
                     chi2TopoXicPlusToPVBefConst, chi2TopoXicPlusToPV,
                     chi2PrimXi, chi2PrimPi0, chi2PrimPi1,
                     chi2DevPi0Pi1, chi2DevPi0Xi, chi2DevPi1Xi,
                     dcaPi0Pi1, dcaPi0Xi, dcaPi1Xi,
                     dcaXYPi0Pi1, dcaXYPi0Xi, dcaXYPi1Xi);
    } // loop over track triplets
  }

  ///////////////////////////////////////////////////////////
  ///                                                     ///
  ///           Process functions with DCAFitter          ///
  ///                                                     ///
  ///////////////////////////////////////////////////////////

  void processNoCentXicplusWithDcaFitter(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                         aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                         CascadesLinked const& cascadesLinked,
                                         CascFull const& cascadesFull,
                                         TracksWCovDcaPidPrPi const& tracks,
                                         aod::BCsWithTimestamps const& bcs)
  {
    runXicplusCreatorWithDcaFitter<o2::hf_centrality::CentralityEstimator::None>(collisions, rowsTrackIndexXicPlus, cascadesLinked, cascadesFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processNoCentXicplusWithDcaFitter, "Run candidate creator with DCAFitter without centrality selection.", true);

  void processCentFT0CXicplusWithDcaFitter(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                           aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                           CascadesLinked const& cascadesLinked,
                                           CascFull const& cascadesFull,
                                           TracksWCovDcaPidPrPi const& tracks,
                                           aod::BCsWithTimestamps const& bcs)
  {
    runXicplusCreatorWithDcaFitter<o2::hf_centrality::CentralityEstimator::FT0C>(collisions, rowsTrackIndexXicPlus, cascadesLinked, cascadesFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCentFT0CXicplusWithDcaFitter, "Run candidate creator with DCAFitter with centrality selection on FT0C.", false);

  void processCentFT0MXicplusWithDcaFitter(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                           aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                           CascadesLinked const& cascadesLinked,
                                           CascFull const& cascadesFull,
                                           TracksWCovDcaPidPrPi const& tracks,
                                           aod::BCsWithTimestamps const& bcs)
  {
    runXicplusCreatorWithDcaFitter<o2::hf_centrality::CentralityEstimator::FT0M>(collisions, rowsTrackIndexXicPlus, cascadesLinked, cascadesFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCentFT0MXicplusWithDcaFitter, "Run candidate creator with DCAFitter with centrality selection on FT0M.", false);

  ///////////////////////////////////////////////////////////
  ///                                                     ///
  ///        Process functions with KFParticle            ///
  ///                                                     ///
  ///////////////////////////////////////////////////////////

  void processNoCentXicplusWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                          aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                          KFCascadesLinked const& kfCascadesLinked,
                                          KFCascFull const& kfCascadesFull,
                                          TracksWCovExtraPidPrPi const& tracks,
                                          aod::BCsWithTimestamps const& bcs)
  {
    runXicplusCreatorWithKFParticle<o2::hf_centrality::CentralityEstimator::None>(collisions, rowsTrackIndexXicPlus, kfCascadesLinked, kfCascadesFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processNoCentXicplusWithKFParticle, "Run candidate creator with KFParticle without centrality selection.", false);

  void processCentFT0CXicplusWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                            aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                            KFCascadesLinked const& kfCascadesLinked,
                                            KFCascFull const& kfCascadesFull,
                                            TracksWCovExtraPidPrPi const& tracks,
                                            aod::BCsWithTimestamps const& bcs)
  {
    runXicplusCreatorWithKFParticle<o2::hf_centrality::CentralityEstimator::FT0C>(collisions, rowsTrackIndexXicPlus, kfCascadesLinked, kfCascadesFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCentFT0CXicplusWithKFParticle, "Run candidate creator with KFParticle with centrality selection on FT0C.", false);

  void processCentFT0MXicplusWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                            aod::HfCascLf3Prongs const& rowsTrackIndexXicPlus,
                                            KFCascadesLinked const& kfCascadesLinked,
                                            KFCascFull const& kfCascadesFull,
                                            TracksWCovExtraPidPrPi const& tracks,
                                            aod::BCsWithTimestamps const& bcs)
  {
    runXicplusCreatorWithKFParticle<o2::hf_centrality::CentralityEstimator::FT0M>(collisions, rowsTrackIndexXicPlus, kfCascadesLinked, kfCascadesFull, tracks, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCentFT0MXicplusWithKFParticle, "Run candidate creator with KFParticle with centrality selection on FT0M.", false);

  ///////////////////////////////////////////////////////////
  ///                                                     ///
  ///   Process functions only for collision monitoring   ///
  ///                                                     ///
  ///////////////////////////////////////////////////////////

  void processCollisions(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCsWithTimestamps const&)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCollisions, "Collision monitoring - no centrality", false);

  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::BCsWithTimestamps const&)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions, aod::BCsWithTimestamps const&)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPi, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);
}; // struct

/// Performs MC matching.
struct HfCandidateCreatorXicToXiPiPiExpressions {
  Spawns<aod::HfCandXicExt> rowCandidateXic;
  Produces<aod::HfCandXicMcRec> rowMcMatchRec;
  Produces<aod::HfCandXicMcGen> rowMcMatchGen;
  Produces<aod::HfCandXicResid> rowResiduals;

  Configurable<bool> fillMcHistograms{"fillMcHistograms", true, "Fill validation plots"};
  Configurable<bool> matchDecayedPions{"matchDecayedPions", true, "Match also candidates with daughter pion tracks that decay with kinked topology"};
  Configurable<bool> matchInteractionsWithMaterial{"matchInteractionsWithMaterial", true, "Match also candidates with daughter tracks that interact with material"};
  Configurable<bool> fillResidualTable{"fillResidualTable", false, "Fill table containing residuals and pulls of PV and SV"};

  HfEventSelectionMc hfEvSelMc;

  enum DebugRec { TotalRec = 0,
                  XicToFinalState,
                  XiToPiPPi,
                  LambdaToPPi };

  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    // add histograms to registry
    if (fillMcHistograms) {
      registry.add("hDecayedPions", "hDecayedPions", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("hInteractionsWithMaterial", "hInteractionsWithMaterial", {HistType::kTH1F, {{6, -0.5, 5.5}}});
      registry.add("hDebugRec", "hDebugRec", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      registry.get<TH1>(HIST("hDebugRec"))->GetXaxis()->SetBinLabel(1 + TotalRec, "total");
      registry.get<TH1>(HIST("hDebugRec"))->GetXaxis()->SetBinLabel(1 + XicToFinalState, "#Xi^{+}_{c} #rightarrow #pi^{#plus} #pi^{#plus} #pi^{#minus} p #pi^{#minus}");
      registry.get<TH1>(HIST("hDebugRec"))->GetXaxis()->SetBinLabel(1 + XiToPiPPi, "#Xi^{#minus} #rightarrow #pi^{#minus} p #pi^{#minus}");
      registry.get<TH1>(HIST("hDebugRec"))->GetXaxis()->SetBinLabel(1 + LambdaToPPi, "#Lambda #rightarrow p #pi^{#minus}");
    }

    // initialize HF event selection helper
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-candidate-creator-xic-to-xi-pi-pi") {
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename McCollisions, typename CollInfos>
  void runMcMatching(aod::TracksWMc const& tracks,
                     aod::McParticles const& mcParticles,
                     McCollisions const& mcCollisions,
                     CollInfos const& collInfos,
                     BCsInfo const&)
  {
    rowCandidateXic->bindExternalIndices(&tracks);

    int indexRec = -1;
    int indexRecXicPlus = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = RecoDecay::OriginType::None;
    float decayLengthGen = -999.f;
    int8_t nPionsDecayed = 0;
    int8_t nInteractionsWithMaterial = 0;
    // for resonance matching
    std::vector<int> arrDaughIndex;
    constexpr std::size_t NDaughtersResonant{2u};
    std::array<int, NDaughtersResonant> arrPDGDaugh{};
    std::array<int, NDaughtersResonant> arrXiResonance = {3324, kPiPlus}; // 3324: Ξ(1530)
    // for non-prompt
    std::vector<int> idxBhadMothers;
    // residuals and pulls
    std::array<float, 2> momentumResiduals{-9999.f};
    std::array<float, 3> pvResiduals{-9999.f};
    std::array<float, 3> pvPulls{-9999.f};
    std::array<float, 3> svResiduals{-9999.f};
    std::array<float, 3> svPulls{-9999.f};

    // Match reconstructed candidates.
    for (const auto& candidate : *rowCandidateXic) {
      sign = 0;
      flag = 0;
      origin = RecoDecay::OriginType::None;
      nPionsDecayed = 0;
      nInteractionsWithMaterial = 0;
      arrDaughIndex.clear();
      if (fillResidualTable) {
        momentumResiduals.fill(-9999.f);
        pvResiduals.fill(-9999.f);
        pvPulls.fill(-9999.f);
        svResiduals.fill(-9999.f);
        svPulls.fill(-9999.f);
      }

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

      if (fillMcHistograms) {
        registry.fill(HIST("hDebugRec"), TotalRec);
      }

      // Xic → pi pi pi pi p
      if (matchDecayedPions && matchInteractionsWithMaterial) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 4, &nPionsDecayed, nullptr, &nInteractionsWithMaterial);
      } else if (matchDecayedPions && !matchInteractionsWithMaterial) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 4, &nPionsDecayed, nullptr, &nInteractionsWithMaterial);
      } else if (!matchDecayedPions && matchInteractionsWithMaterial) {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, false, true>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 4, &nPionsDecayed, nullptr, &nInteractionsWithMaterial);
      } else {
        indexRec = RecoDecay::getMatchedMCRec<false, true, false, false, false>(mcParticles, arrayDaughters, Pdg::kXiCPlus, std::array{+kPiPlus, +kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 4, &nPionsDecayed, nullptr, &nInteractionsWithMaterial);
      }
      indexRecXicPlus = indexRec;
      if (indexRec > -1) {
        if (fillMcHistograms) {
          registry.fill(HIST("hDebugRec"), XicToFinalState);
        }
        // Xi- → pi pi p
        if (matchDecayedPions && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, nullptr, 2);
        } else if (matchDecayedPions && !matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, nullptr, 2);
        } else if (!matchDecayedPions && matchInteractionsWithMaterial) {
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, false, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, nullptr, 2);
        } else {
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, false, false>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, nullptr, 2);
        }
        if (indexRec > -1) {
          if (fillMcHistograms) {
            registry.fill(HIST("hDebugRec"), XiToPiPPi);
          }
          // Lambda → p pi
          if (matchDecayedPions && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true);
          } else if (matchDecayedPions && !matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true);
          } else if (!matchDecayedPions && matchInteractionsWithMaterial) {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true);
          } else {
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, false, false>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true);
          }
          if (indexRec > -1) {
            if (fillMcHistograms) {
              registry.fill(HIST("hDebugRec"), LambdaToPPi);
            }
            auto particleXicPlus = mcParticles.rawIteratorAt(indexRecXicPlus);
            // Check whether XicPlus decays via resonant decay
            RecoDecay::getDaughters(particleXicPlus, &arrDaughIndex, std::array{0}, 1);
            if (arrDaughIndex.size() == NDaughtersResonant) {
              for (auto iProng = 0u; iProng < NDaughtersResonant; ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
                arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
              }
              if ((arrPDGDaugh[0] == arrXiResonance[0] && arrPDGDaugh[1] == arrXiResonance[1]) || (arrPDGDaugh[0] == arrXiResonance[1] && arrPDGDaugh[1] == arrXiResonance[0])) {
                flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi);
              }
            } else {
              flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi);
            }
            // Check whether the charm baryon is non-prompt (from a b quark).
            if (flag != 0) {
              origin = RecoDecay::getCharmHadronOrigin(mcParticles, particleXicPlus, false);
            }
            // Calculate residuals and pulls
            if (flag != 0 && fillResidualTable) {
              auto mcCollision = particleXicPlus.template mcCollision_as<McCollisions>();
              auto particleDaughter0 = mcParticles.rawIteratorAt(arrDaughIndex[0]);

              momentumResiduals[0] = candidate.p() - particleXicPlus.p();
              momentumResiduals[1] = candidate.pt() - particleXicPlus.pt();
              pvResiduals[0] = candidate.posX() - mcCollision.posX();
              pvResiduals[1] = candidate.posY() - mcCollision.posY();
              pvResiduals[2] = candidate.posZ() - mcCollision.posZ();
              svResiduals[0] = candidate.xSecondaryVertex() - particleDaughter0.vx();
              svResiduals[1] = candidate.ySecondaryVertex() - particleDaughter0.vy();
              svResiduals[2] = candidate.zSecondaryVertex() - particleDaughter0.vz();
              try {
                pvPulls[0] = pvResiduals[0] / candidate.xPvErr();
                pvPulls[1] = pvResiduals[1] / candidate.yPvErr();
                pvPulls[2] = pvResiduals[2] / candidate.zPvErr();
                svPulls[0] = svResiduals[0] / candidate.xSvErr();
                svPulls[1] = svResiduals[1] / candidate.ySvErr();
                svPulls[2] = svResiduals[2] / candidate.zSvErr();
              } catch (const std::runtime_error& error) {
                LOG(info) << "Run time error found: " << error.what() << ". Set values of vertex pulls to -9999.9.";
              }
            }
          }
        }
      }

      // Fill histograms
      if (flag != 0 && fillMcHistograms) {
        registry.fill(HIST("hDecayedPions"), nPionsDecayed);
        registry.fill(HIST("hInteractionsWithMaterial"), nInteractionsWithMaterial);
      }

      // Fill tables
      rowMcMatchRec(flag, origin);
      if (flag != 0 && fillResidualTable) {
        rowResiduals(origin, momentumResiduals[0], momentumResiduals[1],
                     pvResiduals[0], pvResiduals[1], pvResiduals[2],
                     pvPulls[0], pvPulls[1], pvPulls[2],
                     svResiduals[0], svResiduals[1], svResiduals[2],
                     svPulls[0], svPulls[1], svPulls[2]);
      }
    } // close loop over candidates

    // Match generated particles.
    for (const auto& mcCollision : mcCollisions) {
      // Slice the particles table to get the particles for the current MC collision
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      // Slice the collisions table to get the collision info for the current MC collision
      float centrality{-1.f};
      o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
      int nSplitColl = 0;
      if constexpr (CentEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (CentEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        nSplitColl = collSlice.size();
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (CentEstimator == o2::hf_centrality::CentralityEstimator::None) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
      }
      hfEvSelMc.fillHistograms<CentEstimator>(mcCollision, rejectionMask, nSplitColl);
      if (rejectionMask != 0) {
        // at least one event selection not satisfied --> reject all particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          rowMcMatchGen(-99, -99, -99, decayLengthGen);
        }
        continue;
      }

      for (const auto& particle : mcParticlesPerMcColl) {
        sign = 0;
        flag = 0;
        origin = RecoDecay::OriginType::None;
        arrDaughIndex.clear();
        idxBhadMothers.clear();

        //  Xic → Xi pi pi
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, Pdg::kXiCPlus, std::array{+kXiMinus, +kPiPlus, +kPiPlus}, true, &sign, 2)) {
          // Xi- -> Lambda pi
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          // Find Xi- from Xi(1530) -> Xi pi in case of resonant decay
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == NDaughtersResonant) {
            auto cascStarMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, cascStarMC, +3324, std::array{+kXiMinus, +kPiPlus}, true)) {
              cascMC = mcParticles.rawIteratorAt(cascStarMC.daughtersIds().front());
            }
          }
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, cascMC, +kXiMinus, std::array{+kLambda0, +kPiMinus}, true)) {
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, v0MC, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
              if (arrDaughIndex.size() == NDaughtersResonant) {
                for (auto iProng = 0u; iProng < NDaughtersResonant; ++iProng) {
                  auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
                  arrPDGDaugh[iProng] = std::abs(daughI.pdgCode());
                }
                if ((arrPDGDaugh[0] == arrXiResonance[0] && arrPDGDaugh[1] == arrXiResonance[1]) || (arrPDGDaugh[0] == arrXiResonance[1] && arrPDGDaugh[1] == arrXiResonance[0])) {
                  flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiResPiToXiPiPi);
                }
              } else {
                flag = sign * (1 << aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi);
              }
            }
          }
        }

        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
          // Calculate the decay length of the generated particle
          auto dau0 = particle.template daughters_as<aod::McParticles>().begin();
          const std::array vtxDau{dau0.vx(), dau0.vy(), dau0.vz()};
          const std::array vtxPV{mcCollision.posX(), mcCollision.posY(), mcCollision.posZ()};
          decayLengthGen = RecoDecay::distance(vtxPV, vtxDau);
        }
        // Fill table
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMcMatchGen(flag, origin, bHadMother.pdgCode(), decayLengthGen);
        } else {
          rowMcMatchGen(flag, origin, 0, decayLengthGen);
        }
      } // close loop over generated particles
    } // close loop over McCollisions
  } // close template function

  void processMc(aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 aod::McCollisions const& mcCollisions,
                 McCollisionsNoCents const& mcCollisionsNoCents,
                 BCsInfo const& bcs)
  {
    runMcMatching<o2::hf_centrality::CentralityEstimator::None>(tracks, mcParticles, mcCollisions, mcCollisionsNoCents, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPiExpressions, processMc, "Perform MC matching with no centrality selection.", true);

  void processMcCentFT0C(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         aod::McCollisions const& mcCollisions,
                         McCollisionsFT0Cs const& mcCollisionsFT0Cs,
                         BCsInfo const& bcs)
  {
    runMcMatching<o2::hf_centrality::CentralityEstimator::FT0C>(tracks, mcParticles, mcCollisions, mcCollisionsFT0Cs, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPiExpressions, processMcCentFT0C, "Perform MC matching with centrality selection on FT0C.", false);

  void processMcCentFT0M(aod::TracksWMc const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsCentFT0Ms const& mcCollisionsCentFT0Ms,
                         McCollisionsFT0Ms const& mcCollisionsFT0Ms,
                         BCsInfo const& bcs)
  {
    runMcMatching<o2::hf_centrality::CentralityEstimator::FT0M>(tracks, mcParticles, mcCollisionsCentFT0Ms, mcCollisionsFT0Ms, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXicToXiPiPiExpressions, processMcCentFT0M, "Perform MC matching with centrality selection on FT0M.", false);
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXicToXiPiPi>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXicToXiPiPiExpressions>(cfgc)};
}
