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

/// \file candidateCreatorXic0XicpToHadronic.cxx
/// \brief Reconstruction of Xic0 and Xicp candiates with hadronic decay chain
///
/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>
#include <KFParticleBase.h>
#include <KFVertex.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::constants::physics;
using namespace o2::hf_evsel;

struct HfCandidateCreatorXic0XicpToHadronic {

  // Cursor to fill tables
  struct : ProducesGroup {
    Produces<aod::HfCandXic0Base> rowCandXic0Base;
    Produces<aod::HfCandXic0KF> rowCandXic0KF;
  } cursors;

  // Configurables
  struct : ConfigurableGroup {
    // Switch for filling histograms
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill validation plots"};
    // Magnetic field setting from CCDB
    Configurable<bool> isRun2{"isRun2", false, "enable Run2 or Run3 GRP objects for magnetic field"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "https://alice-ccdb.cern.ch", "url of the ccdb object"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parameterization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "path of the group file (Run2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run3)"};
    // Cascade preselection - Xic0
    Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
    Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascades"};
    Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in XY plane"};
    // DCAFitter - Xic0
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "Create tracks version propagated to PCA"};
    Configurable<double> maxR{"maxR", 200., "Reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "Reject (if>0) PCA candidate if tracks DZ exceeds this threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "Stop iteration if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "Stop iteration if Chi2/Chi2old > this"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariance, effective only if useAbsDCA is true"};
    // KFParticle - Xic0
    Configurable<bool> useXiMassConstraint{"useXiMassConstraint", true, "use mass constraint for Xi(cascade)"};
    Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks comming from different collisions(effective only for KFParticle w/o derived data)"};
    Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "2 : Daughter particle masses stay fixed in construction process"};
    // KFParticle- Xic0 specific
    Configurable<bool> constrainXic0ToPv{"constrainXic0ToPv", true, "constrain xic0 to pv"};
    Configurable<bool> constrainXiToXic0{"constrainXiToXic0", true, "constrain xi to xic0"};
    // configurbles for histogram binning
    Configurable<int> nBinXiMass{"nBinXiMass", 1000, "nBinXiMass"};
    Configurable<float> xiMassMin{"xiMassMin", 1.0, "xiMassMin"};
    Configurable<float> xiMassMax{"xiMassMax", 2.0, "xiMassMax"};
    Configurable<int> nBinXic0Mass{"nBinXic0Mass", 3000, "nBinXic0Mass"};
    Configurable<float> xic0MassMin{"xic0MassMin", 1.0, "xic0MassMin"};
    Configurable<float> xic0MassMax{"xic0MassMax", 4.0, "xic0MassMax"};
    Configurable<int> nBinCpa2Prong{"nBinCpa2Prong", 240, "nBinCpa2Prong"};
    Configurable<float> cpa2ProngMin{"cpa2ProngMin", -1.2, "cpa2ProngMin"};
    Configurable<float> cpa2ProngMax{"cpa2ProngMax", 1.2, "cpa2ProngMax"};
    Configurable<int> nBinImpParXYXi{"nBinImpParXYXi", 30, "nBinImpParXYXi"};
    Configurable<float> impParXYXiMin{"impParXYXiMin", -1.5, "impParXYXiMin"};
    Configurable<float> impParXYXiMax{"impParXYXiMax", 1.5, "impParXYXiMax"};
    Configurable<int> nBinImpParXYPi{"nBinImpParXYPi", 30, "nBinImpParXYPi"};
    Configurable<float> impParXYPiMin{"impParXYPiMin", -1.5, "impParXYPiMin"};
    Configurable<float> impParXYPiMax{"impParXYPiMax", 1.5, "impParXYPiMax"};
    Configurable<int> nBinPtXi{"nBinPtXi", 100, "nBinPtXi"};
    Configurable<float> ptXiMin{"ptXiMin", 0.0, "ptXiMin"};
    Configurable<float> ptXiMax{"ptXiMax", 20.0, "ptXiMax"};
    Configurable<int> nBinPtPi{"nBinPtPi", 100, "nBinPtPi"};
    Configurable<float> ptPiMin{"ptPiMin", 0.0, "ptPiMin"};
    Configurable<float> ptPiMax{"ptPiMax", 20.0, "ptPiMax"};

  } configs;

  // For magnetic field
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // DCAFitter
  o2::vertexing::DCAFitterN<2> dfXic0;

  int runNumber{0};
  double bz{0.};
  // Mass of Xic0
  float massXiPi{0.};

  enum Xic0CandCounter { All = 0,
                         CascPreSel,
                         VertexFit };

  using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  // For DCAFitter
  using CascadesLinked = soa::Join<aod::Cascades, aod::CascDataLink>;
  using CascFull = soa::Join<aod::CascDatas, aod::CascCovs>;
  // For KFParticle
  using KFCascadesLinked = soa::Join<aod::Cascades, aod::KFCascDataLink>;
  using KFCascFull = soa::Join<aod::KFCascDatas, aod::KFCascCovs>;
  using TracksWCovDcaExtraPidPrPi = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPr, aod::TracksPidPi>;

  HistogramRegistry registry{"hists"};
  HfEventSelection hfEvSel;

  void init(InitContext const&)
  {
    std::vector<bool> processesXic0{doprocessXic0WithDCAFitterNoCent, doprocessXic0WithDCAFitterCentFT0C, doprocessXic0WithDCAFitterCentFT0M,
                                    doprocessXic0WithKFParticleNoCent, doprocessXic0WithKFParticleCentFT0C, doprocessXic0WithKFParticleCentFT0M};
    std::vector<bool> processesCollMonitoring{doprocessCollisionsNoCent, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};

    // Exit if workflow is not configured correctly - More than one process enabled for Xic0
    if (std::accumulate(processesXic0.begin(), processesXic0.end(), 0) > 1) {
      LOGP(fatal, "More than one process fucntion for Xic0 was enabled. Please choose only one process function");
    }

    // Exit if workflow is not configured correctly - More than one process enabled for collision monitoring
    if (std::accumulate(processesCollMonitoring.begin(), processesCollMonitoring.end(), 0) > 1) {
      LOGP(fatal, "More than one process fucntion for CollMonitoring was enabled. Please choose only one process function");
    }

    // Initialization for Xic0
    if (std::accumulate(processesXic0.begin(), processesXic0.end(), 0) == 1) {

      // Set histograms
      if (configs.fillHistograms) {
        // counter
        registry.add("hVertexerType", "Use KF of DCAFitterN;Vertexer type;entries", {kTH1F, {{2, -0.5, 1.5}}});
        registry.add("hCandCounter", "hCandCounter", {kTH1F, {{3, -0.5, 2.5}}});
        registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + All, "Total");
        registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + CascPreSel, "Cascade preselection");
        registry.get<TH1>(HIST("hCandCounter"))->GetXaxis()->SetBinLabel(1 + VertexFit, "Successful vertex fit");

        // QA histograms for Xic0 candidates
        registry.add("hXiMassAfterConstrain", "Xi mass after selection;m(GeV);entries", {kTH1F, {{configs.nBinXiMass, configs.xiMassMin, configs.xiMassMax}}});
        registry.add("hXic0Mass", "Xic0 mass after selection;m(GeV);entries", {kTH1F, {{configs.nBinXic0Mass, configs.xic0MassMin, configs.xic0MassMax}}});
        registry.add("hXic0MassAfterPVConstraint", "Xic0 mass after selection;m(GeV);entries", {kTH1F, {{configs.nBinXic0Mass, configs.xic0MassMin, configs.xic0MassMax}}});
        registry.add("hXic0MassAfterXiConstraint", "Xic0 mass after selection;m(GeV);entries", {kTH1F, {{configs.nBinXic0Mass, configs.xic0MassMin, configs.xic0MassMax}}});

        registry.add("hMassXic0Cand", "2-prong candidates;inv.mass(#Xi #pi)(GeV/#it{c}^{2});entries", {kTH1D, {{500, 2.3, 2.7}}});
        registry.add("hCovPVXX", "2-prong candidates;XX element of cov.matrix of prim.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 1.e-4}}});
        registry.add("hCovSVXX", "2-prong candidates;XX element of cov.matrix of sec.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 0.2}}});
        registry.add("hCovPVYY", "2-prong candidates;YY element of cov.matrix of prim.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 1.e-4}}});
        registry.add("hCovSVYY", "2-prong candidates;YY element of cov.matrix of sec.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 0.2}}});
        registry.add("hCovPVZZ", "2-prong candidates;ZZ element of cov.matrix of prim.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 1.e-4}}});
        registry.add("hCovSVZZ", "2-prong candidates;ZZ element of cov.matrix of sec.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 0.2}}});
        registry.add("hCovPVXZ", "2-prong candidates;XZ element of cov.matrix of prim.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 1.e-4}}});
        registry.add("hCovSVXZ", "2-prong candidates;XZ element of cov.matrix of sec.vts. position(cm^{2});entries", {kTH1D, {{100, 0., 0.2}}});
        registry.add("hDcaXYProngs", "DCAxy of 2-prong candidates;#it{p}_{T} (GeV/#it{c});#it{d}_{xy} (#mum);entries", {kTH2D, {{100, 0., 20.}, {200, -500., 500.}}});
        registry.add("hDcaZProngs", "DCAz of 2-prong candidates;#it{p}_{T} (GeV/#it{c});#it{d}_{xy} (#mum);entries", {kTH2D, {{100, 0., 20.}, {200, -500., 500.}}});

        registry.add("hImpParXiXY", "ImpactParameter of Xi;ImpParXi;entries", {kTH1F, {{configs.nBinImpParXYXi, configs.impParXYXiMin, configs.impParXYXiMax}}});
        registry.add("hImpParPiXY", "ImpactParameter of Pi;ImpParPi;entries", {kTH1F, {{configs.nBinImpParXYPi, configs.impParXYPiMin, configs.impParXYPiMax}}});
        registry.add("hPtXi", "Pt of candidate's Xi;#it{p}_{T};entries", {kTH1F, {{configs.nBinPtPi, configs.ptXiMin, configs.ptXiMax}}});
        registry.add("hPtPi", "Pt of candidate's Pi;#it{p}_{T};entries", {kTH1F, {{configs.nBinPtPi, configs.ptPiMin, configs.ptPiMax}}});

        if (doprocessXic0WithDCAFitterNoCent) {
          registry.get<TH1>(HIST("hVertexerType"))->Fill(aod::hf_cand::VertexerType::DCAFitter);
        }

        if (doprocessXic0WithKFParticleNoCent) {
          registry.get<TH1>(HIST("hVertexerType"))->Fill(aod::hf_cand::VertexerType::KfParticle);
        }
      }

      // initialize 2-prong vertex fitter
      dfXic0.setPropagateToPCA(configs.propagateToPCA);
      dfXic0.setMaxR(configs.maxR);
      dfXic0.setMaxDZIni(configs.maxDZIni);
      dfXic0.setMinParamChange(configs.minParamChange);
      dfXic0.setUseAbsDCA(configs.useAbsDCA);
      dfXic0.setWeightedFinalPCA(configs.useWeightedFinalPCA);

    } // end of Xic0 initilization

    // initialize ccdb
    ccdb->setURL(configs.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(configs.ccdbPathLut));
    runNumber = 0;

    // initailize HF event selection helper
    hfEvSel.init(registry);

  } // end of initialization

  ///////////////////////////////////////////////////////////
  //                                                       //
  //                                                       //
  //                      Xic0                             //
  //                                                       //
  ///////////////////////////////////////////////////////////

  // template function for running xic0 reconstruction via DCAFitter method
  // templated for various centrality estimator usage
  template <o2::hf_centrality::CentralityEstimator centEstimator, typename Colls>
  void runCreatorXic0WithDCAFitter(Colls const&,
                                   aod::HfCascLf2Prongs const& candidates,
                                   CascadesLinked const&,
                                   CascFull const&,
                                   TracksWCovDcaExtraPidPrPi const&,
                                   aod::BCsWithTimestamps const&)
  {
    // Loop over candidate
    for (auto const& cand : candidates) {

      // Fill cascandidates before selection
      if (configs.fillHistograms) {
        registry.get<TH1>(HIST("hCandCounter"))->Fill(All);
      }

      // Apply hfflag selection
      if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
        continue;
      }

      // Event selection
      auto collision = cand.collision_as<Colls>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        continue;
      }

      //------------------------------Cascade pre-selection------------------------------

      // Retrieving skimmed cascade and pion tracks
      // If there is no related tracks, skip
      auto cascAodElement = cand.cascade_as<CascadesLinked>();
      if (!cascAodElement.has_cascData()) {
        continue;
      }
      auto casc = cascAodElement.cascData_as<CascFull>();
      auto trackCharmBachelor = cand.prong0_as<TracksWCovDcaExtraPidPrPi>();
      auto chargeCasc = casc.sign() > 0 ? 1 : -1;

      if (configs.doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > configs.dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - MassXiMinus) > configs.massToleranceCascade) {
          continue;
        }
      }

      if (configs.fillHistograms) {
        registry.fill(HIST("hCandCounter"), CascPreSel);
        registry.fill(HIST("hXiMassAfterConstrain"), casc.mXi());
      }

      //------------------------------Set Magnetic field------------------------------
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>> Current run Number : " << runNumber;
        initCCDB(bc, runNumber, ccdb, configs.isRun2 ? configs.ccdbPathGrp : configs.ccdbPathGrpMag, lut, configs.isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>> Magnetic field: " << bz;
      }
      dfXic0.setBz(bz);

      //------------------------------Info of V0 and cascade tracks from LF table------------------------------
      // -> This quantities are used for physical properties of selected candidates
      // -> Not used for candidate creation
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};

      //------------------------------Create cascade track------------------------------
      constexpr std::size_t NElementsCovMatrix{6u};
      constexpr std::array<int, NElementsCovMatrix> MomInd = {9, 13, 14, 18, 19, 20};
      for (auto i = 0u; i < NElementsCovMatrix; i++) {
        covCasc[i] = casc.positionCovMat()[i];
        covCasc[MomInd[i]] = casc.momentumCovMat()[i];
      }

      o2::track::TrackParCov trackCasc;
      if (chargeCasc < 0) { // Xi-
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else if (chargeCasc > 0) { // Xi+
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else {
        continue;
      }

      trackCasc.setAbsCharge(1);
      trackCasc.setPID(o2::track::PID::XiMinus);

      //------------------------------Fit SV & Create Xic0 track------------------------------
      auto trackParCovCharmBachelor = getTrackParCov(trackCharmBachelor);
      try {
        if (dfXic0.process(trackCasc, trackParCovCharmBachelor) == 0) {
          continue;
        }
      } catch (const std::runtime_error& e) {
        LOG(info) << "Run time error found : " << e.what() << ".DCAFitter cannot work with this candidate. SKIP!";
        continue;
      }

      if (configs.fillHistograms) {
        registry.fill(HIST("hCandCounter"), VertexFit); // -> This candidate has successful SV fit
      }

      //------------------------------Calculate physical properties-----------------------------

      // get SV Properties
      auto const& secondaryVertex = dfXic0.getPCACandidate();
      auto chi2SV = dfXic0.getChi2AtPCACandidate();
      auto covMatrixSV = dfXic0.calcPCACovMatrixFlat();

      // get track momenta
      trackCasc = dfXic0.getTrack(0);
      trackParCovCharmBachelor = dfXic0.getTrack(1);
      std::array<float, 3> pVecXi, pVecPi;
      trackCasc.getPxPyPzGlo(pVecXi);
      trackParCovCharmBachelor.getPxPyPzGlo(pVecPi);

      // get invariant mass of Xic0 candidate
      auto arrMomenta = std::array{pVecXi, pVecPi};
      massXiPi = RecoDecay::m(std::move(arrMomenta), std::array{MassXiMinus, MassPiPlus});

      // get impact parameter
      //! This process modifies track momenta
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      // calculate impact parameter
      o2::dataformats::DCA impactParameterCasc, impactParameterPi;
      trackCasc.propagateToDCA(primaryVertex, bz, &impactParameterCasc);
      trackParCovCharmBachelor.propagateToDCA(primaryVertex, bz, &impactParameterPi);

      // calculate cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      // lambda <- V0
      float cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      // Xi <- Xic0
      float cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float cpaLambdaToXi = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float cpaXYLambdaToXi = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

      // get uncertainty of the decay length
      float phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      //------------------------------Get PID Information-----------------------------
      // Charm Bachelor Pion
      float nSigTpcPiFromXic0 = trackCharmBachelor.tpcNSigmaPi();
      float nSigTofPiFromXic0 = trackCharmBachelor.tofNSigmaPi();
      // Bachelor pion
      auto trackPionFromXi = casc.bachelor_as<TracksWCovDcaExtraPidPrPi>();
      float nSigTpcBachelorPi = trackPionFromXi.tpcNSigmaPi();
      float nSigTofBachelorPi = trackPionFromXi.tofNSigmaPi();
      // V0 daughters
      auto trackPosLambdaDaughter = casc.posTrack_as<TracksWCovDcaExtraPidPrPi>(); // If charge of cascade is negative, this should be track of proton
      auto trackNegLambdaDaughter = casc.negTrack_as<TracksWCovDcaExtraPidPrPi>(); // If charge of cascade is negative, this shoud be track of pion

      const auto& trackProtonFromLambda = chargeCasc < 0 ? trackPosLambdaDaughter : trackNegLambdaDaughter;
      const auto& trackPionFromLambda = chargeCasc < 0 ? trackNegLambdaDaughter : trackPosLambdaDaughter;

      float pPrFromLambda = trackProtonFromLambda.p();
      float nSigTpcPrFromLambda = trackProtonFromLambda.tpcNSigmaPr();
      float nSigTofPrFromLambda = trackProtonFromLambda.tofNSigmaPr();
      float pPiFromLambda = trackPionFromLambda.p();
      float nSigTpcPiFromLambda = trackPionFromLambda.tpcNSigmaPi();
      float nSigTofPiFromLambda = trackPionFromLambda.tofNSigmaPi();

      //------------------------------Fill QA histograms-----------------------------
      if (configs.fillHistograms) {
        registry.fill(HIST("hMassXic0Cand"), massXiPi);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
        registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
        registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
        registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);
        registry.fill(HIST("hCovSVXX"), covMatrixSV[0]);
        registry.fill(HIST("hCovSVYY"), covMatrixSV[2]);
        registry.fill(HIST("hCovSVXZ"), covMatrixSV[3]);
        registry.fill(HIST("hCovSVZZ"), covMatrixSV[5]);
        registry.fill(HIST("hDcaXYProngs"), trackCasc.getPt(), impactParameterCasc.getY());
        registry.fill(HIST("hDcaXYProngs"), trackCharmBachelor.pt(), impactParameterPi.getY());
        registry.fill(HIST("hDcaZProngs"), trackCasc.getPt(), impactParameterCasc.getZ());
        registry.fill(HIST("hImpParXiXY"), impactParameterCasc.getY());
        registry.fill(HIST("hImpParPiXY"), impactParameterPi.getY());
        registry.fill(HIST("hPtXi"), std::sqrt(pVecXi[0] * pVecXi[0] + pVecXi[1] * pVecXi[1])); // pt of Xi
        registry.fill(HIST("hPtPi"), std::sqrt(pVecPi[0] * pVecPi[0] + pVecPi[1] * pVecPi[1])); // pt of Pi
      }

      //------------------------------Fill the table-----------------------------
      cursors.rowCandXic0Base(
        /* Collision informations */
        collision.globalIndex(), collision.posX(), collision.posY(), collision.posZ(),
        std::sqrt(covMatrixPV[0]), std::sqrt(covMatrixPV[2]), std::sqrt(covMatrixPV[5]),
        /* 2-Prong specific columns */
        cand.cascadeId(), cand.prong0Id(),
        casc.bachelorId(), casc.posTrackId(), casc.negTrackId(),
        /* Secondary vertex*/
        secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
        std::sqrt(covMatrixSV[0]), std::sqrt(covMatrixSV[1]), std::sqrt(covMatrixSV[2]),
        /* Decay length error */
        errorDecayLength, errorDecayLengthXY,
        /* Chi2CPA, InvMass, cascade charge */
        chi2SV, massXiPi, chargeCasc,
        /* Cascade, charm bachelor's momentum */
        pVecXi[0], pVecXi[1], pVecXi[2],
        pVecPi[0], pVecPi[1], pVecPi[2],
        /* Impact parameter */
        impactParameterCasc.getY(), impactParameterPi.getY(),
        std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPi.getSigmaY2()),
        /* Cascade specific column */
        trackPionFromXi.p(), pPiFromLambda, pPrFromLambda,
        cpaXi, cpaXYXi,
        cpaLambda, cpaXYLambda,
        cpaLambdaToXi, cpaXYLambdaToXi,
        casc.mXi(), casc.mLambda(),
        /* DCA information*/
        casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(),
        casc.dcaXYCascToPV(), casc.dcaZCascToPV(),
        /* PID information */
        nSigTpcPiFromXic0, nSigTpcBachelorPi, nSigTpcPiFromLambda, nSigTpcPrFromLambda,
        nSigTofPiFromXic0, nSigTofBachelorPi, nSigTofPiFromLambda, nSigTofPrFromLambda);
    } // candidate loop
  }

  // template function for running xic0 reconstruction via KFParticle method
  // templated for various centrality estimator usage
  template <o2::hf_centrality::CentralityEstimator centEstimator, typename Colls>
  void runCreatorXic0WithKfParticle(Colls const&,
                                    aod::HfCascLf2Prongs const& candidates,
                                    KFCascadesLinked const&,
                                    KFCascFull const&,
                                    TracksWCovDcaExtraPidPrPi const&,
                                    aod::BCsWithTimestamps const&)
  {
    // Loop over candidates
    for (auto const& cand : candidates) {

      // Fill cascandidates before selection
      if (configs.fillHistograms) {
        registry.get<TH1>(HIST("hCandCounter"))->Fill(All);
      }

      // Apply hfflag selection
      if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
        continue;
      }

      // Event selection
      auto collision = cand.collision_as<Colls>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        continue;
      }

      //------------------------------Cascade pre-selection------------------------------

      // Retrieving skimmed cascade and pion tracks
      // If there is no related tracks, skip
      auto cascAodElement = cand.cascade_as<KFCascadesLinked>();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<KFCascFull>(); // -> Need to understand this
      auto chargeCasc = casc.sign() > 0 ? 1 : -1;

      if (configs.doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > configs.dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.mXi() - MassXiMinus) > configs.massToleranceCascade) {
          continue;
        }
      }

      if (configs.fillHistograms) {
        registry.fill(HIST("hCandCounter"), CascPreSel);
        registry.fill(HIST("hXiMassAfterConstrain"), casc.mXi());
      }

      //------------------------------Set Magnetic field------------------------------
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>> Current run Number : " << runNumber;
        initCCDB(bc, runNumber, ccdb, configs.isRun2 ? configs.ccdbPathGrp : configs.ccdbPathGrpMag, lut, configs.isRun2);
        bz = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>> Magnetic field: " << bz;
      }
      KFParticle::SetField(bz);

      //------------------------------Info of V0 and cascade tracks from LF table------------------------------
      // -> This quantities are used for physical properties of selected candidates
      // -> Not used for candidate creation
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};

      //------------------------------Create Xic0 as KF Particle object------------------------------

      // initialize primary vertex
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      float covMatrixPV[6];
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle kfPv(kfpVertex); // -> For calculation of DCAs to PV

      // convert charm bachelor pion tracks into KFParticle object
      auto trackCharmBachelor = cand.prong0_as<TracksWCovDcaExtraPidPrPi>();
      KFPTrack kfpTrackCharmBachelor = createKFPTrackFromTrack(trackCharmBachelor);
      KFParticle kfCharmBachelor(kfpTrackCharmBachelor, kPiPlus);

      // create Xi as KFParticle object
      constexpr std::size_t NElementsStateVector{6};
      std::array<float, NElementsStateVector> xyzpxpypz = {casc.x(), casc.y(), casc.z(), casc.px(), casc.py(), casc.pz()};
      float parPosMom[NElementsStateVector];
      std::copy(xyzpxpypz.begin(), xyzpxpypz.end(), parPosMom);

      KFParticle kfXi;
      float massXi = casc.mXi();
      kfXi.Create(parPosMom, casc.kfTrackCovMat(), casc.sign(), massXi);
      if (configs.useXiMassConstraint) {
        kfXi.SetNonlinearMassConstraint(MassXiMinus);
      }

      // create Xic0 as KFParticle object
      KFParticle kfXic0;
      const KFParticle* kfDaughterXic0[2] = {&kfCharmBachelor, &kfXi};
      kfXic0.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfXic0.Construct(kfDaughterXic0, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Xic0 : " << e.what();
        continue;
      }

      if (configs.fillHistograms) {
        registry.fill(HIST("hCandCounter"), VertexFit);
        registry.fill(HIST("hXic0Mass"), kfXic0.GetMass());
      }

      // get geometrical chi2 of xic0
      float chi2GeoXic0 = kfXic0.GetChi2() / kfXic0.GetNDF();

      // topologocal constraint of xic0 to PV
      float chi2topoXic0ToPVBeforeConstraint = kfXic0.GetDeviationFromVertex(kfPv);
      KFParticle kfXic0ToPV = kfXic0;
      kfXic0ToPV.SetProductionVertex(kfPv);
      float chi2topoXic0ToPV = kfXic0ToPV.GetChi2() / kfXic0ToPV.GetNDF();
      if (configs.constrainXic0ToPv) {
        kfXic0 = kfXic0ToPV;             // -> Replacing Xic0 with Xic0 propagated to PV
        kfXic0.TransportToDecayVertex(); // -> What does this do?
      }

      if (configs.fillHistograms) {
        registry.fill(HIST("hXic0MassAfterPVConstraint"), kfXic0.GetMass());
      }

      // topological constraint of Xi to Xic0
      float chi2topoXiToXic0BeforeConstraint = kfXi.GetDeviationFromVertex(kfXic0);
      KFParticle kfXiToXic0 = kfXi;
      kfXiToXic0.SetProductionVertex(kfXic0);
      float chi2topoXiToXic0 = kfXiToXic0.GetChi2() / kfXiToXic0.GetNDF();
      kfXiToXic0.TransportToDecayVertex(); // -> What does this do?
      if (configs.constrainXiToXic0) {
        KFParticle kfXic0WithXiToXic0;
        const KFParticle* kfDaughtersXic0WithXiToXic0[2] = {&kfCharmBachelor, &kfXiToXic0};
        kfXic0WithXiToXic0.SetConstructMethod(configs.kfConstructMethod);
        try {
          kfXic0WithXiToXic0.Construct(kfDaughtersXic0WithXiToXic0, 2);
        } catch (std::runtime_error& e) {
          LOG(debug) << "Failed to construct Xic0 with Xi constrained to Xic0: " << e.what();
          continue;
        }

        kfXic0 = kfXic0WithXiToXic0; // -> Replacing Xic0 propagated to PV with Xic0 with xi constrained to Xic0

        if (configs.fillHistograms) {
          registry.fill(HIST("hXic0MassAfterXiConstraint"), kfXic0.GetMass());
        }
      }

      // Get covariance matrix of xic0
      auto covMatrixXic0 = kfXic0.CovarianceMatrix();

      //------------------------------Calculate physical quantities and fill candidate table------------------------------

      // transport Xic0 daughters to Xic0 decay vertex
      float secondaryVertex[3] = {kfXic0.GetX(), kfXic0.GetY(), kfXic0.GetZ()};
      kfXi.TransportToPoint(secondaryVertex);
      kfCharmBachelor.TransportToPoint(secondaryVertex);

      // impact parameters(daughters~PV) of xic0 daughters
      float impactParameterPiXY = 0., errImpactParameterPiXY = 0.;
      float impactParameterXiXY = 0., errImpactParameterXiXY = 0.;
      kfCharmBachelor.GetDistanceFromVertexXY(kfPv, impactParameterPiXY, errImpactParameterPiXY);
      kfXi.GetDistanceFromVertexXY(kfPv, impactParameterXiXY, errImpactParameterXiXY);

      // calculate cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      float cpaLambda = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYLambda = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaXi = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaXYXi = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float cpaLambdaToXi = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float cpaXYLambdaToXi = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

      // get DCAs of bachelor pion and cascade
      float dcaXYPiXi = kfCharmBachelor.GetDistanceFromParticleXY(kfXi);
      float dcaPiXi = kfCharmBachelor.GetDistanceFromParticle(kfXi);

      // get invariant mass of Xic0 candidate
      float errMassXiPi;
      kfXic0.GetMass(massXiPi, errMassXiPi);

      // decay length of Xic0
      // use Xic0 constrainedto PV (production point must be set before calling GetDecayLength(XY) on KFParticle)
      float kfDecayLength = 0., errKfDecayLength = 0., kfDecayLengthXY = 0., errKfDecayLengthXY = 0.;
      kfXic0ToPV.GetDecayLength(kfDecayLength, errKfDecayLength);
      kfXic0ToPV.GetDecayLengthXY(kfDecayLengthXY, errKfDecayLengthXY);
      float kfDecayLengthNormalised = ldlFromKF(kfXic0, kfPv);
      float kfDecayLengthXYNormalised = ldlXYFromKF(kfXic0, kfPv);

      //-----Get PID information-----
      // Charm Bachelor pion
      float nSigTpcPiFromXic0 = trackCharmBachelor.tpcNSigmaPi();
      float nSigTofPiFromXic0 = trackCharmBachelor.tofNSigmaPi();
      // Bachelor pion(pion from cascade decay)
      auto trackPionFromXi = casc.bachelor_as<TracksWCovDcaExtraPidPrPi>();
      float nSigTpcBachelorPi = trackPionFromXi.tpcNSigmaPi();
      float nSigTofBachelorPi = trackPionFromXi.tofNSigmaPi();
      // V0 daughters
      auto trackPosLambdaDaughter = casc.posTrack_as<TracksWCovDcaExtraPidPrPi>(); // If charge of cascade is negative, this should be track of proton
      auto trackNegLambdaDaughter = casc.negTrack_as<TracksWCovDcaExtraPidPrPi>(); // If charge of cascade is negative, this shoud be track of pion

      const auto& trackProtonFromLambda = chargeCasc < 0 ? trackPosLambdaDaughter : trackNegLambdaDaughter;
      const auto& trackPionFromLambda = chargeCasc < 0 ? trackNegLambdaDaughter : trackPosLambdaDaughter;

      float pPrFromLambda = trackProtonFromLambda.p();
      float nSigTpcPrFromLambda = trackProtonFromLambda.tpcNSigmaPr();
      float nSigTofPrFromLambda = trackProtonFromLambda.tofNSigmaPr();
      float pPiFromLambda = trackPionFromLambda.p();
      float nSigTpcPiFromLambda = trackPionFromLambda.tpcNSigmaPi();
      float nSigTofPiFromLambda = trackPionFromLambda.tofNSigmaPi();

      //------------------------------Calculate physical quantities and fill candidate table------------------------------
      if (configs.fillHistograms) {
        // inv mass
        registry.fill(HIST("hMassXic0Cand"), massXiPi);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);
        registry.fill(HIST("hCovPVYY"), covMatrixPV[2]);
        registry.fill(HIST("hCovPVXZ"), covMatrixPV[3]);
        registry.fill(HIST("hCovPVZZ"), covMatrixPV[5]);
        registry.fill(HIST("hCovSVXX"), covMatrixXic0[0]);
        registry.fill(HIST("hCovSVYY"), covMatrixXic0[2]);
        registry.fill(HIST("hCovSVXZ"), covMatrixXic0[3]);
        registry.fill(HIST("hCovSVZZ"), covMatrixXic0[5]);
        registry.fill(HIST("hDcaXYProngs"), kfXi.GetPt(), impactParameterXiXY);
        registry.fill(HIST("hDcaXYProngs"), kfCharmBachelor.GetPt(), impactParameterPiXY);
        registry.fill(HIST("hImpParXiXY"), impactParameterXiXY);
        registry.fill(HIST("hImpParPiXY"), impactParameterPiXY);
        registry.fill(HIST("hPtXi"), kfXi.GetPt());            // pt of Xi
        registry.fill(HIST("hPtPi"), kfCharmBachelor.GetPt()); // pt of Pi
      }

      //------------------------------Fill the table------------------------------
      cursors.rowCandXic0Base(
        /* Collision information */
        collision.globalIndex(),
        collision.posX(), collision.posY(), collision.posZ(),
        std::sqrt(covMatrixPV[0]), std::sqrt(covMatrixPV[2]), std::sqrt(covMatrixPV[5]),
        /*2-Prong specific columns*/
        cand.cascadeId(), cand.prong0Id(),                       // -> xi-, pi+ from xic0 decay
        casc.bachelorId(), casc.posTrackId(), casc.negTrackId(), // -> pi+ from xi- decay, proton from lambda0 decay, pion- from lambda0 decay
        /*Secondary vertex*/
        secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
        kfXic0.GetErrX(), kfXic0.GetErrY(), kfXic0.GetErrZ(),
        /*DecayLength error*/
        errKfDecayLength, errKfDecayLengthXY,
        /*Chi2 of Geo from KF method, Invmass, cascade charge*/
        chi2GeoXic0, massXiPi, chargeCasc,
        /*Cascade, charm bachelor's momentum*/
        kfXi.GetPx(), kfXi.GetPy(), kfXi.GetPz(),
        kfCharmBachelor.GetPx(), kfCharmBachelor.GetPy(), kfCharmBachelor.GetPz(),
        /*Impact parameter*/
        impactParameterXiXY, impactParameterPiXY,
        errImpactParameterXiXY, errImpactParameterPiXY,
        /*Cascade specific column*/
        trackPionFromXi.p(),
        pPiFromLambda,
        pPrFromLambda,
        cpaXi, cpaXYXi,
        cpaLambda, cpaXYLambda,
        cpaLambdaToXi, cpaXYLambdaToXi,
        massXi, casc.mLambda(),
        /*DCA information*/
        casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(),
        casc.dcaXYCascToPV(), casc.dcaZCascToPV(),
        /*PID information*/
        nSigTpcPiFromXic0, nSigTpcBachelorPi, nSigTpcPiFromLambda, nSigTpcPrFromLambda,
        nSigTofPiFromXic0, nSigTofBachelorPi, nSigTofPiFromLambda, nSigTofPrFromLambda);

      cursors.rowCandXic0KF(
        casc.kfCascadeChi2(), casc.kfV0Chi2(),
        kfDecayLength, kfDecayLengthNormalised, kfDecayLengthXY, kfDecayLengthXYNormalised,
        chi2topoXic0ToPVBeforeConstraint, chi2topoXic0ToPV, chi2topoXiToXic0BeforeConstraint, chi2topoXiToXic0,
        dcaXYPiXi, dcaPiXi);

    } // end candidate loop

  } // end of runCreator

  ////////////////////////////////////////////////////////
  ///                                                  ///
  ///           Process functions with DCAFitter       ///
  ///                                                  ///
  ////////////////////////////////////////////////////////

  void processXic0WithDCAFitterNoCent(SelectedCollisions const& collisions,
                                      aod::HfCascLf2Prongs const& candidates,
                                      CascadesLinked const& cascadesLinked,
                                      CascFull const& cascFull,
                                      TracksWCovDcaExtraPidPrPi const& tracks,
                                      aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorXic0WithDCAFitter<CentralityEstimator::None>(collisions, candidates, cascadesLinked, cascFull, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processXic0WithDCAFitterNoCent, "Xic0 reconstruction via DcaFitter method, no centrality", true);

  void processXic0WithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                        aod::HfCascLf2Prongs const& candidates,
                                        CascadesLinked const& cascadesLinked,
                                        CascFull const& cascFull,
                                        TracksWCovDcaExtraPidPrPi const& tracks,
                                        aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorXic0WithDCAFitter<CentralityEstimator::FT0C>(collisions, candidates, cascadesLinked, cascFull, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processXic0WithDCAFitterCentFT0C, "Xic0 reconstruction via DcaFitter method, centrality selection on FT0C", false);

  void processXic0WithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                        aod::HfCascLf2Prongs const& candidates,
                                        CascadesLinked const& cascadesLinked,
                                        CascFull const& cascFull,
                                        TracksWCovDcaExtraPidPrPi const& tracks,
                                        aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorXic0WithDCAFitter<CentralityEstimator::FT0M>(collisions, candidates, cascadesLinked, cascFull, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processXic0WithDCAFitterCentFT0M, "Xic0 reconstruction via DcaFitter method, centrality selection on FT0M", false);

  ////////////////////////////////////////////////////////
  ///                                                  ///
  ///           Process functions with KFParticle      ///
  ///                                                  ///
  ////////////////////////////////////////////////////////

  void processXic0WithKFParticleNoCent(SelectedCollisions const& collisions,
                                       aod::HfCascLf2Prongs const& candidates,
                                       KFCascadesLinked const& kfCascadesLinked,
                                       KFCascFull const& kfCascFull,
                                       TracksWCovDcaExtraPidPrPi const& tracks,
                                       aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorXic0WithKfParticle<CentralityEstimator::None>(collisions, candidates, kfCascadesLinked, kfCascFull, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processXic0WithKFParticleNoCent, "Xic0 reconstruction via KFParticle method, no centrality", false);

  void processXic0WithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                         aod::HfCascLf2Prongs const& candidates,
                                         KFCascadesLinked const& kfCascadesLinked,
                                         KFCascFull const& kfCascFull,
                                         TracksWCovDcaExtraPidPrPi const& tracks,
                                         aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorXic0WithKfParticle<CentralityEstimator::FT0C>(collisions, candidates, kfCascadesLinked, kfCascFull, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processXic0WithKFParticleCentFT0C, "Xic0 reconstruction via KFParticle method, centrality on FT0C", false);

  void processXic0WithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                         aod::HfCascLf2Prongs const& candidates,
                                         KFCascadesLinked const& kfCascadesLinked,
                                         KFCascFull const& kfCascFull,
                                         TracksWCovDcaExtraPidPrPi const& tracks,
                                         aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorXic0WithKfParticle<CentralityEstimator::FT0M>(collisions, candidates, kfCascadesLinked, kfCascFull, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processXic0WithKFParticleCentFT0M, "Xic0 reconstruction via KFParticle method, centrality on FT0M", false);

  ///////////////////////////////////////////////////////////////////
  ///                                                             ///
  ///           Process functions for Collision monitoring        ///
  ///                                                             ///
  ///////////////////////////////////////////////////////////////////

  void processCollisionsNoCent(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                               aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {

      // bitmask with event selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      // monitor the satisfied event selection
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processCollisionsNoCent, "Collision monitoring - No Centrality", false);

  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                 aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {

      // bitmask with event selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      // monitor the satisfied event selection
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processCollisionsCentFT0C, "Collision monitoring - Centrality selection with FT0C", false);

  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                 aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {

      // bitmask with event selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      // monitor the satisfied event selection
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronic, processCollisionsCentFT0M, "Collision monitoring - Centrality selection with FT0M", false);
};

struct HfCandidateCreatorXic0XicpToHadronicMc {

  Spawns<aod::HfCandXic0Ext> rowCandXic0Ext;

  struct : ProducesGroup {

    Produces<aod::HfCandXic0McRec> rowCandXic0McRec;
    Produces<aod::HfCandXic0McGen> rowCandXic0McGen;

  } cursors;

  struct : ConfigurableGroup {

    Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"}; // -> Used for only Xic0
    Configurable<bool> acceptTrackInteractionWithMaterial{"acceptTrackInteractionWithMaterial", false, "Accept candidates with final daughters interacting with materials"};
    Configurable<bool> fillMcHistograms{"fillMcHistograms", true, "Fill validation plots"};
    Configurable<bool> fillResidualTable{"fillResidualTable", true, "Fill table contaning residuals and pulls of PV and SV"};
    Configurable<bool> matchDecayedPions{"matchedDecayedPions", true, "Match also candidates with daughter pion tracks that decay with kinked toploogy"};

  } configs;

  enum McMatchFlag : uint8_t {
    None = 0,
    CharmBaryonUnmatched,
    CascUnmatched,
    V0Unmatched,
    NumberMcMatchFlag
  };

  // Table aliases
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>; // -> Used for subscription for process functions of centrality with FT0Ms
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId; // -> Why use unsorted??
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};
  HfEventSelectionMc hfEvSelMc;

  void init(InitContext& initContext)
  {
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-xic0-xicp-to-hadronic") == 0) {
        hfEvSelMc.init(device, registry);
        break;
      }
    }

    // Add histograms for QA
    if (configs.fillMcHistograms) {
      registry.add("hDebugStatusRec", "hDebugStatusRec", {HistType::kTH1D, {{NumberMcMatchFlag, 0, NumberMcMatchFlag}}});
      registry.add("hDebugStatusGen", "hDebugStatusGen", {HistType::kTH1D, {{NumberMcMatchFlag, 0, NumberMcMatchFlag}}});
      TString labels[McMatchFlag::NumberMcMatchFlag];
      labels[McMatchFlag::None] = "None";
      labels[McMatchFlag::CharmBaryonUnmatched] = "CharmBaryonUnmatched";
      labels[McMatchFlag::CascUnmatched] = "CascUnmatched";
      labels[McMatchFlag::V0Unmatched] = "V0Unmatched";
      for (int iBin = 0; iBin < McMatchFlag::NumberMcMatchFlag; iBin++) {
        registry.get<TH1>(HIST("hDebugStatusRec"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
        registry.get<TH1>(HIST("hDebugStatusGen"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
      }
    }
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, typename McCollisions, typename CollInfos>
  void runXic0Mc(aod::HfCandXic0Base const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
                 McCollisions const& mcCollisions,
                 CollInfos const& collInfos,
                 BCsInfo const&)
  {
    int indexRec{-1};
    int indexRecXic0{-1};
    int8_t sign{0};
    int8_t signCasc{0};
    int8_t signV0{0};
    int8_t flag{0};
    int8_t origin{0};
    int8_t debug{0};

    ////////////////////////////////////
    // Match reconstructed candidates //
    ////////////////////////////////////

    for (const auto& candidate : candidates) {

      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = McMatchFlag::None;
      std::vector<int> idxBhadMothers{};

      auto arrayDaughters = std::array{candidate.pi_as<aod::TracksWMc>(),        // pi <- Xic0
                                       candidate.bachelor_as<aod::TracksWMc>(),  // pi <- Xi-
                                       candidate.posTrack_as<aod::TracksWMc>(),  // pr <- lambda0 <- xi-
                                       candidate.negTrack_as<aod::TracksWMc>()}; // pi <- lambda0 <- xi-

      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};

      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Reject particles from background events
      if (configs.rejectBackground) {
        bool fromBkg{false};
        for (auto const& daughter : arrayDaughters) {
          if (daughter.has_mcParticle()) {
            auto mcParticle = daughter.mcParticle();
            if (mcParticle.fromBackgroundEvent()) {
              fromBkg = true;
              break;
            }
          }
        } // daughter loop
        if (fromBkg) {
          // fill the tables
          cursors.rowCandXic0McRec(0, McMatchFlag::None, RecoDecay::OriginType::None);
          registry.fill(HIST("hDebugStatusRec"), McMatchFlag::None);
          continue;
        }
      }

      // !Xic0 -> Xi Pi matching
      // Xic0 -> ((pi- p) pi-) pi+
      indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, +kXiC0, std::array{+kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 3);
      indexRecXic0 = indexRec;
      if (indexRec == -1) { // Xic0 not reconstructed
        debug = McMatchFlag::CharmBaryonUnmatched;
        registry.fill(HIST("hDebugStatusRec"), static_cast<double>(McMatchFlag::CharmBaryonUnmatched));
      }
      if (indexRec > -1) {
        // xi- -> pi pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiPlus, +kProton, +kPiMinus}, true, &signCasc, 2);
        if (indexRec == -1) { // Xi- not reconstructed
          debug = McMatchFlag::CascUnmatched;
          registry.fill(HIST("hDebugStatusRec"), static_cast<double>(McMatchFlag::CascUnmatched));
        }
        if (indexRec > -1) {
          // Lambda -> p pi
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);
          if (indexRec == -1) { // V0 not reconstructed
            debug = McMatchFlag::V0Unmatched;
            registry.fill(HIST("hDebugStatusRec"), static_cast<double>(McMatchFlag::V0Unmatched));
          }
          if (indexRec > -1) {
            flag = sign * (1 << aod::hf_cand_xic0_xicp_to_hadronic::DecayTypeXic0::Xic0ToXiPi);
          }
        }
      }

      // Check if Xic0 is from b-hadron decay(prompt vs non-prompt)
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRecXic0);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }

      // Fill the table
      cursors.rowCandXic0McRec(flag, debug, origin);
    } // candidate loop

    ///////////////////////////////
    // Match generated particles //
    ///////////////////////////////

    for (auto const& mcCollision : mcCollisions) {

      auto const mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());

      float centrality{-1.f};
      uint16_t rejectionMask{0};
      int nSplitColl{0};

      if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::None) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
        const auto collSlice = collInfos.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        nSplitColl = collSlice.size();
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      }

      hfEvSelMc.fillHistograms<centEstimator>(mcCollision, rejectionMask, nSplitColl);

      if (rejectionMask != 0) { // none of the event selection was satisfied(?) -> Reject all particles from this event
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          cursors.rowCandXic0McGen(0, McMatchFlag::None, RecoDecay::OriginType::None);
        }
        continue;
      }

      for (auto const& particle : mcParticlesPerMcColl) {
        flag = 0;
        sign = 0;
        debug = 0;
        origin = RecoDecay::OriginType::None;

        // Xic0 -> Xi- pi+
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, +kXiC0, std::array{+kXiMinus, +kPiPlus}, true, &sign)) {

          debug = 1; // -> Matched Xic0

          for (auto const& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != +kXiMinus) {
              continue;
            }
            // Xi -> Lambda + pi
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, +kXiMinus, std::array{+kLambda0, +kPiMinus}, true)) {
              debug = 2; // -> Matched Xi-
              for (auto const& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode() != +kLambda0)) {
                  continue;
                }

                // Lambda -> p + pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
                  debug = 3; // -> Matched Lambda0
                  flag = sign * (1 << o2::aod::hf_cand_xic0_xicp_to_hadronic::DecayTypeXic0::Xic0ToXiPi);
                }
              } // V0 daughter loop
            } // cascade daughter loop
          }
        } // charm daughter loop

        // Check if charm is prompt or non-prompt
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false);
        }
        // Fill the table
        cursors.rowCandXic0McGen(flag, debug, origin);

      } // particle loop

    } // end of collision loop

  } // template run function

  //////////////////////////////////
  ///                            ///
  ///      Process functions     ///
  ///                            ///
  //////////////////////////////////

  void processMcEmpty(aod::Collisions const&)
  {
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronicMc, processMcEmpty, "Empty process function to prevent workflow from getting stuck", true);

  //~~~~Xic0~~~~//
  void processMcXic0(aod::HfCandXic0Base const& candidates,
                     aod::TracksWMc const& tracks,
                     aod::McParticles const& mcParticles,
                     aod::McCollisions const& mcCollisions,
                     McCollisionsNoCents const& mcCollisionsNoCents,
                     BCsInfo const& bcs)
  {
    runXic0Mc<o2::hf_centrality::CentralityEstimator::None>(candidates, tracks, mcParticles, mcCollisions, mcCollisionsNoCents, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronicMc, processMcXic0, "Perform MC matching of Xic0, no cents", false);

  void processMcXic0CentFT0C(aod::HfCandXic0Base const& candidates,
                             aod::TracksWMc const& tracks,
                             aod::McParticles const& mcParticles,
                             aod::McCollisions const& mcCollisions,
                             McCollisionsFT0Cs const& mcCollisionsFT0Cs,
                             BCsInfo const& bcs)
  {
    runXic0Mc<o2::hf_centrality::CentralityEstimator::FT0C>(candidates, tracks, mcParticles, mcCollisions, mcCollisionsFT0Cs, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronicMc, processMcXic0CentFT0C, "Perform MC matching of Xic0, centrality on FT0C", false);

  void processMcXic0CentFT0M(aod::HfCandXic0Base const& candidates,
                             aod::TracksWMc const& tracks,
                             aod::McParticles const& mcParticles,
                             McCollisionsCentFT0Ms const& mcCollisions,
                             McCollisionsFT0Ms const& mcCollisionsFT0Ms,
                             BCsInfo const& bcs)
  {
    runXic0Mc<o2::hf_centrality::CentralityEstimator::FT0M>(candidates, tracks, mcParticles, mcCollisions, mcCollisionsFT0Ms, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0XicpToHadronicMc, processMcXic0CentFT0M, "Perform MC matching of Xic0, centrality on FT0M", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic0XicpToHadronic>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXic0XicpToHadronicMc>(cfgc)};
}
