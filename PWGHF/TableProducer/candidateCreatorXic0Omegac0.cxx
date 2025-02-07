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

/// \file candidateCreatorXic0Omegac0.cxx
/// \brief Reconstruction of Omegac0 and Xic0 decays candidates
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Ruiqi Yin <ruiqi.yin@cern.ch>, Fudan University
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include <iterator>
#include <memory>
#include <string>
#include <vector>

/// includes KFParticle
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/KFparticle/KFUtilities.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::track;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::aod::hf_track_index;
using namespace o2::hf_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_evsel;

// Reconstruction of omegac0 and xic0 candidates
struct HfCandidateCreatorXic0Omegac0 {
  Produces<aod::HfCandToXiPi> rowCandToXiPi;
  Produces<aod::HfCandToOmegaPi> rowCandToOmegaPi;
  Produces<aod::HfCandToOmegaK> rowCandToOmegaK;
  Produces<aod::HfOmegacKf> kfCandidateData;
  Produces<aod::HfCandToXiPiKf> kfCandidateXicData;

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
  Configurable<bool> fillAllHist{"fillAllHist", true, "Fill additional KF histograms to check selector cuts"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // KFParticle process setting
  //  V0 cuts
  Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.0075, "Distance from Lambda mass"};
  // cascade cuts
  Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascade"};
  // for KF particle operation
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
  Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", false, "KF: use Lambda mass constraint"};
  Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: use Cascade mass constraint"};

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter to build the omegac/xic vertex
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{-1};
  double magneticField{0.};

  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>;
  using MyTraCascTable = soa::Join<aod::TraCascDatas, aod::TraCascCovs>; // to use strangeness tracking
  using CascadesLinked = soa::Join<Cascades, CascDataLink>;
  using TraCascadesLinked = soa::Join<Cascades, TraCascDataLink>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MyLFTracksWCov = soa::Join<TracksIU, TracksCovIU>;

  using MyKfTracks = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using MyKfCascTable = soa::Join<KFCascDatas, aod::KFCascCovs>;
  using KFCascadesLinked = soa::Join<aod::Cascades, aod::KFCascDataLink>;

  std::shared_ptr<TH1> hInvMassCharmBaryonToXiPi, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaK, hFitterStatusToXiPi, hFitterStatusToOmegaPi, hFitterStatusToOmegaK, hCandidateCounterToXiPi, hCandidateCounterToOmegaPi, hCandidateCounterToOmegaK, hCascadesCounterToXiPi, hCascadesCounterToOmegaPi, hCascadesCounterToOmegaK;

  HistogramRegistry registry{"registry"};
  // Helper struct to pass  information
  struct {
    float chi2GeoV0;
    float ldlV0;
    float chi2NdfTopoV0ToPv;
    float chi2GeoCasc;
    float ldlCasc;
    float chi2NdfTopoCascToPv;
    float decayLenXYLambda;
    float decayLenXYCasc;
    float cosPaV0ToCasc;
    float cosPaXYV0ToCasc;
    float cosPaV0ToPv;
    float cosPaXYV0ToPv;
    float cosPaCascToOmegac;
    float cosPaXYCascToOmegac;
    float cosPaCascToPv;
    float cosPaXYCascToPv;
    float massV0;
    float massCasc;
    float ptPiFromOmegac;
    float ptOmegac;
    float rapOmegac;
    float massOmegac;
    float cosThetaStarPiFromOmegac;
    float chi2NdfTopoPiFromOmegacToPv;
    float kfDcaXYPiFromOmegac;
    float chi2NdfTopoV0ToCasc;
    float chi2NdfTopoCascToOmegac;
    float decayLenXYOmegac;
    float chi2GeoOmegac;
    float kfDcaV0Dau;
    float kfDcaCascDau;
    float kfDcaOmegacDau;
    float kfDcaXYCascToPv;
    float chi2NdfTopoOmegacToPv;
    float cosPaOmegacToPv;
    float cosPaXYOmegacToPv;
    float ldlOmegac;
    float ctV0;
    float ctCasc;
    float ctOmegac;
    float chi2MassV0;
    float chi2MassCasc;
    float etaOmegac;
    float cascRejectInvmass; // rej
  } kfOmegac0Candidate;

  struct {
    float chi2GeoV0;
    float ldlV0;
    float chi2TopoV0ToPv;
    float chi2GeoCasc;
    float ldlCasc;
    float chi2TopoCascToPv;
    float decayLenXYLambda;
    float decayLenXYCasc;
    float cosPaV0ToCasc;
    float cosPaXYV0ToCasc;
    float cosPaV0ToPv;
    float cosPaXYV0ToPv;
    float cosPaCascToXic;
    float cosPaXYCascToXic;
    float cosPaCascToPv;
    float cosPaXYCascToPv;
    float massV0;
    float massCasc;
    float ptPiFromXic;
    float ptXic;
    float rapXic;
    float massXic;
    float cosThetaStarPiFromXic;
    float chi2TopoPiFromXicToPv;
    float kfDcaXYPiFromXic;
    float chi2TopoV0ToCasc;
    float chi2TopoCascToXic;
    float decayLenXYXic;
    float chi2GeoXic;
    float kfDcaV0Dau;
    float kfDcaCascDau;
    float kfDcaXicDau;
    float kfDcaXYCascToPv;
    float chi2TopoXicToPv;
    float cosPaXicToPv;
    float cosPaXYXicToPv;
    float ldlXic;
    float ctV0;
    float ctCasc;
    float ctXic;
    float ctOmegac;
    float chi2MassV0;
    float chi2MassCasc;
    float etaXic;
  } kfXic0Candidate;
  void init(InitContext const&)
  {
    std::array<bool, 12> allProcesses = {doprocessNoCentToXiPi, doprocessNoCentToXiPiTraCasc, doprocessCentFT0CToXiPi, doprocessCentFT0MToXiPi, doprocessNoCentToOmegaPi, doprocessOmegacToOmegaPiWithKFParticle, doprocessCentFT0CToOmegaPi, doprocessCentFT0MToOmegaPi, doprocessNoCentToOmegaK, doprocessCentFT0CToOmegaK, doprocessCentFT0MToOmegaK, doprocessXicToXiPiWithKFParticle};
    if (std::accumulate(allProcesses.begin(), allProcesses.end(), 0) == 0) {
      LOGP(fatal, "No process function enabled, please select one for at least one channel.");
    }

    std::array<bool, 5> processesToXiPi = {doprocessNoCentToXiPi, doprocessNoCentToXiPiTraCasc, doprocessCentFT0CToXiPi, doprocessCentFT0MToXiPi, doprocessXicToXiPiWithKFParticle};
    if (std::accumulate(processesToXiPi.begin(), processesToXiPi.end(), 0) > 1) {
      LOGP(fatal, "One and only one ToXiPi process function must be enabled at a time.");
    }
    std::array<bool, 4> processesToOmegaPi = {doprocessNoCentToOmegaPi, doprocessCentFT0CToOmegaPi, doprocessCentFT0MToOmegaPi, doprocessOmegacToOmegaPiWithKFParticle};
    if (std::accumulate(processesToOmegaPi.begin(), processesToOmegaPi.end(), 0) > 1) {
      LOGP(fatal, "One and only one process ToOmegaPi function must be enabled at a time.");
    }
    std::array<bool, 3> processesToOmegaK = {doprocessNoCentToOmegaK, doprocessCentFT0CToOmegaK, doprocessCentFT0MToOmegaK};
    if (std::accumulate(processesToOmegaK.begin(), processesToOmegaK.end(), 0) > 1) {
      LOGP(fatal, "One and only one process ToOmegaK function must be enabled at a time.");
    }

    std::array<bool, 3> processesCollisions = {doprocessCollisions, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};
    const int nProcessesCollisions = std::accumulate(processesCollisions.begin(), processesCollisions.end(), 0);
    if (nProcessesCollisions > 1) {
      LOGP(fatal, "At most one process function for collision monitoring can be enabled at a time.");
    }
    if (nProcessesCollisions == 1) {
      if ((doprocessNoCentToXiPi && !doprocessCollisions) || (doprocessNoCentToXiPiTraCasc && !doprocessCollisions) || (doprocessNoCentToOmegaPi && !doprocessCollisions) || (doprocessNoCentToOmegaK && !doprocessCollisions) || (doprocessOmegacToOmegaPiWithKFParticle && !doprocessCollisions) || (doprocessXicToXiPiWithKFParticle && !doprocessCollisions)) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisions\"?");
      }
      if ((doprocessCentFT0CToXiPi && !doprocessCollisionsCentFT0C) || (doprocessCentFT0CToOmegaPi && !doprocessCollisionsCentFT0C) || (doprocessCentFT0CToOmegaK && !doprocessCollisionsCentFT0C)) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0C\"?");
      }
      if ((doprocessCentFT0MToXiPi && !doprocessCollisionsCentFT0M) || (doprocessCentFT0MToOmegaPi && !doprocessCollisionsCentFT0M) || (doprocessCentFT0MToOmegaK && !doprocessCollisionsCentFT0M)) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0M\"?");
      }
    }

    hInvMassCharmBaryonToXiPi = registry.add<TH1>("hInvMassCharmBaryonToXiPi", "Charm baryon invariant mass - #Xi #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2.2, 3.1}}});
    hInvMassCharmBaryonToOmegaPi = registry.add<TH1>("hInvMassCharmBaryonToOmegaPi", "Charm baryon invariant mass - #Omega #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2.2, 3.1}}});
    hInvMassCharmBaryonToOmegaK = registry.add<TH1>("hInvMassCharmBaryonToOmegaK", "Charm baryon invariant mass - #Omega K decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 2.2, 3.1}}});
    hFitterStatusToXiPi = registry.add<TH1>("hFitterStatusToXiPi", "Charm DCAFitter status - #Xi #pi vtx;status;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});                                // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    hFitterStatusToOmegaPi = registry.add<TH1>("hFitterStatusToOmegaPi", "Charm DCAFitter status - #Omega #pi vtx ;status;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});                      // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    hFitterStatusToOmegaK = registry.add<TH1>("hFitterStatusToOmegaK", "Charm DCAFitter status - #Omega K vtx;status;entries", {HistType::kTH1D, {{3, -0.5, 2.5}}});                           // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
    hCandidateCounterToXiPi = registry.add<TH1>("hCandidateCounterToXiPi", "Candidate counter wrt derived data - #Xi #pi decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});          // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    hCandidateCounterToOmegaPi = registry.add<TH1>("hCandidateCounterToOmegaPi", "Candidate counter wrt derived data - #Omega #pi decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}}); // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    hCandidateCounterToOmegaK = registry.add<TH1>("hCandidateCounterToOmegaK", "Candidate counter wrt derived data - #Omega K decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});     // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table
    hCascadesCounterToXiPi = registry.add<TH1>("hCascadesCounterToXiPi", "Cascades counter wrt derived data - #Xi #pi decay;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});             // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table
    hCascadesCounterToOmegaPi = registry.add<TH1>("hCascadesCounterToOmegaPi", "Cascades counter wrt derived data - #Omega #pi decay;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});    // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table
    hCascadesCounterToOmegaK = registry.add<TH1>("hCascadesCounterToOmegaK", "Cascades counter wrt derived data - #Omega K decay;status;entries", {HistType::kTH1D, {{2, -0.5, 1.5}}});        // 0 --> cascades in derived data table (and stored in AOD table), 1 --> cascades in derived data table and also accessible in cascData table

    // KFParticle Variables Histograms
    registry.add("hKFParticleV0TopoChi2", "hKFParticleV0TopoChi2", kTH1D, {{1000, -0.10f, 100.0f}});
    registry.add("hKFParticleCascTopoChi2", "hKFParticleCascTopoChi2", kTH1D, {{1000, -0.1f, 100.0f}});
    registry.add("hKFParticleCascBachTopoChi2", "hKFParticleCascBachTopoChi2", kTH1D, {{1000, -0.1f, 100.0f}});
    registry.add("hKFParticleDcaCharmBaryonDau", "hKFParticleDcaCharmBaryonDau", kTH1D, {{1000, -0.1f, 1.0f}});
    registry.add("hKFParticleDcaXYCascBachToPv", "hKFParticleDcaXYCascBachToPv", kTH1D, {{1000, -0.1f, 15.0f}});
    registry.add("hKfLambda_ldl", "hKfLambda_ldl", kTH1D, {{1000, 0.0f, 1000.0f}});
    registry.add("hKfOmega_ldl", "hKfOmega_ldl", kTH1D, {{1000, 0.0f, 1000.0f}});
    registry.add("hKfXi_ldl", "hKfXi_ldl", kTH1D, {{1000, 0.0f, 1000.0f}});
    registry.add("hKfOmegaC0_ldl", "hKfOmegaC0_ldl", kTH1D, {{1000, 0.0f, 1000.0f}});
    registry.add("hKfXiC0_ldl", "hKfXiC0_ldl", kTH1D, {{1000, 0.0f, 1000.0f}});
    registry.add("hDcaXYCascadeToPVKf", "hDcaXYCascadeToPVKf", kTH1D, {{1000, 0.0f, 2.0f}});
    registry.add("hInvMassOmegaMinus", "hInvMassOmegaMinus", kTH1D, {{1000, 1.6f, 2.0f}});
    registry.add("hInvMassXiMinus", "hInvMassXiMinus", kTH1D, {{1000, 1.25f, 1.65f}});
    registry.add("hInvMassXiMinus_rej", "hInvMassXiMinus_rej", kTH1D, {{1000, 1.25f, 1.65f}});

    // Additional KFParticle Histograms
    registry.add("hKFParticlechi2TopoOmegacToPv", "hKFParticlechi2TopoOmegacToPv", kTH1D, {{1000, -0.1f, 100.0f}});
    registry.add("hKFParticlechi2TopoCascToPv", "hKFParticlechi2TopoCascToPv", kTH1D, {{1000, -0.1f, 100.0f}});
    registry.add("hKFParticleDcaXYV0DauPosToPv", "hKFParticleDcaXYV0DauPosToPv", kTH1D, {{1000, -0.1f, 30.0f}});
    registry.add("hKFParticleDcaXYV0DauNegToPv", "hKFParticleDcaXYV0DauNegToPv", kTH1D, {{1000, -0.1f, 30.0f}});
    registry.add("hEtaV0PosDau", "hEtaV0PosDau", kTH1D, {{1000, -5.0f, 5.0f}});
    registry.add("hEtaV0NegDau", "hEtaV0NegDau", kTH1D, {{1000, -5.0f, 5.0f}});
    registry.add("hEtaKaFromCasc", "hEtaKaFromCasc", kTH1D, {{1000, -5.0f, 5.0f}});
    registry.add("hEtaPiFromCharmBaryon", "hEtaPiFromCharmBaryon", kTH1D, {{1000, -5.0f, 5.0f}});
    registry.add("hCascradius", "hCascradius", kTH1D, {{1000, 0.0f, 50.0f}});
    registry.add("hV0radius", "hV0radius", kTH1D, {{1000, 0.0f, 50.0f}});
    registry.add("hCosPACasc", "hCosPACasc", kTH1D, {{5000, 0.8f, 1.1f}});
    registry.add("hCosPAV0", "hCosPAV0", kTH1D, {{5000, 0.8f, 1.1f}});
    registry.add("hDcaCascDau", "hDcaCascDau", kTH1D, {{1000, -0.1f, 10.0f}});
    registry.add("hDcaV0Dau", "hDcaV0Dau", kTH1D, {{1000, -0.1f, 10.0f}});
    registry.add("hDcaXYToPvKa", "hDcaXYToPvKa", kTH1D, {{1000, -0.1f, 10.0f}});
    registry.add("hImpactParBachFromCharmBaryonXY", "hImpactParBachFromCharmBaryonXY", kTH1D, {{1000, -1.0f, 1.0f}});
    registry.add("hImpactParBachFromCharmBaryonZ", "hImpactParBachFromCharmBaryonZ", kTH1D, {{1000, -2.0f, 2.0f}});
    registry.add("hImpactParCascXY", "hImpactParCascXY", kTH1D, {{1000, -4.0f, 4.0f}});
    registry.add("hImpactParCascZ", "hImpactParCascZ", kTH1D, {{1000, -5.0f, 5.0f}});
    registry.add("hPtKaFromCasc", "hPtKaFromCasc", kTH1D, {{1000, 0.0f, 5.0f}});
    registry.add("hPtPiFromCharmBaryon", "hPtPiFromCharmBaryon", kTH1D, {{1000, 0.0f, 5.0f}});
    registry.add("hCTauOmegac", "hCTauOmegac", kTH1D, {{1000, 0.0f, 0.1f}});
    registry.add("hKFGeoV0Chi2OverNdf", "hKFGeoV0Chi2OverNdf", kTH1D, {{1000, 0.0f, 100.0f}});
    registry.add("hKFGeoCascChi2OverNdf", "hKFGeoCascChi2OverNdf", kTH1D, {{1000, 0.0f, 100.0f}});
    registry.add("hKFGeoCharmbaryonChi2OverNdf", "hKFGeoCharmbaryonChi2OverNdf", kTH1D, {{1000, 0.0f, 100.0f}});
    registry.add("hKFdecayLenXYLambda", "hKFdecayLenXYLambda", kTH1D, {{1000, 0.0f, 50.0f}});
    registry.add("hKFdecayLenXYCasc", "hKFdecayLenXYCasc", kTH1D, {{1000, 0.0f, 50.0f}});
    registry.add("hKFdecayLenXYOmegac", "hKFdecayLenXYOmegac", kTH1D, {{1000, 0.0f, 50.0f}});
    registry.add("hKFcosPaV0ToCasc", "hKFcosPaV0ToCasc", kTH1D, {{5000, 0.8f, 1.1f}});
    registry.add("hKFcosPaCascToOmegac", "hKFcosPaCascToOmegac", kTH1D, {{5000, 0.8f, 1.1f}});

    hfEvSel.addHistograms(registry); // collision monitoring

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

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, int decayChannel, typename Coll, typename Hist, typename TCascTable, typename TCascLinkTable>
  void runXic0Omegac0Creator(Coll const&,
                             aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                             MyLFTracksWCov const& lfTracks,
                             TracksWCovDca const& tracks,
                             TCascTable const&, TCascLinkTable const&,
                             aod::HfCascLf2Prongs const& candidates,
                             Hist& hInvMassCharmBaryon,
                             Hist& hFitterStatus,
                             Hist& hCandidateCounter,
                             Hist& hCascadesCounter)
  {

    if constexpr (decayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi && decayChannel != hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi && decayChannel != hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
      LOGP(fatal, "Decay channel not recognized!");
    }

    for (const auto& cand : candidates) {

      hCandidateCounter->Fill(0);

      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
          continue;
        }
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
          continue;
        }
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK)) {
          continue;
        }
      }

      hCandidateCounter->Fill(1);

      auto collision = cand.collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      // set the magnetic field from CCDB
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << magneticField;
        runNumber = bc.runNumber();
      }
      df.setBz(magneticField);

      auto trackCharmBachelorId = cand.prong0Id();
      auto trackCharmBachelor = tracks.rawIteratorAt(trackCharmBachelorId);

      auto cascAodElement = cand.template cascade_as<TCascLinkTable>();
      hCascadesCounter->Fill(0);
      int v0index = cascAodElement.v0Id();

      // check if the cascade from AO2D has data
      bool hasData = false;
      if constexpr (requires { cascAodElement.cascDataId(); }) { // check if it's the CascDataLink
        if (cascAodElement.has_cascData()) {
          hasData = true;
        }
      }
      if constexpr (requires { cascAodElement.traCascDataId(); }) { // check if it's the TraCascDataLink
        if (cascAodElement.has_traCascData()) {
          hasData = true;
        }
      }
      if (!hasData) {
        continue;
      }

      typename TCascTable::iterator casc;
      if constexpr (requires { cascAodElement.cascDataId(); }) { // check if it's the CascDataLink
        casc = cascAodElement.template cascData_as<TCascTable>();
      }
      if constexpr (requires { cascAodElement.traCascDataId(); }) { // check if it's the TraCascDataLink
        casc = cascAodElement.template traCascData_as<TCascTable>();
      }

      hCascadesCounter->Fill(1);
      auto trackCascDauChargedId = casc.bachelorId();                           // pion <- xi track
      auto trackV0Dau0Id = casc.posTrackId();                                   // V0 positive daughter track
      auto trackV0Dau1Id = casc.negTrackId();                                   // V0 negative daughter track
      auto trackCascDauCharged = lfTracks.rawIteratorAt(trackCascDauChargedId); // pion <- xi track
      auto trackV0Dau0 = lfTracks.rawIteratorAt(trackV0Dau0Id);                 // V0 positive daughter track
      auto trackV0Dau1 = lfTracks.rawIteratorAt(trackV0Dau1Id);                 // V0 negative daughter track

      //-------------------------- V0 info---------------------------
      // pseudorapidity
      float pseudorapV0Dau0 = casc.positiveeta();
      float pseudorapV0Dau1 = casc.negativeeta();

      // info from LF table
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
      std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      float pseudorapCascBachelor = casc.bacheloreta();

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
      if (trackCascDauCharged.sign() > 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else if (trackCascDauCharged.sign() < 0) {
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else {
        continue;
      }
      trackCasc.setAbsCharge(1);
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        trackCasc.setPID(o2::track::PID::XiMinus);
      } else {
        trackCasc.setPID(o2::track::PID::OmegaMinus);
      }

      std::array<float, 3> pVecCascBachelor = {casc.pxbach(), casc.pybach(), casc.pzbach()};

      //------------reconstruct charm baryon decay vtx---------------
      auto trackParVarCharmBachelor = getTrackParCov(trackCharmBachelor); // charm bachelor pion track to be processed with DCAFitter

      // reconstruct charm baryon with DCAFitter
      int nVtxFromFitterCharmBaryon = 0;
      try {
        nVtxFromFitterCharmBaryon = df.process(trackCasc, trackParVarCharmBachelor);
      } catch (...) {
        LOG(error) << "Exception caught in charm DCA fitter process call!";
        hFitterStatus->Fill(1);
        continue;
      }
      if (nVtxFromFitterCharmBaryon == 0) {
        hFitterStatus->Fill(2);
        continue;
      }
      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);
      auto vertexCharmBaryonFromFitter = df.getPCACandidate();
      std::array<float, 3> pVecCascAsD;
      std::array<float, 3> pVecCharmBachelorAsD;
      df.propagateTracksToVertex();
      if (!df.isPropagateTracksToVertexDone()) {
        continue;
      }
      df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
      df.getTrack(1).getPxPyPzGlo(pVecCharmBachelorAsD);
      std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecCharmBachelorAsD[0], pVecCascAsD[1] + pVecCharmBachelorAsD[1], pVecCascAsD[2] + pVecCharmBachelorAsD[2]};

      std::array<float, 3> coordVtxCharmBaryon = df.getPCACandidatePos();
      std::array<float, 6> covVtxCharmBaryon = df.calcPCACovMatrixFlat();

      // pseudorapidity
      float pseudorapCharmBachelor = trackCharmBachelor.eta();

      // primary vertex of the collision
      auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

      // DCAxy and DCAz (computed with propagateToDCABxByBz method)
      o2::dataformats::DCA impactParameterV0Dau0;
      o2::dataformats::DCA impactParameterV0Dau1;
      o2::dataformats::DCA impactParameterCascDauCharged;
      auto trackParVarV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParVarV0Dau1 = getTrackParCov(trackV0Dau1);
      auto trackParVarCascDauCharged = getTrackParCov(trackCascDauCharged);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCascDauCharged, 2.f, matCorr, &impactParameterCascDauCharged);
      float dcaxyV0Dau0 = impactParameterV0Dau0.getY();
      float dcaxyV0Dau1 = impactParameterV0Dau1.getY();
      float dcaxyCascBachelor = impactParameterCascDauCharged.getY();
      float dcazV0Dau0 = impactParameterV0Dau0.getZ();
      float dcazV0Dau1 = impactParameterV0Dau1.getZ();
      float dcazCascBachelor = impactParameterCascDauCharged.getZ();

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterCharmBachelor;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);
      float impactParBachFromCharmBaryonXY = impactParameterCharmBachelor.getY();
      float impactParBachFromCharmBaryonZ = impactParameterCharmBachelor.getZ();

      // invariant mass under the hypothesis of particles ID corresponding to the decay chain
      float mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
      float mCasc = 0.;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        mCasc = casc.mXi();
      } else {
        mCasc = casc.mOmega();
      }
      auto arrMassCharmBaryon = std::array{0., 0.};
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        arrMassCharmBaryon = {MassXiMinus, MassPiPlus};
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        arrMassCharmBaryon = {MassOmegaMinus, MassPiPlus};
      } else {
        arrMassCharmBaryon = {MassOmegaMinus, MassKPlus};
      }
      float mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecCharmBachelorAsD}, arrMassCharmBaryon);

      // computing cosPA
      float cpaV0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      float cpaCasc = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float cpaxyV0 = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      float cpaxyCasc = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);

      // computing decay length and ctau
      float decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      float decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      float decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
      auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
      auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

      float ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassOmegaC0);
      float ctXic = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassXiC0);
      float ctCascade = 0.;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassXiMinus);
      } else {
        ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassOmegaMinus);
      }
      float ctV0 = RecoDecay::ct(pVecV0, decLenV0, MassLambda0);

      // computing eta
      float pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
      float pseudorapCascade = RecoDecay::eta(pVecCasc);
      float pseudorapV0 = RecoDecay::eta(pVecV0);

      // DCA between daughters
      float dcaCascDau = casc.dcacascdaughters();
      float dcaV0Dau = casc.dcaV0daughters();
      float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

      // fill test histograms
      hInvMassCharmBaryon->Fill(mCharmBaryon);
      hCandidateCounter->Fill(3);

      // fill the table
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        rowCandToXiPi(collision.globalIndex(),
                      pvCoord[0], pvCoord[1], pvCoord[2],
                      vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                      vertexCasc[0], vertexCasc[1], vertexCasc[2],
                      vertexV0[0], vertexV0[1], vertexV0[2],
                      trackCascDauCharged.sign(),
                      covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                      pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                      pVecCascAsD[0], pVecCascAsD[1], pVecCascAsD[2],
                      pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                      pVecV0[0], pVecV0[1], pVecV0[2],
                      pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
                      pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                      pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                      impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
                      impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
                      std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                      v0index, casc.posTrackId(), casc.negTrackId(),
                      casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                      mLambda, mCasc, mCharmBaryon,
                      cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                      ctOmegac, ctCascade, ctV0, ctXic,
                      pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                      pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                      dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                      dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
                      dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                      decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);

      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        rowCandToOmegaPi(collision.globalIndex(),
                         pvCoord[0], pvCoord[1], pvCoord[2],
                         vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                         vertexCasc[0], vertexCasc[1], vertexCasc[2],
                         vertexV0[0], vertexV0[1], vertexV0[2],
                         trackCascDauCharged.sign(),
                         covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                         pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                         pVecCascAsD[0], pVecCascAsD[1], pVecCascAsD[2],
                         pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                         pVecV0[0], pVecV0[1], pVecV0[2],
                         pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
                         pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                         pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                         impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
                         impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
                         std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                         v0index, casc.posTrackId(), casc.negTrackId(),
                         casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                         mLambda, mCasc, mCharmBaryon,
                         cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                         ctOmegac, ctCascade, ctV0,
                         pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                         pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                         dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                         dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
                         dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                         decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);

      } else {
        rowCandToOmegaK(
          collision.globalIndex(), pvCoord[0], pvCoord[1], pvCoord[2],
          vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
          vertexCasc[0], vertexCasc[1], vertexCasc[2],
          vertexV0[0], vertexV0[1], vertexV0[2],
          trackCascDauCharged.sign(),
          covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
          pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
          pVecCascAsD[0], pVecCascAsD[1], pVecCascAsD[2],
          pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
          pVecV0[0], pVecV0[1], pVecV0[2],
          pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
          pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
          pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
          impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
          impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
          std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
          v0index, casc.posTrackId(), casc.negTrackId(),
          casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
          mLambda, mCasc, mCharmBaryon,
          cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
          ctOmegac, ctCascade, ctV0,
          pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
          pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
          dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
          dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
          dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
          decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);
      }

    } // loop over LF Cascade-bachelor candidates
  } // end of run function

  template <int decayChannel, typename Coll, typename Hist>
  void runKfOmegac0CreatorWithKFParticle(Coll const&,
                                         aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                                         MyKfTracks const&,
                                         MyKfCascTable const&, KFCascadesLinked const&,
                                         aod::HfCascLf2Prongs const& candidates,
                                         Hist& hInvMassCharmBaryon,
                                         Hist& hFitterStatus,
                                         Hist& hCandidateCounter,
                                         Hist& hCascadesCounter)
  {
    for (const auto& cand : candidates) {
      hCandidateCounter->Fill(1);

      auto collision = cand.collision_as<Coll>();

      // set the magnetic field from CCDB
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << magneticField;
        runNumber = bc.runNumber();
      }
      df.setBz(magneticField);
      KFParticle::SetField(magneticField);
      // bachelor from Omegac0
      auto trackCharmBachelor = cand.prong0_as<MyKfTracks>();

      auto cascAodElement = cand.cascade_as<aod::KFCascadesLinked>();
      hCascadesCounter->Fill(0);
      int v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<MyKfCascTable>();
      hCascadesCounter->Fill(1);
      auto trackCascDauCharged = casc.bachelor_as<MyKfTracks>(); // pion <- xi track
      auto trackV0Dau0 = casc.posTrack_as<MyKfTracks>();         // V0 positive daughter track
      auto trackV0Dau1 = casc.negTrack_as<MyKfTracks>();         // V0 negative daughter track

      auto bachCharge = trackCascDauCharged.signed1Pt() > 0 ? +1 : -1;

      //// pion & p TrackParCov
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);
      // kaon <- casc TrackParCov
      auto omegaDauChargedTrackParCov = getTrackParCov(trackCascDauCharged);
      // convert tracks into KFParticle object
      KFPTrack kfTrack0 = createKFPTrackFromTrack(trackV0Dau0);
      KFPTrack kfTrack1 = createKFPTrackFromTrack(trackV0Dau1);
      KFPTrack kfTrackBach = createKFPTrackFromTrack(trackCascDauCharged);

      KFParticle kfPosPr(kfTrack0, kProton);
      KFParticle kfNegPi(kfTrack1, kPiMinus);
      KFParticle kfNegKa(kfTrackBach, kKMinus);
      KFParticle kfNegPiRej(kfTrackBach, kPiMinus); // rej
      KFParticle kfPosPi(kfTrack0, kPiPlus);
      KFParticle kfNegPr(kfTrack1, kProton);
      KFParticle kfPosKa(kfTrackBach, kKPlus);
      KFParticle kfPosPiRej(kfTrackBach, kPiPlus); // rej

      KFParticle kfBachKaon;
      KFParticle kfPos;
      KFParticle kfNeg;
      KFParticle kfBachPionRej; // rej
      if (bachCharge < 0) {
        kfPos = kfPosPr;
        kfNeg = kfNegPi;
        kfBachKaon = kfNegKa;
        kfBachPionRej = kfNegPiRej; // rej
      } else {
        kfPos = kfPosPi;
        kfNeg = kfNegPr;
        kfBachKaon = kfPosKa;
        kfBachPionRej = kfPosPiRej; // rej
      }

      //__________________________________________
      //*>~<* step 1 : construct V0 with KF
      const KFParticle* v0Daughters[2] = {&kfPos, &kfNeg};
      // construct V0
      KFParticle kfV0;
      kfV0.SetConstructMethod(kfConstructMethod);
      try {
        kfV0.Construct(v0Daughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
        continue;
      }

      // mass window cut on lambda before mass constraint
      float massLam, sigLam;
      kfV0.GetMass(massLam, sigLam);
      if (TMath::Abs(massLam - MassLambda0) > lambdaMassWindow)
        continue;
      // err_mass>0 of Lambda
      if (sigLam <= 0)
        continue;
      // chi2>0 && NDF>0 for selecting Lambda
      if ((kfV0.GetNDF() <= 0 || kfV0.GetChi2() <= 0))
        continue;
      kfOmegac0Candidate.chi2GeoV0 = kfV0.GetChi2();
      KFParticle kfV0MassConstrained = kfV0;
      kfV0MassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassLambda); // set mass constrain to Lambda
      if (kfUseV0MassConstraint) {
        KFParticle kfV0 = kfV0MassConstrained;
      }
      kfV0.TransportToDecayVertex();

      //__________________________________________
      //*>~<* step 2 : reconstruct cascade(Omega) with KF
      const KFParticle* omegaDaugthers[2] = {&kfBachKaon, &kfV0};
      const KFParticle* omegaDaugthersRej[2] = {&kfBachPionRej, &kfV0}; // rej
      // construct cascade
      KFParticle kfOmega;
      KFParticle kfOmegarej; // rej
      kfOmega.SetConstructMethod(kfConstructMethod);
      kfOmegarej.SetConstructMethod(kfConstructMethod); // rej
      try {
        kfOmega.Construct(omegaDaugthers, 2);
        kfOmegarej.Construct(omegaDaugthersRej, 2); // rej
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Omega or Omega_rej from V0 and bachelor track: " << e.what();
        continue;
      }
      float massCasc, sigCasc;
      float massCascrej, sigCascrej;
      kfOmega.GetMass(massCasc, sigCasc);
      kfOmegarej.GetMass(massCascrej, sigCascrej); // rej
      // err_massOmega > 0
      if (sigCasc <= 0)
        continue;
      if (std::abs(massCasc - MassOmegaMinus) > massToleranceCascade)
        continue;
      // chi2>0 && NDF>0
      if (kfOmega.GetNDF() <= 0 || kfOmega.GetChi2() <= 0)
        continue;
      kfOmegac0Candidate.chi2GeoCasc = kfOmega.GetChi2();
      kfOmegac0Candidate.cascRejectInvmass = massCascrej;
      registry.fill(HIST("hInvMassXiMinus_rej"), massCascrej); // rej
      KFParticle kfOmegaMassConstrained = kfOmega;
      kfOmegaMassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassOmegaMinus); // set mass constrain to OmegaMinus
      if (kfUseCascadeMassConstraint) {
        // set mass constraint if requested
        KFParticle kfOmega = kfOmegaMassConstrained;
      }
      registry.fill(HIST("hInvMassOmegaMinus"), massCasc);
      kfOmega.TransportToDecayVertex();
      // rej: Add competing rejection to minimize misidentified Xi impact. Reject if kfBachPionRej is Pion and the constructed cascade has Xi's invariant mass.

      //__________________________________________
      //*>~<* step 3 : reconstruc Omegac0 with KF
      // Create KF charm bach Pion from track
      KFPTrack kfTrackBachPion = createKFPTrackFromTrack(trackCharmBachelor);
      KFParticle kfBachPion(kfTrackBachPion, kPiPlus);
      const KFParticle* omegaC0Daugthers[2] = {&kfBachPion, &kfOmega};

      // construct OmegaC0
      KFParticle kfOmegaC0;
      kfOmegaC0.SetConstructMethod(kfConstructMethod);
      try {
        kfOmegaC0.Construct(omegaC0Daugthers, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct OmegaC0 from Cascade and bachelor pion track: " << e.what();
        continue;
      }
      float massOmegaC0, sigOmegaC0;
      kfOmegaC0.GetMass(massOmegaC0, sigOmegaC0);
      if (sigOmegaC0 <= 0)
        continue;
      // chi2>0 && NDF>0
      if (kfOmegaC0.GetNDF() <= 0 || kfOmegaC0.GetChi2() <= 0)
        continue;
      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);
      kfOmegaC0.TransportToDecayVertex();
      // PV
      KFPVertex kfVertex = createKFPVertexFromCollision(collision);
      KFParticle kfPV(kfVertex);

      // set production vertex;
      kfNeg.SetProductionVertex(kfV0);
      kfPos.SetProductionVertex(kfV0);

      KFParticle kfBachKaonToOmega = kfBachKaon;
      KFParticle kfV0ToCasc = kfV0;
      kfBachKaonToOmega.SetProductionVertex(kfOmega);
      kfV0ToCasc.SetProductionVertex(kfOmega);

      KFParticle kfOmegaToOmegaC = kfOmega;
      KFParticle kfBachPionToOmegaC = kfBachPion;
      kfBachPionToOmegaC.SetProductionVertex(kfOmegaC0);
      kfOmegaToOmegaC.SetProductionVertex(kfOmegaC0);

      // KFParticle to PV
      KFParticle kfV0ToPv = kfV0;
      KFParticle kfOmegaToPv = kfOmega;
      KFParticle kfOmegac0ToPv = kfOmegaC0;
      KFParticle kfPiFromOmegacToPv = kfBachPion;

      kfV0ToPv.SetProductionVertex(kfPV);
      kfOmegaToPv.SetProductionVertex(kfPV);
      kfOmegac0ToPv.SetProductionVertex(kfPV);
      kfPiFromOmegacToPv.SetProductionVertex(kfPV);
      //------------get updated daughter tracks after vertex fit  ---------------
      auto trackParVarCharmBachelor = getTrackParCovFromKFP(kfBachPionToOmegaC, o2::track::PID::Pion, -bachCharge); // chrambaryon bach pion
      trackParVarCharmBachelor.setAbsCharge(1);

      omegaDauChargedTrackParCov = getTrackParCovFromKFP(kfBachKaonToOmega, o2::track::PID::Kaon, bachCharge); // Cascade bach kaon
      omegaDauChargedTrackParCov.setAbsCharge(1);
      o2::track::TrackParCov trackCasc = getTrackParCovFromKFP(kfOmegaToOmegaC, kfOmegaToOmegaC.GetPDG(), bachCharge);
      trackCasc.setAbsCharge(1);

      trackParCovV0Dau0 = getTrackParCovFromKFP(kfPos, kfPos.GetPDG(), 1); // V0 postive daughter
      trackParCovV0Dau0.setAbsCharge(1);
      trackParCovV0Dau1 = getTrackParCovFromKFP(kfNeg, kfNeg.GetPDG(), -1); // V0 negtive daughter
      trackParCovV0Dau1.setAbsCharge(1);

      //-------------------------- V0 info---------------------------
      // pseudorapidity
      float pseudorapV0Dau0 = kfPos.GetEta();
      float pseudorapV0Dau1 = kfNeg.GetEta();

      // info from from KFParticle
      std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      std::array<float, 3> pVecV0Dau0 = {kfPos.GetPx(), kfPos.GetPy(), kfPos.GetPz()};
      std::array<float, 3> pVecV0Dau1 = {kfNeg.GetPx(), kfNeg.GetPy(), kfNeg.GetPz()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      float pseudorapCascBachelor = kfBachKaonToOmega.GetEta();

      // info from KFParticle
      std::array<float, 3> vertexCasc = {kfOmega.GetX(), kfOmega.GetY(), kfOmega.GetZ()};
      std::array<float, 3> pVecCascBachelor = {kfBachKaonToOmega.GetPx(), kfBachKaonToOmega.GetPy(), kfBachKaonToOmega.GetPz()};

      auto primaryVertex = getPrimaryVertex(collision);
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> vertexCharmBaryonFromFitter = {0.0, 0.0, 0.0}; // This variable get from DCAfitter in default process, in KF process it is set as 0.
      std::array<float, 3> pVecCharmBachelorAsD;
      pVecCharmBachelorAsD[0] = kfBachPionToOmegaC.GetPx();
      pVecCharmBachelorAsD[1] = kfBachPionToOmegaC.GetPy();
      pVecCharmBachelorAsD[2] = kfBachPionToOmegaC.GetPz();

      std::array<float, 3> pVecCharmBaryon = {kfOmegaC0.GetPx(), kfOmegaC0.GetPy(), kfOmegaC0.GetPz()};
      std::array<float, 3> coordVtxCharmBaryon = {kfOmegaC0.GetX(), kfOmegaC0.GetY(), kfOmegaC0.GetZ()};
      auto covVtxCharmBaryon = kfOmegaC0.CovarianceMatrix();
      float covMatrixPV[6];
      kfVertex.GetCovarianceMatrix(covMatrixPV);

      // impact parameters
      gpu::gpustd::array<float, 2> impactParameterV0Dau0;
      gpu::gpustd::array<float, 2> impactParameterV0Dau1;
      gpu::gpustd::array<float, 2> impactParameterKaFromCasc;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, omegaDauChargedTrackParCov, 2.f, matCorr, &impactParameterKaFromCasc);
      float dcaxyV0Dau0 = impactParameterV0Dau0[0];
      float dcaxyV0Dau1 = impactParameterV0Dau1[0];
      float dcaxyCascBachelor = impactParameterKaFromCasc[0];
      float dcazV0Dau0 = impactParameterV0Dau0[1];
      float dcazV0Dau1 = impactParameterV0Dau1[1];
      float dcazCascBachelor = impactParameterKaFromCasc[1];

      // pseudorapidity
      float pseudorapCharmBachelor = kfBachPionToOmegaC.GetEta();

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterCharmBachelor;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);
      float impactParBachFromCharmBaryonXY = impactParameterCharmBachelor.getY();
      float impactParBachFromCharmBaryonZ = impactParameterCharmBachelor.getZ();

      // computing decay length and ctau
      float decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      float decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      float decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array{kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
      auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
      auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

      // fill test histograms
      hInvMassCharmBaryon->Fill(massOmegaC0);
      hCandidateCounter->Fill(3);

      //// KFParticle table information
      // KF chi2
      auto v0NDF = kfV0.GetNDF();
      auto v0Chi2OverNdf = kfOmegac0Candidate.chi2GeoV0 / v0NDF;

      auto cascNDF = kfOmega.GetNDF();
      auto cascChi2OverNdf = kfOmegac0Candidate.chi2GeoCasc / cascNDF;

      kfOmegac0Candidate.chi2GeoOmegac = kfOmegaC0.GetChi2();
      auto charmbaryonNDF = kfOmegaC0.GetNDF();
      auto charmbaryonChi2OverNdf = kfOmegac0Candidate.chi2GeoOmegac / charmbaryonNDF;

      kfOmegac0Candidate.chi2MassV0 = kfV0MassConstrained.GetChi2();
      auto v0Ndfm = kfV0MassConstrained.GetNDF();
      auto v0Chi2OverNdfm = kfOmegac0Candidate.chi2MassV0 / v0Ndfm;

      kfOmegac0Candidate.chi2MassCasc = kfOmegaMassConstrained.GetChi2();
      auto cascNdfm = kfOmegaMassConstrained.GetNDF();
      auto cascChi2OverNdfm = kfOmegac0Candidate.chi2MassCasc / cascNdfm;

      // KF topo Chi2 over NDF
      kfOmegac0Candidate.chi2NdfTopoV0ToPv = kfV0ToPv.GetChi2() / kfV0ToPv.GetNDF();
      kfOmegac0Candidate.chi2NdfTopoCascToPv = kfOmegaToPv.GetChi2() / kfOmegaToPv.GetNDF();
      kfOmegac0Candidate.chi2NdfTopoPiFromOmegacToPv = kfPiFromOmegacToPv.GetChi2() / kfPiFromOmegacToPv.GetNDF();
      kfOmegac0Candidate.chi2NdfTopoOmegacToPv = kfOmegac0ToPv.GetChi2() / kfOmegac0ToPv.GetNDF();

      auto cascBachTopoChi2Ndf = kfBachKaonToOmega.GetChi2() / kfBachKaonToOmega.GetNDF();
      kfOmegac0Candidate.chi2NdfTopoV0ToCasc = kfV0ToCasc.GetChi2() / kfV0ToCasc.GetNDF();
      kfOmegac0Candidate.chi2NdfTopoCascToOmegac = kfOmegaToOmegaC.GetChi2() / kfOmegaToOmegaC.GetNDF();

      // KF ldl
      kfOmegac0Candidate.ldlV0 = ldlFromKF(kfV0, kfPV);
      kfOmegac0Candidate.ldlCasc = ldlFromKF(kfOmega, kfPV);
      kfOmegac0Candidate.ldlOmegac = ldlFromKF(kfOmegaC0, kfPV);

      // KF dca
      kfOmegac0Candidate.kfDcaXYPiFromOmegac = kfBachPionToOmegaC.GetDistanceFromVertexXY(kfPV);
      kfOmegac0Candidate.kfDcaV0Dau = kfNeg.GetDistanceFromParticle(kfPos);
      kfOmegac0Candidate.kfDcaCascDau = kfBachKaonToOmega.GetDistanceFromParticle(kfV0ToCasc);
      kfOmegac0Candidate.kfDcaXYCascToPv = kfOmegaToOmegaC.GetDistanceFromVertexXY(kfPV);
      kfOmegac0Candidate.kfDcaOmegacDau = kfBachPionToOmegaC.GetDistanceFromParticle(kfOmegaToOmegaC);

      // KF decay length
      float decayLxyLam, errDecayLxyLam;
      kfV0ToCasc.GetDecayLengthXY(decayLxyLam, errDecayLxyLam);
      kfOmegac0Candidate.decayLenXYLambda = decayLxyLam;

      float decayLxyCasc, errDecayLxyCasc;
      kfOmegaToOmegaC.GetDecayLengthXY(decayLxyCasc, errDecayLxyCasc);
      kfOmegac0Candidate.decayLenXYCasc = decayLxyCasc;

      float decayLxyOmegac0, errDecayLxyOmegac0;
      kfOmegac0ToPv.GetDecayLengthXY(decayLxyOmegac0, errDecayLxyOmegac0);
      kfOmegac0Candidate.decayLenXYOmegac = decayLxyOmegac0;

      // KF cosPA
      kfOmegac0Candidate.cosPaV0ToPv = cpaFromKF(kfV0, kfPV);
      kfOmegac0Candidate.cosPaCascToPv = cpaFromKF(kfOmega, kfPV);
      kfOmegac0Candidate.cosPaOmegacToPv = cpaFromKF(kfOmegaC0, kfPV);
      kfOmegac0Candidate.cosPaXYV0ToPv = cpaXYFromKF(kfV0, kfPV);
      kfOmegac0Candidate.cosPaXYCascToPv = cpaXYFromKF(kfOmega, kfPV);
      kfOmegac0Candidate.cosPaXYOmegacToPv = cpaXYFromKF(kfOmegaC0, kfPV);

      kfOmegac0Candidate.cosPaV0ToCasc = cpaFromKF(kfV0, kfOmega);
      kfOmegac0Candidate.cosPaCascToOmegac = cpaFromKF(kfOmega, kfOmegaC0);
      kfOmegac0Candidate.cosPaXYV0ToCasc = cpaXYFromKF(kfV0, kfOmega);
      kfOmegac0Candidate.cosPaXYCascToOmegac = cpaXYFromKF(kfOmega, kfOmegaC0);
      // KF mass
      kfOmegac0Candidate.massV0 = massLam;
      kfOmegac0Candidate.massCasc = massCasc;
      kfOmegac0Candidate.massOmegac = massOmegaC0;

      // KF pT
      kfOmegac0Candidate.ptPiFromOmegac = kfBachPionToOmegaC.GetPt();
      kfOmegac0Candidate.ptOmegac = kfOmegaC0.GetPt();

      // KF rapidity
      kfOmegac0Candidate.rapOmegac = kfOmegaC0.GetRapidity();

      // KF cosThetaStar
      kfOmegac0Candidate.cosThetaStarPiFromOmegac = cosThetaStarFromKF(0, 4332, 211, 3312, kfBachPionToOmegaC, kfOmegaToOmegaC);

      // KF ct
      kfOmegac0Candidate.ctV0 = kfV0.GetLifeTime();
      kfOmegac0Candidate.ctCasc = kfOmega.GetLifeTime();
      kfOmegac0Candidate.ctOmegac = kfOmegaC0.GetLifeTime();

      // KF eta
      kfOmegac0Candidate.etaOmegac = kfOmegaC0.GetEta();

      // fill KF hist
      registry.fill(HIST("hKFParticleCascBachTopoChi2"), cascBachTopoChi2Ndf);
      registry.fill(HIST("hKFParticleV0TopoChi2"), kfOmegac0Candidate.chi2NdfTopoV0ToCasc);
      registry.fill(HIST("hKFParticleCascTopoChi2"), kfOmegac0Candidate.chi2NdfTopoCascToOmegac);
      registry.fill(HIST("hKFParticlechi2TopoOmegacToPv"), kfOmegac0Candidate.chi2NdfTopoOmegacToPv);
      registry.fill(HIST("hKFParticlechi2TopoCascToPv"), kfOmegac0Candidate.chi2NdfTopoCascToPv);
      registry.fill(HIST("hKFParticleDcaCharmBaryonDau"), kfOmegac0Candidate.kfDcaOmegacDau);
      registry.fill(HIST("hKFParticleDcaXYCascBachToPv"), dcaxyCascBachelor);
      registry.fill(HIST("hKFParticleDcaXYV0DauPosToPv"), dcaxyV0Dau0);
      registry.fill(HIST("hKFParticleDcaXYV0DauNegToPv"), dcaxyV0Dau1);
      registry.fill(HIST("hKfLambda_ldl"), kfOmegac0Candidate.ldlV0);
      registry.fill(HIST("hKfOmega_ldl"), kfOmegac0Candidate.ldlCasc);
      registry.fill(HIST("hKfOmegaC0_ldl"), kfOmegac0Candidate.ldlOmegac);
      registry.fill(HIST("hDcaXYCascadeToPVKf"), kfOmegac0Candidate.kfDcaXYCascToPv);
      // Additional histograms
      if (fillAllHist) {
        registry.fill(HIST("hEtaV0PosDau"), kfPos.GetEta());
        registry.fill(HIST("hEtaV0NegDau"), kfNeg.GetEta());
        registry.fill(HIST("hEtaKaFromCasc"), kfBachKaonToOmega.GetEta());
        registry.fill(HIST("hEtaPiFromCharmBaryon"), kfBachPionToOmegaC.GetEta());
        registry.fill(HIST("hCascradius"), RecoDecay::sqrtSumOfSquares(vertexCasc[0], vertexCasc[1]));
        registry.fill(HIST("hV0radius"), RecoDecay::sqrtSumOfSquares(vertexV0[0], vertexV0[1]));
        registry.fill(HIST("hCosPACasc"), kfOmegac0Candidate.cosPaCascToPv);
        registry.fill(HIST("hCosPAV0"), kfOmegac0Candidate.cosPaV0ToPv);
        registry.fill(HIST("hDcaCascDau"), kfOmegac0Candidate.kfDcaCascDau);
        registry.fill(HIST("hDcaV0Dau"), kfOmegac0Candidate.kfDcaV0Dau);
        registry.fill(HIST("hDcaXYToPvKa"), dcaxyCascBachelor);
        registry.fill(HIST("hImpactParBachFromCharmBaryonXY"), impactParBachFromCharmBaryonXY);
        registry.fill(HIST("hImpactParBachFromCharmBaryonZ"), impactParBachFromCharmBaryonZ);
        registry.fill(HIST("hImpactParCascXY"), impactParameterCasc.getY());
        registry.fill(HIST("hImpactParCascZ"), impactParameterCasc.getZ());
        registry.fill(HIST("hPtKaFromCasc"), RecoDecay::sqrtSumOfSquares(pVecCascBachelor[0], pVecCascBachelor[1]));
        registry.fill(HIST("hPtPiFromCharmBaryon"), RecoDecay::sqrtSumOfSquares(pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1]));
        registry.fill(HIST("hCTauOmegac"), kfOmegac0Candidate.ctOmegac);
        registry.fill(HIST("hKFGeoV0Chi2OverNdf"), v0Chi2OverNdf);
        registry.fill(HIST("hKFGeoCascChi2OverNdf"), cascChi2OverNdf);
        registry.fill(HIST("hKFGeoCharmbaryonChi2OverNdf"), charmbaryonChi2OverNdf);
        registry.fill(HIST("hKFdecayLenXYLambda"), kfOmegac0Candidate.decayLenXYLambda);
        registry.fill(HIST("hKFdecayLenXYCasc"), kfOmegac0Candidate.decayLenXYCasc);
        registry.fill(HIST("hKFdecayLenXYOmegac"), kfOmegac0Candidate.decayLenXYOmegac);
        registry.fill(HIST("hKFcosPaV0ToCasc"), kfOmegac0Candidate.cosPaV0ToCasc);
        registry.fill(HIST("hKFcosPaCascToOmegac"), kfOmegac0Candidate.cosPaCascToOmegac);
      }

      // fill the table
      rowCandToOmegaPi(collision.globalIndex(),
                       pvCoord[0], pvCoord[1], pvCoord[2],
                       vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                       vertexCasc[0], vertexCasc[1], vertexCasc[2],
                       vertexV0[0], vertexV0[1], vertexV0[2],
                       trackCascDauCharged.sign(),
                       covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                       pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                       kfOmegaToOmegaC.GetPx(), kfOmegaToOmegaC.GetPy(), kfOmegaToOmegaC.GetPz(),
                       pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                       pVecV0[0], pVecV0[1], pVecV0[2],
                       pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
                       pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                       pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                       impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
                       impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
                       std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                       v0index, casc.posTrackId(), casc.negTrackId(),
                       casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                       kfOmegac0Candidate.massV0, kfOmegac0Candidate.massCasc, kfOmegac0Candidate.massOmegac,
                       kfOmegac0Candidate.cosPaV0ToPv, kfOmegac0Candidate.cosPaOmegacToPv, kfOmegac0Candidate.cosPaCascToPv, kfOmegac0Candidate.cosPaXYV0ToPv, kfOmegac0Candidate.cosPaXYOmegacToPv, kfOmegac0Candidate.cosPaXYCascToPv,
                       kfOmegac0Candidate.ctOmegac, kfOmegac0Candidate.ctCasc, kfOmegac0Candidate.ctV0,
                       pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                       kfOmegac0Candidate.etaOmegac, kfOmega.GetEta(), kfV0.GetEta(),
                       dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                       dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
                       kfOmegac0Candidate.kfDcaCascDau, kfOmegac0Candidate.kfDcaV0Dau, kfOmegac0Candidate.kfDcaOmegacDau,
                       decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon);
      // fill kf table
      kfCandidateData(kfOmegac0Candidate.kfDcaXYPiFromOmegac, kfOmegac0Candidate.kfDcaXYCascToPv,
                      kfOmegac0Candidate.chi2GeoV0, kfOmegac0Candidate.chi2GeoCasc, kfOmegac0Candidate.chi2GeoOmegac, kfOmegac0Candidate.chi2MassV0, kfOmegac0Candidate.chi2MassCasc,
                      kfOmegac0Candidate.ldlV0, kfOmegac0Candidate.ldlCasc, kfOmegac0Candidate.ldlOmegac,
                      kfOmegac0Candidate.chi2NdfTopoV0ToPv, kfOmegac0Candidate.chi2NdfTopoCascToPv, kfOmegac0Candidate.chi2NdfTopoPiFromOmegacToPv, kfOmegac0Candidate.chi2NdfTopoOmegacToPv,
                      kfOmegac0Candidate.chi2NdfTopoV0ToCasc, kfOmegac0Candidate.chi2NdfTopoCascToOmegac,
                      kfOmegac0Candidate.decayLenXYLambda, kfOmegac0Candidate.decayLenXYCasc, kfOmegac0Candidate.decayLenXYOmegac,
                      kfOmegac0Candidate.cosPaV0ToCasc, kfOmegac0Candidate.cosPaCascToOmegac, kfOmegac0Candidate.cosPaXYV0ToCasc, kfOmegac0Candidate.cosPaXYCascToOmegac,
                      kfOmegac0Candidate.rapOmegac, kfOmegac0Candidate.ptPiFromOmegac, kfOmegac0Candidate.ptOmegac,
                      kfOmegac0Candidate.cosThetaStarPiFromOmegac,
                      v0NDF, cascNDF, charmbaryonNDF, v0Ndfm, cascNdfm,
                      v0Chi2OverNdf, cascChi2OverNdf, charmbaryonChi2OverNdf, v0Chi2OverNdfm, cascChi2OverNdfm, kfOmegac0Candidate.cascRejectInvmass);

    } // loop over LF Cascade-bachelor candidates
  } // end of run function
  //==========================================================
  template <int decayChannel, typename Coll, typename Hist>
  void runKfXic0CreatorWithKFParticle(Coll const&,
                                      aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                                      MyKfTracks const&,
                                      MyKfCascTable const&, KFCascadesLinked const&,
                                      aod::HfCascLf2Prongs const& candidates,
                                      Hist& hInvMassCharmBaryon,
                                      Hist& hFitterStatus,
                                      Hist& hCandidateCounter,
                                      Hist& hCascadesCounter)
  {
    for (const auto& cand : candidates) {
      hCandidateCounter->Fill(1);

      auto collision = cand.collision_as<Coll>();

      // set the magnetic field from CCDB
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << magneticField;
        runNumber = bc.runNumber();
      }
      df.setBz(magneticField);
      KFParticle::SetField(magneticField);
      // bachelor from Xic0
      auto trackCharmBachelor = cand.prong0_as<MyKfTracks>();

      auto cascAodElement = cand.cascade_as<aod::KFCascadesLinked>();
      hCascadesCounter->Fill(0);
      int v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<MyKfCascTable>();
      hCascadesCounter->Fill(1);
      auto trackCascDauCharged = casc.bachelor_as<MyKfTracks>(); // pion <- xi track
      auto trackV0Dau0 = casc.posTrack_as<MyKfTracks>();         // V0 positive daughter track
      auto trackV0Dau1 = casc.negTrack_as<MyKfTracks>();         // V0 negative daughter track

      auto bachCharge = trackCascDauCharged.signed1Pt() > 0 ? +1 : -1;

      //// pion & p TrackParCov
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);
      // pion <- casc TrackParCov
      auto xiDauChargedTrackParCov = getTrackParCov(trackCascDauCharged);

      // convert tracks into KFParticle object
      KFPTrack kfTrack0 = createKFPTrackFromTrack(trackV0Dau0);
      KFPTrack kfTrack1 = createKFPTrackFromTrack(trackV0Dau1);
      KFPTrack kfTrackBach = createKFPTrackFromTrack(trackCascDauCharged);

      KFParticle kfPosPr(kfTrack0, kProton);
      KFParticle kfNegPi(kfTrack1, kPiMinus);
      KFParticle kfNegBachPi(kfTrackBach, kPiMinus);
      KFParticle kfPosPi(kfTrack0, kPiPlus);
      KFParticle kfNegPr(kfTrack1, kProton);
      KFParticle kfPosBachPi(kfTrackBach, kPiPlus);

      KFParticle kfBachPion;
      KFParticle kfPos;
      KFParticle kfNeg;
      if (bachCharge < 0) {
        kfPos = kfPosPr;
        kfNeg = kfNegPi;
        kfBachPion = kfNegBachPi;
      } else {
        kfPos = kfPosPi;
        kfNeg = kfNegPr;
        kfBachPion = kfPosBachPi;
      }

      //__________________________________________
      //*>~<* step 1 : construct V0 with KF
      const KFParticle* v0Daughters[2] = {&kfPos, &kfNeg};
      // construct V0
      KFParticle kfV0;
      kfV0.SetConstructMethod(kfConstructMethod);
      try {
        kfV0.Construct(v0Daughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct cascade V0 from daughter tracks: " << e.what();
        continue;
      }

      // mass window cut on lambda before mass constraint
      float massLam, sigLam;
      kfV0.GetMass(massLam, sigLam);
      if (TMath::Abs(massLam - MassLambda0) > lambdaMassWindow)
        continue;

      // err_mass>0 of Lambda
      if (sigLam <= 0)
        continue;
      // chi2>0 && NDF>0 for selecting Lambda
      if ((kfV0.GetNDF() <= 0 || kfV0.GetChi2() <= 0))
        continue;

      kfXic0Candidate.chi2GeoV0 = kfV0.GetChi2();
      KFParticle kfV0MassConstrained = kfV0;
      kfV0MassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassLambda); // set mass constrain to Lambda
      if (kfUseV0MassConstraint) {
        KFParticle kfV0 = kfV0MassConstrained;
      }
      kfV0.TransportToDecayVertex();

      //__________________________________________
      //*>~<* step 2 : reconstruct cascade(Xi) with KF
      const KFParticle* xiDaugthers[2] = {&kfBachPion, &kfV0};
      // construct cascade
      KFParticle kfXi;
      kfXi.SetConstructMethod(kfConstructMethod);
      try {
        kfXi.Construct(xiDaugthers, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Xi from V0 and bachelor track: " << e.what();
        continue;
      }

      float massCasc, sigCasc;
      kfXi.GetMass(massCasc, sigCasc);
      // err_massXi > 0
      if (sigCasc <= 0)
        continue;

      if (std::abs(massCasc - MassXiMinus) > massToleranceCascade)
        continue;
      // chi2>0 && NDF>0
      if (kfXi.GetNDF() <= 0 || kfXi.GetChi2() <= 0)
        continue;
      kfXic0Candidate.chi2GeoCasc = kfXi.GetChi2();
      KFParticle kfXiMassConstrained = kfXi;
      kfXiMassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassXiMinus); // set mass constrain to XiMinus
      if (kfUseCascadeMassConstraint) {
        // set mass constraint if requested
        KFParticle kfXi = kfXiMassConstrained;
      }
      registry.fill(HIST("hInvMassXiMinus"), massCasc);
      kfXi.TransportToDecayVertex();

      //__________________________________________
      //*>~<* step 3 : reconstruc Xic0 with KF
      // Create KF charm bach Pion from track
      KFPTrack kfTrackBachPion = createKFPTrackFromTrack(trackCharmBachelor);
      KFParticle kfCharmBachPion(kfTrackBachPion, kPiPlus);
      const KFParticle* xiC0Daugthers[2] = {&kfCharmBachPion, &kfXi};

      // construct XiC0
      KFParticle kfXiC0;
      kfXiC0.SetConstructMethod(kfConstructMethod);
      try {
        kfXiC0.Construct(xiC0Daugthers, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct XiC0 from Cascade and bachelor pion track: " << e.what();
        continue;
      }
      float massXiC0, sigXiC0;
      kfXiC0.GetMass(massXiC0, sigXiC0);
      if (sigXiC0 <= 0)
        continue;
      // chi2>0 && NDF>0
      if (kfXiC0.GetNDF() <= 0 || kfXiC0.GetChi2() <= 0)
        continue;

      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);
      kfXiC0.TransportToDecayVertex();
      // PV
      KFPVertex kfVertex = createKFPVertexFromCollision(collision);
      KFParticle kfPV(kfVertex);

      // set production vertex;
      kfNeg.SetProductionVertex(kfV0);
      kfPos.SetProductionVertex(kfV0);

      KFParticle kfBachPionToXi = kfBachPion;
      KFParticle kfV0ToCasc = kfV0;
      kfBachPionToXi.SetProductionVertex(kfXi);
      kfV0ToCasc.SetProductionVertex(kfXi);

      KFParticle kfXiToXiC = kfXi;
      KFParticle kfCharmBachPionToXiC = kfCharmBachPion;
      kfCharmBachPionToXiC.SetProductionVertex(kfXiC0);
      kfXiToXiC.SetProductionVertex(kfXiC0);

      // KFParticle to PV
      KFParticle kfV0ToPv = kfV0;
      KFParticle kfXiToPv = kfXi;
      KFParticle kfXic0ToPv = kfXiC0;
      KFParticle kfPiFromXicToPv = kfCharmBachPion;

      kfV0ToPv.SetProductionVertex(kfPV);
      kfXiToPv.SetProductionVertex(kfPV);
      kfXic0ToPv.SetProductionVertex(kfPV);
      kfPiFromXicToPv.SetProductionVertex(kfPV);
      //------------get updated daughter tracks after vertex fit  ---------------
      auto trackParVarCharmBachelor = getTrackParCovFromKFP(kfCharmBachPionToXiC, o2::track::PID::Pion, -bachCharge); // chrambaryon bach pion
      trackParVarCharmBachelor.setAbsCharge(1);

      xiDauChargedTrackParCov = getTrackParCovFromKFP(kfBachPionToXi, o2::track::PID::Pion, bachCharge); // Cascade bach pion
      xiDauChargedTrackParCov.setAbsCharge(1);
      o2::track::TrackParCov trackCasc = getTrackParCovFromKFP(kfXiToXiC, kfXiToXiC.GetPDG(), bachCharge);
      trackCasc.setAbsCharge(1);

      trackParCovV0Dau0 = getTrackParCovFromKFP(kfPos, kfPos.GetPDG(), 1); // V0 postive daughter
      trackParCovV0Dau0.setAbsCharge(1);
      trackParCovV0Dau1 = getTrackParCovFromKFP(kfNeg, kfNeg.GetPDG(), -1); // V0 negtive daughter
      trackParCovV0Dau1.setAbsCharge(1);

      //-------------------------- V0 info---------------------------
      // pseudorapidity
      float pseudorapV0Dau0 = kfPos.GetEta();
      float pseudorapV0Dau1 = kfNeg.GetEta();

      // info from from KFParticle
      std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      std::array<float, 3> pVecV0Dau0 = {kfPos.GetPx(), kfPos.GetPy(), kfPos.GetPz()};
      std::array<float, 3> pVecV0Dau1 = {kfNeg.GetPx(), kfNeg.GetPy(), kfNeg.GetPz()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      float pseudorapCascBachelor = kfBachPionToXi.GetEta();

      // info from KFParticle
      std::array<float, 3> vertexCasc = {kfXi.GetX(), kfXi.GetY(), kfXi.GetZ()};
      std::array<float, 3> pVecCascBachelor = {kfBachPionToXi.GetPx(), kfBachPionToXi.GetPy(), kfBachPionToXi.GetPz()};

      auto primaryVertex = getPrimaryVertex(collision);
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> vertexCharmBaryonFromFitter = {0.0, 0.0, 0.0}; // This variable get from DCAfitter in default process, in KF process it is set as 0.
      std::array<float, 3> pVecCharmBachelorAsD;
      pVecCharmBachelorAsD[0] = kfCharmBachPionToXiC.GetPx();
      pVecCharmBachelorAsD[1] = kfCharmBachPionToXiC.GetPy();
      pVecCharmBachelorAsD[2] = kfCharmBachPionToXiC.GetPz();

      std::array<float, 3> pVecCharmBaryon = {kfXiC0.GetPx(), kfXiC0.GetPy(), kfXiC0.GetPz()};
      std::array<float, 3> coordVtxCharmBaryon = {kfXiC0.GetX(), kfXiC0.GetY(), kfXiC0.GetZ()};
      auto covVtxCharmBaryon = kfXiC0.CovarianceMatrix();
      float covMatrixPV[6];
      kfVertex.GetCovarianceMatrix(covMatrixPV);

      // impact parameters
      gpu::gpustd::array<float, 2> impactParameterV0Dau0;
      gpu::gpustd::array<float, 2> impactParameterV0Dau1;
      gpu::gpustd::array<float, 2> impactParameterKaFromCasc;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, xiDauChargedTrackParCov, 2.f, matCorr, &impactParameterKaFromCasc);
      float dcaxyV0Dau0 = impactParameterV0Dau0[0];
      float dcaxyV0Dau1 = impactParameterV0Dau1[0];
      float dcaxyCascBachelor = impactParameterKaFromCasc[0];
      float dcazV0Dau0 = impactParameterV0Dau0[1];
      float dcazV0Dau1 = impactParameterV0Dau1[1];
      float dcazCascBachelor = impactParameterKaFromCasc[1];

      // pseudorapidity
      float pseudorapCharmBachelor = kfCharmBachPionToXiC.GetEta();

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterCharmBachelor;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);
      float impactParBachFromCharmBaryonXY = impactParameterCharmBachelor.getY();
      float impactParBachFromCharmBaryonZ = impactParameterCharmBachelor.getZ();

      // computing decay length and ctau
      float decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      float decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      float decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array{kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
      auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
      auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

      // fill test histograms
      hInvMassCharmBaryon->Fill(massXiC0);
      hCandidateCounter->Fill(3);

      //// KFParticle table information
      // KF chi2
      auto v0NDF = kfV0.GetNDF();
      auto v0Chi2OverNdf = kfXic0Candidate.chi2GeoV0 / v0NDF;

      auto cascNDF = kfXi.GetNDF();
      auto cascChi2OverNdf = kfXic0Candidate.chi2GeoCasc / cascNDF;

      kfXic0Candidate.chi2GeoXic = kfXiC0.GetChi2();
      auto charmbaryonNDF = kfXiC0.GetNDF();
      auto charmbaryonChi2OverNdf = kfXic0Candidate.chi2GeoXic / charmbaryonNDF;

      kfXic0Candidate.chi2MassV0 = kfV0MassConstrained.GetChi2();
      auto v0Ndfm = kfV0MassConstrained.GetNDF();
      auto v0Chi2OverNdfm = kfXic0Candidate.chi2MassV0 / v0Ndfm;

      kfXic0Candidate.chi2MassCasc = kfXiMassConstrained.GetChi2();
      auto cascNdfm = kfXiMassConstrained.GetNDF();
      auto cascChi2OverNdfm = kfXic0Candidate.chi2MassCasc / cascNdfm;

      // KF topo Chi2
      kfXic0Candidate.chi2TopoV0ToPv = kfV0ToPv.GetChi2();
      kfXic0Candidate.chi2TopoCascToPv = kfXiToPv.GetChi2();
      kfXic0Candidate.chi2TopoPiFromXicToPv = kfPiFromXicToPv.GetChi2();
      kfXic0Candidate.chi2TopoXicToPv = kfXic0ToPv.GetChi2();

      auto cascBachTopoChi2 = kfBachPionToXi.GetChi2();
      kfXic0Candidate.chi2TopoV0ToCasc = kfV0ToCasc.GetChi2();
      kfXic0Candidate.chi2TopoCascToXic = kfXiToXiC.GetChi2();

      // KF ldl
      kfXic0Candidate.ldlV0 = ldlFromKF(kfV0, kfPV);
      kfXic0Candidate.ldlCasc = ldlFromKF(kfXi, kfPV);
      kfXic0Candidate.ldlXic = ldlFromKF(kfXiC0, kfPV);

      // KF dca
      kfXic0Candidate.kfDcaXYPiFromXic = kfCharmBachPionToXiC.GetDistanceFromVertexXY(kfPV);
      kfXic0Candidate.kfDcaV0Dau = kfNeg.GetDistanceFromParticle(kfPos);
      kfXic0Candidate.kfDcaCascDau = kfBachPionToXi.GetDistanceFromParticle(kfV0ToCasc);
      kfXic0Candidate.kfDcaXYCascToPv = kfXiToXiC.GetDistanceFromVertexXY(kfPV);
      kfXic0Candidate.kfDcaXicDau = kfCharmBachPionToXiC.GetDistanceFromParticle(kfXiToXiC);

      // KF decay length
      float decayLxyLam, errDecayLxyLam;
      kfV0ToCasc.GetDecayLengthXY(decayLxyLam, errDecayLxyLam);
      kfXic0Candidate.decayLenXYLambda = decayLxyLam;

      float decayLxyCasc, errDecayLxyCasc;
      kfXiToXiC.GetDecayLengthXY(decayLxyCasc, errDecayLxyCasc);
      kfXic0Candidate.decayLenXYCasc = decayLxyCasc;

      float decayLxyXic0, errDecayLxyXic0;
      kfXic0ToPv.GetDecayLengthXY(decayLxyXic0, errDecayLxyXic0);
      kfXic0Candidate.decayLenXYXic = decayLxyXic0;

      // KF cosPA
      kfXic0Candidate.cosPaV0ToPv = cpaFromKF(kfV0, kfPV);
      kfXic0Candidate.cosPaCascToPv = cpaFromKF(kfXi, kfPV);
      kfXic0Candidate.cosPaXicToPv = cpaFromKF(kfXiC0, kfPV);
      kfXic0Candidate.cosPaXYV0ToPv = cpaXYFromKF(kfV0, kfPV);
      kfXic0Candidate.cosPaXYCascToPv = cpaXYFromKF(kfXi, kfPV);
      kfXic0Candidate.cosPaXYXicToPv = cpaXYFromKF(kfXiC0, kfPV);

      kfXic0Candidate.cosPaV0ToCasc = cpaFromKF(kfV0, kfXi);
      kfXic0Candidate.cosPaCascToXic = cpaFromKF(kfXi, kfXiC0);
      kfXic0Candidate.cosPaXYV0ToCasc = cpaXYFromKF(kfV0, kfXi);
      kfXic0Candidate.cosPaXYCascToXic = cpaXYFromKF(kfXi, kfXiC0);
      // KF mass
      kfXic0Candidate.massV0 = massLam;
      kfXic0Candidate.massCasc = massCasc;
      kfXic0Candidate.massXic = massXiC0;

      // KF pT
      kfXic0Candidate.ptPiFromXic = kfCharmBachPionToXiC.GetPt();
      kfXic0Candidate.ptXic = kfXiC0.GetPt();

      // KF rapidity
      kfXic0Candidate.rapXic = kfXiC0.GetRapidity();

      // KF cosThetaStar
      kfXic0Candidate.cosThetaStarPiFromXic = cosThetaStarFromKF(0, 4332, 211, 3312, kfCharmBachPionToXiC, kfXiToXiC);

      // KF ct
      kfXic0Candidate.ctV0 = kfV0.GetLifeTime();
      kfXic0Candidate.ctCasc = kfXi.GetLifeTime();
      kfXic0Candidate.ctXic = kfXiC0.GetLifeTime();
      kfXic0Candidate.ctOmegac = kfXiC0.GetLifeTime();

      // KF eta
      kfXic0Candidate.etaXic = kfXiC0.GetEta();

      // fill KF hist
      registry.fill(HIST("hKFParticleCascBachTopoChi2"), cascBachTopoChi2);
      registry.fill(HIST("hKFParticleV0TopoChi2"), kfXic0Candidate.chi2TopoV0ToCasc);
      registry.fill(HIST("hKFParticleCascTopoChi2"), kfXic0Candidate.chi2TopoCascToXic);
      registry.fill(HIST("hKFParticleDcaCharmBaryonDau"), kfXic0Candidate.kfDcaXicDau);
      registry.fill(HIST("hKFParticleDcaXYCascBachToPv"), dcaxyCascBachelor);
      registry.fill(HIST("hKFParticleDcaXYV0DauToPv"), dcaxyV0Dau0);
      registry.fill(HIST("hKfLambda_ldl"), kfXic0Candidate.ldlV0);
      registry.fill(HIST("hKfXi_ldl"), kfXic0Candidate.ldlCasc);
      registry.fill(HIST("hKfXiC0_ldl"), kfXic0Candidate.ldlXic);
      registry.fill(HIST("hDcaXYCascadeToPVKf"), kfXic0Candidate.kfDcaXYCascToPv);

      // fill kf table
      kfCandidateXicData(collision.globalIndex(),
                         pvCoord[0], pvCoord[1], pvCoord[2],
                         vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                         vertexCasc[0], vertexCasc[1], vertexCasc[2],
                         vertexV0[0], vertexV0[1], vertexV0[2],
                         trackCascDauCharged.sign(),
                         covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                         pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                         kfXiToXiC.GetPx(), kfXiToXiC.GetPy(), kfXiToXiC.GetPz(),
                         pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                         pVecV0[0], pVecV0[1], pVecV0[2],
                         pVecCascBachelor[0], pVecCascBachelor[1], pVecCascBachelor[2],
                         pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                         pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                         impactParameterCasc.getY(), impactParBachFromCharmBaryonXY,
                         impactParameterCasc.getZ(), impactParBachFromCharmBaryonZ,
                         std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                         v0index, casc.posTrackId(), casc.negTrackId(),
                         casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                         kfXic0Candidate.massV0, kfXic0Candidate.massCasc, kfXic0Candidate.massXic,
                         kfXic0Candidate.cosPaV0ToPv, kfXic0Candidate.cosPaXicToPv, kfXic0Candidate.cosPaCascToPv, kfXic0Candidate.cosPaXYV0ToPv, kfXic0Candidate.cosPaXYXicToPv, kfXic0Candidate.cosPaXYCascToPv,
                         kfXic0Candidate.ctOmegac, kfXic0Candidate.ctCasc, kfXic0Candidate.ctV0, kfXic0Candidate.ctXic,
                         pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                         kfXic0Candidate.etaXic, kfXi.GetEta(), kfV0.GetEta(),
                         dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                         dcazV0Dau0, dcazV0Dau1, dcazCascBachelor,
                         kfXic0Candidate.kfDcaCascDau, kfXic0Candidate.kfDcaV0Dau, kfXic0Candidate.kfDcaXicDau,
                         decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon,
                         kfXic0Candidate.kfDcaXYPiFromXic, kfXic0Candidate.kfDcaXYCascToPv,
                         kfXic0Candidate.chi2GeoV0, kfXic0Candidate.chi2GeoCasc, kfXic0Candidate.chi2GeoXic, kfXic0Candidate.chi2MassV0, kfXic0Candidate.chi2MassCasc,
                         kfXic0Candidate.ldlV0, kfXic0Candidate.ldlCasc, kfXic0Candidate.ldlXic,
                         kfXic0Candidate.chi2TopoV0ToPv, kfXic0Candidate.chi2TopoCascToPv, kfXic0Candidate.chi2TopoPiFromXicToPv, kfXic0Candidate.chi2TopoXicToPv,
                         kfXic0Candidate.chi2TopoV0ToCasc, kfXic0Candidate.chi2TopoCascToXic,
                         kfXic0Candidate.decayLenXYLambda, kfXic0Candidate.decayLenXYCasc, kfXic0Candidate.decayLenXYXic,
                         kfXic0Candidate.cosPaV0ToCasc, kfXic0Candidate.cosPaCascToXic, kfXic0Candidate.cosPaXYV0ToCasc, kfXic0Candidate.cosPaXYCascToXic,
                         kfXic0Candidate.rapXic, kfXic0Candidate.ptPiFromXic, kfXic0Candidate.ptXic,
                         kfXic0Candidate.cosThetaStarPiFromXic,
                         v0NDF, cascNDF, charmbaryonNDF, v0Ndfm, cascNdfm,
                         v0Chi2OverNdf, cascChi2OverNdf, charmbaryonChi2OverNdf, v0Chi2OverNdfm, cascChi2OverNdfm);

    } // loop over LF Cascade-bachelor candidates
  }
  /// @brief process function w/o centrality selections
  void processNoCentToXiPi(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                           aod::BCsWithTimestamps const& bcWithTimeStamps,
                           TracksWCovDca const& tracks,
                           MyLFTracksWCov const& lfTracks,
                           MyCascTable const& cascades,
                           CascadesLinked const& cascadeLinks,
                           aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToXiPi, "Run candidate creator w/o centrality selections for xi pi decay channel", true);

  void processNoCentToXiPiTraCasc(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                  aod::BCsWithTimestamps const& bcWithTimeStamps,
                                  TracksWCovDca const& tracks,
                                  MyLFTracksWCov const& lfTracks,
                                  MyTraCascTable const& traCascades,
                                  TraCascadesLinked const& traCascadeLinks,
                                  aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, lfTracks, tracks, traCascades, traCascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToXiPiTraCasc, "Run candidate creator w/o centrality selections for xi pi decay channel with tracked cascades", false);

  void processNoCentToOmegaPi(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                              aod::BCsWithTimestamps const& bcWithTimeStamps,
                              TracksWCovDca const& tracks,
                              MyLFTracksWCov const& lfTracks,
                              MyCascTable const& cascades,
                              CascadesLinked const& cascadeLinks,
                              aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToOmegaPi, "Run candidate creator w/o centrality selections for omega pi decay channel", false);

  void processOmegacToOmegaPiWithKFParticle(aod::Collisions const& collisions,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps,
                                            MyKfTracks const& tracks,
                                            MyKfCascTable const& cascades,
                                            KFCascadesLinked const& cascadeLinks,
                                            aod::HfCascLf2Prongs const& candidates)
  {
    runKfOmegac0CreatorWithKFParticle<hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processOmegacToOmegaPiWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega pi decay channel using KFParticle", false);

  void processXicToXiPiWithKFParticle(aod::Collisions const& collisions,
                                      aod::BCsWithTimestamps const& bcWithTimeStamps,
                                      MyKfTracks const& tracks,
                                      MyKfCascTable const& cascades,
                                      KFCascadesLinked const& cascadeLinks,
                                      aod::HfCascLf2Prongs const& candidates)
  {
    runKfXic0CreatorWithKFParticle<hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processXicToXiPiWithKFParticle, "Run candidate creator w/o centrality selections for Xic0 To Xi pi decay channel using KFParticle", false);

  void processNoCentToOmegaK(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                             aod::BCsWithTimestamps const& bcWithTimeStamps,
                             TracksWCovDca const& tracks,
                             MyLFTracksWCov const& lfTracks,
                             MyCascTable const& cascades,
                             CascadesLinked const& cascadeLinks,
                             aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentToOmegaK, "Run candidate creator w/o centrality selections for omega K decay channel", false);

  /// @brief process function w/ FT0C centrality selections
  void processCentFT0CToXiPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                             aod::BCsWithTimestamps const& bcWithTimeStamps,
                             TracksWCovDca const& tracks,
                             MyLFTracksWCov const& lfTracks,
                             MyCascTable const& cascades,
                             CascadesLinked const& cascadeLinks,
                             aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CToXiPi, "Run candidate creator w/ centrality selection on FT0C for xi pi channel", false);

  void processCentFT0CToOmegaPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                TracksWCovDca const& tracks,
                                MyLFTracksWCov const& lfTracks,
                                MyCascTable const& cascades,
                                CascadesLinked const& cascadeLinks,
                                aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CToOmegaPi, "Run candidate creator w/ centrality selection on FT0C for omega pi channel", false);

  void processCentFT0CToOmegaK(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                               aod::BCsWithTimestamps const& bcWithTimeStamps,
                               TracksWCovDca const& tracks,
                               MyLFTracksWCov const& lfTracks,
                               MyCascTable const& cascades,
                               CascadesLinked const& cascadeLinks,
                               aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CToOmegaK, "Run candidate creator w/ centrality selection on FT0C for omega K channel", false);

  /// @brief process function w/ FT0M centrality selections
  void processCentFT0MToXiPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                             aod::BCsWithTimestamps const& bcWithTimeStamps,
                             TracksWCovDca const& tracks,
                             MyLFTracksWCov const& lfTracks,
                             MyCascTable const& cascades,
                             CascadesLinked const& cascadeLinks,
                             aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MToXiPi, "Run candidate creator w/ centrality selection on FT0M for xi pi channel", false);

  void processCentFT0MToOmegaPi(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                TracksWCovDca const& tracks,
                                MyLFTracksWCov const& lfTracks,
                                MyCascTable const& cascades,
                                CascadesLinked const& cascadeLinks,
                                aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MToOmegaPi, "Run candidate creator w/ centrality selection on FT0M for omega pi channel", false);

  void processCentFT0MToOmegaK(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                               aod::BCsWithTimestamps const& bcWithTimeStamps,
                               TracksWCovDca const& tracks,
                               MyLFTracksWCov const& lfTracks,
                               MyCascTable const& cascades,
                               CascadesLinked const& cascadeLinks,
                               aod::HfCascLf2Prongs const& candidates)
  {
    runXic0Omegac0Creator<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, lfTracks, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MToOmegaK, "Run candidate creator w/ centrality selection on FT0M for omega K channel", false);

  ///////////////////////////////////////////////////////////
  ///                                                     ///
  ///   Process functions only for collision monitoring   ///
  ///                                                     ///
  ///////////////////////////////////////////////////////////

  /// @brief process function to monitor collisions - no centrality
  void processCollisions(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCollisions, "Collision monitoring - no centrality", true);

  /// @brief process function to monitor collisions - FT0C centrality
  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCollisionsCentFT0C, "Collision monitoring - FT0C centrality", false);

  /// @brief process function to monitor collisions - FT0M centrality
  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions, aod::BCsWithTimestamps const& /*bcWithTimeStamps*/)
  {
    /// loop over collisions
    for (const auto& collision : collisions) {

      /// bitmask with event. selection info
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);

    } /// end loop over collisions
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCollisionsCentFT0M, "Collision monitoring - FT0M centrality", false);

}; // end of struct

/// Performs MC matching.
struct HfCandidateCreatorXic0Omegac0Mc {
  Produces<aod::HfXicToXiPiMCRec> rowMCMatchRecXicToXiPi;
  Produces<aod::HfXicToXiPiMCGen> rowMCMatchGenXicToXiPi;
  Produces<aod::HfOmegacToXiPiMCRec> rowMCMatchRecOmegacToXiPi;
  Produces<aod::HfOmegacToXiPiMCGen> rowMCMatchGenOmegacToXiPi;
  Produces<aod::HfToOmegaPiMCRec> rowMCMatchRecToOmegaPi;
  Produces<aod::HfToOmegaPiMCGen> rowMCMatchGenToOmegaPi;
  Produces<aod::HfToOmegaKMCRec> rowMCMatchRecToOmegaK;
  Produces<aod::HfToOmegaKMCGen> rowMCMatchGenToOmegaK;

  // Configuration
  o2::framework::Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"};

  using MyTracksWMc = soa::Join<TracksIU, McTrackLabels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  std::shared_ptr<TH1> hGenCharmBaryonPtRapidityTightXicToXiPi, hGenCharmBaryonPtRapidityLooseXicToXiPi, hGenCharmBaryonPtRapidityTightOmegacToXiPi, hGenCharmBaryonPtRapidityLooseOmegacToXiPi, hGenCharmBaryonPtRapidityTightOmegacToOmegaPi, hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi, hGenCharmBaryonPtRapidityTightOmegacToOmegaK, hGenCharmBaryonPtRapidityLooseOmegacToOmegaK;

  HistogramRegistry registry{"registry"};

  // inspect for which zPvPosMax cut was set for reconstructed
  void init(InitContext& initContext)
  {
    std::array<bool, 4> procCollisionsXicToXiPi{doprocessMcXicToXiPi, doprocessMcXicToXiPiFT0m, doprocessMcXicToXiPiFT0c, doprocessMcXicToXiPiKf};
    if (std::accumulate(procCollisionsXicToXiPi.begin(), procCollisionsXicToXiPi.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for XicToXiPi collision study can be enabled at a time.");
    }
    std::array<bool, 3> procCollisionsOmegacToXiPi{doprocessMcOmegacToXiPi, doprocessMcOmegacToXiPiFT0m, doprocessMcOmegacToXiPiFT0c};
    if (std::accumulate(procCollisionsOmegacToXiPi.begin(), procCollisionsOmegacToXiPi.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for OmegacToXiPi collision study can be enabled at a time.");
    }
    std::array<bool, 3> procCollisionsOmegacToOmegaPi{doprocessMcOmegacToOmegaPi, doprocessMcOmegacToOmegaPiFT0m, doprocessMcOmegacToOmegaPiFT0c};
    if (std::accumulate(procCollisionsOmegacToOmegaPi.begin(), procCollisionsOmegacToOmegaPi.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for OmegacToOmegaPi collision study can be enabled at a time.");
    }
    std::array<bool, 3> procCollisionsOmegacToOmegaK{doprocessMcOmegacToOmegaK, doprocessMcOmegacToOmegaKFT0m, doprocessMcOmegacToOmegaKFT0c};
    if (std::accumulate(procCollisionsOmegacToOmegaK.begin(), procCollisionsOmegacToOmegaK.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for OmegacToOmegaK collision study can be enabled at a time.");
    }

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-xic0-omegac0") == 0) {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }
    hfEvSelMc.addHistograms(registry); // particles monitoring

    hGenCharmBaryonPtRapidityTightXicToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityTightXicToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}}); // keep track of generated candidates pt when |y|<0.5
    hGenCharmBaryonPtRapidityLooseXicToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseXicToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}}); // keep track of generated candidates pt when |y|<0.8

    hGenCharmBaryonPtRapidityTightOmegacToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityTightOmegacToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    hGenCharmBaryonPtRapidityLooseOmegacToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseOmegacToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});

    hGenCharmBaryonPtRapidityTightOmegacToOmegaPi = registry.add<TH1>("hGenCharmBaryonPtRapidityTightOmegacToOmegaPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});

    hGenCharmBaryonPtRapidityTightOmegacToOmegaK = registry.add<TH1>("hGenCharmBaryonPtRapidityTightOmegacToOmegaK", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    hGenCharmBaryonPtRapidityLooseOmegacToOmegaK = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseOmegacToOmegaK", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, int decayChannel, typename Colls, typename TMyRecoCand, typename McCollisions>
  void runXic0Omegac0Mc(TMyRecoCand const& candidates,
                        MyTracksWMc const&,
                        aod::McParticles const& mcParticles,
                        Colls const& collsWithMcLabels,
                        McCollisions const& mcCollisions,
                        BCsInfo const&)
  {
    float ptCharmBaryonGen = -999.;
    float rapidityCharmBaryonGen = -999.;
    int indexRec = -1;
    int indexRecCharmBaryon = -1;
    int8_t sign = -9;
    int8_t signCasc = -9;
    int8_t signV0 = -9;
    int8_t flag = 0;
    int8_t origin = 0; // to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenCasc = 0;
    int8_t debugGenLambda = 0;
    bool collisionMatched = false;

    int pdgCodeOmegac0 = Pdg::kOmegaC0;  // 4332
    int pdgCodeXic0 = Pdg::kXiC0;        // 4132
    int pdgCodeXiMinus = kXiMinus;       // 3312
    int pdgCodeOmegaMinus = kOmegaMinus; // 3334
    int pdgCodeLambda = kLambda0;        // 3122
    int pdgCodePiPlus = kPiPlus;         // 211
    int pdgCodePiMinus = kPiMinus;       // -211
    int pdgCodeProton = kProton;         // 2212
    int pdgCodeKaonPlus = kKPlus;        // 321
    int pdgCodeKaonMinus = kKMinus;      // -321

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = 0;
      collisionMatched = false;
      std::vector<int> idxBhadMothers{};

      auto arrayDaughters = std::array{candidate.template bachelorFromCharmBaryon_as<MyTracksWMc>(), // bachelor <- charm baryon
                                       candidate.template bachelor_as<MyTracksWMc>(),                // bachelor <- cascade
                                       candidate.template posTrack_as<MyTracksWMc>(),                // p <- lambda
                                       candidate.template negTrack_as<MyTracksWMc>()};               // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.template bachelor_as<MyTracksWMc>(),
                                           candidate.template posTrack_as<MyTracksWMc>(),
                                           candidate.template negTrack_as<MyTracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.template posTrack_as<MyTracksWMc>(),
                                         candidate.template negTrack_as<MyTracksWMc>()};

      // Check whether the particle is from background events. If so, reject it.
      if (rejectBackground) {
        bool fromBkg{false};
        for (const auto& daughter : arrayDaughters) {
          if (daughter.has_mcParticle()) {
            auto mcParticle = daughter.mcParticle();
            if (mcParticle.fromBackgroundEvent()) {
              fromBkg = true;
              break;
            }
          }
        }
        if (fromBkg) {
          rowMCMatchRecXicToXiPi(flag, debug, origin, collisionMatched, -1.f, 0);
          rowMCMatchRecOmegacToXiPi(flag, debug, origin, collisionMatched, -1.f, 0);
          rowMCMatchRecToOmegaPi(flag, debug, origin, collisionMatched, -1.f, 0);
          rowMCMatchRecToOmegaK(flag, debug, origin, collisionMatched, -1.f, 0);
          continue;
        }
      }

      // Xic0 -> xi pi matching
      if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
        // Xic → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeXic0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi);
              collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecXicToXiPi(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecXicToXiPi(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Xic0ToXiPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) { // Omegac -> xi pi matching
        // Omegac → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeXiMinus, std::array{pdgCodePiMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi);
              collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecOmegacToXiPi(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecOmegacToXiPi(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Omegac0ToXiPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) { // Omegac0 -> omega pi matching
        // Omegac → pi K pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Omega- → K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeOmegaMinus, std::array{pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi);
              collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecToOmegaPi(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecToOmegaPi(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Omegac0ToOmegaPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) { // Omegac0 -> omega K matching
        // Omegac → K K pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodeKaonPlus, pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Omega- → K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgCodeOmegaMinus, std::array{pdgCodeKaonMinus, pdgCodeProton, pdgCodePiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK);
              collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
          rowMCMatchRecToOmegaK(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
        } else {
          rowMCMatchRecToOmegaK(flag, debug, origin, collisionMatched, -1.f, 0);
        }
        if (debug == 2 || debug == 3) {
          LOGF(info, "WARNING: Omegac0ToOmegaK decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      }
    } // close loop over candidates

    for (const auto& mcCollision : mcCollisions) {

      // Slice the particles table to get the particles for the current MC collision
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      // Slice the collisions table to get the collision info for the current MC collision
      float centrality{-1.f};
      uint16_t rejectionMask{0};
      if constexpr (centEstimator == CentralityEstimator::FT0C) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == CentralityEstimator::FT0M) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == CentralityEstimator::None) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      }
      hfEvSelMc.fillHistograms<centEstimator>(mcCollision, rejectionMask);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject all particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
            rowMCMatchGenXicToXiPi(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
            rowMCMatchGenOmegacToXiPi(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
            rowMCMatchGenToOmegaPi(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
            rowMCMatchGenToOmegaK(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          }
        }
        continue;
      }

      // Match generated particles.
      for (const auto& particle : mcParticlesPerMcColl) {
        ptCharmBaryonGen = -999.;
        rapidityCharmBaryonGen = -999.;
        flag = 0;
        sign = -9;
        debugGenCharmBar = 0;
        debugGenCasc = 0;
        debugGenLambda = 0;
        origin = RecoDecay::OriginType::None;
        std::vector<int> idxBhadMothers{};

        // Reject particles from background events
        if (particle.fromBackgroundEvent() && rejectBackground) {
          if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
            rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
            rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
            rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
            rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }
          continue;
        }

        if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
          //  Xic → Xi pi
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeXic0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != pdgCodeXiMinus) {
                continue;
              }
              // Xi -> Lambda pi
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                    debugGenLambda = 1;
                    flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi);
                  }
                }
              }
            }
          }
          // Check whether the charm baryon is non-prompt (from a b quark)
          if (flag != 0) {
            origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
            if (std::abs(rapidityCharmBaryonGen) < 0.5) {
              hGenCharmBaryonPtRapidityTightXicToXiPi->SetBinContent(hGenCharmBaryonPtRapidityTightXicToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightXicToXiPi->GetBinContent(hGenCharmBaryonPtRapidityTightXicToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < 0.8) {
              hGenCharmBaryonPtRapidityLooseXicToXiPi->SetBinContent(hGenCharmBaryonPtRapidityLooseXicToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseXicToXiPi->GetBinContent(hGenCharmBaryonPtRapidityLooseXicToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }

        } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
          //  Omegac → Xi pi
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeXiMinus, pdgCodePiPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != pdgCodeXiMinus) {
                continue;
              }
              // Xi -> Lambda pi
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeXiMinus, std::array{pdgCodeLambda, pdgCodePiMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                    debugGenLambda = 1;
                    flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi);
                  }
                }
              }
            }
          }
          // Check whether the charm baryon is non-prompt (from a b quark)
          if (flag != 0) {
            origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
            if (std::abs(rapidityCharmBaryonGen) < 0.5) {
              hGenCharmBaryonPtRapidityTightOmegacToXiPi->SetBinContent(hGenCharmBaryonPtRapidityTightOmegacToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightOmegacToXiPi->GetBinContent(hGenCharmBaryonPtRapidityTightOmegacToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < 0.8) {
              hGenCharmBaryonPtRapidityLooseOmegacToXiPi->SetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseOmegacToXiPi->GetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }

        } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
          //  Omegac → Omega pi
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeOmegaMinus, pdgCodePiPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != pdgCodeOmegaMinus) {
                continue;
              }
              // Omega -> Lambda K
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeOmegaMinus, std::array{pdgCodeLambda, pdgCodeKaonMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                    debugGenLambda = 1;
                    flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi);
                  }
                }
              }
            }
          }
          // Check whether the charm baryon is non-prompt (from a b quark)
          if (flag != 0) {
            origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
            if (std::abs(rapidityCharmBaryonGen) < 0.5) {
              hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->SetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->GetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < 0.8) {
              hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->SetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->GetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }

        } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
          //  Omegac → Omega K
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeOmegaMinus, pdgCodeKaonPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != pdgCodeOmegaMinus) {
                continue;
              }
              // Omega -> Lambda K
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgCodeOmegaMinus, std::array{pdgCodeLambda, pdgCodeKaonMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != pdgCodeLambda) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
                    debugGenLambda = 1;
                    flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK);
                  }
                }
              }
            }
          }
          // Check whether the charm baryon is non-prompt (from a b quark)
          if (flag != 0) {
            origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
            if (std::abs(rapidityCharmBaryonGen) < 0.5) {
              hGenCharmBaryonPtRapidityTightOmegacToOmegaK->SetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaK->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightOmegacToOmegaK->GetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaK->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < 0.8) {
              hGenCharmBaryonPtRapidityLooseOmegacToOmegaK->SetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToOmegaK->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseOmegacToOmegaK->GetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToOmegaK->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }
        }
      } // close loop on MCParticles
    } // close loop on MCCollisions
  } // close process

  void processDoNoMc(aod::Collisions::iterator const&)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processDoNoMc, "Do not run any MC process function", true);

  void processMcXicToXiPi(aod::HfCandToXiPi const& candidates,
                          MyTracksWMc const& tracks,
                          aod::McParticles const& mcParticles,
                          aod::McCollisions const& mcColls,
                          McCollisionsNoCents const& collsWithMcLabels,
                          BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcXicToXiPi, "Run Xic0 to xi pi MC process function - no centrality", false);

  void processMcXicToXiPiKf(aod::HfCandToXiPiKf const& candidates,
                            MyTracksWMc const& tracks,
                            aod::McParticles const& mcParticles,
                            aod::McCollisions const& mcColls,
                            McCollisionsNoCents const& collsWithMcLabels,
                            BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcXicToXiPiKf, "Run Xic0 to xi pi MC process function - no centrality", false);

  void processMcXicToXiPiFT0m(aod::HfCandToXiPi const& candidates,
                              MyTracksWMc const& tracks,
                              aod::McParticles const& mcParticles,
                              McCollisionsCentFT0Ms const& mcColls,
                              McCollisionsFT0Ms const& collsWithMcLabels,
                              BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcXicToXiPiFT0m, "Run Xic0 to xi pi MC process function - FT0M", false);

  void processMcXicToXiPiFT0c(aod::HfCandToXiPi const& candidates,
                              MyTracksWMc const& tracks,
                              aod::McParticles const& mcParticles,
                              aod::McCollisions const& mcColls,
                              McCollisionsFT0Cs const& collsWithMcLabels,
                              BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcXicToXiPiFT0c, "Run Xic0 to xi pi MC process function - FT0C", false);

  void processMcOmegacToXiPi(aod::HfCandToXiPi const& candidates,
                             MyTracksWMc const& tracks,
                             aod::McParticles const& mcParticles,
                             aod::McCollisions const& mcColls,
                             McCollisionsNoCents const& collsWithMcLabels,
                             BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToXiPi, "Run Omegac0 to xi pi MC process function - no centrality", false);

  void processMcOmegacToXiPiFT0m(aod::HfCandToXiPi const& candidates,
                                 MyTracksWMc const& tracks,
                                 aod::McParticles const& mcParticles,
                                 McCollisionsCentFT0Ms const& mcColls,
                                 McCollisionsFT0Ms const& collsWithMcLabels,
                                 BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToXiPiFT0m, "Run Omegac0 to xi pi MC process function - FT0M", false);

  void processMcOmegacToXiPiFT0c(aod::HfCandToXiPi const& candidates,
                                 MyTracksWMc const& tracks,
                                 aod::McParticles const& mcParticles,
                                 aod::McCollisions const& mcColls,
                                 McCollisionsFT0Cs const& collsWithMcLabels,
                                 BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToXiPiFT0c, "Run Omegac0 to xi pi MC process function - FT0C", false);

  void processMcOmegacToOmegaPi(aod::HfCandToOmegaPi const& candidates,
                                MyTracksWMc const& tracks,
                                aod::McParticles const& mcParticles,
                                aod::McCollisions const& mcColls,
                                McCollisionsNoCents const& collsWithMcLabels,
                                BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaPi, "Run Omegac0 to omega pi MC process function - no centrality", false);

  void processMcOmegacToOmegaPiFT0m(aod::HfCandToOmegaPi const& candidates,
                                    MyTracksWMc const& tracks,
                                    aod::McParticles const& mcParticles,
                                    McCollisionsCentFT0Ms const& mcColls,
                                    McCollisionsFT0Ms const& collsWithMcLabels,
                                    BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaPiFT0m, "Run Omegac0 to omega pi MC process function - FT0M", false);

  void processMcOmegacToOmegaPiFT0c(aod::HfCandToOmegaPi const& candidates,
                                    MyTracksWMc const& tracks,
                                    aod::McParticles const& mcParticles,
                                    aod::McCollisions const& mcColls,
                                    McCollisionsFT0Cs const& collsWithMcLabels,
                                    BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaPiFT0c, "Run Omegac0 to omega pi MC process function - FT0C", false);

  void processMcOmegacToOmegaK(aod::HfCandToOmegaK const& candidates,
                               MyTracksWMc const& tracks,
                               aod::McParticles const& mcParticles,
                               aod::McCollisions const& mcColls,
                               McCollisionsNoCents const& collsWithMcLabels,
                               BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaK, "Run Omegac0 to omega K MC process function - no centrality", false);

  void processMcOmegacToOmegaKFT0m(aod::HfCandToOmegaK const& candidates,
                                   MyTracksWMc const& tracks,
                                   aod::McParticles const& mcParticles,
                                   McCollisionsCentFT0Ms const& mcColls,
                                   McCollisionsFT0Ms const& collsWithMcLabels,
                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaKFT0m, "Run Omegac0 to omega K MC process function - FT0M", false);

  void processMcOmegacToOmegaKFT0c(aod::HfCandToOmegaK const& candidates,
                                   MyTracksWMc const& tracks,
                                   aod::McParticles const& mcParticles,
                                   aod::McCollisions const& mcColls,
                                   McCollisionsFT0Cs const& collsWithMcLabels,
                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcColls, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcOmegacToOmegaKFT0c, "Run Omegac0 to omega K MC process function - FT0C", false);

}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0Mc>(cfgc)};
}
