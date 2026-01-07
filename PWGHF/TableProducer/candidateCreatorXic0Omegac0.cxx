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
/// \author Ran Tu <ran.tu@cern.ch>, Fudan University
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University

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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::track;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::aod::hf_track_index;
using namespace o2::hf_centrality;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_evsel;

enum McMatchFlag : uint8_t {
  None = 0,
  CharmbaryonUnmatched,
  CascUnmatched,
  V0Unmatched
};

// Reconstruction of omegac0 and xic0 candidates
struct HfCandidateCreatorXic0Omegac0 {
  Produces<aod::HfCandToXiPi> rowCandToXiPi;
  Produces<aod::HfCandToOmegaPi> rowCandToOmegaPi;
  Produces<aod::HfCandToOmegaK> rowCandToOmegaK;
  Produces<aod::HfOmegacKf> kfCandidateData;
  Produces<aod::HfCandToXiPiKf> kfCandidateXicData;
  Produces<aod::HfCandToXiPiKfQa> rowKfXic0Qa;
  Produces<aod::HfCandToOmegaKaKf> kfCandidateOmegaKaData;

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
  Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
  Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3.0, "Max cascade DCA to PV in xy plane"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 1.0, "Max DCA of V0 daughter"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 1.0, "Max DCA of cascade daughter"};

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
  Configurable<double> massToleranceCascadeRej{"massToleranceCascadeRej", 0.01, "Invariant mass tolerance for rejected Xi"};
  // for KF particle operation
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
  Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", false, "KF: use Lambda mass constraint"};
  Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: use Cascade mass constraint"};
  Configurable<bool> kfResolutionQA{"kfResolutionQA", false, "KF: KFParticle Quality Assurance"};

  HfEventSelection hfEvSel;        // event selection and monitoring
  o2::vertexing::DCAFitterN<2> df; // 2-prong vertex fitter to build the omegac/xic vertex
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgdb;
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber{-1};
  double magneticField{0.};

  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>;
  using MyTraCascTable = soa::Join<aod::TraCascDatas, aod::TraCascCovs>; // to use strangeness tracking
  using CascadesLinked = soa::Join<Cascades, CascDataLink>;
  using TraCascadesLinked = soa::Join<Cascades, TraCascDataLink>;
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MyLFTracksWCov = soa::Join<TracksIU, TracksCovIU>;

  using MyKfTracksIU = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using MyKfTracks = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using MyKfCascTable = soa::Join<KFCascDatas, aod::KFCascCovs>;
  using KFCascadesLinked = soa::Join<aod::Cascades, aod::KFCascDataLink>;

  std::shared_ptr<TH1> hInvMassCharmBaryonToXiPi, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaK, hFitterStatusToXiPi, hFitterStatusToOmegaPi, hFitterStatusToOmegaK, hCandidateCounterToXiPi, hCandidateCounterToOmegaPi, hCandidateCounterToOmegaK, hCascadesCounterToXiPi, hCascadesCounterToOmegaPi, hCascadesCounterToOmegaK;

  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

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
    float deviationPiFromOmegacToPv;
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
  } kfOmegac0Candidate{};

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
    float cosPaCascToXic;
    float cosPaXYCascToXic;
    float cosPaCascToPv;
    float cosPaXYCascToPv;
    float massV0;
    float massCasc;
    float rapXic;
    float massXic;
    float cosThetaStarPiFromXic;
    float chi2NdfTopoPiFromXicToPv;
    float kfDcaXYPiFromXic;
    float chi2NdfTopoV0ToCasc;
    float chi2NdfTopoCascToXic;
    float decayLenXYXic;
    float chi2GeoXic;
    float kfDcaV0Dau;
    float kfDcaCascDau;
    float kfDcaXicDau;
    float kfDcaXYCascToPv;
    float chi2NdfTopoXicToPv;
    float cosPaXicToPv;
    float cosPaXYXicToPv;
    float ldlXic;
    float ctV0;
    float ctCasc;
    float ctXic;
    float chi2MassV0;
    float chi2MassCasc;
    float etaXic;
  } kfXic0Candidate{};

  void init(InitContext const&)
  {
    std::array<bool, 16> allProcesses = {doprocessNoCentToXiPi, doprocessNoCentToXiPiTraCasc, doprocessCentFT0CToXiPi, doprocessCentFT0MToXiPi, doprocessNoCentToOmegaPi, doprocessNoCentOmegacToOmegaPiWithKFParticle, doprocessCentFT0COmegacToOmegaPiWithKFParticle, doprocessCentFT0MOmegacToOmegaPiWithKFParticle, doprocessCentFT0CToOmegaPi, doprocessCentFT0MToOmegaPi, doprocessNoCentToOmegaK, doprocessCentFT0CToOmegaK, doprocessCentFT0MToOmegaK, doprocessNoCentXicToXiPiWithKFParticle, doprocessCentFT0CXicToXiPiWithKFParticle, doprocessCentFT0MXicToXiPiWithKFParticle};
    if (std::accumulate(allProcesses.begin(), allProcesses.end(), 0) == 0) {
      LOGP(fatal, "No process function enabled, please select one for at least one channel.");
    }

    std::array<bool, 7> processesToXiPi = {doprocessNoCentToXiPi, doprocessNoCentToXiPiTraCasc, doprocessCentFT0CToXiPi, doprocessCentFT0MToXiPi, doprocessNoCentXicToXiPiWithKFParticle, doprocessCentFT0CXicToXiPiWithKFParticle, doprocessCentFT0MXicToXiPiWithKFParticle};
    if (std::accumulate(processesToXiPi.begin(), processesToXiPi.end(), 0) > 1) {
      LOGP(fatal, "One and only one ToXiPi process function must be enabled at a time.");
    }
    std::array<bool, 6> processesToOmegaPi = {doprocessNoCentToOmegaPi, doprocessCentFT0CToOmegaPi, doprocessCentFT0MToOmegaPi, doprocessNoCentOmegacToOmegaPiWithKFParticle, doprocessCentFT0COmegacToOmegaPiWithKFParticle, doprocessCentFT0MOmegacToOmegaPiWithKFParticle};
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
      if ((doprocessNoCentToXiPi && !doprocessCollisions) || (doprocessNoCentToXiPiTraCasc && !doprocessCollisions) || (doprocessNoCentToOmegaPi && !doprocessCollisions) || (doprocessNoCentToOmegaK && !doprocessCollisions) || (doprocessNoCentOmegacToOmegaPiWithKFParticle && !doprocessCollisions) || (doprocessNoCentXicToXiPiWithKFParticle && !doprocessCollisions)) {
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
    registry.add("hKfChi2TopoPiFromCharmBaryon", "hKfChi2TopoPifromCharmBaryon", kTH1F, {{2000, -0.1f, 1000.0f}});
    registry.add("hKfNdfPiFromCharmBaryon", "hKfNDfPifromCharmBaryon", kTH1F, {{2000, -0.1f, 200.0f}});
    registry.add("hKfChi2OverNdfPiFromCharmBaryon", "hKfChi2OverNdfPifromCharmBaryon", kTH1F, {{1000, -0.1f, 200.0f}});
    registry.add("hKFParticlechi2TopoOmegacToPv", "hKFParticlechi2TopoOmegacToPv", kTH1D, {{1000, -0.1f, 100.0f}});
    registry.add("hKfNdfOmegacToPv", "hKfNDfOmegacToPv", kTH1F, {{2000, -0.1f, 200.0f}});
    registry.add("hKfChi2TopoOmegacToPv", "hKfChi2TopoOmegacToPv", kTH1F, {{1000, -0.1f, 200.0f}});
    registry.add("hKfDeviationPiFromCharmBaryon", "hKfDeviationPiFromCharmBaryon", kTH1F, {{1000, -0.1f, 200.0f}});
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
    registry.add("hKFParticlechi2TopoCascToPv", "hKFParticlechi2TopoCascToPv", kTH1D, {{1000, -0.1f, 100.0f}});
    registry.add("hKFParticleDcaXYV0DauPosToPv", "hKFParticleDcaXYV0DauPosToPv", kTH1D, {{1000, -0.1f, 30.0f}});
    registry.add("hKFParticleDcaXYV0DauNegToPv", "hKFParticleDcaXYV0DauNegToPv", kTH1D, {{1000, -0.1f, 30.0f}});

    // Additional KFParticle Histograms
    if (fillAllHist) {
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
    }

    // init HF event selection helper
    hfEvSel.init(registry, &zorroSummary);

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

  template <o2::hf_centrality::CentralityEstimator CentEstimator, int DecayChannel, typename Coll, typename Hist, typename TCascTable, typename TCascLinkTable>
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

    if constexpr (DecayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi && DecayChannel != hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi && DecayChannel != hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
      LOGP(fatal, "Decay channel not recognized!");
    }

    for (const auto& cand : candidates) {

      hCandidateCounter->Fill(0);

      if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
          continue;
        }
      } else if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
          continue;
        }
      } else if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK)) {
          continue;
        }
      }

      hCandidateCounter->Fill(1);

      auto collision = cand.collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      std::array<float, 3> const pVecCasc = {casc.px(), casc.py(), casc.pz()};
      std::array<float, 21> covCasc = {0.};
      constexpr int NumCovElements = 6;
      constexpr int MomInd[NumCovElements] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < NumCovElements; i++) {
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
      if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
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
      std::array<float, 3> pVecCascAsD{};
      std::array<float, 3> pVecCharmBachelorAsD{};
      df.propagateTracksToVertex();
      if (!df.isPropagateTracksToVertexDone()) {
        continue;
      }
      df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
      df.getTrack(1).getPxPyPzGlo(pVecCharmBachelorAsD);
      std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecCharmBachelorAsD[0], pVecCascAsD[1] + pVecCharmBachelorAsD[1], pVecCascAsD[2] + pVecCharmBachelorAsD[2]};

      std::array<float, 3> const coordVtxCharmBaryon = df.getPCACandidatePos();
      std::array<float, 6> covVtxCharmBaryon = df.calcPCACovMatrixFlat();

      // pseudorapidity
      float const pseudorapCharmBachelor = trackCharmBachelor.eta();

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
      float const dcaxyV0Dau0 = impactParameterV0Dau0.getY();
      float const dcaxyV0Dau1 = impactParameterV0Dau1.getY();
      float const dcaxyCascBachelor = impactParameterCascDauCharged.getY();
      float const dcazV0Dau0 = impactParameterV0Dau0.getZ();
      float const dcazV0Dau1 = impactParameterV0Dau1.getZ();
      float const dcazCascBachelor = impactParameterCascDauCharged.getZ();

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterCharmBachelor;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);
      float const impactParBachFromCharmBaryonXY = impactParameterCharmBachelor.getY();
      float const impactParBachFromCharmBaryonZ = impactParameterCharmBachelor.getZ();

      // invariant mass under the hypothesis of particles ID corresponding to the decay chain
      float mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
      float mCasc = 0.;
      if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        mCasc = casc.mXi();
      } else {
        mCasc = casc.mOmega();
      }
      auto arrMassCharmBaryon = std::array{0., 0.};
      if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        arrMassCharmBaryon = {MassXiMinus, MassPiPlus};
      } else if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        arrMassCharmBaryon = {MassOmegaMinus, MassPiPlus};
      } else {
        arrMassCharmBaryon = {MassOmegaMinus, MassKPlus};
      }
      float mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecCharmBachelorAsD}, arrMassCharmBaryon);

      // computing cosPA
      float cpaV0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float const cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      float cpaCasc = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
      float const cpaxyV0 = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float const cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
      float const cpaxyCasc = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);

      // computing decay length and ctau
      float const decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      float const decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      float const decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
      auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
      auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

      float const ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassOmegaC0);
      float const ctXic = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassXiC0);
      float ctCascade = 0.;
      if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassXiMinus);
      } else {
        ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, MassOmegaMinus);
      }
      float const ctV0 = RecoDecay::ct(pVecV0, decLenV0, MassLambda0);

      // computing eta
      float const pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
      float const pseudorapCascade = RecoDecay::eta(pVecCasc);
      float const pseudorapV0 = RecoDecay::eta(pVecV0);

      // DCA between daughters
      float dcaCascDau = casc.dcacascdaughters();
      float dcaV0Dau = casc.dcaV0daughters();
      float const dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

      // fill test histograms
      hInvMassCharmBaryon->Fill(mCharmBaryon);
      hCandidateCounter->Fill(3);

      // fill the table
      if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
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

      } else if constexpr (DecayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
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
                         decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon, cand.hfflag());

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

  template <o2::hf_centrality::CentralityEstimator CentEstimator, int DecayChannel, typename Coll, typename Hist>
  void runKfOmegac0CreatorWithKFParticle(Coll const&,
                                         aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                                         MyKfTracksIU const& tracksIU,
                                         MyKfTracks const& tracks,
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
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      KFParticle::SetField(magneticField);
      // bachelor from Omegac0
      auto trackCharmBachelorId = cand.prong0Id();
      auto trackCharmBachelor = tracks.rawIteratorAt(trackCharmBachelorId);
      auto cascAodElement = cand.cascade_as<aod::KFCascadesLinked>();
      hCascadesCounter->Fill(0);
      int const v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<MyKfCascTable>();
      hCascadesCounter->Fill(1);
      auto trackCascDauChargedId = casc.bachelorId();
      auto trackV0Dau0Id = casc.posTrackId();
      auto trackV0Dau1Id = casc.negTrackId();
      auto trackCascDauCharged = tracksIU.rawIteratorAt(trackCascDauChargedId); // pion <- xi track
      auto trackV0Dau0 = tracksIU.rawIteratorAt(trackV0Dau0Id);                 // V0 positive daughter track
      auto trackV0Dau1 = tracksIU.rawIteratorAt(trackV0Dau1Id);                 // V0 negative daughter track

      auto bachCharge = trackCascDauCharged.signed1Pt() > 0 ? +1 : -1;

      //// pion & p TrackParCov
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);
      // kaon <- casc TrackParCov
      auto omegaDauChargedTrackParCov = getTrackParCov(trackCascDauCharged);
      // convert tracks into KFParticle object
      KFPTrack const kfTrack0 = createKFPTrackFromTrack(trackV0Dau0);
      KFPTrack const kfTrack1 = createKFPTrackFromTrack(trackV0Dau1);
      KFPTrack const kfTrackBach = createKFPTrackFromTrack(trackCascDauCharged);

      KFParticle const kfPosPr(kfTrack0, kProton);
      KFParticle const kfNegPi(kfTrack1, kPiMinus);
      KFParticle const kfNegKa(kfTrackBach, kKMinus);
      KFParticle const kfNegPiRej(kfTrackBach, kPiMinus); // rej
      KFParticle const kfPosPi(kfTrack0, kPiPlus);
      KFParticle const kfNegPr(kfTrack1, kProton);
      KFParticle const kfPosKa(kfTrackBach, kKPlus);
      KFParticle const kfPosPiRej(kfTrackBach, kPiPlus); // rej

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
      if (std::abs(massLam - MassLambda0) > lambdaMassWindow) {
        continue;
      }
      // err_mass>0 of Lambda
      if (sigLam <= 0) {
        continue;
      }
      kfOmegac0Candidate.chi2GeoV0 = kfV0.GetChi2();
      KFParticle kfV0MassConstrained = kfV0;
      kfV0MassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassLambda); // set mass constrain to Lambda
      if (kfUseV0MassConstraint) {
        KFParticle const kfV0 = kfV0MassConstrained;
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
      if (sigCasc <= 0) {
        continue;
      }
      if (std::abs(massCasc - MassOmegaMinus) > massToleranceCascade) {
        continue;
      }

      kfOmegac0Candidate.chi2GeoCasc = kfOmega.GetChi2();
      kfOmegac0Candidate.cascRejectInvmass = massCascrej;
      registry.fill(HIST("hInvMassXiMinus_rej"), massCascrej); // rej
      KFParticle kfOmegaMassConstrained = kfOmega;
      kfOmegaMassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassOmegaMinus); // set mass constrain to OmegaMinus
      if (kfUseCascadeMassConstraint) {
        // set mass constraint if requested
        KFParticle const kfOmega = kfOmegaMassConstrained;
      }
      registry.fill(HIST("hInvMassOmegaMinus"), massCasc);
      kfOmega.TransportToDecayVertex();
      // rej: Add competing rejection to minimize misidentified Xi impact. Reject if kfBachPionRej is Pion and the constructed cascade has Xi's invariant mass.

      //__________________________________________
      //*>~<* step 3 : reconstruc Omegac0 with KF
      // Create KF charm bach Pion from track
      KFPTrack const kfTrackBachPion = createKFPTrackFromTrack(trackCharmBachelor);
      KFParticle const kfBachPion(kfTrackBachPion, kPiPlus);
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
      if (sigOmegaC0 <= 0) {
        continue;
      }

      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);
      kfOmegaC0.TransportToDecayVertex();
      // PV
      KFPVertex const kfVertex = createKFPVertexFromCollision(collision);
      KFParticle const kfPV(kfVertex);

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
      float const pseudorapV0Dau0 = kfPos.GetEta();
      float const pseudorapV0Dau1 = kfNeg.GetEta();

      // info from from KFParticle
      std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      std::array<float, 3> pVecV0Dau0 = {kfPos.GetPx(), kfPos.GetPy(), kfPos.GetPz()};
      std::array<float, 3> pVecV0Dau1 = {kfNeg.GetPx(), kfNeg.GetPy(), kfNeg.GetPz()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      float const pseudorapCascBachelor = kfBachKaonToOmega.GetEta();

      // info from KFParticle
      std::array<float, 3> vertexCasc = {kfOmega.GetX(), kfOmega.GetY(), kfOmega.GetZ()};
      std::array<float, 3> pVecCascBachelor = {kfBachKaonToOmega.GetPx(), kfBachKaonToOmega.GetPy(), kfBachKaonToOmega.GetPz()};

      auto primaryVertex = getPrimaryVertex(collision);
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> vertexCharmBaryonFromFitter = {0.0, 0.0, 0.0}; // This variable get from DCAfitter in default process, in KF process it is set as 0.
      std::array<float, 3> pVecCharmBachelorAsD{};
      pVecCharmBachelorAsD[0] = kfBachPionToOmegaC.GetPx();
      pVecCharmBachelorAsD[1] = kfBachPionToOmegaC.GetPy();
      pVecCharmBachelorAsD[2] = kfBachPionToOmegaC.GetPz();

      std::array<float, 3> pVecCharmBaryon = {kfOmegaC0.GetPx(), kfOmegaC0.GetPy(), kfOmegaC0.GetPz()};
      std::array<float, 3> const coordVtxCharmBaryon = {kfOmegaC0.GetX(), kfOmegaC0.GetY(), kfOmegaC0.GetZ()};
      auto* covVtxCharmBaryon = kfOmegaC0.CovarianceMatrix();
      float covMatrixPV[6];
      kfVertex.GetCovarianceMatrix(covMatrixPV);

      // impact parameters
      std::array<float, 2> impactParameterV0Dau0{};
      std::array<float, 2> impactParameterV0Dau1{};
      std::array<float, 2> impactParameterKaFromCasc{};
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, omegaDauChargedTrackParCov, 2.f, matCorr, &impactParameterKaFromCasc);
      float dcaxyV0Dau0 = impactParameterV0Dau0[0];
      float dcaxyV0Dau1 = impactParameterV0Dau1[0];
      float dcaxyCascBachelor = impactParameterKaFromCasc[0];
      float const dcazV0Dau0 = impactParameterV0Dau0[1];
      float const dcazV0Dau1 = impactParameterV0Dau1[1];
      float const dcazCascBachelor = impactParameterKaFromCasc[1];

      // pseudorapidity
      float const pseudorapCharmBachelor = kfBachPionToOmegaC.GetEta();

      // impact parameters
      o2::dataformats::DCA impactParameterCasc;
      o2::dataformats::DCA impactParameterCharmBachelor;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);
      float impactParBachFromCharmBaryonXY = impactParameterCharmBachelor.getY();
      float impactParBachFromCharmBaryonZ = impactParameterCharmBachelor.getZ();

      // computing decay length and ctau
      float const decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
      float const decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
      float const decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

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
      kfOmegac0Candidate.deviationPiFromOmegacToPv = kfCalculateChi2ToPrimaryVertex(kfOmegaC0, kfPV);

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
      kfOmegac0Candidate.cosThetaStarPiFromOmegac = cosThetaStarFromKF(0, 4332, 211, 3312, kfBachPionToOmegaC, kfOmegaToOmegaC, pdgdb);

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
      registry.fill(HIST("hKfChi2TopoOmegacToPv"), kfOmegac0ToPv.GetChi2());
      registry.fill(HIST("hKfNdfOmegacToPv"), kfOmegac0ToPv.GetNDF());
      registry.fill(HIST("hKFParticlechi2TopoCascToPv"), kfOmegac0Candidate.chi2NdfTopoCascToPv);
      registry.fill(HIST("hKFParticleDcaCharmBaryonDau"), kfOmegac0Candidate.kfDcaOmegacDau);
      registry.fill(HIST("hKFParticleDcaXYCascBachToPv"), dcaxyCascBachelor);
      registry.fill(HIST("hKFParticleDcaXYV0DauPosToPv"), dcaxyV0Dau0);
      registry.fill(HIST("hKFParticleDcaXYV0DauNegToPv"), dcaxyV0Dau1);
      registry.fill(HIST("hKfLambda_ldl"), kfOmegac0Candidate.ldlV0);
      registry.fill(HIST("hKfOmega_ldl"), kfOmegac0Candidate.ldlCasc);
      registry.fill(HIST("hKfOmegaC0_ldl"), kfOmegac0Candidate.ldlOmegac);
      registry.fill(HIST("hDcaXYCascadeToPVKf"), kfOmegac0Candidate.kfDcaXYCascToPv);
      registry.fill(HIST("hKfChi2TopoPiFromCharmBaryon"), kfPiFromOmegacToPv.GetChi2());
      registry.fill(HIST("hKfNdfPiFromCharmBaryon"), kfPiFromOmegacToPv.GetNDF());
      registry.fill(HIST("hKfChi2OverNdfPiFromCharmBaryon"), kfOmegac0Candidate.chi2NdfTopoPiFromOmegacToPv);
      registry.fill(HIST("hKfDeviationPiFromCharmBaryon"), kfOmegac0Candidate.deviationPiFromOmegacToPv);
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
                       decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon, cand.hfflag());
      // fill kf table
      kfCandidateData(kfOmegac0Candidate.kfDcaXYPiFromOmegac, kfOmegac0Candidate.kfDcaXYCascToPv,
                      kfOmegac0Candidate.chi2GeoV0, kfOmegac0Candidate.chi2GeoCasc, kfOmegac0Candidate.chi2GeoOmegac, kfOmegac0Candidate.chi2MassV0, kfOmegac0Candidate.chi2MassCasc,
                      kfOmegac0Candidate.ldlV0, kfOmegac0Candidate.ldlCasc, kfOmegac0Candidate.ldlOmegac,
                      kfOmegac0Candidate.chi2NdfTopoV0ToPv, kfOmegac0Candidate.chi2NdfTopoCascToPv, kfOmegac0Candidate.chi2NdfTopoPiFromOmegacToPv, kfOmegac0Candidate.chi2NdfTopoOmegacToPv, kfOmegac0Candidate.deviationPiFromOmegacToPv,
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
  template <o2::hf_centrality::CentralityEstimator CentEstimator, int DecayChannel, typename Coll, typename Hist>
  void runKfXic0CreatorWithKFParticle(Coll const&,
                                      aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                                      MyKfTracksIU const& tracksIU,
                                      MyKfTracks const& tracks,
                                      MyKfCascTable const&, KFCascadesLinked const&,
                                      aod::HfCascLf2Prongs const& candidates,
                                      Hist& hInvMassCharmBaryon,
                                      Hist& hFitterStatus,
                                      Hist& hCandidateCounter,
                                      Hist& hCascadesCounter)
  {
    for (const auto& cand : candidates) {
      hCandidateCounter->Fill(1);
      if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
        continue;
      }
      auto collision = cand.collision_as<Coll>();

      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
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
      KFParticle::SetField(magneticField);
      // bachelor from Xic0
      auto trackCharmBachelorId = cand.prong0Id();
      auto trackCharmBachelor = tracks.rawIteratorAt(trackCharmBachelorId);
      auto cascAodElement = cand.cascade_as<aod::KFCascadesLinked>();
      hCascadesCounter->Fill(0);
      int const v0index = cascAodElement.v0Id();
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<MyKfCascTable>();
      hCascadesCounter->Fill(1);
      auto trackCascDauChargedId = casc.bachelorId();
      auto trackV0Dau0Id = casc.posTrackId();
      auto trackV0Dau1Id = casc.negTrackId();
      auto trackCascDauCharged = tracksIU.rawIteratorAt(trackCascDauChargedId); // pion <- xi track
      auto trackV0Dau0 = tracksIU.rawIteratorAt(trackV0Dau0Id);                 // V0 positive daughter track
      auto trackV0Dau1 = tracksIU.rawIteratorAt(trackV0Dau1Id);                 // V0 negative daughter track

      auto bachCharge = trackCascDauCharged.signed1Pt() > 0 ? +1 : -1;

      //// pion & p TrackParCov
      auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
      auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);
      // pion <- casc TrackParCov
      auto xiDauChargedTrackParCov = getTrackParCov(trackCascDauCharged);

      // convert tracks into KFParticle object
      KFPTrack const kfTrack0 = createKFPTrackFromTrack(trackV0Dau0);
      KFPTrack const kfTrack1 = createKFPTrackFromTrack(trackV0Dau1);
      KFPTrack const kfTrackBach = createKFPTrackFromTrack(trackCascDauCharged);

      KFParticle const kfPosPr(kfTrack0, kProton);
      KFParticle const kfNegPi(kfTrack1, kPiMinus);
      KFParticle const kfNegBachPi(kfTrackBach, kPiMinus);
      KFParticle const kfPosPi(kfTrack0, kPiPlus);
      KFParticle const kfNegPr(kfTrack1, kProton);
      KFParticle const kfPosBachPi(kfTrackBach, kPiPlus);

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
      if (std::abs(massLam - MassLambda0) > lambdaMassWindow) {
        continue;
      }

      // err_mass>0 of Lambda
      if (sigLam <= 0) {
        continue;
      }
      // chi2>0 && NDF>0 for selecting Lambda
      if ((kfV0.GetNDF() <= 0 || kfV0.GetChi2() <= 0)) {
        continue;
      }

      kfXic0Candidate.chi2GeoV0 = kfV0.GetChi2();
      KFParticle kfV0MassConstrained = kfV0;
      kfV0MassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassLambda); // set mass constrain to Lambda
      if (kfUseV0MassConstraint) {
        kfV0 = kfV0MassConstrained;
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
      if (sigCasc <= 0) {
        continue;
      }

      if (std::abs(massCasc - MassXiMinus) > massToleranceCascade) {
        continue;
      }
      // chi2>0 && NDF>0
      if (kfXi.GetNDF() <= 0 || kfXi.GetChi2() <= 0) {
        continue;
      }
      kfXic0Candidate.chi2GeoCasc = kfXi.GetChi2();
      KFParticle kfXiMassConstrained = kfXi;
      kfXiMassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassXiMinus); // set mass constrain to XiMinus
      if (kfUseCascadeMassConstraint) {
        // set mass constraint if requested
        KFParticle const kfXi = kfXiMassConstrained;
      }
      registry.fill(HIST("hInvMassXiMinus"), massCasc);
      kfXi.TransportToDecayVertex();

      //__________________________________________
      //*>~<* step 3 : reconstruc Xic0 with KF
      // Create KF charm bach Pion from track
      KFPTrack const kfTrackBachPion = createKFPTrackFromTrack(trackCharmBachelor);
      KFParticle const kfCharmBachPion(kfTrackBachPion, kPiPlus);
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
      if (sigXiC0 <= 0) {
        continue;
      }
      // chi2>0 && NDF>0
      if (kfXiC0.GetNDF() <= 0 || kfXiC0.GetChi2() <= 0) {
        continue;
      }

      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);
      kfXiC0.TransportToDecayVertex();
      // PV
      KFPVertex const kfVertex = createKFPVertexFromCollision(collision);
      KFParticle const kfPV(kfVertex);

      KFParticle const kfPosOrigin = kfPos;
      KFParticle const kfNegOrigin = kfNeg;
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
      float const pseudorapV0Dau0 = kfPos.GetEta();
      float const pseudorapV0Dau1 = kfNeg.GetEta();

      // info from from KFParticle
      std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()}; // pVec stands for vector containing the 3-momentum components
      std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      std::array<float, 3> pVecV0Dau0 = {kfPos.GetPx(), kfPos.GetPy(), kfPos.GetPz()};
      std::array<float, 3> pVecV0Dau1 = {kfNeg.GetPx(), kfNeg.GetPy(), kfNeg.GetPz()};

      //-------------------reconstruct cascade track------------------
      // pseudorapidity
      float const pseudorapCascBachelor = kfBachPionToXi.GetEta();

      // info from KFParticle
      std::array<float, 3> vertexCasc = {kfXi.GetX(), kfXi.GetY(), kfXi.GetZ()};
      std::array<float, 3> pVecCascBachelor = {kfBachPionToXi.GetPx(), kfBachPionToXi.GetPy(), kfBachPionToXi.GetPz()};

      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> pVecCharmBachelorAsD{};
      pVecCharmBachelorAsD[0] = kfCharmBachPionToXiC.GetPx();
      pVecCharmBachelorAsD[1] = kfCharmBachPionToXiC.GetPy();
      pVecCharmBachelorAsD[2] = kfCharmBachPionToXiC.GetPz();

      std::array<float, 3> pVecCharmBaryon = {kfXiC0.GetPx(), kfXiC0.GetPy(), kfXiC0.GetPz()};
      auto* covVtxCharmBaryon = kfXiC0.CovarianceMatrix();
      float covMatrixPV[6];
      kfVertex.GetCovarianceMatrix(covMatrixPV);

      // impact parameters
      std::array<float, 2> impactParameterV0Dau0{};
      std::array<float, 2> impactParameterV0Dau1{};
      std::array<float, 2> impactParameterPiFromCasc{};
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, xiDauChargedTrackParCov, 2.f, matCorr, &impactParameterPiFromCasc);
      float const dcaxyV0Dau0 = impactParameterV0Dau0[0];
      float const dcaxyV0Dau1 = impactParameterV0Dau1[0];
      float const dcaxyCascBachelor = impactParameterPiFromCasc[0];

      // pseudorapidity
      float const pseudorapCharmBachelor = kfCharmBachPionToXiC.GetEta();

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
      kfXic0Candidate.chi2NdfTopoV0ToPv = kfV0ToPv.GetChi2() / kfV0ToPv.GetNDF();
      kfXic0Candidate.chi2NdfTopoCascToPv = kfXiToPv.GetChi2() / kfXiToPv.GetNDF();
      kfXic0Candidate.chi2NdfTopoPiFromXicToPv = kfPiFromXicToPv.GetChi2() / kfPiFromXicToPv.GetNDF();
      kfXic0Candidate.chi2NdfTopoXicToPv = kfXic0ToPv.GetChi2() / kfXic0ToPv.GetNDF();

      auto cascBachTopoChi2 = kfBachPionToXi.GetChi2();
      kfXic0Candidate.chi2NdfTopoV0ToCasc = kfV0ToCasc.GetChi2() / kfV0ToCasc.GetNDF();
      kfXic0Candidate.chi2NdfTopoCascToXic = kfXiToXiC.GetChi2() / kfXiToXiC.GetNDF();

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

      // KF rapidity
      kfXic0Candidate.rapXic = kfXiC0.GetRapidity();

      // KF cosThetaStar
      kfXic0Candidate.cosThetaStarPiFromXic = cosThetaStarFromKF(0, 4132, 211, 3312, kfCharmBachPionToXiC, kfXiToXiC, pdgdb);

      // KF ct
      kfXic0Candidate.ctV0 = kfV0ToCasc.GetLifeTime();
      kfXic0Candidate.ctCasc = kfXiToXiC.GetLifeTime();
      kfXic0Candidate.ctXic = kfXic0ToPv.GetLifeTime();
      // KF eta
      kfXic0Candidate.etaXic = kfXiC0.GetEta();

      // fill KF hist
      registry.fill(HIST("hKFParticleCascBachTopoChi2"), cascBachTopoChi2);
      registry.fill(HIST("hKfXiC0_ldl"), kfXic0Candidate.ldlXic);

      // fill kf table
      kfCandidateXicData(collision.globalIndex(),
                         pvCoord[0], pvCoord[1], pvCoord[2],
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
                         v0index, casc.posTrackId(), casc.negTrackId(),
                         casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                         kfXic0Candidate.massV0, kfXic0Candidate.massCasc, kfXic0Candidate.massXic,
                         kfXic0Candidate.cosPaV0ToPv, kfXic0Candidate.cosPaCascToPv,
                         kfXic0Candidate.ctCasc, kfXic0Candidate.ctV0, kfXic0Candidate.ctXic,
                         pseudorapV0Dau0, pseudorapV0Dau1, pseudorapCascBachelor, pseudorapCharmBachelor,
                         kfXic0Candidate.etaXic, kfXi.GetEta(), kfV0.GetEta(),
                         dcaxyV0Dau0, dcaxyV0Dau1, dcaxyCascBachelor,
                         kfXic0Candidate.kfDcaCascDau, kfXic0Candidate.kfDcaV0Dau, kfXic0Candidate.kfDcaXicDau,
                         kfXic0Candidate.kfDcaXYPiFromXic, kfXic0Candidate.kfDcaXYCascToPv,
                         kfXic0Candidate.chi2GeoV0, kfXic0Candidate.chi2GeoCasc, kfXic0Candidate.chi2GeoXic, kfXic0Candidate.chi2MassV0, kfXic0Candidate.chi2MassCasc,
                         kfXic0Candidate.ldlV0, kfXic0Candidate.ldlCasc,
                         kfXic0Candidate.chi2NdfTopoV0ToPv, kfXic0Candidate.chi2NdfTopoCascToPv, kfXic0Candidate.chi2NdfTopoPiFromXicToPv, kfXic0Candidate.chi2NdfTopoXicToPv,
                         kfXic0Candidate.chi2NdfTopoV0ToCasc, kfXic0Candidate.chi2NdfTopoCascToXic,
                         kfXic0Candidate.decayLenXYLambda, kfXic0Candidate.decayLenXYCasc, kfXic0Candidate.decayLenXYXic,
                         kfXic0Candidate.cosPaV0ToCasc, kfXic0Candidate.cosPaCascToXic,
                         kfXic0Candidate.rapXic,
                         kfXic0Candidate.cosThetaStarPiFromXic,
                         v0NDF, cascNDF, charmbaryonNDF, v0Ndfm, cascNdfm,
                         v0Chi2OverNdf, cascChi2OverNdf, charmbaryonChi2OverNdf, v0Chi2OverNdfm, cascChi2OverNdfm);
      // fill QA table
      if (kfResolutionQA) {
        rowKfXic0Qa(massLam, massCasc, massXiC0, sigLam, sigCasc, sigXiC0,
                    collision.globalIndex(), v0index, casc.posTrackId(), casc.negTrackId(), casc.cascadeId(), trackCharmBachelor.globalIndex(), casc.bachelorId(),
                    kfPos.GetX(), kfPos.GetY(), kfPos.GetZ(), kfPos.GetErrX(), kfPos.GetErrY(), kfPos.GetErrZ(), kfPos.GetPt(),
                    kfNeg.GetX(), kfNeg.GetY(), kfNeg.GetZ(), kfNeg.GetErrX(), kfNeg.GetErrY(), kfNeg.GetErrZ(), kfNeg.GetPt(),
                    kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), kfV0.GetErrX(), kfV0.GetErrY(), kfV0.GetErrZ(), kfV0.GetPt(),
                    kfBachPionToXi.GetX(), kfBachPionToXi.GetY(), kfBachPionToXi.GetZ(), kfBachPionToXi.GetErrX(), kfBachPionToXi.GetErrY(), kfBachPionToXi.GetErrZ(), kfBachPionToXi.GetPt(),
                    kfXi.GetX(), kfXi.GetY(), kfXi.GetZ(), kfXi.GetErrX(), kfXi.GetErrY(), kfXi.GetErrZ(), kfXi.GetPt(),
                    kfCharmBachPionToXiC.GetX(), kfCharmBachPionToXiC.GetY(), kfCharmBachPionToXiC.GetZ(), kfCharmBachPionToXiC.GetErrX(), kfCharmBachPionToXiC.GetErrY(), kfCharmBachPionToXiC.GetErrZ(), kfCharmBachPionToXiC.GetPt(),
                    kfXiC0.GetX(), kfXiC0.GetY(), kfXiC0.GetZ(), kfXiC0.GetErrX(), kfXiC0.GetErrY(), kfXiC0.GetErrZ(), kfXiC0.GetPt(),
                    casc.xlambda(), casc.ylambda(), casc.zlambda(), casc.x(), casc.y(), casc.z());
      }
    } // loop over LF Cascade-bachelor candidates
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, int DecayChannel, typename Coll, typename Hist>
  void runOmegac0Xic0ToOmegaKaCreatorWithKFParticle(Coll const&,
                                                    aod::BCsWithTimestamps const& /*bcWithTimeStamps*/,
                                                    MyKfTracksIU const& tracksIU,
                                                    MyKfTracks const& tracks,
                                                    MyKfCascTable const&, KFCascadesLinked const&,
                                                    aod::HfCascLf2Prongs const& candidates,
                                                    Hist& hInvMassCharmBaryon,
                                                    Hist& hFitterStatus,
                                                    Hist& hCandidateCounter,
                                                    Hist& hCascadesCounter)
  {
    for (const auto& cand : candidates) {
      hCandidateCounter->Fill(1);

      //----------------------check if the event is selected-----------------------------
      auto collision = cand.collision_as<Coll>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      //----------------------Set the magnetic field from ccdb-----------------------------
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
        initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>>>> Magnetic field: " << magneticField;
        runNumber = bc.runNumber();
      }
      KFParticle::SetField(magneticField);

      // Retrieve skimmed cascade and pion tracks
      auto cascAodElement = cand.cascade_as<aod::KFCascadesLinked>();
      hCascadesCounter->Fill(0);
      if (!cascAodElement.has_kfCascData()) {
        continue;
      }
      auto casc = cascAodElement.kfCascData_as<MyKfCascTable>();
      hCascadesCounter->Fill(1);

      // convert KaonFromCharm&KaFromOmega&V0DauPos&V0DauNeg tracks into KFParticle object
      auto trackCharmBachelorId = cand.prong0Id();
      auto trackKaFromCharm = tracks.rawIteratorAt(trackCharmBachelorId);
      auto trackCascDauChargedId = casc.bachelorId();
      auto trackV0Dau0Id = casc.posTrackId();
      auto trackV0Dau1Id = casc.negTrackId();
      auto trackKaFromOmega = tracksIU.rawIteratorAt(trackCascDauChargedId); // Ka <- Omega track
      auto trackV0DauPos = tracksIU.rawIteratorAt(trackV0Dau0Id);            // V0 positive daughter track
      auto trackV0DauNeg = tracksIU.rawIteratorAt(trackV0Dau1Id);            // V0 negative daughter track
      auto kaFromOmegaCharge = trackKaFromOmega.sign();
      auto signOmega = casc.sign();

      KFPTrack const kfpTrackKaFromCharm = createKFPTrackFromTrack(trackKaFromCharm);
      KFPTrack const kfpTrackKaFromOmega = createKFPTrackFromTrack(trackKaFromOmega);
      KFPTrack const kfpTrackV0DauPos = createKFPTrackFromTrack(trackV0DauPos);
      KFPTrack const kfpTrackV0DauNeg = createKFPTrackFromTrack(trackV0DauNeg);

      KFParticle kfPrFromV0(kfpTrackV0DauPos, kProton);
      KFParticle kfPiFromV0(kfpTrackV0DauNeg, kPiMinus);
      KFParticle kfKaFromOmega(kfpTrackKaFromOmega, kKMinus);
      KFParticle kfPiFromXiRej(kfpTrackKaFromOmega, kPiMinus); // rej
      KFParticle kfKaFromCharm(kfpTrackKaFromCharm, kKPlus);

      if (signOmega == 0 || kaFromOmegaCharge == 0 || kaFromOmegaCharge != signOmega) {
        continue;
      }
      // convert for Pos and Neg Particles
      if (signOmega > 0) {
        kfPiFromV0 = KFParticle(kfpTrackV0DauPos, kPiPlus);
        kfPrFromV0 = KFParticle(kfpTrackV0DauNeg, -kProton);
        kfKaFromOmega = KFParticle(kfpTrackKaFromOmega, kKPlus);
        kfPiFromXiRej = KFParticle(kfpTrackKaFromOmega, kPiPlus); // rej
        kfKaFromCharm = KFParticle(kfpTrackKaFromCharm, kKMinus);
      }

      if (doCascadePreselection) {
        if (std::abs(casc.dcaXYCascToPV()) > dcaXYToPVCascadeMax) {
          continue;
        }
        if (std::abs(casc.dcaV0daughters()) > dcaV0DaughtersMax) {
          continue;
        }
        if (std::abs(casc.dcacascdaughters()) > dcaCascDaughtersMax) {
          continue;
        }
        if (std::abs(casc.mOmega() - MassOmegaMinus) > massToleranceCascade) {
          continue;
        }
      }

      //----------------------info of V0 and cascade tracks from LF-table------------------
      std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};

      // step 1 : construct V0 with KF
      const KFParticle* v0Daughters[2] = {&kfPrFromV0, &kfPiFromV0};
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
      if (std::abs(massLam - MassLambda0) > lambdaMassWindow) {
        continue;
      }
      // err_mass>0 of Lambda
      if (sigLam <= 0) {
        continue;
      }
      // chi2>0 && NDF>0 for selecting Lambda
      if ((kfV0.GetNDF() <= 0 || kfV0.GetChi2() <= 0)) {
        continue;
      }
      KFParticle kfV0MassConstrained = kfV0;
      kfV0MassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassLambda); // set mass constrain to Lambda
      if (kfUseV0MassConstraint) {
        kfV0 = kfV0MassConstrained;
      }
      kfV0.TransportToDecayVertex();

      // step 2 : reconstruct cascade(Omega) with KF
      const KFParticle* omegaDaugthers[2] = {&kfKaFromOmega, &kfV0};
      const KFParticle* omegaDaugthersRej[2] = {&kfPiFromXiRej, &kfV0}; // rej
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
      // err_massOmega and err_massXiRej > 0
      if (sigCasc <= 0 || sigCascrej <= 0) {
        continue;
      }
      // chi2>0 && NDF>0
      if (kfOmega.GetNDF() <= 0 || kfOmega.GetChi2() <= 0) {
        continue;
      }
      if ((std::abs(massCasc - MassOmegaMinus) > massToleranceCascade) || (std::abs(massCascrej - MassXiMinus) < massToleranceCascadeRej)) {
        continue;
      }
      registry.fill(HIST("hInvMassXiMinus_rej"), massCascrej); // rej: Add competing rejection to minimize misidentified Xi impact. Reject if kfBachPionRej is Pion and the constructed cascade has Xi's invariant mass.
      KFParticle kfOmegaMassConstrained = kfOmega;
      kfOmegaMassConstrained.SetNonlinearMassConstraint(o2::constants::physics::MassOmegaMinus); // set mass constrain to XiMinus
      if (kfUseCascadeMassConstraint) {
        // set mass constraint if requested
        kfOmega = kfOmegaMassConstrained;
      }
      registry.fill(HIST("hInvMassXiMinus"), massCasc);
      kfOmega.TransportToDecayVertex();
      // rej: Add competing rejection to minimize misidentified Xi impact. Reject if kfBachPionRej is Pion and the constructed cascade has Xi's invariant mass.

      // step 3 : reconstruc OmegaKa with KF
      //  Create KF charm bach Pion from track
      const KFParticle* omegaKaDaugthers[2] = {&kfKaFromCharm, &kfOmega};
      // construct Omegac0 or Xic0
      KFParticle kfOmegaKa;
      kfOmegaKa.SetConstructMethod(kfConstructMethod);
      try {
        kfOmegaKa.Construct(omegaKaDaugthers, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct OmegaKa from Cascade and bachelor pion track: " << e.what();
        continue;
      }
      float massOmegaKa, sigOmegaKa;
      kfOmegaKa.GetMass(massOmegaKa, sigOmegaKa);
      if (sigOmegaKa <= 0) {
        continue;
      }
      if (kfOmegaKa.GetNDF() <= 0 || kfOmegaKa.GetChi2() <= 0) {
        continue;
      }
      kfOmegaKa.TransportToDecayVertex();
      hFitterStatus->Fill(0);
      hCandidateCounter->Fill(2);

      // initialize primary vertex
      KFPVertex const kfpVertex = createKFPVertexFromCollision(collision);
      float covMatrixPV[6];
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle const kfPv(kfpVertex); // for calculation of DCAs to PV

      // fill test histograms
      hInvMassCharmBaryon->Fill(massOmegaKa);

      // topological constraint of daughter to mother
      KFParticle kfKaFromCharmToOmegaKa = kfKaFromCharm;
      KFParticle kfOmegaToOmegaKa = kfOmega;
      KFParticle kfV0ToOmega = kfV0;
      KFParticle kfKaToOmega = kfKaFromOmega;
      KFParticle kfPrToV0 = kfPrFromV0;
      KFParticle kfPiToV0 = kfPiFromV0;

      kfPrToV0.SetProductionVertex(kfV0);
      kfPiToV0.SetProductionVertex(kfV0);
      kfV0ToOmega.SetProductionVertex(kfOmega);
      kfKaToOmega.SetProductionVertex(kfOmega);
      kfKaFromCharmToOmegaKa.SetProductionVertex(kfOmegaKa);
      kfOmegaToOmegaKa.SetProductionVertex(kfOmegaKa);

      // topological constraint to PV
      // KFParticle to PV
      KFParticle kfV0ToPv = kfV0;
      KFParticle kfOmegaToPv = kfOmega;
      KFParticle kfCharmToPv = kfOmegaKa;
      KFParticle kfKaFromCharmToPv = kfKaFromCharm;

      kfV0ToPv.SetProductionVertex(kfPv);
      kfOmegaToPv.SetProductionVertex(kfPv);
      kfCharmToPv.SetProductionVertex(kfPv);
      kfKaFromCharmToPv.SetProductionVertex(kfPv);

      //---------------------calculate physical parameters of OmegaKa candidate----------------------

      // transport OmegaKa daughters to decay vertex (secondary vertex)
      float secondaryVertex[3] = {0.};
      secondaryVertex[0] = kfOmegaKa.GetX();
      secondaryVertex[1] = kfOmegaKa.GetY();
      secondaryVertex[2] = kfOmegaKa.GetZ();
      kfKaFromCharm.TransportToPoint(secondaryVertex);
      kfOmega.TransportToPoint(secondaryVertex);

      // get impact parameters of OmegaKa daughters
      float impactParameterKaFromCharmXY = 0., errImpactParameterKaFromCharmXY = 0.;
      float impactParameterOmegaXY = 0., errImpactParameterOmegaXY = 0.;
      kfKaFromCharm.GetDistanceFromVertexXY(kfPv, impactParameterKaFromCharmXY, errImpactParameterKaFromCharmXY);
      kfOmega.GetDistanceFromVertexXY(kfPv, impactParameterOmegaXY, errImpactParameterOmegaXY);

      // calculate cosine of pointing angle
      float const cosPaV0ToPv = cpaFromKF(kfV0, kfPv);
      float const cosPaCascToPv = cpaFromKF(kfOmega, kfPv);
      float const cosPaOmegaKaToPv = cpaFromKF(kfOmegaKa, kfPv);
      float const cosPaXYV0ToPv = cpaXYFromKF(kfV0, kfPv);
      float const cosPaXYCascToPv = cpaXYFromKF(kfOmega, kfPv);
      float const cosPaXYOmegaKaToPv = cpaXYFromKF(kfOmegaKa, kfPv);
      float const cosPaV0ToCasc = cpaFromKF(kfV0, kfOmega);
      float const cosPaCascToOmegaKa = cpaFromKF(kfOmega, kfOmegaKa);
      float const cosPaXYV0ToCasc = cpaXYFromKF(kfV0, kfOmega);
      float const cosPaXYCascToOmegaKa = cpaXYFromKF(kfOmega, kfOmegaKa);

      // Get Chi2Geo/NDF
      float const chi2GeoV0 = kfV0.GetChi2() / kfV0.GetNDF();
      float const chi2GeoCasc = kfOmega.GetChi2() / kfOmega.GetNDF();
      float const chi2GeoOmegaKa = kfOmegaKa.GetChi2() / kfOmegaKa.GetNDF();

      // Get Chi2Topo/NDF
      float const chi2NdfTopoV0ToCasc = kfV0ToOmega.GetChi2() / kfV0ToOmega.GetNDF();
      float const chi2NdfTopoKaToCasc = kfKaToOmega.GetChi2() / kfKaToOmega.GetNDF();
      float const chi2NdfTopoKaFromOmegaKaToOmegaKa = kfKaFromCharmToOmegaKa.GetChi2() / kfKaFromCharmToOmegaKa.GetNDF();
      float const chi2NdfTopoCascToOmegaKa = kfOmegaToOmegaKa.GetChi2() / kfOmegaToOmegaKa.GetNDF();
      float const chi2NdfTopoV0ToPv = kfV0ToPv.GetChi2() / kfV0ToPv.GetNDF();
      float const chi2NdfTopoCascToPv = kfOmegaToPv.GetChi2() / kfOmegaToPv.GetNDF();
      float const chi2NdfTopoOmegaKaToPv = kfCharmToPv.GetChi2() / kfCharmToPv.GetNDF();
      float const chi2NdfTopoKaFromOmegaKaToPv = kfKaFromCharmToPv.GetChi2() / kfKaFromCharmToPv.GetNDF();

      // Get MassChi2/NDF
      auto v0Chi2OverNdfm = kfV0MassConstrained.GetChi2() / kfV0MassConstrained.GetNDF();
      auto cascChi2OverNdfm = kfOmegaMassConstrained.GetChi2() / kfOmegaMassConstrained.GetNDF();

      // KF ldl
      float const ldlV0 = ldlFromKF(kfV0, kfPv);
      float const ldlCasc = ldlFromKF(kfOmega, kfPv);
      float const ldlOmegaKa = ldlFromKF(kfOmegaKa, kfPv);

      // KF decay length
      float decayLxyLam, errDecayLxyLam;
      kfV0ToOmega.GetDecayLengthXY(decayLxyLam, errDecayLxyLam);
      float decayLxyCasc, errDecayLxyCasc;
      kfOmegaToOmegaKa.GetDecayLengthXY(decayLxyCasc, errDecayLxyCasc);
      float decayLxyOmegaKa, errDecayLxyOmegaKa;
      kfCharmToPv.GetDecayLengthXY(decayLxyOmegaKa, errDecayLxyOmegaKa);

      // KF pT
      float const ptOmegaKa = kfOmegaKa.GetPt();
      float const ptKaFromCharm = kfKaFromCharm.GetPt();
      float const ptOmega = kfOmega.GetPt();

      // KF cosThetaStar
      float const cosThetaStarKaFromOmegac = cosThetaStarFromKF(0, 4332, 321, 3334, kfKaFromCharmToOmegaKa, kfOmegaToOmegaKa, pdgdb);
      float const cosThetaStarKaFromXic = cosThetaStarFromKF(0, 4132, 321, 3334, kfKaFromCharmToOmegaKa, kfOmegaToOmegaKa, pdgdb);

      // KF ct
      float const ctV0 = kfV0ToOmega.GetLifeTime();
      float const ctCasc = kfOmegaToOmegaKa.GetLifeTime();
      float const ctOmegaKa = kfCharmToPv.GetLifeTime();

      hCandidateCounter->Fill(3);

      // fill full kf table
      kfCandidateOmegaKaData(collision.globalIndex(),
                             collision.posX(), collision.posY(), collision.posZ(),                                                                                                                                                // PV Coord
                             kfPv.GetX(), kfPv.GetY(), kfPv.GetZ(),                                                                                                                                                               // PV KF
                             vertexV0[0], vertexV0[1], vertexV0[2],                                                                                                                                                               // V0 Vtx from LF-table
                             pVecV0[0], pVecV0[1], pVecV0[2],                                                                                                                                                                     // V0 P from LF-table
                             vertexCasc[0], vertexCasc[1], vertexCasc[2],                                                                                                                                                         // Casc Vtx from LF-table
                             pVecCasc[0], pVecCasc[1], pVecCasc[2],                                                                                                                                                               // Casc P from LF-table
                             kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(),                                                                                                                                                               // V0 Vtx KF
                             kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz(),                                                                                                                                                            // V0 P KF
                             kfOmega.GetX(), kfOmega.GetY(), kfOmega.GetZ(),                                                                                                                                                      // Omega Vtx KF
                             kfOmega.GetPx(), kfOmega.GetPx(), kfOmega.GetPx(),                                                                                                                                                   // Omega Vtx KF
                             kfOmegaKa.GetX(), kfOmegaKa.GetY(), kfOmegaKa.GetZ(),                                                                                                                                                // OmegaKa Vtx KF (SecondaryVertex)
                             kfOmegaKa.GetPx(), kfOmegaKa.GetPx(), kfOmegaKa.GetPx(),                                                                                                                                             // OmegaKa P KF
                             signOmega,                                                                                                                                                                                           // Check Omega sign
                             kfPrFromV0.GetEta(), kfPiFromV0.GetEta(), kfKaFromOmega.GetEta(), kfKaFromCharm.GetEta(), kfV0.GetEta(), kfOmega.GetEta(), kfOmegaKa.GetEta(), kfOmegaKa.GetRapidity(),                              // Eta of daughters and mothers. Rapidity of OmegaKa
                             impactParameterKaFromCharmXY, errImpactParameterKaFromCharmXY, impactParameterOmegaXY, errImpactParameterOmegaXY,                                                                                    // DCAXY of KaFromCharm and Omega
                             kfPrToV0.GetDistanceFromParticle(kfPiToV0), kfV0ToOmega.GetDistanceFromParticle(kfKaToOmega), kfOmegaToOmegaKa.GetDistanceFromParticle(kfKaFromCharmToOmegaKa),                                      // DCA of daughters
                             cosPaV0ToPv, cosPaCascToPv, cosPaOmegaKaToPv, cosPaXYV0ToPv, cosPaXYCascToPv, cosPaXYOmegaKaToPv, cosPaV0ToCasc, cosPaCascToOmegaKa, cosPaXYV0ToCasc, cosPaXYCascToOmegaKa,                          // CosPA of PV and mothers
                             chi2GeoV0, chi2GeoCasc, chi2GeoOmegaKa,                                                                                                                                                              // Chi2Geo/NDF
                             v0Chi2OverNdfm, cascChi2OverNdfm,                                                                                                                                                                    // Chi2Mass/NDF
                             chi2NdfTopoV0ToCasc, chi2NdfTopoKaToCasc, chi2NdfTopoKaFromOmegaKaToOmegaKa, chi2NdfTopoCascToOmegaKa, chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoKaFromOmegaKaToPv, chi2NdfTopoOmegaKaToPv, // Chi2Topo/NDF
                             ldlV0, ldlCasc, ldlOmegaKa,                                                                                                                                                                          // ldl
                             decayLxyLam, decayLxyCasc, decayLxyOmegaKa,                                                                                                                                                          // DecaylengthXY
                             massLam, sigLam, massCasc, sigCasc, massCascrej, sigCascrej, massOmegaKa, sigOmegaKa,                                                                                                                // massKF and masserror
                             ptOmegaKa, ptKaFromCharm, ptOmega,                                                                                                                                                                   // pT
                             cosThetaStarKaFromOmegac, cosThetaStarKaFromXic, ctV0, ctCasc, ctOmegaKa,                                                                                                                            // cosThetaStar & ct
                             cascAodElement.v0Id(), casc.posTrackId(), casc.negTrackId(), casc.cascadeId(), casc.bachelorId(), trackKaFromCharm.globalIndex());
    }
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

  void processNoCentOmegacToOmegaPiWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                  aod::BCsWithTimestamps const& bcWithTimeStamps,
                                                  MyKfTracksIU const& tracksIU,
                                                  MyKfTracks const& tracks,
                                                  MyKfCascTable const& cascades,
                                                  KFCascadesLinked const& cascadeLinks,
                                                  aod::HfCascLf2Prongs const& candidates)
  {
    runKfOmegac0CreatorWithKFParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentOmegacToOmegaPiWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega pi decay channel using KFParticle", false);

  void processNoCentOmegac0Xic0ToOmegaKaCreatorWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                                              aod::BCsWithTimestamps const& bcWithTimeStamps,
                                                              MyKfTracksIU const& tracksIU,
                                                              MyKfTracks const& tracks,
                                                              MyKfCascTable const& cascades,
                                                              KFCascadesLinked const& cascadeLinks,
                                                              aod::HfCascLf2Prongs const& candidates)
  {
    runOmegac0Xic0ToOmegaKaCreatorWithKFParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentOmegac0Xic0ToOmegaKaCreatorWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega ka decay channel using KFParticle", false);

  void processNoCentXicToXiPiWithKFParticle(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                                            aod::BCsWithTimestamps const& bcWithTimeStamps,
                                            MyKfTracksIU const& tracksIU,
                                            MyKfTracks const& tracks,
                                            MyKfCascTable const& cascades,
                                            KFCascadesLinked const& cascadeLinks,
                                            aod::HfCascLf2Prongs const& candidates)
  {
    runKfXic0CreatorWithKFParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processNoCentXicToXiPiWithKFParticle, "Run candidate creator w/o centrality selections for Xic0 To Xi pi decay channel using KFParticle", false);

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

  void processCentFT0COmegacToOmegaPiWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                                    aod::BCsWithTimestamps const& bcWithTimeStamps,
                                                    MyKfTracksIU const& tracksIU,
                                                    MyKfTracks const& tracks,
                                                    MyKfCascTable const& cascades,
                                                    KFCascadesLinked const& cascadeLinks,
                                                    aod::HfCascLf2Prongs const& candidates)
  {
    runKfOmegac0CreatorWithKFParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0COmegacToOmegaPiWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega pi decay channel using KFParticle", false);

  void processCentFT0COmegac0Xic0ToOmegaKaCreatorWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                                                MyKfTracksIU const& tracksIU,
                                                                MyKfTracks const& tracks,
                                                                MyKfCascTable const& cascades,
                                                                KFCascadesLinked const& cascadeLinks,
                                                                aod::HfCascLf2Prongs const& candidates)
  {
    runOmegac0Xic0ToOmegaKaCreatorWithKFParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0COmegac0Xic0ToOmegaKaCreatorWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega ka decay channel using KFParticle", false);

  void processCentFT0CXicToXiPiWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps,
                                              MyKfTracksIU const& tracksIU,
                                              MyKfTracks const& tracks,
                                              MyKfCascTable const& cascades,
                                              KFCascadesLinked const& cascadeLinks,
                                              aod::HfCascLf2Prongs const& candidates)
  {
    runKfXic0CreatorWithKFParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0CXicToXiPiWithKFParticle, "Run candidate creator w FT0C centrality selections for Xic0 To Xi pi decay channel using KFParticle", false);

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

  void processCentFT0MOmegacToOmegaPiWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                                    aod::BCsWithTimestamps const& bcWithTimeStamps,
                                                    MyKfTracksIU const& tracksIU,
                                                    MyKfTracks const& tracks,
                                                    MyKfCascTable const& cascades,
                                                    KFCascadesLinked const& cascadeLinks,
                                                    aod::HfCascLf2Prongs const& candidates)
  {
    runKfOmegac0CreatorWithKFParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaPi, hFitterStatusToOmegaPi, hCandidateCounterToOmegaPi, hCascadesCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MOmegacToOmegaPiWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega pi decay channel using KFParticle", false);

  void processCentFT0MOmegac0Xic0ToOmegaKaCreatorWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                                                aod::BCsWithTimestamps const& bcWithTimeStamps,
                                                                MyKfTracksIU const& tracksIU,
                                                                MyKfTracks const& tracks,
                                                                MyKfCascTable const& cascades,
                                                                KFCascadesLinked const& cascadeLinks,
                                                                aod::HfCascLf2Prongs const& candidates)
  {
    runOmegac0Xic0ToOmegaKaCreatorWithKFParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToOmegaK, hFitterStatusToOmegaK, hCandidateCounterToOmegaK, hCascadesCounterToOmegaK);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MOmegac0Xic0ToOmegaKaCreatorWithKFParticle, "Run candidate creator w/o centrality selections for Omegac0 To omega ka decay channel using KFParticle", false);

  void processCentFT0MXicToXiPiWithKFParticle(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::BCsWithTimestamps const& bcWithTimeStamps,
                                              MyKfTracksIU const& tracksIU,
                                              MyKfTracks const& tracks,
                                              MyKfCascTable const& cascades,
                                              KFCascadesLinked const& cascadeLinks,
                                              aod::HfCascLf2Prongs const& candidates)
  {
    runKfXic0CreatorWithKFParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, bcWithTimeStamps, tracksIU, tracks, cascades, cascadeLinks, candidates, hInvMassCharmBaryonToXiPi, hFitterStatusToXiPi, hCandidateCounterToXiPi, hCascadesCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0, processCentFT0MXicToXiPiWithKFParticle, "Run candidate creator w FT0M centrality selections for Xic0 To Xi pi decay channel using KFParticle", false);

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
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

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
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

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
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, hfEvSel.occEstimator);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      const auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      const auto ir = hfEvSel.getInteractionRate(bc, ccdb); // Hz
      /// monitor the satisfied event selections
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy, ir);

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
  Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"};
  Configurable<bool> acceptTrackIntWithMaterial{"acceptTrackIntWithMaterial", false, " switch to accept candidates with final (i.e. p, K, pi) daughter tracks interacting with material"};

  using MyTracksWMc = soa::Join<TracksIU, McTrackLabels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring

  std::shared_ptr<TH1> hGenCharmBaryonPtRapidityTightXicToXiPi, hGenCharmBaryonPtRapidityLooseXicToXiPi, hGenCharmBaryonPtRapidityTightOmegacToXiPi, hGenCharmBaryonPtRapidityLooseOmegacToXiPi, hGenCharmBaryonPtRapidityTightOmegacToOmegaPi, hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi, hGenCharmBaryonPtRapidityTightOmegacToOmegaK, hGenCharmBaryonPtRapidityLooseOmegacToOmegaK;

  HistogramRegistry registry{"registry"};

  // inspect for which zPvPosMax cut was set for reconstructed
  void init(InitContext& initContext)
  {
    std::array<bool, 5> procCollisionsXicToXiPi{doprocessMcXicToXiPi, doprocessMcXicToXiPiFT0m, doprocessMcXicToXiPiFT0c, doprocessMcXicToXiPiKf, doprocessMcXicToXiPiKfQa};
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
      if (device.name == "hf-candidate-creator-xic0-omegac0") {
        // init HF event selection helper
        hfEvSelMc.init(device, registry);
        break;
      }
    }

    hGenCharmBaryonPtRapidityTightXicToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityTightXicToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}}); // keep track of generated candidates pt when |y|<0.5
    hGenCharmBaryonPtRapidityLooseXicToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseXicToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}}); // keep track of generated candidates pt when |y|<0.8

    hGenCharmBaryonPtRapidityTightOmegacToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityTightOmegacToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    hGenCharmBaryonPtRapidityLooseOmegacToXiPi = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseOmegacToXiPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});

    hGenCharmBaryonPtRapidityTightOmegacToOmegaPi = registry.add<TH1>("hGenCharmBaryonPtRapidityTightOmegacToOmegaPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});

    hGenCharmBaryonPtRapidityTightOmegacToOmegaK = registry.add<TH1>("hGenCharmBaryonPtRapidityTightOmegacToOmegaK", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    hGenCharmBaryonPtRapidityLooseOmegacToOmegaK = registry.add<TH1>("hGenCharmBaryonPtRapidityLooseOmegacToOmegaK", "Generated charm baryon #it{p}_{T};#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});

    // QA
    if (doprocessMcXicToXiPiKfQa) {
      AxisSpec const axisPt{20, 0., 20.};
      AxisSpec const axisDelta{1000, -0.5, 0.5};
      AxisSpec const axisPull{2000, -10., 10.};
      AxisSpec const axisPtRes{400, -0.2, 0.2};
      // mass over pt
      registry.add("hV0MassPullVsPt", "m_{PULL}(V0) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXiMassPullVsPt", "m_{PULL}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXic0MassPullVsPt", "m_{PULL}(#Xic0) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      // delta
      registry.add("hV0DauPosXDelta", "x^{p} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hV0DauPosYDelta", "y^{p} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hV0DauPosZDelta", "z^{p} - z^{MC}", kTH1D, {axisDelta});
      registry.add("hV0DauNegXDelta", "x^{#pi^{-}} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hV0DauNegYDelta", "y^{#pi^{-}} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hV0DauNegZDelta", "z^{#pi^{-}} - z^{MC}", kTH1D, {axisDelta});
      registry.add("hV0XDelta", "x^{#Lambda} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hV0YDelta", "y^{#Lambda} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hV0ZDelta", "z^{#Lambda} - z^{MC}", kTH1D, {axisDelta});

      registry.add("hXiBachelorXDelta", "x^{#pi^{-} from #Xi^{-}} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hXiBachelorYDelta", "y^{#pi^{-} from #Xi^{-}} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hXiBachelorZDelta", "z^{#pi^{-} from #Xi^{-}} - z^{MC}", kTH1D, {axisDelta});

      registry.add("hXiXDelta", "x^{#Xi^{-}} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hXiYDelta", "y^{#Xi^{-}} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hXiZDelta", "z^{#Xi^{-}} - z^{MC}", kTH1D, {axisDelta});

      registry.add("hXic0BachelorXDelta", "x^{#pi^{+} from #Xi_{c}^{0}} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hXic0BachelorYDelta", "y^{#pi^{+} from #Xi_{c}^{0}} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hXic0BachelorZDelta", "z^{#pi^{+} from #Xi_{c}^{0}} - z^{MC}", kTH1D, {axisDelta});

      registry.add("hXic0XDelta", "x^{#Xi_(c)^(0)} - x^{MC}", kTH1D, {axisDelta});
      registry.add("hXic0YDelta", "y^{#Xi_(c)^(0)} - y^{MC}", kTH1D, {axisDelta});
      registry.add("hXic0ZDelta", "z^{#Xi_(c)^(0)} - z^{MC}", kTH1D, {axisDelta});
      // delta over pt
      registry.add("hV0DauPosXDeltaVsPt", "#Delta_{x}(p) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0DauPosYDeltaVsPt", "#Delta_{y}(p) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0DauPosZDeltaVsPt", "#Delta_{z}(p) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0DauNegXDeltaVsPt", "#Delta_{x}(#pi) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0DauNegYDeltaVsPt", "#Delta_{y}(#pi) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0DauNegZDeltaVsPt", "#Delta_{z}(#pi) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0XDeltaVsPt", "#Delta_{x}(#Lambda) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0YDeltaVsPt", "#Delta_{y}(#Lambda) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hV0ZDeltaVsPt", "#Delta_{z}(#Lambda) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});

      registry.add("hXiBachelorXDeltaVsPt", "#Delta_{x}(#pi^{-} from #Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXiBachelorYDeltaVsPt", "#Delta_{y}(#pi^{-} from #Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXiBachelorZDeltaVsPt", "#Delta_{z}(#pi^{-} from #Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});

      registry.add("hXiXDeltaVsPt", "#Delta_{x}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXiYDeltaVsPt", "#Delta_{y}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXiZDeltaVsPt", "#Delta_{z}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});

      registry.add("hXic0BachelorXDeltaVsPt", "#Delta_{x}(#pi^{+} from #Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXic0BachelorYDeltaVsPt", "#Delta_{y}(#pi^{+} from #Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXic0BachelorZDeltaVsPt", "#Delta_{z}(#pi^{+} from #Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});

      registry.add("hXic0XDeltaVsPt", "#Delta_{x}(#Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXic0YDeltaVsPt", "#Delta_{y}(#Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});
      registry.add("hXic0ZDeltaVsPt", "#Delta_{z}(#Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisDelta});

      // pull
      registry.add("hV0DauPosXPull", "x^{PULL}", kTH1D, {axisPull});
      registry.add("hV0DauPosYPull", "y^{PULL}", kTH1D, {axisPull});
      registry.add("hV0DauPosZPull", "z^{PULL}", kTH1D, {axisPull});
      registry.add("hV0DauNegXPull", "x^{PULL}", kTH1D, {axisPull});
      registry.add("hV0DauNegYPull", "y^{PULL}", kTH1D, {axisPull});
      registry.add("hV0DauNegZPull", "z^{PULL}", kTH1D, {axisPull});
      registry.add("hV0XPull", "x^{PULL}(#Lambda)", kTH1D, {axisPull});
      registry.add("hV0YPull", "y^{PULL}(#Lambda)", kTH1D, {axisPull});
      registry.add("hV0ZPull", "z^{PULL}(#Lambda)", kTH1D, {axisPull});

      registry.add("hXiBachelorXPull", "x^{PULL}", kTH1D, {axisPull});
      registry.add("hXiBachelorYPull", "y^{PULL}", kTH1D, {axisPull});
      registry.add("hXiBachelorZPull", "z^{PULL}", kTH1D, {axisPull});

      registry.add("hXiXPull", "x^{PULL}(#Xi^{-})", kTH1D, {axisPull});
      registry.add("hXiYPull", "y^{PULL}(#Xi^{-})", kTH1D, {axisPull});
      registry.add("hXiZPull", "z^{PULL}(#Xi^{-})", kTH1D, {axisPull});

      registry.add("hXic0BachelorXPull", "x^{PULL}", kTH1D, {axisPull});
      registry.add("hXic0BachelorYPull", "y^{PULL}", kTH1D, {axisPull});
      registry.add("hXic0BachelorZPull", "z^{PULL}", kTH1D, {axisPull});

      registry.add("hXic0XPull", "x^{PULL}(#Xi_{c}^{0})", kTH1D, {axisPull});
      registry.add("hXic0YPull", "y^{PULL}(#Xi_{c}^{0})", kTH1D, {axisPull});
      registry.add("hXic0ZPull", "z^{PULL}(#Xi_{c}^{0})", kTH1D, {axisPull});
      // pull over pt
      registry.add("hV0DauPosXPullVsPt", "x_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0DauPosYPullVsPt", "y_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0DauPosZPullVsPt", "z_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0DauNegXPullVsPt", "x_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0DauNegYPullVsPt", "y_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0DauNegZPullVsPt", "z_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0XPullVsPt", "x_{PULL}(#Lambda) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0YPullVsPt", "y_{PULL}(#Lambda) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hV0ZPullVsPt", "z_{PULL}(#Lambda) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});

      registry.add("hXiBachelorXPullVsPt", "x_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXiBachelorYPullVsPt", "y_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXiBachelorZPullVsPt", "z_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});

      registry.add("hXiXPullVsPt", "x_{PULL}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXiYPullVsPt", "y_{PULL}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXiZPullVsPt", "z_{PULL}(#Xi^{-}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});

      registry.add("hXic0BachelorXPullVsPt", "x_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXic0BachelorYPullVsPt", "y_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXic0BachelorZPullVsPt", "z_{PULL} vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});

      registry.add("hXic0XPullVsPt", "x_{PULL}(#Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXic0YPullVsPt", "y_{PULL}(#Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});
      registry.add("hXic0ZPullVsPt", "z_{PULL}(#Xi_{c}^{0}) vs. p_{T}", HistType::kTH2D, {axisPt, axisPull});

      // Defaut delta
      registry.add("hLambdaXDelta", "x^{#Lambda} - x^{MC}(Default)", kTH1D, {axisDelta});
      registry.add("hLambdaYDelta", "y^{#Lambda} - y^{MC}(Default)", kTH1D, {axisDelta});
      registry.add("hLambdaZDelta", "z^{#Lambda} - z^{MC}(Default)", kTH1D, {axisDelta});

      registry.add("hCascXDelta", "x^{#Xi^{-}} - x^{MC}(Default)", kTH1D, {axisDelta});
      registry.add("hCascYDelta", "y^{#Xi^{-}} - y^{MC}(Default)", kTH1D, {axisDelta});
      registry.add("hCascZDelta", "z^{#Xi^{-}} - z^{MC}(Default)", kTH1D, {axisDelta});

      // Pt Resolution
      registry.add("hV0DauPosPtRes", "Pt Resolution (p)", kTH1D, {axisPtRes});
      registry.add("hV0DauNegPtRes", "Pt Resolution (#pi^{-} from #Lambda)", kTH1D, {axisPtRes});
      registry.add("hV0PtRes", "Pt Resolution (V0)", kTH1D, {axisPtRes});
      registry.add("hXiBachelorPtRes", "Pt Resolution (#pi^{-} from #Xi^{-})", kTH1D, {axisPtRes});
      registry.add("hXiPtRes", "Pt Resolution (#Xi^{-})", kTH1D, {axisPtRes});
      registry.add("hXic0BachelorPtRes", "Pt Resolution (#pi^{+} from #Xi_{c}^{0})", kTH1D, {axisPtRes});
      registry.add("hXic0PtRes", "Pt Resolution (#Xi_{c}^{0})", kTH1D, {axisPtRes});
    }
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, int DecayChannel, typename Colls, typename TMyRecoCand, typename McCollisions>
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
    McMatchFlag debug{McMatchFlag::None};
    int8_t debugGenCharmBar = 0;
    int8_t debugGenCasc = 0;
    int8_t debugGenLambda = 0;
    int8_t nPiToMuV0{0}, nPiToMuCasc{0}, nPiToMuOmegac0{0};
    int8_t nKaToPiCasc{0}, nKaToPiOmegac0{0};
    bool collisionMatched = false;

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = McMatchFlag::None;
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
      if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
        // Xic → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, +kXiC0, std::array{+kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = McMatchFlag::CharmbaryonUnmatched;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = McMatchFlag::CascUnmatched;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = McMatchFlag::V0Unmatched;
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
        if (debug == McMatchFlag::CascUnmatched || debug == McMatchFlag::V0Unmatched) {
          LOGF(info, "WARNING: Xic0ToXiPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) { // Omegac -> xi pi matching
        // Omegac → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, +kOmegaC0, std::array{+kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = McMatchFlag::CharmbaryonUnmatched;
        }
        if (indexRec > -1) {
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = McMatchFlag::CascUnmatched;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = McMatchFlag::V0Unmatched;
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
        if (debug == McMatchFlag::CascUnmatched || debug == McMatchFlag::V0Unmatched) {
          LOGF(info, "WARNING: Omegac0ToXiPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) { // Omegac0 -> omega pi matching
        if (acceptTrackIntWithMaterial) {
          // Omegac → pi K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, arrayDaughters, +kOmegaC0, std::array{+kPiPlus, +kKMinus, +kProton, +kPiMinus}, true, &sign, 3, &nPiToMuOmegac0, &nKaToPiOmegac0);
          indexRecCharmBaryon = indexRec;
          if (indexRec == -1) {
            debug = McMatchFlag::CharmbaryonUnmatched;
          }
          if (indexRec > -1) {
            // Omega- → K pi p
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, arrayDaughtersCasc, +kOmegaMinus, std::array{+kKMinus, +kProton, +kPiMinus}, true, &signCasc, 2, &nPiToMuCasc, &nKaToPiCasc);
            if (indexRec == -1) {
              debug = McMatchFlag::CascUnmatched;
            }
            if (indexRec > -1) {
              // Lambda → p pi
              indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1, &nPiToMuV0);
              if (indexRec == -1) {
                debug = McMatchFlag::V0Unmatched;
              }
              if (indexRec > -1 && nPiToMuOmegac0 >= 1 && nKaToPiOmegac0 == 0) {
                flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPiOneMu);
                collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
              } else if (indexRec > -1 && nPiToMuOmegac0 == 0 && nKaToPiOmegac0 == 0) {
                flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi);
                collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
              }
            }
          }
        } else {
          // Omegac → pi K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, arrayDaughters, +kOmegaC0, std::array{+kPiPlus, +kKMinus, +kProton, +kPiMinus}, true, &sign, 3, &nPiToMuOmegac0, &nKaToPiOmegac0);
          indexRecCharmBaryon = indexRec;
          if (indexRec == -1) {
            debug = McMatchFlag::CharmbaryonUnmatched;
          }
          if (indexRec > -1) {
            // Omega- → K pi p
            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, arrayDaughtersCasc, +kOmegaMinus, std::array{+kKMinus, +kProton, +kPiMinus}, true, &signCasc, 2, &nPiToMuCasc, &nKaToPiCasc);
            if (indexRec == -1) {
              debug = McMatchFlag::CascUnmatched;
            }
            if (indexRec > -1) {
              // Lambda → p pi
              indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1, &nPiToMuV0);
              if (indexRec == -1) {
                debug = McMatchFlag::V0Unmatched;
              }
              if (indexRec > -1 && nPiToMuOmegac0 >= 1 && nKaToPiOmegac0 == 0) {
                flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPiOneMu);
                collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
              } else if (indexRec > -1 && nPiToMuOmegac0 == 0 && nKaToPiOmegac0 == 0) {
                flag = sign * (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi);
                collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
              }
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
        if (debug == McMatchFlag::CascUnmatched || debug == McMatchFlag::V0Unmatched) {
          LOGF(info, "WARNING: Omegac0ToOmegaPi decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) { // Omegac0 -> omega K matching
        // Omegac → K K pi p
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, +kOmegaC0, std::array{+kKPlus, +kKMinus, +kProton, +kPiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = McMatchFlag::CharmbaryonUnmatched;
        }
        if (indexRec > -1) {
          // Omega- → K pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, +kOmegaMinus, std::array{+kKMinus, +kProton, +kPiMinus}, true, &signCasc, 2);
          if (indexRec == -1) {
            debug = McMatchFlag::CascUnmatched;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);
            if (indexRec == -1) {
              debug = McMatchFlag::V0Unmatched;
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
        if (debug == McMatchFlag::CascUnmatched || debug == McMatchFlag::V0Unmatched) {
          LOGF(info, "WARNING: Omegac0ToOmegaK decays in the expected final state but the condition on the intermediate states are not fulfilled");
        }
      }
    } // close loop over candidates

    for (const auto& mcCollision : mcCollisions) {

      // Slice the particles table to get the particles for the current MC collision
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      // Slice the collisions table to get the collision info for the current MC collision
      float centrality{-1.f};
      o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
      int nSplitColl = 0;
      if constexpr (CentEstimator == CentralityEstimator::FT0C) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
        nSplitColl = collSlice.size();
      } else if constexpr (CentEstimator == CentralityEstimator::FT0M) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
        nSplitColl = collSlice.size();
      } else if constexpr (CentEstimator == CentralityEstimator::None) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, collSlice, centrality);
        nSplitColl = collSlice.size();
      }
      hfEvSelMc.fillHistograms<CentEstimator>(mcCollision, rejectionMask, nSplitColl);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject all particles from this collision
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
            rowMCMatchGenXicToXiPi(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
            rowMCMatchGenOmegacToXiPi(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
            rowMCMatchGenToOmegaPi(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
          } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
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
        float const kRapidityCutTight = 0.5;
        float const kRapidityCutLoose = 0.8;

        // Reject particles from background events
        if (particle.fromBackgroundEvent() && rejectBackground) {
          if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
            rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
            rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
            rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
            rowMCMatchGenToOmegaK(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }
          continue;
        }

        if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
          //  Xic → Xi pi
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, +kXiC0, std::array{+kXiMinus, +kPiPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != +kXiMinus) {
                continue;
              }
              // Xi -> Lambda pi
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, +kXiMinus, std::array{+kLambda0, +kPiMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != +kLambda0) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
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
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutTight) {
              hGenCharmBaryonPtRapidityTightXicToXiPi->SetBinContent(hGenCharmBaryonPtRapidityTightXicToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightXicToXiPi->GetBinContent(hGenCharmBaryonPtRapidityTightXicToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutLoose) {
              hGenCharmBaryonPtRapidityLooseXicToXiPi->SetBinContent(hGenCharmBaryonPtRapidityLooseXicToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseXicToXiPi->GetBinContent(hGenCharmBaryonPtRapidityLooseXicToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenXicToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }

        } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
          //  Omegac → Xi pi
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, +kOmegaC0, std::array{+kXiMinus, +kPiPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != +kXiMinus) {
                continue;
              }
              // Xi -> Lambda pi
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, +kXiMinus, std::array{+kLambda0, +kPiMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != +kLambda0) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
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
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutTight) {
              hGenCharmBaryonPtRapidityTightOmegacToXiPi->SetBinContent(hGenCharmBaryonPtRapidityTightOmegacToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightOmegacToXiPi->GetBinContent(hGenCharmBaryonPtRapidityTightOmegacToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutLoose) {
              hGenCharmBaryonPtRapidityLooseOmegacToXiPi->SetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToXiPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseOmegacToXiPi->GetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToXiPi->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenOmegacToXiPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }

        } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
          //  Omegac → Omega pi
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, +kOmegaC0, std::array{+kOmegaMinus, +kPiPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != +kOmegaMinus) {
                continue;
              }
              // Omega -> Lambda K
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, +kOmegaMinus, std::array{+kLambda0, +kKMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != +kLambda0) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
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
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutTight) {
              hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->SetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->GetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaPi->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutLoose) {
              hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->SetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->GetBinContent(hGenCharmBaryonPtRapidityLooseOmegacToOmegaPi->FindBin(ptCharmBaryonGen)) + 1);
            }
          }
          if (origin == RecoDecay::OriginType::NonPrompt) {
            rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, idxBhadMothers[0]);
          } else {
            rowMCMatchGenToOmegaPi(flag, debugGenCharmBar, debugGenCasc, debugGenLambda, ptCharmBaryonGen, rapidityCharmBaryonGen, origin, -1);
          }

        } else if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK) {
          //  Omegac → Omega K
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, +kOmegaC0, std::array{+kOmegaMinus, +kKPlus}, true, &sign)) {
            debugGenCharmBar = 1;
            ptCharmBaryonGen = particle.pt();
            rapidityCharmBaryonGen = particle.y();
            for (const auto& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != +kOmegaMinus) {
                continue;
              }
              // Omega -> Lambda K
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, +kOmegaMinus, std::array{+kLambda0, +kKMinus}, true)) {
                debugGenCasc = 1;
                for (const auto& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != +kLambda0) {
                    continue;
                  }
                  // Lambda -> p pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
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
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutTight) {
              hGenCharmBaryonPtRapidityTightOmegacToOmegaK->SetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaK->FindBin(ptCharmBaryonGen), hGenCharmBaryonPtRapidityTightOmegacToOmegaK->GetBinContent(hGenCharmBaryonPtRapidityTightOmegacToOmegaK->FindBin(ptCharmBaryonGen)) + 1);
            }
            if (std::abs(rapidityCharmBaryonGen) < kRapidityCutLoose) {
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

  template <o2::hf_centrality::CentralityEstimator CentEstimator, aod::hf_cand_xic0_omegac0::DecayType DecayChannel, typename TMyRecoCand>
  void runXic0Omegac0McQa(TMyRecoCand const& candidates,
                          MyTracksWMc const&,
                          aod::McParticles const& mcParticles,
                          BCsInfo const&)
  {
    int indexRec = -1;
    int8_t sign = -9;
    int8_t signCasc = -9;
    int8_t signV0 = -9;

    for (const auto& candidate : candidates) {

      auto arrayDaughters = std::array{candidate.template bachelorFromCharmBaryon_as<MyTracksWMc>(), // bachelor <- charm baryon
                                       candidate.template bachelor_as<MyTracksWMc>(),                // bachelor <- cascade
                                       candidate.template posTrack_as<MyTracksWMc>(),                // p <- lambda
                                       candidate.template negTrack_as<MyTracksWMc>()};               // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.template bachelor_as<MyTracksWMc>(),
                                           candidate.template posTrack_as<MyTracksWMc>(),
                                           candidate.template negTrack_as<MyTracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.template posTrack_as<MyTracksWMc>(),
                                         candidate.template negTrack_as<MyTracksWMc>()};

      auto mcV0DauPos = arrayDaughtersV0[0].mcParticle();
      auto mcV0DauNeg = arrayDaughtersV0[1].mcParticle();
      auto mcXiBachelor = arrayDaughtersCasc[0].mcParticle();
      auto mcXic0Bachelor = arrayDaughters[0].mcParticle();

      // Xic0 -> xi pi matching
      if constexpr (DecayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
        // Lambda → p pi
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);
        if (indexRec > -1 && signV0 == 1) {
          auto mcV0 = mcParticles.rawIteratorAt(indexRec - mcParticles.offset());

          float const v0MassPull = (candidate.invMassLambda() - MassLambda0) / candidate.invMassV0Err();
          registry.fill(HIST("hV0MassPullVsPt"), candidate.v0Pt(), v0MassPull);

          float const v0DauPosXDelta = candidate.v0DauPosX() - mcV0DauPos.vx();
          float const v0DauPosYDelta = candidate.v0DauPosY() - mcV0DauPos.vy();
          float const v0DauPosZDelta = candidate.v0DauPosZ() - mcV0DauPos.vz();
          float const v0DauPosPt = mcV0DauPos.pt();
          float const v0DauPosXPull = v0DauPosXDelta / candidate.v0DauPosXError();
          float const v0DauPosYPull = v0DauPosYDelta / candidate.v0DauPosYError();
          float const v0DauPosZPull = v0DauPosZDelta / candidate.v0DauPosZError();

          float const v0DauNegXDelta = candidate.v0DauNegX() - mcV0DauNeg.vx();
          float const v0DauNegYDelta = candidate.v0DauNegY() - mcV0DauNeg.vy();
          float const v0DauNegZDelta = candidate.v0DauNegZ() - mcV0DauNeg.vz();
          float const v0DauNegPt = mcV0DauNeg.pt();
          float const v0DauNegXPull = v0DauNegXDelta / candidate.v0DauNegXError();
          float const v0DauNegYPull = v0DauNegYDelta / candidate.v0DauNegYError();
          float const v0DauNegZPull = v0DauNegZDelta / candidate.v0DauNegZError();

          float const v0XDelta = candidate.v0VtxX() - mcV0DauNeg.vx();
          float const v0YDelta = candidate.v0VtxY() - mcV0DauNeg.vy();
          float const v0ZDelta = candidate.v0VtxZ() - mcV0DauNeg.vz();
          float const v0Pt = mcV0.pt();
          float const v0XPull = v0XDelta / candidate.v0XError();
          float const v0YPull = v0YDelta / candidate.v0YError();
          float const v0ZPull = v0ZDelta / candidate.v0ZError();

          float const lambdaXDelta = candidate.v0X() - mcV0DauNeg.vx();
          float const lambdaYDelta = candidate.v0Y() - mcV0DauNeg.vy();
          float const lambdaZDelta = candidate.v0Z() - mcV0DauNeg.vz();
          registry.fill(HIST("hV0DauPosXDelta"), v0DauPosXDelta);
          registry.fill(HIST("hV0DauPosYDelta"), v0DauPosYDelta);
          registry.fill(HIST("hV0DauPosZDelta"), v0DauPosZDelta);
          registry.fill(HIST("hV0DauPosXDeltaVsPt"), v0DauPosPt, v0DauPosXDelta);
          registry.fill(HIST("hV0DauPosYDeltaVsPt"), v0DauPosPt, v0DauPosYDelta);
          registry.fill(HIST("hV0DauPosZDeltaVsPt"), v0DauPosPt, v0DauPosZDelta);
          registry.fill(HIST("hV0DauPosXPull"), v0DauPosXPull);
          registry.fill(HIST("hV0DauPosYPull"), v0DauPosYPull);
          registry.fill(HIST("hV0DauPosZPull"), v0DauPosZPull);
          registry.fill(HIST("hV0DauPosXPullVsPt"), v0DauPosPt, v0DauPosXPull);
          registry.fill(HIST("hV0DauPosYPullVsPt"), v0DauPosPt, v0DauPosYPull);
          registry.fill(HIST("hV0DauPosZPullVsPt"), v0DauPosPt, v0DauPosZPull);

          registry.fill(HIST("hV0DauNegXDelta"), v0DauNegXDelta);
          registry.fill(HIST("hV0DauNegYDelta"), v0DauNegYDelta);
          registry.fill(HIST("hV0DauNegZDelta"), v0DauNegZDelta);
          registry.fill(HIST("hV0DauNegXDeltaVsPt"), v0DauNegPt, v0DauNegXDelta);
          registry.fill(HIST("hV0DauNegYDeltaVsPt"), v0DauNegPt, v0DauNegYDelta);
          registry.fill(HIST("hV0DauNegZDeltaVsPt"), v0DauNegPt, v0DauNegZDelta);
          registry.fill(HIST("hV0DauNegXPull"), v0DauNegXPull);
          registry.fill(HIST("hV0DauNegYPull"), v0DauNegYPull);
          registry.fill(HIST("hV0DauNegZPull"), v0DauNegZPull);
          registry.fill(HIST("hV0DauNegXPullVsPt"), v0DauNegPt, v0DauNegXPull);
          registry.fill(HIST("hV0DauNegYPullVsPt"), v0DauNegPt, v0DauNegYPull);
          registry.fill(HIST("hV0DauNegZPullVsPt"), v0DauNegPt, v0DauNegZPull);

          registry.fill(HIST("hV0XDelta"), v0XDelta);
          registry.fill(HIST("hV0YDelta"), v0YDelta);
          registry.fill(HIST("hV0ZDelta"), v0ZDelta);
          registry.fill(HIST("hV0XDeltaVsPt"), v0Pt, v0XDelta);
          registry.fill(HIST("hV0YDeltaVsPt"), v0Pt, v0YDelta);
          registry.fill(HIST("hV0ZDeltaVsPt"), v0Pt, v0ZDelta);
          registry.fill(HIST("hV0XPull"), v0XPull);
          registry.fill(HIST("hV0YPull"), v0YPull);
          registry.fill(HIST("hV0ZPull"), v0ZPull);
          registry.fill(HIST("hV0XPullVsPt"), v0Pt, v0XPull);
          registry.fill(HIST("hV0YPullVsPt"), v0Pt, v0YPull);
          registry.fill(HIST("hV0ZPullVsPt"), v0Pt, v0ZPull);

          registry.fill(HIST("hLambdaXDelta"), lambdaXDelta);
          registry.fill(HIST("hLambdaYDelta"), lambdaYDelta);
          registry.fill(HIST("hLambdaZDelta"), lambdaZDelta);

          registry.fill(HIST("hV0DauPosPtRes"), (candidate.v0DauPosPt() - mcV0DauPos.pt()) / candidate.v0DauPosPt());
          registry.fill(HIST("hV0DauNegPtRes"), (candidate.v0DauNegPt() - mcV0DauNeg.pt()) / candidate.v0DauNegPt());
          registry.fill(HIST("hV0PtRes"), (candidate.v0Pt() - mcV0.pt()) / candidate.v0Pt());
          // Xi- → pi pi p
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, +kXiMinus, std::array{+kPiMinus, +kProton, +kPiMinus}, true, &signCasc, 2);
          if (indexRec > -1 && signCasc == 1) {
            // QA
            float const xiMassPull = (candidate.invMassCascade() - MassXiMinus) / candidate.invMassXiErr();
            registry.fill(HIST("hXiMassPullVsPt"), candidate.xiPt(), xiMassPull);

            float const xiBachelorXDelta = candidate.xiBachelorX() - mcXiBachelor.vx();
            float const xiBachelorYDelta = candidate.xiBachelorY() - mcXiBachelor.vy();
            float const xiBachelorZDelta = candidate.xiBachelorZ() - mcXiBachelor.vz();
            float const xiBachelorPt = mcXiBachelor.pt();
            float const xiBachelorXPull = xiBachelorXDelta / candidate.xiBachelorXError();
            float const xiBachelorYPull = xiBachelorYDelta / candidate.xiBachelorYError();
            float const xiBachelorZPull = xiBachelorZDelta / candidate.xiBachelorZError();

            auto mcXi = mcParticles.rawIteratorAt(indexRec - mcParticles.offset());

            float const xiXDelta = candidate.xiX() - mcXiBachelor.vx();
            float const xiYDelta = candidate.xiY() - mcXiBachelor.vy();
            float const xiZDelta = candidate.xiZ() - mcXiBachelor.vz();
            float const xiPt = mcXi.pt();
            float const xiXPull = xiXDelta / candidate.xiXError();
            float const xiYPull = xiYDelta / candidate.xiYError();
            float const xiZPull = xiZDelta / candidate.xiZError();

            float const cascXDelta = candidate.xDecayVtxCascade() - mcXiBachelor.vx();
            float const cascYDelta = candidate.yDecayVtxCascade() - mcXiBachelor.vy();
            float const cascZDelta = candidate.zDecayVtxCascade() - mcXiBachelor.vz();

            registry.fill(HIST("hXiBachelorXDelta"), xiBachelorXDelta);
            registry.fill(HIST("hXiBachelorYDelta"), xiBachelorYDelta);
            registry.fill(HIST("hXiBachelorZDelta"), xiBachelorZDelta);
            registry.fill(HIST("hXiBachelorXDeltaVsPt"), xiBachelorPt, xiBachelorXDelta);
            registry.fill(HIST("hXiBachelorYDeltaVsPt"), xiBachelorPt, xiBachelorYDelta);
            registry.fill(HIST("hXiBachelorZDeltaVsPt"), xiBachelorPt, xiBachelorZDelta);
            registry.fill(HIST("hXiBachelorXPull"), xiBachelorXPull);
            registry.fill(HIST("hXiBachelorYPull"), xiBachelorYPull);
            registry.fill(HIST("hXiBachelorZPull"), xiBachelorZPull);
            registry.fill(HIST("hXiBachelorXPullVsPt"), xiBachelorPt, xiBachelorXPull);
            registry.fill(HIST("hXiBachelorYPullVsPt"), xiBachelorPt, xiBachelorYPull);
            registry.fill(HIST("hXiBachelorZPullVsPt"), xiBachelorPt, xiBachelorZPull);

            registry.fill(HIST("hXiXDelta"), xiXDelta);
            registry.fill(HIST("hXiYDelta"), xiYDelta);
            registry.fill(HIST("hXiZDelta"), xiZDelta);
            registry.fill(HIST("hXiXDeltaVsPt"), xiPt, xiXDelta);
            registry.fill(HIST("hXiYDeltaVsPt"), xiPt, xiYDelta);
            registry.fill(HIST("hXiZDeltaVsPt"), xiPt, xiZDelta);
            registry.fill(HIST("hXiXPull"), xiXPull);
            registry.fill(HIST("hXiYPull"), xiYPull);
            registry.fill(HIST("hXiZPull"), xiZPull);
            registry.fill(HIST("hXiXPullVsPt"), xiPt, xiXPull);
            registry.fill(HIST("hXiYPullVsPt"), xiPt, xiYPull);
            registry.fill(HIST("hXiZPullVsPt"), xiPt, xiZPull);

            registry.fill(HIST("hCascXDelta"), cascXDelta);
            registry.fill(HIST("hCascYDelta"), cascYDelta);
            registry.fill(HIST("hCascZDelta"), cascZDelta);

            registry.fill(HIST("hXiBachelorPtRes"), (candidate.xiBachelorPt() - mcXiBachelor.pt()) / candidate.xiBachelorPt());
            registry.fill(HIST("hXiPtRes"), (candidate.xiPt() - mcXi.pt()) / candidate.xiPt());

            // Xic → pi pi pi p
            indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, +kXiC0, std::array{+kPiPlus, +kPiMinus, +kProton, +kPiMinus}, true, &sign, 3);
            if (indexRec > -1 && sign == 1) {
              auto mcXic0 = mcParticles.rawIteratorAt(indexRec - mcParticles.offset());
              float const xic0MassPull = (candidate.invMassCharmBaryon() - MassXiC0) / candidate.invMassXic0Err();
              registry.fill(HIST("hXic0MassPullVsPt"), candidate.xic0Pt(), xic0MassPull);

              float const xic0BachelorXDelta = candidate.xic0BachelorX() - mcXic0Bachelor.vx();
              float const xic0BachelorYDelta = candidate.xic0BachelorY() - mcXic0Bachelor.vy();
              float const xic0BachelorZDelta = candidate.xic0BachelorZ() - mcXic0Bachelor.vz();
              float const xic0BachelorPt = mcXic0Bachelor.pt();
              float const xic0BachelorXPull = xic0BachelorXDelta / candidate.xic0BachelorXError();
              float const xic0BachelorYPull = xic0BachelorYDelta / candidate.xic0BachelorYError();
              float const xic0BachelorZPull = xic0BachelorZDelta / candidate.xic0BachelorZError();

              float const xic0XDelta = candidate.xDecayVtxCharmBaryon() - mcXic0Bachelor.vx();
              float const xic0YDelta = candidate.yDecayVtxCharmBaryon() - mcXic0Bachelor.vy();
              float const xic0ZDelta = candidate.zDecayVtxCharmBaryon() - mcXic0Bachelor.vz();
              float const xic0Pt = mcXic0.pt();
              float const xic0XPull = xic0XDelta / candidate.xic0XError();
              float const xic0YPull = xic0YDelta / candidate.xic0YError();
              float const xic0ZPull = xic0ZDelta / candidate.xic0ZError();
              registry.fill(HIST("hXic0BachelorXDelta"), xic0BachelorXDelta);
              registry.fill(HIST("hXic0BachelorYDelta"), xic0BachelorYDelta);
              registry.fill(HIST("hXic0BachelorZDelta"), xic0BachelorZDelta);
              registry.fill(HIST("hXic0BachelorXDeltaVsPt"), xic0BachelorPt, xic0BachelorXDelta);
              registry.fill(HIST("hXic0BachelorYDeltaVsPt"), xic0BachelorPt, xic0BachelorYDelta);
              registry.fill(HIST("hXic0BachelorZDeltaVsPt"), xic0BachelorPt, xic0BachelorZDelta);
              registry.fill(HIST("hXic0BachelorXPull"), xic0BachelorXPull);
              registry.fill(HIST("hXic0BachelorYPull"), xic0BachelorYPull);
              registry.fill(HIST("hXic0BachelorZPull"), xic0BachelorZPull);
              registry.fill(HIST("hXic0BachelorXPullVsPt"), xic0BachelorPt, xic0BachelorXPull);
              registry.fill(HIST("hXic0BachelorYPullVsPt"), xic0BachelorPt, xic0BachelorYPull);
              registry.fill(HIST("hXic0BachelorZPullVsPt"), xic0BachelorPt, xic0BachelorZPull);

              registry.fill(HIST("hXic0XDelta"), xic0XDelta);
              registry.fill(HIST("hXic0YDelta"), xic0YDelta);
              registry.fill(HIST("hXic0ZDelta"), xic0ZDelta);
              registry.fill(HIST("hXic0XDeltaVsPt"), xic0Pt, xic0XDelta);
              registry.fill(HIST("hXic0YDeltaVsPt"), xic0Pt, xic0YDelta);
              registry.fill(HIST("hXic0ZDeltaVsPt"), xic0Pt, xic0ZDelta);
              registry.fill(HIST("hXic0XPull"), xic0XPull);
              registry.fill(HIST("hXic0YPull"), xic0YPull);
              registry.fill(HIST("hXic0ZPull"), xic0ZPull);
              registry.fill(HIST("hXic0XPullVsPt"), xic0Pt, xic0XPull);
              registry.fill(HIST("hXic0YPullVsPt"), xic0Pt, xic0YPull);
              registry.fill(HIST("hXic0ZPullVsPt"), xic0Pt, xic0ZPull);

              registry.fill(HIST("hXic0BachelorPtRes"), (candidate.xic0BachelorPt() - mcXic0Bachelor.pt()) / candidate.xic0BachelorPt());
              registry.fill(HIST("hXic0PtRes"), (candidate.xic0Pt() - mcXic0.pt()) / candidate.xic0Pt());
            }
          }
        }
      }
    }
  }

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

  void processMcXicToXiPiKfQa(aod::HfCandToXiPiKfQa const& candidates,
                              MyTracksWMc const& tracks,
                              aod::McParticles const& mcParticles,
                              BCsInfo const& bcs)
  {
    runXic0Omegac0McQa<CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Mc, processMcXicToXiPiKfQa, "Run Xic0 to xi pi MC QA process function - no centrality", false);

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
