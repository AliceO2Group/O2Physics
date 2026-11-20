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

/// \file candidateCreatorXic0Omegac0Qa.cxx
/// \brief Reconstruction of Xic0 and Xicp candiates with hadronic decay chain
///
/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h" // -> Added to test removal of strangeness builder workflow

#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::constants::physics;
using namespace o2::hf_evsel;

constexpr int nMcDecayChannels = hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::NChannelsMain;

struct HfCandidateCreatorXic0Omegac0Qa {
  Produces<HfCandCascBachBaseQa> rowCandCascBatch;
  // Additional information for candidates created with KFParticle
  Produces<HfCandCascBachKfQa> rowCandCascBatchKf;

  struct : ConfigurableGroup {
    std::string prefix = "lfBuilding";
    // Options for internal cascade building - DCAFitter settings
    // ...Initial vaules taken from PWGLF/Utils/strangenessBuilderHelper
    Configurable<float> maxR{"maxR", 200., "Reject PCA's above this radius"};
    Configurable<float> maxDZIni{"maxDZIni", 1e9, "Reject (if>0) PCA candidate if tracks DZ exceeds this threshold"};
    Configurable<float> maxDXYIni{"maxDXYIni", 4.0f, "Reject (if>0) PCA candidate if tracks DXY exceeds this threshold"};
    Configurable<float> minParamChange{"minParamChange", 1.e-3, "Stop iteration if largest change of any X is smaller than this"};
    Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "Stop iteration if Chi2/Chi2old > this"};
    Configurable<float> maxChi2{"maxChi2", 1e9, "Discard vertices with Chi2/Nprongs > this(or sum {DCAi^2}/Nprongs for abs. distance minimization)"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariance, effective only if useAbsDCA is true"};
    Configurable<bool> refitWithMaterialCorrection{"refitWithMaterialCorrection", false, "Do refit after material correction applied"};
    // Options for internal V0 building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    Configurable<int> minCrossedRows{"minCrossedRows", 50, "minimun TPC crossed rows for daughter tracks. Used for internal V0 Building"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "DCV Neg to PV"};
    Configurable<float> dcapostopv{"dcapostopv", 0.1, "DCV Pos To PV"};
    Configurable<float> v0cospa{"v0cospa", 0.95, "V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> v0radius{"v0radius", 0.9, "v0radius"};
    Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximun daughter eta (in abs value)"};
    // Options for internal cascade building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    Configurable<float> dcabachtopv{"dcabachtopv", 0.05, "DCV Bach to PV"};
    Configurable<float> cascradius{"cascradius", 0.9, "DCV Bach to PV"};
    Configurable<float> casccospa{"casccospa", 0.95, "Cascade CosPA"};
    Configurable<float> dcacascdau{"dcacascdau", 1.0, "DCA cascade daughters"};
    Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.01, "Distance from Lambda mass"};
    // Options for LF and charm-baryon building - KF Building specifics
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "0: simplifyed fast mathematics ignoring constraint on fixed mass, 2 : daughter particle masses stay fixed in construction process"};
    Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", true, "KF : Use Lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", true, "KF : Use Cascade mass constraint - WARNING : Not adequate for inv mass analysis of Xi and Omega"};
    Configurable<bool> kfDoDCAFitterPreMinimV0{"kfDoDCAFitterPreMinimV0", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for V0"};
    Configurable<bool> kfDoDCAFitterPreMinimCasc{"kfDoDCAFitterPreMinimCasc", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for Xi and Omega"};
  } lfConfigs;

  struct : ConfigurableGroup {
    // Switch for filling histograms
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill validation plots"};
    // Magnetic field setting from CCDB
    Configurable<bool> isRun2{"isRun2", false, "enable Run2 or Run3 GRP objects for magnetic field"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "https://alice-ccdb.cern.ch", "url of the ccdb object"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parameterization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "path of the group file (Run2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run3)"};
    // Flags to enable creation for different particle species separately
    Configurable<bool> createToXiPi{"createToXiPi", false, "enable XiC0 (OmegaC0) -> Xi pi candidate creation"};
    Configurable<bool> createOmegaCToOmegaPi{"createOmegaCToOmegaPi", false, "enable OmegaC0 -> Omega Pi candidate creation"};
    Configurable<bool> createOmegaCToOmegaK{"createOmegaCToOmegaK", false, "enable OmegaC0 -> Omega K candidate creation"};
    // Cascade pre selection
    Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01, "Invariant mass tolerance for cascades"};
    // Options for DCAFitter
    Configurable<float> maxR{"maxR", 200., "Reject PCA's above this radius"};
    Configurable<float> maxDZIni{"maxDZIni", 4., "Reject (if>0) PCA candidate if tracks DZ exceeds this threshold"};
    Configurable<float> maxDXYIni{"maxDXYIni", 4., "Reject (if>0) PCA candidate if tracks DXY exceeds this threshold"};
    Configurable<float> minParamChange{"minParamChange", 1.e-3, "Stop iteration if largest change of any X is smaller than this"};
    Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "Stop iteration if Chi2/Chi2old > this"};
    Configurable<float> maxChi2{"maxChi2", 100, "Discard vertices with Chi2/Nprongs > this(or sum {DCAi^2}/Nprongs for abs. distance minimization)"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariance, effective only if useAbsDCA is true"};
  } configs;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  Service<o2::framework::O2DatabasePDG> pdgdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // For cascade building using LF strangeness builder
  o2::pwglf::strangenessBuilderHelper straHelper;
  // DCAFitter
  o2::vertexing::DCAFitterN<2> df;

  int runNumber{0};
  double magneticField{0.};

  enum CharmBaryonCandCounter {
    All = 0,
    HfFlagPass,
    CascReconstructed,
    VertexFit
  };

  // Table aliases
  using FilteredHfCascLf2Prongs = soa::Filtered<aod::HfCascLf2Prongs>;
  using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
  using TracksWCovIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
  using TracksWCovDcaExtraPidPrPiKa = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPr, aod::TracksPidPi, aod::TracksPidKa>;
  using TracksWCovExtraPidIU = soa::Join<aod::TracksIU, TracksCovIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  Filter filterSelCascBach = (configs.createToXiPi && (hf_track_index::hfflag & static_cast<uint8_t>(BIT(hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi))) != static_cast<uint8_t>(0)) ||
                             (configs.createOmegaCToOmegaPi && (hf_track_index::hfflag & static_cast<uint8_t>(BIT(hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi))) != static_cast<uint8_t>(0)) ||
                             (configs.createOmegaCToOmegaK && (hf_track_index::hfflag & static_cast<uint8_t>(BIT(hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK))) != static_cast<uint8_t>(0));

  HistogramRegistry registry{"hists"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  HfEventSelection hfEvSel;

  // Pointer of histograms for QA
  std::shared_ptr<TH1> hInvMassCharmBaryonToXiPi, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaKa;
  std::shared_ptr<TH1> hCandidateCounterToXiPi, hCandidateCounterToOmegaPi, hCandidateCounterToOmegaKa;

  void init(InitContext const&)
  {
    // Exit if workflow is not configured correctly
    std::vector<bool> doprocessDCA{doprocessWithDCAFitterNoCent, /*doprocessWithDCAFitterNoCentWithTrackedCasc,*/ doprocessWithDCAFitterCentFT0C, doprocessWithDCAFitterCentFT0M};
    std::vector<bool> doprocessKF{doprocessWithKFParticleNoCent, doprocessWithKFParticleCentFT0C, doprocessWithKFParticleCentFT0M};
    std::vector<bool> processesCollMonitoring{doprocessCollisionsNoCent, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};
    std::vector<bool> creationFlags = {configs.createToXiPi, configs.createOmegaCToOmegaPi, configs.createOmegaCToOmegaK};

    int nProcessesDCA = std::accumulate(doprocessDCA.begin(), doprocessDCA.end(), 0);
    int nProcessesKF = std::accumulate(doprocessKF.begin(), doprocessKF.end(), 0);
    if (nProcessesDCA + nProcessesKF != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

    int nProcessesCollisions = std::accumulate(processesCollMonitoring.begin(), processesCollMonitoring.end(), 0);
    if (nProcessesCollisions > 1) {
      LOGP(fatal, "At most one process function for collision monitoring can be enabled at a time.");
    }
    if (nProcessesCollisions == 1) {
      if ((doprocessWithDCAFitterNoCent || doprocessWithKFParticleNoCent) && !doprocessCollisionsNoCent) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsNoCent\"?");
      }
      if ((doprocessWithDCAFitterCentFT0C || doprocessWithKFParticleCentFT0C) && !doprocessCollisionsCentFT0C) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0C\"?");
      }
      if ((doprocessWithDCAFitterCentFT0M || doprocessWithKFParticleCentFT0M) && !doprocessCollisionsCentFT0M) {
        LOGP(fatal, "Process function for collision monitoring not correctly enabled. Did you enable \"processCollisionsCentFT0M\"?");
      }
    }

    if (std::accumulate(creationFlags.begin(), creationFlags.end(), 0) == 0) {
      LOGP(fatal, "At least one particle specie should be enabled for the creation.");
    }

    LOGF(info, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    LOGF(info, "  Chosen reconstruction info: ToXiPi(%d)/ ToOmegaPi(%d)/ ToOmegaKa(%d)", configs.createToXiPi.value, configs.createOmegaCToOmegaPi.value, configs.createOmegaCToOmegaK.value);
    LOGF(info, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

    // Add histogram to indicate which sv method was used
    registry.add("hVertexerType", "Use KF or DCAFitterN;Vertexer type;entries", {kTH1F, {{2, 0.0, 2.0}}});
    registry.get<TH1>(HIST("hVertexerType"))->GetXaxis()->SetBinLabel(1 + aod::hf_cand::VertexerType::DCAFitter, "DCAFitter");
    registry.get<TH1>(HIST("hVertexerType"))->GetXaxis()->SetBinLabel(1 + aod::hf_cand::VertexerType::KfParticle, "KFParticle");

    // initialize ccdb
    // ---------------
    ccdb->setURL(configs.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(configs.ccdbPathLut));
    straHelper.lut = lut;
    runNumber = 0;

    // Initilization for strangeness builder helper
    // Settings for internal V0 building
    straHelper.v0selections.minCrossedRows = lfConfigs.minCrossedRows;
    straHelper.v0selections.dcanegtopv = lfConfigs.dcanegtopv;
    straHelper.v0selections.dcapostopv = lfConfigs.dcapostopv;
    straHelper.v0selections.v0cospa = lfConfigs.v0cospa;
    straHelper.v0selections.dcav0dau = lfConfigs.dcav0dau;
    straHelper.v0selections.v0radius = lfConfigs.v0radius;
    straHelper.v0selections.maxDaughterEta = lfConfigs.maxDaughterEta;

    // Settings for internal Cascade building
    straHelper.cascadeselections.minCrossedRows = lfConfigs.minCrossedRows;
    straHelper.cascadeselections.dcabachtopv = lfConfigs.dcabachtopv;
    straHelper.cascadeselections.cascradius = lfConfigs.cascradius;
    straHelper.cascadeselections.casccospa = lfConfigs.casccospa;
    straHelper.cascadeselections.dcacascdau = lfConfigs.dcacascdau;
    straHelper.cascadeselections.lambdaMassWindow = lfConfigs.lambdaMassWindow;
    straHelper.cascadeselections.maxDaughterEta = lfConfigs.maxDaughterEta;

    // Fitter setting
    straHelper.fitter.setPropagateToPCA(true);
    straHelper.fitter.setMaxR(lfConfigs.maxR);
    straHelper.fitter.setMaxDZIni(lfConfigs.maxDZIni);
    straHelper.fitter.setMaxDXYIni(lfConfigs.maxDXYIni);
    straHelper.fitter.setMinParamChange(lfConfigs.minParamChange);
    straHelper.fitter.setMinRelChi2Change(lfConfigs.minRelChi2Change);
    straHelper.fitter.setMaxChi2(lfConfigs.maxChi2);
    straHelper.fitter.setUseAbsDCA(lfConfigs.useAbsDCA);
    straHelper.fitter.setWeightedFinalPCA(lfConfigs.useWeightedFinalPCA);
    straHelper.fitter.setRefitWithMatCorr(lfConfigs.refitWithMaterialCorrection);

    // Extra initialization for DCAFitter
    if (nProcessesDCA == 1) {
      registry.get<TH1>(HIST("hVertexerType"))->Fill(aod::hf_cand::VertexerType::DCAFitter);
      df.setPropagateToPCA(true);
      df.setMaxR(configs.maxR);
      df.setMaxDZIni(configs.maxDZIni);
      df.setMaxDXYIni(configs.maxDXYIni);
      df.setMinParamChange(configs.minParamChange);
      df.setMinRelChi2Change(configs.minRelChi2Change);
      df.setMaxChi2(configs.maxChi2);
      df.setUseAbsDCA(configs.useAbsDCA);
      df.setWeightedFinalPCA(configs.useWeightedFinalPCA);
    }

    // Extra initialization for KFParticle
    if (nProcessesKF == 1) {
      registry.get<TH1>(HIST("hVertexerType"))->Fill(aod::hf_cand::VertexerType::KfParticle);
    }

    // initailize HF event selection helper
    hfEvSel.init(registry, &zorroSummary);

    // Histograms for QA
    registry.add("ReconstructedDecayChannel", "DecyayChannel", {kTH1F, {{3, 0.0, 3.0}}});
    registry.get<TH1>(HIST("ReconstructedDecayChannel"))->GetXaxis()->SetBinLabel(1 + hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi, "To #Xi #pi");
    registry.get<TH1>(HIST("ReconstructedDecayChannel"))->GetXaxis()->SetBinLabel(1 + hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi, "To #Omega #pi");
    registry.get<TH1>(HIST("ReconstructedDecayChannel"))->GetXaxis()->SetBinLabel(1 + hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK, "To #Omega K");

    hInvMassCharmBaryonToXiPi = registry.add<TH1>("hInvMassCharmBaryonToXiPi", "Charm baryon invariant mass - #Xi #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1800, 1.8, 3.6}}});
    hCandidateCounterToXiPi = registry.add<TH1>("hCandidateCounterToXiPi", "Candidate counter wrt derived data - #Xi #pi decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});
    registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + All, "Total");
    registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + HfFlagPass, "HfFlagPass");
    registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + CascReconstructed, "CascReconstructed");
    registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + VertexFit, "VertexFit");
    if (configs.createToXiPi) {
      registry.get<TH1>(HIST("ReconstructedDecayChannel"))->Fill(hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi);
    }

    hInvMassCharmBaryonToOmegaPi = registry.add<TH1>("hInvMassCharmBaryonToOmegaPi", "Charm baryon invariant mass - #Omega #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1800, 1.8, 3.6}}});
    hCandidateCounterToOmegaPi = registry.add<TH1>("hCandidateCounterToOmegaPi", "Candidate counter wrt derived data - #Omega #pi decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});
    registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + All, "Total");
    registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + HfFlagPass, "HfFlagPass");
    registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + CascReconstructed, "CascReconstructed");
    registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + VertexFit, "VertexFit");
    if (configs.createOmegaCToOmegaPi) {
      registry.get<TH1>(HIST("ReconstructedDecayChannel"))->Fill(hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi);
    }

    hInvMassCharmBaryonToOmegaKa = registry.add<TH1>("hInvMassCharmBaryonToOmegaKa", "Charm baryon invariant mass - #Omega K decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{1800, 1.8, 3.6}}});
    hCandidateCounterToOmegaKa = registry.add<TH1>("hCandidateCounterToOmegaKa", "Candidate counter wrt derived data - #Omega K decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});
    registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + All, "Total");
    registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + HfFlagPass, "HfFlagPass");
    registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + CascReconstructed, "CascReconstructed");
    registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + VertexFit, "VertexFit");
    if (configs.createOmegaCToOmegaK) {
      registry.get<TH1>(HIST("ReconstructedDecayChannel"))->Fill(hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK);
    }

    registry.add("hCascMass", "Inv mass of reconstructed cascade;Inv mass;Entries", {HistType::kTH1F, {{1000, 1.0, 2.0}}});
    registry.add("hCascPt", "Pt of reconstructed cascade;pT;Entries", {HistType::kTH1F, {{100, 0., 20.0}}});

  } // end of initialization

  // helper function to convert PDG of KFParticle object into appropriate track PID.
  // exclusively used to make input of getTrackParCovFromKFP()
  // \brief kfPdg is Pdg value of KFParticle object
  o2::track::PID::ID convertPDGIntoPID(const int kfPdg)
  {
    switch (std::abs(kfPdg)) {
      case kPiPlus:
        return o2::track::PID::Pion;
      case kKPlus:
        return o2::track::PID::Kaon;
      case kProton:
        return o2::track::PID::Proton;
      case kLambda0:
        return o2::track::PID::Lambda;
      case kXiMinus:
        return o2::track::PID::XiMinus;
      case kOmegaMinus:
        return o2::track::PID::OmegaMinus;
      default:
        LOGF(fatal, "Undefined PDG value from KFParticle given for conversion(%d). Replace it with PID of Pion", kfPdg);
        return o2::track::PID::Pion;
    }
  }

  // template function for running charm-baryon reconstruction with DCAFitter for single candidate and decay channel
  /// \brief decayChannel is for different decay channels
  /// \brief Colls is for collision tables joined with different centraltiy estimators
  /// \brief Hist is for QA histograms
  template <int decayChannel, typename TColl, typename TLfTrack, typename TTrack, typename Hist>
  void buildCandidateWithDCAFitter(TColl const& collision,
                                   TLfTrack const& posTrack,
                                   TLfTrack const& negTrack,
                                   TLfTrack const& bachTrack,
                                   TTrack const& trackCharmBachelor,
                                   Hist& hInvMassCharmBaryon,
                                   Hist& hCandCounter)
  {
    // define constants for decay channel (default to values for Xic0/Omegac0 -> Xi pi)
    int trackPidOfCascade = o2::track::PID::XiMinus;
    double massOfCascade = o2::constants::physics::MassXiMinus;
    double massOfCharmBach = o2::constants::physics::MassPiPlus;

    if (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
      trackPidOfCascade = o2::track::PID::OmegaMinus;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
      massOfCharmBach = o2::constants::physics::MassPiPlus;
    }

    if (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
      trackPidOfCascade = o2::track::PID::OmegaMinus;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
      massOfCharmBach = o2::constants::physics::MassKPlus;
    }

    // Make V0 and cascade starting from daugter tracks
    // If success, fill Cascade and V0 information for reconstruction
    if (!straHelper.buildCascadeCandidate(collision.globalIndex(),
                                          collision.posX(), collision.posY(), collision.posZ(),
                                          posTrack,
                                          negTrack,
                                          bachTrack,
                                          false, // calculateBachelorBaryonVariable
                                          false, // useCascadeMomentumAtPV
                                          true)) {
      return;
    }

    float mCasc = (decayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) ? straHelper.cascade.massOmega : straHelper.cascade.massXi;
    if (std::abs(mCasc - massOfCascade) > configs.cascadeMassWindow) {
      return;
    }

    registry.fill(HIST("hCascMass"), mCasc);
    registry.fill(HIST("hCascPt"), RecoDecay::pt(straHelper.cascade.cascadeMomentum));
    if (configs.fillHistograms) {
      hCandCounter->Fill(CascReconstructed);
    }

    //------------------------------Create cascade track------------------------------
    if (bachTrack.sign() == 0) {
      return;
    }
    const int bachTrackSign = bachTrack.sign() < 0 ? -1 : 1;

    std::array<float, 21> covCasc = {0.};
    std::copy(std::begin(straHelper.cascade.covariance), std::end(straHelper.cascade.covariance), covCasc.begin());
    o2::track::TrackParCov trackCasc = o2::track::TrackParCov(straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum, covCasc, bachTrackSign, true);
    trackCasc.setAbsCharge(1);
    trackCasc.setPID(trackPidOfCascade);

    //------------------------------Fit SV & Create Charm Baryon track------------------------------
    auto trackParCovCharmBachelor = getTrackParCov(trackCharmBachelor);
    if (df.process(trackCasc, trackParCovCharmBachelor) == 0) {
      return;
    }

    if (configs.fillHistograms) {
      hCandCounter->Fill(VertexFit);
    }

    //------------------------------Calculate physical properties-----------------------------
    // get track momenta
    std::array<float, 3> pVecCasc{}, pVecCharmBach{};
    df.getTrack(0).getPxPyPzGlo(pVecCasc);
    df.getTrack(1).getPxPyPzGlo(pVecCharmBach);

    // get PV Properties
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();
    std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

    // get SV Properties
    auto const& secondaryVertex = df.getPCACandidate();
    auto covMatrixSV = df.calcPCACovMatrixFlat();
    std::array<float, 3> svCoord = {static_cast<float>(secondaryVertex[0]),
                                    static_cast<float>(secondaryVertex[1]),
                                    static_cast<float>(secondaryVertex[2])};

    float phi{}, theta{};
    getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
    float errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
    float errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

    // Impact parameters computed with propagatToDCABxByBz method
    auto trackParCovV0DauPos = getTrackParCov(posTrack);
    auto trackParCovV0DauNeg = getTrackParCov(negTrack);
    auto trackParCovBach = getTrackParCov(bachTrack);

    o2::dataformats::DCA impactParameterV0DauPos, impactParameterV0DauNeg, impactParameterCascBach;
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0DauPos, 2.f, matCorr, &impactParameterV0DauPos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0DauNeg, 2.f, matCorr, &impactParameterV0DauNeg);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovBach, 2.f, matCorr, &impactParameterCascBach);

    o2::dataformats::DCA impactParameterCasc, impactParameterCharmBach;
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCharmBachelor, 2.f, matCorr, &impactParameterCharmBach);

    // get Lambda invariant mass - from LF Table
    float mLambda = (straHelper.cascade.charge < 0) ? straHelper.v0.massLambda : straHelper.v0.massAntiLambda;

    // calculate quantities related to cascade and v0
    float cpaCasc = RecoDecay::cpa(pvCoord, straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum);
    float cpaxyCasc = RecoDecay::cpaXY(pvCoord, straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum);
    float cpaV0 = RecoDecay::cpa(pvCoord, straHelper.cascade.v0Position, straHelper.cascade.v0Momentum);
    float cpaxyV0 = RecoDecay::cpaXY(pvCoord, straHelper.cascade.v0Position, straHelper.cascade.v0Momentum);
    float cpaCascToCharmHadron = RecoDecay::cpa(svCoord, straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum);
    float cpaxyCascToCharmHadron = RecoDecay::cpaXY(svCoord, straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum);
    float cpaV0ToCascade = RecoDecay::cpa(straHelper.cascade.cascadePosition, straHelper.cascade.v0Position, straHelper.cascade.v0Momentum);
    float cpaxyV0ToCascade = RecoDecay::cpaXY(straHelper.cascade.cascadePosition, straHelper.cascade.v0Position, straHelper.cascade.v0Momentum);

    float etaCasc = RecoDecay::eta(straHelper.cascade.cascadeMomentum);
    float etaV0 = RecoDecay::eta(straHelper.cascade.v0Momentum);

    // fill control histogram
    if (configs.fillHistograms) {
      float massCharmBaryonCand = RecoDecay::m(std::array{pVecCasc, pVecCharmBach}, std::array{massOfCascade, massOfCharmBach});
      hInvMassCharmBaryon->Fill(massCharmBaryonCand);
    }

    //------------------------------Fill the table-----------------------------
    rowCandCascBatch(collision.globalIndex(),
                     primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                     secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                     straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2],
                     straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                     errorDecayLength, errorDecayLengthXY,
                     pVecCasc[0], pVecCasc[1], pVecCasc[2],
                     pVecCharmBach[0], pVecCharmBach[1], pVecCharmBach[2],
                     impactParameterCasc.getY(), impactParameterCharmBach.getY(),
                     std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBach.getSigmaY2()),
                     impactParameterCasc.getZ(), impactParameterCharmBach.getZ(),
                     std::sqrt(impactParameterCasc.getSigmaZ2()), std::sqrt(impactParameterCharmBach.getSigmaZ2()),
                     impactParameterV0DauPos.getY(), impactParameterV0DauNeg.getY(), impactParameterCascBach.getY(),
                     impactParameterV0DauPos.getZ(), impactParameterV0DauNeg.getZ(), impactParameterCascBach.getZ(),
                     std::sqrt(df.getChi2AtPCACandidate()), straHelper.cascade.cascadeDaughterDCA, straHelper.cascade.v0DaughterDCA,
                     cpaCasc, cpaxyCasc, cpaV0, cpaxyV0,
                     cpaCascToCharmHadron, cpaxyCascToCharmHadron, cpaV0ToCascade, cpaxyV0ToCascade,
                     mCasc, mLambda,
                     etaCasc, etaV0,
                     trackCharmBachelor.globalIndex(), bachTrack.globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                     BIT(decayChannel));
  }

  // template function for running charm-baryon reconstruction with KFParticle for single candidate and decay channel
  /// \brief decayChannel is for different decay channels
  /// \brief Colls is for collision tables joined with different centraltiy estimators
  /// \brief Hist is for QA histograms
  template <int decayChannel, typename TColl, typename TLfTrack, typename TTrack, typename Hist>
  void buildCandidateWithKFParticle(TColl const& collision,
                                    TLfTrack const& posTrack,
                                    TLfTrack const& negTrack,
                                    TLfTrack const& bachTrack,
                                    TTrack const& trackCharmBachelor,
                                    Hist& hInvMassCharmBaryon,
                                    Hist& hCandCounter)
  {
    // define constants for decay channel (default to values for Xic0/Omegac0 -> Xi pi)
    int pdgIdOfV0DauPos = kProton;
    int pdgIdOfV0DauNeg = kPiMinus;
    int pdgIdOfV0 = kLambda0;
    double massOfV0 = o2::constants::physics::MassLambda;

    int pdgIdOfBach = kPiMinus;
    int pdgIdOfCharmBach = kPiPlus;
    int pdgIdOfCascade = kXiMinus;
    int pdgIdOfCharmBaryon = kXiC0;
    double massOfCascade = o2::constants::physics::MassXiMinus;
    bool kfTuneForOmega = false;

    if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
      pdgIdOfBach = kKMinus;
      pdgIdOfCharmBach = kPiPlus;
      pdgIdOfCascade = kOmegaMinus;
      pdgIdOfCharmBaryon = kOmegaC0;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
      kfTuneForOmega = true;
    }

    if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
      pdgIdOfBach = kKMinus;
      pdgIdOfCharmBach = kKPlus;
      pdgIdOfCascade = kOmegaMinus;
      pdgIdOfCharmBaryon = kOmegaC0;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
      kfTuneForOmega = true;
    }

    // Make cascade starting from V0
    // If success, fill Cascade and V0 information for reconstruction
    if (!straHelper.buildCascadeCandidateWithKF(collision.globalIndex(),
                                                collision.posX(), collision.posY(), collision.posZ(),
                                                posTrack,
                                                negTrack,
                                                bachTrack,
                                                false, // calculateBachelorBaryonVariables
                                                lfConfigs.kfConstructMethod,
                                                kfTuneForOmega,
                                                lfConfigs.kfUseV0MassConstraint,
                                                lfConfigs.kfUseCascadeMassConstraint,
                                                lfConfigs.kfDoDCAFitterPreMinimV0,
                                                lfConfigs.kfDoDCAFitterPreMinimCasc)) {
      return;
    }

    //----------Create charm baryon as KF Particle object starting from V0----------
    // TODO: investigate if we can use the KFParticle object of cascade from straHelper instead of creating a new one here
    // Create KFParticle object of V0 Daughter & Bachelor
    const KFPTrack kfTrack0 = createKFPTrackFromTrack(posTrack);
    const KFPTrack kfTrack1 = createKFPTrackFromTrack(negTrack);
    const KFPTrack kfTrackBach = createKFPTrackFromTrack(bachTrack);

    bool isAnti = bachTrack.signed1Pt() > 0;
    int bachCharge = (bachTrack.signed1Pt() > 0 ? +1 : -1);

    KFParticle kfPos(kfTrack0, (isAnti ? -pdgIdOfV0DauNeg : pdgIdOfV0DauPos));
    KFParticle kfNeg(kfTrack1, (isAnti ? -pdgIdOfV0DauPos : pdgIdOfV0DauNeg));
    KFParticle kfBach(kfTrackBach, (isAnti ? -pdgIdOfBach : pdgIdOfBach));

    // ~~~~~~~Construct V0 with KF~~~~~~~
    const KFParticle* v0Daughters[2] = {&kfPos, &kfNeg};
    KFParticle kfV0;
    kfV0.SetConstructMethod(lfConfigs.kfConstructMethod);
    kfV0.Construct(v0Daughters, 2);

    // Require lambda pre-selection before mass constraint
    float massLam{}, sigMassLam{};
    if (kfV0.GetMass(massLam, sigMassLam)) {
      return;
    }
    if (std::abs(massLam - massOfV0) > lfConfigs.lambdaMassWindow) {
      return;
    }
    if (sigMassLam <= 0) {
      return;
    }
    if ((kfV0.GetNDF() <= 0) || (kfV0.GetChi2() <= 0)) {
      return;
    }
    kfV0.SetPDG(pdgIdOfV0);

    // Set mass constraint to lambda
    if (lfConfigs.kfUseV0MassConstraint) {
      kfV0.SetNonlinearMassConstraint(massOfV0);
    }
    kfV0.TransportToDecayVertex();

    //~~~~~~~Construct cascade with KF~~~~~~~
    const KFParticle* cascDaughters[2] = {&kfBach, &kfV0};
    KFParticle kfCasc;
    kfCasc.SetConstructMethod(lfConfigs.kfConstructMethod);
    kfCasc.Construct(cascDaughters, 2);

    float massCasc{}, sigMassCasc{};
    if (kfCasc.GetMass(massCasc, sigMassCasc)) {
      return;
    }
    if (std::abs(massCasc - massOfCascade) > configs.cascadeMassWindow) {
      return;
    }
    if (sigMassCasc <= 0) {
      return;
    }
    if (kfCasc.GetNDF() <= 0 || kfCasc.GetChi2() <= 0) {
      return;
    }
    kfCasc.SetPDG(pdgIdOfCascade);

    // Set mass constraint to cascade
    if (lfConfigs.kfUseCascadeMassConstraint) {
      kfCasc.SetNonlinearMassConstraint(massOfCascade);
    }
    kfCasc.TransportToDecayVertex();

    registry.fill(HIST("hCascMass"), massCasc);
    registry.fill(HIST("hCascPt"), RecoDecay::pt(kfCasc.GetPx(), kfCasc.GetPy()));
    if (configs.fillHistograms) {
      hCandCounter->Fill(CascReconstructed);
    }

    //~~~~~~~Construct Charm Baryon with KF~~~~~~~
    const KFPTrack kfTrackCharmBach = createKFPTrackFromTrack(trackCharmBachelor);
    const KFParticle kfCharmBach(kfTrackCharmBach, (isAnti ? -pdgIdOfCharmBach : pdgIdOfCharmBach));
    const KFParticle* charmBaryonDaughters[2] = {&kfCharmBach, &kfCasc};

    KFParticle kfCharmBaryon;
    kfCharmBaryon.SetConstructMethod(lfConfigs.kfConstructMethod);
    kfCharmBaryon.Construct(charmBaryonDaughters, 2);

    float massCharmBaryon{}, sigMassCharmBaryon{};
    if (kfCharmBaryon.GetMass(massCharmBaryon, sigMassCharmBaryon)) {
      return;
    }
    if (sigMassCharmBaryon <= 0) {
      return;
    }
    if (kfCharmBaryon.GetNDF() <= 0 || kfCharmBaryon.GetChi2() <= 0) {
      return;
    }
    kfCharmBaryon.SetPDG(pdgIdOfCharmBaryon);
    kfCharmBaryon.TransportToDecayVertex();

    if (configs.fillHistograms) {
      hCandCounter->Fill(VertexFit);
      hInvMassCharmBaryon->Fill(massCharmBaryon);
    }

    // Set production vertex
    // PV
    const KFPVertex kfVertex = createKFPVertexFromCollision(collision);
    const KFParticle kfPv(kfVertex);

    // To V0
    kfPos.SetProductionVertex(kfV0);
    kfNeg.SetProductionVertex(kfV0);

    // To Casc
    KFParticle kfBachToCasc = kfBach;
    KFParticle kfV0ToCasc = kfV0;
    kfBachToCasc.SetProductionVertex(kfCasc);
    kfV0ToCasc.SetProductionVertex(kfCasc);

    // To Charm baryon
    KFParticle kfCascToCharmBaryon = kfCasc;
    KFParticle kfCharmBachToCharmBaryon = kfCharmBach;
    kfCascToCharmBaryon.SetProductionVertex(kfCharmBaryon);
    kfCharmBachToCharmBaryon.SetProductionVertex(kfCharmBaryon);

    // To Pv
    KFParticle kfV0ToPv = kfV0;
    KFParticle kfCascToPv = kfCasc;
    KFParticle kfCharmBachToPv = kfCharmBach;
    KFParticle kfCharmBaryonToPv = kfCharmBaryon;
    kfV0ToPv.SetProductionVertex(kfPv);
    kfCascToPv.SetProductionVertex(kfPv);
    kfCharmBachToPv.SetProductionVertex(kfPv);
    kfCharmBaryonToPv.SetProductionVertex(kfPv);

    //----------Calculate physical quantities -------------
    std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
    std::array<float, 3> svCoord = {kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()};
    std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
    std::array<float, 3> vertexCasc = {kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()};
    std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()};
    std::array<float, 3> pVecCasc = {kfCasc.GetPx(), kfCasc.GetPy(), kfCasc.GetPz()};

    auto trackParCovV0DauPos = getTrackParCovFromKFP(kfPos, convertPDGIntoPID(kfPos.GetPDG()), 1);
    auto trackParCovV0DauNeg = getTrackParCovFromKFP(kfNeg, convertPDGIntoPID(kfNeg.GetPDG()), -1);
    auto trackParCovBach = getTrackParCovFromKFP(kfBachToCasc, convertPDGIntoPID(kfBachToCasc.GetPDG()), bachCharge);
    auto trackParCovCharmBach = getTrackParCovFromKFP(kfCharmBachToCharmBaryon, convertPDGIntoPID(kfCharmBachToCharmBaryon.GetPDG()), -bachCharge);
    auto trackParCovCasc = getTrackParCovFromKFP(kfCascToCharmBaryon, convertPDGIntoPID(kfCascToCharmBaryon.GetPDG()), bachCharge);

    trackParCovV0DauPos.setAbsCharge(1);
    trackParCovV0DauNeg.setAbsCharge(1);
    trackParCovBach.setAbsCharge(1);
    trackParCovCharmBach.setAbsCharge(1);
    trackParCovCasc.setAbsCharge(1);

    auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterV0DauPos, impactParameterV0DauNeg, impactParameterCascBach;
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0DauPos, 2.f, matCorr, &impactParameterV0DauPos);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0DauNeg, 2.f, matCorr, &impactParameterV0DauNeg);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovBach, 2.f, matCorr, &impactParameterCascBach);

    o2::dataformats::DCA impactParameterCasc, impactParameterCharmBach;
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCasc, 2.f, matCorr, &impactParameterCasc);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCharmBach, 2.f, matCorr, &impactParameterCharmBach);

    float kfDcaV0Daughters = kfNeg.GetDistanceFromParticle(kfPos);
    float kfDcaCascDaughters = kfBachToCasc.GetDistanceFromParticle(kfV0ToCasc);
    float kfDcaCharmBaryonDaughters = kfCharmBachToCharmBaryon.GetDistanceFromParticle(kfCascToCharmBaryon);

    // Calculate error of decay length using covariance matrix of PV and SV, same as in DCAFitter
    float covMatrixPv[6];
    kfVertex.GetCovarianceMatrix(covMatrixPv);
    auto* covMatrixSV = kfCharmBaryon.CovarianceMatrix();

    double phi{}, theta{};
    getPointDirection(pvCoord, svCoord, phi, theta);
    float errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPv, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
    float errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPv, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

    // Calculate normalised decay length using KFParticle
    float dlV0{}, errDlV0{};
    kfV0ToPv.GetDecayLength(dlV0, errDlV0);

    float dlCasc{}, errDlCasc{};
    kfCascToPv.GetDecayLength(dlCasc, errDlCasc);

    float dlCharmBaryon{}, dlXYCharmBaryon{}, errDlCharmBaryon{}, errDlXYCharmBaryon{};
    kfCharmBaryonToPv.GetDecayLength(dlCharmBaryon, errDlCharmBaryon);
    kfCharmBaryonToPv.GetDecayLengthXY(dlXYCharmBaryon, errDlXYCharmBaryon);

    float cosPaCascToPv = RecoDecay::cpa(pvCoord, vertexCasc, pVecCasc);
    float cosPaXYCascToPv = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
    float cosPaV0ToPv = RecoDecay::cpa(pvCoord, vertexV0, pVecV0);
    float cosPaXYV0ToPv = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);

    float cosPaCascToCharmBaryon = RecoDecay::cpa(svCoord, vertexCasc, pVecCasc);
    float cosPaXYCascToCharmBaryon = RecoDecay::cpaXY(svCoord, vertexCasc, pVecCasc);
    float cosPaV0ToCasc = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
    float cosPaXYV0ToCasc = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);

    //------------------------------Fill the table------------------------------
    rowCandCascBatch(collision.globalIndex(),
                     primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                     kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ(),
                     kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ(),
                     kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(),
                     errorDecayLength, errorDecayLengthXY,
                     kfCascToCharmBaryon.GetPx(), kfCascToCharmBaryon.GetPy(), kfCascToCharmBaryon.GetPz(),
                     kfCharmBachToCharmBaryon.GetPx(), kfCharmBachToCharmBaryon.GetPy(), kfCharmBachToCharmBaryon.GetPz(),
                     impactParameterCasc.getY(), impactParameterCharmBach.getY(),
                     std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBach.getSigmaY2()),
                     impactParameterCasc.getZ(), impactParameterCharmBach.getZ(),
                     std::sqrt(impactParameterCasc.getSigmaZ2()), std::sqrt(impactParameterCharmBach.getSigmaZ2()),
                     impactParameterV0DauPos.getY(), impactParameterV0DauNeg.getY(), impactParameterCascBach.getY(),
                     impactParameterV0DauPos.getZ(), impactParameterV0DauNeg.getZ(), impactParameterCascBach.getZ(),
                     kfDcaCharmBaryonDaughters, kfDcaCascDaughters, kfDcaV0Daughters,
                     cosPaCascToPv, cosPaXYCascToPv, cosPaV0ToPv, cosPaXYV0ToPv,
                     cosPaCascToCharmBaryon, cosPaXYCascToCharmBaryon, cosPaV0ToCasc, cosPaXYV0ToCasc,
                     massCasc, massLam,
                     kfCasc.GetEta(), kfV0.GetEta(),
                     trackCharmBachelor.globalIndex(), bachTrack.globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                     BIT(decayChannel));

    rowCandCascBatchKf(kfCharmBaryon.Chi2(), kfCasc.Chi2(), kfV0.Chi2(),
                       kfCharmBaryon.GetNDF(), kfCasc.GetNDF(), kfV0.GetNDF(),
                       dlCharmBaryon / errDlCharmBaryon, dlCasc / errDlCasc, dlV0 / errDlV0,
                       kfCharmBaryonToPv.GetChi2(), kfCascToPv.GetChi2(), kfV0ToPv.GetChi2(),
                       kfCharmBaryonToPv.GetNDF(), kfCascToPv.GetNDF(), kfV0ToPv.GetNDF(),
                       kfCascToCharmBaryon.GetChi2(), kfV0ToCasc.GetChi2(),
                       kfCascToCharmBaryon.GetNDF(), kfV0ToCasc.GetNDF());
  }

  // template function for running charm baryon reconstruction
  /// \brief recoMethod is for different reconstruction methods (KFParticle or DCAFitter)
  /// \brief centEstimator is for different centrality estimators
  /// \brief Colls is for collision tables joined with different centraltiy estimators
  template <int recoMethod, hf_centrality::CentralityEstimator centEstimator, typename Colls, typename TLfTracks, typename TTracks>
  void runCreator(Colls const&,
                  FilteredHfCascLf2Prongs const& candidates,
                  aod::Cascades const&,
                  aod::V0s const&,
                  TTracks const& tracks,
                  TLfTracks const& lfTracks,
                  aod::BCsWithTimestamps const&)
  {
    // Loop over candidate
    for (auto const& cand : candidates) {

      // Event selection
      auto collision = cand.collision_as<Colls>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        continue;
      }

      //------------------------------Set Magnetic field------------------------------
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (runNumber != bc.runNumber()) {
        LOG(info) << ">>>>>>>>>> Current run Number : " << runNumber;
        initCCDB(bc, runNumber, ccdb, configs.isRun2 ? configs.ccdbPathGrp : configs.ccdbPathGrpMag, lut, configs.isRun2);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
        LOG(info) << ">>>>>>>>>> Magnetic field: " << magneticField;
        runNumber = bc.runNumber();
      }
      // Magnetic field setting for internal cascade building
      straHelper.fitter.setBz(magneticField);
      // Magnetic field setting for charm baryon building
      if constexpr (recoMethod == aod::hf_cand::VertexerType::DCAFitter) {
        df.setBz(magneticField);
      } else {
        KFParticle::SetField(magneticField);
      }

      // Retrieve tracks for candidate building
      auto cascAodElement = cand.cascade_as<aod::Cascades>();
      auto v0AodElement = cascAodElement.v0_as<aod::V0s>();
      auto posTrack = lfTracks.rawIteratorAt(v0AodElement.posTrackId());
      auto negTrack = lfTracks.rawIteratorAt(v0AodElement.negTrackId());
      auto bachTrack = lfTracks.rawIteratorAt(cascAodElement.bachelorId());
      auto trackCharmBachelor = tracks.rawIteratorAt(cand.prong0Id());

      // Fill candidates before selection
      if (configs.fillHistograms) {
        hCandidateCounterToXiPi->Fill(All);
        hCandidateCounterToOmegaPi->Fill(All);
        hCandidateCounterToOmegaKa->Fill(All);
      }

      // Reconstruct charm-baryon candidate for each decay channel it was skimmed for
      // We need to douplicate the candidates here to account for different mass hypotheses of cascades and bachelors
      if (configs.createToXiPi && TESTBIT(cand.hfflag(), hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi)) {
        if (configs.fillHistograms) {
          hCandidateCounterToXiPi->Fill(HfFlagPass);
        }
        if constexpr (recoMethod == aod::hf_cand::VertexerType::DCAFitter) {
          buildCandidateWithDCAFitter<hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collision, posTrack, negTrack, bachTrack, trackCharmBachelor, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
        } else {
          buildCandidateWithKFParticle<hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collision, posTrack, negTrack, bachTrack, trackCharmBachelor, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
        }
      }

      if (configs.createOmegaCToOmegaPi && TESTBIT(cand.hfflag(), hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
        if (configs.fillHistograms) {
          hCandidateCounterToOmegaPi->Fill(HfFlagPass);
        }
        if constexpr (recoMethod == aod::hf_cand::VertexerType::DCAFitter) {
          buildCandidateWithDCAFitter<hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collision, posTrack, negTrack, bachTrack, trackCharmBachelor, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
        } else {
          buildCandidateWithKFParticle<hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collision, posTrack, negTrack, bachTrack, trackCharmBachelor, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
        }
      }

      if (configs.createOmegaCToOmegaK && TESTBIT(cand.hfflag(), hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK)) {
        if (configs.fillHistograms) {
          hCandidateCounterToOmegaKa->Fill(HfFlagPass);
        }
        if constexpr (recoMethod == aod::hf_cand::VertexerType::DCAFitter) {
          buildCandidateWithDCAFitter<hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collision, posTrack, negTrack, bachTrack, trackCharmBachelor, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
        } else {
          buildCandidateWithKFParticle<hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collision, posTrack, negTrack, bachTrack, trackCharmBachelor, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
        }
      }
    } // candidate loop
  } // end of run function

  // Process functions with DCAFitter
  void processWithDCAFitterNoCent(SelectedCollisions const& collisions,
                                  FilteredHfCascLf2Prongs const& candidates,
                                  aod::Cascades const& cascades,
                                  aod::V0s const& v0s,
                                  aod::TracksWCovDca const& tracks,
                                  TracksWCovIU const& lfTracks,
                                  aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::DCAFitter, CentralityEstimator::None>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithDCAFitterNoCent, "Charm-baryon candidate reconstruction with DcaFitter method, no centrality", true);

// TODO: add support for tracked cascades
#if 0
  void processWithDCAFitterNoCentWithTrackedCasc(SelectedCollisions const& collisions,
                                                 FilteredHfCascLf2Prongs const& candidates,
                                                 aod::TrackedCascades const& cascades,
                                                 aod::V0s const& v0s,
                                                 TrackedCascLinked const& trackedCascLinked,
                                                 TracksWCovDcaExtraPidPrPiKa const& tracks,
                                                 aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::DCAFitter, CentralityEstimator::None>(collisions, candidates, trackedCascFull, trackedCascLinked, tracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithDCAFitterNoCentWithTrackedCasc, "Charm-baryon candidate reconstruction with DcaFitter method with tracked cascade, no centrality", false);
#endif

  void processWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                    FilteredHfCascLf2Prongs const& candidates,
                                    aod::Cascades const& cascades,
                                    aod::V0s const& v0s,
                                    aod::TracksWCovDca const& tracks,
                                    TracksWCovIU const& lfTracks,
                                    aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::DCAFitter, CentralityEstimator::FT0C>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithDCAFitterCentFT0C, "Charm-baryon candidate reconstruction with DcaFitter method, centrality selection on FT0C", false);

  void processWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                    FilteredHfCascLf2Prongs const& candidates,
                                    aod::Cascades const& cascades,
                                    aod::V0s const& v0s,
                                    aod::TracksWCovDca const& tracks,
                                    TracksWCovIU const& lfTracks,
                                    aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::DCAFitter, CentralityEstimator::FT0M>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithDCAFitterCentFT0M, "Charm-baryon candidate reconstruction with DcaFitter method, centrality selection on FT0M", false);

  // Process functions with KFParticle
  void processWithKFParticleNoCent(SelectedCollisions const& collisions,
                                   FilteredHfCascLf2Prongs const& candidates,
                                   aod::Cascades const& cascades,
                                   aod::V0s const& v0s,
                                   TracksWCovDcaExtraPidPrPiKa const& tracks,
                                   TracksWCovExtraPidIU const& lfTracks,
                                   aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::KfParticle, CentralityEstimator::None>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithKFParticleNoCent, "Charm-baryon candidate reconstruction with KFParticle method, no centrality", false);

  void processWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                     FilteredHfCascLf2Prongs const& candidates,
                                     aod::Cascades const& cascades,
                                     aod::V0s const& v0s,
                                     TracksWCovDcaExtraPidPrPiKa const& tracks,
                                     TracksWCovExtraPidIU const& lfTracks,
                                     aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::KfParticle, CentralityEstimator::FT0C>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithKFParticleCentFT0C, "Charm-baryon candidate reconstruction with KFParticle method, centrality on FT0C", false);

  void processWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                     FilteredHfCascLf2Prongs const& candidates,
                                     aod::Cascades const& cascades,
                                     aod::V0s const& v0s,
                                     TracksWCovDcaExtraPidPrPiKa const& tracks,
                                     TracksWCovExtraPidIU const& lfTracks,
                                     aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreator<aod::hf_cand::VertexerType::KfParticle, CentralityEstimator::FT0M>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processWithKFParticleCentFT0M, "Charm-baryon candidate reconstruction with KFParticle method, centrality on FT0M", false);

  // Process functions for Collision monitoring
  void processCollisionsNoCent(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults> const& collisions,
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
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processCollisionsNoCent, "Collision monitoring - No Centrality", true);

  void processCollisionsCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                 aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {

      // bitmask with event selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      // monitor the satisfied event selection
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processCollisionsCentFT0C, "Collision monitoring - Centrality selection with FT0C", false);

  void processCollisionsCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                 aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {

      // bitmask with event selection info
      float centrality{-1.f};
      float occupancy = getOccupancyColl(collision, OccupancyEstimator::Its);
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0M, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);

      // monitor the satisfied event selection
      hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processCollisionsCentFT0M, "Collision monitoring - Centrality selection with FT0M", false);
};

struct HfCandidateCreatorXic0Omegac0QaMc {
  Spawns<HfCandCascBachExtQa> rowCandidateCascBach;
  Produces<HfCandCascBachMCRecQa> rowMcMatchRecCascBach;
  Produces<HfCandCascBachMCGenQa> rowMcMatchGenCascBach;

  // Configurables
  struct : ConfigurableGroup {
    Configurable<bool> rejectBackground{"rejectBackground", false, "Reject particles from background events"};
    Configurable<bool> matchKinkedDecayTopology{"matchKinkedDecayTopology", true, "Match also candidates with tracks that decay with kinked topology"};
    Configurable<bool> matchInteractionsWithMaterial{"matchInteractionsWithMaterial", true, "Match also candidates with tracks that interact with material"};
  } configs;

  std::array<int, nMcDecayChannels> pdgOfCharmBaryon{+kXiC0, +kXiC0, +kOmegaC0, +kOmegaC0, +kOmegaC0};
  std::array<int, nMcDecayChannels> pdgOfCascade{+kXiMinus, +kOmegaMinus, +kOmegaMinus, +kOmegaMinus, +kXiMinus};
  std::array<int, nMcDecayChannels> pdgOfCharmBachelor{+kPiPlus, +kKPlus, +kPiPlus, +kKPlus, +kPiPlus};
  std::array<int, nMcDecayChannels> pdgOfBachelor{+kPiMinus, +kKMinus, +kKMinus, +kKMinus, +kPiMinus};

  std::array<int8_t, nMcDecayChannels> decayChannelFlags{hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Xic0ToXiPi,
                                                         hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Xic0ToOmegaK,
                                                         hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToOmegaPi,
                                                         hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToOmegaK,
                                                         hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToXiPi};

  // Table aliases
  using TracksWMcIU = soa::Join<aod::TracksIU, McTrackLabels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;
  using McCollisionsCentFT0Cs = soa::Join<aod::McCollisions, aod::McCentFT0Cs, aod::MultsExtraMC>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<McCollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId; // -> Why use unsorted??
  PresliceUnsorted<McCollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<McCollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;

  HistogramRegistry registry{"registry"};
  HfEventSelectionMc hfEvSelMc;

  void init(InitContext& initContext)
  {
    const std::array<bool, 3> mcProcesses{doprocessMcNoCent, doprocessMcCentFT0C, doprocessMcCentFT0M};
    if (std::accumulate(mcProcesses.begin(), mcProcesses.end(), 0) > 1) {
      LOGP(fatal, "At most one MC process function can be enabled at a time.");
    }

    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-candidate-creator-xic0-omegac0-qa") {
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  template <typename Track, std::size_t N>
  int matchDecay(aod::McParticles const& mcParticles, std::array<Track, N> const& daughters, int pdgMother, std::array<int, N> const& pdgDaughters,
                 int8_t* matchSign, int depth, int8_t* nPiToMu = nullptr, int8_t* nKaToPi = nullptr, int8_t* nMaterial = nullptr)
  {
    if (configs.matchKinkedDecayTopology && configs.matchInteractionsWithMaterial) {
      return RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles, daughters, pdgMother, pdgDaughters, true, matchSign, depth, nPiToMu, nKaToPi, nMaterial);
    }
    if (configs.matchKinkedDecayTopology) {
      return RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles, daughters, pdgMother, pdgDaughters, true, matchSign, depth, nPiToMu, nKaToPi);
    }
    if (configs.matchInteractionsWithMaterial) {
      return RecoDecay::getMatchedMCRec<false, true, false, false, true>(mcParticles, daughters, pdgMother, pdgDaughters, true, matchSign, depth, nullptr, nullptr, nMaterial);
    }
    return RecoDecay::getMatchedMCRec<false, true>(mcParticles, daughters, pdgMother, pdgDaughters, true, matchSign, depth);
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, typename Colls, typename McCollisions>
  void runXic0Omegac0Mc(TracksWMcIU const& tracks,
                        aod::McParticles const& mcParticles,
                        Colls const& collsWithMcLabels,
                        McCollisions const& mcCollisions,
                        BCsInfo const&)
  {
    rowCandidateCascBach->bindExternalIndices(&tracks);

    // Match reconstructed candidates
    for (const auto& candidate : *rowCandidateCascBach) {
      int indexRecCharmBaryon = -1;
      int indexRecCascade = -1;
      int indexRecV0 = -1;
      int8_t sign = 0;
      int8_t signCasc = 0;
      int8_t signV0 = 0;
      int8_t flagChannelMain = 0;
      int8_t flagChannelResonant = 0; // TODO: Implement resonant channel matching
      int8_t origin = 0;
      int8_t nKinkedTracks = 0;
      int8_t nInteractionsWithMaterial = 0;

      std::vector<int> idxBhadMothers{};

      auto arrayDaughters = std::array{candidate.template prong1_as<TracksWMcIU>(),
                                       candidate.template bachelor_as<TracksWMcIU>(),
                                       candidate.template posTrack_as<TracksWMcIU>(),
                                       candidate.template negTrack_as<TracksWMcIU>()};
      auto arrayDaughtersCasc = std::array{candidate.template bachelor_as<TracksWMcIU>(),
                                           candidate.template posTrack_as<TracksWMcIU>(),
                                           candidate.template negTrack_as<TracksWMcIU>()};
      auto arrayDaughtersV0 = std::array{candidate.template posTrack_as<TracksWMcIU>(),
                                         candidate.template negTrack_as<TracksWMcIU>()};

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
        }
        if (fromBkg) {
          rowMcMatchRecCascBach(flagChannelMain, origin, flagChannelResonant, -1.f, 0, 0, 0);
          continue;
        }
      }

      for (int decayChannel = 0; decayChannel < nMcDecayChannels; ++decayChannel) {
        int8_t nPiToMu = 0;
        int8_t nKaToPi = 0;
        int8_t nMaterial = 0;
        // CharmBaryon -> Charm bachelor + Cascade
        indexRecCharmBaryon = matchDecay(mcParticles, arrayDaughters, pdgOfCharmBaryon[decayChannel], std::array{pdgOfCharmBachelor[decayChannel], pdgOfBachelor[decayChannel], +kProton, +kPiMinus}, &sign, 3, &nPiToMu, &nKaToPi, &nMaterial);

        if (indexRecCharmBaryon > -1) {
          // Cascade -> Bachelor + V0
          indexRecCascade = matchDecay(mcParticles, arrayDaughtersCasc, pdgOfCascade[decayChannel], std::array{pdgOfBachelor[decayChannel], +kProton, +kPiMinus}, &signCasc, 2);

          if (indexRecCascade > -1) {
            // V0 -> Pos + Neg
            indexRecV0 = matchDecay(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, &signV0, 1);

            if (indexRecV0 > -1) {
              flagChannelMain = sign * decayChannelFlags[decayChannel];
              nKinkedTracks = nPiToMu + nKaToPi;
              nInteractionsWithMaterial = nMaterial;
              break;
            }
          }
        }
      }

      // Check if charm baryon is from b-hadron decay (prompt vs non-prompt)
      if (flagChannelMain != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        rowMcMatchRecCascBach(flagChannelMain, origin, flagChannelResonant, bHadMother.pt(), bHadMother.pdgCode(), nKinkedTracks, nInteractionsWithMaterial);
      } else {
        rowMcMatchRecCascBach(flagChannelMain, origin, flagChannelResonant, -1.f, 0, nKinkedTracks, nInteractionsWithMaterial);
      }
    } // candidate loop

    // Match generated particles
    for (const auto& mcCollision : mcCollisions) {
      const auto& mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());

      float centrality{-1.f};
      o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
      int nSplitColl{0};

      if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::None) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
        nSplitColl = collSlice.size();
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
        nSplitColl = collSlice.size();
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
        nSplitColl = collSlice.size();
      }

      hfEvSelMc.fillHistograms<centEstimator>(mcCollision, rejectionMask, nSplitColl);

      if (rejectionMask != 0) {
        for (unsigned int i = 0; i < mcParticlesPerMcColl.size(); ++i) {
          rowMcMatchGenCascBach(0, 0, 0, -1);
        }
        continue;
      }

      // Match generated particles
      for (auto const& particle : mcParticlesPerMcColl) {
        int8_t flagChannelMain = 0;
        int8_t flagChannelResonant = 0; // TODO: Implement resonant channel matching
        int8_t origin = RecoDecay::OriginType::None;
        int8_t sign = 0;

        std::vector<int> idxBhadMothers{};

        // Reject particles from background events
        if (particle.fromBackgroundEvent() && configs.rejectBackground) {
          rowMcMatchGenCascBach(flagChannelMain, origin, flagChannelResonant, -1);
          continue;
        }

        for (int decayChannel = 0; decayChannel < nMcDecayChannels; ++decayChannel) {
          // Charm Baryon -> Cascade + Charm bachelor
          if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgOfCharmBaryon[decayChannel], std::array{pdgOfCascade[decayChannel], pdgOfCharmBachelor[decayChannel]}, true, &sign)) {
            for (auto const& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
              if (std::abs(daughterCharm.pdgCode()) != pdgOfCascade[decayChannel]) {
                continue;
              }

              // Cascade -> Lambda + pi
              if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgOfCascade[decayChannel], std::array{+kLambda0, pdgOfBachelor[decayChannel]}, true)) {
                for (auto const& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                  if (std::abs(daughterCascade.pdgCode()) != +kLambda0) {
                    continue;
                  }

                  // Lambda -> p + pi
                  if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
                    flagChannelMain = sign * decayChannelFlags[decayChannel];
                    break;
                  }
                } // V0 daughter loop
              } // cascade daughter loop
            } // charm daughter loop
          }
        }

        // Check if charm is prompt or non-prompt
        if (flagChannelMain != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (origin == RecoDecay::OriginType::NonPrompt) {
          rowMcMatchGenCascBach(flagChannelMain, origin, flagChannelResonant, idxBhadMothers[0]);
        } else {
          rowMcMatchGenCascBach(flagChannelMain, origin, flagChannelResonant, -1);
        }

      } // particle loop
    } // end of collision loop
  }

  // void processMcEmpty(aod::Collisions const&)
  // {
  // }
  // PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcEmpty, "Empty process function to prevent workflow from getting stuck", true);

  void processMcNoCent(TracksWMcIU const& tracks,
                       aod::McParticles const& mcParticles,
                       soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                       McCollisionsNoCents const& collsWithMcLabels,
                       BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None>(tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcNoCent, "Process MC - no centrality", false);

  void processMcCentFT0C(TracksWMcIU const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsCentFT0Cs const& mcCollisions,
                         McCollisionsFT0Cs const& collsWithMcLabels,
                         BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C>(tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcCentFT0C, "Process MC - centrality with FT0C", false);

  void processMcCentFT0M(TracksWMcIU const& tracks,
                         aod::McParticles const& mcParticles,
                         McCollisionsCentFT0Ms const& mcCollisions,
                         McCollisionsFT0Ms const& collsWithMcLabels,
                         BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0M>(tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcCentFT0M, "Process MC - centrality with FT0M", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0Qa>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0QaMc>(cfgc)};
}
