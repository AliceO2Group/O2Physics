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

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTablesQa.h"
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

#include <array>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <stdexcept>
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

struct HfCandidateCreatorXic0Omegac0Qa {
  Produces<aod::HfCandCascBachQa> rowCandCascBatch;
  // Additional information for candidates created with KFParticle
  Produces<aod::HfCandCascBachKfQa> rowCandCascBatchKf;

  // Configurables
  struct : ConfigurableGroup {
    // Options for internal cascade building - DCAFitter settings
    // ...Initial vaules taken from PWGLF/Utils/strangenessBuilderHelper
    // ----------------------------------------------------------------
    Configurable<bool> propagateToPCALF{"propagateToPCALF", true, "Create tracks version propagated to PCA"};
    Configurable<double> maxRLF{"maxRLF", 200., "Reject PCA's above this radius"};
    Configurable<double> maxDZIniLF{"maxDZIniLF", 1e9, "Reject (if>0) PCA candidate if tracks DZ exceeds this threshold"};
    Configurable<double> maxDXYIniLF{"maxDXYIniLF", 4.0f, "Reject (if>0) PCA candidate if tracks DXY exceeds this threshold"};
    Configurable<double> minParamChangeLF{"minParamChangeLF", 1.e-3, "Stop iteration if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2ChangeLF{"minRelChi2ChangeLF", 0.9, "Stop iteration if Chi2/Chi2old > this"};
    Configurable<double> maxChi2LF{"maxChi2LF", 1e9, "Discard vertices with Chi2/Nprongs > this(or sum {DCAi^2}/Nprongs for abs. distance minimization)"};
    Configurable<bool> useAbsDCALF{"useAbsDCALF", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCALF{"useWeightedFinalPCALF", false, "Recalculate vertex position using track covariance, effective only if useAbsDCA is true"};
    Configurable<bool> refitWithMaterialCorrectionLF{"refitWithMaterialCorrectionLF", false, "Do refit after material correction applied"};
    // Options for internal V0 building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // ...Modified according to the configurable in core wagon
    // ---------------------------------------------------------------------
    Configurable<int> minCrossedRowsFromLF{"minCrossedRowsFromLF", 50, "minimun TPC crossed rows for daughter tracks. Used for internal V0 Building"};
    Configurable<float> dcanegtopvFromLF{"dcanegtopvFromLF", .1, "DCV Neg to PV"};
    Configurable<float> dcapostopvFromLF{"dcapostopvFromLF", .1, "DCV Pos To PV"};
    Configurable<double> v0cospaFromLF{"v0cospaFromLF", 0.95, "V0 CosPA"};
    Configurable<float> dcav0dauFromLF{"dcav0dauFromLF", 1.0, "DCA V0 Daughters"};
    Configurable<float> v0radiusFromLF{"v0radiusFromLF", 0.9, "v0radius"};
    Configurable<float> maxDaughterEtaFromLF{"maxDaughterEtaFromLF", 5.0, "Maximun daughter eta (in abs value)"};
    // Options for internal cascade building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // ...Modified according to the configurable in core wagon
    // --------------------------------------------------------------------
    Configurable<float> dcabachtopvFromLF{"dcabachtopvFromLF", .05, "DCV Bach to PV"};
    Configurable<float> cascradiusFromLF{"cascradiusFromLF", .9, "DCV Bach to PV"};
    Configurable<float> casccospaFromLF{"casccospaFromLF", 0.95, "Cascade CosPA"};
    Configurable<float> dcacascdauFromLF{"dcacascdauFromLF", 1.0, "DCA cascade daughters"};
    Configurable<float> lambdaMassWindowFromLF{"lambdaMassWindowFromLF", 0.01, "Distance from Lambda mass(does not apply to KF path)"};
    // Options for internal cascade building - KF Building specifics
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // --------------------------------------------------------------------
    Configurable<bool> kfTuneForOmegaFromLF{"kfTuneForOmegaFromLF", false, "if enabled, take main cascade properties from omega fit instread of Xi fit(=default)"};
    Configurable<int> kfConstructMethodFromLF{"kfConstructMethodFromLF", 2, "2 : Daughter particle masses stay fixed in construction process"};
    Configurable<bool> kfUseV0MassConstraintFromLF{"kfUseV0MassConstraintFromLF", true, "KF : Use Lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraintFromLF{"kfUseCascadeMassConstraintFromLF", false, "KF : Use Cascade mass constraint - WARNING : Not adequate for inv mass analysis of Xi"};
    Configurable<bool> kfDoDCAFitterPreMinimV0FromLF{"kfDoDCAFitterPreMinimV0FromLF", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for V0"};
    Configurable<bool> kfDoDCAFitterPreMinimCascFromLF{"kfDoDCAFitterPreMinimCascFromLF", false, "KF : do DCAFitter pre-optimization before KF fit to include material correction for Xi"};
  } LFConfigs;

  struct : ConfigurableGroup {
    // Switch for filling histograms
    // -----------------------------
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill validation plots"};
    // Magnetic field setting from CCDB
    // --------------------------------
    Configurable<bool> isRun2{"isRun2", false, "enable Run2 or Run3 GRP objects for magnetic field"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "https://alice-ccdb.cern.ch", "url of the ccdb object"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parameterization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "path of the group file (Run2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run3)"};
    // Cascade pre selection
    // ---------------------
    Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
    Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascades"};
    Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in XY plane"};
    Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 1.0, "Max DCA of V0 daughter"};
    Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 1.0, "Max DCA of cascade daughter"};
    // Options for DCAFitter
    // ---------------------
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "Create tracks version propagated to PCA"};
    Configurable<double> maxR{"maxR", 200., "Reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "Reject (if>0) PCA candidate if tracks DZ exceeds this threshold"};
    Configurable<double> maxDXYIni{"maxDXYIni", 4., "Reject (if>0) PCA candidate if tracks DXY exceeds this threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "Stop iteration if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "Stop iteration if Chi2/Chi2old > this"};
    Configurable<double> maxChi2{"maxChi2", 100, "Discard vertices with Chi2/Nprongs > this(or sum {DCAi^2}/Nprongs for abs. distance minimization)"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariance, effective only if useAbsDCA is true"};
    // Options for KFParticle
    // ----------------------
    // V0 cuts
    Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.0075, "Distance from LambdaMass"};
    // For KF Particle operation
    Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct method"};
    Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", false, "KF: Use lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF: Use lambda mass constraint"};
    // Options for QA histogram binning
    // -----------------------------
    // For Cascade
    Configurable<int> nBinMassCasc{"nBinMassCasc", 1000, "nBinCascMass"};
    Configurable<float> minMassCasc{"minMassCasc", 1.0, "xiMassMin"};
    Configurable<float> maxMassCasc{"maxMassCasc", 2.0, "xiMassMax"};
    Configurable<int> nBinPtCasc{"nBinPtCasc", 100, "nBinPtXi"};
    Configurable<float> minPtCasc{"minPtCasc", 0.0, "minimun value of cascade"};
    Configurable<float> maxPtCasc{"maxPtCasc", 20.0, "maximum value of cascade"};
    // For Prong0
    Configurable<int> nBinMassProng0{"nBinMassProng0", 100, "nBinMassPi"};
    Configurable<float> minMassProng0{"minMassProng0", 0.0, "ptPiMin"};
    Configurable<float> maxMassProng0{"maxMassProng0", 20.0, "ptPiMax"};
    Configurable<int> nBinPtProng0{"nBinPtProng0", 100, "nBinPtProng0"};
    Configurable<float> minPtProng0{"minPtProng0", 0.0, "ptPiMin"};
    Configurable<float> maxPtProng0{"maxPtProng0", 20.0, "ptPiMax"};
    // For Charm Baryon
    Configurable<int> nBinMassCharmBaryon{"nBinMassCharmBaryon", 3000, "nBinXic0Mass"};
    Configurable<float> minMassCharmBaryon{"minMassCharmBaryon", 1.0, "xic0MassMin"};
    Configurable<float> maxMassCharmBaryon{"maxMassCharmBaryon", 4.0, "xic0MassMax"};
    Configurable<int> nBinPtCharmBaryon{"nBinPtCharmBaryon", 100, "nBinXic0Pt"};
    Configurable<float> minPtCharmBaryon{"minPtCharmBaryon", 0.0, "xic0PtMin"};
    Configurable<float> maxPtCharmBaryon{"maxPtCharmBaryon", 20.0, "xic0PtMax"};
    // Etc
    Configurable<int> nBinCpa2Prong{"nBinCpa2Prong", 240, "nBinCpa2Prong"};
    Configurable<float> cpa2ProngMin{"cpa2ProngMin", -1.2, "cpa2ProngMin"};
    Configurable<float> cpa2ProngMax{"cpa2ProngMax", 1.2, "cpa2ProngMax"};
    Configurable<int> nBinImpParXYXi{"nBinImpParXYXi", 30, "nBinImpParXYXi"};
    Configurable<float> impParXYXiMin{"impParXYXiMin", -1.5, "impParXYXiMin"};
    Configurable<float> impParXYXiMax{"impParXYXiMax", 1.5, "impParXYXiMax"};
    Configurable<int> nBinImpParXYPi{"nBinImpParXYPi", 30, "nBinImpParXYPi"};
    Configurable<float> impParXYPiMin{"impParXYPiMin", -1.5, "impParXYPiMin"};
    Configurable<float> impParXYPiMax{"impParXYPiMax", 1.5, "impParXYPiMax"};
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

  enum CharmBaryonCandCounter { All = 0,
                                HfFlagPass,
                                CascReconstructed,
                                VertexFit };

  // Table aliases
  using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
  using TracksWCovIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
  using TracksWCovDcaExtraPidPrPiKa = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPr, aod::TracksPidPi, aod::TracksPidKa>;
  using TracksWCovExtraPidIU = soa::Join<aod::TracksIU, TracksCovIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  HistogramRegistry registry{"hists"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  HfEventSelection hfEvSel;

  // PDG Id of daughter tracks & V0s & cascades & charm baryons - Used in KFParticle
  int pdgIdOfV0DauPos{}, pdgIdOfV0DauNeg{}, pdgIdOfBach{}, pdgIdOfBachRej{}, pdgIdOfCharmBach{};
  int pdgIdOfV0{}, pdgIdOfCascade{}, pdgIdOfCharmBaryon{};

  // Track PID: PID value of tracks defined under o2::track::PID namespace
  int trackPidOfCascade{};

  // Mass of daughter tracks & V0s & cascades & charm baryons;
  float massOfV0DauPos{}, massOfV0DauNeg{}, massOfCharmBach{};
  float massOfV0{}, massOfCascade{};

  // Pointer of histograms for QA
  std::shared_ptr<TH1> hInvMassCharmBaryonToXiPi, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaKa;
  std::shared_ptr<TH1> hCandidateCounterToXiPi, hCandidateCounterToOmegaPi, hCandidateCounterToOmegaKa;

  void init(InitContext const&)
  {
    std::vector<bool> processesToXiPiDca{doprocessToXiPiWithDCAFitterNoCent, /*doprocessToXiPiWithDCAFitterNoCentWithTrackedCasc,*/ doprocessToXiPiWithDCAFitterCentFT0C, doprocessToXiPiWithDCAFitterCentFT0M};
    std::vector<bool> processesToOmegaPiDca{doprocessToOmegaPiWithDCAFitterNoCent, doprocessToOmegaPiWithDCAFitterCentFT0C, doprocessToOmegaPiWithDCAFitterCentFT0M};
    std::vector<bool> processesToOmegaKaDca{doprocessToOmegaKaWithDCAFitterNoCent, doprocessToOmegaKaWithDCAFitterCentFT0C, doprocessToOmegaKaWithDCAFitterCentFT0M};
    std::vector<bool> processesToXiPiKf{doprocessToXiPiWithKFParticleNoCent, doprocessToXiPiWithKFParticleCentFT0C, doprocessToXiPiWithKFParticleCentFT0M};
    std::vector<bool> processesToOmegaPiKf{doprocessToOmegaPiWithKFParticleNoCent, doprocessToOmegaPiWithKFParticleCentFT0C, doprocessToOmegaPiWithKFParticleCentFT0M};
    std::vector<bool> processesToOmegaKaKf{doprocessToOmegaKaWithKFParticleNoCent, doprocessToOmegaKaWithKFParticleCentFT0C, doprocessToOmegaKaWithKFParticleCentFT0M};
    std::vector<bool> processesCollMonitoring{doprocessCollisionsNoCent, doprocessCollisionsCentFT0C, doprocessCollisionsCentFT0M};

    int xipiEnabledDca = std::accumulate(processesToXiPiDca.begin(), processesToXiPiDca.end(), 0);
    int xipiEnabledKf = std::accumulate(processesToXiPiKf.begin(), processesToXiPiKf.end(), 0);
    int omegapiEnabledDca = std::accumulate(processesToOmegaPiDca.begin(), processesToOmegaPiDca.end(), 0);
    int omegapiEnabledKf = std::accumulate(processesToOmegaPiKf.begin(), processesToOmegaPiKf.end(), 0);
    int omegakaEnabledDca = std::accumulate(processesToOmegaKaDca.begin(), processesToOmegaKaDca.end(), 0);
    int omegakaEnabledKf = std::accumulate(processesToOmegaKaKf.begin(), processesToOmegaKaKf.end(), 0);

    // Exit if workflow is not configured correctly - More than one process function enabled for candidate to XiPi
    if ((xipiEnabledDca > 0) && (xipiEnabledKf > 0)) {
      LOGP(fatal, "More than one process function enabled for candidte decaying to xi pi");
    }

    // Exit if workflow is not configured correctly - More than one process function enabled for candidate to OmegaPi
    if ((omegapiEnabledDca > 0) && (omegapiEnabledKf > 0)) {
      LOGP(fatal, "More than one process function enabled for candidte decaying to omega pi");
    }

    // Exit if workflow is not configured correctly - More than one process function enabled for candidate to OmegaKa
    if ((omegakaEnabledDca > 0) && (omegakaEnabledKf > 0)) {
      LOGP(fatal, "More than one process function enabled for candidte decaying to omega ka");
    }

    // Exit if workflow is not configured correctly - More than one process enabled for collision monitoring
    if (std::accumulate(processesCollMonitoring.begin(), processesCollMonitoring.end(), 0) > 1) {
      LOGP(fatal, "More than one process fucntion for CollMonitoring was enabled. Please choose only one process function");
    }

    // Assign pdg & mass hypothesis for each decay channel
    if (xipiEnabledDca > 0 || xipiEnabledKf > 0) {
      pdgIdOfV0DauPos = kProton;
      pdgIdOfV0DauNeg = kPiMinus;
      pdgIdOfBach = kPiMinus;
      pdgIdOfBachRej = 0; // -> Not used in XiPi channel
      pdgIdOfCharmBach = kPiPlus;

      pdgIdOfV0 = kLambda0;
      pdgIdOfCascade = kXiMinus;
      pdgIdOfCharmBaryon = kXiC0;

      trackPidOfCascade = o2::track::PID::XiMinus;

      massOfCharmBach = o2::constants::physics::MassPiPlus;
      massOfV0DauPos = o2::constants::physics::MassProton;
      massOfV0DauNeg = o2::constants::physics::MassPionCharged;
      massOfV0 = o2::constants::physics::MassLambda;
      massOfCascade = o2::constants::physics::MassXiMinus;
    } else if (omegapiEnabledDca > 0 || omegapiEnabledKf > 0) {
      pdgIdOfV0DauPos = kProton;
      pdgIdOfV0DauNeg = kPiMinus;
      pdgIdOfBach = kKMinus;
      pdgIdOfBachRej = kPiMinus;
      pdgIdOfCharmBach = kPiPlus;

      pdgIdOfV0 = kLambda0;
      pdgIdOfCascade = kOmegaMinus;
      pdgIdOfCharmBaryon = kOmegaC0;

      trackPidOfCascade = o2::track::PID::OmegaMinus;

      massOfCharmBach = o2::constants::physics::MassPiPlus;
      massOfV0DauPos = o2::constants::physics::MassProton;
      massOfV0DauNeg = o2::constants::physics::MassPionCharged;
      massOfV0 = o2::constants::physics::MassLambda;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
    } else if (omegakaEnabledDca > 0 || omegakaEnabledKf > 0) {
      pdgIdOfV0DauPos = kProton;
      pdgIdOfV0DauNeg = kPiMinus;
      pdgIdOfBach = kKMinus;
      pdgIdOfBachRej = kPiMinus;
      pdgIdOfCharmBach = kKPlus;

      pdgIdOfV0 = kLambda0;
      pdgIdOfCascade = kOmegaMinus;
      pdgIdOfCharmBaryon = kOmegaC0;

      trackPidOfCascade = o2::track::PID::OmegaMinus;

      massOfCharmBach = o2::constants::physics::MassKPlus;
      massOfV0DauPos = o2::constants::physics::MassProton;
      massOfV0DauNeg = o2::constants::physics::MassPionCharged;
      massOfV0 = o2::constants::physics::MassLambda;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
    }
    LOGF(info, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    LOGF(info, "  Chosen reconstruction info: ToXiPi(%d)/ ToOmegaPi(%d)/ ToOmegaKa(%d)", xipiEnabledDca + xipiEnabledKf, omegapiEnabledDca + omegapiEnabledKf, omegakaEnabledDca + omegakaEnabledKf);
    LOGF(info, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    LOGF(info, "  PDG ID of V0 positive daughter: %d", pdgIdOfV0DauPos);
    LOGF(info, "  PDG ID of V0 negative daughter: %d", pdgIdOfV0DauNeg);
    LOGF(info, "  PDG ID of V0: %d", pdgIdOfV0);
    LOGF(info, "  PDG ID of Bachelor: %d", pdgIdOfBach);
    LOGF(info, "  PDG ID of Bachelor Rej: %d", pdgIdOfBachRej);
    LOGF(info, "  PDG ID of Cascade: %d", pdgIdOfCascade);
    LOGF(info, "  PDG ID of Charm Bachelor: %d", pdgIdOfCharmBach);
    LOGF(info, "  PDG ID of Charm Baryon: %d", pdgIdOfCharmBaryon);
    LOGF(info, "-----------------------------------------------------------------------");
    LOGF(info, "  Mass of V0 set as: %f", massOfV0);
    LOGF(info, "  Mass of CharmBach set as: %f", massOfCharmBach);
    LOGF(info, "  Mass of Casc as: %f", massOfCascade);
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
    // --------------------------------------------

    // Settings for internal V0 building
    straHelper.v0selections.minCrossedRows = LFConfigs.minCrossedRowsFromLF;
    straHelper.v0selections.dcanegtopv = LFConfigs.dcanegtopvFromLF;
    straHelper.v0selections.dcapostopv = LFConfigs.dcapostopvFromLF;
    straHelper.v0selections.v0cospa = LFConfigs.v0cospaFromLF;
    straHelper.v0selections.dcav0dau = LFConfigs.dcav0dauFromLF;
    straHelper.v0selections.v0radius = LFConfigs.v0radiusFromLF;
    straHelper.v0selections.maxDaughterEta = LFConfigs.maxDaughterEtaFromLF;

    // Settings for internal Cascade building
    straHelper.cascadeselections.minCrossedRows = LFConfigs.minCrossedRowsFromLF;
    straHelper.cascadeselections.dcabachtopv = LFConfigs.dcabachtopvFromLF;
    straHelper.cascadeselections.cascradius = LFConfigs.cascradiusFromLF;
    straHelper.cascadeselections.casccospa = LFConfigs.casccospaFromLF;
    straHelper.cascadeselections.dcacascdau = LFConfigs.dcacascdauFromLF;
    straHelper.cascadeselections.lambdaMassWindow = LFConfigs.lambdaMassWindowFromLF;
    straHelper.cascadeselections.maxDaughterEta = LFConfigs.maxDaughterEtaFromLF;

    // Fitter setting
    straHelper.fitter.setPropagateToPCA(LFConfigs.propagateToPCALF);
    straHelper.fitter.setMaxR(LFConfigs.maxRLF);
    straHelper.fitter.setMaxDZIni(LFConfigs.maxDZIniLF);
    straHelper.fitter.setMaxDXYIni(LFConfigs.maxDXYIniLF);
    straHelper.fitter.setMinParamChange(LFConfigs.minParamChangeLF);
    straHelper.fitter.setMinRelChi2Change(LFConfigs.minRelChi2ChangeLF);
    straHelper.fitter.setMaxChi2(LFConfigs.maxChi2LF);
    straHelper.fitter.setUseAbsDCA(LFConfigs.useAbsDCALF);
    straHelper.fitter.setWeightedFinalPCA(LFConfigs.useWeightedFinalPCALF);
    straHelper.fitter.setRefitWithMatCorr(LFConfigs.refitWithMaterialCorrectionLF);

    // Extra initialization for DCAFitter
    // ----------------------------------
    if (xipiEnabledDca == 1 || omegapiEnabledDca == 1 || omegakaEnabledDca == 1) {
      registry.get<TH1>(HIST("hVertexerType"))->Fill(aod::hf_cand::VertexerType::DCAFitter);
      df.setPropagateToPCA(configs.propagateToPCA);
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
    // -----------------------------------
    if (xipiEnabledKf == 1 || omegapiEnabledKf == 1 || omegakaEnabledKf == 1) {
      registry.get<TH1>(HIST("hVertexerType"))->Fill(aod::hf_cand::VertexerType::KfParticle);
    }

    // initailize HF event selection helper
    // ------------------------------------
    hfEvSel.init(registry, &zorroSummary);

    // Histograms for QA
    // -----------------
    registry.add("ReconstructedDecayChannel", "DecyayChannel", {kTH1F, {{3, 0.0, 3.0}}});
    registry.get<TH1>(HIST("ReconstructedDecayChannel"))->GetXaxis()->SetBinLabel(1 + hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi, "To #Xi #pi");
    registry.get<TH1>(HIST("ReconstructedDecayChannel"))->GetXaxis()->SetBinLabel(1 + hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi, "To #Omega #pi");
    registry.get<TH1>(HIST("ReconstructedDecayChannel"))->GetXaxis()->SetBinLabel(1 + hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK, "To #Omega K");

    if (xipiEnabledDca != 0 || xipiEnabledKf != 0) {
      hInvMassCharmBaryonToXiPi = registry.add<TH1>("hInvMassCharmBaryonToXiPi", "Charm baryon invariant mass - #Xi #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{configs.nBinMassCharmBaryon, configs.minMassCharmBaryon, configs.maxMassCharmBaryon}}});
      hCandidateCounterToXiPi = registry.add<TH1>("hCandidateCounterToXiPi", "Candidate counter wrt derived data - #Xi #pi decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});
      registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + All, "Total");
      registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + HfFlagPass, "HfFlagPass");
      registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + CascReconstructed, "CascReconstructed");
      registry.get<TH1>(HIST("hCandidateCounterToXiPi"))->GetXaxis()->SetBinLabel(1 + VertexFit, "VertexFit");
      registry.get<TH1>(HIST("ReconstructedDecayChannel"))->Fill(hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi);
    }

    if (omegapiEnabledDca != 0 || omegapiEnabledKf != 0) {
      hInvMassCharmBaryonToOmegaPi = registry.add<TH1>("hInvMassCharmBaryonToOmegaPi", "Charm baryon invariant mass - #Omega #pi decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{configs.nBinMassCharmBaryon, configs.minMassCharmBaryon, configs.maxMassCharmBaryon}}});
      hCandidateCounterToOmegaPi = registry.add<TH1>("hCandidateCounterToOmegaPi", "Candidate counter wrt derived data - #Omega #pi decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});
      registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + All, "Total");
      registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + HfFlagPass, "HfFlagPass");
      registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + CascReconstructed, "CascReconstructed");
      registry.get<TH1>(HIST("hCandidateCounterToOmegaPi"))->GetXaxis()->SetBinLabel(1 + VertexFit, "VertexFit");
      registry.get<TH1>(HIST("ReconstructedDecayChannel"))->Fill(hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi);
    }

    if (omegakaEnabledDca != 0 || omegakaEnabledKf != 0) {
      hInvMassCharmBaryonToOmegaKa = registry.add<TH1>("hInvMassCharmBaryonToOmegaKa", "Charm baryon invariant mass - #Omega K decay;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{configs.nBinMassCharmBaryon, configs.minMassCharmBaryon, configs.maxMassCharmBaryon}}});
      hCandidateCounterToOmegaKa = registry.add<TH1>("hCandidateCounterToOmegaKa", "Candidate counter wrt derived data - #Omega K decay;status;entries", {HistType::kTH1D, {{4, -0.5, 3.5}}});
      registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + All, "Total");
      registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + HfFlagPass, "HfFlagPass");
      registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + CascReconstructed, "CascReconstructed");
      registry.get<TH1>(HIST("hCandidateCounterToOmegaKa"))->GetXaxis()->SetBinLabel(1 + VertexFit, "VertexFit");
      registry.get<TH1>(HIST("ReconstructedDecayChannel"))->Fill(hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK);
    }

    registry.add("hCascMass", "Inv mass of reconstructed cascade;Inv mass;Entries", {HistType::kTH1F, {{configs.nBinMassCasc, configs.minMassCasc, configs.maxMassCasc}}});
    registry.add("hCascPt", "Pt of reconstructed cascade;pT;Entries", {HistType::kTH1F, {{configs.nBinPtCasc, configs.minPtCasc, configs.maxPtCasc}}});

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

  // template function for running charm baryon reconstruction with DCAFitter
  /// \brief centEstimator is for different centrality estimators
  /// \brief decayChannel is for different decay channels. 0 for XiczeroOmegaczeroToXiPi, 1 for OmegaczeroToOmegaPi, 2 for OmegaczeroToOmeagaK
  /// \brief Colls is for collision tables joined with different centraltiy estimators
  /// \brief Hist is for QA histograms
  template <o2::hf_centrality::CentralityEstimator centEstimator, int decayChannel, typename Colls, typename Hist>
  void runCreatorWithDCAFitter(Colls const&,
                               aod::HfCascLf2Prongs const& candidates,
                               aod::Cascades const&, // -> Internal cascade building
                               aod::V0s const&,      // -> Internal v0 building
                               aod::TracksWCovDca const& tracks,
                               TracksWCovIU const& lfTracks,
                               aod::BCsWithTimestamps const&,
                               Hist& hInvMassCharmBaryon,
                               Hist& hCandCounter)
  {
    // arrays which holds values for different decay channel reconstruction

    // Loop over candidate
    for (auto const& cand : candidates) {

      // Fill cascandidates before selection
      if (configs.fillHistograms) {
        hCandCounter->Fill(All);
      }

      // Apply hfflag selection for different candidate reconstruction
      if (!TESTBIT(cand.hfflag(), decayChannel)) {
        continue;
      } else {
        if (configs.fillHistograms) {
          hCandCounter->Fill(HfFlagPass);
        }
      }

      // Event selection
      auto collision = cand.collision_as<Colls>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) { // None of the event selection satisfied -> Reject this candidate
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
      straHelper.fitter.setBz(magneticField); // -> Magnetic field setting for internal cascade building
      df.setBz(magneticField);                // -> Magnetic field setting for charm baryon building

      //------------------Intenal Cascade building------------------
      auto cascAodElement = cand.cascade_as<aod::Cascades>();
      auto v0AodElement = cascAodElement.v0_as<aod::V0s>();
      auto posTrack = lfTracks.rawIteratorAt(v0AodElement.posTrackId());
      auto negTrack = lfTracks.rawIteratorAt(v0AodElement.negTrackId());
      auto bachTrack = lfTracks.rawIteratorAt(cascAodElement.bachelorId());

      // Make cascade starting from V0
      // If success, fill Cascade and V0 information for reconstruction
      if (!straHelper.buildCascadeCandidate(collision.globalIndex(),
                                            collision.posX(), collision.posY(), collision.posZ(),
                                            posTrack,
                                            negTrack,
                                            bachTrack,
                                            false, // calculateBachelorBaryonVariable
                                            false, // useCascadeMomentumAtPV
                                            true)) {
        continue;
      } else {
        float storeMass = (decayChannel == 0) ? straHelper.cascade.massXi : straHelper.cascade.massOmega;
        float storePt = RecoDecay::pt(straHelper.cascade.cascadeMomentum);
        registry.fill(HIST("hCascMass"), storeMass);
        registry.fill(HIST("hCascPt"), storePt);
        if (configs.fillHistograms) {
          hCandCounter->Fill(CascReconstructed);
        }
      }

      //------------------------------Info of V0 and cascade------------------------------
      // Cascade quantities from LF strangeness builder
      // TODO: check if this is necessary or the covariance matrix can be directly used from straHelper.cascade.covariance
      std::array<float, 21> covCasc = {0.};
      constexpr int NumCovElement = 6;
      constexpr int MomInd[NumCovElement] = {9, 13, 14, 18, 19, 20};
      for (int i = 0; i < NumCovElement; i++) {
        covCasc[MomInd[i]] = straHelper.cascade.covariance[MomInd[i]];
        covCasc[i] = straHelper.cascade.covariance[i];
      }

      //------------------------------Create cascade track------------------------------
      o2::track::TrackParCov trackCasc;
      if (bachTrack.sign() < 0) { // Xi- or Omega-
        trackCasc = o2::track::TrackParCov(straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum, covCasc, -1, true);
      } else if (bachTrack.sign() > 0) { // Xi+ or Omega+
        trackCasc = o2::track::TrackParCov(straHelper.cascade.cascadePosition, straHelper.cascade.cascadeMomentum, covCasc, 1, true);
      } else {
        continue;
      }

      trackCasc.setAbsCharge(1);
      trackCasc.setPID(trackPidOfCascade);

      //------------------------------Fit SV & Create Charm Baryon track------------------------------
      auto trackCharmBachelor = tracks.rawIteratorAt(cand.prong0Id());
      auto trackParCovCharmBachelor = getTrackParCov(trackCharmBachelor);
      try {
        if (df.process(trackCasc, trackParCovCharmBachelor) == 0) {
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitter cannot work. Skipping this candidate";
        continue;
      }

      if (configs.fillHistograms) {
        hCandCounter->Fill(VertexFit);
      }

      //------------------------------Calculate physical properties-----------------------------
      // get track momenta
      std::array<float, 3> pVecCasc, pVecCharmBach;
      df.getTrack(0).getPxPyPzGlo(pVecCasc);
      df.getTrack(1).getPxPyPzGlo(pVecCharmBach);

      // get PV Properties
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

      // get SV Properties
      auto const& secondaryVertex = df.getPCACandidate();
      auto covMatrixSV = df.calcPCACovMatrixFlat();
      std::array<float, 3> svCoord = {secondaryVertex[0], secondaryVertex[1], secondaryVertex[2]};

      float phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

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

      // get Casc and Lambda invariant mass - from LF Table
      float mCasc = (decayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) ? straHelper.cascade.massOmega : straHelper.cascade.massXi;
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
                       trackCharmBachelor.globalIndex(), cascAodElement.bachelorId(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                       cand.hfflag());
    } // candidate loop
  } // end of run function

  // template function for running Charm Baryon reconstruction via KFParticle method
  /// \brief centEstimator is for different centrality estimators
  /// \brief decayChannel is for different decay channels. 0 for XiczeroOmegaczeroToXiPi, 1 for OmegaczeroToOmegaPi, 2 for OmegaczeroToOmegaK
  /// \brief Colls is for collision tables joined with different centrality estimators
  /// \brief Hist is for QA histograms
  template <o2::hf_centrality::CentralityEstimator centEstimator, int decayChannel, typename Colls, typename Hist>
  void runCreatorWithKfParticle(Colls const&,
                                aod::HfCascLf2Prongs const& candidates,
                                aod::Cascades const&, // -> Implemented for internal cascade building
                                aod::V0s const&,      // -> Implemented for internal cascade building
                                TracksWCovDcaExtraPidPrPiKa const& tracks,
                                TracksWCovExtraPidIU const& lfTracks,
                                aod::BCsWithTimestamps const&,
                                Hist& hInvMassCharmBaryon,
                                Hist& hCandCounter)
  {
    // Loop over candidates
    for (auto const& cand : candidates) {

      // Fill cascandidates before selection
      if (configs.fillHistograms) {
        hCandCounter->Fill(All);
      }

      // Apply hfflag selection
      if (!TESTBIT(cand.hfflag(), decayChannel)) {
        continue;
      } else {
        if (configs.fillHistograms) {
          hCandCounter->Fill(HfFlagPass);
        }
      }

      // Event selection
      auto collision = cand.collision_as<Colls>();
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, centEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) { // None of the event seletion satisfied -> Reject this candidate
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
      // magnetic field setting for KFParticle
      straHelper.fitter.setBz(magneticField); // -> Manetic field setting for internal cascade building
      KFParticle::SetField(magneticField);    // -> Magnetic field setting for CharmBaryon building

      //------------------Intenal Cascade building------------------
      auto cascAodElement = cand.cascade_as<aod::Cascades>();
      auto v0AodElement = cascAodElement.v0_as<aod::V0s>();
      auto posTrack = lfTracks.rawIteratorAt(v0AodElement.posTrackId());
      auto negTrack = lfTracks.rawIteratorAt(v0AodElement.negTrackId());
      auto bachTrack = lfTracks.rawIteratorAt(cascAodElement.bachelorId());

      // Make cascade starting from V0
      // If success, fill Cascade and V0 information for reconstruction
      if (!straHelper.buildCascadeCandidateWithKF(collision.globalIndex(),
                                                  collision.posX(), collision.posY(), collision.posZ(),
                                                  posTrack,
                                                  negTrack,
                                                  bachTrack,
                                                  false, // calculateBachelorBaryonVariables
                                                  LFConfigs.kfConstructMethodFromLF,
                                                  LFConfigs.kfTuneForOmegaFromLF,
                                                  LFConfigs.kfUseV0MassConstraintFromLF,
                                                  LFConfigs.kfUseCascadeMassConstraintFromLF,
                                                  LFConfigs.kfDoDCAFitterPreMinimV0FromLF,
                                                  LFConfigs.kfDoDCAFitterPreMinimCascFromLF)) {
        continue;
      } else {
        float storeMass = (decayChannel == 0) ? straHelper.cascade.massXi : straHelper.cascade.massOmega;
        float storePt = RecoDecay::pt(straHelper.cascade.cascadeMomentum);
        registry.fill(HIST("hCascMass"), storeMass);
        registry.fill(HIST("hCascPt"), storePt);
        if (configs.fillHistograms) {
          hCandCounter->Fill(CascReconstructed);
        }
      }

      //------------------------------Cascade pre-selection------------------------------
      // ! only for Xic0, Omegac0 -> Omega K
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        if (configs.doCascadePreselection) {
          if (std::abs(straHelper.cascade.cascadeDCAxy) > configs.dcaXYToPVCascadeMax) {
            continue;
          }
          // pre selection on dcaV0Daughters
          if (std::abs(straHelper.cascade.v0DaughterDCA) > configs.dcaV0DaughtersMax) {
            continue;
          }
          // pre selection on dcaCascDaughters
          if (std::abs(straHelper.cascade.cascadeDaughterDCA) > configs.dcaCascDaughtersMax) {
            continue;
          }
          // pre selection on invariantmass
          if (std::abs(straHelper.cascade.massOmega - massOfCascade) > configs.massToleranceCascade) {
            continue;
          }
        }
      }

      //----------Create charm bayron as KF Particle object starting from V0----------
      // Create KFParticle object of V0 Daughter & Bachelor
      const KFPTrack kfTrack0 = createKFPTrackFromTrack(posTrack);
      const KFPTrack kfTrack1 = createKFPTrackFromTrack(negTrack);
      const KFPTrack kfTrackBach = createKFPTrackFromTrack(bachTrack);

      bool isAnti = (bachTrack.signed1Pt() > 0 ? true : false);
      int bachCharge = (bachTrack.signed1Pt() > 0 ? +1 : -1);

      KFParticle kfPos(kfTrack0, (isAnti ? -pdgIdOfV0DauNeg : pdgIdOfV0DauPos));
      KFParticle kfNeg(kfTrack1, (isAnti ? -pdgIdOfV0DauPos : pdgIdOfV0DauNeg));
      KFParticle kfBach(kfTrackBach, (isAnti ? -pdgIdOfBach : pdgIdOfBach));
      KFParticle kfBachRej(kfTrackBach, (isAnti ? -pdgIdOfBachRej : pdgIdOfBachRej)); // Rej -> Used for Omegac0->OmegaPi only

      // ~~~~~~~Construct V0 with KF~~~~~~~
      const KFParticle* v0Daughters[2] = {&kfPos, &kfNeg};
      KFParticle kfV0;
      kfV0.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfV0.Construct(v0Daughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct v0" << e.what();
        continue;
      }

      // Require lambda pre-selection before mass constraint
      float massLam, sigMassLam;
      kfV0.GetMass(massLam, sigMassLam);
      if (std::abs(massLam - MassLambda0) > configs.lambdaMassWindow) {
        continue;
      }
      if (sigMassLam <= 0) {
        continue;
      }
      if ((kfV0.GetNDF() <= 0) || (kfV0.GetChi2() <= 0)) {
        continue;
      }

      // Set mass constraint to lambda
      if (configs.kfUseV0MassConstraint) {
        kfV0.SetNonlinearMassConstraint(massOfV0);
      }
      kfV0.TransportToDecayVertex();
      kfV0.SetPDG(pdgIdOfV0);

      //~~~~~~~Construct cascade with KF~~~~~~~
      const KFParticle* cascDaughters[2] = {&kfBach, &kfV0};
      const KFParticle* cascDaughtersRej[2] = {&kfBachRej, &kfV0};
      KFParticle kfCasc, kfCascRej;

      kfCasc.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfCasc.Construct(cascDaughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Cascade: " << e.what();
        continue;
      }

      float massCasc, sigMassCasc, massCascRej, sigMassCascRej;
      kfCasc.GetMass(massCasc, sigMassCasc);

      if (sigMassCasc <= 0) {
        continue;
      }
      if (std::abs(massCasc - massOfCascade) > configs.massToleranceCascade) {
        continue;
      }
      if (kfCasc.GetNDF() <= 0 || kfCasc.GetChi2() <= 0) {
        continue;
      }
      kfCasc.SetPDG(pdgIdOfCascade);

      // perform cascade building on casc_rej - only for Omega
      if constexpr (decayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        kfCascRej.SetConstructMethod(configs.kfConstructMethod);
        try {
          kfCascRej.Construct(cascDaughtersRej, 2);
        } catch (std::runtime_error& e) {
          LOG(debug) << "Failed to construct Cascade_rej: " << e.what();
          continue;
        }

        kfCascRej.GetMass(massCascRej, sigMassCascRej);
        kfCascRej.SetPDG(pdgIdOfCascade);
      }

      // Set mass constraint to cascade
      if (configs.kfUseCascadeMassConstraint) {
        kfCasc.SetNonlinearMassConstraint(massOfCascade);
      }
      kfCasc.TransportToDecayVertex();

      //~~~~~~~Construct Charm Baryon with KF~~~~~~~
      auto trackCharmBachelor = tracks.rawIteratorAt(cand.prong0Id());
      const KFPTrack kfTrackCharmBach = createKFPTrackFromTrack(trackCharmBachelor);
      const KFParticle kfCharmBach(kfTrackCharmBach, (isAnti ? -pdgIdOfCharmBach : pdgIdOfCharmBach));
      const KFParticle* charmBaryonDaughters[2] = {&kfCharmBach, &kfCasc};

      KFParticle kfCharmBaryon;
      kfCharmBaryon.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfCharmBaryon.Construct(charmBaryonDaughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Charm baryon: " << e.what();
        continue;
      }

      float massCharmBaryon, sigMassCharmBaryon;
      kfCharmBaryon.GetMass(massCharmBaryon, sigMassCharmBaryon);
      if (sigMassCharmBaryon <= 0) {
        continue;
      }
      if (kfCharmBaryon.GetNDF() <= 0 || kfCharmBaryon.GetChi2() <= 0) {
        continue;
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

      // TODO: check if the impact parameter can be calculated from KFParticle directly
      auto primaryVertex = getPrimaryVertex(collision);
      o2::dataformats::DCA impactParameterV0DauPos, impactParameterV0DauNeg, impactParameterCascBach;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0DauPos, 2.f, matCorr, &impactParameterV0DauPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0DauNeg, 2.f, matCorr, &impactParameterV0DauNeg);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovBach, 2.f, matCorr, &impactParameterCascBach);

      o2::dataformats::DCA impactParameterCasc, impactParameterCharmBach;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCharmBach, 2.f, matCorr, &impactParameterCharmBach);

      float kfDcaV0Daughters = kfNeg.GetDistanceFromParticle(kfPos);
      float kfDcaCascDaughters = kfBachToCasc.GetDistanceFromParticle(kfV0ToCasc);
      float kfDcaCharmBaryonDaughters = kfCharmBachToCharmBaryon.GetDistanceFromParticle(kfCascToCharmBaryon);

      float dlV0, errDlV0;
      kfV0ToPv.GetDecayLength(dlV0, errDlV0);

      float dlCasc, errDlCasc;
      kfCascToPv.GetDecayLength(dlCasc, errDlCasc);

      float dlCharmBaryon, dlXYCharmBaryon, errDlCharmBaryon, errDlXYCharmBaryon;
      kfCharmBaryonToPv.GetDecayLength(dlCharmBaryon, errDlCharmBaryon);
      kfCharmBaryonToPv.GetDecayLengthXY(dlXYCharmBaryon, errDlXYCharmBaryon);

      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> svCoord = {kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()};
      std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      std::array<float, 3> vertexCasc = {kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()};
      std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()};
      std::array<float, 3> pVecCasc = {kfCasc.GetPx(), kfCasc.GetPy(), kfCasc.GetPz()};

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
                       errDlCharmBaryon, errDlXYCharmBaryon,
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
                       trackCharmBachelor.globalIndex(), cascAodElement.bachelorId(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                       cand.hfflag());

      rowCandCascBatchKf(kfCharmBaryon.Chi2(), kfCasc.Chi2(), kfV0.Chi2(),
                         kfCharmBaryon.GetNDF(), kfCasc.GetNDF(), kfV0.GetNDF(),
                         dlCharmBaryon / errDlCharmBaryon, dlCasc / errDlCasc, dlV0 / errDlV0,
                         kfCharmBaryonToPv.GetChi2(), kfCascToPv.GetChi2(), kfV0ToPv.GetChi2(),
                         kfCharmBaryonToPv.GetNDF(), kfCascToPv.GetNDF(), kfV0ToPv.GetNDF(),
                         kfCascToCharmBaryon.GetChi2(), kfV0ToCasc.GetChi2(),
                         kfCascToCharmBaryon.GetNDF(), kfV0ToCasc.GetNDF());
    } // end candidate loop
  } // end of runCreator

  /////////////////////////////////////////////////////
  ///                                               ///
  ///        Process functions with DCAFitter       ///
  ///                                               ///
  /////////////////////////////////////////////////////

  /*~~~~~~~~~~~~~~*/
  /*~~~To Xi Pi~~~*/
  /*~~~~~~~~~~~~~~*/
  void processToXiPiWithDCAFitterNoCent(SelectedCollisions const& collisions,
                                        aod::HfCascLf2Prongs const& candidates,
                                        aod::Cascades const& cascades,
                                        aod::V0s const& v0s,
                                        aod::TracksWCovDca const& tracks,
                                        TracksWCovIU const& lfTracks,
                                        aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithDCAFitterNoCent, "Charm candidte reconstruction with Xi Pi via DcaFitter method, no centrality", true);

#if 0
  void processToXiPiWithDCAFitterNoCentWithTrackedCasc(SelectedCollisions const& collisions,
                                                       aod::HfCascLf2Prongs const& candidates,
                                                       aod::TrackedCascades const& cascades,
                                                       aod::V0s const& v0s,
                                                       TrackedCascLinked const& trackedCascLinked,
                                                       TracksWCovDcaExtraPidPrPiKa const& tracks,
                                                       aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, trackedCascFull, trackedCascLinked, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithDCAFitterNoCentWithTrackedCasc, "Charm candidte reconstruction with Xi Pi via DcaFitter method with tracked cascade, no centrality", false);
#endif

  void processToXiPiWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                          aod::HfCascLf2Prongs const& candidates,
                                          aod::Cascades const& cascades,
                                          aod::V0s const& v0s,
                                          aod::TracksWCovDca const& tracks,
                                          TracksWCovIU const& lfTracks,
                                          aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithDCAFitterCentFT0C, "Charm candidate reconstruction with Xi Pi via DcaFitter method, centrality selection on FT0C", false);

  void processToXiPiWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                          aod::HfCascLf2Prongs const& candidates,
                                          aod::Cascades const& cascades,
                                          aod::V0s const& v0s,
                                          aod::TracksWCovDca const& tracks,
                                          TracksWCovIU const& lfTracks,
                                          aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithDCAFitterCentFT0M, "Charm candidate reconstruction with Xi Pi via DcaFitter method, centrality selection on FT0M", false);

  /*~~~~~~~~~~~~~~~~~*/
  /*~~~To Omega Pi~~~*/
  /*~~~~~~~~~~~~~~~~~*/
  void processToOmegaPiWithDCAFitterNoCent(SelectedCollisions const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           aod::TracksWCovDca const& tracks,
                                           TracksWCovIU const& lfTracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithDCAFitterNoCent, "Charm candidte reconstruction with Omega Pi via DcaFitter method, no centrality", false);

  void processToOmegaPiWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             aod::TracksWCovDca const& tracks,
                                             TracksWCovIU const& lfTracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithDCAFitterCentFT0C, "Charm candidate reconstruction with Omega Pi via DcaFitter method, centrality selection on FT0C", false);

  void processToOmegaPiWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             aod::TracksWCovDca const& tracks,
                                             TracksWCovIU const& lfTracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithDCAFitterCentFT0M, "Charm candidate reconstruction with Omega Pi via DcaFitter method, centrality selection on FT0M", false);

  /*~~~~~~~~~~~~~~~~~*/
  /*~~~To Omega Ka~~~*/
  /*~~~~~~~~~~~~~~~~~*/
  void processToOmegaKaWithDCAFitterNoCent(SelectedCollisions const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           aod::TracksWCovDca const& tracks,
                                           TracksWCovIU const& lfTracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithDCAFitterNoCent, "Charm candidte reconstruction with Omega Ka via DcaFitter method, no centrality", false);

  void processToOmegaKaWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             aod::TracksWCovDca const& tracks,
                                             TracksWCovIU const& lfTracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithDCAFitterCentFT0C, "Charm candidate reconstruction with Omega Ka via DcaFitter method, centrality selection on FT0C", false);

  void processToOmegaKaWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             aod::TracksWCovDca const& tracks,
                                             TracksWCovIU const& lfTracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithDCAFitterCentFT0M, "Charm candidate reconstruction with Omega Ka via DcaFitter method, centrality selection on FT0M", false);

  /////////////////////////////////////////////////////
  ///                                               ///
  ///        Process functions with KFParticle      ///
  ///                                               ///
  /////////////////////////////////////////////////////

  /*~~~~~~~~~~~~~~*/
  /*~~~To Xi Pi~~~*/
  /*~~~~~~~~~~~~~~*/
  void processToXiPiWithKFParticleNoCent(SelectedCollisions const& collisions,
                                         aod::HfCascLf2Prongs const& candidates,
                                         aod::Cascades const& cascades,
                                         aod::V0s const& v0s,
                                         TracksWCovDcaExtraPidPrPiKa const& tracks,
                                         TracksWCovExtraPidIU const& lfTracks,
                                         aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithKFParticleNoCent, "Charm Baryon decaying to Xi Pi reconstruction via KFParticle method, no centrality", false);

  void processToXiPiWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           TracksWCovDcaExtraPidPrPiKa const& tracks,
                                           TracksWCovExtraPidIU const& lfTracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithKFParticleCentFT0C, "Charm Baryon decaying to Xi Pi reconstruction via KFParticle method, centrality on FT0C", false);

  void processToXiPiWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           TracksWCovDcaExtraPidPrPiKa const& tracks,
                                           TracksWCovExtraPidIU const& lfTracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithKFParticleCentFT0M, "Charm Baryon decaying to Xi Pireconstruction via KFParticle method, centrality on FT0M", false);

  /*~~~~~~~~~~~~~~~~~*/
  /*~~~To Omega Pi~~~*/
  /*~~~~~~~~~~~~~~~~~*/
  void processToOmegaPiWithKFParticleNoCent(SelectedCollisions const& collisions,
                                            aod::HfCascLf2Prongs const& candidates,
                                            aod::Cascades const& cascades,
                                            aod::V0s const& v0s,
                                            TracksWCovDcaExtraPidPrPiKa const& tracks,
                                            TracksWCovExtraPidIU const& lfTracks,
                                            aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithKFParticleNoCent, "Charm Baryon decaying to Omega Pi reconstruction via KFParticle method, no centrality", false);

  void processToOmegaPiWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              TracksWCovExtraPidIU const& lfTracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithKFParticleCentFT0C, "Charm Baryon decaying to Omega Pi reconstruction via KFParticle method, centrality on FT0C", false);

  void processToOmegaPiWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              TracksWCovExtraPidIU const& lfTracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithKFParticleCentFT0M, "Charm Baryong decaying to Omega Pi reconstruction via KFParticle method, centrality on FT0M", false);

  /*~~~~~~~~~~~~~~~~~*/
  /*~~~To Omega Ka~~~*/
  /*~~~~~~~~~~~~~~~~~*/
  void processToOmegaKaWithKFParticleNoCent(SelectedCollisions const& collisions,
                                            aod::HfCascLf2Prongs const& candidates,
                                            aod::Cascades const& cascades,
                                            aod::V0s const& v0s,
                                            TracksWCovDcaExtraPidPrPiKa const& tracks,
                                            TracksWCovExtraPidIU const& lfTracks,
                                            aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithKFParticleNoCent, "Charm Baryon decaying to Omega Ka reconstruction via KFParticle method, no centrality", false);

  void processToOmegaKaWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Cs> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              TracksWCovExtraPidIU const& lfTracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithKFParticleCentFT0C, "Charm Baryon decaying to Omega Ka reconstruction via KFParticle method, centrality on FT0C", false);

  void processToOmegaKaWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              TracksWCovExtraPidIU const& lfTracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithKFParticleCentFT0M, "Charm Baryong decaying to Omega Ka reconstruction via KFParticle method, centrality on FT0M", false);

  ///////////////////////////////////////////////////////////////
  ///                                                         ///
  ///        Process functions for Collision monitoring       ///
  ///                                                         ///
  ///////////////////////////////////////////////////////////////

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
  Produces<aod::HfCandCascBachMCRecQa> rowMcMatchRecCascBach;
  Produces<aod::HfCandCascBachMCGenQa> rowMcMatchGenCascBach;

  // Configurables
  struct : ConfigurableGroup {
    Configurable<bool> rejectBackground{"rejectBackground", false, "Reject particles from background events"}; // -> Used for only Xic0
    Configurable<bool> acceptTrackInteractionWithMaterial{"acceptTrackInteractionWithMaterial", false, "Accept candidates with final daughters interacting with materials"};
    // Configurable<bool> matchDecayedPions{"matchedDecayedPions", true, "Match also candidates with daughter pion tracks that decay with kinked toploogy"};
  } configs;

  std::array<int, 4> pdgOfCharmBaryon{+kXiC0, +kOmegaC0, +kOmegaC0, +kOmegaC0};
  std::array<int, 4> pdgOfCascade{+kXiMinus, +kXiMinus, +kOmegaMinus, +kOmegaMinus};
  std::array<int, 4> pdgOfCharmBachelor{+kPiPlus, +kPiPlus, +kPiPlus, +kKPlus};
  std::array<int, 4> pdgOfBachelor{+kPiMinus, +kPiMinus, +kKMinus, +kKMinus};

  std::array<int8_t, 4> newDecayChannel{hf_cand_casc_bach_qa::DecayChannelMain::Xic0ToXiPi,
                                        hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToXiPi,
                                        hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToOmegaPi,
                                        hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToOmegaKa};

  // Table aliases
  using TracksWMcIU = soa::Join<aod::TracksIU, McTrackLabels>;
  using McCollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels>;
  using McCollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels, aod::CentFT0Cs>;
  using McCollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>; // -> Used for subscription for process functions of centrality with FT0Ms
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
      if (device.name.compare("hf-candidate-creator-xic0-omegac0-qa") == 0) {
        hfEvSelMc.init(device, registry);
        break;
      }
    }
  }

  template <o2::hf_centrality::CentralityEstimator centEstimator, int decayChannel, typename TRecoCand, typename Colls, typename McCollisions>
  void runXic0Omegac0Mc(TRecoCand const& candidates,
                        TracksWMcIU const&,
                        aod::McParticles const& mcParticles,
                        Colls const& collsWithMcLabels,
                        McCollisions const& mcCollisions,
                        BCsInfo const&)
  {
    // Match reconstructed candidates //
    for (const auto& candidate : candidates) {
      int indexRecCharmBaryon = -1;
      int indexRecCascade = -1;
      int indexRecV0 = -1;
      int8_t sign = 0;
      int8_t signCasc = 0;
      int8_t signV0 = 0;
      int8_t flagChannelMain = 0;
      int8_t flagChannelResonant = 0; // TODO: Implement resonant channel matching
      int8_t origin = 0;
      int8_t nKinkedTracks = 0;             // TODO: Implement kinked track matching
      int8_t nInteractionsWithMaterial = 0; // TODO: Implement interaction with material matching

      std::vector<int> idxBhadMothers{};

      auto arrayDaughters = std::array{candidate.template bachelorFromCharmBaryon_as<TracksWMcIU>(),
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

      // CharmBaryon -> Charm bachelor + Cascade
      indexRecCharmBaryon = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgOfCharmBaryon[decayChannel], std::array{pdgOfCharmBachelor[decayChannel], pdgOfBachelor[decayChannel], +kProton, +kPiMinus}, true, &sign, 3);

      if (indexRecCharmBaryon > -1) {
        // Cascade -> Bachelor + V0
        indexRecCascade = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgOfCascade[decayChannel], std::array{pdgOfBachelor[decayChannel], +kProton, +kPiMinus}, true, &signCasc, 2);

        if (indexRecCascade > -1) {
          // V0 -> Pos + Neg
          indexRecV0 = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);

          if (indexRecV0 > -1) {
            flagChannelMain = sign * newDecayChannel[decayChannel];
          }
        }
      }

      // Check if charm baryon is from b-hadron decay(prompt vs non-prompt)
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

      if (rejectionMask != 0) { // At least one event selection not satisfied --> Reject all particles from this event
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
                  flagChannelMain = sign * newDecayChannel[decayChannel];
                  ;
                }
              } // V0 daughter loop
            } // cascade daughter loop
          } // charm daughter loop
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
  } // template run function

  void processMcEmpty(aod::Collisions const&)
  {
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcEmpty, "Empty process function to prevent workflow from getting stuck", true);

  /////////////////////////////////////////////////////
  ///                                               ///
  ///        Process functions with DCAFitter       ///
  ///                                               ///
  /////////////////////////////////////////////////////

  //~~~~~~~~~~~~~~~~//
  //~~~~To Xi Pi~~~~//
  //~~~~~~~~~~~~~~~~//
  void processMcXicToXiPiWithDCAFitterNoCent(aod::HfCandToXiPi const& candidates,
                                             TracksWMcIU const& tracks,
                                             aod::McParticles const& mcParticles,
                                             soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                             McCollisionsNoCents const& collsWithMcLabels,
                                             BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. No cents", false);

  void processMcXicToXiPiWithDCAFitterCentFT0C(aod::HfCandToXiPi const& candidates,
                                               TracksWMcIU const& tracks,
                                               aod::McParticles const& mcParticles,
                                               soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                               McCollisionsFT0Cs const& collsWithMcLabels,
                                               BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. Cents with FT0C", false);

  void processMcXicToXiPiWithDCAFitterCentFT0M(aod::HfCandToXiPi const& candidates,
                                               TracksWMcIU const& tracks,
                                               aod::McParticles const& mcParticles,
                                               McCollisionsCentFT0Ms const& mcCollisions,
                                               McCollisionsFT0Ms const& collsWithMcLabels,
                                               BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterCentFT0M, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. Cents with FT0M", false);

  void processMcOmegacToXiPiWithDCAFitterNoCent(aod::HfCandToXiPi const& candidates,
                                                TracksWMcIU const& tracks,
                                                aod::McParticles const& mcParticles,
                                                soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                                McCollisionsNoCents const& collsWithMcLabels,
                                                BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToXiPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Xi Pi. No cents", false);

  void processMcOmegacToXiPiWithDCAFitterCentFT0C(aod::HfCandToXiPi const& candidates,
                                                  TracksWMcIU const& tracks,
                                                  aod::McParticles const& mcParticles,
                                                  soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                                  McCollisionsFT0Cs const& collsWithMcLabels,
                                                  BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToXiPiWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Omeagc0 to Xi Pi. Cents with FT0C", false);

  void processMcOmegacToXiPiWithDCAFitterCentFT0M(aod::HfCandToXiPi const& candidates,
                                                  TracksWMcIU const& tracks,
                                                  aod::McParticles const& mcParticles,
                                                  McCollisionsCentFT0Ms const& mcCollisions,
                                                  McCollisionsFT0Ms const& collsWithMcLabels,
                                                  BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToXiPiWithDCAFitterCentFT0M, "Perform MC matching of DCAFitter reconstructed Omegac0 to Xi Pi. Cents with FT0M", false);

  //~~~~~~~~~~~~~~~~~~~//
  //~~~~To Omega Pi~~~~//
  //~~~~~~~~~~~~~~~~~~~//
  void processMcOmegacToOmegaPiWithDCAFitterNoCent(aod::HfCandToOmegaPi const& candidates,
                                                   TracksWMcIU const& tracks,
                                                   aod::McParticles const& mcParticles,
                                                   soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                                   McCollisionsNoCents const& collsWithMcLabels,
                                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Pi. No cents", false);

  void processMcOmegacToOmegaPiWithDCAFitterCentFT0C(aod::HfCandToOmegaPi const& candidates,
                                                     TracksWMcIU const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                                     McCollisionsFT0Cs const& collsWithMcLabels,
                                                     BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaPiWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Pi. Cents with FT0C", false);

  void processMcOmegacToOmegaPiWithDCAFitterCentFT0M(aod::HfCandToOmegaPi const& candidates,
                                                     TracksWMcIU const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     McCollisionsCentFT0Ms const& mcCollisions,
                                                     McCollisionsFT0Ms const& collsWithMcLabels,
                                                     BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaPiWithDCAFitterCentFT0M, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Pi. Cents with FT0M", false);

  //~~~~~~~~~~~~~~~~~~//
  //~~~~To Omega Ka~~~~//
  //~~~~~~~~~~~~~~~~~~//
  void processMcOmegacToOmegaKaWithDCAFitterNoCent(aod::HfCandToOmegaK const& candidates,
                                                   TracksWMcIU const& tracks,
                                                   aod::McParticles const& mcParticles,
                                                   soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                                   McCollisionsNoCents const& collsWithMcLabels,
                                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaKaWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Ka. No cents", false);

  void processMcOmegacToOmegaKaWithDCAFitterCentFT0C(aod::HfCandToOmegaK const& candidates,
                                                     TracksWMcIU const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                                     McCollisionsFT0Cs const& collsWithMcLabels,
                                                     BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaKaWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Ka. Cents with FT0C", false);

  void processMcOmegacToOmegaKaWithDCAFitterCentFT0M(aod::HfCandToOmegaK const& candidates,
                                                     TracksWMcIU const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     McCollisionsCentFT0Ms const& mcCollisions,
                                                     McCollisionsFT0Ms const& collsWithMcLabels,
                                                     BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaKaWithDCAFitterCentFT0M, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Ka. Cents with FT0M", false);

  /////////////////////////////////////////////////////
  ///                                               ///
  ///        Process functions with KFParticle      ///
  ///                                               ///
  /////////////////////////////////////////////////////
  void processMcXicToXiPiWithKFParticleNoCent(aod::HfCandToXiPiKf const& candidates,
                                              TracksWMcIU const& tracks,
                                              aod::McParticles const& mcParticles,
                                              soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions,
                                              McCollisionsNoCents const& collsWithMcLabels,
                                              BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithKFParticleNoCent, "Perform MC matching of KFParticle reconstructed Xic0 to Xi Pi. No cents", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0Qa>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0QaMc>(cfgc)};
}
