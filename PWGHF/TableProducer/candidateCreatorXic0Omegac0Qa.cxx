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

#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h" // -> Added to test removal of strangeness builder workflow

#include "Common/Core/ZorroSummary.h"
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
#include <memory>
#include <string>
#include <utility>
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

  // Cursor to fill tables
  struct : ProducesGroup {
    // Candidates created with DCAFitter
    Produces<aod::HfCandToXiPi> rowCandToXiPi;
    Produces<aod::HfCandToOmegaPi> rowCandToOmegaPi;
    Produces<aod::HfCandToOmegaK> rowCandToOmegaKa;
    // Candidate created with KFParticle
    Produces<aod::HfCandToXiPiKf> rowCandToXiPiKf;
    Produces<aod::HfOmegacKf> rowCandToOmegaPiKf;
    Produces<aod::HfCandToOmegaKaKf> rowCandToOmegaKaKf;
  } cursors;

  // Configurables
  struct : ConfigurableGroup {
    // Options for internal V0 building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
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
    // --------------------------------------------------------------------
    Configurable<float> dcabachtopvFromLF{"dcabachtopvFromLF", .1, "DCV Bach to PV"};
    Configurable<float> cascradiusFromLF{"cascradiusFromLF", .1, "DCV Bach to PV"};
    Configurable<float> casccospaFromLF{"casccospaFromLF", 0.95, "Cascade CosPA"};
    Configurable<float> dcacascdauFromLF{"dcacascdauFromLF", 1.0, "DCA cascade daughters"};
    Configurable<float> lambdaMassWindowFromLF{"lambdaMassWindowFromLF", 0.10, "Distance from Lambda mass(does not apply to KF path)"};

    // Options for internal cascade building - KF Building specifics
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // --------------------------------------------------------------------
    Configurable<bool> kfTuneForOmegaFromLF{"kfTuneForOmegaFromLF", false, "if enabled, take main cascade properties from omega fit instread of Xi fit(=default)"};
    Configurable<int> kfConstructMethodFromLF{"kfConstructMethodFromLF", 2, "2 : Daughter particle masses stay fixed in construction process"};
    Configurable<bool> kfUseV0MassConstraintFromLF{"kfUseV0MassConstraintFromLF", false, "KF : Use Lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraintFromLF{"kfUseCascadeMassConstraintFromLF", false, "KF : Use Cascade mass constraint - WARNING : Not adequate for inv mass analysis of Xi"};
    Configurable<bool> kfDoDCAFitterPreMinimV0FromLF{"kfDoDCAFitterPreMinimV0FromLF", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for V0"};
    Configurable<bool> kfDoDCAFitterPreMinimCascFromLF{"kfDoDCAFitterPreMinimCascFromLF", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for Xi"};
  } LFConfigs;

  // For cascade building using LF strangeness builder
  o2::pwglf::strangenessBuilderHelper straHelper;

  // Configurables
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
    // --------------------
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

  // For magnetic field
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // DCAFitter
  o2::vertexing::DCAFitterN<2> df;

  int runNumber{0};
  double magneticField{0.};

  enum CharmBaryonCandCounter { All = 0,
                                HfFlagPass,
                                CascReconstructed,
                                VertexFit };

  // Table aliases
  using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using TracksWCovIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
  using TracksWCovDcaExtraPidPrPiKa = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPr, aod::TracksPidPi, aod::TracksPidKa>;
  using TracksWCovExtraPidIU = soa::Join<aod::TracksIU, TracksCovIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  HistogramRegistry registry{"hists"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  HfEventSelection hfEvSel;

  // PDG Id of daughter tracks & V0s & cascades & charm baryons - Used in KFParticle
  int pdgIdOfV0DauPos, pdgIdOfV0DauNeg, pdgIdOfBach, pdgIdOfCharmBach;
  int pdgIdOfAntiV0DauPos, pdgIdOfAntiV0DauNeg, pdgIdOfAntiBach, pdgIdOfAntiCharmBach;
  int pdgIdOfV0, pdgIdOfCascade, pdgIdOfCharmBaryon;

  // Track PID - Used in DCAFitter
  int trackPidOfCascade;

  // Mass of daughter tracks & V0s & cascades & charm baryons;
  int massOfV0DauPos, massOfV0DauNeg, massOfBach, massOfCharmBach;
  int massOfV0, massOfCascade, massOfCharmBaryon;

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
    if (xipiEnabledDca || xipiEnabledKf) {
      pdgIdOfV0DauPos = kProton;
      pdgIdOfV0DauNeg = kPiMinus;
      pdgIdOfBach = kPiMinus;
      pdgIdOfCharmBach = kPiPlus;

      pdgIdOfAntiV0DauPos = kPiPlus;
      pdgIdOfAntiV0DauNeg = kProton;
      pdgIdOfAntiBach = kPiPlus;
      pdgIdOfAntiCharmBach = kPiMinus;

      pdgIdOfV0 = kLambda0;
      pdgIdOfCascade = kXiMinus;
      pdgIdOfCharmBaryon = kXiC0;

      trackPidOfCascade = o2::track::PID::XiMinus;
      massOfCharmBach = o2::constants::physics::MassPiPlus;
      massOfV0 = o2::constants::physics::MassLambda;
      massOfCascade = o2::constants::physics::MassXiMinus;
    } else if (omegapiEnabledDca || omegapiEnabledKf) {
      pdgIdOfV0DauPos = kProton;
      pdgIdOfV0DauNeg = kPiMinus;
      pdgIdOfBach = kKMinus;
      pdgIdOfCharmBach = kPiPlus;

      pdgIdOfAntiV0DauPos = kPiPlus;
      pdgIdOfAntiV0DauNeg = kProton;
      pdgIdOfAntiBach = kKPlus;
      pdgIdOfAntiCharmBach = kPiMinus;

      pdgIdOfV0 = kLambda0;
      pdgIdOfCascade = kOmegaMinus;
      pdgIdOfCharmBaryon = kOmegaC0;

      trackPidOfCascade = o2::track::PID::OmegaMinus;
      massOfCharmBach = o2::constants::physics::MassPiPlus;
      massOfV0 = o2::constants::physics::MassLambda;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
    } else if (omegakaEnabledDca || omegakaEnabledKf) {
      pdgIdOfV0DauPos = kProton;
      pdgIdOfV0DauNeg = kPiMinus;
      pdgIdOfBach = kKMinus;
      pdgIdOfCharmBach = kKPlus;

      pdgIdOfAntiV0DauPos = kPiPlus;
      pdgIdOfAntiV0DauNeg = kProton;
      pdgIdOfAntiBach = kKPlus;
      pdgIdOfAntiCharmBach = kKMinus;

      pdgIdOfV0 = kLambda0;
      pdgIdOfCascade = kOmegaMinus;
      pdgIdOfCharmBaryon = kOmegaC0;

      trackPidOfCascade = o2::track::PID::OmegaMinus;
      massOfCharmBach = o2::constants::physics::MassKPlus;
      massOfV0 = o2::constants::physics::MassLambda;
      massOfCascade = o2::constants::physics::MassOmegaMinus;
    }
    LOGF(info, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    LOGF(info, "PDG ID of V0 positive daughter: %d", pdgIdOfV0DauPos);
    LOGF(info, "PDG ID of V0 negative daughter: %d", pdgIdOfV0DauNeg);
    LOGF(info, "PDG ID of Bachelor: %d", pdgIdOfBach);
    LOGF(info, "PDG ID of Charm Bachelor: %d", pdgIdOfCharmBach);
    LOGF(info, "----------");
    LOGF(info, "PDG ID of anti V0 positive daughter: %d", pdgIdOfAntiV0DauPos);
    LOGF(info, "PDG ID of anti V0 negative daughter: %d", pdgIdOfAntiV0DauNeg);
    LOGF(info, "PDG ID of anti Bachelor: %d", pdgIdOfAntiBach);
    LOGF(info, "PDG ID of anti Charm Bachelor: %d", pdgIdOfAntiCharmBach);
    LOGF(info, "----------");
    LOGF(info, "PDG ID of V0: %d", pdgIdOfV0);
    LOGF(info, "PDG ID of Cascade: %d", pdgIdOfCascade);
    LOGF(info, "PDG ID of Charm Baryon: %d", pdgIdOfCharmBaryon);
    LOGF(info, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

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
    straHelper.fitter.setPropagateToPCA(configs.propagateToPCA);
    straHelper.fitter.setMaxR(configs.maxR);
    straHelper.fitter.setMaxDZIni(configs.maxDZIni);
    straHelper.fitter.setMinParamChange(configs.minParamChange);
    straHelper.fitter.setUseAbsDCA(configs.useAbsDCA);
    straHelper.fitter.setWeightedFinalPCA(configs.useWeightedFinalPCA);

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
    registry.add("ReconstructedDecayChannel", "DecayChannel", {kTH1F, {{3, 0.0, 3.0}}});
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

  ////////////////////////////////////////////////////////////
  //                                                        //
  //         Candidate reconstruction with DCAFitter        //
  //                                                        //
  ////////////////////////////////////////////////////////////

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
        // LOG(info) << "!This cascade cannot be rebuilt(cascade ID : " << cand.cascadeId() << "/ collision ID : " << collision.globalIndex();
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
      // V0 quantities from LF strangeness builder
      std::array<float, 3> vertexV0 = {straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2]};
      std::array<float, 3> pVecV0 = {straHelper.cascade.v0Momentum[0], straHelper.cascade.v0Momentum[1], straHelper.cascade.v0Momentum[2]};
      std::array<float, 3> pVecV0DauPos = {straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2]};
      std::array<float, 3> pVecV0DauNeg = {straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2]};

      // pseudo rapidity - V0 daughters
      // float pseudorapV0Dau0 = RecoDecay::eta(pVecV0DauPos);
      // float pseudorapV0Dau1 = RecoDecay::eta(pVecV0DauNeg);

      // Cascade quantities from LF strangeness builder
      int chargeCasc = straHelper.cascade.charge;
      std::array<float, 3> vertexCasc(straHelper.cascade.cascadePosition);
      std::array<float, 3> const pVecCasc(straHelper.cascade.cascadeMomentum);
      std::array<float, 21> covCasc = {0.};
      constexpr int NumCovElement = 6;
      constexpr int MomInd[NumCovElement] = {9, 13, 14, 18, 19, 20};
      for (int i = 0; i < NumCovElement; i++) {
        covCasc[MomInd[i]] = straHelper.cascade.covariance[MomInd[i]];
        covCasc[i] = straHelper.cascade.covariance[i];
      }

      // pseudo rapidity - cascade bachelor
      // float pseudorapCascBachelor = RecoDecay::eta(straHelper.cascade.bachelorMomentum);

      //------------------------------Create cascade track------------------------------

      o2::track::TrackParCov trackCasc;
      if (chargeCasc < 0) { // Xi- or Omega-
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else if (chargeCasc > 0) { // Xi+ or Omega+
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else {
        continue;
      }

      trackCasc.setAbsCharge(1);
      trackCasc.setPID(trackPidOfCascade);

      //------------------------------Fit SV & Create Charm Baryon track------------------------------

      // Perform secondary vertex fitting
      auto trackCharmBachelor = tracks.rawIteratorAt(cand.prong0Id());
      auto trackParCovCharmBachelor = getTrackParCov(trackCharmBachelor);
      try {
        if (df.process(trackCasc, trackParCovCharmBachelor) == 0) {
          continue;
        }
      } catch (std::runtime_error& e) {
        LOG(error) << "Execption caught in charm DCA Fitter process call: " << e.what();
        continue;
      }

      if (configs.fillHistograms) {
        hCandCounter->Fill(VertexFit);
      }

      //------------------------------Calculate physical properties-----------------------------

      // get track momenta
      std::array<float, 3> pVecCascAsD, pVecCharmBachAsD;
      df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
      df.getTrack(1).getPxPyPzGlo(pVecCharmBachAsD);

      std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecCharmBachAsD[0], pVecCascAsD[1] + pVecCharmBachAsD[1], pVecCascAsD[2] + pVecCharmBachAsD[2]};
      std::array<float, 3> pVecBach = {straHelper.cascade.bachelorMomentum[0], straHelper.cascade.bachelorMomentum[1], straHelper.cascade.bachelorMomentum[2]};

      // get PV Properties
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

      // get SV Properties
      auto const& secondaryVertex = df.getPCACandidate();
      auto covMatrixSV = df.calcPCACovMatrixFlat();

      // DCAxy and DCAz. Computed with propagatToDCABxByBz method
      auto trackParCovV0DauPos = getTrackParCov(posTrack);
      auto trackParCovV0DauNeg = getTrackParCov(negTrack);
      auto trackParCovBach = getTrackParCov(bachTrack);

      o2::dataformats::DCA impactParameterV0DauPos, impactParameterV0DauNeg, impactParameterBach;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0DauPos, 2.f, matCorr, &impactParameterV0DauPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0DauNeg, 2.f, matCorr, &impactParameterV0DauNeg);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovBach, 2.f, matCorr, &impactParameterBach);

      float dcaxyV0DauPos = impactParameterV0DauPos.getY();
      float dcaxyV0DauNeg = impactParameterV0DauNeg.getY();
      float dcaxyBach = impactParameterBach.getY();
      float dcazV0DauPos = impactParameterV0DauPos.getZ();
      float dcazV0DauNeg = impactParameterV0DauNeg.getZ();
      float dcazBach = impactParameterBach.getZ();

      // get impact parameter. Compute with propagateToDCABxByBz method
      o2::dataformats::DCA impactParameterCasc, impactParameterCharmBach;
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc); // trackCasc is TrackParCov object
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCharmBachelor, 2.f, matCorr, &impactParameterCharmBach);
      // float impactParCharmBachFromCharmBayonXY = impactParameterCharmBach.getY();
      // float impactParCharmBachFromCharmBayonZ = impactParameterCharmBach.getZ();

      // get v0 invariant mass - from LF Table
      float mLambda = straHelper.v0.massLambda; // from LF Table

      // get Casc mass - from LF Table
      float mCasc = 0.;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        mCasc = straHelper.cascade.massXi;
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        mCasc = straHelper.cascade.massOmega;
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        mCasc = straHelper.cascade.massOmega;
      }

      // get Charm baryon invariant mass
      auto arrMomenta = std::array{pVecCascAsD, pVecCharmBachAsD};
      float massCharmBaryonCand = RecoDecay::m(arrMomenta, std::array{massOfCascade, massOfCharmBach});
      if (configs.fillHistograms) {
        hInvMassCharmBaryon->Fill(massCharmBaryonCand);
      }

      // calculate cosine of pointing angle
      std::array<float, 3> vtxCoordCharmBaryon = df.getPCACandidatePos();
      float cpaV0 = RecoDecay::cpa(pvCoord, vertexV0, pVecV0);
      float cpaCasc = RecoDecay::cpa(pvCoord, vertexCasc, pVecCasc);
      float cpaCharmBaryon = RecoDecay::cpa(pvCoord, vtxCoordCharmBaryon, pVecCharmBaryon);
      float cpaxyV0 = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cpaxyCasc = RecoDecay::cpaXY(pvCoord, vertexCasc, pVecCasc);
      float cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, vtxCoordCharmBaryon, pVecCharmBaryon);

      // calculate decay length
      float decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);
      float decLenCasc = RecoDecay::distance(vtxCoordCharmBaryon, vertexCasc);
      float decLenCharmBaryon = RecoDecay::distance(pvCoord, vtxCoordCharmBaryon);
      float phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      // calcuate ctau
      float ctV0 = RecoDecay::ct(pVecV0, decLenV0, MassLambda0);
      float ctCasc = RecoDecay::ct(pVecCasc, decLenCasc, massOfCascade);
      float ctOmegac0 = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassOmegaC0);
      float ctXic0 = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, MassXiC0);

      // get eta
      float etaV0DauPos = posTrack.eta();
      float etaV0DauNeg = negTrack.eta();
      float etaBach = bachTrack.eta();
      float etaCharmBach = trackCharmBachelor.eta();
      float etaV0 = RecoDecay::eta(pVecV0);
      float etaCasc = RecoDecay::eta(pVecCasc);
      float etaCharmBaryon = RecoDecay::eta(pVecCharmBaryon);

      // DCA between daughters
      float dcaV0Dau = straHelper.cascade.v0DaughterDCA;
      float dcaCascDau = straHelper.cascade.cascadeDaughterDCA;
      float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

      //------------------------------Fill QA histograms-----------------------------

      //------------------------------Fill the table-----------------------------
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        cursors.rowCandToXiPi(collision.globalIndex(),
                              pvCoord[0], pvCoord[1], pvCoord[2],
                              secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                              vertexCasc[0], vertexCasc[1], vertexCasc[2],
                              vertexV0[0], vertexV0[1], vertexV0[2],
                              bachTrack.sign(),
                              covMatrixSV[0], covMatrixSV[1], covMatrixSV[2], covMatrixSV[3], covMatrixSV[4], covMatrixSV[5],
                              pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                              pVecCascAsD[0], pVecCascAsD[1], pVecCascAsD[2],
                              pVecCharmBachAsD[0], pVecCharmBachAsD[1], pVecCharmBachAsD[2],
                              pVecV0[0], pVecV0[1], pVecV0[2],
                              pVecBach[0], pVecBach[1], pVecBach[2],
                              pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                              pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                              impactParameterCasc.getY(), impactParameterCharmBach.getY(),
                              impactParameterCasc.getZ(), impactParameterCharmBach.getZ(),
                              std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBach.getSigmaY2()),
                              cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                              cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                              mLambda, mCasc, massCharmBaryonCand,
                              cpaV0, cpaCharmBaryon, cpaCasc,
                              cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                              ctOmegac0, ctCasc, ctV0, ctXic0,
                              etaV0DauPos, etaV0DauNeg, etaBach, etaCharmBach, etaCharmBaryon, etaCasc, etaV0,
                              dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                              dcazV0DauPos, dcazV0DauNeg, dcazBach,
                              dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                              decLenCharmBaryon, decLenCasc, decLenV0, errorDecayLength, errorDecayLengthXY);
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        cursors.rowCandToOmegaPi(collision.globalIndex(),
                                 pvCoord[0], pvCoord[1], pvCoord[2],
                                 secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                                 vertexCasc[0], vertexCasc[1], vertexCasc[2],
                                 vertexV0[0], vertexV0[1], vertexV0[2],
                                 bachTrack.sign(),
                                 covMatrixSV[0], covMatrixSV[1], covMatrixSV[2], covMatrixSV[3], covMatrixSV[4], covMatrixSV[5],
                                 pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                                 pVecCascAsD[0], pVecCascAsD[1], pVecCascAsD[2],
                                 pVecCharmBachAsD[0], pVecCharmBachAsD[1], pVecCharmBachAsD[2],
                                 pVecV0[0], pVecV0[1], pVecV0[2],
                                 pVecBach[0], pVecBach[1], pVecBach[2],
                                 pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                                 pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                                 impactParameterCasc.getY(), impactParameterCharmBach.getY(),
                                 impactParameterCasc.getZ(), impactParameterCharmBach.getZ(),
                                 std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBach.getSigmaY2()),
                                 cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                                 cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                                 mLambda, mCasc, massCharmBaryonCand,
                                 cpaV0, cpaCharmBaryon, cpaCasc,
                                 cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                                 ctOmegac0, ctCasc, ctV0,
                                 etaV0DauPos, etaV0DauNeg, etaBach, etaCharmBach, etaCharmBaryon, etaCasc, etaV0,
                                 dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                                 dcazV0DauPos, dcazV0DauNeg, dcazBach,
                                 dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                                 decLenCharmBaryon, decLenCasc, decLenV0, errorDecayLength, errorDecayLengthXY, cand.hfflag());
      } else {
        cursors.rowCandToOmegaKa(collision.globalIndex(),
                                 pvCoord[0], pvCoord[1], pvCoord[2],
                                 secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                                 vertexCasc[0], vertexCasc[1], vertexCasc[2],
                                 vertexV0[0], vertexV0[1], vertexV0[2],
                                 bachTrack.sign(),
                                 covMatrixSV[0], covMatrixSV[1], covMatrixSV[2], covMatrixSV[3], covMatrixSV[4], covMatrixSV[5],
                                 pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                                 pVecCascAsD[0], pVecCascAsD[1], pVecCascAsD[2],
                                 pVecCharmBachAsD[0], pVecCharmBachAsD[1], pVecCharmBachAsD[2],
                                 pVecV0[0], pVecV0[1], pVecV0[2],
                                 pVecBach[0], pVecBach[1], pVecBach[2],
                                 pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                                 pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                                 impactParameterCasc.getY(), impactParameterCharmBach.getY(),
                                 impactParameterCasc.getZ(), impactParameterCharmBach.getZ(),
                                 std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBach.getSigmaY2()),
                                 cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                                 cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                                 mLambda, mCasc, massCharmBaryonCand,
                                 cpaV0, cpaCharmBaryon, cpaCasc,
                                 cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                                 ctOmegac0, ctCasc, ctV0,
                                 etaV0DauPos, etaV0DauNeg, etaBach, etaCharmBach, etaCharmBaryon, etaCasc, etaV0,
                                 dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                                 dcazV0DauPos, dcazV0DauNeg, dcazBach,
                                 dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                                 decLenCharmBaryon, decLenCasc, decLenV0, errorDecayLength, errorDecayLengthXY);
      }
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
        LOG(info) << "This cascade cannot be rebuilt";
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

      //----------Create charm bayron as KF Partible object starting from V0----------

      // Create KFParticle object of V0 Daughter & Bachelor
      const KFPTrack kfTrack0 = createKFPTrackFromTrack(posTrack);
      const KFPTrack kfTrack1 = createKFPTrackFromTrack(negTrack);
      const KFPTrack kfTrackBach = createKFPTrackFromTrack(bachTrack);

      bool isAnti = (bachTrack.signed1Pt() > 0 ? true : false);

      KFParticle kfPos(kfTrack0, (isAnti ? pdgIdOfAntiV0DauPos : pdgIdOfV0DauPos));
      KFParticle kfNeg(kfTrack1, (isAnti ? pdgIdOfAntiV0DauNeg : pdgIdOfV0DauNeg));
      KFParticle kfBach(kfTrackBach, (isAnti ? pdgIdOfAntiBach : pdgIdOfBach));
      KFParticle kfBachRej(kfTrackBach, (isAnti ? pdgIdOfAntiBach : pdgIdOfBach)); // Rej -> Used for Omegac0->OmegaPi only

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
      KFParticle kfV0MassConstrained = kfV0;
      kfV0MassConstrained.SetNonlinearMassConstraint(massOfV0);
      if (configs.kfUseV0MassConstraint) {
        kfV0 = kfV0MassConstrained;
      }
      kfV0.TransportToDecayVertex();

      //~~~~~~~Construct cascade with KF~~~~~~~
      const KFParticle* cascDaughters[2] = {&kfBach, &kfV0};
      const KFParticle* cascDaughtersRej[2] = {&kfBachRej, &kfV0};
      KFParticle kfCasc, kfCascRej;

      kfCasc.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfCasc.Construct(cascDaughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Cascade: " << e.what();
      }

      float massCasc, sigMassCasc, massCascRej, sigMassCascRej;
      kfCasc.GetMass(massCasc, sigMassCasc);

      if (sigMassCasc <= 0) {
        continue;
      }
      if (std::abs(massCasc - massCasc) > configs.massToleranceCascade) {
        continue;
      }
      if (kfCasc.GetNDF() <= 0 || kfCasc.GetChi2() <= 0) {
        continue;
      }

      // perform cascade building on casc_rej - only for Omega
      if constexpr (decayChannel != hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        kfCascRej.SetConstructMethod(configs.kfConstructMethod);
        try {
          kfCascRej.Construct(cascDaughtersRej, 2);
        } catch (std::runtime_error& e) {
          LOG(debug) << "Failed to construct Cascade_rej: " << e.what();
        }

        kfCascRej.GetMass(massCascRej, sigMassCascRej);
      }

      // Set mass constraint to cascade
      KFParticle kfCascMassConstrained = kfCasc;
      kfCascMassConstrained.SetNonlinearMassConstraint(massOfCascade);
      if (configs.kfUseCascadeMassConstraint) {
        kfCasc = kfCascMassConstrained;
      }
      kfCasc.TransportToDecayVertex();

      //~~~~~~~Construct Charm Baryon with KF~~~~~~~
      auto trackCharmBachelor = tracks.rawIteratorAt(cand.prong0Id());
      const KFPTrack kfTrackCharmBach = createKFPTrackFromTrack(trackCharmBachelor);
      const KFParticle kfCharmBach(kfTrackCharmBach, (isAnti ? pdgIdOfAntiCharmBach : pdgIdOfCharmBach));
      const KFParticle* charmBaryonDaughters[2] = {&kfCharmBach, &kfCasc};

      KFParticle kfCharmBaryon;
      kfCharmBaryon.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfCharmBaryon.Construct(charmBaryonDaughters, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Charm baryon: " << e.what();
      }

      float massCharmBaryon, sigMassCharmBaryon;
      kfCharmBaryon.GetMass(massCharmBaryon, sigMassCharmBaryon);
      if (sigMassCharmBaryon <= 0) {
        continue;
      }
      if (kfCharmBaryon.GetNDF() <= 0 || kfCharmBaryon.GetChi2() <= 0) {
        continue;
      }
      kfCharmBaryon.TransportToDecayVertex();
      if (configs.fillHistograms) {
        hCandCounter->Fill(VertexFit);
        hInvMassCharmBaryon->Fill(massCharmBaryon);
      }

      // Set production vertex
      // PV
      const KFPVertex kfVertex = createKFPVertexFromCollision(collision);
      const KFParticle kfPv(kfVertex);
      const KFParticle kfPosOrigin = kfPos;
      const KFParticle kfNegOrigin = kfNeg;

      // To V0
      kfPos.SetProductionVertex(kfV0);
      kfNeg.SetProductionVertex(kfV0);

      // To Casc
      KFParticle kfBachToCasc = kfBach;
      KFParticle kfV0ToCasc = kfV0;
      kfBach.SetProductionVertex(kfCasc);
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

      //----------Reconstruct information after vertex fit----------
      std::array<float, 3> vertexV0 = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      std::array<float, 3> vertexCasc = {kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()};

      std::array<float, 3> pVecV0DauPos = {kfPos.GetPx(), kfPos.GetPy(), kfPos.GetPz()};
      std::array<float, 3> pVecV0DauNeg = {kfNeg.GetPx(), kfNeg.GetPy(), kfNeg.GetPz()};
      std::array<float, 3> pVecV0 = {kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz()};
      std::array<float, 3> pVecBach = {kfBachToCasc.GetPx(), kfBachToCasc.GetPy(), kfBachToCasc.GetPz()};
      std::array<float, 3> pVecCharmBachelorAsD = {kfCharmBachToCharmBaryon.GetPx(), kfCharmBachToCharmBaryon.GetPy(), kfCharmBachToCharmBaryon.GetPz()};
      std::array<float, 3> pVecCharmBaryon = {kfCharmBaryon.GetPx(), kfCharmBaryon.GetPy(), kfCharmBaryon.GetPy()};

      auto* covVtxCharmBaryon = kfCharmBaryon.CovarianceMatrix();
      float covMatrixPv[6];
      kfVertex.GetCovarianceMatrix(covMatrixPv);

      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

      auto trackParCovV0DauPos = getTrackParCovFromKFP(kfPos, kfPos.GetPDG(), 1);
      auto trackParCovV0DauNeg = getTrackParCovFromKFP(kfNeg, kfNeg.GetPDG(), -1);
      auto trackParCovBach = getTrackParCovFromKFP(kfBachToCasc, kfBachToCasc.GetPDG(), (isAnti ? 1 : -1));
      auto trackParCovCharmBach = getTrackParCovFromKFP(kfCharmBachToCharmBaryon, kfCharmBachToCharmBaryon.GetPDG(), (isAnti ? -1 : 1));
      auto trackParCovCasc = getTrackParCovFromKFP(kfCascToCharmBaryon, kfCascToCharmBaryon.GetPDG(), (isAnti ? 1 : -1));
      trackParCovV0DauPos.setAbsCharge(1);
      trackParCovV0DauNeg.setAbsCharge(1);
      trackParCovBach.setAbsCharge(1);
      trackParCovCharmBach.setAbsCharge(1);
      trackParCovCasc.setAbsCharge(1);

      //----------Calculate physical quantities and fill candidate table----------

      // impact parameters
      std::array<float, 2> impactParameterV0DauPos;
      std::array<float, 2> impactParameterV0DauNeg;
      std::array<float, 2> impactParameterBach;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0DauPos, 2.f, matCorr, &impactParameterV0DauPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0DauNeg, 2.f, matCorr, &impactParameterV0DauNeg);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovBach, 2.f, matCorr, &impactParameterBach);
      float dcaxyV0DauPos = impactParameterV0DauPos[0];
      float dcaxyV0DauNeg = impactParameterV0DauNeg[0];
      float dcaxyBach = impactParameterBach[0];
      float dcazV0DauPos = impactParameterV0DauPos[1];
      float dcazV0DauNeg = impactParameterV0DauNeg[1];
      float dcazBach = impactParameterBach[1];

      o2::dataformats::DCA impactParameterCasc, impactParameterCharmBachelor;
      auto primaryVertex = getPrimaryVertex(collision);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCasc, 2.f, matCorr, &impactParameterCasc);
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCharmBach, 2.f, matCorr, &impactParameterCharmBachelor);

      // get Chi2Topo/NDF
      float chi2NdfTopoV0ToPv = kfV0ToPv.GetChi2() / kfV0ToPv.GetNDF();
      float chi2NdfTopoCascToPv = kfCascToPv.GetChi2() / kfCascToPv.GetNDF();
      float chi2NdfTopoCharmBachToPv = kfCharmBachToPv.GetChi2() / kfCharmBachToPv.GetNDF();
      float chi2NdfTopoCharmBaryonToPv = kfCharmBaryonToPv.GetChi2() / kfCharmBaryonToPv.GetNDF();
      float chi2NdfTopoBachToCasc = kfBachToCasc.GetChi2() / kfBachToCasc.GetNDF();
      float chi2NdfTopoV0ToCasc = kfV0ToCasc.GetChi2() / kfV0ToCasc.GetNDF();
      float chi2NdfTopoCharmBachToCharmBaryon = kfCharmBachToCharmBaryon.GetChi2() / kfCharmBachToCharmBaryon.GetChi2();
      float chi2NdfTopoCascToCharmBaryon = kfCascToCharmBaryon.GetChi2() / kfCascToCharmBaryon.GetChi2();

      // get ldl
      float ldlV0 = ldlFromKF(kfV0, kfPv);
      float ldlCasc = ldlFromKF(kfCasc, kfPv);
      float ldlCharmBaryon = ldlFromKF(kfCharmBaryon, kfPv);

      // get DCAs
      float kfDcaV0Daughters = kfNeg.GetDistanceFromParticle(kfPos);
      float kfDcaCascDaughters = kfBachToCasc.GetDistanceFromParticle(kfV0ToCasc);
      float kfDcaCharmBaryonDaughters = kfCharmBachToCharmBaryon.GetDistanceFromParticle(kfCascToCharmBaryon);
      float kfDcaXYCharmBachelorToPv = kfCharmBachToCharmBaryon.GetDistanceFromVertexXY(kfPv);
      float kfDcaXYCascToPv = kfCascToCharmBaryon.GetDistanceFromVertexXY(kfPv);

      // get decay length - In XY
      float decayLXYV0, errDecayLXYV0;
      kfV0ToCasc.GetDecayLengthXY(decayLXYV0, errDecayLXYV0);

      float decayLXYCasc, errDecayLXYCasc;
      kfCascToCharmBaryon.GetDecayLengthXY(decayLXYCasc, errDecayLXYCasc);

      float decayLXYCharmBaryon, errDecayLXYCharmBaryon;
      kfCharmBaryonToPv.GetDecayLengthXY(decayLXYCharmBaryon, errDecayLXYCharmBaryon);

      // get decay length - In XYZ
      float decayLV0 = RecoDecay::distance(std::array<float, 3>{kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()}, std::array<float, 3>{kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()});
      float decayLCasc = RecoDecay::distance(std::array<float, 3>{kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()}, std::array<float, 3>{kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()});
      float decayLCharmBaryon = RecoDecay::distance(std::array<float, 3>{collision.posX(), collision.posY(), collision.posZ()}, std::array<float, 3>{kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()});

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array<float, 3>{kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()}, std::array<float, 3>{kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()}, phiCharmBaryon, thetaCharmBaryon);
      float errDecayLCharmBaryon = std::sqrt(getRotatedCovMatrixXX(covMatrixPv, phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));

      // get cosine of pointing angle
      float cosPaV0ToPv = cpaFromKF(kfV0, kfPv);
      float cosPaCascToPv = cpaFromKF(kfCasc, kfPv);
      float cosPaCharmBaryonToPv = cpaFromKF(kfCharmBaryon, kfPv);

      float cosPaXYV0ToPv = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cosPaXYCascToPv = cpaXYFromKF(kfCasc, kfPv);
      float cosPaXYCharmBaryonToPv = cpaXYFromKF(kfCharmBaryon, kfPv);

      float cosPaV0ToCasc = cpaFromKF(kfV0, kfCasc);
      float cosPaCascToCharmBaryon = cpaFromKF(kfCasc, kfCharmBaryon);
      float cosPaXYV0ToCasc = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
      float cosPaXYCascToCharmBaryon = cpaXYFromKF(kfCasc, kfCharmBaryon);

      float deviationCharmBachToPv = kfCalculateChi2ToPrimaryVertex(kfCharmBaryon, kfPv); // -> For Omegac0

      // KF pT, eta
      float ptCasc = kfCascToCharmBaryon.GetPt();
      float ptCharmBachelor = kfCharmBachToCharmBaryon.GetPt();
      float ptCharmBaryon = kfCharmBaryon.GetPt();
      float yCharmBaryon = kfCharmBaryon.GetRapidity();

      // get KF cosThetaStar
      float cosThetaStarCharmBachelor;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        cosThetaStarCharmBachelor = cosThetaStarFromKF(0, pdgIdOfCharmBaryon, pdgIdOfCharmBach, pdgIdOfCascade, kfCharmBachToCharmBaryon, kfCascToCharmBaryon);
      }
      float cosThetaStarKaFromOmegac0, cosThetaStarKaFromXic0; // -> Only requeted to be calculated and filled for Xic0/Omegac0 -> Omega K
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        cosThetaStarKaFromOmegac0 = cosThetaStarFromKF(0, 4332, 321, 3334, kfCharmBachToCharmBaryon, kfCascToCharmBaryon);
        cosThetaStarKaFromXic0 = cosThetaStarFromKF(0, 4132, 321, 3334, kfCharmBachToCharmBaryon, kfCascToCharmBaryon);
      }

      // KF ct
      float ctV0 = kfV0ToCasc.GetLifeTime();
      float ctCasc = kfCascToCharmBaryon.GetLifeTime();
      float ctCharmBaryon = kfCharmBaryonToPv.GetLifeTime();

      //------------------------------Calculate physical quantities and fill candidate table------------------------------

      //------------------------------Fill the table------------------------------
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        cursors.rowCandToXiPiKf(collision.globalIndex(),                     // Global index of collision
                                pvCoord[0], pvCoord[1], pvCoord[2],          // coordination of PV
                                vertexCasc[0], vertexCasc[1], vertexCasc[2], // Decay position of kfCasc
                                vertexV0[0], vertexV0[1], vertexV0[2],
                                bachTrack.sign(),
                                covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                                pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2], // x, y, z momentum of charm baryon
                                kfCascToCharmBaryon.GetPx(), kfCascToCharmBaryon.GetPy(), kfCascToCharmBaryon.GetPz(),
                                pVecCharmBachelorAsD[0], pVecCharmBachelorAsD[1], pVecCharmBachelorAsD[2],
                                pVecV0[0], pVecV0[1], pVecV0[2],
                                pVecBach[0], pVecBach[1], pVecBach[2],
                                pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                                pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                                cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                                cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                                massLam, massCasc, massCharmBaryon,
                                cosPaV0ToPv, cosPaCascToPv,
                                ctCasc, ctV0, ctCharmBaryon,
                                kfPos.GetEta(), kfNeg.GetEta(), kfBach.GetEta(), kfCharmBachToCharmBaryon.GetEta(),
                                kfCharmBaryon.GetEta(), kfCasc.GetEta(), kfV0.GetEta(),
                                dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                                kfDcaCascDaughters, kfDcaV0Daughters, kfDcaCharmBaryonDaughters,
                                kfDcaXYCharmBachelorToPv, kfDcaXYCascToPv,
                                kfV0.GetChi2(), kfCasc.GetChi2(), kfCharmBaryon.GetChi2(), kfV0MassConstrained.GetChi2(), kfCascMassConstrained.GetChi2(),
                                ldlV0, ldlCasc, // ldlCharmBaryon,
                                chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoCharmBachToPv, chi2NdfTopoCharmBaryonToPv,
                                chi2NdfTopoV0ToCasc, chi2NdfTopoCascToCharmBaryon,
                                decayLXYV0, decayLXYCasc, decayLXYCharmBaryon,
                                cosPaV0ToCasc, cosPaCascToCharmBaryon, // cosPaXYV0ToCasc, cosPaXYCascToCharmBaryon,
                                yCharmBaryon,                          // ptCharmBachelor, ptCharmBaryon,
                                cosThetaStarCharmBachelor,
                                kfV0.GetNDF(), kfCasc.GetNDF(), kfCharmBaryon.GetNDF(), kfV0MassConstrained.GetNDF(), kfCascMassConstrained.GetNDF(),
                                kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.GetChi2() / kfCasc.GetNDF(), kfCharmBaryon.GetChi2() / kfCharmBaryon.GetNDF(), kfV0MassConstrained.GetChi2() / kfV0MassConstrained.GetNDF(), kfCascMassConstrained.GetChi2() / kfCascMassConstrained.GetNDF());
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        cursors.rowCandToOmegaPi(collision.globalIndex(),
                                 pvCoord[0], pvCoord[1], pvCoord[2],
                                 0.f, 0.f, 0.f, // -> vertexCharmBaryonFromFitter. For KF, this is 0
                                 vertexCasc[0], vertexCasc[1], vertexCasc[2],
                                 vertexV0[0], vertexV0[1], vertexV0[2],
                                 trackCharmBachelor.sign(),
                                 covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
                                 kfCharmBaryon.GetPx(), kfCharmBaryon.GetPy(), kfCharmBaryon.GetPz(),
                                 kfCascToCharmBaryon.GetPx(), kfCascToCharmBaryon.GetPy(), kfCascToCharmBaryon.GetPz(),
                                 kfCharmBachToCharmBaryon.GetPx(), kfCharmBachToCharmBaryon.GetPy(), kfCharmBachToCharmBaryon.GetPz(),
                                 pVecV0[0], pVecV0[1], pVecV0[2],
                                 pVecBach[0], pVecBach[1], pVecBach[2],
                                 pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                                 pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                                 impactParameterCasc.getY(), impactParameterCharmBachelor.getY(),
                                 impactParameterCasc.getZ(), impactParameterCharmBachelor.getZ(),
                                 std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                                 cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                                 cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                                 massLam, massCasc, massCharmBaryon,
                                 cosPaV0ToPv, cosPaCharmBaryonToPv, cosPaCascToPv, cosPaXYV0ToPv, cosPaXYCharmBaryonToPv, cosPaXYCascToPv,
                                 ctCharmBaryon, ctCasc, ctV0,
                                 kfPos.GetEta(), kfNeg.GetEta(), kfBach.GetEta(), kfCharmBachToCharmBaryon.GetEta(),
                                 kfCharmBaryon.GetEta(), kfCasc.GetEta(), kfV0.GetEta(),
                                 dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                                 dcazV0DauPos, dcazV0DauNeg, dcazBach,
                                 kfDcaCascDaughters, straHelper.cascade.v0DaughterDCA, kfDcaCharmBaryonDaughters,
                                 decayLCharmBaryon, decayLCasc, decayLV0, errDecayLCharmBaryon, errDecayLXYCharmBaryon, cand.hfflag());

        cursors.rowCandToOmegaPiKf(kfDcaXYCharmBachelorToPv, kfDcaXYCascToPv,
                                   kfV0.GetChi2(), kfCasc.GetChi2(), kfCharmBaryon.GetChi2(), kfV0MassConstrained.GetChi2(), kfCascMassConstrained.GetChi2(),
                                   ldlV0, ldlCasc, ldlCharmBaryon,
                                   chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoCharmBachToPv, chi2NdfTopoCharmBaryonToPv, deviationCharmBachToPv,
                                   chi2NdfTopoV0ToCasc, chi2NdfTopoCascToCharmBaryon,
                                   decayLXYV0, decayLXYCasc, decayLXYCharmBaryon,
                                   cosPaV0ToCasc, cosPaCascToCharmBaryon, cosPaXYV0ToCasc, cosPaXYCascToCharmBaryon,
                                   kfCharmBaryon.GetRapidity(), ptCharmBachelor, ptCharmBaryon,
                                   cosThetaStarKaFromOmegac0,
                                   kfV0.GetNDF(), kfCasc.GetNDF(), kfCharmBaryon.GetNDF(), kfV0MassConstrained.GetNDF(), kfCascMassConstrained.GetNDF(),
                                   kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.Chi2() / kfCasc.GetNDF(), kfCharmBaryon.GetChi2() / kfCharmBaryon.GetNDF(),
                                   kfV0MassConstrained.GetChi2() / kfV0.GetNDF(), kfCascMassConstrained.Chi2() / kfCasc.GetNDF(),
                                   massCascRej);

      } else {
        cursors.rowCandToOmegaKaKf(collision.globalIndex(),
                                   collision.posX(), collision.posY(), collision.posZ(),
                                   kfPv.GetX(), kfPv.GetY(), kfPv.GetZ(),
                                   straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2],
                                   straHelper.cascade.v0Momentum[0], straHelper.cascade.v0Momentum[1], straHelper.cascade.v0Momentum[2],
                                   straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2],
                                   straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2],
                                   kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(),
                                   kfV0.GetPx(), kfV0.GetPy(), kfV0.GetPz(),
                                   kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ(),
                                   kfCasc.GetPx(), kfCasc.GetPy(), kfCasc.GetPz(),
                                   kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ(),
                                   kfCharmBaryon.GetPx(), kfCharmBaryon.GetPy(), kfCharmBaryon.GetPz(),
                                   straHelper.cascade.charge,
                                   kfPos.GetEta(), kfNeg.GetEta(), kfBach.GetEta(), kfCharmBach.GetEta(), kfV0.GetEta(), kfCasc.GetEta(), kfCharmBaryon.GetEta(), kfCharmBaryon.GetRapidity(),
                                   impactParameterCharmBachelor.getY(), std::sqrt(impactParameterCharmBachelor.getSigmaY2()), impactParameterCasc.getY(), std::sqrt(impactParameterCasc.getSigmaY2()),
                                   kfDcaV0Daughters, kfDcaCascDaughters, kfDcaCharmBaryonDaughters,
                                   cosPaV0ToPv, cosPaCascToPv, cosPaCharmBaryonToPv, cosPaXYV0ToPv, cosPaXYCascToPv, cosPaXYCharmBaryonToPv, cosPaV0ToCasc, cosPaCascToCharmBaryon, cosPaXYV0ToCasc, cosPaXYCascToCharmBaryon,
                                   kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.GetChi2() / kfCasc.GetNDF(), kfCharmBaryon.GetChi2() / kfCharmBaryon.GetNDF(),
                                   kfV0MassConstrained.GetChi2() / kfV0MassConstrained.GetNDF(), kfCascMassConstrained.GetChi2() / kfCascMassConstrained.GetNDF(),
                                   chi2NdfTopoV0ToCasc, chi2NdfTopoBachToCasc, chi2NdfTopoCharmBachToCharmBaryon, chi2NdfTopoCascToCharmBaryon, // Topological constraints to mother
                                   chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoCharmBachToPv, chi2NdfTopoCharmBaryonToPv,                // Topological constraints to PV
                                   ldlV0, ldlCasc, ldlCharmBaryon,
                                   decayLXYV0, decayLXYCasc, decayLXYCharmBaryon,
                                   massLam, sigMassLam, massCasc, sigMassCasc, massCascRej, sigMassCascRej, massCharmBaryon, sigMassCharmBaryon,
                                   ptCharmBaryon, ptCharmBachelor, ptCasc,
                                   cosThetaStarKaFromOmegac0, cosThetaStarKaFromXic0, ctV0, ctCasc, ctCharmBaryon,
                                   cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(), cand.cascadeId(), cand.prong0Id(), trackCharmBachelor.globalIndex());
      }
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

  void processToXiPiWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
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

  void processToXiPiWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
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

  void processToOmegaPiWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
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

  void processToOmegaPiWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             aod::TracksWCovDca const& tracks,
                                             TracksWCovIU const& lfTracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaPi);
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

  void processToOmegaKaWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
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

  void processToOmegaKaWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             aod::TracksWCovDca const& tracks,
                                             TracksWCovIU const& lfTracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, lfTracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaPi);
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

  void processToXiPiWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
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

  void processToXiPiWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
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

  void processToOmegaPiWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
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

  void processToOmegaPiWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
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

  void processToOmegaKaWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
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

  void processToOmegaKaWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
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
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processCollisionsNoCent, "Collision monitoring - No Centrality", true);

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
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processCollisionsCentFT0C, "Collision monitoring - Centrality selection with FT0C", false);

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
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processCollisionsCentFT0M, "Collision monitoring - Centrality selection with FT0M", false);
};

struct HfCandidateCreatorXic0Omegac0QaMc {

  // Cursor to fill tables
  struct : ProducesGroup {

    Produces<aod::HfXicToXiPiMCRec> rowMcMatchRecXicToXiPi;
    Produces<aod::HfXicToXiPiMCGen> rowMcMatchGenXicToXiPi;
    Produces<aod::HfOmegacToXiPiMCRec> rowMcMatchRecOmegacToXiPi;
    Produces<aod::HfOmegacToXiPiMCGen> rowMcMatchGenOmegacToXiPi;
    Produces<aod::HfToOmegaPiMCRec> rowMcMatchRecToOmegaPi;
    Produces<aod::HfToOmegaPiMCGen> rowMcMatchGenToOmegaPi;
    Produces<aod::HfToOmegaKMCRec> rowMcMatchRecToOmegaKa;
    Produces<aod::HfToOmegaKMCGen> rowMcMatchGenToOmegaKa;

  } cursors;

  // Configurables
  struct : ConfigurableGroup {

    Configurable<bool> rejectBackground{"rejectBackground", true, "Reject particles from background events"}; // -> Used for only Xic0
    Configurable<bool> acceptTrackInteractionWithMaterial{"acceptTrackInteractionWithMaterial", false, "Accept candidates with final daughters interacting with materials"};
    Configurable<bool> fillMcHistograms{"fillMcHistograms", true, "Fill validation plots"};
    Configurable<bool> fillResidualTable{"fillResidualTable", true, "Fill table contaning residuals and pulls of PV and SV"};
    // Configurable<bool> matchDecayedPions{"matchedDecayedPions", true, "Match also candidates with daughter pion tracks that decay with kinked toploogy"};

  } configs;

  enum McMatchFlag : uint8_t {
    None = 0,
    CharmBaryonUnmatched,
    CascUnmatched,
    V0Unmatched,
    NumberMcMatchFlag
  };

  std::array<int, 4> pdgOfCharmBaryon{+kXiC0, +kOmegaC0, +kOmegaC0, +kOmegaC0};
  std::array<int, 4> pdgOfCascade{+kXiMinus, +kXiMinus, +kOmegaMinus, +kOmegaMinus};
  std::array<int, 4> pdgOfCharmBachelor{+kPiPlus, +kPiPlus, +kPiPlus, +kKPlus};
  std::array<int, 4> pdgOfBachelor{+kPiMinus, +kPiMinus, +kKMinus, +kKMinus};

  // Table aliases
  using TracksWMcIU = soa::Join<aod::TracksIU, McTrackLabels>;
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
      if (device.name.compare("hf-candidate-creator-xic0-omegac0-qa") == 0) {
        hfEvSelMc.init(device, registry);
        break;
      }
    }

    // Add histograms for QA
    if (configs.fillMcHistograms) {
      registry.add("hGenCharmBaryonPtRapidityTight", "Generated charm baryon #it{p}_{T};#if{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
      registry.add("hGenCharmBaryonPtRapidityLoose", "Generated charm baryon #it{p}_{T};#if{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{20, 0.0, 20.0}}});
    }
  }

  // Helper function to fill in table for MC Rec
  ///@brief decay Channel is the chosen reconstruction channel
  ///@brief flag
  ///@brief debug
  ///@brief origin
  ///@brief collisionMatched
  ///@brief pt
  ///@brief pdgcode
  template <int decayChannel>
  void fillRecoMcTableByDecayChannel(int8_t flag, int8_t debug, int8_t origin, bool collisionMatched, float pt, int pdgCode)
  {
    if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
      cursors.rowMcMatchRecXicToXiPi(flag, debug, origin, collisionMatched, pt, pdgCode);
    } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
      cursors.rowMcMatchRecOmegacToXiPi(flag, debug, origin, collisionMatched, pt, pdgCode);
    } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
      cursors.rowMcMatchRecToOmegaPi(flag, debug, origin, collisionMatched, pt, pdgCode);
    } else {
      cursors.rowMcMatchRecToOmegaKa(flag, debug, origin, collisionMatched, pt, pdgCode);
    }
  }

  // Helper function to fill in table for MC Gen
  ///@brief decay Channel is the chosen reconstruction channel
  ///@brief flag
  ///@brief debug
  ///@brief origin
  ///@brief collisionMatched
  ///@brief pt
  ///@brief pdgcode
  template <int decayChannel>
  void fillGenMcTableByDecayChannel(int8_t flag, int8_t debugGenCharmBaryon, int8_t debugGenCascade, int8_t debugGenLambda, float ptCharmBaryon, float yCharmBaryon, int8_t origin, int idxMother)
  {
    if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi) {
      cursors.rowMcMatchGenXicToXiPi(flag, debugGenCharmBaryon, debugGenCascade, debugGenLambda, ptCharmBaryon, yCharmBaryon, origin, idxMother);
    } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi) {
      cursors.rowMcMatchGenOmegacToXiPi(flag, debugGenCharmBaryon, debugGenCascade, debugGenLambda, ptCharmBaryon, yCharmBaryon, origin, idxMother);
    } else if constexpr (decayChannel == aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
      cursors.rowMcMatchGenToOmegaPi(flag, debugGenCharmBaryon, debugGenCascade, debugGenLambda, ptCharmBaryon, yCharmBaryon, origin, idxMother);
    } else {
      cursors.rowMcMatchGenToOmegaKa(flag, debugGenCharmBaryon, debugGenCascade, debugGenLambda, ptCharmBaryon, yCharmBaryon, origin, idxMother);
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
    int indexRec{-1};
    int indexRecCharmBaryon{-1};
    int8_t sign{-9};
    int8_t signCasc{-9};
    int8_t signV0{-9};
    int8_t flag{0};
    int8_t origin{0};
    int8_t debug{0};
    int8_t debugGenCharmBaryon{0};
    int8_t debugGenCasc{0};
    int8_t debugGenV0{0};
    bool collisionMatched = false;
    float ptCharmBaryonGen = -999.;
    float yCharmBaryonGen = -999.;

    ////////////////////////////////////
    // Match reconstructed candidates //
    ////////////////////////////////////

    for (const auto& candidate : candidates) {

      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = McMatchFlag::None;
      collisionMatched = false;
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
          // fill the tables
          fillRecoMcTableByDecayChannel<decayChannel>(flag, debug, origin, collisionMatched, -1.f, 0);
          continue;
        }
      }

      // CharmBaryon -> Charm bachelor + Cascade
      indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughters, pdgOfCharmBaryon[decayChannel], std::array{pdgOfCharmBachelor[decayChannel], pdgOfBachelor[decayChannel], +kProton, +kPiMinus}, true, &sign, 3);
      indexRecCharmBaryon = indexRec;
      if (indexRec == -1) { // Xic0 not reconstructed
        debug = McMatchFlag::CharmBaryonUnmatched;
      }
      if (indexRec > -1) {
        // Cascade -> Bachelor + V0
        indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersCasc, pdgOfCascade[decayChannel], std::array{pdgOfBachelor[decayChannel], +kProton, +kPiMinus}, true, &signCasc, 2);
        if (indexRec == -1) { // Xi- not reconstructed
          debug = McMatchFlag::CascUnmatched;
        }
        if (indexRec > -1) {
          // V0 -> Pos + Neg
          indexRec = RecoDecay::getMatchedMCRec<false, true>(mcParticles, arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1);
          if (indexRec == -1) { // V0 not reconstructed
            debug = McMatchFlag::V0Unmatched;
          }
          if (indexRec > -1) {
            flag = sign * (1 << decayChannel);
            collisionMatched = candidate.template collision_as<Colls>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
          }
        }
      }

      // Check if Xic0 is from b-hadron decay(prompt vs non-prompt)
      if (flag != 0) {
        auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
        origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
      }
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
        fillRecoMcTableByDecayChannel<decayChannel>(flag, debug, origin, collisionMatched, bHadMother.pt(), bHadMother.pdgCode());
      } else {
        fillRecoMcTableByDecayChannel<decayChannel>(flag, debug, origin, collisionMatched, -1.f, 0);
      }
      if (debug == McMatchFlag::CascUnmatched || debug == McMatchFlag::V0Unmatched) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fullfilled");
      }
    } // candidate loop

    ///////////////////////////////
    // Match generated particles //
    ///////////////////////////////

    for (const auto& mcCollision : mcCollisions) {

      const auto& mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());

      float centrality{-1.f};
      uint16_t rejectionMask{0};
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
          fillGenMcTableByDecayChannel<decayChannel>(0, 0, 0, 0, -999., -999., RecoDecay::OriginType::None, -1);
        }
        continue;
      }

      // Match generated particles
      for (auto const& particle : mcParticlesPerMcColl) {
        ptCharmBaryonGen = -999.;
        yCharmBaryonGen = -999.;
        flag = 0;
        sign = 0;
        debugGenCharmBaryon = 0;
        debugGenCasc = 0;
        debugGenV0 = 0;
        origin = RecoDecay::OriginType::None;
        std::vector<int> idxBhadMothers{};
        float kYCutTight = 0.5;
        float kYCutLoose = 0.8;

        // Reject particles from background events
        if (particle.fromBackgroundEvent() && configs.rejectBackground) {
          fillGenMcTableByDecayChannel<decayChannel>(flag, debugGenCharmBaryon, debugGenCasc, debugGenV0, ptCharmBaryonGen, yCharmBaryonGen, origin, -1);
          continue;
        }

        // Charm Baryon -> Cascade + Charm bachelor
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, pdgOfCharmBaryon[decayChannel], std::array{pdgOfCascade[decayChannel], pdgOfCharmBachelor[decayChannel]}, true, &sign)) {
          debugGenCharmBaryon = 1;
          ptCharmBaryonGen = particle.pt();
          yCharmBaryonGen = particle.y();
          debug = 1; // -> Matched Xic0

          for (auto const& daughterCharm : particle.template daughters_as<aod::McParticles>()) {
            if (std::abs(daughterCharm.pdgCode()) != pdgOfCascade[decayChannel]) {
              continue;
            }
            // Xi -> Lambda + pi
            if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCharm, pdgOfCascade[decayChannel], std::array{+kLambda0, pdgOfBachelor[decayChannel]}, true)) {
              debugGenCasc = 1; // -> Matched Xi-
              for (auto const& daughterCascade : daughterCharm.template daughters_as<aod::McParticles>()) {
                if (std::abs(daughterCascade.pdgCode()) != +kLambda0) {
                  continue;
                }

                // Lambda -> p + pi
                if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, daughterCascade, +kLambda0, std::array{+kProton, +kPiMinus}, true)) {
                  debugGenV0 = 1; // -> Matched Lambda0
                  flag = sign * (1 << decayChannel);
                }
              } // V0 daughter loop
            } // cascade daughter loop
          }
        } // charm daughter loop

        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
          if (std::abs(yCharmBaryonGen) < kYCutTight) {
            registry.fill(HIST("hGenCharmBaryonPtRapidityTight"), ptCharmBaryonGen);
          }
          if (std::abs(yCharmBaryonGen) < kYCutLoose) {
            registry.fill(HIST("hGenCharmBaryonPtRapidityLoose"), ptCharmBaryonGen);
          }
        }

        // Check if charm is prompt or non-prompt
        if (origin == RecoDecay::OriginType::NonPrompt) {
          fillGenMcTableByDecayChannel<decayChannel>(flag, debugGenCharmBaryon, debugGenCasc, debugGenV0, ptCharmBaryonGen, yCharmBaryonGen, origin, idxBhadMothers[0]);
        } else {
          fillGenMcTableByDecayChannel<decayChannel>(flag, debugGenCharmBaryon, debugGenCasc, debugGenV0, ptCharmBaryonGen, yCharmBaryonGen, origin, -1);
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
                                             aod::McCollisions const& mcCollisions,
                                             McCollisionsNoCents const& collsWithMcLabels,
                                             BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. No cents", false);

  void processMcXicToXiPiWithDCAFitterCentFT0C(aod::HfCandToXiPi const& candidates,
                                               TracksWMcIU const& tracks,
                                               aod::McParticles const& mcParticles,
                                               aod::McCollisions const& mcCollisions,
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
                                                aod::McCollisions const& mcCollisions,
                                                McCollisionsNoCents const& collsWithMcLabels,
                                                BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToXiPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Xi Pi. No cents", false);

  void processMcOmegacToXiPiWithDCAFitterCentFT0C(aod::HfCandToXiPi const& candidates,
                                                  TracksWMcIU const& tracks,
                                                  aod::McParticles const& mcParticles,
                                                  aod::McCollisions const& mcCollisions,
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
                                                   aod::McCollisions const& mcCollisions,
                                                   McCollisionsNoCents const& collsWithMcLabels,
                                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Pi. No cents", false);

  void processMcOmegacToOmegaPiWithDCAFitterCentFT0C(aod::HfCandToOmegaPi const& candidates,
                                                     TracksWMcIU const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     aod::McCollisions const& mcCollisions,
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
                                                   aod::McCollisions const& mcCollisions,
                                                   McCollisionsNoCents const& collsWithMcLabels,
                                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaKaWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Ka. No cents", false);

  void processMcOmegacToOmegaKaWithDCAFitterCentFT0C(aod::HfCandToOmegaK const& candidates,
                                                     TracksWMcIU const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     aod::McCollisions const& mcCollisions,
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
                                              aod::McCollisions const& mcCollisions,
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
