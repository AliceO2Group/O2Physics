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

    // Options for V0 building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // ---------------------------------------------------------------------
    Configurable<int> minCrossedRowsForV0Building{"minCrossedRowsForV0Building", 50, "minimun TPC crossed rows for daughter tracks. Used for internal V0 Building"};
    Configurable<float> dcanegtopvForV0Building{"dcanegtopvForV0Building", .1, "DCV Neg to PV"};
    Configurable<float> dcapostopvForV0Building{"dcapostopvForV0Building", .1, "DCV Pos To PV"};
    Configurable<double> v0cospaForV0Building{"v0cospaForV0Building", 0.95, "V0 CosPA"};
    Configurable<float> dcav0dauForV0Building{"dcav0dauForV0Building", 1.0, "DCA V0 Daughters"};
    Configurable<float> v0radiusForV0Building{"v0radiusForV0Building", 0.9, "v0radius"};
    Configurable<float> maxDaughterEtaForV0Building{"maxDaughterEtaForV0Building", 5.0, "Maximun daughter eta (in abs value)"};

    // Options for internal cascade building
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // --------------------------------------------------------------------
    Configurable<int> minCrossedRowsForCascadeBuilding{"minCrossedRowsForCascadeBuilding", 50, "minimun TPC crossed rows for daughter tracks. Used for internal Cascade Building"};
    Configurable<float> dcabachtopvForCascadeBuilding{"dcabachtopvForCascadeBuilding", .1, "DCV Bach to PV"};
    Configurable<float> cascradiusForCascadeBuilding{"cascradiusForCascadeBuilding", .1, "DCV Bach to PV"};
    Configurable<float> casccospaForCascadeBuilding{"casccospaForCascadeBuilding", 0.95, "Cascade CosPA"};
    Configurable<float> dcacascdauForCascadeBuilding{"dcacascdauForCascadeBuilding", 1.0, "DCA cascade daughters"};
    Configurable<float> lambdaMassWindowForCascadeBuilding{"lambdaMassWindowForCascadeBuilding", 0.10, "Distance from Lambda mass(does not apply to KF path)"};
    Configurable<float> maxDaughterEtaForCascadeBuilding{"maxDaughterEtaForCascadeBuilding", 5.0, "Maximun daughter eta (in abs value)"};

    // Options for internal cascade building - KF Building specifics
    // ...Initial values taken from PWGLF/Utiles/strangenessBuilderModule.h
    // --------------------------------------------------------------------
    Configurable<bool> kfTuneForOmega{"kfTuneForOmega", false, "if enabled, take main cascade properties from omega fit instread of Xi fit(=default)"};
    Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "2 : Daughter particle masses stay fixed in construction process"};
    Configurable<bool> kfUseV0MassConstraint{"kfUseV0MassConstraint", true, "KF : Use Lambda mass constraint"};
    Configurable<bool> kfUseCascadeMassConstraint{"kfUseCascadeMassConstraint", false, "KF : Use Cascade mass constraint - WARNING : Not adequate for inv mass analysis of Xi"};
    Configurable<bool> kfDoDCAFitterPreMinimV0{"kfDoDCAFitterPreMinimV0", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for V0"};
    Configurable<bool> kfDoDCAFitterPreMinimCasc{"kfDoDCAFitterPreMinimCasc", true, "KF : do DCAFitter pre-optimization before KF fit to include material correction for Xi"};

    // Cascade pre selection
    // --------------------
    Configurable<bool> doCascadePreselection{"doCascadePreselection", false, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
    Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascades"};
    Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in XY plane"};
    Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 1.0, "Max DCA of V0 daughter"};
    Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 1.0, "Max DCA of cascade daughter"};

    // Options for DCAFitter
    // ---------------------
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "Create tracks version propagated to PCA"};
    Configurable<double> maxR{"maxR", 200., "Reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "Reject (if>0) PCA candidate if tracks DZ exceeds this threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "Stop iteration if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "Stop iteration if Chi2/Chi2old > this"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariance, effective only if useAbsDCA is true"};

    // Options for KFParticle
    // ----------------------
    Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks comming from different collisions(effective only for KFParticle w/o derived data)"};
    Configurable<bool> kfDoCascadePreselection{"kfDoCascadePreselection", false, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
    // Configurable<bool> kfConstrainTopoV0ToCasc{"kfConstrainTopoV0ToCasc", false, "KF : Use Lambda topo constraint"}; // -> Not sure if this will be used...
    Configurable<bool> kfConstrainInvMassV0{"kfConstrainInvMassV0", true, "use mass constraint for cascade"};
    Configurable<bool> kfConstrainInvMassCasc{"kfConstrainInvMassCasc", true, "use mass constraint for cascade"};
    // Configurable<bool> kfConstrainTopoCascToCharmBaryon{"kfConstrainTopoCascToCharmBaryon", true, "use topo constraint for cascade"};
    Configurable<bool> kfConstrainInvMassCharmBaryon{"kfConstrainInvMassCharmBaryon", false, "constrain invariant mass of charm baryon candidate"};
    // Configurable<bool> constrainTopoCascToChramBaryon{"constrainTopoCascToChramBaryon", false, "constrain Casc to Charm Baryon"};
    // Configurable<bool> kfConstrainTopoCharmBaryon{"kfConstrainTopoCharmBaryon", false, "constrain charm baryon candidate to decay vertex"};

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
  // float massCharmBaryonCand{0.};

  enum CharmBaryonCandCounter { All = 0,
                                HfFlagPass,
                                CascReconstructed,
                                VertexFit };

  // Table aliases
  using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using TracksWCovDcaExtraPidPrPiKa = soa::Join<aod::TracksWCovDcaExtra, aod::TracksPidPr, aod::TracksPidPi, aod::TracksPidKa>;

  HistogramRegistry registry{"hists"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  HfEventSelection hfEvSel;

  // For cascade building
  o2::pwglf::strangenessBuilderHelper straHelper;

  // For candidate reconstruction in different PID hypothesis
  // Each decay channel assgined by following hf_cand_casc_lf_DecayType2Prong enum
  // Index 0: Xic0/Omegac0 -> Xi- Pi+ Index 1: Omegac0 -> Omega- Pi+ Index 2: Omegac0 -> Omega- K+
  // Xi- -> Lambda0 pi-, Omeag- -> Lambda0 K-
  // Lambda0 -> Pr+ Pi-
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfCharmBaryon = {+kXiC0, +kOmegaC0, +kOmegaC0}; // -> This need to be fixed...? +kXiC0, +kOmegaC0, +kOmegaC0, +kOmeagC0
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfCharmBach = {+kPiPlus, +kPiPlus, +kKPlus};
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfCascade = {+kXiMinus, +kOmegaMinus, +kOmegaMinus};
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfBach = {+kPiMinus, +kKMinus, +kKMinus};
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfV0 = {+kLambda0, +kLambda0, +kLambda0};
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfV0DauPos = {+kProton, +kProton, +kProton};
  std::array<int, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> pdgOfV0DauNeg = {+kPiMinus, +kPiMinus, +kPiMinus};

  std::array<o2::track::PID, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> trackPidOfCascade = {o2::track::PID::XiMinus, o2::track::PID::OmegaMinus, o2::track::PID::OmegaMinus};
  std::array<float, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> massOfCascades = {MassXiMinus, MassOmegaMinus, MassOmegaMinus};
  std::array<float, hf_cand_casc_lf::DecayType2Prong::N2ProngDecays> massOfCharmBach = {MassPiPlus, MassPiPlus, MassKPlus};

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
    straHelper.v0selections.minCrossedRows = configs.minCrossedRowsForV0Building;
    straHelper.v0selections.dcanegtopv = configs.dcanegtopvForV0Building;
    straHelper.v0selections.dcapostopv = configs.dcapostopvForV0Building;
    straHelper.v0selections.v0cospa = configs.v0cospaForV0Building;
    straHelper.v0selections.dcav0dau = configs.dcav0dauForV0Building;
    straHelper.v0selections.v0radius = configs.v0radiusForV0Building;
    straHelper.v0selections.maxDaughterEta = configs.maxDaughterEtaForV0Building;

    // Settings for internal Cascade building
    straHelper.cascadeselections.minCrossedRows = configs.minCrossedRowsForCascadeBuilding;
    straHelper.cascadeselections.dcabachtopv = configs.dcabachtopvForCascadeBuilding;
    straHelper.cascadeselections.cascradius = configs.cascradiusForCascadeBuilding;
    straHelper.cascadeselections.casccospa = configs.casccospaForCascadeBuilding;
    straHelper.cascadeselections.dcacascdau = configs.dcacascdauForCascadeBuilding;
    straHelper.cascadeselections.lambdaMassWindow = configs.lambdaMassWindowForCascadeBuilding;
    straHelper.cascadeselections.maxDaughterEta = configs.maxDaughterEtaForCascadeBuilding;

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
      df.setMinParamChange(configs.minParamChange);
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
    hfEvSel.init(registry, zorroSummary);

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
                               TracksWCovDcaExtraPidPrPiKa const&,
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
      auto posTrack = v0AodElement.posTrack_as<TracksWCovDcaExtraPidPrPiKa>();
      auto negTrack = v0AodElement.negTrack_as<TracksWCovDcaExtraPidPrPiKa>();
      auto bachTrack = cascAodElement.bachelor_as<TracksWCovDcaExtraPidPrPiKa>();

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
        LOG(info) << "!This cascade cannot be rebuilt(cascade ID : " << cand.cascadeId() << "/ collision ID : " << collision.globalIndex();
        continue;
      } else {
        float storeMass = (decayChannel == 0) ? straHelper.cascade.massXi : straHelper.cascade.massOmega;
        float storePt = RecoDecay::pt(straHelper.cascade.cascadeMomentum);
        registry.fill(HIST("hCascMass"), storeMass);
        registry.fill(HIST("hCascPt"), storePt);
      }

      //------------------Cascade pre selection------------------
      float massCasc = (decayChannel == 0) ? straHelper.cascade.massXi : straHelper.cascade.massOmega;

      // Perform cascade pre selection
      if (configs.doCascadePreselection) {

        // pre selection of dcaXY
        // FIXME does cascadeDCAxy represents dcaXYtoPV?
        if (std::abs(straHelper.cascade.cascadeDCAxy) > configs.dcaXYToPVCascadeMax) {
          continue;
        }

        // pre selection on invariant mass
        if (std::abs(massCasc - massOfCascades[decayChannel]) > configs.massToleranceCascade) {
          continue;
        }
      }

      if (configs.fillHistograms) {
        hCandCounter->Fill(CascReconstructed);
      }

      //------------------------------Info of V0 and cascade------------------------------
      // -> This quantities are used for physical properties of selected candidates
      // -> Not used for candidate creation
      std::array<float, 3> vertexV0 = {straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2]};
      std::array<float, 3> pVecV0 = {straHelper.cascade.v0Momentum[0], straHelper.cascade.v0Momentum[1], straHelper.cascade.v0Momentum[2]};
      std::array<float, 3> pVecV0DauPos = {straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2]};
      std::array<float, 3> pVecV0DauNeg = {straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2]};

      int chargeCasc = straHelper.cascade.charge > 0 ? 1 : -1;
      std::array<float, 3> vertexCasc = {straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2]};
      std::array<float, 3> pVecCasc = {straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2]};
      std::array<float, 21> covCasc = {0.};

      //------------------------------Create cascade track------------------------------

      constexpr std::size_t NElementsCovMatrix{6u};
      constexpr std::array<int, NElementsCovMatrix> MomInd = {9, 13, 14, 18, 19, 20}; // Momentum index?
      for (auto i = 0u; i < NElementsCovMatrix; i++) {
        covCasc[i] = straHelper.cascade.covariance[i];
        covCasc[MomInd[i]] = straHelper.cascade.covariance[MomInd[i]];
      }

      o2::track::TrackParCov trackCasc;
      if (chargeCasc < 0) { // Xi- or Omega-
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
      } else if (chargeCasc > 0) { // Xi+ or Omega+
        trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
      } else {
        continue;
      }

      trackCasc.setAbsCharge(1);
      trackCasc.setPID(trackPidOfCascade[decayChannel]);

      //------------------------------Fit SV & Create Charm Baryon track------------------------------

      // Perform secondary vertex fitting
      auto trackCharmBachelor = cand.prong0_as<TracksWCovDcaExtraPidPrPiKa>();
      auto trackParCovCharmBachelor = getTrackParCov(trackCharmBachelor);
      try {
        if (df.process(trackCasc, trackParCovCharmBachelor) == 0) {
          continue;
        }
      } catch (const std::runtime_error& e) {
        LOG(info) << "Run time error found : " << e.what() << ".DCAFitter cannot work with this candidate. SKIP!";
        continue;
      }

      // Propagate tracks to vertex
      df.propagateTracksToVertex();
      if (!df.isPropagateTracksToVertexDone()) {
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
      auto chi2SV = df.getChi2AtPCACandidate();
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

      // get v0 invariant mass
      float massLambda = straHelper.v0.massLambda;

      // Get Charm baryon invarian mass
      auto arrMomenta = std::array{pVecCascAsD, pVecCharmBachAsD};
      float massCharmBaryonCand = RecoDecay::m(arrMomenta, std::array{massOfCascades[decayChannel], massOfCharmBach[decayChannel]});
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
      // get uncertainty of the decay length -> Used in previous code...
      float phi, theta;
      getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, secondaryVertex, phi, theta);
      auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSV, phi, theta));
      auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixSV, phi, 0.));

      // calcuate ctau
      float ctV0 = RecoDecay::ct(pVecV0, decLenV0, MassLambda0);
      float ctCasc = RecoDecay::ct(pVecCasc, decLenCasc, massOfCascades[decayChannel]);
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
                              massLambda, massCasc, massCharmBaryonCand,
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
                                 massLambda, massCasc, massCharmBaryonCand,
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
                                 massLambda, massCasc, massCharmBaryonCand,
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
                                TracksWCovDcaExtraPidPrPiKa const&,
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
      auto posTrack = v0AodElement.posTrack_as<TracksWCovDcaExtraPidPrPiKa>();
      auto negTrack = v0AodElement.negTrack_as<TracksWCovDcaExtraPidPrPiKa>();
      auto bachTrack = cascAodElement.bachelor_as<TracksWCovDcaExtraPidPrPiKa>();

      // Make cascade starting from V0
      // If success, fill Cascade and V0 information for reconstruction
      if (!straHelper.buildCascadeCandidateWithKF(collision.globalIndex(),
                                                  collision.posX(), collision.posY(), collision.posZ(),
                                                  posTrack,
                                                  negTrack,
                                                  bachTrack,
                                                  false, // calculateBachelorBaryonVariables
                                                  configs.kfConstructMethod,
                                                  configs.kfTuneForOmega,
                                                  configs.kfUseV0MassConstraint,
                                                  configs.kfUseCascadeMassConstraint,
                                                  configs.kfDoDCAFitterPreMinimV0,
                                                  configs.kfDoDCAFitterPreMinimCasc)) {
        LOG(info) << "This cascade cannot be rebuilt";
        continue;
      } else {
        float storeMass = (decayChannel == 0) ? straHelper.cascade.massXi : straHelper.cascade.massOmega;
        float storePt = RecoDecay::pt(straHelper.cascade.cascadeMomentum);
        registry.fill(HIST("hCascMass"), storeMass);
        registry.fill(HIST("hCascPt"), storePt);
      }

      //------------------------------Cascade pre-selection------------------------------
      float massCasc = (decayChannel == 0) ? straHelper.cascade.massXi : straHelper.cascade.massOmega;

      if (configs.kfDoCascadePreselection) {

        // pre selection of dcaXY
        // FIXME does cascadeDCAxy represents dcaxytoPV?
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
        if (std::abs(massCasc - massOfCascades[decayChannel]) > configs.massToleranceCascade) {
          continue;
        }
      }

      if (configs.fillHistograms) {
        hCandCounter->Fill(CascReconstructed);
      }

      //------------------------------Info of V0 and Cascade------------------------------
      // -> This quantities are used for physical properties of selected candidates
      // -> Not used for candidate creation
      // auto chargeCasc = casc.sign() > 0 ? 1 : -1;
      // std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      // std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
      // std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
      // std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};

      std::array<float, 3> vertexV0 = {straHelper.cascade.v0Position[0], straHelper.cascade.v0Position[1], straHelper.cascade.v0Position[2]};
      std::array<float, 3> pVecV0 = {straHelper.cascade.v0Momentum[0], straHelper.cascade.v0Momentum[1], straHelper.cascade.v0Momentum[2]};
      std::array<float, 3> pVecV0DauPos = {straHelper.cascade.positiveMomentum[0], straHelper.cascade.positiveMomentum[1], straHelper.cascade.positiveMomentum[2]};
      std::array<float, 3> pVecV0DauNeg = {straHelper.cascade.negativeMomentum[0], straHelper.cascade.negativeMomentum[1], straHelper.cascade.negativeMomentum[2]};

      int chargeCasc = straHelper.cascade.charge > 0 ? 1 : -1;
      std::array<float, 3> vertexCasc = {straHelper.cascade.cascadePosition[0], straHelper.cascade.cascadePosition[1], straHelper.cascade.cascadePosition[2]};
      std::array<float, 3> pVecCasc = {straHelper.cascade.cascadeMomentum[0], straHelper.cascade.cascadeMomentum[1], straHelper.cascade.cascadeMomentum[2]};
      std::array<float, 21> covCasc = {0.};

      //------------------------------Create Charm Baryon as KF Particle object------------------------------

      // initialize primary vertex
      KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
      float covMatrixPV[6];
      kfpVertex.GetCovarianceMatrix(covMatrixPV);
      KFParticle kfPv(kfpVertex); // -> For calculation of DCAs to PV

      //~~~~~create charm bachelor track into KFParticle object~~~~~//
      auto trackCharmBachelor = cand.prong0_as<TracksWCovDcaExtraPidPrPiKa>();
      KFPTrack kfpTrackCharmBachelor = createKFPTrackFromTrack(trackCharmBachelor);
      KFParticle kfCharmBachelor(kfpTrackCharmBachelor, pdgOfCharmBach[decayChannel]);

      //~~~~~create Cascade as KFParticle object~~~~~//
      // -> This process is essential in building Charm baryon via KFParticle package.
      //    Since we don't have any information given from other tables about Charm baryon,
      //    only way to construct Charm baryon is to construct Charm baryon from it's daughters.
      constexpr std::size_t NElementsStateVectorCasc{6};
      std::array<float, NElementsStateVectorCasc> xyzpxpypzCasc = {vertexCasc[0], vertexCasc[1], vertexCasc[2], pVecCasc[0], pVecCasc[1], pVecCasc[2]};
      float parPosMomCasc[NElementsStateVectorCasc];
      std::copy(xyzpxpypzCasc.begin(), xyzpxpypzCasc.end(), parPosMomCasc);

      KFParticle kfCasc;
      // kfCasc.Create(parPosMomCasc, casc.kfTrackCovMat(), casc.sign(), massCasc);
      kfCasc.Create(parPosMomCasc, straHelper.cascade.covariance, straHelper.cascade.charge, massCasc);
      if (configs.kfConstrainInvMassCasc) {
        kfCasc.SetNonlinearMassConstraint(massOfCascades[decayChannel]);
      }

      //~~~~~~Create Charm Baryon as KFParticle object~~~~~//
      KFParticle kfCharmBaryon;
      const KFParticle* kfDaughterCharmBaryon[2] = {&kfCharmBachelor, &kfCasc};
      kfCharmBaryon.SetConstructMethod(configs.kfConstructMethod);
      try {
        kfCharmBaryon.Construct(kfDaughterCharmBaryon, 2);
      } catch (std::runtime_error& e) {
        LOG(debug) << "Failed to construct Charm Baryon : " << e.what();
        continue;
      }

      // Get covariance matrix of Charm Baryon
      auto covMatrixCharmBaryon = kfCharmBaryon.CovarianceMatrix();

      // transport Charm Baryon daughters to Charm Baryon to decay vertex
      // FIXME Is this process necessary?
      float secondaryVertex[3] = {kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()};
      kfCasc.TransportToPoint(secondaryVertex);
      kfCharmBachelor.TransportToPoint(secondaryVertex);

      // Fill histogram
      if (configs.fillHistograms) {
        hCandCounter->Fill(VertexFit);
        hInvMassCharmBaryon->Fill(kfCharmBaryon.GetMass());
      }

      // Transport cascade and charm baryon to decay vertex(just to be sure)
      kfCharmBaryon.TransportToDecayVertex();
      kfCasc.TransportToDecayVertex();

      //!~~~~~Extra calculations for QA~~~~~~!//
      /// These quantities were calculated to fill in existing output tables & QA purpose
      /// In the future, these calculation can be skipped? This must be dicussed.

      /// Extra calculation for V0 daughters
      KFPTrack kfpTrackV0Pos = createKFPTrackFromTrack(posTrack);
      KFPTrack kfpTrackV0Neg = createKFPTrackFromTrack(negTrack);

      int pdgV0Pos = (trackCharmBachelor.sign() > 0) ? kProton : kPiPlus;
      int pdgV0Neg = (trackCharmBachelor.sign() > 0) ? kPiMinus : -kProton;

      KFParticle kfV0Pos(kfpTrackV0Pos, pdgV0Pos);
      KFParticle kfV0Neg(kfpTrackV0Neg, pdgV0Neg);

      auto trackParCovV0Pos = getTrackParCovFromKFP(kfV0Pos, kfV0Pos.GetPDG(), 1);
      auto trackParCovV0Neg = getTrackParCovFromKFP(kfV0Neg, kfV0Neg.GetPDG(), -1);

      /// Extra calculations for V0
      KFParticle kfV0;
      float xyzpxpypzV0[6] = {vertexV0[0], vertexV0[1], vertexV0[2], pVecV0[0], pVecV0[1], pVecV0[2]};
      kfV0.Create(xyzpxpypzV0, straHelper.cascade.kfTrackCovarianceV0, 0, MassLambda0);
      if (configs.kfConstrainInvMassV0) {
        kfV0.SetNonlinearMassConstraint(MassLambda);
      }

      float kfMassV0, kfSigMassV0;
      kfV0.GetMass(kfMassV0, kfSigMassV0);

      KFParticle kfV0ConstrainedToPv = kfV0;
      KFParticle kfV0ConstrainedToCasc = kfV0;
      kfV0ConstrainedToPv.SetProductionVertex(kfPv);
      kfV0ConstrainedToCasc.SetProductionVertex(kfCasc);

      /// Extra calculations for bachelor
      KFPTrack kfpTrackBach = createKFPTrackFromTrack(bachTrack);
      KFParticle kfBach(kfpTrackBach, pdgOfBach[decayChannel]);

      KFParticle kfBachConstrainedToPv = kfBach;
      KFParticle kfBachConstrainedToCasc = kfBach;
      kfBachConstrainedToPv.SetProductionVertex(kfPv);
      kfBachConstrainedToCasc.SetProductionVertex(kfCasc);

      auto trackParCovBach = getTrackParCovFromKFP(kfBachConstrainedToCasc, kfBachConstrainedToCasc.GetPDG(), bachTrack.sign());

      /// Extra calculations for casc
      float kfMassCasc, kfSigMassCasc;
      kfCasc.GetMass(kfMassCasc, kfSigMassCasc);

      KFParticle kfCascConstrainedToPv = kfCasc;
      KFParticle kfCascConstrainedToCharmBaryon = kfCasc;
      kfCascConstrainedToPv.SetProductionVertex(kfPv);
      kfCascConstrainedToCharmBaryon.SetProductionVertex(kfCharmBaryon);

      /// Extra calculation for charm bachelor
      KFParticle kfCharmBachelorConstrainedToPv = kfCharmBachelor;
      KFParticle kfCharmBachelorConstrainedToCharmBaryon = kfCharmBachelor;
      kfCharmBachelorConstrainedToPv.SetProductionVertex(kfPv);
      kfCharmBachelorConstrainedToCharmBaryon.SetProductionVertex(kfCharmBaryon);

      /// Extra calculation for charm baryon
      float kfMassCharmBaryon, kfSigMassCharmBaryon;
      kfCharmBaryon.GetMass(kfMassCharmBaryon, kfSigMassCharmBaryon);

      KFParticle kfCharmBaryonConstrainedToPv = kfCharmBaryon;
      kfCharmBaryonConstrainedToPv.SetProductionVertex(kfPv);

      //------------------------------Calculate physical quantities and fill candidate table------------------------------

      // Get updated daughter tracks after vertex fit
      o2::track::TrackParCov trackParCovCasc = getTrackParCovFromKFP(kfCascConstrainedToCharmBaryon, kfCascConstrainedToCharmBaryon.GetPDG(), kfBach.GetQ());
      trackParCovCasc.setAbsCharge(1);
      o2::track::TrackParCov trackParCovCharmBachelor = getTrackParCovFromKFP(kfCharmBachelorConstrainedToCharmBaryon, kfCharmBachelorConstrainedToCharmBaryon.GetPDG(), -kfBach.GetQ());
      trackParCovCharmBachelor.setAbsCharge(1);

      // impact parameters
      std::array<float, 2> impactParameterV0DauPos;
      std::array<float, 2> impactParameterV0DauNeg;
      std::array<float, 2> impactParameterBach;
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Pos, 2.f, matCorr, &impactParameterV0DauPos);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovV0Neg, 2.f, matCorr, &impactParameterV0DauNeg);
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
      o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCharmBachelor, 2.f, matCorr, &impactParameterCharmBachelor);

      // get Chi2Topo/NDF
      float chi2NdfTopoV0ToPv = kfV0ConstrainedToPv.GetChi2() / kfV0ConstrainedToPv.GetNDF();
      float chi2NdfTopoV0ToCasc = kfV0ConstrainedToPv.GetChi2() / kfV0ConstrainedToPv.GetNDF();
      float chi2NdfTopoBachToPv = kfBachConstrainedToPv.GetChi2() / kfBachConstrainedToPv.GetNDF();
      float chi2NdfTopoBachToCasc = kfBachConstrainedToCasc.GetChi2() / kfBachConstrainedToCasc.GetNDF();
      float chi2NdfTopoCascToPv = kfCascConstrainedToPv.GetChi2() / kfCascConstrainedToPv.GetNDF();
      float chi2NdfTopoCascToCharmBaryon = kfCascConstrainedToCharmBaryon.GetChi2() / kfCascConstrainedToCharmBaryon.GetNDF();
      float chi2NdfTopoCharmBachelorToPv = kfCharmBachelorConstrainedToPv.GetChi2() / kfCharmBachelorConstrainedToPv.GetNDF();
      float chi2NdfTopoCharmBachelorToCharmBaryon = kfCharmBachelorConstrainedToCharmBaryon.GetChi2() / kfCharmBachelorConstrainedToCharmBaryon.GetNDF();
      float chi2NdfTopoCharmBaryonToPv = kfCharmBachelorConstrainedToPv.GetChi2() / kfCharmBachelorConstrainedToPv.GetNDF();

      float deviationCharmBachelorToPv = kfCalculateChi2ToPrimaryVertex(kfCharmBachelor, kfPv);
      float chi2NdfBach = kfBach.GetChi2() / kfBach.GetNDF();

      // get ldl
      float ldlV0 = ldlFromKF(kfV0, kfPv);
      float ldlCasc = ldlFromKF(kfCasc, kfPv);
      float ldlCharmBaryon = ldlFromKF(kfCharmBaryon, kfPv);

      // get DCAs
      float kfDcaV0Daughters = kfV0Pos.GetDistanceFromParticle(kfV0Neg);
      float kfDcaCascDaughters = kfBachConstrainedToCasc.GetDistanceFromParticle(kfV0ConstrainedToCasc);
      float kfDcaCharmBaryonDaughters = kfCharmBachelorConstrainedToCharmBaryon.GetDistanceFromParticle(kfCascConstrainedToCharmBaryon);
      float kfDcaXYCharmBachelorToPv = kfCharmBachelorConstrainedToCharmBaryon.GetDistanceFromVertexXY(kfPv);
      float kfDcaXYCascToPv = kfCascConstrainedToCharmBaryon.GetDistanceFromVertexXY(kfPv);

      // get decay length - In XY
      float decayLXYV0, errDecayLXYV0;
      kfV0ConstrainedToCasc.GetDecayLengthXY(decayLXYV0, errDecayLXYV0);

      float decayLXYCasc, errDecayLXYCasc;
      kfCascConstrainedToCharmBaryon.GetDecayLengthXY(decayLXYCasc, decayLXYCasc);

      float decayLXYCharmBaryon, errDecayLXYCharmBaryon;
      kfCharmBaryonConstrainedToPv.GetDecayLengthXY(decayLXYCharmBaryon, errDecayLXYCharmBaryon);

      // get decay length - In XYZ
      float decayLV0 = RecoDecay::distance(std::array<float, 3>{kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()}, std::array<float, 3>{kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()});
      float decayLCasc = RecoDecay::distance(std::array<float, 3>{kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()}, std::array<float, 3>{kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ()});
      float decayLCharmBaryon = RecoDecay::distance(std::array<float, 3>{collision.posX(), collision.posY(), collision.posZ()}, std::array<float, 3>{kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()});

      double phiCharmBaryon, thetaCharmBaryon;
      getPointDirection(std::array<float, 3>{kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()}, std::array<float, 3>{kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ()}, phiCharmBaryon, thetaCharmBaryon);
      float errDecayLCharmBaryon = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covMatrixCharmBaryon, phiCharmBaryon, thetaCharmBaryon));

      // get cosine of pointing angle
      std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};
      // float cosPaV0ToPv = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float cosPaV0ToPv = RecoDecay::cpa(pvCoord, vertexV0, pVecV0);
      float cosPaCascToPv = cpaFromKF(kfCasc, kfPv);
      float cosPaCharmBaryonToPv = cpaFromKF(kfCharmBaryon, kfPv);

      float cosPaXYV0ToPv = RecoDecay::cpaXY(pvCoord, vertexV0, pVecV0);
      float cosPaXYCascToPv = cpaXYFromKF(kfCasc, kfPv);
      float cosPaXYCharmBaryonToPv = cpaXYFromKF(kfCharmBaryon, kfPv);

      float cosPaV0ToCasc = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
      float cosPaCascToCharmBaryon = cpaFromKF(kfCasc, kfCharmBaryon);
      float cosPaXYV0ToCasc = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
      float cosPaXYCascToCharmBaryon = cpaXYFromKF(kfCasc, kfCharmBaryon);

      // KF pT, eta
      float ptCasc = kfCascConstrainedToCharmBaryon.GetPt();
      float ptCharmBachelor = kfCharmBachelorConstrainedToCharmBaryon.GetPt();
      float ptCharmBaryon = kfCharmBaryon.GetPt();
      float yCharmBaryon = kfCharmBaryon.GetRapidity();

      // get KF cosThetaStar(?)
      float cosThetaStarCharmBachelorXic0, cosThetaStarCharmBachelorOmegac0;
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        cosThetaStarCharmBachelorXic0 = cosThetaStarFromKF(0, pdgOfCharmBaryon[decayChannel], pdgOfCascade[decayChannel], pdgOfCharmBach[decayChannel], kfCascConstrainedToCharmBaryon, kfCharmBachelorConstrainedToCharmBaryon);
        cosThetaStarCharmBachelorOmegac0 = cosThetaStarFromKF(0, pdgOfCharmBaryon[decayChannel + 1], pdgOfCascade[decayChannel], pdgOfCharmBach[decayChannel], kfCascConstrainedToCharmBaryon, kfCharmBachelorConstrainedToCharmBaryon);
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        cosThetaStarCharmBachelorOmegac0 = cosThetaStarFromKF(0, pdgOfCharmBaryon[decayChannel + 1], pdgOfCascade[decayChannel], pdgOfCharmBach[decayChannel], kfCascConstrainedToCharmBaryon, kfCharmBachelorConstrainedToCharmBaryon);
      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
        cosThetaStarCharmBachelorXic0 = cosThetaStarFromKF(0, +kXiC0, pdgOfCascade[decayChannel], pdgOfCharmBach[decayChannel], kfCascConstrainedToCharmBaryon, kfCharmBachelorConstrainedToCharmBaryon);
        cosThetaStarCharmBachelorOmegac0 = cosThetaStarFromKF(0, +kOmegaC0, pdgOfCascade[decayChannel], pdgOfCharmBach[decayChannel], kfCascConstrainedToCharmBaryon, kfCharmBachelorConstrainedToCharmBaryon);
      }

      // KF ct
      float ctV0 = kfV0ConstrainedToCasc.GetLifeTime();
      float ctCasc = kfCascConstrainedToCharmBaryon.GetLifeTime();
      float ctCharmBaryon = kfCharmBachelorConstrainedToPv.GetLifeTime();

      //------------------------------Calculate physical quantities and fill candidate table------------------------------

      //------------------------------Fill the table------------------------------
      if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi) {
        cursors.rowCandToXiPiKf(collision.globalIndex(),                     // Global index of collision
                                pvCoord[0], pvCoord[1], pvCoord[2],          // coordination of PV
                                kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ(), // Decay position of kfCasc
                                vertexV0[0], vertexV0[1], vertexV0[2],
                                bachTrack.sign(),
                                covMatrixCharmBaryon[0], covMatrixCharmBaryon[1], covMatrixCharmBaryon[2], covMatrixCharmBaryon[3], covMatrixCharmBaryon[4], covMatrixCharmBaryon[5],
                                kfCharmBaryon.GetPx(), kfCharmBaryon.GetPy(), kfCharmBaryon.GetPz(), // x, y, z momentum of charm baryon
                                kfCascConstrainedToCharmBaryon.GetPx(), kfCascConstrainedToCharmBaryon.GetPy(), kfCascConstrainedToCharmBaryon.GetPz(),
                                kfCharmBachelorConstrainedToCharmBaryon.GetPx(), kfCharmBachelorConstrainedToCharmBaryon.GetPy(), kfCharmBachelorConstrainedToCharmBaryon.GetPz(),
                                pVecV0[0], pVecV0[1], pVecV0[2],
                                kfBachConstrainedToCasc.GetPx(), kfBachConstrainedToCasc.GetPy(), kfBachConstrainedToCasc.GetPz(),
                                pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                                pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                                cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                                cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                                kfMassV0, kfMassCasc, kfMassCharmBaryon,
                                cosPaV0ToPv, cosPaCascToPv,
                                // cosPaCharmBaryonToPv, cosPaXYV0ToPv, cosPaXYCharmBaryonToPv, cosPaXYCascToPv,
                                ctCasc, ctV0, ctCharmBaryon,
                                kfV0Pos.GetEta(), kfV0Neg.GetEta(), kfBach.GetEta(), kfCharmBachelorConstrainedToCharmBaryon.GetEta(),
                                kfCharmBaryon.GetEta(), kfCasc.GetEta(), kfV0.GetEta(),
                                dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                                kfDcaCascDaughters, kfDcaV0Daughters, kfDcaCharmBaryonDaughters,
                                kfDcaXYCharmBachelorToPv, kfDcaXYCascToPv,
                                kfV0.GetChi2(), kfCasc.GetChi2(), kfCharmBaryon.GetChi2(),
                                /*FIXME chi2 of mass constrained, V0 and casc*/ kfV0.GetChi2(), kfCasc.GetChi2(),
                                ldlV0, ldlCasc, // ldlCharmBaryon,
                                chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoCharmBachelorToPv, chi2NdfTopoCharmBaryonToPv,
                                chi2NdfTopoV0ToCasc, chi2NdfTopoCascToCharmBaryon,
                                decayLXYV0, decayLXYCasc, decayLXYCharmBaryon,
                                cosPaV0ToCasc, cosPaCascToCharmBaryon, // cosPaXYV0ToCasc, cosPaXYCascToCharmBaryon,
                                yCharmBaryon,                          // ptCharmBachelor, ptCharmBaryon,
                                cosThetaStarCharmBachelorXic0,
                                kfV0.GetNDF(), kfCasc.GetNDF(), kfCharmBaryon.GetNDF(), /*FIXME chi2, NDF of mass constrained V0/casc*/ kfV0.GetNDF(), kfCasc.GetNDF(),
                                kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.GetChi2() / kfCasc.GetNDF(), kfCharmBaryon.GetChi2() / kfCharmBaryon.GetNDF(), /*FIXME chi2, NDF of mass constrained V0/casc*/ kfV0.GetChi2() / kfV0.GetNDF(), kfV0.GetChi2() / kfV0.GetNDF());

      } else if constexpr (decayChannel == hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
        cursors.rowCandToOmegaPi(collision.globalIndex(),
                                 pvCoord[0], pvCoord[1], pvCoord[2],
                                 /*vertexCharmBaryonFromFitter. For KF, this is 0*/ 0.f, 0.f, 0.f,
                                 kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ(),
                                 vertexV0[0], vertexV0[1], vertexV0[2],
                                 trackCharmBachelor.sign(),
                                 covMatrixCharmBaryon[0], covMatrixCharmBaryon[1], covMatrixCharmBaryon[2], covMatrixCharmBaryon[3], covMatrixCharmBaryon[4], covMatrixCharmBaryon[5],
                                 kfCharmBaryon.GetPx(), kfCharmBaryon.GetPy(), kfCharmBaryon.GetPz(),
                                 kfCascConstrainedToCharmBaryon.GetPx(), kfCascConstrainedToCharmBaryon.GetPy(), kfCascConstrainedToCharmBaryon.GetPz(),
                                 kfCharmBachelorConstrainedToCharmBaryon.GetPx(), kfCharmBachelorConstrainedToCharmBaryon.GetPy(), kfCharmBachelorConstrainedToCharmBaryon.GetPz(),
                                 pVecV0[0], pVecV0[1], pVecV0[2],
                                 kfBachConstrainedToCasc.GetPx(), kfBachConstrainedToCasc.GetPy(), kfBachConstrainedToCasc.GetPz(),
                                 pVecV0DauPos[0], pVecV0DauPos[1], pVecV0DauPos[2],
                                 pVecV0DauNeg[0], pVecV0DauNeg[1], pVecV0DauNeg[2],
                                 impactParameterCasc.getY(), impactParameterCharmBachelor.getY(),
                                 impactParameterCasc.getZ(), impactParameterCharmBachelor.getZ(),
                                 std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterCharmBachelor.getSigmaY2()),
                                 cascAodElement.v0Id(), v0AodElement.posTrackId(), v0AodElement.negTrackId(),
                                 cand.cascadeId(), trackCharmBachelor.globalIndex(), cand.prong0Id(),
                                 kfMassV0, kfMassCasc, kfMassCharmBaryon,
                                 cosPaV0ToPv, cosPaCharmBaryonToPv, cosPaCascToPv, cosPaXYV0ToPv, cosPaXYCharmBaryonToPv, cosPaXYCascToPv,
                                 ctCharmBaryon, ctCasc, ctV0,
                                 kfV0Pos.GetEta(), kfV0Neg.GetEta(), kfBach.GetEta(), kfCharmBachelorConstrainedToCharmBaryon.GetEta(),
                                 kfCharmBaryon.GetEta(), kfCasc.GetEta(), kfV0.GetEta(),
                                 dcaxyV0DauPos, dcaxyV0DauNeg, dcaxyBach,
                                 dcazV0DauPos, dcazV0DauNeg, dcazBach,
                                 kfDcaCascDaughters, straHelper.cascade.v0DaughterDCA, kfDcaCharmBaryonDaughters,
                                 decayLCharmBaryon, decayLCasc, decayLV0, errDecayLCharmBaryon, errDecayLXYCharmBaryon, cand.hfflag());

        cursors.rowCandToOmegaPiKf(kfDcaXYCharmBachelorToPv, kfDcaXYCascToPv,
                                   /*V0 chi2. Taken from LF*/ kfV0.GetChi2(), kfCasc.GetChi2(), kfCharmBaryon.GetChi2(), /*Mass constraint only done when requested*/ kfV0.GetChi2(), /*Mass constraint only done when requested*/ kfCasc.GetChi2(),
                                   ldlV0, ldlCasc, ldlCharmBaryon,
                                   chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoCharmBachelorToPv, chi2NdfTopoCharmBaryonToPv, deviationCharmBachelorToPv,
                                   chi2NdfTopoV0ToCasc, chi2NdfTopoCascToCharmBaryon,
                                   decayLXYV0, decayLXYCasc, decayLXYCharmBaryon,
                                   cosPaV0ToCasc, cosPaCascToCharmBaryon, cosPaXYV0ToCasc, cosPaXYCascToCharmBaryon,
                                   kfCharmBaryon.GetRapidity(), ptCharmBachelor, ptCharmBaryon,
                                   cosThetaStarCharmBachelorOmegac0,
                                   kfV0.GetNDF(), kfCasc.GetNDF(), kfCharmBaryon.GetNDF(), /*Mass constraint only done when requested*/ kfV0.GetNDF(), /*Mass constraint only done when requested*/ kfCasc.GetNDF(),
                                   kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.Chi2() / kfCasc.GetNDF(), kfCharmBaryon.GetChi2() / kfCharmBaryon.GetNDF(),
                                   /*Mass constraint only done when requested*/ kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.Chi2() / kfCasc.GetNDF(),
                                   /*FIXME casc-rej not calculated. For now, fill in mass of KFCasc*/ kfCasc.GetMass());

      } else {
        cursors.rowCandToOmegaKaKf(collision.globalIndex(),
                                   collision.posX(), collision.posY(), collision.posZ(),
                                   kfPv.GetX(), kfPv.GetY(), kfPv.GetZ(),
                                   vertexV0[0], vertexV0[1], vertexV0[2],
                                   pVecV0[0], pVecV0[1], pVecV0[2],
                                   vertexCasc[0], vertexCasc[1], vertexCasc[2],
                                   pVecCasc[0], pVecCasc[1], pVecCasc[2],
                                   vertexV0[0], vertexV0[1], vertexV0[2],
                                   pVecV0[0], pVecV0[1], pVecV0[2],
                                   kfCasc.GetX(), kfCasc.GetY(), kfCasc.GetZ(),
                                   kfCasc.GetPx(), kfCasc.GetPy(), kfCasc.GetPz(),
                                   kfCharmBaryon.GetX(), kfCharmBaryon.GetY(), kfCharmBaryon.GetZ(),
                                   kfCharmBaryon.GetPx(), kfCharmBaryon.GetPy(), kfCharmBaryon.GetPz(),
                                   straHelper.cascade.charge,
                                   kfV0Pos.GetEta(), kfV0Neg.GetEta(), kfBach.GetEta(), kfCharmBachelor.GetEta(), kfV0.GetEta(), kfCasc.GetEta(), kfCharmBaryon.GetEta(), kfCharmBaryon.GetRapidity(),
                                   impactParameterCharmBachelor.getY(), std::sqrt(impactParameterCharmBachelor.getSigmaY2()), impactParameterCasc.getY(), std::sqrt(impactParameterCasc.getSigmaY2()),
                                   kfDcaV0Daughters, kfDcaCascDaughters, kfDcaCharmBaryonDaughters,
                                   cosPaV0ToPv, cosPaCascToPv, cosPaCharmBaryonToPv, cosPaXYV0ToPv, cosPaXYCascToPv, cosPaXYCharmBaryonToPv, cosPaV0ToCasc, cosPaCascToCharmBaryon, cosPaXYV0ToCasc, cosPaXYCascToCharmBaryon,
                                   kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.GetChi2() / kfCasc.GetNDF(), kfCharmBaryon.GetChi2() / kfCharmBaryon.GetNDF(),
                                   /*FIXME mass constraint is optional*/ kfV0.GetChi2() / kfV0.GetNDF(), kfCasc.GetChi2() / kfCasc.GetNDF(),
                                   chi2NdfTopoV0ToCasc, chi2NdfTopoBachToCasc, chi2NdfTopoCharmBachelorToCharmBaryon, chi2NdfTopoCascToCharmBaryon, // Topological constraints to mother
                                   chi2NdfTopoV0ToPv, chi2NdfTopoCascToPv, chi2NdfTopoCharmBachelorToPv, chi2NdfTopoCharmBaryonToPv,                // Topological constraints to PV
                                   ldlV0, ldlCasc, ldlCharmBaryon,
                                   decayLXYV0, decayLXYCasc, decayLXYCharmBaryon,
                                   kfMassV0, kfSigMassV0, kfMassCasc, kfSigMassCasc, /*FIXME no -rej has been made*/ kfMassCasc, /*FIXME no -rej has been made*/ kfSigMassCasc, kfMassCharmBaryon, kfSigMassCharmBaryon,
                                   ptCharmBaryon, ptCharmBachelor, ptCasc,
                                   cosThetaStarCharmBachelorOmegac0, cosThetaStarCharmBachelorXic0, ctV0, ctCasc, ctCharmBaryon,
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
                                        TracksWCovDcaExtraPidPrPiKa const& tracks,
                                        aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
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
                                          TracksWCovDcaExtraPidPrPiKa const& tracks,
                                          aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithDCAFitterCentFT0C, "Charm candidate reconstruction with Xi Pi via DcaFitter method, centrality selection on FT0C", false);

  void processToXiPiWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                          aod::HfCascLf2Prongs const& candidates,
                                          aod::Cascades const& cascades,
                                          aod::V0s const& v0s,
                                          TracksWCovDcaExtraPidPrPiKa const& tracks,
                                          aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithDCAFitterCentFT0M, "Charm candidate reconstruction with Xi Pi via DcaFitter method, centrality selection on FT0M", false);

  /*~~~~~~~~~~~~~~~~~*/
  /*~~~To Omega Pi~~~*/
  /*~~~~~~~~~~~~~~~~~*/
  void processToOmegaPiWithDCAFitterNoCent(SelectedCollisions const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           TracksWCovDcaExtraPidPrPiKa const& tracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithDCAFitterNoCent, "Charm candidte reconstruction with Omega Pi via DcaFitter method, no centrality", false);

  void processToOmegaPiWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             TracksWCovDcaExtraPidPrPiKa const& tracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithDCAFitterCentFT0C, "Charm candidate reconstruction with Omega Pi via DcaFitter method, centrality selection on FT0C", false);

  void processToOmegaPiWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             TracksWCovDcaExtraPidPrPiKa const& tracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hInvMassCharmBaryonToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithDCAFitterCentFT0M, "Charm candidate reconstruction with Omega Pi via DcaFitter method, centrality selection on FT0M", false);

  /*~~~~~~~~~~~~~~~~~*/
  /*~~~To Omega Ka~~~*/
  /*~~~~~~~~~~~~~~~~~*/
  void processToOmegaKaWithDCAFitterNoCent(SelectedCollisions const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           TracksWCovDcaExtraPidPrPiKa const& tracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithDCAFitterNoCent, "Charm candidte reconstruction with Omega Ka via DcaFitter method, no centrality", false);

  void processToOmegaKaWithDCAFitterCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             TracksWCovDcaExtraPidPrPiKa const& tracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithDCAFitterCentFT0C, "Charm candidate reconstruction with Omega Ka via DcaFitter method, centrality selection on FT0C", false);

  void processToOmegaKaWithDCAFitterCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                             aod::HfCascLf2Prongs const& candidates,
                                             aod::Cascades const& cascades,
                                             aod::V0s const& v0s,
                                             TracksWCovDcaExtraPidPrPiKa const& tracks,
                                             aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithDCAFitter<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaPi);
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
                                         aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithKFParticleNoCent, "Charm Baryon decaying to Xi Pi reconstruction via KFParticle method, no centrality", false);

  void processToXiPiWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           TracksWCovDcaExtraPidPrPiKa const& tracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToXiPiWithKFParticleCentFT0C, "Charm Baryon decaying to Xi Pi reconstruction via KFParticle method, centrality on FT0C", false);

  void processToXiPiWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                           aod::HfCascLf2Prongs const& candidates,
                                           aod::Cascades const& cascades,
                                           aod::V0s const& v0s,
                                           TracksWCovDcaExtraPidPrPiKa const& tracks,
                                           aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::XiczeroOmegaczeroToXiPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToXiPi, hCandidateCounterToXiPi);
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
                                            aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithKFParticleNoCent, "Charm Baryon decaying to Omega Pi reconstruction via KFParticle method, no centrality", false);

  void processToOmegaPiWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaPiWithKFParticleCentFT0C, "Charm Baryon decaying to Omega Pi reconstruction via KFParticle method, centrality on FT0C", false);

  void processToOmegaPiWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaPi, hCandidateCounterToOmegaPi);
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
                                            aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::None, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithKFParticleNoCent, "Charm Baryon decaying to Omega Ka reconstruction via KFParticle method, no centrality", false);

  void processToOmegaKaWithKFParticleCentFT0C(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0C, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0Qa, processToOmegaKaWithKFParticleCentFT0C, "Charm Baryon decaying to Omega Ka reconstruction via KFParticle method, centrality on FT0C", false);

  void processToOmegaKaWithKFParticleCentFT0M(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions,
                                              aod::HfCascLf2Prongs const& candidates,
                                              aod::Cascades const& cascades,
                                              aod::V0s const& v0s,
                                              TracksWCovDcaExtraPidPrPiKa const& tracks,
                                              aod::BCsWithTimestamps const& bcsWithTimestamps)
  {
    runCreatorWithKfParticle<CentralityEstimator::FT0M, hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK>(collisions, candidates, cascades, v0s, tracks, bcsWithTimestamps, hInvMassCharmBaryonToOmegaKa, hCandidateCounterToOmegaKa);
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
                        TracksWMc const& tracks,
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

      auto arrayDaughters = std::array{candidate.template bachelorFromCharmBaryon_as<aod::TracksWMc>(),
                                       candidate.template bachelor_as<aod::TracksWMc>(),
                                       candidate.template posTrack_as<aod::TracksWMc>(),
                                       candidate.template negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersCasc = std::array{candidate.template bachelor_as<aod::TracksWMc>(),
                                           candidate.template posTrack_as<aod::TracksWMc>(),
                                           candidate.template negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.template posTrack_as<aod::TracksWMc>(),
                                         candidate.template negTrack_as<aod::TracksWMc>()};

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
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollision, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0C) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      } else if constexpr (centEstimator == o2::hf_centrality::CentralityEstimator::FT0M) {
        const auto collSlice = collsWithMcLabels.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
        rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, centEstimator>(mcCollision, collSlice, centrality);
      }

      hfEvSelMc.fillHistograms<centEstimator>(mcCollision, rejectionMask);

      if (rejectionMask != 0) { // none of the event selection was satisfied(?) -> Reject all particles from this event
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

        // Check if charm is prompt or non-prompt
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
        }
        if (std::abs(yCharmBaryonGen) < kYCutTight) {
          // Fill in some QA histogram. Will be implemented later
        }
        if (std::abs(yCharmBaryonGen) < kYCutLoose) {
          // Fill in some QA histograms. Will be implemented later
        }

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
                                             aod::TracksWMc const& tracks,
                                             aod::McParticles const& mcParticles,
                                             aod::McCollisions const& mcCollisions,
                                             McCollisionsNoCents const& collsWithMcLabels,
                                             BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. No cents", false);

  void processMcXicToXiPiWithDCAFitterCentFT0C(aod::HfCandToXiPi const& candidates,
                                               aod::TracksWMc const& tracks,
                                               aod::McParticles const& mcParticles,
                                               aod::McCollisions const& mcCollisions,
                                               McCollisionsFT0Cs const& collsWithMcLabels,
                                               BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. Cents with FT0C", false);

  void processMcXicToXiPiWithDCAFitterCentFT0M(aod::HfCandToXiPi const& candidates,
                                               aod::TracksWMc const& tracks,
                                               aod::McParticles const& mcParticles,
                                               McCollisionsCentFT0Ms const& mcCollisions,
                                               McCollisionsFT0Ms const& collsWithMcLabels,
                                               BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0M, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithDCAFitterCentFT0M, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. Cents with FT0M", false);

  void processMcOmegacToXiPiWithDCAFitterNoCent(aod::HfCandToXiPi const& candidates,
                                                aod::TracksWMc const& tracks,
                                                aod::McParticles const& mcParticles,
                                                aod::McCollisions const& mcCollisions,
                                                McCollisionsNoCents const& collsWithMcLabels,
                                                BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToXiPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Xi Pi. No cents", false);

  void processMcOmegacToXiPiWithDCAFitterCentFT0C(aod::HfCandToXiPi const& candidates,
                                                  aod::TracksWMc const& tracks,
                                                  aod::McParticles const& mcParticles,
                                                  aod::McCollisions const& mcCollisions,
                                                  McCollisionsFT0Cs const& collsWithMcLabels,
                                                  BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToXiPiWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Omeagc0 to Xi Pi. Cents with FT0C", false);

  void processMcOmegacToXiPiWithDCAFitterCentFT0M(aod::HfCandToXiPi const& candidates,
                                                  aod::TracksWMc const& tracks,
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
                                                   aod::TracksWMc const& tracks,
                                                   aod::McParticles const& mcParticles,
                                                   aod::McCollisions const& mcCollisions,
                                                   McCollisionsNoCents const& collsWithMcLabels,
                                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaPiWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Pi. No cents", false);

  void processMcOmegacToOmegaPiWithDCAFitterCentFT0C(aod::HfCandToOmegaPi const& candidates,
                                                     aod::TracksWMc const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     aod::McCollisions const& mcCollisions,
                                                     McCollisionsFT0Cs const& collsWithMcLabels,
                                                     BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaPiWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Pi. Cents with FT0C", false);

  void processMcOmegacToOmegaPiWithDCAFitterCentFT0M(aod::HfCandToOmegaPi const& candidates,
                                                     aod::TracksWMc const& tracks,
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
                                                   aod::TracksWMc const& tracks,
                                                   aod::McParticles const& mcParticles,
                                                   aod::McCollisions const& mcCollisions,
                                                   McCollisionsNoCents const& collsWithMcLabels,
                                                   BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaKaWithDCAFitterNoCent, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Ka. No cents", false);

  void processMcOmegacToOmegaKaWithDCAFitterCentFT0C(aod::HfCandToOmegaK const& candidates,
                                                     aod::TracksWMc const& tracks,
                                                     aod::McParticles const& mcParticles,
                                                     aod::McCollisions const& mcCollisions,
                                                     McCollisionsFT0Cs const& collsWithMcLabels,
                                                     BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::FT0C, aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaK>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcOmegacToOmegaKaWithDCAFitterCentFT0C, "Perform MC matching of DCAFitter reconstructed Omegac0 to Omega Ka. Cents with FT0C", false);

  void processMcOmegacToOmegaKaWithDCAFitterCentFT0M(aod::HfCandToOmegaK const& candidates,
                                                     aod::TracksWMc const& tracks,
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
                                              aod::TracksWMc const& tracks,
                                              aod::McParticles const& mcParticles,
                                              aod::McCollisions const& mcCollisions,
                                              McCollisionsNoCents const& collsWithMcLabels,
                                              BCsInfo const& bcs)
  {
    runXic0Omegac0Mc<o2::hf_centrality::CentralityEstimator::None, aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi>(candidates, tracks, mcParticles, collsWithMcLabels, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfCandidateCreatorXic0Omegac0QaMc, processMcXicToXiPiWithKFParticleNoCent, "Perform MC matching of DCAFitter reconstructed Xic0 to Xi Pi. No cents", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0Qa>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorXic0Omegac0QaMc>(cfgc)};
}
