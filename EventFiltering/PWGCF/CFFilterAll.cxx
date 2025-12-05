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

/// \file CFFilterAll.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Laura Serksnyte, TU München, laura.serksnyte@cern.ch; Anton Riedel, TU München, anton.riedel@cern.ch; Maximilian Korwieser, TU Munich, maximilian.korwieser@cern.ch

#include "../filterTables.h"

#include "PWGLF/Utils/strangenessBuilderHelper.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include "TMath.h"

#include "fairlogger/Logger.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace cf_trigger
{
// enums
enum CFTriggers {
  kPPP,
  kPPL,
  kPLL,
  kLLL,
  kPPPhi,
  kPPRho,
  kPD,
  kLD,
  kPhiD,
  kRhoD,
  kNTriggers
};

// variables for track selection
const std::vector<std::string> trackNames{"Pion", "Kaon", "Proton", "Deuteron"};
const uint32_t nTrackNames = 4;

const std::vector<std::string> trackSelectionNames{"AbsEtaMax", "TpcClusterMin", "TpcRowMin", "TpcCrossedOverFoundMin", "TpcSharedMax", "TpcFracSharedMax", "ItsClusterMin", "ItsIbClusterMin", "AbsDcaXyMax", "AbsDcaZMax", "Chi2TpcMax", "Chi2ItsMax"};
const uint32_t nTrackSelectionNames = 12;

const float trackSelectionTable[nTrackNames][nTrackSelectionNames] = {
  {0.85, 90, 80, 0.83, 160, 1, 1, 0, 0.15, 0.15, 99, 99}, // Pion
  {0.85, 90, 80, 0.83, 160, 1, 1, 0, 0.15, 0.15, 99, 99}, // Kaon
  {0.85, 90, 80, 0.83, 160, 1, 1, 0, 0.15, 0.15, 99, 99}, // Proton
  {0.85, 90, 80, 0.83, 160, 1, 1, 0, 0.15, 0.15, 99, 99}, // Deuteron
};

const std::vector<std::string> pidSelectionNames{"ItsMin", "ItsMax", "TpcMin", "TpcMax", "TpcTofMax"};
const uint32_t nPidSelectionNames = 5;

const float pidSelectionTable[nTrackNames][nPidSelectionNames] = {
  {-99, 99, -4, 4, 4},         // Pion
  {-99, 99, -4, 4, 4},         // Kaon
  {-99, 99, -4, 4, 4},         // Proton
  {-3.5, 3.5, -3.5, 3.5, 3.5}, // Deuteron
};

const std::vector<std::string> momentumSelectionNames{"PtMin", "PtMax", "PThres", "UseInnerParam"};
const uint32_t nMomentumSelectionNames = 4;

const float momentumSelectionTable[nTrackNames][nMomentumSelectionNames] = {
  {0, 6, 0.4, -1},    // Pion
  {0, 6, 0.5, -1},    // Kaon
  {0.3, 6, 0.75, -1}, // Proton
  {0.4, 2, 1.2, -1},  // Deuteron
};

// variables for triggers
const std::vector<std::string> filterNames{"PPP", "PPL", "PLL", "LLL", "PPPhi", "PPRho", "PD", "LD", "PhiD", "RhoD"};
const uint32_t nFilterNames = 10;

const std::vector<std::string> switches{"Switch"};
const uint32_t nSwitches = 1;

const float filterTable[nSwitches][nFilterNames]{
  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};

const std::vector<std::string> limitNames{"Tight Limit", "Loose Limit"};
const uint32_t nLimitNames = 2;

const float limitTable[nLimitNames][nFilterNames]{
  {0.6f, 0.6f, 0.6f, 0.6f, 0.6f, 0.6f, 0.5f, 0.5f, 0.5f, 0.5f},
  {1.4f, 1.4f, 1.4f, 1.4f, 1.4f, 1.4f, 1.0f, 1.0f, 1.0f, 1.0f}};

using FullCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using FullCollision = FullCollisions::iterator;

using FullTracks = soa::Join<aod::FullTracks, aod::TracksDCA,
                             aod::TracksIU, TracksCovIU, // needed for in house v0 reconstruction
                             aod::pidTPCFullEl, aod::pidTPCFullPi,
                             aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTPCFullDe, aod::pidTOFFullEl,
                             aod::pidTOFFullPi, aod::pidTOFFullKa,
                             aod::pidTOFFullPr, aod::pidTOFFullDe,
                             aod::pidTOFbeta>;

} // namespace cf_trigger

struct CFFilterAll {

  Produces<aod::CFFilters> tags;

  // Configs for events
  struct : ConfigurableGroup {
    std::string prefix = "EventSel";
    Configurable<float> zvtx{"zvtx", 10.f, "Max. z-Vertex (cm)"};
    Configurable<bool> eventSel{"eventSel", true, "Use sel8"};
  } EventSelection;

  // Configs for tracks
  struct : ConfigurableGroup {
    std::string prefix = "TrackSel";
    Configurable<LabeledArray<float>> trackProperties{"trackProperties",
                                                      {cf_trigger::trackSelectionTable[0],
                                                       cf_trigger::nTrackNames,
                                                       cf_trigger::nTrackSelectionNames,
                                                       cf_trigger::trackNames,
                                                       cf_trigger::trackSelectionNames},
                                                      "Track Selections"};

    Configurable<LabeledArray<float>> momentum{"momentum",
                                               {cf_trigger::momentumSelectionTable[0],
                                                cf_trigger::nTrackNames,
                                                cf_trigger::nMomentumSelectionNames,
                                                cf_trigger::trackNames,
                                                cf_trigger::momentumSelectionNames},
                                               "Momentum Selections"};

    Configurable<LabeledArray<float>> pid{"pid",
                                          {cf_trigger::pidSelectionTable[0],
                                           cf_trigger::nTrackNames,
                                           cf_trigger::nPidSelectionNames,
                                           cf_trigger::trackNames,
                                           cf_trigger::pidSelectionNames},
                                          "PID Selections"};
  } TrackSelections;

  // Configs for V0
  struct : ConfigurableGroup {
    std::string prefix = "V0BuilderOpts";
    Configurable<int> minCrossedRows{"minCrossedRows", 70, "minimum TPC crossed rows for daughter tracks"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.04, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", 0.04, "DCA Pos To PV"};
    Configurable<double> v0cospa{"v0cospa", 0.95, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcav0dau{"dcav0dau", 2.0, "DCA V0 Daughters"};
    Configurable<float> v0radius{"v0radius", 0, "v0radius"};
    Configurable<float> maxDaughterEta{"maxDaughterEta", 5, "Maximum daughter eta (in abs value)"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  } V0BuilderOpts;

  struct : ConfigurableGroup {
    std::string prefix = "FitterOpts";
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "Create tracks version propagated to PCA"};
    Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iteraterions if chi2/chi2old > this"};
    Configurable<float> maxDzIni{"maxDzIni", 1.e9, "reject (if>0) PCA candicate if tracks DZ exceeds threshold"};
    Configurable<float> maxDxyIni{"maxDxyIni", 4., "Same as above for DXY"};
    Configurable<float> maxChi2{"maxChi2", 1.e9, "Maximum chi2"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> weightedFinalPCA{"weightedFinalPCA", false, "Weight final PCA"};
  } FitterOpts;

  struct : ConfigurableGroup {
    std::string prefix = "LambdaSel";
    Configurable<float> ptMin{"ptMin", 0.f, "Minimum transverse momentum of V0"};
    Configurable<float> dcaDaughMax{"dcaDaughMax", 2.f, "Maximum DCA between the V0 daughters"};
    Configurable<float> cpaMin{"cpaMin", 0.95f, "Minimum CPA of V0"};
    Configurable<float> tranRadMin{"tranRadMin", 0.f, "Minimum transverse radius"};
    Configurable<float> tranRadMax{"tranRadMax", 100.f, "Maximum transverse radius"};
    Configurable<float> decVtxMax{"decVtxMax", 100.f, "Maximum distance from primary vertex"};
    Configurable<float> invMassLow{"invMassLow", 1.05, "Lower limit of the V0 invariant mass"};
    Configurable<float> invMassUp{"invMassUp", 1.18, "Upper limit of the V0 invariant mass"};
    Configurable<bool> rejectKaons{"rejectKaons", true, "Switch to reject kaons"};
    Configurable<float> invKaonMassLow{"invKaonMassLow", 0.49, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> invKaonMassUp{"invKaonMassUp", 0.505, "Upper limit of the V0 invariant mass for Kaon rejection"};
  } LambdaSelections;

  struct : ConfigurableGroup {
    std::string prefix = "LambdaDaughterSel";
    Configurable<float> absEtaMax{"absEtaMax", 0.85f, "V0 Daugh sel: max eta"};
    Configurable<float> tpcClusterMin{"tpcClusterMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
    Configurable<float> dcaMin{"dcaMin", 0.04f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<float> tpcMax{"tpcMax", 5, "PID selections for Lambda daughters"};
  } LambdaDaughterSelections;

  struct : ConfigurableGroup {
    std::string prefix = "PhiSel";
    Configurable<float> invMassLow{"invMassLow", 1.011461, "Lower limit of Phi invariant mass"};
    Configurable<float> invMassUp{"invMassUp", 1.027461, "Upper limit of Phi invariant mass"};
    Configurable<float> tightInvMassLow{"tightInvMassLow", 1.011461, "Lower tight limit of Phi invariant mass"};
    Configurable<float> tightInvMassUp{"tightInvMassUp", 1.027461, "Upper tight limit of Phi invariant mass"};
  } PhiSelections;

  struct : ConfigurableGroup {
    std::string prefix = "RhoSel";
    Configurable<float> invMassLow{"invMassLow", 0.7, "Lower limit of Rho invariant mass"};
    Configurable<float> invMassUp{"invMassUp", 0.85, "Upper limit of Rho invariant mass"};
    Configurable<float> ptLow{"ptLow", 3, "Lower pt limit for rho"};
    Configurable<float> tightInvMassLow{"tightInvMassLow", 0.73, "Lower tight limit of Rho invariant mass"};
    Configurable<float> tightInvMassUp{"tightInvMassUp", 0.82, "Upper tight limit of Rho invariant mass"};
  } RhoSelections;

  // Trigger selections
  struct : ConfigurableGroup {
    std::string prefix = "Triggers";
    Configurable<LabeledArray<float>> filterSwitches{"filterSwitches",
                                                     {cf_trigger::filterTable[0],
                                                      cf_trigger::nSwitches,
                                                      cf_trigger::nFilterNames,
                                                      cf_trigger::switches,
                                                      cf_trigger::filterNames},
                                                     "Switch for triggers"};
    Configurable<LabeledArray<float>> limits{"limits",
                                             {cf_trigger::limitTable[0],
                                              cf_trigger::nLimitNames,
                                              cf_trigger::nFilterNames,
                                              cf_trigger::limitNames,
                                              cf_trigger::filterNames},
                                             "Limits for trigger. Tight limit without downsampling and loose with downsampling"};
  } TriggerSelections;

  struct : ConfigurableGroup {
    std::string prefix = "Binning";
    ConfigurableAxis multiplicity{"multiplicity", {200, 0, 200}, "Binning Multiplicity"};
    ConfigurableAxis zvtx{"zvtx", {30, -15, 15}, "Binning Zvertex"};

    ConfigurableAxis momentum{"momentum", {600, 0, 6}, "Binning Momentum"};
    ConfigurableAxis eta{"eta", {200, -1, 1}, "Binning eta"};
    ConfigurableAxis phi{"phi", {720, 0, o2::constants::math::TwoPI}, "Binning phi"};
    ConfigurableAxis dca{"dca", {100, -0.2, 0.2}, "Binning Dca"};

    ConfigurableAxis nsigma{"nsigma", {500, -5, 5}, "Binning nsigma"};
    ConfigurableAxis nsigmaComb{"nsigmaComb", {500, 0, 5}, "Binning nsigma comb"};
    ConfigurableAxis tpcSignal{"tpcSignal", {500, 0, 500}, "Binning Tpc Signal"};
    ConfigurableAxis itsSignal{"itsSignal", {150, 0, 15}, "Binning Its Signal"};
    ConfigurableAxis tofSignal{"tofSignal", {120, 0, 1.2}, "Binning Tof Signal"};
    ConfigurableAxis tpcCluster{"tpcCluster", {153, 0, 153}, "Binning Tpc Clusters"};
    ConfigurableAxis tpcChi2{"tpcChi2", {100, 0, 5}, "Binning TPC chi2"};
    ConfigurableAxis itsCluster{"itsCluster", {8, -0.5, 7.5}, "Binning Its Clusters"};
    ConfigurableAxis itsIbCluster{"itsIbCluster", {4, -0.5, 3.5}, "Binning Its Inner Barrel Clusters"};
    ConfigurableAxis itsChi2{"itsChi2", {100, 0, 50}, "Binning ITS chi2"};

    ConfigurableAxis momCor{"momCor", {100, -1, 1}, "Binning Ratios"};
    ConfigurableAxis ratio{"ratio", {200, 0, 2}, "Binning Ratios"};

    ConfigurableAxis invMassLambda{"invMassLambda", {200, 1, 1.2}, "Binning Invariant Mass Lambda"};
    ConfigurableAxis invMassK0short{"invMassK0short", {100, 0.48, 0.52}, "Binning Invariant Mass K0short"};

    ConfigurableAxis dcaDaugh{"dcaDaugh", {200, 0, 2}, "Binning daughter DCA at decay vertex"};
    ConfigurableAxis cpa{"cpa", {100, 0.9, 1}, "Binning CPA"};
    ConfigurableAxis transRad{"transRad", {100, 0, 100}, "Binning Transverse Radius"};
    ConfigurableAxis decayVtx{"decayVtx", {100, 0, 100}, "Binning Decay Vertex"};

    ConfigurableAxis invMassPhi{"invMassPhi", {600, 0.7, 1.3}, "Binning Invariant Mass Phi"};

    ConfigurableAxis invMassRho{"invMassRho", {600, 0.47, 1.07}, "Binning Invariant Mass Rho"};

    ConfigurableAxis q3{"q3", {300, 0, 3}, "Binning Decay Q3"};
    ConfigurableAxis kstar{"kstar", {300, 0, 3}, "Binning Decay Kstar"};

  } Binning;

  // define histogram registry
  // because we have so many histograms, we need to have 2 registries
  HistogramRegistry registryParticleQA{"ParticleQA", {}, OutputObjHandlingPolicy::AnalysisObject}; // for particle histograms
  HistogramRegistry registryTriggerQA{"TriggerQA", {}, OutputObjHandlingPolicy::AnalysisObject};   // for trigger histograms

  // helper object flor building lambdas
  o2::pwglf::strangenessBuilderHelper mStraHelper;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0.;

  // 4vectors for all particles
  std::vector<ROOT::Math::PtEtaPhiMVector> vecProton, vecAntiProton, vecDeuteron, vecAntiDeuteron, vecLambda, vecAntiLambda, vecKaon, vecAntiKaon, vecPhi, vecPion, vecAntiPion, vecRho;
  // indices for all particles
  std::vector<int> idxProton, idxAntiProton, idxDeuteron, idxAntiDeuteron, idxKaon, idxAntiKaon, idxPion, idxAntiPion;
  // indices for lambda daughters
  std::vector<int> idxLambdaDaughProton, idxLambdaDaughPion, idxAntiLambdaDaughProton, idxAntiLambdaDaughPion, idxPhiDaughPos, idxPhiDaughNeg, idxRhoDaughPos, idxRhoDaughNeg;

  // arrays to store found pairs/tripplets and trigger decisions
  std::array<bool, cf_trigger::kNTriggers> keepEventTightLimit;
  std::array<bool, cf_trigger::kNTriggers> keepEventLooseLimit;
  std::array<int, cf_trigger::kNTriggers> signalTightLimit;
  std::array<int, cf_trigger::kNTriggers> signalLooseLimit;

  void init(o2::framework::InitContext&)
  {

    // setup strangeness builder
    mStraHelper.v0selections.minCrossedRows = V0BuilderOpts.minCrossedRows.value;
    mStraHelper.v0selections.dcanegtopv = V0BuilderOpts.dcanegtopv.value;
    mStraHelper.v0selections.dcapostopv = V0BuilderOpts.dcapostopv.value;
    mStraHelper.v0selections.v0cospa = V0BuilderOpts.v0cospa.value;
    mStraHelper.v0selections.dcav0dau = V0BuilderOpts.dcav0dau.value;
    mStraHelper.v0selections.v0radius = V0BuilderOpts.v0radius.value;
    mStraHelper.v0selections.maxDaughterEta = V0BuilderOpts.maxDaughterEta.value;

    mStraHelper.fitter.setPropagateToPCA(FitterOpts.propagateToPCA.value);
    mStraHelper.fitter.setMaxR(FitterOpts.maxR.value);
    mStraHelper.fitter.setMinParamChange(FitterOpts.minParamChange.value);
    mStraHelper.fitter.setMinRelChi2Change(FitterOpts.minRelChi2Change.value);
    mStraHelper.fitter.setMaxDZIni(FitterOpts.maxDzIni.value);
    mStraHelper.fitter.setMaxDXYIni(FitterOpts.maxDxyIni.value);
    mStraHelper.fitter.setMaxChi2(FitterOpts.maxChi2.value);
    mStraHelper.fitter.setUseAbsDCA(FitterOpts.useAbsDCA.value);
    mStraHelper.fitter.setWeightedFinalPCA(FitterOpts.weightedFinalPCA.value);

    // setup histograms
    int allTriggers = 2 * cf_trigger::nFilterNames;
    int prossedEventsBins = 3 + allTriggers;
    std::vector<std::string> triggerTitles = {"ppp_LooseQ3", "ppp_TightQ3",
                                              "ppL_LooseQ3", "ppL_TightQ3",
                                              "pLL_LooseQ3", "pLL_TightQ3",
                                              "LLL_LooseQ3", "LLL_TightQ3",
                                              "ppPhi_LooseQ3", "ppPhi_TightQ3",
                                              "ppRho_LooseQ3", "ppRho_TightQ3",
                                              "pD_LooseKstar", "pD_TightKstar",
                                              "LD_LooseKstar", "LD_TightKstar",
                                              "PhiD_LooseKstar", "PhiD_TightKstar",
                                              "RhoD_LooseKstar", "RhoD_TightKstar"};

    registryTriggerQA.add("fProcessedEvents", "CF - event filtered;;Events", HistType::kTH1F, {{prossedEventsBins, -0.5, prossedEventsBins - 0.5}});
    registryTriggerQA.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(1, "all");
    registryTriggerQA.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(2, "accepted_loose");
    registryTriggerQA.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(3, "accepted_tight");

    registryTriggerQA.add("fTriggerCorrelations", "CF - Trigger correlations", HistType::kTH2F, {{allTriggers, -0.5, allTriggers - 0.5}, {allTriggers, -0.5, allTriggers - 0.5}});

    for (size_t iBin = 0; iBin < triggerTitles.size(); iBin++) {
      registryTriggerQA.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 4, triggerTitles[iBin].data()); // start triggers from 4th bin
      registryTriggerQA.get<TH2>(HIST("fTriggerCorrelations"))->GetXaxis()->SetBinLabel(iBin + 1, triggerTitles[iBin].data());
      registryTriggerQA.get<TH2>(HIST("fTriggerCorrelations"))->GetYaxis()->SetBinLabel(iBin + 1, triggerTitles[iBin].data());
    }

    // event cuts
    registryParticleQA.add("EventQA/Before/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryParticleQA.add("EventQA/Before/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});

    registryParticleQA.add("EventQA/After/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryParticleQA.add("EventQA/After/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});

    // all tracks before cuts
    registryParticleQA.add("TrackQA/Before/Particle/fPt", "Transverse;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/Before/Particle/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/Before/Particle/fPhi", "Azimuthal;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("TrackQA/Before/Particle/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/Before/Particle/fItsSignal", "ITSSignal;p_{TPC} (GeV/c);ITS Signal", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/Before/Particle/fTpcSignal", "TPCSignal;p_{TPC} (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/Before/Particle/fTofSignal", "TOFSignal;p_{TPC} (GeV/c);TOF Signal", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/Before/AntiParticle/fPt", "Transverse momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/Before/AntiParticle/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/Before/AntiParticle/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("TrackQA/Before/AntiParticle/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/Before/AntiParticle/fItsSignal", "ITSSignal;p_{TPC} (GeV/c);ITS Signal", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/Before/AntiParticle/fTpcSignal", "TPCSignal;p_{TPC} (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/Before/AntiParticle/fTofSignal", "TOFSignal;p_{TPC} (GeV/c);TOF Signal", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    // PID vs momentum before cuts
    registryParticleQA.add("TrackQA/Before/Pion/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Pion/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Pion/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Pion/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/AntiPion/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiPion/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiPion/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiPion/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/Kaon/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Kaon/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Kaon/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Kaon/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/AntiKaon/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiKaon/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiKaon/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiKaon/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/Proton/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Proton/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Proton/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Proton/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/AntiProton/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiProton/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiProton/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiProton/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/Deuteron/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Deuteron/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Deuteron/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/Deuteron/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/Before/AntiDeuteron/fNsigmaITS", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiDeuteron/fNsigmaTPC", "NSigmaTPC;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiDeuteron/fNsigmaTOF", "NSigmaTOF;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/Before/AntiDeuteron/fNsigmaTPCTOF", "NsigmaTPCTOF;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    // Pion
    registryParticleQA.add("TrackQA/After/Pion/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Pion/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Pion/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/Pion/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/Pion/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/Pion/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Pion/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Pion/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Pion/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/Pion/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/Pion/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/Pion/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/Pion/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/Pion/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/Pion/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Pion/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Pion/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Pion/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Pion/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Pion/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/Pion/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/Pion/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/Pion/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // antiPion
    registryParticleQA.add("TrackQA/After/AntiPion/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiPion/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiPion/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/AntiPion/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/AntiPion/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/AntiPion/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiPion/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiPion/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiPion/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/AntiPion/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/AntiPion/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/AntiPion/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/AntiPion/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/AntiPion/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/AntiPion/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiPion/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiPion/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiPion/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiPion/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiPion/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/AntiPion/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/AntiPion/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/AntiPion/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // Kaon
    registryParticleQA.add("TrackQA/After/Kaon/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Kaon/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Kaon/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/Kaon/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/Kaon/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/Kaon/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Kaon/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Kaon/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Kaon/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/Kaon/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/Kaon/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/Kaon/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/Kaon/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/Kaon/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/Kaon/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Kaon/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Kaon/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Kaon/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Kaon/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Kaon/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/Kaon/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/Kaon/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/Kaon/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // antiKaon
    registryParticleQA.add("TrackQA/After/AntiKaon/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiKaon/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiKaon/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/AntiKaon/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/AntiKaon/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/AntiKaon/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiKaon/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiKaon/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiKaon/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/AntiKaon/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/AntiKaon/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/AntiKaon/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiKaon/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/AntiKaon/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/AntiKaon/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/AntiKaon/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // proton
    registryParticleQA.add("TrackQA/After/Proton/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Proton/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Proton/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/Proton/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/Proton/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/Proton/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Proton/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Proton/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Proton/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/Proton/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/Proton/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/Proton/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/Proton/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/Proton/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/Proton/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Proton/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Proton/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Proton/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Proton/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Proton/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/Proton/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/Proton/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/Proton/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // antiproton
    registryParticleQA.add("TrackQA/After/AntiProton/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiProton/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiProton/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/AntiProton/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/AntiProton/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/AntiProton/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiProton/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiProton/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiProton/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/AntiProton/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/AntiProton/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/AntiProton/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/AntiProton/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/AntiProton/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/AntiProton/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiProton/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiProton/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiProton/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiProton/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiProton/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/AntiProton/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/AntiProton/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/AntiProton/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // Deuteron
    registryParticleQA.add("TrackQA/After/Deuteron/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Deuteron/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/Deuteron/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/Deuteron/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/Deuteron/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/Deuteron/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Deuteron/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Deuteron/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/Deuteron/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/Deuteron/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/Deuteron/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/Deuteron/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/Deuteron/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/Deuteron/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/Deuteron/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Deuteron/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Deuteron/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/Deuteron/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Deuteron/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/Deuteron/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/Deuteron/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/Deuteron/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/Deuteron/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // AntiDeuteron
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fPTpc", "Momentum at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fMomCor", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {Binning.momentum, Binning.momCor}});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fPhi", "Azimuthal angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    registryParticleQA.add("TrackQA/After/AntiDeuteron/fNsigmaIts", "NSigmaITS;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fNsigmaTpc", "NSigmaTPC;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fNsigmaTof", "NSigmaTOF;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fNsigmaTpcTof", "NSigmaTPCTOF;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {Binning.momentum, Binning.nsigmaComb}});

    registryParticleQA.add("TrackQA/After/AntiDeuteron/fItsSignal", "ITS Signal;p (GeV/c);<cluster size * cos(#lambda)> (cm)", {HistType::kTH2F, {Binning.momentum, Binning.itsSignal}});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcSignal", "TPC Signal;p (GeV/c);TPC Signal", {HistType::kTH2F, {Binning.momentum, Binning.tpcSignal}});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTofBeta", "TOF #beta;p (GeV/c);#beta_{TOF}", {HistType::kTH2F, {Binning.momentum, Binning.tofSignal}});

    registryParticleQA.add("TrackQA/After/AntiDeuteron/fDcaXy", "DCA_{xy};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fDcaZ", "DCA_{z};p_{T} (GeV/c); DCA_{Z};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});

    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcCrossedRows", "TPC Crossed Rows;TPC Crossed Rows;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcSharedClusters", "TPC Shared Clusters;TPC Shared Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcSharedClusterOverClusterss", "TPC Shared Clusters/Clusters;TPC Shared Clusters/Clusters;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcFindableOverRows", "TPC Findabled/Crossed Rows;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {Binning.ratio});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fTpcChi2OverCluster", "TPC #chi^{2}/Cluster;TPC #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.tpcChi2});

    registryParticleQA.add("TrackQA/After/AntiDeuteron/fItsClusters", "ITS Clusters;ITS Clusters;Entries", HistType::kTH1F, {Binning.itsCluster});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fItsIbClusters", "ITS Inner Barrel Clusters;ITS Inner Barrel Clusters;Entries", HistType::kTH1F, {Binning.itsIbCluster});
    registryParticleQA.add("TrackQA/After/AntiDeuteron/fItsChi2OverCluster", "ITS #chi^{2}/Cluster;ITS #chi^{2}/Cluster;Entries", HistType::kTH1F, {Binning.itsChi2});

    // Lambda before
    registryParticleQA.add("LambdaQA/Before/fPt", "Transverse momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/Before/fEta", "Psedurapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/Before/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/Before/fInvMassLambda", "Invariant mass Lambda;M_{#pi p};Entries", HistType::kTH1F, {Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/Before/fInvMassAntiLambda", "Invariant mass AntiLambda;M_{#pi p};Entries", HistType::kTH1F, {Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/Before/fInvMassLambdaVsAntiLambda", "Invariant mass of Lambda vs AntiLambda;M_{#pi p};Entries", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/Before/fInvMassLambdaVsKaon", "Invariant mass of Lambda vs K0;M_{#pi p};;M_{#pi #pi}", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassK0short});
    registryParticleQA.add("LambdaQA/Before/fInvMassAntiLambdaVsKaon", "Invariant mass of AntiLambda vs K0;M_{#pi p};;M_{#pi #pi}", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassK0short});
    registryParticleQA.add("LambdaQA/Before/fDcaDaugh", "DCA_{Daugh};DCA_{daugh};Entries", HistType::kTH1F, {Binning.dcaDaugh});
    registryParticleQA.add("LambdaQA/Before/fCpa", "Cosine of pointing angle;CPA;Entries", HistType::kTH1F, {Binning.cpa});
    registryParticleQA.add("LambdaQA/Before/fTranRad", "Transverse Radisu;TranRad;Entries", HistType::kTH1F, {Binning.transRad});
    registryParticleQA.add("LambdaQA/Before/fDecVtx", "Decay vertex displacement;DecVtx;Entries", HistType::kTH1F, {Binning.decayVtx});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fDcaXy", "DCA_{XY};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fNsigmaTpcProton", "NSigmaTPC Proton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("LambdaQA/Before/PosDaughter/fNsigmaTpcPion", "NSigmaTPC Pion;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fDcaXy", "DCA_{XY};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fNsigmaTpcProton", "NSigmaTPC AnitProton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("LambdaQA/Before/NegDaughter/fNsigmaTpcPion", "NSigmaTPC AntiPion;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});

    // Lambda after
    registryParticleQA.add("LambdaQA/After/Lambda/fPt", "Transverse momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/After/Lambda/fEta", "Psedurapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/After/Lambda/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/After/Lambda/fInvMass", "Invariant mass;M_{#pi p};Entries", HistType::kTH1F, {Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/After/Lambda/fInvMassLambdaVsAntiLambda", "Invariant mass of Lambda vs AntiLambda;M_{#pi p};Entries", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/After/Lambda/fInvMassLambdaVsKaon", "Invariant mass of rejected K0 vs V0s;M_{#pi p};;M_{#pi #pi}", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassK0short});
    registryParticleQA.add("LambdaQA/After/Lambda/fDcaDaugh", "DCA_{Daugh};DCA_{daugh};Entries", HistType::kTH1F, {Binning.dcaDaugh});
    registryParticleQA.add("LambdaQA/After/Lambda/fCpa", "Cosine of pointing angle;CPA;Entries", HistType::kTH1F, {Binning.cpa});
    registryParticleQA.add("LambdaQA/After/Lambda/fTranRad", "Transverse Radisu;TranRad;Entries", HistType::kTH1F, {Binning.transRad});
    registryParticleQA.add("LambdaQA/After/Lambda/fDecVtx", "Decay vertex displacement;DecVtx;Entries", HistType::kTH1F, {Binning.decayVtx});
    registryParticleQA.add("LambdaQA/After/Lambda/PosDaughter/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/After/Lambda/PosDaughter/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/After/Lambda/PosDaughter/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/After/Lambda/PosDaughter/fDcaXy", "DCA_{XY};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("LambdaQA/After/Lambda/PosDaughter/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("LambdaQA/After/Lambda/PosDaughter/fNsigmaTpc", "NSigmaTPC Proton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("LambdaQA/After/Lambda/NegDaughter/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/After/Lambda/NegDaughter/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/After/Lambda/NegDaughter/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/After/Lambda/NegDaughter/fDcaXy", "DCA_{XY};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("LambdaQA/After/Lambda/NegDaughter/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("LambdaQA/After/Lambda/NegDaughter/fNsigmaTpc", "NSigmaTPC AntiPion;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});

    // AntiLambda after
    registryParticleQA.add("LambdaQA/After/AntiLambda/fPt", "Transverse momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fEta", "Psedurapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fInvMass", "Invariant mass;M_{#pi p};Entries", HistType::kTH1F, {Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fInvMassAntiLambdaVsLambda", "Invariant mass of Lambda vs AntiLambda;M_{#pi p};Entries", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassLambda});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fInvMassAntiLambdaVsKaon", "Invariant mass of rejected K0 vs V0s;M_{#pi p};;M_{#pi #pi}", HistType::kTH2F, {Binning.invMassLambda, Binning.invMassK0short});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fDcaDaugh", "DCA_{Daugh};DCA_{daugh};Entries", HistType::kTH1F, {Binning.dcaDaugh});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fCpa", "Cosine of pointing angle;CPA;Entries", HistType::kTH1F, {Binning.cpa});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fTranRad", "Transverse Radisu;TranRad;Entries", HistType::kTH1F, {Binning.transRad});
    registryParticleQA.add("LambdaQA/After/AntiLambda/fDecVtx", "Decay vertex displacement;DecVtx;Entries", HistType::kTH1F, {Binning.decayVtx});
    registryParticleQA.add("LambdaQA/After/AntiLambda/PosDaughter/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/After/AntiLambda/PosDaughter/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/After/AntiLambda/PosDaughter/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/After/AntiLambda/PosDaughter/fDcaXy", "DCA_{XY};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("LambdaQA/After/AntiLambda/PosDaughter/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("LambdaQA/After/AntiLambda/PosDaughter/fNsigmaTpc", "NSigmaTPC Proton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});
    registryParticleQA.add("LambdaQA/After/AntiLambda/NegDaughter/fPt", "Transverse Momentum;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("LambdaQA/After/AntiLambda/NegDaughter/fEta", "Pseudorapidity;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("LambdaQA/After/AntiLambda/NegDaughter/fPhi", "Azimuthal Angle;#varphi;Entries", HistType::kTH1F, {Binning.phi});
    registryParticleQA.add("LambdaQA/After/AntiLambda/NegDaughter/fDcaXy", "DCA_{XY};p_{T} (GeV/c); DCA_{XY};Entries", HistType::kTH2F, {Binning.momentum, Binning.dca});
    registryParticleQA.add("LambdaQA/After/AntiLambda/NegDaughter/fTpcClusters", "TPC Clusters;TPC Clusters;Entries", HistType::kTH1F, {Binning.tpcCluster});
    registryParticleQA.add("LambdaQA/After/AntiLambda/NegDaughter/fNsigmaTpc", "NSigmaTPC AntiPion;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {Binning.momentum, Binning.nsigma}});

    // Phi before
    registryParticleQA.add("PhiQA/Before/fInvMass", "Invariant mass #phi;M_{KK};Entries", HistType::kTH1F, {Binning.invMassPhi});
    registryParticleQA.add("PhiQA/Before/fPt", "Transverse momentum #phi;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("PhiQA/Before/fEta", "Pseudorapidity of V0;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("PhiQA/Before/fPhi", "Azimuthal angle of #phi;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    // Phi after
    registryParticleQA.add("PhiQA/After/fInvMass", "Invariant mass #phi;M_{KK};Entries", HistType::kTH1F, {Binning.invMassPhi});
    registryParticleQA.add("PhiQA/After/fPt", "Transverse momentum #phi;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("PhiQA/After/fEta", "Pseudorapidity of #phi;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("PhiQA/After/fPhi", "Azimuthal angle of #Phi;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    // Rho before
    registryParticleQA.add("RhoQA/Before/fInvMass", "Invariant mass #rho;M_{#pi#pi};Entries", HistType::kTH1F, {Binning.invMassRho});
    registryParticleQA.add("RhoQA/Before/fPt", "Transverse momentum #rho;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("RhoQA/Before/fEta", "Pseudorapidity of #rho;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("RhoQA/Before/fPhi", "Azimuthal angle of #rho;#varphi;Entries", HistType::kTH1F, {Binning.phi});

    // Rho after
    registryParticleQA.add("RhoQA/After/fInvMass", "Invariant mass #rho;M_{#pi#pi};Entries", HistType::kTH1F, {Binning.invMassRho});
    registryParticleQA.add("RhoQA/After/fPt", "Transverse momentum #rho;p_{T} (GeV/c);Entries", HistType::kTH1F, {Binning.momentum});
    registryParticleQA.add("RhoQA/After/fEta", "Pseudorapidity of #rho;#eta;Entries", HistType::kTH1F, {Binning.eta});
    registryParticleQA.add("RhoQA/After/fPhi", "Azimuthal angle of #rho;#phi;Entries", HistType::kTH1F, {Binning.phi});

    // for ppp
    registryTriggerQA.add("PPP/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPP/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPP/all/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPP/all/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPP/all/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPP/all/fAntiProtonQ3VsPt", "Q_{3} vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPP/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPP/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPP/loose/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPP/loose/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPP/loose/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPP/loose/fAntiProtonQ3VsPt", "Q_{3} vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPP/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPP/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPP/tight/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPP/tight/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPP/tight/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPP/tight/fAntiProtonQ3VsPt", "Q_{3} vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});

    // for ppl
    registryTriggerQA.add("PPL/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPL/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPL/all/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPL/all/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPL/all/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/all/fAntiProtonQ3VsPt", "Q_{3} vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/all/fLambdaQ3VsPt", "Q_{3} vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/all/fAntiLambdaQ3VsPt", "Q_{3} vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPL/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPL/loose/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPL/loose/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPL/loose/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/loose/fAntiProtonQ3VsPt", "Q_{3} vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/loose/fLambdaQ3VsPt", "Q_{3} vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/loose/fAntiLambdaQ3VsPt", "Q_{3} vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPL/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPL/tight/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPL/tight/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPL/tight/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/tight/fAntiProtonQ3VsPt", "Q_{3} vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/tight/fLambdaQ3VsPt", "Q_{3} vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PPL/tight/fAntiLambdaQ3VsPt", "Q_{3} vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});

    // for pll
    registryTriggerQA.add("PLL/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PLL/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PLL/all/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PLL/all/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PLL/all/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/all/fAntiProtonQ3VsPt", "Q3 vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/all/fLambdaQ3VsPt", "Q3 vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/all/fAntiLambdaQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PLL/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PLL/loose/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PLL/loose/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PLL/loose/fProtonQ3VsPt", "Q3 vs pT vs Proton p;Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/loose/fAntiProtonQ3VsPt", "Q3 vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/loose/fLambdaQ3VsPt", "Q3 vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/loose/fAntiLambdaQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PLL/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PLL/tight/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PLL/tight/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PLL/tight/fProtonQ3VsPt", "Q3 vs pT vs Proton p;Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/tight/fAntiProtonQ3VsPt", "Q3 vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/tight/fLambdaQ3VsPt", "Q3 vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("PLL/tight/fAntiLambdaQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});

    // for lll
    registryTriggerQA.add("LLL/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("LLL/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("LLL/all/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("LLL/all/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("LLL/all/fLambdaQ3VsPt", "Q3 vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("LLL/all/fAntiLambdaQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("LLL/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("LLL/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("LLL/loose/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("LLL/loose/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("LLL/loose/fLambdaQ3VsPt", "Q3 vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("LLL/loose/fAntiLambdaQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("LLL/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("LLL/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("LLL/tight/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("LLL/tight/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("LLL/tight/fLambdaQ3VsPt", "Q3 vs Lambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});
    registryTriggerQA.add("LLL/tight/fAntiLambdaQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.q3, Binning.momentum}});

    // for ppPhi
    registryTriggerQA.add("PPPhi/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPPhi/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPPhi/all/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPPhi/all/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPPhi/all/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/all/fAntiProtonQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/all/fPhiQ3VsPt", "Q_{3} vs #phi p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/all/fPhiQ3VsInvMass", "Q_{3} vs #phi mass;Q_{3} (GeV/c);M_{K^{+}K^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.q3, Binning.invMassPhi});
    registryTriggerQA.add("PPPhi/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPPhi/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPPhi/loose/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPPhi/loose/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPPhi/loose/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/loose/fAntiProtonQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/loose/fPhiQ3VsPt", "Q_{3} vs #phi p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/loose/fPhiQ3VsInvMass", "Q_{3} vs #phi mass;Q_{3} (GeV/c);M_{K^{+}K^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.q3, Binning.invMassPhi});
    registryTriggerQA.add("PPPhi/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPPhi/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPPhi/tight/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPPhi/tight/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);Entries", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPPhi/tight/fProtonQ3VsPt", "Q_{3} vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/tight/fAntiProtonQ3VsPt", "Q3 vs AntiLambda p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/tight/fPhiQ3VsPt", "Q_{3} vs #phi p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPPhi/tight/fPhiQ3VsInvMass", "Q_{3} vs #phi mass;Q_{3} (GeV/c);M_{K^{+}K^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.q3, Binning.invMassPhi});

    // for ppRho
    registryTriggerQA.add("PPRho/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPRho/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPRho/all/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPRho/all/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPRho/all/fProtonQ3VsPt", "Q3 vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/all/fAntiProtonQ3VsPt", "Q3 vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/all/fRhoQ3VsPt", "Q3 vs #rho p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/all/fRhoQ3VsInvMass", "Q_{3} vs #rho mass;Q_{3} (GeV/c);M_{#pi^{+}#pi^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.q3, Binning.invMassRho});
    registryTriggerQA.add("PPRho/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPRho/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPRho/loose/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPRho/loose/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPRho/loose/fProtonQ3VsPt", "Q3 vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/loose/fAntiProtonQ3VsPt", "Q3 vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/loose/fRhoQ3VsPt", "Q3 vs #rho p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/loose/fRhoQ3VsInvMass", "Q_{3} vs #rho mass;Q_{3} (GeV/c);M_{#pi^{+}#pi^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.q3, Binning.invMassRho});
    registryTriggerQA.add("PPRho/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PPRho/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PPRho/tight/fSE_particle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPRho/tight/fSE_antiparticle", "Same Event distribution;Q_{3} (GeV/c);SE", HistType::kTH1F, {Binning.q3});
    registryTriggerQA.add("PPRho/tight/fProtonQ3VsPt", "Q3 vs Proton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/tight/fAntiProtonQ3VsPt", "Q3 vs AntiProton p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/tight/fRhoQ3VsPt", "Q3 vs #rho p_{T};Q_{3} (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {Binning.q3, Binning.momentum});
    registryTriggerQA.add("PPRho/tight/fRhoQ3VsInvMass", "Q_{3} vs #rho mass;Q_{3} (GeV/c);M_{#pi^{+}#pi^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.q3, Binning.invMassRho});

    // for pd
    registryTriggerQA.add("PD/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PD/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PD/all/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PD/all/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PD/all/fProtonKstarVsPt", "k* vs Proton p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/all/fAntiProtonKstarVsPt", "k* vs AntiProton p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/all/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/all/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PD/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PD/loose/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PD/loose/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PD/loose/fProtonKstarVsPt", "k* vs Proton p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/loose/fAntiProtonKstarVsPt", "k* vs  AntiProton p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/loose/fDeuteronKstarVsPt", "k* Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/loose/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PD/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PD/tight/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PD/tight/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PD/tight/fProtonKstarVsPt", "k* vs Proton p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/tight/fAntiProtonKstarVsPt", "k* vs AntiProton p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/tight/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PD/tight/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});

    // for ld
    registryTriggerQA.add("LD/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("LD/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("LD/all/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("LD/all/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("LD/all/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/all/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/all/fLambdaKstarVsPt", "k* vs Lambda p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/all/fAntiLambdaKstarVsPt", "k* vs AntiLambda p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("LD/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("LD/loose/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("LD/loose/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("LD/loose/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/loose/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/loose/fLambdaKstarVsPt", "k* vs Lambda p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/loose/fAntiLambdaKstarVsPt", "k* vs AntiLambda p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("LD/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("LD/tight/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("LD/tight/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("LD/tight/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/tight/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/tight/fLambdaKstarVsPt", "k* vs Lambda p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("LD/tight/fAntiLambdaKstarVsPt", "k* vs AntiLambda p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});

    // for phid
    registryTriggerQA.add("PhiD/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PhiD/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PhiD/all/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PhiD/all/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PhiD/all/fPhiKstarVsPt", "k* vs Phi p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/all/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/all/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/all/fPhiKstarVsInvMass", "k* vs #phi mass;k* (GeV/c);M_{K^{+}K^{-}} (GeV/c^{2});", HistType::kTH2F, {Binning.kstar, Binning.invMassPhi});
    registryTriggerQA.add("PhiD/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PhiD/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PhiD/loose/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PhiD/loose/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PhiD/loose/fPhiKstarVsPt", "k* vs Phi p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/loose/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/loose/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/loose/fPhiKstarVsInvMass", "k* vs #phi mass;k* (GeV/c);M_{K^{+}K^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.kstar, Binning.invMassPhi});
    registryTriggerQA.add("PhiD/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("PhiD/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("PhiD/tight/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PhiD/tight/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("PhiD/tight/fPhiKstarVsPt", "k* vs Phi p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/tight/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/tight/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("PhiD/tight/fPhiKstarVsInvMass", "k* vs #phi mass;k* (GeV/c);M_{K^{+}K^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.kstar, Binning.invMassPhi});

    // for rhod
    registryTriggerQA.add("RhoD/all/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("RhoD/all/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("RhoD/all/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("RhoD/all/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("RhoD/all/fRhoKstarVsPt", "k* vs Rho p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/all/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/all/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/all/fRhoKstarVsInvMass", "k* vs #rho mass;k* (GeV/c);M_{#pi^{+}#pi^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.kstar, Binning.invMassRho});
    registryTriggerQA.add("RhoD/loose/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("RhoD/loose/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("RhoD/loose/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("RhoD/loose/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("RhoD/loose/fRhoKstarVsPt", "k* vs Rho p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/loose/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/loose/fAntiDeuteronKstarVsPt", "k* vs AntiDeuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/loose/fRhoKstarVsInvMass", "k* vs #rho mass;k* (GeV/c);M_{#pi^{+}#pi^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.kstar, Binning.invMassRho});
    registryTriggerQA.add("RhoD/tight/fMultiplicity", "Multiplicity;Mult;Entries", HistType::kTH1F, {Binning.multiplicity});
    registryTriggerQA.add("RhoD/tight/fZvtx", "Zvtx;Z_{vtx};Entries", HistType::kTH1F, {Binning.zvtx});
    registryTriggerQA.add("RhoD/tight/fSE_particle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("RhoD/tight/fSE_antiparticle", "Same Event distribution;k* (GeV/c);Entries", HistType::kTH1F, {Binning.kstar});
    registryTriggerQA.add("RhoD/tight/fRhoKstarVsPt", "k* vs Rho p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/tight/fDeuteronKstarVsPt", "k* vs Deuteron p_{T};k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/tight/fAntiDeuteronKstarVsPt", "AntiDeuteron p_{T} vs k*;k* (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {Binning.kstar, Binning.momentum}});
    registryTriggerQA.add("RhoD/tight/fRhoKstarVsInvMass", "k* vs #rho mass;k* (GeV/c);M_{#pi^{+}#pi^{-}} (GeV/c^{2})", HistType::kTH2F, {Binning.kstar, Binning.invMassRho});
  }

  void initCCDB(int run)
  {
    if (run != mRunNumber) {
      mRunNumber = run;
      o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", run);
      o2::base::Propagator::initFieldFromGRP(grpmag);
      mBz = static_cast<float>(grpmag->getNominalL3Field());

      mStraHelper.fitter.setBz(mBz);
    }
    if (!mStraHelper.lut) { /// done only once
      ccdb->setURL(V0BuilderOpts.ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(true);
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mStraHelper.lut = lut;
    }
  }

  template <typename T>
  bool checkEvent(T const& col)
  {
    if (std::abs(col.posZ()) > EventSelection.zvtx.value) {
      return false;
    }
    if (EventSelection.eventSel.value && !col.sel8()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool checkTrack(T const& track, std::string trackName)
  {
    if (track.pt() < TrackSelections.momentum->get(trackName.c_str(), "PtMin")) {
      return false;
    }
    if (track.pt() > TrackSelections.momentum->get(trackName.c_str(), "PtMax")) {
      return false;
    }
    if (std::abs(track.eta()) > TrackSelections.trackProperties->get(trackName.c_str(), "AbsEtaMax")) {
      return false;
    }
    if (track.tpcNClsFound() < TrackSelections.trackProperties->get(trackName.c_str(), "TpcClusterMin")) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < TrackSelections.trackProperties->get(trackName.c_str(), "TpcRowMin")) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < TrackSelections.trackProperties->get(trackName.c_str(), "TpcCrossedOverFoundMin")) {
      return false;
    }
    if (track.tpcNClsShared() > TrackSelections.trackProperties->get(trackName.c_str(), "TpcSharedMax")) {
      return false;
    }
    if (track.tpcFractionSharedCls() > TrackSelections.trackProperties->get(trackName.c_str(), "TpcFracSharedMax")) {
      return false;
    }
    if (track.itsNCls() < TrackSelections.trackProperties->get(trackName.c_str(), "ItsClusterMin")) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < TrackSelections.trackProperties->get(trackName.c_str(), "ItsIbClusterMin")) {
      return false;
    }
    if (std::abs(track.dcaXY()) > TrackSelections.trackProperties->get(trackName.c_str(), "AbsDcaXyMax")) {
      return false;
    }
    if (std::abs(track.dcaZ()) > TrackSelections.trackProperties->get(trackName.c_str(), "AbsDcaZMax")) {
      return false;
    }
    if (track.tpcChi2NCl() > TrackSelections.trackProperties->get(trackName.c_str(), "Chi2TpcMax")) {
      return false;
    }
    if (track.itsChi2NCl() > TrackSelections.trackProperties->get(trackName.c_str(), "Chi2ItsMax")) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool checkTrackPid(T const& track, std::string trackName)
  {
    float momentum = -99;

    if (TrackSelections.momentum->get(trackName.c_str(), "UseInnerParam") < 0) {
      momentum = track.p();
    } else {
      momentum = track.tpcInnerParam();
    }

    float nsigmaITS = -99;
    float nsigmaTPC = -99;
    float nsigmaTPCTOF = -99;

    if (trackName == std::string("Pion")) {
      nsigmaITS = track.itsNSigmaPi();
      nsigmaTPC = track.tpcNSigmaPi();
      nsigmaTPCTOF = std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi());
    } else if (trackName == std::string("Kaon")) {
      nsigmaITS = track.itsNSigmaKa();
      nsigmaTPC = track.tpcNSigmaKa();
      nsigmaTPCTOF = std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa());
    } else if (trackName == std::string("Proton")) {
      nsigmaITS = track.itsNSigmaPr();
      nsigmaTPC = track.tpcNSigmaPr();
      nsigmaTPCTOF = std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr());
    } else if (trackName == std::string("Deuteron")) {
      nsigmaITS = track.itsNSigmaDe();
      nsigmaTPC = track.tpcNSigmaDe();
      nsigmaTPCTOF = std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe());
    } else {
      LOG(fatal) << "Unsupported track type";
    }

    if (momentum < TrackSelections.momentum->get(trackName.c_str(), "PThres")) {
      if (nsigmaITS < TrackSelections.pid->get(trackName.c_str(), "ItsMin") || nsigmaITS > TrackSelections.pid->get(trackName.c_str(), "ItsMax")) {
        return false;
      }
      if (nsigmaTPC < TrackSelections.pid->get(trackName.c_str(), "TpcMin") || nsigmaTPC > TrackSelections.pid->get(trackName.c_str(), "TpcMax")) {
        return false;
      }
    } else {
      if (nsigmaTPCTOF > TrackSelections.pid->get(trackName.c_str(), "TpcTofMax")) {
        return false;
      }
    }
    return true;
  }

  bool checkLambda(float lambdaPt, float lambdaDauDca, float lambdaCpa, float lambdaRadius, float lambdaPos, float kaonMass, float lambdaMass)
  {
    if (lambdaPt < LambdaSelections.ptMin) {
      return false;
    }
    if (lambdaDauDca > LambdaSelections.dcaDaughMax) {
      return false;
    }
    if (lambdaCpa < LambdaSelections.cpaMin) {
      return false;
    }
    if (lambdaRadius < LambdaSelections.tranRadMin) {
      return false;
    }
    if (lambdaRadius > LambdaSelections.tranRadMax) {
      return false;
    }
    if (lambdaPos > LambdaSelections.decVtxMax) {
      return false;
    }
    if (LambdaSelections.rejectKaons) {
      if (kaonMass > LambdaSelections.invKaonMassLow && kaonMass < LambdaSelections.invKaonMassUp) {
        return false;
      }
    }
    if (lambdaMass < LambdaSelections.invMassLow) {
      return false;
    }
    if (lambdaMass > LambdaSelections.invMassUp) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool checkLambdaDaughter(T const& track, float eta, float dca, float nSigmaTPC)
  {
    if (std::abs(eta) > LambdaDaughterSelections.absEtaMax.value) {
      return false;
    }
    if (std::abs(dca) < LambdaDaughterSelections.dcaMin.value) {
      return false;
    }
    if (track.tpcNClsFound() < LambdaDaughterSelections.tpcClusterMin.value) {
      return false;
    }
    if (std::abs(nSigmaTPC) > LambdaDaughterSelections.tpcMax.value) {
      return false;
    }
    return true;
  }

  float getkstar(const ROOT::Math::PtEtaPhiMVector part1,
                 const ROOT::Math::PtEtaPhiMVector part2)
  {
    const ROOT::Math::PtEtaPhiMVector trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax =
      beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay =
      beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    ROOT::Math::PxPyPzMVector partOneCMS(part1);
    ROOT::Math::PxPyPzMVector partTwoCMS(part2);
    const ROOT::Math::Boost boostPRF =
      ROOT::Math::Boost(-betax, -betay, -betaz);
    partOneCMS = boostPRF(partOneCMS);
    partTwoCMS = boostPRF(partTwoCMS);
    const ROOT::Math::PxPyPzMVector trackRelK = partOneCMS - partTwoCMS;
    return 0.5 * trackRelK.P();
  }

  ROOT::Math::PxPyPzEVector
    getqij(const ROOT::Math::PtEtaPhiMVector parti, const ROOT::Math::PtEtaPhiMVector partj)
  {
    ROOT::Math::PxPyPzEVector vecparti(parti);
    ROOT::Math::PxPyPzEVector vecpartj(partj);
    ROOT::Math::PxPyPzEVector trackSum = vecparti + vecpartj;
    ROOT::Math::PxPyPzEVector trackDifference = vecparti - vecpartj;
    float scaling = trackDifference.Dot(trackSum) / trackSum.Dot(trackSum);
    return trackDifference - scaling * trackSum;
  }
  float getQ3(const ROOT::Math::PtEtaPhiMVector part1, const ROOT::Math::PtEtaPhiMVector part2, const ROOT::Math::PtEtaPhiMVector part3)
  {
    ROOT::Math::PxPyPzEVector q12 = getqij(part1, part2);
    ROOT::Math::PxPyPzEVector q23 = getqij(part2, part3);
    ROOT::Math::PxPyPzEVector q31 = getqij(part3, part1);
    float q32 = q12.M2() + q23.M2() + q31.M2();
    return std::sqrt(-q32);
  }

  template <typename T>
  float itsSignal(T const& track)
  {
    uint32_t clsizeflag = track.itsClusterSizes();
    auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
    auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
    auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
    auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
    auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
    auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
    auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;
    int numLayers = 7;
    int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
    float cosLamnda = 1. / std::cosh(track.eta());
    return (static_cast<float>(sumClusterSizes) / numLayers) * cosLamnda;
  };

  void process(cf_trigger::FullCollision const& col, aod::BCs const&, cf_trigger::FullTracks const& tracks, o2::aod::V0s const& v0s)
  {

    auto tracksWithItsPid = soa::Attach<cf_trigger::FullTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);

    // reset all arrays
    keepEventTightLimit.fill(false);
    keepEventLooseLimit.fill(false);
    signalTightLimit.fill(0);
    signalLooseLimit.fill(0);

    registryTriggerQA.fill(HIST("fProcessedEvents"), 0);
    registryParticleQA.fill(HIST("EventQA/Before/fMultiplicity"), col.multNTracksPV());
    registryParticleQA.fill(HIST("EventQA/Before/fZvtx"), col.posZ());

    if (!checkEvent(col)) {
      tags(keepEventTightLimit[cf_trigger::kPPP], keepEventLooseLimit[cf_trigger::kPPP],
           keepEventTightLimit[cf_trigger::kPPL], keepEventLooseLimit[cf_trigger::kPPL],
           keepEventTightLimit[cf_trigger::kPLL], keepEventLooseLimit[cf_trigger::kPLL],
           keepEventTightLimit[cf_trigger::kLLL], keepEventLooseLimit[cf_trigger::kLLL],
           keepEventTightLimit[cf_trigger::kPPPhi], keepEventLooseLimit[cf_trigger::kPPPhi],
           keepEventTightLimit[cf_trigger::kPPRho], keepEventLooseLimit[cf_trigger::kPPRho],
           keepEventTightLimit[cf_trigger::kPD], keepEventLooseLimit[cf_trigger::kPD],
           keepEventTightLimit[cf_trigger::kLD], keepEventLooseLimit[cf_trigger::kLD],
           keepEventTightLimit[cf_trigger::kPhiD], keepEventLooseLimit[cf_trigger::kPhiD],
           keepEventTightLimit[cf_trigger::kRhoD], keepEventLooseLimit[cf_trigger::kRhoD]);
      return;
    }

    registryParticleQA.fill(HIST("EventQA/After/fMultiplicity"), col.multNTracksPV());
    registryParticleQA.fill(HIST("EventQA/After/fZvtx"), col.posZ());

    initCCDB(col.bc().runNumber());

    // clear particle vectors
    vecProton.clear();
    vecAntiProton.clear();
    vecDeuteron.clear();
    vecAntiDeuteron.clear();
    vecLambda.clear();
    vecAntiLambda.clear();
    vecKaon.clear();
    vecAntiKaon.clear();
    vecPhi.clear();
    vecPion.clear();
    vecAntiPion.clear();
    vecRho.clear();
    // clear index vectors for all particles
    idxProton.clear();
    idxAntiProton.clear();
    idxDeuteron.clear();
    idxAntiDeuteron.clear();
    idxKaon.clear();
    idxAntiKaon.clear();
    idxPion.clear();
    idxAntiPion.clear();
    // clear index vectors for daughters
    idxLambdaDaughProton.clear();
    idxLambdaDaughPion.clear();
    idxAntiLambdaDaughProton.clear();
    idxAntiLambdaDaughPion.clear();
    idxPhiDaughPos.clear();
    idxPhiDaughNeg.clear();
    idxRhoDaughPos.clear();
    idxRhoDaughNeg.clear();

    for (auto const& track : tracksWithItsPid) {

      // get paritcles
      if (track.sign() > 0) {
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fPt"), track.pt());
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fEta"), track.eta());
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fPhi"), track.phi());
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fItsSignal"), track.p(), itsSignal(track));
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fTpcSignal"), track.p(), track.tpcSignal());
        registryParticleQA.fill(HIST("TrackQA/Before/Particle/fTofSignal"), track.p(), track.beta());

        registryParticleQA.fill(HIST("TrackQA/Before/Pion/fNsigmaITS"), track.p(), track.itsNSigmaPi());
        registryParticleQA.fill(HIST("TrackQA/Before/Pion/fNsigmaTPC"), track.p(), track.tpcNSigmaPi());
        registryParticleQA.fill(HIST("TrackQA/Before/Pion/fNsigmaTOF"), track.p(), track.tofNSigmaPi());
        registryParticleQA.fill(HIST("TrackQA/Before/Pion/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()));

        registryParticleQA.fill(HIST("TrackQA/Before/Kaon/fNsigmaITS"), track.p(), track.itsNSigmaKa());
        registryParticleQA.fill(HIST("TrackQA/Before/Kaon/fNsigmaTPC"), track.p(), track.tpcNSigmaKa());
        registryParticleQA.fill(HIST("TrackQA/Before/Kaon/fNsigmaTOF"), track.p(), track.tofNSigmaKa());
        registryParticleQA.fill(HIST("TrackQA/Before/Kaon/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()));

        registryParticleQA.fill(HIST("TrackQA/Before/Proton/fNsigmaITS"), track.p(), track.itsNSigmaPr());
        registryParticleQA.fill(HIST("TrackQA/Before/Proton/fNsigmaTPC"), track.p(), track.tpcNSigmaPr());
        registryParticleQA.fill(HIST("TrackQA/Before/Proton/fNsigmaTOF"), track.p(), track.tofNSigmaPr());
        registryParticleQA.fill(HIST("TrackQA/Before/Proton/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));

        registryParticleQA.fill(HIST("TrackQA/Before/Deuteron/fNsigmaITS"), track.p(), track.itsNSigmaDe());
        registryParticleQA.fill(HIST("TrackQA/Before/Deuteron/fNsigmaTPC"), track.p(), track.tpcNSigmaDe());
        registryParticleQA.fill(HIST("TrackQA/Before/Deuteron/fNsigmaTOF"), track.p(), track.tofNSigmaDe());
        registryParticleQA.fill(HIST("TrackQA/Before/Deuteron/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe()));

        if (checkTrack(track, std::string("Pion")) && checkTrackPid(track, std::string("Pion"))) {
          vecPion.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassPionCharged);
          idxPion.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/Pion/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/Pion/fNsigmaIts"), track.p(), track.itsNSigmaPi());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fNsigmaTpc"), track.p(), track.tpcNSigmaPi());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fNsigmaTof"), track.p(), track.tofNSigmaPi());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()));

          registryParticleQA.fill(HIST("TrackQA/After/Pion/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/Pion/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/Pion/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/Pion/fItsChi2OverCluster"), track.itsChi2NCl());
        }

        if (checkTrack(track, std::string("Kaon")) && checkTrackPid(track, std::string("Kaon"))) {
          vecKaon.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassKaonCharged);
          idxKaon.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fNsigmaIts"), track.p(), track.itsNSigmaKa());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fNsigmaTpc"), track.p(), track.tpcNSigmaKa());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fNsigmaTof"), track.p(), track.tofNSigmaKa());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()));

          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/Kaon/fItsChi2OverCluster"), track.itsChi2NCl());
        }

        if (checkTrack(track, std::string("Proton")) && checkTrackPid(track, std::string("Proton"))) {
          vecProton.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassProton);
          idxProton.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/Proton/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/Proton/fNsigmaIts"), track.p(), track.itsNSigmaPr());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fNsigmaTpc"), track.p(), track.tpcNSigmaPr());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fNsigmaTof"), track.p(), track.tofNSigmaPr());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));

          registryParticleQA.fill(HIST("TrackQA/After/Proton/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/Proton/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/Proton/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/Proton/fItsChi2OverCluster"), track.itsChi2NCl());
        }

        if (checkTrack(track, std::string("Deuteron")) && checkTrackPid(track, std::string("Deuteron"))) {
          vecDeuteron.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassDeuteron);
          idxDeuteron.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fNsigmaIts"), track.p(), track.itsNSigmaDe());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fNsigmaTpc"), track.p(), track.tpcNSigmaDe());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fNsigmaTof"), track.p(), track.tofNSigmaDe());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe()));

          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/Deuteron/fItsChi2OverCluster"), track.itsChi2NCl());
        }
      }

      if (track.sign() < 0) {
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fPt"), track.pt());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fEta"), track.eta());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fPhi"), track.phi());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fItsSignal"), track.p(), itsSignal(track));
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fTpcSignal"), track.p(), track.tpcSignal());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiParticle/fTofSignal"), track.p(), track.beta());

        registryParticleQA.fill(HIST("TrackQA/Before/AntiPion/fNsigmaITS"), track.p(), track.itsNSigmaPi());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiPion/fNsigmaTPC"), track.p(), track.tpcNSigmaPi());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiPion/fNsigmaTOF"), track.p(), track.tofNSigmaPi());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiPion/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()));

        registryParticleQA.fill(HIST("TrackQA/Before/AntiKaon/fNsigmaITS"), track.p(), track.itsNSigmaKa());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiKaon/fNsigmaTPC"), track.p(), track.tpcNSigmaKa());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiKaon/fNsigmaTOF"), track.p(), track.tofNSigmaKa());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiKaon/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()));

        registryParticleQA.fill(HIST("TrackQA/Before/AntiProton/fNsigmaITS"), track.p(), track.itsNSigmaPr());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiProton/fNsigmaTPC"), track.p(), track.tpcNSigmaPr());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiProton/fNsigmaTOF"), track.p(), track.tofNSigmaPr());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiProton/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));

        registryParticleQA.fill(HIST("TrackQA/Before/AntiDeuteron/fNsigmaITS"), track.p(), track.itsNSigmaDe());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiDeuteron/fNsigmaTPC"), track.p(), track.tpcNSigmaDe());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiDeuteron/fNsigmaTOF"), track.p(), track.tofNSigmaDe());
        registryParticleQA.fill(HIST("TrackQA/Before/AntiDeuteron/fNsigmaTPCTOF"), track.p(), std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe()));

        if (checkTrack(track, std::string("Pion")) && checkTrackPid(track, std::string("Pion"))) {
          vecAntiPion.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassPionCharged);
          idxAntiPion.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fNsigmaIts"), track.p(), track.itsNSigmaPi());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fNsigmaTpc"), track.p(), track.tpcNSigmaPi());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fNsigmaTof"), track.p(), track.tofNSigmaPi());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()));

          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/AntiPion/fItsChi2OverCluster"), track.itsChi2NCl());
        }

        if (checkTrack(track, std::string("Kaon")) && checkTrackPid(track, std::string("Kaon"))) {
          vecAntiKaon.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassKaonCharged);
          idxAntiKaon.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fNsigmaIts"), track.p(), track.itsNSigmaKa());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fNsigmaTpc"), track.p(), track.tpcNSigmaKa());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fNsigmaTof"), track.p(), track.tofNSigmaKa());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()));

          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/AntiKaon/fItsChi2OverCluster"), track.itsChi2NCl());
        }

        if (checkTrack(track, std::string("Proton")) && checkTrackPid(track, std::string("Proton"))) {
          vecAntiProton.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassProton);
          idxAntiProton.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fNsigmaIts"), track.p(), track.itsNSigmaPr());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fNsigmaTpc"), track.p(), track.tpcNSigmaPr());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fNsigmaTof"), track.p(), track.tofNSigmaPr());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));

          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/AntiProton/fItsChi2OverCluster"), track.itsChi2NCl());
        }

        if (checkTrack(track, std::string("Deuteron")) && checkTrackPid(track, std::string("Deuteron"))) {
          vecAntiDeuteron.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassDeuteron);
          idxAntiDeuteron.push_back(track.globalIndex());

          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fPt"), track.pt());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fPTpc"), track.tpcInnerParam());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fMomCor"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fEta"), track.eta());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fPhi"), track.phi());

          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fNsigmaIts"), track.p(), track.itsNSigmaDe());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fNsigmaTpc"), track.p(), track.tpcNSigmaDe());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fNsigmaTof"), track.p(), track.tofNSigmaDe());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fNsigmaTpcTof"), track.p(), std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe()));

          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fItsSignal"), track.p(), itsSignal(track));
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcSignal"), track.p(), track.tpcSignal());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTofBeta"), track.p(), track.beta());

          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fDcaXy"), track.pt(), track.dcaXY());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fDcaZ"), track.pt(), track.dcaZ());

          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcClusters"), track.tpcNClsFound());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcCrossedRows"), track.tpcNClsCrossedRows());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcSharedClusters"), track.tpcNClsShared());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcSharedClusterOverClusterss"), track.tpcFractionSharedCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcFindableOverRows"), track.tpcCrossedRowsOverFindableCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fTpcChi2OverCluster"), track.tpcChi2NCl());

          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fItsClusters"), track.itsNCls());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fItsIbClusters"), track.itsNClsInnerBarrel());
          registryParticleQA.fill(HIST("TrackQA/After/AntiDeuteron/fItsChi2OverCluster"), track.itsChi2NCl());
        }
      }
    }

    // loop over and build v0s
    for (auto const& v0 : v0s) {

      auto posTrack = v0.posTrack_as<cf_trigger::FullTracks>();
      auto negTrack = v0.negTrack_as<cf_trigger::FullTracks>();

      auto posTrackPar = getTrackParCov(posTrack);
      auto negTrackPar = getTrackParCov(negTrack);

      if (!mStraHelper.buildV0Candidate(v0.collisionId(), col.posX(), col.posY(), col.posZ(), posTrack, negTrack, posTrackPar, negTrackPar, false, false, false)) {
        continue;
      }

      float lambdaPt = std::hypot(mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1]);
      float lambdaPos = std::hypot(mStraHelper.v0.position[0] - col.posX(), mStraHelper.v0.position[1] - col.posY(), mStraHelper.v0.position[2] - col.posZ());
      float lambdaRadius = std::hypot(mStraHelper.v0.position[0], mStraHelper.v0.position[1]);
      float lambdaEta = RecoDecay::eta(std::array{mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1], mStraHelper.v0.momentum[2]});
      float lambdaPhi = RecoDecay::phi(mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1]);
      float lambdaCpa = std::cos(mStraHelper.v0.pointingAngle);
      float lambdaDauDca = mStraHelper.v0.daughterDCA;
      float lambdaMass = mStraHelper.v0.massLambda;
      float antiLambdaMass = mStraHelper.v0.massAntiLambda;
      float kaonMass = mStraHelper.v0.massK0Short;

      float posTrackEta = RecoDecay::eta(std::array{mStraHelper.v0.positiveMomentum[0], mStraHelper.v0.positiveMomentum[1], mStraHelper.v0.positiveMomentum[2]});
      float posTrackDca = mStraHelper.v0.positiveDCAxy;
      float negTrackEta = RecoDecay::eta(std::array{mStraHelper.v0.negativeMomentum[0], mStraHelper.v0.negativeMomentum[1], mStraHelper.v0.negativeMomentum[2]});
      float negTrackDca = mStraHelper.v0.negativeDCAxy;

      registryParticleQA.fill(HIST("LambdaQA/Before/fPt"), lambdaPt);
      registryParticleQA.fill(HIST("LambdaQA/Before/fEta"), lambdaEta);
      registryParticleQA.fill(HIST("LambdaQA/Before/fPhi"), lambdaPhi);
      registryParticleQA.fill(HIST("LambdaQA/Before/fInvMassLambda"), lambdaMass);
      registryParticleQA.fill(HIST("LambdaQA/Before/fInvMassAntiLambda"), antiLambdaMass);
      registryParticleQA.fill(HIST("LambdaQA/Before/fInvMassLambdaVsAntiLambda"), lambdaMass, antiLambdaMass);
      registryParticleQA.fill(HIST("LambdaQA/Before/fInvMassLambdaVsKaon"), lambdaMass, kaonMass);
      registryParticleQA.fill(HIST("LambdaQA/Before/fInvMassAntiLambdaVsKaon"), antiLambdaMass, kaonMass);
      registryParticleQA.fill(HIST("LambdaQA/Before/fDcaDaugh"), lambdaDauDca);
      registryParticleQA.fill(HIST("LambdaQA/Before/fCpa"), lambdaCpa);
      registryParticleQA.fill(HIST("LambdaQA/Before/fTranRad"), lambdaRadius);
      registryParticleQA.fill(HIST("LambdaQA/Before/fDecVtx"), lambdaPos);

      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fPt"), posTrack.pt());
      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fEta"), posTrackEta);
      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fPhi"), posTrack.phi());
      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fDcaXy"), posTrack.pt(), posTrackDca);
      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fTpcClusters"), posTrack.tpcNClsFound());
      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fNsigmaTpcProton"), posTrack.p(), posTrack.tpcNSigmaPr());
      registryParticleQA.fill(HIST("LambdaQA/Before/PosDaughter/fNsigmaTpcPion"), posTrack.p(), posTrack.tpcNSigmaPi());

      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fPt"), negTrack.pt());
      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fEta"), negTrackEta);
      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fPhi"), negTrack.phi());
      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fDcaXy"), negTrack.pt(), negTrackDca);
      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fTpcClusters"), negTrack.tpcNClsFound());
      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fNsigmaTpcProton"), negTrack.p(), negTrack.tpcNSigmaPr());
      registryParticleQA.fill(HIST("LambdaQA/Before/NegDaughter/fNsigmaTpcPion"), negTrack.p(), negTrack.tpcNSigmaPi());

      if (checkLambda(lambdaPt, lambdaDauDca, lambdaCpa, lambdaRadius, lambdaPos, kaonMass, lambdaMass) && checkLambdaDaughter(posTrack, posTrackEta, posTrackDca, posTrack.tpcNSigmaPr()) && checkLambdaDaughter(negTrack, negTrackEta, negTrackDca, negTrack.tpcNSigmaPi())) {
        vecLambda.emplace_back(lambdaPt, lambdaEta, lambdaPhi, o2::constants::physics::MassLambda0);
        idxLambdaDaughProton.push_back(posTrack.globalIndex());
        idxLambdaDaughPion.push_back(negTrack.globalIndex());

        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fPt"), lambdaPt);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fEta"), lambdaEta);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fPhi"), lambdaPhi);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fInvMass"), lambdaMass);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fInvMassLambdaVsAntiLambda"), lambdaMass, antiLambdaMass);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fInvMassLambdaVsKaon"), lambdaMass, kaonMass);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fDcaDaugh"), lambdaDauDca);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fCpa"), lambdaCpa);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fTranRad"), lambdaRadius);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/fDecVtx"), lambdaPos);

        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/PosDaughter/fPt"), posTrack.pt());
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/PosDaughter/fEta"), posTrackEta);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/PosDaughter/fPhi"), posTrack.phi());
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/PosDaughter/fDcaXy"), posTrack.pt(), posTrackDca);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/PosDaughter/fTpcClusters"), posTrack.tpcNClsFound());
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/PosDaughter/fNsigmaTpc"), posTrack.p(), posTrack.tpcNSigmaPr());

        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/NegDaughter/fPt"), negTrack.pt());
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/NegDaughter/fEta"), negTrackEta);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/NegDaughter/fPhi"), negTrack.phi());
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/NegDaughter/fDcaXy"), negTrack.pt(), negTrackDca);
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/NegDaughter/fTpcClusters"), negTrack.tpcNClsFound());
        registryParticleQA.fill(HIST("LambdaQA/After/Lambda/NegDaughter/fNsigmaTpc"), negTrack.p(), negTrack.tpcNSigmaPi());
      }

      if (checkLambda(lambdaPt, lambdaDauDca, lambdaCpa, lambdaRadius, lambdaPos, kaonMass, antiLambdaMass) && checkLambdaDaughter(posTrack, posTrackEta, posTrackDca, posTrack.tpcNSigmaPi()) && checkLambdaDaughter(negTrack, negTrackEta, negTrackDca, negTrack.tpcNSigmaPr())) {
        vecAntiLambda.emplace_back(lambdaPt, lambdaEta, lambdaPhi, o2::constants::physics::MassLambda0);

        idxAntiLambdaDaughProton.push_back(negTrack.globalIndex());
        idxAntiLambdaDaughPion.push_back(posTrack.globalIndex());

        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fPt"), lambdaPt);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fEta"), lambdaEta);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fPhi"), lambdaPhi);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fInvMass"), antiLambdaMass);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fInvMassAntiLambdaVsLambda"), antiLambdaMass, lambdaMass);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fInvMassAntiLambdaVsKaon"), antiLambdaMass, kaonMass);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fDcaDaugh"), lambdaDauDca);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fCpa"), lambdaCpa);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fTranRad"), lambdaRadius);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/fDecVtx"), lambdaPos);

        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/PosDaughter/fPt"), posTrack.pt());
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/PosDaughter/fEta"), posTrackEta);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/PosDaughter/fPhi"), posTrack.phi());
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/PosDaughter/fDcaXy"), posTrack.pt(), posTrackDca);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/PosDaughter/fTpcClusters"), posTrack.tpcNClsFound());
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/PosDaughter/fNsigmaTpc"), posTrack.p(), posTrack.tpcNSigmaPr());

        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/NegDaughter/fPt"), negTrack.pt());
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/NegDaughter/fEta"), negTrackEta);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/NegDaughter/fPhi"), negTrack.phi());
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/NegDaughter/fDcaXy"), negTrack.pt(), negTrackDca);
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/NegDaughter/fTpcClusters"), negTrack.tpcNClsFound());
        registryParticleQA.fill(HIST("LambdaQA/After/AntiLambda/NegDaughter/fNsigmaTpc"), negTrack.p(), negTrack.tpcNSigmaPi());
      }
    }

    // build phi candidates
    for (size_t k1 = 0; k1 < vecKaon.size(); k1++) {
      for (size_t k2 = 0; k2 < vecAntiKaon.size(); k2++) {
        ROOT::Math::PtEtaPhiMVector phi = vecKaon.at(k1) + vecAntiKaon.at(k2);

        registryParticleQA.fill(HIST("PhiQA/Before/fInvMass"), phi.M());
        registryParticleQA.fill(HIST("PhiQA/Before/fPt"), phi.Pt());
        registryParticleQA.fill(HIST("PhiQA/Before/fEta"), phi.Eta());
        registryParticleQA.fill(HIST("PhiQA/Before/fPhi"), RecoDecay::constrainAngle(phi.Phi()));

        if ((phi.M() >= PhiSelections.invMassLow.value) &&
            (phi.M() <= PhiSelections.invMassUp.value)) {
          vecPhi.push_back(phi);
          idxPhiDaughPos.push_back(idxKaon.at(k1));
          idxPhiDaughNeg.push_back(idxAntiKaon.at(k2));

          registryParticleQA.fill(HIST("PhiQA/After/fInvMass"), phi.M());
          registryParticleQA.fill(HIST("PhiQA/After/fPt"), phi.Pt());
          registryParticleQA.fill(HIST("PhiQA/After/fEta"), phi.Eta());
          registryParticleQA.fill(HIST("PhiQA/After/fPhi"), RecoDecay::constrainAngle(phi.Phi()));
        }
      }
    }

    // build rho candidates
    for (size_t p1 = 0; p1 < vecPion.size(); p1++) {
      for (size_t p2 = 0; p2 < vecAntiPion.size(); p2++) {
        ROOT::Math::PtEtaPhiMVector rho = vecPion.at(p1) + vecAntiPion.at(p2);

        registryParticleQA.fill(HIST("RhoQA/Before/fInvMass"), rho.M());
        registryParticleQA.fill(HIST("RhoQA/Before/fPt"), rho.Pt());
        registryParticleQA.fill(HIST("RhoQA/Before/fEta"), rho.Eta());
        registryParticleQA.fill(HIST("RhoQA/Before/fPhi"), RecoDecay::constrainAngle(rho.Phi()));

        if (((rho.M() >= RhoSelections.invMassLow.value) && (rho.M() <= RhoSelections.invMassUp.value)) && (rho.Pt() >= RhoSelections.ptLow)) {
          vecRho.push_back(rho);
          idxRhoDaughPos.push_back(idxPion.at(p1));
          idxRhoDaughNeg.push_back(idxAntiPion.at(p2));

          registryParticleQA.fill(HIST("RhoQA/After/fInvMass"), rho.M());
          registryParticleQA.fill(HIST("RhoQA/After/fPt"), rho.Pt());
          registryParticleQA.fill(HIST("RhoQA/After/fEta"), rho.Eta());
          registryParticleQA.fill(HIST("RhoQA/After/fPhi"), RecoDecay::constrainAngle(rho.Phi()));
        }
      }
    }

    float q3 = 999.f, kstar = 999.f;

    // PPP
    if (TriggerSelections.filterSwitches->get("Switch", "PPP") > 0) {
      for (size_t p1 = 0; p1 < vecProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecProton.size(); p2++) {
          for (size_t p3 = p2 + 1; p3 < vecProton.size(); p3++) {
            q3 = getQ3(vecProton.at(p1), vecProton.at(p2), vecProton.at(p3));
            registryTriggerQA.fill(HIST("PPP/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPP/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPP/all/fSE_particle"), q3);
            registryTriggerQA.fill(HIST("PPP/all/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPP/all/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPP/all/fProtonQ3VsPt"), q3, vecProton.at(p3).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPP")) {
              signalLooseLimit[cf_trigger::kPPP] += 1;
              registryTriggerQA.fill(HIST("PPP/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPP/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPP/loose/fSE_particle"), q3);
              registryTriggerQA.fill(HIST("PPP/loose/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPP/loose/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPP/loose/fProtonQ3VsPt"), q3, vecProton.at(p3).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPP")) {
                signalTightLimit[cf_trigger::kPPP] += 1;
                registryTriggerQA.fill(HIST("PPP/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPP/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPP/tight/fSE_particle"), q3);
                registryTriggerQA.fill(HIST("PPP/tight/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPP/tight/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPP/tight/fProtonQ3VsPt"), q3, vecProton.at(p3).Pt());
              }
            }
          }
        }
      }
      for (size_t p1 = 0; p1 < vecAntiProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecAntiProton.size(); p2++) {
          for (size_t p3 = p2 + 1; p3 < vecAntiProton.size(); p3++) {
            q3 = getQ3(vecAntiProton.at(p1), vecAntiProton.at(p2), vecAntiProton.at(p3));
            registryTriggerQA.fill(HIST("PPP/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPP/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPP/all/fSE_antiparticle"), q3);
            registryTriggerQA.fill(HIST("PPP/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPP/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPP/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p3).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPP")) {
              signalLooseLimit[cf_trigger::kPPP] += 1;
              registryTriggerQA.fill(HIST("PPP/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPP/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPP/loose/fSE_antiparticle"), q3);
              registryTriggerQA.fill(HIST("PPP/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPP/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPP/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p3).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPP")) {
                signalTightLimit[cf_trigger::kPPP] += 1;
                registryTriggerQA.fill(HIST("PPP/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPP/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPP/tight/fSE_antiparticle"), q3);
                registryTriggerQA.fill(HIST("PPP/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPP/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPP/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p3).Pt());
              }
            }
          }
        }
      }
    }
    // PPL
    if (TriggerSelections.filterSwitches->get("Switch", "PPL") > 0) {
      for (size_t p1 = 0; p1 < vecProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecProton.size(); p2++) {
          for (size_t l1 = 0; l1 < vecLambda.size(); l1++) {
            if (idxProton.at(p1) == idxLambdaDaughProton.at(l1) || idxProton.at(p2) == idxLambdaDaughProton.at(l1)) {
              continue;
            }
            q3 = getQ3(vecProton.at(p1), vecProton.at(p2), vecLambda.at(l1));
            registryTriggerQA.fill(HIST("PPL/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPL/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPL/all/fSE_particle"), q3);
            registryTriggerQA.fill(HIST("PPL/all/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPL/all/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPL/all/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPL")) {
              signalLooseLimit[cf_trigger::kPPL] += 1;
              registryTriggerQA.fill(HIST("PPL/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPL/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPL/loose/fSE_particle"), q3);
              registryTriggerQA.fill(HIST("PPL/loose/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPL/loose/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPL/loose/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPL")) {
                signalTightLimit[cf_trigger::kPPL] += 1;
                registryTriggerQA.fill(HIST("PPL/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPL/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPL/tight/fSE_particle"), q3);
                registryTriggerQA.fill(HIST("PPL/tight/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPL/tight/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPL/tight/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
              }
            }
          }
        }
      }
      for (size_t p1 = 0; p1 < vecAntiProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecAntiProton.size(); p2++) {
          for (size_t l1 = 0; l1 < vecAntiLambda.size(); l1++) {
            if (idxAntiProton.at(p1) == idxAntiLambdaDaughProton.at(l1) || idxAntiProton.at(p2) == idxAntiLambdaDaughProton.at(l1)) {
              continue;
            }
            q3 = getQ3(vecAntiProton.at(p1), vecAntiProton.at(p2), vecAntiLambda.at(l1));
            registryTriggerQA.fill(HIST("PPL/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPL/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPL/all/fSE_antiparticle"), q3);
            registryTriggerQA.fill(HIST("PPL/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPL/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPL/all/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPL")) {
              signalLooseLimit[cf_trigger::kPPL] += 1;
              registryTriggerQA.fill(HIST("PPL/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPL/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPL/loose/fSE_antiparticle"), q3);
              registryTriggerQA.fill(HIST("PPL/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPL/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPL/loose/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPL")) {
                signalTightLimit[cf_trigger::kPPL] += 1;
                registryTriggerQA.fill(HIST("PPL/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPL/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPL/tight/fSE_antiparticle"), q3);
                registryTriggerQA.fill(HIST("PPL/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPL/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPL/tight/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
              }
            }
          }
        }
      }
    }
    // PLL
    if (TriggerSelections.filterSwitches->get("Switch", "PLL") > 0) {
      for (size_t l1 = 0; l1 < vecLambda.size(); l1++) {
        for (size_t l2 = l1 + 1; l2 < vecLambda.size(); l2++) {
          for (size_t p1 = 0; p1 < vecProton.size(); p1++) {
            if (idxProton.at(p1) == idxLambdaDaughProton.at(l1) || idxProton.at(p1) == idxLambdaDaughProton.at(l2)) {
              continue;
            }
            if (idxLambdaDaughProton.at(l1) == idxLambdaDaughProton.at(l2) || idxLambdaDaughPion.at(l1) == idxLambdaDaughPion.at(l2)) {
              continue;
            }
            q3 = getQ3(vecLambda.at(l1), vecLambda.at(l2), vecProton.at(p1));
            registryTriggerQA.fill(HIST("PLL/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PLL/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PLL/all/fSE_particle"), q3);
            registryTriggerQA.fill(HIST("PLL/all/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
            registryTriggerQA.fill(HIST("PLL/all/fLambdaQ3VsPt"), q3, vecLambda.at(l2).Pt());
            registryTriggerQA.fill(HIST("PLL/all/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PLL")) {
              signalLooseLimit[cf_trigger::kPLL] += 1;
              registryTriggerQA.fill(HIST("PLL/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PLL/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PLL/loose/fSE_particle"), q3);
              registryTriggerQA.fill(HIST("PLL/loose/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
              registryTriggerQA.fill(HIST("PLL/loose/fLambdaQ3VsPt"), q3, vecLambda.at(l2).Pt());
              registryTriggerQA.fill(HIST("PLL/loose/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PLL")) {
                signalTightLimit[cf_trigger::kPLL] += 1;
                registryTriggerQA.fill(HIST("PLL/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PLL/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PLL/tight/fSE_particle"), q3);
                registryTriggerQA.fill(HIST("PLL/tight/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
                registryTriggerQA.fill(HIST("PLL/tight/fLambdaQ3VsPt"), q3, vecLambda.at(l2).Pt());
                registryTriggerQA.fill(HIST("PLL/tight/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
              }
            }
          }
        }
      }
      for (size_t l1 = 0; l1 < vecAntiLambda.size(); l1++) {
        for (size_t l2 = l1 + 1; l2 < vecAntiLambda.size(); l2++) {
          for (size_t p1 = 0; p1 < vecAntiProton.size(); p1++) {
            if (idxAntiProton.at(p1) == idxAntiLambdaDaughProton.at(l1) || idxAntiProton.at(p1) == idxAntiLambdaDaughProton.at(l2)) {
              continue;
            }
            if (idxAntiLambdaDaughProton.at(l1) == idxAntiLambdaDaughProton.at(l2) || idxAntiLambdaDaughPion.at(l1) == idxAntiLambdaDaughPion.at(l2)) {
              continue;
            }
            q3 = getQ3(vecAntiLambda.at(l1), vecAntiLambda.at(l2), vecAntiProton.at(p1));
            registryTriggerQA.fill(HIST("PLL/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PLL/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PLL/all/fSE_antiparticle"), q3);
            registryTriggerQA.fill(HIST("PLL/all/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
            registryTriggerQA.fill(HIST("PLL/all/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l2).Pt());
            registryTriggerQA.fill(HIST("PLL/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PLL")) {
              signalLooseLimit[cf_trigger::kPLL] += 1;
              registryTriggerQA.fill(HIST("PLL/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PLL/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PLL/loose/fSE_antiparticle"), q3);
              registryTriggerQA.fill(HIST("PLL/loose/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
              registryTriggerQA.fill(HIST("PLL/loose/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l2).Pt());
              registryTriggerQA.fill(HIST("PLL/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PLL")) {
                signalTightLimit[cf_trigger::kPLL] += 1;
                registryTriggerQA.fill(HIST("PLL/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PLL/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PLL/tight/fSE_antiparticle"), q3);
                registryTriggerQA.fill(HIST("PLL/tight/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
                registryTriggerQA.fill(HIST("PLL/tight/fAntiLambdaQ3VsPt"), q3, vecAntiLambda.at(l2).Pt());
                registryTriggerQA.fill(HIST("PLL/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
              }
            }
          }
        }
      }
    }
    // LLL
    if (TriggerSelections.filterSwitches->get("Switch", "LLL") > 0) {
      for (size_t l1 = 0; l1 < vecLambda.size(); l1++) {
        for (size_t l2 = l1 + 1; l2 < vecLambda.size(); l2++) {
          for (size_t l3 = l2 + 1; l3 < vecLambda.size(); l3++) {
            if (idxLambdaDaughProton.at(l1) == idxLambdaDaughProton.at(l2) || idxLambdaDaughPion.at(l1) == idxLambdaDaughPion.at(l2) ||
                idxLambdaDaughProton.at(l2) == idxLambdaDaughProton.at(l3) || idxLambdaDaughPion.at(l2) == idxLambdaDaughPion.at(l3) ||
                idxLambdaDaughProton.at(l3) == idxLambdaDaughProton.at(l1) || idxLambdaDaughPion.at(l3) == idxLambdaDaughPion.at(l1)) {
              continue;
            }
            q3 = getQ3(vecLambda.at(l1), vecLambda.at(l2), vecLambda.at(l3));
            registryTriggerQA.fill(HIST("LLL/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("LLL/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("LLL/all/fSE_particle"), q3);
            registryTriggerQA.fill(HIST("LLL/all/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
            registryTriggerQA.fill(HIST("LLL/all/fLambdaQ3VsPt"), q3, vecLambda.at(l2).Pt());
            registryTriggerQA.fill(HIST("LLL/all/fLambdaQ3VsPt"), q3, vecLambda.at(l3).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "LLL")) {
              signalLooseLimit[cf_trigger::kLLL] += 1;
              registryTriggerQA.fill(HIST("LLL/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("LLL/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("LLL/loose/fSE_particle"), q3);
              registryTriggerQA.fill(HIST("LLL/loose/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
              registryTriggerQA.fill(HIST("LLL/loose/fLambdaQ3VsPt"), q3, vecLambda.at(l2).Pt());
              registryTriggerQA.fill(HIST("LLL/loose/fLambdaQ3VsPt"), q3, vecLambda.at(l3).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "LLL")) {
                signalTightLimit[cf_trigger::kLLL] += 1;
                registryTriggerQA.fill(HIST("LLL/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("LLL/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("LLL/tight/fSE_particle"), q3);
                registryTriggerQA.fill(HIST("LLL/tight/fLambdaQ3VsPt"), q3, vecLambda.at(l1).Pt());
                registryTriggerQA.fill(HIST("LLL/tight/fLambdaQ3VsPt"), q3, vecLambda.at(l2).Pt());
                registryTriggerQA.fill(HIST("LLL/tight/fLambdaQ3VsPt"), q3, vecLambda.at(l3).Pt());
              }
            }
          }
        }
      }
      for (size_t l1 = 0; l1 < vecAntiLambda.size(); l1++) {
        for (size_t l2 = l1 + 1; l2 < vecAntiLambda.size(); l2++) {
          for (size_t l3 = l2 + 1; l3 < vecAntiLambda.size(); l3++) {
            if (idxAntiLambdaDaughProton.at(l1) == idxAntiLambdaDaughProton.at(l2) || idxAntiLambdaDaughPion.at(l1) == idxAntiLambdaDaughPion.at(l2) ||
                idxAntiLambdaDaughProton.at(l2) == idxAntiLambdaDaughProton.at(l3) || idxAntiLambdaDaughPion.at(l2) == idxAntiLambdaDaughPion.at(l3) ||
                idxAntiLambdaDaughProton.at(l3) == idxAntiLambdaDaughProton.at(l1) || idxAntiLambdaDaughPion.at(l3) == idxAntiLambdaDaughPion.at(l1)) {
              continue;
            }
            q3 = getQ3(vecAntiLambda.at(l1), vecAntiLambda.at(l2), vecAntiLambda.at(l3));
            registryTriggerQA.fill(HIST("LLL/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("LLL/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("LLL/all/fSE_antiparticle"), q3);
            registryTriggerQA.fill(HIST("LLL/all/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
            registryTriggerQA.fill(HIST("LLL/all/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l2).Pt());
            registryTriggerQA.fill(HIST("LLL/all/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l3).Pt());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "LLL")) {
              signalLooseLimit[cf_trigger::kLLL] += 1;
              registryTriggerQA.fill(HIST("LLL/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("LLL/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("LLL/loose/fSE_antiparticle"), q3);
              registryTriggerQA.fill(HIST("LLL/loose/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
              registryTriggerQA.fill(HIST("LLL/loose/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l2).Pt());
              registryTriggerQA.fill(HIST("LLL/loose/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l3).Pt());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "LLL")) {
                signalTightLimit[cf_trigger::kLLL] += 1;
                registryTriggerQA.fill(HIST("LLL/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("LLL/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("LLL/tight/fSE_antiparticle"), q3);
                registryTriggerQA.fill(HIST("LLL/tight/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l1).Pt());
                registryTriggerQA.fill(HIST("LLL/tight/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l2).Pt());
                registryTriggerQA.fill(HIST("LLL/tight/fLambdaQ3VsPt"), q3, vecAntiLambda.at(l3).Pt());
              }
            }
          }
        }
      }
    }
    // PPPhi
    if (TriggerSelections.filterSwitches->get("Switch", "PPPhi") > 0) {
      for (size_t p1 = 0; p1 < vecProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecProton.size(); p2++) {
          for (size_t phi1 = 0; phi1 < vecPhi.size(); phi1++) {
            if (idxProton.at(p1) == idxPhiDaughPos.at(phi1) || idxProton.at(p2) == idxPhiDaughPos.at(phi1)) {
              continue;
            }
            q3 = getQ3(vecProton.at(p1), vecProton.at(p2), vecPhi.at(phi1));
            registryTriggerQA.fill(HIST("PPPhi/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPPhi/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPPhi/all/fSE_particle"), q3);
            registryTriggerQA.fill(HIST("PPPhi/all/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPPhi/all/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPPhi/all/fPhiQ3VsPt"), q3, vecPhi.at(phi1).Pt());
            registryTriggerQA.fill(HIST("PPPhi/all/fPhiQ3VsInvMass"), q3, vecPhi.at(phi1).M());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPPhi")) {
              signalLooseLimit[cf_trigger::kPPPhi] += 1;
              registryTriggerQA.fill(HIST("PPPhi/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPPhi/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPPhi/loose/fSE_particle"), q3);
              registryTriggerQA.fill(HIST("PPPhi/loose/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPPhi/loose/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPPhi/loose/fPhiQ3VsPt"), q3, vecPhi.at(phi1).Pt());
              registryTriggerQA.fill(HIST("PPPhi/loose/fPhiQ3VsInvMass"), q3, vecPhi.at(phi1).M());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPPhi") &&
                  vecPhi.at(phi1).M() > PhiSelections.tightInvMassLow.value && vecPhi.at(phi1).M() < PhiSelections.tightInvMassUp.value) {
                signalTightLimit[cf_trigger::kPPPhi] += 1;
                registryTriggerQA.fill(HIST("PPPhi/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPPhi/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPPhi/tight/fSE_particle"), q3);
                registryTriggerQA.fill(HIST("PPPhi/tight/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPPhi/tight/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPPhi/tight/fPhiQ3VsPt"), q3, vecPhi.at(phi1).Pt());
                registryTriggerQA.fill(HIST("PPPhi/tight/fPhiQ3VsInvMass"), q3, vecPhi.at(phi1).M());
              }
            }
          }
        }
      }
      for (size_t p1 = 0; p1 < vecAntiProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecAntiProton.size(); p2++) {
          for (size_t phi1 = 0; phi1 < vecPhi.size(); phi1++) {
            if (idxAntiProton.at(p1) == idxPhiDaughNeg.at(phi1) || idxAntiProton.at(p2) == idxPhiDaughNeg.at(phi1)) {
              continue;
            }
            q3 = getQ3(vecAntiProton.at(p1), vecAntiProton.at(p2), vecPhi.at(phi1));
            registryTriggerQA.fill(HIST("PPPhi/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPPhi/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPPhi/all/fSE_antiparticle"), q3);
            registryTriggerQA.fill(HIST("PPPhi/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPPhi/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPPhi/all/fPhiQ3VsPt"), q3, vecPhi.at(phi1).Pt());
            registryTriggerQA.fill(HIST("PPPhi/all/fPhiQ3VsInvMass"), q3, vecPhi.at(phi1).M());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPPhi")) {
              signalLooseLimit[cf_trigger::kPPPhi] += 1;
              registryTriggerQA.fill(HIST("PPPhi/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPPhi/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPPhi/loose/fSE_antiparticle"), q3);
              registryTriggerQA.fill(HIST("PPPhi/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPPhi/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPPhi/loose/fPhiQ3VsPt"), q3, vecPhi.at(phi1).Pt());
              registryTriggerQA.fill(HIST("PPPhi/loose/fPhiQ3VsInvMass"), q3, vecPhi.at(phi1).M());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPPhi") &&
                  vecPhi.at(phi1).M() > PhiSelections.tightInvMassLow.value && vecPhi.at(phi1).M() < PhiSelections.tightInvMassUp.value) {
                signalTightLimit[cf_trigger::kPPPhi] += 1;
                registryTriggerQA.fill(HIST("PPPhi/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPPhi/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPPhi/tight/fSE_antiparticle"), q3);
                registryTriggerQA.fill(HIST("PPPhi/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPPhi/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPPhi/tight/fPhiQ3VsPt"), q3, vecPhi.at(phi1).Pt());
                registryTriggerQA.fill(HIST("PPPhi/tight/fPhiQ3VsInvMass"), q3, vecPhi.at(phi1).M());
              }
            }
          }
        }
      }
    }
    // PPRho
    if (TriggerSelections.filterSwitches->get("Switch", "PPRho") > 0) {
      for (size_t p1 = 0; p1 < vecProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecProton.size(); p2++) {
          for (size_t r1 = 0; r1 < vecRho.size(); r1++) {
            if (idxProton.at(p1) == idxRhoDaughPos.at(r1) || idxProton.at(p2) == idxRhoDaughPos.at(r1)) {
              continue;
            }
            q3 = getQ3(vecProton.at(p1), vecProton.at(p2), vecRho.at(r1));
            registryTriggerQA.fill(HIST("PPRho/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPRho/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPRho/all/fSE_particle"), q3);
            registryTriggerQA.fill(HIST("PPRho/all/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPRho/all/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPRho/all/fRhoQ3VsPt"), q3, vecRho.at(r1).Pt());
            registryTriggerQA.fill(HIST("PPRho/all/fRhoQ3VsInvMass"), q3, vecRho.at(r1).M());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPRho")) {
              signalLooseLimit[cf_trigger::kPPRho] += 1;
              registryTriggerQA.fill(HIST("PPRho/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPRho/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPRho/loose/fSE_particle"), q3);
              registryTriggerQA.fill(HIST("PPRho/loose/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPRho/loose/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPRho/loose/fRhoQ3VsPt"), q3, vecRho.at(r1).Pt());
              registryTriggerQA.fill(HIST("PPRho/loose/fRhoQ3VsInvMass"), q3, vecRho.at(r1).M());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPRho") &&
                  vecRho.at(r1).M() > RhoSelections.tightInvMassLow.value && vecRho.at(r1).M() < RhoSelections.tightInvMassUp.value) {
                signalTightLimit[cf_trigger::kPPRho] += 1;
                registryTriggerQA.fill(HIST("PPRho/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPRho/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPRho/tight/fSE_particle"), q3);
                registryTriggerQA.fill(HIST("PPRho/tight/fProtonQ3VsPt"), q3, vecProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPRho/tight/fProtonQ3VsPt"), q3, vecProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPRho/tight/fRhoQ3VsPt"), q3, vecRho.at(r1).Pt());
                registryTriggerQA.fill(HIST("PPRho/tight/fRhoQ3VsInvMass"), q3, vecRho.at(r1).M());
              }
            }
          }
        }
      }
      for (size_t p1 = 0; p1 < vecAntiProton.size(); p1++) {
        for (size_t p2 = p1 + 1; p2 < vecAntiProton.size(); p2++) {
          for (size_t r1 = 0; r1 < vecRho.size(); r1++) {
            if (idxAntiProton.at(p1) == idxRhoDaughNeg.at(r1) || idxAntiProton.at(p2) == idxRhoDaughNeg.at(r1)) {
              continue;
            }
            q3 = getQ3(vecAntiProton.at(p1), vecAntiProton.at(p2), vecRho.at(r1));
            registryTriggerQA.fill(HIST("PPRho/all/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PPRho/all/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PPRho/all/fSE_antiparticle"), q3);
            registryTriggerQA.fill(HIST("PPRho/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PPRho/all/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
            registryTriggerQA.fill(HIST("PPRho/all/fRhoQ3VsPt"), q3, vecRho.at(r1).Pt());
            registryTriggerQA.fill(HIST("PPRho/all/fRhoQ3VsInvMass"), q3, vecRho.at(r1).M());
            if (q3 < TriggerSelections.limits->get("Loose Limit", "PPRho")) {
              signalLooseLimit[cf_trigger::kPPRho] += 1;
              registryTriggerQA.fill(HIST("PPRho/loose/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PPRho/loose/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PPRho/loose/fSE_antiparticle"), q3);
              registryTriggerQA.fill(HIST("PPRho/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PPRho/loose/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
              registryTriggerQA.fill(HIST("PPRho/loose/fRhoQ3VsPt"), q3, vecRho.at(r1).Pt());
              registryTriggerQA.fill(HIST("PPRho/loose/fRhoQ3VsInvMass"), q3, vecRho.at(r1).M());
              if (q3 < TriggerSelections.limits->get("Tight Limit", "PPRho") &&
                  vecRho.at(r1).M() > RhoSelections.tightInvMassLow.value && vecRho.at(r1).M() < RhoSelections.tightInvMassUp.value) {
                signalTightLimit[cf_trigger::kPPRho] += 1;
                registryTriggerQA.fill(HIST("PPRho/tight/fMultiplicity"), col.multNTracksPV());
                registryTriggerQA.fill(HIST("PPRho/tight/fZvtx"), col.posZ());
                registryTriggerQA.fill(HIST("PPRho/tight/fSE_antiparticle"), q3);
                registryTriggerQA.fill(HIST("PPRho/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p1).Pt());
                registryTriggerQA.fill(HIST("PPRho/tight/fAntiProtonQ3VsPt"), q3, vecAntiProton.at(p2).Pt());
                registryTriggerQA.fill(HIST("PPRho/tight/fRhoQ3VsPt"), q3, vecRho.at(r1).Pt());
                registryTriggerQA.fill(HIST("PPRho/tight/fRhoQ3VsInvMass"), q3, vecRho.at(r1).M());
              }
            }
          }
        }
      }
    }
    // PD
    if (TriggerSelections.filterSwitches->get("Switch", "PD") > 0) {
      for (size_t p1 = 0; p1 < vecProton.size(); p1++) {
        for (size_t d1 = 0; d1 < vecDeuteron.size(); d1++) {
          if (idxProton.at(p1) == idxDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecProton.at(p1), vecDeuteron.at(d1));
          registryTriggerQA.fill(HIST("PD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("PD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("PD/all/fSE_particle"), kstar);
          registryTriggerQA.fill(HIST("PD/all/fProtonKstarVsPt"), kstar, vecProton.at(p1).Pt());
          registryTriggerQA.fill(HIST("PD/all/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "PD")) {
            signalLooseLimit[cf_trigger::kPD] += 1;
            registryTriggerQA.fill(HIST("PD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PD/loose/fSE_particle"), kstar);
            registryTriggerQA.fill(HIST("PD/loose/fProtonKstarVsPt"), kstar, vecProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PD/loose/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "PD")) {
              signalTightLimit[cf_trigger::kPD] += 1;
              registryTriggerQA.fill(HIST("PD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PD/tight/fSE_particle"), kstar);
              registryTriggerQA.fill(HIST("PD/tight/fProtonKstarVsPt"), kstar, vecProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PD/tight/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
            }
          }
        }
      }
      for (size_t p1 = 0; p1 < vecAntiProton.size(); p1++) {
        for (size_t d1 = 0; d1 < vecAntiDeuteron.size(); d1++) {
          if (idxAntiProton.at(p1) == idxAntiDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecAntiProton.at(p1), vecAntiDeuteron.at(d1));
          registryTriggerQA.fill(HIST("PD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("PD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("PD/all/fSE_antiparticle"), kstar);
          registryTriggerQA.fill(HIST("PD/all/fAntiProtonKstarVsPt"), kstar, vecAntiProton.at(p1).Pt());
          registryTriggerQA.fill(HIST("PD/all/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "PD")) {
            signalLooseLimit[cf_trigger::kPD] += 1;
            registryTriggerQA.fill(HIST("PD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PD/loose/fSE_antiparticle"), kstar);
            registryTriggerQA.fill(HIST("PD/loose/fAntiProtonKstarVsPt"), kstar, vecAntiProton.at(p1).Pt());
            registryTriggerQA.fill(HIST("PD/loose/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "PD")) {
              signalTightLimit[cf_trigger::kPD] += 1;
              registryTriggerQA.fill(HIST("PD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PD/tight/fSE_antiparticle"), kstar);
              registryTriggerQA.fill(HIST("PD/tight/fAntiProtonKstarVsPt"), kstar, vecAntiProton.at(p1).Pt());
              registryTriggerQA.fill(HIST("PD/tight/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
            }
          }
        }
      }
    }
    // LD
    if (TriggerSelections.filterSwitches->get("Switch", "LD") > 0) {
      for (size_t l1 = 0; l1 < vecLambda.size(); l1++) {
        for (size_t d1 = 0; d1 < vecDeuteron.size(); d1++) {
          if (idxLambdaDaughProton.at(l1) == idxDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecLambda.at(l1), vecDeuteron.at(d1));
          registryTriggerQA.fill(HIST("LD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("LD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("LD/all/fSE_particle"), kstar);
          registryTriggerQA.fill(HIST("LD/all/fLambdaKstarVsPt"), kstar, vecLambda.at(l1).Pt());
          registryTriggerQA.fill(HIST("LD/all/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "LD")) {
            signalLooseLimit[cf_trigger::kLD] += 1;
            registryTriggerQA.fill(HIST("LD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("LD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("LD/loose/fSE_particle"), kstar);
            registryTriggerQA.fill(HIST("LD/loose/fLambdaKstarVsPt"), kstar, vecLambda.at(l1).Pt());
            registryTriggerQA.fill(HIST("LD/loose/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "LD")) {
              signalTightLimit[cf_trigger::kLD] += 1;
              registryTriggerQA.fill(HIST("LD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("LD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("LD/tight/fSE_particle"), kstar);
              registryTriggerQA.fill(HIST("LD/tight/fLambdaKstarVsPt"), kstar, vecLambda.at(l1).Pt());
              registryTriggerQA.fill(HIST("LD/tight/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
            }
          }
        }
      }
      for (size_t l1 = 0; l1 < vecAntiLambda.size(); l1++) {
        for (size_t d1 = 0; d1 < vecAntiDeuteron.size(); d1++) {
          if (idxAntiLambdaDaughProton.at(l1) == idxAntiDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecAntiLambda.at(l1), vecAntiDeuteron.at(d1));
          registryTriggerQA.fill(HIST("LD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("LD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("LD/all/fSE_antiparticle"), kstar);
          registryTriggerQA.fill(HIST("LD/all/fAntiLambdaKstarVsPt"), kstar, vecAntiLambda.at(l1).Pt());
          registryTriggerQA.fill(HIST("LD/all/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "LD")) {
            signalLooseLimit[cf_trigger::kLD] += 1;
            registryTriggerQA.fill(HIST("LD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("LD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("LD/loose/fSE_antiparticle"), kstar);
            registryTriggerQA.fill(HIST("LD/loose/fAntiLambdaKstarVsPt"), kstar, vecAntiLambda.at(l1).Pt());
            registryTriggerQA.fill(HIST("LD/loose/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "LD")) {
              signalTightLimit[cf_trigger::kLD] += 1;
              registryTriggerQA.fill(HIST("LD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("LD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("LD/tight/fSE_antiparticle"), kstar);
              registryTriggerQA.fill(HIST("LD/tight/fAntiLambdaKstarVsPt"), kstar, vecAntiLambda.at(l1).Pt());
              registryTriggerQA.fill(HIST("LD/tight/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
            }
          }
        }
      }
    }
    // PhiD
    if (TriggerSelections.filterSwitches->get("Switch", "PhiD") > 0) {
      for (size_t phi1 = 0; phi1 < vecPhi.size(); phi1++) {
        for (size_t d1 = 0; d1 < vecDeuteron.size(); d1++) {
          if (idxPhiDaughPos.at(phi1) == idxDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecPhi.at(phi1), vecDeuteron.at(d1));
          registryTriggerQA.fill(HIST("PhiD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("PhiD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("PhiD/all/fSE_particle"), kstar);
          registryTriggerQA.fill(HIST("PhiD/all/fPhiKstarVsPt"), kstar, vecPhi.at(phi1).Pt());
          registryTriggerQA.fill(HIST("PhiD/all/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
          registryTriggerQA.fill(HIST("PhiD/all/fPhiKstarVsInvMass"), kstar, vecPhi.at(phi1).M());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "PhiD")) {
            signalLooseLimit[cf_trigger::kPhiD] += 1;
            registryTriggerQA.fill(HIST("PhiD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PhiD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PhiD/loose/fSE_particle"), kstar);
            registryTriggerQA.fill(HIST("PhiD/loose/fPhiKstarVsPt"), kstar, vecPhi.at(phi1).Pt());
            registryTriggerQA.fill(HIST("PhiD/loose/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
            registryTriggerQA.fill(HIST("PhiD/loose/fPhiKstarVsInvMass"), kstar, vecPhi.at(phi1).M());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "PhiD") &&
                vecPhi.at(phi1).M() > PhiSelections.tightInvMassLow.value && vecPhi.at(phi1).M() < PhiSelections.tightInvMassUp.value) {
              signalTightLimit[cf_trigger::kPhiD] += 1;
              registryTriggerQA.fill(HIST("PhiD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PhiD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PhiD/tight/fSE_particle"), kstar);
              registryTriggerQA.fill(HIST("PhiD/tight/fPhiKstarVsPt"), kstar, vecPhi.at(phi1).Pt());
              registryTriggerQA.fill(HIST("PhiD/tight/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
              registryTriggerQA.fill(HIST("PhiD/tight/fPhiKstarVsInvMass"), kstar, vecPhi.at(phi1).M());
            }
          }
        }
      }
      for (size_t phi1 = 0; phi1 < vecPhi.size(); phi1++) {
        for (size_t d1 = 0; d1 < vecAntiDeuteron.size(); d1++) {
          if (idxPhiDaughNeg.at(phi1) == idxAntiDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecPhi.at(phi1), vecAntiDeuteron.at(d1));
          registryTriggerQA.fill(HIST("PhiD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("PhiD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("PhiD/all/fSE_antiparticle"), kstar);
          registryTriggerQA.fill(HIST("PhiD/all/fPhiKstarVsPt"), kstar, vecPhi.at(phi1).Pt());
          registryTriggerQA.fill(HIST("PhiD/all/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
          registryTriggerQA.fill(HIST("PhiD/all/fPhiKstarVsInvMass"), kstar, vecPhi.at(phi1).M());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "PhiD")) {
            signalLooseLimit[cf_trigger::kPhiD] += 1;
            registryTriggerQA.fill(HIST("PhiD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("PhiD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("PhiD/loose/fSE_antiparticle"), kstar);
            registryTriggerQA.fill(HIST("PhiD/loose/fPhiKstarVsPt"), kstar, vecPhi.at(phi1).Pt());
            registryTriggerQA.fill(HIST("PhiD/loose/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
            registryTriggerQA.fill(HIST("PhiD/loose/fPhiKstarVsInvMass"), kstar, vecPhi.at(phi1).M());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "PhiD") &&
                vecPhi.at(phi1).M() > PhiSelections.tightInvMassLow.value && vecPhi.at(phi1).M() < PhiSelections.tightInvMassUp.value) {
              signalTightLimit[cf_trigger::kPhiD] += 1;
              registryTriggerQA.fill(HIST("PhiD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("PhiD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("PhiD/tight/fSE_antiparticle"), kstar);
              registryTriggerQA.fill(HIST("PhiD/tight/fPhiKstarVsPt"), kstar, vecPhi.at(phi1).Pt());
              registryTriggerQA.fill(HIST("PhiD/tight/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
              registryTriggerQA.fill(HIST("PhiD/tight/fPhiKstarVsInvMass"), kstar, vecPhi.at(phi1).M());
            }
          }
        }
      }
    }
    // RhoD
    if (TriggerSelections.filterSwitches->get("Switch", "RhoD") > 0) {
      for (size_t r1 = 0; r1 < vecRho.size(); r1++) {
        for (size_t d1 = 0; d1 < vecDeuteron.size(); d1++) {
          if (idxRhoDaughPos.at(r1) == idxDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecRho.at(r1), vecDeuteron.at(d1));
          registryTriggerQA.fill(HIST("RhoD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("RhoD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("RhoD/all/fSE_particle"), kstar);
          registryTriggerQA.fill(HIST("RhoD/all/fRhoKstarVsPt"), kstar, vecRho.at(r1).Pt());
          registryTriggerQA.fill(HIST("RhoD/all/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
          registryTriggerQA.fill(HIST("RhoD/all/fRhoKstarVsInvMass"), kstar, vecRho.at(r1).M());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "RhoD")) {
            signalLooseLimit[cf_trigger::kRhoD] += 1;
            registryTriggerQA.fill(HIST("RhoD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("RhoD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("RhoD/loose/fSE_particle"), kstar);
            registryTriggerQA.fill(HIST("RhoD/loose/fRhoKstarVsPt"), kstar, vecRho.at(r1).Pt());
            registryTriggerQA.fill(HIST("RhoD/loose/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
            registryTriggerQA.fill(HIST("RhoD/loose/fRhoKstarVsInvMass"), kstar, vecRho.at(r1).M());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "RhoD") &&
                vecRho.at(r1).M() > RhoSelections.tightInvMassLow.value && vecRho.at(r1).M() < RhoSelections.tightInvMassUp.value) {
              signalTightLimit[cf_trigger::kRhoD] += 1;
              registryTriggerQA.fill(HIST("RhoD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("RhoD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("RhoD/tight/fSE_particle"), kstar);
              registryTriggerQA.fill(HIST("RhoD/tight/fRhoKstarVsPt"), kstar, vecRho.at(r1).Pt());
              registryTriggerQA.fill(HIST("RhoD/tight/fDeuteronKstarVsPt"), kstar, vecDeuteron.at(d1).Pt());
              registryTriggerQA.fill(HIST("RhoD/tight/fRhoKstarVsInvMass"), kstar, vecRho.at(r1).M());
            }
          }
        }
      }
      for (size_t r1 = 0; r1 < vecRho.size(); r1++) {
        for (size_t d1 = 0; d1 < vecAntiDeuteron.size(); d1++) {
          if (idxRhoDaughNeg.at(r1) == idxAntiDeuteron.at(d1)) {
            continue;
          }
          kstar = getkstar(vecRho.at(r1), vecAntiDeuteron.at(d1));
          registryTriggerQA.fill(HIST("RhoD/all/fMultiplicity"), col.multNTracksPV());
          registryTriggerQA.fill(HIST("RhoD/all/fZvtx"), col.posZ());
          registryTriggerQA.fill(HIST("RhoD/all/fSE_antiparticle"), kstar);
          registryTriggerQA.fill(HIST("RhoD/all/fRhoKstarVsPt"), kstar, vecRho.at(r1).Pt());
          registryTriggerQA.fill(HIST("RhoD/all/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
          registryTriggerQA.fill(HIST("RhoD/all/fRhoKstarVsInvMass"), kstar, vecRho.at(r1).M());
          if (kstar < TriggerSelections.limits->get("Loose Limit", "RhoD")) {
            signalLooseLimit[cf_trigger::kRhoD] += 1;
            registryTriggerQA.fill(HIST("RhoD/loose/fMultiplicity"), col.multNTracksPV());
            registryTriggerQA.fill(HIST("RhoD/loose/fZvtx"), col.posZ());
            registryTriggerQA.fill(HIST("RhoD/loose/fSE_antiparticle"), kstar);
            registryTriggerQA.fill(HIST("RhoD/loose/fRhoKstarVsPt"), kstar, vecRho.at(r1).Pt());
            registryTriggerQA.fill(HIST("RhoD/loose/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
            registryTriggerQA.fill(HIST("RhoD/loose/fRhoKstarVsInvMass"), kstar, vecRho.at(r1).M());
            if (kstar < TriggerSelections.limits->get("Tight Limit", "RhoD") &&
                vecRho.at(r1).M() > RhoSelections.tightInvMassLow.value && vecRho.at(r1).M() < RhoSelections.tightInvMassUp.value) {
              signalTightLimit[cf_trigger::kRhoD] += 1;
              registryTriggerQA.fill(HIST("RhoD/tight/fMultiplicity"), col.multNTracksPV());
              registryTriggerQA.fill(HIST("RhoD/tight/fZvtx"), col.posZ());
              registryTriggerQA.fill(HIST("RhoD/tight/fSE_antiparticle"), kstar);
              registryTriggerQA.fill(HIST("RhoD/tight/fRhoKstarVsPt"), kstar, vecRho.at(r1).Pt());
              registryTriggerQA.fill(HIST("RhoD/tight/fAntiDeuteronKstarVsPt"), kstar, vecAntiDeuteron.at(d1).Pt());
              registryTriggerQA.fill(HIST("RhoD/tight/fRhoKstarVsInvMass"), kstar, vecRho.at(r1).M());
            }
          }
        }
      }
    }

    for (int i = 0; i < cf_trigger::kNTriggers; i++) {
      if (signalLooseLimit[i] > 0) {
        registryTriggerQA.fill(HIST("fProcessedEvents"), 3 + 2 * i); // need offset for filling
        keepEventLooseLimit[i] = true;
      }
      if (signalTightLimit[i] > 0) {
        registryTriggerQA.fill(HIST("fProcessedEvents"), 3 + 2 * i + 1); // need offset for filling
        keepEventTightLimit[i] = true;
      }
      for (int j = i; j < cf_trigger::kNTriggers; j++) {
        if (signalLooseLimit[i] > 0 && signalLooseLimit[j]) {
          registryTriggerQA.fill(HIST("fTriggerCorrelations"), 2 * i, 2 * j);
        }
        if (signalLooseLimit[i] > 0 && signalTightLimit[j]) { // only one combination needed, fill only entries above diagonal
          registryTriggerQA.fill(HIST("fTriggerCorrelations"), 2 * i, 2 * j + 1);
        }
        if (signalTightLimit[i] > 0 && signalTightLimit[j]) {
          registryTriggerQA.fill(HIST("fTriggerCorrelations"), 2 * i + 1, 2 * j + 1);
        }
      }
    }

    if (keepEventLooseLimit[cf_trigger::kPPP] ||
        keepEventLooseLimit[cf_trigger::kPPL] ||
        keepEventLooseLimit[cf_trigger::kPLL] ||
        keepEventLooseLimit[cf_trigger::kLLL] ||
        keepEventLooseLimit[cf_trigger::kPPPhi] ||
        keepEventLooseLimit[cf_trigger::kPPRho] ||
        keepEventLooseLimit[cf_trigger::kPD] ||
        keepEventLooseLimit[cf_trigger::kLD] ||
        keepEventLooseLimit[cf_trigger::kPhiD] ||
        keepEventLooseLimit[cf_trigger::kRhoD]) {
      registryTriggerQA.fill(HIST("fProcessedEvents"), 1);
    }

    if (keepEventTightLimit[cf_trigger::kPPP] ||
        keepEventTightLimit[cf_trigger::kPPL] ||
        keepEventTightLimit[cf_trigger::kPLL] ||
        keepEventTightLimit[cf_trigger::kLLL] ||
        keepEventTightLimit[cf_trigger::kPPPhi] ||
        keepEventTightLimit[cf_trigger::kPPRho] ||
        keepEventTightLimit[cf_trigger::kPD] ||
        keepEventTightLimit[cf_trigger::kLD] ||
        keepEventTightLimit[cf_trigger::kPhiD] ||
        keepEventTightLimit[cf_trigger::kRhoD]) {
      registryTriggerQA.fill(HIST("fProcessedEvents"), 2);
    }

    tags(keepEventTightLimit[cf_trigger::kPPP], keepEventLooseLimit[cf_trigger::kPPP],
         keepEventTightLimit[cf_trigger::kPPL], keepEventLooseLimit[cf_trigger::kPPL],
         keepEventTightLimit[cf_trigger::kPLL], keepEventLooseLimit[cf_trigger::kPLL],
         keepEventTightLimit[cf_trigger::kLLL], keepEventLooseLimit[cf_trigger::kLLL],
         keepEventTightLimit[cf_trigger::kPPPhi], keepEventLooseLimit[cf_trigger::kPPPhi],
         keepEventTightLimit[cf_trigger::kPPRho], keepEventLooseLimit[cf_trigger::kPPRho],
         keepEventTightLimit[cf_trigger::kPD], keepEventLooseLimit[cf_trigger::kPD],
         keepEventTightLimit[cf_trigger::kLD], keepEventLooseLimit[cf_trigger::kLD],
         keepEventTightLimit[cf_trigger::kPhiD], keepEventLooseLimit[cf_trigger::kPhiD],
         keepEventTightLimit[cf_trigger::kRhoD], keepEventLooseLimit[cf_trigger::kRhoD]);
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilterAll>(cfg)};
}
