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
///
/// Making modifications to the Strangeness Tutorial code
/// The code is still in development mode
/// Flattenicity part of the code is adopted from
/// https://github.com/AliceO2Group/O2Physics/blob/master/PWGMM/Mult/Tasks/flatenicityFV0.cxx
/// \file lambdak0sflattenicity.cxx
/// \brief V0 task for production of strange hadrons as a function of flattenicity
/// \author Suraj Prasad (suraj.prasad@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>
#include <TProfile2D.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Lambdak0sflattenicity {
  // Histograms are defined with HistogramRegistry
  Service<o2::framework::O2DatabasePDG> pdg{};
  HistogramRegistry rEventSelection{"eventSelection",
                                    {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true,
                                    true};
  HistogramRegistry rKzeroShort{
    "kzeroShort",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry rLambda{
    "lambda",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry rAntiLambda{
    "antilambda",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry rXi{
    "Xi",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry rCommonHist{
    "commonhists",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry rFlattenicity{
    "flattenicity",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry rCharged{
    "charged",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  static constexpr int kNEstimators = 8;
  // forward detector segmentation
  static constexpr int kNChannelsPerT0Sector = 4;
  static constexpr int kNSectorsT0A = 24;
  static constexpr int kNSectorsT0C = 28;
  static constexpr int kNChannelsPerFV0Ring = 8;
  static constexpr int kNFV0EtaRings = 5;
  static constexpr int kOuterFV0RingIndex = kNFV0EtaRings - 1;
  static constexpr int kNSectorsFV0OuterRing = 16;
  // TParticlePDG::Charge() is in units of e/3
  static constexpr float kMinCharge = 0.01f;

  static constexpr std::array<std::string_view, kNEstimators> kHEst = {
    "eGlobaltrack", "eFV0", "e1flatencityFV0", "eFT0",
    "e1flatencityFT0", "eFV0FT0C", "e1flatencityFV0FT0C", "ePtTrig"};
  static constexpr std::array<std::string_view, kNEstimators> kTEst = {
    "GlobalTrk", "FV0", "1-flatencity_FV0", "FT0",
    "1-flatencityFT0", "FV0_FT0C", "1-flatencity_FV0_FT0C", "PtTrig"};
  static constexpr std::array<std::string_view, kNEstimators> kHPtEst = {
    "ptVsGlobaltrack", "ptVsFV0",
    "ptVs1flatencityFV0", "ptVsFT0",
    "ptVs1flatencityFT0", "ptVsFV0FT0C",
    "ptVs1flatencityFV0FT0C", "pTVsPtTrig"};

  // Histogram binning. The pT and 1-rho axes are ConfigurableAxis: pass
  // {nbins, lo, hi} for uniform bins or {VARIABLE_WIDTH, e0, e1, ...} for the
  // edges of the published binning. The 1-rho edges are meant to be the
  // percentile class boundaries measured by the processFlatDist* pass, so that
  // one class is exactly one bin and no rebinning happens downstream.
  struct : ConfigurableGroup {
    std::string prefix = "binning";
    Configurable<int> nBinsVz{"nBinsVz", 100, "N bins in Vz"};
    Configurable<int> nBinsK0sMass{"nBinsK0sMass", 400, "N bins in K0sMass"};
    Configurable<int> nBinsLambdaMass{"nBinsLambdaMass", 400,
                                      "N bins in LambdaMass"};
    Configurable<int> nBinsXiMass{"nBinsXiMass", 400, "N bins in XiMass"};
    Configurable<float> kK0sEPshiftfromMass{"kK0sEPshiftfromMass", 0.1, "distance of K0s Inv mass histogram start and end points from PDG mass"};
    Configurable<float> kLambdaEPshiftfromMass{"kLambdaEPshiftfromMass", 0.05, "distance of Lambda Inv mass histogram start and end points from PDG mass"};
    Configurable<float> kXiEPshiftfromMass{"kXiEPshiftfromMass", 0.05, "distance of Xi Inv mass histogram start and end points from PDG mass"};
    ConfigurableAxis axisPtK0s{"axisPtK0s", {250, 0.0f, 25.0f}, "#it{p}_{T} of K0s"};
    ConfigurableAxis axisPtLambda{"axisPtLambda", {250, 0.0f, 25.0f}, "#it{p}_{T} of Lambda and AntiLambda"};
    ConfigurableAxis axisPtXi{"axisPtXi", {250, 0.0f, 25.0f}, "#it{p}_{T} of Xi"};
    ConfigurableAxis axisPtCharged{"axisPtCharged", {250, 0.0f, 25.0f}, "#it{p}_{T} of charged particles"};
    ConfigurableAxis axisPtPid{"axisPtPid", {250, 0.0f, 25.0f}, "#it{p}_{TPC} of the daughter tracks"};
    ConfigurableAxis axisFlat{"axisFlat", {100, 0.0f, 1.0f}, "1-#rho_{ch} of the spectra"};
    // detector effects move the percentiles, so the closure denominator needs
    // its own class edges, taken from hFlatDistGen
    ConfigurableAxis axisFlatTrue{"axisFlatTrue", {100, 0.0f, 1.0f}, "true 1-#rho_{ch}"};
    ConfigurableAxis axisFlatFine{"axisFlatFine", {2000, 0.0f, 1.0f}, "1-#rho_{ch} used to locate the percentiles"};
    ConfigurableAxis axisCent{"axisCent", {100, 0.0f, 100.0f}, "FT0M percentile"};
    ConfigurableAxis axisNch{"axisNch", {150, -0.5f, 149.5f}, "N_{ch} in the tracking acceptance"};
    ConfigurableAxis axisDcaXy{"axisDcaXy", {200, -1.0f, 1.0f}, "DCA_{xy} (cm)"};
    ConfigurableAxis axisDcaV0ToPv{"axisDcaV0ToPv", {200, 0.0f, 2.0f}, "DCA of the V0 to the PV (cm)"};
    ConfigurableAxis axisPtRes{"axisPtRes", {200, -0.5f, 0.5f}, "(#it{p}_{T}^{rec} - #it{p}_{T}^{gen})/#it{p}_{T}^{gen}"};
  } binning;

  // Event selection
  struct : ConfigurableGroup {
    std::string prefix = "evSel";
    Configurable<bool> applyEvSel{"applyEvSel", true,
                                  "Apply event selection to Data and MCRec"};
    Configurable<bool> issel8{"issel8", true,
                              "Accept events that pass sel8 selection"};
    Configurable<float> cutzvertex{"cutzvertex", 10.0f,
                                   "Accepted z-vertex range (cm)"};
    Configurable<bool> isINELgt0{"isINELgt0", true, "is INEL gt 0"};
    Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", true,
                                           "cut branch crossing at the beginning/end of TF"};
    Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", true,
                                            "cut branch crossing at the beginning/end of ITS ROF"};
    Configurable<bool> isVertexITSTPC{"isVertexITSTPC", false,
                                      "Is Vertex ITSTPC"};
    Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", false,
                                           "Is No Same Bunch Pileup"};
    Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", false,
                                         "Is Good Zvtx FT0 vs PV"};
    Configurable<bool> isTriggerTVX{"isTriggerTVX", true,
                                    "coincidence of a signal in FT0A and FT0C"};
  } evSel;

  // Flattenicity estimation
  struct : ConfigurableGroup {
    std::string prefix = "flatSel";
    Configurable<bool> flattenicityQA{"flattenicityQA", true, "Store Flattenicity QA plots"};
    // Both corrections reshape the lattice and therefore move the flattenicity.
    // With both off nothing is read from CCDB at all.
    Configurable<bool> applyCalibCh{"applyCalibCh", false,
                                    "equalize the per-channel gain, from CCDB"};
    Configurable<bool> applyCalibVtx{"applyCalibVtx", false,
                                     "equalize the per-cell z-vertex dependence, from CCDB"};
    Configurable<bool> applyNorm{"applyNorm", false, "normalization to eta"};
    Configurable<bool> isflattenicitywithFV0{"isflattenicitywithFV0", true,
                                             "Calculate Flattenicity with FV0"};
    Configurable<bool> isflattenicitywithFT0{"isflattenicitywithFT0", true,
                                             "Calculate Flattenicity with FT0"};
    Configurable<bool> isflattenicitywithFV0FT0C{"isflattenicitywithFV0FT0C", true,
                                                 "Calculate Flattenicity with FV0+FT0C"};
    Configurable<int> flattenicityforanalysis{"flattenicityforanalysis", 0,
                                              "Which Flattenicity to be used for analysis, 0 for FV0, 1 for FT0, 2 for FV0+FT0C"};
    Configurable<bool> flattenicityforLossCorrRec{"flattenicityforLossCorrRec", true,
                                                  "Flattenicity from Rec Tracks are used for Signal and Event loss calculations"};
    // same cell convention as the detector lattice, so that gen vs rec is a
    // detector effect and not a difference of definitions
    Configurable<bool> genFlatDetectorLikeNorm{"genFlatDetectorLikeNorm", true,
                                               "Weight the generator-level FV0 cells the way the detector lattice does"};
  } flatSel;

  // Calibration objects, read per run and only when the corresponding switch is
  // on. Same layout as PWGLF/Tasks/GlobalEventProperties/flattenictyPikp.cxx, so
  // the objects are interchangeable between the two tasks: the gain is a
  // std::vector<float> indexed by the raw detector channel, the z-vertex
  // equalization a TProfile2D of (lattice cell, z_vtx).
  struct : ConfigurableGroup {
    std::string prefix = "ccdbConf";
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch",
                                      "url of the ccdb repository"};
    Configurable<std::string> gainEqPath{"gainEqPath", "Users/s/sprasad/flattenicity/GainEq",
                                         "CCDB directory holding FV0, FT0A and FT0C gain vectors"};
    Configurable<std::string> vtxEqPath{"vtxEqPath", "Users/s/sprasad/flattenicity/ZvtxEq",
                                        "CCDB directory holding the FV0, FT0A and FT0C z-vertex maps"};
  } ccdbConf;

  // V0 selection
  struct : ConfigurableGroup {
    std::string prefix = "v0Sel";
    Configurable<float> v0settingDCAv0dau{"v0settingDCAv0dau", 1,
                                          "DCA V0 Daughters"};
    Configurable<float> v0settingDCApostopv{"v0settingDCApostopv", 0.06,
                                            "DCA Pos To PV"};
    Configurable<float> v0settingDCAnegtopv{"v0settingDCAnegtopv", 0.06,
                                            "DCA Neg To PV"};
    Configurable<float> v0settingDCAbactopv{"v0settingDCAbactopv", 0.06,
                                            "DCA Bchelor To PV"};
    Configurable<float> v0settingRapidity{"v0settingRapidity", 0.5,
                                          "V0 rapidity cut"};
    Configurable<float> v0settingCosPAK0s{"v0settingCosPAK0s", 0.97,
                                          "V0 CosPA for K0s"};
    Configurable<float> v0settingRadiusK0s{"v0settingRadiusK0s", 0.5,
                                           "v0radius for K0s"};
    Configurable<float> v0settingcTauK0s{"v0settingcTauK0s", 20,
                                         "v0ctau for K0s"};
    Configurable<float> v0settingMassRejectionK0s{"v0settingMassRejectionK0s", 0.005,
                                                  "Competing Mass Rejection cut for K0s"};
    Configurable<float> v0settingArmePodoK0s{"v0settingArmePodoK0s", 0.2,
                                             "Armenteros-Podolanski cut for K0s"};
    Configurable<float> v0settingCosPALambda{"v0settingCosPALambda", 0.995,
                                             "V0 CosPA for Lambda"};
    Configurable<float> v0settingRadiusLambda{"v0settingRadiusLambda", 0.5,
                                              "v0radius for Lambda"};
    Configurable<float> v0settingcTauLambda{"v0settingcTauLambda", 30,
                                            "v0ctau for Lambda"};
    Configurable<float> v0settingMassRejectionLambda{"v0settingMassRejectionLambda", 0.01,
                                                     "Competing Mass Rejection cut for Lambda"};
    Configurable<float> v0settingNTPCcrossedRows{"v0settingNTPCcrossedRows", 70,
                                                 "Minimum number of TPC crossed pad rows of the V0 daughters, negative: no cut"};
  } v0Sel;

  // Daughter track acceptance and PID
  struct : ConfigurableGroup {
    std::string prefix = "trkPid";
    Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 0.8f,
                                     "Eta range for tracks"};
    Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.0f, "Minimum  pT"};
    Configurable<float> nSigmaTPCPion{"nSigmaTPCPion", 5, "nSigmaTPCPion"};
    Configurable<float> nSigmaTPCProton{"nSigmaTPCProton", 5, "nSigmaTPCProton"};
    Configurable<float> pidQAWindowK0s{"pidQAWindowK0s", 0.05,
                                       "Half width of the K0s mass window used for the PID QA plots"};
    Configurable<float> pidQAWindowLambda{"pidQAWindowLambda", 0.1,
                                          "Half width of the Lambda mass window used for the PID QA plots"};
    // negative values disable a cut, so the defaults reproduce the previous selection
    Configurable<float> minCrossedRowsOverFindable{"minCrossedRowsOverFindable", -1.f,
                                                   "Minimum TPC crossed rows over findable clusters"};
    Configurable<float> maxTpcChi2NCl{"maxTpcChi2NCl", -1.f, "Maximum TPC chi2 per cluster"};
    Configurable<float> maxItsChi2NCl{"maxItsChi2NCl", -1.f, "Maximum ITS chi2 per cluster"};
    Configurable<int> minItsNCls{"minItsNCls", -1, "Minimum number of ITS clusters"};
  } trkPid;

  // Cascade selection
  struct : ConfigurableGroup {
    std::string prefix = "cascSel";
    Configurable<float> nTPCcrossedRows{"nTPCcrossedRows", 52, "Number of TPC crossed pad raws"};
    Configurable<float> cascsettingDCAv0toPV{"cascsettingDCAv0toPV", 0.03, "DCA V0 To PV"};
    Configurable<float> cascsettingDCAv0bach{"cascsettingDCAv0bach", 0.25, "DCA V0 To bachelor"};
    Configurable<float> cascsettingDCAbaryontopv{"cascsettingDCAbaryontopv", 0.05, "DCA of the baryon daughter To PV"};
    Configurable<float> cascsettingDCAmesontopv{"cascsettingDCAmesontopv", 0.1, "DCA of the meson daughter To PV"};
    Configurable<float> cascsettingDCAxybaryonbach{"cascsettingDCAxybaryonbach", 0.02, "DCAxy Bachelor-Baryon To PV, below which the candidate is vetoed"};
    Configurable<float> cascsettingCosPAbaryonbach{"cascsettingCosPAbaryonbach", 0.9999, "CosThetap Bachelor-Baryon, above which the candidate is vetoed"};
    Configurable<float> cascsettingCosPAcascPV{"cascsettingCosPAcascPV", 0.9947, "CosThetap for Cascade to PV"};
    Configurable<float> cascsettingCosPAv0PV{"cascsettingCosPAv0PV", 0.9876, "CosThetap for V0 to PV"};
    Configurable<float> cascsettingv0radius{"cascsettingv0radius", 0.55, "V0 decay radius for cadcades in cm"};
    Configurable<float> cascsettingcascradius{"cascsettingcascradius", 1.01, "Cascade decay radius for cadcades in cm"};
    Configurable<float> cascsettingRapidity{"cascsettingRapidity", 0.5, "Cascade rapidity cut"};
    Configurable<float> cascsettingMassRejectionLambdaXi{"cascsettingMassRejectionLambdaXi", 0.0116, "Casc Mass Rejection cut of Lambda for Xi"};
    Configurable<float> cascsettingMassRejectioOmegaXi{"cascsettingMassRejectioOmegaXi", -1, "Casc Mass Rejection cut of Omega for Xi"};
    Configurable<float> cascsettingproplifetime{"cascsettingproplifetime", 4.6, "Scale for lifetime cut on ctau Xi"};
  } cascSel;

  // values for flattenicityforanalysis
  static constexpr int kFlatFromFV0 = 0;
  static constexpr int kFlatFromFT0 = 1;
  static constexpr int kFlatFromFV0FT0C = 2;

  // FT0M percentile class: run once over the full range for MB and once with a
  // narrow window for HM. Grouped because the framework only decomposes 100
  // struct members (Framework/StructToTuple.h)
  struct : ConfigurableGroup {
    std::string prefix = "eventClass";
    Configurable<bool> applyCentSel{"applyCentSel", false,
                                    "Select events in a FT0M percentile window"};
    Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Minimum FT0M percentile"};
    Configurable<float> cfgCentMax{"cfgCentMax", 100.0f, "Maximum FT0M percentile"};
    // the processFlatDist* pass fills the MB and the HM 1-rho distributions
    // together, so this window is applied to a histogram and not to the event
    Configurable<float> cfgCentMaxHM{"cfgCentMaxHM", 1.0f,
                                     "Upper FT0M percentile of the high-multiplicity class"};

    // keep only primaries in the MC-matched spectra, the rest go to the
    // feed-down histograms
    Configurable<bool> requirePrimaryMC{"requirePrimaryMC", true,
                                        "Require isPhysicalPrimary() on the MC-matched candidate"};
    Configurable<bool> genFlatPrimariesOnly{"genFlatPrimariesOnly", true,
                                            "Use only primaries in the generator-level flattenicity"};
    Configurable<bool> fillChargedQA{"fillChargedQA", true,
                                     "Fill the charged-particle histograms used for <dNch/deta> and Qpp"};
    Configurable<float> cfgEtaChargedCut{"cfgEtaChargedCut", 0.8f,
                                         "Eta window of the charged-particle measurement"};
    // mothers outside the Lambda window still feed it, so the matrix is
    // normalised in a wider one
    Configurable<float> cfgFeedDownMotherRapidity{"cfgFeedDownMotherRapidity", 1.5f,
                                                  "Rapidity window of the generated feed-down mothers"};
  } eventClass;

  // Configurable<float> v0daughter_etacut{"V0DaughterEtaCut", 0.8,
  // "V0DaughterEtaCut"};
  // Configurable<float> v0etacut{"v0etacut", 0.8, "v0etacut"};

  int nbin = 1;
  // hEventsSelected bin for the flattenicity requirement, -1 if not applied
  int nbinFlattenicity = -1;
  // do not fill the detector QA twice when processGenMC runs with a rec-level process
  bool fillFlattenicityQAInGenMC = true;

  // the per-event histograms are shared by every rec-level process function,
  // so only the one designated in init() fills them
  static constexpr int kOwnerNone = -1;
  static constexpr int kOwnerFlatDistData = 0;
  static constexpr int kOwnerFlatDistMC = 1;
  static constexpr int kOwnerRecMCV0 = 2;
  static constexpr int kOwnerDataV0 = 3;
  static constexpr int kOwnerRecMCCasc = 4;
  static constexpr int kOwnerDataCasc = 5;
  int eventHistOwner = kOwnerNone;

  // calibration objects, refreshed when the run changes
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int mRunNumber = -1;
  std::vector<float> gainFV0;
  std::vector<float> gainFT0A;
  std::vector<float> gainFT0C;
  TProfile2D* vtxEqFV0 = nullptr;
  TProfile2D* vtxEqFT0A = nullptr;
  TProfile2D* vtxEqFT0C = nullptr;

  void init(InitContext const&)
  {
    // Axes
    AxisSpec k0sMassAxis = {binning.nBinsK0sMass, o2::constants::physics::MassK0Short - binning.kK0sEPshiftfromMass,
                            o2::constants::physics::MassK0Short + binning.kK0sEPshiftfromMass,
                            "#it{M}_{#pi^{+}#pi^{-}} [GeV/#it{c}^{2}]"};
    AxisSpec lambdaMassAxis = {binning.nBinsLambdaMass, o2::constants::physics::MassLambda0 - binning.kLambdaEPshiftfromMass,
                               o2::constants::physics::MassLambda0 + binning.kLambdaEPshiftfromMass,
                               "#it{M}_{p#pi^{-}} [GeV/#it{c}^{2}]"};
    AxisSpec antilambdaMassAxis = {binning.nBinsLambdaMass, o2::constants::physics::MassLambda0 - binning.kLambdaEPshiftfromMass,
                                   o2::constants::physics::MassLambda0 + binning.kLambdaEPshiftfromMass,
                                   "#it{M}_{#pi^{+}#bar{p}} [GeV/#it{c}^{2}]"};
    AxisSpec xiMassAxis = {binning.nBinsXiMass, o2::constants::physics::MassXiMinus - binning.kXiEPshiftfromMass,
                           o2::constants::physics::MassXiMinus + binning.kXiEPshiftfromMass,
                           "#it{M}_{#Lambda#pi} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {binning.nBinsVz, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptK0sAxis = {binning.axisPtK0s, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptLambdaAxis = {binning.axisPtLambda, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptXiAxis = {binning.axisPtXi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptChargedAxis = {binning.axisPtCharged, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPCAxis = {binning.axisPtPid, "#it{p}_{TPC} (GeV/#it{c})"};
    AxisSpec decayRadiusAxis = {100, 0.0f, 100.0f, "Decay Radius (cm)"};
    AxisSpec flatAxis = {binning.axisFlat, "1-#rho_{ch}"};
    AxisSpec flatTrueAxis = {binning.axisFlatTrue, "true 1-#rho_{ch}"};
    AxisSpec flatFineAxis = {binning.axisFlatFine, "1-#rho_{ch}"};
    AxisSpec centAxis = {binning.axisCent, "FT0M percentile"};
    AxisSpec nchAxis = {binning.axisNch, "#it{N}_{ch}"};
    AxisSpec dcaXyAxis = {binning.axisDcaXy, "DCA_{xy} (cm)"};
    AxisSpec dcaV0ToPvAxis = {binning.axisDcaV0ToPv, "DCA_{V0-PV} (cm)"};
    AxisSpec ptResAxis = {binning.axisPtRes, "(#it{p}_{T}^{rec} - #it{p}_{T}^{gen})/#it{p}_{T}^{gen}"};
    AxisSpec motherAxis = {kNFeedDownMothers, -0.5, kNFeedDownMothers - 0.5, "mother"};

    std::array<int, 8> nBinsEst = {100, 500, 102, 500, 102, 500, 102, 150};
    std::array<float, 8> lowEdgeEst = {-0.5, -0.5, -0.01, -0.5, -0.01, -0.5, -0.01, .0};
    std::array<float, 8> upEdgeEst = {99.5, 49999.5, 1.01, 499.5, 1.01, 499.5, 1.01, 150.0};

    ccdb->setURL(ccdbConf.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZ", "hVertexZ",
                        {HistType::kTH1D, {vertexZAxis}});
    rEventSelection.add("hEventsSelected", "hEventsSelected",
                        {HistType::kTH1D, {{15, 0, 15}}});

    rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "all");
    if (evSel.issel8) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "sel8");
    }

    rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "zvertex");

    if (evSel.isNoTimeFrameBorder) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "TFBorder");
    }
    if (evSel.isNoITSROFrameBorder) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "ITSROFBorder");
    }
    if (evSel.isVertexITSTPC) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "VertexITSTPC");
    }
    if (evSel.isNoSameBunchPileup) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "SameBunchPileup");
    }
    if (evSel.isGoodZvtxFT0vsPV) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "isGoodZvtxFT0vsPV");
    }
    if (evSel.isTriggerTVX) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "TVX");
    }
    if (evSel.isINELgt0) {
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "INEL>0");
    }
    // events without FV0/FT0 are rejected below, they need their own counter
    if (doprocessDataRun3LambdaK0s || doprocessRecMCLambdaK0s ||
        doprocessDataRun3Cascade || doprocessRecMCRun3Cascade ||
        doprocessFlatDistData || doprocessFlatDistMC) {
      nbinFlattenicity = nbin;
      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin++, "flattenicity");
    }

    rEventSelection.add("hFlattenicityDistribution", "hFlattenicityDistribution",
                        {HistType::kTH1D, {flatAxis}});
    rEventSelection.add("hCentFT0M", "hCentFT0M", {HistType::kTH1D, {centAxis}});
    rEventSelection.add("hCentFT0MvsFlattenicity", "hCentFT0MvsFlattenicity",
                        {HistType::kTH2D, {centAxis, flatAxis}});

    // 1-rho on the fine axis, the input to the percentile boundaries of
    // binning.axisFlat. Filled once per accepted collision, by the owner.
    rEventSelection.add("hFlatDistRec", "hFlatDistRec", {HistType::kTH1D, {flatFineAxis}});
    rEventSelection.add("hFlatDistRecINELgt0", "hFlatDistRecINELgt0", {HistType::kTH1D, {flatFineAxis}});
    rEventSelection.add("hFlatDistRecHM", "hFlatDistRecHM", {HistType::kTH1D, {flatFineAxis}});
    rEventSelection.add("hCentFT0MFine", "hCentFT0MFine", {HistType::kTH1D, {{1000, 0.0f, 100.0f, "FT0M percentile"}}});
    if (doprocessRecMCLambdaK0s || doprocessRecMCRun3Cascade || doprocessGenMC ||
        doprocessFlatDistMC) {
      // one entry per generated collision, the sample the closure denominator
      // lives in: binning.axisFlatTrue has to be set from this one, or the
      // numerator and the denominator do not share a class definition
      rEventSelection.add("hFlatDistGen", "hFlatDistGen", {HistType::kTH1D, {flatFineAxis}});
      rEventSelection.add("hFlatDistGenINELgt0", "hFlatDistGenINELgt0", {HistType::kTH1D, {flatFineAxis}});
      rEventSelection.add("hFlatDistGenHM", "hFlatDistGenHM", {HistType::kTH1D, {flatFineAxis}});
      // same, restricted to accepted reconstructed collisions: the ratio to
      // hFlatDistGen is the event selection bias on the class assignment
      rEventSelection.add("hFlatDistGenInRec", "hFlatDistGenInRec", {HistType::kTH1D, {flatFineAxis}});
      rEventSelection.add("hFlatDistGenInRecINELgt0", "hFlatDistGenInRecINELgt0", {HistType::kTH1D, {flatFineAxis}});
      rEventSelection.add("hFlatDistGenInRecHM", "hFlatDistGenInRecHM", {HistType::kTH1D, {flatFineAxis}});
    }
    if (doprocessRecMCLambdaK0s || doprocessRecMCRun3Cascade || doprocessGenMC ||
        doprocessFlatDistMC) {
      rEventSelection.add("hTrueFV0amplvsFlat", "TrueFV0MvsFlat", HistType::kTH2D,
                          {{500, -0.5, +499.5, "True Nch in FV0 region"}, flatTrueAxis});
    }
    if (doprocessRecMCLambdaK0s || doprocessRecMCRun3Cascade || doprocessFlatDistMC) {
      rEventSelection.add("hFlattenicityDistributionMCGen_Rec", "hFlattenicityDistributionMCGen_Rec",
                          {HistType::kTH1D, {flatTrueAxis}});
      rEventSelection.add("hFlattenicity_Corr_Gen_vs_Rec", "hFlattenicity_Corr_Gen_vs_Rec",
                          {HistType::kTH2D, {flatTrueAxis, flatAxis}});
      // migration of the class assignment, needed pT-differentially to unfold
      rEventSelection.add("hFlatGenVsRecFine", "hFlatGenVsRecFine",
                          {HistType::kTH2D, {flatFineAxis, flatFineAxis}});
    }

    // Charged particles in |eta| < cfgEtaChargedCut: <dNch/deta> per class is
    // the scale factor of Qpp, and the DCAxy templates give the secondary
    // contamination the same way the published analysis obtains it.
    if (eventClass.fillChargedQA) {
      rCharged.add("hNchVsFlat", "hNchVsFlat", {HistType::kTH2D, {nchAxis, flatAxis}});
      rCharged.add("hPtChVsFlat", "hPtChVsFlat", {HistType::kTH2D, {ptChargedAxis, flatAxis}});
      if (doprocessRecMCLambdaK0s || doprocessRecMCRun3Cascade || doprocessFlatDistMC) {
        rCharged.add("hPtChVsFlatRecPrim", "hPtChVsFlatRecPrim",
                     {HistType::kTH2D, {ptChargedAxis, flatAxis}});
        rCharged.add("hPtChVsFlatGenInRec", "hPtChVsFlatGenInRec",
                     {HistType::kTH2D, {ptChargedAxis, flatAxis}});
        rCharged.add("hNchVsFlatGenInRec", "hNchVsFlatGenInRec", {HistType::kTH2D, {nchAxis, flatAxis}});
        rCharged.add("hDcaXyPtChPrim", "hDcaXyPtChPrim",
                     {HistType::kTH3D, {dcaXyAxis, ptChargedAxis, flatAxis}});
        rCharged.add("hDcaXyPtChSec", "hDcaXyPtChSec",
                     {HistType::kTH3D, {dcaXyAxis, ptChargedAxis, flatAxis}});
      }
      if (doprocessDataRun3LambdaK0s || doprocessDataRun3Cascade || doprocessFlatDistData) {
        rCharged.add("hDcaXyPtCh", "hDcaXyPtCh",
                     {HistType::kTH3D, {dcaXyAxis, ptChargedAxis, flatAxis}});
      }
      if (doprocessGenMC) {
        rCharged.add("hPtChVsFlatGen", "hPtChVsFlatGen", {HistType::kTH2D, {ptChargedAxis, flatAxis}});
        rCharged.add("hNchVsFlatGen", "hNchVsFlatGen", {HistType::kTH2D, {nchAxis, flatAxis}});
        rCharged.add("hPtChVsTrueFlatGen", "hPtChVsTrueFlatGen", {HistType::kTH2D, {ptChargedAxis, flatTrueAxis}});
        rCharged.add("hNchVsTrueFlatGen", "hNchVsTrueFlatGen", {HistType::kTH2D, {nchAxis, flatTrueAxis}});
      }
    }

    if (doprocessDataRun3LambdaK0s || doprocessRecMCLambdaK0s) {
      // K0s reconstruction
      // Mass
      rKzeroShort.add("hMassK0s", "hMassK0s", {HistType::kTH1D, {k0sMassAxis}});
      rKzeroShort.add("hMassK0sSelected", "hMassK0sSelected",
                      {HistType::kTH1D, {k0sMassAxis}});

      // K0s topological/PID cuts
      rKzeroShort.add("hrapidityK0s", "hrapidityK0s",
                      {HistType::kTH1D, {{40, -2.0f, 2.0f, "y"}}});
      rKzeroShort.add("hctauK0s", "hctauK0s",
                      {HistType::kTH1D, {{40, 0.0f, 40.0f, "c#tau (cm)"}}});
      rKzeroShort.add("h2DdecayRadiusK0s", "h2DdecayRadiusK0s",
                      {HistType::kTH1D, {decayRadiusAxis}});
      rKzeroShort.add("hDCAV0DaughtersK0s", "hDCAV0DaughtersK0s",
                      {HistType::kTH1D, {{55, 0.0f, 2.2f, "DCA Daughters"}}});
      rKzeroShort.add("hV0CosPAK0s", "hV0CosPAK0s",
                      {HistType::kTH1D, {{100, 0.95f, 1.f, "CosPA"}}});
      rKzeroShort.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s",
                      {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rKzeroShort.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s",
                      {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rKzeroShort.add("hMassK0spT", "hMassK0spT",
                      {HistType::kTH2D, {{k0sMassAxis}, {ptK0sAxis}}});
      rKzeroShort.add("hMassK0spTFlat", "hMassK0spTFlat",
                      {HistType::kTH3D, {{k0sMassAxis}, {ptK0sAxis}, {flatAxis}}});
      rKzeroShort.add("hArmPodoAlphavsQTK0sAfterCut", "hArmPodoAlphavsQTK0sAfterCut",
                      {HistType::kTH2D, {{200, -1, 1, "#alpha"}, {70, 0, 0.35, "Q_{T}"}}});

      if (doprocessRecMCLambdaK0s) {
        rKzeroShort.add("Generated_MCRecoCollCheck_INEL_K0Short", "Generated_MCRecoCollCheck_INEL_K0Short",
                        {HistType::kTH2D, {{ptK0sAxis}, {flatAxis}}});
        rKzeroShort.add("Generated_MCRecoCollCheck_INELgt0_K0Short", "Generated_MCRecoCollCheck_INELgt0_K0Short",
                        {HistType::kTH2D, {{ptK0sAxis}, {flatAxis}}});
      }

      // Lambda reconstruction Mass
      rLambda.add("hMassLambda", "hMassLambda",
                  {HistType::kTH1D, {lambdaMassAxis}});
      rLambda.add("hMassLambdaSelected", "hMassLambdaSelected",
                  {HistType::kTH1D, {lambdaMassAxis}});

      // Lambda topological/PID cuts
      rLambda.add("hDCAV0DaughtersLambda", "hDCAV0DaughtersLambda",
                  {HistType::kTH1D, {{55, 0.0f, 2.2f, "DCA Daughters"}}});
      rLambda.add("hV0CosPALambda", "hV0CosPALambda",
                  {HistType::kTH1D, {{100, 0.95f, 1.f, "CosPA"}}});
      rLambda.add("hNSigmaPosProtonFromLambda", "hNSigmaPosProtonFromLambda",
                  {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rLambda.add("hNSigmaNegPionFromLambda", "hNSigmaNegPionFromLambda",
                  {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rLambda.add("hrapidityLambda", "hrapidityLambda",
                  {HistType::kTH1D, {{40, -2.0f, 2.0f, "y"}}});
      rLambda.add("hctauLambda", "hctauLambda",
                  {HistType::kTH1D, {{40, 0.0f, 40.0f, "c#tau (cm)"}}});
      rLambda.add("h2DdecayRadiusLambda", "h2DdecayRadiusLambda",
                  {HistType::kTH1D, {decayRadiusAxis}});
      rLambda.add("hMassLambdapT", "hMassLambdapT",
                  {HistType::kTH2D, {{lambdaMassAxis}, {ptLambdaAxis}}});
      rLambda.add("hMassLambdapTFlat", "hMassLambdapTFlat",
                  {HistType::kTH3D, {{lambdaMassAxis}, {ptLambdaAxis}, {flatAxis}}});
      if (doprocessRecMCLambdaK0s) {
        rLambda.add("Generated_MCRecoCollCheck_INEL_Lambda", "Generated_MCRecoCollCheck_INEL_Lambda",
                    {HistType::kTH2D, {{ptLambdaAxis}, {flatAxis}}});
        rLambda.add("Generated_MCRecoCollCheck_INELgt0_Lambda", "Generated_MCRecoCollCheck_INELgt0_Lambda",
                    {HistType::kTH2D, {{ptLambdaAxis}, {flatAxis}}});
        rLambda.add("hMassFeedDownLambdapTFlat", "hMassFeedDownLambdapTFlat",
                    {HistType::kTH3D, {{lambdaMassAxis}, {ptLambdaAxis}, {flatAxis}}});
        rLambda.add("hFeedDownLambdaPtVsMotherPt", "hFeedDownLambdaPtVsMotherPt",
                    {HistType::kTH2D, {{ptLambdaAxis}, {ptLambdaAxis}}});
        rLambda.add("hFeedDownLambdaMatrix", "hFeedDownLambdaMatrix",
                    {HistType::kTHnSparseF, {{ptLambdaAxis}, {ptLambdaAxis}, {flatAxis}, {motherAxis}}});
        rLambda.add("hFeedDownLambdaMotherPdg", "hFeedDownLambdaMotherPdg",
                    {HistType::kTH1D, {{motherAxis}}});
        rLambda.get<TH1>(HIST("hFeedDownLambdaMotherPdg"))->GetXaxis()->SetBinLabel(1, "#Xi^{-}");
        rLambda.get<TH1>(HIST("hFeedDownLambdaMotherPdg"))->GetXaxis()->SetBinLabel(2, "#Xi^{0}");
        rLambda.get<TH1>(HIST("hFeedDownLambdaMotherPdg"))->GetXaxis()->SetBinLabel(3, "#Omega^{-}");
        rLambda.get<TH1>(HIST("hFeedDownLambdaMotherPdg"))->GetXaxis()->SetBinLabel(4, "other");
        // generated mothers in the same events and classes as
        // hFeedDownLambdaMatrix, so the matrix is a probability per mother that
        // folds with the measured Xi. Xi0 is not measurable and only lives here.
        rLambda.add("hGenFeedDownMotherPtFlat", "hGenFeedDownMotherPtFlat",
                    {HistType::kTHnSparseF, {{ptLambdaAxis}, {flatAxis}, {motherAxis}}});
      }

      // AntiLambda reconstruction
      // Mass
      rAntiLambda.add("hMassAntiLambda", "hMassAntiLambda",
                      {HistType::kTH1D, {antilambdaMassAxis}});
      rAntiLambda.add("hMassAntiLambdaSelected", "hMassAntiLambdaSelected",
                      {HistType::kTH1D, {antilambdaMassAxis}});

      // AntiLambda topological/PID cuts
      rAntiLambda.add("hDCAV0DaughtersAntiLambda", "hDCAV0DaughtersAntiLambda",
                      {HistType::kTH1D, {{55, 0.0f, 2.2f, "DCA Daughters"}}});
      rAntiLambda.add("hV0CosPAAntiLambda", "hV0CosPAAntiLambda",
                      {HistType::kTH1D, {{100, 0.95f, 1.f, "CosPA"}}});
      rAntiLambda.add("hNSigmaPosPionFromAntiLambda",
                      "hNSigmaPosPionFromAntiLambda",
                      {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rAntiLambda.add("hNSigmaNegProtonFromAntiLambda",
                      "hNSigmaNegProtonFromAntiLambda",
                      {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rAntiLambda.add("hrapidityAntiLambda", "hrapidityAntiLambda",
                      {HistType::kTH1D, {{40, -2.0f, 2.0f, "y"}}});
      rAntiLambda.add("hctauAntiLambda", "hctauAntiLambda",
                      {HistType::kTH1D, {{40, 0.0f, 40.0f, "c#tau (cm)"}}});
      rAntiLambda.add("h2DdecayRadiusAntiLambda", "h2DdecayRadiusAntiLambda",
                      {HistType::kTH1D, {decayRadiusAxis}});
      rAntiLambda.add("hMassAntiLambdapT", "hMassAntiLambdapT",
                      {HistType::kTH2D, {{antilambdaMassAxis}, {ptLambdaAxis}}});
      rAntiLambda.add("hMassAntiLambdapTFlat", "hMassAntiLambdapTFlat",
                      {HistType::kTH3D, {{antilambdaMassAxis}, {ptLambdaAxis}, {flatAxis}}});
      if (doprocessRecMCLambdaK0s) {
        rAntiLambda.add("Generated_MCRecoCollCheck_INEL_AntiLambda", "Generated_MCRecoCollCheck_INEL_AntiLambda",
                        {HistType::kTH2D, {{ptLambdaAxis}, {flatAxis}}});
        rAntiLambda.add("Generated_MCRecoCollCheck_INELgt0_AntiLambda", "Generated_MCRecoCollCheck_INELgt0_AntiLambda",
                        {HistType::kTH2D, {{ptLambdaAxis}, {flatAxis}}});
        rAntiLambda.add("hMassFeedDownAntiLambdapTFlat", "hMassFeedDownAntiLambdapTFlat",
                        {HistType::kTH3D, {{antilambdaMassAxis}, {ptLambdaAxis}, {flatAxis}}});
        rAntiLambda.add("hFeedDownAntiLambdaPtVsMotherPt", "hFeedDownAntiLambdaPtVsMotherPt",
                        {HistType::kTH2D, {{ptLambdaAxis}, {ptLambdaAxis}}});
        rAntiLambda.add("hFeedDownAntiLambdaMatrix", "hFeedDownAntiLambdaMatrix",
                        {HistType::kTHnSparseF, {{ptLambdaAxis}, {ptLambdaAxis}, {flatAxis}, {motherAxis}}});
        rAntiLambda.add("hFeedDownAntiLambdaMotherPdg", "hFeedDownAntiLambdaMotherPdg",
                        {HistType::kTH1D, {{motherAxis}}});
        rAntiLambda.get<TH1>(HIST("hFeedDownAntiLambdaMotherPdg"))->GetXaxis()->SetBinLabel(1, "#Xi^{-}");
        rAntiLambda.get<TH1>(HIST("hFeedDownAntiLambdaMotherPdg"))->GetXaxis()->SetBinLabel(2, "#Xi^{0}");
        rAntiLambda.get<TH1>(HIST("hFeedDownAntiLambdaMotherPdg"))->GetXaxis()->SetBinLabel(3, "#Omega^{-}");
        rAntiLambda.get<TH1>(HIST("hFeedDownAntiLambdaMotherPdg"))->GetXaxis()->SetBinLabel(4, "other");
      }

      rCommonHist.add("hArmPodoAlphavsQT", "hArmPodoAlphavsQT",
                      {HistType::kTH2D, {{200, -1, 1, "#alpha"}, {70, 0, 0.35, "Q_{T}"}}});

      // DCA of the V0 to the PV: primaries peak at zero, feed-down does not, so
      // the data distribution can be fitted with the two MC templates
      rLambda.add("hDcaV0ToPVLambda", "hDcaV0ToPVLambda",
                  {HistType::kTH3D, {{dcaV0ToPvAxis}, {ptLambdaAxis}, {flatAxis}}});
      rAntiLambda.add("hDcaV0ToPVAntiLambda", "hDcaV0ToPVAntiLambda",
                      {HistType::kTH3D, {{dcaV0ToPvAxis}, {ptLambdaAxis}, {flatAxis}}});
      if (doprocessRecMCLambdaK0s) {
        rLambda.add("hDcaV0ToPVLambdaPrim", "hDcaV0ToPVLambdaPrim",
                    {HistType::kTH3D, {{dcaV0ToPvAxis}, {ptLambdaAxis}, {flatAxis}}});
        rLambda.add("hDcaV0ToPVLambdaSec", "hDcaV0ToPVLambdaSec",
                    {HistType::kTH3D, {{dcaV0ToPvAxis}, {ptLambdaAxis}, {flatAxis}}});
        rAntiLambda.add("hDcaV0ToPVAntiLambdaPrim", "hDcaV0ToPVAntiLambdaPrim",
                        {HistType::kTH3D, {{dcaV0ToPvAxis}, {ptLambdaAxis}, {flatAxis}}});
        rAntiLambda.add("hDcaV0ToPVAntiLambdaSec", "hDcaV0ToPVAntiLambdaSec",
                        {HistType::kTH3D, {{dcaV0ToPvAxis}, {ptLambdaAxis}, {flatAxis}}});

        rKzeroShort.add("hPtResK0s", "hPtResK0s", {HistType::kTH2D, {{ptK0sAxis}, {ptResAxis}}});
        rLambda.add("hPtResLambda", "hPtResLambda", {HistType::kTH2D, {{ptLambdaAxis}, {ptResAxis}}});
        rAntiLambda.add("hPtResAntiLambda", "hPtResAntiLambda", {HistType::kTH2D, {{ptLambdaAxis}, {ptResAxis}}});

        // numerator classified by the measured 1-rho, denominator by the true
        // one: the pair is the MC closure
        rKzeroShort.add("Generated_MCRecoCollCheck_INELgt0_K0Short_TrueFlat", "Generated_MCRecoCollCheck_INELgt0_K0Short_TrueFlat",
                        {HistType::kTH2D, {{ptK0sAxis}, {flatTrueAxis}}});
        rLambda.add("Generated_MCRecoCollCheck_INELgt0_Lambda_TrueFlat", "Generated_MCRecoCollCheck_INELgt0_Lambda_TrueFlat",
                    {HistType::kTH2D, {{ptLambdaAxis}, {flatTrueAxis}}});
        rAntiLambda.add("Generated_MCRecoCollCheck_INELgt0_AntiLambda_TrueFlat", "Generated_MCRecoCollCheck_INELgt0_AntiLambda_TrueFlat",
                        {HistType::kTH2D, {{ptLambdaAxis}, {flatTrueAxis}}});
        rKzeroShort.add("hMassK0spTTrueFlat", "hMassK0spTTrueFlat",
                        {HistType::kTH3D, {{k0sMassAxis}, {ptK0sAxis}, {flatTrueAxis}}});
        rLambda.add("hMassLambdapTTrueFlat", "hMassLambdapTTrueFlat",
                    {HistType::kTH3D, {{lambdaMassAxis}, {ptLambdaAxis}, {flatTrueAxis}}});
        rAntiLambda.add("hMassAntiLambdapTTrueFlat", "hMassAntiLambdapTTrueFlat",
                        {HistType::kTH3D, {{antilambdaMassAxis}, {ptLambdaAxis}, {flatTrueAxis}}});
      }
    }

    if (doprocessRecMCRun3Cascade || doprocessDataRun3Cascade) {
      rXi.add("hMassXi", "hMassXi", {HistType::kTH1D, {xiMassAxis}});
      rXi.add("hMassXiSelected", "hMassXiSelected",
              {HistType::kTH1D, {xiMassAxis}});

      // Xi topological/PID cuts
      rXi.add("hrapidityXi", "hrapidityXi",
              {HistType::kTH1D, {{40, -2.0f, 2.0f, "y"}}});
      rXi.add("hctauXi", "hctauXi",
              {HistType::kTH1D, {{40, 0.0f, 40.0f, "c#tau (cm)"}}});
      rXi.add("h2DdecayRadiusXi", "h2DdecayRadiusXi",
              {HistType::kTH1D, {decayRadiusAxis}});
      rXi.add("hDCAV0DaughtersXi", "hDCAV0DaughtersXi",
              {HistType::kTH1D, {{55, 0.0f, 2.2f, "DCA Daughters"}}});
      rXi.add("hV0CosPAXi", "hV0CosPAXi",
              {HistType::kTH1D, {{100, 0.95f, 1.f, "CosPA"}}});
      rXi.add("hNSigmaProtonFromXi", "hNSigmaProtonFromXi",
              {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rXi.add("hNSigmaPionFromXi", "hNSigmaPionFromXi",
              {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rXi.add("hNSigmaBachPionFromXi", "hNSigmaBachPionFromXi",
              {HistType::kTH2D, {{100, -5.f, 5.f, "n#sigma_{TPC}"}, {pTPCAxis}}});
      rXi.add("hMassXipT", "hMassXipT",
              {HistType::kTH2D, {{xiMassAxis}, {ptXiAxis}}});
      rXi.add("hMassXipTFlat", "hMassXipTFlat",
              {HistType::kTH3D, {{xiMassAxis}, {ptXiAxis}, {flatAxis}}});
      if (doprocessRecMCRun3Cascade) {
        rXi.add("Generated_MCRecoCollCheck_INEL_Xi", "Generated_MCRecoCollCheck_INEL_Xi",
                {HistType::kTH2D, {{ptXiAxis}, {flatAxis}}});
        rXi.add("Generated_MCRecoCollCheck_INELgt0_Xi", "Generated_MCRecoCollCheck_INELgt0_Xi",
                {HistType::kTH2D, {{ptXiAxis}, {flatAxis}}});
        rXi.add("hMassFeedDownXipTFlat", "hMassFeedDownXipTFlat",
                {HistType::kTH3D, {{xiMassAxis}, {ptXiAxis}, {flatAxis}}});
        rXi.add("hFeedDownXiPtVsMotherPt", "hFeedDownXiPtVsMotherPt",
                {HistType::kTH2D, {{ptXiAxis}, {ptXiAxis}}});
        rXi.add("hFeedDownXiMatrix", "hFeedDownXiMatrix",
                {HistType::kTHnSparseF, {{ptXiAxis}, {ptXiAxis}, {flatAxis}, {motherAxis}}});
        rXi.add("hFeedDownXiMotherPdg", "hFeedDownXiMotherPdg",
                {HistType::kTH1D, {{motherAxis}}});
        rXi.get<TH1>(HIST("hFeedDownXiMotherPdg"))->GetXaxis()->SetBinLabel(1, "#Xi^{-}");
        rXi.get<TH1>(HIST("hFeedDownXiMotherPdg"))->GetXaxis()->SetBinLabel(2, "#Xi^{0}");
        rXi.get<TH1>(HIST("hFeedDownXiMotherPdg"))->GetXaxis()->SetBinLabel(3, "#Omega^{-}");
        rXi.get<TH1>(HIST("hFeedDownXiMotherPdg"))->GetXaxis()->SetBinLabel(4, "other");

        rXi.add("hPtResXi", "hPtResXi", {HistType::kTH2D, {{ptXiAxis}, {ptResAxis}}});
        rXi.add("Generated_MCRecoCollCheck_INELgt0_Xi_TrueFlat", "Generated_MCRecoCollCheck_INELgt0_Xi_TrueFlat",
                {HistType::kTH2D, {{ptXiAxis}, {flatTrueAxis}}});
        rXi.add("hMassXipTTrueFlat", "hMassXipTTrueFlat",
                {HistType::kTH3D, {{xiMassAxis}, {ptXiAxis}, {flatTrueAxis}}});
      }
    }
    if (doprocessGenMC) {

      rEventSelection.get<TH1>(HIST("hEventsSelected"))->GetXaxis()->SetBinLabel(nbin, "Applied selection");

      rEventSelection.add("hVertexZGen", "hVertexZGen",
                          {HistType::kTH1D, {vertexZAxis}});

      rEventSelection.add("hFlattenicityDistributionMCGen", "hFlattenicityDistributionMCGen",
                          {HistType::kTH1D, {flatTrueAxis}});

      rEventSelection.add("hFlattenicityDistributionRecMCGen", "hFlattenicityDistributionRecMCGen",
                          {HistType::kTH1D, {flatAxis}});

      rEventSelection.add("hFlat_RecoColl_MC", "hFlat_RecoColl_MC", {HistType::kTH1D, {flatAxis}});
      rEventSelection.add("hFlat_RecoColl_MC_INELgt0", "hFlat_RecoColl_MC_INELgt0", {HistType::kTH1D, {flatAxis}});
      rEventSelection.add("hFlat_GenRecoColl_MC", "hFlat_GenRecoColl_MC", {HistType::kTH1D, {flatAxis}});
      rEventSelection.add("hFlat_GenRecoColl_MC_INELgt0", "hFlat_GenRecoColl_MC_INELgt0", {HistType::kTH1D, {flatAxis}});
      rEventSelection.add("hFlat_GenColl_MC", "hFlat_GenColl_MC", {HistType::kTH1D, {flatAxis}});
      rEventSelection.add("hFlat_GenColl_MC_INELgt0", "hFlat_GenColl_MC_INELgt0", {HistType::kTH1D, {flatAxis}});
      rEventSelection.add("hNEventsMCGen", "hNEventsMCGen", {HistType::kTH1D, {{4, 0.f, 4.f}}});
      rEventSelection.get<TH1>(HIST("hNEventsMCGen"))->GetXaxis()->SetBinLabel(1, "all");
      rEventSelection.get<TH1>(HIST("hNEventsMCGen"))->GetXaxis()->SetBinLabel(2, "zvertex_true");
      rEventSelection.get<TH1>(HIST("hNEventsMCGen"))->GetXaxis()->SetBinLabel(3, "INELgt0_true");
      rEventSelection.add("hNEventsMCGenReco", "hNEventsMCGenReco", {HistType::kTH1D, {{2, 0.f, 2.f}}});
      rEventSelection.get<TH1>(HIST("hNEventsMCGenReco"))->GetXaxis()->SetBinLabel(1, "INEL");
      rEventSelection.get<TH1>(HIST("hNEventsMCGenReco"))->GetXaxis()->SetBinLabel(2, "INELgt0");
      rEventSelection.add("hNEventsMCReco", "hNEventsMCReco", {HistType::kTH1D, {{3, 0.f, 3.f}}});
      rEventSelection.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(1, "all");
      rEventSelection.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(2, "pass ev sel");
      rEventSelection.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      rKzeroShort.add("pGen_MCGenRecoColl_INEL_K0Short", "pGen_MCGenRecoColl_INEL_K0Short",
                      {HistType::kTH2D, {ptK0sAxis, flatAxis}});
      rKzeroShort.add("Generated_MCRecoColl_INEL_K0Short", "Generated_MCRecoColl_INEL_K0Short",
                      {HistType::kTH2D, {ptK0sAxis, flatAxis}});
      rKzeroShort.add("pGen_MCGenColl_INEL_K0Short", "pGen_MCGenColl_INEL_K0Short",
                      {HistType::kTH2D, {ptK0sAxis, flatAxis}});
      rKzeroShort.add("pGen_MCGenRecoColl_INELgt0_K0Short", "pGen_MCGenRecoColl_INELgt0_K0Short",
                      {HistType::kTH2D, {ptK0sAxis, flatAxis}});
      rKzeroShort.add("Generated_MCRecoColl_INELgt0_K0Short", "Generated_MCRecoColl_INELgt0_K0Short",
                      {HistType::kTH2D, {ptK0sAxis, flatAxis}});
      rKzeroShort.add("pGen_MCGenColl_INELgt0_K0Short", "pGen_MCGenColl_INELgt0_K0Short",
                      {HistType::kTH2D, {ptK0sAxis, flatAxis}});

      rLambda.add("pGen_MCGenRecoColl_INEL_Lambda", "pGen_MCGenRecoColl_INEL_Lambda",
                  {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rLambda.add("Generated_MCRecoColl_INEL_Lambda", "Generated_MCRecoColl_INEL_Lambda",
                  {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rLambda.add("pGen_MCGenColl_INEL_Lambda", "pGen_MCGenColl_INEL_Lambda",
                  {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rLambda.add("pGen_MCGenRecoColl_INELgt0_Lambda", "pGen_MCGenRecoColl_INELgt0_Lambda",
                  {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rLambda.add("Generated_MCRecoColl_INELgt0_Lambda", "Generated_MCRecoColl_INELgt0_Lambda",
                  {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rLambda.add("pGen_MCGenColl_INELgt0_Lambda", "pGen_MCGenColl_INELgt0_Lambda",
                  {HistType::kTH2D, {ptLambdaAxis, flatAxis}});

      rAntiLambda.add("pGen_MCGenRecoColl_INEL_AntiLambda", "pGen_MCGenRecoColl_INEL_AntiLambda",
                      {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rAntiLambda.add("Generated_MCRecoColl_INEL_AntiLambda", "Generated_MCRecoColl_INEL_AntiLambda",
                      {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rAntiLambda.add("pGen_MCGenColl_INEL_AntiLambda", "pGen_MCGenColl_INEL_AntiLambda",
                      {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rAntiLambda.add("pGen_MCGenRecoColl_INELgt0_AntiLambda", "pGen_MCGenRecoColl_INELgt0_AntiLambda",
                      {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rAntiLambda.add("Generated_MCRecoColl_INELgt0_AntiLambda", "Generated_MCRecoColl_INELgt0_AntiLambda",
                      {HistType::kTH2D, {ptLambdaAxis, flatAxis}});
      rAntiLambda.add("pGen_MCGenColl_INELgt0_AntiLambda", "pGen_MCGenColl_INELgt0_AntiLambda",
                      {HistType::kTH2D, {ptLambdaAxis, flatAxis}});

      rXi.add("pGen_MCGenRecoColl_INEL_Xi", "pGen_MCGenRecoColl_INEL_Xi",
              {HistType::kTH2D, {ptXiAxis, flatAxis}});
      rXi.add("Generated_MCRecoColl_INEL_Xi", "Generated_MCRecoColl_INEL_Xi",
              {HistType::kTH2D, {ptXiAxis, flatAxis}});
      rXi.add("pGen_MCGenColl_INEL_Xi", "pGen_MCGenColl_INEL_Xi",
              {HistType::kTH2D, {ptXiAxis, flatAxis}});
      rXi.add("pGen_MCGenRecoColl_INELgt0_Xi", "pGen_MCGenRecoColl_INELgt0_Xi",
              {HistType::kTH2D, {ptXiAxis, flatAxis}});
      rXi.add("Generated_MCRecoColl_INELgt0_Xi", "Generated_MCRecoColl_INELgt0_Xi",
              {HistType::kTH2D, {ptXiAxis, flatAxis}});
      rXi.add("pGen_MCGenColl_INELgt0_Xi", "pGen_MCGenColl_INELgt0_Xi",
              {HistType::kTH2D, {ptXiAxis, flatAxis}});

      // The same INEL>0 counters against the true 1-rho. Together with the
      // reconstructed ones above they give the closure without a second job.
      rEventSelection.add("hFlat_RecoColl_MC_INELgt0_TrueFlat", "hFlat_RecoColl_MC_INELgt0_TrueFlat",
                          {HistType::kTH1D, {flatTrueAxis}});
      rEventSelection.add("hFlat_GenRecoColl_MC_INELgt0_TrueFlat", "hFlat_GenRecoColl_MC_INELgt0_TrueFlat",
                          {HistType::kTH1D, {flatTrueAxis}});
      rEventSelection.add("hFlat_GenColl_MC_INELgt0_TrueFlat", "hFlat_GenColl_MC_INELgt0_TrueFlat",
                          {HistType::kTH1D, {flatTrueAxis}});
      rKzeroShort.add("pGen_MCGenRecoColl_INELgt0_K0Short_TrueFlat", "pGen_MCGenRecoColl_INELgt0_K0Short_TrueFlat",
                      {HistType::kTH2D, {ptK0sAxis, flatTrueAxis}});
      rKzeroShort.add("Generated_MCRecoColl_INELgt0_K0Short_TrueFlat", "Generated_MCRecoColl_INELgt0_K0Short_TrueFlat",
                      {HistType::kTH2D, {ptK0sAxis, flatTrueAxis}});
      rKzeroShort.add("pGen_MCGenColl_INELgt0_K0Short_TrueFlat", "pGen_MCGenColl_INELgt0_K0Short_TrueFlat",
                      {HistType::kTH2D, {ptK0sAxis, flatTrueAxis}});
      rLambda.add("pGen_MCGenRecoColl_INELgt0_Lambda_TrueFlat", "pGen_MCGenRecoColl_INELgt0_Lambda_TrueFlat",
                  {HistType::kTH2D, {ptLambdaAxis, flatTrueAxis}});
      rLambda.add("Generated_MCRecoColl_INELgt0_Lambda_TrueFlat", "Generated_MCRecoColl_INELgt0_Lambda_TrueFlat",
                  {HistType::kTH2D, {ptLambdaAxis, flatTrueAxis}});
      rLambda.add("pGen_MCGenColl_INELgt0_Lambda_TrueFlat", "pGen_MCGenColl_INELgt0_Lambda_TrueFlat",
                  {HistType::kTH2D, {ptLambdaAxis, flatTrueAxis}});
      rAntiLambda.add("pGen_MCGenRecoColl_INELgt0_AntiLambda_TrueFlat", "pGen_MCGenRecoColl_INELgt0_AntiLambda_TrueFlat",
                      {HistType::kTH2D, {ptLambdaAxis, flatTrueAxis}});
      rAntiLambda.add("Generated_MCRecoColl_INELgt0_AntiLambda_TrueFlat", "Generated_MCRecoColl_INELgt0_AntiLambda_TrueFlat",
                      {HistType::kTH2D, {ptLambdaAxis, flatTrueAxis}});
      rAntiLambda.add("pGen_MCGenColl_INELgt0_AntiLambda_TrueFlat", "pGen_MCGenColl_INELgt0_AntiLambda_TrueFlat",
                      {HistType::kTH2D, {ptLambdaAxis, flatTrueAxis}});
      rXi.add("pGen_MCGenRecoColl_INELgt0_Xi_TrueFlat", "pGen_MCGenRecoColl_INELgt0_Xi_TrueFlat",
              {HistType::kTH2D, {ptXiAxis, flatTrueAxis}});
      rXi.add("Generated_MCRecoColl_INELgt0_Xi_TrueFlat", "Generated_MCRecoColl_INELgt0_Xi_TrueFlat",
              {HistType::kTH2D, {ptXiAxis, flatTrueAxis}});
      rXi.add("pGen_MCGenColl_INELgt0_Xi_TrueFlat", "pGen_MCGenColl_INELgt0_Xi_TrueFlat",
              {HistType::kTH2D, {ptXiAxis, flatTrueAxis}});
    }

    if (flatSel.flattenicityQA) {
      rFlattenicity.add("hEv", "Ev", HistType::kTH1D,
                        {{6, -0.5, 5.5, "index activated detector"}});
      rFlattenicity.add("hFV0amplRing1to4", "FV01to4", HistType::kTH1D,
                        {{4000, -0.5, +49999.5, "FV0 amplitude"}});
      rFlattenicity.add("hFT0Aampl", "FTAampl", HistType::kTH1D,
                        {{2000, -0.5, +199999.5, "Summed FT0A amplitude"}});
      rFlattenicity.add("hFT0Campl", "FTCampl", HistType::kTH1D,
                        {{2000, -0.5, +99999.5, "Summed FT0C amplitude"}});
      rFlattenicity.add("hFT0C", "FT0C", HistType::kTH1D,
                        {{2000, -0.5, 1999.5, "FT0C amplitude per channel"}});
      rFlattenicity.add("hFT0A", "FT0A", HistType::kTH1D,
                        {{2000, -0.5, 1999.5, "FT0A amplitude per channel"}});
      rFlattenicity.add("hFV0amplvsFlat", "FV0MvsFlat", HistType::kTH2D,
                        {{4000, -0.5, +49999.5, "FV0 amplitude"}, flatAxis});

      // estimators
      for (int iEe = 0; iEe < kNEstimators; ++iEe) {
        rFlattenicity.add(
          kHEst[iEe].data(), "", HistType::kTH2D,
          {{nBinsEst[iEe], lowEdgeEst[iEe], upEdgeEst[iEe], kTEst[iEe].data()},
           {100, -0.5, +99.5, "Global track"}});
      }

      // vs pT
      for (int iEe = 0; iEe < kNEstimators; ++iEe) {
        rFlattenicity.add(
          kHPtEst[iEe].data(), "", HistType::kTProfile,
          {{nBinsEst[iEe], lowEdgeEst[iEe], upEdgeEst[iEe], kTEst[iEe].data()}});
      }

      rFlattenicity.add("fMultFv0", "FV0 amp", HistType::kTH1D,
                        {{5000, -0.5, +199999.5, "FV0 amplitude"}});
      rFlattenicity.add(
        "hAmpV0VsCh", "", HistType::kTH2D,
        {{48, -0.5, 47.5, "channel"}, {500, -0.5, +19999.5, "FV0 amplitude"}});
      rFlattenicity.add(
        "hAmpV0VsChBeforeCalibration", "", HistType::kTH2D,
        {{48, -0.5, 47.5, "channel"}, {500, -0.5, +19999.5, "FV0 amplitude"}});

      rFlattenicity.add(
        "hAmpT0AVsChBeforeCalibration", "", HistType::kTH2D,
        {{24, -0.5, 23.5, "channel"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
      rFlattenicity.add(
        "hAmpT0CVsChBeforeCalibration", "", HistType::kTH2D,
        {{28, -0.5, 27.5, "channel"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

      rFlattenicity.add(
        "hAmpT0AVsCh", "", HistType::kTH2D,
        {{24, -0.5, 23.5, "channel"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
      rFlattenicity.add(
        "hAmpT0CVsCh", "", HistType::kTH2D,
        {{28, -0.5, 27.5, "channel"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

      rFlattenicity.add("hFlatFT0CvsFlatFT0A", "", HistType::kTH2D,
                        {{20, -0.01, +1.01, "flatenicity (FT0C)"},
                         {20, -0.01, +1.01, "flatenicity (FT0A)"}});
      rFlattenicity.add(
        "fEtaPhiFv0", "eta vs phi", HistType::kTH2D,
        {{8, 0.0, constants::math::TwoPI, "#phi (rad)"}, {5, 2.2, 5.1, "#eta"}});

      rFlattenicity.add("hAmpV0vsVtxBeforeCalibration", "", HistType::kTH2D,
                        {{30, -15.0, +15.0, "Trk mult"},
                         {1000, -0.5, +39999.5, "FV0 amplitude"}});
      rFlattenicity.add(
        "hAmpT0AvsVtxBeforeCalibration", "", HistType::kTH2D,
        {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
      rFlattenicity.add(
        "hAmpT0CvsVtxBeforeCalibration", "", HistType::kTH2D,
        {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

      rFlattenicity.add("hAmpV0vsVtx", "", HistType::kTH2D,
                        {{30, -15.0, +15.0, "Trk mult"},
                         {1000, -0.5, +39999.5, "FV0 amplitude"}});
      rFlattenicity.add(
        "hAmpT0AvsVtx", "", HistType::kTH2D,
        {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
      rFlattenicity.add(
        "hAmpT0CvsVtx", "", HistType::kTH2D,
        {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0C amplitude"}});
    }

    if ((doprocessDataRun3LambdaK0s || doprocessRecMCLambdaK0s) && (doprocessDataRun3Cascade || doprocessRecMCRun3Cascade)) {
      LOGF(fatal, "Can not run both LambdaK0s and Cascade process functions simulatenously. Try one at a time.");
    }

    if ((doprocessDataRun3LambdaK0s || doprocessDataRun3Cascade) && doprocessGenMC) {
      LOGF(fatal, "Can not run MCGen and Data process functions together. Try one of these at a time");
    }

    // the spectra passes come first, their event counts normalise the yields
    if (doprocessRecMCLambdaK0s) {
      eventHistOwner = kOwnerRecMCV0;
    } else if (doprocessDataRun3LambdaK0s) {
      eventHistOwner = kOwnerDataV0;
    } else if (doprocessRecMCRun3Cascade) {
      eventHistOwner = kOwnerRecMCCasc;
    } else if (doprocessDataRun3Cascade) {
      eventHistOwner = kOwnerDataCasc;
    } else if (doprocessFlatDistMC) {
      eventHistOwner = kOwnerFlatDistMC;
    } else if (doprocessFlatDistData) {
      eventHistOwner = kOwnerFlatDistData;
    }

    // tied to the owner so a new process function cannot bring the double fill back
    fillFlattenicityQAInGenMC = (eventHistOwner == kOwnerNone);

    // the generated 1-rho covers the FV0 only, the other estimators would
    // classify rec and gen with two different observables
    if ((doprocessRecMCLambdaK0s || doprocessRecMCRun3Cascade || doprocessGenMC ||
         doprocessFlatDistMC) &&
        flatSel.flattenicityforanalysis != kFlatFromFV0) {
      LOGF(fatal,
           "flattenicityforanalysis!=0 has no generator-level counterpart: "
           "estimateFlattenicityFV0MC covers the FV0 acceptance only");
    }

    // the estimator used for the analysis has to be computed
    if (flatSel.flattenicityforanalysis == kFlatFromFV0 && !flatSel.isflattenicitywithFV0 && !flatSel.isflattenicitywithFV0FT0C) {
      LOGF(fatal, "flattenicityforanalysis=0 (FV0) needs isflattenicitywithFV0 or isflattenicitywithFV0FT0C enabled");
    }
    if (flatSel.flattenicityforanalysis == kFlatFromFT0 && !flatSel.isflattenicitywithFT0) {
      LOGF(fatal, "flattenicityforanalysis=1 (FT0) needs isflattenicitywithFT0 enabled");
    }
    if (flatSel.flattenicityforanalysis == kFlatFromFV0FT0C && !(flatSel.isflattenicitywithFV0FT0C || (flatSel.isflattenicitywithFV0 && flatSel.isflattenicitywithFT0))) {
      LOGF(fatal, "flattenicityforanalysis=2 (FV0+FT0C) needs isflattenicitywithFV0FT0C enabled");
    }
  }

  // FT0 sectors group kNChannelsPerT0Sector consecutive channels
  int getT0ASector(int iCh)
  {
    const int iSec = iCh / kNChannelsPerT0Sector;
    return (iCh >= 0 && iSec < kNSectorsT0A) ? iSec : -1;
  }

  int getT0CSector(int iCh)
  {
    const int iSec = iCh / kNChannelsPerT0Sector;
    return (iCh >= 0 && iSec < kNSectorsT0C) ? iSec : -1;
  }

  // four inner rings of kNChannelsPerFV0Ring channels each, the rest is the outer ring
  int getFV0Ring(int iCh)
  {
    const int iRing = iCh / kNChannelsPerFV0Ring;
    if (iRing < 0) {
      return 0;
    }
    return (iRing < kOuterFV0RingIndex) ? iRing : kOuterFV0RingIndex;
  }

  // FV0 channel -> phi-ordered index, a fixed permutation of the 48 channels
  int getFV0IndexPhi(int iCh)
  {
    if (iCh < 0 || iCh >= kNCells) {
      return -1;
    }
    return kFV0PhiIndex[iCh];
  }

  float getFlatenicity(std::span<float> signals)
  {
    int entries = signals.size();
    float flat = 9999;
    float mRho = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      mRho += 1.0 * signals[iCell];
    }
    // average activity per cell
    mRho /= (1.0 * entries);
    // get sigma
    float sRhoTmp = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      sRhoTmp += std::pow(1.0 * signals[iCell] - mRho, 2);
    }
    sRhoTmp /= (1.0 * entries * entries);
    float sRho = std::sqrt(sRhoTmp);
    if (mRho > 0) {
      flat = sRho / mRho;
    }
    return flat;
  }
  // V0A signal and flatenicity calculation

  static constexpr int kNeta5 = 2; // FT0C + FT0A
  static constexpr std::array<float, kNeta5> kWeigthsEta5 = {0.0490638, 0.010958415};
  static constexpr std::array<float, kNeta5> kDeltaEeta5 = {1.1, 1.2};

  static constexpr int kNeta6 = 2; // FT0C + FV0
  static constexpr std::array<float, kNeta6> kWeigthsEta6 = {0.0490638, 0.00353962};
  static constexpr std::array<float, kNeta6> kDeltaEeta6 = {1.1, 2.9};

  static constexpr int kInnerFV0 = 32;
  static constexpr float kMaxEtaFV0 = 5.1;
  static constexpr float kMinEtaFV0 = 2.2;
  static constexpr float kDetaFV0 = (kMaxEtaFV0 - kMinEtaFV0) / 5.0;

  // no FV0/FT0 information
  static constexpr float kInvalidFlattenicity = -1.f;

  static constexpr int kNCells = 48; // 48 sectors in FV0
  static constexpr std::array<int, kNCells> kFV0PhiIndex = {
    0, 1, 2, 3, 7, 6, 5, 4,
    8, 9, 10, 11, 15, 14, 13, 12,
    16, 17, 18, 19, 23, 22, 21, 20,
    24, 25, 26, 27, 31, 30, 29, 28,
    32, 34, 36, 38, 47, 45, 43, 41,
    33, 35, 37, 39, 46, 44, 42, 40};
  std::array<float, kNCells> rhoLattice{};
  std::array<float, kNCells> rhoLatticeFV0AMC{};
  std::array<float, kNCells> ampchannel{};
  std::array<float, kNCells> ampchannelBefore{};
  static constexpr int kNCellsT0A = 24;
  std::array<float, kNCellsT0A> rhoLatticeT0A{};
  static constexpr int kNCellsT0C = 28;
  std::array<float, kNCellsT0C> rhoLatticeT0C{};

  std::array<float, kNEstimators> estimator{};

  // fillCounter=false runs the same cuts without filling hEventsSelected
  template <typename TCollision>
  bool isEventSelected(TCollision const& collision, bool fillCounter = true)
  {
    float nbinev = 0.5;
    if (fillCounter) {
      rEventSelection.fill(HIST("hEventsSelected"), nbinev);
    }

    if (evSel.issel8 && !collision.sel8()) {
      return false;
    }
    if (evSel.issel8) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }

    if (std::abs(collision.posZ()) > evSel.cutzvertex) {
      return false;
    }

    nbinev++;
    if (fillCounter) {
      rEventSelection.fill(HIST("hEventsSelected"), nbinev);
    }

    if (evSel.isNoTimeFrameBorder &&
        !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (evSel.isNoTimeFrameBorder) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }

    if (evSel.isNoITSROFrameBorder &&
        !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (evSel.isNoITSROFrameBorder) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }
    if (evSel.isVertexITSTPC &&
        !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (evSel.isVertexITSTPC) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }

    if (evSel.isNoSameBunchPileup &&
        !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (evSel.isNoSameBunchPileup) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }

    if (evSel.isGoodZvtxFT0vsPV &&
        !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (evSel.isGoodZvtxFT0vsPV) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }
    if (evSel.isTriggerTVX &&
        !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (evSel.isTriggerTVX) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }

    if (eventClass.applyCentSel && (collision.centFT0M() < eventClass.cfgCentMin || collision.centFT0M() > eventClass.cfgCentMax)) {
      return false;
    }

    if (evSel.isINELgt0 && (collision.isInelGt0() == false)) {
      return false;
    }
    if (evSel.isINELgt0) {
      nbinev++;
      if (fillCounter) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinev);
      }
    }

    return true;
  }

  // Daughter track acceptance and quality. The quality cuts are off by default
  // and exist so a systematic variation can move them from the JSON.
  template <typename TTrack>
  bool isSelectedDaughterTrack(TTrack const& track, float minCrossedRows)
  {
    if (std::abs(track.eta()) > trkPid.cfgTrkEtaCut || track.pt() < trkPid.cfgTrkLowPtCut) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < minCrossedRows) {
      return false;
    }
    if (trkPid.minCrossedRowsOverFindable > 0.f &&
        track.tpcCrossedRowsOverFindableCls() < trkPid.minCrossedRowsOverFindable) {
      return false;
    }
    if (trkPid.maxTpcChi2NCl > 0.f && track.tpcChi2NCl() > trkPid.maxTpcChi2NCl) {
      return false;
    }
    if (trkPid.maxItsChi2NCl > 0.f && track.itsChi2NCl() > trkPid.maxItsChi2NCl) {
      return false;
    }
    if (trkPid.minItsNCls > 0 && static_cast<int>(track.itsNCls()) < trkPid.minItsNCls) {
      return false;
    }
    return true;
  }

  // 1-rho on the fine axis for the percentile boundaries, inclusive, INEL>0 and
  // in the top FT0M class
  void fillFlatDistRec(float flattenicity, float centFT0M, bool isInelGt0)
  {
    rEventSelection.fill(HIST("hFlatDistRec"), flattenicity);
    rEventSelection.fill(HIST("hCentFT0MFine"), centFT0M);
    if (isInelGt0) {
      rEventSelection.fill(HIST("hFlatDistRecINELgt0"), flattenicity);
      if (centFT0M < eventClass.cfgCentMaxHM) {
        rEventSelection.fill(HIST("hFlatDistRecHM"), flattenicity);
      }
    }
  }

  // processGenMC only. A negative centFT0M means no accepted reconstructed
  // counterpart, so the collision has no FT0M class and stays out of the HM one.
  void fillFlatDistGen(float flattenicityGen, float centFT0M, bool isInelGt0)
  {
    rEventSelection.fill(HIST("hFlatDistGen"), flattenicityGen);
    if (isInelGt0) {
      rEventSelection.fill(HIST("hFlatDistGenINELgt0"), flattenicityGen);
      if (centFT0M >= 0.f && centFT0M < eventClass.cfgCentMaxHM) {
        rEventSelection.fill(HIST("hFlatDistGenHM"), flattenicityGen);
      }
    }
  }

  // the same from the rec-level process functions, diagnostic only
  void fillFlatDistGenInRec(float flattenicityGen, float centFT0M, bool isInelGt0)
  {
    rEventSelection.fill(HIST("hFlatDistGenInRec"), flattenicityGen);
    if (isInelGt0) {
      rEventSelection.fill(HIST("hFlatDistGenInRecINELgt0"), flattenicityGen);
      if (centFT0M >= 0.f && centFT0M < eventClass.cfgCentMaxHM) {
        rEventSelection.fill(HIST("hFlatDistGenInRecHM"), flattenicityGen);
      }
    }
  }

  // Charged tracks in the tracking acceptance: <dNch/deta> per class is the
  // scale factor of Qpp, and the DCAxy distribution carries the secondary
  // contamination.
  template <bool isMC, typename TTracks>
  int fillChargedRec(TTracks const& tracks, float flattenicity)
  {
    if (!eventClass.fillChargedQA) {
      return 0;
    }
    int nch = 0;
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack() || std::abs(track.eta()) > eventClass.cfgEtaChargedCut) {
        continue;
      }
      nch++;
      rCharged.fill(HIST("hPtChVsFlat"), track.pt(), flattenicity);
      if constexpr (isMC) {
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto& mcParticle = track.mcParticle();
        if (mcParticle.isPhysicalPrimary()) {
          rCharged.fill(HIST("hPtChVsFlatRecPrim"), track.pt(), flattenicity);
          rCharged.fill(HIST("hDcaXyPtChPrim"), track.dcaXY(), track.pt(), flattenicity);
        } else {
          rCharged.fill(HIST("hDcaXyPtChSec"), track.dcaXY(), track.pt(), flattenicity);
        }
      } else {
        rCharged.fill(HIST("hDcaXyPtCh"), track.dcaXY(), track.pt(), flattenicity);
      }
    }
    rCharged.fill(HIST("hNchVsFlat"), nch, flattenicity);
    return nch;
  }

  template <typename TMcParticles>
  int fillChargedGen(TMcParticles const& mcParticles, float flattenicity, bool trueFlat)
  {
    if (!eventClass.fillChargedQA) {
      return 0;
    }
    int nch = 0;
    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.eta()) > eventClass.cfgEtaChargedCut) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
      if (!(pdgParticle && std::abs(pdgParticle->Charge()) > kMinCharge)) {
        continue;
      }
      nch++;
      if (trueFlat) {
        rCharged.fill(HIST("hPtChVsTrueFlatGen"), mcParticle.pt(), flattenicity);
      } else {
        rCharged.fill(HIST("hPtChVsFlatGen"), mcParticle.pt(), flattenicity);
      }
    }
    if (trueFlat) {
      rCharged.fill(HIST("hNchVsTrueFlatGen"), nch, flattenicity);
    } else {
      rCharged.fill(HIST("hNchVsFlatGen"), nch, flattenicity);
    }
    return nch;
  }

  // ================= Calibration objects from CCDB ==================== //
  // A missing or malformed object falls back to unity, so a run without
  // calibration is left uncorrected instead of being equalized with somebody
  // else's constants.

  std::vector<float> fetchGainEq(const std::string& path, int run, std::size_t nChannels)
  {
    const auto* obj = ccdb->getForRun<std::vector<float>>(path, run);
    if (!obj || obj->size() != nChannels) {
      LOGF(warning, "No gain equalization of size %zu at %s for run %d, using unity",
           nChannels, path.c_str(), run);
      return std::vector<float>(nChannels, 1.f);
    }
    return *obj;
  }

  TProfile2D* fetchVtxEq(const std::string& path, int run)
  {
    auto* obj = ccdb->getForRun<TProfile2D>(path, run);
    if (!obj) {
      LOGF(warning, "No z-vertex equalization at %s for run %d, lattice left uncorrected",
           path.c_str(), run);
    }
    return obj;
  }

  // called per collision, does nothing unless a correction is switched on
  template <typename TBC>
  void initCcdb(TBC const& bc)
  {
    const int run = bc.runNumber();
    if (run == mRunNumber) {
      return;
    }
    mRunNumber = run;
    if (flatSel.applyCalibCh) {
      gainFV0 = fetchGainEq(ccdbConf.gainEqPath.value + "/FV0", run, kNCells);
      gainFT0A = fetchGainEq(ccdbConf.gainEqPath.value + "/FT0A", run, kNCellsT0A);
      gainFT0C = fetchGainEq(ccdbConf.gainEqPath.value + "/FT0C", run, kNCellsT0C);
    }
    if (flatSel.applyCalibVtx) {
      vtxEqFV0 = fetchVtxEq(ccdbConf.vtxEqPath.value + "/FV0", run);
      vtxEqFT0A = fetchVtxEq(ccdbConf.vtxEqPath.value + "/FT0A", run);
      vtxEqFT0C = fetchVtxEq(ccdbConf.vtxEqPath.value + "/FT0C", run);
    }
  }

  // the gain is divided out, the convention of flattenictyPikp
  static float gainEqFactor(std::vector<float> const& gain, std::size_t channel)
  {
    if (channel >= gain.size() || !(gain[channel] > 0.f)) {
      return 1.f;
    }
    return gain[channel];
  }

  // 1 for a missing map or an empty bin, so a bad bin never zeroes a cell
  static float vtxEqFactor(TProfile2D const* map, int cell, float vtxZ)
  {
    if (!map) {
      return 1.f;
    }
    const float factor = map->GetBinContent(map->GetXaxis()->FindBin(cell),
                                            map->GetYaxis()->FindBin(vtxZ));
    return (factor > 0.f) ? factor : 1.f;
  }

  // ============== Flattenicity estimation begins  ===================== //
  // fillQA=false skips the QA registry, for callers that would fill it twice
  template <typename TCollision, typename Tracks>
  float estimateFlattenicity(TCollision const& collision, Tracks const& tracks, bool fillQA = true)
  {
    const bool flattenicityQAhere = flatSel.flattenicityQA && fillQA;
    if (flatSel.applyCalibCh || flatSel.applyCalibVtx) {
      initCcdb(collision.template bc_as<soa::Join<aod::BCs, aod::Timestamps>>());
    }
    std::array<float, kNeta5> ampl5 = {0, 0};
    std::array<float, kNeta6> ampl6 = {0, 0};

    auto vtxZ = collision.posZ();

    float sumAmpFV0 = 0;
    float sumAmpFV01to4Ch = 0;
    // gain-equalized but not yet vertex-equalized, the input to the z-vertex QA
    float sumAmpFV0BeforeVtx = 0;

    ampchannel.fill(0.0);
    ampchannelBefore.fill(0.0);
    rhoLattice.fill(0);

    if ((flatSel.isflattenicitywithFV0 || flatSel.isflattenicitywithFV0FT0C) &&
        collision.has_foundFV0()) {

      auto fv0 = collision.foundFV0();
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
        float phiv0 = -999.0;
        float etav0 = -999.0;
        int channelv0 = fv0.channel()[ich];
        float amplCh = fv0.amplitude()[ich];
        int ringindex = getFV0Ring(channelv0);
        int channelv0phi = getFV0IndexPhi(channelv0);
        etav0 = kMaxEtaFV0 - (kDetaFV0 / 2.0) * (2.0 * ringindex + 1);
        if (channelv0 < kInnerFV0) {
          phiv0 = (2.0 * (channelv0phi - 8 * ringindex) + 1) * constants::math::PI / (8.0);
        } else {
          phiv0 = ((2.0 * channelv0phi) + 1 - 64.0) * constants::math::TwoPI / (32.0);
        }
        ampchannelBefore[channelv0phi] = amplCh;
        if (flatSel.applyCalibCh) {
          amplCh /= gainEqFactor(gainFV0, channelv0);
        }
        sumAmpFV0BeforeVtx += amplCh;
        if (flatSel.applyCalibVtx) {
          amplCh *= vtxEqFactor(vtxEqFV0, channelv0phi, vtxZ);
        }
        sumAmpFV0 += amplCh;

        if (channelv0 >= kNChannelsPerFV0Ring) { // exclude the 1st ring, eta 2.2,4.52
          sumAmpFV01to4Ch += amplCh;
        }
        if (flattenicityQAhere) {
          rFlattenicity.fill(HIST("fEtaPhiFv0"), phiv0, etav0, amplCh);
        }
        ampchannel[channelv0phi] = amplCh;
        if (channelv0 < kInnerFV0) {
          rhoLattice[channelv0phi] = amplCh;
        } else {
          rhoLattice[channelv0phi] = amplCh / 2.0; // two channels per bin
        }
      }

      if (flattenicityQAhere) {
        rFlattenicity.fill(HIST("hAmpV0vsVtxBeforeCalibration"), vtxZ, sumAmpFV0BeforeVtx);
        rFlattenicity.fill(HIST("hAmpV0vsVtx"), vtxZ, sumAmpFV0);
      }
    }

    float flattenicityfv0 = 9999;
    if (flatSel.isflattenicitywithFV0 || flatSel.isflattenicitywithFV0FT0C) {
      flattenicityfv0 = getFlatenicity({rhoLattice.data(), rhoLattice.size()});
    }

    // global tracks
    float ptT = 0.;
    int multGlob = 0;
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
      if (track.pt() > ptT) {
        ptT = track.pt();
      }
      multGlob++;
    }

    // FT0
    float sumAmpFT0A = 0.f;
    float sumAmpFT0C = 0.f;
    float sumAmpFT0ABeforeVtx = 0.f;
    float sumAmpFT0CBeforeVtx = 0.f;

    rhoLatticeT0A.fill(0);
    rhoLatticeT0C.fill(0);

    if ((flatSel.isflattenicitywithFT0 || flatSel.isflattenicitywithFV0FT0C) &&
        collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      if (flatSel.isflattenicitywithFT0) {
        for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
          float amplitude = ft0.amplitudeA()[i_a];
          uint8_t channel = ft0.channelA()[i_a];
          int sector = getT0ASector(channel);
          float amplitudeBeforeVtx = amplitude;
          if (sector >= 0 && sector < kNCellsT0A) {
            if (flattenicityQAhere) {
              rFlattenicity.fill(HIST("hAmpT0AVsChBeforeCalibration"), sector,
                                 amplitude);
            }
            if (flatSel.applyCalibCh) {
              amplitude /= gainEqFactor(gainFT0A, sector);
            }
            if (flattenicityQAhere) {
              rFlattenicity.fill(HIST("hAmpT0AVsCh"), sector, amplitude);
            }
            amplitudeBeforeVtx = amplitude;
            if (flatSel.applyCalibVtx) {
              amplitude *= vtxEqFactor(vtxEqFT0A, sector, vtxZ);
            }
            rhoLatticeT0A[sector] += amplitude;
          }
          sumAmpFT0A += amplitude;
          sumAmpFT0ABeforeVtx += amplitudeBeforeVtx;
          if (flattenicityQAhere) {
            rFlattenicity.fill(HIST("hFT0A"), amplitude);
          }
        }
      }

      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        uint8_t channel = ft0.channelC()[i_c];
        int sector = getT0CSector(channel);
        float amplitudeBeforeVtx = amplitude;
        if (sector >= 0 && sector < kNCellsT0C) {
          if (flattenicityQAhere) {
            rFlattenicity.fill(HIST("hAmpT0CVsChBeforeCalibration"), sector,
                               amplitude);
          }
          if (flatSel.applyCalibCh) {
            amplitude /= gainEqFactor(gainFT0C, sector);
          }
          if (flattenicityQAhere) {
            rFlattenicity.fill(HIST("hAmpT0CVsCh"), sector, amplitude);
          }
          amplitudeBeforeVtx = amplitude;
          if (flatSel.applyCalibVtx) {
            amplitude *= vtxEqFactor(vtxEqFT0C, sector, vtxZ);
          }
          rhoLatticeT0C[sector] += amplitude;
        }
        sumAmpFT0C += amplitude;
        sumAmpFT0CBeforeVtx += amplitudeBeforeVtx;
        if (flattenicityQAhere) {
          rFlattenicity.fill(HIST("hFT0C"), amplitude);
        }
      }
      if (flattenicityQAhere) {
        rFlattenicity.fill(HIST("hAmpT0AvsVtxBeforeCalibration"), vtxZ,
                           sumAmpFT0ABeforeVtx);
        rFlattenicity.fill(HIST("hAmpT0CvsVtxBeforeCalibration"), vtxZ,
                           sumAmpFT0CBeforeVtx);
        rFlattenicity.fill(HIST("hAmpT0AvsVtx"), vtxZ, sumAmpFT0A);
        rFlattenicity.fill(HIST("hAmpT0CvsVtx"), vtxZ, sumAmpFT0C);
      }
    }
    float flatenicityT0a = 9999;
    if (flatSel.isflattenicitywithFT0) {
      flatenicityT0a =
        getFlatenicity({rhoLatticeT0A.data(), rhoLatticeT0A.size()});
    }
    float flatenicityT0c = 9999;
    if (flatSel.isflattenicitywithFT0 || flatSel.isflattenicitywithFV0FT0C) {
      flatenicityT0c =
        getFlatenicity({rhoLatticeT0C.data(), rhoLatticeT0C.size()});
    }

    bool isOKEstimator5 = false;
    bool isOKEstimator6 = false;
    float combinedEstimator5 = 0;
    float combinedEstimator6 = 0;

    for (int iEe = 0; iEe < kNEstimators; ++iEe) {
      estimator[iEe] = 0;
    }

    if (!collision.has_foundFV0() || !collision.has_foundFT0()) {
      // no FV0/FT0, flattenicity undefined
      return kInvalidFlattenicity;
    }

    float allWeights = 0;
    // option 5
    ampl5[0] = sumAmpFT0C;
    ampl5[1] = sumAmpFT0A;
    if (sumAmpFT0C > 0 && sumAmpFT0A > 0) {
      isOKEstimator5 = true;
    }
    if (isOKEstimator5) {
      if (flatSel.applyNorm) {
        allWeights = 0;
        for (int i5 = 0; i5 < kNeta5; ++i5) {
          combinedEstimator5 +=
            ampl5[i5] * kWeigthsEta5[i5] / kDeltaEeta5[i5];
          allWeights += kWeigthsEta5[i5];
        }
        combinedEstimator5 /= allWeights;
      } else {
        for (int i5 = 0; i5 < kNeta5; ++i5) {
          combinedEstimator5 += ampl5[i5] * kWeigthsEta5[i5];
        }
      }
    }
    // option 6: FT0C + FV0
    ampl6[0] = sumAmpFT0C;
    ampl6[1] = sumAmpFV0;
    if (sumAmpFT0C > 0 && sumAmpFV0 > 0) {
      isOKEstimator6 = true;
    }
    if (isOKEstimator6) {
      if (flatSel.applyNorm) {
        allWeights = 0;
        for (int i6 = 0; i6 < kNeta6; ++i6) {
          combinedEstimator6 +=
            ampl6[i6] * kWeigthsEta6[i6] / kDeltaEeta6[i6];
          allWeights += kWeigthsEta6[i6];
        }
        combinedEstimator6 /= allWeights;
      } else {
        for (int i6 = 0; i6 < kNeta6; ++i6) {
          combinedEstimator6 += ampl6[i6] * kWeigthsEta6[i6];
        }
      }
    }
    if (flattenicityQAhere) {
      rFlattenicity.fill(HIST("hFT0Aampl"), sumAmpFT0A);
      rFlattenicity.fill(HIST("hFT0Campl"), sumAmpFT0C);
      rFlattenicity.fill(HIST("hFV0amplRing1to4"), sumAmpFV01to4Ch);
      rFlattenicity.fill(HIST("hEv"), 4);
    }
    estimator[0] = multGlob;
    estimator[1] = sumAmpFV0;
    estimator[2] = 1.0 - flattenicityfv0;
    estimator[3] = combinedEstimator5;
    float flatenicityFT0 = (flatenicityT0a + flatenicityT0c) / 2.0;
    estimator[4] = 1.0 - flatenicityFT0;
    estimator[5] = combinedEstimator6;
    float flatenicityFT0v0 = 0.5 * flattenicityfv0 + 0.5 * flatenicityT0c;
    estimator[6] = 1.0 - flatenicityFT0v0;
    estimator[7] = ptT;
    if (flattenicityQAhere) {
      rFlattenicity.fill(HIST(kHEst[0]), estimator[0], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[1]), estimator[1], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[2]), estimator[2], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[3]), estimator[3], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[4]), estimator[4], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[5]), estimator[5], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[6]), estimator[6], estimator[0]);
      rFlattenicity.fill(HIST(kHEst[7]), estimator[7], estimator[0]);

      // plot pt vs estimators
      for (const auto& track : tracks) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        float pt = track.pt();
        rFlattenicity.fill(HIST(kHPtEst[0]), estimator[0], pt);
        rFlattenicity.fill(HIST(kHPtEst[1]), estimator[1], pt);
        rFlattenicity.fill(HIST(kHPtEst[2]), estimator[2], pt);
        rFlattenicity.fill(HIST(kHPtEst[3]), estimator[3], pt);
        rFlattenicity.fill(HIST(kHPtEst[4]), estimator[4], pt);
        rFlattenicity.fill(HIST(kHPtEst[5]), estimator[5], pt);
        rFlattenicity.fill(HIST(kHPtEst[6]), estimator[6], pt);
        rFlattenicity.fill(HIST(kHPtEst[7]), estimator[7], pt);
      }

      if (flatSel.isflattenicitywithFV0) {
        for (int iCh = 0; iCh < kNCells; ++iCh) {
          rFlattenicity.fill(HIST("hAmpV0VsCh"), iCh, ampchannel[iCh]);
          rFlattenicity.fill(HIST("hAmpV0VsChBeforeCalibration"), iCh,
                             ampchannelBefore[iCh]);
        }
      }

      rFlattenicity.fill(HIST("fMultFv0"), sumAmpFV0);
      rFlattenicity.fill(HIST("hFlatFT0CvsFlatFT0A"), flatenicityT0c,
                         flatenicityT0a);
    }
    float finalflattenicity = estimator[2];
    if (flattenicityQAhere) {
      rFlattenicity.fill(HIST("hFV0amplvsFlat"), sumAmpFV0, estimator[2]);
    }

    if (flatSel.flattenicityforanalysis == kFlatFromFT0) {
      finalflattenicity = estimator[4];
    }
    if (flatSel.flattenicityforanalysis == kFlatFromFV0FT0C) {
      finalflattenicity = estimator[6];
    }
    return finalflattenicity;
  }

  // which weak decay produced a non-primary candidate, bins of hFeedDown*MotherPdg
  static constexpr int kFdXiMinus = 0;
  static constexpr int kFdXiZero = 1;
  static constexpr int kFdOmegaMinus = 2;
  static constexpr int kFdOther = 3;
  static constexpr int kNFeedDownMothers = 4;

  // -1 for a species that does not feed the measured hadrons
  static int feedDownSpecies(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kXiMinus:
        return kFdXiMinus;
      case o2::constants::physics::Pdg::kXi0:
        return kFdXiZero;
      case PDG_t::kOmegaMinus:
        return kFdOmegaMinus;
      default:
        return -1;
    }
  }

  template <typename TMcParticle>
  int getFeedDownMother(TMcParticle const& mcParticle, float& motherPt)
  {
    motherPt = -1.f;
    if (!mcParticle.has_mothers()) {
      return kFdOther;
    }
    for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
      motherPt = mother.pt();
      const int species = feedDownSpecies(mother.pdgCode());
      return (species >= 0) ? species : kFdOther;
    }
    return kFdOther;
  }

  // fillQA=false skips hTrueFV0amplvsFlat, for callers that would fill it twice
  template <typename McParticles>
  float estimateFlattenicityFV0MC(McParticles const& mcParticles, bool fillQA = true)
  {
    rhoLatticeFV0AMC.fill(0);
    int multFV0 = 0;

    for (const auto& mcParticle : mcParticles) {
      // the measured 1-rho also sees secondaries; dropping the primary
      // requirement is the handle on the leading MC non-closure source
      if (eventClass.genFlatPrimariesOnly && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (!(mcParticle.pt() > 0)) {
        continue;
      }

      auto pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
      if (!(pdgParticle && std::abs(pdgParticle->Charge()) > kMinCharge)) {
        continue;
      }

      const float etap = mcParticle.eta();
      const float phip = mcParticle.phi();
      int isegment = 0;

      for (int ieta = 0; ieta < kNFV0EtaRings; ieta++) {
        // the outer ring ends exactly at kMinEtaFV0, the rest follow from the ring width
        const float etamax = kMaxEtaFV0 - ieta * kDetaFV0;
        const float etamin = (ieta == kOuterFV0RingIndex) ? kMinEtaFV0
                                                          : kMaxEtaFV0 - (ieta + 1) * kDetaFV0;
        const int nsectors = (ieta == kOuterFV0RingIndex) ? kNSectorsFV0OuterRing
                                                          : kNChannelsPerFV0Ring;
        for (int iphi = 0; iphi < nsectors; iphi++) {
          const float minphi = iphi * constants::math::TwoPI / nsectors;
          const float maxphi = (iphi + 1) * constants::math::TwoPI / nsectors;
          const float dphi = std::abs(maxphi - minphi);
          if (etap >= etamin && etap < etamax && phip >= minphi && phip < maxphi) {
            // yield per cell with the outer ring halved, as the amplitude is in
            // estimateFlattenicity; the alternative normalises to the cell area
            rhoLatticeFV0AMC[isegment] +=
              flatSel.genFlatDetectorLikeNorm
                ? ((ieta == kOuterFV0RingIndex) ? 0.5f : 1.0f)
                : 1.0 / std::abs(dphi * kDetaFV0);
            multFV0++;
          }
          isegment++;
        }
      }
    }

    const float flattenicity =
      1.0 - getFlatenicity({rhoLatticeFV0AMC.data(), rhoLatticeFV0AMC.size()});
    if (fillQA) {
      rEventSelection.fill(HIST("hTrueFV0amplvsFlat"), multFV0, flattenicity);
    }
    return flattenicity;
  }
  // ====================== Flattenicity estimation ends =====================

  // Filters on V0s
  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between
  // daughters only
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0Sel.v0settingDCApostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0Sel.v0settingDCAnegtopv &&
                        aod::v0data::dcaV0daughters < v0Sel.v0settingDCAv0dau);

  Filter trackFilter =
    (nabs(aod::track::eta) < trkPid.cfgTrkEtaCut && aod::track::pt > trkPid.cfgTrkLowPtCut);

  using TrackCandidates = soa::Filtered<
    soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA,
              aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>>;

  void processDataRun3LambdaK0s(
    soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms,
              aod::PVMults>::iterator const& collision,
    soa::Filtered<aod::V0Datas> const& V0s, TrackCandidates const& tracks,
    soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
    aod::FV0As const& /*fv0s*/)
  {
    const bool own = (eventHistOwner == kOwnerDataV0);
    if (evSel.applyEvSel &&
        !(isEventSelected(collision, own))) { // Checking if the event passes the
                                              // selection criteria
      return;
    }

    auto vtxZ = collision.posZ();
    auto vtxY = collision.posY();
    auto vtxX = collision.posX();

    float flattenicity = estimateFlattenicity(collision, tracks);
    if (flattenicity < 0.f) {
      return;
    }
    if (own) {
      rEventSelection.fill(HIST("hEventsSelected"), nbinFlattenicity - 0.5);

      rEventSelection.fill(HIST("hVertexZ"), vtxZ);
      rEventSelection.fill(HIST("hFlattenicityDistribution"), flattenicity);
      rEventSelection.fill(HIST("hCentFT0M"), collision.centFT0M());
      rEventSelection.fill(HIST("hCentFT0MvsFlattenicity"), collision.centFT0M(), flattenicity);
      fillFlatDistRec(flattenicity, collision.centFT0M(), collision.isInelGt0());
      fillChargedRec<false>(tracks, flattenicity);
    }

    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<TrackCandidates>();
      const auto& negDaughterTrack = v0.negTrack_as<TrackCandidates>();

      if (!isSelectedDaughterTrack(posDaughterTrack, v0Sel.v0settingNTPCcrossedRows) ||
          !isSelectedDaughterTrack(negDaughterTrack, v0Sel.v0settingNTPCcrossedRows)) {
        continue;
      }
      float massK0s = v0.mK0Short();
      float massLambda = v0.mLambda();
      float massAntiLambda = v0.mAntiLambda();

      rKzeroShort.fill(HIST("hMassK0s"), massK0s);
      rLambda.fill(HIST("hMassLambda"), massLambda);
      rAntiLambda.fill(HIST("hMassAntiLambda"), massAntiLambda);

      float decayvtxX = v0.x();
      float decayvtxY = v0.y();
      float decayvtxZ = v0.z();

      float decaylength = std::sqrt(std::pow(decayvtxX - vtxX, 2) +
                                    std::pow(decayvtxY - vtxY, 2) +
                                    std::pow(decayvtxZ - vtxZ, 2));
      float v0p = std::sqrt(v0.pt() * v0.pt() + v0.pz() * v0.pz());

      float ctauK0s = decaylength * o2::constants::physics::MassK0Short / v0p;
      // same PDG mass for Lambda and AntiLambda
      float ctauLambda = decaylength * o2::constants::physics::MassLambda0 / v0p;

      float alpha = v0.alpha();
      float qtarm = v0.qtarm();

      // Cut on dynamic columns for K0s
      rCommonHist.fill(HIST("hArmPodoAlphavsQT"), alpha, qtarm);

      if (v0.v0cosPA() >= v0Sel.v0settingCosPAK0s &&
          v0.v0radius() >= v0Sel.v0settingRadiusK0s &&
          std::abs(posDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
          std::abs(negDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
          ctauK0s < v0Sel.v0settingcTauK0s &&
          std::abs(v0.rapidity(0)) <= v0Sel.v0settingRapidity &&
          std::abs(massLambda - o2::constants::physics::MassLambda0) > v0Sel.v0settingMassRejectionK0s &&
          std::abs(massAntiLambda - o2::constants::physics::MassLambda0) >
            v0Sel.v0settingMassRejectionK0s &&
          qtarm > v0Sel.v0settingArmePodoK0s * std::abs(alpha)) {

        rKzeroShort.fill(HIST("hMassK0sSelected"), massK0s);
        rKzeroShort.fill(HIST("hDCAV0DaughtersK0s"), v0.dcaV0daughters());
        rKzeroShort.fill(HIST("hV0CosPAK0s"), v0.v0cosPA());
        rKzeroShort.fill(HIST("hrapidityK0s"), v0.rapidity(0));
        rKzeroShort.fill(HIST("hctauK0s"), ctauK0s);
        rKzeroShort.fill(HIST("h2DdecayRadiusK0s"), v0.v0radius());
        rKzeroShort.fill(HIST("hMassK0spT"), massK0s, v0.pt());
        rKzeroShort.fill(HIST("hMassK0spTFlat"), massK0s, v0.pt(), flattenicity);
        rKzeroShort.fill(HIST("hArmPodoAlphavsQTK0sAfterCut"), alpha, qtarm);

        // Filling the PID of the V0 daughters in the region of the K0s peak
        if (std::abs(massK0s - o2::constants::physics::MassK0Short) < trkPid.pidQAWindowK0s) {
          rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"),
                           posDaughterTrack.tpcNSigmaPi(),
                           posDaughterTrack.tpcInnerParam());
          rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"),
                           negDaughterTrack.tpcNSigmaPi(),
                           negDaughterTrack.tpcInnerParam());
        }
      }

      // Cut on dynamic columns for Lambda
      if (v0.v0cosPA() >= v0Sel.v0settingCosPALambda &&
          v0.v0radius() >= v0Sel.v0settingRadiusLambda &&
          std::abs(posDaughterTrack.tpcNSigmaPr()) <= trkPid.nSigmaTPCProton &&
          std::abs(negDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
          ctauLambda < v0Sel.v0settingcTauLambda &&
          std::abs(v0.rapidity(1)) <= v0Sel.v0settingRapidity &&
          std::abs(massK0s - o2::constants::physics::MassK0Short) > v0Sel.v0settingMassRejectionLambda) {

        rLambda.fill(HIST("hMassLambdaSelected"), massLambda);
        rLambda.fill(HIST("hDCAV0DaughtersLambda"), v0.dcaV0daughters());
        rLambda.fill(HIST("hV0CosPALambda"), v0.v0cosPA());
        rLambda.fill(HIST("hrapidityLambda"), v0.rapidity(1));
        rLambda.fill(HIST("hctauLambda"), ctauLambda);
        rLambda.fill(HIST("h2DdecayRadiusLambda"), v0.v0radius());
        rLambda.fill(HIST("hMassLambdapT"), massLambda, v0.pt());
        rLambda.fill(HIST("hMassLambdapTFlat"), massLambda, v0.pt(), flattenicity);
        rLambda.fill(HIST("hDcaV0ToPVLambda"), v0.dcav0topv(), v0.pt(), flattenicity);

        // Filling the PID of the V0 daughters in the region of the Lambda peak
        if (std::abs(massLambda - o2::constants::physics::MassLambda0) < trkPid.pidQAWindowLambda) {
          rLambda.fill(HIST("hNSigmaPosProtonFromLambda"),
                       posDaughterTrack.tpcNSigmaPr(),
                       posDaughterTrack.tpcInnerParam());
          rLambda.fill(HIST("hNSigmaNegPionFromLambda"),
                       negDaughterTrack.tpcNSigmaPi(),
                       negDaughterTrack.tpcInnerParam());
        }
      }

      // Cut on dynamic columns for AntiLambda
      if (v0.v0cosPA() >= v0Sel.v0settingCosPALambda &&
          v0.v0radius() >= v0Sel.v0settingRadiusLambda &&
          std::abs(posDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
          std::abs(negDaughterTrack.tpcNSigmaPr()) <= trkPid.nSigmaTPCProton &&
          ctauLambda < v0Sel.v0settingcTauLambda &&
          std::abs(v0.rapidity(2)) <= v0Sel.v0settingRapidity &&
          std::abs(massK0s - o2::constants::physics::MassK0Short) > v0Sel.v0settingMassRejectionLambda) {

        rAntiLambda.fill(HIST("hMassAntiLambdaSelected"), massAntiLambda);
        rAntiLambda.fill(HIST("hDCAV0DaughtersAntiLambda"),
                         v0.dcaV0daughters());
        rAntiLambda.fill(HIST("hV0CosPAAntiLambda"), v0.v0cosPA());
        rAntiLambda.fill(HIST("hrapidityAntiLambda"), v0.rapidity(2));
        rAntiLambda.fill(HIST("hctauAntiLambda"), ctauLambda);
        rAntiLambda.fill(HIST("h2DdecayRadiusAntiLambda"), v0.v0radius());
        rAntiLambda.fill(HIST("hMassAntiLambdapT"), massAntiLambda, v0.pt());

        rAntiLambda.fill(HIST("hMassAntiLambdapTFlat"), massAntiLambda, v0.pt(), flattenicity);
        rAntiLambda.fill(HIST("hDcaV0ToPVAntiLambda"), v0.dcav0topv(), v0.pt(), flattenicity);
        // Filling the PID of the V0 daughters in the region of the AntiLambda
        // peak
        if (std::abs(massAntiLambda - o2::constants::physics::MassLambda0) < trkPid.pidQAWindowLambda) {
          rAntiLambda.fill(HIST("hNSigmaPosPionFromAntiLambda"),
                           posDaughterTrack.tpcNSigmaPi(),
                           posDaughterTrack.tpcInnerParam());
          rAntiLambda.fill(HIST("hNSigmaNegProtonFromAntiLambda"),
                           negDaughterTrack.tpcNSigmaPr(),
                           negDaughterTrack.tpcInnerParam());
        }
      }
    }
  }

  using TrackCandidatesMC =
    soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA,
                            aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr,
                            aod::McTrackLabels>>;

  Preslice<soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>>> perCol = aod::v0data::collisionId;
  Preslice<TrackCandidatesMC> perColTracksMC = aod::track::collisionId;
  // needed by the sliceByCached calls below, else they throw "Disabled cache"
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache1;

  void processRecMCLambdaK0s(
    soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms,
              aod::PVMults, aod::McCollisionLabels> const& collisions,
    soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>> const& V0s, aod::McCollisions const&, TrackCandidatesMC const& tracks,
    soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
    aod::FV0As const& /*fv0s*/, aod::McParticles const& mcParticles)
  {
    const bool own = (eventHistOwner == kOwnerRecMCV0);
    for (const auto& collision : collisions) {
      if (evSel.applyEvSel &&
          !(isEventSelected(collision, own))) { // Checking if the event passes the
                                                // selection criteria
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto vtxZ = collision.posZ();
      auto vtxY = collision.posY();
      auto vtxX = collision.posX();

      auto tracksThisCollision = tracks.sliceBy(perColTracksMC, collision.globalIndex());
      float flattenicity = estimateFlattenicity(collision, tracksThisCollision);
      if (flattenicity < 0.f) {
        continue;
      }
      if (own) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinFlattenicity - 0.5);

        rEventSelection.fill(HIST("hVertexZ"), vtxZ);
        rEventSelection.fill(HIST("hFlattenicityDistribution"), flattenicity);
        rEventSelection.fill(HIST("hCentFT0M"), collision.centFT0M());
        rEventSelection.fill(HIST("hCentFT0MvsFlattenicity"), collision.centFT0M(), flattenicity);
        fillFlatDistRec(flattenicity, collision.centFT0M(), collision.isInelGt0());
        fillChargedRec<true>(tracksThisCollision, flattenicity);
      }

      auto v0sThisCollision = V0s.sliceBy(perCol, collision.globalIndex());
      const auto& mcCollision = collision.mcCollision_as<aod::McCollisions>();

      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache1);
      const float flattenicityMCGen = estimateFlattenicityFV0MC(particlesInCollision, own);
      if (own) {
        rEventSelection.fill(HIST("hFlattenicityDistributionMCGen_Rec"), flattenicityMCGen);
        rEventSelection.fill(HIST("hFlattenicity_Corr_Gen_vs_Rec"), flattenicityMCGen, flattenicity);
        rEventSelection.fill(HIST("hFlatGenVsRecFine"), flattenicityMCGen, flattenicity);
        fillFlatDistGenInRec(flattenicityMCGen, collision.centFT0M(), collision.isInelGt0());
      }

      for (const auto& v0 : v0sThisCollision) {

        const auto& posDaughterTrack = v0.posTrack_as<TrackCandidatesMC>();
        const auto& negDaughterTrack = v0.negTrack_as<TrackCandidatesMC>();

        if (!isSelectedDaughterTrack(posDaughterTrack, v0Sel.v0settingNTPCcrossedRows) ||
            !isSelectedDaughterTrack(negDaughterTrack, v0Sel.v0settingNTPCcrossedRows)) {
          continue;
        }

        if (!v0.has_mcParticle()) {
          continue;
        }

        float massK0s = v0.mK0Short();
        float massLambda = v0.mLambda();
        float massAntiLambda = v0.mAntiLambda();

        rKzeroShort.fill(HIST("hMassK0s"), massK0s);
        rLambda.fill(HIST("hMassLambda"), massLambda);
        rAntiLambda.fill(HIST("hMassAntiLambda"), massAntiLambda);

        float decayvtxX = v0.x();
        float decayvtxY = v0.y();
        float decayvtxZ = v0.z();

        float decaylength = std::sqrt(std::pow(decayvtxX - vtxX, 2) +
                                      std::pow(decayvtxY - vtxY, 2) +
                                      std::pow(decayvtxZ - vtxZ, 2));
        float v0p = std::sqrt(v0.pt() * v0.pt() + v0.pz() * v0.pz());

        float ctauK0s = decaylength * o2::constants::physics::MassK0Short / v0p;
        // same PDG mass for Lambda and AntiLambda
        float ctauLambda = decaylength * o2::constants::physics::MassLambda0 / v0p;

        float alpha = v0.alpha();
        float qtarm = v0.qtarm();
        rCommonHist.fill(HIST("hArmPodoAlphavsQT"), alpha, qtarm);

        auto v0mcParticle = v0.mcParticle();
        const bool isPrimaryV0 = v0mcParticle.isPhysicalPrimary();
        const bool keepForEfficiency = isPrimaryV0 || !eventClass.requirePrimaryMC;
        float motherPt = -1.f;
        // Cut on dynamic columns for K0s

        if (v0mcParticle.pdgCode() == PDG_t::kK0Short && keepForEfficiency &&
            v0.v0cosPA() >= v0Sel.v0settingCosPAK0s &&
            v0.v0radius() >= v0Sel.v0settingRadiusK0s &&
            std::abs(posDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
            std::abs(negDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
            ctauK0s < v0Sel.v0settingcTauK0s &&
            std::abs(v0.rapidity(0)) <= v0Sel.v0settingRapidity &&
            std::abs(massLambda - o2::constants::physics::MassLambda0) > v0Sel.v0settingMassRejectionK0s &&
            std::abs(massAntiLambda - o2::constants::physics::MassLambda0) >
              v0Sel.v0settingMassRejectionK0s &&
            qtarm > v0Sel.v0settingArmePodoK0s * std::abs(alpha)) {

          rKzeroShort.fill(HIST("hMassK0sSelected"), massK0s);
          rKzeroShort.fill(HIST("hDCAV0DaughtersK0s"), v0.dcaV0daughters());
          rKzeroShort.fill(HIST("hV0CosPAK0s"), v0.v0cosPA());
          rKzeroShort.fill(HIST("hrapidityK0s"), v0.rapidity(0));
          rKzeroShort.fill(HIST("hctauK0s"), ctauK0s);
          rKzeroShort.fill(HIST("h2DdecayRadiusK0s"), v0.v0radius());
          rKzeroShort.fill(HIST("hMassK0spT"), massK0s, v0.pt());
          rKzeroShort.fill(HIST("hMassK0spTFlat"), massK0s, v0.pt(), flattenicity);
          rKzeroShort.fill(HIST("hMassK0spTTrueFlat"), massK0s, v0.pt(), flattenicityMCGen);
          rKzeroShort.fill(HIST("hPtResK0s"), v0mcParticle.pt(),
                           (v0.pt() - v0mcParticle.pt()) / v0mcParticle.pt());
          rKzeroShort.fill(HIST("hArmPodoAlphavsQTK0sAfterCut"), alpha, qtarm);

          // Filling the PID of the V0 daughters in the region of the K0s peak
          if (std::abs(massK0s - o2::constants::physics::MassK0Short) < trkPid.pidQAWindowK0s) {
            rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"),
                             posDaughterTrack.tpcNSigmaPi(),
                             posDaughterTrack.tpcInnerParam());
            rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"),
                             negDaughterTrack.tpcNSigmaPi(),
                             negDaughterTrack.tpcInnerParam());
          }
        }

        // Cut on dynamic columns for Lambda
        if (v0mcParticle.pdgCode() == PDG_t::kLambda0 &&
            v0.v0cosPA() >= v0Sel.v0settingCosPALambda &&
            v0.v0radius() >= v0Sel.v0settingRadiusLambda &&
            std::abs(posDaughterTrack.tpcNSigmaPr()) <= trkPid.nSigmaTPCProton &&
            std::abs(negDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
            ctauLambda < v0Sel.v0settingcTauLambda &&
            std::abs(v0.rapidity(1)) <= v0Sel.v0settingRapidity &&
            std::abs(massK0s - o2::constants::physics::MassK0Short) > v0Sel.v0settingMassRejectionLambda) {

          // secondaries kept out of the efficiency numerator, booked for the feed-down
          if (!isPrimaryV0) {
            const int motherIndex = getFeedDownMother(v0mcParticle, motherPt);
            rLambda.fill(HIST("hMassFeedDownLambdapTFlat"), massLambda, v0.pt(), flattenicity);
            rLambda.fill(HIST("hFeedDownLambdaPtVsMotherPt"), v0.pt(), motherPt);
            rLambda.fill(HIST("hFeedDownLambdaMatrix"), v0.pt(), motherPt, flattenicity, motherIndex);
            rLambda.fill(HIST("hFeedDownLambdaMotherPdg"), motherIndex);
            rLambda.fill(HIST("hDcaV0ToPVLambdaSec"), v0.dcav0topv(), v0.pt(), flattenicity);
          } else {
            rLambda.fill(HIST("hDcaV0ToPVLambdaPrim"), v0.dcav0topv(), v0.pt(), flattenicity);
          }
          if (keepForEfficiency) {

            rLambda.fill(HIST("hMassLambdaSelected"), massLambda);
            rLambda.fill(HIST("hDCAV0DaughtersLambda"), v0.dcaV0daughters());
            rLambda.fill(HIST("hV0CosPALambda"), v0.v0cosPA());
            rLambda.fill(HIST("hrapidityLambda"), v0.rapidity(1));
            rLambda.fill(HIST("hctauLambda"), ctauLambda);
            rLambda.fill(HIST("h2DdecayRadiusLambda"), v0.v0radius());
            rLambda.fill(HIST("hMassLambdapT"), massLambda, v0.pt());
            rLambda.fill(HIST("hMassLambdapTFlat"), massLambda, v0.pt(), flattenicity);
            rLambda.fill(HIST("hMassLambdapTTrueFlat"), massLambda, v0.pt(), flattenicityMCGen);
            rLambda.fill(HIST("hDcaV0ToPVLambda"), v0.dcav0topv(), v0.pt(), flattenicity);
            rLambda.fill(HIST("hPtResLambda"), v0mcParticle.pt(),
                         (v0.pt() - v0mcParticle.pt()) / v0mcParticle.pt());

            // Filling the PID of the V0 daughters in the region of the Lambda peak
            if (std::abs(massLambda - o2::constants::physics::MassLambda0) < trkPid.pidQAWindowLambda) {
              rLambda.fill(HIST("hNSigmaPosProtonFromLambda"),
                           posDaughterTrack.tpcNSigmaPr(),
                           posDaughterTrack.tpcInnerParam());
              rLambda.fill(HIST("hNSigmaNegPionFromLambda"),
                           negDaughterTrack.tpcNSigmaPi(),
                           negDaughterTrack.tpcInnerParam());
            }
          }
        }

        // Cut on dynamic columns for AntiLambda
        if (v0mcParticle.pdgCode() == PDG_t::kLambda0Bar &&
            v0.v0cosPA() >= v0Sel.v0settingCosPALambda &&
            v0.v0radius() >= v0Sel.v0settingRadiusLambda &&
            std::abs(posDaughterTrack.tpcNSigmaPi()) <= trkPid.nSigmaTPCPion &&
            std::abs(negDaughterTrack.tpcNSigmaPr()) <= trkPid.nSigmaTPCProton &&
            ctauLambda < v0Sel.v0settingcTauLambda &&
            std::abs(v0.rapidity(2)) <= v0Sel.v0settingRapidity &&
            std::abs(massK0s - o2::constants::physics::MassK0Short) > v0Sel.v0settingMassRejectionLambda) {

          // secondaries kept out of the efficiency numerator, booked for the feed-down
          if (!isPrimaryV0) {
            const int motherIndex = getFeedDownMother(v0mcParticle, motherPt);
            rAntiLambda.fill(HIST("hMassFeedDownAntiLambdapTFlat"), massAntiLambda, v0.pt(), flattenicity);
            rAntiLambda.fill(HIST("hFeedDownAntiLambdaPtVsMotherPt"), v0.pt(), motherPt);
            rAntiLambda.fill(HIST("hFeedDownAntiLambdaMatrix"), v0.pt(), motherPt, flattenicity, motherIndex);
            rAntiLambda.fill(HIST("hFeedDownAntiLambdaMotherPdg"), motherIndex);
            rAntiLambda.fill(HIST("hDcaV0ToPVAntiLambdaSec"), v0.dcav0topv(), v0.pt(), flattenicity);
          } else {
            rAntiLambda.fill(HIST("hDcaV0ToPVAntiLambdaPrim"), v0.dcav0topv(), v0.pt(), flattenicity);
          }
          if (keepForEfficiency) {

            rAntiLambda.fill(HIST("hMassAntiLambdaSelected"), massAntiLambda);
            rAntiLambda.fill(HIST("hDCAV0DaughtersAntiLambda"),
                             v0.dcaV0daughters());
            rAntiLambda.fill(HIST("hV0CosPAAntiLambda"), v0.v0cosPA());
            rAntiLambda.fill(HIST("hrapidityAntiLambda"), v0.rapidity(2));
            rAntiLambda.fill(HIST("hctauAntiLambda"), ctauLambda);
            rAntiLambda.fill(HIST("h2DdecayRadiusAntiLambda"), v0.v0radius());
            rAntiLambda.fill(HIST("hMassAntiLambdapT"), massAntiLambda, v0.pt());
            rAntiLambda.fill(HIST("hMassAntiLambdapTFlat"), massAntiLambda, v0.pt(), flattenicity);
            rAntiLambda.fill(HIST("hMassAntiLambdapTTrueFlat"), massAntiLambda, v0.pt(), flattenicityMCGen);
            rAntiLambda.fill(HIST("hDcaV0ToPVAntiLambda"), v0.dcav0topv(), v0.pt(), flattenicity);
            rAntiLambda.fill(HIST("hPtResAntiLambda"), v0mcParticle.pt(),
                             (v0.pt() - v0mcParticle.pt()) / v0mcParticle.pt());

            // Filling the PID of the V0 daughters in the region of the AntiLambda
            // peak
            if (std::abs(massAntiLambda - o2::constants::physics::MassLambda0) < trkPid.pidQAWindowLambda) {
              rAntiLambda.fill(HIST("hNSigmaPosPionFromAntiLambda"),
                               posDaughterTrack.tpcNSigmaPi(),
                               posDaughterTrack.tpcInnerParam());
              rAntiLambda.fill(HIST("hNSigmaNegProtonFromAntiLambda"),
                               negDaughterTrack.tpcNSigmaPr(),
                               negDaughterTrack.tpcInnerParam());
            }
          }
        }
      }

      const bool isInelGt0Rec = collision.isInelGt0();

      int nchGenInRec = 0;
      for (const auto& mcParticle : particlesInCollision) {
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }

        if (eventClass.fillChargedQA && std::abs(mcParticle.eta()) <= eventClass.cfgEtaChargedCut) {
          auto pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
          if (pdgParticle && std::abs(pdgParticle->Charge()) > kMinCharge) {
            nchGenInRec++;
            rCharged.fill(HIST("hPtChVsFlatGenInRec"), mcParticle.pt(), flattenicity);
          }
        }

        // normalisation of hFeedDownLambdaMatrix, before the Lambda rapidity cut
        if (std::abs(mcParticle.y()) <= eventClass.cfgFeedDownMotherRapidity) {
          const int motherSpecies = feedDownSpecies(mcParticle.pdgCode());
          if (motherSpecies >= 0) {
            rLambda.fill(HIST("hGenFeedDownMotherPtFlat"), mcParticle.pt(), flattenicity, motherSpecies);
          }
        }

        if (std::abs(mcParticle.y()) > v0Sel.v0settingRapidity) {
          continue;
        }

        if (mcParticle.pdgCode() == PDG_t::kK0Short) {
          rKzeroShort.fill(HIST("Generated_MCRecoCollCheck_INEL_K0Short"), mcParticle.pt(), flattenicity);
          if (isInelGt0Rec) {
            rKzeroShort.fill(HIST("Generated_MCRecoCollCheck_INELgt0_K0Short"), mcParticle.pt(), flattenicity);
            rKzeroShort.fill(HIST("Generated_MCRecoCollCheck_INELgt0_K0Short_TrueFlat"), mcParticle.pt(), flattenicityMCGen);
          }
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0) {
          rLambda.fill(HIST("Generated_MCRecoCollCheck_INEL_Lambda"), mcParticle.pt(), flattenicity);
          if (isInelGt0Rec) {
            rLambda.fill(HIST("Generated_MCRecoCollCheck_INELgt0_Lambda"), mcParticle.pt(), flattenicity);
            rLambda.fill(HIST("Generated_MCRecoCollCheck_INELgt0_Lambda_TrueFlat"), mcParticle.pt(), flattenicityMCGen);
          }
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0Bar) {
          rAntiLambda.fill(HIST("Generated_MCRecoCollCheck_INEL_AntiLambda"), mcParticle.pt(), flattenicity);
          if (isInelGt0Rec) {
            rAntiLambda.fill(HIST("Generated_MCRecoCollCheck_INELgt0_AntiLambda"), mcParticle.pt(), flattenicity);
            rAntiLambda.fill(HIST("Generated_MCRecoCollCheck_INELgt0_AntiLambda_TrueFlat"), mcParticle.pt(), flattenicityMCGen);
          }
        }
      }
      if (eventClass.fillChargedQA) {
        rCharged.fill(HIST("hNchVsFlatGenInRec"), nchGenInRec, flattenicity);
      }
    }
  }

  // Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < evSel.cutzvertex);
  void processGenMC(
    o2::aod::McCollision const& mcCollision, const soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels, o2::aod::CentFT0Ms, aod::PVMults>>& collisions, TrackCandidatesMC const& tracks,
    soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
    aod::FV0As const& /*fv0s*/, o2::aod::McParticles const& mcParticles)
  {

    // Both estimates are always computed: the reconstructed one classifies the
    // numerator of the closure, the true one its denominator. Without a
    // reconstructed counterpart the sentinel is kept, it falls in the underflow
    // so the loss counters stay complete.
    float flattenicityRec = kInvalidFlattenicity;
    float centFT0M = -1.f;
    for (const auto& collision : collisions) {
      if (evSel.applyEvSel && !isEventSelected(collision, false)) {
        continue;
      }
      auto tracksThisCollision = tracks.sliceBy(perColTracksMC, collision.globalIndex());
      flattenicityRec = estimateFlattenicity(collision, tracksThisCollision, fillFlattenicityQAInGenMC);
      if (flattenicityRec >= 0.f) {
        centFT0M = collision.centFT0M();
        break;
      }
    }
    const float flattenicityTrue = estimateFlattenicityFV0MC(mcParticles, fillFlattenicityQAInGenMC);
    rEventSelection.fill(HIST("hFlattenicityDistributionRecMCGen"), flattenicityRec);
    rEventSelection.fill(HIST("hFlattenicityDistributionMCGen"), flattenicityTrue);

    // which of the two drives the loss corrections written without a suffix
    const float flattenicity = flatSel.flattenicityforLossCorrRec ? flattenicityRec : flattenicityTrue;

    //====================================
    //===== Event Loss Denominator =======
    //====================================

    rEventSelection.fill(HIST("hNEventsMCGen"), 0.5);

    if (std::abs(mcCollision.posZ()) > evSel.cutzvertex) {
      return;
    }
    rEventSelection.fill(HIST("hNEventsMCGen"), 1.5);

    // The FT0M class has no generator-level counterpart, so with a window
    // requested the collision only belongs to it through an accepted
    // reconstructed one, which is what leaves centFT0M non-negative above.
    // Otherwise the generated denominators stay MB while the numerator is HM.
    if (eventClass.applyCentSel && centFT0M < 0.f) {
      return;
    }

    rEventSelection.fill(HIST("hFlat_GenColl_MC"), flattenicity);

    const bool isINELgt0true = pwglf::isINELgtNmc(mcParticles, 0, pdg);

    if (isINELgt0true) {
      rEventSelection.fill(HIST("hNEventsMCGen"), 2.5);
      rEventSelection.fill(HIST("hFlat_GenColl_MC_INELgt0"), flattenicity);
      rEventSelection.fill(HIST("hFlat_GenColl_MC_INELgt0_TrueFlat"), flattenicityTrue);
    }

    fillFlatDistGen(flattenicityTrue, centFT0M, isINELgt0true);
    fillChargedGen(mcParticles, flattenicity, false);
    fillChargedGen(mcParticles, flattenicityTrue, true);

    //=====================================
    //===== Signal Loss Denominator =======
    //=====================================

    for (const auto& mcParticle : mcParticles) {

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const bool inV0Rapidity = std::abs(mcParticle.y()) <= v0Sel.v0settingRapidity;
      const bool inCascRapidity = std::abs(mcParticle.y()) <= cascSel.cascsettingRapidity;

      if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kK0Short) {
        rKzeroShort.fill(HIST("pGen_MCGenColl_INEL_K0Short"), mcParticle.pt(), flattenicity);
        if (isINELgt0true) {
          rKzeroShort.fill(HIST("pGen_MCGenColl_INELgt0_K0Short"), mcParticle.pt(), flattenicity);
          rKzeroShort.fill(HIST("pGen_MCGenColl_INELgt0_K0Short_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
      if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kLambda0) {
        rLambda.fill(HIST("pGen_MCGenColl_INEL_Lambda"), mcParticle.pt(), flattenicity);
        if (isINELgt0true) {
          rLambda.fill(HIST("pGen_MCGenColl_INELgt0_Lambda"), mcParticle.pt(), flattenicity);
          rLambda.fill(HIST("pGen_MCGenColl_INELgt0_Lambda_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
      if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kLambda0Bar) {
        rAntiLambda.fill(HIST("pGen_MCGenColl_INEL_AntiLambda"), mcParticle.pt(), flattenicity);
        if (isINELgt0true) {
          rAntiLambda.fill(HIST("pGen_MCGenColl_INELgt0_AntiLambda"), mcParticle.pt(), flattenicity);
          rAntiLambda.fill(HIST("pGen_MCGenColl_INELgt0_AntiLambda_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
      if (inCascRapidity && std::abs(mcParticle.pdgCode()) == PDG_t::kXiMinus) {
        rXi.fill(HIST("pGen_MCGenColl_INEL_Xi"), mcParticle.pt(), flattenicity);
        if (isINELgt0true) {
          rXi.fill(HIST("pGen_MCGenColl_INELgt0_Xi"), mcParticle.pt(), flattenicity);
          rXi.fill(HIST("pGen_MCGenColl_INELgt0_Xi_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
    }

    int recoCollIndexINEL = 0;
    int recoCollIndexINELgt0 = 0;
    for (const auto& collision : collisions) { // loop on reconstructed collisions

      //=====================================
      //====== Event Split Numerator ========
      //=====================================

      rEventSelection.fill(HIST("hNEventsMCReco"), 0.5);
      // the per-cut ladder of hEventsSelected belongs to the reconstructed-level
      // process functions; filling it here as well would double count every
      // collision whenever processGenMC runs alongside one of them, which is the
      // normal MC configuration. hNEventsMCReco already records all/passed here.
      if (evSel.applyEvSel && !isEventSelected(collision, false)) {
        continue;
      }
      rEventSelection.fill(HIST("hEventsSelected"), nbin - 0.5);
      rEventSelection.fill(HIST("hNEventsMCReco"), 1.5);
      rEventSelection.fill(HIST("hFlat_RecoColl_MC"), flattenicity);

      recoCollIndexINEL++;

      if (collision.isInelGt0() && isINELgt0true) {
        rEventSelection.fill(HIST("hNEventsMCReco"), 2.5);
        rEventSelection.fill(HIST("hFlat_RecoColl_MC_INELgt0"), flattenicity);
        rEventSelection.fill(HIST("hFlat_RecoColl_MC_INELgt0_TrueFlat"), flattenicityTrue);

        recoCollIndexINELgt0++;
      }

      //=====================================
      //======== Sgn Split Numerator ========
      //=====================================

      for (const auto& mcParticle : mcParticles) {

        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        const bool inV0Rapidity = std::abs(mcParticle.y()) <= v0Sel.v0settingRapidity;
        const bool inCascRapidity = std::abs(mcParticle.y()) <= cascSel.cascsettingRapidity;

        if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kK0Short) {
          rKzeroShort.fill(HIST("Generated_MCRecoColl_INEL_K0Short"), mcParticle.pt(), flattenicity);
          if (recoCollIndexINELgt0 > 0) {
            rKzeroShort.fill(HIST("Generated_MCRecoColl_INELgt0_K0Short"), mcParticle.pt(), flattenicity);
            rKzeroShort.fill(HIST("Generated_MCRecoColl_INELgt0_K0Short_TrueFlat"), mcParticle.pt(), flattenicityTrue);
          }
        }
        if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kLambda0) {
          rLambda.fill(HIST("Generated_MCRecoColl_INEL_Lambda"), mcParticle.pt(), flattenicity);
          if (recoCollIndexINELgt0 > 0) {
            rLambda.fill(HIST("Generated_MCRecoColl_INELgt0_Lambda"), mcParticle.pt(), flattenicity);
            rLambda.fill(HIST("Generated_MCRecoColl_INELgt0_Lambda_TrueFlat"), mcParticle.pt(), flattenicityTrue);
          }
        }
        if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kLambda0Bar) {
          rAntiLambda.fill(HIST("Generated_MCRecoColl_INEL_AntiLambda"), mcParticle.pt(), flattenicity);
          if (recoCollIndexINELgt0 > 0) {
            rAntiLambda.fill(HIST("Generated_MCRecoColl_INELgt0_AntiLambda"), mcParticle.pt(), flattenicity);
            rAntiLambda.fill(HIST("Generated_MCRecoColl_INELgt0_AntiLambda_TrueFlat"), mcParticle.pt(), flattenicityTrue);
          }
        }
        if (inCascRapidity && std::abs(mcParticle.pdgCode()) == PDG_t::kXiMinus) {
          rXi.fill(HIST("Generated_MCRecoColl_INEL_Xi"), mcParticle.pt(), flattenicity);
          if (recoCollIndexINELgt0 > 0) {
            rXi.fill(HIST("Generated_MCRecoColl_INELgt0_Xi"), mcParticle.pt(), flattenicity);
            rXi.fill(HIST("Generated_MCRecoColl_INELgt0_Xi_TrueFlat"), mcParticle.pt(), flattenicityTrue);
          }
        }
      }
    }

    // From now on keep only mc collisions with at least one reconstructed collision (INEL)
    if (recoCollIndexINEL < 1) {
      return;
    }

    //=====================================
    //====== Event Loss Numerator =========
    //=====================================

    rEventSelection.fill(HIST("hNEventsMCGenReco"), 0.5);
    rEventSelection.fill(HIST("hFlat_GenRecoColl_MC"), flattenicity);

    if (recoCollIndexINELgt0 > 0) {
      rEventSelection.fill(HIST("hNEventsMCGenReco"), 1.5);
      rEventSelection.fill(HIST("hFlat_GenRecoColl_MC_INELgt0"), flattenicity);
      rEventSelection.fill(HIST("hFlat_GenRecoColl_MC_INELgt0_TrueFlat"), flattenicityTrue);
    }

    //=====================================
    //===== Signal Loss Numerator =========
    //=====================================

    for (const auto& mcParticle : mcParticles) {

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const bool inV0Rapidity = std::abs(mcParticle.y()) <= v0Sel.v0settingRapidity;
      const bool inCascRapidity = std::abs(mcParticle.y()) <= cascSel.cascsettingRapidity;

      if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kK0Short) {
        rKzeroShort.fill(HIST("pGen_MCGenRecoColl_INEL_K0Short"), mcParticle.pt(), flattenicity);
        if (recoCollIndexINELgt0 > 0) {
          rKzeroShort.fill(HIST("pGen_MCGenRecoColl_INELgt0_K0Short"), mcParticle.pt(), flattenicity);
          rKzeroShort.fill(HIST("pGen_MCGenRecoColl_INELgt0_K0Short_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
      if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kLambda0) {
        rLambda.fill(HIST("pGen_MCGenRecoColl_INEL_Lambda"), mcParticle.pt(), flattenicity);
        if (recoCollIndexINELgt0 > 0) {
          rLambda.fill(HIST("pGen_MCGenRecoColl_INELgt0_Lambda"), mcParticle.pt(), flattenicity);
          rLambda.fill(HIST("pGen_MCGenRecoColl_INELgt0_Lambda_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
      if (inV0Rapidity && mcParticle.pdgCode() == PDG_t::kLambda0Bar) {
        rAntiLambda.fill(HIST("pGen_MCGenRecoColl_INEL_AntiLambda"), mcParticle.pt(), flattenicity);
        if (recoCollIndexINELgt0 > 0) {
          rAntiLambda.fill(HIST("pGen_MCGenRecoColl_INELgt0_AntiLambda"), mcParticle.pt(), flattenicity);
          rAntiLambda.fill(HIST("pGen_MCGenRecoColl_INELgt0_AntiLambda_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
      if (inCascRapidity && std::abs(mcParticle.pdgCode()) == PDG_t::kXiMinus) {
        rXi.fill(HIST("pGen_MCGenRecoColl_INEL_Xi"), mcParticle.pt(), flattenicity);
        if (recoCollIndexINELgt0 > 0) {
          rXi.fill(HIST("pGen_MCGenRecoColl_INELgt0_Xi"), mcParticle.pt(), flattenicity);
          rXi.fill(HIST("pGen_MCGenRecoColl_INELgt0_Xi_TrueFlat"), mcParticle.pt(), flattenicityTrue);
        }
      }
    }
  }
  // Cascade Analysis Starts here
  using DauTracks = soa::Join<aod::TracksIU, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using LabeledDauTracks = soa::Join<DauTracks, aod::McTrackLabels>;
  using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

  static constexpr float kCtauXi = 4.91; // Xi lifetime, in cm

  // for Xi- the baryon daughter is the positive one, the other way around for Xi+
  template <typename TCollision, typename TCascade, typename TTrack>
  bool isSelectedXi(TCollision const& collision, TCascade const& casc,
                    TTrack const& posDaughterTrack, TTrack const& negDaughterTrack,
                    TTrack const& bacDaughterTrack, float ctauXi)
  {
    const bool isXiMinus = casc.sign() < 0;
    const auto& protonDaughter = isXiMinus ? posDaughterTrack : negDaughterTrack;
    const auto& pionDaughter = isXiMinus ? negDaughterTrack : posDaughterTrack;

    // PID, one hypothesis per daughter
    if (std::abs(protonDaughter.tpcNSigmaPr()) > trkPid.nSigmaTPCProton ||
        std::abs(pionDaughter.tpcNSigmaPi()) > trkPid.nSigmaTPCPion ||
        std::abs(bacDaughterTrack.tpcNSigmaPi()) > trkPid.nSigmaTPCPion) {
      return false;
    }

    // track quality
    if (!isSelectedDaughterTrack(posDaughterTrack, cascSel.nTPCcrossedRows) ||
        !isSelectedDaughterTrack(negDaughterTrack, cascSel.nTPCcrossedRows) ||
        !isSelectedDaughterTrack(bacDaughterTrack, cascSel.nTPCcrossedRows)) {
      return false;
    }

    // topology
    const float dcav0pv = casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ());
    const float cosPAcasc = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
    const float cosPAv0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());

    // DCA to PV per daughter role, not per charge
    const float dcaBaryonToPV = isXiMinus ? casc.dcapostopv() : casc.dcanegtopv();
    const float dcaMesonToPV = isXiMinus ? casc.dcanegtopv() : casc.dcapostopv();

    if (std::abs(dcaBaryonToPV) < cascSel.cascsettingDCAbaryontopv || std::abs(dcaMesonToPV) < cascSel.cascsettingDCAmesontopv ||
        std::abs(casc.dcabachtopv()) < v0Sel.v0settingDCAbactopv || casc.dcaV0daughters() > v0Sel.v0settingDCAv0dau ||
        dcav0pv < cascSel.cascsettingDCAv0toPV || casc.dcacascdaughters() > cascSel.cascsettingDCAv0bach) {
      return false;
    }
    // bachelor-baryon veto
    if (casc.bachBaryonCosPA() > cascSel.cascsettingCosPAbaryonbach ||
        std::abs(casc.bachBaryonDCAxyToPV()) < cascSel.cascsettingDCAxybaryonbach) {
      return false;
    }
    if (cosPAcasc < cascSel.cascsettingCosPAcascPV || cosPAv0 < cascSel.cascsettingCosPAv0PV ||
        casc.cascradius() < cascSel.cascsettingcascradius || casc.v0radius() < cascSel.cascsettingv0radius) {
      return false;
    }

    // kinematics and competing species
    if (std::abs(casc.yXi()) > cascSel.cascsettingRapidity || ctauXi > kCtauXi * cascSel.cascsettingproplifetime) {
      return false;
    }
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascSel.cascsettingMassRejectionLambdaXi ||
        std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < cascSel.cascsettingMassRejectioOmegaXi) {
      return false;
    }

    return true;
  }

  void processDataRun3Cascade(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms,
                                        aod::PVMults>::iterator const& collision,
                              aod::CascDataExt const& Cascades,
                              aod::V0Datas const&, DauTracks const& tracks,
                              soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
                              aod::FV0As const& /*fv0s*/)
  {
    const bool own = (eventHistOwner == kOwnerDataCasc);
    if (evSel.applyEvSel &&
        !(isEventSelected(collision, own))) { // Checking if the event passes the
                                              // selection criteria
      return;
    }

    auto vtxZ = collision.posZ();
    auto vtxY = collision.posY();
    auto vtxX = collision.posX();

    float flattenicity = estimateFlattenicity(collision, tracks);
    if (flattenicity < 0.f) {
      return;
    }
    if (own) {
      rEventSelection.fill(HIST("hEventsSelected"), nbinFlattenicity - 0.5);

      rEventSelection.fill(HIST("hVertexZ"), vtxZ);
      rEventSelection.fill(HIST("hFlattenicityDistribution"), flattenicity);
      rEventSelection.fill(HIST("hCentFT0M"), collision.centFT0M());
      rEventSelection.fill(HIST("hCentFT0MvsFlattenicity"), collision.centFT0M(), flattenicity);
      fillFlatDistRec(flattenicity, collision.centFT0M(), collision.isInelGt0());
      fillChargedRec<false>(tracks, flattenicity);
    }

    for (const auto& casc : Cascades) {

      auto posDaughterTrack = casc.posTrack_as<DauTracks>();
      auto negDaughterTrack = casc.negTrack_as<DauTracks>();
      auto bacDaughterTrack = casc.bachelor_as<DauTracks>();

      float cascPos = std::hypot(casc.x() - vtxX, casc.y() - vtxY, casc.z() - vtxZ);
      float cascTotMom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctauXi = o2::constants::physics::MassXiMinus * cascPos / (cascTotMom + 1e-13);
      float cosPAv0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
      float massXi = casc.mXi();
      rXi.fill(HIST("hMassXi"), massXi);
      // Cascade
      if (isSelectedXi(collision, casc, posDaughterTrack, negDaughterTrack, bacDaughterTrack, ctauXi)) {

        const bool isXiMinus = casc.sign() < 0;
        const auto& protonDaughter = isXiMinus ? posDaughterTrack : negDaughterTrack;
        const auto& pionDaughter = isXiMinus ? negDaughterTrack : posDaughterTrack;

        rXi.fill(HIST("hMassXiSelected"), massXi);
        rXi.fill(HIST("hDCAV0DaughtersXi"), casc.dcaV0daughters());
        rXi.fill(HIST("hV0CosPAXi"), cosPAv0);
        rXi.fill(HIST("hrapidityXi"), casc.yXi());
        rXi.fill(HIST("hctauXi"), ctauXi);
        rXi.fill(HIST("h2DdecayRadiusXi"), casc.cascradius());
        rXi.fill(HIST("hMassXipT"), massXi, casc.pt());
        rXi.fill(HIST("hMassXipTFlat"), massXi, casc.pt(), flattenicity);
        rXi.fill(HIST("hNSigmaProtonFromXi"), protonDaughter.tpcNSigmaPr(), protonDaughter.tpcInnerParam());
        rXi.fill(HIST("hNSigmaPionFromXi"), pionDaughter.tpcNSigmaPi(), pionDaughter.tpcInnerParam());
        rXi.fill(HIST("hNSigmaBachPionFromXi"), bacDaughterTrack.tpcNSigmaPi(), bacDaughterTrack.tpcInnerParam());
      }
    }
  }
  Preslice<LabeledCascades> perColCasc = aod::cascade::collisionId;
  Preslice<LabeledDauTracks> perColDauTracksMC = aod::track::collisionId;
  SliceCache cacheCasc;

  void processRecMCRun3Cascade(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms,
                                         aod::PVMults, aod::McCollisionLabels> const& collisions,
                               LabeledCascades const& Cascades,
                               aod::V0Datas const&, LabeledDauTracks const& tracks,
                               soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
                               aod::FV0As const& /*fv0s*/, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    const bool own = (eventHistOwner == kOwnerRecMCCasc);
    for (const auto& collision : collisions) {
      if (evSel.applyEvSel &&
          !(isEventSelected(collision, own))) { // Checking if the event passes the
                                                // selection criteria
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto vtxZ = collision.posZ();
      auto vtxY = collision.posY();
      auto vtxX = collision.posX();

      auto tracksThisCollision = tracks.sliceBy(perColDauTracksMC, collision.globalIndex());
      float flattenicity = estimateFlattenicity(collision, tracksThisCollision);
      if (flattenicity < 0.f) {
        continue;
      }
      if (own) {
        rEventSelection.fill(HIST("hEventsSelected"), nbinFlattenicity - 0.5);

        rEventSelection.fill(HIST("hVertexZ"), vtxZ);
        rEventSelection.fill(HIST("hFlattenicityDistribution"), flattenicity);
        rEventSelection.fill(HIST("hCentFT0M"), collision.centFT0M());
        rEventSelection.fill(HIST("hCentFT0MvsFlattenicity"), collision.centFT0M(), flattenicity);
        fillFlatDistRec(flattenicity, collision.centFT0M(), collision.isInelGt0());
        fillChargedRec<true>(tracksThisCollision, flattenicity);
      }

      auto cascsThisCollision = Cascades.sliceBy(perColCasc, collision.globalIndex());
      const auto& mcCollision = collision.mcCollision_as<aod::McCollisions>();

      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cacheCasc);
      const float flattenicityMCGen = estimateFlattenicityFV0MC(particlesInCollision, own);
      if (own) {
        rEventSelection.fill(HIST("hFlattenicityDistributionMCGen_Rec"), flattenicityMCGen);
        rEventSelection.fill(HIST("hFlattenicity_Corr_Gen_vs_Rec"), flattenicityMCGen, flattenicity);
        rEventSelection.fill(HIST("hFlatGenVsRecFine"), flattenicityMCGen, flattenicity);
        fillFlatDistGenInRec(flattenicityMCGen, collision.centFT0M(), collision.isInelGt0());
      }

      for (const auto& casc : cascsThisCollision) {

        // MC truth matching, else the spectra keep the combinatorial background
        if (!casc.has_mcParticle()) {
          continue;
        }
        const auto& cascMcParticle = casc.mcParticle();

        auto posDaughterTrack = casc.posTrack_as<LabeledDauTracks>();
        auto negDaughterTrack = casc.negTrack_as<LabeledDauTracks>();
        auto bacDaughterTrack = casc.bachelor_as<LabeledDauTracks>();

        float cascPos = std::hypot(casc.x() - vtxX, casc.y() - vtxY, casc.z() - vtxZ);
        float cascTotMom = std::hypot(casc.px(), casc.py(), casc.pz());
        float ctauXi = o2::constants::physics::MassXiMinus * cascPos / (cascTotMom + 1e-13);
        float cosPAv0 = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
        float massXi = casc.mXi();
        rXi.fill(HIST("hMassXi"), massXi);
        // Cascade
        if (std::abs(cascMcParticle.pdgCode()) == PDG_t::kXiMinus &&
            isSelectedXi(collision, casc, posDaughterTrack, negDaughterTrack, bacDaughterTrack, ctauXi)) {

          const bool isPrimaryCasc = cascMcParticle.isPhysicalPrimary();
          if (!isPrimaryCasc) {
            float motherPt = -1.f;
            const int motherIndex = getFeedDownMother(cascMcParticle, motherPt);
            rXi.fill(HIST("hMassFeedDownXipTFlat"), massXi, casc.pt(), flattenicity);
            rXi.fill(HIST("hFeedDownXiPtVsMotherPt"), casc.pt(), motherPt);
            rXi.fill(HIST("hFeedDownXiMatrix"), casc.pt(), motherPt, flattenicity, motherIndex);
            rXi.fill(HIST("hFeedDownXiMotherPdg"), motherIndex);
          }
          if (!isPrimaryCasc && eventClass.requirePrimaryMC) {
            continue;
          }

          const bool isXiMinus = casc.sign() < 0;
          const auto& protonDaughter = isXiMinus ? posDaughterTrack : negDaughterTrack;
          const auto& pionDaughter = isXiMinus ? negDaughterTrack : posDaughterTrack;

          rXi.fill(HIST("hMassXiSelected"), massXi);
          rXi.fill(HIST("hDCAV0DaughtersXi"), casc.dcaV0daughters());
          rXi.fill(HIST("hV0CosPAXi"), cosPAv0);
          rXi.fill(HIST("hrapidityXi"), casc.yXi());
          rXi.fill(HIST("hctauXi"), ctauXi);
          rXi.fill(HIST("h2DdecayRadiusXi"), casc.cascradius());
          rXi.fill(HIST("hMassXipT"), massXi, casc.pt());
          rXi.fill(HIST("hMassXipTFlat"), massXi, casc.pt(), flattenicity);
          rXi.fill(HIST("hMassXipTTrueFlat"), massXi, casc.pt(), flattenicityMCGen);
          rXi.fill(HIST("hPtResXi"), cascMcParticle.pt(),
                   (casc.pt() - cascMcParticle.pt()) / cascMcParticle.pt());
          rXi.fill(HIST("hNSigmaProtonFromXi"), protonDaughter.tpcNSigmaPr(), protonDaughter.tpcInnerParam());
          rXi.fill(HIST("hNSigmaPionFromXi"), pionDaughter.tpcNSigmaPi(), pionDaughter.tpcInnerParam());
          rXi.fill(HIST("hNSigmaBachPionFromXi"), bacDaughterTrack.tpcNSigmaPi(), bacDaughterTrack.tpcInnerParam());
        }
      }
      const bool isInelGt0Rec = collision.isInelGt0();

      for (const auto& mcParticle : particlesInCollision) {
        if (mcParticle.isPhysicalPrimary() && std::abs(mcParticle.y()) <= cascSel.cascsettingRapidity &&
            std::abs(mcParticle.pdgCode()) == PDG_t::kXiMinus) {
          rXi.fill(HIST("Generated_MCRecoCollCheck_INEL_Xi"), mcParticle.pt(), flattenicity);
          if (isInelGt0Rec) {
            rXi.fill(HIST("Generated_MCRecoCollCheck_INELgt0_Xi"), mcParticle.pt(), flattenicity);
            rXi.fill(HIST("Generated_MCRecoCollCheck_INELgt0_Xi_TrueFlat"), mcParticle.pt(), flattenicityMCGen);
          }
        }
      }
    }
  }

  // ================== Percentile determination pass ====================== //
  //
  // Event selection and flattenicity only, so a short run gives the 1-rho
  // distribution the percentile boundaries are read off. The boundaries then go
  // back into binning.axisFlat for the spectra pass, where one class is exactly
  // one bin.

  void processFlatDistData(
    soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision,
    TrackCandidates const& tracks,
    soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
    aod::FV0As const& /*fv0s*/)
  {
    const bool own = (eventHistOwner == kOwnerFlatDistData);
    if (evSel.applyEvSel && !isEventSelected(collision, own)) {
      return;
    }
    const float flattenicity = estimateFlattenicity(collision, tracks);
    if (flattenicity < 0.f) {
      return;
    }
    if (own) {
      rEventSelection.fill(HIST("hEventsSelected"), nbinFlattenicity - 0.5);
      rEventSelection.fill(HIST("hVertexZ"), collision.posZ());
      rEventSelection.fill(HIST("hFlattenicityDistribution"), flattenicity);
      rEventSelection.fill(HIST("hCentFT0M"), collision.centFT0M());
      rEventSelection.fill(HIST("hCentFT0MvsFlattenicity"), collision.centFT0M(), flattenicity);
      fillFlatDistRec(flattenicity, collision.centFT0M(), collision.isInelGt0());
      fillChargedRec<false>(tracks, flattenicity);
    }
  }

  void processFlatDistMC(
    soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults,
              aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&, TrackCandidatesMC const& tracks,
    soa::Join<aod::BCs, aod::Timestamps> const& /*bcs*/, aod::FT0s const& /*ft0s*/,
    aod::FV0As const& /*fv0s*/, aod::McParticles const& mcParticles)
  {
    const bool own = (eventHistOwner == kOwnerFlatDistMC);
    for (const auto& collision : collisions) {
      if (evSel.applyEvSel && !isEventSelected(collision, own)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto tracksThisCollision = tracks.sliceBy(perColTracksMC, collision.globalIndex());
      const float flattenicity = estimateFlattenicity(collision, tracksThisCollision);
      if (flattenicity < 0.f) {
        continue;
      }
      if (!own) {
        continue;
      }
      rEventSelection.fill(HIST("hEventsSelected"), nbinFlattenicity - 0.5);
      rEventSelection.fill(HIST("hVertexZ"), collision.posZ());
      rEventSelection.fill(HIST("hFlattenicityDistribution"), flattenicity);
      rEventSelection.fill(HIST("hCentFT0M"), collision.centFT0M());
      rEventSelection.fill(HIST("hCentFT0MvsFlattenicity"), collision.centFT0M(), flattenicity);
      fillFlatDistRec(flattenicity, collision.centFT0M(), collision.isInelGt0());
      fillChargedRec<true>(tracksThisCollision, flattenicity);

      const auto& mcCollision = collision.mcCollision_as<aod::McCollisions>();
      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache1);
      const float flattenicityMCGen = estimateFlattenicityFV0MC(particlesInCollision);
      rEventSelection.fill(HIST("hFlattenicityDistributionMCGen_Rec"), flattenicityMCGen);
      rEventSelection.fill(HIST("hFlattenicity_Corr_Gen_vs_Rec"), flattenicityMCGen, flattenicity);
      rEventSelection.fill(HIST("hFlatGenVsRecFine"), flattenicityMCGen, flattenicity);
      fillFlatDistGenInRec(flattenicityMCGen, collision.centFT0M(), collision.isInelGt0());
    }
  }

  PROCESS_SWITCH(Lambdak0sflattenicity, processFlatDistData, "Flattenicity distribution only, data", false);
  PROCESS_SWITCH(Lambdak0sflattenicity, processFlatDistMC, "Flattenicity distribution only, MC", false);
  PROCESS_SWITCH(Lambdak0sflattenicity, processDataRun3LambdaK0s, "Process Run 3 Data LambdaK0s", false);
  PROCESS_SWITCH(Lambdak0sflattenicity, processRecMCLambdaK0s, "Process Run 3 MC reconstructed LambdaK0s", false);
  PROCESS_SWITCH(Lambdak0sflattenicity, processGenMC, "Process Run 3 MC generated", false);
  PROCESS_SWITCH(Lambdak0sflattenicity, processDataRun3Cascade, "Process Run 3 Data Cascade", false);
  PROCESS_SWITCH(Lambdak0sflattenicity, processRecMCRun3Cascade, "Process Run 3 mc Rec Cascade", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Lambdak0sflattenicity>(cfgc)};
}
