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
/// \author Laura Serksnyte, TU München, laura.serksnyte@cern.ch; Anton Riedel, TU München, anton.riedel@cern.ch

#include <Framework/Configurable.h>
#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TMath.h>
#include <iostream>
#include <iterator>
#include <string>

#include "../filterTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "CommonConstants/MathConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace CFTrigger
{
// enums
enum CFThreeBodyTriggers { kPPP,
                           kPPL,
                           kPLL,
                           kLLL,
                           kNThreeBodyTriggers };
enum CFTwoBodyTriggers { kPD,
                         kLD,
                         kNTwoBodyTriggers
};
enum ParticleSpecies {
  kProton,
  kDeuteron,
  kLambda,
  kNParticleSpecies
};
enum V0Daughters {
  kDaughPion,
  kDaughProton,
  kNV0Daughters
};
enum ParticleRejection { kRejProton,
                         kRejPion,
                         kRejElectron,
                         kNParticleRejection
};
enum PIDLimits { kTPCMin,
                 kTPCMax,
                 kTOFMin,
                 kTOFMax,
                 kTPCTOF,
                 kNPIDLimits
};

// For configurable tables
static const std::vector<std::string> CFTriggerNames3B{"ppp", "ppL", "pLL", "LLL"};
static const std::vector<std::string> CFTriggerNames2B{"pd", "Ld"};
static const std::vector<std::string> CFTriggerNamesALL{"ppp", "ppL", "pLL", "LLL", "pd", "Ld"};
static const std::vector<std::string> SpeciesNameAll{"Proton", "Deuteron", "Lambda"};
static const std::vector<std::string> SpeciesName{"Proton", "Deuteron"};
static const std::vector<std::string> SpeciesNameAnti{"AntiProton", "AntiDeuteron"};
static const std::vector<std::string> SpeciesV0DaughterName{"Pion", "Proton"};
static const std::vector<std::string> SpeciesRejectionName{"Proton", "Pion", "Electron"};
static const std::vector<std::string> TPCCutName{"TPC min", "TPC max"};
static const std::vector<std::string> SpeciesAvgName{"Proton", "Deuteron"};
static const std::vector<std::string> TPCTOFAvgName{"TPC Avg", "TOF Avg"};
static const std::vector<std::string> PidCutsName{"TPC min", "TPC max", "TOF min", "TOF max", "TPCTOF max"};
static const std::vector<std::string> PtCutsName{"Pt min", "Pt max", "P thres"};
static const std::vector<std::string> ThreeBodyFilterNames{"PPP", "PPL", "PLL", "LLL"};
static const std::vector<std::string> TwoBodyFilterNames{"PD", "LD"};

static const int nPidRejection = 2;
static const int nTracks = 2;
static const int nPidAvg = 2;
static const int nPidCutsDaughers = 2;
static const int nPtCuts = 3;
static const int nAllTriggers = 6;

static const float pidcutsTable[nTracks][kNPIDLimits]{
  {-6.f, 6.f, -6.f, 6.f, 6.f},
  {-6.f, 6.f, -99.f, 99.f, 99.f}};
static const float pidcutsTableAnti[nTracks][kNPIDLimits]{
  {-6.f, 6.f, -6.f, 6.f, 6.f},
  {-6.f, 6.f, -99.f, 99.f, 99.f}};
static const float pidRejectionTable[kNParticleRejection][nPidRejection]{
  {-2.f, 2.f},
  {-2.f, 2.f}};
static const float pidTPCTOFAvgTable[nPidAvg][nPidAvg]{
  {0.f, 0.f},
  {0.f, 0.f}};
static const float pidcutsV0DaughterTable[kNV0Daughters][nPidCutsDaughers]{
  {-6.f, 6.f},
  {-6.f, 6.f}};
static const float ptcutsTable[kNParticleRejection][nPtCuts]{
  {0.35f, 6.f, 0.75f},
  {0.35f, 1.6f, 99.f},
  {0.35f, 6.f, 99.f}};

static const float triggerSwitches[1][nAllTriggers]{
  {1, 1, 1, 1, 1, 1}};

static const float Q3Limits[1][kNThreeBodyTriggers]{
  {0.6f, 0.6f, 0.6f, 0.6f}};
static const float KstarLimits[1][kNTwoBodyTriggers]{
  {1.2f, 1.2f}};
static const float NClustersMin[1][nTracks]{
  {60.0f, 60.0f}};
} // namespace CFTrigger

namespace o2::aod
{
using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;

using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA,
            aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
            aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
} // namespace o2::aod

struct CFFilter {

  Produces<aod::CFFilters> tags3N;
  Produces<aod::CFFiltersTwoN> tags2N;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // Configs for events
  Configurable<bool> ConfIsRun3{
    "ConfIsRun3",
    true,
    "Is Run3"};

  Configurable<bool> ConfEvtSelectZvtx{
    "ConfEvtSelectZvtx",
    true,
    "Event selection includes max. z-Vertex"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx",
                                  10.f,
                                  "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtOfflineCheck{
    "ConfEvtOfflineCheck",
    false,
    "Evt sel: check for offline selection"};
  Configurable<bool> ConfAutocorRejection{
    "ConfAutocorRejection",
    true,
    "Rejection autocorrelation pL pairs"};

  // Configs for tracks
  Configurable<bool> ConfDeuteronThPVMom{
    "ConfDeuteronThPVMom",
    false,
    "True: use momentum at PV instead of TPCinnerparameter for threshold"};

  Configurable<bool> ConfUseManualPIDproton{
    "ConfUseManualPIDproton",
    false,
    "True: use home-made PID solution for proton "};
  Configurable<std::string> ConfPIDBBProton{
    "ConfPIDBBProton",
    "Users/l/lserksny/PIDProton",
    "Path to the CCDB ocject for proton BB param"};
  Configurable<std::string> ConfPIDBBAntiProton{
    "ConfPIDBBAntiProton",
    "Users/l/lserksny/PIDAntiProton",
    "Path to the CCDB ocject for antiproton BB param"};

  Configurable<bool> ConfUseManualPIDdeuteron{
    "ConfUseManualPIDdeuteron",
    false,
    "True: use home-made PID solution for deuteron "};
  Configurable<std::string> ConfPIDBBDeuteron{
    "ConfPIDBBDeuteron",
    "Users/l/lserksny/PIDDeuteron",
    "Path to the CCDB ocject for Deuteron BB param"};
  Configurable<std::string> ConfPIDBBAntiDeuteron{
    "ConfPIDBBAntiDeuteron",
    "Users/l/lserksny/PIDAntiDeuteron",
    "Path to the CCDB ocject for antiDeuteron BB param"};

  Configurable<bool> ConfUseManualPIDpion{
    "ConfUseManualPIDpion",
    false,
    "True: use home-made PID solution for pions"};
  Configurable<std::string> ConfPIDBBPion{
    "ConfPIDBBPion",
    "Users/l/lserksny/PIDPion",
    "Path to the CCDB ocject for Pion BB param"};
  Configurable<std::string> ConfPIDBBAntiPion{
    "ConfPIDBBAntiPion",
    "Users/l/lserksny/PIDAntiPion",
    "Path to the CCDB ocject for antiPion BB param"};

  Configurable<bool> ConfUseManualPIDel{
    "ConfUseManualPIDel",
    false,
    "True: use home-made PID solution for electron"};
  Configurable<std::string> ConfPIDBBElectron{
    "ConfPIDBBElectron",
    "Users/l/lserksny/PIDElectron",
    "Path to the CCDB ocject for Electron BB param"};
  Configurable<std::string> ConfPIDBBAntiElectron{
    "ConfPIDBBAntiElectron",
    "Users/l/lserksny/PIDAntiElectron",
    "Path to the CCDB ocject for antiElectron BB param"};

  Configurable<bool> ConfUseManualPIDdaughterPion{
    "ConfUseManualPIDdaughterPion",
    false,
    "True: use home-made PID solution for pion from V0"};
  Configurable<bool> ConfUseManualPIDdaughterProton{
    "ConfUseManualPIDdaughterProton",
    false,
    "True: use home-made PID solution for proton from V0"};
  Configurable<bool> ConfRejectNotPropagatedTracks{
    "ConfRejectNotPropagatedTracks",
    false,
    "True: reject not propagated tracks"};
  Configurable<float> ConfTrkEta{
    "ConfTrkEta",
    0.85,
    "Eta"};
  Configurable<LabeledArray<float>> ConfTPCNClustersMin{
    "ConfTPCNClustersMin",
    {CFTrigger::NClustersMin[0], 1, CFTrigger::nTracks, std::vector<std::string>{"TPCNClusMin"}, CFTrigger::SpeciesAvgName},
    "kstar limit for two body trigger"};
  Configurable<float> ConfTrkTPCfCls{
    "ConfTrkTPCfCls",
    0.83,
    "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> ConfTrkTPCcRowsMin{
    "ConfTrkTPCcRowsMin",
    70,
    "Minimum number of crossed TPC rows"};
  Configurable<float> ConfTrkTPCsClsMax{
    "ConfTrkTPCsClsMax",
    160,
    "Maximum number of shared TPC clusters"};
  Configurable<float> ConfTrkITSnclsMin{
    "ConfTrkITSnclsMin",
    0,
    "Minimum number of ITS clusters"};
  Configurable<float> ConfTrkITSnclsIbMin{
    "ConfTrkITSnclsIbMin",
    0,
    "Minimum number of ITS clusters in the inner barrel"};
  Configurable<float> ConfTrkDCAxyMax{
    "ConfTrkDCAxyMax",
    0.15,
    "Maximum DCA_xy"};
  Configurable<float> ConfTrkDCAzMax{
    "ConfTrkDCAzMax",
    0.3,
    "Maximum DCA_z"};
  // Checks taken from global track definition
  Configurable<bool> ConfTrkRequireChi2MaxTPC{
    "ConfTrkRequireChi2MaxTPC", false,
    "True: require max chi2 per TPC cluster"};
  Configurable<bool> ConfTrkRequireChi2MaxITS{
    "ConfTrkRequireChi2MaxITS", false,
    "True: require max chi2 per ITS cluster"};
  Configurable<float>
    ConfTrkMaxChi2PerClusterTPC{
      "ConfTrkMaxChi2PerClusterTPC",
      4.0f,
      "Minimal track selection: max allowed chi2 per TPC cluster"}; // 4.0 is default of
                                                                    // global tracks
                                                                    // on 20.01.2023
  Configurable<float>
    ConfTrkMaxChi2PerClusterITS{
      "ConfTrkMaxChi2PerClusterITS",
      36.0f,
      "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default of
                                                                    // global tracks
                                                                    // on 20.01.2023
  Configurable<bool> ConfTrkTPCRefit{
    "ConfTrkTPCRefit",
    false,
    "True: require TPC refit"};
  Configurable<bool> ConfTrkITSRefit{
    "ConfTrkITSRefit",
    false,
    "True: require ITS refit"};
  // PID selections

  Configurable<LabeledArray<float>> ConfPIDCuts{
    "ConfPIDCuts",
    {CFTrigger::pidcutsTable[0], CFTrigger::nTracks, CFTrigger::kNPIDLimits, CFTrigger::SpeciesName, CFTrigger::PidCutsName},
    "Particle PID selections"};
  Configurable<LabeledArray<float>> ConfPIDCutsAnti{
    "ConfPIDCutsAnti",
    {CFTrigger::pidcutsTableAnti[0], CFTrigger::nTracks, CFTrigger::kNPIDLimits, CFTrigger::SpeciesNameAnti, CFTrigger::PidCutsName},
    "Particle PID selections for antiparticles; perfect case scenario identical to particles"};
  Configurable<LabeledArray<float>> ConfPtCuts{
    "ConfPtCuts",
    {CFTrigger::ptcutsTable[0], CFTrigger::kNParticleRejection, CFTrigger::nPtCuts, CFTrigger::SpeciesNameAll, CFTrigger::PtCutsName},
    "Particle Momentum selections"};
  Configurable<bool> ConfRejectNOTDeuteron{
    "ConfRejectNOTDeuteron",
    false,
    "Reject deuteron candidates if they are compatible with electron, pion, proton"};
  Configurable<LabeledArray<float>> ConfPIDRejection{
    "ConfPIDRejection",
    {CFTrigger::pidRejectionTable[0], CFTrigger::kNParticleRejection, CFTrigger::nPidRejection, CFTrigger::SpeciesRejectionName, CFTrigger::TPCCutName},
    "Particle PID Rejection selections (Deuteron candidates only)"};
  Configurable<LabeledArray<float>> ConfPIDTPCTOFAvg{
    "ConfPIDTPCTOFAvg",
    {CFTrigger::pidTPCTOFAvgTable[0], CFTrigger::nPidAvg, CFTrigger::nPidAvg, CFTrigger::SpeciesAvgName, CFTrigger::TPCTOFAvgName},
    "Average expected nSigma of TPC and TOF, which is substracted in calculation of combined TPC and TOF nSigma"};

  // Configs for V0
  Configurable<float> ConfV0PtMin{
    "ConfV0PtMin",
    0.f,
    "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0DCADaughMax{
    "ConfV0DCADaughMax",
    1.8f,
    "Maximum DCA between the V0 daughters"};
  Configurable<float> ConfV0CPAMin{
    "ConfV0CPAMin",
    0.985f,
    "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{
    "ConfV0TranRadV0Min",
    0.2f,
    "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{
    "ConfV0TranRadV0Max",
    100.f,
    "Maximum transverse radius"};
  Configurable<float> ConfV0DecVtxMax{"ConfV0DecVtxMax",
                                      100.f,
                                      "Maximum distance from primary vertex"};
  Configurable<float> ConfV0InvMassLowLimit{
    "ConfV0InvMassLowLimit",
    1.05,
    "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfV0InvMassUpLimit{
    "ConfV0InvMassUpLimit",
    1.18,
    "Upper limit of the V0 invariant mass"};

  Configurable<bool> ConfV0RejectKaons{"ConfV0RejectKaons",
                                       true,
                                       "Switch to reject kaons"};
  Configurable<float> ConfV0InvKaonMassLowLimit{
    "ConfV0InvKaonMassLowLimit",
    0.49,
    "Lower limit of the V0 invariant mass for Kaon rejection"};
  Configurable<float> ConfV0InvKaonMassUpLimit{
    "ConfV0InvKaonMassUpLimit",
    0.505,
    "Upper limit of the V0 invariant mass for Kaon rejection"};

  // config for V0 daughters
  Configurable<float> ConfDaughEta{
    "ConfDaughEta",
    0.85f,
    "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{
    "ConfDaughTPCnclsMin",
    60.f,
    "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughDCAMin{
    "ConfDaughDCAMin",
    0.04f,
    "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<LabeledArray<float>> ConfDaughPIDCuts{
    "ConfDaughPIDCuts",
    {CFTrigger::pidcutsV0DaughterTable[0], CFTrigger::kNV0Daughters, CFTrigger::nPidCutsDaughers, CFTrigger::SpeciesV0DaughterName, CFTrigger::TPCCutName},
    "PID selections for Lambda daughters"};

  // Trigger selections
  Configurable<LabeledArray<float>> ConfTriggerSwitches{
    "ConfTriggerSwitches",
    {CFTrigger::triggerSwitches[0], 1, CFTrigger::nAllTriggers, std::vector<std::string>{"Switch"}, CFTrigger::CFTriggerNamesALL},
    "Turn on specific trigger"};

  Configurable<LabeledArray<float>> ConfQ3Limits{
    "ConfQ3Limits",
    {CFTrigger::Q3Limits[0], 1, CFTrigger::kNThreeBodyTriggers, std::vector<std::string>{"Limit"}, CFTrigger::ThreeBodyFilterNames},
    "Q3 limits for three body trigger"};

  Configurable<LabeledArray<float>> ConfKstarLimits{
    "ConfKstarLimits",
    {CFTrigger::KstarLimits[0], 1, CFTrigger::kNTwoBodyTriggers, std::vector<std::string>{"Limit"}, CFTrigger::TwoBodyFilterNames},
    "kstar limit for two body trigger"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<double> BBProton, BBAntiproton, BBDeuteron, BBAntideuteron, BBPion, BBAntipion, BBElectron, BBAntielectron;
  void init(o2::framework::InitContext&)
  {

    // init the ccdb
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // global histograms
    registry.add("fProcessedEvents", "CF - event filtered;;events", HistType::kTH1F, {{8, -0.5, 7.5}});
    std::vector<std::string> eventTitles = {"all", "rejected", "ppp", "ppL", "pLL", "LLL", "pD", "LD"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    // event cuts
    registry.add("EventCuts/fMultiplicityBefore", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fMultiplicityAfter", "Multiplicity after event cuts", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fZvtxBefore", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("EventCuts/fZvtxAfter", "Zvtx after event cuts", HistType::kTH1F, {{1000, -15, 15}});

    // all tracks
    registry.add("TrackCuts/TracksBefore/fPtTrackBefore", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/TracksBefore/fEtaTrackBefore", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/TracksBefore/fPhiTrackBefore", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/TracksBefore/fMomCorrelation", "fMomCorrelation", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.0f, 10.0f}}});

    // PID vs momentum before cuts
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPProtonBefore", "NSigmaTPC Proton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTOFvsPProtonBefore", "NSigmaTOF Proton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCTOFvsPProtonBefore", "NSigmaTPCTOF Proton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPAntiProtonBefore", "NSigmaTPC AntiProton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTOFvsPAntiProtonBefore", "NSigmaTOF AntiProton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCTOFvsPAntiProtonBefore", "NSigmaTPCTOF AntiProton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPDeuteronBefore", "NSigmaTPC Deuteron Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPAntiDeuteronBefore", "NSigmaTPC AntiDeuteron Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // TPC signal
    registry.add("TrackCuts/TPCSignal/fTPCSignal", "TPCSignal", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalALLCUTS", "TPCSignalALLCUTS", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    // TPC signal anti
    registry.add("TrackCuts/TPCSignal/fTPCSignalAnti", "TPCSignal", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalAntiALLCUTS", "TPCSignalALLCUTS", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});

    registry.add("TrackCuts/TPCSignal/fTPCSignalProton", "fTPCSignalProton", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalAntiProton", "fTPCSignalAntiProton", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalDeuteron", "fTPCSignalDeuteron", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalAntiDeuteron", "fTPCSignalAntiDeuteron", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalPionMinusV0Daughter", "fTPCSignalPionMinusV0Daughter", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalPionPlusV0Daughter", "fTPCSignalPionPlusV0Daughter", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalProtonMinusV0Daughter", "fTPCSignalProtonMinusV0Daughter", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/TPCSignal/fTPCSignalProtonPlusV0Daughter", "fTPCSignalProtonPlusV0Daughter", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});

    // PID vs momentum before cuts daughters
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPProtonV0DaughBefore", "NSigmaTPC Proton V0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPPionMinusV0DaughBefore", "NSigmaTPC AntiPion V0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPAntiProtonAntiV0DaughBefore", "NSigmaTPC AntiProton antiV0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/NSigmaBefore/fNsigmaTPCvsPPionPlusAntiV0DaughBefore", "NSigmaTPC Pion antiV0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // proton
    // TEST P TPC
    registry.add("TrackCuts/Proton/fPProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/fPTPCProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});

    registry.add("TrackCuts/Proton/fPtProton", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/fEtaProton", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Proton/fPhiProton", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/Proton/fNsigmaTPCvsPProton", "NSigmaTPC Proton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/fNsigmaTOFvsPProton", "NSigmaTOF Proton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/fNsigmaTPCTOFvsPProton", "NSigmaTPCTOF Proton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    registry.add("TrackCuts/Proton/fDCAxyProton", "fDCAxy Proton", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Proton/fDCAzProton", "fDCAz Proton", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Proton/fTPCsClsProton", "fTPCsCls Proton", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Proton/fTPCcRowsProton", "fTPCcRows Proton", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Proton/fTrkTPCfClsProton", "fTrkTPCfCls Proton", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/Proton/fTPCnclsProton", "fTPCncls Proton", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Proton/fNsigmaTPCvsPDProtonP", "NSigmaTPC Proton vd P", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // antiproton
    registry.add("TrackCuts/AntiProton/fPtAntiProton", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/fEtaAntiProton", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiProton/fPhiAntiProton", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProton", "NSigmaTPC AntiProton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProton", "NSigmaTOF AntiProton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProton", "NSigmaTPCTOF AntiProton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    registry.add("TrackCuts/AntiProton/fDCAxyAntiProton", "fDCAxy AntiProton", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiProton/fDCAzAntiProton", "fDCAz AntiProton", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiProton/fTPCsClsAntiProton", "fTPCsCls AntiProton", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiProton/fTPCcRowsAntiProton", "fTPCcRows AntiProton", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiProton/fTrkTPCfClsAntiProton", "fTrkTPCfCls AntiProton", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/AntiProton/fTPCnclsAntiProton", "fTPCncls AntiProton", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // deuteron
    registry.add("TrackCuts/Deuteron/fPtDeuteron", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Deuteron/fEtaDeuteron", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Deuteron/fPhiDeuteron", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteron", "NSigmaTPC Deuteron", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteronP", "NSigmaTPC Deuteron vd P", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    registry.add("TrackCuts/Deuteron/fDCAxyDeuteron", "fDCAxy Deuteron", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Deuteron/fDCAzDeuteron", "fDCAz Deuteron", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/Deuteron/fTPCsClsDeuteron", "fTPCsCls Deuteron", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Deuteron/fTPCcRowsDeuteron", "fTPCcRows Deuteron", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Deuteron/fTrkTPCfClsDeuteron", "fTrkTPCfCls Deuteron", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/Deuteron/fTPCnclsDeuteron", "fTPCncls Deuteron", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // antideuteron
    registry.add("TrackCuts/AntiDeuteron/fPtAntiDeuteron", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiDeuteron/fEtaAntiDeuteron", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiDeuteron/fPhiAntiDeuteron", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteron", "NSigmaTPC AntiDeuteron", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteronP", "NSigmaTPC AntiDeuteron vd P", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    registry.add("TrackCuts/AntiDeuteron/fDCAxyAntiDeuteron", "fDCAxy AntiDeuteron", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiDeuteron/fDCAzAntiDeuteron", "fDCAz AntiDeuteron", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    registry.add("TrackCuts/AntiDeuteron/fTPCsClsAntiDeuteron", "fTPCsCls AntiDeuteron", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiDeuteron/fTPCcRowsAntiDeuteron", "fTPCcRows AntiDeuteron", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiDeuteron/fTrkTPCfClsAntiDeuteron", "fTrkTPCfCls AntiDeuteron", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    registry.add("TrackCuts/AntiDeuteron/fTPCnclsAntiDeuteron", "fTPCncls AntiDeuteron", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // lambda before selections
    registry.add("TrackCuts/V0Before/fPtLambdaBefore", "Transverse momentum of all processed V0s before cuts", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/V0Before/fInvMassLambdaBefore", "Invariant mass of all processed V0s (Lambda) before cuts", HistType::kTH1F, {{1000, 0.7, 1.5}});

    // lambda
    registry.add("TrackCuts/Lambda/fPtLambda", "Transverse momentum of all selected V0s", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Lambda/fInvMassLambda", "Invariant mass of all selected V0s (Lambda)", HistType::kTH1F, {{1000, 0.7, 1.5}});
    registry.add("TrackCuts/Lambda/fInvMassLambdaKaonvsLambda", "Invariant mass of rejected K0 vs V0s (Lambda)", HistType::kTH2F, {{1000, 0.7, 1.5}, {1000, 0.3, 0.6}});
    registry.add("TrackCuts/Lambda/fV0DCADaugh", "V0DCADaugh", HistType::kTH1F, {{1000, -4, 4}});
    registry.add("TrackCuts/Lambda/fV0CPA", "V0 CPA", HistType::kTH1F, {{1000, 0.7, 1}});
    registry.add("TrackCuts/Lambda/fV0TranRad", "V0 TranRad", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("TrackCuts/Lambda/f0DecVtxX", "V0 DecVtxX", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("TrackCuts/Lambda/f0DecVtxY", "V0 DecVtxY", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("TrackCuts/Lambda/f0DecVtxZ", "V0 DecVtxZ", HistType::kTH1F, {{1000, 0, 1000}});

    // Lambda daughter
    registry.add("TrackCuts/Lambda/PosDaughter/Eta", "Lambda Pos Daugh Eta", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Lambda/PosDaughter/DCAXY", "Lambda Pos Daugh DCAXY", HistType::kTH1F, {{1000, -2.5f, 2.5f}});
    registry.add("TrackCuts/Lambda/PosDaughter/fTPCncls", "Lambda Pos Daugh TPCncls", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Lambda/NegDaughter/Eta", "Lambda Neg Daugh Eta", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Lambda/NegDaughter/DCAXY", "Lambda Neg Daugh DCAXY", HistType::kTH1F, {{1000, -2.5f, 2.5f}});
    registry.add("TrackCuts/Lambda/NegDaughter/fTPCncls", "Lambda Neg Daugh TPCncls", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/Lambda/PosDaughter/fNsigmaTPCvsPProtonV0Daugh", "NSigmaTPC Proton V0Daught ", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Lambda/NegDaughter/fNsigmaTPCvsPPionMinusV0Daugh", "NSigmaTPC AntiPion V0Daught ", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // antilambda
    registry.add("TrackCuts/AntiLambda/fPtAntiLambda", "Transverse momentum of all selected V0s", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiLambda/fInvMassAntiLambda", "Invariant mass of all selected V0s (Lambda)", HistType::kTH1F, {{1000, 0.7, 1.5}});
    registry.add("TrackCuts/AntiLambda/fInvMassAntiLambdaKaonvsAntiLambda", "Invariant mass of rejected K0 vs V0s (Lambda)", HistType::kTH2F, {{1000, 0.7, 1.5}, {1000, 0.3, 0.6}});
    registry.add("TrackCuts/AntiLambda/fV0DCADaugh", "V0DCADaugh", HistType::kTH1F, {{1000, -4, 4}});
    registry.add("TrackCuts/AntiLambda/fV0CPA", "V0 CPA", HistType::kTH1F, {{1000, 0.7, 1}});
    registry.add("TrackCuts/AntiLambda/fV0TranRad", "V0 TranRad", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("TrackCuts/AntiLambda/f0DecVtxX", "V0 DecVtxX", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("TrackCuts/AntiLambda/f0DecVtxY", "V0 DecVtxY", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("TrackCuts/AntiLambda/f0DecVtxZ", "V0 DecVtxZ", HistType::kTH1F, {{1000, 0, 1000}});

    // AntiLambda daughter
    registry.add("TrackCuts/AntiLambda/PosDaughter/Eta", "AntiLambda Pos Daugh Eta", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiLambda/PosDaughter/DCAXY", "AntiLambda Pos Daugh DCAXY", HistType::kTH1F, {{1000, -2.5f, 2.5f}});
    registry.add("TrackCuts/AntiLambda/PosDaughter/fTPCncls", "AntiLambda Pos Daugh TPCncls", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiLambda/NegDaughter/Eta", "AntiLambda Neg Daugh Eta", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiLambda/NegDaughter/DCAXY", "AntiLambda Neg Daugh DCAXY", HistType::kTH1F, {{1000, -2.5f, 2.5f}});
    registry.add("TrackCuts/AntiLambda/NegDaughter/fTPCncls", "AntiLambda Neg Daugh TPCncls", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    registry.add("TrackCuts/AntiLambda/NegDaughter/fNsigmaTPCvsPAntiProtonAntiV0Daugh", "NSigmaTPC AntiProton antiV0Daught ", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiLambda/PosDaughter/fNsigmaTPCvsPPionPlusAntiV0Daugh", "NSigmaTPC Pion antiV0Daught ", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // for ppp
    registry.add("ppp/fMultiplicity", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("ppp/fZvtx", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("ppp/fSE_particle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppp/fSE_antiparticle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // for ppl
    registry.add("ppl/fMultiplicity", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("ppl/fZvtx", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("ppl/fSE_particle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ppl/fSE_antiparticle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // for pll
    registry.add("pll/fMultiplicity", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("pll/fZvtx", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("pll/fSE_particle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("pll/fSE_antiparticle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // for pll
    registry.add("lll/fMultiplicity", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("lll/fZvtx", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("lll/fSE_particle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("lll/fSE_antiparticle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // for pd
    registry.add("pd/fMultiplicity", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("pd/fZvtx", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("pd/fSE_particle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("pd/fSE_antiparticle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // for ld
    registry.add("ld/fMultiplicity", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("ld/fZvtx", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("ld/fSE_particle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    registry.add("ld/fSE_antiparticle", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
  }

  float mMassElectron = o2::constants::physics::MassElectron;
  float mMassPion = o2::constants::physics::MassPionCharged;
  float mMassProton = o2::constants::physics::MassProton;
  float mMassLambda = o2::constants::physics::MassLambda;
  float mMassDeuteron = o2::constants::physics::MassDeuteron;
  int currentRunNumber = -999;
  int lastRunNumber = -999;

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (ConfEvtSelectZvtx && std::abs(col.posZ()) > ConfEvtZvtx) {
      return false;
    }
    if (ConfEvtOfflineCheck && !col.sel8()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedTrack(T const& track, CFTrigger::ParticleSpecies partSpecies)
  {
    const auto pT = track.pt();
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto tpcRClsC = track.tpcCrossedRowsOverFindableCls();
    const auto tpcNClsC = track.tpcNClsCrossedRows();
    const auto tpcNClsS = track.tpcNClsShared();
    const auto itsNCls = track.itsNCls();
    const auto itsNClsIB = track.itsNClsInnerBarrel();
    const auto dcaXY = track.dcaXY();
    const auto dcaZ = track.dcaZ();

    if (pT < ConfPtCuts->get(partSpecies, "Pt min")) {
      return false;
    }
    if (pT > ConfPtCuts->get(partSpecies, "Pt max")) {
      return false;
    }
    if (std::abs(eta) > ConfTrkEta) {
      return false;
    }
    if (tpcNClsF < ConfTPCNClustersMin->get("TPCNClusMin", partSpecies)) {
      return false;
    }
    if (tpcRClsC < ConfTrkTPCfCls) {
      return false;
    }
    if (tpcNClsC < ConfTrkTPCcRowsMin) {
      return false;
    }
    if (tpcNClsS > ConfTrkTPCsClsMax) {
      return false;
    }
    if (itsNCls < ConfTrkITSnclsMin) {
      return false;
    }
    if (itsNClsIB < ConfTrkITSnclsIbMin) {
      return false;
    }
    if (std::abs(dcaXY) > ConfTrkDCAxyMax) {
      return false;
    }
    if (std::abs(dcaZ) > ConfTrkDCAzMax) {
      return false;
    }
    // TODO: which dca, put dcaxy for now
    if (ConfRejectNotPropagatedTracks && std::abs(dcaXY) > 1e3) {
      return false;
    }
    if (ConfTrkRequireChi2MaxTPC && track.tpcChi2NCl() >= ConfTrkMaxChi2PerClusterTPC) {
      return false;
    }
    if (ConfTrkRequireChi2MaxITS && track.itsChi2NCl() >= ConfTrkMaxChi2PerClusterITS) {
      return false;
    }
    if (ConfTrkTPCRefit && !track.hasTPC()) {
      return false;
    }
    if (ConfTrkITSRefit && !track.hasITS()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge, CFTrigger::V0Daughters species, double nSigmaTPCDaug[2])
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
    const auto sign = track.sign();
    double nSigmaTPC = -999.f;

    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
    }
    if (std::abs(eta) > ConfDaughEta) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (std::abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }

    switch (species) {
      case CFTrigger::kDaughPion:
        nSigmaTPC = nSigmaTPCDaug[1];
        break;
      case CFTrigger::kDaughProton:
        nSigmaTPC = nSigmaTPCDaug[0];
        break;
      default:
        LOG(fatal) << "Particle species for V0 daughters not found";
    }

    if (nSigmaTPC < ConfDaughPIDCuts->get(species, "TPC min") ||
        nSigmaTPC > ConfDaughPIDCuts->get(species, "TPC max")) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedTrackPID(T const& track, CFTrigger::ParticleSpecies partSpecies, bool Rejection, double nSigmaTPC[2], int charge)
  {
    // nSigma should have entries [proton, deuteron]
    bool isSelected = false;
    bool pThres = true;
    float nSigma = -999.;

    // check momentum threshold
    if (track.tpcInnerParam() <= ConfPtCuts->get(partSpecies, "P thres")) {
      pThres = true;
    } else {
      pThres = false;
    }
    if (CFTrigger::kDeuteron == partSpecies && ConfDeuteronThPVMom) {
      if (track.p() <= ConfPtCuts->get(partSpecies, "P thres")) {
        pThres = true;
      } else {
        pThres = false;
      }
    }

    // compute nsigma
    switch (partSpecies) {
      case CFTrigger::kProton:
        nSigma = (pThres) ? nSigmaTPC[0]
                          : std::sqrt(std::pow(nSigmaTPC[0] - ConfPIDTPCTOFAvg->get("Proton", "TPC Avg"), 2) +
                                      std::pow(track.tofNSigmaPr() - ConfPIDTPCTOFAvg->get("Proton", "TOF Avg"), 2));
        break;
      case CFTrigger::kDeuteron:
        nSigma = (pThres) ? nSigmaTPC[1]
                          : std::sqrt(std::pow(nSigmaTPC[1] - ConfPIDTPCTOFAvg->get("Deuteron", "TPC Avg"), 2) +
                                      std::pow(track.tofNSigmaDe() - ConfPIDTPCTOFAvg->get("Deuteron", "TOF Avg"), 2));
        break;
      case CFTrigger::kLambda:
        LOG(fatal) << "No PID selection for Lambdas";
        break;
      default:
        LOG(fatal) << "Particle species not known";
    }
    // check if track is selected

    auto TPCmin = charge > 0 ? ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMin)
                             : ConfPIDCutsAnti->get(partSpecies, CFTrigger::kTPCMin);

    auto TPCmax = charge > 0 ? ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMax)
                             : ConfPIDCutsAnti->get(partSpecies, CFTrigger::kTPCMax);

    auto TPCTOFmax = charge > 0 ? ConfPIDCuts->get(partSpecies, CFTrigger::kTPCTOF)
                                : ConfPIDCutsAnti->get(partSpecies, CFTrigger::kTPCTOF);

    if (pThres) {
      if (nSigma > TPCmin &&
          nSigma < TPCmax) {
        isSelected = true;
      }
    } else {
      if (nSigma < TPCTOFmax) {
        isSelected = true;
      }
    }
    // for deuterons normally, we want to reject tracks that have a high
    // probablilty of being another particle
    if (Rejection) {
      double nSigmaPi = track.tpcNSigmaPi();
      double nSigmaEl = track.tpcNSigmaEl();
      if (ConfUseManualPIDpion) {
        auto bgScalingPion = 1 / mMassPion; // momentum scaling?
        if (BBPion.size() == 6 && charge > 0)
          nSigmaPi = updatePID(track, bgScalingPion, BBPion);
        if (BBAntipion.size() == 6 && charge < 0)
          nSigmaPi = updatePID(track, bgScalingPion, BBAntipion);
      }
      if (ConfUseManualPIDel) {
        auto bgScalingElectron = 1 / mMassElectron; // momentum scaling?
        if (BBElectron.size() == 6 && charge < 0)
          nSigmaEl = updatePID(track, bgScalingElectron, BBElectron);
        if (BBAntielectron.size() == 6 && charge > 0)
          nSigmaEl = updatePID(track, bgScalingElectron, BBAntielectron);
      }
      if ((ConfPIDRejection->get(CFTrigger::kRejProton, CFTrigger::kTPCMin) < nSigmaTPC[0] &&
           ConfPIDRejection->get(CFTrigger::kRejProton, CFTrigger::kTPCMax) > nSigmaTPC[0]) ||
          (ConfPIDRejection->get(CFTrigger::kRejPion, CFTrigger::kTPCMin) < nSigmaPi &&
           ConfPIDRejection->get(CFTrigger::kRejPion, CFTrigger::kTPCMax) > nSigmaPi) ||
          (ConfPIDRejection->get(CFTrigger::kRejElectron, CFTrigger::kTPCMin) < nSigmaEl &&
           ConfPIDRejection->get(CFTrigger::kRejElectron, CFTrigger::kTPCMax) > nSigmaEl)) {
        return false;
      }
    }
    return isSelected;
  }

  template <typename C, typename V, typename T>
  bool isSelectedMinimalV0(C const& col, V const& v0, T const& posTrack,
                           T const& negTrack, float charge, double nSigmaTPCPos[2], double nSigmaTPCNeg[2])
  {
    const auto signPos = posTrack.sign();
    const auto signNeg = negTrack.sign();
    if (signPos < 0 || signNeg > 0) {
      LOG(info) << "Something wrong in isSelectedMinimal";
      LOG(info) << "ERROR - Wrong sign for V0 daughters";
    }
    const float pT = v0.pt();
    const std::vector<float> decVtx = {v0.x(), v0.y(), v0.z()};
    const float tranRad = v0.v0radius();
    const float dcaDaughv0 = v0.dcaV0daughters();
    const float cpav0 = v0.v0cosPA(col.posX(), col.posY(), col.posZ());

    const float invMassLambda = v0.mLambda();
    const float invMassAntiLambda = v0.mAntiLambda();

    if ((invMassLambda < ConfV0InvMassLowLimit || invMassLambda > ConfV0InvMassUpLimit) &&
        (invMassAntiLambda < ConfV0InvMassLowLimit || invMassAntiLambda > ConfV0InvMassUpLimit)) {
      return false;
    }
    if (ConfV0RejectKaons) {
      const float invMassKaon = v0.mK0Short();
      if (invMassKaon > ConfV0InvKaonMassLowLimit && invMassKaon < ConfV0InvKaonMassUpLimit) {
        return false;
      }
    }
    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    for (size_t i = 0; i < decVtx.size(); i++) {
      if (decVtx.at(i) > ConfV0DecVtxMax) {
        return false;
      }
    }
    if (charge > 0) {
      if (!isSelectedV0Daughter(posTrack, 1, CFTrigger::kDaughProton, nSigmaTPCPos)) {
        return false;
      }
      if (!isSelectedV0Daughter(negTrack, -1, CFTrigger::kDaughPion, nSigmaTPCNeg)) {
        return false;
      }
    }
    if (charge < 0) {
      if (!isSelectedV0Daughter(posTrack, 1, CFTrigger::kDaughPion, nSigmaTPCPos)) {
        return false;
      }
      if (!isSelectedV0Daughter(negTrack, -1, CFTrigger::kDaughProton, nSigmaTPCNeg)) {
        return false;
      }
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
    ROOT::Math::PxPyPzMVector PartOneCMS(part1);
    ROOT::Math::PxPyPzMVector PartTwoCMS(part2);
    const ROOT::Math::Boost boostPRF =
      ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    const ROOT::Math::PxPyPzMVector trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  ROOT::Math::PxPyPzEVector getqij(const ROOT::Math::PtEtaPhiMVector parti,
                                   const ROOT::Math::PtEtaPhiMVector partj)
  {
    ROOT::Math::PxPyPzEVector vecparti(parti);
    ROOT::Math::PxPyPzEVector vecpartj(partj);
    ROOT::Math::PxPyPzEVector trackSum = vecparti + vecpartj;
    ROOT::Math::PxPyPzEVector trackDifference = vecparti - vecpartj;
    float scaling = trackDifference.Dot(trackSum) / trackSum.Dot(trackSum);
    return trackDifference - scaling * trackSum;
  }
  float getQ3(const ROOT::Math::PtEtaPhiMVector part1,
              const ROOT::Math::PtEtaPhiMVector part2,
              const ROOT::Math::PtEtaPhiMVector part3)
  {
    ROOT::Math::PxPyPzEVector q12 = getqij(part1, part2);
    ROOT::Math::PxPyPzEVector q23 = getqij(part2, part3);
    ROOT::Math::PxPyPzEVector q31 = getqij(part3, part1);
    float Q32 = q12.M2() + q23.M2() + q31.M2();
    return sqrt(-Q32);
  }

  std::vector<double> setValuesBB(aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string ccdbPath)
  {
    map<string, string> metadata;
    auto h = ccdbApi.retrieveFromTFileAny<TH1F>(ccdbPath, metadata, bunchCrossing.timestamp());
    // auto h = ccdb->getForTimeStamp<TH1F>(ccdbPath, bunchCrossing.timestamp()); //check if possible to use this without getting fatal
    if (!h) {
      std::vector<double> dummy;
      LOG(info) << "File from CCDB in path " << ccdbPath << " was not found for run " << bunchCrossing.runNumber() << ". Will use default PID task values!";
      return dummy;
    }
    TAxis* axis = h->GetXaxis();
    std::vector<double> v{static_cast<double>(h->GetBinContent(axis->FindBin("bb1"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb2"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb3"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb4"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb5"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("Resolution")))};
    return v;
  }

  template <typename T>
  double updatePID(T const& track, double bgScaling, std::vector<double> BB)
  {
    double expBethe = tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScaling), BB[0], BB[1], BB[2], BB[3], BB[4]);
    double expSigma = expBethe * BB[5];
    return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
  }

  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, o2::aod::V0Datas const& fullV0s)
  {

    if (!ConfIsRun3) {
      LOG(fatal) << "Run 2 processing is not implemented!";
    }
    if (ConfUseManualPIDproton || ConfUseManualPIDdeuteron) {
      currentRunNumber = col.bc_as<aod::BCsWithTimestamps>().runNumber();
      if (currentRunNumber != lastRunNumber) {
        auto bc = col.bc_as<aod::BCsWithTimestamps>();
        if (ConfUseManualPIDproton || ConfUseManualPIDdaughterProton) {
          BBProton = setValuesBB(bc, ConfPIDBBProton);
          BBAntiproton = setValuesBB(bc, ConfPIDBBAntiProton);
        }
        if (ConfUseManualPIDdeuteron) {
          BBDeuteron = setValuesBB(bc, ConfPIDBBDeuteron);
          BBAntideuteron = setValuesBB(bc, ConfPIDBBAntiDeuteron);
        }
        if (ConfUseManualPIDpion || ConfUseManualPIDdaughterPion) {
          BBPion = setValuesBB(bc, ConfPIDBBPion);
          BBAntipion = setValuesBB(bc, ConfPIDBBAntiPion);
        }
        if (ConfUseManualPIDpion) {
          BBElectron = setValuesBB(bc, ConfPIDBBElectron);
          BBAntielectron = setValuesBB(bc, ConfPIDBBAntiElectron);
        }
        lastRunNumber = currentRunNumber;
      }
    }

    registry.fill(HIST("fProcessedEvents"), 0);
    registry.fill(HIST("EventCuts/fMultiplicityBefore"), col.multNTracksPV());
    registry.fill(HIST("EventCuts/fZvtxBefore"), col.posZ());

    bool keepEvent3N[CFTrigger::kNThreeBodyTriggers] = {false, false, false, false};
    int lowQ3Triplets[CFTrigger::kNThreeBodyTriggers] = {0, 0, 0, 0};

    bool keepEvent2N[CFTrigger::kNTwoBodyTriggers] = {false, false};
    int lowKstarPairs[CFTrigger::kNTwoBodyTriggers] = {0, 0};

    if (isSelectedEvent(col)) {

      registry.fill(HIST("EventCuts/fMultiplicityAfter"), col.multNTracksPV());
      registry.fill(HIST("EventCuts/fZvtxAfter"), col.posZ());

      // keep track of proton indices
      std::vector<int> ProtonIndex = {};
      std::vector<int> AntiProtonIndex = {};

      // Prepare vectors for different species
      std::vector<ROOT::Math::PtEtaPhiMVector> protons, antiprotons, deuterons, antideuterons, lambdas, antilambdas;

      // create deuteron and proton vectors (and corresponding antiparticles) for pair and triplet creation
      for (auto& track : tracks) {

        double nTPCSigmaP[2]{track.tpcNSigmaPr(), track.tpcNSigmaDe()};
        double nTPCSigmaN[2]{track.tpcNSigmaPr(), track.tpcNSigmaDe()};

        if (ConfUseManualPIDproton) {
          auto bgScalingProton = 1 / mMassProton; // momentum scaling?
          if (BBProton.size() == 6)
            nTPCSigmaP[0] = updatePID(track, bgScalingProton, BBProton);
          if (BBAntiproton.size() == 6)
            nTPCSigmaN[0] = updatePID(track, bgScalingProton, BBAntiproton);
        }
        if (ConfUseManualPIDdeuteron) {
          auto bgScalingDeuteron = 1 / mMassDeuteron; // momentum scaling?
          if (BBDeuteron.size() == 6)
            nTPCSigmaP[1] = updatePID(track, bgScalingDeuteron, BBDeuteron);
          if (BBAntideuteron.size() == 6)
            nTPCSigmaN[1] = updatePID(track, bgScalingDeuteron, BBAntideuteron);
        }

        registry.fill(HIST("TrackCuts/TracksBefore/fPtTrackBefore"), track.pt());
        registry.fill(HIST("TrackCuts/TracksBefore/fEtaTrackBefore"), track.eta());
        registry.fill(HIST("TrackCuts/TracksBefore/fPhiTrackBefore"), track.phi());

        if (track.sign() > 0) {
          // Fill PID info
          registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignal"), track.tpcInnerParam(), track.tpcSignal());
          if (isSelectedTrack(track, CFTrigger::kProton)) {
            registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalALLCUTS"), track.tpcInnerParam(), track.tpcSignal());
          }
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPProtonBefore"), track.tpcInnerParam(), nTPCSigmaP[0]);
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTOFvsPProtonBefore"), track.tpcInnerParam(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCTOFvsPProtonBefore"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaP[0] - ConfPIDTPCTOFAvg->get("Proton", "TPC Avg"), 2) + std::pow(track.tofNSigmaPr() - ConfPIDTPCTOFAvg->get("Proton", "TOF Avg"), 2)));
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPDeuteronBefore"), track.tpcInnerParam(), nTPCSigmaP[1]);
          registry.fill(HIST("TrackCuts/TracksBefore/fMomCorrelation"), track.p(), track.tpcInnerParam());
        }
        if (track.sign() < 0) {

          registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAnti"), track.tpcInnerParam(), track.tpcSignal());
          if (isSelectedTrack(track, CFTrigger::kProton)) {
            registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiALLCUTS"), track.tpcInnerParam(), track.tpcSignal());
          }

          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPAntiProtonBefore"), track.tpcInnerParam(), nTPCSigmaN[0]);
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTOFvsPAntiProtonBefore"), track.tpcInnerParam(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCTOFvsPAntiProtonBefore"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaN[0] - ConfPIDTPCTOFAvg->get("Proton", "TPC Avg"), 2) + std::pow(track.tofNSigmaPr() - ConfPIDTPCTOFAvg->get("Proton", "TOF Avg"), 2)));
          registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPAntiDeuteronBefore"), track.tpcInnerParam(), nTPCSigmaN[1]);
        }

        // get protons
        if (isSelectedTrack(track, CFTrigger::kProton)) {
          ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), mMassProton);
          if (track.sign() > 0 && isSelectedTrackPID(track, CFTrigger::kProton, false, nTPCSigmaP, 1)) {
            protons.push_back(temp);
            ProtonIndex.push_back(track.globalIndex());
            registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalProton"), track.tpcInnerParam(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/Proton/fPProton"), track.p());
            registry.fill(HIST("TrackCuts/Proton/fPTPCProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/Proton/fPtProton"), track.pt());
            registry.fill(HIST("TrackCuts/Proton/fEtaProton"), track.eta());
            registry.fill(HIST("TrackCuts/Proton/fPhiProton"), track.phi());
            registry.fill(HIST("TrackCuts/Proton/fNsigmaTPCvsPProton"), track.tpcInnerParam(), nTPCSigmaP[0]);
            registry.fill(HIST("TrackCuts/Proton/fNsigmaTOFvsPProton"), track.tpcInnerParam(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/Proton/fNsigmaTPCTOFvsPProton"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaP[0] - ConfPIDTPCTOFAvg->get("Proton", "TPC Avg"), 2) + std::pow(track.tofNSigmaPr() - ConfPIDTPCTOFAvg->get("Proton", "TOF Avg"), 2)));

            registry.fill(HIST("TrackCuts/Proton/fDCAxyProton"), track.dcaXY());
            registry.fill(HIST("TrackCuts/Proton/fDCAzProton"), track.dcaZ());
            registry.fill(HIST("TrackCuts/Proton/fTPCsClsProton"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/Proton/fTPCcRowsProton"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/Proton/fTrkTPCfClsProton"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/Proton/fTPCnclsProton"), track.tpcNClsFound());
          }
          if (track.sign() < 0 && isSelectedTrackPID(track, CFTrigger::kProton, false, nTPCSigmaN, -1)) {
            antiprotons.push_back(temp);
            AntiProtonIndex.push_back(track.globalIndex());
            registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiProton"), track.tpcInnerParam(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/AntiProton/fPtAntiProton"), track.pt());
            registry.fill(HIST("TrackCuts/AntiProton/fEtaAntiProton"), track.eta());
            registry.fill(HIST("TrackCuts/AntiProton/fPhiAntiProton"), track.phi());
            registry.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProton"), track.tpcInnerParam(), nTPCSigmaN[0]);
            registry.fill(HIST("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProton"), track.tpcInnerParam(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProton"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaN[0] - ConfPIDTPCTOFAvg->get("Proton", "TPC Avg"), 2) + std::pow(track.tofNSigmaPr() - ConfPIDTPCTOFAvg->get("Proton", "TOF Avg"), 2)));

            registry.fill(HIST("TrackCuts/AntiProton/fDCAxyAntiProton"), track.dcaXY());
            registry.fill(HIST("TrackCuts/AntiProton/fDCAzAntiProton"), track.dcaZ());
            registry.fill(HIST("TrackCuts/AntiProton/fTPCsClsAntiProton"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/AntiProton/fTPCcRowsAntiProton"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/AntiProton/fTrkTPCfClsAntiProton"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/AntiProton/fTPCnclsAntiProton"), track.tpcNClsFound());
          }
        }
        // get deuterons
        if (isSelectedTrack(track, CFTrigger::kDeuteron)) {
          ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), mMassDeuteron);
          if (track.sign() > 0 && isSelectedTrackPID(track, CFTrigger::kDeuteron, ConfRejectNOTDeuteron.value, nTPCSigmaP, 1)) {
            deuterons.push_back(temp);

            registry.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteronP"), track.p(), nTPCSigmaP[1]);
            registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalDeuteron"), track.tpcInnerParam(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/Deuteron/fPtDeuteron"), track.pt());
            registry.fill(HIST("TrackCuts/Deuteron/fEtaDeuteron"), track.eta());
            registry.fill(HIST("TrackCuts/Deuteron/fPhiDeuteron"), track.phi());
            registry.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteron"), track.tpcInnerParam(), nTPCSigmaP[1]);

            registry.fill(HIST("TrackCuts/Deuteron/fDCAxyDeuteron"), track.dcaXY());
            registry.fill(HIST("TrackCuts/Deuteron/fDCAzDeuteron"), track.dcaZ());
            registry.fill(HIST("TrackCuts/Deuteron/fTPCsClsDeuteron"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/Deuteron/fTPCcRowsDeuteron"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/Deuteron/fTrkTPCfClsDeuteron"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/Deuteron/fTPCnclsDeuteron"), track.tpcNClsFound());
          }
          if (track.sign() < 0 && isSelectedTrackPID(track, CFTrigger::kDeuteron, ConfRejectNOTDeuteron.value, nTPCSigmaN, -1)) {
            antideuterons.push_back(temp);
            registry.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteronP"), track.p(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiDeuteron"), track.tpcInnerParam(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fPtAntiDeuteron"), track.pt());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fEtaAntiDeuteron"), track.eta());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fPhiAntiDeuteron"), track.phi());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteron"), track.tpcInnerParam(), nTPCSigmaN[1]);

            registry.fill(HIST("TrackCuts/AntiDeuteron/fDCAxyAntiDeuteron"), track.dcaXY());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fDCAzAntiDeuteron"), track.dcaZ());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fTPCsClsAntiDeuteron"), track.tpcNClsShared());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fTPCcRowsAntiDeuteron"), track.tpcNClsCrossedRows());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fTrkTPCfClsAntiDeuteron"), track.tpcCrossedRowsOverFindableCls());
            registry.fill(HIST("TrackCuts/AntiDeuteron/fTPCnclsAntiDeuteron"), track.tpcNClsFound());
          }
        }
      }

      // keep track of daugher indices to avoid selfcorrelations
      std::vector<int> LambdaPosDaughIndex = {};
      std::vector<int> LambdaNegDaughIndex = {};
      std::vector<int> AntiLambdaPosDaughIndex = {};
      std::vector<int> AntiLambdaNegDaughIndex = {};

      for (auto& v0 : fullV0s) {

        registry.fill(HIST("TrackCuts/V0Before/fPtLambdaBefore"), v0.pt());
        registry.fill(HIST("TrackCuts/V0Before/fInvMassLambdaBefore"), v0.mLambda());

        auto postrack = v0.template posTrack_as<aod::FemtoFullTracks>();
        auto negtrack = v0.template negTrack_as<aod::FemtoFullTracks>();
        double nTPCSigmaPos[2]{postrack.tpcNSigmaPr(), postrack.tpcNSigmaPi()};
        double nTPCSigmaNeg[2]{negtrack.tpcNSigmaPr(), negtrack.tpcNSigmaPi()};
        if (ConfUseManualPIDdaughterPion) {
          auto bgScalingPion = 1 / mMassPion; // momentum scaling?
          if (BBPion.size() == 6)
            nTPCSigmaPos[1] = updatePID(postrack, bgScalingPion, BBPion);
          if (BBAntipion.size() == 6)
            nTPCSigmaNeg[1] = updatePID(negtrack, bgScalingPion, BBAntipion);
        }
        if (ConfUseManualPIDdaughterProton) {
          auto bgScalingProton = 1 / mMassProton; // momentum scaling?
          if (BBProton.size() == 6)
            nTPCSigmaPos[0] = updatePID(postrack, bgScalingProton, BBProton);
          if (BBAntiproton.size() == 6)
            nTPCSigmaNeg[0] = updatePID(negtrack, bgScalingProton, BBAntiproton);
        }

        registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPProtonV0DaughBefore"), postrack.tpcInnerParam(), nTPCSigmaPos[0]);
        registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPPionPlusAntiV0DaughBefore"), postrack.tpcInnerParam(), nTPCSigmaNeg[1]);

        registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPPionMinusV0DaughBefore"), negtrack.tpcInnerParam(), nTPCSigmaNeg[1]);
        registry.fill(HIST("TrackCuts/NSigmaBefore/fNsigmaTPCvsPAntiProtonAntiV0DaughBefore"), negtrack.tpcInnerParam(), nTPCSigmaNeg[0]);

        if (isSelectedMinimalV0(col, v0, postrack, negtrack, 1, nTPCSigmaPos, nTPCSigmaNeg)) {
          ROOT::Math::PtEtaPhiMVector temp(v0.pt(), v0.eta(), v0.phi(), mMassLambda);
          lambdas.push_back(temp);
          LambdaPosDaughIndex.push_back(postrack.globalIndex());
          LambdaNegDaughIndex.push_back(negtrack.globalIndex());
          registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalProtonPlusV0Daughter"), postrack.tpcInnerParam(), postrack.tpcSignal());
          registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalPionMinusV0Daughter"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
          registry.fill(HIST("TrackCuts/Lambda/fPtLambda"), v0.pt());
          registry.fill(HIST("TrackCuts/Lambda/fInvMassLambda"), v0.mLambda());
          registry.fill(HIST("TrackCuts/Lambda/fInvMassLambdaKaonvsLambda"), v0.mK0Short(), v0.mLambda());
          registry.fill(HIST("TrackCuts/Lambda/fV0DCADaugh"), v0.dcaV0daughters());
          registry.fill(HIST("TrackCuts/Lambda/fV0CPA"), v0.v0cosPA(col.posX(), col.posY(), col.posZ()));
          registry.fill(HIST("TrackCuts/Lambda/fV0TranRad"), v0.v0radius());
          registry.fill(HIST("TrackCuts/Lambda/f0DecVtxX"), v0.x());
          registry.fill(HIST("TrackCuts/Lambda/f0DecVtxY"), v0.y());
          registry.fill(HIST("TrackCuts/Lambda/f0DecVtxZ"), v0.z());

          registry.fill(HIST("TrackCuts/Lambda/PosDaughter/Eta"), postrack.eta());
          registry.fill(HIST("TrackCuts/Lambda/PosDaughter/DCAXY"), postrack.dcaXY());
          registry.fill(HIST("TrackCuts/Lambda/PosDaughter/fTPCncls"), postrack.tpcNClsFound());
          registry.fill(HIST("TrackCuts/Lambda/NegDaughter/Eta"), negtrack.eta());
          registry.fill(HIST("TrackCuts/Lambda/NegDaughter/DCAXY"), negtrack.dcaXY());
          registry.fill(HIST("TrackCuts/Lambda/NegDaughter/fTPCncls"), negtrack.tpcNClsFound());
          registry.fill(HIST("TrackCuts/Lambda/PosDaughter/fNsigmaTPCvsPProtonV0Daugh"), postrack.tpcInnerParam(), nTPCSigmaPos[0]);
          registry.fill(HIST("TrackCuts/Lambda/NegDaughter/fNsigmaTPCvsPPionMinusV0Daugh"), negtrack.tpcInnerParam(), nTPCSigmaNeg[1]);
        }
        if (isSelectedMinimalV0(col, v0, postrack, negtrack, -1, nTPCSigmaPos, nTPCSigmaNeg)) {
          ROOT::Math::PtEtaPhiMVector temp(v0.pt(), v0.eta(), v0.phi(), mMassLambda);
          antilambdas.push_back(temp);
          AntiLambdaPosDaughIndex.push_back(postrack.globalIndex());
          AntiLambdaNegDaughIndex.push_back(negtrack.globalIndex());
          registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalPionPlusV0Daughter"), postrack.tpcInnerParam(), postrack.tpcSignal());
          registry.fill(HIST("TrackCuts/TPCSignal/fTPCSignalProtonMinusV0Daughter"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
          registry.fill(HIST("TrackCuts/AntiLambda/fPtAntiLambda"), v0.pt());
          registry.fill(HIST("TrackCuts/AntiLambda/fInvMassAntiLambda"), v0.mAntiLambda());
          registry.fill(HIST("TrackCuts/AntiLambda/fInvMassAntiLambdaKaonvsAntiLambda"), v0.mK0Short(), v0.mAntiLambda());
          registry.fill(HIST("TrackCuts/AntiLambda/fV0DCADaugh"), v0.dcaV0daughters());
          registry.fill(HIST("TrackCuts/AntiLambda/fV0CPA"), v0.v0cosPA(col.posX(), col.posY(), col.posZ()));
          registry.fill(HIST("TrackCuts/AntiLambda/fV0TranRad"), v0.v0radius());
          registry.fill(HIST("TrackCuts/AntiLambda/f0DecVtxX"), v0.x());
          registry.fill(HIST("TrackCuts/AntiLambda/f0DecVtxY"), v0.y());
          registry.fill(HIST("TrackCuts/AntiLambda/f0DecVtxZ"), v0.z());

          registry.fill(HIST("TrackCuts/AntiLambda/PosDaughter/Eta"), postrack.eta());
          registry.fill(HIST("TrackCuts/AntiLambda/PosDaughter/DCAXY"), postrack.dcaXY());
          registry.fill(HIST("TrackCuts/AntiLambda/PosDaughter/fTPCncls"), postrack.tpcNClsFound());
          registry.fill(HIST("TrackCuts/AntiLambda/NegDaughter/Eta"), negtrack.eta());
          registry.fill(HIST("TrackCuts/AntiLambda/NegDaughter/DCAXY"), negtrack.dcaXY());
          registry.fill(HIST("TrackCuts/AntiLambda/NegDaughter/fTPCncls"), negtrack.tpcNClsFound());
          registry.fill(HIST("TrackCuts/AntiLambda/NegDaughter/fNsigmaTPCvsPAntiProtonAntiV0Daugh"), negtrack.tpcInnerParam(), nTPCSigmaNeg[0]);
          registry.fill(HIST("TrackCuts/AntiLambda/PosDaughter/fNsigmaTPCvsPPionPlusAntiV0Daugh"), postrack.tpcInnerParam(), nTPCSigmaPos[1]);
        }
      }

      float Q3 = 999.f, kstar = 999.f;
      // if(ConfTriggerSwitches->get(static_cast<uint>(0), CFTrigger::)>0.){
      if (ConfTriggerSwitches->get("Switch", "ppp") > 0.) {
        // ppp trigger
        for (auto iProton1 = protons.begin(); iProton1 != protons.end(); ++iProton1) {
          auto iProton2 = iProton1 + 1;
          for (; iProton2 != protons.end(); ++iProton2) {
            auto iProton3 = iProton2 + 1;
            for (; iProton3 != protons.end(); ++iProton3) {
              Q3 = getQ3(*iProton1, *iProton2, *iProton3);
              registry.fill(HIST("ppp/fSE_particle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kPPP)) {
                lowQ3Triplets[CFTrigger::kPPP] += 1;
              }
            }
          }
        }
        for (auto iAntiProton1 = antiprotons.begin(); iAntiProton1 != antiprotons.end(); ++iAntiProton1) {
          auto iAntiProton2 = iAntiProton1 + 1;
          for (; iAntiProton2 != antiprotons.end(); ++iAntiProton2) {
            auto iAntiProton3 = iAntiProton2 + 1;
            for (; iAntiProton3 != antiprotons.end(); ++iAntiProton3) {
              Q3 = getQ3(*iAntiProton1, *iAntiProton2, *iAntiProton3);
              registry.fill(HIST("ppp/fSE_antiparticle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kPPP)) {
                lowQ3Triplets[CFTrigger::kPPP] += 1;
              }
            }
          }
        }
      }
      if (ConfTriggerSwitches->get("Switch", "ppL") > 0.) {
        // ppl trigger
        for (auto iProton1 = protons.begin(); iProton1 != protons.end(); ++iProton1) {
          auto iProton2 = iProton1 + 1;
          auto i1 = std::distance(protons.begin(), iProton1);
          for (; iProton2 != protons.end(); ++iProton2) {
            auto i2 = std::distance(protons.begin(), iProton2);
            for (auto iLambda1 = lambdas.begin(); iLambda1 != lambdas.end(); ++iLambda1) {
              auto i3 = std::distance(lambdas.begin(), iLambda1);
              if (ConfAutocorRejection.value &&
                  (ProtonIndex.at(i1) == LambdaPosDaughIndex.at(i3) ||
                   ProtonIndex.at(i2) == LambdaPosDaughIndex.at(i3))) {
                continue;
              }
              Q3 = getQ3(*iProton1, *iProton2, *iLambda1);
              registry.fill(HIST("ppl/fSE_particle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kPPL)) {
                lowQ3Triplets[CFTrigger::kPPL] += 1;
              }
            }
          }
        }
        for (auto iAntiProton1 = antiprotons.begin(); iAntiProton1 != antiprotons.end(); ++iAntiProton1) {
          auto iAntiProton2 = iAntiProton1 + 1;
          auto i1 = std::distance(antiprotons.begin(), iAntiProton1);
          for (; iAntiProton2 != antiprotons.end(); ++iAntiProton2) {
            auto i2 = std::distance(antiprotons.begin(), iAntiProton2);
            for (auto iAntiLambda1 = antilambdas.begin(); iAntiLambda1 != antilambdas.end(); ++iAntiLambda1) {
              auto i3 = std::distance(antilambdas.begin(), iAntiLambda1);
              if (ConfAutocorRejection.value &&
                  (AntiProtonIndex.at(i1) == AntiLambdaNegDaughIndex.at(i3) ||
                   AntiProtonIndex.at(i2) == AntiLambdaNegDaughIndex.at(i3))) {
                continue;
              }
              Q3 = getQ3(*iAntiProton1, *iAntiProton2, *iAntiLambda1);
              registry.fill(HIST("ppl/fSE_antiparticle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kPPL)) {
                lowQ3Triplets[CFTrigger::kPPL] += 1;
              }
            }
          }
        }
      }
      if (ConfTriggerSwitches->get("Switch", "pLL") > 0.) {
        // pll trigger
        for (auto iLambda1 = lambdas.begin(); iLambda1 != lambdas.end(); ++iLambda1) {
          auto iLambda2 = iLambda1 + 1;
          auto i1 = std::distance(lambdas.begin(), iLambda1);
          for (; iLambda2 != lambdas.end(); ++iLambda2) {
            auto i2 = std::distance(lambdas.begin(), iLambda2);
            if (ConfAutocorRejection.value &&
                (LambdaPosDaughIndex.at(i1) == LambdaPosDaughIndex.at(i2) ||
                 LambdaNegDaughIndex.at(i1) == LambdaNegDaughIndex.at(i2))) {
              continue;
            }
            for (auto iProton1 = protons.begin(); iProton1 != protons.end(); ++iProton1) {
              auto i3 = std::distance(protons.begin(), iProton1);
              if (ConfAutocorRejection.value &&
                  (LambdaPosDaughIndex.at(i1) == ProtonIndex.at(i3) ||
                   LambdaPosDaughIndex.at(i2) == ProtonIndex.at(i3))) {
                continue;
              }
              Q3 = getQ3(*iLambda1, *iLambda2, *iProton1);
              registry.fill(HIST("pll/fSE_particle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kPLL)) {
                lowQ3Triplets[CFTrigger::kPLL] += 1;
              }
            }
          }
        }
        for (auto iAntiLambda1 = antilambdas.begin(); iAntiLambda1 != antilambdas.end(); ++iAntiLambda1) {
          auto iAntiLambda2 = iAntiLambda1 + 1;
          auto i1 = std::distance(antilambdas.begin(), iAntiLambda1);
          for (; iAntiLambda2 != antilambdas.end(); ++iAntiLambda2) {
            auto i2 = std::distance(antilambdas.begin(), iAntiLambda2);
            if (ConfAutocorRejection.value &&
                (AntiLambdaPosDaughIndex.at(i1) == AntiLambdaPosDaughIndex.at(i2) ||
                 AntiLambdaNegDaughIndex.at(i1) == AntiLambdaNegDaughIndex.at(i2))) {
              continue;
            }
            for (auto iAntiProton1 = antiprotons.begin(); iAntiProton1 != antiprotons.end(); ++iAntiProton1) {
              auto i3 = std::distance(antiprotons.begin(), iAntiProton1);
              if (ConfAutocorRejection.value &&
                  (AntiLambdaNegDaughIndex.at(i1) == AntiProtonIndex.at(i3) ||
                   AntiLambdaNegDaughIndex.at(i2) == AntiProtonIndex.at(i3))) {
                continue;
              }
              Q3 = getQ3(*iAntiLambda1, *iAntiLambda2, *iAntiProton1);
              registry.fill(HIST("pll/fSE_antiparticle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kPLL)) {
                lowQ3Triplets[CFTrigger::kPLL] += 1;
              }
            }
          }
        }
      }
      if (ConfTriggerSwitches->get("Switch", "LLL") > 0.) {
        // lll trigger
        for (auto iLambda1 = lambdas.begin(); iLambda1 != lambdas.end(); ++iLambda1) {
          auto iLambda2 = iLambda1 + 1;
          auto i1 = std::distance(lambdas.begin(), iLambda1);
          for (; iLambda2 != lambdas.end(); ++iLambda2) {
            auto i2 = std::distance(lambdas.begin(), iLambda2);
            if (ConfAutocorRejection.value &&
                (LambdaPosDaughIndex.at(i1) == LambdaPosDaughIndex.at(i2) ||
                 LambdaNegDaughIndex.at(i1) == LambdaNegDaughIndex.at(i2))) {
              continue;
            }
            auto iLambda3 = iLambda2 + 1;
            for (; iLambda3 != lambdas.end(); ++iLambda3) {
              auto i3 = std::distance(lambdas.begin(), iLambda3);
              if (ConfAutocorRejection.value &&
                  (LambdaPosDaughIndex.at(i1) == LambdaPosDaughIndex.at(i3) ||
                   LambdaNegDaughIndex.at(i1) == LambdaNegDaughIndex.at(i3) ||
                   LambdaPosDaughIndex.at(i2) == LambdaPosDaughIndex.at(i3) ||
                   LambdaNegDaughIndex.at(i2) == LambdaNegDaughIndex.at(i3))) {
                continue;
              }
              Q3 = getQ3(*iLambda1, *iLambda2, *iLambda3);
              registry.fill(HIST("lll/fSE_particle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kLLL)) {
                lowQ3Triplets[CFTrigger::kLLL] += 1;
              }
            }
          }
        }
        for (auto iAntiLambda1 = antilambdas.begin(); iAntiLambda1 != antilambdas.end(); ++iAntiLambda1) {
          auto iAntiLambda2 = iAntiLambda1 + 1;
          auto i1 = std::distance(antilambdas.begin(), iAntiLambda1);
          for (; iAntiLambda2 != antilambdas.end(); ++iAntiLambda2) {
            auto i2 = std::distance(antilambdas.begin(), iAntiLambda2);
            if (ConfAutocorRejection.value &&
                (AntiLambdaPosDaughIndex.at(i1) == AntiLambdaPosDaughIndex.at(i2) ||
                 AntiLambdaNegDaughIndex.at(i1) == AntiLambdaNegDaughIndex.at(i2))) {
              continue;
            }
            auto iAntiLambda3 = iAntiLambda2 + 1;
            for (; iAntiLambda3 != antilambdas.end(); ++iAntiLambda3) {
              auto i3 = std::distance(antilambdas.begin(), iAntiLambda3);
              if (ConfAutocorRejection.value &&
                  (AntiLambdaPosDaughIndex.at(i1) == AntiLambdaPosDaughIndex.at(i3) ||
                   AntiLambdaNegDaughIndex.at(i1) == AntiLambdaNegDaughIndex.at(i3) ||
                   AntiLambdaPosDaughIndex.at(i2) == AntiLambdaPosDaughIndex.at(i3) ||
                   AntiLambdaNegDaughIndex.at(i2) == AntiLambdaNegDaughIndex.at(i3))) {
                continue;
              }
              Q3 = getQ3(*iAntiLambda1, *iAntiLambda2, *iAntiLambda3);
              registry.fill(HIST("lll/fSE_antiparticle"), Q3);
              if (Q3 < ConfQ3Limits->get(static_cast<uint>(0), CFTrigger::kLLL)) {
                lowQ3Triplets[CFTrigger::kLLL] += 1;
              }
            }
          }
        }
      }
      if (ConfTriggerSwitches->get("Switch", "pd") > 0.) {
        // pd trigger
        for (auto iProton = protons.begin(); iProton != protons.end(); ++iProton) {
          for (auto iDeuteron = deuterons.begin(); iDeuteron != deuterons.end(); ++iDeuteron) {
            kstar = getkstar(*iProton, *iDeuteron);
            registry.fill(HIST("pd/fSE_particle"), kstar);
            if (kstar < ConfKstarLimits->get(static_cast<uint>(0), CFTrigger::kPD)) {
              lowKstarPairs[CFTrigger::kPD] += 1;
            }
          }
        }
        for (auto iAntiProton = antiprotons.begin(); iAntiProton != antiprotons.end(); ++iAntiProton) {
          for (auto iAntiDeuteron = antideuterons.begin(); iAntiDeuteron != antideuterons.end(); ++iAntiDeuteron) {
            kstar = getkstar(*iAntiProton, *iAntiDeuteron);
            registry.fill(HIST("pd/fSE_antiparticle"), kstar);
            if (kstar < ConfKstarLimits->get(static_cast<uint>(0), CFTrigger::kPD)) {
              lowKstarPairs[CFTrigger::kPD] += 1;
            }
          }
        }
      }
      if (ConfTriggerSwitches->get("Switch", "Ld") > 0.) {
        // ld trigger
        for (auto iDeuteron = deuterons.begin(); iDeuteron != deuterons.end(); ++iDeuteron) {
          for (auto iLambda = lambdas.begin(); iLambda != lambdas.end(); ++iLambda) {
            kstar = getkstar(*iDeuteron, *iLambda);
            if (kstar < ConfKstarLimits->get(static_cast<uint>(0), CFTrigger::kLD)) {
              registry.fill(HIST("ld/fSE_particle"), kstar);
              lowKstarPairs[CFTrigger::kLD] += 1;
            }
          }
        }
        for (auto iAntiDeuteron = antideuterons.begin(); iAntiDeuteron != antideuterons.end(); ++iAntiDeuteron) {
          for (auto iAntiLambda = antilambdas.begin(); iAntiLambda != antilambdas.end(); ++iAntiLambda) {
            kstar = getkstar(*iAntiDeuteron, *iAntiLambda);
            if (kstar < ConfKstarLimits->get(static_cast<uint>(0), CFTrigger::kLD)) {
              registry.fill(HIST("ld/fSE_antiparticle"), kstar);
              lowKstarPairs[CFTrigger::kLD] += 1;
            }
          }
        }
      }

    } // if(isSelectedEvent)

    // create tags for three body triggers
    if (lowQ3Triplets[CFTrigger::kPPP] > 0) {
      keepEvent3N[CFTrigger::kPPP] = true;
      registry.fill(HIST("fProcessedEvents"), 2);
      registry.fill(HIST("ppp/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("ppp/fZvtx"), col.posZ());
    }
    if (lowQ3Triplets[CFTrigger::kPPL] > 0) {
      keepEvent3N[CFTrigger::kPPL] = true;
      registry.fill(HIST("fProcessedEvents"), 3);
      registry.fill(HIST("ppl/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("ppl/fZvtx"), col.posZ());
    }
    if (lowQ3Triplets[CFTrigger::kPLL] > 0) {
      keepEvent3N[CFTrigger::kPLL] = true;
      registry.fill(HIST("fProcessedEvents"), 4);
      registry.fill(HIST("pll/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("pll/fZvtx"), col.posZ());
    }
    if (lowQ3Triplets[CFTrigger::kLLL] > 0) {
      keepEvent3N[CFTrigger::kLLL] = true;
      registry.fill(HIST("fProcessedEvents"), 5);
      registry.fill(HIST("lll/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("lll/fZvtx"), col.posZ());
    }

    tags3N(keepEvent3N[CFTrigger::kPPP],
           keepEvent3N[CFTrigger::kPPL],
           keepEvent3N[CFTrigger::kPLL],
           keepEvent3N[CFTrigger::kLLL]);

    // create tags for two body triggers
    if (lowKstarPairs[CFTrigger::kPD] > 0) {
      keepEvent2N[CFTrigger::kPD] = true;
      registry.fill(HIST("fProcessedEvents"), 6);
      registry.fill(HIST("pd/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("pd/fZvtx"), col.posZ());
    }
    if (lowKstarPairs[CFTrigger::kLD] > 0) {
      keepEvent2N[CFTrigger::kLD] = true;
      registry.fill(HIST("fProcessedEvents"), 7);
      registry.fill(HIST("ld/fMultiplicity"), col.multNTracksPV());
      registry.fill(HIST("ld/fZvtx"), col.posZ());
    }
    tags2N(keepEvent2N[CFTrigger::kPD],
           keepEvent2N[CFTrigger::kLD]);

    if (!keepEvent3N[CFTrigger::kPPP] && !keepEvent3N[CFTrigger::kPPL] && !keepEvent3N[CFTrigger::kPLL] && !keepEvent3N[CFTrigger::kLLL] &&
        !keepEvent2N[CFTrigger::kPD] && !keepEvent2N[CFTrigger::kLD]) {
      registry.fill(HIST("fProcessedEvents"), 1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilter>(cfg)};
}
