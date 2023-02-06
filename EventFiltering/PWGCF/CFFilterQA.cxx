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
#include <cstdint>
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
#include "PWGCF/DataModel/FemtoDerived.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace CFTrigger
{
// enums
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
static const std::vector<std::string> nSpeciesName{"Proton", "Deuteron", "Lambda"};
static const std::vector<std::string> nSpeciesV0DaughterName{"Pion", "Proton"};
static const std::vector<std::string> nSpeciesRejectionName{"Proton", "Pion", "Electron"};
static const std::vector<std::string> nTPCCutName{"TPC min", "TPC max"};
static const std::vector<std::string> nPidCutsName{"TPC min", "TPC max", "TOF min", "TOF max", "TPCTOF max"};
static const std::vector<std::string> nPtCutsName{"Pt min", "Pt max", "P thres"};
static const std::vector<std::string> nThreeBodyFilterNames{"PPP", "PPL", "PLL", "LLL"};
static const std::vector<std::string> nTwoBodyFilterNames{"PD", "LD"};

static const int nPidRejection = 2;
static const int nPidCutsDaughers = 2;
static const int nPtCuts = 3;
static const int nAllTriggers = 6;

static const float pidcutsTable[kNParticleSpecies][kNPIDLimits]{
  {-6.f, 6.f, -6.f, 6.f, 6.f},
  {-6.f, 6.f, -99.f, 99.f, 99.f},
  {-6.f, 6.f, -99.f, 99.f, 99.f}};
static const float pidRejectionTable[kNParticleRejection][nPidRejection]{
  {-2.f, 2.f},
  {-2.f, 2.f},
  {-2.f, 2.f}};
static const float pidcutsV0DaughterTable[kNV0Daughters][nPidCutsDaughers]{
  {-6.f, 6.f},
  {-6.f, 6.f}};
static const float ptcutsTable[kNParticleSpecies][nPtCuts]{
  {0.35f, 6.f, 0.75f},
  {0.35f, 1.6f, 99.f},
  {0.35f, 6.f, 99.f}};

static const float triggerSwitches[1][nAllTriggers]{
  {1, 1, 1, 1, 1, 1}};

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

  Produces<aod::FemtoDreamCollisions> outputCollision;
  Produces<aod::FemtoDreamParticles> outputParts;

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

  Configurable<int> ConfCutBitPart{
    "ConfCutBitPart",
    8190,
    "Cutbit for particle (charge +1)"};
  Configurable<int> ConfCutBitAntiPart{
    "ConfCutBitAntiPart",
    8189,
    "Cutbit for antiparticle"};

  Configurable<int> ConfPidBitProton{
    "ConfPidBitProton",
    1,
    "Pidbit for proton"};
  // Configurable<int> ConfPidBitDeuteron{
  //   "ConfPidBitDeuteron",
  //   4,
  //   "Pidbit for proton"};

  // Configs for tracks
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
  Configurable<float> ConfTrkTPCnclsMin{
    "ConfTrkTPCnclsMin",
    65,
    "Minimum number of TPC clusters"};
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
  Configurable<float>
    ConfTPCAvg{
      "ConfTPCAvg",
      0.0f,
      "TPC center value, which is substracted in calculation of combined TPC and TOF nSigma"};
  Configurable<float>
    ConfTOFAvg{
      "ConfTOFAvg",
      0.0f,
      "TOF center value, which is substracted in calculation of combined TPC and TOF nSigma"};

  Configurable<LabeledArray<float>> ConfPIDCuts{
    "ConfPIDCuts",
    {CFTrigger::pidcutsTable[0], CFTrigger::kNParticleSpecies, CFTrigger::kNPIDLimits, CFTrigger::nSpeciesName, CFTrigger::nPidCutsName},
    "Particle PID selections"};
  Configurable<LabeledArray<float>> ConfPtCuts{
    "ConfPtCuts",
    {CFTrigger::ptcutsTable[0], CFTrigger::kNParticleSpecies, CFTrigger::nPtCuts, CFTrigger::nSpeciesName, CFTrigger::nPtCutsName},
    "Particle Momentum selections"};
  Configurable<bool> ConfRejectNOTDeuteron{
    "ConfRejectNOTDeuteron",
    false,
    "Reject deuteron candidates if they are compatible with electron, pion, proton"};
  Configurable<LabeledArray<float>> ConfPIDRejection{
    "ConfPIDRejection",
    {CFTrigger::pidRejectionTable[0], CFTrigger::kNParticleRejection, CFTrigger::nPidRejection, CFTrigger::nSpeciesRejectionName, CFTrigger::nTPCCutName},
    "Particle PID Rejection selections (Deuteron candidates only)"};

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
  Configurable<float> ConfV0DaughPIDnSigmaMax{
    "ConfV0DaughPIDnSigmaMax",
    6.f,
    "V0 Daugh sel: Max. PID nSigma TPC"};
  Configurable<LabeledArray<float>> ConfDaughPIDCuts{
    "ConfDaughPIDCuts",
    {CFTrigger::pidcutsV0DaughterTable[0], CFTrigger::kNV0Daughters, CFTrigger::nPidCutsDaughers, CFTrigger::nSpeciesV0DaughterName, CFTrigger::nTPCCutName},
    "Q3 limits for three body trigger"};

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

    // event cuts
    registry.add("EventCuts/fMultiplicityBefore", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fMultiplicityAfter", "Multiplicity after event cuts", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fZvtxBefore", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("EventCuts/fZvtxAfter", "Zvtx after event cuts", HistType::kTH1F, {{1000, -15, 15}});

    // all tracks
    registry.add("TrackCuts/fPtTrackBefore", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaTrackBefore", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiTrackBefore", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // PID vs momentum before cuts
    registry.add("TrackCuts/fNsigmaTPCvsPProtonBefore", "NSigmaTPC Proton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTOFvsPProtonBefore", "NSigmaTOF Proton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCTOFvsPProtonBefore", "NSigmaTPCTOF Proton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPAntiProtonBefore", "NSigmaTPC AntiProton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTOFvsPAntiProtonBefore", "NSigmaTOF AntiProton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCTOFvsPAntiProtonBefore", "NSigmaTPCTOF AntiProton Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    // registry.add("TrackCuts/fNsigmaTPCvsPDeuteronBefore", "NSigmaTPC Deuteron Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    // registry.add("TrackCuts/fNsigmaTPCvsPAntiDeuteronBefore", "NSigmaTPC AntiDeuteron Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // TPC signal
    registry.add("TrackCuts/fTPCSignal", "TPCSignal", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/fTPCSignalALLCUTS", "TPCSignalALLCUTS", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -1000.f, 1000.f}}});
    // TPC signal anti
    registry.add("TrackCuts/fTPCSignalAnti", "TPCSignal", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/fTPCSignalAntiALLCUTS", "TPCSignalALLCUTS", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});

    registry.add("TrackCuts/fTPCSignalProton", "fTPCSignalProton", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    registry.add("TrackCuts/fTPCSignalAntiProton", "fTPCSignalAntiProton", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    // registry.add("TrackCuts/fTPCSignalDeuteron", "fTPCSignalDeuteron", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});
    // registry.add("TrackCuts/fTPCSignalAntiDeuteron", "fTPCSignalAntiDeuteron", {HistType::kTH2F, {{1000, 0.0f, 6.0f}, {20000, -100.f, 1000.f}}});

    // PID vs momentum before cuts daughters
    registry.add("TrackCuts/fNsigmaTPCvsPProtonV0DaughBefore", "NSigmaTPC Proton V0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPPionMinusV0DaughBefore", "NSigmaTPC AntiPion V0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPAntiProtonAntiV0DaughBefore", "NSigmaTPC AntiProton antiV0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPPionPlusAntiV0DaughBefore", "NSigmaTPC Pion antiV0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // proton
    registry.add("TrackCuts/fPProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fPTPCProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fNsigmaTPCvsPTPCProton", "NSigmaTPC Proton P TPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    // registry.add("TrackCuts/fNsigmaTPCvsPTPCDeuteron", "NSigmaTPC Deuterons P TPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    registry.add("TrackCuts/fPtProton", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaProton", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiProton", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/fNsigmaTPCvsPProton", "NSigmaTPC Proton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTOFvsPProton", "NSigmaTOF Proton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCTOFvsPProton", "NSigmaTPCTOF Proton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    // antiproton
    registry.add("TrackCuts/fPtAntiProton", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaAntiProton", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiAntiProton", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/fNsigmaTPCvsPAntiProton", "NSigmaTPC AntiProton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTOFvsPAntiProton", "NSigmaTOF AntiProton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCTOFvsPAntiProton", "NSigmaTPCTOF AntiProton", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    // deuteron
    // registry.add("TrackCuts/fPtDeuteron", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    // registry.add("TrackCuts/fEtaDeuteron", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    // registry.add("TrackCuts/fPhiDeuteron", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    // registry.add("TrackCuts/fNsigmaTPCvsPDeuteron", "NSigmaTPC Deuteron", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // antideuteron
    // registry.add("TrackCuts/fPtAntiDeuteron", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    // registry.add("TrackCuts/fEtaAntiDeuteron", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    // registry.add("TrackCuts/fPhiAntiDeuteron", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    // registry.add("TrackCuts/fNsigmaTPCvsPAntiDeuteron", "NSigmaTPC AntiDeuteron", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // lambda before selections
    // registry.add("TrackCuts/fPtLambdaBefore", "Transverse momentum of all processed V0s before cuts", HistType::kTH1F, {{1000, 0, 10}});
    // registry.add("TrackCuts/fInvMassLambdaBefore", "Invariant mass of all processed V0s (Lambda) before cuts", HistType::kTH1F, {{1000, 0.7, 1.5}});

    // lambda
    // registry.add("TrackCuts/fPtLambda", "Transverse momentum of all selected V0s", HistType::kTH1F, {{1000, 0, 10}});
    // registry.add("TrackCuts/fInvMassLambda", "Invariant mass of all selected V0s (Lambda)", HistType::kTH1F, {{1000, 0.7, 1.5}});

    // antilambda
    // registry.add("TrackCuts/fPtAntiLambda", "Transverse momentum of all selected V0s", HistType::kTH1F, {{1000, 0, 10}});
    // registry.add("TrackCuts/fInvMassAntiLambda", "Invariant mass of all selected V0s (Lambda)", HistType::kTH1F, {{1000, 0.7, 1.5}});
  }

  float mMassProton = o2::constants::physics::MassProton;
  float mMassLambda = o2::constants::physics::MassLambda;
  float mMassDeuteron = o2::constants::physics::MassDeuteron;
  int currentRunNumber = -999;
  int lastRunNumber = -999;

  template <typename T>
  int getRowDaughters(int daughID, T const& vecID)
  {
    int rowInPrimaryTrackTableDaugh = -1;
    for (size_t i = 0; i < vecID.size(); i++) {
      if (vecID.at(i) == daughID) {
        rowInPrimaryTrackTableDaugh = i;
        break;
      }
    }
    return rowInPrimaryTrackTableDaugh;
  }

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
    if (tpcNClsF < ConfTrkTPCnclsMin) {
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
  bool isSelectedTrackPID(T const& track, CFTrigger::ParticleSpecies partSpecies, bool Rejection, double nSigmaTPC[2])
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
    // compute nsigma
    switch (partSpecies) {
      case CFTrigger::kProton:
        nSigma = (pThres) ? nSigmaTPC[0]
                          : std::sqrt(std::pow(nSigmaTPC[0] - ConfTPCAvg, 2) +
                                      std::pow(track.tofNSigmaPr() - ConfTOFAvg, 2));
        break;
      case CFTrigger::kDeuteron:
        nSigma = (pThres) ? nSigmaTPC[1]
                          : std::sqrt(std::pow(nSigmaTPC[1] - ConfTPCAvg, 2) +
                                      std::pow(track.tofNSigmaDe() - ConfTOFAvg, 2));
        break;
      case CFTrigger::kLambda:
        LOG(fatal) << "No PID selection for Lambdas";
        break;
      default:
        LOG(fatal) << "Particle species not known";
    }
    // check if track is selected
    if (pThres) {
      if (nSigma > ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMin) &&
          nSigma < ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMax)) {
        isSelected = true;
      }
    } else {
      if (nSigma < ConfPIDCuts->get(partSpecies, CFTrigger::kTPCTOF)) {
        isSelected = true;
      }
    }
    // for deuterons normally, we want to reject tracks that have a high
    // probablilty of being another particle
    if (Rejection) {
      double nSigmaPi{track.tpcNSigmaPi()};
      double nSigmaEl{track.tpcNSigmaEl()};
      if (ConfUseManualPIDpion) {
        // TO BE IMPLEMENTED
      }
      if (ConfUseManualPIDel) {
        // TO BE IMPLEMENTED
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

  // void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, o2::aod::V0Datas const& fullV0s)
  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {

    if (!ConfIsRun3) {
      LOG(fatal) << "Run 2 processing is not implemented!";
    }
    if (ConfUseManualPIDproton || ConfUseManualPIDdeuteron) {
      currentRunNumber = col.bc_as<aod::BCsWithTimestamps>().runNumber();
      if (currentRunNumber != lastRunNumber) {
        auto bc = col.bc_as<aod::BCsWithTimestamps>();
        if (ConfUseManualPIDproton) {
          BBProton = setValuesBB(bc, ConfPIDBBProton);
          BBAntiproton = setValuesBB(bc, ConfPIDBBAntiProton);
        }
        if (ConfUseManualPIDdeuteron) {
          BBDeuteron = setValuesBB(bc, ConfPIDBBDeuteron);
          BBAntideuteron = setValuesBB(bc, ConfPIDBBAntiDeuteron);
        }
        if (ConfUseManualPIDpion) {
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

    registry.fill(HIST("EventCuts/fMultiplicityBefore"), col.multNTracksPV());
    registry.fill(HIST("EventCuts/fZvtxBefore"), col.posZ());

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    if (isSelectedEvent(col)) {

      outputCollision(col.posZ(), col.multFV0M(), col.multNTracksPV(), -2, -2);

      registry.fill(HIST("EventCuts/fMultiplicityAfter"), col.multNTracksPV());
      registry.fill(HIST("EventCuts/fZvtxAfter"), col.posZ());

      // keep track of proton indices
      std::vector<int> ProtonIndex = {};
      std::vector<int> AntiProtonIndex = {};

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

        registry.fill(HIST("TrackCuts/fPtTrackBefore"), track.pt());
        registry.fill(HIST("TrackCuts/fEtaTrackBefore"), track.eta());
        registry.fill(HIST("TrackCuts/fPhiTrackBefore"), track.phi());

        if (track.sign() > 0) {
          // Fill PID info
          registry.fill(HIST("TrackCuts/fTPCSignal"), track.tpcInnerParam(), track.tpcSignal());
          if (isSelectedTrack(track, CFTrigger::kProton)) {
            registry.fill(HIST("TrackCuts/fTPCSignalALLCUTS"), track.tpcInnerParam(), track.tpcSignal());
          }
          registry.fill(HIST("TrackCuts/fNsigmaTPCvsPProtonBefore"), track.tpcInnerParam(), nTPCSigmaP[0]);
          registry.fill(HIST("TrackCuts/fNsigmaTOFvsPProtonBefore"), track.tpcInnerParam(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/fNsigmaTPCTOFvsPProtonBefore"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaP[0] - ConfTPCAvg, 2) + std::pow(track.tofNSigmaPr() - ConfTOFAvg, 2)));
          // registry.fill(HIST("TrackCuts/fNsigmaTPCvsPDeuteronBefore"), track.tpcInnerParam(), nTPCSigmaP[1]);
        }
        if (track.sign() < 0) {

          registry.fill(HIST("TrackCuts/fTPCSignalAnti"), track.tpcInnerParam(), track.tpcSignal());
          if (isSelectedTrack(track, CFTrigger::kProton)) {
            registry.fill(HIST("TrackCuts/fTPCSignalAntiALLCUTS"), track.tpcInnerParam(), track.tpcSignal());
          }

          registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiProtonBefore"), track.tpcInnerParam(), nTPCSigmaN[0]);
          registry.fill(HIST("TrackCuts/fNsigmaTOFvsPAntiProtonBefore"), track.tpcInnerParam(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/fNsigmaTOFvsPAntiProtonBefore"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaN[0] - ConfTPCAvg, 2) + std::pow(track.tofNSigmaPr() - ConfTOFAvg, 2)));
          // registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiDeuteronBefore"), track.tpcInnerParam(), nTPCSigmaN[1]);
        }

        // get protons
        if (isSelectedTrack(track, CFTrigger::kProton)) {
          if (track.sign() > 0 && isSelectedTrackPID(track, CFTrigger::kProton, false, nTPCSigmaP)) {
            outputParts(outputCollision.lastIndex(),
                        track.pt(),
                        track.eta(),
                        track.phi(),
                        aod::femtodreamparticle::ParticleType::kTrack,
                        static_cast<uint32_t>(ConfCutBitPart.value), // cutbit for particle
                        static_cast<uint32_t>(ConfPidBitProton.value),
                        track.dcaXY(),
                        childIDs,
                        0.f,
                        0.f);
            tmpIDtrack.push_back(track.globalIndex());
            registry.fill(HIST("TrackCuts/fTPCSignalProton"), track.tpcInnerParam(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/fPProton"), track.p());
            registry.fill(HIST("TrackCuts/fPTPCProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/fPtProton"), track.pt());
            registry.fill(HIST("TrackCuts/fEtaProton"), track.eta());
            registry.fill(HIST("TrackCuts/fPhiProton"), track.phi());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPProton"), track.tpcInnerParam(), nTPCSigmaP[0]);
            registry.fill(HIST("TrackCuts/fNsigmaTOFvsPProton"), track.tpcInnerParam(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/fNsigmaTPCTOFvsPProton"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaP[0] - ConfTPCAvg, 2) + std::pow(track.tofNSigmaPr() - ConfTOFAvg, 2)));
          }
          if (track.sign() < 0 && isSelectedTrackPID(track, CFTrigger::kProton, false, nTPCSigmaN)) {
            outputParts(outputCollision.lastIndex(),
                        track.pt(),
                        track.eta(),
                        track.phi(),
                        aod::femtodreamparticle::ParticleType::kTrack,
                        static_cast<uint32_t>(ConfCutBitAntiPart.value), // cutbit for antiparticle
                        static_cast<uint32_t>(ConfPidBitProton.value),
                        track.dcaXY(),
                        childIDs,
                        0.f,
                        0.f);
            tmpIDtrack.push_back(track.globalIndex());
            registry.fill(HIST("TrackCuts/fTPCSignalAntiProton"), track.tpcInnerParam(), track.tpcSignal());
            registry.fill(HIST("TrackCuts/fPtAntiProton"), track.pt());
            registry.fill(HIST("TrackCuts/fEtaAntiProton"), track.eta());
            registry.fill(HIST("TrackCuts/fPhiAntiProton"), track.phi());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiProton"), track.tpcInnerParam(), nTPCSigmaN[0]);
            registry.fill(HIST("TrackCuts/fNsigmaTOFvsPAntiProton"), track.tpcInnerParam(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/fNsigmaTPCTOFvsPAntiProton"), track.tpcInnerParam(), std::sqrt(std::pow(nTPCSigmaN[0] - ConfTPCAvg, 2) + std::pow(track.tofNSigmaPr() - ConfTOFAvg, 2)));
          }
        }
        // get deuterons
        // if (isSelectedTrack(track, CFTrigger::kDeuteron)) {
        //   if (track.sign() > 0 && isSelectedTrackPID(track, CFTrigger::kDeuteron, ConfRejectNOTDeuteron.value, nTPCSigmaP)) {
        //     outputParts(outputCollision.lastIndex(),
        //                 track.pt(),
        //                 track.eta(),
        //                 track.phi(),
        //                 aod::femtodreamparticle::ParticleType::kTrack,
        //                 static_cast<uint32_t>(ConfCutBitPart.value), // cutbit for particle
        //                 static_cast<uint32_t>(ConfPidBitDeuteron.value),
        //                 track.dcaXY(),
        //                 childIDs,
        //                 0.f,
        //                 0.f);
        //     registry.fill(HIST("TrackCuts/fTPCSignalDeuteron"), track.tpcInnerParam(), track.tpcSignal());
        //     registry.fill(HIST("TrackCuts/fPtDeuteron"), track.pt());
        //     registry.fill(HIST("TrackCuts/fEtaDeuteron"), track.eta());
        //     registry.fill(HIST("TrackCuts/fPhiDeuteron"), track.phi());
        //     registry.fill(HIST("TrackCuts/fNsigmaTPCvsPDeuteron"), track.tpcInnerParam(), nTPCSigmaP[1]);
        //   }
        //   if (track.sign() < 0 && isSelectedTrackPID(track, CFTrigger::kDeuteron, ConfRejectNOTDeuteron.value, nTPCSigmaN)) {
        //     outputParts(outputCollision.lastIndex(),
        //                 track.pt(),
        //                 track.eta(),
        //                 track.phi(),
        //                 aod::femtodreamparticle::ParticleType::kTrack,
        //                 static_cast<uint32_t>(ConfCutBitAntiPart.value), // cutbit for antiparticle
        //                 static_cast<uint32_t>(ConfPidBitDeuteron.value),
        //                 track.dcaXY(),
        //                 childIDs,
        //                 0.f,
        //                 0.f);
        //     registry.fill(HIST("TrackCuts/fTPCSignalAntiDeuteron"), track.tpcInnerParam(), track.tpcSignal());
        //     registry.fill(HIST("TrackCuts/fPtAntiDeuteron"), track.pt());
        //     registry.fill(HIST("TrackCuts/fEtaAntiDeuteron"), track.eta());
        //     registry.fill(HIST("TrackCuts/fPhiAntiDeuteron"), track.phi());
        //     registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiDeuteron"), track.tpcInnerParam(), nTPCSigmaN[1]);
        //   }
        // }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilter>(cfg)};
}
