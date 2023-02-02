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

static const float Q3Limits[1][kNThreeBodyTriggers]{
  {0.6f, 0.6f, 0.6f, 0.6f}};
static const float KstarLimits[1][kNTwoBodyTriggers]{
  {1.2f, 1.2f}};
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

  // Trigger selections
  Configurable<LabeledArray<float>> ConfTriggerSwitches{
    "ConfTriggerSwitches",
    {CFTrigger::triggerSwitches[0], 1, CFTrigger::nAllTriggers, std::vector<std::string>{"Switch"}, CFTrigger::CFTriggerNamesALL},
    "Turn on specific trigger"};

  Configurable<LabeledArray<float>> ConfQ3Limits{
    "ConfQ3Limits",
    {CFTrigger::Q3Limits[0], 1, CFTrigger::kNThreeBodyTriggers, std::vector<std::string>{"Limit"}, CFTrigger::nThreeBodyFilterNames},
    "Q3 limits for three body trigger"};

  Configurable<LabeledArray<float>> ConfKstarLimits{
    "ConfKstarLimits",
    {CFTrigger::KstarLimits[0], 1, CFTrigger::kNTwoBodyTriggers, std::vector<std::string>{"Limit"}, CFTrigger::nTwoBodyFilterNames},
    "kstar limit for two body trigger"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  void init(o2::framework::InitContext&)
  {
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
    registry.add("TrackCuts/fNsigmaTPCvsPDeuteronBefore", "NSigmaTPC Deuteron Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPAntiDeuteronBefore", "NSigmaTPC AntiDeuteron Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // PID vs momentum before cuts daughters
    registry.add("TrackCuts/fNsigmaTPCvsPProtonV0DaughBefore", "NSigmaTPC Proton V0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPPionMinusV0DaughBefore", "NSigmaTPC AntiPion V0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPAntiProtonAntiV0DaughBefore", "NSigmaTPC AntiProton antiV0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/fNsigmaTPCvsPPionPlusAntiV0DaughBefore", "NSigmaTPC Pion antiV0Daught Before", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // proton
    // TEST P TPC
    registry.add("TrackCuts/fPProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fPTPCProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fNsigmaTPCvsPTPCProton", "NSigmaTPC Proton P TPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

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
    registry.add("TrackCuts/fPtDeuteron", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaDeuteron", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiDeuteron", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/fNsigmaTPCvsPDeuteron", "NSigmaTPC Deuteron", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // antideuteron
    registry.add("TrackCuts/fPtAntiDeuteron", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaAntiDeuteron", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiAntiDeuteron", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    registry.add("TrackCuts/fNsigmaTPCvsPAntiDeuteron", "NSigmaTPC AntiDeuteron", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    // lambda before selections
    registry.add("TrackCuts/fPtLambdaBefore", "Transverse momentum of all processed V0s before cuts", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fInvMassLambdaBefore", "Invariant mass of all processed V0s (Lambda) before cuts", HistType::kTH1F, {{1000, 0.7, 1.5}});

    // lambda
    registry.add("TrackCuts/fPtLambda", "Transverse momentum of all selected V0s", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fInvMassLambda", "Invariant mass of all selected V0s (Lambda)", HistType::kTH1F, {{1000, 0.7, 1.5}});

    // antilambda
    registry.add("TrackCuts/fPtAntiLambda", "Transverse momentum of all selected V0s", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fInvMassAntiLambda", "Invariant mass of all selected V0s (Lambda)", HistType::kTH1F, {{1000, 0.7, 1.5}});

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

  float mMassProton = o2::constants::physics::MassProton;
  float mMassLambda = o2::constants::physics::MassLambda;
  float mMassDeuteron = o2::constants::physics::MassDeuteron;

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
  bool isSelectedV0Daughter(T const& track, float charge, CFTrigger::V0Daughters species)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
    const auto sign = track.sign();
    float nSigmaTPC = -999.f;

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
        nSigmaTPC = track.tpcNSigmaPi();
        break;
      case CFTrigger::kDaughProton:
        nSigmaTPC = track.tpcNSigmaPr();
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
  bool isSelectedTrackPID(T const& track, CFTrigger::ParticleSpecies partSpecies, bool Rejection = false)
  {
    bool isSelected = false;
    bool pThres = true;
    float nSigma = -999.;
    float TPCAvg = 0; //(ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMax) + ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMin)) / 2.;
    float TOFAvg = 0; //(ConfPIDCuts->get(partSpecies, CFTrigger::kTOFMax) + ConfPIDCuts->get(partSpecies, CFTrigger::kTOFMin)) / 2.;

    // check momentum threshold
    if (track.p() <= ConfPtCuts->get(partSpecies, "P thres")) {
      pThres = true;
    } else {
      pThres = false;
    }
    // compute nsigma
    switch (partSpecies) {
      case CFTrigger::kProton:
        nSigma = (pThres) ? track.tpcNSigmaPr()
                          : std::sqrt(std::pow(track.tpcNSigmaPr() - TPCAvg, 2) +
                                      std::pow(track.tofNSigmaPr() - TOFAvg, 2));
        break;
      case CFTrigger::kDeuteron:
        nSigma = (pThres) ? track.tpcNSigmaDe()
                          : std::sqrt(std::pow(track.tpcNSigmaDe() - TPCAvg, 2) +
                                      std::pow(track.tofNSigmaDe() - TOFAvg, 2));
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
      if ((ConfPIDRejection->get(CFTrigger::kRejProton, CFTrigger::kTPCMin) < track.tpcNSigmaPr() &&
           ConfPIDRejection->get(CFTrigger::kRejProton, CFTrigger::kTPCMax) > track.tpcNSigmaPr()) ||
          (ConfPIDRejection->get(CFTrigger::kRejPion, CFTrigger::kTPCMin) < track.tpcNSigmaPi() &&
           ConfPIDRejection->get(CFTrigger::kRejPion, CFTrigger::kTPCMax) > track.tpcNSigmaPi()) ||
          (ConfPIDRejection->get(CFTrigger::kRejElectron, CFTrigger::kTPCMin) < track.tpcNSigmaEl() &&
           ConfPIDRejection->get(CFTrigger::kRejElectron, CFTrigger::kTPCMax) > track.tpcNSigmaEl())) {
        return false;
      }
    }
    return isSelected;
  }

  template <typename C, typename V, typename T>
  bool isSelectedMinimalV0(C const& col, V const& v0, T const& posTrack,
                           T const& negTrack, float charge)
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
      if (!isSelectedV0Daughter(posTrack, 1, CFTrigger::kDaughProton)) {
        return false;
      }
      if (!isSelectedV0Daughter(negTrack, -1, CFTrigger::kDaughPion)) {
        return false;
      }
    }
    if (charge < 0) {
      if (!isSelectedV0Daughter(posTrack, 1, CFTrigger::kDaughPion)) {
        return false;
      }
      if (!isSelectedV0Daughter(negTrack, -1, CFTrigger::kDaughProton)) {
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

  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks, o2::aod::V0Datas const& fullV0s)
  {

    if (!ConfIsRun3) {
      LOG(fatal) << "Run 2 processing is not implemented!";
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

        registry.fill(HIST("TrackCuts/fPtTrackBefore"), track.pt());
        registry.fill(HIST("TrackCuts/fEtaTrackBefore"), track.eta());
        registry.fill(HIST("TrackCuts/fPhiTrackBefore"), track.phi());

        if (track.sign() > 0) {
          // Fill PID info
          registry.fill(HIST("TrackCuts/fNsigmaTPCvsPProtonBefore"), track.p(), track.tpcNSigmaPr());
          registry.fill(HIST("TrackCuts/fNsigmaTOFvsPProtonBefore"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/fNsigmaTPCTOFvsPProtonBefore"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));
          registry.fill(HIST("TrackCuts/fNsigmaTPCvsPDeuteronBefore"), track.p(), track.tpcNSigmaDe());
        }
        if (track.sign() < 0) {
          registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiProtonBefore"), track.p(), track.tpcNSigmaPr());
          registry.fill(HIST("TrackCuts/fNsigmaTOFvsPAntiProtonBefore"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/fNsigmaTOFvsPAntiProtonBefore"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));
          registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiDeuteronBefore"), track.p(), track.tpcNSigmaDe());
        }

        // get protons
        if (isSelectedTrack(track, CFTrigger::kProton) && isSelectedTrackPID(track, CFTrigger::kProton, false)) {
          ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), mMassProton);
          if (track.sign() > 0) {
            protons.push_back(temp);
            ProtonIndex.push_back(track.globalIndex());
            registry.fill(HIST("TrackCuts/fPProton"), track.p());
            registry.fill(HIST("TrackCuts/fPTPCProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPTPCProton"), track.tpcInnerParam(), track.tpcNSigmaPr());

            registry.fill(HIST("TrackCuts/fPtProton"), track.pt());
            registry.fill(HIST("TrackCuts/fEtaProton"), track.eta());
            registry.fill(HIST("TrackCuts/fPhiProton"), track.phi());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPProton"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/fNsigmaTOFvsPProton"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/fNsigmaTPCTOFvsPProton"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));
          }
          if (track.sign() < 0) {
            antiprotons.push_back(temp);
            AntiProtonIndex.push_back(track.globalIndex());
            registry.fill(HIST("TrackCuts/fPtAntiProton"), track.pt());
            registry.fill(HIST("TrackCuts/fEtaAntiProton"), track.eta());
            registry.fill(HIST("TrackCuts/fPhiAntiProton"), track.phi());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiProton"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/fNsigmaTOFvsPAntiProton"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/fNsigmaTPCTOFvsPAntiProton"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));
          }
        }
        // get deuterons
        if (isSelectedTrack(track, CFTrigger::kDeuteron) && isSelectedTrackPID(track, CFTrigger::kDeuteron, ConfRejectNOTDeuteron.value)) {
          ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), mMassDeuteron);
          if (track.sign() > 0) {
            deuterons.push_back(temp);
            registry.fill(HIST("TrackCuts/fPtDeuteron"), track.pt());
            registry.fill(HIST("TrackCuts/fEtaDeuteron"), track.eta());
            registry.fill(HIST("TrackCuts/fPhiDeuteron"), track.phi());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPDeuteron"), track.p(), track.tpcNSigmaDe());
          }
          if (track.sign() < 0) {
            antideuterons.push_back(temp);
            registry.fill(HIST("TrackCuts/fPtAntiDeuteron"), track.pt());
            registry.fill(HIST("TrackCuts/fEtaAntiDeuteron"), track.eta());
            registry.fill(HIST("TrackCuts/fPhiAntiDeuteron"), track.phi());
            registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiDeuteron"), track.p(), track.tpcNSigmaDe());
          }
        }
      }

      // keep track of daugher indices to avoid selfcorrelations
      std::vector<int> LambdaPosDaughIndex = {};
      std::vector<int> LambdaNegDaughIndex = {};
      std::vector<int> AntiLambdaPosDaughIndex = {};
      std::vector<int> AntiLambdaNegDaughIndex = {};

      for (auto& v0 : fullV0s) {

        registry.fill(HIST("TrackCuts/fPtLambdaBefore"), v0.pt());
        registry.fill(HIST("TrackCuts/fInvMassLambdaBefore"), v0.mLambda());

        auto postrack = v0.template posTrack_as<aod::FemtoFullTracks>();
        auto negtrack = v0.template negTrack_as<aod::FemtoFullTracks>();

        registry.fill(HIST("TrackCuts/fNsigmaTPCvsPProtonV0DaughBefore"), postrack.p(), postrack.tpcNSigmaPr());
        registry.fill(HIST("TrackCuts/fNsigmaTPCvsPPionPlusAntiV0DaughBefore"), postrack.p(), postrack.tpcNSigmaPi());

        registry.fill(HIST("TrackCuts/fNsigmaTPCvsPPionMinusV0DaughBefore"), negtrack.p(), negtrack.tpcNSigmaPi());
        registry.fill(HIST("TrackCuts/fNsigmaTPCvsPAntiProtonAntiV0DaughBefore"), negtrack.p(), negtrack.tpcNSigmaPr());

        if (isSelectedMinimalV0(col, v0, postrack, negtrack, 1)) {
          ROOT::Math::PtEtaPhiMVector temp(v0.pt(), v0.eta(), v0.phi(), mMassLambda);
          lambdas.push_back(temp);
          LambdaPosDaughIndex.push_back(postrack.globalIndex());
          LambdaNegDaughIndex.push_back(negtrack.globalIndex());
          registry.fill(HIST("TrackCuts/fPtLambda"), v0.pt());
          registry.fill(HIST("TrackCuts/fInvMassLambda"), v0.mLambda());
        }
        if (isSelectedMinimalV0(col, v0, postrack, negtrack, -1)) {
          ROOT::Math::PtEtaPhiMVector temp(v0.pt(), v0.eta(), v0.phi(), mMassLambda);
          antilambdas.push_back(temp);
          AntiLambdaPosDaughIndex.push_back(postrack.globalIndex());
          AntiLambdaNegDaughIndex.push_back(negtrack.globalIndex());
          registry.fill(HIST("TrackCuts/fPtAntiLambda"), v0.pt());
          registry.fill(HIST("TrackCuts/fInvMassAntiLambda"), v0.mAntiLambda());
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
    //     for (int i = 0; i < CFTrigger::kNThreeBodyTriggers; i++) {
    //       if (keepEvent3N[i]) {
    //         registry.fill(HIST("fProcessedEvents"), i + 2);
    //       }
    //     }
    //   }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilter>(cfg)};
}
