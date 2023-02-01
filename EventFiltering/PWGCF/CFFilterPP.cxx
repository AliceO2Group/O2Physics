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

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <vector>

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using namespace o2::analysis::femtoDream;

namespace CFTrigger
{
// enums
enum PIDLimits { kTPCMin,
                 kTPCMax,
                 kTOFMin,
                 kTOFMax,
                 kTPCTOF,
                 kNPIDLimits
};

// For configurable tables
static const std::vector<std::string> nTPCCutName{"TPC min", "TPC max"};
static const std::vector<std::string> nPidCutsName{"TPC min", "TPC max", "TOF min", "TOF max", "TPCTOF max"};
static const std::vector<std::string> nPtCutsName{"Pt min", "Pt max", "P thres"};

static const int nPtCuts = 3;

static const float pidcutsTable[1][kNPIDLimits]{
  {-6.f, 6.f, -6.f, 6.f, 6.f}};
static const float ptcutsTable[1][nPtCuts]{
  {0.35f, 6.f, 0.75f}};
} // namespace CFTrigger

namespace o2::aod
{
using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;

using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTOFFullPr>;
} // namespace o2::aod

struct CFFilter {

  Produces<aod::FemtoDreamCollisions> outputCollision;
  Produces<aod::FemtoDreamParticles> outputParts;

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
    {CFTrigger::pidcutsTable[0], 1, CFTrigger::kNPIDLimits, std::vector<std::string>{"Proton"}, CFTrigger::nPidCutsName},
    "Particle PID selections"};
  Configurable<LabeledArray<float>> ConfPtCuts{
    "ConfPtCuts",
    {CFTrigger::ptcutsTable[0], 1, CFTrigger::nPtCuts, std::vector<std::string>{""}, CFTrigger::nPtCutsName},
    "Particle Momentum selections"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  void init(o2::framework::InitContext&)
  {
    // event cuts
    registry.add("EventCuts/fMultiplicityBefore", "Multiplicity of all processed events", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fMultiplicityAfter", "Multiplicity after event cuts", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("EventCuts/fZvtxBefore", "Zvtx of all processed events", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("EventCuts/fZvtxAfter", "Zvtx after event cuts", HistType::kTH1F, {{1000, -15, 15}});

    // all tracks
    registry.add("TrackCuts/fPtTrackBefore", "Transverse momentum of all processed tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/fEtaTrackBefore", "Pseudorapidity of all processed tracks", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/fPhiTrackBefore", "Azimuthal angle of all processed tracks", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // proton
    // before cuts
    registry.add("TrackCuts/Proton/Before/fPProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/Before/fPTPCProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/Before/fPtProton", "Transverse momentum", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/Before/fEtaProton", "Pseudorapidity", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Proton/Before/fPhiProton", "Azimuthal angle", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    registry.add("TrackCuts/Proton/Before/fNsigmaTPCvsPProton", "NSigmaTPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/Before/fNsigmaTPCvsPTPCProton", "NSigmaTPC p_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/Before/fNsigmaTOFvsPProton", "NSigmaTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/Before/fNsigmaTPCTOFvsPProton", "NSigmaTPCTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    // after cuts
    registry.add("TrackCuts/Proton/After/fPProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/After/fPTPCProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/After/fPtProton", "Transverse momentum", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/Proton/After/fEtaProton", "Pseudorapidity", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/Proton/After/fPhiProton", "Azimuthal angle", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    registry.add("TrackCuts/Proton/After/fNsigmaTPCvsPProton", "NSigmaTPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/After/fNsigmaTPCvsPTPCProton", "NSigmaTPC p_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/After/fNsigmaTOFvsPProton", "NSigmaTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/Proton/After/fNsigmaTPCTOFvsPProton", "NSigmaTPCTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    // antiproton
    // before cuts
    registry.add("TrackCuts/AntiProton/Before/fPAntiProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/Before/fPTPCAntiProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/Before/fPtAntiProton", "Transverse momentum", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/Before/fEtaAntiProton", "Pseudorapidity", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiProton/Before/fPhiAntiProton", "Azimuthal angle", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    registry.add("TrackCuts/AntiProton/Before/fNsigmaTPCvsPAntiProton", "NSigmaTPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/Before/fNsigmaTPCvsPTPCAntiProton", "NSigmaTPC p_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/Before/fNsigmaTOFvsPAntiProton", "NSigmaTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/Before/fNsigmaTPCTOFvsPAntiProton", "NSigmaTPCTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    // after cuts
    registry.add("TrackCuts/AntiProton/After/fPAntiProton", "Momentum of protons at PV", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/After/fPTPCAntiProton", "Momentum of protons at TPC inner wall", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/After/fPtAntiProton", "Transverse momentum", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("TrackCuts/AntiProton/After/fEtaAntiProton", "Pseudorapidity", HistType::kTH1F, {{1000, -2, 2}});
    registry.add("TrackCuts/AntiProton/After/fPhiAntiProton", "Azimuthal angle", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    registry.add("TrackCuts/AntiProton/After/fNsigmaTPCvsPAntiProton", "NSigmaTPC", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/After/fNsigmaTPCvsPTPCAntiProton", "NSigmaTPC p_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/After/fNsigmaTOFvsPAntiProton", "NSigmaTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    registry.add("TrackCuts/AntiProton/After/fNsigmaTPCTOFvsPAntiProton", "NSigmaTPCTOF", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
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
  bool isSelectedTrack(T const& track)
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

    if (pT < ConfPtCuts->get(0u, "Pt min")) {
      return false;
    }
    if (pT > ConfPtCuts->get(0u, "Pt max")) {
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
  bool isSelectedTrackPID(T const& track)
  {
    bool isSelected = false;
    bool pThres = true;
    float nSigma = -999.;
    float TPCAvg = (ConfPIDCuts->get(0u, CFTrigger::kTPCMax) + ConfPIDCuts->get(0u, CFTrigger::kTPCMin)) / 2.;
    float TOFAvg = (ConfPIDCuts->get(0u, CFTrigger::kTOFMax) + ConfPIDCuts->get(0u, CFTrigger::kTOFMin)) / 2.;

    // check momentum threshold
    if (track.p() <= ConfPtCuts->get(0u, "P thres")) {
      pThres = true;
    } else {
      pThres = false;
    }
    // compute nsigma
    nSigma = (pThres) ? track.tpcNSigmaPr()
                      : std::sqrt(std::pow(track.tpcNSigmaPr() - TPCAvg, 2) +
                                  std::pow(track.tofNSigmaPr() - TOFAvg, 2));
    // check if track is selected
    if (pThres) {
      if (nSigma > ConfPIDCuts->get(0u, CFTrigger::kTPCMin) &&
          nSigma < ConfPIDCuts->get(0u, CFTrigger::kTPCMax)) {
        isSelected = true;
      }
    } else {
      if (nSigma < ConfPIDCuts->get(0u, CFTrigger::kTPCTOF)) {
        isSelected = true;
      }
    }
    return isSelected;
  }

  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {

    if (!ConfIsRun3) {
      LOG(fatal) << "Run 2 processing is not implemented!";
    }

    registry.fill(HIST("EventCuts/fMultiplicityBefore"), col.multNTracksPV());
    registry.fill(HIST("EventCuts/fZvtxBefore"), col.posZ());

    int childIDs[2] = {0, 0};
    uint32_t pidCuts = 1; // set to 0

    if (isSelectedEvent(col)) {

      registry.fill(HIST("EventCuts/fMultiplicityAfter"), col.multNTracksPV());
      registry.fill(HIST("EventCuts/fZvtxAfter"), col.posZ());

      outputCollision(col.posZ(), col.multFV0M(), col.multNTracksPV(), -2, -2);

      for (auto& track : tracks) {

        registry.fill(HIST("TrackCuts/fPtTrackBefore"), track.pt());
        registry.fill(HIST("TrackCuts/fEtaTrackBefore"), track.eta());
        registry.fill(HIST("TrackCuts/fPhiTrackBefore"), track.phi());

        if (track.sign() > 0) {
          registry.fill(HIST("TrackCuts/Proton/Before/fPProton"), track.p());
          registry.fill(HIST("TrackCuts/Proton/Before/fPTPCProton"), track.tpcInnerParam());
          registry.fill(HIST("TrackCuts/Proton/Before/fPtProton"), track.pt());
          registry.fill(HIST("TrackCuts/Proton/Before/fEtaProton"), track.eta());
          registry.fill(HIST("TrackCuts/Proton/Before/fPhiProton"), track.phi());
          registry.fill(HIST("TrackCuts/Proton/Before/fNsigmaTPCvsPProton"), track.p(), track.tpcNSigmaPr());
          registry.fill(HIST("TrackCuts/Proton/Before/fNsigmaTPCvsPTPCProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
          registry.fill(HIST("TrackCuts/Proton/Before/fNsigmaTOFvsPProton"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/Proton/Before/fNsigmaTPCTOFvsPProton"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));
        }

        if (track.sign() < 0) {
          registry.fill(HIST("TrackCuts/AntiProton/Before/fPAntiProton"), track.p());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fPTPCAntiProton"), track.tpcInnerParam());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fPtAntiProton"), track.pt());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fEtaAntiProton"), track.eta());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fPhiAntiProton"), track.phi());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fNsigmaTPCvsPAntiProton"), track.p(), track.tpcNSigmaPr());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fNsigmaTPCvsPTPCAntiProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fNsigmaTOFvsPAntiProton"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("TrackCuts/AntiProton/Before/fNsigmaTPCTOFvsPAntiProton"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));
        }

        // get protons
        if (isSelectedTrack(track) && isSelectedTrackPID(track)) {
          if (track.sign() > 0) {
            registry.fill(HIST("TrackCuts/Proton/After/fPProton"), track.p());
            registry.fill(HIST("TrackCuts/Proton/After/fPTPCProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/Proton/After/fPtProton"), track.pt());
            registry.fill(HIST("TrackCuts/Proton/After/fEtaProton"), track.eta());
            registry.fill(HIST("TrackCuts/Proton/After/fPhiProton"), track.phi());
            registry.fill(HIST("TrackCuts/Proton/After/fNsigmaTPCvsPProton"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/Proton/After/fNsigmaTPCvsPTPCProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/Proton/After/fNsigmaTOFvsPProton"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/Proton/After/fNsigmaTPCTOFvsPProton"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));

            outputParts(outputCollision.lastIndex(),
                        track.pt(),
                        track.eta(),
                        track.phi(),
                        0u, // aod::femtodreamparticle::ParticleType::kTrack,
                        0u, // charge 0 for proton
                        pidCuts,
                        track.dcaXY(),
                        childIDs,
                        0.f,
                        0.f);
          }
          if (track.sign() < 0) {
            registry.fill(HIST("TrackCuts/AntiProton/After/fPAntiProton"), track.p());
            registry.fill(HIST("TrackCuts/AntiProton/After/fPTPCAntiProton"), track.tpcInnerParam());
            registry.fill(HIST("TrackCuts/AntiProton/After/fPtAntiProton"), track.pt());
            registry.fill(HIST("TrackCuts/AntiProton/After/fEtaAntiProton"), track.eta());
            registry.fill(HIST("TrackCuts/AntiProton/After/fPhiAntiProton"), track.phi());
            registry.fill(HIST("TrackCuts/AntiProton/After/fNsigmaTPCvsPAntiProton"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/AntiProton/After/fNsigmaTPCvsPTPCAntiProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
            registry.fill(HIST("TrackCuts/AntiProton/After/fNsigmaTOFvsPAntiProton"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("TrackCuts/AntiProton/After/fNsigmaTPCTOFvsPAntiProton"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)));

            outputParts(outputCollision.lastIndex(),
                        track.pt(),
                        track.eta(),
                        track.phi(),
                        0u, // aod::femtodreamparticle::ParticleType::kTrack,
                        1u, // 1 for antiproton
                        pidCuts,
                        track.dcaXY(),
                        childIDs,
                        0.f,
                        0.f);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilter>(cfg)};
}
