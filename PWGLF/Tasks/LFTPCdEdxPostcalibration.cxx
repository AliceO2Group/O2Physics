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
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since June 27, 2023

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta,
  aod::pidTOFmass, aod::TrackSelection, aod::TrackSelectionExtension,
  aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
  aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTOFFullPi, aod::pidTOFFullKa,
  aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct LFTPCdEdxPostcalibration {

  // dE/dx for all charged particles
  HistogramRegistry registryCh{
    "registryCh",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  // dE/dx and nsigma_{TPC} for different hadron species
  HistogramRegistry registryPi{
    "registryPi",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryKa{
    "registryKa",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryPr{
    "registryPr",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryDe{
    "registryDe",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryTr{
    "registryTr",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryHe{
    "registryHe",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // External Parameters
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 70.0f,
                                      "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{
    "minNCrossedRowsTPC", 70.0f, "min number of found TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f,
                                 "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f,
                                 "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8f, "etaMax"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.998, "Minimum V0 CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f,
                                      "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100.0f,
                                      "Maximum V0 Radius"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f,
                                        "Maximum DCA Daughters"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", 3.0f, "Maximum nsigma TOF"};
  Configurable<float> minMassK0s{"minMassK0s", 0.4f, "Minimum Mass K0s"};
  Configurable<float> maxMassK0s{"maxMassK0s", 0.5f, "Maximum Mass K0s"};
  Configurable<float> minMassLambda{"minMassLambda", 1.0f,
                                    "Minimum Mass Lambda"};
  Configurable<float> maxMassLambda{"maxMassLambda", 1.1f,
                                    "Maximum Mass Lambda"};
  Configurable<float> minReqClusterITS{
    "minReqClusterITS", 4.0f, "min number of clusters required in ITS"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.1f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.1f, "maxDCAz"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<bool> useTOFpi{"useTOFpi", true, "use TOF for pion ID"};
  Configurable<bool> useTOFpr{"useTOFpr", true, "use TOF for proton ID"};

  void init(InitContext const&)
  {
    // Raw dE/dx vs. TPC momentum
    registryCh.add(
      "dEdx_vs_Momentum", "dE/dx", HistType::kTH2F,
      {{200, -10.0, 10.0, "z#cdot p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});
    registryPi.add(
      "dEdx_vs_Momentum_Pi", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});
    registryKa.add(
      "dEdx_vs_Momentum_Ka", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});
    registryPr.add(
      "dEdx_vs_Momentum_Pr", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});
    registryDe.add(
      "dEdx_vs_Momentum_De", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});
    registryTr.add(
      "dEdx_vs_Momentum_Tr", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});
    registryHe.add(
      "dEdx_vs_Momentum_He", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {1400, 0, 1400, "dE/dx (a. u.)"}});

    // nsigma_(TPC) vs. TPC momentum
    registryPi.add(
      "nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {200, -5, 5, "n#sigma_{TPC}"}});
    registryKa.add(
      "nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {200, -5, 5, "n#sigma_{TPC}"}});
    registryPr.add(
      "nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {200, -5, 5, "n#sigma_{TPC}"}});
    registryDe.add(
      "nsigmaTPC_vs_Momentum_De", "nsigmaTPC", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {200, -5, 5, "n#sigma_{TPC}"}});
    registryTr.add(
      "nsigmaTPC_vs_Momentum_Tr", "nsigmaTPC", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {200, -5, 5, "n#sigma_{TPC}"}});
    registryHe.add(
      "nsigmaTPC_vs_Momentum_He", "nsigmaTPC", HistType::kTH2F,
      {{100, 0.0, 10.0, "p (GeV/c)"}, {200, -5, 5, "n#sigma_{TPC}"}});

    // Event Counter
    registryCh.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20.0, +20.0, "z_{vtx} (cm)"}});
  }

  // Single-Track Selection
  template <typename T1, typename C>
  bool passedSingleTrackSelection(const T1& track, const C& collision)
  {
    // Single-Track Selections
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;

    return true;
  }

  // General V0 Selections
  template <typename T1, typename C>
  bool passedV0Selection(const T1& v0, const C& collision)
  {
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) <
        v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    return true;
  }

  // K0s Selections
  template <typename T1, typename T2, typename C>
  bool passedK0Selection(const T1& v0, const T2& ntrack, const T2& ptrack,
                         const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;
    if (useTOFpi && (!ptrack.hasTOF()))
      return false;
    if (useTOFpi && (!ntrack.hasTOF()))
      return false;
    if (TMath::Abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
      return false;
    if (TMath::Abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
      return false;

    // Invariant-Mass Selection
    if (v0.mK0Short() < minMassK0s || v0.mK0Short() > maxMassK0s)
      return false;

    return true;
  }

  // Lambda Selections
  template <typename T1, typename T2, typename C>
  bool passedLambdaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                             const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;
    if (useTOFpr && (!ptrack.hasTOF()))
      return false;
    if (useTOFpi && (!ntrack.hasTOF()))
      return false;
    if (TMath::Abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
      return false;
    if (TMath::Abs(ptrack.tofNSigmaPr()) > nsigmaTOFmax)
      return false;

    // Invariant-Mass Selection
    if (v0.mLambda() < minMassLambda || v0.mLambda() > maxMassLambda)
      return false;

    return true;
  }

  // AntiLambda Selections
  template <typename T1, typename T2, typename C>
  bool passedAntiLambdaSelection(const T1& v0, const T2& ntrack,
                                 const T2& ptrack, const C& collision)
  {

    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;
    if (useTOFpr && (!ntrack.hasTOF()))
      return false;
    if (useTOFpi && (!ptrack.hasTOF()))
      return false;
    if (TMath::Abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
      return false;
    if (TMath::Abs(ntrack.tofNSigmaPr()) > nsigmaTOFmax)
      return false;

    // Invariant-Mass Selection
    if (v0.mAntiLambda() < minMassLambda || v0.mAntiLambda() > maxMassLambda)
      return false;

    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision,
               aod::V0Datas const& fullV0s, PIDTracks const& tracks)
  {
    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter
    registryCh.fill(HIST("histRecVtxZData"), collision.posZ());

    // Kaons and nuclei
    for (auto& trk : tracks) {

      if (!passedSingleTrackSelection(trk, collision))
        continue;
      if (!trk.passedITSRefit())
        continue;
      if (!trk.passedTPCRefit())
        continue;
      if (trk.itsNCls() < minReqClusterITS)
        continue;
      if (TMath::Abs(trk.dcaXY()) > maxDCAxy)
        continue;
      if (TMath::Abs(trk.dcaZ()) > maxDCAz)
        continue;
      if (trk.itsChi2NCl() > maxChi2ITS)
        continue;

      // Charged Particles
      registryCh.fill(HIST("dEdx_vs_Momentum"), trk.sign() * trk.tpcInnerParam(), trk.tpcSignal());

      // Kaons
      if (trk.hasTOF() && TMath::Abs(trk.tofNSigmaKa()) < 2.0) {
        registryKa.fill(HIST("dEdx_vs_Momentum_Ka"), trk.tpcInnerParam(),
                        trk.tpcSignal());
        registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), trk.tpcInnerParam(),
                        trk.tpcNSigmaKa());
      }

      // Selection of high dE/dx Objects
      if (trk.tpcNSigmaDe() < -4.0)
        continue;

      // Deuterons
      if (trk.hasTOF() && TMath::Abs(trk.tofNSigmaDe()) < 3.0) {
        registryDe.fill(HIST("dEdx_vs_Momentum_De"), trk.tpcInnerParam(),
                        trk.tpcSignal());
        registryDe.fill(HIST("nsigmaTPC_vs_Momentum_De"), trk.tpcInnerParam(),
                        trk.tpcNSigmaDe());
      }

      // Heavier Nuclei
      registryTr.fill(HIST("dEdx_vs_Momentum_Tr"), trk.tpcInnerParam(),
                      trk.tpcSignal());
      registryTr.fill(HIST("nsigmaTPC_vs_Momentum_Tr"), trk.tpcInnerParam(),
                      trk.tpcNSigmaTr());
      registryHe.fill(HIST("dEdx_vs_Momentum_He"), 2.0 * trk.tpcInnerParam(),
                      trk.tpcSignal());
      registryHe.fill(HIST("nsigmaTPC_vs_Momentum_He"),
                      2.0 * trk.tpcInnerParam(), trk.tpcNSigmaHe());
    }

    // Loop over Reconstructed V0s
    for (auto& v0 : fullV0s) {

      // Standard V0 Selections
      if (!passedV0Selection(v0, collision)) {
        continue;
      }

      if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
        continue;
      }

      // Positive and Negative Tracks
      const auto& posTrack = v0.posTrack_as<PIDTracks>();
      const auto& negTrack = v0.negTrack_as<PIDTracks>();

      if (!posTrack.passedTPCRefit())
        continue;
      if (!negTrack.passedTPCRefit())
        continue;

      // K0s Selection
      if (passedK0Selection(v0, negTrack, posTrack, collision)) {
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), negTrack.tpcInnerParam(),
                        negTrack.tpcSignal());
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), posTrack.tpcInnerParam(),
                        posTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        negTrack.tpcInnerParam(), negTrack.tpcNSigmaPi());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        posTrack.tpcInnerParam(), posTrack.tpcNSigmaPi());
      }

      // Lambda Selection
      if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {
        registryPr.fill(HIST("dEdx_vs_Momentum_Pr"), posTrack.tpcInnerParam(),
                        posTrack.tpcSignal());
        registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                        posTrack.tpcInnerParam(), posTrack.tpcNSigmaPr());
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), negTrack.tpcInnerParam(),
                        negTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        negTrack.tpcInnerParam(), negTrack.tpcNSigmaPi());
      }

      // AntiLambda Selection
      if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), posTrack.tpcInnerParam(),
                        posTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        posTrack.tpcInnerParam(), posTrack.tpcNSigmaPi());
        registryPr.fill(HIST("dEdx_vs_Momentum_Pr"), negTrack.tpcInnerParam(),
                        negTrack.tpcSignal());
        registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                        negTrack.tpcInnerParam(), negTrack.tpcNSigmaPr());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LFTPCdEdxPostcalibration>(cfgc)};
}
