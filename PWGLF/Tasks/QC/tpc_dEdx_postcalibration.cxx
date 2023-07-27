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

struct tpc_dEdx_postcalibration {

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
  Configurable<float> minNClsTPCdEdx{
    "minNClsTPCdEdx", 50.0f, "min number of TPC clusters for PID"};
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
  Configurable<float> maxMassK0s{"maxMassK0s", 0.6f, "Maximum Mass K0s"};
  Configurable<float> minMassLambda{"minMassLambda", 1.1f,
                                    "Minimum Mass Lambda"};
  Configurable<float> maxMassLambda{"maxMassLambda", 1.2f,
                                    "Maximum Mass Lambda"};
  Configurable<float> minReqClusterITS{
    "minReqClusterITS", 4.0f, "min number of clusters required in ITS"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.1f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.1f, "maxDCAz"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<bool> useTOFpi{"useTOFpi", true, "use TOF for pion ID"};
  Configurable<bool> useTOFpr{"useTOFpr", true, "use TOF for proton ID"};
  Configurable<bool> doContaminations{"doContaminations", true, "Flag to produce the plots for the contaminations"};
  Configurable<bool> usePt{"usePt", true, "Flag to use PT instead of the TPCInnerParam"};
  ConfigurableAxis pBins{"pBins", {200, -10.0, 10.0}, "Binning in p"};
  ConfigurableAxis dEdxBins{"dEdxBins", {3000, 0.f, 1500.f}, "Binning in dE/dx"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {200, -5, 5}, "Binning in nsigma"};

  void init(InitContext const&)
  {
    AxisSpec pAxis{pBins, "z#upoint p (GeV/c)"};
    if (usePt) {
      pAxis.title = "#it{p}_{T} (GeV/c)";
    }
    const AxisSpec dedxAxis{dEdxBins, "d#it{E}/d#it{x} Arb. units"};
    const AxisSpec nsigmaAxis{nsigmaBins, "n#sigma_{TPC}"};

    // Raw dE/dx vs. TPC momentum
    registryCh.add(
      "dEdx_vs_Momentum", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryPi.add(
      "dEdx_vs_Momentum_Pi", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryKa.add(
      "dEdx_vs_Momentum_Ka", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryPr.add(
      "dEdx_vs_Momentum_Pr", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryDe.add(
      "dEdx_vs_Momentum_De", "dE/dx", HistType::kTH2F,
      {{2000, -3.0, 3.0, "z#upoint p (GeV/c)"}, dEdxBins});
    registryTr.add(
      "dEdx_vs_Momentum_Tr", "dE/dx", HistType::kTH2F,
      {{200, -3.0, 3.0, "z#upoint p (GeV/c)"}, dEdxBins});
    registryHe.add(
      "dEdx_vs_Momentum_He", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});

    // nsigma_(TPC) vs. TPC momentum
    registryPi.add(
      "nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryKa.add(
      "nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryPr.add(
      "nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryDe.add(
      "nsigmaTPC_vs_Momentum_De", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryTr.add(
      "nsigmaTPC_vs_Momentum_Tr", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryHe.add(
      "nsigmaTPC_vs_Momentum_He", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    if (doContaminations) {
      registryPi.add(
        "nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F,
        {pAxis, nsigmaAxis});
      registryPi.add(
        "nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F,
        {pAxis, nsigmaAxis});

      registryKa.add(
        "nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F,
        {pAxis, nsigmaAxis});
      registryKa.add(
        "nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F,
        {pAxis, nsigmaAxis});

      registryPr.add(
        "nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F,
        {pAxis, nsigmaAxis});
      registryPr.add(
        "nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F,
        {pAxis, nsigmaAxis});
    }
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
    // if (track.tpcSignalN() < minNClsTPCdEdx)
    // return false;
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
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospaMin)
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

    if (ptrack.tpcInnerParam() > 0.6) {
      if (!ptrack.hasTOF())
        return false;
      if (TMath::Abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

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

    if (ptrack.tpcInnerParam() > 0.6) {
      if (!ptrack.hasTOF())
        return false;
      if (TMath::Abs(ptrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

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

    if (ptrack.tpcInnerParam() > 0.6) {
      if (!ptrack.hasTOF())
        return false;
      if (TMath::Abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

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
      if (!trk.passedTPCRefit())
        continue;
      float signedP = trk.sign() * trk.tpcInnerParam();
      if (usePt) {
        signedP = trk.sign() * trk.pt();
      }

      // Charged Particles
      registryCh.fill(HIST("dEdx_vs_Momentum"), signedP,
                      trk.tpcSignal());

      // Kaons
      if (trk.tpcInnerParam() > 0.4 && trk.hasTOF() && TMath::Abs(trk.tofNSigmaKa()) < 2.0) {
        registryKa.fill(HIST("dEdx_vs_Momentum_Ka"), signedP,
                        trk.tpcSignal());
        registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedP,
                        trk.tpcNSigmaKa());

        if (doContaminations) {
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pi"), signedP,
                          trk.tpcNSigmaPi());
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedP,
                          trk.tpcNSigmaPr());
        }
      }

      if (trk.tpcInnerParam() < 0.4) {
        registryKa.fill(HIST("dEdx_vs_Momentum_Ka"), signedP,
                        trk.tpcSignal());
        registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedP,
                        trk.tpcNSigmaKa());

        if (doContaminations) {
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pi"), signedP,
                          trk.tpcNSigmaPi());
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedP,
                          trk.tpcNSigmaPr());
        }
      }

      // Selection of high dE/dx Objects
      if (trk.tpcNSigmaDe() < -4.0)
        continue;

      // Deuterons
      if (trk.tpcInnerParam() > 1.0 && trk.hasTOF() && TMath::Abs(trk.tofNSigmaDe()) < 3.0) {
        registryDe.fill(HIST("dEdx_vs_Momentum_De"), signedP,
                        trk.tpcSignal());
        registryDe.fill(HIST("nsigmaTPC_vs_Momentum_De"), signedP,
                        trk.tpcNSigmaDe());
      }

      if (trk.tpcInnerParam() < 1.0) {
        registryDe.fill(HIST("dEdx_vs_Momentum_De"), signedP,
                        trk.tpcSignal());
        registryDe.fill(HIST("nsigmaTPC_vs_Momentum_De"), signedP,
                        trk.tpcNSigmaDe());
      }

      // Heavier Nuclei
      registryTr.fill(HIST("dEdx_vs_Momentum_Tr"), signedP,
                      trk.tpcSignal());
      registryTr.fill(HIST("nsigmaTPC_vs_Momentum_Tr"), signedP,
                      trk.tpcNSigmaTr());
      registryHe.fill(HIST("dEdx_vs_Momentum_He"), 2.0 * signedP,
                      trk.tpcSignal());
      registryHe.fill(HIST("nsigmaTPC_vs_Momentum_He"),
                      2.0 * signedP, trk.tpcNSigmaHe());
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

      float signedPpos = posTrack.sign() * posTrack.tpcInnerParam();
      float signedPneg = negTrack.sign() * negTrack.tpcInnerParam();
      if (usePt) {
        signedPpos = posTrack.sign() * posTrack.pt();
        signedPneg = negTrack.sign() * negTrack.pt();
      }

      // K0s Selection
      if (passedK0Selection(v0, negTrack, posTrack, collision)) {
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPneg, negTrack.tpcSignal());
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPpos, posTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPneg, negTrack.tpcNSigmaPi());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPpos, posTrack.tpcNSigmaPi());
        if (doContaminations) {
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"),
                          signedPneg, negTrack.tpcNSigmaKa());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                          signedPneg, negTrack.tpcNSigmaPr());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"),
                          signedPpos, posTrack.tpcNSigmaKa());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                          signedPpos, posTrack.tpcNSigmaPr());
        }
      }

      // Lambda Selection
      if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {
        registryPr.fill(HIST("dEdx_vs_Momentum_Pr"), signedPpos, posTrack.tpcSignal());
        registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                        signedPpos, posTrack.tpcNSigmaPr());
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPneg, negTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPneg, negTrack.tpcNSigmaPi());

        if (doContaminations) {
          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                          signedPpos, posTrack.tpcNSigmaPi());
          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Ka"),
                          signedPpos, posTrack.tpcNSigmaKa());

          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                          signedPneg, negTrack.tpcNSigmaPr());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"),
                          signedPneg, negTrack.tpcNSigmaKa());
        }
      }

      // AntiLambda Selection
      if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPpos, posTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPpos, posTrack.tpcNSigmaPi());
        registryPr.fill(HIST("dEdx_vs_Momentum_Pr"), signedPneg, negTrack.tpcSignal());
        registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                        signedPneg, negTrack.tpcNSigmaPr());

        if (doContaminations) {
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                          signedPpos, posTrack.tpcNSigmaPi());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"),
                          signedPpos, posTrack.tpcNSigmaKa());

          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                          signedPneg, negTrack.tpcNSigmaPr());
          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Ka"),
                          signedPneg, negTrack.tpcNSigmaKa());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tpc_dEdx_postcalibration>(cfgc)};
}
