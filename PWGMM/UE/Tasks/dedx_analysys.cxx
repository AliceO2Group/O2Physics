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
/// \author Paola Vargas
/// \since January 8, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta,
  aod::pidTOFmass, aod::TrackSelection, aod::TrackSelectionExtension,
  aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
  aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTOFFullPi, aod::pidTOFFullKa,
  aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullEl>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;

struct dedx_analysys {

  // dE/dx for all charged particles
  HistogramRegistry registryDeDx{
    "registryDeDx",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Corrections
  // constexpr double Correc[8]{54.3344, 55.1277, 56.0811, 56.7974, 56.9533, 56.4622, 55.8873, 55.1449};
  // constexpr double Correc[9]{-0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8};

  // Configurable Parameters
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
  Configurable<float> v0cospaMin{"v0cospaMin", 0.998f, "Minimum V0 CosPA"};
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
  Configurable<float> minMassGamma{"minMassGamma", 0.000922f,
                                   "Minimum Mass Gamma"};
  Configurable<float> maxMassGamma{"maxMassGamma", 0.002022f,
                                   "Maximum Mass Gamma"};
  Configurable<float> minReqClusterITS{
    "minReqClusterITS", 4.0f, "min number of clusters required in ITS"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.1f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.1f, "maxDCAz"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  // Histograms names
  static constexpr std::string_view Armenteros[5] = {"Armenteros", "Armenteros_K0S", "Armenteros_Lambda", "Armenteros_AntiLambda", "Armenteros_Gamma"};
  static constexpr std::string_view Qt_vs_alpha[5] = {"Qt_vs_alpha", "Qt_vs_alpha_K0S", "Qt_vs_alpha_Lambda", "Qt_vs_alpha_AntiLambda", "Qt_vs_alpha_Gamma"};
  static constexpr std::string_view dEdx_vs_Momentum[5] = {"dEdx_vs_Momentum_all", "dEdx_vs_Momentum_02", "dEdx_vs_Momentum_0204", "dEdx_vs_Momentum_0406", "dEdx_vs_Momentum_0608"};
  static constexpr std::string_view dEdx_vs_Momentum_neg[4] = {"dEdx_vs_Momentum_02neg", "dEdx_vs_Momentum_0204neg", "dEdx_vs_Momentum_0406neg", "dEdx_vs_Momentum_0608neg"};
  static constexpr std::string_view dEdx_vs_Momentum_pos[4] = {"dEdx_vs_Momentum_02pos", "dEdx_vs_Momentum_0204pos", "dEdx_vs_Momentum_0406pos", "dEdx_vs_Momentum_0608pos"};
  static constexpr std::string_view dEdx_vs_Momentum_v0[3] = {"dEdx_vs_Momentum_Pi_v0", "dEdx_vs_Momentum_Pr_v0", "dEdx_vs_Momentum_El_v0"};
  static constexpr std::string_view hInvMass[4] = {"InvMass_K0S", "InvMass_Lambda", "InvMass_AntiLambda", "InvMass_Gamma"};
  static constexpr double EtaCut[9] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  static constexpr double Correction[8] = {54.3344, 55.1277, 56.0811, 56.7974, 56.9533, 56.4622, 55.8873, 55.1449};

  void init(InitContext const&)
  {

    // MIP for pions
    registryDeDx.add(
      "dEdxMIP_vs_eta", "dE/dx", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0.0, 600.0, "dE/dx MIP (a. u.)"}});

    registryDeDx.add(
      "dEdxMIP_vs_phi", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 6.4, "#phi"}, {100, 0.0, 600.0, "dE/dx MIP (a. u.)"}});

    registryDeDx.add(
      "dEdxMIP_vs_eta_AfterCorr", "dE/dx", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0.0, 600.0, "dE/dx MIP (a. u.)"}});

    ////////////////////////////////
    registryDeDx.add(hInvMass[0].data(), "mass", HistType::kTH1F,
                     {{100, 400, 600, "m (MeV/c)"}});
    registryDeDx.add(hInvMass[1].data(), "mass", HistType::kTH1F,
                     {{100, 1.08, 1.25, "m (GeV/c)"}});
    registryDeDx.add(hInvMass[2].data(), "mass", HistType::kTH1F,
                     {{100, 1.08, 1.25, "m (GeV/c)"}});
    registryDeDx.add(hInvMass[3].data(), "mass", HistType::kTH1F,
                     {{100, 0, 2, "m (MeV/c)"}});
    // Armenteros plot and De/Dx for eta cut inclusive
    for (int i = 0; i < 5; ++i) {
      registryDeDx.add(Armenteros[i].data(), Qt_vs_alpha[i].data(), HistType::kTH2F,
                       {{100, -1., 1., "#alpha (a. u.)"}, {100, 0.0, 0.3, "q_T (GeV/c)"}});
      registryDeDx.add(dEdx_vs_Momentum[i].data(), "dE/dx", HistType::kTH2F,
                       {{100, -20.0, 20.0, "#it{p}/Z (GeV/c)"}, {100, 0.0, 600.0, "dE/dx (a. u.)"}});
    }
    // De/Dx for eta cut negative and positive
    for (int i = 0; i < 4; ++i) {

      registryDeDx.add(dEdx_vs_Momentum_neg[i].data(), "dE/dx", HistType::kTH2F,
                       {{100, -20.0, 20.0, "#it{p}/Z (GeV/c)"}, {100, 0.0, 600.0, "dE/dx (a. u.)"}});

      registryDeDx.add(dEdx_vs_Momentum_pos[i].data(), "dE/dx", HistType::kTH2F,
                       {{100, -20.0, 20.0, "#it{p}/Z (GeV/c)"}, {100, 0.0, 600.0, "dE/dx (a. u.)"}});
    }

    // De/Dx for v0 particles
    for (int i = 0; i < 3; ++i) {

      registryDeDx.add(dEdx_vs_Momentum_v0[i].data(), "dE/dx", HistType::kTH2F,
                       {{100, -20.0, 20.0, "#it{p}/Z (GeV/c)"}, {100, 0.0, 600.0, "dE/dx (a. u.)"}});
    }

    // Event Counter
    registryDeDx.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{100, -20.0, +20.0, "z_{vtx} (cm)"}});
  }

  // Single-Track Selection
  template <typename T1, typename C>
  bool passedSingleTrackSelection(const T1& track, const C& /*collision*/)
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
  bool passedV0Selection(const T1& v0, const C& /*collision*/)
  {
    if (v0.v0cosPA() < v0cospaMin)
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

  // Gamma Selections
  template <typename T1, typename T2, typename C>
  bool passedGammaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
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
      if (TMath::Abs(ptrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mGamma() < minMassGamma || v0.mGamma() > maxMassGamma)
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
    registryDeDx.fill(HIST("histRecVtxZData"), collision.posZ());

    // Centrality
    float centrality = collision.centFT0C();
    if (centrality < 0.0 || centrality > 100.0)
      centrality = 1.0;

    // Kaons
    for (auto& trk : tracks) {

      if (!passedSingleTrackSelection(trk, collision))
        continue;
      if (!trk.passedTPCRefit())
        continue;
      float signedP = trk.sign() * trk.tpcInnerParam();

      // DeDx all particles
      registryDeDx.fill(HIST(dEdx_vs_Momentum[0]), signedP, trk.tpcSignal());

      ////////////////////////////////

      // MIP for pions

      // TODO: 0.25 TO 0.35
      // range 0.2 to 0.6
      // float pion_m = 0.139; //GeV
      if (trk.tpcInnerParam() >= 0.25 && trk.tpcInnerParam() <= 0.35) {
        registryDeDx.fill(HIST("dEdxMIP_vs_eta"), trk.eta(), trk.tpcSignal());
        registryDeDx.fill(HIST("dEdxMIP_vs_phi"), trk.phi(), trk.tpcSignal());
        // After calibration

        for (int i = 0; i < 8; ++i) {
          if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
            registryDeDx.fill(HIST("dEdxMIP_vs_eta_AfterCorr"), trk.eta(), trk.tpcSignal() * 50 / Correction[i]);
          }
        }
      }

      // After calibration
      for (int i = 0; i < 8; ++i) {
        if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
          if (i == 0) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_neg[0]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[4]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 1) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_neg[1]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[3]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 2) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_neg[2]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[2]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 3) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_neg[3]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[1]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 4) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_pos[0]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[1]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 5) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_pos[1]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[2]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 6) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_pos[2]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[3]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          } else if (i == 7) {
            registryDeDx.fill(HIST(dEdx_vs_Momentum_pos[3]), signedP, trk.tpcSignal() * 50 / Correction[i]);
            registryDeDx.fill(HIST(dEdx_vs_Momentum[4]), signedP, trk.tpcSignal() * 50 / Correction[i]);
          }
        }
      }
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

      float px_v0 = v0.px();
      float py_v0 = v0.py();
      float pz_v0 = v0.pz();
      float p_v0 = sqrt(px_v0 * px_v0 + py_v0 * py_v0 + pz_v0 * pz_v0);

      float px_pos = posTrack.px();
      float py_pos = posTrack.py();
      float pz_pos = posTrack.pz();

      float px_neg = negTrack.px();
      float py_neg = negTrack.py();
      float pz_neg = negTrack.pz();

      const float protonMass = 0.938; // GeV/c^2
      const float piMass = 0.13957;   // GeV/c^2
      const float K0SMass = 493.677;  // MeV/c^2
      const float LambdaMass = 1.115; // GeV/c^2
      const float GammaMass = 1.022;  // MeV/c^2
      const float eMass = 0.511;      // MeV/c^2

      // float E_pos_man = std::sqrt(px_pos * px_pos + py_pos * py_pos + pz_pos * pz_pos + piMass * piMass);

      //-------------------Armenteros plots--------

      float pl_pos = (px_pos * px_v0 + py_pos * py_v0 + pz_pos * pz_v0) / p_v0;
      float pl_neg = (px_neg * px_v0 + py_neg * py_v0 + pz_neg * pz_v0) / p_v0;

      float alpha = (pl_pos - pl_neg) / (pl_pos + pl_neg);
      float p_pos = sqrt(px_pos * px_pos + py_pos * py_pos + pz_pos * pz_pos);
      float qt = sqrt(p_pos * p_pos - pl_pos * pl_pos);

      registryDeDx.fill(HIST(Armenteros[0]), alpha, qt);

      //-------------------------------------------

      // K0s Selection
      if (passedK0Selection(v0, negTrack, posTrack, collision)) {
        float E_pos_pi = posTrack.energy(piMass);
        float E_neg_pi = negTrack.energy(piMass);

        float InvMass = sqrt((E_neg_pi + E_pos_pi) * (E_neg_pi + E_pos_pi) - ((px_neg + px_pos) * (px_neg + px_pos) + (py_neg + py_pos) * (py_neg + py_pos) + (pz_neg + pz_pos) * (pz_neg + pz_pos)));

        if (TMath::Abs(InvMass * 1000 - K0SMass) > 10) {
          continue;
        }

        registryDeDx.fill(HIST(hInvMass[0]), InvMass * 1000);
        registryDeDx.fill(HIST(Armenteros[1]), alpha, qt);

        for (int i = 0; i < 8; ++i) {
          if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[0]), signedPneg, negTrack.tpcSignal() * 50 / Correction[i]);
          }
          if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[0]), signedPpos, posTrack.tpcSignal() * 50 / Correction[i]);
          }
        }
      }

      // Lambda Selection
      if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {
        float E_pos_p = posTrack.energy(protonMass);
        float E_neg_pi = negTrack.energy(piMass);

        float InvMass = sqrt((E_neg_pi + E_pos_p) * (E_neg_pi + E_pos_p) - ((px_neg + px_pos) * (px_neg + px_pos) + (py_neg + py_pos) * (py_neg + py_pos) + (pz_neg + pz_pos) * (pz_neg + pz_pos)));

        if (TMath::Abs(InvMass - LambdaMass) > 0.01) {
          continue;
        }

        registryDeDx.fill(HIST(hInvMass[1]), InvMass);
        registryDeDx.fill(HIST(Armenteros[2]), alpha, qt);

        for (int i = 0; i < 8; ++i) {
          if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[0]), signedPneg, negTrack.tpcSignal() * 50 / Correction[i]);
          }
          if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[1]), signedPpos, posTrack.tpcSignal() * 50 / Correction[i]);
          }
        }
      }

      // AntiLambda Selection
      if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {

        float E_pos_pi = posTrack.energy(piMass);
        float E_neg_p = negTrack.energy(protonMass);

        float InvMass = sqrt((E_neg_p + E_pos_pi) * (E_neg_p + E_pos_pi) - ((px_neg + px_pos) * (px_neg + px_pos) + (py_neg + py_pos) * (py_neg + py_pos) + (pz_neg + pz_pos) * (pz_neg + pz_pos)));

        if (TMath::Abs(InvMass - LambdaMass) > 0.01) {
          continue;
        }

        registryDeDx.fill(HIST(hInvMass[2]), InvMass);
        registryDeDx.fill(HIST(Armenteros[3]), alpha, qt);

        for (int i = 0; i < 8; ++i) {
          if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[1]), signedPneg, negTrack.tpcSignal() * 50 / Correction[i]);
          }
          if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[0]), signedPpos, posTrack.tpcSignal() * 50 / Correction[i]);
          }
        }
      }

      // Gamma Selection
      if (passedGammaSelection(v0, negTrack, posTrack, collision)) {

        float E_pos_e = posTrack.energy(eMass);
        float E_neg_e = negTrack.energy(eMass);

        float InvMass = sqrt((E_neg_e + E_pos_e) * (E_neg_e + E_pos_e) - ((px_neg + px_pos) * (px_neg + px_pos) + (py_neg + py_pos) * (py_neg + py_pos) + (pz_neg + pz_pos) * (pz_neg + pz_pos)));

        if (TMath::Abs(InvMass - GammaMass) > 10) {
          continue;
        }

        registryDeDx.fill(HIST(hInvMass[3]), InvMass);
        registryDeDx.fill(HIST(Armenteros[4]), alpha, qt);

        for (int i = 0; i < 8; ++i) {
          if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[2]), signedPneg, negTrack.tpcSignal() * 50 / Correction[i]);
          }
          if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {

            registryDeDx.fill(HIST(dEdx_vs_Momentum_v0[2]), signedPpos, posTrack.tpcSignal() * 50 / Correction[i]);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<dedx_analysys>(cfgc)};
}
