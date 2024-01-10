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
/// \since September 26, 2023

#include <cmath>
#include <TMath.h>
#include <TPDGCode.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

using MCTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::McTrackLabels>;

struct vzero_cascade_absorption {

  // QC Histograms
  HistogramRegistry registryQC{
    "registryQC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: Data
  HistogramRegistry registryData{
    "registryData",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: MC
  HistogramRegistry registryMC{
    "registryMC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Configurable Parameters
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of found TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> etaMin{"etaMin", -0.5f, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.5f, "etaMax"};
  Configurable<float> pMin_k0postrack{"pMin_k0postrack", 0.3f, "Min Momentum K0 pos track"};
  Configurable<float> pMax_k0postrack{"pMax_k0postrack", 5.0f, "Max Momentum K0 pos track"};
  Configurable<float> pMin_k0negtrack{"pMin_k0negtrack", 0.3f, "Min Momentum K0 neg track"};
  Configurable<float> pMax_k0negtrack{"pMax_k0negtrack", 5.0f, "Max Momentum K0 neg track"};
  Configurable<float> pMin_Lambda_proton{"pMin_Lambda_proton", 0.3f, "Min Momentum proton from Lambda"};
  Configurable<float> pMax_Lambda_proton{"pMax_Lambda_proton", 5.0f, "Max Momentum proton from Lambda"};
  Configurable<float> pMin_Lambda_pion{"pMin_Lambda_pion", 0.2f, "Min Momentum pion from Lambda"};
  Configurable<float> pMax_Lambda_pion{"pMax_Lambda_pion", 1.0f, "Max Momentum pion from Lambda"};
  Configurable<bool> requirehitsITS{"requirehitsITS", true, "require ITS hits for daughters"};
  Configurable<std::vector<float>> hit_requirement_before_target{"hit_requirement_before_target", {0, 0, 0, 1, 1, 1, 1}, "ITS Hits before target"};
  Configurable<std::vector<float>> hit_requirement_after_target{"hit_requirement_after_target", {0, 0, 0, 0, 0, 1, 1}, "ITS Hits after target"};
  Configurable<bool> useTOF{"useTOF", true, "use TOF for PID"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 40.0f, "Maximum V0 Radius"};
  Configurable<float> minimumCascRadius{"minimumCascRadius", 0.5f, "Minimum Cascade Radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 40.0f, "Maximum Cascade Radius"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<float> dcaBachelorToPVmin{"dcaBachelorToPVmin", 0.1f, "Minimum DCA bachelor To PV"};
  Configurable<float> dcaV0ToPVmin{"dcaV0ToPVmin", 0.1f, "Minimum DCA V0 To PV for cascades"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 CosPA"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum Cascade CosPA"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA Daughters"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.5f, "Maximum DCA Cascade Daughters"};
  Configurable<float> Rmin_beforeAbs{"Rmin_beforeAbs", 10.0f, "Rmin before target"};
  Configurable<float> Rmax_beforeAbs{"Rmax_beforeAbs", 15.0f, "Rmax before target"};
  Configurable<float> Rmin_afterAbs{"Rmin_afterAbs", 26.0f, "Rmin after target"};
  Configurable<float> Rmax_afterAbs{"Rmax_afterAbs", 31.0f, "Rmax after target"};

  void init(InitContext const&)
  {
    // Histograms
    registryQC.add("event_counter", "event counter", HistType::kTH1F, {{4, 0, 4, "number of events in data (first 2 bins) and mc (last 2 bins)"}});

    // K0 short
    registryData.add("K0_before_target_data", "K0 before target data", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {240, 0.44, 0.56, "m (GeV/c^{2})"}});
    registryData.add("K0_after_target_data", "K0 after target data", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {240, 0.44, 0.56, "m (GeV/c^{2})"}});
    registryMC.add("K0_before_target_mc", "K0 before target mc", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {240, 0.44, 0.56, "m (GeV/c^{2})"}});
    registryMC.add("K0_after_target_mc", "K0 after target mc", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {240, 0.44, 0.56, "m (GeV/c^{2})"}});

    // Lambda and Antilambda
    registryData.add("Lambda_before_target_data", "Lambda before target data", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryData.add("Lambda_after_target_data", "Lambda after target data", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryData.add("AntiLambda_before_target_data", "AntiLambda before target data", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryData.add("AntiLambda_after_target_data", "AntiLambda after target data", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryMC.add("Lambda_before_target_mc", "Lambda before target mc", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryMC.add("Lambda_after_target_mc", "Lambda after target mc", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryMC.add("AntiLambda_before_target_mc", "AntiLambda before target mc", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});
    registryMC.add("AntiLambda_after_target_mc", "AntiLambda after target mc", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.09, 1.14, "m (GeV/c^{2})"}});

    // Generated Distributions
    registryMC.add("Lambda_before_target_mc_gen", "Lambda before target mc gen", HistType::kTH1F, {{200, 0.0, 10.0, "p (GeV/c)"}});
    registryMC.add("Lambda_after_target_mc_gen", "Lambda after target mc gen", HistType::kTH1F, {{200, 0.0, 10.0, "p (GeV/c)"}});
    registryMC.add("AntiLambda_before_target_mc_gen", "AntiLambda before target mc gen", HistType::kTH1F, {{200, 0.0, 10.0, "p (GeV/c)"}});
    registryMC.add("AntiLambda_after_target_mc_gen", "AntiLambda after target mc gen", HistType::kTH1F, {{200, 0.0, 10.0, "p (GeV/c)"}});

    // Resolution
    registryMC.add("K0_Rresolution_before_target", "K0 Rresolution before target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {500, -5, 5, "#Delta R (cm)"}});
    registryMC.add("K0_Rresolution_after_target", "K0 Rresolution after target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {500, -5, 5, "#Delta R (cm)"}});
    registryMC.add("Lambda_Rresolution_before_target", "Lambda Rresolution before target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {500, -5, 5, "#Delta R (cm)"}});
    registryMC.add("Lambda_Rresolution_after_target", "Lambda Rresolution after target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {500, -5, 5, "#Delta R (cm)"}});
    registryMC.add("AntiLambda_Rresolution_before_target", "AntiLambda Rresolution before target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {500, -5, 5, "#Delta R (cm)"}});
    registryMC.add("AntiLambda_Rresolution_after_target", "AntiLambda Rresolution after target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {500, -5, 5, "#Delta R (cm)"}});

    /*
    registryData.add("Csi_before_target", "Csi_before_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.305, 1.34, "m (GeV/c)"}});
    registryData.add("Csi_after_target", "Csi_after_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.305, 1.34, "m (GeV/c)"}});
    registryData.add("AntiCsi_before_target", "AntiCsi_before_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.3, 1.35, "m (GeV/c)"}});
    registryData.add("AntiCsi_after_target", "AntiCsi_after_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.3, 1.35, "m (GeV/c)"}});
    registryData.add("Omega_before_target", "Omega_before_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.6, 1.75, "m (GeV/c)"}});
    registryData.add("Omega_after_target", "Omega_after_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.6, 1.75, "m (GeV/c)"}});
    registryData.add("AntiOmega_before_target", "AntiOmega_before_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.6, 1.75, "m (GeV/c)"}});
    registryData.add("AntiOmega_after_target", "AntiOmega_after_target", HistType::kTH2F, {{200, 0.0, 10.0, "p (GeV/c)"}, {200, 1.6, 1.75, "m (GeV/c)"}});*/
  }

  bool hasHitOnITSlayer(uint8_t itsClsmap, int layer)
  {

    unsigned char test_bit = 1 << layer;
    return (itsClsmap & test_bit);
  }

  // Single-Track Selection
  template <typename T1, typename C>
  bool passedSingleTrackSelection(const T1& track, const C& collision)
  {
    // Single-Track Selections
    if (requirehitsITS && (!track.hasITS()))
      return false;
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
    if (useTOF && (!track.hasTOF()))
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

    // Momentum K0 Daughters
    float p_k0postrack = TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos());
    float p_k0negtrack = TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg());

    // Momentum Interval
    if (p_k0postrack < pMin_k0postrack || p_k0postrack > pMax_k0postrack)
      return false;
    if (p_k0negtrack < pMin_k0negtrack || p_k0negtrack > pMax_k0negtrack)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (v0.dcapostopv() < dcapostoPVmin)
      return false;
    if (v0.dcanegtopv() < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (useTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

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

    // Momentum Lambda Daughters
    float p_Lambdapostrack = TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos());
    float p_Lambdanegtrack = TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg());

    // Momentum Interval
    if (p_Lambdapostrack < pMin_Lambda_proton || p_Lambdapostrack > pMax_Lambda_proton)
      return false;
    if (p_Lambdanegtrack < pMin_Lambda_pion || p_Lambdanegtrack > pMax_Lambda_pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (v0.dcapostopv() < dcapostoPVmin)
      return false;
    if (v0.dcanegtopv() < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (useTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    return true;
  }

  // AntiLambda Selections
  template <typename T1, typename T2, typename C>
  bool passedAntiLambdaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                                 const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    // Momentum Lambda Daughters
    float p_Lambdapostrack = TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos());
    float p_Lambdanegtrack = TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg());

    // Momentum Interval
    if (p_Lambdapostrack < pMin_Lambda_pion || p_Lambdapostrack > pMax_Lambda_pion)
      return false;
    if (p_Lambdanegtrack < pMin_Lambda_proton || p_Lambdanegtrack > pMax_Lambda_proton)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (v0.dcapostopv() < dcapostoPVmin)
      return false;
    if (v0.dcanegtopv() < dcanegtoPVmin)
      return false;

    // PID Selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (useTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }

    return true;
  }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, aod::V0Datas const& fullV0s, FullTracks const& tracks)
  {

    // Event Counter (before event sel)
    registryQC.fill(HIST("event_counter"), 0.5);

    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter (after event sel)
    registryQC.fill(HIST("event_counter"), 1.5);

    auto hit_ITS_before_target = static_cast<std::vector<float>>(hit_requirement_before_target);
    auto hit_ITS_after_target = static_cast<std::vector<float>>(hit_requirement_after_target);

    // Loop over Reconstructed V0s
    for (auto& v0 : fullV0s) {

      // Positive and Negative Tracks
      const auto& posTrack = v0.posTrack_as<FullTracks>();
      const auto& negTrack = v0.negTrack_as<FullTracks>();

      // Require TPC Refit
      if (!posTrack.passedTPCRefit())
        continue;
      if (!negTrack.passedTPCRefit())
        continue;

      // ITS Requirement
      bool satisfyITSreq = true;
      if (requirehitsITS) {
        for (int i = 0; i < 7; i++) {
          if (hit_ITS_before_target[i] > 0 && !hasHitOnITSlayer(posTrack.itsClusterMap(), i)) {
            satisfyITSreq = false;
            break;
          }
          if (hit_ITS_after_target[i] > 0 && !hasHitOnITSlayer(negTrack.itsClusterMap(), i)) {
            satisfyITSreq = false;
            break;
          }
        }
      }

      if (requirehitsITS && (!satisfyITSreq))
        continue;

      // K0 Short
      if (passedK0Selection(v0, negTrack, posTrack, collision)) {

        // Before Target
        if (v0.v0radius() > Rmin_beforeAbs && v0.v0radius() < Rmax_beforeAbs)
          registryData.fill(HIST("K0_before_target_data"), v0.p(), v0.mK0Short());

        // After Target
        if (v0.v0radius() > Rmin_afterAbs && v0.v0radius() < Rmax_afterAbs)
          registryData.fill(HIST("K0_after_target_data"), v0.p(), v0.mK0Short());
      }

      // Lambda
      if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {

        // Before Target
        if (v0.v0radius() > Rmin_beforeAbs && v0.v0radius() < Rmax_beforeAbs)
          registryData.fill(HIST("Lambda_before_target_data"), v0.p(), v0.mLambda());

        // After Target
        if (v0.v0radius() > Rmin_afterAbs && v0.v0radius() < Rmax_afterAbs)
          registryData.fill(HIST("Lambda_after_target_data"), v0.p(), v0.mLambda());
      }

      // AntiLambda
      if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {

        // Before Target
        if (v0.v0radius() > Rmin_beforeAbs && v0.v0radius() < Rmax_beforeAbs)
          registryData.fill(HIST("AntiLambda_before_target_data"), v0.p(), v0.mAntiLambda());

        // After Target
        if (v0.v0radius() > Rmin_afterAbs && v0.v0radius() < Rmax_afterAbs)
          registryData.fill(HIST("AntiLambda_after_target_data"), v0.p(), v0.mAntiLambda());
      }

    } // end loop on V0s
  }   // end processData
  PROCESS_SWITCH(vzero_cascade_absorption, processData, "Process data", true);

  // Process MC
  void processMC(soa::Join<SelectedCollisions, aod::McCollisionLabels>::iterator const& collision, aod::V0Datas const& fullV0s, MCTracks const& tracks, aod::McParticles& mcParticles, aod::McCollisions const& mcCollisions)
  {

    // Event Counter (before event sel)
    registryQC.fill(HIST("event_counter"), 2.5);

    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter (after event sel)
    registryQC.fill(HIST("event_counter"), 3.5);

    // Loop over Reconstructed V0s
    for (auto& v0 : fullV0s) {

      // Positive and Negative Tracks
      const auto& posTrack = v0.posTrack_as<MCTracks>();
      const auto& negTrack = v0.negTrack_as<MCTracks>();

      // Require TPC Refit
      if (!posTrack.passedTPCRefit())
        continue;
      if (!negTrack.passedTPCRefit())
        continue;

      auto hit_ITS_before_target = static_cast<std::vector<float>>(hit_requirement_before_target);
      auto hit_ITS_after_target = static_cast<std::vector<float>>(hit_requirement_after_target);

      // ITS Requirement
      bool satisfyITSreq = true;
      if (requirehitsITS) {
        for (int i = 0; i < 7; i++) {
          if (hit_ITS_before_target[i] > 0 && !hasHitOnITSlayer(posTrack.itsClusterMap(), i)) {
            satisfyITSreq = false;
            break;
          }
          if (hit_ITS_after_target[i] > 0 && !hasHitOnITSlayer(negTrack.itsClusterMap(), i)) {
            satisfyITSreq = false;
            break;
          }
        }
      }

      if (requirehitsITS && (!satisfyITSreq))
        continue;

      // MC Particles
      if (!posTrack.has_mcParticle())
        continue;
      if (!negTrack.has_mcParticle())
        continue;

      auto posParticle = posTrack.mcParticle_as<aod::McParticles>();
      auto negParticle = negTrack.mcParticle_as<aod::McParticles>();
      if (!posParticle.has_mothers() || !negParticle.has_mothers()) {
        continue;
      }

      bool isK0s = false;
      bool isLambda = false;
      bool isAntiLambda = false;

      for (auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
        for (auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 310)
            isK0s = true;

          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == +3122)
            isLambda = true;
          if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == -3122)
            isAntiLambda = true;
        }
      }

      // Resolution in radial position
      float Rgen = TMath::Sqrt(posParticle.vx() * posParticle.vx() + posParticle.vy() * posParticle.vy());
      float Rrec = v0.v0radius();
      float deltaR = Rgen - Rrec;

      // K0 Short
      if (passedK0Selection(v0, negTrack, posTrack, collision) && isK0s) {

        // Before Target
        if (v0.v0radius() > Rmin_beforeAbs && v0.v0radius() < Rmax_beforeAbs) {
          registryMC.fill(HIST("K0_before_target_mc"), v0.p(), v0.mK0Short());
          registryMC.fill(HIST("K0_Rresolution_before_target"), v0.p(), deltaR);
        }

        // After Target
        if (v0.v0radius() > Rmin_afterAbs && v0.v0radius() < Rmax_afterAbs) {
          registryMC.fill(HIST("K0_after_target_mc"), v0.p(), v0.mK0Short());
          registryMC.fill(HIST("K0_Rresolution_after_target"), v0.p(), deltaR);
        }
      }

      // Lambda
      if (passedLambdaSelection(v0, negTrack, posTrack, collision) && isLambda) {

        // Before Target
        if (v0.v0radius() > Rmin_beforeAbs && v0.v0radius() < Rmax_beforeAbs) {
          registryMC.fill(HIST("Lambda_before_target_mc"), v0.p(), v0.mLambda());
          registryMC.fill(HIST("Lambda_Rresolution_before_target"), v0.p(), deltaR);
        }

        // After Target
        if (v0.v0radius() > Rmin_afterAbs && v0.v0radius() < Rmax_afterAbs) {
          registryMC.fill(HIST("Lambda_after_target_mc"), v0.p(), v0.mLambda());
          registryMC.fill(HIST("Lambda_Rresolution_after_target"), v0.p(), deltaR);
        }
      }

      // AntiLambda
      if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision) && isAntiLambda) {

        // Before Target
        if (v0.v0radius() > Rmin_beforeAbs && v0.v0radius() < Rmax_beforeAbs) {
          registryMC.fill(HIST("AntiLambda_before_target_mc"), v0.p(), v0.mAntiLambda());
          registryMC.fill(HIST("AntiLambda_Rresolution_before_target"), v0.p(), deltaR);
        }

        // After Target
        if (v0.v0radius() > Rmin_afterAbs && v0.v0radius() < Rmax_afterAbs) {
          registryMC.fill(HIST("AntiLambda_after_target_mc"), v0.p(), v0.mAntiLambda());
          registryMC.fill(HIST("AntiLambda_Rresolution_after_target"), v0.p(), deltaR);
        }
      }

    } // end loop on V0s
  }   // end processMC
  PROCESS_SWITCH(vzero_cascade_absorption, processMC, "Process mc", false);

  void processMCgen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles)
  {
    for (auto& mcParticle : mcParticles) {

      if (mcParticle.eta() < etaMin || mcParticle.eta() > etaMax)
        continue;

      float R = TMath::Sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy());
      float p = mcParticle.p();

      // Lambda
      if (mcParticle.pdgCode() == +3122) {
        if (R > Rmin_beforeAbs && R < Rmax_beforeAbs)
          registryMC.fill(HIST("Lambda_before_target_mc_gen"), p);
        if (R > Rmin_afterAbs && R < Rmax_afterAbs)
          registryMC.fill(HIST("Lambda_after_target_mc_gen"), p);
      }

      // AntiLambda
      if (mcParticle.pdgCode() == -3122) {
        if (R > Rmin_beforeAbs && R < Rmax_beforeAbs)
          registryMC.fill(HIST("AntiLambda_before_target_mc_gen"), p);
        if (R > Rmin_afterAbs && R < Rmax_afterAbs)
          registryMC.fill(HIST("AntiLambda_after_target_mc_gen"), p);
      }
    }
  }
  PROCESS_SWITCH(vzero_cascade_absorption, processMCgen, "Process generated MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vzero_cascade_absorption>(cfgc)};
}
