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
/// \author Dushmanta Sahu (Dushmanta.Sahu@cern.ch)
/// \since October 8, 2025
/// \file multiplicitypt.cxx
/// \brief  Analysis to do PID

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TH1.h>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace constants::physics;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection,
  aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFbeta>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;

struct MultiplicityPt {

  // dE/dx for all charged particles
  HistogramRegistry registryDeDx{
    "registryDeDx",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Configurable Parameters
  // Tracks cuts
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 70.0f,
                                      "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of found TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f,
                                 "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f,
                                 "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8f, "etaMax"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxDCAz{"maxDCAz", 2.f, "maxDCAz"};
  // v0 cuts
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
  Configurable<bool> calibrationMode{"calibrationMode", false, "calibration mode"};
  // pT spectra and multiplicity parameters
  Configurable<float> minPtCut{"minPtCut", 0.15f, "Minimum pT cut for primary particles"};
  Configurable<float> maxPtCut{"maxPtCut", 10.0f, "Maximum pT cut for primary particles"};
  Configurable<bool> enablePrimarySpectra{"enablePrimarySpectra", true, "Enable primary particle spectra analysis"};
  // Normalization parameters
  Configurable<float> rapidityRange{"rapidityRange", 0.8f, "Rapidity range for normalization (|y| < y_max)"};
  Configurable<bool> enableNormalization{"enableNormalization", true, "Enable pT spectrum normalization"};

  // Multiplicity classes
  Configurable<std::vector<double>> multiplicityClasses{"multiplicityClasses",
                                                        {0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0},
                                                        "Multiplicity class boundaries in percent"};

  // Histograms names
  static constexpr std::string_view kDedxvsMomentumPos[4] = {"dEdx_vs_Momentum_all_Pos", "dEdx_vs_Momentum_Pi_v0_Pos", "dEdx_vs_Momentum_Pr_v0_Pos", "dEdx_vs_Momentum_El_v0_Pos"};
  static constexpr std::string_view kDedxvsMomentumNeg[4] = {"dEdx_vs_Momentum_all_Neg", "dEdx_vs_Momentum_Pi_v0_Neg", "dEdx_vs_Momentum_Pr_v0_Neg", "dEdx_vs_Momentum_El_v0_Neg"};
  static constexpr double EtaCut[9] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  Configurable<std::vector<float>> calibrationFactorNeg{"calibrationFactorNeg", {50.4011, 50.4764, 50.186, 49.2955, 48.8222, 49.4273, 49.9292, 50.0556}, "negative calibration factors"};
  Configurable<std::vector<float>> calibrationFactorPos{"calibrationFactorPos", {50.5157, 50.6359, 50.3198, 49.3345, 48.9197, 49.4931, 50.0188, 50.1406}, "positive calibration factors"};
  ConfigurableAxis binP{"binP", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, ""};

  // Store multiplicity class names and histogram references
  std::vector<std::string> multClassNames;
  std::vector<std::string> multHistNames;
  int nClasses; // Number of multiplicity classes

  TrackSelection myTrackSelection()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(etaMin, etaMax);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    selectedTracks.SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC);
    selectedTracks.SetMaxChi2PerClusterTPC(maxChi2TPC);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1});
    selectedTracks.SetMaxChi2PerClusterITS(maxChi2ITS);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / std::pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(maxDCAz);

    return selectedTracks;
  }

  TrackSelection mySelectionPrim;

  void init(InitContext const&)
  {
    AxisSpec dedxAxis{100, 0.0, 100.0, "dE/dx MIP (a. u.)"};
    AxisSpec etaAxis{8, -0.8, 0.8, "#eta"};
    AxisSpec pAxis = {binP, "#it{p}/Z (GeV/c)"};

    // pT spectra and multiplicity histograms
    AxisSpec ptAxis{100, 0.0, 10.0, "#it{p}_{T} (GeV/c)"};
    AxisSpec multiplicityAxis{200, -0.5, 199.5, "N_{tracks}"};
    AxisSpec centralityAxis{100, 0.0, 100.0, "Centrality (%)"};
    AxisSpec etaFullAxis{50, -2.5, 2.5, "#eta"};
    AxisSpec phiAxis{100, 0.0, 6.283, "#phi"};

    // Add pT spectra and multiplicity histograms
    registryDeDx.add("hPtSpectrum", "pT spectrum of primary charged particles", HistType::kTH1F, {ptAxis});
    registryDeDx.add("hEtaDistribution", "eta distribution", HistType::kTH1F, {etaFullAxis});
    registryDeDx.add("hPhiDistribution", "phi distribution", HistType::kTH1F, {phiAxis});
    registryDeDx.add("hMultiplicity", "Track multiplicity", HistType::kTH1F, {multiplicityAxis});
    registryDeDx.add("hMultiplicityVsCentrality", "Multiplicity vs centrality", HistType::kTH2F, {centralityAxis, multiplicityAxis});
    registryDeDx.add("hPtVsEta", "pT vs eta", HistType::kTH2F, {etaFullAxis, ptAxis});

    // Add event counter histograms for normalization
    registryDeDx.add("hEventCounter", "Event counter", HistType::kTH1D, {{3, 0.5, 3.5}});
    registryDeDx.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(1, "All Events");
    registryDeDx.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(2, "Events with Tracks");
    registryDeDx.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(3, "Tracks in |y|<y_max");

    // Create multiplicity class names and histograms
    nClasses = multiplicityClasses->size() - 1;
    for (int i = 0; i < nClasses; ++i) {
      std::string className = Form("%.0f-%.0f%%", multiplicityClasses->at(i), multiplicityClasses->at(i + 1));
      multClassNames.push_back(className);

      std::string histName = Form("hPtSpectrum_mult%d", static_cast<int>(i));
      multHistNames.push_back(histName);

      // Create pT spectra for each multiplicity class
      registryDeDx.add(histName.c_str(),
                       Form("pT spectrum for multiplicity %s", className.c_str()),
                       HistType::kTH1F, {ptAxis});
    }

    // ADDED: Event counters for each multiplicity class
    registryDeDx.add("hEventsPerMultClass", "Number of events per multiplicity class",
                     HistType::kTH1D, {{nClasses, 0.5, nClasses + 0.5}});

    // ADDED: Event counters with tracks for each multiplicity class
    registryDeDx.add("hEventsWithTracksPerMultClass", "Number of events with tracks per multiplicity class",
                     HistType::kTH1D, {{nClasses, 0.5, nClasses + 0.5}});

    // Additional histogram for multiplicity distribution
    registryDeDx.add("hMultiplicityClasses", "Multiplicity class distribution",
                     HistType::kTH1D, {{nClasses, 0.5, nClasses + 0.5}});

    // Set bin labels for multiplicity classes
    auto hMultiplicityClasses = registryDeDx.get<TH1>(HIST("hMultiplicityClasses"));
    auto hEventsPerMultClass = registryDeDx.get<TH1>(HIST("hEventsPerMultClass"));
    auto hEventsWithTracksPerMultClass = registryDeDx.get<TH1>(HIST("hEventsWithTracksPerMultClass"));

    for (int i = 0; i < nClasses; ++i) {
      hMultiplicityClasses->GetXaxis()->SetBinLabel(i + 1, multClassNames[i].c_str());
      hEventsPerMultClass->GetXaxis()->SetBinLabel(i + 1, multClassNames[i].c_str());
      hEventsWithTracksPerMultClass->GetXaxis()->SetBinLabel(i + 1, multClassNames[i].c_str());
    }

    if (calibrationMode) {
      // MIP for pions
      registryDeDx.add(
        "hdEdx_vs_eta_Neg_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_Pos_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});
      // MIP for electrons
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      // Pions from TOF
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

    } else {
      // MIP for pions
      registryDeDx.add(
        "hdEdx_vs_eta_Neg_calibrated_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_Pos_calibrated_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});

      // MIP for electrons
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_calibrated_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_calibrated_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      // Pions from TOF
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_calibrated_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_calibrated_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      // De/Dx for ch and v0 particles
      for (int i = 0; i < 4; ++i) {
        registryDeDx.add(kDedxvsMomentumPos[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxis}, {dedxAxis}, {etaAxis}});
        registryDeDx.add(kDedxvsMomentumNeg[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxis}, {dedxAxis}, {etaAxis}});
      }
    }

    registryDeDx.add(
      "hdEdx_vs_phi", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 6.4, "#phi"}, {dedxAxis}});

    registryDeDx.add(
      "hbeta_vs_p_Neg", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    registryDeDx.add(
      "hbeta_vs_p_Pos", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});
    // Event Counter
    registryDeDx.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{100, -20.0, +20.0, "z_{vtx} (cm)"}});

    mySelectionPrim = myTrackSelection();
  }

  // Function to determine multiplicity class
  int getMultiplicityClass(float centrality)
  {
    for (size_t i = 0; i < multiplicityClasses->size() - 1; ++i) {
      if (centrality >= multiplicityClasses->at(i) && centrality < multiplicityClasses->at(i + 1)) {
        return i;
      }
    }
    return -1; // Outside defined classes
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
      if (std::abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
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
      if (std::abs(ptrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
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
      if (std::abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPr()) > nsigmaTOFmax)
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
      if (std::abs(ptrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaEl()) > nsigmaTOFmax)
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
    registryDeDx.fill(HIST("hEventCounter"), 1); // Count all events

    // Centrality
    float centrality = collision.centFT0C();
    if (centrality < 0.0 || centrality > 100.0)
      centrality = 1.0;

    // Determine multiplicity class
    int multClass = getMultiplicityClass(centrality);

    // Count event in multiplicity class (regardless of track selection)
    if (multClass >= 0) {
      registryDeDx.fill(HIST("hMultiplicityClasses"), multClass + 1);
      registryDeDx.fill(HIST("hEventsPerMultClass"), multClass + 1); // Count ALL events in this class
    }

    // Multiplicity counter for primary charged particles
    int nPrimaryTracks = 0;
    int nTracksInRapidity = 0;
    bool hasGoodTracks = false;

    for (const auto& trk : tracks) {

      // track Selection
      if (!passedSingleTrackSelection(trk, collision))
        continue;

      if (!mySelectionPrim.IsSelected(trk))
        continue;

      hasGoodTracks = true;

      // Fill basic track distributions for primary particles
      registryDeDx.fill(HIST("hEtaDistribution"), trk.eta());
      registryDeDx.fill(HIST("hPhiDistribution"), trk.phi());
      registryDeDx.fill(HIST("hPtVsEta"), trk.eta(), trk.pt());

      // Fill pT spectrum for primary charged particles
      if (enablePrimarySpectra && trk.pt() >= minPtCut && trk.pt() <= maxPtCut) {
        registryDeDx.fill(HIST("hPtSpectrum"), trk.pt());
        nPrimaryTracks++;

        // Fill pT spectrum for specific multiplicity class
        if (multClass >= 0 && multClass < nClasses) {
          switch (multClass) {
            case 0:
              registryDeDx.fill(HIST("hPtSpectrum_mult0"), trk.pt());
              break;
            case 1:
              registryDeDx.fill(HIST("hPtSpectrum_mult1"), trk.pt());
              break;
            case 2:
              registryDeDx.fill(HIST("hPtSpectrum_mult2"), trk.pt());
              break;
            case 3:
              registryDeDx.fill(HIST("hPtSpectrum_mult3"), trk.pt());
              break;
            case 4:
              registryDeDx.fill(HIST("hPtSpectrum_mult4"), trk.pt());
              break;
            case 5:
              registryDeDx.fill(HIST("hPtSpectrum_mult5"), trk.pt());
              break;
            case 6:
              registryDeDx.fill(HIST("hPtSpectrum_mult6"), trk.pt());
              break;
            case 7:
              registryDeDx.fill(HIST("hPtSpectrum_mult7"), trk.pt());
              break;
            default:
              break;
          }
        }

        // Count tracks in rapidity range for normalization
        if (enableNormalization && std::abs(trk.eta()) < rapidityRange) {
          nTracksInRapidity++;
        }
      }

      float signedP = trk.sign() * trk.tpcInnerParam();

      // MIP calibration for pions
      if (trk.tpcInnerParam() >= 0.35 && trk.tpcInnerParam() <= 0.45) {
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_Pi"), trk.eta(), trk.tpcSignal());
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_Pi"), trk.eta(), trk.tpcSignal());
          }

        } else {
          for (int i = 0; i < 8; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_calibrated_Pi"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_calibrated_Pi"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i));
              }
            }
          }
        }
      }
      // Beta from TOF
      if (signedP < 0) {
        registryDeDx.fill(HIST("hbeta_vs_p_Neg"), std::abs(signedP), trk.beta());
      } else {
        registryDeDx.fill(HIST("hbeta_vs_p_Pos"), signedP, trk.beta());
      }
      // Electrons from TOF
      if (std::abs(trk.beta() - 1) < 0.1) { // beta cut
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_El"), trk.eta(), trk.tpcSignal(), std::abs(signedP));
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_El"), trk.eta(), trk.tpcSignal(), signedP);
          }
        } else {
          for (int i = 0; i < 8; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_calibrated_El"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), std::abs(signedP));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_calibrated_El"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i), signedP);
              }
            }
          }
        }
      }
      // pions from TOF
      if (trk.beta() > 1. && trk.beta() < 1.05) { // beta cut
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_TOF"), trk.eta(), trk.tpcSignal(), std::abs(signedP));
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_TOF"), trk.eta(), trk.tpcSignal(), signedP);
          }
        } else {
          for (int i = 0; i < 8; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_calibrated_TOF"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), std::abs(signedP));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_calibrated_TOF"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i), signedP);
              }
            }
          }
        }
      }

      registryDeDx.fill(HIST("hdEdx_vs_phi"), trk.phi(), trk.tpcSignal());

      if (!calibrationMode) {
        for (int i = 0; i < 8; ++i) {
          if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
            if (signedP > 0) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[0]), signedP, trk.tpcSignal() * 50 / calibrationFactorPos->at(i), trk.eta());
            } else {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[0]), std::abs(signedP), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), trk.eta());
            }
          }
        }
      }
    }

    // Count events with tracks for normalization
    if (hasGoodTracks) {
      registryDeDx.fill(HIST("hEventCounter"), 2); // Events with tracks

      // ADDED: Count events with tracks per multiplicity class
      if (multClass >= 0) {
        registryDeDx.fill(HIST("hEventsWithTracksPerMultClass"), multClass + 1);
      }

      if (nTracksInRapidity > 0) {
        registryDeDx.fill(HIST("hEventCounter"), 3, nTracksInRapidity); // Tracks in rapidity range
      }
    }

    // Fill multiplicity histograms
    registryDeDx.fill(HIST("hMultiplicity"), nPrimaryTracks);
    registryDeDx.fill(HIST("hMultiplicityVsCentrality"), centrality, nPrimaryTracks);

    // Loop over Reconstructed V0s
    if (!calibrationMode) {
      for (const auto& v0 : fullV0s) {

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

        float pxPos = posTrack.px();
        float pyPos = posTrack.py();
        float pzPos = posTrack.pz();

        float pxNeg = negTrack.px();
        float pyNeg = negTrack.py();
        float pzNeg = negTrack.pz();

        const float gammaMass = 2 * MassElectron; // GeV/c^2

        // K0s Selection
        if (passedK0Selection(v0, negTrack, posTrack, collision)) {
          float ePosPi = posTrack.energy(MassPionCharged);
          float eNegPi = negTrack.energy(MassPionCharged);

          float invMass = std::sqrt((eNegPi + ePosPi) * (eNegPi + ePosPi) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassK0Short) > 0.01) {
            continue;
          }

          for (int i = 0; i < 8; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Lambda Selection
        if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {

          float ePosPr = posTrack.energy(MassProton);
          float eNegPi = negTrack.energy(MassPionCharged);

          float invMass = std::sqrt((eNegPi + ePosPr) * (eNegPi + ePosPr) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassLambda) > 0.01) {
            continue;
          }

          for (int i = 0; i < 8; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[2]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // AntiLambda Selection
        if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {

          float ePosPi = posTrack.energy(MassPionCharged);
          float eNegPr = negTrack.energy(MassProton);

          float invMass = std::sqrt((eNegPr + ePosPi) * (eNegPr + ePosPi) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassLambda) > 0.01) {
            continue;
          }

          for (int i = 0; i < 8; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[2]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Gamma Selection
        if (passedGammaSelection(v0, negTrack, posTrack, collision)) {

          float ePosEl = posTrack.energy(MassElectron);
          float eNegEl = negTrack.energy(MassElectron);

          float invMass = std::sqrt((eNegEl + ePosEl) * (eNegEl + ePosEl) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - gammaMass) > 0.0015) {
            continue;
          }

          for (int i = 0; i < 8; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[3]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[3]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityPt>(cfgc)};
}
