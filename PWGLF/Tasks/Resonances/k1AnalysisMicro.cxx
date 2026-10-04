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
/// \file k1AnalysisMicro.cxx
/// \brief Reconstruction of track-track decay resonance candidates
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#include "PWGLF/Core/K1AnalysisMicroCore.h"
#include "PWGLF/Core/ResoAnalysisSelectionCore.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>
#include <cmath>
#include <tuple>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::analysis::resonance;
using namespace o2::analysis::k1micro;

/// Histogram binning, QA and debug output
struct HistogramOptions : ConfigurableGroup {
  Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<int> cfgTruthDebug{"cfgTruthDebug", 0, "Maximum logged matched candidates per truth channel"};
};

enum BinAnti : unsigned int {
  kNormal = 0,
  kAnti,
  kNAEnd
};

enum BinType : unsigned int {
  kK1P = 0,
  kK1N,
  kK1P_Mix,
  kK1N_Mix,
  kK1P_GenINEL10,
  kK1N_GenINEL10,
  kK1P_GenINELgt10,
  kK1N_GenINELgt10,
  kK1P_GenTrig10,
  kK1N_GenTrig10,
  kK1P_GenEvtSel,
  kK1N_GenEvtSel,
  kK1P_Rec,
  kK1N_Rec,
  kTYEnd
};

enum class QAFolder {
  Before, // QA/*: before the candidate cuts
  After,  // QAcut/*: after the candidate cuts
  MC      // QAMC/*: matched K1 truth candidates
};

struct K1AnalysisMicro {
  // Module-initializer v001 tables; full tracks keep their unversioned schema as a fallback.
  using ResoCollisions = aod::ResoCollisions_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;
  using ResoTracks = aod::ResoTracks; // no v001 exists; K1 does not need ResoTrackTracks (trackId unused)
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCTracks = soa::Join<ResoTracks, aod::ResoMCTracks>;
  using ResoMCMicroTracks = soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001>;
  using ResoMCParents = aod::ResoMCParents_001;

  SliceCache cache;
  // Registered only to enable the slice cache that SameKindPair (event mixing) needs, as in Xi1820Analysis
  Preslice<ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<ResoMicroTracks> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Selection shared with the K1 training-table task (plain JSON keys, no group prefix)
  EventCuts eventCuts;
  TrackCuts trackCuts;
  PionPidCuts pionPID;
  KaonPidCuts kaonPID;
  SecondaryCuts secondaryCuts;
  CandidateCuts candidateCuts;
  HistogramOptions histogramOptions;

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};

  K1AnalysisMicroCore core;
  std::array<int, NTruthChannels> truthDebugCounts{};

  void init(InitContext&)
  {
    const int sameEventModes = static_cast<int>(doprocessResoTracks) + static_cast<int>(doprocessResoMicroTracks) +
                               static_cast<int>(doprocessMC) + static_cast<int>(doprocessMCMicro);
    const int mixedEventModes = static_cast<int>(doprocessME) + static_cast<int>(doprocessMEMicro);
    if (sameEventModes > 1 || mixedEventModes > 1) {
      LOG(fatal) << "Enable at most one same-event mode and one mixing mode";
    }

    ProcessModes modes;
    modes.microTracks = doprocessResoMicroTracks || doprocessMCMicro || doprocessMEMicro;
    modes.mcReco = doprocessMC || doprocessMCMicro;
    modes.mcRecoMicro = doprocessMCMicro;
    modes.mcGen = doprocessMCTrue;
    core.init(histos, eventCuts, trackCuts, pionPID, kaonPID, secondaryCuts, candidateCuts, modes);
    registerHistograms(modes);

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in K1 Analysis Task";
    histos.print();
  }

  void registerHistograms(ProcessModes const& modes)
  {
    const int nBinsDiv = histogramOptions.cNbinsDiv;
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 200.};
    AxisSpec centAxis = {centBinning, "T0M (%)"};
    AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaxyAxis = {300, 0, 3, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {500, 0, 5, "DCA_{#it{z}} (cm)"};
    AxisSpec invMassAxisK892 = {1400 / nBinsDiv, 0.6, 2.0, "Invariant Mass (GeV/#it{c}^2)"};   // K(892)0
    AxisSpec invMassAxisRho = {2000 / nBinsDiv, 0.0, 2.0, "Invariant Mass (GeV/#it{c}^2)"};    // rho
    AxisSpec invMassAxisReso = {1600 / nBinsDiv, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"}; // K1
    AxisSpec pidQAAxis = {130, -6.5, 6.5};

    // THnSparse
    AxisSpec axisAnti = {BinAnti::kNAEnd, 0, BinAnti::kNAEnd, "Type of bin: Normal or Anti"};
    AxisSpec axisType = {BinType::kTYEnd, 0, BinType::kTYEnd, "Type of bin with charge and mix"};

    // DCA QA
    // Primary pion
    histos.add("QA/trkppionDCAxy", "DCAxy disstribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkppionDCAz", "DCAz disstribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkppionDCAxy", "DCAxy distribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkppionDCAz", "DCAz distribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // Secondary pion
    histos.add("QA/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // Kaon
    histos.add("QA/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // K1
    histos.add("QA/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
    histos.add("QA/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
    histos.add("QA/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
    histos.add("QA/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
    histos.add("QA/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});

    histos.add("QAcut/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
    histos.add("QAcut/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QAcut/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
    histos.add("QAcut/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
    histos.add("QAcut/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
    histos.add("QAcut/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});

    // Invariant mass
    histos.add("hInvmass_K1", "Invariant mass of K1(1270) (US)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K1_LS", "Invariant mass of K1(1270) (LS)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K1_Mix", "Invariant mass of K1(1270) (ME)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    // Mass QA (quick check)
    histos.add("k1invmass", "Invariant mass of K1(1270) (US)", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_LS", "Invariant mass of K1(1270) (LS)", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_Mix", "Invariant mass of K1(1270) (ME)", HistType::kTH1F, {invMassAxisReso});

    // MC
    if (modes.mcReco) {
      AxisSpec channelAxis = {3, -0.5, 2.5, "0: non-K1, 1: rho K, 2: K* pi"};
      histos.add("MCReco/collisions", "Selected reconstructed MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/microTracks", "Input micro tracks in selected MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/channel", "All selected pi-pi-K combinations by truth channel", HistType::kTH1D, {channelAxis});
      histos.add("MCReco/mass", "Reconstructed mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisReso});
      histos.add("MCReco/pt", "Reconstructed pT by truth channel", HistType::kTH2D, {channelAxis, ptAxis});
      histos.add("MCReco/piPiMass", "pi-pi mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisRho});
      histos.add("MCReco/pi1KMass", "First pion-kaon mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisK892});
      histos.add("MCReco/pi2KMass", "Second pion-kaon mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisK892});
      histos.add("k1invmass_MC", "Invariant mass of K1(1270)", HistType::kTH1F, {invMassAxisReso});
      histos.add("k1invmass_MC_noK1", "Invariant mass of K1(1270)", HistType::kTH1F, {invMassAxisReso});

      histos.add("QAMC/trkppionDCAxy", "DCAxy distribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkppionDCAz", "DCAz distribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
      histos.add("QAMC/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
      histos.add("QAMC/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
      histos.add("QAMC/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
      histos.add("QAMC/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
      histos.add("QAMC/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});
    } // mcReco
    if (modes.mcGen) {
      AxisSpec channelAxis = {3, -0.5, 2.5, "0: other/unresolved, 1: rho K, 2: K* pi"};
      histos.add("MCGen/chargeChannel", "K1 parents in selected reconstructed events, inside the K1 rapidity window", HistType::kTH2D, {{2, -1.5, 1.5, "K1 charge"}, channelAxis});
      histos.add("MCGen/ptChannel", "Generated K1 pT by immediate decay channel", HistType::kTH2D, {channelAxis, ptAxis});
    }
  }

  // Track QA of a pion; isPrimary selects the trkppion (first) or trkspion (second) histograms
  template <QAFolder Folder, typename TrackType>
  void fillPionQA(const TrackType& track, bool isPrimary)
  {
    const bool hasTOF = track.hasTOF();
    if (isPrimary) {
      if constexpr (Folder == QAFolder::Before) {
        histos.fill(HIST("QA/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QA/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QA/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QA/trkppionpT"), track.pt());
        histos.fill(HIST("QA/trkppionDCAxy"), track.dcaXY());
        histos.fill(HIST("QA/trkppionDCAz"), track.dcaZ());
      } else if constexpr (Folder == QAFolder::After) {
        histos.fill(HIST("QAcut/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAcut/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAcut/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAcut/trkppionpT"), track.pt());
        histos.fill(HIST("QAcut/trkppionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAcut/trkppionDCAz"), track.dcaZ());
      } else {
        histos.fill(HIST("QAMC/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAMC/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAMC/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAMC/trkppionpT"), track.pt());
        histos.fill(HIST("QAMC/trkppionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAMC/trkppionDCAz"), track.dcaZ());
      }
    } else {
      if constexpr (Folder == QAFolder::Before) {
        histos.fill(HIST("QA/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QA/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QA/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QA/trkspionpT"), track.pt());
        histos.fill(HIST("QA/trkspionDCAxy"), track.dcaXY());
        histos.fill(HIST("QA/trkspionDCAz"), track.dcaZ());
      } else if constexpr (Folder == QAFolder::After) {
        histos.fill(HIST("QAcut/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAcut/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAcut/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAcut/trkspionpT"), track.pt());
        histos.fill(HIST("QAcut/trkspionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAcut/trkspionDCAz"), track.dcaZ());
      } else {
        histos.fill(HIST("QAMC/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAMC/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAMC/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAMC/trkspionpT"), track.pt());
        histos.fill(HIST("QAMC/trkspionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAMC/trkspionDCAz"), track.dcaZ());
      }
    }
  }

  // Track QA of the bachelor kaon
  template <QAFolder Folder, typename TrackType>
  void fillKaonQA(const TrackType& track)
  {
    const bool hasTOF = track.hasTOF();
    if constexpr (Folder == QAFolder::Before) {
      histos.fill(HIST("QA/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("QA/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QA/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
      histos.fill(HIST("QA/trkkaonpT"), track.pt());
      histos.fill(HIST("QA/trkkaonDCAxy"), track.dcaXY());
      histos.fill(HIST("QA/trkkaonDCAz"), track.dcaZ());
    } else if constexpr (Folder == QAFolder::After) {
      histos.fill(HIST("QAcut/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("QAcut/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QAcut/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
      histos.fill(HIST("QAcut/trkkaonpT"), track.pt());
      histos.fill(HIST("QAcut/trkkaonDCAxy"), track.dcaXY());
      histos.fill(HIST("QAcut/trkkaonDCAz"), track.dcaZ());
    } else {
      histos.fill(HIST("QAMC/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("QAMC/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QAMC/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
      histos.fill(HIST("QAMC/trkkaonpT"), track.pt());
      histos.fill(HIST("QAMC/trkkaonDCAxy"), track.dcaXY());
      histos.fill(HIST("QAMC/trkkaonDCAz"), track.dcaZ());
    }
  }

  // Histograms of the selected pion pairs and (pion, pion, kaon) candidates of one collision
  // (or one mixed pair of collisions). The selection itself is the shared K1 core.
  // dTracks1: bachelor kaons, dTracks2: pions.
  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    const bool fillQA = !IsMix && histogramOptions.additionalQAplots;
    const auto multiplicity = collision.cent();

    // Pion pair passing the pion selection; trk1 is the pion with the lower index
    auto onPair = [&](auto const& trk1, auto const& trk2, ROOT::Math::PxPyPzMVector const& secondary, bool passesPairPt) {
      if (fillQA) {
        fillPionQA<QAFolder::Before>(trk1, true);
        fillPionQA<QAFolder::Before>(trk2, false);
      }
      if (!passesPairPt) {
        return;
      }
      if (fillQA) {
        histos.fill(HIST("QA/hInvmassSecon"), secondary.M());
      }
      if constexpr (IsMC) {
        histos.fill(HIST("QAMC/hpT_Secondary"), secondary.Pt());
      }
    };

    // Candidate passing the pion, pair and kaon selection; pion1 is the K*0 partner in the unlike-sign case
    auto onCandidate = [&](auto const& kaon, auto const& pion1, auto const& pion2, K1CandidateValues const& c) {
      if (fillQA) {
        fillKaonQA<QAFolder::Before>(kaon);
      }
      if (!c.inRapidity) {
        return;
      }

      // QA histogram before the candidate cuts
      if (fillQA) {
        histos.fill(HIST("QA/K1OA"), c.angle);
        histos.fill(HIST("QA/K1PairAsym"), c.pairAsym);
        histos.fill(HIST("QA/hInvmassK892_Rho"), c.mass13, c.secondary.M());
        histos.fill(HIST("QA/hInvmassSecon_PiKa"), c.secondary.M(), c.mass23);
        histos.fill(HIST("QA/hpT_Secondary"), c.secondary.Pt());
      }
      if (!c.passesCandidateCuts) {
        return;
      }

      // QA histograms after the candidate cuts
      if (fillQA) {
        fillPionQA<QAFolder::After>(pion1, true);
        fillPionQA<QAFolder::After>(pion2, false);
        fillKaonQA<QAFolder::After>(kaon);
        histos.fill(HIST("QAcut/K1OA"), c.angle);
        histos.fill(HIST("QAcut/K1PairAsym"), c.pairAsym);
        histos.fill(HIST("QAcut/hInvmassK892_Rho"), c.mass13, c.secondary.M());
        histos.fill(HIST("QAcut/hInvmassSecon_PiKa"), c.secondary.M(), c.mass23);
        histos.fill(HIST("QAcut/hInvmassSecon"), c.secondary.M());
        histos.fill(HIST("QAcut/hpT_Secondary"), c.secondary.Pt());
      }

      const unsigned int typeNormal = BinAnti::kNormal;
      if constexpr (IsMix) {
        const unsigned int typeK1 = kaon.sign() > 0 ? BinType::kK1P_Mix : BinType::kK1N_Mix;
        histos.fill(HIST("hInvmass_K1_Mix"), typeNormal, typeK1, multiplicity, c.k1.Pt(), c.k1.M());
        histos.fill(HIST("k1invmass_Mix"), c.k1.M());
        return;
      }
      unsigned int typeK1 = kaon.sign() > 0 ? BinType::kK1P : BinType::kK1N;
      if (c.isUnlikeSign) {
        histos.fill(HIST("k1invmass"), c.k1.M());
        histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, c.k1.Pt(), c.k1.M());
      } else {
        histos.fill(HIST("k1invmass_LS"), c.k1.M());
        histos.fill(HIST("hInvmass_K1_LS"), typeNormal, typeK1, multiplicity, c.k1.Pt(), c.k1.M());
      }

      if constexpr (IsMC) {
        const auto channel = classifyK1Truth(pion1, pion2, kaon);
        const int channelBin = static_cast<int>(channel);
        histos.fill(HIST("MCReco/channel"), channelBin);
        histos.fill(HIST("MCReco/mass"), channelBin, c.k1.M());
        histos.fill(HIST("MCReco/pt"), channelBin, c.k1.Pt());
        histos.fill(HIST("MCReco/piPiMass"), channelBin, c.secondary.M());
        histos.fill(HIST("MCReco/pi1KMass"), channelBin, c.mass13);
        histos.fill(HIST("MCReco/pi2KMass"), channelBin, c.mass23);
        if (channel == K1TruthChannel::None) {
          histos.fill(HIST("k1invmass_MC_noK1"), c.k1.M());
          return;
        }
        if (truthDebugCounts[channelBin] < histogramOptions.cfgTruthDebug) {
          ++truthDebugCounts[channelBin];
          LOGP(info, "K1Truth channel={} collision={} tracks=({},{},{}) pdg=({},{},{}) mothers=({},{},{}) motherPDG=({},{},{}) siblings=(({},{}),({},{}),({},{}))",
               channelBin, collision.globalIndex(),
               pion1.globalIndex(), pion2.globalIndex(), kaon.globalIndex(),
               pion1.pdgCode(), pion2.pdgCode(), kaon.pdgCode(), pion1.motherId(), pion2.motherId(), kaon.motherId(),
               pion1.motherPDG(), pion2.motherPDG(), kaon.motherPDG(),
               pion1.siblingIds()[0], pion1.siblingIds()[1], pion2.siblingIds()[0], pion2.siblingIds()[1], kaon.siblingIds()[0], kaon.siblingIds()[1]);
        }
        typeK1 = kaon.sign() > 0 ? BinType::kK1P_Rec : BinType::kK1N_Rec;
        histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, c.k1.Pt(), c.k1.M());
        histos.fill(HIST("k1invmass_MC"), c.k1.M());
        histos.fill(HIST("QAMC/K1OA"), c.angle);
        histos.fill(HIST("QAMC/K1PairAsym"), c.pairAsym);
        histos.fill(HIST("QAMC/hInvmassK892_Rho"), c.mass13, c.secondary.M());
        histos.fill(HIST("QAMC/hInvmassSecon_PiKa"), c.secondary.M(), c.mass23);
        histos.fill(HIST("QAMC/hInvmassSecon"), c.secondary.M());
        histos.fill(HIST("QAMC/hpT_Secondary"), c.secondary.Pt());

        // PID QA primary and secondary pion
        fillPionQA<QAFolder::MC>(pion1, true);
        fillPionQA<QAFolder::MC>(pion2, false);
        fillKaonQA<QAFolder::MC>(kaon);
      }
    };

    core.forEachCandidate<IsMC, IsMix, IsResoMicrotrack>(histos, collision, dTracks1, dTracks2, IsMC || fillQA, onPair, onCandidate);
  }

  void processResoTracks(ResoCollisions::iterator const& collision,
                         ResoTracks const& resotracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    fillHistograms<false, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoTracks, "Process ResoTracks", false);

  void processResoMicroTracks(ResoCollisions::iterator const& collision,
                              ResoMicroTracks const& resomicrotracks)
  {
    if (!core.passesEventCuts(collision)) {
      return;
    }
    fillHistograms<false, false, true>(collision, resomicrotracks, resomicrotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoMicroTracks, "Process ResoMicroTracks", true);

  void processMC(ResoMCCols::iterator const& collision,
                 ResoMCTracks const& resotracks)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    fillHistograms<true, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMC, "Process Event for MC", false);

  void processMCMicro(ResoMCCols::iterator const& collision, ResoMCMicroTracks const& tracks)
  {
    // The modular producer already selected these reconstructed collisions.
    // Apply precisely the same reconstruction loop as the data baseline.
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    histos.fill(HIST("MCReco/microTracks"), 0.5, tracks.size());
    fillHistograms<true, false, true>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCMicro, "Process reconstructed MC with micro v001 tables", false);

  void processMCTrue(ResoMCCols::iterator const& collision, ResoMCParents const& resoParents)
  {
    if (!core.passesEventCuts(collision) || !core.passesMCEventCuts(collision)) {
      return;
    }
    // Keep other/unresolved immediate decays too; never require both pairs.
    core.forEachGeneratedK1(histos, resoParents, [&](auto const& part, K1TruthChannel channel) {
      const int charge = part.pdgCode() > 0 ? 1 : -1;
      histos.fill(HIST("MCGen/chargeChannel"), charge, static_cast<int>(channel));
      histos.fill(HIST("MCGen/ptChannel"), static_cast<int>(channel), part.pt());
    });
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCTrue, "Process generated K1 in selected events with v001 parents", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(ResoCollisions const& collisions, ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<ResoCollisions, ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!core.passesEventCuts(collision1) || !core.passesEventCuts(collision2)) {
        continue;
      }
      fillHistograms<false, true, false>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processME, "Process EventMixing light without partition", false);

  // Processing Event Mixing -- Micro
  void processMEMicro(ResoCollisions const& collisions, ResoMicroTracks const& resomicrotracks)
  {
    auto tracksTuple = std::make_tuple(resomicrotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<ResoCollisions, ResoMicroTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!core.passesEventCuts(collision1) || !core.passesEventCuts(collision2)) {
        continue;
      }
      fillHistograms<false, true, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processMEMicro, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K1AnalysisMicro>(cfgc)};
}
