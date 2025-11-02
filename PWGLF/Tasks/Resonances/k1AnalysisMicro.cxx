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

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TDatabasePDG.h> // FIXME
#include <TLorentzVector.h>
#include <TPDGCode.h> // FIXME

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::constants::math;
;

struct K1AnalysisMicro {
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
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;

  //// Configurables
  Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  /// Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};

  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  /// PID Selections
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};  // Combined
  Configurable<bool> cTOFVeto{"cTOFVeto", true, "TOF Veto, if false, TOF is nessessary for PID selection"};   // TOF Veto
  Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};                // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};                // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", -999, "Combined nSigma cut for Kaon"};  // Combined
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalEvsel{"additionalEvsel", true, "Additional event selcection"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};

  // Secondary selection
  Configurable<double> cMinSecondaryPtCut{"cMinSecondaryPtCut", 0.5, "Min pT cut for secondary selection"};
  /*
  Configurable<bool> cfgModeK892orRho{"cfgModeK892orRho", false, "Secondary scenario for K892 (true) or Rho (false)"};
  Configurable<double> cSecondaryMasswindow{"cSecondaryMasswindow", 0.1, "Secondary inv mass selection window"};
  Configurable<double> cMinAnotherSecondaryMassCut{"cMinAnotherSecondaryMassCut", 0, "Min inv. mass selection of another secondary scenario"};
  Configurable<double> cMaxAnotherSecondaryMassCut{"cMaxAnotherSecondaryMassCut", 999, "MAx inv. mass selection of another secondary scenario"};
  Configurable<double> cMinPiKaMassCut{"cMinPiKaMassCut", 0, "bPion-Kaon pair inv mass selection minimum"};
  Configurable<double> cMaxPiKaMassCut{"cMaxPiKaMassCut", 999, "bPion-Kaon pair inv mass selection maximum"};
  Configurable<double> cMinAngle{"cMinAngle", 0, "Minimum angle between K(892)0 and bachelor pion"};
  Configurable<double> cMaxAngle{"cMaxAngle", 4, "Maximum angle between K(892)0 and bachelor pion"};
  Configurable<double> cMinPairAsym{"cMinPairAsym", -1, "Minimum pair asymmetry"};
  Configurable<double> cMaxPairAsym{"cMaxPairAsym", 1, "Maximum pair asymmetry"};
*/

  // K1 selection
  Configurable<double> cK1MaxRap{"cK1MaxRap", 0.5, "K1 maximum rapidity"};
  Configurable<double> cK1MinRap{"cK1MinRap", -0.5, "K1 minimum rapidity"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 200.};
    AxisSpec centAxis = {centBinning, "T0M (%)"};
    AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaxyAxis = {300, 0, 3, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {500, 0, 5, "DCA_{#it{z}} (cm)"};
    AxisSpec invMassAxisK892 = {1400 / cNbinsDiv, 0.6, 2.0, "Invariant Mass (GeV/#it{c}^2)"};   // K(892)0
    AxisSpec invMassAxisRho = {2000 / cNbinsDiv, 0.0, 2.0, "Invariant Mass (GeV/#it{c}^2)"};    // rho
    AxisSpec invMassAxisReso = {1600 / cNbinsDiv, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"}; // K1
    AxisSpec invMassAxisScan = {250, 0, 2.5, "Invariant Mass (GeV/#it{c}^2)"};                  // For selection
    AxisSpec pidQAAxis = {130, -6.5, 6.5};
    AxisSpec dataTypeAxis = {9, 0, 9, "Histogram types"};
    AxisSpec mcTypeAxis = {4, 0, 4, "Histogram types"};

    // THnSparse
    AxisSpec axisAnti = {BinAnti::kNAEnd, 0, BinAnti::kNAEnd, "Type of bin: Normal or Anti"};
    AxisSpec axisType = {BinType::kTYEnd, 0, BinType::kTYEnd, "Type of bin with charge and mix"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

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
    if (doprocessMC) {
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
    } // doprocessMC
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in K1 Analysis Task";
    histos.print();
  } // init

  // PDG code
  int kPDGRho770 = 113;
  int kK1Plus = 10323;

  template <bool IsResoMicrotrack, typename TrackType>
  bool trackCut(const TrackType track)
  {
    if constexpr (!IsResoMicrotrack) {
      // basic track cuts
      if (std::abs(track.pt()) < cMinPtcut)
        return false;
      if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
        return false;
      if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
        return false;
      if (track.tpcNClsFound() < cfgTPCcluster)
        return false;
      if (cfgHasTOF && !track.hasTOF())
        return false;
      if (cfgUseITSRefit && !track.passedITSRefit())
        return false;
      if (cfgUseTPCRefit && !track.passedTPCRefit())
        return false;
      if (cfgPVContributor && !track.isPVContributor())
        return false;
      if (cfgPrimaryTrack && !track.isPrimaryTrack())
        return false;
      if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
        return false;
      if (cfgGlobalTrack && !track.isGlobalTrack())
        return false;
    } else {
      if (std::abs(track.pt()) < cMinPtcut)
        return false;
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()) > cMaxDCArToPVcut - Epsilon)
        return false;
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags()) > cMaxDCAzToPVcut - Epsilon)
        return false;
      if (cfgPrimaryTrack && !track.isPrimaryTrack())
        return false;
      if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
        return false;
      if (cfgPVContributor && !track.isPVContributor())
        return false;
    }
    return true;
  }

  // Pion PID selection tools
  template <bool IsResoMicrotrack, typename T>
  bool selectionPIDpion(const T& candidate)
  {
    if constexpr (!IsResoMicrotrack) {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
        tpcPIDPassed = true;
      } else {
        return false;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedPion > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
          tofPIDPassed = true;
        }
      } else {
        if (!cTOFVeto) {
          return false;
        }
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      tpcPIDPassed = std::abs(o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaPiFlag())) < cMaxTPCnSigmaPion + Epsilon;
      tofPIDPassed = candidate.hasTOF() ? std::abs(o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaPiFlag())) < cMaxTOFnSigmaPion + Epsilon : true;
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return false;
  }

  // Kaon PID selection tools
  template <bool IsResoMicrotrack, typename T>
  bool selectionPIDkaon(const T& candidate)
  {
    if constexpr (!IsResoMicrotrack) {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) {
        tpcPIDPassed = true;
      } else {
        return false;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedKaon > 0) && (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa() < nsigmaCutCombinedKaon * nsigmaCutCombinedKaon)) {
          tofPIDPassed = true;
        }
      } else {
        if (!cTOFVeto) {
          return false;
        }
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      tpcPIDPassed = std::abs(o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaKaFlag())) < cMaxTPCnSigmaKaon + Epsilon;
      tofPIDPassed = candidate.hasTOF() ? std::abs(o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaKaFlag())) < cMaxTOFnSigmaKaon + Epsilon : true;
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return false;
  }

  template <typename T, typename T2>
  bool isTrueK1(const T& trk1, const T& trk2, const T2& bTrack)
  {
    if (std::abs(trk1.pdgCode()) != kPiPlus || std::abs(trk2.pdgCode()) != kPiPlus)
      return false;
    if (std::abs(bTrack.pdgCode()) != kKPlus)
      return false;
    auto mother1 = trk1.motherId();
    auto mother2 = trk2.motherId();
    if (mother1 != mother2)
      return false;
    if (((std::abs(trk1.motherPDG()) && std::abs(trk2.motherPDG()) != kPDGRho770) && (std::abs(bTrack.motherPDG()) != kK1Plus)) || (std::abs(trk1.motherPDG()) && std::abs(bTrack.motherPDG()) != kK0Star892 && (std::abs(trk2.motherPDG()) != kK1Plus)) || (std::abs(trk2.motherPDG()) && std::abs(bTrack.motherPDG()) != kK0Star892 && (std::abs(trk1.motherPDG()) != kK1Plus)))
      return false;
    auto siblings = bTrack.siblingIds();
    if (siblings[0] != mother1 && siblings[1] != mother2)
      return false;
    return true;
  } // isTrueK1

  template <typename T>
  bool isTrueK892(const T& trk1, const T& trk2)
  {
    if (std::abs(trk1.pdgCode()) != kPiPlus || std::abs(trk2.pdgCode()) != kKPlus)
      return false;
    auto mother1 = trk1.motherId();
    auto mother2 = trk2.motherId();
    if (mother1 != mother2)
      return false;
    if (std::abs(trk1.motherPDG()) != kK0Star892)
      return false;
    return true;
  }

  template <typename T>
  bool isTrueRho(const T& trk1, const T& trk2)
  {
    if (std::abs(trk1.pdgCode()) != kPiPlus || std::abs(trk2.pdgCode()) != kPiPlus)
      return false;
    auto mother1 = trk1.motherId();
    auto mother2 = trk2.motherId();
    if (mother1 != mother2)
      return false;
    if (std::abs(trk1.motherPDG()) != kPDGRho770)
      return false;
    return true;
  }
  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.cent();
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonanceSecondary, lDecayDaughter_bach, lResonanceK1;
    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks2, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs too. But the same id pairs are not needed.
      // trk1: pion, trk2: pion, bTrack: kaon
      if (!trackCut<IsResoMicrotrack>(trk1) || !trackCut<IsResoMicrotrack>(trk2))
        continue;

      auto trk1pt = trk1.pt();
      auto trk2pt = trk2.pt();
      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();

      if constexpr (!IsResoMicrotrack) {
        auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
        auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;
        auto trk2NSigmaPiTPC = trk2.tpcNSigmaPi();
        auto trk2NSigmaPiTOF = (isTrk2hasTOF) ? trk2.tofNSigmaPi() : -999.;

        if (cUseOnlyTOFTrackPi && !isTrk1hasTOF)
          continue;
        if (!selectionPIDpion<IsResoMicrotrack>(trk1) || !selectionPIDpion<IsResoMicrotrack>(trk2))
          continue;

        if constexpr (!IsMix) {

          histos.fill(HIST("QA/trkppionTPCPID"), trk1pt, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/trkppionTOFPID"), trk1pt, trk1NSigmaPiTOF);
            histos.fill(HIST("QA/trkppionTPCTOFPID"), trk1NSigmaPiTPC, trk1NSigmaPiTOF);
          }
          histos.fill(HIST("QA/trkppionpT"), trk1pt);
          histos.fill(HIST("QA/trkppionDCAxy"), trk1.dcaXY());
          histos.fill(HIST("QA/trkppionDCAz"), trk1.dcaZ());

          histos.fill(HIST("QA/trkspionTPCPID"), trk2pt, trk2NSigmaPiTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/trkspionTOFPID"), trk2pt, trk2NSigmaPiTOF);
            histos.fill(HIST("QA/trkspionTPCTOFPID"), trk2NSigmaPiTPC, trk2NSigmaPiTOF);
          }
          histos.fill(HIST("QA/trkspionpT"), trk2pt);
          histos.fill(HIST("QA/trkspionDCAxy"), trk2.dcaXY());
          histos.fill(HIST("QA/trkspionDCAz"), trk2.dcaZ());
        }
      } else {
        histos.fill(HIST("QA/trkppionTPCPID"), trk1pt, o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk1.pidNSigmaPiFlag()));
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/trkppionTOFPID"), trk1pt, o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk1.pidNSigmaPiFlag()));
          histos.fill(HIST("QA/trkppionTPCTOFPID"), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk1.pidNSigmaPiFlag()), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk1.pidNSigmaPiFlag()));
        }
        histos.fill(HIST("QA/trkppionpT"), trk1pt);
        histos.fill(HIST("QA/trkppionDCAxy"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(trk1.trackSelectionFlags()));
        histos.fill(HIST("QA/trkppionDCAz"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(trk1.trackSelectionFlags()));

        histos.fill(HIST("QA/trkspionTPCPID"), trk2pt, o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk2.pidNSigmaPiFlag()));
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/trkspionTOFPID"), trk2pt, o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk2.pidNSigmaPiFlag()));
          histos.fill(HIST("QA/trkspionTPCTOFPID"), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk2.pidNSigmaPiFlag()), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk2.pidNSigmaPiFlag()));
        }
        histos.fill(HIST("QA/trkspionpT"), trk2pt);
        histos.fill(HIST("QA/trkspionDCAxy"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(trk2.trackSelectionFlags()));
        histos.fill(HIST("QA/trkspionDCAz"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(trk2.trackSelectionFlags()));
      }

      // Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), MassPionCharged);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), MassPionCharged);
      lResonanceSecondary = lDecayDaughter1 + lDecayDaughter2;

      if (lResonanceSecondary.Pt() < cMinSecondaryPtCut)
        continue;

      if constexpr (!IsMix) {
        histos.fill(HIST("QA/hInvmassSecon"), lResonanceSecondary.M());
      }
      if constexpr (IsMC) {
        /*
           if (isTrueK892(trk1, trk2))
           histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());
           } else {
           if (isTrueRho(trk1, trk2))
           histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());
           }
         */
        histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());
      }
      // Mass Window cut is removed

      for (const auto& bTrack : dTracks1) {
        if (bTrack.index() == trk1.index() || bTrack.index() == trk2.index())
          continue;
        if (!trackCut<IsResoMicrotrack>(bTrack))
          continue;
        if (!selectionPIDkaon<IsResoMicrotrack>(bTrack))
          continue;

        // K1 reconstruction
        lDecayDaughter_bach.SetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), MassKaonCharged);
        lResonanceK1 = lResonanceSecondary + lDecayDaughter_bach;

        // Cuts
        if (lResonanceK1.Rapidity() > cK1MaxRap || lResonanceK1.Rapidity() < cK1MinRap)
          continue;

        auto lK1Angle = lResonanceSecondary.Angle(lDecayDaughter_bach.Vect());
        auto lPairAsym = (lResonanceSecondary.E() - lDecayDaughter_bach.E()) / (lResonanceSecondary.E() + lDecayDaughter_bach.E());

        TLorentzVector temp13 = lDecayDaughter1 + lDecayDaughter_bach;
        TLorentzVector temp23 = lDecayDaughter2 + lDecayDaughter_bach;

        // QA histogram
        if constexpr (!IsMix) {
          histos.fill(HIST("QA/K1OA"), lK1Angle);
          histos.fill(HIST("QA/K1PairAsym"), lPairAsym);
          histos.fill(HIST("QA/hInvmassK892_Rho"), temp13.M(), lResonanceSecondary.M());
          histos.fill(HIST("QA/hInvmassSecon_PiKa"), lResonanceSecondary.M(), temp23.M());
          histos.fill(HIST("QA/hpT_Secondary"), lResonanceSecondary.Pt());
        }
        // Selection cuts are removed
        // QA histograms after the cuts are removed as no cuts are applied

        if constexpr (!IsMix) {
          unsigned int typeK1 = bTrack.sign() > 0 ? BinType::kK1P : BinType::kK1N;
          unsigned int typeNormal = BinAnti::kNormal;
          if (trk1.sign() * trk2.sign() < 0) {
            histos.fill(HIST("k1invmass"), lResonanceK1.M());
            histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          } else {
            histos.fill(HIST("k1invmass_LS"), lResonanceK1.M());
            histos.fill(HIST("hInvmass_K1_LS"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          }

          if constexpr (IsMC) {
            if (isTrueK1(trk1, trk2, bTrack)) {
              typeK1 = bTrack.sign() > 0 ? BinType::kK1P_Rec : BinType::kK1N_Rec;
              histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
              histos.fill(HIST("k1invmass_MC"), lResonanceK1.M());
              histos.fill(HIST("QAMC/K1OA"), lK1Angle);
              histos.fill(HIST("QAMC/K1PairAsym"), lPairAsym);
              histos.fill(HIST("QAMC/hInvmassK892_Rho"), temp13.M(), lResonanceSecondary.M());
              histos.fill(HIST("QAMC/hInvmassSecon_PiKa"), lResonanceSecondary.M(), temp23.M());
              histos.fill(HIST("QAMC/hInvmassSecon"), lResonanceSecondary.M());
              histos.fill(HIST("QAMC/hpT_Seocondary"), lResonanceSecondary.Pt());

              if constexpr (!IsResoMicrotrack) {

                auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
                auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;
                auto trk2NSigmaPiTPC = trk2.tpcNSigmaPi();
                auto trk2NSigmaPiTOF = (isTrk2hasTOF) ? trk2.tofNSigmaPi() : -999.;

                // PID QA primary pion
                histos.fill(HIST("QAMC/trkppionTPCPID"), trk1pt, trk1NSigmaPiTPC);
                if (isTrk1hasTOF) {
                  histos.fill(HIST("QAMC/trkppionTOFPID"), trk1pt, trk1NSigmaPiTOF);
                  histos.fill(HIST("QAMC/trkppionTPCTOFPID"), trk1NSigmaPiTPC, trk1NSigmaPiTOF);
                }
                histos.fill(HIST("QAMC/trkppionpT"), trk1pt);
                histos.fill(HIST("QAMC/trkppionDCAxy"), trk1.dcaXY());
                histos.fill(HIST("QAMC/trkppionDCAz"), trk1.dcaZ());

                // PID QA secondary pion
                histos.fill(HIST("QAMC/trkspionTPCPID"), trk2pt, trk2NSigmaPiTPC);
                if (isTrk2hasTOF) {
                  histos.fill(HIST("QAMC/trkspionTOFPID"), trk2pt, trk2NSigmaPiTOF);
                  histos.fill(HIST("QAMC/trkspionTPCTOFPID"), trk2NSigmaPiTPC, trk2NSigmaPiTOF);
                }
                histos.fill(HIST("QAMC/trkspionpT"), trk2pt);
                histos.fill(HIST("QAMC/trkspionDCAxy"), trk2.dcaXY());
                histos.fill(HIST("QAMC/trkspionDCAz"), trk2.dcaZ());

              } else {

                histos.fill(HIST("QAMC/trkppionTPCPID"), trk1pt, o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk1.pidNSigmaSelectionFlags()));
                if (isTrk1hasTOF) {
                  histos.fill(HIST("QAMC/trkppionTOFPID"), trk1pt, o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk1.pidNSigmaSelectionFlags()));
                  histos.fill(HIST("QAMC/trkppionTPCTOFPID"), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk1.pidNSigmaSelectionFlags()), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk1.pidNSigmaSelectionFlags()));
                }
                histos.fill(HIST("QAMC/trkppionpT"), trk1pt);
                histos.fill(HIST("QAMC/trkppionDCAxy"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(trk1.trackSelectionFlags()));
                histos.fill(HIST("QAMC/trkppionDCAz"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(trk1.trackSelectionFlags()));

                // PID QA secondary pion
                histos.fill(HIST("QAMC/trkspionTPCPID"), trk2pt, o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk2.pidNSigmaSelectionFlags()));
                if (isTrk2hasTOF) {
                  histos.fill(HIST("QAMC/trkspionTOFPID"), trk2pt, o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk2.pidNSigmaSelectionFlags()));
                  histos.fill(HIST("QAMC/trkspionTPCTOFPID"), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk2.pidNSigmaSelectionFlags()), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk2.pidNSigmaSelectionFlags()));
                }
                histos.fill(HIST("QAMC/trkspionpT"), trk2pt);
                histos.fill(HIST("QAMC/trkspionDCAxy"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(trk2.trackSelectionFlags()));
                histos.fill(HIST("QAMC/trkspionDCAz"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(trk2.trackSelectionFlags()));
              }
            } else {
              histos.fill(HIST("k1invmass_MC_noK1"), lResonanceK1.M());
            }
          } // IsMC
        } else {
          unsigned int typeK1 = bTrack.sign() > 0 ? BinType::kK1P_Mix : BinType::kK1N_Mix;
          unsigned int typeNormal = BinAnti::kNormal;
          histos.fill(HIST("hInvmass_K1_Mix"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          histos.fill(HIST("k1invmass_Mix"), lResonanceK1.M());
        }
      } // bTrack
    }
  } // fillHistograms

  void processResoTracks(aod::ResoCollision const& collision,
                         aod::ResoTracks const& resotracks)
  {
    fillHistograms<false, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoTracks, "Process ResoTracks", false);

  void processResoMicroTracks(aod::ResoCollision const& collision,
                              aod::ResoMicroTracks const& resomicrotracks)
  {
    fillHistograms<false, false, true>(collision, resomicrotracks, resomicrotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoMicroTracks, "Process ResoMicroTracks", true);

  void processMC(aod::ResoCollision const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    fillHistograms<true, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMC, "Process Event for MC", false);

  void processMCTrue(ResoMCCols::iterator const& collision, aod::ResoMCParents const& resoParents)
  {
    auto multiplicity = collision.cent();
    for (const auto& part : resoParents) {
      if (std::abs(part.pdgCode()) != kK1Plus)
        continue;
      if (std::abs(part.y()) > 0.5) {
        continue;
      }
      bool pass1 = false;
      bool pass2 = false;
      bool pass3 = false;
      bool pass4 = false;
      if (std::abs(part.daughterPDG1()) == 313 || std::abs(part.daughterPDG2()) == 313) { // At least one decay into K892
        pass2 = true;
      }
      if (std::abs(part.daughterPDG1()) == kPiPlus || std::abs(part.daughterPDG2()) == kPiPlus) { // At lest one decay into pion
        pass1 = true;
      }
      if (std::abs(part.daughterPDG1()) == kPDGRho770 || std::abs(part.daughterPDG2()) == kPDGRho770) {
        pass4 = true;
      }
      if (std::abs(part.daughterPDG1()) == kKPlus || std::abs(part.daughterPDG2()) == kKPlus) {
        pass3 = true;
      }
      if (!pass1 || !pass2 || !pass3 || !pass4) // If we have both decay products
        continue;
      auto typeNormal = part.pdgCode() > 0 ? BinAnti::kNormal : BinAnti::kAnti;
      if (collision.isVtxIn10()) // INEL>10
      {
        auto typeK1 = part.pdgCode() > 0 ? BinType::kK1P_GenINEL10 : BinType::kK1N_GenINEL10;
        histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, part.pt(), 1);
      }
      if (collision.isVtxIn10() && collision.isInSel8()) // INEL>10, vtx10
      {
        auto typeK1 = part.pdgCode() > 0 ? BinType::kK1P_GenINELgt10 : BinType::kK1N_GenINELgt10;
        histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, part.pt(), 1);
      }
      if (collision.isVtxIn10() && collision.isTriggerTVX()) // vtx10, TriggerTVX
      {
        auto typeK1 = part.pdgCode() > 0 ? BinType::kK1P_GenTrig10 : BinType::kK1N_GenTrig10;
        histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, part.pt(), 1);
      }
      if (collision.isInAfterAllCuts()) // after all event selection
      {
        auto typeK1 = part.pdgCode() > 0 ? BinType::kK1P_GenEvtSel : BinType::kK1N_GenEvtSel;
        histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, part.pt(), 1);
      }
    }
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCTrue, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(o2::aod::ResoCollisions const& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      fillHistograms<false, true, false>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processME, "Process EventMixing light without partition", false);

  // Processing Event Mixing -- Micro
  // using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMEMicro(o2::aod::ResoCollisions const& collisions, aod::ResoMicroTracks const& resomicrotracks)
  {
    auto tracksTuple = std::make_tuple(resomicrotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoMicroTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      fillHistograms<false, true, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processMEMicro, "Process EventMixing light without partition", true);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K1AnalysisMicro>(cfgc)};
}
