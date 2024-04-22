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

/// \file phiInJets.cxx
/// \brief Reconstruction of Phi yield through track-track Minv correlations for resonance hadrochemistry analysis.
///
///
/// \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch>

#include <TLorentzVector.h>
#include <TVector2.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phiInJets {
  SliceCache cache;
  HistogramRegistry JEhistos{"JEhistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry registry{"registry",
                             {{"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_matched_REC_jet_pt", "matched_REC level jet pT;#it{p}_{T,jet part} (GeV/#it{c});Delta", {HistType::kTH2F, {{200, 0., 200.}, {400, -20., 20.}}}},
                              {"h_matched_REC_jet_eta", "matched_REC level jet #eta;#eta_{jet part};Delta", {HistType::kTH2F, {{100, -1.0, 1.0}, {400, -20., 20.}}}},
                              {"h_matched_REC_jet_phi", "matched_REC level jet #phi;#phi_{jet part};Delta", {HistType::kTH2F, {{80, -1.0, 7.}, {400, -20., 20.}}}},
                              {"h_matched_GEN_jet_pt", "matched_GEN level jet pT;#it{p}_{T,jet part} (GeV/#it{c});Delta", {HistType::kTH2F, {{200, 0., 200.}, {400, -20., 20.}}}},
                              {"h_matched_GEN_jet_eta", "matched_GEN level jet #eta;#eta_{jet part};Delta", {HistType::kTH2F, {{100, -1.0, 1.0}, {400, -20., 20.}}}},
                              {"h_matched_GEN_jet_phi", "matched_GEN level jet #phi;#phi_{jet part};Delta", {HistType::kTH2F, {{80, -1.0, 7.}, {400, -20., 20.}}}},
                              {"h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}}}};

  Configurable<std::string> cfgeventSelections{"cfgeventSelections", "sel8", "choose event selection"};
  Configurable<std::string> cfgtrackSelections{"cfgtrackSelections", "globalTracks", "set track selections"};

  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<int> cfgnTPCPID{"cfgnTPCPID", 4, "nTPC PID"};
  Configurable<int> cfgnTOFPID{"cfgnTOFPID", 4, "nTOF PID"};
  Configurable<float> cfgjetPtMin{"cfgjetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> cfgjetR{"cfgjetR", 0.4, "jet resolution parameter"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};
  // CONFIG DONE
  /////////////////////////////////////////  //INIT

  int eventSelection = -1;

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec axisPt{200, 0, +200, "#pt"};
    const AxisSpec MinvAxis = {400, 0.95, 1.35};
    const AxisSpec PtAxis = {200, 0, 20.0};
    const AxisSpec MultAxis = {100, 0, 100};
    const AxisSpec dRAxis = {100, 0, 3.6};

    JEhistos.add("ptJEHistogramPion", "ptJEHistogramPion", kTH1F, {PtAxis});
    JEhistos.add("ptJEHistogramKaon", "ptJEHistogramKaon", kTH1F, {PtAxis});
    JEhistos.add("ptJEHistogramProton", "ptJEHistogramProton", kTH1F, {PtAxis});
    JEhistos.add("ptJEHistogramPhi", "ptJEHistogramPhi", kTH1F, {PtAxis});
    JEhistos.add("minvJEHistogramPhi", "minvJEHistogramPhi", kTH1F, {MinvAxis});

    JEhistos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {PtAxis});
    JEhistos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {PtAxis});
    JEhistos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {PtAxis});
    JEhistos.add("ptGeneratedPhi", "ptGeneratedPhi", kTH1F, {PtAxis});
    JEhistos.add("mGeneratedPhi", "mGeneratedPhi", kTH1F, {MinvAxis});

    JEhistos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    JEhistos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
    JEhistos.add("FJetaHistogram", "FJetaHistogram", kTH1F, {axisEta});
    JEhistos.add("FJphiHistogram", "FJphiHistogram", kTH1F, {axisPhi});
    JEhistos.add("FJptHistogram", "FJptHistogram", kTH1F, {axisPt});
    JEhistos.add("FJetaHistogram_MCRec", "FJetaHistogram_MCRec", kTH1F, {axisEta});
    JEhistos.add("FJphiHistogram_MCRec", "FJphiHistogram_MCRec", kTH1F, {axisPhi});
    JEhistos.add("FJptHistogram_MCRec", "FJptHistogram_MCRec", kTH1F, {axisPt});
    JEhistos.add("FJetaHistogram_MCTrue", "FJetaHistogram_MCTrue", kTH1F, {axisEta});
    JEhistos.add("FJphiHistogram_MCTrue", "FJphiHistogram_MCTrue", kTH1F, {axisPhi});
    JEhistos.add("FJptHistogram_MCTrue", "FJptHistogram_MCTrue", kTH1F, {axisPt});

    JEhistos.add("FJnchHistogram", "FJnchHistogram", kTH1F, {MultAxis});
    JEhistos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
    JEhistos.add("nEvents_MCRec", "nEvents_MCRec", kTH1F, {{4, 0.0, 4.0}});
    JEhistos.add("nEvents_MCGen", "nEvents_MCGen", kTH1F, {{4, 0.0, 4.0}});

    JEhistos.add("nJetsPerEvent", "nJetsPerEvent", kTH1F, {{10, 0.0, 10.0}});

    JEhistos.add("hDCArToPv", "DCArToPv", kTH1F, {{300, 0.0, 3.0}});
    JEhistos.add("hDCAzToPv", "DCAzToPv", kTH1F, {{300, 0.0, 3.0}});
    JEhistos.add("rawpT", "rawpT", kTH1F, {{1000, 0.0, 10.0}});
    JEhistos.add("rawDpT", "rawDpT", kTH2F, {{1000, 0.0, 10.0}, {300, -1.5, 1.5}});
    JEhistos.add("hIsPrim", "hIsPrim", kTH1F, {{2, -0.5, +1.5}});
    JEhistos.add("hIsGood", "hIsGood", kTH1F, {{2, -0.5, +1.5}});
    JEhistos.add("hIsPrimCont", "hIsPrimCont", kTH1F, {{2, -0.5, +1.5}});
    JEhistos.add("hFindableTPCClusters", "hFindableTPCClusters", kTH1F, {{200, 0, 200}});
    JEhistos.add("hFindableTPCRows", "hFindableTPCRows", kTH1F, {{200, 0, 200}});
    JEhistos.add("hClustersVsRows", "hClustersVsRows", kTH1F, {{200, 0, 2}});
    JEhistos.add("hTPCChi2", "hTPCChi2", kTH1F, {{200, 0, 100}});
    JEhistos.add("hITSChi2", "hITSChi2", kTH1F, {{200, 0, 100}});

    JEhistos.add("hUSS", "hUSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hUSS_1D", "hUSS_1D", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_1D_2_3", "hUSS_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hLSS", "hLSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hLSS_1D", "hLSS_1D", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_1D_2_3", "hLSS_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_INSIDE", "hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hUSS_INSIDE_1D", "hUSS_INSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_INSIDE_1D_2_3", "hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_INSIDE", "hLSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hLSS_INSIDE_1D", "hLSS_INSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_INSIDE_1D_2_3", "hLSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_OUTSIDE", "hUSS_OUTSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hUSS_OUTSIDE_1D", "hUSS_OUTSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hUSS_OUTSIDE_1D_2_3", "hUSS_OUTSIDE_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_OUTSIDE", "hLSS_OUTSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hLSS_OUTSIDE_1D", "hLSS_OUTSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hLSS_OUTSIDE_1D_2_3", "hLSS_OUTSIDE_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hMCTrue_hUSS_INSIDE", "hMCTrue_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hMCTrue_hUSS_INSIDE_1D", "hMCTrue_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hMCTrue_hUSS_INSIDE_1D_2_3", "hMCTrue_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hMCTrue_hUSS_OUTSIDE", "hMCTrue_hUSS_OUTSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hMCTrue_hUSS_OUTSIDE_1D", "hMCTrue_hUSS_OUTSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hMCTrue_hUSS_OUTSIDE_1D_2_3", "hMCTrue_hUSS_OUTSIDE_1D_2_3", kTH1F, {MinvAxis});

    JEhistos.add("hMCRec_hUSS", "hMCRec_hUSS", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hMCRec_hUSS_1D", "hMCRec_hUSS_1D", kTH1F, {MinvAxis});
    JEhistos.add("hMCRec_hUSS_1D_2_3", "hMCRec_hUSS_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hMCRec_hUSS_INSIDE", "hMCRec_hUSS_INSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hMCRec_hUSS_INSIDE_1D", "hMCRec_hUSS_INSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hMCRec_hUSS_INSIDE_1D_2_3", "hMCRec_hUSS_INSIDE_1D_2_3", kTH1F, {MinvAxis});
    JEhistos.add("hMCRec_hUSS_OUTSIDE", "hMCRec_hUSS_OUTSIDE", kTH3F, {dRAxis, PtAxis, MinvAxis});
    JEhistos.add("hMCRec_hUSS_OUTSIDE_1D", "hMCRec_hUSS_OUTSIDE_1D", kTH1F, {MinvAxis});
    JEhistos.add("hMCRec_hUSS_OUTSIDE_1D_2_3", "hMCRec_hUSS_OUTSIDE_1D_2_3", kTH1F, {MinvAxis});

    // EVENT SELECTION
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(cfgeventSelections));

  } // end of init

  double massKa = o2::constants::physics::MassKPlus;

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>; // , aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa>;
  Filter jetCuts = aod::jet::pt > cfgjetPtMin&& aod::jet::r == nround(cfgjetR.node() * 100.0f);

  // Function for track quality cuts
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <typename TrackType>
  bool trackSelection(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cfgtrkMinPt)
      return false;

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (track.tpcNClsFindable() < cfgnFindableTPCClusters)
      return false;

    if (std::abs(track.tpcNClsCrossedRows()) < cfgnTPCCrossedRows)
      return false;

    if (std::abs(track.tpcCrossedRowsOverFindableCls()) > cfgnRowsOverFindable)
      return false;

    if (std::abs(track.tpcChi2NCl()) > cfgnTPCChi2)
      return false;

    if (std::abs(track.itsChi2NCl()) > cfgnITSChi2)
      return false;

    if (cfgConnectedToPV && !track.isPVContributor())
      return false;
    return true;
  };
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <typename T>
  bool trackPID(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < cfgnTPCPID)
      tpcPIDPassed = true;

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cfgnTOFPID) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <bool IsMC, bool IsMix, typename TracksType, typename JetType>
  void minvReconstruction(double mult, const TracksType& trk1, const TracksType& trk2, const JetType& jets)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    if (!trackSelection(trk1) || !trackSelection(trk2))
      return;

    if (!trackPID(trk1) || !trackPID(trk2))
      return;

    if (trk1.index() >= trk2.index())
      return; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

    lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
    lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Rapidity()) > 0.5)
      return;

    /////////////////////////////////////////////////////////////////////////////
    // Fill Global Event Minv
    if (trk1.sign() * trk2.sign() < 0) {
      if (!IsMC) {
        JEhistos.fill(HIST("hUSS_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hUSS_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hUSS"), mult, lResonance.Pt(), lResonance.M());
      }

      if (IsMC) {
        JEhistos.fill(HIST("hMCRec_hUSS_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hMCRec_hUSS_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hMCRec_hUSS"), mult, lResonance.Pt(), lResonance.M());
      }

    } else if (trk1.sign() * trk2.sign() > 0) {

      JEhistos.fill(HIST("hLSS_1D"), lResonance.M());
      if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
        JEhistos.fill(HIST("hLSS_1D_2_3"), lResonance.M());
      JEhistos.fill(HIST("hLSS"), mult, lResonance.Pt(), lResonance.M());
    }
    /////////////////////////////////////////////////////////////////////////////

    bool jetFlag = false;
    for (auto const& jet : jets) {

      double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
      double etadiff = jet.eta() - lResonance.Eta();
      double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
      if (R < cfgjetR)
        jetFlag = true;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Fill inside Jet
    if (jetFlag) {
      if (trk1.sign() * trk2.sign() < 0) {
        if (!IsMC) {
          JEhistos.fill(HIST("hUSS_INSIDE_1D"), lResonance.M());
          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hUSS_INSIDE_1D_2_3"), lResonance.M());
          JEhistos.fill(HIST("hUSS_INSIDE"), mult, lResonance.Pt(), lResonance.M());
        }

        if (IsMC) {
          JEhistos.fill(HIST("hMCRec_hUSS_INSIDE_1D"), lResonance.M());
          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hMCRec_hUSS_INSIDE_1D_2_3"), lResonance.M());
          JEhistos.fill(HIST("hMCRec_hUSS_INSIDE"), mult, lResonance.Pt(), lResonance.M());
        }

      } else if (trk1.sign() * trk2.sign() > 0) {

        JEhistos.fill(HIST("hLSS_INSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hLSS_INSIDE_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hLSS_INSIDE"), mult, lResonance.Pt(), lResonance.M());
      }
    } // jetflag
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    // Fill outside Jet
    if (!jetFlag) {
      if (trk1.sign() * trk2.sign() < 0) {
        if (!IsMC) {
          JEhistos.fill(HIST("hUSS_OUTSIDE_1D"), lResonance.M());
          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hUSS_OUTSIDE_1D_2_3"), lResonance.M());
          JEhistos.fill(HIST("hUSS_OUTSIDE"), mult, lResonance.Pt(), lResonance.M());
        }

        if (IsMC) {
          JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_1D"), lResonance.M());
          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_1D_2_3"), lResonance.M());
          JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE"), mult, lResonance.Pt(), lResonance.M());
        }

      } else if (trk1.sign() * trk2.sign() > 0) {

        JEhistos.fill(HIST("hLSS_OUTSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hLSS_OUTSIDE_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hLSS_OUTSIDE"), mult, lResonance.Pt(), lResonance.M());
      }
    } //! jetflag
    /////////////////////////////////////////////////////////////////////////////
  } // MinvReconstruction

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  int nEvents = 0;
  void processJetTracks(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& chargedjets, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, TrackCandidates const&)
  {
    if (cDebugLevel > 0) {
      nEvents++;
      if ((nEvents + 1) % 10000 == 0)
        std::cout << nEvents << std::endl;
    }
    JEhistos.fill(HIST("nEvents"), 0.5);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;
    if (fabs(collision.posZ()) > 10)
      return;

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      auto trk1 = track1.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>();
      auto trk2 = track2.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>();
      minvReconstruction<false, false>(1.0, trk1, trk2, chargedjets);
    }
    int nJets = 0;
    for (auto chargedjet : chargedjets) {
      JEhistos.fill(HIST("FJetaHistogram"), chargedjet.eta());
      JEhistos.fill(HIST("FJphiHistogram"), chargedjet.phi());
      JEhistos.fill(HIST("FJptHistogram"), chargedjet.pt());
      JEhistos.fill(HIST("FJnchHistogram"), chargedjet.tracks().size());
      nJets++;
    }
    JEhistos.fill(HIST("nJetsPerEvent"), nJets);

    JEhistos.fill(HIST("nEvents"), 1.5);

    for (auto const& track : tracks) {
      auto originalTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>();
      JEhistos.fill(HIST("hDCArToPv"), originalTrack.dcaXY());
      JEhistos.fill(HIST("hDCAzToPv"), originalTrack.dcaZ());
      JEhistos.fill(HIST("rawpT"), originalTrack.pt());
      JEhistos.fill(HIST("rawDpT"), track.pt(), track.pt() - originalTrack.pt());
      JEhistos.fill(HIST("hIsPrim"), originalTrack.isPrimaryTrack());
      JEhistos.fill(HIST("hIsGood"), originalTrack.isGlobalTrackWoDCA());
      JEhistos.fill(HIST("hIsPrimCont"), originalTrack.isPVContributor());
      JEhistos.fill(HIST("hFindableTPCClusters"), originalTrack.tpcNClsFindable());
      JEhistos.fill(HIST("hFindableTPCRows"), originalTrack.tpcNClsCrossedRows());
      JEhistos.fill(HIST("hClustersVsRows"), originalTrack.tpcCrossedRowsOverFindableCls());
      JEhistos.fill(HIST("hTPCChi2"), originalTrack.tpcChi2NCl());
      JEhistos.fill(HIST("hITSChi2"), originalTrack.itsChi2NCl());

      if (!trackSelection(originalTrack))
        continue;

      JEhistos.fill(HIST("etaHistogram"), track.eta());
      JEhistos.fill(HIST("phiHistogram"), track.phi());
    } // JTrack Loop

  }; // Process Switch
  PROCESS_SWITCH(phiInJets, processJetTracks, "process JE Framework", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////MC STUFF////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using myCompleteJetTracks = soa::Join<aod::JTracks, aod::JTrackPIs, aod::McTrackLabels>;
  int nJEEvents = 0;
  void processRec(o2::aod::JCollision const& collision, myCompleteJetTracks const& tracks, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, aod::McParticles const&, myCompleteTracks const&)
  {
    if (cDebugLevel > 0) {
      nJEEvents++;
      if ((nJEEvents + 1) % 10000 == 0)
        std::cout << "JEvents: " << nJEEvents << std::endl;
    }

    JEhistos.fill(HIST("nEvents_MCRec"), 0.5);
    if (fabs(collision.posZ()) > 10)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;
    JEhistos.fill(HIST("nEvents_MCRec"), 1.5);

    // Track Eff
    for (const auto& track : tracks) {
      auto originalTrack = track.track_as<myCompleteTracks>();
      if (!trackSelection(originalTrack))
        continue;
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y()) <= 0.5) { // do this in the context of the track ! (context matters!!!)
          if (abs(mcParticle.pdgCode()) == 211)
            JEhistos.fill(HIST("ptJEHistogramPion"), mcParticle.pt());
          if (abs(mcParticle.pdgCode()) == 321)
            JEhistos.fill(HIST("ptJEHistogramKaon"), mcParticle.pt());
          if (abs(mcParticle.pdgCode()) == 2212)
            JEhistos.fill(HIST("ptJEHistogramProton"), mcParticle.pt());
        }
      }
    }

    // Jet Eff
    for (auto& mcdjet : mcdjets) {
      registry.fill(HIST("h_jet_pt"), mcdjet.pt());
      registry.fill(HIST("h_jet_eta"), mcdjet.eta());
      registry.fill(HIST("h_jet_phi"), mcdjet.phi());
    }
    // Phi Eff

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      auto trk1 = track1.track_as<myCompleteTracks>();
      auto trk2 = track2.track_as<myCompleteTracks>();
      if (trk1.index() >= trk2.index())
        continue;
      if (fabs(trk1.eta()) > 0.8 || fabs(trk2.eta()) > 0.8)
        continue;

      if (!trackSelection(trk1))
        continue;
      if (!trackSelection(trk2))
        continue;
      if ((trk1.sign() * trk2.sign()) > 0)
        continue; // Not K+K-
      if (!track1.has_mcParticle())
        continue;
      if (!track2.has_mcParticle())
        continue;

      auto part1 = track1.mcParticle();
      auto part2 = track2.mcParticle();
      if (fabs(part1.pdgCode()) != 321)
        continue; // Not Kaon
      if (fabs(part2.pdgCode()) != 321)
        continue; // Not Kaon
      if (!part1.has_mothers())
        continue; // Not decaying Kaon
      if (!part2.has_mothers())
        continue; // Not decaying Kaon

      std::vector<int> mothers1{};
      std::vector<int> mothers1PDG{};
      for (auto& part1_mom : part1.mothers_as<aod::McParticles>()) {
        mothers1.push_back(part1_mom.globalIndex());
        mothers1PDG.push_back(part1_mom.pdgCode());
      }
      if (mothers1.size() > 1)
        continue; // Strange multi-mother decay

      std::vector<int> mothers2{};
      std::vector<int> mothers2PDG{};
      for (auto& part2_mom : part2.mothers_as<aod::McParticles>()) {
        mothers2.push_back(part2_mom.globalIndex());
        mothers2PDG.push_back(part2_mom.pdgCode());
      }

      if (mothers2.size() > 1)
        continue; // Strange multi-mother decay
      if (mothers1PDG[0] != 333)
        continue; // mother not phi
      if (mothers2PDG[0] != 333)
        continue; // mother not phi
      if (mothers1[0] != mothers2[0])
        continue; // Kaons not from the same phi
      TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      if (lResonance.Rapidity() > 0.5)
        continue;
      JEhistos.fill(HIST("ptJEHistogramPhi"), lResonance.Pt());
      JEhistos.fill(HIST("minvJEHistogramPhi"), lResonance.M());

    } // end of phi check
  }
  PROCESS_SWITCH(phiInJets, processRec, "pikp detector level MC JE", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void processSim(o2::aod::JMcCollision const& collision, aod::JMcParticles const& mcParticles, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  {
    JEhistos.fill(HIST("nEvents_MCGen"), 0.5);

    if (fabs(collision.posZ()) > 10)
      return;
    JEhistos.fill(HIST("nEvents_MCGen"), 1.5);

    // Check pikp and phi
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && fabs(mcParticle.y()) <= 0.5) { // watch out for context!!!
        if (abs(mcParticle.pdgCode()) == 211)
          JEhistos.fill(HIST("ptGeneratedPion"), mcParticle.pt());
        if (abs(mcParticle.pdgCode()) == 321)
          JEhistos.fill(HIST("ptGeneratedKaon"), mcParticle.pt());
        if (abs(mcParticle.pdgCode()) == 2212)
          JEhistos.fill(HIST("ptGeneratedProton"), mcParticle.pt());
      }
      if (fabs(mcParticle.y()) <= 0.5) { // watch out for context!!!
        TLorentzVector lResonance;
        lResonance.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
        if (abs(mcParticle.pdgCode()) == 333)
          JEhistos.fill(HIST("ptGeneratedPhi"), mcParticle.pt());
        if (abs(mcParticle.pdgCode()) == 333)
          JEhistos.fill(HIST("mGeneratedPhi"), lResonance.M());
      }
    }
    // Check jets
    for (auto& mcpjet : mcpjets) {
      registry.fill(HIST("h_part_jet_pt"), mcpjet.pt());
      registry.fill(HIST("h_part_jet_eta"), mcpjet.eta());
      registry.fill(HIST("h_part_jet_phi"), mcpjet.phi());
    }
  }
  PROCESS_SWITCH(phiInJets, processSim, "pikp particle level MC", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  using JetMCDTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;

  // void processMatchedGen(o2::aod::JMcCollision const& collision, aod::JMcParticles const& mcParticles, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  void processMatchedGen(aod::JMcCollision const&,
                         JetMCDTable const&,
                         JetMCPTable const& mcpjets,
                         aod::JMcParticles const& mcParticles)

  {
    // if(fabs(collision.posZ())>10)
    //   return;

    std::vector<double> mcd_pt{};
    std::vector<double> mcd_phi{};
    std::vector<double> mcd_eta{};
    std::vector<double> mcp_pt{};
    std::vector<double> mcp_phi{};
    std::vector<double> mcp_eta{};

    for (auto& mcpjet : mcpjets) {
      if (!mcpjet.has_matchedJetGeo())
        continue;

      for (auto& mcdjet : mcpjet.template matchedJetGeo_as<JetMCDTable>()) {
        if (!mcdjet.has_matchedJetGeo())
          continue;

        registry.fill(HIST("h_matched_GEN_jet_pt"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        registry.fill(HIST("h_matched_GEN_jet_phi"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        registry.fill(HIST("h_matched_GEN_jet_eta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        mcd_pt.push_back(mcdjet.pt());
        mcd_eta.push_back(mcdjet.eta());
        mcd_phi.push_back(mcdjet.phi());
        mcp_pt.push_back(mcpjet.pt());
        mcp_eta.push_back(mcpjet.eta());
        mcp_phi.push_back(mcpjet.phi());
      } // mcpjets
    }   // mcdjets

    // First we do GEN part
    for (const auto& mcParticle : mcParticles) {
      if (fabs(mcParticle.y()) > 0.5)
        continue;
      if (fabs(mcParticle.eta() > 0.8))
        continue;
      if (abs(mcParticle.pdgCode()) == 333) {

        TLorentzVector lResonance;
        lResonance.SetPxPyPzE(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());
        bool jetFlag = false;
        for (int i = 0; i < mcp_pt.size(); i++) {
          double phidiff = TVector2::Phi_mpi_pi(mcp_pt[i] - lResonance.Phi());
          double etadiff = mcp_eta[i] - lResonance.Eta();
          double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
          if (R < cfgjetR)
            jetFlag = true;
        }

        if (jetFlag) {
          JEhistos.fill(HIST("hMCTrue_hUSS_INSIDE_1D"), lResonance.M());
          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hMCTrue_hUSS_INSIDE_1D_2_3"), lResonance.M());
          JEhistos.fill(HIST("hMCTrue_hUSS_INSIDE"), 1.0, lResonance.Pt(), lResonance.M());

        } // jetflag
        if (!jetFlag) {
          JEhistos.fill(HIST("hMCTrue_hUSS_OUTSIDE_1D"), lResonance.M());

          if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
            JEhistos.fill(HIST("hMCTrue_hUSS_OUTSIDE_1D_2_3"), lResonance.M());

          JEhistos.fill(HIST("hMCTrue_hUSS_OUTSIDE"), 1.0, lResonance.Pt(), lResonance.M());

        } //! jetflag
      }   // chech for phi
    }     // MC Particles
  }       // main fcn
  PROCESS_SWITCH(phiInJets, processMatchedGen, "phi matched level MC", true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //  void processMatchedRec(o2::aod::JCollision const& collision, myCompleteJetTracks const& tracks, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, aod::McParticles const&, myCompleteTracks const& originalTracks)
  void processMatchedRec(aod::JCollision const& collision,
                         JetMCDTable const& mcdjets,
                         JetMCPTable const&,
                         myCompleteJetTracks const& tracks,
                         myCompleteTracks const&,
                         aod::McParticles const&)
  {
    if (fabs(collision.posZ()) > 10)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    std::vector<double> mcd_pt{};
    std::vector<double> mcd_phi{};
    std::vector<double> mcd_eta{};
    std::vector<double> mcp_pt{};
    std::vector<double> mcp_phi{};
    std::vector<double> mcp_eta{};

    // std::vector<double> mcp_ID{};
    // std::vector<double> mcp_ID{};

    for (auto& mcdjet : mcdjets) {
      if (!mcdjet.has_matchedJetGeo())
        continue;
      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {
        if (!mcpjet.has_matchedJetGeo())
          continue;

        registry.fill(HIST("h_matched_REC_jet_pt"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        registry.fill(HIST("h_matched_REC_jet_phi"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        registry.fill(HIST("h_matched_REC_jet_eta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        mcd_pt.push_back(mcdjet.pt());
        mcd_eta.push_back(mcdjet.eta());
        mcd_phi.push_back(mcdjet.phi());
        mcp_pt.push_back(mcpjet.pt());
        mcp_eta.push_back(mcpjet.eta());
        mcp_phi.push_back(mcpjet.phi());
      } // mcpjets
    }   // mcdjets

    // Now we do REC part
    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      auto trk1 = track1.track_as<myCompleteTracks>();
      auto trk2 = track2.track_as<myCompleteTracks>();
      if (trk1.index() >= trk2.index())
        continue;
      if (fabs(trk1.eta()) > 0.8 || fabs(trk2.eta()) > 0.8)
        continue;
      if (!trackSelection(trk1))
        continue;
      if (!trackSelection(trk2))
        continue;
      if ((trk1.sign() * trk2.sign()) > 0)
        continue; // Not K+K-
      if (!track1.has_mcParticle())
        continue;
      if (!track2.has_mcParticle())
        continue;

      auto part1 = track1.mcParticle();
      auto part2 = track2.mcParticle();
      if (fabs(part1.pdgCode()) != 321)
        continue; // Not Kaon
      if (fabs(part2.pdgCode()) != 321)
        continue; // Not Kaon
      if (!part1.has_mothers())
        continue; // Not decaying Kaon
      if (!part2.has_mothers())
        continue; // Not decaying Kaon

      std::vector<int> mothers1{};
      std::vector<int> mothers1PDG{};
      for (auto& part1_mom : part1.mothers_as<aod::McParticles>()) {
        mothers1.push_back(part1_mom.globalIndex());
        mothers1PDG.push_back(part1_mom.pdgCode());
      }
      if (mothers1.size() > 1)
        continue; // Strange multi-mother decay

      std::vector<int> mothers2{};
      std::vector<int> mothers2PDG{};
      for (auto& part2_mom : part2.mothers_as<aod::McParticles>()) {
        mothers2.push_back(part2_mom.globalIndex());
        mothers2PDG.push_back(part2_mom.pdgCode());
      }
      if (mothers2.size() > 1)
        continue; // Strange multi-mother decay
      if (mothers1PDG[0] != 333)
        continue; // mother not phi
      if (mothers2PDG[0] != 333)
        continue; // mother not phi
      if (mothers1[0] != mothers2[0])
        continue; // Kaons not from the same phi

      TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      if (lResonance.Rapidity() > 0.5)
        continue;

      bool jetFlag = false;
      for (int i = 0; i < mcd_pt.size(); i++) {
        double phidiff = TVector2::Phi_mpi_pi(mcd_pt[i] - lResonance.Phi());
        double etadiff = mcd_eta[i] - lResonance.Eta();
        double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
        if (R < cfgjetR)
          jetFlag = true;
      }

      if (jetFlag) {
        JEhistos.fill(HIST("hMCRec_hUSS_INSIDE_1D"), lResonance.M());
        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hMCRec_hUSS_INSIDE_1D_2_3"), lResonance.M());
        JEhistos.fill(HIST("hMCRec_hUSS_INSIDE"), 1.0, lResonance.Pt(), lResonance.M());

      } // jetflag
      if (!jetFlag) {
        JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_1D"), lResonance.M());

        if (lResonance.Pt() > 2.0 && lResonance.Pt() < 3)
          JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE_1D_2_3"), lResonance.M());

        JEhistos.fill(HIST("hMCRec_hUSS_OUTSIDE"), 1.0, lResonance.Pt(), lResonance.M());

      } //! jetflag
    }   // tracks
  }     // main fcn
  PROCESS_SWITCH(phiInJets, processMatchedRec, "phi matched Rec level MC", true);

}; // end of main struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiInJets>(cfgc)};
};
