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

/// \file kstarqa.cxx
/// \brief Code for Kstar resonance without resonance initializer
/// \author prottay das, sawan
/// \since 13/03/2024

// #include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <vector>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct Kstarqa {

  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hInvMass{"hInvMass", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hPID{"hPID", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hOthers{"hOthers", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Confugrable for QA histograms
  Configurable<bool> calcLikeSign{"calcLikeSign", true, "Calculate Like Sign"};
  Configurable<bool> calcRotational{"calcRotational", true, "Calculate Rotational"};
  Configurable<bool> cQAplots{"cQAplots", true, "cQAplots"};
  Configurable<bool> cQAevents{"cQAevents", true, "Multiplicity dist, DCAxy, DCAz"};
  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
  Configurable<bool> onlyTPC{"onlyTPC", true, "only TPC tracks"};

  // Configurables for track selections
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgRCRFC{"cfgRCRFC", 0.8f, "Crossed Rows to Findable Clusters"};
  Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 36.0, "ITS Chi2/NCl"};
  Configurable<float> cfgTPCChi2NCl{"cfgTPCChi2NCl", 4.0, "TPC Chi2/NCl"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                       // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<float> cBetaCutTOF{"cBetaCutTOF", 0.0, "cut TOF beta"};
  Configurable<float> cFakeTrackCutKa{"cFakeTrackCutKa", 0.5, "Cut based on momentum difference in global and TPC tracks for Kaons"};
  Configurable<float> cFakeTrackCutPi{"cFakeTrackCutPi", 0.5, "Cut based on momentum difference in global and TPC tracks for Pions"};
  Configurable<bool> cFakeTrack{"cFakeTrack", true, "Fake track selection"};

  // PID selections
  Configurable<float> nsigmaCutTPCPi{"nsigmaCutTPCPi", 3.0, "TPC Nsigma cut for pions"};
  Configurable<float> nsigmaCutTPCKa{"nsigmaCutTPCKa", 3.0, "TPC Nsigma cut for kaons"};
  Configurable<float> nsigmaCutTOFPi{"nsigmaCutTOFPi", 3.0, "TOF Nsigma cut for pions"};
  Configurable<float> nsigmaCutTOFKa{"nsigmaCutTOFKa", 3.0, "TOF Nsigma cut for kaons"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Combined Nsigma cut"};

  // Event selection configurables
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> cTVXEvsel{"cTVXEvsel", false, "Triggger selection"};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  // Configurable<bool> cMID{"cMID", false, "Misidentification of tracks"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "No. of bins in Vz distribution"};
  Configurable<int> nBinsinvMass{"nBinsinvMass", 180, "N bins in invMass hInvMass"};
  Configurable<float> invMassbinlow{"invMassbinlow", 0.6, "invMass bin low"};
  Configurable<float> invMassbinhigh{"invMassbinhigh", 1.5, "invMass bin high"};
  Configurable<int> nBinspT{"nBinspT", 200, "N bins in pT hInvMass"};
  Configurable<float> pTbinlow{"pTbinlow", 0.0, "pT bin low"};
  Configurable<float> pTbinhigh{"pTbinhigh", 20.0, "pT bin high"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", true, "avoid split track in MC"};
  Configurable<bool> cAllGenCollisions{"cAllGenCollisions", false, "To fill all generated collisions for the signal loss calculations"};
  ConfigurableAxis binsMultPlot{"binsMultPlot", {201, -0.5f, 200.5f}, "centrality axis bins"};
  ConfigurableAxis axisdEdx{"axisdEdx", {1, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {1, 0, 20}, "pT (GeV/c)"};
  ConfigurableAxis axisMultdist{"axisMultdist", {1, 0, 70000}, "Multiplicity distribution"};

  // Event plane configurables
  Configurable<bool> boostDaugter1{"boostDaugter1", false, "Boost daughter Kaon in the COM frame"};
  Configurable<bool> boostDaugter2{"boostDaugter2", true, "Boost daughter Pion in the COM frame"};
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
  TRandom* rn = new TRandom();

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm] for plots"};
    AxisSpec ptAxis = {nBinspT, pTbinlow, pTbinhigh, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invmassAxis = {nBinsinvMass, invMassbinlow, invMassbinhigh, "Invariant mass (GeV/#it{c}^{2})"};
    AxisSpec thnAxisPOL{configThnAxisPOL, "cos(#theta)"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Multiplicity percentile", kTH1F, {{binsMultPlot}});

    // for primary tracks
    if (cQAplots) {
      hOthers.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
      hOthers.add("hphi", "Phi distribution", kTH1F, {{65, 0, 6.5}});
      hOthers.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hCRFC_after", "CRFC after distribution", kTH1F, {{100, 0.0f, 10.0f}});
      hOthers.add("hCRFC_before", "CRFC before distribution", kTH1F, {{100, 0.0f, 10.0f}});

      hPID.add("Before/hNsigmaTPC_Ka_before", "N #sigma Kaon TPC before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("Before/hNsigmaTOF_Ka_before", "N #sigma Kaon TOF before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("Before/hNsigmaTPC_Pi_before", "N #sigma Pion TPC before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("Before/hNsigmaTOF_Pi_before", "N #sigma Pion TOF before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("Before/hNsigma_TPC_TOF_Ka_before", "N #sigma Kaon TOF before", kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      hPID.add("Before/hNsigma_TPC_TOF_Pi_before", "N #sigma Pion TOF before", kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      hPID.add("h1PID_TPC_kaon_data", "Kaon PID distribution in data", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TPC_pion_data", "Pion PID distribution in data", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TPC_kaon_MC", "Kaon PID distribution in MC", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TPC_pion_MC", "Pion PID distribution in MC", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TOF_kaon_data", "Kaon PID distribution in data", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TOF_pion_data", "Pion PID distribution in data", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TOF_kaon_MC", "Kaon PID distribution in MC", kTH1F, {{100, -10.0f, 10.0f}});
      hPID.add("h1PID_TOF_pion_MC", "Pion PID distribution in MC", kTH1F, {{100, -10.0f, 10.0f}});

      hPID.add("After/hNsigmaPionTPC_after", "N #Pi TPC after", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("After/hNsigmaPionTOF_after", "N #Pi TOF after", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("After/hNsigmaKaonTPC_after", "N #sigma Kaon TPC after", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("After/hNsigmaKaonTOF_after", "N #sigma Kaon TOF after", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      hPID.add("After/hNsigma_TPC_TOF_Ka_after", "N #sigma Kaon TOF after", kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      hPID.add("After/hNsigma_TPC_TOF_Pi_after", "N #sigma Pion TOF after", kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
    }

    // KStar histograms
    hInvMass.add("h3KstarInvMassUnlikeSign", "kstar Unlike Sign", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL});
    hInvMass.add("h3KstarInvMassMixed", "kstar Mixed", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL});
    if (calcLikeSign)
      hInvMass.add("h3KstarInvMasslikeSign", "kstar like Sign", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL});
    if (calcRotational)
      hInvMass.add("h3KstarInvMassRotated", "kstar rotated", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL});

    // MC generated histograms
    hInvMass.add("hk892GenpT", "pT distribution of True MC K(892)0", kTH2F, {ptAxis, binsMultPlot});
    // hInvMass.add("hk892GenpTAnti", "pT distribution of True MC Anti-K(892)0", kTH2F, {ptAxis}, {binsMultPlot});
    // Reconstructed MC histogram
    hInvMass.add("h1KstarRecMass", "Invariant mass of kstar meson", kTH1F, {invmassAxis});
    hInvMass.add("h2KstarRecpt1", "pT of kstar meson", kTH2F, {ptAxis, binsMultPlot});
    hInvMass.add("h2KstarRecpt2", "pT of generated kstar meson", kTH2F, {ptAxis, binsMultPlot});
    hInvMass.add("h1genmass", "Invariant mass of generated kstar meson", kTH1F, {invmassAxis});
    rEventSelection.add("events_check_data", "No. of events in the data", kTH1I, {{20, 0, 20}});
    rEventSelection.add("events_check", "No. of events in the generated MC", kTH1I, {{20, 0, 20}});
    rEventSelection.add("events_checkrec", "No. of events in the reconstructed MC", kTH1I, {{20, 0, 20}});
    hInvMass.add("h1KSRecsplit", "KS meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});

    // Multplicity distribution
    if (cQAevents) {
      rEventSelection.add("multdist_FT0M", "FT0M Multiplicity distribution", kTH1F, {axisMultdist});
      // hInvMass.add("multdist_FT0A", "FT0A Multiplicity distribution", kTH1F, {axisMultdist});
      // hInvMass.add("multdist_FT0C", "FT0C Multiplicity distribution", kTH1F, {axisMultdist});
      // hInvMass.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
      rEventSelection.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
      rEventSelection.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    }
  }

  // double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass(); // FIXME: Get from the common header
  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster)) {
      return false;
    } else if (!ismanualDCAcut) {
      if (std::abs(candidate.pt()) < cfgCutPT)
        return false;
      if (std::abs(candidate.dcaXY()) > cfgCutDCAxy)
        return false;
      if (std::abs(candidate.dcaZ()) > cfgCutDCAz)
        return false;
      if (candidate.tpcCrossedRowsOverFindableCls() < cfgRCRFC)
        return false;
      if (candidate.itsNCls() < cfgITScluster)
        return false;
      if (candidate.tpcNClsFound() < cfgTPCcluster)
        return false;
      if (candidate.itsChi2NCl() >= cfgITSChi2NCl)
        return false;
      if (candidate.tpcChi2NCl() >= cfgTPCChi2NCl)
        return false;
      if (cfgPVContributor && !candidate.isPVContributor())
        return false;
      if (cfgPrimaryTrack && !candidate.isPrimaryTrack())
        return false;
      if (cfgGlobalWoDCATrack && !candidate.isGlobalTrackWoDCA())
        return false;
      if (cfgGlobalTrack && !candidate.isGlobalTrack())
        return false;
    }

    return true;
  }

  template <typename T>
  bool isFakeTrack(const T& track, int PID)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (PID == 0 && std::abs(pglobal - ptpc) > cFakeTrackCutPi) {
      return true;
    }
    if (PID == 1 && std::abs(pglobal - ptpc) > cFakeTrackCutKa) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi && candidate.beta() > cBetaCutTOF) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (nsigmaCutCombined * nsigmaCutCombined) && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa && candidate.beta() > cBetaCutTOF) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined) && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      }
    }
    return false;
  }

  // template <typename T>
  // bool cMIDselectionPID(const T& candidate, int PID)
  // {
  //   if (PID == 0) {
  //     if (onlyTOF) {
  //       if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < 3.0) {
  //         return true;
  //       }
  //     } else if (onlyTOFHIT) {
  //       if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < 3.0) {
  //         return true;
  //       }
  //       if (!candidate.hasTOF() &&
  //           std::abs(candidate.tpcNSigmaPi()) < 3.0) {
  //         return true;
  //       }
  //     } else if (onlyTPC) {
  //       if (std::abs(candidate.tpcNSigmaPi()) < 3.0) {
  //         return true;
  //       }
  //     } else {
  //       if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (3.0 * 3.0)) {
  //         return true;
  //       }
  //       if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < 3.0) {
  //         return true;
  //       }
  //     }
  //   } else if (PID == 1) {
  //     if (onlyTOF) {
  //       if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < 3.0) {
  //         return true;
  //       }
  //     } else if (onlyTOFHIT) {
  //       if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < 3.0) {
  //         return true;
  //       }
  //       if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
  //         return true;
  //       }
  //     } else if (onlyTPC) {
  //       if (std::abs(candidate.tpcNSigmaKa()) < 3.0) {
  //         return true;
  //       }
  //     } else {
  //       if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (3.0 * 3.0)) {
  //         return true;
  //       }
  //       if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
  //         return true;
  //       }
  //     }
  //   } else if (PID == 2) {
  //     if (onlyTOF) {
  //       if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < 3.0) {
  //         return true;
  //       }
  //     } else if (onlyTOFHIT) {
  //       if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < 3.0) {
  //         return true;
  //       }
  //       if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
  //         return true;
  //       }
  //     } else if (onlyTPC) {
  //       if (std::abs(candidate.tpcNSigmaPr()) < 3.0) {
  //         return true;
  //       }
  //     } else {
  //       if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (3.0 * 3.0)) {
  //         return true;
  //       }
  //       if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
  //         return true;
  //       }
  //     }
  //   }
  //   return false;
  // }

  std::array<float, 3> pvec0;
  std::array<float, 3> pvec1;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter fDCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTOFbeta>>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;

  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::McTrackLabels, aod::pidTOFbeta>>;

  //*********Varibles declaration***************
  TLorentzVector lv1, lv2, lv3, lv4, lv5;
  float multiplicity = 0.0f;
  float theta2;
  ROOT::Math::PxPyPzMVector daughter1, daughter2, daughterSelected, fourVecDau1, fourVecMother, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, helicityVec, randomVec, beamVec, normalVec;
  bool isMix = false;

  template <typename T1, typename T2, typename T3, typename T4>
  void fillInvMass(const T1& track1, const T2& track2, const T3& lv2, const T4& lv3, float multiplicity, bool isMix)
  {
    daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa); // Kaon
    daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi); // Pion
    daughterSelected = (boostDaugter1) ? daughter1 : daughter2;
    auto selectedDauMass = (boostDaugter1) ? massKa : massPi;

    // polarization calculations

    fourVecDau1 = ROOT::Math::PxPyPzMVector(daughterSelected.Px(), daughterSelected.Py(), daughterSelected.Pz(), selectedDauMass); // Kaon or Pion

    fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M()); // mass of KshortKshort pair
    ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                               // boost mother to center of mass frame
    fourVecDauCM = boost(fourVecDau1);                                                // boost the frame of daughter same as mother
    threeVecDauCM = fourVecDauCM.Vect();                                              // get the 3 vector of daughter in the frame of mother

    if (std::abs(lv3.Rapidity()) < 0.5) {
      if (activateTHnSparseCosThStarHelicity) {
        helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
        auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);

            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (calcRotational)
                hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarHelicity);
            }
          } else {
            hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);
          }
        } else {
          if (!isMix) {
            if (calcLikeSign)
              hInvMass.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);
          }
        }

      } else if (activateTHnSparseCosThStarProduction) {
        normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
        auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(0, o2::constants::math::PI);
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (calcRotational)
                hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarProduction);
            }
          } else {
            hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
          }
        } else {
          if (!isMix) {
            if (calcLikeSign)
              hInvMass.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
          }
        }
      } else if (activateTHnSparseCosThStarBeam) {
        beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
        auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(0, o2::constants::math::PI);
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (calcRotational)
                hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarBeam);
            }
          } else {
            hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
          }
        } else {
          if (calcLikeSign)
            hInvMass.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
        }
      } else if (activateTHnSparseCosThStarRandom) {
        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

        randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(0, o2::constants::math::PI);
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (calcRotational)
                hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarRandom);
            }
          } else {
            hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
          }
        } else {
          if (!isMix) {
            if (calcLikeSign)
              hInvMass.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
          }
        }
      }
    }
  }

  // int counter = 0;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)

  {
    rEventSelection.fill(HIST("events_check_data"), 0.5);

    if (cTVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }
    rEventSelection.fill(HIST("events_check_data"), 1.5);

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    rEventSelection.fill(HIST("events_check_data"), 2.5);

    if (!collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("events_check_data"), 3.5);

    multiplicity = collision.centFT0M();

    // Fill the event counter
    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hmult"), multiplicity);
      rEventSelection.fill(HIST("multdist_FT0M"), collision.multFT0M());
      // rEventSelection.fill(HIST("multdist_FT0A"), collision.multFT0A());
      // rEventSelection.fill(HIST("multdist_FT0C"), collision.multFT0C());
      // rEventSelection.fill(HIST("hNcontributor"), collision.numContrib());
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (cQAplots) {
        hPID.fill(HIST("Before/hNsigmaTPC_Ka_before"), track1.pt(), track1.tpcNSigmaKa());
        hPID.fill(HIST("Before/hNsigmaTOF_Ka_before"), track1.pt(), track1.tofNSigmaKa());
        hPID.fill(HIST("Before/hNsigmaTPC_Pi_before"), track2.pt(), track2.tpcNSigmaPi());
        hPID.fill(HIST("Before/hNsigmaTOF_Pi_before"), track2.pt(), track2.tofNSigmaPi());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_before"), track2.tpcNSigmaPi(), track2.tofNSigmaPi());
        hPID.fill(HIST("h1PID_TPC_kaon_data"), track1.tpcNSigmaKa());
        hPID.fill(HIST("h1PID_TPC_pion_data"), track2.tpcNSigmaPi());
        hPID.fill(HIST("h1PID_TOF_kaon_data"), track1.tofNSigmaKa());
        hPID.fill(HIST("h1PID_TOF_pion_data"), track2.tofNSigmaPi());

        hOthers.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());
        hOthers.fill(HIST("hphi"), track1.phi());
      }
      rEventSelection.fill(HIST("events_check_data"), 4.5);

      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionTrack(track2)) {
        continue;
      }
      rEventSelection.fill(HIST("events_check_data"), 5.5);
      // if (counter < 1e4)
      //   std::cout << "TOF beta value is " << track1.beta() << std::endl;
      // counter++;

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy"), track1.dcaXY());
        rEventSelection.fill(HIST("hDcaz"), track1.dcaZ());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if (!selectionPID(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (!selectionPID(track2, 0)) // Track 2 is checked with Pion
        continue;

      rEventSelection.fill(HIST("events_check_data"), 6.5);

      if (cFakeTrack && isFakeTrack(track1, 1)) // Kaon
        continue;
      if (cFakeTrack && isFakeTrack(track2, 0)) // Pion
        continue;

      // if (cMID) {
      //   if (cMIDselectionPID(track1, 0)) // Kaon misidentified as pion
      //     continue;
      //   if (cMIDselectionPID(track1, 2)) // Kaon misidentified as proton
      //     continue;
      //   if (cMIDselectionPID(track2, 1)) // Pion misidentified as kaon
      //     continue;
      // }

      rEventSelection.fill(HIST("events_check_data"), 7.5);

      if (cQAplots) {
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        hPID.fill(HIST("After/hNsigmaKaonTPC_after"), track1.pt(), track1.tpcNSigmaKa());
        hPID.fill(HIST("After/hNsigmaKaonTOF_after"), track1.pt(), track1.tofNSigmaKa());
        hPID.fill(HIST("After/hNsigmaPionTPC_after"), track2.pt(), track2.tpcNSigmaPi());
        hPID.fill(HIST("After/hNsigmaPionTOF_after"), track2.pt(), track2.tofNSigmaPi());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_after"), track1.tpcNSigmaKa(), track1.tofNSigmaKa());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_after"), track2.tpcNSigmaPi(), track2.tofNSigmaPi());
      }

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      rEventSelection.fill(HIST("events_check_data"), 8.5);

      lv3.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
      lv1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), massKa);
      lv2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), massPi);
      lv3 = lv1 + lv2;
      isMix = false;
      fillInvMass(track1, track2, lv2, lv3, multiplicity, isMix);
    }
  }

  PROCESS_SWITCH(Kstarqa, processSE, "Process Same event", true);

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for ME mixing"};
  // ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for ME mixing"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity axis for ME mixing"};

  // using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
  // BinningTypeCentralityM binningOnCentrality{{axisVertex, axisMultiplicity}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair1{binningOnPositions, cfgNoMixedEvents, -1, &cache};
  // SameKindPair<EventCandidates, TrackCandidates, BinningTypeCentralityM> pair2{binningOnCentrality, cfgNoMixedEvents, -1, &cache};

  void processME(EventCandidates const&, TrackCandidates const&)
  {
    for (const auto& [c1, tracks1, c2, tracks2] : pair1) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }

      if (cTVXEvsel && (!c1.selection_bit(aod::evsel::kIsTriggerTVX) || !c2.selection_bit(aod::evsel::kIsTriggerTVX))) {
        return;
      }

      multiplicity = c1.centFT0M();

      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(t1)) // Kaon
          continue;
        if (!selectionTrack(t2)) // Pion
          continue;
        if (!selectionPID(t1, 1)) // Kaon
          continue;
        if (!selectionPID(t2, 0)) // Pion
          continue;

        // if (cMID) {
        //   if (cMIDselectionPID(t1, 0)) // misidentified as pion
        //     continue;
        //   if (cMIDselectionPID(t1, 2)) // misidentified as proton
        //     continue;
        //   if (cMIDselectionPID(t2, 1)) // misidentified as kaon
        //     continue;
        // }

        // TLorentzVector vKAON;
        // vKAON.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
        // TLorentzVector vPION;
        // vPION.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massPi);
        lv3.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
        lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massPi);

        // TLorentzVector kstar = vKAON + vPION;
        lv3 = lv1 + lv2;
        isMix = true;

        // if (std::abs(kstar.Rapidity()) < 0.5) {
        //   fillInvMass(t1, t2, vPION, kstar, multiplicity, isMix);

        if (std::abs(lv3.Rapidity()) < 0.5) {
          fillInvMass(t1, t2, lv2, lv3, multiplicity, isMix);
        }
      }
    }
  }

  PROCESS_SWITCH(Kstarqa, processME, "Process Mixed event", true);

  void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    rEventSelection.fill(HIST("events_check"), 0.5);
    if (std::abs(mcCollision.posZ()) < cutzvertex) {
      rEventSelection.fill(HIST("events_check"), 1.5);
    }

    int nChInel = 0;
    for (const auto& mcParticle : mcParticles) {
      auto pdgcode = std::abs(mcParticle.pdgCode());
      if (mcParticle.isPhysicalPrimary() && (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11 || pdgcode == 13)) {
        if (std::abs(mcParticle.eta()) < 1.0) {
          nChInel = nChInel + 1;
        }
      }
    }
    if (nChInel > 0 && std::abs(mcCollision.posZ()) < cutzvertex)
      rEventSelection.fill(HIST("events_check"), 2.5);

    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;

    multiplicity = 0;
    for (const auto& collision : collisions) {
      // if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cutzvertex) {
      if (std::abs(collision.mcCollision().posZ()) > cutzvertex) {
        continue;
      }

      if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      if (cTVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        continue;
      }
      multiplicity = collision.centFT0M();
      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);
    rEventSelection.fill(HIST("events_check"), 3.5);

    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();

    if (!cAllGenCollisions && !evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    rEventSelection.fill(HIST("events_check"), 4.5);

    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      rEventSelection.fill(HIST("events_check"), 5.5);

      if (std::abs(mcParticle.pdgCode()) != 313) {
        continue;
      }
      rEventSelection.fill(HIST("events_check"), 6.5);

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2) {
        continue;
      }
      rEventSelection.fill(HIST("events_check"), 7.5);

      auto passkaon = false;
      auto passpion = false;
      for (const auto& kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        rEventSelection.fill(HIST("events_check"), 8.5);

        if (std::abs(kCurrentDaughter.pdgCode()) == 321) {
          // if (kCurrentDaughter.pdgCode() == +321) {
          passkaon = true;
          rEventSelection.fill(HIST("events_check"), 9.5);

        } else if (std::abs(kCurrentDaughter.pdgCode()) == 211) {
          //} else if (kCurrentDaughter.pdgCode() == -321) {
          passpion = true;
          // rEventSelection.fill(HIST("events_check"), 10.5);
        }
      }
      if (passkaon && passpion) {
        // if (mcParticle.pdgCode() > 0)
        hInvMass.fill(HIST("hk892GenpT"), mcParticle.pt(), multiplicity);
        // else
        //   hInvMass.fill(HIST("hk892GenpTAnti"), mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processGen, "Process Generated", false);

  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {

    // TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    multiplicity = collision.centFT0M();

    rEventSelection.fill(HIST("events_checkrec"), 0.5);

    if (!collision.has_mcCollision()) {
      return;
    }
    rEventSelection.fill(HIST("events_checkrec"), 1.5);

    // if (std::abs(collision.mcCollision().posZ()) > cutzvertex || !collision.sel8()) {
    if (std::abs(collision.mcCollision().posZ()) > cutzvertex) {
      return;
    }
    rEventSelection.fill(HIST("events_checkrec"), 2.5);

    if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    rEventSelection.fill(HIST("events_checkrec"), 3.5);

    if (cTVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }
    rEventSelection.fill(HIST("events_checkrec"), 4.5);

    auto oldindex = -999;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      rEventSelection.fill(HIST("events_checkrec"), 5.5);

      if (!track1.has_mcParticle()) {
        continue;
      }
      rEventSelection.fill(HIST("events_checkrec"), 6.5);

      auto track1ID = track1.index();
      for (const auto& track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 7.5);

        if (!selectionTrack(track2)) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 8.5);

        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 9.5);

        if (track1.sign() * track2.sign() >= 0) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 10.5);

        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());

        if (cQAplots && track1PDG == 211) {
          hPID.fill(HIST("h1PID_TPC_kaon_MC"), track1.tpcNSigmaKa());
          hPID.fill(HIST("h1PID_TOF_kaon_MC"), track1.tofNSigmaKa());
        }
        if (cQAplots && track1PDG == 321) {
          hPID.fill(HIST("h1PID_TPC_pion_MC"), track1.tpcNSigmaPi());
          hPID.fill(HIST("h1PID_TOF_pion_MC"), track1.tofNSigmaPi());
        }

        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 11.5);

        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 12.5);

        // if (!(track1PDG == 321 && track2PDG == 211)) {
        //   continue;
        // }
        if (!(track1PDG == 211) && !(track1PDG == 321)) {
          continue;
        }
        if (!(track2PDG == 211) && !(track2PDG == 321)) {
          continue;
        }
        rEventSelection.fill(HIST("events_checkrec"), 13.5);

        if (track1PDG == 211) {
          if (!(selectionPID(track1, 0) && selectionPID(track2, 1))) { // pion and kaon
            continue;
          }
        } else {
          if (!(selectionPID(track1, 1) && selectionPID(track2, 0))) { // kaon and pion
            continue;
          }
        }
        rEventSelection.fill(HIST("events_checkrec"), 14.5);

        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            rEventSelection.fill(HIST("events_checkrec"), 15.5);

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            rEventSelection.fill(HIST("events_checkrec"), 16.5);

            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            rEventSelection.fill(HIST("events_checkrec"), 17.5);

            if (std::abs(mothertrack1.y()) >= 0.5) {
              continue;
            }
            rEventSelection.fill(HIST("events_checkrec"), 18.5);

            if (std::abs(mothertrack1.pdgCode()) != 313) {
              continue;
            }

            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              hInvMass.fill(HIST("h1KSRecsplit"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();
            pvec0 = std::array{track1.px(), track1.py(), track1.pz()};
            pvec1 = std::array{track2.px(), track2.py(), track2.pz()};
            auto arrMomrec = std::array{pvec0, pvec1};
            auto motherP = mothertrack1.p();
            auto motherE = mothertrack1.e();
            auto genMass = std::sqrt(motherE * motherE - motherP * motherP);
            auto recMass = RecoDecay::m(arrMomrec, std::array{massKa, massPi});
            auto recpt = std::sqrt((track1.px() + track2.px()) * (track1.px() + track2.px()) + (track1.py() + track2.py()) * (track1.py() + track2.py()));
            //// Resonance reconstruction
            // lDecayDaughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), massKa);
            // lDecayDaughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), massPi);
            // lResonance = lDecayDaughter1 + lDecayDaughter2;

            hInvMass.fill(HIST("h1KstarRecMass"), recMass);
            hInvMass.fill(HIST("h1genmass"), genMass);
            hInvMass.fill(HIST("h2KstarRecpt1"), mothertrack1.pt(), multiplicity);
            hInvMass.fill(HIST("h2KstarRecpt2"), recpt, multiplicity);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processRec, "Process Reconstructed", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Kstarqa>(cfgc)};
}
