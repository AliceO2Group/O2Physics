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
/// \brief this is a code for the elliptic flow of identified hadrons
/// \author prottay das, preet
/// \since 29/05/2024

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <iostream>

#include <array>
#include <cmath>
#include <cstdlib>

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
using std::array;

struct v2ellip {

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Confugrable for QA histograms
  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
  bool onlyTPC = true;

  // Configurables for track selections
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPCPi{"nsigmacutTPCPi", 3.0, "Value of the TPC Nsigma cut for pions"};
  Configurable<float> nsigmaCutTPCKa{"nsigmacutTPCKa", 3.0, "Value of the TPC Nsigma cut for kaons"};
  Configurable<float> nsigmaCutTPCPr{"nsigmacutTPCPr", 3.0, "Value of the TPC Nsigma cut for protons"};
  Configurable<float> nsigmaCutTOFPi{"nsigmacutTOFPi", 3.0, "Value of the TOF Nsigma cut for pions"};
  Configurable<float> nsigmaCutTOFKa{"nsigmacutTOFKa", 3.0, "Value of the TOF Nsigma cut for kaons"};
  Configurable<float> nsigmaCutTOFPr{"nsigmacutTOFPr", 3.0, "Value of the TOF Nsigma cut for protons"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};

  // Event selection configurables
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> TVXEvsel{"TVXEvsel", false, "Triggger selection"};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  ConfigurableAxis binsMultPlot{"binsCent", {201, -0.5f, 200.5f}, "Binning of the centrality axis for plots"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm] for plots"};
    AxisSpec axisv2ref = {10, 0, 10, "v2_{ref}"};
    AxisSpec axisv2diff = {14, 0, 14, "v2_{def}"};
    AxisSpec axisphi = {700, 0, 7, "#phi"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Centrality distribution", kTH1F, {{binsMultPlot}});

    // v2 tprofiles for reference and differential flow
    histos.add("profv2ref", "profv2ref", kTProfile, {axisv2ref});
    histos.add("profv2diff_pr_10_20", "profv2diff_pr_10_20", kTProfile, {axisv2diff});
    histos.add("profv2diff_pi_10_20", "profv2diff_pi_10_20", kTProfile, {axisv2diff});
    histos.add("profv2diff_k_10_20", "profv2diff_k_10_20", kTProfile, {axisv2diff});

    histos.add("profv2diff_pr_20_30", "profv2diff_pr_20_30", kTProfile, {axisv2diff});
    histos.add("profv2diff_pi_20_30", "profv2diff_pi_20_30", kTProfile, {axisv2diff});
    histos.add("profv2diff_k_20_30", "profv2diff_k_20_30", kTProfile, {axisv2diff});

    histos.add("profv2diff_pr_30_40", "profv2diff_pr_30_40", kTProfile, {axisv2diff});
    histos.add("profv2diff_pi_30_40", "profv2diff_pi_30_40", kTProfile, {axisv2diff});
    histos.add("profv2diff_k_30_40", "profv2diff_k_30_40", kTProfile, {axisv2diff});

    histos.add("profv2diff_pr_40_50", "profv2diff_pr_40_50", kTProfile, {axisv2diff});
    histos.add("profv2diff_pi_40_50", "profv2diff_pi_40_50", kTProfile, {axisv2diff});
    histos.add("profv2diff_k_40_50", "profv2diff_k_40_50", kTProfile, {axisv2diff});

    // histogram for phi distribution
    histos.add("hphi", "hphi", kTH1F, {axisphi});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (candidate.pt() > 0.4) {
      onlyTPC = false;
    }

    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi) {
          return true;
        }
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa) {
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
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      }
    } else if (PID == 2) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOFPr) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOFPr) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCPr) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCPr) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPCPr) {
          return true;
        }
      }
    }
    return false;
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements

  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelectionExtension>>;

  // Defining partitions for subevents for eta-gap method
  Partition<TrackCandidates> Atracks = (aod::track::eta > 0.4f) && (aod::track::eta < 0.8f);   // partition for subevent A
  Partition<TrackCandidates> Btracks = (aod::track::eta < -0.4f) && (aod::track::eta > -0.8f); // partition for subevent B

  array<float, 15> ptbins = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0};

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)

  {
    float sum_sinA = 0.0, sum_cosA = 0.0, sum_sinB = 0.0, sum_cosB = 0.0;
    int multA = 0, multB = 0;

    // Q vector elements
    array<float, 14> sum_sindsA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // sin component of Q vector for subevent A
    array<float, 14> sum_cosdsA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // cos component of Q vector for subevent A
    array<float, 14> sum_sindsB = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // sin component of Q vector for subevent B
    array<float, 14> sum_cosdsB = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // cos component of Q vector for subevent B

    // p vector definitions for subevent A
    array<float, 14> pn_sumsinA_pr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for proton
    array<float, 14> pn_sumcosA_pr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for proton

    array<float, 14> pn_sumsinA_pi = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for pion
    array<float, 14> pn_sumcosA_pi = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for pion

    array<float, 14> pn_sumsinA_k = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for kaon
    array<float, 14> pn_sumcosA_k = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for kaon

    // p vector definitions for subevent B
    array<float, 14> pn_sumsinB_pr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for proton
    array<float, 14> pn_sumcosB_pr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for proton

    array<float, 14> pn_sumsinB_pi = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for pion
    array<float, 14> pn_sumcosB_pi = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for pion

    array<float, 14> pn_sumsinB_k = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for kaon
    array<float, 14> pn_sumcosB_k = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // to store ptwise pn vector components for kaon

    // POI multiplicities
    array<int, 14> mpA_pr = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // proton multiplicity for subevent A
    array<int, 14> mpB_pr = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // proton multiplicity for subevent B

    array<int, 14> mpA_pi = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // pion multiplicity for subevent A
    array<int, 14> mpB_pi = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // pion multiplicity for subevent B

    array<int, 14> mpA_k = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // kaon multiplicity for subevent A
    array<int, 14> mpB_k = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // kaon multiplicity for subevent B

    if (!collision.sel8()) {
      return;
    }

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }

    float multiplicity = 0.0f;
    multiplicity = collision.centFT0C();

    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    rEventSelection.fill(HIST("hmult"), multiplicity);

    auto atrack = Atracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto btrack = Btracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto track : tracks) {
      if (!selectionTrack(track)) {
        continue;
      }
      if (selectionPID(track, 0) || selectionPID(track, 1) || selectionPID(track, 2)) { // If track pion, kaon or proton
        histos.fill(HIST("hphi"), track.phi());
      } else {
        continue;
      }
    } // end of track loop

    for (auto track1 : atrack) {
      if (!selectionTrack(track1)) {
        continue;
      }

      sum_sinA += TMath::Sin(2.0 * track1.phi()); // sum of sin components of Q vector
      sum_cosA += TMath::Cos(2.0 * track1.phi()); // sum of cos components of Q vector
      multA++;                                    // charged particle multiplicity

      if (selectionPID(track1, 0) || selectionPID(track1, 1) || selectionPID(track1, 2)) { // If track  pion, kaon or proton
        // pt loop for component sums of p vector, POI multiplicities pt wise
        for (auto pt = 0; pt < 14; pt++) {
          sum_sindsA[pt] += TMath::Sin(2 * track1.phi());
          sum_cosdsA[pt] += TMath::Cos(2 * track1.phi());

          if (track1.pt() > ptbins[pt] && track1.pt() <= ptbins[pt + 1] && selectionPID(track1, 0)) { // for pion
            pn_sumsinA_pi[pt] += TMath::Sin(2 * track1.phi());
            pn_sumcosA_pi[pt] += TMath::Cos(2 * track1.phi());
            mpA_pi[pt]++;
          } else if (track1.pt() > ptbins[pt] && track1.pt() <= ptbins[pt + 1] && selectionPID(track1, 1)) { // for kaon
            pn_sumsinA_k[pt] += TMath::Sin(2 * track1.phi());
            pn_sumcosA_k[pt] += TMath::Cos(2 * track1.phi());
            mpA_k[pt]++;
          } else if (track1.pt() > ptbins[pt] && track1.pt() <= ptbins[pt + 1] && selectionPID(track1, 2)) { // for proton
            pn_sumsinA_pr[pt] += TMath::Sin(2 * track1.phi());
            pn_sumcosA_pr[pt] += TMath::Cos(2 * track1.phi());
            mpA_pr[pt]++;
          } else {
            continue;
          }
        } // end of pt loop
      } else {
        continue;
      }
    } // track loop ends

    for (auto track2 : btrack) {
      if (!selectionTrack(track2)) {
        continue;
      }

      sum_sinB += TMath::Sin(2.0 * track2.phi()); // sum of sin components of Q vector
      sum_cosB += TMath::Cos(2.0 * track2.phi()); // sum of cos components of Q vector
      multB++;                                    // charged particle multiplicity

      if (selectionPID(track2, 0) || selectionPID(track2, 1) || selectionPID(track2, 2)) { // If track  pion, kaon or proton
        // pt loop for component sums of p vector, POI multiplicities pt wise
        for (auto pt = 0; pt < 14; pt++) {
          sum_sindsB[pt] += TMath::Sin(2 * track2.phi());
          sum_cosdsB[pt] += TMath::Cos(2 * track2.phi());

          if (track2.pt() > ptbins[pt] && track2.pt() <= ptbins[pt + 1] && selectionPID(track2, 0)) { // for pion
            pn_sumsinB_pi[pt] += TMath::Sin(2 * track2.phi());
            pn_sumcosB_pi[pt] += TMath::Cos(2 * track2.phi());
            mpB_pi[pt]++;
          } else if (track2.pt() > ptbins[pt] && track2.pt() <= ptbins[pt + 1] && selectionPID(track2, 1)) { // for kaon
            pn_sumsinB_k[pt] += TMath::Sin(2 * track2.phi());
            pn_sumcosB_k[pt] += TMath::Cos(2 * track2.phi());
            mpB_k[pt]++;
          } else if (track2.pt() > ptbins[pt] && track2.pt() <= ptbins[pt + 1] && selectionPID(track2, 2)) { // for proton
            pn_sumsinB_pr[pt] += TMath::Sin(2 * track2.phi());
            pn_sumcosB_pr[pt] += TMath::Cos(2 * track2.phi());
            mpB_pr[pt]++;
          } else {
            continue;
          }
        } // end of pt loop
      } else {
        continue;
      }
    } // track loop ends

    if (10.0 < multiplicity && multiplicity <= 20.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 1, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }

      // pt wise differential flow
      for (auto pt = 0; pt < 14; pt++) {
        if ((mpA_pr[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pr_10_20"), pt + 1, ((pn_sumcosA_pr[pt] * sum_cosB + pn_sumsinA_pr[pt] * sum_sinB) / (mpA_pr[pt] * multB)), mpA_pr[pt] * multB);
        } // for proton
        if ((mpA_pi[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pi_10_20"), pt + 1, ((pn_sumcosA_pi[pt] * sum_cosB + pn_sumsinA_pi[pt] * sum_sinB) / (mpA_pi[pt] * multB)), mpA_pi[pt] * multB);
        } // for pion
        if ((mpA_k[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_k_10_20"), pt + 1, ((pn_sumcosA_k[pt] * sum_cosB + pn_sumsinA_k[pt] * sum_sinB) / (mpA_k[pt] * multB)), mpA_k[pt] * multB);
        } // for kaon
      }
    } // 10 to 20 percent centrality

    if (20.0 < multiplicity && multiplicity <= 30.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 2, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }

      // pt wise differential flow
      for (auto pt = 0; pt < 14; pt++) {
        if ((mpA_pr[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pr_20_30"), pt + 1, ((pn_sumcosA_pr[pt] * sum_cosB + pn_sumsinA_pr[pt] * sum_sinB) / (mpA_pr[pt] * multB)), mpA_pr[pt] * multB);
        } // for proton
        if ((mpA_pi[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pi_20_30"), pt + 1, ((pn_sumcosA_pi[pt] * sum_cosB + pn_sumsinA_pi[pt] * sum_sinB) / (mpA_pi[pt] * multB)), mpA_pi[pt] * multB);
        } // for pion
        if ((mpA_k[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_k_20_30"), pt + 1, ((pn_sumcosA_k[pt] * sum_cosB + pn_sumsinA_k[pt] * sum_sinB) / (mpA_k[pt] * multB)), mpA_k[pt] * multB);
        } // for kaon
      }
    } // 20 to 30 percent centrality

    if (30.0 < multiplicity && multiplicity <= 40.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 3, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }

      // pt wise differential flow
      for (auto pt = 0; pt < 14; pt++) {
        if ((mpA_pr[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pr_30_40"), pt + 1, ((pn_sumcosA_pr[pt] * sum_cosB + pn_sumsinA_pr[pt] * sum_sinB) / (mpA_pr[pt] * multB)), mpA_pr[pt] * multB);
        } // for proton
        if ((mpA_pi[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pi_30_40"), pt + 1, ((pn_sumcosA_pi[pt] * sum_cosB + pn_sumsinA_pi[pt] * sum_sinB) / (mpA_pi[pt] * multB)), mpA_pi[pt] * multB);
        } // for pion
        if ((mpA_k[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_k_30_40"), pt + 1, ((pn_sumcosA_k[pt] * sum_cosB + pn_sumsinA_k[pt] * sum_sinB) / (mpA_k[pt] * multB)), mpA_k[pt] * multB);
        } // for kaon
      }
    } // 30 to 40 percent centrality

    if (40.0 < multiplicity && multiplicity <= 50.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 4, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }

      // pt wise differential flow
      for (auto pt = 0; pt < 14; pt++) {
        if ((mpA_pr[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pr_40_50"), pt + 1, ((pn_sumcosA_pr[pt] * sum_cosB + pn_sumsinA_pr[pt] * sum_sinB) / (mpA_pr[pt] * multB)), mpA_pr[pt] * multB);
        } // for proton
        if ((mpA_pi[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_pi_40_50"), pt + 1, ((pn_sumcosA_pi[pt] * sum_cosB + pn_sumsinA_pi[pt] * sum_sinB) / (mpA_pi[pt] * multB)), mpA_pi[pt] * multB);
        } // for pion
        if ((mpA_k[pt] * multB) != 0) {
          histos.fill(HIST("profv2diff_k_40_50"), pt + 1, ((pn_sumcosA_k[pt] * sum_cosB + pn_sumsinA_k[pt] * sum_sinB) / (mpA_k[pt] * multB)), mpA_k[pt] * multB);
        } // for kaon
      }
    } // 40 to 50 percent centrality

    if (50.0 < multiplicity && multiplicity <= 60.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 5, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }
    } // 50 to 60 percent centrality

    if (60.0 < multiplicity && multiplicity <= 70.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 6, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }
    } // 60 to 70 percent centrality

    if (70.0 < multiplicity && multiplicity <= 80.0) {
      // reference flow
      if ((multA * multB) != 0) {
        histos.fill(HIST("profv2ref"), 7, ((sum_cosA * sum_cosB + sum_sinA * sum_sinB) / (multA * multB)), multA * multB);
      }
    } // 70 to 80 percent centrality

  } // end of process

  PROCESS_SWITCH(v2ellip, processSE, "Process Same event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<v2ellip>(cfgc)};
}
