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

/// \file pidFeatureExtractor.cxx
/// \brief Task to extract particle identification features from ALICE AO2D data for machine learning workflows
/// \author Robert Forynski

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <memory>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// PIDFeatureExtractor task for extracting particle identification features from AO2D files
struct PIDFeatureExtractor {
  // ============================================================================
  // OUTPUT OBJECTS - File and data structures for feature storage
  // ============================================================================
  /// Output ROOT file for storing the TTree with extracted features
  std::unique_ptr<TFile> outputFile;

  /// TTree storing all extracted features for each track
  std::unique_ptr<TTree> featureTree;

  /// CSV output stream for exporting features in comma-separated format
  std::ofstream csvFile;

  // ============================================================================
  // KINEMATIC VARIABLES - Track momentum and position information
  // ============================================================================
  int eventId;      /// Unique identifier for each collision event
  int trackId;      /// Track index within the event

  // Momentum components (in GeV/c)
  float px, py, pz; /// Cartesian momentum components
  float pt, p;      /// Transverse momentum and total momentum

  // Angular variables
  float eta;   /// Pseudorapidity
  float phi;   /// Azimuthal angle
  float theta; /// Polar angle (calculated from eta)

  // Track properties
  int charge;    /// Track charge (+1 or -1)
  int trackType; /// Type of track (e.g., 0=global, 1=TPC-only, etc.)

  // ============================================================================
  // TPC VARIABLES - Time Projection Chamber PID information
  // ============================================================================
  float tpcSignal; /// dE/dx energy loss in TPC (specific ionization)

  // n-sigma values: standard deviations from expected energy loss for each particle
  float tpcNsigmaPi; /// n-sigma for pion (π)
  float tpcNsigmaKa; /// n-sigma for kaon (K)
  float tpcNsigmaPr; /// n-sigma for proton (p)
  float tpcNsigmaEl; /// n-sigma for electron (e)

  // Track quality variables
  int tpcNclusters; /// Number of TPC clusters used in track fit
  float tpcChi2;    /// Chi-square per degree of freedom of TPC fit

  // ============================================================================
  // TOF VARIABLES - Time-Of-Flight PID information
  // ============================================================================
  float tofBeta; /// β = v/c (velocity over speed of light)
  float tofMass; /// Reconstructed mass from TOF measurement

  // n-sigma values for TOF detection
  float tofNsigmaPi; /// n-sigma for pion in TOF
  float tofNsigmaKa; /// n-sigma for kaon in TOF
  float tofNsigmaPr; /// n-sigma for proton in TOF
  float tofNsigmaEl; /// n-sigma for electron in TOF

  // ============================================================================
  // BAYESIAN PID VARIABLES - Combined PID probabilities
  // ============================================================================
  /// Bayesian probability that track is a pion (probability sum = 1.0)
  float bayesProbPi;
  /// Bayesian probability that track is a kaon
  float bayesProbKa;
  /// Bayesian probability that track is a proton
  float bayesProbPr;
  /// Bayesian probability that track is an electron
  float bayesProbEl;

  // ============================================================================
  // MONTE CARLO TRUTH INFORMATION - For simulated data
  // ============================================================================
  int mcPdg;              /// PDG code of true particle (0 if no MC match)
  float mcPx, mcPy, mcPz; /// True momentum components from simulation

  // ============================================================================
  // DETECTOR AVAILABILITY FLAGS
  // ============================================================================
  bool hasTpc; /// Flag: track has TPC information
  bool hasTof; /// Flag: track has TOF information

  // ============================================================================
  // TRACK IMPACT PARAMETERS - Quality and background rejection
  // ============================================================================
  float dcaXy; /// Distance of closest approach in xy-plane
  float dcaZ;  /// Distance of closest approach in z-direction

  // ============================================================================
  // HISTOGRAM REGISTRY - Quality control histograms
  // ============================================================================
  /// Registry for quality control histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ============================================================================
  // CONFIGURABLE PARAMETERS - User-adjustable settings
  // ============================================================================
  /// Base path and filename for output files (without extension)
  Configurable<std::string> outputPath{"outputPath", "pid_features", "Output file base"};

  /// Enable CSV export of features
  Configurable<bool> exportCSV{"exportCSV", true, "Export CSV"};

  /// Enable ROOT file export of features
  Configurable<bool> exportROOT{"exportROOT", true, "Export ROOT"};

  /// Minimum pseudorapidity cut for track selection
  Configurable<float> etaMin{"etaMin", -1.5f, "Minimum eta"};

  /// Maximum pseudorapidity cut for track selection
  Configurable<float> etaMax{"etaMax", 1.5f, "Maximum eta"};

  /// Minimum transverse momentum cut (GeV/c)
  Configurable<float> ptMin{"ptMin", 0.1f, "Minimum pT"};

  /// Maximum transverse momentum cut (GeV/c)
  Configurable<float> ptMax{"ptMax", 20.0f, "Maximum pT"};

  // ============================================================================
  // CONSTANTS
  // ============================================================================
  static constexpr int KNumSpecies = 4;
  static constexpr float KPriorPi = 1.0f;
  static constexpr float KPriorKa = 0.2f;
  static constexpr float KPriorPr = 0.1f;
  static constexpr float KPriorEl = 0.05f;
  static constexpr float KSentinelValue = -999.0f;

  // ============================================================================
  // INITIALIZATION FUNCTION
  // ============================================================================
  /// Initialize output files and histograms
  void init(InitContext const&)
  {
    std::string base = outputPath.value;

    // ROOT OUTPUT SETUP
    if (exportROOT) {
      outputFile = std::make_unique<TFile>((base + ".root").c_str(), "RECREATE");
      featureTree = std::make_unique<TTree>("pid_features", "PID features");

      // KINEMATIC VARIABLES
      featureTree->Branch("eventId", &eventId);
      featureTree->Branch("trackId", &trackId);
      featureTree->Branch("px", &px);
      featureTree->Branch("py", &py);
      featureTree->Branch("pz", &pz);
      featureTree->Branch("pt", &pt);
      featureTree->Branch("p", &p);
      featureTree->Branch("eta", &eta);
      featureTree->Branch("phi", &phi);
      featureTree->Branch("theta", &theta);
      featureTree->Branch("charge", &charge);
      featureTree->Branch("trackType", &trackType);

      // TPC VARIABLES
      featureTree->Branch("tpcSignal", &tpcSignal);
      featureTree->Branch("tpcNsigmaPi", &tpcNsigmaPi);
      featureTree->Branch("tpcNsigmaKa", &tpcNsigmaKa);
      featureTree->Branch("tpcNsigmaPr", &tpcNsigmaPr);
      featureTree->Branch("tpcNsigmaEl", &tpcNsigmaEl);
      featureTree->Branch("tpcNclusters", &tpcNclusters);
      featureTree->Branch("tpcChi2", &tpcChi2);

      // TOF VARIABLES
      featureTree->Branch("tofBeta", &tofBeta);
      featureTree->Branch("tofMass", &tofMass);
      featureTree->Branch("tofNsigmaPi", &tofNsigmaPi);
      featureTree->Branch("tofNsigmaKa", &tofNsigmaKa);
      featureTree->Branch("tofNsigmaPr", &tofNsigmaPr);
      featureTree->Branch("tofNsigmaEl", &tofNsigmaEl);

      // BAYESIAN PID VARIABLES
      featureTree->Branch("bayesProbPi", &bayesProbPi);
      featureTree->Branch("bayesProbKa", &bayesProbKa);
      featureTree->Branch("bayesProbPr", &bayesProbPr);
      featureTree->Branch("bayesProbEl", &bayesProbEl);

      // MONTE CARLO TRUTH
      featureTree->Branch("mcPdg", &mcPdg);
      featureTree->Branch("mcPx", &mcPx);
      featureTree->Branch("mcPy", &mcPy);
      featureTree->Branch("mcPz", &mcPz);

      // DETECTOR FLAGS
      featureTree->Branch("hasTpc", &hasTpc);
      featureTree->Branch("hasTof", &hasTof);

      // IMPACT PARAMETERS
      featureTree->Branch("dcaXy", &dcaXy);
      featureTree->Branch("dcaZ", &dcaZ);
    }

    // CSV OUTPUT SETUP
    if (exportCSV) {
      csvFile.open((base + ".csv").c_str());
      csvFile << "eventId,trackId,px,py,pz,pt,p,eta,phi,theta,charge,trackType,"
                 "tpcSignal,tpcNsigmaPi,tpcNsigmaKa,tpcNsigmaPr,tpcNsigmaEl,"
                 "tpcNclusters,tpcChi2,"
                 "tofBeta,tofMass,tofNsigmaPi,tofNsigmaKa,tofNsigmaPr,tofNsigmaEl,"
                 "bayesProbPi,bayesProbKa,bayesProbPr,bayesProbEl,"
                 "mcPdg,mcPx,mcPy,mcPz,hasTpc,hasTof,dcaXy,dcaZ\n";
    }

    // HISTOGRAM SETUP
    const AxisSpec axisPt{200, 0, 10, "pT"};
    const AxisSpec axisEta{60, -1.5, 1.5, "eta"};
    const AxisSpec axisdEdx{300, 0, 300, "dE/dx"};
    const AxisSpec axisBeta{120, 0, 1.2, "beta"};
    const AxisSpec axisMass{100, -0.2, 2.0, "mass"};

    histos.add("QC/nTracks", "Tracks", kTH1F, {{10000, 0, 100000}});
    histos.add("QC/pt", "pT", kTH1F, {axisPt});
    histos.add("QC/eta", "eta", kTH1F, {axisEta});
    histos.add("QC/tpc_dEdx_vs_pt", "dE/dx vs pT", kTH2F, {axisPt, axisdEdx});
    histos.add("QC/tof_beta_vs_p", "beta vs p", kTH2F, {axisPt, axisBeta});
    histos.add("QC/mass_vs_p", "mass vs p", kTH2F, {axisPt, axisMass});
  }

  // ============================================================================
  // BAYESIAN PID CALCULATION FUNCTION
  // ============================================================================
  /// Compute Bayesian probabilities combining TPC and TOF information
  void computeBayesianPID(const float nsTPC[KNumSpecies], const float nsTOF[KNumSpecies], const float pri[KNumSpecies], float out[KNumSpecies])
  {
    float sum = 0;

    for (int i = 0; i < KNumSpecies; i++) {
      float l = std::exp(-0.5f * (nsTPC[i] * nsTPC[i] +
                                  (std::isfinite(nsTOF[i]) ? nsTOF[i] * nsTOF[i] : 0.0f)));

      out[i] = l * pri[i];
      sum += out[i];
    }

    for (int i = 0; i < KNumSpecies; i++) {
      out[i] = sum > 0 ? out[i] / sum : 0.0f;
    }
  }

  // ============================================================================
  // MAIN PROCESSING FUNCTION
  // ============================================================================
  /// Process collision and track data, extract PID features
  void process(
    aod::Collision const& collision,
    soa::Join<
      aod::Tracks,
      aod::TracksExtra,
      aod::TracksDCA,
      aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
      aod::pidTPCEl,
      aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr,
      aod::pidTOFEl,
      aod::pidTOFmass, aod::pidTOFbeta,
      aod::McTrackLabels> const& tracks,
    aod::McParticles const& mcParticles)
  {
    static int eventCounter = 0;
    eventId = eventCounter++;
    int idx = 0;

    for (const auto& t : tracks) {
      if (t.pt() < ptMin || t.pt() > ptMax)
        continue;
      if (t.eta() < etaMin || t.eta() > etaMax)
        continue;

      trackId = idx++;

      // Kinematics
      px = t.px();
      py = t.py();
      pz = t.pz();
      pt = t.pt();
      p = t.p();
      eta = t.eta();
      phi = t.phi();
      theta = 2.0f * std::atanf(std::expf(-eta));
      charge = t.sign();
      trackType = t.trackType();

      // TPC info
      hasTpc = t.hasTPC();
      if (hasTpc) {
        tpcSignal = t.tpcSignal();
        tpcNsigmaPi = t.tpcNSigmaPi();
        tpcNsigmaKa = t.tpcNSigmaKa();
        tpcNsigmaPr = t.tpcNSigmaPr();
        tpcNsigmaEl = t.tpcNSigmaEl();
        tpcNclusters = t.tpcNClsFound();
        tpcChi2 = t.tpcChi2NCl();
      } else {
        tpcSignal = tpcNsigmaPi = tpcNsigmaKa = tpcNsigmaPr = tpcNsigmaEl = KSentinelValue;
        tpcNclusters = 0;
        tpcChi2 = KSentinelValue;
      }

      // TOF info
      hasTof = t.hasTOF();
      if (hasTof) {
        tofBeta = t.beta();
        tofMass = t.mass();
        tofNsigmaPi = t.tofNSigmaPi();
        tofNsigmaKa = t.tofNSigmaKa();
        tofNsigmaPr = t.tofNSigmaPr();
        tofNsigmaEl = t.tofNSigmaEl();
      } else {
        tofBeta = tofMass = KSentinelValue;
        tofNsigmaPi = tofNsigmaKa = tofNsigmaPr = tofNsigmaEl = KSentinelValue;
      }

      // Impact parameters
      dcaXy = t.dcaXY();
      dcaZ = t.dcaZ();

      // Bayesian PID calculation
      float arrTPC[KNumSpecies] = {tpcNsigmaPi, tpcNsigmaKa, tpcNsigmaPr, tpcNsigmaEl};
      float arrTOF[KNumSpecies] = {tofNsigmaPi, tofNsigmaKa, tofNsigmaPr, tofNsigmaEl};
      float priors[KNumSpecies] = {KPriorPi, KPriorKa, KPriorPr, KPriorEl};
      float probs[KNumSpecies];

      computeBayesianPID(arrTPC, arrTOF, priors, probs);
      bayesProbPi = probs[0];
      bayesProbKa = probs[1];
      bayesProbPr = probs[2];
      bayesProbEl = probs[3];

      // MC truth
      if (t.has_mcParticle()) {
        auto mc = t.mcParticle();
        mcPdg = mc.pdgCode();
        mcPx = mc.px();
        mcPy = mc.py();
        mcPz = mc.pz();
      } else {
        mcPdg = 0;
        mcPx = mcPy = mcPz = 0;
      }

      // Write outputs
      if (exportROOT)
        featureTree->Fill();
      if (exportCSV) {
        csvFile << eventId << "," << trackId << ","
                << px << "," << py << "," << pz << ","
                << pt << "," << p << ","
                << eta << "," << phi << "," << theta << ","
                << charge << "," << trackType << ","
                << tpcSignal << "," << tpcNsigmaPi << "," << tpcNsigmaKa << "," << tpcNsigmaPr << "," << tpcNsigmaEl << ","
                << tpcNclusters << "," << tpcChi2 << ","
                << tofBeta << "," << tofMass << "," << tofNsigmaPi << "," << tofNsigmaKa << "," << tofNsigmaPr << "," << tofNsigmaEl << ","
                << bayesProbPi << "," << bayesProbKa << "," << bayesProbPr << "," << bayesProbEl << ","
                << mcPdg << "," << mcPx << "," << mcPy << "," << mcPz << ","
                << hasTpc << "," << hasTof << ","
                << dcaXy << "," << dcaZ << "\n";
      }

      // Fill QC histograms
      histos.fill(HIST("QC/nTracks"), 1);
      histos.fill(HIST("QC/pt"), pt);
      histos.fill(HIST("QC/eta"), eta);
      if (hasTpc)
        histos.fill(HIST("QC/tpc_dEdx_vs_pt"), pt, tpcSignal);
      if (hasTof) {
        histos.fill(HIST("QC/tof_beta_vs_p"), p, tofBeta);
        histos.fill(HIST("QC/mass_vs_p"), p, tofMass);
      }
    }
  }

  // ============================================================================
  // FINALIZATION FUNCTION
  // ============================================================================
  /// Clean up and finalize output files
  void finalize()
  {
    if (exportROOT) {
      outputFile->cd();
      featureTree->Write();
      outputFile->Close();
    }
    if (exportCSV) {
      csvFile.close();
    }
  }
};

// ============================================================================
// WORKFLOW DEFINITION
// ============================================================================
/// Define the O2Physics workflow
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PIDFeatureExtractor>(cfgc)};
}
