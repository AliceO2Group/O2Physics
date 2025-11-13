#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/**
 * @struct PIDFeatureExtractor
 * @brief O2Physics task for extracting particle identification features from AO2D files
 *
 * This task processes track data from the ALICE experiment and extracts comprehensive
 * PID (Particle Identification) features for machine learning applications.
 * It combines TPC and TOF information to compute Bayesian probabilities and saves
 * features to both ROOT TTree and CSV formats.
 */
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
  
  int event_id;      /// Unique identifier for each collision event
  int track_id;      /// Track index within the event
  
  // Momentum components (in GeV/c)
  float px, py, pz;  /// Cartesian momentum components
  float pt, p;       /// Transverse momentum and total momentum
  
  // Angular variables
  float eta;         /// Pseudorapidity
  float phi;         /// Azimuthal angle
  float theta;       /// Polar angle (calculated from eta)
  
  // Track properties
  int charge;        /// Track charge (+1 or -1)
  int track_type;    /// Type of track (e.g., 0=global, 1=TPC-only, etc.)

  // ============================================================================
  // TPC VARIABLES - Time Projection Chamber PID information
  // ============================================================================
  
  float tpc_signal;               /// dE/dx energy loss in TPC (specific ionization)
  
  // n-sigma values: standard deviations from expected energy loss for each particle
  float tpc_nsigma_pi;            /// n-sigma for pion (π)
  float tpc_nsigma_ka;            /// n-sigma for kaon (K)
  float tpc_nsigma_pr;            /// n-sigma for proton (p)
  float tpc_nsigma_el;            /// n-sigma for electron (e)
  
  // Track quality variables
  int tpc_nclusters;              /// Number of TPC clusters used in track fit
  float tpc_chi2;                 /// Chi-square per degree of freedom of TPC fit

  // ============================================================================
  // TOF VARIABLES - Time-Of-Flight PID information
  // ============================================================================
  
  float tof_beta;                 /// β = v/c (velocity over speed of light)
  float tof_mass;                 /// Reconstructed mass from TOF measurement
  
  // n-sigma values for TOF detection
  float tof_nsigma_pi;            /// n-sigma for pion in TOF
  float tof_nsigma_ka;            /// n-sigma for kaon in TOF
  float tof_nsigma_pr;            /// n-sigma for proton in TOF
  float tof_nsigma_el;            /// n-sigma for electron in TOF

  // ============================================================================
  // BAYESIAN PID VARIABLES - Combined PID probabilities
  // ============================================================================
  
  /// Bayesian probability that track is a pion (probability sum = 1.0)
  float bayes_prob_pi;
  /// Bayesian probability that track is a kaon
  float bayes_prob_ka;
  /// Bayesian probability that track is a proton
  float bayes_prob_pr;
  /// Bayesian probability that track is an electron
  float bayes_prob_el;

  // ============================================================================
  // MONTE CARLO TRUTH INFORMATION - For simulated data
  // ============================================================================
  
  int mc_pdg;                     /// PDG code of true particle (0 if no MC match)
  float mc_px, mc_py, mc_pz;      /// True momentum components from simulation

  // ============================================================================
  // DETECTOR AVAILABILITY FLAGS
  // ============================================================================
  
  bool has_tpc;                   /// Flag: track has TPC information
  bool has_tof;                   /// Flag: track has TOF information

  // ============================================================================
  // TRACK IMPACT PARAMETERS - Quality and background rejection
  // ============================================================================
  
  float dca_xy;                   /// Distance of closest approach in xy-plane
  float dca_z;                    /// Distance of closest approach in z-direction

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
  // INITIALIZATION FUNCTION
  // ============================================================================
  
  /**
   * @brief Initialize output files and histograms
   *
   * Called once at task startup. Creates ROOT TTree and CSV file headers,
   * and initializes all quality control histograms.
   */
  void init(InitContext const&) {
    std::string base = outputPath.value;
    
    // ========================================================================
    // ROOT OUTPUT SETUP
    // ========================================================================
    if (exportROOT) {
      // Create ROOT file for storing the TTree
      outputFile = std::make_unique<TFile>((base + ".root").c_str(), "RECREATE");
      
      // Create TTree with descriptive name and title
      featureTree = std::make_unique<TTree>("pid_features", "PID features");

      // Create branches for KINEMATIC VARIABLES
      featureTree->Branch("event_id", &event_id);
      featureTree->Branch("track_id", &track_id);
      featureTree->Branch("px", &px);
      featureTree->Branch("py", &py);
      featureTree->Branch("pz", &pz);
      featureTree->Branch("pt", &pt);
      featureTree->Branch("p", &p);
      featureTree->Branch("eta", &eta);
      featureTree->Branch("phi", &phi);
      featureTree->Branch("theta", &theta);
      featureTree->Branch("charge", &charge);
      featureTree->Branch("track_type", &track_type);

      // Create branches for TPC VARIABLES
      featureTree->Branch("tpc_signal", &tpc_signal);
      featureTree->Branch("tpc_nsigma_pi", &tpc_nsigma_pi);
      featureTree->Branch("tpc_nsigma_ka", &tpc_nsigma_ka);
      featureTree->Branch("tpc_nsigma_pr", &tpc_nsigma_pr);
      featureTree->Branch("tpc_nsigma_el", &tpc_nsigma_el);
      featureTree->Branch("tpc_nclusters", &tpc_nclusters);
      featureTree->Branch("tpc_chi2", &tpc_chi2);

      // Create branches for TOF VARIABLES
      featureTree->Branch("tof_beta", &tof_beta);
      featureTree->Branch("tof_mass", &tof_mass);
      featureTree->Branch("tof_nsigma_pi", &tof_nsigma_pi);
      featureTree->Branch("tof_nsigma_ka", &tof_nsigma_ka);
      featureTree->Branch("tof_nsigma_pr", &tof_nsigma_pr);
      featureTree->Branch("tof_nsigma_el", &tof_nsigma_el);

      // Create branches for BAYESIAN PID VARIABLES
      featureTree->Branch("bayes_prob_pi", &bayes_prob_pi);
      featureTree->Branch("bayes_prob_ka", &bayes_prob_ka);
      featureTree->Branch("bayes_prob_pr", &bayes_prob_pr);
      featureTree->Branch("bayes_prob_el", &bayes_prob_el);

      // Create branches for MONTE CARLO TRUTH (simulated data only)
      featureTree->Branch("mc_pdg", &mc_pdg);
      featureTree->Branch("mc_px", &mc_px);
      featureTree->Branch("mc_py", &mc_py);
      featureTree->Branch("mc_pz", &mc_pz);

      // Create branches for DETECTOR FLAGS
      featureTree->Branch("has_tpc", &has_tpc);
      featureTree->Branch("has_tof", &has_tof);
      
      // Create branches for IMPACT PARAMETERS
      featureTree->Branch("dca_xy", &dca_xy);
      featureTree->Branch("dca_z", &dca_z);
    }

    // ========================================================================
    // CSV OUTPUT SETUP
    // ========================================================================
    if (exportCSV) {
      csvFile.open((base + ".csv").c_str());
      // Write CSV header with all column names
      csvFile <<
        "event_id,track_id,px,py,pz,pt,p,eta,phi,theta,charge,track_type,"
        "tpc_signal,tpc_nsigma_pi,tpc_nsigma_ka,tpc_nsigma_pr,tpc_nsigma_el,"
        "tpc_nclusters,tpc_chi2,"
        "tof_beta,tof_mass,tof_nsigma_pi,tof_nsigma_ka,tof_nsigma_pr,tof_nsigma_el,"
        "bayes_prob_pi,bayes_prob_ka,bayes_prob_pr,bayes_prob_el,"
        "mc_pdg,mc_px,mc_py,mc_pz,has_tpc,has_tof,dca_xy,dca_z\n";
    }

    // ========================================================================
    // HISTOGRAM SETUP - Quality Control Plots
    // ========================================================================
    
    // Define histogram axes with binning
    const AxisSpec axisPt{200, 0, 10, "pT"};               // 200 bins, 0-10 GeV/c
    const AxisSpec axisEta{60, -1.5, 1.5, "eta"};         // 60 bins, -1.5 to 1.5
    const AxisSpec axisdEdx{300, 0, 300, "dE/dx"};        // 300 bins, 0-300
    const AxisSpec axisBeta{120, 0, 1.2, "beta"};         // 120 bins, 0 to 1.2
    const AxisSpec axisMass{100, -0.2, 2.0, "mass"};      // 100 bins, -0.2 to 2.0 GeV/c²

    // Add histograms to registry
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
  
  /**
   * @brief Compute Bayesian probabilities combining TPC and TOF information
   *
   * Uses Gaussian likelihood in n-sigma space and Bayesian inference to combine
   * TPC dE/dx and TOF mass measurements.
   *
   * @param[in] nsTPC[4]  n-sigma values for [pion, kaon, proton, electron] from TPC
   * @param[in] nsTOF[4]  n-sigma values for [pion, kaon, proton, electron] from TOF
   * @param[in] pri[4]    Prior probabilities for each particle hypothesis
   * @param[out] out[4]   Output Bayesian probabilities (normalized to sum=1)
   *
   * Formula: P(particle|TPC,TOF) ∝ P(TPC|particle) * P(TOF|particle) * P(particle)
   *
   * Likelihood: L_i = exp(-0.5 * (ns_TPC_i² + ns_TOF_i²))
   */
  void computeBayesianPID(float nsTPC[4], float nsTOF[4], float pri[4], float out[4]) {
    float sum = 0;
    
    // Calculate likelihood for each particle species
    for (int i = 0; i < 4; i++) {
      // Gaussian likelihood: exp(-0.5 * chi²)
      // Handle invalid TOF values (NaN) by replacing with 0 contribution
      float l = std::exp(-0.5f * (nsTPC[i]*nsTPC[i] + 
                                  (std::isfinite(nsTOF[i]) ? nsTOF[i]*nsTOF[i] : 0.f)));
      
      // Apply prior probability and accumulate
      out[i] = l * pri[i];
      sum += out[i];
    }
    
    // Normalize probabilities so they sum to 1.0
    for (int i = 0; i < 4; i++) {
      out[i] = sum > 0 ? out[i] / sum : 0.f;
    }
  }

  // ============================================================================
  // MAIN PROCESSING FUNCTION
  // ============================================================================
  
  /**
   * @brief Process collision and track data, extract PID features
   *
   * Called for each collision event in the input data. Applies track selections,
   * extracts features from TPC and TOF detectors, computes Bayesian PID,
   * and writes output to ROOT and/or CSV.
   *
   * @param collision  Collision event data
   * @param tracks     Table of tracks with all associated PID information
   * @param mcParticles Monte Carlo particle information (for simulated data)
   */
  void process(
    aod::Collision const& collision,
    soa::Join<
      aod::Tracks,                                          // Base track properties
      aod::TracksExtra,                                     // Extended track info
      aod::TracksDCA,                                       // Impact parameters (DCA)
      aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,         // TPC PID for pion, kaon, proton
      aod::pidTPCEl,                                        // TPC PID for electron
      aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr,         // TOF PID for pion, kaon, proton
      aod::pidTOFEl,                                        // TOF PID for electron
      aod::pidTOFmass, aod::pidTOFbeta,                     // TOF mass and beta
      aod::McTrackLabels                                    // MC truth matching
    > const& tracks,
    aod::McParticles const& mcParticles)
  {
    // Use static counter to maintain event numbering across process calls
    static int eventCounter = 0;
    event_id = eventCounter++;
    int idx = 0;

    // ======================================================================
    // TRACK LOOP - Process each track in the event
    // ======================================================================
    for (auto& t : tracks) {
      
      // ====================================================================
      // TRACK SELECTION - Apply kinematic cuts
      // ====================================================================
      if (t.pt() < ptMin || t.pt() > ptMax) continue;       // Apply pT cut
      if (t.eta() < etaMin || t.eta() > etaMax) continue;   // Apply eta cut
      
      track_id = idx++;
      
      // ====================================================================
      // EXTRACT KINEMATIC VARIABLES
      // ====================================================================
      px = t.px();
      py = t.py();
      pz = t.pz();
      pt = t.pt();
      p = t.p();
      eta = t.eta();
      phi = t.phi();
      // Calculate polar angle from pseudorapidity: θ = 2*arctan(exp(-η))
      theta = 2.f * atanf(expf(-eta));
      charge = t.sign();           // Track charge
      track_type = t.trackType();  // Track categorization

      // ====================================================================
      // EXTRACT TPC INFORMATION
      // ====================================================================
      has_tpc = t.hasTPC();
      if (has_tpc) {
        // TPC has valid measurement
        tpc_signal = t.tpcSignal();           // dE/dx specific ionization
        tpc_nsigma_pi = t.tpcNSigmaPi();      // Deviation from pion hypothesis
        tpc_nsigma_ka = t.tpcNSigmaKa();      // Deviation from kaon hypothesis
        tpc_nsigma_pr = t.tpcNSigmaPr();      // Deviation from proton hypothesis
        tpc_nsigma_el = t.tpcNSigmaEl();      // Deviation from electron hypothesis
        tpc_nclusters = t.tpcNClsFound();     // Quality: number of clusters
        tpc_chi2 = t.tpcChi2NCl();            // Quality: fit chi-square
      } else {
        // TPC has no valid measurement - set sentinel values
        tpc_signal = tpc_nsigma_pi = tpc_nsigma_ka = tpc_nsigma_pr = tpc_nsigma_el = -999;
        tpc_nclusters = 0;
        tpc_chi2 = -999;
      }

      // ====================================================================
      // EXTRACT TOF INFORMATION
      // ====================================================================
      has_tof = t.hasTOF();
      if (has_tof) {
        // TOF has valid measurement
        tof_beta = t.beta();                  // Velocity over c
        tof_mass = t.mass();                  // Reconstructed mass
        tof_nsigma_pi = t.tofNSigmaPi();      // Deviation from pion hypothesis
        tof_nsigma_ka = t.tofNSigmaKa();      // Deviation from kaon hypothesis
        tof_nsigma_pr = t.tofNSigmaPr();      // Deviation from proton hypothesis
        tof_nsigma_el = t.tofNSigmaEl();      // Deviation from electron hypothesis
      } else {
        // TOF has no valid measurement - set sentinel values
        tof_beta = tof_mass = -999;
        tof_nsigma_pi = tof_nsigma_ka = tof_nsigma_pr = tof_nsigma_el = -999;
      }

      // ====================================================================
      // EXTRACT IMPACT PARAMETERS (track quality)
      // ====================================================================
      dca_xy = t.dcaXY();  // Distance of closest approach in transverse plane
      dca_z = t.dcaZ();    // Distance of closest approach along beam axis

      // ====================================================================
      // COMPUTE BAYESIAN PID
      // ====================================================================
      float arrTPC[4] = {tpc_nsigma_pi, tpc_nsigma_ka, tpc_nsigma_pr, tpc_nsigma_el};
      float arrTOF[4] = {tof_nsigma_pi, tof_nsigma_ka, tof_nsigma_pr, tof_nsigma_el};
      float priors[4] = {1.f, 0.2f, 0.1f, 0.05f};  // Prior prob: π, K, p, e
      float probs[4];
      
      // Compute combined PID probabilities
      computeBayesianPID(arrTPC, arrTOF, priors, probs);
      bayes_prob_pi = probs[0];
      bayes_prob_ka = probs[1];
      bayes_prob_pr = probs[2];
      bayes_prob_el = probs[3];

      // ====================================================================
      // EXTRACT MONTE CARLO TRUTH (if available)
      // ====================================================================
      // Safely access MC particle information with existence check
      if (t.has_mcParticle()) {
        auto mc = t.mcParticle();
        mc_pdg = mc.pdgCode();     // Particle identifier code
        mc_px = mc.px();           // True momentum components
        mc_py = mc.py();
        mc_pz = mc.pz();
      } else {
        // No MC match - set sentinel values
        mc_pdg = 0;
        mc_px = mc_py = mc_pz = 0;
      }

      // ====================================================================
      // WRITE OUTPUT
      // ====================================================================
      
      // Write to ROOT TTree
      if (exportROOT) featureTree->Fill();
      
      // Write to CSV file
      if (exportCSV) {
        csvFile << event_id << "," << track_id << ","
                << px << "," << py << "," << pz << ","
                << pt << "," << p << ","
                << eta << "," << phi << "," << theta << ","
                << charge << "," << track_type << ","
                << tpc_signal << "," << tpc_nsigma_pi << "," << tpc_nsigma_ka << "," << tpc_nsigma_pr << "," << tpc_nsigma_el << ","
                << tpc_nclusters << "," << tpc_chi2 << ","
                << tof_beta << "," << tof_mass << "," << tof_nsigma_pi << "," << tof_nsigma_ka << "," << tof_nsigma_pr << "," << tof_nsigma_el << ","
                << bayes_prob_pi << "," << bayes_prob_ka << "," << bayes_prob_pr << "," << bayes_prob_el << ","
                << mc_pdg << "," << mc_px << "," << mc_py << "," << mc_pz << ","
                << has_tpc << "," << has_tof << ","
                << dca_xy << "," << dca_z << "\n";
      }

      // ====================================================================
      // FILL QUALITY CONTROL HISTOGRAMS
      // ====================================================================
      histos.fill(HIST("QC/nTracks"), 1);      // Count total tracks processed
      histos.fill(HIST("QC/pt"), pt);          // pT distribution
      histos.fill(HIST("QC/eta"), eta);        // eta distribution
      
      // TPC dE/dx vs pT (only if TPC measurement exists)
      if (has_tpc) histos.fill(HIST("QC/tpc_dEdx_vs_pt"), pt, tpc_signal);
      
      // TOF beta and mass vs momentum (only if TOF measurement exists)
      if (has_tof) {
        histos.fill(HIST("QC/tof_beta_vs_p"), p, tof_beta);
        histos.fill(HIST("QC/mass_vs_p"), p, tof_mass);
      }
    }
  }

  // ============================================================================
  // FINALIZATION FUNCTION
  // ============================================================================
  
  /**
   * @brief Clean up and finalize output files
   *
   * Called at task completion. Writes TTree to file and closes all output files.
   */
  void finalize() {
    if (exportROOT) {
      // Write TTree to ROOT file and close
      outputFile->cd();
      featureTree->Write();
      outputFile->Close();
    }
    if (exportCSV) {
      // Close CSV file
      csvFile.close();
    }
  }
};

// ============================================================================
// WORKFLOW DEFINITION
// ============================================================================

/**
 * @brief Define the O2Physics workflow
 *
 * This function creates and registers the PIDFeatureExtractor task
 * into the O2 data processing workflow.
 */
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{adaptAnalysisTask<PIDFeatureExtractor>(cfgc)};
}