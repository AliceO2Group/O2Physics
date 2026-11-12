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

// Fragmentation function task for Jpsi mesons in Jets.
//
/// \author Andreas Ulvseth Mikalsen <Andreas.ulveseth.mikalsen@cern.ch>
//
#include "PWGJE/Core/JetCandidateUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DqFragmentationFunctionTask {
  // Jet Pt bins for multidimension histograms, not affected by masterswitch unlike other configurables.
  // The first variable is not a bin edge, but a flag for variable width bins, the rest are bin edges. LEAVE AT 0 IN CONFIG!
  ConfigurableAxis jetPtBins{"jetPtBins", {VARIABLE_WIDTH, 0., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 25., 30., 35., 40.}, "Jet pT bin edges (GeV/c), used for multidimensional histograms"};

  // Master switch
  Configurable<bool> useCuts{"useCuts", true, "Enable or disable: eventselection, jetselection and dielectron selection cuts"};
  Configurable<bool> useMCPCuts{"useMCPCuts", false, "Enable or disable (truth only): eventselection, jetselection and dielectron selection cuts"};

  // Event selection
  Configurable<std::string> eventSelections{"eventSelections", "NoTimeFrameBorder+NoITSROFrameBorder", "List of available event cuts JetDerivedDataUtilities.h line 113"};
  Configurable<bool> skipMBgapEvents{"skipMBgapEvents", true, "Skip MB gap events"};
  Configurable<bool> applyRCTflag{"applyRCTflag", true, "Apply RCT flag"};
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT label to use for event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.f, "Maximum absolute primary vertex z (cm)"};

  // Electron kinematics
  Configurable<float> electronPMin{"electronPMin", 1.f, "Minimum electron momentum (GeV/c)"};
  Configurable<float> electronPairPtMin{"electronPairPtMin", 0.f, "Minimum electron pair transverse momentum (GeV/c)"};
  Configurable<float> electronEtaMax{"electronEtaMax", 0.8f, "Maximum absolute electron pseudorapidity"};
  Configurable<float> electronPairEtaMax{"electronPairEtaMax", 0.8f, "Maximum absolute electron pair pseudorapidity"};

  // Electron TPC PID
  Configurable<float> nSigmaElLowPMin{"nSigmaElLowPMin", -2.5f, "Minimum electron n-sigma for p below threshold"};
  Configurable<float> nSigmaElHighPMin{"nSigmaElHighPMin", -1.5f, "Minimum electron n-sigma for p above threshold"};
  Configurable<float> nSigmaElMax{"nSigmaElMax", 3.f, "Maximum electron n-sigma"};
  Configurable<float> electronPIDMomentumThreshold{"electronPIDMomentumThreshold", 5.f, "Momentum threshold for electron PID (GeV/c)"};

  // Pion and proton rejection (low momentum only)
  Configurable<float> pionNSigmaAbsMin{"pionNSigmaAbsMin", 2.5f, "Minimum absolute pion TPC n-sigma"};
  Configurable<float> protonNSigmaAbsMin{"protonNSigmaAbsMin", 3.f, "Minimum absolute proton TPC n-sigma"};

  // Track DCA
  Configurable<float> dcaXYMax{"dcaXYMax", 0.5f, "Maximum absolute DCAxy (cm)"};
  Configurable<float> dcaZMax{"dcaZMax", 0.5f, "Maximum absolute DCAz (cm)"};

  // Track quality
  Configurable<float> tpcChi2Max{"tpcChi2Max", 4.f, "Maximum TPC chi2 per cluster"};
  Configurable<int> tpcClustersMin{"tpcClustersMin", 70, "Minimum number of TPC clusters"};
  Configurable<int> tpcClustersMax{"tpcClustersMax", 161, "Maximum number of TPC clusters"};

  // ITS inner-barrel selection
  Configurable<int> nITSHitsFirst3{"nITSHitsFirst3", 1, "Minimum number of ITS inner-barrel hits (layers 0-2)"};
  Configurable<int> nITSHitsAll7{"nITSHitsAll7", 1, "Minimum number of ITS inner-barrel hits (layers 0-6)"};
  Configurable<bool> requireITSFirstLayer{"requireITSFirstLayer", false, "Require a hit in the first ITS layer (bit 0=1)"};

  // Jet selection
  Configurable<float> jetPtMin{"jetPtMin", 0.f, "Minimum jet transverse momentum (GeV/c)"};
  Configurable<float> jetPtMax{"jetPtMax", 1.e6f, "Maximum jet transverse momentum (GeV/c)"};
  Configurable<float> jetEtaMax{"jetEtaMax", 1.e6f, "Maximum absolute jet pseudorapidity"};

  HistogramRegistry registry{"registry",
                             {{"Data/Jet/h_dielectron_jet_pt", ";p_{T,jet} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"Data/Jet/h_constituent_multiplicity_vs_jet_pt", "Jet constituent multiplicity vs p_{T};p_{T,jet} (GeV/c);N_{constituents}", {HistType::kTH2F, {{200, 0., 80.}, {101, -0.5, 100.5}}}},
                              {"Data/Jet/h_sum_constituent_pt_vs_jet_pt", "Scalar constituent p_{T} sum vs jet p_{T};p_{T,jet} (GeV/c);#Sigma p_{T,constituent} (GeV/c)", {HistType::kTH2F, {{100, 0., 40.}, {100, 0., 40.}}}},
                              {"Data/Jet/h_dielectron_jet_eta", ";#eta_{jet};Counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"Data/Jet/h_dielectron_jet_phi", ";#phi_{jet};Counts", {HistType::kTH1F, {{128, 0., o2::constants::math::TwoPI}}}},
                              {"Data/Jet/h_dielectron_jet_Fragmentation_z", ";z^{e^{+}e^{-},jet};dN/dz^{e^{+}e^{-},jet}", {HistType::kTH1F, {{100, 0., 1.}}}},
                              {"Data/Jet/h_dielectron_jet_distance", ";#DeltaR_{e^{+}e^{-},jet};Counts", {HistType::kTH1F, {{100, 0., 1.}}}},
                              {"Data/Jet/h_dielectron_jet_distance_vs_projection", ";#DeltaR_{e^{+}e^{-},jet};z^{e^{+}e^{-},jet}", {HistType::kTH2F, {{100, 0., 1.}, {100, 0., 1.}}}},
                              {"Data/Jet/h_electron_pairs_per_jet", "Number of accepted electron pairs per jet;Electron pairs per jet;Jets", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
                              {"Data/Dielectron/h_dielectron_pt", ";p_{T,e^{+}e^{-}};dN/dp_{T,e^{+}e^{-}}", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"Data/Dielectron/h_dielectron_pt_vs_jet_pt", "Dielectron p_{T} vs jet p_{T};p_{T,jet} (GeV/c);p_{T,e^{+}e^{-}} (GeV/c)", {HistType::kTH2F, {{200, 0., 80.}, {200, 0., 80.}}}},
                              {"Data/Dielectron/h_dielectron_eta", ";#eta_{e^{+}e^{-}};Counts", {HistType::kTH1F, {{160, -0.8, 0.8}}}},
                              {"Data/Dielectron/h_dielectron_phi", ";#phi_{e^{+}e^{-}};Counts", {HistType::kTH1F, {{128, 0., o2::constants::math::TwoPI}}}},
                              {"Data/Dielectron/h_dielectron_mass", ";m_{e^{+}e^{-}} (GeV/c^{2});dN/dm_{e^{+}e^{-}}", {HistType::kTH1F, {{1000, 0., 10.}}}},
                              {"Data/Dielectron/h_dielectron_dcaxy", "Electron daughter DCA_{xy};DCA_{xy} (cm);Tracks", {HistType::kTH1F, {{200, -1.5, 1.5}}}},
                              {"Data/Dielectron/h_dielectron_dcaz", "Electron daughter DCA_{z};DCA_{z} (cm);Tracks", {HistType::kTH1F, {{200, -1.5, 1.5}}}},
                              {"Data/Dielectron/h_dielectron_Lxy", ";L_{xy};dN/dL_{xy}", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"Data/Dielectron/h_dielectron_Lz", ";L_{z};dN/dL_{z}", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"Data/Dielectron/h_dielectron_Tau_xy", ";#tau_{xy};dN/d#tau_{xy}", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"Data/Dielectron/h_dielectron_Tau_z", ";#tau_{z};dN/d#tau_{z}", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"Data/SingleTrack/h_electron_pt", "Electron daughter p_{T};p_{T,e} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 20.}}}},
                              {"Data/SingleTrack/h_electron_eta", "Electron daughter #eta;#eta_{e};Counts", {HistType::kTH1F, {{160, -0.8, 0.8}}}},
                              {"Data/SingleTrack/h_electron_phi", "Electron daughter #phi;#phi_{e};Counts", {HistType::kTH1F, {{128, 0., o2::constants::math::TwoPI}}}},
                              {"Data/SingleTrack/h_electron_tpc_nsigma", "Electron daughter TPC PID;TPC n#sigma_{e};Counts", {HistType::kTH1F, {{160, -8., 8.}}}},
                              {"Data/SingleTrack/h_electron_tpc_nsigma_vs_pt", "Electron daughter TPC PID;p_{T,e} (GeV/c);TPC n#sigma_{e}", {HistType::kTH2F, {{200, 0., 20.}, {160, -8., 8.}}}},
                              {"Data/SingleTrack/h_electron_tpc_ncls_found", "Electron daughter TPC clusters found;N_{TPC clusters};Tracks", {HistType::kTH1F, {{162, -0.5, 161.5}}}},
                              {"Data/SingleTrack/h_electron_tpc_chi2_per_cluster", "Electron daughter TPC #chi^{2}/cluster;TPC #chi^{2}/cluster;Tracks", {HistType::kTH1F, {{200, 0., 10.}}}},
                              {"Data/SingleTrack/h_electron_its_chi2_per_cluster", "Electron daughter ITS #chi^{2}/cluster;ITS #chi^{2}/cluster;Tracks", {HistType::kTH1F, {{200, 0., 10.}}}},
                              {"Data/SingleTrack/h_jet_track_pt", "Jet constituent p_{T};p_{T,track} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"Data/SingleTrack/h_jet_track_eta", "Jet constituent #eta;#eta_{track};Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
                              {"Data/SingleTrack/h_n_tracks_per_jet", "Charged constituents per jet;N_{constituents};Counts", {HistType::kTH1F, {{101, -0.5, 100.5}}}},

                              {"Matching/h_collision_counter", "Collision counter;Stage;Counts", {HistType::kTH1F, {{4, 0., 4.}}}},
                              {"Matching/h_jet_counter", "Jet counter;Stage;Counts", {HistType::kTH1F, {{4, 0., 4.}}}},
                              {"Matching/h_candidate_counter", "Candidate matching stages;Stage;Counts", {HistType::kTH1F, {{4, 0., 4.}}}},
                              {"Matching/h_2D_prompt_z", "Prompt z comparison;z^{MCD};z^{MCP}", {HistType::kTH2F, {{20, 0., 1.}, {20, 0., 1.}}}},
                              {"Matching/h_2D_prompt_jetpt", "Prompt jetpt;p_{T,jet}^{MCD};p_{T,jet}^{MCP};", {HistType::kTH2F, {{40, 0., 40.}, {40, 0., 40.}}}},
                              {"Matching/h_2D_prompt_jpsipt", "Prompt J/#psi pT comparison;p_{T,J/#psi}^{MCD} (GeV/c);p_{T,J/#psi}^{MCP} (GeV/c)", {HistType::kTH2F, {{40, 0., 40.}, {40, 0., 40.}}}},
                              {"Matching/h_2D_prompt_tauxy", "Prompt dielectron #tau_{xy};#tau_{xy}^{MCD};#tau_{xy}^{MCP}", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
                              {"Matching/h_2D_nonprompt_z", "Nonprompt z comparison;z^{MCD};z^{MCP}", {HistType::kTH2F, {{20, 0., 1.}, {20, 0., 1.}}}},
                              {"Matching/h_2D_nonprompt_jetpt", "Nonprompt jetpt;p_{T,jet}^{MCD};p_{T,jet}^{MCP}", {HistType::kTH2F, {{40, 0., 40.}, {40, 0., 40.}}}},
                              {"Matching/h_2D_nonprompt_jpsipt", "Nonprompt J/#psi pt;p_{T,J/#psi}^{MCD};p_{T,J/#psi}^{MCP}", {HistType::kTH2F, {{40, 0., 40.}, {40, 0., 40.}}}},
                              {"Matching/h_2D_nonprompt_tauxy", "Nonprompt dielectron #tau_{xy};#tau_{xy}^{MCD};#tau_{xy}^{MCP}", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},

                              {"MCD/Jet/h_MCD_jet_pt", "MCD jet p_{T};p_{T,jet} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"MCD/Jet/h_constituent_multiplicity_vs_jet_pt", "MCD jet constituent multiplicity vs p_{T};p_{T,jet} (GeV/c);N_{constituents}", {HistType::kTH2F, {{200, 0., 80.}, {101, -0.5, 100.5}}}},
                              {"MCD/Jet/h_sum_constituent_pt_vs_jet_pt", "MCD scalar constituent p_{T} sum vs jet p_{T};p_{T,jet} (GeV/c);#Sigma p_{T,constituent} (GeV/c)", {HistType::kTH2F, {{100, 0., 40.}, {100, 0., 40.}}}},
                              {"MCD/Jet/h_MCD_jet_eta", "MCD jet #eta;#eta_{jet};Counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"MCD/Jet/h_MCD_jet_phi", "MCD jet #phi;#phi_{jet};Counts", {HistType::kTH1F, {{128, 0., o2::constants::math::TwoPI}}}},
                              {"MCD/Jet/h_MCD_z", "MCD z;z;Counts", {HistType::kTH1F, {{100, 0., 1.}}}},
                              {"MCD/Jet/h_MCD_jet_distance", "MCD J/#psi-jet distance;#DeltaR_{e^{+}e^{-},jet};Counts", {HistType::kTH1F, {{100, 0., 1.}}}},
                              {"MCD/Jet/h_MCD_jet_distance_vs_projection", ";#DeltaR_{e^{+}e^{-},jet};z^{e^{+}e^{-},jet}", {HistType::kTH2F, {{100, 0., 1.}, {100, 0., 1.}}}},
                              {"MCD/Jet/h_electron_pairs_per_jet", "Number of accepted electron pairs per jet;Electron pairs per jet;Jets", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_pt", "MCD J/#psi p_{T};p_{T,J/#psi} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"MCD/Dielectron/h_dielectron_pt_vs_jet_pt", "MCD dielectron p_{T} vs jet p_{T};p_{T,jet} (GeV/c);p_{T,e^{+}e^{-}} (GeV/c)", {HistType::kTH2F, {{200, 0., 80.}, {200, 0., 80.}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_eta", "MCD dielectron #eta;#eta_{ee};Counts", {HistType::kTH1F, {{160, -0.8, 0.8}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_mass", "MCD J/#psi mass;m_{ee} (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0., 10.}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_dcaxy", "MCD electron daughter DCA_{xy};DCA_{xy} (cm);Tracks", {HistType::kTH1F, {{200, -1.5, 1.5}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_dcaz", "MCD electron daughter DCA_{z};DCA_{z} (cm);Tracks", {HistType::kTH1F, {{200, -1.5, 1.5}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_origin", "MCD J/#psi origin;origin;Counts", {HistType::kTH1F, {{3, 0., 3.}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_Lxy", "MCD J/#psi L_{xy};L_{xy} (cm);Counts", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_Lz", "MCD J/#psi L_{z};L_{z};Counts", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_tau_xy", "MCD J/#psi #tau_{xy};#tau_{xy};Counts", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"MCD/Dielectron/h_MCD_dielectron_tau_z", "MCD J/#psi #tau_{z};#tau_{z};Counts", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"MCD/SingleTrack/h_MCD_electron_pt", "MCD electron daughter p_{T};p_{T,e} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 20.}}}},
                              {"MCD/SingleTrack/h_MCD_electron_eta", "MCD electron daughter #eta;#eta_{e};Counts", {HistType::kTH1F, {{160, -0.8, 0.8}}}},
                              {"MCD/SingleTrack/h_MCD_electron_phi", "MCD electron daughter #phi;#phi_{e};Counts", {HistType::kTH1F, {{128, 0., o2::constants::math::TwoPI}}}},
                              {"MCD/SingleTrack/h_MCD_electron_tpc_nsigma", "MCD electron daughter TPC PID;TPC n#sigma_{e};Counts", {HistType::kTH1F, {{160, -8., 8.}}}},
                              {"MCD/SingleTrack/h_MCD_electron_tpc_nsigma_vs_pt", "MCD electron daughter TPC PID;p_{T,e} (GeV/c);TPC n#sigma_{e}", {HistType::kTH2F, {{200, 0., 20.}, {160, -8., 8.}}}},
                              {"MCD/SingleTrack/h_MCD_electron_tpc_ncls_found", "MCD electron daughter TPC clusters found;N_{TPC clusters};Tracks", {HistType::kTH1F, {{162, -0.5, 161.5}}}},
                              {"MCD/SingleTrack/h_MCD_electron_tpc_chi2_per_cluster", "MCD electron daughter TPC #chi^{2}/cluster;TPC #chi^{2}/cluster;Tracks", {HistType::kTH1F, {{200, 0., 10.}}}},
                              {"MCD/SingleTrack/h_MCD_electron_its_chi2_per_cluster", "MCD electron daughter ITS #chi^{2}/cluster;ITS #chi^{2}/cluster;Tracks", {HistType::kTH1F, {{200, 0., 10.}}}},
                              {"MCD/SingleTrack/h_MCD_jet_track_pt", "MCD jet constituent p_{T};p_{T,track} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"MCD/SingleTrack/h_MCD_jet_track_eta", "MCD jet constituent #eta;#eta_{track};Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
                              {"MCD/SingleTrack/h_MCD_n_tracks_per_jet", "MCD charged constituents per jet;N_{constituents};Counts", {HistType::kTH1F, {{101, -0.5, 100.5}}}},

                              {"MCP/Jet/h_MCP_jet_pt", "MCP jet p_{T};p_{T,jet} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"MCP/Jet/h_constituent_multiplicity_vs_jet_pt", "MCP jet constituent multiplicity vs p_{T};p_{T,jet} (GeV/c);N_{constituents}", {HistType::kTH2F, {{200, 0., 80.}, {101, -0.5, 100.5}}}},
                              {"MCP/Jet/h_sum_constituent_pt_vs_jet_pt", "MCP scalar constituent p_{T} sum vs jet p_{T};p_{T,jet} (GeV/c);#Sigma p_{T,constituent} (GeV/c)", {HistType::kTH2F, {{100, 0., 40.}, {100, 0., 40.}}}},
                              {"MCP/Jet/h_MCP_jet_eta", "MCP jet #eta;#eta_{jet};Counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"MCP/Jet/h_MCP_jet_phi", "MCP jet #phi;#phi_{jet};Counts", {HistType::kTH1F, {{128, 0., o2::constants::math::TwoPI}}}},
                              {"MCP/Jet/h_MCP_z", "MCP z;z;Counts", {HistType::kTH1F, {{100, 0., 1.}}}},
                              {"MCP/Jet/h_MCP_jet_distance", "MCP J/#psi-jet distance;#DeltaR_{e^{+}e^{-},jet};Counts", {HistType::kTH1F, {{100, 0., 1.}}}},
                              {"MCP/Jet/h_MCP_jet_distance_vs_projection", ";#DeltaR_{e^{+}e^{-},jet};z^{e^{+}e^{-},jet}", {HistType::kTH2F, {{100, 0., 1.}, {100, 0., 1.}}}},
                              {"MCP/Jet/h_electron_pairs_per_jet", "Number of accepted electron pairs per jet;Electron pairs per jet;Jets", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
                              {"MCP/Dielectron/h_MCP_jpsi_pt", "MCP J/#psi p_{T};p_{T,J/#psi} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"MCP/Dielectron/h_dielectron_pt_vs_jet_pt", "MCP dielectron p_{T} vs jet p_{T};p_{T,jet} (GeV/c);p_{T,e^{+}e^{-}} (GeV/c)", {HistType::kTH2F, {{200, 0., 80.}, {200, 0., 80.}}}},
                              {"MCP/Dielectron/h_MCP_dielectron_eta", "MCP dielectron #eta;#eta_{ee};Counts", {HistType::kTH1F, {{160, -0.8, 0.8}}}},
                              {"MCP/Dielectron/h_MCP_dielectron_mass", "MCP dielectron mass;m_{ee} (GeV/c^{2});Counts", {HistType::kTH1F, {{1000, 0., 10.}}}},
                              {"MCP/Dielectron/h_MCP_dielectron_origin", "MCP J/#psi origin;origin;Counts", {HistType::kTH1F, {{3, 0., 3.}}}},
                              {"MCP/Dielectron/h_MCP_jpsi_lxy", "MCP L_{xy};L_{xy} (cm);Counts", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"MCP/Dielectron/h_MCP_jpsi_tauxy", "MCP #tau_{xy};#tau_{xy};Counts", {HistType::kTH1F, {{1000, -2., 2.}}}},
                              {"MCP/SingleTrack/h_MCP_jet_particle_pt_ElectronPairInfoMissing", "MCP jet constituent p_{T};p_{T,particle} (GeV/c);Counts", {HistType::kTH1F, {{200, 0., 80.}}}},
                              {"MCP/SingleTrack/h_MCP_jet_particle_eta_ElectronPairInfoMissing", "MCP jet constituent #eta;#eta_{particle};Counts", {HistType::kTH1F, {{180, -0.9, 0.9}}}},
                              {"MCP/SingleTrack/h_MCP_n_particles_per_jet", "MCP constituents per jet;N_{constituents};Counts", {HistType::kTH1F, {{101, -0.5, 100.5}}}}}};

  std::vector<int> eventSelectionBits;

  void init(InitContext const&)
  {
    // Add histograms using configurable jet pT binning, placed in init so jetPtBins are applied correctly.
    registry.add("Data/Multidimensional/h_4d_histogram", "Dielectron fragmentation;p_{T,jet} (GeV/c);z;m_{ee} (GeV/c^{2});#tau_{xy}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {150, 2.5, 3.7}, {250, -1., 1.}});
    registry.add("Matching/h_response_prompt", "Prompt J/#psi response;p_{T,jet}^{MCD};z^{MCD};p_{T,jet}^{MCP};z^{MCP}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, AxisSpec{jetPtBins}, {20, 0., 1.}});
    registry.add("Matching/h_response_nonprompt", "Nonprompt J/#psi response;p_{T,jet}^{MCD};z^{MCD};p_{T,jet}^{MCP};z^{MCP}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, AxisSpec{jetPtBins}, {20, 0., 1.}});
    registry.add("MCD/Multidimensional/h_MCD_5d_allParticles", "MCD J/#psi fragmentation;p_{T,jet};z;m_{ee};Tau_{xy};origin", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {150, 2.5, 3.7}, {250, -1., 1.}, {3, 0., 3.}});
    registry.add("MCD/Multidimensional/h_MCD_4d_PromptJpsi", "MCD J/#psi fragmentation promptJpsi;p_{T,jet};z;m_{ee};Tau_{xy}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {150, 2.5, 3.7}, {250, -1., 1.}});
    registry.add("MCD/Multidimensional/h_MCD_4d_NonPromptJpsi", "MCD J/#psi fragmentation NonPromptJpsi;p_{T,jet};z;m_{ee};Tau_{xy}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {150, 2.5, 3.7}, {250, -1., 1.}});
    registry.add("MCP/Multidimensional/h_MCP_4d_allParticles", "MCP J/#psi fragmentation;p_{T,jet};z;Tau_{xy};origin", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {250, -1., 1.}, {3, 0., 3.}});
    registry.add("MCP/Multidimensional/h_MCP_4d_PromptJpsi", "MCP J/#psi fragmentation promptJpsi;p_{T,jet};z;m_{ee};Tau_{xy}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {150, 2.5, 3.7}, {250, -1., 1.}});
    registry.add("MCP/Multidimensional/h_MCP_4d_NonPromptJpsi", "MCP J/#psi fragmentation NonPromptJpsi;p_{T,jet};z;m_{ee};Tau_{xy}", HistType::kTHnF, {AxisSpec{jetPtBins}, {20, 0., 1.}, {150, 2.5, 3.7}, {250, -1., 1.}});

    // Initialize event selection
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(
      static_cast<std::string>(eventSelections));

    // Set matching counter labels
    auto collisionCounter = registry.get<TH1>(HIST("Matching/h_collision_counter"));
    collisionCounter->GetXaxis()->SetBinLabel(1, "MCP collisions");
    collisionCounter->GetXaxis()->SetBinLabel(2, "MCP collisions after collision cuts (if enabled)");
    collisionCounter->GetXaxis()->SetBinLabel(3, "matched MCDcollisions to MCP collisions");
    collisionCounter->GetXaxis()->SetBinLabel(4, "matched MCD collisions after cut (if enabled)");

    auto jetCounter = registry.get<TH1>(HIST("Matching/h_jet_counter"));
    jetCounter->GetXaxis()->SetBinLabel(1, "MCD jets");
    jetCounter->GetXaxis()->SetBinLabel(2, "MCD jets after jet cuts");
    jetCounter->GetXaxis()->SetBinLabel(3, "MCD jets matched to MCP");
    jetCounter->GetXaxis()->SetBinLabel(4, "MCD jets matched to MCP after MCP cut (if enabled)");

    auto candidateCounter = registry.get<TH1>(HIST("Matching/h_candidate_counter"));
    candidateCounter->GetXaxis()->SetBinLabel(1, "MCD candidates in matched jets");
    candidateCounter->GetXaxis()->SetBinLabel(2, "Candidates passing z and optional e cuts");
    candidateCounter->GetXaxis()->SetBinLabel(3, "MCD candidates matched to MCP");
    candidateCounter->GetXaxis()->SetBinLabel(4, "MCD candidates matched to MCP after MCP cut (if enabled)");
  }

  using DielectronMCDJets = soa::Join<
    aod::DielectronChargedMCDetectorLevelJets,
    aod::DielectronChargedMCDetectorLevelJetConstituents,
    aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>;

  using DielectronMCPJets = soa::Join<
    aod::DielectronChargedMCParticleLevelJets,
    aod::DielectronChargedMCParticleLevelJetConstituents,
    aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>;

  // Preslice detector-level collisions belonging to each MC collision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  // Preslice detector-level jets belonging to each detector-level collision
  Preslice<DielectronMCDJets> mcdJetsPerCollisionPreslice = aod::jet::collisionId;

  template <typename JetType, typename JpsiType>
  double calculateZ(
    JetType const& jet,
    JpsiType const& jpsi)
  {
    double z = (jet.pt() > 0.) ? (jpsi.pt()) / jet.pt() : -1.;
    // avoiding z=1.0 overflowing into next bin.
    if (z > 0.999) {
      z = 0.9991;
    }
    return z;
  }

  template <typename Candidate>
  bool passesElectronCuts(const Candidate& candidate, int electronIndex)
  {
    // electronIndex = 1: first electron
    // electronIndex = 2: second electron
    const double pt = electronIndex == 1
                        ? candidate.pt1()
                        : candidate.pt2();

    const double eta = electronIndex == 1
                         ? candidate.eta1()
                         : candidate.eta2();

    // Momentum magnitude: p = pT * cosh(eta) (because dqtrigger 2025 uses p not pt)
    const double p = pt * std::cosh(eta);

    const double nSigmaEl = electronIndex == 1
                              ? candidate.tpcNSigmaEl1()
                              : candidate.tpcNSigmaEl2();

    const double tpcChi2NCl = electronIndex == 1
                                ? candidate.tpcChi2NCl1()
                                : candidate.tpcChi2NCl2();

    const int tpcNClsFound = electronIndex == 1
                               ? candidate.tpcNClsFound1()
                               : candidate.tpcNClsFound2();

    const double dcaXY = electronIndex == 1
                           ? candidate.dcaXY1()
                           : candidate.dcaXY2();

    const double dcaZ = electronIndex == 1
                          ? candidate.dcaZ1()
                          : candidate.dcaZ2();

    const double nSigmaPi = electronIndex == 1
                              ? candidate.tpcNSigmaPi1()
                              : candidate.tpcNSigmaPi2();

    const double nSigmaPr = electronIndex == 1
                              ? candidate.tpcNSigmaPr1()
                              : candidate.tpcNSigmaPr2();

    const unsigned char itsClusterMap = electronIndex == 1
                                          ? candidate.itsClusterMap1()
                                          : candidate.itsClusterMap2();

    // Electron kinematic selection
    if (p < electronPMin.value ||
        std::abs(eta) > electronEtaMax.value) {
      return false;
    }

    // Momentum-dependent electron TPC PID
    const double nSigmaElMin =
      p < electronPIDMomentumThreshold.value
        ? nSigmaElLowPMin.value
        : nSigmaElHighPMin.value;

    if (nSigmaEl < nSigmaElMin ||
        nSigmaEl > nSigmaElMax.value) {
      return false;
    }

    // Pion and proton rejection: low momentum only
    if (p < electronPIDMomentumThreshold.value) {

      if (std::abs(nSigmaPi) < pionNSigmaAbsMin.value) {
        return false;
      }

      if (std::abs(nSigmaPr) < protonNSigmaAbsMin.value) {
        return false;
      }
    }

    // TPC track quality
    if (tpcChi2NCl > tpcChi2Max.value) {
      return false;
    }

    if (tpcNClsFound < tpcClustersMin.value ||
        tpcNClsFound > tpcClustersMax.value) {
      return false;
    }

    // Distance of closest approach
    if (std::abs(dcaXY) > dcaXYMax.value) {
      return false;
    }

    if (std::abs(dcaZ) > dcaZMax.value) {
      return false;
    }

    // Hits in ITS
    // Count hits in the first three ITS layers (bits 0, 1, 2).
    const int nITSIBHits =
      ((itsClusterMap & (1u << 0)) != 0) +
      ((itsClusterMap & (1u << 1)) != 0) +
      ((itsClusterMap & (1u << 2)) != 0);

    if (nITSIBHits < nITSHitsFirst3.value) {
      return false;
    }
    // Count hits in the all 7 ITS layers (bits 0, 1, 2, 3, 4, 5, 6).
    const int nITSIBHitsAllLayers =
      ((itsClusterMap & (1u << 0)) != 0) +
      ((itsClusterMap & (1u << 1)) != 0) +
      ((itsClusterMap & (1u << 2)) != 0) +
      ((itsClusterMap & (1u << 3)) != 0) +
      ((itsClusterMap & (1u << 4)) != 0) +
      ((itsClusterMap & (1u << 5)) != 0) +
      ((itsClusterMap & (1u << 6)) != 0);

    if (nITSIBHitsAllLayers < nITSHitsAll7.value) {
      return false;
    }

    // Require a hit in the first ITS layer (bit 0).
    if (requireITSFirstLayer.value &&
        (itsClusterMap & (1u << 0)) == 0) {
      return false;
    }

    // All cuts passed
    return true;
  }
  template <typename Candidate>
  bool passesPairCuts(const Candidate& candidate)
  {
    if (candidate.pt() < electronPairPtMin.value ||
        std::abs(candidate.eta()) > electronPairEtaMax.value) {
      return false;
    }

    return true;
  }

  void FragmentationFunctionData(aod::JetCollision const& collision,
                                 soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents> const& jets,
                                 aod::CandidatesDielectronData const&,
                                 aod::JetTracks const&)
  {
    if (useCuts && (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBgapEvents.value, applyRCTflag.value, rctLabel.value) || !(std::abs(collision.posZ()) < vertexZCut))) {
      return;
    }
    for (const auto& jet : jets) {

      if (useCuts) {
        if (jet.pt() < jetPtMin || jet.pt() > jetPtMax ||
            std::abs(jet.eta()) > jetEtaMax) {
          continue;
        }
      }

      const auto jetTracks = jet.template tracks_as<aod::JetTracks>();
      const auto dielectronCandidates = jet.template candidates_as<aod::CandidatesDielectronData>();
      double sumConstituentPt = 0.;
      int electronPairs = 0;
      for (const auto& dielectronCandidate : dielectronCandidates) {
        if (useCuts &&
            (!passesPairCuts(dielectronCandidate) ||
             !passesElectronCuts(dielectronCandidate, 1) ||
             !passesElectronCuts(dielectronCandidate, 2))) {
          continue;
        }

        const double z = calculateZ(jet, dielectronCandidate);
        if (z < 0.) {
          continue;
        }
        ++electronPairs;
        sumConstituentPt += dielectronCandidate.pt1();
        sumConstituentPt += dielectronCandidate.pt2();

        // calculating angular distance in eta-phi plane
        const double axisDistance =
          jetutilities::deltaR(jet, dielectronCandidate);

        // filling candidate histograms
        registry.fill(HIST("Data/SingleTrack/h_jet_track_pt"), dielectronCandidate.pt1());
        registry.fill(HIST("Data/SingleTrack/h_jet_track_pt"), dielectronCandidate.pt2());
        registry.fill(HIST("Data/SingleTrack/h_jet_track_eta"), dielectronCandidate.eta1());
        registry.fill(HIST("Data/SingleTrack/h_jet_track_eta"), dielectronCandidate.eta2());
        registry.fill(HIST("Data/Jet/h_dielectron_jet_Fragmentation_z"), z);
        registry.fill(HIST("Data/Jet/h_dielectron_jet_distance_vs_projection"), axisDistance, z);
        registry.fill(HIST("Data/Jet/h_dielectron_jet_distance"), axisDistance);

        registry.fill(HIST("Data/Dielectron/h_dielectron_mass"), dielectronCandidate.mass());
        registry.fill(HIST("Data/Dielectron/h_dielectron_dcaxy"), dielectronCandidate.dcaXY1());
        registry.fill(HIST("Data/Dielectron/h_dielectron_dcaxy"), dielectronCandidate.dcaXY2());
        registry.fill(HIST("Data/Dielectron/h_dielectron_dcaz"), dielectronCandidate.dcaZ1());
        registry.fill(HIST("Data/Dielectron/h_dielectron_dcaz"), dielectronCandidate.dcaZ2());
        registry.fill(HIST("Data/Dielectron/h_dielectron_pt"), dielectronCandidate.pt());
        registry.fill(HIST("Data/Dielectron/h_dielectron_pt_vs_jet_pt"), jet.pt(), dielectronCandidate.pt());
        registry.fill(HIST("Data/Dielectron/h_dielectron_eta"), dielectronCandidate.eta());
        registry.fill(HIST("Data/Dielectron/h_dielectron_phi"), dielectronCandidate.phi());
        registry.fill(HIST("Data/Dielectron/h_dielectron_Lz"), dielectronCandidate.lz());
        registry.fill(HIST("Data/Dielectron/h_dielectron_Lxy"), dielectronCandidate.lxy());
        registry.fill(HIST("Data/Dielectron/h_dielectron_Tau_z"), dielectronCandidate.tauz());
        registry.fill(HIST("Data/Dielectron/h_dielectron_Tau_xy"), dielectronCandidate.tauxy());
        registry.fill(HIST("Data/SingleTrack/h_electron_pt"), dielectronCandidate.pt1());
        registry.fill(HIST("Data/SingleTrack/h_electron_pt"), dielectronCandidate.pt2());
        registry.fill(HIST("Data/SingleTrack/h_electron_eta"), dielectronCandidate.eta1());
        registry.fill(HIST("Data/SingleTrack/h_electron_eta"), dielectronCandidate.eta2());
        registry.fill(HIST("Data/SingleTrack/h_electron_phi"), dielectronCandidate.phi1());
        registry.fill(HIST("Data/SingleTrack/h_electron_phi"), dielectronCandidate.phi2());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_nsigma"), dielectronCandidate.tpcNSigmaEl1());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_nsigma"), dielectronCandidate.tpcNSigmaEl2());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_nsigma_vs_pt"), dielectronCandidate.pt1(), dielectronCandidate.tpcNSigmaEl1());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_nsigma_vs_pt"), dielectronCandidate.pt2(), dielectronCandidate.tpcNSigmaEl2());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_ncls_found"), dielectronCandidate.tpcNClsFound1());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_ncls_found"), dielectronCandidate.tpcNClsFound2());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_chi2_per_cluster"), dielectronCandidate.tpcChi2NCl1());
        registry.fill(HIST("Data/SingleTrack/h_electron_tpc_chi2_per_cluster"), dielectronCandidate.tpcChi2NCl2());
        registry.fill(HIST("Data/SingleTrack/h_electron_its_chi2_per_cluster"), dielectronCandidate.itsChi2NCl1());
        registry.fill(HIST("Data/SingleTrack/h_electron_its_chi2_per_cluster"), dielectronCandidate.itsChi2NCl2());
        registry.fill(HIST("Data/Multidimensional/h_4d_histogram"), jet.pt(), z, dielectronCandidate.mass(), dielectronCandidate.tauxy());
      } // end of loop over dielectron candidates
      registry.fill(HIST("Data/Jet/h_electron_pairs_per_jet"), electronPairs);
      if (electronPairs == 0) {
        continue;
      }
      for (const auto& track : jetTracks) {
        sumConstituentPt += track.pt();
        registry.fill(HIST("Data/SingleTrack/h_jet_track_pt"), track.pt());
        registry.fill(HIST("Data/SingleTrack/h_jet_track_eta"), track.eta());
      }
      registry.fill(HIST("Data/Jet/h_dielectron_jet_pt"), jet.pt());
      registry.fill(HIST("Data/Jet/h_dielectron_jet_eta"), jet.eta());
      registry.fill(HIST("Data/Jet/h_dielectron_jet_phi"), jet.phi());
      const int constituentMultiplicity = static_cast<int>(jetTracks.size()) + 2 * electronPairs; // adding 2 for the dielectron candidate
      registry.fill(HIST("Data/Jet/h_constituent_multiplicity_vs_jet_pt"), jet.pt(), constituentMultiplicity);
      registry.fill(HIST("Data/Jet/h_sum_constituent_pt_vs_jet_pt"), jet.pt(), sumConstituentPt);
      registry.fill(HIST("Data/SingleTrack/h_n_tracks_per_jet"), constituentMultiplicity);
    } // end of loop over jets
  }

  PROCESS_SWITCH(DqFragmentationFunctionTask, FragmentationFunctionData, "fragmentation function data", false);

  void FragmentationFunctionMCD(
    aod::JetCollisionMCD const& collision,
    soa::Join<aod::DielectronChargedMCDetectorLevelJets,
              aod::DielectronChargedMCDetectorLevelJetConstituents> const& jets,
    aod::CandidatesDielectronMCD const&,
    aod::JetTracksMCD const&)
  {
    // Event selection
    if (useCuts && (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBgapEvents.value, applyRCTflag.value, rctLabel.value) || !(std::abs(collision.posZ()) < vertexZCut))) {
      return;
    }

    // Loop over detector-level jets
    for (const auto& jet : jets) {
      // Optional jet selection
      if (useCuts) {
        if (jet.pt() < jetPtMin || jet.pt() > jetPtMax ||
            std::abs(jet.eta()) > jetEtaMax) {
          continue;
        }
      }
      const auto jetTracks = jet.template tracks_as<aod::JetTracksMCD>();
      const auto dielectronCandidates = jet.template candidates_as<aod::CandidatesDielectronMCD>();
      double sumConstituentPt = 0.;
      int electronPairs = 0;
      for (const auto& dielectronCandidate : dielectronCandidates) {
        if (useCuts &&
            (!passesPairCuts(dielectronCandidate) ||
             !passesElectronCuts(dielectronCandidate, 1) ||
             !passesElectronCuts(dielectronCandidate, 2))) {
          continue;
        }
        const double z = calculateZ(jet, dielectronCandidate);
        if (z < 0.) {
          continue;
        }
        sumConstituentPt += dielectronCandidate.pt1();
        sumConstituentPt += dielectronCandidate.pt2();
        ++electronPairs;

        const double axisDistance = jetutilities::deltaR(jet, dielectronCandidate);

        // Fill candidate histograms
        registry.fill(HIST("MCD/SingleTrack/h_MCD_jet_track_pt"), dielectronCandidate.pt1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_jet_track_pt"), dielectronCandidate.pt2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_jet_track_eta"), dielectronCandidate.eta1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_jet_track_eta"), dielectronCandidate.eta2());
        registry.fill(HIST("MCD/Jet/h_MCD_z"), z);
        registry.fill(HIST("MCD/Jet/h_MCD_jet_distance"), axisDistance);
        registry.fill(HIST("MCD/Jet/h_MCD_jet_distance_vs_projection"), axisDistance, z);
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_mass"), dielectronCandidate.mass());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_dcaxy"), dielectronCandidate.dcaXY1());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_dcaxy"), dielectronCandidate.dcaXY2());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_dcaz"), dielectronCandidate.dcaZ1());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_dcaz"), dielectronCandidate.dcaZ2());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_pt"), dielectronCandidate.pt());
        registry.fill(HIST("MCD/Dielectron/h_dielectron_pt_vs_jet_pt"), jet.pt(), dielectronCandidate.pt());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_eta"), dielectronCandidate.eta());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_Lxy"), dielectronCandidate.lxy());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_Lz"), dielectronCandidate.lz());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_tau_xy"), dielectronCandidate.tauxy());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_tau_z"), dielectronCandidate.tauz());
        registry.fill(HIST("MCD/Dielectron/h_MCD_dielectron_origin"), dielectronCandidate.mcDecision());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_pt"), dielectronCandidate.pt1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_pt"), dielectronCandidate.pt2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_eta"), dielectronCandidate.eta1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_eta"), dielectronCandidate.eta2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_phi"), dielectronCandidate.phi1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_phi"), dielectronCandidate.phi2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_nsigma"), dielectronCandidate.tpcNSigmaEl1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_nsigma"), dielectronCandidate.tpcNSigmaEl2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_nsigma_vs_pt"), dielectronCandidate.pt1(), dielectronCandidate.tpcNSigmaEl1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_nsigma_vs_pt"), dielectronCandidate.pt2(), dielectronCandidate.tpcNSigmaEl2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_ncls_found"), dielectronCandidate.tpcNClsFound1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_ncls_found"), dielectronCandidate.tpcNClsFound2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_chi2_per_cluster"), dielectronCandidate.tpcChi2NCl1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_tpc_chi2_per_cluster"), dielectronCandidate.tpcChi2NCl2());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_its_chi2_per_cluster"), dielectronCandidate.itsChi2NCl1());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_electron_its_chi2_per_cluster"), dielectronCandidate.itsChi2NCl2());
        registry.fill(HIST("MCD/Multidimensional/h_MCD_5d_allParticles"), jet.pt(), z, dielectronCandidate.mass(), dielectronCandidate.tauxy(), dielectronCandidate.mcDecision());
        if (dielectronCandidate.mcDecision() == 1) {
          registry.fill(HIST("MCD/Multidimensional/h_MCD_4d_PromptJpsi"), jet.pt(), z, dielectronCandidate.mass(), dielectronCandidate.tauxy());
        } else if (dielectronCandidate.mcDecision() == 2) {
          registry.fill(HIST("MCD/Multidimensional/h_MCD_4d_NonPromptJpsi"), jet.pt(), z, dielectronCandidate.mass(), dielectronCandidate.tauxy());
        }
      } // end of loop over dielectron candidates
      registry.fill(HIST("MCD/Jet/h_electron_pairs_per_jet"), electronPairs);
      if (electronPairs == 0) {
        continue;
      }
      for (const auto& track : jetTracks) {
        sumConstituentPt += track.pt();
        registry.fill(HIST("MCD/SingleTrack/h_MCD_jet_track_pt"), track.pt());
        registry.fill(HIST("MCD/SingleTrack/h_MCD_jet_track_eta"), track.eta());
      }
      registry.fill(HIST("MCD/Jet/h_MCD_jet_pt"), jet.pt());
      registry.fill(HIST("MCD/Jet/h_MCD_jet_eta"), jet.eta());
      registry.fill(HIST("MCD/Jet/h_MCD_jet_phi"), jet.phi());
      const int constituentMultiplicity = static_cast<int>(jetTracks.size()) + 2 * electronPairs; // adding 2 for the dielectron candidate
      registry.fill(HIST("MCD/Jet/h_constituent_multiplicity_vs_jet_pt"), jet.pt(), constituentMultiplicity);
      registry.fill(HIST("MCD/Jet/h_sum_constituent_pt_vs_jet_pt"), jet.pt(), sumConstituentPt);
      registry.fill(HIST("MCD/SingleTrack/h_MCD_n_tracks_per_jet"), constituentMultiplicity);
    } // end of loop over jets
  }

  PROCESS_SWITCH(DqFragmentationFunctionTask, FragmentationFunctionMCD, "Jpsi fragmentation function detector-level MC", false);

  void FragmentationFunctionMCP(
    aod::JetMcCollision const& mccollision,
    soa::Join<aod::DielectronChargedMCParticleLevelJets,
              aod::DielectronChargedMCParticleLevelJetConstituents> const& jets,
    aod::CandidatesDielectronMCP const&,
    aod::JetParticles const&)
  {
    // Event selection
    if (useMCPCuts && (!jetderiveddatautilities::selectCollision(mccollision, eventSelectionBits, skipMBgapEvents.value, applyRCTflag.value, rctLabel.value) || !(std::abs(mccollision.posZ()) < vertexZCut))) {
      return;
    }
    for (const auto& jet : jets) {
      if (useMCPCuts) {
        if (jet.pt() < jetPtMin || jet.pt() > jetPtMax ||
            std::abs(jet.eta()) > jetEtaMax) {
          continue;
        }
      }
      const auto jetParticles = jet.template tracks_as<aod::JetParticles>();
      const auto dielectronCandidates = jet.template candidates_as<aod::CandidatesDielectronMCP>();
      int electronPairs = 0;
      double sumConstituentPt = 0.;
      for (const auto& dielectronCandidate : dielectronCandidates) {
        if (useMCPCuts && !passesPairCuts(dielectronCandidate)) {
          continue;
        }
        const double z = calculateZ(jet, dielectronCandidate);
        if (z < 0.) {
          continue;
        }
        ++electronPairs;

        // Because I dont have access to the individual electron pt, I will use the dielectron candidate pt as a proxy for the sum of the two electrons' pt.
        sumConstituentPt += dielectronCandidate.pt();

        const double axisDistance = jetutilities::deltaR(jet, dielectronCandidate);

        // Calculate Lxy for the J/psi candidate
        const double dx = dielectronCandidate.vx() - mccollision.posX();
        const double dy = dielectronCandidate.vy() - mccollision.posY();

        const double pt = std::sqrt(
          dielectronCandidate.px() * dielectronCandidate.px() +
          dielectronCandidate.py() * dielectronCandidate.py());

        double lxyMC = 0.;
        double tauxyProjected = 0.;

        if (pt > 0.) {
          lxyMC = (dx * dielectronCandidate.px() +
                   dy * dielectronCandidate.py()) /
                  pt;

          tauxyProjected = lxyMC * dielectronCandidate.m() / pt;
        }
        registry.fill(HIST("MCP/Dielectron/h_MCP_jpsi_pt"), dielectronCandidate.pt());
        registry.fill(HIST("MCP/Dielectron/h_dielectron_pt_vs_jet_pt"), jet.pt(), dielectronCandidate.pt());
        registry.fill(HIST("MCP/Dielectron/h_MCP_jpsi_lxy"), lxyMC);
        registry.fill(HIST("MCP/Dielectron/h_MCP_jpsi_tauxy"), tauxyProjected);
        registry.fill(HIST("MCP/Jet/h_MCP_z"), z);
        registry.fill(HIST("MCP/Jet/h_MCP_jet_distance"), axisDistance);
        registry.fill(HIST("MCP/Jet/h_MCP_jet_distance_vs_projection"), axisDistance, z);
        registry.fill(HIST("MCP/Dielectron/h_MCP_dielectron_mass"), dielectronCandidate.m());
        registry.fill(HIST("MCP/Dielectron/h_MCP_dielectron_eta"), dielectronCandidate.eta());
        registry.fill(HIST("MCP/Dielectron/h_MCP_dielectron_origin"), dielectronCandidate.origin());
        registry.fill(HIST("MCP/Multidimensional/h_MCP_4d_allParticles"), jet.pt(), z, tauxyProjected, dielectronCandidate.origin());
        if (dielectronCandidate.origin() == 1) {
          registry.fill(HIST("MCP/Multidimensional/h_MCP_4d_PromptJpsi"), jet.pt(), z, dielectronCandidate.m(), tauxyProjected);
        } else if (dielectronCandidate.origin() == 2) {
          registry.fill(HIST("MCP/Multidimensional/h_MCP_4d_NonPromptJpsi"), jet.pt(), z, dielectronCandidate.m(), tauxyProjected);
        }
      } // end of loop over dielectron candidates
      registry.fill(HIST("MCP/Jet/h_electron_pairs_per_jet"), electronPairs);
      if (electronPairs == 0) {
        continue;
      }
      for (const auto& particle : jetParticles) {
        sumConstituentPt += particle.pt();
        registry.fill(HIST("MCP/SingleTrack/h_MCP_jet_particle_pt_ElectronPairInfoMissing"), particle.pt()); // named ElectronPairInfoMissing, because for MCP, I don't have access to the individual electron pt.
        registry.fill(HIST("MCP/SingleTrack/h_MCP_jet_particle_eta_ElectronPairInfoMissing"), particle.eta());
      }
      const int constituentMultiplicity = static_cast<int>(jetParticles.size()) + 2 * electronPairs;
      registry.fill(HIST("MCP/Jet/h_constituent_multiplicity_vs_jet_pt"), jet.pt(), constituentMultiplicity);
      registry.fill(HIST("MCP/Jet/h_sum_constituent_pt_vs_jet_pt"), jet.pt(), sumConstituentPt);
      registry.fill(HIST("MCP/Jet/h_MCP_jet_pt"), jet.pt());
      registry.fill(HIST("MCP/Jet/h_MCP_jet_eta"), jet.eta());
      registry.fill(HIST("MCP/Jet/h_MCP_jet_phi"), jet.phi());
      registry.fill(HIST("MCP/SingleTrack/h_MCP_n_particles_per_jet"), constituentMultiplicity);
    } // end of loop over jets
  }

  PROCESS_SWITCH(DqFragmentationFunctionTask, FragmentationFunctionMCP, "Jpsi fragmentation function particle-level MC", false);

  void processMcMatching(
    aod::JetMcCollisions const& mccollisions,
    aod::JetCollisionsMCD const& collisions,
    DielectronMCDJets const& mcdjets,
    DielectronMCPJets const&,
    aod::CandidatesDielectronMCD const&,
    aod::CandidatesDielectronMCP const&,
    aod::JetTracksMCD const& tracksMCD,
    aod::JetParticles const& particlesMCP)
  {
    for (const auto& mccollision : mccollisions) {
      registry.fill(HIST("Matching/h_collision_counter"), 0.5);

      if (useMCPCuts && (!jetderiveddatautilities::selectCollision(mccollision, eventSelectionBits, skipMBgapEvents.value, applyRCTflag.value, rctLabel.value) || !(std::abs(mccollision.posZ()) < vertexZCut))) {
        continue;
      }
      registry.fill(HIST("Matching/h_collision_counter"), 1.5);

      const auto collisionsPerMCCollision =
        collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {
        registry.fill(HIST("Matching/h_collision_counter"), 2.5);
        // Event selection
        if (useCuts && (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBgapEvents.value, applyRCTflag.value, rctLabel.value) || !(std::abs(collision.posZ()) < vertexZCut))) {
          continue;
        }
        registry.fill(HIST("Matching/h_collision_counter"), 3.5);
        const auto mcdJetsPerCollision =
          mcdjets.sliceBy(mcdJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : mcdJetsPerCollision) {
          registry.fill(HIST("Matching/h_jet_counter"), 0.5);
          if (useCuts) {
            if (mcdjet.pt() < jetPtMin || mcdjet.pt() > jetPtMax ||
                std::abs(mcdjet.eta()) > jetEtaMax) {
              continue;
            }
          }
          registry.fill(HIST("Matching/h_jet_counter"), 1.5);

          if (!mcdjet.has_matchedJetCand()) {
            continue;
          }
          registry.fill(HIST("Matching/h_jet_counter"), 2.5);

          const auto matchedMCPJets =
            mcdjet.template matchedJetCand_as<DielectronMCPJets>();
          for (const auto& mcpjet : matchedMCPJets) {
            if (useMCPCuts) {
              if (mcpjet.pt() < jetPtMin || mcpjet.pt() > jetPtMax ||
                  std::abs(mcpjet.eta()) > jetEtaMax) {
                continue;
              }
            }
            registry.fill(HIST("Matching/h_jet_counter"), 3.5);
            for (const auto& mcdCandidate :
                 mcdjet.template candidates_as<aod::CandidatesDielectronMCD>()) {
              registry.fill(HIST("Matching/h_candidate_counter"), 0.5);

              const double zMCD = calculateZ(mcdjet, mcdCandidate);
              if (zMCD < 0.) {
                continue;
              }
              if (useCuts &&
                  (!passesPairCuts(mcdCandidate) ||
                   !passesElectronCuts(mcdCandidate, 1) ||
                   !passesElectronCuts(mcdCandidate, 2))) {
                continue;
              }
              registry.fill(HIST("Matching/h_candidate_counter"), 1.5);

              const auto matchedParticleId =
                jetcandidateutilities::matchedParticleId(
                  mcdCandidate, tracksMCD, particlesMCP);

              for (const auto& mcpCandidate :
                   mcpjet.template candidates_as<aod::CandidatesDielectronMCP>()) {
                if (mcpCandidate.mcParticleId() != matchedParticleId) {
                  continue;
                }

                const double zMCP = calculateZ(mcpjet, mcpCandidate);
                if (zMCP < 0.) {
                  continue;
                }
                registry.fill(HIST("Matching/h_candidate_counter"), 2.5);
                if (useMCPCuts && !passesPairCuts(mcpCandidate)) {
                  continue;
                }
                registry.fill(HIST("Matching/h_candidate_counter"), 3.5);
                const double mcpCandidatePt = std::sqrt(mcpCandidate.px() * mcpCandidate.px() + mcpCandidate.py() * mcpCandidate.py());
                double tauxyMCP = 0.;
                if (mcpCandidatePt > 0.) {
                  const double dxMCP = mcpCandidate.vx() - mccollision.posX();
                  const double dyMCP = mcpCandidate.vy() - mccollision.posY();
                  const double lxyMCP = (dxMCP * mcpCandidate.px() + dyMCP * mcpCandidate.py()) / mcpCandidatePt;
                  tauxyMCP = lxyMCP * mcpCandidate.m() / mcpCandidatePt;
                }

                if (mcpCandidate.origin() == 1) {
                  registry.fill(HIST("Matching/h_response_prompt"), mcdjet.pt(), zMCD, mcpjet.pt(), zMCP);
                  registry.fill(HIST("Matching/h_2D_prompt_z"), zMCD, zMCP);
                  registry.fill(HIST("Matching/h_2D_prompt_jetpt"), mcdjet.pt(), mcpjet.pt());
                  registry.fill(HIST("Matching/h_2D_prompt_jpsipt"), mcdCandidate.pt(), mcpCandidate.pt());
                  registry.fill(HIST("Matching/h_2D_prompt_tauxy"), mcdCandidate.tauxy(), tauxyMCP);

                } else if (mcpCandidate.origin() == 2) {
                  registry.fill(HIST("Matching/h_response_nonprompt"), mcdjet.pt(), zMCD, mcpjet.pt(), zMCP);
                  registry.fill(HIST("Matching/h_2D_nonprompt_z"), zMCD, zMCP);
                  registry.fill(HIST("Matching/h_2D_nonprompt_jetpt"), mcdjet.pt(), mcpjet.pt());
                  registry.fill(HIST("Matching/h_2D_nonprompt_jpsipt"), mcdCandidate.pt(), mcpCandidate.pt());
                  registry.fill(HIST("Matching/h_2D_nonprompt_tauxy"), mcdCandidate.tauxy(), tauxyMCP);
                }
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DqFragmentationFunctionTask, processMcMatching, "MC collision, jet and candidate matching", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DqFragmentationFunctionTask>(cfgc, TaskName{"dq-fragmentation-function"})};
}
