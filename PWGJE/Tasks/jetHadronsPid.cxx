// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetHadronsPid.cxx
/// \brief Analysis of hadrons in jets
/// \author Leonard Lorenc, WUT Warsaw, leonard.lorenc@cern.ch
/// \author Aleksandra Mulewicz, WUT Warsaw, aleksandra.mulewicz@cern.ch
/// \author Hubert Zalewski, WUT Warsaw, hubert.zalewski@cern.ch
/// \author Janik Małgorzata malgorzata.janik@cern.ch
/// \author Daniela Ruggiano daniela.ruggiano@cern.ch

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TVector3.h>

#include <Rtypes.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;

using StandardEvents = soa::Join<aod::Collisions, aod::EvSels>;

using JetEvents = soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionPIs>;
using MCPJetEvents = soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>;

using ChargedMCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
using ChargedMCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

using HadronTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,
                               aod::pidTPCFullPi, aod::pidTOFFullPi,
                               aod::pidTPCFullKa, aod::pidTOFFullKa,
                               aod::pidTPCFullPr, aod::pidTOFFullPr,
                               aod::pidTOFbeta>;
using HadronTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA,
                                 aod::pidTPCFullPi, aod::pidTOFFullPi,
                                 aod::pidTPCFullKa, aod::pidTOFFullKa,
                                 aod::pidTPCFullPr, aod::pidTOFFullPr,
                                 aod::pidTOFbeta, aod::McTrackLabels>;

struct JetHadronsPid {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<std::string> eventSelections{"eventSelections", "sel8+NoITSROFrameBorder+NoTimeFrameBorder+NoSameBunchPileup+IsGoodZvtxFT0vsPV", "choose event selection"};
  std::vector<int> eventSelectionBits;

  Configurable<bool> isppRefAnalysis{"isppRefAnalysis", true, "Is ppRef analysis"};
  Configurable<double> cfgEtaJetMax{"cfgEtaJetMax", 0.5, "max jet eta"};

  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "Reject events near the ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "Reject events near the TF border"};
  Configurable<bool> requireVtxITSTPC{"requireVtxITSTPC", true, "Require at least one ITS-TPC matched track"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "Reject events with same-bunch pileup collisions"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "Require consistent FT0 vs PV z-vertex"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "Require at least one vertex track matched to TOF"};

  Configurable<double> minJetPt{"minJetPt", 10, "Minimum pt of the jet after bkg subtraction"};
  Configurable<double> maxJetPt{"maxJetPt", 1e+06, "Maximum pt of the jet after bkg subtraction"};
  Configurable<double> rJet{"rJet", 0.4, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<bool> applyAreaCut{"applyAreaCut", false, "apply area cut"};
  Configurable<double> minNormalizedJetArea{"minNormalizedJetArea", 0.6, "Minimum normalized area cut to reject fake jets"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};

  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 70, "minimum number of TPC crossed pad rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.05, "minimum pt of the tracks"};
  Configurable<double> maxPt{"maxPt", 4.0, "maximum pt of the tracks for PID analysis"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.2, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.1, "Maximum DCAz"};
  Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", true, "set MC default parameters"};

  struct : ConfigurableGroup {
    Configurable<int> pidMethod{"pidMethod", 2, "Variable for choosing pid method"};
    Configurable<float> rejectionSigma{"rejectionSigma", 3.0, "Rejection sigma for tof and tpc pid"};
    Configurable<float> ptThresholdPion{"ptThresholdPion", 0.75, "Threshold for pt for different pid methods"};
    Configurable<float> ptThresholdKaon{"ptThresholdKaon", 0.5, "Threshold for pt for different pid methods"};
    Configurable<float> ptThresholdProton{"ptThresholdProton", 0.75, "Threshold for pt for different pid methods"};
    Configurable<float> nSigmaCutPi{"nSigmaCutPi", 3.0, "Acceptence cut for pid methods"};
    Configurable<float> nSigmaCutKa{"nSigmaCutKa", 3.0, "Acceptence cut for pid methods"};
    Configurable<float> nSigmaCutPr{"nSigmaCutPr", 3.0, "Acceptence cut for pid methods"};
    Configurable<float> minPtPion{"minPtPion", 0.2, "Minimal Pion pt"};
    Configurable<float> maxPtPion{"maxPtPion", 4.0, "Maximum Pion pt"};
    Configurable<float> minPtKaon{"minPtKaon", 0.3, "Minimal Kaon pt"};
    Configurable<float> maxPtKaon{"maxPtKaon", 3.0, "Maximum Kaon pt"};
    Configurable<float> minPtProton{"minPtProton", 0.5, "Minimal Proton pt"};
    Configurable<float> maxPtProton{"maxPtProton", 4.0, "Maximum Proton pt"};
  } cfg;

  struct : ConfigurableGroup {
    // ----------------- Kinematics -----------------
    ConfigurableAxis axisPt{"axisPt", {120, 0.0, 4.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisEta{"axisEta", {100, -1.0, 1.0}, "#eta"};
    ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, TwoPI}, "#varphi"};
    ConfigurableAxis axisMomentum{"axisMomentum", {300, 0.05, 30.0}, "#it{p} (GeV/#it{c})"};

    // ----------------- Spatial / DCA -----------------
    ConfigurableAxis axisDcaXY{"axisDcaXY", {300, -0.5, 0.5}, "DCA_{xy} (cm)"};
    ConfigurableAxis axisDcaZ{"axisDcaZ", {300, -0.5, 0.5}, "DCA_{z} (cm)"};
    ConfigurableAxis axisZVtx{"axisZVtx", {200, -20.0, 20.0}, "z_{vtx} (cm)"};
    ConfigurableAxis axisVtxXY{"axisVtxXY", {200, -1.0, 1.0}, "vtx_{xy} (cm)"};

    // ----------------- PID -----------------
    ConfigurableAxis axisTPCSignal{"axisTPCSignal", {400, 20.0, 2000.0}, "TPC d#it{E}/d#it{x} (a.u.)"};
    ConfigurableAxis axisTOFBeta{"axisTOFBeta", {240, 0.0, 1.2}, "#beta_{TOF}"};
    ConfigurableAxis axisNSigmaTPC{"axisNSigmaTPC", {200, -5.0, 5.0}, "n#sigma_{TPC}"};
    ConfigurableAxis axisNSigmaTOF{"axisNSigmaTOF", {200, -5.0, 5.0}, "n#sigma_{TOF}"};
    ConfigurableAxis axisNSigmaCOMB{"axisNSigmaCOMB", {200, 0.0, 15.0}, "n#sigma_{COMB}"};
    ConfigurableAxis axisParticleSpecies{"axisParticleSpecies", {3, 0.5, 3.5}, "Particle species"};

    // ----------------- Cone / Jet Dynamics -----------------
    ConfigurableAxis axisDeltaEtaCone{"axisDeltaEtaCone", {80, -0.8, 0.8}, "#Delta#eta"};
    ConfigurableAxis axisDeltaPhiCone{"axisDeltaPhiCone", {80, -0.8, 0.8}, "#Delta#varphi"};
    ConfigurableAxis axisPtTrackCone{"axisPtTrackCone", {100, 0.0, 10.0}, "#it{p}_{T}^{track} (GeV/#it{c})"};
    ConfigurableAxis axisJetPt{"axisJetPt", {200, 0.0, 100.0}, "#it{p}_{T}^{jet} (GeV/#it{c})"};
    ConfigurableAxis axisJetSubPt{"axisJetSubPt", {200, 0.0, 10.0}, "#it{p}_{T}^{sub} (GeV/#it{c})"};
    ConfigurableAxis axisJetArea{"axisJetArea", {100, 0.0, 1.5}, "A_{jet}"};
    ConfigurableAxis axisJetNConst{"axisJetNConst", {101, -0.5, 100.5}, "N_{constituents}"};

    // ----------------- Multiplicities / Counting -----------------
    ConfigurableAxis axisEventMultiplicity{"axisEventMultiplicity", {201, -0.5, 200.5}, "N_{particles/contributors}"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {31, -0.5, 30.5}, "N per jet"};
    ConfigurableAxis axisNPions{"axisNPions", {21, -0.5, 20.5}, "N_{#pi} per jet"};
    ConfigurableAxis axisNKaons{"axisNKaons", {11, -0.5, 10.5}, "N_{K} per jet"};
    ConfigurableAxis axisNProtons{"axisNProtons", {11, -0.5, 10.5}, "N_{p} per jet"};

    // ----------------- Others -----------------
    ConfigurableAxis axisMCPdg{"axisMCPdg", {6001, -3000.5, 3000.5}, "PDG Code"};
    ConfigurableAxis axisTrackCategory{"axisTrackCategory", {2, 0.5, 2.5}, "Track category"};
    ConfigurableAxis axisEventCounter{"axisEventCounter", {1, 0.5, 1.5}, "N_{events}"};

  } axis;

  o2::aod::ITSResponse itsResponse;

  Preslice<HadronTracks> tracksPerCollision = aod::track::collisionId;
  Preslice<aod::McParticles> mcParticlesPerCollision = aod::mcparticle::mcCollisionId;

  void init(InitContext const&)
  {
    if (setMCDefaultItsParams) {
      itsResponse.setMCDefaultParameters();
    }

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    // --- DATA: PURE TRACKS ---
    registryData.add("data/pure/general/tpc_signal_vs_p", "TPC n#sigma_{#pi} vs transverse momentum ;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/pure/general/tof_beta_vs_p", "TOF #beta vs momentum;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});

    // Pions
    registryData.add("data/pure/pions/pion_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected pion candidates;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/pure/pions/pion_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected pion candidates;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/pure/pions/pion_pure_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/pions/pion_pure_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/pions/pion_pure_combined_tof_tpc", "Combined TOF+TPC n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{COMB}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaCOMB});

    registryData.add("data/pure/pions/pion_pure_pt", "Transverse-momentum distribution of PID-selected pion candidates;#it{p}_{T} (GeV/#it{c});N_{#pi}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/pions/pion_pure_eta", "Pseudorapidity distribution of PID-selected pion candidates;#eta;N_{#pi}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/pions/pion_pure_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected pion candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/pions/pion_pure_dcaz", "DCA_{z} vs transverse momentum for PID-selected pion candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/pions/pion_pure_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected pion candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    registryData.add("data/pure/pions/pos/pion_pure_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/pions/pos/pion_pure_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/pions/pos/pion_pure_pos_pt", "Transverse-momentum distribution of PID-selected #pi^{+} candidates;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/pions/pos/pion_pure_pos_eta", "Pseudorapidity distribution of PID-selected #pi^{+} candidates;#eta;N_{#pi^{+}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/pions/pos/pion_pure_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #pi^{+} candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/pions/pos/pion_pure_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected #pi^{+} candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/pions/pos/pion_pure_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected #pi^{+} candidates;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryData.add("data/pure/pions/pos/pion_pure_pos_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected #pi^{+} candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    registryData.add("data/pure/pions/neg/pion_pure_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/pions/neg/pion_pure_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/pions/neg/pion_pure_neg_pt", "Transverse-momentum distribution of PID-selected #pi^{-} candidates;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/pions/neg/pion_pure_neg_eta", "Pseudorapidity distribution of PID-selected #pi^{-} candidates;#eta;N_{#pi^{-}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/pions/neg/pion_pure_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #pi^{-} candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/pions/neg/pion_pure_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected #pi^{-} candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/pions/neg/pion_pure_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected #pi^{-} candidates;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryData.add("data/pure/pions/neg/pion_pure_neg_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected #pi^{-} candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    // Kaons
    registryData.add("data/pure/kaons/kaon_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected kaon candidates;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/pure/kaons/kaon_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected kaon candidates;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/pure/kaons/kaon_pure_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected kaon candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/kaons/kaon_pure_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected kaon candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/kaons/kaon_pure_combined_tof_tpc", "Combined TOF+TPC n#sigma_{K} vs transverse momentum for PID-selected kaon candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{COMB}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaCOMB});

    registryData.add("data/pure/kaons/kaon_pure_pt", "Transverse-momentum distribution of PID-selected kaon candidates;#it{p}_{T} (GeV/#it{c});N_{K}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/kaons/kaon_pure_eta", "Pseudorapidity distribution of PID-selected kaon candidates;#eta;N_{K}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/kaons/kaon_pure_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected kaon candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/kaons/kaon_pure_dcaz", "DCA_{z} vs transverse momentum for PID-selected kaon candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/kaons/kaon_pure_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected kaon candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    registryData.add("data/pure/kaons/pos/kaon_pure_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_pt", "Transverse-momentum distribution of PID-selected K^{+} candidates;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_eta", "Pseudorapidity distribution of PID-selected K^{+} candidates;#eta;N_{K^{+}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected K^{+} candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected K^{+} candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected K^{+} candidates;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryData.add("data/pure/kaons/pos/kaon_pure_pos_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected K^{+} candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    registryData.add("data/pure/kaons/neg/kaon_pure_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_pt", "Transverse-momentum distribution of PID-selected K^{-} candidates;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_eta", "Pseudorapidity distribution of PID-selected K^{-} candidates;#eta;N_{K^{-}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected K^{-} candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected K^{-} candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected K^{-} candidates;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryData.add("data/pure/kaons/neg/kaon_pure_neg_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected K^{-} candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    // Protons
    registryData.add("data/pure/protons/proton_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected proton candidates;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/pure/protons/proton_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected proton candidates;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/pure/protons/proton_pure_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/protons/proton_pure_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/protons/proton_pure_combined_tof_tpc", "Combined TOF+TPC n#sigma_{p} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{COMB}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaCOMB});

    registryData.add("data/pure/protons/proton_pure_pt", "Transverse-momentum distribution of PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/protons/proton_pure_eta", "Pseudorapidity distribution of PID-selected proton candidates;#eta;N_{p}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/protons/proton_pure_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/protons/proton_pure_dcaz", "DCA_{z} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/protons/proton_pure_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected proton candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    registryData.add("data/pure/protons/pos/proton_pure_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/protons/pos/proton_pure_pos_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/protons/pos/proton_pure_pos_pt", "Transverse-momentum distribution of PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/protons/pos/proton_pure_pos_eta", "Pseudorapidity distribution of PID-selected proton candidates;#eta;N_{p}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/protons/pos/proton_pure_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/protons/pos/proton_pure_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/protons/pos/proton_pure_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected proton candidates;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryData.add("data/pure/protons/pos/proton_pure_pos_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected proton candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    registryData.add("data/pure/protons/neg/proton_pure_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected antiproton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/pure/protons/neg/proton_pure_neg_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected antiproton candidates;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/pure/protons/neg/proton_pure_neg_pt", "Transverse-momentum distribution of PID-selected antiproton candidates;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/pure/protons/neg/proton_pure_neg_eta", "Pseudorapidity distribution of PID-selected antiproton candidates;#eta;N_{#bar{p}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/pure/protons/neg/proton_pure_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected antiproton candidates;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/pure/protons/neg/proton_pure_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected antiproton candidates;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/pure/protons/neg/proton_pure_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected antiproton candidates;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryData.add("data/pure/protons/neg/proton_pure_neg_eta_phi", "Azimuthal angle vs pseudorapidity for PID-selected antiproton candidates;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});

    // --- DATA: COLLISIONS ---
    registryData.add("data/pure/collisions/z_vertex", "Collision vertex;z_{vtx} (cm);N_{collisions}", HistType::kTH1F, {axis.axisZVtx});
    registryData.add("data/pure/collisions/x_vertex", "Collision vertex;x_{vtx} (cm);N_{collisions}", HistType::kTH1F, {axis.axisVtxXY});
    registryData.add("data/pure/collisions/y_vertex", "Collision vertex;y_{vtx} (cm);N_{collisions}", HistType::kTH1F, {axis.axisVtxXY});
    registryData.add("data/pure/collisions/xy_vertex", "Collision vertex;x_{vtx} (cm);y_{vtx} (cm)", HistType::kTH2F, {axis.axisVtxXY, axis.axisVtxXY});
    registryData.add("data/pure/collisions/n_contributors", "PV contributors;N_{contributors};N_{collisions}", HistType::kTH1I, {axis.axisEventMultiplicity});

    // --- DATA: JETS GENERAL ---
    registryData.add("data/jets/general/n_events", "Event counter", HistType::kTH1F, {axis.axisEventCounter});
    registryData.add("data/jets/general/n_events_raw", "All events", HistType::kTH1F, {axis.axisEventCounter});
    registryData.add("data/jets/general/collision_multiplicity", "Number of Jets in Collison;Number of Jets;Number of Collisions", HistType::kTH1F, {axis.axisMultiplicity});
    registryData.add("data/jets/general/z_vtx", "Z-Vertex Distribution", HistType::kTH1F, {axis.axisZVtx});

    registryData.add("data/jets/general/jet_pt", "Jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetPt});
    registryData.add("data/jets/general/jet_pt_subtracted", "Jet #it{p}_{T} subtracted;#it{p}_{T}^{sub} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetSubPt});
    registryData.add("data/jets/general/jet_pt_raw_vs_sub", "Raw vs sub jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T}^{sub} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisJetSubPt});
    registryData.add("data/jets/general/jet_eta", "Jet #eta;#eta_{jet};N_{jets}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/general/jet_phi", "Jet #varphi;#varphi_{jet};N_{jets}", HistType::kTH1F, {axis.axisPhi});
    registryData.add("data/jets/general/jet_area", "Jet area;A_{jet};N_{jets}", HistType::kTH1F, {axis.axisJetArea});
    registryData.add("data/jets/general/jet_n_constituents", "Jet multiplicity;N_{constituents};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});
    registryData.add("data/jets/general/jet_n_selected_constituents", "Jet multiplicity of selected tracks;N_{constituents}^{selected};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});

    registryData.add("data/jets/general/tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for tracks in selected jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/jets/general/tof_beta_vs_p", "TOF #beta vs momentum for tracks in selected jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/jets/general/total_tracks_pure_vs_jets", "Selected tracks and tracks used in jets;Track category;N_{tracks}", HistType::kTH1F, {axis.axisTrackCategory});

    auto hTrackTotals = registryData.get<TH1>(HIST("data/jets/general/total_tracks_pure_vs_jets"));
    hTrackTotals->GetXaxis()->SetBinLabel(1, "Pure tracks");
    hTrackTotals->GetXaxis()->SetBinLabel(2, "Tracks in jets");

    // --- DATA: JET CONE ---
    registryData.add("data/jets/cone/all_tracks_2d", "Track distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});

    // Pions
    registryData.add("data/jets/cone/pions/pion_pos_2d", "PID-selected #pi^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryData.add("data/jets/cone/pions/pion_neg_2d", "PID-selected #pi^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});

    // Kaons
    registryData.add("data/jets/cone/kaons/kaon_pos_2d", "PID-selected K^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryData.add("data/jets/cone/kaons/kaon_neg_2d", "PID-selected K^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});

    // Protons
    registryData.add("data/jets/cone/protons/proton_pos_2d", "PID-selected p distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryData.add("data/jets/cone/protons/proton_neg_2d", "PID-selected #bar{p} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet}", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});

    // 3D Distributions
    registryData.add("data/jets/cone/all_tracks", "Track distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryData.add("data/jets/cone/pions/pion_pos", "PID-selected #pi^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryData.add("data/jets/cone/pions/pion_neg", "PID-selected #pi^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryData.add("data/jets/cone/kaons/kaon_pos", "PID-selected K^{+} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryData.add("data/jets/cone/kaons/kaon_neg", "PID-selected K^{-} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryData.add("data/jets/cone/protons/proton_pos", "PID-selected p distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryData.add("data/jets/cone/protons/proton_neg", "PID-selected #bar{p} distribution inside selected jet cones;#Delta#eta_{track,jet};#Delta#varphi_{track,jet};#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});

    // --- DATA: JET COMPOSITION ---

    // Pions
    registryData.add("data/jets/pions/n_pions_per_jet", "Multiplicity of PID-selected #pi^{#pm} in selected jets;N_{#pi};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryData.add("data/jets/pions/n_pions_vs_jet_pt", "Multiplicity of PID-selected #pi^{#pm} vs jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});N_{#pi}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryData.add("data/jets/pions/pion_pt_vs_jet_pt", "PID-selected #pi^{#pm} #it{p}_{T} vs jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryData.add("data/jets/pions/pos/n_pions_per_jet", "Multiplicity of PID-selected #pi^{+} in selected jets;N_{#pi^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryData.add("data/jets/pions/neg/n_pions_per_jet", "Multiplicity of PID-selected #pi^{-} in selected jets;N_{#pi^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});

    registryData.add("data/jets/pions/pion_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected pion candidates in jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/jets/pions/pion_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected pion candidates in jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/jets/pions/pion_jet_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/pions/pion_jet_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});

    registryData.add("data/jets/pions/pos/pion_jet_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/pions/pos/pion_jet_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/jets/pions/pos/pion_jet_pos_pt", "Transverse-momentum distribution of PID-selected #pi^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/jets/pions/pos/pion_jet_pos_eta", "Pseudorapidity distribution of PID-selected #pi^{+} candidates in jets;#eta;N_{#pi^{+}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/pions/pos/pion_jet_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #pi^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/jets/pions/pos/pion_jet_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected #pi^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/jets/pions/pos/pion_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected #pi^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryData.add("data/jets/pions/neg/pion_jet_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/pions/neg/pion_jet_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/jets/pions/neg/pion_jet_neg_pt", "Transverse-momentum distribution of PID-selected #pi^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/jets/pions/neg/pion_jet_neg_eta", "Pseudorapidity distribution of PID-selected #pi^{-} candidates in jets;#eta;N_{#pi^{-}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/pions/neg/pion_jet_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #pi^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/jets/pions/neg/pion_jet_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected #pi^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/jets/pions/neg/pion_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected #pi^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // Kaons
    registryData.add("data/jets/kaons/n_kaons_per_jet", "Multiplicity of PID-selected K^{#pm} in selected jets;N_{K};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryData.add("data/jets/kaons/n_kaons_vs_jet_pt", "Multiplicity of PID-selected K^{#pm} vs jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});N_{K}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryData.add("data/jets/kaons/kaon_pt_vs_jet_pt", "PID-selected K^{#pm} #it{p}_{T} vs jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryData.add("data/jets/kaons/pos/n_kaons_per_jet", "Multiplicity of PID-selected K^{+} in selected jets;N_{K^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryData.add("data/jets/kaons/neg/n_kaons_per_jet", "Multiplicity of PID-selected K^{-} in selected jets;N_{K^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});

    registryData.add("data/jets/kaons/kaon_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected kaon candidates in jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/jets/kaons/kaon_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected kaon candidates in jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/jets/kaons/kaon_jet_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/kaons/kaon_jet_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});

    registryData.add("data/jets/kaons/pos/kaon_jet_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_pt", "Transverse-momentum distribution of PID-selected K^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_eta", "Pseudorapidity distribution of PID-selected K^{+} candidates in jets;#eta;N_{K^{+}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected K^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected K^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/jets/kaons/pos/kaon_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected K^{+} candidates in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryData.add("data/jets/kaons/neg/kaon_jet_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_pt", "Transverse-momentum distribution of PID-selected K^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_eta", "Pseudorapidity distribution of PID-selected K^{-} candidates in jets;#eta;N_{K^{-}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected K^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected K^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/jets/kaons/neg/kaon_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected K^{-} candidates in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // Protons
    registryData.add("data/jets/protons/n_protons_per_jet", "Multiplicity of PID-selected p+#bar{p} in selected jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryData.add("data/jets/protons/n_protons_vs_jet_pt", "Multiplicity of PID-selected p+#bar{p} vs jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryData.add("data/jets/protons/proton_pt_vs_jet_pt", "PID-selected p+#bar{p} #it{p}_{T} vs jet #it{p}_{T};#it{p}_{T}^{jet} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryData.add("data/jets/protons/pos/n_protons_per_jet", "Multiplicity of PID-selected p in selected jets;N_{p};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryData.add("data/jets/protons/neg/n_protons_per_jet", "Multiplicity of PID-selected #bar{p} in selected jets;N_{#bar{p}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});

    registryData.add("data/jets/protons/proton_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected proton candidates in jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryData.add("data/jets/protons/proton_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected proton candidates in jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryData.add("data/jets/protons/proton_jet_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected p+#bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/protons/proton_jet_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected p+#bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});

    registryData.add("data/jets/protons/pos/proton_jet_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected p candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/protons/pos/proton_jet_pos_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected p candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/jets/protons/pos/proton_jet_pos_pt", "Transverse-momentum distribution of PID-selected p candidates in jets;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/jets/protons/pos/proton_jet_pos_eta", "Pseudorapidity distribution of PID-selected p candidates in jets;#eta;N_{p}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/protons/pos/proton_jet_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected p candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/jets/protons/pos/proton_jet_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected p candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/jets/protons/pos/proton_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected p candidates in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryData.add("data/jets/protons/neg/proton_jet_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected #bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/jets/protons/neg/proton_jet_neg_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected #bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/jets/protons/neg/proton_jet_neg_pt", "Transverse-momentum distribution of PID-selected #bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/jets/protons/neg/proton_jet_neg_eta", "Pseudorapidity distribution of PID-selected #bar{p} candidates in jets;#eta;N_{#bar{p}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/jets/protons/neg/proton_jet_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/jets/protons/neg/proton_jet_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected #bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryData.add("data/jets/protons/neg/proton_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for PID-selected #bar{p} candidates in jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryData.add("data/jets/composition/particle_multiplicity_per_jet", "Identified-particle multiplicity in selected jets;Particle species;N_{particles} per jet", HistType::kTH2I, {axis.axisParticleSpecies, axis.axisMultiplicity});

    auto hComposition = registryData.get<TH2>(HIST("data/jets/composition/particle_multiplicity_per_jet"));
    hComposition->GetXaxis()->SetBinLabel(1, "#pi");
    hComposition->GetXaxis()->SetBinLabel(2, "K");
    hComposition->GetXaxis()->SetBinLabel(3, "p");
    hComposition->SetOption("COLZ");
    hComposition->SetStats(false);

    // --- DATA: UNDERLYING EVENT (UE) ---

    // Pions
    registryData.add("data/ue/pions/pion_ue_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/pions/pion_ue_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected pion candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/pions/pion_ue_pt", "Transverse-momentum distribution of PID-selected pion candidates in UE;#it{p}_{T} (GeV/#it{c});N_{#pi}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/pions/pion_ue_eta", "Pseudorapidity distribution of PID-selected pion candidates in UE;#eta;N_{#pi}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/pions/pion_ue_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected pion candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/pions/pion_ue_dcaz", "DCA_{z} vs transverse momentum for PID-selected pion candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/pions/pos/pion_ue_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/pions/pos/pion_ue_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/pions/pos/pion_ue_pos_pt", "Transverse-momentum distribution of PID-selected #pi^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/pions/pos/pion_ue_pos_eta", "Pseudorapidity distribution of PID-selected #pi^{+} candidates in UE;#eta;N_{#pi^{+}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/pions/pos/pion_ue_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #pi^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/pions/pos/pion_ue_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected #pi^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/pions/neg/pion_ue_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/pions/neg/pion_ue_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/pions/neg/pion_ue_neg_pt", "Transverse-momentum distribution of PID-selected #pi^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/pions/neg/pion_ue_neg_eta", "Pseudorapidity distribution of PID-selected #pi^{-} candidates in UE;#eta;N_{#pi^{-}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/pions/neg/pion_ue_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #pi^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/pions/neg/pion_ue_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected #pi^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    // Kaons
    registryData.add("data/ue/kaons/kaon_ue_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected kaon candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/kaons/kaon_ue_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected kaon candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/kaons/kaon_ue_pt", "Transverse-momentum distribution of PID-selected kaon candidates in UE;#it{p}_{T} (GeV/#it{c});N_{K}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/kaons/kaon_ue_eta", "Pseudorapidity distribution of PID-selected kaon candidates in UE;#eta;N_{K}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/kaons/kaon_ue_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected kaon candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/kaons/kaon_ue_dcaz", "DCA_{z} vs transverse momentum for PID-selected kaon candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/kaons/pos/kaon_ue_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_pt", "Transverse-momentum distribution of PID-selected K^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_eta", "Pseudorapidity distribution of PID-selected K^{+} candidates in UE;#eta;N_{K^{+}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected K^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/kaons/pos/kaon_ue_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected K^{+} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/kaons/neg/kaon_ue_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_pt", "Transverse-momentum distribution of PID-selected K^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_eta", "Pseudorapidity distribution of PID-selected K^{-} candidates in UE;#eta;N_{K^{-}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected K^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/kaons/neg/kaon_ue_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected K^{-} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    // Protons
    registryData.add("data/ue/protons/proton_ue_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected proton candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/protons/proton_ue_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected proton candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/protons/proton_ue_pt", "Transverse-momentum distribution of PID-selected proton candidates in UE;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/protons/proton_ue_eta", "Pseudorapidity distribution of PID-selected proton candidates in UE;#eta;N_{p}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/protons/proton_ue_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected proton candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/protons/proton_ue_dcaz", "DCA_{z} vs transverse momentum for PID-selected proton candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/protons/pos/proton_ue_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected p candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/protons/pos/proton_ue_pos_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected p candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/protons/pos/proton_ue_pos_pt", "Transverse-momentum distribution of PID-selected p candidates in UE;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/protons/pos/proton_ue_pos_eta", "Pseudorapidity distribution of PID-selected p candidates in UE;#eta;N_{p}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/protons/pos/proton_ue_pos_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected p candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/protons/pos/proton_ue_pos_dcaz", "DCA_{z} vs transverse momentum for PID-selected p candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/protons/neg/proton_ue_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected #bar{p} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryData.add("data/ue/protons/neg/proton_ue_neg_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected #bar{p} candidates in UE;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryData.add("data/ue/protons/neg/proton_ue_neg_pt", "Transverse-momentum distribution of PID-selected #bar{p} candidates in UE;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryData.add("data/ue/protons/neg/proton_ue_neg_eta", "Pseudorapidity distribution of PID-selected #bar{p} candidates in UE;#eta;N_{#bar{p}}", HistType::kTH1F, {axis.axisEta});
    registryData.add("data/ue/protons/neg/proton_ue_neg_dcaxy", "DCA_{xy} vs transverse momentum for PID-selected #bar{p} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryData.add("data/ue/protons/neg/proton_ue_neg_dcaz", "DCA_{z} vs transverse momentum for PID-selected #bar{p} candidates in UE;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryData.add("data/ue/general/tracks_n_in_ue", "Number of tracks in UE;N_{tracks};Counts", HistType::kTH1I, {{100, 0, 100, "N_{tracks}"}});

    // =========================================================================
    // --- MC REGISTRY ---
    // =========================================================================

    // --- MC: TRUTH ---

    registryMC.add("mc/truth/general/n_events", "Selected particle-level events;Event counter;N_{events}", HistType::kTH1I, {axis.axisEventCounter});
    registryMC.add("mc/truth/general/z_vtx", "Particle-level Z-vertex distribution;z_{vtx} (cm);N_{events}", HistType::kTH1F, {axis.axisZVtx});

    // Pions
    registryMC.add("mc/truth/pions/mc_gen_pion_pt", "Transverse-momentum distribution of generated primary #pi^{#pm};#it{p}_{T} (GeV/#it{c});N_{#pi}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/truth/pions/pos/mc_gen_pion_pos_pt", "Transverse-momentum distribution of generated primary #pi^{+};#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/truth/pions/neg/mc_gen_pion_neg_pt", "Transverse-momentum distribution of generated primary #pi^{-};#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});

    // Kaons
    registryMC.add("mc/truth/kaons/mc_gen_kaon_pt", "Transverse-momentum distribution of generated primary K^{#pm};#it{p}_{T} (GeV/#it{c});N_{K}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/truth/kaons/pos/mc_gen_kaon_pos_pt", "Transverse-momentum distribution of generated primary K^{+};#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/truth/kaons/neg/mc_gen_kaon_neg_pt", "Transverse-momentum distribution of generated primary K^{-};#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});

    // Protons
    registryMC.add("mc/truth/protons/mc_gen_proton_pt", "Transverse-momentum distribution of generated primary p+#bar{p};#it{p}_{T} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/truth/protons/pos/mc_gen_proton_pos_pt", "Transverse-momentum distribution of generated primary p;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/truth/protons/neg/mc_gen_proton_neg_pt", "Transverse-momentum distribution of generated primary #bar{p};#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});

    // Hadrons
    registryMC.add("mc/truth/hadrons/mc_gen_hadron_pt", "Transverse-momentum distribution of generated primary hadrons;#it{p}_{T} (GeV/#it{c});N_{hadrons}", HistType::kTH1F, {axis.axisPt});

    // --- MC: MATCHED RECONSTRUCTION ---
    registryMC.add("mc/matched/general/n_events", "Event counter;Selection;N_{events}", HistType::kTH1F, {axis.axisEventCounter});
    registryMC.add("mc/matched/general/z_vertex", "Collision vertex;z_{vtx} (cm);N_{collisions}", HistType::kTH1F, {axis.axisZVtx});
    registryMC.add("mc/matched/general/x_vertex", "Collision vertex;x_{vtx} (cm);N_{collisions}", HistType::kTH1F, {axis.axisVtxXY});
    registryMC.add("mc/matched/general/y_vertex", "Collision vertex;y_{vtx} (cm);N_{collisions}", HistType::kTH1F, {axis.axisVtxXY});
    registryMC.add("mc/matched/general/xy_vertex", "Collision vertex;x_{vtx} (cm);y_{vtx} (cm)", HistType::kTH2F, {axis.axisVtxXY, axis.axisVtxXY});
    registryMC.add("mc/matched/general/n_contributors", "PV contributors;N_{contributors};N_{collisions}", HistType::kTH1I, {axis.axisEventMultiplicity});

    // General Hadrons
    registryMC.add("mc/matched/reconstruction/hadrons/rec_hadron_all", "All reconstructed hadrons;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/hadrons/rec_hadron_tof_match", "TOF-matched reconstructed hadrons;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/hadrons/mc_rec_hadron_pt", "True primary reconstructed hadrons;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/hadrons/mc_sec_hadron_pt", "Secondary/fake reconstructed hadrons;#it{p}_{T} (GeV/#it{c});N_{tracks}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/general/tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for reconstructed tracks;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/matched/reconstruction/general/tof_beta_vs_p", "TOF #beta vs momentum for reconstructed tracks;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});

    // Pions
    registryMC.add("mc/matched/reconstruction/pions/rec_pion_all", "All PID #pi;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/contamination_matrix_pion", "PID #pi contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/mc_rec_pion_pt", "True primary #pi;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/mc_sec_pion_pt", "Secondary/fake #pi;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/pion_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for reconstructed #pi;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/matched/reconstruction/pions/pion_tof_beta_vs_p", "TOF #beta vs momentum for reconstructed #pi;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryMC.add("mc/matched/reconstruction/pions/rec_pion_eta", "Pseudorapidity distribution of reconstructed #pi;#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/pions/rec_pion_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed #pi;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/pions/rec_pion_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed #pi;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/pions/rec_pion_dcaz", "DCA_{z} vs transverse momentum for reconstructed #pi;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_all", "All PID #pi^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/pos/contamination_matrix_pion_pos", "PID #pi^{+} contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_tpc", "TPC n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{+};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_tof", "TOF n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{+};#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_tof_matched", "TOF-matched reconstructed #pi^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/pos/mc_rec_pion_pos_pt", "True primary #pi^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/pos/mc_sec_pion_pos_pt", "Secondary/fake #pi^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_eta", "Pseudorapidity distribution of reconstructed #pi^{+};#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed #pi^{+};#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed #pi^{+};#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_dcaz", "DCA_{z} vs transverse momentum for reconstructed #pi^{+};#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/matched/reconstruction/pions/pos/rec_pion_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for reconstructed #pi^{+};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_all", "All PID #pi^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/neg/contamination_matrix_pion_neg", "PID #pi^{-} contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_tpc", "TPC n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{-};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_tof", "TOF n#sigma_{#pi} vs transverse momentum for reconstructed #pi^{-};#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_tof_matched", "TOF-matched reconstructed #pi^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/neg/mc_rec_pion_neg_pt", "True primary #pi^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/neg/mc_sec_pion_neg_pt", "Secondary/fake #pi^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_eta", "Pseudorapidity distribution of reconstructed #pi^{-};#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed #pi^{-};#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed #pi^{-};#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_dcaz", "DCA_{z} vs transverse momentum for reconstructed #pi^{-};#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/matched/reconstruction/pions/neg/rec_pion_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for reconstructed #pi^{-};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // Kaons
    registryMC.add("mc/matched/reconstruction/kaons/rec_kaon_all", "All PID K;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/contamination_matrix_kaon", "PID K contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/mc_rec_kaon_pt", "True primary K;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/mc_sec_kaon_pt", "Secondary/fake K;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/kaon_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for reconstructed K;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/matched/reconstruction/kaons/kaon_tof_beta_vs_p", "TOF #beta vs momentum for reconstructed K;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryMC.add("mc/matched/reconstruction/kaons/rec_kaon_eta", "Pseudorapidity distribution of reconstructed K;#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/kaons/rec_kaon_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed K;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/kaons/rec_kaon_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed K;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/kaons/rec_kaon_dcaz", "DCA_{z} vs transverse momentum for reconstructed K;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_all", "All PID K^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/pos/contamination_matrix_kaon_pos", "PID K^{+} contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_tpc", "TPC n#sigma_{K} vs transverse momentum for reconstructed K^{+};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_tof", "TOF n#sigma_{K} vs transverse momentum for reconstructed K^{+};#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_tof_matched", "TOF-matched reconstructed K^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/pos/mc_rec_kaon_pos_pt", "True primary K^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/pos/mc_sec_kaon_pos_pt", "Secondary/fake K^{+};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_eta", "Pseudorapidity distribution of reconstructed K^{+};#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed K^{+};#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed K^{+};#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_dcaz", "DCA_{z} vs transverse momentum for reconstructed K^{+};#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for reconstructed K^{+};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_all", "All PID K^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/neg/contamination_matrix_kaon_neg", "PID K^{-} contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_tpc", "TPC n#sigma_{K} vs transverse momentum for reconstructed K^{-};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_tof", "TOF n#sigma_{K} vs transverse momentum for reconstructed K^{-};#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_tof_matched", "TOF-matched reconstructed K^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/neg/mc_rec_kaon_neg_pt", "True primary K^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/neg/mc_sec_kaon_neg_pt", "Secondary/fake K^{-};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_eta", "Pseudorapidity distribution of reconstructed K^{-};#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed K^{-};#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed K^{-};#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_dcaz", "DCA_{z} vs transverse momentum for reconstructed K^{-};#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for reconstructed K^{-};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // Protons
    registryMC.add("mc/matched/reconstruction/protons/rec_proton_all", "All PID p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/contamination_matrix_proton", "PID p contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/mc_rec_proton_pt", "True primary p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/mc_sec_proton_pt", "Secondary/fake p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/proton_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for reconstructed p+#bar{p};#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/matched/reconstruction/protons/proton_tof_beta_vs_p", "TOF #beta vs momentum for reconstructed p+#bar{p};#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryMC.add("mc/matched/reconstruction/protons/rec_proton_eta", "Pseudorapidity distribution of reconstructed p+#bar{p};#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/protons/rec_proton_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed p+#bar{p};#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/protons/rec_proton_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed p+#bar{p};#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/protons/rec_proton_dcaz", "DCA_{z} vs transverse momentum for reconstructed p+#bar{p};#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});

    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_all", "All PID p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/pos/contamination_matrix_proton_pos", "PID p contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_tpc", "TPC n#sigma_{p} vs transverse momentum for reconstructed p;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_tof", "TOF n#sigma_{p} vs transverse momentum for reconstructed p;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_tof_matched", "TOF-matched reconstructed p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/pos/mc_rec_proton_pos_pt", "True primary p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/pos/mc_sec_proton_pos_pt", "Secondary/fake p;#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_eta", "Pseudorapidity distribution of reconstructed p;#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed p;#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed p;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_dcaz", "DCA_{z} vs transverse momentum for reconstructed p;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/matched/reconstruction/protons/pos/rec_proton_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for reconstructed p;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_all", "All PID #bar{p};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/neg/contamination_matrix_proton_neg", "PID #bar{p} contamination matrix;PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_tpc", "TPC n#sigma_{p} vs transverse momentum for reconstructed #bar{p};#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_tof", "TOF n#sigma_{p} vs transverse momentum for reconstructed #bar{p};#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_tof_matched", "TOF-matched reconstructed #bar{p};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/neg/mc_rec_proton_neg_pt", "True primary #bar{p};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/neg/mc_sec_proton_neg_pt", "Secondary/fake #bar{p};#it{p}_{T} (GeV/#it{c});N", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_eta", "Pseudorapidity distribution of reconstructed #bar{p};#eta;N", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_eta_phi", "Azimuthal angle vs pseudorapidity for reconstructed #bar{p};#varphi;#eta", HistType::kTH2F, {axis.axisPhi, axis.axisEta});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_dcaxy", "DCA_{xy} vs transverse momentum for reconstructed #bar{p};#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_dcaz", "DCA_{z} vs transverse momentum for reconstructed #bar{p};#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/matched/reconstruction/protons/neg/rec_proton_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for reconstructed #bar{p};#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // --- MC: MATCHED TRUTH ---

    // Pions
    registryMC.add("mc/matched/truth/pions/mc_gen_pion_pt", "Transverse-momentum distribution of generated primary #pi^{#pm};#it{p}_{T} (GeV/#it{c});N_{#pi}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/truth/pions/pos/mc_gen_pion_pos_pt", "Transverse-momentum distribution of generated primary #pi^{+};#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/truth/pions/neg/mc_gen_pion_neg_pt", "Transverse-momentum distribution of generated primary #pi^{-};#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});

    // Kaons
    registryMC.add("mc/matched/truth/kaons/mc_gen_kaon_pt", "Transverse-momentum distribution of generated primary K^{#pm};#it{p}_{T} (GeV/#it{c});N_{K}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/truth/kaons/pos/mc_gen_kaon_pos_pt", "Transverse-momentum distribution of generated primary K^{+};#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/truth/kaons/neg/mc_gen_kaon_neg_pt", "Transverse-momentum distribution of generated primary K^{-};#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});

    // Protons
    registryMC.add("mc/matched/truth/protons/mc_gen_proton_pt", "Transverse-momentum distribution of generated primary p+#bar{p};#it{p}_{T} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/truth/protons/pos/mc_gen_proton_pos_pt", "Transverse-momentum distribution of generated primary p;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/matched/truth/protons/neg/mc_gen_proton_neg_pt", "Transverse-momentum distribution of generated primary #bar{p};#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});

    // Hadrons
    registryMC.add("mc/matched/truth/hadrons/mc_gen_hadron_pt", "Transverse-momentum distribution of generated primary hadrons;#it{p}_{T} (GeV/#it{c});N_{hadrons}", HistType::kTH1F, {axis.axisPt});

    // --- MC: TRUTH JETS (PARTICLE LEVEL) ---

    // General Jet QA
    registryMC.add("mc/jets/truth/general/n_events", "Selected particle-level jet events;Event counter;N_{events}", HistType::kTH1I, {axis.axisEventCounter});
    registryMC.add("mc/jets/truth/general/n_events_raw", "All events", HistType::kTH1I, {axis.axisEventCounter});
    registryMC.add("mc/jets/truth/general/collision_multiplicity", "Number of particle-level jets in collision;Number of Jets;Number of Collisions", HistType::kTH1F, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/general/jet_pt", "Transverse-momentum distribution of particle-level charged jets;#it{p}_{T,jet}^{part} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetPt});
    registryMC.add("mc/jets/truth/general/jet_eta", "Pseudorapidity distribution of particle-level charged jets;#eta_{jet}^{part};N_{jets}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/truth/general/jet_phi", "Azimuthal-angle distribution of particle-level charged jets;#varphi_{jet}^{part};N_{jets}", HistType::kTH1F, {axis.axisPhi});
    registryMC.add("mc/jets/truth/general/jet_n_constituents", "Constituent multiplicity of particle-level charged jets;N_{constituents}^{part};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});
    registryMC.add("mc/jets/truth/general/jet_area", "Area distribution of particle-level charged jets;A_{jet}^{part};N_{jets}", HistType::kTH1F, {axis.axisJetArea});
    registryMC.add("mc/jets/truth/general/jet_n_selected_constituents", "Selected-particle multiplicity in particle-level charged jets;N_{constituents}^{selected};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});
    registryMC.add("mc/jets/truth/general/z_vtx", "Particle-level Z-vertex distribution;z_{vtx} (cm);N_{events}", HistType::kTH1F, {axis.axisZVtx});
    registryMC.add("mc/jets/truth/general/jet_pt_subtracted", "Particle-level jet #it{p}_{T} after background subtraction;#it{p}_{T,jet}^{part,sub} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetSubPt});
    registryMC.add("mc/jets/truth/general/jet_pt_raw_vs_sub", "Particle-level raw vs subtracted jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet}^{part,sub} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisJetSubPt});

    // Pions
    registryMC.add("mc/jets/truth/pions/n_pions_per_jet", "Multiplicity of generated #pi^{#pm} in particle-level jets;N_{#pi};N_{jets}", HistType::kTH1D, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/pions/n_pions_vs_jet_pt", "Multiplicity of generated #pi^{#pm} vs particle-level jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{#pi}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/pions/pion_pt_vs_jet_pt", "Generated #pi^{#pm} #it{p}_{T} vs particle-level jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/truth/pions/mc_gen_pion_pt", "Transverse-momentum distribution of generated primary #pi^{#pm} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{#pi}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/pions/pos/n_pions_per_jet", "Multiplicity of generated #pi^{+} in particle-level jets;N_{#pi^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/pions/pos/mc_gen_pion_pos_pt", "Transverse-momentum distribution of generated primary #pi^{+} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/pions/neg/n_pions_per_jet", "Multiplicity of generated #pi^{-} in particle-level jets;N_{#pi^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/pions/neg/mc_gen_pion_neg_pt", "Transverse-momentum distribution of generated primary #pi^{-} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});

    // Kaons
    registryMC.add("mc/jets/truth/kaons/n_kaons_per_jet", "Multiplicity of generated K^{#pm} in particle-level jets;N_{K};N_{jets}", HistType::kTH1D, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/kaons/n_kaons_vs_jet_pt", "Multiplicity of generated K^{#pm} vs particle-level jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{K}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/kaons/kaon_pt_vs_jet_pt", "Generated K^{#pm} #it{p}_{T} vs particle-level jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/truth/kaons/mc_gen_kaon_pt", "Transverse-momentum distribution of generated primary K^{#pm} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{K}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/kaons/pos/n_kaons_per_jet", "Multiplicity of generated K^{+} in particle-level jets;N_{K^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/kaons/pos/mc_gen_kaon_pos_pt", "Transverse-momentum distribution of generated primary K^{+} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/kaons/neg/n_kaons_per_jet", "Multiplicity of generated K^{-} in particle-level jets;N_{K^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/kaons/neg/mc_gen_kaon_neg_pt", "Transverse-momentum distribution of generated primary K^{-} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});

    // Protons
    registryMC.add("mc/jets/truth/protons/n_protons_per_jet", "Multiplicity of generated p+#bar{p} in particle-level jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1D, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/protons/n_protons_vs_jet_pt", "Multiplicity of generated p+#bar{p} vs particle-level jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/protons/proton_pt_vs_jet_pt", "Generated p+#bar{p} #it{p}_{T} vs particle-level jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/truth/protons/mc_gen_proton_pt", "Transverse-momentum distribution of generated primary p+#bar{p} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/protons/pos/n_protons_per_jet", "Multiplicity of generated p in particle-level jets;N_{p};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/protons/pos/mc_gen_proton_pos_pt", "Transverse-momentum distribution of generated primary p in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/protons/neg/n_protons_per_jet", "Multiplicity of generated #bar{p} in particle-level jets;N_{#bar{p}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/truth/protons/neg/mc_gen_proton_neg_pt", "Transverse-momentum distribution of generated primary #bar{p} in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/truth/composition/particle_multiplicity_per_jet", "Identified-particle multiplicity in particle-level jets;Particle species;N_{particles} per jet", HistType::kTH2I, {axis.axisParticleSpecies, axis.axisMultiplicity});

    // General Hadrons
    registryMC.add("mc/jets/truth/hadrons/mc_gen_hadron_pt", "Transverse-momentum distribution of generated primary hadrons in particle-level jets;#it{p}_{T} (GeV/#it{c});N_{hadrons}", HistType::kTH1F, {axis.axisPt});

    // --- MC: MATCHED RECONSTRUCTED JETS ---

    // General Jet QA
    registryMC.add("mc/jets/matched/reconstructed/general/n_events", "Selected matched events with reconstructed charged jets;Event counter;N_{events}", HistType::kTH1I, {axis.axisEventCounter});
    registryMC.add("mc/jets/matched/reconstructed/general/n_events_raw", "All events", HistType::kTH1I, {axis.axisEventCounter});
    registryMC.add("mc/jets/matched/reconstructed/general/collision_multiplicity", "Number of matched jets in collision;Number of Jets;Number of Collisions", HistType::kTH1F, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_pt", "Transverse-momentum distribution of matched reconstructed charged jets;#it{p}_{T,jet}^{rec} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetPt});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_pt_subtracted", "Matched reconstructed jet #it{p}_{T} after background subtraction;#it{p}_{T,jet}^{rec,sub} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetSubPt});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_pt_raw_vs_sub", "Raw vs subtracted matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});#it{p}_{T,jet}^{rec,sub} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisJetSubPt});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_eta", "Pseudorapidity distribution of matched reconstructed charged jets;#eta_{jet}^{rec};N_{jets}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_phi", "Azimuthal-angle distribution of matched reconstructed charged jets;#varphi_{jet}^{rec};N_{jets}", HistType::kTH1F, {axis.axisPhi});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_n_constituents", "Constituent multiplicity of matched reconstructed charged jets;N_{constituents}^{rec};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_area", "Area distribution of matched reconstructed charged jets;A_{jet}^{rec};N_{jets}", HistType::kTH1F, {axis.axisJetArea});
    registryMC.add("mc/jets/matched/reconstructed/general/jet_n_selected_constituents", "Selected-track multiplicity in matched reconstructed charged jets;N_{constituents}^{selected};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});

    // Hadrons in jets
    registryMC.add("mc/jets/matched/reconstructed/hadrons/rec_hadron_all", "All reconstructed hadrons in matched jets;#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/hadrons/mc_rec_hadron_pt", "True primary reconstructed hadrons in matched jets;#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/hadrons/mc_sec_hadron_pt", "Secondary/fake reconstructed hadrons in matched jets;#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});

    // Jet Cone
    registryMC.add("mc/jets/matched/reconstructed/cone/all_tracks_2d", "Selected tracks in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/all_tracks_3d", "Selected tracks in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/pions_2d", "Reconstructed #pi in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/pions_3d", "Reconstructed #pi in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/pions/pion_pos_2d", "Reconstructed #pi^{+} in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/pions/pion_pos", "Reconstructed #pi^{+} in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/pions/pion_neg_2d", "Reconstructed #pi^{-} in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/pions/pion_neg", "Reconstructed #pi^{-} in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/kaons_2d", "Reconstructed K in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/kaons_3d", "Reconstructed K in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/kaons/kaon_pos_2d", "Reconstructed K^{+} in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/kaons/kaon_pos", "Reconstructed K^{+} in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/kaons/kaon_neg_2d", "Reconstructed K^{-} in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/kaons/kaon_neg", "Reconstructed K^{-} in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/protons_2d", "Reconstructed p/#bar{p} in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/protons_3d", "Reconstructed p/#bar{p} in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3F, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/protons/proton_pos_2d", "Reconstructed p in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/protons/proton_pos", "Reconstructed p in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/protons/proton_neg_2d", "Reconstructed #bar{p} in MCD matched jet cone;#Delta#eta;#Delta#varphi", HistType::kTH2D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone});
    registryMC.add("mc/jets/matched/reconstructed/cone/protons/proton_neg", "Reconstructed #bar{p} in MCD matched jet cone;#Delta#eta;#Delta#varphi;#it{p}_{T}^{track} (GeV/#it{c})", HistType::kTH3D, {axis.axisDeltaEtaCone, axis.axisDeltaPhiCone, axis.axisPtTrackCone});

    // Composition
    registryMC.add("mc/jets/matched/reconstructed/composition/n_pions_per_jet", "Multiplicity of PID-selected #pi^{#pm} in matched MC jets;N_{#pi};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/composition/n_kaons_per_jet", "Multiplicity of PID-selected K^{#pm} in matched MC jets;N_{K};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/composition/n_protons_per_jet", "Multiplicity of PID-selected p+#bar{p} in matched MC jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/composition/particle_multiplicity_per_jet", "Identified-particle multiplicity in matched MC jets;Particle species;N_{particles} per jet", HistType::kTH2I, {axis.axisParticleSpecies, axis.axisMultiplicity});

    // Truth vs Reconstructed Matrices
    registryMC.add("mc/jets/matched/reconstructed/composition/n_pions_matrix", "Truth vs reconstructed pion multiplicity;N_{truth #pi};N_{reconstructed #pi};N_{jets}", HistType::kTH2I, {axis.axisMultiplicity, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/composition/n_kaons_matrix", "Truth vs reconstructed kaon multiplicity;N_{truth K};N_{reconstructed K};N_{jets}", HistType::kTH2I, {axis.axisMultiplicity, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/composition/n_protons_matrix", "Truth vs reconstructed proton multiplicity;N_{truth p};N_{reconstructed p};N_{jets}", HistType::kTH2I, {axis.axisMultiplicity, axis.axisMultiplicity});

    // Pions
    registryMC.add("mc/jets/matched/reconstructed/pions/n_pions_vs_jet_pt", "Multiplicity of PID-selected #pi^{#pm} vs matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});N_{#pi}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/pions/pion_pt_vs_jet_pt", "Reconstructed #pi^{#pm} #it{p}_{T} vs matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/n_pions_per_jet", "Multiplicity of PID-selected #pi^{+} in matched MC jets;N_{#pi^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/n_pions_per_jet", "Multiplicity of PID-selected #pi^{-} in matched MC jets;N_{#pi^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/pions/rec_pion_all", "All reconstructed pions #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/contamination_matrix_pion", "Pion contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/mc_rec_pion_pt", "MC true primary pions #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/mc_sec_pion_pt", "MC secondary/fake pions #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pion_jet_tpc", "TPC n#sigma_{#pi} vs transverse momentum for PID-selected #pi candidates in matched jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/pions/pion_jet_tof", "TOF n#sigma_{#pi} vs transverse momentum for PID-selected #pi candidates in matched jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/pions/pion_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected pion candidates in matched jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/jets/matched/reconstructed/pions/pion_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected pion candidates in matched jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/rec_pion_pos_all", "All reconstructed #pi^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/contamination_matrix_pion_pos", "#pi^{+} contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_tpc", "#pi^{+} TPC n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC n#sigma_{#pi}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/rec_pion_pos_tof_matched", "TOF-matched reconstructed #pi^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_tof", "#pi^{+} TOF n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF n#sigma_{#pi}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/mc_rec_pion_pos_pt", "MC true primary #pi^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/mc_sec_pion_pos_pt", "MC secondary/fake #pi^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_pt", "Transverse-momentum distribution of #pi^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_eta", "Pseudorapidity distribution of #pi^{+} candidates in matched MC jets;#eta;N_{#pi^{+}}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_dcaxy", "DCA_{xy} vs transverse momentum for #pi^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_dcaz", "DCA_{z} vs transverse momentum for #pi^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for #pi^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/rec_pion_neg_all", "All reconstructed #pi^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/contamination_matrix_pion_neg", "#pi^{-} contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_tpc", "#pi^{-} TPC n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC n#sigma_{#pi}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/rec_pion_neg_tof_matched", "TOF-matched reconstructed #pi^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_tof", "#pi^{-} TOF n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF n#sigma_{#pi}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/mc_rec_pion_neg_pt", "MC true primary #pi^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/mc_sec_pion_neg_pt", "MC secondary/fake #pi^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_pt", "Transverse-momentum distribution of #pi^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_eta", "Pseudorapidity distribution of #pi^{-} candidates in matched MC jets;#eta;N_{#pi^{-}}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_dcaxy", "DCA_{xy} vs transverse momentum for #pi^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_dcaz", "DCA_{z} vs transverse momentum for #pi^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for #pi^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // Kaons
    registryMC.add("mc/jets/matched/reconstructed/kaons/n_kaons_vs_jet_pt", "Multiplicity of PID-selected K^{#pm} vs matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});N_{K}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/kaons/kaon_pt_vs_jet_pt", "Reconstructed K^{#pm} #it{p}_{T} vs matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/n_kaons_per_jet", "Multiplicity of PID-selected K^{+} in matched MC jets;N_{K^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/n_kaons_per_jet", "Multiplicity of PID-selected K^{-} in matched MC jets;N_{K^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/kaons/rec_kaon_all", "All reconstructed kaons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/contamination_matrix_kaon", "Kaon contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/mc_rec_kaon_pt", "MC true primary kaons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/mc_sec_kaon_pt", "MC secondary/fake kaons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/kaon_jet_tpc", "TPC n#sigma_{K} vs transverse momentum for PID-selected K candidates in matched jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/kaons/kaon_jet_tof", "TOF n#sigma_{K} vs transverse momentum for PID-selected K candidates in matched jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/kaons/kaon_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected kaon candidates in matched jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/jets/matched/reconstructed/kaons/kaon_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected kaon candidates in matched jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/rec_kaon_pos_all", "All reconstructed K^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/contamination_matrix_kaon_pos", "K^{+} contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_tpc", "K^{+} TPC n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC n#sigma_{K}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/rec_kaon_pos_tof_matched", "TOF-matched reconstructed K^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_tof", "K^{+} TOF n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF n#sigma_{K}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/mc_rec_kaon_pos_pt", "MC true primary K^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/mc_sec_kaon_pos_pt", "MC secondary/fake K^{+} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_pt", "Transverse-momentum distribution of K^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_eta", "Pseudorapidity distribution of K^{+} candidates in matched MC jets;#eta;N_{K^{+}}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_dcaxy", "DCA_{xy} vs transverse momentum for K^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_dcaz", "DCA_{z} vs transverse momentum for K^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for K^{+} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_all", "All reconstructed K^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/contamination_matrix_kaon_neg", "K^{-} contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_tpc", "K^{-} TPC n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC n#sigma_{K}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_tof_matched", "TOF-matched reconstructed K^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_tof", "K^{-} TOF n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF n#sigma_{K}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/mc_rec_kaon_neg_pt", "MC true primary K^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/mc_sec_kaon_neg_pt", "MC secondary/fake K^{-} #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_pt", "Transverse-momentum distribution of K^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_eta", "Pseudorapidity distribution of K^{-} candidates in matched MC jets;#eta;N_{K^{-}}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_dcaxy", "DCA_{xy} vs transverse momentum for K^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_dcaz", "DCA_{z} vs transverse momentum for K^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for K^{-} candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // Protons
    registryMC.add("mc/jets/matched/reconstructed/protons/n_protons_vs_jet_pt", "Multiplicity of PID-selected p+#bar{p} vs matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/protons/proton_pt_vs_jet_pt", "Reconstructed p+#bar{p} #it{p}_{T} vs matched reconstructed jet #it{p}_{T};#it{p}_{T,jet}^{rec} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/n_protons_per_jet", "Multiplicity of PID-selected p in matched MC jets;N_{p};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/n_protons_per_jet", "Multiplicity of PID-selected #bar{p} in matched MC jets;N_{#bar{p}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/reconstructed/protons/rec_proton_all", "All reconstructed protons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/contamination_matrix_proton", "Proton contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/mc_rec_proton_pt", "MC true primary protons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/mc_sec_proton_pt", "MC secondary/fake protons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/proton_jet_tpc", "TPC n#sigma_{p} vs transverse momentum for PID-selected p+#bar{p} candidates in matched jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/protons/proton_jet_tof", "TOF n#sigma_{p} vs transverse momentum for PID-selected p+#bar{p} candidates in matched jets;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", HistType::kTH2F, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/protons/proton_tpc_signal_vs_p", "TPC d#it{E}/d#it{x} signal vs momentum for PID-selected proton candidates in matched jets;#it{p} (GeV/#it{c});TPC d#it{E}/d#it{x} signal (a.u.)", HistType::kTH2F, {axis.axisMomentum, axis.axisTPCSignal});
    registryMC.add("mc/jets/matched/reconstructed/protons/proton_tof_beta_vs_p", "TOF #beta vs momentum for PID-selected proton candidates in matched jets;#it{p} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {axis.axisMomentum, axis.axisTOFBeta});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/rec_proton_pos_all", "All reconstructed protons (positive) #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/contamination_matrix_proton_pos", "Proton (positive) contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_tpc", "Proton TPC n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC n#sigma_{p}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/rec_proton_pos_tof_matched", "TOF-matched reconstructed protons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_tof", "Proton TOF n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF n#sigma_{p}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/mc_rec_proton_pos_pt", "MC true primary protons (positive) #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/mc_sec_proton_pos_pt", "MC secondary/fake protons (positive) #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_pt", "Transverse-momentum distribution of proton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_eta", "Pseudorapidity distribution of proton candidates in matched MC jets;#eta;N_{p}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_dcaxy", "DCA_{xy} vs transverse momentum for proton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_dcaz", "DCA_{z} vs transverse momentum for proton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_phi_vs_pt", "Azimuthal angle vs transverse momentum for proton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/rec_proton_neg_all", "All reconstructed antiprotons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/contamination_matrix_proton_neg", "Antiproton contamination matrix;True PDG Code;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {axis.axisMCPdg, axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_tpc", "Antiproton TPC n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TPC n#sigma_{p}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTPC});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/rec_proton_neg_tof_matched", "TOF-matched reconstructed antiprotons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_tof", "Antiproton TOF n#sigma vs #it{p}_{T};#it{p}_{T} (GeV/#it{c});TOF n#sigma_{p}", HistType::kTH2D, {axis.axisPt, axis.axisNSigmaTOF});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/mc_rec_proton_neg_pt", "MC true primary antiprotons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/mc_sec_proton_neg_pt", "MC secondary/fake antiprotons #it{p}_{T};#it{p}_{T} (GeV/#it{c});Counts", HistType::kTH1D, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_pt", "Transverse-momentum distribution of antiproton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_eta", "Pseudorapidity distribution of antiproton candidates in matched MC jets;#eta;N_{#bar{p}}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_dcaxy", "DCA_{xy} vs transverse momentum for antiproton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaXY});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_dcaz", "DCA_{z} vs transverse momentum for antiproton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2F, {axis.axisPt, axis.axisDcaZ});
    registryMC.add("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_phi_vs_pt", "Azimuthal angle vs transverse momentum for antiproton candidates in matched MC jets;#it{p}_{T} (GeV/#it{c});#varphi", HistType::kTH2F, {axis.axisPt, axis.axisPhi});

    // --- MC: MATCHED TRUTH JETS ---

    // General Jet QA
    registryMC.add("mc/jets/matched/truth/general/jet_pt", "Transverse-momentum distribution of particle-level matched charged jets;#it{p}_{T,jet}^{part} (GeV/#it{c});N_{jets}", HistType::kTH1F, {axis.axisJetPt});
    registryMC.add("mc/jets/matched/truth/general/jet_eta", "Pseudorapidity distribution of particle-level matched charged jets;#eta_{jet}^{part};N_{jets}", HistType::kTH1F, {axis.axisEta});
    registryMC.add("mc/jets/matched/truth/general/jet_phi", "Azimuthal-angle distribution of particle-level matched charged jets;#varphi_{jet}^{part};N_{jets}", HistType::kTH1F, {axis.axisPhi});
    registryMC.add("mc/jets/matched/truth/general/jet_n_constituents", "Constituent multiplicity of particle-level matched charged jets;N_{constituents}^{part};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});
    registryMC.add("mc/jets/matched/truth/general/jet_area", "Area distribution of particle-level matched charged jets;A_{jet}^{part};N_{jets}", HistType::kTH1F, {axis.axisJetArea});
    registryMC.add("mc/jets/matched/truth/general/jet_n_selected_constituents", "Selected-particle multiplicity in particle-level matched charged jets;N_{constituents}^{selected};N_{jets}", HistType::kTH1I, {axis.axisJetNConst});

    // Pions
    registryMC.add("mc/jets/matched/truth/pions/n_pions_per_jet", "Multiplicity of generated #pi^{#pm} in particle-level matched jets;N_{#pi};N_{jets}", HistType::kTH1D, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/pions/n_pions_vs_jet_pt", "Multiplicity of generated #pi^{#pm} vs particle-level matched jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{#pi}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/pions/pion_pt_vs_jet_pt", "Generated #pi^{#pm} #it{p}_{T} vs particle-level matched jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/matched/truth/pions/mc_gen_pion_pt", "Transverse-momentum distribution of generated primary #pi^{#pm} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{#pi}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/truth/pions/pos/n_pions_per_jet", "Multiplicity of generated #pi^{+} in particle-level matched jets;N_{#pi^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/pions/pos/mc_gen_pion_pos_pt", "Transverse-momentum distribution of generated primary #pi^{+} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/truth/pions/neg/n_pions_per_jet", "Multiplicity of generated #pi^{-} in particle-level matched jets;N_{#pi^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/pions/neg/mc_gen_pion_neg_pt", "Transverse-momentum distribution of generated primary #pi^{-} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{#pi^{-}}", HistType::kTH1F, {axis.axisPt});

    // Kaons
    registryMC.add("mc/jets/matched/truth/kaons/n_kaons_per_jet", "Multiplicity of generated K^{#pm} in particle-level matched jets;N_{K};N_{jets}", HistType::kTH1D, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/kaons/n_kaons_vs_jet_pt", "Multiplicity of generated K^{#pm} vs particle-level matched jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{K}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/kaons/kaon_pt_vs_jet_pt", "Generated K^{#pm} #it{p}_{T} vs particle-level matched jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/matched/truth/kaons/mc_gen_kaon_pt", "Transverse-momentum distribution of generated primary K^{#pm} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{K}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/truth/kaons/pos/n_kaons_per_jet", "Multiplicity of generated K^{+} in particle-level matched jets;N_{K^{+}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/kaons/pos/mc_gen_kaon_pos_pt", "Transverse-momentum distribution of generated primary K^{+} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{K^{+}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/truth/kaons/neg/n_kaons_per_jet", "Multiplicity of generated K^{-} in particle-level matched jets;N_{K^{-}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/kaons/neg/mc_gen_kaon_neg_pt", "Transverse-momentum distribution of generated primary K^{-} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{K^{-}}", HistType::kTH1F, {axis.axisPt});

    // Protons
    registryMC.add("mc/jets/matched/truth/protons/n_protons_per_jet", "Multiplicity of generated p+#bar{p} in particle-level matched jets;N_{p+#bar{p}};N_{jets}", HistType::kTH1D, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/protons/n_protons_vs_jet_pt", "Multiplicity of generated p+#bar{p} vs particle-level matched jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH2I, {axis.axisJetPt, axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/protons/proton_pt_vs_jet_pt", "Generated p+#bar{p} #it{p}_{T} vs particle-level matched jet #it{p}_{T};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {axis.axisJetPt, axis.axisPt});
    registryMC.add("mc/jets/matched/truth/protons/mc_gen_proton_pt", "Transverse-momentum distribution of generated primary p+#bar{p} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{p+#bar{p}}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/truth/protons/pos/n_protons_per_jet", "Multiplicity of generated p in particle-level matched jets;N_{p};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/protons/pos/mc_gen_proton_pos_pt", "Transverse-momentum distribution of generated primary p in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{p}", HistType::kTH1F, {axis.axisPt});
    registryMC.add("mc/jets/matched/truth/protons/neg/n_protons_per_jet", "Multiplicity of generated #bar{p} in particle-level matched jets;N_{#bar{p}};N_{jets}", HistType::kTH1I, {axis.axisMultiplicity});
    registryMC.add("mc/jets/matched/truth/protons/neg/mc_gen_proton_neg_pt", "Transverse-momentum distribution of generated primary #bar{p} in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{#bar{p}}", HistType::kTH1F, {axis.axisPt});

    // General Hadrons
    registryMC.add("mc/jets/matched/truth/hadrons/mc_gen_hadron_pt", "Transverse-momentum distribution of generated primary hadrons in particle-level matched jets;#it{p}_{T} (GeV/#it{c});N_{hadrons}", HistType::kTH1F, {axis.axisPt});
  }

  void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2)
  {
    double const threshold = 1e-9;

    if (p.Pt() < threshold) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }

    double pt = p.Pt();
    double eta = p.Eta();
    double phi = p.Phi();

    u1.SetPtEtaPhi(pt, eta, phi + PIHalf);
    u2.SetPtEtaPhi(pt, eta, phi - PIHalf);
  }

  template <typename TrackIts>
  bool hasITSLayerHit(const TrackIts& track, int layer)
  {
    return TESTBIT(track.itsClusterMap(), layer - 1);
  }

  struct PidResult {
    bool isPion, isKaon, isProton;
  };

  template <typename TrackType>
  PidResult getPid(const TrackType& track)
  {
    constexpr int ClosestMatch = 0;
    constexpr int ExclusiveMatch = 1;
    constexpr int RejectionBased = 2;

    double const buffer = 999.0;
    double pt = track.pt();

    double dForPionsPi = 0;
    double dForPionsKa = 0;
    double dForPionsPr = 0;

    double dForKaonsPi = 0;
    double dForKaonsKa = 0;
    double dForKaonsPr = 0;

    double dForProtonsPi = 0;
    double dForProtonsKa = 0;
    double dForProtonsPr = 0;

    if (pt < cfg.ptThresholdPion) {
      dForPionsPi = std::abs(track.tpcNSigmaPi());
      dForPionsKa = std::abs(track.tpcNSigmaKa());
      dForPionsPr = std::abs(track.tpcNSigmaPr());
    } else {
      if (track.hasTOF()) {
        dForPionsPi = std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi());
        dForPionsKa = std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa());
        dForPionsPr = std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr());
      } else {
        dForPionsPi = buffer;
        dForPionsKa = buffer;
        dForPionsPr = buffer;
      }
    }

    if (pt < cfg.ptThresholdKaon) {
      dForKaonsPi = std::abs(track.tpcNSigmaPi());
      dForKaonsKa = std::abs(track.tpcNSigmaKa());
      dForKaonsPr = std::abs(track.tpcNSigmaPr());
    } else {
      if (track.hasTOF()) {
        dForKaonsPi = std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi());
        dForKaonsKa = std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa());
        dForKaonsPr = std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr());
      } else {
        dForKaonsPi = buffer;
        dForKaonsKa = buffer;
        dForKaonsPr = buffer;
      }
    }

    if (pt < cfg.ptThresholdProton) {
      dForProtonsPi = std::abs(track.tpcNSigmaPi());
      dForProtonsKa = std::abs(track.tpcNSigmaKa());
      dForProtonsPr = std::abs(track.tpcNSigmaPr());
    } else {
      if (track.hasTOF()) {
        dForProtonsPi = std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi());
        dForProtonsKa = std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa());
        dForProtonsPr = std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr());
      } else {
        dForProtonsPi = buffer;
        dForProtonsKa = buffer;
        dForProtonsPr = buffer;
      }
    }

    bool isPiMatch = (dForPionsPi <= cfg.nSigmaCutPi);
    bool isKaMatch = (dForKaonsKa <= cfg.nSigmaCutKa);
    bool isPrMatch = (dForProtonsPr <= cfg.nSigmaCutPr);

    PidResult res{.isPion = false, .isKaon = false, .isProton = false};

    if (cfg.pidMethod == ClosestMatch) {
      if (isPiMatch && dForPionsPi < dForPionsKa && dForPionsPi < dForPionsPr) {
        res.isPion = true;
      }
      if (isKaMatch && dForKaonsKa < dForKaonsPi && dForKaonsKa < dForKaonsPr) {
        res.isKaon = true;
      }
      if (isPrMatch && dForProtonsPr < dForProtonsPi && dForProtonsPr < dForProtonsKa) {
        res.isProton = true;
      }
    } else if (cfg.pidMethod == ExclusiveMatch) {
      if (isPiMatch && !isKaMatch && !isPrMatch) {
        res.isPion = true;
      }
      if (isKaMatch && !isPiMatch && !isPrMatch) {
        res.isKaon = true;
      }
      if (isPrMatch && !isPiMatch && !isKaMatch) {
        res.isProton = true;
      }
    } else if (cfg.pidMethod == RejectionBased) {
      if (isKaMatch && dForKaonsPi > cfg.rejectionSigma && dForKaonsPr > cfg.rejectionSigma) {
        res.isKaon = true;
      }
      if (isPiMatch && dForPionsKa > cfg.rejectionSigma && dForPionsPr > cfg.rejectionSigma) {
        res.isPion = true;
      }
      if (isPrMatch && dForProtonsPi > cfg.rejectionSigma && dForProtonsKa > cfg.rejectionSigma) {
        res.isProton = true;
      }
    }

    if (res.isPion && (pt < cfg.minPtPion || pt > cfg.maxPtPion)) {
      res.isPion = false;
    }
    if (res.isKaon && (pt < cfg.minPtKaon || pt > cfg.maxPtKaon)) {
      res.isKaon = false;
    }
    if (res.isProton && (pt < cfg.minPtProton || pt > cfg.maxPtProton)) {
      res.isProton = false;
    }

    return res;
  }

  template <typename TrackType>
  bool passedTrackSelection(const TrackType& track)
  {
    if (requirePvContributor && !(track.isPVContributor())) {
      return false;
    }
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }
    if ((!hasITSLayerHit(track, 1)) && (!hasITSLayerHit(track, 2)) && (!hasITSLayerHit(track, 3))) {
      return false;
    }
    if (track.itsNCls() < minItsNclusters) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows) {
      return false;
    }
    if (track.tpcChi2NCl() < minChiSquareTpc || track.tpcChi2NCl() > maxChiSquareTpc) {
      return false;
    }
    if (track.itsChi2NCl() > maxChiSquareIts) {
      return false;
    }
    if (track.eta() < minEta || track.eta() > maxEta) {
      return false;
    }
    if (track.pt() < minPt || track.pt() > maxPt) {
      return false;
    }
    if (std::abs(track.dcaXY()) > maxDcaxy || std::abs(track.dcaZ()) > maxDcaz) {
      return false;
    }
    return true;
  }

  void processPureTracks(StandardEvents::iterator const& collision, HadronTracks const& globalTracks)
  {
    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx) {
      return;
    }
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }

    int nSelectedTracks = 0;

    registryData.fill(HIST("data/pure/collisions/z_vertex"), collision.posZ());
    registryData.fill(HIST("data/pure/collisions/x_vertex"), collision.posX());
    registryData.fill(HIST("data/pure/collisions/y_vertex"), collision.posY());
    registryData.fill(HIST("data/pure/collisions/xy_vertex"), collision.posX(), collision.posY());
    registryData.fill(HIST("data/pure/collisions/n_contributors"), collision.numContrib());

    for (auto const& track : globalTracks) {
      if (!passedTrackSelection(track)) {
        continue;
      }

      ++nSelectedTracks;
      double pt = track.pt();
      double eta = track.eta();
      double dcaxy = track.dcaXY();
      double dcaz = track.dcaZ();
      int charge = track.sign();
      double phi = track.phi();

      const double p = track.p();
      registryData.fill(HIST("data/pure/general/tpc_signal_vs_p"), p, track.tpcSignal());
      if (track.hasTOF()) {
        registryData.fill(HIST("data/pure/general/tof_beta_vs_p"), p, track.beta());
      }

      PidResult pid = getPid(track);

      if (pid.isPion) {
        registryData.fill(HIST("data/pure/pions/pion_tpc_signal_vs_p"), p, track.tpcSignal());
        registryData.fill(HIST("data/pure/pions/pion_pure_tpc"), pt, track.tpcNSigmaPi());
        registryData.fill(HIST("data/pure/pions/pion_pure_eta_phi"), phi, eta);
        if (track.hasTOF()) {
          registryData.fill(
            HIST("data/pure/pions/pion_tof_beta_vs_p"),
            p,
            track.beta());
          registryData.fill(HIST("data/pure/pions/pion_pure_tof"), pt, track.tofNSigmaPi());
          registryData.fill(HIST("data/pure/pions/pion_pure_combined_tof_tpc"), pt, std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi()));
        }
        registryData.fill(HIST("data/pure/pions/pion_pure_pt"), pt);
        registryData.fill(HIST("data/pure/pions/pion_pure_eta"), eta);
        registryData.fill(HIST("data/pure/pions/pion_pure_dcaxy"), pt, dcaxy);
        registryData.fill(HIST("data/pure/pions/pion_pure_dcaz"), pt, dcaz);

        if (charge > 0) {
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_tpc"), pt, track.tpcNSigmaPi());
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_tof"), pt, track.tofNSigmaPi());
          }
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_pt"), pt);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_eta"), eta);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/pure/pions/pos/pion_pure_pos_dcaz"), pt, dcaz);
        } else {
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_tpc"), pt, track.tpcNSigmaPi());
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_tof"), pt, track.tofNSigmaPi());
          }
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_pt"), pt);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_eta"), eta);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/pure/pions/neg/pion_pure_neg_dcaz"), pt, dcaz);
        }
      }

      if (pid.isKaon) {
        registryData.fill(HIST("data/pure/kaons/kaon_tpc_signal_vs_p"), track.p(), track.tpcSignal());
        registryData.fill(HIST("data/pure/kaons/kaon_pure_tpc"), pt, track.tpcNSigmaKa());
        registryData.fill(HIST("data/pure/kaons/kaon_pure_eta_phi"), phi, eta);
        if (track.hasTOF()) {
          registryData.fill(
            HIST("data/pure/kaons/kaon_tof_beta_vs_p"),
            p,
            track.beta());
          registryData.fill(HIST("data/pure/kaons/kaon_pure_tof"), pt, track.tofNSigmaKa());
          registryData.fill(HIST("data/pure/kaons/kaon_pure_combined_tof_tpc"), pt, std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa()));
        }
        registryData.fill(HIST("data/pure/kaons/kaon_pure_pt"), pt);
        registryData.fill(HIST("data/pure/kaons/kaon_pure_eta"), eta);
        registryData.fill(HIST("data/pure/kaons/kaon_pure_dcaxy"), pt, dcaxy);
        registryData.fill(HIST("data/pure/kaons/kaon_pure_dcaz"), pt, dcaz);

        if (charge > 0) {
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_tpc"), pt, track.tpcNSigmaKa());
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_tof"), pt, track.tofNSigmaKa());
          }
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_pt"), pt);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_eta"), eta);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/pure/kaons/pos/kaon_pure_pos_dcaz"), pt, dcaz);
        } else {
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_tpc"), pt, track.tpcNSigmaKa());
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_tof"), pt, track.tofNSigmaKa());
          }
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_pt"), pt);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_eta"), eta);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/pure/kaons/neg/kaon_pure_neg_dcaz"), pt, dcaz);
        }
      }

      if (pid.isProton) {
        registryData.fill(HIST("data/pure/protons/proton_tpc_signal_vs_p"), track.p(), track.tpcSignal());
        registryData.fill(HIST("data/pure/protons/proton_pure_tpc"), pt, track.tpcNSigmaPr());
        registryData.fill(HIST("data/pure/protons/proton_pure_eta_phi"), phi, eta);
        if (track.hasTOF()) {
          registryData.fill(
            HIST("data/pure/protons/proton_tof_beta_vs_p"),
            p,
            track.beta());
          registryData.fill(HIST("data/pure/protons/proton_pure_tof"), pt, track.tofNSigmaPr());
          registryData.fill(HIST("data/pure/protons/proton_pure_combined_tof_tpc"), pt, std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr()));
        }
        registryData.fill(HIST("data/pure/protons/proton_pure_pt"), pt);
        registryData.fill(HIST("data/pure/protons/proton_pure_eta"), eta);
        registryData.fill(HIST("data/pure/protons/proton_pure_dcaxy"), pt, dcaxy);
        registryData.fill(HIST("data/pure/protons/proton_pure_dcaz"), pt, dcaz);

        if (charge > 0) {
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_tpc"), pt, track.tpcNSigmaPr());
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_tof"), pt, track.tofNSigmaPr());
          }
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_pt"), pt);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_eta"), eta);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/pure/protons/pos/proton_pure_pos_dcaz"), pt, dcaz);
        } else {
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_tpc"), pt, track.tpcNSigmaPr());
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_phi_vs_pt"), pt, phi);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_eta_phi"), phi, eta);

          if (track.hasTOF()) {
            registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_tof"), pt, track.tofNSigmaPr());
          }
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_pt"), pt);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_eta"), eta);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/pure/protons/neg/proton_pure_neg_dcaz"), pt, dcaz);
        }
      }
    }
    registryData.fill(HIST("data/jets/general/total_tracks_pure_vs_jets"), 1.0, static_cast<double>(nSelectedTracks));
  }
  PROCESS_SWITCH(JetHadronsPid, processPureTracks, "Pure Tracks Analysis", true);

  void processMC(soa::Join<StandardEvents, aod::McCollisionLabels>::iterator const& collision, HadronTracksMC const& tracks, aod::McParticles const& particles)
  {
    if (!collision.sel8() || std::abs(collision.posZ()) > zVtx) {
      return;
    }
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }

    registryMC.fill(HIST("mc/matched/general/n_events"), 1);
    registryMC.fill(HIST("mc/matched/general/z_vertex"), collision.posZ());
    registryMC.fill(HIST("mc/matched/general/x_vertex"), collision.posX());
    registryMC.fill(HIST("mc/matched/general/y_vertex"), collision.posY());
    registryMC.fill(HIST("mc/matched/general/xy_vertex"), collision.posX(), collision.posY());
    registryMC.fill(HIST("mc/matched/general/n_contributors"), collision.numContrib());

    if (!collision.has_mcCollision()) {
      return;
    }

    auto particlesPerCollision = particles.sliceBy(mcParticlesPerCollision, collision.mcCollisionId());

    for (auto const& track : tracks) {

      if (!passedTrackSelection(track)) {
        continue;
      }

      double pt = track.pt();
      double eta = track.eta();
      double phi = track.phi();
      double dcaxy = track.dcaXY();
      double dcaz = track.dcaZ();
      const double p = track.p();

      if (!track.has_mcParticle()) {
        continue;
      }
      auto const& trueParticle = track.mcParticle();
      int pdg = std::abs(trueParticle.pdgCode());

      int realpdg = trueParticle.pdgCode();
      bool isPrimary = trueParticle.isPhysicalPrimary();
      int charge = track.sign();

      registryMC.fill(HIST("mc/matched/reconstruction/general/tpc_signal_vs_p"), p, track.tpcSignal());
      if (track.hasTOF()) {
        registryMC.fill(HIST("mc/matched/reconstruction/general/tof_beta_vs_p"), p, track.beta());
      }

      registryMC.fill(HIST("mc/matched/reconstruction/hadrons/rec_hadron_all"), pt);

      if (track.hasTOF()) {
        registryMC.fill(HIST("mc/matched/reconstruction/hadrons/rec_hadron_tof_match"), pt);
      }

      if (isPrimary) {
        registryMC.fill(HIST("mc/matched/reconstruction/hadrons/mc_rec_hadron_pt"), pt);
      } else {
        registryMC.fill(HIST("mc/matched/reconstruction/hadrons/mc_sec_hadron_pt"), pt);
      }

      PidResult pid = getPid(track);

      if (pid.isPion) {
        registryMC.fill(HIST("mc/matched/reconstruction/pions/rec_pion_all"), pt);
        registryMC.fill(HIST("mc/matched/reconstruction/pions/contamination_matrix_pion"), realpdg, pt);
        registryMC.fill(HIST("mc/matched/reconstruction/pions/pion_tpc_signal_vs_p"), p, track.tpcSignal());
        registryMC.fill(HIST("mc/matched/reconstruction/pions/rec_pion_eta"), eta);
        registryMC.fill(HIST("mc/matched/reconstruction/pions/rec_pion_eta_phi"), phi, eta);
        registryMC.fill(HIST("mc/matched/reconstruction/pions/rec_pion_dcaxy"), pt, dcaxy);
        registryMC.fill(HIST("mc/matched/reconstruction/pions/rec_pion_dcaz"), pt, dcaz);
        if (track.hasTOF()) {
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pion_tof_beta_vs_p"), p, track.beta());
        }

        if (charge > 0) {

          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_all"), pt);

          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_tof_matched"), pt);
          }
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/contamination_matrix_pion_pos"), realpdg, pt);

          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_tof"), pt, track.tofNSigmaPi());
          }
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_eta"), eta);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_eta_phi"), phi, eta);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_dcaxy"), pt, dcaxy);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_dcaz"), pt, dcaz);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/rec_pion_pos_phi_vs_pt"), pt, phi);
        } else {
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_all"), pt);

          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_tof_matched"), pt);
          }
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/contamination_matrix_pion_neg"), realpdg, pt);

          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_tof"), pt, track.tofNSigmaPi());
          }
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_eta"), eta);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_eta_phi"), phi, eta);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_dcaxy"), pt, dcaxy);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_dcaz"), pt, dcaz);
          registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/rec_pion_neg_phi_vs_pt"), pt, phi);
        }

        if (isPrimary) {
          registryMC.fill(HIST("mc/matched/reconstruction/pions/mc_rec_pion_pt"), pt);
          if (charge > 0) {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/mc_rec_pion_pos_pt"), pt);
          } else {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/mc_rec_pion_neg_pt"), pt);
          }
        } else {
          if (pdg == PDG_t::kPiPlus) {
            registryMC.fill(HIST("mc/matched/reconstruction/pions/mc_sec_pion_pt"), pt);
            if (charge > 0) {
              registryMC.fill(HIST("mc/matched/reconstruction/pions/pos/mc_sec_pion_pos_pt"), pt);
            } else {
              registryMC.fill(HIST("mc/matched/reconstruction/pions/neg/mc_sec_pion_neg_pt"), pt);
            }
          }
        }
      }

      if (pid.isKaon) {
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/rec_kaon_all"), pt);
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/contamination_matrix_kaon"), realpdg, pt);
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/kaon_tpc_signal_vs_p"), p, track.tpcSignal());
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/rec_kaon_eta"), eta);
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/rec_kaon_eta_phi"), phi, eta);
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/rec_kaon_dcaxy"), pt, dcaxy);
        registryMC.fill(HIST("mc/matched/reconstruction/kaons/rec_kaon_dcaz"), pt, dcaz);
        if (track.hasTOF()) {
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/kaon_tof_beta_vs_p"), p, track.beta());
        }

        if (charge > 0) {

          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_all"), pt);

          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_tof_matched"), pt);
          }
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/contamination_matrix_kaon_pos"), realpdg, pt);

          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_tof"), pt, track.tofNSigmaKa());
          }
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_eta"), eta);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_eta_phi"), phi, eta);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_dcaxy"), pt, dcaxy);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_dcaz"), pt, dcaz);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/rec_kaon_pos_phi_vs_pt"), pt, phi);
        } else {

          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_all"), pt);

          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_tof_matched"), pt);
          }
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/contamination_matrix_kaon_neg"), realpdg, pt);

          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_tof"), pt, track.tofNSigmaKa());
          }
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_eta"), eta);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_eta_phi"), phi, eta);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_dcaxy"), pt, dcaxy);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_dcaz"), pt, dcaz);
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/rec_kaon_neg_phi_vs_pt"), pt, phi);
        }

        if (isPrimary) {
          registryMC.fill(HIST("mc/matched/reconstruction/kaons/mc_rec_kaon_pt"), pt);
          if (charge > 0) {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/mc_rec_kaon_pos_pt"), pt);
          } else {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/mc_rec_kaon_neg_pt"), pt);
          }
        } else {
          if (pdg == PDG_t::kKPlus) {
            registryMC.fill(HIST("mc/matched/reconstruction/kaons/mc_sec_kaon_pt"), pt);
            if (charge > 0) {
              registryMC.fill(HIST("mc/matched/reconstruction/kaons/pos/mc_sec_kaon_pos_pt"), pt);
            } else {
              registryMC.fill(HIST("mc/matched/reconstruction/kaons/neg/mc_sec_kaon_neg_pt"), pt);
            }
          }
        }
      }

      if (pid.isProton) {
        registryMC.fill(HIST("mc/matched/reconstruction/protons/rec_proton_all"), pt);
        registryMC.fill(HIST("mc/matched/reconstruction/protons/contamination_matrix_proton"), realpdg, pt);
        registryMC.fill(HIST("mc/matched/reconstruction/protons/proton_tpc_signal_vs_p"), p, track.tpcSignal());
        registryMC.fill(HIST("mc/matched/reconstruction/protons/rec_proton_eta"), eta);
        registryMC.fill(HIST("mc/matched/reconstruction/protons/rec_proton_eta_phi"), phi, eta);
        registryMC.fill(HIST("mc/matched/reconstruction/protons/rec_proton_dcaxy"), pt, dcaxy);
        registryMC.fill(HIST("mc/matched/reconstruction/protons/rec_proton_dcaz"), pt, dcaz);
        if (track.hasTOF()) {
          registryMC.fill(HIST("mc/matched/reconstruction/protons/proton_tof_beta_vs_p"), p, track.beta());
        }

        if (charge > 0) {
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_all"), pt);
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_tof_matched"), pt);
          }
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/contamination_matrix_proton_pos"), realpdg, pt);

          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_tof"), pt, track.tofNSigmaPr());
          }
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_eta"), eta);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_eta_phi"), phi, eta);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_dcaxy"), pt, dcaxy);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_dcaz"), pt, dcaz);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/rec_proton_pos_phi_vs_pt"), pt, phi);
        } else {

          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_all"), pt);
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_tof_matched"), pt);
          }
          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/contamination_matrix_proton_neg"), realpdg, pt);

          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_tof"), pt, track.tofNSigmaPr());
          }
          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_eta"), eta);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_eta_phi"), phi, eta);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_dcaxy"), pt, dcaxy);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_dcaz"), pt, dcaz);
          registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/rec_proton_neg_phi_vs_pt"), pt, phi);
        }

        if (isPrimary) {
          registryMC.fill(HIST("mc/matched/reconstruction/protons/mc_rec_proton_pt"), pt);
          if (charge > 0) {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/mc_rec_proton_pos_pt"), pt);
          } else {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/mc_rec_proton_neg_pt"), pt);
          }
        } else {
          if (pdg == PDG_t::kProton) {
            registryMC.fill(HIST("mc/matched/reconstruction/protons/mc_sec_proton_pt"), pt);
            if (charge > 0) {
              registryMC.fill(HIST("mc/matched/reconstruction/protons/pos/mc_sec_proton_pos_pt"), pt);
            } else {
              registryMC.fill(HIST("mc/matched/reconstruction/protons/neg/mc_sec_proton_neg_pt"), pt);
            }
          }
        }
      }
    }

    for (auto const& mcpart : particlesPerCollision) {

      int realpdg = mcpart.pdgCode();
      int pdg = std::abs(realpdg);
      bool isPrimary = mcpart.isPhysicalPrimary();
      double ptTrue = mcpart.pt();

      if (!isPrimary) {
        continue;
      }

      if (mcpart.eta() < minEta || mcpart.eta() > maxEta) {
        continue;
      }
      if (ptTrue < minPt || ptTrue > maxPt) {
        continue;
      }

      registryMC.fill(HIST("mc/matched/truth/hadrons/mc_gen_hadron_pt"), ptTrue);

      if (pdg == PDG_t::kPiPlus) {
        if (ptTrue >= cfg.minPtPion && ptTrue <= cfg.maxPtPion) {
          registryMC.fill(HIST("mc/matched/truth/pions/mc_gen_pion_pt"), ptTrue);
          if (realpdg > 0) {
            registryMC.fill(HIST("mc/matched/truth/pions/pos/mc_gen_pion_pos_pt"), ptTrue);
          } else {
            registryMC.fill(HIST("mc/matched/truth/pions/neg/mc_gen_pion_neg_pt"), ptTrue);
          }
        }
      } else if (pdg == PDG_t::kKPlus) {
        if (ptTrue >= cfg.minPtKaon && ptTrue <= cfg.maxPtKaon) {
          registryMC.fill(HIST("mc/matched/truth/kaons/mc_gen_kaon_pt"), ptTrue);
          if (realpdg > 0) {
            registryMC.fill(HIST("mc/matched/truth/kaons/pos/mc_gen_kaon_pos_pt"), ptTrue);
          } else {
            registryMC.fill(HIST("mc/matched/truth/kaons/neg/mc_gen_kaon_neg_pt"), ptTrue);
          }
        }
      } else if (pdg == PDG_t::kProton) {
        if (ptTrue >= cfg.minPtProton && ptTrue <= cfg.maxPtProton) {
          registryMC.fill(HIST("mc/matched/truth/protons/mc_gen_proton_pt"), ptTrue);
          if (realpdg > 0) {
            registryMC.fill(HIST("mc/matched/truth/protons/pos/mc_gen_proton_pos_pt"), ptTrue);
          } else {
            registryMC.fill(HIST("mc/matched/truth/protons/neg/mc_gen_proton_neg_pt"), ptTrue);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetHadronsPid, processMC, "Run on Monte Carlo", false);

  void processMCTruth(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {

    if (std::abs(mcCollision.posZ()) > zVtx) {
      return;
    }

    registryMC.fill(HIST("mc/truth/general/n_events"), 1);
    registryMC.fill(HIST("mc/truth/general/z_vtx"), mcCollision.posZ());

    for (auto const& mcpart : mcParticles) {

      if (!mcpart.isPhysicalPrimary()) {
        continue;
      }

      if (mcpart.eta() < minEta || mcpart.eta() > maxEta) {
        continue;
      }
      if (mcpart.pt() < minPt || mcpart.pt() > maxPt) {
        continue;
      }

      int originalPdg = mcpart.pdgCode();
      int pdg = std::abs(originalPdg);
      double pt = mcpart.pt();

      registryMC.fill(HIST("mc/truth/hadrons/mc_gen_hadron_pt"), pt);

      if (pdg == PDG_t::kPiPlus) {
        if (pt >= cfg.minPtPion && pt <= cfg.maxPtPion) {
          registryMC.fill(HIST("mc/truth/pions/mc_gen_pion_pt"), pt);
          if (originalPdg > 0) {
            registryMC.fill(HIST("mc/truth/pions/pos/mc_gen_pion_pos_pt"), pt);
          } else {
            registryMC.fill(HIST("mc/truth/pions/neg/mc_gen_pion_neg_pt"), pt);
          }
        }
      } else if (pdg == PDG_t::kKPlus) {
        if (pt >= cfg.minPtKaon && pt <= cfg.maxPtKaon) {
          registryMC.fill(HIST("mc/truth/kaons/mc_gen_kaon_pt"), pt);
          if (originalPdg > 0) {
            registryMC.fill(HIST("mc/truth/kaons/pos/mc_gen_kaon_pos_pt"), pt);
          } else {
            registryMC.fill(HIST("mc/truth/kaons/neg/mc_gen_kaon_neg_pt"), pt);
          }
        }
      } else if (pdg == PDG_t::kProton) {
        if (pt >= cfg.minPtProton && pt <= cfg.maxPtProton) {
          registryMC.fill(HIST("mc/truth/protons/mc_gen_proton_pt"), pt);
          if (originalPdg > 0) {
            registryMC.fill(HIST("mc/truth/protons/pos/mc_gen_proton_pos_pt"), pt);
          } else {
            registryMC.fill(HIST("mc/truth/protons/neg/mc_gen_proton_neg_pt"), pt);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetHadronsPid, processMCTruth, "Run on Monte Carlo (Pure Truth)", false);

  void processJets(JetEvents::iterator const& collision,
                   soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                   soa::Join<aod::JetTracks, aod::JTrackPIs> const&,
                   HadronTracks const& globalTracks)
  {
    registryData.fill(HIST("data/jets/general/n_events_raw"), 1);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    float zVertex = collision.posZ();
    if (std::abs(zVertex) > zVtx) {
      return;
    }

    registryData.fill(HIST("data/jets/general/n_events"), 1);
    registryData.fill(HIST("data/jets/general/z_vtx"), zVertex);

    double centralRho = collision.rho();
    int jetsInCollision = 0;

    auto collTracks = globalTracks.sliceBy(tracksPerCollision, collision.collisionId());
    int tracksInJetsCollision = 0;

    for (auto const& jet : jets) {

      if (!isppRefAnalysis && ((std::abs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))) {
        continue;
      }
      if (isppRefAnalysis && std::abs(jet.eta()) > cfgEtaJetMax) {
        continue;
      }

      if (isppRefAnalysis && (jet.pt() < minJetPt || jet.pt() > maxJetPt)) {
        continue;
      }

      double ptSub = jet.pt() - (centralRho * jet.area());
      if (ptSub < 0) {
        ptSub = 0.0;
      }

      if (!isppRefAnalysis && (ptSub < minJetPt || ptSub > maxJetPt)) {
        continue;
      }

      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && normalizedJetArea < minNormalizedJetArea) {
        continue;
      }

      jetsInCollision++;

      registryData.fill(HIST("data/jets/general/jet_pt_subtracted"), ptSub);
      registryData.fill(HIST("data/jets/general/jet_pt_raw_vs_sub"), jet.pt(), ptSub);

      registryData.fill(HIST("data/jets/general/jet_eta"), jet.eta());
      registryData.fill(HIST("data/jets/general/jet_phi"), jet.phi());
      registryData.fill(HIST("data/jets/general/jet_area"), jet.area());
      registryData.fill(HIST("data/jets/general/jet_pt"), jet.pt());
      registryData.fill(HIST("data/jets/general/jet_n_constituents"), jet.tracksIds().size());

      const double magnitudeThreshold = 1e-9;
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);

      if (ueAxis1.Mag() < magnitudeThreshold || ueAxis2.Mag() < magnitudeThreshold) {
        continue;
      }

      int selectedConstituentCount = 0;

      int nPionsInJet = 0;
      int nKaonsInJet = 0;
      int nProtonsInJet = 0;

      int nPionsInJetPos = 0;
      int nKaonsInJetPos = 0;
      int nProtonsInJetPos = 0;

      int nPionsInJetNeg = 0;
      int nKaonsInJetNeg = 0;
      int nProtonsInJetNeg = 0;

      for (auto const& jtrack : jet.tracks_as<soa::Join<aod::JetTracks, aod::JTrackPIs>>()) {

        auto track = jtrack.track_as<HadronTracks>();

        if (!passedTrackSelection(track)) {
          continue;
        }
        tracksInJetsCollision++;

        selectedConstituentCount++;
        const double p = track.p();
        registryData.fill(HIST("data/jets/general/tpc_signal_vs_p"), p, track.tpcSignal());
        if (track.hasTOF()) {
          registryData.fill(HIST("data/jets/general/tof_beta_vs_p"), p, track.beta());
        }

        double pt = track.pt();
        double eta = track.eta();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();
        int charge = track.sign();
        double phi = track.phi();

        double deltaEtaJet = track.eta() - jet.eta();
        double deltaPhiJet = RecoDecay::constrainAngle(track.phi() - jet.phi(), -PI);

        registryData.fill(HIST("data/jets/cone/all_tracks"), deltaEtaJet, deltaPhiJet, pt);
        registryData.fill(HIST("data/jets/cone/all_tracks_2d"), deltaEtaJet, deltaPhiJet);

        PidResult pid = getPid(track);

        if (pid.isPion) {

          registryData.fill(HIST("data/jets/pions/pion_tpc_signal_vs_p"), p, track.tpcSignal());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/jets/pions/pion_jet_tof"), pt, track.tofNSigmaPi());
            registryData.fill(HIST("data/jets/pions/pion_tof_beta_vs_p"), p, track.beta());
          }
          registryData.fill(HIST("data/jets/pions/pion_jet_tpc"), pt, track.tpcNSigmaPi());
          registryData.fill(HIST("data/jets/pions/pion_pt_vs_jet_pt"), jet.pt(), pt);

          if (charge > 0) {
            nPionsInJetPos++;

            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_tpc"), pt, track.tpcNSigmaPi());
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/cone/pions/pion_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/pions/pion_pos_2d"), deltaEtaJet, deltaPhiJet);

            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_pt"), pt);
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_eta"), eta);
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/jets/pions/pos/pion_jet_pos_dcaz"), pt, dcaz);
          } else {
            nPionsInJetNeg++;

            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_tpc"), pt, track.tpcNSigmaPi());
            registryData.fill(HIST("data/jets/cone/pions/pion_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/pions/pion_neg_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_pt"), pt);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_eta"), eta);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/jets/pions/neg/pion_jet_neg_dcaz"), pt, dcaz);
          }
        }

        if (pid.isKaon) {

          registryData.fill(HIST("data/jets/kaons/kaon_tpc_signal_vs_p"), p, track.tpcSignal());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/jets/kaons/kaon_jet_tof"), pt, track.tofNSigmaKa());
            registryData.fill(HIST("data/jets/kaons/kaon_tof_beta_vs_p"), p, track.beta());
          }
          registryData.fill(HIST("data/jets/kaons/kaon_jet_tpc"), pt, track.tpcNSigmaKa());
          registryData.fill(HIST("data/jets/kaons/kaon_pt_vs_jet_pt"), jet.pt(), pt);

          if (charge > 0) {
            nKaonsInJetPos++;

            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_tpc"), pt, track.tpcNSigmaKa());
            registryData.fill(HIST("data/jets/cone/kaons/kaon_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/kaons/kaon_pos_2d"), deltaEtaJet, deltaPhiJet);

            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_pt"), pt);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_eta"), eta);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/jets/kaons/pos/kaon_jet_pos_dcaz"), pt, dcaz);
          } else {
            nKaonsInJetNeg++;

            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_tpc"), pt, track.tpcNSigmaKa());
            registryData.fill(HIST("data/jets/cone/kaons/kaon_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/kaons/kaon_neg_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_pt"), pt);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_eta"), eta);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/jets/kaons/neg/kaon_jet_neg_dcaz"), pt, dcaz);
          }
        }

        if (pid.isProton) {

          registryData.fill(HIST("data/jets/protons/proton_tpc_signal_vs_p"), p, track.tpcSignal());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/jets/protons/proton_jet_tof"), pt, track.tofNSigmaPr());
            registryData.fill(HIST("data/jets/protons/proton_tof_beta_vs_p"), p, track.beta());
          }
          registryData.fill(HIST("data/jets/protons/proton_jet_tpc"), pt, track.tpcNSigmaPr());
          registryData.fill(HIST("data/jets/protons/proton_pt_vs_jet_pt"), jet.pt(), pt);

          if (charge > 0) {
            nProtonsInJetPos++;

            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_tpc"), pt, track.tpcNSigmaPr());
            registryData.fill(HIST("data/jets/cone/protons/proton_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/protons/proton_pos_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_pt"), pt);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_eta"), eta);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/jets/protons/pos/proton_jet_pos_dcaz"), pt, dcaz);
          } else {
            nProtonsInJetNeg++;
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_phi_vs_pt"), pt, phi);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_tpc"), pt, track.tpcNSigmaPr());
            registryData.fill(HIST("data/jets/cone/protons/proton_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryData.fill(HIST("data/jets/cone/protons/proton_neg_2d"), deltaEtaJet, deltaPhiJet);
            if (track.hasTOF()) {
              registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_pt"), pt);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_eta"), eta);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/jets/protons/neg/proton_jet_neg_dcaz"), pt, dcaz);
          }
        }
      }
      registryData.fill(HIST("data/jets/general/jet_n_selected_constituents"), selectedConstituentCount);

      registryData.fill(HIST("data/jets/pions/pos/n_pions_per_jet"), nPionsInJetPos);
      registryData.fill(HIST("data/jets/kaons/pos/n_kaons_per_jet"), nKaonsInJetPos);
      registryData.fill(HIST("data/jets/protons/pos/n_protons_per_jet"), nProtonsInJetPos);

      registryData.fill(HIST("data/jets/pions/neg/n_pions_per_jet"), nPionsInJetNeg);
      registryData.fill(HIST("data/jets/kaons/neg/n_kaons_per_jet"), nKaonsInJetNeg);
      registryData.fill(HIST("data/jets/protons/neg/n_protons_per_jet"), nProtonsInJetNeg);

      nPionsInJet = nPionsInJetPos + nPionsInJetNeg;
      nProtonsInJet = nProtonsInJetNeg + nProtonsInJetPos;
      nKaonsInJet = nKaonsInJetNeg + nKaonsInJetPos;

      registryData.fill(HIST("data/jets/pions/n_pions_per_jet"), nPionsInJet);
      registryData.fill(HIST("data/jets/kaons/n_kaons_per_jet"), nKaonsInJet);
      registryData.fill(HIST("data/jets/protons/n_protons_per_jet"), nProtonsInJet);

      registryData.fill(HIST("data/jets/pions/n_pions_vs_jet_pt"), jet.pt(), nPionsInJet);
      registryData.fill(HIST("data/jets/kaons/n_kaons_vs_jet_pt"), jet.pt(), nKaonsInJet);
      registryData.fill(HIST("data/jets/protons/n_protons_vs_jet_pt"), jet.pt(), nProtonsInJet);

      registryData.fill(HIST("data/jets/composition/particle_multiplicity_per_jet"), 1.0, nPionsInJet);
      registryData.fill(HIST("data/jets/composition/particle_multiplicity_per_jet"), 2.0, nKaonsInJet);
      registryData.fill(HIST("data/jets/composition/particle_multiplicity_per_jet"), 3.0, nProtonsInJet);

      int nTracksOut = 0;

      for (auto const& track : collTracks) {

        if (!passedTrackSelection(track)) {
          continue;
        }

        double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
        double deltaPhiUe1 = RecoDecay::constrainAngle(track.phi() - ueAxis1.Phi(), -PI);
        double deltaRUe1 = std::hypot(deltaEtaUe1, deltaPhiUe1);

        double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
        double deltaPhiUe2 = RecoDecay::constrainAngle(track.phi() - ueAxis2.Phi(), -PI);
        double deltaRUe2 = std::hypot(deltaEtaUe2, deltaPhiUe2);

        if (deltaRUe1 > rJet && deltaRUe2 > rJet) {
          continue;
        }

        double pt = track.pt();
        double eta = track.eta();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();
        int charge = track.sign();

        PidResult pid = getPid(track);
        nTracksOut++;

        if (pid.isPion) {
          registryData.fill(HIST("data/ue/pions/pion_ue_tpc"), pt, track.tpcNSigmaPi());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/ue/pions/pion_ue_tof"), pt, track.tofNSigmaPi());
          }

          registryData.fill(HIST("data/ue/pions/pion_ue_pt"), pt);
          registryData.fill(HIST("data/ue/pions/pion_ue_eta"), eta);
          registryData.fill(HIST("data/ue/pions/pion_ue_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/ue/pions/pion_ue_dcaz"), pt, dcaz);

          if (charge > 0) {
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_tpc"), pt, track.tpcNSigmaPi());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_pt"), pt);
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_eta"), eta);
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/ue/pions/pos/pion_ue_pos_dcaz"), pt, dcaz);
          } else {
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_tpc"), pt, track.tpcNSigmaPi());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_tof"), pt, track.tofNSigmaPi());
            }
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_pt"), pt);
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_eta"), eta);
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/ue/pions/neg/pion_ue_neg_dcaz"), pt, dcaz);
          }
        }

        if (pid.isKaon) {
          registryData.fill(HIST("data/ue/kaons/kaon_ue_tpc"), pt, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/ue/kaons/kaon_ue_tof"), pt, track.tofNSigmaKa());
          }

          registryData.fill(HIST("data/ue/kaons/kaon_ue_pt"), pt);
          registryData.fill(HIST("data/ue/kaons/kaon_ue_eta"), eta);
          registryData.fill(HIST("data/ue/kaons/kaon_ue_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/ue/kaons/kaon_ue_dcaz"), pt, dcaz);

          if (charge > 0) {
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_tpc"), pt, track.tpcNSigmaKa());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_pt"), pt);
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_eta"), eta);
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/ue/kaons/pos/kaon_ue_pos_dcaz"), pt, dcaz);
          } else {
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_tpc"), pt, track.tpcNSigmaKa());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_tof"), pt, track.tofNSigmaKa());
            }
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_pt"), pt);
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_eta"), eta);
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/ue/kaons/neg/kaon_ue_neg_dcaz"), pt, dcaz);
          }
        }

        if (pid.isProton) {
          registryData.fill(HIST("data/ue/protons/proton_ue_tpc"), pt, track.tpcNSigmaPr());
          if (track.hasTOF()) {
            registryData.fill(HIST("data/ue/protons/proton_ue_tof"), pt, track.tofNSigmaPr());
          }

          registryData.fill(HIST("data/ue/protons/proton_ue_pt"), pt);
          registryData.fill(HIST("data/ue/protons/proton_ue_eta"), eta);
          registryData.fill(HIST("data/ue/protons/proton_ue_dcaxy"), pt, dcaxy);
          registryData.fill(HIST("data/ue/protons/proton_ue_dcaz"), pt, dcaz);

          if (charge > 0) {
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_tpc"), pt, track.tpcNSigmaPr());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_pt"), pt);
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_eta"), eta);
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/ue/protons/pos/proton_ue_pos_dcaz"), pt, dcaz);
          } else {
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_tpc"), pt, track.tpcNSigmaPr());
            if (track.hasTOF()) {
              registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_tof"), pt, track.tofNSigmaPr());
            }
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_pt"), pt);
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_eta"), eta);
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_dcaxy"), pt, dcaxy);
            registryData.fill(HIST("data/ue/protons/neg/proton_ue_neg_dcaz"), pt, dcaz);
          }
        }
      }
      registryData.fill(HIST("data/ue/general/tracks_n_in_ue"), nTracksOut);
    }
    registryData.fill(HIST("data/jets/general/total_tracks_pure_vs_jets"), 2.0, static_cast<double>(tracksInJetsCollision));
    registryData.fill(HIST("data/jets/general/collision_multiplicity"), jetsInCollision);
  }
  PROCESS_SWITCH(JetHadronsPid, processJets, "Jets Analysis", true);

  void processMatchJets(JetEvents::iterator const& collision,
                        soa::Join<ChargedMCDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& detectorLevelJets,
                        soa::Join<aod::JetTracks, aod::JTrackPIs> const&,
                        ChargedMCPJets const&,
                        aod::JetParticles const&,
                        HadronTracksMC const&, aod::McParticles const&)
  {

    registryMC.fill(HIST("mc/jets/matched/reconstructed/general/n_events_raw"), 1);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    if (std::abs(collision.posZ()) > zVtx) {
      return;
    }

    registryMC.fill(HIST("mc/jets/matched/reconstructed/general/n_events"), 1);

    double centralRho = collision.rho();
    int jetsInCollision = 0;

    for (auto const& djet : detectorLevelJets) {

      if (!isppRefAnalysis && ((std::abs(djet.eta()) + rJet) > (maxEta - deltaEtaEdge))) {
        continue;
      }

      if (isppRefAnalysis && std::abs(djet.eta()) > cfgEtaJetMax) {
        continue;
      }

      if (isppRefAnalysis && (djet.pt() < minJetPt || djet.pt() > maxJetPt)) {
        continue;
      }

      double ptSub = djet.pt() - (centralRho * djet.area());
      if (ptSub < 0) {
        ptSub = 0.0;
      }

      if (!isppRefAnalysis && (ptSub < minJetPt || ptSub > maxJetPt)) {
        continue;
      }

      const double normalizedJetArea = djet.area() / (PI * rJet * rJet);

      if (applyAreaCut && normalizedJetArea < minNormalizedJetArea) {
        continue;
      }

      if (!djet.has_matchedJetGeo() || !djet.has_matchedJetPt()) {
        continue;
      }

      auto geoMatches = djet.matchedJetGeo_as<ChargedMCPJets>();
      auto ptMatches = djet.matchedJetPt_as<ChargedMCPJets>();
      auto matchedIt = geoMatches.end();

      for (auto itGeo = geoMatches.begin(); itGeo != geoMatches.end(); ++itGeo) {
        bool found = false;
        for (auto itPt = ptMatches.begin(); itPt != ptMatches.end(); ++itPt) {
          if (itGeo->globalIndex() == itPt->globalIndex()) {
            matchedIt = itGeo;
            found = true;
            break;
          }
        }
        if (found) {
          break;
        }
      }

      if (matchedIt == geoMatches.end()) {
        continue;
      }

      auto const& pjet = *matchedIt;

      ++jetsInCollision;

      int nSelectedConstituents = 0;

      int nPionsInJet = 0;
      int nKaonsInJet = 0;
      int nProtonsInJet = 0;

      int nPionsInJetPos = 0;
      int nKaonsInJetPos = 0;
      int nProtonsInJetPos = 0;

      int nPionsInJetNeg = 0;
      int nKaonsInJetNeg = 0;
      int nProtonsInJetNeg = 0;

      int nPionsInJetTruth = 0;
      int nKaonsInJetTruth = 0;
      int nProtonsInJetTruth = 0;

      int nPionsInJetPosTruth = 0;
      int nKaonsInJetPosTruth = 0;
      int nProtonsInJetPosTruth = 0;

      int nPionsInJetNegTruth = 0;
      int nKaonsInJetNegTruth = 0;
      int nProtonsInJetNegTruth = 0;

      int nSelectedConstituentsTruth = 0;

      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_pt"), djet.pt());
      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_pt_subtracted"), ptSub);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_pt_raw_vs_sub"), djet.pt(), ptSub);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_eta"), djet.eta());
      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_phi"), djet.phi());
      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_area"), djet.area());
      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_n_constituents"), djet.tracksIds().size());

      registryMC.fill(HIST("mc/jets/matched/truth/general/jet_pt"), pjet.pt());
      registryMC.fill(HIST("mc/jets/matched/truth/general/jet_eta"), pjet.eta());
      registryMC.fill(HIST("mc/jets/matched/truth/general/jet_phi"), pjet.phi());
      registryMC.fill(HIST("mc/jets/matched/truth/general/jet_area"), pjet.area());
      registryMC.fill(HIST("mc/jets/matched/truth/general/jet_n_constituents"), pjet.tracksIds().size());

      for (auto const& jtrack : djet.tracks_as<soa::Join<aod::JetTracks, aod::JTrackPIs>>()) {

        auto track = jtrack.track_as<HadronTracksMC>();

        if (!passedTrackSelection(track)) {
          continue;
        }

        ++nSelectedConstituents;

        const double deltaEtaJet = track.eta() - djet.eta();
        const double deltaPhiJet = RecoDecay::constrainAngle(track.phi() - djet.phi(), -PI);
        const double pt = track.pt();
        const double eta = track.eta();
        const double dcaxy = track.dcaXY();
        const double dcaz = track.dcaZ();
        const double phi = track.phi();
        const int charge = track.sign();
        const double p = track.p();

        int realpdg = 0;
        int pdg = 0;
        bool isPrimary = false;

        if (track.has_mcParticle()) {
          auto const& trueParticle = track.mcParticle();
          realpdg = trueParticle.pdgCode();
          pdg = std::abs(realpdg);
          isPrimary = trueParticle.isPhysicalPrimary();
        }

        registryMC.fill(HIST("mc/jets/matched/reconstructed/hadrons/rec_hadron_all"), pt);

        if (isPrimary) {
          registryMC.fill(HIST("mc/jets/matched/reconstructed/hadrons/mc_rec_hadron_pt"), pt);
        } else {
          registryMC.fill(HIST("mc/jets/matched/reconstructed/hadrons/mc_sec_hadron_pt"), pt);
        }

        registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/all_tracks_2d"), deltaEtaJet, deltaPhiJet);
        registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/all_tracks_3d"), deltaEtaJet, deltaPhiJet, pt);

        const PidResult pid = getPid(track);

        if (pid.isPion) {
          registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/pions_2d"), deltaEtaJet, deltaPhiJet);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/pions_3d"), deltaEtaJet, deltaPhiJet, pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pion_pt_vs_jet_pt"), djet.pt(), pt);

          registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/rec_pion_all"), pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/contamination_matrix_pion"), realpdg, pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pion_jet_tpc"), pt, track.tpcNSigmaPi());
          registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pion_tpc_signal_vs_p"), p, track.tpcSignal());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pion_jet_tof"), pt, track.tofNSigmaPi());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pion_tof_beta_vs_p"), p, track.beta());
          }

          if (charge > 0) {
            nPionsInJetPos++;

            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/rec_pion_pos_all"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/contamination_matrix_pion_pos"), realpdg, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_tpc"), pt, track.tpcNSigmaPi());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_phi_vs_pt"), pt, phi);

            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/pions/pion_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/pions/pion_pos_2d"), deltaEtaJet, deltaPhiJet);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_eta"), eta);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_dcaxy"), pt, dcaxy);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_dcaz"), pt, dcaz);

            if (track.hasTOF()) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/rec_pion_pos_tof_matched"), pt);
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/pion_jet_pos_tof"), pt, track.tofNSigmaPi());
            }
          } else {
            nPionsInJetNeg++;

            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/rec_pion_neg_all"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/contamination_matrix_pion_neg"), realpdg, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_tpc"), pt, track.tpcNSigmaPi());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_phi_vs_pt"), pt, phi);

            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/pions/pion_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/pions/pion_neg_2d"), deltaEtaJet, deltaPhiJet);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_eta"), eta);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_dcaxy"), pt, dcaxy);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_dcaz"), pt, dcaz);

            if (track.hasTOF()) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/rec_pion_neg_tof_matched"), pt);
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/pion_jet_neg_tof"), pt, track.tofNSigmaPi());
            }
          }

          if (isPrimary) {
            registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/mc_rec_pion_pt"), pt);
            if (charge > 0) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/mc_rec_pion_pos_pt"), pt);
            } else {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/mc_rec_pion_neg_pt"), pt);
            }
          } else {
            if (pdg == PDG_t::kPiPlus) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/mc_sec_pion_pt"), pt);
              if (charge > 0) {
                registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/mc_sec_pion_pos_pt"), pt);
              } else {
                registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/mc_sec_pion_neg_pt"), pt);
              }
            }
          }
        }

        // --- KAONS ---
        if (pid.isKaon) {
          registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/kaons_2d"), deltaEtaJet, deltaPhiJet);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/kaons_3d"), deltaEtaJet, deltaPhiJet, pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/kaon_pt_vs_jet_pt"), djet.pt(), pt);

          registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/rec_kaon_all"), pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/contamination_matrix_kaon"), realpdg, pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/kaon_jet_tpc"), pt, track.tpcNSigmaKa());
          registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/kaon_tpc_signal_vs_p"), p, track.tpcSignal());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/kaon_jet_tof"), pt, track.tofNSigmaKa());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/kaon_tof_beta_vs_p"), p, track.beta());
          }

          if (charge > 0) {
            nKaonsInJetPos++;

            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/rec_kaon_pos_all"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/contamination_matrix_kaon_pos"), realpdg, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_tpc"), pt, track.tpcNSigmaKa());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_phi_vs_pt"), pt, phi);

            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/kaons/kaon_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/kaons/kaon_pos_2d"), deltaEtaJet, deltaPhiJet);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_eta"), eta);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_dcaxy"), pt, dcaxy);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_dcaz"), pt, dcaz);

            if (track.hasTOF()) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/rec_kaon_pos_tof_matched"), pt);
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/kaon_jet_pos_tof"), pt, track.tofNSigmaKa());
            }
          } else {
            nKaonsInJetNeg++;

            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_all"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/contamination_matrix_kaon_neg"), realpdg, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_tpc"), pt, track.tpcNSigmaKa());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_phi_vs_pt"), pt, phi);

            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/kaons/kaon_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/kaons/kaon_neg_2d"), deltaEtaJet, deltaPhiJet);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_eta"), eta);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_dcaxy"), pt, dcaxy);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/kaon_jet_neg_dcaz"), pt, dcaz);

            if (track.hasTOF()) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_tof_matched"), pt);
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/rec_kaon_neg_tof"), pt, track.tofNSigmaKa());
            }
          }

          if (isPrimary) {
            registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/mc_rec_kaon_pt"), pt);
            if (charge > 0) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/mc_rec_kaon_pos_pt"), pt);
            } else {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/mc_rec_kaon_neg_pt"), pt);
            }
          } else {
            if (pdg == PDG_t::kKPlus) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/mc_sec_kaon_pt"), pt);
              if (charge > 0) {
                registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/mc_sec_kaon_pos_pt"), pt);
              } else {
                registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/mc_sec_kaon_neg_pt"), pt);
              }
            }
          }
        }

        // --- PROTONS ---
        if (pid.isProton) {
          registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/protons_2d"), deltaEtaJet, deltaPhiJet);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/protons_3d"), deltaEtaJet, deltaPhiJet, pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/proton_pt_vs_jet_pt"), djet.pt(), pt);

          registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/rec_proton_all"), pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/contamination_matrix_proton"), realpdg, pt);
          registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/proton_jet_tpc"), pt, track.tpcNSigmaPr());
          registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/proton_tpc_signal_vs_p"), p, track.tpcSignal());
          if (track.hasTOF()) {
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/proton_jet_tof"), pt, track.tofNSigmaPr());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/proton_tof_beta_vs_p"), p, track.beta());
          }

          if (charge > 0) {
            nProtonsInJetPos++;

            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/rec_proton_pos_all"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/contamination_matrix_proton_pos"), realpdg, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_tpc"), pt, track.tpcNSigmaPr());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_phi_vs_pt"), pt, phi);

            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/protons/proton_pos"), deltaEtaJet, deltaPhiJet, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/protons/proton_pos_2d"), deltaEtaJet, deltaPhiJet);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_eta"), eta);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_dcaxy"), pt, dcaxy);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_dcaz"), pt, dcaz);

            if (track.hasTOF()) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/rec_proton_pos_tof_matched"), pt);
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/proton_jet_pos_tof"), pt, track.tofNSigmaPr());
            }
          } else {
            nProtonsInJetNeg++;

            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/rec_proton_neg_all"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/contamination_matrix_proton_neg"), realpdg, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_tpc"), pt, track.tpcNSigmaPr());
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_phi_vs_pt"), pt, phi);

            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/protons/proton_neg"), deltaEtaJet, deltaPhiJet, pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/cone/protons/proton_neg_2d"), deltaEtaJet, deltaPhiJet);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_eta"), eta);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_dcaxy"), pt, dcaxy);
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_dcaz"), pt, dcaz);

            if (track.hasTOF()) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/rec_proton_neg_tof_matched"), pt);
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/proton_jet_neg_tof"), pt, track.tofNSigmaPr());
            }
          }

          if (isPrimary) {
            registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/mc_rec_proton_pt"), pt);
            if (charge > 0) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/mc_rec_proton_pos_pt"), pt);
            } else {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/mc_rec_proton_neg_pt"), pt);
            }
          } else {
            if (pdg == PDG_t::kProton) {
              registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/mc_sec_proton_pt"), pt);
              if (charge > 0) {
                registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/mc_sec_proton_pos_pt"), pt);
              } else {
                registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/mc_sec_proton_neg_pt"), pt);
              }
            }
          }
        }
      }

      registryMC.fill(HIST("mc/jets/matched/reconstructed/general/jet_n_selected_constituents"), nSelectedConstituents);

      nPionsInJet = nPionsInJetPos + nPionsInJetNeg;
      nKaonsInJet = nKaonsInJetPos + nKaonsInJetNeg;
      nProtonsInJet = nProtonsInJetPos + nProtonsInJetNeg;

      registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/pos/n_pions_per_jet"), nPionsInJetPos);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/pos/n_kaons_per_jet"), nKaonsInJetPos);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/pos/n_protons_per_jet"), nProtonsInJetPos);

      registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/neg/n_pions_per_jet"), nPionsInJetNeg);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/neg/n_kaons_per_jet"), nKaonsInJetNeg);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/neg/n_protons_per_jet"), nProtonsInJetNeg);

      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/n_pions_per_jet"), nPionsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/n_kaons_per_jet"), nKaonsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/n_protons_per_jet"), nProtonsInJet);

      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/particle_multiplicity_per_jet"), 1.0, nPionsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/particle_multiplicity_per_jet"), 2.0, nKaonsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/particle_multiplicity_per_jet"), 3.0, nProtonsInJet);

      registryMC.fill(HIST("mc/jets/matched/reconstructed/pions/n_pions_vs_jet_pt"), djet.pt(), nPionsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/kaons/n_kaons_vs_jet_pt"), djet.pt(), nKaonsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/protons/n_protons_vs_jet_pt"), djet.pt(), nProtonsInJet);

      for (auto const& particle : pjet.tracks_as<aod::JetParticles>()) {

        if (!particle.isPhysicalPrimary()) {
          continue;
        }

        if (particle.eta() < minEta || particle.eta() > maxEta) {
          continue;
        }

        double pt = particle.pt();
        if (pt < minPt || pt > maxPt) {
          continue;
        }

        ++nSelectedConstituentsTruth;

        int originalPdg = particle.pdgCode();
        int pdg = std::abs(originalPdg);

        registryMC.fill(HIST("mc/jets/matched/truth/hadrons/mc_gen_hadron_pt"), pt);

        if (pdg == PDG_t::kPiPlus) {
          if (pt >= cfg.minPtPion && pt <= cfg.maxPtPion) {
            registryMC.fill(HIST("mc/jets/matched/truth/pions/mc_gen_pion_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/truth/pions/pion_pt_vs_jet_pt"), pjet.pt(), pt);

            if (originalPdg > 0) {
              nPionsInJetPosTruth++;
              registryMC.fill(HIST("mc/jets/matched/truth/pions/pos/mc_gen_pion_pos_pt"), pt);
            } else {
              nPionsInJetNegTruth++;
              registryMC.fill(HIST("mc/jets/matched/truth/pions/neg/mc_gen_pion_neg_pt"), pt);
            }
          }
        } else if (pdg == PDG_t::kKPlus) {
          if (pt >= cfg.minPtKaon && pt <= cfg.maxPtKaon) {
            registryMC.fill(HIST("mc/jets/matched/truth/kaons/mc_gen_kaon_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/truth/kaons/kaon_pt_vs_jet_pt"), pjet.pt(), pt);

            if (originalPdg > 0) {
              nKaonsInJetPosTruth++;
              registryMC.fill(HIST("mc/jets/matched/truth/kaons/pos/mc_gen_kaon_pos_pt"), pt);
            } else {
              nKaonsInJetNegTruth++;
              registryMC.fill(HIST("mc/jets/matched/truth/kaons/neg/mc_gen_kaon_neg_pt"), pt);
            }
          }
        } else if (pdg == PDG_t::kProton) {
          if (pt >= cfg.minPtProton && pt <= cfg.maxPtProton) {
            registryMC.fill(HIST("mc/jets/matched/truth/protons/mc_gen_proton_pt"), pt);
            registryMC.fill(HIST("mc/jets/matched/truth/protons/proton_pt_vs_jet_pt"), pjet.pt(), pt);
            if (originalPdg > 0) {
              nProtonsInJetPosTruth++;
              registryMC.fill(HIST("mc/jets/matched/truth/protons/pos/mc_gen_proton_pos_pt"), pt);
            } else {
              nProtonsInJetNegTruth++;
              registryMC.fill(HIST("mc/jets/matched/truth/protons/neg/mc_gen_proton_neg_pt"), pt);
            }
          }
        }
      }

      registryMC.fill(HIST("mc/jets/matched/truth/general/jet_n_selected_constituents"), nSelectedConstituentsTruth);

      registryMC.fill(HIST("mc/jets/matched/truth/pions/pos/n_pions_per_jet"), nPionsInJetPosTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/kaons/pos/n_kaons_per_jet"), nKaonsInJetPosTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/protons/pos/n_protons_per_jet"), nProtonsInJetPosTruth);

      registryMC.fill(HIST("mc/jets/matched/truth/pions/neg/n_pions_per_jet"), nPionsInJetNegTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/kaons/neg/n_kaons_per_jet"), nKaonsInJetNegTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/protons/neg/n_protons_per_jet"), nProtonsInJetNegTruth);

      nPionsInJetTruth = nPionsInJetPosTruth + nPionsInJetNegTruth;
      nKaonsInJetTruth = nKaonsInJetPosTruth + nKaonsInJetNegTruth;
      nProtonsInJetTruth = nProtonsInJetPosTruth + nProtonsInJetNegTruth;

      registryMC.fill(HIST("mc/jets/matched/truth/pions/n_pions_per_jet"), nPionsInJetTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/kaons/n_kaons_per_jet"), nKaonsInJetTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/protons/n_protons_per_jet"), nProtonsInJetTruth);

      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/n_pions_matrix"), nPionsInJetTruth, nPionsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/n_kaons_matrix"), nKaonsInJetTruth, nKaonsInJet);
      registryMC.fill(HIST("mc/jets/matched/reconstructed/composition/n_protons_matrix"), nProtonsInJetTruth, nProtonsInJet);

      registryMC.fill(HIST("mc/jets/matched/truth/pions/n_pions_vs_jet_pt"), pjet.pt(), nPionsInJetTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/kaons/n_kaons_vs_jet_pt"), pjet.pt(), nKaonsInJetTruth);
      registryMC.fill(HIST("mc/jets/matched/truth/protons/n_protons_vs_jet_pt"), pjet.pt(), nProtonsInJetTruth);
    }
    registryMC.fill(HIST("mc/jets/matched/reconstructed/general/collision_multiplicity"), jetsInCollision);
  }
  PROCESS_SWITCH(JetHadronsPid, processMatchJets, "Process dector level jets matched to particle level", false);

  void processJetsMCP(MCPJetEvents::iterator const& collision,
                      ChargedMCPJets const& particleLevelJets, aod::JetParticles const&)
  {

    registryMC.fill(HIST("mc/jets/truth/general/n_events_raw"), 1);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    double zVertex = collision.posZ();
    if (std::abs(zVertex) > zVtx) {
      return;
    }

    registryMC.fill(HIST("mc/jets/truth/general/n_events"), 1);
    registryMC.fill(HIST("mc/jets/truth/general/z_vtx"), zVertex);

    double centralRho = collision.rho();
    int jetsInCollision = 0;

    for (auto const& jet : particleLevelJets) {

      if (!isppRefAnalysis && ((std::abs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))) {
        continue;
      }

      if (isppRefAnalysis && std::abs(jet.eta()) > cfgEtaJetMax) {
        continue;
      }

      if (isppRefAnalysis && (jet.pt() < minJetPt || jet.pt() > maxJetPt)) {
        continue;
      }

      double ptSub = jet.pt() - (centralRho * jet.area());
      if (ptSub < 0) {
        ptSub = 0.0;
      }

      if (!isppRefAnalysis && (ptSub < minJetPt || ptSub > maxJetPt)) {
        continue;
      }

      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && normalizedJetArea < minNormalizedJetArea) {
        continue;
      }

      ++jetsInCollision;

      registryMC.fill(HIST("mc/jets/truth/general/jet_pt_subtracted"), ptSub);
      registryMC.fill(HIST("mc/jets/truth/general/jet_pt_raw_vs_sub"), jet.pt(), ptSub);
      registryMC.fill(HIST("mc/jets/truth/general/jet_eta"), jet.eta());
      registryMC.fill(HIST("mc/jets/truth/general/jet_phi"), jet.phi());
      registryMC.fill(HIST("mc/jets/truth/general/jet_area"), jet.area());
      registryMC.fill(HIST("mc/jets/truth/general/jet_pt"), jet.pt());
      registryMC.fill(HIST("mc/jets/truth/general/jet_n_constituents"), jet.tracksIds().size());

      int nPionsInJetPos = 0;
      int nKaonsInJetPos = 0;
      int nProtonsInJetPos = 0;
      int nPionsInJetNeg = 0;
      int nKaonsInJetNeg = 0;
      int nProtonsInJetNeg = 0;
      int selectedConstituentCount = 0;

      for (auto const& particle : jet.tracks_as<aod::JetParticles>()) {

        if (!particle.isPhysicalPrimary()) {
          continue;
        }

        if (particle.eta() < minEta || particle.eta() > maxEta) {
          continue;
        }

        double pt = particle.pt();
        if (pt < minPt || pt > maxPt) {
          continue;
        }

        selectedConstituentCount++;

        int originalPdg = particle.pdgCode();
        int pdg = std::abs(originalPdg);

        registryMC.fill(HIST("mc/jets/truth/hadrons/mc_gen_hadron_pt"), pt);
        if (pdg == PDG_t::kPiPlus) {
          if (pt >= cfg.minPtPion && pt <= cfg.maxPtPion) {
            registryMC.fill(HIST("mc/jets/truth/pions/mc_gen_pion_pt"), pt);
            registryMC.fill(HIST("mc/jets/truth/pions/pion_pt_vs_jet_pt"), jet.pt(), pt);

            if (originalPdg > 0) {
              nPionsInJetPos++;
              registryMC.fill(HIST("mc/jets/truth/pions/pos/mc_gen_pion_pos_pt"), pt);
            } else {
              nPionsInJetNeg++;
              registryMC.fill(HIST("mc/jets/truth/pions/neg/mc_gen_pion_neg_pt"), pt);
            }
          }
        } else if (pdg == PDG_t::kKPlus) {
          if (pt >= cfg.minPtKaon && pt <= cfg.maxPtKaon) {
            registryMC.fill(HIST("mc/jets/truth/kaons/mc_gen_kaon_pt"), pt);
            registryMC.fill(HIST("mc/jets/truth/kaons/kaon_pt_vs_jet_pt"), jet.pt(), pt);

            if (originalPdg > 0) {
              nKaonsInJetPos++;
              registryMC.fill(HIST("mc/jets/truth/kaons/pos/mc_gen_kaon_pos_pt"), pt);
            } else {
              nKaonsInJetNeg++;
              registryMC.fill(HIST("mc/jets/truth/kaons/neg/mc_gen_kaon_neg_pt"), pt);
            }
          }
        } else if (pdg == PDG_t::kProton) {
          if (pt >= cfg.minPtProton && pt <= cfg.maxPtProton) {
            registryMC.fill(HIST("mc/jets/truth/protons/mc_gen_proton_pt"), pt);
            registryMC.fill(HIST("mc/jets/truth/protons/proton_pt_vs_jet_pt"), jet.pt(), pt);
            if (originalPdg > 0) {
              nProtonsInJetPos++;
              registryMC.fill(HIST("mc/jets/truth/protons/pos/mc_gen_proton_pos_pt"), pt);
            } else {
              nProtonsInJetNeg++;
              registryMC.fill(HIST("mc/jets/truth/protons/neg/mc_gen_proton_neg_pt"), pt);
            }
          }
        }
      }

      registryMC.fill(HIST("mc/jets/truth/general/jet_n_selected_constituents"), selectedConstituentCount);

      registryMC.fill(HIST("mc/jets/truth/pions/pos/n_pions_per_jet"), nPionsInJetPos);
      registryMC.fill(HIST("mc/jets/truth/kaons/pos/n_kaons_per_jet"), nKaonsInJetPos);
      registryMC.fill(HIST("mc/jets/truth/protons/pos/n_protons_per_jet"), nProtonsInJetPos);
      registryMC.fill(HIST("mc/jets/truth/pions/neg/n_pions_per_jet"), nPionsInJetNeg);
      registryMC.fill(HIST("mc/jets/truth/kaons/neg/n_kaons_per_jet"), nKaonsInJetNeg);
      registryMC.fill(HIST("mc/jets/truth/protons/neg/n_protons_per_jet"), nProtonsInJetNeg);

      double nPionsInJet = nPionsInJetPos + nPionsInJetNeg;
      double nKaonsInJet = nKaonsInJetPos + nKaonsInJetNeg;
      double nProtonsInJet = nProtonsInJetPos + nProtonsInJetNeg;

      registryMC.fill(HIST("mc/jets/truth/pions/n_pions_per_jet"), nPionsInJet);
      registryMC.fill(HIST("mc/jets/truth/kaons/n_kaons_per_jet"), nKaonsInJet);
      registryMC.fill(HIST("mc/jets/truth/protons/n_protons_per_jet"), nProtonsInJet);

      registryMC.fill(HIST("mc/jets/truth/composition/particle_multiplicity_per_jet"), 1.0, nPionsInJet);
      registryMC.fill(HIST("mc/jets/truth/composition/particle_multiplicity_per_jet"), 2.0, nKaonsInJet);
      registryMC.fill(HIST("mc/jets/truth/composition/particle_multiplicity_per_jet"), 3.0, nProtonsInJet);

      registryMC.fill(HIST("mc/jets/truth/pions/n_pions_vs_jet_pt"), jet.pt(), nPionsInJet);
      registryMC.fill(HIST("mc/jets/truth/kaons/n_kaons_vs_jet_pt"), jet.pt(), nKaonsInJet);
      registryMC.fill(HIST("mc/jets/truth/protons/n_protons_vs_jet_pt"), jet.pt(), nProtonsInJet);
    }
    registryMC.fill(HIST("mc/jets/truth/general/collision_multiplicity"), jetsInCollision);
  }
  PROCESS_SWITCH(JetHadronsPid, processJetsMCP, "Particle-level charged jets in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetHadronsPid>(cfgc)};
}
