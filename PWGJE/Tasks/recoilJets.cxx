// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \author Kotliarov Artem <artem.kotliarov@cern.ch>
/// \file recoilJets.cxx
/// \brief hadron-jet correlation analysis

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include "TRandom3.h"
#include <TH1.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Shorthand notations

// --- Collisions (+ rho)
using CollDataIt = soa::Filtered<aod::JetCollisions>::iterator;
using CollRhoDataIt = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
using CollRhoOutlierDetIt = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos, aod::JCollisionOutliers>>::iterator;
using CollRhoDetIt = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>>::iterator;

using CollPartIt = soa::Filtered<aod::JetMcCollisions>::iterator;
using CollRhoPartTbl = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>;
using CollRhoPartIt = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>>::iterator;
using CollRhoOutlierPartIt = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>>::iterator;
using CollRhoOutlierPartTbl = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>>;

// --- Event multiplicity (+ ZDC etc.)
using EvMultZDCDataIt = soa::Filtered<soa::Join<aod::JetCollisions, aod::ZDCMults>>::iterator;
using EvMultOutlierZDCDetIt = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers, aod::ZDCMults>>::iterator;
using EvMultOutlierPartIt = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::JMcCollisionOutliers>>::iterator;

// --- Tracks / Particles
using TrackTbl = soa::Filtered<aod::JetTracks>;
using PartTbl = soa::Filtered<aod::JetParticles>;

// --- Jets (with constituents)
using JetsDataTbl = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using JetsDetTbl = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>>;
using JetsPartTbl = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>;

// --- Matched jets (det <-> part)
using MatchedJetsDetToPartTbl = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using MatchedJetsPartToDetTbl = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

// --- O2 collisions event selection (not JCollisions)
using CollEvSelExtendedIt = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FT0MultZeqs, aod::MultsExtra, aod::PVMults>>::iterator;
using BCsRun3Tbl = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>; // aod::Run3MatchedToBCExclusive

struct RecoilJets {

  // List of configurable parameters

  // ---------- Event selection ----------
  struct EvCfg : ConfigurableGroup {
    std::string prefix = "event";
    Configurable<std::string> sel{"sel", "sel8", "Choose event selection"},
      triggerMasks{"triggerMasks", "", "Relevant trigger masks: fTrackLowPt,fTrackHighPt"};

    Configurable<float> vertexZCut{"vertexZCut", 10., "Accepted z-vertex range"};
    Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false,
                                       "Flag to choose to reject min. bias gap events; jet-level rejection "
                                       "applied at the jet finder level, here rejection is applied for "
                                       "collision and track process functions"};
  } ev;

  // ---------- RCT / flag-based selections ----------
  struct Rct : ConfigurableGroup {
    std::string prefix = "rct";
    Configurable<std::string> label{"label", "CBT_hadronPID", "Apply rct flag"}; // CBT

    Configurable<bool> enable{"enable", true, "Apply RCT selections"},
      requireZDC{"requireZDC", false, "Require ZDC flag"},
      rejectLimitedAcceptance{"rejectLimitedAcceptance", false, "Reject LimitedAcceptance flag"};
  } rct;

  // ---------- Track selection ----------
  struct Trk : ConfigurableGroup {
    std::string prefix = "track";
    Configurable<std::string> sel{"sel", "globalTracks", "Set track selection"};

    Configurable<float> ptMin{"ptMin", 0.15, "Minimum pT of acceptanced tracks"},
      ptMax{"ptMax", 100., "Maximum pT of acceptanced tracks"},
      etaCut{"etaCut", 0.9, "Eta acceptance of TPC"};
  } trk;

  // ---------- Jet reconstruction ----------
  struct Jet : ConfigurableGroup {
    std::string prefix = "jet";
    Configurable<float> constituentPtMax{"constituentPtMax", 100., "Remove jets with constituent above this pT cut"},
      radius{"radius", 0.4, "Jet cone radius"};
  } jet;

  // ---------- Background tools ----------
  struct Bkgd : ConfigurableGroup {
    std::string prefix = "bkgd";

    // Random cone method
    Configurable<float> randomConeR{"randomConeR", 0.4, "Size of random cone for estimating background fluctuations"},
      minDeltaRToJet{"minDeltaRToJet", 0.0,
                     "Min dR between random cone axis and lead./sublead. jet axis; if 0 -> use R_jet + R_rc"};
  } bkgd;

  // ---------- Normalization FT0 by means ----------
  struct FT0A : ConfigurableGroup {
    std::string prefix = "ft0a";
    Configurable<float> mean{"mean", -1., "Mean FT0A signal"},
      meanPartLevel{"meanPartLevel", -1., "Mean Nch (part level) within FT0A acceptance"},
      meanZeq{"meanZeq", -1., "Mean equalized FT0A signal"};
  } ft0a;

  struct FT0C : ConfigurableGroup {
    std::string prefix = "ft0c";
    Configurable<float> mean{"mean", -1., "Mean FT0C signal"},
      meanPartLevel{"meanPartLevel", -1., "Mean Nch (part level) within FT0C acceptance"},
      meanZeq{"meanZeq", -1., "Mean equalized FT0C signal"};
  } ft0c;

  // ---------- TT / recoil ----------
  struct TT : ConfigurableGroup {
    std::string prefix = "triggerTrack";
    Configurable<float> fracSig{"fracSig", 0.9, "Fraction of events used for signal TT"};
    Configurable<float> recoilRegion{"recoilRegion", 0.6, "Width of recoil acceptance"};

    Configurable<std::vector<float>> refPtRange{"refPtRange", {5., 7}, "Reference TT pT range [min,max] (GeV/c)"},
      sigPtRange{"sigPtRange", {20., 50}, "Signal TT pT range [min,max] (GeV/c)"};

    Configurable<std::vector<float>> phiRestr{"phiRestr", {0., 6.3}, "Phi restriction [min,max] (rad) for TT search"};
  } tt;

  // ---------- Two particle correlations ----------
  struct TwoPartCorrel : ConfigurableGroup {
    std::string prefix = "twoPartCorrel";
    Configurable<std::vector<float>> leadPtRange{"leadPtRange", {4., 6.}, "Leading track pT range [min,max] (GeV/c)"};
    Configurable<float> associatTrackPtMin{"associatTrackPtMin", 2., "Associated track minimum pT (GeV/c)"};
  } twoPartCorrel;

  // ---------- Histogram settings ----------
  struct Hist : ConfigurableGroup {
    std::string prefix = "hist";
    Configurable<bool> sumw2{"sumw2", false, "Enable Sumw2() for histograms"};

    Configurable<uint16_t> jetPtMax{"jetPtMax", 100, "Maximum jet pT stored"},
      multNBins{"multNBins", 600, "Number of bins for scaled FT0M multiplicity"},
      zdcTimeNBins{"zdcTimeNBins", 240, "Number of bins for ZDC timing histograms"};

    ConfigurableAxis multFT0CThresh{"multFT0CThresh", {VARIABLE_WIDTH, 0.0, 0.133, 0.233, 0.367, 0.567, 0.767, 1.067, 1.4, 1.867, 2.5, 3.9, 5.4, 6.9, 20.}, "Percentiles of scaled FT0C: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"}; // default values for raw data
    ConfigurableAxis multFT0CThreshPartLevel{"multFT0CThreshPartLevel", {VARIABLE_WIDTH, 0.0, 0.133, 0.233, 0.367, 0.567, 0.767, 1.067, 1.4, 1.867, 2.5, 3.9, 5.4, 6.9, 20.}, "Percentiles of scaled FT0C: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"};

    ConfigurableAxis multFT0MThresh{"multFT0MThresh", {VARIABLE_WIDTH, 0.0, 0.167, 0.267, 0.4, 0.567, 0.8, 1.067, 1.4, 1.833, 2.433, 3.667, 5.1, 6.433, 20.}, "Percentiles of scaled FT0M: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"}; // default values for raw data
    ConfigurableAxis multFT0MThreshPartLevel{"multFT0MThreshPartLevel", {VARIABLE_WIDTH, 0.0, 0.167, 0.267, 0.4, 0.567, 0.8, 1.067, 1.4, 1.833, 2.433, 3.667, 5.1, 6.433, 20.}, "Percentiles of scaled FT0M: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"};
  } hist;

  // Auxiliary variables
  std::unique_ptr<TRandom3> randGen = std::make_unique<TRandom3>(0);

  // Declare filter on collision Z vertex
  Filter jCollisionFilter = nabs(aod::jcollision::posZ) < ev.vertexZCut.node();
  Filter jCollisionFilterMC = nabs(aod::jmccollision::posZ) < ev.vertexZCut.node();
  Filter collisionFilter = nabs(aod::collision::posZ) < ev.vertexZCut.node();

  // Declare filters on accepted tracks and MC particles (settings for jet reco are provided in the jet finder wagon)
  Filter trackFilter = aod::jtrack::pt > trk.ptMin.node() && aod::jtrack::pt < trk.ptMax.node() && nabs(aod::jtrack::eta) < trk.etaCut.node();
  Filter partFilter = nabs(aod::jmcparticle::eta) < trk.etaCut.node();

  // Declare filter on jets
  Filter jetRadiusFilter = aod::jet::r == nround(jet.radius.node() * 100.);
  Filter jetEtaFilter = nabs(aod::jet::eta) < trk.etaCut.node() - jet.radius.node(); // 0.5 in our analysis

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<MatchedJetsPartToDetTbl> partJetsPerMcCollision = aod::jet::mcCollisionId;

  template <typename AxisObject>
  struct AxisDesc {
    AxisDesc(const char* label_, const AxisObject& axis_, const std::string& axisName_ = "")
      : label(label_), axis(axis_), axisName(axisName_) {}

    const char* label;      // "FT0C" / "FT0M"
    const AxisObject& axis; // AxisSpec or ConfigurableAxis
    std::string axisName;   // Empty for AxisSpec
  };

  void init(InitContext const&)
  {
    // Initialize histogram axes: configurable
    AxisSpec pT{hist.jetPtMax, 0.0, hist.jetPtMax * 1., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPTcorr{hist.jetPtMax + 20, -20., hist.jetPtMax * 1.0, "#it{p}_{T, jet}^{ch, corr} (GeV/#it{c})"};
    AxisSpec scaledFT0A{hist.multNBins, 0.0, 20., "FT0A / #LT FT0A #GT"};
    AxisSpec scaledFT0C{hist.multNBins, 0.0, 20., "FT0C / #LT FT0C #GT"};
    AxisSpec scaledFT0M{hist.multNBins, 0.0, 20., "FT0M^{*}"};
    AxisSpec zdcTiming{hist.zdcTimeNBins, -30., 30., ""};

    // Fixed size histo
    AxisSpec unscaledFT0A{2000, 0.0, 40000., "FT0A"};
    AxisSpec unscaledFT0C{2000, 0.0, 40000., "FT0C"};
    AxisSpec unscaledFT0M{3000, 0.0, 60000., "FT0M (FT0A + FT0C)"};

    AxisSpec zdcNeutronA{1000, 0.0, 5000., "ZNA"};
    AxisSpec zdcNeutronC{1000, 0.0, 5000., "ZNC"};
    AxisSpec zdcNeutronM{4000, 0.0, 8000., "ZNM (ZNA + ZNC)"};

    AxisSpec zdcProtonA{1000, 0.0, 5000., "ZPA"};
    AxisSpec zdcProtonC{1000, 0.0, 5000., "ZPC"};
    AxisSpec zdcProtonM{4000, 0.0, 8000., "ZPM (ZPA + ZPC)"};

    AxisSpec phiAngle{40, 0.0, constants::math::TwoPI, "#it{#varphi} (rad)"};
    AxisSpec deltaPhiAngle{52, 0.0, constants::math::PI, "#Delta#it{#varphi} (rad)"};
    AxisSpec pseudorap{40, -1., 1., "#it{#eta}"};
    AxisSpec pseudorapJets{20, -0.5, 0.5, "#it{#eta}_{jet}"};
    AxisSpec jetArea{50, 0.0, 5., "Area_{jet}"};
    AxisSpec rho{50, 0.0, 50., "#it{#rho}"};

    std::string nameFT0Caxis = "FT0C / #LT FT0C #GT";
    std::string nameFT0Maxis = "FT0M^{*}";

    std::array<AxisDesc<AxisSpec>, 2> arrAxisSpecScaledEA = {{{"FT0C", scaledFT0C},
                                                              {"FT0M", scaledFT0M}}};

    std::array<AxisDesc<AxisSpec>, 3> arrAxisSpecUnscaledEA = {{{"FT0A", unscaledFT0A},
                                                                {"FT0C", unscaledFT0C},
                                                                {"FT0M", unscaledFT0M}}};

    std::array<AxisDesc<ConfigurableAxis>, 2> arrConfigurableAxis = {{{"FT0C", hist.multFT0CThresh, nameFT0Caxis},
                                                                      {"FT0M", hist.multFT0MThresh, nameFT0Maxis}}};

    std::array<AxisDesc<ConfigurableAxis>, 2> arrConfigurableAxisPartLevel = {{{"FT0C", hist.multFT0CThreshPartLevel, nameFT0Caxis},
                                                                               {"FT0M", hist.multFT0MThreshPartLevel, nameFT0Maxis}}};

    // Zero-degree calorimeter
    std::array<AxisDesc<AxisSpec>, 3> arrAxisSpecZDCNeutron = {{{"ZNA", zdcNeutronA},
                                                                {"ZNC", zdcNeutronC},
                                                                {"ZNM", zdcNeutronM}}};

    std::array<AxisDesc<AxisSpec>, 3> arrAxisSpecZDCProton = {{{"ZPA", zdcProtonA},
                                                               {"ZPC", zdcProtonC},
                                                               {"ZPM", zdcProtonM}}};

    // Convert configurable strings to std::string
    std::string evSelToString = static_cast<std::string>(ev.sel);
    std::string trkSelToString = static_cast<std::string>(trk.sel);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(evSelToString);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trkSelToString);
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(ev.triggerMasks);

    const auto phiMin = tt.phiRestr->at(0);
    const auto phiMax = tt.phiRestr->at(1);

    // List of raw and MC det. distributions
    if (doprocessData || doprocessMCDetLevel || doprocessMCDetLevelWeighted) {
      spectra.add("hEventSelectionCount", "Count # of events in the analysis", kTH1F, {{5, 0.0, 5.}});
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of events");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of events after sel. %s", evSelToString.data()));
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(3, "# of events w. outlier");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(4, "# of events w/o assoc MC.");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(5, "# of selected events");

      spectra.add("hTrackSelectionCount", "Count # of tracks in the analysis", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hTrackSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of tracks");
      spectra.get<TH1>(HIST("hTrackSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of tracks after sel. %s", trkSelToString.data()));

      spectra.add("hTTSig_pT", "pT spectrum of all found TT_{Sig} cand.", kTH1F, {{40, 10., 50.}}); // needed to distinguish merged data from diff. wagons

      spectra.add("hJetPtEtaPhiRhoArea", "Charact. of inclusive jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea}, hist.sumw2);
      spectra.add("hJetArea_JetPt_Rho_TTRef", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho}, hist.sumw2);
      spectra.add("hJetArea_JetPt_Rho_TTSig", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho}, hist.sumw2);

      for (const auto& eaAxis : arrConfigurableAxis) {
        spectra.add(Form("hScaled%s_vertexZ", eaAxis.label),
                    "Z vertex of collisions",
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, {60, -12., 12., "#it{z}_{vertex}"}}, hist.sumw2);

        spectra.add(Form("hScaled%sTrackPtEtaPhi", eaAxis.label),
                    "Charact. of tracks",
                    kTHnSparseF, {{eaAxis.axis, eaAxis.axisName}, pT, pseudorap, phiAngle}, hist.sumw2);

        auto tmpHistPointer = spectra.add<TH2>(Form("hScaled%s_Ntrig", eaAxis.label),
                                               Form("Total number of selected triggers per class vs scaled %s", eaAxis.label),
                                               kTH2F, {{eaAxis.axis, eaAxis.axisName}, {2, 0.0, 2.}});
        tmpHistPointer->GetYaxis()->SetBinLabel(1, "TT_{Ref}");
        tmpHistPointer->GetYaxis()->SetBinLabel(2, "TT_{Sig}");

        spectra.add(Form("hScaled%s_TTRef_per_event", eaAxis.label),
                    Form("Number of TT_{Ref} per event vs scaled %s", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, {15, 0.5, 15.5, "# of TT_{Ref}"}});

        spectra.add(Form("hScaled%s_TTSig_per_event", eaAxis.label),
                    Form("Number of TT_{Sig} per event vs scaled %s", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, {10, 0.5, 10.5, "# of TT_{Sig}"}});

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTRef", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTSig", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_TTRef", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_TTSig", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_Corr_TTRef", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_Corr_TTSig", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_TTRef", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_TTSig", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_Rho", eaAxis.label),
                    Form("Scaled %s & #rho", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, rho}, hist.sumw2);

        spectra.add(Form("hScaled%s_Rho_TTRef", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #rho", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, rho}, hist.sumw2);

        spectra.add(Form("hScaled%s_Rho_TTSig", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #rho", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, rho}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTRef_RectrictedPhi", eaAxis.label),
                    Form("Events w. TT_{Ref} #in #varphi (%.2f, %.2f): scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiMin, phiMax, eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTSig_RectrictedPhi", eaAxis.label),
                    Form("Events w. TT_{Sig} #in #varphi (%.2f, %.2f): scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiMin, phiMax, eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);
      }

      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaled%s_TTRef", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s", eaAxis.label),
                    kTH1F, {eaAxis.axis});

        spectra.add(Form("hScaled%s_TTSig", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s", eaAxis.label),
                    kTH1F, {eaAxis.axis});
      }
    }

    // List of MC particle level distributions
    if (doprocessMCPartLevel || doprocessMCPartLevelWeighted) {

      spectra.add("hEventSelectionCountPartLevel", "Count # of events in the part. level analysis", kTH1F, {{4, 0.0, 4.}});
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(1, "Total # of events");
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(2, Form("# of events after sel. %s", evSelToString.data()));
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(3, "# of events w. outlier");
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(4, "# of selected events");

      spectra.add("ptHat", "Distribution of pT hat", kTH1F, {{2000, 0.0, 1000.}});

      spectra.add("hJetPtEtaPhiRhoArea_Part", "Charact. of inclusive part. level jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea}, hist.sumw2);
      spectra.add("hJetArea_JetPt_Rho_TTRef_Part", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho}, hist.sumw2);
      spectra.add("hJetArea_JetPt_Rho_TTSig_Part", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho}, hist.sumw2);

      for (const auto& eaAxis : arrConfigurableAxisPartLevel) {
        spectra.add(Form("hScaled%s_vertexZMC", eaAxis.label),
                    "Z vertex of MC collision",
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, {60, -12., 12., "#it{z}_{vertex}"}}, hist.sumw2);

        auto tmpHistPointer = spectra.add<TH2>(Form("hScaled%s_Ntrig_Part", eaAxis.label),
                                               Form("Total number of selected triggers per class vs scaled %s", eaAxis.label),
                                               kTH2F, {{eaAxis.axis, eaAxis.axisName}, {2, 0.0, 2.}});
        tmpHistPointer->GetYaxis()->SetBinLabel(1, "TT_{Ref}");
        tmpHistPointer->GetYaxis()->SetBinLabel(2, "TT_{Sig}");

        spectra.add(Form("hScaled%sPartPtEtaPhi", eaAxis.label),
                    "Charact. of particles",
                    kTHnSparseF, {{eaAxis.axis, eaAxis.axisName}, pT, pseudorap, phiAngle}, hist.sumw2);

        spectra.add(Form("hScaled%s_TTRef_per_event_Part", eaAxis.label),
                    Form("Number of TT_{Ref} per event vs scaled %s", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, {15, 0.5, 15.5, "# of TT_{Ref}"}});

        spectra.add(Form("hScaled%s_TTSig_per_event_Part", eaAxis.label),
                    Form("Number of TT_{Sig} per event vs scaled %s", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, {10, 0.5, 10.5, "# of TT_{Sig}"}});

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTRef_Part", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTSig_Part", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_TTRef_Part", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_TTSig_Part", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_Corr_TTRef_Part", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_Corr_TTSig_Part", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_TTRef_Part", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_Recoil_JetPt_TTSig_Part", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #it{p}_{T} of recoil jets", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, pT}, hist.sumw2);

        spectra.add(Form("hScaled%s_Rho_Part", eaAxis.label),
                    Form("Scaled %s & #rho", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, rho}, hist.sumw2);

        spectra.add(Form("hScaled%s_Rho_TTRef_Part", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s & #rho", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, rho}, hist.sumw2);

        spectra.add(Form("hScaled%s_Rho_TTSig_Part", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s & #rho", eaAxis.label),
                    kTH2F, {{eaAxis.axis, eaAxis.axisName}, rho}, hist.sumw2);

        // Rectricted phi range for TT selection
        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTRef_RectrictedPhi_Part", eaAxis.label),
                    Form("Events w. TT_{Ref} #in #varphi (%.2f, %.2f): scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiMin, phiMax, eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);

        spectra.add(Form("hScaled%s_DPhi_JetPt_Corr_TTSig_RectrictedPhi_Part", eaAxis.label),
                    Form("Events w. TT_{Sig} #in #varphi (%.2f, %.2f): scaled %s & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiMin, phiMax, eaAxis.label),
                    kTH3F, {{eaAxis.axis, eaAxis.axisName}, deltaPhiAngle, jetPTcorr}, hist.sumw2);
      }

      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaled%s_TTRef_Part", eaAxis.label),
                    Form("Events w. TT_{Ref}: scaled %s", eaAxis.label),
                    kTH1F, {eaAxis.axis});

        spectra.add(Form("hScaled%s_TTSig_Part", eaAxis.label),
                    Form("Events w. TT_{Sig}: scaled %s", eaAxis.label),
                    kTH1F, {eaAxis.axis});
      }
    }

    // Jet matching
    if (doprocessJetsGeoMatching || doprocessJetsGeoMatchingWeighted || doprocessJetsGeoPtMatchingWeighted || doprocessJetsGeoPtMatching) {
      AxisSpec detJetPt{200, 0.0, 200., "#it{p}_{T, det} (GeV/#it{c})"};
      AxisSpec detJetPtCorr{220, -20.0, 200., "#it{p}_{T, det}^{corr.} (GeV/#it{c})"};

      AxisSpec partJetPt{200, 0.0, 200., "#it{p}_{T, part} (GeV/#it{c})"};
      AxisSpec partJetPtCorr{220, -20.0, 200., "#it{p}_{T, part}^{corr.} (GeV/#it{c})"};

      AxisSpec relJetSmearPt{100, -5., 1., "(#it{p}_{T, part} - #it{p}_{T, det}) / #it{p}_{T, part}"};

      //====================================================================================
      // Part. level jets
      spectra.add("hPartLevelInclusiveJetsPt",
                  "All part. level inclusive jets",
                  kTH1F, {partJetPt}, hist.sumw2);

      spectra.add("hPartLevelInclusiveJetsPtCorr",
                  "All part. level inclusive jets",
                  kTH1F, {partJetPtCorr}, hist.sumw2);

      spectra.add("hPartLevelRecoilJetsPt",
                  "All part. level recoil jets",
                  kTH1F, {partJetPt}, hist.sumw2);

      spectra.add("hPartLevelRecoilJetsPtCorr",
                  "All part. level recoil jets",
                  kTH1F, {partJetPtCorr}, hist.sumw2);

      spectra.add("hMissedInclusiveJetsPt",
                  "Part. level inclusive jets w/o matched pair",
                  kTH1F, {partJetPt}, hist.sumw2);

      spectra.add("hMissedInclusiveJetsPtCorr",
                  "Part. level inclusive jets w/o matched pair",
                  kTH1F, {partJetPtCorr}, hist.sumw2);

      spectra.add("hMissedRecoilJetsPt",
                  "Part. level recoil jets w/o matched pair",
                  kTH1F, {partJetPt}, hist.sumw2);

      spectra.add("hMissedRecoilJetsPtCorr",
                  "Part. level recoil jets w/o matched pair",
                  kTH1F, {partJetPtCorr}, hist.sumw2);

      //====================================================================================
      // Det. level jets
      spectra.add("hDetLevelInclusiveJetsPt",
                  "All reconstructed inclusive jets",
                  kTH1F, {detJetPt}, hist.sumw2);

      spectra.add("hDetLevelInclusiveJetsPtCorr",
                  "All reconstructed inclusive jets",
                  kTH1F, {detJetPtCorr}, hist.sumw2);

      spectra.add("hDetLevelRecoilJetsPt",
                  "All reconstructed recoil jets",
                  kTH1F, {detJetPt}, hist.sumw2);

      spectra.add("hDetLevelRecoilJetsPtCorr",
                  "All reconstructed recoil jets",
                  kTH1F, {detJetPtCorr}, hist.sumw2);

      spectra.add("hFakeInclusiveJetsPt",
                  "Det. level inclusive jets w/o matched pair",
                  kTH1F, {detJetPt}, hist.sumw2);

      spectra.add("hFakeInclusiveJetsPtCorr",
                  "Det. level inclusive jets w/o matched pair",
                  kTH1F, {detJetPtCorr}, hist.sumw2);

      spectra.add("hFakeRecoilJetsPt",
                  "Det. level recoil jets w/o matched pair",
                  kTH1F, {detJetPt}, hist.sumw2);

      spectra.add("hFakeRecoilJetsPtCorr",
                  "Det. level recoil jets w/o matched pair",
                  kTH1F, {detJetPtCorr}, hist.sumw2);

      //====================================================================================
      // Response matices with inclusive jets
      spectra.add("hResponseMatrixInclusiveJetsPt",
                  "Correlation inclusive #it{p}_{T, det.} vs. #it{p}_{T, part}",
                  kTH2F, {detJetPt, partJetPt}, hist.sumw2);

      spectra.add("hResponseMatrixInclusiveJetsPtCorr",
                  "Correlation inclusive #it{p}_{T, det.} vs. #it{p}_{T, part}",
                  kTH2F, {detJetPtCorr, partJetPtCorr}, hist.sumw2);

      //====================================================================================
      // Response matices with recoil jets
      spectra.add("hResponseMatrixRecoilJetsPt",
                  "Correlation recoil #it{p}_{T, det.} vs. #it{p}_{T, part}",
                  kTH2F, {detJetPt, partJetPt}, hist.sumw2);

      spectra.add("hResponseMatrixRecoilJetsPtCorr",
                  "Correlation recoil #it{p}_{T, det.} vs. #it{p}_{T, part}",
                  kTH2F, {detJetPtCorr, partJetPtCorr}, hist.sumw2);

      //====================================================================================
      // Jet energy scale and resolution: pT and phi
      spectra.add("hInclusiveJESPt",
                  "ES of inclusive jets vs. #it{p}_{T, part}",
                  kTH2F, {relJetSmearPt, partJetPt}, hist.sumw2);

      spectra.add("hInclusiveJESPtCorr",
                  "ES of inclusive jets vs. #it{p}_{T, part}",
                  kTH2F, {relJetSmearPt, partJetPtCorr}, hist.sumw2);

      spectra.add("hRecoilJESPt",
                  "ES of recoil jets vs. #it{p}_{T, part}",
                  kTH2F, {relJetSmearPt, partJetPt}, hist.sumw2);

      spectra.add("hRecoilJESPtCorr",
                  "ES of recoil jets vs. #it{p}_{T, part}",
                  kTH2F, {relJetSmearPt, partJetPtCorr}, hist.sumw2);

      AxisSpec jetSmearPhi{100, -0.5, 0.5, "#it{#varphi}_{part} - #it{#varphi}_{det}"};
      spectra.add("hInclusiveJESPhi",
                  "#varphi resolution as a func. of jet #it{p}_{T, part}",
                  kTH2F, {jetSmearPhi, partJetPt}, hist.sumw2);

      spectra.add("hRecoilJESPhi",
                  "#varphi resolution as a func. of jet #it{p}_{T, part}",
                  kTH2F, {jetSmearPhi, partJetPt}, hist.sumw2);

      //====================================================================================
      // Utility histograms
      spectra.add("hNumberMatchedInclusiveDetJetsPerOnePartJet",
                  "# of det. level inclusive jets per 1 part. level jet vs. #it{p}_{T, det.} vs. #it{p}_{T, part.}",
                  kTH3F, {{6, 0.5, 6.5, "# of matched det. level jets"}, {200, 0.0, 200., "#it{p}_{T, det.}"}, {200, 0.0, 200., "#it{p}_{T, part.}"}}, hist.sumw2);

      spectra.add("hNumberMatchedRecoilDetJetsPerOnePartJet",
                  "# of det. level recoil jets per 1 part. level jet vs. #it{p}_{T, det.} vs. #it{p}_{T, part.}",
                  kTH3F, {{6, 0.5, 6.5, "# of matched det. level jets"}, {200, 0.0, 200., "#it{p}_{T, det.}"}, {200, 0.0, 200., "#it{p}_{T, part.}"}}, hist.sumw2);
    }

    // Multiplicity for raw data and detector level MC
    if (doprocessEventActivityOO || doprocessEventActivityMCDetLevelWeightedOO) {

      //====================================================================================
      // FIT data
      for (const auto& eaAxis : arrAxisSpecUnscaledEA) {
        spectra.add(Form("hMult%s", eaAxis.label),
                    Form("Mult. signal %s", eaAxis.label),
                    kTH1F, {eaAxis.axis}, hist.sumw2);
      }

      spectra.add("hScaleMultFT0A", "Scaled FTOA signal", kTH1F, {scaledFT0A}, hist.sumw2);
      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaleMult%s", eaAxis.label),
                    Form("Scaled %s signal", eaAxis.label),
                    kTH1F, {eaAxis.axis}, hist.sumw2);
      }

      //====================================================================================
      // Zero-degree calorimeter data
      for (size_t i = 0; i < arrAxisSpecZDCNeutron.size(); ++i) {
        spectra.add(Form("hMult%s", arrAxisSpecZDCNeutron[i].label),
                    Form("Mult. signal from %s", arrAxisSpecZDCNeutron[i].label),
                    kTH1F, {arrAxisSpecZDCNeutron[i].axis}, hist.sumw2);

        spectra.add(Form("hMult%s", arrAxisSpecZDCProton[i].label),
                    Form("Mult. signal from %s", arrAxisSpecZDCProton[i].label),
                    kTH1F, {arrAxisSpecZDCProton[i].axis}, hist.sumw2);

        // Correlation
        spectra.add(Form("h%s_vs_%s", arrAxisSpecZDCProton[i].label, arrAxisSpecZDCNeutron[i].label),
                    Form("Correlation of signals %s vs %s", arrAxisSpecZDCProton[i].label, arrAxisSpecZDCNeutron[i].label),
                    kTH2F, {{arrAxisSpecZDCProton[i].axis}, {arrAxisSpecZDCNeutron[i].axis}}, hist.sumw2);
      }

      //====================================================================================
      // FT0 vs. ZDC correlation
      for (size_t i = 0; i < arrAxisSpecUnscaledEA.size(); ++i) {
        spectra.add(Form("hMult%s_vs_%s", arrAxisSpecUnscaledEA[i].label, arrAxisSpecZDCNeutron[i].label),
                    Form("Correlation of signals %s vs %s", arrAxisSpecUnscaledEA[i].label, arrAxisSpecZDCNeutron[i].label),
                    kTH2F, {{arrAxisSpecUnscaledEA[i].axis}, {arrAxisSpecZDCNeutron[i].axis}}, hist.sumw2);

        if (i != 0) {
          spectra.add(Form("hScaleMult%s_vs_%s", arrAxisSpecScaledEA[i - 1].label, arrAxisSpecZDCNeutron[i].label),
                      Form("Correlation of signals scaled %s vs %s", arrAxisSpecScaledEA[i - 1].label, arrAxisSpecZDCNeutron[i].label),
                      kTH2F, {{arrAxisSpecScaledEA[i - 1].axis}, {arrAxisSpecZDCNeutron[i].axis}}, hist.sumw2);

          spectra.add(Form("hScaleMult%s_vs_%s", arrAxisSpecScaledEA[i - 1].label, arrAxisSpecZDCProton[i].label),
                      Form("Correlation of signals scaled %s vs %s", arrAxisSpecScaledEA[i - 1].label, arrAxisSpecZDCProton[i].label),
                      kTH2F, {{arrAxisSpecScaledEA[i - 1].axis}, {arrAxisSpecZDCProton[i].axis}}, hist.sumw2);
        } else {
          spectra.add("hScaleMultFT0A_vs_ZNA", "Correlation of signals scaled FT0A vs ZNA", kTH2F, {{scaledFT0A}, {arrAxisSpecZDCNeutron[i].axis}}, hist.sumw2);
          spectra.add("hScaleMultFT0A_vs_ZPA", "Correlation of signals scaled FT0A vs ZPA", kTH2F, {{scaledFT0A}, {arrAxisSpecZDCProton[i].axis}}, hist.sumw2);
        }
      }

      spectra.add("hScaleMultFT0M_vs_ZNA_vs_ZNC",
                  "Correlation of signals FT0M^{*} vs ZNA vs ZNC",
                  kTH3F, {{scaledFT0M}, {600, 0.0, 3000., "ZNA"}, {600, 0.0, 3000., "ZNC"}}, hist.sumw2);

      spectra.add("hScaleMultFT0M_vs_ZPA_vs_ZPC",
                  "Correlation of signals FT0M^{*} vs ZPA vs ZPC",
                  kTH3F, {{scaledFT0M}, {600, 0.0, 3000., "ZPA"}, {600, 0.0, 3000., "ZPC"}}, hist.sumw2);
    }

    // Multiplicity for particle level MC
    if (doprocessEventActivityMCPartLevel || doprocessEventActivityMCPartLevelWeighted) {
      spectra.add("hMultFT0APartLevel", "# of primary particles within FTOA acceptance", kTH1F, {{2000, 0.0, 1000., "FT0A"}}, hist.sumw2);
      spectra.add("hMultFT0CPartLevel", "# of primary particles within FTOC acceptance", kTH1F, {{2000, 0.0, 1000., "FT0C"}}, hist.sumw2);
      spectra.add("hMultFT0MPartLevel", "Total # of primary particles from FT0A & FTOC", kTH1F, {{4000, 0.0, 2000., "FT0M"}}, hist.sumw2);

      spectra.add("hScaleMultFT0APartLevel", "Scaled # of primary particles within FTOA acceptance", kTH1F, {{scaledFT0A}}, hist.sumw2);
      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaleMult%sPartLevel", eaAxis.label),
                    Form("Scaled # of primary particles within %s acceptance", eaAxis.label),
                    kTH1F, {{eaAxis.axis}}, hist.sumw2);
      }
    }

    if (doprocessEventActivityQA) {
      spectra.add("hEventSelectionCount", "Count # of events in the analysis", kTH1F, {{5, 0.0, 5.}});
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(1, "sel8");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(2, "IsGoodZvtxFT0vsPV");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(3, "NoSameBunchPileup");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(4, "NoCollInTimeRangeStandard");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(5, "All flags");

      //====================================================================================
      // ZNA vs. ZNC correlation
      spectra.add("hTimeCorrZnaZnc",
                  "Correlat. #it{t}_{ZNA} - #it{t}_{ZNC} vs. #it{t}_{ZNA} + #it{t}_{ZNC}",
                  kTH2F, {{500, -10., 10., "#it{t}_{ZNA} - #it{t}_{ZNC} (ns)"}, {500, -10., 10., "#it{t}_{ZNA} + #it{t}_{ZNC} (ns)"}}, hist.sumw2);
      for (const auto& eaAxis : arrAxisSpecScaledEA) {

        // ZNA vs. ZNC vs. scaled FIT signal
        spectra.add(Form("hTimeZnaVsZncVs%s", eaAxis.label),
                    Form("Correlat. #it{t}_{ZNA} (ns) vs. #it{t}_{ZNC} (ns) vs. scaled %s", eaAxis.label),
                    kTH3F, {{zdcTiming}, {zdcTiming}, {eaAxis.axis}}, hist.sumw2);

        // Number of tracks from PV within acceptance |eta| < 0.8
        spectra.add(Form("hScaled%s_TracksPV", eaAxis.label),
                    Form("Correlat. scaled %s vs. PV tracks", eaAxis.label),
                    kTH2F, {{eaAxis.axis}, {700, 0., 700.}}, hist.sumw2);

        // ITS-only tracks
        spectra.add(Form("hScaled%s_ITStracks", eaAxis.label),
                    Form("Correlat. scaled %s vs. number of ITS tracks", eaAxis.label),
                    kTH2F, {{eaAxis.axis}, {700, 0., 700.}}, hist.sumw2);
      }

      //====================================================================================
      // EA equalized for the vertex position with FT0 detector
      for (const auto& eaAxis : arrAxisSpecUnscaledEA) {
        spectra.add(Form("hMultZeq%s", eaAxis.label),
                    Form("Equalized mult. %s", eaAxis.label),
                    kTH1F, {{eaAxis.axis}}, hist.sumw2);
      }

      // Scaled EA
      spectra.add("hScaledZeqFT0A", "Equalized scaled FT0A", kTH1F, {{scaledFT0A}}, hist.sumw2);
      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaledZeq%s", eaAxis.label),
                    Form("Equalized scaled %s", eaAxis.label),
                    kTH1F, {{eaAxis.axis}}, hist.sumw2);
      }

      //====================================================================================
      // Run-by-run study of EA
      std::vector<const char*> runNumbersOO = {
        "564356", "564359", "564373", "564374", "564387", "564400", "564414", "564430", "564445"};
      const int nRunsOO = runNumbersOO.size();

      std::vector<const char*> evSelFlags = {
        "sel8", "sel8 + IsGoodZvtxFT0vsPV", "sel8 + NoSameBunchPileup", "NoCollInTimeRangeStandard", "sel8 && IsGoodZvtxFT0vsPV && NoSameBunchPileup && NoCollInTimeRangeStandard"};
      const int nEvSelFlags = evSelFlags.size();

      // Scaled FT0 signal; Run-by-run QA
      spectra.add("hScaledFT0APerRunPerSetOfFlags",
                  "Scaled FT0A signal per run per set of ev. sel. flags",
                  kTH3F, {{scaledFT0A}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}}, hist.sumw2);
      setBinLablesYZaxes(spectra.get<TH3>(HIST("hScaledFT0APerRunPerSetOfFlags")), runNumbersOO, evSelFlags);

      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        auto tmpHistPointer = spectra.add<TH3>(Form("hScaled%sPerRunPerSetOfFlags", eaAxis.label),
                                               Form("Scaled %s signal per run per set of ev. sel. flags", eaAxis.label),
                                               kTH3F, {{eaAxis.axis}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}}, hist.sumw2);
        setBinLablesYZaxes(tmpHistPointer, runNumbersOO, evSelFlags);
      }

      // Unscaled FT0 signal; check whether mean value is the same for all runs
      for (const auto& eaAxis : arrAxisSpecUnscaledEA) {
        auto tmpHistPointer = spectra.add<TH3>(Form("h%sPerRunPerSetOfFlags", eaAxis.label),
                                               Form("%s signal per run per set of ev. sel. flags", eaAxis.label),
                                               kTH3F, {{eaAxis.axis}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}}, hist.sumw2);
        setBinLablesYZaxes(tmpHistPointer, runNumbersOO, evSelFlags);
      }

      // Check whether each BC has FT0 signal
      spectra.add("hIsFT0SignalComeFromCollPerRun", "", kTH2F, {{4, 0., 4.}, {nRunsOO, 0., nRunsOO * 1.}});
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(1, "BC has FT0");
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(2, "BC has not FT0");
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(3, "Coll. w. BC");
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(4, "Coll. w/o BC");
      setBinLablesYZaxes(spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun")), runNumbersOO, {});

      // FT0 signal for the case when there is no associated BC
      spectra.add("hScaledFT0AsignalWithoutBC", "", kTH2F, {{scaledFT0A}, {nRunsOO, 0., nRunsOO * 1.}});
      setBinLablesYZaxes(spectra.get<TH2>(HIST("hScaledFT0AsignalWithoutBC")), runNumbersOO, {});

      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        auto tmpHistPointer = spectra.add<TH2>(Form("hScaled%ssignalWithoutBC", eaAxis.label),
                                               "",
                                               kTH2F, {{eaAxis.axis}, {nRunsOO, 0., nRunsOO * 1.}});
        setBinLablesYZaxes(tmpHistPointer, runNumbersOO, {});
      }
    }

    // Di-hadron correlation
    if (doprocessLeadingAndAssociatedTracksTask) {
      const auto pTLeadTrackMin = twoPartCorrel.leadPtRange->at(0);
      const auto pTLeadTrackMax = twoPartCorrel.leadPtRange->at(1);

      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaled%s_NleadTracks", eaAxis.label),
                    Form("Total number of selected leading tracks vs scaled %s", eaAxis.label),
                    kTH2F, {{eaAxis.axis}, {1, 0.0, 1.}}, hist.sumw2);

        spectra.add(Form("hScaled%s_Correlation_LeadTrack_AssociatTracks", eaAxis.label),
                    Form("Leading track #it{p}_{T} #in (%.2f, %.2f); Associated track #it{p}_{T} #in (%.2f, #it{p}_{T, lead. trk})", pTLeadTrackMin, pTLeadTrackMax, twoPartCorrel.associatTrackPtMin.value),
                    kTH2F, {{eaAxis.axis}, {160, -1.28, 5.0, "#it{#varphi} (rad)"}}, hist.sumw2);
      }
    }

    // Bkgd fluctuations in raw and MC det. level data
    const auto ptTTsigMin = tt.sigPtRange->at(0);
    const auto ptTTsigMax = tt.sigPtRange->at(1);

    if (doprocessBkgdFluctuations || doprocessBkgdFluctuationsMCDetLevel || doprocessBkgdFluctuationsMCDetLevelWeighted) {
      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaled%s_deltaPtRandomCone", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f vs. EA", bkgd.randomConeR.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtRandomConeAvoidLeadJet", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f avoid lead jet in vic. %.1f vs. EA", bkgd.randomConeR.value, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtPerpConeAvoidLeadJet", eaAxis.label),
                    Form("Bkgd fluctuations PC with #it{R} = %.1f avoid lead jet in vic. %.1f vs. EA", bkgd.randomConeR.value, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtRandomConeAvoidLeadAndSubleadJet", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f avoid lead, sublead jet in vic. %.1f vs. EA", bkgd.randomConeR.value, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtRandomConeInEventTTSig", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f in events with TT{%.0f, %.0f} in vic. %.1f vs. EA", bkgd.randomConeR.value, ptTTsigMin, ptTTsigMax, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);
      }
    }

    // Bkgd fluctuations in MC part. level data
    if (doprocessBkgdFluctuationsMCPartLevel || doprocessBkgdFluctuationsMCPartLevelWeighted) {
      for (const auto& eaAxis : arrAxisSpecScaledEA) {
        spectra.add(Form("hScaled%s_deltaPtRandomCone_PartLevel", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f vs. EA", bkgd.randomConeR.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtRandomConeAvoidLeadJet_PartLevel", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f avoid lead jet in vic. %.1f vs. EA", bkgd.randomConeR.value, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtPerpConeAvoidLeadJet_PartLevel", eaAxis.label),
                    Form("Bkgd fluctuations PC with #it{R} = %.1f avoid lead jet in vic. %.1f vs. EA", bkgd.randomConeR.value, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtRandomConeAvoidLeadAndSubleadJet_PartLevel", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f avoid lead, sublead jet in vic. %.1f vs. EA", bkgd.randomConeR.value, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);

        spectra.add(Form("hScaled%s_deltaPtRandomConeInEventTTSig_PartLevel", eaAxis.label),
                    Form("Bkgd fluctuations RC with #it{R} = %.1f in events with TT{%.0f, %.0f} in vic. %.1f vs. EA", bkgd.randomConeR.value, ptTTsigMin, ptTTsigMax, bkgd.minDeltaRToJet.value),
                    kTH2F, {{eaAxis.axis}, {400, -40., 60., "#delta#it{p}_{T} (GeV/#it{c})"}}, hist.sumw2);
      }
    }
    spectra.add("hResponseLoopDet", "", kTH2F, {{100, 0., 100., "det. level"}, {100, 0., 100., "part. level"}});
    spectra.add("hResponseLoopPart", "", kTH2F, {{100, 0., 100., "det. level"}, {100, 0., 100., "part. level"}});
  }

  //=============================================================================
  //  Recoil jet analysis
  //=============================================================================

  // Fill histograms with raw or MC det. level data
  template <typename JCollision, typename Jets, typename JTracks>
  void fillHistograms(JCollision const& collision,
                      Jets const& jets,
                      JTracks const& tracks,
                      float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float rho = collision.rho();
    float scaledFT0C = getScaledFT0(collision.multFT0C(), ft0c.mean);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), ft0a.mean), scaledFT0C);

    auto dice = randGen->Rndm();
    if (dice < tt.fracSig)
      bSigEv = true;

    spectra.fill(HIST("hScaledFT0C_vertexZ"), scaledFT0C, collision.posZ(), weight);
    spectra.fill(HIST("hScaledFT0M_vertexZ"), scaledFT0M, collision.posZ(), weight);

    spectra.fill(HIST("hScaledFT0C_Rho"), scaledFT0C, rho, weight);
    spectra.fill(HIST("hScaledFT0M_Rho"), scaledFT0M, rho, weight);

    for (const auto& track : tracks) {
      spectra.fill(HIST("hTrackSelectionCount"), 0.5);

      if (skipTrack(track))
        continue;

      float trackPt = track.pt();
      float trackPhi = track.phi();

      spectra.fill(HIST("hTrackSelectionCount"), 1.5);
      spectra.fill(HIST("hScaledFT0CTrackPtEtaPhi"), scaledFT0C, trackPt, track.eta(), trackPhi, weight);
      spectra.fill(HIST("hScaledFT0MTrackPtEtaPhi"), scaledFT0M, trackPt, track.eta(), trackPhi, weight);

      // Search for TT candidate
      const auto ptTTsigMin = tt.sigPtRange->at(0);
      const auto ptTTsigMax = tt.sigPtRange->at(1);
      if (bSigEv && (trackPt > ptTTsigMin && trackPt < ptTTsigMax)) {
        vPhiOfTT.push_back(trackPhi);
        spectra.fill(HIST("hTTSig_pT"), trackPt, weight);
        ++nTT;
      }

      const auto ptTTrefMin = tt.refPtRange->at(0);
      const auto ptTTrefMax = tt.refPtRange->at(1);
      if (!bSigEv && (trackPt > ptTTrefMin && trackPt < ptTTrefMax)) {
        vPhiOfTT.push_back(trackPhi);
        ++nTT;
      }
    }

    if (nTT > 0) { // at least 1 TT

      phiTT = getPhiTT(vPhiOfTT);

      if (bSigEv) {
        spectra.fill(HIST("hScaledFT0C_Ntrig"), scaledFT0C, 1.5, weight);
        spectra.fill(HIST("hScaledFT0M_Ntrig"), scaledFT0M, 1.5, weight);
        spectra.fill(HIST("hScaledFT0C_TTSig_per_event"), scaledFT0C, nTT, weight);
        spectra.fill(HIST("hScaledFT0M_TTSig_per_event"), scaledFT0M, nTT, weight);

        spectra.fill(HIST("hScaledFT0C_TTSig"), scaledFT0C, weight);
        spectra.fill(HIST("hScaledFT0M_TTSig"), scaledFT0M, weight);

        spectra.fill(HIST("hScaledFT0C_Rho_TTSig"), scaledFT0C, rho, weight);
        spectra.fill(HIST("hScaledFT0M_Rho_TTSig"), scaledFT0M, rho, weight);
      } else {
        spectra.fill(HIST("hScaledFT0C_Ntrig"), scaledFT0C, 0.5, weight);
        spectra.fill(HIST("hScaledFT0M_Ntrig"), scaledFT0M, 0.5, weight);
        spectra.fill(HIST("hScaledFT0C_TTRef_per_event"), scaledFT0C, nTT, weight);
        spectra.fill(HIST("hScaledFT0M_TTRef_per_event"), scaledFT0M, nTT, weight);

        spectra.fill(HIST("hScaledFT0C_TTRef"), scaledFT0C, weight);
        spectra.fill(HIST("hScaledFT0M_TTRef"), scaledFT0M, weight);

        spectra.fill(HIST("hScaledFT0C_Rho_TTRef"), scaledFT0C, rho, weight);
        spectra.fill(HIST("hScaledFT0M_Rho_TTRef"), scaledFT0M, rho, weight);
      }
    }

    for (const auto& jet : jets) {
      // skip jets which have a constituent with pT above specified cut
      if (isJetWithHighPtConstituent<JTracks>(jet))
        continue;

      float jetPt = jet.pt();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea"), jetPt, jet.eta(), jet.phi(), rho, jetArea, weight);

      if (nTT > 0) {
        const auto phiMin = tt.phiRestr->at(0);
        const auto phiMax = tt.phiRestr->at(1);

        auto [dphi, bRecoilJet] = isRecoilJet(jet, phiTT);

        if (bSigEv) {
          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTSig"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTSig"), scaledFT0M, dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig"), jetArea, jetPt, rho, weight);

          if (phiTT > phiMin && phiTT < phiMax) {
            spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig_RectrictedPhi"), scaledFT0C, dphi, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig_RectrictedPhi"), scaledFT0M, dphi, jetPtCorr, weight);
          }

          if (bRecoilJet) {
            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_Corr_TTSig"), scaledFT0C, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_Corr_TTSig"), scaledFT0M, jetPtCorr, weight);

            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_TTSig"), scaledFT0C, jetPt, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_TTSig"), scaledFT0M, jetPt, weight);
          }
        } else {
          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTRef"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTRef"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTRef"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTRef"), scaledFT0M, dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTRef"), jetArea, jetPt, rho, weight);

          if (phiTT > phiMin && phiTT < phiMax) {
            spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTRef_RectrictedPhi"), scaledFT0C, dphi, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTRef_RectrictedPhi"), scaledFT0M, dphi, jetPtCorr, weight);
          }

          if (bRecoilJet) {
            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_Corr_TTRef"), scaledFT0C, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_Corr_TTRef"), scaledFT0M, jetPtCorr, weight);

            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_TTRef"), scaledFT0C, jetPt, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_TTRef"), scaledFT0M, jetPt, weight);
          }
        }
      }
    }
  }

  template <typename JCollision, typename Jets, typename JParticles>
  void fillHistogramsMCPartLevel(JCollision const& collision,
                                 Jets const& jets,
                                 JParticles const& particles,
                                 float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float rho = collision.rho();
    float scaledFT0C = getScaledFT0(collision.multFT0C(), ft0c.meanPartLevel);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), ft0a.meanPartLevel), scaledFT0C);

    auto dice = randGen->Rndm();
    if (dice < tt.fracSig)
      bSigEv = true;

    spectra.fill(HIST("hScaledFT0C_vertexZMC"), scaledFT0C, collision.posZ(), weight);
    spectra.fill(HIST("hScaledFT0M_vertexZMC"), scaledFT0M, collision.posZ(), weight);

    spectra.fill(HIST("hScaledFT0C_Rho_Part"), scaledFT0C, rho, weight);
    spectra.fill(HIST("hScaledFT0M_Rho_Part"), scaledFT0M, rho, weight);

    for (const auto& particle : particles) {
      if (skipParticle(particle))
        continue;

      float particlePt = particle.pt();
      float particlePhi = particle.phi();

      spectra.fill(HIST("hScaledFT0CPartPtEtaPhi"), scaledFT0C, particlePt, particle.eta(), particlePhi, weight);
      spectra.fill(HIST("hScaledFT0MPartPtEtaPhi"), scaledFT0M, particlePt, particle.eta(), particlePhi, weight);

      // Search for TT candidate
      const auto ptTTsigMin = tt.sigPtRange->at(0);
      const auto ptTTsigMax = tt.sigPtRange->at(1);
      if (bSigEv && (particlePt > ptTTsigMin && particlePt < ptTTsigMax)) {
        vPhiOfTT.push_back(particlePhi);
        ++nTT;
      }

      const auto ptTTrefMin = tt.refPtRange->at(0);
      const auto ptTTrefMax = tt.refPtRange->at(1);
      if (!bSigEv && (particlePt > ptTTrefMin && particlePt < ptTTrefMax)) {
        vPhiOfTT.push_back(particlePhi);
        ++nTT;
      }
    }

    if (nTT > 0) {

      phiTT = getPhiTT(vPhiOfTT);

      if (bSigEv) {
        spectra.fill(HIST("hScaledFT0C_Ntrig_Part"), scaledFT0C, 1.5, weight);
        spectra.fill(HIST("hScaledFT0M_Ntrig_Part"), scaledFT0M, 1.5, weight);
        spectra.fill(HIST("hScaledFT0C_TTSig_per_event_Part"), scaledFT0C, nTT, weight);
        spectra.fill(HIST("hScaledFT0M_TTSig_per_event_Part"), scaledFT0M, nTT, weight);

        spectra.fill(HIST("hScaledFT0C_TTSig_Part"), scaledFT0C, weight);
        spectra.fill(HIST("hScaledFT0M_TTSig_Part"), scaledFT0M, weight);

        spectra.fill(HIST("hScaledFT0C_Rho_TTSig_Part"), scaledFT0C, rho, weight);
        spectra.fill(HIST("hScaledFT0M_Rho_TTSig_Part"), scaledFT0M, rho, weight);
      } else {
        spectra.fill(HIST("hScaledFT0C_Ntrig_Part"), scaledFT0C, 0.5, weight);
        spectra.fill(HIST("hScaledFT0M_Ntrig_Part"), scaledFT0M, 0.5, weight);
        spectra.fill(HIST("hScaledFT0C_TTRef_per_event_Part"), scaledFT0C, nTT, weight);
        spectra.fill(HIST("hScaledFT0M_TTRef_per_event_Part"), scaledFT0M, nTT, weight);

        spectra.fill(HIST("hScaledFT0C_TTRef_Part"), scaledFT0C, weight);
        spectra.fill(HIST("hScaledFT0M_TTRef_Part"), scaledFT0M, weight);

        spectra.fill(HIST("hScaledFT0C_Rho_TTRef_Part"), scaledFT0C, rho, weight);
        spectra.fill(HIST("hScaledFT0M_Rho_TTRef_Part"), scaledFT0M, rho, weight);
      }
    }

    for (const auto& jet : jets) {
      float jetPt = jet.pt();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea_Part"), jetPt, jet.eta(), jet.phi(), rho, jetArea, weight);

      if (nTT > 0) {
        const auto phiMin = tt.phiRestr->at(0);
        const auto phiMax = tt.phiRestr->at(1);

        auto [dphi, bRecoilJet] = isRecoilJet(jet, phiTT);

        if (bSigEv) {
          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig_Part"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig_Part"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTSig_Part"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTSig_Part"), scaledFT0M, dphi, jetPt, weight);

          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig_Part"), jetArea, jetPt, rho, weight);

          if (phiTT > phiMin && phiTT < phiMax) {
            spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig_RectrictedPhi_Part"), scaledFT0C, dphi, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig_RectrictedPhi_Part"), scaledFT0M, dphi, jetPtCorr, weight);
          }

          if (bRecoilJet) {
            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_Corr_TTSig_Part"), scaledFT0C, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_Corr_TTSig_Part"), scaledFT0M, jetPtCorr, weight);

            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_TTSig_Part"), scaledFT0C, jetPt, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_TTSig_Part"), scaledFT0M, jetPt, weight);
          }

        } else {

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTRef_Part"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTRef_Part"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTRef_Part"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTRef_Part"), scaledFT0M, dphi, jetPt, weight);

          spectra.fill(HIST("hJetArea_JetPt_Rho_TTRef_Part"), jetArea, jetPt, rho, weight);

          if (phiTT > phiMin && phiTT < phiMax) {
            spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTRef_RectrictedPhi_Part"), scaledFT0C, dphi, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTRef_RectrictedPhi_Part"), scaledFT0M, dphi, jetPtCorr, weight);
          }

          if (bRecoilJet) {
            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_Corr_TTRef_Part"), scaledFT0C, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_Corr_TTRef_Part"), scaledFT0M, jetPtCorr, weight);

            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_TTRef_Part"), scaledFT0C, jetPt, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_TTRef_Part"), scaledFT0M, jetPt, weight);
          }
        }
      }
    }
  }

  //=============================================================================
  // Construction of response matrix
  //=============================================================================
  template <typename JTracks, typename JetsPart, typename JetsDet>
  void fillMatchedGeoHistograms(JetsPart const& jetsPart,
                                JetsDet const& jetsDet,
                                JTracks const& tracks,
                                float partLevelCollRho,
                                float detLevelCollRho,
                                float weight = 1.)
  {
    /// TODO: consider TT cand. on both MC levels simultaneously

    // Utilities for recoil jets; TTsis at det. level MC
    std::vector<double> vPhiOfTT = getPhiOfAllTTsigCandidates(tracks);
    double phiTTSig = 0.0;
    bool bIsThereTTSig = vPhiOfTT.size() > 0;

    if (bIsThereTTSig)
      phiTTSig = getPhiTT(vPhiOfTT);

    for (const auto& jetPart : jetsPart) {
      auto partJetPt = jetPart.pt();
      auto partJetPtCorr = partJetPt - partLevelCollRho * jetPart.area();

      bool bIsPartJetRecoil =
        get<1>(isRecoilJet(jetPart, phiTTSig)) && bIsThereTTSig;

      // Distribution of all part. level jets
      spectra.fill(HIST("hPartLevelInclusiveJetsPt"), partJetPt, weight);
      spectra.fill(HIST("hPartLevelInclusiveJetsPtCorr"), partJetPtCorr, weight);

      if (bIsPartJetRecoil) {
        spectra.fill(HIST("hPartLevelRecoilJetsPt"), partJetPt, weight);
        spectra.fill(HIST("hPartLevelRecoilJetsPtCorr"), partJetPtCorr, weight);
      }

      if (jetPart.has_matchedJetGeo()) {
        const auto& jetsDetMatched = jetPart.template matchedJetGeo_as<JetsDet>();

        for (const auto& jetDetMatched : jetsDetMatched) {
          // Skip matches where detector level jets have a constituent with pT above specified cut
          const bool skipMatchedDetJet = isJetWithHighPtConstituent<JTracks>(jetDetMatched);

          if (skipMatchedDetJet) {
            // Miss jets
            spectra.fill(HIST("hMissedInclusiveJetsPt"), partJetPt, weight);
            spectra.fill(HIST("hMissedInclusiveJetsPtCorr"), partJetPtCorr, weight);

            if (bIsPartJetRecoil) {
              spectra.fill(HIST("hMissedRecoilJetsPt"), partJetPt, weight);
              spectra.fill(HIST("hMissedRecoilJetsPtCorr"), partJetPtCorr, weight);
            }

          } else {
            auto detJetPt = jetDetMatched.pt();
            auto detJetPtCorr = detJetPt - detLevelCollRho * jetDetMatched.area();

            spectra.fill(HIST("hNumberMatchedInclusiveDetJetsPerOnePartJet"), jetsDetMatched.size(), detJetPt, partJetPt, weight);

            spectra.fill(HIST("hResponseMatrixInclusiveJetsPt"), detJetPt, partJetPt, weight);
            spectra.fill(HIST("hResponseMatrixInclusiveJetsPtCorr"), detJetPtCorr, partJetPtCorr, weight);

            spectra.fill(HIST("hInclusiveJESPt"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
            spectra.fill(HIST("hInclusiveJESPtCorr"), (partJetPtCorr - detJetPtCorr) / partJetPtCorr, partJetPtCorr, weight);

            spectra.fill(HIST("hInclusiveJESPhi"), jetPart.phi() - jetDetMatched.phi(), partJetPt, weight);

            if (bIsPartJetRecoil) {
              spectra.fill(HIST("hNumberMatchedRecoilDetJetsPerOnePartJet"), jetsDetMatched.size(), detJetPt, partJetPt, weight);

              spectra.fill(HIST("hResponseMatrixRecoilJetsPt"), detJetPt, partJetPt, weight);
              spectra.fill(HIST("hResponseMatrixRecoilJetsPtCorr"), detJetPtCorr, partJetPtCorr, weight);

              spectra.fill(HIST("hRecoilJESPt"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
              spectra.fill(HIST("hRecoilJESPtCorr"), (partJetPtCorr - detJetPtCorr) / partJetPtCorr, partJetPtCorr, weight);

              spectra.fill(HIST("hRecoilJESPhi"), jetPart.phi() - jetDetMatched.phi(), partJetPt, weight);
            }
          }
        }
      } else {
        // Miss jets
        spectra.fill(HIST("hMissedInclusiveJetsPt"), partJetPt, weight);
        spectra.fill(HIST("hMissedInclusiveJetsPtCorr"), partJetPtCorr, weight);

        if (bIsPartJetRecoil) {
          spectra.fill(HIST("hMissedRecoilJetsPt"), partJetPt, weight);
          spectra.fill(HIST("hMissedRecoilJetsPtCorr"), partJetPtCorr, weight);
        }
      }
    }

    // Reconstructed jets
    for (const auto& jetDet : jetsDet) {
      if (isJetWithHighPtConstituent<JTracks>(jetDet))
        continue;

      auto detJetPt = jetDet.pt();
      auto detJetPtCorr = detJetPt - detLevelCollRho * jetDet.area();

      bool bIsJetRecoil =
        get<1>(isRecoilJet(jetDet, phiTTSig)) && bIsThereTTSig;

      // Distribution of all det. level jets
      spectra.fill(HIST("hDetLevelInclusiveJetsPt"), detJetPt, weight);
      spectra.fill(HIST("hDetLevelInclusiveJetsPtCorr"), detJetPtCorr, weight);

      if (bIsJetRecoil) {
        spectra.fill(HIST("hDetLevelRecoilJetsPt"), detJetPt, weight);
        spectra.fill(HIST("hDetLevelRecoilJetsPtCorr"), detJetPtCorr, weight);
      }

      if (!jetDet.has_matchedJetGeo()) {
        spectra.fill(HIST("hFakeInclusiveJetsPt"), detJetPt, weight);
        spectra.fill(HIST("hFakeInclusiveJetsPtCorr"), detJetPtCorr, weight);

        if (bIsJetRecoil) {
          spectra.fill(HIST("hFakeRecoilJetsPt"), detJetPt, weight);
          spectra.fill(HIST("hFakeRecoilJetsPtCorr"), detJetPtCorr, weight);
        }
      }
    }
  }

  /// TODO: think how to implement matching matching geo/geo+pT in one function
  template <typename JTracks, typename JetsPart, typename JetsDet>
  void fillMatchedGeoPtHistograms(JetsPart const& jetsPart,
                                  JetsDet const& jetsDet,
                                  JTracks const& tracks,
                                  float partLevelCollRho,
                                  float detLevelCollRho,
                                  float weight = 1.)
  {
    // Utilities for recoil jets; TTsis at det. level MC
    std::vector<double> vPhiOfTT = getPhiOfAllTTsigCandidates(tracks);
    double phiTTSig = 0.0;
    bool bIsThereTTSig = vPhiOfTT.size() > 0;

    if (bIsThereTTSig)
      phiTTSig = getPhiTT(vPhiOfTT);

    for (const auto& jetPart : jetsPart) {
      auto partJetPt = jetPart.pt();
      auto partJetPtCorr = partJetPt - partLevelCollRho * jetPart.area();

      bool bIsPartJetRecoil =
        get<1>(isRecoilJet(jetPart, phiTTSig)) && bIsThereTTSig;

      // Distribution of all part. level jets
      spectra.fill(HIST("hPartLevelInclusiveJetsPt"), partJetPt, weight);
      spectra.fill(HIST("hPartLevelInclusiveJetsPtCorr"), partJetPtCorr, weight);

      if (bIsPartJetRecoil) {
        spectra.fill(HIST("hPartLevelRecoilJetsPt"), partJetPt, weight);
        spectra.fill(HIST("hPartLevelRecoilJetsPtCorr"), partJetPtCorr, weight);
      }

      // Geo + pT matching
      if (jetPart.has_matchedJetGeo() && jetPart.has_matchedJetPt()) {
        const auto& jetsDetMatched = jetPart.template matchedJetGeo_as<JetsDet>();
        auto both = intersectMatchIds(jetPart.matchedJetGeoIds(), jetPart.matchedJetPtIds());

        for (const auto& jetDetMatched : jetsDetMatched) {
          // Skip matches where detector level jets have a constituent with pT above specified cut
          const bool skipMatchedDetJet = isJetWithHighPtConstituent<JTracks>(jetDetMatched);

          // If indecies for geo and pT matching do not coinside, reject detector level jet
          const bool skipNotBothId = !both.contains(jetDetMatched.globalIndex());

          if (skipMatchedDetJet || skipNotBothId) {
            // Miss jets
            spectra.fill(HIST("hMissedInclusiveJetsPt"), partJetPt, weight);
            spectra.fill(HIST("hMissedInclusiveJetsPtCorr"), partJetPtCorr, weight);

            if (bIsPartJetRecoil) {
              spectra.fill(HIST("hMissedRecoilJetsPt"), partJetPt, weight);
              spectra.fill(HIST("hMissedRecoilJetsPtCorr"), partJetPtCorr, weight);
            }
          } else {
            auto detJetPt = jetDetMatched.pt();
            auto detJetPtCorr = detJetPt - detLevelCollRho * jetDetMatched.area();

            spectra.fill(HIST("hNumberMatchedInclusiveDetJetsPerOnePartJet"), jetsDetMatched.size(), detJetPt, partJetPt, weight);

            spectra.fill(HIST("hResponseMatrixInclusiveJetsPt"), detJetPt, partJetPt, weight);
            spectra.fill(HIST("hResponseMatrixInclusiveJetsPtCorr"), detJetPtCorr, partJetPtCorr, weight);

            spectra.fill(HIST("hInclusiveJESPt"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
            spectra.fill(HIST("hInclusiveJESPtCorr"), (partJetPtCorr - detJetPtCorr) / partJetPtCorr, partJetPtCorr, weight);

            spectra.fill(HIST("hInclusiveJESPhi"), jetPart.phi() - jetDetMatched.phi(), partJetPt, weight);

            if (bIsPartJetRecoil) {
              spectra.fill(HIST("hNumberMatchedRecoilDetJetsPerOnePartJet"), jetsDetMatched.size(), detJetPt, partJetPt, weight);

              spectra.fill(HIST("hResponseMatrixRecoilJetsPt"), detJetPt, partJetPt, weight);
              spectra.fill(HIST("hResponseMatrixRecoilJetsPtCorr"), detJetPtCorr, partJetPtCorr, weight);

              spectra.fill(HIST("hRecoilJESPt"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
              spectra.fill(HIST("hRecoilJESPtCorr"), (partJetPtCorr - detJetPtCorr) / partJetPtCorr, partJetPtCorr, weight);

              spectra.fill(HIST("hRecoilJESPhi"), jetPart.phi() - jetDetMatched.phi(), partJetPt, weight);
            }
          }
        }
      } else {
        // Miss jets
        spectra.fill(HIST("hMissedInclusiveJetsPt"), partJetPt, weight);
        spectra.fill(HIST("hMissedInclusiveJetsPtCorr"), partJetPtCorr, weight);

        if (bIsPartJetRecoil) {
          spectra.fill(HIST("hMissedRecoilJetsPt"), partJetPt, weight);
          spectra.fill(HIST("hMissedRecoilJetsPtCorr"), partJetPtCorr, weight);
        }
      }
    }

    // Reconstructed jets
    for (const auto& jetDet : jetsDet) {
      if (isJetWithHighPtConstituent<JTracks>(jetDet))
        continue;

      auto detJetPt = jetDet.pt();
      auto detJetPtCorr = detJetPt - detLevelCollRho * jetDet.area();

      bool bIsJetRecoil =
        get<1>(isRecoilJet(jetDet, phiTTSig)) && bIsThereTTSig;

      // Distribution of all det. level jets
      spectra.fill(HIST("hDetLevelInclusiveJetsPt"), detJetPt, weight);
      spectra.fill(HIST("hDetLevelInclusiveJetsPtCorr"), detJetPtCorr, weight);

      if (bIsJetRecoil) {
        spectra.fill(HIST("hDetLevelRecoilJetsPt"), detJetPt, weight);
        spectra.fill(HIST("hDetLevelRecoilJetsPtCorr"), detJetPtCorr, weight);
      }

      if (!jetDet.has_matchedJetGeo()) {
        spectra.fill(HIST("hFakeInclusiveJetsPt"), detJetPt, weight);
        spectra.fill(HIST("hFakeInclusiveJetsPtCorr"), detJetPtCorr, weight);

        if (bIsJetRecoil) {
          spectra.fill(HIST("hFakeRecoilJetsPt"), detJetPt, weight);
          spectra.fill(HIST("hFakeRecoilJetsPtCorr"), detJetPtCorr, weight);
        }
      }
    }
  }

  //=============================================================================
  // Event Activity analysis
  //=============================================================================
  template <typename JCollision>
  void fillMultiplicityHistogramsOO(JCollision const& collision,
                                    float weight = 1.)
  {
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multFT0M = collision.multFT0M();
    float scaledFT0A = getScaledFT0(multFT0A, ft0a.mean);
    float scaledFT0C = getScaledFT0(multFT0C, ft0c.mean);
    float scaledFT0M = getScaledFT0M(scaledFT0A, scaledFT0C);

    float multZNA = collision.multZNA();
    float multZNC = collision.multZNC();
    float multZNM = multZNA + multZNC;

    float multZPA = collision.multZPA();
    float multZPC = collision.multZPC();
    float multZPM = multZPA + multZPC;

    // Individual distributions
    spectra.fill(HIST("hMultFT0A"), multFT0A, weight);
    spectra.fill(HIST("hMultFT0C"), multFT0C, weight);
    spectra.fill(HIST("hMultFT0M"), multFT0M, weight);

    spectra.fill(HIST("hScaleMultFT0A"), scaledFT0A, weight);
    spectra.fill(HIST("hScaleMultFT0C"), scaledFT0C, weight);
    spectra.fill(HIST("hScaleMultFT0M"), scaledFT0M, weight);

    spectra.fill(HIST("hMultZNA"), multZNA, weight);
    spectra.fill(HIST("hMultZNC"), multZNC, weight);
    spectra.fill(HIST("hMultZNM"), multZNM, weight);

    spectra.fill(HIST("hMultZPA"), multZPA, weight);
    spectra.fill(HIST("hMultZPC"), multZPC, weight);
    spectra.fill(HIST("hMultZPM"), multZPM, weight);

    // Correlations
    spectra.fill(HIST("hZPA_vs_ZNA"), multZPA, multZNA, weight);
    spectra.fill(HIST("hZPC_vs_ZNC"), multZPC, multZNC, weight);

    spectra.fill(HIST("hMultFT0A_vs_ZNA"), multFT0A, multZNA, weight);
    spectra.fill(HIST("hMultFT0C_vs_ZNC"), multFT0C, multZNC, weight);
    spectra.fill(HIST("hMultFT0M_vs_ZNM"), multFT0M, multZNM, weight);

    spectra.fill(HIST("hScaleMultFT0A_vs_ZNA"), scaledFT0A, multZNA, weight);
    spectra.fill(HIST("hScaleMultFT0C_vs_ZNC"), scaledFT0C, multZNC, weight);
    spectra.fill(HIST("hScaleMultFT0M_vs_ZNM"), scaledFT0M, multZNM, weight);

    spectra.fill(HIST("hScaleMultFT0A_vs_ZPA"), scaledFT0A, multZPA, weight);
    spectra.fill(HIST("hScaleMultFT0C_vs_ZPC"), scaledFT0C, multZPC, weight);
    spectra.fill(HIST("hScaleMultFT0M_vs_ZPM"), scaledFT0M, multZPM, weight);

    spectra.fill(HIST("hScaleMultFT0M_vs_ZNA_vs_ZNC"), scaledFT0M, multZNA, multZNC, weight);
    spectra.fill(HIST("hScaleMultFT0M_vs_ZPA_vs_ZPC"), scaledFT0M, multZPA, multZPC, weight);
  }

  template <typename JCollisionMC>
  void fillMultiplicityHistogramsMCPartLevel(JCollisionMC const& collision,
                                             float weight = 1.)
  {
    spectra.fill(HIST("hMultFT0APartLevel"), collision.multFT0A(), weight);
    spectra.fill(HIST("hMultFT0CPartLevel"), collision.multFT0C(), weight);
    spectra.fill(HIST("hMultFT0MPartLevel"), collision.multFT0A() + collision.multFT0C(), weight);

    auto scaledFT0A = getScaledFT0(collision.multFT0A(), ft0a.meanPartLevel);
    auto scaledFT0C = getScaledFT0(collision.multFT0C(), ft0c.meanPartLevel);
    spectra.fill(HIST("hScaleMultFT0APartLevel"), scaledFT0A, weight);
    spectra.fill(HIST("hScaleMultFT0CPartLevel"), scaledFT0C, weight);
    spectra.fill(HIST("hScaleMultFT0MPartLevel"), getScaledFT0M(scaledFT0A, scaledFT0C), weight);
  }

  //=============================================================================
  // Event Activity QA analysis in OO collisions (raw and MC detector level (no weight; MB events))
  //=============================================================================
  template <typename BC, typename Collision, typename ZDC>
  void fillMultiplicityQA(Collision const& collision,
                          BC const&,
                          ZDC const&,
                          float weight = 1.)
  {
    int runNumber = collision.multRunNumber();
    int fillNumber = getBinNumberOnYaxisForGivenRun(spectra.get<TH3>(HIST("hScaledFT0CPerRunPerSetOfFlags")), runNumber) - 0.5; // Same for FT0M distrib.

    // FT0 Signal
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multFT0M = collision.multFT0M();
    float scaledFT0A = getScaledFT0(multFT0A, ft0a.mean);
    float scaledFT0C = getScaledFT0(multFT0C, ft0c.mean);
    float scaledFT0M = getScaledFT0M(scaledFT0A, scaledFT0C);

    // Event with flag Sel8
    spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 0.5, weight);
    spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 0.5, weight);
    spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 0.5, weight);

    spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 0.5, weight);
    spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 0.5, weight);
    spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 0.5, weight);

    spectra.fill(HIST("hEventSelectionCount"), 0.5);

    bool isGoodZvtxFT0vsPV = collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV);
    if (isGoodZvtxFT0vsPV) {
      spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 1.5, weight);
      spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 1.5, weight);
      spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 1.5, weight);

      spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 1.5, weight);
      spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 1.5, weight);
      spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 1.5, weight);

      spectra.fill(HIST("hEventSelectionCount"), 1.5);
    }

    bool isNoSameBunchPileup = collision.selection_bit(aod::evsel::kNoSameBunchPileup);
    if (isNoSameBunchPileup) {
      spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 2.5, weight);
      spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 2.5, weight);
      spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 2.5, weight);

      spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 2.5, weight);
      spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 2.5, weight);
      spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 2.5, weight);

      spectra.fill(HIST("hEventSelectionCount"), 2.5);
    }

    bool isNoCollInTimeRangeStandard = collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard);
    if (isNoCollInTimeRangeStandard) {

      spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 3.5, weight);
      spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 3.5, weight);
      spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 3.5, weight);

      spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 3.5, weight);
      spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 3.5, weight);
      spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 3.5, weight);

      spectra.fill(HIST("hEventSelectionCount"), 3.5);
    }

    if (!(isGoodZvtxFT0vsPV && isNoSameBunchPileup & isNoCollInTimeRangeStandard))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 4.5); // All accepted events after 4 flags cut

    spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 4.5, weight);
    spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 4.5, weight);
    spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 4.5, weight);

    spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 4.5, weight);
    spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 4.5, weight);
    spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 4.5, weight);

    //____________________________________________________________________________________
    // Investigate other EA variables

    // Multiplicity equalized for the vertex position with FT0 detector
    float multZeqFT0A = collision.multZeqFT0A();
    float multZeqFT0C = collision.multZeqFT0C();
    float multZeqFT0M = multZeqFT0A + multZeqFT0C;
    float scaledZeqFT0A = getScaledFT0(multZeqFT0A, ft0a.meanZeq);
    float scaledZeqFT0C = getScaledFT0(multZeqFT0C, ft0c.meanZeq);
    float scaledZeqFT0M = getScaledFT0M(scaledZeqFT0A, scaledZeqFT0C);

    spectra.fill(HIST("hMultZeqFT0A"), multZeqFT0A, weight);
    spectra.fill(HIST("hMultZeqFT0C"), multZeqFT0C, weight);
    spectra.fill(HIST("hMultZeqFT0M"), multZeqFT0M, weight);
    spectra.fill(HIST("hScaledZeqFT0A"), scaledZeqFT0A, weight);
    spectra.fill(HIST("hScaledZeqFT0C"), scaledZeqFT0C, weight);
    spectra.fill(HIST("hScaledZeqFT0M"), scaledZeqFT0M, weight);

    // ZDC timing info
    auto const& foundBC = collision.template foundBC_as<BC>();
    float timeZNA = foundBC.has_zdc() ? foundBC.zdc().timeZNA() : -999.f;
    float timeZNC = foundBC.has_zdc() ? foundBC.zdc().timeZNC() : -999.f;
    float timeDiffZDC = timeZNA - timeZNC;
    float timeSumZDC = timeZNA + timeZNC;

    spectra.fill(HIST("hTimeCorrZnaZnc"), timeDiffZDC, timeSumZDC, weight);
    spectra.fill(HIST("hTimeZnaVsZncVsFT0C"), timeZNA, timeZNC, scaledFT0C, weight);
    spectra.fill(HIST("hTimeZnaVsZncVsFT0M"), timeZNA, timeZNC, scaledFT0M, weight);

    // ITS only tracks
    int nITSonly = collision.multNTracksITSOnly();

    spectra.fill(HIST("hScaledFT0C_ITStracks"), scaledFT0C, nITSonly, weight);
    spectra.fill(HIST("hScaledFT0M_ITStracks"), scaledFT0M, nITSonly, weight);

    // Global tracks from PV within |eta| < 0.8
    int multNContribs = collision.multNTracksPV();

    spectra.fill(HIST("hScaledFT0C_TracksPV"), scaledFT0C, multNContribs, weight);
    spectra.fill(HIST("hScaledFT0M_TracksPV"), scaledFT0M, multNContribs, weight);

    if (foundBC.foundFT0Id() > 0) // -1 if does not
    {
      spectra.fill(HIST("hIsFT0SignalComeFromCollPerRun"), 0.5, fillNumber, weight);
    } else {
      spectra.fill(HIST("hIsFT0SignalComeFromCollPerRun"), 1.5, fillNumber, weight);
      spectra.fill(HIST("hScaledFT0AsignalWithoutBC"), scaledFT0A, fillNumber, weight);
      spectra.fill(HIST("hScaledFT0CsignalWithoutBC"), scaledFT0C, fillNumber, weight);
      spectra.fill(HIST("hScaledFT0MsignalWithoutBC"), scaledFT0M, fillNumber, weight);
    }

    if (collision.foundBCId() > 0)
      spectra.fill(HIST("hIsFT0SignalComeFromCollPerRun"), 2.5, fillNumber, weight);
    else
      spectra.fill(HIST("hIsFT0SignalComeFromCollPerRun"), 3.5, fillNumber, weight);
  }

  //=============================================================================
  // Di-hadron azimuthal correlation in raw and MC det. level (no weight; MB events) data
  //=============================================================================
  template <typename JCollision, typename JTracks>
  void fillLeadingAndAssociatedTracksTask(JCollision const& collision, JTracks const& tracks, float weight = 1.)
  {
    std::vector<double> vPhiOfLeadingTracks;
    std::vector<double> vPtOfLeadingTracks;
    std::vector<double> vPhiOfAssociatedTracks;

    float scaledFT0C = getScaledFT0(collision.multFT0C(), ft0c.mean);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), ft0a.mean), scaledFT0C);

    // Search for leading tracks
    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      float trackPt = track.pt();

      if (trackPt > twoPartCorrel.leadPtRange->at(0) && trackPt < twoPartCorrel.leadPtRange->at(1)) {
        vPhiOfLeadingTracks.push_back(track.phi());
        vPtOfLeadingTracks.push_back(trackPt);
      }
    }

    int nLeadingTracks = vPhiOfLeadingTracks.size();

    if (nLeadingTracks > 0) {
      auto indexLeadTrack = randGen->Integer(nLeadingTracks);

      double phiLeadingTrack = vPhiOfLeadingTracks[indexLeadTrack];
      double pTLeadingTrack = vPtOfLeadingTracks[indexLeadTrack];

      spectra.fill(HIST("hScaledFT0C_NleadTracks"), scaledFT0C, 0.5, weight);
      spectra.fill(HIST("hScaledFT0M_NleadTracks"), scaledFT0M, 0.5, weight);

      for (const auto& track : tracks) {
        if (skipTrack(track))
          continue;

        float trackPt = track.pt();
        float trackPhi = track.phi();

        // Search for associated tracks
        if (trackPt > twoPartCorrel.associatTrackPtMin && trackPt < pTLeadingTrack) {
          double dphi = RecoDecay::constrainAngle(phiLeadingTrack - trackPhi, -1.3);
          spectra.fill(HIST("hScaledFT0C_Correlation_LeadTrack_AssociatTracks"), scaledFT0C, dphi, weight);
          spectra.fill(HIST("hScaledFT0M_Correlation_LeadTrack_AssociatTracks"), scaledFT0M, dphi, weight);
        }
      }
    }
  }

  //=============================================================================
  // Estimation of bkgd fluctuations
  //=============================================================================

  // Background fluctuations in raw data and MC det. level
  template <typename JCollision, typename JTracks, typename Jets>
  void fillBkgdFluctuations(JCollision const& collision,
                            Jets const& jets,
                            JTracks const& tracks,
                            float weight = 1.)
  {
    //----------------------------------------------------------
    float rho = collision.rho();
    float scaledFT0C = getScaledFT0(collision.multFT0C(), ft0c.mean);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), ft0a.mean), scaledFT0C);

    //----------------------------------------------------------
    // Study bkgd fluctuations in events with TTsig
    std::vector<uint64_t> vCandForTT;

    //----------------------------------------------------------
    // Place absolutely random cone (anywhere in an event)
    float randomConeEta = randGen->Uniform(-trk.etaCut + bkgd.randomConeR, trk.etaCut - bkgd.randomConeR);
    float randomConePhi = randGen->Uniform(0.0, constants::math::TwoPI);
    float radiusRC2 = std::pow(bkgd.randomConeR, 2);
    float areaRC = constants::math::PI * radiusRC2;
    float randomConePt = 0.0;

    uint64_t index = 0;
    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      float dEta = std::pow(randomConeEta - track.eta(), 2);
      float dPhi = std::pow(RecoDecay::constrainAngle(randomConePhi - track.phi(), -constants::math::PI), 2);

      if ((dEta + dPhi) < radiusRC2) // inside RC
      {
        randomConePt += track.pt();
      }

      // Search for TT_Sig candidate
      const auto ptTTsigMin = tt.sigPtRange->at(0);
      const auto ptTTsigMax = tt.sigPtRange->at(1);
      if (track.pt() > ptTTsigMin && track.pt() < ptTTsigMax) {
        vCandForTT.emplace_back(index);
      }
      ++index;
    }
    spectra.fill(HIST("hScaledFT0C_deltaPtRandomCone"), scaledFT0C, randomConePt - areaRC * rho, weight);
    spectra.fill(HIST("hScaledFT0M_deltaPtRandomCone"), scaledFT0M, randomConePt - areaRC * rho, weight);

    //----------------------------------------------------------
    // Avoid leading jet (JE jet reconstruction sorts jets by pT)

    // square of distance to accept RC placement in events with leading jet
    float dMinR2 = std::pow(jet.radius + bkgd.randomConeR + bkgd.minDeltaRToJet, 2);

    // max # of attempts to find a place for RC; to avoid possibility with infinite loop in While cycle
    const int maxAttempts = 15000;

    if (jets.size() > 0) // at least 1 jet
    {
      float leadJetEta = jets.iteratorAt(0).eta();
      float leadJetPhi = jets.iteratorAt(0).phi();

      float dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
      float dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);

      bool isTherePlaceForRC = false;
      for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        if ((dPhiLeadJet + dEtaLeadJet) > dMinR2) {
          isTherePlaceForRC = true;
          break;
        }

        randomConeEta = randGen->Uniform(-trk.etaCut + bkgd.randomConeR, trk.etaCut - bkgd.randomConeR);
        randomConePhi = randGen->Uniform(0.0, constants::math::TwoPI);

        dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
        dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);
      }

      if (isTherePlaceForRC) {
        randomConePt = 0.0;
        for (const auto& track : tracks) {
          if (skipTrack(track))
            continue;

          float dEta = std::pow(randomConeEta - track.eta(), 2);
          float dPhi = std::pow(RecoDecay::constrainAngle(randomConePhi - track.phi(), -constants::math::PI), 2);

          if ((dEta + dPhi) < radiusRC2) // inside RC
          {
            randomConePt += track.pt();
          }
        }
        spectra.fill(HIST("hScaledFT0C_deltaPtRandomConeAvoidLeadJet"), scaledFT0C, randomConePt - areaRC * rho, weight);
        spectra.fill(HIST("hScaledFT0M_deltaPtRandomConeAvoidLeadJet"), scaledFT0M, randomConePt - areaRC * rho, weight);
      }

      //----------------------------------------------------------
      // Place cone perpendicular to the leading jet (perpendicular cone)
      float perpConeEta = leadJetEta;
      float perpConePhi = leadJetPhi + constants::math::PIHalf;

      float perpConePt = 0.0;
      for (const auto& track : tracks) {
        if (skipTrack(track))
          continue;

        float dEta = std::pow(perpConeEta - track.eta(), 2);
        float dPhi = std::pow(RecoDecay::constrainAngle(perpConePhi - track.phi(), -constants::math::PI), 2);

        if ((dEta + dPhi) < radiusRC2) // inside RC
        {
          perpConePt += track.pt();
        }
      }
      spectra.fill(HIST("hScaledFT0C_deltaPtPerpConeAvoidLeadJet"), scaledFT0C, perpConePt - areaRC * rho, weight);
      spectra.fill(HIST("hScaledFT0M_deltaPtPerpConeAvoidLeadJet"), scaledFT0M, perpConePt - areaRC * rho, weight);
    }

    //----------------------------------------------------------
    // Avoid leading and subleading jets
    if (jets.size() > 1) // at least 2 jets in an event
    {

      // Leading jet
      float leadJetEta = jets.iteratorAt(0).eta();
      float leadJetPhi = jets.iteratorAt(0).phi();
      float dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
      float dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);

      // Subleading jet
      float subleadJetEta = jets.iteratorAt(1).eta();
      float subleadJetPhi = jets.iteratorAt(1).phi();
      float dEtaSubleadJet = std::pow(subleadJetEta - randomConeEta, 2);
      float dPhiSubleadJet = std::pow(RecoDecay::constrainAngle(subleadJetPhi - randomConePhi, -constants::math::PI), 2);

      // Try to add events with TTsig
      bool keepEventWithTT = false;
      if (vCandForTT.size() > 0) // at least 1 TT
      {
        auto randIndexTrack = randGen->Integer(vCandForTT.size());
        auto objTT = tracks.iteratorAt(vCandForTT[randIndexTrack]);

        // Skip events where TT is not a part of leading or subleading jets (mutlijet event, difficult to place RC and avoid hard jets)
        if (isTrackInJet(jets.iteratorAt(0), objTT) || isTrackInJet(jets.iteratorAt(1), objTT)) {
          keepEventWithTT = true;
        }
      }

      //----------------------------------------------------------
      bool isTherePlaceForRC = false;
      for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        if ((dPhiLeadJet + dEtaLeadJet) > dMinR2 && (dEtaSubleadJet + dPhiSubleadJet) > dMinR2) {
          isTherePlaceForRC = true;
          break;
        }

        randomConeEta = randGen->Uniform(-trk.etaCut + bkgd.randomConeR, trk.etaCut - bkgd.randomConeR);
        randomConePhi = randGen->Uniform(0.0, constants::math::TwoPI);

        dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
        dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);

        dEtaSubleadJet = std::pow(subleadJetEta - randomConeEta, 2);
        dPhiSubleadJet = std::pow(RecoDecay::constrainAngle(subleadJetPhi - randomConePhi, -constants::math::PI), 2);
      }

      if (isTherePlaceForRC) {
        randomConePt = 0.0;
        for (const auto& track : tracks) {
          if (skipTrack(track))
            continue;

          float dEta = std::pow(randomConeEta - track.eta(), 2);
          float dPhi = std::pow(RecoDecay::constrainAngle(randomConePhi - track.phi(), -constants::math::PI), 2);

          if ((dEta + dPhi) < radiusRC2) // inside RC
          {
            randomConePt += track.pt();
          }
        }
        spectra.fill(HIST("hScaledFT0C_deltaPtRandomConeAvoidLeadAndSubleadJet"), scaledFT0C, randomConePt - areaRC * rho, weight);
        spectra.fill(HIST("hScaledFT0M_deltaPtRandomConeAvoidLeadAndSubleadJet"), scaledFT0M, randomConePt - areaRC * rho, weight);

        if (keepEventWithTT) {
          spectra.fill(HIST("hScaledFT0C_deltaPtRandomConeInEventTTSig"), scaledFT0C, randomConePt - areaRC * rho, weight);
          spectra.fill(HIST("hScaledFT0M_deltaPtRandomConeInEventTTSig"), scaledFT0M, randomConePt - areaRC * rho, weight);
        }
      }
    }
  }

  template <typename JCollision, typename JParticles, typename Jets>
  void fillBkgdFluctuationsMCPartLevel(JCollision const& collision,
                                       Jets const& jets,
                                       JParticles const& particles,
                                       float weight = 1.)
  {
    //----------------------------------------------------------
    float rho = collision.rho();
    float scaledFT0C = getScaledFT0(collision.multFT0C(), ft0c.meanPartLevel);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), ft0a.meanPartLevel), scaledFT0C);

    //----------------------------------------------------------
    // Study bkgd fluctuations in events with TTsig
    std::vector<uint64_t> vCandForTT;

    //----------------------------------------------------------
    // Place absolutely random cone (anywhere in an event)
    float randomConeEta = randGen->Uniform(-trk.etaCut + bkgd.randomConeR, trk.etaCut - bkgd.randomConeR);
    float randomConePhi = randGen->Uniform(0.0, constants::math::TwoPI);
    float radiusRC2 = std::pow(bkgd.randomConeR, 2);
    float areaRC = constants::math::PI * radiusRC2;
    float randomConePt = 0.0;

    uint64_t index = 0;
    for (const auto& particle : particles) {
      if (skipParticle(particle))
        continue;

      float dEta = std::pow(randomConeEta - particle.eta(), 2);
      float dPhi = std::pow(RecoDecay::constrainAngle(randomConePhi - particle.phi(), -constants::math::PI), 2);

      if ((dEta + dPhi) < radiusRC2) // inside RC
      {
        randomConePt += particle.pt();
      }

      // Search for TT_Sig candidate
      const auto ptTTsigMin = tt.sigPtRange->at(0);
      const auto ptTTsigMax = tt.sigPtRange->at(1);
      if (particle.pt() > ptTTsigMin && particle.pt() < ptTTsigMax) {
        vCandForTT.emplace_back(index);
      }
      ++index;
    }
    spectra.fill(HIST("hScaledFT0C_deltaPtRandomCone_PartLevel"), scaledFT0C, randomConePt - areaRC * rho, weight);
    spectra.fill(HIST("hScaledFT0M_deltaPtRandomCone_PartLevel"), scaledFT0M, randomConePt - areaRC * rho, weight);

    //----------------------------------------------------------
    // Avoid leading jet (JE jet reconstruction sorts jets by pT)

    // square of distance to accept RC placement in events with leading jet
    float dMinR2 = std::pow(jet.radius + bkgd.randomConeR + bkgd.minDeltaRToJet, 2);

    // max # of attempts to find a place for RC; to avoid possibility with infinite loop in While cycle
    const int maxAttempts = 15000;

    if (jets.size() > 0) // at least 1 jet
    {
      float leadJetEta = jets.iteratorAt(0).eta();
      float leadJetPhi = jets.iteratorAt(0).phi();

      float dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
      float dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);

      bool isTherePlaceForRC = false;
      for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        if ((dPhiLeadJet + dEtaLeadJet) > dMinR2) {
          isTherePlaceForRC = true;
          break;
        }

        randomConeEta = randGen->Uniform(-trk.etaCut + bkgd.randomConeR, trk.etaCut - bkgd.randomConeR);
        randomConePhi = randGen->Uniform(0.0, constants::math::TwoPI);

        dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
        dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);
      }

      if (isTherePlaceForRC) {
        randomConePt = 0.0;
        for (const auto& particle : particles) {
          if (skipParticle(particle))
            continue;

          float dEta = std::pow(randomConeEta - particle.eta(), 2);
          float dPhi = std::pow(RecoDecay::constrainAngle(randomConePhi - particle.phi(), -constants::math::PI), 2);

          if ((dEta + dPhi) < radiusRC2) // inside RC
          {
            randomConePt += particle.pt();
          }
        }
        spectra.fill(HIST("hScaledFT0C_deltaPtRandomConeAvoidLeadJet_PartLevel"), scaledFT0C, randomConePt - areaRC * rho, weight);
        spectra.fill(HIST("hScaledFT0M_deltaPtRandomConeAvoidLeadJet_PartLevel"), scaledFT0M, randomConePt - areaRC * rho, weight);
      }

      //----------------------------------------------------------
      // Place cone perpendicular to the leading jet (perpendicular cone)
      float perpConeEta = leadJetEta;
      float perpConePhi = leadJetPhi + constants::math::PIHalf;

      float perpConePt = 0.0;
      for (const auto& particle : particles) {
        if (skipParticle(particle))
          continue;

        float dEta = std::pow(perpConeEta - particle.eta(), 2);
        float dPhi = std::pow(RecoDecay::constrainAngle(perpConePhi - particle.phi(), -constants::math::PI), 2);

        if ((dEta + dPhi) < radiusRC2) // inside RC
        {
          perpConePt += particle.pt();
        }
      }
      spectra.fill(HIST("hScaledFT0C_deltaPtPerpConeAvoidLeadJet_PartLevel"), scaledFT0C, perpConePt - areaRC * rho, weight);
      spectra.fill(HIST("hScaledFT0M_deltaPtPerpConeAvoidLeadJet_PartLevel"), scaledFT0M, perpConePt - areaRC * rho, weight);
    }

    //----------------------------------------------------------
    // Avoid leading and subleading jets
    if (jets.size() > 1) // at least 2 jets in an event
    {

      // Leading jet
      float leadJetEta = jets.iteratorAt(0).eta();
      float leadJetPhi = jets.iteratorAt(0).phi();
      float dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
      float dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);

      // Subleading jet
      float subleadJetEta = jets.iteratorAt(1).eta();
      float subleadJetPhi = jets.iteratorAt(1).phi();
      float dEtaSubleadJet = std::pow(subleadJetEta - randomConeEta, 2);
      float dPhiSubleadJet = std::pow(RecoDecay::constrainAngle(subleadJetPhi - randomConePhi, -constants::math::PI), 2);

      // Try to add events with TTsig
      bool keepEventWithTT = false;
      if (vCandForTT.size() > 0) // at least 1 TT
      {
        auto randIndexParticle = randGen->Integer(vCandForTT.size());
        auto objTT = particles.iteratorAt(vCandForTT[randIndexParticle]);

        // Skip events where TT is not a part of leading or subleading jets (mutlijet event, difficult to place RC and avoid hard jets)
        if (isTrackInJet(jets.iteratorAt(0), objTT) || isTrackInJet(jets.iteratorAt(1), objTT)) {
          keepEventWithTT = true;
        }
      }

      //----------------------------------------------------------
      bool isTherePlaceForRC = false;
      for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        if ((dPhiLeadJet + dEtaLeadJet) > dMinR2 && (dEtaSubleadJet + dPhiSubleadJet) > dMinR2) {
          isTherePlaceForRC = true;
          break;
        }

        randomConeEta = randGen->Uniform(-trk.etaCut + bkgd.randomConeR, trk.etaCut - bkgd.randomConeR);
        randomConePhi = randGen->Uniform(0.0, constants::math::TwoPI);

        dEtaLeadJet = std::pow(leadJetEta - randomConeEta, 2);
        dPhiLeadJet = std::pow(RecoDecay::constrainAngle(leadJetPhi - randomConePhi, -constants::math::PI), 2);

        dEtaSubleadJet = std::pow(subleadJetEta - randomConeEta, 2);
        dPhiSubleadJet = std::pow(RecoDecay::constrainAngle(subleadJetPhi - randomConePhi, -constants::math::PI), 2);
      }

      if (isTherePlaceForRC) {
        randomConePt = 0.0;
        for (const auto& particle : particles) {
          if (skipParticle(particle))
            continue;

          float dEta = std::pow(randomConeEta - particle.eta(), 2);
          float dPhi = std::pow(RecoDecay::constrainAngle(randomConePhi - particle.phi(), -constants::math::PI), 2);

          if ((dEta + dPhi) < radiusRC2) // inside RC
          {
            randomConePt += particle.pt();
          }
        }
        spectra.fill(HIST("hScaledFT0C_deltaPtRandomConeAvoidLeadAndSubleadJet_PartLevel"), scaledFT0C, randomConePt - areaRC * rho, weight);
        spectra.fill(HIST("hScaledFT0M_deltaPtRandomConeAvoidLeadAndSubleadJet_PartLevel"), scaledFT0M, randomConePt - areaRC * rho, weight);

        if (keepEventWithTT) {
          spectra.fill(HIST("hScaledFT0C_deltaPtRandomConeInEventTTSig_PartLevel"), scaledFT0C, randomConePt - areaRC * rho, weight);
          spectra.fill(HIST("hScaledFT0M_deltaPtRandomConeInEventTTSig_PartLevel"), scaledFT0M, randomConePt - areaRC * rho, weight);
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  // Block of Process Functions

  //=============================================================================
  //  Recoil jet analysis
  //=============================================================================
  void processData(CollRhoDataIt const& collision,
                   TrackTbl const& tracksPerColl,
                   JetsDataTbl const& jetsPerColl)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);

    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5); // number of events selected for analysis

    fillHistograms(collision, jetsPerColl, tracksPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processData, "process raw data", true);

  //____________________
  void processMCDetLevel(CollRhoDataIt const& collision,
                         TrackTbl const& tracksPerColl,
                         JetsDetTbl const& jetsPerColl)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    spectra.fill(HIST("hEventSelectionCount"), 4.5); // number of events selected for analysis
    fillHistograms(collision, jetsPerColl, tracksPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevel, "process MC det. level data (no weight; MB events)", false);

  //____________________________
  void processMCDetLevelWeighted(CollRhoOutlierDetIt const& collision,
                                 aod::JetMcCollisions const&,
                                 TrackTbl const& tracksPerColl,
                                 JetsDetTbl const& jetsPerColl)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    if (collision.isOutlier()) {
      spectra.fill(HIST("hEventSelectionCount"), 2.5);
      return;
    }

    if (!collision.has_mcCollision()) {
      spectra.fill(HIST("hEventSelectionCount"), 3.5);
      return;
    }

    spectra.fill(HIST("hEventSelectionCount"), 4.5); // number of events selected for analysis
    auto weight = collision.mcCollision().weight();
    fillHistograms(collision, jetsPerColl, tracksPerColl, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevelWeighted, "process MC det. level data (weighted JJ)", false);

  //_____________________
  void processMCPartLevel(CollRhoPartIt const& collision,
                          PartTbl const& particlesPerColl,
                          JetsPartTbl const& jetsPerColl)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);

    if (skipMCEvent(collision)) {
      spectra.fill(HIST("hEventSelectionCountPartLevel"), 1.5);
      return;
    }

    spectra.fill(HIST("hEventSelectionCountPartLevel"), 3.5); // number of events selected for analysis
    fillHistogramsMCPartLevel(collision, jetsPerColl, particlesPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevel, "process MC part. level data (no weight; MB events)", false);

  //_____________________________
  void processMCPartLevelWeighted(CollRhoOutlierPartIt const& collision,
                                  PartTbl const& particlesPerColl,
                                  JetsPartTbl const& jetsPerColl)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);

    if (skipMCEvent(collision)) {
      spectra.fill(HIST("hEventSelectionCountPartLevel"), 1.5);
      return;
    }

    if (collision.isOutlier()) {
      spectra.fill(HIST("hEventSelectionCountPartLevel"), 2.5);
      return;
    }

    spectra.fill(HIST("hEventSelectionCountPartLevel"), 3.5); // number of events selected for analysis

    auto weight = collision.weight();
    spectra.fill(HIST("ptHat"), collision.ptHard(), weight);
    fillHistogramsMCPartLevel(collision, jetsPerColl, particlesPerColl, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevelWeighted, "process MC part. level data (weighted JJ)", false);

  //=============================================================================
  // Construction of response matrix
  //=============================================================================
  void processJetsGeoMatching(CollRhoDetIt const& collision,
                              CollRhoPartTbl const& mcCollisions,
                              TrackTbl const& tracksPerColl,
                              MatchedJetsDetToPartTbl const& mcDetJetsPerColl,
                              MatchedJetsPartToDetTbl const& mcPartJets)
  {
    // Skip detector level collisions
    if (skipEvent(collision) || !collision.has_mcCollision())
      return;

    auto mcCollisionId = collision.mcCollisionId();
    auto detLevelCollRho = collision.rho();
    auto partLevelCollRho = mcCollisions.iteratorAt(mcCollisionId).rho();

    // Slice for mc part level jets associated to a given mcCollisionId
    auto mcPartJetsPerMcCollision = mcPartJets.sliceBy(partJetsPerMcCollision, mcCollisionId); // signature: (__column to slice___, __index__)

    fillMatchedGeoHistograms(mcPartJetsPerMcCollision, mcDetJetsPerColl, tracksPerColl, partLevelCollRho, detLevelCollRho);
  }
  PROCESS_SWITCH(RecoilJets, processJetsGeoMatching, "process matching of MC jets using Geo criterion (no weight; MB events)", false);

  //___________________________
  void processJetsGeoPtMatching(CollRhoDetIt const& collision,
                                CollRhoPartTbl const& mcCollisions,
                                TrackTbl const& tracksPerColl,
                                MatchedJetsDetToPartTbl const& mcDetJetsPerColl,
                                MatchedJetsPartToDetTbl const& mcPartJets)
  {
    // Skip detector level collisions
    if (skipEvent(collision) || !collision.has_mcCollision())
      return;

    auto mcCollisionId = collision.mcCollisionId();
    auto detLevelCollRho = collision.rho();
    auto partLevelCollRho = mcCollisions.iteratorAt(mcCollisionId).rho();

    // Slice for mc part level jets associated to a given mcCollisionId
    auto mcPartJetsPerMcCollision = mcPartJets.sliceBy(partJetsPerMcCollision, mcCollisionId); // signature: (__column to slice___, __index__)

    fillMatchedGeoPtHistograms(mcPartJetsPerMcCollision, mcDetJetsPerColl, tracksPerColl, partLevelCollRho, detLevelCollRho);
  }
  PROCESS_SWITCH(RecoilJets, processJetsGeoPtMatching, "process matching of MC jets using Geo+Pt criteria (no weight; MB events)", false);

  //_________________________________
  void processJetsGeoMatchingWeighted(CollRhoOutlierDetIt const& collision,
                                      CollRhoOutlierPartTbl const& mcCollisions,
                                      TrackTbl const& tracksPerColl,
                                      MatchedJetsDetToPartTbl const& mcDetJetsPerColl,
                                      MatchedJetsPartToDetTbl const& mcPartJets)
  {
    // Skip detector level collisions
    if (skipEvent(collision) || collision.isOutlier() || !collision.has_mcCollision())
      return;

    auto mcCollisionId = collision.mcCollisionId();
    auto detLevelCollRho = collision.rho();
    auto partLevelCollRho = mcCollisions.iteratorAt(mcCollisionId).rho();
    auto weight = mcCollisions.iteratorAt(mcCollisionId).weight();

    // Slice for mc part level jets associated to a given mcCollisionId
    auto mcPartJetsPerMcCollision = mcPartJets.sliceBy(partJetsPerMcCollision, mcCollisionId); // signature: (__column to slice___, __index__)

    fillMatchedGeoHistograms(mcPartJetsPerMcCollision, mcDetJetsPerColl, tracksPerColl, partLevelCollRho, detLevelCollRho, weight);
  }
  PROCESS_SWITCH(RecoilJets, processJetsGeoMatchingWeighted, "process matching of MC jets using Geo criterion (weighted JJ)", false);

  //___________________________________
  void processJetsGeoPtMatchingWeighted(CollRhoOutlierDetIt const& collision,
                                        CollRhoOutlierPartTbl const& mcCollisions,
                                        TrackTbl const& tracksPerColl,
                                        MatchedJetsDetToPartTbl const& mcDetJetsPerColl,
                                        MatchedJetsPartToDetTbl const& mcPartJets)
  {
    // Skip detector level collisions
    if (skipEvent(collision) || collision.isOutlier() || !collision.has_mcCollision())
      return;

    auto mcCollisionId = collision.mcCollisionId();
    auto detLevelCollRho = collision.rho();
    auto partLevelCollRho = mcCollisions.iteratorAt(mcCollisionId).rho();
    auto weight = mcCollisions.iteratorAt(mcCollisionId).weight();

    // Slice for mc part level jets associated to a given mcCollisionId
    auto mcPartJetsPerMcCollision = mcPartJets.sliceBy(partJetsPerMcCollision, mcCollisionId); // signature: (__column to slice___, __index__)

    fillMatchedGeoPtHistograms(mcPartJetsPerMcCollision, mcDetJetsPerColl, tracksPerColl, partLevelCollRho, detLevelCollRho, weight);
  }
  PROCESS_SWITCH(RecoilJets, processJetsGeoPtMatchingWeighted, "process matching of MC jets using Geo+Pt criteria (weighted JJ)", false);

  //=============================================================================
  // Event Activity analysis in OO collisions (raw and MC detector level (no weight; MB events))
  //=============================================================================
  void processEventActivityOO(EvMultZDCDataIt const& collision)
  {
    if (skipEvent(collision))
      return;

    fillMultiplicityHistogramsOO(collision);
  }
  PROCESS_SWITCH(RecoilJets, processEventActivityOO, "process event activity in OO collisions and MC det. level (no weight; MB events)", false);

  //___________________________________________
  void processEventActivityMCDetLevelWeightedOO(EvMultOutlierZDCDetIt const& collision,
                                                aod::JetMcCollisions const&)
  {
    if (skipEvent(collision) || collision.isOutlier() || !collision.has_mcCollision())
      return;

    auto weight = collision.mcCollision().weight();
    fillMultiplicityHistogramsOO(collision, weight);
  }
  PROCESS_SWITCH(RecoilJets, processEventActivityMCDetLevelWeightedOO, "process event activity in MC det. level OO collisions (weighted JJ)", false);

  //=============================================================================
  // Event Activity analysis in OO and pp collisions at Particle level
  //=============================================================================
  void processEventActivityMCPartLevel(CollPartIt const& collision)
  {
    if (skipMCEvent(collision))
      return;

    fillMultiplicityHistogramsMCPartLevel(collision);
  }
  PROCESS_SWITCH(RecoilJets, processEventActivityMCPartLevel, "process event activity in MC part. level events (no weight; MB events)", false);

  //__________________________________________
  void processEventActivityMCPartLevelWeighted(EvMultOutlierPartIt const& collision)
  {
    if (skipMCEvent(collision) || collision.isOutlier())
      return;

    auto weight = collision.weight();
    fillMultiplicityHistogramsMCPartLevel(collision, weight);
  }
  PROCESS_SWITCH(RecoilJets, processEventActivityMCPartLevelWeighted, "process event activity in MC part. level events (weighted JJ)", false);

  //=============================================================================
  // Event Activity QA analysis in OO collisions (raw and MC detector level (no weight; MB events))
  //=============================================================================
  void processEventActivityQA(CollEvSelExtendedIt const& collision,
                              BCsRun3Tbl const& BCs,
                              aod::Zdcs const& ZDCs)
  {
    // Base flag for event selection
    if (!collision.sel8())
      return;

    fillMultiplicityQA(collision, BCs, ZDCs);
  }
  PROCESS_SWITCH(RecoilJets, processEventActivityQA, "process function for EA QA purposes in raw and MC det. level (no weight; MB events) data", false);

  //=============================================================================
  // Di-hadron azimuthal correlation in raw and MC det. level (no weight; MB events) data
  //=============================================================================
  void processLeadingAndAssociatedTracksTask(CollDataIt const& collision,
                                             TrackTbl const& tracksPerColl)
  {
    if (skipEvent(collision))
      return;
    fillLeadingAndAssociatedTracksTask(collision, tracksPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processLeadingAndAssociatedTracksTask, "process di-hadron azimuthal correlation in raw and MC det. level (no weight; MB events) data", false);

  //=============================================================================
  // Estimation of bkgd fluctuations
  //=============================================================================
  void processBkgdFluctuations(CollRhoDataIt const& collision,
                               TrackTbl const& tracksPerColl,
                               JetsDataTbl const& jetsPerColl)
  {
    if (skipEvent(collision))
      return;

    fillBkgdFluctuations(collision, jetsPerColl, tracksPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processBkgdFluctuations, "process raw data to estimate bkgd fluctuations", false);

  //____________________________________
  void processBkgdFluctuationsMCDetLevel(CollRhoDataIt const& collision,
                                         TrackTbl const& tracksPerColl,
                                         JetsDetTbl const& jetsPerColl)
  {
    if (skipEvent(collision))
      return;

    fillBkgdFluctuations(collision, jetsPerColl, tracksPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processBkgdFluctuationsMCDetLevel, "process MC det. level (no weight; MB events) data to estimate bkgd fluctuations", false);

  //____________________________________________
  void processBkgdFluctuationsMCDetLevelWeighted(CollRhoOutlierDetIt const& collision,
                                                 aod::JetMcCollisions const&,
                                                 TrackTbl const& tracksPerColl,
                                                 JetsDetTbl const& jetsPerColl)
  {
    if (skipEvent(collision) || collision.isOutlier() || !collision.has_mcCollision())
      return;

    auto weight = collision.mcCollision().weight();
    fillBkgdFluctuations(collision, jetsPerColl, tracksPerColl, weight);
  }
  PROCESS_SWITCH(RecoilJets, processBkgdFluctuationsMCDetLevelWeighted, "process MC det. level (weighted JJ) data to estimate bkgd fluctuations", false);

  //_____________________________________
  void processBkgdFluctuationsMCPartLevel(CollRhoPartIt const& collision,
                                          PartTbl const& particlesPerColl,
                                          JetsPartTbl const& jetsPerColl)
  {
    if (skipMCEvent(collision))
      return;

    fillBkgdFluctuationsMCPartLevel(collision, jetsPerColl, particlesPerColl);
  }
  PROCESS_SWITCH(RecoilJets, processBkgdFluctuationsMCPartLevel, "process MC part. level (no weight; MB events) data to estimate bkgd fluctuations", false);

  //_____________________________________________
  void processBkgdFluctuationsMCPartLevelWeighted(CollRhoOutlierPartIt const& collision,
                                                  PartTbl const& particlesPerColl,
                                                  JetsPartTbl const& jetsPerColl)
  {
    if (skipMCEvent(collision) || collision.isOutlier())
      return;

    auto weight = collision.weight();
    fillBkgdFluctuationsMCPartLevel(collision, jetsPerColl, particlesPerColl, weight);
  }
  PROCESS_SWITCH(RecoilJets, processBkgdFluctuationsMCPartLevelWeighted, "process MC part. level (weighted JJ) data to estimate bkgd fluctuations", false);

  //------------------------------------------------------------------------------
  // Auxiliary functions
  template <typename Collision>
  bool skipEvent(const Collision& coll)
  {
    /// \brief: trigger cut is needed for pp data
    return !jetderiveddatautilities::selectCollision(coll, eventSelectionBits, ev.skipMBGapEvents, rct.enable, rct.label, rct.rejectLimitedAcceptance, rct.requireZDC) || !jetderiveddatautilities::selectTrigger(coll, triggerMaskBits);
  }

  template <typename Collision>
  bool skipMCEvent(const Collision& coll)
  {
    return !jetderiveddatautilities::selectMcCollision(coll, ev.skipMBGapEvents, rct.enable, rct.label, rct.rejectLimitedAcceptance, rct.requireZDC);
  }

  template <typename Track>
  bool skipTrack(const Track& track)
  {
    return !jetderiveddatautilities::selectTrack(track, trackSelection);
  }

  template <typename Particle>
  bool skipParticle(const Particle& particle)
  {
    auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      return true;
    }
    const bool bParticleNeutral = (static_cast<int8_t>(pdgParticle->Charge()) == 0);
    return bParticleNeutral || !particle.isPhysicalPrimary();
  }

  template <typename Jet>
  std::tuple<double, bool> isRecoilJet(const Jet& jet, double phiTT)
  {
    double dphi = std::fabs(RecoDecay::constrainAngle(jet.phi() - phiTT, -constants::math::PI));
    return {dphi, (constants::math::PI - tt.recoilRegion) < dphi};
  }

  double getPhiTT(const std::vector<double>& vPhiOfTT)
  {
    auto iTrig = randGen->Integer(vPhiOfTT.size());
    return vPhiOfTT[iTrig];
  }

  float getScaledFT0(const float& multFT0, const float& meanFT0)
  {
    return multFT0 / meanFT0;
  }

  float getScaledFT0M(const float& scaledMultFT0A, const float& scaledMultFT0C)
  {
    return 0.5 * (scaledMultFT0A + scaledMultFT0C);
  }

  template <typename Tracks, typename Jet>
  bool isJetWithHighPtConstituent(Jet const& chJet)
  {
    bool bIsJetWithHighPtConstituent = false;
    for (const auto& chJetConstituent : chJet.template tracks_as<Tracks>()) {
      if (chJetConstituent.pt() > jet.constituentPtMax) {
        bIsJetWithHighPtConstituent = true;
        break;
      }
    }
    return bIsJetWithHighPtConstituent;
  }

  template <typename Jet, typename Track>
  bool isTrackInJet(Jet const& jet, Track const& track)
  {
    for (auto const& constituentId : jet.tracksIds()) {
      if (constituentId == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  template <typename histo>
  int getBinNumberOnYaxisForGivenRun(const std::shared_ptr<histo>& histogram, int runNumber)
  {
    int nBins = histogram->GetYaxis()->GetNbins();
    int binNumber = -1;

    for (int iBin = 1; iBin <= nBins; ++iBin) {
      const char* binLabel = histogram->GetYaxis()->GetBinLabel(iBin);
      if (std::stoi(binLabel) == runNumber) {
        binNumber = iBin;
        break;
      }
    }

    if (binNumber == -1) // No bin found
      return 0;

    return binNumber;
  }

  template <typename JTracksTable>
  std::vector<double> getPhiOfAllTTsigCandidates(JTracksTable const& tracks)
  {
    std::vector<double> vPhiOfTT;

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      // Search for TT_Sig candidate
      const auto ptTTsigMin = tt.sigPtRange->at(0);
      const auto ptTTsigMax = tt.sigPtRange->at(1);
      if (track.pt() > ptTTsigMin && track.pt() < ptTTsigMax) {
        vPhiOfTT.emplace_back(track.phi());
      }
    }
    return vPhiOfTT;
  }

  template <typename typeHist>
  void setBinLablesYZaxes(const std::shared_ptr<typeHist>& histPointer,
                          const std::vector<const char*>& yAxis,
                          const std::vector<const char*>& zAxis)
  {
    const int nRunsOO = yAxis.size();
    const int nEvSelFlags = zAxis.size();
    for (int iRun = 0; iRun < nRunsOO; ++iRun) {
      histPointer->GetYaxis()->SetBinLabel(iRun + 1, yAxis[iRun]);
    }

    const int dimOf2DHist = 2;
    bool isHist2D = histPointer->GetDimension() == dimOf2DHist;
    if (isHist2D)
      return;

    for (int iFlag = 0; iFlag < nEvSelFlags; ++iFlag) {
      histPointer->GetZaxis()->SetBinLabel(iFlag + 1, zAxis[iFlag]);
    }
  }

  template <typename MassiveA, typename MassiveB>
  std::unordered_set<int32_t> intersectMatchIds(MassiveA const& geoIds, MassiveB const& ptIds)
  {
    std::unordered_set<int32_t> geoSet;
    geoSet.reserve(geoIds.size());
    for (const auto& id : geoIds) {
      if (id >= 0) {
        geoSet.insert(id);
      }
    }

    std::unordered_set<int32_t> bothSet;
    bothSet.reserve(std::min(geoIds.size(), ptIds.size()));
    for (const auto& id : ptIds) {
      if (id >= 0 && geoSet.contains(id)) {
        bothSet.insert(id);
      }
    }
    return bothSet;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<RecoilJets>(cfgc)};
}
