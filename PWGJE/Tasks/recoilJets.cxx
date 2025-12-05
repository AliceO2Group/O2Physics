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

#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Shorthand notations
using FilteredColl = soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
using FilteredCollPartLevel = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>>::iterator;
using FilteredCollDetLevelGetWeight = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos, aod::JCollisionOutliers>>::iterator;
using FilteredEventMultiplicity = soa::Filtered<soa::Join<aod::JetCollisions, aod::ZDCMults>>::iterator;
using FilteredEventMultiplicityDetLevelGetWeight = soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers, aod::ZDCMults>>::iterator;
using FilteredEventMultiplicityPartLevel = soa::Filtered<soa::Join<aod::JetMcCollisions, aod::JMcCollisionOutliers>>::iterator;

using FilteredJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using FilteredJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>>;
using FilteredJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>;

using FilteredMatchedJetsDetLevel = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using FilteredMatchedJetsPartLevel = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

using FilteredTracks = soa::Filtered<aod::JetTracks>;
using FilteredParticles = soa::Filtered<aod::JetParticles>;

using ColEvSelEA = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FT0MultZeqs, aod::MultsExtra, aod::PVMults>>::iterator;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>; // aod::Run3MatchedToBCExclusive

struct RecoilJets {

  // List of configurable parameters
  Configurable<std::string> evSel{"evSel", "sel8", "Choose event selection"};
  Configurable<std::string> trkSel{"trkSel", "globalTracks", "Set track selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10., "Accepted z-vertex range"};
  Configurable<float> fracSig{"fracSig", 0.9, "Fraction of events to use for signal TT"};

  Configurable<float> trkPtMin{"trkPtMin", 0.15, "Minimum pT of acceptanced tracks"};
  Configurable<float> trkPtMax{"trkPtMax", 100., "Maximum pT of acceptanced tracks"};

  Configurable<float> trkEtaCut{"trkEtaCut", 0.9, "Eta acceptance of TPC"};
  Configurable<float> jetR{"jetR", 0.4, "Jet cone radius"};
  Configurable<float> maxJetConstituentPt{"maxJetConstituentPt", 100., "Remove jets with constituent above this pT cut"};

  Configurable<std::string> triggerMasks{"triggerMasks", "", "Relevant trigger masks: fTrackLowPt,fTrackHighPt"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false,
                                     "flag to choose to reject min. bias gap events; jet-level rejection "
                                     "applied at the jet finder level, here rejection is applied for "
                                     "collision and track process functions"};

  Configurable<float> meanFT0A{"meanFT0A", -1., "Mean value of FT0A signal"};
  Configurable<float> meanFT0C{"meanFT0C", -1., "Mean value of FT0C signal"};

  Configurable<float> meanZeqFT0A{"meanZeqFT0A", -1., "Mean value of equalized FT0A signal"};
  Configurable<float> meanZeqFT0C{"meanZeqFT0C", -1., "Mean value of equalized FT0C signal"};

  Configurable<float> meanFT0APartLevel{"meanFT0APartLevel", -1., "Mean number of charged part. within FT0A acceptance"};
  Configurable<float> meanFT0CPartLevel{"meanFT0CPartLevel", -1., "Mean number of charged part. within FT0C acceptance"};

  // Parameters for recoil jet selection
  Configurable<std::vector<float>> ptTTref{"ptTTref", {5., 7}, "Transverse momentum (min, max) range for reference TT"};
  Configurable<std::vector<float>> ptTTsig{"ptTTsig", {20., 50}, "Transverse momentum (min, max) range for signal TT"};
  Configurable<float> recoilRegion{"recoilRegion", 0.6, "Width of recoil acceptance"};
  Configurable<std::vector<float>> phiRestrTTSelection{"phiRestrTTSelection", {0., 6.3}, "Restriction on phi angle (min, max) to search for TT"};

  // Leading track and associated track
  Configurable<std::vector<float>> pTLeadTrack{"pTLeadTrack", {4., 6.}, "Transverse momenturm range (min, max) for leading tracks"};
  Configurable<float> pTAssociatTrackMin{"pTAssociatTrackMin", 2., "Min transverse momenturm for associated tracks"};

  // List of configurable parameters for histograms
  Configurable<uint16_t> histJetPt{"histJetPt", 100, "Maximum value of jet pT stored in histograms"};
  Configurable<uint16_t> histMultBins{"histMultBins", 600, "Number of bins for scaled FT0M multiplicity"};
  Configurable<uint16_t> histZDCTimeBins{"histZDCTimeBins", 240, "Number of bins for ZDC timing histograms"};

  // Axes specification
  ConfigurableAxis multFT0CThresh{"multFT0CThresh", {VARIABLE_WIDTH, 0.0, 0.133, 0.233, 0.367, 0.567, 0.767, 1.067, 1.4, 1.867, 2.5, 3.9, 5.4, 6.9, 20.}, "Percentiles of scaled FT0C: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"};   // default values for raw data
  ConfigurableAxis multFT0MThresh{"multFT0MThresh", {VARIABLE_WIDTH, 0.0, 0.167, 0.267, 0.4, 0.567, 0.8, 1.067, 1.4, 1.833, 2.433, 3.667, 5.1, 6.433, 20.}, "Percentiles of scaled FT0M: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"}; // default values for raw data

  ConfigurableAxis multFT0CThreshPartLevel{"multFT0CThreshPartLevel", {VARIABLE_WIDTH, 0.0, 0.133, 0.233, 0.367, 0.567, 0.767, 1.067, 1.4, 1.867, 2.5, 3.9, 5.4, 6.9, 20.}, "Percentiles of scaled FT0C: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"};
  ConfigurableAxis multFT0MThreshPartLevel{"multFT0MThreshPartLevel", {VARIABLE_WIDTH, 0.0, 0.167, 0.267, 0.4, 0.567, 0.8, 1.067, 1.4, 1.833, 2.433, 3.667, 5.1, 6.433, 20.}, "Percentiles of scaled FT0M: 100%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 20%, 10%, 1%, 0.1%, 0.01%"};

  // Auxiliary variables
  TRandom3* rand = new TRandom3(0);

  // Declare filter on collision Z vertex
  Filter jCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter jCollisionFilterMC = nabs(aod::jmccollision::posZ) < vertexZCut;
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;

  // Declare filters on accepted tracks and MC particles (settings for jet reco are provided in the jet finder wagon)
  Filter trackFilter = aod::jtrack::pt > trkPtMin&& aod::jtrack::pt < trkPtMax&& nabs(aod::jtrack::eta) < trkEtaCut;
  Filter partFilter = nabs(aod::jmcparticle::eta) < trkEtaCut;

  // Declare filter on jets
  Filter jetRadiusFilter = aod::jet::r == nround(jetR.node() * 100.);

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<FilteredMatchedJetsPartLevel> partJetsPerCollision = aod::jet::mcCollisionId;

  void init(InitContext const&)
  {
    // Initialize histogram axes
    AxisSpec pT{histJetPt, 0.0, histJetPt * 1., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPTcorr{histJetPt + 20, -20., histJetPt * 1.0, "#it{p}_{T, jet}^{ch, corr} (GeV/#it{c})"};
    AxisSpec scaledFT0A{histMultBins, 0.0, 20., "FT0A / #LT FT0A #GT"};
    AxisSpec scaledFT0C{histMultBins, 0.0, 20., "FT0C / #LT FT0C #GT"};
    AxisSpec scaledFT0M{histMultBins, 0.0, 20., "FT0M^{*}"};
    AxisSpec zdcTiming{histZDCTimeBins, -30., 30., ""};

    AxisSpec phiAngle{40, 0.0, constants::math::TwoPI, "#it{#varphi} (rad)"};
    AxisSpec deltaPhiAngle{52, 0.0, constants::math::PI, "#Delta#it{#varphi} (rad)"};
    AxisSpec pseudorap{40, -1., 1., "#it{#eta}"};
    AxisSpec pseudorapJets{20, -0.5, 0.5, "#it{#eta}_{jet}"};
    AxisSpec jetArea{50, 0.0, 5., "Area_{jet}"};
    AxisSpec rho{50, 0.0, 50., "#it{#rho}"};

    std::string nameFT0Caxis = "FT0C / #LT FT0C #GT";
    std::string nameFT0Maxis = "FT0M^{*}";

    // Convert configurable strings to std::string
    std::string evSelToString = static_cast<std::string>(evSel);
    std::string trkSelToString = static_cast<std::string>(trkSel);

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(evSelToString);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trkSelToString);
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);

    // List of raw and MC det. distributions
    if (doprocessData || doprocessMCDetLevel || doprocessMCDetLevelWeighted) {
      spectra.add("hEventSelectionCount", "Count # of events in the analysis", kTH1F, {{6, 0.0, 6.}});
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of events");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of events after sel. %s", evSelToString.data()));
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(3, "# of events skipMBGap");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(4, "# of events w. outlier");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(5, "# of events w/o assoc MC.");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(6, "# of selected events");

      spectra.add("hScaledFT0C_vertexZ", "Z vertex of collisions", kTH2F, {{multFT0CThresh, nameFT0Caxis}, {60, -12., 12., "#it{z}_{vertex}"}});
      spectra.add("hScaledFT0M_vertexZ", "Z vertex of collisions", kTH2F, {{multFT0MThresh, nameFT0Maxis}, {60, -12., 12., "#it{z}_{vertex}"}});

      spectra.add("hTrackSelectionCount", "Count # of tracks in the analysis", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hTrackSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of tracks");
      spectra.get<TH1>(HIST("hTrackSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of tracks after sel. %s", trkSelToString.data()));

      spectra.add("hScaledFT0CTrackPtEtaPhi", "Charact. of tracks", kTHnSparseF, {{multFT0CThresh, nameFT0Caxis}, pT, pseudorap, phiAngle});
      spectra.add("hScaledFT0MTrackPtEtaPhi", "Charact. of tracks", kTHnSparseF, {{multFT0MThresh, nameFT0Maxis}, pT, pseudorap, phiAngle});
      spectra.add("hTTSig_pT", "pT spectrum of all found TT_{Sig} cand.", kTH1F, {{40, 10., 50.}}); // needed to distinguish merged data from diff. wagons

      spectra.add("hScaledFT0C_Ntrig", "Total number of selected triggers per class vs scaled FT0C", kTH2F, {{multFT0CThresh, nameFT0Caxis}, {2, 0.0, 2.}});
      spectra.get<TH2>(HIST("hScaledFT0C_Ntrig"))->GetYaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH2>(HIST("hScaledFT0C_Ntrig"))->GetYaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hScaledFT0M_Ntrig", "Total number of selected triggers per class vs scaled FT0M", kTH2F, {{multFT0MThresh, nameFT0Maxis}, {2, 0.0, 2.}});
      spectra.get<TH2>(HIST("hScaledFT0M_Ntrig"))->GetYaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH2>(HIST("hScaledFT0M_Ntrig"))->GetYaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hScaledFT0C_TTRef_per_event", "Number of TT_{Ref} per event vs scaled FT0C", kTH2F, {{multFT0CThresh, nameFT0Caxis}, {15, 0.5, 15.5, "# of TT_{Ref}"}});
      spectra.add("hScaledFT0M_TTRef_per_event", "Number of TT_{Ref} per event vs scaled FT0M", kTH2F, {{multFT0MThresh, nameFT0Maxis}, {15, 0.5, 15.5, "# of TT_{Ref}"}});

      spectra.add("hScaledFT0C_TTSig_per_event", "Number of TT_{Sig} per event vs scaled FT0C", kTH2F, {{multFT0CThresh, nameFT0Caxis}, {10, 0.5, 10.5, "# of TT_{Sig}"}});
      spectra.add("hScaledFT0M_TTSig_per_event", "Number of TT_{Sig} per event vs scaled FT0M", kTH2F, {{multFT0MThresh, nameFT0Maxis}, {10, 0.5, 10.5, "# of TT_{Sig}"}});

      spectra.add("hJetPtEtaPhiRhoArea", "Charact. of inclusive jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThresh, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThresh, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThresh, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThresh, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThresh, nameFT0Caxis}, deltaPhiAngle, pT});
      spectra.add("hScaledFT0M_DPhi_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThresh, nameFT0Maxis}, deltaPhiAngle, pT});

      spectra.add("hScaledFT0C_DPhi_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThresh, nameFT0Caxis}, deltaPhiAngle, pT});
      spectra.add("hScaledFT0M_DPhi_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThresh, nameFT0Maxis}, deltaPhiAngle, pT});

      spectra.add("hScaledFT0C_Recoil_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThresh, nameFT0Caxis}, jetPTcorr});
      spectra.add("hScaledFT0M_Recoil_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThresh, nameFT0Maxis}, jetPTcorr});

      spectra.add("hScaledFT0C_Recoil_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThresh, nameFT0Caxis}, jetPTcorr});
      spectra.add("hScaledFT0M_Recoil_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThresh, nameFT0Maxis}, jetPTcorr});

      spectra.add("hScaledFT0C_Recoil_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThresh, nameFT0Caxis}, pT});
      spectra.add("hScaledFT0M_Recoil_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThresh, nameFT0Maxis}, pT});

      spectra.add("hScaledFT0C_Recoil_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThresh, nameFT0Caxis}, pT});
      spectra.add("hScaledFT0M_Recoil_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThresh, nameFT0Maxis}, pT});

      spectra.add("hJetArea_JetPt_Rho_TTRef", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
      spectra.add("hJetArea_JetPt_Rho_TTSig", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});

      spectra.add("hScaledFT0C_Rho", "Scaled FT0C & #rho", kTH2F, {{multFT0CThresh, nameFT0Caxis}, rho});
      spectra.add("hScaledFT0M_Rho", "Scaled FT0M & #rho", kTH2F, {{multFT0MThresh, nameFT0Maxis}, rho});

      spectra.add("hScaledFT0C_Rho_TTRef", "Events w. TT_{Ref}: scaled FT0C & #rho", kTH2F, {{multFT0CThresh, nameFT0Caxis}, rho});
      spectra.add("hScaledFT0M_Rho_TTRef", "Events w. TT_{Ref}: scaled FT0M & #rho", kTH2F, {{multFT0MThresh, nameFT0Maxis}, rho});

      spectra.add("hScaledFT0C_Rho_TTSig", "Events w. TT_{Sig}: scaled FT0C & #rho", kTH2F, {{multFT0CThresh, nameFT0Caxis}, rho});
      spectra.add("hScaledFT0M_Rho_TTSig", "Events w. TT_{Sig}: scaled FT0M & #rho", kTH2F, {{multFT0MThresh, nameFT0Maxis}, rho});

      spectra.add("hScaledFT0C_TTRef", "Events w. TT_{Ref}: scaled FT0C", kTH1F, {scaledFT0C});
      spectra.add("hScaledFT0M_TTRef", "Events w. TT_{Ref}: scaled FT0M", kTH1F, {scaledFT0M});

      spectra.add("hScaledFT0C_TTSig", "Events w. TT_{Sig}: scaled FT0C", kTH1F, {scaledFT0C});
      spectra.add("hScaledFT0M_TTSig", "Events w. TT_{Sig}: scaled FT0M", kTH1F, {scaledFT0M});

      // Rectricted phi range for TT selection
      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTRef_RectrictedPhi", Form("Events w. TT_{Ref} #in #varphi (%.2f, %.2f): scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0CThresh, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTRef_RectrictedPhi", Form("Events w. TT_{Ref} #in #varphi (%.2f, %.2f): scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0MThresh, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTSig_RectrictedPhi", Form("Events w. TT_{Sig} #in #varphi (%.2f, %.2f): scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0CThresh, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTSig_RectrictedPhi", Form("Events w. TT_{Sig} #in #varphi (%.2f, %.2f): scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0MThresh, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});
    }

    // List of MC particle level distributions
    if (doprocessMCPartLevel || doprocessMCPartLevelWeighted) {
      spectra.add("hScaledFT0C_vertexZMC", "Z vertex of MCcollision", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, {60, -12., 12., "#it{z}_{vertex}"}});
      spectra.add("hScaledFT0M_vertexZMC", "Z vertex of MCcollision", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, {60, -12., 12., "#it{z}_{vertex}"}});
      spectra.add("ptHat", "Distribution of pT hat", kTH1F, {{5000, 0.0, 1000.}});

      spectra.add("hEventSelectionCountPartLevel", "Count # of events in the part. level analysis", kTH1F, {{4, 0.0, 4.}});
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(1, "Total # of events");
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(2, "# of events skipMB gap");
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(3, "# of events w. outlier");
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(4, "# of selected events");

      spectra.add("hScaledFT0CPartPtEtaPhi", "Charact. of particles", kTHnSparseF, {{multFT0CThreshPartLevel, nameFT0Caxis}, pT, pseudorap, phiAngle});
      spectra.add("hScaledFT0MPartPtEtaPhi", "Charact. of particles", kTHnSparseF, {{multFT0MThreshPartLevel, nameFT0Maxis}, pT, pseudorap, phiAngle});

      spectra.add("hScaledFT0C_Ntrig_Part", "Total number of selected triggers per class vs scaled FT0C", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, {2, 0.0, 2.}});
      spectra.get<TH2>(HIST("hScaledFT0C_Ntrig_Part"))->GetXaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH2>(HIST("hScaledFT0C_Ntrig_Part"))->GetXaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hScaledFT0M_Ntrig_Part", "Total number of selected triggers per class vs scaled FT0M", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, {2, 0.0, 2.}});
      spectra.get<TH2>(HIST("hScaledFT0M_Ntrig_Part"))->GetXaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH2>(HIST("hScaledFT0M_Ntrig_Part"))->GetXaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hScaledFT0C_TTRef_per_event_Part", "Number of TT_{Ref} per event vs scaled FT0C", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, {15, 0.5, 15.5, "# of TT_{Ref}"}});
      spectra.add("hScaledFT0M_TTRef_per_event_Part", "Number of TT_{Ref} per event vs scaled FT0M", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, {15, 0.5, 15.5, "# of TT_{Ref}"}});

      spectra.add("hScaledFT0C_TTSig_per_event_Part", "Number of TT_{Sig} per event vs scaled FT0C", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, {10, 0.5, 10.5, "# of TT_{Sig}"}});
      spectra.add("hScaledFT0M_TTSig_per_event_Part", "Number of TT_{Sig} per event vs scaled FT0M", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, {10, 0.5, 10.5, "# of TT_{Sig}"}});

      spectra.add("hJetPtEtaPhiRhoArea_Part", "Charact. of inclusive part. level jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThreshPartLevel, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThreshPartLevel, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThreshPartLevel, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThreshPartLevel, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_TTRef_Part", "Events w. TT_{Ref}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThreshPartLevel, nameFT0Caxis}, deltaPhiAngle, pT});
      spectra.add("hScaledFT0M_DPhi_JetPt_TTRef_Part", "Events w. TT_{Ref}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThreshPartLevel, nameFT0Maxis}, deltaPhiAngle, pT});

      spectra.add("hScaledFT0C_DPhi_JetPt_TTSig_Part", "Events w. TT_{Sig}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0CThreshPartLevel, nameFT0Caxis}, deltaPhiAngle, pT});
      spectra.add("hScaledFT0M_DPhi_JetPt_TTSig_Part", "Events w. TT_{Sig}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {{multFT0MThreshPartLevel, nameFT0Maxis}, deltaPhiAngle, pT});

      spectra.add("hScaledFT0C_Recoil_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, jetPTcorr});
      spectra.add("hScaledFT0M_Recoil_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, jetPTcorr});

      spectra.add("hScaledFT0C_Recoil_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, jetPTcorr});
      spectra.add("hScaledFT0M_Recoil_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, jetPTcorr});

      spectra.add("hScaledFT0C_Recoil_JetPt_TTRef_Part", "Events w. TT_{Ref}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, pT});
      spectra.add("hScaledFT0M_Recoil_JetPt_TTRef_Part", "Events w. TT_{Ref}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, pT});

      spectra.add("hScaledFT0C_Recoil_JetPt_TTSig_Part", "Events w. TT_{Sig}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, pT});
      spectra.add("hScaledFT0M_Recoil_JetPt_TTSig_Part", "Events w. TT_{Sig}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, pT});

      spectra.add("hJetArea_JetPt_Rho_TTRef_Part", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
      spectra.add("hJetArea_JetPt_Rho_TTSig_Part", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});

      spectra.add("hScaledFT0C_Rho_Part", "Scaled FT0C & #rho", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, rho});
      spectra.add("hScaledFT0M_Rho_Part", "Scaled FT0M & #rho", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, rho});

      spectra.add("hScaledFT0C_Rho_TTRef_Part", "Events w. TT_{Ref}: scaled FT0C & #rho", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, rho});
      spectra.add("hScaledFT0M_Rho_TTRef_Part", "Events w. TT_{Ref}: scaled FT0M & #rho", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, rho});

      spectra.add("hScaledFT0C_Rho_TTSig_Part", "Events w. TT_{Sig}: scaled FT0C & #rho", kTH2F, {{multFT0CThreshPartLevel, nameFT0Caxis}, rho});
      spectra.add("hScaledFT0M_Rho_TTSig_Part", "Events w. TT_{Sig}: scaled FT0M & #rho", kTH2F, {{multFT0MThreshPartLevel, nameFT0Maxis}, rho});

      spectra.add("hScaledFT0C_TTRef_Part", "Events w. TT_{Ref}: scaled FT0C", kTH1F, {scaledFT0C});
      spectra.add("hScaledFT0M_TTRef_Part", "Events w. TT_{Ref}: scaled FT0M", kTH1F, {scaledFT0M});

      spectra.add("hScaledFT0C_TTSig_Part", "Events w. TT_{Sig}: scaled FT0C", kTH1F, {scaledFT0C});
      spectra.add("hScaledFT0M_TTSig_Part", "Events w. TT_{Sig}: scaled FT0M", kTH1F, {scaledFT0M});

      // Rectricted phi range for TT selection
      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTRef_RectrictedPhi_Part", Form("Events w. TT_{Ref} #in #varphi (%.2f, %.2f): scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0CThreshPartLevel, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTRef_RectrictedPhi_Part", Form("Events w. TT_{Ref} #in #varphi (%.2f, %.2f): scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0MThreshPartLevel, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTSig_RectrictedPhi_Part", Form("Events w. TT_{Sig} #in #varphi (%.2f, %.2f): scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0CThreshPartLevel, nameFT0Caxis}, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTSig_RectrictedPhi_Part", Form("Events w. TT_{Sig} #in #varphi (%.2f, %.2f): scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", phiRestrTTSelection->at(0), phiRestrTTSelection->at(1)), kTH3F, {{multFT0MThreshPartLevel, nameFT0Maxis}, deltaPhiAngle, jetPTcorr});
    }

    // Jet matching: part. vs. det.
    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      spectra.add("hJetPt_DetLevel_vs_PartLevel", "Correlation jet pT at det. vs. part. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}});
      // spectra.add("hJetPt_Corr_PartLevel_vs_DetLevel", "Correlation jet pT at
      // part. vs. det. levels", kTH2F, {jetPTcorr, jetPTcorr});
      spectra.add("hJetPt_DetLevel_vs_PartLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}});
      // spectra.add("hJetPt_Corr_PartLevel_vs_DetLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {jetPTcorr, jetPTcorr});

      spectra.add("hMissedJets_pT", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hMissedJets_Corr_pT", "Part. level jets w/o matched pair", kTH1F, {jetPTcorr});
      spectra.add("hMissedJets_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hMissedJets_Corr_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {jetPTcorr});

      spectra.add("hFakeJets_pT", "Det. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hFakeJets_Corr_pT", "Det. level jets w/o matched pair", kTH1F, {jetPTcorr});
      spectra.add("hFakeJets_pT_RecoilJets", "Det. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hFakeJets_Corr_pT_RecoilJets", "Det. level jets w/o matched pair", kTH1F, {jetPTcorr});

      spectra.add("hJetPt_resolution", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT});
      spectra.add("hJetPt_resolution_RecoilJets", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT});

      spectra.add("hJetPhi_resolution", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT});
      spectra.add("hJetPhi_resolution_RecoilJets", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT});

      spectra.add("hNumberMatchedJetsPerOneBaseJet", "# of tagged jets per 1 base jet vs. jet pT", kTH2F, {{10, 0.5, 10.5}, {100, 0.0, 100.}});
    }

    // Multiplicity for raw data and detector level MC
    if (doprocessMultiplicityOO || doprocessMultiplicityMCDetLevelWeightedOO) {
      spectra.add("hMultFT0A", "Mult. signal from FTOA", kTH1F, {{2000, 0.0, 40000., "FT0A"}});
      spectra.add("hMultFT0C", "Mult. signal from FTOC", kTH1F, {{2000, 0.0, 40000., "FT0C"}});
      spectra.add("hMultFT0M", "Total mult. signal from FT0A & FTOC", kTH1F, {{3000, 0.0, 60000., "FT0M"}});

      spectra.add("hScaleMultFT0A", "Scaled mult. signal from FTOA", kTH1F, {scaledFT0A});
      spectra.add("hScaleMultFT0C", "Scaled mult. signal from FTOC", kTH1F, {scaledFT0C});
      spectra.add("hScaleMultFT0M", "Scaled total mult. signal from FT0A & FTOC", kTH1F, {scaledFT0M});

      spectra.add("hMultZNA", "Mult. signal from ZDC A-side", kTH1F, {{1000, 0.0, 5000., "ZNA"}});
      spectra.add("hMultZNC", "Mult. signal from ZDC C-side", kTH1F, {{1000, 0.0, 5000., "ZNC"}});
      spectra.add("hMultZNM", "Total mult. signal from ZDCs for neutrons", kTH1F, {{4000, 0.0, 8000., "ZNM"}});

      spectra.add("hMultZPA", "Mult. signal from ZDC A-side", kTH1F, {{1000, 0.0, 5000., "ZPA"}});
      spectra.add("hMultZPC", "Mult. signal from ZDC C-side", kTH1F, {{1000, 0.0, 5000., "ZPC"}});
      spectra.add("hMultZPM", "Total mult. signal from ZDCs for protons", kTH1F, {{4000, 0.0, 8000., "ZPM"}});

      // Correlations
      spectra.add("hZPA_vs_ZNA", "Correlation of signals ZPA vs ZNA", kTH2F, {{1000, 0.0, 5000., "ZPA"}, {1000, 0.0, 5000., "ZNA"}});
      spectra.add("hZPC_vs_ZNC", "Correlation of signals ZPC vs ZNC", kTH2F, {{1000, 0.0, 5000., "ZPC"}, {1000, 0.0, 5000., "ZNC"}});

      spectra.add("hMultFT0A_vs_ZNA", "Correlation of signals FTOA vs ZNA", kTH2F, {{2000, 0.0, 40000., "FT0A"}, {1000, 0.0, 5000., "ZNA"}});
      spectra.add("hMultFT0C_vs_ZNC", "Correlation of signals FTOC vs ZNC", kTH2F, {{2000, 0.0, 40000., "FT0C"}, {1000, 0.0, 5000., "ZNC"}});
      spectra.add("hMultFT0M_vs_ZNM", "Correlation of signals FTOM vs ZNM", kTH2F, {{3000, 0.0, 60000., "FT0M"}, {4000, 0.0, 8000., "ZNM"}});

      spectra.add("hScaleMultFT0A_vs_ZNA", "Correlation of signals FT0A/meanFT0A vs ZNA", kTH2F, {{scaledFT0A}, {1000, 0.0, 5000., "ZNA"}});
      spectra.add("hScaleMultFT0C_vs_ZNC", "Correlation of signals FT0C/meanFT0C vs ZNC", kTH2F, {{scaledFT0C}, {1000, 0.0, 5000., "ZNC"}});
      spectra.add("hScaleMultFT0M_vs_ZNM", "Correlation of signals FT0M^{*} vs ZNM", kTH2F, {{scaledFT0M}, {4000, 0.0, 8000., "ZNM"}});

      spectra.add("hScaleMultFT0A_vs_ZPA", "Correlation of signals FT0A/meanFT0A vs ZPA", kTH2F, {{scaledFT0A}, {1000, 0.0, 5000., "ZPA"}});
      spectra.add("hScaleMultFT0C_vs_ZPC", "Correlation of signals FT0C/meanFT0C vs ZPC", kTH2F, {{scaledFT0C}, {1000, 0.0, 5000., "ZPC"}});
      spectra.add("hScaleMultFT0M_vs_ZPM", "Correlation of signals FT0M^{*} vs ZPM", kTH2F, {{scaledFT0M}, {4000, 0.0, 8000., "ZPM"}});

      spectra.add("hScaleMultFT0M_vs_ZNA_vs_ZNC", "Correlation of signals FT0M^{*} vs ZNA vs ZNC", kTH3F, {{scaledFT0M}, {600, 0.0, 3000., "ZNA"}, {600, 0.0, 3000., "ZNC"}});
      spectra.add("hScaleMultFT0M_vs_ZPA_vs_ZPC", "Correlation of signals FT0M^{*} vs ZPA vs ZPC", kTH3F, {{scaledFT0M}, {600, 0.0, 3000., "ZPA"}, {600, 0.0, 3000., "ZPC"}});
    }

    // Multiplicity for particle level MC
    if (doprocessMultiplicityPartLevelMC || doprocessMultiplicityPartLevelMCWeighted) {
      spectra.add("hMultFT0APartLevel", "# of primary particles within FTOA acceptance", kTH1F, {{2000, 0.0, 500., "FT0A"}});
      spectra.add("hMultFT0CPartLevel", "# of primary particles within FTOC acceptance", kTH1F, {{2000, 0.0, 500., "FT0C"}});
      spectra.add("hMultFT0MPartLevel", "Total # of primary particles from FT0A & FTOC", kTH1F, {{4000, 0.0, 1000., "FT0M"}});

      spectra.add("hScaleMultFT0APartLevel", "Scaled # of primary particles within FTOA acceptance", kTH1F, {scaledFT0A});
      spectra.add("hScaleMultFT0CPartLevel", "Scaled # of primary particles within FTOC acceptance", kTH1F, {scaledFT0C});
      spectra.add("hScaleMultFT0MPartLevel", "Scaled total # of primary particles from FT0A & FTOC", kTH1F, {scaledFT0M});
    }

    if (doprocessMultiplicityQA) {

      // ZNC timing QA
      spectra.add("hTimeCorrZnaZnc", "Correlat. #it{t}_{ZNA} - #it{t}_{ZNC} vs. #it{t}_{ZNA} + #it{t}_{ZNC}", kTH2F, {{1000, -10., 10., "#it{t}_{ZNA} - #it{t}_{ZNC} (ns)"}, {1000, -10., 10., "#it{t}_{ZNA} + #it{t}_{ZNC} (ns)"}});
      spectra.add("hTimeZnaVsZncVsFT0C", "Correlat. #it{t}_{ZNA} (ns) vs. #it{t}_{ZNC} (ns) vs. FT0C/meanFT0C", kTH3F, {{zdcTiming}, {zdcTiming}, {scaledFT0C}});
      spectra.add("hTimeZnaVsZncVsFT0M", "Correlat. #it{t}_{ZNA} (ns) vs. #it{t}_{ZNC} (ns) vs. FT0M^{*}", kTH3F, {{zdcTiming}, {zdcTiming}, {scaledFT0M}});

      // Number of tracks from PV within acceptance |eta| < 0.8
      spectra.add("hScaledFT0C_TracksPV", "Correlat. FT0C/meanFT0C vs. PV tracks", kTH2F, {{scaledFT0C}, {5000, 0., 5000.}});
      spectra.add("hScaledFT0M_TracksPV", "Correlat. FT0M^{*} vs. PV tracks", kTH2F, {{scaledFT0M}, {5000, 0., 5000.}});

      // ITS-only tracks
      spectra.add("hScaledFT0C_ITStracks", "Correlat. FT0C/meanFT0C vs. number of ITS tracks", kTH2F, {{scaledFT0C}, {5000, 0., 5000.}});
      spectra.add("hScaledFT0M_ITStracks", "Correlat. FT0M^{*} vs. number of ITS tracks", kTH2F, {{scaledFT0M}, {5000, 0., 5000.}});

      // Multiplicity equalized for the vertex position with FT0 detector
      spectra.add("hMultZeqFT0A", "Equalized mult. FT0A", kTH1F, {{{2000, 0.0, 40000., "FT0A"}}});
      spectra.add("hMultZeqFT0C", "Equalized mult. FT0C", kTH1F, {{{2000, 0.0, 40000., "FT0C"}}});
      spectra.add("hMultZeqFT0M", "Equalized mult. FT0M", kTH1F, {{{3000, 0.0, 60000., "FT0M"}}});

      spectra.add("hScaledZeqFT0A", "Equalized FT0A/meanFT0A", kTH1F, {{scaledFT0A}});
      spectra.add("hScaledZeqFT0C", "Equalized FT0C/meanFT0C", kTH1F, {{scaledFT0C}});
      spectra.add("hScaledZeqFT0M", "Equalized FT0M^{*}", kTH1F, {{scaledFT0M}});

      // Run-by-run study of EA
      std::vector<const char*> runNumbersOO = {
        "564356", "564359", "564373", "564374", "564387", "564400", "564414", "564430", "564445"};
      const int nRunsOO = runNumbersOO.size();

      std::vector<const char*> evSelFlags = {
        "sel8", "sel8 + IsGoodZvtxFT0vsPV", "sel8 + NoSameBunchPileup", "sel8 + IsGoodZvtxFT0vsPV + NoSameBunchPileup"};
      const int nEvSelFlags = evSelFlags.size();

      // Scaled FT0 signal
      spectra.add("hScaledFT0APerRunPerSetOfFlags", "FT0A/meanFT0A signal per run per set of ev. sel. flags", kTH3F, {{scaledFT0A}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}});
      spectra.add("hScaledFT0CPerRunPerSetOfFlags", "FT0C/meanFT0C signal per run per set of ev. sel. flags", kTH3F, {{scaledFT0C}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}});
      spectra.add("hScaledFT0MPerRunPerSetOfFlags", "FT0M^{*} signal per run per set of ev. sel. flags", kTH3F, {{scaledFT0M}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}});

      // Unscaled FT0 signal; check whether mean value is the same for all runs
      spectra.add("hFT0APerRunPerSetOfFlags", "FT0A signal per run per set of ev. sel. flags", kTH3F, {{2000, 0.0, 40000., "FT0A"}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}});
      spectra.add("hFT0CPerRunPerSetOfFlags", "FT0C signal per run per set of ev. sel. flags", kTH3F, {{2000, 0.0, 40000., "FT0C"}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}});
      spectra.add("hFT0MPerRunPerSetOfFlags", "FT0M signal per run per set of ev. sel. flags", kTH3F, {{3000, 0.0, 60000., "FT0M"}, {nRunsOO, 0., nRunsOO * 1.}, {nEvSelFlags, 0., nEvSelFlags * 1.}});

      // Check whether each BC has FT0 signal
      spectra.add("hIsFT0SignalComeFromCollPerRun", "", kTH2F, {{4, 0., 4.}, {nRunsOO, 0., nRunsOO * 1.}});
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(1, "BC has FT0");
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(2, "BC has not FT0");
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(3, "Coll. w. BC");
      spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetXaxis()->SetBinLabel(4, "Coll. w/o BC");

      // FT0 signal for the case when there is no associated BC
      spectra.add("hFT0AsignalWithoutBC", "", kTH2F, {{2000, 0.0, 40000., "FT0A"}, {nRunsOO, 0., nRunsOO * 1.}});
      spectra.add("hFT0CsignalWithoutBC", "", kTH2F, {{2000, 0.0, 40000., "FT0C"}, {nRunsOO, 0., nRunsOO * 1.}});
      spectra.add("hFT0MsignalWithoutBC", "", kTH2F, {{2000, 0.0, 40000., "FT0M"}, {nRunsOO, 0., nRunsOO * 1.}});

      // Rename Y axis with Run numbers
      for (int iRun = 0; iRun < nRunsOO; ++iRun) {
        spectra.get<TH3>(HIST("hScaledFT0APerRunPerSetOfFlags"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
        spectra.get<TH3>(HIST("hScaledFT0CPerRunPerSetOfFlags"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
        spectra.get<TH3>(HIST("hScaledFT0MPerRunPerSetOfFlags"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);

        spectra.get<TH3>(HIST("hFT0APerRunPerSetOfFlags"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
        spectra.get<TH3>(HIST("hFT0CPerRunPerSetOfFlags"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
        spectra.get<TH3>(HIST("hFT0MPerRunPerSetOfFlags"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);

        spectra.get<TH2>(HIST("hIsFT0SignalComeFromCollPerRun"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);

        spectra.get<TH2>(HIST("hFT0AsignalWithoutBC"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
        spectra.get<TH2>(HIST("hFT0CsignalWithoutBC"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
        spectra.get<TH2>(HIST("hFT0MsignalWithoutBC"))->GetYaxis()->SetBinLabel(iRun + 1, runNumbersOO[iRun]);
      }

      // Rename Z axis with event selection flags
      for (int iFlag = 0; iFlag < nEvSelFlags; ++iFlag) {
        spectra.get<TH3>(HIST("hScaledFT0APerRunPerSetOfFlags"))->GetZaxis()->SetBinLabel(iFlag + 1, evSelFlags[iFlag]);
        spectra.get<TH3>(HIST("hScaledFT0CPerRunPerSetOfFlags"))->GetZaxis()->SetBinLabel(iFlag + 1, evSelFlags[iFlag]);
        spectra.get<TH3>(HIST("hScaledFT0MPerRunPerSetOfFlags"))->GetZaxis()->SetBinLabel(iFlag + 1, evSelFlags[iFlag]);

        spectra.get<TH3>(HIST("hFT0APerRunPerSetOfFlags"))->GetZaxis()->SetBinLabel(iFlag + 1, evSelFlags[iFlag]);
        spectra.get<TH3>(HIST("hFT0CPerRunPerSetOfFlags"))->GetZaxis()->SetBinLabel(iFlag + 1, evSelFlags[iFlag]);
        spectra.get<TH3>(HIST("hFT0MPerRunPerSetOfFlags"))->GetZaxis()->SetBinLabel(iFlag + 1, evSelFlags[iFlag]);
      }
    }

    if (doprocessLeadingAndAssociatedTracksTask) {

      spectra.add("hScaledFT0C_NleadTracks", "Total number of selected leading tracks vs scaled FT0C", kTH2F, {{scaledFT0C}, {1, 0.0, 1.}});
      spectra.add("hScaledFT0M_NleadTracks", "Total number of selected leading tracks vs scaled FT0M", kTH2F, {{scaledFT0M}, {1, 0.0, 1.}});

      spectra.add("hScaledFT0C_Correlation_LeadTrack_AssociatTracks", Form("Leading track #it{p}_{T} #in (%.2f, %.2f); Associated track #it{p}_{T} #in (%.2f, #it{p}_{T, lead. trk})", pTLeadTrack->at(0), pTLeadTrack->at(1), pTAssociatTrackMin.value), kTH2F, {{scaledFT0C}, {160, -1.28, 5.0, "#it{#varphi} (rad)"}});
      spectra.add("hScaledFT0M_Correlation_LeadTrack_AssociatTracks", Form("Leading track #it{p}_{T} #in (%.2f, %.2f); Associated track #it{p}_{T} #in (%.2f, #it{p}_{T, lead. trk})", pTLeadTrack->at(0), pTLeadTrack->at(1), pTAssociatTrackMin.value), kTH2F, {{scaledFT0M}, {160, -1.28, 5.0, "#it{#varphi} (rad)"}});
    }
  }

  // Fill histograms with raw or MC det. level data
  template <typename JCollision, typename Jets, typename JTracks>
  void fillHistograms(JCollision const& collision, Jets const& jets,
                      JTracks const& tracks, float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float rho = collision.rho();
    float scaledFT0C = getScaledFT0(collision.multFT0C(), meanFT0C);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), meanFT0A), scaledFT0C);

    auto dice = rand->Rndm();
    if (dice < fracSig)
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
      if (bSigEv && (trackPt > ptTTsig->at(0) && trackPt < ptTTsig->at(1))) {
        vPhiOfTT.push_back(trackPhi);
        spectra.fill(HIST("hTTSig_pT"), trackPt, weight);
        ++nTT;
      }

      if (!bSigEv && (trackPt > ptTTref->at(0) && trackPt < ptTTref->at(1))) {
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
      if (isJetWithHighPtConstituent(jet, tracks))
        continue;

      float jetPt = jet.pt();
      float jetArea = jet.area();
      float jetPtCorr = jetPt - rho * jetArea;

      spectra.fill(HIST("hJetPtEtaPhiRhoArea"), jetPt, jet.eta(), jet.phi(), rho, jetArea, weight);

      if (nTT > 0) {
        auto [dphi, bRecoilJet] = isRecoilJet(jet, phiTT);

        if (bSigEv) {
          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTSig"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTSig"), scaledFT0M, dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig"), jetArea, jetPt, rho, weight);

          if (phiTT > phiRestrTTSelection->at(0) && phiTT < phiRestrTTSelection->at(1)) {
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

          if (phiTT > phiRestrTTSelection->at(0) && phiTT < phiRestrTTSelection->at(1)) {
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
  void fillMCPHistograms(JCollision const& collision, Jets const& jets,
                         JParticles const& particles, float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float rho = collision.rho();
    float scaledFT0C = getScaledFT0(collision.multFT0C(), meanFT0CPartLevel);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), meanFT0APartLevel), scaledFT0C);

    auto dice = rand->Rndm();
    if (dice < fracSig)
      bSigEv = true;

    spectra.fill(HIST("hScaledFT0C_vertexZMC"), scaledFT0C, collision.posZ(), weight);
    spectra.fill(HIST("hScaledFT0M_vertexZMC"), scaledFT0M, collision.posZ(), weight);

    spectra.fill(HIST("hScaledFT0C_Rho_Part"), scaledFT0C, rho, weight);
    spectra.fill(HIST("hScaledFT0M_Rho_Part"), scaledFT0M, rho, weight);

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle)
        continue;

      float particlePt = particle.pt();
      float particlePhi = particle.phi();

      // Need charge and physical primary particles
      bool bParticleNeutral = (static_cast<int8_t>(pdgParticle->Charge()) == 0);
      if (bParticleNeutral || !particle.isPhysicalPrimary())
        continue;

      spectra.fill(HIST("hScaledFT0CPartPtEtaPhi"), scaledFT0C, particlePt, particle.eta(), particlePhi, weight);
      spectra.fill(HIST("hScaledFT0MPartPtEtaPhi"), scaledFT0M, particlePt, particle.eta(), particlePhi, weight);

      if (bSigEv && (particlePt > ptTTsig->at(0) && particlePt < ptTTsig->at(1))) {
        vPhiOfTT.push_back(particlePhi);
        ++nTT;
      }

      if (!bSigEv && (particlePt > ptTTref->at(0) && particlePt < ptTTref->at(1))) {
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

        auto [dphi, bRecoilJet] = isRecoilJet(jet, phiTT);

        if (bSigEv) {
          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig_Part"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig_Part"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTSig_Part"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTSig_Part"), scaledFT0M, dphi, jetPt, weight);

          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig_Part"), jetArea, jetPt, rho, weight);

          if (phiTT > phiRestrTTSelection->at(0) && phiTT < phiRestrTTSelection->at(1)) {
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

          if (phiTT > phiRestrTTSelection->at(0) && phiTT < phiRestrTTSelection->at(1)) {
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

  template <typename JTracksTable, typename JetsBase, typename JetsTag>
  void fillMatchedHistograms(JTracksTable const& tracks,
                             JetsBase const& jetsBase, JetsTag const& jetsTag,
                             float weight = 1.)
  {
    std::vector<double> vPhiOfTT;
    double phiTTSig = 0.;

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      if (track.pt() > ptTTsig->at(0) && track.pt() < ptTTsig->at(1)) {
        vPhiOfTT.push_back(track.phi());
      }
    }

    bool bIsThereTTSig = vPhiOfTT.size() > 0;

    if (bIsThereTTSig)
      phiTTSig = getPhiTT(vPhiOfTT);

    for (const auto& jetBase : jetsBase) {
      bool bIsBaseJetRecoil =
        get<1>(isRecoilJet(jetBase, phiTTSig)) && bIsThereTTSig;
      dataForUnfolding(jetBase, jetsTag, bIsBaseJetRecoil, tracks, weight);
    }
  }

  template <typename JCollision>
  void fillMultiplicityHistogramsOO(JCollision const& collision,
                                    float weight = 1.)
  {
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multFT0M = collision.multFT0M();
    float scaledFT0A = getScaledFT0(multFT0A, meanFT0A);
    float scaledFT0C = getScaledFT0(multFT0C, meanFT0C);
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
  void fillMultiplicityHistogramsPartLevelMC(JCollisionMC const& collision,
                                             float weight = 1.)
  {
    spectra.fill(HIST("hMultFT0APartLevel"), collision.multFT0A(), weight);
    spectra.fill(HIST("hMultFT0CPartLevel"), collision.multFT0C(), weight);
    spectra.fill(HIST("hMultFT0MPartLevel"), collision.multFT0A() + collision.multFT0C(), weight);

    auto scaledFT0A = getScaledFT0(collision.multFT0A(), meanFT0APartLevel);
    auto scaledFT0C = getScaledFT0(collision.multFT0C(), meanFT0CPartLevel);
    spectra.fill(HIST("hScaleMultFT0APartLevel"), scaledFT0A, weight);
    spectra.fill(HIST("hScaleMultFT0CPartLevel"), scaledFT0C, weight);
    spectra.fill(HIST("hScaleMultFT0MPartLevel"), getScaledFT0M(scaledFT0A, scaledFT0C), weight);
  }

  template <typename BC, typename Collision, typename ZDC>
  void fillMultiplicityQA(Collision const& collision, BC const&,
                          ZDC const&, float weight = 1.)
  {
    int runNumber = collision.multRunNumber();
    int fillNumber = getBinNumberOnYaxisForGivenRun(spectra.get<TH3>(HIST("hScaledFT0CPerRunPerSetOfFlags")), runNumber) - 0.5; // Same for FT0M distrib.

    // FT0 Signal
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multFT0M = collision.multFT0M();
    float scaledFT0A = getScaledFT0(multFT0A, meanFT0A);
    float scaledFT0C = getScaledFT0(multFT0C, meanFT0C);
    float scaledFT0M = getScaledFT0M(scaledFT0A, scaledFT0C);

    // Event with flag Sel8
    spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 0.5, weight);
    spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 0.5, weight);
    spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 0.5, weight);

    spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 0.5, weight);
    spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 0.5, weight);
    spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 0.5, weight);

    bool isGoodZvtxFT0vsPV = collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV);
    if (isGoodZvtxFT0vsPV) {
      spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 1.5, weight);
      spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 1.5, weight);
      spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 1.5, weight);

      spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 1.5, weight);
      spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 1.5, weight);
      spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 1.5, weight);
    }

    bool isNoSameBunchPileup = collision.selection_bit(aod::evsel::kNoSameBunchPileup);
    if (isNoSameBunchPileup) {
      spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 2.5, weight);
      spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 2.5, weight);
      spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 2.5, weight);

      spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 2.5, weight);
      spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 2.5, weight);
      spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 2.5, weight);
    }

    if (isGoodZvtxFT0vsPV && isNoSameBunchPileup) {
      spectra.fill(HIST("hFT0APerRunPerSetOfFlags"), multFT0A, fillNumber, 3.5, weight);
      spectra.fill(HIST("hFT0CPerRunPerSetOfFlags"), multFT0C, fillNumber, 3.5, weight);
      spectra.fill(HIST("hFT0MPerRunPerSetOfFlags"), multFT0M, fillNumber, 3.5, weight);

      spectra.fill(HIST("hScaledFT0APerRunPerSetOfFlags"), scaledFT0A, fillNumber, 3.5, weight);
      spectra.fill(HIST("hScaledFT0CPerRunPerSetOfFlags"), scaledFT0C, fillNumber, 3.5, weight);
      spectra.fill(HIST("hScaledFT0MPerRunPerSetOfFlags"), scaledFT0M, fillNumber, 3.5, weight);
    } else {
      return;
    }

    // Investigate other EA variables
    //____________________________________________________________________________________

    // Multiplicity equalized for the vertex position with FT0 detector
    float multZeqFT0A = collision.multZeqFT0A();
    float multZeqFT0C = collision.multZeqFT0C();
    float multZeqFT0M = multZeqFT0A + multZeqFT0C;
    float scaledZeqFT0A = getScaledFT0(multZeqFT0A, meanZeqFT0A);
    float scaledZeqFT0C = getScaledFT0(multZeqFT0C, meanZeqFT0C);
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
      spectra.fill(HIST("hFT0AsignalWithoutBC"), multFT0A, fillNumber, weight);
      spectra.fill(HIST("hFT0CsignalWithoutBC"), multFT0C, fillNumber, weight);
      spectra.fill(HIST("hFT0MsignalWithoutBC"), multFT0M, fillNumber, weight);
    }

    if (collision.foundBCId() > 0)
      spectra.fill(HIST("hIsFT0SignalComeFromCollPerRun"), 2.5, fillNumber, weight);
    else
      spectra.fill(HIST("hIsFT0SignalComeFromCollPerRun"), 3.5, fillNumber, weight);
  }

  template <typename JCollision, typename JTracks>
  void fillLeadingAndAssociatedTracksTask(JCollision const& collision, JTracks const& tracks, float weight = 1.)
  {
    std::vector<double> vPhiOfLeadingTracks;
    std::vector<double> vPtOfLeadingTracks;
    std::vector<double> vPhiOfAssociatedTracks;

    float scaledFT0C = getScaledFT0(collision.multFT0C(), meanFT0C);
    float scaledFT0M = getScaledFT0M(getScaledFT0(collision.multFT0A(), meanFT0A), scaledFT0C);

    // Search for leading tracks
    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      float trackPt = track.pt();

      if (trackPt > pTLeadTrack->at(0) && trackPt < pTLeadTrack->at(1)) {
        vPhiOfLeadingTracks.push_back(track.phi());
        vPtOfLeadingTracks.push_back(trackPt);
      }
    }

    int nLeadingTracks = vPhiOfLeadingTracks.size();

    if (nLeadingTracks > 0) {
      auto indexLeadTrack = rand->Integer(nLeadingTracks);

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
        if (trackPt > pTAssociatTrackMin && trackPt < pTLeadingTrack) {
          double dphi = RecoDecay::constrainAngle(phiLeadingTrack - trackPhi, -1.3);
          spectra.fill(HIST("hScaledFT0C_Correlation_LeadTrack_AssociatTracks"), scaledFT0C, dphi, weight);
          spectra.fill(HIST("hScaledFT0M_Correlation_LeadTrack_AssociatTracks"), scaledFT0M, dphi, weight);
        }
      }
    }
  }

  //------------------------------------------------------------------------------
  // Process functions
  void processData(FilteredColl const& collision, FilteredTracks const& tracks,
                   FilteredJets const& jets)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);

    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5); // number of events selected for analysis

    fillHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(RecoilJets, processData, "process raw data", true);

  void processMCDetLevel(FilteredColl const& collision,
                         FilteredTracks const& tracks,
                         FilteredJetsDetLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    if (skipMBGapEvent(collision)) {
      spectra.fill(HIST("hEventSelectionCount"), 2.5);
      return;
    }

    spectra.fill(HIST("hEventSelectionCount"), 5.5); // number of events selected for analysis
    fillHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevel, "process MC detector level data (no weight)", false);

  void processMCDetLevelWeighted(FilteredCollDetLevelGetWeight const& collision,
                                 aod::JetMcCollisions const&,
                                 FilteredTracks const& tracks,
                                 FilteredJetsDetLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);
    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    if (skipMBGapEvent(collision)) {
      spectra.fill(HIST("hEventSelectionCount"), 2.5);
      return;
    }

    if (collision.isOutlier()) {
      spectra.fill(HIST("hEventSelectionCount"), 3.5);
      return;
    }

    if (!collision.has_mcCollision()) {
      spectra.fill(HIST("hEventSelectionCount"), 4.5);
      return;
    }

    spectra.fill(HIST("hEventSelectionCount"), 5.5); // number of events selected for analysis
    auto weight = collision.mcCollision().weight();
    fillHistograms(collision, jets, tracks, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevelWeighted, "process MC detector level data (weighted)", false);

  void processMCPartLevel(FilteredCollPartLevel const& collision,
                          FilteredParticles const& particles,
                          FilteredJetsPartLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);

    if (skipMBGapEvent(collision)) {
      spectra.fill(HIST("hEventSelectionCountPartLevel"), 1.5);
      return;
    }

    spectra.fill(HIST("hEventSelectionCountPartLevel"), 3.5); // number of events selected for analysis
    fillMCPHistograms(collision, jets, particles);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevel, "process MC particle level data (no weight)", false);

  void processMCPartLevelWeighted(FilteredCollPartLevel const& collision,
                                  FilteredParticles const& particles,
                                  FilteredJetsPartLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);

    if (skipMBGapEvent(collision)) {
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
    fillMCPHistograms(collision, jets, particles, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevelWeighted, "process MC particle level data (weighted)", false);

  void processJetsMatched(FilteredCollDetLevelGetWeight const& collision,
                          aod::JetMcCollisions const&,
                          FilteredTracks const& tracks,
                          FilteredMatchedJetsDetLevel const& mcdjets,
                          FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision) || skipMBGapEvent(collision) || collision.isOutlier())
      return;

    auto mcpjetsPerMCCollision = mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());

    fillMatchedHistograms(tracks, mcpjetsPerMCCollision, mcdjets);
  }
  PROCESS_SWITCH(RecoilJets, processJetsMatched, "process matching of MC jets (no weight)", false);

  void processJetsMatchedWeighted(FilteredCollDetLevelGetWeight const& collision,
                                  aod::JetMcCollisions const&,
                                  FilteredTracks const& tracks,
                                  FilteredMatchedJetsDetLevel const& mcdjets,
                                  FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision) || skipMBGapEvent(collision) || collision.isOutlier())
      return;

    auto mcpjetsPerMCCollision = mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    auto weight = collision.mcCollision().weight();

    fillMatchedHistograms(tracks, mcpjetsPerMCCollision, mcdjets, weight);
  }
  PROCESS_SWITCH(RecoilJets, processJetsMatchedWeighted, "process matching of MC jets (weighted)", false);

  void processMultiplicityOO(FilteredEventMultiplicity const& collision)
  {
    if (skipEvent(collision))
      return;

    fillMultiplicityHistogramsOO(collision);
  }
  PROCESS_SWITCH(RecoilJets, processMultiplicityOO, "process multiplicity for OO collisions and MC detector level (no weight)", false);

  void processMultiplicityMCDetLevelWeightedOO(FilteredEventMultiplicityDetLevelGetWeight const& collision)
  {
    if (skipEvent(collision) || collision.isOutlier() || !collision.has_mcCollision())
      return;

    auto weight = collision.mcCollision().weight();
    fillMultiplicityHistogramsOO(collision, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMultiplicityMCDetLevelWeightedOO, "process multiplicity for MC detector level OO collisions (weighted)", false);

  void processMultiplicityPartLevelMC(FilteredEventMultiplicityPartLevel const& collision)
  {
    if (skipMBGapEvent(collision))
      return;

    fillMultiplicityHistogramsPartLevelMC(collision);
  }
  PROCESS_SWITCH(RecoilJets, processMultiplicityPartLevelMC, "process multiplicity for MC particle level events (no weight)", false);

  void processMultiplicityPartLevelMCWeighted(FilteredEventMultiplicityPartLevel const& collision)
  {
    if (skipMBGapEvent(collision) || collision.isOutlier())
      return;

    auto weight = collision.weight();
    fillMultiplicityHistogramsPartLevelMC(collision, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMultiplicityPartLevelMCWeighted, "process multiplicity for MC particle level events (weighted)", false);

  void processMultiplicityQA(ColEvSelEA const& collision,
                             BCsRun3 const& BCs,
                             aod::Zdcs const& ZDCs)
  {
    // Base flag for event selection
    if (!collision.sel8())
      return;

    fillMultiplicityQA(collision, BCs, ZDCs);
  }
  PROCESS_SWITCH(RecoilJets, processMultiplicityQA, "process function for EA QA purposes", false);

  void processLeadingAndAssociatedTracksTask(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                             soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (skipEvent(collision))
      return;
    fillLeadingAndAssociatedTracksTask(collision, tracks);
  }
  PROCESS_SWITCH(RecoilJets, processLeadingAndAssociatedTracksTask, "process function for correlation between leading and associated tracks", false);

  //------------------------------------------------------------------------------
  // Auxiliary functions
  template <typename Collision>
  bool skipEvent(const Collision& coll)
  {
    /// \brief: trigger cut is needed for pp data
    return !jetderiveddatautilities::selectCollision(coll, eventSelectionBits) || !jetderiveddatautilities::selectTrigger(coll, triggerMaskBits);
  }

  template <typename Collision>
  bool skipMBGapEvent(const Collision& coll)
  {
    return skipMBGapEvents && coll.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap;
  }

  template <typename Track>
  bool skipTrack(const Track& track)
  {
    return !jetderiveddatautilities::selectTrack(track, trackSelection);
  }

  template <typename Jet>
  std::tuple<double, bool> isRecoilJet(const Jet& jet, double phiTT)
  {
    double dphi = std::fabs(RecoDecay::constrainAngle(jet.phi() - phiTT, -constants::math::PI));
    return {dphi, (constants::math::PI - recoilRegion) < dphi};
  }

  double getPhiTT(const std::vector<double>& vPhiOfTT)
  {
    auto iTrig = rand->Integer(vPhiOfTT.size());
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

  template <typename Jet, typename Tracks>
  bool isJetWithHighPtConstituent(Jet const& jet, Tracks const&)
  {
    bool bIsJetWithHighPtConstituent = false;
    for (const auto& jetConstituent : jet.template tracks_as<Tracks>()) {
      if (jetConstituent.pt() > maxJetConstituentPt) {
        bIsJetWithHighPtConstituent = true;
        break;
      }
    }
    return bIsJetWithHighPtConstituent;
  }

  template <typename histo>
  int getBinNumberOnYaxisForGivenRun(std::shared_ptr<histo> histogram, int runNumber)
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

  template <typename PartJet, typename DetJet, typename TracksTable>
  void dataForUnfolding(PartJet const& partJet, DetJet const& detJets,
                        bool bIsBaseJetRecoil, TracksTable const& tracks, float weight = 1.)
  {

    float partJetPt = partJet.pt();
    bool bIsThereMatchedJet = partJet.has_matchedJetGeo();

    if (bIsThereMatchedJet) {
      const auto& jetsMatched = partJet.template matchedJetGeo_as<std::decay_t<DetJet>>();

      for (const auto& jetMatched : jetsMatched) {

        // skip matches where detector level jets have a constituent with pT above specified cut
        bool skipMatchedDetJet = isJetWithHighPtConstituent(jetMatched, tracks);

        if (skipMatchedDetJet) {
          // Miss jets
          spectra.fill(HIST("hMissedJets_pT"), partJetPt, weight);
          if (bIsBaseJetRecoil)
            spectra.fill(HIST("hMissedJets_pT_RecoilJets"), partJetPt, weight);
        } else {
          float detJetPt = jetMatched.pt();

          spectra.fill(HIST("hNumberMatchedJetsPerOneBaseJet"), jetsMatched.size(), detJetPt, weight);
          spectra.fill(HIST("hJetPt_DetLevel_vs_PartLevel"), detJetPt, partJetPt, weight);
          spectra.fill(HIST("hJetPt_resolution"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
          spectra.fill(HIST("hJetPhi_resolution"), partJet.phi() - jetMatched.phi(), partJetPt, weight);

          if (bIsBaseJetRecoil) {
            spectra.fill(HIST("hJetPt_DetLevel_vs_PartLevel_RecoilJets"), detJetPt, partJetPt, weight);
            spectra.fill(HIST("hJetPt_resolution_RecoilJets"), (partJetPt - detJetPt) / partJetPt, partJetPt, weight);
            spectra.fill(HIST("hJetPhi_resolution_RecoilJets"), partJet.phi() - jetMatched.phi(), partJetPt, weight);
          }
        }
      }
    } else {
      // Miss jets
      spectra.fill(HIST("hMissedJets_pT"), partJetPt, weight);
      if (bIsBaseJetRecoil)
        spectra.fill(HIST("hMissedJets_pT_RecoilJets"), partJetPt, weight);
    }

    // Fake jets
    for (const auto& detJet : detJets) {
      if (isJetWithHighPtConstituent(detJet, tracks))
        continue;

      bIsThereMatchedJet = detJet.has_matchedJetGeo();
      if (!bIsThereMatchedJet) {
        spectra.fill(HIST("hFakeJets_pT"), detJet.pt(), weight);
        if (bIsBaseJetRecoil)
          spectra.fill(HIST("hFakeJets_pT_RecoilJets"), detJet.pt(), weight);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<RecoilJets>(cfgc)};
}
