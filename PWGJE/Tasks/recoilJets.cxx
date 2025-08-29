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
using FilteredColl =
  soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator;
using FilteredCollPartLevel =
  soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos,
                          aod::JMcCollisionOutliers>>::iterator;
using FilteredCollDetLevelGetWeight =
  soa::Filtered<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos,
                          aod::JCollisionOutliers>>::iterator;
using FilteredEventMultiplicity =
  soa::Filtered<soa::Join<aod::JetCollisions, aod::ZDCMults>>::iterator;

using FilteredJets =
  soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using FilteredJetsDetLevel =
  soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets,
                          aod::ChargedMCDetectorLevelJetConstituents>>;
using FilteredJetsPartLevel =
  soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets,
                          aod::ChargedMCParticleLevelJetConstituents>>;

using FilteredMatchedJetsDetLevel = soa::Filtered<soa::Join<
  aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents,
  aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using FilteredMatchedJetsPartLevel = soa::Filtered<soa::Join<
  aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents,
  aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

using FilteredTracks = soa::Filtered<aod::JetTracks>;
using FilteredParticles = soa::Filtered<aod::JetParticles>;

struct RecoilJets {

  // List of configurable parameters
  Configurable<std::string> evSel{"evSel", "sel8", "Choose event selection"};
  Configurable<std::string> trkSel{"trkSel", "globalTracks",
                                   "Set track selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10., "Accepted z-vertex range"};
  Configurable<float> fracSig{"fracSig", 0.9,
                              "Fraction of events to use for signal TT"};

  Configurable<float> trkPtMin{"trkPtMin", 0.15,
                               "Minimum pT of acceptanced tracks"};
  Configurable<float> trkPtMax{"trkPtMax", 100.,
                               "Maximum pT of acceptanced tracks"};

  Configurable<float> trkEtaCut{"trkEtaCut", 0.9, "Eta acceptance of TPC"};
  Configurable<float> jetR{"jetR", 0.4, "Jet cone radius"};

  Configurable<std::string> triggerMasks{"triggerMasks", "",
                                         "Relevant trigger masks: fTrackLowPt,fTrackHighPt"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false,
                                     "flag to choose to reject min. bias gap events; jet-level rejection "
                                     "applied at the jet finder level, here rejection is applied for "
                                     "collision and track process functions"};

  Configurable<float> meanFT0A{"meanFT0A", -1.0, "Mean value of FT0A"};
  Configurable<float> meanFT0C{"meanFT0C", -1.0, "Mean value of FT0C"};

  // List of configurable parameters for MC
  Configurable<float> pTHatExponent{"pTHatExponent", 4.0,
                                    "Exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMax{"pTHatMax", 999.0,
                               "Maximum fraction of hard scattering for jet acceptance in MC"};

  // Parameters for recoil jet selection
  Configurable<float> ptTTrefMin{"ptTTrefMin", 5.,
                                 "Minimum pT of reference TT"};
  Configurable<float> ptTTrefMax{"ptTTrefMax", 7.,
                                 "Maximum pT of reference TT"};
  Configurable<float> ptTTsigMin{"ptTTsigMin", 20., "Minimum pT of signal TT"};
  Configurable<float> ptTTsigMax{"ptTTsigMax", 50., "Maximum pT of signal TT"};
  Configurable<float> recoilRegion{"recoilRegion", 0.6,
                                   "Width of recoil acceptance"};

  Configurable<float> maxJetConstituentPt{"maxJetConstituentPt", 100.,
                                          "Remove jets with constituent above this pt cut"};

  // List of configurable parameters for histograms
  Configurable<uint16_t> histJetPt{"histJetPt", 100,
                                   "Maximum value of jet pT shown in histograms"};

  Configurable<uint16_t> histMultBins{"histMultBins", 1000,
                                      "Number of bins for scaled FT0M multiplicity"};

  // Axes specification
  AxisSpec pT{histJetPt, 0.0, histJetPt * 1.0, "#it{p}_{T} (GeV/#it{c})"};
  AxisSpec jetPTcorr{histJetPt + 20, -20., histJetPt * 1.0,
                     "#it{p}_{T, jet}^{ch, corr} (GeV/#it{c})"};
  AxisSpec phiAngle{40, 0.0, constants::math::TwoPI, "#it{#varphi} (rad)"};
  AxisSpec deltaPhiAngle{52, 0.0, constants::math::PI,
                         "#Delta#it{#varphi} (rad)"};
  AxisSpec pseudorap{40, -1., 1., "#it{#eta}"};
  AxisSpec pseudorapJets{20, -0.5, 0.5, "#it{#eta}_{jet}"};
  AxisSpec jetArea{50, 0.0, 5., "Area_{jet}"};
  AxisSpec rhoArea{60, 0.0, 60., "#it{#rho} #times Area_{jet}"};
  AxisSpec rho{50, 0.0, 50., "#it{#rho}"};
  AxisSpec scaledFT0C{histMultBins, 0.0, 20., "FT0C / #LT FT0C #GT"};
  AxisSpec scaledFT0M{histMultBins, 0.0, 20., "FT0M^{*}"};
  ConfigurableAxis multFT0CPercentile{"multFT0CPercentile", {VARIABLE_WIDTH, 0, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.4, 1.8, 2.4, 3.6, 5.0}, "Percentiles of scaled FT0C: 100-90%, 90-80%, 80-70%, 70-60%, 60-50%, 50-40%, 40-30%, 30-20%, 20-10%, 10-1%, 1-0.1%"}; // to adjust the boarders
  ConfigurableAxis multFT0MPercentile{"multFT0MPercentile", {VARIABLE_WIDTH, 0, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.4, 1.8, 2.4, 3.6, 5.0}, "Percentiles of scaled FT0M: 100-90%, 90-80%, 80-70%, 70-60%, 60-50%, 50-40%, 40-30%, 30-20%, 20-10%, 10-1%, 1-0.1%"};

  Preslice<FilteredMatchedJetsPartLevel> partJetsPerCollision = aod::jet::mcCollisionId;

  TRandom3* rand = new TRandom3(0);

  // Declare filter on collision Z vertex
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter collisionFilterMC = nabs(aod::jmccollision::posZ) < vertexZCut;

  // Declare filters on accepted tracks and MC particles (settings for jet reco
  // are provided in the jet finder wagon)
  Filter trackFilter = aod::jtrack::pt > trkPtMin&& aod::jtrack::pt <
                       trkPtMax&& nabs(aod::jtrack::eta) < trkEtaCut;
  Filter partFilter = nabs(aod::jmcparticle::eta) < trkEtaCut;

  // Declare filter on jets
  Filter jetRadiusFilter = aod::jet::r == nround(jetR.node() * 100.);

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  Service<o2::framework::O2DatabasePDG> pdg;

  void init(InitContext const&)
  {
    std::string evSelToString = static_cast<std::string>(evSel);
    std::string trkSelToString = static_cast<std::string>(trkSel);

    eventSelectionBits =
      jetderiveddatautilities::initialiseEventSelectionBits(evSelToString);
    trackSelection =
      jetderiveddatautilities::initialiseTrackSelection(trkSelToString);
    triggerMaskBits =
      jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);

    // List of raw and MC det. distributions
    if (doprocessData || doprocessMCDetLevel || doprocessMCDetLevelWeighted) {
      spectra.add("hEventSelectionCount", "Count # of events in the analysis", kTH1F, {{3, 0.0, 3.}});
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of events");
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of events after sel. %s", evSelToString.data()));
      spectra.get<TH1>(HIST("hEventSelectionCount"))->GetXaxis()->SetBinLabel(3, "# of events w. outlier");

      spectra.add("vertexZ", "Z vertex of collisions", kTH1F, {{60, -12., 12.}});
      spectra.add("hHasAssocMcCollision", "Has det. level coll. associat. MC coll.", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hHasAssocMcCollision"))->GetXaxis()->SetBinLabel(1, "Yes");
      spectra.get<TH1>(HIST("hHasAssocMcCollision"))->GetXaxis()->SetBinLabel(2, "No");

      spectra.add("hTrackSelectionCount", "Count # of tracks in the analysis", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hTrackSelectionCount"))->GetXaxis()->SetBinLabel(1, "Total # of tracks");
      spectra.get<TH1>(HIST("hTrackSelectionCount"))->GetXaxis()->SetBinLabel(2, Form("# of tracks after sel. %s", trkSelToString.data()));

      spectra.add("hTrackPtEtaPhi", "Charact. of tracks", kTH3F, {pT, pseudorap, phiAngle});
      spectra.add("hTTSig_pT", "pT spectrum of all found TT_{Sig} cand.", kTH1F, {{40, 10., 50.}}); // needed to distinguish merged data from diff. wagons

      spectra.add("hScaledFT0C_vs_Ntrig", "Total number of selected triggers per class vs scaled FT0C", kTH2F, {{multFT0CPercentile}, {2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hScaledFT0C_vs_Ntrig"))->GetYaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH1>(HIST("hScaledFT0C_vs_Ntrig"))->GetYaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hScaledFT0M_vs_Ntrig", "Total number of selected triggers per class vs scaled FT0M", kTH2F, {{multFT0MPercentile}, {2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hScaledFT0M_vs_Ntrig"))->GetYaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH1>(HIST("hScaledFT0M_vs_Ntrig"))->GetYaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hScaledFT0C_vs_TTRef_per_event", "Number of TT_{Ref} per event vs scaled FT0C", kTH2F, {{multFT0CPercentile}, {15, 0.5, 15.5}});
      spectra.add("hScaledFT0M_vs_TTRef_per_event", "Number of TT_{Ref} per event vs scaled FT0M", kTH2F, {{multFT0MPercentile}, {15, 0.5, 15.5}});

      spectra.add("hScaledFT0C_vs_TTSig_per_event", "Number of TT_{Sig} per event vs scaled FT0C", kTH2F, {{multFT0CPercentile}, {10, 0.5, 10.5}});
      spectra.add("hScaledFT0M_vs_TTSig_per_event", "Number of TT_{Sig} per event vs scaled FT0M", kTH2F, {{multFT0MPercentile}, {10, 0.5, 10.5}});

      spectra.add("hJetPtEtaPhiRhoArea", "Charact. of inclusive jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0CPercentile, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0MPercentile, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0CPercentile, deltaPhiAngle, jetPTcorr});
      spectra.add("hScaledFT0M_DPhi_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0MPercentile, deltaPhiAngle, jetPTcorr});

      spectra.add("hScaledFT0C_DPhi_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0CPercentile, deltaPhiAngle, pT});
      spectra.add("hScaledFT0M_DPhi_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0MPercentile, deltaPhiAngle, pT});

      spectra.add("hScaledFT0C_DPhi_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0C & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0CPercentile, deltaPhiAngle, pT});
      spectra.add("hScaledFT0M_DPhi_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0M & #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH3F, {multFT0MPercentile, deltaPhiAngle, pT});

      spectra.add("hScaledFT0C_Recoil_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {multFT0CPercentile, jetPTcorr});
      spectra.add("hScaledFT0M_Recoil_JetPt_Corr_TTRef", "Events w. TT_{Ref}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {multFT0MPercentile, jetPTcorr});

      spectra.add("hScaledFT0C_Recoil_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {multFT0CPercentile, jetPTcorr});
      spectra.add("hScaledFT0M_Recoil_JetPt_Corr_TTSig", "Events w. TT_{Sig}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {multFT0MPercentile, jetPTcorr});

      spectra.add("hScaledFT0C_Recoil_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {multFT0CPercentile, pT});
      spectra.add("hScaledFT0M_Recoil_JetPt_TTRef", "Events w. TT_{Ref}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {multFT0MPercentile, pT});

      spectra.add("hScaledFT0C_Recoil_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0C & #it{p}_{T} of recoil jets", kTH2F, {multFT0CPercentile, pT});
      spectra.add("hScaledFT0M_Recoil_JetPt_TTSig", "Events w. TT_{Sig}: scaled FT0M & #it{p}_{T} of recoil jets", kTH2F, {multFT0MPercentile, pT});

      spectra.add("hJetArea_JetPt_Rho_TTRef", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
      spectra.add("hJetArea_JetPt_Rho_TTSig", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});

      spectra.add("hScaledFT0C_Rho_TTRef", "Events w. TT_{Ref}: scaled FT0C & #rho", kTH2F, {multFT0CPercentile, rho});
      spectra.add("hScaledFT0M_Rho_TTRef", "Events w. TT_{Ref}: scaled FT0M & #rho", kTH2F, {multFT0MPercentile, rho});

      spectra.add("hScaledFT0C_Rho_TTSig", "Events w. TT_{Sig}: scaled FT0C & #rho", kTH2F, {multFT0CPercentile, rho});
      spectra.add("hScaledFT0M_Rho_TTSig", "Events w. TT_{Sig}: scaled FT0M & #rho", kTH2F, {multFT0MPercentile, rho});

      spectra.add("hScaledFT0C_TTRef", "Events w. TT_{Ref}: scaled FT0C", kTH1F, {scaledFT0C});
      spectra.add("hScaledFT0M_TTRef", "Events w. TT_{Ref}: scaled FT0M", kTH1F, {scaledFT0M});

      spectra.add("hScaledFT0C_TTSig", "Events w. TT_{Sig}: scaled FT0C", kTH1F, {scaledFT0C});
      spectra.add("hScaledFT0M_TTSig", "Events w. TT_{Sig}: scaled FT0M", kTH1F, {scaledFT0M});
    }

    // List of MC particle level distributions
    if (doprocessMCPartLevel || doprocessMCPartLevelWeighted) {
      spectra.add("vertexZMC", "Z vertex of jmccollision", kTH1F, {{60, -12., 12.}});
      spectra.add("ptHat", "Distribution of pT hat", kTH1F, {{500, 0.0, 100.}});

      spectra.add("hEventSelectionCountPartLevel", "Count # of events in the part. level analysis", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(1, "Total # of events");
      spectra.get<TH1>(HIST("hEventSelectionCountPartLevel"))->GetXaxis()->SetBinLabel(2, "# of events w. outlier");

      spectra.add("hCountNumberOutliersFrameWork", "Count # of outlier events based on flag from JE fw", kTH1F, {{1, 0.0, 1.}});
      spectra.get<TH1>(HIST("hCountNumberOutliersFrameWork"))->GetXaxis()->SetBinLabel(1, "Outlier flag true");

      spectra.add("hPartPtEtaPhi", "Charact. of particles", kTH3F, {pT, pseudorap, phiAngle});
      spectra.add("hNtrig_Part", "Total number of selected triggers per class", kTH1F, {{2, 0.0, 2.}});
      spectra.get<TH1>(HIST("hNtrig_Part"))->GetXaxis()->SetBinLabel(1, "TT_{ref}");
      spectra.get<TH1>(HIST("hNtrig_Part"))->GetXaxis()->SetBinLabel(2, "TT_{sig}");

      spectra.add("hTTRef_per_event_Part", "Number of TT_{Ref} per event", kTH1F, {{15, 0.5, 15.5}});
      spectra.add("hTTSig_per_event_Part", "Number of TT_{Sig} per event", kTH1F, {{10, 0.5, 10.5}});

      spectra.add("hJetPtEtaPhiRhoArea_Part", "Charact. of inclusive part. level jets", kTHnSparseF, {pT, pseudorapJets, phiAngle, rho, jetArea});

      spectra.add("hDPhi_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, jetPTcorr});
      spectra.add("hDPhi_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, jetPTcorr});
      spectra.add("hDPhi_JetPt_TTRef_Part", "Events w. TT_{Ref}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, pT});
      spectra.add("hDPhi_JetPt_TTSig_Part", "Events w. TT_{Sig}: #Delta#varphi & #it{p}_{T, jet}^{ch}", kTH2F, {deltaPhiAngle, pT});

      spectra.add("hRecoil_JetPt_Corr_TTRef_Part", "Events w. TT_{Ref}: #it{p}_{T} of recoil jets", kTH1F, {jetPTcorr});
      spectra.add("hRecoil_JetPt_Corr_TTSig_Part", "Events w. TT_{Sig}: #it{p}_{T} of recoil jets", kTH1F, {jetPTcorr});
      spectra.add("hRecoil_JetPt_TTRef_Part", "Events w. TT_{Ref}: #it{p}_{T} of recoil jets", kTH1F, {pT});
      spectra.add("hRecoil_JetPt_TTSig_Part", "Events w. TT_{Sig}: #it{p}_{T} of recoil jets", kTH1F, {pT});

      spectra.add("hJetArea_JetPt_Rho_TTRef_Part", "Events w. TT_{Ref}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});
      spectra.add("hJetArea_JetPt_Rho_TTSig_Part", "Events w. TT_{Sig}: A_{jet} & jet pT & #rho", kTH3F, {jetArea, pT, rho});

      spectra.add("hDiffInOutlierRemove", "Difference between pT hat from code and fw", kTH1F, {{502, -0.2, 50.}});
    }

    // Jet matching: part. vs. det.
    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      spectra.add("hJetPt_DetLevel_vs_PartLevel", "Correlation jet pT at det. vs. part. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}});
      // spectra.add("hJetPt_Corr_PartLevel_vs_DetLevel", "Correlation jet pT at
      // part. vs. det. levels", kTH2F, {jetPTcorr, jetPTcorr});
      spectra.add("hJetPt_DetLevel_vs_PartLevel_RecoilJets", "Correlation recoil jet pT at part. vs. det. levels", kTH2F, {{200, 0.0, 200.}, {200, 0.0, 200.}});
      // spectra.add("hJetPt_Corr_PartLevel_vs_DetLevel_RecoilJets",
      // "Correlation recoil jet pT at part. vs. det. levels", kTH2F,
      // {jetPTcorr, jetPTcorr});

      spectra.add("hMissedJets_pT", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hMissedJets_Corr_pT", "Part. level jets w/o matched pair",
      // kTH1F, {jetPTcorr});
      spectra.add("hMissedJets_pT_RecoilJets", "Part. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hMissedJets_Corr_pT_RecoilJets", "Part. level jets w/o
      // matched pair", kTH1F, {jetPTcorr});
      spectra.add("hFakeJets_pT", "Det. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hFakeJets_Corr_pT", "Det. level jets w/o matched pair",
      // kTH1F, {jetPTcorr});
      spectra.add("hFakeJets_pT_RecoilJets", "Det. level jets w/o matched pair", kTH1F, {{200, 0.0, 200.}});
      // spectra.add("hFakeJets_Corr_pT_RecoilJets", "Det. level jets w/o
      // matched pair", kTH1F, {jetPTcorr});

      spectra.add("hJetPt_resolution", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT});
      spectra.add("hJetPt_resolution_RecoilJets", "Jet p_{T} relative resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{100, -5., 5.}, pT});

      spectra.add("hJetPhi_resolution", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT});
      spectra.add("hJetPhi_resolution_RecoilJets", "#varphi resolution as a func. of jet #it{p}_{T, part}", kTH2F, {{40, -1., 1.}, pT});

      spectra.add("hNumberMatchedJetsPerOneBaseJet", "# of tagged jets per 1 base jet vs. jet pT", kTH2F, {{10, 0.5, 10.5}, {100, 0.0, 100.}});
    }

    if (doprocessMultiplicity) {
      spectra.add("hMultFT0A", "Mult. signal from FTOA", kTH1F, {{2000, 0.0, 40000.}});
      spectra.add("hMultFT0C", "Mult. signal from FTOC", kTH1F, {{2000, 0.0, 40000.}});
      spectra.add("hMultFT0M", "Total mult. signal from FT0A & FTOC", kTH1F, {{3000, 0.0, 60000.}});

      spectra.add("hScaleMultFT0A", "Scaled mult. signal from FTOA", kTH1F, {{200, 0.0, 20., "FT0A / #LT FT0A #GT"}});
      spectra.add("hScaleMultFT0C", "Scaled mult. signal from FTOC", kTH1F, {scaledFT0C});
      spectra.add("hScaleMultFT0M", "Scaled total mult. signal from FT0A & FTOC", kTH1F, {scaledFT0M});

      spectra.add("hMultZNA", "Mult. signal from ZDC A-side", kTH1F, {{1000, 0.0, 5000.}});
      spectra.add("hMultZNC", "Mult. signal from ZDC C-side", kTH1F, {{1000, 0.0, 5000.}});
      spectra.add("hMultZNM", "Total mult. signal from ZDCs", kTH1F, {{4000, 0.0, 8000.}});

      // Correlations
      spectra.add("hMultFT0A_vs_ZNA", "Correlation of signals FTOA vs ZNA", kTH2F, {{2000, 0.0, 40000.}, {1000, 0.0, 5000.}});
      spectra.add("hMultFT0C_vs_ZNC", "Correlation of signals FTOC vs ZNC", kTH2F, {{2000, 0.0, 40000.}, {1000, 0.0, 5000.}});
      spectra.add("hMultFT0M_vs_ZNM", "Correlation of signals FTOM vs ZNM", kTH2F, {{3000, 0.0, 60000.}, {4000, 0.0, 8000.}});

      spectra.add("hScaleMultFT0A_vs_ZNA", "Correlation of signals FT0A/meanFT0A vs ZNA", kTH2F, {{200, 0.0, 20., "FT0A / #LT FT0A #GT"}, {1000, 0.0, 5000.}});
      spectra.add("hScaleMultFT0C_vs_ZNC", "Correlation of signals FT0C/meanFT0C vs ZNC", kTH2F, {{scaledFT0C}, {1000, 0.0, 5000.}});
      spectra.add("hScaleMultFT0M_vs_ZNM", "Correlation of signals FT0M^{*} vs ZNM", kTH2F, {{scaledFT0M}, {4000, 0.0, 8000.}});
    }
  }

  // Fill histograms with raw or MC det. level data
  template <typename Collision, typename Jets, typename Tracks>
  void fillHistograms(Collision const& collision, Jets const& jets,
                      Tracks const& tracks, float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float pTHat = getPtHat(weight);
    float rho = collision.rho();
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float scaledFT0C = getScaledFT0C(multFT0C);
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);

    auto dice = rand->Rndm();
    if (dice < fracSig)
      bSigEv = true;

    // Remove whole event if jet passes the outlier removal condition
    for (const auto& jet : jets) {
      if (jet.pt() > pTHatMax * pTHat) {
        spectra.fill(HIST("hEventSelectionCount"), 2.5);
        return;
      }
    }

    for (const auto& track : tracks) {
      spectra.fill(HIST("hTrackSelectionCount"), 0.5);

      if (skipTrack(track))
        continue;

      float trackPt = track.pt();

      spectra.fill(HIST("hTrackSelectionCount"), 1.5);
      spectra.fill(HIST("hTrackPtEtaPhi"), trackPt, track.eta(), track.phi(), weight);

      // Search for TT candidate
      if (bSigEv && (trackPt > ptTTsigMin && trackPt < ptTTsigMax)) {
        vPhiOfTT.push_back(track.phi());
        spectra.fill(HIST("hTTSig_pT"), trackPt, weight);
        ++nTT;
      }

      if (!bSigEv && (trackPt > ptTTrefMin && trackPt < ptTTrefMax)) {
        vPhiOfTT.push_back(track.phi());
        ++nTT;
      }
    }

    if (nTT > 0) { // at least 1 TT

      phiTT = getPhiTT(vPhiOfTT);

      if (bSigEv) {
        spectra.fill(HIST("hScaledFT0C_vs_Ntrig"), scaledFT0C, 1.5, weight);
        spectra.fill(HIST("hScaledFT0M_vs_Ntrig"), scaledFT0M, 1.5, weight);
        spectra.fill(HIST("hScaledFT0C_vs_TTSig_per_event"), scaledFT0C, nTT, weight);
        spectra.fill(HIST("hScaledFT0M_vs_TTSig_per_event"), scaledFT0M, nTT, weight);

        spectra.fill(HIST("hScaledFT0C_TTSig"), scaledFT0C, weight);
        spectra.fill(HIST("hScaledFT0M_TTSig"), scaledFT0M, weight);
      } else {
        spectra.fill(HIST("hScaledFT0C_vs_Ntrig"), scaledFT0C, 0.5, weight);
        spectra.fill(HIST("hScaledFT0M_vs_Ntrig"), scaledFT0M, 0.5, weight);
        spectra.fill(HIST("hScaledFT0C_vs_TTRef_per_event"), scaledFT0C, nTT, weight);
        spectra.fill(HIST("hScaledFT0M_vs_TTRef_per_event"), scaledFT0M, nTT, weight);

        spectra.fill(HIST("hScaledFT0C_TTRef"), scaledFT0C, weight);
        spectra.fill(HIST("hScaledFT0M_TTRef"), scaledFT0M, weight);
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
          spectra.fill(HIST("hScaledFT0C_Rho_TTSig"), scaledFT0C, rho, weight);
          spectra.fill(HIST("hScaledFT0M_Rho_TTSig"), scaledFT0M, rho, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTSig"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTSig"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTSig"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTSig"), scaledFT0M, dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig"), jetArea, jetPt, rho, weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_Corr_TTSig"), scaledFT0C, jetPtCorr, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_Corr_TTSig"), scaledFT0M, jetPtCorr, weight);

            spectra.fill(HIST("hScaledFT0C_Recoil_JetPt_TTSig"), scaledFT0C, jetPt, weight);
            spectra.fill(HIST("hScaledFT0M_Recoil_JetPt_TTSig"), scaledFT0M, jetPt, weight);
          }

        } else {
          spectra.fill(HIST("hScaledFT0C_Rho_TTRef"), scaledFT0C, rho, weight);
          spectra.fill(HIST("hScaledFT0M_Rho_TTRef"), scaledFT0M, rho, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_Corr_TTRef"), scaledFT0C, dphi, jetPtCorr, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_Corr_TTRef"), scaledFT0M, dphi, jetPtCorr, weight);

          spectra.fill(HIST("hScaledFT0C_DPhi_JetPt_TTRef"), scaledFT0C, dphi, jetPt, weight);
          spectra.fill(HIST("hScaledFT0M_DPhi_JetPt_TTRef"), scaledFT0M, dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTRef"), jetArea, jetPt, rho, weight);

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

  template <typename Collision, typename Jets, typename Particles>
  void fillMCPHistograms(Collision const& collision, Jets const& jets,
                         Particles const& particles, float weight = 1.)
  {
    bool bSigEv = false;
    std::vector<double> vPhiOfTT;
    double phiTT = 0.;
    int nTT = 0;
    float pTHat = getPtHat(weight);
    float rho = collision.rho();

    spectra.fill(HIST("ptHat"), pTHat, weight);

    auto dice = rand->Rndm();
    if (dice < fracSig)
      bSigEv = true;

    for (const auto& jet : jets) {
      if (jet.pt() > pTHatMax * pTHat) {
        spectra.fill(HIST("hEventSelectionCountPartLevel"), 1.5);
        return;
      }
    }

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle)
        continue;

      // Need charge and physical primary particles
      bool bParticleNeutral = (static_cast<int8_t>(pdgParticle->Charge()) == 0);
      if (bParticleNeutral || !particle.isPhysicalPrimary())
        continue;

      spectra.fill(HIST("hPartPtEtaPhi"), particle.pt(), particle.eta(),
                   particle.phi(), weight);

      if (bSigEv &&
          (particle.pt() > ptTTsigMin && particle.pt() < ptTTsigMax)) {
        vPhiOfTT.push_back(particle.phi());
        ++nTT;
      }

      if (!bSigEv &&
          (particle.pt() > ptTTrefMin && particle.pt() < ptTTrefMax)) {
        vPhiOfTT.push_back(particle.phi());
        ++nTT;
      }
    }

    if (nTT > 0) {

      phiTT = getPhiTT(vPhiOfTT);

      if (bSigEv) {
        spectra.fill(HIST("hNtrig_Part"), 1.5, weight);
        spectra.fill(HIST("hTTSig_per_event_Part"), nTT, weight);
      } else {
        spectra.fill(HIST("hNtrig_Part"), 0.5, weight);
        spectra.fill(HIST("hTTRef_per_event_Part"), nTT, weight);
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

          spectra.fill(HIST("hDPhi_JetPt_Corr_TTSig_Part"), dphi, jetPtCorr, weight);
          spectra.fill(HIST("hDPhi_JetPt_TTSig_Part"), dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTSig_Part"), jetArea, jetPt, rho, weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hRecoil_JetPt_Corr_TTSig_Part"), jetPtCorr, weight);
            spectra.fill(HIST("hRecoil_JetPt_TTSig_Part"), jetPt, weight);
          }

        } else {

          spectra.fill(HIST("hDPhi_JetPt_Corr_TTRef_Part"), dphi, jetPtCorr, weight);
          spectra.fill(HIST("hDPhi_JetPt_TTRef_Part"), dphi, jetPt, weight);
          spectra.fill(HIST("hJetArea_JetPt_Rho_TTRef_Part"), jetArea, jetPt, rho, weight);

          if (bRecoilJet) {
            spectra.fill(HIST("hRecoil_JetPt_Corr_TTRef_Part"), jetPtCorr, weight);
            spectra.fill(HIST("hRecoil_JetPt_TTRef_Part"), jetPt, weight);
          }
        }
      }
    }
  }

  template <typename TracksTable, typename JetsBase, typename JetsTag>
  void fillMatchedHistograms(TracksTable const& tracks,
                             JetsBase const& jetsBase, JetsTag const& jetsTag,
                             float weight = 1.)
  {
    std::vector<double> vPhiOfTT;
    double phiTTSig = 0.;
    float pTHat = getPtHat(weight);

    for (const auto& jetBase : jetsBase) {
      if (jetBase.pt() > pTHatMax * pTHat)
        return;
    }

    for (const auto& track : tracks) {
      if (skipTrack(track))
        continue;

      if (track.pt() > ptTTsigMin && track.pt() < ptTTsigMax) {
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

  template <typename Collision>
  void fillMultiplicityHistograms(Collision const& collision,
                                  float weight = 1.)
  {
    float multFT0A = collision.multFT0A();
    float multFT0C = collision.multFT0C();
    float multFT0M = collision.multFT0M();
    float scaledFT0A = getScaledFT0A(multFT0A);
    float scaledFT0C = getScaledFT0C(multFT0C);
    float scaledFT0M = getScaledFT0M(multFT0A, multFT0C);

    float multZNA = collision.multZNA();
    float multZNC = collision.multZNC();
    float multZNM = collision.multZNA() + collision.multZNC();

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

    // Correlations
    spectra.fill(HIST("hMultFT0A_vs_ZNA"), multFT0A, multZNA, weight);
    spectra.fill(HIST("hMultFT0C_vs_ZNC"), multFT0C, multZNC, weight);
    spectra.fill(HIST("hMultFT0M_vs_ZNM"), multFT0M, multZNM, weight);

    spectra.fill(HIST("hScaleMultFT0A_vs_ZNA"), scaledFT0A, multZNA, weight);
    spectra.fill(HIST("hScaleMultFT0C_vs_ZNC"), scaledFT0C, multZNC, weight);
    spectra.fill(HIST("hScaleMultFT0M_vs_ZNM"), scaledFT0M, multZNM, weight);
  }

  //------------------------------------------------------------------------------
  // Process functions
  void processData(FilteredColl const& collision, FilteredTracks const& tracks,
                   FilteredJets const& jets)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);

    if (skipEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    spectra.fill(HIST("vertexZ"), collision.posZ());
    fillHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(RecoilJets, processData, "process data", true);

  void processMCDetLevel(FilteredColl const& collision,
                         FilteredTracks const& tracks,
                         FilteredJetsDetLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);
    if (skipEvent(collision) || skipMBGapEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    spectra.fill(HIST("vertexZ"), collision.posZ());
    fillHistograms(collision, jets, tracks);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevel, "process MC detector level",
                 false);

  void processMCDetLevelWeighted(FilteredCollDetLevelGetWeight const& collision,
                                 aod::JetMcCollisions const&,
                                 FilteredTracks const& tracks,
                                 FilteredJetsDetLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCount"), 0.5);
    if (skipEvent(collision) || skipMBGapEvent(collision))
      return;

    spectra.fill(HIST("hEventSelectionCount"), 1.5);

    auto weight = collision.mcCollision().weight();
    spectra.fill(HIST("vertexZ"), collision.posZ(), weight);

    if (collision.has_mcCollision()) {
      spectra.fill(HIST("hHasAssocMcCollision"), 0.5, weight);
    } else {
      spectra.fill(HIST("hHasAssocMcCollision"), 1.5, weight);
    }

    fillHistograms(collision, jets, tracks, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCDetLevelWeighted,
                 "process MC detector level with event weight", false);

  void processMCPartLevel(FilteredCollPartLevel const& collision,
                          FilteredParticles const& particles,
                          FilteredJetsPartLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);
    if (skipMBGapEvent(collision))
      return;

    spectra.fill(HIST("vertexZMC"), collision.posZ());
    fillMCPHistograms(collision, jets, particles);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevel, "process MC particle level",
                 false);

  void processMCPartLevelWeighted(FilteredCollPartLevel const& collision,
                                  FilteredParticles const& particles,
                                  FilteredJetsPartLevel const& jets)
  {
    spectra.fill(HIST("hEventSelectionCountPartLevel"), 0.5);
    if (skipMBGapEvent(collision))
      return;

    auto weight = collision.weight();

    auto calcPtHat = getPtHat(weight);
    auto pThatFromFW = collision.ptHard();
    spectra.fill(HIST("hDiffInOutlierRemove"), calcPtHat - pThatFromFW);
    if (collision.isOutlier())
      spectra.fill(HIST("hCountNumberOutliersFrameWork"), 0.5);

    // LOG(debug) << "Difference between pT hat: " << calcPtHat - pThatFromFW;

    spectra.fill(HIST("vertexZMC"), collision.posZ(), weight);
    fillMCPHistograms(collision, jets, particles, weight);
  }
  PROCESS_SWITCH(RecoilJets, processMCPartLevelWeighted,
                 "process MC particle level with event weight", false);

  void processJetsMatched(FilteredCollDetLevelGetWeight const& collision,
                          aod::JetMcCollisions const&,
                          FilteredTracks const& tracks,
                          FilteredMatchedJetsDetLevel const& mcdjets,
                          FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision) || skipMBGapEvent(collision))
      return;

    auto mcpjetsPerMCCollision =
      mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());

    fillMatchedHistograms(tracks, mcpjetsPerMCCollision, mcdjets);
  }
  PROCESS_SWITCH(RecoilJets, processJetsMatched,
                 "process matching of MC jets (no weight)", false);

  void processJetsMatchedWeighted(FilteredCollDetLevelGetWeight const& collision,
                                  aod::JetMcCollisions const&,
                                  FilteredTracks const& tracks,
                                  FilteredMatchedJetsDetLevel const& mcdjets,
                                  FilteredMatchedJetsPartLevel const& mcpjets)
  {
    if (skipEvent(collision) || skipMBGapEvent(collision))
      return;

    auto mcpjetsPerMCCollision =
      mcpjets.sliceBy(partJetsPerCollision, collision.mcCollisionId());
    auto weight = collision.mcCollision().weight();

    fillMatchedHistograms(tracks, mcpjetsPerMCCollision, mcdjets, weight);
  }
  PROCESS_SWITCH(RecoilJets, processJetsMatchedWeighted,
                 "process matching of MC jets (weighted)", false);

  void processMultiplicity(FilteredEventMultiplicity const& collision)
  {
    if (skipEvent(collision))
      return;

    fillMultiplicityHistograms(collision);
  }
  PROCESS_SWITCH(RecoilJets, processMultiplicity, "process multiplicity", false);

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

  float getPtHat(float weight)
  {
    return 10. / (std::pow(weight, 1.0 / pTHatExponent));
  }

  float getScaledFT0A(const float multFT0A)
  {
    return multFT0A / meanFT0A;
  }

  float getScaledFT0C(const float multFT0C)
  {
    return multFT0C / meanFT0C;
  }

  float getScaledFT0M(const float multFT0A, const float multFT0C)
  {
    return 0.5 * (getScaledFT0A(multFT0A) + getScaledFT0C(multFT0C));
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
