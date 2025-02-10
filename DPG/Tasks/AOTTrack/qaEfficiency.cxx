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
/// \file   qaEfficiency.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to analyse both data and MC to produce efficiency vs pT, eta and phi.
///         In MC the efficiency for particles is computed according to the PDG code (sign included and not charge).
///

// O2 includes
#include <memory>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/RecoDecay.h"

// ROOT includes
#include "TPDGCode.h"
#include "TEfficiency.h"
#include "THashList.h"

using namespace o2;
using namespace o2::framework;

// Indices for the track cut histogram
static constexpr int trkCutIdxTrkRead = 1;
static constexpr int trkCutIdxHasMcPart = 2;
static constexpr int trkCutIdxPassedPt = 3;
static constexpr int trkCutIdxPassedEta = 4;
static constexpr int trkCutIdxPassedPhi = 5;
static constexpr int trkCutIdxPassedY = 6;
static constexpr int trkCutIdxPassedFake = 7;
static constexpr int trkCutIdxHasCollision = 8;
static constexpr int trkCutIdxPassedTrkType = 9;
static constexpr int trkCutIdxPassedPtRange = 10;
static constexpr int trkCutIdxPassedEtaRange = 11;
static constexpr int trkCutIdxPassedDcaXYMax = 12;
static constexpr int trkCutIdxPassedDcaXYMin = 13;
static constexpr int trkCutIdxPassedDcaZMax = 14;
static constexpr int trkCutIdxPassedDcaZMin = 15;
static constexpr int trkCutIdxPassedGoldenChi2 = 16;
static constexpr int trkCutIdxPassedIsPvCont = 17;
static constexpr int trkCutIdxPassedITSPartial = 18;
static constexpr int trkCutIdxPassedTPCPartial = 19;
static constexpr int trkCutIdxPassedTOFPartial = 20;
static constexpr int trkCutIdxPassedGlobal = 21;
static constexpr int trkCutSameColl = 22;
static constexpr int trkCutIdxN = 23;

// Particle information
static constexpr int nSpecies = o2::track::PID::NIDs; // One per PDG
static constexpr int nCharges = 2;
static constexpr int nParticles = nSpecies * nCharges;
static constexpr const char* particleTitle[nParticles] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha",
                                                          "e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
static constexpr int PDGs[nParticles] = {11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040,
                                         -11, -13, -211, -321, -2212, -1000010020, -1000010030, -1000020030, -1000020040};

// Pt
std::array<std::shared_ptr<TH1>, nParticles> hPtIts;
std::array<std::shared_ptr<TH1>, nParticles> hPtTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTof;
std::array<std::shared_ptr<TH1>, nParticles> hPtTpcTof;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTof;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTrdTof;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTrd;
std::array<std::shared_ptr<TH1>, nParticles> hPtTrkItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPtGenerated;
std::array<std::shared_ptr<TH1>, nParticles> hPtGeneratedRecoEv;

// Pt for primaries
std::array<std::shared_ptr<TH1>, nParticles> hPtItsPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPtTrkItsTpcPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTofPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPtTrkItsTpcTofPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPtGeneratedPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPtGeneratedPrmRecoEv;

// Pt for secondaries from weak decay
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcStr;
std::array<std::shared_ptr<TH1>, nParticles> hPtTrkItsTpcStr;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTofStr;
std::array<std::shared_ptr<TH1>, nParticles> hPtGeneratedStr;
std::array<std::shared_ptr<TH1>, nParticles> hPtmotherGenerated; // histogram to store pT of mother
std::array<std::shared_ptr<TH1>, nParticles> hdecaylengthmother; // histogram to store decaylength of mother

// Pt for secondaries from material
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcMat;
std::array<std::shared_ptr<TH1>, nParticles> hPtTrkItsTpcMat;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTofMat;
std::array<std::shared_ptr<TH1>, nParticles> hPtGeneratedMat;

// Pt for tertiaries from secondary weak decay
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTer;
std::array<std::shared_ptr<TH1>, nParticles> hPtTrkItsTpcTer;
std::array<std::shared_ptr<TH1>, nParticles> hPtItsTpcTofTer;
std::array<std::shared_ptr<TH1>, nParticles> hPtGeneratedTer;

// P
std::array<std::shared_ptr<TH1>, nParticles> hPItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPTrkItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPItsTpcTof;
std::array<std::shared_ptr<TH1>, nParticles> hPGenerated;

// Eta
std::array<std::shared_ptr<TH1>, nParticles> hEtaItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hEtaTrkItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hEtaItsTpcTof;
std::array<std::shared_ptr<TH1>, nParticles> hEtaGenerated;

// Eta for primaries
std::array<std::shared_ptr<TH1>, nParticles> hEtaItsTpcPrm;
std::array<std::shared_ptr<TH1>, nParticles> hEtaTrkItsTpcPrm;
std::array<std::shared_ptr<TH1>, nParticles> hEtaItsTpcTofPrm;
std::array<std::shared_ptr<TH1>, nParticles> hEtaGeneratedPrm;

// Y
std::array<std::shared_ptr<TH1>, nParticles> hYItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hYItsTpcTof;
std::array<std::shared_ptr<TH1>, nParticles> hYGenerated;

// Phi
std::array<std::shared_ptr<TH1>, nParticles> hPhiItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPhiTrkItsTpc;
std::array<std::shared_ptr<TH1>, nParticles> hPhiItsTpcTof;
std::array<std::shared_ptr<TH1>, nParticles> hPhiGenerated;

// Phi for primaries
std::array<std::shared_ptr<TH1>, nParticles> hPhiItsTpcPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPhiTrkItsTpcPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPhiItsTpcTofPrm;
std::array<std::shared_ptr<TH1>, nParticles> hPhiGeneratedPrm;

// 2D
std::array<std::shared_ptr<TH2>, nParticles> hPtEtaItsTpc;
std::array<std::shared_ptr<TH2>, nParticles> hPtEtaTrkItsTpc;
std::array<std::shared_ptr<TH2>, nParticles> hPtEtaItsTpcTof;
std::array<std::shared_ptr<TH2>, nParticles> hPtEtaGenerated;
// 2D  Pt vs Radius
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpc;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusTrkItsTpc;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcTof;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusGenerated;

std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcPrm;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusTrkItsTpcPrm;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcTofPrm;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusGeneratedPrm;

std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcStr;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusTrkItsTpcStr;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcTofStr;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusGeneratedStr;

std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcTer;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusTrkItsTpcTer;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusItsTpcTofTer;
std::array<std::shared_ptr<TH2>, nParticles> hPtRadiusGeneratedTer;

struct QaEfficiency {
  // Track/particle selection
  Configurable<bool> numSameCollision{"numSameCollision", false, "Flag to ask that the numerator is in the same collision as the denominator"};
  Configurable<bool> noFakesHits{"noFakesHits", false, "Flag to reject tracks that have fake hits"};
  Configurable<bool> skipEventsWithoutTPCTracks{"skipEventsWithoutTPCTracks", false, "Flag to reject events that have no tracks reconstructed in the TPC"};
  Configurable<float> maxProdRadius{"maxProdRadius", 9999.f, "Maximum production radius of the particle under study"};
  Configurable<float> nsigmaTPCDe{"nsigmaTPCDe", 3.f, "Value of the Nsigma TPC cut for deuterons PID"};
  // Charge selection
  Configurable<bool> doPositivePDG{"doPositivePDG", false, "Flag to fill histograms for positive PDG codes."};
  Configurable<bool> doNegativePDG{"doNegativePDG", false, "Flag to fill histograms for negative PDG codes."};
  // Particle only selection
  Configurable<bool> doEl{"do-el", false, "Flag to run with the PDG code of electrons"};
  Configurable<bool> doMu{"do-mu", false, "Flag to run with the PDG code of muons"};
  Configurable<bool> doPi{"do-pi", false, "Flag to run with the PDG code of pions"};
  Configurable<bool> doKa{"do-ka", false, "Flag to run with the PDG code of kaons"};
  Configurable<bool> doPr{"do-pr", false, "Flag to run with the PDG code of protons"};
  Configurable<bool> doDe{"do-de", false, "Flag to run with the PDG code of deuterons"};
  Configurable<bool> doTr{"do-tr", false, "Flag to run with the PDG code of tritons"};
  Configurable<bool> doHe{"do-he", false, "Flag to run with the PDG code of helium 3"};
  Configurable<bool> doAl{"do-al", false, "Flag to run with the PDG code of helium 4"};
  // Selection on mothers
  Configurable<bool> checkForMothers{"checkForMothers", false, "Flag to use the array of mothers to check if the particle of interest come from any of those particles"};
  Configurable<std::vector<int>> mothersPDGs{"mothersPDGs", std::vector<int>{3312, -3312}, "PDGs of origin of the particle under study"};
  Configurable<bool> keepOnlyHfParticles{"keepOnlyHfParticles", false, "Flag to decide wheter to consider only HF particles"};
  Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "Flag to check specific event generator (for HF): -1 -> no check, 0 -> MB events, 4 -> charm triggered, 5 -> beauty triggered"};
  // Track only selection, options to select only specific tracks
  Configurable<bool> trackSelection{"trackSelection", true, "Local track selection"};
  Configurable<int> globalTrackSelection{"globalTrackSelection", 0, "Global track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks, 6 -> custom track cuts via Configurable"};
  // Event selection
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the primary vertez in Z (cm)"};
  Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the primary vertez in Z (cm)"};
  Configurable<bool> applyPvZCutGenColl{"applyPvZCutGenColl", false, "Flag to enable the cut on the generated vertex z coordinate"};
  Configurable<bool> applyPvZCutInProcessMcWoColl{"applyPvZCutInProcessMcWoColl", false, "Flag to enable the cut on the vertex z coordinate (reco. & gen.) also in processMCWithoutCollisions"};
  // Histogram configuration
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  Configurable<int> logPt{"log-pt", 0, "Flag to use a logarithmic pT axis"};
  ConfigurableAxis etaBins{"etaBins", {200, -3.f, 3.f}, "Eta binning"};
  ConfigurableAxis phiBins{"phiBins", {200, 0.f, 6.284f}, "Phi binning"};
  ConfigurableAxis yBins{"yBins", {200, -0.5f, 0.5f}, "Y binning"};
  ConfigurableAxis occBins{"occBins", {100, 0.f, 14000.f}, "Occupancy binning"};
  ConfigurableAxis centBins{"centBins", {110, 0.f, 110.f}, "Centrality binning"};
  ConfigurableAxis radiusBins{"radiusBins", {200, 0.f, 100.f}, "Radius binning"};
  // Task configuration
  Configurable<bool> makeEff{"make-eff", false, "Flag to produce the efficiency with TEfficiency"};
  Configurable<bool> doPtEta{"doPtEta", false, "Flag to produce the efficiency vs pT and Eta"};
  Configurable<bool> doPtRadius{"doPtRadius", false, "Flag to produce the efficiency vs pT and Radius"};
  Configurable<int> applyEvSel{"applyEvSel", 0, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  // Custom track cuts for debug purposes
  TrackSelection customTrackCuts;
  struct : ConfigurableGroup {
    Configurable<bool> tracksIU{"tracksIU", false, "Additional cut for IU tracks"};
    Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
    Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
    Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
    Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
    Configurable<int> minITScl{"minITScl", 4, "Additional cut on the ITS cluster"};
    Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
    Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
    Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
    Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
    Configurable<float> maxDcaXY{"maxDcaXY", 10000.f, "Additional cut on the maximum abs value of the DCA xy"};
    Configurable<float> maxDcaZ{"maxDcaZ", 2.f, "Additional cut on the maximum abs value of the DCA z"};
    Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  } cfgCustomTrackCuts;

  Configurable<bool> doPVContributorCut{"doPVContributorCut", false, "Select tracks used for primary vertex recostruction (isPVContributor)"};
  Configurable<float> minDcaZ{"minDcaZ", -2.f, "Additional cut on the minimum abs value of the DCA z"};
  Configurable<float> minDcaXY{"minDcaXY", -1.f, "Additional cut on the minimum abs value of the DCA xy"};

  Configurable<bool> doOccupancy{"doOccupancyStudy", false, "Flag to store Occupancy-related information"};
  Configurable<bool> useFT0OccEstimator{"useFT0OccEstimator", false, "Flag to adopt FT0c to estimate occupancy instead of ITS"};

  // Output objects for TEfficiency
  OutputObj<THashList> listEfficiencyMC{"EfficiencyMC"};
  OutputObj<THashList> listEfficiencyData{"EfficiencyData"};

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, aod::CentFT0Cs>;
  using CollisionCandidatesMC = o2::soa::Join<CollisionCandidates, o2::aod::McCollisionLabels>;
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension, o2::aod::TracksDCA>;
  using TrackCandidatesMC = o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels>;

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static const char* particleName(int pdgSign, o2::track::PID::ID id)
  {
    if (pdgSign == 0) { // Positive PDG
      return particleTitle[id];
    }
    // Negative PDG
    return particleTitle[id + o2::track::PID::NIDs];
  }

  template <int pdgSign, o2::track::PID::ID id>
  void makeMCHistograms(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    if constexpr (pdgSign == 0) {
      if (!doPositivePDG) { // Positive
        return;
      }
    } else if constexpr (pdgSign == 1) {
      if (!doNegativePDG) { // Negative
        return;
      }
    } else {
      LOG(fatal) << "Can't interpret pdgSign " << pdgSign;
    }

    AxisSpec axisDecayLength{100, 0.0, 10.0, "Decay Length (cm)"};
    AxisSpec axisPt{ptBins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisP{ptBins, "#it{p} (GeV/#it{c})"};
    if (logPt) {
      axisPt.makeLogarithmic();
      axisP.makeLogarithmic();
    }
    const AxisSpec axisEta{etaBins, "#it{#eta}"};
    const AxisSpec axisY{yBins, "#it{y}"};
    const AxisSpec axisPhi{phiBins, "#it{#varphi} (rad)"};
    const AxisSpec axisRadius{radiusBins, "Radius (cm)"};
    const AxisSpec axisOcc{occBins, "Occupancy"};

    const char* partName = particleName(pdgSign, id);
    LOG(info) << "Preparing histograms for particle: " << partName << " pdgSign " << pdgSign;

    const TString tagPt = Form("%s #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                               partName,
                               etaMin, etaMax,
                               yMin, yMax,
                               phiMin, phiMax);

    const TString tagEta = Form("%s #it{p}_{T} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                                partName,
                                ptMin, ptMax,
                                yMin, yMax,
                                phiMin, phiMax);

    const TString tagY = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                              partName,
                              ptMin, ptMax,
                              etaMin, etaMax,
                              phiMin, phiMax);

    const TString tagPhi = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
                                partName,
                                ptMin, ptMax,
                                etaMin, etaMax,
                                yMin, yMax);

    const TString tagPtEta = Form("%s #it{#varphi} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
                                  partName,
                                  phiMin, phiMax,
                                  yMin, yMax);
    const int histogramIndex = id + pdgSign * nSpecies;

    // Pt
    hPtIts[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/its", PDGs[histogramIndex]), "ITS tracks " + tagPt, kTH1D, {axisPt});
    hPtTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/tpc", PDGs[histogramIndex]), "TPC tracks " + tagPt, kTH1D, {axisPt});
    hPtItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPt, kTH1D, {axisPt});
    hPtItsTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/its_tof", PDGs[histogramIndex]), "ITS-TOF tracks " + tagPt, kTH1D, {axisPt});
    hPtTpcTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/tpc_tof", PDGs[histogramIndex]), "TPC-TOF tracks " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTrdTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/its_tpc_trd_tof", PDGs[histogramIndex]), "ITS-TPC-TRD-TOF tracks " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTrd[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/its_tpc_trd", PDGs[histogramIndex]), "ITS-TPC-TRD tracks " + tagPt, kTH1D, {axisPt});
    hPtTrkItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC track (reco) " + tagPt, kTH1D, {axisPt});
    hPtGenerated[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/generated", PDGs[histogramIndex]), "Generated " + tagPt, kTH1D, {axisPt});
    hPtGeneratedRecoEv[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/generated_reco_ev", PDGs[histogramIndex]), "Generated Reco Ev. " + tagPt, kTH1D, {axisPt});

    // Prm
    hPtItsPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/its", PDGs[histogramIndex]), "ITS tracks (primaries) " + tagPt, kTH1D, {axisPt});
    hPtItsTpcPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (primaries) " + tagPt, kTH1D, {axisPt});
    hPtTrkItsTpcPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco primaries) " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTofPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (primaries) " + tagPt, kTH1D, {axisPt});
    hPtTrkItsTpcTofPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/trk/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (reco primaries) " + tagPt, kTH1D, {axisPt});
    hPtGeneratedPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/generated", PDGs[histogramIndex]), "Generated (primaries) " + tagPt, kTH1D, {axisPt});
    hPtGeneratedPrmRecoEv[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/prm/generated_reco_ev", PDGs[histogramIndex]), "Generated Reco Ev. " + tagPt, kTH1D, {axisPt});

    // Str
    hPtItsTpcStr[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/str/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (from weak decays) " + tagPt, kTH1D, {axisPt});
    hPtTrkItsTpcStr[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/str/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco from weak decays) " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTofStr[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/str/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (from weak decays) " + tagPt, kTH1D, {axisPt});
    hPtGeneratedStr[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/str/generated", PDGs[histogramIndex]), "Generated (from weak decays) " + tagPt, kTH1D, {axisPt});
    hPtmotherGenerated[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/str/generated_mother", PDGs[histogramIndex]), "Generated Mother " + tagPt, kTH1D, {axisPt});
    hdecaylengthmother[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/str/decayLength", PDGs[histogramIndex]), "Decay Length of mother particle" + tagPt, kTH1D, {axisDecayLength});

    // Ter
    hPtItsTpcTer[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/ter/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (from secondary weak decays) " + tagPt, kTH1D, {axisPt});
    hPtTrkItsTpcTer[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/ter/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco from secondary weak decays) " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTofTer[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/ter/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (from secondary weak decays) " + tagPt, kTH1D, {axisPt});
    hPtGeneratedTer[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/ter/generated", PDGs[histogramIndex]), "Generated (from secondary weak decays) " + tagPt, kTH1D, {axisPt});

    // Mat
    hPtItsTpcMat[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/mat/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (from material)" + tagPt, kTH1D, {axisPt});
    hPtTrkItsTpcMat[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/mat/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco from material) " + tagPt, kTH1D, {axisPt});
    hPtItsTpcTofMat[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/mat/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (from material) " + tagPt, kTH1D, {axisPt});
    hPtGeneratedMat[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/pt/mat/generated", PDGs[histogramIndex]), "Generated ( from material) " + tagPt, kTH1D, {axisPt});

    // P
    hPItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/p/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPt, kTH1D, {axisP});
    hPTrkItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/p/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco) " + tagPt, kTH1D, {axisP});
    hPItsTpcTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/p/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPt, kTH1D, {axisP});
    hPGenerated[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/p/generated", PDGs[histogramIndex]), "Generated " + tagPt, kTH1D, {axisP});

    // Eta
    hEtaItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagEta, kTH1D, {axisEta});
    hEtaTrkItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco) " + tagEta, kTH1D, {axisEta});
    hEtaItsTpcTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagEta, kTH1D, {axisEta});
    hEtaGenerated[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/generated", PDGs[histogramIndex]), "Generated " + tagEta, kTH1D, {axisEta});

    // Prm
    hEtaItsTpcPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/prm/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (primaries) " + tagEta, kTH1D, {axisEta});
    hEtaTrkItsTpcPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/prm/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco primaries) " + tagEta, kTH1D, {axisEta});
    hEtaItsTpcTofPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/prm/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (primaries) " + tagEta, kTH1D, {axisEta});
    hEtaGeneratedPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/eta/prm/generated", PDGs[histogramIndex]), "Generated (primaries) " + tagEta, kTH1D, {axisEta});

    // Y
    hYItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/y/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagY, kTH1D, {axisY});
    hYItsTpcTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/y/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagY, kTH1D, {axisY});
    hYGenerated[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/y/generated", PDGs[histogramIndex]), "Generated " + tagY, kTH1D, {axisY});

    // Phi
    hPhiItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPhi, kTH1D, {axisPhi});
    hPhiTrkItsTpc[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco) " + tagPhi, kTH1D, {axisPhi});
    hPhiItsTpcTof[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPhi, kTH1D, {axisPhi});
    hPhiGenerated[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/generated", PDGs[histogramIndex]), "Generated " + tagPhi, kTH1D, {axisPhi});

    // Phi prm
    hPhiItsTpcPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/prm/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (primaries) " + tagPhi, kTH1D, {axisPhi});
    hPhiTrkItsTpcPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/prm/trk/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks (reco primaries) " + tagPhi, kTH1D, {axisPhi});
    hPhiItsTpcTofPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/prm/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks (primaries) " + tagPhi, kTH1D, {axisPhi});
    hPhiGeneratedPrm[histogramIndex] = histos.add<TH1>(Form("MC/pdg%i/phi/prm/generated", PDGs[histogramIndex]), "Generated (primaries) " + tagPhi, kTH1D, {axisPhi});

    if (doPtEta) {
      hPtEtaItsTpc[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/asd", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPtEta, kTH2D, {axisPt, axisEta});
      hPtEtaTrkItsTpc[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/asd", PDGs[histogramIndex]), "ITS-TPC tracks (reco) " + tagPtEta, kTH2D, {axisPt, axisEta});
      hPtEtaItsTpcTof[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/asd", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPtEta, kTH2D, {axisPt, axisEta});
      hPtEtaGenerated[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/asd", PDGs[histogramIndex]), "Generated " + tagPtEta, kTH2D, {axisPt, axisEta});
    }

    if (doPtRadius) {
      hPtRadiusItsTpc[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/radius/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusItsTpcTof[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/radius/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusGenerated[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/radius/generated", PDGs[histogramIndex]), "Generated " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});

      hPtRadiusItsTpcPrm[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/prm/radius/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusItsTpcTofPrm[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/prm/radius/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusGeneratedPrm[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/prm/radius/generated", PDGs[histogramIndex]), "Generated " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});

      hPtRadiusItsTpcStr[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/str/radius/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusItsTpcTofStr[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/str/radius/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusGeneratedStr[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/str/radius/generated", PDGs[histogramIndex]), "Generated " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});

      hPtRadiusItsTpcTer[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/ter/radius/its_tpc", PDGs[histogramIndex]), "ITS-TPC tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusItsTpcTofTer[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/ter/radius/its_tpc_tof", PDGs[histogramIndex]), "ITS-TPC-TOF tracks " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
      hPtRadiusGeneratedTer[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/ter/radius/generated", PDGs[histogramIndex]), "Generated " + tagPt + " vs Radius", kTH2D, {axisPt, axisRadius});
    }

    LOG(info) << "Done with making histograms for particle: " << partName;
  }

  template <int pdgSign, o2::track::PID::ID id>
  void makeMCEfficiency(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    if (!makeEff) {
      return;
    }

    if constexpr (pdgSign == 0) {
      if (!doPositivePDG) { // Positive
        return;
      }
    } else if constexpr (pdgSign == 1) {
      if (!doNegativePDG) { // Negative
        return;
      }
    } else {
      LOG(fatal) << "Can't interpret pdgSign index";
    }

    const TString partName = particleName(pdgSign, id);
    LOG(info) << "Making TEfficiency for MC for particle " << partName;
    THashList* subList = new THashList();
    subList->SetName(partName);
    listEfficiencyMC->Add(subList);

    auto makeEfficiency = [&](const TString effname, auto h) { // 1D efficiencies
      LOG(debug) << " Making 1D TEfficiency " << effname << " from " << h->GetName();
      const TAxis* axis = h->GetXaxis();
      TString efftitle = h->GetTitle();
      efftitle.ReplaceAll("Numerator", "").Strip(TString::kBoth);
      efftitle = Form("%s;%s;Efficiency", efftitle.Data(), axis->GetTitle());
      if (axis->IsVariableBinSize()) {
        subList->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
      } else {
        subList->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
      }
    };

    const int histogramIndex = id + pdgSign * nSpecies;

    makeEfficiency("ITS_vsPt", hPtIts[histogramIndex]);
    makeEfficiency("TPC_vsPt", hPtTpc[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt", hPtItsTpc[histogramIndex]);
    makeEfficiency("ITS-TOF_vsPt", hPtItsTof[histogramIndex]);
    makeEfficiency("TPC-TOF_vsPt", hPtTpcTof[histogramIndex]);
    makeEfficiency("ITS-TPC-TRD_vsPt", hPtItsTpcTrd[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPt", hPtItsTpcTof[histogramIndex]);
    makeEfficiency("ITS-TPC-TRD-TOF_vsPt", hPtItsTpcTrdTof[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Trk", hPtTrkItsTpc[histogramIndex]);

    makeEfficiency("ITS-TPC_vsPt_RecoEv", hPtItsTpc[histogramIndex]);

    makeEfficiency("ITS_vsPt_Prm", hPtItsPrm[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Prm", hPtItsTpcPrm[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Prm_Trk", hPtTrkItsTpcPrm[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPt_Prm", hPtItsTpcTofPrm[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPt_Prm_Trk", hPtTrkItsTpcTofPrm[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Prm_RecoEv", hPtItsTpcPrm[histogramIndex]);

    makeEfficiency("ITS-TPC_vsPt_Str", hPtItsTpcStr[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Str_Trk", hPtTrkItsTpcStr[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPt_Str", hPtItsTpcTofStr[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Mat", hPtItsTpcMat[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Mat_Trk", hPtTrkItsTpcMat[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPt_Mat", hPtItsTpcTofMat[histogramIndex]);

    makeEfficiency("ITS-TPC_vsPt_Ter", hPtItsTpcTer[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPt_Ter_Trk", hPtTrkItsTpcTer[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPt_Ter", hPtItsTpcTofTer[histogramIndex]);

    makeEfficiency("ITS-TPC_vsP", hPItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC_vsP_Trk", hPTrkItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsP", hPItsTpcTof[histogramIndex]);
    makeEfficiency("ITS-TPC_vsEta", hEtaItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC_vsEta_Trk", hEtaTrkItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsEta", hEtaItsTpcTof[histogramIndex]);
    makeEfficiency("ITS-TPC_vsEta_Prm", hEtaItsTpcPrm[histogramIndex]);
    makeEfficiency("ITS-TPC_vsEta_Prm_Trk", hEtaTrkItsTpcPrm[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsEta_Prm", hEtaItsTpcTofPrm[histogramIndex]);
    makeEfficiency("ITS-TPC_vsY", hYItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsY", hYItsTpcTof[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPhi", hPhiItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPhi_Trk", hPhiTrkItsTpc[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPhi", hPhiItsTpcTof[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPhi_Prm", hPhiItsTpcPrm[histogramIndex]);
    makeEfficiency("ITS-TPC_vsPhi_Prm_Trk", hPhiTrkItsTpcPrm[histogramIndex]);
    makeEfficiency("ITS-TPC-TOF_vsPhi_Prm", hPhiItsTpcTofPrm[histogramIndex]);

    auto makeEfficiency2D = [&](const TString effname, auto h) { // 2D efficiencies
      LOG(debug) << " Making 2D TEfficiency " << effname << " from " << h->GetName();
      const TAxis* axisX = h->GetXaxis();
      const TAxis* axisY = h->GetYaxis();
      TString efftitle = h->GetTitle();
      efftitle.ReplaceAll("Numerator", "").Strip(TString::kBoth);
      efftitle = Form("%s;%s;%s;Efficiency", efftitle.Data(), axisX->GetTitle(), axisY->GetTitle());
      if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
        subList->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
      } else {
        subList->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
      }
    };

    if (doPtEta) {
      makeEfficiency2D("ITS-TPC_vsPt_vsEta", hPtEtaItsTpc[histogramIndex]);
      makeEfficiency2D("ITS-TPC_vsPt_vsEta_Trk", hPtEtaTrkItsTpc[histogramIndex]);
      makeEfficiency2D("ITS-TPC-TOF_vsPt_vsEta", hPtEtaItsTpcTof[histogramIndex]);
    }
    if (doPtRadius) {
      makeEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpc[histogramIndex]);
      makeEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTof[histogramIndex]);
      makeEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpcPrm[histogramIndex]);
      makeEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTofPrm[histogramIndex]);
      makeEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpcStr[histogramIndex]);
      makeEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTofStr[histogramIndex]);
      makeEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpcTer[histogramIndex]);
      makeEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTofTer[histogramIndex]);
    }

    LOG(info) << "Done with making histograms for particle: " << partName << " for efficiencies";
  }

  void initMC(const AxisSpec& axisSel)
  {
    if (!doprocessMC && !doprocessMCWithoutCollisions) {
      return;
    }
    if (doprocessMC && doprocessMCWithoutCollisions) {
      LOG(fatal) << "Both processMC and processMCWithoutCollisions are set to true. Please set only one of them to true.";
    }

    auto h = histos.add<TH1>("MC/trackSelection", "Track Selection", kTH1D, {axisSel});
    h->GetXaxis()->SetBinLabel(trkCutIdxTrkRead, "Tracks read");
    h->GetXaxis()->SetBinLabel(trkCutIdxHasMcPart, "Passed has MC part.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedPt, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedEta, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedPhi, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedY, "Passed y");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedFake, "Passed Fake");
    h->GetXaxis()->SetBinLabel(trkCutIdxHasCollision, "Passed has collision");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedTrkType, "passedTrackType");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedPtRange, "passedPtRange");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedEtaRange, "passedEtaRange");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaXYMax, "passedDCAxy max.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaXYMin, "passedDCAxy min.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaZMax, "passedDCAz max.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaZMin, "passedDCAz min.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedGoldenChi2, "passedGoldenChi2");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedIsPvCont, "passed isPVContributor");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedITSPartial, "passedITS (partial)");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedTPCPartial, "passedTPC (partial)");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedTOFPartial, "passedTOF (partial)");
    switch (globalTrackSelection) {
      case 0:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "No extra selection");
        break;
      case 1:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isGlobalTrack");
        break;
      case 2:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isGlobalTrackWoPtEta");
        break;
      case 3:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isGlobalTrackWoDCA");
        break;
      case 4:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isQualityTrack");
        break;
      case 5:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isInAcceptanceTrack");
        break;
      case 6:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "customTrackSelection");
        break;
      default:
        LOG(fatal) << "Can't interpret track asked selection " << globalTrackSelection;
    }
    h->GetXaxis()->SetBinLabel(trkCutSameColl, "passedSameColl");

    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(trkCutIdxN + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1D, {{1, 0, 1}});

    h = histos.add<TH1>("MC/particleSelection", "Particle Selection", kTH1D, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Particles read");
    h->GetXaxis()->SetBinLabel(2, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(5, "Passed y");
    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(6 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/eventMultiplicity", "Event Selection", kTH1D, {{1000, 0, 5000}});

    histos.add("MC/trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    listEfficiencyMC.setObject(new THashList);

    if (doOccupancy) {
      const AxisSpec axisPt{ptBins, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec axisOcc{occBins, "Occupancy"};
      const AxisSpec axisCent{centBins, "Centrality"};

      histos.add("MC/occ_cent/gen/pos", "Generated Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("MC/occ_cent/gen/neg", "Generated Negative ", kTH3D, {axisOcc, axisCent, axisPt});

      histos.add("MC/occ_cent/reco/pos/its_tpc_tof", "ITS-TPC-TOF Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("MC/occ_cent/reco/neg/its_tpc_tof", "ITS-TPC-TOF Negative ", kTH3D, {axisOcc, axisCent, axisPt});

      histos.add("MC/occ_cent/reco/pos/its_tpc", "ITS-TPC Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("MC/occ_cent/reco/neg/its_tpc", "ITS-TPC Negative ", kTH3D, {axisOcc, axisCent, axisPt});

      histos.add("MC/occ_cent/reco/pos/its", "ITS Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("MC/occ_cent/reco/neg/its", "ITS Negative ", kTH3D, {axisOcc, axisCent, axisPt});
    }

    static_for<0, 1>([&](auto pdgSign) {
      makeMCHistograms<pdgSign, o2::track::PID::Electron>(doEl);
      makeMCHistograms<pdgSign, o2::track::PID::Muon>(doMu);
      makeMCHistograms<pdgSign, o2::track::PID::Pion>(doPi);
      makeMCHistograms<pdgSign, o2::track::PID::Kaon>(doKa);
      makeMCHistograms<pdgSign, o2::track::PID::Proton>(doPr);
      makeMCHistograms<pdgSign, o2::track::PID::Deuteron>(doDe);
      makeMCHistograms<pdgSign, o2::track::PID::Triton>(doTr);
      makeMCHistograms<pdgSign, o2::track::PID::Helium3>(doHe);
      makeMCHistograms<pdgSign, o2::track::PID::Alpha>(doAl);

      makeMCEfficiency<pdgSign, o2::track::PID::Electron>(doEl);
      makeMCEfficiency<pdgSign, o2::track::PID::Muon>(doMu);
      makeMCEfficiency<pdgSign, o2::track::PID::Pion>(doPi);
      makeMCEfficiency<pdgSign, o2::track::PID::Kaon>(doKa);
      makeMCEfficiency<pdgSign, o2::track::PID::Proton>(doPr);
      makeMCEfficiency<pdgSign, o2::track::PID::Deuteron>(doDe);
      makeMCEfficiency<pdgSign, o2::track::PID::Triton>(doTr);
      makeMCEfficiency<pdgSign, o2::track::PID::Helium3>(doHe);
      makeMCEfficiency<pdgSign, o2::track::PID::Alpha>(doAl);
    });
  }

  // Efficiencies
  TEfficiency* effITSTPCMatchingVsPt = nullptr;
  TEfficiency* effTPCITSMatchingVsPt = nullptr;
  TEfficiency* effTPCTOFMatchingVsPt = nullptr;
  TEfficiency* effTPCTOFMatchingVsP = nullptr;
  TEfficiency* effTPCTOFMatchingVsEta = nullptr;
  TEfficiency* effTPCTOFMatchingVsPhi = nullptr;

  // 2D
  TEfficiency* effTPCTOFMatchingVsPtVsEta = nullptr;
  TEfficiency* effTPCTOFMatchingVsPtVsPhi = nullptr;

  void initData(const AxisSpec& axisSel)
  {
    if (!doprocessData && !doprocessDataWithPID) {
      return;
    }

    if (doprocessData == true && doprocessDataWithPID == true) {
      LOG(fatal) << "Can't enable processData and doprocessDataWithPID in the same time, pick one!";
    }

    auto h = histos.add<TH1>("Data/trackSelection", "Track Selection", kTH1D, {axisSel});
    h->GetXaxis()->SetBinLabel(trkCutIdxTrkRead, "Tracks read");
    h->GetXaxis()->SetBinLabel(trkCutIdxHasMcPart, "");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedPt, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedEta, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedPhi, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedY, "");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedFake, "");
    h->GetXaxis()->SetBinLabel(trkCutIdxHasCollision, "Passed has collision");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedTrkType, "passedTrackType");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedPtRange, "passedPtRange");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedEtaRange, "passedEtaRange");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaXYMax, "passedDCAxy max.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaXYMin, "passedDCAxy min.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaZMax, "passedDCAz max.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedDcaZMin, "passedDCAz min.");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedGoldenChi2, "passedGoldenChi2");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedIsPvCont, "passed isPVContributor");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedITSPartial, "passedITS (partial)");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedTPCPartial, "passedTPC (partial)");
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedTOFPartial, "passedTOF (partial)");
    switch (globalTrackSelection) {
      case 0:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "No extra selection");
        break;
      case 1:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isGlobalTrack");
        break;
      case 2:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isGlobalTrackWoPtEta");
        break;
      case 3:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isGlobalTrackWoDCA");
        break;
      case 4:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isQualityTrack");
        break;
      case 5:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "isInAcceptanceTrack");
        break;
      case 6:
        h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "customTrackSelection");
        break;
      default:
        LOG(fatal) << "Can't interpret track asked selection " << globalTrackSelection;
    }

    const TString tagPt = Form("#it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                               etaMin, etaMax,
                               phiMin, phiMax);
    AxisSpec axisPt{ptBins, "#it{p}_{T} (GeV/#it{c})"};
    if (logPt) {
      axisPt.makeLogarithmic();
    }

    const TString tagEta = Form("#it{p}_{T} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                                ptMin, ptMax,
                                phiMin, phiMax);
    const AxisSpec axisEta{etaBins, "#it{#eta}"};
    const AxisSpec axisOcc{occBins, "Occupancy"};
    const AxisSpec axisCent{centBins, "Centrality"};

    const TString tagPhi = Form("#it{#eta} [%.2f,%.2f] #it{p}_{T} [%.2f,%.2f]",
                                etaMin, etaMax,
                                ptMin, ptMax);
    const AxisSpec axisPhi{phiBins, "#it{#varphi} (rad)"};

    const TString tagEtaPhi = Form("#it{p}_{T} [%.2f,%.2f]",
                                   ptMin, ptMax);

    histos.add("Data/trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    if (doOccupancy) {
      histos.add("Data/occ_cent/pos/its_tpc_tof", "ITS-TPC-TOF Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("Data/occ_cent/neg/its_tpc_tof", "ITS-TPC-TOF Negative ", kTH3D, {axisOcc, axisCent, axisPt});

      histos.add("Data/occ_cent/pos/its_tpc", "ITS-TPC Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("Data/occ_cent/neg/its_tpc", "ITS-TPC Negative ", kTH3D, {axisOcc, axisCent, axisPt});

      histos.add("Data/occ_cent/pos/tpc", "TPC Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("Data/occ_cent/neg/tpc", "TPC Negative ", kTH3D, {axisOcc, axisCent, axisPt});

      histos.add("Data/occ_cent/pos/its", "ITS Positive ", kTH3D, {axisOcc, axisCent, axisPt});
      histos.add("Data/occ_cent/neg/its", "ITS Negative ", kTH3D, {axisOcc, axisCent, axisPt});
    }

    // ITS-TPC-TOF
    histos.add("Data/pos/pt/its_tpc_tof", "ITS-TPC-TOF Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/its_tpc_tof", "ITS-TPC-TOF Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/pos/eta/its_tpc_tof", "ITS-TPC-TOF Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/its_tpc_tof", "ITS-TPC-TOF Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/pos/phi/its_tpc_tof", "ITS-TPC-TOF Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/its_tpc_tof", "ITS-TPC-TOF Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/pos/etaphi/its_tpc_tof", "ITS-TPC-TOF Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc_tof", "ITS-TPC-TOF Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // ITS-TPC
    histos.add("Data/pos/pt/its_tpc", "ITS-TPC Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/its_tpc", "ITS-TPC Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/pos/eta/its_tpc", "ITS-TPC Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/its_tpc", "ITS-TPC Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/pos/phi/its_tpc", "ITS-TPC Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/its_tpc", "ITS-TPC Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/pos/etaphi/its_tpc", "ITS-TPC Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc", "ITS-TPC Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // TPC
    histos.add("Data/pos/pt/tpc", "TPC Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/tpc", "TPC Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/pos/eta/tpc", "TPC Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/tpc", "TPC Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/pos/phi/tpc", "TPC Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/tpc", "TPC Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/pos/etaphi/tpc", "TPC Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/tpc", "TPC Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // ITS
    histos.add("Data/pos/pt/its", "ITS Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/its", "ITS Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/pos/eta/its", "ITS Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/its", "ITS Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/pos/phi/its", "ITS Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/its", "ITS Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/pos/etaphi/its", "ITS Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its", "ITS Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // HMPID
    if (doprocessHmpid) {
      histos.add("Data/pos/hmpidMomDiff", "HMPID Positive Momentum difference", kTH1D, {{100, -10, 10, "#it{p} - #it{p}_{HMPID} (GeV/#it{c})"}});
      histos.add("Data/neg/hmpidMomDiff", "HMPID Negative Momentum difference", kTH1D, {{100, -10, 10, "#it{p} - #it{p}_{HMPID} (GeV/#it{c})"}});
      histos.add("Data/pos/pt/hmpid", "HMPID Positive " + tagPt, kTH1D, {axisPt});
      histos.add("Data/neg/pt/hmpid", "HMPID Negative " + tagPt, kTH1D, {axisPt});

      histos.add("Data/pos/eta/hmpid", "HMPID Positive " + tagEta, kTH1D, {axisEta});
      histos.add("Data/neg/eta/hmpid", "HMPID Negative " + tagEta, kTH1D, {axisEta});

      histos.add("Data/pos/phi/hmpid", "HMPID Positive " + tagPhi, kTH1D, {axisPhi});
      histos.add("Data/neg/phi/hmpid", "HMPID Negative " + tagPhi, kTH1D, {axisPhi});

      histos.add("Data/pos/etaphi/hmpid", "HMPID Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
      histos.add("Data/neg/etaphi/hmpid", "HMPID Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    }

    listEfficiencyData.setObject(new THashList);
    if (makeEff) {
      LOG(debug) << "Making TEfficiency for Data";
      auto makeEfficiency = [&](TString effname, TString efftitle, auto templateHisto, TEfficiency*& eff) {
        TAxis* axis = histos.get<TH1>(templateHisto)->GetXaxis();
        if (axis->IsVariableBinSize()) {
          eff = new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray());
        } else {
          eff = new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
        }
        listEfficiencyData->Add(eff);
      };

      makeEfficiency("ITSTPCMatchingEfficiencyVsPt",
                     "ITS-TPC M.E. (ITS-TPC)/ITS in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"),
                     effITSTPCMatchingVsPt);
      makeEfficiency("TPCITSMatchingEfficiencyVsPt",
                     "ITS-TPC M.E. (ITS-TPC)/TPC in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"),
                     effTPCITSMatchingVsPt);
      makeEfficiency("TPCTOFMatchingEfficiencyVsPt",
                     "TPC-TOF M.E. in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"),
                     effTPCTOFMatchingVsPt);
      makeEfficiency("TPCTOFMatchingEfficiencyVsP",
                     "TPC-TOF M.E. in data " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"),
                     effTPCTOFMatchingVsP);
      makeEfficiency("TPCTOFMatchingEfficiencyVsEta",
                     "TPC-TOF M.E. in data " + tagEta + ";#it{#eta};Efficiency", HIST("Data/pos/eta/its_tpc_tof"),
                     effTPCTOFMatchingVsEta);
      makeEfficiency("TPCTOFMatchingEfficiencyVsPhi",
                     "TPC-TOF M.E. in data " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST("Data/pos/phi/its_tpc_tof"),
                     effTPCTOFMatchingVsPhi);

      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHistoX, auto templateHistoY, TEfficiency*& eff) {
        TAxis* axisX = histos.get<TH1>(templateHistoX)->GetXaxis();
        TAxis* axisY = histos.get<TH1>(templateHistoY)->GetYaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          eff = new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray());
        } else {
          eff = new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax());
        }
        listEfficiencyData->Add(eff);
      };

      makeEfficiency2D("TPCTOFMatchingEfficiencyVsPtVsEta",
                       Form("TPC-TOF M.E. in data #it{#varphi} [%.2f,%.2f];%s;%s;Efficiency", phiMin, phiMax, "#it{p}_{T} (GeV/#it{c})", "#it{#eta}"), HIST("Data/pos/pt/its_tpc_tof"), HIST("Data/pos/eta/its_tpc_tof"),
                       effTPCTOFMatchingVsPtVsEta);
      makeEfficiency2D("TPCTOFMatchingEfficiencyVsPtVsPhi",
                       Form("TPC-TOF M.E. in data #it{#eta} [%.2f,%.2f];%s;%s;Efficiency", etaMin, etaMax, "#it{p}_{T} (GeV/#it{c})", "#it{#varphi} (rad)"), HIST("Data/pos/pt/its_tpc_tof"), HIST("Data/pos/phi/its_tpc_tof"),
                       effTPCTOFMatchingVsPtVsPhi);
    }
  }

  // Selection cuts defined from the binning
  double ptMin, ptMax;
  double etaMin, etaMax;
  double phiMin, phiMax;
  double yMin, yMax;

  void init(InitContext&)
  {

    // Printing configuration
    LOG(info) << "Printing configuration";
    LOG(info) << "Set noFakesHits to: " << (noFakesHits ? "true" : "false");
    LOG(info) << "Set doPositivePDG to: " << (doPositivePDG ? "true" : "false");
    LOG(info) << "Set doNegativePDG to: " << (doNegativePDG ? "true" : "false");
    LOG(info) << "Set doEl to: " << (doEl ? "true" : "false");
    LOG(info) << "Set doMu to: " << (doMu ? "true" : "false");
    LOG(info) << "Set doPi to: " << (doPi ? "true" : "false");
    LOG(info) << "Set doKa to: " << (doKa ? "true" : "false");
    LOG(info) << "Set doPr to: " << (doPr ? "true" : "false");
    LOG(info) << "Set doDe to: " << (doDe ? "true" : "false");
    LOG(info) << "Set doTr to: " << (doTr ? "true" : "false");
    LOG(info) << "Set doHe to: " << (doHe ? "true" : "false");
    LOG(info) << "Set doAl to: " << (doAl ? "true" : "false");
    LOG(info) << "Set trackSelection to: " << (trackSelection ? "true" : "false");
    LOG(info) << "Set makeEff to: " << (makeEff ? "true" : "false");
    LOG(info) << "Set doPtEta to: " << (doPtEta ? "true" : "false");

    auto doLimits = [&](double& min, double& max, const ConfigurableAxis& binning) {
      const AxisSpec a{binning, "dummy"};
      min = a.binEdges[0];
      max = a.binEdges[a.binEdges.size() - 1];
      LOG(info) << "Making limits from " << min << ", " << max << " size " << a.getNbins();
    };

    doLimits(ptMin, ptMax, ptBins);
    doLimits(etaMin, etaMax, etaBins);
    doLimits(phiMin, phiMax, phiBins);
    doLimits(yMin, yMax, yBins);

    histos.add("eventSelection", "Event Selection", kTH1D, {{10, 0.5, 10.5, "Selection"}});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    if (applyEvSel == 0) {
      histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel. (no ev. sel)");
    } else if (applyEvSel == 1) {
      histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel. (sel7)");
    } else if (applyEvSel == 2) {
      histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel. (sel8)");
    } else {
      LOG(fatal) << "Can't interpret event selection asked " << applyEvSel << " (0: no event selection, 1: sel7, 2: sel8)";
    }

    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Position");

    if (doprocessMC) {
      histos.add("MC/generatedCollisions", "Generated Collisions", kTH1D, {{10, 0.5, 10.5, "Generated collisions"}});
      histos.get<TH1>(HIST("MC/generatedCollisions"))->GetXaxis()->SetBinLabel(1, "Gen. coll");
      histos.get<TH1>(HIST("MC/generatedCollisions"))->GetXaxis()->SetBinLabel(2, "At least 1 reco");
      histos.get<TH1>(HIST("MC/generatedCollisions"))->GetXaxis()->SetBinLabel(3, "At least 1 TPC track");
      histos.get<TH1>(HIST("MC/generatedCollisions"))->GetXaxis()->SetBinLabel(4, "Reco. coll.");
      histos.get<TH1>(HIST("MC/generatedCollisions"))->GetXaxis()->SetBinLabel(5, "Reco. good coll.");
    }

    const AxisSpec axisSel{40, 0.5, 40.5, "Selection"};
    initData(axisSel);
    initMC(axisSel);

    // Custom track cuts
    if (globalTrackSelection.value == 6) {
      customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(cfgCustomTrackCuts.itsPattern);
      LOG(info) << "Customizing track cuts:";
      if (cfgCustomTrackCuts.tracksIU.value) {
        customTrackCuts.SetTrackType(o2::aod::track::TrackTypeEnum::TrackIU);
      }
      customTrackCuts.SetRequireITSRefit(cfgCustomTrackCuts.requireITS);
      customTrackCuts.SetRequireTPCRefit(cfgCustomTrackCuts.requireTPC);
      customTrackCuts.SetRequireGoldenChi2(cfgCustomTrackCuts.requireGoldenChi2);
      customTrackCuts.SetRequireHitsInITSLayers(cfgCustomTrackCuts.minITScl.value, {0, 1, 2, 3, 4, 5, 6});
      customTrackCuts.SetMaxChi2PerClusterTPC(cfgCustomTrackCuts.maxChi2PerClusterTPC);
      customTrackCuts.SetMaxChi2PerClusterITS(cfgCustomTrackCuts.maxChi2PerClusterITS);
      customTrackCuts.SetMinNCrossedRowsTPC(cfgCustomTrackCuts.minNCrossedRowsTPC);
      customTrackCuts.SetMinNClustersTPC(cfgCustomTrackCuts.minTPCNClsFound);
      customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(cfgCustomTrackCuts.minNCrossedRowsOverFindableClustersTPC);
      customTrackCuts.SetMaxDcaXYPtDep([&](float /*pt*/) { return cfgCustomTrackCuts.maxDcaXY; }); // No DCAxy cut will be used, this is done via the member function of the task
      customTrackCuts.SetMaxDcaZ(cfgCustomTrackCuts.maxDcaZ);
      customTrackCuts.print();
    }
  }

  template <int pdgSign, o2::track::PID::ID id, typename particleType>
  bool isPdgSelected(const particleType& mcParticle)
  {
    static_assert(pdgSign == 0 || pdgSign == 1);
    static_assert(id > 0 || id < nSpecies);
    constexpr int index = id + pdgSign * nSpecies;
    return mcParticle.pdgCode() == PDGs[index];
  }

  bool isPhysicalPrimary(const o2::aod::McParticles::iterator& mcParticle)
  {
    if (maxProdRadius < 999.f) {
      if ((mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()) > maxProdRadius * maxProdRadius) {
        return false;
      }
    }
    return mcParticle.isPhysicalPrimary();
  }
  bool isFinal(const o2::aod::McParticles::iterator& mcParticle)
  {
    if (!mcParticle.has_daughters() && !mcParticle.isPhysicalPrimary() && mcParticle.getProcess() == 4) {
      auto mothers = mcParticle.mothers_as<o2::aod::McParticles>();
      for (const auto& mother : mothers) {
        if (!mother.isPhysicalPrimary() && mother.getProcess() == 4) {
          return true;
        }
      }
    }
    return false; // Otherwise, not considered a tertiary particle
  }
  template <int pdgSign, o2::track::PID::ID id>
  void fillMCTrackHistograms(const TrackCandidatesMC::iterator& track, const bool doMakeHistograms)
  {
    static_assert(pdgSign == 0 || pdgSign == 1);
    if (!doMakeHistograms) {
      return;
    }

    if constexpr (pdgSign == 0) {
      if (!doPositivePDG) {
        return;
      }
    } else {
      if (!doNegativePDG) {
        return;
      }
    }
    constexpr int histogramIndex = id + pdgSign * nSpecies;
    LOG(debug) << "fillMCTrackHistograms for pdgSign '" << pdgSign << "' and id '" << static_cast<int>(id) << "' " << particleName(pdgSign, id) << " with index " << histogramIndex;
    const o2::aod::McParticles::iterator& mcParticle = track.mcParticle();
    const CollisionCandidatesMC::iterator& collision = track.collision_as<CollisionCandidatesMC>();
    float radius = std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy());
    if (numSameCollision) {
      if (!collision.has_mcCollision()) {
        return;
      }
      if (mcParticle.mcCollision().globalIndex() != collision.mcCollision().globalIndex()) {
        return;
      }
    }
    histos.fill(HIST("MC/trackSelection"), trkCutSameColl);

    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("MC/trackSelection"), trkCutIdxN + id);

    if (passedITS) {
      hPtIts[histogramIndex]->Fill(mcParticle.pt());
    }
    if (passedTPC) {
      hPtTpc[histogramIndex]->Fill(mcParticle.pt());
    }
    if (passedITS && passedTPC) {
      hPItsTpc[histogramIndex]->Fill(mcParticle.p());
      hPtItsTpc[histogramIndex]->Fill(mcParticle.pt());
      hEtaItsTpc[histogramIndex]->Fill(mcParticle.eta());
      hYItsTpc[histogramIndex]->Fill(mcParticle.y());
      hPhiItsTpc[histogramIndex]->Fill(mcParticle.phi());

      hPTrkItsTpc[histogramIndex]->Fill(track.p());
      hPtTrkItsTpc[histogramIndex]->Fill(track.pt());
      hEtaTrkItsTpc[histogramIndex]->Fill(track.eta());
      hPhiTrkItsTpc[histogramIndex]->Fill(track.phi());

      if (doPtEta) {
        hPtEtaItsTpc[histogramIndex]->Fill(mcParticle.pt(), mcParticle.eta());
        hPtEtaTrkItsTpc[histogramIndex]->Fill(track.pt(), track.eta());
        if (passedTOF) {
          hPtEtaItsTpcTof[histogramIndex]->Fill(mcParticle.pt(), mcParticle.eta());
        }
      }
      if (doPtRadius) {
        hPtRadiusItsTpc[histogramIndex]->Fill(mcParticle.pt(), radius);
        if (passedTOF) {
          hPtRadiusItsTpcTof[histogramIndex]->Fill(mcParticle.pt(), radius);
        }
      }
    }
    if (passedITS && passedTOF) {
      hPtItsTof[histogramIndex]->Fill(mcParticle.pt());
    }
    if (passedTPC && passedTOF) {
      hPtTpcTof[histogramIndex]->Fill(mcParticle.pt());
    }
    if (passedITS && passedTPC && passedTRD) {
      hPtItsTpcTrd[histogramIndex]->Fill(mcParticle.p());
    }
    if (passedITS && passedTPC && passedTOF) {
      hPItsTpcTof[histogramIndex]->Fill(mcParticle.p());
      hPtItsTpcTof[histogramIndex]->Fill(mcParticle.pt());
      hEtaItsTpcTof[histogramIndex]->Fill(mcParticle.eta());
      hYItsTpcTof[histogramIndex]->Fill(mcParticle.y());
      hPhiItsTpcTof[histogramIndex]->Fill(mcParticle.phi());
    }
    if (passedITS && passedTPC && passedTRD && passedTOF) {
      hPtItsTpcTrdTof[histogramIndex]->Fill(mcParticle.p());
    }
    if (isPhysicalPrimary(mcParticle)) {
      if (passedITS) {
        hPtItsPrm[histogramIndex]->Fill(mcParticle.pt());
      }
      if (passedITS && passedTPC) {
        hPtItsTpcPrm[histogramIndex]->Fill(mcParticle.pt());
        hPtTrkItsTpcPrm[histogramIndex]->Fill(track.pt());
        hEtaItsTpcPrm[histogramIndex]->Fill(mcParticle.eta());
        hEtaTrkItsTpcPrm[histogramIndex]->Fill(track.eta());
        hPhiItsTpcPrm[histogramIndex]->Fill(mcParticle.phi());
        hPhiTrkItsTpcPrm[histogramIndex]->Fill(track.phi());
        if (passedTOF) {
          hPtItsTpcTofPrm[histogramIndex]->Fill(mcParticle.pt());
          hPtTrkItsTpcTofPrm[histogramIndex]->Fill(track.pt());
          hEtaItsTpcTofPrm[histogramIndex]->Fill(mcParticle.eta());
          hPhiItsTpcTofPrm[histogramIndex]->Fill(mcParticle.phi());
        }
        if (doPtRadius) {
          hPtRadiusItsTpcPrm[histogramIndex]->Fill(mcParticle.pt(), radius);
          if (passedTOF) {
            hPtRadiusItsTpcTofPrm[histogramIndex]->Fill(mcParticle.pt(), radius);
          }
        }
      }
    } else if (mcParticle.getProcess() == 4) { // Particle decay
      // Checking mothers
      bool motherIsAccepted = true;
      if (checkForMothers.value && mothersPDGs.value.size() > 0 && mcParticle.has_mothers()) {
        motherIsAccepted = false;
        auto mothers = mcParticle.mothers_as<o2::aod::McParticles>();
        for (const auto& mother : mothers) {
          for (const auto& pdgToCheck : mothersPDGs.value) {
            if (mother.pdgCode() == pdgToCheck) {
              motherIsAccepted = true;
              // Calculate the decay length
              double decayLength = std::sqrt(std::pow(mother.vx() - mother.mcCollision().posX(), 2) + std::pow(mother.vy() - mother.mcCollision().posY(), 2) + std::pow(mother.vz() - mother.mcCollision().posZ(), 2));
              hdecaylengthmother[histogramIndex]->Fill(decayLength);
              break;
            }
            if (motherIsAccepted) {
              break;
            }
          }
        }
      }
      if (passedITS && passedTPC && motherIsAccepted) {
        hPtItsTpcStr[histogramIndex]->Fill(mcParticle.pt());
        hPtTrkItsTpcStr[histogramIndex]->Fill(track.pt());
        if (passedTOF) {
          hPtItsTpcTofStr[histogramIndex]->Fill(mcParticle.pt());
        }
        if (doPtRadius) {
          hPtRadiusItsTpcStr[histogramIndex]->Fill(mcParticle.pt(), radius);
          if (passedTOF) {
            hPtRadiusItsTpcTofStr[histogramIndex]->Fill(mcParticle.pt(), radius);
          }
        }
      }
      if (isFinal(mcParticle)) {
        if (passedITS && passedTPC && motherIsAccepted) {
          hPtItsTpcTer[histogramIndex]->Fill(mcParticle.pt());
          hPtTrkItsTpcTer[histogramIndex]->Fill(track.pt());
          if (passedTOF) {
            hPtItsTpcTofTer[histogramIndex]->Fill(mcParticle.pt());
          }
          if (doPtRadius) {
            hPtRadiusItsTpcTer[histogramIndex]->Fill(mcParticle.pt(), radius);
            if (passedTOF) {
              hPtRadiusItsTpcTofTer[histogramIndex]->Fill(mcParticle.pt(), radius);
            }
          }
        }
      }
    } else { // Material
      if (passedITS && passedTPC) {
        hPtItsTpcMat[histogramIndex]->Fill(mcParticle.pt());
        hPtTrkItsTpcMat[histogramIndex]->Fill(track.pt());
        if (passedTOF) {
          hPtItsTpcTofMat[histogramIndex]->Fill(mcParticle.pt());
        }
      }
    }
  }

  template <int pdgSign, o2::track::PID::ID id, bool recoEv = false>
  void fillMCParticleHistograms(const o2::aod::McParticles::iterator& mcParticle, const bool doMakeHistograms)
  {
    static_assert(pdgSign == 0 || pdgSign == 1);
    if (!doMakeHistograms) {
      return;
    }

    if constexpr (pdgSign == 0) {
      if (!doPositivePDG) {
        return;
      }
    } else {
      if (!doNegativePDG) {
        return;
      }
    }

    constexpr int histogramIndex = id + pdgSign * nSpecies;
    LOG(debug) << "fillMCParticleHistograms for pdgSign '" << pdgSign << "' and id '" << static_cast<int>(id) << "' " << particleName(pdgSign, id) << " with index " << histogramIndex;
    float radius = std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy());
    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }
    if constexpr (recoEv) {
      hPtGeneratedRecoEv[histogramIndex]->Fill(mcParticle.pt());
      if (isPhysicalPrimary(mcParticle)) {
        hPtGeneratedPrmRecoEv[histogramIndex]->Fill(mcParticle.pt());
      }
      return;
    }
    histos.fill(HIST("MC/particleSelection"), 6 + id);

    hPGenerated[histogramIndex]->Fill(mcParticle.p());
    hPtGenerated[histogramIndex]->Fill(mcParticle.pt());

    if (isPhysicalPrimary(mcParticle)) {
      hPtGeneratedPrm[histogramIndex]->Fill(mcParticle.pt());
      hEtaGeneratedPrm[histogramIndex]->Fill(mcParticle.eta());
      hPhiGeneratedPrm[histogramIndex]->Fill(mcParticle.phi());
      if (doPtRadius) {
        hPtRadiusGeneratedPrm[histogramIndex]->Fill(mcParticle.pt(), radius);
      }
    } else {
      if (mcParticle.getProcess() == 4) { // Particle decay
        bool motherIsAccepted = true;
        // Check for mothers if needed
        if (checkForMothers.value && mothersPDGs.value.size() > 0 && mcParticle.has_mothers()) {
          motherIsAccepted = false;
          auto mothers = mcParticle.mothers_as<o2::aod::McParticles>();
          // Loop over mother particles
          for (const auto& mother : mothers) {
            for (const auto& pdgToCheck : mothersPDGs.value) {
              if (mother.pdgCode() == pdgToCheck) {
                motherIsAccepted = true; // Mother matches the list of specified PDGs
                hPtmotherGenerated[histogramIndex]->Fill(mother.pt()); // Fill generated pT for mother
                break;
              }
              if (motherIsAccepted) {
                break;
              }
            }
          }
        }
        if (motherIsAccepted) {
          hPtGeneratedStr[histogramIndex]->Fill(mcParticle.pt());
          if (doPtRadius) {
            hPtRadiusGeneratedStr[histogramIndex]->Fill(mcParticle.pt(), radius);
          }
          if (isFinal(mcParticle)) {
            hPtGeneratedTer[histogramIndex]->Fill(mcParticle.pt());
            if (doPtRadius) {
              hPtRadiusGeneratedTer[histogramIndex]->Fill(mcParticle.pt(), radius);
            }
          }
        }
      } else { // Material
        hPtGeneratedMat[histogramIndex]->Fill(mcParticle.pt());
      }
    }
    hEtaGenerated[histogramIndex]->Fill(mcParticle.eta());
    hYGenerated[histogramIndex]->Fill(mcParticle.y());
    hPhiGenerated[histogramIndex]->Fill(mcParticle.phi());
    if (doPtEta) {
      hPtEtaGenerated[histogramIndex]->Fill(mcParticle.pt(), mcParticle.eta());
    }
    if (doPtRadius) {
      hPtRadiusGenerated[histogramIndex]->Fill(mcParticle.pt(), radius);
    }
  }

  template <int pdgSign, o2::track::PID::ID id>
  void fillMCEfficiency(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    static_assert(pdgSign == 0 || pdgSign == 1);
    if constexpr (pdgSign == 0) {
      if (!doPositivePDG) {
        return;
      }
    } else {
      if (!doNegativePDG) {
        return;
      }
    }
    if (!makeEff) {
      return;
    }

    constexpr int histogramIndex = id + pdgSign * nSpecies;

    const char* partName = particleName(pdgSign, id);
    LOG(debug) << "Filling efficiency for particle " << static_cast<int>(id) << " " << partName;
    THashList* subList = static_cast<THashList*>(listEfficiencyMC->FindObject(partName));
    if (!subList) {
      LOG(warning) << "Cannot find list of efficiency objects for particle " << partName;
      return;
    }

    // Filling 1D efficiencies
    auto doFillEfficiency = [&](const TString effname, auto num, auto den) {
      TEfficiency* eff = static_cast<TEfficiency*>(subList->FindObject(effname));
      if (!eff) {
        LOG(warning) << "Cannot find TEfficiency " << effname;
        return;
      }
      eff->SetTotalHistogram(*den, "f");
      eff->SetPassedHistogram(*num, "f");
    };

    doFillEfficiency("ITS_vsPt", hPtIts[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("TPC_vsPt", hPtTpc[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt", hPtItsTpc[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("ITS-TOF_vsPt", hPtItsTof[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("TPC-TOF_vsPt", hPtTpcTof[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TRD_vsPt", hPtItsTpcTrd[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPt", hPtItsTpcTof[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TRD-TOF_vsPt", hPtItsTpcTrdTof[histogramIndex], hPtGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt_Trk", hPtTrkItsTpc[histogramIndex], hPtGenerated[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPt_RecoEv", hPtItsTpcPrm[histogramIndex], hPtGeneratedRecoEv[histogramIndex]);

    doFillEfficiency("ITS_vsPt_Prm", hPtItsPrm[histogramIndex], hPtGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt_Prm", hPtItsTpcPrm[histogramIndex], hPtGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt_Prm_Trk", hPtTrkItsTpcPrm[histogramIndex], hPtGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPt_Prm", hPtItsTpcTofPrm[histogramIndex], hPtGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPt_Prm_Trk", hPtTrkItsTpcTofPrm[histogramIndex], hPtGeneratedPrm[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPt_Prm_RecoEv", hPtItsTpcPrm[histogramIndex], hPtGeneratedPrmRecoEv[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPt_Str", hPtItsTpcStr[histogramIndex], hPtGeneratedStr[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt_Str_Trk", hPtTrkItsTpcStr[histogramIndex], hPtGeneratedStr[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPt_Str", hPtItsTpcTofStr[histogramIndex], hPtGeneratedStr[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPt_Mat", hPtItsTpcMat[histogramIndex], hPtGeneratedMat[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt_Mat_Trk", hPtTrkItsTpcMat[histogramIndex], hPtGeneratedMat[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPt_Mat", hPtItsTpcTofMat[histogramIndex], hPtGeneratedMat[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPt_Ter", hPtItsTpcTer[histogramIndex], hPtGeneratedTer[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPt_Ter_Trk", hPtTrkItsTpcTer[histogramIndex], hPtGeneratedTer[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPt_Ter", hPtItsTpcTofTer[histogramIndex], hPtGeneratedTer[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsP", hPItsTpc[histogramIndex], hPGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsP_Trk", hPTrkItsTpc[histogramIndex], hPGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsP", hPItsTpcTof[histogramIndex], hPGenerated[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsEta", hEtaItsTpc[histogramIndex], hEtaGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsEta_Trk", hEtaTrkItsTpc[histogramIndex], hEtaGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsEta", hEtaItsTpcTof[histogramIndex], hEtaGenerated[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsEta_Prm", hEtaItsTpcPrm[histogramIndex], hEtaGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsEta_Prm_Trk", hEtaTrkItsTpcPrm[histogramIndex], hEtaGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsEta_Prm", hEtaItsTpcTofPrm[histogramIndex], hEtaGeneratedPrm[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPhi_Prm", hPhiItsTpcPrm[histogramIndex], hPhiGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPhi_Prm_Trk", hPhiTrkItsTpcPrm[histogramIndex], hPhiGeneratedPrm[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPhi_Prm", hPhiItsTpcTofPrm[histogramIndex], hPhiGeneratedPrm[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsY", hYItsTpc[histogramIndex], hYGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsY", hYItsTpcTof[histogramIndex], hYGenerated[histogramIndex]);

    doFillEfficiency("ITS-TPC_vsPhi", hPhiItsTpc[histogramIndex], hPhiGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC_vsPhi_Trk", hPhiTrkItsTpc[histogramIndex], hPhiGenerated[histogramIndex]);
    doFillEfficiency("ITS-TPC-TOF_vsPhi", hPhiItsTpcTof[histogramIndex], hPhiGenerated[histogramIndex]);

    if (!doPtEta) {
      return;
    }

    // Filling 2D efficiencies
    auto fillEfficiency2D = [&](const TString effname, auto num, auto den) {
      TEfficiency* eff = static_cast<TEfficiency*>(subList->FindObject(effname));
      if (!eff) {
        LOG(warning) << "Cannot find TEfficiency " << effname;
        return;
      }
      eff->SetTotalHistogram(*den, "f");
      eff->SetPassedHistogram(*num, "f");
    };
    fillEfficiency2D("ITS-TPC_vsPt_vsEta", hPtEtaItsTpc[histogramIndex], hPtEtaGenerated[histogramIndex]);
    fillEfficiency2D("ITS-TPC_vsPt_vsEta_Trk", hPtEtaTrkItsTpc[histogramIndex], hPtEtaGenerated[histogramIndex]);
    fillEfficiency2D("ITS-TPC-TOF_vsPt_vsEta", hPtEtaItsTpcTof[histogramIndex], hPtEtaGenerated[histogramIndex]);
    if (!doPtRadius) {
      return;
    }
    fillEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpc[histogramIndex], hPtRadiusGenerated[histogramIndex]);
    fillEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTof[histogramIndex], hPtRadiusGenerated[histogramIndex]);
    fillEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpcPrm[histogramIndex], hPtRadiusGeneratedPrm[histogramIndex]);
    fillEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTofPrm[histogramIndex], hPtRadiusGeneratedPrm[histogramIndex]);
    fillEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpcStr[histogramIndex], hPtRadiusGeneratedStr[histogramIndex]);
    fillEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTofStr[histogramIndex], hPtRadiusGeneratedStr[histogramIndex]);
    fillEfficiency2D("ITS-TPC_vsPt_vsRadius", hPtRadiusItsTpcTer[histogramIndex], hPtRadiusGeneratedTer[histogramIndex]);
    fillEfficiency2D("ITS-TPC-TOF_vsPt_vsRadius", hPtRadiusItsTpcTofTer[histogramIndex], hPtRadiusGeneratedTer[histogramIndex]);
  }
  template <bool doFillHistograms, typename CollType>
  bool isCollisionSelected(const CollType& collision)
  {
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 1);
    }
    if (applyEvSel == 1 && !collision.sel7()) {
      return false;
    } else if (applyEvSel == 2 && !collision.sel8()) {
      return false;
    }
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 2);
    }
    if (collision.numContrib() < nMinNumberOfContributors) {
      return false;
    }
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 3);
    }
    if ((collision.posZ() < vertexZMin || collision.posZ() > vertexZMax)) {
      return false;
    }
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 4);
    }
    return true;
  }

  // Global process
  void process(CollisionCandidates::iterator const& collision) { isCollisionSelected<true>(collision); }

  // Function to apply particle selection
  template <bool isMC = true, bool doFillHisto = true, typename particleType, typename histoType = int>
  bool isInAcceptance(const particleType& particle, const histoType& countingHisto = 0, const int offset = 0)
  {
    if (particle.pt() < ptMin || particle.pt() > ptMax) { // Check pt
      return false;
    }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, 1 + offset);
    }
    if (particle.eta() < etaMin || particle.eta() > etaMax) { // Check eta
      return false;
    }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, 2 + offset);
    }
    if (particle.phi() < phiMin || particle.phi() > phiMax) { // Check phi
      return false;
    }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, 3 + offset);
    }
    if constexpr (isMC) {
      if (particle.y() < yMin || particle.y() > yMax) { // Check rapidity
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, 4 + offset);
      }
    }

    return true;
  }

  // Function to apply track selection
  bool passedITS = false;
  bool passedTPC = false;
  bool passedTRD = false;
  bool passedTOF = false;
  template <bool isMC = true, bool doFillHisto = true, typename trackType, typename histoType = int>
  bool isTrackSelected(trackType& track, const histoType& countingHisto = 0)
  {
    // Reset selections
    passedITS = false;
    passedTPC = false;
    passedTRD = false;
    passedTOF = false;

    if constexpr (doFillHisto) {
      histos.fill(countingHisto, trkCutIdxTrkRead); // Read tracks
    }

    if constexpr (isMC) { // MC only
      if (!track.has_mcParticle()) {
        histos.fill(HIST("MC/fakeTrackNoiseHits"), 0.5);
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxHasMcPart); // Tracks with particles (i.e. no fakes)
      }
      const auto mcParticle = track.mcParticle();
      if (!isInAcceptance<true, doFillHisto>(mcParticle, countingHisto, trkCutIdxHasMcPart)) {
        // 3: pt cut 4: eta cut 5: phi cut 6: y cut
        return false;
      }

      if (noFakesHits) { // Selecting tracks with no fake hits
        bool hasFakeHit = false;
        for (int i = 0; i < 10; i++) { // From ITS to TPC
          if (track.mcMask() & 1 << i) {
            hasFakeHit = true;
            break;
          }
        }
        if (hasFakeHit) {
          return false;
        }
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedFake);
      }
    } else { // Data only
      if (!isInAcceptance<false, doFillHisto>(track, countingHisto, trkCutIdxHasMcPart)) {
        // 3: pt cut 4: eta cut 5: phi cut 6: y cut
        return false;
      }
    }

    if (!track.has_collision()) {
      return false;
    }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, trkCutIdxHasCollision);
    }

    if (trackSelection) { // Check general cuts
      if (!track.passedTrackType()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedTrkType);
      }
      if (!track.passedPtRange()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedPtRange);
      }
      if (!track.passedEtaRange()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedEtaRange);
      }
      if (!track.passedDCAxy()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedDcaXYMax);
      }
      if (std::abs(track.dcaXY()) < minDcaXY) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedDcaXYMin);
      }
      if (!track.passedDCAz()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedDcaZMax);
      }
      if (std::abs(track.dcaZ()) < minDcaZ) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedDcaZMin);
      }
      if (!track.passedGoldenChi2()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedGoldenChi2);
      }
      if (doPVContributorCut && !track.isPVContributor()) {
        return false;
      }
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedIsPvCont);
      }

      passedITS = track.passedITSNCls() &&
                  track.passedITSChi2NDF() &&
                  track.passedITSRefit() &&
                  track.passedITSHits() &&
                  track.hasITS();

      passedTPC = track.passedTPCNCls() &&
                  track.passedTPCCrossedRows() &&
                  track.passedTPCCrossedRowsOverNCls() &&
                  track.passedTPCChi2NDF() &&
                  track.passedTPCRefit() &&
                  track.hasTPC();
      passedTRD = track.hasTRD();
      passedTOF = track.hasTOF();
    } else {
      passedITS = track.hasITS();
      passedTPC = track.hasTPC();
      passedTRD = track.hasTRD();
      passedTOF = track.hasTOF();
    }

    if (passedITS) { // Partial
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedITSPartial);
      }
    }

    if (passedTPC) { // Partial
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedTPCPartial);
      }
    }

    if (passedTOF) { // Partial
      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxPassedTOFPartial);
      }
    }
    bool isTrackSelectedAfteAll = false;
    switch (globalTrackSelection) {
      case 0:
        isTrackSelectedAfteAll = true;
        break;
      case 1:
        isTrackSelectedAfteAll = track.isGlobalTrack();
        break;
      case 2:
        isTrackSelectedAfteAll = track.isGlobalTrackWoPtEta();
        break;
      case 3:
        isTrackSelectedAfteAll = track.isGlobalTrackWoDCA();
        break;
      case 4:
        isTrackSelectedAfteAll = track.isQualityTrack();
        break;
      case 5:
        isTrackSelectedAfteAll = track.isInAcceptanceTrack();
        break;
      case 6:
        isTrackSelectedAfteAll = customTrackCuts.IsSelected(track);
        break;
      default:
        LOG(fatal) << "Can't interpret track asked selection " << globalTrackSelection;
    }
    if (!isTrackSelectedAfteAll) {
      return false;
    }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, trkCutIdxPassedGlobal);
    }

    return true;
  }

  /// \brief Function to get MC collision occupancy
  /// \param collSlice collection of reconstructed collisions
  /// \return collision occupancy
  template <typename CCs>
  int getOccupancyColl(CCs const& collSlice)
  {
    float multiplicity{0.f};
    int occupancy = 0;
    for (const auto& collision : collSlice) {
      float collMult{0.f};
      collMult = collision.numContrib();

      if (collMult > multiplicity) {
        if (useFT0OccEstimator) {
          /// occupancy estimator (FT0c signal amplitudes in +-10us from current collision)
          occupancy = collision.ft0cOccupancyInTimeRange();
        } else {
          /// occupancy estimator (ITS tracks with at least 5 clusters in +-10us from current collision)
          occupancy = static_cast<float>(collision.trackOccupancyInTimeRange());
        }
        multiplicity = collMult;
      }
    } // end loop over collisions

    return occupancy;
  }

  /// \brief Function to get MC collision centrality
  /// \param collSlice collection of reconstructed collisions
  /// \return collision centrality
  template <typename CCs>
  int getCentralityColl(CCs const& collSlice)
  {
    float multiplicity{0.f};
    int centrality = 0;
    for (const auto& collision : collSlice) {
      float collMult{0.f};
      collMult = collision.numContrib();

      if (collMult > multiplicity) {
        centrality = collision.centFT0C();
        multiplicity = collMult;
      }
    } // end loop over collisions

    return centrality;
  }

  // MC process
  // Single-track efficiency calculated only for MC collisions with at least 1 reco. collision
  SliceCache cache;
  Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  Preslice<o2::aod::McParticles> perCollisionMc = o2::aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CollisionCandidatesMC> collPerCollMc = o2::aod::mccollisionlabel::mcCollisionId;
  void processMC(o2::aod::McCollisions const& mcCollisions,
                 // o2::soa::SmallGroups<CollisionCandidatesMC> const& collisions,
                 CollisionCandidatesMC const& collisions,
                 TrackCandidatesMC const& tracks,
                 o2::aod::McParticles const& mcParticles)
  {

    /// loop over generated collisions
    for (const auto& mcCollision : mcCollisions) {
      histos.fill(HIST("MC/generatedCollisions"), 1);

      const auto groupedCollisions = collisions.sliceBy(collPerCollMc, mcCollision.globalIndex());
      const auto groupedMcParticles = mcParticles.sliceBy(perCollisionMc, mcCollision.globalIndex());

      // LOG(info) << "groupedCollisions.size() " << groupedCollisions.size();

      if (groupedCollisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
        continue;
      }
      float centrality = -1.;
      float occupancy = -1.;
      if (doOccupancy) {
        centrality = getCentralityColl(groupedCollisions);
        occupancy = getOccupancyColl(groupedCollisions);
      }
      histos.fill(HIST("MC/generatedCollisions"), 2);
      if (skipEventsWithoutTPCTracks) {
        int nTPCTracks = 0;
        for (const auto& collision : groupedCollisions) {
          const auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
          for (const auto& track : groupedTracks) {
            if (track.hasTPC()) {
              nTPCTracks++;
              break;
            }
          }
        }
        if (nTPCTracks == 0) {
          LOG(info) << "Skipping event with no TPC tracks";
          continue;
        }
      }
      histos.fill(HIST("MC/generatedCollisions"), 3);

      if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
        LOG(debug) << "Skipping event with different type of generator than the one requested";
        continue;
      }

      /// loop over reconstructed collisions
      for (const auto& collision : groupedCollisions) {
        histos.fill(HIST("MC/generatedCollisions"), 4);
        if (!isCollisionSelected<false>(collision)) {
          continue;
        }
        histos.fill(HIST("MC/generatedCollisions"), 5);
        if (useFT0OccEstimator) {
          /// occupancy estimator (FT0c signal amplitudes in +-10us from current collision)
          occupancy = collision.ft0cOccupancyInTimeRange();
        } else {
          /// occupancy estimator (ITS tracks with at least 5 clusters in +-10us from current collision)
          occupancy = static_cast<float>(collision.trackOccupancyInTimeRange());
        }
        centrality = collision.centFT0C();

        const auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());

        // Track loop
        for (const auto& track : groupedTracks) {
          if (!isTrackSelected(track, HIST("MC/trackSelection"))) {
            continue;
          }

          // search for particles from HF decays
          // no need to check if track.has_mcParticle() == true, this is done already in isTrackSelected
          const auto& particle = track.mcParticle();
          if (keepOnlyHfParticles && !RecoDecay::getCharmHadronOrigin(mcParticles, particle, /*searchUpToQuark*/ true)) {
            continue;
          }

          if (doOccupancy) {
            float trackPt = track.pt();
            float trackSign = track.sign();
            if (trackSign > 0) {
              if (passedTOF && passedTPC && passedITS) {
                histos.fill(HIST("MC/occ_cent/reco/pos/its_tpc_tof"), occupancy, centrality, trackPt);
              }
              if (passedTPC && passedITS) {
                histos.fill(HIST("MC/occ_cent/reco/pos/its_tpc"), occupancy, centrality, trackPt);
              }
              if (passedITS) {
                histos.fill(HIST("MC/occ_cent/reco/pos/its"), occupancy, centrality, trackPt);
              }
            }
            if (trackSign < 0) {
              if (passedTOF && passedTPC && passedITS) {
                histos.fill(HIST("MC/occ_cent/reco/neg/its_tpc_tof"), occupancy, centrality, trackPt);
              }
              if (passedTPC && passedITS) {
                histos.fill(HIST("MC/occ_cent/reco/neg/its_tpc"), occupancy, centrality, trackPt);
              }
              if (passedITS) {
                histos.fill(HIST("MC/occ_cent/reco/neg/its"), occupancy, centrality, trackPt);
              }
            }
          }

          // Filling variable histograms
          histos.fill(HIST("MC/trackLength"), track.length());
          static_for<0, 1>([&](auto pdgSign) {
            fillMCTrackHistograms<pdgSign, o2::track::PID::Electron>(track, doEl);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Muon>(track, doMu);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Pion>(track, doPi);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Kaon>(track, doKa);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Proton>(track, doPr);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Deuteron>(track, doDe);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Triton>(track, doTr);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Helium3>(track, doHe);
            fillMCTrackHistograms<pdgSign, o2::track::PID::Alpha>(track, doAl);
          });
        }

        // Skipping collisions without the generated collisions
        // Actually this should never happen, since we group per MC collision
        if (!collision.has_mcCollision()) {
          continue;
        } else {
          // skip generated collisions outside the allowed vtx-z range
          // putting this condition here avoids the particle loop a few lines below
          if (applyPvZCutGenColl) {
            const float genPvZ = mcCollision.posZ();
            if (genPvZ < vertexZMin || genPvZ > vertexZMax) {
              continue;
            }
          }
        }

        /// only to fill denominator of ITS-TPC matched primary tracks only in MC events with at least 1 reco. vtx
        for (const auto& particle : groupedMcParticles) { // Particle loop

          /// require generated particle in acceptance
          if (!isInAcceptance<true, false>(particle, nullptr)) {
            continue;
          }

          // search for particles from HF decays
          if (keepOnlyHfParticles && !RecoDecay::getCharmHadronOrigin(mcParticles, particle, /*searchUpToQuark*/ true)) {
            continue;
          }

          if (doOccupancy) {
            float partSign = particle.pdgCode();
            float mcPartPt = particle.pt();
            if (partSign > 0) {
              histos.fill(HIST("MC/occ_cent/gen/pos"), occupancy, centrality, mcPartPt);
            }
            if (partSign < 0) {
              histos.fill(HIST("MC/occ_cent/gen/neg"), occupancy, centrality, mcPartPt);
            }
          }

          static_for<0, 1>([&](auto pdgSign) {
            fillMCParticleHistograms<pdgSign, o2::track::PID::Electron, true>(particle, doEl);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Muon, true>(particle, doMu);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Pion, true>(particle, doPi);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Kaon, true>(particle, doKa);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Proton, true>(particle, doPr);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Deuteron, true>(particle, doDe);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Triton, true>(particle, doTr);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Helium3, true>(particle, doHe);
            fillMCParticleHistograms<pdgSign, o2::track::PID::Alpha, true>(particle, doAl);
          });
        }
      } /// end loop over reconstructed collisions

      // skip generated collisions outside the allowed vtx-z range
      // putting this condition here avoids the particle loop a few lines below
      if (applyPvZCutGenColl) {
        const float genPvZ = mcCollision.posZ();
        if (genPvZ < vertexZMin || genPvZ > vertexZMax) {
          continue;
        }
      }

      // Loop on particles to fill the denominator
      float dNdEta = 0; // Multiplicity
      for (const auto& mcParticle : groupedMcParticles) {
        if (TMath::Abs(mcParticle.eta()) <= 2.f && !mcParticle.has_daughters()) {
          dNdEta += 1.f;
        }
        if (!isInAcceptance(mcParticle, HIST("MC/particleSelection"))) {
          continue;
        }

        // search for particles from HF decays
        if (keepOnlyHfParticles && !RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, /*searchUpToQuark*/ true)) {
          continue;
        }

        static_for<0, 1>([&](auto pdgSign) {
          fillMCParticleHistograms<pdgSign, o2::track::PID::Electron>(mcParticle, doEl);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Muon>(mcParticle, doMu);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Pion>(mcParticle, doPi);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Kaon>(mcParticle, doKa);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Proton>(mcParticle, doPr);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Deuteron>(mcParticle, doDe);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Triton>(mcParticle, doTr);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Helium3>(mcParticle, doHe);
          fillMCParticleHistograms<pdgSign, o2::track::PID::Alpha>(mcParticle, doAl);
        });
      }
      histos.fill(HIST("MC/eventMultiplicity"), dNdEta * 0.5f / 2.f);

      // Fill TEfficiencies
      static_for<0, 1>([&](auto pdgSign) {
        fillMCEfficiency<pdgSign, o2::track::PID::Electron>(doEl);
        fillMCEfficiency<pdgSign, o2::track::PID::Muon>(doMu);
        fillMCEfficiency<pdgSign, o2::track::PID::Pion>(doPi);
        fillMCEfficiency<pdgSign, o2::track::PID::Kaon>(doKa);
        fillMCEfficiency<pdgSign, o2::track::PID::Proton>(doPr);
        fillMCEfficiency<pdgSign, o2::track::PID::Deuteron>(doDe);
        fillMCEfficiency<pdgSign, o2::track::PID::Triton>(doTr);
        fillMCEfficiency<pdgSign, o2::track::PID::Helium3>(doHe);
        fillMCEfficiency<pdgSign, o2::track::PID::Alpha>(doAl);
      });
    } /// end loop over generated collisions
  }
  PROCESS_SWITCH(QaEfficiency, processMC, "process MC", false);

  // MC process without the collision association
  // Single-track efficiency calculated:
  //  - considering also MC collisions without any reco. collision
  //  - considering also tracks not associated to any collision
  //  - ignoring the track-to-collision association
  void processMCWithoutCollisions(TrackCandidatesMC const& tracks,
                                  o2::aod::Collisions const&,
                                  o2::aod::McParticles const& mcParticles,
                                  o2::aod::McCollisions const&)
  {
    // Track loop
    for (const auto& track : tracks) {
      if (!isTrackSelected(track, HIST("MC/trackSelection"))) {
        continue;
      }

      /// checking the PV z coordinate, if the track has been assigned to any collision
      if (applyPvZCutInProcessMcWoColl && track.has_collision()) {
        const auto collision = track.collision();
        const float posZ = collision.posZ();
        if (posZ < vertexZMin || posZ > vertexZMax) {
          continue;
        }
      }

      // search for particles from HF decays
      // no need to check if track.has_mcParticle() == true, this is done already in isTrackSelected
      const auto& particle = track.mcParticle();
      if (keepOnlyHfParticles && !RecoDecay::getCharmHadronOrigin(mcParticles, particle, /*searchUpToQuark*/ true)) {
        continue;
      }

      // Filling variable histograms
      histos.fill(HIST("MC/trackLength"), track.length());
      static_for<0, 1>([&](auto pdgSign) {
        fillMCTrackHistograms<pdgSign, o2::track::PID::Electron>(track, doEl);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Muon>(track, doMu);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Pion>(track, doPi);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Kaon>(track, doKa);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Proton>(track, doPr);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Deuteron>(track, doDe);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Triton>(track, doTr);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Helium3>(track, doHe);
        fillMCTrackHistograms<pdgSign, o2::track::PID::Alpha>(track, doAl);
      });
    }

    for (const auto& mcParticle : mcParticles) {
      if (!isInAcceptance(mcParticle, HIST("MC/particleSelection"))) {
        continue;
      }

      /// checking the PV z coordinate for the generated collision
      if (applyPvZCutInProcessMcWoColl && applyPvZCutGenColl) {
        const auto mcCollision = mcParticle.mcCollision();
        const float posZ = mcCollision.posZ();
        if (posZ < vertexZMin || posZ > vertexZMax) {
          continue;
        }
      }

      // search for particles from HF decays
      if (keepOnlyHfParticles && !RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, /*searchUpToQuark*/ true)) {
        continue;
      }

      static_for<0, 1>([&](auto pdgSign) {
        fillMCParticleHistograms<pdgSign, o2::track::PID::Electron>(mcParticle, doEl);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Muon>(mcParticle, doMu);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Pion>(mcParticle, doPi);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Kaon>(mcParticle, doKa);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Proton>(mcParticle, doPr);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Deuteron>(mcParticle, doDe);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Triton>(mcParticle, doTr);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Helium3>(mcParticle, doHe);
        fillMCParticleHistograms<pdgSign, o2::track::PID::Alpha>(mcParticle, doAl);
      });
    }

    // Fill TEfficiencies
    static_for<0, 1>([&](auto pdgSign) {
      fillMCEfficiency<pdgSign, o2::track::PID::Electron>(doEl);
      fillMCEfficiency<pdgSign, o2::track::PID::Muon>(doMu);
      fillMCEfficiency<pdgSign, o2::track::PID::Pion>(doPi);
      fillMCEfficiency<pdgSign, o2::track::PID::Kaon>(doKa);
      fillMCEfficiency<pdgSign, o2::track::PID::Proton>(doPr);
      fillMCEfficiency<pdgSign, o2::track::PID::Deuteron>(doDe);
      fillMCEfficiency<pdgSign, o2::track::PID::Triton>(doTr);
      fillMCEfficiency<pdgSign, o2::track::PID::Helium3>(doHe);
      fillMCEfficiency<pdgSign, o2::track::PID::Alpha>(doAl);
    });
  }
  PROCESS_SWITCH(QaEfficiency, processMCWithoutCollisions, "process MC without the collision association", false);

  void processData(CollisionCandidates::iterator const& collision,
                   TrackCandidates const& tracks)
  {

    if (!isCollisionSelected<false>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      if (!isTrackSelected<false>(track, HIST("Data/trackSelection"))) {
        continue;
      }

      histos.fill(HIST("Data/trackLength"), track.length());

      float trackPt = track.pt();
      float trackEta = track.eta();
      float trackPhi = track.phi();
      float trackSign = track.sign();
      float occupancy{};
      float centrality{};
      if (doOccupancy) {
        centrality = collision.centFT0C();
        if (useFT0OccEstimator) {
          /// occupancy estimator (FT0c signal amplitudes in +-10us from current collision)
          occupancy = collision.ft0cOccupancyInTimeRange();
        } else {
          /// occupancy estimator (ITS tracks with at least 5 clusters in +-10us from current collision)
          occupancy = static_cast<float>(collision.trackOccupancyInTimeRange());
        }
      }
      if (passedITS) {
        if (trackSign > 0) {
          histos.fill(HIST("Data/pos/pt/its"), trackPt);
          histos.fill(HIST("Data/pos/eta/its"), trackEta);
          histos.fill(HIST("Data/pos/phi/its"), trackPhi);
          histos.fill(HIST("Data/pos/etaphi/its"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/pos/its"), occupancy, centrality, trackPt);
          }
        } else {
          histos.fill(HIST("Data/neg/pt/its"), trackPt);
          histos.fill(HIST("Data/neg/eta/its"), trackEta);
          histos.fill(HIST("Data/neg/phi/its"), trackPhi);
          histos.fill(HIST("Data/neg/etaphi/its"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/neg/its"), occupancy, centrality, trackPt);
          }
        }
      }

      if (passedTPC) {
        if (trackSign > 0) {
          histos.fill(HIST("Data/pos/pt/tpc"), trackPt);
          histos.fill(HIST("Data/pos/eta/tpc"), trackEta);
          histos.fill(HIST("Data/pos/phi/tpc"), trackPhi);
          histos.fill(HIST("Data/pos/etaphi/tpc"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/pos/tpc"), occupancy, centrality, trackPt);
          }
        } else {
          histos.fill(HIST("Data/neg/pt/tpc"), trackPt);
          histos.fill(HIST("Data/neg/eta/tpc"), trackEta);
          histos.fill(HIST("Data/neg/phi/tpc"), trackPhi);
          histos.fill(HIST("Data/neg/etaphi/tpc"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/neg/tpc"), occupancy, centrality, trackPt);
          }
        }
      }

      if (passedITS && passedTPC) {
        if (trackSign > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc"), trackPt);
          histos.fill(HIST("Data/pos/eta/its_tpc"), trackEta);
          histos.fill(HIST("Data/pos/phi/its_tpc"), trackPhi);
          histos.fill(HIST("Data/pos/etaphi/its_tpc"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/pos/its_tpc"), occupancy, centrality, trackPt);
          }
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc"), trackPt);
          histos.fill(HIST("Data/neg/eta/its_tpc"), trackEta);
          histos.fill(HIST("Data/neg/phi/its_tpc"), trackPhi);
          histos.fill(HIST("Data/neg/etaphi/its_tpc"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/neg/its_tpc"), occupancy, centrality, trackPt);
          }
        }
      }

      if (passedITS && passedTPC && passedTOF) {
        if (trackSign > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_tof"), trackPt);
          histos.fill(HIST("Data/pos/eta/its_tpc_tof"), trackEta);
          histos.fill(HIST("Data/pos/phi/its_tpc_tof"), trackPhi);
          histos.fill(HIST("Data/pos/etaphi/its_tpc_tof"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/pos/its_tpc_tof"), occupancy, centrality, trackPt);
          }
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_tof"), trackPt);
          histos.fill(HIST("Data/neg/eta/its_tpc_tof"), trackEta);
          histos.fill(HIST("Data/neg/phi/its_tpc_tof"), trackPhi);
          histos.fill(HIST("Data/neg/etaphi/its_tpc_tof"), trackEta, trackPhi);

          if (doOccupancy) {
            histos.fill(HIST("Data/occ_cent/neg/its"), occupancy, centrality, trackPt);
          }
        }
      }

      if (makeEff) {
        if (passedITS) {
          effITSTPCMatchingVsPt->Fill(passedTPC, trackPt);
        }
        if (passedTPC) {
          effTPCITSMatchingVsPt->Fill(passedITS, trackPt);
        }
        if (passedITS && passedTPC) {
          effTPCTOFMatchingVsPt->Fill(passedTOF, trackPt);
          effTPCTOFMatchingVsP->Fill(passedTOF, track.p());
          effTPCTOFMatchingVsEta->Fill(passedTOF, trackEta);
          effTPCTOFMatchingVsPhi->Fill(passedTOF, trackPhi);
          effTPCTOFMatchingVsPtVsEta->Fill(passedTOF, trackPt, trackEta);
          effTPCTOFMatchingVsPtVsPhi->Fill(passedTOF, trackPt, trackPhi);
        }
      }
    }
  }
  PROCESS_SWITCH(QaEfficiency, processData, "process data", true);

  void processDataWithPID(CollisionCandidates::iterator const& collision,
                          o2::soa::Join<TrackCandidates, o2::aod::pidTPCLfFullDe> const& tracks)
  {

    if (!isCollisionSelected<false>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      if (!isTrackSelected<false>(track, HIST("Data/trackSelection"))) {
        continue;
      }
      if (std::abs(track.tpcNSigmaDe()) > nsigmaTPCDe) {
        continue;
      }
      histos.fill(HIST("Data/trackLength"), track.length());

      if (passedITS) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its"), track.pt());
          histos.fill(HIST("Data/pos/eta/its"), track.eta());
          histos.fill(HIST("Data/pos/phi/its"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/its"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/pt/its"), track.pt());
          histos.fill(HIST("Data/neg/eta/its"), track.eta());
          histos.fill(HIST("Data/neg/phi/its"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/its"), track.eta(), track.phi());
        }
      }

      if (passedTPC) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc"), track.pt());
          histos.fill(HIST("Data/pos/eta/tpc"), track.eta());
          histos.fill(HIST("Data/pos/phi/tpc"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/tpc"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc"), track.pt());
          histos.fill(HIST("Data/neg/eta/tpc"), track.eta());
          histos.fill(HIST("Data/neg/phi/tpc"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/tpc"), track.eta(), track.phi());
        }
      }

      if (passedITS && passedTPC) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc"), track.pt());
          histos.fill(HIST("Data/pos/eta/its_tpc"), track.eta());
          histos.fill(HIST("Data/pos/phi/its_tpc"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/its_tpc"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc"), track.pt());
          histos.fill(HIST("Data/neg/eta/its_tpc"), track.eta());
          histos.fill(HIST("Data/neg/phi/its_tpc"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/its_tpc"), track.eta(), track.phi());
        }
      }

      if (passedITS && passedTPC && passedTOF) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_tof"), track.pt());
          histos.fill(HIST("Data/pos/eta/its_tpc_tof"), track.eta());
          histos.fill(HIST("Data/pos/phi/its_tpc_tof"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/its_tpc_tof"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_tof"), track.pt());
          histos.fill(HIST("Data/neg/eta/its_tpc_tof"), track.eta());
          histos.fill(HIST("Data/neg/phi/its_tpc_tof"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/its_tpc_tof"), track.eta(), track.phi());
        }
      }

      if (makeEff) {
        if (passedITS) {
          effITSTPCMatchingVsPt->Fill(passedTPC, track.pt());
        }
        if (passedTPC) {
          effTPCITSMatchingVsPt->Fill(passedITS, track.pt());
        }
        if (passedITS && passedTPC) {
          effTPCTOFMatchingVsPt->Fill(passedTOF, track.pt());
          effTPCTOFMatchingVsP->Fill(passedTOF, track.p());
          effTPCTOFMatchingVsEta->Fill(passedTOF, track.eta());
          effTPCTOFMatchingVsPhi->Fill(passedTOF, track.phi());
          effTPCTOFMatchingVsPtVsEta->Fill(passedTOF, track.pt(), track.eta());
          effTPCTOFMatchingVsPtVsPhi->Fill(passedTOF, track.pt(), track.phi());
        }
      }
    }
  }
  PROCESS_SWITCH(QaEfficiency, processDataWithPID, "process data with PID", false);

  void processHmpid(CollisionCandidates::iterator const& collision, TrackCandidates const&,
                    o2::aod::HMPIDs const& hmpids)
  {

    if (!isCollisionSelected<false>(collision)) {
      return;
    }

    for (const auto& hmpid : hmpids) {
      const auto& track = hmpid.track_as<TrackCandidates>();
      if (!isTrackSelected<false, false>(track)) {
        continue;
      }

      if (track.sign() > 0) {
        histos.fill(HIST("Data/pos/hmpidMomDiff"), track.p() - hmpid.hmpidMom());
        histos.fill(HIST("Data/pos/pt/hmpid"), track.pt());
        histos.fill(HIST("Data/pos/eta/hmpid"), track.eta());
        histos.fill(HIST("Data/pos/phi/hmpid"), track.phi());
        histos.fill(HIST("Data/pos/etaphi/hmpid"), track.eta(), track.phi());
      } else {
        histos.fill(HIST("Data/neg/hmpidMomDiff"), track.p() - hmpid.hmpidMom());
        histos.fill(HIST("Data/neg/pt/hmpid"), track.pt());
        histos.fill(HIST("Data/neg/eta/hmpid"), track.eta());
        histos.fill(HIST("Data/neg/phi/hmpid"), track.phi());
        histos.fill(HIST("Data/neg/etaphi/hmpid"), track.eta(), track.phi());
      }
    }
  }
  PROCESS_SWITCH(QaEfficiency, processHmpid, "process HMPID matching", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<QaEfficiency>(cfgc)}; }
