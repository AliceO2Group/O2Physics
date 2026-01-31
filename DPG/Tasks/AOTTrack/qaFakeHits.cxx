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
/// \file   qaFakeHits.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2024-04-08
/// \brief  Task to analyze the fraction of the true and fake hits depending on where the fake hits are picked
///

#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>

using namespace o2::framework;
// Particle information
static constexpr int nSpecies = o2::track::PID::NIDs; // One per PDG
static constexpr int nCharges = 2;
static constexpr int nParticles = nSpecies * nCharges;
static constexpr const char* particleTitle[nParticles] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha",
                                                          "e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
static constexpr int PDGs[nParticles] = {11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040,
                                         -11, -13, -211, -321, -2212, -1000010020, -1000010030, -1000020030, -1000020040};
std::array<std::shared_ptr<TH2>, nParticles> hPtAll;
std::array<std::shared_ptr<TH2>, nParticles> hPtITS;
std::array<std::shared_ptr<TH2>, nParticles> hPtITSAB;
std::array<std::shared_ptr<TH2>, nParticles> hPtITSTPC;
std::array<std::shared_ptr<TH2>, nParticles> hPtTPC;
std::array<std::shared_ptr<TH2>, nParticles> hPtTRD;
std::array<std::shared_ptr<TH2>, nParticles> hPtTOF;
std::array<std::shared_ptr<TH2>, nParticles> hPtOverall;
std::array<std::shared_ptr<TH2>, nParticles> hPtNoise;
std::array<std::shared_ptr<TH2>, nParticles> hPtNoMismatch;

struct QaFakeHits {
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
  // Track only selection, options to select only specific tracks
  Configurable<bool> trackSelection{"trackSelection", true, "Local track selection"};
  Configurable<int> globalTrackSelection{"globalTrackSelection", 0, "Global track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks, 6 -> custom track cuts via Configurable"};
  // Event selection
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the generated vertez in Z (cm)"};
  Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the generated vertez in Z (cm)"};
  // Histogram configuration
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  // Histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static const char* particleName(int pdgSign, o2::track::PID::ID id)
  {
    if (pdgSign == 0) { // Positive PDG
      return particleTitle[id];
    }
    // Negative PDG
    return particleTitle[id + o2::track::PID::NIDs];
  }

  void makeMCHistograms(const bool doMakeHistograms,
                        const int pdgSign,
                        o2::track::PID::ID id)
  {
    if (!doMakeHistograms) {
      return;
    }

    switch (pdgSign) {
      case 0: // Positive
        if (!doPositivePDG) {
          return;
        }
        break;
      case 1: // Negative
        if (!doNegativePDG) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Can't interpret pdgSign " << pdgSign;
    }

    const AxisSpec axisPt{ptBins, "#it{p}_{T}^{gen} (GeV/#it{c})"};
    const AxisSpec axisPtDiff{100, -1, 1, "#it{p}_{T}^{gen} - #it{p}_{T}^{rec} (GeV/#it{c})"};

    const char* partName = particleName(pdgSign, id);
    LOG(info) << "Preparing histograms for particle: " << partName << " pdgSign " << pdgSign;
    const int histogramIndex = id + pdgSign * nSpecies;

    const TString tagPt = Form("%s ", partName);
    hPtAll[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/all", PDGs[histogramIndex]), "All tracks " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtITS[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/its", PDGs[histogramIndex]), "ITS mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtITSAB[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/itsab", PDGs[histogramIndex]), "ITSAB mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtTPC[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/tpc", PDGs[histogramIndex]), "TPC mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtITSTPC[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/itstpc", PDGs[histogramIndex]), "ITSTPC mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtTRD[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/trd", PDGs[histogramIndex]), "TRD mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtTOF[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/tof", PDGs[histogramIndex]), "TOF mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtOverall[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/overall", PDGs[histogramIndex]), "Overall mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtNoise[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/noise", PDGs[histogramIndex]), "Noise " + tagPt, kTH2D, {axisPt, axisPtDiff});
    hPtNoMismatch[histogramIndex] = histos.add<TH2>(Form("MC/pdg%i/pt/mismatched/no", PDGs[histogramIndex]), "No mismatch " + tagPt, kTH2D, {axisPt, axisPtDiff});

    LOG(info) << "Done with particle: " << partName;
  }

  void init(InitContext&)
  {
    for (int pdgSign = 0; pdgSign < 2; pdgSign++) {
      makeMCHistograms(doEl, pdgSign, o2::track::PID::Electron);
      makeMCHistograms(doMu, pdgSign, o2::track::PID::Muon);
      makeMCHistograms(doPi, pdgSign, o2::track::PID::Pion);
      makeMCHistograms(doKa, pdgSign, o2::track::PID::Kaon);
      makeMCHistograms(doPr, pdgSign, o2::track::PID::Proton);
      makeMCHistograms(doDe, pdgSign, o2::track::PID::Deuteron);
      makeMCHistograms(doTr, pdgSign, o2::track::PID::Triton);
      makeMCHistograms(doHe, pdgSign, o2::track::PID::Helium3);
      makeMCHistograms(doAl, pdgSign, o2::track::PID::Alpha);
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

  template <typename TrackType>
  bool isMismatched(TrackType const& track, int layer)
  {
    return (track.mcMask() & 1 << layer);
  }

  template <int pdgSign, o2::track::PID::ID id, typename trackType>
  void fillMCTrackHistograms(const trackType& track, const bool doMakeHistograms)
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
    if (!track.has_mcParticle()) {
      return;
    }
    const auto& mcParticle = track.mcParticle();

    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }
    if (!track.isGlobalTrack()) {
      return;
    }
    const float ptDiff = mcParticle.pt() - track.pt();
    hPtAll[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    bool mismatchInITS = false;
    for (int i = 0; i < 7; i++) {
      if (isMismatched(track, i)) {
        mismatchInITS = true;
        break;
      }
    }
    if (mismatchInITS) {
      hPtITS[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }

    bool mismatchInTPC = false;
    for (int i = 7; i < 10; i++) {
      if (isMismatched(track, i)) {
        mismatchInTPC = true;
        break;
      }
    }
    if (mismatchInTPC) {
      hPtTPC[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }

    if (track.mcMask() == 0) {
      hPtNoMismatch[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
    const bool mismatchInTRD = isMismatched(track, 10); // Not set
    if (mismatchInTRD) {
      hPtTRD[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
    const bool mismatchInTOF = isMismatched(track, 15) && track.hasTOF();
    if (mismatchInTOF) {
      hPtTOF[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
    const bool itsAfterBurnerMismatch = isMismatched(track, 12);
    if (itsAfterBurnerMismatch) {
      hPtITSAB[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
    const bool itstpcMismatch = isMismatched(track, 13);
    if (itstpcMismatch) {
      hPtITSTPC[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
    const bool isNoise = isMismatched(track, 14);
    if (isNoise) {
      hPtNoise[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
    const bool overallMismatch = isMismatched(track, 15);
    if (overallMismatch) {
      hPtOverall[histogramIndex]->Fill(mcParticle.pt(), ptDiff);
    }
  }

  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension, o2::aod::TracksDCA>;
  using TrackCandidatesMC = o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels>;

  // MC process
  void process(o2::aod::Collision const& /*collision*/,
               TrackCandidatesMC const& tracks,
               o2::aod::McCollisions const&,
               o2::aod::McParticles const&)
  {
    // Track loop
    for (const auto& track : tracks) {

      LOG(debug) << "Processing track with global index " << track.globalIndex() << " and mask " << track.mcMask();
      LOG(debug) << "mcMask bits: " << std::bitset<16>(track.mcMask());

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
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QaFakeHits>(cfgc)};
}
