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
/// \file   qaDcaMC.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to analyse the DCA distributions in the MC according to PDG code and status
///

// O2 includes
#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

// ROOT includes
#include "TEfficiency.h"
#include "THashList.h"
#include "TPDGCode.h"

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
static constexpr int trkCutIdxN = 22;

// Particle information
static constexpr int nSpecies = o2::track::PID::NIDs; // One per PDG
static constexpr int nCharges = 2;                    // Positive and negative
static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
static constexpr const char* particleNames[nSpecies] = {"el", "mu", "pi", "ka", "pr", "de", "tr", "he", "al"};
static constexpr const char* chargeNames[nSpecies] = {"pos_pdg", "neg_pdg"};
static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040};

// Histograms
// Pt
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtPrm;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtStr;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtMat;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtMatSM;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtMatSM1Dau;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtMotherOfPDG1;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtMotherOfPDG2;
std::array<std::array<std::shared_ptr<TH3>, nCharges>, nSpecies> hPtMotherOfPDG3;

struct QaDcaMc {
  // Track/particle selection
  Configurable<float> maxProdRadius{"maxProdRadius", 9999.f, "Maximum production radius of the particle under study"};
  Configurable<int> motherPDG{"motherPDG", 0, "PDG code of the mother particle"};
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
  Configurable<int> minNClustersITS{"minNClustersITS", -1, "Minimum required number of ITS clusters"};
  Configurable<bool> enableTrackSelection{"enableTrackSelection", true, "Enable the track selection"};

  // Event selection
  Configurable<bool> enableEventSelection{"enableEventSelection", true, "Enable the event selection"};
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the generated vertez in Z (cm)"};
  Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the generated vertez in Z (cm)"};
  // Histogram configuration
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  Configurable<int> logPt{"log-pt", 0, "Flag to use a logarithmic pT axis"};
  ConfigurableAxis etaBins{"etaBins", {200, -3.f, 3.f}, "Eta binning"};
  ConfigurableAxis phiBins{"phiBins", {200, 0.f, 6.284f}, "Phi binning"};
  ConfigurableAxis yBins{"yBins", {200, -0.5f, 0.5f}, "Y binning"};
  ConfigurableAxis dcaBinsxy{"dcaBinsxy", {500, -1.f, 1.f}, "DCAxy binning"};
  ConfigurableAxis dcaBinsz{"dcaBinsz", {100, -0.1f, 0.1f}, "DCAz binning"};
  Configurable<bool> doPVContributorCut{"doPVContributorCut", false, "Select tracks used for primary vertex recostruction (isPVContributor)"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static const char* particleName(int pdgSign, o2::track::PID::ID id)
  {
    return Form("%s %s", pdgSign == 0 ? "Positive PDG" : "Negative PDG", o2::track::PID::getName(id));
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

    AxisSpec axisPt{ptBins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisP{ptBins, "#it{p} (GeV/#it{c})"};
    if (logPt) {
      axisPt.makeLogarithmic();
      axisP.makeLogarithmic();
    }
    const AxisSpec axisEta{etaBins, "#it{#eta}"};
    const AxisSpec axisY{yBins, "#it{y}"};
    const AxisSpec axisPhi{phiBins, "#it{#varphi} (rad)"};
    const AxisSpec axisDCAxy{dcaBinsxy, "DCA_{xy} (cm)"};
    const AxisSpec axisDCAz{dcaBinsz, "DCA_{z} (cm)"};

    const char* partName = particleName(pdgSign, id);
    LOG(info) << "Preparing histograms for particle: " << partName << " pdgSign " << pdgSign;

    const TString tagPt = Form("%s #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                               partName,
                               etaMin, etaMax,
                               yMin, yMax,
                               phiMin, phiMax);

    hPtPrm[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/prm", particleNames[id], chargeNames[pdgSign]), "DCA Prm. " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
    hPtStr[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/str", particleNames[id], chargeNames[pdgSign]), "DCA Str. " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
    hPtMat[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/mat", particleNames[id], chargeNames[pdgSign]), "DCA Mat. " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
    hPtMatSM[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/sm", particleNames[id], chargeNames[pdgSign]), "DCA Mat. SM " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
    hPtMatSM1Dau[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/sm1dau", particleNames[id], chargeNames[pdgSign]), "DCA Mat. SM 1 Dau " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
    if (motherPDG.value != 0) {
      hPtMotherOfPDG1[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/motherpdg1", particleNames[id], chargeNames[pdgSign]), "DCA mother pdg " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
      hPtMotherOfPDG2[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/motherpdg2", particleNames[id], chargeNames[pdgSign]), "DCA mother pdg " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
      hPtMotherOfPDG3[id][pdgSign] = histos.add<TH3>(Form("MC/%s/%s/dcaxyz/pt/motherpdg3", particleNames[id], chargeNames[pdgSign]), "DCA mother pdg " + tagPt, kTH3F, {axisPt, axisDCAxy, axisDCAz});
    }

    LOG(info) << "Done with making histograms for particle: " << partName;
  }

  void initMC(const AxisSpec& axisSel)
  {
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
    h->GetXaxis()->SetBinLabel(trkCutIdxPassedGlobal, "No extra selection");

    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(trkCutIdxN + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1D, {{1, 0, 1}});

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
    });

    histos.print();
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
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel. (sel8)");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Position");

    const AxisSpec axisSel{40, 0.5, 40.5, "Selection"};
    initMC(axisSel);
  }

  template <int pdgSign, o2::track::PID::ID id>
  bool isPdgSelected(const o2::aod::McParticles::iterator& mcParticle)
  {
    static_assert(pdgSign == 0 || pdgSign == 1);
    static_assert(id > 0 || id < nSpecies);

    // Selecting a specific PDG
    if constexpr (pdgSign == 0) {
      return mcParticle.pdgCode() == PDGs[id];
    } else {
      return mcParticle.pdgCode() == -PDGs[id];
    }
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

    const o2::aod::McParticles::iterator& mcParticle = track.mcParticle();

    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("MC/trackSelection"), trkCutIdxN + id);

    if (motherPDG.value != 0) {
      int motherPdgCounter = 0;
      if (mcParticle.has_mothers()) {
        const auto& mothers = mcParticle.mothers_as<o2::aod::McParticles>();
        for (const auto& mother : mothers) {
          if (!mother.has_mothers()) {
            continue;
          }
          const auto& mothers2 = mother.mothers_as<o2::aod::McParticles>();
          for (const auto& mother2 : mothers2) {
            if (abs(mother2.pdgCode()) == motherPDG.value) {
              motherPdgCounter++;
            }
          }
        }
      }
      if (motherPdgCounter == 1) {
        hPtMotherOfPDG1[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
      } else if (motherPdgCounter == 2) {
        hPtMotherOfPDG2[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
      } else if (motherPdgCounter > 0) {
        hPtMotherOfPDG3[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
      }
    }

    if (isPhysicalPrimary(mcParticle)) {
      hPtPrm[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
    } else if (mcParticle.getProcess() == 4) { // Particle decay
      hPtStr[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
    } else { // Material
      hPtMat[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
      if (mcParticle.has_mothers()) {
        if (mcParticle.mothers_as<o2::aod::McParticles>()[0].pdgCode() == mcParticle.pdgCode()) {
          hPtMatSM[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
          if (mcParticle.mothers_as<o2::aod::McParticles>().size() == 1) {
            hPtMatSM1Dau[id][pdgSign]->Fill(mcParticle.pt(), track.dcaXY(), track.dcaZ());
          }
        }
      }
    }
  }

  template <bool doFillHistograms, typename CollType>
  bool isCollisionSelected(const CollType& collision)
  {
    if constexpr (doFillHistograms) {
      histos.fill(HIST("eventSelection"), 1);
    }
    if (!enableEventSelection.value) {
      return true;
    }
    if (!collision.sel8()) {
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
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension, o2::aod::TracksDCA>;
  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>::iterator const& collision,
               o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels> const& tracks,
               o2::aod::McParticles const&)
  {
    if (!isCollisionSelected<true>(collision)) {
      return;
    }

    // Track loop
    for (const auto& track : tracks) {
      if (!isTrackSelected(track, HIST("MC/trackSelection"))) {
        continue;
      }
      // Filling variable histograms
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
      if (!enableTrackSelection.value) {
        return true;
      }

      if constexpr (doFillHisto) {
        histos.fill(countingHisto, trkCutIdxHasMcPart); // Tracks with particles (i.e. no fakes)
      }
      const auto mcParticle = track.mcParticle();
      if (!isInAcceptance<true, doFillHisto>(mcParticle, countingHisto, trkCutIdxHasMcPart)) {
        // 3: pt cut 4: eta cut 5: phi cut 6: y cut
        return false;
      }

      // if (noFakesHits) {              // Selecting tracks with no fake hits
      //   for (int i = 0; i < 7; i++) { // ITS
      //     if (track.mcMismatchInITS(i)) {
      //       return false;
      //     }
      //   }
      //   for (int i = 7; i < 10; i++) { // TPC
      //     if (track.mcMismatchInTPC(i)) {
      //       return false;
      //     }
      //   }
      // }
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
    // if (!track.passedDCAxy()) {
    //   return false;
    // }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, trkCutIdxPassedDcaXYMax);
    }
    // if (std::abs(track.dcaXY()) < minDcaXY) {
    //   return false;
    // }
    if (track.itsNCls() < minNClustersITS) {
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
    // if (std::abs(track.dcaZ()) < minDcaZ) {
    //   return false;
    // }
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
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if constexpr (doFillHisto) {
      histos.fill(countingHisto, trkCutIdxPassedGlobal);
    }

    return true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<QaDcaMc>(cfgc)}; }
