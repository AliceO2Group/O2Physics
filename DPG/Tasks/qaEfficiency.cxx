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
/// \author Peter Hristov <Peter.Hristov@cern.ch>, CERN
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Henrique J C Zanoli <henrique.zanoli@cern.ch>, Utrecht University
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

// ROOT includes
#include "TPDGCode.h"
#include "TEfficiency.h"
#include "TList.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"eff-data", VariantType::Int, 1, {"Efficiency for the data"}},
    {"eff-el", VariantType::Int, 0, {"Efficiency for the Electron PDG code"}},
    {"eff-mu", VariantType::Int, 0, {"Efficiency for the Muon PDG code"}},
    {"eff-pi", VariantType::Int, 0, {"Efficiency for the Pion PDG code"}},
    {"eff-ka", VariantType::Int, 0, {"Efficiency for the Kaon PDG code"}},
    {"eff-pr", VariantType::Int, 0, {"Efficiency for the Proton PDG code"}},
    {"eff-de", VariantType::Int, 0, {"Efficiency for the Deuteron PDG code"}},
    {"eff-tr", VariantType::Int, 0, {"Efficiency for the Triton PDG code"}},
    {"eff-he", VariantType::Int, 0, {"Efficiency for the Helium3 PDG code"}},
    {"eff-el-pos", VariantType::Int, 0, {"Efficiency for the Electron Positive PDG code"}},
    {"eff-mu-pos", VariantType::Int, 0, {"Efficiency for the Muon Positive PDG code"}},
    {"eff-pi-pos", VariantType::Int, 0, {"Efficiency for the Pion Positive PDG code"}},
    {"eff-ka-pos", VariantType::Int, 0, {"Efficiency for the Kaon Positive PDG code"}},
    {"eff-pr-pos", VariantType::Int, 0, {"Efficiency for the Proton Positive PDG code"}},
    {"eff-de-pos", VariantType::Int, 0, {"Efficiency for the Deuteron Positive PDG code"}},
    {"eff-tr-pos", VariantType::Int, 0, {"Efficiency for the Triton Positive PDG code"}},
    {"eff-he-pos", VariantType::Int, 0, {"Efficiency for the Helium3 Positive PDG code"}},
    {"eff-el-neg", VariantType::Int, 0, {"Efficiency for the Electron Negative PDG code"}},
    {"eff-mu-neg", VariantType::Int, 0, {"Efficiency for the Muon Negative PDG code"}},
    {"eff-pi-neg", VariantType::Int, 0, {"Efficiency for the Pion Negative PDG code"}},
    {"eff-ka-neg", VariantType::Int, 0, {"Efficiency for the Kaon Negative PDG code"}},
    {"eff-pr-neg", VariantType::Int, 0, {"Efficiency for the Proton Negative PDG code"}},
    {"eff-de-neg", VariantType::Int, 0, {"Efficiency for the Deuteron Negative PDG code"}},
    {"eff-tr-neg", VariantType::Int, 0, {"Efficiency for the Triton Negative PDG code"}},
    {"eff-he-neg", VariantType::Int, 0, {"Efficiency for the Helium3 Negative PDG code"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to QA the efficiency of a particular particle defined by its pdg code
template <o2::track::pid_constants::ID particle, int pdgSign>
struct QaTrackingEfficiency {
  static constexpr int nSpecies = 8;
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030};
  static_assert(particle < nSpecies && "Maximum of particles reached");
  static constexpr int pdg = PDGs[particle];
  // Particle selection
  Configurable<float> etaMin{"eta-min", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 3.f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 6.284f, "Upper limit in phi"};
  Configurable<float> yMin{"y-min", -0.5f, "Lower limit in y"};
  Configurable<float> yMax{"y-max", 0.5f, "Upper limit in y"};
  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 5.f, "Upper limit in pT"};
  Configurable<int> selPrim{"sel-prim", 1, "1 select primaries, 0 select all particles"};
  Configurable<bool> noFakes{"noFakes", false, "Flag to reject tracks that have fake hits"};
  // Event selection
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the generated vertez in Z (cm)"};
  Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the generated vertez in Z (cm)"};
  // Histogram configuration
  Configurable<int> ptBins{"pt-bins", 500, "Number of pT bins"};
  Configurable<int> logPt{"log-pt", 0, "Flag to use a logarithmic pT axis"};
  Configurable<int> etaBins{"eta-bins", 500, "Number of eta bins"};
  Configurable<int> yBins{"y-bins", 500, "Number of eta bins"};
  Configurable<int> phiBins{"phi-bins", 500, "Number of phi bins"};
  // Task configuration
  Configurable<bool> makeEff{"make-eff", false, "Flag to produce the efficiency with TEfficiency"};

  OutputObj<TList> list{"Efficiency"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    if constexpr (pdgSign != 0 && pdgSign != 1 && pdgSign != -1) {
      LOG(fatal) << "Provide pdgSign as 0, 1, -1. Provided: " << pdgSign;
    }
    AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} (GeV/#it{c})"};
    if (logPt) {
      axisPt.makeLogaritmic();
    }
    AxisSpec axisP{ptBins, ptMin, ptMax, "#it{p} (GeV/#it{c})"};
    if (logPt) {
      axisP.makeLogaritmic();
    }
    const AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};
    const AxisSpec axisY{yBins, yMin, yMax, "#it{y}"};
    const AxisSpec axisPhi{phiBins, phiMin, phiMax, "#it{#varphi} (rad)"};

    const AxisSpec axisSel{20, 0.5, 20.5, "Selection"};
    histos.add("eventSelection", "Event Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Position");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Ev. Sel.");

    histos.add("trackSelection", "Track Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(1, "Tracks read");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Reco.");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(6, "Passed y");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(7, "Passed Prim.");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(8, Form("Passed PDG %i", pdg));
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(9, "Passed Prim. MC");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(10, "Passed Fake");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(11, "Passed standard quality cuts");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(12, "Passed has collision");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(13, "Passed hasTOF");
    histos.add("fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1D, {{1, 0, 1}});

    histos.add("partSelection", "Particle Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(1, "Particles read");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Reco.");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(6, "Passed y");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(7, "Passed Prim.");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(8, Form("Passed PDG %i", pdg));
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(9, "Passed Prim. MC");

    histos.add("eventMultiplicity", "Event Selection", kTH1D, {{1000, 0, 5000}});
    histos.add("trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    const TString tagPt = Form("%s #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f] Prim %i",
                               o2::track::pid_constants::sNames[particle],
                               etaMin.value, etaMax.value,
                               yMin.value, yMax.value,
                               phiMin.value, phiMax.value,
                               selPrim.value);
    histos.add("pt/num", "Numerator " + tagPt, kTH1D, {axisPt});
    histos.add("pt/numtof", "Numerator TOF " + tagPt, kTH1D, {axisPt});
    histos.add("pt/den", "Denominator " + tagPt, kTH1D, {axisPt});

    histos.add("p/num", "Numerator " + tagPt, kTH1D, {axisP});
    histos.add("p/numtof", "Numerator TOF " + tagPt, kTH1D, {axisP});
    histos.add("p/den", "Denominator " + tagPt, kTH1D, {axisP});

    const TString tagEta = Form("%s #it{p}_{T} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f] Prim %i",
                                o2::track::pid_constants::sNames[particle],
                                ptMin.value, ptMax.value,
                                yMin.value, yMax.value,
                                phiMin.value, phiMax.value,
                                selPrim.value);
    histos.add("eta/num", "Numerator " + tagEta, kTH1D, {axisEta});
    histos.add("eta/numtof", "Numerator TOF " + tagEta, kTH1D, {axisEta});
    histos.add("eta/den", "Denominator " + tagEta, kTH1D, {axisEta});

    const TString tagY = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f] Prim %i",
                              o2::track::pid_constants::sNames[particle],
                              ptMin.value, ptMax.value,
                              etaMin.value, etaMax.value,
                              phiMin.value, phiMax.value,
                              selPrim.value);
    histos.add("y/num", "Numerator " + tagY, kTH1D, {axisY});
    histos.add("y/numtof", "Numerator TOF " + tagY, kTH1D, {axisY});
    histos.add("y/den", "Denominator " + tagY, kTH1D, {axisY});

    const TString tagPhi = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] Prim %i",
                                o2::track::pid_constants::sNames[particle],
                                ptMin.value, ptMax.value,
                                etaMin.value, etaMax.value,
                                yMin.value, yMax.value,
                                selPrim.value);
    histos.add("phi/num", "Numerator " + tagPhi, kTH1D, {axisPhi});
    histos.add("phi/numtof", "Numerator TOF " + tagPhi, kTH1D, {axisPhi});
    histos.add("phi/den", "Denominator " + tagPhi, kTH1D, {axisPhi});

    const TString tagPtEta = Form("%s #it{#varphi} [%.2f,%.2f] #it{y} [%.2f,%.2f] Prim %i",
                                  o2::track::pid_constants::sNames[particle],
                                  phiMin.value, phiMax.value,
                                  yMin.value, yMax.value,
                                  selPrim.value);
    histos.add("pteta/num", "Numerator " + tagPtEta, kTH2D, {axisPt, axisEta});
    histos.add("pteta/den", "Denominator " + tagPtEta, kTH2D, {axisPt, axisEta});

    list.setObject(new TList);
    if (makeEff) {
      auto makeEfficiency = [&](TString effname, TString efftitle, auto templateHisto) {
        const TAxis* axis = histos.get<TH1>(templateHisto)->GetXaxis();
        if (axis->IsVariableBinSize()) {
          list->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
        } else {
          list->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
        }
      };
      makeEfficiency("efficiencyVsPt", "Efficiency " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("pt/num"));
      makeEfficiency("efficiencyVsP", "Efficiency " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST("pt/num"));
      makeEfficiency("efficiencyVsEta", "Efficiency " + tagEta + ";#it{#eta};Efficiency", HIST("eta/num"));
      makeEfficiency("efficiencyVsPhi", "Efficiency " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST("phi/num"));

      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHisto) {
        const TAxis* axisX = histos.get<TH2>(templateHisto)->GetXaxis();
        const TAxis* axisY = histos.get<TH2>(templateHisto)->GetXaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          list->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
        } else {
          list->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
        }
      };
      makeEfficiency2D("efficiencyVsPtVsEta", "Efficiency " + tagPtEta + ";#it{#varphi} (rad);Efficiency", HIST("pteta/num"));
    }
  }

  void process(const o2::aod::McParticles& mcParticles,
               const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>& collisions,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::McTrackLabels, o2::aod::TrackSelection>& tracks,
               const o2::aod::McCollisions&)
  {

    std::vector<int64_t> recoEvt(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      histos.fill(HIST("eventSelection"), 1);
      if (collision.numContrib() < nMinNumberOfContributors) {
        continue;
      }
      histos.fill(HIST("eventSelection"), 2);
      const auto mcCollision = collision.mcCollision();
      if ((mcCollision.posZ() < vertexZMin || mcCollision.posZ() > vertexZMax)) {
        continue;
      }
      histos.fill(HIST("eventSelection"), 3);
      if (!collision.sel8()) {
        continue;
      }
      histos.fill(HIST("eventSelection"), 4);
      recoEvt[nevts++] = mcCollision.globalIndex();
    }
    recoEvt.resize(nevts);

    auto rejectParticle = [&](const auto& p, auto h) {
      histos.fill(h, 1);
      const auto evtReconstructed = std::find(recoEvt.begin(), recoEvt.end(), p.mcCollision().globalIndex()) != recoEvt.end();
      if (!evtReconstructed) { // Check that the event is reconstructed
        return true;
      }

      histos.fill(h, 2);
      if ((p.pt() < ptMin || p.pt() > ptMax)) { // Check pt
        return true;
      }
      histos.fill(h, 3);
      if ((p.eta() < etaMin || p.eta() > etaMax)) { // Check eta
        return true;
      }
      histos.fill(h, 4);
      if ((p.phi() < phiMin || p.phi() > phiMax)) { // Check phi
        return true;
      }
      histos.fill(h, 5);
      if ((p.y() < yMin || p.y() > yMax)) { // Check rapidity
        return true;
      }
      histos.fill(h, 6);
      if ((selPrim == 1) && !p.isPhysicalPrimary()) { // Requiring is physical primary
        return true;
      }
      histos.fill(h, 7);

      // Selecting PDG code
      switch (pdgSign) {
        case 0:
          if (abs(p.pdgCode()) != pdg) {
            return true;
          }
          break;
        case 1:
          if (p.pdgCode() != pdg) {
            return true;
          }
          break;
        case -1:
          if (p.pdgCode() != -pdg) {
            return true;
          }
          break;
        default:
          LOG(fatal) << "Provide pdgSign as 0, 1, -1. Provided: " << pdgSign;
          break;
      }
      histos.fill(h, 8);
      // Select primaries based on position
      const float dx = p.vx() - p.mcCollision().posX();
      const float dy = p.vy() - p.mcCollision().posY();
      const float dz = p.vz() - p.mcCollision().posZ();
      if (sqrt(dx * dx + dy * dy + dz * dz) > 0.0001) {
        return true;
      }
      histos.fill(h, 9);

      return false;
    };

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        histos.fill(HIST("fakeTrackNoiseHits"), 0.5);
        continue;
      }
      const auto mcParticle = track.mcParticle();
      if (rejectParticle(mcParticle, HIST("trackSelection"))) {
        continue;
      }

      if (noFakes) { // Selecting tracks with no fake hits
        bool hasFake = false;
        for (int i = 0; i < 10; i++) { // From ITS to TPC
          if (track.mcMask() & 1 << i) {
            hasFake = true;
            break;
          }
        }
        if (hasFake) {
          continue;
        }
      }

      histos.fill(HIST("trackSelection"), 10);
      if (!track.isGlobalTrack()) { // Check general cuts
        continue;
      }
      histos.fill(HIST("trackSelection"), 11);
      if (!track.has_collision()) {
        continue;
      }
      histos.fill(HIST("trackSelection"), 12);
      histos.fill(HIST("trackLength"), track.length());
      histos.fill(HIST("p/num"), mcParticle.p());
      histos.fill(HIST("pt/num"), mcParticle.pt());
      histos.fill(HIST("eta/num"), mcParticle.eta());
      histos.fill(HIST("y/num"), mcParticle.y());
      histos.fill(HIST("phi/num"), mcParticle.phi());
      histos.fill(HIST("pteta/num"), mcParticle.pt(), mcParticle.eta());
      if (!track.hasTOF()) {
        continue;
      }
      histos.fill(HIST("trackSelection"), 13);
      histos.fill(HIST("p/numtof"), mcParticle.p());
      histos.fill(HIST("pt/numtof"), mcParticle.pt());
      histos.fill(HIST("eta/numtof"), mcParticle.eta());
      histos.fill(HIST("y/numtof"), mcParticle.y());
      histos.fill(HIST("phi/numtof"), mcParticle.phi());
    }

    float dNdEta = 0;
    for (const auto& mcParticle : mcParticles) {
      if (TMath::Abs(mcParticle.eta()) <= 2.f && !mcParticle.has_daughters()) {
        dNdEta += 1.f;
      }
      if (rejectParticle(mcParticle, HIST("partSelection"))) {
        continue;
      }

      histos.fill(HIST("p/den"), mcParticle.p());
      histos.fill(HIST("pt/den"), mcParticle.pt());
      histos.fill(HIST("eta/den"), mcParticle.eta());
      histos.fill(HIST("y/den"), mcParticle.y());
      histos.fill(HIST("phi/den"), mcParticle.phi());
      histos.fill(HIST("pteta/den"), mcParticle.pt(), mcParticle.eta());
    }
    histos.fill(HIST("eventMultiplicity"), dNdEta * 0.5f / 2.f);

    if (makeEff) {
      auto fillEfficiency = [&](int index, auto num, auto den) {
        static_cast<TEfficiency*>(list->At(index))->SetTotalHistogram(*histos.get<TH1>(den).get(), "f");
        static_cast<TEfficiency*>(list->At(index))->SetPassedHistogram(*histos.get<TH1>(num).get(), "f");
      };
      fillEfficiency(0, HIST("pt/num"), HIST("pt/den"));
      fillEfficiency(1, HIST("p/num"), HIST("p/den"));
      fillEfficiency(2, HIST("eta/num"), HIST("eta/den"));
      fillEfficiency(3, HIST("phi/num"), HIST("phi/den"));

      auto fillEfficiency2D = [&](int index, auto num, auto den) {
        static_cast<TEfficiency*>(list->At(index))->SetTotalHistogram(*histos.get<TH2>(den).get(), "f");
        static_cast<TEfficiency*>(list->At(index))->SetPassedHistogram(*histos.get<TH2>(num).get(), "f");
      };
      fillEfficiency2D(4, HIST("pteta/num"), HIST("pteta/den"));
    }
  }
};

/// Task to QA the efficiency of a particular particle defined by its pdg code
struct QaTrackingEfficiencyData {
  // Track selection
  Configurable<float> etaMin{"eta-min", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 3.f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 6.284f, "Upper limit in phi"};
  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 5.f, "Upper limit in pT"};
  Configurable<bool> trackSelection{"track-sel", false, "Flag to use the standard track selection when selecting numerator and denominator"};
  // Event selection
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the generated vertez in Z (cm)"};
  Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the generated vertez in Z (cm)"};
  // Histogram configuration
  Configurable<int> ptBins{"pt-bins", 500, "Number of pT bins"};
  Configurable<int> logPt{"log-pt", 0, "Flag to use a logarithmic pT axis"};
  Configurable<int> etaBins{"eta-bins", 500, "Number of eta bins"};
  Configurable<int> phiBins{"phi-bins", 500, "Number of phi bins"};
  // Task configuration
  Configurable<bool> makeEff{"make-eff", false, "Flag to produce the efficiency with TEfficiency"};

  OutputObj<TList> list{"Efficiency"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const TString tagPt = Form("#it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                               etaMin.value, etaMax.value,
                               phiMin.value, phiMax.value);
    AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} (GeV/#it{c})"};
    if (logPt) {
      axisPt.makeLogaritmic();
    }

    const TString tagEta = Form("#it{p}_{T} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                                ptMin.value, ptMax.value,
                                phiMin.value, phiMax.value);
    const AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};

    const TString tagPhi = Form("#it{#eta} [%.2f,%.2f] #it{p}_{T} [%.2f,%.2f]",
                                etaMin.value, etaMax.value,
                                ptMin.value, ptMax.value);
    const AxisSpec axisPhi{phiBins, phiMin, phiMax, "#it{#varphi} (rad)"};

    const AxisSpec axisSel{9, 0.5, 9.5, "Selection"};
    histos.add("eventSelection", "Event Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Position");

    histos.add("trackSelection", "Track Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(1, "Tracks read");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(2, "Passed #it{p}_{T}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(3, "Passed #it{#eta}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(4, "Passed #it{#varphi}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(5, "Passed quality cuts");

    histos.add("trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    histos.add("pt/num", "Numerator " + tagPt, kTH1D, {axisPt});
    histos.add("pt/den", "Denominator " + tagPt, kTH1D, {axisPt});

    histos.add("eta/num", "Numerator " + tagEta, kTH1D, {axisEta});
    histos.add("eta/den", "Denominator " + tagEta, kTH1D, {axisEta});

    histos.add("phi/num", "Numerator " + tagPhi, kTH1D, {axisPhi});
    histos.add("phi/den", "Denominator " + tagPhi, kTH1D, {axisPhi});

    // Pos
    histos.add("ptpos/num", "Numerator Positve " + tagPt, kTH1D, {axisPt});
    histos.add("ptpos/den", "Denominator Positve " + tagPt, kTH1D, {axisPt});

    histos.add("etapos/num", "Numerator Positve " + tagEta, kTH1D, {axisEta});
    histos.add("etapos/den", "Denominator Positve " + tagEta, kTH1D, {axisEta});

    histos.add("phipos/num", "Numerator Positve " + tagPhi, kTH1D, {axisPhi});
    histos.add("phipos/den", "Denominator Positve " + tagPhi, kTH1D, {axisPhi});

    // Neg
    histos.add("ptneg/num", "Numerator Negative " + tagPt, kTH1D, {axisPt});
    histos.add("ptneg/den", "Denominator Negative " + tagPt, kTH1D, {axisPt});

    histos.add("etaneg/num", "Numerator Negative " + tagEta, kTH1D, {axisEta});
    histos.add("etaneg/den", "Denominator Negative " + tagEta, kTH1D, {axisEta});

    histos.add("phineg/num", "Numerator Negative " + tagPhi, kTH1D, {axisPhi});
    histos.add("phineg/den", "Denominator Negative " + tagPhi, kTH1D, {axisPhi});

    list.setObject(new TList);
    if (makeEff) {
      auto makeEfficiency = [&](TString effname, TString efftitle, auto templateHisto) {
        TAxis* axis = histos.get<TH1>(templateHisto)->GetXaxis();
        if (axis->IsVariableBinSize()) {
          list->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
        } else {
          list->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
        }
      };
      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHistoX, auto templateHistoY) {
        TAxis* axisX = histos.get<TH1>(templateHistoX)->GetXaxis();
        TAxis* axisY = histos.get<TH1>(templateHistoY)->GetXaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          list->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
        } else {
          list->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
        }
      };
      makeEfficiency("efficiencyVsPt", "Efficiency in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("pt/num"));
      makeEfficiency("efficiencyVsP", "Efficiency in data " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST("pt/num"));
      makeEfficiency("efficiencyVsEta", "Efficiency in data " + tagEta + ";#it{#eta};Efficiency", HIST("eta/num"));
      makeEfficiency("efficiencyVsPhi", "Efficiency in data " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST("phi/num"));

      makeEfficiency2D("efficiencyVsPtVsEta", Form("Efficiency in data #it{#varphi} [%.2f,%.2f];%s;%s;Efficiency", phiMin.value, phiMax.value, "#it{p}_{T} (GeV/#it{c})", "#it{#eta}"), HIST("pt/num"), HIST("eta/num"));
      makeEfficiency2D("efficiencyVsPtVsPhi", Form("Efficiency in data #it{#eta} [%.2f,%.2f];%s;%s;Efficiency", etaMin.value, etaMax.value, "#it{p}_{T} (GeV/#it{c})", "#it{#varphi} (rad)"), HIST("pt/num"), HIST("phi/num"));
    }
  }

  void process(const o2::aod::Collision& collision,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection>& tracks)
  {

    histos.fill(HIST("eventSelection"), 1);
    if (collision.numContrib() < nMinNumberOfContributors) {
      return;
    }
    histos.fill(HIST("eventSelection"), 2);
    if ((collision.posZ() < vertexZMin || collision.posZ() > vertexZMax)) {
      return;
    }
    histos.fill(HIST("eventSelection"), 3);

    for (const auto& track : tracks) {
      histos.fill(HIST("trackSelection"), 1);
      if ((track.pt() < ptMin || track.pt() > ptMax)) { // Check pt
        continue;
      }
      histos.fill(HIST("trackSelection"), 2);
      if ((track.eta() < etaMin || track.eta() > etaMax)) { // Check eta
        continue;
      }
      histos.fill(HIST("trackSelection"), 3);
      if ((track.phi() < phiMin || track.phi() > phiMax)) { // Check phi
        continue;
      }
      histos.fill(HIST("trackSelection"), 4);
      if (!track.isGlobalTrack()) { // Check general cuts
        continue;
      }
      histos.fill(HIST("trackSelection"), 5);

      histos.fill(HIST("trackLength"), track.length());
      histos.fill(HIST("pt/den"), track.pt());
      histos.fill(HIST("eta/den"), track.eta());
      histos.fill(HIST("phi/den"), track.phi());
      if (track.hasTOF()) {
        histos.fill(HIST("pt/num"), track.pt());
        histos.fill(HIST("eta/num"), track.eta());
        histos.fill(HIST("phi/num"), track.phi());
        if (track.sign() > 0) {
          histos.fill(HIST("ptpos/num"), track.pt());
          histos.fill(HIST("etapos/num"), track.eta());
          histos.fill(HIST("phipos/num"), track.phi());
        } else {
          histos.fill(HIST("ptneg/num"), track.pt());
          histos.fill(HIST("etaneg/num"), track.eta());
          histos.fill(HIST("phineg/num"), track.phi());
        }
      }

      if (makeEff) {
        static_cast<TEfficiency*>(list->At(0))->Fill(track.hasTOF(), track.pt());
        static_cast<TEfficiency*>(list->At(1))->Fill(track.hasTOF(), track.p());
        static_cast<TEfficiency*>(list->At(2))->Fill(track.hasTOF(), track.eta());
        static_cast<TEfficiency*>(list->At(3))->Fill(track.hasTOF(), track.phi());
        static_cast<TEfficiency*>(list->At(4))->Fill(track.hasTOF(), track.pt(), track.eta());
        static_cast<TEfficiency*>(list->At(5))->Fill(track.hasTOF(), track.pt(), track.phi());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  // Data
  if (cfgc.options().get<int>("eff-data")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiencyData>(cfgc));
  }
  // Sign blind
  if (cfgc.options().get<int>("eff-el")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Electron, 0>>(cfgc, TaskName{"qa-tracking-efficiency-electron"}));
  }
  if (cfgc.options().get<int>("eff-mu")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Muon, 0>>(cfgc, TaskName{"qa-tracking-efficiency-muon"}));
  }
  if (cfgc.options().get<int>("eff-pi")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Pion, 0>>(cfgc, TaskName{"qa-tracking-efficiency-pion"}));
  }
  if (cfgc.options().get<int>("eff-ka")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Kaon, 0>>(cfgc, TaskName{"qa-tracking-efficiency-kaon"}));
  }
  if (cfgc.options().get<int>("eff-pr")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Proton, 0>>(cfgc, TaskName{"qa-tracking-efficiency-proton"}));
  }
  if (cfgc.options().get<int>("eff-de")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Deuteron, 0>>(cfgc, TaskName{"qa-tracking-efficiency-deuteron"}));
  }
  if (cfgc.options().get<int>("eff-tr")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Triton, 0>>(cfgc, TaskName{"qa-tracking-efficiency-triton"}));
  }
  if (cfgc.options().get<int>("eff-he")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Helium3, 0>>(cfgc, TaskName{"qa-tracking-efficiency-helium3"}));
  }
  // Pos
  if (cfgc.options().get<int>("eff-el-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Electron, 1>>(cfgc, TaskName{"qa-tracking-efficiency-electron-pos"}));
  }
  if (cfgc.options().get<int>("eff-mu-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Muon, 1>>(cfgc, TaskName{"qa-tracking-efficiency-muon-pos"}));
  }
  if (cfgc.options().get<int>("eff-pi-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Pion, 1>>(cfgc, TaskName{"qa-tracking-efficiency-pion-pos"}));
  }
  if (cfgc.options().get<int>("eff-ka-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Kaon, 1>>(cfgc, TaskName{"qa-tracking-efficiency-kaon-pos"}));
  }
  if (cfgc.options().get<int>("eff-pr-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Proton, 1>>(cfgc, TaskName{"qa-tracking-efficiency-proton-pos"}));
  }
  if (cfgc.options().get<int>("eff-de-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Deuteron, 1>>(cfgc, TaskName{"qa-tracking-efficiency-deuteron-pos"}));
  }
  if (cfgc.options().get<int>("eff-tr-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Triton, 1>>(cfgc, TaskName{"qa-tracking-efficiency-triton-pos"}));
  }
  if (cfgc.options().get<int>("eff-he-pos")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Helium3, 1>>(cfgc, TaskName{"qa-tracking-efficiency-helium3-pos"}));
  }
  // Neg
  if (cfgc.options().get<int>("eff-el-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Electron, -1>>(cfgc, TaskName{"qa-tracking-efficiency-electron-neg"}));
  }
  if (cfgc.options().get<int>("eff-mu-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Muon, -1>>(cfgc, TaskName{"qa-tracking-efficiency-muon-neg"}));
  }
  if (cfgc.options().get<int>("eff-pi-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Pion, -1>>(cfgc, TaskName{"qa-tracking-efficiency-pion-neg"}));
  }
  if (cfgc.options().get<int>("eff-ka-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Kaon, -1>>(cfgc, TaskName{"qa-tracking-efficiency-kaon-neg"}));
  }
  if (cfgc.options().get<int>("eff-pr-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Proton, -1>>(cfgc, TaskName{"qa-tracking-efficiency-proton-neg"}));
  }
  if (cfgc.options().get<int>("eff-de-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Deuteron, -1>>(cfgc, TaskName{"qa-tracking-efficiency-deuteron-neg"}));
  }
  if (cfgc.options().get<int>("eff-tr-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Triton, -1>>(cfgc, TaskName{"qa-tracking-efficiency-triton-neg"}));
  }
  if (cfgc.options().get<int>("eff-he-neg")) {
    w.push_back(adaptAnalysisTask<QaTrackingEfficiency<o2::track::PID::Helium3, -1>>(cfgc, TaskName{"qa-tracking-efficiency-helium3-neg"}));
  }
  return w;
}
