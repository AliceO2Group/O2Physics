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
    {"eff-mc", VariantType::Int, 0, {"Efficiency for the MC"}},
    {"eff-mc-pos", VariantType::Int, 0, {"Efficiency for the Electron Positive PDG code"}},
    {"eff-mc-neg", VariantType::Int, 0, {"Efficiency for the Electron Negative PDG code"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to QA the efficiency of a particular particle defined by its pdg code
// template <o2::track::pid_constants::ID particle, int pdgSign>
template <int pdgSign>
struct QaEfficiencyMc {
  static constexpr int nSpecies = 9;
  static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040};
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
  Configurable<bool> doEl{"do-el", false, "Flag to run with the PDG code of pions"};
  Configurable<bool> doMu{"do-mu", false, "Flag to run with the PDG code of muons"};
  Configurable<bool> doPi{"do-pi", false, "Flag to run with the PDG code of pions"};
  Configurable<bool> doKa{"do-ka", false, "Flag to run with the PDG code of kaons"};
  Configurable<bool> doPr{"do-pr", false, "Flag to run with the PDG code of protons"};
  Configurable<bool> doDe{"do-de", false, "Flag to run with the PDG code of deuterons"};
  Configurable<bool> doTr{"do-tr", false, "Flag to run with the PDG code of tritons"};
  Configurable<bool> doHe{"do-he", false, "Flag to run with the PDG code of helium 3"};
  Configurable<bool> doAl{"do-al", false, "Flag to run with the PDG code of helium 4"};
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

  OutputObj<TList> listEfficiency{"Efficiency"};
  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Pt
  static constexpr std::string_view hptnum[nSpecies] = {"el/pt/num", "mu/pt/num", "pi/pt/num",
                                                        "ka/pt/num", "pr/pt/num", "de/pt/num",
                                                        "tr/pt/num", "he/pt/num", "al/pt/num"};
  static constexpr std::string_view hptnumtrk[nSpecies] = {"el/pt/numtrk", "mu/pt/numtrk", "pi/pt/numtrk",
                                                           "ka/pt/numtrk", "pr/pt/numtrk", "de/pt/numtrk",
                                                           "tr/pt/numtrk", "he/pt/numtrk", "al/pt/numtrk"};
  static constexpr std::string_view hptnumtof[nSpecies] = {"el/pt/numtof", "mu/pt/numtof", "pi/pt/numtof",
                                                           "ka/pt/numtof", "pr/pt/numtof", "de/pt/numtof",
                                                           "tr/pt/numtof", "he/pt/numtof", "al/pt/numtof"};
  static constexpr std::string_view hptden[nSpecies] = {"el/pt/den", "mu/pt/den", "pi/pt/den",
                                                        "ka/pt/den", "pr/pt/den", "de/pt/den",
                                                        "tr/pt/den", "he/pt/den", "al/pt/den"};
  // P
  static constexpr std::string_view hpnum[nSpecies] = {"el/p/num", "mu/p/num", "pi/p/num",
                                                       "ka/p/num", "pr/p/num", "de/p/num",
                                                       "tr/p/num", "he/p/num", "al/p/num"};
  static constexpr std::string_view hpnumtrk[nSpecies] = {"el/p/numtrk", "mu/p/numtrk", "pi/p/numtrk",
                                                          "ka/p/numtrk", "pr/p/numtrk", "de/p/numtrk",
                                                          "tr/p/numtrk", "he/p/numtrk", "al/p/numtrk"};
  static constexpr std::string_view hpnumtof[nSpecies] = {"el/p/numtof", "mu/p/numtof", "pi/p/numtof",
                                                          "ka/p/numtof", "pr/p/numtof", "de/p/numtof",
                                                          "tr/p/numtof", "he/p/numtof", "al/p/numtof"};
  static constexpr std::string_view hpden[nSpecies] = {"el/p/den", "mu/p/den", "pi/p/den",
                                                       "ka/p/den", "pr/p/den", "de/p/den",
                                                       "tr/p/den", "he/p/den", "al/p/den"};
  // Eta
  static constexpr std::string_view hetanum[nSpecies] = {"el/eta/num", "mu/eta/num", "pi/eta/num",
                                                         "ka/eta/num", "pr/eta/num", "de/eta/num",
                                                         "tr/eta/num", "he/eta/num", "al/eta/num"};
  static constexpr std::string_view hetanumtrk[nSpecies] = {"el/eta/numtrk", "mu/eta/numtrk", "pi/eta/numtrk",
                                                            "ka/eta/numtrk", "pr/eta/numtrk", "de/eta/numtrk",
                                                            "tr/eta/numtrk", "he/eta/numtrk", "al/eta/numtrk"};
  static constexpr std::string_view hetanumtof[nSpecies] = {"el/eta/numtof", "mu/eta/numtof", "pi/eta/numtof",
                                                            "ka/eta/numtof", "pr/eta/numtof", "de/eta/numtof",
                                                            "tr/eta/numtof", "he/eta/numtof", "al/eta/numtof"};
  static constexpr std::string_view hetaden[nSpecies] = {"el/eta/den", "mu/eta/den", "pi/eta/den",
                                                         "ka/eta/den", "pr/eta/den", "de/eta/den",
                                                         "tr/eta/den", "he/eta/den", "al/eta/den"};
  // Y
  static constexpr std::string_view hynum[nSpecies] = {"el/y/num", "mu/y/num", "pi/y/num",
                                                       "ka/y/num", "pr/y/num", "de/y/num",
                                                       "tr/y/num", "he/y/num", "al/y/num"};
  static constexpr std::string_view hynumtof[nSpecies] = {"el/y/numtof", "mu/y/numtof", "pi/y/numtof",
                                                          "ka/y/numtof", "pr/y/numtof", "de/y/numtof",
                                                          "tr/y/numtof", "he/y/numtof", "al/y/numtof"};
  static constexpr std::string_view hyden[nSpecies] = {"el/y/den", "mu/y/den", "pi/y/den",
                                                       "ka/y/den", "pr/y/den", "de/y/den",
                                                       "tr/y/den", "he/y/den", "al/y/den"};
  // Phi
  static constexpr std::string_view hphinum[nSpecies] = {"el/phi/num", "mu/phi/num", "pi/phi/num",
                                                         "ka/phi/num", "pr/phi/num", "de/phi/num",
                                                         "tr/phi/num", "he/phi/num", "al/phi/num"};
  static constexpr std::string_view hphinumtrk[nSpecies] = {"el/phi/numtrk", "mu/phi/numtrk", "pi/phi/numtrk",
                                                            "ka/phi/numtrk", "pr/phi/numtrk", "de/phi/numtrk",
                                                            "tr/phi/numtrk", "he/phi/numtrk", "al/phi/numtrk"};
  static constexpr std::string_view hphinumtof[nSpecies] = {"el/phi/numtof", "mu/phi/numtof", "pi/phi/numtof",
                                                            "ka/phi/numtof", "pr/phi/numtof", "de/phi/numtof",
                                                            "tr/phi/numtof", "he/phi/numtof", "al/phi/numtof"};
  static constexpr std::string_view hphiden[nSpecies] = {"el/phi/den", "mu/phi/den", "pi/phi/den",
                                                         "ka/phi/den", "pr/phi/den", "de/phi/den",
                                                         "tr/phi/den", "he/phi/den", "al/phi/den"};
  // Eta-Phi
  static constexpr std::string_view hptetanum[nSpecies] = {"el/pteta/num", "mu/pteta/num", "pi/pteta/num",
                                                           "ka/pteta/num", "pr/pteta/num", "de/pteta/num",
                                                           "tr/pteta/num", "he/pteta/num", "al/pteta/num"};
  static constexpr std::string_view hptetaden[nSpecies] = {"el/pteta/den", "mu/pteta/den", "pi/pteta/den",
                                                           "ka/pteta/den", "pr/pteta/den", "de/pteta/den",
                                                           "tr/pteta/den", "he/pteta/den", "al/pteta/den"};

  template <o2::track::PID::ID particle>
  void makeHistograms()
  {
    AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisP{ptBins, ptMin, ptMax, "#it{p} (GeV/#it{c})"};
    if (logPt) {
      axisPt.makeLogaritmic();
      axisP.makeLogaritmic();
    }
    const AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};
    const AxisSpec axisY{yBins, yMin, yMax, "#it{y}"};
    const AxisSpec axisPhi{phiBins, phiMin, phiMax, "#it{#varphi} (rad)"};

    const TString tagPt = Form("%s #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f] Prim %i",
                               o2::track::PID::getName(particle),
                               etaMin.value, etaMax.value,
                               yMin.value, yMax.value,
                               phiMin.value, phiMax.value,
                               selPrim.value);
    histos.add(hptnum[particle].data(), "Numerator " + tagPt, kTH1D, {axisPt});
    histos.add(hptnumtrk[particle].data(), "Numerator Track " + tagPt, kTH1D, {axisPt});
    histos.add(hptnumtof[particle].data(), "Numerator TOF " + tagPt, kTH1D, {axisPt});
    histos.add(hptden[particle].data(), "Denominator " + tagPt, kTH1D, {axisPt});

    histos.add(hpnum[particle].data(), "Numerator " + tagPt, kTH1D, {axisP});
    histos.add(hpnumtrk[particle].data(), "Numerator Track " + tagPt, kTH1D, {axisP});
    histos.add(hpnumtof[particle].data(), "Numerator TOF " + tagPt, kTH1D, {axisP});
    histos.add(hpden[particle].data(), "Denominator " + tagPt, kTH1D, {axisP});

    const TString tagEta = Form("%s #it{p}_{T} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f] Prim %i",
                                o2::track::PID::getName(particle),
                                ptMin.value, ptMax.value,
                                yMin.value, yMax.value,
                                phiMin.value, phiMax.value,
                                selPrim.value);
    histos.add(hetanum[particle].data(), "Numerator " + tagEta, kTH1D, {axisEta});
    histos.add(hetanumtrk[particle].data(), "Numerator Track " + tagEta, kTH1D, {axisEta});
    histos.add(hetanumtof[particle].data(), "Numerator TOF " + tagEta, kTH1D, {axisEta});
    histos.add(hetaden[particle].data(), "Denominator " + tagEta, kTH1D, {axisEta});

    const TString tagY = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f] Prim %i",
                              o2::track::PID::getName(particle),
                              ptMin.value, ptMax.value,
                              etaMin.value, etaMax.value,
                              phiMin.value, phiMax.value,
                              selPrim.value);
    histos.add(hynum[particle].data(), "Numerator " + tagY, kTH1D, {axisY});
    histos.add(hynumtof[particle].data(), "Numerator TOF " + tagY, kTH1D, {axisY});
    histos.add(hyden[particle].data(), "Denominator " + tagY, kTH1D, {axisY});

    const TString tagPhi = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] Prim %i",
                                o2::track::PID::getName(particle),
                                ptMin.value, ptMax.value,
                                etaMin.value, etaMax.value,
                                yMin.value, yMax.value,
                                selPrim.value);
    histos.add(hphinum[particle].data(), "Numerator " + tagPhi, kTH1D, {axisPhi});
    histos.add(hphinumtrk[particle].data(), "Numerator Track " + tagPhi, kTH1D, {axisPhi});
    histos.add(hphinumtof[particle].data(), "Numerator TOF " + tagPhi, kTH1D, {axisPhi});
    histos.add(hphiden[particle].data(), "Denominator " + tagPhi, kTH1D, {axisPhi});

    const TString tagPtEta = Form("%s #it{#varphi} [%.2f,%.2f] #it{y} [%.2f,%.2f] Prim %i",
                                  o2::track::PID::getName(particle),
                                  phiMin.value, phiMax.value,
                                  yMin.value, yMax.value,
                                  selPrim.value);
    histos.add(hptetanum[particle].data(), "Numerator " + tagPtEta, kTH2D, {axisPt, axisEta});
    histos.add(hptetaden[particle].data(), "Denominator " + tagPtEta, kTH2D, {axisPt, axisEta});

    if (makeEff) {
      TList* subList = new TList();
      subList->SetName(o2::track::PID::getName(particle));
      listEfficiency->Add(subList);
      auto makeEfficiency = [&](TString effname, TString efftitle, auto templateHisto) {
        effname = o2::track::PID::getName(particle) + effname;
        const TAxis* axis = histos.get<TH1>(templateHisto)->GetXaxis();
        if (axis->IsVariableBinSize()) {
          subList->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
        } else {
          subList->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
        }
      };
      makeEfficiency("efficiencyVsPt", "Efficiency " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST(hptnum[particle]));
      makeEfficiency("efficiencyVsP", "Efficiency " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST(hpnum[particle]));
      makeEfficiency("efficiencyVsEta", "Efficiency " + tagEta + ";#it{#eta};Efficiency", HIST(hetanum[particle]));
      makeEfficiency("efficiencyVsPhi", "Efficiency " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST(hphinum[particle]));

      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHisto) {
        effname = o2::track::PID::getName(particle) + effname;
        const TAxis* axisX = histos.get<TH2>(templateHisto)->GetXaxis();
        const TAxis* axisY = histos.get<TH2>(templateHisto)->GetYaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          subList->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
        } else {
          subList->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
        }
      };
      makeEfficiency2D("efficiencyVsPtVsEta", "Efficiency " + tagPtEta + ";#it{#varphi} (rad);Efficiency", HIST(hptetanum[particle]));
    }
  }

  void init(InitContext&)
  {
    const AxisSpec axisSel{30, 0.5, 30.5, "Selection"};
    histos.add("eventSelection", "Event Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Position");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Ev. Sel.");

    histos.add("trackSelection", "Track Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(1, "Tracks read");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(2, "Passed has MC part.");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(3, "Passed Ev. Reco.");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(4, "Passed #it{p}_{T}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(5, "Passed #it{#eta}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(6, "Passed #it{#varphi}");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(7, "Passed y");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(8, "Passed Prim.");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(9, "Passed Prim. MC");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(10, "Passed Fake");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(11, "Passed standard quality cuts");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(12, "Passed has collision");
    for (int i = 0; i < nSpecies; i++) {
      histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(13 + i, Form("Passed PDG %i", PDGs[i]));
    }
    histos.add("fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1D, {{1, 0, 1}});

    histos.add("partSelection", "Particle Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(1, "Particles read");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Reco.");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(6, "Passed y");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(7, "Passed Prim.");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(8, "Passed Prim. MC");
    for (int i = 0; i < nSpecies; i++) {
      histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(9 + i, Form("Passed PDG %i", PDGs[i]));
    }

    histos.add("eventMultiplicity", "Event Selection", kTH1D, {{1000, 0, 5000}});
    histos.add("trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    listEfficiency.setObject(new TList);
    if (doEl) {
      makeHistograms<o2::track::PID::Electron>();
    }
    if (doMu) {
      makeHistograms<o2::track::PID::Muon>();
    }
    if (doPi) {
      makeHistograms<o2::track::PID::Pion>();
    }
    if (doKa) {
      makeHistograms<o2::track::PID::Kaon>();
    }
    if (doPr) {
      makeHistograms<o2::track::PID::Proton>();
    }
    if (doDe) {
      makeHistograms<o2::track::PID::Deuteron>();
    }
    if (doTr) {
      makeHistograms<o2::track::PID::Triton>();
    }
    if (doHe) {
      makeHistograms<o2::track::PID::Helium3>();
    }
    if (doAl) {
      makeHistograms<o2::track::PID::Alpha>();
    }
  }

  template <o2::track::PID::ID particle, typename particleType>
  bool isPdgSelected(particleType mcParticle)
  {
    // Selecting PDG code
    if constexpr (pdgSign == 0) {
      if (abs(mcParticle.pdgCode()) != PDGs[particle]) {
        return false;
      }
    } else if constexpr (pdgSign == 1) {
      if (mcParticle.pdgCode() != PDGs[particle]) {
        return false;
      }
    } else if constexpr (pdgSign == -1) {
      if (mcParticle.pdgCode() != -PDGs[particle]) {
        return false;
      }
    }
    return true;
  }

  template <o2::track::PID::ID particle, typename trackType>
  void fillTrackHistograms(const trackType& track)
  {
    const auto mcParticle = track.mcParticle();

    if (!isPdgSelected<particle>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("trackSelection"), 13 + particle);

    histos.fill(HIST(hpnum[particle]), mcParticle.p());
    histos.fill(HIST(hptnum[particle]), mcParticle.pt());
    histos.fill(HIST(hetanum[particle]), mcParticle.eta());
    histos.fill(HIST(hynum[particle]), mcParticle.y());
    histos.fill(HIST(hphinum[particle]), mcParticle.phi());
    histos.fill(HIST(hptetanum[particle]), mcParticle.pt(), mcParticle.eta());

    histos.fill(HIST(hpnumtrk[particle]), track.p());
    histos.fill(HIST(hptnumtrk[particle]), track.pt());
    histos.fill(HIST(hetanumtrk[particle]), track.eta());
    histos.fill(HIST(hphinumtrk[particle]), track.phi());
    if (!track.hasTOF()) {
      return;
    }
    histos.fill(HIST(hpnumtof[particle]), mcParticle.p());
    histos.fill(HIST(hptnumtof[particle]), mcParticle.pt());
    histos.fill(HIST(hetanumtof[particle]), mcParticle.eta());
    histos.fill(HIST(hynumtof[particle]), mcParticle.y());
    histos.fill(HIST(hphinumtof[particle]), mcParticle.phi());
  }

  template <o2::track::PID::ID particle, typename particleType>
  void fillParticleHistograms(const particleType& mcParticle)
  {
    if (!isPdgSelected<particle>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST(hpden[particle]), mcParticle.p());
    histos.fill(HIST(hptden[particle]), mcParticle.pt());
    histos.fill(HIST(hetaden[particle]), mcParticle.eta());
    histos.fill(HIST(hyden[particle]), mcParticle.y());
    histos.fill(HIST(hphiden[particle]), mcParticle.phi());
    histos.fill(HIST(hptetaden[particle]), mcParticle.pt(), mcParticle.eta());
  }

  template <o2::track::PID::ID particle>
  void fillEfficiency()
  {

    TList* subList = static_cast<TList*>(listEfficiency->FindObject(o2::track::PID::getName(particle)));
    if (!subList) {
      return;
    }

    auto fillEfficiency = [&](int index, auto num, auto den) {
      static_cast<TEfficiency*>(subList->At(index))->SetTotalHistogram(*histos.get<TH1>(den).get(), "f");
      static_cast<TEfficiency*>(subList->At(index))->SetPassedHistogram(*histos.get<TH1>(num).get(), "f");
    };
    fillEfficiency(0, HIST(hptnum[particle]), HIST(hptden[particle]));
    fillEfficiency(1, HIST(hpnum[particle]), HIST(hpden[particle]));
    fillEfficiency(2, HIST(hetanum[particle]), HIST(hetaden[particle]));
    fillEfficiency(3, HIST(hphinum[particle]), HIST(hphiden[particle]));

    auto fillEfficiency2D = [&](int index, auto num, auto den) {
      static_cast<TEfficiency*>(subList->At(index))->SetTotalHistogram(*histos.get<TH2>(den).get(), "f");
      static_cast<TEfficiency*>(subList->At(index))->SetPassedHistogram(*histos.get<TH2>(num).get(), "f");
    };
    fillEfficiency2D(4, HIST(hptetanum[particle]), HIST(hptetaden[particle]));
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

    auto rejectParticle = [&](const auto& p, auto h, const int& offset = 0) {
      histos.fill(h, 1 + offset);
      const auto evtReconstructed = std::find(recoEvt.begin(), recoEvt.end(), p.mcCollision().globalIndex()) != recoEvt.end();
      if (!evtReconstructed) { // Check that the event is reconstructed
        return true;
      }

      histos.fill(h, 2 + offset);
      if ((p.pt() < ptMin || p.pt() > ptMax)) { // Check pt
        return true;
      }
      histos.fill(h, 3 + offset);
      if ((p.eta() < etaMin || p.eta() > etaMax)) { // Check eta
        return true;
      }
      histos.fill(h, 4 + offset);
      if ((p.phi() < phiMin || p.phi() > phiMax)) { // Check phi
        return true;
      }
      histos.fill(h, 5 + offset);
      if ((p.y() < yMin || p.y() > yMax)) { // Check rapidity
        return true;
      }
      histos.fill(h, 6 + offset);
      if ((selPrim == 1) && !p.isPhysicalPrimary()) { // Requiring is physical primary
        return true;
      }
      histos.fill(h, 7 + offset);

      // Select primaries based on position
      const float dx = p.vx() - p.mcCollision().posX();
      const float dy = p.vy() - p.mcCollision().posY();
      const float dz = p.vz() - p.mcCollision().posZ();
      if (sqrt(dx * dx + dy * dy + dz * dz) > 0.0001) {
        return true;
      }
      histos.fill(h, 8 + offset);

      return false;
    };

    // Track loop
    for (const auto& track : tracks) {
      histos.fill(HIST("trackSelection"), 1);
      if (!track.has_mcParticle()) {
        histos.fill(HIST("fakeTrackNoiseHits"), 0.5);
        continue;
      }
      const auto mcParticle = track.mcParticle();
      if (rejectParticle(mcParticle, HIST("trackSelection"), 1)) {
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
      // Filling variable histograms
      histos.fill(HIST("trackLength"), track.length());
      if (doEl) {
        fillTrackHistograms<o2::track::PID::Electron>(track);
      }
      if (doMu) {
        fillTrackHistograms<o2::track::PID::Muon>(track);
      }
      if (doPi) {
        fillTrackHistograms<o2::track::PID::Pion>(track);
      }
      if (doKa) {
        fillTrackHistograms<o2::track::PID::Kaon>(track);
      }
      if (doPr) {
        fillTrackHistograms<o2::track::PID::Proton>(track);
      }
      if (doDe) {
        fillTrackHistograms<o2::track::PID::Deuteron>(track);
      }
      if (doTr) {
        fillTrackHistograms<o2::track::PID::Triton>(track);
      }
      if (doHe) {
        fillTrackHistograms<o2::track::PID::Helium3>(track);
      }
      if (doAl) {
        fillTrackHistograms<o2::track::PID::Alpha>(track);
      }
    }

    float dNdEta = 0;
    for (const auto& mcParticle : mcParticles) {
      if (TMath::Abs(mcParticle.eta()) <= 2.f && !mcParticle.has_daughters()) {
        dNdEta += 1.f;
      }
      if (rejectParticle(mcParticle, HIST("partSelection"))) {
        continue;
      }

      if (doEl) {
        fillParticleHistograms<o2::track::PID::Electron>(mcParticle);
      }
      if (doMu) {
        fillParticleHistograms<o2::track::PID::Muon>(mcParticle);
      }
      if (doPi) {
        fillParticleHistograms<o2::track::PID::Pion>(mcParticle);
      }
      if (doKa) {
        fillParticleHistograms<o2::track::PID::Kaon>(mcParticle);
      }
      if (doPr) {
        fillParticleHistograms<o2::track::PID::Proton>(mcParticle);
      }
      if (doDe) {
        fillParticleHistograms<o2::track::PID::Deuteron>(mcParticle);
      }
      if (doTr) {
        fillParticleHistograms<o2::track::PID::Triton>(mcParticle);
      }
      if (doHe) {
        fillParticleHistograms<o2::track::PID::Helium3>(mcParticle);
      }
      if (doAl) {
        fillParticleHistograms<o2::track::PID::Alpha>(mcParticle);
      }
    }
    histos.fill(HIST("eventMultiplicity"), dNdEta * 0.5f / 2.f);

    if (makeEff) {
      if (doEl) {
        fillEfficiency<o2::track::PID::Electron>();
      }
      if (doMu) {
        fillEfficiency<o2::track::PID::Muon>();
      }
      if (doPi) {
        fillEfficiency<o2::track::PID::Pion>();
      }
      if (doKa) {
        fillEfficiency<o2::track::PID::Kaon>();
      }
      if (doPr) {
        fillEfficiency<o2::track::PID::Proton>();
      }
      if (doDe) {
        fillEfficiency<o2::track::PID::Deuteron>();
      }
      if (doTr) {
        fillEfficiency<o2::track::PID::Triton>();
      }
      if (doHe) {
        fillEfficiency<o2::track::PID::Helium3>();
      }
      if (doAl) {
        fillEfficiency<o2::track::PID::Alpha>();
      }
    }
  }
};

/// Task to QA the efficiency of a particular particle defined by its pdg code
struct QaEfficiencyData {
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

  OutputObj<TList> listEfficiency{"Efficiency"};
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

    const TString tagEtaPhi = Form("#it{p}_{T} [%.2f,%.2f]",
                                   ptMin.value, ptMax.value);

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

    // Sum
    histos.add("sum/pt/num", "Numerator " + tagPt, kTH1D, {axisPt});
    histos.add("sum/pt/den", "Denominator " + tagPt, kTH1D, {axisPt});

    histos.add("sum/eta/num", "Numerator " + tagEta, kTH1D, {axisEta});
    histos.add("sum/eta/den", "Denominator " + tagEta, kTH1D, {axisEta});

    histos.add("sum/phi/num", "Numerator " + tagPhi, kTH1D, {axisPhi});
    histos.add("sum/phi/den", "Denominator " + tagPhi, kTH1D, {axisPhi});

    histos.add("sum/etaphi/num", "Numerator " + tagPhi, kTH2D, {axisEta, axisPhi});
    histos.add("sum/etaphi/den", "Denominator " + tagPhi, kTH2D, {axisEta, axisPhi});

    // Pos
    histos.add("pos/pt/num", "Numerator Positve " + tagPt, kTH1D, {axisPt});
    histos.add("pos/pt/den", "Denominator Positve " + tagPt, kTH1D, {axisPt});

    histos.add("pos/eta/num", "Numerator Positve " + tagEta, kTH1D, {axisEta});
    histos.add("pos/eta/den", "Denominator Positve " + tagEta, kTH1D, {axisEta});

    histos.add("pos/phi/num", "Numerator Positve " + tagPhi, kTH1D, {axisPhi});
    histos.add("pos/phi/den", "Denominator Positve " + tagPhi, kTH1D, {axisPhi});

    histos.add("pos/etaphi/num", "Numerator Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("pos/etaphi/den", "Denominator Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // Neg
    histos.add("neg/pt/num", "Numerator Negative " + tagPt, kTH1D, {axisPt});
    histos.add("neg/pt/den", "Denominator Negative " + tagPt, kTH1D, {axisPt});

    histos.add("neg/eta/num", "Numerator Negative " + tagEta, kTH1D, {axisEta});
    histos.add("neg/eta/den", "Denominator Negative " + tagEta, kTH1D, {axisEta});

    histos.add("neg/phi/num", "Numerator Negative " + tagPhi, kTH1D, {axisPhi});
    histos.add("neg/phi/den", "Denominator Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("neg/etaphi/num", "Numerator Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("neg/etaphi/den", "Denominator Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    listEfficiency.setObject(new TList);
    if (makeEff) {
      auto makeEfficiency = [&](TString effname, TString efftitle, auto templateHisto) {
        TAxis* axis = histos.get<TH1>(templateHisto)->GetXaxis();
        if (axis->IsVariableBinSize()) {
          listEfficiency->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
        } else {
          listEfficiency->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
        }
      };
      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHistoX, auto templateHistoY) {
        TAxis* axisX = histos.get<TH1>(templateHistoX)->GetXaxis();
        TAxis* axisY = histos.get<TH1>(templateHistoY)->GetYaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          listEfficiency->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
        } else {
          listEfficiency->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
        }
      };
      makeEfficiency("efficiencyVsPt", "Efficiency in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("sum/pt/num"));
      makeEfficiency("efficiencyVsP", "Efficiency in data " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST("sum/pt/num"));
      makeEfficiency("efficiencyVsEta", "Efficiency in data " + tagEta + ";#it{#eta};Efficiency", HIST("sum/eta/num"));
      makeEfficiency("efficiencyVsPhi", "Efficiency in data " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST("sum/phi/num"));

      makeEfficiency2D("efficiencyVsPtVsEta", Form("Efficiency in data #it{#varphi} [%.2f,%.2f];%s;%s;Efficiency", phiMin.value, phiMax.value, "#it{p}_{T} (GeV/#it{c})", "#it{#eta}"), HIST("sum/pt/num"), HIST("sum/eta/num"));
      makeEfficiency2D("efficiencyVsPtVsPhi", Form("Efficiency in data #it{#eta} [%.2f,%.2f];%s;%s;Efficiency", etaMin.value, etaMax.value, "#it{p}_{T} (GeV/#it{c})", "#it{#varphi} (rad)"), HIST("sum/pt/num"), HIST("sum/phi/num"));
    }
  }

  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
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
      histos.fill(HIST("sum/pt/den"), track.pt());
      histos.fill(HIST("sum/eta/den"), track.eta());
      histos.fill(HIST("sum/phi/den"), track.phi());
      histos.fill(HIST("sum/etaphi/den"), track.eta(), track.phi());
      if (track.sign() > 0) {
        histos.fill(HIST("pos/pt/den"), track.pt());
        histos.fill(HIST("pos/eta/den"), track.eta());
        histos.fill(HIST("pos/phi/den"), track.phi());
        histos.fill(HIST("pos/etaphi/den"), track.eta(), track.phi());
      } else {
        histos.fill(HIST("neg/pt/den"), track.pt());
        histos.fill(HIST("neg/eta/den"), track.eta());
        histos.fill(HIST("neg/phi/den"), track.phi());
        histos.fill(HIST("neg/etaphi/den"), track.eta(), track.phi());
      }
      if (track.hasTOF()) {
        histos.fill(HIST("sum/pt/num"), track.pt());
        histos.fill(HIST("sum/eta/num"), track.eta());
        histos.fill(HIST("sum/phi/num"), track.phi());
        histos.fill(HIST("sum/etaphi/num"), track.eta(), track.phi());
        if (track.sign() > 0) {
          histos.fill(HIST("pos/pt/num"), track.pt());
          histos.fill(HIST("pos/eta/num"), track.eta());
          histos.fill(HIST("pos/phi/num"), track.phi());
          histos.fill(HIST("pos/etaphi/num"), track.eta(), track.phi());
        } else {
          histos.fill(HIST("neg/pt/num"), track.pt());
          histos.fill(HIST("neg/eta/num"), track.eta());
          histos.fill(HIST("neg/phi/num"), track.phi());
          histos.fill(HIST("neg/etaphi/num"), track.eta(), track.phi());
        }
      }

      if (makeEff) {
        static_cast<TEfficiency*>(listEfficiency->At(0))->Fill(track.hasTOF(), track.pt());
        static_cast<TEfficiency*>(listEfficiency->At(1))->Fill(track.hasTOF(), track.p());
        static_cast<TEfficiency*>(listEfficiency->At(2))->Fill(track.hasTOF(), track.eta());
        static_cast<TEfficiency*>(listEfficiency->At(3))->Fill(track.hasTOF(), track.phi());
        static_cast<TEfficiency*>(listEfficiency->At(4))->Fill(track.hasTOF(), track.pt(), track.eta());
        static_cast<TEfficiency*>(listEfficiency->At(5))->Fill(track.hasTOF(), track.pt(), track.phi());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  // Data
  if (cfgc.options().get<int>("eff-data")) {
    w.push_back(adaptAnalysisTask<QaEfficiencyData>(cfgc));
  }
  if (cfgc.options().get<int>("eff-mc")) { // Sign blind
    w.push_back(adaptAnalysisTask<QaEfficiencyMc<0>>(cfgc, TaskName{"qa-efficiency"}));
  }
  if (cfgc.options().get<int>("eff-mc-pos")) { // Pos
    w.push_back(adaptAnalysisTask<QaEfficiencyMc<1>>(cfgc, TaskName{"qa-efficiency-pos"}));
  }
  if (cfgc.options().get<int>("eff-mc-neg")) { // Neg
    w.push_back(adaptAnalysisTask<QaEfficiencyMc<-1>>(cfgc, TaskName{"qa-efficiency-neg"}));
  }
  return w;
}
