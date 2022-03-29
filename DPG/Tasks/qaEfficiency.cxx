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
///

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
template <int pdgSign>
struct QaEfficiencyMc {
  static constexpr int nSpecies = o2::track::PID::NIDs + 1;
  static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha", "All"};
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040, 0};
  // Particle selection
  Configurable<float> etaMin{"eta-min", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 3.f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 6.284f, "Upper limit in phi"};
  Configurable<float> yMin{"y-min", -0.5f, "Lower limit in y"};
  Configurable<float> yMax{"y-max", 0.5f, "Upper limit in y"};
  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 5.f, "Upper limit in pT"};
  Configurable<bool> noFakes{"noFakes", false, "Flag to reject tracks that have fake hits"};
  Configurable<bool> doUnId{"do-un-id", true, "Flag to run without PDG code"};
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
  static constexpr std::string_view hPtNum[nSpecies] = {"el/pt/num", "mu/pt/num", "pi/pt/num",
                                                        "ka/pt/num", "pr/pt/num", "de/pt/num",
                                                        "tr/pt/num", "he/pt/num", "al/pt/num",
                                                        "all/pt/num"};
  static constexpr std::string_view hPtNumTrk[nSpecies] = {"el/pt/numtrk", "mu/pt/numtrk", "pi/pt/numtrk",
                                                           "ka/pt/numtrk", "pr/pt/numtrk", "de/pt/numtrk",
                                                           "tr/pt/numtrk", "he/pt/numtrk", "al/pt/numtrk",
                                                           "all/pt/numtrk"};
  static constexpr std::string_view hPtNumTof[nSpecies] = {"el/pt/numtof", "mu/pt/numtof", "pi/pt/numtof",
                                                           "ka/pt/numtof", "pr/pt/numtof", "de/pt/numtof",
                                                           "tr/pt/numtof", "he/pt/numtof", "al/pt/numtof",
                                                           "all/pt/numtof"};
  static constexpr std::string_view hPtDen[nSpecies] = {"el/pt/den", "mu/pt/den", "pi/pt/den",
                                                        "ka/pt/den", "pr/pt/den", "de/pt/den",
                                                        "tr/pt/den", "he/pt/den", "al/pt/den",
                                                        "all/pt/den"};
  // Pt for primaries
  static constexpr std::string_view hPtPrmNum[nSpecies] = {"el/ptprm/num", "mu/ptprm/num", "pi/ptprm/num",
                                                           "ka/ptprm/num", "pr/ptprm/num", "de/ptprm/num",
                                                           "tr/ptprm/num", "he/ptprm/num", "al/ptprm/num",
                                                           "all/ptprm/num"};
  static constexpr std::string_view hPtPrmNumTrk[nSpecies] = {"el/ptprm/numtrk", "mu/ptprm/numtrk", "pi/ptprm/numtrk",
                                                              "ka/ptprm/numtrk", "pr/ptprm/numtrk", "de/ptprm/numtrk",
                                                              "tr/ptprm/numtrk", "he/ptprm/numtrk", "al/ptprm/numtrk",
                                                              "all/ptprm/numtrk"};
  static constexpr std::string_view hPtPrmNumTof[nSpecies] = {"el/ptprm/numtof", "mu/ptprm/numtof", "pi/ptprm/numtof",
                                                              "ka/ptprm/numtof", "pr/ptprm/numtof", "de/ptprm/numtof",
                                                              "tr/ptprm/numtof", "he/ptprm/numtof", "al/ptprm/numtof",
                                                              "all/ptprm/numtof"};
  static constexpr std::string_view hPtPrmDen[nSpecies] = {"el/ptprm/den", "mu/ptprm/den", "pi/ptprm/den",
                                                           "ka/ptprm/den", "pr/ptprm/den", "de/ptprm/den",
                                                           "tr/ptprm/den", "he/ptprm/den", "al/ptprm/den",
                                                           "all/ptprm/den"};
  // Pt for secondaries from weak decay
  static constexpr std::string_view hPtDecNum[nSpecies] = {"el/ptdec/num", "mu/ptdec/num", "pi/ptdec/num",
                                                           "ka/ptdec/num", "pr/ptdec/num", "de/ptdec/num",
                                                           "tr/ptdec/num", "he/ptdec/num", "al/ptdec/num",
                                                           "all/ptdec/num"};
  static constexpr std::string_view hPtDecNumTrk[nSpecies] = {"el/ptdec/numtrk", "mu/ptdec/numtrk", "pi/ptdec/numtrk",
                                                              "ka/ptdec/numtrk", "pr/ptdec/numtrk", "de/ptdec/numtrk",
                                                              "tr/ptdec/numtrk", "he/ptdec/numtrk", "al/ptdec/numtrk",
                                                              "all/ptdec/numtrk"};
  static constexpr std::string_view hPtDecNumTof[nSpecies] = {"el/ptdec/numtof", "mu/ptdec/numtof", "pi/ptdec/numtof",
                                                              "ka/ptdec/numtof", "pr/ptdec/numtof", "de/ptdec/numtof",
                                                              "tr/ptdec/numtof", "he/ptdec/numtof", "al/ptdec/numtof",
                                                              "all/ptdec/numtof"};
  static constexpr std::string_view hPtDecDen[nSpecies] = {"el/ptdec/den", "mu/ptdec/den", "pi/ptdec/den",
                                                           "ka/ptdec/den", "pr/ptdec/den", "de/ptdec/den",
                                                           "tr/ptdec/den", "he/ptdec/den", "al/ptdec/den",
                                                           "all/ptdec/den"};
  // Pt for secondaries from material
  static constexpr std::string_view hPtMatNum[nSpecies] = {"el/ptmat/num", "mu/ptmat/num", "pi/ptmat/num",
                                                           "ka/ptmat/num", "pr/ptmat/num", "de/ptmat/num",
                                                           "tr/ptmat/num", "he/ptmat/num", "al/ptmat/num",
                                                           "all/ptmat/num"};
  static constexpr std::string_view hPtMatNumTrk[nSpecies] = {"el/ptmat/numtrk", "mu/ptmat/numtrk", "pi/ptmat/numtrk",
                                                              "ka/ptmat/numtrk", "pr/ptmat/numtrk", "de/ptmat/numtrk",
                                                              "tr/ptmat/numtrk", "he/ptmat/numtrk", "al/ptmat/numtrk",
                                                              "all/ptmat/numtrk"};
  static constexpr std::string_view hPtMatNumTof[nSpecies] = {"el/ptmat/numtof", "mu/ptmat/numtof", "pi/ptmat/numtof",
                                                              "ka/ptmat/numtof", "pr/ptmat/numtof", "de/ptmat/numtof",
                                                              "tr/ptmat/numtof", "he/ptmat/numtof", "al/ptmat/numtof",
                                                              "all/ptmat/numtof"};
  static constexpr std::string_view hPtMatDen[nSpecies] = {"el/ptmat/den", "mu/ptmat/den", "pi/ptmat/den",
                                                           "ka/ptmat/den", "pr/ptmat/den", "de/ptmat/den",
                                                           "tr/ptmat/den", "he/ptmat/den", "al/ptmat/den",
                                                           "all/ptmat/den"};
  // P
  static constexpr std::string_view hPNum[nSpecies] = {"el/p/num", "mu/p/num", "pi/p/num",
                                                       "ka/p/num", "pr/p/num", "de/p/num",
                                                       "tr/p/num", "he/p/num", "al/p/num",
                                                       "all/p/num"};
  static constexpr std::string_view hPNumTrk[nSpecies] = {"el/p/numtrk", "mu/p/numtrk", "pi/p/numtrk",
                                                          "ka/p/numtrk", "pr/p/numtrk", "de/p/numtrk",
                                                          "tr/p/numtrk", "he/p/numtrk", "al/p/numtrk",
                                                          "all/p/numtrk"};
  static constexpr std::string_view hPNumTof[nSpecies] = {"el/p/numtof", "mu/p/numtof", "pi/p/numtof",
                                                          "ka/p/numtof", "pr/p/numtof", "de/p/numtof",
                                                          "tr/p/numtof", "he/p/numtof", "al/p/numtof",
                                                          "all/p/numtof"};
  static constexpr std::string_view hPDen[nSpecies] = {"el/p/den", "mu/p/den", "pi/p/den",
                                                       "ka/p/den", "pr/p/den", "de/p/den",
                                                       "tr/p/den", "he/p/den", "al/p/den",
                                                       "all/p/den"};
  // Eta
  static constexpr std::string_view hEtaNum[nSpecies] = {"el/eta/num", "mu/eta/num", "pi/eta/num",
                                                         "ka/eta/num", "pr/eta/num", "de/eta/num",
                                                         "tr/eta/num", "he/eta/num", "al/eta/num",
                                                         "all/eta/num"};
  static constexpr std::string_view hEtaNumTrk[nSpecies] = {"el/eta/numtrk", "mu/eta/numtrk", "pi/eta/numtrk",
                                                            "ka/eta/numtrk", "pr/eta/numtrk", "de/eta/numtrk",
                                                            "tr/eta/numtrk", "he/eta/numtrk", "al/eta/numtrk",
                                                            "all/eta/numtrk"};
  static constexpr std::string_view hEtaNumTof[nSpecies] = {"el/eta/numtof", "mu/eta/numtof", "pi/eta/numtof",
                                                            "ka/eta/numtof", "pr/eta/numtof", "de/eta/numtof",
                                                            "tr/eta/numtof", "he/eta/numtof", "al/eta/numtof",
                                                            "all/eta/numtof"};
  static constexpr std::string_view hEtaDen[nSpecies] = {"el/eta/den", "mu/eta/den", "pi/eta/den",
                                                         "ka/eta/den", "pr/eta/den", "de/eta/den",
                                                         "tr/eta/den", "he/eta/den", "al/eta/den",
                                                         "all/eta/den"};
  // Y
  static constexpr std::string_view hYNum[nSpecies] = {"el/y/num", "mu/y/num", "pi/y/num",
                                                       "ka/y/num", "pr/y/num", "de/y/num",
                                                       "tr/y/num", "he/y/num", "al/y/num",
                                                       "all/y/num"};
  static constexpr std::string_view hYNumTof[nSpecies] = {"el/y/numtof", "mu/y/numtof", "pi/y/numtof",
                                                          "ka/y/numtof", "pr/y/numtof", "de/y/numtof",
                                                          "tr/y/numtof", "he/y/numtof", "al/y/numtof",
                                                          "all/y/numtof"};
  static constexpr std::string_view hYDen[nSpecies] = {"el/y/den", "mu/y/den", "pi/y/den",
                                                       "ka/y/den", "pr/y/den", "de/y/den",
                                                       "tr/y/den", "he/y/den", "al/y/den",
                                                       "all/y/den"};
  // Phi
  static constexpr std::string_view hPhiNum[nSpecies] = {"el/phi/num", "mu/phi/num", "pi/phi/num",
                                                         "ka/phi/num", "pr/phi/num", "de/phi/num",
                                                         "tr/phi/num", "he/phi/num", "al/phi/num",
                                                         "all/phi/num"};
  static constexpr std::string_view hPhiNumTrk[nSpecies] = {"el/phi/numtrk", "mu/phi/numtrk", "pi/phi/numtrk",
                                                            "ka/phi/numtrk", "pr/phi/numtrk", "de/phi/numtrk",
                                                            "tr/phi/numtrk", "he/phi/numtrk", "al/phi/numtrk",
                                                            "all/phi/numtrk"};
  static constexpr std::string_view hPhiNumTof[nSpecies] = {"el/phi/numtof", "mu/phi/numtof", "pi/phi/numtof",
                                                            "ka/phi/numtof", "pr/phi/numtof", "de/phi/numtof",
                                                            "tr/phi/numtof", "he/phi/numtof", "al/phi/numtof",
                                                            "all/phi/numtof"};
  static constexpr std::string_view hPhiDen[nSpecies] = {"el/phi/den", "mu/phi/den", "pi/phi/den",
                                                         "ka/phi/den", "pr/phi/den", "de/phi/den",
                                                         "tr/phi/den", "he/phi/den", "al/phi/den",
                                                         "all/phi/den"};
  // Eta-Phi
  static constexpr std::string_view hPtEtaNum[nSpecies] = {"el/pteta/num", "mu/pteta/num", "pi/pteta/num",
                                                           "ka/pteta/num", "pr/pteta/num", "de/pteta/num",
                                                           "tr/pteta/num", "he/pteta/num", "al/pteta/num",
                                                           "all/pteta/num"};
  static constexpr std::string_view hPtEtaDen[nSpecies] = {"el/pteta/den", "mu/pteta/den", "pi/pteta/den",
                                                           "ka/pteta/den", "pr/pteta/den", "de/pteta/den",
                                                           "tr/pteta/den", "he/pteta/den", "al/pteta/den",
                                                           "all/pteta/den"};

  template <o2::track::PID::ID id>
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

    const char* partName = id == o2::track::PID::NIDs ? "All" : o2::track::PID::getName(id);
    LOG(debug) << "Preparing histograms for particle: " << partName;

    const TString tagPt = Form("%s #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                               partName,
                               etaMin.value, etaMax.value,
                               yMin.value, yMax.value,
                               phiMin.value, phiMax.value);
    histos.add(hPtNum[id].data(), "Numerator " + tagPt, kTH1D, {axisPt});
    histos.add(hPtNumTrk[id].data(), "Numerator Track " + tagPt, kTH1D, {axisPt});
    histos.add(hPtNumTof[id].data(), "Numerator TOF " + tagPt, kTH1D, {axisPt});
    histos.add(hPtDen[id].data(), "Denominator " + tagPt, kTH1D, {axisPt});

    histos.add(hPtPrmNum[id].data(), "Numerator " + tagPt + " Primaries", kTH1D, {axisPt});
    histos.add(hPtPrmNumTrk[id].data(), "Numerator Track " + tagPt + " Primaries", kTH1D, {axisPt});
    histos.add(hPtPrmNumTof[id].data(), "Numerator TOF " + tagPt + " Primaries", kTH1D, {axisPt});
    histos.add(hPtPrmDen[id].data(), "Denominator " + tagPt + " Primaries", kTH1D, {axisPt});

    histos.add(hPtDecNum[id].data(), "Numerator " + tagPt + " Sec. from decays", kTH1D, {axisPt});
    histos.add(hPtDecNumTrk[id].data(), "Numerator Track " + tagPt + " Sec. from decays", kTH1D, {axisPt});
    histos.add(hPtDecNumTof[id].data(), "Numerator TOF " + tagPt + " Sec. from decays", kTH1D, {axisPt});
    histos.add(hPtDecDen[id].data(), "Denominator " + tagPt + " Sec. from decays", kTH1D, {axisPt});

    histos.add(hPtMatNum[id].data(), "Numerator " + tagPt + " Sec. from material", kTH1D, {axisPt});
    histos.add(hPtMatNumTrk[id].data(), "Numerator Track " + tagPt + " Sec. from material", kTH1D, {axisPt});
    histos.add(hPtMatNumTof[id].data(), "Numerator TOF " + tagPt + " Sec. from material", kTH1D, {axisPt});
    histos.add(hPtMatDen[id].data(), "Denominator " + tagPt + " Sec. from material", kTH1D, {axisPt});

    histos.add(hPNum[id].data(), "Numerator " + tagPt, kTH1D, {axisP});
    histos.add(hPNumTrk[id].data(), "Numerator Track " + tagPt, kTH1D, {axisP});
    histos.add(hPNumTof[id].data(), "Numerator TOF " + tagPt, kTH1D, {axisP});
    histos.add(hPDen[id].data(), "Denominator " + tagPt, kTH1D, {axisP});

    const TString tagEta = Form("%s #it{p}_{T} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                                partName,
                                ptMin.value, ptMax.value,
                                yMin.value, yMax.value,
                                phiMin.value, phiMax.value);
    histos.add(hEtaNum[id].data(), "Numerator " + tagEta, kTH1D, {axisEta});
    histos.add(hEtaNumTrk[id].data(), "Numerator Track " + tagEta, kTH1D, {axisEta});
    histos.add(hEtaNumTof[id].data(), "Numerator TOF " + tagEta, kTH1D, {axisEta});
    histos.add(hEtaDen[id].data(), "Denominator " + tagEta, kTH1D, {axisEta});

    const TString tagY = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                              partName,
                              ptMin.value, ptMax.value,
                              etaMin.value, etaMax.value,
                              phiMin.value, phiMax.value);
    histos.add(hYNum[id].data(), "Numerator " + tagY, kTH1D, {axisY});
    histos.add(hYNumTof[id].data(), "Numerator TOF " + tagY, kTH1D, {axisY});
    histos.add(hYDen[id].data(), "Denominator " + tagY, kTH1D, {axisY});

    const TString tagPhi = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
                                partName,
                                ptMin.value, ptMax.value,
                                etaMin.value, etaMax.value,
                                yMin.value, yMax.value);
    histos.add(hPhiNum[id].data(), "Numerator " + tagPhi, kTH1D, {axisPhi});
    histos.add(hPhiNumTrk[id].data(), "Numerator Track " + tagPhi, kTH1D, {axisPhi});
    histos.add(hPhiNumTof[id].data(), "Numerator TOF " + tagPhi, kTH1D, {axisPhi});
    histos.add(hPhiDen[id].data(), "Denominator " + tagPhi, kTH1D, {axisPhi});

    const TString tagPtEta = Form("%s #it{#varphi} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
                                  partName,
                                  phiMin.value, phiMax.value,
                                  yMin.value, yMax.value);
    histos.add(hPtEtaNum[id].data(), "Numerator " + tagPtEta, kTH2D, {axisPt, axisEta});
    histos.add(hPtEtaDen[id].data(), "Denominator " + tagPtEta, kTH2D, {axisPt, axisEta});

    if (makeEff) {
      LOG(debug) << "Making TEfficiency";
      TList* subList = new TList();
      subList->SetName(partName);
      listEfficiency->Add(subList);
      auto makeEfficiency = [&](TString effname, auto templateHisto) {
        effname = partName + effname;
        LOG(debug) << " - " << effname;
        const auto h = histos.get<TH1>(templateHisto);
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
      makeEfficiency("efficiencyVsPt", HIST(hPtNum[id]));
      makeEfficiency("efficiencyVsPtPrm", HIST(hPtPrmNum[id]));
      makeEfficiency("efficiencyVsPtDec", HIST(hPtDecNum[id]));
      makeEfficiency("efficiencyVsPtMat", HIST(hPtMatNum[id]));
      makeEfficiency("efficiencyVsP", HIST(hPNum[id]));
      makeEfficiency("efficiencyVsEta", HIST(hEtaNum[id]));
      makeEfficiency("efficiencyVsPhi", HIST(hPhiNum[id]));

      auto makeEfficiency2D = [&](TString effname, auto templateHisto) {
        effname = partName + effname;
        LOG(debug) << " - " << effname;
        const auto h = histos.get<TH2>(templateHisto);
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
      makeEfficiency2D("efficiencyVsPtVsEta", HIST(hPtEtaNum[id]));
    }
    LOG(debug) << "Done with particle: " << partName;
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
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(8, "Passed Fake");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(9, "Passed standard quality cuts");
    histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(10, "Passed has collision");
    for (int i = 0; i < nSpecies; i++) {
      histos.get<TH1>(HIST("trackSelection"))->GetXaxis()->SetBinLabel(11 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1D, {{1, 0, 1}});

    histos.add("partSelection", "Particle Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(1, "Particles read");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Reco.");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(6, "Passed y");
    for (int i = 0; i < nSpecies; i++) {
      histos.get<TH1>(HIST("partSelection"))->GetXaxis()->SetBinLabel(7 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
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
    if (doUnId) {
      makeHistograms<o2::track::PID::NIDs>();
    }
  }

  template <o2::track::PID::ID id, typename particleType>
  bool isPdgSelected(particleType mcParticle)
  {
    // Selecting PDG code
    if constexpr (PDGs[id] == 0) {
      if constexpr (pdgSign == 0) {
        return true;
      } else if constexpr (pdgSign == 1) {
        return mcParticle.pdgCode() > 0;
      }
      if constexpr (pdgSign == -1) {
        return mcParticle.pdgCode() < 0;
      }
    }
    if constexpr (pdgSign == 0) {
      if (abs(mcParticle.pdgCode()) != PDGs[id]) {
        return false;
      }
    } else if constexpr (pdgSign == 1) {
      if (mcParticle.pdgCode() != PDGs[id]) {
        return false;
      }
    } else if constexpr (pdgSign == -1) {
      if (mcParticle.pdgCode() != -PDGs[id]) {
        return false;
      }
    }
    return true;
  }

  template <o2::track::PID::ID id, typename trackType>
  void fillTrackHistograms(const trackType& track)
  {
    LOG(debug) << "Filling track histograms for id " << static_cast<int>(id);
    const auto mcParticle = track.mcParticle();

    if (!isPdgSelected<id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("trackSelection"), 11 + id);

    histos.fill(HIST(hPNum[id]), mcParticle.p());
    histos.fill(HIST(hPtNum[id]), mcParticle.pt());
    histos.fill(HIST(hEtaNum[id]), mcParticle.eta());
    histos.fill(HIST(hYNum[id]), mcParticle.y());
    histos.fill(HIST(hPhiNum[id]), mcParticle.phi());
    histos.fill(HIST(hPtEtaNum[id]), mcParticle.pt(), mcParticle.eta());

    histos.fill(HIST(hPNumTrk[id]), track.p());
    histos.fill(HIST(hPtNumTrk[id]), track.pt());
    histos.fill(HIST(hEtaNumTrk[id]), track.eta());
    histos.fill(HIST(hPhiNumTrk[id]), track.phi());

    if (mcParticle.isPhysicalPrimary()) {
      histos.fill(HIST(hPtPrmNum[id]), mcParticle.pt());
      histos.fill(HIST(hPtPrmNumTrk[id]), track.pt());
      if (track.hasTOF()) {
        histos.fill(HIST(hPtPrmNumTof[id]), mcParticle.pt());
      }
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        histos.fill(HIST(hPtDecNum[id]), mcParticle.pt());
        histos.fill(HIST(hPtDecNumTrk[id]), track.pt());
        if (track.hasTOF()) {
          histos.fill(HIST(hPtDecNumTof[id]), mcParticle.pt());
        }
      } else { // Material
        histos.fill(HIST(hPtMatNum[id]), mcParticle.pt());
        histos.fill(HIST(hPtMatNumTrk[id]), track.pt());
        if (track.hasTOF()) {
          histos.fill(HIST(hPtMatNumTof[id]), mcParticle.pt());
        }
      }
    }
    if (!track.hasTOF()) {
      return;
    }
    histos.fill(HIST(hPNumTof[id]), mcParticle.p());
    histos.fill(HIST(hPtNumTof[id]), mcParticle.pt());
    histos.fill(HIST(hEtaNumTof[id]), mcParticle.eta());
    histos.fill(HIST(hYNumTof[id]), mcParticle.y());
    histos.fill(HIST(hPhiNumTof[id]), mcParticle.phi());
  }

  template <o2::track::PID::ID id, typename particleType>
  void fillParticleHistograms(const particleType& mcParticle)
  {
    LOG(debug) << "Filling particle histograms for id " << static_cast<int>(id);
    if (!isPdgSelected<id>(mcParticle)) { // Selecting PDG code
      return;
    }
    histos.fill(HIST("partSelection"), 7 + id);

    histos.fill(HIST(hPDen[id]), mcParticle.p());
    histos.fill(HIST(hPtDen[id]), mcParticle.pt());

    if (mcParticle.isPhysicalPrimary()) {
      histos.fill(HIST(hPtPrmDen[id]), mcParticle.pt());
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        histos.fill(HIST(hPtDecDen[id]), mcParticle.pt());
      } else { // Material
        histos.fill(HIST(hPtMatDen[id]), mcParticle.pt());
      }
    }

    histos.fill(HIST(hEtaDen[id]), mcParticle.eta());
    histos.fill(HIST(hYDen[id]), mcParticle.y());
    histos.fill(HIST(hPhiDen[id]), mcParticle.phi());
    histos.fill(HIST(hPtEtaDen[id]), mcParticle.pt(), mcParticle.eta());
  }

  template <o2::track::PID::ID id>
  void fillEfficiency()
  {
    if (!makeEff) {
      return;
    }
    const char* partName = id == o2::track::PID::NIDs ? "All" : o2::track::PID::getName(id);
    LOG(debug) << "Filling efficiency for particle " << static_cast<int>(id) << " " << partName;
    TList* subList = static_cast<TList*>(listEfficiency->FindObject(partName));
    if (!subList) {
      LOG(warning) << "Cannot find list of efficiency objects for particle " << partName;
      return;
    }

    auto doFillEfficiency = [&](TString effname, auto num, auto den) {
      effname = partName + effname;
      TEfficiency* eff = static_cast<TEfficiency*>(subList->FindObject(effname));
      if (!eff) {
        LOG(warning) << "Cannot find TEfficiency " << effname;
        return;
      }
      eff->SetTotalHistogram(*histos.get<TH1>(den).get(), "f");
      eff->SetPassedHistogram(*histos.get<TH1>(num).get(), "f");
    };

    doFillEfficiency("efficiencyVsPt", HIST(hPtNum[id]), HIST(hPtDen[id]));
    doFillEfficiency("efficiencyVsPtPrm", HIST(hPtPrmNum[id]), HIST(hPtPrmDen[id]));
    doFillEfficiency("efficiencyVsPtDec", HIST(hPtDecNum[id]), HIST(hPtDecDen[id]));
    doFillEfficiency("efficiencyVsPtMat", HIST(hPtMatNum[id]), HIST(hPtMatDen[id]));
    doFillEfficiency("efficiencyVsP", HIST(hPNum[id]), HIST(hPDen[id]));
    doFillEfficiency("efficiencyVsEta", HIST(hEtaNum[id]), HIST(hEtaDen[id]));
    doFillEfficiency("efficiencyVsPhi", HIST(hPhiNum[id]), HIST(hPhiDen[id]));
    auto fillEfficiency2D = [&](TString effname, auto num, auto den) {
      effname = partName + effname;
      TEfficiency* eff = static_cast<TEfficiency*>(subList->FindObject(effname));
      if (!eff) {
        LOG(warning) << "Cannot find TEfficiency " << effname;
        return;
      }
      eff->SetTotalHistogram(*histos.get<TH2>(den).get(), "f");
      eff->SetPassedHistogram(*histos.get<TH2>(num).get(), "f");
    };
    fillEfficiency2D("efficiencyVsPtVsEta", HIST(hPtEtaNum[id]), HIST(hPtEtaDen[id]));
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
      if (doUnId) {
        fillTrackHistograms<o2::track::PID::NIDs>(track);
      }
    }

    float dNdEta = 0;
    for (const auto& mcParticle : mcParticles) { // Loop on particles to fill the denominator
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
      if (doUnId) {
        fillParticleHistograms<o2::track::PID::NIDs>(mcParticle);
      }
    }
    histos.fill(HIST("eventMultiplicity"), dNdEta * 0.5f / 2.f);

    // Fill TEfficiencies
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
    if (doUnId) {
      fillEfficiency<o2::track::PID::NIDs>();
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
  Configurable<int> applyEvSel{"applyEvSel", 0, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};

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
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Position");

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
    if (applyEvSel == 1 && !collision.sel7()) {
      return;
    } else if (applyEvSel == 2 && !collision.sel8()) {
      return;
    }
    histos.fill(HIST("eventSelection"), 2);
    if (collision.numContrib() < nMinNumberOfContributors) {
      return;
    }
    histos.fill(HIST("eventSelection"), 3);
    if ((collision.posZ() < vertexZMin || collision.posZ() > vertexZMax)) {
      return;
    }
    histos.fill(HIST("eventSelection"), 4);

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
    w.push_back(adaptAnalysisTask<QaEfficiencyMc<0>>(cfgc, TaskName{"qa-efficiency-mc"}));
  }
  if (cfgc.options().get<int>("eff-mc-pos")) { // Pos
    w.push_back(adaptAnalysisTask<QaEfficiencyMc<1>>(cfgc, TaskName{"qa-efficiency-mc-pos"}));
  }
  if (cfgc.options().get<int>("eff-mc-neg")) { // Neg
    w.push_back(adaptAnalysisTask<QaEfficiencyMc<-1>>(cfgc, TaskName{"qa-efficiency-mc-neg"}));
  }
  return w;
}
