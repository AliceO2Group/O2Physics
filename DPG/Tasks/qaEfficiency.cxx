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
    {"eff", VariantType::Int, 1, {"Efficiency for the sum"}},
    {"eff-pos", VariantType::Int, 0, {"Efficiency for the Positive"}},
    {"eff-neg", VariantType::Int, 0, {"Efficiency for the Negative"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to QA the efficiency of a particular particle defined by its pdg code
template <int pdgSign>
struct QaEfficiency {
  // Particle information
  static constexpr int nSpecies = o2::track::PID::NIDs + 1;
  static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha", "All"};
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040, 0};
  // Track/particle selection
  Configurable<float> etaMin{"eta-min", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 3.f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 6.284f, "Upper limit in phi"};
  Configurable<float> yMin{"y-min", -0.5f, "Lower limit in y"};
  Configurable<float> yMax{"y-max", 0.5f, "Upper limit in y"};
  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 5.f, "Upper limit in pT"};
  // Particle only selection
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
  // Track only selection
  Configurable<bool> applyTrackSelection{"applyTrackSelection", false, "Flag to use the standard track selection when selecting numerator and denominator"};
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
  Configurable<int> applyEvSel{"applyEvSel", 0, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};

  OutputObj<TList> listEfficiencyMC{"EfficiencyMC"};
  OutputObj<TList> listEfficiencyData{"EfficiencyData"};
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
  void makeMCHistograms()
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
      LOG(debug) << "Making TEfficiency for MC";
      TList* subList = new TList();
      subList->SetName(partName);
      listEfficiencyMC->Add(subList);
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

  void initMC(const AxisSpec& axisSel)
  {
    if (!doprocessMC) {
      return;
    }

    auto h = histos.add<TH1>("MC/trackSelection", "Track Selection", kTH1D, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Passed has MC part.");
    h->GetXaxis()->SetBinLabel(3, "Passed Ev. Reco.");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(5, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(6, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(7, "Passed y");
    h->GetXaxis()->SetBinLabel(8, "Passed Fake");
    h->GetXaxis()->SetBinLabel(9, "Passed standard quality cuts");
    h->GetXaxis()->SetBinLabel(10, "Passed has collision");
    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(11 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1D, {{1, 0, 1}});

    h = histos.add<TH1>("MC/particleSelection", "Particle Selection", kTH1D, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Particles read");
    h->GetXaxis()->SetBinLabel(2, "Passed Ev. Reco.");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(6, "Passed y");
    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(7 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/eventMultiplicity", "Event Selection", kTH1D, {{1000, 0, 5000}});

    histos.add("MC/trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    listEfficiencyMC.setObject(new TList);
    if (doEl) {
      makeMCHistograms<o2::track::PID::Electron>();
    }
    if (doMu) {
      makeMCHistograms<o2::track::PID::Muon>();
    }
    if (doPi) {
      makeMCHistograms<o2::track::PID::Pion>();
    }
    if (doKa) {
      makeMCHistograms<o2::track::PID::Kaon>();
    }
    if (doPr) {
      makeMCHistograms<o2::track::PID::Proton>();
    }
    if (doDe) {
      makeMCHistograms<o2::track::PID::Deuteron>();
    }
    if (doTr) {
      makeMCHistograms<o2::track::PID::Triton>();
    }
    if (doHe) {
      makeMCHistograms<o2::track::PID::Helium3>();
    }
    if (doAl) {
      makeMCHistograms<o2::track::PID::Alpha>();
    }
    if (doUnId) {
      makeMCHistograms<o2::track::PID::NIDs>();
    }
  }

  void initData(const AxisSpec& axisSel)
  {
    if (!doprocessData) {
      return;
    }

    auto h = histos.add<TH1>("Data/trackSelection", "Track Selection", kTH1D, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(5, "Passed TrackType");
    h->GetXaxis()->SetBinLabel(6, "Passed PtRange");
    h->GetXaxis()->SetBinLabel(7, "Passed EtaRange");
    h->GetXaxis()->SetBinLabel(8, "Passed DCAxy");
    h->GetXaxis()->SetBinLabel(9, "Passed DCAz");
    h->GetXaxis()->SetBinLabel(10, "Passed GoldenChi2");
    h->GetXaxis()->SetBinLabel(11, "Passed quality cuts");

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

    histos.add("Data/trackLength", "Track length;Track length (cm)", kTH1D, {{2000, -1000, 1000}});

    // ITS-TPC-TOF
    histos.add("Data/sum/pt/its_tpc_tof", "ITS-TPC-TOF " + tagPt, kTH1D, {axisPt});
    histos.add("Data/pos/pt/its_tpc_tof", "ITS-TPC-TOF Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/its_tpc_tof", "ITS-TPC-TOF Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/sum/eta/its_tpc_tof", "ITS-TPC-TOF " + tagEta, kTH1D, {axisEta});
    histos.add("Data/pos/eta/its_tpc_tof", "ITS-TPC-TOF Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/its_tpc_tof", "ITS-TPC-TOF Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/sum/phi/its_tpc_tof", "ITS-TPC-TOF " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/pos/phi/its_tpc_tof", "ITS-TPC-TOF Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/its_tpc_tof", "ITS-TPC-TOF Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/sum/etaphi/its_tpc_tof", "ITS-TPC-TOF " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/pos/etaphi/its_tpc_tof", "ITS-TPC-TOF Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc_tof", "ITS-TPC-TOF Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // ITS-TPC
    histos.add("Data/sum/pt/its_tpc", "ITS-TPC " + tagPt, kTH1D, {axisPt});
    histos.add("Data/pos/pt/its_tpc", "ITS-TPC Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/its_tpc", "ITS-TPC Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/sum/eta/its_tpc", "ITS-TPC " + tagEta, kTH1D, {axisEta});
    histos.add("Data/pos/eta/its_tpc", "ITS-TPC Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/its_tpc", "ITS-TPC Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/sum/phi/its_tpc", "ITS-TPC " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/pos/phi/its_tpc", "ITS-TPC Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/its_tpc", "ITS-TPC Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/sum/etaphi/its_tpc", "ITS-TPC " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/pos/etaphi/its_tpc", "ITS-TPC Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc", "ITS-TPC Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // TPC
    histos.add("Data/sum/pt/tpc", "TPC " + tagPt, kTH1D, {axisPt});
    histos.add("Data/pos/pt/tpc", "TPC Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/tpc", "TPC Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/sum/eta/tpc", "TPC " + tagEta, kTH1D, {axisEta});
    histos.add("Data/pos/eta/tpc", "TPC Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/tpc", "TPC Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/sum/phi/tpc", "TPC " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/pos/phi/tpc", "TPC Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/tpc", "TPC Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/sum/etaphi/tpc", "TPC " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/pos/etaphi/tpc", "TPC Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/tpc", "TPC Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // ITS
    histos.add("Data/sum/pt/its", "ITS " + tagPt, kTH1D, {axisPt});
    histos.add("Data/pos/pt/its", "ITS Positive " + tagPt, kTH1D, {axisPt});
    histos.add("Data/neg/pt/its", "ITS Negative " + tagPt, kTH1D, {axisPt});

    histos.add("Data/sum/eta/its", "ITS " + tagEta, kTH1D, {axisEta});
    histos.add("Data/pos/eta/its", "ITS Positive " + tagEta, kTH1D, {axisEta});
    histos.add("Data/neg/eta/its", "ITS Negative " + tagEta, kTH1D, {axisEta});

    histos.add("Data/sum/phi/its", "ITS " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/pos/phi/its", "ITS Positive " + tagPhi, kTH1D, {axisPhi});
    histos.add("Data/neg/phi/its", "ITS Negative " + tagPhi, kTH1D, {axisPhi});

    histos.add("Data/sum/etaphi/its", "ITS " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/pos/etaphi/its", "ITS Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its", "ITS Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    listEfficiencyData.setObject(new TList);
    if (makeEff) {
      LOG(debug) << "Making TEfficiency for Data";
      auto makeEfficiency = [&](TString effname, TString efftitle, auto templateHisto) {
        TAxis* axis = histos.get<TH1>(templateHisto)->GetXaxis();
        if (axis->IsVariableBinSize()) {
          listEfficiencyData->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
        } else {
          listEfficiencyData->Add(new TEfficiency(effname, efftitle, axis->GetNbins(), axis->GetXmin(), axis->GetXmax()));
        }
      };
      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHistoX, auto templateHistoY) {
        TAxis* axisX = histos.get<TH1>(templateHistoX)->GetXaxis();
        TAxis* axisY = histos.get<TH1>(templateHistoY)->GetYaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          listEfficiencyData->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
        } else {
          listEfficiencyData->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
        }
      };
      makeEfficiency("ITSTPCMatchingEfficiencyVsPt", "ITS-TPC M.E. in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/sum/pt/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsPt", "TPC-TOF M.E. in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/sum/pt/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsP", "TPC-TOF M.E. in data " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST("Data/sum/pt/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsEta", "TPC-TOF M.E. in data " + tagEta + ";#it{#eta};Efficiency", HIST("Data/sum/eta/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsPhi", "TPC-TOF M.E. in data " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST("Data/sum/phi/its_tpc_tof"));

      makeEfficiency2D("TPCTOFMatchingEfficiencyVsPtVsEta", Form("TPC-TOF M.E. in data #it{#varphi} [%.2f,%.2f];%s;%s;Efficiency", phiMin.value, phiMax.value, "#it{p}_{T} (GeV/#it{c})", "#it{#eta}"), HIST("Data/sum/pt/its_tpc_tof"), HIST("Data/sum/eta/its_tpc_tof"));
      makeEfficiency2D("TPCTOFMatchingEfficiencyVsPtVsPhi", Form("TPC-TOF M.E. in data #it{#eta} [%.2f,%.2f];%s;%s;Efficiency", etaMin.value, etaMax.value, "#it{p}_{T} (GeV/#it{c})", "#it{#varphi} (rad)"), HIST("Data/sum/pt/its_tpc_tof"), HIST("Data/sum/phi/its_tpc_tof"));
    }
  }

  void init(InitContext&)
  {
    const AxisSpec axisSel{30, 0.5, 30.5, "Selection"};
    histos.add("eventSelection", "Event Selection", kTH1D, {axisSel});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Position");

    initData(axisSel);
    initMC(axisSel);
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
  void fillMCTrackHistograms(const trackType& track)
  {
    LOG(debug) << "Filling track histograms for id " << static_cast<int>(id);
    const auto mcParticle = track.mcParticle();

    if (!isPdgSelected<id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("MC/trackSelection"), 11 + id);

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
  void fillMCParticleHistograms(const particleType& mcParticle)
  {
    LOG(debug) << "Filling particle histograms for id " << static_cast<int>(id);
    if (!isPdgSelected<id>(mcParticle)) { // Selecting PDG code
      return;
    }
    histos.fill(HIST("MC/particleSelection"), 7 + id);

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
  void fillMCEfficiency()
  {
    if (!makeEff) {
      return;
    }
    const char* partName = id == o2::track::PID::NIDs ? "All" : o2::track::PID::getName(id);
    LOG(debug) << "Filling efficiency for particle " << static_cast<int>(id) << " " << partName;
    TList* subList = static_cast<TList*>(listEfficiencyMC->FindObject(partName));
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
  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection>& tracks)
  {
    isCollisionSelected<true>(collision);
  }

  // MC process
  void processMC(const o2::aod::McParticles& mcParticles,
                 const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>& collisions,
                 const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::McTrackLabels, o2::aod::TrackSelection>& tracks,
                 const o2::aod::McCollisions&)
  {

    std::vector<int64_t> recoEvt(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (!isCollisionSelected<false>(collision)) {
        continue;
      }
      recoEvt[nevts++] = collision.mcCollision().globalIndex();
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
      histos.fill(HIST("MC/trackSelection"), 1);
      if (!track.has_mcParticle()) {
        histos.fill(HIST("MC/fakeTrackNoiseHits"), 0.5);
        continue;
      }
      const auto mcParticle = track.mcParticle();
      if (rejectParticle(mcParticle, HIST("MC/trackSelection"), 1)) {
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

      histos.fill(HIST("MC/trackSelection"), 10);
      if (applyTrackSelection && !track.isGlobalTrack()) { // Check general cuts
        continue;
      }
      histos.fill(HIST("MC/trackSelection"), 11);
      if (!track.has_collision()) {
        continue;
      }
      histos.fill(HIST("MC/trackSelection"), 12);
      // Filling variable histograms
      histos.fill(HIST("MC/trackLength"), track.length());
      if (doEl) {
        fillMCTrackHistograms<o2::track::PID::Electron>(track);
      }
      if (doMu) {
        fillMCTrackHistograms<o2::track::PID::Muon>(track);
      }
      if (doPi) {
        fillMCTrackHistograms<o2::track::PID::Pion>(track);
      }
      if (doKa) {
        fillMCTrackHistograms<o2::track::PID::Kaon>(track);
      }
      if (doPr) {
        fillMCTrackHistograms<o2::track::PID::Proton>(track);
      }
      if (doDe) {
        fillMCTrackHistograms<o2::track::PID::Deuteron>(track);
      }
      if (doTr) {
        fillMCTrackHistograms<o2::track::PID::Triton>(track);
      }
      if (doHe) {
        fillMCTrackHistograms<o2::track::PID::Helium3>(track);
      }
      if (doAl) {
        fillMCTrackHistograms<o2::track::PID::Alpha>(track);
      }
      if (doUnId) {
        fillMCTrackHistograms<o2::track::PID::NIDs>(track);
      }
    }

    float dNdEta = 0;
    for (const auto& mcParticle : mcParticles) { // Loop on particles to fill the denominator
      if (TMath::Abs(mcParticle.eta()) <= 2.f && !mcParticle.has_daughters()) {
        dNdEta += 1.f;
      }
      if (rejectParticle(mcParticle, HIST("MC/particleSelection"))) {
        continue;
      }

      if (doEl) {
        fillMCParticleHistograms<o2::track::PID::Electron>(mcParticle);
      }
      if (doMu) {
        fillMCParticleHistograms<o2::track::PID::Muon>(mcParticle);
      }
      if (doPi) {
        fillMCParticleHistograms<o2::track::PID::Pion>(mcParticle);
      }
      if (doKa) {
        fillMCParticleHistograms<o2::track::PID::Kaon>(mcParticle);
      }
      if (doPr) {
        fillMCParticleHistograms<o2::track::PID::Proton>(mcParticle);
      }
      if (doDe) {
        fillMCParticleHistograms<o2::track::PID::Deuteron>(mcParticle);
      }
      if (doTr) {
        fillMCParticleHistograms<o2::track::PID::Triton>(mcParticle);
      }
      if (doHe) {
        fillMCParticleHistograms<o2::track::PID::Helium3>(mcParticle);
      }
      if (doAl) {
        fillMCParticleHistograms<o2::track::PID::Alpha>(mcParticle);
      }
      if (doUnId) {
        fillMCParticleHistograms<o2::track::PID::NIDs>(mcParticle);
      }
    }
    histos.fill(HIST("MC/eventMultiplicity"), dNdEta * 0.5f / 2.f);

    // Fill TEfficiencies
    if (doEl) {
      fillMCEfficiency<o2::track::PID::Electron>();
    }
    if (doMu) {
      fillMCEfficiency<o2::track::PID::Muon>();
    }
    if (doPi) {
      fillMCEfficiency<o2::track::PID::Pion>();
    }
    if (doKa) {
      fillMCEfficiency<o2::track::PID::Kaon>();
    }
    if (doPr) {
      fillMCEfficiency<o2::track::PID::Proton>();
    }
    if (doDe) {
      fillMCEfficiency<o2::track::PID::Deuteron>();
    }
    if (doTr) {
      fillMCEfficiency<o2::track::PID::Triton>();
    }
    if (doHe) {
      fillMCEfficiency<o2::track::PID::Helium3>();
    }
    if (doAl) {
      fillMCEfficiency<o2::track::PID::Alpha>();
    }
    if (doUnId) {
      fillMCEfficiency<o2::track::PID::NIDs>();
    }
  }
  PROCESS_SWITCH(QaEfficiency, processMC, "process MC", false);

  void processData(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
                   const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection>& tracks)
  {

    if (!isCollisionSelected<false>(collision)) {
      return;
    }

    bool passedITSCuts = true;
    bool passedTPCCuts = true;
    for (const auto& track : tracks) {
      histos.fill(HIST("Data/trackSelection"), 1);
      if ((track.pt() < ptMin || track.pt() > ptMax)) { // Check pt
        continue;
      }
      histos.fill(HIST("Data/trackSelection"), 2);
      if ((track.eta() < etaMin || track.eta() > etaMax)) { // Check eta
        continue;
      }
      histos.fill(HIST("Data/trackSelection"), 3);
      if ((track.phi() < phiMin || track.phi() > phiMax)) { // Check phi
        continue;
      }
      histos.fill(HIST("Data/trackSelection"), 4);
      if (applyTrackSelection) { // Check general cuts
        if (!track.passedTrackType()) {
          continue;
        }
        histos.fill(HIST("Data/trackSelection"), 5);
        if (!track.passedPtRange()) {
          continue;
        }
        histos.fill(HIST("Data/trackSelection"), 6);
        if (!track.passedEtaRange()) {
          continue;
        }
        histos.fill(HIST("Data/trackSelection"), 7);
        if (!track.passedDCAxy()) {
          continue;
        }
        histos.fill(HIST("Data/trackSelection"), 8);
        if (!track.passedDCAz()) {
          continue;
        }
        histos.fill(HIST("Data/trackSelection"), 9);
        if (!track.passedGoldenChi2()) {
          continue;
        }
        histos.fill(HIST("Data/trackSelection"), 10);

        passedITSCuts = track.passedITSNCls() &&
                        track.passedITSChi2NDF() &&
                        track.passedITSRefit() &&
                        track.passedITSHits() &&
                        track.hasITS();

        passedTPCCuts = track.passedTPCNCls() &&
                        track.passedTPCCrossedRows() &&
                        track.passedTPCCrossedRowsOverNCls() &&
                        track.passedTPCChi2NDF() &&
                        track.passedTPCRefit() &&
                        track.hasTPC();
      } else {
        passedITSCuts = track.hasITS();
        passedTPCCuts = track.hasTPC();
      }

      histos.fill(HIST("Data/trackSelection"), 11);

      histos.fill(HIST("Data/trackLength"), track.length());

      if (passedITSCuts) {
        histos.fill(HIST("Data/sum/pt/its"), track.pt());
        histos.fill(HIST("Data/sum/eta/its"), track.eta());
        histos.fill(HIST("Data/sum/phi/its"), track.phi());
        histos.fill(HIST("Data/sum/etaphi/its"), track.eta(), track.phi());
      }
      if (passedTPCCuts) {
        histos.fill(HIST("Data/sum/pt/tpc"), track.pt());
        histos.fill(HIST("Data/sum/eta/tpc"), track.eta());
        histos.fill(HIST("Data/sum/phi/tpc"), track.phi());
        histos.fill(HIST("Data/sum/etaphi/tpc"), track.eta(), track.phi());
      }
      if (passedITSCuts && passedTPCCuts) {
        histos.fill(HIST("Data/sum/pt/its_tpc"), track.pt());
        histos.fill(HIST("Data/sum/eta/its_tpc"), track.eta());
        histos.fill(HIST("Data/sum/phi/its_tpc"), track.phi());
        histos.fill(HIST("Data/sum/etaphi/its_tpc"), track.eta(), track.phi());
      }
      if (passedITSCuts && passedTPCCuts && track.hasTOF()) {
        histos.fill(HIST("Data/sum/pt/its_tpc_tof"), track.pt());
        histos.fill(HIST("Data/sum/eta/its_tpc_tof"), track.eta());
        histos.fill(HIST("Data/sum/phi/its_tpc_tof"), track.phi());
        histos.fill(HIST("Data/sum/etaphi/its_tpc_tof"), track.eta(), track.phi());
      }

      if (track.sign() > 0) {
        if (passedITSCuts) {
          histos.fill(HIST("Data/pos/pt/its"), track.pt());
          histos.fill(HIST("Data/pos/eta/its"), track.eta());
          histos.fill(HIST("Data/pos/phi/its"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/its"), track.eta(), track.phi());
        }
        if (passedTPCCuts) {
          histos.fill(HIST("Data/pos/pt/tpc"), track.pt());
          histos.fill(HIST("Data/pos/eta/tpc"), track.eta());
          histos.fill(HIST("Data/pos/phi/tpc"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/tpc"), track.eta(), track.phi());
        }
        if (passedITSCuts && passedTPCCuts) {
          histos.fill(HIST("Data/pos/pt/its_tpc"), track.pt());
          histos.fill(HIST("Data/pos/eta/its_tpc"), track.eta());
          histos.fill(HIST("Data/pos/phi/its_tpc"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/its_tpc"), track.eta(), track.phi());
        }
        if (track.hasITS() && passedTPCCuts && track.hasTOF()) {
          histos.fill(HIST("Data/pos/pt/its_tpc_tof"), track.pt());
          histos.fill(HIST("Data/pos/eta/its_tpc_tof"), track.eta());
          histos.fill(HIST("Data/pos/phi/its_tpc_tof"), track.phi());
          histos.fill(HIST("Data/pos/etaphi/its_tpc_tof"), track.eta(), track.phi());
        }
      } else {
        if (passedITSCuts) {
          histos.fill(HIST("Data/neg/pt/its"), track.pt());
          histos.fill(HIST("Data/neg/eta/its"), track.eta());
          histos.fill(HIST("Data/neg/phi/its"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/its"), track.eta(), track.phi());
        }
        if (passedTPCCuts) {
          histos.fill(HIST("Data/neg/pt/tpc"), track.pt());
          histos.fill(HIST("Data/neg/eta/tpc"), track.eta());
          histos.fill(HIST("Data/neg/phi/tpc"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/tpc"), track.eta(), track.phi());
        }
        if (passedITSCuts && passedTPCCuts) {
          histos.fill(HIST("Data/neg/pt/its_tpc"), track.pt());
          histos.fill(HIST("Data/neg/eta/its_tpc"), track.eta());
          histos.fill(HIST("Data/neg/phi/its_tpc"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/its_tpc"), track.eta(), track.phi());
        }
        if (passedITSCuts && passedTPCCuts && track.hasTOF()) {
          histos.fill(HIST("Data/neg/pt/its_tpc_tof"), track.pt());
          histos.fill(HIST("Data/neg/eta/its_tpc_tof"), track.eta());
          histos.fill(HIST("Data/neg/phi/its_tpc_tof"), track.phi());
          histos.fill(HIST("Data/neg/etaphi/its_tpc_tof"), track.eta(), track.phi());
        }
      }

      if (makeEff) {
        if (passedITSCuts) {
          static_cast<TEfficiency*>(listEfficiencyData->At(0))->Fill(passedTPCCuts, track.pt());
        }
        if (passedITSCuts && passedTPCCuts) {
          static_cast<TEfficiency*>(listEfficiencyData->At(1))->Fill(track.hasTOF(), track.pt());
          static_cast<TEfficiency*>(listEfficiencyData->At(2))->Fill(track.hasTOF(), track.p());
          static_cast<TEfficiency*>(listEfficiencyData->At(3))->Fill(track.hasTOF(), track.eta());
          static_cast<TEfficiency*>(listEfficiencyData->At(4))->Fill(track.hasTOF(), track.phi());
          static_cast<TEfficiency*>(listEfficiencyData->At(5))->Fill(track.hasTOF(), track.pt(), track.eta());
          static_cast<TEfficiency*>(listEfficiencyData->At(6))->Fill(track.hasTOF(), track.pt(), track.phi());
        }
      }
    }
  }
  PROCESS_SWITCH(QaEfficiency, processData, "process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  // Data
  if (cfgc.options().get<int>("eff")) { // Sign blind
    w.push_back(adaptAnalysisTask<QaEfficiency<0>>(cfgc, TaskName{"qa-efficiency"}));
  }
  if (cfgc.options().get<int>("eff-pos")) { // Pos
    w.push_back(adaptAnalysisTask<QaEfficiency<1>>(cfgc, TaskName{"qa-efficiency-pos"}));
  }
  if (cfgc.options().get<int>("eff-neg")) { // Neg
    w.push_back(adaptAnalysisTask<QaEfficiency<-1>>(cfgc, TaskName{"qa-efficiency-neg"}));
  }
  return w;
}
