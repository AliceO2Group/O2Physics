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
#include "Framework/runDataProcessing.h"
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

struct QaEfficiency {
  // Particle information
  static constexpr int nSpecies = o2::track::PID::NIDs + 1; // One per PDG + 1 for unidentified
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
  Configurable<bool> noFakes{"noFakes", false, "Flag to reject tracks that have fake hits"};
  // Charge selection
  Configurable<bool> doSumPDG{"doSumPDG", true, "Flag to fill histograms for summed PDG codes. Required to fill the efficiencies"};
  Configurable<bool> doPositivePDG{"doPositivePDG", false, "Flag to fill histograms for positive PDG codes."};
  Configurable<bool> doNegativePDG{"doNegativePDG", false, "Flag to fill histograms for negative PDG codes."};
  // Particle only selection
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
  static constexpr int nHistograms = nSpecies * 3;
  // Pt
  static constexpr std::string_view hPtNum[nHistograms] = {"MC/el/sum/pt/num", "MC/mu/sum/pt/num", "MC/pi/sum/pt/num",
                                                           "MC/ka/sum/pt/num", "MC/pr/sum/pt/num", "MC/de/sum/pt/num",
                                                           "MC/tr/sum/pt/num", "MC/he/sum/pt/num", "MC/al/sum/pt/num",
                                                           "MC/all/sum/pt/num",
                                                           "MC/el/pos/pt/num", "MC/mu/pos/pt/num", "MC/pi/pos/pt/num",
                                                           "MC/ka/pos/pt/num", "MC/pr/pos/pt/num", "MC/de/pos/pt/num",
                                                           "MC/tr/pos/pt/num", "MC/he/pos/pt/num", "MC/al/pos/pt/num",
                                                           "MC/all/pos/pt/num",
                                                           "MC/el/neg/pt/num", "MC/mu/neg/pt/num", "MC/pi/neg/pt/num",
                                                           "MC/ka/neg/pt/num", "MC/pr/neg/pt/num", "MC/de/neg/pt/num",
                                                           "MC/tr/neg/pt/num", "MC/he/neg/pt/num", "MC/al/neg/pt/num",
                                                           "MC/all/neg/pt/num"};
  static constexpr std::string_view hPtNumTrk[nHistograms] = {"MC/el/sum/pt/numtrk", "MC/mu/sum/pt/numtrk", "MC/pi/sum/pt/numtrk",
                                                              "MC/ka/sum/pt/numtrk", "MC/pr/sum/pt/numtrk", "MC/de/sum/pt/numtrk",
                                                              "MC/tr/sum/pt/numtrk", "MC/he/sum/pt/numtrk", "MC/al/sum/pt/numtrk",
                                                              "MC/all/sum/pt/numtrk",
                                                              "MC/el/pos/pt/numtrk", "MC/mu/pos/pt/numtrk", "MC/pi/pos/pt/numtrk",
                                                              "MC/ka/pos/pt/numtrk", "MC/pr/pos/pt/numtrk", "MC/de/pos/pt/numtrk",
                                                              "MC/tr/pos/pt/numtrk", "MC/he/pos/pt/numtrk", "MC/al/pos/pt/numtrk",
                                                              "MC/all/pos/pt/numtrk",
                                                              "MC/el/neg/pt/numtrk", "MC/mu/neg/pt/numtrk", "MC/pi/neg/pt/numtrk",
                                                              "MC/ka/neg/pt/numtrk", "MC/pr/neg/pt/numtrk", "MC/de/neg/pt/numtrk",
                                                              "MC/tr/neg/pt/numtrk", "MC/he/neg/pt/numtrk", "MC/al/neg/pt/numtrk",
                                                              "MC/all/neg/pt/numtrk"};
  static constexpr std::string_view hPtNumTof[nHistograms] = {"MC/el/sum/pt/numtof", "MC/mu/sum/pt/numtof", "MC/pi/sum/pt/numtof",
                                                              "MC/ka/sum/pt/numtof", "MC/pr/sum/pt/numtof", "MC/de/sum/pt/numtof",
                                                              "MC/tr/sum/pt/numtof", "MC/he/sum/pt/numtof", "MC/al/sum/pt/numtof",
                                                              "MC/all/sum/pt/numtof",
                                                              "MC/el/pos/pt/numtof", "MC/mu/pos/pt/numtof", "MC/pi/pos/pt/numtof",
                                                              "MC/ka/pos/pt/numtof", "MC/pr/pos/pt/numtof", "MC/de/pos/pt/numtof",
                                                              "MC/tr/pos/pt/numtof", "MC/he/pos/pt/numtof", "MC/al/pos/pt/numtof",
                                                              "MC/all/pos/pt/numtof",
                                                              "MC/el/sum/pt/numtof", "MC/mu/sum/pt/numtof", "MC/pi/sum/pt/numtof",
                                                              "MC/ka/sum/pt/numtof", "MC/pr/sum/pt/numtof", "MC/de/sum/pt/numtof",
                                                              "MC/tr/sum/pt/numtof", "MC/he/sum/pt/numtof", "MC/al/sum/pt/numtof",
                                                              "MC/all/sum/pt/numtof"};
  static constexpr std::string_view hPtDen[nHistograms] = {"MC/el/sum/pt/den", "MC/mu/sum/pt/den", "MC/pi/sum/pt/den",
                                                           "MC/ka/sum/pt/den", "MC/pr/sum/pt/den", "MC/de/sum/pt/den",
                                                           "MC/tr/sum/pt/den", "MC/he/sum/pt/den", "MC/al/sum/pt/den",
                                                           "MC/all/sum/pt/den",
                                                           "MC/el/pos/pt/den", "MC/mu/pos/pt/den", "MC/pi/pos/pt/den",
                                                           "MC/ka/pos/pt/den", "MC/pr/pos/pt/den", "MC/de/pos/pt/den",
                                                           "MC/tr/pos/pt/den", "MC/he/pos/pt/den", "MC/al/pos/pt/den",
                                                           "MC/all/pos/pt/den",
                                                           "MC/el/neg/pt/den", "MC/mu/neg/pt/den", "MC/pi/neg/pt/den",
                                                           "MC/ka/neg/pt/den", "MC/pr/neg/pt/den", "MC/de/neg/pt/den",
                                                           "MC/tr/neg/pt/den", "MC/he/neg/pt/den", "MC/al/neg/pt/den",
                                                           "MC/all/neg/pt/den"};
  // Pt for primaries
  static constexpr std::string_view hPtPrmNum[nHistograms] = {"MC/el/sum/prm/pt/num", "MC/mu/sum/prm/pt/num", "MC/pi/sum/prm/pt/num",
                                                              "MC/ka/sum/prm/pt/num", "MC/pr/sum/prm/pt/num", "MC/de/sum/prm/pt/num",
                                                              "MC/tr/sum/prm/pt/num", "MC/he/sum/prm/pt/num", "MC/al/sum/prm/pt/num",
                                                              "MC/all/sum/prm/pt/num",
                                                              "MC/el/pos/prm/pt/num", "MC/mu/pos/prm/pt/num", "MC/pi/pos/prm/pt/num",
                                                              "MC/ka/pos/prm/pt/num", "MC/pr/pos/prm/pt/num", "MC/de/pos/prm/pt/num",
                                                              "MC/tr/pos/prm/pt/num", "MC/he/pos/prm/pt/num", "MC/al/pos/prm/pt/num",
                                                              "MC/all/pos/prm/pt/num",
                                                              "MC/el/den/prm/pt/num", "MC/mu/den/prm/pt/num", "MC/pi/den/prm/pt/num",
                                                              "MC/ka/den/prm/pt/num", "MC/pr/den/prm/pt/num", "MC/de/den/prm/pt/num",
                                                              "MC/tr/den/prm/pt/num", "MC/he/den/prm/pt/num", "MC/al/den/prm/pt/num",
                                                              "MC/all/den/prm/pt/num"};
  static constexpr std::string_view hPtPrmNumTrk[nHistograms] = {"MC/el/sum/prm/pt/numtrk", "MC/mu/sum/prm/pt/numtrk", "MC/pi/sum/prm/pt/numtrk",
                                                                 "MC/ka/sum/prm/pt/numtrk", "MC/pr/sum/prm/pt/numtrk", "MC/de/sum/prm/pt/numtrk",
                                                                 "MC/tr/sum/prm/pt/numtrk", "MC/he/sum/prm/pt/numtrk", "MC/al/sum/prm/pt/numtrk",
                                                                 "MC/all/sum/prm/pt/numtrk",
                                                                 "MC/el/pos/prm/pt/numtrk", "MC/mu/pos/prm/pt/numtrk", "MC/pi/pos/prm/pt/numtrk",
                                                                 "MC/ka/pos/prm/pt/numtrk", "MC/pr/pos/prm/pt/numtrk", "MC/de/pos/prm/pt/numtrk",
                                                                 "MC/tr/pos/prm/pt/numtrk", "MC/he/pos/prm/pt/numtrk", "MC/al/pos/prm/pt/numtrk",
                                                                 "MC/all/pos/prm/pt/numtrk",
                                                                 "MC/el/neg/prm/pt/numtrk", "MC/mu/neg/prm/pt/numtrk", "MC/pi/neg/prm/pt/numtrk",
                                                                 "MC/ka/neg/prm/pt/numtrk", "MC/pr/neg/prm/pt/numtrk", "MC/de/neg/prm/pt/numtrk",
                                                                 "MC/tr/neg/prm/pt/numtrk", "MC/he/neg/prm/pt/numtrk", "MC/al/neg/prm/pt/numtrk",
                                                                 "MC/all/neg/prm/pt/numtrk"};
  static constexpr std::string_view hPtPrmNumTof[nHistograms] = {"MC/el/sum/prm/pt/numtof", "MC/mu/sum/prm/pt/numtof", "MC/pi/sum/prm/pt/numtof",
                                                                 "MC/ka/sum/prm/pt/numtof", "MC/pr/sum/prm/pt/numtof", "MC/de/sum/prm/pt/numtof",
                                                                 "MC/tr/sum/prm/pt/numtof", "MC/he/sum/prm/pt/numtof", "MC/al/sum/prm/pt/numtof",
                                                                 "MC/all/sum/prm/pt/numtof",
                                                                 "MC/el/pos/prm/pt/numtof", "MC/mu/pos/prm/pt/numtof", "MC/pi/pos/prm/pt/numtof",
                                                                 "MC/ka/pos/prm/pt/numtof", "MC/pr/pos/prm/pt/numtof", "MC/de/pos/prm/pt/numtof",
                                                                 "MC/tr/pos/prm/pt/numtof", "MC/he/pos/prm/pt/numtof", "MC/al/pos/prm/pt/numtof",
                                                                 "MC/all/pos/prm/pt/numtof",
                                                                 "MC/el/neg/prm/pt/numtof", "MC/mu/neg/prm/pt/numtof", "MC/pi/neg/prm/pt/numtof",
                                                                 "MC/ka/neg/prm/pt/numtof", "MC/pr/neg/prm/pt/numtof", "MC/de/neg/prm/pt/numtof",
                                                                 "MC/tr/neg/prm/pt/numtof", "MC/he/neg/prm/pt/numtof", "MC/al/neg/prm/pt/numtof",
                                                                 "MC/all/neg/prm/pt/numtof"};
  static constexpr std::string_view hPtPrmDen[nHistograms] = {"MC/el/sum/prm/pt/den", "MC/mu/sum/prm/pt/den", "MC/pi/sum/prm/pt/den",
                                                              "MC/ka/sum/prm/pt/den", "MC/pr/sum/prm/pt/den", "MC/de/sum/prm/pt/den",
                                                              "MC/tr/sum/prm/pt/den", "MC/he/sum/prm/pt/den", "MC/al/sum/prm/pt/den",
                                                              "MC/all/sum/prm/pt/den",
                                                              "MC/el/pos/prm/pt/den", "MC/mu/pos/prm/pt/den", "MC/pi/pos/prm/pt/den",
                                                              "MC/ka/pos/prm/pt/den", "MC/pr/pos/prm/pt/den", "MC/de/pos/prm/pt/den",
                                                              "MC/tr/pos/prm/pt/den", "MC/he/pos/prm/pt/den", "MC/al/pos/prm/pt/den",
                                                              "MC/all/pos/prm/pt/den",
                                                              "MC/el/neg/prm/pt/den", "MC/mu/neg/prm/pt/den", "MC/pi/neg/prm/pt/den",
                                                              "MC/ka/neg/prm/pt/den", "MC/pr/neg/prm/pt/den", "MC/de/neg/prm/pt/den",
                                                              "MC/tr/neg/prm/pt/den", "MC/he/neg/prm/pt/den", "MC/al/neg/prm/pt/den",
                                                              "MC/all/neg/prm/pt/den"};
  // Pt for secondaries from weak decay
  static constexpr std::string_view hPtDecNum[nHistograms] = {"MC/el/sum/dec/pt/num", "MC/mu/sum/dec/pt/num", "MC/pi/sum/dec/pt/num",
                                                              "MC/ka/sum/dec/pt/num", "MC/pr/sum/dec/pt/num", "MC/de/sum/dec/pt/num",
                                                              "MC/tr/sum/dec/pt/num", "MC/he/sum/dec/pt/num", "MC/al/sum/dec/pt/num",
                                                              "MC/all/sum/dec/pt/num",
                                                              "MC/el/pos/dec/pt/num", "MC/mu/pos/dec/pt/num", "MC/pi/pos/dec/pt/num",
                                                              "MC/ka/pos/dec/pt/num", "MC/pr/pos/dec/pt/num", "MC/de/pos/dec/pt/num",
                                                              "MC/tr/pos/dec/pt/num", "MC/he/pos/dec/pt/num", "MC/al/pos/dec/pt/num",
                                                              "MC/all/pos/dec/pt/num",
                                                              "MC/el/den/dec/pt/num", "MC/mu/den/dec/pt/num", "MC/pi/den/dec/pt/num",
                                                              "MC/ka/den/dec/pt/num", "MC/pr/den/dec/pt/num", "MC/de/den/dec/pt/num",
                                                              "MC/tr/den/dec/pt/num", "MC/he/den/dec/pt/num", "MC/al/den/dec/pt/num",
                                                              "MC/all/den/dec/pt/num"};
  static constexpr std::string_view hPtDecNumTrk[nHistograms] = {"MC/el/sum/dec/pt/numtrk", "MC/mu/sum/dec/pt/numtrk", "MC/pi/sum/dec/pt/numtrk",
                                                                 "MC/ka/sum/dec/pt/numtrk", "MC/pr/sum/dec/pt/numtrk", "MC/de/sum/dec/pt/numtrk",
                                                                 "MC/tr/sum/dec/pt/numtrk", "MC/he/sum/dec/pt/numtrk", "MC/al/sum/dec/pt/numtrk",
                                                                 "MC/all/sum/dec/pt/numtrk",
                                                                 "MC/el/pos/dec/pt/numtrk", "MC/mu/pos/dec/pt/numtrk", "MC/pi/pos/dec/pt/numtrk",
                                                                 "MC/ka/pos/dec/pt/numtrk", "MC/pr/pos/dec/pt/numtrk", "MC/de/pos/dec/pt/numtrk",
                                                                 "MC/tr/pos/dec/pt/numtrk", "MC/he/pos/dec/pt/numtrk", "MC/al/pos/dec/pt/numtrk",
                                                                 "MC/all/pos/dec/pt/numtrk",
                                                                 "MC/el/neg/dec/pt/numtrk", "MC/mu/neg/dec/pt/numtrk", "MC/pi/neg/dec/pt/numtrk",
                                                                 "MC/ka/neg/dec/pt/numtrk", "MC/pr/neg/dec/pt/numtrk", "MC/de/neg/dec/pt/numtrk",
                                                                 "MC/tr/neg/dec/pt/numtrk", "MC/he/neg/dec/pt/numtrk", "MC/al/neg/dec/pt/numtrk",
                                                                 "MC/all/neg/dec/pt/numtrk"};
  static constexpr std::string_view hPtDecNumTof[nHistograms] = {"MC/el/sum/dec/pt/numtof", "MC/mu/sum/dec/pt/numtof", "MC/pi/sum/dec/pt/numtof",
                                                                 "MC/ka/sum/dec/pt/numtof", "MC/pr/sum/dec/pt/numtof", "MC/de/sum/dec/pt/numtof",
                                                                 "MC/tr/sum/dec/pt/numtof", "MC/he/sum/dec/pt/numtof", "MC/al/sum/dec/pt/numtof",
                                                                 "MC/all/sum/dec/pt/numtof",
                                                                 "MC/el/pos/dec/pt/numtof", "MC/mu/pos/dec/pt/numtof", "MC/pi/pos/dec/pt/numtof",
                                                                 "MC/ka/pos/dec/pt/numtof", "MC/pr/pos/dec/pt/numtof", "MC/de/pos/dec/pt/numtof",
                                                                 "MC/tr/pos/dec/pt/numtof", "MC/he/pos/dec/pt/numtof", "MC/al/pos/dec/pt/numtof",
                                                                 "MC/all/pos/dec/pt/numtof",
                                                                 "MC/el/neg/dec/pt/numtof", "MC/mu/neg/dec/pt/numtof", "MC/pi/neg/dec/pt/numtof",
                                                                 "MC/ka/neg/dec/pt/numtof", "MC/pr/neg/dec/pt/numtof", "MC/de/neg/dec/pt/numtof",
                                                                 "MC/tr/neg/dec/pt/numtof", "MC/he/neg/dec/pt/numtof", "MC/al/neg/dec/pt/numtof",
                                                                 "MC/all/neg/dec/pt/numtof"};
  static constexpr std::string_view hPtDecDen[nHistograms] = {"MC/el/sum/dec/pt/den", "MC/mu/sum/dec/pt/den", "MC/pi/sum/dec/pt/den",
                                                              "MC/ka/sum/dec/pt/den", "MC/pr/sum/dec/pt/den", "MC/de/sum/dec/pt/den",
                                                              "MC/tr/sum/dec/pt/den", "MC/he/sum/dec/pt/den", "MC/al/sum/dec/pt/den",
                                                              "MC/all/sum/dec/pt/den",
                                                              "MC/el/pos/dec/pt/den", "MC/mu/pos/dec/pt/den", "MC/pi/pos/dec/pt/den",
                                                              "MC/ka/pos/dec/pt/den", "MC/pr/pos/dec/pt/den", "MC/de/pos/dec/pt/den",
                                                              "MC/tr/pos/dec/pt/den", "MC/he/pos/dec/pt/den", "MC/al/pos/dec/pt/den",
                                                              "MC/all/pos/dec/pt/den",
                                                              "MC/el/neg/dec/pt/den", "MC/mu/neg/dec/pt/den", "MC/pi/neg/dec/pt/den",
                                                              "MC/ka/neg/dec/pt/den", "MC/pr/neg/dec/pt/den", "MC/de/neg/dec/pt/den",
                                                              "MC/tr/neg/dec/pt/den", "MC/he/neg/dec/pt/den", "MC/al/neg/dec/pt/den",
                                                              "MC/all/neg/dec/pt/den"};
  // Pt for secondaries from material
  static constexpr std::string_view hPtMatNum[nHistograms] = {"MC/el/sum/mat/pt/num", "MC/mu/sum/mat/pt/num", "MC/pi/sum/mat/pt/num",
                                                              "MC/ka/sum/mat/pt/num", "MC/pr/sum/mat/pt/num", "MC/de/sum/mat/pt/num",
                                                              "MC/tr/sum/mat/pt/num", "MC/he/sum/mat/pt/num", "MC/al/sum/mat/pt/num",
                                                              "MC/all/sum/mat/pt/num",
                                                              "MC/el/pos/mat/pt/num", "MC/mu/pos/mat/pt/num", "MC/pi/pos/mat/pt/num",
                                                              "MC/ka/pos/mat/pt/num", "MC/pr/pos/mat/pt/num", "MC/de/pos/mat/pt/num",
                                                              "MC/tr/pos/mat/pt/num", "MC/he/pos/mat/pt/num", "MC/al/pos/mat/pt/num",
                                                              "MC/all/pos/mat/pt/num",
                                                              "MC/el/den/mat/pt/num", "MC/mu/den/mat/pt/num", "MC/pi/den/mat/pt/num",
                                                              "MC/ka/den/mat/pt/num", "MC/pr/den/mat/pt/num", "MC/de/den/mat/pt/num",
                                                              "MC/tr/den/mat/pt/num", "MC/he/den/mat/pt/num", "MC/al/den/mat/pt/num",
                                                              "MC/all/den/mat/pt/num"};
  static constexpr std::string_view hPtMatNumTrk[nHistograms] = {"MC/el/sum/mat/pt/numtrk", "MC/mu/sum/mat/pt/numtrk", "MC/pi/sum/mat/pt/numtrk",
                                                                 "MC/ka/sum/mat/pt/numtrk", "MC/pr/sum/mat/pt/numtrk", "MC/de/sum/mat/pt/numtrk",
                                                                 "MC/tr/sum/mat/pt/numtrk", "MC/he/sum/mat/pt/numtrk", "MC/al/sum/mat/pt/numtrk",
                                                                 "MC/all/sum/mat/pt/numtrk",
                                                                 "MC/el/pos/mat/pt/numtrk", "MC/mu/pos/mat/pt/numtrk", "MC/pi/pos/mat/pt/numtrk",
                                                                 "MC/ka/pos/mat/pt/numtrk", "MC/pr/pos/mat/pt/numtrk", "MC/de/pos/mat/pt/numtrk",
                                                                 "MC/tr/pos/mat/pt/numtrk", "MC/he/pos/mat/pt/numtrk", "MC/al/pos/mat/pt/numtrk",
                                                                 "MC/all/pos/mat/pt/numtrk",
                                                                 "MC/el/neg/mat/pt/numtrk", "MC/mu/neg/mat/pt/numtrk", "MC/pi/neg/mat/pt/numtrk",
                                                                 "MC/ka/neg/mat/pt/numtrk", "MC/pr/neg/mat/pt/numtrk", "MC/de/neg/mat/pt/numtrk",
                                                                 "MC/tr/neg/mat/pt/numtrk", "MC/he/neg/mat/pt/numtrk", "MC/al/neg/mat/pt/numtrk",
                                                                 "MC/all/neg/mat/pt/numtrk"};
  static constexpr std::string_view hPtMatNumTof[nHistograms] = {"MC/el/sum/mat/pt/numtof", "MC/mu/sum/mat/pt/numtof", "MC/pi/sum/mat/pt/numtof",
                                                                 "MC/ka/sum/mat/pt/numtof", "MC/pr/sum/mat/pt/numtof", "MC/de/sum/mat/pt/numtof",
                                                                 "MC/tr/sum/mat/pt/numtof", "MC/he/sum/mat/pt/numtof", "MC/al/sum/mat/pt/numtof",
                                                                 "MC/all/sum/mat/pt/numtof",
                                                                 "MC/el/pos/mat/pt/numtof", "MC/mu/pos/mat/pt/numtof", "MC/pi/pos/mat/pt/numtof",
                                                                 "MC/ka/pos/mat/pt/numtof", "MC/pr/pos/mat/pt/numtof", "MC/de/pos/mat/pt/numtof",
                                                                 "MC/tr/pos/mat/pt/numtof", "MC/he/pos/mat/pt/numtof", "MC/al/pos/mat/pt/numtof",
                                                                 "MC/all/pos/mat/pt/numtof",
                                                                 "MC/el/neg/mat/pt/numtof", "MC/mu/neg/mat/pt/numtof", "MC/pi/neg/mat/pt/numtof",
                                                                 "MC/ka/neg/mat/pt/numtof", "MC/pr/neg/mat/pt/numtof", "MC/de/neg/mat/pt/numtof",
                                                                 "MC/tr/neg/mat/pt/numtof", "MC/he/neg/mat/pt/numtof", "MC/al/neg/mat/pt/numtof",
                                                                 "MC/all/neg/mat/pt/numtof"};
  static constexpr std::string_view hPtMatDen[nHistograms] = {"MC/el/sum/mat/pt/den", "MC/mu/sum/mat/pt/den", "MC/pi/sum/mat/pt/den",
                                                              "MC/ka/sum/mat/pt/den", "MC/pr/sum/mat/pt/den", "MC/de/sum/mat/pt/den",
                                                              "MC/tr/sum/mat/pt/den", "MC/he/sum/mat/pt/den", "MC/al/sum/mat/pt/den",
                                                              "MC/all/sum/mat/pt/den",
                                                              "MC/el/pos/mat/pt/den", "MC/mu/pos/mat/pt/den", "MC/pi/pos/mat/pt/den",
                                                              "MC/ka/pos/mat/pt/den", "MC/pr/pos/mat/pt/den", "MC/de/pos/mat/pt/den",
                                                              "MC/tr/pos/mat/pt/den", "MC/he/pos/mat/pt/den", "MC/al/pos/mat/pt/den",
                                                              "MC/all/pos/mat/pt/den",
                                                              "MC/el/neg/mat/pt/den", "MC/mu/neg/mat/pt/den", "MC/pi/neg/mat/pt/den",
                                                              "MC/ka/neg/mat/pt/den", "MC/pr/neg/mat/pt/den", "MC/de/neg/mat/pt/den",
                                                              "MC/tr/neg/mat/pt/den", "MC/he/neg/mat/pt/den", "MC/al/neg/mat/pt/den",
                                                              "MC/all/neg/mat/pt/den"};
  // P
  static constexpr std::string_view hPNum[nHistograms] = {"MC/el/sum/p/num", "MC/mu/sum/p/num", "MC/pi/sum/p/num",
                                                          "MC/ka/sum/p/num", "MC/pr/sum/p/num", "MC/de/sum/p/num",
                                                          "MC/tr/sum/p/num", "MC/he/sum/p/num", "MC/al/sum/p/num",
                                                          "MC/all/sum/p/num",
                                                          "MC/el/pos/p/num", "MC/mu/pos/p/num", "MC/pi/pos/p/num",
                                                          "MC/ka/pos/p/num", "MC/pr/pos/p/num", "MC/de/pos/p/num",
                                                          "MC/tr/pos/p/num", "MC/he/pos/p/num", "MC/al/pos/p/num",
                                                          "MC/all/pos/p/num",
                                                          "MC/el/neg/p/num", "MC/mu/neg/p/num", "MC/pi/neg/p/num",
                                                          "MC/ka/neg/p/num", "MC/pr/neg/p/num", "MC/de/neg/p/num",
                                                          "MC/tr/neg/p/num", "MC/he/neg/p/num", "MC/al/neg/p/num",
                                                          "MC/all/neg/p/num"};
  static constexpr std::string_view hPNumTrk[nHistograms] = {"MC/el/sum/p/numtrk", "MC/mu/sum/p/numtrk", "MC/pi/sum/p/numtrk",
                                                             "MC/ka/sum/p/numtrk", "MC/pr/sum/p/numtrk", "MC/de/sum/p/numtrk",
                                                             "MC/tr/sum/p/numtrk", "MC/he/sum/p/numtrk", "MC/al/sum/p/numtrk",
                                                             "MC/all/sum/p/numtrk",
                                                             "MC/el/pos/p/numtrk", "MC/mu/pos/p/numtrk", "MC/pi/pos/p/numtrk",
                                                             "MC/ka/pos/p/numtrk", "MC/pr/pos/p/numtrk", "MC/de/pos/p/numtrk",
                                                             "MC/tr/pos/p/numtrk", "MC/he/pos/p/numtrk", "MC/al/pos/p/numtrk",
                                                             "MC/all/pos/p/numtrk",
                                                             "MC/el/neg/p/numtrk", "MC/mu/neg/p/numtrk", "MC/pi/neg/p/numtrk",
                                                             "MC/ka/neg/p/numtrk", "MC/pr/neg/p/numtrk", "MC/de/neg/p/numtrk",
                                                             "MC/tr/neg/p/numtrk", "MC/he/neg/p/numtrk", "MC/al/neg/p/numtrk",
                                                             "MC/all/neg/p/numtrk"};
  static constexpr std::string_view hPNumTof[nHistograms] = {"MC/el/sum/p/numtof", "MC/mu/sum/p/numtof", "MC/pi/sum/p/numtof",
                                                             "MC/ka/sum/p/numtof", "MC/pr/sum/p/numtof", "MC/de/sum/p/numtof",
                                                             "MC/tr/sum/p/numtof", "MC/he/sum/p/numtof", "MC/al/sum/p/numtof",
                                                             "MC/all/sum/p/numtof",
                                                             "MC/el/pos/p/numtof", "MC/mu/pos/p/numtof", "MC/pi/pos/p/numtof",
                                                             "MC/ka/pos/p/numtof", "MC/pr/pos/p/numtof", "MC/de/pos/p/numtof",
                                                             "MC/tr/pos/p/numtof", "MC/he/pos/p/numtof", "MC/al/pos/p/numtof",
                                                             "MC/all/pos/p/numtof",
                                                             "MC/el/sum/p/numtof", "MC/mu/sum/p/numtof", "MC/pi/sum/p/numtof",
                                                             "MC/ka/sum/p/numtof", "MC/pr/sum/p/numtof", "MC/de/sum/p/numtof",
                                                             "MC/tr/sum/p/numtof", "MC/he/sum/p/numtof", "MC/al/sum/p/numtof",
                                                             "MC/all/sum/p/numtof"};
  static constexpr std::string_view hPDen[nHistograms] = {"MC/el/sum/p/den", "MC/mu/sum/p/den", "MC/pi/sum/p/den",
                                                          "MC/ka/sum/p/den", "MC/pr/sum/p/den", "MC/de/sum/p/den",
                                                          "MC/tr/sum/p/den", "MC/he/sum/p/den", "MC/al/sum/p/den",
                                                          "MC/all/sum/p/den",
                                                          "MC/el/pos/p/den", "MC/mu/pos/p/den", "MC/pi/pos/p/den",
                                                          "MC/ka/pos/p/den", "MC/pr/pos/p/den", "MC/de/pos/p/den",
                                                          "MC/tr/pos/p/den", "MC/he/pos/p/den", "MC/al/pos/p/den",
                                                          "MC/all/pos/p/den",
                                                          "MC/el/neg/p/den", "MC/mu/neg/p/den", "MC/pi/neg/p/den",
                                                          "MC/ka/neg/p/den", "MC/pr/neg/p/den", "MC/de/neg/p/den",
                                                          "MC/tr/neg/p/den", "MC/he/neg/p/den", "MC/al/neg/p/den",
                                                          "MC/all/neg/p/den"};
  // Eta
  static constexpr std::string_view hEtaNum[nHistograms] = {"MC/el/sum/eta/num", "MC/mu/sum/eta/num", "MC/pi/sum/eta/num",
                                                            "MC/ka/sum/eta/num", "MC/pr/sum/eta/num", "MC/de/sum/eta/num",
                                                            "MC/tr/sum/eta/num", "MC/he/sum/eta/num", "MC/al/sum/eta/num",
                                                            "MC/all/sum/eta/num",
                                                            "MC/el/pos/eta/num", "MC/mu/pos/eta/num", "MC/pi/pos/eta/num",
                                                            "MC/ka/pos/eta/num", "MC/pr/pos/eta/num", "MC/de/pos/eta/num",
                                                            "MC/tr/pos/eta/num", "MC/he/pos/eta/num", "MC/al/pos/eta/num",
                                                            "MC/all/pos/eta/num",
                                                            "MC/el/neg/eta/num", "MC/mu/neg/eta/num", "MC/pi/neg/eta/num",
                                                            "MC/ka/neg/eta/num", "MC/pr/neg/eta/num", "MC/de/neg/eta/num",
                                                            "MC/tr/neg/eta/num", "MC/he/neg/eta/num", "MC/al/neg/eta/num",
                                                            "MC/all/neg/eta/num"};
  static constexpr std::string_view hEtaNumTrk[nHistograms] = {"MC/el/sum/eta/numtrk", "MC/mu/sum/eta/numtrk", "MC/pi/sum/eta/numtrk",
                                                               "MC/ka/sum/eta/numtrk", "MC/pr/sum/eta/numtrk", "MC/de/sum/eta/numtrk",
                                                               "MC/tr/sum/eta/numtrk", "MC/he/sum/eta/numtrk", "MC/al/sum/eta/numtrk",
                                                               "MC/all/sum/eta/numtrk",
                                                               "MC/el/pos/eta/numtrk", "MC/mu/pos/eta/numtrk", "MC/pi/pos/eta/numtrk",
                                                               "MC/ka/pos/eta/numtrk", "MC/pr/pos/eta/numtrk", "MC/de/pos/eta/numtrk",
                                                               "MC/tr/pos/eta/numtrk", "MC/he/pos/eta/numtrk", "MC/al/pos/eta/numtrk",
                                                               "MC/all/pos/eta/numtrk",
                                                               "MC/el/neg/eta/numtrk", "MC/mu/neg/eta/numtrk", "MC/pi/neg/eta/numtrk",
                                                               "MC/ka/neg/eta/numtrk", "MC/pr/neg/eta/numtrk", "MC/de/neg/eta/numtrk",
                                                               "MC/tr/neg/eta/numtrk", "MC/he/neg/eta/numtrk", "MC/al/neg/eta/numtrk",
                                                               "MC/all/neg/eta/numtrk"};
  static constexpr std::string_view hEtaNumTof[nHistograms] = {"MC/el/sum/eta/numtof", "MC/mu/sum/eta/numtof", "MC/pi/sum/eta/numtof",
                                                               "MC/ka/sum/eta/numtof", "MC/pr/sum/eta/numtof", "MC/de/sum/eta/numtof",
                                                               "MC/tr/sum/eta/numtof", "MC/he/sum/eta/numtof", "MC/al/sum/eta/numtof",
                                                               "MC/all/sum/eta/numtof",
                                                               "MC/el/pos/eta/numtof", "MC/mu/pos/eta/numtof", "MC/pi/pos/eta/numtof",
                                                               "MC/ka/pos/eta/numtof", "MC/pr/pos/eta/numtof", "MC/de/pos/eta/numtof",
                                                               "MC/tr/pos/eta/numtof", "MC/he/pos/eta/numtof", "MC/al/pos/eta/numtof",
                                                               "MC/all/pos/eta/numtof",
                                                               "MC/el/sum/eta/numtof", "MC/mu/sum/eta/numtof", "MC/pi/sum/eta/numtof",
                                                               "MC/ka/sum/eta/numtof", "MC/pr/sum/eta/numtof", "MC/de/sum/eta/numtof",
                                                               "MC/tr/sum/eta/numtof", "MC/he/sum/eta/numtof", "MC/al/sum/eta/numtof",
                                                               "MC/all/sum/eta/numtof"};
  static constexpr std::string_view hEtaDen[nHistograms] = {"MC/el/sum/eta/den", "MC/mu/sum/eta/den", "MC/pi/sum/eta/den",
                                                            "MC/ka/sum/eta/den", "MC/pr/sum/eta/den", "MC/de/sum/eta/den",
                                                            "MC/tr/sum/eta/den", "MC/he/sum/eta/den", "MC/al/sum/eta/den",
                                                            "MC/all/sum/eta/den",
                                                            "MC/el/pos/eta/den", "MC/mu/pos/eta/den", "MC/pi/pos/eta/den",
                                                            "MC/ka/pos/eta/den", "MC/pr/pos/eta/den", "MC/de/pos/eta/den",
                                                            "MC/tr/pos/eta/den", "MC/he/pos/eta/den", "MC/al/pos/eta/den",
                                                            "MC/all/pos/eta/den",
                                                            "MC/el/neg/eta/den", "MC/mu/neg/eta/den", "MC/pi/neg/eta/den",
                                                            "MC/ka/neg/eta/den", "MC/pr/neg/eta/den", "MC/de/neg/eta/den",
                                                            "MC/tr/neg/eta/den", "MC/he/neg/eta/den", "MC/al/neg/eta/den",
                                                            "MC/all/neg/eta/den"};
  // Y
  static constexpr std::string_view hYNum[nHistograms] = {"MC/el/sum/y/num", "MC/mu/sum/y/num", "MC/pi/sum/y/num",
                                                          "MC/ka/sum/y/num", "MC/pr/sum/y/num", "MC/de/sum/y/num",
                                                          "MC/tr/sum/y/num", "MC/he/sum/y/num", "MC/al/sum/y/num",
                                                          "MC/all/sum/y/num",
                                                          "MC/el/pos/y/num", "MC/mu/pos/y/num", "MC/pi/pos/y/num",
                                                          "MC/ka/pos/y/num", "MC/pr/pos/y/num", "MC/de/pos/y/num",
                                                          "MC/tr/pos/y/num", "MC/he/pos/y/num", "MC/al/pos/y/num",
                                                          "MC/all/pos/y/num",
                                                          "MC/el/neg/y/num", "MC/mu/neg/y/num", "MC/pi/neg/y/num",
                                                          "MC/ka/neg/y/num", "MC/pr/neg/y/num", "MC/de/neg/y/num",
                                                          "MC/tr/neg/y/num", "MC/he/neg/y/num", "MC/al/neg/y/num",
                                                          "MC/all/neg/y/num"};
  static constexpr std::string_view hYNumTof[nHistograms] = {"MC/el/sum/y/numtof", "MC/mu/sum/y/numtof", "MC/pi/sum/y/numtof",
                                                             "MC/ka/sum/y/numtof", "MC/pr/sum/y/numtof", "MC/de/sum/y/numtof",
                                                             "MC/tr/sum/y/numtof", "MC/he/sum/y/numtof", "MC/al/sum/y/numtof",
                                                             "MC/all/sum/y/numtof",
                                                             "MC/el/pos/y/numtof", "MC/mu/pos/y/numtof", "MC/pi/pos/y/numtof",
                                                             "MC/ka/pos/y/numtof", "MC/pr/pos/y/numtof", "MC/de/pos/y/numtof",
                                                             "MC/tr/pos/y/numtof", "MC/he/pos/y/numtof", "MC/al/pos/y/numtof",
                                                             "MC/all/pos/y/numtof",
                                                             "MC/el/sum/y/numtof", "MC/mu/sum/y/numtof", "MC/pi/sum/y/numtof",
                                                             "MC/ka/sum/y/numtof", "MC/pr/sum/y/numtof", "MC/de/sum/y/numtof",
                                                             "MC/tr/sum/y/numtof", "MC/he/sum/y/numtof", "MC/al/sum/y/numtof",
                                                             "MC/all/sum/y/numtof"};
  static constexpr std::string_view hYDen[nHistograms] = {"MC/el/sum/y/den", "MC/mu/sum/y/den", "MC/pi/sum/y/den",
                                                          "MC/ka/sum/y/den", "MC/pr/sum/y/den", "MC/de/sum/y/den",
                                                          "MC/tr/sum/y/den", "MC/he/sum/y/den", "MC/al/sum/y/den",
                                                          "MC/all/sum/y/den",
                                                          "MC/el/pos/y/den", "MC/mu/pos/y/den", "MC/pi/pos/y/den",
                                                          "MC/ka/pos/y/den", "MC/pr/pos/y/den", "MC/de/pos/y/den",
                                                          "MC/tr/pos/y/den", "MC/he/pos/y/den", "MC/al/pos/y/den",
                                                          "MC/all/pos/y/den",
                                                          "MC/el/neg/y/den", "MC/mu/neg/y/den", "MC/pi/neg/y/den",
                                                          "MC/ka/neg/y/den", "MC/pr/neg/y/den", "MC/de/neg/y/den",
                                                          "MC/tr/neg/y/den", "MC/he/neg/y/den", "MC/al/neg/y/den",
                                                          "MC/all/neg/y/den"};
  // Phi
  static constexpr std::string_view hPhiNum[nHistograms] = {"MC/el/sum/phi/num", "MC/mu/sum/phi/num", "MC/pi/sum/phi/num",
                                                            "MC/ka/sum/phi/num", "MC/pr/sum/phi/num", "MC/de/sum/phi/num",
                                                            "MC/tr/sum/phi/num", "MC/he/sum/phi/num", "MC/al/sum/phi/num",
                                                            "MC/all/sum/phi/num",
                                                            "MC/el/pos/phi/num", "MC/mu/pos/phi/num", "MC/pi/pos/phi/num",
                                                            "MC/ka/pos/phi/num", "MC/pr/pos/phi/num", "MC/de/pos/phi/num",
                                                            "MC/tr/pos/phi/num", "MC/he/pos/phi/num", "MC/al/pos/phi/num",
                                                            "MC/all/pos/phi/num",
                                                            "MC/el/neg/phi/num", "MC/mu/neg/phi/num", "MC/pi/neg/phi/num",
                                                            "MC/ka/neg/phi/num", "MC/pr/neg/phi/num", "MC/de/neg/phi/num",
                                                            "MC/tr/neg/phi/num", "MC/he/neg/phi/num", "MC/al/neg/phi/num",
                                                            "MC/all/neg/phi/num"};
  static constexpr std::string_view hPhiNumTrk[nHistograms] = {"MC/el/sum/phi/numtrk", "MC/mu/sum/phi/numtrk", "MC/pi/sum/phi/numtrk",
                                                               "MC/ka/sum/phi/numtrk", "MC/pr/sum/phi/numtrk", "MC/de/sum/phi/numtrk",
                                                               "MC/tr/sum/phi/numtrk", "MC/he/sum/phi/numtrk", "MC/al/sum/phi/numtrk",
                                                               "MC/all/sum/phi/numtrk",
                                                               "MC/el/pos/phi/numtrk", "MC/mu/pos/phi/numtrk", "MC/pi/pos/phi/numtrk",
                                                               "MC/ka/pos/phi/numtrk", "MC/pr/pos/phi/numtrk", "MC/de/pos/phi/numtrk",
                                                               "MC/tr/pos/phi/numtrk", "MC/he/pos/phi/numtrk", "MC/al/pos/phi/numtrk",
                                                               "MC/all/pos/phi/numtrk",
                                                               "MC/el/neg/phi/numtrk", "MC/mu/neg/phi/numtrk", "MC/pi/neg/phi/numtrk",
                                                               "MC/ka/neg/phi/numtrk", "MC/pr/neg/phi/numtrk", "MC/de/neg/phi/numtrk",
                                                               "MC/tr/neg/phi/numtrk", "MC/he/neg/phi/numtrk", "MC/al/neg/phi/numtrk",
                                                               "MC/all/neg/phi/numtrk"};
  static constexpr std::string_view hPhiNumTof[nHistograms] = {"MC/el/sum/phi/numtof", "MC/mu/sum/phi/numtof", "MC/pi/sum/phi/numtof",
                                                               "MC/ka/sum/phi/numtof", "MC/pr/sum/phi/numtof", "MC/de/sum/phi/numtof",
                                                               "MC/tr/sum/phi/numtof", "MC/he/sum/phi/numtof", "MC/al/sum/phi/numtof",
                                                               "MC/all/sum/phi/numtof",
                                                               "MC/el/pos/phi/numtof", "MC/mu/pos/phi/numtof", "MC/pi/pos/phi/numtof",
                                                               "MC/ka/pos/phi/numtof", "MC/pr/pos/phi/numtof", "MC/de/pos/phi/numtof",
                                                               "MC/tr/pos/phi/numtof", "MC/he/pos/phi/numtof", "MC/al/pos/phi/numtof",
                                                               "MC/all/pos/phi/numtof",
                                                               "MC/el/sum/phi/numtof", "MC/mu/sum/phi/numtof", "MC/pi/sum/phi/numtof",
                                                               "MC/ka/sum/phi/numtof", "MC/pr/sum/phi/numtof", "MC/de/sum/phi/numtof",
                                                               "MC/tr/sum/phi/numtof", "MC/he/sum/phi/numtof", "MC/al/sum/phi/numtof",
                                                               "MC/all/sum/phi/numtof"};
  static constexpr std::string_view hPhiDen[nHistograms] = {"MC/el/sum/phi/den", "MC/mu/sum/phi/den", "MC/pi/sum/phi/den",
                                                            "MC/ka/sum/phi/den", "MC/pr/sum/phi/den", "MC/de/sum/phi/den",
                                                            "MC/tr/sum/phi/den", "MC/he/sum/phi/den", "MC/al/sum/phi/den",
                                                            "MC/all/sum/phi/den",
                                                            "MC/el/pos/phi/den", "MC/mu/pos/phi/den", "MC/pi/pos/phi/den",
                                                            "MC/ka/pos/phi/den", "MC/pr/pos/phi/den", "MC/de/pos/phi/den",
                                                            "MC/tr/pos/phi/den", "MC/he/pos/phi/den", "MC/al/pos/phi/den",
                                                            "MC/all/pos/phi/den",
                                                            "MC/el/neg/phi/den", "MC/mu/neg/phi/den", "MC/pi/neg/phi/den",
                                                            "MC/ka/neg/phi/den", "MC/pr/neg/phi/den", "MC/de/neg/phi/den",
                                                            "MC/tr/neg/phi/den", "MC/he/neg/phi/den", "MC/al/neg/phi/den",
                                                            "MC/all/neg/phi/den"};
  // Pt-Eta
  static constexpr std::string_view hPtEtaNum[nHistograms] = {"MC/el/sum/pteta/num", "MC/mu/sum/pteta/num", "MC/pi/sum/pteta/num",
                                                              "MC/ka/sum/pteta/num", "MC/pr/sum/pteta/num", "MC/de/sum/pteta/num",
                                                              "MC/tr/sum/pteta/num", "MC/he/sum/pteta/num", "MC/al/sum/pteta/num",
                                                              "MC/all/sum/pteta/num",
                                                              "MC/el/pos/pteta/num", "MC/mu/pos/pteta/num", "MC/pi/pos/pteta/num",
                                                              "MC/ka/pos/pteta/num", "MC/pr/pos/pteta/num", "MC/de/pos/pteta/num",
                                                              "MC/tr/pos/pteta/num", "MC/he/pos/pteta/num", "MC/al/pos/pteta/num",
                                                              "MC/all/pos/pteta/num",
                                                              "MC/el/neg/pteta/num", "MC/mu/neg/pteta/num", "MC/pi/neg/pteta/num",
                                                              "MC/ka/neg/pteta/num", "MC/pr/neg/pteta/num", "MC/de/neg/pteta/num",
                                                              "MC/tr/neg/pteta/num", "MC/he/neg/pteta/num", "MC/al/neg/pteta/num",
                                                              "MC/all/neg/pteta/num"};
  static constexpr std::string_view hPtEtaNumTrk[nHistograms] = {"MC/el/sum/pteta/numtrk", "MC/mu/sum/pteta/numtrk", "MC/pi/sum/pteta/numtrk",
                                                                 "MC/ka/sum/pteta/numtrk", "MC/pr/sum/pteta/numtrk", "MC/de/sum/pteta/numtrk",
                                                                 "MC/tr/sum/pteta/numtrk", "MC/he/sum/pteta/numtrk", "MC/al/sum/pteta/numtrk",
                                                                 "MC/all/sum/pteta/numtrk",
                                                                 "MC/el/pos/pteta/numtrk", "MC/mu/pos/pteta/numtrk", "MC/pi/pos/pteta/numtrk",
                                                                 "MC/ka/pos/pteta/numtrk", "MC/pr/pos/pteta/numtrk", "MC/de/pos/pteta/numtrk",
                                                                 "MC/tr/pos/pteta/numtrk", "MC/he/pos/pteta/numtrk", "MC/al/pos/pteta/numtrk",
                                                                 "MC/all/pos/pteta/numtrk",
                                                                 "MC/el/neg/pteta/numtrk", "MC/mu/neg/pteta/numtrk", "MC/pi/neg/pteta/numtrk",
                                                                 "MC/ka/neg/pteta/numtrk", "MC/pr/neg/pteta/numtrk", "MC/de/neg/pteta/numtrk",
                                                                 "MC/tr/neg/pteta/numtrk", "MC/he/neg/pteta/numtrk", "MC/al/neg/pteta/numtrk",
                                                                 "MC/all/neg/pteta/numtrk"};
  static constexpr std::string_view hPtEtaNumTof[nHistograms] = {"MC/el/sum/pteta/numtof", "MC/mu/sum/pteta/numtof", "MC/pi/sum/pteta/numtof",
                                                                 "MC/ka/sum/pteta/numtof", "MC/pr/sum/pteta/numtof", "MC/de/sum/pteta/numtof",
                                                                 "MC/tr/sum/pteta/numtof", "MC/he/sum/pteta/numtof", "MC/al/sum/pteta/numtof",
                                                                 "MC/all/sum/pteta/numtof",
                                                                 "MC/el/pos/pteta/numtof", "MC/mu/pos/pteta/numtof", "MC/pi/pos/pteta/numtof",
                                                                 "MC/ka/pos/pteta/numtof", "MC/pr/pos/pteta/numtof", "MC/de/pos/pteta/numtof",
                                                                 "MC/tr/pos/pteta/numtof", "MC/he/pos/pteta/numtof", "MC/al/pos/pteta/numtof",
                                                                 "MC/all/pos/pteta/numtof",
                                                                 "MC/el/sum/pteta/numtof", "MC/mu/sum/pteta/numtof", "MC/pi/sum/pteta/numtof",
                                                                 "MC/ka/sum/pteta/numtof", "MC/pr/sum/pteta/numtof", "MC/de/sum/pteta/numtof",
                                                                 "MC/tr/sum/pteta/numtof", "MC/he/sum/pteta/numtof", "MC/al/sum/pteta/numtof",
                                                                 "MC/all/sum/pteta/numtof"};
  static constexpr std::string_view hPtEtaDen[nHistograms] = {"MC/el/sum/pteta/den", "MC/mu/sum/pteta/den", "MC/pi/sum/pteta/den",
                                                              "MC/ka/sum/pteta/den", "MC/pr/sum/pteta/den", "MC/de/sum/pteta/den",
                                                              "MC/tr/sum/pteta/den", "MC/he/sum/pteta/den", "MC/al/sum/pteta/den",
                                                              "MC/all/sum/pteta/den",
                                                              "MC/el/pos/pteta/den", "MC/mu/pos/pteta/den", "MC/pi/pos/pteta/den",
                                                              "MC/ka/pos/pteta/den", "MC/pr/pos/pteta/den", "MC/de/pos/pteta/den",
                                                              "MC/tr/pos/pteta/den", "MC/he/pos/pteta/den", "MC/al/pos/pteta/den",
                                                              "MC/all/pos/pteta/den",
                                                              "MC/el/neg/pteta/den", "MC/mu/neg/pteta/den", "MC/pi/neg/pteta/den",
                                                              "MC/ka/neg/pteta/den", "MC/pr/neg/pteta/den", "MC/de/neg/pteta/den",
                                                              "MC/tr/neg/pteta/den", "MC/he/neg/pteta/den", "MC/al/neg/pteta/den",
                                                              "MC/all/neg/pteta/den"};

  template <o2::track::PID::ID id>
  void makeMCHistograms(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }
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

    const TString tagEta = Form("%s #it{p}_{T} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                                partName,
                                ptMin.value, ptMax.value,
                                yMin.value, yMax.value,
                                phiMin.value, phiMax.value);

    const TString tagY = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
                              partName,
                              ptMin.value, ptMax.value,
                              etaMin.value, etaMax.value,
                              phiMin.value, phiMax.value);

    const TString tagPhi = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
                                partName,
                                ptMin.value, ptMax.value,
                                etaMin.value, etaMax.value,
                                yMin.value, yMax.value);

    const TString tagPtEta = Form("%s #it{#varphi} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
                                  partName,
                                  phiMin.value, phiMax.value,
                                  yMin.value, yMax.value);

    auto makeHistogramsPerCharge = [&](const int chargeIndex) {
      const int histogramIndex = id + chargeIndex * nSpecies;
      histos.add(hPtNum[histogramIndex].data(), "Numerator " + tagPt, kTH1D, {axisPt});
      histos.add(hPtNumTrk[histogramIndex].data(), "Numerator Track " + tagPt, kTH1D, {axisPt});
      histos.add(hPtNumTof[histogramIndex].data(), "Numerator TOF " + tagPt, kTH1D, {axisPt});
      histos.add(hPtDen[histogramIndex].data(), "Denominator " + tagPt, kTH1D, {axisPt});

      histos.add(hPtPrmNum[histogramIndex].data(), "Numerator " + tagPt + " Primaries", kTH1D, {axisPt});
      histos.add(hPtPrmNumTrk[histogramIndex].data(), "Numerator Track " + tagPt + " Primaries", kTH1D, {axisPt});
      histos.add(hPtPrmNumTof[histogramIndex].data(), "Numerator TOF " + tagPt + " Primaries", kTH1D, {axisPt});
      histos.add(hPtPrmDen[histogramIndex].data(), "Denominator " + tagPt + " Primaries", kTH1D, {axisPt});

      histos.add(hPtDecNum[histogramIndex].data(), "Numerator " + tagPt + " Sec. from decays", kTH1D, {axisPt});
      histos.add(hPtDecNumTrk[histogramIndex].data(), "Numerator Track " + tagPt + " Sec. from decays", kTH1D, {axisPt});
      histos.add(hPtDecNumTof[histogramIndex].data(), "Numerator TOF " + tagPt + " Sec. from decays", kTH1D, {axisPt});
      histos.add(hPtDecDen[histogramIndex].data(), "Denominator " + tagPt + " Sec. from decays", kTH1D, {axisPt});

      histos.add(hPtMatNum[histogramIndex].data(), "Numerator " + tagPt + " Sec. from material", kTH1D, {axisPt});
      histos.add(hPtMatNumTrk[histogramIndex].data(), "Numerator Track " + tagPt + " Sec. from material", kTH1D, {axisPt});
      histos.add(hPtMatNumTof[histogramIndex].data(), "Numerator TOF " + tagPt + " Sec. from material", kTH1D, {axisPt});
      histos.add(hPtMatDen[histogramIndex].data(), "Denominator " + tagPt + " Sec. from material", kTH1D, {axisPt});

      histos.add(hPNum[histogramIndex].data(), "Numerator " + tagPt, kTH1D, {axisP});
      histos.add(hPNumTrk[histogramIndex].data(), "Numerator Track " + tagPt, kTH1D, {axisP});
      histos.add(hPNumTof[histogramIndex].data(), "Numerator TOF " + tagPt, kTH1D, {axisP});
      histos.add(hPDen[histogramIndex].data(), "Denominator " + tagPt, kTH1D, {axisP});

      histos.add(hEtaNum[histogramIndex].data(), "Numerator " + tagEta, kTH1D, {axisEta});
      histos.add(hEtaNumTrk[histogramIndex].data(), "Numerator Track " + tagEta, kTH1D, {axisEta});
      histos.add(hEtaNumTof[histogramIndex].data(), "Numerator TOF " + tagEta, kTH1D, {axisEta});
      histos.add(hEtaDen[histogramIndex].data(), "Denominator " + tagEta, kTH1D, {axisEta});

      histos.add(hYNum[histogramIndex].data(), "Numerator " + tagY, kTH1D, {axisY});
      histos.add(hYNumTof[histogramIndex].data(), "Numerator TOF " + tagY, kTH1D, {axisY});
      histos.add(hYDen[histogramIndex].data(), "Denominator " + tagY, kTH1D, {axisY});

      histos.add(hPhiNum[histogramIndex].data(), "Numerator " + tagPhi, kTH1D, {axisPhi});
      histos.add(hPhiNumTrk[histogramIndex].data(), "Numerator Track " + tagPhi, kTH1D, {axisPhi});
      histos.add(hPhiNumTof[histogramIndex].data(), "Numerator TOF " + tagPhi, kTH1D, {axisPhi});
      histos.add(hPhiDen[histogramIndex].data(), "Denominator " + tagPhi, kTH1D, {axisPhi});

      histos.add(hPtEtaNum[histogramIndex].data(), "Numerator " + tagPtEta, kTH2D, {axisPt, axisEta});
      histos.add(hPtEtaDen[histogramIndex].data(), "Denominator " + tagPtEta, kTH2D, {axisPt, axisEta});
    };

    if (doSumPDG) { // Sum
      makeHistogramsPerCharge(0);
    }
    if (doPositivePDG) { // Positive
      makeHistogramsPerCharge(1);
    }
    if (doNegativePDG) { // Negative
      makeHistogramsPerCharge(2);
    }

    if (makeEff && doSumPDG) {
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
    makeMCHistograms<o2::track::PID::Electron>(doEl);
    makeMCHistograms<o2::track::PID::Muon>(doMu);
    makeMCHistograms<o2::track::PID::Pion>(doPi);
    makeMCHistograms<o2::track::PID::Kaon>(doKa);
    makeMCHistograms<o2::track::PID::Proton>(doPr);
    makeMCHistograms<o2::track::PID::Deuteron>(doDe);
    makeMCHistograms<o2::track::PID::Triton>(doTr);
    makeMCHistograms<o2::track::PID::Helium3>(doHe);
    makeMCHistograms<o2::track::PID::Alpha>(doAl);
    makeMCHistograms<o2::track::PID::NIDs>(doUnId);
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

  template <int pdgSign, o2::track::PID::ID id, typename particleType>
  bool isPdgSelected(particleType mcParticle)
  {
    static_assert(pdgSign == 0 || pdgSign == 1 || pdgSign == 2);
    // Selecting PDG code
    if constexpr (PDGs[id] == 0) {
      if constexpr (pdgSign == 0) {
        return true;
      } else if constexpr (pdgSign == 1) {
        return mcParticle.pdgCode() > 0;
      }
      if constexpr (pdgSign == 2) {
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
    } else if constexpr (pdgSign == 2) {
      if (mcParticle.pdgCode() != -PDGs[id]) {
        return false;
      }
    }
    return true;
  }

  template <int pdgSign, o2::track::PID::ID id, typename trackType>
  void fillMCTrackHistograms(const trackType& track)
  {
    static_assert(pdgSign == 0 || pdgSign == 1 || pdgSign == 2);
    if constexpr (pdgSign == 0) {
      if (!doSumPDG) {
        return;
      }
    } else if constexpr (pdgSign == 1) {
      if (!doPositivePDG) {
        return;
      }
    } else if constexpr (pdgSign == 2) {
      if (!doNegativePDG) {
        return;
      }
    }
    constexpr int histogramIndex = id + pdgSign * nSpecies;
    LOG(debug) << "Filling track histograms for id " << static_cast<int>(id);
    const auto mcParticle = track.mcParticle();

    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("MC/trackSelection"), 11 + id);

    histos.fill(HIST(hPNum[histogramIndex]), mcParticle.p());
    histos.fill(HIST(hPtNum[histogramIndex]), mcParticle.pt());
    histos.fill(HIST(hEtaNum[histogramIndex]), mcParticle.eta());
    histos.fill(HIST(hYNum[histogramIndex]), mcParticle.y());
    histos.fill(HIST(hPhiNum[histogramIndex]), mcParticle.phi());
    histos.fill(HIST(hPtEtaNum[histogramIndex]), mcParticle.pt(), mcParticle.eta());

    histos.fill(HIST(hPNumTrk[histogramIndex]), track.p());
    histos.fill(HIST(hPtNumTrk[histogramIndex]), track.pt());
    histos.fill(HIST(hEtaNumTrk[histogramIndex]), track.eta());
    histos.fill(HIST(hPhiNumTrk[histogramIndex]), track.phi());

    if (mcParticle.isPhysicalPrimary()) {
      histos.fill(HIST(hPtPrmNum[histogramIndex]), mcParticle.pt());
      histos.fill(HIST(hPtPrmNumTrk[histogramIndex]), track.pt());
      if (track.hasTOF()) {
        histos.fill(HIST(hPtPrmNumTof[histogramIndex]), mcParticle.pt());
      }
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        histos.fill(HIST(hPtDecNum[histogramIndex]), mcParticle.pt());
        histos.fill(HIST(hPtDecNumTrk[histogramIndex]), track.pt());
        if (track.hasTOF()) {
          histos.fill(HIST(hPtDecNumTof[histogramIndex]), mcParticle.pt());
        }
      } else { // Material
        histos.fill(HIST(hPtMatNum[histogramIndex]), mcParticle.pt());
        histos.fill(HIST(hPtMatNumTrk[histogramIndex]), track.pt());
        if (track.hasTOF()) {
          histos.fill(HIST(hPtMatNumTof[histogramIndex]), mcParticle.pt());
        }
      }
    }
    if (!track.hasTOF()) {
      return;
    }
    histos.fill(HIST(hPNumTof[histogramIndex]), mcParticle.p());
    histos.fill(HIST(hPtNumTof[histogramIndex]), mcParticle.pt());
    histos.fill(HIST(hEtaNumTof[histogramIndex]), mcParticle.eta());
    histos.fill(HIST(hYNumTof[histogramIndex]), mcParticle.y());
    histos.fill(HIST(hPhiNumTof[histogramIndex]), mcParticle.phi());
  }

  template <int pdgSign, o2::track::PID::ID id, typename particleType>
  void fillMCParticleHistograms(const particleType& mcParticle)
  {
    static_assert(pdgSign == 0 || pdgSign == 1 || pdgSign == 2);
    if constexpr (pdgSign == 0) {
      if (!doSumPDG) {
        return;
      }
    } else if constexpr (pdgSign == 1) {
      if (!doPositivePDG) {
        return;
      }
    } else if constexpr (pdgSign == 2) {
      if (!doNegativePDG) {
        return;
      }
    }
    constexpr int histogramIndex = id + pdgSign * nSpecies;
    LOG(debug) << "Filling particle histograms for id " << static_cast<int>(id);
    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }
    histos.fill(HIST("MC/particleSelection"), 7 + id);

    histos.fill(HIST(hPDen[histogramIndex]), mcParticle.p());
    histos.fill(HIST(hPtDen[histogramIndex]), mcParticle.pt());

    if (mcParticle.isPhysicalPrimary()) {
      histos.fill(HIST(hPtPrmDen[histogramIndex]), mcParticle.pt());
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        histos.fill(HIST(hPtDecDen[histogramIndex]), mcParticle.pt());
      } else { // Material
        histos.fill(HIST(hPtMatDen[histogramIndex]), mcParticle.pt());
      }
    }

    histos.fill(HIST(hEtaDen[histogramIndex]), mcParticle.eta());
    histos.fill(HIST(hYDen[histogramIndex]), mcParticle.y());
    histos.fill(HIST(hPhiDen[histogramIndex]), mcParticle.phi());
    histos.fill(HIST(hPtEtaDen[histogramIndex]), mcParticle.pt(), mcParticle.eta());
  }

  template <o2::track::PID::ID id>
  void fillMCEfficiency()
  {
    if (!makeEff) {
      return;
    }
    if (!doSumPDG) {
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
        fillMCTrackHistograms<0, o2::track::PID::Electron>(track);
        fillMCTrackHistograms<1, o2::track::PID::Electron>(track);
        fillMCTrackHistograms<2, o2::track::PID::Electron>(track);
      }
      if (doMu) {
        fillMCTrackHistograms<0, o2::track::PID::Muon>(track);
        fillMCTrackHistograms<1, o2::track::PID::Muon>(track);
        fillMCTrackHistograms<2, o2::track::PID::Muon>(track);
      }
      if (doPi) {
        fillMCTrackHistograms<0, o2::track::PID::Pion>(track);
        fillMCTrackHistograms<1, o2::track::PID::Pion>(track);
        fillMCTrackHistograms<2, o2::track::PID::Pion>(track);
      }
      if (doKa) {
        fillMCTrackHistograms<0, o2::track::PID::Kaon>(track);
        fillMCTrackHistograms<1, o2::track::PID::Kaon>(track);
        fillMCTrackHistograms<2, o2::track::PID::Kaon>(track);
      }
      if (doPr) {
        fillMCTrackHistograms<0, o2::track::PID::Proton>(track);
        fillMCTrackHistograms<1, o2::track::PID::Proton>(track);
        fillMCTrackHistograms<2, o2::track::PID::Proton>(track);
      }
      if (doDe) {
        fillMCTrackHistograms<0, o2::track::PID::Deuteron>(track);
        fillMCTrackHistograms<1, o2::track::PID::Deuteron>(track);
        fillMCTrackHistograms<2, o2::track::PID::Deuteron>(track);
      }
      if (doTr) {
        fillMCTrackHistograms<0, o2::track::PID::Triton>(track);
        fillMCTrackHistograms<1, o2::track::PID::Triton>(track);
        fillMCTrackHistograms<2, o2::track::PID::Triton>(track);
      }
      if (doHe) {
        fillMCTrackHistograms<0, o2::track::PID::Helium3>(track);
        fillMCTrackHistograms<1, o2::track::PID::Helium3>(track);
        fillMCTrackHistograms<2, o2::track::PID::Helium3>(track);
      }
      if (doAl) {
        fillMCTrackHistograms<0, o2::track::PID::Alpha>(track);
        fillMCTrackHistograms<1, o2::track::PID::Alpha>(track);
        fillMCTrackHistograms<2, o2::track::PID::Alpha>(track);
      }
      if (doUnId) {
        fillMCTrackHistograms<0, o2::track::PID::NIDs>(track);
        fillMCTrackHistograms<1, o2::track::PID::NIDs>(track);
        fillMCTrackHistograms<2, o2::track::PID::NIDs>(track);
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
        fillMCParticleHistograms<0, o2::track::PID::Electron>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Electron>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Electron>(mcParticle);
      }
      if (doMu) {
        fillMCParticleHistograms<0, o2::track::PID::Muon>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Muon>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Muon>(mcParticle);
      }
      if (doPi) {
        fillMCParticleHistograms<0, o2::track::PID::Pion>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Pion>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Pion>(mcParticle);
      }
      if (doKa) {
        fillMCParticleHistograms<0, o2::track::PID::Kaon>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Kaon>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Kaon>(mcParticle);
      }
      if (doPr) {
        fillMCParticleHistograms<0, o2::track::PID::Proton>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Proton>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Proton>(mcParticle);
      }
      if (doDe) {
        fillMCParticleHistograms<0, o2::track::PID::Deuteron>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Deuteron>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Deuteron>(mcParticle);
      }
      if (doTr) {
        fillMCParticleHistograms<0, o2::track::PID::Triton>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Triton>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Triton>(mcParticle);
      }
      if (doHe) {
        fillMCParticleHistograms<0, o2::track::PID::Helium3>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Helium3>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Helium3>(mcParticle);
      }
      if (doAl) {
        fillMCParticleHistograms<0, o2::track::PID::Alpha>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Alpha>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::Alpha>(mcParticle);
      }
      if (doUnId) {
        fillMCParticleHistograms<0, o2::track::PID::NIDs>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::NIDs>(mcParticle);
        fillMCParticleHistograms<2, o2::track::PID::NIDs>(mcParticle);
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
  return WorkflowSpec{adaptAnalysisTask<QaEfficiency>(cfgc)};
}
