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
#include "Framework/StaticFor.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

// ROOT includes
#include "TPDGCode.h"
#include "TEfficiency.h"
#include "THashList.h"

using namespace o2::framework;

struct QaEfficiency {
  // Particle information
  static constexpr int nSpecies = o2::track::PID::NIDs + 1; // One per PDG + 1 for unidentified
  static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha", "All"};
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040, 0};
  // Track/particle selection
  Configurable<bool> noFakesHits{"noFakesHits", false, "Flag to reject tracks that have fake hits"};
  // Charge selection
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
  // Track only selection, options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  // Event selection
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the generated vertez in Z (cm)"};
  Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the generated vertez in Z (cm)"};
  // Histogram configuration
  ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  Configurable<int> logPt{"log-pt", 0, "Flag to use a logarithmic pT axis"};
  ConfigurableAxis etaBins{"etaBins", {200, -3.f, 3.f}, "Eta binning"};
  ConfigurableAxis phiBins{"phiBins", {200, 0.f, 6.284f}, "Phi binning"};
  ConfigurableAxis yBins{"yBins", {200, -0.5f, 0.5f}, "Y binning"};
  // Task configuration
  Configurable<bool> makeEff{"make-eff", false, "Flag to produce the efficiency with TEfficiency"};
  Configurable<bool> doPtEta{"doPtEta", false, "Flag to produce the efficiency vs pT and Eta"};
  Configurable<int> applyEvSel{"applyEvSel", 0, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};

  OutputObj<THashList> listEfficiencyMC{"EfficiencyMC"};
  OutputObj<THashList> listEfficiencyData{"EfficiencyData"};
  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  static constexpr int nHistograms = nSpecies * 2;

  // Pt
  static constexpr std::string_view hPtIts[nHistograms] = {"MC/el/pos/pt/its", "MC/mu/pos/pt/its", "MC/pi/pos/pt/its",
                                                           "MC/ka/pos/pt/its", "MC/pr/pos/pt/its", "MC/de/pos/pt/its",
                                                           "MC/tr/pos/pt/its", "MC/he/pos/pt/its", "MC/al/pos/pt/its",
                                                           "MC/all/pos/pt/its",
                                                           "MC/el/neg/pt/its", "MC/mu/neg/pt/its", "MC/pi/neg/pt/its",
                                                           "MC/ka/neg/pt/its", "MC/pr/neg/pt/its", "MC/de/neg/pt/its",
                                                           "MC/tr/neg/pt/its", "MC/he/neg/pt/its", "MC/al/neg/pt/its",
                                                           "MC/all/neg/pt/its"};
  static constexpr std::string_view hPtTpc[nHistograms] = {"MC/el/pos/pt/tpc", "MC/mu/pos/pt/tpc", "MC/pi/pos/pt/tpc",
                                                           "MC/ka/pos/pt/tpc", "MC/pr/pos/pt/tpc", "MC/de/pos/pt/tpc",
                                                           "MC/tr/pos/pt/tpc", "MC/he/pos/pt/tpc", "MC/al/pos/pt/tpc",
                                                           "MC/all/pos/pt/tpc",
                                                           "MC/el/neg/pt/tpc", "MC/mu/neg/pt/tpc", "MC/pi/neg/pt/tpc",
                                                           "MC/ka/neg/pt/tpc", "MC/pr/neg/pt/tpc", "MC/de/neg/pt/tpc",
                                                           "MC/tr/neg/pt/tpc", "MC/he/neg/pt/tpc", "MC/al/neg/pt/tpc",
                                                           "MC/all/neg/pt/tpc"};
  static constexpr std::string_view hPtItsTpc[nHistograms] = {"MC/el/pos/pt/its_tpc", "MC/mu/pos/pt/its_tpc", "MC/pi/pos/pt/its_tpc",
                                                              "MC/ka/pos/pt/its_tpc", "MC/pr/pos/pt/its_tpc", "MC/de/pos/pt/its_tpc",
                                                              "MC/tr/pos/pt/its_tpc", "MC/he/pos/pt/its_tpc", "MC/al/pos/pt/its_tpc",
                                                              "MC/all/pos/pt/its_tpc",
                                                              "MC/el/neg/pt/its_tpc", "MC/mu/neg/pt/its_tpc", "MC/pi/neg/pt/its_tpc",
                                                              "MC/ka/neg/pt/its_tpc", "MC/pr/neg/pt/its_tpc", "MC/de/neg/pt/its_tpc",
                                                              "MC/tr/neg/pt/its_tpc", "MC/he/neg/pt/its_tpc", "MC/al/neg/pt/its_tpc",
                                                              "MC/all/neg/pt/its_tpc"};
  static constexpr std::string_view hPtItsTof[nHistograms] = {"MC/el/pos/pt/its_tof", "MC/mu/pos/pt/its_tof", "MC/pi/pos/pt/its_tof",
                                                              "MC/ka/pos/pt/its_tof", "MC/pr/pos/pt/its_tof", "MC/de/pos/pt/its_tof",
                                                              "MC/tr/pos/pt/its_tof", "MC/he/pos/pt/its_tof", "MC/al/pos/pt/its_tof",
                                                              "MC/all/pos/pt/its_tof",
                                                              "MC/el/neg/pt/its_tof", "MC/mu/neg/pt/its_tof", "MC/pi/neg/pt/its_tof",
                                                              "MC/ka/neg/pt/its_tof", "MC/pr/neg/pt/its_tof", "MC/de/neg/pt/its_tof",
                                                              "MC/tr/neg/pt/its_tof", "MC/he/neg/pt/its_tof", "MC/al/neg/pt/its_tof",
                                                              "MC/all/neg/pt/its_tof"};
  static constexpr std::string_view hPtTpcTof[nHistograms] = {"MC/el/pos/pt/tpc_tof", "MC/mu/pos/pt/tpc_tof", "MC/pi/pos/pt/tpc_tof",
                                                              "MC/ka/pos/pt/tpc_tof", "MC/pr/pos/pt/tpc_tof", "MC/de/pos/pt/tpc_tof",
                                                              "MC/tr/pos/pt/tpc_tof", "MC/he/pos/pt/tpc_tof", "MC/al/pos/pt/tpc_tof",
                                                              "MC/all/pos/pt/tpc_tof",
                                                              "MC/el/neg/pt/tpc_tof", "MC/mu/neg/pt/tpc_tof", "MC/pi/neg/pt/tpc_tof",
                                                              "MC/ka/neg/pt/tpc_tof", "MC/pr/neg/pt/tpc_tof", "MC/de/neg/pt/tpc_tof",
                                                              "MC/tr/neg/pt/tpc_tof", "MC/he/neg/pt/tpc_tof", "MC/al/neg/pt/tpc_tof",
                                                              "MC/all/neg/pt/tpc_tof"};
  static constexpr std::string_view hPtItsTpcTof[nHistograms] = {"MC/el/pos/pt/its_tpc_tof", "MC/mu/pos/pt/its_tpc_tof", "MC/pi/pos/pt/its_tpc_tof",
                                                                 "MC/ka/pos/pt/its_tpc_tof", "MC/pr/pos/pt/its_tpc_tof", "MC/de/pos/pt/its_tpc_tof",
                                                                 "MC/tr/pos/pt/its_tpc_tof", "MC/he/pos/pt/its_tpc_tof", "MC/al/pos/pt/its_tpc_tof",
                                                                 "MC/all/pos/pt/its_tpc_tof",
                                                                 "MC/el/neg/pt/its_tpc_tof", "MC/mu/neg/pt/its_tpc_tof", "MC/pi/neg/pt/its_tpc_tof",
                                                                 "MC/ka/neg/pt/its_tpc_tof", "MC/pr/neg/pt/its_tpc_tof", "MC/de/neg/pt/its_tpc_tof",
                                                                 "MC/tr/neg/pt/its_tpc_tof", "MC/he/neg/pt/its_tpc_tof", "MC/al/neg/pt/its_tpc_tof",
                                                                 "MC/all/neg/pt/its_tpc_tof"};
  static constexpr std::string_view hPtTrkItsTpc[nHistograms] = {"MC/el/pos/pt/its_tpc/trk", "MC/mu/pos/pt/its_tpc/trk", "MC/pi/pos/pt/its_tpc/trk",
                                                                 "MC/ka/pos/pt/its_tpc/trk", "MC/pr/pos/pt/its_tpc/trk", "MC/de/pos/pt/its_tpc/trk",
                                                                 "MC/tr/pos/pt/its_tpc/trk", "MC/he/pos/pt/its_tpc/trk", "MC/al/pos/pt/its_tpc/trk",
                                                                 "MC/all/pos/pt/its_tpc/trk",
                                                                 "MC/el/neg/pt/its_tpc/trk", "MC/mu/neg/pt/its_tpc/trk", "MC/pi/neg/pt/its_tpc/trk",
                                                                 "MC/ka/neg/pt/its_tpc/trk", "MC/pr/neg/pt/its_tpc/trk", "MC/de/neg/pt/its_tpc/trk",
                                                                 "MC/tr/neg/pt/its_tpc/trk", "MC/he/neg/pt/its_tpc/trk", "MC/al/neg/pt/its_tpc/trk",
                                                                 "MC/all/neg/pt/its_tpc/trk"};
  static constexpr std::string_view hPtGenerated[nHistograms] = {"MC/el/pos/pt/generated", "MC/mu/pos/pt/generated", "MC/pi/pos/pt/generated",
                                                                 "MC/ka/pos/pt/generated", "MC/pr/pos/pt/generated", "MC/de/pos/pt/generated",
                                                                 "MC/tr/pos/pt/generated", "MC/he/pos/pt/generated", "MC/al/pos/pt/generated",
                                                                 "MC/all/pos/pt/generated",
                                                                 "MC/el/neg/pt/generated", "MC/mu/neg/pt/generated", "MC/pi/neg/pt/generated",
                                                                 "MC/ka/neg/pt/generated", "MC/pr/neg/pt/generated", "MC/de/neg/pt/generated",
                                                                 "MC/tr/neg/pt/generated", "MC/he/neg/pt/generated", "MC/al/neg/pt/generated",
                                                                 "MC/all/neg/pt/generated"};

  // Pt for primaries
  static constexpr std::string_view hPtItsTpcPrm[nHistograms] = {"MC/el/pos/pt/its_tpc/prm", "MC/mu/pos/pt/its_tpc/prm", "MC/pi/pos/pt/its_tpc/prm",
                                                                 "MC/ka/pos/pt/its_tpc/prm", "MC/pr/pos/pt/its_tpc/prm", "MC/de/pos/pt/its_tpc/prm",
                                                                 "MC/tr/pos/pt/its_tpc/prm", "MC/he/pos/pt/its_tpc/prm", "MC/al/pos/pt/its_tpc/prm",
                                                                 "MC/all/pos/pt/its_tpc/prm",
                                                                 "MC/el/neg/pt/its_tpc/prm", "MC/mu/neg/pt/its_tpc/prm", "MC/pi/neg/pt/its_tpc/prm",
                                                                 "MC/ka/neg/pt/its_tpc/prm", "MC/pr/neg/pt/its_tpc/prm", "MC/de/neg/pt/its_tpc/prm",
                                                                 "MC/tr/neg/pt/its_tpc/prm", "MC/he/neg/pt/its_tpc/prm", "MC/al/neg/pt/its_tpc/prm",
                                                                 "MC/all/neg/pt/its_tpc/prm"};
  static constexpr std::string_view hPtTrkItsTpcPrm[nHistograms] = {"MC/el/pos/pt/its_tpc/trk/prm", "MC/mu/pos/pt/its_tpc/trk/prm", "MC/pi/pos/pt/its_tpc/trk/prm",
                                                                    "MC/ka/pos/pt/its_tpc/trk/prm", "MC/pr/pos/pt/its_tpc/trk/prm", "MC/de/pos/pt/its_tpc/trk/prm",
                                                                    "MC/tr/pos/pt/its_tpc/trk/prm", "MC/he/pos/pt/its_tpc/trk/prm", "MC/al/pos/pt/its_tpc/trk/prm",
                                                                    "MC/all/pos/pt/its_tpc/trk/prm",
                                                                    "MC/el/neg/pt/its_tpc/trk/prm", "MC/mu/neg/pt/its_tpc/trk/prm", "MC/pi/neg/pt/its_tpc/trk/prm",
                                                                    "MC/ka/neg/pt/its_tpc/trk/prm", "MC/pr/neg/pt/its_tpc/trk/prm", "MC/de/neg/pt/its_tpc/trk/prm",
                                                                    "MC/tr/neg/pt/its_tpc/trk/prm", "MC/he/neg/pt/its_tpc/trk/prm", "MC/al/neg/pt/its_tpc/trk/prm",
                                                                    "MC/all/neg/pt/its_tpc/trk/prm"};
  static constexpr std::string_view hPtItsTpcTofPrm[nHistograms] = {"MC/el/pos/pt/its_tpc_tof/prm", "MC/mu/pos/pt/its_tpc_tof/prm", "MC/pi/pos/pt/its_tpc_tof/prm",
                                                                    "MC/ka/pos/pt/its_tpc_tof/prm", "MC/pr/pos/pt/its_tpc_tof/prm", "MC/de/pos/pt/its_tpc_tof/prm",
                                                                    "MC/tr/pos/pt/its_tpc_tof/prm", "MC/he/pos/pt/its_tpc_tof/prm", "MC/al/pos/pt/its_tpc_tof/prm",
                                                                    "MC/all/pos/pt/its_tpc_tof/prm",
                                                                    "MC/el/neg/pt/its_tpc_tof/prm", "MC/mu/neg/pt/its_tpc_tof/prm", "MC/pi/neg/pt/its_tpc_tof/prm",
                                                                    "MC/ka/neg/pt/its_tpc_tof/prm", "MC/pr/neg/pt/its_tpc_tof/prm", "MC/de/neg/pt/its_tpc_tof/prm",
                                                                    "MC/tr/neg/pt/its_tpc_tof/prm", "MC/he/neg/pt/its_tpc_tof/prm", "MC/al/neg/pt/its_tpc_tof/prm",
                                                                    "MC/all/neg/pt/its_tpc_tof/prm"};
  static constexpr std::string_view hPtGeneratedPrm[nHistograms] = {"MC/el/pos/pt/generated/prm", "MC/mu/pos/pt/generated/prm", "MC/pi/pos/pt/generated/prm",
                                                                    "MC/ka/pos/pt/generated/prm", "MC/pr/pos/pt/generated/prm", "MC/de/pos/pt/generated/prm",
                                                                    "MC/tr/pos/pt/generated/prm", "MC/he/pos/pt/generated/prm", "MC/al/pos/pt/generated/prm",
                                                                    "MC/all/pos/pt/generated/prm",
                                                                    "MC/el/neg/pt/generated/prm", "MC/mu/neg/pt/generated/prm", "MC/pi/neg/pt/generated/prm",
                                                                    "MC/ka/neg/pt/generated/prm", "MC/pr/neg/pt/generated/prm", "MC/de/neg/pt/generated/prm",
                                                                    "MC/tr/neg/pt/generated/prm", "MC/he/neg/pt/generated/prm", "MC/al/neg/pt/generated/prm",
                                                                    "MC/all/neg/pt/generated/prm"};

  // Pt for secondaries from weak decay
  static constexpr std::string_view hPtItsTpcStr[nHistograms] = {"MC/el/pos/pt/its_tpc/str", "MC/mu/pos/pt/its_tpc/str", "MC/pi/pos/pt/its_tpc/str",
                                                                 "MC/ka/pos/pt/its_tpc/str", "MC/pr/pos/pt/its_tpc/str", "MC/de/pos/pt/its_tpc/str",
                                                                 "MC/tr/pos/pt/its_tpc/str", "MC/he/pos/pt/its_tpc/str", "MC/al/pos/pt/its_tpc/str",
                                                                 "MC/all/pos/pt/its_tpc/str",
                                                                 "MC/el/neg/pt/its_tpc/str", "MC/mu/neg/pt/its_tpc/str", "MC/pi/neg/pt/its_tpc/str",
                                                                 "MC/ka/neg/pt/its_tpc/str", "MC/pr/neg/pt/its_tpc/str", "MC/de/neg/pt/its_tpc/str",
                                                                 "MC/tr/neg/pt/its_tpc/str", "MC/he/neg/pt/its_tpc/str", "MC/al/neg/pt/its_tpc/str",
                                                                 "MC/all/neg/pt/its_tpc/str"};
  static constexpr std::string_view hPtTrkItsTpcStr[nHistograms] = {"MC/el/pos/pt/its_tpc/trk/str", "MC/mu/pos/pt/its_tpc/trk/str", "MC/pi/pos/pt/its_tpc/trk/str",
                                                                    "MC/ka/pos/pt/its_tpc/trk/str", "MC/pr/pos/pt/its_tpc/trk/str", "MC/de/pos/pt/its_tpc/trk/str",
                                                                    "MC/tr/pos/pt/its_tpc/trk/str", "MC/he/pos/pt/its_tpc/trk/str", "MC/al/pos/pt/its_tpc/trk/str",
                                                                    "MC/all/pos/pt/its_tpc/trk/str",
                                                                    "MC/el/neg/pt/its_tpc/trk/str", "MC/mu/neg/pt/its_tpc/trk/str", "MC/pi/neg/pt/its_tpc/trk/str",
                                                                    "MC/ka/neg/pt/its_tpc/trk/str", "MC/pr/neg/pt/its_tpc/trk/str", "MC/de/neg/pt/its_tpc/trk/str",
                                                                    "MC/tr/neg/pt/its_tpc/trk/str", "MC/he/neg/pt/its_tpc/trk/str", "MC/al/neg/pt/its_tpc/trk/str",
                                                                    "MC/all/neg/pt/its_tpc/trk/str"};
  static constexpr std::string_view hPtItsTpcTofStr[nHistograms] = {"MC/el/pos/pt/its_tpc_tof/str", "MC/mu/pos/pt/its_tpc_tof/str", "MC/pi/pos/pt/its_tpc_tof/str",
                                                                    "MC/ka/pos/pt/its_tpc_tof/str", "MC/pr/pos/pt/its_tpc_tof/str", "MC/de/pos/pt/its_tpc_tof/str",
                                                                    "MC/tr/pos/pt/its_tpc_tof/str", "MC/he/pos/pt/its_tpc_tof/str", "MC/al/pos/pt/its_tpc_tof/str",
                                                                    "MC/all/pos/pt/its_tpc_tof/str",
                                                                    "MC/el/neg/pt/its_tpc_tof/str", "MC/mu/neg/pt/its_tpc_tof/str", "MC/pi/neg/pt/its_tpc_tof/str",
                                                                    "MC/ka/neg/pt/its_tpc_tof/str", "MC/pr/neg/pt/its_tpc_tof/str", "MC/de/neg/pt/its_tpc_tof/str",
                                                                    "MC/tr/neg/pt/its_tpc_tof/str", "MC/he/neg/pt/its_tpc_tof/str", "MC/al/neg/pt/its_tpc_tof/str",
                                                                    "MC/all/neg/pt/its_tpc_tof/str"};
  static constexpr std::string_view hPtGeneratedStr[nHistograms] = {"MC/el/pos/pt/generated/str", "MC/mu/pos/pt/generated/str", "MC/pi/pos/pt/generated/str",
                                                                    "MC/ka/pos/pt/generated/str", "MC/pr/pos/pt/generated/str", "MC/de/pos/pt/generated/str",
                                                                    "MC/tr/pos/pt/generated/str", "MC/he/pos/pt/generated/str", "MC/al/pos/pt/generated/str",
                                                                    "MC/all/pos/pt/generated/str",
                                                                    "MC/el/neg/pt/generated/str", "MC/mu/neg/pt/generated/str", "MC/pi/neg/pt/generated/str",
                                                                    "MC/ka/neg/pt/generated/str", "MC/pr/neg/pt/generated/str", "MC/de/neg/pt/generated/str",
                                                                    "MC/tr/neg/pt/generated/str", "MC/he/neg/pt/generated/str", "MC/al/neg/pt/generated/str",
                                                                    "MC/all/neg/pt/generated/str"};

  // Pt for secondaries from material
  static constexpr std::string_view hPtItsTpcMat[nHistograms] = {"MC/el/pos/pt/its_tpc/mat", "MC/mu/pos/pt/its_tpc/mat", "MC/pi/pos/pt/its_tpc/mat",
                                                                 "MC/ka/pos/pt/its_tpc/mat", "MC/pr/pos/pt/its_tpc/mat", "MC/de/pos/pt/its_tpc/mat",
                                                                 "MC/tr/pos/pt/its_tpc/mat", "MC/he/pos/pt/its_tpc/mat", "MC/al/pos/pt/its_tpc/mat",
                                                                 "MC/all/pos/pt/its_tpc/mat",
                                                                 "MC/el/neg/pt/its_tpc/mat", "MC/mu/neg/pt/its_tpc/mat", "MC/pi/neg/pt/its_tpc/mat",
                                                                 "MC/ka/neg/pt/its_tpc/mat", "MC/pr/neg/pt/its_tpc/mat", "MC/de/neg/pt/its_tpc/mat",
                                                                 "MC/tr/neg/pt/its_tpc/mat", "MC/he/neg/pt/its_tpc/mat", "MC/al/neg/pt/its_tpc/mat",
                                                                 "MC/all/neg/pt/its_tpc/mat"};
  static constexpr std::string_view hPtTrkItsTpcMat[nHistograms] = {"MC/el/pos/pt/its_tpc/trk/mat", "MC/mu/pos/pt/its_tpc/trk/mat", "MC/pi/pos/pt/its_tpc/trk/mat",
                                                                    "MC/ka/pos/pt/its_tpc/trk/mat", "MC/pr/pos/pt/its_tpc/trk/mat", "MC/de/pos/pt/its_tpc/trk/mat",
                                                                    "MC/tr/pos/pt/its_tpc/trk/mat", "MC/he/pos/pt/its_tpc/trk/mat", "MC/al/pos/pt/its_tpc/trk/mat",
                                                                    "MC/all/pos/pt/its_tpc/trk/mat",
                                                                    "MC/el/neg/pt/its_tpc/trk/mat", "MC/mu/neg/pt/its_tpc/trk/mat", "MC/pi/neg/pt/its_tpc/trk/mat",
                                                                    "MC/ka/neg/pt/its_tpc/trk/mat", "MC/pr/neg/pt/its_tpc/trk/mat", "MC/de/neg/pt/its_tpc/trk/mat",
                                                                    "MC/tr/neg/pt/its_tpc/trk/mat", "MC/he/neg/pt/its_tpc/trk/mat", "MC/al/neg/pt/its_tpc/trk/mat",
                                                                    "MC/all/neg/pt/its_tpc/trk/mat"};
  static constexpr std::string_view hPtItsTpcTofMat[nHistograms] = {"MC/el/pos/pt/its_tpc_tof/mat", "MC/mu/pos/pt/its_tpc_tof/mat", "MC/pi/pos/pt/its_tpc_tof/mat",
                                                                    "MC/ka/pos/pt/its_tpc_tof/mat", "MC/pr/pos/pt/its_tpc_tof/mat", "MC/de/pos/pt/its_tpc_tof/mat",
                                                                    "MC/tr/pos/pt/its_tpc_tof/mat", "MC/he/pos/pt/its_tpc_tof/mat", "MC/al/pos/pt/its_tpc_tof/mat",
                                                                    "MC/all/pos/pt/its_tpc_tof/mat",
                                                                    "MC/el/neg/pt/its_tpc_tof/mat", "MC/mu/neg/pt/its_tpc_tof/mat", "MC/pi/neg/pt/its_tpc_tof/mat",
                                                                    "MC/ka/neg/pt/its_tpc_tof/mat", "MC/pr/neg/pt/its_tpc_tof/mat", "MC/de/neg/pt/its_tpc_tof/mat",
                                                                    "MC/tr/neg/pt/its_tpc_tof/mat", "MC/he/neg/pt/its_tpc_tof/mat", "MC/al/neg/pt/its_tpc_tof/mat",
                                                                    "MC/all/neg/pt/its_tpc_tof/mat"};
  static constexpr std::string_view hPtGeneratedMat[nHistograms] = {"MC/el/pos/pt/generated/mat", "MC/mu/pos/pt/generated/mat", "MC/pi/pos/pt/generated/mat",
                                                                    "MC/ka/pos/pt/generated/mat", "MC/pr/pos/pt/generated/mat", "MC/de/pos/pt/generated/mat",
                                                                    "MC/tr/pos/pt/generated/mat", "MC/he/pos/pt/generated/mat", "MC/al/pos/pt/generated/mat",
                                                                    "MC/all/pos/pt/generated/mat",
                                                                    "MC/el/neg/pt/generated/mat", "MC/mu/neg/pt/generated/mat", "MC/pi/neg/pt/generated/mat",
                                                                    "MC/ka/neg/pt/generated/mat", "MC/pr/neg/pt/generated/mat", "MC/de/neg/pt/generated/mat",
                                                                    "MC/tr/neg/pt/generated/mat", "MC/he/neg/pt/generated/mat", "MC/al/neg/pt/generated/mat",
                                                                    "MC/all/neg/pt/generated/mat"};

  // P
  static constexpr std::string_view hPItsTpc[nHistograms] = {"MC/el/pos/p/its_tpc", "MC/mu/pos/p/its_tpc", "MC/pi/pos/p/its_tpc",
                                                             "MC/ka/pos/p/its_tpc", "MC/pr/pos/p/its_tpc", "MC/de/pos/p/its_tpc",
                                                             "MC/tr/pos/p/its_tpc", "MC/he/pos/p/its_tpc", "MC/al/pos/p/its_tpc",
                                                             "MC/all/pos/p/its_tpc",
                                                             "MC/el/neg/p/its_tpc", "MC/mu/neg/p/its_tpc", "MC/pi/neg/p/its_tpc",
                                                             "MC/ka/neg/p/its_tpc", "MC/pr/neg/p/its_tpc", "MC/de/neg/p/its_tpc",
                                                             "MC/tr/neg/p/its_tpc", "MC/he/neg/p/its_tpc", "MC/al/neg/p/its_tpc",
                                                             "MC/all/neg/p/its_tpc"};
  static constexpr std::string_view hPTrkItsTpc[nHistograms] = {"MC/el/pos/p/its_tpc/trk", "MC/mu/pos/p/its_tpc/trk", "MC/pi/pos/p/its_tpc/trk",
                                                                "MC/ka/pos/p/its_tpc/trk", "MC/pr/pos/p/its_tpc/trk", "MC/de/pos/p/its_tpc/trk",
                                                                "MC/tr/pos/p/its_tpc/trk", "MC/he/pos/p/its_tpc/trk", "MC/al/pos/p/its_tpc/trk",
                                                                "MC/all/pos/p/its_tpc/trk",
                                                                "MC/el/neg/p/its_tpc/trk", "MC/mu/neg/p/its_tpc/trk", "MC/pi/neg/p/its_tpc/trk",
                                                                "MC/ka/neg/p/its_tpc/trk", "MC/pr/neg/p/its_tpc/trk", "MC/de/neg/p/its_tpc/trk",
                                                                "MC/tr/neg/p/its_tpc/trk", "MC/he/neg/p/its_tpc/trk", "MC/al/neg/p/its_tpc/trk",
                                                                "MC/all/neg/p/its_tpc/trk"};
  static constexpr std::string_view hPItsTpcTof[nHistograms] = {"MC/el/pos/p/its_tpc_tof", "MC/mu/pos/p/its_tpc_tof", "MC/pi/pos/p/its_tpc_tof",
                                                                "MC/ka/pos/p/its_tpc_tof", "MC/pr/pos/p/its_tpc_tof", "MC/de/pos/p/its_tpc_tof",
                                                                "MC/tr/pos/p/its_tpc_tof", "MC/he/pos/p/its_tpc_tof", "MC/al/pos/p/its_tpc_tof",
                                                                "MC/all/pos/p/its_tpc_tof",
                                                                "MC/el/neg/p/its_tpc_tof", "MC/mu/neg/p/its_tpc_tof", "MC/pi/neg/p/its_tpc_tof",
                                                                "MC/ka/neg/p/its_tpc_tof", "MC/pr/neg/p/its_tpc_tof", "MC/de/neg/p/its_tpc_tof",
                                                                "MC/tr/neg/p/its_tpc_tof", "MC/he/neg/p/its_tpc_tof", "MC/al/neg/p/its_tpc_tof",
                                                                "MC/all/neg/p/its_tpc_tof"};
  static constexpr std::string_view hPGenerated[nHistograms] = {"MC/el/pos/p/generated", "MC/mu/pos/p/generated", "MC/pi/pos/p/generated",
                                                                "MC/ka/pos/p/generated", "MC/pr/pos/p/generated", "MC/de/pos/p/generated",
                                                                "MC/tr/pos/p/generated", "MC/he/pos/p/generated", "MC/al/pos/p/generated",
                                                                "MC/all/pos/p/generated",
                                                                "MC/el/neg/p/generated", "MC/mu/neg/p/generated", "MC/pi/neg/p/generated",
                                                                "MC/ka/neg/p/generated", "MC/pr/neg/p/generated", "MC/de/neg/p/generated",
                                                                "MC/tr/neg/p/generated", "MC/he/neg/p/generated", "MC/al/neg/p/generated",
                                                                "MC/all/neg/p/generated"};

  // Eta
  static constexpr std::string_view hEtaItsTpc[nHistograms] = {"MC/el/pos/eta/its_tpc", "MC/mu/pos/eta/its_tpc", "MC/pi/pos/eta/its_tpc",
                                                               "MC/ka/pos/eta/its_tpc", "MC/pr/pos/eta/its_tpc", "MC/de/pos/eta/its_tpc",
                                                               "MC/tr/pos/eta/its_tpc", "MC/he/pos/eta/its_tpc", "MC/al/pos/eta/its_tpc",
                                                               "MC/all/pos/eta/its_tpc",
                                                               "MC/el/neg/eta/its_tpc", "MC/mu/neg/eta/its_tpc", "MC/pi/neg/eta/its_tpc",
                                                               "MC/ka/neg/eta/its_tpc", "MC/pr/neg/eta/its_tpc", "MC/de/neg/eta/its_tpc",
                                                               "MC/tr/neg/eta/its_tpc", "MC/he/neg/eta/its_tpc", "MC/al/neg/eta/its_tpc",
                                                               "MC/all/neg/eta/its_tpc"};
  static constexpr std::string_view hEtaTrkItsTpc[nHistograms] = {"MC/el/pos/eta/its_tpc/trk", "MC/mu/pos/eta/its_tpc/trk", "MC/pi/pos/eta/its_tpc/trk",
                                                                  "MC/ka/pos/eta/its_tpc/trk", "MC/pr/pos/eta/its_tpc/trk", "MC/de/pos/eta/its_tpc/trk",
                                                                  "MC/tr/pos/eta/its_tpc/trk", "MC/he/pos/eta/its_tpc/trk", "MC/al/pos/eta/its_tpc/trk",
                                                                  "MC/all/pos/eta/its_tpc/trk",
                                                                  "MC/el/neg/eta/its_tpc/trk", "MC/mu/neg/eta/its_tpc/trk", "MC/pi/neg/eta/its_tpc/trk",
                                                                  "MC/ka/neg/eta/its_tpc/trk", "MC/pr/neg/eta/its_tpc/trk", "MC/de/neg/eta/its_tpc/trk",
                                                                  "MC/tr/neg/eta/its_tpc/trk", "MC/he/neg/eta/its_tpc/trk", "MC/al/neg/eta/its_tpc/trk",
                                                                  "MC/all/neg/eta/its_tpc/trk"};
  static constexpr std::string_view hEtaItsTpcTof[nHistograms] = {"MC/el/pos/eta/its_tpc_tof", "MC/mu/pos/eta/its_tpc_tof", "MC/pi/pos/eta/its_tpc_tof",
                                                                  "MC/ka/pos/eta/its_tpc_tof", "MC/pr/pos/eta/its_tpc_tof", "MC/de/pos/eta/its_tpc_tof",
                                                                  "MC/tr/pos/eta/its_tpc_tof", "MC/he/pos/eta/its_tpc_tof", "MC/al/pos/eta/its_tpc_tof",
                                                                  "MC/all/pos/eta/its_tpc_tof",
                                                                  "MC/el/neg/eta/its_tpc_tof", "MC/mu/neg/eta/its_tpc_tof", "MC/pi/neg/eta/its_tpc_tof",
                                                                  "MC/ka/neg/eta/its_tpc_tof", "MC/pr/neg/eta/its_tpc_tof", "MC/de/neg/eta/its_tpc_tof",
                                                                  "MC/tr/neg/eta/its_tpc_tof", "MC/he/neg/eta/its_tpc_tof", "MC/al/neg/eta/its_tpc_tof",
                                                                  "MC/all/neg/eta/its_tpc_tof"};
  static constexpr std::string_view hEtaGenerated[nHistograms] = {"MC/el/pos/eta/generated", "MC/mu/pos/eta/generated", "MC/pi/pos/eta/generated",
                                                                  "MC/ka/pos/eta/generated", "MC/pr/pos/eta/generated", "MC/de/pos/eta/generated",
                                                                  "MC/tr/pos/eta/generated", "MC/he/pos/eta/generated", "MC/al/pos/eta/generated",
                                                                  "MC/all/pos/eta/generated",
                                                                  "MC/el/neg/eta/generated", "MC/mu/neg/eta/generated", "MC/pi/neg/eta/generated",
                                                                  "MC/ka/neg/eta/generated", "MC/pr/neg/eta/generated", "MC/de/neg/eta/generated",
                                                                  "MC/tr/neg/eta/generated", "MC/he/neg/eta/generated", "MC/al/neg/eta/generated",
                                                                  "MC/all/neg/eta/generated"};

  // Y
  static constexpr std::string_view hYItsTpc[nHistograms] = {"MC/el/pos/y/its_tpc", "MC/mu/pos/y/its_tpc", "MC/pi/pos/y/its_tpc",
                                                             "MC/ka/pos/y/its_tpc", "MC/pr/pos/y/its_tpc", "MC/de/pos/y/its_tpc",
                                                             "MC/tr/pos/y/its_tpc", "MC/he/pos/y/its_tpc", "MC/al/pos/y/its_tpc",
                                                             "MC/all/pos/y/its_tpc",
                                                             "MC/el/neg/y/its_tpc", "MC/mu/neg/y/its_tpc", "MC/pi/neg/y/its_tpc",
                                                             "MC/ka/neg/y/its_tpc", "MC/pr/neg/y/its_tpc", "MC/de/neg/y/its_tpc",
                                                             "MC/tr/neg/y/its_tpc", "MC/he/neg/y/its_tpc", "MC/al/neg/y/its_tpc",
                                                             "MC/all/neg/y/its_tpc"};
  static constexpr std::string_view hYItsTpcTof[nHistograms] = {"MC/el/pos/y/its_tpc_tof", "MC/mu/pos/y/its_tpc_tof", "MC/pi/pos/y/its_tpc_tof",
                                                                "MC/ka/pos/y/its_tpc_tof", "MC/pr/pos/y/its_tpc_tof", "MC/de/pos/y/its_tpc_tof",
                                                                "MC/tr/pos/y/its_tpc_tof", "MC/he/pos/y/its_tpc_tof", "MC/al/pos/y/its_tpc_tof",
                                                                "MC/all/pos/y/its_tpc_tof",
                                                                "MC/el/neg/y/its_tpc_tof", "MC/mu/neg/y/its_tpc_tof", "MC/pi/neg/y/its_tpc_tof",
                                                                "MC/ka/neg/y/its_tpc_tof", "MC/pr/neg/y/its_tpc_tof", "MC/de/neg/y/its_tpc_tof",
                                                                "MC/tr/neg/y/its_tpc_tof", "MC/he/neg/y/its_tpc_tof", "MC/al/neg/y/its_tpc_tof",
                                                                "MC/all/neg/y/its_tpc_tof"};
  static constexpr std::string_view hYGenerated[nHistograms] = {"MC/el/pos/y/generated", "MC/mu/pos/y/generated", "MC/pi/pos/y/generated",
                                                                "MC/ka/pos/y/generated", "MC/pr/pos/y/generated", "MC/de/pos/y/generated",
                                                                "MC/tr/pos/y/generated", "MC/he/pos/y/generated", "MC/al/pos/y/generated",
                                                                "MC/all/pos/y/generated",
                                                                "MC/el/neg/y/generated", "MC/mu/neg/y/generated", "MC/pi/neg/y/generated",
                                                                "MC/ka/neg/y/generated", "MC/pr/neg/y/generated", "MC/de/neg/y/generated",
                                                                "MC/tr/neg/y/generated", "MC/he/neg/y/generated", "MC/al/neg/y/generated",
                                                                "MC/all/neg/y/generated"};

  // Phi
  static constexpr std::string_view hPhiItsTpc[nHistograms] = {"MC/el/pos/phi/its_tpc", "MC/mu/pos/phi/its_tpc", "MC/pi/pos/phi/its_tpc",
                                                               "MC/ka/pos/phi/its_tpc", "MC/pr/pos/phi/its_tpc", "MC/de/pos/phi/its_tpc",
                                                               "MC/tr/pos/phi/its_tpc", "MC/he/pos/phi/its_tpc", "MC/al/pos/phi/its_tpc",
                                                               "MC/all/pos/phi/its_tpc",
                                                               "MC/el/neg/phi/its_tpc", "MC/mu/neg/phi/its_tpc", "MC/pi/neg/phi/its_tpc",
                                                               "MC/ka/neg/phi/its_tpc", "MC/pr/neg/phi/its_tpc", "MC/de/neg/phi/its_tpc",
                                                               "MC/tr/neg/phi/its_tpc", "MC/he/neg/phi/its_tpc", "MC/al/neg/phi/its_tpc",
                                                               "MC/all/neg/phi/its_tpc"};
  static constexpr std::string_view hPhiTrkItsTpc[nHistograms] = {"MC/el/pos/phi/its_tpc/trk", "MC/mu/pos/phi/its_tpc/trk", "MC/pi/pos/phi/its_tpc/trk",
                                                                  "MC/ka/pos/phi/its_tpc/trk", "MC/pr/pos/phi/its_tpc/trk", "MC/de/pos/phi/its_tpc/trk",
                                                                  "MC/tr/pos/phi/its_tpc/trk", "MC/he/pos/phi/its_tpc/trk", "MC/al/pos/phi/its_tpc/trk",
                                                                  "MC/all/pos/phi/its_tpc/trk",
                                                                  "MC/el/neg/phi/its_tpc/trk", "MC/mu/neg/phi/its_tpc/trk", "MC/pi/neg/phi/its_tpc/trk",
                                                                  "MC/ka/neg/phi/its_tpc/trk", "MC/pr/neg/phi/its_tpc/trk", "MC/de/neg/phi/its_tpc/trk",
                                                                  "MC/tr/neg/phi/its_tpc/trk", "MC/he/neg/phi/its_tpc/trk", "MC/al/neg/phi/its_tpc/trk",
                                                                  "MC/all/neg/phi/its_tpc/trk"};
  static constexpr std::string_view hPhiItsTpcTof[nHistograms] = {"MC/el/pos/phi/its_tpc_tof", "MC/mu/pos/phi/its_tpc_tof", "MC/pi/pos/phi/its_tpc_tof",
                                                                  "MC/ka/pos/phi/its_tpc_tof", "MC/pr/pos/phi/its_tpc_tof", "MC/de/pos/phi/its_tpc_tof",
                                                                  "MC/tr/pos/phi/its_tpc_tof", "MC/he/pos/phi/its_tpc_tof", "MC/al/pos/phi/its_tpc_tof",
                                                                  "MC/all/pos/phi/its_tpc_tof",
                                                                  "MC/el/neg/phi/its_tpc_tof", "MC/mu/neg/phi/its_tpc_tof", "MC/pi/neg/phi/its_tpc_tof",
                                                                  "MC/ka/neg/phi/its_tpc_tof", "MC/pr/neg/phi/its_tpc_tof", "MC/de/neg/phi/its_tpc_tof",
                                                                  "MC/tr/neg/phi/its_tpc_tof", "MC/he/neg/phi/its_tpc_tof", "MC/al/neg/phi/its_tpc_tof",
                                                                  "MC/all/neg/phi/its_tpc_tof"};
  static constexpr std::string_view hPhiGenerated[nHistograms] = {"MC/el/pos/phi/generated", "MC/mu/pos/phi/generated", "MC/pi/pos/phi/generated",
                                                                  "MC/ka/pos/phi/generated", "MC/pr/pos/phi/generated", "MC/de/pos/phi/generated",
                                                                  "MC/tr/pos/phi/generated", "MC/he/pos/phi/generated", "MC/al/pos/phi/generated",
                                                                  "MC/all/pos/phi/generated",
                                                                  "MC/el/neg/phi/generated", "MC/mu/neg/phi/generated", "MC/pi/neg/phi/generated",
                                                                  "MC/ka/neg/phi/generated", "MC/pr/neg/phi/generated", "MC/de/neg/phi/generated",
                                                                  "MC/tr/neg/phi/generated", "MC/he/neg/phi/generated", "MC/al/neg/phi/generated",
                                                                  "MC/all/neg/phi/generated"};

  // Pt-Eta
  static constexpr std::string_view hPtEtaItsTpc[nHistograms] = {"MC/el/pos/pteta/its_tpc", "MC/mu/pos/pteta/its_tpc", "MC/pi/pos/pteta/its_tpc",
                                                                 "MC/ka/pos/pteta/its_tpc", "MC/pr/pos/pteta/its_tpc", "MC/de/pos/pteta/its_tpc",
                                                                 "MC/tr/pos/pteta/its_tpc", "MC/he/pos/pteta/its_tpc", "MC/al/pos/pteta/its_tpc",
                                                                 "MC/all/pos/pteta/its_tpc",
                                                                 "MC/el/neg/pteta/its_tpc", "MC/mu/neg/pteta/its_tpc", "MC/pi/neg/pteta/its_tpc",
                                                                 "MC/ka/neg/pteta/its_tpc", "MC/pr/neg/pteta/its_tpc", "MC/de/neg/pteta/its_tpc",
                                                                 "MC/tr/neg/pteta/its_tpc", "MC/he/neg/pteta/its_tpc", "MC/al/neg/pteta/its_tpc",
                                                                 "MC/all/neg/pteta/its_tpc"};
  static constexpr std::string_view hPtEtaTrkItsTpc[nHistograms] = {"MC/el/pos/pteta/its_tpc/trk", "MC/mu/pos/pteta/its_tpc/trk", "MC/pi/pos/pteta/its_tpc/trk",
                                                                    "MC/ka/pos/pteta/its_tpc/trk", "MC/pr/pos/pteta/its_tpc/trk", "MC/de/pos/pteta/its_tpc/trk",
                                                                    "MC/tr/pos/pteta/its_tpc/trk", "MC/he/pos/pteta/its_tpc/trk", "MC/al/pos/pteta/its_tpc/trk",
                                                                    "MC/all/pos/pteta/its_tpc/trk",
                                                                    "MC/el/neg/pteta/its_tpc/trk", "MC/mu/neg/pteta/its_tpc/trk", "MC/pi/neg/pteta/its_tpc/trk",
                                                                    "MC/ka/neg/pteta/its_tpc/trk", "MC/pr/neg/pteta/its_tpc/trk", "MC/de/neg/pteta/its_tpc/trk",
                                                                    "MC/tr/neg/pteta/its_tpc/trk", "MC/he/neg/pteta/its_tpc/trk", "MC/al/neg/pteta/its_tpc/trk",
                                                                    "MC/all/neg/pteta/its_tpc/trk"};
  static constexpr std::string_view hPtEtaItsTpcTof[nHistograms] = {"MC/el/pos/pteta/its_tpc_tof", "MC/mu/pos/pteta/its_tpc_tof", "MC/pi/pos/pteta/its_tpc_tof",
                                                                    "MC/ka/pos/pteta/its_tpc_tof", "MC/pr/pos/pteta/its_tpc_tof", "MC/de/pos/pteta/its_tpc_tof",
                                                                    "MC/tr/pos/pteta/its_tpc_tof", "MC/he/pos/pteta/its_tpc_tof", "MC/al/pos/pteta/its_tpc_tof",
                                                                    "MC/all/pos/pteta/its_tpc_tof",
                                                                    "MC/el/neg/pteta/its_tpc_tof", "MC/mu/neg/pteta/its_tpc_tof", "MC/pi/neg/pteta/its_tpc_tof",
                                                                    "MC/ka/neg/pteta/its_tpc_tof", "MC/pr/neg/pteta/its_tpc_tof", "MC/de/neg/pteta/its_tpc_tof",
                                                                    "MC/tr/neg/pteta/its_tpc_tof", "MC/he/neg/pteta/its_tpc_tof", "MC/al/neg/pteta/its_tpc_tof",
                                                                    "MC/all/neg/pteta/its_tpc_tof"};
  static constexpr std::string_view hPtEtaGenerated[nHistograms] = {"MC/el/pos/pteta/generated", "MC/mu/pos/pteta/generated", "MC/pi/pos/pteta/generated",
                                                                    "MC/ka/pos/pteta/generated", "MC/pr/pos/pteta/generated", "MC/de/pos/pteta/generated",
                                                                    "MC/tr/pos/pteta/generated", "MC/he/pos/pteta/generated", "MC/al/pos/pteta/generated",
                                                                    "MC/all/pos/pteta/generated",
                                                                    "MC/el/neg/pteta/generated", "MC/mu/neg/pteta/generated", "MC/pi/neg/pteta/generated",
                                                                    "MC/ka/neg/pteta/generated", "MC/pr/neg/pteta/generated", "MC/de/neg/pteta/generated",
                                                                    "MC/tr/neg/pteta/generated", "MC/he/neg/pteta/generated", "MC/al/neg/pteta/generated",
                                                                    "MC/all/neg/pteta/generated"};

  template <o2::track::PID::ID id>
  void makeMCHistograms(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
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

    const char* partName = id == o2::track::PID::NIDs ? "All" : o2::track::PID::getName(id);
    LOG(debug) << "Preparing histograms for particle: " << partName;

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

    auto makeHistogramsPerCharge = [&](const int chargeIndex) {
      const int histogramIndex = id + chargeIndex * nSpecies;

      histos.add(hPtIts[histogramIndex].data(), "ITS tracks " + tagPt, kTH1F, {axisPt});
      histos.add(hPtTpc[histogramIndex].data(), "TPC tracks " + tagPt, kTH1F, {axisPt});
      histos.add(hPtItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPt, kTH1F, {axisPt});
      histos.add(hPtItsTof[histogramIndex].data(), "ITS-TOF tracks " + tagPt, kTH1F, {axisPt});
      histos.add(hPtTpcTof[histogramIndex].data(), "TPC-TOF tracks " + tagPt, kTH1F, {axisPt});
      histos.add(hPtItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPt, kTH1F, {axisPt});
      histos.add(hPtTrkItsTpc[histogramIndex].data(), "ITS-TPC track (reco) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtGenerated[histogramIndex].data(), "Generated " + tagPt, kTH1F, {axisPt});

      histos.add(hPtItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (primaries) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtTrkItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (reco primaries) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtItsTpcTofPrm[histogramIndex].data(), "ITS-TPC-TOF tracks (primaries) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtGeneratedPrm[histogramIndex].data(), "Generated (primaries) " + tagPt, kTH1F, {axisPt});

      histos.add(hPtItsTpcStr[histogramIndex].data(), "ITS-TPC tracks (from weak decays) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtTrkItsTpcStr[histogramIndex].data(), "ITS-TPC tracks (reco from weak decays) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtItsTpcTofStr[histogramIndex].data(), "ITS-TPC-TOF tracks (from weak decays) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtGeneratedStr[histogramIndex].data(), "Generated (from weak decays) " + tagPt, kTH1F, {axisPt});

      histos.add(hPtItsTpcMat[histogramIndex].data(), "ITS-TPC tracks (from material)" + tagPt, kTH1F, {axisPt});
      histos.add(hPtTrkItsTpcMat[histogramIndex].data(), "ITS-TPC tracks (reco from material) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtItsTpcTofMat[histogramIndex].data(), "ITS-TPC-TOF tracks ( from material) " + tagPt, kTH1F, {axisPt});
      histos.add(hPtGeneratedMat[histogramIndex].data(), "Generated ( from material) " + tagPt, kTH1F, {axisPt});

      histos.add(hPItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPt, kTH1F, {axisP});
      histos.add(hPTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagPt, kTH1F, {axisP});
      histos.add(hPItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPt, kTH1F, {axisP});
      histos.add(hPGenerated[histogramIndex].data(), "Generated " + tagPt, kTH1F, {axisP});

      histos.add(hEtaItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagEta, kTH1F, {axisEta});
      histos.add(hEtaTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagEta, kTH1F, {axisEta});
      histos.add(hEtaItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagEta, kTH1F, {axisEta});
      histos.add(hEtaGenerated[histogramIndex].data(), "Generated " + tagEta, kTH1F, {axisEta});

      histos.add(hYItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagY, kTH1F, {axisY});
      histos.add(hYItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagY, kTH1F, {axisY});
      histos.add(hYGenerated[histogramIndex].data(), "Generated " + tagY, kTH1F, {axisY});

      histos.add(hPhiItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPhi, kTH1F, {axisPhi});
      histos.add(hPhiTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagPhi, kTH1F, {axisPhi});
      histos.add(hPhiItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPhi, kTH1F, {axisPhi});
      histos.add(hPhiGenerated[histogramIndex].data(), "Generated " + tagPhi, kTH1F, {axisPhi});

      if (doPtEta) {
        histos.add(hPtEtaItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPtEta, kTH2D, {axisPt, axisEta});
        histos.add(hPtEtaTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagPtEta, kTH2D, {axisPt, axisEta});
        histos.add(hPtEtaItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPtEta, kTH2D, {axisPt, axisEta});
        histos.add(hPtEtaGenerated[histogramIndex].data(), "Generated " + tagPtEta, kTH2D, {axisPt, axisEta});
      }
    };

    if (doPositivePDG) { // Positive
      makeHistogramsPerCharge(0);
    }
    if (doNegativePDG) { // Negative
      makeHistogramsPerCharge(1);
    }
    LOG(debug) << "Done with particle: " << partName;
  }

  template <o2::track::PID::ID id, int chargeIndex>
  void makeMCEfficiency(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    if (!makeEff) {
      return;
    }

    if constexpr (chargeIndex == 0) {
      if (!doPositivePDG) { // Positive
        return;
      }
    } else if constexpr (chargeIndex == 1) {
      if (!doNegativePDG) { // Negative
        return;
      }
    } else {
      LOG(fatal) << "Can't interpret charge index";
    }

    const char* partName = id == o2::track::PID::NIDs ? "All" : o2::track::PID::getName(id);
    LOG(debug) << "Making TEfficiency for MC for particle " << partName;
    THashList* subList = new THashList();
    subList->SetName(Form("%s_%s", chargeIndex == 0 ? "Positive" : "Negative", partName));
    listEfficiencyMC->Add(subList);

    auto makeEfficiency = [&](TString effname, auto templateHisto) { // 1D efficiencies
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

    const int histogramIndex = id + chargeIndex * nSpecies;

    makeEfficiency("ITS_vsPt", HIST(hPtIts[histogramIndex]));
    makeEfficiency("TPC_vsPt", HIST(hPtTpc[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt", HIST(hPtItsTpc[histogramIndex]));
    makeEfficiency("ITS-TOF_vsPt", HIST(hPtItsTof[histogramIndex]));
    makeEfficiency("Tpc-TOF_vsPt", HIST(hPtTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPt", HIST(hPtItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Trk", HIST(hPtTrkItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Prm", HIST(hPtItsTpcPrm[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Prm_Trk", HIST(hPtTrkItsTpcPrm[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPt_Prm", HIST(hPtItsTpcTofPrm[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Str", HIST(hPtItsTpcStr[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Str_Trk", HIST(hPtTrkItsTpcStr[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPt_Str", HIST(hPtItsTpcTofStr[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Mat", HIST(hPtItsTpcMat[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Mat_Trk", HIST(hPtTrkItsTpcMat[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPt_Mat", HIST(hPtItsTpcTofMat[histogramIndex]));
    makeEfficiency("ITS-TPC_vsP", HIST(hPItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC_vsP_Trk", HIST(hPTrkItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsP", HIST(hPItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsEta", HIST(hEtaItsTpc[histogramIndex]));
    makeEfficiency("ITS-Tpc_vsEta_Trk", HIST(hEtaTrkItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsEta", HIST(hEtaItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsY", HIST(hYItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsY", HIST(hYItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPhi", HIST(hPhiItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPhi_Trk", HIST(hPhiTrkItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPhi", HIST(hPhiItsTpcTof[histogramIndex]));

    auto makeEfficiency2D = [&](TString effname, auto templateHisto) { // 2D efficiencies
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

    if (doPtEta) {
      makeEfficiency2D("ITS-TPC_vsPt_vsEta", HIST(hPtEtaItsTpc[histogramIndex]));
    }
  }

  void initMC(const AxisSpec& axisSel)
  {
    if (!doprocessMC) {
      return;
    }

    auto h = histos.add<TH1>("MC/trackSelection", "Track Selection", kTH1F, {axisSel});
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
    histos.add("MC/fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1F, {{1, 0, 1}});

    h = histos.add<TH1>("MC/particleSelection", "Particle Selection", kTH1F, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Particles read");
    h->GetXaxis()->SetBinLabel(2, "Passed Ev. Reco.");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(6, "Passed y");
    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(7 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/eventMultiplicity", "Event Selection", kTH1F, {{1000, 0, 5000}});

    histos.add("MC/trackLength", "Track length;Track length (cm)", kTH1F, {{2000, -1000, 1000}});

    listEfficiencyMC.setObject(new THashList);
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

    static_for<0, 1>([&](auto i) {
      makeMCEfficiency<o2::track::PID::Electron, i>(doEl);
      makeMCEfficiency<o2::track::PID::Muon, i>(doMu);
      makeMCEfficiency<o2::track::PID::Pion, i>(doPi);
      makeMCEfficiency<o2::track::PID::Kaon, i>(doKa);
      makeMCEfficiency<o2::track::PID::Proton, i>(doPr);
      makeMCEfficiency<o2::track::PID::Deuteron, i>(doDe);
      makeMCEfficiency<o2::track::PID::Triton, i>(doTr);
      makeMCEfficiency<o2::track::PID::Helium3, i>(doHe);
      makeMCEfficiency<o2::track::PID::Alpha, i>(doAl);
      makeMCEfficiency<o2::track::PID::NIDs, i>(doUnId);
    });
  }

  void initData(const AxisSpec& axisSel)
  {
    if (!doprocessData) {
      return;
    }

    auto h = histos.add<TH1>("Data/trackSelection", "Track Selection", kTH1F, {axisSel});
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

    const TString tagPhi = Form("#it{#eta} [%.2f,%.2f] #it{p}_{T} [%.2f,%.2f]",
                                etaMin, etaMax,
                                ptMin, ptMax);
    const AxisSpec axisPhi{phiBins, "#it{#varphi} (rad)"};

    const TString tagEtaPhi = Form("#it{p}_{T} [%.2f,%.2f]",
                                   ptMin, ptMax);

    histos.add("Data/trackLength", "Track length;Track length (cm)", kTH1F, {{2000, -1000, 1000}});

    // ITS-TPC-TOF
    histos.add("Data/pos/pt/its_tpc_tof", "ITS-TPC-TOF Positive " + tagPt, kTH1F, {axisPt});
    histos.add("Data/neg/pt/its_tpc_tof", "ITS-TPC-TOF Negative " + tagPt, kTH1F, {axisPt});

    histos.add("Data/pos/eta/its_tpc_tof", "ITS-TPC-TOF Positive " + tagEta, kTH1F, {axisEta});
    histos.add("Data/neg/eta/its_tpc_tof", "ITS-TPC-TOF Negative " + tagEta, kTH1F, {axisEta});

    histos.add("Data/pos/phi/its_tpc_tof", "ITS-TPC-TOF Positive " + tagPhi, kTH1F, {axisPhi});
    histos.add("Data/neg/phi/its_tpc_tof", "ITS-TPC-TOF Negative " + tagPhi, kTH1F, {axisPhi});

    histos.add("Data/pos/etaphi/its_tpc_tof", "ITS-TPC-TOF Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc_tof", "ITS-TPC-TOF Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // ITS-TPC
    histos.add("Data/pos/pt/its_tpc", "ITS-TPC Positive " + tagPt, kTH1F, {axisPt});
    histos.add("Data/neg/pt/its_tpc", "ITS-TPC Negative " + tagPt, kTH1F, {axisPt});

    histos.add("Data/pos/eta/its_tpc", "ITS-TPC Positive " + tagEta, kTH1F, {axisEta});
    histos.add("Data/neg/eta/its_tpc", "ITS-TPC Negative " + tagEta, kTH1F, {axisEta});

    histos.add("Data/pos/phi/its_tpc", "ITS-TPC Positive " + tagPhi, kTH1F, {axisPhi});
    histos.add("Data/neg/phi/its_tpc", "ITS-TPC Negative " + tagPhi, kTH1F, {axisPhi});

    histos.add("Data/pos/etaphi/its_tpc", "ITS-TPC Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its_tpc", "ITS-TPC Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // TPC
    histos.add("Data/pos/pt/tpc", "TPC Positive " + tagPt, kTH1F, {axisPt});
    histos.add("Data/neg/pt/tpc", "TPC Negative " + tagPt, kTH1F, {axisPt});

    histos.add("Data/pos/eta/tpc", "TPC Positive " + tagEta, kTH1F, {axisEta});
    histos.add("Data/neg/eta/tpc", "TPC Negative " + tagEta, kTH1F, {axisEta});

    histos.add("Data/pos/phi/tpc", "TPC Positive " + tagPhi, kTH1F, {axisPhi});
    histos.add("Data/neg/phi/tpc", "TPC Negative " + tagPhi, kTH1F, {axisPhi});

    histos.add("Data/pos/etaphi/tpc", "TPC Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/tpc", "TPC Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    // ITS
    histos.add("Data/pos/pt/its", "ITS Positive " + tagPt, kTH1F, {axisPt});
    histos.add("Data/neg/pt/its", "ITS Negative " + tagPt, kTH1F, {axisPt});

    histos.add("Data/pos/eta/its", "ITS Positive " + tagEta, kTH1F, {axisEta});
    histos.add("Data/neg/eta/its", "ITS Negative " + tagEta, kTH1F, {axisEta});

    histos.add("Data/pos/phi/its", "ITS Positive " + tagPhi, kTH1F, {axisPhi});
    histos.add("Data/neg/phi/its", "ITS Negative " + tagPhi, kTH1F, {axisPhi});

    histos.add("Data/pos/etaphi/its", "ITS Positive " + tagEtaPhi, kTH2D, {axisEta, axisPhi});
    histos.add("Data/neg/etaphi/its", "ITS Negative " + tagEtaPhi, kTH2D, {axisEta, axisPhi});

    listEfficiencyData.setObject(new THashList);
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

      makeEfficiency("ITSTPCMatchingEfficiencyVsPt", "ITS-TPC M.E. in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsPt", "TPC-TOF M.E. in data " + tagPt + ";#it{p}_{T} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsP", "TPC-TOF M.E. in data " + tagPt + ";#it{p} (GeV/#it{c});Efficiency", HIST("Data/pos/pt/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsEta", "TPC-TOF M.E. in data " + tagEta + ";#it{#eta};Efficiency", HIST("Data/pos/eta/its_tpc_tof"));
      makeEfficiency("TPCTOFMatchingEfficiencyVsPhi", "TPC-TOF M.E. in data " + tagPhi + ";#it{#varphi} (rad);Efficiency", HIST("Data/pos/phi/its_tpc_tof"));

      auto makeEfficiency2D = [&](TString effname, TString efftitle, auto templateHistoX, auto templateHistoY) {
        TAxis* axisX = histos.get<TH1>(templateHistoX)->GetXaxis();
        TAxis* axisY = histos.get<TH1>(templateHistoY)->GetYaxis();
        if (axisX->IsVariableBinSize() || axisY->IsVariableBinSize()) {
          listEfficiencyData->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXbins()->GetArray(), axisY->GetNbins(), axisY->GetXbins()->GetArray()));
        } else {
          listEfficiencyData->Add(new TEfficiency(effname, efftitle, axisX->GetNbins(), axisX->GetXmin(), axisX->GetXmax(), axisY->GetNbins(), axisY->GetXmin(), axisY->GetXmax()));
        }
      };

      makeEfficiency2D("TPCTOFMatchingEfficiencyVsPtVsEta", Form("TPC-TOF M.E. in data #it{#varphi} [%.2f,%.2f];%s;%s;Efficiency", phiMin, phiMax, "#it{p}_{T} (GeV/#it{c})", "#it{#eta}"), HIST("Data/pos/pt/its_tpc_tof"), HIST("Data/pos/eta/its_tpc_tof"));
      makeEfficiency2D("TPCTOFMatchingEfficiencyVsPtVsPhi", Form("TPC-TOF M.E. in data #it{#eta} [%.2f,%.2f];%s;%s;Efficiency", etaMin, etaMax, "#it{p}_{T} (GeV/#it{c})", "#it{#varphi} (rad)"), HIST("Data/pos/pt/its_tpc_tof"), HIST("Data/pos/phi/its_tpc_tof"));
    }
  }

  // Selection cuts defined from the binning
  double ptMin, ptMax;
  double etaMin, etaMax;
  double phiMin, phiMax;
  double yMin, yMax;

  void init(InitContext&)
  {
    auto doLimits = [&](double& min, double& max, const ConfigurableAxis& binning) {
      const AxisSpec a{binning, "dummy"};
      min = a.binEdges[1];
      max = a.binEdges[a.getNbins() - 1];
    };

    doLimits(ptMin, ptMax, ptBins);
    doLimits(etaMin, etaMax, etaBins);
    doLimits(phiMin, phiMax, phiBins);
    doLimits(yMin, yMax, yBins);

    const AxisSpec axisSel{30, 0.5, 30.5, "Selection"};
    histos.add("eventSelection", "Event Selection", kTH1F, {axisSel});
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
    static_assert(pdgSign == 0 || pdgSign == 1);

    // Selecting PDG code
    if constexpr (PDGs[id] == 0) { // All PDGs
      if constexpr (pdgSign == 0) {
        return mcParticle.pdgCode() > 0; // Positive
      } else {
        return mcParticle.pdgCode() < 0; // Negative
      }
    }
    // Specific PDGs
    if constexpr (pdgSign == 0) {
      return mcParticle.pdgCode() == PDGs[id];
    } else {
      return mcParticle.pdgCode() == -PDGs[id];
    }
  }

  template <int pdgSign, o2::track::PID::ID id, typename trackType>
  void fillMCTrackHistograms(const trackType& track)
  {
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

    constexpr int histogramIndex = id + pdgSign * nSpecies;
    LOG(debug) << "Filling track histograms for id " << static_cast<int>(id);
    const auto mcParticle = track.mcParticle();

    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("MC/trackSelection"), 11 + id);

    histos.fill(HIST(hPItsTpc[histogramIndex]), mcParticle.p());
    histos.fill(HIST(hPtItsTpc[histogramIndex]), mcParticle.pt());
    histos.fill(HIST(hEtaItsTpc[histogramIndex]), mcParticle.eta());
    histos.fill(HIST(hYItsTpc[histogramIndex]), mcParticle.y());
    histos.fill(HIST(hPhiItsTpc[histogramIndex]), mcParticle.phi());
    if (doPtEta) {
      histos.fill(HIST(hPtEtaItsTpc[histogramIndex]), mcParticle.pt(), mcParticle.eta());
    }

    histos.fill(HIST(hPTrkItsTpc[histogramIndex]), track.p());
    histos.fill(HIST(hPtTrkItsTpc[histogramIndex]), track.pt());
    histos.fill(HIST(hEtaTrkItsTpc[histogramIndex]), track.eta());
    histos.fill(HIST(hPhiTrkItsTpc[histogramIndex]), track.phi());

    if (mcParticle.isPhysicalPrimary()) {
      histos.fill(HIST(hPtItsTpcPrm[histogramIndex]), mcParticle.pt());
      histos.fill(HIST(hPtTrkItsTpcPrm[histogramIndex]), track.pt());
      if (track.hasTOF()) {
        histos.fill(HIST(hPtItsTpcTofPrm[histogramIndex]), mcParticle.pt());
      }
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        histos.fill(HIST(hPtItsTpcStr[histogramIndex]), mcParticle.pt());
        histos.fill(HIST(hPtTrkItsTpcStr[histogramIndex]), track.pt());
        if (track.hasTOF()) {
          histos.fill(HIST(hPtItsTpcTofStr[histogramIndex]), mcParticle.pt());
        }
      } else { // Material
        histos.fill(HIST(hPtItsTpcMat[histogramIndex]), mcParticle.pt());
        histos.fill(HIST(hPtTrkItsTpcMat[histogramIndex]), track.pt());
        if (track.hasTOF()) {
          histos.fill(HIST(hPtItsTpcTofMat[histogramIndex]), mcParticle.pt());
        }
      }
    }
    if (!track.hasTOF()) {
      return;
    }
    histos.fill(HIST(hPItsTpcTof[histogramIndex]), mcParticle.p());
    histos.fill(HIST(hPtItsTpcTof[histogramIndex]), mcParticle.pt());
    histos.fill(HIST(hEtaItsTpcTof[histogramIndex]), mcParticle.eta());
    histos.fill(HIST(hYItsTpcTof[histogramIndex]), mcParticle.y());
    histos.fill(HIST(hPhiItsTpcTof[histogramIndex]), mcParticle.phi());
  }

  template <int pdgSign, o2::track::PID::ID id, typename particleType>
  void fillMCParticleHistograms(const particleType& mcParticle)
  {
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

    constexpr int histogramIndex = id + pdgSign * nSpecies;
    LOG(debug) << "Filling particle histograms for id " << static_cast<int>(id);
    if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
      return;
    }
    histos.fill(HIST("MC/particleSelection"), 7 + id);

    histos.fill(HIST(hPGenerated[histogramIndex]), mcParticle.p());
    histos.fill(HIST(hPtGenerated[histogramIndex]), mcParticle.pt());

    if (mcParticle.isPhysicalPrimary()) {
      histos.fill(HIST(hPtGeneratedPrm[histogramIndex]), mcParticle.pt());
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        histos.fill(HIST(hPtGeneratedStr[histogramIndex]), mcParticle.pt());
      } else { // Material
        histos.fill(HIST(hPtGeneratedMat[histogramIndex]), mcParticle.pt());
      }
    }

    histos.fill(HIST(hEtaGenerated[histogramIndex]), mcParticle.eta());
    histos.fill(HIST(hYGenerated[histogramIndex]), mcParticle.y());
    histos.fill(HIST(hPhiGenerated[histogramIndex]), mcParticle.phi());
    if (doPtEta) {
      histos.fill(HIST(hPtEtaGenerated[histogramIndex]), mcParticle.pt(), mcParticle.eta());
    }
  }

  template <int pdgSign, o2::track::PID::ID id>
  void fillMCEfficiency()
  {
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

    const char* partName = id == o2::track::PID::NIDs ? "All" : o2::track::PID::getName(id);
    LOG(debug) << "Filling efficiency for particle " << static_cast<int>(id) << " " << partName;
    THashList* subList = static_cast<THashList*>(listEfficiencyMC->FindObject(partName));
    if (!subList) {
      LOG(warning) << "Cannot find list of efficiency objects for particle " << partName;
      return;
    }

    // Filling 1D efficiencies
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

    doFillEfficiency("efficiencyVsPt", HIST(hPtItsTpc[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("efficiencyVsPtPrm", HIST(hPtItsTpcPrm[histogramIndex]), HIST(hPtGeneratedPrm[histogramIndex]));
    doFillEfficiency("efficiencyVsPtDec", HIST(hPtItsTpcStr[histogramIndex]), HIST(hPtGeneratedStr[histogramIndex]));
    doFillEfficiency("efficiencyVsPtMat", HIST(hPtItsTpcMat[histogramIndex]), HIST(hPtGeneratedMat[histogramIndex]));
    doFillEfficiency("efficiencyVsP", HIST(hPItsTpc[histogramIndex]), HIST(hPGenerated[histogramIndex]));
    doFillEfficiency("efficiencyVsEta", HIST(hEtaItsTpc[histogramIndex]), HIST(hEtaGenerated[histogramIndex]));
    doFillEfficiency("efficiencyVsPhi", HIST(hPhiItsTpc[histogramIndex]), HIST(hPhiGenerated[histogramIndex]));

    if (!doPtEta) {
      return;
    }

    // Filling 2D efficiencies
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
    fillEfficiency2D("efficiencyVsPtVsEta", HIST(hPtEtaItsTpc[histogramIndex]), HIST(hPtEtaGenerated[histogramIndex]));
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

  // Function to apply particle selection
  template <bool isMC = true, typename particleType, typename histoType>
  bool isInAcceptance(const particleType& particle, const histoType& countingHisto, const int offset = 0)
  {
    histos.fill(countingHisto, 1 + offset);
    if ((particle.pt() < ptMin || particle.pt() > ptMax)) { // Check pt
      return false;
    }
    histos.fill(countingHisto, 2 + offset);
    if ((particle.eta() < etaMin || particle.eta() > etaMax)) { // Check eta
      return false;
    }
    histos.fill(countingHisto, 3 + offset);
    if ((particle.phi() < phiMin || particle.phi() > phiMax)) { // Check phi
      return false;
    }
    histos.fill(countingHisto, 4 + offset);
    if constexpr (isMC) {
      if ((particle.y() < yMin || particle.y() > yMax)) { // Check rapidity
        return false;
      }
      histos.fill(countingHisto, 5 + offset);
    }

    return true;
  }

  // Function to apply track selection
  bool passedITS = false;
  bool passedTPC = false;
  bool passedTOF = false;
  template <bool isMC = true, typename trackType, typename histoType>
  bool isTrackSelected(trackType& track, const histoType& countingHisto)
  {
    // Reset selections
    passedITS = false;
    passedTPC = false;
    passedTOF = false;

    histos.fill(countingHisto, 1); // Read tracks

    if constexpr (isMC) { // MC only
      if (!track.has_mcParticle()) {
        histos.fill(HIST("MC/fakeTrackNoiseHits"), 0.5);
        return false;
      }
      histos.fill(countingHisto, 2); // Tracks with particles (i.e. no fakes)
      const auto mcParticle = track.mcParticle();
      if (!isInAcceptance(mcParticle, countingHisto, 2)) {
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
      histos.fill(countingHisto, 8);
    } else { // Data only
      if (!isInAcceptance<false>(track, countingHisto, 2)) {
        return false;
      }
    }

    if (!track.has_collision()) {
      return false;
    }
    histos.fill(countingHisto, 9);

    if (trackSelection.value > 0) { // Check general cuts
      if (!track.passedTrackType()) {
        return false;
      }
      histos.fill(countingHisto, 5);
      if (!track.passedPtRange()) {
        return false;
      }
      histos.fill(countingHisto, 6);
      if (!track.passedEtaRange()) {
        return false;
      }
      histos.fill(countingHisto, 7);
      if (!track.passedDCAxy()) {
        return false;
      }
      histos.fill(countingHisto, 8);
      if (!track.passedDCAz()) {
        return false;
      }
      histos.fill(countingHisto, 9);
      if (!track.passedGoldenChi2()) {
        return false;
      }
      histos.fill(countingHisto, 10);

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
      passedTOF = track.hasTOF();
    } else {
      passedITS = track.hasITS();
      passedTPC = track.hasTPC();
      passedTOF = track.hasTOF();
    }

    switch (trackSelection.value) {
      case 0:
        return true;
      case 1:
        return track.isGlobalTrack();
      case 2:
        return track.isGlobalTrackWoPtEta();
      case 3:
        return track.isGlobalTrackWoDCA();
      case 4:
        return track.isQualityTrack();
      case 5:
        return track.isInAcceptanceTrack();
      default:
        LOG(fatal) << "Can't interpret track asked selection";
    }
    return false;
  }

  // MC process
  Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  void processMC(o2::aod::McCollision const& mcCollision,
                 o2::soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>> const& collisions,
                 o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::McTrackLabels, o2::aod::TrackSelection> const& tracks,
                 o2::aod::McParticles const& mcParticles)
  {
    if (collisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
      return;
    }

    for (const auto& collision : collisions) {
      if (!isCollisionSelected<false>(collision)) {
        continue;
      }
      const auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());

      // Track loop
      for (const auto& track : groupedTracks) {
        if (!isTrackSelected(track, HIST("MC/trackSelection"))) {
          continue;
        }
        // Filling variable histograms
        histos.fill(HIST("MC/trackLength"), track.length());
        if (doEl) {
          fillMCTrackHistograms<0, o2::track::PID::Electron>(track);
          fillMCTrackHistograms<1, o2::track::PID::Electron>(track);
        }
        if (doMu) {
          fillMCTrackHistograms<0, o2::track::PID::Muon>(track);
          fillMCTrackHistograms<1, o2::track::PID::Muon>(track);
        }
        if (doPi) {
          fillMCTrackHistograms<0, o2::track::PID::Pion>(track);
          fillMCTrackHistograms<1, o2::track::PID::Pion>(track);
        }
        if (doKa) {
          fillMCTrackHistograms<0, o2::track::PID::Kaon>(track);
          fillMCTrackHistograms<1, o2::track::PID::Kaon>(track);
        }
        if (doPr) {
          fillMCTrackHistograms<0, o2::track::PID::Proton>(track);
          fillMCTrackHistograms<1, o2::track::PID::Proton>(track);
        }
        if (doDe) {
          fillMCTrackHistograms<0, o2::track::PID::Deuteron>(track);
          fillMCTrackHistograms<1, o2::track::PID::Deuteron>(track);
        }
        if (doTr) {
          fillMCTrackHistograms<0, o2::track::PID::Triton>(track);
          fillMCTrackHistograms<1, o2::track::PID::Triton>(track);
        }
        if (doHe) {
          fillMCTrackHistograms<0, o2::track::PID::Helium3>(track);
          fillMCTrackHistograms<1, o2::track::PID::Helium3>(track);
        }
        if (doAl) {
          fillMCTrackHistograms<0, o2::track::PID::Alpha>(track);
          fillMCTrackHistograms<1, o2::track::PID::Alpha>(track);
        }
        if (doUnId) {
          fillMCTrackHistograms<0, o2::track::PID::NIDs>(track);
          fillMCTrackHistograms<1, o2::track::PID::NIDs>(track);
        }
      }
    }

    // Loop on particles to fill the denominator
    float dNdEta = 0; // Multiplicity
    for (const auto& mcParticle : mcParticles) {
      if (TMath::Abs(mcParticle.eta()) <= 2.f && !mcParticle.has_daughters()) {
        dNdEta += 1.f;
      }
      if (!isInAcceptance(mcParticle, HIST("MC/particleSelection"))) {
        continue;
      }

      if (doEl) {
        fillMCParticleHistograms<0, o2::track::PID::Electron>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Electron>(mcParticle);
      }
      if (doMu) {
        fillMCParticleHistograms<0, o2::track::PID::Muon>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Muon>(mcParticle);
      }
      if (doPi) {
        fillMCParticleHistograms<0, o2::track::PID::Pion>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Pion>(mcParticle);
      }
      if (doKa) {
        fillMCParticleHistograms<0, o2::track::PID::Kaon>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Kaon>(mcParticle);
      }
      if (doPr) {
        fillMCParticleHistograms<0, o2::track::PID::Proton>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Proton>(mcParticle);
      }
      if (doDe) {
        fillMCParticleHistograms<0, o2::track::PID::Deuteron>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Deuteron>(mcParticle);
      }
      if (doTr) {
        fillMCParticleHistograms<0, o2::track::PID::Triton>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Triton>(mcParticle);
      }
      if (doHe) {
        fillMCParticleHistograms<0, o2::track::PID::Helium3>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Helium3>(mcParticle);
      }
      if (doAl) {
        fillMCParticleHistograms<0, o2::track::PID::Alpha>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::Alpha>(mcParticle);
      }
      if (doUnId) {
        fillMCParticleHistograms<0, o2::track::PID::NIDs>(mcParticle);
        fillMCParticleHistograms<1, o2::track::PID::NIDs>(mcParticle);
      }
    }
    histos.fill(HIST("MC/eventMultiplicity"), dNdEta * 0.5f / 2.f);

    // Fill TEfficiencies
    if (doEl) {
      fillMCEfficiency<0, o2::track::PID::Electron>();
      fillMCEfficiency<1, o2::track::PID::Electron>();
    }
    if (doMu) {
      fillMCEfficiency<0, o2::track::PID::Muon>();
      fillMCEfficiency<1, o2::track::PID::Muon>();
    }
    if (doPi) {
      fillMCEfficiency<0, o2::track::PID::Pion>();
      fillMCEfficiency<1, o2::track::PID::Pion>();
    }
    if (doKa) {
      fillMCEfficiency<0, o2::track::PID::Kaon>();
      fillMCEfficiency<1, o2::track::PID::Kaon>();
    }
    if (doPr) {
      fillMCEfficiency<0, o2::track::PID::Proton>();
      fillMCEfficiency<1, o2::track::PID::Proton>();
    }
    if (doDe) {
      fillMCEfficiency<0, o2::track::PID::Deuteron>();
      fillMCEfficiency<1, o2::track::PID::Deuteron>();
    }
    if (doTr) {
      fillMCEfficiency<0, o2::track::PID::Triton>();
      fillMCEfficiency<1, o2::track::PID::Triton>();
    }
    if (doHe) {
      fillMCEfficiency<0, o2::track::PID::Helium3>();
      fillMCEfficiency<1, o2::track::PID::Helium3>();
    }
    if (doAl) {
      fillMCEfficiency<0, o2::track::PID::Alpha>();
      fillMCEfficiency<1, o2::track::PID::Alpha>();
    }
    if (doUnId) {
      fillMCEfficiency<0, o2::track::PID::NIDs>();
      fillMCEfficiency<1, o2::track::PID::NIDs>();
    }
  }
  PROCESS_SWITCH(QaEfficiency, processMC, "process MC", false);

  void processData(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
                   const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection>& tracks)
  {

    if (!isCollisionSelected<false>(collision)) {
      return;
    }

    for (const auto& track : tracks) {
      if (!isTrackSelected<false>(track, HIST("Data/trackSelection"))) {
        continue;
      }
      if (trackSelection.value > 0) { // Check general cuts
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
      } else {
        passedITS = track.hasITS();
        passedTPC = track.hasTPC();
      }

      histos.fill(HIST("Data/trackSelection"), 11);

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

        if (track.hasITS() && passedTPC && track.hasTOF()) {
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
      }

      if (makeEff) {
        if (passedITS) {
          static_cast<TEfficiency*>(listEfficiencyData->At(0))->Fill(passedTPC, track.pt());
        }
        if (passedITS && passedTPC) {
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
