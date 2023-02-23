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
  static constexpr int nSpecies = o2::track::PID::NIDs; // One per PDG
  static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040};
  // Track/particle selection
  Configurable<bool> noFakesHits{"noFakesHits", false, "Flag to reject tracks that have fake hits"};
  Configurable<float> maxProdRadius{"maxProdRadius", 9999.f, "Maximum production radius of the particle under study"};
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
  Configurable<int> logPt{"log-pt", 0, "Flag to use a logarithmic pT axis"};
  ConfigurableAxis etaBins{"etaBins", {200, -3.f, 3.f}, "Eta binning"};
  ConfigurableAxis phiBins{"phiBins", {200, 0.f, 6.284f}, "Phi binning"};
  ConfigurableAxis yBins{"yBins", {200, -0.5f, 0.5f}, "Y binning"};
  // Task configuration
  Configurable<bool> makeEff{"make-eff", false, "Flag to produce the efficiency with TEfficiency"};
  Configurable<bool> doPtEta{"doPtEta", false, "Flag to produce the efficiency vs pT and Eta"};
  Configurable<int> applyEvSel{"applyEvSel", 0, "Flag to apply event selection: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  // Custom track cuts for debug purposes
  TrackSelection customTrackCuts;
  Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 1.f, "Additional cut on the maximum value of the DCA xy (multiplicative factor)"};
  Configurable<float> maxDcaZ{"maxDcaZ", 2.f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};

  OutputObj<THashList> listEfficiencyMC{"EfficiencyMC"};
  OutputObj<THashList> listEfficiencyData{"EfficiencyData"};
  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosPos{"HistosPos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosNeg{"HistosNeg", {}, OutputObjHandlingPolicy::AnalysisObject};
  static constexpr int nHistograms = nSpecies * 2;

  // Pt
  static constexpr std::string_view hPtIts[nHistograms] = {"MC/el/pos/pt/its", "MC/mu/pos/pt/its", "MC/pi/pos/pt/its",
                                                           "MC/ka/pos/pt/its", "MC/pr/pos/pt/its", "MC/de/pos/pt/its",
                                                           "MC/tr/pos/pt/its", "MC/he/pos/pt/its", "MC/al/pos/pt/its",
                                                           "MC/el/neg/pt/its", "MC/mu/neg/pt/its", "MC/pi/neg/pt/its",
                                                           "MC/ka/neg/pt/its", "MC/pr/neg/pt/its", "MC/de/neg/pt/its",
                                                           "MC/tr/neg/pt/its", "MC/he/neg/pt/its", "MC/al/neg/pt/its"};
  static constexpr std::string_view hPtTpc[nHistograms] = {"MC/el/pos/pt/tpc", "MC/mu/pos/pt/tpc", "MC/pi/pos/pt/tpc",
                                                           "MC/ka/pos/pt/tpc", "MC/pr/pos/pt/tpc", "MC/de/pos/pt/tpc",
                                                           "MC/tr/pos/pt/tpc", "MC/he/pos/pt/tpc", "MC/al/pos/pt/tpc",
                                                           "MC/el/neg/pt/tpc", "MC/mu/neg/pt/tpc", "MC/pi/neg/pt/tpc",
                                                           "MC/ka/neg/pt/tpc", "MC/pr/neg/pt/tpc", "MC/de/neg/pt/tpc",
                                                           "MC/tr/neg/pt/tpc", "MC/he/neg/pt/tpc", "MC/al/neg/pt/tpc"};
  static constexpr std::string_view hPtItsTpc[nHistograms] = {"MC/el/pos/pt/its_tpc", "MC/mu/pos/pt/its_tpc", "MC/pi/pos/pt/its_tpc",
                                                              "MC/ka/pos/pt/its_tpc", "MC/pr/pos/pt/its_tpc", "MC/de/pos/pt/its_tpc",
                                                              "MC/tr/pos/pt/its_tpc", "MC/he/pos/pt/its_tpc", "MC/al/pos/pt/its_tpc",
                                                              "MC/el/neg/pt/its_tpc", "MC/mu/neg/pt/its_tpc", "MC/pi/neg/pt/its_tpc",
                                                              "MC/ka/neg/pt/its_tpc", "MC/pr/neg/pt/its_tpc", "MC/de/neg/pt/its_tpc",
                                                              "MC/tr/neg/pt/its_tpc", "MC/he/neg/pt/its_tpc", "MC/al/neg/pt/its_tpc"};
  static constexpr std::string_view hPtItsTof[nHistograms] = {"MC/el/pos/pt/its_tof", "MC/mu/pos/pt/its_tof", "MC/pi/pos/pt/its_tof",
                                                              "MC/ka/pos/pt/its_tof", "MC/pr/pos/pt/its_tof", "MC/de/pos/pt/its_tof",
                                                              "MC/tr/pos/pt/its_tof", "MC/he/pos/pt/its_tof", "MC/al/pos/pt/its_tof",
                                                              "MC/el/neg/pt/its_tof", "MC/mu/neg/pt/its_tof", "MC/pi/neg/pt/its_tof",
                                                              "MC/ka/neg/pt/its_tof", "MC/pr/neg/pt/its_tof", "MC/de/neg/pt/its_tof",
                                                              "MC/tr/neg/pt/its_tof", "MC/he/neg/pt/its_tof", "MC/al/neg/pt/its_tof"};
  static constexpr std::string_view hPtTpcTof[nHistograms] = {"MC/el/pos/pt/tpc_tof", "MC/mu/pos/pt/tpc_tof", "MC/pi/pos/pt/tpc_tof",
                                                              "MC/ka/pos/pt/tpc_tof", "MC/pr/pos/pt/tpc_tof", "MC/de/pos/pt/tpc_tof",
                                                              "MC/tr/pos/pt/tpc_tof", "MC/he/pos/pt/tpc_tof", "MC/al/pos/pt/tpc_tof",
                                                              "MC/el/neg/pt/tpc_tof", "MC/mu/neg/pt/tpc_tof", "MC/pi/neg/pt/tpc_tof",
                                                              "MC/ka/neg/pt/tpc_tof", "MC/pr/neg/pt/tpc_tof", "MC/de/neg/pt/tpc_tof",
                                                              "MC/tr/neg/pt/tpc_tof", "MC/he/neg/pt/tpc_tof", "MC/al/neg/pt/tpc_tof"};
  static constexpr std::string_view hPtItsTpcTof[nHistograms] = {"MC/el/pos/pt/its_tpc_tof", "MC/mu/pos/pt/its_tpc_tof", "MC/pi/pos/pt/its_tpc_tof",
                                                                 "MC/ka/pos/pt/its_tpc_tof", "MC/pr/pos/pt/its_tpc_tof", "MC/de/pos/pt/its_tpc_tof",
                                                                 "MC/tr/pos/pt/its_tpc_tof", "MC/he/pos/pt/its_tpc_tof", "MC/al/pos/pt/its_tpc_tof",
                                                                 "MC/el/neg/pt/its_tpc_tof", "MC/mu/neg/pt/its_tpc_tof", "MC/pi/neg/pt/its_tpc_tof",
                                                                 "MC/ka/neg/pt/its_tpc_tof", "MC/pr/neg/pt/its_tpc_tof", "MC/de/neg/pt/its_tpc_tof",
                                                                 "MC/tr/neg/pt/its_tpc_tof", "MC/he/neg/pt/its_tpc_tof", "MC/al/neg/pt/its_tpc_tof"};
  static constexpr std::string_view hPtItsTpcTrdTof[nHistograms] = {"MC/el/pos/pt/its_tpc_trd_tof", "MC/mu/pos/pt/its_tpc_trd_tof", "MC/pi/pos/pt/its_tpc_trd_tof",
                                                                    "MC/ka/pos/pt/its_tpc_trd_tof", "MC/pr/pos/pt/its_tpc_trd_tof", "MC/de/pos/pt/its_tpc_trd_tof",
                                                                    "MC/tr/pos/pt/its_tpc_trd_tof", "MC/he/pos/pt/its_tpc_trd_tof", "MC/al/pos/pt/its_tpc_trd_tof",
                                                                    "MC/el/neg/pt/its_tpc_trd_tof", "MC/mu/neg/pt/its_tpc_trd_tof", "MC/pi/neg/pt/its_tpc_trd_tof",
                                                                    "MC/ka/neg/pt/its_tpc_trd_tof", "MC/pr/neg/pt/its_tpc_trd_tof", "MC/de/neg/pt/its_tpc_trd_tof",
                                                                    "MC/tr/neg/pt/its_tpc_trd_tof", "MC/he/neg/pt/its_tpc_trd_tof", "MC/al/neg/pt/its_tpc_trd_tof"};
  static constexpr std::string_view hPtItsTpcTrd[nHistograms] = {"MC/el/pos/pt/its_tpc_trd", "MC/mu/pos/pt/its_tpc_trd", "MC/pi/pos/pt/its_tpc_trd",
                                                                 "MC/ka/pos/pt/its_tpc_trd", "MC/pr/pos/pt/its_tpc_trd", "MC/de/pos/pt/its_tpc_trd",
                                                                 "MC/tr/pos/pt/its_tpc_trd", "MC/he/pos/pt/its_tpc_trd", "MC/al/pos/pt/its_tpc_trd",
                                                                 "MC/el/neg/pt/its_tpc_trd", "MC/mu/neg/pt/its_tpc_trd", "MC/pi/neg/pt/its_tpc_trd",
                                                                 "MC/ka/neg/pt/its_tpc_trd", "MC/pr/neg/pt/its_tpc_trd", "MC/de/neg/pt/its_tpc_trd",
                                                                 "MC/tr/neg/pt/its_tpc_trd", "MC/he/neg/pt/its_tpc_trd", "MC/al/neg/pt/its_tpc_trd"};
  static constexpr std::string_view hPtTrkItsTpc[nHistograms] = {"MC/el/pos/pt/trk/its_tpc", "MC/mu/pos/pt/trk/its_tpc", "MC/pi/pos/pt/trk/its_tpc",
                                                                 "MC/ka/pos/pt/trk/its_tpc", "MC/pr/pos/pt/trk/its_tpc", "MC/de/pos/pt/trk/its_tpc",
                                                                 "MC/tr/pos/pt/trk/its_tpc", "MC/he/pos/pt/trk/its_tpc", "MC/al/pos/pt/trk/its_tpc",
                                                                 "MC/el/neg/pt/trk/its_tpc", "MC/mu/neg/pt/trk/its_tpc", "MC/pi/neg/pt/trk/its_tpc",
                                                                 "MC/ka/neg/pt/trk/its_tpc", "MC/pr/neg/pt/trk/its_tpc", "MC/de/neg/pt/trk/its_tpc",
                                                                 "MC/tr/neg/pt/trk/its_tpc", "MC/he/neg/pt/trk/its_tpc", "MC/al/neg/pt/trk/its_tpc"};
  static constexpr std::string_view hPtGenerated[nHistograms] = {"MC/el/pos/pt/generated", "MC/mu/pos/pt/generated", "MC/pi/pos/pt/generated",
                                                                 "MC/ka/pos/pt/generated", "MC/pr/pos/pt/generated", "MC/de/pos/pt/generated",
                                                                 "MC/tr/pos/pt/generated", "MC/he/pos/pt/generated", "MC/al/pos/pt/generated",
                                                                 "MC/el/neg/pt/generated", "MC/mu/neg/pt/generated", "MC/pi/neg/pt/generated",
                                                                 "MC/ka/neg/pt/generated", "MC/pr/neg/pt/generated", "MC/de/neg/pt/generated",
                                                                 "MC/tr/neg/pt/generated", "MC/he/neg/pt/generated", "MC/al/neg/pt/generated"};

  // Pt for primaries
  static constexpr std::string_view hPtItsPrm[nHistograms] = {"MC/el/pos/pt/prm/its", "MC/mu/pos/pt/prm/its", "MC/pi/pos/pt/prm/its",
                                                              "MC/ka/pos/pt/prm/its", "MC/pr/pos/pt/prm/its", "MC/de/pos/pt/prm/its",
                                                              "MC/tr/pos/pt/prm/its", "MC/he/pos/pt/prm/its", "MC/al/pos/pt/prm/its",
                                                              "MC/el/neg/pt/prm/its", "MC/mu/neg/pt/prm/its", "MC/pi/neg/pt/prm/its",
                                                              "MC/ka/neg/pt/prm/its", "MC/pr/neg/pt/prm/its", "MC/de/neg/pt/prm/its",
                                                              "MC/tr/neg/pt/prm/its", "MC/he/neg/pt/prm/its", "MC/al/neg/pt/prm/its"};
  static constexpr std::string_view hPtItsTpcPrm[nHistograms] = {"MC/el/pos/pt/prm/its_tpc", "MC/mu/pos/pt/prm/its_tpc", "MC/pi/pos/pt/prm/its_tpc",
                                                                 "MC/ka/pos/pt/prm/its_tpc", "MC/pr/pos/pt/prm/its_tpc", "MC/de/pos/pt/prm/its_tpc",
                                                                 "MC/tr/pos/pt/prm/its_tpc", "MC/he/pos/pt/prm/its_tpc", "MC/al/pos/pt/prm/its_tpc",
                                                                 "MC/el/neg/pt/prm/its_tpc", "MC/mu/neg/pt/prm/its_tpc", "MC/pi/neg/pt/prm/its_tpc",
                                                                 "MC/ka/neg/pt/prm/its_tpc", "MC/pr/neg/pt/prm/its_tpc", "MC/de/neg/pt/prm/its_tpc",
                                                                 "MC/tr/neg/pt/prm/its_tpc", "MC/he/neg/pt/prm/its_tpc", "MC/al/neg/pt/prm/its_tpc"};
  static constexpr std::string_view hPtTrkItsTpcPrm[nHistograms] = {"MC/el/pos/pt/prm/trk/its_tpc", "MC/mu/pos/pt/prm/trk/its_tpc", "MC/pi/pos/pt/prm/trk/its_tpc",
                                                                    "MC/ka/pos/pt/prm/trk/its_tpc", "MC/pr/pos/pt/prm/trk/its_tpc", "MC/de/pos/pt/prm/trk/its_tpc",
                                                                    "MC/tr/pos/pt/prm/trk/its_tpc", "MC/he/pos/pt/prm/trk/its_tpc", "MC/al/pos/pt/prm/trk/its_tpc",
                                                                    "MC/el/neg/pt/prm/trk/its_tpc", "MC/mu/neg/pt/prm/trk/its_tpc", "MC/pi/neg/pt/prm/trk/its_tpc",
                                                                    "MC/ka/neg/pt/prm/trk/its_tpc", "MC/pr/neg/pt/prm/trk/its_tpc", "MC/de/neg/pt/prm/trk/its_tpc",
                                                                    "MC/tr/neg/pt/prm/trk/its_tpc", "MC/he/neg/pt/prm/trk/its_tpc", "MC/al/neg/pt/prm/trk/its_tpc"};
  static constexpr std::string_view hPtItsTpcTofPrm[nHistograms] = {"MC/el/pos/pt/prm/its_tpc_tof", "MC/mu/pos/pt/prm/its_tpc_tof", "MC/pi/pos/pt/prm/its_tpc_tof",
                                                                    "MC/ka/pos/pt/prm/its_tpc_tof", "MC/pr/pos/pt/prm/its_tpc_tof", "MC/de/pos/pt/prm/its_tpc_tof",
                                                                    "MC/tr/pos/pt/prm/its_tpc_tof", "MC/he/pos/pt/prm/its_tpc_tof", "MC/al/pos/pt/prm/its_tpc_tof",
                                                                    "MC/el/neg/pt/prm/its_tpc_tof", "MC/mu/neg/pt/prm/its_tpc_tof", "MC/pi/neg/pt/prm/its_tpc_tof",
                                                                    "MC/ka/neg/pt/prm/its_tpc_tof", "MC/pr/neg/pt/prm/its_tpc_tof", "MC/de/neg/pt/prm/its_tpc_tof",
                                                                    "MC/tr/neg/pt/prm/its_tpc_tof", "MC/he/neg/pt/prm/its_tpc_tof", "MC/al/neg/pt/prm/its_tpc_tof"};
  static constexpr std::string_view hPtGeneratedPrm[nHistograms] = {"MC/el/pos/pt/prm/generated", "MC/mu/pos/pt/prm/generated", "MC/pi/pos/pt/prm/generated",
                                                                    "MC/ka/pos/pt/prm/generated", "MC/pr/pos/pt/prm/generated", "MC/de/pos/pt/prm/generated",
                                                                    "MC/tr/pos/pt/prm/generated", "MC/he/pos/pt/prm/generated", "MC/al/pos/pt/prm/generated",
                                                                    "MC/el/neg/pt/prm/generated", "MC/mu/neg/pt/prm/generated", "MC/pi/neg/pt/prm/generated",
                                                                    "MC/ka/neg/pt/prm/generated", "MC/pr/neg/pt/prm/generated", "MC/de/neg/pt/prm/generated",
                                                                    "MC/tr/neg/pt/prm/generated", "MC/he/neg/pt/prm/generated", "MC/al/neg/pt/prm/generated"};

  // Pt for secondaries from weak decay
  static constexpr std::string_view hPtItsTpcStr[nHistograms] = {"MC/el/pos/pt/str/its_tpc", "MC/mu/pos/pt/str/its_tpc", "MC/pi/pos/pt/str/its_tpc",
                                                                 "MC/ka/pos/pt/str/its_tpc", "MC/pr/pos/pt/str/its_tpc", "MC/de/pos/pt/str/its_tpc",
                                                                 "MC/tr/pos/pt/str/its_tpc", "MC/he/pos/pt/str/its_tpc", "MC/al/pos/pt/str/its_tpc",
                                                                 "MC/el/neg/pt/str/its_tpc", "MC/mu/neg/pt/str/its_tpc", "MC/pi/neg/pt/str/its_tpc",
                                                                 "MC/ka/neg/pt/str/its_tpc", "MC/pr/neg/pt/str/its_tpc", "MC/de/neg/pt/str/its_tpc",
                                                                 "MC/tr/neg/pt/str/its_tpc", "MC/he/neg/pt/str/its_tpc", "MC/al/neg/pt/str/its_tpc"};
  static constexpr std::string_view hPtTrkItsTpcStr[nHistograms] = {"MC/el/pos/pt/str/trk/its_tpc", "MC/mu/pos/pt/str/trk/its_tpc", "MC/pi/pos/pt/str/trk/its_tpc",
                                                                    "MC/ka/pos/pt/str/trk/its_tpc", "MC/pr/pos/pt/str/trk/its_tpc", "MC/de/pos/pt/str/trk/its_tpc",
                                                                    "MC/tr/pos/pt/str/trk/its_tpc", "MC/he/pos/pt/str/trk/its_tpc", "MC/al/pos/pt/str/trk/its_tpc",
                                                                    "MC/el/neg/pt/str/trk/its_tpc", "MC/mu/neg/pt/str/trk/its_tpc", "MC/pi/neg/pt/str/trk/its_tpc",
                                                                    "MC/ka/neg/pt/str/trk/its_tpc", "MC/pr/neg/pt/str/trk/its_tpc", "MC/de/neg/pt/str/trk/its_tpc",
                                                                    "MC/tr/neg/pt/str/trk/its_tpc", "MC/he/neg/pt/str/trk/its_tpc", "MC/al/neg/pt/str/trk/its_tpc"};
  static constexpr std::string_view hPtItsTpcTofStr[nHistograms] = {"MC/el/pos/pt/str/its_tpc_tof", "MC/mu/pos/pt/str/its_tpc_tof", "MC/pi/pos/pt/str/its_tpc_tof",
                                                                    "MC/ka/pos/pt/str/its_tpc_tof", "MC/pr/pos/pt/str/its_tpc_tof", "MC/de/pos/pt/str/its_tpc_tof",
                                                                    "MC/tr/pos/pt/str/its_tpc_tof", "MC/he/pos/pt/str/its_tpc_tof", "MC/al/pos/pt/str/its_tpc_tof",
                                                                    "MC/el/neg/pt/str/its_tpc_tof", "MC/mu/neg/pt/str/its_tpc_tof", "MC/pi/neg/pt/str/its_tpc_tof",
                                                                    "MC/ka/neg/pt/str/its_tpc_tof", "MC/pr/neg/pt/str/its_tpc_tof", "MC/de/neg/pt/str/its_tpc_tof",
                                                                    "MC/tr/neg/pt/str/its_tpc_tof", "MC/he/neg/pt/str/its_tpc_tof", "MC/al/neg/pt/str/its_tpc_tof"};
  static constexpr std::string_view hPtGeneratedStr[nHistograms] = {"MC/el/pos/pt/str/generated", "MC/mu/pos/pt/str/generated", "MC/pi/pos/pt/str/generated",
                                                                    "MC/ka/pos/pt/str/generated", "MC/pr/pos/pt/str/generated", "MC/de/pos/pt/str/generated",
                                                                    "MC/tr/pos/pt/str/generated", "MC/he/pos/pt/str/generated", "MC/al/pos/pt/str/generated",
                                                                    "MC/el/neg/pt/str/generated", "MC/mu/neg/pt/str/generated", "MC/pi/neg/pt/str/generated",
                                                                    "MC/ka/neg/pt/str/generated", "MC/pr/neg/pt/str/generated", "MC/de/neg/pt/str/generated",
                                                                    "MC/tr/neg/pt/str/generated", "MC/he/neg/pt/str/generated", "MC/al/neg/pt/str/generated"};

  // Pt for secondaries from material
  static constexpr std::string_view hPtItsTpcMat[nHistograms] = {"MC/el/pos/pt/mat/its_tpc", "MC/mu/pos/pt/mat/its_tpc", "MC/pi/pos/pt/mat/its_tpc",
                                                                 "MC/ka/pos/pt/mat/its_tpc", "MC/pr/pos/pt/mat/its_tpc", "MC/de/pos/pt/mat/its_tpc",
                                                                 "MC/tr/pos/pt/mat/its_tpc", "MC/he/pos/pt/mat/its_tpc", "MC/al/pos/pt/mat/its_tpc",
                                                                 "MC/el/neg/pt/mat/its_tpc", "MC/mu/neg/pt/mat/its_tpc", "MC/pi/neg/pt/mat/its_tpc",
                                                                 "MC/ka/neg/pt/mat/its_tpc", "MC/pr/neg/pt/mat/its_tpc", "MC/de/neg/pt/mat/its_tpc",
                                                                 "MC/tr/neg/pt/mat/its_tpc", "MC/he/neg/pt/mat/its_tpc", "MC/al/neg/pt/mat/its_tpc"};
  static constexpr std::string_view hPtTrkItsTpcMat[nHistograms] = {"MC/el/pos/pt/mat/trk/its_tpc", "MC/mu/pos/pt/mat/trk/its_tpc", "MC/pi/pos/pt/mat/trk/its_tpc",
                                                                    "MC/ka/pos/pt/mat/trk/its_tpc", "MC/pr/pos/pt/mat/trk/its_tpc", "MC/de/pos/pt/mat/trk/its_tpc",
                                                                    "MC/tr/pos/pt/mat/trk/its_tpc", "MC/he/pos/pt/mat/trk/its_tpc", "MC/al/pos/pt/mat/trk/its_tpc",
                                                                    "MC/el/neg/pt/mat/trk/its_tpc", "MC/mu/neg/pt/mat/trk/its_tpc", "MC/pi/neg/pt/mat/trk/its_tpc",
                                                                    "MC/ka/neg/pt/mat/trk/its_tpc", "MC/pr/neg/pt/mat/trk/its_tpc", "MC/de/neg/pt/mat/trk/its_tpc",
                                                                    "MC/tr/neg/pt/mat/trk/its_tpc", "MC/he/neg/pt/mat/trk/its_tpc", "MC/al/neg/pt/mat/trk/its_tpc"};
  static constexpr std::string_view hPtItsTpcTofMat[nHistograms] = {"MC/el/pos/pt/mat/its_tpc_tof", "MC/mu/pos/pt/mat/its_tpc_tof", "MC/pi/pos/pt/mat/its_tpc_tof",
                                                                    "MC/ka/pos/pt/mat/its_tpc_tof", "MC/pr/pos/pt/mat/its_tpc_tof", "MC/de/pos/pt/mat/its_tpc_tof",
                                                                    "MC/tr/pos/pt/mat/its_tpc_tof", "MC/he/pos/pt/mat/its_tpc_tof", "MC/al/pos/pt/mat/its_tpc_tof",
                                                                    "MC/el/neg/pt/mat/its_tpc_tof", "MC/mu/neg/pt/mat/its_tpc_tof", "MC/pi/neg/pt/mat/its_tpc_tof",
                                                                    "MC/ka/neg/pt/mat/its_tpc_tof", "MC/pr/neg/pt/mat/its_tpc_tof", "MC/de/neg/pt/mat/its_tpc_tof",
                                                                    "MC/tr/neg/pt/mat/its_tpc_tof", "MC/he/neg/pt/mat/its_tpc_tof", "MC/al/neg/pt/mat/its_tpc_tof"};
  static constexpr std::string_view hPtGeneratedMat[nHistograms] = {"MC/el/pos/pt/mat/generated", "MC/mu/pos/pt/mat/generated", "MC/pi/pos/pt/mat/generated",
                                                                    "MC/ka/pos/pt/mat/generated", "MC/pr/pos/pt/mat/generated", "MC/de/pos/pt/mat/generated",
                                                                    "MC/tr/pos/pt/mat/generated", "MC/he/pos/pt/mat/generated", "MC/al/pos/pt/mat/generated",
                                                                    "MC/el/neg/pt/mat/generated", "MC/mu/neg/pt/mat/generated", "MC/pi/neg/pt/mat/generated",
                                                                    "MC/ka/neg/pt/mat/generated", "MC/pr/neg/pt/mat/generated", "MC/de/neg/pt/mat/generated",
                                                                    "MC/tr/neg/pt/mat/generated", "MC/he/neg/pt/mat/generated", "MC/al/neg/pt/mat/generated"};

  // P
  static constexpr std::string_view hPItsTpc[nHistograms] = {"MC/el/pos/p/its_tpc", "MC/mu/pos/p/its_tpc", "MC/pi/pos/p/its_tpc",
                                                             "MC/ka/pos/p/its_tpc", "MC/pr/pos/p/its_tpc", "MC/de/pos/p/its_tpc",
                                                             "MC/tr/pos/p/its_tpc", "MC/he/pos/p/its_tpc", "MC/al/pos/p/its_tpc",
                                                             "MC/el/neg/p/its_tpc", "MC/mu/neg/p/its_tpc", "MC/pi/neg/p/its_tpc",
                                                             "MC/ka/neg/p/its_tpc", "MC/pr/neg/p/its_tpc", "MC/de/neg/p/its_tpc",
                                                             "MC/tr/neg/p/its_tpc", "MC/he/neg/p/its_tpc", "MC/al/neg/p/its_tpc"};
  static constexpr std::string_view hPTrkItsTpc[nHistograms] = {"MC/el/pos/p/trk/its_tpc", "MC/mu/pos/p/trk/its_tpc", "MC/pi/pos/p/trk/its_tpc",
                                                                "MC/ka/pos/p/trk/its_tpc", "MC/pr/pos/p/trk/its_tpc", "MC/de/pos/p/trk/its_tpc",
                                                                "MC/tr/pos/p/trk/its_tpc", "MC/he/pos/p/trk/its_tpc", "MC/al/pos/p/trk/its_tpc",
                                                                "MC/el/neg/p/trk/its_tpc", "MC/mu/neg/p/trk/its_tpc", "MC/pi/neg/p/trk/its_tpc",
                                                                "MC/ka/neg/p/trk/its_tpc", "MC/pr/neg/p/trk/its_tpc", "MC/de/neg/p/trk/its_tpc",
                                                                "MC/tr/neg/p/trk/its_tpc", "MC/he/neg/p/trk/its_tpc", "MC/al/neg/p/trk/its_tpc"};
  static constexpr std::string_view hPItsTpcTof[nHistograms] = {"MC/el/pos/p/its_tpc_tof", "MC/mu/pos/p/its_tpc_tof", "MC/pi/pos/p/its_tpc_tof",
                                                                "MC/ka/pos/p/its_tpc_tof", "MC/pr/pos/p/its_tpc_tof", "MC/de/pos/p/its_tpc_tof",
                                                                "MC/tr/pos/p/its_tpc_tof", "MC/he/pos/p/its_tpc_tof", "MC/al/pos/p/its_tpc_tof",
                                                                "MC/el/neg/p/its_tpc_tof", "MC/mu/neg/p/its_tpc_tof", "MC/pi/neg/p/its_tpc_tof",
                                                                "MC/ka/neg/p/its_tpc_tof", "MC/pr/neg/p/its_tpc_tof", "MC/de/neg/p/its_tpc_tof",
                                                                "MC/tr/neg/p/its_tpc_tof", "MC/he/neg/p/its_tpc_tof", "MC/al/neg/p/its_tpc_tof"};
  static constexpr std::string_view hPGenerated[nHistograms] = {"MC/el/pos/p/generated", "MC/mu/pos/p/generated", "MC/pi/pos/p/generated",
                                                                "MC/ka/pos/p/generated", "MC/pr/pos/p/generated", "MC/de/pos/p/generated",
                                                                "MC/tr/pos/p/generated", "MC/he/pos/p/generated", "MC/al/pos/p/generated",
                                                                "MC/el/neg/p/generated", "MC/mu/neg/p/generated", "MC/pi/neg/p/generated",
                                                                "MC/ka/neg/p/generated", "MC/pr/neg/p/generated", "MC/de/neg/p/generated",
                                                                "MC/tr/neg/p/generated", "MC/he/neg/p/generated", "MC/al/neg/p/generated"};

  // Eta
  static constexpr std::string_view hEtaItsTpc[nHistograms] = {"MC/el/pos/eta/its_tpc", "MC/mu/pos/eta/its_tpc", "MC/pi/pos/eta/its_tpc",
                                                               "MC/ka/pos/eta/its_tpc", "MC/pr/pos/eta/its_tpc", "MC/de/pos/eta/its_tpc",
                                                               "MC/tr/pos/eta/its_tpc", "MC/he/pos/eta/its_tpc", "MC/al/pos/eta/its_tpc",
                                                               "MC/el/neg/eta/its_tpc", "MC/mu/neg/eta/its_tpc", "MC/pi/neg/eta/its_tpc",
                                                               "MC/ka/neg/eta/its_tpc", "MC/pr/neg/eta/its_tpc", "MC/de/neg/eta/its_tpc",
                                                               "MC/tr/neg/eta/its_tpc", "MC/he/neg/eta/its_tpc", "MC/al/neg/eta/its_tpc"};
  static constexpr std::string_view hEtaTrkItsTpc[nHistograms] = {"MC/el/pos/eta/trk/its_tpc", "MC/mu/pos/eta/trk/its_tpc", "MC/pi/pos/eta/trk/its_tpc",
                                                                  "MC/ka/pos/eta/trk/its_tpc", "MC/pr/pos/eta/trk/its_tpc", "MC/de/pos/eta/trk/its_tpc",
                                                                  "MC/tr/pos/eta/trk/its_tpc", "MC/he/pos/eta/trk/its_tpc", "MC/al/pos/eta/trk/its_tpc",
                                                                  "MC/el/neg/eta/trk/its_tpc", "MC/mu/neg/eta/trk/its_tpc", "MC/pi/neg/eta/trk/its_tpc",
                                                                  "MC/ka/neg/eta/trk/its_tpc", "MC/pr/neg/eta/trk/its_tpc", "MC/de/neg/eta/trk/its_tpc",
                                                                  "MC/tr/neg/eta/trk/its_tpc", "MC/he/neg/eta/trk/its_tpc", "MC/al/neg/eta/trk/its_tpc"};
  static constexpr std::string_view hEtaItsTpcTof[nHistograms] = {"MC/el/pos/eta/its_tpc_tof", "MC/mu/pos/eta/its_tpc_tof", "MC/pi/pos/eta/its_tpc_tof",
                                                                  "MC/ka/pos/eta/its_tpc_tof", "MC/pr/pos/eta/its_tpc_tof", "MC/de/pos/eta/its_tpc_tof",
                                                                  "MC/tr/pos/eta/its_tpc_tof", "MC/he/pos/eta/its_tpc_tof", "MC/al/pos/eta/its_tpc_tof",
                                                                  "MC/el/neg/eta/its_tpc_tof", "MC/mu/neg/eta/its_tpc_tof", "MC/pi/neg/eta/its_tpc_tof",
                                                                  "MC/ka/neg/eta/its_tpc_tof", "MC/pr/neg/eta/its_tpc_tof", "MC/de/neg/eta/its_tpc_tof",
                                                                  "MC/tr/neg/eta/its_tpc_tof", "MC/he/neg/eta/its_tpc_tof", "MC/al/neg/eta/its_tpc_tof"};
  static constexpr std::string_view hEtaGenerated[nHistograms] = {"MC/el/pos/eta/generated", "MC/mu/pos/eta/generated", "MC/pi/pos/eta/generated",
                                                                  "MC/ka/pos/eta/generated", "MC/pr/pos/eta/generated", "MC/de/pos/eta/generated",
                                                                  "MC/tr/pos/eta/generated", "MC/he/pos/eta/generated", "MC/al/pos/eta/generated",
                                                                  "MC/el/neg/eta/generated", "MC/mu/neg/eta/generated", "MC/pi/neg/eta/generated",
                                                                  "MC/ka/neg/eta/generated", "MC/pr/neg/eta/generated", "MC/de/neg/eta/generated",
                                                                  "MC/tr/neg/eta/generated", "MC/he/neg/eta/generated", "MC/al/neg/eta/generated"};

  // Eta for primaries
  static constexpr std::string_view hEtaItsTpcPrm[nHistograms] = {"MC/el/pos/eta/prm/its_tpc", "MC/mu/pos/eta/prm/its_tpc", "MC/pi/pos/eta/prm/its_tpc",
                                                                  "MC/ka/pos/eta/prm/its_tpc", "MC/pr/pos/eta/prm/its_tpc", "MC/de/pos/eta/prm/its_tpc",
                                                                  "MC/tr/pos/eta/prm/its_tpc", "MC/he/pos/eta/prm/its_tpc", "MC/al/pos/eta/prm/its_tpc",
                                                                  "MC/el/neg/eta/prm/its_tpc", "MC/mu/neg/eta/prm/its_tpc", "MC/pi/neg/eta/prm/its_tpc",
                                                                  "MC/ka/neg/eta/prm/its_tpc", "MC/pr/neg/eta/prm/its_tpc", "MC/de/neg/eta/prm/its_tpc",
                                                                  "MC/tr/neg/eta/prm/its_tpc", "MC/he/neg/eta/prm/its_tpc", "MC/al/neg/eta/prm/its_tpc"};
  static constexpr std::string_view hEtaTrkItsTpcPrm[nHistograms] = {"MC/el/pos/eta/prm/trk/its_tpc", "MC/mu/pos/eta/prm/trk/its_tpc", "MC/pi/pos/eta/prm/trk/its_tpc",
                                                                     "MC/ka/pos/eta/prm/trk/its_tpc", "MC/pr/pos/eta/prm/trk/its_tpc", "MC/de/pos/eta/prm/trk/its_tpc",
                                                                     "MC/tr/pos/eta/prm/trk/its_tpc", "MC/he/pos/eta/prm/trk/its_tpc", "MC/al/pos/eta/prm/trk/its_tpc",
                                                                     "MC/el/neg/eta/prm/trk/its_tpc", "MC/mu/neg/eta/prm/trk/its_tpc", "MC/pi/neg/eta/prm/trk/its_tpc",
                                                                     "MC/ka/neg/eta/prm/trk/its_tpc", "MC/pr/neg/eta/prm/trk/its_tpc", "MC/de/neg/eta/prm/trk/its_tpc",
                                                                     "MC/tr/neg/eta/prm/trk/its_tpc", "MC/he/neg/eta/prm/trk/its_tpc", "MC/al/neg/eta/prm/trk/its_tpc"};
  static constexpr std::string_view hEtaItsTpcTofPrm[nHistograms] = {"MC/el/pos/eta/prm/its_tpc_tof", "MC/mu/pos/eta/prm/its_tpc_tof", "MC/pi/pos/eta/prm/its_tpc_tof",
                                                                     "MC/ka/pos/eta/prm/its_tpc_tof", "MC/pr/pos/eta/prm/its_tpc_tof", "MC/de/pos/eta/prm/its_tpc_tof",
                                                                     "MC/tr/pos/eta/prm/its_tpc_tof", "MC/he/pos/eta/prm/its_tpc_tof", "MC/al/pos/eta/prm/its_tpc_tof",
                                                                     "MC/el/neg/eta/prm/its_tpc_tof", "MC/mu/neg/eta/prm/its_tpc_tof", "MC/pi/neg/eta/prm/its_tpc_tof",
                                                                     "MC/ka/neg/eta/prm/its_tpc_tof", "MC/pr/neg/eta/prm/its_tpc_tof", "MC/de/neg/eta/prm/its_tpc_tof",
                                                                     "MC/tr/neg/eta/prm/its_tpc_tof", "MC/he/neg/eta/prm/its_tpc_tof", "MC/al/neg/eta/prm/its_tpc_tof"};
  static constexpr std::string_view hEtaGeneratedPrm[nHistograms] = {"MC/el/pos/eta/prm/generated", "MC/mu/pos/eta/prm/generated", "MC/pi/pos/eta/prm/generated",
                                                                     "MC/ka/pos/eta/prm/generated", "MC/pr/pos/eta/prm/generated", "MC/de/pos/eta/prm/generated",
                                                                     "MC/tr/pos/eta/prm/generated", "MC/he/pos/eta/prm/generated", "MC/al/pos/eta/prm/generated",
                                                                     "MC/el/neg/eta/prm/generated", "MC/mu/neg/eta/prm/generated", "MC/pi/neg/eta/prm/generated",
                                                                     "MC/ka/neg/eta/prm/generated", "MC/pr/neg/eta/prm/generated", "MC/de/neg/eta/prm/generated",
                                                                     "MC/tr/neg/eta/prm/generated", "MC/he/neg/eta/prm/generated", "MC/al/neg/eta/prm/generated"};

  // Y
  static constexpr std::string_view hYItsTpc[nHistograms] = {"MC/el/pos/y/its_tpc", "MC/mu/pos/y/its_tpc", "MC/pi/pos/y/its_tpc",
                                                             "MC/ka/pos/y/its_tpc", "MC/pr/pos/y/its_tpc", "MC/de/pos/y/its_tpc",
                                                             "MC/tr/pos/y/its_tpc", "MC/he/pos/y/its_tpc", "MC/al/pos/y/its_tpc",
                                                             "MC/el/neg/y/its_tpc", "MC/mu/neg/y/its_tpc", "MC/pi/neg/y/its_tpc",
                                                             "MC/ka/neg/y/its_tpc", "MC/pr/neg/y/its_tpc", "MC/de/neg/y/its_tpc",
                                                             "MC/tr/neg/y/its_tpc", "MC/he/neg/y/its_tpc", "MC/al/neg/y/its_tpc"};
  static constexpr std::string_view hYItsTpcTof[nHistograms] = {"MC/el/pos/y/its_tpc_tof", "MC/mu/pos/y/its_tpc_tof", "MC/pi/pos/y/its_tpc_tof",
                                                                "MC/ka/pos/y/its_tpc_tof", "MC/pr/pos/y/its_tpc_tof", "MC/de/pos/y/its_tpc_tof",
                                                                "MC/tr/pos/y/its_tpc_tof", "MC/he/pos/y/its_tpc_tof", "MC/al/pos/y/its_tpc_tof",
                                                                "MC/el/neg/y/its_tpc_tof", "MC/mu/neg/y/its_tpc_tof", "MC/pi/neg/y/its_tpc_tof",
                                                                "MC/ka/neg/y/its_tpc_tof", "MC/pr/neg/y/its_tpc_tof", "MC/de/neg/y/its_tpc_tof",
                                                                "MC/tr/neg/y/its_tpc_tof", "MC/he/neg/y/its_tpc_tof", "MC/al/neg/y/its_tpc_tof"};
  static constexpr std::string_view hYGenerated[nHistograms] = {"MC/el/pos/y/generated", "MC/mu/pos/y/generated", "MC/pi/pos/y/generated",
                                                                "MC/ka/pos/y/generated", "MC/pr/pos/y/generated", "MC/de/pos/y/generated",
                                                                "MC/tr/pos/y/generated", "MC/he/pos/y/generated", "MC/al/pos/y/generated",
                                                                "MC/el/neg/y/generated", "MC/mu/neg/y/generated", "MC/pi/neg/y/generated",
                                                                "MC/ka/neg/y/generated", "MC/pr/neg/y/generated", "MC/de/neg/y/generated",
                                                                "MC/tr/neg/y/generated", "MC/he/neg/y/generated", "MC/al/neg/y/generated"};

  // Phi
  static constexpr std::string_view hPhiItsTpc[nHistograms] = {"MC/el/pos/phi/its_tpc", "MC/mu/pos/phi/its_tpc", "MC/pi/pos/phi/its_tpc",
                                                               "MC/ka/pos/phi/its_tpc", "MC/pr/pos/phi/its_tpc", "MC/de/pos/phi/its_tpc",
                                                               "MC/tr/pos/phi/its_tpc", "MC/he/pos/phi/its_tpc", "MC/al/pos/phi/its_tpc",
                                                               "MC/el/neg/phi/its_tpc", "MC/mu/neg/phi/its_tpc", "MC/pi/neg/phi/its_tpc",
                                                               "MC/ka/neg/phi/its_tpc", "MC/pr/neg/phi/its_tpc", "MC/de/neg/phi/its_tpc",
                                                               "MC/tr/neg/phi/its_tpc", "MC/he/neg/phi/its_tpc", "MC/al/neg/phi/its_tpc"};
  static constexpr std::string_view hPhiTrkItsTpc[nHistograms] = {"MC/el/pos/phi/trk/its_tpc", "MC/mu/pos/phi/trk/its_tpc", "MC/pi/pos/phi/trk/its_tpc",
                                                                  "MC/ka/pos/phi/trk/its_tpc", "MC/pr/pos/phi/trk/its_tpc", "MC/de/pos/phi/trk/its_tpc",
                                                                  "MC/tr/pos/phi/trk/its_tpc", "MC/he/pos/phi/trk/its_tpc", "MC/al/pos/phi/trk/its_tpc",
                                                                  "MC/el/neg/phi/trk/its_tpc", "MC/mu/neg/phi/trk/its_tpc", "MC/pi/neg/phi/trk/its_tpc",
                                                                  "MC/ka/neg/phi/trk/its_tpc", "MC/pr/neg/phi/trk/its_tpc", "MC/de/neg/phi/trk/its_tpc",
                                                                  "MC/tr/neg/phi/trk/its_tpc", "MC/he/neg/phi/trk/its_tpc", "MC/al/neg/phi/trk/its_tpc"};
  static constexpr std::string_view hPhiItsTpcTof[nHistograms] = {"MC/el/pos/phi/its_tpc_tof", "MC/mu/pos/phi/its_tpc_tof", "MC/pi/pos/phi/its_tpc_tof",
                                                                  "MC/ka/pos/phi/its_tpc_tof", "MC/pr/pos/phi/its_tpc_tof", "MC/de/pos/phi/its_tpc_tof",
                                                                  "MC/tr/pos/phi/its_tpc_tof", "MC/he/pos/phi/its_tpc_tof", "MC/al/pos/phi/its_tpc_tof",
                                                                  "MC/el/neg/phi/its_tpc_tof", "MC/mu/neg/phi/its_tpc_tof", "MC/pi/neg/phi/its_tpc_tof",
                                                                  "MC/ka/neg/phi/its_tpc_tof", "MC/pr/neg/phi/its_tpc_tof", "MC/de/neg/phi/its_tpc_tof",
                                                                  "MC/tr/neg/phi/its_tpc_tof", "MC/he/neg/phi/its_tpc_tof", "MC/al/neg/phi/its_tpc_tof"};
  static constexpr std::string_view hPhiGenerated[nHistograms] = {"MC/el/pos/phi/generated", "MC/mu/pos/phi/generated", "MC/pi/pos/phi/generated",
                                                                  "MC/ka/pos/phi/generated", "MC/pr/pos/phi/generated", "MC/de/pos/phi/generated",
                                                                  "MC/tr/pos/phi/generated", "MC/he/pos/phi/generated", "MC/al/pos/phi/generated",
                                                                  "MC/el/neg/phi/generated", "MC/mu/neg/phi/generated", "MC/pi/neg/phi/generated",
                                                                  "MC/ka/neg/phi/generated", "MC/pr/neg/phi/generated", "MC/de/neg/phi/generated",
                                                                  "MC/tr/neg/phi/generated", "MC/he/neg/phi/generated", "MC/al/neg/phi/generated"};

  // Phi for primaries
  static constexpr std::string_view hPhiItsTpcPrm[nHistograms] = {"MC/el/pos/phi/prm/its_tpc", "MC/mu/pos/phi/prm/its_tpc", "MC/pi/pos/phi/prm/its_tpc",
                                                                  "MC/ka/pos/phi/prm/its_tpc", "MC/pr/pos/phi/prm/its_tpc", "MC/de/pos/phi/prm/its_tpc",
                                                                  "MC/tr/pos/phi/prm/its_tpc", "MC/he/pos/phi/prm/its_tpc", "MC/al/pos/phi/prm/its_tpc",
                                                                  "MC/el/neg/phi/prm/its_tpc", "MC/mu/neg/phi/prm/its_tpc", "MC/pi/neg/phi/prm/its_tpc",
                                                                  "MC/ka/neg/phi/prm/its_tpc", "MC/pr/neg/phi/prm/its_tpc", "MC/de/neg/phi/prm/its_tpc",
                                                                  "MC/tr/neg/phi/prm/its_tpc", "MC/he/neg/phi/prm/its_tpc", "MC/al/neg/phi/prm/its_tpc"};
  static constexpr std::string_view hPhiTrkItsTpcPrm[nHistograms] = {"MC/el/pos/phi/prm/trk/its_tpc", "MC/mu/pos/phi/prm/trk/its_tpc", "MC/pi/pos/phi/prm/trk/its_tpc",
                                                                     "MC/ka/pos/phi/prm/trk/its_tpc", "MC/pr/pos/phi/prm/trk/its_tpc", "MC/de/pos/phi/prm/trk/its_tpc",
                                                                     "MC/tr/pos/phi/prm/trk/its_tpc", "MC/he/pos/phi/prm/trk/its_tpc", "MC/al/pos/phi/prm/trk/its_tpc",
                                                                     "MC/el/neg/phi/prm/trk/its_tpc", "MC/mu/neg/phi/prm/trk/its_tpc", "MC/pi/neg/phi/prm/trk/its_tpc",
                                                                     "MC/ka/neg/phi/prm/trk/its_tpc", "MC/pr/neg/phi/prm/trk/its_tpc", "MC/de/neg/phi/prm/trk/its_tpc",
                                                                     "MC/tr/neg/phi/prm/trk/its_tpc", "MC/he/neg/phi/prm/trk/its_tpc", "MC/al/neg/phi/prm/trk/its_tpc"};
  static constexpr std::string_view hPhiItsTpcTofPrm[nHistograms] = {"MC/el/pos/phi/prm/its_tpc_tof", "MC/mu/pos/phi/prm/its_tpc_tof", "MC/pi/pos/phi/prm/its_tpc_tof",
                                                                     "MC/ka/pos/phi/prm/its_tpc_tof", "MC/pr/pos/phi/prm/its_tpc_tof", "MC/de/pos/phi/prm/its_tpc_tof",
                                                                     "MC/tr/pos/phi/prm/its_tpc_tof", "MC/he/pos/phi/prm/its_tpc_tof", "MC/al/pos/phi/prm/its_tpc_tof",
                                                                     "MC/el/neg/phi/prm/its_tpc_tof", "MC/mu/neg/phi/prm/its_tpc_tof", "MC/pi/neg/phi/prm/its_tpc_tof",
                                                                     "MC/ka/neg/phi/prm/its_tpc_tof", "MC/pr/neg/phi/prm/its_tpc_tof", "MC/de/neg/phi/prm/its_tpc_tof",
                                                                     "MC/tr/neg/phi/prm/its_tpc_tof", "MC/he/neg/phi/prm/its_tpc_tof", "MC/al/neg/phi/prm/its_tpc_tof"};
  static constexpr std::string_view hPhiGeneratedPrm[nHistograms] = {"MC/el/pos/phi/prm/generated", "MC/mu/pos/phi/prm/generated", "MC/pi/pos/phi/prm/generated",
                                                                     "MC/ka/pos/phi/prm/generated", "MC/pr/pos/phi/prm/generated", "MC/de/pos/phi/prm/generated",
                                                                     "MC/tr/pos/phi/prm/generated", "MC/he/pos/phi/prm/generated", "MC/al/pos/phi/prm/generated",
                                                                     "MC/el/neg/phi/prm/generated", "MC/mu/neg/phi/prm/generated", "MC/pi/neg/phi/prm/generated",
                                                                     "MC/ka/neg/phi/prm/generated", "MC/pr/neg/phi/prm/generated", "MC/de/neg/phi/prm/generated",
                                                                     "MC/tr/neg/phi/prm/generated", "MC/he/neg/phi/prm/generated", "MC/al/neg/phi/prm/generated"};

  // Pt-Eta
  static constexpr std::string_view hPtEtaItsTpc[nHistograms] = {"MC/el/pos/pteta/its_tpc", "MC/mu/pos/pteta/its_tpc", "MC/pi/pos/pteta/its_tpc",
                                                                 "MC/ka/pos/pteta/its_tpc", "MC/pr/pos/pteta/its_tpc", "MC/de/pos/pteta/its_tpc",
                                                                 "MC/tr/pos/pteta/its_tpc", "MC/he/pos/pteta/its_tpc", "MC/al/pos/pteta/its_tpc",
                                                                 "MC/el/neg/pteta/its_tpc", "MC/mu/neg/pteta/its_tpc", "MC/pi/neg/pteta/its_tpc",
                                                                 "MC/ka/neg/pteta/its_tpc", "MC/pr/neg/pteta/its_tpc", "MC/de/neg/pteta/its_tpc",
                                                                 "MC/tr/neg/pteta/its_tpc", "MC/he/neg/pteta/its_tpc", "MC/al/neg/pteta/its_tpc"};
  static constexpr std::string_view hPtEtaTrkItsTpc[nHistograms] = {"MC/el/pos/pteta/trk/its_tpc", "MC/mu/pos/pteta/trk/its_tpc", "MC/pi/pos/pteta/trk/its_tpc",
                                                                    "MC/ka/pos/pteta/trk/its_tpc", "MC/pr/pos/pteta/trk/its_tpc", "MC/de/pos/pteta/trk/its_tpc",
                                                                    "MC/tr/pos/pteta/trk/its_tpc", "MC/he/pos/pteta/trk/its_tpc", "MC/al/pos/pteta/trk/its_tpc",
                                                                    "MC/el/neg/pteta/trk/its_tpc", "MC/mu/neg/pteta/trk/its_tpc", "MC/pi/neg/pteta/trk/its_tpc",
                                                                    "MC/ka/neg/pteta/trk/its_tpc", "MC/pr/neg/pteta/trk/its_tpc", "MC/de/neg/pteta/trk/its_tpc",
                                                                    "MC/tr/neg/pteta/trk/its_tpc", "MC/he/neg/pteta/trk/its_tpc", "MC/al/neg/pteta/trk/its_tpc"};
  static constexpr std::string_view hPtEtaItsTpcTof[nHistograms] = {"MC/el/pos/pteta/its_tpc_tof", "MC/mu/pos/pteta/its_tpc_tof", "MC/pi/pos/pteta/its_tpc_tof",
                                                                    "MC/ka/pos/pteta/its_tpc_tof", "MC/pr/pos/pteta/its_tpc_tof", "MC/de/pos/pteta/its_tpc_tof",
                                                                    "MC/tr/pos/pteta/its_tpc_tof", "MC/he/pos/pteta/its_tpc_tof", "MC/al/pos/pteta/its_tpc_tof",
                                                                    "MC/el/neg/pteta/its_tpc_tof", "MC/mu/neg/pteta/its_tpc_tof", "MC/pi/neg/pteta/its_tpc_tof",
                                                                    "MC/ka/neg/pteta/its_tpc_tof", "MC/pr/neg/pteta/its_tpc_tof", "MC/de/neg/pteta/its_tpc_tof",
                                                                    "MC/tr/neg/pteta/its_tpc_tof", "MC/he/neg/pteta/its_tpc_tof", "MC/al/neg/pteta/its_tpc_tof"};
  static constexpr std::string_view hPtEtaGenerated[nHistograms] = {"MC/el/pos/pteta/generated", "MC/mu/pos/pteta/generated", "MC/pi/pos/pteta/generated",
                                                                    "MC/ka/pos/pteta/generated", "MC/pr/pos/pteta/generated", "MC/de/pos/pteta/generated",
                                                                    "MC/tr/pos/pteta/generated", "MC/he/pos/pteta/generated", "MC/al/pos/pteta/generated",
                                                                    "MC/el/neg/pteta/generated", "MC/mu/neg/pteta/generated", "MC/pi/neg/pteta/generated",
                                                                    "MC/ka/neg/pteta/generated", "MC/pr/neg/pteta/generated", "MC/de/neg/pteta/generated",
                                                                    "MC/tr/neg/pteta/generated", "MC/he/neg/pteta/generated", "MC/al/neg/pteta/generated"};

  static const char* particleName(int charge, o2::track::PID::ID id)
  {
    return Form("%s %s", charge == 0 ? "Positive" : "Negative", o2::track::PID::getName(id));
  }

  template <int charge, o2::track::PID::ID id>
  void makeMCHistograms(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    if constexpr (charge == 0) {
      if (!doPositivePDG) { // Positive
        return;
      }
    } else if constexpr (charge == 1) {
      if (!doNegativePDG) { // Negative
        return;
      }
    } else {
      LOG(fatal) << "Can't interpret charge " << charge;
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

    const char* partName = particleName(charge, id);
    LOG(info) << "Preparing histograms for particle: " << partName << " charge " << charge;

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

    const int histogramIndex = id + charge * nSpecies;
    HistogramRegistry* registry = &histosPos;
    if (charge == 1) {
      registry = &histosNeg;
    }

    registry->add(hPtIts[histogramIndex].data(), "ITS tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtTpc[histogramIndex].data(), "TPC tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTof[histogramIndex].data(), "ITS-TOF tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtTpcTof[histogramIndex].data(), "TPC-TOF tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcTrd[histogramIndex].data(), "ITS-TPC-TRD tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcTrdTof[histogramIndex].data(), "ITS-TPC-TRD-TOF tracks " + tagPt, kTH1F, {axisPt});
    registry->add(hPtTrkItsTpc[histogramIndex].data(), "ITS-TPC track (reco) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtGenerated[histogramIndex].data(), "Generated " + tagPt, kTH1F, {axisPt});

    registry->add(hPtItsPrm[histogramIndex].data(), "ITS tracks (primaries) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (primaries) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtTrkItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (reco primaries) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcTofPrm[histogramIndex].data(), "ITS-TPC-TOF tracks (primaries) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtGeneratedPrm[histogramIndex].data(), "Generated (primaries) " + tagPt, kTH1F, {axisPt});

    registry->add(hPtItsTpcStr[histogramIndex].data(), "ITS-TPC tracks (from weak decays) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtTrkItsTpcStr[histogramIndex].data(), "ITS-TPC tracks (reco from weak decays) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcTofStr[histogramIndex].data(), "ITS-TPC-TOF tracks (from weak decays) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtGeneratedStr[histogramIndex].data(), "Generated (from weak decays) " + tagPt, kTH1F, {axisPt});

    registry->add(hPtItsTpcMat[histogramIndex].data(), "ITS-TPC tracks (from material)" + tagPt, kTH1F, {axisPt});
    registry->add(hPtTrkItsTpcMat[histogramIndex].data(), "ITS-TPC tracks (reco from material) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtItsTpcTofMat[histogramIndex].data(), "ITS-TPC-TOF tracks ( from material) " + tagPt, kTH1F, {axisPt});
    registry->add(hPtGeneratedMat[histogramIndex].data(), "Generated ( from material) " + tagPt, kTH1F, {axisPt});

    registry->add(hPItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPt, kTH1F, {axisP});
    registry->add(hPTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagPt, kTH1F, {axisP});
    registry->add(hPItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPt, kTH1F, {axisP});
    registry->add(hPGenerated[histogramIndex].data(), "Generated " + tagPt, kTH1F, {axisP});

    registry->add(hEtaItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagEta, kTH1F, {axisEta});
    registry->add(hEtaTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagEta, kTH1F, {axisEta});
    registry->add(hEtaItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagEta, kTH1F, {axisEta});
    registry->add(hEtaGenerated[histogramIndex].data(), "Generated " + tagEta, kTH1F, {axisEta});

    registry->add(hEtaItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (primaries) " + tagEta, kTH1F, {axisEta});
    registry->add(hEtaTrkItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (reco primaries) " + tagEta, kTH1F, {axisEta});
    registry->add(hEtaItsTpcTofPrm[histogramIndex].data(), "ITS-TPC-TOF tracks (primaries) " + tagEta, kTH1F, {axisEta});
    registry->add(hEtaGeneratedPrm[histogramIndex].data(), "Generated (primaries) " + tagEta, kTH1F, {axisEta});

    registry->add(hYItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagY, kTH1F, {axisY});
    registry->add(hYItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagY, kTH1F, {axisY});
    registry->add(hYGenerated[histogramIndex].data(), "Generated " + tagY, kTH1F, {axisY});

    registry->add(hPhiItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPhi, kTH1F, {axisPhi});
    registry->add(hPhiTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagPhi, kTH1F, {axisPhi});
    registry->add(hPhiItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPhi, kTH1F, {axisPhi});
    registry->add(hPhiGenerated[histogramIndex].data(), "Generated " + tagPhi, kTH1F, {axisPhi});

    registry->add(hPhiItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (primaries) " + tagPhi, kTH1F, {axisPhi});
    registry->add(hPhiTrkItsTpcPrm[histogramIndex].data(), "ITS-TPC tracks (reco primaries) " + tagPhi, kTH1F, {axisPhi});
    registry->add(hPhiItsTpcTofPrm[histogramIndex].data(), "ITS-TPC-TOF tracks (primaries) " + tagPhi, kTH1F, {axisPhi});
    registry->add(hPhiGeneratedPrm[histogramIndex].data(), "Generated (primaries) " + tagPhi, kTH1F, {axisPhi});

    if (doPtEta) {
      registry->add(hPtEtaItsTpc[histogramIndex].data(), "ITS-TPC tracks " + tagPtEta, kTH2D, {axisPt, axisEta});
      registry->add(hPtEtaTrkItsTpc[histogramIndex].data(), "ITS-TPC tracks (reco) " + tagPtEta, kTH2D, {axisPt, axisEta});
      registry->add(hPtEtaItsTpcTof[histogramIndex].data(), "ITS-TPC-TOF tracks " + tagPtEta, kTH2D, {axisPt, axisEta});
      registry->add(hPtEtaGenerated[histogramIndex].data(), "Generated " + tagPtEta, kTH2D, {axisPt, axisEta});
    }

    LOG(info) << "Done with particle: " << partName;
  }

  template <int charge, o2::track::PID::ID id>
  void makeMCEfficiency(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    if (!makeEff) {
      return;
    }

    if constexpr (charge == 0) {
      if (!doPositivePDG) { // Positive
        return;
      }
    } else if constexpr (charge == 1) {
      if (!doNegativePDG) { // Negative
        return;
      }
    } else {
      LOG(fatal) << "Can't interpret charge index";
    }

    const TString partName = particleName(charge, id);
    LOG(info) << "Making TEfficiency for MC for particle " << partName;
    THashList* subList = new THashList();
    subList->SetName(partName);
    listEfficiencyMC->Add(subList);

    HistogramRegistry* registry = &histosPos;
    if (charge == 1) {
      registry = &histosNeg;
    }

    auto makeEfficiency = [&](const TString effname, auto templateHisto) { // 1D efficiencies
      const auto h = registry->get<TH1>(templateHisto);
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

    const int histogramIndex = id + charge * nSpecies;

    makeEfficiency("ITS_vsPt", HIST(hPtIts[histogramIndex]));
    makeEfficiency("TPC_vsPt", HIST(hPtTpc[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt", HIST(hPtItsTpc[histogramIndex]));
    makeEfficiency("ITS-TOF_vsPt", HIST(hPtItsTof[histogramIndex]));
    makeEfficiency("Tpc-TOF_vsPt", HIST(hPtTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC-TRD_vsPt", HIST(hPtItsTpcTrd[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPt", HIST(hPtItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC-TRD-TOF_vsPt", HIST(hPtItsTpcTrdTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPt_Trk", HIST(hPtTrkItsTpc[histogramIndex]));

    makeEfficiency("ITS_vsPt_Prm", HIST(hPtItsPrm[histogramIndex]));
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
    makeEfficiency("ITS-TPC_vsEta_Trk", HIST(hEtaTrkItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsEta", HIST(hEtaItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsEta_Prm", HIST(hEtaItsTpcPrm[histogramIndex]));
    makeEfficiency("ITS-TPC_vsEta_Prm_Trk", HIST(hEtaTrkItsTpcPrm[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsEta_Prm", HIST(hEtaItsTpcTofPrm[histogramIndex]));
    makeEfficiency("ITS-TPC_vsY", HIST(hYItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsY", HIST(hYItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPhi", HIST(hPhiItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPhi_Trk", HIST(hPhiTrkItsTpc[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPhi", HIST(hPhiItsTpcTof[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPhi_Prm", HIST(hPhiItsTpcPrm[histogramIndex]));
    makeEfficiency("ITS-TPC_vsPhi_Prm_Trk", HIST(hPhiTrkItsTpcPrm[histogramIndex]));
    makeEfficiency("ITS-TPC-TOF_vsPhi_Prm", HIST(hPhiItsTpcTofPrm[histogramIndex]));

    auto makeEfficiency2D = [&](const TString effname, auto templateHisto) { // 2D efficiencies
      const auto h = registry->get<TH2>(templateHisto);
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
      makeEfficiency2D("ITS-TPC_vsPt_vsEta", HIST(hPtEtaItsTpc[histogramIndex]));
      makeEfficiency2D("ITS-TPC_vsPt_vsEta_Trk", HIST(hPtEtaTrkItsTpc[histogramIndex]));
      makeEfficiency2D("ITS-TPC-TOF_vsPt_vsEta", HIST(hPtEtaItsTpcTof[histogramIndex]));
    }

    LOG(info) << "Done with particle: " << partName << " for efficiencies";
  }

  void initMC(const AxisSpec& axisSel)
  {
    if (!doprocessMC && !doprocessMCWithoutCollisions) {
      return;
    }
    if (doprocessMC && doprocessMCWithoutCollisions) {
      LOG(fatal) << "Both processMC and processMCWithoutCollisions are set to true. Please set only one of them to true.";
    }

    auto h = histos.add<TH1>("MC/trackSelection", "Track Selection", kTH1F, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Passed has MC part.");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(6, "Passed y");
    h->GetXaxis()->SetBinLabel(7, "Passed Fake");
    h->GetXaxis()->SetBinLabel(8, "Passed has collision");
    h->GetXaxis()->SetBinLabel(9, "passedTrackType");
    h->GetXaxis()->SetBinLabel(10, "passedPtRange");
    h->GetXaxis()->SetBinLabel(11, "passedEtaRange");
    h->GetXaxis()->SetBinLabel(12, "passedDCAxy");
    h->GetXaxis()->SetBinLabel(13, "passedDCAz");
    h->GetXaxis()->SetBinLabel(14, "passedGoldenChi2");
    h->GetXaxis()->SetBinLabel(15, "passedITS (partial)");
    h->GetXaxis()->SetBinLabel(16, "passedTPC (partial)");
    h->GetXaxis()->SetBinLabel(17, "passedTOF (partial)");
    switch (globalTrackSelection) {
      case 0:
        h->GetXaxis()->SetBinLabel(18, "No extra selection");
        break;
      case 1:
        h->GetXaxis()->SetBinLabel(18, "isGlobalTrack");
        break;
      case 2:
        h->GetXaxis()->SetBinLabel(18, "isGlobalTrackWoPtEta");
        break;
      case 3:
        h->GetXaxis()->SetBinLabel(18, "isGlobalTrackWoDCA");
        break;
      case 4:
        h->GetXaxis()->SetBinLabel(18, "isQualityTrack");
        break;
      case 5:
        h->GetXaxis()->SetBinLabel(18, "isInAcceptanceTrack");
        break;
      case 6:
        h->GetXaxis()->SetBinLabel(18, "customTrackSelection");
        break;
      default:
        LOG(fatal) << "Can't interpret track asked selection";
    }

    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(19 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/fakeTrackNoiseHits", "Fake tracks from noise hits", kTH1F, {{1, 0, 1}});

    h = histos.add<TH1>("MC/particleSelection", "Particle Selection", kTH1F, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Particles read");
    h->GetXaxis()->SetBinLabel(2, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(5, "Passed y");
    for (int i = 0; i < nSpecies; i++) {
      h->GetXaxis()->SetBinLabel(6 + i, Form("Passed PDG %i %s", PDGs[i], particleTitle[i]));
    }
    histos.add("MC/eventMultiplicity", "Event Selection", kTH1F, {{1000, 0, 5000}});

    histos.add("MC/trackLength", "Track length;Track length (cm)", kTH1F, {{2000, -1000, 1000}});

    listEfficiencyMC.setObject(new THashList);

    static_for<0, 1>([&](auto charge) {
      makeMCHistograms<charge, o2::track::PID::Electron>(doEl);
      makeMCHistograms<charge, o2::track::PID::Muon>(doMu);
      makeMCHistograms<charge, o2::track::PID::Pion>(doPi);
      makeMCHistograms<charge, o2::track::PID::Kaon>(doKa);
      makeMCHistograms<charge, o2::track::PID::Proton>(doPr);
      makeMCHistograms<charge, o2::track::PID::Deuteron>(doDe);
      makeMCHistograms<charge, o2::track::PID::Triton>(doTr);
      makeMCHistograms<charge, o2::track::PID::Helium3>(doHe);
      makeMCHistograms<charge, o2::track::PID::Alpha>(doAl);

      makeMCEfficiency<charge, o2::track::PID::Electron>(doEl);
      makeMCEfficiency<charge, o2::track::PID::Muon>(doMu);
      makeMCEfficiency<charge, o2::track::PID::Pion>(doPi);
      makeMCEfficiency<charge, o2::track::PID::Kaon>(doKa);
      makeMCEfficiency<charge, o2::track::PID::Proton>(doPr);
      makeMCEfficiency<charge, o2::track::PID::Deuteron>(doDe);
      makeMCEfficiency<charge, o2::track::PID::Triton>(doTr);
      makeMCEfficiency<charge, o2::track::PID::Helium3>(doHe);
      makeMCEfficiency<charge, o2::track::PID::Alpha>(doAl);
    });
  }

  void initData(const AxisSpec& axisSel)
  {
    if (!doprocessData) {
      return;
    }

    auto h = histos.add<TH1>("Data/trackSelection", "Track Selection", kTH1F, {axisSel});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "");
    h->GetXaxis()->SetBinLabel(3, "Passed #it{p}_{T}");
    h->GetXaxis()->SetBinLabel(4, "Passed #it{#eta}");
    h->GetXaxis()->SetBinLabel(5, "Passed #it{#varphi}");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, "");
    h->GetXaxis()->SetBinLabel(8, "Passed has collision");
    h->GetXaxis()->SetBinLabel(9, "passedTrackType");
    h->GetXaxis()->SetBinLabel(10, "passedPtRange");
    h->GetXaxis()->SetBinLabel(11, "passedEtaRange");
    h->GetXaxis()->SetBinLabel(12, "passedDCAxy");
    h->GetXaxis()->SetBinLabel(13, "passedDCAz");
    h->GetXaxis()->SetBinLabel(14, "passedGoldenChi2");
    h->GetXaxis()->SetBinLabel(15, "passedITS (partial)");
    h->GetXaxis()->SetBinLabel(16, "passedTPC (partial)");
    h->GetXaxis()->SetBinLabel(17, "passedTOF (partial)");
    h->GetXaxis()->SetBinLabel(18, "Passed globalCut");

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

    const AxisSpec axisSel{40, 0.5, 40.5, "Selection"};
    histos.add("eventSelection", "Event Selection", kTH1F, {axisSel});
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(1, "Events read");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(2, "Passed Ev. Sel.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(3, "Passed Contrib.");
    histos.get<TH1>(HIST("eventSelection"))->GetXaxis()->SetBinLabel(4, "Passed Position");

    initData(axisSel);
    initMC(axisSel);

    // Custom track cuts
    if (globalTrackSelection.value == 6) {
      customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
      LOG(info) << "Customizing track cuts:";
      customTrackCuts.SetRequireITSRefit(requireITS.value);
      customTrackCuts.SetRequireTPCRefit(requireTPC.value);
      customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
      customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
      customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
      customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
      customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
      customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
      customTrackCuts.SetMaxDcaXYPtDep([](float pt) { return 10000.f; }); // No DCAxy cut will be used, this is done via the member function of the task
      customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
      customTrackCuts.print();
    }
  }

  template <int charge, o2::track::PID::ID id, typename particleType>
  bool isPdgSelected(const particleType& mcParticle)
  {
    static_assert(charge == 0 || charge == 1);
    static_assert(id > 0 || id < nSpecies);

    // Selecting a specific PDG
    if constexpr (charge == 0) {
      return mcParticle.pdgCode() == PDGs[id];
    } else {
      return mcParticle.pdgCode() == -PDGs[id];
    }
  }

  template <int charge, o2::track::PID::ID id, typename trackType>
  void fillMCTrackHistograms(const trackType& track, const bool doMakeHistograms)
  {
    static_assert(charge == 0 || charge == 1);
    if (!doMakeHistograms) {
      return;
    }

    if constexpr (charge == 0) {
      if (!doPositivePDG) {
        return;
      }
    } else {
      if (!doNegativePDG) {
        return;
      }
    }

    HistogramRegistry* h = &histosPos;
    if constexpr (charge == 1) {
      h = &histosNeg;
    }

    constexpr int histogramIndex = id + charge * nSpecies;
    LOG(debug) << "fillMCTrackHistograms for charge '" << charge << "' and id '" << static_cast<int>(id) << "' " << particleName(charge, id) << " with index " << histogramIndex;
    const auto mcParticle = track.mcParticle();

    if (!isPdgSelected<charge, id>(mcParticle)) { // Selecting PDG code
      return;
    }

    histos.fill(HIST("MC/trackSelection"), 19 + id);

    if (passedITS) {
      h->fill(HIST(hPtIts[histogramIndex]), mcParticle.pt());
    }
    if (passedTPC) {
      h->fill(HIST(hPtTpc[histogramIndex]), mcParticle.pt());
    }
    if (passedITS && passedTPC) {
      h->fill(HIST(hPItsTpc[histogramIndex]), mcParticle.p());
      h->fill(HIST(hPtItsTpc[histogramIndex]), mcParticle.pt());
      h->fill(HIST(hEtaItsTpc[histogramIndex]), mcParticle.eta());
      h->fill(HIST(hYItsTpc[histogramIndex]), mcParticle.y());
      h->fill(HIST(hPhiItsTpc[histogramIndex]), mcParticle.phi());

      h->fill(HIST(hPTrkItsTpc[histogramIndex]), track.p());
      h->fill(HIST(hPtTrkItsTpc[histogramIndex]), track.pt());
      h->fill(HIST(hEtaTrkItsTpc[histogramIndex]), track.eta());
      h->fill(HIST(hPhiTrkItsTpc[histogramIndex]), track.phi());

      if (doPtEta) {
        h->fill(HIST(hPtEtaItsTpc[histogramIndex]), mcParticle.pt(), mcParticle.eta());
        h->fill(HIST(hPtEtaTrkItsTpc[histogramIndex]), track.pt(), track.eta());
        if (passedTOF) {
          h->fill(HIST(hPtEtaItsTpcTof[histogramIndex]), mcParticle.pt(), mcParticle.eta());
        }
      }
    }
    if (passedITS && passedTOF) {
      h->fill(HIST(hPtItsTof[histogramIndex]), mcParticle.pt());
    }
    if (passedTPC && passedTOF) {
      h->fill(HIST(hPtTpcTof[histogramIndex]), mcParticle.pt());
    }
    if (passedITS && passedTPC && passedTRD) {
      h->fill(HIST(hPtItsTpcTrd[histogramIndex]), mcParticle.p());
    }
    if (passedITS && passedTPC && passedTOF) {
      h->fill(HIST(hPItsTpcTof[histogramIndex]), mcParticle.p());
      h->fill(HIST(hPtItsTpcTof[histogramIndex]), mcParticle.pt());
      h->fill(HIST(hEtaItsTpcTof[histogramIndex]), mcParticle.eta());
      h->fill(HIST(hYItsTpcTof[histogramIndex]), mcParticle.y());
      h->fill(HIST(hPhiItsTpcTof[histogramIndex]), mcParticle.phi());
    }
    if (passedITS && passedTPC && passedTRD && passedTOF) {
      h->fill(HIST(hPtItsTpcTrdTof[histogramIndex]), mcParticle.p());
    }

    bool isPhysicalPrimary = mcParticle.isPhysicalPrimary();
    if (maxProdRadius < 999.f) {
      if ((mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()) > maxProdRadius * maxProdRadius) {
        isPhysicalPrimary = false;
      }
    }

    if (isPhysicalPrimary) {
      if (passedITS) {
        h->fill(HIST(hPtItsPrm[histogramIndex]), mcParticle.pt());
      }
      if (passedITS && passedTPC) {
        h->fill(HIST(hPtItsTpcPrm[histogramIndex]), mcParticle.pt());
        h->fill(HIST(hPtTrkItsTpcPrm[histogramIndex]), track.pt());
        h->fill(HIST(hEtaItsTpcPrm[histogramIndex]), mcParticle.eta());
        h->fill(HIST(hEtaTrkItsTpcPrm[histogramIndex]), track.eta());
        h->fill(HIST(hPhiItsTpcPrm[histogramIndex]), mcParticle.phi());
        h->fill(HIST(hPhiTrkItsTpcPrm[histogramIndex]), track.phi());
        if (passedTOF) {
          h->fill(HIST(hPtItsTpcTofPrm[histogramIndex]), mcParticle.pt());
          h->fill(HIST(hEtaItsTpcTofPrm[histogramIndex]), mcParticle.eta());
          h->fill(HIST(hPhiItsTpcTofPrm[histogramIndex]), mcParticle.phi());
        }
      }
    } else if (mcParticle.getProcess() == 4) { // Particle decay
      if (passedITS && passedTPC) {
        h->fill(HIST(hPtItsTpcStr[histogramIndex]), mcParticle.pt());
        h->fill(HIST(hPtTrkItsTpcStr[histogramIndex]), track.pt());
        if (passedTOF) {
          h->fill(HIST(hPtItsTpcTofStr[histogramIndex]), mcParticle.pt());
        }
      }
    } else { // Material
      if (passedITS && passedTPC) {
        h->fill(HIST(hPtItsTpcMat[histogramIndex]), mcParticle.pt());
        h->fill(HIST(hPtTrkItsTpcMat[histogramIndex]), track.pt());
        if (passedTOF) {
          h->fill(HIST(hPtItsTpcTofMat[histogramIndex]), mcParticle.pt());
        }
      }
    }
  }

  template <int charge, o2::track::PID::ID id, typename particleType>
  void fillMCParticleHistograms(const particleType& mcParticle, const bool doMakeHistograms)
  {
    static_assert(charge == 0 || charge == 1);
    if (!doMakeHistograms) {
      return;
    }

    if constexpr (charge == 0) {
      if (!doPositivePDG) {
        return;
      }
    } else {
      if (!doNegativePDG) {
        return;
      }
    }

    HistogramRegistry* h = &histosPos;
    if (charge == 1) {
      h = &histosNeg;
    }

    constexpr int histogramIndex = id + charge * nSpecies;
    LOG(debug) << "fillMCParticleHistograms for charge '" << charge << "' and id '" << static_cast<int>(id) << "' " << particleName(charge, id) << " with index " << histogramIndex;
    if (!isPdgSelected<charge, id>(mcParticle)) { // Selecting PDG code
      return;
    }
    histos.fill(HIST("MC/particleSelection"), 7 + id);

    h->fill(HIST(hPGenerated[histogramIndex]), mcParticle.p());
    h->fill(HIST(hPtGenerated[histogramIndex]), mcParticle.pt());

    bool isPhysicalPrimary = mcParticle.isPhysicalPrimary();
    if (maxProdRadius < 999.f) {
      if ((mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()) > maxProdRadius * maxProdRadius) {
        isPhysicalPrimary = false;
      }
    }

    if (isPhysicalPrimary) {
      h->fill(HIST(hPtGeneratedPrm[histogramIndex]), mcParticle.pt());
      h->fill(HIST(hEtaGeneratedPrm[histogramIndex]), mcParticle.eta());
      h->fill(HIST(hPhiGeneratedPrm[histogramIndex]), mcParticle.phi());
    } else {
      if (mcParticle.getProcess() == 4) { // Particle deday
        h->fill(HIST(hPtGeneratedStr[histogramIndex]), mcParticle.pt());
      } else { // Material
        h->fill(HIST(hPtGeneratedMat[histogramIndex]), mcParticle.pt());
      }
    }

    h->fill(HIST(hEtaGenerated[histogramIndex]), mcParticle.eta());
    h->fill(HIST(hYGenerated[histogramIndex]), mcParticle.y());
    h->fill(HIST(hPhiGenerated[histogramIndex]), mcParticle.phi());
    if (doPtEta) {
      h->fill(HIST(hPtEtaGenerated[histogramIndex]), mcParticle.pt(), mcParticle.eta());
    }
  }

  template <int charge, o2::track::PID::ID id>
  void fillMCEfficiency(const bool doMakeHistograms)
  {
    if (!doMakeHistograms) {
      return;
    }

    static_assert(charge == 0 || charge == 1);
    if constexpr (charge == 0) {
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

    HistogramRegistry* registry = &histosPos;
    if (charge == 1) {
      registry = &histosNeg;
    }

    constexpr int histogramIndex = id + charge * nSpecies;

    const char* partName = particleName(charge, id);
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
      eff->SetTotalHistogram(*registry->get<TH1>(den).get(), "f");
      eff->SetPassedHistogram(*registry->get<TH1>(num).get(), "f");
    };

    doFillEfficiency("ITS_vsPt", HIST(hPtIts[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("TPC_vsPt", HIST(hPtTpc[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPt", HIST(hPtItsTpc[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("ITS-TOF_vsPt", HIST(hPtItsTof[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("Tpc-TOF_vsPt", HIST(hPtTpcTof[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TRD_vsPt", HIST(hPtItsTpcTrd[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsPt", HIST(hPtItsTpcTof[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TRD-TOF_vsPt", HIST(hPtItsTpcTrdTof[histogramIndex]), HIST(hPtGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPt_Trk", HIST(hPtTrkItsTpc[histogramIndex]), HIST(hPtGenerated[histogramIndex]));

    doFillEfficiency("ITS_vsPt_Prm", HIST(hPtItsPrm[histogramIndex]), HIST(hPtGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPt_Prm", HIST(hPtItsTpcPrm[histogramIndex]), HIST(hPtGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPt_Prm_Trk", HIST(hPtTrkItsTpcPrm[histogramIndex]), HIST(hPtGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsPt_Prm", HIST(hPtItsTpcTofPrm[histogramIndex]), HIST(hPtGeneratedPrm[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsPt_Str", HIST(hPtItsTpcStr[histogramIndex]), HIST(hPtGeneratedStr[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPt_Str_Trk", HIST(hPtTrkItsTpcStr[histogramIndex]), HIST(hPtGeneratedStr[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsPt_Str", HIST(hPtItsTpcTofStr[histogramIndex]), HIST(hPtGeneratedStr[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsPt_Mat", HIST(hPtItsTpcMat[histogramIndex]), HIST(hPtGeneratedMat[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPt_Mat_Trk", HIST(hPtTrkItsTpcMat[histogramIndex]), HIST(hPtGeneratedMat[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsPt_Mat", HIST(hPtItsTpcTofMat[histogramIndex]), HIST(hPtGeneratedMat[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsP", HIST(hPItsTpc[histogramIndex]), HIST(hPGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsP_Trk", HIST(hPTrkItsTpc[histogramIndex]), HIST(hPGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsP", HIST(hPItsTpcTof[histogramIndex]), HIST(hPGenerated[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsEta", HIST(hEtaItsTpc[histogramIndex]), HIST(hEtaGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsEta_Trk", HIST(hEtaTrkItsTpc[histogramIndex]), HIST(hEtaGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsEta", HIST(hEtaItsTpcTof[histogramIndex]), HIST(hEtaGenerated[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsEta_Prm", HIST(hEtaItsTpcPrm[histogramIndex]), HIST(hEtaGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsEta_Prm_Trk", HIST(hEtaTrkItsTpcPrm[histogramIndex]), HIST(hEtaGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsEta_Prm", HIST(hEtaItsTpcTofPrm[histogramIndex]), HIST(hEtaGeneratedPrm[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsPhi_Prm", HIST(hPhiItsTpcPrm[histogramIndex]), HIST(hPhiGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPhi_Prm_Trk", HIST(hPhiTrkItsTpcPrm[histogramIndex]), HIST(hPhiGeneratedPrm[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsPhi_Prm", HIST(hPhiItsTpcTofPrm[histogramIndex]), HIST(hPhiGeneratedPrm[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsY", HIST(hYItsTpc[histogramIndex]), HIST(hYGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsY", HIST(hYItsTpcTof[histogramIndex]), HIST(hYGenerated[histogramIndex]));

    doFillEfficiency("ITS-TPC_vsPhi", HIST(hPhiItsTpc[histogramIndex]), HIST(hPhiGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC_vsPhi_Trk", HIST(hPhiTrkItsTpc[histogramIndex]), HIST(hPhiGenerated[histogramIndex]));
    doFillEfficiency("ITS-TPC-TOF_vsPhi", HIST(hPhiItsTpcTof[histogramIndex]), HIST(hPhiGenerated[histogramIndex]));

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
      eff->SetTotalHistogram(*registry->get<TH2>(den).get(), "f");
      eff->SetPassedHistogram(*registry->get<TH2>(num).get(), "f");
    };
    fillEfficiency2D("ITS-TPC_vsPt_vsEta", HIST(hPtEtaItsTpc[histogramIndex]), HIST(hPtEtaGenerated[histogramIndex]));
    fillEfficiency2D("ITS-TPC_vsPt_vsEta_Trk", HIST(hPtEtaTrkItsTpc[histogramIndex]), HIST(hPtEtaGenerated[histogramIndex]));
    fillEfficiency2D("ITS-TPC-TOF_vsPt_vsEta", HIST(hPtEtaItsTpcTof[histogramIndex]), HIST(hPtEtaGenerated[histogramIndex]));
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
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension, o2::aod::TracksDCA>;
  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
               TrackCandidates const& tracks)
  {
    isCollisionSelected<true>(collision);
  }

  // Function to apply particle selection
  template <bool isMC = true, typename particleType, typename histoType>
  bool isInAcceptance(const particleType& particle, const histoType& countingHisto, const int offset = 0)
  {
    if (particle.pt() < ptMin || particle.pt() > ptMax) { // Check pt
      return false;
    }
    histos.fill(countingHisto, 1 + offset);
    if (particle.eta() < etaMin || particle.eta() > etaMax) { // Check eta
      return false;
    }
    histos.fill(countingHisto, 2 + offset);
    if (particle.phi() < phiMin || particle.phi() > phiMax) { // Check phi
      return false;
    }
    histos.fill(countingHisto, 3 + offset);
    if constexpr (isMC) {
      if (particle.y() < yMin || particle.y() > yMax) { // Check rapidity
        return false;
      }
      histos.fill(countingHisto, 4 + offset);
    }

    return true;
  }

  // Function to apply track selection
  bool passedITS = false;
  bool passedTPC = false;
  bool passedTRD = false;
  bool passedTOF = false;
  template <bool isMC = true, typename trackType, typename histoType>
  bool isTrackSelected(trackType& track, const histoType& countingHisto)
  {
    // Reset selections
    passedITS = false;
    passedTPC = false;
    passedTRD = false;
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
        // 2: pt cut 3: eta cut 4: phi cut 5: y cut
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
      histos.fill(countingHisto, 7);
    } else { // Data only
      if (!isInAcceptance<false>(track, countingHisto, 2)) {
        return false;
      }
    }

    if (!track.has_collision()) {
      return false;
    }
    histos.fill(countingHisto, 8);

    if (trackSelection) { // Check general cuts
      if (!track.passedTrackType()) {
        return false;
      }
      histos.fill(countingHisto, 9);
      if (!track.passedPtRange()) {
        return false;
      }
      histos.fill(countingHisto, 10);
      if (!track.passedEtaRange()) {
        return false;
      }
      histos.fill(countingHisto, 11);
      if (!track.passedDCAxy()) {
        return false;
      }
      histos.fill(countingHisto, 12);
      if (!track.passedDCAz()) {
        return false;
      }
      histos.fill(countingHisto, 13);
      if (!track.passedGoldenChi2()) {
        return false;
      }
      histos.fill(countingHisto, 14);

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
      histos.fill(countingHisto, 15);
    }

    if (passedTPC) { // Partial
      histos.fill(countingHisto, 16);
    }

    if (passedTOF) { // Partial
      histos.fill(countingHisto, 17);
    }

    switch (globalTrackSelection) {
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
      case 6:
        return customTrackCuts.IsSelected(track);
      default:
        LOG(fatal) << "Can't interpret track asked selection";
    }
    histos.fill(countingHisto, 18);

    return false;
  }

  // MC process
  Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  void processMC(o2::aod::McCollision const& mcCollision,
                 o2::soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>> const& collisions,
                 o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels> const& tracks,
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
        static_for<0, 1>([&](auto charge) {
          fillMCTrackHistograms<charge, o2::track::PID::Electron>(track, doEl);
          fillMCTrackHistograms<charge, o2::track::PID::Muon>(track, doMu);
          fillMCTrackHistograms<charge, o2::track::PID::Pion>(track, doPi);
          fillMCTrackHistograms<charge, o2::track::PID::Kaon>(track, doKa);
          fillMCTrackHistograms<charge, o2::track::PID::Proton>(track, doPr);
          fillMCTrackHistograms<charge, o2::track::PID::Deuteron>(track, doDe);
          fillMCTrackHistograms<charge, o2::track::PID::Triton>(track, doTr);
          fillMCTrackHistograms<charge, o2::track::PID::Helium3>(track, doHe);
          fillMCTrackHistograms<charge, o2::track::PID::Alpha>(track, doAl);
        });
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

      static_for<0, 1>([&](auto charge) {
        fillMCParticleHistograms<charge, o2::track::PID::Electron>(mcParticle, doEl);
        fillMCParticleHistograms<charge, o2::track::PID::Muon>(mcParticle, doMu);
        fillMCParticleHistograms<charge, o2::track::PID::Pion>(mcParticle, doPi);
        fillMCParticleHistograms<charge, o2::track::PID::Kaon>(mcParticle, doKa);
        fillMCParticleHistograms<charge, o2::track::PID::Proton>(mcParticle, doPr);
        fillMCParticleHistograms<charge, o2::track::PID::Deuteron>(mcParticle, doDe);
        fillMCParticleHistograms<charge, o2::track::PID::Triton>(mcParticle, doTr);
        fillMCParticleHistograms<charge, o2::track::PID::Helium3>(mcParticle, doHe);
        fillMCParticleHistograms<charge, o2::track::PID::Alpha>(mcParticle, doAl);
      });
    }
    histos.fill(HIST("MC/eventMultiplicity"), dNdEta * 0.5f / 2.f);

    // Fill TEfficiencies
    static_for<0, 1>([&](auto charge) {
      fillMCEfficiency<charge, o2::track::PID::Electron>(doEl);
      fillMCEfficiency<charge, o2::track::PID::Muon>(doMu);
      fillMCEfficiency<charge, o2::track::PID::Pion>(doPi);
      fillMCEfficiency<charge, o2::track::PID::Kaon>(doKa);
      fillMCEfficiency<charge, o2::track::PID::Proton>(doPr);
      fillMCEfficiency<charge, o2::track::PID::Deuteron>(doDe);
      fillMCEfficiency<charge, o2::track::PID::Triton>(doTr);
      fillMCEfficiency<charge, o2::track::PID::Helium3>(doHe);
      fillMCEfficiency<charge, o2::track::PID::Alpha>(doAl);
    });
  }
  PROCESS_SWITCH(QaEfficiency, processMC, "process MC", false);

  // MC process without the collision association
  void processMCWithoutCollisions(o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels> const& tracks,
                                  o2::aod::McParticles const& mcParticles)
  {
    // Track loop
    for (const auto& track : tracks) {
      if (!isTrackSelected(track, HIST("MC/trackSelection"))) {
        continue;
      }
      // Filling variable histograms
      histos.fill(HIST("MC/trackLength"), track.length());
      static_for<0, 1>([&](auto charge) {
        fillMCTrackHistograms<charge, o2::track::PID::Electron>(track, doEl);
        fillMCTrackHistograms<charge, o2::track::PID::Muon>(track, doMu);
        fillMCTrackHistograms<charge, o2::track::PID::Pion>(track, doPi);
        fillMCTrackHistograms<charge, o2::track::PID::Kaon>(track, doKa);
        fillMCTrackHistograms<charge, o2::track::PID::Proton>(track, doPr);
        fillMCTrackHistograms<charge, o2::track::PID::Deuteron>(track, doDe);
        fillMCTrackHistograms<charge, o2::track::PID::Triton>(track, doTr);
        fillMCTrackHistograms<charge, o2::track::PID::Helium3>(track, doHe);
        fillMCTrackHistograms<charge, o2::track::PID::Alpha>(track, doAl);
      });
    }

    for (const auto& mcParticle : mcParticles) {
      if (!isInAcceptance(mcParticle, HIST("MC/particleSelection"))) {
        continue;
      }

      static_for<0, 1>([&](auto charge) {
        fillMCParticleHistograms<charge, o2::track::PID::Electron>(mcParticle, doEl);
        fillMCParticleHistograms<charge, o2::track::PID::Muon>(mcParticle, doMu);
        fillMCParticleHistograms<charge, o2::track::PID::Pion>(mcParticle, doPi);
        fillMCParticleHistograms<charge, o2::track::PID::Kaon>(mcParticle, doKa);
        fillMCParticleHistograms<charge, o2::track::PID::Proton>(mcParticle, doPr);
        fillMCParticleHistograms<charge, o2::track::PID::Deuteron>(mcParticle, doDe);
        fillMCParticleHistograms<charge, o2::track::PID::Triton>(mcParticle, doTr);
        fillMCParticleHistograms<charge, o2::track::PID::Helium3>(mcParticle, doHe);
        fillMCParticleHistograms<charge, o2::track::PID::Alpha>(mcParticle, doAl);
      });
    }

    // Fill TEfficiencies
    static_for<0, 1>([&](auto charge) {
      fillMCEfficiency<charge, o2::track::PID::Electron>(doEl);
      fillMCEfficiency<charge, o2::track::PID::Muon>(doMu);
      fillMCEfficiency<charge, o2::track::PID::Pion>(doPi);
      fillMCEfficiency<charge, o2::track::PID::Kaon>(doKa);
      fillMCEfficiency<charge, o2::track::PID::Proton>(doPr);
      fillMCEfficiency<charge, o2::track::PID::Deuteron>(doDe);
      fillMCEfficiency<charge, o2::track::PID::Triton>(doTr);
      fillMCEfficiency<charge, o2::track::PID::Helium3>(doHe);
      fillMCEfficiency<charge, o2::track::PID::Alpha>(doAl);
    });
  }
  PROCESS_SWITCH(QaEfficiency, processMCWithoutCollisions, "process MC without the collision association", false);

  void processData(o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator const& collision,
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
          static_cast<TEfficiency*>(listEfficiencyData->At(0))->Fill(passedTPC, track.pt());
        }
        if (passedITS && passedTPC) {
          static_cast<TEfficiency*>(listEfficiencyData->At(1))->Fill(passedTOF, track.pt());
          static_cast<TEfficiency*>(listEfficiencyData->At(2))->Fill(passedTOF, track.p());
          static_cast<TEfficiency*>(listEfficiencyData->At(3))->Fill(passedTOF, track.eta());
          static_cast<TEfficiency*>(listEfficiencyData->At(4))->Fill(passedTOF, track.phi());
          static_cast<TEfficiency*>(listEfficiencyData->At(5))->Fill(passedTOF, track.pt(), track.eta());
          static_cast<TEfficiency*>(listEfficiencyData->At(6))->Fill(passedTOF, track.pt(), track.phi());
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
