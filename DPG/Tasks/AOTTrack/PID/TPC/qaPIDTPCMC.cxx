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
/// \file   qaPIDTPCMC.cxx
/// \author Nicolò Jacazio
/// \brief  Task to produce QA output of the PID with TPC running on the MC.
///

// O2 includes
#include "Common/DataModel/PIDResponseTPC.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TPC QA plots
struct pidTpcQaMc {
  SliceCache cache;

  static constexpr int Np = 9;
  static constexpr int NpNp = Np * Np;
  static constexpr std::string_view hparticlept[Np] = {"particlept/El", "particlept/Mu", "particlept/Pi",
                                                       "particlept/Ka", "particlept/Pr", "particlept/De",
                                                       "particlept/Tr", "particlept/He", "particlept/Al"};
  static constexpr std::string_view hparticlep[Np] = {"particlep/El", "particlep/Mu", "particlep/Pi",
                                                      "particlep/Ka", "particlep/Pr", "particlep/De",
                                                      "particlep/Tr", "particlep/He", "particlep/Al"};
  static constexpr std::string_view hparticleeta[Np] = {"particleeta/El", "particleeta/Mu", "particleeta/Pi",
                                                        "particleeta/Ka", "particleeta/Pr", "particleeta/De",
                                                        "particleeta/Tr", "particleeta/He", "particleeta/Al"};
  static constexpr std::string_view htrackpt[Np] = {"trackpt/El", "trackpt/Mu", "trackpt/Pi",
                                                    "trackpt/Ka", "trackpt/Pr", "trackpt/De",
                                                    "trackpt/Tr", "trackpt/He", "trackpt/Al"};
  static constexpr std::string_view htrackp[Np] = {"trackp/El", "trackp/Mu", "trackp/Pi",
                                                   "trackp/Ka", "trackp/Pr", "trackp/De",
                                                   "trackp/Tr", "trackp/He", "trackp/Al"};
  static constexpr std::string_view htracketa[Np] = {"tracketa/El", "tracketa/Mu", "tracketa/Pi",
                                                     "tracketa/Ka", "tracketa/Pr", "tracketa/De",
                                                     "tracketa/Tr", "tracketa/He", "tracketa/Al"};
  static constexpr std::string_view htracklength[Np] = {"tracklength/El", "tracklength/Mu", "tracklength/Pi",
                                                        "tracklength/Ka", "tracklength/Pr", "tracklength/De",
                                                        "tracklength/Tr", "tracklength/He", "tracklength/Al"};
  // Signal
  static constexpr std::string_view hsignalMC[Np] = {"signalMC/El", "signalMC/Mu", "signalMC/Pi",
                                                     "signalMC/Ka", "signalMC/Pr", "signalMC/De",
                                                     "signalMC/Tr", "signalMC/He", "signalMC/Al"};
  static constexpr std::string_view hsignalMCprm[Np] = {"signalMCprm/El", "signalMCprm/Mu", "signalMCprm/Pi",
                                                        "signalMCprm/Ka", "signalMCprm/Pr", "signalMCprm/De",
                                                        "signalMCprm/Tr", "signalMCprm/He", "signalMCprm/Al"};
  static constexpr std::string_view hsignalMCstr[Np] = {"signalMCstr/El", "signalMCstr/Mu", "signalMCstr/Pi",
                                                        "signalMCstr/Ka", "signalMCstr/Pr", "signalMCstr/De",
                                                        "signalMCstr/Tr", "signalMCstr/He", "signalMCstr/Al"};
  static constexpr std::string_view hsignalMCmat[Np] = {"signalMCmat/El", "signalMCmat/Mu", "signalMCmat/Pi",
                                                        "signalMCmat/Ka", "signalMCmat/Pr", "signalMCmat/De",
                                                        "signalMCmat/Tr", "signalMCmat/He", "signalMCmat/Al"};
  // Nsigma
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigmaprm[Np] = {"nsigmaprm/El", "nsigmaprm/Mu", "nsigmaprm/Pi",
                                                      "nsigmaprm/Ka", "nsigmaprm/Pr", "nsigmaprm/De",
                                                      "nsigmaprm/Tr", "nsigmaprm/He", "nsigmaprm/Al"};
  static constexpr std::string_view hnsigmastr[Np] = {"nsigmastr/El", "nsigmastr/Mu", "nsigmastr/Pi",
                                                      "nsigmastr/Ka", "nsigmastr/Pr", "nsigmastr/De",
                                                      "nsigmastr/Tr", "nsigmastr/He", "nsigmastr/Al"};
  static constexpr std::string_view hnsigmamat[Np] = {"nsigmamat/El", "nsigmamat/Mu", "nsigmamat/Pi",
                                                      "nsigmamat/Ka", "nsigmamat/Pr", "nsigmamat/De",
                                                      "nsigmamat/Tr", "nsigmamat/He", "nsigmamat/Al"};
  static constexpr std::string_view hnsigmaMC[NpNp] = {"nsigmaMC/El/El", "nsigmaMC/El/Mu", "nsigmaMC/El/Pi",
                                                       "nsigmaMC/El/Ka", "nsigmaMC/El/Pr", "nsigmaMC/El/De",
                                                       "nsigmaMC/El/Tr", "nsigmaMC/El/He", "nsigmaMC/El/Al",
                                                       "nsigmaMC/Mu/El", "nsigmaMC/Mu/Mu", "nsigmaMC/Mu/Pi",
                                                       "nsigmaMC/Mu/Ka", "nsigmaMC/Mu/Pr", "nsigmaMC/Mu/De",
                                                       "nsigmaMC/Mu/Tr", "nsigmaMC/Mu/He", "nsigmaMC/Mu/Al",
                                                       "nsigmaMC/Pi/El", "nsigmaMC/Pi/Mu", "nsigmaMC/Pi/Pi",
                                                       "nsigmaMC/Pi/Ka", "nsigmaMC/Pi/Pr", "nsigmaMC/Pi/De",
                                                       "nsigmaMC/Pi/Tr", "nsigmaMC/Pi/He", "nsigmaMC/Pi/Al",
                                                       "nsigmaMC/Ka/El", "nsigmaMC/Ka/Mu", "nsigmaMC/Ka/Pi",
                                                       "nsigmaMC/Ka/Ka", "nsigmaMC/Ka/Pr", "nsigmaMC/Ka/De",
                                                       "nsigmaMC/Ka/Tr", "nsigmaMC/Ka/He", "nsigmaMC/Ka/Al",
                                                       "nsigmaMC/Pr/El", "nsigmaMC/Pr/Mu", "nsigmaMC/Pr/Pi",
                                                       "nsigmaMC/Pr/Ka", "nsigmaMC/Pr/Pr", "nsigmaMC/Pr/De",
                                                       "nsigmaMC/Pr/Tr", "nsigmaMC/Pr/He", "nsigmaMC/Pr/Al",
                                                       "nsigmaMC/De/El", "nsigmaMC/De/Mu", "nsigmaMC/De/Pi",
                                                       "nsigmaMC/De/Ka", "nsigmaMC/De/Pr", "nsigmaMC/De/De",
                                                       "nsigmaMC/De/Tr", "nsigmaMC/De/He", "nsigmaMC/De/Al",
                                                       "nsigmaMC/Tr/El", "nsigmaMC/Tr/Mu", "nsigmaMC/Tr/Pi",
                                                       "nsigmaMC/Tr/Ka", "nsigmaMC/Tr/Pr", "nsigmaMC/Tr/De",
                                                       "nsigmaMC/Tr/Tr", "nsigmaMC/Tr/He", "nsigmaMC/Tr/Al",
                                                       "nsigmaMC/He/El", "nsigmaMC/He/Mu", "nsigmaMC/He/Pi",
                                                       "nsigmaMC/He/Ka", "nsigmaMC/He/Pr", "nsigmaMC/He/De",
                                                       "nsigmaMC/He/Tr", "nsigmaMC/He/He", "nsigmaMC/He/Al",
                                                       "nsigmaMC/Al/El", "nsigmaMC/Al/Mu", "nsigmaMC/Al/Pi",
                                                       "nsigmaMC/Al/Ka", "nsigmaMC/Al/Pr", "nsigmaMC/Al/De",
                                                       "nsigmaMC/Al/Tr", "nsigmaMC/Al/He", "nsigmaMC/Al/Al"};
  static constexpr std::string_view hnsigmaMCstr[NpNp] = {"nsigmaMCstr/El/El", "nsigmaMCstr/El/Mu", "nsigmaMCstr/El/Pi",
                                                          "nsigmaMCstr/El/Ka", "nsigmaMCstr/El/Pr", "nsigmaMCstr/El/De",
                                                          "nsigmaMCstr/El/Tr", "nsigmaMCstr/El/He", "nsigmaMCstr/El/Al",
                                                          "nsigmaMCstr/Mu/El", "nsigmaMCstr/Mu/Mu", "nsigmaMCstr/Mu/Pi",
                                                          "nsigmaMCstr/Mu/Ka", "nsigmaMCstr/Mu/Pr", "nsigmaMCstr/Mu/De",
                                                          "nsigmaMCstr/Mu/Tr", "nsigmaMCstr/Mu/He", "nsigmaMCstr/Mu/Al",
                                                          "nsigmaMCstr/Pi/El", "nsigmaMCstr/Pi/Mu", "nsigmaMCstr/Pi/Pi",
                                                          "nsigmaMCstr/Pi/Ka", "nsigmaMCstr/Pi/Pr", "nsigmaMCstr/Pi/De",
                                                          "nsigmaMCstr/Pi/Tr", "nsigmaMCstr/Pi/He", "nsigmaMCstr/Pi/Al",
                                                          "nsigmaMCstr/Ka/El", "nsigmaMCstr/Ka/Mu", "nsigmaMCstr/Ka/Pi",
                                                          "nsigmaMCstr/Ka/Ka", "nsigmaMCstr/Ka/Pr", "nsigmaMCstr/Ka/De",
                                                          "nsigmaMCstr/Ka/Tr", "nsigmaMCstr/Ka/He", "nsigmaMCstr/Ka/Al",
                                                          "nsigmaMCstr/Pr/El", "nsigmaMCstr/Pr/Mu", "nsigmaMCstr/Pr/Pi",
                                                          "nsigmaMCstr/Pr/Ka", "nsigmaMCstr/Pr/Pr", "nsigmaMCstr/Pr/De",
                                                          "nsigmaMCstr/Pr/Tr", "nsigmaMCstr/Pr/He", "nsigmaMCstr/Pr/Al",
                                                          "nsigmaMCstr/De/El", "nsigmaMCstr/De/Mu", "nsigmaMCstr/De/Pi",
                                                          "nsigmaMCstr/De/Ka", "nsigmaMCstr/De/Pr", "nsigmaMCstr/De/De",
                                                          "nsigmaMCstr/De/Tr", "nsigmaMCstr/De/He", "nsigmaMCstr/De/Al",
                                                          "nsigmaMCstr/Tr/El", "nsigmaMCstr/Tr/Mu", "nsigmaMCstr/Tr/Pi",
                                                          "nsigmaMCstr/Tr/Ka", "nsigmaMCstr/Tr/Pr", "nsigmaMCstr/Tr/De",
                                                          "nsigmaMCstr/Tr/Tr", "nsigmaMCstr/Tr/He", "nsigmaMCstr/Tr/Al",
                                                          "nsigmaMCstr/He/El", "nsigmaMCstr/He/Mu", "nsigmaMCstr/He/Pi",
                                                          "nsigmaMCstr/He/Ka", "nsigmaMCstr/He/Pr", "nsigmaMCstr/He/De",
                                                          "nsigmaMCstr/He/Tr", "nsigmaMCstr/He/He", "nsigmaMCstr/He/Al",
                                                          "nsigmaMCstr/Al/El", "nsigmaMCstr/Al/Mu", "nsigmaMCstr/Al/Pi",
                                                          "nsigmaMCstr/Al/Ka", "nsigmaMCstr/Al/Pr", "nsigmaMCstr/Al/De",
                                                          "nsigmaMCstr/Al/Tr", "nsigmaMCstr/Al/He", "nsigmaMCstr/Al/Al"};
  static constexpr std::string_view hnsigmaMCmat[NpNp] = {"nsigmaMCmat/El/El", "nsigmaMCmat/El/Mu", "nsigmaMCmat/El/Pi",
                                                          "nsigmaMCmat/El/Ka", "nsigmaMCmat/El/Pr", "nsigmaMCmat/El/De",
                                                          "nsigmaMCmat/El/Tr", "nsigmaMCmat/El/He", "nsigmaMCmat/El/Al",
                                                          "nsigmaMCmat/Mu/El", "nsigmaMCmat/Mu/Mu", "nsigmaMCmat/Mu/Pi",
                                                          "nsigmaMCmat/Mu/Ka", "nsigmaMCmat/Mu/Pr", "nsigmaMCmat/Mu/De",
                                                          "nsigmaMCmat/Mu/Tr", "nsigmaMCmat/Mu/He", "nsigmaMCmat/Mu/Al",
                                                          "nsigmaMCmat/Pi/El", "nsigmaMCmat/Pi/Mu", "nsigmaMCmat/Pi/Pi",
                                                          "nsigmaMCmat/Pi/Ka", "nsigmaMCmat/Pi/Pr", "nsigmaMCmat/Pi/De",
                                                          "nsigmaMCmat/Pi/Tr", "nsigmaMCmat/Pi/He", "nsigmaMCmat/Pi/Al",
                                                          "nsigmaMCmat/Ka/El", "nsigmaMCmat/Ka/Mu", "nsigmaMCmat/Ka/Pi",
                                                          "nsigmaMCmat/Ka/Ka", "nsigmaMCmat/Ka/Pr", "nsigmaMCmat/Ka/De",
                                                          "nsigmaMCmat/Ka/Tr", "nsigmaMCmat/Ka/He", "nsigmaMCmat/Ka/Al",
                                                          "nsigmaMCmat/Pr/El", "nsigmaMCmat/Pr/Mu", "nsigmaMCmat/Pr/Pi",
                                                          "nsigmaMCmat/Pr/Ka", "nsigmaMCmat/Pr/Pr", "nsigmaMCmat/Pr/De",
                                                          "nsigmaMCmat/Pr/Tr", "nsigmaMCmat/Pr/He", "nsigmaMCmat/Pr/Al",
                                                          "nsigmaMCmat/De/El", "nsigmaMCmat/De/Mu", "nsigmaMCmat/De/Pi",
                                                          "nsigmaMCmat/De/Ka", "nsigmaMCmat/De/Pr", "nsigmaMCmat/De/De",
                                                          "nsigmaMCmat/De/Tr", "nsigmaMCmat/De/He", "nsigmaMCmat/De/Al",
                                                          "nsigmaMCmat/Tr/El", "nsigmaMCmat/Tr/Mu", "nsigmaMCmat/Tr/Pi",
                                                          "nsigmaMCmat/Tr/Ka", "nsigmaMCmat/Tr/Pr", "nsigmaMCmat/Tr/De",
                                                          "nsigmaMCmat/Tr/Tr", "nsigmaMCmat/Tr/He", "nsigmaMCmat/Tr/Al",
                                                          "nsigmaMCmat/He/El", "nsigmaMCmat/He/Mu", "nsigmaMCmat/He/Pi",
                                                          "nsigmaMCmat/He/Ka", "nsigmaMCmat/He/Pr", "nsigmaMCmat/He/De",
                                                          "nsigmaMCmat/He/Tr", "nsigmaMCmat/He/He", "nsigmaMCmat/He/Al",
                                                          "nsigmaMCmat/Al/El", "nsigmaMCmat/Al/Mu", "nsigmaMCmat/Al/Pi",
                                                          "nsigmaMCmat/Al/Ka", "nsigmaMCmat/Al/Pr", "nsigmaMCmat/Al/De",
                                                          "nsigmaMCmat/Al/Tr", "nsigmaMCmat/Al/He", "nsigmaMCmat/Al/Al"};
  static constexpr std::string_view hnsigmaMCprm[NpNp] = {"nsigmaMCprm/El/El", "nsigmaMCprm/El/Mu", "nsigmaMCprm/El/Pi",
                                                          "nsigmaMCprm/El/Ka", "nsigmaMCprm/El/Pr", "nsigmaMCprm/El/De",
                                                          "nsigmaMCprm/El/Tr", "nsigmaMCprm/El/He", "nsigmaMCprm/El/Al",
                                                          "nsigmaMCprm/Mu/El", "nsigmaMCprm/Mu/Mu", "nsigmaMCprm/Mu/Pi",
                                                          "nsigmaMCprm/Mu/Ka", "nsigmaMCprm/Mu/Pr", "nsigmaMCprm/Mu/De",
                                                          "nsigmaMCprm/Mu/Tr", "nsigmaMCprm/Mu/He", "nsigmaMCprm/Mu/Al",
                                                          "nsigmaMCprm/Pi/El", "nsigmaMCprm/Pi/Mu", "nsigmaMCprm/Pi/Pi",
                                                          "nsigmaMCprm/Pi/Ka", "nsigmaMCprm/Pi/Pr", "nsigmaMCprm/Pi/De",
                                                          "nsigmaMCprm/Pi/Tr", "nsigmaMCprm/Pi/He", "nsigmaMCprm/Pi/Al",
                                                          "nsigmaMCprm/Ka/El", "nsigmaMCprm/Ka/Mu", "nsigmaMCprm/Ka/Pi",
                                                          "nsigmaMCprm/Ka/Ka", "nsigmaMCprm/Ka/Pr", "nsigmaMCprm/Ka/De",
                                                          "nsigmaMCprm/Ka/Tr", "nsigmaMCprm/Ka/He", "nsigmaMCprm/Ka/Al",
                                                          "nsigmaMCprm/Pr/El", "nsigmaMCprm/Pr/Mu", "nsigmaMCprm/Pr/Pi",
                                                          "nsigmaMCprm/Pr/Ka", "nsigmaMCprm/Pr/Pr", "nsigmaMCprm/Pr/De",
                                                          "nsigmaMCprm/Pr/Tr", "nsigmaMCprm/Pr/He", "nsigmaMCprm/Pr/Al",
                                                          "nsigmaMCprm/De/El", "nsigmaMCprm/De/Mu", "nsigmaMCprm/De/Pi",
                                                          "nsigmaMCprm/De/Ka", "nsigmaMCprm/De/Pr", "nsigmaMCprm/De/De",
                                                          "nsigmaMCprm/De/Tr", "nsigmaMCprm/De/He", "nsigmaMCprm/De/Al",
                                                          "nsigmaMCprm/Tr/El", "nsigmaMCprm/Tr/Mu", "nsigmaMCprm/Tr/Pi",
                                                          "nsigmaMCprm/Tr/Ka", "nsigmaMCprm/Tr/Pr", "nsigmaMCprm/Tr/De",
                                                          "nsigmaMCprm/Tr/Tr", "nsigmaMCprm/Tr/He", "nsigmaMCprm/Tr/Al",
                                                          "nsigmaMCprm/He/El", "nsigmaMCprm/He/Mu", "nsigmaMCprm/He/Pi",
                                                          "nsigmaMCprm/He/Ka", "nsigmaMCprm/He/Pr", "nsigmaMCprm/He/De",
                                                          "nsigmaMCprm/He/Tr", "nsigmaMCprm/He/He", "nsigmaMCprm/He/Al",
                                                          "nsigmaMCprm/Al/El", "nsigmaMCprm/Al/Mu", "nsigmaMCprm/Al/Pi",
                                                          "nsigmaMCprm/Al/Ka", "nsigmaMCprm/Al/Pr", "nsigmaMCprm/Al/De",
                                                          "nsigmaMCprm/Al/Tr", "nsigmaMCprm/Al/He", "nsigmaMCprm/Al/Al"};

  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr int PDGs[Np] = {11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> checkPrimaries{"checkPrimaries", 1,
                                   "Whether to check physical primary and secondaries particles for the resolution."};
  Configurable<int> pdgSign{"pdgSign", 0, "Sign of the PDG, -1 0 or 1"};
  Configurable<int> doEl{"doEl", 0, "Process electrons"};
  Configurable<int> doMu{"doMu", 0, "Process muons"};
  Configurable<int> doPi{"doPi", 0, "Process pions"};
  Configurable<int> doKa{"doKa", 0, "Process kaons"};
  Configurable<int> doPr{"doPr", 0, "Process protons"};
  Configurable<int> doDe{"doDe", 0, "Process deuterons"};
  Configurable<int> doTr{"doTr", 0, "Process tritons"};
  Configurable<int> doHe{"doHe", 0, "Process helium3"};
  Configurable<int> doAl{"doAl", 0, "Process alpha"};
  ConfigurableAxis binsPt{"binsPt", {2000, 0.f, 20.f}, "Binning of the pT axis"};
  ConfigurableAxis binsNsigma{"binsNsigma", {2000, -30.f, 30.f}, "Binning of the NSigma axis"};
  ConfigurableAxis binsSignal{"binsSignal", {6000, 0, 2000}, "Binning of the TPC signal axis"};
  ConfigurableAxis binsLength{"binsLength", {1000, 0, 3000}, "Binning of the Length axis"};
  ConfigurableAxis binsEta{"binsEta", {100, -4, 4}, "Binning of the Eta axis"};
  Configurable<float> minEta{"minEta", -0.8, "Minimum eta in range"};
  Configurable<float> maxEta{"maxEta", 0.8, "Maximum eta in range"};
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the vertex"};
  Configurable<int> logAxis{"logAxis", 0, "Flag to use a logarithmic pT axis, in this case the pT limits are the expontents"};

  template <uint8_t mcID, uint8_t massID>
  void addParticleHistos(const AxisSpec& ptAxis, const AxisSpec& pAxis, const AxisSpec& signalAxis)
  {
    switch (mcID) {
      case 0:
        if (!doEl) {
          return;
        }
        break;
      case 1:
        if (!doMu) {
          return;
        }
        break;
      case 2:
        if (!doPi) {
          return;
        }
        break;
      case 3:
        if (!doKa) {
          return;
        }
        break;
      case 4:
        if (!doPr) {
          return;
        }
        break;
      case 5:
        if (!doDe) {
          return;
        }
        break;
      case 6:
        if (!doTr) {
          return;
        }
        break;
      case 7:
        if (!doHe) {
          return;
        }
        break;
      case 8:
        if (!doAl) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Can't interpret index";
    }

    const AxisSpec lengthAxis{binsLength, "Track length (cm)"};
    const AxisSpec etaAxis{binsEta, "#it{#eta}"};
    const AxisSpec nSigmaAxis{binsNsigma, Form("N_{#sigma}^{TPC}(%s)", pT[massID])};

    // Particle info
    if constexpr (mcID == massID) {
      histos.add(hparticlept[mcID].data(), "", kTH1F, {pAxis});
      histos.add(hparticlep[mcID].data(), "", kTH1F, {pAxis});
      histos.add(hparticleeta[mcID].data(), "", kTH1F, {pAxis});
    }
    // Track info
    if constexpr (mcID == massID) {
      histos.add(htrackpt[mcID].data(), "", kTH1F, {pAxis});
      histos.add(htrackp[mcID].data(), "", kTH1F, {pAxis});
      histos.add(htracketa[mcID].data(), "", kTH1F, {pAxis});
      histos.add(htracklength[mcID].data(), "", kTH1F, {pAxis});
    }
    // NSigma
    if constexpr (mcID == massID) {
      histos.add(hnsigma[mcID].data(), pT[mcID], HistType::kTH2F, {ptAxis, nSigmaAxis});
    }
    histos.add(hnsigmaMC[mcID * Np + massID].data(), Form("True %s", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    if (!checkPrimaries) {
      return;
    }
    if constexpr (mcID == massID) {
      histos.add(hnsigmaprm[mcID].data(), Form("Primary %s", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigmastr[mcID].data(), Form("Secondary %s from decay", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigmamat[mcID].data(), Form("Secondary %s from material", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    }
    histos.add(hnsigmaMCprm[mcID * Np + massID].data(), Form("True Primary %s", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmaMCstr[mcID * Np + massID].data(), Form("True Secondary %s from decay", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmaMCmat[mcID * Np + massID].data(), Form("True Secondary %s from material", pT[mcID]), HistType::kTH2F, {ptAxis, nSigmaAxis});

    if constexpr (mcID == massID) {
      histos.add(hsignalMC[mcID].data(), Form("%s", pT[mcID]), HistType::kTH2F, {pAxis, signalAxis});
      histos.add(hsignalMCprm[mcID].data(), Form("Primary %s", pT[mcID]), HistType::kTH2F, {pAxis, signalAxis});
      histos.add(hsignalMCstr[mcID].data(), Form("Secondary %s from decay", pT[mcID]), HistType::kTH2F, {pAxis, signalAxis});
      histos.add(hsignalMCmat[mcID].data(), Form("Secondary %s from material", pT[mcID]), HistType::kTH2F, {pAxis, signalAxis});
    }
  }

  void init(o2::framework::InitContext&)
  {
    AxisSpec pAxis{binsPt, "#it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogarithmic();
      ptAxis.makeLogarithmic();
    }
    const AxisSpec signalAxis{binsSignal, "TPC d#it{E}/d#it{x} A.U."};

    histos.add("event/vertexz", ";Vtx_{z} (cm);Entries", kTH1F, {{100, -20, 20}});

    static_for<0, 8>([&](auto i) {
      static_for<0, 8>([&](auto j) {
        addParticleHistos<i, j>(ptAxis, pAxis, signalAxis);
      });
    });

    histos.add("event/tpcsignal", "All", HistType::kTH2F, {pAxis, signalAxis});
    if (checkPrimaries) {
      histos.add("event/tpcsignalPrm", "Primaries", HistType::kTH2F, {pAxis, signalAxis});
      histos.add("event/tpcsignalStr", "Secondaries from weak decays", HistType::kTH2F, {pAxis, signalAxis});
      histos.add("event/tpcsignalMat", "Secondaries from material", HistType::kTH2F, {pAxis, signalAxis});
    }
  }

  template <uint8_t mcID, typename T>
  void fillParticleInfoForPdg(const T& particle)
  {
    switch (pdgSign.value) {
      case 0:
        if (abs(particle.pdgCode()) != PDGs[mcID]) {
          return;
        }
        break;
      case 1:
        if (particle.pdgCode() != PDGs[mcID]) {
          return;
        }
        break;
      case 2:
        if (particle.pdgCode() != -PDGs[mcID]) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Can't interpret pdgSign";
    }
    switch (mcID) {
      case 0:
        if (!doEl) {
          return;
        }
        break;
      case 1:
        if (!doMu) {
          return;
        }
        break;
      case 2:
        if (!doPi) {
          return;
        }
        break;
      case 3:
        if (!doKa) {
          return;
        }
        break;
      case 4:
        if (!doPr) {
          return;
        }
        break;
      case 5:
        if (!doDe) {
          return;
        }
        break;
      case 6:
        if (!doTr) {
          return;
        }
        break;
      case 7:
        if (!doHe) {
          return;
        }
        break;
      case 8:
        if (!doAl) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Can't interpret index";
    }

    histos.fill(HIST(hparticlep[mcID]), particle.p());
    histos.fill(HIST(hparticlept[mcID]), particle.pt());
    histos.fill(HIST(hparticleeta[mcID]), particle.eta());
  }

  template <uint8_t mcID, typename T, typename TT>
  void fillTrackInfoForPdg(const T& track, const TT& particle)
  {

    switch (mcID) {
      case 0:
        if (!doEl) {
          return;
        }
        break;
      case 1:
        if (!doMu) {
          return;
        }
        break;
      case 2:
        if (!doPi) {
          return;
        }
        break;
      case 3:
        if (!doKa) {
          return;
        }
        break;
      case 4:
        if (!doPr) {
          return;
        }
        break;
      case 5:
        if (!doDe) {
          return;
        }
        break;
      case 6:
        if (!doTr) {
          return;
        }
        break;
      case 7:
        if (!doHe) {
          return;
        }
        break;
      case 8:
        if (!doAl) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Can't interpret index";
    }

    const float nsigma = o2::aod::pidutils::tpcNSigma<mcID>(track);

    // Fill for all
    histos.fill(HIST(hnsigma[mcID]), track.pt(), nsigma);

    if (checkPrimaries) {
      if (!particle.isPhysicalPrimary()) {
        if (particle.getProcess() == 4) {
          histos.fill(HIST(hnsigmastr[mcID]), track.pt(), nsigma);
        } else {
          histos.fill(HIST(hnsigmamat[mcID]), track.pt(), nsigma);
        }
      } else {
        histos.fill(HIST(hnsigmaprm[mcID]), track.pt(), nsigma);
      }
    }

    switch (pdgSign.value) {
      case 0:
        if (abs(particle.pdgCode()) != PDGs[mcID]) {
          return;
        }
        break;
      case 1:
        if (particle.pdgCode() != PDGs[mcID]) {
          return;
        }
        break;
      case 2:
        if (particle.pdgCode() != -PDGs[mcID]) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Can't interpret pdgSign";
    }

    // Track info
    histos.fill(HIST(htrackp[mcID]), track.p());
    histos.fill(HIST(htrackpt[mcID]), track.pt());
    histos.fill(HIST(htracketa[mcID]), track.eta());
    histos.fill(HIST(htracklength[mcID]), track.length());

    // PID info
    histos.fill(HIST(hsignalMC[mcID]), track.tpcInnerParam(), track.tpcSignal());
    histos.fill(HIST(hnsigmaMC[mcID]), track.pt(), nsigma);
    if (!particle.isPhysicalPrimary()) {
      if (particle.getProcess() == 4) {
        histos.fill(HIST(hsignalMCstr[mcID]), track.tpcInnerParam(), track.tpcSignal());
        histos.fill(HIST(hnsigmaMCstr[mcID]), track.pt(), nsigma);
      } else {
        histos.fill(HIST(hsignalMCmat[mcID]), track.tpcInnerParam(), track.tpcSignal());
        histos.fill(HIST(hnsigmaMCmat[mcID]), track.pt(), nsigma);
      }
    } else {
      histos.fill(HIST(hsignalMCprm[mcID]), track.tpcInnerParam(), track.tpcSignal());
      histos.fill(HIST(hnsigmaMCprm[mcID]), track.pt(), nsigma);
    }
  }

  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
               soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                         aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                         aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                         aod::McTrackLabels>& tracks,
               aod::McParticles& mcParticles,
               aod::McCollisions&)
  {
    for (const auto& collision : collisions) {
      if (collision.numContrib() < nMinNumberOfContributors) {
        return;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto tracksInCollision = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache);

      for (const auto& p : particlesInCollision) {
        static_for<0, 8>([&](auto i) {
          fillParticleInfoForPdg<i>(p);
        });
      }

      for (const auto& t : tracksInCollision) {
        if (t.eta() < minEta || t.eta() > maxEta) {
          continue;
        }

        // Fill for all
        histos.fill(HIST("event/tpcsignal"), t.tpcInnerParam(), t.tpcSignal());
        if (!t.has_mcParticle()) {
          continue;
        }

        const auto& particle = t.mcParticle();

        if (!particle.isPhysicalPrimary()) {
          if (particle.getProcess() == 4) {
            histos.fill(HIST("event/tpcsignalStr"), t.tpcInnerParam(), t.tpcSignal());
          } else {
            histos.fill(HIST("event/tpcsignalMat"), t.tpcInnerParam(), t.tpcSignal());
          }
        } else {
          histos.fill(HIST("event/tpcsignalPrm"), t.tpcInnerParam(), t.tpcSignal());
        }

        // Fill with PDG codes
        static_for<0, 8>([&](auto i) {
          fillTrackInfoForPdg<i>(t, particle);
        });
      }
      histos.fill(HIST("event/vertexz"), collision.posZ());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<pidTpcQaMc>(cfgc)}; }
