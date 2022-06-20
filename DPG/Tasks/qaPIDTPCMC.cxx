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
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"qa-el", VariantType::Int, 0, {"Produce PID information for the electron mass hypothesis"}},
    {"qa-mu", VariantType::Int, 0, {"Produce PID information for the muon mass hypothesis"}},
    {"qa-pikapr", VariantType::Int, 1, {"Produce PID information for the Pion, Kaon, Proton mass hypothesis"}},
    {"qa-nuclei", VariantType::Int, 0, {"Produce PID information for the Deuteron, Triton, Alpha mass hypothesis"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

template <o2::track::PID::ID pid_type>
struct pidTPCTaskQA {

  static constexpr int Np = 9;
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigmaprm[Np] = {"nsigmaprm/El", "nsigmaprm/Mu", "nsigmaprm/Pi",
                                                      "nsigmaprm/Ka", "nsigmaprm/Pr", "nsigmaprm/De",
                                                      "nsigmaprm/Tr", "nsigmaprm/He", "nsigmaprm/Al"};
  static constexpr std::string_view hnsigmasec[Np] = {"nsigmasec/El", "nsigmasec/Mu", "nsigmasec/Pi",
                                                      "nsigmasec/Ka", "nsigmasec/Pr", "nsigmasec/De",
                                                      "nsigmasec/Tr", "nsigmasec/He", "nsigmasec/Al"};
  static constexpr std::string_view hnsigmaMC[Np] = {"nsigmaMC/El", "nsigmaMC/Mu", "nsigmaMC/Pi",
                                                     "nsigmaMC/Ka", "nsigmaMC/Pr", "nsigmaMC/De",
                                                     "nsigmaMC/Tr", "nsigmaMC/He", "nsigmaMC/Al"};
  static constexpr std::string_view hnsigmaMCsec[Np] = {"nsigmaMCsec/El", "nsigmaMCsec/Mu", "nsigmaMCsec/Pi",
                                                        "nsigmaMCsec/Ka", "nsigmaMCsec/Pr", "nsigmaMCsec/De",
                                                        "nsigmaMCsec/Tr", "nsigmaMCsec/He", "nsigmaMCsec/Al"};
  static constexpr std::string_view hnsigmaMCprm[Np] = {"nsigmaMCprm/El", "nsigmaMCprm/Mu", "nsigmaMCprm/Pi",
                                                        "nsigmaMCprm/Ka", "nsigmaMCprm/Pr", "nsigmaMCprm/De",
                                                        "nsigmaMCprm/Tr", "nsigmaMCprm/He", "nsigmaMCprm/Al"};
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr int PDGs[Np] = {11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> checkPrimaries{"checkPrimaries", 1,
                                   "Whether to check physical primary and secondaries particles for the resolution."};
  Configurable<int> nBinsP{"nBinsP", 2000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20.f, "Maximum momentum in range"};
  Configurable<int> nBinsNsigma{"nBinsNsigma", 2000, "Number of bins for the momentum"};
  Configurable<float> minNsigma{"minNsigma", -30.f, "Minimum momentum in range"};
  Configurable<float> maxNsigma{"maxNsigma", 30.f, "Maximum momentum in range"};
  Configurable<float> minEta{"minEta", -0.8, "Minimum eta in range"};
  Configurable<float> maxEta{"maxEta", 0.8, "Maximum eta in range"};
  Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the vertex"};
  Configurable<int> logAxis{"logAxis", 0, "Flag to use a logarithmic pT axis, in this case the pT limits are the expontents"};

  template <uint8_t i>
  void addParticleHistos()
  {

    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaAxis{nBinsNsigma, minNsigma, maxNsigma, Form("N_{#sigma}^{TPC}(%s)", pT[pid_type])};
    if (logAxis) {
      ptAxis.makeLogaritmic();
    }

    // NSigma
    histos.add(hnsigmaMC[i].data(), Form("True %s", pT[i]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    if (!checkPrimaries) {
      return;
    }
    histos.add(hnsigmaMCprm[i].data(), Form("True Primary %s", pT[i]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmaMCsec[i].data(), Form("True Secondary %s", pT[i]), HistType::kTH2F, {ptAxis, nSigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogaritmic();
      ptAxis.makeLogaritmic();
    }
    const AxisSpec nSigmaAxis{nBinsNsigma, minNsigma, maxNsigma, Form("N_{#sigma}^{TPC}(%s)", pT[pid_type])};
    const AxisSpec signalAxis{6000, 0, 2000, "TPC d#it{E}/d#it{x} A.U."};

    histos.add(hnsigma[pid_type].data(), pT[pid_type], HistType::kTH2F, {ptAxis, nSigmaAxis});
    if (checkPrimaries) {
      histos.add(hnsigmaprm[pid_type].data(), Form("Primary %s", pT[pid_type]), HistType::kTH2F, {ptAxis, nSigmaAxis});
      histos.add(hnsigmasec[pid_type].data(), Form("Secondary %s", pT[pid_type]), HistType::kTH2F, {ptAxis, nSigmaAxis});
    }
    const AxisSpec lengthAxis{1000, 0, 3000, "Track length (cm)"};
    const AxisSpec etaAxis{100, -4, 4, "#it{#eta}"};

    histos.add("event/vertexz", ";Vtx_{z} (cm);Entries", kTH1F, {{100, -20, 20}});
    histos.add("particle/p", "", kTH1F, {pAxis});
    histos.add("particle/pt", "", kTH1F, {ptAxis});
    histos.add("particle/eta", "", kTH1F, {etaAxis});
    histos.add("tracks/p", "", kTH1F, {pAxis});
    histos.add("tracks/pt", "", kTH1F, {ptAxis});
    histos.add("tracks/eta", "", kTH1F, {etaAxis});
    histos.add("tracks/length", "", kTH1F, {lengthAxis});

    addParticleHistos<0>();
    addParticleHistos<1>();
    addParticleHistos<2>();
    addParticleHistos<3>();
    addParticleHistos<4>();
    addParticleHistos<5>();
    addParticleHistos<6>();
    addParticleHistos<7>();
    addParticleHistos<8>();
    histos.add("event/tpcsignal", "All", HistType::kTH2F, {pAxis, signalAxis});
    histos.add("event/tpcsignalMC", pT[pid_type], HistType::kTH2F, {pAxis, signalAxis});
    if (checkPrimaries) {
      histos.add("event/tpcsignalPrm", "Primaries", HistType::kTH2F, {pAxis, signalAxis});
      histos.add("event/tpcsignalSec", "Secondaries", HistType::kTH2F, {pAxis, signalAxis});
      histos.add("event/tpcsignalMCPrm", Form("Primary %s", pT[pid_type]), HistType::kTH2F, {pAxis, signalAxis});
      histos.add("event/tpcsignalMCSec", Form("Secondary %s", pT[pid_type]), HistType::kTH2F, {pAxis, signalAxis});
    }
  }

  template <uint8_t pidIndex, typename T>
  void fillNsigma(const T& track, const float& nsigma)
  {
    const auto particle = track.mcParticle();
    if (abs(particle.pdgCode()) == PDGs[pidIndex]) {

      histos.fill(HIST(hnsigmaMC[pidIndex]), track.pt(), nsigma);
      // Selecting primaries
      if (particle.isPhysicalPrimary()) {
        histos.fill(HIST(hnsigmaMCprm[pidIndex]), track.pt(), nsigma);
      } else {
        histos.fill(HIST(hnsigmaMCsec[pidIndex]), track.pt(), nsigma);
      }
    }
  }

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
      const auto tracksInCollision = tracks.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex());
      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex());

      for (const auto& p : particlesInCollision) {
        histos.fill(HIST("particle/p"), p.p());
        histos.fill(HIST("particle/pt"), p.pt());
        histos.fill(HIST("particle/eta"), p.eta());
      }

      for (const auto& t : tracksInCollision) {
        //
        if (t.eta() < minEta || t.eta() > maxEta) {
          continue;
        }

        histos.fill(HIST("tracks/p"), t.p());
        histos.fill(HIST("tracks/pt"), t.pt());
        histos.fill(HIST("tracks/eta"), t.eta());
        histos.fill(HIST("tracks/length"), t.length());

        const float nsigma = o2::aod::pidutils::tpcNSigma<pid_type>(t);

        // Fill for all
        histos.fill(HIST(hnsigma[pid_type]), t.pt(), nsigma);
        histos.fill(HIST("event/tpcsignal"), t.p(), t.tpcSignal());
        const auto particle = t.mcParticle();
        if (particle.isPhysicalPrimary()) { // Selecting primaries
          histos.fill(HIST(hnsigmaprm[pid_type]), t.pt(), nsigma);
          histos.fill(HIST("event/tpcsignalPrm"), t.p(), t.tpcSignal());
        } else {
          histos.fill(HIST(hnsigmasec[pid_type]), t.pt(), nsigma);
          histos.fill(HIST("event/tpcsignalSec"), t.p(), t.tpcSignal());
        }
        if (abs(particle.pdgCode()) == PDGs[pid_type]) { // Checking the PDG code
          histos.fill(HIST("event/tpcsignalMC"), t.pt(), t.tpcSignal());
          if (particle.isPhysicalPrimary()) {
            histos.fill(HIST("event/tpcsignalMCPrm"), t.pt(), t.tpcSignal());
          } else {
            histos.fill(HIST("event/tpcsignalMCSec"), t.pt(), t.tpcSignal());
          }
        }
        // Fill with PDG codes
        fillNsigma<0>(t, nsigma);
        fillNsigma<1>(t, nsigma);
        fillNsigma<2>(t, nsigma);
        fillNsigma<3>(t, nsigma);
        fillNsigma<4>(t, nsigma);
        fillNsigma<5>(t, nsigma);
        fillNsigma<6>(t, nsigma);
        fillNsigma<7>(t, nsigma);
        fillNsigma<8>(t, nsigma);
      }
      histos.fill(HIST("event/vertexz"), collision.posZ());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{};
  if (cfgc.options().get<int>("qa-el")) {
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Electron>>(cfgc, TaskName{"pidTPC-qa-mc-El"}));
  }
  if (cfgc.options().get<int>("qa-mu")) {
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Muon>>(cfgc, TaskName{"pidTPC-qa-mc-Mu"}));
  }
  if (cfgc.options().get<int>("qa-pikapr")) {
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Pion>>(cfgc, TaskName{"pidTPC-qa-mc-Pi"}));
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Kaon>>(cfgc, TaskName{"pidTPC-qa-mc-Ka"}));
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Proton>>(cfgc, TaskName{"pidTPC-qa-mc-Pr"}));
  }
  if (cfgc.options().get<int>("qa-nuclei")) {
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Deuteron>>(cfgc, TaskName{"pidTPC-qa-mc-De"}));
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Triton>>(cfgc, TaskName{"pidTPC-qa-mc-Tr"}));
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Helium3>>(cfgc, TaskName{"pidTPC-qa-mc-He"}));
    workflow.push_back(adaptAnalysisTask<pidTPCTaskQA<PID::Alpha>>(cfgc, TaskName{"pidTPC-qa-mc-Al"}));
  }
  return workflow;
}
