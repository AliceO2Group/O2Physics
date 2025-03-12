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
/// \file   alice3-pidTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce PID tables for TOF split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/PIDResponse.h"
#include "ALICE3/Core/TOFResoALICE3.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TOF PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

struct ALICE3tofSignal { /// Task that produces the TOF signal from the trackTime
  Produces<o2::aod::TOFSignal> table;
  bool enableTable = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        const std::string table = "TOFSignal";
        if (input.matcher.binding == table) {
          enableTable = true;
        }
      }
    }
  }
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
  void process(Trks const& tracks)
  {
    if (!enableTable) {
      return;
    }
    table.reserve(tracks.size());
    for (auto& t : tracks) {
      table(t.trackTime() * 1000.f);
    }
  }
};

struct ALICE3pidTOFTask {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
  using Coll = aod::Collisions;
  // Tables to produce
  Produces<o2::aod::pidTOFFullEl> tablePIDEl;
  Produces<o2::aod::pidTOFFullMu> tablePIDMu;
  Produces<o2::aod::pidTOFFullPi> tablePIDPi;
  Produces<o2::aod::pidTOFFullKa> tablePIDKa;
  Produces<o2::aod::pidTOFFullPr> tablePIDPr;
  Produces<o2::aod::pidTOFFullDe> tablePIDDe;
  Produces<o2::aod::pidTOFFullTr> tablePIDTr;
  Produces<o2::aod::pidTOFFullHe> tablePIDHe;
  Produces<o2::aod::pidTOFFullAl> tablePIDAl;
  Parameters resoParameters{1};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if empty the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFResoALICE3", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdb->setTimestamp(timestamp.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    //
    const std::vector<float> p = {24.5};
    resoParameters.SetParameters(p);
    const std::string fname = paramfile.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading parametrization from file" << fname << ", using param: " << sigmaname;
      resoParameters.LoadParamFromFile(fname.data(), sigmaname.value.data());
    } else { // Loading it from CCDB
      const std::string path = "Analysis/ALICE3/PID/TOF/Parameters";
      resoParameters.SetParameters(ccdb->getForTimeStamp<Parameters>(path + "/" + sigmaname.value, timestamp.value));
    }
  }

  template <o2::track::PID::ID id>
  float sigma(Trks::iterator track)
  {
    return o2::pid::tof::TOFResoALICE3ParamTrack<id>(track, resoParameters);
  }
  template <o2::track::PID::ID id>
  float nsigma(Trks::iterator track)
  {
    if (!track.hasTOF()) {
      return -999.f;
    }
    return ((track.trackTime() - track.collision().collisionTime()) * 1000.f - o2::pid::tof::ExpTimes<Trks::iterator, id>::ComputeExpectedTime(track.tofExpMom() * o2::constants::physics::invLightSpeedCm2PS,
                                                                                                                                               track.length())) /
           sigma<id>(track);
  }
  void process(Trks const& tracks, Coll const&)
  {
    tablePIDEl.reserve(tracks.size());
    tablePIDMu.reserve(tracks.size());
    tablePIDPi.reserve(tracks.size());
    tablePIDKa.reserve(tracks.size());
    tablePIDPr.reserve(tracks.size());
    tablePIDDe.reserve(tracks.size());
    tablePIDTr.reserve(tracks.size());
    tablePIDHe.reserve(tracks.size());
    tablePIDAl.reserve(tracks.size());
    for (auto const& trk : tracks) {
      tablePIDEl(sigma<PID::Electron>(trk), nsigma<PID::Electron>(trk));
      tablePIDMu(sigma<PID::Muon>(trk), nsigma<PID::Muon>(trk));
      tablePIDPi(sigma<PID::Pion>(trk), nsigma<PID::Pion>(trk));
      tablePIDKa(sigma<PID::Kaon>(trk), nsigma<PID::Kaon>(trk));
      tablePIDPr(sigma<PID::Proton>(trk), nsigma<PID::Proton>(trk));
      tablePIDDe(sigma<PID::Deuteron>(trk), nsigma<PID::Deuteron>(trk));
      tablePIDTr(sigma<PID::Triton>(trk), nsigma<PID::Triton>(trk));
      tablePIDHe(sigma<PID::Helium3>(trk), nsigma<PID::Helium3>(trk));
      tablePIDAl(sigma<PID::Alpha>(trk), nsigma<PID::Alpha>(trk));
    }
  }
};

struct ALICE3pidTOFTaskQA {

  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hexpected_diff[Np] = {"expected_diff/El", "expected_diff/Mu", "expected_diff/Pi",
                                                          "expected_diff/Ka", "expected_diff/Pr", "expected_diff/De",
                                                          "expected_diff/Tr", "expected_diff/He", "expected_diff/Al"};
  static constexpr std::string_view hexpsigma[Np] = {"expsigma/El", "expsigma/Mu", "expsigma/Pi",
                                                     "expsigma/Ka", "expsigma/Pr", "expsigma/De",
                                                     "expsigma/Tr", "expsigma/He", "expsigma/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20.f, "Maximum momentum in range"};
  Configurable<int> nBinsDelta{"nBinsDelta", 200, "Number of bins for the Delta"};
  Configurable<float> minDelta{"minDelta", -1000.f, "Minimum Delta in range"};
  Configurable<float> maxDelta{"maxDelta", 1000.f, "Maximum Delta in range"};
  Configurable<int> nBinsExpSigma{"nBinsExpSigma", 200, "Number of bins for the ExpSigma"};
  Configurable<float> minExpSigma{"minExpSigma", 0.f, "Minimum ExpSigma in range"};
  Configurable<float> maxExpSigma{"maxExpSigma", 200.f, "Maximum ExpSigma in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Maximum NSigma in range"};
  Configurable<int> minMult{"minMult", 1, "Minimum track multiplicity with TOF"};

  template <uint8_t i>
  void addParticleHistos(const AxisSpec& pAxis, const AxisSpec& /*ptAxis*/)
  {
    // Exp signal
    const AxisSpec expAxis{1000, 0, 2e6, Form("t_{exp}(%s)", pT[i])};
    histos.add(hexpected[i].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("(t-t_{evt}-t_{exp}(%s))", pT[i])};
    histos.add(hexpected_diff[i].data(), "", kTH2F, {pAxis, deltaAxis});

    // Exp Sigma
    const AxisSpec expSigmaAxis{nBinsExpSigma, minExpSigma, maxExpSigma, Form("Exp_{#sigma}^{TOF}(%s)", pT[i])};
    histos.add(hexpsigma[i].data(), "", kTH2F, {pAxis, expSigmaAxis});

    // NSigma
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, Form("N_{#sigma}^{TOF}(%s)", pT[i])};
    histos.add(hnsigma[i].data(), "", kTH2F, {pAxis, nSigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {

    const AxisSpec multAxis{100, 0, 100, "TOF multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec tofAxis{10000, 0, 2e6, "TOF Signal"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec colTimeAxis{100, -2000, 2000, "Collision time (ps)"};
    const AxisSpec colTimeResoAxis{100, 0, 1000, "#sigma_{Collision time} (ps)"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec ptResoAxis{100, 0, 0.1, "#sigma_{#it{p}_{T}}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    AxisSpec pExpAxis{nBinsP, minP, maxP, "#it{p}_{Exp. TOF} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogarithmic();
      pAxis.makeLogarithmic();
      pExpAxis.makeLogarithmic();
    }

    // Event properties
    histos.add("event/vertexz", "", kTH1F, {vtxZAxis});
    histos.add("event/tofmultiplicity", "", kTH1F, {multAxis});
    histos.add("event/colltime", "", kTH1F, {colTimeAxis});
    histos.add("event/colltimereso", "", kTH2F, {multAxis, colTimeResoAxis});
    histos.add("event/tofsignal", "", kTH2F, {pAxis, tofAxis});
    histos.add("event/pexp", "", kTH2F, {pAxis, pExpAxis});
    histos.add("event/eta", "", kTH1F, {etaAxis});
    histos.add("event/length", "", kTH1F, {lAxis});
    histos.add("event/pt", "", kTH1F, {ptAxis});
    histos.add("event/p", "", kTH1F, {pAxis});
    histos.add("event/ptreso", "", kTH2F, {pAxis, ptResoAxis});

    static_for<0, 8>([&](auto i) {
      addParticleHistos<i>(pAxis, ptAxis);
    });
  }

  template <o2::track::PID::ID i, typename T>
  void fillParticleHistos(const T& t, const float& tof, const float& exp_diff, const float& expsigma)
  {
    histos.fill(HIST(hexpected[i]), t.p(), tof - exp_diff);
    histos.fill(HIST(hexpected_diff[i]), t.p(), exp_diff);
    histos.fill(HIST(hexpsigma[i]), t.p(), expsigma);
    histos.fill(HIST(hnsigma[i]), t.p(), o2::aod::pidutils::tofNSigma(i, t));
  }

  void process(aod::Collision const& collision,
               soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                         aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                         aod::TrackSelection> const& tracks)
  {
    // Computing Multiplicity first
    int mult = 0;
    for (auto t : tracks) {
      //
      if (!t.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      mult++;
    }
    if (mult < minMult) { // Cutting on low multiplicity events
      return;
    }

    const float collisionTime_ps = collision.collisionTime() * 1000.f;
    histos.fill(HIST("event/vertexz"), collision.posZ());
    histos.fill(HIST("event/colltime"), collisionTime_ps);
    histos.fill(HIST("event/tofmultiplicity"), mult);
    histos.fill(HIST("event/colltimereso"), mult, collision.collisionTimeRes() * 1000.f);

    for (auto t : tracks) {
      //
      if (!t.hasTOF()) { // Skipping tracks without TOF
        continue;
      }
      if (!t.isGlobalTrack()) {
        continue;
      }

      const float tofSignal = t.trackTime() * 1e3f;
      // const float tof = t.tofSignal() - collisionTime_ps;
      const float tof = tofSignal - collisionTime_ps;

      //
      // histos.fill(HIST("event/tofsignal"), t.p(), t.tofSignal());
      histos.fill(HIST("event/tofsignal"), t.p(), tofSignal);
      histos.fill(HIST("event/pexp"), t.p(), t.tofExpMom() * o2::constants::physics::invLightSpeedCm2PS);
      histos.fill(HIST("event/eta"), t.eta());
      histos.fill(HIST("event/length"), t.length());
      histos.fill(HIST("event/pt"), t.pt());
      histos.fill(HIST("event/ptreso"), t.p(), t.sigma1Pt() * t.pt() * t.pt());
      //
      fillParticleHistos<o2::track::PID::Electron>(t, tof, t.tofExpSignalDiffEl(), t.tofExpSigmaEl());
      fillParticleHistos<o2::track::PID::Muon>(t, tof, t.tofExpSignalDiffMu(), t.tofExpSigmaMu());
      fillParticleHistos<o2::track::PID::Pion>(t, tof, t.tofExpSignalDiffPi(), t.tofExpSigmaPi());
      fillParticleHistos<o2::track::PID::Kaon>(t, tof, t.tofExpSignalDiffKa(), t.tofExpSigmaKa());
      fillParticleHistos<o2::track::PID::Proton>(t, tof, t.tofExpSignalDiffPr(), t.tofExpSigmaPr());
      fillParticleHistos<o2::track::PID::Deuteron>(t, tof, t.tofExpSignalDiffDe(), t.tofExpSigmaDe());
      fillParticleHistos<o2::track::PID::Triton>(t, tof, t.tofExpSignalDiffTr(), t.tofExpSigmaTr());
      fillParticleHistos<o2::track::PID::Helium3>(t, tof, t.tofExpSignalDiffHe(), t.tofExpSigmaHe());
      fillParticleHistos<o2::track::PID::Alpha>(t, tof, t.tofExpSignalDiffAl(), t.tofExpSigmaAl());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<ALICE3pidTOFTask>(cfgc, TaskName{"alice3-pidTOF-task"}),
                               adaptAnalysisTask<ALICE3tofSignal>(cfgc, TaskName{"alice3-tof-signal-task"})};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<ALICE3pidTOFTaskQA>(cfgc, TaskName{"alice3-pidTOFQA-task"}));
  }
  return workflow;
}
