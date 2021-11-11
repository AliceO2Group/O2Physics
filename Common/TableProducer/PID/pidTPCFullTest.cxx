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
/// \file   pidTPCFullTest.cxx
/// \author Annalena Kalteyer, Christian Sonnabend
/// \brief  Task to produce PID tables for TPC split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

#include "TFile.h"

//  // O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TPC PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

struct tpcPidFull {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  using Coll = aod::Collisions;
  // Tables to produce
  Produces<o2::aod::pidTPCFullEl> tablePIDEl;
  Produces<o2::aod::pidTPCFullMu> tablePIDMu;
  Produces<o2::aod::pidTPCFullPi> tablePIDPi;
  Produces<o2::aod::pidTPCFullKa> tablePIDKa;
  Produces<o2::aod::pidTPCFullPr> tablePIDPr;
  Produces<o2::aod::pidTPCFullDe> tablePIDDe;
  Produces<o2::aod::pidTPCFullTr> tablePIDTr;
  Produces<o2::aod::pidTPCFullHe> tablePIDHe;
  Produces<o2::aod::pidTPCFullAl> tablePIDAl;
  // TPC PID Response
  o2::pid::tpc::Response* response = nullptr;
  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<long> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<int> pidEl{"pid-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidMu{"pid-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPi{"pid-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidKa{"pid-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPr{"pid-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidDe{"pid-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTr{"pid-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidHe{"pid-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidAl{"pid-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};

  void init(o2::framework::InitContext& initContext)
  {
    // Checking the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec device : workflows.devices) {
      for (auto input : device.inputs) {
        auto enableFlag = [&input](const std::string particle, Configurable<int>& flag) {
          const std::string table = "pidTPCFull" + particle;
          if (input.matcher.binding == table) {
            if (flag < 0) {
              flag.value = 1;
              LOG(info) << "Auto-enabling table: " + table;
            } else if (flag > 0) {
              flag.value = 1;
              LOG(info) << "Table enabled: " + table;
            } else {
              LOG(info) << "Table disabled: " + table;
            }
          }
        };
        enableFlag("El", pidEl);
        enableFlag("Mu", pidMu);
        enableFlag("Pi", pidPi);
        enableFlag("Ka", pidKa);
        enableFlag("Pr", pidPr);
        enableFlag("De", pidDe);
        enableFlag("Tr", pidTr);
        enableFlag("He", pidHe);
        enableFlag("Al", pidAl);
      }
    }

    const TString fname = paramfile.value;
    if (fname != "") { // Loading the parametrization from file
      LOGP(info, "Loading TPC response from file {}", fname);
      try {
        std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
        f->GetObject("Response", response);
      } catch (...) {
        LOGP(info, "Loading the TPC PID Response from file {} failed!", fname);
      };
    } else {
      const std::string path = ccdbPath.value;
      const auto time = timestamp.value;
      ccdb->setURL(url.value);
      ccdb->setTimestamp(time);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, time);
      LOGP(info, "Loading TPC response from CCDB, using path: {} for timestamp {}", path, time);
    }
  }

  void process(Coll const& collisions, Trks const& tracks)
  {
    // LOG(info) << "Custom TPCPID: MIP value from object is " << response->GetMIP();

    // Check and fill enabled tables
    auto makeTable = [&tracks](const Configurable<int>& flag, auto& table, const o2::pid::tpc::Response& response, const o2::track::PID::ID pid) {
      if (flag.value == 1) {
        // Prepare memory for enabled tables
        table.reserve(tracks.size());
        for (auto const& trk : tracks) { // Loop on Tracks
          table(response.GetExpectedSigma(trk, pid),
                response.GetNumberOfSigma(trk, pid));
        }
      }
    };
    //const o2::pid::tpc::Response& response;
    makeTable(pidEl, tablePIDEl, *response, o2::track::PID::Electron);
    makeTable(pidMu, tablePIDMu, *response, o2::track::PID::Muon);
    makeTable(pidPi, tablePIDPi, *response, o2::track::PID::Pion);
    makeTable(pidKa, tablePIDKa, *response, o2::track::PID::Kaon);
    makeTable(pidPr, tablePIDPr, *response, o2::track::PID::Proton);
    makeTable(pidDe, tablePIDDe, *response, o2::track::PID::Deuteron);
    makeTable(pidTr, tablePIDTr, *response, o2::track::PID::Triton);
    makeTable(pidHe, tablePIDHe, *response, o2::track::PID::Helium3);
    makeTable(pidAl, tablePIDAl, *response, o2::track::PID::Alpha);
  }
};

struct tpcPidFullQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hexpected[Np] = {"expected/El", "expected/Mu", "expected/Pi",
                                                     "expected/Ka", "expected/Pr", "expected/De",
                                                     "expected/Tr", "expected/He", "expected/Al"};
  static constexpr std::string_view hexpected_diff[Np] = {"expected_diff/El", "expected_diff/Mu", "expected_diff/Pi",
                                                          "expected_diff/Ka", "expected_diff/Pr", "expected_diff/De",
                                                          "expected_diff/Tr", "expected_diff/He", "expected_diff/Al"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 400, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> nBinsDelta{"nBinsDelta", 200, "Number of bins for the Delta"};
  Configurable<float> minDelta{"minDelta", -1000.f, "Minimum Delta in range"};
  Configurable<float> maxDelta{"maxDelta", 1000.f, "Maximum Delta in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Maximum NSigma in range"};

  template <uint8_t i>
  void addParticleHistos()
  {
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogaritmic();
    }

    // Exp signal
    const AxisSpec expAxis{1000, 0, 1000, Form("d#it{E}/d#it{x}_(%s) A.U.", pT[i])};
    histos.add(hexpected[i].data(), "", kTH2F, {pAxis, expAxis});

    // Signal - Expected signal
    const AxisSpec deltaAxis{nBinsDelta, minDelta, maxDelta, Form("d#it{E}/d#it{x} - d#it{E}/d#it{x}(%s)", pT[i])};
    histos.add(hexpected_diff[i].data(), "", kTH2F, {pAxis, deltaAxis});

    // NSigma
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, Form("N_{#sigma}^{TPC}(%s)", pT[i])};
    histos.add(hnsigma[i].data(), "", kTH2F, {pAxis, nSigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {

    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogaritmic();
    }
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec dedxAxis{1000, 0, 1000, "d#it{E}/d#it{x} A.U."};

    // Event properties
    histos.add("event/vertexz", "", kTH1F, {vtxZAxis});
    histos.add("event/tpcsignal", "", kTH2F, {pAxis, dedxAxis});

    addParticleHistos<0>();
    addParticleHistos<1>();
    addParticleHistos<2>();
    addParticleHistos<3>();
    addParticleHistos<4>();
    addParticleHistos<5>();
    addParticleHistos<6>();
    addParticleHistos<7>();
    addParticleHistos<8>();
  }

  template <uint8_t i, typename T>
  void fillParticleHistos(const T& t, const float mom, const float exp_diff, const float nsigma)
  {
    histos.fill(HIST(hexpected[i]), mom, t.tpcSignal() - exp_diff);
    histos.fill(HIST(hexpected_diff[i]), mom, exp_diff);
    histos.fill(HIST(hnsigma[i]), t.p(), nsigma);
  }

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra,
                                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                          aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                                          aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                                                          aod::TrackSelection> const& tracks)
  {
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (auto t : tracks) {
      // const float mom = t.p();
      const float mom = t.tpcInnerParam();
      histos.fill(HIST("event/tpcsignal"), mom, t.tpcSignal());
      //
      fillParticleHistos<0>(t, mom, t.tpcExpSignalDiffEl(), t.tpcNSigmaEl());
      fillParticleHistos<1>(t, mom, t.tpcExpSignalDiffMu(), t.tpcNSigmaMu());
      fillParticleHistos<2>(t, mom, t.tpcExpSignalDiffPi(), t.tpcNSigmaPi());
      fillParticleHistos<3>(t, mom, t.tpcExpSignalDiffKa(), t.tpcNSigmaKa());
      fillParticleHistos<4>(t, mom, t.tpcExpSignalDiffPr(), t.tpcNSigmaPr());
      fillParticleHistos<5>(t, mom, t.tpcExpSignalDiffDe(), t.tpcNSigmaDe());
      fillParticleHistos<6>(t, mom, t.tpcExpSignalDiffTr(), t.tpcNSigmaTr());
      fillParticleHistos<7>(t, mom, t.tpcExpSignalDiffHe(), t.tpcNSigmaHe());
      fillParticleHistos<8>(t, mom, t.tpcExpSignalDiffAl(), t.tpcNSigmaAl());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tpcPidFull>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tpcPidFullQa>(cfgc));
  }
  return workflow;
}