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
/// \file   pidTPC.cxx
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce PID tables for TPC split for each particle with only the Nsigma information.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

#include "TFile.h"

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/StaticFor.h"
#include "Common/TableProducer/PID/pidTPCML.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::pid::tpc;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TPC PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the response table
struct tpcPid {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  using Coll = soa::Join<aod::Collisions, aod::Mults>;
  // Tables to produce
  Produces<o2::aod::pidTPCEl> tablePIDEl;
  Produces<o2::aod::pidTPCMu> tablePIDMu;
  Produces<o2::aod::pidTPCPi> tablePIDPi;
  Produces<o2::aod::pidTPCKa> tablePIDKa;
  Produces<o2::aod::pidTPCPr> tablePIDPr;
  Produces<o2::aod::pidTPCDe> tablePIDDe;
  Produces<o2::aod::pidTPCTr> tablePIDTr;
  Produces<o2::aod::pidTPCHe> tablePIDHe;
  Produces<o2::aod::pidTPCAl> tablePIDAl;
  // TPC PID Response
  o2::pid::tpc::Response response;
  o2::pid::tpc::Response* responseptr = nullptr;
  // Network correction for TPC PID response
  Network network;

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<long> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};
  // Parameters for loading network from a file / downloading the file
  Configurable<int> useNetworkCorrection{"useNetworkCorrection", 0, "Using the network correction for the TPC dE/dx signal"};
  Configurable<int> downloadNetworkFromAlien{"downloadNetworkFromAlien", 0, "Download network from AliEn (1) or use a local file (filepath must be provided by --networkPathLocally /path/to/file) (0)"};
  Configurable<std::string> networkPathAlien{"networkPathAlien", "alien:///alice/cern.ch/user/c/csonnabe/tpc_network_testing/net_onnx_0.onnx", "Path to .onnx file containing the network on AliEn"};
  Configurable<std::string> networkPathLocally{"networkPathLocally", "network.onnx", "Path to local .onnx file containing the network"};
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

  // Paramatrization configuration
  bool useCCDBParam = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      const std::string table = "pidTPC" + particle;
      if (isTableRequiredInWorkflow(initContext, table)) {
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

    const TString fname = paramfile.value;
    if (fname != "") { // Loading the parametrization from file
      LOGP(info, "Loading TPC response from file {}", fname);
      try {
        std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
        f->GetObject("Response", responseptr);
        response.SetParameters(responseptr);
      } catch (...) {
        LOGF(fatal, "Loading the TPC PID Response from file {} failed!", fname);
      };
    } else {
      useCCDBParam = true;
      const std::string path = ccdbPath.value;
      const auto time = ccdbTimestamp.value;
      ccdb->setURL(url.value);
      ccdb->setTimestamp(time);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      response.SetParameters(ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, time));
      LOGP(info, "Loading TPC response from CCDB, using path: {} for ccdbTimestamp {}", path, time);
      response.PrintAll();
    }

    if (!useNetworkCorrection) {
      return;
    } else {
      Network temp_net(networkPathLocally.value,
                       downloadNetworkFromAlien.value,
                       networkPathAlien.value,
                       true);
      network = temp_net;
    }
  }

  void process(Coll const& collisions, Trks const& tracks,
               aod::BCsWithTimestamps const&)
  {
    auto reserveTable = [&tracks](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(tracks.size());
    };
    // Prepare memory for enabled tables
    reserveTable(pidEl, tablePIDEl);
    reserveTable(pidMu, tablePIDMu);
    reserveTable(pidPi, tablePIDPi);
    reserveTable(pidKa, tablePIDKa);
    reserveTable(pidPr, tablePIDPr);
    reserveTable(pidDe, tablePIDDe);
    reserveTable(pidTr, tablePIDTr);
    reserveTable(pidHe, tablePIDHe);
    reserveTable(pidAl, tablePIDAl);

    std::vector<float> network_prediction;

    if (useNetworkCorrection) {

      auto start_overhead = std::chrono::high_resolution_clock::now();
      std::vector<float> track_properties;
      for (int i = 0; i < 9; i++) { // Loop over particle number for which network correction is used
        for (auto const& trk : tracks) {
          std::vector<float> net_tensor = network.createInputFromTrack(trk, i);
          for (auto value : net_tensor) {
            track_properties.push_back(value);
          }
        }
      }
      const unsigned long track_prop_size = tracks.size() * 9;
      auto stop_overhead = std::chrono::high_resolution_clock::now();
      float duration_overhead = std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_overhead - start_overhead).count();
      float time_per_track_overhead = duration_overhead / track_prop_size; // There are n (typically n=7) variables in each track which are being extracted in track_properties. Each network evaluation takes time_per_track_overhead/9 nano-seconds
      LOG(info) << "Time per track (overhead): " << time_per_track_overhead << "ns ; Overhead total: " << duration_overhead / 1000000000 << "s";

      auto start_network = std::chrono::high_resolution_clock::now();
      float* output_network = network.evalNetwork(track_properties);
      for (unsigned long i = 0; i < track_prop_size; i++) {
        network_prediction.push_back(output_network[i]);
      }
      track_properties.clear();
      auto stop_network = std::chrono::high_resolution_clock::now();
      float duration_network = std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network - start_network).count();
      float time_per_track_net = duration_network / track_prop_size;
      LOG(info) << "Time per track (net): " << time_per_track_net << "ns ; Network total: " << duration_network / 1000000000 << "s"; // The time per track but with 9 particle mass hypotheses: So actual time per track is (time_per_track_net / 9)

      // Uncomment if you want to check example-outputs of the netwwork:
      // for(int i=0; i<100; i++){
      //   LOG(info) << "Output " << i << ": " << network_prediction[i] << " ; Input: [" << track_properties[7*i + 0] << ", " << track_properties[7*i + 1] << ", " << track_properties[7*i + 2] << ", " << track_properties[7*i + 3] << ", " << track_properties[7*i + 4] << ", " << track_properties[7*i + 5] << ", " << track_properties[7*i + 6] << "]";
      // }
    }

    int lastCollisionId = -1; // Last collision ID analysed
    unsigned long count_tracks = 0;
    const int tracks_size = tracks.size();

    for (auto const& trk : tracks) {                                                                                 // Loop on Tracks
      if (useCCDBParam && ccdbTimestamp.value == 0 && trk.has_collision() && trk.collisionId() != lastCollisionId) { // Updating parametrization only if the initial timestamp is 0
        lastCollisionId = trk.collisionId();
        const auto& bc = collisions.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>();
        response.SetParameters(ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp()));
      }
      // Check and fill enabled tables
      auto makeTable = [&trk, &collisions, &network_prediction, &count_tracks, &tracks_size, this](const Configurable<int>& flag, auto& table, const o2::track::PID::ID pid) {
        if (flag.value != 1) {
          return;
        }

        if (useNetworkCorrection) {
          aod::pidutils::packInTable<aod::pidtpc_tiny::binning>((trk.tpcSignal() - (network_prediction[count_tracks + tracks_size * pid]) * response.GetExpectedSignal(trk, pid)) / response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid), table);
        } else {
          aod::pidutils::packInTable<aod::pidtpc_tiny::binning>(response.GetNumberOfSigma(collisions.iteratorAt(trk.collisionId()), trk, pid), table);
        }
      };

      makeTable(pidEl, tablePIDEl, o2::track::PID::Electron);
      makeTable(pidMu, tablePIDMu, o2::track::PID::Muon);
      makeTable(pidPi, tablePIDPi, o2::track::PID::Pion);
      makeTable(pidKa, tablePIDKa, o2::track::PID::Kaon);
      makeTable(pidPr, tablePIDPr, o2::track::PID::Proton);
      makeTable(pidDe, tablePIDDe, o2::track::PID::Deuteron);
      makeTable(pidTr, tablePIDTr, o2::track::PID::Triton);
      makeTable(pidHe, tablePIDHe, o2::track::PID::Helium3);
      makeTable(pidAl, tablePIDAl, o2::track::PID::Alpha);

      count_tracks++;
    }
  }
};

/// Task to produce the TPC QA plots
struct tpcPidQa {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hnsigma[Np] = {"nsigma/El", "nsigma/Mu", "nsigma/Pi",
                                                   "nsigma/Ka", "nsigma/Pr", "nsigma/De",
                                                   "nsigma/Tr", "nsigma/He", "nsigma/Al"};
  static constexpr std::string_view hnsigmapt[Np] = {"nsigmapt/El", "nsigmapt/Mu", "nsigmapt/Pi",
                                                     "nsigmapt/Ka", "nsigmapt/Pr", "nsigmapt/De",
                                                     "nsigmapt/Tr", "nsigmapt/He", "nsigmapt/Al"};
  static constexpr std::string_view hnsigmapospt[Np] = {"nsigmapospt/El", "nsigmapospt/Mu", "nsigmapospt/Pi",
                                                        "nsigmapospt/Ka", "nsigmapospt/Pr", "nsigmapospt/De",
                                                        "nsigmapospt/Tr", "nsigmapospt/He", "nsigmapospt/Al"};
  static constexpr std::string_view hnsigmanegpt[Np] = {"nsigmanegpt/El", "nsigmanegpt/Mu", "nsigmanegpt/Pi",
                                                        "nsigmanegpt/Ka", "nsigmanegpt/Pr", "nsigmanegpt/De",
                                                        "nsigmanegpt/Tr", "nsigmanegpt/He", "nsigmanegpt/Al"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 200, "Number of bins for the NSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Minimum NSigma in range"};
  Configurable<float> maxNSigma{"maxNSigma", 10.f, "Maximum NSigma in range"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<bool> applyTrackCut{"applyTrackCut", false, "Flag to apply standard track cuts"};
  Configurable<bool> applyRapidityCut{"applyRapidityCut", false, "Flag to apply rapidity cut"};

  template <uint8_t i>
  void addParticleHistos(const AxisSpec& pAxis, const AxisSpec& ptAxis)
  {
    // NSigma
    const char* axisTitle = Form("N_{#sigma}^{TPC}(%s)", pT[i]);
    const AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, axisTitle};
    histos.add(hnsigma[i].data(), axisTitle, kTH2F, {pAxis, nSigmaAxis});
    histos.add(hnsigmapt[i].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmapospt[i].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
    histos.add(hnsigmanegpt[i].data(), axisTitle, kTH2F, {ptAxis, nSigmaAxis});
  }

  void init(o2::framework::InitContext&)
  {

    const AxisSpec multAxis{1000, 0.f, 1000.f, "Track multiplicity"};
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec etaAxis{100, -2, 2, "#it{#eta}"};
    const AxisSpec lAxis{100, 0, 500, "Track length (cm)"};
    const AxisSpec pAxisPosNeg{nBinsP, -maxP, maxP, "Signed #it{p} (GeV/#it{c})"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    if (logAxis) {
      ptAxis.makeLogaritmic();
      pAxis.makeLogaritmic();
    }
    const AxisSpec dedxAxis{1000, 0, 1000, "d#it{E}/d#it{x} A.U."};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1F, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed mult.");
    h->GetXaxis()->SetBinLabel(4, "Passed vtx Z");

    histos.add("event/vertexz", "", kTH1F, {vtxZAxis});
    h = histos.add<TH1>("event/particlehypo", "", kTH1F, {{10, 0, 10, "PID in tracking"}});
    for (int i = 0; i < 9; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PID::getName(i));
    }
    histos.add("event/trackmultiplicity", "", kTH1F, {multAxis});
    histos.add("event/tpcsignal", "", kTH2F, {pAxis, dedxAxis});
    histos.add("event/signedtpcsignal", "", kTH2F, {pAxisPosNeg, dedxAxis});
    histos.add("event/eta", "", kTH1F, {etaAxis});
    histos.add("event/length", "", kTH1F, {lAxis});
    histos.add("event/pt", "", kTH1F, {ptAxis});
    histos.add("event/p", "", kTH1F, {pAxis});

    static_for<0, 8>([&](auto i) {
      addParticleHistos<i>(pAxis, ptAxis);
    });
  }

  template <o2::track::PID::ID id, typename T>
  void fillParticleHistos(const T& t)
  {
    if (applyRapidityCut) {
      const float y = TMath::ASinH(t.pt() / TMath::Sqrt(PID::getMass2(id) + t.pt() * t.pt()) * TMath::SinH(t.eta()));
      if (abs(y) > 0.5) {
        return;
      }
    }

    // Fill histograms
    const auto& nsigma = o2::aod::pidutils::tpcNSigma<id>(t);
    histos.fill(HIST(hnsigma[id]), t.p(), nsigma);
    histos.fill(HIST(hnsigmapt[id]), t.pt(), nsigma);
    if (t.sign() > 0) {
      histos.fill(HIST(hnsigmapospt[id]), t.pt(), nsigma);
    } else {
      histos.fill(HIST(hnsigmanegpt[id]), t.pt(), nsigma);
    }
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                         aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                         aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl,
                         aod::TrackSelection>;
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               Trks const& tracks)
  {

    histos.fill(HIST("event/evsel"), 1);
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return;
      }
    }

    histos.fill(HIST("event/evsel"), 2);

    // Computing Multiplicity first
    float ntracks = 0;
    for (auto t : tracks) {
      if (applyTrackCut && !t.isGlobalTrack()) {
        continue;
      }
      ntracks += 1;
    }
    // if (0 && ntracks < 1) {
    //   return;
    // }
    histos.fill(HIST("event/evsel"), 3);
    if (abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("event/evsel"), 4);
    histos.fill(HIST("event/vertexz"), collision.posZ());
    histos.fill(HIST("event/trackmultiplicity"), ntracks);

    for (auto t : tracks) {
      if (applyTrackCut && !t.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("event/particlehypo"), t.pidForTracking());
      histos.fill(HIST("event/tpcsignal"), t.tpcInnerParam(), t.tpcSignal());
      histos.fill(HIST("event/signedtpcsignal"), t.tpcInnerParam() * t.sign(), t.tpcSignal());
      histos.fill(HIST("event/eta"), t.eta());
      histos.fill(HIST("event/length"), t.length());
      histos.fill(HIST("event/pt"), t.pt());
      histos.fill(HIST("event/p"), t.p());
      //
      fillParticleHistos<PID::Electron>(t);
      fillParticleHistos<PID::Muon>(t);
      fillParticleHistos<PID::Pion>(t);
      fillParticleHistos<PID::Kaon>(t);
      fillParticleHistos<PID::Proton>(t);
      fillParticleHistos<PID::Deuteron>(t);
      fillParticleHistos<PID::Triton>(t);
      fillParticleHistos<PID::Helium3>(t);
      fillParticleHistos<PID::Alpha>(t);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tpcPid>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tpcPidQa>(cfgc));
  }
  return workflow;
}