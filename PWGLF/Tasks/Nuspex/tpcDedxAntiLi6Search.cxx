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

/// \file tpcDedxAntiLi6Search.cxx
/// \brief TPC dE/dx vs. signed momentum histograms (linear and logarithmic momentum bins)
///        and optional search for anti-lithium-6 candidates along a user-given Bethe-Bloch curve.
///
/// The candidates are written to a flat TTree (one branch per feature) together with all
/// markers needed to locate the corresponding clusters in the CTF files (run number, global BC,
/// orbit, timestamp, data frame ID, collision/track times, TPC sector, ...).
/// The width of the selection window can either come from the default parametrisation
/// (relative resolution * sqrt(NClNorm / NCl)) or from the neural network used in pidTPCModule.h,
/// evaluated with the proton mass hypothesis.
///
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/ML/model.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/PID.h>

#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct TpcDedxAntiLi6Search {
  // ---------------------------------------------------------------- histogram configuration
  Configurable<float> pMax{"pMax", 10.f, "Maximum |p/Z| (GeV/c) of the dE/dx histograms"};
  Configurable<float> pMinLog{"pMinLog", 0.1f, "Minimum |p/Z| (GeV/c) of the logarithmic binning"};
  Configurable<int> nBinsPLinear{"nBinsPLinear", 1000, "Number of p/Z bins (total, both signs) of the linear histogram"};
  Configurable<int> nBinsPLogPerSign{"nBinsPLogPerSign", 300, "Number of p/Z bins per charge sign of the logarithmic histogram"};
  Configurable<float> dEdxMax{"dEdxMax", 1300.f, "Maximum dE/dx (arb. units) of the histograms"};
  Configurable<int> nBinsDedx{"nBinsDedx", 650, "Number of dE/dx bins"};

  // ---------------------------------------------------------------- anti-Li6 search
  Configurable<bool> doSearch{"doSearch", false, "Enable the candidate extraction (see mode) and the output tree"};
  Configurable<std::string> mode{"mode", "post-bb-fit", "pre-bb-fit: extract Li6 candidates in [-preFitNSigmaBelow, +preFitNSigmaAbove] sigma of the standard curve (bbParamsStd) for the offline fit; post-bb-fit: select anti-Li6 with the fitted curve (bbParams)"};
  Configurable<std::vector<float>> bbParamsStd{"bbParamsStd", std::vector<float>{0.0320981f, 19.9768f, 2.52666e-16f, 2.72123f, 6.08092f}, "[pre-bb-fit] Fallback standard Bethe-Bloch (ALEPH) parameters p0..p4, only used if useCcdbStdParams=false"};
  Configurable<int> preFitSign{"preFitSign", 1, "[pre-bb-fit] Required charge sign of the Li6 candidates (1: Li6, -1: anti-Li6, 0: both)"};
  Configurable<float> preFitNSigmaBelow{"preFitNSigmaBelow", 5.f, "[pre-bb-fit] Candidates may lie at most this many sigma below the standard curve"};
  Configurable<float> preFitNSigmaAbove{"preFitNSigmaAbove", 10.f, "[pre-bb-fit] Candidates may lie at most this many sigma above the standard curve"};
  Configurable<std::vector<float>> bbParams{"bbParams", std::vector<float>{0.0320981f, 19.9768f, 2.52666e-16f, 2.72123f, 6.08092f}, "[post-bb-fit] Fitted Bethe-Bloch (ALEPH) parameters p0..p4 used for the anti-Li6 selection"};
  Configurable<bool> useCcdbStdParams{"useCcdbStdParams", true, "[pre-bb-fit] Take the standard Bethe-Bloch parameters, MIP, charge factor and NCl normalisation from the TPC response object on the CCDB (per run)"};
  Configurable<std::string> ccdbPathResponse{"ccdbPathResponse", "Analysis/PID/TPC/Response", "[pre-bb-fit] Path of the TPC response object on the CCDB"};
  Configurable<float> bbMip{"bbMip", 50.f, "dE/dx of a MIP (scaling of the Bethe-Bloch curve). post-bb-fit: used as is; pre-bb-fit: only if useCcdbStdParams=false"};
  Configurable<float> bbChargeFactor{"bbChargeFactor", 2.3f, "Exponent of the charge scaling Z^x of the Bethe-Bloch curve. post-bb-fit: used as is; pre-bb-fit: only if useCcdbStdParams=false"};
  Configurable<float> massLi6{"massLi6", 5.60152f, "Mass of Li6 (GeV/c^2)"};
  Configurable<float> chargeLi6{"chargeLi6", 3.f, "Charge number Z of Li6"};
  Configurable<bool> rigidityToMomentum{"rigidityToMomentum", true, "If true: beta*gamma = Z * tpcInnerParam / mass (tpcInnerParam is treated as rigidity p/Z), otherwise tpcInnerParam / mass"};
  Configurable<float> nSigmaBelow{"nSigmaBelow", 5.f, "[post-bb-fit] Candidates may lie at most this many sigma below the curve (nSigma > -nSigmaBelow)"};
  Configurable<float> nSigmaAbove{"nSigmaAbove", 5.f, "[post-bb-fit] Candidates may lie at most this many sigma above the curve (nSigma < +nSigmaAbove)"};
  Configurable<int> minTpcNCls{"minTpcNCls", 60, "Candidates must have more than this number of TPC clusters"};
  Configurable<float> maxAbsEta{"maxAbsEta", 1.f, "Candidates must have |eta| below this value"};
  Configurable<int> requiredSign{"requiredSign", -1, "[post-bb-fit] Required track charge sign (0: both)"};
  Configurable<float> relResolution{"relResolution", 0.07f, "Default relative dE/dx resolution (used if the network sigma is off)"};
  Configurable<float> nClNorm{"nClNorm", 152.f, "Normalisation of the cluster scaling sqrt(nClNorm / NCl) (default sigma and network input)"};
  Configurable<int> tfLengthOrbits{"tfLengthOrbits", 128, "Number of orbits per time frame (used for the time frame estimate)"};

  // ---------------------------------------------------------------- neural network sigma
  Configurable<bool> useNetworkSigma{"useNetworkSigma", false, "Use the sigma from the TPC PID neural network (proton hypothesis) instead of the default parametrisation"};
  Configurable<bool> autofetchNetworks{"autofetchNetworks", true, "Fetch the network from the CCDB for the timestamp of the data frame"};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB"};
  Configurable<std::string> networkPathCCDB{"networkPathCCDB", "Analysis/PID/TPC/ML", "Path of the network on the CCDB"};
  Configurable<std::string> networkPathLocally{"networkPathLocally", "network.onnx", "Local .onnx file (also the download location when autofetching)"};
  Configurable<std::string> networkVersionLocal{"networkVersionLocal", "5", "NN-Version of the local network file (ignored when autofetching)"};
  Configurable<std::string> recoPass{"recoPass", "", "Reconstruction pass name for the CCDB network query (empty: latest)"};
  Configurable<bool> enableNetworkOptimizations{"enableNetworkOptimizations", true, "Enable ORT_ENABLE_EXTENDED in the ONNX session"};
  Configurable<int> networkSetNumThreads{"networkSetNumThreads", 0, "Number of threads for the ONNX session (0: default)"};
  Configurable<std::string> cfgPathGrpLhcIf{"ccdb-path-grplhcif", "GLO/Config/GRPLHCIF", "Path of the GRPLHCIF object on the CCDB"};
  Configurable<float> phiEntranceCoeff1{"phiEntranceCoeff1", 1.026f, "phiEntrance = phi + c1 * LightSpeedDm2PS * 0.5 * c2 * 1/pT (network version 5)"};
  Configurable<float> phiEntranceCoeff2{"phiEntranceCoeff2", 85.f, "phiEntrance = phi + c1 * LightSpeedDm2PS * 0.5 * c2 * 1/pT (network version 5)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TTree> tree{"candidates", OutputObjHandlingPolicy::AnalysisObject};

  using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra>;

  // network state
  ml::OnnxModel network;
  std::map<std::string, std::string> metadata;
  std::map<std::string, std::string> headers;
  std::string networkVersion;
  std::string irSource;
  o2::common::core::CollisionSystemType::collType collsys = o2::common::core::CollisionSystemType::kCollSysUndef;
  ctpRateFetcher rateFetcher;
  bool networkReady = false;

  // selection settings resolved from the mode in init()
  std::array<float, 5> selParams{};
  int selSign = -1;
  float selBelow = 5.f;
  float selAbove = 5.f;
  float selMip = 50.f;
  float selChargeFactor = 2.3f;
  float selNClNorm = 152.f;
  bool preFitMode = false;

  // tree branches
  struct Output {
    // identification / CTF markers
    int64_t trackId = -1, collisionId = -1;
    uint64_t dataframeId = 0, globalBC = 0, orbit = 0, timestamp = 0, timeFrameEstimate = 0;
    int runNumber = 0, bcInOrbit = 0, trackType = 0, tpcSector = -1;
    float collisionTime = 0, collisionTimeRes = 0, collisionPosX = 0, collisionPosY = 0, collisionPosZ = 0;
    int collisionNContrib = 0;
    float trackTime = 0, trackTimeRes = 0, alpha = 0, x = 0, y = 0, z = 0;
    // kinematics
    float p = 0, pt = 0, tpcInnerParam = 0, signedP = 0, eta = 0, phi = 0, tgl = 0, signed1Pt = 0;
    int sign = 0;
    // TPC
    float tpcSignal = 0, tpcExpSignal = 0, tpcSigma = 0, tpcNSigma = 0, tpcRelSigma = 0, tpcChi2NCl = 0;
    int tpcNClsFound = 0, tpcNClsPID = -1, tpcNClsCrossedRows = 0, tpcNClsFindable = 0, tpcNClsShared = 0;
    // other detectors / flags
    int hasITS = 0, hasTRD = 0, hasTOF = 0, itsNCls = 0;
    float itsChi2NCl = 0, trdChi2 = 0, tofChi2 = 0;
    uint32_t flags = 0;
  } out;

  std::vector<double> makeSignedLogEdges()
  {
    std::vector<double> edges;
    const int n = nBinsPLogPerSign;
    const double lo = std::log10(static_cast<double>(pMinLog));
    const double hi = std::log10(static_cast<double>(pMax));
    for (int i = n; i >= 0; --i) {
      edges.push_back(-std::pow(10., lo + (hi - lo) * i / n));
    }
    for (int i = 0; i <= n; ++i) {
      edges.push_back(std::pow(10., lo + (hi - lo) * i / n));
    }
    return edges;
  }

  void init(InitContext&)
  {
    const AxisSpec axisDedx{static_cast<int>(nBinsDedx), 0., static_cast<double>(dEdxMax), "d#it{E}/d#it{x} (arb. units)"};
    const char* const pTitle = "#it{p} / #it{Z} (GeV/#it{c})";
    const AxisSpec axisPLinear{static_cast<int>(nBinsPLinear), -static_cast<double>(pMax), static_cast<double>(pMax), pTitle};
    const AxisSpec axisPLog{makeSignedLogEdges(), pTitle};
    histos.add("dEdxVsPLinear", "TPC d#it{E}/d#it{x} vs. signed momentum (linear bins)", HistType::kTH2F, {axisPLinear, axisDedx});
    histos.add("dEdxVsPLog", "TPC d#it{E}/d#it{x} vs. signed momentum (logarithmic bins)", HistType::kTH2F, {axisPLog, axisDedx});

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (!doSearch) {
      tree.setObject(new TTree("candidates", "candidates (disabled)"));
      return;
    }
    if (mode.value != "pre-bb-fit" && mode.value != "post-bb-fit") {
      LOG(fatal) << "mode must be 'pre-bb-fit' or 'post-bb-fit' (got '" << mode.value << "')";
    }
    const bool preFit = mode.value == "pre-bb-fit";
    const auto& par = preFit ? bbParamsStd.value : bbParams.value;
    if (par.size() != 5) {
      LOG(fatal) << "The Bethe-Bloch parameter list needs exactly 5 entries (got " << par.size() << ")";
    }
    std::copy(par.begin(), par.end(), selParams.begin());
    selSign = preFit ? preFitSign.value : requiredSign.value;
    selBelow = preFit ? preFitNSigmaBelow.value : nSigmaBelow.value;
    selAbove = preFit ? preFitNSigmaAbove.value : nSigmaAbove.value;
    selMip = bbMip.value;
    selChargeFactor = bbChargeFactor.value;
    selNClNorm = nClNorm.value;
    preFitMode = preFit;
    LOG(info) << "Mode " << mode.value << ": sign " << selSign << ", window [-" << selBelow << ", +" << selAbove << "] sigma, BB parameters " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4];

    auto* t = new TTree("candidates", preFit ? "Li6 candidates around the standard Bethe-Bloch curve (pre-bb-fit)" : "anti-Li6 candidates around the fitted Bethe-Bloch curve (post-bb-fit)");
    tree.setObject(t);
    t->Branch("trackId", &out.trackId);
    t->Branch("collisionId", &out.collisionId);
    t->Branch("dataframeId", &out.dataframeId);
    t->Branch("runNumber", &out.runNumber);
    t->Branch("globalBC", &out.globalBC);
    t->Branch("orbit", &out.orbit);
    t->Branch("bcInOrbit", &out.bcInOrbit);
    t->Branch("timestamp", &out.timestamp);
    t->Branch("timeFrameEstimate", &out.timeFrameEstimate);
    t->Branch("collisionTime", &out.collisionTime);
    t->Branch("collisionTimeRes", &out.collisionTimeRes);
    t->Branch("collisionPosX", &out.collisionPosX);
    t->Branch("collisionPosY", &out.collisionPosY);
    t->Branch("collisionPosZ", &out.collisionPosZ);
    t->Branch("collisionNContrib", &out.collisionNContrib);
    t->Branch("trackType", &out.trackType);
    t->Branch("trackTime", &out.trackTime);
    t->Branch("trackTimeRes", &out.trackTimeRes);
    t->Branch("tpcSector", &out.tpcSector);
    t->Branch("alpha", &out.alpha);
    t->Branch("x", &out.x);
    t->Branch("y", &out.y);
    t->Branch("z", &out.z);
    t->Branch("p", &out.p);
    t->Branch("pt", &out.pt);
    t->Branch("tpcInnerParam", &out.tpcInnerParam);
    t->Branch("signedP", &out.signedP);
    t->Branch("eta", &out.eta);
    t->Branch("phi", &out.phi);
    t->Branch("tgl", &out.tgl);
    t->Branch("signed1Pt", &out.signed1Pt);
    t->Branch("sign", &out.sign);
    t->Branch("tpcSignal", &out.tpcSignal);
    t->Branch("tpcExpSignal", &out.tpcExpSignal);
    t->Branch("tpcSigma", &out.tpcSigma);
    t->Branch("tpcNSigma", &out.tpcNSigma);
    t->Branch("tpcRelSigma", &out.tpcRelSigma);
    t->Branch("tpcChi2NCl", &out.tpcChi2NCl);
    t->Branch("tpcNClsFound", &out.tpcNClsFound);
    t->Branch("tpcNClsPID", &out.tpcNClsPID);
    t->Branch("tpcNClsCrossedRows", &out.tpcNClsCrossedRows);
    t->Branch("tpcNClsFindable", &out.tpcNClsFindable);
    t->Branch("tpcNClsShared", &out.tpcNClsShared);
    t->Branch("hasITS", &out.hasITS);
    t->Branch("hasTRD", &out.hasTRD);
    t->Branch("hasTOF", &out.hasTOF);
    t->Branch("itsNCls", &out.itsNCls);
    t->Branch("itsChi2NCl", &out.itsChi2NCl);
    t->Branch("trdChi2", &out.trdChi2);
    t->Branch("tofChi2", &out.tofChi2);
    t->Branch("flags", &out.flags);

    if (useNetworkSigma && !autofetchNetworks) {
      if (networkPathLocally.value.empty()) {
        LOG(fatal) << "networkPathLocally must be set when the network is not autofetched";
      }
      networkVersion = networkVersionLocal.value;
      network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value);
      std::vector<float> dummy(network.getNumInputNodes(), 1.f);
      network.evalModel(dummy); // initialisation, reduces the overhead of the first call
      networkReady = true;
    }
  }

  /// Fetch the network (and the collision system) for the timestamp of the data frame
  void prepareNetwork(uint64_t timestamp)
  {
    if (!autofetchNetworks) {
      return;
    }
    if (irSource.empty()) {
      if (auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>(cfgPathGrpLhcIf.value, timestamp)) {
        collsys = o2::common::core::CollisionSystemType::getCollisionTypeFromGrp(grpo);
        irSource = collsys == o2::common::core::CollisionSystemType::kCollSyspp ? "T0VTX" : "ZNC hadronic";
      } else {
        LOG(warning) << "No GRPLHCIF object found, hadronic rate will be 0";
      }
    }
    if (networkReady && timestamp >= network.getValidityFrom() && timestamp <= network.getValidityUntil()) {
      return;
    }
    metadata.clear();
    if (!recoPass.value.empty()) {
      metadata["RecoPassName"] = recoPass.value;
    }
    LOG(info) << "Fetching network for timestamp " << timestamp;
    const bool ok = ccdb->getCCDBAccessor().retrieveBlob(networkPathCCDB.value, ".", metadata, timestamp, false, networkPathLocally.value, "", "", &headers);
    if (!ok) {
      LOG(fatal) << "No valid TPC PID network found on the CCDB for timestamp " << timestamp;
    }
    networkVersion = headers["NN-Version"];
    network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), nullptr, 0), strtoul(headers["Valid-Until"].c_str(), nullptr, 0));
    std::vector<float> dummy(network.getNumInputNodes(), 1.f);
    network.evalModel(dummy);
    networkReady = true;
  }

  /// pre-bb-fit: take the standard curve from the TPC response object on the CCDB
  void updateStdResponse(uint64_t timestamp)
  {
    if (!preFitMode || !useCcdbStdParams) {
      return;
    }
    std::map<std::string, std::string> md;
    std::map<std::string, std::string> hd;
    if (!recoPass.value.empty()) {
      md["RecoPassName"] = recoPass.value;
    }
    auto* response = ccdb->getSpecific<o2::pid::tpc::Response>(ccdbPathResponse.value, timestamp, md, &hd);
    if (!response) {
      LOG(warning) << "No TPC response for pass '" << recoPass.value << "', falling back to the latest object";
      response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPathResponse.value, timestamp, &hd);
    }
    if (!response) {
      LOG(fatal) << "No TPC response object found on the CCDB at " << ccdbPathResponse.value << " for timestamp " << timestamp;
    }
    const auto par = response->GetBetheBlochParams();
    if (par != selParams || response->GetMIP() != selMip || response->GetChargeFactor() != selChargeFactor || response->GetNClNormalization() != selNClNorm) {
      selParams = par;
      selMip = response->GetMIP();
      selChargeFactor = response->GetChargeFactor();
      selNClNorm = response->GetNClNormalization();
      LOG(info) << "Standard curve from CCDB (period " << hd["LPMProductionTag"] << ", pass " << hd["RecoPassName"] << "): BB " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << ", MIP " << selMip << ", chargeFactor " << selChargeFactor << ", nClNorm " << selNClNorm;
    }
  }

  float expectedSignal(float tpcInnerParam)
  {
    const float bg = (rigidityToMomentum ? chargeLi6 * tpcInnerParam : tpcInnerParam) / massLi6;
    const auto& par = selParams;
    return selMip * o2::common::BetheBlochAleph(bg, par[0], par[1], par[2], par[3], par[4]) * std::pow(chargeLi6.value, selChargeFactor);
  }

  template <typename TOrigins>
  void run(Colls const& collisions, TrackCandidates const& tracks, BCsWithTimestamps const& bcs, [[maybe_unused]] TOrigins const* origins)
  {
    constexpr float OneToKilo = 1.e-3f;
    constexpr float MultiplicityNorm = 11000.f;
    constexpr float Ft0cOccupancyNorm = 60000.f;
    constexpr float HadronicRateNormPp = 1500.f;
    constexpr float HadronicRateNormAa = 50.f;
    constexpr float LightSpeedDm2PS = o2::constants::physics::LightSpeedCm2PS / 10.f;

    // ---------------------------------------------------------------- pass 1: histograms and candidate preselection
    std::vector<int32_t> candidateSlot;
    std::vector<size_t> candidateRow;
    std::vector<float> nnInput;
    std::vector<int> nTpcTracks;
    int nnVersion = 0;
    int nFeatures = 0;
    float hadronicRateNorm = HadronicRateNormPp;
    std::vector<float> hadronicRate;

    if (doSearch) {
      if (bcs.size() > 0) {
        updateStdResponse(bcs.begin().timestamp());
      }
      candidateSlot.assign(tracks.size(), -1);
      if (useNetworkSigma) {
        if (bcs.size() == 0) {
          return;
        }
        prepareNetwork(bcs.begin().timestamp());
        nFeatures = network.getNumInputNodes();
        struct VersionEntry {
          std::string_view name;
          int nFeatures;
          int version;
        };
        constexpr std::array<VersionEntry, 6> Versions{{{"", 6, 1}, {"1", 6, 1}, {"2", 7, 2}, {"3", 8, 3}, {"4", 9, 4}, {"5", 9, 5}}};
        for (const auto& v : Versions) {
          if (networkVersion == v.name && nFeatures == v.nFeatures) {
            nnVersion = v.version;
            break;
          }
        }
        if (nnVersion == 0) {
          LOG(fatal) << "Network version '" << networkVersion << "' with " << nFeatures << " input features is not supported";
        }
        hadronicRateNorm = collsys == o2::common::core::CollisionSystemType::kCollSyspp ? HadronicRateNormPp : HadronicRateNormAa;
        hadronicRate.assign(collisions.size(), -1.f);
        nTpcTracks.assign(collisions.size(), 0);
        for (const auto& trk : tracks) {
          if (trk.hasTPC() && trk.has_collision()) {
            ++nTpcTracks[trk.collisionId()];
          }
        }
      }
    }

    size_t row = 0;
    for (const auto& trk : tracks) {
      const size_t thisRow = row++;
      if (!trk.hasTPC() || trk.tpcSignal() <= 0.f) {
        continue;
      }
      const float signedP = trk.sign() * trk.tpcInnerParam();
      histos.fill(HIST("dEdxVsPLinear"), signedP, trk.tpcSignal());
      histos.fill(HIST("dEdxVsPLog"), signedP, trk.tpcSignal());

      if (!doSearch || (selSign != 0 && trk.sign() != selSign) || trk.tpcNClsFound() <= minTpcNCls || std::abs(trk.eta()) >= maxAbsEta) {
        continue;
      }
      candidateSlot[thisRow] = static_cast<int32_t>(candidateRow.size());
      candidateRow.push_back(thisRow);

      if (useNetworkSigma) {
        // Same features as createNetworkPrediction() in pidTPCModule.h, with the proton mass hypothesis
        const bool hasColl = trk.has_collision() && collisions.size() > 0;
        const auto base = nnInput.size();
        nnInput.resize(base + nFeatures, 1.f);
        nnInput[base + 0] = trk.tpcInnerParam();
        nnInput[base + 1] = trk.tgl();
        nnInput[base + 2] = trk.signed1Pt();
        nnInput[base + 3] = o2::track::pid_constants::sMasses[o2::track::PID::Proton];
        nnInput[base + 4] = hasColl ? nTpcTracks[trk.collisionId()] / MultiplicityNorm : 1.f;
        nnInput[base + 5] = std::sqrt(selNClNorm / trk.tpcNClsFound());
        if (nnVersion >= 2) { // FT0C occupancy
          nnInput[base + 6] = hasColl ? collisions.rawIteratorAt(trk.collisionId()).ft0cOccupancyInTimeRange() / Ft0cOccupancyNorm : 1.f;
        }
        if (nnVersion >= 3) { // hadronic rate
          float rate = 0.f;
          if (hasColl && !irSource.empty()) {
            float& cached = hadronicRate[trk.collisionId()];
            if (cached < 0.f) {
              const auto bc = collisions.rawIteratorAt(trk.collisionId()).template bc_as<BCsWithTimestamps>();
              cached = rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * OneToKilo;
            }
            rate = cached;
          }
          nnInput[base + 7] = rate / hadronicRateNorm;
        }
        if (nnVersion >= 4) { // phi modulo sector
          float phi = trk.phi();
          if (nnVersion == 5) {
            phi += phiEntranceCoeff1 * LightSpeedDm2PS * 0.5f * phiEntranceCoeff2 * trk.signed1Pt();
          }
          nnInput[base + 8] = RecoDecay::constrainAngle(phi, 0.f, o2::constants::math::NSectors);
        }
      }
    }

    if (!doSearch || candidateRow.empty()) {
      return;
    }

    // ---------------------------------------------------------------- network evaluation (one batch)
    // Output layout per track: 1 node: mean only (no sigma -> default sigma), 2 nodes: {mean, mean + sigma},
    // 3 nodes: {mean, mean + sigma_up, mean - sigma_down}. Only the width is used, the proton mean correction is not applied.
    std::vector<float> nnOutput;
    int nOutputs = 0;
    if (useNetworkSigma) {
      network.evalModel(nnInput, nnOutput);
      nOutputs = network.getNumOutputNodes();
      if (nnOutput.size() != static_cast<size_t>(nOutputs) * candidateRow.size()) {
        LOG(fatal) << "Network output size " << nnOutput.size() << " does not match " << nOutputs << " x " << candidateRow.size();
      }
    }

    // ---------------------------------------------------------------- pass 2: selection and output
    uint64_t dataframeId = 0; // 0: not available (use processWithOrigins)
    if constexpr (!std::is_same_v<TOrigins, std::nullptr_t>) {
      if (origins != nullptr && origins->size() > 0) {
        dataframeId = origins->begin().dataframeID();
      }
    }
    row = 0;
    for (const auto& trk : tracks) {
      const int32_t slot = candidateSlot[row++];
      if (slot < 0) {
        continue;
      }
      const float expSignal = expectedSignal(trk.tpcInnerParam());
      if (!(expSignal > 0.f)) {
        continue;
      }
      float relSigma = relResolution * std::sqrt(selNClNorm / trk.tpcNClsFound());
      if (useNetworkSigma && nOutputs >= 2) {
        const float* o = nnOutput.data() + static_cast<size_t>(nOutputs) * slot;
        const float sigmaUp = o[1] - o[0];
        const float sigmaDown = nOutputs == 3 ? o[0] - o[2] : sigmaUp;
        relSigma = trk.tpcSignal() >= expSignal ? sigmaUp : sigmaDown;
      }
      const float sigma = relSigma * expSignal;
      if (!(sigma > 0.f)) {
        continue;
      }
      const float nSigma = (trk.tpcSignal() - expSignal) / sigma;
      if (nSigma <= -selBelow || nSigma >= selAbove) {
        continue;
      }

      out = Output{};
      out.trackId = trk.globalIndex();
      out.dataframeId = dataframeId;
      out.trackType = trk.trackType();
      out.trackTime = trk.trackTime();
      out.trackTimeRes = trk.trackTimeRes();
      out.alpha = trk.alpha();
      out.x = trk.x();
      out.y = trk.y();
      out.z = trk.z();
      const float alpha = out.alpha < 0.f ? out.alpha + o2::constants::math::TwoPI : out.alpha;
      out.tpcSector = static_cast<int>(alpha / (o2::constants::math::TwoPI / o2::constants::math::NSectors)) % o2::constants::math::NSectors;
      if (trk.has_collision()) {
        const auto coll = collisions.rawIteratorAt(trk.collisionId());
        out.collisionId = coll.globalIndex();
        out.collisionTime = coll.collisionTime();
        out.collisionTimeRes = coll.collisionTimeRes();
        out.collisionPosX = coll.posX();
        out.collisionPosY = coll.posY();
        out.collisionPosZ = coll.posZ();
        out.collisionNContrib = coll.numContrib();
        const auto bc = coll.template bc_as<BCsWithTimestamps>();
        out.runNumber = bc.runNumber();
        out.globalBC = bc.globalBC();
        out.timestamp = bc.timestamp();
        out.orbit = out.globalBC / o2::constants::lhc::LHCMaxBunches;
        out.bcInOrbit = static_cast<int>(out.globalBC % o2::constants::lhc::LHCMaxBunches);
        out.timeFrameEstimate = out.orbit / tfLengthOrbits;
      }
      out.p = trk.p();
      out.pt = trk.pt();
      out.tpcInnerParam = trk.tpcInnerParam();
      out.signedP = trk.sign() * trk.tpcInnerParam();
      out.eta = trk.eta();
      out.phi = trk.phi();
      out.tgl = trk.tgl();
      out.signed1Pt = trk.signed1Pt();
      out.sign = trk.sign();
      out.tpcSignal = trk.tpcSignal();
      out.tpcExpSignal = expSignal;
      out.tpcSigma = sigma;
      out.tpcNSigma = nSigma;
      out.tpcRelSigma = relSigma;
      out.tpcChi2NCl = trk.tpcChi2NCl();
      out.tpcNClsFound = trk.tpcNClsFound();
      if constexpr (requires { trk.tpcNClsPID(); }) {
        out.tpcNClsPID = trk.tpcNClsPID();
      }
      out.tpcNClsCrossedRows = trk.tpcNClsCrossedRows();
      out.tpcNClsFindable = trk.tpcNClsFindable();
      out.tpcNClsShared = trk.tpcNClsShared();
      out.hasITS = trk.hasITS();
      out.hasTRD = trk.hasTRD();
      out.hasTOF = trk.hasTOF();
      out.itsNCls = trk.itsNCls();
      out.itsChi2NCl = trk.itsChi2NCl();
      out.trdChi2 = trk.trdChi2();
      out.tofChi2 = trk.tofChi2();
      out.flags = trk.flags();
      tree->Fill();
    }
  }

  void process(Colls const& collisions, TrackCandidates const& tracks, BCsWithTimestamps const& bcs)
  {
    run<std::nullptr_t>(collisions, tracks, bcs, nullptr);
  }

  /// Also stores the data frame ID (needs the ORIGIN table in the AO2D)
  void processWithOrigins(Colls const& collisions, TrackCandidates const& tracks, BCsWithTimestamps const& bcs, aod::Origins const& origins)
  {
    run(collisions, tracks, bcs, &origins);
  }
  PROCESS_SWITCH(TpcDedxAntiLi6Search, processWithOrigins, "Also store the data frame ID (requires the ORIGIN table)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TpcDedxAntiLi6Search>(cfgc)};
}
