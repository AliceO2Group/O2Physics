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

/// \file phosPi0.cxx
/// \brief PHOS pi0/eta analysis
/// \author Dmitri Peresunko <Dmitri.Peresunko@cern.ch>
///

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhosPi0 {
  Configurable<bool> skimmedProcessing{"skimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<std::string> trigName{"trigName", "fPHOSPhoton", "name of offline trigger"};
  Configurable<std::string> zorroCCDBpath{"zorroCCDBpath", "/Users/m/mpuccio/EventFiltering/OTS/", "path to the zorro ccdb objects"};
  Configurable<int> evSelTrig{"evSelTrig", aod::evsel::kIsTriggerTVX, "Select events with this trigger"};
  Configurable<bool> isMC{"isMC", false, "to fill MC histograms"};
  Configurable<float> minCluE{"minCluE", 0.3, "Minimum cluster energy for analysis"};
  Configurable<float> minCluTime{"minCluTime", -25.e-9, "Min. cluster time"};
  Configurable<float> maxCluTime{"maxCluTime", 25.e-9, "Max. cluster time"};
  Configurable<int> minCluNcell{"minCluNcell", 2, "min cells in cluster"};
  Configurable<float> minM02{"minM02", 0.2, "Min disp M02 cut"};
  Configurable<float> cpvCut{"cpvCut", 2., "Min distance to track"};
  Configurable<int> nMixedEvents{"nMixedEvents", 10, "number of events to mix"};
  Configurable<bool> fillQC{"fillQC", true, "Fill QC histos"};
  Configurable<float> minOccE{"minOccE", 0.5, "Min. cluster energy of occupancy plots"};
  Configurable<float> nonlinA{"nonlinA", 1., "nonlinsrity param A (scale)"};
  Configurable<float> nonlinB{"nonlinB", 0., "nonlinsrity param B (a+b*exp(-e/c))"};
  Configurable<float> nonlinC{"nonlinC", 1., "nonlinsrity param C (a+b*exp(-e/c))"};
  Configurable<int> tofEffParam{"tofEffParam", 0, "parameterization of TOF cut efficiency"};
  Configurable<float> timeOffset{"timeOffset", 0., "time offset to compensate imperfection of time calibration"};

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using SelCollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using BCsWithBcSels = soa::Join<aod::BCs, aod::BcSels>;
  using McClusters = soa::Join<aod::CaloClusters, aod::PHOSCluLabels>;
  using McAmbClusters = soa::Join<aod::CaloAmbiguousClusters, aod::PHOSAmbCluLabels>;

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  HistogramRegistry mHistManager{"PHOSPi0Histograms"};

  // class to keep photon candidate parameters
  class Photon
  {
   public:
    Photon() = default;
    Photon(double x, double y, double z, double ee, double t, int m, bool isDispOK, bool isCPVOK, int mcLabel) : px(x), py(y), pz(z), e(ee), time(t), mod(m), mPID(isDispOK << 1 | isCPVOK << 2), label(mcLabel) {}
    ~Photon() = default;

    bool isCPVOK() const { return (mPID >> 2) & 1; }
    bool isDispOK() const { return (mPID >> 1) & 1; }
    double pt() const { return std::sqrt(px * px + py * py); }

   public:
    double px = 0.;   // px
    double py = 0.;   // py
    double pz = 0.;   // pz
    double e = 0.;    // energy
    double time = 0.; // time
    int mod = 0;      // module
    int mPID = 0;     // store PID bits
    int label = -1;   // label of MC particle
  };

  int mRunNumber = 0;    // Current run number
  int mixedEventBin = 0; // Which list of Mixed use for mixing
  std::vector<Photon> mCurEvent;
  static constexpr int kMaxMixBins = 20; // maximal number of kinds of events for mixing
  std::array<std::deque<std::vector<Photon>>, kMaxMixBins> mMixedEvents;
  std::array<std::deque<std::vector<Photon>>, kMaxMixBins> mAmbMixedEvents;

  int mPrevMCColId = -1; // mark MC collissions already scanned
  // fast access to histos
  TH1* hColl;
  TH3 *hReMod, *hMiMod;
  TH2 *hReAll, *hReDisp, *hReCPV, *hReBoth, *hSignalAll, *hPi0SignalAll, *hPi0SignalCPV, *hPi0SignalDisp,
    *hPi0SignalBoth, *hMiAll, *hMiDisp, *hMiCPV, *hMiBoth;
  TH2 *hReOneAll, *hReOneDisp, *hReOneCPV, *hReOneBoth, *hMiOneAll, *hMiOneDisp, *hMiOneCPV, *hMiOneBoth;
  TH2 *hReTime12, *hReTime30, *hReTime50, *hReTime100;

  std::vector<double> pt = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1,
                            1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.5, 4.6, 4.8, 5.0,
                            5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28.,
                            30., 34., 38., 42., 46., 50., 55., 60., 70., 75., 80., 85., 90., 95., 100., 110., 120., 130., 140., 150., 160., 180., 200.};

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS pi0 analysis task ...";
    zorroSummary.setObject(zorro.getZorroSummary());
    zorro.setBaseCCDBPath(zorroCCDBpath.value);

    const AxisSpec ptAxis{pt, "p_{T} (GeV/c)"},
      mggAxis{625, 0., 1.25, "m_{#gamma#gamma} (GeV/c^{2})"},
      timeAxis{100, -500.e-9, 500.e-9, "t (s)"},
      m02Axis{100, 0., 20., "M02 (cm^{2})"},
      m20Axis{100, 0., 20., "M20 (cm^{2})"},
      nCellsAxis{100, 0., 100., "N_{cell}"},
      zAxis{56, -63., 63., "z", "z (cm)"},
      phiAxis{64, -72., 72., "x", "x (cm)"},
      modAxis{4, 1., 5., "module", "Module"},
      multAxis{100, 0., 100.},
      vertexAxis{100, -20., 20., "z", "z (cm)"},
      modCombAxis{10, 0., 10.},
      centAxis{10, 0., 10.},
      centralityAxis{100, 0., 100., "centrality", "centrality"};

    hColl = std::get<std::shared_ptr<TH1>>(mHistManager.add("eventsCol", "Number of events", HistType::kTH1F, {{10, 0., 10.}})).get();
    hColl->GetXaxis()->SetBinLabel(1, "All");
    hColl->GetXaxis()->SetBinLabel(2, "SwTr");
    hColl->GetXaxis()->SetBinLabel(3, "T0a||T0c");
    hColl->GetXaxis()->SetBinLabel(4, "T0a&&T0c");
    hColl->GetXaxis()->SetBinLabel(5, "kTVXinPHOS");
    hColl->GetXaxis()->SetBinLabel(6, "kIsTriggerTVX");
    hColl->GetXaxis()->SetBinLabel(7, "PHOSClu");
    hColl->GetXaxis()->SetBinLabel(8, "PHOSClu&&kTVXinPHOS");
    hColl->GetXaxis()->SetBinLabel(9, "Accepted");

    auto h2{std::get<std::shared_ptr<TH1>>(mHistManager.add("eventsBC", "Number of events per trigger", HistType::kTH1F, {{8, 0., 8.}}))};
    h2->GetXaxis()->SetBinLabel(1, "All");
    h2->GetXaxis()->SetBinLabel(2, "T0a||T0c");
    h2->GetXaxis()->SetBinLabel(3, "T0a&&T0c");
    h2->GetXaxis()->SetBinLabel(4, "alias kTVXinPHOS");
    h2->GetXaxis()->SetBinLabel(5, "kIsTriggerTVX");
    h2->GetXaxis()->SetBinLabel(6, "kTVXinPHOS");
    h2->GetXaxis()->SetBinLabel(7, "PHOSClu");
    h2->GetXaxis()->SetBinLabel(8, "PHOSClu&&Trig");

    mHistManager.add("contributors", "Contributors per collision", HistType::kTH2F, {{10, 0., 10.}, {10, 0., 100.}});
    mHistManager.add("vertex", "vertex", HistType::kTH1F, {vertexAxis});

    if (fillQC) {
      // QC histograms for normal collisions
      mHistManager.add("cluSp", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("cluSpDisp", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("cluSpCPV", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("cluSpBoth", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("hM02Clu", "(M02,M20) in clusters", HistType::kTH2F, {ptAxis, m02Axis});
      mHistManager.add("hM20Clu", "(M02,M20) in clusters", HistType::kTH2F, {ptAxis, m20Axis});
      mHistManager.add("hNcellClu", "Number of cells in clusters", HistType::kTH3F, {ptAxis, nCellsAxis, modAxis});
      mHistManager.add("cluETime", "Cluster time vs E",
                       HistType::kTH3F, {ptAxis, timeAxis, modAxis});
      mHistManager.add("cluOcc", "Cluster occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
      mHistManager.add("cluCPVOcc", "Cluster with CPV occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
      mHistManager.add("cluDispOcc", "Cluster with Disp occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
      mHistManager.add("cluBothOcc", "Cluster with Both occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    }

    hReMod = std::get<std::shared_ptr<TH3>>(mHistManager.add("mggReModComb", "inv mass per module",
                                                             HistType::kTH3F, {mggAxis, ptAxis, modCombAxis}))
               .get();
    hReAll = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReAll", "inv mass for centrality",
                                                             HistType::kTH2F, {mggAxis, ptAxis}))
               .get();
    hReCPV = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReCPV", "inv mass for centrality",
                                                             HistType::kTH2F, {mggAxis, ptAxis}))
               .get();
    hReDisp = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReDisp", "inv mass for centrality",
                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                .get();
    hReBoth = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReBoth", "inv mass for centrality",
                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                .get();
    hReOneAll = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReOneAll", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hReOneCPV = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReOneCPV", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hReOneDisp = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReOneDisp", "inv mass for centrality",
                                                                 HistType::kTH2F, {mggAxis, ptAxis}))
                   .get();
    hReOneBoth = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReOneBoth", "inv mass for centrality",
                                                                 HistType::kTH2F, {mggAxis, ptAxis}))
                   .get();

    hReTime12 = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReTime12", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hReTime30 = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReTime30", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hReTime50 = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReTime50", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hReTime100 = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggReTime100", "inv mass for centrality",
                                                                 HistType::kTH2F, {mggAxis, ptAxis}))
                   .get();

    if (isMC) {
      hSignalAll = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggSignal", "inv mass for correlated pairs",
                                                                   HistType::kTH2F, {mggAxis, ptAxis}))
                     .get();
      hPi0SignalAll = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggPi0SignalAll", "inv mass for pi0 pairs",
                                                                      HistType::kTH2F, {mggAxis, ptAxis}))
                        .get();
      hPi0SignalCPV = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggPi0SignalCPV", "inv mass for pi0 pairs",
                                                                      HistType::kTH2F, {mggAxis, ptAxis}))
                        .get();
      hPi0SignalDisp = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggPi0SignalDisp", "inv mass for pi0 pairs",
                                                                       HistType::kTH2F, {mggAxis, ptAxis}))
                         .get();
      hPi0SignalBoth = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggPi0SignalBoth", "inv mass for pi0 pairs",
                                                                       HistType::kTH2F, {mggAxis, ptAxis}))
                         .get();
    }

    hMiMod = std::get<std::shared_ptr<TH3>>(mHistManager.add("mggMiModComb", "inv mass for centrality",
                                                             HistType::kTH3F, {mggAxis, ptAxis, modCombAxis}))
               .get();
    hMiAll = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiAll", "inv mass for centrality",
                                                             HistType::kTH2F, {mggAxis, ptAxis}))
               .get();
    hMiCPV = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiCPV", "inv mass for centrality",
                                                             HistType::kTH2F, {mggAxis, ptAxis}))
               .get();
    hMiDisp = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiDisp", "inv mass for centrality",
                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                .get();
    hMiBoth = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiBoth", "inv mass for centrality",
                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                .get();
    hMiOneAll = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiOneAll", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hMiOneCPV = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiOneCPV", "inv mass for centrality",
                                                                HistType::kTH2F, {mggAxis, ptAxis}))
                  .get();
    hMiOneDisp = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiOneDisp", "inv mass for centrality",
                                                                 HistType::kTH2F, {mggAxis, ptAxis}))
                   .get();
    hMiOneBoth = std::get<std::shared_ptr<TH2>>(mHistManager.add("mggMiOneBoth", "inv mass for centrality",
                                                                 HistType::kTH2F, {mggAxis, ptAxis}))
                   .get();
    if (isMC) {
      mHistManager.add("hMCPi0SpAll", "pi0 spectrum inclusive", HistType::kTH1F, {ptAxis});
      mHistManager.add("hMCPi0SpPrim", "pi0 spectrum Primary", HistType::kTH1F, {ptAxis});
      mHistManager.add("hMCPi0RapPrim", "pi0 rapidity primary", HistType::kTH1F, {{100, -1., 1., "Rapidity"}});
      mHistManager.add("hMCPi0PhiPrim", "pi0 phi primary", HistType::kTH1F, {{100, 0., o2::constants::math::TwoPI, "#phi (rad)"}});
      mHistManager.add("hMCPi0SecVtx", "pi0 secondary", HistType::kTH2F, {{100, 0., 500., "R (cm)"}, {100, -o2::constants::math::PI, o2::constants::math::PI, "#phi (rad)"}});
    }
  }

  // template <typename Tcoll>
  // float getCentrality(Tcoll const& collision)
  // {
  //   float centrality = 1.;
  //   if constexpr (std::is_same<Tcoll, CollWithCent>::value || std::is_same<Tcoll, CollWithEP>::value || std::is_same<Tcoll, CollWithQvec>::value) {
  //     if (cfgCentralityEstimator == nuclei::centDetectors::kFV0A) {
  //       centrality = collision.centFV0A();
  //     } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0M) {
  //       centrality = collision.centFT0M();
  //     } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0A) {
  //       centrality = collision.centFT0A();
  //     } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0C) {
  //       centrality = collision.centFT0C();
  //     } else {
  //       LOG(warning) << "Centrality estimator not valid. Possible values: (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3). Centrality set to 1.";
  //     }
  //   }
  //   return centrality;
  // }

  /// \brief Process PHOS data
  void processData(SelCollisions::iterator const& col,
                   aod::CaloClusters const& clusters,
                   aod::BCsWithTimestamps const&)
  {
    aod::McParticles const* mcPart = nullptr;
    scanAll<false>(col, clusters, mcPart);
  }
  PROCESS_SWITCH(PhosPi0, processData, "processData", true);
  void processMC(SelCollisionsMC::iterator const& col,
                 McClusters const& clusters,
                 aod::McParticles const& mcPart,
                 aod::McCollisions const& /*mcCol*/,
                 aod::BCsWithTimestamps const&)
  {
    scanAll<true>(col, clusters, &mcPart);
  }
  PROCESS_SWITCH(PhosPi0, processMC, "processMC", false);

  template <bool isMC, typename TCollision, typename TClusters>
  void scanAll(TCollision& col,
               TClusters& clusters,
               aod::McParticles const* mcPart)
  {
    mixedEventBin = 0;
    hColl->Fill(0.5);
    if (skimmedProcessing) {
      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      if (mRunNumber != bc.runNumber()) {
        zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), trigName);
        zorro.populateHistRegistry(mHistManager, bc.runNumber());
        mRunNumber = bc.runNumber();
      }

      if (!zorro.isSelected(bc.globalBC())) {
        return; ///
      }
    } else {
      if (!col.selection_bit(evSelTrig)) {
        return;
      }
    }

    hColl->Fill(1.5);
    const double vtxCut = 10.;
    double vtxZ = col.posZ();
    mHistManager.fill(HIST("vertex"), vtxZ);
    bool isColSelected = false;
    if constexpr (isMC) {
      isColSelected = (col.selection_bit(kIsTriggerTVX) && (clusters.size() > 0));
    } else {
      isColSelected = col.selection_bit(evSelTrig) && std::abs(vtxZ) < vtxCut; // col.alias_bit(evSelTrig)
      //               collision.selection_bit(aod::evsel::kNoTimeFrameBorder);
    }

    if (col.selection_bit(kIsBBT0A) || col.selection_bit(kIsBBT0C)) {
      hColl->Fill(2.5);
    }
    if (col.selection_bit(kIsBBT0A) && col.selection_bit(kIsBBT0C)) {
      hColl->Fill(3.5);
    }
    if (col.alias_bit(kTVXinPHOS)) {
      hColl->Fill(4.5);
    }
    if (col.selection_bit(kIsTriggerTVX)) {
      hColl->Fill(5.5);
    }
    if (clusters.size() > 0) {
      hColl->Fill(6.5);
      if (col.alias_bit(kTVXinPHOS)) {
        hColl->Fill(7.5);
      }
    }
    // //Event Plane| jet orientation
    //   if (flag & (kProton | kDeuteron | kTriton | kHe3 | kHe4) || doprocessMC) { /// ignore PID pre-selections for the MC
    //     if constexpr (std::is_same<Tcoll, CollWithEP>::value) {
    //       nuclei::candidates_flow.emplace_back(NucleusCandidateFlow{
    //         collision.centFV0A(),
    //         collision.centFT0M(),
    //         collision.centFT0A(),
    //         collision.centFT0C(),
    //         collision.psiFT0A(),
    //         collision.multFT0A(),
    //         collision.psiFT0C(),
    //         collision.multFT0C(),
    //         collision.psiTPC(),
    //         collision.psiTPCL(),
    //         collision.psiTPCR(),
    //         collision.multTPC()});
    //     } else if constexpr (std::is_same<Tcoll, CollWithQvec>::value) {
    //       nuclei::candidates_flow.emplace_back(NucleusCandidateFlow{
    //         collision.centFV0A(),
    //         collision.centFT0M(),
    //         collision.centFT0A(),
    //         collision.centFT0C(),
    //         0.5 * std::atan2(collision.qvecFT0AIm(), collision.qvecFT0ARe()),
    //         collision.multFT0A(),
    //         0.5 * std::atan2(collision.qvecFT0CIm(), collision.qvecFT0CRe()),
    //         collision.multFT0C(),
    //         -999.,
    //         0.5 * std::atan2(collision.qvecBNegIm(), collision.qvecBNegRe()),
    //         0.5 * std::atan2(collision.qvecBPosIm(), collision.qvecBPosRe()),
    //         collision.multTPC()});
    //     }

    int mult = 1.; // multiplicity TODO!!!
    mixedEventBin = findMixedEventBin(vtxZ, mult);

    if (!isColSelected) {
      return;
    }
    hColl->Fill(8.5);

    // Fill MC distributions
    // pion rapidity, pt, phi
    // secondary pi0s
    if constexpr (isMC) {
      // check current collision Id for clusters
      int cluMcBCId = -1;
      for (const auto& clu : clusters) {
        auto mcList = clu.labels(); // const std::vector<int>
        int nParents = mcList.size();
        for (int iParent = 0; iParent < nParents; iParent++) { // Not found nbar parent yiet
          int label = mcList[iParent];
          if (label > -1) {
            auto parent = mcPart->iteratorAt(label);
            cluMcBCId = parent.mcCollision().bcId();
            break;
          }
        }
        if (cluMcBCId > -1) {
          break;
        }
      }
      if (mcPart->begin() != mcPart->end()) {
        if (mcPart->begin().mcCollisionId() != mPrevMCColId) {
          mPrevMCColId = mcPart->begin().mcCollisionId(); // to avoid scanning full MC table each BC
          for (const auto& part : *mcPart) {
            if (part.mcCollision().bcId() != col.bcId()) {
              continue;
            }
            if (part.pdgCode() == 111) {
              double r = std::sqrt(std::pow(part.vx(), 2) + std::pow(part.vy(), 2));
              if (r < 0.5) {
                mHistManager.fill(HIST("hMCPi0RapPrim"), part.y());
              }
              if (std::abs(part.y()) < .5) {
                double pt = part.pt();
                mHistManager.fill(HIST("hMCPi0SpAll"), pt);
                double phiVtx = std::atan2(part.vy(), part.vx());
                if (r > 0.5) {
                  mHistManager.fill(HIST("hMCPi0SecVtx"), r, phiVtx);
                } else {
                  mHistManager.fill(HIST("hMCPi0SpPrim"), pt);
                  mHistManager.fill(HIST("hMCPi0PhiPrim"), part.phi());
                }
              }
            }
          }
        }
      }
    }

    // Fill invariant mass distributions
    mCurEvent.clear();
    for (const auto& clu : clusters) {
      // Fill QC histos
      if (fillQC) {
        mHistManager.fill(HIST("hM02Clu"), clu.e(), clu.m02());
        mHistManager.fill(HIST("hM20Clu"), clu.e(), clu.m20());
        mHistManager.fill(HIST("hNcellClu"), clu.e(), clu.ncell(), clu.mod());
        mHistManager.fill(HIST("cluETime"), clu.e(), clu.time(), clu.mod());
      }
      if (clu.e() < minCluE ||
          clu.ncell() < minCluNcell ||
          clu.time() > maxCluTime || clu.time() < minCluTime ||
          clu.m02() < minM02) {
        continue;
      }
      if (fillQC) {
        mHistManager.fill(HIST("cluSp"), clu.e(), clu.mod());
        if (clu.e() > minOccE) {
          mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
          if (clu.trackdist() > 2.) {
            mHistManager.fill(HIST("cluCPVOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpCPV"), clu.e(), clu.mod());
            if (testLambda(clu.e(), clu.m02(), clu.m20())) {
              mHistManager.fill(HIST("cluBothOcc"), clu.x(), clu.z(), clu.mod());
              mHistManager.fill(HIST("cluSpBoth"), clu.e(), clu.mod());
            }
          }
          if (testLambda(clu.e(), clu.m02(), clu.m20())) {
            mHistManager.fill(HIST("cluDispOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpDisp"), clu.e(), clu.mod());
          }
        }
      }

      int mcLabel = -1;
      if constexpr (isMC) {         // test parent
        auto mcList = clu.labels(); // const std::vector<int>
        if (mcList.size() > 0) {
          mcLabel = mcList[0];
        }
      }
      double enCorr = 1;
      if constexpr (isMC) { // correct MC energy
        enCorr = nonlinearity(clu.e());
      }
      Photon ph1(clu.px() * enCorr, clu.py() * enCorr, clu.pz() * enCorr, clu.e() * enCorr, clu.time(), clu.mod(), testLambda(clu.e(), clu.m02(), clu.m20()), clu.trackdist() > cpvCut, mcLabel);
      // Mix with other photons added to stack
      for (const auto& ph2 : mCurEvent) {
        double m = std::pow(ph1.e + ph2.e, 2) - std::pow(ph1.px + ph2.px, 2) -
                   std::pow(ph1.py + ph2.py, 2) - std::pow(ph1.pz + ph2.pz, 2);
        if (m > 0) {
          m = std::sqrt(m);
        }
        double pt = std::sqrt(std::pow(ph1.px + ph2.px, 2) +
                              std::pow(ph1.py + ph2.py, 2));
        int modComb = moduleCombination(ph1.mod, ph2.mod);
        double w = 1.;
        if constexpr (isMC) { // correct MC energy
          w = tofCutEff(ph1.e) * tofCutEff(ph2.e);
        }
        hReMod->Fill(m, pt, modComb, w);
        hReAll->Fill(m, pt, w);
        hReOneAll->Fill(m, ph1.pt(), w);
        hReOneAll->Fill(m, ph2.pt(), w);
        if (ph1.isCPVOK()) {
          hReOneCPV->Fill(m, ph1.pt(), w);
        }
        if (ph2.isCPVOK()) {
          hReOneCPV->Fill(m, ph2.pt(), w);
        }
        if (ph1.isDispOK()) {
          hReOneDisp->Fill(m, ph1.pt(), w);
          if (ph1.isCPVOK()) {
            hReOneBoth->Fill(m, ph1.pt(), w);
          }
        }
        if (ph2.isDispOK()) {
          hReOneDisp->Fill(m, ph2.pt(), w);
          if (ph2.isCPVOK()) {
            hReOneBoth->Fill(m, ph2.pt(), w);
          }
        }
        // Test time eff
        if (std::abs(ph1.time - timeOffset) < 12.5e-9) { // strict cut on first photon
          if (std::abs(ph2.time - timeOffset) < 100.e-9) {
            hReTime100->Fill(m, ph2.pt());
            if (std::abs(ph2.time - timeOffset) < 50.e-9) {
              hReTime50->Fill(m, ph2.pt());
              if (std::abs(ph2.time - timeOffset) < 30.e-9) {
                hReTime30->Fill(m, ph2.pt());
                if (std::abs(ph2.time - timeOffset) < 12.5e-9) {
                  hReTime12->Fill(m, ph2.pt());
                }
              }
            }
          }
        }
        if (std::abs(ph2.time - timeOffset) < 12.5e-9) { // strict cut on first photon
          if (std::abs(ph1.time - timeOffset) < 100.e-9) {
            hReTime100->Fill(m, ph1.pt());
            if (std::abs(ph1.time - timeOffset) < 50.e-9) {
              hReTime50->Fill(m, ph1.pt());
              if (std::abs(ph1.time - timeOffset) < 30.e-9) {
                hReTime30->Fill(m, ph1.pt());
                if (std::abs(ph1.time - timeOffset) < 12.5e-9) {
                  hReTime12->Fill(m, ph1.pt());
                }
              }
            }
          }
        }

        bool isPi0 = false;
        if constexpr (isMC) { // test parent
          int cp = commonParentPDG(ph1.label, ph2.label, mcPart);
          if (cp != 0) {
            hSignalAll->Fill(m, pt, w);
            if (cp == 111) {
              isPi0 = true;
              hPi0SignalAll->Fill(m, pt, w);
            }
          }
        }

        if (ph1.isCPVOK() && ph2.isCPVOK()) {
          hReCPV->Fill(m, pt, w);
          if (isPi0) {
            hPi0SignalCPV->Fill(m, pt, w);
          }
        }
        if (ph1.isDispOK() && ph2.isDispOK()) {
          hReDisp->Fill(m, pt, w);
          if (isPi0) {
            hPi0SignalDisp->Fill(m, pt, w);
          }
          if (ph1.isCPVOK() && ph2.isCPVOK()) {
            hReBoth->Fill(m, pt, w);
            if (isPi0) {
              hPi0SignalBoth->Fill(m, pt, w);
            }
          }
        }
      }

      // Add photon to stack
      mCurEvent.emplace_back(ph1);
    }

    // Mixed
    for (const auto& ph1 : mCurEvent) {
      for (const auto& mixEvent : mMixedEvents[mixedEventBin]) {
        for (const auto& ph2 : mixEvent) {
          double m = std::pow(ph1.e + ph2.e, 2) - std::pow(ph1.px + ph2.px, 2) -
                     std::pow(ph1.py + ph2.py, 2) - std::pow(ph1.pz + ph2.pz, 2);
          if (m > 0) {
            m = std::sqrt(m);
          }
          double pt = std::sqrt(std::pow(ph1.px + ph2.px, 2) +
                                std::pow(ph1.py + ph2.py, 2));
          int modComb = moduleCombination(ph1.mod, ph2.mod);
          double w = 1.;
          if constexpr (isMC) { // correct MC energy
            w = tofCutEff(ph1.e) * tofCutEff(ph2.e);
          }
          hMiMod->Fill(m, pt, modComb, w);
          hMiAll->Fill(m, pt, w);
          hMiOneAll->Fill(m, ph1.pt(), w);
          hMiOneAll->Fill(m, ph2.pt(), w);
          if (ph1.isCPVOK()) {
            hMiOneCPV->Fill(m, ph1.pt(), w);
          }
          if (ph2.isCPVOK()) {
            hMiOneCPV->Fill(m, ph2.pt(), w);
          }
          if (ph1.isDispOK()) {
            hMiOneDisp->Fill(m, ph1.pt(), w);
            if (ph1.isCPVOK()) {
              hMiOneBoth->Fill(m, ph1.pt(), w);
            }
          }
          if (ph2.isDispOK()) {
            hMiOneDisp->Fill(m, ph2.pt(), w);
            if (ph2.isCPVOK()) {
              hMiOneBoth->Fill(m, ph2.pt(), w);
            }
          }
          if (ph1.isCPVOK() && ph2.isCPVOK()) {
            hMiCPV->Fill(m, pt, w);
          }
          if (ph1.isDispOK() && ph2.isDispOK()) {
            hMiDisp->Fill(m, pt, w);
            if (ph1.isCPVOK() && ph2.isCPVOK()) {
              hMiBoth->Fill(m, pt, w);
            }
          }
        }
      }
    }

    // Fill events to store and remove oldest to keep buffer size
    if (mCurEvent.size() > 0) {
      mMixedEvents[mixedEventBin].emplace_back(mCurEvent);
      if (mMixedEvents[mixedEventBin].size() > static_cast<size_t>(nMixedEvents)) {
        mMixedEvents[mixedEventBin].pop_front();
      }
    }
  }
  void processBC(aod::BCs::iterator const& /*bc*/,
                 aod::CaloAmbiguousClusters const& clusters)
  {
    bool isSelected = true;
    mixedEventBin = 0;

    mHistManager.fill(HIST("eventsBC"), 0.);
    double vtxZ = 0; // no vtx info
    int mult = 1.;   // multiplicity TODO!!!
    mixedEventBin = findMixedEventBin(vtxZ, mult);

    if (!isSelected) {
      return;
    }

    // Fill invariant mass distributions
    mCurEvent.clear();
    for (const auto& clu : clusters) {
      // Fill QC histos
      if (fillQC) {
        mHistManager.fill(HIST("hM02Clu"), clu.e(), clu.m02());
        mHistManager.fill(HIST("hM20Clu"), clu.e(), clu.m20());
        mHistManager.fill(HIST("hNcellClu"), clu.e(), clu.ncell(), clu.mod());
        mHistManager.fill(HIST("cluETime"), clu.e(), clu.time(), clu.mod());
      }
      if (clu.e() < minCluE ||
          clu.ncell() < minCluNcell ||
          clu.time() > maxCluTime || clu.time() < minCluTime ||
          clu.m02() < minM02) {
        continue;
      }
      if (fillQC) {
        mHistManager.fill(HIST("cluSp"), clu.e(), clu.mod());
        if (clu.e() > minOccE) {
          mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
          if (clu.trackdist() > 2.) {
            mHistManager.fill(HIST("cluCPVOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpCPV"), clu.e(), clu.mod());
            if (testLambda(clu.e(), clu.m02(), clu.m20())) {
              mHistManager.fill(HIST("cluBothOcc"), clu.x(), clu.z(), clu.mod());
              mHistManager.fill(HIST("cluSpBoth"), clu.e(), clu.mod());
            }
          }
          if (testLambda(clu.e(), clu.m02(), clu.m20())) {
            mHistManager.fill(HIST("cluDispOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpDisp"), clu.e(), clu.mod());
          }
        }
      }

      int mcLabel = -1;
      double enCorr = 1;
      if (isMC) { // correct MC energy
        enCorr = nonlinearity(clu.e());
      }
      Photon ph1(clu.px() * enCorr, clu.py() * enCorr, clu.pz() * enCorr, clu.e() * enCorr, clu.time(), clu.mod(), testLambda(clu.e(), clu.m02(), clu.m20()), clu.trackdist() > cpvCut, mcLabel);

      // Mix with other photons added to stack
      for (const auto& ph2 : mCurEvent) {
        double m = std::pow(ph1.e + ph2.e, 2) - std::pow(ph1.px + ph2.px, 2) -
                   std::pow(ph1.py + ph2.py, 2) - std::pow(ph1.pz + ph2.pz, 2);
        if (m > 0) {
          m = std::sqrt(m);
        }
        double pt = std::sqrt(std::pow(ph1.px + ph2.px, 2) +
                              std::pow(ph1.py + ph2.py, 2));
        int modComb = moduleCombination(ph1.mod, ph2.mod);
        hReMod->Fill(m, pt, modComb);
        hReAll->Fill(m, pt);

        if (ph1.isCPVOK() && ph2.isCPVOK()) {
          hReCPV->Fill(m, pt);
        }
        if (ph1.isDispOK() && ph2.isDispOK()) {
          hReDisp->Fill(m, pt);
          if (ph1.isCPVOK() && ph2.isCPVOK()) {
            hReBoth->Fill(m, pt);
          }
        }
      }

      // Add photon to stack
      mCurEvent.emplace_back(ph1);
    }

    // Mixed
    for (const auto& ph1 : mCurEvent) {
      for (const auto& mixEvent : mMixedEvents[mixedEventBin]) {
        for (const auto& ph2 : mixEvent) {
          double m = std::pow(ph1.e + ph2.e, 2) - std::pow(ph1.px + ph2.px, 2) -
                     std::pow(ph1.py + ph2.py, 2) - std::pow(ph1.pz + ph2.pz, 2);
          if (m > 0) {
            m = std::sqrt(m);
          }
          double pt = std::sqrt(std::pow(ph1.px + ph2.px, 2) +
                                std::pow(ph1.py + ph2.py, 2));
          int modComb = moduleCombination(ph1.mod, ph2.mod);
          hMiMod->Fill(m, pt, modComb);
          hMiAll->Fill(m, pt);
          if (ph1.isCPVOK() && ph2.isCPVOK()) {
            hMiCPV->Fill(m, pt);
          }
          if (ph1.isDispOK() && ph2.isDispOK()) {
            hMiDisp->Fill(m, pt);
            if (ph1.isCPVOK() && ph2.isCPVOK()) {
              hMiBoth->Fill(m, pt);
            }
          }
        }
      }
    }

    // Fill events to store and remove oldest to keep buffer size
    if (mCurEvent.size() > 0) {
      mMixedEvents[mixedEventBin].emplace_back(mCurEvent);
      if (mMixedEvents[mixedEventBin].size() > static_cast<size_t>(nMixedEvents)) {
        mMixedEvents[mixedEventBin].pop_front();
      }
    }
  }
  PROCESS_SWITCH(PhosPi0, processBC, "processBC", false);

  //_____________________________________________________________________________
  int moduleCombination(int m1, int m2)
  {
    // enumerates possible module combinations
    // (1,1)=0, (2,2)=1, (3,3)=2, (4,4)=3, (1,2)=(2,1)=4, (2,3)=(3,2)=5, (3,4)=(4,3)=6, (1,3)=(3,1)=7,
    // (2,4)=(4,2)=8, (1,4)=(4,1)=9
    int d = std::abs(m1 - m2);
    if (d == 0) {
      return m1 - 1;
    }
    if (d == 1) {
      return 3 + std::min(m1, m2);
    }
    if (d == 2) {
      return 6 + std::min(m1, m2);
    }
    return 9;
  }
  //_____________________________________________________________________________
  bool testLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    // Parameterizatino for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * std::exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * std::exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
  //_____________________________________________________________________________
  int findMixedEventBin(double vtxZ, double /*mult */)
  {
    // calculate index for event mixing
    const double zwidth = 1.; // Width of zvtx bin
    int res = static_cast<int>((vtxZ + 10.) / zwidth);

    if (res < 0)
      return 0;
    if (res >= kMaxMixBins)
      return kMaxMixBins - 1;
    return res;
  }
  //----------------------------------------
  int commonParentPDG(int lab1, int lab2, aod::McParticles const* mcParticles)
  {
    // Tests if two labels contain common ancestor
    // return 0 if no common parent
    //        PDG code if found
    int iparent1 = lab1;
    while (iparent1 > -1) {
      int iparent2 = lab2;
      while (iparent2 > -1) {
        if (iparent1 == iparent2) {
          return mcParticles->iteratorAt(iparent1).pdgCode();
        }
        auto parent2 = mcParticles->iteratorAt(iparent2);
        if (parent2.mothersIds().size() == 0 || parent2.pdgCode() == 21 || std::abs(parent2.pdgCode()) < 11 || std::abs(parent2.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
          break;
        }
        iparent2 = parent2.mothersIds()[0];
      }
      auto parent1 = mcParticles->iteratorAt(iparent1);
      if (parent1.mothersIds().size() == 0 || parent1.pdgCode() == 21 || std::abs(parent1.pdgCode()) < 11 || std::abs(parent1.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
        return 0;
      }
      iparent1 = parent1.mothersIds()[0];
    }
    return 0; // nothing found
  }
  double nonlinearity(double e)
  {
    return nonlinA + nonlinB * std::exp(-e / nonlinC);
  }
  double tofCutEff(double en)
  {
    if (tofEffParam == 0) {
      return 1.;
    }
    if (tofEffParam == 1) { // Run2 100 ns
      // parameterization 01.08.2020
      if (en > 1.1)
        en = 1.1;
      if (en < 0.11)
        en = 0.11;
      return std::exp((-1.15295e+05 + 2.26754e+05 * en - 1.26063e+05 * en * en + en * en * en) /
                      (1. - 3.16443e+05 * en + 3.68044e+06 * en * en + en * en * en));
    }
    if (tofEffParam == 2) { // Run2 30 ns
      if (en > 1.6)
        en = 1.6;
      return 1. / (1. + std::exp((4.83230e+01 - 8.89758e+01 * en + 1.10897e+03 * en * en - 5.73755e+03 * en * en * en -
                                  1.43777e+03 * en * en * en * en) /
                                 (1. - 1.23667e+02 * en + 1.07255e+03 * en * en + 5.87221e+02 * en * en * en)));
    }
    if (tofEffParam == 2) { // Run2 12.5 ns
      if (en < 4.6) {
        return std::exp(3.64952e-03 *
                        (-5.80032e+01 - 1.53442e+02 * en + 1.30994e+02 * en * en + -3.53094e+01 * en * en * en + en * en * en * en) /
                        (-7.75638e-02 + 8.64761e-01 * en + 1.22320e+00 * en * en - 1.00177e+00 * en * en * en + en * en * en * en));
      } else {
        return 0.63922783 * (1. - 1.63273e-01 * std::tanh((en - 7.94528e+00) / 1.28997e+00)) *
               (-4.39257e+00 * en + 2.25503e+00 * en * en + en * en * en) / (2.37160e+00 * en - 6.93786e-01 * en * en + en * en * en);
      }
    }
    return 1.;
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<PhosPi0>(cfgc)};
}
