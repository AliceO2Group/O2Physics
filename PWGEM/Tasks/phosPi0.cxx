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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
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

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

/// \struct PHOS pi0 analysis
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since Nov, 2022
///

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosPi0 {
  Configurable<bool> mIsMC{"isMC", false, "to fill MC histograms"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<float> mMinCluE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"};
  Configurable<float> mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"};
  Configurable<float> mMaxCluTime{"maxCluTime", 25.e-9, "Max. cluster time"};
  Configurable<int> mMinCluNcell{"minCluNcell", 2, "min cells in cluster"};
  Configurable<float> mMinM02{"minM02", 0.2, "Min disp M02 cut"};
  Configurable<float> mCPVCut{"CPVCut", 2., "Min distance to track"};
  Configurable<int> mNMixedEvents{"nMixedEvents", 10, "number of events to mix"};
  Configurable<bool> mSelectOneCollPerBC{"selectOneColPerBC", true, "skip multiple coll. per bc"};
  Configurable<bool> mFillQC{"fillQC", true, "Fill QC histos"};
  Configurable<float> mOccE{"minOccE", 0.5, "Min. cluster energy of occupancy plots"};

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using SelCollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using BCsWithBcSels = soa::Join<aod::BCs, aod::BcSels>;
  using mcClusters = soa::Join<aod::CaloClusters, aod::PHOSCluLabels>;
  using mcAmbClusters = soa::Join<aod::CaloAmbiguousClusters, aod::PHOSAmbCluLabels>;

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry mHistManager{"phosPi0Histograms"};

  // class to keep photon candidate parameters
  class photon
  {
   public:
    photon() = default;
    photon(double x, double y, double z, double ee, int m, bool isDispOK, bool isCPVOK, int mcLabel) : px(x), py(y), pz(z), e(ee), mod(m), mPID(isDispOK << 1 | isCPVOK << 2), label(mcLabel) {}
    ~photon() = default;

    bool isCPVOK() { return (mPID >> 2) & 1; }
    bool isDispOK() { return (mPID >> 1) & 1; }

   public:
    double px = 0.; // px
    double py = 0.; // py
    double pz = 0.; // pz
    double e = 0.;  // energy
    int mod = 0;    // module
    int mPID = 0;   // store PID bits
    int label = -1; // label of MC particle
  };

  int mixedEventBin = 0; // Which list of Mixed use for mixing
  std::vector<photon> mCurEvent;
  static constexpr int nMaxMixBins = 20; // maximal number of kinds of events for mixing
  std::array<std::deque<std::vector<photon>>, nMaxMixBins> mMixedEvents;
  std::array<std::deque<std::vector<photon>>, nMaxMixBins> mAmbMixedEvents;

  int mPrevMCColId = -1; // mark MC collissions already scanned
  // fast access to histos
  TH1* hColl;
  TH3 *hReMod, *hMiMod;
  TH2 *hReAll, *hReDisp, *hReCPV, *hReBoth, *hSignalAll, *hPi0SignalAll, *hPi0SignalCPV, *hPi0SignalDisp,
    *hPi0SignalBoth, *hMiAll, *hMiDisp, *hMiCPV, *hMiBoth;

  std::vector<double> pt = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2,
                            1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,
                            6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 26., 28.,
                            30., 34., 38., 42., 46., 50., 55., 60., 70., 80., 90., 100., 110., 120., 150.};

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS pi0 analysis task ...";

    const AxisSpec
      ptAxis{pt, "p_{T} (GeV/c)"},
      mggAxis{625, 0., 1.25, "m_{#gamma#gamma} (GeV/c^{2})"},
      timeAxis{100, -500.e-9, 500.e-9, "t (s)"},
      M02Axis{100, 0., 20., "M02 (cm^{2})"},
      M20Axis{100, 0., 20., "M20 (cm^{2})"},
      nCellsAxis{100, 0., 100., "N_{cell}"},
      zAxis{56, -63., 63., "z", "z (cm)"},
      phiAxis{64, -72., 72., "x", "x (cm)"},
      modAxis{4, 1., 5., "module", "Module"},
      multAxis{100, 0., 100.},
      vertexAxis{100, -20., 20., "z", "z (cm)"},
      modCombAxis{10, 0., 10.},
      centAxis{10, 0., 10.},
      centralityAxis{100, 0., 100., "centrality", "centrality"};

    hColl = std::get<std::shared_ptr<TH1>>(mHistManager.add("eventsCol", "Number of events", HistType::kTH1F, {{9, 0., 9.}})).get();
    hColl->GetXaxis()->SetBinLabel(1, "All");
    hColl->GetXaxis()->SetBinLabel(2, "T0a||T0c");
    hColl->GetXaxis()->SetBinLabel(3, "T0a&&T0c");
    hColl->GetXaxis()->SetBinLabel(4, "kTVXinPHOS");
    hColl->GetXaxis()->SetBinLabel(5, "kIsTriggerTVX");
    hColl->GetXaxis()->SetBinLabel(6, "PHOSClu");
    hColl->GetXaxis()->SetBinLabel(7, "PHOSClu&&kTVXinPHOS");
    hColl->GetXaxis()->SetBinLabel(8, "Accepted");

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

    if (mFillQC) {
      // QC histograms for normal collisions
      mHistManager.add("cluSp", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("cluSpDisp", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("cluSpCPV", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("cluSpBoth", "Cluster spectrum per module", HistType::kTH2F, {ptAxis, modAxis});
      mHistManager.add("hM02Clu", "(M02,M20) in clusters", HistType::kTH2F, {ptAxis, M02Axis});
      mHistManager.add("hM20Clu", "(M02,M20) in clusters", HistType::kTH2F, {ptAxis, M20Axis});
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

    if (mIsMC) {
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
    if (mIsMC) {
      mHistManager.add("hMCPi0SpAll", "pi0 spectrum inclusive", HistType::kTH1F, {ptAxis});
      mHistManager.add("hMCPi0SpPrim", "pi0 spectrum Primary", HistType::kTH1F, {ptAxis});
      mHistManager.add("hMCPi0RapPrim", "pi0 rapidity primary", HistType::kTH1F, {{100, -1., 1., "Rapidity"}});
      mHistManager.add("hMCPi0PhiPrim", "pi0 phi primary", HistType::kTH1F, {{100, 0., TMath::TwoPi(), "#phi (rad)"}});
      mHistManager.add("hMCPi0SecVtx", "pi0 secondary", HistType::kTH2F, {{100, 0., 500., "R (cm)"}, {100, -TMath::Pi(), TMath::Pi(), "#phi (rad)"}});
    }
  }

  /// \brief Process PHOS data
  void processData(SelCollisions::iterator const& col,
                   aod::CaloClusters const& clusters)
  {
    aod::McParticles const* mcPart = nullptr;
    scanAll<false>(col, clusters, mcPart);
  }
  PROCESS_SWITCH(phosPi0, processData, "processData", true);
  void processMC(SelCollisionsMC::iterator const& col,
                 mcClusters const& clusters,
                 aod::McParticles const& mcPart,
                 aod::McCollisions const& /*mcCol*/)
  {
    scanAll<true>(col, clusters, &mcPart);
  }
  PROCESS_SWITCH(phosPi0, processMC, "processMC", false);

  template <bool isMC, typename TCollision, typename TClusters>
  void scanAll(TCollision& col,
               TClusters& clusters,
               aod::McParticles const* mcPart)
  {
    bool isColSelected = false;
    mixedEventBin = 0;

    hColl->Fill(0.5);
    if (col.selection_bit(kIsBBT0A) || col.selection_bit(kIsBBT0C)) {
      hColl->Fill(1.5);
    }
    if (col.selection_bit(kIsBBT0A) && col.selection_bit(kIsBBT0C)) {
      hColl->Fill(2.5);
    }
    if (col.alias_bit(kTVXinPHOS)) {
      hColl->Fill(3.5);
    }
    if (col.selection_bit(kIsTriggerTVX)) {
      hColl->Fill(4.5);
    }
    if (clusters.size() > 0) {
      hColl->Fill(5.5);
      if (col.alias_bit(kTVXinPHOS)) {
        hColl->Fill(6.5);
      }
    }
    isColSelected = false;
    if constexpr (isMC) {
      isColSelected = (col.selection_bit(kIsTriggerTVX) && (clusters.size() > 0));
    } else {
      isColSelected = col.alias_bit(mEvSelTrig);
    }
    double vtxZ = col.posZ();
    mHistManager.fill(HIST("vertex"), vtxZ);
    int mult = 1.; // multiplicity TODO!!!
    mixedEventBin = findMixedEventBin(vtxZ, mult);

    if (!isColSelected) {
      return;
    }
    hColl->Fill(7.5);

    // Fill MC distributions
    // pion rapidity, pt, phi
    // secondary pi0s
    if constexpr (isMC) {
      // check current collision Id for clusters
      int cluMcBCId = -1;
      for (auto clu : clusters) {
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
          for (auto part : *mcPart) {
            if (part.mcCollision().bcId() != cluMcBCId) {
              continue;
            }
            if (part.pdgCode() == 111) {
              double r = sqrt(pow(part.vx(), 2) + pow(part.vy(), 2));
              if (r < 0.5) {
                mHistManager.fill(HIST("hMCPi0RapPrim"), part.y());
              }
              if (abs(part.y()) < .5) {
                double pt = part.pt();
                mHistManager.fill(HIST("hMCPi0SpAll"), pt);
                double phiVtx = atan2(part.vy(), part.vx());
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
      if (mFillQC) {
        mHistManager.fill(HIST("hM02Clu"), clu.e(), clu.m02());
        mHistManager.fill(HIST("hM20Clu"), clu.e(), clu.m20());
        mHistManager.fill(HIST("hNcellClu"), clu.e(), clu.ncell(), clu.mod());
        mHistManager.fill(HIST("cluETime"), clu.e(), clu.time(), clu.mod());
      }
      if (clu.e() < mMinCluE ||
          clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime ||
          clu.m02() < mMinM02) {
        continue;
      }
      if (mFillQC) {
        mHistManager.fill(HIST("cluSp"), clu.e(), clu.mod());
        if (clu.e() > mOccE) {
          mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
          if (clu.trackdist() > 2.) {
            mHistManager.fill(HIST("cluCPVOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpCPV"), clu.e(), clu.mod());
            if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
              mHistManager.fill(HIST("cluBothOcc"), clu.x(), clu.z(), clu.mod());
              mHistManager.fill(HIST("cluSpBoth"), clu.e(), clu.mod());
            }
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
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
      photon ph1(clu.px(), clu.py(), clu.pz(), clu.e(), clu.mod(), TestLambda(clu.e(), clu.m02(), clu.m20()), clu.trackdist() > mCPVCut, mcLabel);
      // Mix with other photons added to stack
      for (auto ph2 : mCurEvent) {
        double m = pow(ph1.e + ph2.e, 2) - pow(ph1.px + ph2.px, 2) -
                   pow(ph1.py + ph2.py, 2) - pow(ph1.pz + ph2.pz, 2);
        if (m > 0) {
          m = sqrt(m);
        }
        double pt = sqrt(pow(ph1.px + ph2.px, 2) +
                         pow(ph1.py + ph2.py, 2));
        int modComb = ModuleCombination(ph1.mod, ph2.mod);
        hReMod->Fill(m, pt, modComb);
        hReAll->Fill(m, pt);
        bool isPi0 = false;
        if constexpr (isMC) { // test parent
          int cp = commonParentPDG(ph1.label, ph2.label, mcPart);
          if (cp != 0) {
            hSignalAll->Fill(m, pt);
            if (cp == 111) {
              isPi0 = true;
              hPi0SignalAll->Fill(m, pt);
            }
          }
        }

        if (ph1.isCPVOK() && ph2.isCPVOK()) {
          hReCPV->Fill(m, pt);
          if (isPi0) {
            hPi0SignalCPV->Fill(m, pt);
          }
        }
        if (ph1.isDispOK() && ph2.isDispOK()) {
          hReDisp->Fill(m, pt);
          if (isPi0) {
            hPi0SignalDisp->Fill(m, pt);
          }
          if (ph1.isCPVOK() && ph2.isCPVOK()) {
            hReBoth->Fill(m, pt);
            if (isPi0) {
              hPi0SignalBoth->Fill(m, pt);
            }
          }
        }
      }

      // Add photon to stack
      mCurEvent.emplace_back(ph1);
    }

    // Mixed
    for (auto ph1 : mCurEvent) {
      for (auto mixEvent : mMixedEvents[mixedEventBin]) {
        for (auto ph2 : mixEvent) {
          double m = pow(ph1.e + ph2.e, 2) - pow(ph1.px + ph2.px, 2) -
                     pow(ph1.py + ph2.py, 2) - pow(ph1.pz + ph2.pz, 2);
          if (m > 0) {
            m = sqrt(m);
          }
          double pt = sqrt(pow(ph1.px + ph2.px, 2) +
                           pow(ph1.py + ph2.py, 2));
          int modComb = ModuleCombination(ph1.mod, ph2.mod);
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
      if (mMixedEvents[mixedEventBin].size() > static_cast<size_t>(mNMixedEvents)) {
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
      if (mFillQC) {
        mHistManager.fill(HIST("hM02Clu"), clu.e(), clu.m02());
        mHistManager.fill(HIST("hM20Clu"), clu.e(), clu.m20());
        mHistManager.fill(HIST("hNcellClu"), clu.e(), clu.ncell(), clu.mod());
        mHistManager.fill(HIST("cluETime"), clu.e(), clu.time(), clu.mod());
      }
      if (clu.e() < mMinCluE ||
          clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime ||
          clu.m02() < mMinM02) {
        continue;
      }
      if (mFillQC) {
        mHistManager.fill(HIST("cluSp"), clu.e(), clu.mod());
        if (clu.e() > mOccE) {
          mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
          if (clu.trackdist() > 2.) {
            mHistManager.fill(HIST("cluCPVOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpCPV"), clu.e(), clu.mod());
            if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
              mHistManager.fill(HIST("cluBothOcc"), clu.x(), clu.z(), clu.mod());
              mHistManager.fill(HIST("cluSpBoth"), clu.e(), clu.mod());
            }
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
            mHistManager.fill(HIST("cluDispOcc"), clu.x(), clu.z(), clu.mod());
            mHistManager.fill(HIST("cluSpDisp"), clu.e(), clu.mod());
          }
        }
      }

      int mcLabel = -1;
      photon ph1(clu.px(), clu.py(), clu.pz(), clu.e(), clu.mod(), TestLambda(clu.e(), clu.m02(), clu.m20()), clu.trackdist() > mCPVCut, mcLabel);
      // Mix with other photons added to stack
      for (auto ph2 : mCurEvent) {
        double m = pow(ph1.e + ph2.e, 2) - pow(ph1.px + ph2.px, 2) -
                   pow(ph1.py + ph2.py, 2) - pow(ph1.pz + ph2.pz, 2);
        if (m > 0) {
          m = sqrt(m);
        }
        double pt = sqrt(pow(ph1.px + ph2.px, 2) +
                         pow(ph1.py + ph2.py, 2));
        int modComb = ModuleCombination(ph1.mod, ph2.mod);
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
    for (auto ph1 : mCurEvent) {
      for (auto mixEvent : mMixedEvents[mixedEventBin]) {
        for (auto ph2 : mixEvent) {
          double m = pow(ph1.e + ph2.e, 2) - pow(ph1.px + ph2.px, 2) -
                     pow(ph1.py + ph2.py, 2) - pow(ph1.pz + ph2.pz, 2);
          if (m > 0) {
            m = sqrt(m);
          }
          double pt = sqrt(pow(ph1.px + ph2.px, 2) +
                           pow(ph1.py + ph2.py, 2));
          int modComb = ModuleCombination(ph1.mod, ph2.mod);
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
      if (mMixedEvents[mixedEventBin].size() > static_cast<size_t>(mNMixedEvents)) {
        mMixedEvents[mixedEventBin].pop_front();
      }
    }
  }
  PROCESS_SWITCH(phosPi0, processBC, "processBC", false);

  //_____________________________________________________________________________
  int ModuleCombination(int m1, int m2)
  {
    // enumerates possible module combinations
    // (1,1)=0, (2,2)=1, (3,3)=2, (4,4)=3, (1,2)=(2,1)=4, (2,3)=(3,2)=5, (3,4)=(4,3)=6, (1,3)=(3,1)=7,
    // (2,4)=(4,2)=8, (1,4)=(4,1)=9
    int d = TMath::Abs(m1 - m2);
    if (d == 0) {
      return m1 - 1;
    }
    if (d == 1) {
      return 3 + TMath::Min(m1, m2);
    }
    if (d == 2) {
      return 6 + TMath::Min(m1, m2);
    }
    return 9;
  }
  //_____________________________________________________________________________
  bool TestLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    // Parameterizatino for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * TMath::Exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * TMath::Exp(-0.390730 * pt);

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
    if (res >= nMaxMixBins)
      return nMaxMixBins - 1;
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
        if (parent2.mothersIds().size() == 0 || parent2.pdgCode() == 21 || abs(parent2.pdgCode()) < 11 || abs(parent2.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
          break;
        }
        iparent2 = parent2.mothersIds()[0];
      }
      auto parent1 = mcParticles->iteratorAt(iparent1);
      if (parent1.mothersIds().size() == 0 || parent1.pdgCode() == 21 || abs(parent1.pdgCode()) < 11 || abs(parent1.pdgCode()) > 5000) { // no parents, parent not quark/gluon, strings
        return 0;
      }
      iparent1 = parent1.mothersIds()[0];
    }
    return 0; // nothing found
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosPi0>(cfgc)};
}
