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
//
/// \brief (Anti)triton analysis
/// \author Nikita Gladin (ngladin@cern.ch)
/// Based on PWGCF/Femto3D/DataModel/singletrackselector.h

#include "PWGCF/Femto3D/Core/femto3dPairTask.h"
#include "PWGCF/Femto3D/DataModel/PIDutils.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponseITS.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AntitritonAnalysis {
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryTrue{"registryTrue", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> isMC{"isMC", false, ""};

  Configurable<bool> removeSameBunchPileup{"removeSameBunchPileup", false, ""};
  Configurable<bool> requestGoodZvtxFT0vsPV{"requestGoodZvtxFT0vsPV", false, ""};
  Configurable<bool> requestVertexITSTPC{"requestVertexITSTPC", false, ""};
  Configurable<int> requestVertexTOFmatched{"requestVertexTOFmatched", 0, "0 -> no selection; 1 -> vertex is matched to TOF or TRD; 2 -> matched to both;"};
  Configurable<bool> requestNoCollInTimeRangeStandard{"requestNoCollInTimeRangeStandard", false, ""};
  Configurable<bool> requestIsGoodITSLayersAll{"requestIsGoodITSLayersAll", false, "cut time intervals with dead ITS staves"};
  Configurable<std::pair<float, float>> irCut{"irCut", std::pair<float, float>{0.f, 100.f}, "[min., max.] IR range to keep events within"};
  Configurable<std::pair<int, int>> occupancyCut{"occupancyCut", std::pair<int, int>{0, 10000}, "[min., max.] occupancy range to keep events within"};

  Configurable<float> vertexZ{"vertexZ", 10.0, "abs vertexZ value limit"};
  Configurable<float> minP{"minP", 0.0, "lower mometum limit"};
  Configurable<float> maxP{"maxP", 100.0, "upper mometum limit"};
  Configurable<float> eta{"eta", 100.0, "abs eta value limit"};
  Configurable<std::vector<float>> dcaXY{"dcaXY", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaXY value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<std::vector<float>> dcaZ{"dcaZ", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaZ value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<int16_t> minTpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<float> maxTpcFractionSharedCls{"maxTpcFractionSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> minItsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters for a track"};
  Configurable<float> itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters for a track"};
  Configurable<int> particlePDG{"particlePDG", o2::constants::physics::Pdg::kTriton, "PDG code of a particle to perform PID for"};
  Configurable<std::vector<float>> tpcNSigma{"tpcNSigma", std::vector<float>{-4.0f, 4.0f}, "Nsigma range in TPC before the TOF is used"};
  Configurable<std::vector<float>> itsNSigma{"itsNSigma", std::vector<float>{-10.0f, 10.0f}, "Nsigma range in ITS to use along with TPC"};
  Configurable<float> pidTrshld{"pidTrshld", 0.0, "value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<float> massCuttrshld{"massCuttrshld", 0.0, "value of momentum from which the mass cut is done with TOF"};
  Configurable<std::vector<float>> itsNsigNL{"itsNsigNL", std::vector<float>{-100.0f, 100.0f}, "Nsigma p-dependent cut in ITS ; formula: [0] + [1]*pT"};

  Configurable<std::vector<float>> tofNSigma{"tofNSigma", std::vector<float>{-4.0f, 4.0f}, "Nsigma range in TOF"};
  Configurable<std::vector<float>> tpcNSigmaResidual{"tpcNSigmaResidual", std::vector<float>{-5.0f, 5.0f}, "residual TPC Nsigma cut to use with the TOF"};

  Configurable<std::vector<int>> particlePDGtoReject{"particlePDGtoReject", std::vector<int>{PDG_t::kProton}, "PDG codes of perticles that will be rejected with TPC"};
  Configurable<std::vector<float>> rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for particles specified with PDG to be rejected"};
  Configurable<std::vector<float>> rejectWithinNsigmaTPC{"rejectWithinNsigmaTPC", std::vector<float>{-0.0f, 0.0f}, "TPC rejection Nsigma range for particles specified with PDG to be rejected"};

  Configurable<std::pair<float, float>> centCut{"centCut", std::pair<float, float>{0.f, 100.f}, "[min., max.] centrality range to keep tracks within"};

  Configurable<std::vector<float>> dcaBinning{"dcaBinning", std::vector<float>{501, 0.5f, 1}, "setup for variable binning (geometric progression is used): 1st (int) -- N_bins (must be odd, otherwise will be increased by 1); 2nd (float) -- abs value of the edge of axises in histos (-2nd, +2nd); 3d (int) -- desired ratio between w_bin at the edges and at 0;"};
  Configurable<std::vector<float>> tofMass{"tofMass", std::vector<float>{3.0f, 3.0f}, "TOF Mass bounds"};

  std::pair<int, std::vector<float>> ITScuts;
  std::pair<int, std::vector<float>> TPCcuts;
  std::pair<int, std::vector<float>> TOFcuts;

  static constexpr float kMassTriton = o2::constants::physics::MassTriton; // GeV/c^2

  Filter pFilter = o2::aod::singletrackselector::p > minP&& o2::aod::singletrackselector::p < maxP;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < eta;
  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= minTpcNClsFound && o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) < tpcChi2NCl && o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) > tpcCrossedRowsOverFindableCls;
  Filter itsTrkFilter = o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) < itsChi2NCl;
  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < vertexZ;

  struct HistSet {
    std::shared_ptr<TH1> eta;
    std::shared_ptr<TH2> eta_to_y;
    std::shared_ptr<TH1> y;
    std::shared_ptr<TH1> phi;
    std::shared_ptr<TH1> p;
    std::shared_ptr<TH1> pt;
    std::shared_ptr<TH2> dcaxy_to_p;
    std::shared_ptr<TH2> dcaxy_to_pt;
    std::shared_ptr<TH2> dcaz_to_p;
    std::shared_ptr<TH2> dcaz_to_pt;
    std::shared_ptr<TH1> TPCClusters;
    std::shared_ptr<TH1> ITSClusters;
    std::shared_ptr<TH2> TPCSignal;
    std::shared_ptr<TH2> TOFSignal;
    std::shared_ptr<TH2> mTOF;
    std::shared_ptr<TH2> nsigmaITS;
    std::shared_ptr<TH2> nsigmaTPC;
    std::shared_ptr<TH2> nsigmaTOF;
    std::shared_ptr<TH2> TOFoverTPC;
    std::shared_ptr<TH2> innerParamToP;
    std::shared_ptr<TH2> PtGenPtRec;
    std::shared_ptr<TH3> dcaxy_dcaz_to_pt;
    std::shared_ptr<TH2> dcaxy_to_dcaz;
    std::shared_ptr<TH1> origin;
  };

  static inline const std::vector<std::string> stageDirs = {
    "init",
    "No_Cuts",
    "ITS_Cuts",
    "TPC_Cuts",
    "TOF_Cuts",
    "ITSTPC_Cut",
    "PID_Cuts",
    "TOFmass_Cut",
    "Rejection_Cut",
    "TPCLinear_Cut",
  };

  // signDirs[0] = "t" for sign > 0 (triton), signDirs[1] = "at" for sign < 0 (antitriton).
  static inline const std::array<std::string, 2> signDirs = {"t", "at"};

  std::vector<std::array<HistSet, 2>> hSets;
  std::vector<std::array<HistSet, 2>> hSetsTrue;

  HistSet makeHistSet(HistogramRegistry& reg, const std::string& dir, int pdgForPid = 0)
  {
    HistSet h;

    h.eta = reg.add<TH1>((dir + "/eta").c_str(), "#eta;#eta;Entries", kTH1F, {{200, -2.5, 2.5}});
    h.y = reg.add<TH1>((dir + "/y").c_str(), "y;y;Entries", kTH1F, {{200, -2.5, 2.5}});
    h.eta_to_y = reg.add<TH2>((dir + "/eta_to_y").c_str(), "#eta;y;Entries", kTH2F, {{200, -2.5, 2.5}, {200, -2.5, 2.5}});
    h.phi = reg.add<TH1>((dir + "/phi").c_str(), "#phi;#phi;Entries", kTH1F, {{200, 0., 2. * M_PI}});
    h.p = reg.add<TH1>((dir + "/p").c_str(), "p;p (GeV/#it{c});Entries", kTH1F, {{100, 0., 5.}});
    h.pt = reg.add<TH1>((dir + "/pt").c_str(), "p_{T};p_{T} (GeV/#it{c});Entries", kTH1F, {{100, 0., 5.}});
    h.dcaxy_to_p = reg.add<TH2>((dir + "/dcaxy_to_p").c_str(), "DCA_{xy} vs p;p (GeV/#it{c});DCA_{xy} (cm)", kTH2F, {{100, 0., 5.0}, {501, -0.5, 0.5}});
    h.dcaxy_to_pt = reg.add<TH2>((dir + "/dcaxy_to_pt").c_str(), "DCA_{xy} vs p_{T};p_{T} (GeV/#it{c});DCA_{xy} (cm)", kTH2F, {{100, 0., 5.0}, {501, -0.5, 0.5}});
    h.dcaz_to_p = reg.add<TH2>((dir + "/dcaz_to_p").c_str(), "DCA_{z} vs p;p (GeV/#it{c});DCA_{z} (cm)", kTH2F, {{100, 0., 5.0}, {501, -0.5, 0.5}});
    h.dcaz_to_pt = reg.add<TH2>((dir + "/dcaz_to_pt").c_str(), "DCA_{z} vs p_{T};p_{T} (GeV/#it{c});DCA_{z} (cm)", kTH2F, {{100, 0., 5.0}, {501, -0.5, 0.5}});
    h.TPCClusters = reg.add<TH1>((dir + "/TPCClusters").c_str(), "TPC Clusters;N_{TPC clusters};Entries", kTH1F, {{163, -0.5, 162.5}});
    h.ITSClusters = reg.add<TH1>((dir + "/ITSClusters").c_str(), "ITS Clusters;N_{ITS clusters};Entries", kTH1F, {{10, -0.5, 9.5}});
    h.TPCSignal = reg.add<TH2>((dir + "/TPCSignal").c_str(), "TPC Signal;#it{p}_{inner} (GeV/#it{c});dE/dx in TPC (arbitrary units)", kTH2F, {{200, 0., 5.0}, {1000, 0., 1000.0}});
    h.TOFSignal = reg.add<TH2>((dir + "/TOFSignal").c_str(), "TOF Signal;#it{p} (GeV/#it{c});#beta", kTH2F, {{200, 0., 5.0}, {100, 0., 1.5}});
    h.mTOF = reg.add<TH2>((dir + "/mTOF").c_str(), "m_{TOF};p_{T} (GeV/#it{c});m_{TOF}(GeV/#it{c}^{2})", kTH2F, {{120, 0., 6.}, {1000, -1.5, 5.0}});
    h.innerParamToP = reg.add<TH2>((dir + "/ip_to_p").c_str(), ";#it{p} (GeV/#it{c});ip_to_p", kTH2F, {{100, 0., 5.0}, {100, 0.0, 2.0}});
    h.PtGenPtRec = reg.add<TH2>((dir + "/ptgenptrec").c_str(), ";#it{p_{T, rec}} (GeV/#it{c});#it{p_{T, rec}} - #it{p_{T, gen}}(GeV/#it{c})", kTH2F, {{100, 0., 5.0}, {100, -2.5, 2.5}});
    h.dcaxy_dcaz_to_pt = reg.add<TH3>(
      (dir + "/dcaxy_dcaz_to_pt").c_str(),
      "DCA_{xy} vs DCA_{z} vs p_{T};DCA_{xy} (cm);DCA_{z} (cm);p_{T} (GeV/#it{c})",
      kTH3F,
      {
        {101, -0.5, 0.5}, // x axis: DCA_xy
        {101, -0.5, 0.5}, // y axis: DCA_z
        {20, 0., 5.0}     // z axis: pT
      });
    h.dcaxy_to_dcaz = reg.add<TH2>(
      (dir + "/dcaxy_to_dcaz").c_str(),
      "DCA_{z} vs DCA_{xy};DCA_{xy} (cm);DCA_{z} (cm);Counts",
      kTH2F,
      {
        {501, -0.5, 0.5}, // x axis: DCAxy
        {501, -0.5, 0.5}  // y axis: DCAz
      });
    h.origin = reg.add<TH1>((dir + "/origin").c_str(), "MC origin;origin (-1: unmatched, 0: primary, 1: weak decay, 2: material);Entries", kTH1F, {{4, -1.5, 2.5}});

    if (pdgForPid != 0) {
      h.TOFoverTPC = reg.add<TH2>((dir + Form("/TOFoverTPC_PDG%i", pdgForPid)).c_str(), "n#sigma_{TOF};n#sigma_{TPC};", kTH2F, {{100, -10., 10.}, {100, -10., 10.}});
      h.nsigmaITS = reg.add<TH2>((dir + Form("/nsigmaITS_PDG%i", pdgForPid)).c_str(), Form("n#sigma_{ITS};p (GeV/#it{c});n#sigma_{ITS}"), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
      h.nsigmaTPC = reg.add<TH2>((dir + Form("/nsigmaTPC_PDG%i", pdgForPid)).c_str(), Form("n#sigma_{TPC};p (GeV/#it{c});n#sigma_{TPC}"), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
      h.nsigmaTOF = reg.add<TH2>((dir + Form("/nsigmaTOF_PDG%i", pdgForPid)).c_str(), Form("n#sigma_{TOF};p (GeV/#it{c});n#sigma_{TOF}"), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    }

    return h;
  }

  template <typename TrackType>
  static float getMTOF(const TrackType& track)
  {
    const float beta = track.beta();
    if (beta <= 0.f || beta >= 1.f) {
      return -1.f;
    }
    return track.p() * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  template <typename TrackType>
  void fillHistSet(HistSet& h, const TrackType& track, int pdg)
  {
    const float mass = kMassTriton;
    const auto y = static_cast<float>(RecoDecay::y(std::array{track.px(), track.py(), track.pz()}, mass));
    h.eta->Fill(track.eta());
    h.eta_to_y->Fill(track.eta(), y);
    h.y->Fill(y);
    h.phi->Fill(track.phi());
    h.p->Fill(track.p());
    h.pt->Fill(track.pt());
    h.dcaxy_to_p->Fill(track.p(), track.dcaXY());
    h.dcaxy_to_pt->Fill(track.pt(), track.dcaXY());
    h.dcaz_to_p->Fill(track.p(), track.dcaZ());
    h.dcaz_to_pt->Fill(track.pt(), track.dcaZ());
    h.dcaxy_dcaz_to_pt->Fill(track.dcaXY(), track.dcaZ(), track.pt());
    h.dcaxy_to_dcaz->Fill(track.dcaXY(), track.dcaZ());
    h.TPCClusters->Fill(track.tpcNClsFound());
    h.ITSClusters->Fill(track.itsNCls());
    h.nsigmaITS->Fill(track.p(), o2::aod::singletrackselector::getITSNsigma(track, pdg));
    h.nsigmaTPC->Fill(track.p(), o2::aod::singletrackselector::getTPCNsigma(track, pdg));
    h.nsigmaTOF->Fill(track.p(), o2::aod::singletrackselector::getTOFNsigma(track, pdg));
    h.TOFoverTPC->Fill(o2::aod::singletrackselector::getTPCNsigma(track, pdg), o2::aod::singletrackselector::getTOFNsigma(track, pdg));
  }

  template <bool HasMC, typename TrackType>
  void fillHistSetExtra(HistSet& h, const TrackType& track)
  {
    h.TPCSignal->Fill(track.tpcInnerParam(), track.tpcSignal());
    h.TOFSignal->Fill(track.p(), track.beta());
    h.mTOF->Fill(track.pt(), getMTOF(track));
    h.innerParamToP->Fill(track.p(), track.tpcInnerParam() / track.p());
    if constexpr (HasMC)
      h.PtGenPtRec->Fill(track.pt(), track.pt() - track.pt_MC());
  }

  template <bool FillExtra, bool HasMC, typename TrackType>
  void fillSet(HistSet& Set, HistSet& mcSet, const TrackType& track, int pdg)
  {
    fillHistSet(Set, track, pdg);
    if constexpr (FillExtra) {
      fillHistSetExtra<false>(Set, track);
    }

    if constexpr (HasMC) {
      if ((std::abs(track.pdgCode()) == pdg) && (track.pdgCode() / std::abs(track.pdgCode()) == track.sign())) {
        mcSet.origin->Fill(track.origin());
        fillHistSet(mcSet, track, pdg);
        if constexpr (FillExtra) {
          fillHistSetExtra<true>(mcSet, track);
        }
      }
    }
  }

  void init(o2::framework::InitContext&)
  {

    if (isMC.value)
      o2::aod::ITSResponse::setMCDefaultParameters(); // set MC parametrisation for the ITS PID

    ITScuts = std::make_pair(particlePDG.value, itsNSigma);
    TPCcuts = std::make_pair(particlePDG.value, tpcNSigma);
    TOFcuts = std::make_pair(particlePDG.value, tofNSigma);

    int N = static_cast<int>(dcaBinning.value[0]); // number of bins -- must be odd otherwise will be increased by 1
    if (N % 2 != 1) {
      N += 1;
    }

    std::unique_ptr<double[]> dca_bins;
    if (static_cast<int>(dcaBinning.value[2]) != 1.0) {
      dca_bins = calc_var_bins(N + 1, dcaBinning.value[1], static_cast<int>(dcaBinning.value[2]));
    } else {
      dca_bins = calc_const_bins(N, -dcaBinning.value[1], dcaBinning.value[1]);
    }
    auto const_bins_p = calc_const_bins(100, 0., 5.0);

    hSets.resize(stageDirs.size());
    hSetsTrue.resize(stageDirs.size());
    for (size_t i = 0; i < stageDirs.size(); i++) {
      for (size_t s = 0; s < signDirs.size(); s++) {
        const std::string dir = stageDirs[i] + "/" + signDirs[s];
        hSets[i][s] = makeHistSet(registry, dir, particlePDG.value);
        hSetsTrue[i][s] = makeHistSet(registryTrue, dir + "/true", particlePDG.value);
      }
    }

    const AxisSpec axisMult{5001, -0.5, 5000.5, "mult."};
    const AxisSpec axisPerc{101, -0.5, 100.5, "percentile"};
    registry.add("posZ", "posZ", kTH1F, {{300, -16., 16., "posZ"}});
    registry.add("mult", "mult", kTH1F, {axisMult});
    registry.add("MultVsCent", "MultVsCent", kTH2F, {axisMult, axisPerc});
    registry.add("IRvsOccupancy", "IRvsOccupancy", kTH2I, {{10000, 0, 10000, "Occupancy"}, {50, 0, 50, "IR, kHz"}});
  }

  template <bool FillExtra, bool HasMC, typename ColsType, typename TracksType>
  void fillHistograms(ColsType const& collisions, TracksType const& tracks)
  {
    for (const auto& track : tracks) {
      const int signIdx = (-track.sign() + 1) / 2;                                                   // 0: triton, 1: antitriton
      fillSet<FillExtra, HasMC>(hSets[1][signIdx], hSetsTrue[1][signIdx], track, particlePDG.value); // nocuts
    }
    for (const auto& collision : collisions) {
      if (removeSameBunchPileup && !collision.isNoSameBunchPileup())
        continue;
      if (requestGoodZvtxFT0vsPV && !collision.isGoodZvtxFT0vsPV())
        continue;
      if (requestVertexITSTPC && !collision.isVertexITSTPC())
        continue;
      if (requestVertexTOFmatched > collision.isVertexTOForTRDmatched())
        continue;
      if (requestNoCollInTimeRangeStandard && !collision.noCollInTimeRangeStandard())
        continue;
      if (requestIsGoodITSLayersAll && !collision.isGoodITSLayersAll())
        continue;
      if (collision.multPerc() < centCut.value.first || collision.multPerc() >= centCut.value.second)
        continue;
      if (collision.hadronicRate() < irCut.value.first || collision.hadronicRate() >= irCut.value.second)
        continue;
      if (collision.occupancy() < occupancyCut.value.first || collision.occupancy() >= occupancyCut.value.second)
        continue;

      registry.fill(HIST("posZ"), collision.posZ());
      registry.fill(HIST("mult"), collision.mult());
      registry.fill(HIST("MultVsCent"), collision.mult(), collision.multPerc());
      registry.fill(HIST("IRvsOccupancy"), collision.occupancy(), collision.hadronicRate());
    }

    for (const auto& track : tracks) {
      const auto& coll = track.template singleCollSel_as<ColsType>();
      const int signIdx = (-track.sign() + 1) / 2; // 0: triton, 1: antitriton

      if (removeSameBunchPileup && !coll.isNoSameBunchPileup())
        continue;
      if (requestGoodZvtxFT0vsPV && !coll.isGoodZvtxFT0vsPV())
        continue;
      if (requestVertexITSTPC && !coll.isVertexITSTPC())
        continue;
      if (requestVertexTOFmatched > coll.isVertexTOForTRDmatched())
        continue;
      if (requestNoCollInTimeRangeStandard && !coll.noCollInTimeRangeStandard())
        continue;
      if (requestIsGoodITSLayersAll && !coll.isGoodITSLayersAll())
        continue;
      if (std::fabs(coll.posZ()) > vertexZ)
        continue;
      if (coll.multPerc() < centCut.value.first || coll.multPerc() >= centCut.value.second)
        continue;
      if (coll.hadronicRate() < irCut.value.first || coll.hadronicRate() >= irCut.value.second)
        continue;
      if (coll.occupancy() < occupancyCut.value.first || coll.occupancy() >= occupancyCut.value.second)
        continue;
      if (track.tpcFractionSharedCls() > maxTpcFractionSharedCls)
        continue;
      if (track.itsNCls() < minItsNCls)
        continue;
      if (std::fabs(track.dcaXY()) > dcaXY.value[0] + dcaXY.value[1] * std::pow(track.pt(), dcaXY.value[2]))
        continue;
      if (std::fabs(track.dcaZ()) > dcaZ.value[0] + dcaZ.value[1] * std::pow(track.pt(), dcaZ.value[2]))
        continue;

      if (o2::aod::singletrackselector::ITSselection(track, ITScuts)) {
        fillSet<FillExtra, HasMC>(hSets[2][signIdx], hSetsTrue[2][signIdx], track, particlePDG.value); // ITSCuts
      }
      if (o2::aod::singletrackselector::TPCselection<false>(track, TPCcuts, itsNSigma.value)) {
        fillSet<FillExtra, HasMC>(hSets[3][signIdx], hSetsTrue[3][signIdx], track, particlePDG.value); // TPCCuts
      }
      if (o2::aod::singletrackselector::TOFselection(track, TOFcuts, tpcNSigmaResidual.value)) {
        fillSet<FillExtra, HasMC>(hSets[4][signIdx], hSetsTrue[4][signIdx], track, particlePDG.value); // TOFCuts
      }
      if (o2::aod::singletrackselector::TPCselection<true>(track, TPCcuts, itsNSigma.value)) {
        fillSet<FillExtra, HasMC>(hSets[5][signIdx], hSetsTrue[5][signIdx], track, particlePDG.value); // ITSTPC
      }

      const bool pidOk =
        (track.p() < pidTrshld)
          ? o2::aod::singletrackselector::TPCselection<true>(track, TPCcuts, itsNSigma.value)
          : o2::aod::singletrackselector::TPCselection<false>(track, std::make_pair(particlePDG.value, tpcNSigmaResidual.value), tpcNSigmaResidual.value);
      if (!pidOk)
        continue;
      if (std::fabs(static_cast<float>(RecoDecay::y(std::array{track.px(), track.py(), track.pz()}, kMassTriton))) > eta)
        continue;

      fillSet<FillExtra, HasMC>(hSets[6][signIdx], hSetsTrue[6][signIdx], track, particlePDG.value); // PIDcuts
      if constexpr (FillExtra) {
        if (track.p() >= massCuttrshld && (getMTOF(track) <= tofMass.value[0] || getMTOF(track) >= tofMass.value[1])) {
          continue;
        }
      }
      fillSet<FillExtra, HasMC>(hSets[7][signIdx], hSetsTrue[7][signIdx], track, particlePDG.value); // masscut

      bool rejectedByPDG = false;
      for (const int& pdgToReject : particlePDGtoReject.value) {
        if (o2::aod::singletrackselector::TPCselection<false>(track, std::make_pair(pdgToReject, rejectWithinNsigmaTPC.value))) {
          rejectedByPDG = true;
          break;
        }
      }
      if (rejectedByPDG)
        continue;
      fillSet<FillExtra, HasMC>(hSets[8][signIdx], hSetsTrue[8][signIdx], track, particlePDG.value); // rd

      if (getITSNsigma(track, particlePDG.value) <= itsNsigNL.value[0] + itsNsigNL.value[1] * track.p()) {
        continue;
      }
      fillSet<FillExtra, HasMC>(hSets[9][signIdx], hSetsTrue[9][signIdx], track, particlePDG.value); // nl
    }
  }

  void processDefault(soa::Filtered<soa::Join<aod::SingleCollSels, aod::SingleCollExtras>> const& collisions, soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes>> const& tracks)
  {
    fillHistograms<false, false>(collisions, tracks);
  }
  PROCESS_SWITCH(AntitritonAnalysis, processDefault, "process default", true);

  void processExtra(soa::Filtered<soa::Join<aod::SingleCollSels, aod::SingleCollExtras>> const& collisions, soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkExtras, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes>> const& tracks)
  {
    fillHistograms<true, false>(collisions, tracks);
  }
  PROCESS_SWITCH(AntitritonAnalysis, processExtra, "process extra", false);

  void processMC(soa::Filtered<soa::Join<aod::SingleCollSels, aod::SingleCollExtras>> const& collisions, soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkExtras, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes, aod::SingleTrkMCs>> const& tracks)
  {
    fillHistograms<true, true>(collisions, tracks);
  }
  PROCESS_SWITCH(AntitritonAnalysis, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntitritonAnalysis>(cfgc)};
}
