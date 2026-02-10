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

/// \file V0PhotonCut.h
/// \brief Header of class for V0 photon selection.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EmMlResponsePCM.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCandidate.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/Array2D.h>
#include <Framework/HistogramRegistry.h>

#include <TMath.h>
#include <TNamed.h>

#include <fairlogger/Logger.h>

#include <Rtypes.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace o2::analysis
{

// namespace per channel
namespace em_cuts_ml
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int NBins = 12;

static constexpr int NBinsPt = 12;
static constexpr int NCutScores = 2;
// default values for the pT bin edges, offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  0.25,
  0.5,
  0.75,
  1.,
  1.5,
  2.,
  4.,
  6.,
  10.,
  20.,
  50.,
  100.};
const auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};
static constexpr int NBinsCent = 11;
constexpr double BinsCent[NBinsCent + 1] = {
  0.,
  5,
  10,
  20,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100.};
const auto vecBinsCent = std::vector<double>{BinsCent, BinsCent + NBinsCent + 1};

// default values for the ML model paths, one model per pT bin
static const std::vector<std::string> modelPaths = {
  ""};

// default values for the cut directions
constexpr int CutDir[NCutScores] = {CutGreater, CutSmaller};
const auto vecCutDir = std::vector<int>{CutDir, CutDir + NCutScores};

// default values for the cuts
constexpr double Cuts[NBins][NCutScores] = {
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5}};

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};
// labels
static const std::vector<std::string> labelsCent = {
  "Cent bin 0",
  "Cent bin 1",
  "Cent bin 2",
  "Cent bin 3",
  "Cent bin 4",
  "Cent bin 5",
  "Cent bin 6",
  "Cent bin 7",
  "Cent bin 8",
  "Cent bin 9",
  "Cent bin 10"};

// column labels
static const std::vector<std::string> labelsCutScore = {"score primary photons", "score background"};
} // namespace em_cuts_ml

} // namespace o2::analysis

class V0PhotonCut : public TNamed
{
 public:
  V0PhotonCut() = default;
  V0PhotonCut(const char* name, const char* title) : TNamed(name, title) {}
  ~V0PhotonCut() override { delete mEmMlResponse; };

  enum class V0PhotonCuts : int {
    // v0 cut
    kMee = 0,
    kV0PtRange,
    kV0EtaRange,
    kAP,
    kPsiPair,
    kPhiV,
    kRxy,
    kCosPA,
    kPCA,
    kChi2KF,
    kRZLine,
    kOnWwireIB,
    kOnWwireOB,
    // leg cut
    kTrackPtRange,
    kTrackEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCFracSharedClusters,
    kTPCChi2NDF,
    kTPCNsigmaEl,
    kTPCNsigmaPi,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kITSClusterSize,
    kRequireITSTPC,
    kRequireITSonly,
    kRequireTPConly,
    kRequireTPCTRD,
    kRequireTPCTOF,
    kNCuts
  };

  /// \brief add histograms to registry
  /// \param fRegistry pointer to histogram registry
  void addQAHistograms(o2::framework::HistogramRegistry* fRegistry = nullptr) const
  {
    if (mDoQA && fRegistry != nullptr) {
      const o2::framework::AxisSpec thAxispT{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"};
      const o2::framework::AxisSpec thAxisMomentum{250, 0., 25., "#it{p} (GeV/#it{c})"};
      const o2::framework::AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
      const o2::framework::AxisSpec thAxisPhi{500, 0, o2::constants::math::TwoPI, "#varphi (rad)"};
      const o2::framework::AxisSpec thAxisNSigmaE{500, -10, 10, "#it{N#sigma}_{e} (a.u.)"};
      const o2::framework::AxisSpec thAxisNSigmaPi{500, -10, 10, "#it{N#sigma}_{#pi} (a.u.)"};
      const o2::framework::AxisSpec thAxisNClusTPC{152, 0., 152., "#it{N}_{cl., TPC}"};
      const o2::framework::AxisSpec thAxisNCrossedTPC{152, 0., 152., "#it{N}_{cr., TPC}"};
      const o2::framework::AxisSpec thAxisIDV0{1000, 0., 1000., "ID_{V0}"};
      const o2::framework::AxisSpec thAxisIDLeg{1000, 0., 1000., "ID_{Leg}"};
      const o2::framework::AxisSpec thAxisAlpha{250, -1., 1., "#alpha=(#it{p}^{+}_{L}-#it{p}^{-}_{L})/(#it{p}^{+}_{L}+#it{p}^{-}_{L})"};
      const o2::framework::AxisSpec thAxisQt{200, 0., 0.1, "#it{q}_{T} (GeV/#it{c})"};
      const o2::framework::AxisSpec thAxisConvX{320, -160., +160, "X (cm)"};
      const o2::framework::AxisSpec thAxisConvY{320, -160., +160, "Y (cm)"};
      const o2::framework::AxisSpec thAxisConvZ{180, -90., +90, "Z (cm)"};
      const o2::framework::AxisSpec thAxisConvR{200, 0., +100, "R (cm)"};
      const o2::framework::AxisSpec thAxisChi2{100, 0., +50, "#chi^{2}_{KF}/ndf"};
      const o2::framework::AxisSpec thAxisPsiPair{200, -0.1, +0.1, "#Psi_{pair}"};

      fRegistry->add("QA/V0Photon/before/hE", "p_{T};#it{p}_{T} (GeV/#it{c});#it{N}_{#gamma}", o2::framework::kTH1D, {thAxispT}, true);
      fRegistry->add("QA/V0Photon/before/hPt", "Transverse momenta of clusters;#it{p}_{T} (GeV/c);#it{N}_{#gamma}", o2::framework::kTH1D, {thAxispT}, true);
      fRegistry->add("QA/V0Photon/before/hNgamma", "Number of #gamma candidates per collision;#it{N}_{#gamma} per collision;#it{N}_{collisions}", o2::framework::kTH1D, {{1001, -0.5f, 1000.5f}}, true);
      fRegistry->add("QA/V0Photon/before/hEtaPhi", "#eta vs #varphi;#eta;#varphi (rad.)", o2::framework::kTH2F, {thAxisEta, thAxisPhi}, true);
      fRegistry->add("QA/V0Photon/before/hAP", "Armenteros-Podolanski #alpha vs qT", o2::framework::kTH2F, {thAxisAlpha, thAxisQt}, true);
      fRegistry->add("QA/V0Photon/before/hConvXY", "Conversion point XY", o2::framework::kTH2F, {thAxisConvX, thAxisConvY}, true);
      fRegistry->add("QA/V0Photon/before/hConvZR", "Conversion point ZR", o2::framework::kTH2F, {thAxisConvZ, thAxisConvR}, true);
      fRegistry->add("QA/V0Photon/before/hChi2", "Chi2/ndf from KFParticle;#chi^{2}_{KF}/ndf;counts", o2::framework::kTH1D, {thAxisChi2}, true);

      // TODO: add psi_pair once available
      // fRegistry->add("QA/V0Photon/before/hPsiPair", "Psi pair;#Psi_{pair};counts", o2::framework::kTH1D, {thAxisPsiPair}, true);

      fRegistry->add("QA/V0Photon/before/Pos/NSigmaE", "NSigmaE of pos leg vs momentum", o2::framework::kTH2F, {thAxisMomentum, thAxisNSigmaE}, true);
      fRegistry->add("QA/V0Photon/before/Pos/NSigmaPi", "NSigmaE of pos leg vs momentum", o2::framework::kTH2F, {thAxisMomentum, thAxisNSigmaPi}, true);
      fRegistry->add("QA/V0Photon/before/Pos/hEtaPhi", "eta vs phi of pos leg", o2::framework::kTH2F, {thAxisEta, thAxisPhi}, true);
      fRegistry->add("QA/V0Photon/before/Pos/hTPCHits", "NCluster vs NFindable TPC", o2::framework::kTH2F, {thAxisNClusTPC, thAxisNCrossedTPC}, true);
      fRegistry->add("QA/V0Photon/before/Neg/NSigmaE", "NSigmaE of neg leg vs momentum", o2::framework::kTH2F, {thAxisMomentum, thAxisNSigmaE}, true);
      fRegistry->add("QA/V0Photon/before/Neg/NSigmaPi", "NSigmaE of neg leg vs momentum", o2::framework::kTH2F, {thAxisMomentum, thAxisNSigmaPi}, true);
      fRegistry->add("QA/V0Photon/before/Neg/hEtaPhi", "eta vs phi of neg leg", o2::framework::kTH2F, {thAxisEta, thAxisPhi}, true);
      fRegistry->add("QA/V0Photon/before/Neg/hTPCHits", "NCluster vs NFindable TPC", o2::framework::kTH2F, {thAxisNClusTPC, thAxisNCrossedTPC}, true);

      fRegistry->addClone("QA/V0Photon/before/", "QA/V0Photon/after/");

      auto hPhotonQualityCuts = fRegistry->add<TH2>("QA/V0Photon/hPhotonQualityCuts", "pT at which v0 photons are removed by a given cut", o2::framework::kTH2F, {{static_cast<int>(V0PhotonCut::V0PhotonCuts::kNCuts) + 2, -0.5, static_cast<double>(V0PhotonCut::V0PhotonCuts::kNCuts) + 1.5}, thAxispT}, true);
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(1, "In");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(2, "#it{M}_{ee}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(3, "#it{p}_{T}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(4, "#it{#eta}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(5, "AP");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(6, "#Psi_{pair}"); // currently not implemented!
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(7, "#Phi_{v}");    // currently not implemented!
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(8, "#it{R}_{xy}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(9, "CosPA");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(10, "PCA");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(11, "#chi^{2}_{KF}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(12, "RZ_{line}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(13, "Wire_{IB}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(14, "Wire_{OB}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(15, "#it{p}_{T,leg}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(16, "#it{#eta}_{leg}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(17, "#it{N}_{cl,TPC}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(18, "#it{N}_{cr,TPC}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(19, "#it{N}_{cr,TPC}/#it{N}_{cl,TPC}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(20, "FracSharedCl");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(21, "#chi^{2}_{TPC}/NDF");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(22, "#it{N#sigma}_{e,TPC}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(23, "#it{N#sigma}_{#pi,TPC}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(24, "DCA_{xy}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(25, "DCA_{z}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(26, "#it{N}_{cl,ITS}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(27, "#chi^{2}_{ITS}/NDF");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(28, "size_{ITS}");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(29, "ITSTPC");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(30, "ITSOnly");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(31, "TPCOnly");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(32, "TPCTRD");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(33, "TPCTOF");
      hPhotonQualityCuts->GetXaxis()->SetBinLabel(34, "Out");
    }
  }

  template <o2::soa::is_iterator TV0, o2::soa::is_iterator TLeg1, o2::soa::is_iterator TLeg2>
  void fillBeforePhotonHistogram(TV0 const& v0, TLeg1 const& pos, TLeg2 const& ele, o2::framework::HistogramRegistry* fRegistry = nullptr) const
  {

    if (mDoQA == false || fRegistry == nullptr) {
      return;
    }

    fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), 0, v0.pt());
    fRegistry->fill(HIST("QA/V0Photon/before/hE"), v0.e());
    fRegistry->fill(HIST("QA/V0Photon/before/hPt"), v0.pt());
    fRegistry->fill(HIST("QA/V0Photon/before/hEtaPhi"), v0.eta(), v0.phi());
    fRegistry->fill(HIST("QA/V0Photon/before/hAP"), v0.alpha(), v0.qtarm());
    fRegistry->fill(HIST("QA/V0Photon/before/hConvXY"), v0.vx(), v0.vy());
    fRegistry->fill(HIST("QA/V0Photon/before/hConvZR"), v0.vz(), v0.v0radius());
    fRegistry->fill(HIST("QA/V0Photon/before/hChi2"), v0.chiSquareNDF());

    // TODO: add psi_pair once available
    // fRegistry->fill(HIST("QA/V0Photon/before/hPsiPair"), v0.psiPair());

    fRegistry->fill(HIST("QA/V0Photon/before/Pos/NSigmaE"), pos.p(), pos.tpcNSigmaEl());
    fRegistry->fill(HIST("QA/V0Photon/before/Pos/NSigmaPi"), pos.p(), pos.tpcNSigmaPi());
    fRegistry->fill(HIST("QA/V0Photon/before/Pos/hEtaPhi"), pos.eta(), pos.phi());
    fRegistry->fill(HIST("QA/V0Photon/before/Pos/hTPCHits"), pos.tpcNClsFound(), pos.tpcNClsCrossedRows());
    fRegistry->fill(HIST("QA/V0Photon/before/Neg/NSigmaE"), ele.p(), ele.tpcNSigmaEl());
    fRegistry->fill(HIST("QA/V0Photon/before/Neg/NSigmaPi"), ele.p(), ele.tpcNSigmaPi());
    fRegistry->fill(HIST("QA/V0Photon/before/Neg/hEtaPhi"), ele.eta(), ele.phi());
    fRegistry->fill(HIST("QA/V0Photon/before/Neg/hTPCHits"), ele.tpcNClsFound(), ele.tpcNClsCrossedRows());
  }

  template <o2::soa::is_iterator TV0, o2::soa::is_iterator TLeg1, o2::soa::is_iterator TLeg2>
  void fillAfterPhotonHistogram(TV0 const& v0, TLeg1 const& pos, TLeg2 const& ele, o2::framework::HistogramRegistry* fRegistry = nullptr) const
  {

    if (mDoQA == false || fRegistry == nullptr) {
      return;
    }

    fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kNCuts) + 1, v0.pt());
    fRegistry->fill(HIST("QA/V0Photon/after/hE"), v0.e());
    fRegistry->fill(HIST("QA/V0Photon/after/hPt"), v0.pt());
    fRegistry->fill(HIST("QA/V0Photon/after/hEtaPhi"), v0.eta(), v0.phi());
    fRegistry->fill(HIST("QA/V0Photon/after/hAP"), v0.alpha(), v0.qtarm());
    fRegistry->fill(HIST("QA/V0Photon/after/hConvXY"), v0.vx(), v0.vy());
    fRegistry->fill(HIST("QA/V0Photon/after/hConvZR"), v0.vz(), v0.v0radius());
    fRegistry->fill(HIST("QA/V0Photon/after/hChi2"), v0.chiSquareNDF());

    // TODO: add psi_pair once available
    // fRegistry->fill(HIST("QA/V0Photon/after/hPsiPair"), v0.psiPair());

    fRegistry->fill(HIST("QA/V0Photon/after/Pos/NSigmaE"), pos.p(), pos.tpcNSigmaEl());
    fRegistry->fill(HIST("QA/V0Photon/after/Pos/NSigmaPi"), pos.p(), pos.tpcNSigmaPi());
    fRegistry->fill(HIST("QA/V0Photon/after/Pos/hEtaPhi"), pos.eta(), pos.phi());
    fRegistry->fill(HIST("QA/V0Photon/after/Pos/hTPCHits"), pos.tpcNClsFound(), pos.tpcNClsCrossedRows());
    fRegistry->fill(HIST("QA/V0Photon/after/Neg/NSigmaE"), ele.p(), ele.tpcNSigmaEl());
    fRegistry->fill(HIST("QA/V0Photon/after/Neg/NSigmaPi"), ele.p(), ele.tpcNSigmaPi());
    fRegistry->fill(HIST("QA/V0Photon/after/Neg/hEtaPhi"), ele.eta(), ele.phi());
    fRegistry->fill(HIST("QA/V0Photon/after/Neg/hTPCHits"), ele.tpcNClsFound(), ele.tpcNClsCrossedRows());
  }

  /// \brief check if given v0 photon survives all cuts
  /// \param flags EMBitFlags where results will be stored
  /// \param v0s v0 photon table to check
  template <o2::soa::is_table TV0, typename TLeg>
  void AreSelectedRunning(EMBitFlags& flags, TV0 const& v0s, o2::framework::HistogramRegistry* fRegistry = nullptr) const
  {
    // auto legIter = legs.begin();
    // auto legEnd = legs.end();
    size_t iV0 = 0;

    const bool doQA = mDoQA && fRegistry != nullptr;

    uint nTotV0PerColl = 0;
    currentCollID = v0s.iteratorAt(0).emeventId();

    for (const auto& v0 : v0s) {
      const auto collID = v0.emeventId();
      if (!IsSelected<decltype(v0), TLeg>(v0, fRegistry)) {
        flags.set(iV0);
      }
      if (collID != currentCollID) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/before/hNgamma"), nTotV0PerColl);
        }
        nTotV0PerColl = 0;
        currentCollID = collID;
      }
      ++nTotV0PerColl;
      ++iV0;
    }
  }

  template <o2::soa::is_iterator TV0, typename TLeg>
  bool IsSelected(TV0 const& v0, o2::framework::HistogramRegistry* fRegistry = nullptr) const
  {
    auto pos = v0.template posTrack_as<TLeg>();
    auto ele = v0.template negTrack_as<TLeg>();

    const float v0Pt = v0.pt();

    const auto doQA = mDoQA && fRegistry != nullptr;

    if (doQA) {
      fillBeforePhotonHistogram(v0, pos, ele, fRegistry);
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kMee)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kMee) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kV0PtRange)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kV0PtRange) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kV0EtaRange)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kV0EtaRange) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kAP)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kAP) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPsiPair)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kPsiPair) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPhiV)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kPhiV) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kRxy)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRxy) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kCosPA)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kCosPA) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPCA)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kPCA) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kChi2KF)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kChi2KF) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kRZLine)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRZLine) + 1, v0Pt);
      }
      return false;
    }

    if (mIsOnWwireIB && mIsOnWwireOB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB) && !IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kOnWwireIB) + 1, v0Pt);
        }
        return false;
      }
    } else if (mIsOnWwireIB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kOnWwireIB) + 1, v0Pt);
        }
        return false;
      }
    } else if (mIsOnWwireOB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kOnWwireOB) + 1, v0Pt);
        }
        return false;
      }
    }

    for (const auto& track : {pos, ele}) {
      if (!IsSelectedTrack(track, V0PhotonCuts::kTrackPtRange)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTrackPtRange) + 1, v0Pt);
        }
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kTrackEtaRange)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTrackEtaRange) + 1, v0Pt);
        }
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAxy)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kDCAxy) + 1, v0Pt);
        }
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAz)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kDCAz) + 1, v0Pt);
        }
        return false;
      }
      if (!track.hasITS() && !track.hasTPC()) { // track has to be ITSonly or TPConly or ITS-TPC
        return false;
      }
      if (mDisableITSonly && o2::pwgem::photonmeson::isITSonlyTrack(track)) {
        return false;
      }
      if (mDisableTPConly && o2::pwgem::photonmeson::isTPConlyTrack(track)) {
        return false;
      }

      if (mRejectITSib) {
        auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        auto hits_ob = std::count_if(its_ob_Requirement.second.begin(), its_ob_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only = (hits_ib <= its_ib_Requirement.first) && (hits_ob >= its_ob_Requirement.first);
        if (o2::pwgem::photonmeson::isITSonlyTrack(track) && !its_ob_only) { // ITSonly tracks should not have any ITSib hits.
          return false;
        }

        auto hits_ob_itstpc = std::count_if(its_ob_Requirement_ITSTPC.second.begin(), its_ob_Requirement_ITSTPC.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only_itstpc = (hits_ib <= its_ib_Requirement.first) && (hits_ob_itstpc >= its_ob_Requirement_ITSTPC.first);
        if (o2::pwgem::photonmeson::isITSTPCTrack(track) && !its_ob_only_itstpc) { // ITSTPC tracks should not have any ITSib hits.
          return false;
        }
      }

      if (track.hasITS() && !CheckITSCuts(track, fRegistry, v0Pt)) {
        return false;
      }
      if (track.hasTPC() && !CheckTPCCuts(track, fRegistry, v0Pt)) {
        return false;
      }
      if (track.hasITS() && !track.hasTPC() && (track.hasTRD() || track.hasTOF())) { // remove ITS-TRD, ITS-TOF, ITS-TRD-TOF that are unrealistic tracks.
        return false;
      }

      if (mIsOnWwireIB && !CheckITSCuts(track)) { // photon conversion on ibw requires ITS hits.
        return false;
      }

      if (mRequireITSonly && !IsSelectedTrack(track, V0PhotonCuts::kRequireITSonly)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRequireITSonly) + 1, v0Pt);
        }
        return false;
      }
      if (mRequireITSTPC && !IsSelectedTrack(track, V0PhotonCuts::kRequireITSTPC)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRequireITSTPC) + 1, v0Pt);
        }
        return false;
      }
      if (mRequireTPConly && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPConly)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRequireTPConly) + 1, v0Pt);
        }
        return false;
      }
      if (mRequireTPCTRD && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTRD)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRequireTPCTRD) + 1, v0Pt);
        }
        return false;
      }
      if (mRequireTPCTOF && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTOF)) {
        if (doQA) {
          fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kRequireTPCTOF) + 1, v0Pt);
        }
        return false;
      }
    }
    if (mApplyMlCuts) {
      if (!mEmMlResponse) {
        LOG(error) << "EM ML Response is not initialized!";
        return false;
      }
      bool mIsSelectedMl = false;
      std::vector<float> mOutputML;
      V0PhotonCandidate v0photoncandidate(v0, pos, ele, mCentFT0A, mCentFT0C, mCentFT0M, mD_Bz);
      std::vector<float> mlInputFeatures = mEmMlResponse->getInputFeatures(v0photoncandidate, pos, ele);
      if (mUse2DBinning) {
        if (mCentralityTypeMl == "CentFT0C") {
          mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCentFT0C(), mOutputML);
        } else if (mCentralityTypeMl == "CentFT0A") {
          mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCentFT0A(), mOutputML);
        } else if (mCentralityTypeMl == "CentFT0M") {
          mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCentFT0M(), mOutputML);
        } else {
          LOG(fatal) << "Unsupported centTypePCMMl: " << mCentralityTypeMl << " , please choose from CentFT0C, CentFT0A, CentFT0M.";
        }
      } else {
        mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), mOutputML);
      }
      if (!mIsSelectedMl) {
        return false;
      }
    }
    if (doQA) {
      fillAfterPhotonHistogram(v0, pos, ele, fRegistry);
      if (v0.emeventId() != currentCollID) {
        fRegistry->fill(HIST("QA/V0Photon/after/hNgamma"), nAccV0PerColl);
        nAccV0PerColl = 0;
      }
      ++nAccV0PerColl;
    }
    return true;
  }

  template <typename T>
  bool CheckITSCuts(T const& track, o2::framework::HistogramRegistry* fRegistry = nullptr, const float v0Pt = 0.f) const
  {
    const auto doQA = mDoQA && fRegistry != nullptr;
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSNCls)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kITSNCls) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSChi2NDF)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kITSChi2NDF) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSClusterSize)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kITSClusterSize) + 1, v0Pt);
      }
      return false;
    }
    return true;
  }

  template <typename T>
  bool CheckTPCCuts(T const& track, o2::framework::HistogramRegistry* fRegistry = nullptr, const float v0Pt = 0.f) const
  {
    const auto doQA = mDoQA && fRegistry != nullptr;
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNCls)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCNCls) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCCrossedRows)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCCrossedRows) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCCrossedRowsOverNCls)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCCrossedRowsOverNCls) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCFracSharedClusters)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCFracSharedClusters) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCChi2NDF)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCChi2NDF) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNsigmaEl)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCNsigmaEl) + 1, v0Pt);
      }
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNsigmaPi)) {
      if (doQA) {
        fRegistry->fill(HIST("QA/V0Photon/hPhotonQualityCuts"), static_cast<int>(V0PhotonCuts::kTPCNsigmaPi) + 1, v0Pt);
      }
      return false;
    }
    return true;
  }

  template <o2::soa::is_iterator T>
  bool IsSelectedV0(T const& v0, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kV0PtRange:
        return v0.pt() >= mMinV0Pt && v0.pt() <= mMaxV0Pt;

      case V0PhotonCuts::kV0EtaRange:
        return v0.eta() >= mMinV0Eta && v0.eta() <= mMaxV0Eta;

      case V0PhotonCuts::kMee:
        return v0.mGamma() < mMaxQt * 2.f;

      case V0PhotonCuts::kAP:
        return std::pow(v0.alpha() / mMaxAlpha, 2) + std::pow(v0.qtarm() / mMaxQt, 2) < 1.0;

      case V0PhotonCuts::kPsiPair:
        return true;

      case V0PhotonCuts::kPhiV:
        return true;

      case V0PhotonCuts::kRxy: {
        if (v0.v0radius() < mMinRxy || mMaxRxy < v0.v0radius()) {
          return false;
        }
        return true;
      }

      case V0PhotonCuts::kCosPA:
        return v0.cospa() >= mMinCosPA;

      case V0PhotonCuts::kPCA:
        return v0.pca() <= mMaxPCA;

      case V0PhotonCuts::kChi2KF:
        return v0.chiSquareNDF() <= mMaxChi2KF;

      case V0PhotonCuts::kRZLine:
        return v0.v0radius() > std::fabs(v0.vz()) * std::tan(2 * std::atan(std::exp(-mMaxV0Eta))) - mMaxMarginZ;

      case V0PhotonCuts::kOnWwireIB: {
        const float margin_xy = 1.0; // cm
        // const float margin_z = 20.0; // cm
        // const float rxy_min = 5.506; // cm
        // const float rxy_max = 14.846;         // cm
        // const float z_min = -17.56; // cm
        // const float z_max = +31.15;           // cm
        float x = std::fabs(v0.vx()); // cm, measured secondary vertex of gamma->ee
        float y = v0.vy();            // cm, measured secondary vertex of gamma->ee
        float z = v0.vz();            // cm, measured secondary vertex of gamma->ee

        float rxy = std::sqrt(x * x + y * y);
        if (rxy < 7.0 || 14.0 < rxy) {
          return false;
        }

        // r = 0.192 * z + 8.88 (cm) expected wire position in RZ plane.TMath::Tan(10.86 * o2::constants::math::Deg2Rad) = 0.192
        if (rxy > 0.192 * z + 14.0) { // upper limit
          return false;
        }

        float dxy = std::fabs(1.0 * y - x * std::tan(-8.52f * o2::constants::math::Deg2Rad)) / std::sqrt(std::pow(1.0, 2) + std::pow(std::tan(-8.52f * o2::constants::math::Deg2Rad), 2));
        return !(dxy > margin_xy);
      }
      case V0PhotonCuts::kOnWwireOB: {
        const float margin_xy = 1.0;                                                  // cm
        const float rxy_exp = 30.8;                                                   // cm
        const float x_exp = rxy_exp * std::cos(-1.3f * o2::constants::math::Deg2Rad); // cm, expected position x of W wire
        const float y_exp = rxy_exp * std::sin(-1.3f * o2::constants::math::Deg2Rad); // cm, expected position y of W wire
        // const float z_min = -47.0;                                          // cm
        // const float z_max = +47.0;                                          // cm
        float x = v0.vx(); // cm, measured secondary vertex of gamma->ee
        float y = v0.vy(); // cm, measured secondary vertex of gamma->ee
        // float z = v0.vz();                                                  // cm, measured secondary vertex of gamma->ee

        // float rxy = sqrt(x * x + y * y);
        // if (rxy < 28.0 || 33.0 < rxy) {
        //   return false;
        // }

        float dxy = std::sqrt(std::pow(x - x_exp, 2) + std::pow(y - y_exp, 2));
        return !(dxy > margin_xy);
      }
      default:
        return false;
    }
  }

  // Temporary function to check if track passes a given selection criteria. To be replaced by framework filters.
  template <typename T>
  bool IsSelectedTrack(T const& track, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kTrackPtRange:
        return track.pt() > mMinTrackPt && track.pt() < mMaxTrackPt;

      case V0PhotonCuts::kTrackEtaRange:
        return track.eta() > mMinTrackEta && track.eta() < mMaxTrackEta;

      case V0PhotonCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case V0PhotonCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case V0PhotonCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case V0PhotonCuts::kTPCFracSharedClusters:
        return track.tpcFractionSharedCls() < mMaxFracSharedClustersTPC;

      case V0PhotonCuts::kTPCChi2NDF:
        return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

      case V0PhotonCuts::kTPCNsigmaEl:
        return track.tpcNSigmaEl() > mMinTPCNsigmaEl && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;

      case V0PhotonCuts::kTPCNsigmaPi:
        return track.tpcNSigmaPi() > mMinTPCNsigmaPi && track.tpcNSigmaPi() < mMaxTPCNsigmaPi;

      case V0PhotonCuts::kDCAxy:
        return std::fabs(track.dcaXY()) < ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case V0PhotonCuts::kDCAz:
        return std::fabs(track.dcaZ()) < mMaxDcaZ;

      case V0PhotonCuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case V0PhotonCuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      case V0PhotonCuts::kITSClusterSize: {
        if (!o2::pwgem::photonmeson::isITSonlyTrack(track)) {
          return true;
        }
        return mMinMeanClusterSizeITS < track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())) && track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())) < mMaxMeanClusterSizeITS;
      }

      case V0PhotonCuts::kRequireITSTPC:
        return o2::pwgem::photonmeson::isITSTPCTrack(track);

      case V0PhotonCuts::kRequireITSonly:
        return o2::pwgem::photonmeson::isITSonlyTrack(track);

      case V0PhotonCuts::kRequireTPConly:
        return o2::pwgem::photonmeson::isTPConlyTrack(track);

      case V0PhotonCuts::kRequireTPCTRD:
        return o2::pwgem::photonmeson::isTPCTRDTrack(track);

      case V0PhotonCuts::kRequireTPCTOF:
        return o2::pwgem::photonmeson::isTPCTOFTrack(track);

      default:
        return false;
    }
  }

  void initV0MlModels(o2::ccdb::CcdbApi& ccdbApi)
  {
    if (!mEmMlResponse) {
      mEmMlResponse = new o2::analysis::EmMlResponsePCM<float>();
    }
    if (mUse2DBinning) {
      int binsNPt = static_cast<int>(mBinsPtMl.size()) - 1;
      int binsNCent = static_cast<int>(mBinsCentMl.size()) - 1;
      int binsN = binsNPt * binsNCent;
      if (binsN * static_cast<int>(mCutDirMl.size()) != static_cast<int>(mCutsMlFlat.size())) {
        LOG(fatal) << "Mismatch in number of bins and cuts provided for 2D ML application: binsN * mCutDirMl: " << int(binsN) * int(mCutDirMl.size()) << " bins vs. mCutsMlFlat: " << mCutsMlFlat.size() << " cuts";
      }
      if (binsN != static_cast<int>(mOnnxFileNames.size())) {
        LOG(fatal) << "Mismatch in number of bins and ONNX files provided for 2D ML application: binsN " << binsN << " bins vs. mOnnxFileNames: " << mOnnxFileNames.size() << " ONNX files";
      }
      if (binsN != static_cast<int>(mLabelsBinsMl.size())) {
        LOG(fatal) << "Mismatch in number of bins and labels provided for 2D ML application: binsN:" << binsN << " bins vs. mLabelsBinsMl: " << mLabelsBinsMl.size() << " labels";
      }
      if (static_cast<int>(mCutDirMl.size()) != mNClassesMl) {
        LOG(fatal) << "Mismatch in number of classes and cut directions provided for 2D ML application: mNClassesMl: " << mNClassesMl << " classes vs. mCutDirMl: " << mCutDirMl.size() << " cut directions";
      }
      if (static_cast<int>(mLabelsCutScoresMl.size()) != mNClassesMl) {
        LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for 2D ML application: mNClassesMl: " << mNClassesMl << " classes vs. mLabelsCutScoresMl: " << mLabelsCutScoresMl.size() << " labels";
      }
      o2::framework::LabeledArray<double> mCutsMl(mCutsMlFlat.data(), binsN, mNClassesMl, mLabelsBinsMl, mLabelsCutScoresMl);
      mEmMlResponse->configure2D(mBinsPtMl, mBinsCentMl, mCutsMl, mCutDirMl, mNClassesMl);
    } else {
      int binsNPt = static_cast<int>(mBinsPtMl.size()) - 1;
      if (binsNPt * static_cast<int>(mCutDirMl.size()) != static_cast<int>(mCutsMlFlat.size())) {
        LOG(fatal) << "Mismatch in number of pT bins and cuts provided for ML application: binsNPt * mCutDirMl:" << binsNPt * mCutDirMl.size() << " bins vs. mCutsMlFlat: " << mCutsMlFlat.size() << " cuts";
      }
      if (binsNPt != static_cast<int>(mOnnxFileNames.size())) {
        LOG(fatal) << "Mismatch in number of pT bins and ONNX files provided for ML application: binsNPt " << binsNPt << " bins vs. mOnnxFileNames: " << mOnnxFileNames.size() << " ONNX files";
      }
      if (binsNPt != static_cast<int>(mLabelsBinsMl.size())) {
        LOG(fatal) << "Mismatch in number of pT bins and labels provided for ML application: binsNPt:" << binsNPt << " bins vs. mLabelsBinsMl: " << mLabelsBinsMl.size() << " labels";
      }
      if (mNClassesMl != static_cast<int>(mCutDirMl.size())) {
        LOG(fatal) << "Mismatch in number of classes and cut directions provided for ML application: mNClassesMl: " << mNClassesMl << " classes vs. mCutDirMl: " << mCutDirMl.size() << " cut directions";
      }
      if (static_cast<int>(mLabelsCutScoresMl.size()) != mNClassesMl) {
        LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for ML application: mNClassesMl:" << mNClassesMl << " classes vs. mLabelsCutScoresMl: " << mLabelsCutScoresMl.size() << " labels";
      }
      o2::framework::LabeledArray<double> mCutsMl(mCutsMlFlat.data(), binsNPt, mNClassesMl, mLabelsBinsMl, mLabelsCutScoresMl);
      mEmMlResponse->configure(mBinsPtMl, mCutsMl, mCutDirMl, mNClassesMl);
    }
    if (mLoadMlModelsFromCCDB) {
      ccdbApi.init(mCcdbUrl);
      mEmMlResponse->setModelPathsCCDB(mOnnxFileNames, ccdbApi, mModelPathsCCDB, mTimestampCCDB);
    } else {
      mEmMlResponse->setModelPathsLocal(mOnnxFileNames);
    }
    mEmMlResponse->cacheInputFeaturesIndices(mNamesInputFeatures);
    mEmMlResponse->init();
  }

  template <o2::soa::is_iterator TMCPhoton>
  bool IsConversionPointInAcceptance(TMCPhoton const& mcphoton) const
  {

    float rGenXY = std::sqrt(std::pow(mcphoton.vx(), 2) + std::pow(mcphoton.vy(), 2));

    // eta cut
    if (mcphoton.eta() >= mMinV0Eta && mcphoton.eta() <= mMaxV0Eta) {
      return false;
    }

    // radius cut
    if (rGenXY < mMinRxy || mMaxRxy < rGenXY) {
      return false;
    }

    // line cut
    if (rGenXY < std::abs(mcphoton.vz()) * std::tan(2 * std::atan(std::exp(-mMaxV0Eta))) - mMaxMarginZ) {
      return false;
    }

    return true;
  }

  // Setters
  void SetV0PtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetV0EtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMeeRange(float min = 0.f, float max = 0.1);
  void SetPsiPairRange(float min = -3.15, float max = +3.15);
  void SetPhivPairRange(float min = 0.f, float max = +3.15);
  void SetAPRange(float max_alpha = 0.95, float max_qt = 0.05); // Armenteros Podolanski
  void SetRxyRange(float min = 0.f, float max = 180.f);
  void SetMinCosPA(float min = 0.95);
  void SetMaxPCA(float max = 2.f);
  void SetMaxChi2KF(float max = 1e+10);
  void SetMaxMarginZ(float max = 7.f);
  void SetMaxMeePsiPairDep(std::function<float(float)> psiDepCut);
  void SetOnWwireIB(bool flag = false);
  void SetOnWwireOB(bool flag = false);
  void RejectITSib(bool flag = false);

  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetMaxFracSharedClustersTPC(float max);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);
  void SetMeanClusterSizeITSob(float min, float max);

  void SetTPCNsigmaElRange(float min = -3, float max = +3);
  void SetTPCNsigmaPiRange(float min = -1e+10, float max = 1e+10);

  void SetMaxDcaXY(float maxDcaXY);
  void SetMaxDcaZ(float maxDcaZ);
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void SetIsWithinBeamPipe(bool flag);
  void SetRequireITSTPC(bool flag);
  void SetRequireITSonly(bool flag);
  void SetRequireTPConly(bool flag);
  void SetRequireTPCTRD(bool flag);
  void SetRequireTPCTOF(bool flag);
  void SetDisableITSonly(bool flag);
  void SetDisableTPConly(bool flag);

  void SetApplyMlCuts(bool flag = false);
  void SetUse2DBinning(bool flag = true);
  void SetLoadMlModelsFromCCDB(bool flag = true);
  void SetNClassesMl(int nClasses);
  void SetMlTimestampCCDB(int timestamp);
  void SetCentrality(float centFT0A, float centFT0C, float centFT0M);
  void SetD_Bz(float d_bz);
  void SetCcdbUrl(const std::string& url = "http://alice-ccdb.cern.ch");
  void SetCentralityTypeMl(const std::string& centType);
  void SetCutDirMl(const std::vector<int>& cutDirMl);
  void SetMlModelPathsCCDB(const std::vector<std::string>& modelPaths);
  void SetMlOnnxFileNames(const std::vector<std::string>& onnxFileNamesVec);
  void SetLabelsBinsMl(const std::vector<std::string>& labelsBins);
  void SetLabelsCutScoresMl(const std::vector<std::string>& labelsCutScores);
  void SetBinsPtMl(const std::vector<double>& binsPt);
  void SetBinsCentMl(const std::vector<double>& binsCent);
  void SetCutsMl(const std::vector<double>& cutsMlFlat);
  void SetNamesInputFeatures(const std::vector<std::string>& namesInputFeaturesVec);

  void setDoQA(bool flag = false);

 private:
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ob_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ob_Requirement_ITSTPC;
  // v0 cuts
  float mMinMee{0.f}, mMaxMee{0.1f};
  float mMinV0Pt{0.f}, mMaxV0Pt{1e10f};      // range in pT
  float mMinV0Eta{-1e10f}, mMaxV0Eta{1e10f}; // range in eta
  float mMaxAlpha{0.95}, mMaxQt{0.05};
  float mMinPsiPair{-3.15}, mMaxPsiPair{+3.15};
  float mMinPhivPair{0.f}, mMaxPhivPair{+3.15};
  float mMinRxy{0.f}, mMaxRxy{180.f};
  float mMinCosPA{0.95};
  float mMaxPCA{2.f};
  float mMaxChi2KF{1e+10};
  float mMaxMarginZ{7.f};
  std::function<float(float)> mMaxMeePsiPairDep{}; // max mee as a function of psipair
  bool mIsOnWwireIB{false};
  bool mIsOnWwireOB{false};
  bool mRejectITSib{false};

  // ML cuts
  bool mApplyMlCuts{false};
  bool mUse2DBinning{true};
  bool mLoadMlModelsFromCCDB{true};
  int mTimestampCCDB{-1};
  int mNClassesMl{static_cast<int>(o2::analysis::em_cuts_ml::NCutScores)};
  float mCentFT0A{0.f};
  float mCentFT0C{0.f};
  float mCentFT0M{0.f};
  float mD_Bz{0.f};
  std::string mCcdbUrl{"http://alice-ccdb.cern.ch"};
  std::string mCentralityTypeMl{"CentFT0C"};
  std::vector<int> mCutDirMl{std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}};
  std::vector<std::string> mModelPathsCCDB{std::vector<std::string>{"path_ccdb/BDT_PCM/"}};
  std::vector<std::string> mOnnxFileNames{std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}};
  std::vector<std::string> mNamesInputFeatures{std::vector<std::string>{"feature1", "feature2"}};
  std::vector<std::string> mLabelsBinsMl{std::vector<std::string>{"bin 0", "bin 1"}};
  std::vector<std::string> mLabelsCutScoresMl{std::vector<std::string>{"score primary photons", "score background"}};
  std::vector<double> mBinsPtMl{std::vector<double>{o2::analysis::em_cuts_ml::vecBinsPt}};
  std::vector<double> mBinsCentMl{std::vector<double>{o2::analysis::em_cuts_ml::vecBinsCent}};
  std::vector<double> mCutsMlFlat{std::vector<double>{0.5}};
  o2::analysis::EmMlResponsePCM<float>* mEmMlResponse{nullptr};

  // pid cuts
  float mMinTPCNsigmaEl{-5}, mMaxTPCNsigmaEl{+5};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta

  // track quality cuts
  int mMinNClustersTPC{0};                                             // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                                          // min number of crossed rows in TPC
  float mMinChi2PerClusterTPC{-1e10f}, mMaxChi2PerClusterTPC{1e10f};   // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f};                  // min ratio crossed rows / findable clusters
  float mMaxFracSharedClustersTPC{999.f};                              // max ratio shared clusters / clusters in TPC
  int mMinNClustersITS{0}, mMaxNClustersITS{7};                        // range in number of ITS clusters
  float mMinChi2PerClusterITS{-1e10f}, mMaxChi2PerClusterITS{1e10f};   // max its fit chi2 per ITS cluster
  float mMinMeanClusterSizeITS{-1e10f}, mMaxMeanClusterSizeITS{1e10f}; // max <its cluster size> x cos(Lmabda)

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mRequireITSTPC{false};
  bool mRequireITSonly{false};
  bool mRequireTPConly{false};
  bool mRequireTPCTRD{false};
  bool mRequireTPCTOF{false};
  bool mDisableITSonly{false};
  bool mDisableTPConly{false};

  bool mDoQA{false};             ///< flag to decide if QA should be done or not
  mutable uint nAccV0PerColl{0}; ///< running number of accepted v0 photons per collision used for QA
  mutable int currentCollID{-1}; ///< running collision ID of v0 photon used for QA

  ClassDef(V0PhotonCut, 5);
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
