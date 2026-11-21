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

/// \file flowGfwDerivedProducer.cxx
/// \brief Produce the GFW derived dataset used by flowGfwV02.
///
/// Keeps sel8 collisions and global tracks that pass the ITS cluster and ITS chi2
/// cuts. Stores the event mask and the combined track and PID bitmask. PID bits
/// stay unset above the TOF pT cut when the track has no TOF match. Run it after
/// event selection, track selection, multiplicity, centrality, and the TPC/TOF PID
/// helpers (pid-tpc, pid-tof-full, pid-tof-beta).
///
/// \author Maxim Virta, NBI, maxim.virta@cern.ch

#include "PWGCF/GenericFramework/DataModel/GfwDerived.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr float VertexResMax = 0.25f;
constexpr int VertexContribMin = 20;
constexpr float DcaShapeMinPt = 0.001f;
constexpr float DcaShapeMaxPt = 100.f;
constexpr float DcaZReferenceNsigma = 7.f;
constexpr float NsigmaScaleNominal = 1.f;
constexpr float NsigmaScaleTight = 2.f / 3.f;
constexpr float NsigmaScaleLoose = 4.f / 3.f;
constexpr uint32_t NsigmaRows = 6;
constexpr uint32_t NsigmaCols = 3;
constexpr int NsigmaUpPi = 0;
constexpr int NsigmaUpKa = 1;
constexpr int NsigmaUpPr = 2;
constexpr int NsigmaLowPi = 3;
constexpr int NsigmaLowKa = 4;
constexpr int NsigmaLowPr = 5;
constexpr int NsigmaColTpc = 0;
constexpr int NsigmaColTof = 1;
constexpr std::array<std::array<float, NsigmaCols>, NsigmaRows> NsigmaDefault = {{{{1.1f, 2.1f, 3.1f}}, {{1.2f, 2.2f, 3.2f}}, {{1.3f, 2.3f, 3.3f}}, {{-1.1f, -2.1f, -3.1f}}, {{-1.2f, -2.2f, -3.2f}}, {{-1.3f, -2.3f, -3.3f}}}};
} // namespace

struct FlowGfwDerivedProducer {
  using CollisionsWithInfo = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultsGlobal, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>;
  using TracksWithPid = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

  enum CentEstimator : int {
    CentFT0C = 0,
    CentFT0CVariant1,
    CentFT0M,
    CentFV0A,
    CentNTPV,
    CentNGlobal,
    CentMFT
  };

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.f, "Maximum |z| of collisions written to the derived dataset (cm)"};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum track pT written to the derived dataset (GeV/c)"};
  Configurable<float> cfgPtMax{"cfgPtMax", 10.f, "Maximum track pT written to the derived dataset (GeV/c)"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum |eta| of tracks written to the derived dataset"};
  Configurable<bool> cfgSkipEmptyEvents{"cfgSkipEmptyEvents", true, "Do not write collisions that keep no tracks"};
  Configurable<bool> cfgRequireLooseTrackBits{"cfgRequireLooseTrackBits", false, "Keep only tracks that pass the loose DCA and TPC-crossed-rows bits"};
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FV0A; 4:NTPV; 5:NGlobal; 6:MFT"};
  Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 2000, "Maximum track occupancy for the occupancy bit; negative disables the cut"};

  Configurable<float> cfgDCAxyNSigma{"cfgDCAxyNSigma", 7.f, "Nominal pT-dependent DCAxy nSigma. Non-positive disables the nominal bit"};
  Configurable<float> cfgDCAxyNSigmaTight{"cfgDCAxyNSigmaTight", 5.f, "Tight pT-dependent DCAxy nSigma"};
  Configurable<float> cfgDCAxyNSigmaLoose{"cfgDCAxyNSigmaLoose", 10.f, "Loose pT-dependent DCAxy nSigma"};
  Configurable<std::string> cfgDCAxy{"cfgDCAxy", "(0.0026+0.005/(x^1.01))", "1-sigma shape of the pT-dependent DCAxy cut; a bit passes when |DCAxy| <= nSigma*shape(pT)"};
  Configurable<float> cfgDCAzNSigma{"cfgDCAzNSigma", 7.f, "Nominal DCAz nSigma relative to the 7-sigma pT shape. Non-positive disables the nominal bit"};
  Configurable<float> cfgDCAzNSigmaTight{"cfgDCAzNSigmaTight", 5.f, "Tight DCAz nSigma"};
  Configurable<float> cfgDCAzNSigmaLoose{"cfgDCAzNSigmaLoose", 10.f, "Loose DCAz nSigma"};
  Configurable<std::string> cfgDCAz{"cfgDCAz", "(0.0105+0.0350/(x^1.1))", "7-sigma pT-dependent DCAz shape; a bit passes when |DCAz| <= (nSigma/7)*shape(pT)"};
  Configurable<float> cfgNTPCXrows{"cfgNTPCXrows", 70.f, "Nominal minimum number of TPC crossed rows"};
  Configurable<float> cfgNTPCXrowsTight{"cfgNTPCXrowsTight", 80.f, "Tight minimum number of TPC crossed rows"};
  Configurable<float> cfgNTPCXrowsLoose{"cfgNTPCXrowsLoose", 60.f, "Loose minimum number of TPC crossed rows"};
  Configurable<float> cfgNTPCCls{"cfgNTPCCls", 50.f, "Nominal minimum number of found TPC clusters"};
  Configurable<float> cfgNTPCClsTight{"cfgNTPCClsTight", 60.f, "Tight minimum number of found TPC clusters"};
  Configurable<float> cfgNTPCClsLoose{"cfgNTPCClsLoose", 40.f, "Loose minimum number of found TPC clusters"};
  Configurable<float> cfgMinNITSCls{"cfgMinNITSCls", 5.f, "Minimum number of ITS clusters required to keep a track"};
  Configurable<float> cfgChi2PrTPCCls{"cfgChi2PrTPCCls", 2.5f, "Nominal maximum chi2 per TPC cluster"};
  Configurable<float> cfgChi2PrTPCClsLoose{"cfgChi2PrTPCClsLoose", 4.f, "Loose maximum chi2 per TPC cluster"};
  Configurable<float> cfgChi2PrITSCls{"cfgChi2PrITSCls", 36.f, "Maximum chi2 per ITS cluster required to keep a track"};

  Configurable<std::vector<double>> cfgMultGlobalCutPars{"cfgMultGlobalCutPars", std::vector<double>{2272.16, -76.6932, 1.01204, -0.00631545, 1.59868e-05, 136.336, -4.97006, 0.121199, -0.0015921, 7.66197e-06}, "Global vs FT0C multiplicity cut parameters"};
  Configurable<std::vector<double>> cfgMultPVCutPars{"cfgMultPVCutPars", std::vector<double>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}, "PV vs FT0C multiplicity cut parameters"};
  Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.223013, 0.715849, 0.664242, 0.0829653, -0.000503733, 1.21185e-06}, "Global vs PV multiplicity cut parameters"};
  Configurable<std::string> cfgMultCorrHighCutFunction{"cfgMultCorrHighCutFunction", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Upper multiplicity-correlation cut"};
  Configurable<std::string> cfgMultCorrLowCutFunction{"cfgMultCorrLowCutFunction", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Lower multiplicity-correlation cut"};
  Configurable<std::string> cfgMultGlobalPVCorrCutFunction{"cfgMultGlobalPVCorrCutFunction", "[0] + [1]*x + 3*([2] + [3]*x + [4]*x*x + [5]*x*x*x)", "Global vs PV multiplicity-correlation cut"};
  Configurable<float> cfgTofPtCut{"cfgTofPtCut", 0.5f, "Above this pT, PID bits stay unset when the track has no TOF match"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {NsigmaDefault.front().data(), NsigmaRows, NsigmaCols, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Nominal nSigma windows for TPC, TOF, and ITS. ITS is accepted so the analysis JSON can be reused and is not stored as a bit"};
  Configurable<float> cfgNsigmaScaleTight{"cfgNsigmaScaleTight", static_cast<float>(NsigmaScaleTight), "Scale applied to the nominal nSigma bounds for the tight PID bits"};
  Configurable<float> cfgNsigmaScaleLoose{"cfgNsigmaScaleLoose", static_cast<float>(NsigmaScaleLoose), "Scale applied to the nominal nSigma bounds for the loose PID bits"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (aod::track::pt > cfgPtMin) && (aod::track::pt < cfgPtMax) && (nabs(aod::track::eta) < cfgEtaMax);

  Produces<aod::GfwCollisions> outputCollisions;
  Produces<aod::GfwTracks> outputTracks;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::unique_ptr<TF1> fDcaXyShape;
  std::unique_ptr<TF1> fDcaZShape;
  std::array<float, NsigmaRows> tpcNsigmaCut{};
  std::array<float, NsigmaRows> tofNsigmaCut{};
  std::unique_ptr<TF1> fMultPVCutLow;
  std::unique_ptr<TF1> fMultPVCutHigh;
  std::unique_ptr<TF1> fMultCutLow;
  std::unique_ptr<TF1> fMultCutHigh;
  std::unique_ptr<TF1> fMultPVGlobalCutHigh;

  void init(InitContext const&)
  {
    const AxisSpec axisCounter{1, 0., 1., ""};
    registry.add("events", "Written derived collisions", HistType::kTH1F, {axisCounter});
    registry.add("tracks", "Written derived tracks", HistType::kTH1F, {axisCounter});

    const std::string dcaXyFormula = "[0]*" + cfgDCAxy.value;
    fDcaXyShape = std::make_unique<TF1>("gfwDerivedDcaXy", dcaXyFormula.c_str(), DcaShapeMinPt, DcaShapeMaxPt);
    fDcaXyShape->SetParameter(0, 1.f);
    fDcaZShape = std::make_unique<TF1>("gfwDerivedDcaZ", cfgDCAz.value.c_str(), DcaShapeMinPt, DcaShapeMaxPt);

    const auto nsigma = nSigmas->getData();
    for (uint32_t row = 0; row < NsigmaRows; ++row) {
      tpcNsigmaCut[row] = nsigma(row, NsigmaColTpc);
      tofNsigmaCut[row] = nsigma(row, NsigmaColTof);
    }

    constexpr double MultCorrMaxCent = 100.;
    constexpr double MultCorrMaxNch = 3000.;
    fMultPVCutLow = std::make_unique<TF1>("gfwDerivedMultPVCutLow", cfgMultCorrLowCutFunction.value.c_str(), 0., MultCorrMaxCent);
    fMultPVCutLow->SetParameters(cfgMultPVCutPars.value.data());
    fMultPVCutHigh = std::make_unique<TF1>("gfwDerivedMultPVCutHigh", cfgMultCorrHighCutFunction.value.c_str(), 0., MultCorrMaxCent);
    fMultPVCutHigh->SetParameters(cfgMultPVCutPars.value.data());
    fMultCutLow = std::make_unique<TF1>("gfwDerivedMultCutLow", cfgMultCorrLowCutFunction.value.c_str(), 0., MultCorrMaxCent);
    fMultCutLow->SetParameters(cfgMultGlobalCutPars.value.data());
    fMultCutHigh = std::make_unique<TF1>("gfwDerivedMultCutHigh", cfgMultCorrHighCutFunction.value.c_str(), 0., MultCorrMaxCent);
    fMultCutHigh->SetParameters(cfgMultGlobalCutPars.value.data());
    fMultPVGlobalCutHigh = std::make_unique<TF1>("gfwDerivedMultPVGlobalCutHigh", cfgMultGlobalPVCorrCutFunction.value.c_str(), 0., MultCorrMaxNch);
    fMultPVGlobalCutHigh->SetParameters(cfgMultGlobalPVCutPars.value.data());
  }

  template <typename TCollision>
  float centrality(TCollision const& collision) const
  {
    switch (cfgCentEstimator.value) {
      case CentFT0C:
        return collision.centFT0C();
      case CentFT0CVariant1:
        return collision.centFT0CVariant1();
      case CentFT0M:
        return collision.centFT0M();
      case CentFV0A:
        return collision.centFV0A();
      case CentNTPV:
        return collision.centNTPV();
      case CentNGlobal:
        return collision.centNGlobal();
      case CentMFT:
        return collision.centMFT();
      default:
        return collision.centFT0C();
    }
  }

  template <typename TCollision>
  bool vertexQuality(TCollision const& collision) const
  {
    if (collision.numContrib() <= 1) {
      return false;
    }
    const float zRes = std::sqrt(collision.covZZ());
    return !(zRes > VertexResMax && collision.numContrib() < VertexContribMin);
  }

  template <typename TCollision>
  bool multiplicityCorrelated(TCollision const& collision, float cent, float multTrk) const
  {
    const float multPV = collision.multNTracksPV();
    if (multPV < fMultPVCutLow->Eval(cent) || multPV > fMultPVCutHigh->Eval(cent)) {
      return false;
    }
    if (multTrk < fMultCutLow->Eval(cent) || multTrk > fMultCutHigh->Eval(cent)) {
      return false;
    }
    return multTrk <= fMultPVGlobalCutHigh->Eval(multPV);
  }

  template <typename TCollision>
  uint16_t packEventSel(TCollision const& collision, float cent, float multTrk) const
  {
    uint16_t bits = 0;
    if (collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      bits |= aod::gfwcollision::kNoSameBunchPileup;
    }
    if (collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      bits |= aod::gfwcollision::kIsGoodZvtxFT0vsPV;
    }
    if (collision.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) {
      bits |= aod::gfwcollision::kNoCollInTimeRangeStandard;
    }
    if (collision.selection_bit(aod::evsel::kNoCollInRofStandard)) {
      bits |= aod::gfwcollision::kNoCollInRofStandard;
    }
    if (collision.selection_bit(aod::evsel::kNoHighMultCollInPrevRof)) {
      bits |= aod::gfwcollision::kNoHighMultCollInPrevRof;
    }
    if (collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      bits |= aod::gfwcollision::kIsVertexITSTPC;
    }
    if (collision.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      bits |= aod::gfwcollision::kIsGoodITSLayersAll;
    }
    if (collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      bits |= aod::gfwcollision::kNoTimeFrameBorder;
    }
    if (collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      bits |= aod::gfwcollision::kNoITSROFrameBorder;
    }
    if (cfgOccupancyMax.value < 0 || collision.trackOccupancyInTimeRange() <= cfgOccupancyMax.value) {
      bits |= aod::gfwcollision::kOccupancy;
    }
    if (!collision.alias_bit(kTVXinTRD)) {
      bits |= aod::gfwcollision::kNotTVXinTRD;
    }
    if (vertexQuality(collision)) {
      bits |= aod::gfwcollision::kGoodVertex;
    }
    if (multiplicityCorrelated(collision, cent, multTrk)) {
      bits |= aod::gfwcollision::kMultCorr;
    }
    return bits;
  }

  template <typename TTrack>
  bool withinDcaXy(TTrack const& track, float nSigma) const
  {
    return std::abs(track.dcaXY()) <= nSigma * fDcaXyShape->Eval(track.pt());
  }

  template <typename TTrack>
  bool withinDcaZ(TTrack const& track, float nSigma) const
  {
    return std::abs(track.dcaZ()) <= (nSigma / DcaZReferenceNsigma) * fDcaZShape->Eval(track.pt());
  }

  template <typename TTrack>
  bool passProducerTrackCuts(TTrack const& track) const
  {
    const bool global = track.isGlobalTrack() || track.isGlobalTrackSDD() != 0;
    return global && track.itsNCls() >= cfgMinNITSCls.value && track.itsChi2NCl() <= cfgChi2PrITSCls.value;
  }

  template <typename TTrack>
  uint32_t packTrackSel(TTrack const& track) const
  {
    uint32_t bits = 0;
    if (cfgDCAxyNSigma.value <= 0.f || withinDcaXy(track, cfgDCAxyNSigma.value)) {
      bits |= aod::gfwtrack::kDCAxy;
    }
    if (withinDcaXy(track, cfgDCAxyNSigmaTight.value)) {
      bits |= aod::gfwtrack::kDCAxyTight;
    }
    if (withinDcaXy(track, cfgDCAxyNSigmaLoose.value)) {
      bits |= aod::gfwtrack::kDCAxyLoose;
    }
    if (cfgDCAzNSigma.value <= 0.f || withinDcaZ(track, cfgDCAzNSigma.value)) {
      bits |= aod::gfwtrack::kDCAz;
    }
    if (withinDcaZ(track, cfgDCAzNSigmaTight.value)) {
      bits |= aod::gfwtrack::kDCAzTight;
    }
    if (withinDcaZ(track, cfgDCAzNSigmaLoose.value)) {
      bits |= aod::gfwtrack::kDCAzLoose;
    }
    if (track.tpcNClsCrossedRows() >= cfgNTPCXrows.value) {
      bits |= aod::gfwtrack::kTPCXrows;
    }
    if (track.tpcNClsCrossedRows() >= cfgNTPCXrowsTight.value) {
      bits |= aod::gfwtrack::kTPCXrowsTight;
    }
    if (track.tpcNClsCrossedRows() >= cfgNTPCXrowsLoose.value) {
      bits |= aod::gfwtrack::kTPCXrowsLoose;
    }
    if (track.tpcNClsFound() >= cfgNTPCCls.value) {
      bits |= aod::gfwtrack::kTPCCls;
    }
    if (track.tpcNClsFound() >= cfgNTPCClsTight.value) {
      bits |= aod::gfwtrack::kTPCClsTight;
    }
    if (track.tpcNClsFound() >= cfgNTPCClsLoose.value) {
      bits |= aod::gfwtrack::kTPCClsLoose;
    }
    if (track.tpcChi2NCl() <= cfgChi2PrTPCCls.value) {
      bits |= aod::gfwtrack::kChi2TPC;
    }
    if (track.tpcChi2NCl() <= cfgChi2PrTPCClsLoose.value) {
      bits |= aod::gfwtrack::kChi2TPCLoose;
    }
    return bits;
  }

  static bool insideNsigma(float nsigma, float up, float low, float scale)
  {
    return nsigma < up * scale && nsigma > low * scale;
  }

  template <typename TTrack>
  uint32_t packPidWindow(TTrack const& track, float scale) const
  {
    if (track.pt() > cfgTofPtCut.value && !track.hasTOF()) {
      return 0;
    }
    uint32_t bits = 0;
    if (insideNsigma(track.tpcNSigmaPi(), tpcNsigmaCut[NsigmaUpPi], tpcNsigmaCut[NsigmaLowPi], scale)) {
      bits |= aod::gfwtrack::kTpcPi;
    }
    if (insideNsigma(track.tpcNSigmaKa(), tpcNsigmaCut[NsigmaUpKa], tpcNsigmaCut[NsigmaLowKa], scale)) {
      bits |= aod::gfwtrack::kTpcKa;
    }
    if (insideNsigma(track.tpcNSigmaPr(), tpcNsigmaCut[NsigmaUpPr], tpcNsigmaCut[NsigmaLowPr], scale)) {
      bits |= aod::gfwtrack::kTpcPr;
    }
    if (insideNsigma(track.tofNSigmaPi(), tofNsigmaCut[NsigmaUpPi], tofNsigmaCut[NsigmaLowPi], scale)) {
      bits |= aod::gfwtrack::kTofPi;
    }
    if (insideNsigma(track.tofNSigmaKa(), tofNsigmaCut[NsigmaUpKa], tofNsigmaCut[NsigmaLowKa], scale)) {
      bits |= aod::gfwtrack::kTofKa;
    }
    if (insideNsigma(track.tofNSigmaPr(), tofNsigmaCut[NsigmaUpPr], tofNsigmaCut[NsigmaLowPr], scale)) {
      bits |= aod::gfwtrack::kTofPr;
    }
    return bits;
  }

  template <typename TTrack>
  uint32_t packPidSel(TTrack const& track) const
  {
    const uint32_t nominal = packPidWindow(track, NsigmaScaleNominal);
    const uint32_t tight = packPidWindow(track, cfgNsigmaScaleTight.value);
    const uint32_t loose = packPidWindow(track, cfgNsigmaScaleLoose.value);
    // Nominal TPC bits start at bit 14 and nominal TOF bits at bit 23. Shifting by one window moves each detector into its tight or loose block.
    return nominal | (tight << (aod::gfwtrack::PidWindowTight * aod::gfwtrack::PidWindowBits)) | (loose << (aod::gfwtrack::PidWindowLoose * aod::gfwtrack::PidWindowBits));
  }

  template <typename TTrack>
  void fillTrack(TTrack const& track, uint32_t bits)
  {
    outputTracks(outputCollisions.lastIndex(), track.pt(), track.eta(), track.phi(),
                 track.hasTOF(), bits | packPidSel(track));
    registry.fill(HIST("tracks"), 0.5);
  }

  void process(soa::Filtered<CollisionsWithInfo>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<TracksWithPid> const& tracks)
  {
    if (!collision.sel8()) {
      return;
    }

    int nKept = 0;
    int nForMultCorr = 0;
    for (const auto& track : tracks) {
      if (!passProducerTrackCuts(track)) {
        continue;
      }
      const uint32_t bits = packTrackSel(track);
      const bool passLoose = !cfgRequireLooseTrackBits || ((bits & aod::gfwtrack::LooseSkimMask) == aod::gfwtrack::LooseSkimMask);
      if (!passLoose) {
        continue;
      }
      ++nKept;
      if ((bits & aod::gfwtrack::MultCorrTrackMask) == aod::gfwtrack::MultCorrTrackMask) {
        ++nForMultCorr;
      }
    }
    if (nKept == 0 && cfgSkipEmptyEvents) {
      return;
    }

    const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const float cent = centrality(collision);
    outputCollisions(bc.runNumber(), collision.posZ(), cent, static_cast<float>(collision.multNTracksGlobal()), bc.timestamp(), packEventSel(collision, cent, static_cast<float>(nForMultCorr)));
    registry.fill(HIST("events"), 0.5);

    for (const auto& track : tracks) {
      if (!passProducerTrackCuts(track)) {
        continue;
      }
      const uint32_t bits = packTrackSel(track);
      const bool passLoose = !cfgRequireLooseTrackBits || ((bits & aod::gfwtrack::LooseSkimMask) == aod::gfwtrack::LooseSkimMask);
      if (!passLoose) {
        continue;
      }
      fillTrack(track, bits);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwDerivedProducer>(cfgc),
  };
}
