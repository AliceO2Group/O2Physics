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

/// \file piHypertritonFemto.cxx
/// \brief LF pi-hypertriton femto pair table producer Shifted from PWGCF
/// \author zhengqing and meiyi
/// \date 2026-09-17

#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/DataModel/LFPiHypertritonFemtoTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/CallbackService.h>
#include <Framework/Configurable.h>
#include <Framework/EndOfStreamContext.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <Math/GenVector/Boost.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using HyperCandidates = aod::DataHypCandsWColl;
using HyperCandidatesMC = aod::MCHypCandsWColl;
using HadHyperCollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using HadHyperCollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection,
                                  aod::pidTPCFullPi, aod::pidTOFFullPi,
                                  aod::TOFSignal, aod::TOFEvTime>;
using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullPi, aod::pidTOFFullPi,
                                    aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;

namespace
{
constexpr std::array<float, 3> InvalidMomentum{-1.f, -1.f, -1.f};
constexpr int HyperTritonPDG = o2::constants::physics::Pdg::kHyperTriton;
// FemtoDream average-phi-star reference radii, in cm.
constexpr std::array<float, 9> CPRTPCRadii{85.f, 105.f, 125.f, 145.f, 165.f, 185.f, 205.f, 225.f, 245.f};
using PairLorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>;
} // namespace

struct PiHypertritonFemto {
  // First aggregate constructor argument, supplied by defineDataProcessing.
  std::shared_ptr<std::function<void(EndOfStreamContext&)>> mEndOfStreamCallback;
  Produces<aod::PiHypertritonFemtoTable> mOutputDataTable;
  Produces<aod::PiHypertritonFemtoTableMC> mOutputMCTable;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"eventMixing"};
    Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
    Configurable<int> settingNoMixedEvents{"settingNoMixedEvents", 5, "Number of mixed events per event"};
  } eventMixing;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"eventSelection"};
    Configurable<bool> disableITSROFCut{"disableITSROFCut", false, "Disable the ITS readout-frame border cut for data, as in hyperRecoTask"};
    Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", false, "Reject collisions sharing the same found-by-T0 bunch crossing"};
    Configurable<bool> cfgEvSelkIsGoodZvtxFT0vsPV{"cfgEvSelkIsGoodZvtxFT0vsPV", false, "Require consistent primary vertex z positions from tracking and FT0 timing"};
  } eventSelection;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"pionTrack"};
    Configurable<float> settingCutEta{"settingCutEta", 0.8f, "Maximum pion track |eta|"};
    Configurable<float> settingPtMin{"settingPtMin", 0.14f, "Minimum pion pT"};
    Configurable<float> settingPtMax{"settingPtMax", 2.5f, "Maximum pion pT"};
    Configurable<int> settingITSInnerBarrelMin{"settingITSInnerBarrelMin", 3, "Minimum pion ITS inner-barrel clusters"};
    Configurable<int> settingITSNClsMin{"settingITSNClsMin", 7, "Minimum pion ITS clusters"};
    Configurable<int> settingTPCNClsFoundMin{"settingTPCNClsFoundMin", 80, "Minimum pion found TPC clusters"};
    Configurable<int> settingTPCCrossedRowsMin{"settingTPCCrossedRowsMin", 90, "Minimum pion crossed TPC rows"};
    Configurable<float> settingDCAxyOffset{"settingDCAxyOffset", 0.004f, "Pion DCAxy offset"};
    Configurable<float> settingDCAxyPtCoeff{"settingDCAxyPtCoeff", 0.013f, "Pion DCAxy 1/pT coefficient"};
    Configurable<float> settingDCAzOffset{"settingDCAzOffset", 0.004f, "Pion DCAz offset"};
    Configurable<float> settingDCAzPtCoeff{"settingDCAzPtCoeff", 0.013f, "Pion DCAz 1/pT coefficient"};
  } pionTrack;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"pionPid"};
    Configurable<float> settingMomCombMin{"settingMomCombMin", 0.5f, "Minimum momentum to use combined TPC+TOF pion PID"};
    Configurable<float> settingTPCNsigMax{"settingTPCNsigMax", 3.0f, "Maximum pion TPC n-sigma below the TOF threshold"};
    Configurable<float> settingCombNsigMax{"settingCombNsigMax", 3.0f, "Maximum combined pion TPC+TOF n-sigma"};
    Configurable<bool> settingReqSingleNsig{"settingReqSingleNsig", false, "Also require individual TPC and TOF n-sigma cuts in combined PID"};
  } pionPid;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"CPR"};
    Configurable<float> settingClosePairDeltaPhiMax{"settingClosePairDeltaPhiMax", 0.01f, "Average delta-phi-star ellipse semiaxis for the offline close-pair flag"};
    Configurable<float> settingClosePairDeltaEtaMax{"settingClosePairDeltaEtaMax", 0.01f, "Delta-eta ellipse semiaxis for the offline close-pair flag"};
    Configurable<float> settingClosePairDistanceMax{"settingClosePairDistanceMax", 8.f, "Average 3D TPC separation below which a same-sign daughter pair is tagged (cm)"};
  } CPR;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"mc"};
    Configurable<bool> settingRequireSel8{"settingRequireSel8", false, "Additionally require sel8 for reconstructed MC; TVX and timeframe border cuts always apply"};
    Configurable<bool> settingRequireRecoMCCollisionMatch{"settingRequireRecoMCCollisionMatch", true, "Require reconstructed collision MC labels"};
  } mc;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"hypertriton"};
    Configurable<float> settingHypMassMin{"settingHypMassMin", 2.94f, "Minimum hypertriton invariant mass"};
    Configurable<float> settingHypMassMax{"settingHypMassMax", 3.10f, "Maximum hypertriton invariant mass"};
  } hypertriton;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"output"};
    Configurable<bool> settingFillTable{"settingFillTable", false, "Enable output table filling"};
  } output;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"hadHyper"};
    Configurable<bool> enableMixing{"enableMixing", true, "Build mixed-event pion-hypertriton pairs"};
    Configurable<float> maxOutputKstar{"maxOutputKstar", -1.f, "Maximum pair k* (GeV/c); negative saves all selected pairs"};
  } hadHyper;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"zorro"};
    Configurable<bool> settingSkimmedProcessing{"settingSkimmedProcessing", false, "Skimmed dataset processing"};
    Configurable<std::string> settingTriggerMask{"settingTriggerMask", "fPiHypertritonFemto", "Zorro trigger mask"};
  } zorro;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"ccdb"};
    Configurable<std::string> settingCcdburl{"settingCcdburl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository used by Zorro and the magnetic field"};
    Configurable<std::string> settingGrpmagPath{"settingGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the Run 3 magnetic field"};
  } ccdb;

  Preslice<TrackCandidates> mPerCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> mPerColMC = aod::track::collisionId;
  PresliceUnsorted<HyperCandidates> hypPerCol = o2::aod::hyperrec::collisionId;
  PresliceUnsorted<HyperCandidatesMC> hypPerColMC = o2::aod::hyperrec::collisionId;

  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "Binning for mixed-event vertex z"};
  ConfigurableAxis axisCentrality{"axisCentrality", {40, 0, 100}, "Binning for mixed-event centrality"};
  ConfigurableAxis axisPairMultiplicity{"axisPairMultiplicity", {1501, -0.5, 1500.5}, "Binning for accepted pairs per hypertriton candidate (SE and ME)"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  using HadHyperTrackInfo = std::tuple<float, float, float, int8_t, float, uint8_t, uint8_t, float, uint32_t, float, bool, float, float>;
  using HadHyperEventInfo = std::tuple<float, float, int, float, float, float, float, float>;
  using HadHyperCandidateInfo = std::tuple<bool, float, float, float, float, float, float, float, float, float, float, uint8_t, uint8_t, float, float, uint16_t, uint16_t, uint32_t, uint32_t, float, float, float>;
  using HadHyperDataInfo = decltype(std::tuple_cat(std::declval<std::tuple<bool, uint8_t, int>>(),
                                                   std::declval<HadHyperEventInfo>(),
                                                   std::declval<HadHyperCandidateInfo>(),
                                                   std::declval<HadHyperTrackInfo>()));
  using HadHyperMCInfo = std::tuple<float, float, float, float, float, float, bool, bool, bool, bool, bool, int16_t,
                                    float, bool, bool, float, float, float, float, float, float, bool, int16_t,
                                    bool, bool, bool, bool, bool, bool>;

  struct HadHyperParticleTruth {
    int64_t particleId{-1};
    int64_t collisionId{-1};
    int32_t pdgCode{0};
    float pt{-1.f};
    float eta{-999.f};
    float phi{-999.f};
    bool isPhysicalPrimary{false};
    int16_t statusCode{0};
    int16_t process{0};
    std::array<float, 3> momentum{InvalidMomentum};
  };

  struct ClosePairTrack {
    float eta{0.f};
    int8_t sign{0};
    std::array<float, CPRTPCRadii.size()> phiStar{};
    std::array<bool, CPRTPCRadii.size()> valid{};
    std::array<std::array<float, 3>, CPRTPCRadii.size()> positions{}; // Relative to the source collision PV.
    std::array<bool, CPRTPCRadii.size()> positionValid{};
  };

  struct ClosePairResult {
    float deltaEta{0.f};
    float deltaPhiStar{0.f};
    int validRadii{0};
    float averageDistance{0.f};
    int validPositionRadii{0};
    bool angularRejected{false};
    bool distanceRejected{false};
    bool angularUnavailable{false};
    bool distanceUnavailable{false};
  };

  struct ClosePairQA {
    std::shared_ptr<TH2> before;
    std::shared_ptr<TH2> after;
    std::shared_ptr<TH1> validRadii;
    std::shared_ptr<TH1> distanceBefore;
    std::shared_ptr<TH1> distanceAfter;
    std::shared_ptr<TH1> validPositionRadii;
    std::shared_ptr<TH2> decisions;
  };

  struct HadHyperCandidate {
    HadHyperCandidateInfo info{};
    HadHyperParticleTruth hyperTruth{};
    HadHyperParticleTruth heTruth{};
    HadHyperParticleTruth decayPionTruth{};
    int64_t heTrackId{-1};
    int64_t pionTrackId{-1};
    int16_t statusCode{0};
    bool isReco{true};
    bool isSignal{false};
    bool isRecoMCCollision{false};
    bool isSurvEvSel{false};
    bool isTwoBodyDecay{false};
    float genPt{-1.f};
    float genEta{-1.f};
    float genPhi{-1.f};
    float genPtHe3{-1.f};
    std::array<float, 3> genDecVtx{InvalidMomentum};
    std::array<float, 3> momentum{};
    std::array<float, 3> genMomentum{InvalidMomentum};
    float mass{0.f};
    bool isMatter{false};
    std::array<ClosePairTrack, 2> daughterClosePairTracks{}; // He3, decay pion.
    uint64_t mixedPairCount{0};                              // Accumulated over the complete residence in the mixing pool.

    [[nodiscard]] float pt() const { return std::hypot(momentum[0], momentum[1]); }
    [[nodiscard]] float eta() const
    {
      const float transverseMomentum = pt();
      return transverseMomentum > 0.f ? std::asinh(momentum[2] / transverseMomentum) : 999.f;
    }
    [[nodiscard]] float phi() const { return std::atan2(momentum[1], momentum[0]); }
    [[nodiscard]] int8_t sign() const { return isMatter ? 1 : -1; }
  };

  struct HadHyperHadron {
    HadHyperTrackInfo info{};
    HadHyperParticleTruth truth{};
    std::vector<int64_t> motherIds;
    int64_t sourceId{-1};
    std::array<float, 3> momentum{};
    int8_t signValue{0};
    ClosePairTrack closePairTrack{};

    [[nodiscard]] float eta() const { return std::get<1>(info); }
    [[nodiscard]] float phi() const { return std::get<2>(info); }
    [[nodiscard]] int8_t sign() const { return signValue; }
  };

  struct HadHyperEvent {
    HadHyperEventInfo info{};
    uint64_t sourceFrameId{0}; // Task-local scope for track and MC indices stored in the mixing pool.
    int64_t mcCollisionId{-1};
    bool hasMCCollision{false};
    float centrality{0.f};
    std::vector<HadHyperCandidate> candidates;
    std::vector<HadHyperHadron> hadrons;
  };

  std::unordered_map<int, std::deque<HadHyperEvent>> mMixingPools;
  std::unordered_map<int, std::deque<HadHyperEvent>> mMCMixingPools;
  int mMixingRunNumber{-1};
  int mMCMixingRunNumber{-1};
  uint64_t mSourceFrameId{0};
  int mRunNumber{-1};
  float mMagneticFieldTesla{0.f};
  std::array<std::array<ClosePairQA, 2>, 2> mClosePairQA{}; // SE/ME, He3/decay pion.
  Service<o2::ccdb::BasicCCDBManager> mCcdb{};
  Zorro mZorro;
  OutputObj<ZorroSummary> mZorroSummary{"zorroSummary"};

  HistogramRegistry mQaRegistry{
    "QA",
    {{"hEvents", "Event selection cut-flow;Selection step;Events", {HistType::kTH1F, {{11, -0.5, 10.5}}}},
     {"hVtxZ", "Vertex distribution in Z;Z (cm);Entries", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"hNcontributor", "Number of primary vertex contributors;contributors;Entries", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
     {"hPionTPCNSigmaPreselection", "Pion TPC n#sigma before PID;signed #it{p}_{TPC};n#sigma_{TPC}^{#pi}", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {400, -10.0f, 10.0f}}}},
     {"hPionTOFNSigmaPreselection", "Pion TOF n#sigma before PID;signed #it{p}_{T};n#sigma_{TOF}^{#pi}", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"hPionCombNSigmaPreselection", "Pion combined n#sigma before PID;signed #it{p}_{T};n#sigma_{comb}^{#pi}", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"hPionTPCNSigma", "Selected pion TPC n#sigma;signed #it{p}_{T};n#sigma_{TPC}^{#pi}", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"hPionTOFNSigma", "Selected pion TOF n#sigma;signed #it{p}_{T};n#sigma_{TOF}^{#pi}", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"hPionCombNSigma", "Selected pion combined n#sigma;signed #it{p}_{T};n#sigma_{comb}^{#pi}", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"hPionTPC", "Selected pion TPC d#it{E}/d#it{x};signed #it{p}_{TPC};TPC signal", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"hHypHe3TPCNSigma", "Selected hypertriton candidates;#it{p}_{T}^{^{3}He};n#sigma_{TPC}^{^{3}He}", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"hHypHe3TPCMomPreselection", "Hypertriton candidate He3 TPC momentum before mass cut;#it{p}_{TPC};Entries", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"hHypHe3TPCMom", "Hypertriton candidate He3 TPC momentum after mass cut;#it{p}_{TPC};Entries", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  HistogramRegistry hadHyperRegistry{
    "hadHyperRegistry",
    {{"hSE", "Raw same-event pairs;k*;Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"hME", "Raw mixed-event pairs;k*;Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"hMass", "Unique selected candidates;mass He3-pi;Entries", {HistType::kTH1D, {{160, 2.94, 3.10}}}},
     {"hSelfPairs", "Rejected daughter reuse;reason (0=He3,1=pion,2=truth);k*", {HistType::kTH2F, {{3, -0.5, 2.5}, {300, 0., 3.}}}},
     {"hSameEventSelfPairs", "Same-event rejected daughter reuse;reason (0=He3,1=pion,2=truth);k*", {HistType::kTH2F, {{3, -0.5, 2.5}, {300, 0., 3.}}}},
     {"hMixEventDeltaPosZVsCent", "Mixed-event #Delta z_{vtx} vs hypertriton centrality;hypertriton CentFT0C;#Delta z_{vtx}", {HistType::kTH2F, {{100, 0., 100.}, {120, -30., 30.}}}},
     {"hMixEventDeltaCentFT0CVsCent", "Mixed-event #Delta CentFT0C vs hypertriton centrality;hypertriton CentFT0C;#Delta CentFT0C", {HistType::kTH2F, {{100, 0., 100.}, {200, -100., 100.}}}},
     {"hMixingDepth", "Available partner events;events;anchor events", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
     {"hPoolFlow", "Mixing pool selection;0=selected,1=in range;events", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"hDaughterHeTPC", "He3 daughter QA per candidate;TPC rigidity;TPC signal", {HistType::kTH2F, {{100, 0., 10.}, {300, 0., 1500.}}}},
     {"hDaughterPiTPC", "Pion daughter QA per candidate;TPC rigidity;TPC signal", {HistType::kTH2F, {{100, 0., 10.}, {300, 0., 1500.}}}},
     {"hPionTPC", "Pair pion QA per selected track;TPC rigidity;TPC signal", {HistType::kTH2F, {{100, 0., 10.}, {300, 0., 1500.}}}},
     {"MC/SE/hKstarRecVsGenHyperReco", "Pion-hypertriton k* response for reconstructed true hypertritons;generated k* (GeV/#it{c});reconstructed k* (GeV/#it{c})", {HistType::kTH2F, {{300, 0., 3.}, {300, 0., 3.}}}},
     {"MC/SE/hKstarResolutionHyperReco", "Pion-hypertriton k* resolution for reconstructed true hypertritons;reconstructed k* (GeV/#it{c});k*_{reco}-k*_{gen} (GeV/#it{c})", {HistType::kTH2F, {{300, 0., 3.}, {200, -0.2, 0.2}}}},
     {"MC/SE/hPrimaryPionVsKstarDen", "Truth-matched selected pion pair-weighted denominator;k* (GeV/#it{c});Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"MC/SE/hPrimaryPionVsKstarNum", "Physical-primary selected pion pair-weighted numerator;k* (GeV/#it{c});Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"MC/SE/hPrimaryPionVsCentDen", "Truth-matched selected pion pair-weighted denominator;hypertriton CentFT0C;Entries", {HistType::kTH1D, {{100, 0., 100.}}}},
     {"MC/SE/hPrimaryPionVsCentNum", "Physical-primary selected pion pair-weighted numerator;hypertriton CentFT0C;Entries", {HistType::kTH1D, {{100, 0., 100.}}}}}};

  void init(o2::framework::InitContext&);
  template <bool isMC>
  void initCCDB(const aod::BCsWithTimestamps::iterator& bc);
  template <bool isMC, typename Tcollision>
  bool selectCollision(const Tcollision& collision, const aod::BCsWithTimestamps& bcs);
  template <typename Ttrack>
  bool selectPionTrack(const Ttrack& candidate) const;
  template <typename Ttrack>
  bool selectionPIDPion(const Ttrack& candidate);
  template <typename Tcandidate>
  bool selectHyperCandidate(const Tcandidate& candidate);
  template <bool isMC, typename Tcollisions, typename Ttracks, typename Tcandidates>
  void processPairs(const Tcollisions& collisions, const Ttracks& tracks, const Tcandidates& candidates,
                    const aod::BCsWithTimestamps& bcs,
                    std::unordered_map<int, std::deque<HadHyperEvent>>& mixingPools,
                    int& mixingRunNumber);

  template <typename Tcandidate>
  float computeHyperCandidateMass(const Tcandidate& candidate) const
  {
    const std::array<float, 3> heMomentum{candidate.ptHe3() * std::cos(candidate.phiHe3()), candidate.ptHe3() * std::sin(candidate.phiHe3()), candidate.ptHe3() * std::sinh(candidate.etaHe3())};
    const std::array<float, 3> pionMomentum{candidate.ptPi() * std::cos(candidate.phiPi()), candidate.ptPi() * std::sin(candidate.phiPi()), candidate.ptPi() * std::sinh(candidate.etaPi())};
    return RecoDecay::m(std::array{heMomentum, pionMomentum},
                        std::array{static_cast<float>(o2::constants::physics::MassHelium3),
                                   static_cast<float>(o2::constants::physics::MassPiPlus)});
  }

  float wrapDeltaPhi(float dphi) const { return std::atan2(std::sin(dphi), std::cos(dphi)); }

  template <typename Ttrack>
  ClosePairTrack makeClosePairTrack(const Ttrack& track, const std::array<float, 3>& primaryVertex) const
  {
    ClosePairTrack result;
    result.eta = track.eta();
    result.sign = track.sign();
    if (!std::isfinite(result.eta) || !std::isfinite(track.phi()) || !std::isfinite(track.signed1Pt()) || track.signed1Pt() == 0.f) {
      return result;
    }
    // FemtoDream helix approximation, using the original daughter track, not the mother.
    // AO2D signed1Pt is signed inverse rigidity: He3's |Z|=2 is already included.
    // Do not multiply it by two again. Cache each source event's field-dependent values.
    constexpr float CurvatureFactor = 0.3f * 0.01f / 2.f; // B in T, radius in cm.
    for (std::size_t i = 0; i < CPRTPCRadii.size(); ++i) {
      const float argument = CurvatureFactor * mMagneticFieldTesla * CPRTPCRadii[i] * track.signed1Pt();
      if (std::isfinite(argument) && std::abs(argument) < 1.f) {
        result.phiStar[i] = track.phi() - std::asin(argument);
        result.valid[i] = true;
      }
    }
    // Uniform-field extrapolation from the measured track state, retaining daughter displacement.
    // This is not material-aware transport. getTrackPar keeps q/pT (including He3 rigidity).
    const auto trackPar = getTrackPar(track);
    constexpr float TeslaToKilogauss = 10.f;
    constexpr float TPCHalfLength = 250.f; // cm
    const float fieldKilogauss = TeslaToKilogauss * mMagneticFieldTesla;
    for (std::size_t i = 0; i < CPRTPCRadii.size(); ++i) {
      float localX = 0.f;
      if (!trackPar.getXatLabR(CPRTPCRadii[i], localX, fieldKilogauss, o2::track::DirOutward)) {
        continue;
      }
      bool propagated = false;
      const auto point = trackPar.getXYZGloAt(localX, fieldKilogauss, propagated);
      if (!propagated || !std::isfinite(point.X()) || !std::isfinite(point.Y()) || !std::isfinite(point.Z()) || std::abs(point.Z()) > TPCHalfLength) {
        continue;
      }
      // Apply the same source-PV convention to SE and ME; in SE the common PV cancels.
      result.positions[i] = {point.X() - primaryVertex[0], point.Y() - primaryVertex[1], point.Z() - primaryVertex[2]};
      result.positionValid[i] = true;
    }
    return result;
  }

  ClosePairResult evaluateClosePair(const ClosePairTrack& pion, const ClosePairTrack& daughter) const
  {
    ClosePairResult result;
    result.deltaEta = pion.eta - daughter.eta;
    for (std::size_t i = 0; i < CPRTPCRadii.size(); ++i) {
      if (pion.valid[i] && daughter.valid[i]) {
        result.deltaPhiStar += wrapDeltaPhi(pion.phiStar[i] - daughter.phiStar[i]);
        ++result.validRadii;
      }
      if (pion.positionValid[i] && daughter.positionValid[i]) {
        const auto& pionPosition = pion.positions[i];
        const auto& daughterPosition = daughter.positions[i];
        result.averageDistance += std::hypot(pionPosition[0] - daughterPosition[0],
                                             pionPosition[1] - daughterPosition[1],
                                             pionPosition[2] - daughterPosition[2]);
        ++result.validPositionRadii;
      }
    }
    const bool sameSign = pion.sign * daughter.sign > 0;
    result.angularUnavailable = sameSign && result.validRadii == 0;
    result.distanceUnavailable = sameSign && result.validPositionRadii == 0;
    // Undefined comparisons set separate status bits, never a rejection bit.
    if (result.validRadii > 0) {
      result.deltaPhiStar /= result.validRadii;
      const float scaledEta = result.deltaEta / CPR.settingClosePairDeltaEtaMax.value;
      const float scaledPhi = result.deltaPhiStar / CPR.settingClosePairDeltaPhiMax.value;
      result.angularRejected = sameSign && scaledEta * scaledEta + scaledPhi * scaledPhi < 1.f;
    }
    if (result.validPositionRadii > 0) {
      result.averageDistance /= result.validPositionRadii;
      result.distanceRejected = sameSign && result.averageDistance < CPR.settingClosePairDistanceMax.value;
    }
    return result;
  }

  uint8_t flagClosePair(const HadHyperCandidate& candidate, const HadHyperHadron& pion, bool mixed)
  {
    const std::array<ClosePairResult, 2> results{
      evaluateClosePair(pion.closePairTrack, candidate.daughterClosePairTracks[0]),
      evaluateClosePair(pion.closePairTrack, candidate.daughterClosePairTracks[1])};
    uint8_t rejectionFlags = 0;
    for (const auto& result : results) {
      if (result.angularRejected) {
        rejectionFlags |= aod::pihypertritonfemto::ClosePairAngular;
      }
      if (result.distanceRejected) {
        rejectionFlags |= aod::pihypertritonfemto::ClosePairDistance;
      }
      if (result.angularUnavailable) {
        rejectionFlags |= aod::pihypertritonfemto::ClosePairAngularUnavailable;
      }
      if (result.distanceUnavailable) {
        rejectionFlags |= aod::pihypertritonfemto::ClosePairDistanceUnavailable;
      }
    }
    for (std::size_t i = 0; i < results.size(); ++i) {
      const auto& result = results[i];
      const auto& qa = mClosePairQA[mixed ? 1 : 0][i];
      qa.validRadii->Fill(result.validRadii);
      qa.validPositionRadii->Fill(result.validPositionRadii);
      if (result.validRadii > 0 && result.validPositionRadii > 0) {
        qa.decisions->Fill(result.angularRejected, result.distanceRejected);
      }
      if (result.validPositionRadii > 0) {
        qa.distanceBefore->Fill(result.averageDistance);
        if (rejectionFlags == 0) {
          qa.distanceAfter->Fill(result.averageDistance);
        }
      }
      if (result.validRadii > 0) {
        qa.before->Fill(result.deltaEta, result.deltaPhiStar);
        if (rejectionFlags == 0) {
          qa.after->Fill(result.deltaEta, result.deltaPhiStar);
        }
      }
    }
    return rejectionFlags;
  }

  float computePairKstar(const std::array<float, 3>& momPion, const std::array<float, 3>& momHyper) const
  {
    const PairLorentzVector vecPion(momPion[0], momPion[1], momPion[2], o2::constants::physics::MassPiPlus);
    const PairLorentzVector vecHyper(momHyper[0], momHyper[1], momHyper[2], o2::constants::physics::MassHyperTriton);
    const PairLorentzVector trackSum = vecPion + vecHyper;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    PairLorentzVector pionCMS(vecPion);
    PairLorentzVector hyperCMS(vecHyper);
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    pionCMS = boostPRF(pionCMS);
    hyperCMS = boostPRF(hyperCMS);
    return 0.5f * (pionCMS - hyperCMS).P();
  }

  template <typename Ttrack>
  HadHyperTrackInfo makePionInfo(const Ttrack& track)
  {
    constexpr float InvalidPID = -999.f;
    return {track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign()),
            track.dcaXY(), track.tpcNClsCrossedRows(), track.tpcNClsPID(), track.tpcChi2NCl(),
            track.itsClusterSizes(), track.itsChi2NCl(), track.hasTOF(),
            track.tpcNSigmaPi(), track.hasTOF() ? track.tofNSigmaPi() : InvalidPID};
  }

  template <typename Tcollision>
  HadHyperEventInfo makeEventInfo(const Tcollision& collision)
  {
    return {collision.posZ(), collision.centFT0C(),
            collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(),
            collision.multFT0C(),
            collision.posX(), collision.posY(), collision.posZ()};
  }

  template <typename Tparticle>
  HadHyperParticleTruth makeParticleTruth(const Tparticle& particle)
  {
    return {particle.globalIndex(), particle.mcCollisionId(), particle.pdgCode(), particle.pt(), particle.eta(), particle.phi(), particle.isPhysicalPrimary(), static_cast<int16_t>(particle.statusCode()), static_cast<int16_t>(particle.getProcess()), {particle.px(), particle.py(), particle.pz()}};
  }

  template <bool isMC, typename Ttrack>
  HadHyperHadron makePion(const Ttrack& track, const std::array<float, 3>& primaryVertex)
  {
    HadHyperHadron pion{makePionInfo(track), {}, {}, track.globalIndex(), {track.px(), track.py(), track.pz()}, static_cast<int8_t>(track.sign())};
    pion.closePairTrack = makeClosePairTrack(track, primaryVertex);
    if constexpr (isMC) {
      if (track.has_mcParticle()) {
        const auto particle = track.template mcParticle_as<aod::McParticles>();
        pion.truth = makeParticleTruth(particle);
        if (particle.has_mothers()) {
          for (const auto& mother : particle.template mothers_as<aod::McParticles>()) {
            pion.motherIds.push_back(mother.globalIndex());
          }
        }
      }
    }
    return pion;
  }

  template <bool isMC, typename Tcandidate, typename Ttracks>
  HadHyperCandidate makeHyperCandidate(const Tcandidate& candidate, const Ttracks& tracks, const std::array<float, 3>& primaryVertex)
  {
    const auto he = tracks.rawIteratorAt(candidate.heTrackId() - tracks.offset());
    const auto decayPion = tracks.rawIteratorAt(candidate.piTrackId() - tracks.offset());
    const std::array<float, 3> heMomentum{candidate.ptHe3() * std::cos(candidate.phiHe3()), candidate.ptHe3() * std::sin(candidate.phiHe3()), candidate.ptHe3() * std::sinh(candidate.etaHe3())};
    const std::array<float, 3> decayPionMomentum{candidate.ptPi() * std::cos(candidate.phiPi()), candidate.ptPi() * std::sin(candidate.phiPi()), candidate.ptPi() * std::sinh(candidate.etaPi())};

    HadHyperCandidate result;
    for (std::size_t i = 0; i < result.momentum.size(); ++i) {
      result.momentum[i] = heMomentum[i] + decayPionMomentum[i];
    }
    result.mass = computeHyperCandidateMass(candidate);
    result.heTrackId = candidate.heTrackId();
    result.pionTrackId = candidate.piTrackId();
    result.isMatter = candidate.isMatter();
    result.daughterClosePairTracks = {makeClosePairTrack(he, primaryVertex), makeClosePairTrack(decayPion, primaryVertex)};

    if constexpr (isMC) {
      if (he.has_mcParticle()) {
        result.heTruth = makeParticleTruth(he.template mcParticle_as<aod::McParticles>());
      }
      if (decayPion.has_mcParticle()) {
        result.decayPionTruth = makeParticleTruth(decayPion.template mcParticle_as<aod::McParticles>());
      }
      if (he.has_mcParticle() && decayPion.has_mcParticle()) {
        const auto heParticle = he.template mcParticle_as<aod::McParticles>();
        const auto decayPionParticle = decayPion.template mcParticle_as<aod::McParticles>();
        if (heParticle.has_mothers() && decayPionParticle.has_mothers()) {
          for (const auto& heMother : heParticle.template mothers_as<aod::McParticles>()) {
            for (const auto& pionMother : decayPionParticle.template mothers_as<aod::McParticles>()) {
              if (heMother.globalIndex() == pionMother.globalIndex() && std::abs(heMother.pdgCode()) == HyperTritonPDG) {
                result.hyperTruth = makeParticleTruth(heMother);
              }
            }
          }
        }
      }
    }

    result.info = {candidate.isMatter(), candidate.ptHe3(), candidate.etaHe3(), candidate.phiHe3(),
                   candidate.ptPi(), candidate.etaPi(), candidate.phiPi(),
                   candidate.dcaV0Daug(), candidate.dcaHe(), candidate.dcaPi(),
                   candidate.nSigmaHe(),
                   candidate.nTPCCrossedRowsHe(), candidate.nTPCCrossedRowsPi(),
                   candidate.tpcMomHe(), candidate.tpcMomPi(), candidate.tpcSignalHe(), candidate.tpcSignalPi(),
                   candidate.itsClusterSizesHe(), candidate.itsClusterSizesPi(),
                   candidate.xDecVtx(), candidate.yDecVtx(), candidate.zDecVtx()};

    if constexpr (isMC) {
      result.statusCode = static_cast<int16_t>(candidate.statusCode());
      result.isReco = candidate.isReco();
      result.isSignal = candidate.isSignal();
      result.isRecoMCCollision = candidate.isRecoMCCollision();
      result.isSurvEvSel = candidate.isSurvEvSel();
      result.isTwoBodyDecay = candidate.isTwoBodyDecay();
      result.genPt = candidate.genPt();
      result.genEta = candidate.genEta();
      result.genPhi = candidate.genPhi();
      result.genPtHe3 = candidate.genPtHe3();
      result.genDecVtx = {candidate.genXDecVtx(), candidate.genYDecVtx(), candidate.genZDecVtx()};
      const float absGenPt = std::abs(candidate.genPt());
      result.genMomentum = {absGenPt * std::cos(candidate.genPhi()), absGenPt * std::sin(candidate.genPhi()), absGenPt * std::sinh(candidate.genEta())};
    }

    hadHyperRegistry.fill(HIST("hMass"), result.mass);
    hadHyperRegistry.fill(HIST("hDaughterHeTPC"), he.tpcInnerParam(), he.tpcSignal());
    hadHyperRegistry.fill(HIST("hDaughterPiTPC"), decayPion.tpcInnerParam(), decayPion.tpcSignal());
    return result;
  }

  HadHyperDataInfo makeDataInfo(const HadHyperCandidate& candidate, const HadHyperHadron& pion, const HadHyperEvent& hyperEvent, bool mixed, uint8_t closePairRejected, int mixingDepth)
  {
    return std::tuple_cat(std::make_tuple(mixed, closePairRejected, mixingDepth), hyperEvent.info, candidate.info, pion.info);
  }

  bool hasTruthMother(const HadHyperHadron& pion, int64_t motherId) const
  {
    return motherId >= 0 && std::find(pion.motherIds.begin(), pion.motherIds.end(), motherId) != pion.motherIds.end();
  }

  bool isRecoSelfCorrelation(const HadHyperCandidate& candidate, const HadHyperHadron& pion,
                             const HadHyperEvent& hyperEvent, const HadHyperEvent& pionEvent) const
  {
    return hyperEvent.sourceFrameId == pionEvent.sourceFrameId && pion.sourceId >= 0 &&
           (pion.sourceId == candidate.heTrackId || pion.sourceId == candidate.pionTrackId);
  }

  bool isTruthSelfCorrelation(const HadHyperCandidate& candidate, const HadHyperHadron& pion,
                              const HadHyperEvent& hyperEvent, const HadHyperEvent& pionEvent) const
  {
    return hyperEvent.sourceFrameId == pionEvent.sourceFrameId && pion.truth.particleId >= 0 &&
           (pion.truth.particleId == candidate.heTruth.particleId ||
            pion.truth.particleId == candidate.decayPionTruth.particleId ||
            hasTruthMother(pion, candidate.hyperTruth.particleId));
  }

  HadHyperMCInfo makeMCInfo(const HadHyperCandidate& candidate, const HadHyperHadron& pion,
                            const HadHyperEvent& hyperEvent, const HadHyperEvent& pionEvent) const
  {
    const bool sameMCCollision = hyperEvent.sourceFrameId == pionEvent.sourceFrameId &&
                                 candidate.hyperTruth.collisionId >= 0 && pion.truth.collisionId == candidate.hyperTruth.collisionId;
    const bool matchesHypRecoMCCollision = hyperEvent.hasMCCollision && candidate.hyperTruth.collisionId == hyperEvent.mcCollisionId;
    const bool matchesPairRecoMCCollision = pionEvent.hasMCCollision && pion.truth.collisionId == pionEvent.mcCollisionId;
    const bool truthSelfCorrelation = isTruthSelfCorrelation(candidate, pion, hyperEvent, pionEvent);
    const bool hadIsTruePion = pion.truth.particleId >= 0 && std::abs(pion.truth.pdgCode) == static_cast<int>(PDG_t::kPiPlus);
    const bool isTruePrimaryPiHyperPair = candidate.isSignal &&
                                          std::abs(candidate.hyperTruth.pdgCode) == HyperTritonPDG &&
                                          hadIsTruePion &&
                                          pion.truth.isPhysicalPrimary && sameMCCollision && !truthSelfCorrelation;
    return {candidate.genPt, candidate.genEta, candidate.genPhi,
            candidate.genDecVtx[0], candidate.genDecVtx[1], candidate.genDecVtx[2],
            candidate.isReco, candidate.isSignal, candidate.isRecoMCCollision, candidate.isSurvEvSel, candidate.isTwoBodyDecay, candidate.statusCode,
            candidate.heTruth.pt, candidate.heTruth.isPhysicalPrimary, candidate.decayPionTruth.isPhysicalPrimary,
            std::get<0>(pion.info), std::get<1>(pion.info), std::get<2>(pion.info),
            pion.truth.pt, pion.truth.eta, pion.truth.phi, pion.truth.isPhysicalPrimary, pion.truth.process,
            sameMCCollision, matchesHypRecoMCCollision, matchesPairRecoMCCollision, truthSelfCorrelation, isTruePrimaryPiHyperPair, hadIsTruePion};
  }

  void fillMCQA(const HadHyperCandidate& candidate, const HadHyperHadron& pion, const HadHyperEvent& hyperEvent, float kstar, bool mixed)
  {
    if (!candidate.isReco || !candidate.isSignal || std::abs(candidate.hyperTruth.pdgCode) != HyperTritonPDG || std::abs(pion.truth.pdgCode) != static_cast<int>(PDG_t::kPiPlus)) {
      return;
    }

    const float kstarMC = computePairKstar(pion.truth.momentum, candidate.genMomentum);
    const auto fill = [&](const auto& prefix) {
      if (std::isfinite(kstarMC)) {
        hadHyperRegistry.fill(prefix + HIST("hKstarRecVsGenHyperReco"), kstarMC, kstar);
        hadHyperRegistry.fill(prefix + HIST("hKstarResolutionHyperReco"), kstar, kstar - kstarMC);
      }

      hadHyperRegistry.fill(prefix + HIST("hPrimaryPionVsKstarDen"), kstar);
      hadHyperRegistry.fill(prefix + HIST("hPrimaryPionVsCentDen"), hyperEvent.centrality);
      if (pion.truth.isPhysicalPrimary) {
        hadHyperRegistry.fill(prefix + HIST("hPrimaryPionVsKstarNum"), kstar);
        hadHyperRegistry.fill(prefix + HIST("hPrimaryPionVsCentNum"), hyperEvent.centrality);
      }
    };
    if (mixed) {
      fill(HIST("MC/ME/"));
    } else {
      fill(HIST("MC/SE/"));
    }
  }

  template <bool isMC>
  bool fillPair(const HadHyperCandidate& candidate, const HadHyperHadron& pion, const HadHyperEvent& hyperEvent, const HadHyperEvent& pionEvent, bool mixed, int mixingDepth)
  {
    const float kstar = computePairKstar(pion.momentum, candidate.momentum);
    if (isRecoSelfCorrelation(candidate, pion, hyperEvent, pionEvent)) {
      const int reason = pion.sourceId == candidate.heTrackId ? 0 : 1;
      hadHyperRegistry.fill(HIST("hSelfPairs"), reason, kstar);
      if (!mixed) {
        hadHyperRegistry.fill(HIST("hSameEventSelfPairs"), reason, kstar);
      }
      return false;
    }
    if constexpr (isMC) {
      if (isTruthSelfCorrelation(candidate, pion, hyperEvent, pionEvent)) {
        hadHyperRegistry.fill(HIST("hSelfPairs"), 2, kstar);
        if (!mixed) {
          hadHyperRegistry.fill(HIST("hSameEventSelfPairs"), 2, kstar);
        }
        return false;
      }
    }

    if (!std::isfinite(kstar)) {
      return false;
    }
    // CPR is an offline flag only: do not remove rows or change raw-pair/MC QA counts.
    const uint8_t closePairRejected = flagClosePair(candidate, pion, mixed);

    if (mixed) {
      hadHyperRegistry.fill(HIST("hME"), kstar);
    } else {
      hadHyperRegistry.fill(HIST("hSE"), kstar);
    }
    if constexpr (isMC) {
      fillMCQA(candidate, pion, hyperEvent, kstar, mixed);
    }

    if (!output.settingFillTable || (hadHyper.maxOutputKstar.value > 0.f && kstar >= hadHyper.maxOutputKstar.value)) {
      return true;
    }

    const auto dataInfo = makeDataInfo(candidate, pion, hyperEvent, mixed, closePairRejected, mixingDepth);
    if constexpr (isMC) {
      std::apply([this](const auto&... columns) { mOutputMCTable(columns...); }, std::tuple_cat(dataInfo, makeMCInfo(candidate, pion, hyperEvent, pionEvent)));
    } else {
      std::apply([this](const auto&... columns) { mOutputDataTable(columns...); }, dataInfo);
    }
    return true;
  }

  template <typename Ttracks>
  bool hasValidHyperDaughterIndices(const Ttracks& tracks, int64_t heTrackId, int64_t pionTrackId) const
  {
    const auto first = static_cast<int64_t>(tracks.offset());
    const auto last = first + static_cast<int64_t>(tracks.size());
    return heTrackId >= first && heTrackId < last && pionTrackId >= first && pionTrackId < last;
  }

  template <bool isMC, typename Ttracks>
  void collectPions(HadHyperEvent& event, const Ttracks& eventTracks)
  {
    const std::array<float, 3> primaryVertex{std::get<5>(event.info), std::get<6>(event.info), std::get<7>(event.info)};
    for (const auto& track : eventTracks) {
      if (!selectPionTrack(track) || !selectionPIDPion(track)) {
        continue;
      }
      event.hadrons.push_back(makePion<isMC>(track, primaryVertex));
      hadHyperRegistry.fill(HIST("hPionTPC"), track.tpcInnerParam(), track.tpcSignal());
    }
  }

  template <bool isMC, typename Tcandidates, typename Ttracks>
  void collectHyperCandidates(HadHyperEvent& event, const Tcandidates& eventCandidates, const Ttracks& tracks)
  {
    const std::array<float, 3> primaryVertex{std::get<5>(event.info), std::get<6>(event.info), std::get<7>(event.info)};
    for (const auto& candidate : eventCandidates) {
      if constexpr (isMC) {
        if (!candidate.isReco()) {
          continue;
        }
      }
      if (!hasValidHyperDaughterIndices(tracks, candidate.heTrackId(), candidate.piTrackId())) {
        if constexpr (isMC) {
          continue;
        } else {
          LOG(fatal) << "Hypertriton daughter indices must reference the input Tracks table";
        }
      }
      if (!selectHyperCandidate(candidate)) {
        continue;
      }
      event.candidates.push_back(makeHyperCandidate<isMC>(candidate, tracks, primaryVertex));
    }
  }

  template <bool isMC, typename Tcollision, typename Ttracks, typename Tcandidates>
  HadHyperEvent buildEvent(const Tcollision& collision, const Ttracks& tracks, const Tcandidates& candidates)
  {
    HadHyperEvent event;
    event.info = makeEventInfo(collision);
    event.centrality = collision.centFT0C();
    if constexpr (isMC) {
      event.hasMCCollision = collision.has_mcCollision();
      event.mcCollisionId = collision.has_mcCollision() ? collision.mcCollisionId() : -1;
    }

    if constexpr (isMC) {
      const auto eventTracks = tracks.sliceBy(mPerColMC, collision.globalIndex());
      collectPions<isMC>(event, eventTracks);
      const auto eventCandidates = candidates.sliceBy(hypPerColMC, collision.globalIndex());
      collectHyperCandidates<isMC>(event, eventCandidates, tracks);
    } else {
      const auto eventTracks = tracks.sliceBy(mPerCol, collision.globalIndex());
      collectPions<isMC>(event, eventTracks);
      const auto eventCandidates = candidates.sliceBy(hypPerCol, collision.globalIndex());
      collectHyperCandidates<isMC>(event, eventCandidates, tracks);
    }
    return event;
  }

  template <bool isMC>
  void fillSameEventPairs(const HadHyperEvent& event)
  {
    for (const auto& candidate : event.candidates) {
      int acceptedPairs = 0;
      for (const auto& pion : event.hadrons) {
        if (fillPair<isMC>(candidate, pion, event, event, false, 0)) {
          ++acceptedPairs;
        }
      }
      hadHyperRegistry.fill(HIST("hCandidatePairMultiplicitySE"), acceptedPairs, event.centrality);
    }
  }

  template <bool isMC>
  void fillMixedEventPairs(HadHyperEvent& currentEvent, std::deque<HadHyperEvent>& pool)
  {
    const int depth = static_cast<int>(pool.size());
    hadHyperRegistry.fill(HIST("hMixingDepth"), depth);
    for (std::size_t partnerIndex = 0; partnerIndex < pool.size(); ++partnerIndex) {
      auto& partner = pool[partnerIndex];
      const float currentPosZ = std::get<0>(currentEvent.info);
      const float partnerPosZ = std::get<0>(partner.info);
      hadHyperRegistry.fill(HIST("hMixEventDeltaPosZVsCent"), currentEvent.centrality, currentPosZ - partnerPosZ);
      hadHyperRegistry.fill(HIST("hMixEventDeltaCentFT0CVsCent"), currentEvent.centrality, currentEvent.centrality - partner.centrality);
      hadHyperRegistry.fill(HIST("hMixEventDeltaPosZVsCent"), partner.centrality, partnerPosZ - currentPosZ);
      hadHyperRegistry.fill(HIST("hMixEventDeltaCentFT0CVsCent"), partner.centrality, partner.centrality - currentEvent.centrality);
      for (std::size_t candidateIndex = 0; candidateIndex < currentEvent.candidates.size(); ++candidateIndex) {
        auto& candidate = currentEvent.candidates[candidateIndex];
        for (const auto& pion : partner.hadrons) {
          if (fillPair<isMC>(candidate, pion, currentEvent, partner, true, depth)) {
            ++candidate.mixedPairCount;
          }
        }
      }
      for (std::size_t candidateIndex = 0; candidateIndex < partner.candidates.size(); ++candidateIndex) {
        auto& candidate = partner.candidates[candidateIndex];
        for (const auto& pion : currentEvent.hadrons) {
          if (fillPair<isMC>(candidate, pion, partner, currentEvent, true, depth)) {
            ++candidate.mixedPairCount;
          }
        }
      }
    }
  }

  void flushMixedEventMultiplicity(const HadHyperEvent& event)
  {
    // One entry per candidate, including zero partners, after all its actual ME pairings.
    for (const auto& candidate : event.candidates) {
      hadHyperRegistry.fill(HIST("hCandidatePairMultiplicityME"), candidate.mixedPairCount, event.centrality);
    }
  }

  void flushMixingPools(std::unordered_map<int, std::deque<HadHyperEvent>>& pools)
  {
    for (const auto& [bin, pool] : pools) {
      for (const auto& event : pool) {
        flushMixedEventMultiplicity(event);
      }
    }
    pools.clear();
  }

  void storeEventInPool(std::deque<HadHyperEvent>& pool, HadHyperEvent&& event)
  {
    const int requestedMixingDepth = eventMixing.settingNoMixedEvents.value;
    if (requestedMixingDepth <= 0) {
      return;
    }
    const auto mixingDepth = static_cast<std::size_t>(requestedMixingDepth);
    if (pool.size() >= mixingDepth) {
      flushMixedEventMultiplicity(pool.front());
      pool.pop_front();
    }
    pool.push_back(std::move(event));
  }

  void endOfStream(o2::framework::EndOfStreamContext&)
  {
    flushMixingPools(mMixingPools);
    flushMixingPools(mMCMixingPools);
  }

  void processHyper(const HadHyperCollisionsFull& collisions, const TrackCandidates& tracks,
                    const HyperCandidates& candidates, const aod::BCsWithTimestamps& bcs)
  {
    processPairs</*isMC*/ false>(collisions, tracks, candidates, bcs, mMixingPools, mMixingRunNumber);
  }
  PROCESS_SWITCH(PiHypertritonFemto, processHyper, "Process pion-hypertriton same-event and mixed-event pairs", false);

  void processMCHyper(const HadHyperCollisionsFullMC& collisions, const TrackCandidatesMC& tracks,
                      const HyperCandidatesMC& candidates, const aod::McParticles&, const aod::BCsWithTimestamps& bcs)
  {
    processPairs</*isMC*/ true>(collisions, tracks, candidates, bcs, mMCMixingPools, mMCMixingRunNumber);
  }
  PROCESS_SWITCH(PiHypertritonFemto, processMCHyper, "Process MC pion-hypertriton same-event and mixed-event pairs", false);
};

void PiHypertritonFemto::init(o2::framework::InitContext&)
{
  *mEndOfStreamCallback = [this](EndOfStreamContext& context) { endOfStream(context); };
  hadHyperRegistry.addClone("MC/SE/", "MC/ME/");
  const AxisSpec pairMultiplicityAxis{axisPairMultiplicity, "pairs"};
  hadHyperRegistry.add("hCandidatePairMultiplicitySE", "Accepted same-event pion pairs per hypertriton candidate;pairs;centrality",
                       HistType::kTH2F, {pairMultiplicityAxis, {10, 0., 100.}});
  hadHyperRegistry.add("hCandidatePairMultiplicityME", "Total mixed-event pion pairs per hypertriton candidate over its pool lifetime;pairs;centrality",
                       HistType::kTH2F, {pairMultiplicityAxis, {10, 0., 100.}});

  const std::array<std::string, 2> eventNames{"SE", "ME"};
  const std::array<std::string, 2> daughterNames{"He3", "DecayPi"};
  for (std::size_t eventIndex = 0; eventIndex < eventNames.size(); ++eventIndex) {
    for (std::size_t daughterIndex = 0; daughterIndex < daughterNames.size(); ++daughterIndex) {
      const std::string prefix = "CPR/" + eventNames[eventIndex] + "/" + daughterNames[daughterIndex];
      auto& qa = mClosePairQA[eventIndex][daughterIndex];
      qa.before = hadHyperRegistry.add<TH2>((prefix + "/hBefore").c_str(), "All pairs;#Delta#eta;#LT#Delta#varphi*#GT",
                                            HistType::kTH2F, {{200, -0.1, 0.1}, {200, -0.1, 0.1}});
      qa.after = hadHyperRegistry.add<TH2>((prefix + "/hAfter").c_str(), "Pairs with offline CPR flags zero;#Delta#eta;#LT#Delta#varphi*#GT",
                                           HistType::kTH2F, {{200, -0.1, 0.1}, {200, -0.1, 0.1}});
      qa.validRadii = hadHyperRegistry.add<TH1>((prefix + "/hValidRadii").c_str(), "Common valid TPC radii (0=undefined CPR);radii;Comparisons",
                                                HistType::kTH1F, {{10, -0.5, 9.5}});
      qa.distanceBefore = hadHyperRegistry.add<TH1>((prefix + "/hDistanceBefore").c_str(), "All pairs;Average TPC separation (cm);Comparisons",
                                                    HistType::kTH1F, {{500, 0., 100.}});
      qa.distanceAfter = hadHyperRegistry.add<TH1>((prefix + "/hDistanceAfter").c_str(), "Pairs with offline CPR flags zero;Average TPC separation (cm);Comparisons",
                                                   HistType::kTH1F, {{500, 0., 100.}});
      qa.validPositionRadii = hadHyperRegistry.add<TH1>((prefix + "/hValidPositionRadii").c_str(), "Common extrapolated TPC radii (0=undefined distance);radii;Comparisons",
                                                        HistType::kTH1F, {{10, -0.5, 9.5}});
      qa.decisions = hadHyperRegistry.add<TH2>((prefix + "/hDecisions").c_str(), "Valid angular and distance comparisons;Angular flag;Distance flag",
                                               HistType::kTH2F, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    }
  }

  const bool processHyperPairs = doprocessHyper || doprocessMCHyper;
  if (processHyperPairs && hadHyper.maxOutputKstar.value == 0.f) {
    LOG(fatal) << "Pion-hypertriton mode requires a nonzero output k* range";
  }
  if (processHyperPairs && hadHyper.enableMixing.value && eventMixing.settingNoMixedEvents.value <= 0) {
    LOG(fatal) << "Pion-hypertriton mixed-event mode requires a positive mixing depth";
  }
  if (processHyperPairs && hypertriton.settingHypMassMin.value >= hypertriton.settingHypMassMax.value) {
    LOG(fatal) << "Pion-hypertriton mode requires settingHypMassMin < settingHypMassMax";
  }
  if (!std::isfinite(CPR.settingClosePairDeltaEtaMax.value) || !std::isfinite(CPR.settingClosePairDeltaPhiMax.value) ||
      CPR.settingClosePairDeltaEtaMax.value <= 0.f || CPR.settingClosePairDeltaPhiMax.value <= 0.f) {
    LOG(fatal) << "Close-pair rejection requires positive delta-eta and delta-phi limits";
  }
  if (!std::isfinite(CPR.settingClosePairDistanceMax.value) || CPR.settingClosePairDistanceMax.value <= 0.f) {
    LOG(fatal) << "Close-pair tagging requires a finite positive distance threshold";
  }
  mZorroSummary.setObject(mZorro.getZorroSummary());
  mRunNumber = -1;
  mCcdb->setURL(ccdb.settingCcdburl);
  mCcdb->setCaching(true);
  mCcdb->setLocalObjectValidityChecking();
  mCcdb->setFatalWhenNull(false);

  const std::array<std::string, 11> eventsLabels = {
    "All",
    "kNoITSROFrameBorder (data, optional)",
    "kIsTriggerTVX",
    "kNoTimeFrameBorder",
    "z_{vtx}",
    "kNoSameBunchPileup (optional)",
    "kIsGoodZvtxFT0vsPV (optional)",
    "sel8 (MC, optional)",
    "Zorro (data, optional)",
    "MC collision label (MC, optional)",
    "Mixing bin (if enabled)"};
  for (std::size_t i = 0; i < eventsLabels.size(); i++) {
    mQaRegistry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(i + 1, eventsLabels[i].c_str());
  }
}

template <bool isMC>
void PiHypertritonFemto::initCCDB(const aod::BCsWithTimestamps::iterator& bc)
{
  if (mRunNumber == bc.runNumber()) {
    return;
  }
  auto* magneticField = mCcdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdb.settingGrpmagPath.value, bc.timestamp());
  if (!magneticField) {
    LOG(fatal) << "Cannot load GRPMagField for daughter close-pair tagging at timestamp " << bc.timestamp();
    return;
  }
  constexpr float KilogaussToTesla = 0.1f;
  mMagneticFieldTesla = KilogaussToTesla * magneticField->getNominalL3Field();
  if constexpr (!isMC) {
    if (zorro.settingSkimmedProcessing.value) {
      mZorro.initCCDB(mCcdb.service, bc.runNumber(), bc.timestamp(), zorro.settingTriggerMask.value.c_str());
      mZorro.populateHistRegistry(mQaRegistry, bc.runNumber());
    }
  }
  mRunNumber = bc.runNumber();
}

template <bool isMC, typename Tcollision>
bool PiHypertritonFemto::selectCollision(const Tcollision& collision, const aod::BCsWithTimestamps&)
{
  auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
  initCCDB<isMC>(bc);
  mQaRegistry.fill(HIST("hEvents"), 0);

  // Match hyperRecoTask: the ITS readout-frame border cut applies only to data.
  if constexpr (!isMC) {
    if (!eventSelection.disableITSROFCut.value && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
  }
  // Disabled or inapplicable cuts pass through, keeping the cut-flow cumulative.
  mQaRegistry.fill(HIST("hEvents"), 1);
  if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
    return false;
  }
  mQaRegistry.fill(HIST("hEvents"), 2);
  if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
    return false;
  }
  mQaRegistry.fill(HIST("hEvents"), 3);
  if (std::abs(collision.posZ()) > eventMixing.settingCutVertex.value) {
    return false;
  }
  mQaRegistry.fill(HIST("hEvents"), 4);
  if (eventSelection.cfgEvSelkNoSameBunchPileup.value && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
    return false;
  }
  mQaRegistry.fill(HIST("hEvents"), 5);
  if (eventSelection.cfgEvSelkIsGoodZvtxFT0vsPV.value && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
    return false;
  }
  mQaRegistry.fill(HIST("hEvents"), 6);
  if constexpr (isMC) {
    if (mc.settingRequireSel8.value && !collision.sel8()) {
      return false;
    }
  }

  mQaRegistry.fill(HIST("hEvents"), 7);
  if constexpr (!isMC) {
    if (zorro.settingSkimmedProcessing.value) {
      if (!mZorro.isSelected(bc.globalBC())) {
        return false;
      }
    }
  }
  mQaRegistry.fill(HIST("hEvents"), 8);
  mQaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
  mQaRegistry.fill(HIST("hVtxZ"), collision.posZ());
  return true;
}

template <typename Ttrack>
bool PiHypertritonFemto::selectPionTrack(const Ttrack& candidate) const
{
  const float absPt = std::abs(candidate.pt());
  if (std::abs(candidate.eta()) > pionTrack.settingCutEta.value || absPt < pionTrack.settingPtMin.value || absPt > pionTrack.settingPtMax.value || absPt <= 0.f) {
    return false;
  }
  if (candidate.itsNClsInnerBarrel() < pionTrack.settingITSInnerBarrelMin.value ||
      candidate.itsNCls() < pionTrack.settingITSNClsMin.value ||
      candidate.tpcNClsFound() < pionTrack.settingTPCNClsFoundMin.value ||
      candidate.tpcNClsCrossedRows() < pionTrack.settingTPCCrossedRowsMin.value) {
    return false;
  }
  const float pionDCAxyMax = pionTrack.settingDCAxyOffset.value + pionTrack.settingDCAxyPtCoeff.value / absPt;
  const float pionDCAzMax = pionTrack.settingDCAzOffset.value + pionTrack.settingDCAzPtCoeff.value / absPt;
  return !(std::abs(candidate.dcaXY()) > pionDCAxyMax || std::abs(candidate.dcaZ()) > pionDCAzMax);
}

template <typename Ttrack>
bool PiHypertritonFemto::selectionPIDPion(const Ttrack& candidate)
{
  const float tpcNSigmaPi = candidate.tpcNSigmaPi();
  const float absP = std::abs(candidate.p());
  mQaRegistry.fill(HIST("hPionTPCNSigmaPreselection"), candidate.sign() * candidate.tpcInnerParam(), tpcNSigmaPi);
  if (absP <= pionPid.settingMomCombMin.value) {
    if (std::abs(tpcNSigmaPi) > pionPid.settingTPCNsigMax.value) {
      return false;
    }
    mQaRegistry.fill(HIST("hPionTPCNSigma"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
    mQaRegistry.fill(HIST("hPionTPC"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    return true;
  }
  if (!candidate.hasTOF()) {
    return false;
  }
  const float tofNSigmaPi = candidate.tofNSigmaPi();
  const float combNsigma = std::sqrt(tofNSigmaPi * tofNSigmaPi + tpcNSigmaPi * tpcNSigmaPi);
  mQaRegistry.fill(HIST("hPionTOFNSigmaPreselection"), candidate.sign() * candidate.pt(), tofNSigmaPi);
  mQaRegistry.fill(HIST("hPionCombNSigmaPreselection"), candidate.sign() * candidate.pt(), combNsigma);
  if (combNsigma > pionPid.settingCombNsigMax.value) {
    return false;
  }
  if (pionPid.settingReqSingleNsig.value && (std::abs(tpcNSigmaPi) > pionPid.settingCombNsigMax.value || std::abs(tofNSigmaPi) > pionPid.settingCombNsigMax.value)) {
    return false;
  }
  mQaRegistry.fill(HIST("hPionTPCNSigma"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
  mQaRegistry.fill(HIST("hPionTOFNSigma"), candidate.sign() * candidate.pt(), tofNSigmaPi);
  mQaRegistry.fill(HIST("hPionCombNSigma"), candidate.sign() * candidate.pt(), combNsigma);
  mQaRegistry.fill(HIST("hPionTPC"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
  return true;
}

template <typename Tcandidate>
bool PiHypertritonFemto::selectHyperCandidate(const Tcandidate& candidate)
{
  mQaRegistry.fill(HIST("hHypHe3TPCMomPreselection"), candidate.tpcMomHe());
  const float mass = computeHyperCandidateMass(candidate);
  if (!std::isfinite(mass) || mass < hypertriton.settingHypMassMin.value || mass > hypertriton.settingHypMassMax.value) {
    return false;
  }
  mQaRegistry.fill(HIST("hHypHe3TPCMom"), candidate.tpcMomHe());
  mQaRegistry.fill(HIST("hHypHe3TPCNSigma"), candidate.ptHe3(), candidate.nSigmaHe());
  return true;
}

template <bool isMC, typename Tcollisions, typename Ttracks, typename Tcandidates>
void PiHypertritonFemto::processPairs(const Tcollisions& collisions, const Ttracks& tracks, const Tcandidates& candidates,
                                      const aod::BCsWithTimestamps& bcs,
                                      std::unordered_map<int, std::deque<HadHyperEvent>>& mixingPools,
                                      int& mixingRunNumber)
{
  // Each invocation receives a new input-table scope; row indices can repeat in later frames.
  const uint64_t sourceFrameId = ++mSourceFrameId;
  const BinningType configuredBinningPolicy{{axisVertex, axisCentrality}, true};
  for (const auto& collision : collisions) {
    if (!selectCollision<isMC>(collision, bcs)) {
      continue;
    }
    if constexpr (isMC) {
      if (mc.settingRequireRecoMCCollisionMatch.value && !collision.has_mcCollision()) {
        continue;
      }
    }
    mQaRegistry.fill(HIST("hEvents"), 9);
    hadHyperRegistry.fill(HIST("hPoolFlow"), 0);
    int poolBin = -1;
    if (hadHyper.enableMixing.value) {
      poolBin = configuredBinningPolicy.getBin(std::make_tuple(collision.posZ(), collision.centFT0C()));
      if (poolBin < 0) {
        continue;
      }
      hadHyperRegistry.fill(HIST("hPoolFlow"), 1);
    }
    mQaRegistry.fill(HIST("hEvents"), 10);
    auto event = buildEvent<isMC>(collision, tracks, candidates);
    event.sourceFrameId = sourceFrameId;
    fillSameEventPairs<isMC>(event);
    if (!hadHyper.enableMixing.value) {
      continue;
    }
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (mixingRunNumber != bc.runNumber()) {
      flushMixingPools(mixingPools);
      mixingRunNumber = bc.runNumber();
    }
    auto& pool = mixingPools[poolBin];
    fillMixedEventPairs<isMC>(event, pool);
    storeEventInPool(pool, std::move(event));
  }
}

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  auto finishMixing = std::make_shared<std::function<void(EndOfStreamContext&)>>();
  auto task = adaptAnalysisTask<PiHypertritonFemto>(cfgc, finishMixing);
  auto initializeTask = task.algorithm.onInit;
  task.algorithm.onInit = [initializeTask, finishMixing](InitContext& context) {
    // CallbackService preserves registration order. Register before adaptAnalysisTask's
    // initializer so the pools are flushed BEFORE its histogram snapshot/cleanup callback.
    context.services().get<CallbackService>().set<CallbackService::Id::EndOfStream>(
      [finishMixing](EndOfStreamContext& eosContext) {
        if (*finishMixing) {
          (*finishMixing)(eosContext);
        }
      });
    return initializeTask(context);
  };
  return WorkflowSpec{std::move(task)};
}
