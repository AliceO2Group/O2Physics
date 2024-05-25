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

/// \file tagAndProbeDmesons.cxx
/// \brief Task for tracking efficiency studies with tag-and-probe using 3-prong D-meson decays
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace tagandprobe
{
enum TagChannels : uint8_t {
  DplusToKPiPi = 0,
  DsOrDplusToKKPi,
  DstarPlusToDzeroPi,
  DstarMinusToDzeroBarPi,
  DstarToDzeroToKK,
  NTagChannels
};

enum TrackTypes : uint8_t {
  GlobalWoDca = 0,
  GlobalWoDcaWoItsIb,
  GlobalWoDcaWoIts,
  GlobalWoDcaWoTpc,
  NTrackTypes
};

enum SignalFlags : uint8_t {
  Bkg = 0,
  Prompt,
  NonPrompt,
  Resonant
};

static constexpr int nBinsPt = 7;
static constexpr int nCutVars = 6;
static constexpr int nCutVarsDzero = 9;
constexpr float binsPt[nBinsPt + 1] = {0., 1., 2., 4., 6., 10., 20., 1000.};
auto vecBinsPt = std::vector<float>{binsPt, binsPt + nBinsPt + 1};

// default values for the cuts
constexpr float cuts[nBinsPt][nCutVars] = {{0.1f, 1.5f, 0.01f, 0.01f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.01f, 0.01f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.02f, 0.02f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.02f, 0.02f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.04f, 0.04f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.04f, 0.04f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.06f, 0.06f, 2.f, 2.f}};

constexpr float cutsDzero[nBinsPt][nCutVarsDzero] = {{1.815f, 1.915f, 0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                                     {1.815f, 1.915f, 0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                                     {1.815f, 1.915f, 0.06f, 0.06f, 2.f, 2.f, 0.f, 0.95f, 0.95f}};

static const std::vector<std::string> labelsPt{};
static const std::vector<std::string> labelsCutVar = {"minMass", "maxMass", "decayLength", "decayLengthXY", "normDecayLength", "normDecayLengthXY"};
static const std::vector<std::string> labelsCutVarDzero = {"minMass", "maxMass", "decayLength", "decayLengthXY", "normDecayLength", "normDecayLengthXY", "impParProd", "cosPointing", "cosPointingXY"};

DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, Tracks, "_0"); //! Index to first track
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, Tracks, "_1"); //! Index to second track
// Topological variables
DECLARE_SOA_COLUMN(TagPt, tagPt, float);                                     //! Tag's pT
DECLARE_SOA_COLUMN(TagInvMass, tagInvMass2, float);                          //! Tag's invMass
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                         //! Decay length of the tag (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                     //! Transverse decay length of the tag (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);     //! Normalised decay length of the tag
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float); //! Normalised transverse decay length of the tag
DECLARE_SOA_COLUMN(TrackDcaXY0, trackDcaXY0, float);                         //! DCAxy of the first tag track
DECLARE_SOA_COLUMN(TrackDcaXY1, trackDcaXY1, float);                         //! DCAxy of the second tag track
DECLARE_SOA_COLUMN(ProductTrackDcaXY, productTrackDcaXY, float);             //! Product of DCAxy of the two tag tracks
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                         //! Cosine pointing angle of the tag
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                     //! Cosine of the pointing angle in XY of the tag
DECLARE_SOA_COLUMN(IsSignal, isSignal, uint8_t);                             //! Flag for a signal
DECLARE_SOA_COLUMN(DecChannel, decChannel, uint8_t);                         //! Flag the selected decay channel
} // namespace tagandprobe

DECLARE_SOA_TABLE(PiPiFromDpTags, "AOD", "PIPIFROMDPTAG", //! Table for same sign 2-pion vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id);
DECLARE_SOA_TABLE(KaKaFromDspTags, "AOD", "KAKAFROMDSPTAG", //! Table for opposite sign 2-kaon vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(PiKaFromDzTags, "AOD", "PIKAFROMDZTAG", //! Table for opposite sign pion(+)-kaon(-) vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  soa::Marker<2>);
DECLARE_SOA_TABLE(KaPiFromDzTags, "AOD", "KAPIFROMDZTAG", //! Table for opposite sign kaon(+)-pion(-) vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  soa::Marker<3>);
DECLARE_SOA_TABLE(TagTopoVariables, "AOD", "TAGTOPOVARIABLE", //! Table for the Tags' Topological variables
                  aod::tagandprobe::TagPt,
                  aod::tagandprobe::TagInvMass,
                  aod::tagandprobe::DecayLength,
                  aod::tagandprobe::DecayLengthXY,
                  aod::tagandprobe::DecayLengthNormalised,
                  aod::tagandprobe::DecayLengthXYNormalised,
                  aod::tagandprobe::TrackDcaXY0,
                  aod::tagandprobe::TrackDcaXY1,
                  aod::tagandprobe::ProductTrackDcaXY,
                  aod::tagandprobe::Cpa,
                  aod::tagandprobe::CpaXY,
                  aod::tagandprobe::IsSignal,
                  aod::tagandprobe::DecChannel);
} // namespace o2::aod

/// Reconstruction of 2-prong displaced vertices (very good quality and purity)
/// 1) K∓K± for φ from Ds± or D± → φπ± decays
/// 2) π±π± for D± → K∓π±π± decays
/// 3) K∓π± for D0 from D±* → D0π± decays
struct TagTwoProngDisplacedVertices {

  Produces<aod::PiPiFromDpTags> tagPiPiTable;
  Produces<aod::KaKaFromDspTags> tagKaKaTable;
  Produces<aod::KaPiFromDzTags> tagKaPiTable;
  Produces<aod::PiKaFromDzTags> tagPiKaTable;
  Produces<aod::TagTopoVariables> tagVarsTable;
  SliceCache cache;
  Configurable<int> fillTagTable{"fillTagTable", 0, "flag to fill tag table with topological variables (0 -> disabled, 1 -> signal only, 2 -> bkg only, 3 -> both)"};
  Configurable<bool> applyTofPid{"applyTofPid", true, "flag to enable TOF PID selection"};
  Configurable<bool> studyDzeroReflections{"studyDzeroReflections", false, "flag to study Dzero reflections"};
  Configurable<float> trackNumSigmaTof{"trackNumSigmaTof", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackNumSigmaTpc{"trackNumSigmaTpc", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackDcaXyMin{"trackDcaXyMin", 0.002f, "minimum DCAxy for tracks with pT < 2 GeV/c"};
  Configurable<float> trackPtMin{"trackPtMin", 0.4f, "minimum track pT"};

  Configurable<std::vector<float>> binsPtPiPiFromDplus{"binsPtPiPiFromDplus", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for pipi pairs from D+ decays"};
  Configurable<std::vector<float>> binsKaKaFromDsOrDplus{"binsKaKaFromDsOrDplus", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from Ds or D+ decays"};
  Configurable<std::vector<float>> binsPtDzeroFromDstar{"binsPtDzeroFromDstar", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for Kpi pairs from D0 <- D*+ decays"};
  Configurable<std::vector<float>> binsPtDzeroKaKaFromDstar{"binsPtDzeroKaKaFromDstar", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from D0 <- D*+ decays"};

  Configurable<LabeledArray<float>> cutsPiPiFromDplus{"cutsPiPiFromDplus", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVar}, "Selections for pipi pairs from D+ decays"};
  Configurable<LabeledArray<float>> cutsKaKaFromDsOrDplus{"cutsKaKaFromDsOrDplus", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVar}, "Selections for KK pairs from Ds or D+ decays"};
  Configurable<LabeledArray<float>> cutsDzeroFromDstar{"cutsDzeroFromDstar", {aod::tagandprobe::cutsDzero[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVarsDzero, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVarDzero}, "Selections for Kpi pairs from D0 <- D*+ decays"};
  Configurable<LabeledArray<float>> cutsDzeroKaKaFromDstar{"cutsDzeroKaKaFromDstar", {aod::tagandprobe::cutsDzero[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVarsDzero, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVarDzero}, "Selections for Kpi pairs from D0 <- D*+ decays"};

  using TracksWithSelAndDca = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using TracksWithSelAndDcaMc = soa::Join<TracksWithSelAndDca, aod::McTrackLabels>;
  using CollisionsWithEvSel = soa::Join<aod::Collisions, aod::EvSels>;

  Filter evSelFilter = aod::evsel::sel8 == true;                                                                                                              // simple event selection
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;                                                                                                 // simple event selection
  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && aod::track::pt > trackPtMin && (nabs(aod::track::dcaXY) > trackDcaXyMin || aod::track::pt > 2.f); // for the tag, we only consider global tracks with large dcaXY (low pT only)
  using TracksWithSelAndDcaFiltered = soa::Filtered<TracksWithSelAndDca>;
  using TracksWithSelAndDcaMcFiltered = soa::Filtered<TracksWithSelAndDcaMc>;
  using CollisionsFiltered = soa::Filtered<CollisionsWithEvSel>;

  // in the partition we only apply TPC PID
  Preslice<TracksWithSelAndDcaFiltered> perCollision = aod::track::collisionId;
  Partition<TracksWithSelAndDcaFiltered> positivePions = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativePions = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> positiveKaons = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativeKaons = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> positivePionsMc = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> negativePionsMc = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> positiveKaonsMc = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> negativeKaonsMc = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;

  ccdb::CcdbApi ccdbApi;
  Service<ccdb::BasicCCDBManager> ccdb;
  vertexing::DCAFitterN<2> vertexer;
  int runNumber{0};

  std::array<std::array<double, 2>, aod::tagandprobe::TagChannels::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged}};

  std::array<LabeledArray<float>, aod::tagandprobe::TagChannels::NTagChannels> topologicalCuts{};
  std::array<std::vector<float>, aod::tagandprobe::TagChannels::NTagChannels> ptBinsForTopologicalCuts{};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::string ccdbUrl = "http://alice-ccdb.cern.ch";
    ccdb->setURL(ccdbUrl.data());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(ccdbUrl.data());

    vertexer.setPropagateToPCA(true);
    vertexer.setMaxR(200.f);
    vertexer.setMaxDZIni(4.f);
    vertexer.setMinParamChange(1.e-3);
    vertexer.setMinRelChi2Change(0.9f);
    vertexer.setUseAbsDCA(false);

    topologicalCuts = {cutsPiPiFromDplus, cutsKaKaFromDsOrDplus, cutsDzeroFromDstar, cutsDzeroFromDstar, cutsDzeroKaKaFromDstar};
    ptBinsForTopologicalCuts = {binsPtPiPiFromDplus, binsKaKaFromDsOrDplus, binsPtDzeroFromDstar, binsPtDzeroFromDstar, binsPtDzeroKaKaFromDstar};

    const AxisSpec axisPt{250, 0.f, 50.f};
    const AxisSpec axisPtDzeroRefl{{0.f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f, 24.0f, 30.0f, 40.0f}};
    const AxisSpec axisMassPiPi{250, 0.f, 2.5f};
    const AxisSpec axisReflFlag{3, 0.5f, 3.5f};
    const AxisSpec axisMassKaKa{200, constants::physics::MassPhi - 0.05f, constants::physics::MassPhi + 0.05f};
    const AxisSpec axisMassKaPi{400, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f};

    if (doprocessPiPiFromDplus) {
      registry.add<TH2>("hMassPiPiVsPt", ";#it{p}_{T}(#pi#pi) (GeV/#it{c}); #it{M}(#pi#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassPiPi});
    }
    if (doprocessKaKaFromDsOrDplus) {
      registry.add<TH2>("hMassKaKaVsPt", ";#it{p}_{T}(KK) (GeV/#it{c}); #it{M}(KK) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaKa});
    }
    if (doprocessKaPiFromDstar) {
      if (!studyDzeroReflections) {
        registry.add<TH2>("hMassKaPiVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaPi});
      } else {
        registry.add<THnSparse>("hMassKaPiVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2}); #it{M}(#piK) (GeV/#it{c}^{2}); ReflFag", HistType::kTHnSparseF, {axisPtDzeroRefl, axisMassKaPi, axisMassKaPi, axisReflFlag});
      }
    }
    if (doprocessKaKaFromDzero) {
      registry.add<TH2>("hMassDzeroKaKaVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaPi});
    }
  }

  /// Finds pT bin in an array.
  /// \param bins  array of pT bins
  /// \param value  pT
  /// \return index of the pT bin
  template <typename T1, typename T2>
  int findBin(T1 const& binsPt, T2 value)
  {
    if (value < binsPt->front()) {
      return -1;
    }
    if (value >= binsPt->back()) {
      return -1;
    }
    return std::distance(binsPt->begin(), std::upper_bound(binsPt->begin(), binsPt->end(), value)) - 1;
  }

  /// Fill a vector with the Mothers pdg codes
  /// \param pdgMother vector with the pdg codes
  /// \param channel decay channel
  void pdgMothersDecayChannel(std::vector<int>& pdgDecayMothers, std::vector<int>& pdgResonances, const uint8_t channel)
  {
    pdgDecayMothers.clear();
    pdgResonances.clear();
    if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
      pdgDecayMothers.push_back(constants::physics::Pdg::kDPlus);
      pdgResonances.push_back(313);    // K*(892)0
      pdgResonances.push_back(10313);  // K1(1270)0
      pdgResonances.push_back(100313); // K*(1410)0
      pdgResonances.push_back(10311);  // K*0(1430)0
      pdgResonances.push_back(100311); // K*(1460)0
      pdgResonances.push_back(20313);  // K1(1400)0
      pdgResonances.push_back(30313);  // K*(1680)0
    } else if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
      pdgDecayMothers.push_back(constants::physics::Pdg::kDPlus);
      pdgDecayMothers.push_back(constants::physics::Pdg::kDS);
    } else if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi || channel == aod::tagandprobe::TagChannels::DstarToDzeroToKK) {
      pdgDecayMothers.push_back(constants::physics::Pdg::kDStar);
    }
  }

  /// Check if the given tag tracks belong to a D meson
  /// \param firstTrack candidate
  /// \param SecondTrack candidate
  /// \param particlesMc McParticles table
  /// \param channel decay channel
  /// \return a flag that contains the information of MC truth (see aod::tagandprobe::SignalFlags)
  template <typename PParticles, typename TTrack>
  uint8_t getTagOrigin(TTrack const& firsTrack,
                       TTrack const& secondTrack,
                       PParticles const& particlesMc,
                       const uint8_t channel,
                       std::vector<int>& pdgDecayMothers,
                       std::vector<int>& pdgResonances)
  {
    int pdgTagMother{0};
    int pdgProbeParticle{-1};
    uint8_t signalFlag = 0;

    if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
      pdgTagMother = constants::physics::Pdg::kDPlus;
      pdgProbeParticle = 321; // Ka
    } else if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
      pdgTagMother = constants::physics::Pdg::kPhi;
      pdgProbeParticle = 211; // Pi
    } else {
      pdgTagMother = constants::physics::Pdg::kD0;
      pdgProbeParticle = 211; // Pi
    }

    if (!firsTrack.has_mcParticle() || !secondTrack.has_mcParticle()) {
      return BIT(aod::tagandprobe::SignalFlags::Bkg);
    } else {
      auto firstMcTrack = firsTrack.template mcParticle_as<PParticles>();
      auto secondMcTrack = secondTrack.template mcParticle_as<PParticles>();
      auto firstTrackMotherId = RecoDecay::getMother(particlesMc, firstMcTrack, pdgTagMother, true);
      auto secondTrackMotherId = RecoDecay::getMother(particlesMc, secondMcTrack, pdgTagMother, true);

      bool isTaggedAsSignal{false}, isResonant{false};
      if ((firstTrackMotherId == secondTrackMotherId) && (firstTrackMotherId != -1)) {
        auto particleMother = particlesMc.rawIteratorAt(firstTrackMotherId);

        /// π±π± for D± → K∓π±π± decays
        if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
          auto particleMother = particlesMc.rawIteratorAt(firstTrackMotherId);
          auto daughters = particleMother.template daughters_as<PParticles>();

          // Check if the probe is within the mother's particle daughters
          if (daughters.size() == 3) { // non-resonant decay
            for (auto& daughter : daughters) {
              if (std::abs(daughter.pdgCode()) == pdgProbeParticle) {
                isTaggedAsSignal = true;
                break;
              }
            }
          } else if (daughters.size() == 2) { // resonant decay
            for (auto& daughter : daughters) {
              auto absPdg = std::abs(daughter.pdgCode());
              if (std::find(pdgResonances.begin(), pdgResonances.end(), absPdg) != pdgResonances.end()) {
                isTaggedAsSignal = true;
                isResonant = true;
                break;
              }
            }
          }
        } else {
          ///  K∓K± for φ from Ds± or D± → φπ± decays
          ///  K∓π± for D0 from D±* → D0π± decays
          for (auto pdgGrandMother : pdgDecayMothers) {
            auto grandMotherId = RecoDecay::getMother(particlesMc, particleMother, pdgGrandMother, true);
            if (grandMotherId != -1) {
              auto particleGrandMother = particlesMc.rawIteratorAt(grandMotherId);
              auto daughters = particleGrandMother.template daughters_as<PParticles>();
              // Check if the probe is within the GrandMother's particle daughters
              if (daughters.size() == 2) { // exclude undesired decays, such as Ds± → φπ±π±π∓
                for (auto& daughter : daughters) {
                  if (std::abs(daughter.pdgCode()) == pdgProbeParticle) {
                    isTaggedAsSignal = true;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      // check if it is non-prompt from beauty
      if (isTaggedAsSignal) {
        if (RecoDecay::getCharmHadronOrigin(particlesMc, particlesMc.rawIteratorAt(firstTrackMotherId)) == RecoDecay::OriginType::NonPrompt) {
          SETBIT(signalFlag, aod::tagandprobe::SignalFlags::NonPrompt);
        } else {
          SETBIT(signalFlag, aod::tagandprobe::SignalFlags::Prompt);
        }
        if (isResonant) {
          SETBIT(signalFlag, aod::tagandprobe::SignalFlags::Resonant);
        }
        return signalFlag;
      }

      return BIT(aod::tagandprobe::SignalFlags::Bkg);
    }
  }

  template <typename Pvec>
  bool isSelectedInvariantMass(const Pvec& pVecTrackFirst,
                               const Pvec& pVecTrackSecond,
                               const uint8_t channel,
                               float& invMass2,
                               const int& ptBin)
  {
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond};

    auto invMassMin = topologicalCuts[channel].get(ptBin, 0u);
    auto invMassMax = topologicalCuts[channel].get(ptBin, 1u);
    invMass2 = RecoDecay::m2(arrMomentum, masses[channel]);
    if (invMass2 > invMassMax * invMassMax || invMass2 < invMassMin * invMassMin) {
      return false;
    }
    return true;
  }

  /// Calculate all the topology variables and store them in the Topology table
  /// \param doMc 1 for the Mc and 0 for the data
  /// \param primVtx  primary vertex
  /// \param secVtx  secondary vertex
  /// \param trackDcaXy array with the Tags' TrackDCAXY
  /// \param channel decay channel
  /// \param firstTrack candidate
  /// \param SecondTrack candidate
  /// \param particlesMc McParticle table
  template <bool doMc, typename PV, typename SV, typename CovMatSV, typename PVec, typename TTrack, typename PParticles>
  void getTagInfo(const PV& primVtx,
                  const SV& secVtx,
                  const CovMatSV& covMatrixSecVtx,
                  const PVec& pVec,
                  std::array<float, 2>& trackDcaXy,
                  const uint8_t channel,
                  const TTrack& firstTrack,
                  const TTrack& secondTrack,
                  float& invMass2,
                  std::vector<int>& pdgDecayMothers,
                  std::vector<int>& pdgResonances,
                  const PParticles& particlesMc)
  {
    auto covMatrixPV = primVtx.getCov();
    float phi, theta;
    std::array<float, 3> pvCoord = {primVtx.getX(), primVtx.getY(), primVtx.getZ()};
    getPointDirection(pvCoord, secVtx, phi, theta);

    auto decLen = RecoDecay::distance(pvCoord, secVtx);
    auto errorDecLen = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, theta));
    auto errorDecLenXy = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.f) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, 0.f));
    auto decLenXy = RecoDecay::distanceXY(pvCoord, secVtx);
    auto cpa = RecoDecay::cpa(pvCoord, secVtx, pVec);
    auto cpaXy = RecoDecay::cpaXY(pvCoord, secVtx, pVec);
    auto normDecLen = decLen / errorDecLen;
    auto normDecLenXy = decLenXy / errorDecLenXy;
    auto tagsPt = RecoDecay::pt(pVec);
    auto invMass = std::sqrt(invMass2);

    uint8_t isSignal = 0; // default value for data (no bkg, no signal)

    bool fillTable{true};
    if constexpr (doMc) {
      isSignal = getTagOrigin(firstTrack, secondTrack, particlesMc, channel, pdgDecayMothers, pdgResonances);
      if (fillTagTable == 1 && !(TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Prompt) || TESTBIT(isSignal, aod::tagandprobe::SignalFlags::NonPrompt))) { // only signal
        fillTable = false;
      } else if (fillTagTable == 2 && !TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Bkg)) { // only background
        fillTable = false;
      }
    }
    if (fillTable) {
      tagVarsTable(tagsPt, invMass, decLen, decLenXy, normDecLen, normDecLenXy, trackDcaXy[0], trackDcaXy[1], trackDcaXy[0] * trackDcaXy[1], cpa, cpaXy, isSignal, channel);
    }
  }

  template <typename TTrack>
  bool isSelectedPidTof(const TTrack& track,
                        const uint8_t channel)
  {
    if (!track.hasTOF()) { // TOF not forced anyway
      return true;
    }

    switch (channel) {
      case aod::tagandprobe::TagChannels::DplusToKPiPi: {
        if (std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DsOrDplusToKKPi: {
        if (std::abs(track.tofNSigmaKa()) < trackNumSigmaTof) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DstarPlusToDzeroPi: {
        if ((track.signed1Pt() > 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() < 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi: {
        if ((track.signed1Pt() < 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() > 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DstarToDzeroToKK: {
        if (std::abs(track.tofNSigmaKa()) < trackNumSigmaTof) {
          return true;
        }
        break;
      }
    }
    return false;
  }

  template <typename PV, typename SV, typename CovMatSV, typename PVec>
  bool isSelectedTopology(const PV& primVtx,
                          const SV& secVtx,
                          const CovMatSV& covMatrixSecVtx,
                          const PVec& pVec,
                          std::array<float, 2>& trackDcaXy,
                          const uint8_t channel,
                          const int& ptBin)
  {
    std::array<float, 3> pvCoord = {primVtx.getX(), primVtx.getY(), primVtx.getZ()};
    auto decLen = RecoDecay::distance(pvCoord, secVtx);
    if (decLen < topologicalCuts[channel].get(ptBin, 2u)) {
      return false;
    }

    auto covMatrixPV = primVtx.getCov();

    auto decLenXy = RecoDecay::distanceXY(pvCoord, secVtx);
    if (decLenXy < topologicalCuts[channel].get(ptBin, 3u)) {
      return false;
    }

    float phi, theta;
    getPointDirection(pvCoord, secVtx, phi, theta);
    auto errorDecLen = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, theta));
    if (decLen / errorDecLen < topologicalCuts[channel].get(ptBin, 4u)) {
      return false;
    }

    auto errorDecLenXy = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.f) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, 0.f));
    if (decLenXy / errorDecLenXy < topologicalCuts[channel].get(ptBin, 5u)) {
      return false;
    }

    // only for D0 meson
    if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi || channel == aod::tagandprobe::TagChannels::DstarToDzeroToKK) {
      if (trackDcaXy[0] * trackDcaXy[1] > topologicalCuts[channel].get(ptBin, 6u)) {
        return false;
      }
      auto cpa = RecoDecay::cpa(pvCoord, secVtx, pVec);
      if (cpa < topologicalCuts[channel].get(ptBin, 7u)) {
        return false;
      }
      auto cpaXy = RecoDecay::cpaXY(pvCoord, secVtx, pVec);
      if (cpaXy < topologicalCuts[channel].get(ptBin, 8u)) {
        return false;
      }
    }

    return true;
  }

  template <bool doMc, typename CCollision, typename TTracks, typename PParticles>
  void computeCombinatorialSameCharge(CCollision const& collision,
                                      TTracks const& tracks, // pool of tracks
                                      const uint8_t channel,
                                      float& /*bz*/,
                                      std::vector<int>& pdgDecayMothers,
                                      std::vector<int>& pdgResonances,
                                      PParticles const& particlesMc)
  {
    for (auto trackFirst = tracks.begin(); trackFirst != tracks.end(); ++trackFirst) {

      if (applyTofPid && !isSelectedPidTof(trackFirst, channel)) {
        continue;
      }

      for (auto trackSecond = trackFirst + 1; trackSecond != tracks.end(); ++trackSecond) {

        if (applyTofPid && !isSelectedPidTof(trackSecond, channel)) {
          continue;
        }

        float invMass2{0.f};
        std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
        std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};

        auto pVec = RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond);
        auto ptBin = findBin(&ptBinsForTopologicalCuts[channel], RecoDecay::pt(pVec));
        if (ptBin == -1) {
          continue;
        }

        if (!isSelectedInvariantMass(pVecTrackFirst, pVecTrackSecond, channel, invMass2, ptBin)) {
          continue;
        }

        auto trackParCovFirst = getTrackParCov(trackFirst);
        auto trackParCovSecond = getTrackParCov(trackSecond);

        int nVertices{0};
        try {
          nVertices = vertexer.process(trackParCovFirst, trackParCovSecond);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices != 1) {
          continue;
        }

        auto primVtx = getPrimaryVertex(collision);
        const auto& secVtx = vertexer.getPCACandidate();
        const auto& covMatrixPCA = vertexer.calcPCACovMatrixFlat();
        std::array<float, 2> trackDcaXy{trackFirst.dcaXY(), trackSecond.dcaXY()};
        if (fillTagTable) {
          getTagInfo<doMc>(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel, trackFirst, trackSecond, invMass2, pdgDecayMothers, pdgResonances, particlesMc);
        } else {
          if (!isSelectedTopology(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel, ptBin)) {
            continue;
          }
          registry.fill(HIST("hMassPiPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2)); // only channel with same sign tracks for the moment
          tagPiPiTable(trackFirst.collisionId(), trackFirst.globalIndex(), trackSecond.globalIndex());
        }
      }
    }
  }

  template <bool doMc, typename CCollision, typename TTracks, typename PParticles>
  void computeCombinatorialOppositeCharge(CCollision const& collision,
                                          TTracks const& tracksPos,
                                          TTracks const& tracksNeg,
                                          const uint8_t channel,
                                          float& /*bz*/,
                                          std::vector<int>& pdgDecayMothers,
                                          std::vector<int>& pdgResonances,
                                          PParticles const& particlesMc)
  {
    for (const auto& trackPos : tracksPos) {

      if (applyTofPid && !isSelectedPidTof(trackPos, channel)) {
        continue;
      }

      for (const auto& trackNeg : tracksNeg) {

        if (applyTofPid && !isSelectedPidTof(trackNeg, channel)) {
          continue;
        }

        float invMass2{0.f};
        std::array<float, 3> pVecTrackPos{trackPos.px(), trackPos.py(), trackPos.pz()};
        std::array<float, 3> pVecTrackNeg{trackNeg.px(), trackNeg.py(), trackNeg.pz()};

        auto pVec = RecoDecay::pVec(pVecTrackPos, pVecTrackNeg);
        auto ptBin = findBin(&ptBinsForTopologicalCuts[channel], RecoDecay::pt(pVec));
        if (ptBin == -1) {
          continue;
        }

        if (!isSelectedInvariantMass(pVecTrackPos, pVecTrackNeg, channel, invMass2, ptBin)) {
          continue;
        }

        auto trackParCovPos = getTrackParCov(trackPos);
        auto trackParCovNeg = getTrackParCov(trackNeg);

        int nVertices{0};
        try {
          nVertices = vertexer.process(trackParCovPos, trackParCovNeg);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices != 1) {
          continue;
        }

        auto primVtx = getPrimaryVertex(collision);
        const auto& secVtx = vertexer.getPCACandidate();
        const auto& covMatrixPCA = vertexer.calcPCACovMatrixFlat();
        std::array<float, 2> trackDcaXy{trackPos.dcaXY(), trackNeg.dcaXY()};
        if (fillTagTable) {
          getTagInfo<doMc>(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel, trackPos, trackNeg, invMass2, pdgDecayMothers, pdgResonances, particlesMc);
        } else {
          if (!isSelectedTopology(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel, ptBin)) {
            continue;
          }
          if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
            registry.fill(HIST("hMassKaKaVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
            tagKaKaTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex());
          } else if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) {
            if (!studyDzeroReflections) {
              registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
            } else {
              float invMassrefl{0.f};
              int isDzero = 1;
              if (std::abs(trackPos.tpcNSigmaKa()) < trackNumSigmaTpc && (std::abs(trackNeg.tpcNSigmaPi()) < trackNumSigmaTpc)) {
                isDzero = 3;
                if (applyTofPid) {
                  if (!isSelectedPidTof(trackNeg, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) || !isSelectedPidTof(trackPos, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi))
                    isDzero = 1;
                }
              }
              if (isDzero == 3) {
                auto arrMomentum = std::array{pVecTrackNeg, pVecTrackPos};
                invMassrefl = std::sqrt(RecoDecay::m2(arrMomentum, masses[channel]));
              }
              registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2), invMassrefl, isDzero);
            }
            tagPiKaTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex());
          } else if (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
            if (!studyDzeroReflections) {
              registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
            } else {
              float invMassrefl{0.f};
              int isDzero = 2;
              if (std::abs(trackNeg.tpcNSigmaKa()) < trackNumSigmaTpc && (std::abs(trackPos.tpcNSigmaPi()) < trackNumSigmaTpc)) {
                isDzero = 3;
                if (applyTofPid) {
                  if (!isSelectedPidTof(trackNeg, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) || !isSelectedPidTof(trackPos, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi))
                    isDzero = 2;
                }
              }
              if (isDzero == 3) {
                auto arrMomentum = std::array{pVecTrackNeg, pVecTrackPos};
                invMassrefl = std::sqrt(RecoDecay::m2(arrMomentum, masses[channel]));
              }
              registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2), invMassrefl, isDzero);
            }
            tagKaPiTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex());
          } else if (channel == aod::tagandprobe::TagChannels::DstarToDzeroToKK) {
            registry.fill(HIST("hMassDzeroKaKaVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
          }
        }
      }
    }
  }

  void processPiPiFromDplusMc(CollisionsFiltered::iterator const& collision,
                              TracksWithSelAndDcaMcFiltered const&,
                              aod::BCsWithTimestamps const&,
                              aod::McParticles const& particlesMc)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DplusToKPiPi);

    auto groupPositive = positivePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<true>(collision, groupPositive, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, particlesMc);

    auto groupNegative = negativePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<true>(collision, groupNegative, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, particlesMc);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processPiPiFromDplusMc, "Process pipi combinatorial to tag pion pairs from D+ decays Mc", false);

  void processPiPiFromDplus(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DplusToKPiPi);

    auto groupPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<false>(collision, groupPositive, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, tracks);

    auto groupNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<false>(collision, groupNegative, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processPiPiFromDplus, "Process pipi combinatorial to tag pion pairs from D+ decays", false);

  void processKaKaFromDsOrDplusMc(CollisionsFiltered::iterator const& collision,
                                  TracksWithSelAndDcaMcFiltered const&,
                                  aod::BCsWithTimestamps const&,
                                  aod::McParticles const& particlesMc)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DsOrDplusToKKPi);

    auto groupPositive = positiveKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<true>(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DsOrDplusToKKPi, bz, pdgDecayMothers, pdgResonances, particlesMc);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDsOrDplusMc, "Process KK combinatorial to tag kaon pairs from Ds+/D+ decays Mc", false);

  void processKaKaFromDsOrDplus(CollisionsFiltered::iterator const& collision,
                                TracksWithSelAndDcaFiltered const& tracks,
                                aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DsOrDplusToKKPi);

    auto groupPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<false>(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DsOrDplusToKKPi, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDsOrDplus, "Process KK combinatorial to tag kaon pairs from Ds+/D+ decays", false);

  void processKaKaFromDzero(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DstarToDzeroToKK);

    auto groupPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<false>(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DstarToDzeroToKK, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDzero, "Process KK combinatorial to tag kaon pairs from Dzero decays", false);

  void processKaPiFromDstar(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi);

    auto groupPionPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupPionNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<false>(collision, groupPionPositive, groupKaonNegative, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, bz, pdgDecayMothers, pdgResonances, tracks);
    computeCombinatorialOppositeCharge<false>(collision, groupKaonPositive, groupPionNegative, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaPiFromDstar, "Process Kpi combinatorial to tag D0 from D*+ decays", false);

  void processKaPiFromDstarMc(CollisionsFiltered::iterator const& collision,
                              TracksWithSelAndDcaMcFiltered const&,
                              aod::BCsWithTimestamps const&,
                              aod::McParticles const& particlesMc)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi);

    auto groupPionPositive = positivePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupPionNegative = negativePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonPositive = positiveKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonNegative = negativeKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<true>(collision, groupPionPositive, groupKaonNegative, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, bz, pdgDecayMothers, pdgResonances, particlesMc);
    computeCombinatorialOppositeCharge<true>(collision, groupKaonPositive, groupPionNegative, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, bz, pdgDecayMothers, pdgResonances, particlesMc);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaPiFromDstarMc, "Process Kpi combinatorial to tag D0 from D*+ decays", false);
};

/// Probe third track reconstruction efficiency with different selections
struct ProbeThirdTrack {

  using TracksWithDca = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>;

  Preslice<aod::PiPiFromDpTags> tagsPiPiPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::KaKaFromDspTags> tagsKaKaPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::PiKaFromDzTags> tagsPiKaPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::KaPiFromDzTags> tagsKaPiPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  std::array<std::array<double, 3>, aod::tagandprobe::TagChannels::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged}};

  std::array<TrackSelection, aod::tagandprobe::TrackTypes::NTrackTypes> trackSelector{}; // define the track selectors

  std::array<std::array<std::shared_ptr<THnSparse>, aod::tagandprobe::TrackTypes::NTrackTypes>, aod::tagandprobe::TagChannels::NTagChannels> histos{};
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    // ITS-TPC tracks (global tracks)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireITSRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMinNCrossedRowsTPC(50);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxChi2PerClusterTPC(10.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any IB layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxDcaZ(2.f);

    // TPC tracks (global tracks without ITS IB requirement)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMinNCrossedRowsTPC(50);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMaxChi2PerClusterTPC(10.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetRequireHitsInITSLayers(3, {0, 1, 2, 3, 4, 5, 6}); // at least three hits in whatever layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMaxDcaZ(2.f);

    // TPC tracks (global tracks without ITS requirements)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMinNCrossedRowsTPC(50);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMaxChi2PerClusterTPC(10.f);

    // ITS tracks (global tracks without TPC requirements)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetRequireITSRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any SPD layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetMaxDcaZ(2.f);

    const AxisSpec axisPtProbe{{0.05f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.f, 12.f, 15.f, 20.f, 25.f, 30.f}};
    const AxisSpec axisPtTag{{0.05f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.f, 12.f, 15.f, 20.f, 25.f, 30.f}};
    const AxisSpec axisPtD{{0.f, 0.5f, 1.f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 5.5f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 8.5f, 9.0f, 9.5f, 10.f, 11.f, 12.f, 14.f, 16.f, 20.f, 24.f, 36.f, 50.f}};
    const AxisSpec axisEtaProbe{20, -1.f, 1.f};
    const AxisSpec axisNumCrossRowTpc{51, 49.5f, 100.5f};
    const AxisSpec axisTpcChi2PerClus{8, 2.f, 10.f};
    const AxisSpec axisNumCluIts{5, 2.5f, 7.5f};
    std::array<AxisSpec, aod::tagandprobe::TagChannels::NTagChannels> axisMass = {AxisSpec{225, 1.65f, 2.10f}, AxisSpec{225, 1.65f, 2.10f}, AxisSpec{350, 0.135f, 0.17f}, AxisSpec{350, 0.135f, 0.17f}, AxisSpec{350, 0.135f, 0.17f}};
    std::array<AxisSpec, aod::tagandprobe::TagChannels::NTagChannels> axisMassTag = {AxisSpec{125, 0.f, 2.5f}, AxisSpec{100, constants::physics::MassPhi - 0.05f, constants::physics::MassPhi + 0.05f}, AxisSpec{200, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f}, AxisSpec{200, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f}, AxisSpec{200, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f}};

    std::string trackTypes[aod::tagandprobe::TrackTypes::NTrackTypes] = {"ItsTpc", "ItsTpcNoIb", "Tpc", "Its"};
    std::string tagChannels[aod::tagandprobe::TagChannels::NTagChannels] = {"DplusToKPiPi", "DsOrDplusToKKPi", "DstarPlusToDzeroPi", "DstarMinusToDzeroBarPi", "DstarChargedToDzeroToKK"};

    for (int iChannel{0}; iChannel < aod::tagandprobe::TagChannels::NTagChannels; ++iChannel) {
      for (int iTrackType{0}; iTrackType < aod::tagandprobe::TrackTypes::NTrackTypes; ++iTrackType) {
        histos[iChannel][iTrackType] = registry.add<THnSparse>(Form("h%sVsPtProbeTag_%s", tagChannels[iChannel].data(), trackTypes[iTrackType].data()),
                                                               "; #it{p}_{T}(D) (GeV/#it{c}); #it{p}_{T}(tag) (GeV/#it{c}); #it{p}_{T}(probe) (GeV/#it{c}); #it{p}_{T}^{TPC in}(probe) (GeV/#it{c}); #it{M}(D) (GeV/#it{c}^{2}); #it{M}(tag) (GeV/#it{c}^{2}); #it{#eta}(probe); #it{N}_{cross rows}^{TPC}(probe); #chi^{2}/#it{N}_{clusters}^{TPC}(probe); #it{N}_{clusters}^{ITS}(probe);",
                                                               HistType::kTHnSparseF, {axisPtD, axisPtTag, axisPtProbe, axisPtProbe, axisMass[iChannel], axisMassTag[iChannel], axisEtaProbe, axisNumCrossRowTpc, axisTpcChi2PerClus, axisNumCluIts});
      }
    }
  }

  template <typename TTrack>
  void computeInvariantMass(TTrack const& trackFirst, TTrack const& trackSecond, TTrack const& trackThird, const uint8_t channel, float& ptTag, float& invMassTag, float& ptD, float& invMass)
  {
    std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
    std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};
    std::array<float, 3> pVecTrackThird{trackThird.px(), trackThird.py(), trackThird.pz()};
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond, pVecTrackThird};
    auto arrMomentumTag = std::array{pVecTrackFirst, pVecTrackSecond};
    ptTag = RecoDecay::pt(RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond));
    ptD = RecoDecay::pt(RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond, pVecTrackThird));
    invMass = RecoDecay::m(arrMomentum, masses[channel]);
    auto massesTagDau = std::array{masses[channel][0], masses[channel][1]};
    invMassTag = RecoDecay::m(arrMomentumTag, massesTagDau);

    if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
      invMass -= invMassTag;
    }
  }

  template <typename TTrackIndices, typename TTrack, typename TTracks>
  void loopOverThirdTrack(TTrackIndices const& groupedTrackThirdIndices, TTracks const& /*tracks*/, TTrack const& trackFirst, TTrack const& trackSecond, const uint8_t channel)
  {
    for (const auto& trackIndex : groupedTrackThirdIndices) {
      auto trackThird = trackIndex.template track_as<TTracks>();
      if (trackThird.globalIndex() == trackFirst.globalIndex() || trackThird.globalIndex() == trackSecond.globalIndex()) {
        continue;
      }
      if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi && trackThird.signed1Pt() * trackFirst.signed1Pt() > 0.) { // must be opposite sign
        continue;
      }
      if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi && trackThird.signed1Pt() < 0.) { // must be positive
        continue;
      }
      if (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi && trackThird.signed1Pt() > 0.) { // must be negative
        continue;
      }
      auto ptTrackThird = trackThird.pt();
      auto ptTpcInnerTrackThird = trackThird.tpcInnerParam() / std::sqrt(1.f + trackThird.tgl() * trackThird.tgl());
      auto etaTrackThird = trackThird.eta();
      auto numTpcCrossRowTrackThird = trackThird.tpcNClsCrossedRows();
      auto numTpcChi2NumCluTrackThird = trackThird.tpcChi2NCl();
      auto numItsCluTrackThird = trackThird.itsNCls();
      float invMass{-1.f}, invMassTag{-1.f}, ptTag{-1.f}, ptD{-1.f};
      computeInvariantMass(trackFirst, trackSecond, trackThird, channel, ptTag, invMassTag, ptD, invMass);
      if ((channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) && invMass > 0.17f) {
        continue;
      } else if ((channel == aod::tagandprobe::TagChannels::DplusToKPiPi || channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) && (invMass < 1.65f || invMass > 2.10f)) {
        continue;
      }
      for (int iTrackType{0}; iTrackType < aod::tagandprobe::TrackTypes::NTrackTypes; ++iTrackType) {
        if (trackSelector[iTrackType].IsSelected(trackThird)) {
          histos[channel][iTrackType]->Fill(ptD, ptTag, ptTrackThird, ptTpcInnerTrackThird, invMass, invMassTag, etaTrackThird, numTpcCrossRowTrackThird, numTpcChi2NumCluTrackThird, numItsCluTrackThird);
        }
      }
    }
  }

  void processCombinatorialDplusToKaPiPi(aod::Collisions const& collisions,
                                         aod::PiPiFromDpTags const& tagsPiPi,
                                         aod::TrackAssoc const& trackIndices,
                                         TracksWithDca const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D+ -> pi+pi+K- and c.c.
      auto groupedTagsPiPi = tagsPiPi.sliceBy(tagsPiPiPerCollision, thisCollId);
      for (const auto& tagPiPi : groupedTagsPiPi) {
        auto trackFirst = tagPiPi.track0_as<TracksWithDca>();
        auto trackSecond = tagPiPi.track1_as<TracksWithDca>();
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DplusToKPiPi);
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDplusToKaPiPi, "Process combinatorial of tagged 2-pion vertices with additional track", true);

  void processCombinatorialDsToPhiPi(aod::Collisions const& collisions,
                                     aod::KaKaFromDspTags const& tagsKaKa,
                                     aod::TrackAssoc const& trackIndices,
                                     TracksWithDca const& tracks)
  {
    for (const auto& collision : collisions) {
      // Ds+/D+ -> phi(->K+K-)pi+ and c.c.
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto groupedTagsKaKa = tagsKaKa.sliceBy(tagsKaKaPerCollision, thisCollId);
      for (const auto& tagKaKa : groupedTagsKaKa) {
        auto trackFirst = tagKaKa.track0_as<TracksWithDca>();
        auto trackSecond = tagKaKa.track1_as<TracksWithDca>();
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DsOrDplusToKKPi);
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDsToPhiPi, "Process combinatorial of tagged 2-kaon (phi) vertices with additional track", true);

  void processCombinatorialDstarToDzeroPi(aod::Collisions const& collisions,
                                          aod::PiKaFromDzTags const& tagsPiKa,
                                          aod::KaPiFromDzTags const& tagsKaPi,
                                          aod::TrackAssoc const& trackIndices,
                                          TracksWithDca const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D*+ -> D0(->pi+K-)pi+
      auto groupedTagsPiKa = tagsPiKa.sliceBy(tagsPiKaPerCollision, thisCollId);
      for (const auto& tagPiKa : groupedTagsPiKa) {
        auto trackFirst = tagPiKa.track0_as<TracksWithDca>();  // positive --> pion
        auto trackSecond = tagPiKa.track1_as<TracksWithDca>(); // negative --> kaon
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi);
      }
      // D*- -> D0bar(->K+pi-)pi-
      auto groupedTagsKaPi = tagsKaPi.sliceBy(tagsKaPiPerCollision, thisCollId);
      for (const auto& tagKaPi : groupedTagsKaPi) {
        auto trackFirst = tagKaPi.track0_as<TracksWithDca>();  // positive --> kaon
        auto trackSecond = tagKaPi.track1_as<TracksWithDca>(); // negative --> pion
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi);
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDstarToDzeroPi, "Process combinatorial of tagged pion-kaon (D0) vertices with additional track", true);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(ProbeThirdTrack, processDummy, "Dummy process function that does nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TagTwoProngDisplacedVertices>(cfgc));
  workflow.push_back(adaptAnalysisTask<ProbeThirdTrack>(cfgc));
  return workflow;
}
