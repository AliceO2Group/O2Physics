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

/// \file HadNucleiFemto.cxx
/// \brief Analysis task for Nuclei-Hadron femto analysis
/// \author CMY
/// \date 2025-04-10

#include "PWGCF/Femto/FemtoNuclei/DataModel/HadronNucleiTables.h"
#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/GenVector/Boost.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using CollBracket = o2::math_utils::Bracket<int>;
using HyperCandidates = aod::DataHypCandsWColl;
using HyperCandidatesMC = aod::MCHypCands;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using HadHyperCollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::FT0Mults>;
using HadHyperCollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::TOFSignal, aod::TOFEvTime>;
using TrackCandidatesMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;

namespace
{
constexpr std::array<double, 12> BetheBlochDefault{-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09,
                                                   -321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09};
const std::vector<std::string> betheBlochParticleNames{"De", "He3"};
const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr std::array<float, 9> tmpRadiiTPC{{85.f, 105.f, 125.f, 145.f, 165.f, 185.f, 205.f, 225.f, 245.f}};
constexpr int DeuteronPDG = o2::constants::physics::Pdg::kDeuteron;
constexpr int TritonPDG = o2::constants::physics::Pdg::kTriton;
constexpr int He3PDG = o2::constants::physics::Pdg::kHelium3;
constexpr int HyperTritonPDG = o2::constants::physics::Pdg::kHyperTriton;
constexpr float He3TPCChi2NClMin = 0.5f;
using PairLorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>;

enum Selections {
  kNoCuts = 0,
  kTrackCuts,
  kPID,
  kAll
};

} // namespace

struct HadNucandidate {

  [[nodiscard]] float recoPtNu() const { return signNu * std::hypot(momNu[0], momNu[1]); }
  [[nodiscard]] float recoPhiNu() const { return std::atan2(momNu[1], momNu[0]); }
  [[nodiscard]] float recoEtaNu() const { return std::asinh(momNu[2] / std::abs(recoPtNu())); }
  [[nodiscard]] float recoPtHad() const { return signHad * std::hypot(momHad[0], momHad[1]); }
  [[nodiscard]] float recoPhiHad() const { return std::atan2(momHad[1], momHad[0]); }
  [[nodiscard]] float recoEtaHad() const { return std::asinh(momHad[2] / std::abs(recoPtHad())); }

  std::array<float, 3> momNu = {99.f, 99.f, 99.f};
  std::array<float, 3> momHad = {99.f, 99.f, 99.f};

  float ptHe3 = 1.f;
  float etaHe3 = 1.f;
  float signNu = 1.f;
  float signHad = 1.f;
  float invMass = -10.f;
  float dcaxyNu = -10.f;
  float dcazNu = -10.f;
  float dcaxyHad = -10.f;
  float dcazHad = -10.f;
  float dcaPair = -10.f;

  uint16_t tpcSignalNu = 0u;
  uint16_t tpcSignalHad = 0u;
  float momNuTPC = -99.f;
  float momHadTPC = -99.f;
  uint8_t nTPCClustersNu = 0u;
  uint8_t nTPCClustersHad = 0u;
  uint8_t nTPCCrossedRowsNu = 0u;
  uint8_t nTPCCrossedRowsHad = 0u;
  uint8_t sharedClustersNu = 0u;
  uint8_t sharedClustersHad = 0u;
  float chi2TPCNu = -10.f;
  float chi2TPCHad = -10.f;
  float nSigmaNu = -10.f;
  float nSigmaHad = -10.f;
  float nSigmaTOFNu = -10.f;
  float nSigmaITSNu = -10.f;
  float nSigmaTOFHad = -10.f;
  float nSigmaITSHad = -10.f;
  float nSigmaTPCHadPi = -10.f;
  float nSigmaTPCHadKa = -10.f;
  float nSigmaTPCHadPr = -10.f;
  float nSigmaTOFHadPi = -10.f;
  float nSigmaTOFHadKa = -10.f;
  float nSigmaTOFHadPr = -10.f;
  float tpcPrnsigma = -10.f;
  float tofPrnsigma = -10.f;
  uint32_t pidTrkNu = 0xFFFFF; // PID in tracking
  uint32_t pidTrkHad = 0xFFFFF;
  float massTOFNu = -10;
  float massTOFHad = -10;
  uint32_t itsClSizeNu = 0u;
  uint32_t itsClSizeHad = 0u;

  uint8_t nClsItsNu = 0u;
  uint8_t nClsItsHad = 0u;

  bool isBkgUS = false; // unlike sign
  bool isBkgEM = false; // event mixing

  int trackIDNu = -1;
  int trackIDHad = -1;

  float deltaEta = -99.f;
  float deltaPhi = -99.f;
  float kstar = 1.f;
  float mT = 1.f;

  // collision information
  int32_t collisionID = 0;
  float cent = 1.f;
};

struct BufferedTrack {
  [[nodiscard]] float pt() const { return ptValue; }
  [[nodiscard]] float eta() const { return etaValue; }
  [[nodiscard]] float phi() const { return phiValue; }
  [[nodiscard]] int8_t sign() const { return signValue; }

  std::array<float, 3> momentum{};
  float ptValue{0.f};
  float etaValue{0.f};
  float phiValue{0.f};
  int8_t signValue{0};
  float dcaXY{-10.f};
  float dcaZ{-10.f};
  uint16_t tpcSignal{0u};
  float tpcInnerParam{-99.f};
  uint8_t tpcNClsFound{0u};
  uint8_t tpcNClsCrossedRows{0u};
  uint8_t tpcNClsShared{0u};
  uint8_t itsNCls{0u};
  float tpcChi2NCl{-10.f};
  float nSigmaTPC{-10.f};
  float nSigmaTOF{-10.f};
  float nSigmaITS{-10.f};
  float nSigmaTPCHadPi{-10.f};
  float nSigmaTPCHadKa{-10.f};
  float nSigmaTPCHadPr{-10.f};
  float nSigmaTOFHadPi{-10.f};
  float nSigmaTOFHadKa{-10.f};
  float nSigmaTOFHadPr{-10.f};
  float massTOF{-10.f};
  uint32_t pidForTracking{0xFFFFFu};
  uint32_t itsClusterSizes{0u};
  int64_t trackId{-1};
};

struct BufferedCollision {
  int64_t eventId{-1};
  float posZ{0.f};
  uint16_t numContrib{0u};
  float centFT0C{0.f};
  float multFT0C{0.f};
  std::vector<BufferedTrack> nuclei;
  std::vector<BufferedTrack> hadrons;
};

struct HadNucleiFemto {

  Produces<aod::HadronNucleiTable> mOutputDataTable;
  Produces<aod::HadronNucleiTableMC> mOutputMCTable;
  Produces<aod::HadronHyperTable> mOutputHadHyperDataTable;
  Produces<aod::HadronHyperTableMC> mOutputHadHyperMCTable;
  Produces<aod::HadronNucleiMult> mOutputMultiplicityTable;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"species"};
    // Particle species configuration
    Configurable<int> settingNuPDGCode{"settingNuPDGCode", static_cast<int>(DeuteronPDG), "Nucleus - PDG code"};
    Configurable<int> settingHadPDGCode{"settingHadPDGCode", 211, "Hadron - PDG code"};
  } species;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"eventMixing"};
    // Event selection and mixing configuration
    Configurable<float> settingCutVertex{"settingCutVertex", 10.0f, "Accepted z-vertex range"};
    Configurable<int> settingNoMixedEvents{"settingNoMixedEvents", 5, "Number of mixed events per event"};
    Configurable<bool> settingEnableBkgUS{"settingEnableBkgUS", false, "Enable US background"};
    Configurable<bool> settingSaveUSandLS{"settingSaveUSandLS", true, "Save All Pairs"};
  } eventMixing;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"trackCut"};
    // Common track-quality cuts
    Configurable<float> settingCutEta{"settingCutEta", 0.8f, "Eta cut on daughter track"};
    Configurable<float> settingCutNCls{"settingCutNCls", 5.0f, "Minimum ITS Ncluster for tracks"};
    Configurable<float> settingCutChi2tpcLow{"settingCutChi2tpcLow", 0.5f, "Low cut on TPC chi2"};
    Configurable<float> settingCutChi2tpcHigh{"settingCutChi2tpcHigh", 4.f, "High cut on TPC chi2"};
    Configurable<float> settingCutChi2NClITS{"settingCutChi2NClITS", 36.f, "Maximum ITS Chi2 for tracks"};
  } trackCut;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"hadronPid"};
    // Hadron purity and PID cuts
    Configurable<float> settingCutPinMinTOFHad{"settingCutPinMinTOFHad", 0.5f, "Minimum Pin to apply the TOF cut on hadrons"};
    Configurable<float> settingCutNsigmaTPCHad{"settingCutNsigmaTPCHad", 3.0f, "Value of the TPC Nsigma cut on Had"};
    Configurable<float> settingCutNsigmaTOFHad{"settingCutNsigmaTOFHad", 3.0f, "Value of the hsdron TOF Nsigma cut"};
    Configurable<float> settingCutNsigTPCPrMin{"settingCutNsigTPCPrMin", 3.0f, "Minimum TPC Pr Nsigma cut for rejection"};
    Configurable<float> settingCutNsigTPCPiMin{"settingCutNsigTPCPiMin", 3.0f, "Minimum TPC Pi Nsigma cut for rejection"};
    Configurable<float> settingCutNsigTOFPrMin{"settingCutNsigTOFPrMin", 3.0f, "Minimum TOF Pr Nsigma cut for rejection"};
    Configurable<float> settingCutNsigTOFPiMin{"settingCutNsigTOFPiMin", 3.0f, "Minimum TOF Pi Nsigma cut for rejection"};
    Configurable<float> settingHadptMin{"settingHadptMin", 0.14f, "Minimum pT for the reference pion track cuts"};
    Configurable<float> settingHadptMax{"settingHadptMax", 2.5f, "Maximum pT for the reference pion track cuts"};
    Configurable<int> settingPionITSInnerBarrelMin{"settingPionITSInnerBarrelMin", 3, "Minimum ITS inner barrel clusters for the reference pion track cuts"};
    Configurable<int> settingPionITSNClsMin{"settingPionITSNClsMin", 7, "Minimum ITS clusters for the reference pion track cuts"};
    Configurable<int> settingPionTPCNClsFoundMin{"settingPionTPCNClsFoundMin", 80, "Minimum found TPC clusters for the reference pion track cuts"};
    Configurable<int> settingPionTPCCrossedRowsMin{"settingPionTPCCrossedRowsMin", 90, "Minimum crossed TPC rows for the reference pion track cuts"};
    Configurable<float> settingPionDCAxyOffset{"settingPionDCAxyOffset", 0.004f, "DCAxy offset for the reference pion track cuts"};
    Configurable<float> settingPionDCAxyPtCoeff{"settingPionDCAxyPtCoeff", 0.013f, "DCAxy 1/pT coefficient for the reference pion track cuts"};
    Configurable<float> settingPionDCAzOffset{"settingPionDCAzOffset", 0.004f, "DCAz offset for the reference pion track cuts"};
    Configurable<float> settingPionDCAzPtCoeff{"settingPionDCAzPtCoeff", 0.013f, "DCAz 1/pT coefficient for the reference pion track cuts"};
    Configurable<float> settingPionMomCombMin{"settingPionMomCombMin", 0.5f, "Minimum momentum to use combined TPC+TOF PID for reference pions"};
    Configurable<float> settingPionTPCNsigMax{"settingPionTPCNsigMax", 3.0f, "Maximum TPC n-sigma for reference pions below the TOF threshold"};
    Configurable<float> settingPionCombNsigMax{"settingPionCombNsigMax", 3.0f, "Maximum combined TPC+TOF n-sigma for reference pions"};
  } hadronPid;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"deuteronPid"};
    // Deuteron purity and PID cuts
    Configurable<float> settingCutPinMinDe{"settingCutPinMinDe", 0.0f, "Minimum Pin for De"};
    Configurable<float> settingCutDeptMin{"settingCutDeptMin", 0.6f, "Minimum PT cut on De"};
    Configurable<float> settingCutDeptMax{"settingCutDeptMax", 1.6f, "Maximum PT cut on De"};
    Configurable<float> settingCutPinMinTOFITSDe{"settingCutPinMinTOFITSDe", 1.2f, "Minimum p to apply the TOF ITS cut on De"};
    Configurable<float> settingCutNsigmaTPCDe{"settingCutNsigmaTPCDe", 2.5f, "Value of the TPC Nsigma cut on De"};
    Configurable<float> settingCutNsigmaITSDe{"settingCutNsigmaITSDe", 2.5f, "Value of the ITD Nsigma cut on De"};
    Configurable<float> settingCutNsigmaTOFTPCDe{"settingCutNsigmaTOFTPCDe", 2.5f, "Value of the De TOF TPC combNsigma cut"};
    Configurable<bool> settingReqSingleNsig{"settingReqSingleNsig", false, "If true, also require individual TPC and TOF n-sigma cuts in branches using combined TPC+TOF PID"};
    Configurable<bool> settingUseProtonMassForKstarMt{"settingUseProtonMassForKstarMt", false, "If true, use proton mass instead of deuteron mass for kstar and mT"};
  } deuteronPid;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"helium3Pid"};
    Configurable<float> settingRigidityMinHe3{"settingRigidityMinHe3", 0.8f, "Minimum He3 TPC rigidity"};
    Configurable<int> settingTPCNClsFoundMinHe3{"settingTPCNClsFoundMinHe3", 110, "Minimum found TPC clusters for He3"};
    Configurable<int> settingTPCCrossedRowsMinHe3{"settingTPCCrossedRowsMinHe3", 70, "Minimum crossed TPC rows for He3"};
    Configurable<float> settingTPCNSigmaMaxHe3{"settingTPCNSigmaMaxHe3", 3.f, "Maximum absolute TPC n-sigma for He3"};
    Configurable<float> settingITSNSigmaMinHe3{"settingITSNSigmaMinHe3", -1.5f, "Minimum ITS n-sigma for He3"};
  } helium3Pid;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"tritonPid"};
    // Triton track-quality, purity and PID cuts
    Configurable<float> settingPIDMomentumSplitTr{"settingPIDMomentumSplitTr", 2.f, "Momentum separating low- and high-momentum triton TPC PID"};
    Configurable<float> settingCutTPCNsigmaLowPTr{"settingCutTPCNsigmaLowPTr", 3.f, "Maximum absolute TPC n-sigma for tritons below the momentum split"};
    Configurable<float> settingCutITSNsigmaLowPTr{"settingCutITSNsigmaLowPTr", 3.f, "Maximum absolute ITS n-sigma for tritons below the momentum split"};
    Configurable<float> settingCutTPCNsigmaHighPMinTr{"settingCutTPCNsigmaHighPMinTr", -2.f, "Minimum TPC n-sigma for tritons above the momentum split"};
    Configurable<float> settingCutTPCNsigmaHighPMaxTr{"settingCutTPCNsigmaHighPMaxTr", 3.f, "Maximum TPC n-sigma for tritons above the momentum split"};
    Configurable<float> settingTOFMassMomentumMinTr{"settingTOFMassMomentumMinTr", 1.2f, "Minimum momentum to apply the triton TOF mass cut"};
    Configurable<float> settingTOFMassMinTr{"settingTOFMassMinTr", 2.5f, "Minimum triton TOF mass"};
    Configurable<float> settingTOFMassMaxTr{"settingTOFMassMaxTr", 3.4f, "Maximum triton TOF mass"};
    Configurable<float> settingTPCRejectNsig{"settingTPCRejectNsig", 3.f, "Minimum absolute TPC n-sigma from deuteron, proton and pion hypotheses"};
  } tritonPid;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"CPR"};
    // Close pair rejection controls
    Configurable<bool> settingEnableClosePairRejection{"settingEnableClosePairRejection", false, "Enable close pair rejection for nucleus-hadron track pairs"};
    Configurable<float> settingClosePairDeltaPhiMax{"settingClosePairDeltaPhiMax", 0.01f, "Maximum delta phi star for close pair rejection"};
    Configurable<float> settingClosePairDeltaEtaMax{"settingClosePairDeltaEtaMax", 0.01f, "Maximum delta eta for close pair rejection"};
    Configurable<int> settingClosePairRadiusMode{"settingClosePairRadiusMode", 1, "Close pair rejection mode: 0 = PV, 1 = average phi star, 2 = specific TPC radius"};
    Configurable<float> settingClosePairSpecificRadius{"settingClosePairSpecificRadius", 85.f, "TPC radius in cm used when close pair rejection mode is 2"};
  } CPR;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"mc"};
    Configurable<bool> settingRequireSel8{"settingRequireSel8", true, "Apply the same sel8 event selection to reconstructed MC as to data"};
    Configurable<bool> settingRequireTruthSpecies{"settingRequireTruthSpecies", true, "Store only truth-matched pion-nucleus pairs"};
    Configurable<bool> settingRequireSameMCCollision{"settingRequireSameMCCollision", true, "Require the two truth particles to come from the same MC collision"};
    Configurable<bool> settingRequireRecoMCCollisionMatch{"settingRequireRecoMCCollisionMatch", true, "Require both truth particles to match the reconstructed collision MC label"};
    Configurable<bool> settingRequirePhysicalPrimaries{"settingRequirePhysicalPrimaries", false, "Store only pairs in which both truth particles are physical primaries"};
  } mc;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"hypertriton"};
    // Hypertriton-specific cuts
    Configurable<float> settingHypMassMin{"settingHypMassMin", 2.94f, "Minimum hypertriton invariant mass"};
    Configurable<float> settingHypMassMax{"settingHypMassMax", 3.10f, "Maximum hypertriton invariant mass"};
  } hypertriton;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"hadHyper"};
    Configurable<bool> enableMixing{"enableMixing", true, "Build mixed-event pion-hypertriton pairs"};
    Configurable<float> maxOutputKstar{"maxOutputKstar", -1.f, "Maximum pair k* (GeV/c); negative saves all selected pairs"};
  } hadHyper;

  HistogramRegistry hadHyperRegistry{
    "hadHyperRegistry",
    {{"hSE", "Raw same-event pairs;k*;Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"hME", "Raw mixed-event pairs;k*;Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"hMass", "Unique selected candidates;mass He3-pi;Entries", {HistType::kTH1D, {{160, 2.94, 3.10}}}},
     {"hSelfPairs", "Rejected daughter reuse;reason (0=He3,1=pion,2=truth);k*", {HistType::kTH2F, {{3, -0.5, 2.5}, {300, 0., 3.}}}},
     {"hSameEventSelfPairs", "Same-event rejected daughter reuse;reason (0=He3,1=pion,2=truth);k*", {HistType::kTH2F, {{3, -0.5, 2.5}, {300, 0., 3.}}}},
     {"hCandidatePairMultiplicitySE", "Accepted same-event hadron pairs per hypertriton candidate;pairs;centrality", {HistType::kTH2F, {{200, -0.5, 199.5}, {10, 0., 100.}}}},
     {"hCandidatePairMultiplicityME", "Accepted mixed-event hadron pairs per hypertriton candidate;pairs;centrality", {HistType::kTH2F, {{200, -0.5, 199.5}, {10, 0., 100.}}}},
     {"hMixEventDeltaPosZVsCent", "Mixed-event #Delta z_{vtx} vs hypertriton centrality;hypertriton CentFT0C;#Delta z_{vtx}", {HistType::kTH2F, {{100, 0., 100.}, {120, -30., 30.}}}},
     {"hMixEventDeltaCentFT0CVsCent", "Mixed-event #Delta CentFT0C vs hypertriton centrality;hypertriton CentFT0C;#Delta CentFT0C", {HistType::kTH2F, {{100, 0., 100.}, {200, -100., 100.}}}},
     {"hMixingDepth", "Available partner events;events;anchor events", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
     {"hPoolFlow", "Mixing pool selection;0=selected,1=in range;events", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
     {"hClosePairBeforePV", "Hadron-hypertriton PV close-pair QA before rejection;#Delta#eta;#Delta#varphi", {HistType::kTH2F, {{200, -0.1, 0.1}, {200, -0.1, 0.1}}}},
     {"hClosePairAfterPV", "Hadron-hypertriton PV close-pair QA after rejection;#Delta#eta;#Delta#varphi", {HistType::kTH2F, {{200, -0.1, 0.1}, {200, -0.1, 0.1}}}},
     {"hDaughterHeTPC", "He3 daughter QA per candidate;TPC rigidity;TPC signal", {HistType::kTH2F, {{100, 0., 10.}, {300, 0., 1500.}}}},
     {"hDaughterPiTPC", "Pion daughter QA per candidate;TPC rigidity;TPC signal", {HistType::kTH2F, {{100, 0., 10.}, {300, 0., 1500.}}}},
     {"hHadronTPC", "Pair hadron QA per selected track;TPC rigidity;TPC signal", {HistType::kTH2F, {{100, 0., 10.}, {300, 0., 1500.}}}},
     {"MC/hKstarRecVsGenHyperReco", "Hadron-hypertriton k* response for reconstructed true hypertritons;generated k* (GeV/#it{c});reconstructed k* (GeV/#it{c})", {HistType::kTH2F, {{300, 0., 3.}, {300, 0., 3.}}}},
     {"MC/hKstarResolutionHyperReco", "Hadron-hypertriton k* resolution for reconstructed true hypertritons;reconstructed k* (GeV/#it{c});k*_{reco}-k*_{gen} (GeV/#it{c})", {HistType::kTH2F, {{300, 0., 3.}, {200, -0.2, 0.2}}}},
     {"MC/hPrimaryHadronVsKstarDen", "Truth-matched selected hadron denominator;k* (GeV/#it{c});Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"MC/hPrimaryHadronVsKstarNum", "Physical-primary selected hadron numerator;k* (GeV/#it{c});Entries", {HistType::kTH1D, {{300, 0., 3.}}}},
     {"MC/hPrimaryHadronVsCentDen", "Truth-matched selected hadron denominator;hypertriton CentFT0C;Entries", {HistType::kTH1D, {{100, 0., 100.}}}},
     {"MC/hPrimaryHadronVsCentNum", "Physical-primary selected hadron numerator;hypertriton CentFT0C;Entries", {HistType::kTH1D, {{100, 0., 100.}}}}}};

  // Tuple order follows the column blocks in HadronNucleiTables.h.
  using HadHyperTrackInfo = std::tuple<float, float, float, int8_t, float, uint8_t, uint8_t, float, uint32_t, float, bool, float, float, float, float, float, float>;
  using HadHyperEventInfo = std::tuple<int32_t, float, float, float, float, int, float, uint16_t, float, float, float, float>;
  using HadHyperCandidateInfo = std::tuple<bool, float, float, float, float, float, float, float, float, float, float, float, uint8_t, uint8_t, float, float, uint16_t, uint16_t, float, float, float, uint32_t, uint32_t, float, float, float>;
  using HadHyperDataInfo = decltype(std::tuple_cat(std::declval<std::tuple<bool, int>>(),
                                                   std::declval<HadHyperEventInfo>(),
                                                   std::declval<HadHyperCandidateInfo>(),
                                                   std::declval<HadHyperTrackInfo>()));
  using HadHyperMCInfo = std::tuple<float, float, float, float, float, float, bool, bool, bool, bool, bool, int16_t,
                                    float, bool, bool, float, float, float, float, float, float, bool, int16_t,
                                    bool, bool, bool, bool, bool>;

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
    std::array<float, 3> momentum{0.f, 0.f, 0.f};
  };

  struct HadHyperCandidate {
    HadHyperCandidateInfo info{};
    HadHyperParticleTruth hyperTruth{};
    HadHyperParticleTruth heTruth{};
    HadHyperParticleTruth decayPionTruth{};
    int64_t sourceCandidateId{-1};
    int64_t heTrackId{-1};
    int64_t pionTrackId{-1};
    int16_t statusCode{0};
    bool isReco{true};
    bool isSignal{false};
    bool isRecoMCCollision{false};
    bool isSurvEvSel{false};
    bool isTwoBodyDecay{false};
    uint8_t isFakeHeOnITSLayer{0u};
    float genPt{-1.f};
    float genEta{-1.f};
    float genPhi{-1.f};
    float genPtHe3{-1.f};
    std::array<float, 3> genDecVtx{-1.f, -1.f, -1.f};
    std::array<float, 3> momentum{};
    std::array<float, 3> genMomentum{};
    float mass{0.f};
    float etaHe{0.f};
    float phiHe{0.f};
    float etaPi{0.f};
    float phiPi{0.f};
    bool isMatter{false};

    float pt() const { return std::hypot(momentum[0], momentum[1]); }
    float eta() const
    {
      const float transverseMomentum = pt();
      return transverseMomentum > 0.f ? std::asinh(momentum[2] / transverseMomentum) : 999.f;
    }
    float phi() const { return std::atan2(momentum[1], momentum[0]); }
    int8_t sign() const { return isMatter ? 1 : -1; }
  };

  struct HadHyperHadron {
    HadHyperTrackInfo info{};
    HadHyperParticleTruth truth{};
    std::vector<int64_t> motherIds;
    int64_t sourceId{-1};
    std::array<float, 3> momentum{};
    float tpcInnerParam{0.f};
    float tpcSignal{0.f};
    float etaValue{0.f};
    float phiValue{0.f};
    int8_t signValue{0};

    float pt() const { return std::hypot(momentum[0], momentum[1]); }
    float eta() const { return etaValue; }
    float phi() const { return phiValue; }
    int8_t sign() const { return signValue; }
  };

  struct HadHyperEvent {
    HadHyperEventInfo info{};
    int64_t mcCollisionId{-1};
    bool hasMCCollision{false};
    float centrality{0.f};
    std::vector<HadHyperCandidate> candidates;
    std::vector<HadHyperHadron> hadrons;
  };

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"output"};
    // Output and QA controls
    Configurable<bool> settingFillTable{"settingFillTable", false, "Enable table filling"};
    Configurable<bool> settingFillTableLowKstarOnly{"settingFillTableLowKstarOnly", false, "If true, apply the configured low-kstar threshold to pion-deuteron and pion-triton pairs"};
    Configurable<float> settingFillTablePiDeKstarMax{"settingFillTablePiDeKstarMax", 0.1f, "Maximum kstar for pion-deuteron pairs written to the output table"};
    Configurable<float> settingFillTablePiTrKstarMax{"settingFillTablePiTrKstarMax", 0.1f, "Maximum kstar for pion-triton pairs written to the output table"};
    Configurable<bool> settingFillMultiplicity{"settingFillMultiplicity", false, "Fill multiplicity table"};
    Configurable<bool> settingUseBBcomputeDeNsigma{"settingUseBBcomputeDeNsigma", false, "Use BB params to compute De TPC Nsigma"};
  } output;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"zorro"};
    // Zorro
    Configurable<bool> settingSkimmedProcessing{"settingSkimmedProcessing", false, "Skimmed dataset processing"};
  } zorro;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"ccdb"};
    // CCDB options
    Configurable<double> settingDbz{"settingDbz", -999, "bz field, -999 is automatic"};
    Configurable<std::string> settingCcdburl{"settingCcdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> settingGrpPath{"settingGrpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> settingGrpmagPath{"settingGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> settingLutPath{"settingLutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> settingGeoPath{"settingGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> settingPidPath{"settingPidPath", "", "Path to the PID response object"};
  } ccdb;

  struct : o2::framework::ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"pidCalibration"};
    // PID calibration and material correction controls
    Configurable<LabeledArray<double>> settingBetheBlochParams{"settingBetheBlochParams", {BetheBlochDefault.data(), 2, 6, betheBlochParticleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for the selected nucleus"};
    Configurable<bool> settingCompensatePIDinTracking{"settingCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
    Configurable<int> settingMaterialCorrection{"settingMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Material correction type"};
  } pidCalibration;

  Preslice<TrackCandidates> mPerCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> mPerColMC = aod::track::collisionId;
  PresliceUnsorted<HyperCandidates> hypPerCol = o2::aod::hyperrec::collisionId;
  PresliceUnsorted<HyperCandidatesMC> hypPerColMC = o2::aod::hyperrec::collisionId;

  // binning for EM background
  ConfigurableAxis axisVertex{"axisVertex", {30, -10, 10}, "Binning for vtxz"};
  ConfigurableAxis axisCentrality{"axisCentrality", {40, 0, 100}, "Binning for centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  std::array<float, 6> mBBparamsNucleus{};
  float mMassHad{0.f};
  std::vector<bool> mGoodCollisions;
  std::vector<SVCand> mTrackPairs;
  std::unordered_map<int, std::deque<BufferedCollision>> mMixingPools;
  std::unordered_map<int, std::deque<HadHyperEvent>> mHyperMixingPools;
  std::unordered_map<int, std::deque<HadHyperEvent>> mHyperMCMixingPools;
  int64_t mNextMixedEventId{0};
  int mMixingRunNumber{-1};
  int mHyperMixingRunNumber{-1};
  int mHyperMCMixingRunNumber{-1};
  o2::vertexing::DCAFitterN<2> mFitter;

  int mRunNumber{0};
  float mDbz{0.f};
  Service<o2::ccdb::BasicCCDBManager> mCcdb{};
  Zorro mZorro;
  OutputObj<ZorroSummary> mZorroSummary{"zorroSummary"};

  HistogramRegistry mQaRegistry{
    "QA",
    {// Event-level
     {"hVtxZ", "Vertex distribution in Z;Z (cm)", {HistType::kTH1F, {{400, -20.0, 20.0}}}},
     {"hNcontributor", "Number of primary vertex contributor", {HistType::kTH1F, {{2000, 0.0f, 2000.0f}}}},
     {"hCentrality", "Centrality", {HistType::kTH1F, {{100, 0.0f, 100.0f}}}},
     {"hSkipReasons", "Why storedEvent skipped;Reason;Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
     {"hEvents", "; Events;", {HistType::kTH1F, {{3, -0.5, 2.5}}}},
     {"hMixedEventSelections", "Mixed-event collision selection;Step;Counts", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hMixingPoolOccupancy", "Events already available in the mixing pool;Events;Counts", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
     {"hMixedNucleiPerEvent", "Selected nuclei in mixed-event input;Candidates;Events", {HistType::kTH1F, {{201, -0.5, 200.5}}}},
     {"hMixedHadronsPerEvent", "Selected hadrons in mixed-event input;Candidates;Events", {HistType::kTH1F, {{501, -0.5, 500.5}}}},

     // Candidate topology and kinematics
     {"hTrackSel", "Accepted hadron tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"hTrackSelNu", "Accepted nucleus tracks", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"hNuPairFlow", "Nucleus pair-building flow;step;counts", {HistType::kTH1F, {{3, -0.5, 2.5}}}},

     {"hdcaxyNu", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazNu", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNClsNuITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hNuPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hSingleNuPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hNuPin", "#it{p} distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hSingleNuPin", "#it{p} distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{240, -6.0f, 6.0f}}}},
     {"hNuEta", "eta distribution; #eta(Nu)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNuPhi", "phi distribution; phi(Nu)", {HistType::kTH1F, {{600, -4.0f, 4.0f}}}},

     {"hdcaxyHad", ";DCA_{xy} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hdcazHad", ";DCA_{z} (cm)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hNClsHadITS", ";N_{ITS} Cluster", {HistType::kTH1F, {{20, -10.0f, 10.0f}}}},
     {"hHadPt", "Pt distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hSingleHadPt", "#it{p}_{T} distribution; #it{p}_{T} (GeV/#it{c})", {HistType::kTH1F, {{280, -7.0f, 7.0f}}}},
     {"hHadPin", "P distribution; #it{p} (GeV/#it{c})", {HistType::kTH1F, {{120, -4.0f, 4.0f}}}},
     {"hHadEta", "eta distribution; #eta(had)", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
     {"hHadPhi", "phi distribution; phi(had)", {HistType::kTH1F, {{600, -4.0f, 4.0f}}}},
     {"h2CPRBefore", "Close pair rejection before cut; #Delta#eta; #Delta#phi^{*}", {HistType::kTH2F, {{300, -0.15f, 0.15f}, {400, -0.2f, 0.2f}}}},
     {"h2CPRAfter", "Close pair rejection after cut; #Delta#eta; #Delta#phi^{*}", {HistType::kTH2F, {{300, -0.15f, 0.15f}, {400, -0.2f, 0.2f}}}},

     // Reconstructed MC pair QA
     {"MC/hPairFlow", "MC pair flow;step;counts", {HistType::kTH1F, {{7, -0.5f, 6.5f}}}},
     {"MC/hKstarRecVsGen", "Reconstructed versus generated k*;generated k* (GeV/c);reconstructed k* (GeV/c)", {HistType::kTH2F, {{300, 0.f, 3.f}, {300, 0.f, 3.f}}}},
     {"MC/hPtNuRecVsGen", "Reconstructed versus generated signed nucleus pT;generated pT (GeV/c);reconstructed pT (GeV/c)", {HistType::kTH2F, {{280, -7.f, 7.f}, {280, -7.f, 7.f}}}},
     {"MC/hPtHadRecVsGen", "Reconstructed versus generated signed pion pT;generated pT (GeV/c);reconstructed pT (GeV/c)", {HistType::kTH2F, {{280, -7.f, 7.f}, {280, -7.f, 7.f}}}},

     // dE/dx
     {"h2dEdxNucandidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2dEdxHadcandidates", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},
     {"h2dEdx", "dEdx distribution; #it{p} (GeV/#it{c}); dE/dx (a.u.)", {HistType::kTH2F, {{200, -5.0f, 5.0f}, {100, 0.0f, 2000.0f}}}},

     // Nucleus PID
     {"h2NsigmaNuTPC", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaNuComb", "NsigmaNu TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaNuTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaNuTPC_preselecComp", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NSigmaNuITS_preselection", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NSigmaNuITS", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {120, -3.0f, 3.0f}}}},
     {"h2NsigmaNuTOF", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaNuTOF_preselection", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2MassTOFTr", "Triton TOF mass; signed #it{p} (GeV/#it{c}); m_{TOF} (GeV/#it{c}^{2})", {HistType::kTH2F, {{240, -6.0f, 6.0f}, {200, 0.0f, 5.0f}}}},

     // Hadron PID
     {"h2NsigmaHadComb", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaHadTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaHadTOF", "NsigmaHad TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadTOF_preselection", "NsigmaHad TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"h2NsigmaHadComb_preselection", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"h2NsigmaHadPrTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPiTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(pi)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadKaTPC", "NsigmaHad TPC distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(K)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPrTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(p)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadPiTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(pi)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"h2NsigmaHadKaTOF", "NsigmaHad TOF distribution; #it{p}_{T}(GeV/#it{c}); n#sigma_{TPC}(K)", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},

     // Purity
     {"purity/h2NsigmaNuTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaNuTPC_preselecComp", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NSigmaNuITS_preselection", "NsigmaNu ITS distribution; signed #it{p}_{T} (GeV/#it{c}); n#sigma_{ITS} Nu", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {120, -3.0f, 3.0f}}}},
     {"purity/h2NsigmaNuTOF_preselection", "NsigmaNu TOF distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaNuComb_preselection", "NsigmaNu TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},
     {"purity/h2NsigmaHadTPC_preselection", "NsigmaNu TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(Nu)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaHadTOF_preselection", "NsigmaHad TOF distribution; #iit{p}_{T} (GeV/#it{c}); n#sigma_{TOF}(p)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {400, -10.0f, 10.0f}}}},
     {"purity/h2NsigmaHadComb_preselection", "NsigmaHad TPCTOF comb distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{comb}(had)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {100, 0.0f, 5.0f}}}},

     // Hypertriton
     {"hHe3TPCnsigma", "NsigmaHe3 TPC distribution; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}(He3)", {HistType::kTH2F, {{280, -7.0f, 7.0f}, {200, -5.0f, 5.0f}}}},
     {"hHe3P", "Pin distribution; p (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},
     {"hHe3P_preselected", "Pin distribution_preselected; p (GeV/#it{c})", {HistType::kTH1F, {{120, -3.0f, 3.0f}}}},

     // Correlation observables
     {"hkStar_LS_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_LS_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_M", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStar_US_A", ";kStar (GeV/c)", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
     {"hkStaVsmT_LS_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_LS_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_US_M", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsmT_US_A", ";kStar (GeV/c);mT (GeV/#it{c}^{2})", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {2000, 0.8, 2.0}}}},
     {"hkStaVsCent_LS_M", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsCent_LS_A", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsCent_US_M", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hkStaVsCent_US_A", ";kStar (GeV/c);Centrality", {HistType::kTH2F, {{300, 0.0f, 3.0f}, {100, 0.0f, 100.0f}}}},
     {"hNuHadtInvMass", "; M(Nu + had) (GeV/#it{c}^{2})", {HistType::kTH1F, {{500, 2.5f, 4.5f}}}},

     // Mixed-event
     {"hisBkgEM", "; isBkgEM;", {HistType::kTH1F, {{3, -1, 2}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  void init(o2::framework::InitContext&)
  {
    const bool processHyperPairs = doprocessHyper || doprocessMCHyper;
    if (processHyperPairs && hadHyper.maxOutputKstar.value == 0.f) {
      LOG(fatal) << "Hadron-hypertriton mode requires a nonzero output k* range";
    }
    if (processHyperPairs && hadHyper.enableMixing.value && eventMixing.settingNoMixedEvents.value <= 0) {
      LOG(fatal) << "Hadron-hypertriton mixed-event mode requires enabled mixing and positive mixing depth";
    }
    if (processHyperPairs && hypertriton.settingHypMassMin.value >= hypertriton.settingHypMassMax.value) {
      LOG(fatal) << "Hadron-hypertriton mode requires settingHypMassMin < settingHypMassMax";
    }
    constexpr int closePairRadiusModePv = 0;
    constexpr int closePairRadiusModeSpecificTpc = 2;
    if (CPR.settingEnableClosePairRejection.value) {
      if (CPR.settingClosePairDeltaEtaMax.value <= 0.f || CPR.settingClosePairDeltaPhiMax.value <= 0.f) {
        LOG(fatal) << "Close-pair rejection requires positive delta-eta and delta-phi-star limits";
      }
      if (CPR.settingClosePairRadiusMode.value < closePairRadiusModePv || CPR.settingClosePairRadiusMode.value > closePairRadiusModeSpecificTpc) {
        LOG(fatal) << "Invalid close-pair radius mode " << CPR.settingClosePairRadiusMode.value << "; expected 0, 1, or 2";
      }
      if (CPR.settingClosePairRadiusMode.value == closePairRadiusModeSpecificTpc && CPR.settingClosePairSpecificRadius.value <= 0.f) {
        LOG(fatal) << "Close-pair rejection at a specific TPC radius requires a positive radius";
      }
    }

    mZorroSummary.setObject(mZorro.getZorroSummary());
    mRunNumber = 0;

    mCcdb->setURL(ccdb.settingCcdburl);
    mCcdb->setCaching(true);
    mCcdb->setLocalObjectValidityChecking();
    mCcdb->setFatalWhenNull(false);

    mFitter.setPropagateToPCA(true);
    mFitter.setMaxR(200.);
    mFitter.setMinParamChange(1e-3);
    mFitter.setMinRelChi2Change(0.9);
    mFitter.setMaxDZIni(1e9);
    mFitter.setMaxChi2(1e9);
    mFitter.setUseAbsDCA(true);
    int mat{static_cast<int>(pidCalibration.settingMaterialCorrection)};
    mFitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    if (!useTritonNucleus()) {
      const int numParticles = 5;
      const char* betheBlochLabel = nucleusBetheBlochLabel();
      for (int i = 0; i < numParticles; i++) {
        mBBparamsNucleus[i] = pidCalibration.settingBetheBlochParams->get(betheBlochLabel, Form("p%i", i));
      }
      mBBparamsNucleus[5] = pidCalibration.settingBetheBlochParams->get(betheBlochLabel, "resolution");
    }

    std::vector<std::string> selectionLabels = {"All", "Track selection", "PID"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hTrackSel"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i].c_str());
      mQaRegistry.get<TH1>(HIST("hTrackSelNu"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i].c_str());
    }

    std::vector<std::string> eventsLabels = {"All", "Selected", "Zorro selected events"};
    for (int i = 0; i < Selections::kAll; i++) {
      mQaRegistry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(i + 1, eventsLabels[i].c_str());
    }

    const std::array<std::string, 4> mixedEventLabels = {"All collisions", "Event selection", "Mixing pool", "Mixed combinations"};
    for (size_t i = 0; i < mixedEventLabels.size(); i++) {
      mQaRegistry.get<TH1>(HIST("hMixedEventSelections"))->GetXaxis()->SetBinLabel(i + 1, mixedEventLabels[i].c_str());
    }
  }

  template <bool isMC>
  void initCCDB(const aod::BCsWithTimestamps::iterator& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if constexpr (!isMC) {
      if (zorro.settingSkimmedProcessing) {
        mZorro.initCCDB(mCcdb.service, bc.runNumber(), bc.timestamp(), zorroTriggerMask());
        mZorro.populateHistRegistry(mQaRegistry, bc.runNumber());
      }
    }
    mRunNumber = bc.runNumber();
    const float defaultBzValue = -999.0f;

    // A fixed field is sufficient for CPR and DCAFitter when material
    // corrections are disabled. This also makes local MC tests independent
    // of an AliEn token when the CCDB payload is stored on Grid.
    if (ccdb.settingDbz > defaultBzValue && pidCalibration.settingMaterialCorrection == 0) {
      mDbz = ccdb.settingDbz;
      mFitter.setBz(mDbz);
      LOG(info) << "Using configured magnetic field of " << mDbz << " kZG";
      return;
    }

    auto run3GrpTimestamp = bc.timestamp();
    auto* grpo = mCcdb->getForTimeStamp<o2::parameters::GRPObject>(ccdb.settingGrpPath, run3GrpTimestamp);
    o2::parameters::GRPMagField* grpmag = nullptr;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (ccdb.settingDbz <= defaultBzValue) {
        // Fetch magnetic field from ccdb for current collision
        mDbz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3GrpTimestamp << " with magnetic field of " << mDbz << " kZG";
      } else {
        mDbz = ccdb.settingDbz;
      }
    } else {
      grpmag = mCcdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdb.settingGrpmagPath, run3GrpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdb.settingGrpmagPath << " of object GRPMagField and " << ccdb.settingGrpPath << " of object GRPObject for timestamp " << run3GrpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (ccdb.settingDbz <= defaultBzValue) {
        // Fetch magnetic field from ccdb for current collision
        mDbz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3GrpTimestamp << " with magnetic field of " << mDbz << " kZG";
      } else {
        mDbz = ccdb.settingDbz;
      }
    }
    mFitter.setBz(mDbz);
  }

  // ==================================================================================================================

  template <bool isMC, typename Tcollision>
  bool passesEventSelection(const Tcollision& collision)
  {
    // CPR uses phi* and therefore needs the magnetic field for MC as well.
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB<isMC>(bc);

    if constexpr (isMC) {
      if ((mc.settingRequireSel8.value && !collision.sel8()) || std::abs(collision.posZ()) > eventMixing.settingCutVertex) {
        return false;
      }
    } else {
      if (!collision.sel8() || std::abs(collision.posZ()) > eventMixing.settingCutVertex) {
        return false;
      }
    }

    return true;
  }

  template <typename Tcollision>
  bool passesZorroSelection(const Tcollision& collision)
  {
    if (!zorro.settingSkimmedProcessing) {
      return true;
    }
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    return mZorro.isSelected(bc.globalBC());
  }

  template <bool isMC, typename Tcollision>
  bool selectCollision(const Tcollision& collision, const aod::BCsWithTimestamps&)
  {
    mQaRegistry.fill(HIST("hEvents"), 0);

    if (!passesEventSelection<isMC>(collision)) {
      return false;
    }

    mQaRegistry.fill(HIST("hEvents"), 1);

    if constexpr (!isMC) {
      if (zorro.settingSkimmedProcessing) {
        if (!passesZorroSelection(collision)) {
          return false;
        }
        mQaRegistry.fill(HIST("hEvents"), 2);
      }
    }

    mQaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
    mQaRegistry.fill(HIST("hVtxZ"), collision.posZ());
    return true;
  }

  template <typename Ttrack>
  bool selectTrack(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > trackCut.settingCutEta) {
      return false;
    }
    const int minTPCNClsFound = 90;
    const int minTPCNClsCrossedRows = 70;
    const float crossedRowsToFindableRatio = 0.83f;
    return !(candidate.itsNCls() < trackCut.settingCutNCls ||
             candidate.tpcNClsFound() < minTPCNClsFound ||
             candidate.tpcNClsCrossedRows() < minTPCNClsCrossedRows ||
             candidate.tpcNClsCrossedRows() < crossedRowsToFindableRatio * candidate.tpcNClsFindable() ||
             candidate.tpcChi2NCl() > trackCut.settingCutChi2tpcHigh ||
             candidate.tpcChi2NCl() < trackCut.settingCutChi2tpcLow ||
             candidate.itsChi2NCl() > trackCut.settingCutChi2NClITS);
  }

  template <typename Ttrack>
  bool selectTrackPion(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > trackCut.settingCutEta) {
      return false;
    }

    const float absPt = std::abs(candidate.pt());
    if (absPt < hadronPid.settingHadptMin || absPt > hadronPid.settingHadptMax) {
      return false;
    }

    if (candidate.itsNClsInnerBarrel() < hadronPid.settingPionITSInnerBarrelMin ||
        candidate.itsNCls() < hadronPid.settingPionITSNClsMin ||
        candidate.tpcNClsFound() < hadronPid.settingPionTPCNClsFoundMin ||
        candidate.tpcNClsCrossedRows() < hadronPid.settingPionTPCCrossedRowsMin) {
      return false;
    }

    if (absPt <= 0.f) {
      return false;
    }

    const float pionDCAxyMax = hadronPid.settingPionDCAxyOffset + hadronPid.settingPionDCAxyPtCoeff / absPt;
    const float pionDCAzMax = hadronPid.settingPionDCAzOffset + hadronPid.settingPionDCAzPtCoeff / absPt;
    return !(std::abs(candidate.dcaXY()) > pionDCAxyMax || std::abs(candidate.dcaZ()) > pionDCAzMax);
  }

  template <typename Ttrack>
  bool selectTrackHadron(const Ttrack& candidate)
  {
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return selectTrackProton(candidate);
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return selectTrackPion(candidate);
    }
    return selectTrack(candidate);
  }

  template <typename Ttrack>
  bool selectTrackProton(const Ttrack& candidate)
  {
    constexpr float protonEtaMax = 0.8f;
    constexpr int protonTPCNClsFoundMin = 90;
    constexpr int protonTPCCrossedRowsMin = 80;
    constexpr float protonDCAzMax = 0.2f;

    if (std::abs(candidate.eta()) >= protonEtaMax) {
      return false;
    }

    const float absPt = std::abs(candidate.pt());
    if (absPt <= 0.f) {
      return false;
    }

    if (candidate.tpcNClsFound() <= protonTPCNClsFoundMin ||
        candidate.tpcNClsCrossedRows() <= protonTPCCrossedRowsMin) {
      return false;
    }

    const float prDCAxyMax = 105.e-3f + 30.5e-3f / std::pow(absPt, 1.1f);
    return !(std::abs(candidate.dcaXY()) >= prDCAxyMax || std::abs(candidate.dcaZ()) >= protonDCAzMax);
  }

  template <typename Ttrack>
  bool selectTrackDe(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > trackCut.settingCutEta) {
      return false;
    }

    constexpr int minTPCNClsFound = 110;
    constexpr int minTPCNClsCrossedRows = 100;
    constexpr float minTPCCrossedRowsOverFound = 0.f;
    constexpr int maxTPCNClsShared = 160;
    constexpr float maxSharedTPCFraction = 1.f;
    constexpr int minITSNClsInnerBarrel = 1;
    const float tpcCrossedRowsOverFound = candidate.tpcNClsFound() > 0 ? static_cast<float>(candidate.tpcNClsCrossedRows()) / candidate.tpcNClsFound() : 0.f;

    return !(candidate.tpcNClsFound() < minTPCNClsFound ||
             candidate.tpcNClsCrossedRows() < minTPCNClsCrossedRows ||
             tpcCrossedRowsOverFound < minTPCCrossedRowsOverFound ||
             candidate.tpcNClsShared() > maxTPCNClsShared ||
             candidate.tpcFractionSharedCls() > maxSharedTPCFraction ||
             candidate.itsNCls() < trackCut.settingCutNCls ||
             candidate.itsNClsInnerBarrel() < minITSNClsInnerBarrel);
  }

  template <typename Ttrack>
  bool selectTrackTr(const Ttrack& candidate)
  {
    constexpr float maxAbsEta = 0.8f;
    constexpr float maxAbsDcaXY = 0.2f;
    constexpr float maxAbsDcaZ = 0.2f;
    constexpr int minTPCCrossedRows = 70;
    constexpr float maxTPCChi2NCl = 5.f;
    constexpr float maxTPCFractionSharedCls = 0.3f;
    constexpr int minITSNCls = 5;
    constexpr float maxITSChi2NCl = 10.f;

    return !(std::abs(candidate.eta()) >= maxAbsEta ||
             std::abs(candidate.dcaXY()) >= maxAbsDcaXY ||
             std::abs(candidate.dcaZ()) >= maxAbsDcaZ ||
             candidate.tpcNClsCrossedRows() < minTPCCrossedRows ||
             candidate.tpcChi2NCl() >= maxTPCChi2NCl ||
             candidate.tpcFractionSharedCls() >= maxTPCFractionSharedCls ||
             candidate.itsNCls() < minITSNCls ||
             candidate.itsChi2NCl() >= maxITSChi2NCl);
  }

  bool useDeuteronNucleus() const
  {
    return species.settingNuPDGCode.value == DeuteronPDG;
  }

  bool useTritonNucleus() const
  {
    return species.settingNuPDGCode.value == TritonPDG;
  }

  bool useHelium3Nucleus() const
  {
    return species.settingNuPDGCode.value == He3PDG;
  }

  const char* zorroTriggerMask() const
  {
    if (useTritonNucleus()) {
      return "fTritonFemto";
    }
    return useHelium3Nucleus() ? "fHe" : "fDe";
  }

  const char* nucleusBetheBlochLabel() const
  {
    return useHelium3Nucleus() ? "He3" : "De";
  }

  float nucleusChargeFactor() const
  {
    return useHelium3Nucleus() ? 2.f : 1.f;
  }

  float nucleusMass() const
  {
    if (useHelium3Nucleus()) {
      return static_cast<float>(o2::constants::physics::MassHelium3);
    }
    if (useTritonNucleus()) {
      return static_cast<float>(o2::constants::physics::MassTriton);
    }
    return static_cast<float>(o2::constants::physics::MassDeuteron);
  }

  template <typename Ttrack>
  bool selectTrackHe3(const Ttrack& candidate)
  {
    if (std::abs(candidate.eta()) > trackCut.settingCutEta) {
      return false;
    }

    constexpr float crossedRowsToFindableRatio = 0.8f;
    return !(candidate.itsNCls() < trackCut.settingCutNCls ||
             candidate.tpcNClsFound() < helium3Pid.settingTPCNClsFoundMinHe3 ||
             candidate.tpcNClsCrossedRows() < helium3Pid.settingTPCCrossedRowsMinHe3 ||
             candidate.tpcNClsCrossedRows() < crossedRowsToFindableRatio * candidate.tpcNClsFindable() ||
             candidate.tpcChi2NCl() > trackCut.settingCutChi2tpcHigh ||
             candidate.tpcChi2NCl() < He3TPCChi2NClMin ||
             candidate.itsChi2NCl() > trackCut.settingCutChi2NClITS);
  }

  template <typename Ttrack>
  bool selectTrackNu(const Ttrack& candidate)
  {
    if (useHelium3Nucleus()) {
      return selectTrackHe3(candidate);
    }
    if (useTritonNucleus()) {
      return selectTrackTr(candidate);
    }
    if (useDeuteronNucleus()) {
      return selectTrackDe(candidate);
    }
    LOG(info) << "invalid nucleus PDG code";
    return false;
  }

  void fillNucleusTrackSelection(const Selections selection)
  {
    mQaRegistry.fill(HIST("hTrackSelNu"), selection);
  }

  void fillNucleusPairFlow(const int step)
  {
    mQaRegistry.fill(HIST("hNuPairFlow"), step);
  }

  bool isPionDeuteronPair() const
  {
    return species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus) && useDeuteronNucleus();
  }

  bool isPionTritonPair() const
  {
    return species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus) && useTritonNucleus();
  }

  bool shouldFillOutputTable(const HadNucandidate& hadNucand) const
  {
    if (!output.settingFillTableLowKstarOnly.value) {
      return true;
    }
    if (isPionDeuteronPair()) {
      return hadNucand.kstar < output.settingFillTablePiDeKstarMax.value;
    }
    if (isPionTritonPair()) {
      return hadNucand.kstar < output.settingFillTablePiTrKstarMax.value;
    }
    return true;
  }

  template <typename Ttrack>
  float phiAtSpecificRadiiTPC(const Ttrack& track, float radius) const
  {
    const float absPt = std::abs(track.pt());
    if (absPt <= 0.f) {
      return 999.f;
    }
    const float arg = 0.3f * static_cast<float>(track.sign()) * 0.1f * mDbz * radius * 0.01f / (2.f * absPt);
    if (std::fabs(arg) >= 1.f) {
      return 999.f;
    }
    return track.phi() - std::asin(arg);
  }

  float wrapDeltaPhi(float dphi) const
  {
    return std::atan2(std::sin(dphi), std::cos(dphi));
  }

  template <typename Ttrack1, typename Ttrack2>
  float averagePhiStar(const Ttrack1& firstTrack, const Ttrack2& secondTrack) const
  {
    constexpr float invalidPhiStar = 999.f;
    float dPhiAvg = 0.f;
    int meaningfulEntries = 0;
    for (const auto& radius : tmpRadiiTPC) {
      const float phi1 = phiAtSpecificRadiiTPC(firstTrack, radius);
      const float phi2 = phiAtSpecificRadiiTPC(secondTrack, radius);
      if (phi1 == invalidPhiStar || phi2 == invalidPhiStar) {
        continue;
      }
      dPhiAvg += wrapDeltaPhi(phi1 - phi2);
      meaningfulEntries++;
    }
    if (meaningfulEntries == 0) {
      return invalidPhiStar;
    }
    return dPhiAvg / static_cast<float>(meaningfulEntries);
  }

  template <typename Ttrack1, typename Ttrack2>
  bool computeClosePairDeltas(const Ttrack1& firstTrack, const Ttrack2& secondTrack, float& deltaEta, float& deltaPhi) const
  {
    constexpr int closePairRadiusModePv = 0;
    constexpr int closePairRadiusModeSpecificTpc = 2;
    constexpr float invalidPhiStar = 999.f;
    constexpr float invalidOutputDelta = -99.f;

    deltaEta = invalidOutputDelta;
    deltaPhi = invalidOutputDelta;
    if (!CPR.settingEnableClosePairRejection.value) {
      return false;
    }

    float selectedDeltaPhi = averagePhiStar(firstTrack, secondTrack);
    if (CPR.settingClosePairRadiusMode.value == closePairRadiusModePv) {
      selectedDeltaPhi = wrapDeltaPhi(firstTrack.phi() - secondTrack.phi());
    } else if (CPR.settingClosePairRadiusMode.value == closePairRadiusModeSpecificTpc) {
      const float firstPhi = phiAtSpecificRadiiTPC(firstTrack, CPR.settingClosePairSpecificRadius.value);
      const float secondPhi = phiAtSpecificRadiiTPC(secondTrack, CPR.settingClosePairSpecificRadius.value);
      if (firstPhi == invalidPhiStar || secondPhi == invalidPhiStar) {
        return false;
      }
      selectedDeltaPhi = wrapDeltaPhi(firstPhi - secondPhi);
    }

    if (selectedDeltaPhi == invalidPhiStar) {
      return false;
    }

    deltaEta = firstTrack.eta() - secondTrack.eta();
    deltaPhi = selectedDeltaPhi;
    return true;
  }

  template <typename Ttrack1, typename Ttrack2>
  bool isClosePair(const Ttrack1& firstTrack, const Ttrack2& secondTrack, bool fillQA)
  {
    if (!CPR.settingEnableClosePairRejection.value) {
      return false;
    }
    if (firstTrack.sign() != secondTrack.sign()) {
      return false;
    }

    float deltaEta = -99.f;
    float deltaPhi = -99.f;
    if (!computeClosePairDeltas(firstTrack, secondTrack, deltaEta, deltaPhi)) {
      return false;
    }

    if (fillQA) {
      mQaRegistry.fill(HIST("h2CPRBefore"), deltaEta, deltaPhi);
    }
    const bool isRejected = std::pow(deltaPhi, 2.f) / std::pow(CPR.settingClosePairDeltaPhiMax.value, 2.f) +
                              std::pow(deltaEta, 2.f) / std::pow(CPR.settingClosePairDeltaEtaMax.value, 2.f) <
                            1.f;
    if (fillQA && !isRejected) {
      mQaRegistry.fill(HIST("h2CPRAfter"), deltaEta, deltaPhi);
    }
    return isRejected;
  }

  bool isCloseHadHyperPairAtPV(const HadHyperCandidate& candidate, const HadHyperHadron& hadron, bool fillQA)
  {
    if (!CPR.settingEnableClosePairRejection.value) {
      return false;
    }
    constexpr int closePairRadiusModePv = 0;
    if (CPR.settingClosePairRadiusMode.value != closePairRadiusModePv) {
      return false;
    }
    if (candidate.sign() != hadron.sign()) {
      return false;
    }

    const float deltaEta = candidate.eta() - hadron.eta();
    const float deltaPhi = wrapDeltaPhi(candidate.phi() - hadron.phi());
    if (fillQA) {
      hadHyperRegistry.fill(HIST("hClosePairBeforePV"), deltaEta, deltaPhi);
    }
    const bool isRejected = std::pow(deltaPhi, 2.f) / std::pow(CPR.settingClosePairDeltaPhiMax.value, 2.f) +
                              std::pow(deltaEta, 2.f) / std::pow(CPR.settingClosePairDeltaEtaMax.value, 2.f) <
                            1.f;
    if (fillQA && !isRejected) {
      hadHyperRegistry.fill(HIST("hClosePairAfterPV"), deltaEta, deltaPhi);
    }
    return isRejected;
  }

  template <typename Ttrack>
  bool selectionPIDProton(const Ttrack& candidate)
  {
    constexpr float protonPtMin = 0.5f;
    constexpr float protonPtMax = 3.0f;
    constexpr float protonPCombMin = 0.75f;
    constexpr float protonTPCNsigmaMax = 3.0f;
    constexpr float protonCombNsigmaMax = 3.0f;

    const float tpcNSigmaPr = candidate.tpcNSigmaPr();
    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.sign() * candidate.tpcInnerParam(), tpcNSigmaPr);

    if (std::abs(candidate.pt()) <= protonPtMin || std::abs(candidate.pt()) >= protonPtMax) {
      return false;
    }

    const float absPin = std::abs(candidate.tpcInnerParam());
    if (absPin < protonPCombMin) {
      if (std::abs(tpcNSigmaPr) > protonTPCNsigmaMax) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPr);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    const float tofNSigmaPr = candidate.tofNSigmaPr();
    const float combNsigma = std::sqrt(tpcNSigmaPr * tpcNSigmaPr + tofNSigmaPr * tofNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.sign() * candidate.pt(), tofNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadComb_preselection"), candidate.sign() * candidate.pt(), combNsigma);
    if (combNsigma > protonCombNsigmaMax) {
      return false;
    }
    if (deuteronPid.settingReqSingleNsig.value &&
        (std::abs(tpcNSigmaPr) > protonCombNsigmaMax || std::abs(tofNSigmaPr) > protonCombNsigmaMax)) {
      return false;
    }

    mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaPr);
    mQaRegistry.fill(HIST("h2NsigmaHadComb"), candidate.sign() * candidate.pt(), combNsigma);
    mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    return true;
  }

  template <typename Ttrack>
  bool selectionPIDKaon(const Ttrack& candidate)
  {
    auto tpcNSigmaKa = candidate.tpcNSigmaKa();
    float DeDCAxyMin = 0.004 + (0.013 / candidate.pt());
    float DeDCAzMin = 0.004 + (0.013 / candidate.pt());
    if (std::abs(candidate.dcaXY()) > DeDCAxyMin || std::abs(candidate.dcaZ()) > DeDCAzMin) {
      return false;
    }

    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.tpcInnerParam(), tpcNSigmaKa);
    if (std::abs(candidate.pt()) < hadronPid.settingHadptMin || std::abs(candidate.pt()) > hadronPid.settingHadptMax) {
      return false;
    }

    // reject protons and pions
    if (std::abs(candidate.tpcNSigmaPr()) < hadronPid.settingCutNsigTPCPrMin || std::abs(candidate.tpcNSigmaPi()) < hadronPid.settingCutNsigTPCPiMin) {
      return false;
    }
    mQaRegistry.fill(HIST("h2NsigmaHadPrTPC"), candidate.tpcNSigmaPr());
    mQaRegistry.fill(HIST("h2NsigmaHadPiTPC"), candidate.tpcNSigmaPi());
    if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < hadronPid.settingCutNsigTOFPrMin) {
      return false;
    }
    if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < hadronPid.settingCutNsigTOFPiMin) {
      return false;
    }
    mQaRegistry.fill(HIST("h2NsigmaHadPrTOF"), candidate.tofNSigmaPr());
    mQaRegistry.fill(HIST("h2NsigmaHadPiTOF"), candidate.tofNSigmaPi());
    // rejection end

    if (candidate.hasTOF() && candidate.tpcInnerParam() >= hadronPid.settingCutPinMinTOFHad) {
      auto tofNSigmaKa = candidate.tofNSigmaKa();

      mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.pt(), tofNSigmaKa);
      if (std::abs(tofNSigmaKa) > hadronPid.settingCutNsigmaTOFHad) {
        return false;
      }
      if (std::abs(tpcNSigmaKa) > hadronPid.settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaKa);
      mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaKa);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }
    if (candidate.tpcInnerParam() < hadronPid.settingCutPinMinTOFHad) {
      if (std::abs(tpcNSigmaKa) > hadronPid.settingCutNsigmaTPCHad) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaKa);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  bool selectionPIDPion(const Ttrack& candidate)
  {
    const float tpcNSigmaPi = candidate.tpcNSigmaPi();
    const float absP = std::abs(candidate.p());
    mQaRegistry.fill(HIST("h2NsigmaHadTPC_preselection"), candidate.sign() * candidate.tpcInnerParam(), tpcNSigmaPi);

    if (absP <= hadronPid.settingPionMomCombMin) {
      if (std::abs(tpcNSigmaPi) > hadronPid.settingPionTPCNsigMax) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
      mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
      return true;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    const float tofNSigmaPi = candidate.tofNSigmaPi();
    const float combNsigma = std::sqrt(tofNSigmaPi * tofNSigmaPi + tpcNSigmaPi * tpcNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF_preselection"), candidate.sign() * candidate.pt(), tofNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadComb_preselection"), candidate.sign() * candidate.pt(), combNsigma);
    if (combNsigma > hadronPid.settingPionCombNsigMax) {
      return false;
    }
    if (deuteronPid.settingReqSingleNsig.value &&
        (std::abs(tpcNSigmaPi) > hadronPid.settingPionCombNsigMax || std::abs(tofNSigmaPi) > hadronPid.settingPionCombNsigMax)) {
      return false;
    }

    mQaRegistry.fill(HIST("h2NsigmaHadTPC"), candidate.sign() * candidate.pt(), tpcNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadTOF"), candidate.sign() * candidate.pt(), tofNSigmaPi);
    mQaRegistry.fill(HIST("h2NsigmaHadComb"), candidate.sign() * candidate.pt(), combNsigma);
    mQaRegistry.fill(HIST("h2dEdxHadcandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    return true;
  }

  template <typename Ttrack>
  bool selectionPIDHadron(const Ttrack& candidate)
  {
    bool PID = false;
    if (species.settingHadPDGCode == PDG_t::kPiPlus) {
      PID = selectionPIDPion(candidate);
      mMassHad = o2::constants::physics::MassPiPlus;
    } else if (species.settingHadPDGCode == PDG_t::kKPlus) {
      PID = selectionPIDKaon(candidate);
      mMassHad = o2::constants::physics::MassKPlus;
    } else if (species.settingHadPDGCode == PDG_t::kProton) {
      PID = selectionPIDProton(candidate);
      mMassHad = o2::constants::physics::MassProton;
    } else {
      LOG(info) << "invalid PDG code";
    }
    return PID;
  }

  template <typename Ttrack>
  float getHadronTPCNSigma(const Ttrack& candidate) const
  {
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return candidate.tpcNSigmaPi();
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kKPlus)) {
      return candidate.tpcNSigmaKa();
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return candidate.tpcNSigmaPr();
    }
    return -10.f;
  }

  template <typename Ttrack>
  float getHadronTOFNSigma(const Ttrack& candidate) const
  {
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return candidate.tofNSigmaPi();
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kKPlus)) {
      return candidate.tofNSigmaKa();
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return candidate.tofNSigmaPr();
    }
    return -10.f;
  }

  template <typename Ttrack>
  float getHadronITSNSigma(const Ttrack& candidate) const
  {
    o2::aod::ITSResponse itsResponse;
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return itsResponse.nSigmaITS<o2::track::PID::Pion>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kKPlus)) {
      return itsResponse.nSigmaITS<o2::track::PID::Kaon>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return itsResponse.nSigmaITS<o2::track::PID::Proton>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    }
    return -10.f;
  }

  template <typename Ttrack>
  float computeNSigmaDe(const Ttrack& candidate)
  {
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() / constants::physics::MassDeuteron), mBBparamsNucleus[0], mBBparamsNucleus[1], mBBparamsNucleus[2], mBBparamsNucleus[3], mBBparamsNucleus[4]);
    double resoTPC{expTPCSignal * mBBparamsNucleus[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  float correctedTPCInnerParamHe3(const Ttrack& candidate) const
  {
    const bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    return (heliumPID && pidCalibration.settingCompensatePIDinTracking.value) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
  }

  template <typename Ttrack>
  float computeNSigmaHe3(const Ttrack& candidate)
  {
    const float correctedTPCinnerParam = correctedTPCInnerParamHe3(candidate);
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), mBBparamsNucleus[0], mBBparamsNucleus[1], mBBparamsNucleus[2], mBBparamsNucleus[3], mBBparamsNucleus[4]);
    double resoTPC{expTPCSignal * mBBparamsNucleus[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename Ttrack>
  bool selectionPIDDe(const Ttrack& candidate)
  {
    float tpcInnerParam = candidate.tpcInnerParam();
    mQaRegistry.fill(HIST("h2dEdx"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());

    if (std::abs(tpcInnerParam) < deuteronPid.settingCutPinMinDe) {
      return false;
    }
    float tpcNSigmaDe = 0.f;
    if (output.settingUseBBcomputeDeNsigma) {
      tpcNSigmaDe = computeNSigmaDe(candidate);
    } else {
      tpcNSigmaDe = candidate.tpcNSigmaDe();
    }

    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselection"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselecComp"), candidate.sign() * candidate.pt(), candidate.tpcNSigmaDe());
    if (std::abs(candidate.pt()) < deuteronPid.settingCutDeptMin || std::abs(candidate.pt()) > deuteronPid.settingCutDeptMax) {
      return false;
    }
    const float absPt = std::abs(candidate.pt());
    if (absPt <= 0.f) {
      return false;
    }
    const float deDCAxyMax = 0.004f + 0.013f / absPt;
    const float deDCAzMax = 0.004f + 0.013f / absPt;
    if (std::abs(candidate.dcaXY()) > deDCAxyMax || std::abs(candidate.dcaZ()) > deDCAzMax) {
      return false;
    }

    if (candidate.hasTOF() && candidate.tpcInnerParam() > deuteronPid.settingCutPinMinTOFITSDe) {
      auto tofNSigmaDe = candidate.tofNSigmaDe();
      auto combNsigma = std::sqrt(tofNSigmaDe * tofNSigmaDe + tpcNSigmaDe * tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NsigmaNuTOF_preselection"), candidate.sign() * candidate.pt(), tofNSigmaDe);
      if (combNsigma > deuteronPid.settingCutNsigmaTOFTPCDe) {
        return false;
      }
      if (deuteronPid.settingReqSingleNsig.value &&
          (std::abs(tpcNSigmaDe) > deuteronPid.settingCutNsigmaTOFTPCDe || std::abs(tofNSigmaDe) > deuteronPid.settingCutNsigmaTOFTPCDe)) {
        return false;
      }
      mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());
      mQaRegistry.fill(HIST("h2NsigmaNuComb"), candidate.sign() * candidate.pt(), combNsigma);
      mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NsigmaNuTOF"), candidate.sign() * candidate.pt(), tofNSigmaDe);
      return true;
    }
    if (candidate.tpcInnerParam() <= deuteronPid.settingCutPinMinTOFITSDe) {
      if (std::abs(tpcNSigmaDe) > deuteronPid.settingCutNsigmaTPCDe) {
        return false;
      }
      o2::aod::ITSResponse itsResponse;
      auto itsnSigmaDe = itsResponse.nSigmaITS<o2::track::PID::Deuteron>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
      mQaRegistry.fill(HIST("h2NSigmaNuITS_preselection"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      if (std::abs(itsnSigmaDe) > deuteronPid.settingCutNsigmaITSDe) {
        return false;
      }
      mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * candidate.pt(), tpcNSigmaDe);
      mQaRegistry.fill(HIST("h2NSigmaNuITS"), candidate.sign() * candidate.pt(), itsnSigmaDe);
      mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * tpcInnerParam, candidate.tpcSignal());
      return true;
    }
    return false;
  }

  template <typename Ttrack>
  bool selectionPIDHe3(const Ttrack& candidate)
  {
    const float correctedTPCinnerParam = correctedTPCInnerParamHe3(candidate);
    mQaRegistry.fill(HIST("h2dEdx"), candidate.sign() * correctedTPCinnerParam, candidate.tpcSignal());

    if (correctedTPCinnerParam < helium3Pid.settingRigidityMinHe3) {
      return false;
    }

    const float nSigmaHe3 = computeNSigmaHe3(candidate);
    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselection"), candidate.sign() * 2.f * candidate.pt(), nSigmaHe3);
    if (std::abs(nSigmaHe3) > helium3Pid.settingTPCNSigmaMaxHe3) {
      return false;
    }

    o2::aod::ITSResponse itsResponse;
    const float itsNsigmaHe3 = itsResponse.nSigmaITS<o2::track::PID::Helium3>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    mQaRegistry.fill(HIST("h2NSigmaNuITS_preselection"), candidate.sign() * 2.f * candidate.pt(), itsNsigmaHe3);
    if (itsNsigmaHe3 < helium3Pid.settingITSNSigmaMinHe3) {
      return false;
    }

    mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * correctedTPCinnerParam, candidate.tpcSignal());
    mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * 2.f * candidate.pt(), nSigmaHe3);
    mQaRegistry.fill(HIST("h2NSigmaNuITS"), candidate.sign() * 2.f * candidate.pt(), itsNsigmaHe3);
    return true;
  }

  template <typename Ttrack>
  float computeTOFMass(const Ttrack& candidate) const
  {
    if (!candidate.hasTOF()) {
      return -1.f;
    }
    float beta = o2::pid::tof::Beta::GetBeta(candidate);
    beta = std::clamp(beta, 1.e-4f, 1.f - 1.e-6f);
    return std::abs(candidate.tpcInnerParam()) * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  template <typename Ttrack>
  bool selectionPIDTr(const Ttrack& candidate)
  {
    const float tpcNSigmaTr = candidate.tpcNSigmaTr();
    const float absP = std::abs(candidate.p());
    o2::aod::ITSResponse itsResponse;
    const float itsNSigmaTr = itsResponse.nSigmaITS<o2::track::PID::Triton>(candidate.itsClusterSizes(), absP, candidate.eta());

    mQaRegistry.fill(HIST("h2dEdx"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    mQaRegistry.fill(HIST("h2NsigmaNuTPC_preselection"), candidate.sign() * candidate.pt(), tpcNSigmaTr);
    mQaRegistry.fill(HIST("h2NSigmaNuITS_preselection"), candidate.sign() * candidate.pt(), itsNSigmaTr);

    if (std::abs(candidate.tpcNSigmaDe()) < tritonPid.settingTPCRejectNsig ||
        std::abs(candidate.tpcNSigmaPr()) < tritonPid.settingTPCRejectNsig ||
        std::abs(candidate.tpcNSigmaPi()) < tritonPid.settingTPCRejectNsig) {
      return false;
    }

    if (absP < tritonPid.settingPIDMomentumSplitTr) {
      if (std::abs(tpcNSigmaTr) >= tritonPid.settingCutTPCNsigmaLowPTr ||
          std::abs(itsNSigmaTr) >= tritonPid.settingCutITSNsigmaLowPTr) {
        return false;
      }
    } else if (tpcNSigmaTr <= tritonPid.settingCutTPCNsigmaHighPMinTr ||
               tpcNSigmaTr >= tritonPid.settingCutTPCNsigmaHighPMaxTr) {
      return false;
    }

    if (absP > tritonPid.settingTOFMassMomentumMinTr) {
      if (!candidate.hasTOF()) {
        return false;
      }
      const float massTOFTr = computeTOFMass(candidate);
      mQaRegistry.fill(HIST("h2MassTOFTr"), candidate.sign() * absP, massTOFTr);
      if (massTOFTr <= tritonPid.settingTOFMassMinTr || massTOFTr >= tritonPid.settingTOFMassMaxTr) {
        return false;
      }
    }

    mQaRegistry.fill(HIST("h2dEdxNucandidates"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
    mQaRegistry.fill(HIST("h2NsigmaNuTPC"), candidate.sign() * candidate.pt(), tpcNSigmaTr);
    mQaRegistry.fill(HIST("h2NSigmaNuITS"), candidate.sign() * candidate.pt(), itsNSigmaTr);
    return true;
  }

  template <typename Ttrack>
  bool selectionPIDNu(const Ttrack& candidate)
  {
    if (useHelium3Nucleus()) {
      return selectionPIDHe3(candidate);
    }
    if (useTritonNucleus()) {
      return selectionPIDTr(candidate);
    }
    if (useDeuteronNucleus()) {
      return selectionPIDDe(candidate);
    }
    return false;
  }

  template <typename Ttrack>
  float getNucleusTPCNSigma(const Ttrack& candidate)
  {
    if (useHelium3Nucleus()) {
      return computeNSigmaHe3(candidate);
    }
    if (useTritonNucleus()) {
      return candidate.tpcNSigmaTr();
    }
    return computeNSigmaDe(candidate);
  }

  template <typename Ttrack>
  float getNucleusTOFNSigma(const Ttrack& candidate) const
  {
    if (useHelium3Nucleus()) {
      return candidate.tofNSigmaHe();
    }
    if (useTritonNucleus()) {
      return candidate.tofNSigmaTr();
    }
    if (useDeuteronNucleus()) {
      return candidate.tofNSigmaDe();
    }
    return -10.f;
  }

  template <typename Ttrack>
  float getNucleusITSNSigma(const Ttrack& candidate) const
  {
    o2::aod::ITSResponse itsResponse;
    if (useHelium3Nucleus()) {
      return itsResponse.nSigmaITS<o2::track::PID::Helium3>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    }
    if (useTritonNucleus()) {
      return itsResponse.nSigmaITS<o2::track::PID::Triton>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    }
    if (useDeuteronNucleus()) {
      return itsResponse.nSigmaITS<o2::track::PID::Deuteron>(candidate.itsClusterSizes(), candidate.p(), candidate.eta());
    }
    return -10.f;
  }

  template <typename Tcandidate>
  float computeHyperCandidateMass(const Tcandidate& candidate) const
  {
    const std::array<float, 3> heMomentum{
      candidate.ptHe3() * std::cos(candidate.phiHe3()),
      candidate.ptHe3() * std::sin(candidate.phiHe3()),
      candidate.ptHe3() * std::sinh(candidate.etaHe3())};
    const std::array<float, 3> pionMomentum{
      candidate.ptPi() * std::cos(candidate.phiPi()),
      candidate.ptPi() * std::sin(candidate.phiPi()),
      candidate.ptPi() * std::sinh(candidate.etaPi())};
    return RecoDecay::m(std::array{heMomentum, pionMomentum},
                        std::array{static_cast<float>(o2::constants::physics::MassHelium3),
                                   static_cast<float>(o2::constants::physics::MassPiPlus)});
  }

  template <typename Tcandidate>
  bool selectHyperCandidate(const Tcandidate& candidate)
  {
    mQaRegistry.fill(HIST("hHe3P_preselected"), candidate.tpcMomHe());
    const float mass = computeHyperCandidateMass(candidate);
    if (!std::isfinite(mass) || mass < hypertriton.settingHypMassMin || mass > hypertriton.settingHypMassMax) {
      return false;
    }
    mQaRegistry.fill(HIST("hHe3P"), candidate.tpcMomHe());
    mQaRegistry.fill(HIST("hHe3TPCnsigma"), candidate.ptHe3(), candidate.nSigmaHe());

    return true;
  }

  float computePairKstar(const std::array<float, 3>& momHad, const float massHad, const std::array<float, 3>& momNu, const float massNu) const
  {
    const PairLorentzVector vecHad(momHad[0], momHad[1], momHad[2], massHad);
    const PairLorentzVector vecNu(momNu[0], momNu[1], momNu[2], massNu);
    const PairLorentzVector trackSum = vecHad + vecNu;

    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());

    PairLorentzVector partHadCMS(vecHad);
    PairLorentzVector partNuCMS(vecNu);

    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    partHadCMS = boostPRF(partHadCMS);
    partNuCMS = boostPRF(partNuCMS);

    const PairLorentzVector trackRelK = partHadCMS - partNuCMS;
    return 0.5f * trackRelK.P();
  }

  float computePairMT(const std::array<float, 3>& momHad, const float massHad, const std::array<float, 3>& momNu, const float massNu) const
  {
    const PairLorentzVector vecHad(momHad[0], momHad[1], momHad[2], massHad);
    const PairLorentzVector vecNu(momNu[0], momNu[1], momNu[2], massNu);
    const PairLorentzVector trackSum = vecHad + vecNu;
    const float kT = 0.5f * trackSum.Pt();
    return std::sqrt(kT * kT + std::pow(0.5f * (massHad + massNu), 2.f));
  }

  // ==================================================================================================================

  template <typename Ttrack, typename Tcollisions, typename Ttracks>
  bool fillCandidateInfo(const Ttrack& trackDe, const Ttrack& trackHad, const CollBracket& collBracket, const Tcollisions& collisions, HadNucandidate& hadNucand, const Ttracks& /*trackTable*/, bool isMixedEvent)
  {
    const int numCoordinates = 3;
    if (!isMixedEvent) {
      auto trackCovDe = getTrackParCov(trackDe);
      auto trackCovHad = getTrackParCov(trackHad);
      int nCand = 0;
      try {
        nCand = mFitter.process(trackCovDe, trackCovHad);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        mQaRegistry.fill(HIST("hSkipReasons"), 0);
        return false;
      }
      if (nCand == 0) {
        mQaRegistry.fill(HIST("hSkipReasons"), 1);
        return false;
      }
      hadNucand.dcaPair = std::sqrt(std::abs(mFitter.getChi2AtPCACandidate()));

      // associate collision id as the one that minimises the distance between the vertex and the PCAs of the daughters
      double distanceMin = -1;
      unsigned int collIdxMin = 0;

      for (int collIdx = collBracket.getMin(); collIdx <= collBracket.getMax(); collIdx++) {
        auto collision = collisions.rawIteratorAt(collIdx);
        std::array<float, 3> collVtx = {collision.posX(), collision.posY(), collision.posZ()};
        const auto& pca = mFitter.getPCACandidate();
        float distance = 0;
        for (int i = 0; i < numCoordinates; i++) {
          distance += (pca[i] - collVtx[i]) * (pca[i] - collVtx[i]);
        }
        if (distanceMin < 0 || distance < distanceMin) {
          distanceMin = distance;
          collIdxMin = collIdx;
        }
      }

      if (!mGoodCollisions[collIdxMin]) {
        mQaRegistry.fill(HIST("hSkipReasons"), 2);
        return false;
      }
      hadNucand.collisionID = collIdxMin;
    } else {
      hadNucand.collisionID = collBracket.getMin();
    }

    const float nuChargeFactor = nucleusChargeFactor();
    hadNucand.momNu = std::array{trackDe.px(), trackDe.py(), trackDe.pz()};
    for (auto i = 0u; i < hadNucand.momNu.size(); ++i) {
      hadNucand.momNu[i] *= nuChargeFactor;
    }
    hadNucand.momHad = std::array{trackHad.px(), trackHad.py(), trackHad.pz()};
    float invMass = 0;
    invMass = RecoDecay::m(std::array<std::array<float, 3>, 2>{hadNucand.momNu, hadNucand.momHad}, std::array<float, 2>{nucleusMass(), mMassHad});

    hadNucand.signNu = trackDe.sign();
    hadNucand.signHad = trackHad.sign();
    computeClosePairDeltas(trackDe, trackHad, hadNucand.deltaEta, hadNucand.deltaPhi);

    hadNucand.dcaxyNu = trackDe.dcaXY();
    hadNucand.dcaxyHad = trackHad.dcaXY();

    hadNucand.dcazNu = trackDe.dcaZ();
    hadNucand.dcazHad = trackHad.dcaZ();

    hadNucand.tpcSignalNu = trackDe.tpcSignal();
    hadNucand.momNuTPC = useHelium3Nucleus() ? correctedTPCInnerParamHe3(trackDe) : trackDe.tpcInnerParam();
    hadNucand.tpcSignalHad = trackHad.tpcSignal();
    hadNucand.momHadTPC = trackHad.tpcInnerParam();

    hadNucand.nTPCClustersNu = trackDe.tpcNClsFound();
    hadNucand.nTPCClustersHad = trackHad.tpcNClsFound();
    hadNucand.nTPCCrossedRowsNu = trackDe.tpcNClsCrossedRows();
    hadNucand.nTPCCrossedRowsHad = trackHad.tpcNClsCrossedRows();
    hadNucand.nSigmaNu = getNucleusTPCNSigma(trackDe);
    hadNucand.nSigmaHad = getHadronTPCNSigma(trackHad);
    hadNucand.nSigmaTOFNu = getNucleusTOFNSigma(trackDe);
    hadNucand.nSigmaITSNu = getNucleusITSNSigma(trackDe);
    hadNucand.nSigmaTOFHad = getHadronTOFNSigma(trackHad);
    hadNucand.nSigmaITSHad = getHadronITSNSigma(trackHad);
    hadNucand.nSigmaTPCHadPi = trackHad.tpcNSigmaPi();
    hadNucand.nSigmaTPCHadKa = trackHad.tpcNSigmaKa();
    hadNucand.nSigmaTPCHadPr = trackHad.tpcNSigmaPr();
    hadNucand.nSigmaTOFHadPi = trackHad.tofNSigmaPi();
    hadNucand.nSigmaTOFHadKa = trackHad.tofNSigmaKa();
    hadNucand.nSigmaTOFHadPr = trackHad.tofNSigmaPr();

    hadNucand.chi2TPCNu = trackDe.tpcChi2NCl();
    hadNucand.chi2TPCHad = trackHad.tpcChi2NCl();

    hadNucand.pidTrkNu = trackDe.pidForTracking();
    hadNucand.pidTrkHad = trackHad.pidForTracking();

    hadNucand.itsClSizeNu = trackDe.itsClusterSizes();
    hadNucand.itsClSizeHad = trackHad.itsClusterSizes();

    hadNucand.nClsItsNu = trackDe.itsNCls();
    hadNucand.nClsItsHad = trackHad.itsNCls();

    hadNucand.sharedClustersNu = trackDe.tpcNClsShared();
    hadNucand.sharedClustersHad = trackHad.tpcNClsShared();

    hadNucand.isBkgUS = trackDe.sign() * trackHad.sign() < 0;
    hadNucand.isBkgEM = isMixedEvent;

    hadNucand.invMass = invMass;

    hadNucand.trackIDNu = trackDe.globalIndex();
    hadNucand.trackIDHad = trackHad.globalIndex();

    if (trackDe.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackDe);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      hadNucand.massTOFNu = hadNucand.momNuTPC * nuChargeFactor * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    if (trackHad.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(trackHad);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      hadNucand.massTOFHad = trackHad.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
    }

    float massLightNucleusForKstarMt = nucleusMass();
    if (useDeuteronNucleus() && deuteronPid.settingUseProtonMassForKstarMt) {
      massLightNucleusForKstarMt = static_cast<float>(o2::constants::physics::MassProton);
    }
    hadNucand.kstar = computePairKstar(hadNucand.momHad, mMassHad, hadNucand.momNu, massLightNucleusForKstarMt);
    hadNucand.mT = computePairMT(hadNucand.momHad, mMassHad, hadNucand.momNu, massLightNucleusForKstarMt);

    return true;
  }

  template <typename Ttrack>
  void pairTracksSameEvent(const Ttrack& tracks, float /*cent*/)
  {
    // LOG(info) << "Number of tracks: " << tracks.size();
    for (const auto& track0 : tracks) {

      mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
      fillNucleusTrackSelection(Selections::kNoCuts);

      if (!selectTrackNu(track0)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);
      fillNucleusTrackSelection(Selections::kTrackCuts);

      if (!selectionPIDNu(track0)) {
        continue;
      }
      mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);
      fillNucleusTrackSelection(Selections::kPID);
      mQaRegistry.fill(HIST("hSingleNuPt"), track0.pt() * track0.sign() * nucleusChargeFactor());
      mQaRegistry.fill(HIST("hSingleNuPin"), (useHelium3Nucleus() ? correctedTPCInnerParamHe3(track0) : track0.tpcInnerParam()) * track0.sign());
      fillNucleusPairFlow(0);

      bool hasHadronSelected = false;
      bool hasStoredPair = false;

      for (const auto& track1 : tracks) {
        if (track0 == track1) {
          continue;
        }

        if (!eventMixing.settingSaveUSandLS) {
          if (!eventMixing.settingEnableBkgUS && (track0.sign() * track1.sign() < 0)) {
            continue;
          }
          if (eventMixing.settingEnableBkgUS && (track0.sign() * track1.sign() > 0)) {
            continue;
          }
        }

        if (!selectTrackHadron(track1) || !selectionPIDHadron(track1)) {
          continue;
        }
        hasHadronSelected = true;
        if (isClosePair(track0, track1, /*fillQA*/ true)) {
          continue;
        }

        SVCand trackPair;
        trackPair.tr0Idx = track0.globalIndex();
        trackPair.tr1Idx = track1.globalIndex();
        const int collIdx = track0.collisionId();
        CollBracket collBracket{collIdx, collIdx};
        trackPair.collBracket = collBracket;
        mTrackPairs.push_back(trackPair);
        hasStoredPair = true;
      }

      if (hasHadronSelected) {
        fillNucleusPairFlow(1);
      }
      if (hasStoredPair) {
        fillNucleusPairFlow(2);
      }
    }
  }

  template <typename Ttrack>
  BufferedTrack makeBufferedTrack(const Ttrack& track, bool isNucleus)
  {
    BufferedTrack bufferedTrack;
    bufferedTrack.momentum = {track.px(), track.py(), track.pz()};
    bufferedTrack.ptValue = track.pt();
    bufferedTrack.etaValue = track.eta();
    bufferedTrack.phiValue = track.phi();
    bufferedTrack.signValue = track.sign();
    bufferedTrack.dcaXY = track.dcaXY();
    bufferedTrack.dcaZ = track.dcaZ();
    bufferedTrack.tpcSignal = track.tpcSignal();
    bufferedTrack.tpcInnerParam = isNucleus && useHelium3Nucleus() ? correctedTPCInnerParamHe3(track) : track.tpcInnerParam();
    bufferedTrack.tpcNClsFound = track.tpcNClsFound();
    bufferedTrack.tpcNClsCrossedRows = track.tpcNClsCrossedRows();
    bufferedTrack.tpcNClsShared = track.tpcNClsShared();
    bufferedTrack.itsNCls = track.itsNCls();
    bufferedTrack.tpcChi2NCl = track.tpcChi2NCl();
    bufferedTrack.pidForTracking = track.pidForTracking();
    bufferedTrack.itsClusterSizes = track.itsClusterSizes();
    bufferedTrack.trackId = track.globalIndex();

    if (isNucleus) {
      bufferedTrack.nSigmaTPC = getNucleusTPCNSigma(track);
      bufferedTrack.nSigmaTOF = getNucleusTOFNSigma(track);
      bufferedTrack.nSigmaITS = getNucleusITSNSigma(track);
    } else {
      bufferedTrack.nSigmaTPC = getHadronTPCNSigma(track);
      bufferedTrack.nSigmaTOF = getHadronTOFNSigma(track);
      bufferedTrack.nSigmaITS = getHadronITSNSigma(track);
      bufferedTrack.nSigmaTPCHadPi = track.tpcNSigmaPi();
      bufferedTrack.nSigmaTPCHadKa = track.tpcNSigmaKa();
      bufferedTrack.nSigmaTPCHadPr = track.tpcNSigmaPr();
      bufferedTrack.nSigmaTOFHadPi = track.tofNSigmaPi();
      bufferedTrack.nSigmaTOFHadKa = track.tofNSigmaKa();
      bufferedTrack.nSigmaTOFHadPr = track.tofNSigmaPr();
    }

    if (track.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(track);
      beta = std::clamp(beta, 1.e-4f, 1.f - 1.e-6f);
      const float momentumForMass = isNucleus ? bufferedTrack.tpcInnerParam * nucleusChargeFactor() : bufferedTrack.tpcInnerParam;
      bufferedTrack.massTOF = momentumForMass * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    return bufferedTrack;
  }

  float configuredHadronMass() const
  {
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kPiPlus)) {
      return static_cast<float>(o2::constants::physics::MassPiPlus);
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kKPlus)) {
      return static_cast<float>(o2::constants::physics::MassKPlus);
    }
    if (species.settingHadPDGCode.value == static_cast<int>(PDG_t::kProton)) {
      return static_cast<float>(o2::constants::physics::MassProton);
    }
    return 0.f;
  }

  void fillBufferedCandidateInfo(const BufferedTrack& nucleus, const BufferedTrack& hadron, HadNucandidate& hadNucand)
  {
    const float nuChargeFactor = nucleusChargeFactor();
    hadNucand.momNu = nucleus.momentum;
    for (size_t iComponent = 0; iComponent < hadNucand.momNu.size(); ++iComponent) {
      hadNucand.momNu[iComponent] *= nuChargeFactor;
    }
    hadNucand.momHad = hadron.momentum;

    const float massHadron = configuredHadronMass();
    hadNucand.invMass = RecoDecay::m(std::array<std::array<float, 3>, 2>{hadNucand.momNu, hadNucand.momHad}, std::array<float, 2>{nucleusMass(), massHadron});
    hadNucand.signNu = nucleus.sign();
    hadNucand.signHad = hadron.sign();
    computeClosePairDeltas(nucleus, hadron, hadNucand.deltaEta, hadNucand.deltaPhi);
    hadNucand.dcaxyNu = nucleus.dcaXY;
    hadNucand.dcazNu = nucleus.dcaZ;
    hadNucand.dcaxyHad = hadron.dcaXY;
    hadNucand.dcazHad = hadron.dcaZ;
    hadNucand.tpcSignalNu = nucleus.tpcSignal;
    hadNucand.tpcSignalHad = hadron.tpcSignal;
    hadNucand.momNuTPC = nucleus.tpcInnerParam;
    hadNucand.momHadTPC = hadron.tpcInnerParam;
    hadNucand.nTPCClustersNu = nucleus.tpcNClsFound;
    hadNucand.nTPCClustersHad = hadron.tpcNClsFound;
    hadNucand.nTPCCrossedRowsNu = nucleus.tpcNClsCrossedRows;
    hadNucand.nTPCCrossedRowsHad = hadron.tpcNClsCrossedRows;
    hadNucand.sharedClustersNu = nucleus.tpcNClsShared;
    hadNucand.sharedClustersHad = hadron.tpcNClsShared;
    hadNucand.nClsItsNu = nucleus.itsNCls;
    hadNucand.nClsItsHad = hadron.itsNCls;
    hadNucand.chi2TPCNu = nucleus.tpcChi2NCl;
    hadNucand.chi2TPCHad = hadron.tpcChi2NCl;
    hadNucand.nSigmaNu = nucleus.nSigmaTPC;
    hadNucand.nSigmaHad = hadron.nSigmaTPC;
    hadNucand.nSigmaTOFNu = nucleus.nSigmaTOF;
    hadNucand.nSigmaITSNu = nucleus.nSigmaITS;
    hadNucand.nSigmaTOFHad = hadron.nSigmaTOF;
    hadNucand.nSigmaITSHad = hadron.nSigmaITS;
    hadNucand.nSigmaTPCHadPi = hadron.nSigmaTPCHadPi;
    hadNucand.nSigmaTPCHadKa = hadron.nSigmaTPCHadKa;
    hadNucand.nSigmaTPCHadPr = hadron.nSigmaTPCHadPr;
    hadNucand.nSigmaTOFHadPi = hadron.nSigmaTOFHadPi;
    hadNucand.nSigmaTOFHadKa = hadron.nSigmaTOFHadKa;
    hadNucand.nSigmaTOFHadPr = hadron.nSigmaTOFHadPr;
    hadNucand.massTOFNu = nucleus.massTOF;
    hadNucand.massTOFHad = hadron.massTOF;
    hadNucand.pidTrkNu = nucleus.pidForTracking;
    hadNucand.pidTrkHad = hadron.pidForTracking;
    hadNucand.itsClSizeNu = nucleus.itsClusterSizes;
    hadNucand.itsClSizeHad = hadron.itsClusterSizes;
    hadNucand.trackIDNu = static_cast<int>(nucleus.trackId);
    hadNucand.trackIDHad = static_cast<int>(hadron.trackId);
    hadNucand.isBkgUS = nucleus.sign() * hadron.sign() < 0;
    hadNucand.isBkgEM = true;

    float massLightNucleusForKstarMt = nucleusMass();
    if (useDeuteronNucleus() && deuteronPid.settingUseProtonMassForKstarMt) {
      massLightNucleusForKstarMt = static_cast<float>(o2::constants::physics::MassProton);
    }
    hadNucand.kstar = computePairKstar(hadNucand.momHad, massHadron, hadNucand.momNu, massLightNucleusForKstarMt);
    hadNucand.mT = computePairMT(hadNucand.momHad, massHadron, hadNucand.momNu, massLightNucleusForKstarMt);
  }

  void fillDataTable(const HadNucandidate& hadNucand)
  {
    mOutputDataTable(
      hadNucand.recoPtNu(),
      hadNucand.recoEtaNu(),
      hadNucand.recoPhiNu(),
      hadNucand.recoPtHad(),
      hadNucand.recoEtaHad(),
      hadNucand.recoPhiHad(),
      hadNucand.dcaxyNu,
      hadNucand.dcazNu,
      hadNucand.dcaxyHad,
      hadNucand.dcazHad,
      hadNucand.dcaPair,
      hadNucand.tpcSignalNu,
      hadNucand.momNuTPC,
      hadNucand.tpcSignalHad,
      hadNucand.momHadTPC,
      hadNucand.nTPCClustersNu,
      hadNucand.nTPCClustersHad,
      hadNucand.nTPCCrossedRowsNu,
      hadNucand.nTPCCrossedRowsHad,
      hadNucand.nSigmaNu,
      hadNucand.nSigmaTPCHadPi,
      hadNucand.nSigmaTPCHadKa,
      hadNucand.nSigmaTPCHadPr,
      hadNucand.nSigmaTOFHadPi,
      hadNucand.nSigmaTOFHadKa,
      hadNucand.nSigmaTOFHadPr,
      hadNucand.chi2TPCNu,
      hadNucand.chi2TPCHad,
      hadNucand.massTOFNu,
      hadNucand.massTOFHad,
      hadNucand.pidTrkNu,
      hadNucand.pidTrkHad,
      hadNucand.itsClSizeNu,
      hadNucand.itsClSizeHad,
      hadNucand.sharedClustersNu,
      hadNucand.sharedClustersHad,
      hadNucand.deltaEta,
      hadNucand.deltaPhi,
      hadNucand.nSigmaHad,
      hadNucand.nSigmaTOFNu,
      hadNucand.nSigmaITSNu,
      hadNucand.nSigmaTOFHad,
      hadNucand.nSigmaITSHad);
  }

  template <typename Tcoll>
  void fillTable(const HadNucandidate& hadNucand, const Tcoll& collision)
  {
    fillDataTable(hadNucand);
    if (output.settingFillMultiplicity) {
      mOutputMultiplicityTable(
        collision.globalIndex(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0C(),
        collision.multFT0C());
    }
  }

  void fillTable(const HadNucandidate& hadNucand, const BufferedCollision& collision)
  {
    fillDataTable(hadNucand);
    if (output.settingFillMultiplicity) {
      mOutputMultiplicityTable(
        collision.eventId,
        collision.posZ,
        collision.numContrib,
        collision.centFT0C,
        collision.multFT0C);
    }
  }

  template <typename TparticleNu, typename TparticleHad>
  void fillMCTable(const HadNucandidate& hadNucand, const TparticleNu& particleNu, const TparticleHad& particleHad, bool sameMCCollision, bool matchesRecoMCCollision)
  {
    const float signedPtNuMC = particleNu.pdgCode() >= 0 ? particleNu.pt() : -particleNu.pt();
    const float signedPtHadMC = particleHad.pdgCode() >= 0 ? particleHad.pt() : -particleHad.pt();
    const std::array<float, 3> momentumNuMC{particleNu.px(), particleNu.py(), particleNu.pz()};
    const std::array<float, 3> momentumHadMC{particleHad.px(), particleHad.py(), particleHad.pz()};
    const float kstarMC = computePairKstar(momentumHadMC, mMassHad, momentumNuMC, nucleusMass());

    mOutputMCTable(
      signedPtNuMC,
      particleNu.eta(),
      particleNu.phi(),
      signedPtHadMC,
      particleHad.eta(),
      particleHad.phi(),
      kstarMC,
      particleNu.pdgCode(),
      particleHad.pdgCode(),
      particleNu.isPhysicalPrimary(),
      particleHad.isPhysicalPrimary(),
      sameMCCollision,
      matchesRecoMCCollision);

    mQaRegistry.fill(HIST("MC/hKstarRecVsGen"), kstarMC, hadNucand.kstar);
    mQaRegistry.fill(HIST("MC/hPtNuRecVsGen"), signedPtNuMC, hadNucand.recoPtNu());
    mQaRegistry.fill(HIST("MC/hPtHadRecVsGen"), signedPtHadMC, hadNucand.recoPtHad());
  }

  void fillHistograms(const HadNucandidate& hadNucand)
  {
    mQaRegistry.fill(HIST("hNuPt"), hadNucand.recoPtNu());
    mQaRegistry.fill(HIST("hHadPt"), hadNucand.recoPtHad());
    mQaRegistry.fill(HIST("hNuPin"), hadNucand.momNuTPC * hadNucand.signNu);
    mQaRegistry.fill(HIST("hHadPin"), hadNucand.momHadTPC * hadNucand.signHad);
    mQaRegistry.fill(HIST("hNuEta"), hadNucand.recoEtaNu());
    mQaRegistry.fill(HIST("hHadEta"), hadNucand.recoEtaHad());
    mQaRegistry.fill(HIST("hNuPhi"), hadNucand.recoPhiNu());
    mQaRegistry.fill(HIST("hHadPhi"), hadNucand.recoPhiHad());
    mQaRegistry.fill(HIST("hNuHadtInvMass"), hadNucand.invMass);
    mQaRegistry.fill(HIST("hdcaxyNu"), hadNucand.dcaxyNu);
    mQaRegistry.fill(HIST("hdcazNu"), hadNucand.dcazNu);
    mQaRegistry.fill(HIST("hdcaxyHad"), hadNucand.dcaxyHad);
    mQaRegistry.fill(HIST("hdcazHad"), hadNucand.dcazHad);
    mQaRegistry.fill(HIST("hNClsNuITS"), hadNucand.nClsItsNu);
    mQaRegistry.fill(HIST("hNClsHadITS"), hadNucand.nClsItsHad);
    mQaRegistry.fill(HIST("hisBkgEM"), hadNucand.isBkgEM);
  }

  void fillKstarAtCentrality(const HadNucandidate& hadNucand, float centrality)
  {
    if (!hadNucand.isBkgUS) {
      if (hadNucand.recoPtNu() > 0) {
        mQaRegistry.fill(HIST("hkStar_LS_M"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_LS_M"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_LS_M"), hadNucand.kstar, centrality);
      } else {
        mQaRegistry.fill(HIST("hkStar_LS_A"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_LS_A"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_LS_A"), hadNucand.kstar, centrality);
      }
    } else {
      if (hadNucand.recoPtNu() > 0) {
        mQaRegistry.fill(HIST("hkStar_US_M"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_US_M"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_US_M"), hadNucand.kstar, centrality);
      } else {
        mQaRegistry.fill(HIST("hkStar_US_A"), hadNucand.kstar);
        mQaRegistry.fill(HIST("hkStaVsmT_US_A"), hadNucand.kstar, hadNucand.mT);
        mQaRegistry.fill(HIST("hkStaVsCent_US_A"), hadNucand.kstar, centrality);
      }
    }
  }

  template <typename Tcoll>
  void fillKstar(const HadNucandidate& hadNucand, const Tcoll& collision)
  {
    fillKstarAtCentrality(hadNucand, collision.centFT0C());
  }

  void fillMixedPair(const BufferedTrack& nucleus, const BufferedTrack& hadron, const BufferedCollision& nucleusCollision)
  {
    if (isClosePair(nucleus, hadron, /*fillQA*/ true)) {
      return;
    }

    HadNucandidate hadNucand;
    fillBufferedCandidateInfo(nucleus, hadron, hadNucand);
    fillKstarAtCentrality(hadNucand, nucleusCollision.centFT0C);
    fillHistograms(hadNucand);

    if (output.settingFillTable && shouldFillOutputTable(hadNucand)) {
      fillTable(hadNucand, nucleusCollision);
    }
  }

  // ==================================================================================================================

  template <typename Tcollisions, typename Ttracks>
  void fillPairs(const Tcollisions& collisions, const Ttracks& tracks, const bool isMixedEvent)
  {
    for (const auto& trackPair : mTrackPairs) {

      auto deTrack = tracks.rawIteratorAt(trackPair.tr0Idx);
      auto hadTrack = tracks.rawIteratorAt(trackPair.tr1Idx);
      auto collBracket = trackPair.collBracket;

      HadNucandidate hadNucand;
      if (!fillCandidateInfo(deTrack, hadTrack, collBracket, collisions, hadNucand, tracks, isMixedEvent)) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(hadNucand.collisionID);
      fillKstar(hadNucand, collision);
      fillHistograms(hadNucand);

      if (output.settingFillTable && shouldFillOutputTable(hadNucand)) {
        fillTable(hadNucand, collision);
      }
    }
  }

  template <typename Ttrack>
  HadHyperTrackInfo makeHadHyperTrackInfo(const Ttrack& track)
  {
    constexpr float InvalidPID = -999.f;
    return {track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign()),
            track.dcaXY(), track.tpcNClsCrossedRows(), track.tpcNClsPID(), track.tpcChi2NCl(),
            track.itsClusterSizes(), track.itsChi2NCl(), track.hasTOF(),
            track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
            track.hasTOF() ? track.tofNSigmaPi() : InvalidPID,
            track.hasTOF() ? track.tofNSigmaKa() : InvalidPID,
            track.hasTOF() ? track.tofNSigmaPr() : InvalidPID};
  }

  template <typename Tcollision>
  HadHyperEventInfo makeHadHyperEventInfo(const Tcollision& collision)
  {
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    return {bc.runNumber(), collision.posZ(), collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
            collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange(),
            static_cast<uint16_t>(collision.numContrib()), collision.multFT0C(),
            collision.posX(), collision.posY(), collision.posZ()};
  }

  template <typename Tparticle>
  HadHyperParticleTruth makeHadHyperParticleTruth(const Tparticle& particle)
  {
    return {particle.globalIndex(), particle.mcCollisionId(), particle.pdgCode(), particle.pt(), particle.eta(), particle.phi(), particle.isPhysicalPrimary(), static_cast<int16_t>(particle.statusCode()), static_cast<int16_t>(particle.getProcess()), {particle.px(), particle.py(), particle.pz()}};
  }

  template <bool isMC, typename Ttrack>
  HadHyperHadron makeHadHyperHadron(const Ttrack& track)
  {
    HadHyperHadron hadron{makeHadHyperTrackInfo(track), {}, {}, track.globalIndex(), {track.px(), track.py(), track.pz()}, track.tpcInnerParam(), track.tpcSignal(), track.eta(), track.phi(), static_cast<int8_t>(track.sign())};
    if constexpr (isMC) {
      if (track.has_mcParticle()) {
        const auto particle = track.template mcParticle_as<aod::McParticles>();
        hadron.truth = makeHadHyperParticleTruth(particle);
        if (particle.has_mothers()) {
          for (const auto& mother : particle.template mothers_as<aod::McParticles>()) {
            hadron.motherIds.push_back(mother.globalIndex());
          }
        }
      }
    }
    return hadron;
  }

  template <bool isMC, typename Tcandidate, typename Ttracks>
  HadHyperCandidate makeHadHyperCandidate(const Tcandidate& candidate, const Ttracks& tracks)
  {
    const auto he = tracks.rawIteratorAt(candidate.heTrackId() - tracks.offset());
    const auto pion = tracks.rawIteratorAt(candidate.piTrackId() - tracks.offset());
    const std::array<float, 3> heMomentum{
      candidate.ptHe3() * std::cos(candidate.phiHe3()),
      candidate.ptHe3() * std::sin(candidate.phiHe3()),
      candidate.ptHe3() * std::sinh(candidate.etaHe3())};
    const std::array<float, 3> pionMomentum{
      candidate.ptPi() * std::cos(candidate.phiPi()),
      candidate.ptPi() * std::sin(candidate.phiPi()),
      candidate.ptPi() * std::sinh(candidate.etaPi())};
    HadHyperCandidate result;
    for (size_t i = 0; i < result.momentum.size(); ++i) {
      result.momentum[i] = heMomentum[i] + pionMomentum[i];
    }
    result.mass = computeHyperCandidateMass(candidate);
    result.heTrackId = candidate.heTrackId();
    result.pionTrackId = candidate.piTrackId();
    result.sourceCandidateId = candidate.globalIndex();
    result.etaHe = candidate.etaHe3();
    result.phiHe = candidate.phiHe3();
    result.etaPi = candidate.etaPi();
    result.phiPi = candidate.phiPi();
    result.isMatter = candidate.isMatter();
    if constexpr (isMC) {
      if (he.has_mcParticle()) {
        result.heTruth = makeHadHyperParticleTruth(he.template mcParticle_as<aod::McParticles>());
      }
      if (pion.has_mcParticle()) {
        result.decayPionTruth = makeHadHyperParticleTruth(pion.template mcParticle_as<aod::McParticles>());
      }
      if (he.has_mcParticle() && pion.has_mcParticle()) {
        const auto heParticle = he.template mcParticle_as<aod::McParticles>();
        const auto pionParticle = pion.template mcParticle_as<aod::McParticles>();
        if (heParticle.has_mothers() && pionParticle.has_mothers()) {
          for (const auto& heMother : heParticle.template mothers_as<aod::McParticles>()) {
            for (const auto& pionMother : pionParticle.template mothers_as<aod::McParticles>()) {
              if (heMother.globalIndex() == pionMother.globalIndex() &&
                  std::abs(heMother.pdgCode()) == HyperTritonPDG) {
                result.hyperTruth = makeHadHyperParticleTruth(heMother);
              }
            }
          }
        }
      }
    }
    result.info = {candidate.isMatter(),
                   candidate.ptHe3(), candidate.etaHe3(), candidate.phiHe3(),
                   candidate.ptPi(), candidate.etaPi(), candidate.phiPi(),
                   candidate.dcaV0Daug(), candidate.dcaHe(), candidate.dcaPi(),
                   candidate.nSigmaHe(), candidate.tofMass(),
                   candidate.nTPCCrossedRowsHe(), candidate.nTPCCrossedRowsPi(),
                   candidate.tpcMomHe(), candidate.tpcMomPi(),
                   candidate.tpcSignalHe(), candidate.tpcSignalPi(),
                   candidate.tpcChi2He(), candidate.itsChi2He(), candidate.itsChi2Pi(),
                   candidate.itsClusterSizesHe(), candidate.itsClusterSizesPi(),
                   candidate.xDecVtx(), candidate.yDecVtx(), candidate.zDecVtx()};
    if constexpr (isMC) {
      result.statusCode = static_cast<int16_t>(candidate.statusCode());
      result.isReco = candidate.isReco();
      result.isSignal = candidate.isSignal();
      result.isRecoMCCollision = candidate.isRecoMCCollision();
      result.isSurvEvSel = candidate.isSurvEvSel();
      result.isTwoBodyDecay = candidate.isTwoBodyDecay();
      result.isFakeHeOnITSLayer = candidate.isFakeHeOnITSLayer();
      result.genPt = candidate.genPt();
      result.genEta = candidate.genEta();
      result.genPhi = candidate.genPhi();
      result.genPtHe3 = candidate.genPtHe3();
      result.genDecVtx = {candidate.genXDecVtx(), candidate.genYDecVtx(), candidate.genZDecVtx()};
      const float absGenPt = std::abs(candidate.genPt());
      result.genMomentum = {absGenPt * std::cos(candidate.genPhi()),
                            absGenPt * std::sin(candidate.genPhi()),
                            absGenPt * std::sinh(candidate.genEta())};
    }
    hadHyperRegistry.fill(HIST("hMass"), result.mass);
    hadHyperRegistry.fill(HIST("hDaughterHeTPC"), he.tpcInnerParam(), he.tpcSignal());
    hadHyperRegistry.fill(HIST("hDaughterPiTPC"), pion.tpcInnerParam(), pion.tpcSignal());
    return result;
  }

  HadHyperDataInfo makeHadHyperDataInfo(const HadHyperCandidate& candidate, const HadHyperHadron& hadron,
                                        const HadHyperEvent& hyperEvent, const HadHyperEvent& hadronEvent,
                                        bool mixed, int mixingDepth)
  {
    (void)hadronEvent;
    return std::tuple_cat(std::make_tuple(mixed, mixingDepth),
                          hyperEvent.info, candidate.info,
                          hadron.info);
  }

  bool hasTruthMother(const HadHyperHadron& hadron, int64_t motherId) const
  {
    return motherId >= 0 && std::find(hadron.motherIds.begin(), hadron.motherIds.end(), motherId) != hadron.motherIds.end();
  }

  bool isHadHyperTruthSelfCorrelation(const HadHyperCandidate& candidate, const HadHyperHadron& hadron) const
  {
    return hadron.truth.particleId >= 0 &&
           (hadron.truth.particleId == candidate.heTruth.particleId ||
            hadron.truth.particleId == candidate.decayPionTruth.particleId ||
            hasTruthMother(hadron, candidate.hyperTruth.particleId));
  }

  HadHyperMCInfo makeHadHyperMCInfo(const HadHyperCandidate& candidate, const HadHyperHadron& hadron,
                                    const HadHyperEvent& hyperEvent, const HadHyperEvent& hadronEvent) const
  {
    const bool sameMCCollision = candidate.hyperTruth.collisionId >= 0 &&
                                 hadron.truth.collisionId == candidate.hyperTruth.collisionId;
    const bool matchesHypRecoMCCollision = hyperEvent.hasMCCollision &&
                                           candidate.hyperTruth.collisionId == hyperEvent.mcCollisionId;
    const bool matchesPairRecoMCCollision = hadronEvent.hasMCCollision &&
                                            hadron.truth.collisionId == hadronEvent.mcCollisionId;
    const bool isTruthSelfCorrelation = isHadHyperTruthSelfCorrelation(candidate, hadron);
    const bool isTruePrimaryHadHyperPair = candidate.isSignal &&
                                           std::abs(candidate.hyperTruth.pdgCode) == HyperTritonPDG &&
                                           std::abs(hadron.truth.pdgCode) == std::abs(species.settingHadPDGCode.value) &&
                                           hadron.truth.isPhysicalPrimary && sameMCCollision && !isTruthSelfCorrelation;
    return {candidate.genPt, candidate.genEta, candidate.genPhi,
            candidate.genDecVtx[0], candidate.genDecVtx[1], candidate.genDecVtx[2],
            candidate.isReco, candidate.isSignal,
            candidate.isRecoMCCollision, candidate.isSurvEvSel, candidate.isTwoBodyDecay, candidate.statusCode,
            candidate.heTruth.pt, candidate.heTruth.isPhysicalPrimary,
            candidate.decayPionTruth.isPhysicalPrimary,
            std::get<0>(hadron.info), std::get<1>(hadron.info), std::get<2>(hadron.info),
            hadron.truth.pt, hadron.truth.eta, hadron.truth.phi, hadron.truth.isPhysicalPrimary,
            hadron.truth.process,
            sameMCCollision, matchesHypRecoMCCollision, matchesPairRecoMCCollision,
            isTruthSelfCorrelation, isTruePrimaryHadHyperPair};
  }

  void fillHadHyperMCQA(const HadHyperCandidate& candidate, const HadHyperHadron& hadron,
                        const HadHyperEvent& hyperEvent, float kstar)
  {
    if (!candidate.isReco || !candidate.isSignal ||
        std::abs(candidate.hyperTruth.pdgCode) != HyperTritonPDG) {
      return;
    }

    const bool isSelectedHadronSpecies = std::abs(hadron.truth.pdgCode) == std::abs(species.settingHadPDGCode.value);
    if (!isSelectedHadronSpecies) {
      return;
    }

    const float kstarMC = computePairKstar(hadron.truth.momentum, configuredHadronMass(),
                                           candidate.genMomentum, o2::constants::physics::MassHyperTriton);
    if (std::isfinite(kstarMC)) {
      hadHyperRegistry.fill(HIST("MC/hKstarRecVsGenHyperReco"), kstarMC, kstar);
      hadHyperRegistry.fill(HIST("MC/hKstarResolutionHyperReco"), kstar, kstar - kstarMC);
    }

    hadHyperRegistry.fill(HIST("MC/hPrimaryHadronVsKstarDen"), kstar);
    hadHyperRegistry.fill(HIST("MC/hPrimaryHadronVsCentDen"), hyperEvent.centrality);
    if (hadron.truth.isPhysicalPrimary) {
      hadHyperRegistry.fill(HIST("MC/hPrimaryHadronVsKstarNum"), kstar);
      hadHyperRegistry.fill(HIST("MC/hPrimaryHadronVsCentNum"), hyperEvent.centrality);
    }
  }

  template <bool isMC>
  bool fillHadHyperPair(const HadHyperCandidate& candidate, const HadHyperHadron& hadron,
                        const HadHyperEvent& hyperEvent, const HadHyperEvent& hadronEvent,
                        bool mixed, int mixingDepth)
  {
    const float kstar = computePairKstar(hadron.momentum, configuredHadronMass(),
                                         candidate.momentum, o2::constants::physics::MassHyperTriton);
    // All source indices here belong to the same timeframe. Also clean mixed
    // pairs in case LF reassigned a daughter to a different collision.
    if (hadron.sourceId == candidate.heTrackId || hadron.sourceId == candidate.pionTrackId) {
      const int reason = hadron.sourceId == candidate.heTrackId ? 0 : 1;
      hadHyperRegistry.fill(HIST("hSelfPairs"), reason, kstar);
      if (!mixed) {
        hadHyperRegistry.fill(HIST("hSameEventSelfPairs"), reason, kstar);
      }
      return false;
    }
    if constexpr (isMC) {
      if (isHadHyperTruthSelfCorrelation(candidate, hadron)) {
        hadHyperRegistry.fill(HIST("hSelfPairs"), 2, kstar);
        if (!mixed) {
          hadHyperRegistry.fill(HIST("hSameEventSelfPairs"), 2, kstar);
        }
        return false;
      }
    }
    const bool unlikeSign = candidate.sign() * hadron.sign() < 0;
    if (!eventMixing.settingSaveUSandLS && unlikeSign != eventMixing.settingEnableBkgUS.value) {
      return false;
    }
    if (!std::isfinite(kstar)) {
      return false;
    }
    if (isCloseHadHyperPairAtPV(candidate, hadron, /*fillQA*/ true)) {
      return false;
    }
    if (mixed) {
      hadHyperRegistry.fill(HIST("hME"), kstar);
    } else {
      hadHyperRegistry.fill(HIST("hSE"), kstar);
    }

    if constexpr (isMC) {
      fillHadHyperMCQA(candidate, hadron, hyperEvent, kstar);
    }

    if (!output.settingFillTable || (hadHyper.maxOutputKstar.value > 0.f && kstar >= hadHyper.maxOutputKstar.value)) {
      return true;
    }
    const auto dataInfo = makeHadHyperDataInfo(candidate, hadron, hyperEvent, hadronEvent, mixed, mixingDepth);
    if constexpr (isMC) {
      std::apply([this](const auto&... columns) { mOutputHadHyperMCTable(columns...); },
                 std::tuple_cat(dataInfo, makeHadHyperMCInfo(candidate, hadron, hyperEvent, hadronEvent)));
    } else {
      std::apply([this](const auto&... columns) { mOutputHadHyperDataTable(columns...); }, dataInfo);
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
  void collectHyperHadronTracks(HadHyperEvent& event, const Ttracks& eventTracks)
  {
    for (const auto& track : eventTracks) {
      if (!selectTrackHadron(track) || !selectionPIDHadron(track)) {
        continue;
      }
      event.hadrons.push_back(makeHadHyperHadron<isMC>(track));
      hadHyperRegistry.fill(HIST("hHadronTPC"), track.tpcInnerParam(), track.tpcSignal());
    }
  }

  template <bool isMC, typename Tcandidates, typename Ttracks>
  void collectHyperCandidates(HadHyperEvent& event, const Tcandidates& eventCandidates, const Ttracks& tracks)
  {
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
      event.candidates.push_back(makeHadHyperCandidate<isMC>(candidate, tracks));
    }
  }

  template <bool isMC, typename Tcollision, typename Ttracks, typename Tcandidates>
  HadHyperEvent buildHyperEvent(const Tcollision& collision, const Ttracks& tracks, const Tcandidates& candidates)
  {
    HadHyperEvent event;
    event.info = makeHadHyperEventInfo(collision);
    event.centrality = collision.centFT0C();
    if constexpr (isMC) {
      event.hasMCCollision = collision.has_mcCollision();
      event.mcCollisionId = collision.has_mcCollision() ? collision.mcCollisionId() : -1;
    }

    if constexpr (isMC) {
      const auto eventTracks = tracks.sliceBy(mPerColMC, collision.globalIndex());
      collectHyperHadronTracks<isMC>(event, eventTracks);

      const auto eventCandidates = candidates.sliceBy(hypPerColMC, collision.globalIndex());
      collectHyperCandidates<isMC>(event, eventCandidates, tracks);
    } else {
      const auto eventTracks = tracks.sliceBy(mPerCol, collision.globalIndex());
      collectHyperHadronTracks<isMC>(event, eventTracks);

      const auto eventCandidates = candidates.sliceBy(hypPerCol, collision.globalIndex());
      collectHyperCandidates<isMC>(event, eventCandidates, tracks);
    }

    return event;
  }

  template <bool isMC>
  void fillSameEventHyperPairs(const HadHyperEvent& event)
  {
    for (const auto& candidate : event.candidates) {
      int acceptedPairs = 0;
      for (const auto& hadron : event.hadrons) {
        if (fillHadHyperPair<isMC>(candidate, hadron, event, event, false, 0)) {
          ++acceptedPairs;
        }
      }
      hadHyperRegistry.fill(HIST("hCandidatePairMultiplicitySE"), acceptedPairs, event.centrality);
    }
  }

  template <bool isMC>
  void fillMixedEventHyperPairs(const HadHyperEvent& currentEvent, const std::deque<HadHyperEvent>& pool)
  {
    const int depth = static_cast<int>(pool.size());
    hadHyperRegistry.fill(HIST("hMixingDepth"), depth);
    std::vector<int> currentCandidatePairCounts(currentEvent.candidates.size(), 0);
    for (const auto& partner : pool) {
      const float currentPosZ = std::get<1>(currentEvent.info);
      const float partnerPosZ = std::get<1>(partner.info);
      hadHyperRegistry.fill(HIST("hMixEventDeltaPosZVsCent"), currentEvent.centrality, currentPosZ - partnerPosZ);
      hadHyperRegistry.fill(HIST("hMixEventDeltaCentFT0CVsCent"), currentEvent.centrality, currentEvent.centrality - partner.centrality);
      size_t candidateIndex = 0;
      for (const auto& candidate : currentEvent.candidates) {
        for (const auto& hadron : partner.hadrons) {
          if (fillHadHyperPair<isMC>(candidate, hadron, currentEvent, partner, true, depth)) {
            ++currentCandidatePairCounts[candidateIndex];
          }
        }
        ++candidateIndex;
      }
      for (const auto& candidate : partner.candidates) {
        int acceptedPairs = 0;
        for (const auto& hadron : currentEvent.hadrons) {
          if (fillHadHyperPair<isMC>(candidate, hadron, partner, currentEvent, true, depth)) {
            ++acceptedPairs;
          }
        }
        hadHyperRegistry.fill(HIST("hCandidatePairMultiplicityME"), acceptedPairs, partner.centrality);
      }
    }
    for (const auto& acceptedPairs : currentCandidatePairCounts) {
      hadHyperRegistry.fill(HIST("hCandidatePairMultiplicityME"), acceptedPairs, currentEvent.centrality);
    }
  }

  void storeHyperEventInPool(std::deque<HadHyperEvent>& pool, HadHyperEvent&& event)
  {
    const int requestedMixingDepth = eventMixing.settingNoMixedEvents.value;
    if (requestedMixingDepth <= 0) {
      return;
    }
    const auto mixingDepth = static_cast<size_t>(requestedMixingDepth);
    if (pool.size() >= mixingDepth) {
      pool.pop_front();
    }
    pool.push_back(std::move(event));
  }

  template <bool isMC, typename Tcollisions, typename Ttracks, typename Tcandidates>
  void processHyperPairs(const Tcollisions& collisions, const Ttracks& tracks, const Tcandidates& candidates,
                         const aod::BCsWithTimestamps& bcs,
                         std::unordered_map<int, std::deque<HadHyperEvent>>& mixingPools,
                         int& mixingRunNumber)
  {
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

      hadHyperRegistry.fill(HIST("hPoolFlow"), 0);
      int poolBin = -1;
      if (hadHyper.enableMixing.value) {
        poolBin = configuredBinningPolicy.getBin(std::make_tuple(collision.posZ(), collision.centFT0C()));
        if (poolBin < 0) {
          continue;
        }
        hadHyperRegistry.fill(HIST("hPoolFlow"), 1);
      }

      auto event = buildHyperEvent<isMC>(collision, tracks, candidates);
      fillSameEventHyperPairs<isMC>(event);

      if (!hadHyper.enableMixing.value) {
        continue;
      }

      const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (mixingRunNumber != bc.runNumber()) {
        mixingPools.clear();
        mixingRunNumber = bc.runNumber();
      }

      auto& pool = mixingPools[poolBin];
      fillMixedEventHyperPairs<isMC>(event, pool);
      storeHyperEventInPool(pool, std::move(event));
    }
  }

  // ==================================================================================================================

  void processMC(const CollisionsFullMC& collisions, const TrackCandidatesMC& tracks, const aod::McParticles&, const aod::BCsWithTimestamps& bcs)
  {
    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);

    for (const auto& collision : collisions) {
      mTrackPairs.clear();

      if (!selectCollision</*isMC*/ true>(collision, bcs)) {
        continue;
      }
      if (mc.settingRequireRecoMCCollisionMatch.value && !collision.has_mcCollision()) {
        continue;
      }

      const uint64_t collIdx = collision.globalIndex();
      mGoodCollisions[collIdx] = true;
      auto tracksThisCollision = tracks.sliceBy(mPerColMC, collIdx);
      tracksThisCollision.bindExternalIndices(&tracks);

      // This is deliberately the same pair builder as for data. It applies
      // nucleus/pion selections and CPR before truth matching.
      pairTracksSameEvent(tracksThisCollision, collision.centFT0C());

      for (const auto& trackPair : mTrackPairs) {
        mQaRegistry.fill(HIST("MC/hPairFlow"), 0);
        auto trackNu = tracks.rawIteratorAt(trackPair.tr0Idx);
        auto trackHad = tracks.rawIteratorAt(trackPair.tr1Idx);

        if (!trackNu.has_mcParticle() || !trackHad.has_mcParticle()) {
          continue;
        }
        mQaRegistry.fill(HIST("MC/hPairFlow"), 1);

        const auto particleNu = trackNu.template mcParticle_as<aod::McParticles>();
        const auto particleHad = trackHad.template mcParticle_as<aod::McParticles>();
        const bool truthSpeciesMatch = std::abs(particleNu.pdgCode()) == std::abs(species.settingNuPDGCode.value) &&
                                       std::abs(particleHad.pdgCode()) == std::abs(species.settingHadPDGCode.value);
        if (mc.settingRequireTruthSpecies.value && !truthSpeciesMatch) {
          continue;
        }
        mQaRegistry.fill(HIST("MC/hPairFlow"), 2);

        const bool sameMCCollision = particleNu.mcCollisionId() == particleHad.mcCollisionId();
        if (mc.settingRequireSameMCCollision.value && !sameMCCollision) {
          continue;
        }
        mQaRegistry.fill(HIST("MC/hPairFlow"), 3);

        const bool matchesRecoMCCollision = collision.has_mcCollision() &&
                                            particleNu.mcCollisionId() == collision.mcCollisionId() &&
                                            particleHad.mcCollisionId() == collision.mcCollisionId();
        if (mc.settingRequireRecoMCCollisionMatch.value && !matchesRecoMCCollision) {
          continue;
        }
        mQaRegistry.fill(HIST("MC/hPairFlow"), 4);

        if (mc.settingRequirePhysicalPrimaries.value &&
            (!particleNu.isPhysicalPrimary() || !particleHad.isPhysicalPrimary())) {
          continue;
        }
        mQaRegistry.fill(HIST("MC/hPairFlow"), 5);

        HadNucandidate hadNucand;
        if (!fillCandidateInfo(trackNu, trackHad, trackPair.collBracket, collisions, hadNucand, tracks, /*isMixedEvent*/ false)) {
          continue;
        }
        mQaRegistry.fill(HIST("MC/hPairFlow"), 6);
        auto selectedCollision = collisions.rawIteratorAt(hadNucand.collisionID);
        fillKstar(hadNucand, selectedCollision);
        fillHistograms(hadNucand);

        if (output.settingFillTable && shouldFillOutputTable(hadNucand)) {
          fillTable(hadNucand, selectedCollision);
          fillMCTable(hadNucand, particleNu, particleHad, sameMCCollision, matchesRecoMCCollision);
        }
      }
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processMC, "Process reconstructed MC same-event pairs", false);

  void processSameEvent(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::BCsWithTimestamps& bcs)
  {
    mGoodCollisions.clear();
    mGoodCollisions.resize(collisions.size(), false);

    for (const auto& collision : collisions) {

      mTrackPairs.clear();

      if (!selectCollision</*isMC*/ false>(collision, bcs)) {
        continue;
      }

      mGoodCollisions[collision.globalIndex()] = true;
      const uint64_t collIdx = collision.globalIndex();
      auto trackTableThisCollision = tracks.sliceBy(mPerCol, collIdx);
      trackTableThisCollision.bindExternalIndices(&tracks);

      pairTracksSameEvent(trackTableThisCollision, collision.centFT0C());

      if (mTrackPairs.empty()) {
        continue;
      }

      fillPairs(collisions, tracks, /*isMixedEvent*/ false);
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processSameEvent, "Process Same event", false);

  void processMixedEvent(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::BCsWithTimestamps&)
  {
    LOG(debug) << "Processing mixed event";
    const BinningType configuredBinningPolicy{{axisVertex, axisCentrality}, true};

    for (const auto& collision : collisions) {
      mQaRegistry.fill(HIST("hMixedEventSelections"), 0);
      if (!passesEventSelection</*isMC*/ false>(collision)) {
        continue;
      }
      mQaRegistry.fill(HIST("hMixedEventSelections"), 1);
      if (!passesZorroSelection(collision)) {
        continue;
      }
      mQaRegistry.fill(HIST("hMixedEventSelections"), 2);

      const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (mMixingRunNumber != bc.runNumber()) {
        mMixingPools.clear();
        mMixingRunNumber = bc.runNumber();
      }

      BufferedCollision currentCollision;
      currentCollision.eventId = mNextMixedEventId++;
      currentCollision.posZ = collision.posZ();
      currentCollision.numContrib = collision.numContrib();
      currentCollision.centFT0C = collision.centFT0C();
      currentCollision.multFT0C = collision.multFT0C();

      const uint64_t collIdx = collision.globalIndex();
      auto tracksThisCollision = tracks.sliceBy(mPerCol, collIdx);
      tracksThisCollision.bindExternalIndices(&tracks);
      for (const auto& track : tracksThisCollision) {
        if (selectTrackNu(track) && selectionPIDNu(track)) {
          currentCollision.nuclei.push_back(makeBufferedTrack(track, /*isNucleus*/ true));
        }
        if (selectTrackHadron(track) && selectionPIDHadron(track)) {
          currentCollision.hadrons.push_back(makeBufferedTrack(track, /*isNucleus*/ false));
        }
      }

      mQaRegistry.fill(HIST("hMixedNucleiPerEvent"), currentCollision.nuclei.size());
      mQaRegistry.fill(HIST("hMixedHadronsPerEvent"), currentCollision.hadrons.size());

      const int poolBin = configuredBinningPolicy.getBin(std::make_tuple(collision.posZ(), collision.centFT0C()));
      auto& pool = mMixingPools[poolBin];
      mQaRegistry.fill(HIST("hMixingPoolOccupancy"), pool.size());

      for (const auto& bufferedCollision : pool) {
        mQaRegistry.fill(HIST("hMixedEventSelections"), 3);
        mQaRegistry.fill(HIST("hNcontributor"), collision.numContrib());
        mQaRegistry.fill(HIST("hVtxZ"), collision.posZ());

        for (const auto& nucleus : currentCollision.nuclei) {
          for (const auto& hadron : bufferedCollision.hadrons) {
            fillMixedPair(nucleus, hadron, currentCollision);
          }
        }
        for (const auto& nucleus : bufferedCollision.nuclei) {
          for (const auto& hadron : currentCollision.hadrons) {
            fillMixedPair(nucleus, hadron, bufferedCollision);
          }
        }
      }

      const int requestedMixingDepth = eventMixing.settingNoMixedEvents.value;
      if (requestedMixingDepth <= 0) {
        continue;
      }
      const auto mixingDepth = static_cast<size_t>(requestedMixingDepth);
      if (pool.size() >= mixingDepth) {
        pool.pop_front();
      }
      pool.push_back(std::move(currentCollision));
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processMixedEvent, "Process Mixed event", false);

  void processPurity(const CollisionsFull& collisions, const TrackCandidates& tracks, const aod::BCsWithTimestamps& bcs)
  {
    for (const auto& collision : collisions) {
      if (!selectCollision</*isMC*/ false>(collision, bcs)) {
        continue;
      }

      const uint64_t collIdx = collision.globalIndex();
      auto trackTableThisCollision = tracks.sliceBy(mPerCol, collIdx);
      trackTableThisCollision.bindExternalIndices(&tracks);

      for (const auto& track : trackTableThisCollision) {
        const bool passTrackHad = selectTrackHadron(track);
        const bool passTrackNu = selectTrackNu(track);

        mQaRegistry.fill(HIST("hTrackSel"), Selections::kNoCuts);
        if (passTrackHad) {
          mQaRegistry.fill(HIST("hTrackSel"), Selections::kTrackCuts);
        }

        fillNucleusTrackSelection(Selections::kNoCuts);
        if (passTrackNu) {
          fillNucleusTrackSelection(Selections::kTrackCuts);
        }

        if (passTrackHad && species.settingHadPDGCode == PDG_t::kPiPlus) {
          const float tpcNSigmaHad = track.tpcNSigmaPi();
          mQaRegistry.fill(HIST("purity/h2NsigmaHadTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);
          if (track.hasTOF() && std::abs(track.p()) > hadronPid.settingPionMomCombMin) {
            const float tofNSigmaHad = track.tofNSigmaPi();
            const float combNsigmaHad = std::sqrt(tofNSigmaHad * tofNSigmaHad + tpcNSigmaHad * tpcNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadComb_preselection"), track.sign() * track.pt(), combNsigmaHad);
          }
        } else if (passTrackHad && species.settingHadPDGCode == PDG_t::kKPlus) {
          const float tpcNSigmaHad = track.tpcNSigmaKa();
          mQaRegistry.fill(HIST("purity/h2NsigmaHadTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);
          if (track.hasTOF() && track.tpcInnerParam() >= hadronPid.settingCutPinMinTOFHad) {
            const float tofNSigmaHad = track.tofNSigmaKa();
            const float combNsigmaHad = std::sqrt(tofNSigmaHad * tofNSigmaHad + tpcNSigmaHad * tpcNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadComb_preselection"), track.sign() * track.pt(), combNsigmaHad);
          }
        } else if (passTrackHad && species.settingHadPDGCode == PDG_t::kProton) {
          constexpr float protonPCombMin = 0.75f;
          const float tpcNSigmaHad = track.tpcNSigmaPr();
          mQaRegistry.fill(HIST("purity/h2NsigmaHadTPC_preselection"), track.sign() * track.pt(), tpcNSigmaHad);
          if (track.hasTOF() && std::abs(track.tpcInnerParam()) >= protonPCombMin) {
            const float tofNSigmaHad = track.tofNSigmaPr();
            const float combNsigmaHad = std::sqrt(tofNSigmaHad * tofNSigmaHad + tpcNSigmaHad * tpcNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadTOF_preselection"), track.sign() * track.pt(), tofNSigmaHad);
            mQaRegistry.fill(HIST("purity/h2NsigmaHadComb_preselection"), track.sign() * track.pt(), combNsigmaHad);
          }
        }

        if (passTrackNu && useDeuteronNucleus()) {
          const float tpcNSigmaDe = output.settingUseBBcomputeDeNsigma ? computeNSigmaDe(track) : track.tpcNSigmaDe();
          const float absTPCInnerParam = std::abs(track.tpcInnerParam());
          if (absTPCInnerParam > deuteronPid.settingCutPinMinTOFITSDe) {
            mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselection"), track.sign() * track.pt(), tpcNSigmaDe);
            mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselecComp"), track.sign() * track.pt(), track.tpcNSigmaDe());
            if (track.hasTOF()) {
              const float tofNSigmaDe = track.tofNSigmaDe();
              const float combNsigmaDe = std::sqrt(tofNSigmaDe * tofNSigmaDe + tpcNSigmaDe * tpcNSigmaDe);
              mQaRegistry.fill(HIST("purity/h2NsigmaNuTOF_preselection"), track.sign() * track.pt(), tofNSigmaDe);
              mQaRegistry.fill(HIST("purity/h2NsigmaNuComb_preselection"), track.sign() * track.pt(), combNsigmaDe);
            }
          } else {
            o2::aod::ITSResponse itsResponse;
            const float itsNSigmaDe = itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track.itsClusterSizes(), track.p(), track.eta());
            mQaRegistry.fill(HIST("purity/h2NSigmaNuITS_preselection"), track.sign() * track.pt(), itsNSigmaDe);
            if (std::abs(itsNSigmaDe) <= deuteronPid.settingCutNsigmaITSDe) {
              mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselection"), track.sign() * track.pt(), tpcNSigmaDe);
              mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselecComp"), track.sign() * track.pt(), track.tpcNSigmaDe());
            }
          }
        } else if (passTrackNu && useHelium3Nucleus()) {
          const float tpcNSigmaHe3 = computeNSigmaHe3(track);
          const float signedPtHe3 = track.sign() * 2.f * track.pt();
          mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselection"), signedPtHe3, tpcNSigmaHe3);
          o2::aod::ITSResponse itsResponse;
          const float itsNSigmaHe3 = itsResponse.nSigmaITS<o2::track::PID::Helium3>(track.itsClusterSizes(), track.p(), track.eta());
          mQaRegistry.fill(HIST("purity/h2NSigmaNuITS_preselection"), signedPtHe3, itsNSigmaHe3);
        } else if (passTrackNu && useTritonNucleus()) {
          const float tpcNSigmaTr = track.tpcNSigmaTr();
          o2::aod::ITSResponse itsResponse;
          const float itsNSigmaTr = itsResponse.nSigmaITS<o2::track::PID::Triton>(track.itsClusterSizes(), std::abs(track.p()), track.eta());
          mQaRegistry.fill(HIST("purity/h2NsigmaNuTPC_preselection"), track.sign() * track.pt(), tpcNSigmaTr);
          mQaRegistry.fill(HIST("purity/h2NSigmaNuITS_preselection"), track.sign() * track.pt(), itsNSigmaTr);
        }

        const bool isHadronSelected = passTrackHad && selectionPIDHadron(track);
        const bool isNucleusSelected = passTrackNu && selectionPIDNu(track);
        if (!isHadronSelected && !isNucleusSelected) {
          continue;
        }

        if (isHadronSelected) {
          mQaRegistry.fill(HIST("hTrackSel"), Selections::kPID);
          mQaRegistry.fill(HIST("hSingleHadPt"), track.pt() * track.sign());
        }

        if (isNucleusSelected) {
          fillNucleusTrackSelection(Selections::kPID);
          mQaRegistry.fill(HIST("hSingleNuPt"), track.pt() * track.sign() * nucleusChargeFactor());
          mQaRegistry.fill(HIST("hSingleNuPin"), (useHelium3Nucleus() ? correctedTPCInnerParamHe3(track) : track.tpcInnerParam()) * track.sign());
        }
      }
    }
  }
  PROCESS_SWITCH(HadNucleiFemto, processPurity, "Process for hadron and nucleus purity QA", false);

  void processHyper(const HadHyperCollisionsFull& collisions, const TrackCandidates& tracks,
                    const HyperCandidates& candidates, const aod::BCsWithTimestamps& bcs)
  {
    processHyperPairs</*isMC*/ false>(collisions, tracks, candidates, bcs, mHyperMixingPools, mHyperMixingRunNumber);
  }
  PROCESS_SWITCH(HadNucleiFemto, processHyper, "Process same-event and mixed-event hadron-hypertriton pairs", false);

  void processMCHyper(const HadHyperCollisionsFullMC& collisions, const TrackCandidatesMC& tracks,
                      const HyperCandidatesMC& candidates, const aod::McParticles&, const aod::BCsWithTimestamps& bcs)
  {
    processHyperPairs</*isMC*/ true>(collisions, tracks, candidates, bcs, mHyperMCMixingPools, mHyperMCMixingRunNumber);
  }
  PROCESS_SWITCH(HadNucleiFemto, processMCHyper, "Process MC same-event and mixed-event hadron-hypertriton pairs", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HadNucleiFemto>(cfgc)};
}
