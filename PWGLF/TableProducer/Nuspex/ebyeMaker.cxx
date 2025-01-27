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

#include <vector>
#include <map>
#include <utility>
#include <random>
#include <string>
#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "CCDB/CcdbApi.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFEbyeTables.h"

#include "TDatabasePDG.h"
#include "TFormula.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime>;
using TracksFullPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::pidTOFPr>;
using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime>;
using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

namespace
{
constexpr int kNpart = 2;
constexpr float trackSels[12]{/* 60, */ 80, 100, 2, 3, /* 4,  */ 0.05, 0.1, /* 0.15,  */ 0.5, 1, /* 1.5, */ 2, 3 /* , 4 */, 2, 3, /*, 4 */};
constexpr float dcaSels[3]{10., 10., 10.};
constexpr double betheBlochDefault[kNpart][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}, {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
constexpr double betheBlochDefaultITS[6]{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32};
constexpr double estimatorsCorrelationCoef[2]{-0.669108, 1.04489};
constexpr double estimatorsSigmaPars[4]{0.933321, 0.0416976, -0.000936344, 8.92179e-06};
constexpr double deltaEstimatorNsigma[2]{5.5, 5.};
constexpr double partMass[kNpart]{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron};
constexpr double partPdg[kNpart]{2212, o2::constants::physics::kDeuteron};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNamesPar{"p", "d"};
static const std::vector<std::string> trackSelsNames{"tpcClsMid", "tpcClsTight", "chi2TpcTight", "chi2TpcMid", "dcaxyTight", "dcaxyMid", "dcazTight", "dcazMid", "tpcNsigmaTight", "tpcNsigmaMid", "itsNsigmaTight", "itsNsigmaMid"};
static const std::vector<std::string> dcaSelsNames{"dcaxy", "dcaz", "dca"};
static const std::vector<std::string> particleName{"p"};
std::array<std::shared_ptr<TH3>, kNpart> tofMass;
void momTotXYZ(std::array<float, 3>& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC)
{
  for (int i = 0; i < 3; ++i) {
    momA[i] = momB[i] + momC[i];
  }
}
float invMass2Body(std::array<float, 3> const& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC, float const& massB, float const& massC)
{
  float p2B = momB[0] * momB[0] + momB[1] * momB[1] + momB[2] * momB[2];
  float p2C = momC[0] * momC[0] + momC[1] * momC[1] + momC[2] * momC[2];
  float eB = std::sqrt(p2B + massB * massB);
  float eC = std::sqrt(p2C + massC * massC);
  float eA = eB + eC;
  float massA = std::sqrt(eA * eA - momA[0] * momA[0] - momA[1] * momA[1] - momA[2] * momA[2]);
  return massA;
}
float alphaAP(std::array<float, 3> const& momA, std::array<float, 3> const& momB, std::array<float, 3> const& momC)
{
  float momTot = std::sqrt(std::pow(momA[0], 2.) + std::pow(momA[1], 2.) + std::pow(momA[2], 2.));
  float lQlPos = (momB[0] * momA[0] + momB[1] * momA[1] + momB[2] * momA[2]) / momTot;
  float lQlNeg = (momC[0] * momA[0] + momC[1] * momA[1] + momC[2] * momA[2]) / momTot;
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}
float etaFromMom(std::array<float, 3> const& momA, std::array<float, 3> const& momB)
{
  if (std::sqrt((1.f * momA[0] + 1.f * momB[0]) * (1.f * momA[0] + 1.f * momB[0]) +
                (1.f * momA[1] + 1.f * momB[1]) * (1.f * momA[1] + 1.f * momB[1]) +
                (1.f * momA[2] + 1.f * momB[2]) * (1.f * momA[2] + 1.f * momB[2])) -
        (1.f * momA[2] + 1.f * momB[2]) <
      static_cast<float>(1e-7)) {
    if ((1.f * momA[2] + 1.f * momB[2]) < 0.f)
      return -100.f;
    return 100.f;
  }
  return 0.5f * std::log((std::sqrt((1.f * momA[0] + 1.f * momB[0]) * (1.f * momA[0] + 1.f * momB[0]) +
                                    (1.f * momA[1] + 1.f * momB[1]) * (1.f * momA[1] + 1.f * momB[1]) +
                                    (1.f * momA[2] + 1.f * momB[2]) * (1.f * momA[2] + 1.f * momB[2])) +
                          (1.f * momA[2] + 1.f * momB[2])) /
                         (std::sqrt((1.f * momA[0] + 1.f * momB[0]) * (1.f * momA[0] + 1.f * momB[0]) +
                                    (1.f * momA[1] + 1.f * momB[1]) * (1.f * momA[1] + 1.f * momB[1]) +
                                    (1.f * momA[2] + 1.f * momB[2]) * (1.f * momA[2] + 1.f * momB[2])) -
                          (1.f * momA[2] + 1.f * momB[2])));
}
float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
{
  return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
}
} // namespace

struct CandidateV0 {
  float pt = -999.f;
  float eta = -999.f;
  float mass = -999.f;
  float cpa = -999.f;
  float dcav0daugh = -999.f;
  float dcanegpv = -999.f;
  float dcapospv = -999.f;
  float dcav0pv = -999.f;
  float tpcnsigmaneg = -999.f;
  float tpcnsigmapos = -999.f;
  float genpt = -999.f;
  float geneta = -999.f;
  int pdgcode = -999;
  bool isreco = 0;
  int64_t mcIndex = -999;
  int64_t globalIndexPos = -999;
  int64_t globalIndexNeg = -999;
};

struct CandidateTrack {
  float pt = -999.f;
  float eta = -999.f;
  uint8_t mass = 100;
  float dcapv = 0;
  float dcaxypv = 0;
  float dcazpv = 0;
  uint8_t tpcncls = 0;
  float tpcchi2 = 0;
  float tpcnsigma = -999.f;
  float itsnsigma = -999.f;
  float tofmass = -999.f;
  float outerPID = -999.f;
  float genpt = -999.f;
  float geneta = -999.f;
  int pdgcode = -999;
  int pdgcodemoth = -999;
  bool isreco = 0;
  int64_t mcIndex = -999;
  int64_t globalIndex = -999;
};

enum selBits {
  kTPCclsTight = BIT(0),
  kTPCclsMid = BIT(1),
  kChi2TPCTight = BIT(2),
  kChi2TPCMid = BIT(3),
  kDCAxyTight = BIT(4),
  kDCAxyMid = BIT(5),
  kDCAzTight = BIT(6),
  kDCAzMid = BIT(7),
  kITSPIDTight = BIT(8),
  kITSPIDMid = BIT(9),
  kTPCPIDTight = BIT(10),
  kTPCPIDMid = BIT(11)
};

enum PartTypes {
  kLa = BIT(20),
  kSig = BIT(21),
  kPhysPrim = BIT(22)
};

struct tagRun2V0MCalibration {
  bool mCalibrationStored = false;
  TH1* mhVtxAmpCorrV0A = nullptr;
  TH1* mhVtxAmpCorrV0C = nullptr;
  TH1* mhMultSelCalib = nullptr;
  float mMCScalePars[6] = {0.0};
  TFormula* mMCScale = nullptr;
} Run2V0MInfo;

struct tagRun2CL0Calibration {
  bool mCalibrationStored = false;
  TH1* mhVtxAmpCorr = nullptr;
  TH1* mhMultSelCalib = nullptr;
} Run2CL0Info;

struct ebyeMaker {
  Produces<aod::CollEbyeTable> collisionEbyeTable;
  Produces<aod::MiniCollTable> miniCollTable;
  Produces<aod::NucleiEbyeTable> nucleiEbyeTable;
  Produces<aod::LambdaEbyeTable> lambdaEbyeTable;
  Produces<aod::MiniTrkTable> miniTrkTable;
  Produces<aod::McNucleiEbyeTable> mcNucleiEbyeTable;
  Produces<aod::McLambdaEbyeTable> mcLambdaEbyeTable;
  Produces<aod::McMiniTrkTable> mcMiniTrkTable;
  std::mt19937 gen32;
  std::vector<CandidateV0> candidateV0s;
  std::array<std::vector<CandidateTrack>, 2> candidateTracks;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;
  std::vector<int> classIds;

  int mRunNumber;
  float d_bz;
  uint8_t nTrackletsColl;
  // o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 2, 6, particleNamesPar, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for deuteron"};
  Configurable<LabeledArray<double>> cfgBetheBlochParamsITS{"cfgBetheBlochParamsITS", {betheBlochDefaultITS, 1, 6, particleName, betheBlochParNames}, "ITS Bethe-Bloch parameterisation for deuteron"};

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis multAxis{"multAxis", {100, 0, 10000}, "Binning for the multiplicity axis"};
  ConfigurableAxis multFt0Axis{"multFt0Axis", {100, 0, 100000}, "Binning for the ft0 multiplicity axis"};

  // binning of (anti)lambda mass QA histograms
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.03f, o2::constants::physics::MassLambda0 + 0.03f}, "binning for the lambda invariant-mass"};

  // binning of PID QA histograms
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};
  ConfigurableAxis tpcAxis{"tpcAxis", {4.e2, 0.f, 4.e3f}, "tpc signal axis binning"};
  ConfigurableAxis tofMassAxis{"tofMassAxis", {1000, 0., 3.f}, "tof mass axis"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
  Configurable<float> etaMaxV0dau{"etaMaxV0dau", 0.8f, "maximum eta V0 daughters"};
  Configurable<float> outerPIDMin{"outerPIDMin", -4.f, "minimum outer PID"};

  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "fill histograms only for true signal candidates (MC)"};
  Configurable<std::string> genName{"genname", "", "Genearator name: HIJING, PYTHIA8, ... Default: \"\""};

  Configurable<uint8_t> triggerCut{"triggerCut", 0x0, "trigger cut to select"};
  Configurable<bool> kINT7Intervals{"kINT7Intervals", false, "toggle kINT7 trigger selection in the 10-30% and 50-90% centrality intervals (2018 Pb-Pb)"};
  Configurable<bool> kUseTPCPileUpCut{"kUseTPCPileUpCut", false, "toggle strong correlation cuts (Run 2)"};
  Configurable<bool> kUseEstimatorsCorrelationCut{"kUseEstimatorsCorrelationCut", false, "toggle cut on the correlation between centrality estimators (2018 Pb-Pb)"};

  Configurable<float> antidPtMin{"antidPtMin", 0.6f, "minimum antideuteron pT (GeV/c)"};
  Configurable<float> antidPtTof{"antidPtTof", 1.0f, "antideuteron pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antidPtMax{"antidPtMax", 1.8f, "maximum antideuteron pT (GeV/c)"};

  Configurable<float> antipPtMin{"antipPtMin", 0.4f, "minimum antiproton pT (GeV/c)"};
  Configurable<float> antipPtTof{"antipPtTof", 0.6f, "antiproton pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antipPtMax{"antipPtMax", 0.9f, "maximum antiproton pT (GeV/c)"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 1.f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 4.f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> trackNcrossedRows{"trackNcrossedRows", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> trackNclusItsCut{"trackNclusITScut", 2, "Minimum number of ITS clusters"};
  Configurable<float> trackNclusTpcCut{"trackNclusTPCcut", 60, "Minimum number of TPC clusters"};
  Configurable<float> trackChi2Cut{"trackChi2Cut", 4.f, "Maximum chi2/ncls in TPC"};
  Configurable<LabeledArray<float>> cfgDcaSels{"cfgDcaSels", {dcaSels, 1, 3, particleName, dcaSelsNames}, "DCA selections"};

  Configurable<float> v0trackNcrossedRows{"v0trackNcrossedRows", 100, "Minimum number of crossed TPC rows for V0 daughter"};
  Configurable<float> v0trackNclusItsCut{"v0trackNclusITScut", 0, "Minimum number of ITS clusters for V0 daughter"};
  Configurable<float> v0trackNclusTpcCut{"v0trackNclusTPCcut", 100, "Minimum number of TPC clusters for V0 daughter"};
  Configurable<float> v0trackNsharedClusTpc{"v0trackNsharedClusTpc", 5, "Maximum number of shared TPC clusters for V0 daughter"};
  Configurable<bool> v0requireITSrefit{"v0requireITSrefit", false, "require ITS refit for V0 daughter"};
  Configurable<float> vetoMassK0Short{"vetoMassK0Short", 0.01f, "veto for V0 compatible with K0s mass"};
  Configurable<float> v0radiusMax{"v0radiusMax", 100.f, "maximum V0 radius eccepted"};

  Configurable<float> antidNsigmaTpcCutLow{"antidNsigmaTpcCutLow", -4.f, "TPC PID cut low"};
  Configurable<float> antidNsigmaTpcCutUp{"antidNsigmaTpcCutUp", 4.f, "TPC PID cut up"};
  Configurable<float> antidTpcInnerParamMax{"tpcInnerParamMax", 0.f, "(temporary) tpc inner param cut"};
  Configurable<float> antidTofMassMax{"tofMassMax", 0.3f, "(temporary) tof mass cut"};

  Configurable<float> antipNsigmaTpcCutLow{"antipNsigmaTpcCutLow", -4.f, "TPC PID cut low"};
  Configurable<float> antipNsigmaTpcCutUp{"antipNsigmaTpcCutUp", 4.f, "TPC PID cut up"};
  Configurable<float> antipTpcInnerParamMax{"antipTpcInnerParamMax", 0.f, "(temporary) tpc inner param cut"};
  Configurable<float> antipTofMassMax{"antipTofMassMax", 0.3f, "(temporary) tof mass cut"};
  Configurable<float> tofMassMaxQA{"tofMassMaxQA", 0.6f, "(temporary) tof mass cut (for QA histograms)"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 0.5f, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcav0pv{"v0setting_dcav0pv", 1.f, "DCA V0 to Pv"};
  Configurable<float> v0setting_dcadaughtopv{"v0setting_dcadaughtopv", 0.1f, "DCA Pos To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.99f, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 5.f, "v0radius"};
  Configurable<float> v0setting_lifetime{"v0setting_lifetime", 40.f, "v0 lifetime cut"};
  Configurable<float> v0setting_nsigmatpc{"v0setting_nsigmatpc", 4.f, "nsigmatpc"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.02f, "maximum deviation from PDG mass (for QA histograms)"};

  Configurable<bool> constDCASel{"constDCASel", true, "use DCA selections independent of pt"};

  Configurable<float> antidItsClsSizeCut{"antidItsClsSizeCut", 1.e-10f, "cluster size cut for antideuterons"};
  Configurable<float> antidPtItsClsSizeCut{"antidPtItsClsSizeCut", 10.f, "pt for cluster size cut for antideuterons"};

  Configurable<LabeledArray<float>> cfgTrackSels{"cfgTrackSels", {trackSels, 1, 12, particleName, trackSelsNames}, "Track selections"};

  std::array<float, kNpart> ptMin;
  std::array<float, kNpart> ptTof;
  std::array<float, kNpart> ptMax;
  std::array<float, kNpart> nSigmaTpcCutLow;
  std::array<float, kNpart> nSigmaTpcCutUp;
  std::array<float, kNpart> tpcInnerParamMax;
  std::array<float, kNpart> tofMassMax;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Preslice<TracksFull> perCollisionTracksFull = o2::aod::track::collisionId;
  Preslice<TracksFullPID> perCollisionTracksFullPID = o2::aod::track::collisionId;
  Preslice<TracksFullIU> perCollisionTracksFullIU = o2::aod::track::collisionId;
  Preslice<aod::V0s> perCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::McParticles> perCollisionMcParts = o2::aod::mcparticle::mcCollisionId;

  template <class P>
  int getPartTypeMother(P const& mcPart)
  {
    for (auto& mother : mcPart.template mothers_as<aod::McParticles>()) {
      if (!mother.isPhysicalPrimary())
        return -1;
      int pdgCode = mother.pdgCode();
      switch (std::abs(pdgCode)) {
        case 3122: {
          int foundPi = 0;
          for (auto& mcDaught : mother.template daughters_as<aod::McParticles>()) {
            if (std::abs(mcDaught.pdgCode()) == 211) {
              foundPi = mcDaught.pdgCode();
              break;
            }
          }
          if (foundPi * mcPart.pdgCode() < -0.5)
            return PartTypes::kLa;
          return -1;
        }
        // case 3222:
        //   return PartTypes::kSig;
        // case 3112:
        //   return PartTypes::kSig;
        default:
          return -1;
      }
    }
    return -1;
  }

  template <class T>
  bool selectV0Daughter(T const& track)
  {
    if (std::abs(track.eta()) > etaMaxV0dau) {
      return false;
    }
    if (track.itsNCls() < v0trackNclusItsCut ||
        track.tpcNClsFound() < v0trackNclusTpcCut ||
        track.tpcNClsCrossedRows() < v0trackNclusTpcCut ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > v0trackNsharedClusTpc) {
      return false;
    }
    if (doprocessRun2 || doprocessMiniRun2 || doprocessMcRun2 || doprocessMiniMcRun2) {
      if (!(track.trackType() & o2::aod::track::Run2Track) ||
          !(track.flags() & o2::aod::track::TPCrefit)) {
        return false;
      }
      if (v0requireITSrefit && !(track.flags() & o2::aod::track::ITSrefit)) {
        return false;
      }
    }
    return true;
  }

  template <class T>
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (!(track.itsClusterMap() & 0x01) && !(track.itsClusterMap() & 0x02)) {
      return false;
    }
    if (track.itsNCls() < trackNclusItsCut ||
        track.tpcNClsFound() < trackNclusTpcCut ||
        track.tpcNClsCrossedRows() < trackNcrossedRows ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcChi2NCl() > trackChi2Cut ||
        track.itsChi2NCl() > 36.f) {
      return false;
    }
    if (doprocessRun2 || doprocessMiniRun2 || doprocessMcRun2 || doprocessMiniMcRun2) {
      if (!(track.trackType() & o2::aod::track::Run2Track) ||
          !(track.flags() & o2::aod::track::TPCrefit) ||
          !(track.flags() & o2::aod::track::ITSrefit)) {
        return false;
      }
    }
    return true;
  }

  template <class T>
  float getITSClSize(T const& track)
  {
    float sum{0.f};
    for (int iL{0}; iL < 6; ++iL) {
      sum += (track.itsClusterSizes() >> (iL * 4)) & 0xf;
    }
    return sum / track.itsNCls();
  }

  float dcaSigma(float const& pt)
  {
    return 0.0105 + 0.0350 / std::pow(std::abs(pt), 1.1);
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    classIds.clear();
    auto timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (doprocessRun2 || doprocessMcRun2 || doprocessMiniRun2 || doprocessMiniMcRun2) {
      auto grpPath{"GLO/GRP/GRP"};
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (!grpo) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpPath << " of object GRPObject for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      TList* callst = ccdb->getForTimeStamp<TList>("Centrality/Estimators", bc.timestamp());
      if (callst != nullptr) {
        auto getccdb = [callst](const char* ccdbhname) {
          TH1* h = reinterpret_cast<TH1*>(callst->FindObject(ccdbhname));
          return h;
        };
        auto getformulaccdb = [callst](const char* ccdbhname) {
          TFormula* f = reinterpret_cast<TFormula*>(callst->FindObject(ccdbhname));
          return f;
        };
        Run2V0MInfo.mhVtxAmpCorrV0A = getccdb("hVtx_fAmplitude_V0A_Normalized");
        Run2V0MInfo.mhVtxAmpCorrV0C = getccdb("hVtx_fAmplitude_V0C_Normalized");
        Run2V0MInfo.mhMultSelCalib = getccdb("hMultSelCalib_V0M");
        Run2V0MInfo.mMCScale = getformulaccdb(TString::Format("%s-V0M", genName->c_str()).Data());
        if ((Run2V0MInfo.mhVtxAmpCorrV0A != nullptr) && (Run2V0MInfo.mhVtxAmpCorrV0C != nullptr) && (Run2V0MInfo.mhMultSelCalib != nullptr)) {
          if (genName->length() != 0) {
            if (Run2V0MInfo.mMCScale != nullptr) {
              for (int ixpar = 0; ixpar < 6; ++ixpar) {
                Run2V0MInfo.mMCScalePars[ixpar] = Run2V0MInfo.mMCScale->GetParameter(ixpar);
              }
            } else {
              LOGF(fatal, "MC Scale information from V0M for run %d not available", bc.runNumber());
            }
          }
          Run2V0MInfo.mCalibrationStored = true;
        } else {
          LOGF(fatal, "Calibration information from V0M for run %d corrupted", bc.runNumber());
        }
        if (doprocessRun2) {
          Run2CL0Info.mhVtxAmpCorr = getccdb("hVtx_fnSPDClusters0_Normalized");
          Run2CL0Info.mhMultSelCalib = getccdb("hMultSelCalib_CL0");
          if ((Run2CL0Info.mhVtxAmpCorr != nullptr) && (Run2CL0Info.mhMultSelCalib != nullptr)) {
            Run2CL0Info.mCalibrationStored = true;
          } else {
            LOGF(fatal, "Calibration information from CL0 multiplicity for run %d corrupted", bc.runNumber());
          }
        }
      }
    } else {
      auto grpmagPath{"GLO/Config/GRPMagField"};
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
    }
    // Fetch magnetic field from ccdb for current collision
    d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << d_bz << " kG";
    mRunNumber = bc.runNumber();
    if (doprocessMiniRun2) {
      o2::ccdb::CcdbApi ccdbApi;
      ccdbApi.init("http://alice-ccdb.cern.ch");
      std::map<std::string, std::string> metadata;
      std::map<std::string, int>* classNameToIndexMap = ccdbApi.retrieveFromTFileAny<std::map<std::string, int>>("CTP/ClassNameToIndexMap", metadata, mRunNumber);
      for (const auto& classToIndexPair : *classNameToIndexMap) {
        bool hasClassName = classToIndexPair.first.find("HMV0M") < classToIndexPair.first.length();
        int classId = hasClassName ? classToIndexPair.second - 1 : -1;
        if (classId < 0) {
          continue;
        }
        classIds.push_back(classId);
      }
    }
    fitter.setBz(d_bz);

    // o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  template <class T>
  std::pair<float, float> getITSSignal(T const& track, aod::Run2TrackExtras const& trackExtraRun2)
  {
    if ((doprocessMiniRun2 || doprocessMiniMcRun2) && track.hasITS()) {
      auto extra = trackExtraRun2.rawIteratorAt(track.globalIndex());
      double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.p() / partMass[0]), cfgBetheBlochParamsITS->get("p0"), cfgBetheBlochParamsITS->get("p1"), cfgBetheBlochParamsITS->get("p2"), cfgBetheBlochParamsITS->get("p3"), cfgBetheBlochParamsITS->get("p4"))};
      double expSigma{expBethe * cfgBetheBlochParamsITS->get("resolution")};
      auto nSigmaITS = static_cast<float>((extra.itsSignal() - expBethe) / expSigma);
      return std::make_pair(extra.itsSignal(), nSigmaITS);
    }
    return std::make_pair(-999.f, -999.f);
  }

  template <class T>
  float getOuterPID(T const& track)
  {
    if ((doprocessMiniRun2 || doprocessMiniMcRun2) && track.hasTOF() && track.pt() > antipPtTof)
      return track.tofNSigmaPr();
    return -999.f;
  }

  float getV0M(int64_t const id, float const zvtx, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs)
  {
    auto fv0a = fv0as.rawIteratorAt(id);
    auto fv0c = fv0cs.rawIteratorAt(id);
    float multFV0A = 0;
    float multFV0C = 0;
    for (float amplitude : fv0a.amplitude()) {
      multFV0A += amplitude;
    }

    for (float amplitude : fv0c.amplitude()) {
      multFV0C += amplitude;
    }

    float v0m = -1;
    auto scaleMC = [](float x, float pars[6]) {
      return pow(((pars[0] + pars[1] * pow(x, pars[2])) - pars[3]) / pars[4], 1.0f / pars[5]);
    };

    if (Run2V0MInfo.mMCScale != nullptr) {
      float multFV0M = multFV0A + multFV0C;
      v0m = scaleMC(multFV0M, Run2V0MInfo.mMCScalePars);
      LOGF(debug, "Unscaled v0m: %f, scaled v0m: %f", multFV0M, v0m);
    } else if (Run2V0MInfo.mCalibrationStored) {
      v0m = multFV0A * Run2V0MInfo.mhVtxAmpCorrV0A->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0A->FindFixBin(zvtx)) +
            multFV0C * Run2V0MInfo.mhVtxAmpCorrV0C->GetBinContent(Run2V0MInfo.mhVtxAmpCorrV0C->FindFixBin(zvtx));
    }
    return v0m;
  }

  template <class T>
  int getTrackSelMask(T const& track)
  {
    int mask = 0x0;
    if (track.tpcncls > cfgTrackSels->get("tpcClsTight"))
      mask |= kTPCclsTight;
    else if (track.tpcncls > cfgTrackSels->get("tpcClsMid"))
      mask |= kTPCclsMid;
    if (track.tpcchi2 < cfgTrackSels->get("chi2TpcTight"))
      mask |= kChi2TPCTight;
    else if (track.tpcchi2 < cfgTrackSels->get("chi2TpcMid"))
      mask |= kChi2TPCMid;
    if (std::abs(track.dcaxypv) < cfgTrackSels->get("dcaxyTight") * (constDCASel ? 1. : dcaSigma(track.pt)))
      mask |= kDCAxyTight;
    else if (std::abs(track.dcaxypv) < cfgTrackSels->get("dcaxyMid") * (constDCASel ? 1. : dcaSigma(track.pt)))
      mask |= kDCAxyMid;
    if (std::abs(track.dcazpv) < cfgTrackSels->get("dcazTight"))
      mask |= kDCAzTight;
    else if (std::abs(track.dcazpv) < cfgTrackSels->get("dcazMid"))
      mask |= kDCAzMid;
    if (std::abs(track.tpcnsigma) < cfgTrackSels->get("tpcNsigmaTight"))
      mask |= kTPCPIDTight;
    else if (std::abs(track.tpcnsigma) < cfgTrackSels->get("tpcNsigmaMid"))
      mask |= kTPCPIDMid;
    if (std::abs(track.itsnsigma) < cfgTrackSels->get("itsNsigmaTight"))
      mask |= kITSPIDTight;
    else if (std::abs(track.itsnsigma) < cfgTrackSels->get("itsNsigmaMid"))
      mask |= kITSPIDMid;
    return mask;
  }

  void init(o2::framework::InitContext&)
  {

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    // lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(4);
    fitter.setMaxDXYIni(1);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
    int mat{static_cast<int>(cfgMaterialCorrection)};
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

    // event QA
    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    if (doprocessRun3) {
      histos.add<TH2>("QA/PvMultVsCent", ";Centrality T0C (%);#it{N}_{PV contributors};", HistType::kTH2F, {centAxis, multAxis});
      histos.add<TH2>("QA/MultVsCent", ";Centrality T0C (%);Multiplicity T0C;", HistType::kTH2F, {centAxis, multFt0Axis});
    } else if (doprocessRun2 || doprocessMiniRun2 || doprocessMcRun2 || doprocessMiniMcRun2) {
      histos.add<TH2>("QA/V0MvsCL0", ";Centrality CL0 (%);Centrality V0M (%)", HistType::kTH2F, {centAxis, centAxis});
      histos.add<TH2>("QA/trackletsVsV0M", ";Centrality CL0 (%);Centrality V0M (%)", HistType::kTH2F, {centAxis, multAxis});
      histos.add<TH2>("QA/nTrklCorrelation", ";Tracklets |#eta| < 0.6; Tracklets |#eta| > 0.7", HistType::kTH2D, {{201, -0.5, 200.5}, {201, -0.5, 200.5}});
      histos.add<TH1>("QA/TrklEta", ";Tracklets #eta; Entries", HistType::kTH1D, {{100, -3., 3.}});
    }

    // v0 QA
    histos.add<TH3>("QA/massLambda", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massLambdaAxis});

    // antid and antip QA
    histos.add<TH2>("QA/tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcAxis});
    histos.add<TH2>("QA/tpcSignalPr", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcAxis});
    histos.add<TH2>("QA/itsSignal", ";#it{p}_{glo} (GeV/#it{c});d#it{E}/d#it{x}_{ITS} (a.u.)", HistType::kTH2F, {momAxis, tpcAxis});
    tofMass[0] = histos.add<TH3>("QA/tofMass_p", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});
    tofMass[1] = histos.add<TH3>("QA/tofMass_d", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});

    ptMin = std::array<float, kNpart>{antipPtMin, antidPtMin};
    ptMax = std::array<float, kNpart>{antipPtMax, antidPtMax};
    ptTof = std::array<float, kNpart>{antipPtTof, antidPtTof};

    nSigmaTpcCutLow = std::array<float, kNpart>{antipNsigmaTpcCutLow, antidNsigmaTpcCutLow};
    nSigmaTpcCutUp = std::array<float, kNpart>{antipNsigmaTpcCutUp, antidNsigmaTpcCutUp};
    tpcInnerParamMax = std::array<float, kNpart>{antipTpcInnerParamMax, antidTpcInnerParamMax};
    tofMassMax = std::array<float, kNpart>{antipTofMassMax, antidTofMassMax};
  }

  template <class T>
  auto tracksSlice(T const& tracksAll, uint64_t const& collId)
  {
    if (doprocessRun3 || doprocessMcRun3)
      return tracksAll.sliceBy(perCollisionTracksFullIU, collId);
    else if (doprocessRun2 || doprocessMcRun2)
      return tracksAll.sliceBy(perCollisionTracksFull, collId);
    else
      return tracksAll.sliceBy(perCollisionTracksFullPID, collId);
  }

  template <class C, class T>
  void fillRecoEvent(C const& collision, T const& tracksAll, aod::V0s const& V0s, float const& centrality)
  {
    auto tracks = tracksSlice(tracksAll, collision.globalIndex());
    candidateTracks[0].clear();
    candidateTracks[1].clear();
    candidateV0s.clear();

    gpu::gpustd::array<float, 2> dcaInfo;
    uint8_t nTracklets[2]{0, 0};
    for (const auto& track : tracks) {

      if (track.trackType() == 255 && std::abs(track.eta()) < 1.2) { // tracklet
        if (std::abs(track.eta()) < 0.6)
          nTracklets[0]++;
        else if (std::abs(track.eta()) > 0.7)
          nTracklets[1]++;
      }

      if (!selectTrack(track)) {
        continue;
      }

      auto trackParCov = getTrackParCov(track);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      auto dca = std::hypot(dcaInfo[0], dcaInfo[1]);
      auto trackPt = trackParCov.getPt();
      auto trackEta = trackParCov.getEta();
      if (dca > cfgDcaSels->get("dca")) { // dca
        continue;
      }
      if (std::abs(dcaInfo[1]) > cfgDcaSels->get("dcaz")) { // dcaz
        continue;
      }
      if (std::abs(dcaInfo[0]) > cfgDcaSels->get("dcaxy") * (constDCASel ? 1. : dcaSigma(track.pt()))) { // dcaxy
        continue;
      }
      histos.fill(HIST("QA/tpcSignal"), track.tpcInnerParam(), track.tpcSignal());

      for (int iP{0}; iP < kNpart; ++iP) {
        if (trackPt < ptMin[iP] || trackPt > ptMax[iP]) {
          continue;
        }

        if (doprocessRun3 || doprocessMcRun3) {
          float cosL = 1 / std::sqrt(1.f + track.tgl() * track.tgl());
          if (iP && getITSClSize(track) * cosL < antidItsClsSizeCut && trackPt < antidPtItsClsSizeCut) {
            continue;
          }
        }

        double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / partMass[iP]), cfgBetheBlochParams->get(iP, "p0"), cfgBetheBlochParams->get(iP, "p1"), cfgBetheBlochParams->get(iP, "p2"), cfgBetheBlochParams->get(iP, "p3"), cfgBetheBlochParams->get(iP, "p4"))};
        double expSigma{expBethe * cfgBetheBlochParams->get(iP, "resolution")};
        auto nSigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);

        float beta{track.hasTOF() ? track.length() / (track.tofSignal() - track.tofEvTime()) * o2::constants::physics::invLightSpeedCm2PS : -999.f};
        beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta));
        float mass{track.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f)};
        bool hasTof = track.hasTOF() && track.tofChi2() < 3;

        if (trackPt <= ptTof[iP] || (trackPt > ptTof[iP] && hasTof && std::abs(mass - partMass[iP]) < tofMassMaxQA)) { // for QA histograms
          if (nSigmaTPC > nSigmaTpcCutLow[iP] && nSigmaTPC < nSigmaTpcCutUp[iP]) {
            tofMass[iP]->Fill(centrality, trackPt, mass);
          }
        }

        if (nSigmaTPC < nSigmaTpcCutLow[iP] || nSigmaTPC > nSigmaTpcCutUp[iP]) {
          continue;
        }

        // temporary cut to reject fake matches (run 3)
        if (track.tpcInnerParam() < tpcInnerParamMax[iP]) {
          continue;
        }
        if (trackPt > ptTof[iP] && !hasTof) {
          continue;
        }

        if (trackPt <= ptTof[iP] || (trackPt > ptTof[iP] && hasTof && std::abs(mass - partMass[iP]) < tofMassMax[iP])) {
          CandidateTrack candTrack;
          candTrack.pt = track.sign() > 0. ? trackPt : -trackPt;
          candTrack.eta = trackEta;
          candTrack.mass = iP;
          candTrack.dcapv = dca;
          candTrack.dcaxypv = dcaInfo[0];
          candTrack.dcazpv = dcaInfo[1];
          candTrack.tpcchi2 = track.tpcChi2NCl();
          candTrack.tpcncls = track.tpcNClsFound();
          candTrack.tpcnsigma = nSigmaTPC;
          candTrack.tofmass = hasTof ? mass : -999.f;
          candTrack.globalIndex = track.globalIndex();
          candTrack.outerPID = nSigmaTPC;
          candidateTracks[iP].push_back(candTrack);
        }
      }
    }
    if (doprocessRun2 || doprocessMcRun2 || doprocessMiniRun2 || doprocessMiniMcRun2) {
      histos.fill(HIST("QA/nTrklCorrelation"), nTracklets[0], nTracklets[1]);
      nTrackletsColl = nTracklets[1];
    }

    if (lambdaPtMax > lambdaPtMin) {
      std::vector<int64_t> trkId;
      for (const auto& v0 : V0s) {
        auto posTrack = v0.posTrack_as<T>();
        auto negTrack = v0.negTrack_as<T>();

        bool posSelect = selectV0Daughter(posTrack);
        bool negSelect = selectV0Daughter(negTrack);
        if (!posSelect || !negSelect)
          continue;

        if (doprocessRun2 || doprocessMiniRun2 || doprocessMcRun2 || doprocessMiniMcRun2) {
          bool checkPosPileUp = posTrack.hasTOF() || (posTrack.flags() & o2::aod::track::ITSrefit);
          bool checkNegPileUp = negTrack.hasTOF() || (negTrack.flags() & o2::aod::track::ITSrefit);
          if (!checkPosPileUp && !checkNegPileUp) {
            continue;
          }
        }

        auto posTrackCov = getTrackParCov(posTrack);
        auto negTrackCov = getTrackParCov(negTrack);

        int nCand = 0;
        try {
          nCand = fitter.process(posTrackCov, negTrackCov);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nCand == 0) {
          continue;
        }

        auto& posPropTrack = fitter.getTrack(0);
        auto& negPropTrack = fitter.getTrack(1);

        std::array<float, 3> momPos;
        std::array<float, 3> momNeg;
        std::array<float, 3> momV0;
        posPropTrack.getPxPyPzGlo(momPos);
        negPropTrack.getPxPyPzGlo(momNeg);
        momTotXYZ(momV0, momPos, momNeg);

        auto ptV0 = std::hypot(momV0[0], momV0[1]);
        if (ptV0 < lambdaPtMin || ptV0 > lambdaPtMax) {
          continue;
        }

        auto etaV0 = etaFromMom(momPos, momNeg);
        if (std::abs(etaV0) > etaMax) {
          continue;
        }

        auto alpha = alphaAP(momV0, momPos, momNeg);
        bool matter = alpha > 0;
        auto massPos = matter ? o2::constants::physics::MassProton : o2::constants::physics::MassPionCharged;
        auto massNeg = matter ? o2::constants::physics::MassPionCharged : o2::constants::physics::MassProton;
        auto mLambda = invMass2Body(momV0, momPos, momNeg, massPos, massNeg);
        auto mK0Short = invMass2Body(momV0, momPos, momNeg, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged);

        // pid selections
        double expBethePos{tpc::BetheBlochAleph(static_cast<double>(posTrack.tpcInnerParam() / massPos), cfgBetheBlochParams->get("p0"), cfgBetheBlochParams->get("p1"), cfgBetheBlochParams->get("p2"), cfgBetheBlochParams->get("p3"), cfgBetheBlochParams->get("p4"))};
        double expSigmaPos{expBethePos * cfgBetheBlochParams->get("resolution")};
        auto nSigmaTPCPos = static_cast<float>((posTrack.tpcSignal() - expBethePos) / expSigmaPos);
        double expBetheNeg{tpc::BetheBlochAleph(static_cast<double>(negTrack.tpcInnerParam() / massNeg), cfgBetheBlochParams->get("p0"), cfgBetheBlochParams->get("p1"), cfgBetheBlochParams->get("p2"), cfgBetheBlochParams->get("p3"), cfgBetheBlochParams->get("p4"))};
        double expSigmaNeg{expBetheNeg * cfgBetheBlochParams->get("resolution")};
        auto nSigmaTPCNeg = static_cast<float>((negTrack.tpcSignal() - expBetheNeg) / expSigmaNeg);
        float tpcSigPr = matter ? posTrack.tpcSignal() : negTrack.tpcSignal();

        if (std::abs(nSigmaTPCPos) > v0setting_nsigmatpc || std::abs(nSigmaTPCNeg) > v0setting_nsigmatpc) {
          continue;
        }

        // veto on K0s mass
        if (std::abs(mK0Short - o2::constants::physics::MassK0Short) < vetoMassK0Short) {
          continue;
        }

        float dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
        if (dcaV0dau > v0setting_dcav0dau) {
          continue;
        }

        std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
        const auto& vtx = fitter.getPCACandidate();

        float radiusV0 = std::hypot(vtx[0], vtx[1]);
        if (radiusV0 < v0setting_radius || radiusV0 > v0radiusMax) {
          continue;
        }

        float dcaV0Pv = CalculateDCAStraightToPV(
          vtx[0], vtx[1], vtx[2],
          momPos[0] + momNeg[0],
          momPos[1] + momNeg[1],
          momPos[2] + momNeg[2],
          collision.posX(), collision.posY(), collision.posZ());
        if (std::abs(dcaV0Pv) > v0setting_dcav0pv) {
          continue;
        }

        double cosPA = RecoDecay::cpa(primVtx, vtx, momV0);
        if (cosPA < v0setting_cospa) {
          continue;
        }

        auto ptotal = RecoDecay::sqrtSumOfSquares(momV0[0], momV0[1], momV0[2]);
        auto lengthTraveled = RecoDecay::sqrtSumOfSquares(vtx[0] - primVtx[0], vtx[1] - primVtx[1], vtx[2] - primVtx[2]);
        float ML2P_Lambda = o2::constants::physics::MassLambda * lengthTraveled / ptotal;
        if (ML2P_Lambda > v0setting_lifetime) {
          continue;
        }

        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, posTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
        auto posDcaToPv = std::hypot(dcaInfo[0], dcaInfo[1]);
        if (posDcaToPv < v0setting_dcadaughtopv && std::abs(dcaInfo[0]) < v0setting_dcadaughtopv) {
          continue;
        }

        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, negTrackCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
        auto negDcaToPv = std::hypot(dcaInfo[0], dcaInfo[1]);
        if (negDcaToPv < v0setting_dcadaughtopv && std::abs(dcaInfo[0]) < v0setting_dcadaughtopv) {
          continue;
        }

        if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) { // for QA histograms
          continue;
        }
        histos.fill(HIST("QA/massLambda"), centrality, ptV0, mLambda);

        histos.fill(HIST("QA/tpcSignalPr"), matter > 0. ? posTrack.tpcInnerParam() : negTrack.tpcInnerParam(), tpcSigPr);

        CandidateV0 candV0;
        candV0.pt = matter > 0. ? ptV0 : -ptV0;
        candV0.eta = etaV0;
        candV0.mass = mLambda;
        candV0.cpa = cosPA;
        candV0.dcav0daugh = dcaV0dau;
        candV0.dcav0pv = dcaV0Pv;
        candV0.dcanegpv = negDcaToPv;
        candV0.dcapospv = posDcaToPv;
        candV0.tpcnsigmaneg = nSigmaTPCNeg;
        candV0.tpcnsigmapos = nSigmaTPCPos;
        candV0.globalIndexPos = posTrack.globalIndex();
        candV0.globalIndexNeg = negTrack.globalIndex();
        candidateV0s.push_back(candV0);
      }
    }
  }

  template <class C, class T>
  void fillMcEvent(C const& collision, T const& tracks, aod::V0s const& V0s, float const& centrality, aod::McParticles const& particlesMC, aod::McTrackLabels const& mcLabels)
  {
    fillRecoEvent<C, T>(collision, tracks, V0s, centrality);

    for (int iP{0}; iP < kNpart; ++iP) {
      for (auto& candidateTrack : candidateTracks[iP]) {
        candidateTrack.isreco = true;

        auto mcLab = mcLabels.rawIteratorAt(candidateTrack.globalIndex);

        if (mcLab.mcParticleId() < -1 || mcLab.mcParticleId() >= particlesMC.size()) {
          continue;
        }
        if (mcLab.has_mcParticle()) {
          auto mcTrack = mcLab.template mcParticle_as<aod::McParticles>();
          if (std::abs(mcTrack.pdgCode()) != partPdg[iP])
            continue;
          if (((mcTrack.flags() & 0x8) && (doprocessMcRun2 || doprocessMiniMcRun2)) || (mcTrack.flags() & 0x2) || ((mcTrack.flags() & 0x1) && !doprocessMiniMcRun2))
            continue;

          if (!mcTrack.isPhysicalPrimary() && !doprocessMiniMcRun2)
            continue;
          if (mcTrack.isPhysicalPrimary())
            candidateTrack.pdgcodemoth = PartTypes::kPhysPrim;
          else if (mcTrack.has_mothers() && iP == 0)
            candidateTrack.pdgcodemoth = getPartTypeMother(mcTrack);

          auto genPt = std::hypot(mcTrack.px(), mcTrack.py());
          candidateTrack.pdgcode = mcTrack.pdgCode();
          candidateTrack.genpt = genPt;
          candidateTrack.geneta = mcTrack.eta();
          candidateTrack.mcIndex = mcTrack.globalIndex();
        }
      }
    }
    for (auto& candidateV0 : candidateV0s) {
      candidateV0.isreco = true;
      auto mcLabPos = mcLabels.rawIteratorAt(candidateV0.globalIndexPos);
      auto mcLabNeg = mcLabels.rawIteratorAt(candidateV0.globalIndexNeg);

      if (mcLabPos.has_mcParticle() && mcLabNeg.has_mcParticle()) {
        auto mcTrackPos = mcLabPos.template mcParticle_as<aod::McParticles>();
        auto mcTrackNeg = mcLabNeg.template mcParticle_as<aod::McParticles>();
        if (mcTrackPos.has_mothers() && mcTrackNeg.has_mothers()) {
          for (auto& negMother : mcTrackNeg.template mothers_as<aod::McParticles>()) {
            for (auto& posMother : mcTrackPos.template mothers_as<aod::McParticles>()) {
              if (posMother.globalIndex() != negMother.globalIndex())
                continue;
              if (!((mcTrackPos.pdgCode() == 2212 && mcTrackNeg.pdgCode() == -211) || (mcTrackPos.pdgCode() == 211 && mcTrackNeg.pdgCode() == -2212)))
                continue;
              if (std::abs(posMother.pdgCode()) != 3122) {
                continue;
              }
              if (!posMother.isPhysicalPrimary() && !posMother.has_mothers())
                continue;
              if (((posMother.flags() & 0x8) && (doprocessMcRun2 || doprocessMiniMcRun2)) || (posMother.flags() & 0x2) || (posMother.flags() & 0x1))
                continue;

              auto genPt = std::hypot(posMother.px(), posMother.py());
              candidateV0.pdgcode = posMother.pdgCode();
              candidateV0.genpt = genPt;
              candidateV0.geneta = posMother.eta();
              candidateV0.mcIndex = posMother.globalIndex();
            }
          }
        }
      }
    }
  }

  void fillMcGen(aod::McParticles const& mcParticles, aod::McTrackLabels const& /*mcLab*/, uint64_t const& collisionId)
  {
    auto mcParticles_thisCollision = mcParticles.sliceBy(perCollisionMcParts, collisionId);
    for (auto& mcPart : mcParticles_thisCollision) {
      auto genEta = mcPart.eta();
      if (std::abs(genEta) > etaMax) {
        continue;
      }
      if (((mcPart.flags() & 0x8) && (doprocessMcRun2 || doprocessMiniMcRun2)) || (mcPart.flags() & 0x2) || ((mcPart.flags() & 0x1) && !doprocessMiniMcRun2))
        continue;
      auto pdgCode = mcPart.pdgCode();
      if (std::abs(pdgCode) == 3122) {
        if (!mcPart.isPhysicalPrimary() && !mcPart.has_mothers())
          continue;
        bool foundPr = false;
        for (auto& mcDaught : mcPart.daughters_as<aod::McParticles>()) {
          if (std::abs(mcDaught.pdgCode()) == 2212) {
            foundPr = true;
            break;
          }
        }
        if (!foundPr) {
          continue;
        }
        auto genPt = std::hypot(mcPart.px(), mcPart.py());
        CandidateV0 candV0;
        candV0.genpt = genPt;
        candV0.geneta = mcPart.eta();
        candV0.pdgcode = pdgCode;
        auto it = find_if(candidateV0s.begin(), candidateV0s.end(), [&](CandidateV0 v0) { return v0.mcIndex == mcPart.globalIndex(); });
        if (it != candidateV0s.end()) {
          continue;
        } else {
          LOGF(debug, "not found!");
          candidateV0s.emplace_back(candV0);
        }
      } else if (std::abs(pdgCode) == partPdg[0] || std::abs(pdgCode) == partPdg[1]) {
        int iP = 1;
        if (std::abs(pdgCode) == partPdg[0]) {
          iP = 0;
        }
        if ((!mcPart.isPhysicalPrimary() && !doprocessMiniMcRun2))
          continue;
        auto genPt = std::hypot(mcPart.px(), mcPart.py());
        CandidateTrack candTrack;
        candTrack.genpt = genPt;
        candTrack.geneta = mcPart.eta();
        candTrack.pdgcode = pdgCode;
        if (mcPart.isPhysicalPrimary())
          candTrack.pdgcodemoth = PartTypes::kPhysPrim;
        else if (mcPart.has_mothers() && iP == 0)
          candTrack.pdgcodemoth = getPartTypeMother(mcPart);

        auto it = find_if(candidateTracks[iP].begin(), candidateTracks[iP].end(), [&](CandidateTrack trk) { return trk.mcIndex == mcPart.globalIndex(); });
        if (it != candidateTracks[iP].end()) {
          continue;
        } else {
          candidateTracks[iP].emplace_back(candTrack);
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults> const& collisions, TracksFullIU const& tracks, aod::V0s const& V0s, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      auto multiplicity = collision.multFT0C();
      auto centrality = collision.centFT0C();
      fillRecoEvent(collision, tracks, V0Table_thisCollision, centrality);

      histos.fill(HIST("QA/PvMultVsCent"), centrality, collision.numContrib());
      histos.fill(HIST("QA/MultVsCent"), centrality, multiplicity);

      collisionEbyeTable(centrality, collision.posZ());

      for (auto& candidateV0 : candidateV0s) {
        lambdaEbyeTable(
          collisionEbyeTable.lastIndex(),
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.dcav0pv,
          // candidateV0.dcanegpv,
          // candidateV0.dcapospv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          // candidateV0.tpcnsigmaneg,
          // candidateV0.tpcnsigmapos,
          candidateV0.globalIndexNeg,
          candidateV0.globalIndexPos);
      }

      for (int iP{0}; iP < kNpart; ++iP) {
        for (auto& candidateTrack : candidateTracks[iP]) { // deuterons + protons
          nucleiEbyeTable(
            collisionEbyeTable.lastIndex(),
            candidateTrack.pt,
            candidateTrack.eta,
            candidateTrack.mass,
            candidateTrack.dcapv,
            candidateTrack.tpcncls,
            candidateTrack.tpcnsigma,
            candidateTrack.tofmass);
        }
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processRun3, "process (Run 3)", false);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::TrackletMults> const& collisions, TracksFull const& tracks, aod::V0s const& V0s, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs, BCsWithRun2Info const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      if (kUseTPCPileUpCut && !(bc.eventCuts() & BIT(aod::Run2EventCuts::kTPCPileUp)))
        continue;

      float v0m = getV0M(bc.globalIndex(), collision.posZ(), fv0as, fv0cs);
      float cV0M = 105.f;
      if (Run2V0MInfo.mCalibrationStored) {
        cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
        if (!(collision.sel7() && collision.alias_bit(kINT7)) && (!kINT7Intervals || (kINT7Intervals && ((cV0M >= 10 && cV0M < 30) || cV0M > 50))))
          continue;
      }

      auto centralityCl0 = 105.0f;
      if (Run2CL0Info.mCalibrationStored) {
        float cl0m = bc.spdClustersL0() * Run2CL0Info.mhVtxAmpCorr->GetBinContent(Run2CL0Info.mhVtxAmpCorr->FindFixBin(collision.posZ()));
        centralityCl0 = Run2CL0Info.mhMultSelCalib->GetBinContent(Run2CL0Info.mhMultSelCalib->FindFixBin(cl0m));
      }
      if (kUseEstimatorsCorrelationCut) {
        const auto& x = centralityCl0;
        const double center = estimatorsCorrelationCoef[0] + estimatorsCorrelationCoef[1] * x;
        const double sigma = estimatorsSigmaPars[0] + estimatorsSigmaPars[1] * x + estimatorsSigmaPars[2] * std::pow(x, 2) + estimatorsSigmaPars[3] * std::pow(x, 3);
        if (cV0M < center - deltaEstimatorNsigma[0] * sigma || cV0M > center + deltaEstimatorNsigma[1] * sigma) {
          continue;
        }
      }

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      auto multTracklets = collision.multTracklets();
      fillRecoEvent(collision, tracks, V0Table_thisCollision, cV0M);

      histos.fill(HIST("QA/V0MvsCL0"), centralityCl0, cV0M);
      histos.fill(HIST("QA/trackletsVsV0M"), cV0M, multTracklets);

      collisionEbyeTable(cV0M, collision.posZ());

      for (auto& candidateV0 : candidateV0s) {
        lambdaEbyeTable(
          collisionEbyeTable.lastIndex(),
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.dcav0pv,
          // candidateV0.dcanegpv,
          // candidateV0.dcapospv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          // candidateV0.tpcnsigmaneg,
          // candidateV0.tpcnsigmapos,
          candidateV0.globalIndexNeg,
          candidateV0.globalIndexPos);
      }

      for (int iP{0}; iP < kNpart; ++iP) {
        for (auto& candidateTrack : candidateTracks[iP]) { // deuterons + protons
          nucleiEbyeTable(
            collisionEbyeTable.lastIndex(),
            candidateTrack.pt,
            candidateTrack.eta,
            candidateTrack.mass,
            candidateTrack.dcapv,
            candidateTrack.tpcncls,
            candidateTrack.tpcnsigma,
            candidateTrack.tofmass);
        }
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processRun2, "process (Run 2)", false);

  void processMiniRun2(soa::Join<aod::Collisions, aod::EvSels> const& collisions, TracksFullPID const& tracks, aod::Run2TrackExtras const& trackExtraRun2, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs, aod::V0s const& V0s, BCsWithRun2Info const&)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kINELgtZERO)))
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      float v0m = getV0M(bc.globalIndex(), collision.posZ(), fv0as, fv0cs);
      float cV0M = 105.f;
      if (Run2V0MInfo.mCalibrationStored) {
        cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillRecoEvent(collision, tracks, V0Table_thisCollision, cV0M);

      uint8_t trigger = collision.alias_bit(kINT7) ? 0x1 : 0x0;
      for (auto& classId : classIds) {
        if (bc.triggerMask() & BIT(classId)) {
          trigger |= 0x2;
          cV0M = cV0M < 104.f ? cV0M * 100. : cV0M;
          break;
        }
      }
      if (trigger == 0x0) {
        continue;
      }
      if (triggerCut != 0x0 && (trigger & triggerCut) != triggerCut) {
        continue;
      }
      miniCollTable(static_cast<int8_t>(collision.posZ() * 10), trigger, nTrackletsColl, cV0M);

      for (auto& candidateTrack : candidateTracks[0]) { // protons
        auto tk = tracks.rawIteratorAt(candidateTrack.globalIndex);
        float outerPID = getOuterPID(tk);
        auto [itsSignal, nSigmaITS] = getITSSignal(tk, trackExtraRun2);
        histos.fill(HIST("QA/itsSignal"), tk.p(), itsSignal);
        candidateTrack.itsnsigma = nSigmaITS;
        candidateTrack.outerPID = tk.pt() < antipPtTof ? candidateTrack.outerPID : outerPID;
        int selMask = getTrackSelMask(candidateTrack);
        if (candidateTrack.outerPID < outerPIDMin)
          continue;
        miniTrkTable(
          miniCollTable.lastIndex(),
          candidateTrack.pt,
          static_cast<int8_t>(candidateTrack.eta * 100),
          selMask,
          candidateTrack.outerPID);
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processMiniRun2, "process mini tables(Run 2)", false);

  void processMcRun3(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::McCollisions const& /*mcCollisions*/, TracksFullIU const& tracks, aod::V0s const& V0s, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      auto centrality = collision.centFT0C();

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillMcEvent(collision, tracks, V0Table_thisCollision, centrality, mcParticles, mcLab);
      fillMcGen(mcParticles, mcLab, collision.mcCollisionId());

      collisionEbyeTable(centrality, collision.posZ());

      for (auto& candidateV0 : candidateV0s) {
        mcLambdaEbyeTable(
          collisionEbyeTable.lastIndex(),
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.dcav0pv,
          // candidateV0.dcanegpv,
          // candidateV0.dcapospv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          // candidateV0.tpcnsigmaneg,
          // candidateV0.tpcnsigmapos,
          candidateV0.globalIndexNeg,
          candidateV0.globalIndexPos,
          candidateV0.genpt,
          candidateV0.geneta,
          candidateV0.pdgcode,
          candidateV0.isreco);
      }

      for (int iP{0}; iP < kNpart; ++iP) {
        for (auto& candidateTrack : candidateTracks[iP]) { // deuterons + protons
          mcNucleiEbyeTable(
            collisionEbyeTable.lastIndex(),
            candidateTrack.pt,
            candidateTrack.eta,
            candidateTrack.mass,
            candidateTrack.dcapv,
            candidateTrack.tpcncls,
            candidateTrack.tpcnsigma,
            candidateTrack.tofmass,
            candidateTrack.genpt,
            candidateTrack.geneta,
            candidateTrack.pdgcode,
            candidateTrack.isreco);
        }
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processMcRun3, "process MC (Run 3)", false);

  void processMcRun2(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const& /*mcCollisions*/, TracksFull const& tracks, aod::V0s const& V0s, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, BCsWithRun2Info const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      float v0m = getV0M(bc.globalIndex(), collision.posZ(), fv0as, fv0cs);
      float cV0M = 105.f;
      if (Run2V0MInfo.mCalibrationStored) {
        cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillMcEvent(collision, tracks, V0Table_thisCollision, cV0M, mcParticles, mcLab);
      fillMcGen(mcParticles, mcLab, collision.mcCollisionId());

      collisionEbyeTable(cV0M, collision.posZ());

      for (auto& candidateV0 : candidateV0s) {
        mcLambdaEbyeTable(
          collisionEbyeTable.lastIndex(),
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.dcav0pv,
          // candidateV0.dcanegpv,
          // candidateV0.dcapospv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          // candidateV0.tpcnsigmaneg,
          // candidateV0.tpcnsigmapos,
          candidateV0.globalIndexNeg,
          candidateV0.globalIndexPos,
          candidateV0.genpt,
          candidateV0.geneta,
          candidateV0.pdgcode,
          candidateV0.isreco);
      }

      for (int iP{0}; iP < kNpart; ++iP) {
        for (auto& candidateTrack : candidateTracks[iP]) { // deuterons + protons
          mcNucleiEbyeTable(
            collisionEbyeTable.lastIndex(),
            candidateTrack.pt,
            candidateTrack.eta,
            candidateTrack.mass,
            candidateTrack.dcapv,
            candidateTrack.tpcncls,
            candidateTrack.tpcnsigma,
            candidateTrack.tofmass,
            candidateTrack.genpt,
            candidateTrack.geneta,
            candidateTrack.pdgcode,
            candidateTrack.isreco);
        }
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processMcRun2, "process MC (Run 2)", false);

  void processMiniMcRun2(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const& /*mcCollisions*/, TracksFullPID const& tracks, aod::Run2TrackExtras const& trackExtraRun2, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs, aod::V0s const& V0s, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, BCsWithRun2Info const&)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      float v0m = getV0M(bc.globalIndex(), collision.posZ(), fv0as, fv0cs);
      float cV0M = 105.f;
      if (Run2V0MInfo.mCalibrationStored) {
        cV0M = Run2V0MInfo.mhMultSelCalib->GetBinContent(Run2V0MInfo.mhMultSelCalib->FindFixBin(v0m));
      }

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillMcEvent(collision, tracks, V0Table_thisCollision, cV0M, mcParticles, mcLab);
      fillMcGen(mcParticles, mcLab, collision.mcCollisionId());

      miniCollTable(static_cast<int8_t>(collision.posZ() * 10), 0x0, nTrackletsColl, cV0M);

      for (auto& candidateTrack : candidateTracks[0]) { // protons
        int selMask = -1;
        if (candidateTrack.isreco) {
          auto tk = tracks.rawIteratorAt(candidateTrack.globalIndex);
          float outerPID = getOuterPID(tk);
          auto [itsSignal, nSigmaITS] = getITSSignal(tk, trackExtraRun2);
          histos.fill(HIST("QA/itsSignal"), tk.p(), itsSignal);
          candidateTrack.itsnsigma = nSigmaITS;
          candidateTrack.outerPID = tk.pt() < antipPtTof ? candidateTrack.outerPID : outerPID;
          selMask = getTrackSelMask(candidateTrack);
          // if (candidateTrack.outerPID < -4)
          //   continue;
          if (candidateTrack.pdgcodemoth > 0)
            selMask |= candidateTrack.pdgcodemoth;
        } else if (candidateTrack.pdgcodemoth > 0) {
          selMask = candidateTrack.pdgcodemoth;
        }
        if (selMask < 0)
          continue;
        mcMiniTrkTable(
          miniCollTable.lastIndex(),
          candidateTrack.pt,
          static_cast<int8_t>(candidateTrack.eta * 100),
          selMask,
          candidateTrack.outerPID,
          candidateTrack.pdgcode > 0 ? candidateTrack.genpt : -candidateTrack.genpt,
          static_cast<int8_t>(candidateTrack.geneta * 100),
          candidateTrack.isreco);
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processMiniMcRun2, "process mini tables for mc(Run 2)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ebyeMaker>(cfgc)};
}
