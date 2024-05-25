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
#include <utility>
#include <random>

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

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFEbyeTables.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime>;
using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TOFSignal, aod::TOFEvTime>;
using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

namespace
{
constexpr int kNpart = 2;
constexpr double betheBlochDefault[kNpart][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}, {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
constexpr double estimatorsCorrelationCoef[2]{-0.669108, 1.04489};
constexpr double estimatorsSigmaPars[4]{0.933321, 0.0416976, -0.000936344, 8.92179e-06};
constexpr double deltaEstimatorNsigma[2]{5.5, 5.};
constexpr double partMass[kNpart]{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron};
constexpr double partPdg[kNpart]{2212, o2::constants::physics::kDeuteron};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNamesBB{"p", "d"};
std::array<std::shared_ptr<TH2>, kNpart> tempTracks;
std::shared_ptr<TH2> tempAntiLambda;
std::shared_ptr<TH2> tempLambda;
std::array<std::shared_ptr<TH3>, kNpart> recTracks;
std::array<std::shared_ptr<TH3>, kNpart> recAntiTracks;
std::array<std::shared_ptr<TH3>, kNpart> genTracks;
std::array<std::shared_ptr<TH3>, kNpart> genAntiTracks;
std::array<std::shared_ptr<TH2>, kNpart> tpcNsigma;
std::array<std::shared_ptr<TH3>, kNpart> tpcNsigmaGlo;
std::array<std::shared_ptr<TH3>, kNpart> tofMass;
std::array<std::shared_ptr<TH2>, kNpart> tofSignal;
std::array<std::shared_ptr<TH2>, kNpart> tofSignal_glo;
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
  bool matter = 0;
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
  float mass = -999.f;
  bool matter = -999.f;
  float dcapv = 0;
  uint8_t tpcncls = 0;
  float tpcnsigma = -999.f;
  float tofmass = -999.f;
  float genpt = -999.f;
  float geneta = -999.f;
  int pdgcode = -999;
  bool isreco = 0;
  int64_t mcIndex = -999;
  int64_t globalIndex = -999;
};

struct ebyeMaker {
  Produces<aod::NucleiEbyeTable> nucleiEbyeTable;
  Produces<aod::LambdaEbyeTable> lambdaEbyeTable;
  Produces<aod::McNucleiEbyeTable> mcNucleiEbyeTable;
  Produces<aod::McLambdaEbyeTable> mcLambdaEbyeTable;
  std::mt19937 gen32;
  std::vector<CandidateV0> candidateV0s;
  std::array<std::vector<CandidateTrack>, 2> candidateTracks;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;

  int nSubsamples;
  int mRunNumber;
  float d_bz;
  // o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 2, 6, particleNamesBB, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for deuteron"};

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis subsampleAxis{"subsampleAxis", {30, 0, 30}, "binning of the subsample axis"};
  ConfigurableAxis deltaEtaAxis{"deltaEtaAxis", {4, 0, 0.8}, "binning of the delta eta axis"};
  ConfigurableAxis ptAntidAxis{"ptAntidAxis", {VARIABLE_WIDTH, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f}, "binning of the antideuteron pT axis (GeV/c)"};
  ConfigurableAxis ptAntipAxis{"ptAntipAxis", {VARIABLE_WIDTH, 0.4f, 0.6f, 0.7f, 0.8f, 0.9f}, "binning of the antiproton pT axis (GeV/c)"};
  ConfigurableAxis ptLambdaAxis{"ptLambdaAxis", {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f}, "binning of the (anti)lambda pT axis (GeV/c)"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis multAxis{"multAxis", {100, 0, 10000}, "Binning for the multiplicity axis"};
  ConfigurableAxis multFt0Axis{"multFt0Axis", {100, 0, 100000}, "Binning for the ft0 multiplicity axis"};
  ConfigurableAxis nGenRecAxis{"nGenRecAxis", {20, 0, 20}, "binning for the number of reconstructed or generated candidates per event"};

  // binning of (anti)lambda QA histograms
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.03f, o2::constants::physics::MassLambda0 + 0.03f}, "binning for the lambda invariant-mass"};
  ConfigurableAxis cosPaAxis{"cosPaAxis", {1e3, 0.95f, 1.00f}, "binning for the cosPa axis"};
  ConfigurableAxis radiusAxis{"radiusAxis", {1e3, 0.f, 100.f}, "binning for the radius axis"};
  ConfigurableAxis dcaV0daughAxis{"dcaV0daughAxis", {2e2, 0.f, 2.f}, "binning for the dca of V0 daughters"};
  ConfigurableAxis dcaDaughPvAxis{"dcaDaughPvAxis", {1e3, -10.f, 10.f}, "binning for the dca of positive daughter to PV"};

  // binning of deuteron QA histograms
  ConfigurableAxis tpcNsigmaAxis{"tpcNsigmaAxis", {100, -5.f, 5.f}, "tpc nsigma axis"};
  ConfigurableAxis tofMassAxis{"tofMassAxis", {1000, 0., 3.f}, "tof mass axis"};
  ConfigurableAxis momAxis{"momAxis", {60., 0.f, 3.f}, "momentum axis binning"};
  ConfigurableAxis momAxisFine{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};
  ConfigurableAxis momResAxis{"momResAxis", {1.e2, -1.f, 1.f}, "momentum resolution binning"};
  ConfigurableAxis tpcAxis{"tpcAxis", {4.e2, 0.f, 4.e3f}, "tpc signal axis binning"};
  ConfigurableAxis tofAxis{"tofAxis", {1.e3, 0.f, 1.f}, "tof signal axis binning"};
  ConfigurableAxis tpcClsAxis{"tpcClsAxis", {160, 0.f, 160.f}, "tpc n clusters binning"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};
  Configurable<float> etaMaxV0dau{"etaMaxV0dau", 0.8f, "maximum eta V0 daughters"};

  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "fill histograms only for true signal candidates (MC)"};

  Configurable<bool> kINT7Intervals{"kINT7Intervals", false, "toggle kINT7 trigger selection in the 10-30% and 50-90% centrality intervals (2018 Pb-Pb)"};
  Configurable<bool> kUseTPCPileUpCut{"kUseTPCPileUpCut", false, "toggle strong correlation cuts (Run 2)"};
  Configurable<bool> kUseEstimatorsCorrelationCut{"kUseEstimatorsCorrelationCut", false, "toggle cut on the correlation between centrality estimators (2018 Pb-Pb)"};

  Configurable<float> antidPtMin{"antidPtMin", 0.8f, "minimum antideuteron pT (GeV/c)"};
  Configurable<float> antidPtTof{"antidPtTof", 1.0f, "antideuteron pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antidPtMax{"antidPtMax", 1.8f, "maximum antideuteron pT (GeV/c)"};

  Configurable<float> antipPtMin{"antipPtMin", 0.4f, "minimum antiproton pT (GeV/c)"};
  Configurable<float> antipPtTof{"antipPtTof", 0.6f, "antiproton pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antipPtMax{"antipPtMax", 0.9f, "maximum antiproton pT (GeV/c)"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 3.0f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> trackNcrossedRows{"trackNcrossedRows", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> trackNclusItsCut{"trackNclusITScut", 5, "Minimum number of ITS clusters"};
  Configurable<float> trackNclusTpcCut{"trackNclusTPCcut", 70, "Minimum number of TPC clusters"};
  Configurable<float> trackDcaCut{"trackDcaCut", 0.1f, "DCA antid to PV"};

  Configurable<float> v0trackNcrossedRows{"v0trackNcrossedRows", 70, "Minimum number of crossed TPC rows for V0 daughter"};
  Configurable<float> v0trackNclusItsCut{"v0trackNclusITScut", 1, "Minimum number of ITS clusters for V0 daughter"};
  Configurable<float> v0trackNclusTpcCut{"v0trackNclusTPCcut", 70, "Minimum number of TPC clusters for V0 daughter"};
  Configurable<float> v0trackNsharedClusTpc{"v0trackNsharedClusTpc", 10, "Maximum number of shared TPC clusters for V0 daughter"};
  Configurable<bool> v0requireITSrefit{"v0requireITSrefit", false, "require ITS refit for V0 daughter"};
  Configurable<float> vetoMassK0Short{"vetoMassK0Short", -999.f, "veto for V0 compatible with K0s mass"};
  Configurable<float> v0radiusMax{"v0radiusMax", 100.f, "maximum V0 radius eccepted"};

  Configurable<float> antidNsigmaTpcCutLow{"antidNsigmaTpcCutLow", 4.f, "TPC PID cut low"};
  Configurable<float> antidNsigmaTpcCutUp{"antidNsigmaTpcCutUp", 4.f, "TPC PID cut up"};
  Configurable<float> antidTpcInnerParamMax{"tpcInnerParamMax", 0.6f, "(temporary) tpc inner param cut"};
  Configurable<float> antidTofMassMax{"tofMassMax", 0.3f, "(temporary) tof mass cut"};

  Configurable<float> antipNsigmaTpcCutLow{"antipNsigmaTpcCutLow", 4.f, "TPC PID cut low"};
  Configurable<float> antipNsigmaTpcCutUp{"antipNsigmaTpcCutUp", 4.f, "TPC PID cut up"};
  Configurable<float> antipTpcInnerParamMax{"antipTpcInnerParamMax", 0.6f, "(temporary) tpc inner param cut"};
  Configurable<float> antipTofMassMax{"antipTofMassMax", 0.3f, "(temporary) tof mass cut"};
  Configurable<float> tofMassMaxQA{"tofMassMaxQA", 0.6f, "(temporary) tof mass cut (for QA histograms)"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcav0pv{"v0setting_dcav0pv", 1, "DCA V0 to Pv"};
  Configurable<float> v0setting_dcadaughtopv{"v0setting_dcadaughtopv", 0.1f, "DCA Pos To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5f, "v0radius"};
  Configurable<float> v0setting_lifetime{"v0setting_lifetime", 40.f, "v0 lifetime cut"};
  Configurable<float> v0setting_nsigmatpc{"v0setting_nsigmatpc", 4.f, "nsigmatpc"};
  Configurable<float> lambdaMassCutQA{"lambdaMassCutQA", 0.02f, "maximum deviation from PDG mass (for QA histograms)"};

  Configurable<float> antidItsClsSizeCut{"antidItsClsSizeCut", 2.f, "cluster size cut for antideuterons"};
  Configurable<float> antidPtItsClsSizeCut{"antidPtItsClsSizeCut", 1.f, "pt for cluster size cut for antideuterons"};

  std::array<float, kNpart> ptMin;
  std::array<float, kNpart> ptTof;
  std::array<float, kNpart> ptMax;
  std::array<float, kNpart> nSigmaTpcCutLow;
  std::array<float, kNpart> nSigmaTpcCutUp;
  std::array<float, kNpart> tpcInnerParamMax;
  std::array<float, kNpart> tofMassMax;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry tempHistos{"tempHistos", {}, OutputObjHandlingPolicy::TransientObject};

  Preslice<TracksFull> perCollisionTracksFull = o2::aod::track::collisionId;
  Preslice<TracksFullIU> perCollisionTracksFullIU = o2::aod::track::collisionId;
  Preslice<aod::V0s> perCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::McParticles> perCollisionMcParts = o2::aod::mcparticle::mcCollisionId;

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
    if (doprocessRun2 || doprocessMcRun2) {
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
        track.tpcChi2NCl() > 4.f ||
        track.itsChi2NCl() > 36.f) {
      return false;
    }
    if (doprocessRun2 || doprocessMcRun2) {
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

  void fillHistoN(std::shared_ptr<THnSparse> hFull, std::shared_ptr<TH2> const& hTmp, int const subsample, int const centrality)
  {
    for (int iEta{1}; iEta < hTmp->GetNbinsX() + 1; ++iEta) {
      for (int iPt{1}; iPt < hTmp->GetNbinsY() + 1; ++iPt) {
        auto eta = hTmp->GetXaxis()->GetBinCenter(iEta);
        auto pt = hTmp->GetYaxis()->GetBinCenter(iPt);
        auto num = hTmp->Integral(1, iEta, iPt, iPt);

        hFull->Fill(subsample, centrality, eta, pt, num);
      }
    }
  }

  void fillHistoN(std::shared_ptr<THnSparse> hFull, std::shared_ptr<TH2> const& hTmpA, std::shared_ptr<TH2> const& hTmpB, int const subsample, int const centrality)
  {
    for (int iEta{1}; iEta < hTmpA->GetNbinsX() + 1; ++iEta) {
      auto eta = hTmpA->GetXaxis()->GetBinCenter(iEta);
      for (int iPtA{1}; iPtA < hTmpA->GetNbinsY() + 1; ++iPtA) {
        for (int iPtB{1}; iPtB < hTmpB->GetNbinsY() + 1; ++iPtB) {
          auto ptA = hTmpA->GetYaxis()->GetBinCenter(iPtA);
          auto ptB = hTmpB->GetYaxis()->GetBinCenter(iPtB);
          auto numA = hTmpA->Integral(1, iEta, iPtA, iPtA);
          auto numB = hTmpB->Integral(1, iEta, iPtB, iPtB);

          hFull->Fill(subsample, centrality, eta, ptA, ptB, numA * numB);
        }
      }
    }
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (doprocessRun2 || doprocessMcRun2) {
      auto grpPath{"GLO/GRP/GRP"};
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (!grpo) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpPath << " of object GRPObject for timestamp " << timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
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
    fitter.setBz(d_bz);

    // o2::base::Propagator::Instance()->setMatLUT(lut);
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

    uint32_t randomSeed = static_cast<uint32_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    gen32.seed(randomSeed);

    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    auto hNev = histos.add<THnSparse>("nEv", ";Subsample;Centrality (%);", HistType::kTHnSparseD, {subsampleAxis, centAxis});
    nSubsamples = hNev->GetAxis(0)->GetNbins();

    histos.add<TH2>("QA/nRecPerEvAntid", ";Centrality (%);#it{N}_{#bar{d}};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});
    histos.add<TH2>("QA/nRecPerEvAntip", ";Centrality (%);#it{N}_{#bar{p}};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});
    histos.add<TH2>("QA/nRecPerEvAntiL", ";Centrality (%);#it{N}_{#bar{#Lambda}};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});
    histos.add<TH2>("QA/nRecPerEvL", ";Centrality (%);#it{N}_{#Lambda};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});

    // event QA
    if (doprocessRun3) {
      histos.add<TH2>("QA/PvMultVsCent", ";Centrality T0C (%);#it{N}_{PV contributors};", HistType::kTH2F, {centAxis, multAxis});
      histos.add<TH2>("QA/MultVsCent", ";Centrality T0C (%);Multiplicity T0C;", HistType::kTH2F, {centAxis, multFt0Axis});
    } else if (doprocessRun2) {
      histos.add<TH2>("QA/V0MvsCL0", ";Centrality CL0 (%);Centrality V0M (%)", HistType::kTH2F, {centAxis, centAxis});
      histos.add<TH2>("QA/trackletsVsV0M", ";Centrality CL0 (%);Centrality V0M (%)", HistType::kTH2F, {centAxis, multAxis});
    }

    // v0 QA
    histos.add<TH3>("QA/massLambda", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massLambdaAxis});
    histos.add<TH1>("QA/cosPa", ";cosPa;Entries", HistType::kTH1F, {cosPaAxis});
    histos.add<TH1>("QA/cosPaSig", ";cosPa;Entries", HistType::kTH1F, {cosPaAxis});
    histos.add<TH1>("QA/cosPaBkg", ";cosPa;Entries", HistType::kTH1F, {cosPaAxis});
    histos.add<TH1>("QA/dcaV0daughSig", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0daughBkg", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0PvSig", ";dcaV0Pv;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0PvBkg", ";dcaV0Pv;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH2>("QA/cosPaDcaV0daughSig", ";cosPa;dcaV0daugh", HistType::kTH2F, {cosPaAxis, dcaV0daughAxis});
    histos.add<TH2>("QA/cosPaDcaV0daughBkg", ";cosPa;dcaV0daugh", HistType::kTH2F, {cosPaAxis, dcaV0daughAxis});
    histos.add<TH3>("QA/massLambdaEvRej", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massLambdaAxis});
    histos.add<TH3>("QA/massLambdaEvRejSig", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massLambdaAxis});
    histos.add<TH3>("QA/massLambdaEvRejBkg", ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, massLambdaAxis});
    histos.add<TH1>("QA/radius", ";radius;Entries", HistType::kTH1F, {radiusAxis});
    histos.add<TH1>("QA/dcaV0daugh", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0Pv", ";dcaV0Pv;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaPosPv", ";dcaPosPv;Entries", HistType::kTH1F, {dcaDaughPvAxis});
    histos.add<TH1>("QA/dcaNegPv", ";dcaNegPv;Entries", HistType::kTH1F, {dcaDaughPvAxis});
    histos.add<TH1>("QA/cosPaBeforeCut", ";cosPa;Entries", HistType::kTH1F, {cosPaAxis});
    histos.add<TH1>("QA/radiusBeforeCut", ";radius;Entries", HistType::kTH1F, {radiusAxis});
    histos.add<TH1>("QA/dcaV0daughBeforeCut", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0PvBeforeCut", ";dcaV0Pv;Entries", HistType::kTH1F, {dcaV0daughAxis});

    // d QA
    histos.add<TH2>("QA/dcaPv", ";#it{p}_{T} (GeV/#it{c});dcaPv;Entries", HistType::kTH2F, {momAxis, dcaDaughPvAxis});
    histos.add<TH1>("QA/nClsTPC", ";tpcCls;Entries", HistType::kTH1F, {tpcClsAxis});
    histos.add<TH1>("QA/nCrossedRowsTPC", ";nCrossedRowsTPC;Entries", HistType::kTH1F, {tpcClsAxis});
    histos.add<TH2>("QA/dcaPvBefore", ";#it{p}_{T} (GeV/#it{c});dcaPv;Entries", HistType::kTH2F, {momAxis, dcaDaughPvAxis});
    histos.add<TH1>("QA/nClsTPCBeforeCut", ";tpcCls;Entries", HistType::kTH1F, {tpcClsAxis});
    histos.add<TH1>("QA/nCrossedRowsTPCBeforeCut", ";nCrossedRowsTPC;Entries", HistType::kTH1F, {tpcClsAxis});

    // antid and antip QA
    histos.add<TH2>("QA/tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", HistType::kTH2F, {momAxisFine, tpcAxis});
    histos.add<TH2>("QA/tpcSignal_glo", ";#it{p}_{glo} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.);", HistType::kTH2F, {momAxisFine, tpcAxis});

    tpcNsigma[0] = histos.add<TH2>("QA/tpcNsigma_p", ";#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcNsigmaAxis});
    tpcNsigmaGlo[0] = histos.add<TH3>("QA/tpcNsigmaGlo_p", ";Centrality (%);#it{p}_{T} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH3F, {centAxis, momAxis, tpcNsigmaAxis});
    tofMass[0] = histos.add<TH3>("QA/tofMass_p", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});
    tofSignal[0] = histos.add<TH2>("QA/tofSignal_p", ";#it{p}_{TPC} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {momAxisFine, tofAxis});
    tofSignal_glo[0] = histos.add<TH2>("QA/tofSignal_glo_p", ";#it{p}_{T} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {momAxisFine, tofAxis});

    tpcNsigma[1] = histos.add<TH2>("QA/tpcNsigma_d", ";#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcNsigmaAxis});
    tpcNsigmaGlo[1] = histos.add<TH3>("QA/tpcNsigmaGlo_d", ";Centrality (%);#it{p}_{T} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH3F, {centAxis, momAxis, tpcNsigmaAxis});
    tofMass[1] = histos.add<TH3>("QA/tofMass_d", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});
    tofSignal[1] = histos.add<TH2>("QA/tofSignal_d", ";#it{p}_{TPC} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {momAxisFine, tofAxis});
    tofSignal_glo[1] = histos.add<TH2>("QA/tofSignal_glo_d", ";#it{p}_{T} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {momAxisFine, tofAxis});

    // mc histograms
    if (doprocessMcRun3 || doprocessMcRun2) {
      histos.add<TH3>("recL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      histos.add<TH3>("recAntiL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      recTracks[0] = histos.add<TH3>("recP", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntipAxis, deltaEtaAxis});
      recTracks[1] = histos.add<TH3>("recD", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
      recAntiTracks[0] = histos.add<TH3>("recAntip", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntipAxis, deltaEtaAxis});
      recAntiTracks[1] = histos.add<TH3>("recAntid", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
      histos.add<TH3>("genL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      histos.add<TH3>("genAntiL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      genTracks[0] = histos.add<TH3>("genP", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntipAxis, deltaEtaAxis});
      genTracks[1] = histos.add<TH3>("genD", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
      genAntiTracks[0] = histos.add<TH3>("genAntip", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntipAxis, deltaEtaAxis});
      genAntiTracks[1] = histos.add<TH3>("genAntid", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
    }

    // temporary histograms
    tempTracks[0] = tempHistos.add<TH2>("tempAntip", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptAntipAxis});
    tempTracks[1] = tempHistos.add<TH2>("tempAntid", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptAntidAxis});
    tempLambda = tempHistos.add<TH2>("tempLambda", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptLambdaAxis});
    tempAntiLambda = tempHistos.add<TH2>("tempAntiLambda", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptLambdaAxis});

    ptMin = std::array<float, kNpart>{antipPtMin, antidPtMin};
    ptMax = std::array<float, kNpart>{antipPtMax, antidPtMax};
    ptTof = std::array<float, kNpart>{antipPtTof, antidPtTof};

    nSigmaTpcCutLow = std::array<float, kNpart>{antipNsigmaTpcCutLow, antidNsigmaTpcCutLow};
    nSigmaTpcCutUp = std::array<float, kNpart>{antipNsigmaTpcCutUp, antidNsigmaTpcCutUp};
    tpcInnerParamMax = std::array<float, kNpart>{antipTpcInnerParamMax, antidTpcInnerParamMax};
    tofMassMax = std::array<float, kNpart>{antipTofMassMax, antidTofMassMax};
  }

  template <class C, class T>
  int fillRecoEvent(C const& collision, T const& tracksAll, aod::V0s const& V0s, float const& centrality)
  {
    auto tracks = (doprocessRun3 || doprocessMcRun3) ? tracksAll.sliceBy(perCollisionTracksFullIU, collision.globalIndex()) : tracksAll.sliceBy(perCollisionTracksFull, collision.globalIndex());
    candidateTracks[0].clear();
    candidateTracks[1].clear();
    candidateV0s.clear();

    auto rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
    auto subsample = static_cast<int>(rnd * nSubsamples);
    LOGF(info, "%d", subsample);

    gpu::gpustd::array<float, 2> dcaInfo;
    for (const auto& track : tracks) {

      histos.fill(HIST("QA/nClsTPCBeforeCut"), track.tpcNClsFound());
      histos.fill(HIST("QA/nCrossedRowsTPCBeforeCut"), track.tpcNClsCrossedRows());

      if (!selectTrack(track)) {
        continue;
      }

      auto trackParCov = getTrackParCov(track);
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCov, 2.f, fitter.getMatCorrType(), &dcaInfo);
      auto dca = std::hypot(dcaInfo[0], dcaInfo[1]);
      auto trackPt = trackParCov.getPt();
      auto trackEta = trackParCov.getEta();
      histos.fill(HIST("QA/dcaPvBefore"), trackPt, dca);
      if (dca > trackDcaCut) {
        continue;
      }
      histos.fill(HIST("QA/dcaPv"), trackPt, dca);

      histos.fill(HIST("QA/nClsTPC"), track.tpcNClsFound());
      histos.fill(HIST("QA/nCrossedRowsTPC"), track.tpcNClsCrossedRows());
      histos.fill(HIST("QA/tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
      histos.fill(HIST("QA/tpcSignal_glo"), track.p(), track.tpcSignal());

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

        float beta{track.hasTOF() ? track.length() / (track.tofSignal() - track.tofEvTime()) * o2::pid::tof::kCSPEDDInv : -999.f};
        beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta));
        float mass{track.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f)};
        bool hasTof = track.hasTOF() && track.tofChi2() < 3;

        if (trackPt <= ptTof[iP] || (trackPt > ptTof[iP] && hasTof && std::abs(mass - partMass[iP]) < tofMassMaxQA)) { // for QA histograms
          tpcNsigmaGlo[iP]->Fill(centrality, trackPt, nSigmaTPC);
          if (nSigmaTPC > nSigmaTpcCutLow[iP] && nSigmaTPC < nSigmaTpcCutUp[iP]) {
            tofMass[iP]->Fill(centrality, trackPt, mass);
          }
        }

        if (nSigmaTPC < nSigmaTpcCutLow[iP] || nSigmaTPC > nSigmaTpcCutUp[iP]) {
          continue;
        }

        tpcNsigma[iP]->Fill(track.tpcInnerParam(), nSigmaTPC);
        if (trackPt > ptTof[iP] && hasTof) {
          tofSignal_glo[iP]->Fill(track.p(), beta);
          tofSignal[iP]->Fill(track.tpcInnerParam(), beta);
        }

        // temporary cut to reject fake matches (run 3)
        if (track.tpcInnerParam() < tpcInnerParamMax[iP]) {
          continue;
        }
        if (trackPt > ptTof[iP] && !hasTof) {
          continue;
        }

        if (trackPt <= ptTof[iP] || (trackPt > ptTof[iP] && hasTof && std::abs(mass - partMass[iP]) < tofMassMax[iP])) {
          tempTracks[iP]->Fill(std::abs(trackEta), trackPt);
          CandidateTrack candTrack;
          candTrack.pt = trackPt;
          candTrack.eta = trackEta;
          candTrack.mass = partMass[iP];
          candTrack.matter = track.sign() > 0.;
          candTrack.dcapv = dca;
          candTrack.tpcncls = track.tpcNClsFound();
          candTrack.tpcnsigma = nSigmaTPC;
          candTrack.tofmass = hasTof ? mass : -999.f;
          candTrack.globalIndex = track.globalIndex();
          candidateTracks[iP].push_back(candTrack);
        }
      }
    }

    std::vector<int64_t> trkId;
    for (const auto& v0 : V0s) {
      auto posTrack = v0.posTrack_as<T>();
      auto negTrack = v0.negTrack_as<T>();

      bool posSelect = selectV0Daughter(posTrack);
      bool negSelect = selectV0Daughter(negTrack);
      if (!posSelect || !negSelect)
        continue;

      if (doprocessRun2 || doprocessMcRun2) {
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

      if (std::abs(nSigmaTPCPos) > v0setting_nsigmatpc || std::abs(nSigmaTPCNeg) > v0setting_nsigmatpc) {
        continue;
      }

      // veto on K0s mass
      if (std::abs(mK0Short - o2::constants::physics::MassK0Short) < vetoMassK0Short) {
        continue;
      }

      float dcaV0dau = std::sqrt(fitter.getChi2AtPCACandidate());
      histos.fill(HIST("QA/dcaV0daughBeforeCut"), dcaV0dau);
      if (dcaV0dau > v0setting_dcav0dau) {
        continue;
      }

      std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
      const auto& vtx = fitter.getPCACandidate();

      float radiusV0 = std::hypot(vtx[0], vtx[1]);
      histos.fill(HIST("QA/radiusBeforeCut"), radiusV0);
      if (radiusV0 < v0setting_radius || radiusV0 > v0radiusMax) {
        continue;
      }

      float dcaV0Pv = CalculateDCAStraightToPV(
        vtx[0], vtx[1], vtx[2],
        momPos[0] + momNeg[0],
        momPos[1] + momNeg[1],
        momPos[2] + momNeg[2],
        collision.posX(), collision.posY(), collision.posZ());
      histos.fill(HIST("QA/dcaV0PvBeforeCut"), dcaV0Pv);
      if (std::abs(dcaV0Pv) > v0setting_dcav0pv) {
        continue;
      }

      double cosPA = RecoDecay::cpa(primVtx, vtx, momV0);
      histos.fill(HIST("QA/cosPaBeforeCut"), cosPA);
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

      if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCutQA) { // for QA histograms
        continue;
      }
      histos.fill(HIST("QA/massLambda"), centrality, ptV0, mLambda);

      // if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) {
      //   continue;
      // }
      histos.fill(HIST("QA/cosPa"), cosPA);
      histos.fill(HIST("QA/radius"), radiusV0);
      histos.fill(HIST("QA/dcaV0daugh"), dcaV0dau);
      histos.fill(HIST("QA/dcaPosPv"), posDcaToPv);
      histos.fill(HIST("QA/dcaNegPv"), negDcaToPv);
      histos.fill(HIST("QA/dcaV0Pv"), dcaV0Pv);

      if (matter) {
        tempHistos.fill(HIST("tempLambda"), std::abs(etaV0), ptV0);
      } else {
        tempHistos.fill(HIST("tempAntiLambda"), std::abs(etaV0), ptV0);
      }

      trkId.emplace_back(posTrack.globalIndex());
      trkId.emplace_back(negTrack.globalIndex());

      CandidateV0 candV0;
      candV0.pt = ptV0;
      candV0.eta = etaV0;
      candV0.mass = mLambda;
      candV0.cpa = cosPA;
      candV0.dcav0daugh = dcaV0dau;
      candV0.dcav0pv = dcaV0Pv;
      candV0.dcanegpv = negDcaToPv;
      candV0.dcapospv = posDcaToPv;
      candV0.tpcnsigmaneg = nSigmaTPCNeg;
      candV0.tpcnsigmapos = nSigmaTPCPos;
      candV0.matter = matter;
      candV0.globalIndexPos = posTrack.globalIndex();
      candV0.globalIndexNeg = negTrack.globalIndex();
      candidateV0s.push_back(candV0);
    }

    for (auto& candidateV0 : candidateV0s) {
      histos.fill(HIST("QA/massLambdaEvRej"), centrality, candidateV0.pt, candidateV0.mass);
    }

    histos.fill(HIST("nEv"), subsample, centrality);

    if ((doprocessMcRun3 || doprocessMcRun2) && fillOnlySignal)
      return subsample;

    histos.fill(HIST("QA/nRecPerEvAntip"), centrality, tempTracks[0]->GetEntries());
    histos.fill(HIST("QA/nRecPerEvAntid"), centrality, tempTracks[1]->GetEntries());
    histos.fill(HIST("QA/nRecPerEvAntiL"), centrality, tempAntiLambda->GetEntries());
    histos.fill(HIST("QA/nRecPerEvL"), centrality, tempLambda->GetEntries());

    return 0;
  }

  template <class C, class T>
  void fillMcEvent(C const& collision, T const& tracks, aod::V0s const& V0s, float const& centrality, aod::McParticles const&, aod::McTrackLabels const& mcLabels)
  {
    /* int subsample =  */fillRecoEvent<C, T>(collision, tracks, V0s, centrality);
    // subsample += 1;
    // if (candidateV0s.size() == 1 && candidateV0s[0].pt < -998.f && candidateV0s[0].eta < -998.f && candidateV0s[0].globalIndexPos == -999 && candidateV0s[0].globalIndexPos == -999) {
    //   return;
    // }

    for (int iP{0}; iP < kNpart; ++iP) {
      for (auto& candidateTrack : candidateTracks[iP]) {
        candidateTrack.isreco = true;

        auto mcLab = mcLabels.rawIteratorAt(candidateTrack.globalIndex);
        if (mcLab.has_mcParticle()) {
          auto mcTrack = mcLab.template mcParticle_as<aod::McParticles>();
          if (std::abs(mcTrack.pdgCode()) != partPdg[iP])
            continue;
          if (((mcTrack.flags() & 0x8) && doprocessMcRun2) || (mcTrack.flags() & 0x2) || (mcTrack.flags() & 0x1))
            continue;
          if (!mcTrack.isPhysicalPrimary())
            continue;

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
                histos.fill(HIST("QA/cosPaBkg"), candidateV0.cpa);
                histos.fill(HIST("QA/dcaV0daughBkg"), candidateV0.dcav0daugh);
                histos.fill(HIST("QA/dcaV0PvBkg"), candidateV0.dcav0pv);
                histos.fill(HIST("QA/cosPaDcaV0daughBkg"), candidateV0.cpa, candidateV0.dcav0daugh);
                histos.fill(HIST("QA/massLambdaEvRejBkg"), centrality, candidateV0.pt, candidateV0.mass);
                continue;
              }
              if (!posMother.isPhysicalPrimary() && !posMother.has_mothers())
                continue;
              if (((posMother.flags() & 0x8) && doprocessMcRun2) || (posMother.flags() & 0x2) || (posMother.flags() & 0x1))
                continue;
              histos.fill(HIST("QA/cosPaSig"), candidateV0.cpa);
              histos.fill(HIST("QA/dcaV0daughSig"), candidateV0.dcav0daugh);
              histos.fill(HIST("QA/dcaV0PvSig"), candidateV0.dcav0pv);
              histos.fill(HIST("QA/cosPaDcaV0daughSig"), candidateV0.cpa, candidateV0.dcav0daugh);
              histos.fill(HIST("QA/massLambdaEvRejSig"), centrality, candidateV0.pt, candidateV0.mass);

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

    if (fillOnlySignal) {
      histos.fill(HIST("QA/nRecPerEvAntip"), centrality, tempTracks[0]->GetEntries());
      histos.fill(HIST("QA/nRecPerEvAntid"), centrality, tempTracks[1]->GetEntries());
      histos.fill(HIST("QA/nRecPerEvAntiL"), centrality, tempAntiLambda->GetEntries());
      histos.fill(HIST("QA/nRecPerEvL"), centrality, tempLambda->GetEntries());
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
      if (((mcPart.flags() & 0x8) && doprocessMcRun2) || (mcPart.flags() & 0x2) || (mcPart.flags() & 0x1))
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
        candV0.pdgcode = mcPart.pdgCode();
        auto it = find_if(candidateV0s.begin(), candidateV0s.end(), [&](CandidateV0 v0){return v0.mcIndex == mcPart.globalIndex();});
        if (it == candidateV0s.end()) {
          candidateV0s.emplace_back(candV0);
        }
      } else if (std::abs(pdgCode) == partPdg[0] || std::abs(pdgCode) == partPdg[1]) {
        int iP = 1;
        if (std::abs(pdgCode) == partPdg[0]) {
          iP = 0;
        }
        if (!mcPart.isPhysicalPrimary() && !mcPart.has_mothers())
          continue;
        auto genPt = std::hypot(mcPart.px(), mcPart.py());
        CandidateTrack candTrack;
        candTrack.genpt = genPt;
        candTrack.geneta = mcPart.eta();
        candTrack.pdgcode = mcPart.pdgCode();
        auto it = find_if(candidateTracks[iP].begin(), candidateTracks[iP].end(), [&](CandidateTrack trk){return trk.mcIndex == mcPart.globalIndex();});
        if (it == candidateTracks[iP].end()) {
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

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        continue;

      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;

      if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;

      if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
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
    }
  }
  PROCESS_SWITCH(ebyeMaker, processRun3, "process (Run 3)", false);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::TrackletMults> const& collisions, TracksFull const& tracks, aod::V0s const& V0s, BCsWithRun2Info const&)
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

      auto centrality = collision.centRun2V0M();
      if (!(collision.sel7() && collision.alias_bit(kINT7)) && (!kINT7Intervals || (kINT7Intervals && ((centrality >= 10 && centrality < 30) || centrality > 50))))
        continue;

      auto centralityCl0 = collision.centRun2CL0();
      if (kUseEstimatorsCorrelationCut) {
        const auto& x = centralityCl0;
        const double center = estimatorsCorrelationCoef[0] + estimatorsCorrelationCoef[1] * x;
        const double sigma = estimatorsSigmaPars[0] + estimatorsSigmaPars[1] * x + estimatorsSigmaPars[2] * std::pow(x, 2) + estimatorsSigmaPars[3] * std::pow(x, 3);
        if (centrality < center - deltaEstimatorNsigma[0] * sigma || centrality > center + deltaEstimatorNsigma[1] * sigma) {
          continue;
        }
      }

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      auto multTracklets = collision.multTracklets();
      fillRecoEvent(collision, tracks, V0Table_thisCollision, centrality);

      histos.fill(HIST("QA/V0MvsCL0"), centralityCl0, centrality);
      histos.fill(HIST("QA/trackletsVsV0M"), centrality, multTracklets);

      for (auto& candidateV0 : candidateV0s) {
        // LOG(info) << candidateV0.pt;
        lambdaEbyeTable(
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.matter,
          candidateV0.dcav0pv,
          candidateV0.dcanegpv,
          candidateV0.dcapospv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          candidateV0.tpcnsigmaneg,
          candidateV0.tpcnsigmapos,
          candidateV0.globalIndexNeg,
          candidateV0.globalIndexPos);
      }

      for (auto& candidateTrack : candidateTracks[0]) { // deuterons
        // LOG(info) << candidateV0.pt;
        nucleiEbyeTable(
          candidateTrack.pt,
          candidateTrack.eta,
          candidateTrack.mass,
          candidateTrack.matter,
          candidateTrack.dcapv,
          candidateTrack.tpcncls,
          candidateTrack.tpcnsigma,
          candidateTrack.tofmass);
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processRun2, "process (Run 2)", false);

  void processMcRun3(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::McCollisions const& mcCollisions, TracksFullIU const& tracks, aod::V0s const& V0s, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, aod::BCsWithTimestamps const&)
  {
    std::vector<std::pair<bool, float>> goodCollisions(mcCollisions.size(), std::make_pair(false, -999.));
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.sel8())
        continue;

      if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        continue;

      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;

      if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;

      if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      auto centrality = collision.centFT0C();
      goodCollisions[collision.mcCollisionId()].first = true;
      goodCollisions[collision.mcCollisionId()].second = centrality;

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillMcEvent(collision, tracks, V0Table_thisCollision, centrality, mcParticles, mcLab);
      if (candidateV0s.size() == 1 && candidateV0s[0].pt < -998.f && candidateV0s[0].eta < -998.f && candidateV0s[0].globalIndexPos == -999 && candidateV0s[0].globalIndexPos == -999) {
        goodCollisions[collision.mcCollisionId()].first = false;
      }
    }
  }
  PROCESS_SWITCH(ebyeMaker, processMcRun3, "process MC (Run 3)", false);

  void processMcRun2(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentRun2V0Ms> const& collisions, aod::McCollisions const& mcCollisions, TracksFull const& tracks, aod::V0s const& V0s, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, BCsWithRun2Info const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      auto centrality = collision.centRun2V0M();

      histos.fill(HIST("QA/zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      fillMcEvent(collision, tracks, V0Table_thisCollision, centrality, mcParticles, mcLab);
      fillMcGen(mcParticles, mcLab, collision.mcCollisionId());

      for (auto& candidateV0 : candidateV0s) {
        // LOG(info) << candidateV0.pt;
        mcLambdaEbyeTable(
          candidateV0.pt,
          candidateV0.eta,
          candidateV0.mass,
          candidateV0.matter,
          candidateV0.dcav0pv,
          candidateV0.dcanegpv,
          candidateV0.dcapospv,
          candidateV0.dcav0daugh,
          candidateV0.cpa,
          candidateV0.tpcnsigmaneg,
          candidateV0.tpcnsigmapos,
          candidateV0.globalIndexNeg,
          candidateV0.globalIndexPos,
          candidateV0.genpt,
          candidateV0.geneta,
          candidateV0.pdgcode,
          candidateV0.isreco);
      }

      for (auto& candidateTrack : candidateTracks[0]) { // deuterons
        // LOG(info) << candidateV0.pt;
        mcNucleiEbyeTable(
          candidateTrack.pt,
          candidateTrack.eta,
          candidateTrack.mass,
          candidateTrack.matter,
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
  PROCESS_SWITCH(ebyeMaker, processMcRun2, "process MC (Run 2)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ebyeMaker>(cfgc)};
}
