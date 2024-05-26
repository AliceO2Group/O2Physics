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

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr int kNpart = 2;
constexpr double partMass[kNpart]{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron};
constexpr double partPdg[kNpart]{2212, o2::constants::physics::kDeuteron};
std::shared_ptr<THnSparse> nAntid;
std::shared_ptr<THnSparse> nAntip;
std::shared_ptr<THnSparse> nAntiL;
std::shared_ptr<THnSparse> nL;
std::shared_ptr<THnSparse> nSqAntid;
std::shared_ptr<THnSparse> nSqAntip;
std::shared_ptr<THnSparse> nSqAntiL;
std::shared_ptr<THnSparse> nSqL;
std::shared_ptr<THnSparse> nAntipAntid;
std::shared_ptr<THnSparse> nLantiL;
std::shared_ptr<THnSparse> nLantid;
std::shared_ptr<THnSparse> nAntiLantid;
std::shared_ptr<THnSparse> nGenAntid;
std::shared_ptr<THnSparse> nGenAntip;
std::shared_ptr<THnSparse> nGenAntiL;
std::shared_ptr<THnSparse> nGenL;
std::shared_ptr<THnSparse> nGenSqAntid;
std::shared_ptr<THnSparse> nGenSqAntip;
std::shared_ptr<THnSparse> nGenSqAntiL;
std::shared_ptr<THnSparse> nGenSqL;
std::shared_ptr<THnSparse> nGenAntipAntid;
std::shared_ptr<THnSparse> nGenLantiL;
std::shared_ptr<THnSparse> nGenLantid;
std::shared_ptr<THnSparse> nGenAntiLantid;
std::shared_ptr<TH2> tempAntiLambda;
std::shared_ptr<TH2> tempLambda;
std::array<std::shared_ptr<TH3>, kNpart> recTracks;
std::array<std::shared_ptr<TH3>, kNpart> recAntiTracks;
std::array<std::shared_ptr<TH3>, kNpart> genTracks;
std::array<std::shared_ptr<TH3>, kNpart> genAntiTracks;
std::array<std::shared_ptr<TH2>, kNpart> tpcNsigma;
std::array<std::shared_ptr<TH3>, kNpart> tofMass;
} // namespace

struct nucleiEbye {
  std::mt19937 gen32;
  int nSubsamples;

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis subsampleAxis{"subsampleAxis", {30, 0, 30}, "binning of the subsample axis"};
  ConfigurableAxis deltaEtaAxis{"deltaEtaAxis", {4, 0, 0.8}, "binning of the delta eta axis"};
  ConfigurableAxis ptAntidAxis{"ptAntidAxis", {VARIABLE_WIDTH, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f}, "binning of the antideuteron pT axis (GeV/c)"};
  ConfigurableAxis ptAntipAxis{"ptAntipAxis", {VARIABLE_WIDTH, 0.4f, 0.6f, 0.7f, 0.8f, 0.9f}, "binning of the antiproton pT axis (GeV/c)"};
  ConfigurableAxis ptLambdaAxis{"ptLambdaAxis", {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f}, "binning of the (anti)lambda pT axis (GeV/c)"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
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

  Configurable<float> antidPtMin{"antidPtMin", 0.8f, "minimum antideuteron pT (GeV/c)"};
  Configurable<float> antidPtTof{"antidPtTof", 1.0f, "antideuteron pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antidPtMax{"antidPtMax", 1.8f, "maximum antideuteron pT (GeV/c)"};

  Configurable<float> antipPtMin{"antipPtMin", 0.4f, "minimum antiproton pT (GeV/c)"};
  Configurable<float> antipPtTof{"antipPtTof", 0.6f, "antiproton pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antipPtMax{"antipPtMax", 0.9f, "maximum antiproton pT (GeV/c)"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 3.0f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> trackNclusTpcCut{"trackNclusTPCcut", 70, "Minimum number of TPC clusters"};
  Configurable<float> trackDcaCut{"trackDcaCut", 0.1f, "DCA antid to PV"};

  Configurable<float> antidNsigmaTpcCutLow{"antidNsigmaTpcCutLow", 4.f, "TPC PID cut low"};
  Configurable<float> antidNsigmaTpcCutUp{"antidNsigmaTpcCutUp", 4.f, "TPC PID cut up"};
  Configurable<float> antidTofMassMax{"tofMassMax", 0.3f, "(temporary) tof mass cut"};

  Configurable<float> antipNsigmaTpcCutLow{"antipNsigmaTpcCutLow", 4.f, "TPC PID cut low"};
  Configurable<float> antipNsigmaTpcCutUp{"antipNsigmaTpcCutUp", 4.f, "TPC PID cut up"};
  Configurable<float> antipTofMassMax{"antipTofMassMax", 0.3f, "(temporary) tof mass cut"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcav0pv{"v0setting_dcav0pv", 1, "DCA V0 to Pv"};
  Configurable<float> v0setting_dcadaughtopv{"v0setting_dcadaughtopv", 0.1f, "DCA Pos To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_nsigmatpc{"v0setting_nsigmatpc", 4.f, "nsigmatpc"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.005f, "maximum deviation from PDG mass"};
  Configurable<float> lambdaMassCutQA{"lambdaMassCutQA", 0.02f, "maximum deviation from PDG mass (for QA histograms)"};

  std::array<float, kNpart> ptMin;
  std::array<float, kNpart> ptTof;
  std::array<float, kNpart> ptMax;
  std::array<float, kNpart> nSigmaTpcCutLow;
  std::array<float, kNpart> nSigmaTpcCutUp;
  std::array<float, kNpart> tofMassMax;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry tempHistos{"tempHistos", {}, OutputObjHandlingPolicy::TransientObject};

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

  void init(o2::framework::InitContext&)
  {

    uint32_t randomSeed = static_cast<uint32_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    gen32.seed(randomSeed);

    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    auto hNev = histos.add<THnSparse>("nEv", ";Subsample;Centrality (%);", HistType::kTHnSparseD, {subsampleAxis, centAxis});
    nSubsamples = hNev->GetAxis(0)->GetNbins();

    histos.add<TH2>("QA/nRecPerEvAntid", ";Centrality (%);#it{N}_{#bar{d}};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});
    histos.add<TH2>("QA/nRecPerEvAntip", ";Centrality (%);#it{N}_{#bar{p}};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});
    histos.add<TH2>("QA/nRecPerEvAntiL", ";Centrality (%);#it{N}_{#bar{#Lambda}};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});
    histos.add<TH2>("QA/nRecPerEvL", ";Centrality (%);#it{N}_{#Lambda};#it{N}_{ev}", HistType::kTH2D, {centAxis, nGenRecAxis});

    nAntid = histos.add<THnSparse>("nAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis});
    nAntip = histos.add<THnSparse>("nAntip", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{p}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntipAxis});
    nAntiL = histos.add<THnSparse>("nAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});
    nL = histos.add<THnSparse>("nL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});

    nSqAntid = histos.add<THnSparse>("nSqAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis, ptAntidAxis});
    nSqAntip = histos.add<THnSparse>("nSqAntip", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{p}) (GeV/#it{c});#it{p}_{T}(#bar{p}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntipAxis, ptAntipAxis});
    nSqAntiL = histos.add<THnSparse>("nSqAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    nSqL = histos.add<THnSparse>("nSqL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#Lambda) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});

    nAntipAntid = histos.add<THnSparse>("nAntipAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{p}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntipAxis, ptAntidAxis});
    nLantiL = histos.add<THnSparse>("nLantiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    nLantid = histos.add<THnSparse>("nLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});
    nAntiLantid = histos.add<THnSparse>("nAntiLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});

    // mc generated
    nGenAntid = histos.add<THnSparse>("nGenAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis});
    nGenAntip = histos.add<THnSparse>("nGenAntip", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{p}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntipAxis});
    nGenAntiL = histos.add<THnSparse>("nGenAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});
    nGenL = histos.add<THnSparse>("nGenL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});

    nGenSqAntid = histos.add<THnSparse>("nGenSqAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis, ptAntidAxis});
    nGenSqAntip = histos.add<THnSparse>("nGenSqAntip", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{p}) (GeV/#it{c});#it{p}_{T}(#bar{p}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntipAxis, ptAntipAxis});
    nGenSqAntiL = histos.add<THnSparse>("nGenSqAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    nGenSqL = histos.add<THnSparse>("nGenSqL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#Lambda) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});

    nGenAntipAntid = histos.add<THnSparse>("nGenAntipAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{p}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntipAxis, ptAntidAxis});
    nGenLantiL = histos.add<THnSparse>("nGenLantiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    nGenLantid = histos.add<THnSparse>("nGenLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});
    nGenAntiLantid = histos.add<THnSparse>("nGenAntiLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});

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
    histos.add<TH1>("QA/dcaV0daugh", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0Pv", ";dcaV0Pv;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaPosPv", ";dcaPosPv;Entries", HistType::kTH1F, {dcaDaughPvAxis});
    histos.add<TH1>("QA/dcaNegPv", ";dcaNegPv;Entries", HistType::kTH1F, {dcaDaughPvAxis});
    histos.add<TH1>("QA/cosPaBeforeCut", ";cosPa;Entries", HistType::kTH1F, {cosPaAxis});
    histos.add<TH1>("QA/dcaV0daughBeforeCut", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("QA/dcaV0PvBeforeCut", ";dcaV0Pv;Entries", HistType::kTH1F, {dcaV0daughAxis});

    // d QA
    histos.add<TH2>("QA/dcaPv", ";#it{p}_{T} (GeV/#it{c});dcaPv;Entries", HistType::kTH2F, {momAxis, dcaDaughPvAxis});
    histos.add<TH1>("QA/nClsTPC", ";tpcCls;Entries", HistType::kTH1F, {tpcClsAxis});
    histos.add<TH2>("QA/dcaPvBefore", ";#it{p}_{T} (GeV/#it{c});dcaPv;Entries", HistType::kTH2F, {momAxis, dcaDaughPvAxis});
    histos.add<TH1>("QA/nClsTPCBeforeCut", ";tpcCls;Entries", HistType::kTH1F, {tpcClsAxis});

    tpcNsigma[0] = histos.add<TH2>("QA/tpcNsigma_p", ";#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcNsigmaAxis});
    tofMass[0] = histos.add<TH3>("QA/tofMass_p", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});

    tpcNsigma[1] = histos.add<TH2>("QA/tpcNsigma_d", ";#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcNsigmaAxis});
    tofMass[1] = histos.add<TH3>("QA/tofMass_d", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});

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
    tofMassMax = std::array<float, kNpart>{antipTofMassMax, antidTofMassMax};
  }

  template <class C, class T>
  int fillRecoEvent(C const& collision, T const& tracksAll, aod::V0s const& V0s, float const& centrality)
  {
    auto tracks = (doprocessRun3 || doprocessMcRun3) ? tracksAll.sliceBy(perCollisionTracksFullIU, collision.globalIndex()) : tracksAll.sliceBy(perCollisionTracksFull, collision.globalIndex());
    candidateTracks[0].clear();
    candidateTracks[1].clear();
    candidateV0s.clear();

    tempTracks[0]->Reset();
    tempTracks[1]->Reset();
    tempLambda->Reset();
    tempAntiLambda->Reset();
    auto rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
    auto subsample = static_cast<int>(rnd * nSubsamples);

    gpu::gpustd::array<float, 2> dcaInfo;
    for (const auto& track : tracks) {

      histos.fill(HIST("QA/nClsTPCBeforeCut"), track.tpcNClsFound());
      histos.fill(HIST("QA/nCrossedRowsTPCBeforeCut"), track.tpcNClsCrossedRows());

      if (!selectTrack(track)) {
        continue;
      }

      if (track.sign() > 0.) {
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

      if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) {
        continue;
      }
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
      candV0.globalIndexPos = posTrack.globalIndex();
      candV0.globalIndexNeg = negTrack.globalIndex();
      candidateV0s.push_back(candV0);
    }

    // reject events having multiple v0s from same tracks (TODO: also across collisions?)
    std::sort(trkId.begin(), trkId.end());
    if (std::adjacent_find(trkId.begin(), trkId.end()) != trkId.end()) {
      candidateV0s.clear();

      CandidateV0 candV0;
      candV0.pt = -999.f;
      candV0.eta = -999.f;
      candV0.globalIndexPos = -999;
      candV0.globalIndexNeg = -999;
      candidateV0s.push_back(candV0);
      return -1;
    }
    for (auto& candidateV0 : candidateV0s) {
      histos.fill(HIST("QA/massLambdaEvRej"), centrality, candidateV0.pt, candidateV0.mass);
    }

    histos.fill(HIST("nEv"), subsample, centrality);

    if ((doprocessMcRun3 || doprocessMcRun2) && fillOnlySignal)
      return subsample;

    fillHistoN(nAntip, tempTracks[0], subsample, centrality);
    fillHistoN(nAntid, tempTracks[1], subsample, centrality);
    fillHistoN(nAntiL, tempAntiLambda, subsample, centrality);
    fillHistoN(nL, tempLambda, subsample, centrality);

    fillHistoN(nSqAntip, tempTracks[0], tempTracks[0], subsample, centrality);
    fillHistoN(nSqAntid, tempTracks[1], tempTracks[1], subsample, centrality);
    fillHistoN(nSqAntiL, tempAntiLambda, tempAntiLambda, subsample, centrality);
    fillHistoN(nSqL, tempLambda, tempLambda, subsample, centrality);

    fillHistoN(nAntipAntid, tempTracks[0], tempTracks[1], subsample, centrality);
    fillHistoN(nLantid, tempLambda, tempTracks[1], subsample, centrality);
    fillHistoN(nLantiL, tempLambda, tempAntiLambda, subsample, centrality);
    fillHistoN(nAntiLantid, tempAntiLambda, tempTracks[1], subsample, centrality);

    histos.fill(HIST("QA/nRecPerEvAntip"), centrality, tempTracks[0]->GetEntries());
    histos.fill(HIST("QA/nRecPerEvAntid"), centrality, tempTracks[1]->GetEntries());
    histos.fill(HIST("QA/nRecPerEvAntiL"), centrality, tempAntiLambda->GetEntries());
    histos.fill(HIST("QA/nRecPerEvL"), centrality, tempLambda->GetEntries());

    return 0;
  }

  template <class C, class T>
  void fillMcEvent(C const& collision, T const& tracks, aod::V0s const& V0s, float const& centrality, aod::McParticles const&, aod::McTrackLabels const& mcLabels)
  {
    int subsample = fillRecoEvent<C, T>(collision, tracks, V0s, centrality);
    if (candidateV0s.size() == 1 && candidateV0s[0].pt < -998.f && candidateV0s[0].eta < -998.f && candidateV0s[0].globalIndexPos == -999 && candidateV0s[0].globalIndexPos == -999) {
      return;
    }

    if (fillOnlySignal) {
      tempTracks[0]->Reset();
      tempTracks[1]->Reset();
      tempLambda->Reset();
      tempAntiLambda->Reset();
    }

    for (int iP{0}; iP < kNpart; ++iP) {
      for (auto& candidateTrack : candidateTracks[iP]) {
        auto mcLab = mcLabels.rawIteratorAt(candidateTrack.globalIndex);
        if (mcLab.has_mcParticle()) {
          auto mcTrack = mcLab.template mcParticle_as<aod::McParticles>();
          if (std::abs(mcTrack.pdgCode()) != partPdg[iP])
            continue;
          if (((mcTrack.flags() & 0x8) && doprocessMcRun2) || (mcTrack.flags() & 0x2) || (mcTrack.flags() & 0x1))
            continue;
          if (!mcTrack.isPhysicalPrimary())
            continue;
          if (mcTrack.pdgCode() > 0) {
            recTracks[iP]->Fill(centrality, candidateTrack.pt, std::abs(candidateTrack.eta));
          } else {
            recAntiTracks[iP]->Fill(centrality, candidateTrack.pt, std::abs(candidateTrack.eta));
            if (fillOnlySignal)
              tempTracks[iP]->Fill(std::abs(candidateTrack.eta), candidateTrack.pt);
          }
        }
      }
    }
    for (auto& candidateV0 : candidateV0s) {
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
              if (posMother.pdgCode() > 0) {
                histos.fill(HIST("recL"), centrality, candidateV0.pt, std::abs(candidateV0.eta));
                if (fillOnlySignal)
                  tempLambda->Fill(std::abs(candidateV0.eta), candidateV0.pt);
              } else {
                histos.fill(HIST("recAntiL"), centrality, candidateV0.pt, std::abs(candidateV0.eta));
                if (fillOnlySignal)
                  tempAntiLambda->Fill(std::abs(candidateV0.eta), candidateV0.pt);
              }
            }
          }
        }
      }
    }

    if (fillOnlySignal) {
      fillHistoN(nAntip, tempTracks[0], subsample, centrality);
      fillHistoN(nAntid, tempTracks[1], subsample, centrality);
      fillHistoN(nAntiL, tempAntiLambda, subsample, centrality);
      fillHistoN(nL, tempLambda, subsample, centrality);

      fillHistoN(nSqAntip, tempTracks[0], tempTracks[0], subsample, centrality);
      fillHistoN(nSqAntid, tempTracks[1], tempTracks[1], subsample, centrality);
      fillHistoN(nSqAntiL, tempAntiLambda, tempAntiLambda, subsample, centrality);
      fillHistoN(nSqL, tempLambda, tempLambda, subsample, centrality);

      fillHistoN(nAntipAntid, tempTracks[0], tempTracks[1], subsample, centrality);
      fillHistoN(nLantid, tempLambda, tempTracks[1], subsample, centrality);
      fillHistoN(nLantiL, tempLambda, tempAntiLambda, subsample, centrality);
      fillHistoN(nAntiLantid, tempAntiLambda, tempTracks[1], subsample, centrality);

      histos.fill(HIST("QA/nRecPerEvAntip"), centrality, tempTracks[0]->GetEntries());
      histos.fill(HIST("QA/nRecPerEvAntid"), centrality, tempTracks[1]->GetEntries());
      histos.fill(HIST("QA/nRecPerEvAntiL"), centrality, tempAntiLambda->GetEntries());
      histos.fill(HIST("QA/nRecPerEvL"), centrality, tempLambda->GetEntries());
    }
  }

  void fillMcGen(aod::McParticles const& mcParticles, aod::McTrackLabels const& /*mcLab*/, std::vector<std::pair<bool, float>> const& goodCollisions)
  {
    for (uint64_t iC{0}; iC < goodCollisions.size(); ++iC) {
      if (goodCollisions[iC].first == false) {
        continue;
      }

      tempTracks[0]->Reset();
      tempTracks[1]->Reset();
      tempLambda->Reset();
      tempAntiLambda->Reset();

      auto centrality = goodCollisions[iC].second;
      auto rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
      auto subsample = static_cast<int>(rnd * nSubsamples);
      auto mcParticles_thisCollision = mcParticles.sliceBy(perCollisionMcParts, iC);
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

          if (pdgCode > 0) {
            histos.fill(HIST("genL"), centrality, genPt, std::abs(genEta));
            tempHistos.fill(HIST("tempLambda"), std::abs(genEta), genPt);
          } else {
            histos.fill(HIST("genAntiL"), centrality, genPt, std::abs(genEta));
            tempHistos.fill(HIST("tempAntiLambda"), std::abs(genEta), genPt);
          }
        } else if (std::abs(pdgCode) == partPdg[0] || std::abs(pdgCode) == partPdg[1]) {
          int iP = 1;
          if (std::abs(pdgCode) == partPdg[0]) {
            iP = 0;
          }
          if (!mcPart.isPhysicalPrimary() && !mcPart.has_mothers())
            continue;
          auto genPt = std::hypot(mcPart.px(), mcPart.py());
          if (pdgCode > 0) {
            genTracks[iP]->Fill(centrality, genPt, std::abs(genEta));
          } else {
            genAntiTracks[iP]->Fill(centrality, genPt, std::abs(genEta));
            tempTracks[iP]->Fill(std::abs(genEta), genPt);
          }
        }
      }

      fillHistoN(nGenAntip, tempTracks[0], subsample, centrality);
      fillHistoN(nGenAntid, tempTracks[1], subsample, centrality);
      fillHistoN(nGenAntiL, tempAntiLambda, subsample, centrality);
      fillHistoN(nGenL, tempLambda, subsample, centrality);

      fillHistoN(nGenSqAntip, tempTracks[0], tempTracks[0], subsample, centrality);
      fillHistoN(nGenSqAntid, tempTracks[1], tempTracks[1], subsample, centrality);
      fillHistoN(nGenSqAntiL, tempAntiLambda, tempAntiLambda, subsample, centrality);
      fillHistoN(nGenSqL, tempLambda, tempLambda, subsample, centrality);

      fillHistoN(nGenAntipAntid, tempTracks[0], tempTracks[1], subsample, centrality);
      fillHistoN(nGenLantid, tempLambda, tempTracks[1], subsample, centrality);
      fillHistoN(nGenLantiL, tempLambda, tempAntiLambda, subsample, centrality);
      fillHistoN(nGenAntiLantid, tempAntiLambda, tempTracks[1], subsample, centrality);
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
  PROCESS_SWITCH(nucleiEbye, processRun3, "process (Run 3)", false);

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
    }
  }
  PROCESS_SWITCH(nucleiEbye, processRun2, "process (Run 2)", false);

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

    fillMcGen(mcParticles, mcLab, goodCollisions);
  }
  PROCESS_SWITCH(nucleiEbye, processMcRun3, "process MC (Run 3)", false);

  void processMcRun2(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentRun2V0Ms> const& collisions, aod::McCollisions const& mcCollisions, TracksFull const& tracks, aod::V0s const& V0s, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab, BCsWithRun2Info const&)
  {
    std::vector<std::pair<bool, float>> goodCollisions(mcCollisions.size(), std::make_pair(false, -999.));
    for (auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
        continue;

      auto centrality = collision.centRun2V0M();
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

    fillMcGen(mcParticles, mcLab, goodCollisions);
  }
  PROCESS_SWITCH(nucleiEbye, processMcRun2, "process MC (Run 2)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiEbye>(cfgc)};
}
