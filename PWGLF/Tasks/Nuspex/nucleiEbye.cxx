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
#include <iostream>
#include <memory>
#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "PWGLF/DataModel/LFEbyeTables.h"

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
std::array<std::shared_ptr<TH2>, kNpart> tempTracks;
std::array<std::shared_ptr<TH3>, kNpart> recTracks;
std::array<std::shared_ptr<TH3>, kNpart> recAntiTracks;
std::array<std::shared_ptr<TH3>, kNpart> genTracks;
std::array<std::shared_ptr<TH3>, kNpart> genAntiTracks;
std::array<std::shared_ptr<TH3>, kNpart> tpcNsigmaGlo;
std::array<std::shared_ptr<TH3>, kNpart> tofMass;
} // namespace

struct CandidateV0 {
  float pt;
  float eta;
  float mass;
  float cpa;
  float dcav0daugh;
  float dcav0pv;
  int64_t globalIndexPos = -999;
  int64_t globalIndexNeg = -999;
  int64_t globalIndex = -999;
};

struct CandidateTrack {
  float pt;
  float eta;
  int64_t globalIndex = -999;
};

struct nucleiEbye {
  std::mt19937 gen32;
  std::vector<CandidateV0> candidateV0s;
  std::array<std::vector<CandidateTrack>, 2> candidateTracks;
  PresliceUnsorted<aod::McNucleiEbyeTables> perCollTrack = o2::aod::LFEbyeTable::collEbyeTableId;
  PresliceUnsorted<aod::McLambdaEbyeTables> perCollV0s = o2::aod::LFEbyeTable::collEbyeTableId;
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
  Configurable<float> tofMassMaxQA{"tofMassMaxQA", 0.6f, "(temporary) tof mass cut (for QA histograms)"};

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
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (track.tpcNcls() < trackNclusTpcCut) {
      return false;
    }
    return true;
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

    tpcNsigmaGlo[0] = histos.add<TH3>("QA/tpcNsigmaGlo_p", ";Centrality (%);#it{p}_{T} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH3F, {centAxis, momAxis, tpcNsigmaAxis});
    tofMass[0] = histos.add<TH3>("QA/tofMass_p", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});

    tpcNsigmaGlo[1] = histos.add<TH3>("QA/tpcNsigmaGlo_d", ";Centrality (%);#it{p}_{T} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH3F, {centAxis, momAxis, tpcNsigmaAxis});
    tofMass[1] = histos.add<TH3>("QA/tofMass_d", ";Centrality (%);#it{p}_{T} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH3F, {centAxis, momAxis, tofMassAxis});

    // mc histograms
    if (doprocessMc) {
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

  template <class C, class T, class V0>
  int fillRecoEvent(C const& /*collision*/, T const& tracks, V0 const& V0s, float const& centrality)
  {
    candidateTracks[0].clear();
    candidateTracks[1].clear();
    candidateV0s.clear();

    tempTracks[0]->Reset();
    tempTracks[1]->Reset();
    tempLambda->Reset();
    tempAntiLambda->Reset();
    auto rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
    auto subsample = static_cast<int>(rnd * nSubsamples);

    for (const auto& track : tracks) {

      histos.fill(HIST("QA/nClsTPCBeforeCut"), track.tpcNcls());

      if (!selectTrack(track)) {
        continue;
      }

      if (track.pt() > 0.) {
        continue;
      }

      auto dca = track.dcaPV();
      auto trackPt = std::abs(track.pt());
      auto trackEta = track.eta();
      histos.fill(HIST("QA/dcaPvBefore"), trackPt, dca);
      if (dca > trackDcaCut) {
        continue;
      }
      histos.fill(HIST("QA/dcaPv"), trackPt, dca);
      histos.fill(HIST("QA/nClsTPC"), track.tpcNcls());

      for (int iP{0}; iP < kNpart; ++iP) {
        if (track.mass() != iP)
          continue;
        if (trackPt < ptMin[iP] || trackPt > ptMax[iP]) {
          continue;
        }

        auto nSigmaTPC = track.tpcNsigma();
        float mass{track.tofMass()};
        bool hasTof = track.tofMass() > 0;

        if (trackPt <= ptTof[iP] || (trackPt > ptTof[iP] && hasTof && std::abs(mass - partMass[iP]) < tofMassMaxQA)) { // for QA histograms
          tpcNsigmaGlo[iP]->Fill(centrality, trackPt, nSigmaTPC);
          if (nSigmaTPC > nSigmaTpcCutLow[iP] && nSigmaTPC < nSigmaTpcCutUp[iP]) {
            tofMass[iP]->Fill(centrality, trackPt, mass);
          }
        }

        if (nSigmaTPC < nSigmaTpcCutLow[iP] || nSigmaTPC > nSigmaTpcCutUp[iP]) {
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

      auto ptV0 = std::abs(v0.pt());
      if (ptV0 < lambdaPtMin || ptV0 > lambdaPtMax) {
        continue;
      }

      auto etaV0 = v0.eta();
      if (std::abs(etaV0) > etaMax) {
        continue;
      }

      bool matter = v0.pt() > 0;
      auto mLambda = v0.mass();

      // pid selections
      // auto nSigmaTPCPos = v0.tpcNsigmaPos();
      // auto nSigmaTPCNeg = v0.tpcNsigmaNeg();

      // if (std::abs(nSigmaTPCPos) > v0setting_nsigmatpc || std::abs(nSigmaTPCNeg) > v0setting_nsigmatpc) {
      //   continue;
      // }

      float dcaV0dau = v0.dcaV0tracks();
      histos.fill(HIST("QA/dcaV0daughBeforeCut"), dcaV0dau);
      if (dcaV0dau > v0setting_dcav0dau) {
        continue;
      }

      float dcaV0Pv = v0.dcaV0Pv();
      histos.fill(HIST("QA/dcaV0PvBeforeCut"), dcaV0Pv);
      if (std::abs(dcaV0Pv) > v0setting_dcav0pv) {
        continue;
      }

      double cosPA = v0.cosPa();
      histos.fill(HIST("QA/cosPaBeforeCut"), cosPA);
      if (cosPA < v0setting_cospa) {
        continue;
      }

      // auto posDcaToPv = v0.dcaPosPv();
      // if (posDcaToPv < v0setting_dcadaughtopv) {
      //   continue;
      // }

      // auto negDcaToPv = v0.dcaNegPv();
      // if (negDcaToPv < v0setting_dcadaughtopv) {
      //   continue;
      // }

      if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCutQA) { // for QA histograms
        continue;
      }
      histos.fill(HIST("QA/massLambda"), centrality, ptV0, mLambda);

      if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) {
        continue;
      }
      histos.fill(HIST("QA/cosPa"), cosPA);
      histos.fill(HIST("QA/dcaV0daugh"), dcaV0dau);
      // histos.fill(HIST("QA/dcaPosPv"), posDcaToPv);
      // histos.fill(HIST("QA/dcaNegPv"), negDcaToPv);
      histos.fill(HIST("QA/dcaV0Pv"), dcaV0Pv);

      if (matter) {
        tempHistos.fill(HIST("tempLambda"), std::abs(etaV0), ptV0);
      } else {
        tempHistos.fill(HIST("tempAntiLambda"), std::abs(etaV0), ptV0);
      }

      trkId.emplace_back(v0.idNeg());
      trkId.emplace_back(v0.idPos());

      CandidateV0 candV0;
      candV0.pt = ptV0;
      candV0.eta = etaV0;
      candV0.mass = mLambda;
      candV0.cpa = cosPA;
      candV0.dcav0daugh = dcaV0dau;
      candV0.dcav0pv = dcaV0Pv;
      candV0.globalIndexPos = v0.idPos();
      candV0.globalIndexNeg = v0.idNeg();
      candV0.globalIndex = v0.globalIndex();
      candidateV0s.push_back(candV0);
    }

    // reject events having multiple v0s from same tracks (TODO: also across collisions?)
    std::sort(trkId.begin(), trkId.end());
    if (std::adjacent_find(trkId.begin(), trkId.end()) != trkId.end()) {
      candidateV0s.clear();

      CandidateV0 candV0;
      candV0.pt = -999.f;
      candV0.eta = -999.f;
      candV0.mass = -999.f;
      candV0.cpa = -999.f;
      candV0.dcav0daugh = -999.f;
      candV0.dcav0pv = -999.f;
      candV0.globalIndexPos = -999;
      candV0.globalIndexNeg = -999;
      candidateV0s.push_back(candV0);
      return -1;
    }
    for (auto& candidateV0 : candidateV0s) {
      histos.fill(HIST("QA/massLambdaEvRej"), centrality, candidateV0.pt, candidateV0.mass);
    }

    histos.fill(HIST("nEv"), subsample, centrality);

    if (doprocessMc && fillOnlySignal)
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

  template <class C, class T, class V0>
  void fillMcEvent(C const& collision, T const& tracks, V0 const& v0s, float const& centrality)
  {
    int subsample = fillRecoEvent<C, T, V0>(collision, tracks, v0s, centrality);
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
        auto mcTrack = tracks.rawIteratorAt(candidateTrack.globalIndex);
        if (std::abs(mcTrack.pdgCode()) != partPdg[iP])
          continue;
        if (!mcTrack.isReco())
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
    for (auto& candidateV0 : candidateV0s) {
      auto mcTrack = v0s.rawIteratorAt(candidateV0.globalIndex);
      if (!mcTrack.isReco())
        continue;
      if (std::abs(mcTrack.pdgCode()) != 3122) {
        histos.fill(HIST("QA/cosPaBkg"), candidateV0.cpa);
        histos.fill(HIST("QA/dcaV0daughBkg"), candidateV0.dcav0daugh);
        histos.fill(HIST("QA/dcaV0PvBkg"), candidateV0.dcav0pv);
        histos.fill(HIST("QA/cosPaDcaV0daughBkg"), candidateV0.cpa, candidateV0.dcav0daugh);
        histos.fill(HIST("QA/massLambdaEvRejBkg"), centrality, candidateV0.pt, candidateV0.mass);
        continue;
      }
      histos.fill(HIST("QA/cosPaSig"), candidateV0.cpa);
      histos.fill(HIST("QA/dcaV0daughSig"), candidateV0.dcav0daugh);
      histos.fill(HIST("QA/dcaV0PvSig"), candidateV0.dcav0pv);
      histos.fill(HIST("QA/cosPaDcaV0daughSig"), candidateV0.cpa, candidateV0.dcav0daugh);
      histos.fill(HIST("QA/massLambdaEvRejSig"), centrality, candidateV0.pt, candidateV0.mass);
      if (mcTrack.pdgCode() > 0) {
        histos.fill(HIST("recL"), centrality, candidateV0.pt, std::abs(candidateV0.eta));
        if (fillOnlySignal)
          tempLambda->Fill(std::abs(candidateV0.eta), candidateV0.pt);
      } else {
        histos.fill(HIST("recAntiL"), centrality, candidateV0.pt, std::abs(candidateV0.eta));
        if (fillOnlySignal)
          tempAntiLambda->Fill(std::abs(candidateV0.eta), candidateV0.pt);
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

  template <class C, class T, class V0>
  void fillMcGen(C const& /*collision*/, T const& tracks, V0 const& v0s, float const& centrality)
  {

    tempTracks[0]->Reset();
    tempTracks[1]->Reset();
    tempLambda->Reset();
    tempAntiLambda->Reset();

    auto rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
    auto subsample = static_cast<int>(rnd * nSubsamples);
    for (auto& mcPart : v0s) {
      auto genEta = mcPart.genEta();
      if (std::abs(genEta) > etaMax) {
        continue;
      }
      auto pdgCode = mcPart.pdgCode();
      if (std::abs(pdgCode) == 3122) {
        auto genPt = mcPart.genPt();
        if (pdgCode > 0) {
          histos.fill(HIST("genL"), centrality, genPt, std::abs(genEta));
          tempHistos.fill(HIST("tempLambda"), std::abs(genEta), genPt);
        } else {
          histos.fill(HIST("genAntiL"), centrality, genPt, std::abs(genEta));
          tempHistos.fill(HIST("tempAntiLambda"), std::abs(genEta), genPt);
        }
      }
    }
    for (auto& mcPart : tracks) {
      auto genEta = mcPart.genEta();
      if (std::abs(genEta) > etaMax) {
        continue;
      }
      auto pdgCode = mcPart.pdgCode();
      if (std::abs(pdgCode) == partPdg[0] || std::abs(pdgCode) == partPdg[1]) {
        int iP = 1;
        if (std::abs(pdgCode) == partPdg[0]) {
          iP = 0;
        }
        auto genPt = mcPart.genPt();
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

  void processData(aod::CollEbyeTable const& collision, aod::NucleiEbyeTables const& tracks, aod::LambdaEbyeTables const& v0s)
  {
    if (std::abs(collision.zvtx()) > zVtxMax)
      return;
    histos.fill(HIST("QA/zVtx"), collision.zvtx());
    fillRecoEvent(collision, tracks, v0s, collision.centrality());
  }
  PROCESS_SWITCH(nucleiEbye, processData, "process data", false);

  void processMc(aod::CollEbyeTables const& collisions, aod::McNucleiEbyeTables const& tracksTot, aod::McLambdaEbyeTables const& v0sTot)
  {
    for (auto& collision : collisions) {
      if (std::abs(collision.zvtx()) > zVtxMax)
        continue;
      auto tracks = tracksTot.sliceBy(perCollTrack, collision.globalIndex());
      auto v0s = v0sTot.sliceBy(perCollV0s, collision.globalIndex());
      histos.fill(HIST("QA/zVtx"), collision.zvtx());
      fillMcEvent(collision, tracks, v0s, collision.centrality());
      fillMcGen(collision, tracks, v0s, collision.centrality());
    }
  }
  PROCESS_SWITCH(nucleiEbye, processMc, "process Mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiEbye>(cfgc)};
}
