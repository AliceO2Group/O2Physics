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

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCDe, aod::pidTPCPr, aod::pidTPCPi, aod::TOFSignal, aod::TOFEvTime>;

namespace
{
constexpr double betheBlochDefault[1][6]{{-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNamesBB{"d"};
static const std::vector<std::string> pidHypotheses{"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "He3", "Alpha", "Pion0", "Photon", "K0", "Lambda", "HyperTriton", "Hyperhydrog4", "XiMinus", "OmegaMinus"};
std::shared_ptr<TH2> tempAntid;
std::shared_ptr<TH2> tempAntiLambda;
std::shared_ptr<TH2> tempLambda;
std::shared_ptr<THnSparse> nAntid;
std::shared_ptr<THnSparse> nAntiL;
std::shared_ptr<THnSparse> nL;
std::shared_ptr<THnSparse> nSqAntid;
std::shared_ptr<THnSparse> nSqAntiL;
std::shared_ptr<THnSparse> nSqL;
std::shared_ptr<THnSparse> nLantiL;
std::shared_ptr<THnSparse> nLantid;
std::shared_ptr<THnSparse> nAntiLantid;
} // namespace

struct CandidateV0 {
  float pt;
  float eta;
  int64_t globalIndexPos = -999;
  int64_t globalIndexNeg = -999;
};

struct CandidateAntid {
  float pt;
  float eta;
  int64_t globalIndex = -999;
};

struct antidLambdaEbye {
  o2::pid::tof::Beta<TracksFull::iterator> responseBeta;
  std::mt19937 gen32;
  std::vector<CandidateV0> candidateV0s;
  std::vector<CandidateAntid> candidateAntids;

  int nSubsamples;

  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNamesBB, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for deuteron"};

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis subsampleAxis{"subsampleAxis", {30, 0, 30}, "binning of the subsample axis"};
  ConfigurableAxis deltaEtaAxis{"deltaEtaAxis", {4, 0, 0.8}, "binning of the delta eta axis"};
  ConfigurableAxis ptAntidAxis{"ptAntidAxis", {VARIABLE_WIDTH, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f}, "binning of the antideuteron pT axis (GeV/c)"};
  ConfigurableAxis ptLambdaAxis{"ptLambdaAxis", {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f}, "binning of the (anti)lambda pT axis (GeV/c)"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};

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
  ConfigurableAxis trackingPidAxis{"trackingPidAxis", {static_cast<double>(pidHypotheses.size()), 0, static_cast<double>(pidHypotheses.size())}, "tracking pid hypothesis binning"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.8f, "maximum eta"};

  Configurable<float> antidPtMin{"antidPtMin", 0.8f, "minimum antideuteron pT (GeV/c)"};
  Configurable<float> antidPtTof{"antidPtTof", 1.0f, "antideuteron pT to switch to TOF pid (GeV/c) "};
  Configurable<float> antidPtMax{"antidPtMax", 1.8f, "maximum antideuteron pT (GeV/c)"};

  Configurable<float> lambdaPtMin{"lambdaPtMin", 0.5f, "minimum (anti)lambda pT (GeV/c)"};
  Configurable<float> lambdaPtMax{"lambdaPtMax", 3.0f, "maximum (anti)lambda pT (GeV/c)"};

  Configurable<float> antidNcrossedRows{"antidNcrossedRows", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> antidNclusItsCut{"antidNclusITScut", 5, "Minimum number of ITS clusters"};
  Configurable<float> antidNclusTpcCut{"antidNclusTPCcut", 70, "Minimum number of TPC clusters"};
  Configurable<float> antidNsigmaTpcCut{"antidNsigmaTpcCut", 4.f, "TPC PID cut"};
  Configurable<float> antidNsigmaTofCut{"antidNsigmaTofCut", 4.f, "TOF PID cut"};
  Configurable<float> antidDcaCut{"antidDcaCut", 0.1f, "DCA antid to PV"};
  Configurable<float> tpcInnerParamMax{"tpcInnerParamMax", 0.6f, "(temporary) tpc inner param cut"};
  Configurable<float> tofMassMax{"tofMassMax", 0.3f, "(temporary) tof mass cut"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1f, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1f, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5f, "v0radius"};
  Configurable<float> lambdaMassCut{"lambdaMassCut", 0.005f, "maximum deviation from PDG mass"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry tempHistos{"tempHistos", {}, OutputObjHandlingPolicy::TransientObject};

  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  Preslice<TracksFull> perCollisionTracksFull = o2::aod::track::collisionId;

  template <class RecV0>
  bool selectLambda(RecV0 const& v0) // TODO: apply ML
  {
    if (std::abs(v0.eta()) > etaMax ||
        v0.v0cosPA() < v0setting_cospa ||
        v0.v0radius() < v0setting_radius) {
      return false;
    }
    auto mLambda = v0.alpha() > 0 ? v0.mLambda() : v0.mAntiLambda();
    if (std::abs(mLambda - o2::constants::physics::MassLambda0) > lambdaMassCut) {
      return false;
    }
    return true;
  }

  template <class T>
  bool selectAntid(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (track.sign() > 0.) {
      return false;
    }
    if (track.itsNCls() < antidNclusItsCut ||
        track.tpcNClsFound() < antidNclusTpcCut ||
        track.tpcNClsCrossedRows() < antidNcrossedRows ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcChi2NCl() > 4.f ||
        track.itsChi2NCl() > 36.f) {
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

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    auto hNev = histos.add<THnSparse>("nEv", ";Subsample;Centrality (%);", HistType::kTHnSparseD, {subsampleAxis, centAxis});
    nSubsamples = hNev->GetAxis(0)->GetNbins();

    nAntid = histos.add<THnSparse>("nAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis});
    nAntiL = histos.add<THnSparse>("nAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});
    nL = histos.add<THnSparse>("nL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});

    nSqAntid = histos.add<THnSparse>("nSqAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis, ptAntidAxis});
    nSqAntiL = histos.add<THnSparse>("nSqAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    nSqL = histos.add<THnSparse>("nSqL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#Lambda) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});

    nLantiL = histos.add<THnSparse>("nLantiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    nLantid = histos.add<THnSparse>("nLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});
    nAntiLantid = histos.add<THnSparse>("nAntiLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", HistType::kTHnSparseD, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});

    // v0 QA
    histos.add<TH1>("massLambda", ";#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH1F, {massLambdaAxis});
    histos.add<TH1>("cosPa", ";cosPa;Entries", HistType::kTH1F, {cosPaAxis});
    histos.add<TH1>("radius", ";radius;Entries", HistType::kTH1F, {radiusAxis});
    histos.add<TH1>("dcaV0daugh", ";dcaV0daugh;Entries", HistType::kTH1F, {dcaV0daughAxis});
    histos.add<TH1>("dcaPosPv", ";dcaPosPv;Entries", HistType::kTH1F, {dcaDaughPvAxis});
    histos.add<TH1>("dcaNegPv", ";dcaNegPv;Entries", HistType::kTH1F, {dcaDaughPvAxis});

    // antid QA
    histos.add<TH2>("tpcNsigma", ";#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcNsigmaAxis});
    histos.add<TH2>("tpcNsigmaGlo", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC} (a.u.)", HistType::kTH2F, {momAxis, tpcNsigmaAxis});
    histos.add<TH2>("tofMass", ";#it{p}_{glo} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, tofMassAxis});
    auto hmomCorr = histos.add<TH3>("momCorr", ";#it{p}_{glo} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c});", HistType::kTH3F, {momAxisFine, momResAxis, trackingPidAxis});
    histos.add<TH2>("tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", HistType::kTH2F, {momAxisFine, tpcAxis});
    histos.add<TH2>("tpcSignalBkg", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", HistType::kTH2F, {momAxisFine, tpcAxis});
    histos.add<TH2>("tpcSignal_glo", ";#it{p}_{glo} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.);", HistType::kTH2F, {momAxisFine, tpcAxis});
    histos.add<TH2>("tofSignal", ";#it{p}_{TPC} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {momAxisFine, tofAxis});
    histos.add<TH2>("tofSignal_glo", ";#it{p}_{T} (GeV/#it{c});#beta_{TOF}", HistType::kTH2F, {momAxisFine, tofAxis});

    if (doprocessMcRun3) {
      histos.add<TH3>("recL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      histos.add<TH3>("recAntiL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      histos.add<TH3>("recD", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
      histos.add<TH3>("recAntid", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
      histos.add<TH3>("genL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      histos.add<TH3>("genAntiL", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptLambdaAxis, deltaEtaAxis});
      histos.add<TH3>("genD", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
      histos.add<TH3>("genAntid", ";Centrality (%); #it{p}_{T} (GeV/#it{c});#Delta#eta", HistType::kTH3D, {centAxis, ptAntidAxis, deltaEtaAxis});
    }

    for (int i{1}; i < hmomCorr->GetNbinsZ() + 1; ++i) {
      hmomCorr->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
    }

    // temporary histograms
    tempAntid = tempHistos.add<TH2>("tempAntid", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptAntidAxis});
    tempLambda = tempHistos.add<TH2>("tempLambda", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptLambdaAxis});
    tempAntiLambda = tempHistos.add<TH2>("tempAntiLambda", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", HistType::kTH2D, {deltaEtaAxis, ptLambdaAxis});
  }

  void fillRecoEvent(TracksFull const& tracks, aod::V0Datas const& V0s, float const& centrality)
  {
    candidateAntids.clear();
    candidateV0s.clear();

    tempAntid->Reset();
    tempLambda->Reset();
    tempAntiLambda->Reset();
    auto rnd = static_cast<float>(gen32()) / static_cast<float>(gen32.max());
    auto subsample = static_cast<int>(rnd * nSubsamples);

    for (const auto& track : tracks) {
      if (!selectAntid(track)) {
        continue;
      }

      if (std::hypot(track.dcaXY(), track.dcaZ()) > antidDcaCut) {
        continue;
      }

      histos.fill(HIST("tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
      histos.fill(HIST("tpcSignal_glo"), track.p(), track.tpcSignal());

      if (track.pt() < antidPtMin || track.pt() > antidPtMax) {
        continue;
      }

      double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / constants::physics::MassDeuteron), cfgBetheBlochParams->get("d", "p0"), cfgBetheBlochParams->get("d", "p1"), cfgBetheBlochParams->get("d", "p2"), cfgBetheBlochParams->get("d", "p3"), cfgBetheBlochParams->get("d", "p4"))};
      double expSigma{expBethe * cfgBetheBlochParams->get("d", "resolution")};
      auto nSigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);

      float beta{track.hasTOF() ? responseBeta.GetBeta(track) : -999.f};
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta));
      float mass{track.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f)};
      bool hasTof = track.hasTOF() && track.tofChi2() < 3;

      if (std::abs(nSigmaTPC) > antidNsigmaTpcCut) {
        continue;
      }
      histos.fill(HIST("tpcNsigma"), track.tpcInnerParam(), nSigmaTPC);
      histos.fill(HIST("momCorr"), track.p(), track.tpcInnerParam() - track.p(), track.pidForTracking());
      // check contamination
      if (track.tpcInnerParam() < tpcInnerParamMax) {
        histos.fill(HIST("tpcSignalBkg"), track.tpcInnerParam(), track.tpcSignal());
      }

      if (track.pt() > antidPtTof && hasTof) {
        histos.fill(HIST("tofSignal_glo"), track.p(), beta);
        histos.fill(HIST("tofSignal"), track.tpcInnerParam(), beta);
      }

      // temporary cut to reject fake matches
      if (track.tpcInnerParam() < tpcInnerParamMax) {
        continue;
      }
      if (track.pt() > antidPtTof && !hasTof) {
        continue;
      }
      histos.fill(HIST("tofMass"), track.pt(), mass);

      if (track.pt() <= antidPtTof || (track.pt() > antidPtTof && hasTof && std::abs(mass - o2::constants::physics::MassDeuteron) < tofMassMax)) {
        histos.fill(HIST("tpcNsigmaGlo"), track.pt(), nSigmaTPC);
        tempHistos.fill(HIST("tempAntid"), std::abs(track.eta()), track.pt());
        CandidateAntid candAntid;
        candAntid.pt = track.pt();
        candAntid.eta = track.eta();
        candAntid.globalIndex = track.globalIndex();
        candidateAntids.push_back(candAntid);
      }
    }

    std::vector<int64_t> trkId;
    for (const auto& v0 : V0s) {
      if (v0.pt() < lambdaPtMin || v0.pt() > lambdaPtMax) {
        continue;
      }

      if (!selectLambda(v0)) {
        continue;
      }

      auto pos = v0.template posTrack_as<TracksFull>();
      auto neg = v0.template negTrack_as<TracksFull>();
      if (std::abs(pos.eta()) > etaMax || std::abs(neg.eta()) > etaMax) {
        continue;
      }

      bool matter = v0.alpha() > 0;

      histos.fill(HIST("massLambda"), matter ? v0.mLambda() : v0.mAntiLambda());
      histos.fill(HIST("cosPa"), v0.v0cosPA());
      histos.fill(HIST("radius"), v0.v0radius());
      histos.fill(HIST("dcaV0daugh"), v0.dcaV0daughters());
      histos.fill(HIST("dcaPosPv"), v0.dcapostopv());
      histos.fill(HIST("dcaNegPv"), v0.dcanegtopv());

      if (matter) {
        tempHistos.fill(HIST("tempLambda"), std::abs(v0.eta()), v0.pt());
      } else {
        tempHistos.fill(HIST("tempAntiLambda"), std::abs(v0.eta()), v0.pt());
      }

      trkId.emplace_back(pos.globalIndex());
      trkId.emplace_back(neg.globalIndex());

      CandidateV0 candV0;
      candV0.pt = v0.pt();
      candV0.eta = v0.eta();
      candV0.globalIndexPos = pos.globalIndex();
      candV0.globalIndexNeg = neg.globalIndex();
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
      return;
    }

    histos.fill(HIST("nEv"), subsample, centrality);

    fillHistoN(nAntid, tempAntid, subsample, centrality);
    fillHistoN(nAntiL, tempAntiLambda, subsample, centrality);
    fillHistoN(nL, tempLambda, subsample, centrality);

    fillHistoN(nSqAntid, tempAntid, tempAntid, subsample, centrality);
    fillHistoN(nSqAntiL, tempAntiLambda, tempAntiLambda, subsample, centrality);
    fillHistoN(nSqL, tempLambda, tempLambda, subsample, centrality);

    fillHistoN(nLantid, tempLambda, tempAntid, subsample, centrality);
    fillHistoN(nLantiL, tempLambda, tempAntiLambda, subsample, centrality);
    fillHistoN(nAntiLantid, tempAntiLambda, tempAntid, subsample, centrality);
  }

  void fillMcEvent(TracksFull const& tracks, aod::V0Datas const& V0s, float const& centrality, aod::McParticles const&, aod::McTrackLabels const& mcLabels)
  {
    fillRecoEvent(tracks, V0s, centrality);
    if (candidateV0s.size() == 1 && candidateV0s[0].pt < -998.f && candidateV0s[0].eta < -998.f && candidateV0s[0].globalIndexPos == -999 && candidateV0s[0].globalIndexPos == -999) {
      return;
    }
    for (auto& candidateAntid : candidateAntids) {
      auto mcLab = mcLabels.rawIteratorAt(candidateAntid.globalIndex);
      if (mcLab.has_mcParticle()) {
        auto mcTrack = mcLab.template mcParticle_as<aod::McParticles>();
        if (std::abs(mcTrack.pdgCode()) != o2::constants::physics::kDeuteron)
          continue;
        if (!mcTrack.isPhysicalPrimary())
          continue;
        if (mcTrack.pdgCode() > 0) {
          histos.fill(HIST("recD"), centrality, candidateAntid.pt, candidateAntid.eta);
        } else {
          histos.fill(HIST("recAntid"), centrality, candidateAntid.pt, candidateAntid.eta);
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
              if (std::abs(posMother.pdgCode()) != 3122)
                continue;
              if (!posMother.has_mothers())
                continue;

              if (posMother.pdgCode() > 0) {
                histos.fill(HIST("recL"), centrality, candidateV0.pt, candidateV0.eta);
              } else {
                histos.fill(HIST("recAntiL"), centrality, candidateV0.pt, candidateV0.eta);
              }
            }
          }
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions, TracksFull const& tracks, soa::Filtered<aod::V0Datas> const& V0s)
  {
    for (const auto& collision : collisions) {
      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      histos.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perCollisionTracksFull, collIdx);
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      auto centrality = collision.centFT0C();
      fillRecoEvent(TrackTable_thisCollision, V0Table_thisCollision, centrality);
    }
  }
  PROCESS_SWITCH(antidLambdaEbye, processRun3, "process (Run 3)", false);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::FV0Mults, aod::CentRun2V0Ms>::iterator const& collision, TracksFull const& tracks, soa::Filtered<aod::V0Datas> const& V0s)
  {
    if (!collision.sel7())
      return;

    if (!collision.alias_bit(kINT7))
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    histos.fill(HIST("zVtx"), collision.posZ());
    auto centrality = collision.centRun2V0M();
    fillRecoEvent(tracks, V0s, centrality);
  }
  PROCESS_SWITCH(antidLambdaEbye, processRun2, "process (Run 2)", false);

  void processMcRun3(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs> const& collisions, aod::McCollisions const& mcCollisions, TracksFull const& tracks, soa::Filtered<aod::V0Datas> const& V0s, aod::McParticles const& mcParticles, aod::McTrackLabels const& mcLab)
  {
    std::vector<std::pair<bool, float>> goodCollisions(mcCollisions.size(), std::make_pair(false, -999.));
    for (auto& collision : collisions) {
      if (!collision.sel8())
        continue;

      if (std::abs(collision.posZ()) > zVtxMax)
        continue;

      auto centrality = collision.centFT0C();
      goodCollisions[collision.mcCollisionId()].first = true;
      goodCollisions[collision.mcCollisionId()].second = centrality;

      histos.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(perCollisionTracksFull, collIdx);
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionV0, collIdx);
      V0Table_thisCollision.bindExternalIndices(&tracks);

      // fillMC(tracks, V0s, centrality);
      fillMcEvent(TrackTable_thisCollision, V0Table_thisCollision, centrality, mcParticles, mcLab);
      if (candidateV0s.size() == 1 && candidateV0s[0].pt < -998.f && candidateV0s[0].eta < -998.f && candidateV0s[0].globalIndexPos == -999 && candidateV0s[0].globalIndexPos == -999) {
        goodCollisions[collision.mcCollisionId()].first = false;
      }
    }

    for (auto& mcPart : mcParticles) {
      auto genEta = mcPart.eta();
      if (std::abs(genEta) > etaMax) {
        continue;
      }
      if (goodCollisions[mcPart.mcCollisionId()].first == false) {
        continue;
      }
      auto centrality = goodCollisions[mcPart.mcCollisionId()].second;
      auto pdgCode = mcPart.pdgCode();

      if (std::abs(pdgCode) == 3122) {
        if (!mcPart.has_mothers())
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
        } else {
          histos.fill(HIST("genAntiL"), centrality, genPt, std::abs(genEta));
        }
      } else if (std::abs(pdgCode) == o2::constants::physics::kDeuteron) {
        if (!mcPart.isPhysicalPrimary())
          continue;
        std::cout << pdgCode << std::endl;
        auto genPt = std::hypot(mcPart.px(), mcPart.py());
        if (pdgCode > 0) {
          histos.fill(HIST("genD"), centrality, genPt, std::abs(genEta));
        } else {
          histos.fill(HIST("genAntid"), centrality, genPt, std::abs(genEta));
        }
      }
    }
  }
  PROCESS_SWITCH(antidLambdaEbye, processMcRun3, "process MC (Run 3)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<antidLambdaEbye>(cfgc)};
}
