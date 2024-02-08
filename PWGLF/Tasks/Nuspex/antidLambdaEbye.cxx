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
} // namespace

struct antidLambdaEbye {
  o2::pid::tof::Beta<TracksFull::iterator> responseBeta;
  std::mt19937 gen32;

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

  std::shared_ptr<TH2> tempAntid;
  std::shared_ptr<TH2> tempAntiLambda;
  std::shared_ptr<TH2> tempLambda;

  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);

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

  void init(o2::framework::InitContext&)
  {

    uint32_t randomSeed = static_cast<uint32_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    gen32.seed(randomSeed);

    histos.add<TH1>("zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});

    auto hNev = histos.add<THnSparse>("nEv", ";Subsample;Centrality (%);", {HistType::kTHnSparseD}, {subsampleAxis, centAxis});
    nSubsamples = hNev->GetAxis(0)->GetNbins();

    histos.add<THnSparse>("nAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis});
    histos.add<THnSparse>("nSqAntid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{d}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptAntidAxis});

    histos.add<THnSparse>("nAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});
    histos.add<THnSparse>("nSqAntiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});

    histos.add<THnSparse>("nL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});
    histos.add<THnSparse>("nSqL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis});

    histos.add<THnSparse>("nLantiL", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptLambdaAxis});
    histos.add<THnSparse>("nLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#Lambda) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});
    histos.add<THnSparse>("nAntiLantid", ";Subsample;Centrality (%);#Delta#eta;#it{p}_{T}(#bar{#Lambda}) (GeV/#it{c});#it{p}_{T}(#bar{d}) (GeV/#it{c});", {HistType::kTHnSparseD}, {subsampleAxis, centAxis, deltaEtaAxis, ptLambdaAxis, ptAntidAxis});

    // v0 QA
    histos.add<TH1>("massLambda", ";#it{M}(p + #pi^{-}) (GeV/#it{c}^{2});Entries", {HistType::kTH1F, {massLambdaAxis}});
    histos.add<TH1>("cosPa", ";cosPa;Entries", {HistType::kTH1F}, {cosPaAxis});
    histos.add<TH1>("radius", ";radius;Entries", {HistType::kTH1F}, {radiusAxis});
    histos.add<TH1>("dcaV0daugh", ";dcaV0daugh;Entries", {HistType::kTH1F}, {dcaV0daughAxis});
    histos.add<TH1>("dcaPosPv", ";dcaPosPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});
    histos.add<TH1>("dcaNegPv", ";dcaNegPv;Entries", {HistType::kTH1F}, {dcaDaughPvAxis});

    // antid QA
    histos.add<TH2>("tpcNsigma", ";#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC} (a.u.)", {HistType::kTH2F}, {momAxis, tpcNsigmaAxis});
    histos.add<TH2>("tpcNsigmaGlo", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC} (a.u.)", {HistType::kTH2F}, {momAxis, tpcNsigmaAxis});
    histos.add<TH2>("tofMass", ";#it{p}_{glo} (GeV/#it{c});Mass (GeV/#it{c}^{2});Entries", {HistType::kTH2F}, {momAxis, tofMassAxis});
    auto hmomCorr = histos.add<TH3>("momCorr", ";#it{p}_{glo} (GeV/#it{c});#it{p}_{TPC} - #it{p}_{glo} (GeV/#it{c});", {HistType::kTH3F}, {momAxisFine, momResAxis, trackingPidAxis});
    histos.add<TH2>("tpcSignal", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", {HistType::kTH2F}, {momAxisFine, tpcAxis});
    histos.add<TH2>("tpcSignalBkg", ";#it{p}_{TPC} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.)", {HistType::kTH2F}, {momAxisFine, tpcAxis});
    histos.add<TH2>("tpcSignal_glo", ";#it{p}_{glo} (GeV/#it{c});d#it{E}/d#it{x}_{TPC} (a.u.);", {HistType::kTH2F}, {momAxisFine, tpcAxis});
    histos.add<TH2>("tofSignal", ";#it{p}_{TPC} (GeV/#it{c});#beta_{TOF}", {HistType::kTH2F}, {momAxisFine, tofAxis});
    histos.add<TH2>("tofSignal_glo", ";#it{p}_{T} (GeV/#it{c});#beta_{TOF}", {HistType::kTH2F}, {momAxisFine, tofAxis});

    for (int i{1}; i < hmomCorr->GetNbinsZ() + 1; ++i) {
      hmomCorr->GetZaxis()->SetBinLabel(i, pidHypotheses[i - 1].data());
    }

    // temporary histograms
    tempAntid = tempHistos.add<TH2>("tempAntid", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2D}, {deltaEtaAxis, ptAntidAxis});
    tempLambda = tempHistos.add<TH2>("tempLambda", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2D}, {deltaEtaAxis, ptLambdaAxis});
    tempAntiLambda = tempHistos.add<TH2>("tempAntiLambda", ";#Delta#eta;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2D}, {deltaEtaAxis, ptLambdaAxis});
  }

  void fillEvent(TracksFull const& tracks, soa::Filtered<aod::V0Datas> const& V0s, float const& centrality)
  {
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
    }

    // reject events having multiple v0s from same tracks (TODO: also across collisions?)
    std::sort(trkId.begin(), trkId.end());
    if (std::adjacent_find(trkId.begin(), trkId.end()) != trkId.end()) {
      return;
    }

    histos.fill(HIST("nEv"), subsample, centrality);

    for (int iEta{1}; iEta < tempAntid->GetNbinsX() + 1; ++iEta) {
      for (int iPtAntid{1}; iPtAntid < tempAntid->GetNbinsY() + 1; ++iPtAntid) {
        for (int iPtL{1}; iPtL < tempLambda->GetNbinsY() + 1; ++iPtL) {
          for (int iPtAntiL{1}; iPtAntiL < tempAntiLambda->GetNbinsY() + 1; ++iPtAntiL) {

            auto eta = tempAntid->GetXaxis()->GetBinCenter(iEta);
            auto antidPt = tempAntid->GetYaxis()->GetBinCenter(iPtAntid);
            auto lPt = tempLambda->GetYaxis()->GetBinCenter(iPtL);
            auto antilPt = tempAntiLambda->GetYaxis()->GetBinCenter(iPtAntiL);
            auto antidNum = tempAntid->Integral(1, iEta, iPtAntid, iPtAntid);
            auto lNum = tempLambda->Integral(1, iEta, iPtL, iPtL);
            auto antilNum = tempAntiLambda->Integral(1, iEta, iPtAntiL, iPtAntiL);

            if (iPtL == 1) {
              histos.fill(HIST("nAntiLantid"), subsample, centrality, eta, antilPt, antidPt, antilNum * antidNum);
              if (iPtAntiL == 1) {
                histos.fill(HIST("nAntid"), subsample, centrality, eta, antidPt, antidNum);
                histos.fill(HIST("nSqAntid"), subsample, centrality, eta, antidPt, std::pow(antidNum, 2));
              }
            }
            if (iPtAntiL == 1) {
              histos.fill(HIST("nLantid"), subsample, centrality, eta, lPt, antidPt, lNum * antidNum);
              if (iPtAntid == 1) {
                histos.fill(HIST("nL"), subsample, centrality, eta, lPt, lNum);
                histos.fill(HIST("nSqL"), subsample, centrality, eta, lPt, std::pow(lNum, 2));
              }
            }
            if (iPtAntid == 1) {
              histos.fill(HIST("nLantiL"), subsample, centrality, eta, lPt, antilPt, lNum * antilNum);
              if (iPtL == 1) {
                histos.fill(HIST("nAntiL"), subsample, centrality, eta, antilPt, antilNum);
                histos.fill(HIST("nSqAntiL"), subsample, centrality, eta, antilPt, std::pow(antilNum, 2));
              }
            }
          }
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>::iterator const& collision, TracksFull const& tracks, soa::Filtered<aod::V0Datas> const& V0s)
  {
    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    histos.fill(HIST("zVtx"), collision.posZ());
    auto centrality = collision.centFT0C();
    fillEvent(tracks, V0s, centrality);
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
    fillEvent(tracks, V0s, centrality);
  }
  PROCESS_SWITCH(antidLambdaEbye, processRun2, "process (Run 2)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<antidLambdaEbye>(cfgc)};
}
