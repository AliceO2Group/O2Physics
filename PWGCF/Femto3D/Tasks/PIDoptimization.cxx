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
/// \brief optimization of particle identification for femtoscopic analysis.
/// \author Sofia Tomassini
/// \since July 2025

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "TRandom3.h"
#include <TParameter.h>

#include "fairlogger/Logger.h"

#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
namespace o2::aod
{
using SelectedTracks = soa::Join<aod::FullTracks, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
}
struct PidOptimization {

  HistogramRegistry histos{"Histos"};

  Configurable<bool> _removeSameBunchPileup{"removeSameBunchPileup", false, ""};
  Configurable<bool> _requestGoodZvtxFT0vsPV{"requestGoodZvtxFT0vsPV", false, ""};
  Configurable<bool> _requestVertexITSTPC{"requestVertexITSTPC", false, ""};
  Configurable<int> _requestVertexTOForTRDmatched{"requestVertexTOFmatched", 0, "0 -> no selectio; 1 -> vertex is matched to TOF or TRD; 2 -> matched to both;"};
  Configurable<bool> _requestNoCollInTimeRangeStandard{"requestNoCollInTimeRangeStandard", false, ""};
  Configurable<std::pair<float, float>> _IRcut{"IRcut", std::pair<float, float>{0.f, 100.f}, "[min., max.] IR range to keep events within"};
  Configurable<std::pair<int, int>> _OccupancyCut{"OccupancyCut", std::pair<int, int>{0, 10000}, "[min., max.] occupancy range to keep events within"};

  Configurable<int> _sign{"sign", 1, "sign of a track"};
  Configurable<float> _vertexZ{"VertexZ", 20.0, "abs vertexZ value limit"};
  Configurable<float> _min_P{"min_P", 0.0, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.0, "upper mometum limit"};
  Configurable<float> _eta{"eta", 100.0, "abs eta value limit"};

  Configurable<std::vector<float>> _dcaXY{"dcaXY", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaXY value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<std::vector<float>> _dcaZ{"dcaZ", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaZ value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<int16_t> _tpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> _tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> _tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<float> _tpcFractionSharedCls{"maxTpcFractionSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> _itsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters for a track"};
  Configurable<float> _itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters for a track"};
  Configurable<int> _particlePDG{"particlePDG", 2212, "PDG code of a particle to perform PID for (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::pair<float, float>> _tpcNSigma{"tpcNSigma", std::pair<float, float>{-100, 100}, "Nsigma range in TPC before the TOF is used"};
  Configurable<std::pair<float, float>> _itsNSigma{"itsNSigma", std::pair<float, float>{-100, 100}, "Nsigma range in ITS to use along with TPC"};
  Configurable<float> _PIDtrshld{"PIDtrshld", 10.0, "value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::pair<float, float>> _tofNSigma{"tofNSigma", std::pair<float, float>{-100, 100}, "Nsigma range in TOF"};
  Configurable<std::pair<float, float>> _tpcNSigmaResidual{"tpcNSigmaResidual", std::pair<float, float>{-10, 10}, "residual TPC Nsigma cut to use with the TOF"};
  Configurable<int> _particlePDGtoReject{"particlePDGtoReject", 211, "PDG codes of perticles that will be rejected with TOF (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::pair<float, float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::pair<float, float>{-10, 10}, "TOF rejection Nsigma range for particles specified with PDG to be rejected"};

  std::shared_ptr<TH2> ITShisto;
  std::shared_ptr<TH2> TPChisto;
  std::shared_ptr<TH2> TOFhisto;
  std::shared_ptr<TH2> ITSvsTPChisto;
  std::shared_ptr<TH2> dcaxy_p;
  std::shared_ptr<TH2> dcaxy_pt;
  std::shared_ptr<TH2> dcaz_p;
  std::shared_ptr<TH2> dcaz_pt;

  void init(o2::framework::InitContext& context)
  {
    o2::aod::ITSResponse::setParameters(context);
    histos.add("vtz", "vtz", kTH1F, {{100, -20., 20., "vtxz"}});
    histos.add("eta", "eta", kTH1F, {{200, -2.5, 2.5, "eta"}});
    histos.add("phi", "phi", kTH1F, {{200, 0., 2. * M_PI, "phi"}});
    histos.add("px", "px", kTH1F, {{100, 0., 5., "px"}});
    histos.add("py", "py", kTH1F, {{100, 0., 5., "py"}});
    histos.add("pz", "pz", kTH1F, {{100, 0., 5., "pz"}});
    histos.add("p", "p", kTH1F, {{100, 0., 5., "p"}});
    histos.add("pt", "pt", kTH1F, {{100, 0., 5., "pt"}});
    histos.add("sign", "sign", kTH1F, {{3, -1.5, 1.5, "sign"}});
    histos.add("TPCClusters", "TPCClusters", kTH1F, {{163, -0.5, 162.5, "NTPCClust"}});
    histos.add("TPCCrossedRowsOverFindableCls", "TPCCrossedRowsOverFindableCls", kTH1F, {{100, 0.0, 10.0, "NcrossedRowsOverFindable"}});
    histos.add("TPCFractionSharedCls", "TPCFractionSharedCls", kTH1F, {{100, 0.0, 1.0, "TPCsharedFraction"}});
    histos.add("ITSClusters", "ITSClusters", kTH1F, {{10, -0.5, 9.5, "NITSClust"}});
    histos.add("ITSchi2", "ITSchi2", kTH1F, {{100, 0.0, 40., "ITSchi2"}});
    histos.add("TPCchi2", "TPCchi2", kTH1F, {{100, 0.0, 6., "TPCchi2"}});

    dcaxy_p = histos.add<TH2>("dcaxy_p", "dcaxy_p", kTH2F, {{100, 0., 5.0, "p"}, {500, -0.5, 0.5, "dcaxy"}});
    dcaxy_pt = histos.add<TH2>("dcaxy_pt", "dcaxy_pt", kTH2F, {{100, 0., 5.0, "pt"}, {500, -0.5, 0.5, "dcaxy"}});
    dcaz_p = histos.add<TH2>("dcaz_p", "dcaz_p", kTH2F, {{100, 0., 5.0, "p"}, {500, -0.5, 0.5, "dcaz"}});
    dcaz_pt = histos.add<TH2>("dcaz_pt", "dcaz_pt", kTH2F, {{100, 0., 5.0, "pt"}, {500, -0.5, 0.5, "dcaxy"}});

    ITShisto = histos.add<TH2>(Form("nsigmaITS_PDG%i", _particlePDG.value), Form("nsigmaITS_PDG%i", _particlePDG.value), kTH2F, {{100, 0., 10.}, {1000, -50., 50.}});
    TPChisto = histos.add<TH2>(Form("nsigmaTPC_PDG%i", _particlePDG.value), Form("nsigmaTPC_PDG%i", _particlePDG.value), kTH2F, {{100, 0., 10.}, {1000, -50., 50.}});
    TOFhisto = histos.add<TH2>(Form("nsigmaTOF_PDG%i", _particlePDG.value), Form("nsigmaTOF_PDG%i", _particlePDG.value), kTH2F, {{100, 0., 10.}, {2000, -100., 100.}});

    ITSvsTPChisto = histos.add<TH2>(Form("nsigmaITSvsTPC_PDG%i", _particlePDG.value), Form("nsigmaITSvsTPC_PDG%i", _particlePDG.value), kTH2F, {{1000, -50., 50.}, {1000, -50., 50.}});
  }

  template <typename TrackType>
  inline float getITSNsigma(TrackType const& track, int const& PDG)
  {
    switch (PDG) {
      case 211:
        return track.itsNSigmaPi();
      case 321:
        return track.itsNSigmaKa();
      case 2212:
        return track.itsNSigmaPr();
      case 1000010020:
        return track.itsNSigmaDe();
      case 0:
        return -1000.0;
      default:
        LOG(fatal) << "Cannot interpret PDG for ITS selection: " << PDG;
        return -1000.0;
    }
  }
  template <typename TrackType>
  inline float getTPCNsigma(TrackType const& track, int const& PDG)
  {
    switch (PDG) {
      case 211:
        return track.tpcNSigmaPi();
      case 321:
        return track.tpcNSigmaKa();
      case 2212:
        return track.tpcNSigmaPr();
      case 1000010020:
        return track.tpcNSigmaDe();
      case 0:
        return -1000.0;
      default:
        LOG(fatal) << "Cannot interpret PDG for TPC selection: " << PDG;
        return -1000.0;
    }
  }
  template <typename TrackType>
  inline float getTOFNsigma(TrackType const& track, int const& PDG)
  {
    switch (PDG) {
      case 211:
        return track.tofNSigmaPi();
      case 321:
        return track.tofNSigmaKa();
      case 2212:
        return track.tofNSigmaPr();
      case 1000010020:
        return track.tofNSigmaDe();
      case 0:
        return -1000.0;
      default:
        LOG(fatal) << "Cannot interpret PDG for TOF selection: " << PDG;
        return -1000.0;
    }
  }

  bool isInRange(float value, std::pair<float, float> range)
  {
    return value > range.first && value < range.second;
  }

  void process(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::SelectedTracks const& tracks)
  {
    auto tracksWithItsPid = soa::Attach<aod::SelectedTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);

    for (auto& collision : collisions) {
      if (!collision.sel8() || abs(collision.posZ()) > _vertexZ)
        continue;
      if (_requestGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX))
        continue;
      if (_removeSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        continue;
      if (_requestGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      if (_requestVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
        continue;
      // if (_requestVertexTOForTRDmatched > collision.selection_bit(o2::aod::evsel::kisVertexTOForTRDmatched()))
      //      continue;
      if (_requestNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
        continue;
      //   if (collision.multPerc() < _centCut.value.first || collision.multPerc() >= _centCut.value.second)
      //     continue;
      //   if (collision.hadronicRate() < _IRcut.value.first || collision.hadronicRate() >= _IRcut.value.second)
      //     continue;
      if (collision.trackOccupancyInTimeRange() < _OccupancyCut.value.first || collision.trackOccupancyInTimeRange() >= _OccupancyCut.value.second)
        continue;

      histos.fill(HIST("vtz"), collision.posZ());
    }

    for (auto& track : tracksWithItsPid) {
      if (track.sign() != _sign)
        continue;
      if (track.p() < _min_P || track.p() > _max_P || abs(track.eta()) > _eta)
        continue;
      if ((track.itsChi2NCl() >= _itsChi2NCl) || (track.itsChi2NCl() <= 0.f) || (track.tpcChi2NCl() <= 0.f) || (track.tpcChi2NCl() >= _tpcChi2NCl))
        continue;
      if ((track.tpcFractionSharedCls()) > _tpcFractionSharedCls || (track.tpcNClsFound()) < _tpcNClsFound || (track.tpcCrossedRowsOverFindableCls()) < _tpcCrossedRowsOverFindableCls || (track.itsNCls()) < _itsNCls)
        continue;
      if (std::fabs(track.dcaXY()) > _dcaXY.value[0] + _dcaXY.value[1] * std::pow(track.pt(), _dcaXY.value[2]) || std::fabs(track.dcaZ()) > _dcaZ.value[0] + _dcaZ.value[1] * std::pow(track.pt(), _dcaZ.value[2]))
        continue;

      bool belowThreshold = track.p() < _PIDtrshld;
      float tpcSigma = getTPCNsigma(track, _particlePDG);
      float itsSigma = getITSNsigma(track, _particlePDG);
      float tofSigma = getTOFNsigma(track, _particlePDG);
      float tofRejection = getTOFNsigma(track, _particlePDGtoReject);

      bool passTPC = belowThreshold ? isInRange(tpcSigma, _tpcNSigma.value) : isInRange(tpcSigma, _tpcNSigmaResidual.value);

      bool passITS = belowThreshold ? isInRange(itsSigma, _itsNSigma.value) : true;

      bool passTOF = belowThreshold ? true : (isInRange(tofSigma, _tofNSigma.value) && !isInRange(tofRejection, _rejectWithinNsigmaTOF.value));

      if (passTPC && passITS && passTOF) {
        histos.fill(HIST("eta"), track.eta());
        histos.fill(HIST("phi"), track.phi());
        histos.fill(HIST("px"), track.px());
        histos.fill(HIST("py"), track.py());
        histos.fill(HIST("pz"), track.pz());
        histos.fill(HIST("p"), track.p());
        histos.fill(HIST("pt"), track.pt());
        histos.fill(HIST("sign"), track.sign());
        histos.fill(HIST("dcaxy_p"), track.p(), track.dcaXY());
        histos.fill(HIST("dcaxy_pt"), track.pt(), track.dcaXY());
        histos.fill(HIST("dcaz_p"), track.p(), track.dcaZ());
        histos.fill(HIST("dcaz_pt"), track.pt(), track.dcaZ());
        histos.fill(HIST("TPCClusters"), track.tpcNClsFound());
        histos.fill(HIST("TPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        histos.fill(HIST("TPCFractionSharedCls"), track.tpcFractionSharedCls());
        histos.fill(HIST("ITSClusters"), track.itsNCls());
        histos.fill(HIST("ITSchi2"), track.itsChi2NCl());
        histos.fill(HIST("TPCchi2"), track.tpcChi2NCl());

        switch (_particlePDG) {
          case 211:
            ITShisto->Fill(track.p(), track.itsNSigmaPi());
            TPChisto->Fill(track.p(), track.tpcNSigmaPi());
            TOFhisto->Fill(track.p(), track.tofNSigmaPi());
            ITSvsTPChisto->Fill(track.itsNSigmaPi(), track.tpcNSigmaPi());
            break;
          case 321:
            ITShisto->Fill(track.p(), track.itsNSigmaKa());
            TPChisto->Fill(track.p(), track.tpcNSigmaKa());
            TOFhisto->Fill(track.p(), track.tofNSigmaKa());
            ITSvsTPChisto->Fill(track.itsNSigmaKa(), track.tpcNSigmaKa());
            break;
          case 2212:
            ITShisto->Fill(track.p(), track.itsNSigmaPr());
            TPChisto->Fill(track.p(), track.tpcNSigmaPr());
            TOFhisto->Fill(track.p(), track.tofNSigmaPr());
            ITSvsTPChisto->Fill(track.itsNSigmaPr(), track.tpcNSigmaPr());
            break;
          case 1000010020:
            ITShisto->Fill(track.p(), track.itsNSigmaDe());
            TPChisto->Fill(track.p(), track.tpcNSigmaDe());
            TOFhisto->Fill(track.p(), track.tofNSigmaDe());
            ITSvsTPChisto->Fill(track.itsNSigmaDe(), track.tpcNSigmaDe());
            break;
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidOptimization>(cfgc)};
}
