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
/// \brief Femto3D QA task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include <vector>
#include <memory>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <TParameter.h>
#include <TH1F.h>

#include "Framework/ASoA.h"
#include "MathUtils/Utils.h"
#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"

#include "Framework/StaticFor.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"
#include "PWGCF/Femto3D/Core/femto3dPairTask.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct QAHistograms {
  // using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> _removeSameBunchPileup{"removeSameBunchPileup", false, ""};
  Configurable<bool> _requestGoodZvtxFT0vsPV{"requestGoodZvtxFT0vsPV", false, ""};
  Configurable<bool> _requestVertexITSTPC{"requestVertexITSTPC", false, ""};
  Configurable<int> _requestVertexTOForTRDmatched{"requestVertexTOFmatched", 0, "0 -> no selectio; 1 -> vertex is matched to TOF or TRD; 2 -> matched to both;"};
  Configurable<bool> _requestNoCollInTimeRangeStandard{"requestNoCollInTimeRangeStandard", false, ""};
  Configurable<std::pair<float, float>> _IRcut{"IRcut", std::pair<float, float>{0.f, 100.f}, "[min., max.] IR range to keep events within"};
  Configurable<std::pair<int, int>> _OccupancyCut{"OccupancyCut", std::pair<int, int>{0, 10000}, "[min., max.] occupancy range to keep events within"};

  Configurable<int> _sign{"sign", 1, "sign of a track"};
  Configurable<float> _vertexZ{"VertexZ", 10.0, "abs vertexZ value limit"};
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
  Configurable<std::vector<float>> _tpcNSigma{"tpcNSigma", std::vector<float>{-4.0f, 4.0f}, "Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld{"PIDtrshld", 10.0, "value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma{"tofNSigma", std::vector<float>{-4.0f, 4.0f}, "Nsigma range in TOF"};
  Configurable<std::vector<float>> _tpcNSigmaResidual{"tpcNSigmaResidual", std::vector<float>{-5.0f, 5.0f}, "residual TPC Nsigma cut to use with the TOF"};

  Configurable<int> _particlePDGtoReject{"particlePDGtoReject", 211, "PDG codes of perticles that will be rejected with TOF (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for particles specified with PDG to be rejected"};

  Configurable<std::pair<float, float>> _centCut{"centCut", std::pair<float, float>{0.f, 100.f}, "[min., max.] centrality range to keep tracks within"};

  Configurable<std::vector<float>> _dcaBinning{"dcaBinning", std::vector<float>{501, 0.5f, 1}, "setup for variable binning (geometric progression is used): 1st (int) -- N_bins (must be odd, otherwise will be increased by 1); 2nd (float) -- abs value of the edge of axises in histos (-2nd, +2nd); 3d (int) -- desired ratio between w_bin at the edges and at 0;"};

  std::pair<int, std::vector<float>> TPCcuts;
  std::pair<int, std::vector<float>> TOFcuts;

  Filter signFilter = o2::aod::singletrackselector::sign == _sign;
  Filter pFilter = o2::aod::singletrackselector::p > _min_P&& o2::aod::singletrackselector::p < _max_P;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < _eta;

  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= _tpcNClsFound &&
                        o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) < _tpcChi2NCl &&
                        o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) > _tpcCrossedRowsOverFindableCls;

  Filter itsTrkFilter = o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) < _itsChi2NCl;

  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < _vertexZ;

  void init(o2::framework::InitContext&)
  {
    TPCcuts = std::make_pair(_particlePDG, _tpcNSigma);
    TOFcuts = std::make_pair(_particlePDG, _tofNSigma);

    int N = _dcaBinning.value[0]; // number of bins -- must be odd otherwise will be increased by 1
    if (N % 2 != 1) {
      N += 1;
    }

    std::unique_ptr<double[]> dca_bins;
    if (static_cast<int>(_dcaBinning.value[2]) != 1.0) {
      dca_bins = calc_var_bins(N + 1, _dcaBinning.value[1], static_cast<int>(_dcaBinning.value[2]));
    } else {
      dca_bins = calc_const_bins(N, -_dcaBinning.value[1], _dcaBinning.value[1]);
    }
    auto const_bins_p = calc_const_bins(100, 0., 5.0);

    auto DCA_XY_p = registry.add<TH2>("dcaxy_to_p", "dcaxy_to_p", kTH2F, {{100, 0., 5.0, "p"}, {501, -0.5, 0.5, "dcaxy"}});
    auto DCA_XY_pt = registry.add<TH2>("dcaxy_to_pt", "dcaxy_to_pt", kTH2F, {{100, 0., 5.0, "pt"}, {501, -0.5, 0.5, "dcaxy"}});

    auto DCA_Z_p = registry.add<TH2>("dcaz_to_p", "dcaz_to_p", kTH2F, {{100, 0., 5.0, "p"}, {501, -0.5, 0.5, "dcaz"}});
    auto DCA_Z_pt = registry.add<TH2>("dcaz_to_pt", "dcaz_to_pt", kTH2F, {{100, 0., 5.0, "pt"}, {501, -0.5, 0.5, "dcaz"}});

    DCA_XY_p->SetBins(100, &const_bins_p[0], N, &dca_bins[0]);  // set variable bins in Y and Z axis; constant on X
    DCA_XY_pt->SetBins(100, &const_bins_p[0], N, &dca_bins[0]); // set variable bins in Y and Z axis; constant on X
    DCA_Z_p->SetBins(100, &const_bins_p[0], N, &dca_bins[0]);   // set variable bins in Y and Z axis; constant on X
    DCA_Z_pt->SetBins(100, &const_bins_p[0], N, &dca_bins[0]);  // set variable bins in Y and Z axis; constant on X

    registry.add("TPCSignal_nocuts", "TPC signal without cuts", kTH2F, {{{200, 0., 5.0, "#it{p}_{inner} (GeV/#it{c})"}, {1000, 0., 1000.0, "dE/dx in TPC (arbitrary units)"}}});
    registry.add("TOFSignal_nocuts", "TOF signal without cuts", kTH2F, {{{200, 0., 5.0, "#it{p} (GeV/#it{c})"}, {100, 0., 1.5, "#beta"}}});

    registry.add("eta", "eta", kTH1F, {{200, -2.5, 2.5, "eta"}});
    registry.add("phi", "phi", kTH1F, {{200, 0., 2. * M_PI, "phi"}});
    registry.add("px", "px", kTH1F, {{100, 0., 5., "px"}});
    registry.add("py", "py", kTH1F, {{100, 0., 5., "py"}});
    registry.add("pz", "pz", kTH1F, {{100, 0., 5., "pz"}});
    registry.add("p", "p", kTH1F, {{100, 0., 5., "p"}});
    registry.add("pt", "pt", kTH1F, {{100, 0., 5., "pt"}});
    registry.add("sign", "sign", kTH1F, {{3, -1.5, 1.5, "sign"}});
    registry.add("TPCClusters", "TPCClusters", kTH1F, {{163, -0.5, 162.5, "NTPCClust"}});
    registry.add("TPCCrossedRowsOverFindableCls", "TPCCrossedRowsOverFindableCls", kTH1F, {{100, 0.0, 10.0, "NcrossedRowsOverFindable"}});
    registry.add("TPCFractionSharedCls", "TPCFractionSharedCls", kTH1F, {{100, 0.0, 1.0, "TPCsharedFraction"}});
    registry.add("ITSClusters", "ITSClusters", kTH1F, {{10, -0.5, 9.5, "NITSClust"}});
    registry.add("ITSchi2", "ITSchi2", kTH1F, {{100, 0.0, 40., "ITSchi2"}});
    registry.add("TPCchi2", "TPCchi2", kTH1F, {{100, 0.0, 6., "TPCchi2"}});

    registry.add("TPCSignal", "TPC Signal", kTH2F, {{{200, 0., 5.0, "#it{p}_{inner} (GeV/#it{c})"}, {1000, 0., 1000.0, "dE/dx in TPC (arbitrary units)"}}});
    registry.add("TOFSignal", "TOF Signal", kTH2F, {{200, 0., 5.0, "#it{p} (GeV/#it{c})"}, {100, 0., 1.5, "#beta"}});

    switch (_particlePDG) {
      case 211:
        registry.add("nsigmaTOFPi", "nsigmaTOFPi", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaTPCPi", "nsigmaTPCPi", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaITSPi", "nsigmaITSPi", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        break;
      case 321:
        registry.add("nsigmaTOFKa", "nsigmaTOFKa", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaTPCKa", "nsigmaTPCKa", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaITSKa", "nsigmaITSKa", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        break;
      case 2212:
        registry.add("nsigmaTOFPr", "nsigmaTOFPr", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaTPCPr", "nsigmaTPCPr", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaITSPr", "nsigmaITSPr", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        break;
      case 1000010020:
        registry.add("nsigmaTOFDe", "nsigmaTOFDe", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaTPCDe", "nsigmaTPCDe", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaITSDe", "nsigmaITSDe", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        break;
      case 1000010030:
        registry.add("nsigmaTOFTr", "nsigmaTOFTr", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaTPCTr", "nsigmaTPCTr", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaITSTr", "nsigmaITSTr", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        break;
      case 1000020030:
        registry.add("nsigmaTOFHe", "nsigmaTOFHe", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaTPCHe", "nsigmaTPCHe", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        registry.add("nsigmaITSHe", "nsigmaITSHe", kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
        break;
      default:
        break;
    }

    const AxisSpec axisMult{5001, -0.5, 5000.5, "mult."};
    const AxisSpec axisPerc{101, -0.5, 100.5, "percentile"};
    registry.add("posZ", "posZ", kTH1F, {{300, -16., 16., "posZ"}});
    registry.add("mult", "mult", kTH1F, {axisMult});
    registry.add("MultVsCent", "MultVsCent", kTH2F, {axisMult, axisPerc});
    registry.add("IRvsOccupancy", "IRvsOccupancy", kTH2I, {{10000, 0, 10000, "Occupancy"}, {50, 0, 50, "IR, kHz"}});
  }

  template <bool FillExtra, typename ColsType, typename TracksType>
  void fillHistograms(ColsType const& collisions, TracksType const& tracks)
  {
    for (auto& collision : collisions) {
      if (_removeSameBunchPileup && !collision.isNoSameBunchPileup())
        continue;
      if (_requestGoodZvtxFT0vsPV && !collision.isGoodZvtxFT0vsPV())
        continue;
      if (_requestVertexITSTPC && !collision.isVertexITSTPC())
        continue;
      if (_requestVertexTOForTRDmatched > collision.isVertexTOForTRDmatched())
        continue;
      if (_requestNoCollInTimeRangeStandard && !collision.noCollInTimeRangeStandard())
        continue;
      if (collision.multPerc() < _centCut.value.first || collision.multPerc() >= _centCut.value.second)
        continue;
      if (collision.hadronicRate() < _IRcut.value.first || collision.hadronicRate() >= _IRcut.value.second)
        continue;
      if (collision.occupancy() < _OccupancyCut.value.first || collision.occupancy() >= _OccupancyCut.value.second)
        continue;

      registry.fill(HIST("posZ"), collision.posZ());
      registry.fill(HIST("mult"), collision.mult());
      registry.fill(HIST("MultVsCent"), collision.mult(), collision.multPerc());
      registry.fill(HIST("IRvsOccupancy"), collision.occupancy(), collision.hadronicRate());
    }

    for (auto& track : tracks) {

      if (_removeSameBunchPileup && !track.template singleCollSel_as<ColsType>().isNoSameBunchPileup())
        continue;
      if (_requestGoodZvtxFT0vsPV && !track.template singleCollSel_as<ColsType>().isGoodZvtxFT0vsPV())
        continue;
      if (_requestVertexITSTPC && !track.template singleCollSel_as<ColsType>().isVertexITSTPC())
        continue;
      if (_requestVertexTOForTRDmatched > track.template singleCollSel_as<ColsType>().isVertexTOForTRDmatched())
        continue;
      if (_requestNoCollInTimeRangeStandard && !track.template singleCollSel_as<ColsType>().noCollInTimeRangeStandard())
        continue;

      if (std::fabs(track.template singleCollSel_as<ColsType>().posZ()) > _vertexZ)
        continue;
      if (track.template singleCollSel_as<ColsType>().multPerc() < _centCut.value.first || track.template singleCollSel_as<ColsType>().multPerc() >= _centCut.value.second)
        continue;
      if (track.template singleCollSel_as<ColsType>().hadronicRate() < _IRcut.value.first || track.template singleCollSel_as<ColsType>().hadronicRate() >= _IRcut.value.second)
        continue;
      if (track.template singleCollSel_as<ColsType>().occupancy() < _OccupancyCut.value.first || track.template singleCollSel_as<ColsType>().occupancy() >= _OccupancyCut.value.second)
        continue;
      if ((track.tpcFractionSharedCls()) > _tpcFractionSharedCls || (track.itsNCls()) < _itsNCls)
        continue;
      if (std::fabs(track.dcaXY()) > _dcaXY.value[0] + _dcaXY.value[1] * std::pow(track.pt(), _dcaXY.value[2]) || std::fabs(track.dcaZ()) > _dcaZ.value[0] + _dcaZ.value[1] * std::pow(track.pt(), _dcaZ.value[2]))
        continue;

      if constexpr (FillExtra) {
        registry.fill(HIST("TPCSignal_nocuts"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("TOFSignal_nocuts"), track.p(), track.beta());
      }

      if (!TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF)) && (track.p() < _PIDtrshld ? o2::aod::singletrackselector::TPCselection(track, TPCcuts) : o2::aod::singletrackselector::TOFselection(track, TOFcuts, _tpcNSigmaResidual.value))) {
        registry.fill(HIST("eta"), track.eta());
        registry.fill(HIST("phi"), track.phi());
        registry.fill(HIST("px"), track.px());
        registry.fill(HIST("py"), track.py());
        registry.fill(HIST("pz"), track.pz());
        registry.fill(HIST("p"), track.p());
        registry.fill(HIST("pt"), track.pt());
        registry.fill(HIST("sign"), track.sign());
        registry.fill(HIST("dcaxy_to_p"), track.p(), track.dcaXY());
        registry.fill(HIST("dcaxy_to_pt"), track.pt(), track.dcaXY());
        registry.fill(HIST("dcaz_to_p"), track.p(), track.dcaZ());
        registry.fill(HIST("dcaz_to_pt"), track.pt(), track.dcaZ());
        registry.fill(HIST("TPCClusters"), track.tpcNClsFound());
        registry.fill(HIST("TPCCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        registry.fill(HIST("TPCFractionSharedCls"), track.tpcFractionSharedCls());
        registry.fill(HIST("ITSClusters"), track.itsNCls());
        registry.fill(HIST("ITSchi2"), track.itsChi2NCl());
        registry.fill(HIST("TPCchi2"), track.tpcChi2NCl());

        switch (_particlePDG) {
          case 211:
            registry.fill(HIST("nsigmaTOFPi"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("nsigmaTPCPi"), track.p(), track.tpcNSigmaPi());
            registry.fill(HIST("nsigmaITSPi"), track.p(), track.itsNSigmaPi());
            break;
          case 321:
            registry.fill(HIST("nsigmaTOFKa"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("nsigmaTPCKa"), track.p(), track.tpcNSigmaKa());
            registry.fill(HIST("nsigmaITSKa"), track.p(), track.itsNSigmaKa());
            break;
          case 2212:
            registry.fill(HIST("nsigmaTOFPr"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPCPr"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("nsigmaITSPr"), track.p(), track.itsNSigmaPr());
            break;
          case 1000010020:
            registry.fill(HIST("nsigmaTOFDe"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPCDe"), track.p(), track.tpcNSigmaDe());
            registry.fill(HIST("nsigmaITSDe"), track.p(), track.itsNSigmaDe());
            break;
          case 1000010030:
            registry.fill(HIST("nsigmaTOFTr"), track.p(), track.tofNSigmaTr());
            registry.fill(HIST("nsigmaTPCTr"), track.p(), track.tpcNSigmaTr());
            registry.fill(HIST("nsigmaITSTr"), track.p(), track.itsNSigmaTr());
            break;
          case 1000020030:
            registry.fill(HIST("nsigmaTOFHe"), track.p(), track.tofNSigmaHe());
            registry.fill(HIST("nsigmaTPCHe"), track.p(), track.tpcNSigmaHe());
            registry.fill(HIST("nsigmaITSHe"), track.p(), track.itsNSigmaHe());
            break;
          default:
            break;
        }

        if constexpr (FillExtra) {
          registry.fill(HIST("TPCSignal"), track.tpcInnerParam(), track.tpcSignal());
          registry.fill(HIST("TOFSignal"), track.p(), track.beta());
        }
      }
    }
  }

  void processDefault(soa::Filtered<soa::Join<aod::SingleCollSels, aod::SingleCollExtras>> const& collisions, soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SinglePIDPis, aod::SinglePIDKas, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes>> const& tracks)
  {
    fillHistograms<false>(collisions, tracks);
  }
  PROCESS_SWITCH(QAHistograms, processDefault, "process default", true);

  void processExtra(soa::Filtered<soa::Join<aod::SingleCollSels, aod::SingleCollExtras>> const& collisions, soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SingleTrkExtras, aod::SinglePIDPis, aod::SinglePIDKas, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes>> const& tracks)
  {
    fillHistograms<true>(collisions, tracks);
  }
  PROCESS_SWITCH(QAHistograms, processExtra, "process extra", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QAHistograms>(cfgc)};
}
